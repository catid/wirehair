#include "Wh2FrozenTrace.h"
#include "Wh2K2ProjectiveCodec.h"
#include "Wh2NativeCodec.h"

#include <wirehair/wirehair.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <climits>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#if defined(__linux__)
#include <sched.h>
#include <unistd.h>
#endif

namespace {

using wirehair_wh2_bench::K2ProjectiveCodec;
using wirehair_wh2_bench::K2ProjectiveDirection;
using wirehair_wh2_bench::NativeArm;
using wirehair_wh2_bench::NativeArmSpec;
using wirehair_wh2_bench::NativeEncoderFixture;
using wirehair_wh2_bench::NativeReceiveFixture;
using wirehair_wh2_bench::TimedArmResult;

static const uint32_t kK = 2u;
static const uint32_t kLossPpm = 100000u;
static const uint32_t kReceiveOverheadCap = 256u;
static const uint32_t kWarmupInvocations = 1u;
static const uint32_t kObservationCount = 32u;
static const uint32_t kInvocationsPerObservation = 8u;
static const uint32_t kBlockBytes[] = { 64u, 1280u };
static const uint64_t kTrainingRoots[] = {
    UINT64_C(0xd1b54a32d192ed03),
    UINT64_C(0x94d049bb133111eb),
    UINT64_C(0x8538ecb5bd456ea3)
};

static const char kDevelopmentDomainSha256[] =
    "f97f28c211428cd77aed97160073b192d93014cb4a61a844bc7d76375ac61b77";
static const char kK2CanonicalCellStreamSha256[] =
    "85195a20cd8eed3c60b78487e5436e3a4574e58fac27a084de22d0e550a9143f";
static const char kK2TraceStreamSha256[] =
    "3ed52a503688c75cbb3e6a4735c8a03449eec8aa351ca73c8546187c6186da52";
static const char kBoundaryDirectionSha256[] =
    "d842960b16718ac7afae1c6ed536bb2abaf1d3924e333f0783ad2d12c55117cc";
static const char kHighDirectionSha256[] =
    "e7f0e80ebbc2c2e93c8193da6fbfb8ae48012cbf05adddf26e4f5b60239ab10f";
static const char kCurrentWh2Comparison[] =
    "candidate/current-wh2-selected-attempt-diagnostic-only";

template <bool Measure>
class ScreenTimer;

template <>
class ScreenTimer<true>
{
public:
    ScreenTimer()
        : Start(std::chrono::steady_clock::now())
    {
    }

    uint64_t Stop() const
    {
        const std::chrono::steady_clock::time_point finish =
            std::chrono::steady_clock::now();
        const std::chrono::nanoseconds elapsed =
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                finish - Start);
        return elapsed.count() > 0 ?
            static_cast<uint64_t>(elapsed.count()) : 0u;
    }

private:
    std::chrono::steady_clock::time_point Start;
};

template <>
class ScreenTimer<false>
{
public:
    uint64_t Stop() const
    {
        return 0u;
    }
};

enum class Scope : uint32_t
{
    Encoder = 0u,
    IidReceive = 1u,
    TwoRepairReceive = 2u
};

enum ArmIndex : uint32_t
{
    CandidateArm = 0u,
    Wirehair1Arm = 1u,
    CurrentWh2Arm = 2u,
    Wirehair1AaArm = 3u,
    ArmCount = 4u
};

enum class HarnessOutcome : uint32_t
{
    Complete = 0u,
    CandidateReject,
    Invalid
};

const char* ScopeName(Scope scope)
{
    switch (scope)
    {
    case Scope::Encoder: return "encoder-init-systematic";
    case Scope::IidReceive: return "iid10-receive-to-success";
    case Scope::TwoRepairReceive: return "two-repair-receive-to-success";
    }
    return "invalid";
}

const char* ArmName(ArmIndex arm)
{
    switch (arm)
    {
    case CandidateArm: return "candidate";
    case Wirehair1Arm: return "wirehair1";
    case CurrentWh2Arm: return "current-wh2-selected-attempt";
    case Wirehair1AaArm: return "wirehair1-aa";
    case ArmCount: break;
    }
    return "invalid";
}

bool CheckedPacketStorageBytes(
    size_t packet_count,
    uint32_t block_bytes,
    size_t& bytes_out)
{
    bytes_out = 0u;
    if (block_bytes == 0u ||
        packet_count > std::numeric_limits<size_t>::max() / block_bytes)
    {
        return false;
    }
    bytes_out = packet_count * block_bytes;
    return true;
}

void AppendLe16(std::string& bytes, uint16_t value)
{
    bytes.push_back(static_cast<char>(value & 0xffu));
    bytes.push_back(static_cast<char>((value >> 8) & 0xffu));
}

void AppendLe32(std::string& bytes, uint32_t value)
{
    for (unsigned shift = 0u; shift != 32u; shift += 8u) {
        bytes.push_back(static_cast<char>((value >> shift) & 0xffu));
    }
}

void AppendLe64(std::string& bytes, uint64_t value)
{
    for (unsigned shift = 0u; shift != 64u; shift += 8u) {
        bytes.push_back(static_cast<char>((value >> shift) & 0xffu));
    }
}

bool SameDirection(
    const K2ProjectiveDirection& left,
    const K2ProjectiveDirection& right)
{
    return left.Point == right.Point &&
        left.Alpha == right.Alpha && left.Beta == right.Beta;
}

bool IndependentDirections(uint32_t first_id, uint32_t second_id)
{
    return wirehair_wh2_bench::K2ProjectiveDirectionForPacketId(first_id).
               Point !=
        wirehair_wh2_bench::K2ProjectiveDirectionForPacketId(second_id).
               Point;
}

struct AliasGroup
{
    uint16_t Point = 0u;
    std::vector<uint32_t> Ids;
};

bool FindFirstHighAliasGroup(size_t required_ids, AliasGroup& group_out)
{
    if (required_ids < 2u) return false;
    std::map<uint16_t, std::vector<uint32_t> > groups;
    for (uint32_t id = 257u; id < 100000u; ++id)
    {
        const K2ProjectiveDirection direction =
            wirehair_wh2_bench::K2ProjectiveDirectionForPacketId(id);
        std::vector<uint32_t>& ids = groups[direction.Point];
        ids.push_back(id);
        if (ids.size() >= required_ids)
        {
            AliasGroup result;
            result.Point = direction.Point;
            result.Ids.assign(ids.begin(), ids.begin() + required_ids);
            group_out = result;
            return true;
        }
    }
    return false;
}

bool EncodeCandidatePacket(
    const K2ProjectiveCodec& encoder,
    uint32_t block_id,
    std::vector<uint8_t>& packet_out,
    uint32_t& data_bytes_out)
{
    std::vector<uint8_t> packet(encoder.BlockBytes(), uint8_t{0xa5});
    uint32_t data_bytes = UINT32_MAX;
    const WirehairResult result = encoder.EncodeResult(
        block_id,
        packet.data(),
        static_cast<uint32_t>(packet.size()),
        &data_bytes);
    if (result != Wirehair_Success || data_bytes > packet.size()) {
        return false;
    }
    packet.resize(data_bytes);
    packet_out.swap(packet);
    data_bytes_out = data_bytes;
    return true;
}

class CandidateEncoderFixture
{
public:
    bool Initialize(
        uint32_t block_bytes,
        const std::vector<uint8_t>& source)
    {
        if (block_bytes == 0u ||
            source.size() != static_cast<size_t>(kK) * block_bytes)
        {
            return false;
        }
        BlockBytes = block_bytes;
        Source = source;
        Initialized = true;
        return true;
    }

    TimedArmResult Preflight() const
    {
        return RunImpl<false>();
    }

    TimedArmResult Run() const
    {
        return RunImpl<true>();
    }

private:
    template <bool Measure>
    TimedArmResult RunImpl() const
    {
        TimedArmResult result;
        if (!Initialized) return result;
        try
        {
            std::vector<uint8_t> owned_source(Source);
            std::vector<uint8_t> encoded(Source.size(), uint8_t{0});
            K2ProjectiveCodec encoder;
            WirehairResult invocation_result = Wirehair_Error;
            const ScreenTimer<Measure> timer;
            invocation_result =
                encoder.InitializeEncoderOwnedSourceAfterGlobalInit(
                    std::move(owned_source), BlockBytes);
            if (invocation_result == Wirehair_Success)
            {
                for (uint32_t id = 0u; id < kK; ++id)
                {
                    uint32_t written = 0u;
                    invocation_result = encoder.EncodeResult(
                        id,
                        encoded.data() +
                            static_cast<size_t>(id) * BlockBytes,
                        BlockBytes,
                        &written);
                    if (invocation_result != Wirehair_Success ||
                        written != BlockBytes)
                    {
                        invocation_result = Wirehair_Error;
                        break;
                    }
                }
            }
            result.ElapsedNanoseconds = timer.Stop();
            result.Result = invocation_result;
            result.BytesVerified = encoded == Source;
            if (result.Result == Wirehair_Success &&
                !result.BytesVerified) {
                result.Result = Wirehair_Error;
            }
            if (result.Result != Wirehair_Success) {
                result.BytesVerified = false;
            }
            return result;
        }
        catch (const std::bad_alloc&) {
            result.Result = Wirehair_OOM;
            return result;
        }
        catch (const std::length_error&) {
            result.Result = Wirehair_OOM;
            return result;
        }
    }

    uint32_t BlockBytes = 0u;
    std::vector<uint8_t> Source;
    bool Initialized = false;
};

class CandidateReceiveFixture
{
public:
    WirehairResult Initialize(
        uint32_t block_bytes,
        const std::vector<uint8_t>& source,
        const std::vector<uint32_t>& packet_ids)
    {
        if (block_bytes == 0u ||
            source.size() != static_cast<size_t>(kK) * block_bytes ||
            packet_ids.size() < kK)
        {
            return Wirehair_InvalidInput;
        }
        size_t storage_bytes = 0u;
        if (!CheckedPacketStorageBytes(
                packet_ids.size(), block_bytes, storage_bytes))
        {
            return Wirehair_OOM;
        }
        try
        {
            K2ProjectiveCodec encoder;
            const WirehairResult initialize_result = encoder.InitializeEncoder(
                source.data(), source.size(), block_bytes);
            if (initialize_result != Wirehair_Success) {
                return initialize_result;
            }
            std::vector<uint8_t> packets(storage_bytes, uint8_t{0});
            std::vector<uint32_t> packet_bytes(packet_ids.size(), 0u);
            for (size_t i = 0u; i < packet_ids.size(); ++i)
            {
                uint32_t written = 0u;
                const WirehairResult encode_result = encoder.EncodeResult(
                    packet_ids[i],
                    packets.data() + i * block_bytes,
                    block_bytes,
                    &written);
                if (encode_result != Wirehair_Success) {
                    return encode_result;
                }
                if (written == 0u || written > block_bytes)
                {
                    return Wirehair_Error;
                }
                packet_bytes[i] = written;
            }
            BlockBytes = block_bytes;
            Source = source;
            PacketIds = packet_ids;
            PacketStorage.swap(packets);
            PacketBytes.swap(packet_bytes);
            Initialized = true;
            return Wirehair_Success;
        }
        catch (const std::bad_alloc&) {
            return Wirehair_OOM;
        }
        catch (const std::length_error&) {
            return Wirehair_OOM;
        }
    }

    TimedArmResult Preflight() const
    {
        return RunImpl<false>();
    }

    TimedArmResult Run() const
    {
        return RunImpl<true>();
    }

private:
    template <bool Measure>
    TimedArmResult RunImpl() const
    {
        TimedArmResult result;
        if (!Initialized) return result;
        try
        {
            std::vector<uint8_t> recovered(Source.size(), uint8_t{0});
            K2ProjectiveCodec decoder;
            const WirehairResult init_result = decoder.InitializeDecoder(
                Source.size(), BlockBytes);
            if (init_result != Wirehair_Success) {
                result.Result = init_result;
                return result;
            }

            WirehairResult receive_result = Wirehair_NeedMore;
            uint32_t success_count = 0u;
            const ScreenTimer<Measure> timer;
            for (size_t i = 0u; i < PacketIds.size(); ++i)
            {
                receive_result = decoder.DecodeResult(
                    PacketIds[i],
                    PacketStorage.data() + i * BlockBytes,
                    PacketBytes[i]);
                if (receive_result == Wirehair_Success)
                {
                    success_count = static_cast<uint32_t>(i + 1u);
                    receive_result = decoder.RecoverResult(
                        recovered.data(), recovered.size());
                    break;
                }
                if (receive_result != Wirehair_NeedMore) break;
            }
            result.ElapsedNanoseconds = timer.Stop();
            result.Result = receive_result;
            result.BytesVerified = recovered == Source;
            if (success_count >= kK) {
                result.DecodedOverhead = success_count - kK;
            }
            if (result.Result == Wirehair_Success &&
                (!result.BytesVerified ||
                 result.DecodedOverhead == UINT32_MAX))
            {
                result.Result = Wirehair_Error;
            }
            if (result.Result != Wirehair_Success)
            {
                result.BytesVerified = false;
                result.DecodedOverhead = UINT32_MAX;
            }
            return result;
        }
        catch (const std::bad_alloc&) {
            result.Result = Wirehair_OOM;
            return result;
        }
        catch (const std::length_error&) {
            result.Result = Wirehair_OOM;
            return result;
        }
    }

    uint32_t BlockBytes = 0u;
    std::vector<uint8_t> Source;
    std::vector<uint32_t> PacketIds;
    std::vector<uint8_t> PacketStorage;
    std::vector<uint32_t> PacketBytes;
    bool Initialized = false;
};

struct ReceiveFixtures
{
    CandidateReceiveFixture Candidate;
    NativeReceiveFixture Wirehair1;
    NativeReceiveFixture CurrentWh2;
};

struct RootCell
{
    uint32_t BlockBytes = 0u;
    uint32_t CurrentWh2Attempt = UINT32_MAX;
    uint64_t Root = 0u;
    std::vector<uint8_t> Source;
    std::vector<uint32_t> IidIds;
    std::vector<uint32_t> RepairIds;
    CandidateEncoderFixture CandidateEncoder;
    NativeEncoderFixture Wirehair1Encoder;
    NativeEncoderFixture CurrentWh2Encoder;
    NativeArm Wirehair1State;
    NativeArm CurrentWh2State;
    ReceiveFixtures IidReceive;
    ReceiveFixtures RepairReceive;
};

bool CurrentWh2SelectedAttempt(
    uint32_t block_bytes,
    uint32_t& attempt_out)
{
    if (block_bytes == 64u) {
        attempt_out = 7u;
        return true;
    }
    if (block_bytes == 1280u) {
        attempt_out = 1u;
        return true;
    }
    return false;
}

bool ExpectedTwoRepairOverhead(
    uint32_t block_bytes,
    ArmIndex arm,
    uint32_t& overhead_out)
{
    if (block_bytes != 64u && block_bytes != 1280u) return false;
    overhead_out = arm == CurrentWh2Arm && block_bytes == 1280u ? 1u : 0u;
    return true;
}

bool SemanticallyValidResult(
    uint32_t block_bytes,
    Scope scope,
    ArmIndex arm,
    const TimedArmResult& result)
{
    if (result.Result != Wirehair_Success || !result.BytesVerified) {
        return false;
    }
    if (scope == Scope::Encoder) {
        return result.DecodedOverhead == UINT32_MAX;
    }
    if (scope == Scope::TwoRepairReceive) {
        uint32_t expected_overhead = 0u;
        return ExpectedTwoRepairOverhead(
                   block_bytes, arm, expected_overhead) &&
            result.DecodedOverhead == expected_overhead;
    }
    return result.DecodedOverhead != UINT32_MAX &&
        result.DecodedOverhead <= kReceiveOverheadCap;
}

HarnessOutcome ClassifyInvocation(
    uint32_t block_bytes,
    Scope scope,
    ArmIndex arm,
    const TimedArmResult& result,
    bool measured,
    const char* phase,
    std::string& diagnostic)
{
    if (!SemanticallyValidResult(block_bytes, scope, arm, result))
    {
        std::ostringstream message;
        message << ScopeName(scope) << ' ' << ArmName(arm) << ' ' << phase
                << " semantic failure: result="
                << static_cast<uint32_t>(result.Result)
                << " bytes_verified=" << (result.BytesVerified ? 1 : 0)
                << " decoded_overhead=" << result.DecodedOverhead;
        if (scope == Scope::TwoRepairReceive)
        {
            uint32_t expected_overhead = UINT32_MAX;
            if (ExpectedTwoRepairOverhead(
                    block_bytes, arm, expected_overhead))
            {
                message << " expected_decoded_overhead="
                        << expected_overhead;
            }
        }
        diagnostic = message.str();
        if (arm != CandidateArm || result.Result == Wirehair_OOM ||
            result.Result == Wirehair_UnsupportedPlatform)
        {
            return HarnessOutcome::Invalid;
        }
        return HarnessOutcome::CandidateReject;
    }
    const bool elapsed_valid = measured ?
        result.ElapsedNanoseconds > 0u :
        result.ElapsedNanoseconds == 0u;
    if (!elapsed_valid)
    {
        std::ostringstream message;
        message << ScopeName(scope) << ' ' << ArmName(arm) << ' ' << phase
                << " clock-envelope failure: elapsed_ns="
                << result.ElapsedNanoseconds;
        diagnostic = message.str();
        return HarnessOutcome::Invalid;
    }
    return HarnessOutcome::Complete;
}

HarnessOutcome BuildReceiveFixtures(
    RootCell& cell,
    const std::vector<uint32_t>& ids,
    const char* trace_name,
    ReceiveFixtures& fixtures,
    std::string& diagnostic)
{
    const WirehairResult candidate_result = fixtures.Candidate.Initialize(
        cell.BlockBytes, cell.Source, ids);
    if (candidate_result != Wirehair_Success)
    {
        std::ostringstream message;
        message << "candidate " << trace_name
                << " packet construction failed: result="
                << static_cast<uint32_t>(candidate_result);
        diagnostic = message.str();
        return candidate_result == Wirehair_OOM ||
                candidate_result == Wirehair_UnsupportedPlatform ?
            HarnessOutcome::Invalid : HarnessOutcome::CandidateReject;
    }
    const WirehairResult wirehair1_result = fixtures.Wirehair1.Initialize(
        cell.Wirehair1State, ids, kReceiveOverheadCap);
    if (wirehair1_result != Wirehair_Success)
    {
        std::ostringstream message;
        message << "Wirehair1 " << trace_name
                << " fixture construction failed: result="
                << static_cast<uint32_t>(wirehair1_result);
        diagnostic = message.str();
        return HarnessOutcome::Invalid;
    }
    const WirehairResult current_wh2_result =
        fixtures.CurrentWh2.Initialize(
            cell.CurrentWh2State, ids, kReceiveOverheadCap);
    if (current_wh2_result != Wirehair_Success)
    {
        std::ostringstream message;
        message << "current WH2 " << trace_name
                << " fixture construction failed: result="
                << static_cast<uint32_t>(current_wh2_result);
        diagnostic = message.str();
        return HarnessOutcome::Invalid;
    }
    return HarnessOutcome::Complete;
}

HarnessOutcome BuildRootCell(
    uint32_t block_bytes,
    uint64_t root,
    RootCell& cell_out,
    std::string& diagnostic)
{
    diagnostic.clear();
    std::unique_ptr<RootCell> cell(new RootCell);
    cell->BlockBytes = block_bytes;
    cell->Root = root;
    if (!CurrentWh2SelectedAttempt(
            block_bytes, cell->CurrentWh2Attempt))
    {
        diagnostic = "block width has no frozen selected WH2 attempt";
        return HarnessOutcome::Invalid;
    }
    if (!wirehair_wh2_bench::MakeDeterministicSource(
            kK, block_bytes, root, cell->Source))
    {
        diagnostic = "deterministic source generation failed";
        return HarnessOutcome::Invalid;
    }

    wirehair::wh2_benchmark::FrozenPacketTrace iid_trace;
    if (wirehair::wh2_benchmark::GenerateFrozenPacketTrace(
            kK,
            block_bytes,
            kLossPpm,
            root,
            wirehair::wh2_benchmark::FrozenSchedule::Iid,
            kReceiveOverheadCap,
            iid_trace) !=
            wirehair::wh2_benchmark::FrozenTraceStatus::Complete ||
        iid_trace.delivered_ids.size() !=
            static_cast<size_t>(kK) + kReceiveOverheadCap)
    {
        diagnostic = "IID10 K+256 trace generation failed";
        return HarnessOutcome::Invalid;
    }
    cell->IidIds = iid_trace.delivered_ids;
    cell->RepairIds.resize(static_cast<size_t>(kK) + kReceiveOverheadCap);
    for (size_t i = 0u; i < cell->RepairIds.size(); ++i) {
        cell->RepairIds[i] = static_cast<uint32_t>(i) + 2u;
    }

    const NativeArmSpec wh1_spec = wirehair_wh2_bench::MakeWirehair1Arm();
    const NativeArmSpec wh2_spec =
        wirehair_wh2_bench::MakeCertifiedWh2Arm(cell->CurrentWh2Attempt);
    if (cell->Wirehair1State.Initialize(
            wh1_spec, kK, block_bytes, cell->Source) != Wirehair_Success)
    {
        diagnostic = "Wirehair1 exact state construction failed";
        return HarnessOutcome::Invalid;
    }
    if (cell->CurrentWh2State.Initialize(
            wh2_spec, kK, block_bytes, cell->Source) != Wirehair_Success)
    {
        diagnostic = "current WH2 exact state construction failed";
        return HarnessOutcome::Invalid;
    }
    if (cell->CurrentWh2State.ConstructionAttempt() !=
        cell->CurrentWh2Attempt)
    {
        diagnostic = "current WH2 state did not retain selected attempt";
        return HarnessOutcome::Invalid;
    }
    if (!cell->CandidateEncoder.Initialize(block_bytes, cell->Source) ||
        cell->Wirehair1Encoder.Initialize(
            wh1_spec, kK, block_bytes, cell->Source) != Wirehair_Success ||
        cell->CurrentWh2Encoder.Initialize(
            wh2_spec, kK, block_bytes, cell->Source) != Wirehair_Success)
    {
        diagnostic = "encoder fixture construction failed";
        return HarnessOutcome::Invalid;
    }
    const HarnessOutcome iid_outcome = BuildReceiveFixtures(
        *cell, cell->IidIds, "IID", cell->IidReceive, diagnostic);
    if (iid_outcome != HarnessOutcome::Complete) {
        return iid_outcome;
    }
    const HarnessOutcome repair_outcome = BuildReceiveFixtures(
        *cell,
        cell->RepairIds,
        "two-repair",
        cell->RepairReceive,
        diagnostic);
    if (repair_outcome != HarnessOutcome::Complete)
    {
        return repair_outcome;
    }

    cell_out = std::move(*cell);
    return HarnessOutcome::Complete;
}

TimedArmResult RunArm(
    const RootCell& cell,
    Scope scope,
    ArmIndex arm)
{
    if (scope == Scope::Encoder)
    {
        if (arm == CandidateArm) return cell.CandidateEncoder.Run();
        if (arm == CurrentWh2Arm) return cell.CurrentWh2Encoder.Run();
        return cell.Wirehair1Encoder.Run();
    }
    const ReceiveFixtures& fixtures = scope == Scope::IidReceive ?
        cell.IidReceive : cell.RepairReceive;
    if (arm == CandidateArm) return fixtures.Candidate.Run();
    if (arm == CurrentWh2Arm) return fixtures.CurrentWh2.Run();
    return fixtures.Wirehair1.Run();
}

TimedArmResult PreflightArm(
    const RootCell& cell,
    Scope scope,
    ArmIndex arm)
{
    if (scope == Scope::Encoder)
    {
        if (arm == CandidateArm) return cell.CandidateEncoder.Preflight();
        if (arm == CurrentWh2Arm) {
            return cell.CurrentWh2Encoder.Preflight();
        }
        return cell.Wirehair1Encoder.Preflight();
    }
    const ReceiveFixtures& fixtures = scope == Scope::IidReceive ?
        cell.IidReceive : cell.RepairReceive;
    if (arm == CandidateArm) return fixtures.Candidate.Preflight();
    if (arm == CurrentWh2Arm) return fixtures.CurrentWh2.Preflight();
    return fixtures.Wirehair1.Preflight();
}

HarnessOutcome PreflightRootCell(
    const RootCell& cell,
    std::string& diagnostic)
{
    diagnostic.clear();
    const Scope scopes[] = {
        Scope::Encoder, Scope::IidReceive, Scope::TwoRepairReceive
    };
    for (Scope scope : scopes)
    {
        for (uint32_t arm_value = 0u;
             arm_value < static_cast<uint32_t>(ArmCount);
             ++arm_value)
        {
            const ArmIndex arm = static_cast<ArmIndex>(arm_value);
            const TimedArmResult result = PreflightArm(cell, scope, arm);
            const HarnessOutcome outcome = ClassifyInvocation(
                cell.BlockBytes,
                scope,
                arm,
                result,
                false,
                "preflight",
                diagnostic);
            if (outcome != HarnessOutcome::Complete) return outcome;
        }
    }
    return HarnessOutcome::Complete;
}

struct LogRatioStats
{
    size_t Count = 0u;
    double GeometricMean = 0.0;
    double Lower95 = 0.0;
    double Upper95 = 0.0;
};

double StudentT975(size_t sample_count)
{
    if (sample_count == 32u) return 2.039513446;
    if (sample_count == 96u) return 1.985251004;
    if (sample_count == 576u) return 1.964095;
    return 2.0;
}

bool ComputeLogRatioStats(
    const std::vector<double>& logs,
    LogRatioStats& stats_out)
{
    if (logs.size() < 2u) return false;
    long double sum = 0.0L;
    for (double value : logs) {
        if (!std::isfinite(value)) return false;
        sum += static_cast<long double>(value);
    }
    const long double mean = sum / logs.size();
    long double squares = 0.0L;
    for (double value : logs)
    {
        const long double difference =
            static_cast<long double>(value) - mean;
        squares += difference * difference;
    }
    const long double variance = squares / (logs.size() - 1u);
    const long double standard_error =
        std::sqrt(variance / logs.size());
    const long double half_width =
        StudentT975(logs.size()) * standard_error;
    LogRatioStats stats;
    stats.Count = logs.size();
    stats.GeometricMean = static_cast<double>(std::exp(mean));
    stats.Lower95 = static_cast<double>(std::exp(mean - half_width));
    stats.Upper95 = static_cast<double>(std::exp(mean + half_width));
    if (!std::isfinite(stats.GeometricMean) ||
        !std::isfinite(stats.Lower95) ||
        !std::isfinite(stats.Upper95) || stats.Lower95 <= 0.0 ||
        stats.Lower95 > stats.GeometricMean ||
        stats.GeometricMean > stats.Upper95)
    {
        return false;
    }
    stats_out = stats;
    return true;
}

struct MeasurementLogs
{
    std::vector<double> CandidateVsWirehair1;
    std::vector<double> CandidateVsCurrentWh2;
    std::vector<double> Wirehair1Aa;
    std::vector<std::array<uint64_t, ArmCount> > Observations;
};

struct PanelMeasurement
{
    size_t WidthIndex = 0u;
    size_t RootIndex = 0u;
    Scope ScopeValue = Scope::Encoder;
    MeasurementLogs Logs;
};

void ReserveMeasurementLogs(MeasurementLogs& logs)
{
    logs.CandidateVsWirehair1.reserve(kObservationCount);
    logs.CandidateVsCurrentWh2.reserve(kObservationCount);
    logs.Wirehair1Aa.reserve(kObservationCount);
    logs.Observations.reserve(kObservationCount);
}

static const ArmIndex kCounterbalancedOrders[4][4] = {
    { CandidateArm, Wirehair1Arm, Wirehair1AaArm, CurrentWh2Arm },
    { Wirehair1Arm, CurrentWh2Arm, CandidateArm, Wirehair1AaArm },
    { CurrentWh2Arm, Wirehair1AaArm, Wirehair1Arm, CandidateArm },
    { Wirehair1AaArm, CandidateArm, CurrentWh2Arm, Wirehair1Arm }
};

bool AddElapsed(uint64_t value, uint64_t& total)
{
    if (value > std::numeric_limits<uint64_t>::max() - total) return false;
    total += value;
    return true;
}

HarnessOutcome WarmScope(
    const RootCell& cell,
    Scope scope,
    uint32_t order_offset,
    std::string& diagnostic)
{
    for (uint32_t invocation = 0u;
         invocation < kWarmupInvocations;
         ++invocation)
    {
        const uint32_t row = (invocation + order_offset) % 4u;
        for (uint32_t position = 0u; position < 4u; ++position)
        {
            const ArmIndex arm = kCounterbalancedOrders[row][position];
            const HarnessOutcome outcome = ClassifyInvocation(
                cell.BlockBytes,
                scope,
                arm,
                RunArm(cell, scope, arm),
                true,
                "warmup",
                diagnostic);
            if (outcome != HarnessOutcome::Complete) return outcome;
        }
    }
    return HarnessOutcome::Complete;
}

HarnessOutcome MeasureScope(
    const RootCell& cell,
    Scope scope,
    uint32_t order_offset,
    MeasurementLogs& logs_out,
    std::string& diagnostic)
{
    if (!logs_out.CandidateVsWirehair1.empty() ||
        !logs_out.CandidateVsCurrentWh2.empty() ||
        !logs_out.Wirehair1Aa.empty() || !logs_out.Observations.empty())
    {
        diagnostic = "measurement storage was not empty";
        return HarnessOutcome::Invalid;
    }
    const HarnessOutcome warmup =
        WarmScope(cell, scope, order_offset, diagnostic);
    if (warmup != HarnessOutcome::Complete) return warmup;

    for (uint32_t observation = 0u;
         observation < kObservationCount;
         ++observation)
    {
        uint64_t totals[ArmCount] = {};
        for (uint32_t invocation = 0u;
             invocation < kInvocationsPerObservation;
             ++invocation)
        {
            const uint32_t row =
                (observation * kInvocationsPerObservation + invocation +
                    order_offset) % 4u;
            for (uint32_t position = 0u; position < 4u; ++position)
            {
                const ArmIndex arm = kCounterbalancedOrders[row][position];
                const TimedArmResult result = RunArm(cell, scope, arm);
                const HarnessOutcome outcome = ClassifyInvocation(
                    cell.BlockBytes,
                    scope,
                    arm,
                    result,
                    true,
                    "measurement",
                    diagnostic);
                if (outcome != HarnessOutcome::Complete) return outcome;
                if (!AddElapsed(
                        result.ElapsedNanoseconds,
                        totals[static_cast<uint32_t>(arm)]))
                {
                    diagnostic = "aggregate duration overflow";
                    return HarnessOutcome::Invalid;
                }
            }
        }
        for (uint32_t arm = 0u; arm < static_cast<uint32_t>(ArmCount); ++arm)
        {
            if (totals[arm] == 0u) {
                diagnostic = "zero aggregate duration";
                return HarnessOutcome::Invalid;
            }
        }
        std::array<uint64_t, ArmCount> record = {};
        for (uint32_t arm = 0u; arm < static_cast<uint32_t>(ArmCount); ++arm) {
            record[arm] = totals[arm];
        }
        logs_out.Observations.push_back(record);
    }
    return HarnessOutcome::Complete;
}

bool FinalizeMeasurementLogs(
    MeasurementLogs& logs,
    std::string& diagnostic)
{
    if (logs.Observations.size() != kObservationCount ||
        !logs.CandidateVsWirehair1.empty() ||
        !logs.CandidateVsCurrentWh2.empty() ||
        !logs.Wirehair1Aa.empty())
    {
        diagnostic = "measurement observation cardinality invalid";
        return false;
    }
    for (const std::array<uint64_t, ArmCount>& totals : logs.Observations)
    {
        const double candidate_vs_wh1 = std::log(
            static_cast<double>(totals[CandidateArm]) /
            static_cast<double>(totals[Wirehair1Arm]));
        const double candidate_vs_wh2 = std::log(
            static_cast<double>(totals[CandidateArm]) /
            static_cast<double>(totals[CurrentWh2Arm]));
        const double aa = std::log(
            static_cast<double>(totals[Wirehair1Arm]) /
            static_cast<double>(totals[Wirehair1AaArm]));
        if (!std::isfinite(candidate_vs_wh1) ||
            !std::isfinite(candidate_vs_wh2) || !std::isfinite(aa))
        {
            diagnostic = "non-finite paired log ratio";
            return false;
        }
        logs.CandidateVsWirehair1.push_back(candidate_vs_wh1);
        logs.CandidateVsCurrentWh2.push_back(candidate_vs_wh2);
        logs.Wirehair1Aa.push_back(aa);
    }
    return true;
}

void AppendLogs(
    const std::vector<double>& source,
    std::vector<double>& destination)
{
    destination.insert(destination.end(), source.begin(), source.end());
}

void PrintObservations(
    const RootCell& cell,
    Scope scope,
    const MeasurementLogs& logs)
{
    for (size_t observation = 0u;
         observation < logs.Observations.size();
         ++observation)
    {
        const std::array<uint64_t, ArmCount>& totals =
            logs.Observations[observation];
        std::cout << "OBS,scope=" << ScopeName(scope)
                  << ",block_bytes=" << cell.BlockBytes
                  << ",root=0x" << std::hex << cell.Root << std::dec
                  << ",observation=" << observation
                  << ",candidate_ns=" << totals[CandidateArm]
                  << ",wirehair1_ns=" << totals[Wirehair1Arm]
                  << ",current_wh2_ns=" << totals[CurrentWh2Arm]
                  << ",current_wh2_attempt=" << cell.CurrentWh2Attempt
                  << ",current_wh2_diagnostic_only=1"
                  << ",wirehair1_aa_ns=" << totals[Wirehair1AaArm]
                  ;
        if (scope == Scope::TwoRepairReceive)
        {
            uint32_t current_wh2_expected_overhead = UINT32_MAX;
            const bool valid_expected_overhead = ExpectedTwoRepairOverhead(
                cell.BlockBytes,
                CurrentWh2Arm,
                current_wh2_expected_overhead);
            std::cout << ",candidate_expected_overhead=0"
                      << ",wirehair1_expected_overhead=0"
                      << ",current_wh2_expected_overhead="
                      << (valid_expected_overhead ?
                          current_wh2_expected_overhead : UINT32_MAX)
                      << ",wirehair1_aa_expected_overhead=0";
        }
        std::cout << '\n';
    }
}

void PrintStats(
    const char* level,
    Scope scope,
    uint32_t block_bytes,
    const char* comparison,
    const LogRatioStats& stats,
    uint64_t root = 0u)
{
    std::cout << std::setprecision(12)
              << "RESULT,level=" << level
              << ",scope=" << ScopeName(scope)
              << ",block_bytes=" << block_bytes;
    if (root != 0u) {
        std::cout << ",root=0x" << std::hex << root << std::dec;
    }
    std::cout << ",comparison=" << comparison
              << ",n=" << stats.Count
              << ",gm=" << stats.GeometricMean
              << ",lower95=" << stats.Lower95
              << ",upper95=" << stats.Upper95 << '\n';
}

bool AaGate(const LogRatioStats& stats)
{
    return stats.Lower95 >= 0.98 && stats.Lower95 <= 1.0 &&
        stats.Upper95 >= 1.0 && stats.Upper95 <= 1.02;
}

bool CandidateCellGate(const LogRatioStats& stats)
{
    return stats.GeometricMean < 1.0 && stats.Upper95 < 0.80;
}

#if defined(__linux__)
bool IsExactSingletonCpuSet(const cpu_set_t& set, int target_cpu)
{
    if (target_cpu < 0 || target_cpu >= CPU_SETSIZE ||
        !CPU_ISSET(target_cpu, &set))
    {
        return false;
    }
    for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu)
    {
        if (cpu != target_cpu && CPU_ISSET(cpu, &set)) return false;
    }
    return true;
}
#endif

bool DefaultCpu(int& cpu_out, std::string& diagnostic)
{
#if defined(__linux__)
    cpu_set_t allowed;
    CPU_ZERO(&allowed);
    if (sched_getaffinity(0, sizeof(allowed), &allowed) != 0)
    {
        diagnostic = std::string("sched_getaffinity failed: ") +
            std::strerror(errno);
        return false;
    }
    for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu)
    {
        if (CPU_ISSET(cpu, &allowed)) {
            cpu_out = cpu;
            return true;
        }
    }
    diagnostic = "no CPU is available to the process";
    return false;
#else
    (void)cpu_out;
    diagnostic = "K2 projective timing screen requires Linux affinity";
    return false;
#endif
}

bool PinCpu(int target_cpu, std::string& diagnostic)
{
#if defined(__linux__)
    if (target_cpu < 0 || target_cpu >= CPU_SETSIZE)
    {
        diagnostic = "CPU is outside CPU_SETSIZE";
        return false;
    }
    cpu_set_t requested;
    CPU_ZERO(&requested);
    CPU_SET(target_cpu, &requested);
    if (sched_setaffinity(0, sizeof(requested), &requested) != 0)
    {
        diagnostic = std::string("sched_setaffinity failed: ") +
            std::strerror(errno);
        return false;
    }
    return true;
#else
    (void)target_cpu;
    diagnostic = "K2 projective timing screen requires Linux affinity";
    return false;
#endif
}

bool VerifyCpu(int target_cpu, std::string& diagnostic)
{
#if defined(__linux__)
    if (target_cpu < 0 || target_cpu >= CPU_SETSIZE)
    {
        diagnostic = "CPU is outside CPU_SETSIZE";
        return false;
    }
    cpu_set_t observed;
    CPU_ZERO(&observed);
    if (sched_getaffinity(0, sizeof(observed), &observed) != 0)
    {
        diagnostic = std::string("sched_getaffinity verification failed: ") +
            std::strerror(errno);
        return false;
    }
    if (!IsExactSingletonCpuSet(observed, target_cpu))
    {
        diagnostic = "observed affinity is not the requested singleton";
        return false;
    }
    const int observed_cpu = sched_getcpu();
    if (observed_cpu != target_cpu)
    {
        std::ostringstream message;
        message << "sched_getcpu observed " << observed_cpu
                << " instead of " << target_cpu;
        diagnostic = message.str();
        return false;
    }
    return true;
#else
    (void)target_cpu;
    diagnostic = "K2 projective timing screen requires Linux affinity";
    return false;
#endif
}

bool ParseCpu(const char* text, int& cpu_out)
{
    if (!text || !*text) return false;
    errno = 0;
    char* end = nullptr;
    const long value = std::strtol(text, &end, 10);
    if (errno != 0 || !end || *end != '\0' || value < 0 ||
        value > INT_MAX)
    {
        return false;
    }
    cpu_out = static_cast<int>(value);
    return true;
}

bool CheckDirectionHashes()
{
    static const uint32_t boundary_ids[] = {
        0u, 1u, 2u, 255u, 256u, 257u, 258u, 265u, 270u,
        UINT32_MAX - 1u, UINT32_MAX
    };
    std::string boundary_bytes;
    for (uint32_t id : boundary_ids)
    {
        const K2ProjectiveDirection direction =
            wirehair_wh2_bench::K2ProjectiveDirectionForPacketId(id);
        AppendLe32(boundary_bytes, id);
        AppendLe16(boundary_bytes, direction.Point);
        boundary_bytes.push_back(static_cast<char>(direction.Alpha));
        boundary_bytes.push_back(static_cast<char>(direction.Beta));
    }
    if (wirehair::wh2_benchmark::Sha256Hex(boundary_bytes) !=
        kBoundaryDirectionSha256)
    {
        return false;
    }

    std::string high_bytes;
    high_bytes.reserve(static_cast<size_t>(65536u - 257u + 1u) * 4u);
    for (uint32_t id = 257u; id <= 65536u; ++id)
    {
        const K2ProjectiveDirection direction =
            wirehair_wh2_bench::K2ProjectiveDirectionForPacketId(id);
        AppendLe16(high_bytes, direction.Point);
        high_bytes.push_back(static_cast<char>(direction.Alpha));
        high_bytes.push_back(static_cast<char>(direction.Beta));
    }
    return wirehair::wh2_benchmark::Sha256Hex(high_bytes) ==
        kHighDirectionSha256;
}

bool TestDirectionBoundariesAndAliases()
{
    if (wirehair_wh2_bench::K2ProjectiveHighIdSalt() != 0u ||
        !CheckDirectionHashes())
    {
        return false;
    }
    const K2ProjectiveDirection zero =
        wirehair_wh2_bench::K2ProjectiveDirectionForPacketId(0u);
    const K2ProjectiveDirection one =
        wirehair_wh2_bench::K2ProjectiveDirectionForPacketId(1u);
    const K2ProjectiveDirection two =
        wirehair_wh2_bench::K2ProjectiveDirectionForPacketId(2u);
    const K2ProjectiveDirection last_low =
        wirehair_wh2_bench::K2ProjectiveDirectionForPacketId(256u);
    if (zero.Point != 0u || zero.Alpha != 1u || zero.Beta != 0u ||
        one.Point != 1u || one.Alpha != 0u || one.Beta != 1u ||
        two.Point != 2u || two.Alpha != 1u || two.Beta != 1u ||
        last_low.Point != 256u || last_low.Alpha != 1u ||
        last_low.Beta != 255u)
    {
        return false;
    }
    const uint32_t high_boundaries[] = {
        257u, 258u, UINT32_MAX - 1u, UINT32_MAX
    };
    for (uint32_t id : high_boundaries)
    {
        const K2ProjectiveDirection direction =
            wirehair_wh2_bench::K2ProjectiveDirectionForPacketId(id);
        if (direction.Point < 2u || direction.Point > 256u ||
            direction.Alpha != 1u || direction.Beta == 0u ||
            direction.Point != static_cast<uint16_t>(direction.Beta) + 1u ||
            !SameDirection(
                direction,
                wirehair_wh2_bench::K2ProjectiveDirectionForPacketId(id)))
        {
            return false;
        }
    }

    AliasGroup pair;
    if (!FindFirstHighAliasGroup(2u, pair) || pair.Ids.size() != 2u ||
        pair.Ids[0] != 265u || pair.Ids[1] != 270u ||
        pair.Point != 191u)
    {
        return false;
    }
    return SameDirection(
        wirehair_wh2_bench::K2ProjectiveDirectionForPacketId(pair.Ids[0]),
        wirehair_wh2_bench::K2ProjectiveDirectionForPacketId(pair.Ids[1]));
}

bool TestAllLowDirectionPairs()
{
    static const uint32_t block_bytes = 8u;
    std::vector<uint8_t> source(static_cast<size_t>(kK) * block_bytes);
    for (size_t i = 0u; i < source.size(); ++i) {
        source[i] = static_cast<uint8_t>(i * 29u + 7u);
    }
    K2ProjectiveCodec encoder;
    if (encoder.InitializeEncoder(
            source.data(), source.size(), block_bytes) != Wirehair_Success)
    {
        return false;
    }
    std::vector<std::vector<uint8_t> > packets(257u);
    for (uint32_t id = 0u; id <= 256u; ++id)
    {
        uint32_t bytes = 0u;
        if (!EncodeCandidatePacket(encoder, id, packets[id], bytes) ||
            bytes != block_bytes)
        {
            return false;
        }
    }

    uint32_t pairs = 0u;
    for (uint32_t first = 0u; first <= 256u; ++first)
    {
        for (uint32_t second = first + 1u; second <= 256u; ++second)
        {
            if (!IndependentDirections(first, second)) return false;
            K2ProjectiveCodec decoder;
            if (decoder.InitializeDecoder(
                    source.size(), block_bytes) != Wirehair_Success ||
                decoder.DecodeResult(
                    first,
                    packets[first].data(),
                    static_cast<uint32_t>(packets[first].size())) !=
                    Wirehair_NeedMore ||
                decoder.DecodeResult(
                    second,
                    packets[second].data(),
                    static_cast<uint32_t>(packets[second].size())) !=
                    Wirehair_Success ||
                decoder.Rank() != 2u)
            {
                return false;
            }
            std::vector<uint8_t> recovered(source.size(), uint8_t{0});
            if (decoder.RecoverResult(
                    recovered.data(), recovered.size()) != Wirehair_Success ||
                recovered != source)
            {
                return false;
            }
            ++pairs;
        }
    }
    return pairs == 32896u;
}

bool TestRawIdAndAliasSemantics()
{
    static const uint32_t block_bytes = 32u;
    std::vector<uint8_t> source(static_cast<size_t>(kK) * block_bytes);
    for (size_t i = 0u; i < source.size(); ++i) {
        source[i] = static_cast<uint8_t>(i * 17u + 3u);
    }
    AliasGroup aliases;
    if (!FindFirstHighAliasGroup(3u, aliases) || aliases.Ids.size() != 3u) {
        return false;
    }
    K2ProjectiveCodec encoder;
    if (encoder.InitializeEncoder(
            source.data(), source.size(), block_bytes) != Wirehair_Success)
    {
        return false;
    }
    std::vector<uint8_t> alias_packets[3];
    uint32_t alias_bytes[3] = {};
    for (size_t i = 0u; i < 3u; ++i)
    {
        if (!EncodeCandidatePacket(
                encoder, aliases.Ids[i], alias_packets[i], alias_bytes[i]) ||
            alias_bytes[i] != block_bytes)
        {
            return false;
        }
    }
    if (alias_packets[0] != alias_packets[1] ||
        alias_packets[0] != alias_packets[2])
    {
        return false;
    }
    std::vector<uint8_t> conflict(alias_packets[0]);
    conflict[0] ^= 1u;
    K2ProjectiveCodec decoder;
    if (decoder.InitializeDecoder(
            source.size(), block_bytes) != Wirehair_Success ||
        decoder.DecodeResult(
            aliases.Ids[0], alias_packets[0].data(), block_bytes) !=
            Wirehair_NeedMore ||
        decoder.Rank() != 1u || decoder.AcceptedIdCount() != 1u ||
        decoder.DecodeResult(
            aliases.Ids[0], alias_packets[0].data(), block_bytes) !=
            Wirehair_NeedMore ||
        decoder.AcceptedIdCount() != 1u ||
        decoder.DecodeResult(
            aliases.Ids[0], conflict.data(), block_bytes) !=
            Wirehair_InvalidInput ||
        decoder.Rank() != 1u || decoder.AcceptedIdCount() != 1u ||
        decoder.DecodeResult(
            aliases.Ids[1], alias_packets[1].data(), block_bytes) !=
            Wirehair_NeedMore ||
        decoder.Rank() != 1u || decoder.AcceptedIdCount() != 2u ||
        decoder.DecodeResult(
            aliases.Ids[2], conflict.data(), block_bytes) != Wirehair_Error ||
        decoder.Rank() != 1u || decoder.AcceptedIdCount() != 2u)
    {
        return false;
    }

    std::vector<uint8_t> systematic;
    uint32_t systematic_bytes = 0u;
    if (!EncodeCandidatePacket(
            encoder, 0u, systematic, systematic_bytes) ||
        systematic_bytes != block_bytes ||
        decoder.DecodeResult(
            0u, systematic.data(), systematic_bytes) != Wirehair_Success ||
        decoder.Rank() != 2u || decoder.AcceptedIdCount() != 3u)
    {
        return false;
    }
    std::vector<uint8_t> recovered(source.size(), uint8_t{0});
    if (decoder.RecoverResult(
            recovered.data(), recovered.size()) != Wirehair_Success ||
        recovered != source ||
        decoder.DecodeResult(
            aliases.Ids[0], alias_packets[0].data(), block_bytes) !=
            Wirehair_Success ||
        decoder.DecodeResult(
            aliases.Ids[0], conflict.data(), block_bytes) !=
            Wirehair_InvalidInput ||
        decoder.DecodeResult(
            aliases.Ids[2], conflict.data(), block_bytes) != Wirehair_Error)
    {
        return false;
    }
    return decoder.Rank() == 2u && decoder.AcceptedIdCount() == 3u;
}

bool TestPartialTailPaddingAndRetry()
{
    static const uint32_t block_bytes = 31u;
    static const uint32_t tail_bytes = 7u;
    const size_t message_bytes = block_bytes + tail_bytes;
    std::vector<uint8_t> source(message_bytes);
    for (size_t i = 0u; i < source.size(); ++i) {
        source[i] = static_cast<uint8_t>(i * 11u + 5u);
    }
    AliasGroup aliases;
    if (!FindFirstHighAliasGroup(2u, aliases)) return false;
    K2ProjectiveCodec encoder;
    if (encoder.InitializeEncoder(
            source.data(), source.size(), block_bytes) != Wirehair_Success)
    {
        return false;
    }
    std::vector<uint8_t> tail;
    std::vector<uint8_t> high_first;
    std::vector<uint8_t> high_second;
    std::vector<uint8_t> systematic_zero;
    uint32_t tail_size = 0u;
    uint32_t high_first_size = 0u;
    uint32_t high_second_size = 0u;
    uint32_t systematic_zero_size = 0u;
    if (!EncodeCandidatePacket(encoder, 1u, tail, tail_size) ||
        !EncodeCandidatePacket(
            encoder, aliases.Ids[0], high_first, high_first_size) ||
        !EncodeCandidatePacket(
            encoder, aliases.Ids[1], high_second, high_second_size) ||
        !EncodeCandidatePacket(
            encoder, 0u, systematic_zero, systematic_zero_size) ||
        tail_size != tail_bytes || high_first_size != block_bytes ||
        high_second_size != block_bytes ||
        systematic_zero_size != block_bytes || high_first != high_second)
    {
        return false;
    }

    K2ProjectiveCodec alias_decoder;
    if (alias_decoder.InitializeDecoder(
            source.size(), block_bytes) != Wirehair_Success ||
        alias_decoder.DecodeResult(
            aliases.Ids[0], high_first.data(), high_first_size) !=
            Wirehair_NeedMore ||
        alias_decoder.DecodeResult(
            aliases.Ids[1], high_second.data(), high_second_size) !=
            Wirehair_NeedMore ||
        alias_decoder.Rank() != 1u ||
        alias_decoder.DecodeResult(1u, tail.data(), tail_size) !=
            Wirehair_Success)
    {
        return false;
    }
    std::vector<uint8_t> recovered(source.size(), uint8_t{0});
    if (alias_decoder.RecoverResult(
            recovered.data(), recovered.size()) != Wirehair_Success ||
        recovered != source)
    {
        return false;
    }

    std::vector<uint8_t> corrupt(high_first);
    corrupt[tail_bytes] ^= 1u;
    K2ProjectiveCodec retry_decoder;
    if (retry_decoder.InitializeDecoder(
            source.size(), block_bytes) != Wirehair_Success ||
        retry_decoder.DecodeResult(
            0u, systematic_zero.data(), systematic_zero_size) !=
            Wirehair_NeedMore ||
        retry_decoder.DecodeResult(
            aliases.Ids[0], corrupt.data(), high_first_size) !=
            Wirehair_Error ||
        retry_decoder.Rank() != 1u ||
        retry_decoder.AcceptedIdCount() != 1u ||
        retry_decoder.DecodeResult(
            aliases.Ids[0], high_first.data(), high_first_size) !=
            Wirehair_Success)
    {
        return false;
    }
    std::fill(recovered.begin(), recovered.end(), uint8_t{0});
    return retry_decoder.RecoverResult(
               recovered.data(), recovered.size()) == Wirehair_Success &&
        recovered == source;
}

bool TestUint32MaxEndToEnd()
{
    static const uint32_t block_bytes = 24u;
    std::vector<uint8_t> source(static_cast<size_t>(kK) * block_bytes);
    for (size_t i = 0u; i < source.size(); ++i) {
        source[i] = static_cast<uint8_t>(i * 23u + 13u);
    }
    K2ProjectiveCodec encoder;
    if (encoder.InitializeEncoder(
            source.data(), source.size(), block_bytes) != Wirehair_Success)
    {
        return false;
    }
    std::vector<uint8_t> high_packet;
    std::vector<uint8_t> systematic_packet;
    uint32_t high_bytes = 0u;
    uint32_t systematic_bytes = 0u;
    if (!EncodeCandidatePacket(
            encoder, UINT32_MAX, high_packet, high_bytes) ||
        !EncodeCandidatePacket(
            encoder, 0u, systematic_packet, systematic_bytes) ||
        high_bytes != block_bytes || systematic_bytes != block_bytes ||
        !IndependentDirections(UINT32_MAX, 0u))
    {
        return false;
    }
    std::vector<uint8_t> conflict(high_packet);
    conflict[0] ^= 1u;
    K2ProjectiveCodec decoder;
    if (decoder.InitializeDecoder(
            source.size(), block_bytes) != Wirehair_Success ||
        decoder.DecodeResult(
            UINT32_MAX, high_packet.data(), high_bytes) != Wirehair_NeedMore ||
        decoder.DecodeResult(
            UINT32_MAX, high_packet.data(), high_bytes) != Wirehair_NeedMore ||
        decoder.AcceptedIdCount() != 1u || decoder.Rank() != 1u ||
        decoder.DecodeResult(
            UINT32_MAX, conflict.data(), high_bytes) !=
            Wirehair_InvalidInput ||
        decoder.AcceptedIdCount() != 1u || decoder.Rank() != 1u ||
        decoder.DecodeResult(
            0u, systematic_packet.data(), systematic_bytes) !=
            Wirehair_Success ||
        decoder.AcceptedIdCount() != 2u || decoder.Rank() != 2u)
    {
        return false;
    }
    std::vector<uint8_t> recovered(source.size(), uint8_t{0});
    if (decoder.RecoverResult(
            recovered.data(), recovered.size()) != Wirehair_Success ||
        recovered != source ||
        decoder.DecodeResult(
            UINT32_MAX, high_packet.data(), high_bytes) != Wirehair_Success ||
        decoder.DecodeResult(
            UINT32_MAX, conflict.data(), high_bytes) != Wirehair_InvalidInput)
    {
        return false;
    }
    std::fill(recovered.begin(), recovered.end(), uint8_t{0});
    return decoder.RecoverResult(
               recovered.data(), recovered.size()) == Wirehair_Success &&
        recovered == source;
}

bool TestPostSuccessRawIdLedgerCap()
{
    static const uint32_t block_bytes = 8u;
    std::vector<uint8_t> source(static_cast<size_t>(kK) * block_bytes);
    for (size_t i = 0u; i < source.size(); ++i) {
        source[i] = static_cast<uint8_t>(i * 37u + 9u);
    }
    K2ProjectiveCodec encoder;
    K2ProjectiveCodec decoder;
    if (encoder.InitializeEncoder(
            source.data(), source.size(), block_bytes) != Wirehair_Success ||
        decoder.InitializeDecoder(
            source.size(), block_bytes) != Wirehair_Success)
    {
        return false;
    }

    std::vector<uint8_t> packet;
    uint32_t packet_bytes = 0u;
    if (!EncodeCandidatePacket(encoder, 0u, packet, packet_bytes) ||
        decoder.DecodeResult(0u, packet.data(), packet_bytes) !=
            Wirehair_NeedMore ||
        !EncodeCandidatePacket(encoder, 1u, packet, packet_bytes) ||
        decoder.DecodeResult(1u, packet.data(), packet_bytes) !=
            Wirehair_Success ||
        decoder.AcceptedIdCount() != 2u)
    {
        return false;
    }

    uint32_t next_id = 2u;
    while (decoder.AcceptedIdCount() <
        K2ProjectiveCodec::kMaxAcceptedPacketIds)
    {
        if (!EncodeCandidatePacket(
                encoder, next_id, packet, packet_bytes) ||
            decoder.DecodeResult(
                next_id, packet.data(), packet_bytes) != Wirehair_Success)
        {
            return false;
        }
        ++next_id;
    }
    if (decoder.AcceptedIdCount() !=
        K2ProjectiveCodec::kMaxAcceptedPacketIds)
    {
        return false;
    }

    std::vector<uint8_t> systematic_zero;
    uint32_t systematic_zero_bytes = 0u;
    if (!EncodeCandidatePacket(
            encoder, 0u, systematic_zero, systematic_zero_bytes) ||
        decoder.DecodeResult(
            0u,
            systematic_zero.data(),
            systematic_zero_bytes) != Wirehair_Success)
    {
        return false;
    }
    std::vector<uint8_t> conflict(systematic_zero);
    conflict[0] ^= 1u;
    if (decoder.DecodeResult(
            0u, conflict.data(), systematic_zero_bytes) !=
            Wirehair_InvalidInput ||
        !EncodeCandidatePacket(encoder, next_id, packet, packet_bytes))
    {
        return false;
    }
    std::vector<uint8_t> new_id_conflict(packet);
    new_id_conflict[0] ^= 1u;
    if (decoder.DecodeResult(
            next_id, new_id_conflict.data(), packet_bytes) != Wirehair_Error ||
        decoder.DecodeResult(
            next_id, packet.data(), packet_bytes) !=
            Wirehair_ExtraInsufficient ||
        decoder.AcceptedIdCount() !=
            K2ProjectiveCodec::kMaxAcceptedPacketIds)
    {
        return false;
    }
    std::vector<uint8_t> recovered(source.size(), uint8_t{0});
    return decoder.RecoverResult(
               recovered.data(), recovered.size()) == Wirehair_Success &&
        recovered == source;
}

#if defined(WIREHAIR_WH2_K2_PROJECTIVE_TEST_HOOKS)
bool TestAllocationTransactionality()
{
    static const uint32_t block_bytes = 16u;
    std::vector<uint8_t> source(static_cast<size_t>(kK) * block_bytes);
    std::vector<uint8_t> replacement(source.size());
    for (size_t i = 0u; i < source.size(); ++i)
    {
        source[i] = static_cast<uint8_t>(i + 1u);
        replacement[i] = static_cast<uint8_t>(255u - i);
    }
    K2ProjectiveCodec encoder;
    if (encoder.InitializeEncoder(
            source.data(), source.size(), block_bytes) != Wirehair_Success)
    {
        return false;
    }
    std::vector<uint8_t> expected;
    uint32_t expected_bytes = 0u;
    if (!EncodeCandidatePacket(
            encoder, 2u, expected, expected_bytes)) return false;
    wirehair_wh2_bench::
        SetK2ProjectiveAllocationFailureCountdownForTesting(0);
    const WirehairResult encoder_oom = encoder.InitializeEncoder(
        replacement.data(), replacement.size(), block_bytes);
    wirehair_wh2_bench::
        SetK2ProjectiveAllocationFailureCountdownForTesting(-1);
    std::vector<uint8_t> retained;
    uint32_t retained_bytes = 0u;
    if (encoder_oom != Wirehair_OOM ||
        !EncodeCandidatePacket(encoder, 2u, retained, retained_bytes) ||
        retained_bytes != expected_bytes || retained != expected)
    {
        return false;
    }

    K2ProjectiveCodec decoder;
    if (decoder.InitializeDecoder(
            source.size(), block_bytes) != Wirehair_Success ||
        decoder.DecodeResult(
            2u, expected.data(), expected_bytes) != Wirehair_NeedMore)
    {
        return false;
    }
    for (int64_t countdown = 0; countdown <= 1; ++countdown)
    {
        wirehair_wh2_bench::
            SetK2ProjectiveAllocationFailureCountdownForTesting(countdown);
        const WirehairResult decoder_oom = decoder.InitializeDecoder(
            replacement.size(), block_bytes);
        wirehair_wh2_bench::
            SetK2ProjectiveAllocationFailureCountdownForTesting(-1);
        if (decoder_oom != Wirehair_OOM || decoder.Rank() != 1u ||
            decoder.AcceptedIdCount() != 1u)
        {
            return false;
        }
    }
    std::vector<uint8_t> systematic;
    uint32_t systematic_bytes = 0u;
    if (!EncodeCandidatePacket(
            encoder, 0u, systematic, systematic_bytes) ||
        decoder.DecodeResult(
            0u, systematic.data(), systematic_bytes) != Wirehair_Success)
    {
        return false;
    }
    std::vector<uint8_t> recovered(source.size(), uint8_t{0});
    return decoder.RecoverResult(
               recovered.data(), recovered.size()) == Wirehair_Success &&
        recovered == source;
}
#endif

bool PrefixHasProjectiveRankTwo(
    const std::vector<uint32_t>& ids,
    size_t count)
{
    if (count < 2u || count > ids.size()) return false;
    const uint16_t first_point =
        wirehair_wh2_bench::K2ProjectiveDirectionForPacketId(ids[0]).Point;
    for (size_t i = 1u; i < count; ++i)
    {
        if (wirehair_wh2_bench::K2ProjectiveDirectionForPacketId(ids[i]).
                Point != first_point)
        {
            return true;
        }
    }
    return false;
}

struct FrozenK2CorpusEntry
{
    wirehair::wh2_benchmark::FrozenRecoveryCell Cell;
    std::vector<uint32_t> DeliveredIds;
};

bool BuildValidatedFrozenK2Corpus(
    std::vector<FrozenK2CorpusEntry>& corpus_out,
    std::string& diagnostic)
{
    corpus_out.clear();
    diagnostic.clear();
    if (wirehair::wh2_benchmark::DevelopmentRecoveryDomainSha256() !=
        kDevelopmentDomainSha256)
    {
        diagnostic = "development recovery domain hash mismatch";
        return false;
    }
    const std::vector<wirehair::wh2_benchmark::FrozenRecoveryCell> cells =
        wirehair::wh2_benchmark::EnumerateDevelopmentRecoveryCells();
    std::string canonical_stream;
    std::string trace_stream;
    uint32_t k2_cell_count = 0u;
    for (const wirehair::wh2_benchmark::FrozenRecoveryCell& cell : cells)
    {
        if (cell.K != kK) continue;
        ++k2_cell_count;
        canonical_stream +=
            wirehair::wh2_benchmark::CanonicalRecoveryCellJson(cell);
        canonical_stream.push_back('\n');
        wirehair::wh2_benchmark::FrozenPacketTrace trace;
        if (wirehair::wh2_benchmark::GenerateFrozenPacketTrace(cell, trace) !=
                wirehair::wh2_benchmark::FrozenTraceStatus::Complete ||
            trace.delivered_ids.size() != 6u)
        {
            diagnostic = "frozen K2 trace generation/cardinality mismatch";
            return false;
        }
        AppendLe64(trace_stream, static_cast<uint64_t>(cell.ordinal));
        for (uint32_t id : trace.delivered_ids) AppendLe32(trace_stream, id);

        FrozenK2CorpusEntry entry;
        entry.Cell = cell;
        entry.DeliveredIds.swap(trace.delivered_ids);
        corpus_out.push_back(std::move(entry));
    }
    if (k2_cell_count != 12u || corpus_out.size() != 12u)
    {
        diagnostic = "frozen K2 cell cardinality mismatch";
        return false;
    }
    if (wirehair::wh2_benchmark::Sha256Hex(canonical_stream) !=
            kK2CanonicalCellStreamSha256)
    {
        diagnostic = "frozen K2 canonical stream hash mismatch";
        return false;
    }
    if (wirehair::wh2_benchmark::Sha256Hex(trace_stream) !=
            kK2TraceStreamSha256)
    {
        diagnostic = "frozen K2 trace stream hash mismatch";
        return false;
    }
    return true;
}

HarnessOutcome CheckCandidateFrozenK2Recovery(
    const std::vector<FrozenK2CorpusEntry>& corpus,
    std::string& diagnostic)
{
    diagnostic.clear();
    if (corpus.size() != 12u)
    {
        diagnostic = "validated frozen K2 corpus cardinality changed";
        return HarnessOutcome::Invalid;
    }
    for (const FrozenK2CorpusEntry& entry : corpus)
    {
        const wirehair::wh2_benchmark::FrozenRecoveryCell& cell = entry.Cell;
        const std::vector<uint32_t>& delivered_ids = entry.DeliveredIds;
        bool seen_points[257] = {};
        for (uint32_t id : delivered_ids)
        {
            const uint16_t point =
                wirehair_wh2_bench::K2ProjectiveDirectionForPacketId(id).
                    Point;
            if (point > 256u || seen_points[point])
            {
                diagnostic =
                    "candidate frozen K2 trace has a repeated direction";
                return HarnessOutcome::CandidateReject;
            }
            seen_points[point] = true;
        }

        std::vector<uint8_t> source;
        if (!wirehair_wh2_bench::MakeDeterministicSource(
                kK,
                cell.block_bytes,
                cell.loss_seed ^ UINT64_C(0x6a09e667f3bcc909),
                source))
        {
            diagnostic = "frozen K2 deterministic source generation failed";
            return HarnessOutcome::Invalid;
        }
        K2ProjectiveCodec encoder;
        const WirehairResult encoder_result = encoder.InitializeEncoder(
            source.data(), source.size(), cell.block_bytes);
        if (encoder_result != Wirehair_Success)
        {
            diagnostic = "candidate frozen K2 encoder initialization failed";
            return encoder_result == Wirehair_OOM ||
                    encoder_result == Wirehair_UnsupportedPlatform ?
                HarnessOutcome::Invalid : HarnessOutcome::CandidateReject;
        }
        std::vector<std::vector<uint8_t> > packets(6u);
        std::vector<uint32_t> packet_bytes(6u, 0u);
        for (size_t i = 0u; i < 6u; ++i)
        {
            if (!EncodeCandidatePacket(
                    encoder,
                    delivered_ids[i],
                    packets[i],
                    packet_bytes[i]))
            {
                diagnostic = "candidate frozen K2 packet encoding failed";
                return HarnessOutcome::CandidateReject;
            }
        }
        for (uint32_t overhead = 0u; overhead <= 4u; ++overhead)
        {
            const size_t prefix = static_cast<size_t>(kK + overhead);
            if (!PrefixHasProjectiveRankTwo(delivered_ids, prefix))
            {
                diagnostic = "candidate frozen K2 prefix lacks rank two";
                return HarnessOutcome::CandidateReject;
            }
            K2ProjectiveCodec decoder;
            const WirehairResult decoder_result = decoder.InitializeDecoder(
                source.size(), cell.block_bytes);
            if (decoder_result != Wirehair_Success)
            {
                diagnostic =
                    "candidate frozen K2 decoder initialization failed";
                return decoder_result == Wirehair_OOM ||
                        decoder_result == Wirehair_UnsupportedPlatform ?
                    HarnessOutcome::Invalid : HarnessOutcome::CandidateReject;
            }
            WirehairResult result = Wirehair_NeedMore;
            for (size_t i = 0u; i < prefix; ++i)
            {
                result = decoder.DecodeResult(
                    delivered_ids[i],
                    packets[i].data(),
                    packet_bytes[i]);
                if ((i == 0u && result != Wirehair_NeedMore) ||
                    (i > 0u && result != Wirehair_Success))
                {
                    diagnostic = "candidate frozen K2 decode failed";
                    return result == Wirehair_OOM ||
                            result == Wirehair_UnsupportedPlatform ?
                        HarnessOutcome::Invalid :
                        HarnessOutcome::CandidateReject;
                }
            }
            std::vector<uint8_t> recovered(source.size(), uint8_t{0});
            const WirehairResult recover_result = decoder.RecoverResult(
                recovered.data(), recovered.size());
            if (!decoder.IsDecoded() || decoder.Rank() != 2u ||
                recover_result != Wirehair_Success || recovered != source)
            {
                diagnostic = "candidate frozen K2 recovery mismatch";
                return recover_result == Wirehair_OOM ||
                        recover_result == Wirehair_UnsupportedPlatform ?
                    HarnessOutcome::Invalid : HarnessOutcome::CandidateReject;
            }
        }
    }
    return HarnessOutcome::Complete;
}

bool TestFrozenK2Recovery()
{
    std::vector<FrozenK2CorpusEntry> corpus;
    std::string diagnostic;
    return BuildValidatedFrozenK2Corpus(corpus, diagnostic) &&
        CheckCandidateFrozenK2Recovery(corpus, diagnostic) ==
            HarnessOutcome::Complete;
}

bool TestStatisticsAndCounterbalance()
{
    std::vector<double> identity_logs(kObservationCount, 0.0);
    LogRatioStats identity;
    if (!ComputeLogRatioStats(identity_logs, identity) ||
        identity.Count != kObservationCount ||
        identity.GeometricMean != 1.0 || identity.Lower95 != 1.0 ||
        identity.Upper95 != 1.0 || !AaGate(identity))
    {
        return false;
    }
    for (uint32_t row = 0u; row < 4u; ++row)
    {
        bool seen[ArmCount] = {};
        for (uint32_t position = 0u; position < 4u; ++position)
        {
            const uint32_t arm =
                static_cast<uint32_t>(kCounterbalancedOrders[row][position]);
            if (arm >= static_cast<uint32_t>(ArmCount) || seen[arm]) {
                return false;
            }
            seen[arm] = true;
        }
    }
    for (uint32_t position = 0u; position < 4u; ++position)
    {
        bool seen[ArmCount] = {};
        for (uint32_t row = 0u; row < 4u; ++row)
        {
            const uint32_t arm =
                static_cast<uint32_t>(kCounterbalancedOrders[row][position]);
            if (seen[arm]) return false;
            seen[arm] = true;
        }
    }
    return true;
}

bool TestOutcomeClassification()
{
    std::string diagnostic;
    TimedArmResult encoder_success;
    encoder_success.Result = Wirehair_Success;
    encoder_success.ElapsedNanoseconds = 1u;
    encoder_success.BytesVerified = true;
    if (ClassifyInvocation(
            64u,
            Scope::Encoder,
            CandidateArm,
            encoder_success,
            true,
            "selftest",
            diagnostic) != HarnessOutcome::Complete)
    {
        return false;
    }
    encoder_success.ElapsedNanoseconds = 0u;
    if (ClassifyInvocation(
            64u,
            Scope::Encoder,
            CandidateArm,
            encoder_success,
            true,
            "selftest",
            diagnostic) != HarnessOutcome::Invalid ||
        ClassifyInvocation(
            64u,
            Scope::Encoder,
            CandidateArm,
            encoder_success,
            false,
            "selftest",
            diagnostic) != HarnessOutcome::Complete)
    {
        return false;
    }
    encoder_success.ElapsedNanoseconds = 1u;
    if (ClassifyInvocation(
            64u,
            Scope::Encoder,
            CandidateArm,
            encoder_success,
            false,
            "selftest",
            diagnostic) != HarnessOutcome::Invalid)
    {
        return false;
    }

    TimedArmResult candidate_failure;
    candidate_failure.Result = Wirehair_Error;
    candidate_failure.ElapsedNanoseconds = 1u;
    if (ClassifyInvocation(
            64u,
            Scope::Encoder,
            CandidateArm,
            candidate_failure,
            true,
            "selftest",
            diagnostic) != HarnessOutcome::CandidateReject)
    {
        return false;
    }
    candidate_failure.Result = Wirehair_OOM;
    if (ClassifyInvocation(
            64u,
            Scope::Encoder,
            CandidateArm,
            candidate_failure,
            true,
            "selftest",
            diagnostic) != HarnessOutcome::Invalid)
    {
        return false;
    }
    candidate_failure.Result = Wirehair_UnsupportedPlatform;
    if (ClassifyInvocation(
            64u,
            Scope::Encoder,
            CandidateArm,
            candidate_failure,
            true,
            "selftest",
            diagnostic) != HarnessOutcome::Invalid)
    {
        return false;
    }
    candidate_failure.Result = Wirehair_Error;
    if (ClassifyInvocation(
            64u,
            Scope::Encoder,
            Wirehair1Arm,
            candidate_failure,
            true,
            "selftest",
            diagnostic) != HarnessOutcome::Invalid)
    {
        return false;
    }

    TimedArmResult receive_success;
    receive_success.Result = Wirehair_Success;
    receive_success.ElapsedNanoseconds = 1u;
    receive_success.BytesVerified = true;
    receive_success.DecodedOverhead = 0u;
    if (ClassifyInvocation(
            64u,
            Scope::TwoRepairReceive,
            Wirehair1AaArm,
            receive_success,
            true,
            "selftest",
            diagnostic) != HarnessOutcome::Complete)
    {
        return false;
    }
    receive_success.DecodedOverhead = 1u;
    if (ClassifyInvocation(
            64u,
            Scope::TwoRepairReceive,
            CandidateArm,
            receive_success,
            true,
            "selftest",
            diagnostic) != HarnessOutcome::CandidateReject ||
        ClassifyInvocation(
            64u,
            Scope::TwoRepairReceive,
            Wirehair1AaArm,
            receive_success,
            true,
            "selftest",
            diagnostic) != HarnessOutcome::Invalid ||
        ClassifyInvocation(
            1280u,
            Scope::TwoRepairReceive,
            CurrentWh2Arm,
            receive_success,
            true,
            "selftest",
            diagnostic) != HarnessOutcome::Complete ||
        ClassifyInvocation(
            64u,
            Scope::TwoRepairReceive,
            CurrentWh2Arm,
            receive_success,
            true,
            "selftest",
            diagnostic) != HarnessOutcome::Invalid ||
        ClassifyInvocation(
            64u,
            Scope::IidReceive,
            CandidateArm,
            receive_success,
            true,
            "selftest",
            diagnostic) != HarnessOutcome::Complete)
    {
        return false;
    }
    receive_success.DecodedOverhead = UINT32_MAX;
    if (ClassifyInvocation(
            64u,
            Scope::IidReceive,
            CandidateArm,
            receive_success,
            true,
            "selftest",
            diagnostic) != HarnessOutcome::CandidateReject)
    {
        return false;
    }

    std::vector<FrozenK2CorpusEntry> corpus;
    if (CheckCandidateFrozenK2Recovery(corpus, diagnostic) !=
            HarnessOutcome::Invalid ||
        !BuildValidatedFrozenK2Corpus(corpus, diagnostic) || corpus.empty())
    {
        return false;
    }
    corpus[0].DeliveredIds[1] = corpus[0].DeliveredIds[0];
    return CheckCandidateFrozenK2Recovery(corpus, diagnostic) ==
        HarnessOutcome::CandidateReject;
}

bool TestClockFreeFixturePreflight()
{
    CandidateEncoderFixture empty_candidate_encoder;
    CandidateReceiveFixture empty_candidate_receive;
    NativeEncoderFixture empty_native_encoder;
    NativeReceiveFixture empty_native_receive;
    const TimedArmResult empty_results[] = {
        empty_candidate_encoder.Preflight(),
        empty_candidate_receive.Preflight(),
        empty_native_encoder.Preflight(),
        empty_native_receive.Preflight()
    };
    for (const TimedArmResult& result : empty_results)
    {
        if (result.Result == Wirehair_Success ||
            result.ElapsedNanoseconds != 0u || result.BytesVerified)
        {
            return false;
        }
    }

    std::string diagnostic;
    for (uint32_t block_bytes : kBlockBytes)
    {
        uint32_t expected_attempt = UINT32_MAX;
        if (!CurrentWh2SelectedAttempt(block_bytes, expected_attempt)) {
            return false;
        }
        for (uint64_t root : kTrainingRoots)
        {
            RootCell cell;
            if (BuildRootCell(block_bytes, root, cell, diagnostic) !=
                    HarnessOutcome::Complete ||
                cell.CurrentWh2Attempt != expected_attempt ||
                PreflightRootCell(cell, diagnostic) !=
                    HarnessOutcome::Complete)
            {
                std::cerr << "clock-free fixture validation failed for B="
                          << block_bytes << " root=0x" << std::hex << root
                          << std::dec << ": "
                          << diagnostic << '\n';
                return false;
            }
        }
    }
    return true;
}

bool RunSelfTest()
{
    const WirehairResult init_result = wirehair_init();
    if (init_result != Wirehair_Success)
    {
        std::cerr << "K2 projective selftest global init failed\n";
        return false;
    }
    struct Test
    {
        const char* Name;
        bool (*Run)();
    };
    const Test tests[] = {
        { "direction-boundaries-aliases", TestDirectionBoundariesAndAliases },
        { "all-low-direction-pairs", TestAllLowDirectionPairs },
        { "raw-id-alias-semantics", TestRawIdAndAliasSemantics },
        { "partial-tail-padding-retry", TestPartialTailPaddingAndRetry },
        { "uint32-max-end-to-end", TestUint32MaxEndToEnd },
        { "post-success-raw-id-ledger-cap", TestPostSuccessRawIdLedgerCap },
#if defined(WIREHAIR_WH2_K2_PROJECTIVE_TEST_HOOKS)
        { "allocation-transactionality", TestAllocationTransactionality },
#endif
        { "frozen-k2-recovery", TestFrozenK2Recovery },
        { "statistics-counterbalance", TestStatisticsAndCounterbalance },
        { "outcome-classification", TestOutcomeClassification },
        { "clock-free-fixture-preflight", TestClockFreeFixturePreflight }
    };
    for (const Test& test : tests)
    {
        if (!test.Run())
        {
            std::cerr << "K2 projective selftest failed: "
                      << test.Name << '\n';
            return false;
        }
        std::cout << "SELFTEST,name=" << test.Name << ",status=PASS\n";
    }
    std::cout << "K2_PROJECTIVE_SELFTEST=PASS\n";
    return true;
}

int RunScreen(int requested_cpu, bool cpu_was_explicit)
{
    std::string diagnostic;
    int target_cpu = requested_cpu;
    if (!cpu_was_explicit && !DefaultCpu(target_cpu, diagnostic))
    {
        std::cerr << "K2 projective screen CPU selection failed: "
                  << diagnostic << '\n';
        return 1;
    }
    if (wirehair_init() != Wirehair_Success)
    {
        std::cerr << "K2 projective screen global initialization failed\n";
        return 1;
    }
    std::vector<FrozenK2CorpusEntry> frozen_k2_corpus;
    if (!BuildValidatedFrozenK2Corpus(frozen_k2_corpus, diagnostic))
    {
        std::cerr << "K2 projective frozen-control preflight invalid: "
                  << diagnostic << '\n';
        std::cout << "K2_PROJECTIVE_SCREEN=INVALID\n";
        return 1;
    }
    if (!TestDirectionBoundariesAndAliases())
    {
        std::cerr << "K2 projective screen mapping preflight failed\n";
        std::cout << "K2_PROJECTIVE_SCREEN=REJECT\n";
        return 2;
    }
    const HarnessOutcome frozen_recovery_outcome =
        CheckCandidateFrozenK2Recovery(frozen_k2_corpus, diagnostic);
    if (frozen_recovery_outcome != HarnessOutcome::Complete)
    {
        std::cerr << "K2 projective candidate frozen-recovery preflight "
                  << "failed: " << diagnostic << '\n';
        if (frozen_recovery_outcome == HarnessOutcome::CandidateReject)
        {
            std::cout << "K2_PROJECTIVE_SCREEN=REJECT\n";
            return 2;
        }
        std::cout << "K2_PROJECTIVE_SCREEN=INVALID\n";
        return 1;
    }

    std::vector<std::unique_ptr<RootCell> > cells;
    for (uint32_t block_bytes : kBlockBytes)
    {
        for (uint64_t root : kTrainingRoots)
        {
            std::unique_ptr<RootCell> cell(new RootCell);
            const HarnessOutcome build_outcome =
                BuildRootCell(block_bytes, root, *cell, diagnostic);
            if (build_outcome != HarnessOutcome::Complete)
            {
                std::cerr << "K2 projective screen construction failed for B="
                          << block_bytes << " root=0x" << std::hex << root
                          << std::dec << ": " << diagnostic << '\n';
                if (build_outcome == HarnessOutcome::CandidateReject) {
                    std::cout << "K2_PROJECTIVE_SCREEN=REJECT\n";
                    return 2;
                }
                std::cout << "K2_PROJECTIVE_SCREEN=INVALID\n";
                return 1;
            }
            cells.push_back(std::move(cell));
        }
    }
    for (size_t cell_index = 0u; cell_index < cells.size(); ++cell_index)
    {
        const HarnessOutcome outcome =
            PreflightRootCell(*cells[cell_index], diagnostic);
        if (outcome != HarnessOutcome::Complete)
        {
            std::cerr << "K2 projective screen clock-free preflight failed "
                      << "for B=" << cells[cell_index]->BlockBytes
                      << " root=0x" << std::hex << cells[cell_index]->Root
                      << std::dec << ": " << diagnostic << '\n';
            if (outcome == HarnessOutcome::CandidateReject) {
                std::cout << "K2_PROJECTIVE_SCREEN=REJECT\n";
                return 2;
            }
            std::cout << "K2_PROJECTIVE_SCREEN=INVALID\n";
            return 1;
        }
    }

    std::vector<PanelMeasurement> panel;
    const size_t expected_panel_size =
        (sizeof(kBlockBytes) / sizeof(kBlockBytes[0])) * 3u *
        (sizeof(kTrainingRoots) / sizeof(kTrainingRoots[0]));
    panel.reserve(expected_panel_size);
    for (size_t width_index = 0u;
         width_index < sizeof(kBlockBytes) / sizeof(kBlockBytes[0]);
         ++width_index)
    {
        for (uint32_t scope_value = 0u; scope_value < 3u; ++scope_value)
        {
            for (size_t root_index = 0u;
                 root_index <
                    sizeof(kTrainingRoots) / sizeof(kTrainingRoots[0]);
                 ++root_index)
            {
                PanelMeasurement measurement;
                measurement.WidthIndex = width_index;
                measurement.RootIndex = root_index;
                measurement.ScopeValue = static_cast<Scope>(scope_value);
                ReserveMeasurementLogs(measurement.Logs);
                panel.push_back(std::move(measurement));
            }
        }
    }
    if (panel.size() != expected_panel_size)
    {
        std::cerr << "K2 projective screen panel preparation invalid\n";
        std::cout << "K2_PROJECTIVE_SCREEN=INVALID\n";
        return 1;
    }

    std::cout << "CONFIG,K=2,cpu=" << target_cpu
              << ",loss_ppm=" << kLossPpm
              << ",delivered_overhead=" << kReceiveOverheadCap
              << ",observations=" << kObservationCount
              << ",invocations_per_observation="
              << kInvocationsPerObservation
              << ",warmup_invocations=" << kWarmupInvocations
              << ",high_id_salt="
              << wirehair_wh2_bench::K2ProjectiveHighIdSalt()
              << ",current_wh2_selected_attempt_b64=7"
              << ",current_wh2_selected_attempt_b1280=1"
              << ",current_wh2_diagnostic_only=1"
              << ",two_repair_raw_id_start=2"
              << ",two_repair_candidate_expected_overhead=0"
              << ",two_repair_wirehair1_expected_overhead=0"
              << ",two_repair_wirehair1_aa_expected_overhead=0"
              << ",two_repair_current_wh2_b64_expected_overhead=0"
              << ",two_repair_current_wh2_b1280_expected_overhead=1\n";

    if (!PinCpu(target_cpu, diagnostic))
    {
        std::cerr << "K2 projective screen affinity pin failed: "
                  << diagnostic << '\n';
        std::cout << "K2_PROJECTIVE_SCREEN=INVALID\n";
        return 1;
    }
    if (!VerifyCpu(target_cpu, diagnostic))
    {
        std::cerr << "K2 projective screen pre-panel affinity failed: "
                  << diagnostic << '\n';
        std::cout << "K2_PROJECTIVE_SCREEN=INVALID\n";
        return 1;
    }

    HarnessOutcome panel_outcome = HarnessOutcome::Complete;
    for (PanelMeasurement& measurement : panel)
    {
        const size_t cell_index = measurement.WidthIndex * 3u +
            measurement.RootIndex;
        const uint32_t order_offset = static_cast<uint32_t>(
            measurement.WidthIndex * 7u + measurement.RootIndex * 5u +
            static_cast<uint32_t>(measurement.ScopeValue) * 3u);
        panel_outcome = MeasureScope(
            *cells[cell_index],
            measurement.ScopeValue,
            order_offset,
            measurement.Logs,
            diagnostic);
        if (panel_outcome != HarnessOutcome::Complete) break;
    }

    const bool post_panel_affinity_valid = VerifyCpu(target_cpu, diagnostic);
    if (!post_panel_affinity_valid)
    {
        std::cerr << "K2 projective screen post-panel affinity failed: "
                  << diagnostic << '\n';
        std::cout << "K2_PROJECTIVE_SCREEN=INVALID\n";
        return 1;
    }
    if (panel_outcome != HarnessOutcome::Complete)
    {
        std::cerr << "K2 projective screen panel stopped: "
                  << diagnostic << '\n';
        if (panel_outcome == HarnessOutcome::CandidateReject) {
            std::cout << "K2_PROJECTIVE_SCREEN=REJECT\n";
            return 2;
        }
        std::cout << "K2_PROJECTIVE_SCREEN=INVALID\n";
        return 1;
    }

    if (panel.size() != expected_panel_size)
    {
        std::cerr << "K2 projective screen panel cardinality invalid\n";
        std::cout << "K2_PROJECTIVE_SCREEN=INVALID\n";
        return 1;
    }
    for (PanelMeasurement& measurement : panel)
    {
        if (!FinalizeMeasurementLogs(measurement.Logs, diagnostic))
        {
            std::cerr << "K2 projective screen log finalization failed: "
                      << diagnostic << '\n';
            std::cout << "K2_PROJECTIVE_SCREEN=INVALID\n";
            return 1;
        }
    }

    bool candidate_pass = true;
    bool controls_valid = true;
    size_t panel_index = 0u;
    MeasurementLogs pooled;
    for (size_t width_index = 0u;
         width_index < sizeof(kBlockBytes) / sizeof(kBlockBytes[0]);
         ++width_index)
    {
        const uint32_t block_bytes = kBlockBytes[width_index];
        for (uint32_t scope_value = 0u; scope_value < 3u; ++scope_value)
        {
            const Scope scope = static_cast<Scope>(scope_value);
            MeasurementLogs cell_logs;
            for (size_t root_index = 0u;
                 root_index <
                    sizeof(kTrainingRoots) / sizeof(kTrainingRoots[0]);
                 ++root_index)
            {
                const PanelMeasurement& measurement = panel[panel_index++];
                if (measurement.WidthIndex != width_index ||
                    measurement.RootIndex != root_index ||
                    measurement.ScopeValue != scope)
                {
                    std::cerr <<
                        "K2 projective screen panel ordering invalid\n";
                    return 1;
                }
                const size_t cell_index = width_index * 3u + root_index;
                PrintObservations(
                    *cells[cell_index], scope, measurement.Logs);
                LogRatioStats candidate_wh1;
                LogRatioStats candidate_wh2;
                LogRatioStats aa;
                if (!ComputeLogRatioStats(
                        measurement.Logs.CandidateVsWirehair1,
                        candidate_wh1) ||
                    !ComputeLogRatioStats(
                        measurement.Logs.CandidateVsCurrentWh2,
                        candidate_wh2) ||
                    !ComputeLogRatioStats(
                        measurement.Logs.Wirehair1Aa, aa))
                {
                    std::cerr << "K2 projective root statistics invalid\n";
                    return 1;
                }
                PrintStats(
                    "root", scope, block_bytes, "candidate/wirehair1",
                    candidate_wh1, kTrainingRoots[root_index]);
                PrintStats(
                    "root", scope, block_bytes,
                    kCurrentWh2Comparison,
                    candidate_wh2, kTrainingRoots[root_index]);
                PrintStats(
                    "root", scope, block_bytes, "wirehair1/wirehair1-aa",
                    aa, kTrainingRoots[root_index]);
                AppendLogs(
                    measurement.Logs.CandidateVsWirehair1,
                    cell_logs.CandidateVsWirehair1);
                AppendLogs(
                    measurement.Logs.CandidateVsCurrentWh2,
                    cell_logs.CandidateVsCurrentWh2);
                AppendLogs(
                    measurement.Logs.Wirehair1Aa,
                    cell_logs.Wirehair1Aa);
            }
            LogRatioStats candidate_wh1;
            LogRatioStats candidate_wh2;
            LogRatioStats aa;
            if (!ComputeLogRatioStats(
                    cell_logs.CandidateVsWirehair1, candidate_wh1) ||
                !ComputeLogRatioStats(
                    cell_logs.CandidateVsCurrentWh2, candidate_wh2) ||
                !ComputeLogRatioStats(cell_logs.Wirehair1Aa, aa))
            {
                std::cerr << "K2 projective cell statistics invalid\n";
                return 1;
            }
            PrintStats(
                "cell", scope, block_bytes, "candidate/wirehair1",
                candidate_wh1);
            PrintStats(
                "cell", scope, block_bytes,
                kCurrentWh2Comparison,
                candidate_wh2);
            PrintStats(
                "cell", scope, block_bytes, "wirehair1/wirehair1-aa", aa);
            candidate_pass = candidate_pass &&
                CandidateCellGate(candidate_wh1);
            controls_valid = controls_valid && AaGate(aa);
            AppendLogs(
                cell_logs.CandidateVsWirehair1,
                pooled.CandidateVsWirehair1);
            AppendLogs(
                cell_logs.CandidateVsCurrentWh2,
                pooled.CandidateVsCurrentWh2);
            AppendLogs(cell_logs.Wirehair1Aa, pooled.Wirehair1Aa);
        }
    }

    LogRatioStats pooled_candidate_wh1;
    LogRatioStats pooled_candidate_wh2;
    LogRatioStats pooled_aa;
    if (!ComputeLogRatioStats(
            pooled.CandidateVsWirehair1, pooled_candidate_wh1) ||
        !ComputeLogRatioStats(
            pooled.CandidateVsCurrentWh2, pooled_candidate_wh2) ||
        !ComputeLogRatioStats(pooled.Wirehair1Aa, pooled_aa))
    {
        std::cerr << "K2 projective pooled statistics invalid\n";
        return 1;
    }
    std::cout << std::setprecision(12)
              << "RESULT,level=pooled,comparison=candidate/wirehair1"
              << ",n=" << pooled_candidate_wh1.Count
              << ",gm=" << pooled_candidate_wh1.GeometricMean
              << ",lower95=" << pooled_candidate_wh1.Lower95
              << ",upper95=" << pooled_candidate_wh1.Upper95 << '\n'
              << "RESULT,level=pooled,comparison="
              << kCurrentWh2Comparison
              << ",n=" << pooled_candidate_wh2.Count
              << ",gm=" << pooled_candidate_wh2.GeometricMean
              << ",lower95=" << pooled_candidate_wh2.Lower95
              << ",upper95=" << pooled_candidate_wh2.Upper95 << '\n'
              << "RESULT,level=pooled,comparison=wirehair1/wirehair1-aa"
              << ",n=" << pooled_aa.Count
              << ",gm=" << pooled_aa.GeometricMean
              << ",lower95=" << pooled_aa.Lower95
              << ",upper95=" << pooled_aa.Upper95 << '\n';
    candidate_pass = candidate_pass &&
        pooled_candidate_wh1.GeometricMean < 1.0 &&
        pooled_candidate_wh1.Upper95 < 0.70;
    controls_valid = controls_valid && AaGate(pooled_aa);
    if (!controls_valid)
    {
        std::cout << "K2_PROJECTIVE_SCREEN=INVALID\n";
        return 1;
    }
    std::cout << "K2_PROJECTIVE_SCREEN="
              << (candidate_pass ? "PASS" : "REJECT") << '\n';
    return candidate_pass ? 0 : 2;
}

void PrintUsage(const char* program)
{
    std::cerr << "usage: " << program
              << " --selftest | --screen [--cpu N]\n";
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
#if defined(WIREHAIR_WH2_K2_PROJECTIVE_TEST_HOOKS)
        if (argc >= 2 && std::strcmp(argv[1], "--screen") == 0)
        {
            std::cerr <<
                "hook-enabled K2 projective selftest refuses --screen\n";
            return 1;
        }
#endif
        if (argc == 2 && std::strcmp(argv[1], "--selftest") == 0) {
            return RunSelfTest() ? 0 : 1;
        }
        if (argc >= 2 && std::strcmp(argv[1], "--screen") == 0)
        {
            int requested_cpu = -1;
            bool cpu_was_explicit = false;
            if (argc == 4 && std::strcmp(argv[2], "--cpu") == 0 &&
                ParseCpu(argv[3], requested_cpu))
            {
                cpu_was_explicit = true;
            }
            else if (argc != 2)
            {
                PrintUsage(argv[0]);
                return 1;
            }
            return RunScreen(requested_cpu, cpu_was_explicit);
        }
        PrintUsage(argv[0]);
        return 1;
    }
    catch (const std::bad_alloc&) {
        std::cerr << "K2 projective harness allocation failed\n";
        return 1;
    }
    catch (const std::length_error&) {
        std::cerr << "K2 projective harness allocation length failed\n";
        return 1;
    }
    catch (...) {
        std::cerr << "K2 projective harness failed unexpectedly\n";
        return 1;
    }
}
