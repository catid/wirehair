#include "Wh2FrozenTrace.h"
#include "Wh2NativeCodec.h"

#include "../codec/WirehairV2Plan.h"

#include <wirehair/wirehair.h>

#include <algorithm>
#include <cstdint>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <new>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifndef WIREHAIR_WH2_SOURCE_GIT_COMMIT
#define WIREHAIR_WH2_SOURCE_GIT_COMMIT "0000000000000000000000000000000000000000"
#endif

namespace {

using wirehair::wh2_benchmark::FrozenPacketTrace;
using wirehair::wh2_benchmark::FrozenSchedule;
using wirehair::wh2_benchmark::FrozenTraceStatus;
using wirehair_wh2_bench::NativeArm;
using wirehair_wh2_bench::NativeArmSpec;
using wirehair_wh2_bench::NativeRecoveryFixture;
using wirehair_wh2_bench::RecoveryCellResult;
using wirehair_wh2_bench::RecoveryOutcome;
using wirehair_wh2_bench::ResolvedNativeWh2Configuration;

static const char kWorkerSchema[] =
    "wirehair.wh2.mix2-seed-repair-worker.v3";
static const char kDescriptionSchema[] =
    "wirehair.wh2.mix2-seed-repair-worker-description.v3";
static const char kDerivationSchema[] =
    "wirehair.wh2.mix2-seed-repair-derivation-record.v3";
static const char kValidationSchema[] =
    "wirehair.wh2.mix2-seed-repair-validation-record.v3";
static const char kContractSchema[] =
    "wirehair.wh2.mix2-seed-repair-contract.v7";
static const char kValidationRosterSchema[] =
    "wirehair.wh2.mix2-seed-repair-validation-roster.v1";
static const char kValidationRosterSha256[] =
    "030bb1c51e21777266edd4c2349d4a81ccf6e79e2fe4ed9eb75856e16f3387c7";
static const char kProfileSchema[] =
    "wirehair.wh2.mix2-production-profile.v1";
static const char kCandidateArm[] =
    "wirehair2_two07_mix2_graph_b2_v1";
static const char kSeedScheduleCanonical[] =
    "wirehair:wh2:two07-mix2:graph-b2:seed-schedule-v1;"
    "profile=SelectSeedProfile(K,2);"
    "precode=MatrixSeedFromProfile(profile,0,0x763263707265636f);"
    "packet=PacketPeelSeedFromProfile(profile,0x76327265636f7665);"
    "precode_stride=0x9e3779b97f4a7c15;"
    "packet_stride=0x9e3779b9";
static const char kSeedScheduleSha256[] =
    "fe8101781024b4a30797d66cb1512640e6a939586ad5e32f73c5b7ff6f411294";
static const char kCandidateProfileSha256[] =
    "90233b44a0893f96c1a18c19aa61ada052c935a48c6bf7d6a2813065856651f0";
static const uint32_t kBlockBytes = 2u;
static const uint32_t kMaximumK = 64000u;
static const uint32_t kAttemptCount = 256u;
static const uint32_t kScheduleCount = 3u;
static const uint32_t kSelectionRootCount = 6u;
static const uint32_t kValidationRootCount = 3u;
static const uint32_t kSelectionCellCount =
    kSelectionRootCount * kScheduleCount;
static const uint32_t kValidationCellCount =
    kValidationRootCount * kScheduleCount;
static const uint64_t kSourceSeedBase = UINT64_C(0x6d69783273656564);

// Selection retains the v5 root order exactly: its three former training
// roots followed by its three former validation roots.  The cell ordinal is
// global across this concatenated six-root roster.
static const uint64_t kSelectionRoots[kSelectionRootCount] = {
    UINT64_C(0xd1b54a32d192ed03),
    UINT64_C(0x94d049bb133111eb),
    UINT64_C(0x8538ecb5bd456ea3),
    UINT64_C(0xc0ac29b7c97c50dd),
    UINT64_C(0x3f84d5b5b5470917),
    UINT64_C(0x9216d5d98979fb1b)
};
static const uint64_t kValidationRoots[kValidationRootCount] = {
    // Fresh v7 all-K holdout.  Ordinary worker tests must not execute V.
    UINT64_C(0xb501025fdce63900),
    UINT64_C(0x7fb960494dece7de),
    UINT64_C(0x6ad0017d0069e483)
};
static const FrozenSchedule kSchedules[kScheduleCount] = {
    FrozenSchedule::Burst,
    FrozenSchedule::Adversarial,
    FrozenSchedule::RepairOnly
};

struct CellReceipt
{
    uint32_t CellOrdinal = 0u;
    uint32_t RootIndex = 0u;
    uint64_t LossRoot = 0u;
    FrozenSchedule Schedule = FrozenSchedule::Burst;
    uint64_t AttemptedCandidates = 0u;
    std::string TraceSha256;
    std::string Outcome;
    bool HasDecodedExtra = false;
    uint32_t DecodedExtra = 0u;
    int ResultCode = 0;
};

struct AttemptProbe
{
    bool AllOh0 = false;
    uint64_t BasePrecodeSeed = 0u;
    uint32_t BasePacketSeed = 0u;
    std::vector<CellReceipt> Receipts;
    std::string SourceSha256;
    ResolvedNativeWh2Configuration Resolved;
};

std::string Hex64(uint64_t value)
{
    static const char digits[] = "0123456789abcdef";
    std::string text("0x0000000000000000");
    for (unsigned i = 0u; i < 16u; ++i) {
        text[17u - i] = digits[value & 15u];
        value >>= 4;
    }
    return text;
}

std::string Hex32(uint32_t value)
{
    static const char digits[] = "0123456789abcdef";
    std::string text("0x00000000");
    for (unsigned i = 0u; i < 8u; ++i) {
        text[9u - i] = digits[value & 15u];
        value >>= 4;
    }
    return text;
}

std::string ValidationRosterJson()
{
    // Keep keys in canonical JSON order.  Values are derived from the exact
    // arrays used by ValidationRecord(), so a root or schedule permutation
    // changes the preflight identity before any V command is accepted.
    std::string json;
    json.reserve(384u);
    json += "{\"cell_count\":";
    json += std::to_string(kValidationCellCount);
    json += ",\"cell_order\":\"root-major-then-schedule\",\"root_count\":";
    json += std::to_string(kValidationRootCount);
    json += ",\"roots\":[";
    for (uint32_t i = 0u; i < kValidationRootCount; ++i)
    {
        if (i != 0u) json += ',';
        json += '\"';
        json += Hex64(kValidationRoots[i]);
        json += '\"';
    }
    json += "],\"schedule_count\":";
    json += std::to_string(kScheduleCount);
    json += ",\"schedules\":[";
    for (uint32_t i = 0u; i < kScheduleCount; ++i)
    {
        if (i != 0u) json += ',';
        json += '\"';
        json += wirehair::wh2_benchmark::FrozenScheduleName(kSchedules[i]);
        json += '\"';
    }
    json += "],\"schema\":\"";
    json += kValidationRosterSchema;
    json += "\"}";
    return json;
}

std::string ValidationRosterSha256()
{
    return wirehair::wh2_benchmark::Sha256Hex(ValidationRosterJson());
}

bool IsLowerHex(const std::string& value, std::size_t digits)
{
    if (value.size() != digits) return false;
    return std::find_if(
        value.begin(), value.end(), [](char c) {
            return !((c >= '0' && c <= '9') || (c >= 'a' && c <= 'f'));
        }) == value.end();
}

std::string CandidateProfileJson()
{
    return
        "{\"arm\":\"wirehair2_two07_mix2_graph_b2_v1\","
        "\"binary_dense_rows\":12,\"codec\":\"wirehair2_candidate\","
        "\"construction_seed_basis\":\"production-normalized-b2-v1\","
        "\"dense_anchor_layout\":\"two07\","
        "\"dense_identity_corner\":false,\"field\":\"GF(256)\","
        "\"gf256_heavy_rows\":12,\"graph_seed_block_bytes\":2,"
        "\"heavy_family\":\"periodic-cauchy\",\"mix_count\":2,"
        "\"packet_seed_salt\":\"0x76327265636f7665\","
        "\"precode_seed_salt\":\"0x763263707265636f\","
        "\"schema\":\"wirehair.wh2.mix2-production-profile.v1\","
        "\"seed_schedule_sha256\":"
        "\"fe8101781024b4a30797d66cb1512640e6a939586ad5e32f73c5b7ff6f411294\","
        "\"source_hits\":\"certified-by-K\"}";
}

std::string CandidateProfileSha256()
{
    return wirehair::wh2_benchmark::Sha256Hex(CandidateProfileJson());
}

bool SameCandidateConfiguration(
    const wirehair_v2::PrecodeParams& left_params,
    const wirehair_v2::PacketRowConfig& left_packet,
    const wirehair_v2::PrecodeParams& right_params,
    const wirehair_v2::PacketRowConfig& right_packet)
{
    return left_params.BlockCount == right_params.BlockCount &&
        left_params.Staircase == right_params.Staircase &&
        left_params.DenseRows == right_params.DenseRows &&
        left_params.HeavyRows == right_params.HeavyRows &&
        left_params.SourceHits == right_params.SourceHits &&
        left_params.DenseIdentityCorner ==
            right_params.DenseIdentityCorner &&
        left_params.HeavyFamily == right_params.HeavyFamily &&
        left_params.DenseAnchors == right_params.DenseAnchors &&
        left_params.Seed == right_params.Seed &&
        left_packet.PeelSeed == right_packet.PeelSeed &&
        left_packet.MixCount == right_packet.MixCount;
}

bool MakeNormalizedCandidateBase(
    uint32_t K,
    wirehair_v2::PrecodeParams& params,
    wirehair_v2::PacketRowConfig& packet)
{
    if (K < 2u || K > kMaximumK ||
        kBlockBytes != wirehair_v2::kTwo07Mix2GraphSeedBlockBytes ||
        wirehair_v2::kMessagePrecodeSeedSalt !=
            UINT64_C(0x763263707265636f) ||
        wirehair_v2::kMessageRecoveryRowSeedSalt !=
            UINT64_C(0x76327265636f7665) ||
        wirehair::wh2_benchmark::Sha256Hex(kSeedScheduleCanonical) !=
            kSeedScheduleSha256)
    {
        return false;
    }

    // Derive the unstepped graph explicitly from SelectSeedProfile(K, 2).
    // Do not reconstruct it from a payload-width profile or from an already
    // stepped resolved configuration: the frozen map is indexed only by K.
    const wirehair_v2::SeedProfile graph_profile =
        wirehair_v2::SelectSeedProfile(K, kBlockBytes);
    if (graph_profile.BlockCount != K ||
        graph_profile.BlockBytes != kBlockBytes ||
        graph_profile.DenseCount == 0u)
    {
        return false;
    }

    params = wirehair_v2::MakeCertifiedParams(
        K,
        wirehair_v2::MatrixSeedFromProfile(
            graph_profile,
            0u,
            wirehair_v2::kMessagePrecodeSeedSalt));
    params.Staircase = graph_profile.DenseCount;
    params.DenseAnchors = wirehair_v2::DenseAnchorLayout::Two07;
    params.DenseIdentityCorner = false;
    packet.PeelSeed = wirehair_v2::PacketPeelSeedFromProfile(
        graph_profile, wirehair_v2::kMessageRecoveryRowSeedSalt);
    packet.MixCount = 2u;

    return params.BlockCount == K &&
        params.Staircase == graph_profile.DenseCount &&
        params.DenseRows == 12u &&
        params.HeavyRows == 12u &&
        params.SourceHits != 0u &&
        !params.DenseIdentityCorner &&
        params.HeavyFamily ==
            wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy &&
        params.DenseAnchors == wirehair_v2::DenseAnchorLayout::Two07 &&
        packet.MixCount == 2u;
}

bool ReadSelfBinary(std::string& bytes)
{
#if defined(__linux__)
    std::ifstream input("/proc/self/exe", std::ios::binary);
    if (!input) return false;
    bytes.assign(
        std::istreambuf_iterator<char>(input),
        std::istreambuf_iterator<char>());
    return !input.bad() && !bytes.empty();
#else
    (void)bytes;
    return false;
#endif
}

bool SelfBinarySha256(std::string& digest)
{
    try
    {
        std::string bytes;
        if (!ReadSelfBinary(bytes)) return false;
        digest = wirehair::wh2_benchmark::Sha256Hex(bytes);
        return IsLowerHex(digest, 64u);
    }
    catch (const std::bad_alloc&) {
        return false;
    }
    catch (const std::length_error&) {
        return false;
    }
}

bool CandidateTransform(
    uint32_t block_count,
    uint32_t block_bytes,
    wirehair_v2::PrecodeParams& params,
    wirehair_v2::PacketRowConfig& packet,
    void*)
{
    if (block_count < 2u || block_count > kMaximumK ||
        block_bytes != kBlockBytes || params.BlockCount != block_count)
    {
        return false;
    }
    wirehair_v2::PrecodeParams normalized_params;
    wirehair_v2::PacketRowConfig normalized_packet;
    return MakeNormalizedCandidateBase(
            block_count, normalized_params, normalized_packet) &&
        SameCandidateConfiguration(
            params, packet, normalized_params, normalized_packet);
}

NativeArmSpec CandidateSpec(uint32_t attempt)
{
    wirehair_v2::MessagePrecodeEncoderOptions options;
    options.DenseAnchors = wirehair_v2::DenseAnchorLayout::Two07;
    options.RecoveryMixCount = 2u;
    return wirehair_wh2_bench::MakeExperimentalWh2Arm(
        attempt,
        &CandidateTransform,
        nullptr,
        &options,
        wirehair_wh2_bench::NativeWh2BaseKind::ProductionProfile);
}

bool ValidateResolved(
    uint32_t K,
    uint32_t attempt,
    const ResolvedNativeWh2Configuration& resolved)
{
    wirehair_v2::PrecodeParams base_params;
    wirehair_v2::PacketRowConfig base_packet;
    if (!MakeNormalizedCandidateBase(K, base_params, base_packet)) {
        return false;
    }
    const wirehair_v2::PrecodeParams expected_params =
        wirehair_v2::PrecodeParamsForAttempt(base_params, attempt);
    const wirehair_v2::PacketRowConfig expected_packet =
        wirehair_v2::PacketConfigForAttempt(base_packet, attempt);
    return resolved.PrecodeAttempt == attempt &&
        resolved.PacketAttempt == attempt &&
        SameCandidateConfiguration(
            resolved.Params,
            resolved.PacketConfig,
            expected_params,
            expected_packet);
}

uint64_t SourceSeed(uint32_t K)
{
    (void)K;
    return kSourceSeedBase;
}

bool MakeTrace(
    uint32_t K,
    uint32_t cell_ordinal,
    const uint64_t* roots,
    uint32_t root_count,
    FrozenPacketTrace& trace)
{
    if (roots == nullptr || root_count == 0u ||
        cell_ordinal >= root_count * kScheduleCount)
    {
        return false;
    }
    const uint32_t root_index = cell_ordinal / kScheduleCount;
    const uint32_t schedule_index = cell_ordinal % kScheduleCount;
    return wirehair::wh2_benchmark::GenerateFrozenPacketTrace(
        K,
        kBlockBytes,
        500000u,
        roots[root_index],
        kSchedules[schedule_index],
        4u,
        trace) == FrozenTraceStatus::Complete &&
        trace.delivered_ids.size() == static_cast<std::size_t>(K) + 4u &&
        IsLowerHex(trace.trace_sha256, 64u);
}

CellReceipt MakeReceipt(
    uint32_t cell_ordinal,
    const uint64_t* roots,
    const FrozenPacketTrace& trace,
    const RecoveryCellResult& result)
{
    CellReceipt receipt;
    receipt.CellOrdinal = cell_ordinal;
    receipt.RootIndex = cell_ordinal / kScheduleCount;
    receipt.LossRoot = roots[receipt.RootIndex];
    receipt.Schedule = kSchedules[cell_ordinal % kScheduleCount];
    receipt.AttemptedCandidates = trace.attempted_candidates;
    receipt.TraceSha256 = trace.trace_sha256;
    receipt.ResultCode = static_cast<int>(result.Result);
    if (result.Outcome == RecoveryOutcome::Success)
    {
        receipt.HasDecodedExtra = true;
        receipt.DecodedExtra = result.FirstOverhead;
        receipt.Outcome = result.FirstOverhead == 0u ?
            "success" : "need_more_at_oh0";
    }
    else if (result.Outcome == RecoveryOutcome::NeedMoreAtCap) {
        receipt.Outcome = "need_more_at_cap";
    }
    else if (result.Outcome == RecoveryOutcome::ConstructFailed) {
        receipt.Outcome = "construct_failed";
    }
    else if (result.Outcome == RecoveryOutcome::Unsupported) {
        receipt.Outcome = "unsupported";
    }
    else {
        receipt.Outcome = "fatal";
    }
    return receipt;
}

bool IsOh0Success(const CellReceipt& receipt)
{
    return receipt.Outcome == "success" &&
        receipt.HasDecodedExtra && receipt.DecodedExtra == 0u &&
        receipt.ResultCode == static_cast<int>(Wirehair_Success);
}

bool ProbeAttempt(
    uint32_t K,
    uint32_t attempt,
    const uint64_t* roots,
    uint32_t root_count,
    bool stop_at_first_failure,
    AttemptProbe& output,
    std::string& error)
{
    if (K < 2u || K > kMaximumK || attempt >= kAttemptCount ||
        roots == nullptr || root_count == 0u ||
        root_count > std::numeric_limits<uint32_t>::max() / kScheduleCount)
    {
        error = "probe coordinate is outside the frozen domain";
        return false;
    }
    const uint32_t cell_count = root_count * kScheduleCount;
    std::vector<uint8_t> source;
    if (!wirehair_wh2_bench::MakeDeterministicSource(
            K, kBlockBytes, SourceSeed(K), source))
    {
        error = "cannot construct deterministic source";
        return false;
    }
    output.SourceSha256 = wirehair::wh2_benchmark::Sha256Hex(
        source.data(), source.size());
    if (!IsLowerHex(output.SourceSha256, 64u)) {
        error = "cannot hash deterministic source";
        return false;
    }
    const NativeArmSpec spec = CandidateSpec(attempt);
    wirehair_v2::PrecodeParams base_params;
    wirehair_v2::PacketRowConfig base_packet;
    if (!wirehair_wh2_bench::ResolveNativeWh2Configuration(
            spec, K, kBlockBytes, output.Resolved) ||
        !ValidateResolved(K, attempt, output.Resolved) ||
        !MakeNormalizedCandidateBase(K, base_params, base_packet))
    {
        error = "candidate configuration differs from the frozen profile";
        return false;
    }
    output.BasePrecodeSeed = base_params.Seed;
    output.BasePacketSeed = base_packet.PeelSeed;

    NativeArm arm;
    const WirehairResult construction = arm.Initialize(
        spec, K, kBlockBytes, source);
    if (construction != Wirehair_Success &&
        construction != Wirehair_BadPeelSeed)
    {
        error = "candidate construction returned unsupported/fatal result";
        return false;
    }

    output.Receipts.clear();
    for (uint32_t cell = 0u; cell < cell_count; ++cell)
    {
        FrozenPacketTrace trace;
        if (!MakeTrace(K, cell, roots, root_count, trace)) {
            error = "cannot generate frozen packet trace";
            return false;
        }
        RecoveryCellResult result;
        if (construction == Wirehair_BadPeelSeed)
        {
            result.Outcome = RecoveryOutcome::ConstructFailed;
            result.Result = construction;
        }
        else
        {
            NativeRecoveryFixture fixture;
            const WirehairResult initialized = fixture.Initialize(
                arm, trace.delivered_ids);
            if (initialized != Wirehair_Success) {
                error = "cannot initialize native recovery fixture";
                return false;
            }
            result = fixture.RunNested();
            if (result.Outcome == RecoveryOutcome::Fatal ||
                result.Outcome == RecoveryOutcome::Unsupported)
            {
                error = "native recovery returned unsupported/fatal result";
                return false;
            }
        }
        const CellReceipt receipt = MakeReceipt(cell, roots, trace, result);
        if (IsOh0Success(receipt)) {
            output.Receipts.push_back(receipt);
            continue;
        }
        if (stop_at_first_failure) {
            output.Receipts.clear();
            output.Receipts.push_back(receipt);
            output.AllOh0 = false;
            return true;
        }
        output.Receipts.push_back(receipt);
    }
    output.AllOh0 = std::all_of(
        output.Receipts.begin(), output.Receipts.end(), IsOh0Success);
    return output.Receipts.size() == cell_count;
}

std::string CellReceiptJson(const CellReceipt& receipt, uint32_t attempt)
{
    std::string json;
    json.reserve(384u);
    json += "{\"attempted_candidates\":";
    json += std::to_string(receipt.AttemptedCandidates);
    json += ",\"cell_ordinal\":";
    json += std::to_string(receipt.CellOrdinal);
    json += ",\"construction_attempt\":";
    json += std::to_string(attempt);
    json += ",\"decoded_extra\":";
    json += receipt.HasDecodedExtra ?
        std::to_string(receipt.DecodedExtra) : "null";
    json += ",\"loss_ppm\":500000,\"loss_root\":\"";
    json += Hex64(receipt.LossRoot);
    json += "\",\"outcome\":\"";
    json += receipt.Outcome;
    json += "\",\"result_code\":";
    json += std::to_string(receipt.ResultCode);
    json += ",\"root_index\":";
    json += std::to_string(receipt.RootIndex);
    json += ",\"schedule\":\"";
    json += wirehair::wh2_benchmark::FrozenScheduleName(receipt.Schedule);
    json += "\",\"trace_sha256\":\"";
    json += receipt.TraceSha256;
    json += "\"}";
    return json;
}

void AppendReceiptArray(
    std::string& json,
    const std::vector<CellReceipt>& receipts,
    uint32_t attempt)
{
    json += '[';
    for (std::size_t i = 0u; i < receipts.size(); ++i)
    {
        if (i != 0u) json += ',';
        json += CellReceiptJson(receipts[i], attempt);
    }
    json += ']';
}

std::string DeriveRecord(
    std::size_t ordinal,
    uint32_t K,
    const std::string& binary_sha256,
    std::string& error)
{
    std::vector<CellReceipt> lower_witnesses;
    AttemptProbe selected;
    uint32_t selected_attempt = kAttemptCount;
    for (uint32_t attempt = 0u; attempt < kAttemptCount; ++attempt)
    {
        AttemptProbe probe;
        if (!ProbeAttempt(
                K,
                attempt,
                kSelectionRoots,
                kSelectionRootCount,
                true,
                probe,
                error))
        {
            return std::string();
        }
        if (probe.AllOh0)
        {
            if (probe.Receipts.size() != kSelectionCellCount) {
                error = "selected attempt omitted selection successes";
                return std::string();
            }
            selected_attempt = attempt;
            selected = std::move(probe);
            break;
        }
        if (probe.Receipts.size() != 1u ||
            IsOh0Success(probe.Receipts[0]))
        {
            error = "lower attempt did not retain one failure witness";
            return std::string();
        }
        lower_witnesses.push_back(probe.Receipts[0]);
    }
    if (selected_attempt >= kAttemptCount ||
        lower_witnesses.size() != selected_attempt)
    {
        error = "no uint8 construction attempt passes the selection census";
        return std::string();
    }

    std::string json;
    json.reserve(4096u + lower_witnesses.size() * 384u);
    json += "{\"K\":";
    json += std::to_string(K);
    json += ",\"base_packet_seed\":\"";
    json += Hex32(selected.BasePacketSeed);
    json += "\",\"base_precode_seed\":\"";
    json += Hex64(selected.BasePrecodeSeed);
    json += "\",\"candidate_profile_sha256\":\"";
    json += CandidateProfileSha256();
    json += "\",\"effective_packet_seed\":\"";
    json += Hex32(selected.Resolved.PacketConfig.PeelSeed);
    json += "\",\"effective_precode_seed\":\"";
    json += Hex64(selected.Resolved.Params.Seed);
    json += "\",\"lower_attempt_failure_witnesses\":[";
    for (std::size_t i = 0u; i < lower_witnesses.size(); ++i)
    {
        if (i != 0u) json += ',';
        json += CellReceiptJson(lower_witnesses[i], static_cast<uint32_t>(i));
    }
    json += "],\"mode\":\"derive\",\"ordinal\":";
    json += std::to_string(ordinal);
    json += ",\"schema\":\"";
    json += kDerivationSchema;
    json += "\",\"selected_attempt\":";
    json += std::to_string(selected_attempt);
    json += ",\"selected_successes\":";
    AppendReceiptArray(json, selected.Receipts, selected_attempt);
    json += ",\"source_sha256\":\"";
    json += selected.SourceSha256;
    json += "\",\"worker_binary_sha256\":\"";
    json += binary_sha256;
    json += "\"}";
    return json;
}

std::string ValidationRecord(
    std::size_t ordinal,
    uint32_t K,
    uint32_t attempt,
    const std::string& binary_sha256,
    std::string& error)
{
    AttemptProbe probe;
    if (!ProbeAttempt(
            K,
            attempt,
            kValidationRoots,
            kValidationRootCount,
            false,
            probe,
            error) ||
        probe.Receipts.size() != kValidationCellCount)
    {
        return std::string();
    }
    std::string json;
    json.reserve(4096u);
    json += "{\"K\":";
    json += std::to_string(K);
    json += ",\"all_oh0_success\":";
    json += probe.AllOh0 ? "true" : "false";
    json += ",\"base_packet_seed\":\"";
    json += Hex32(probe.BasePacketSeed);
    json += "\",\"base_precode_seed\":\"";
    json += Hex64(probe.BasePrecodeSeed);
    json += "\",\"candidate_profile_sha256\":\"";
    json += CandidateProfileSha256();
    json += "\",\"cells\":";
    AppendReceiptArray(json, probe.Receipts, attempt);
    json += ",\"construction_attempt\":";
    json += std::to_string(attempt);
    json += ",\"effective_packet_seed\":\"";
    json += Hex32(probe.Resolved.PacketConfig.PeelSeed);
    json += "\",\"effective_precode_seed\":\"";
    json += Hex64(probe.Resolved.Params.Seed);
    json += "\",\"mode\":\"validate\",\"ordinal\":";
    json += std::to_string(ordinal);
    json += ",\"schema\":\"";
    json += kValidationSchema;
    json += "\",\"source_sha256\":\"";
    json += probe.SourceSha256;
    json += "\",\"worker_binary_sha256\":\"";
    json += binary_sha256;
    json += "\"}";
    return json;
}

std::string Description(const std::string& binary_sha256)
{
    std::string json;
    json.reserve(1536u);
    json += "{\"binary_sha256\":\"";
    json += binary_sha256;
    json += "\",\"candidate_profile\":";
    json += CandidateProfileJson();
    json += ",\"candidate_profile_sha256\":\"";
    json += CandidateProfileSha256();
    json += "\",\"contract_schema\":\"";
    json += kContractSchema;
    json += "\",\"derivation_schema\":\"";
    json += kDerivationSchema;
    json += "\",\"protocol\":\"D ordinal K | V ordinal K attempt | Q\","
        "\"schema\":\"";
    json += kDescriptionSchema;
    json += "\",\"source_git_commit\":\"";
    json += WIREHAIR_WH2_SOURCE_GIT_COMMIT;
    json += "\",\"validation_roster_schema\":\"";
    json += kValidationRosterSchema;
    json += "\",\"validation_roster_sha256\":\"";
    json += ValidationRosterSha256();
    json += "\",\"validation_schema\":\"";
    json += kValidationSchema;
    json += "\",\"worker_schema\":\"";
    json += kWorkerSchema;
    json += "\"}";
    return json;
}

bool ParseUnsigned(const std::string& text, uint64_t maximum, uint64_t& value)
{
    if (text.empty() || (text.size() > 1u && text[0] == '0')) return false;
    uint64_t parsed = 0u;
    for (char c : text)
    {
        if (c < '0' || c > '9') return false;
        const uint64_t digit = static_cast<uint64_t>(c - '0');
        if (parsed > (maximum - digit) / 10u) return false;
        parsed = parsed * 10u + digit;
    }
    value = parsed;
    return true;
}

static const std::size_t kMaximumCommandBytes = 64u;

enum class CommandReadStatus
{
    Line,
    End,
    TooLong,
    Error
};

CommandReadStatus ReadCommandLine(std::string& line)
{
    char buffer[kMaximumCommandBytes];
    std::size_t length = 0u;
    for (;;)
    {
        const std::char_traits<char>::int_type next = std::cin.get();
        if (next == std::char_traits<char>::eof())
        {
            if (std::cin.bad() || !std::cin.eof()) {
                return CommandReadStatus::Error;
            }
            if (length == 0u) return CommandReadStatus::End;
            line.assign(buffer, length);
            return CommandReadStatus::Line;
        }
        if (next == '\n')
        {
            line.assign(buffer, length);
            return CommandReadStatus::Line;
        }
        if (length >= kMaximumCommandBytes) {
            return CommandReadStatus::TooLong;
        }
        buffer[length++] = std::char_traits<char>::to_char_type(next);
    }
}

bool RunProtocol(const std::string& binary_sha256)
{
    std::string line;
    for (;;)
    {
        const CommandReadStatus read_status = ReadCommandLine(line);
        if (read_status == CommandReadStatus::End)
        {
            std::cerr << "mix2 seed-repair worker: EOF before Q\n";
            return false;
        }
        if (read_status == CommandReadStatus::TooLong)
        {
            std::cerr << "mix2 seed-repair worker: command exceeds "
                << kMaximumCommandBytes << "-byte limit\n";
            return false;
        }
        if (read_status == CommandReadStatus::Error)
        {
            std::cerr << "mix2 seed-repair worker: command read failed\n";
            return false;
        }
        if (line == "Q") return true;
        std::istringstream input(line);
        std::string kind;
        std::string ordinal_text;
        std::string K_text;
        std::string attempt_text;
        std::string extra;
        if (!(input >> kind >> ordinal_text >> K_text) ||
            (kind != "D" && kind != "V") ||
            (kind == "V" && !(input >> attempt_text)) ||
            (input >> extra))
        {
            std::cerr << "mix2 seed-repair worker: malformed command\n";
            return false;
        }
        uint64_t ordinal_wide = 0u;
        uint64_t K_wide = 0u;
        uint64_t attempt_wide = 0u;
        if (!ParseUnsigned(
                ordinal_text,
                static_cast<uint64_t>(std::numeric_limits<std::size_t>::max()),
                ordinal_wide) ||
            !ParseUnsigned(K_text, kMaximumK, K_wide) ||
            K_wide < 2u ||
            ordinal_wide != K_wide - 2u ||
            (kind == "V" &&
             !ParseUnsigned(attempt_text, kAttemptCount - 1u, attempt_wide)))
        {
            std::cerr << "mix2 seed-repair worker: command coordinate invalid\n";
            return false;
        }
        std::string error;
        const std::string result = kind == "D" ?
            DeriveRecord(
                static_cast<std::size_t>(ordinal_wide),
                static_cast<uint32_t>(K_wide), binary_sha256, error) :
            ValidationRecord(
                static_cast<std::size_t>(ordinal_wide),
                static_cast<uint32_t>(K_wide),
                static_cast<uint32_t>(attempt_wide), binary_sha256, error);
        if (result.empty())
        {
            std::cerr << "mix2 seed-repair worker: " << error << '\n';
            return false;
        }
        std::cout << result << '\n' << std::flush;
        if (!std::cout) return false;
    }
}

bool SelfTest(const std::string& binary_sha256)
{
    if (CandidateProfileSha256() != kCandidateProfileSha256 ||
        wirehair::wh2_benchmark::Sha256Hex(kSeedScheduleCanonical) !=
            kSeedScheduleSha256 ||
        ValidationRosterSha256() != kValidationRosterSha256 ||
        CandidateProfileJson().find(kProfileSchema) == std::string::npos ||
        CandidateProfileJson().find(kCandidateArm) == std::string::npos)
    {
        return false;
    }
    wirehair_v2::PrecodeParams base_params;
    wirehair_v2::PacketRowConfig base_packet;
    if (!MakeNormalizedCandidateBase(8u, base_params, base_packet) ||
        base_params.Seed != UINT64_C(0x1dcd06e5a961880b) ||
        base_packet.PeelSeed != UINT32_C(0x77d69c8c))
    {
        return false;
    }
    const uint32_t payload_widths[] = { 2u, 17u, 1280u };
    for (uint32_t payload_width : payload_widths)
    {
        const wirehair_v2::SeedProfile payload_profile =
            wirehair_v2::SelectSeedProfile(8u, payload_width);
        wirehair_v2::MessagePrecodeEncoderOptions options;
        options.DenseAnchors = wirehair_v2::DenseAnchorLayout::Two07;
        options.RecoveryMixCount = 2u;
        wirehair_v2::PrecodeParams production_params;
        wirehair_v2::PacketRowConfig production_packet;
        wirehair_v2::MakeMessagePrecodeConfiguration(
            payload_profile, options, production_params, production_packet);
        production_params.Staircase = payload_profile.DenseCount;
        if (!SameCandidateConfiguration(
                production_params,
                production_packet,
                base_params,
                base_packet))
        {
            return false;
        }
    }
    ResolvedNativeWh2Configuration resolved;
    if (!wirehair_wh2_bench::ResolveNativeWh2Configuration(
            CandidateSpec(255u), 8u, kBlockBytes, resolved) ||
        !ValidateResolved(8u, 255u, resolved))
    {
        return false;
    }
    std::string error;
    const std::string record = DeriveRecord(0u, 2u, binary_sha256, error);
    return !record.empty() && error.empty() &&
        record.find("\"K\":2") != std::string::npos &&
        record.find("\"base_packet_seed\":\"0x9799409a\"") !=
            std::string::npos &&
        record.find("\"base_precode_seed\":\"0xecbb2b8a0e18da1e\"") !=
            std::string::npos &&
        record.find("\"effective_packet_seed\":\"0x278c881b\"") !=
            std::string::npos &&
        record.find("\"effective_precode_seed\":\"0x7cae730f87b736db\"") !=
            std::string::npos &&
        record.find("\"ordinal\":0") != std::string::npos &&
        record.find("\"selected_attempt\":9") != std::string::npos &&
        record.find(kDerivationSchema) != std::string::npos;
}

int Usage()
{
    std::cerr << "usage: wirehair_wh2_mix2_seed_repair_worker "
        "--describe | --worker | --self-test\n";
    return 2;
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        std::string binary_sha256;
        if (!SelfBinarySha256(binary_sha256)) {
            std::cerr << "mix2 seed-repair worker: cannot hash self binary\n";
            return 1;
        }
        if (argc == 2 && std::string(argv[1]) == "--describe") {
            std::cout << Description(binary_sha256) << '\n';
            return std::cout ? 0 : 1;
        }
        if (argc == 2 && std::string(argv[1]) == "--worker") {
            return RunProtocol(binary_sha256) ? 0 : 1;
        }
        if (argc == 2 && std::string(argv[1]) == "--self-test") {
            return SelfTest(binary_sha256) ? 0 : 1;
        }
        return Usage();
    }
    catch (const std::bad_alloc&)
    {
        std::cerr << "mix2 seed-repair worker: allocation failed\n";
    }
    catch (const std::length_error&)
    {
        std::cerr << "mix2 seed-repair worker: allocation length invalid\n";
    }
    catch (const std::exception&)
    {
        std::cerr << "mix2 seed-repair worker: standard exception\n";
    }
    catch (...)
    {
        std::cerr << "mix2 seed-repair worker: unexpected exception\n";
    }
    return 1;
}
