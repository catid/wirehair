#include "V2TinyDenseOracle.h"
#include "WirehairV2Plan.h"
#include "WirehairV2Precode.h"
#include "WirehairV2Solve.h"

#include <wirehair/wirehair.h>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace {

class TinyDirectStateReset
{
public:
    ~TinyDirectStateReset()
    {
        wirehair_v2::SetTinyDirectSolveForTesting(false);
        wirehair_v2::SetTinyDirectSolveAllocationFailureForTesting(
            -1, false);
        wirehair_v2::SetOddPacketPeelSeedXorForTesting(0u);
        wirehair_v2::SetPacketRowSeedAvalancheForTesting(false);
        (void)wirehair_v2::SetPacketRowSeedMultiplierForTesting(1u);
    }
};

struct Fixture
{
    uint32_t K = 0u;
    uint32_t BlockBytes = 0u;
    wirehair_v2::PrecodeSystem System;
    wirehair_v2::PacketRowConfig Config;
    std::vector<uint8_t> Message;
    std::vector<wirehair_v2::SolvePacket> Packets;
};

bool MakeFixtureFromParams(
    const wirehair_v2::PrecodeParams& params,
    uint32_t block_bytes,
    Fixture& out,
    uint32_t mix_count = wirehair_v2::kCertifiedPacketMixCount)
{
    wirehair_v2::PacketRowConfig base_config;
    base_config.PeelSeed =
        UINT32_C(0x6d2b79f5) ^ params.BlockCount ^ block_bytes;
    base_config.MixCount = mix_count;

    out = Fixture();
    out.K = params.BlockCount;
    out.BlockBytes = block_bytes;
    if (wirehair_v2::SelectSystematicConfiguration(
            params,
            base_config,
            out.System,
            out.Config) != Wirehair_Success)
    {
        std::fprintf(stderr,
            "tiny-direct: configuration failed K=%u bb=%u\n",
            params.BlockCount, block_bytes);
        return false;
    }

    out.Message.resize((size_t)out.K * block_bytes);
    for (size_t i = 0u; i < out.Message.size(); ++i) {
        out.Message[i] = (uint8_t)(
            i * 73u + (i >> 5) + out.K * 19u + block_bytes * 7u);
    }
    out.Packets.resize(out.K);
    for (uint32_t id = 0u; id < out.K; ++id)
    {
        out.Packets[id].BlockId = id;
        out.Packets[id].Data =
            out.Message.data() + (size_t)id * block_bytes;
    }
    return true;
}

bool MakeFixture(
    uint32_t K,
    uint32_t block_bytes,
    wirehair_v2::DenseAnchorLayout layout,
    Fixture& out)
{
    wirehair_v2::PrecodeParams params =
        wirehair_v2::MakeCertifiedParams(
            K,
            UINT64_C(0x74696e7964697265) ^
                ((uint64_t)K << 32) ^ (uint32_t)layout);
    params.DenseAnchors = layout;
    return MakeFixtureFromParams(params, block_bytes, out);
}

bool SameUntimedStats(
    const wirehair_v2::PrecodeSolveStats& a,
    const wirehair_v2::PrecodeSolveStats& b)
{
    return a.PacketRows == b.PacketRows &&
        a.PeeledColumns == b.PeeledColumns &&
        a.InactivatedColumns == b.InactivatedColumns &&
        a.ResidualRows == b.ResidualRows &&
        a.ResidualRank == b.ResidualRank &&
        a.BinaryResidualRank == b.BinaryResidualRank &&
        a.BinaryRowReferences == b.BinaryRowReferences &&
        a.BinaryRowStorageBytes == b.BinaryRowStorageBytes &&
        a.BinaryAdjacencyStorageBytes == b.BinaryAdjacencyStorageBytes &&
        a.BinaryRowStorageAllocations == b.BinaryRowStorageAllocations &&
        a.BinaryAdjacencyStorageAllocations ==
            b.BinaryAdjacencyStorageAllocations &&
        a.BlockXors == b.BlockXors &&
        a.BlockMulAdds == b.BlockMulAdds &&
        a.PacketSeedAttempt == b.PacketSeedAttempt;
}

WirehairResult RunWithMode(
    const Fixture& fixture,
    bool direct,
    const std::vector<wirehair_v2::SolvePacket>& packets,
    std::vector<uint8_t>& output,
    wirehair_v2::PrecodeSolveStats* stats = nullptr,
    wirehair_v2::PrecodeSolveResumeState* resume = nullptr)
{
    wirehair_v2::SetTinyDirectSolveForTesting(direct);
    const WirehairResult result = wirehair_v2::SolvePrecodeSystem(
        fixture.System,
        fixture.Config,
        packets,
        fixture.BlockBytes,
        output,
        stats,
        resume);
    wirehair_v2::SetTinyDirectSolveForTesting(false);
    return result;
}

bool SolveWithMode(
    const Fixture& fixture,
    bool direct,
    const std::vector<wirehair_v2::SolvePacket>& packets,
    std::vector<uint8_t>& output,
    wirehair_v2::PrecodeSolveStats* stats = nullptr,
    wirehair_v2::PrecodeSolveResumeState* resume = nullptr)
{
    return RunWithMode(
        fixture, direct, packets, output, stats, resume) == Wirehair_Success;
}

bool CheckForcedOracleMatrix()
{
    static const wirehair_v2::DenseAnchorLayout kLayouts[] = {
        wirehair_v2::DenseAnchorLayout::Disabled,
        wirehair_v2::DenseAnchorLayout::Two07
    };

    uint32_t checked = 0u;
    const auto check_case = [&](
        uint32_t K,
        uint32_t block_bytes,
        wirehair_v2::DenseAnchorLayout layout) -> bool {
        Fixture fixture;
        if (!MakeFixture(K, block_bytes, layout, fixture)) {
            return false;
        }
        std::vector<uint8_t> general;
        wirehair_v2::PrecodeSolveStats general_stats;
        wirehair_v2::ResetTinyDirectSolveCountersForTesting();
        if (!SolveWithMode(
                fixture,
                false,
                fixture.Packets,
                general,
                &general_stats) ||
            wirehair_v2::TinyDirectSolveAttemptsForTesting() != 0u)
        {
            std::fprintf(stderr,
                "tiny-direct: forced-off contamination K=%u bb=%u "
                "layout=%u\n",
                K, block_bytes, (unsigned)layout);
            return false;
        }

        std::vector<uint8_t> oracle;
        if (wirehair_v2::test::SolvePrecodeSystemTinyDenseOracle(
                fixture.System,
                fixture.Config,
                fixture.Packets,
                block_bytes,
                oracle) != Wirehair_Success ||
            oracle != general)
        {
            std::fprintf(stderr,
                "tiny-direct: baseline/oracle mismatch K=%u bb=%u "
                "layout=%u\n",
                K, block_bytes, (unsigned)layout);
            return false;
        }

        std::vector<uint8_t> direct;
        wirehair_v2::PrecodeSolveStats direct_stats;
        wirehair_v2::ResetTinyDirectSolveCountersForTesting();
        if (!SolveWithMode(
                fixture,
                true,
                fixture.Packets,
                direct,
                &direct_stats) ||
            direct != general ||
            wirehair_v2::TinyDirectSolveAttemptsForTesting() != 1u ||
            wirehair_v2::TinyDirectSolveCompletionsForTesting() != 1u ||
            wirehair_v2::TinyDirectSolveFallbacksForTesting() != 0u)
        {
            std::fprintf(stderr,
                "tiny-direct: forced-on/oracle mismatch K=%u bb=%u "
                "layout=%u attempts=%llu complete=%llu fallback=%llu\n",
                K, block_bytes, (unsigned)layout,
                (unsigned long long)
                    wirehair_v2::TinyDirectSolveAttemptsForTesting(),
                (unsigned long long)
                    wirehair_v2::TinyDirectSolveCompletionsForTesting(),
                (unsigned long long)
                    wirehair_v2::TinyDirectSolveFallbacksForTesting());
            return false;
        }
        const uint32_t L = K + fixture.System.Params.Staircase +
            fixture.System.Params.DenseRows +
            fixture.System.Params.HeavyRows;
        if (direct_stats.PacketRows != K ||
            direct_stats.PeeledColumns != 0u ||
            direct_stats.InactivatedColumns != L ||
            direct_stats.ResidualRows != L ||
            direct_stats.ResidualRank != L ||
            direct_stats.BinaryRowStorageAllocations != 1u ||
            direct_stats.BinaryAdjacencyStorageAllocations != 0u ||
            !wirehair_v2::VerifyPrecodeSolution(
                fixture.System,
                fixture.Config,
                fixture.Packets,
                direct.data(),
                block_bytes))
        {
            std::fprintf(stderr,
                "tiny-direct: stats/solution mismatch K=%u bb=%u "
                "layout=%u rows=%u rank=%u L=%u\n",
                K, block_bytes, (unsigned)layout,
                direct_stats.ResidualRows,
                direct_stats.ResidualRank,
                L);
            return false;
        }
        ++checked;
        return true;
    };

    for (wirehair_v2::DenseAnchorLayout layout : kLayouts) {
        for (uint32_t K = 2u; K <= 100u; ++K) {
            if (!check_case(K, 2u, layout)) {
                return false;
            }
        }
    }
    static const uint32_t kWideBlockCounts[] = {2u, 8u, 32u, 100u};
    static const uint32_t kWideBlockBytes[] = {64u, 1280u, 4096u};
    for (wirehair_v2::DenseAnchorLayout layout : kLayouts) {
        for (uint32_t K : kWideBlockCounts) {
            for (uint32_t block_bytes : kWideBlockBytes) {
                if (!check_case(K, block_bytes, layout)) {
                    return false;
                }
            }
        }
    }

    // Exercise the cached-runtime/trusted-system seam explicitly.
    Fixture cached;
    if (!MakeFixture(
            64u, 13u, wirehair_v2::DenseAnchorLayout::Two07, cached))
    {
        return false;
    }
    const uint32_t P = cached.System.Params.Staircase +
        cached.System.Params.DenseRows + cached.System.Params.HeavyRows;
    wirehair_v2::PacketRowRuntime runtime;
    if (!runtime.Initialize(cached.K, P, cached.Config.MixCount)) {
        return false;
    }
    std::vector<uint8_t> cached_output;
    wirehair_v2::ResetTinyDirectSolveCountersForTesting();
    wirehair_v2::SetTinyDirectSolveForTesting(true);
    const WirehairResult cached_result =
        wirehair_v2::SolvePrecodeSystemForValidatedSystemWithRuntime(
            cached.System,
            cached.Config,
            runtime,
            cached.Packets,
            cached.BlockBytes,
            cached_output);
    wirehair_v2::SetTinyDirectSolveForTesting(false);
    if (cached_result != Wirehair_Success ||
        wirehair_v2::TinyDirectSolveAttemptsForTesting() != 1u ||
        wirehair_v2::TinyDirectSolveCompletionsForTesting() != 1u ||
        !wirehair_v2::VerifyPrecodeSolution(
            cached.System,
            cached.Config,
            cached.Packets,
            cached_output.data(),
            cached.BlockBytes))
    {
        std::fprintf(stderr,
            "tiny-direct: cached-runtime path failed result=%d\n",
            (int)cached_result);
        return false;
    }

    std::printf(
        "tiny-direct forced/off/on dense-oracle matrix: %u cases PASS\n",
        checked);
    return true;
}

uint64_t Mix64(uint64_t value)
{
    value += UINT64_C(0x9e3779b97f4a7c15);
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

bool MakeEncodedPacketCorpus(
    const Fixture& fixture,
    const std::vector<uint8_t>& intermediate,
    uint32_t packet_count,
    uint64_t seed,
    std::vector<uint8_t>& storage,
    std::vector<wirehair_v2::SolvePacket>& packets)
{
    const uint32_t P = fixture.System.Params.Staircase +
        fixture.System.Params.DenseRows + fixture.System.Params.HeavyRows;
    wirehair_v2::PacketRowRuntime runtime;
    if (!runtime.Initialize(fixture.K, P, fixture.Config.MixCount)) {
        return false;
    }
    storage.assign((size_t)packet_count * fixture.BlockBytes, 0u);
    packets.clear();
    packets.reserve(packet_count);
    for (uint32_t index = 0u; index < packet_count; ++index)
    {
        uint32_t block_id = 0u;
        bool duplicate = false;
        do
        {
            block_id = (uint32_t)Mix64(
                seed + (uint64_t)index * UINT64_C(0x632be59bd9b4e019));
            seed = Mix64(seed ^ block_id);
            duplicate = false;
            for (const wirehair_v2::SolvePacket& packet : packets) {
                if (packet.BlockId == block_id) {
                    duplicate = true;
                    break;
                }
            }
        } while (duplicate);

        uint8_t* const data = storage.data() +
            (size_t)index * fixture.BlockBytes;
        if (!wirehair_v2::EvaluatePacketBlockForValidatedSystemWithRuntime(
                fixture.System,
                fixture.Config,
                runtime,
                intermediate.data(),
                fixture.BlockBytes,
                block_id,
                data))
        {
            return false;
        }
        wirehair_v2::SolvePacket packet;
        packet.BlockId = block_id;
        packet.Data = data;
        packets.push_back(packet);
    }
    return true;
}

bool ComparePacketCorpus(
    const Fixture& fixture,
    const std::vector<wirehair_v2::SolvePacket>& packets,
    const char* label)
{
    std::vector<uint8_t> general(7u, 0xa5u);
    std::vector<uint8_t> oracle = general;
    std::vector<uint8_t> direct = general;
    const WirehairResult general_result = RunWithMode(
        fixture, false, packets, general);
    const WirehairResult oracle_result =
        wirehair_v2::test::SolvePrecodeSystemTinyDenseOracle(
            fixture.System,
            fixture.Config,
            packets,
            fixture.BlockBytes,
            oracle);
    wirehair_v2::ResetTinyDirectSolveCountersForTesting();
    const WirehairResult direct_result = RunWithMode(
        fixture, true, packets, direct);
    if (direct_result != general_result || direct_result != oracle_result ||
        direct != general || direct != oracle ||
        wirehair_v2::TinyDirectSolveAttemptsForTesting() != 1u ||
        wirehair_v2::TinyDirectSolveCompletionsForTesting() != 1u ||
        wirehair_v2::TinyDirectSolveFallbacksForTesting() != 0u)
    {
        std::fprintf(stderr,
            "tiny-direct: packet corpus mismatch %s K=%u bb=%u "
            "general=%d oracle=%d direct=%d\n",
            label,
            fixture.K,
            fixture.BlockBytes,
            (int)general_result,
            (int)oracle_result,
            (int)direct_result);
        return false;
    }
    return direct_result != Wirehair_Success ||
        wirehair_v2::VerifyPrecodeSolution(
            fixture.System,
            fixture.Config,
            packets,
            direct.data(),
            fixture.BlockBytes);
}

bool CheckRandomDecoderShapes()
{
    static const uint32_t kBlockCounts[] = {2u, 8u, 32u, 64u, 100u};
    static const uint32_t kBlockBytes[] = {2u, 64u};
    static const wirehair_v2::DenseAnchorLayout kLayouts[] = {
        wirehair_v2::DenseAnchorLayout::Disabled,
        wirehair_v2::DenseAnchorLayout::Two07
    };
    static const uint32_t kOverheads[] = {0u, 2u, 4u};
    uint32_t comparisons = 0u;
    for (wirehair_v2::DenseAnchorLayout layout : kLayouts)
    {
        for (uint32_t K : kBlockCounts)
        {
            for (uint32_t block_bytes : kBlockBytes)
            {
                Fixture fixture;
                if (!MakeFixture(K, block_bytes, layout, fixture)) {
                    return false;
                }
                std::vector<uint8_t> intermediate;
                if (!SolveWithMode(
                        fixture,
                        false,
                        fixture.Packets,
                        intermediate))
                {
                    return false;
                }
                std::vector<uint8_t> storage;
                std::vector<wirehair_v2::SolvePacket> packets;
                if (!MakeEncodedPacketCorpus(
                        fixture,
                        intermediate,
                        K + 4u,
                        UINT64_C(0x72616e646f6d6964) ^
                            ((uint64_t)K << 32) ^ block_bytes ^
                            (uint32_t)layout,
                        storage,
                        packets))
                {
                    return false;
                }
                for (uint32_t overhead : kOverheads)
                {
                    const std::vector<wirehair_v2::SolvePacket> prefix(
                        packets.begin(), packets.begin() + K + overhead);
                    if (!ComparePacketCorpus(fixture, prefix, "encoded")) {
                        return false;
                    }

                    std::vector<uint8_t> arbitrary_storage(
                        (size_t)(K + overhead) * block_bytes);
                    for (size_t i = 0u; i < arbitrary_storage.size(); ++i) {
                        arbitrary_storage[i] = (uint8_t)Mix64(
                            i ^ ((uint64_t)K << 40) ^
                            ((uint64_t)block_bytes << 8) ^ overhead);
                    }
                    std::vector<wirehair_v2::SolvePacket> arbitrary = prefix;
                    for (uint32_t row = 0u; row < K + overhead; ++row) {
                        arbitrary[row].Data = arbitrary_storage.data() +
                            (size_t)row * block_bytes;
                    }
                    if (!ComparePacketCorpus(fixture, arbitrary, "arbitrary")) {
                        return false;
                    }
                    comparisons += 2u;
                }
            }
        }
    }
    std::printf(
        "tiny-direct random loss/repair/RHS differential: %u cases PASS\n",
        comparisons);
    return true;
}

bool CheckPacketRowHookAgreement()
{
    Fixture fixture;
    if (!MakeFixture(
            32u, 13u, wirehair_v2::DenseAnchorLayout::Two07, fixture))
    {
        return false;
    }
    std::vector<uint8_t> intermediate;
    if (!SolveWithMode(
            fixture, false, fixture.Packets, intermediate) ||
        !wirehair_v2::SetPacketRowSeedMultiplierForTesting(
            UINT32_C(0x9e3779b1)))
    {
        return false;
    }
    wirehair_v2::SetPacketRowSeedAvalancheForTesting(true);
    wirehair_v2::SetOddPacketPeelSeedXorForTesting(UINT32_C(0xa53c91e7));

    std::vector<uint8_t> storage;
    std::vector<wirehair_v2::SolvePacket> packets;
    if (!MakeEncodedPacketCorpus(
            fixture,
            intermediate,
            fixture.K + 4u,
            UINT64_C(0x686f6f6b63616368),
            storage,
            packets))
    {
        return false;
    }
    const bool matched = ComparePacketCorpus(
        fixture, packets, "cached-row-hooks");
    const uint32_t P = fixture.System.Params.Staircase +
        fixture.System.Params.DenseRows + fixture.System.Params.HeavyRows;
    wirehair_v2::PacketRowRuntime runtime;
    std::vector<uint8_t> cached_output;
    wirehair_v2::ResetTinyDirectSolveCountersForTesting();
    wirehair_v2::SetTinyDirectSolveForTesting(true);
    const WirehairResult cached_result = runtime.Initialize(
            fixture.K, P, fixture.Config.MixCount) ?
        wirehair_v2::SolvePrecodeSystemForValidatedSystemWithRuntime(
            fixture.System,
            fixture.Config,
            runtime,
            packets,
            fixture.BlockBytes,
            cached_output) : Wirehair_InvalidInput;
    wirehair_v2::SetTinyDirectSolveForTesting(false);
    const bool cached_verified = cached_result == Wirehair_Success &&
        wirehair_v2::VerifyPrecodeSolution(
            fixture.System,
            fixture.Config,
            packets,
            cached_output.data(),
            fixture.BlockBytes);
    wirehair_v2::SetOddPacketPeelSeedXorForTesting(0u);
    wirehair_v2::SetPacketRowSeedAvalancheForTesting(false);
    (void)wirehair_v2::SetPacketRowSeedMultiplierForTesting(1u);
    if (!matched || cached_result != Wirehair_Success ||
        wirehair_v2::TinyDirectSolveAttemptsForTesting() != 1u ||
        wirehair_v2::TinyDirectSolveCompletionsForTesting() != 1u ||
        !cached_verified)
    {
        std::fprintf(stderr,
            "tiny-direct: cached packet-hook runtime failed result=%d\n",
            (int)cached_result);
        return false;
    }
    std::printf("tiny-direct cached packet-row hooks: PASS\n");
    return true;
}

bool CheckDuplicateConflictAndResume()
{
    Fixture fixture;
    if (!MakeFixture(
            64u, 17u, wirehair_v2::DenseAnchorLayout::Two07, fixture))
    {
        return false;
    }

    std::vector<uint8_t> expected;
    if (!SolveWithMode(
            fixture, false, fixture.Packets, expected))
    {
        return false;
    }

    std::vector<wirehair_v2::SolvePacket> duplicate = fixture.Packets;
    duplicate.push_back(fixture.Packets[0]);
    std::vector<uint8_t> output;
    wirehair_v2::ResetTinyDirectSolveCountersForTesting();
    if (!SolveWithMode(fixture, true, duplicate, output) ||
        output != expected ||
        wirehair_v2::TinyDirectSolveCompletionsForTesting() != 1u)
    {
        std::fprintf(stderr, "tiny-direct: exact surplus duplicate failed\n");
        return false;
    }

    std::vector<uint8_t> corrupt(
        fixture.Message.begin(),
        fixture.Message.begin() + fixture.BlockBytes);
    corrupt[0] ^= 1u;
    duplicate.back().Data = corrupt.data();
    output.assign(11u, 0xa5u);
    const std::vector<uint8_t> output_before = output;
    wirehair_v2::PrecodeSolveResumeState conflict_resume;
    wirehair_v2::ResetTinyDirectSolveCountersForTesting();
    wirehair_v2::SetTinyDirectSolveForTesting(true);
    const WirehairResult conflict_result = wirehair_v2::SolvePrecodeSystem(
        fixture.System,
        fixture.Config,
        duplicate,
        fixture.BlockBytes,
        output,
        nullptr,
        &conflict_resume);
    wirehair_v2::SetTinyDirectSolveForTesting(false);
    if (conflict_result != Wirehair_Error || output != output_before ||
        conflict_resume.Active ||
        wirehair_v2::TinyDirectSolveAttemptsForTesting() != 1u ||
        wirehair_v2::TinyDirectSolveCompletionsForTesting() != 1u ||
        wirehair_v2::TinyDirectSolveFallbacksForTesting() != 0u)
    {
        std::fprintf(stderr,
            "tiny-direct: surplus conflict classification failed result=%d\n",
            (int)conflict_result);
        return false;
    }

    std::vector<wirehair_v2::SolvePacket> deficient(
        fixture.K, fixture.Packets[0]);
    output.assign(7u, 0x6du);
    const std::vector<uint8_t> deficient_before = output;
    wirehair_v2::PrecodeSolveStats direct_stats;
    wirehair_v2::ResetTinyDirectSolveCountersForTesting();
    wirehair_v2::SetTinyDirectSolveForTesting(true);
    const WirehairResult direct_need_more = wirehair_v2::SolvePrecodeSystem(
        fixture.System,
        fixture.Config,
        deficient,
        fixture.BlockBytes,
        output,
        &direct_stats);
    wirehair_v2::SetTinyDirectSolveForTesting(false);
    const uint32_t L = fixture.K + fixture.System.Params.Staircase +
        fixture.System.Params.DenseRows + fixture.System.Params.HeavyRows;
    if (direct_need_more != Wirehair_NeedMore ||
        output != deficient_before ||
        direct_stats.InactivatedColumns != L ||
        direct_stats.ResidualRank >= L ||
        wirehair_v2::TinyDirectSolveCompletionsForTesting() != 1u ||
        wirehair_v2::TinyDirectSolveFallbacksForTesting() != 0u)
    {
        std::fprintf(stderr,
            "tiny-direct: direct NeedMore classification failed result=%d\n",
            (int)direct_need_more);
        return false;
    }

    wirehair_v2::PrecodeSolveResumeState resume;
    wirehair_v2::PrecodeSolveStats resume_stats;
    output = deficient_before;
    wirehair_v2::ResetTinyDirectSolveCountersForTesting();
    wirehair_v2::SetTinyDirectSolveForTesting(true);
    const WirehairResult resumed_need_more = wirehair_v2::SolvePrecodeSystem(
        fixture.System,
        fixture.Config,
        deficient,
        fixture.BlockBytes,
        output,
        &resume_stats,
        &resume);
    wirehair_v2::SetTinyDirectSolveForTesting(false);
    if (resumed_need_more != Wirehair_NeedMore || !resume.Active ||
        output != deficient_before ||
        resume_stats.PeeledColumns + resume_stats.InactivatedColumns != L ||
        wirehair_v2::TinyDirectSolveAttemptsForTesting() != 1u ||
        wirehair_v2::TinyDirectSolveCompletionsForTesting() != 0u ||
        wirehair_v2::TinyDirectSolveFallbacksForTesting() != 1u)
    {
        std::fprintf(stderr,
            "tiny-direct: resumable NeedMore did not fall through result=%d\n",
            (int)resumed_need_more);
        return false;
    }

    WirehairResult completion_result = Wirehair_NeedMore;
    for (uint32_t id = 1u;
         id < fixture.K && completion_result == Wirehair_NeedMore;
         ++id)
    {
        completion_result = wirehair_v2::ResumePrecodeSystem(
            fixture.System,
            fixture.Config,
            fixture.Packets[id].BlockId,
            fixture.Packets[id].Data,
            fixture.BlockBytes,
            resume,
            output,
            &resume_stats,
            true);
    }
    if (completion_result != Wirehair_Success || resume.Active ||
        output != expected ||
        !wirehair_v2::VerifyPrecodeSolution(
            fixture.System,
            fixture.Config,
            fixture.Packets,
            output.data(),
            fixture.BlockBytes))
    {
        std::fprintf(stderr,
            "tiny-direct: general resume completion failed result=%d\n",
            (int)completion_result);
        return false;
    }

    std::printf("tiny-direct duplicate/conflict/resume: PASS\n");
    return true;
}

bool CheckAllocationFallbacks()
{
    Fixture fixture;
    if (!MakeFixture(
            17u, 13u, wirehair_v2::DenseAnchorLayout::Disabled, fixture))
    {
        return false;
    }
    std::vector<uint8_t> expected;
    wirehair_v2::PrecodeSolveStats expected_stats;
    if (!SolveWithMode(
            fixture,
            false,
            fixture.Packets,
            expected,
            &expected_stats))
    {
        return false;
    }

    // Checkpoint zero is the augmented arena; checkpoint one is the
    // transactional output allocation after full rank.
    for (uint32_t exception = 0u; exception < 2u; ++exception)
    {
        for (int64_t checkpoint = 0; checkpoint < 2; ++checkpoint)
        {
            std::vector<uint8_t> output(9u, 0x3cu);
            wirehair_v2::PrecodeSolveStats stats;
            wirehair_v2::ResetTinyDirectSolveCountersForTesting();
            wirehair_v2::SetTinyDirectSolveAllocationFailureForTesting(
                checkpoint, exception != 0u);
            wirehair_v2::SetTinyDirectSolveForTesting(true);
            const WirehairResult result = wirehair_v2::SolvePrecodeSystem(
                fixture.System,
                fixture.Config,
                fixture.Packets,
                fixture.BlockBytes,
                output,
                &stats);
            wirehair_v2::SetTinyDirectSolveForTesting(false);
            wirehair_v2::SetTinyDirectSolveAllocationFailureForTesting(
                -1, false);
            if (result != Wirehair_Success || output != expected ||
                !SameUntimedStats(stats, expected_stats) ||
                wirehair_v2::TinyDirectSolveAttemptsForTesting() != 1u ||
                wirehair_v2::TinyDirectSolveCompletionsForTesting() != 0u ||
                wirehair_v2::TinyDirectSolveFallbacksForTesting() != 1u)
            {
                std::fprintf(stderr,
                    "tiny-direct: allocation fallback failed exception=%u "
                    "checkpoint=%lld result=%d\n",
                    exception,
                    (long long)checkpoint,
                    (int)result);
                return false;
            }
        }
    }
    std::printf("tiny-direct bad_alloc/length_error fallbacks: PASS\n");
    return true;
}

bool CheckEligibilityBoundaries()
{
    struct Boundary
    {
        uint32_t K;
        uint32_t BlockBytes;
        uint32_t Surplus;
        const char* Name;
    };
    static const Boundary kDeclines[] = {
        {101u, 13u, 0u, "K"},
        {2u, 4097u, 0u, "block-bytes"},
        {64u, 13u, 5u, "surplus"}
    };
    for (const Boundary& boundary : kDeclines)
    {
        Fixture fixture;
        if (!MakeFixture(
                boundary.K,
                boundary.BlockBytes,
                wirehair_v2::DenseAnchorLayout::Disabled,
                fixture))
        {
            return false;
        }
        std::vector<wirehair_v2::SolvePacket> packets = fixture.Packets;
        for (uint32_t i = 0u; i < boundary.Surplus; ++i) {
            packets.push_back(fixture.Packets[0]);
        }
        std::vector<uint8_t> output;
        wirehair_v2::ResetTinyDirectSolveCountersForTesting();
        if (!SolveWithMode(fixture, true, packets, output) ||
            wirehair_v2::TinyDirectSolveAttemptsForTesting() != 0u)
        {
            std::fprintf(stderr,
                "tiny-direct: ineligible %s did not decline\n",
                boundary.Name);
            return false;
        }
    }

    Fixture exact_surplus;
    if (!MakeFixture(
            100u,
            13u,
            wirehair_v2::DenseAnchorLayout::Two07,
            exact_surplus))
    {
        return false;
    }
    for (uint32_t i = 0u; i < 4u; ++i) {
        exact_surplus.Packets.push_back(exact_surplus.Packets[0]);
    }
    std::vector<uint8_t> output;
    wirehair_v2::ResetTinyDirectSolveCountersForTesting();
    if (!SolveWithMode(
            exact_surplus,
            true,
            exact_surplus.Packets,
            output) ||
        wirehair_v2::TinyDirectSolveAttemptsForTesting() != 1u ||
        wirehair_v2::TinyDirectSolveCompletionsForTesting() != 1u)
    {
        std::fprintf(stderr,
            "tiny-direct: exact K+4 surplus boundary was not selected\n");
        return false;
    }

    struct GeometryDecline
    {
        wirehair_v2::PrecodeParams Params;
        uint32_t MixCount;
        const char* Name;
    };
    const wirehair_v2::PrecodeParams certified =
        wirehair_v2::MakeCertifiedParams(
            64u, UINT64_C(0x637573746f6d5331));
    std::vector<GeometryDecline> geometry_declines;
    GeometryDecline decline = {certified, 3u, "S"};
    ++decline.Params.Staircase;
    geometry_declines.push_back(decline);
    decline = GeometryDecline{certified, 3u, "D2"};
    --decline.Params.DenseRows;
    geometry_declines.push_back(decline);
    decline = GeometryDecline{certified, 3u, "H"};
    --decline.Params.HeavyRows;
    geometry_declines.push_back(decline);
    decline = GeometryDecline{certified, 3u, "N1"};
    decline.Params.SourceHits = 1u;
    geometry_declines.push_back(decline);
    decline = GeometryDecline{certified, 2u, "mix"};
    geometry_declines.push_back(decline);
    decline = GeometryDecline{certified, 3u, "heavy-family"};
    decline.Params.HeavyFamily =
        wirehair_v2::HeavyCoefficientFamily::HashedNonzero;
    geometry_declines.push_back(decline);
    decline = GeometryDecline{certified, 3u, "identity"};
    decline.Params.DenseIdentityCorner = true;
    geometry_declines.push_back(decline);
    decline = GeometryDecline{certified, 3u, "four-anchor"};
    decline.Params.DenseAnchors =
        wirehair_v2::DenseAnchorLayout::Four0369;
    geometry_declines.push_back(decline);

    for (const GeometryDecline& geometry : geometry_declines)
    {
        Fixture custom_fixture;
        if (!MakeFixtureFromParams(
                geometry.Params,
                13u,
                custom_fixture,
                geometry.MixCount))
        {
            return false;
        }
        wirehair_v2::ResetTinyDirectSolveCountersForTesting();
        if (!SolveWithMode(
                custom_fixture,
                true,
                custom_fixture.Packets,
                output) ||
            wirehair_v2::TinyDirectSolveAttemptsForTesting() != 0u)
        {
            std::fprintf(stderr,
                "tiny-direct: ineligible %s geometry was not declined\n",
                geometry.Name);
            return false;
        }
    }

    Fixture malformed;
    if (!MakeFixture(
            17u,
            13u,
            wirehair_v2::DenseAnchorLayout::Disabled,
            malformed))
    {
        return false;
    }
    const uint32_t malformed_P = malformed.System.Params.Staircase +
        malformed.System.Params.DenseRows +
        malformed.System.Params.HeavyRows;
    const uint32_t malformed_L = malformed.K + malformed_P;
    malformed.System.StaircaseRows[0].back() = malformed_L;
    wirehair_v2::PacketRowRuntime runtime;
    if (!runtime.Initialize(
            malformed.K, malformed_P, malformed.Config.MixCount))
    {
        return false;
    }
    output.assign(15u, 0xb7u);
    const std::vector<uint8_t> before = output;
    wirehair_v2::PrecodeSolveStats stats;
    stats.PacketRows = UINT32_C(0xa5a5a5a5);
    wirehair_v2::ResetTinyDirectSolveCountersForTesting();
    wirehair_v2::SetTinyDirectSolveForTesting(true);
    const WirehairResult malformed_result =
        wirehair_v2::SolvePrecodeSystemForValidatedSystemWithRuntime(
            malformed.System,
            malformed.Config,
            runtime,
            malformed.Packets,
            malformed.BlockBytes,
            output,
            &stats);
    wirehair_v2::SetTinyDirectSolveForTesting(false);
    if (malformed_result != Wirehair_InvalidInput || output != before ||
        stats.PacketRows != UINT32_C(0xa5a5a5a5) ||
        wirehair_v2::TinyDirectSolveAttemptsForTesting() != 1u ||
        wirehair_v2::TinyDirectSolveCompletionsForTesting() != 1u ||
        wirehair_v2::TinyDirectSolveFallbacksForTesting() != 0u)
    {
        std::fprintf(stderr,
            "tiny-direct: trusted malformed row did not fail safely result=%d\n",
            (int)malformed_result);
        return false;
    }

    std::printf("tiny-direct eligibility/trusted-row boundaries: PASS\n");
    return true;
}

bool TimeSolveInvocation(
    const Fixture& fixture,
    const std::vector<wirehair_v2::SolvePacket>& packets,
    const std::vector<uint8_t>& expected,
    bool direct,
    uint64_t& nanoseconds)
{
    std::vector<uint8_t> output;
    wirehair_v2::ResetTinyDirectSolveCountersForTesting();
    wirehair_v2::SetTinyDirectSolveForTesting(direct);
    const std::chrono::steady_clock::time_point start =
        std::chrono::steady_clock::now();
    const WirehairResult result = wirehair_v2::SolvePrecodeSystem(
        fixture.System,
        fixture.Config,
        packets,
        fixture.BlockBytes,
        output);
    const std::chrono::steady_clock::time_point finish =
        std::chrono::steady_clock::now();
    wirehair_v2::SetTinyDirectSolveForTesting(false);
    const std::chrono::nanoseconds elapsed =
        std::chrono::duration_cast<std::chrono::nanoseconds>(finish - start);
    nanoseconds = elapsed.count() > 0 ? (uint64_t)elapsed.count() : 0u;
    const uint64_t attempts =
        wirehair_v2::TinyDirectSolveAttemptsForTesting();
    const uint64_t completions =
        wirehair_v2::TinyDirectSolveCompletionsForTesting();
    const uint64_t fallbacks =
        wirehair_v2::TinyDirectSolveFallbacksForTesting();
    const bool identity_matches = direct ?
        attempts == 1u && completions == 1u && fallbacks == 0u :
        attempts == 0u && completions == 0u && fallbacks == 0u;
    return result == Wirehair_Success && nanoseconds != 0u &&
        output == expected && identity_matches;
}

bool RunTimingPrefilter(uint32_t repetitions)
{
    static const uint32_t kBlockCounts[] = {
        2u, 3u, 4u, 5u, 6u, 7u, 8u
    };
    static const uint32_t kBlockBytes[] = {
        2u, 64u, 1280u, 4096u
    };
    static const wirehair_v2::DenseAnchorLayout kLayouts[] = {
        wirehair_v2::DenseAnchorLayout::Disabled,
        wirehair_v2::DenseAnchorLayout::Two07
    };
    struct Panel
    {
        const char* Name;
        bool ADirect;
        bool BDirect;
    };
    static const Panel kPanels[] = {
        {"off_off", false, false},
        {"on_on", true, true},
        {"on_off", true, false}
    };

    std::printf(
        "# tiny-direct directional public-SolvePrecodeSystem prefilter; "
        "includes runtime construction and system validation; 8 alternating "
        "panel-specific untimed warmup invocations per cell/panel; measured "
        "times exclude mode setter, counter identity checks, and byte "
        "verification\n");
    for (wirehair_v2::DenseAnchorLayout layout : kLayouts)
    {
        for (uint32_t K : kBlockCounts)
        {
            for (uint32_t block_bytes : kBlockBytes)
            {
                Fixture fixture;
                if (!MakeFixture(K, block_bytes, layout, fixture)) {
                    return false;
                }
                std::vector<uint8_t> expected;
                if (!SolveWithMode(
                        fixture,
                        false,
                        fixture.Packets,
                        expected))
                {
                    return false;
                }
                std::vector<uint8_t> storage;
                std::vector<wirehair_v2::SolvePacket> packets;
                if (!MakeEncodedPacketCorpus(
                        fixture,
                        expected,
                        K + 4u,
                        UINT64_C(0x70726566696c7465) ^
                            ((uint64_t)K << 32) ^ block_bytes ^
                            (uint32_t)layout,
                        storage,
                        packets))
                {
                    return false;
                }

                for (const Panel& panel : kPanels)
                {
                    // Warm the exact arm or pair immediately before its
                    // panel.  AA panels get eight invocations of their one
                    // mode; AB gets four of each mode in alternating order.
                    for (uint32_t warmup = 0u; warmup < 8u; ++warmup)
                    {
                        uint64_t ignored_nanoseconds = 0u;
                        const bool warmup_direct = (warmup & 1u) == 0u ?
                            panel.ADirect : panel.BDirect;
                        if (!TimeSolveInvocation(
                                fixture,
                                packets,
                                expected,
                                warmup_direct,
                                ignored_nanoseconds))
                        {
                            return false;
                        }
                    }
                    for (uint32_t repetition = 0u;
                         repetition < repetitions;
                         ++repetition)
                    {
                        const bool abba = (repetition & 1u) == 0u;
                        const bool order[4] = {
                            abba ? panel.ADirect : panel.BDirect,
                            abba ? panel.BDirect : panel.ADirect,
                            abba ? panel.BDirect : panel.ADirect,
                            abba ? panel.ADirect : panel.BDirect
                        };
                        uint64_t a_nanoseconds = 0u;
                        uint64_t b_nanoseconds = 0u;
                        for (uint32_t invocation = 0u;
                             invocation < 4u;
                             ++invocation)
                        {
                            uint64_t elapsed = 0u;
                            if (!TimeSolveInvocation(
                                    fixture,
                                    packets,
                                    expected,
                                    order[invocation],
                                    elapsed))
                            {
                                return false;
                            }
                            const bool is_a = abba ?
                                invocation == 0u || invocation == 3u :
                                invocation == 1u || invocation == 2u;
                            if (is_a) a_nanoseconds += elapsed;
                            else b_nanoseconds += elapsed;
                        }
                        std::printf(
                            "tiny_direct_prefilter,layout=%s,K=%u,bb=%u,"
                            "panel=%s,rep=%u,order=%s,a_direct=%u,b_direct=%u,"
                            "a_ns=%llu,b_ns=%llu\n",
                            layout == wirehair_v2::DenseAnchorLayout::Two07 ?
                                "two07" : "disabled",
                            K,
                            block_bytes,
                            panel.Name,
                            repetition,
                            abba ? "ABBA" : "BAAB",
                            panel.ADirect ? 1u : 0u,
                            panel.BDirect ? 1u : 0u,
                            (unsigned long long)a_nanoseconds,
                            (unsigned long long)b_nanoseconds);
                    }
                }
            }
        }
    }
    return true;
}

} // namespace

int main(int argc, char** argv)
{
    TinyDirectStateReset reset;
    if (argc == 3 &&
        std::strcmp(argv[1], "--timing-prefilter") == 0)
    {
        char* end = nullptr;
        const unsigned long parsed = std::strtoul(argv[2], &end, 10);
        if (parsed == 0u || parsed > 10000u || !end || *end != '\0') {
            return 2;
        }
        return RunTimingPrefilter((uint32_t)parsed) ? 0 : 1;
    }
    if (argc != 1) {
        return 2;
    }
    const bool ok = CheckForcedOracleMatrix() &&
        CheckRandomDecoderShapes() &&
        CheckPacketRowHookAgreement() &&
        CheckDuplicateConflictAndResume() &&
        CheckAllocationFallbacks() &&
        CheckEligibilityBoundaries();
    return ok ? 0 : 1;
}
