#include "WirehairV2PrecodeDecode.h"
#include "WirehairV2PrecodeEncode.h"
#include "../gf256.h"

#include <wirehair/wirehair.h>

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <limits>
#include <memory>
#include <thread>
#include <vector>

namespace {

using namespace wirehair_v2;

uint64_t NextRandom(uint64_t& state)
{
    state ^= state >> 12;
    state ^= state << 25;
    state ^= state >> 27;
    return state * UINT64_C(0x2545f4914f6cdd1d);
}

std::vector<uint8_t> MakeMessage(uint64_t bytes, uint64_t seed)
{
    std::vector<uint8_t> message((size_t)bytes);
    for (size_t i = 0u; i < message.size(); ++i) {
        message[i] = (uint8_t)NextRandom(seed);
    }
    return message;
}

bool SamePeelingCodec(const PeelingCodec& a, const PeelingCodec& b)
{
    return a.Solver == b.Solver && a.Structure == b.Structure &&
        a.Family == b.Family && a.MinDegree == b.MinDegree &&
        a.MaxDegree == b.MaxDegree &&
        a.SolverCandidateLimit == b.SolverCandidateLimit &&
        a.Degree1Mass == b.Degree1Mass &&
        a.Degree2Mass == b.Degree2Mass && a.RobustC == b.RobustC &&
        a.RobustDelta == b.RobustDelta &&
        a.FullyRandomRows == b.FullyRandomRows &&
        a.UseWirehairRowDistribution == b.UseWirehairRowDistribution;
}

bool SameProfile(const SeedProfile& a, const SeedProfile& b)
{
    return a.BlockCount == b.BlockCount && a.BlockBytes == b.BlockBytes &&
        a.Policy.Solver == b.Policy.Solver &&
        a.Policy.Structure == b.Policy.Structure &&
        a.Policy.ByteClass == b.Policy.ByteClass &&
        a.Policy.CountBand == b.Policy.CountBand &&
        SamePeelingCodec(a.Policy.Codec, b.Policy.Codec) &&
        a.DenseCount == b.DenseCount && a.PeelSeed == b.PeelSeed &&
        a.DenseSeed == b.DenseSeed &&
        a.PeelSeedBucket == b.PeelSeedBucket &&
        a.UsedPeelFixup == b.UsedPeelFixup &&
        a.UsedDenseFixup == b.UsedDenseFixup &&
        a.V2SeedSelected == b.V2SeedSelected &&
        a.V2SeedAttempt == b.V2SeedAttempt &&
        a.V2PrecodeContractVersion == b.V2PrecodeContractVersion &&
        a.V2PacketRowContractVersion == b.V2PacketRowContractVersion &&
        a.V2StaircaseCount == b.V2StaircaseCount &&
        a.V2DenseRowCount == b.V2DenseRowCount &&
        a.V2HeavyRowCount == b.V2HeavyRowCount &&
        a.V2SourceHits == b.V2SourceHits &&
        a.V2PrecodeSeed == b.V2PrecodeSeed &&
        a.V2PacketPeelSeed == b.V2PacketPeelSeed &&
        a.V2RecoveryMixCount == b.V2RecoveryMixCount &&
        a.V2DenseIdentityCorner == b.V2DenseIdentityCorner &&
        a.V2PrecodeSeedSalt == b.V2PrecodeSeedSalt &&
        a.V2RecoveryRowSeedSalt == b.V2RecoveryRowSeedSalt &&
        a.Tuned == b.Tuned &&
        a.TuningResidualMean == b.TuningResidualMean &&
        a.TuningResidualColumns == b.TuningResidualColumns &&
        a.TuningXorCost == b.TuningXorCost &&
        a.TuningCandidatesRequested == b.TuningCandidatesRequested &&
        a.TuningCandidatesUnique == b.TuningCandidatesUnique &&
        a.TuningCandidatesCompleted == b.TuningCandidatesCompleted &&
        a.TuningTrials == b.TuningTrials;
}

bool SameOptions(
    const MessagePrecodeEncoderOptions& a,
    const MessagePrecodeEncoderOptions& b)
{
    return a.RecoveryMixCount == b.RecoveryMixCount &&
        a.DenseIdentityCorner == b.DenseIdentityCorner &&
        a.PrecodeSeedSalt == b.PrecodeSeedSalt &&
        a.RecoveryRowSeedSalt == b.RecoveryRowSeedSalt &&
        a.CacheSystematicSource == b.CacheSystematicSource &&
        a.CacheReceivedSystematicPackets ==
            b.CacheReceivedSystematicPackets;
}

bool SameParams(const PrecodeParams& a, const PrecodeParams& b)
{
    return a.BlockCount == b.BlockCount && a.Staircase == b.Staircase &&
        a.DenseRows == b.DenseRows && a.HeavyRows == b.HeavyRows &&
        a.SourceHits == b.SourceHits &&
        a.DenseIdentityCorner == b.DenseIdentityCorner &&
        a.HeavyFamily == b.HeavyFamily &&
        a.DenseAnchors == b.DenseAnchors && a.Seed == b.Seed;
}

bool EmptyPlanSolveStats(
    const PrecodeSolveStats& stats,
    uint32_t packet_seed_attempt)
{
    return stats.PacketRows == 0u && stats.PeeledColumns == 0u &&
        stats.InactivatedColumns == 0u && stats.ResidualRows == 0u &&
        stats.ResidualRank == 0u && stats.BinaryResidualRank == 0u &&
        stats.BinaryRowReferences == 0u &&
        stats.BinaryRowStorageBytes == 0u &&
        stats.BinaryAdjacencyStorageBytes == 0u &&
        stats.BinaryRowStorageAllocations == 0u &&
        stats.BinaryAdjacencyStorageAllocations == 0u &&
        stats.BlockXors == 0u && stats.BlockMulAdds == 0u &&
        stats.BuildNanoseconds == 0u && stats.PeelNanoseconds == 0u &&
        stats.ProjectNanoseconds == 0u &&
        stats.ResidualNanoseconds == 0u &&
        stats.BackSubNanoseconds == 0u &&
        stats.PacketSeedAttempt == packet_seed_attempt;
}

bool CompareOutputs(
    const MessagePrecodeEncoder& candidate,
    const MessagePrecodeEncoder& legacy,
    uint32_t K,
    uint32_t B)
{
    std::vector<uint8_t> candidate_block(B, uint8_t{0xa5});
    std::vector<uint8_t> legacy_block(B, uint8_t{0x5a});
    for (uint32_t id = 0u; id < 3u * K + 9u; ++id)
    {
        std::fill(candidate_block.begin(), candidate_block.end(), uint8_t{0xa5});
        std::fill(legacy_block.begin(), legacy_block.end(), uint8_t{0x5a});
        uint32_t candidate_bytes = UINT32_MAX;
        uint32_t legacy_bytes = UINT32_MAX;
        uint64_t candidate_ops = UINT64_MAX;
        uint64_t legacy_ops = UINT64_MAX;
        const WirehairResult candidate_result = candidate.EncodeResult(
            id,
            candidate_block.data(),
            B,
            &candidate_bytes,
            &candidate_ops);
        const WirehairResult legacy_result = legacy.EncodeResult(
            id,
            legacy_block.data(),
            B,
            &legacy_bytes,
            &legacy_ops);
        if (candidate_result != legacy_result ||
            candidate_result != Wirehair_Success ||
            candidate_bytes != legacy_bytes || candidate_ops != legacy_ops ||
            std::memcmp(
                candidate_block.data(),
                legacy_block.data(),
                candidate_bytes) != 0)
        {
            std::fprintf(stderr, "generator plan output mismatch id=%u\n", id);
            return false;
        }
    }
    return true;
}

bool CompareEncoderState(
    const MessagePrecodeEncoder& candidate,
    const MessagePrecodeEncoder& legacy)
{
    if (!candidate.IsInitialized() || !legacy.IsInitialized() ||
        candidate.MessageBytes() != legacy.MessageBytes() ||
        candidate.SourceBlockCount() != legacy.SourceBlockCount() ||
        candidate.BlockBytes() != legacy.BlockBytes() ||
        !SameProfile(candidate.Profile(), legacy.Profile()) ||
        !SameOptions(candidate.Options(), legacy.Options()) ||
        candidate.HasSystematicSourceCache() !=
            legacy.HasSystematicSourceCache() ||
        candidate.SystematicSourceCacheBytes() !=
            legacy.SystematicSourceCacheBytes() ||
        candidate.BlockEncoder().HasCompleteSystem() ||
        legacy.BlockEncoder().HasCompleteSystem() ||
        !SameParams(
            candidate.BlockEncoder().System().Params,
            legacy.BlockEncoder().System().Params) ||
        candidate.BlockEncoder().RecoveryRowSeed() !=
            legacy.BlockEncoder().RecoveryRowSeed() ||
        candidate.BlockEncoder().RecoveryMixCount() !=
            legacy.BlockEncoder().RecoveryMixCount() ||
        candidate.SolveStats().PacketSeedAttempt !=
            legacy.SolveStats().PacketSeedAttempt)
    {
        return false;
    }
    const uint32_t K = candidate.SourceBlockCount();
    const uint32_t B = candidate.BlockBytes();
    const PrecodeParams& params = candidate.BlockEncoder().System().Params;
    const size_t intermediate_bytes =
        (size_t)(K + params.Staircase + params.DenseRows + params.HeavyRows) * B;
    return candidate.IntermediateBlocks() && legacy.IntermediateBlocks() &&
        std::memcmp(
            candidate.IntermediateBlocks(),
            legacy.IntermediateBlocks(),
            intermediate_bytes) == 0 &&
        CompareOutputs(candidate, legacy, K, B);
}

bool PrepareSelected(
    const std::vector<uint8_t>& message,
    uint32_t B,
    const MessagePrecodeEncoderOptions& requested_options,
    MessagePrecodeEncoder& encoder)
{
    return encoder.InitializeResult(
        message.data(),
        message.size(),
        B,
        nullptr,
        &requested_options) == Wirehair_Success;
}

bool CheckSemanticCell(
    uint32_t K,
    uint32_t B,
    uint32_t tail,
    bool& saw_nonzero_attempt)
{
    const uint64_t message_bytes = (uint64_t)(K - 1u) * B + tail;
    std::vector<uint8_t> message = MakeMessage(
        message_bytes,
        UINT64_C(0x243f6a8885a308d3) ^ ((uint64_t)K << 32) ^ B ^ tail);
    MessagePrecodeEncoderOptions requested_options;
    MessagePrecodeEncoder legacy;
    if (!PrepareSelected(message, B, requested_options, legacy)) {
        return false;
    }
    saw_nonzero_attempt = saw_nonzero_attempt ||
        legacy.Profile().V2SeedAttempt != 0u;

    std::shared_ptr<const SystematicGeneratorPlan> plan;
    if (SystematicGeneratorPlan::Build(
            legacy.Profile(), &legacy.Options(), plan) != Wirehair_Success ||
        !plan || plan->Key().SourceBlockCount != K ||
        plan->Key().BlockBytes != B ||
        plan->Key().PacketSeedAttempt != legacy.Profile().V2SeedAttempt ||
        plan->BuildStats().GeneratorBytes !=
            (uint64_t)plan->ColumnBlockCount() * K ||
        plan->BuildStats().Solve.PacketSeedAttempt !=
            legacy.Profile().V2SeedAttempt ||
        !SameParams(
            plan->SystemParams(), legacy.BlockEncoder().System().Params))
    {
        return false;
    }
    PrecodeSystem expected_system;
    uint64_t expected_fingerprint0 = 0u;
    uint64_t expected_fingerprint1 = 0u;
    if (!BuildPrecodeSystem(plan->SystemParams(), expected_system)) {
        return false;
    }
    internal::PrecodeSystemFingerprint(
        expected_system, expected_fingerprint0, expected_fingerprint1);
    const PacketRowEquationIdentity expected_packet_identity =
        internal::PacketRowEquationIdentitySnapshot();
    if (plan->GraphFingerprint0() != expected_fingerprint0 ||
        plan->GraphFingerprint1() != expected_fingerprint1 ||
        plan->PacketEquationIdentity().BlockIdMultiplier !=
            expected_packet_identity.BlockIdMultiplier ||
        plan->PacketEquationIdentity().BlockIdAvalanche !=
            expected_packet_identity.BlockIdAvalanche ||
        plan->PacketEquationIdentity().OddPeelSeedXor !=
            expected_packet_identity.OddPeelSeedXor ||
        !plan->PacketRuntime().IsValidFor(
            K,
            plan->ColumnBlockCount() - K,
            plan->PacketConfig().MixCount))
    {
        return false;
    }

    MessagePrecodeEncoder candidate;
    bool used_plan = false;
    SystematicGeneratorPlanReplayStats replay_stats;
    if (candidate.InitializeResultWithSystematicGeneratorPlanForTesting(
            plan.get(),
            message.data(),
            message.size(),
            B,
            &legacy.Profile(),
            &legacy.Options(),
            &used_plan,
            &replay_stats) != Wirehair_Success ||
        !used_plan || replay_stats.SourceBlockCount != K ||
        replay_stats.ColumnBlockCount != plan->ColumnBlockCount() ||
        replay_stats.BlockBytes != B ||
        replay_stats.OutputBytes !=
            (uint64_t)plan->ColumnBlockCount() * B ||
        replay_stats.LinearCombinations != plan->ColumnBlockCount() ||
        replay_stats.SourceTerms !=
            (uint64_t)plan->ColumnBlockCount() * K ||
        !EmptyPlanSolveStats(
            candidate.SolveStats(), legacy.Profile().V2SeedAttempt) ||
        !CompareEncoderState(candidate, legacy))
    {
        return false;
    }

    // Both encoders own their complete intermediate values after return.
    std::fill(message.begin(), message.end(), uint8_t{0x69});
    if (!CompareEncoderState(candidate, legacy)) {
        return false;
    }

    // The initialized encoder retains no reference to the caller-owned plan.
    plan.reset();
    return CompareEncoderState(candidate, legacy);
}

bool RunSemanticGrid()
{
    bool saw_nonzero_attempt = false;
    // Every 64-byte-aligned eligible width, with one-byte, near-full, and
    // complete final source blocks: 7 * 17 * 3 = 357 message shapes.
    for (uint32_t K = 2u; K <= 8u; ++K)
    {
        for (uint32_t B = 1024u; B <= 2048u; B += 64u)
        {
            const uint32_t tails[] = { 1u, B - 7u, B };
            for (uint32_t tail : tails)
            {
                if (!CheckSemanticCell(
                        K, B, tail, saw_nonzero_attempt)) {
                    std::fprintf(
                        stderr,
                        "generator plan semantic cell failed K=%u B=%u tail=%u\n",
                        K,
                        B,
                        tail);
                    return false;
                }
            }
        }
    }
    if (!saw_nonzero_attempt) {
        std::fprintf(
            stderr,
            "generator plan semantic grid missed nonzero seed attempts\n");
        return false;
    }
    return true;
}

bool CompareFallback(
    const SystematicGeneratorPlan* plan,
    const std::vector<uint8_t>& message,
    uint32_t B,
    const SeedProfile* profile,
    const MessagePrecodeEncoderOptions* options)
{
    MessagePrecodeEncoder legacy;
    if (legacy.InitializeResult(
            message.data(), message.size(), B, profile, options) !=
        Wirehair_Success)
    {
        return false;
    }
    MessagePrecodeEncoder candidate;
    bool used_plan = true;
    SystematicGeneratorPlanReplayStats replay_stats;
    replay_stats.OutputBytes = UINT64_MAX;
    if (candidate.InitializeResultWithSystematicGeneratorPlanForTesting(
            plan,
            message.data(),
            message.size(),
            B,
            profile,
            options,
            &used_plan,
            &replay_stats) != Wirehair_Success ||
        used_plan || replay_stats.OutputBytes != 0u)
    {
        return false;
    }
    return CompareEncoderState(candidate, legacy);
}

bool RunUnsupportedGFNIFallback()
{
    const uint32_t K = 8u;
    const uint32_t B = 1280u;
    std::vector<uint8_t> message = MakeMessage(
        (uint64_t)K * B - 3u, UINT64_C(0x9e3779b97f4a7c15));
    MessagePrecodeEncoderOptions options;
    MessagePrecodeEncoder selected;
    if (!PrepareSelected(message, B, options, selected)) {
        return false;
    }
    std::shared_ptr<const SystematicGeneratorPlan> plan;
    if (SystematicGeneratorPlan::Build(
            selected.Profile(), &selected.Options(), plan) !=
            Wirehair_UnsupportedPlatform ||
        plan)
    {
        return false;
    }
    return CompareFallback(
        nullptr,
        message,
        B,
        &selected.Profile(),
        &selected.Options());
}

bool RunMismatchAndFallback()
{
    const uint32_t K = 8u;
    const uint32_t B = 1280u;
    std::vector<uint8_t> message = MakeMessage(
        (uint64_t)K * B - 11u, UINT64_C(0x13198a2e03707344));
    MessagePrecodeEncoderOptions options;
    MessagePrecodeEncoder selected;
    if (!PrepareSelected(message, B, options, selected)) {
        return false;
    }
    std::shared_ptr<const SystematicGeneratorPlan> plan;
    if (SystematicGeneratorPlan::Build(
            selected.Profile(), &selected.Options(), plan) !=
            Wirehair_Success ||
        !plan)
    {
        return false;
    }

    if (!CompareFallback(
            nullptr,
            message,
            B,
            &selected.Profile(),
            &selected.Options()))
    {
        return false;
    }

    const uint32_t mismatch_counts[] = { 7u, 9u };
    for (uint32_t mismatch_K : mismatch_counts)
    {
        std::vector<uint8_t> other = MakeMessage(
            (uint64_t)mismatch_K * B - 5u,
            UINT64_C(0xa4093822299f31d0) ^ mismatch_K);
        MessagePrecodeEncoder other_selected;
        if (!PrepareSelected(other, B, options, other_selected) ||
            !CompareFallback(
                plan.get(),
                other,
                B,
                &other_selected.Profile(),
                &other_selected.Options()))
        {
            return false;
        }
    }

    const uint32_t mismatch_widths[] = { 1023u, 1344u };
    for (uint32_t other_B : mismatch_widths)
    {
        std::vector<uint8_t> other = MakeMessage(
            (uint64_t)K * other_B - 5u,
            UINT64_C(0x082efa98ec4e6c89) ^ other_B);
        MessagePrecodeEncoder other_selected;
        if (!PrepareSelected(other, other_B, options, other_selected) ||
            !CompareFallback(
                plan.get(),
                other,
                other_B,
                &other_selected.Profile(),
                &other_selected.Options()))
        {
            return false;
        }
    }

    SeedProfile unselected = SelectSeedProfile(K, B);
    if (!CompareFallback(plan.get(), message, B, &unselected, &options)) {
        return false;
    }

    MessagePrecodeEncoderOptions salted;
    salted.PrecodeSeedSalt ^= UINT64_C(0xd6e8feb86659fd93);
    salted.RecoveryRowSeedSalt ^= UINT64_C(0xa0761d6478bd642f);
    MessagePrecodeEncoder salted_selected;
    if (!PrepareSelected(message, B, salted, salted_selected) ||
        !CompareFallback(
            plan.get(),
            message,
            B,
            &salted_selected.Profile(),
            &salted_selected.Options()))
    {
        return false;
    }

    // SeedProfile diagnostics are not equation identity and must not make a
    // semantically identical selected profile miss the plan.
    SeedProfile diagnostic_only = selected.Profile();
    diagnostic_only.Tuned = !diagnostic_only.Tuned;
    diagnostic_only.TuningTrials += 17u;
    MessagePrecodeEncoder diagnostic_legacy;
    if (diagnostic_legacy.InitializeResult(
            message.data(),
            message.size(),
            B,
            &diagnostic_only,
            &selected.Options()) != Wirehair_Success)
    {
        return false;
    }
    MessagePrecodeEncoder diagnostic_candidate;
    bool diagnostic_used = false;
    SystematicGeneratorPlanReplayStats diagnostic_stats;
    if (diagnostic_candidate.
            InitializeResultWithSystematicGeneratorPlanForTesting(
                plan.get(),
                message.data(),
                message.size(),
                B,
                &diagnostic_only,
                &selected.Options(),
                &diagnostic_used,
                &diagnostic_stats) != Wirehair_Success ||
        !diagnostic_used ||
        !CompareEncoderState(diagnostic_candidate, diagnostic_legacy))
    {
        return false;
    }

    if (!SetPacketRowSeedMultiplierForTesting(3u)) {
        return false;
    }
    const bool packet_identity_fallback = CompareFallback(
        plan.get(),
        message,
        B,
        &selected.Profile(),
        &selected.Options());
    if (!SetPacketRowSeedMultiplierForTesting(1u) ||
        !packet_identity_fallback)
    {
        return false;
    }
    return true;
}

bool RunEligibilityAndBuildOOM()
{
    const uint32_t K = 8u;
    const uint32_t B = 1280u;
    const MessagePrecodeEncoderOptions options;
    std::vector<uint8_t> message = MakeMessage(
        (uint64_t)K * B, UINT64_C(0x452821e638d01377));
    MessagePrecodeEncoder selected;
    if (!PrepareSelected(message, B, options, selected)) {
        return false;
    }
    std::shared_ptr<const SystematicGeneratorPlan> plan;
    if (SystematicGeneratorPlan::Build(
            selected.Profile(), &selected.Options(), plan) !=
            Wirehair_Success ||
        !plan)
    {
        return false;
    }

    std::shared_ptr<const SystematicGeneratorPlan> preserved = plan;
    SetAllocationFailureCountdownForTesting(0);
    const WirehairResult unselected_result = SystematicGeneratorPlan::Build(
        SelectSeedProfile(K, B), &options, preserved);
    SetAllocationFailureCountdownForTesting(-1);
    if (unselected_result != Wirehair_InvalidInput ||
        preserved.get() != plan.get())
    {
        return false;
    }

    const struct IneligibleCell { uint32_t K; uint32_t B; } cells[] = {
        { 9u, 1280u }, { 8u, 1023u }, { 8u, 2112u }
    };
    for (const IneligibleCell& cell : cells)
    {
        std::vector<uint8_t> other = MakeMessage(
            (uint64_t)cell.K * cell.B,
            UINT64_C(0xbe5466cf34e90c6c) ^
                ((uint64_t)cell.K << 32) ^ cell.B);
        MessagePrecodeEncoder other_selected;
        if (!PrepareSelected(other, cell.B, options, other_selected)) {
            return false;
        }
        preserved = plan;
        SetAllocationFailureCountdownForTesting(0);
        const WirehairResult result = SystematicGeneratorPlan::Build(
            other_selected.Profile(), &other_selected.Options(), preserved);
        SetAllocationFailureCountdownForTesting(-1);
        if (result != Wirehair_InvalidInput || preserved.get() != plan.get()) {
            return false;
        }
    }

    bool reached_success = false;
    for (int64_t countdown = 0; countdown < 12; ++countdown)
    {
        preserved = plan;
        SetAllocationFailureCountdownForTesting(countdown);
        const WirehairResult result = SystematicGeneratorPlan::Build(
            selected.Profile(), &selected.Options(), preserved);
        SetAllocationFailureCountdownForTesting(-1);
        if (result == Wirehair_Success)
        {
            reached_success = preserved && preserved.get() != plan.get();
            break;
        }
        if (result != Wirehair_OOM || preserved.get() != plan.get()) {
            return false;
        }
    }
    return reached_success;
}

bool RunReplayOOMAndTransaction()
{
    const uint32_t K = 8u;
    const uint32_t B = 1280u;
    std::vector<uint8_t> message = MakeMessage(
        (uint64_t)K * B, UINT64_C(0xc0ac29b7c97c50dd));
    MessagePrecodeEncoderOptions options;
    MessagePrecodeEncoder legacy;
    if (!PrepareSelected(message, B, options, legacy)) {
        return false;
    }
    std::shared_ptr<const SystematicGeneratorPlan> plan;
    if (SystematicGeneratorPlan::Build(
            legacy.Profile(), &legacy.Options(), plan) != Wirehair_Success ||
        !plan)
    {
        return false;
    }

    const test::EncodeAllocationFailureException exceptions[] = {
        test::EncodeAllocationFailureException::BadAlloc,
        test::EncodeAllocationFailureException::LengthError
    };
    for (test::EncodeAllocationFailureException exception : exceptions)
    {
        test::SetEncodeAllocationFailurePointForTesting(
            test::EncodeAllocationFailurePoint::
                GeneratorPlanReplayIntermediate,
            exception);
        MessagePrecodeEncoder candidate;
        bool used_plan = true;
        SystematicGeneratorPlanReplayStats stats;
        stats.OutputBytes = UINT64_MAX;
        const WirehairResult result =
            candidate.InitializeResultWithSystematicGeneratorPlanForTesting(
                plan.get(),
                message.data(),
                message.size(),
                B,
                &legacy.Profile(),
                &legacy.Options(),
                &used_plan,
                &stats);
        test::SetEncodeAllocationFailurePointForTesting(
            test::EncodeAllocationFailurePoint::None,
            test::EncodeAllocationFailureException::BadAlloc);
        if (result != Wirehair_Success || used_plan ||
            stats.OutputBytes != 0u ||
            !CompareEncoderState(candidate, legacy))
        {
            return false;
        }
    }

    MessagePrecodeEncoderOptions cached_options = legacy.Options();
    cached_options.CacheSystematicSource = true;
    MessagePrecodeEncoder cached_legacy;
    if (cached_legacy.InitializeResult(
            message.data(),
            message.size(),
            B,
            &legacy.Profile(),
            &cached_options) != Wirehair_Success)
    {
        return false;
    }
    for (test::EncodeAllocationFailureException exception : exceptions)
    {
        test::SetEncodeAllocationFailurePointForTesting(
            test::EncodeAllocationFailurePoint::
                GeneratorPlanSystematicCache,
            exception);
        MessagePrecodeEncoder cache_fallback;
        bool cache_used = true;
        SystematicGeneratorPlanReplayStats cache_stats;
        cache_stats.OutputBytes = UINT64_MAX;
        const WirehairResult cache_result =
            cache_fallback.
                InitializeResultWithSystematicGeneratorPlanForTesting(
                    plan.get(),
                    message.data(),
                    message.size(),
                    B,
                    &legacy.Profile(),
                    &cached_options,
                    &cache_used,
                    &cache_stats);
        test::SetEncodeAllocationFailurePointForTesting(
            test::EncodeAllocationFailurePoint::None,
            test::EncodeAllocationFailureException::BadAlloc);
        if (cache_result != Wirehair_Success || cache_used ||
            cache_stats.OutputBytes != 0u ||
            !cache_fallback.HasSystematicSourceCache() ||
            !CompareEncoderState(cache_fallback, cached_legacy))
        {
            return false;
        }
    }

    MessagePrecodeEncoder transactional;
    bool initially_used = false;
    SystematicGeneratorPlanReplayStats initial_stats;
    if (transactional.InitializeResultWithSystematicGeneratorPlanForTesting(
            plan.get(),
            message.data(),
            message.size(),
            B,
            &legacy.Profile(),
            &legacy.Options(),
            &initially_used,
            &initial_stats) != Wirehair_Success ||
        !initially_used)
    {
        return false;
    }
    std::vector<uint8_t> before(B);
    std::vector<uint8_t> after(B);
    uint32_t before_bytes = 0u;
    uint32_t after_bytes = 0u;
    if (transactional.EncodeResult(
            K + 3u, before.data(), B, &before_bytes) != Wirehair_Success)
    {
        return false;
    }
    bool failed_used = true;
    SystematicGeneratorPlanReplayStats failed_stats;
    failed_stats.OutputBytes = UINT64_MAX;
    SetAllocationFailureCountdownForTesting(0);
    const WirehairResult failed_result =
        transactional.InitializeResultWithSystematicGeneratorPlanForTesting(
            plan.get(),
            message.data(),
            message.size(),
            B,
            &legacy.Profile(),
            &legacy.Options(),
            &failed_used,
            &failed_stats);
    SetAllocationFailureCountdownForTesting(-1);
    if (failed_result != Wirehair_OOM || !failed_used ||
        failed_stats.OutputBytes != UINT64_MAX ||
        transactional.EncodeResult(
            K + 3u, after.data(), B, &after_bytes) != Wirehair_Success ||
        before_bytes != after_bytes || before != after)
    {
        return false;
    }
    return true;
}

bool RunInputAlias()
{
    const uint32_t K = 8u;
    const uint32_t B = 1280u;
    std::vector<uint8_t> message = MakeMessage(
        (uint64_t)K * B, UINT64_C(0x3f84d5b5b5470917));
    MessagePrecodeEncoderOptions options;
    MessagePrecodeEncoder encoder;
    if (!PrepareSelected(message, B, options, encoder)) {
        return false;
    }
    std::shared_ptr<const SystematicGeneratorPlan> plan;
    if (SystematicGeneratorPlan::Build(
            encoder.Profile(), &encoder.Options(), plan) != Wirehair_Success ||
        !plan)
    {
        return false;
    }
    const SeedProfile profile = encoder.Profile();
    const MessagePrecodeEncoderOptions selected_options = encoder.Options();
    const uint8_t* const alias = encoder.IntermediateBlocks();
    std::vector<uint8_t> copied(alias, alias + (size_t)K * B);
    MessagePrecodeEncoder legacy;
    if (legacy.InitializeResult(
            copied.data(),
            copied.size(),
            B,
            &profile,
            &selected_options) != Wirehair_Success)
    {
        return false;
    }
    bool used = false;
    SystematicGeneratorPlanReplayStats stats;
    if (encoder.InitializeResultWithSystematicGeneratorPlanForTesting(
            plan.get(),
            alias,
            (uint64_t)K * B,
            B,
            &profile,
            &selected_options,
            &used,
            &stats) != Wirehair_Success ||
        !used)
    {
        return false;
    }
    return CompareEncoderState(encoder, legacy);
}

bool RunConcurrentReuse()
{
    const uint32_t K = 8u;
    const uint32_t B = 1280u;
    const uint64_t message_bytes = (uint64_t)K * B - 13u;
    MessagePrecodeEncoderOptions options;
    std::vector<uint8_t> first = MakeMessage(
        message_bytes, UINT64_C(0x9216d5d98979fb1b));
    MessagePrecodeEncoder selected;
    if (!PrepareSelected(first, B, options, selected)) {
        return false;
    }
    const SeedProfile profile = selected.Profile();
    const MessagePrecodeEncoderOptions selected_options = selected.Options();
    std::shared_ptr<const SystematicGeneratorPlan> plan;
    if (SystematicGeneratorPlan::Build(
            profile, &selected_options, plan) != Wirehair_Success || !plan)
    {
        return false;
    }

    std::atomic<uint32_t> failures(0u);
    std::vector<std::thread> threads;
    for (uint32_t worker = 0u; worker < 8u; ++worker)
    {
        threads.push_back(std::thread([&, worker]() {
            std::vector<uint8_t> message = MakeMessage(
                message_bytes,
                UINT64_C(0xd1310ba698dfb5ac) ^
                    ((uint64_t)worker << 40));
            MessagePrecodeEncoder candidate;
            MessagePrecodeEncoder legacy;
            bool used = false;
            SystematicGeneratorPlanReplayStats stats;
            if (candidate.
                    InitializeResultWithSystematicGeneratorPlanForTesting(
                        plan.get(),
                        message.data(),
                        message.size(),
                        B,
                        &profile,
                        &selected_options,
                        &used,
                        &stats) != Wirehair_Success ||
                !used ||
                legacy.InitializeResult(
                    message.data(),
                    message.size(),
                    B,
                    &profile,
                    &selected_options) != Wirehair_Success ||
                !CompareEncoderState(candidate, legacy))
            {
                failures.fetch_add(1u, std::memory_order_relaxed);
            }
        }));
    }
    for (std::thread& thread : threads) {
        thread.join();
    }
    return failures.load(std::memory_order_relaxed) == 0u;
}

bool RunRecovery()
{
    const uint32_t K = 8u;
    const uint32_t B = 1280u;
    const uint64_t message_bytes = (uint64_t)K * B - 17u;
    std::vector<uint8_t> message = MakeMessage(
        message_bytes, UINT64_C(0x2ffd72dbd01adfb7));
    MessagePrecodeEncoderOptions options;
    MessagePrecodeEncoder selected;
    if (!PrepareSelected(message, B, options, selected)) {
        return false;
    }
    std::shared_ptr<const SystematicGeneratorPlan> plan;
    if (SystematicGeneratorPlan::Build(
            selected.Profile(), &selected.Options(), plan) !=
            Wirehair_Success ||
        !plan)
    {
        return false;
    }
    MessagePrecodeEncoder encoder;
    bool used = false;
    SystematicGeneratorPlanReplayStats stats;
    if (encoder.InitializeResultWithSystematicGeneratorPlanForTesting(
            plan.get(),
            message.data(),
            message.size(),
            B,
            &selected.Profile(),
            &selected.Options(),
            &used,
            &stats) != Wirehair_Success ||
        !used)
    {
        return false;
    }
    MessagePrecodeDecoder decoder;
    if (decoder.InitializeResult(
            message_bytes, B, &encoder.Profile(), &encoder.Options()) !=
        Wirehair_Success)
    {
        return false;
    }
    std::vector<uint8_t> block(B);
    WirehairResult decode_result = Wirehair_NeedMore;
    for (uint32_t id = K;
         decode_result == Wirehair_NeedMore && id < 4u * K + 64u;
         ++id)
    {
        uint32_t bytes = 0u;
        if (encoder.EncodeResult(
                id, block.data(), B, &bytes) != Wirehair_Success)
        {
            return false;
        }
        decode_result = decoder.DecodeResult(id, block.data(), bytes);
    }
    std::vector<uint8_t> recovered((size_t)message_bytes, uint8_t{0xcc});
    return decode_result == Wirehair_Success &&
        decoder.RecoverResult(recovered.data(), recovered.size()) ==
            Wirehair_Success &&
        recovered == message;
}

// Non-promotional same-binary timing screen.  This deliberately lives in the
// focused test executable: the production codec exposes no automatic plan
// cache or policy.  Every comparison starts from one already-selected exact
// profile so seed search is outside all four measured arms.
enum class PreliminaryTimingPhase
{
    FirstK,
    RepairK
};

const char* PreliminaryTimingPhaseName(PreliminaryTimingPhase phase)
{
    return phase == PreliminaryTimingPhase::FirstK ?
        "init_first_k" : "init_repair_k";
}

volatile uint64_t PreliminaryTimingSink = 0u;

uint64_t ConsumePreliminaryTimingOutput(const uint8_t* data, size_t bytes)
{
    uint64_t value = 0u;
    for (size_t offset = 0u; offset < bytes; offset += 31u) {
        value = value * UINT64_C(0x9e3779b97f4a7c15) + data[offset];
    }
    return value;
}

bool IsPreliminaryTimingCell(uint32_t K, uint32_t B)
{
    const bool supported_count =
        (K >= 2u && K <= 8u) || K == 9u || K == 32u || K == 100u;
    return supported_count && B >= 1024u && B <= 2048u &&
        (B & 63u) == 0u;
}

struct PreliminaryTimingContext
{
    uint32_t K = 0u;
    uint32_t B = 0u;
    std::vector<uint8_t> Message;
    SeedProfile Profile = {};
    MessagePrecodeEncoderOptions Options = {};
    std::shared_ptr<const SystematicGeneratorPlan> WarmPlan;
    bool PlanEligible = false;
};

bool InitializePreliminaryPlanArm(
    const PreliminaryTimingContext& context,
    bool cold,
    MessagePrecodeEncoder& encoder)
{
    std::shared_ptr<const SystematicGeneratorPlan> cold_plan;
    const SystematicGeneratorPlan* plan = context.WarmPlan.get();
    if (cold)
    {
        const WirehairResult build_result = SystematicGeneratorPlan::Build(
            context.Profile, &context.Options, cold_plan);
        if (context.PlanEligible)
        {
            if (build_result != Wirehair_Success || !cold_plan) {
                return false;
            }
        }
        else if (build_result != Wirehair_InvalidInput || cold_plan) {
            return false;
        }
        plan = cold_plan.get();
    }

    bool used_plan = !context.PlanEligible;
    SystematicGeneratorPlanReplayStats replay_stats;
    replay_stats.OutputBytes = UINT64_MAX;
    const WirehairResult result =
        encoder.InitializeResultWithSystematicGeneratorPlanForTesting(
            plan,
            context.Message.data(),
            context.Message.size(),
            context.B,
            &context.Profile,
            &context.Options,
            &used_plan,
            &replay_stats);
    if (result != Wirehair_Success || used_plan != context.PlanEligible) {
        return false;
    }
    if (context.PlanEligible)
    {
        return replay_stats.SourceBlockCount == context.K &&
            replay_stats.BlockBytes == context.B &&
            replay_stats.ColumnBlockCount == plan->ColumnBlockCount() &&
            replay_stats.OutputBytes ==
                (uint64_t)plan->ColumnBlockCount() * context.B &&
            EmptyPlanSolveStats(
                encoder.SolveStats(), context.Profile.V2SeedAttempt);
    }
    return replay_stats.SourceBlockCount == 0u &&
        replay_stats.ColumnBlockCount == 0u &&
        replay_stats.BlockBytes == 0u && replay_stats.OutputBytes == 0u;
}

bool EncodePreliminaryWH2Range(
    MessagePrecodeEncoder& encoder,
    const PreliminaryTimingContext& context,
    PreliminaryTimingPhase phase)
{
    const uint32_t first =
        phase == PreliminaryTimingPhase::FirstK ? 0u : context.K;
    std::vector<uint8_t> output(context.B);
    for (uint32_t offset = 0u; offset < context.K; ++offset)
    {
        uint32_t written = 0u;
        if (encoder.EncodeResult(
                first + offset,
                output.data(),
                context.B,
                &written) != Wirehair_Success ||
            written != context.B)
        {
            return false;
        }
        PreliminaryTimingSink ^=
            ConsumePreliminaryTimingOutput(output.data(), written);
    }
    return true;
}

bool RunPreliminaryLegacyAction(
    const PreliminaryTimingContext& context,
    PreliminaryTimingPhase phase)
{
    MessagePrecodeEncoder encoder;
    return encoder.InitializeResult(
               context.Message.data(),
               context.Message.size(),
               context.B,
               &context.Profile,
               &context.Options) == Wirehair_Success &&
        EncodePreliminaryWH2Range(encoder, context, phase);
}

bool RunPreliminaryPlanAction(
    const PreliminaryTimingContext& context,
    PreliminaryTimingPhase phase,
    bool cold)
{
    MessagePrecodeEncoder encoder;
    return InitializePreliminaryPlanArm(context, cold, encoder) &&
        EncodePreliminaryWH2Range(encoder, context, phase);
}

bool RunPreliminaryWH1Action(
    const PreliminaryTimingContext& context,
    PreliminaryTimingPhase phase)
{
    WirehairCodec encoder = wirehair_encoder_create(
        nullptr,
        context.Message.data(),
        context.Message.size(),
        context.B);
    if (!encoder) {
        return false;
    }
    const uint32_t first =
        phase == PreliminaryTimingPhase::FirstK ? 0u : context.K;
    std::vector<uint8_t> output(context.B);
    bool success = true;
    for (uint32_t offset = 0u; offset < context.K; ++offset)
    {
        uint32_t written = 0u;
        if (wirehair_encode(
                encoder,
                first + offset,
                output.data(),
                context.B,
                &written) != Wirehair_Success ||
            written != context.B)
        {
            success = false;
            break;
        }
        PreliminaryTimingSink ^=
            ConsumePreliminaryTimingOutput(output.data(), written);
    }
    wirehair_free(encoder);
    return success;
}

bool VerifyPreliminaryWH1Systematic(
    const PreliminaryTimingContext& context)
{
    WirehairCodec encoder = wirehair_encoder_create(
        nullptr,
        context.Message.data(),
        context.Message.size(),
        context.B);
    if (!encoder) {
        return false;
    }
    std::vector<uint8_t> output(context.B);
    bool success = true;
    for (uint32_t block_id = 0u; block_id < context.K; ++block_id)
    {
        uint32_t written = 0u;
        if (wirehair_encode(
                encoder,
                block_id,
                output.data(),
                context.B,
                &written) != Wirehair_Success ||
            written != context.B ||
            std::memcmp(
                output.data(),
                context.Message.data() + (size_t)block_id * context.B,
                context.B) != 0)
        {
            success = false;
            break;
        }
    }
    wirehair_free(encoder);
    return success;
}

bool PreparePreliminaryTimingContext(
    uint32_t K,
    uint32_t B,
    PreliminaryTimingContext& context)
{
    if (!IsPreliminaryTimingCell(K, B)) {
        return false;
    }
    PreliminaryTimingContext next;
    next.K = K;
    next.B = B;
    next.PlanEligible = K <= 8u;
    next.Message = MakeMessage(
        (uint64_t)K * B,
        UINT64_C(0x6a09e667f3bcc909) ^ ((uint64_t)K << 32) ^ B);

    MessagePrecodeEncoder selected;
    MessagePrecodeEncoderOptions requested_options;
    if (!PrepareSelected(next.Message, B, requested_options, selected)) {
        return false;
    }
    next.Profile = selected.Profile();
    next.Options = selected.Options();

    const WirehairResult build_result = SystematicGeneratorPlan::Build(
        next.Profile, &next.Options, next.WarmPlan);
    if (next.PlanEligible)
    {
        if (build_result != Wirehair_Success || !next.WarmPlan ||
            next.WarmPlan->Key().SourceBlockCount != K ||
            next.WarmPlan->Key().BlockBytes != B)
        {
            return false;
        }
    }
    else if (build_result != Wirehair_InvalidInput || next.WarmPlan) {
        return false;
    }

    MessagePrecodeEncoder legacy;
    if (legacy.InitializeResult(
            next.Message.data(),
            next.Message.size(),
            B,
            &next.Profile,
            &next.Options) != Wirehair_Success)
    {
        return false;
    }
    MessagePrecodeEncoder cold;
    MessagePrecodeEncoder warm;
    if (!InitializePreliminaryPlanArm(next, true, cold) ||
        !InitializePreliminaryPlanArm(next, false, warm) ||
        !CompareEncoderState(cold, legacy) ||
        !CompareEncoderState(warm, legacy) ||
        !VerifyPreliminaryWH1Systematic(next) ||
        !RunPreliminaryWH1Action(
            next, PreliminaryTimingPhase::RepairK))
    {
        return false;
    }
    context = std::move(next);
    return true;
}

double MedianPreliminaryTiming(std::vector<double>& values)
{
    if (values.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    std::sort(values.begin(), values.end());
    const size_t middle = values.size() / 2u;
    return (values.size() & 1u) != 0u ? values[middle] :
        (values[middle - 1u] + values[middle]) * .5;
}

double TimePreliminaryAction(
    const std::function<bool()>& action,
    uint32_t batch)
{
    const std::chrono::steady_clock::time_point start =
        std::chrono::steady_clock::now();
    for (uint32_t iteration = 0u; iteration < batch; ++iteration) {
        if (!action()) {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
    const std::chrono::steady_clock::time_point end =
        std::chrono::steady_clock::now();
    return (double)std::chrono::duration_cast<std::chrono::nanoseconds>(
        end - start).count() / batch;
}

struct PreliminaryPairTiming
{
    double A = 0.;
    double B = 0.;
    bool Valid = false;
};

PreliminaryPairTiming MeasurePreliminaryPair(
    const std::function<bool()>& action_a,
    const std::function<bool()>& action_b,
    uint32_t replicates,
    uint32_t batch)
{
    PreliminaryPairTiming result;
    for (uint32_t warmup = 0u; warmup < 2u; ++warmup) {
        if (!action_a() || !action_b()) {
            return result;
        }
    }
    std::vector<double> a_values;
    std::vector<double> b_values;
    a_values.reserve(replicates);
    b_values.reserve(replicates);
    for (uint32_t replicate = 0u; replicate < replicates; ++replicate)
    {
        double a_first = 0.;
        double a_second = 0.;
        double b_first = 0.;
        double b_second = 0.;
        if ((replicate & 1u) == 0u)
        {
            a_first = TimePreliminaryAction(action_a, batch);
            b_first = TimePreliminaryAction(action_b, batch);
            b_second = TimePreliminaryAction(action_b, batch);
            a_second = TimePreliminaryAction(action_a, batch);
        }
        else
        {
            b_first = TimePreliminaryAction(action_b, batch);
            a_first = TimePreliminaryAction(action_a, batch);
            a_second = TimePreliminaryAction(action_a, batch);
            b_second = TimePreliminaryAction(action_b, batch);
        }
        if (!std::isfinite(a_first) || !std::isfinite(a_second) ||
            !std::isfinite(b_first) || !std::isfinite(b_second))
        {
            return result;
        }
        a_values.push_back((a_first + a_second) * .5);
        b_values.push_back((b_first + b_second) * .5);
    }
    result.A = MedianPreliminaryTiming(a_values);
    result.B = MedianPreliminaryTiming(b_values);
    result.Valid = std::isfinite(result.A) && result.A > 0. &&
        std::isfinite(result.B) && result.B > 0.;
    return result;
}

bool PrintPreliminaryTimingPhase(
    const PreliminaryTimingContext& context,
    PreliminaryTimingPhase phase)
{
    const std::function<bool()> wh1 = [&]() {
        return RunPreliminaryWH1Action(context, phase);
    };
    const std::function<bool()> legacy = [&]() {
        return RunPreliminaryLegacyAction(context, phase);
    };
    const std::function<bool()> cold = [&]() {
        return RunPreliminaryPlanAction(context, phase, true);
    };
    const std::function<bool()> warm = [&]() {
        return RunPreliminaryPlanAction(context, phase, false);
    };
    const uint32_t replicates = 4u;
    const uint32_t batch = context.K <= 8u ? 2u : 1u;
    const PreliminaryPairTiming aa = MeasurePreliminaryPair(
        legacy, legacy, replicates, batch);
    const PreliminaryPairTiming wh1_pair = MeasurePreliminaryPair(
        wh1, legacy, replicates, batch);
    const PreliminaryPairTiming cold_pair = MeasurePreliminaryPair(
        legacy, cold, replicates, batch);
    const PreliminaryPairTiming warm_pair = MeasurePreliminaryPair(
        legacy, warm, replicates, batch);
    if (!aa.Valid || !wh1_pair.Valid || !cold_pair.Valid ||
        !warm_pair.Valid)
    {
        return false;
    }
    std::printf(
        "PRELIM_TIMING\tK=%u\tB=%u\tphase=%s\tplan_eligible=%u"
        "\twh1_ns=%.3f\tlegacy_wh1_pair_ns=%.3f"
        "\tlegacy_over_wh1=%.9f"
        "\tlegacy_cold_pair_ns=%.3f\tcold_ns=%.3f"
        "\tcold_over_legacy=%.9f"
        "\tlegacy_warm_pair_ns=%.3f\twarm_ns=%.3f"
        "\twarm_over_legacy=%.9f"
        "\taa_a_ns=%.3f\taa_b_ns=%.3f\taa_ratio=%.9f"
        "\tsemantics=PASS\n",
        context.K,
        context.B,
        PreliminaryTimingPhaseName(phase),
        context.PlanEligible ? 1u : 0u,
        wh1_pair.A,
        wh1_pair.B,
        wh1_pair.B / wh1_pair.A,
        cold_pair.A,
        cold_pair.B,
        cold_pair.B / cold_pair.A,
        warm_pair.A,
        warm_pair.B,
        warm_pair.B / warm_pair.A,
        aa.A,
        aa.B,
        aa.B / aa.A);
    return true;
}

bool RunPreliminaryTimingCell(uint32_t K, uint32_t B)
{
    PreliminaryTimingContext context;
    if (!PreparePreliminaryTimingContext(K, B, context)) {
        std::fprintf(
            stderr,
            "PRELIM_FAIL\tK=%u\tB=%u\treason=semantic_setup\n",
            K,
            B);
        return false;
    }
    return PrintPreliminaryTimingPhase(
               context, PreliminaryTimingPhase::FirstK) &&
        PrintPreliminaryTimingPhase(
               context, PreliminaryTimingPhase::RepairK);
}

void PrintPreliminaryTimingProtocol()
{
    std::printf(
        "PRELIM_PROTOCOL\tsame_binary=1\tselected_profile=precomputed"
        "\tpair_order=ABBA_BAAB\treplicates=4"
        "\tlegacy_aa=1\tsemantic_preflight=1"
        "\trole=non_promotional_preliminary\n");
}

bool RunPreliminaryTimingGrid()
{
    for (uint32_t K = 2u; K <= 8u; ++K) {
        for (uint32_t B = 1024u; B <= 2048u; B += 64u) {
            if (!RunPreliminaryTimingCell(K, B)) return false;
        }
    }
    const uint32_t controls[] = { 9u, 32u, 100u };
    for (uint32_t K : controls) {
        for (uint32_t B = 1024u; B <= 2048u; B += 64u) {
            if (!RunPreliminaryTimingCell(K, B)) return false;
        }
    }
    return true;
}

enum class TimingPanelArm : uint8_t
{
    Legacy,
    Cold,
    Warm,
    Wirehair1
};

static const uint32_t kTimingPanelBatch = 8192u;

const char* TimingPanelArmName(TimingPanelArm arm)
{
    switch (arm)
    {
    case TimingPanelArm::Legacy: return "legacy";
    case TimingPanelArm::Cold: return "cold";
    case TimingPanelArm::Warm: return "warm";
    case TimingPanelArm::Wirehair1: return "wh1";
    }
    return "unknown";
}

struct TimingPanelDefinition
{
    const char* Name;
    TimingPanelArm A;
    TimingPanelArm B;
};

static const TimingPanelDefinition TimingPanels[] = {
    { "cold_legacy", TimingPanelArm::Cold, TimingPanelArm::Legacy },
    { "warm_legacy", TimingPanelArm::Warm, TimingPanelArm::Legacy },
    { "warm_wh1", TimingPanelArm::Warm, TimingPanelArm::Wirehair1 },
    { "legacy_legacy", TimingPanelArm::Legacy, TimingPanelArm::Legacy },
    { "cold_cold", TimingPanelArm::Cold, TimingPanelArm::Cold },
    { "warm_warm", TimingPanelArm::Warm, TimingPanelArm::Warm },
    { "wh1_wh1", TimingPanelArm::Wirehair1, TimingPanelArm::Wirehair1 }
};

bool RunTimingPanelArmAction(
    const PreliminaryTimingContext& context,
    PreliminaryTimingPhase phase,
    TimingPanelArm arm)
{
    switch (arm)
    {
    case TimingPanelArm::Legacy:
        return RunPreliminaryLegacyAction(context, phase);
    case TimingPanelArm::Cold:
        return RunPreliminaryPlanAction(context, phase, true);
    case TimingPanelArm::Warm:
        return RunPreliminaryPlanAction(context, phase, false);
    case TimingPanelArm::Wirehair1:
        return RunPreliminaryWH1Action(context, phase);
    }
    return false;
}

bool TimeTimingPanelArm(
    const PreliminaryTimingContext& context,
    PreliminaryTimingPhase phase,
    TimingPanelArm arm,
    uint64_t& elapsed_ns)
{
    const std::chrono::steady_clock::time_point start =
        std::chrono::steady_clock::now();
    for (uint32_t iteration = 0u;
         iteration < kTimingPanelBatch;
         ++iteration)
    {
        if (!RunTimingPanelArmAction(context, phase, arm)) {
            return false;
        }
    }
    const std::chrono::steady_clock::time_point end =
        std::chrono::steady_clock::now();
    elapsed_ns = (uint64_t)
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            end - start).count();
    return elapsed_ns != 0u;
}

struct TimingPanelSlotRecord
{
    uint32_t Lane = 0u;
    uint32_t Round = 0u;
    uint32_t Chronology = 0u;
    uint32_t Cell = 0u;
    uint32_t K = 0u;
    uint32_t B = 0u;
    PreliminaryTimingPhase Phase = PreliminaryTimingPhase::FirstK;
    uint32_t Panel = 0u;
    uint32_t Slot = 0u;
    uint32_t HalfSlot = 0u;
    bool IsABBA = true;
    char Pattern = 'A';
    TimingPanelArm Arm = TimingPanelArm::Legacy;
    uint64_t ElapsedNanoseconds = 0u;
};

char TimingPanelPattern(
    uint32_t slot,
    bool swap_halves,
    bool& is_abba,
    uint32_t& half_slot)
{
    static const char ABBAPattern[4] = {
        'A', 'B', 'B', 'A'
    };
    static const char BAABPattern[4] = {
        'B', 'A', 'A', 'B'
    };
    if (slot >= 8u) {
        return '?';
    }
    is_abba = (slot < 4u) != swap_halves;
    half_slot = (slot & 3u);
    return is_abba ? ABBAPattern[half_slot] : BAABPattern[half_slot];
}

bool RunTimingPanelDurationSmoke()
{
    PreliminaryTimingContext context;
    if (!PreparePreliminaryTimingContext(5u, 1280u, context) ||
        !context.PlanEligible || !context.WarmPlan)
    {
        std::fprintf(stderr, "PANEL_DURATION_SMOKE_FAIL\treason=preflight\n");
        return false;
    }
    uint64_t summed_elapsed_ns = 0u;
    uint32_t slot_count = 0u;
    const std::chrono::steady_clock::time_point wall_start =
        std::chrono::steady_clock::now();
    for (uint32_t phase_index = 0u; phase_index < 2u; ++phase_index)
    {
        const PreliminaryTimingPhase phase = phase_index == 0u ?
            PreliminaryTimingPhase::FirstK : PreliminaryTimingPhase::RepairK;
        for (uint32_t panel_index = 0u;
             panel_index < sizeof(TimingPanels) / sizeof(TimingPanels[0]);
             ++panel_index)
        {
            const TimingPanelDefinition& panel = TimingPanels[panel_index];
            const bool swap_halves = ((phase_index ^ panel_index) & 1u) != 0u;
            for (uint32_t slot = 0u; slot < 8u; ++slot)
            {
                bool is_abba = false;
                uint32_t half_slot = 0u;
                const char pattern = TimingPanelPattern(
                    slot, swap_halves, is_abba, half_slot);
                if (pattern != 'A' && pattern != 'B') {
                    return false;
                }
                const TimingPanelArm arm =
                    pattern == 'A' ? panel.A : panel.B;
                uint64_t elapsed_ns = 0u;
                if (!TimeTimingPanelArm(context, phase, arm, elapsed_ns)) {
                    return false;
                }
                summed_elapsed_ns += elapsed_ns;
                ++slot_count;
                PreliminaryTimingSink ^= (uint64_t)is_abba ^ half_slot;
            }
        }
    }
    const std::chrono::steady_clock::time_point wall_end =
        std::chrono::steady_clock::now();
    if (slot_count != 2u * 7u * 8u) {
        return false;
    }
    const uint64_t wall_ns = (uint64_t)
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            wall_end - wall_start).count();
    const uint64_t projected_timed_ns = summed_elapsed_ns * 9u * 12u;
    const uint64_t projected_wall_ns = wall_ns * 9u * 12u;
    std::printf(
        "PANEL_DURATION_SMOKE\tK=5\tB=1280\tphases=2\tpanels=7"
        "\tslots=%u\tbatch=%u\tsemantic_preflight=PASS"
        "\tsummed_elapsed_ns=%llu\twall_ns=%llu"
        "\tprojection_cells=9\tprojection_rounds=12"
        "\tprojected_timed_wall_ns=%llu\tprojected_total_wall_ns=%llu"
        "\tsink=%llu\tstatus=PASS\n",
        slot_count,
        kTimingPanelBatch,
        (unsigned long long)summed_elapsed_ns,
        (unsigned long long)wall_ns,
        (unsigned long long)projected_timed_ns,
        (unsigned long long)projected_wall_ns,
        (unsigned long long)PreliminaryTimingSink);
    return true;
}

bool RunTimingPanel(uint32_t lane, uint32_t round)
{
    if (lane > 1u || round >= 12u) {
        return false;
    }
    const uint32_t counts[] = { 2u, 5u, 8u };
    const uint32_t widths[] = { 1024u, 1280u, 2048u };
    std::vector<PreliminaryTimingContext> contexts;
    contexts.reserve(9u);
    for (uint32_t K : counts)
    {
        for (uint32_t B : widths)
        {
            PreliminaryTimingContext context;
            if (!PreparePreliminaryTimingContext(K, B, context) ||
                !context.PlanEligible || !context.WarmPlan)
            {
                std::fprintf(
                    stderr,
                    "PANEL_FAIL\tlane=%u\tround=%u\tK=%u\tB=%u"
                    "\treason=semantic_preflight\n",
                    lane,
                    round,
                    K,
                    B);
                return false;
            }
            contexts.push_back(std::move(context));
        }
    }

    const bool reverse = ((lane ^ round) & 1u) != 0u;
    const bool swap_halves =
        ((lane ^ (round >> 1u)) & 1u) != 0u;
    std::printf(
        "PANEL_PROTOCOL\tschema=wirehair.wh2.generator_plan_panel.v2"
        "\tlane=%u\tround=%u\tchronology=%s\tbatch=%u"
        "\tpatterns=ABBA_BAAB\tphysical_half_order=%s\thalf_swap=%u"
        "\tsemantic_preflight=PASS"
        "\trole=internal_gate\n",
        lane,
        round,
        reverse ? "reverse" : "forward",
        kTimingPanelBatch,
        swap_halves ? "BAAB_ABBA" : "ABBA_BAAB",
        swap_halves ? 1u : 0u);
    for (uint32_t cell = 0u; cell < contexts.size(); ++cell)
    {
        const PreliminaryTimingContext& context = contexts[cell];
        std::printf(
            "PANEL_SEMANTIC\tschema=wirehair.wh2.generator_plan_panel.v2"
            "\tlane=%u\tround=%u\tcell=%u\tK=%u\tB=%u"
            "\tstatus=PASS\n",
            lane,
            round,
            cell,
            context.K,
            context.B);
    }
    std::fflush(stdout);

    static const uint32_t kExpectedRecords = 9u * 2u * 7u * 8u;
    std::vector<TimingPanelSlotRecord> records;
    records.reserve(kExpectedRecords);
    uint32_t chronology = 0u;
    for (uint32_t cell_order = 0u; cell_order < contexts.size(); ++cell_order)
    {
        const uint32_t cell = reverse ?
            (uint32_t)contexts.size() - 1u - cell_order : cell_order;
        const PreliminaryTimingContext& context = contexts[cell];
        for (uint32_t phase_order = 0u; phase_order < 2u; ++phase_order)
        {
            const uint32_t phase_index = reverse ? 1u - phase_order :
                phase_order;
            const PreliminaryTimingPhase phase = phase_index == 0u ?
                PreliminaryTimingPhase::FirstK :
                PreliminaryTimingPhase::RepairK;
            for (uint32_t panel_order = 0u;
                 panel_order < sizeof(TimingPanels) / sizeof(TimingPanels[0]);
                 ++panel_order)
            {
                const uint32_t panel_index = reverse ?
                    (uint32_t)(sizeof(TimingPanels) /
                        sizeof(TimingPanels[0])) - 1u - panel_order :
                    panel_order;
                const TimingPanelDefinition& panel =
                    TimingPanels[panel_index];
                for (uint32_t slot = 0u; slot < 8u; ++slot)
                {
                    bool is_abba = false;
                    uint32_t half_slot = 0u;
                    const char pattern = TimingPanelPattern(
                        slot, swap_halves, is_abba, half_slot);
                    const TimingPanelArm arm =
                        pattern == 'A' ? panel.A : panel.B;
                    uint64_t elapsed_ns = 0u;
                    if (!TimeTimingPanelArm(
                            context, phase, arm, elapsed_ns))
                    {
                        std::fprintf(
                            stderr,
                            "PANEL_FAIL\tlane=%u\tround=%u\tK=%u\tB=%u"
                            "\tphase=%s\tpanel=%s\tslot=%u"
                            "\treason=timed_action\n",
                            lane,
                            round,
                            context.K,
                            context.B,
                            PreliminaryTimingPhaseName(phase),
                            panel.Name,
                            slot + 1u);
                        return false;
                    }
                    TimingPanelSlotRecord record;
                    record.Lane = lane;
                    record.Round = round;
                    record.Chronology = chronology++;
                    record.Cell = cell;
                    record.K = context.K;
                    record.B = context.B;
                    record.Phase = phase;
                    record.Panel = panel_index;
                    record.Slot = slot + 1u;
                    record.HalfSlot = half_slot + 1u;
                    record.IsABBA = is_abba;
                    record.Pattern = pattern;
                    record.Arm = arm;
                    record.ElapsedNanoseconds = elapsed_ns;
                    records.push_back(record);
                }
            }
        }
    }
    if (records.size() != kExpectedRecords) {
        return false;
    }
    for (const TimingPanelSlotRecord& record : records)
    {
        const TimingPanelDefinition& panel = TimingPanels[record.Panel];
        std::printf(
            "PANEL_SLOT\tschema=wirehair.wh2.generator_plan_panel.v2"
            "\tlane=%u\tround=%u\tchronology=%u\tcell=%u"
            "\tK=%u\tB=%u\tphase=%s\tpanel=%s\tslot=%u"
            "\thalf=%s\thalf_slot=%u\tpattern=%c\tarm=%s"
            "\tbatch=%u\telapsed_ns=%llu"
            "\tns_per_action=%.9f\n",
            record.Lane,
            record.Round,
            record.Chronology,
            record.Cell,
            record.K,
            record.B,
            PreliminaryTimingPhaseName(record.Phase),
            panel.Name,
            record.Slot,
            record.IsABBA ? "ABBA" : "BAAB",
            record.HalfSlot,
            record.Pattern,
            TimingPanelArmName(record.Arm),
            kTimingPanelBatch,
            (unsigned long long)record.ElapsedNanoseconds,
            (double)record.ElapsedNanoseconds / kTimingPanelBatch);
    }
    std::printf(
        "PANEL_DONE\tschema=wirehair.wh2.generator_plan_panel.v2"
        "\tlane=%u\tround=%u\trecords=%u\tsink=%llu\tstatus=PASS\n",
        lane,
        round,
        (unsigned)records.size(),
        (unsigned long long)PreliminaryTimingSink);
    return true;
}

bool ParsePreliminaryU32(const char* text, uint32_t& value)
{
    if (!text || !*text || *text == '-' || *text == '+') {
        return false;
    }
    errno = 0;
    char* end = nullptr;
    const unsigned long parsed = std::strtoul(text, &end, 0);
    if (errno != 0 || !end || end == text || *end != '\0' ||
        parsed > UINT32_MAX)
    {
        return false;
    }
    value = (uint32_t)parsed;
    return true;
}

} // namespace

int main(int argc, char** argv)
{
    if (wirehair_init() != Wirehair_Success) {
        std::fprintf(stderr, "generator plan: wirehair_init failed\n");
        return 1;
    }
    gf256_x86_cpu_features features = {};
    gf256_get_active_x86_cpu_features(&features);
    if (!features.GFNI) {
        if (!RunUnsupportedGFNIFallback()) {
            std::fprintf(
                stderr,
                "systematic generator plan: non-GFNI fallback failed\n");
            return 2;
        }
        if (argc > 1) {
            std::fprintf(
                stderr,
                "generator plan timing requires the active GFNI backend\n");
            return 2;
        }
        std::printf(
            "systematic generator plan: SKIP (active GFNI unavailable)\n");
        return 0;
    }
    if (argc == 4 && std::strcmp(argv[1], "--timing-cell") == 0)
    {
        uint32_t K = 0u;
        uint32_t B = 0u;
        if (!ParsePreliminaryU32(argv[2], K) ||
            !ParsePreliminaryU32(argv[3], B) ||
            !IsPreliminaryTimingCell(K, B))
        {
            std::fprintf(stderr, "invalid preliminary timing cell\n");
            return 2;
        }
        PrintPreliminaryTimingProtocol();
        return RunPreliminaryTimingCell(K, B) ? 0 : 2;
    }
    if (argc == 4 && std::strcmp(argv[1], "--timing-panel") == 0)
    {
        uint32_t lane = 0u;
        uint32_t round = 0u;
        if (!ParsePreliminaryU32(argv[2], lane) ||
            !ParsePreliminaryU32(argv[3], round) ||
            lane > 1u || round >= 12u)
        {
            std::fprintf(stderr, "invalid timing panel lane/round\n");
            return 2;
        }
        return RunTimingPanel(lane, round) ? 0 : 2;
    }
    if (argc == 2 &&
        std::strcmp(argv[1], "--timing-panel-duration-smoke") == 0)
    {
        return RunTimingPanelDurationSmoke() ? 0 : 2;
    }
    if (argc == 2 &&
        std::strcmp(argv[1], "--timing-preliminary") == 0)
    {
        PrintPreliminaryTimingProtocol();
        return RunPreliminaryTimingGrid() ? 0 : 2;
    }
    if (argc != 1)
    {
        std::fprintf(
            stderr,
            "usage: %s [--timing-cell K B | --timing-panel lane round | "
            "--timing-panel-duration-smoke | --timing-preliminary]\n",
            argv[0]);
        return 2;
    }
    if (!RunSemanticGrid()) {
        return 2;
    }
    if (!RunMismatchAndFallback()) {
        std::fprintf(stderr, "generator plan: mismatch/fallback failed\n");
        return 3;
    }
    if (!RunEligibilityAndBuildOOM()) {
        std::fprintf(stderr, "generator plan: eligibility/build OOM failed\n");
        return 4;
    }
    if (!RunReplayOOMAndTransaction()) {
        std::fprintf(stderr, "generator plan: replay OOM/transaction failed\n");
        return 5;
    }
    if (!RunInputAlias()) {
        std::fprintf(stderr, "generator plan: input alias failed\n");
        return 6;
    }
    if (!RunConcurrentReuse()) {
        std::fprintf(stderr, "generator plan: concurrent reuse failed\n");
        return 7;
    }
    if (!RunRecovery()) {
        std::fprintf(stderr, "generator plan: recovery failed\n");
        return 8;
    }
    std::printf("systematic generator plan: PASS\n");
    return 0;
}
