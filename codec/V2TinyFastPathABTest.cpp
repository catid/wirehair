/*
    Tiny-K mixed fast-path A/B equality oracle.

    The tiny fast path claims to be equation-preserving: forcing it on or
    off must never change a solve result code, an intermediate block byte,
    or (therefore) an encoder output byte.  This test proves that directly:

      * scalar helper identities (exhaustive GF(2^16) inverse, fused row
        kernels versus the scalar reference kernels);
      * full encode A/B: for a grid of (profile, K, construction seed,
        payload size), byte-compare the complete intermediate block vector
        and K generated repair symbols with the fast path forced on versus
        forced off;
      * decode A/B: deterministic systematic loss patterns repaired from the
        encoded repair stream, including rank-deficient (NeedMore) and
        corrupted (Error) classifications, must agree exactly.

    Default mode is a bounded grid suitable for ctest; --full runs the
    dense verification ladder (K stepped densely through 1024, 8
    construction seeds).
*/

#include "WirehairV2Precode.h"
#include "WirehairV2Solve.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using namespace wirehair_v2;

namespace {

uint64_t SplitMix64(uint64_t& state)
{
    state += UINT64_C(0x9E3779B97F4A7C15);
    uint64_t z = state;
    z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
    return z ^ (z >> 31);
}

bool CheckTinyMixedFastPathModeControl()
{
    if (!SetTinyMixedFastPathModeForTesting(0) ||
        TinyMixedFastPathModeForTesting() != 0)
    {
        return false;
    }
    static const int kModes[] = {-1, 1, 0};
    for (int mode : kModes)
    {
        if (!SetTinyMixedFastPathModeForTesting(mode) ||
            TinyMixedFastPathModeForTesting() != mode)
        {
            return false;
        }
    }
    if (SetTinyMixedFastPathModeForTesting(-2) ||
        TinyMixedFastPathModeForTesting() != 0 ||
        SetTinyMixedFastPathModeForTesting(2) ||
        TinyMixedFastPathModeForTesting() != 0)
    {
        return false;
    }
    return true;
}

struct ProfileConfig
{
    const char* Name;
    bool Mixed;
    // ic-d8-9x2-p48x style knobs (zero/false = profile defaults).
    bool SharedX;
    uint32_t Period;
    uint32_t GF256Rows;
    uint32_t DenseRows;
    bool IdentityCorner;
    uint32_t MinK;
    bool CanonicalDispatch;
    MixedResidueBucketMode BucketMode;
};

bool ApplyProfileHooks(const ProfileConfig& profile)
{
    if (!SetMixedCoefficientGeometryForTesting(
            profile.SharedX ?
                MixedCoefficientGeometry::SharedCauchyX :
                MixedCoefficientGeometry::FrozenPowerX))
    {
        return false;
    }
    if (!SetMixedGF16RowsForTesting(2u)) {
        return false;
    }
    if (!SetMixedCoefficientPeriodForTesting(
            profile.Period != 0u ? profile.Period :
                kMixedCoefficientPeriod))
    {
        return false;
    }
    if (!SetMixedGF256RowsForTesting(
            profile.GF256Rows != 0u ? profile.GF256Rows :
                kMixedGF256Rows))
    {
        return false;
    }
    return SetMixedResidueBucketModeForTesting(profile.BucketMode);
}

struct SolveOutcome
{
    WirehairResult Result = Wirehair_Error;
    std::vector<uint8_t> Blocks;
    PrecodeSolveStats Stats;
};

SolveOutcome RunSolve(
    int fast_path_mode,
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const std::vector<SolvePacket>& packets,
    uint32_t block_bytes)
{
    SolveOutcome outcome;
    if (!SetTinyMixedFastPathModeForTesting(fast_path_mode) ||
        TinyMixedFastPathModeForTesting() != fast_path_mode)
    {
        outcome.Result = Wirehair_Error;
        return outcome;
    }
    outcome.Result = SolvePrecodeSystem(
        system, config, packets, block_bytes, outcome.Blocks,
        &outcome.Stats);
    if (!SetTinyMixedFastPathModeForTesting(0) ||
        TinyMixedFastPathModeForTesting() != 0)
    {
        outcome.Result = Wirehair_Error;
        outcome.Blocks.clear();
        outcome.Stats = PrecodeSolveStats();
    }
    return outcome;
}

struct Counters
{
    uint64_t EncodeCases = 0;
    uint64_t DecodeCases = 0;
    uint64_t SuccessCases = 0;
    uint64_t NeedMoreCases = 0;
    uint64_t ErrorCases = 0;
    uint64_t SuccessfulForcedOnMixedSolves = 0;
    uint64_t ForcedOnMixedAcceptances = 0;
    uint64_t ForcedOffAcceptances = 0;
    uint64_t ForcedOnNonMixedAcceptances = 0;
    uint64_t SuccessfulAutomaticEligibleMixedSolves = 0;
    uint64_t AutomaticMixedAcceptances = 0;
    uint64_t AutomaticIneligibleAcceptances = 0;
    uint64_t CanonicalDispatchEncodeCases = 0;
    uint32_t AutomaticExclusionSeamCases = 0;
    bool CanonicalDispatchWidth2Success[101] = {};
    bool CanonicalDispatchWidth4Success[101] = {};
};

bool SamePathIndependentStats(
    const PrecodeSolveStats& a,
    const PrecodeSolveStats& b)
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
        a.PeelAdjacencyVisits == b.PeelAdjacencyVisits &&
        a.PeelRowScanSteps == b.PeelRowScanSteps &&
        a.PeelHeapOperations == b.PeelHeapOperations &&
        a.ProjectionWordXors == b.ProjectionWordXors &&
        a.ResidualCoeffWordXors == b.ResidualCoeffWordXors &&
        a.ResidualCoeffByteOps == b.ResidualCoeffByteOps &&
        a.PacketSeedAttempt == b.PacketSeedAttempt;
}

bool CompareOutcomes(
    const char* what,
    const ProfileConfig& profile,
    uint32_t K,
    uint64_t seed,
    uint32_t block_bytes,
    const SolveOutcome& off,
    const SolveOutcome& on,
    const SolveOutcome& automatic,
    Counters& counters)
{
    const uint32_t off_acceptances =
        off.Stats.TinyMixedFastPathAcceptances;
    const uint32_t on_acceptances =
        on.Stats.TinyMixedFastPathAcceptances;
    counters.ForcedOffAcceptances += off_acceptances;
    if (off_acceptances != 0u)
    {
        std::fprintf(stderr,
            "FAIL %s %s K=%u seed=0x%llx bb=%u forced-off accepted "
            "tiny path %u times\n",
            what, profile.Name, K, (unsigned long long)seed, block_bytes,
            off_acceptances);
        return false;
    }
    if (!profile.Mixed)
    {
        counters.ForcedOnNonMixedAcceptances += on_acceptances;
        if (on_acceptances != 0u)
        {
            std::fprintf(stderr,
                "FAIL %s %s K=%u seed=0x%llx bb=%u nonmixed solve "
                "accepted tiny path %u times\n",
                what, profile.Name, K, (unsigned long long)seed, block_bytes,
                on_acceptances);
            return false;
        }
    }
    else
    {
        counters.ForcedOnMixedAcceptances += on_acceptances;
        if (on_acceptances > 1u)
        {
            std::fprintf(stderr,
                "FAIL %s %s K=%u seed=0x%llx bb=%u forced-on mixed "
                "solve accepted tiny path %u times\n",
                what, profile.Name, K, (unsigned long long)seed, block_bytes,
                on_acceptances);
            return false;
        }
        if (on.Result == Wirehair_Success)
        {
            ++counters.SuccessfulForcedOnMixedSolves;
            if (on_acceptances != 1u)
            {
                std::fprintf(stderr,
                    "FAIL %s %s K=%u seed=0x%llx bb=%u successful "
                    "forced-on mixed solve did not accept tiny path\n",
                    what, profile.Name, K, (unsigned long long)seed,
                    block_bytes);
                return false;
            }
        }
    }
    const uint32_t automatic_acceptances =
        automatic.Stats.TinyMixedFastPathAcceptances;
    if (profile.Mixed) {
        counters.AutomaticMixedAcceptances += automatic_acceptances;
    }
    if (automatic_acceptances > 1u)
    {
        std::fprintf(stderr,
            "FAIL %s %s K=%u seed=0x%llx bb=%u automatic solve "
            "accepted tiny path %u times\n",
            what, profile.Name, K, (unsigned long long)seed, block_bytes,
            automatic_acceptances);
        return false;
    }
    const bool automatic_eligible =
        profile.Mixed &&
        K <= TinyMixedFastPathMaxSourceBlocks() &&
        block_bytes <= TinyMixedFastPathMaxBlockBytes() &&
        (uint64_t)K * block_bytes <=
            TinyMixedFastPathMaxProductBytes() &&
        profile.BucketMode == MixedResidueBucketMode::Automatic;
    if (!automatic_eligible)
    {
        counters.AutomaticIneligibleAcceptances += automatic_acceptances;
        if (automatic_acceptances != 0u)
        {
            std::fprintf(stderr,
                "FAIL %s %s K=%u seed=0x%llx bb=%u ineligible "
                "automatic solve accepted tiny path %u times\n",
                what, profile.Name, K, (unsigned long long)seed, block_bytes,
                automatic_acceptances);
            return false;
        }
    }
    else
    {
        if (automatic_acceptances != on_acceptances)
        {
            std::fprintf(stderr,
                "FAIL %s %s K=%u seed=0x%llx bb=%u eligible automatic "
                "acceptance %u differs from forced-on %u\n",
                what, profile.Name, K, (unsigned long long)seed,
                block_bytes, automatic_acceptances, on_acceptances);
            return false;
        }
        if (automatic.Result == Wirehair_Success)
        {
            ++counters.SuccessfulAutomaticEligibleMixedSolves;
            if (automatic_acceptances != 1u)
            {
                std::fprintf(stderr,
                    "FAIL %s %s K=%u seed=0x%llx bb=%u eligible "
                    "successful automatic solve did not accept tiny path\n",
                    what, profile.Name, K, (unsigned long long)seed,
                    block_bytes);
                return false;
            }
        }
    }
    if (!SamePathIndependentStats(off.Stats, on.Stats) ||
        !SamePathIndependentStats(off.Stats, automatic.Stats))
    {
        std::fprintf(stderr,
            "FAIL %s %s K=%u seed=0x%llx bb=%u path-independent "
            "stats differ: off rank=%u/%u on=%u/%u auto=%u/%u\n",
            what, profile.Name, K, (unsigned long long)seed, block_bytes,
            off.Stats.BinaryResidualRank, off.Stats.ResidualRank,
            on.Stats.BinaryResidualRank, on.Stats.ResidualRank,
            automatic.Stats.BinaryResidualRank,
            automatic.Stats.ResidualRank);
        return false;
    }
    if (off.Result != on.Result || off.Result != automatic.Result)
    {
        std::fprintf(stderr,
            "FAIL %s %s K=%u seed=0x%llx bb=%u result "
            "off=%d on=%d auto=%d\n",
            what, profile.Name, K, (unsigned long long)seed, block_bytes,
            (int)off.Result, (int)on.Result, (int)automatic.Result);
        return false;
    }
    if (off.Result == Wirehair_Success)
    {
        ++counters.SuccessCases;
        if (off.Blocks != on.Blocks || off.Blocks != automatic.Blocks)
        {
            std::fprintf(stderr,
                "FAIL %s %s K=%u seed=0x%llx bb=%u intermediate "
                "blocks differ across off/on/automatic\n",
                what, profile.Name, K, (unsigned long long)seed,
                block_bytes);
            return false;
        }
    }
    else if (off.Result == Wirehair_NeedMore) {
        ++counters.NeedMoreCases;
    }
    else {
        ++counters.ErrorCases;
    }
    return true;
}

bool RunCase(
    const ProfileConfig& profile,
    uint32_t K,
    uint64_t construction_seed,
    uint32_t block_bytes,
    bool decode_variants,
    Counters& counters)
{
    if (profile.Mixed && !ApplyProfileHooks(profile)) {
        std::fprintf(stderr, "FAIL %s hooks\n", profile.Name);
        return false;
    }
    PrecodeParams params = profile.Mixed ?
        MakeMixedParams(K, construction_seed) :
        MakeCertifiedParams(K, construction_seed);
    if (profile.DenseRows != 0u) {
        params.DenseRows = profile.DenseRows;
    }
    params.DenseIdentityCorner = profile.IdentityCorner;
    if (profile.CanonicalDispatch &&
        (!profile.Mixed || profile.SharedX ||
         profile.Period != kMixedCoefficientPeriod ||
         profile.GF256Rows != kMixedGF256Rows ||
         profile.DenseRows != 4u || profile.IdentityCorner ||
         profile.MinK != 2u ||
         !IsCanonicalMixedCompletionState() ||
         !IsCanonicalStableTargetPacketRowState() ||
         ActiveMixedResidueBucketModeForTesting() !=
            MixedResidueBucketMode::Automatic ||
         params.BlockCount != K ||
         params.Staircase != SmallBandStaircaseCount(K) ||
         params.DenseRows != 4u ||
         params.HeavyRows != kMixedGF256Rows + kMixedGF16Rows ||
         params.Field != CompletionField::MixedGF256GF16 ||
         params.DenseIdentityCorner || params.DenseTwoAnchor ||
         params.DenseTwoAnchorPhase != 0u ||
         params.DegreeBalancedStaircase ||
         params.StaircaseDegreeScale != kStaircaseDegreeScaleUnset ||
         params.Seed != construction_seed))
    {
        std::fprintf(stderr,
            "FAIL %s K=%u seed=0x%llx canonical dispatch params drifted\n",
            profile.Name, K, (unsigned long long)construction_seed);
        return false;
    }
    PrecodeSystem system;
    if (!BuildPrecodeSystem(params, system))
    {
        std::fprintf(stderr,
            "FAIL %s K=%u seed=0x%llx construction rejected\n",
            profile.Name, K, (unsigned long long)construction_seed);
        return false;
    }
    PacketRowConfig config;
    config.PeelSeed = (uint32_t)construction_seed ^
        (uint32_t)(construction_seed >> 32);
    config.MixCount = kCertifiedPacketMixCount;
    if (profile.CanonicalDispatch &&
        (config.PeelSeed != RawUniformPacketPeelSeed(construction_seed) ||
         config.MixCount != kCertifiedPacketMixCount))
    {
        std::fprintf(stderr,
            "FAIL %s K=%u seed=0x%llx canonical packet config drifted\n",
            profile.Name, K, (unsigned long long)construction_seed);
        return false;
    }

    // Deterministic payload.
    uint64_t payload_state = construction_seed ^
        ((uint64_t)K << 32) ^ block_bytes;
    std::vector<uint8_t> source((size_t)K * block_bytes);
    for (uint8_t& byte : source) {
        byte = (uint8_t)SplitMix64(payload_state);
    }
    std::vector<SolvePacket> packets(K);
    for (uint32_t id = 0; id < K; ++id)
    {
        packets[id].BlockId = id;
        packets[id].Data = source.data() + (size_t)id * block_bytes;
    }

    ++counters.EncodeCases;
    if (profile.CanonicalDispatch) {
        ++counters.CanonicalDispatchEncodeCases;
    }
    const SolveOutcome off =
        RunSolve(-1, system, config, packets, block_bytes);
    const SolveOutcome on =
        RunSolve(1, system, config, packets, block_bytes);
    const SolveOutcome automatic =
        RunSolve(0, system, config, packets, block_bytes);
    if (!CompareOutcomes(
            "encode", profile, K, construction_seed, block_bytes,
            off, on, automatic, counters))
    {
        return false;
    }
    if (profile.CanonicalDispatch && K <= 100u &&
        on.Result == Wirehair_Success &&
        on.Stats.TinyMixedFastPathAcceptances == 1u)
    {
        if (block_bytes == 2u) {
            counters.CanonicalDispatchWidth2Success[K] = true;
        }
        else if (block_bytes == 4u) {
            counters.CanonicalDispatchWidth4Success[K] = true;
        }
    }

    std::vector<uint8_t> repair_off;
    std::vector<uint8_t> repair_on;
    if (off.Result == Wirehair_Success)
    {
        // Full encoder output: systematic symbols are the source data by
        // construction; generate K repair symbols from each intermediate
        // vector and byte-compare.
        repair_off.resize((size_t)K * block_bytes);
        repair_on.resize((size_t)K * block_bytes);
        for (uint32_t repair = 0; repair < K; ++repair)
        {
            const uint32_t block_id = K + repair;
            if (!EvaluatePacketBlock(
                    system, config, off.Blocks.data(), block_bytes,
                    block_id,
                    repair_off.data() + (size_t)repair * block_bytes) ||
                !EvaluatePacketBlock(
                    system, config, on.Blocks.data(), block_bytes,
                    block_id,
                    repair_on.data() + (size_t)repair * block_bytes))
            {
                std::fprintf(stderr,
                    "FAIL %s K=%u seed=0x%llx repair evaluation\n",
                    profile.Name, K,
                    (unsigned long long)construction_seed);
                return false;
            }
        }
        if (repair_off != repair_on)
        {
            std::fprintf(stderr,
                "FAIL %s K=%u seed=0x%llx repair symbols differ\n",
                profile.Name, K, (unsigned long long)construction_seed);
            return false;
        }
    }

    if (!decode_variants || off.Result != Wirehair_Success) {
        return true;
    }

    // Decode A/B: drop ~10% of systematic symbols deterministically and
    // replace them with repair symbols; also probe a rank-deficient
    // partial set (NeedMore) and a corrupted duplicate (Error).
    uint64_t loss_state = construction_seed ^ 0xD1CEu ^ K;
    std::vector<uint8_t> kept(K, 1u);
    uint32_t dropped = 0u;
    for (uint32_t id = 0; id < K; ++id)
    {
        if ((SplitMix64(loss_state) % 10u) == 0u) {
            kept[id] = 0u;
            ++dropped;
        }
    }
    if (dropped == 0u && K > 2u)
    {
        kept[SplitMix64(loss_state) % K] = 0u;
        dropped = 1u;
    }
    const uint32_t extra = 2u;
    std::vector<SolvePacket> received;
    received.reserve((size_t)K + extra);
    for (uint32_t id = 0; id < K; ++id)
    {
        if (kept[id]) {
            received.push_back(packets[id]);
        }
    }
    for (uint32_t repair = 0;
         repair < dropped + extra && repair < K; ++repair)
    {
        SolvePacket packet;
        packet.BlockId = K + repair;
        packet.Data = repair_off.data() + (size_t)repair * block_bytes;
        received.push_back(packet);
    }
    ++counters.DecodeCases;
    const SolveOutcome decode_off =
        RunSolve(-1, system, config, received, block_bytes);
    const SolveOutcome decode_on =
        RunSolve(1, system, config, received, block_bytes);
    const SolveOutcome decode_automatic =
        RunSolve(0, system, config, received, block_bytes);
    if (!CompareOutcomes(
            "decode", profile, K, construction_seed, block_bytes,
            decode_off, decode_on, decode_automatic, counters))
    {
        return false;
    }
    if (decode_off.Result == Wirehair_Success &&
        decode_off.Blocks != off.Blocks)
    {
        std::fprintf(stderr,
            "FAIL %s K=%u seed=0x%llx decode disagrees with encode\n",
            profile.Name, K, (unsigned long long)construction_seed);
        return false;
    }

    // Rank-deficient partial set: strictly fewer than K packets after the
    // mandatory count check requires more symbols in both modes.
    if (received.size() > K)
    {
        std::vector<SolvePacket> deficient(
            received.begin(), received.begin() + K);
        // Duplicate one packet so the K received rows cannot have full rank.
        deficient[0] = deficient[K > 1u ? 1u : 0u];
        ++counters.DecodeCases;
        const SolveOutcome deficient_off =
            RunSolve(-1, system, config, deficient, block_bytes);
        const SolveOutcome deficient_on =
            RunSolve(1, system, config, deficient, block_bytes);
        const SolveOutcome deficient_automatic =
            RunSolve(0, system, config, deficient, block_bytes);
        if (!CompareOutcomes(
                "deficient", profile, K, construction_seed, block_bytes,
                deficient_off, deficient_on, deficient_automatic, counters))
        {
            return false;
        }
    }

    // Corrupted duplicate: a repeated block id with different payload makes
    // the system inconsistent; classification must agree byte-for-byte.
    {
        std::vector<uint8_t> corrupted(block_bytes);
        std::memcpy(
            corrupted.data(), received.back().Data, block_bytes);
        corrupted[0] ^= 0x5Au;
        std::vector<SolvePacket> conflicted = received;
        SolvePacket duplicate = received.back();
        duplicate.Data = corrupted.data();
        conflicted.push_back(duplicate);
        ++counters.DecodeCases;
        const SolveOutcome conflict_off =
            RunSolve(-1, system, config, conflicted, block_bytes);
        const SolveOutcome conflict_on =
            RunSolve(1, system, config, conflicted, block_bytes);
        const SolveOutcome conflict_automatic =
            RunSolve(0, system, config, conflicted, block_bytes);
        if (!CompareOutcomes(
                "conflict", profile, K, construction_seed, block_bytes,
                conflict_off, conflict_on, conflict_automatic, counters))
        {
            return false;
        }
    }
    return true;
}

std::vector<uint32_t> BuildKList(bool full)
{
    std::vector<uint32_t> k_values;
    if (full)
    {
        for (uint32_t k = 2; k <= 128u; ++k) {
            k_values.push_back(k);
        }
        for (uint32_t k = 130u; k <= 256u; k += 2u) {
            k_values.push_back(k);
        }
        for (uint32_t k = 260u; k <= 512u; k += 4u) {
            k_values.push_back(k);
        }
        // Fine steps across the product-gate seam (K*bb = 3072 excludes
        // K=513 at 6-byte blocks): odd K right at the first excluded cells.
        for (uint32_t k = 513u; k <= 519u; k += 2u) {
            k_values.push_back(k);
        }
        for (uint32_t k = 520u; k <= 1024u; k += 8u) {
            k_values.push_back(k);
        }
    }
    else
    {
        for (uint32_t k = 2; k <= 40u; ++k) {
            k_values.push_back(k);
        }
        static const uint32_t kSparse[] = {
            48u, 56u, 64u, 80u, 96u, 128u, 192u, 256u, 384u, 512u,
            768u, 1024u};
        for (uint32_t k : kSparse) {
            k_values.push_back(k);
        }
    }
    return k_values;
}

} // namespace

namespace {

// Nonzero-payload solve timing: unlike the bench's zero-payload timing
// loop, every peeled constant is nonzero here, so the fast path's
// zero-block skip cannot flatter it.  Prints one median per (K, bb, mode).
// The optional profile selector times either the frozen baseline or the
// dispatched mixed-p244-d4 configuration (dense rows 4).
int RunTiming(
    uint32_t K,
    uint32_t block_bytes,
    int mode,
    uint32_t reps,
    const char* profile_name)
{
    ProfileConfig profile =
        {"mixed-p244-baseline", true, false, 0u, 0u, 0u, false, 2u,
            false, MixedResidueBucketMode::Automatic};
    if (profile_name && !std::strcmp(profile_name, "d4")) {
        profile.Name = "mixed-p244-d4";
        profile.Period = kMixedCoefficientPeriod;
        profile.GF256Rows = kMixedGF256Rows;
        profile.DenseRows = 4u;
        profile.CanonicalDispatch = true;
    }
    else if (profile_name && std::strcmp(profile_name, "baseline")) {
        std::fprintf(stderr,
            "unknown timing profile %s (want baseline or d4)\n",
            profile_name);
        return 1;
    }
    if (!ApplyProfileHooks(profile)) {
        return 1;
    }
    PrecodeParams params = MakeMixedParams(K, UINT64_C(0xC0FFEE));
    if (profile.DenseRows != 0u) {
        params.DenseRows = profile.DenseRows;
    }
    PrecodeSystem system;
    if (!BuildPrecodeSystem(params, system)) {
        return 1;
    }
    PacketRowConfig config;
    config.PeelSeed = 0xC0FFEEu;
    config.MixCount = kCertifiedPacketMixCount;
    uint64_t payload_state = K ^ ((uint64_t)block_bytes << 32);
    std::vector<uint8_t> source((size_t)K * block_bytes);
    for (uint8_t& byte : source) {
        byte = (uint8_t)(SplitMix64(payload_state) | 1u);
    }
    std::vector<SolvePacket> packets(K);
    for (uint32_t id = 0; id < K; ++id)
    {
        packets[id].BlockId = id;
        packets[id].Data = source.data() + (size_t)id * block_bytes;
    }
    if (!SetTinyMixedFastPathModeForTesting(mode)) {
        return 1;
    }
    std::vector<uint64_t> rep_ns;
    rep_ns.reserve(reps);
    for (uint32_t rep = 0; rep < reps; ++rep)
    {
        std::vector<uint8_t> blocks;
        const std::chrono::steady_clock::time_point t0 =
            std::chrono::steady_clock::now();
        const WirehairResult result = SolvePrecodeSystem(
            system, config, packets, block_bytes, blocks);
        const std::chrono::steady_clock::time_point t1 =
            std::chrono::steady_clock::now();
        if (result != Wirehair_Success && result != Wirehair_NeedMore) {
            return 1;
        }
        rep_ns.push_back((uint64_t)
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                t1 - t0).count());
    }
    std::sort(rep_ns.begin(), rep_ns.end());
    std::printf("timing profile=%s K=%u bb=%u mode=%d median_solve_ns=%llu\n",
        profile.Name, K, block_bytes, mode,
        (unsigned long long)rep_ns[rep_ns.size() / 2u]);
    return 0;
}

} // namespace

int main(int argc, char** argv)
{
    bool full = false;
    for (int i = 1; i < argc; ++i)
    {
        if (!std::strcmp(argv[i], "--full")) {
            full = true;
        }
        else if (!std::strcmp(argv[i], "--timing") && i + 4 < argc)
        {
            const uint32_t K = (uint32_t)std::strtoul(
                argv[i + 1], nullptr, 0);
            const uint32_t block_bytes = (uint32_t)std::strtoul(
                argv[i + 2], nullptr, 0);
            const int mode = (int)std::strtol(argv[i + 3], nullptr, 0);
            const uint32_t reps = (uint32_t)std::strtoul(
                argv[i + 4], nullptr, 0);
            const char* profile_name = i + 5 < argc ? argv[i + 5] : nullptr;
            return RunTiming(K, block_bytes, mode, reps, profile_name);
        }
        else
        {
            std::fprintf(stderr, "unknown argument: %s\n", argv[i]);
            return 1;
        }
    }

    if (!CheckTinyMixedScalarHelpersForTesting())
    {
        std::fprintf(stderr,
            "FAIL tiny scalar helper verification\n");
        return 1;
    }
    if (!CheckTinyMixedFastPathModeControl())
    {
        std::fprintf(stderr,
            "FAIL tiny fast-path mode control verification\n");
        return 1;
    }
    std::printf("tiny scalar helpers: exhaustive inverse and row kernels "
                "verified; mode control verified\n");

    static const ProfileConfig kProfiles[] = {
        // mixed-p244-baseline: frozen geometry defaults
        {"mixed-p244-baseline", true, false, 0u, 0u, 0u, false, 2u,
            false, MixedResidueBucketMode::Automatic},
        // Exact dispatch-v1 completion band: frozen P244, 10 GF(256) plus
        // two GF(2^16) rows, and four binary dense rows.
        {"mixed-p244-d4", true, false, kMixedCoefficientPeriod,
            kMixedGF256Rows, 4u, false, 2u, true,
            MixedResidueBucketMode::Automatic},
        // ic-d8-9x2-p48x: shared-x, period 48, 9 GF256 rows, D2=8,
        // identity corner.  The small-band staircase gives S=2 at K=5, so
        // K+S=7 cannot carry the eight distinct corner anchors; K=6 is the
        // first valid construction.
        {"ic-d8-9x2-p48x", true, true, 48u, 9u, 8u, true, 6u,
            false, MixedResidueBucketMode::Automatic},
        // certified control: never touches the mixed fast path
        {"cert-baseline", false, false, 0u, 0u, 0u, false, 2u,
            false, MixedResidueBucketMode::Automatic},
    };

    const std::vector<uint32_t> k_values = BuildKList(full);
    const uint32_t seed_count = full ? 8u : 3u;
    Counters counters;
    for (const ProfileConfig& profile : kProfiles)
    {
        for (uint32_t k_index = 0; k_index < k_values.size(); ++k_index)
        {
            const uint32_t K = k_values[k_index];
            if (K < profile.MinK) {
                continue;
            }
            for (uint32_t seed_index = 0;
                 seed_index < seed_count; ++seed_index)
            {
                const uint64_t construction_seed =
                    UINT64_C(0xC0FFEE) +
                    seed_index * UINT64_C(0x9E3779B97F4A7C15);
                // Payload widths: harness width, small even widths crossing
                // the scalar/kernel regimes, and an MTU-scale width kept to
                // bounded K so the forced scalar path stays fast enough.
                std::vector<uint32_t> widths;
                widths.push_back(2u);
                if (profile.Mixed) {
                    widths.push_back(4u);
                }
                widths.push_back(6u);
                // 8-byte blocks reach the widened product gate (K*bb <=
                // 3072 engages bb=8 through K=384); cover the flat-K
                // region plus the first excluded cells above it.
                if (K <= 512u) {
                    widths.push_back(8u);
                }
                if ((k_index + seed_index) % 2u == 0u) {
                    widths.push_back(32u);
                }
                else {
                    widths.push_back(64u);
                }
                if (K <= 96u && seed_index == 0u) {
                    widths.push_back(1354u);
                }
                for (uint32_t block_bytes : widths)
                {
                    const bool decode_variants =
                        block_bytes <= 64u &&
                        (full || block_bytes != 6u);
                    if (!RunCase(
                            profile, K, construction_seed, block_bytes,
                            decode_variants, counters))
                    {
                        return 1;
                    }
                }
            }
        }
    }
    // Directly straddle every automatic-dispatch exclusion that the regular
    // width/K ladder could otherwise skip: the first valid even width above
    // bb=8, the first K above K*bb=3072 at bb=8, and an explicitly requested
    // residue-bucket implementation at an otherwise eligible small cell.
    const uint64_t seam_seed = UINT64_C(0xC0FFEE);
    if (!RunCase(
            kProfiles[1], 2u, seam_seed,
            TinyMixedFastPathMaxBlockBytes() + 2u, false, counters) ||
        !RunCase(
            kProfiles[1],
            TinyMixedFastPathMaxProductBytes() /
                TinyMixedFastPathMaxBlockBytes() + 1u,
            seam_seed, TinyMixedFastPathMaxBlockBytes(), false, counters))
    {
        return 1;
    }
    counters.AutomaticExclusionSeamCases += 2u;
    ProfileConfig explicit_buckets = kProfiles[1];
    explicit_buckets.Name = "mixed-p244-d4-explicit-buckets";
    explicit_buckets.CanonicalDispatch = false;
    explicit_buckets.BucketMode = MixedResidueBucketMode::Separate;
    if (!RunCase(
            explicit_buckets, 16u, seam_seed, 2u, false, counters) ||
        !SetMixedResidueBucketModeForTesting(
            MixedResidueBucketMode::Automatic))
    {
        return 1;
    }
    ++counters.AutomaticExclusionSeamCases;
    const uint64_t expected_encode_cases = full ? 47139u : 2998u;
    const uint64_t expected_dispatch_cases = full ? 12473u : 805u;
    if (counters.EncodeCases != expected_encode_cases ||
        counters.CanonicalDispatchEncodeCases != expected_dispatch_cases)
    {
        std::fprintf(stderr,
            "FAIL grid cardinality encode=%llu expected=%llu "
            "dispatch=%llu expected_dispatch=%llu\n",
            (unsigned long long)counters.EncodeCases,
            (unsigned long long)expected_encode_cases,
            (unsigned long long)counters.CanonicalDispatchEncodeCases,
            (unsigned long long)expected_dispatch_cases);
        return 1;
    }
    if (counters.AutomaticExclusionSeamCases != 3u)
    {
        std::fprintf(stderr,
            "FAIL automatic exclusion seam count=%u expected=3\n",
            counters.AutomaticExclusionSeamCases);
        return 1;
    }
    uint32_t required_dispatch_k = 0u;
    uint32_t covered_dispatch_k = 0u;
    for (uint32_t K : k_values)
    {
        if (K > 100u) {
            continue;
        }
        ++required_dispatch_k;
        if (!counters.CanonicalDispatchWidth2Success[K] ||
            !counters.CanonicalDispatchWidth4Success[K])
        {
            std::fprintf(stderr,
                "FAIL canonical dispatch lacks accepted successful "
                "bb=2/4 coverage at K=%u\n", K);
            return 1;
        }
        ++covered_dispatch_k;
    }
    std::printf(
        "tiny fast path A/B: %llu encode cases, %llu decode cases "
        "(%llu success, %llu need-more, %llu error) byte-identical; "
        "tiny_acceptances=%llu successful_forced_on_mixed=%llu "
        "forced_off_acceptances=%llu "
        "forced_on_nonmixed_acceptances=%llu "
        "automatic_mixed_acceptances=%llu "
        "successful_automatic_eligible_mixed=%llu "
        "automatic_ineligible_acceptances=%llu "
        "automatic_exclusion_seams=%u/3 "
        "canonical_dispatch_encode=%llu "
        "canonical_dispatch_k_covered=%u/%u\n",
        (unsigned long long)counters.EncodeCases,
        (unsigned long long)counters.DecodeCases,
        (unsigned long long)counters.SuccessCases,
        (unsigned long long)counters.NeedMoreCases,
        (unsigned long long)counters.ErrorCases,
        (unsigned long long)counters.ForcedOnMixedAcceptances,
        (unsigned long long)counters.SuccessfulForcedOnMixedSolves,
        (unsigned long long)counters.ForcedOffAcceptances,
        (unsigned long long)counters.ForcedOnNonMixedAcceptances,
        (unsigned long long)counters.AutomaticMixedAcceptances,
        (unsigned long long)
            counters.SuccessfulAutomaticEligibleMixedSolves,
        (unsigned long long)counters.AutomaticIneligibleAcceptances,
        counters.AutomaticExclusionSeamCases,
        (unsigned long long)counters.CanonicalDispatchEncodeCases,
        covered_dispatch_k, required_dispatch_k);
    return 0;
}
