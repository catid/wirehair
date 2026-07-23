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
    return true;
}

struct SolveOutcome
{
    WirehairResult Result = Wirehair_Error;
    std::vector<uint8_t> Blocks;
};

SolveOutcome RunSolve(
    int fast_path_mode,
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const std::vector<SolvePacket>& packets,
    uint32_t block_bytes)
{
    SolveOutcome outcome;
    if (!SetTinyMixedFastPathModeForTesting(fast_path_mode)) {
        outcome.Result = Wirehair_Error;
        return outcome;
    }
    outcome.Result = SolvePrecodeSystem(
        system, config, packets, block_bytes, outcome.Blocks);
    (void)SetTinyMixedFastPathModeForTesting(0);
    return outcome;
}

struct Counters
{
    uint64_t EncodeCases = 0;
    uint64_t DecodeCases = 0;
    uint64_t SuccessCases = 0;
    uint64_t NeedMoreCases = 0;
    uint64_t ErrorCases = 0;
    uint64_t FastPathEngaged = 0;
};

bool CompareOutcomes(
    const char* what,
    const ProfileConfig& profile,
    uint32_t K,
    uint64_t seed,
    uint32_t block_bytes,
    const SolveOutcome& off,
    const SolveOutcome& on,
    Counters& counters)
{
    if (off.Result != on.Result)
    {
        std::fprintf(stderr,
            "FAIL %s %s K=%u seed=0x%llx bb=%u result off=%d on=%d\n",
            what, profile.Name, K, (unsigned long long)seed, block_bytes,
            (int)off.Result, (int)on.Result);
        return false;
    }
    if (off.Result == Wirehair_Success)
    {
        ++counters.SuccessCases;
        if (off.Blocks != on.Blocks)
        {
            std::fprintf(stderr,
                "FAIL %s %s K=%u seed=0x%llx bb=%u intermediate "
                "blocks differ\n",
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
    const SolveOutcome off =
        RunSolve(-1, system, config, packets, block_bytes);
    const SolveOutcome on =
        RunSolve(1, system, config, packets, block_bytes);
    if (!CompareOutcomes(
            "encode", profile, K, construction_seed, block_bytes,
            off, on, counters))
    {
        return false;
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
    if (!CompareOutcomes(
            "decode", profile, K, construction_seed, block_bytes,
            decode_off, decode_on, counters))
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
        if (!CompareOutcomes(
                "deficient", profile, K, construction_seed, block_bytes,
                deficient_off, deficient_on, counters))
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
        if (!CompareOutcomes(
                "conflict", profile, K, construction_seed, block_bytes,
                conflict_off, conflict_on, counters))
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
        {"mixed-p244-baseline", true, false, 0u, 0u, 0u, false, 2u};
    if (profile_name && !std::strcmp(profile_name, "d4")) {
        profile.Name = "mixed-p244-d4";
        profile.DenseRows = 4u;
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
    std::printf("tiny scalar helpers: exhaustive inverse and row kernels "
                "verified\n");

    static const ProfileConfig kProfiles[] = {
        // mixed-p244-baseline: frozen geometry defaults
        {"mixed-p244-baseline", true, false, 0u, 0u, 0u, false, 2u},
        // ic-d8-9x2-p48x: shared-x, period 48, 9 GF256 rows, D2=8,
        // identity corner (structurally unavailable below K=5)
        {"ic-d8-9x2-p48x", true, true, 48u, 9u, 8u, true, 5u},
        // certified control: never touches the mixed fast path
        {"cert-baseline", false, false, 0u, 0u, 0u, false, 2u},
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
    std::printf(
        "tiny fast path A/B: %llu encode cases, %llu decode cases "
        "(%llu success, %llu need-more, %llu error) byte-identical\n",
        (unsigned long long)counters.EncodeCases,
        (unsigned long long)counters.DecodeCases,
        (unsigned long long)counters.SuccessCases,
        (unsigned long long)counters.NeedMoreCases,
        (unsigned long long)counters.ErrorCases);
    return 0;
}
