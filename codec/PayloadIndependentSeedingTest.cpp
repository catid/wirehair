// Unit tests for the experiment-only payload-independent WH2 equation
// seeding profile (wirehair-sxvz.16.1.5.2).
//
// The profile is a test-hook-gated, thread-local knob: while enabled,
// MatrixSeedFromProfile() hashes the seed profile normalized to
// kPayloadIndependentEquationSeedBlockBytes (= 2) instead of the true
// payload width, while memory layout and kernels keep the true BlockBytes.
// These tests pin the pure normalization contract, the cross-width seed
// invariance, the bb=2 identity that makes K-indexed fixups portable, and a
// cheap end-to-end roundtrip through the real message-precode codec.

#include "WirehairV2Codec.h"
#include "WirehairV2Plan.h"
#include "WirehairV2Policy.h"
#include "WirehairV2Seeds.h"
#include "WirehairV2Solve.h"

#include <cstdio>
#include <cstring>
#include <vector>

namespace {

int Failures = 0;

void Check(bool condition, const char* message)
{
    if (!condition) {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++Failures;
    }
}

const uint32_t kBlockCounts[] = {
    2u, 16u, 96u, 500u, 1000u, 3200u, 10000u, 32000u, 64000u
};

// The matched evaluation strata from the cross-width holdout.
const uint32_t kBlockByteStrata[] = { 2u, 64u, 256u, 1280u, 4096u };

const uint64_t kSalts[] = {
    0u,
    UINT64_C(0x5eedf411),
    wirehair_v2::kMessagePrecodeSeedSalt,
    wirehair_v2::kMessageRecoveryRowSeedSalt
};

bool SamePeelingCodec(
    const wirehair_v2::PeelingCodec& a,
    const wirehair_v2::PeelingCodec& b)
{
    return a.Solver == b.Solver &&
        a.Structure == b.Structure &&
        a.Family == b.Family &&
        a.MinDegree == b.MinDegree &&
        a.MaxDegree == b.MaxDegree &&
        a.SolverCandidateLimit == b.SolverCandidateLimit &&
        a.Degree1Mass == b.Degree1Mass &&
        a.Degree2Mass == b.Degree2Mass &&
        a.RobustC == b.RobustC &&
        a.RobustDelta == b.RobustDelta &&
        a.FullyRandomRows == b.FullyRandomRows &&
        a.UseWirehairRowDistribution == b.UseWirehairRowDistribution;
}

bool SamePolicy(
    const wirehair_v2::PeelPolicy& a,
    const wirehair_v2::PeelPolicy& b)
{
    return a.Solver == b.Solver &&
        a.Structure == b.Structure &&
        a.ByteClass == b.ByteClass &&
        a.CountBand == b.CountBand &&
        SamePeelingCodec(a.Codec, b.Codec);
}

void CheckCanonicalConstant()
{
    Check(wirehair_v2::kPayloadIndependentEquationSeedBlockBytes == 2u,
        "the canonical equation-seed width must stay pinned to the "
        "validated bb=2 normalization");
}

void CheckNormalizationIsPure()
{
    for (uint32_t K : kBlockCounts)
    {
        const wirehair_v2::SeedProfile canonical =
            wirehair_v2::SelectSeedProfile(
                K, wirehair_v2::kPayloadIndependentEquationSeedBlockBytes);
        for (uint32_t bb : kBlockByteStrata)
        {
            const wirehair_v2::SeedProfile profile =
                wirehair_v2::SelectSeedProfile(K, bb);
            const wirehair_v2::SeedProfile normalized =
                wirehair_v2::NormalizeProfileForEquationSeeding(
                    profile,
                    wirehair_v2::kPayloadIndependentEquationSeedBlockBytes);

            Check(profile.BlockBytes == bb,
                "normalization must not mutate its input profile");
            Check(normalized.BlockBytes ==
                    wirehair_v2::kPayloadIndependentEquationSeedBlockBytes,
                "normalized profile must carry the canonical block bytes");
            Check(SamePolicy(normalized.Policy, canonical.Policy),
                "normalized policy must equal the canonical-width policy");

            // Every other equation-seed input derives from BlockCount alone
            // and must be preserved bit for bit.
            Check(normalized.BlockCount == profile.BlockCount &&
                    normalized.DenseCount == profile.DenseCount &&
                    normalized.PeelSeed == profile.PeelSeed &&
                    normalized.DenseSeed == profile.DenseSeed &&
                    normalized.PeelSeedBucket == profile.PeelSeedBucket &&
                    normalized.UsedPeelFixup == profile.UsedPeelFixup &&
                    normalized.UsedDenseFixup == profile.UsedDenseFixup,
                "normalization must preserve all K-derived profile fields");

            // The canonical-width profile has identical K-derived fields, so
            // the normalized copy must hash exactly like it.
            for (uint64_t salt : kSalts)
            {
                Check(wirehair_v2::MatrixSeedFromProfile(
                        normalized, 0u, salt) ==
                        wirehair_v2::MatrixSeedFromProfile(
                            canonical, 0u, salt),
                    "normalized profile must hash like the canonical-width "
                    "profile at row_count=0");
                Check(wirehair_v2::MatrixSeedFromProfile(
                        normalized, 100u, salt) ==
                        wirehair_v2::MatrixSeedFromProfile(
                            canonical, 100u, salt),
                    "normalized profile must hash like the canonical-width "
                    "profile at nonzero row_count");
            }

            // Zero canonical width is the documented no-op ("off").
            const wirehair_v2::SeedProfile untouched =
                wirehair_v2::NormalizeProfileForEquationSeeding(profile, 0u);
            Check(untouched.BlockBytes == profile.BlockBytes &&
                    SamePolicy(untouched.Policy, profile.Policy),
                "zero canonical width must be a no-op");
        }
    }
}

void CheckBaselinePayloadDependence()
{
    // Documents the effect under study: without the experiment profile, the
    // graph seed depends on the payload width even at matched K.  The hash
    // is deterministic, so these expectations are stable.
    for (uint32_t K : kBlockCounts)
    {
        const uint64_t canonical_seed = wirehair_v2::MatrixSeedFromProfile(
            wirehair_v2::SelectSeedProfile(
                K, wirehair_v2::kPayloadIndependentEquationSeedBlockBytes),
            0u, wirehair_v2::kMessagePrecodeSeedSalt);
        for (uint32_t bb : kBlockByteStrata)
        {
            if (bb == wirehair_v2::kPayloadIndependentEquationSeedBlockBytes) {
                continue;
            }
            Check(wirehair_v2::MatrixSeedFromProfile(
                    wirehair_v2::SelectSeedProfile(K, bb), 0u,
                    wirehair_v2::kMessagePrecodeSeedSalt) != canonical_seed,
                "baseline seeding is expected to depend on the payload "
                "width at matched K");
        }
    }
}

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)

struct ScopedPayloadIndependentSeeding
{
    explicit ScopedPayloadIndependentSeeding(bool enabled)
    {
        Saved = wirehair_v2::PayloadIndependentEquationSeedingForTesting();
        wirehair_v2::SetPayloadIndependentEquationSeedingForTesting(enabled);
    }
    ~ScopedPayloadIndependentSeeding()
    {
        wirehair_v2::SetPayloadIndependentEquationSeedingForTesting(Saved);
    }
    bool Saved;
};

void CheckHookNormalizesAcrossWidths()
{
    Check(!wirehair_v2::PayloadIndependentEquationSeedingForTesting(),
        "the payload-independent seeding hook must default to off");

    for (uint32_t K : kBlockCounts)
    {
        // Baseline values with the hook off.
        const wirehair_v2::SeedProfile canonical =
            wirehair_v2::SelectSeedProfile(
                K, wirehair_v2::kPayloadIndependentEquationSeedBlockBytes);
        const uint64_t canonical_matrix_seed =
            wirehair_v2::MatrixSeedFromProfile(
                canonical, 0u, wirehair_v2::kMessagePrecodeSeedSalt);
        const uint64_t canonical_row_seed =
            wirehair_v2::MatrixSeedFromProfile(
                canonical, 7u, wirehair_v2::kMessagePrecodeSeedSalt);
        const uint32_t canonical_packet_seed =
            wirehair_v2::PacketPeelSeedFromProfile(
                canonical, wirehair_v2::kMessageRecoveryRowSeedSalt);
        const uint64_t baseline_wide_seed =
            wirehair_v2::MatrixSeedFromProfile(
                wirehair_v2::SelectSeedProfile(K, 1280u), 0u,
                wirehair_v2::kMessagePrecodeSeedSalt);

        {
            ScopedPayloadIndependentSeeding enabled(true);
            for (uint32_t bb : kBlockByteStrata)
            {
                const wirehair_v2::SeedProfile profile =
                    wirehair_v2::SelectSeedProfile(K, bb);
                // Layout inputs are untouched by the hook: the profile still
                // carries the true payload width.
                Check(profile.BlockBytes == bb,
                    "hook must not change SelectSeedProfile block bytes");
                Check(wirehair_v2::MatrixSeedFromProfile(
                        profile, 0u,
                        wirehair_v2::kMessagePrecodeSeedSalt) ==
                        canonical_matrix_seed,
                    "hooked matrix seed must match the canonical-width "
                    "baseline at every stratum");
                Check(wirehair_v2::MatrixSeedFromProfile(
                        profile, 7u,
                        wirehair_v2::kMessagePrecodeSeedSalt) ==
                        canonical_row_seed,
                    "hooked matrix seed must match the canonical-width "
                    "baseline at nonzero row_count");
                Check(wirehair_v2::PacketPeelSeedFromProfile(
                        profile,
                        wirehair_v2::kMessageRecoveryRowSeedSalt) ==
                        canonical_packet_seed,
                    "hooked packet peel seed must match the canonical-width "
                    "baseline at every stratum");
            }
        }

        // The scope guard restored the hook: baseline behavior returns.
        Check(!wirehair_v2::PayloadIndependentEquationSeedingForTesting(),
            "scope exit must restore the disabled hook");
        Check(wirehair_v2::MatrixSeedFromProfile(
                wirehair_v2::SelectSeedProfile(K, 1280u), 0u,
                wirehair_v2::kMessagePrecodeSeedSalt) == baseline_wide_seed,
            "disabling the hook must restore payload-dependent seeding");
    }
}

bool RunHookedRoundTrip(
    uint32_t K,
    uint32_t block_bytes,
    bool mixed_profile,
    uint64_t& precode_seed_out,
    uint32_t& packet_peel_seed_out,
    uint32_t& seed_attempt_out)
{
    const uint64_t message_bytes = (uint64_t)K * block_bytes;
    std::vector<uint8_t> message((size_t)message_bytes);
    for (size_t i = 0; i < message.size(); ++i) {
        message[i] = (uint8_t)(i * 131u + K * 17u + block_bytes * 3u + 1u);
    }

    wirehair_v2::Codec encoder;
    wirehair_v2::Codec decoder;
    wirehair_v2::MessagePrecodeEncoderOptions options;
    if (mixed_profile) {
        options.Completion = wirehair_v2::CompletionField::MixedGF256GF16;
    }
    if (encoder.InitializePrecodeEncoder(
            message.data(), message_bytes, block_bytes,
            nullptr, &options) != Wirehair_Success ||
        decoder.InitializePrecodeDecoder(
            message_bytes, block_bytes,
            &encoder.Profile()) != Wirehair_Success)
    {
        std::fprintf(stderr,
            "FAIL: hooked roundtrip init K=%u bb=%u mixed=%d\n",
            K, block_bytes, mixed_profile ? 1 : 0);
        return false;
    }
    precode_seed_out = encoder.Profile().V2PrecodeSeed;
    packet_peel_seed_out = encoder.Profile().V2PacketPeelSeed;
    seed_attempt_out = encoder.Profile().V2SeedAttempt;

    WirehairResult decode_result = Wirehair_NeedMore;
    std::vector<uint8_t> block(block_bytes);
    // Drop every fifth source packet, then repair until solved.
    for (uint32_t id = 0; id < K && decode_result == Wirehair_NeedMore; ++id)
    {
        if (id % 5u == 0u) {
            continue;
        }
        uint32_t bytes = 0u;
        if (encoder.Encode(id, block.data(), block_bytes, &bytes) !=
            Wirehair_Success)
        {
            std::fprintf(stderr,
                "FAIL: hooked roundtrip encode id=%u K=%u bb=%u\n",
                id, K, block_bytes);
            return false;
        }
        decode_result = decoder.Decode(id, block.data(), bytes);
        if (decode_result != Wirehair_NeedMore &&
            decode_result != Wirehair_Success)
        {
            std::fprintf(stderr,
                "FAIL: hooked roundtrip decode id=%u result=%d\n",
                id, (int)decode_result);
            return false;
        }
    }
    for (uint32_t id = K;
        decode_result == Wirehair_NeedMore && id < K + 4u * K + 64u; ++id)
    {
        uint32_t bytes = 0u;
        if (encoder.Encode(id, block.data(), block_bytes, &bytes) !=
            Wirehair_Success)
        {
            std::fprintf(stderr,
                "FAIL: hooked roundtrip repair encode id=%u\n", id);
            return false;
        }
        decode_result = decoder.Decode(id, block.data(), bytes);
        if (decode_result != Wirehair_NeedMore &&
            decode_result != Wirehair_Success)
        {
            std::fprintf(stderr,
                "FAIL: hooked roundtrip repair decode id=%u result=%d\n",
                id, (int)decode_result);
            return false;
        }
    }
    if (decode_result != Wirehair_Success) {
        std::fprintf(stderr,
            "FAIL: hooked roundtrip never completed K=%u bb=%u\n",
            K, block_bytes);
        return false;
    }

    std::vector<uint8_t> recovered((size_t)message_bytes, 0xccu);
    if (decoder.Recover(recovered.data(), message_bytes) !=
            Wirehair_Success ||
        std::memcmp(recovered.data(), message.data(),
            (size_t)message_bytes) != 0)
    {
        std::fprintf(stderr,
            "FAIL: hooked roundtrip recover mismatch K=%u bb=%u\n",
            K, block_bytes);
        return false;
    }
    return true;
}

void CheckHookedCodecRoundTrip()
{
    const uint32_t K = 64u;
    for (int mixed = 0; mixed < 2; ++mixed)
    {
        // Baseline equation state at the canonical width with the hook OFF.
        uint64_t base_precode_seed = 0u;
        uint32_t base_packet_seed = 0u;
        uint32_t base_attempt = 0u;
        Check(RunHookedRoundTrip(
                K, wirehair_v2::kPayloadIndependentEquationSeedBlockBytes,
                mixed != 0, base_precode_seed, base_packet_seed,
                base_attempt),
            "baseline canonical-width roundtrip must succeed");

        ScopedPayloadIndependentSeeding enabled(true);
        // Both an even non-canonical width (mixed-compatible) and the
        // canonical width itself.
        const uint32_t widths[] = {
            wirehair_v2::kPayloadIndependentEquationSeedBlockBytes, 64u
        };
        for (uint32_t bb : widths)
        {
            uint64_t precode_seed = 0u;
            uint32_t packet_seed = 0u;
            uint32_t attempt = 0u;
            Check(RunHookedRoundTrip(
                    K, bb, mixed != 0, precode_seed, packet_seed, attempt),
                "hooked roundtrip must decode at every width");
            Check(precode_seed == base_precode_seed &&
                    packet_seed == base_packet_seed &&
                    attempt == base_attempt,
                "hooked equation state must equal the canonical-width "
                "baseline, making K-indexed behavior portable");
        }
    }
}

#endif // WIREHAIR_V2_ENABLE_TEST_HOOKS

} // namespace

int main()
{
    CheckCanonicalConstant();
    CheckNormalizationIsPure();
    CheckBaselinePayloadDependence();
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    CheckHookNormalizesAcrossWidths();
    CheckHookedCodecRoundTrip();
#endif
    if (Failures != 0) {
        std::fprintf(stderr,
            "payload-independent seeding test failures: %d\n", Failures);
        return 1;
    }
    std::printf("payload-independent seeding tests passed\n");
    return 0;
}
