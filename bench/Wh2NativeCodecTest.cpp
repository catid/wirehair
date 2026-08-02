#include "Wh2NativeCodec.h"

#include <wirehair/wirehair.h>

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <numeric>
#include <vector>

namespace {

using wirehair_wh2_bench::IsolatedSolveResult;
using wirehair_wh2_bench::NativeArm;
using wirehair_wh2_bench::NativeArmKind;
using wirehair_wh2_bench::NativeArmSpec;
using wirehair_wh2_bench::NativeEncoderFixture;
using wirehair_wh2_bench::NativeReceiveFixture;
using wirehair_wh2_bench::NativeRecoveryFixture;
using wirehair_wh2_bench::NativeSolveFixture;
using wirehair_wh2_bench::RecoveryCellResult;
using wirehair_wh2_bench::RecoveryOutcome;

int Failures = 0;

bool Check(bool condition, const char* message)
{
    if (!condition)
    {
        std::fprintf(stderr, "FAIL: %s\n", message);
        ++Failures;
    }
    return condition;
}

std::vector<uint32_t> ConsecutiveIds(uint32_t K, uint32_t first = 0u)
{
    std::vector<uint32_t> ids((size_t)K + 4u);
    std::iota(ids.begin(), ids.end(), first);
    return ids;
}

bool EncodePublicWh2(
    const std::vector<uint8_t>& source,
    uint32_t block_bytes,
    uint32_t attempt,
    const std::vector<uint32_t>& packet_ids,
    std::vector<uint8_t>& packets_out,
    WirehairV2Result& create_result_out)
{
    WirehairV2Profile profile = {};
    profile.struct_bytes = (uint32_t)sizeof(profile);
    profile.profile_version = WIREHAIR_V2_PROFILE_VERSION;
    profile.profile_id = WIREHAIR_V2_PROFILE_CERTIFIED_2026_07;
    profile.message_bytes = source.size();
    profile.block_bytes = block_bytes;
    profile.seed_attempt = (uint8_t)attempt;

    uint8_t serialized[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES] = {};
    uint32_t serialized_bytes = 0u;
    if (wirehair_v2_profile_serialize(
            &profile,
            serialized,
            sizeof(serialized),
            &serialized_bytes) != WirehairV2_Success ||
        serialized_bytes != sizeof(serialized))
    {
        return false;
    }

    WirehairV2Codec encoder = nullptr;
    create_result_out = wirehair_v2_encoder_create_profile(
        source.data(), serialized, sizeof(serialized), &encoder);
    if (create_result_out != WirehairV2_Success) {
        return true;
    }

    std::vector<uint8_t> packets(packet_ids.size() * block_bytes, 0u);
    bool valid = true;
    for (size_t i = 0u; i < packet_ids.size(); ++i)
    {
        uint32_t written = 0u;
        if (wirehair_v2_encode(
                encoder,
                packet_ids[i],
                packets.data() + i * block_bytes,
                block_bytes,
                &written) != WirehairV2_Success ||
            written != block_bytes)
        {
            valid = false;
            break;
        }
    }
    wirehair_v2_free(encoder);
    if (valid) {
        packets_out.swap(packets);
    }
    return valid;
}

void CheckDeterministicSource()
{
    std::vector<uint8_t> source;
    const uint8_t expected[] = {
        0x67, 0xf8, 0x00, 0xa1, 0x2a,
        0x6e, 0xe2, 0xbb, 0x29, 0xb7,
        0xae, 0x62, 0x33, 0xf7, 0x68
    };
    Check(wirehair_wh2_bench::MakeDeterministicSource(
              3u, 5u, UINT64_C(0x0123456789abcdef), source),
        "deterministic source generation failed");
    Check(source.size() == sizeof(expected) &&
              std::memcmp(source.data(), expected, sizeof(expected)) == 0,
        "deterministic source golden changed");

    const std::vector<uint8_t> before = source;
    Check(!wirehair_wh2_bench::MakeDeterministicSource(
              1u, 5u, 7u, source),
        "invalid deterministic source dimensions accepted");
    Check(source == before,
        "failed deterministic source generation modified output");
}

void CheckCertifiedMatchesPublicProfile()
{
    const uint32_t block_counts[] = { 2u, 3u, 4u, 5u, 17u, 48u, 257u };
    uint32_t total_successful_attempts = 0u;
    for (uint32_t K : block_counts)
    {
        const uint32_t block_bytes = 17u + K % 23u;
        std::vector<uint8_t> source;
        Check(wirehair_wh2_bench::MakeDeterministicSource(
                  K,
                  block_bytes,
                  UINT64_C(0x824d3b297c6e19a1) ^ K,
                  source),
            "public comparison source generation failed");
        const uint32_t packet_id_values[] = {
            0u, 1u, K - 1u, K, K + 7u, UINT32_MAX
        };
        const std::vector<uint32_t> packet_ids(
            packet_id_values,
            packet_id_values + sizeof(packet_id_values) /
                sizeof(packet_id_values[0]));

        for (uint32_t attempt = 0u; attempt < 3u; ++attempt)
        {
            NativeArm native;
            const WirehairResult native_result = native.Initialize(
                wirehair_wh2_bench::MakeCertifiedWh2Arm(attempt),
                K,
                block_bytes,
                source);

            std::vector<uint8_t> public_packets;
            WirehairV2Result public_result = WirehairV2_Error;
            Check(EncodePublicWh2(
                      source,
                      block_bytes,
                      attempt,
                      packet_ids,
                      public_packets,
                      public_result),
                "public explicit-profile packet generation failed");

            if (public_result != WirehairV2_Success &&
                public_result != WirehairV2_BadSeed &&
                public_result != WirehairV2_Error)
            {
                std::fprintf(
                    stderr,
                    "public exact profile K=%u attempt=%u returned %s; native=%s\n",
                    K,
                    attempt,
                    wirehair_v2_result_string(public_result),
                    wirehair_result_string(native_result));
            }

            if (public_result == WirehairV2_BadSeed)
            {
                Check(native_result == Wirehair_BadPeelSeed ||
                          native_result == Wirehair_BadDenseSeed,
                    "native exact-attempt rank failure disagrees with public profile");
                continue;
            }
            if (public_result == WirehairV2_Error)
            {
                Check(native_result == Wirehair_Error,
                    "native exact-attempt fatal failure disagrees with public profile");
                continue;
            }
            Check(public_result == WirehairV2_Success,
                "public exact profile returned an unexpected result");
            if (!Check(native_result == Wirehair_Success,
                    "native exact profile rejected a public-success attempt"))
            {
                continue;
            }
            ++total_successful_attempts;

            std::vector<uint8_t> native_packets;
            for (uint32_t block_id : packet_ids)
            {
                std::vector<uint8_t> packet;
                Check(native.Encode(block_id, packet) == Wirehair_Success,
                    "native exact profile encode failed");
                native_packets.insert(
                    native_packets.end(), packet.begin(), packet.end());
            }
            Check(native_packets == public_packets,
                "native exact profile packets differ from serialized public profile");
        }
    }
    Check(total_successful_attempts != 0u,
        "no public-success attempt was available for exact-profile comparison");
}

void CheckSourceIndependentRankClassification()
{
    const uint32_t K = 8u;
    const uint32_t block_bytes = 64u;
    const uint32_t attempt = 2u;
    std::vector<uint8_t> timing_source;
    Check(wirehair_wh2_bench::MakeDeterministicSource(
              K,
              block_bytes,
              UINT64_C(0xb686a5d9065d521e),
              timing_source),
        "rank-classification timing source generation failed");
    const std::vector<uint8_t> zero_source(
        (size_t)K * block_bytes, 0u);
    const std::vector<std::vector<uint8_t> > sources = {
        timing_source,
        zero_source
    };

    for (const std::vector<uint8_t>& source : sources)
    {
        NativeArm native;
        Check(native.Initialize(
                  wirehair_wh2_bench::MakeCertifiedWh2Arm(attempt),
                  K,
                  block_bytes,
                  source) == Wirehair_BadPeelSeed,
            "native exact rank failure depended on source RHS");

        std::vector<uint8_t> public_packets;
        WirehairV2Result public_result = WirehairV2_Error;
        Check(EncodePublicWh2(
                  source,
                  block_bytes,
                  attempt,
                  std::vector<uint32_t>(),
                  public_packets,
                  public_result),
            "public exact rank-classification setup failed");
        Check(public_result == WirehairV2_BadSeed,
            "public exact rank failure depended on source RHS");
    }
}

struct TransformState
{
    uint32_t Calls = 0u;
};

bool RuntimeCandidateTransform(
    uint32_t block_count,
    uint32_t block_bytes,
    wirehair_v2::PrecodeParams& params,
    wirehair_v2::PacketRowConfig& packet_config,
    void* context)
{
    TransformState* state = static_cast<TransformState*>(context);
    if (!state || block_count != params.BlockCount || block_bytes == 0u) {
        return false;
    }
    ++state->Calls;
    params.DenseIdentityCorner = true;
    params.HeavyFamily =
        wirehair_v2::HeavyCoefficientFamily::HashedNonzero;
    params.Seed ^= UINT64_C(0x53a9d2e17c4086bf);
    packet_config.PeelSeed ^= UINT32_C(0x6d2b79f5);
    return true;
}

bool InvalidCandidateTransform(
    uint32_t,
    uint32_t,
    wirehair_v2::PrecodeParams& params,
    wirehair_v2::PacketRowConfig&,
    void*)
{
    params.Staircase = 0u;
    return true;
}

void CheckArmsAndNestedRecovery()
{
    const uint32_t K = 32u;
    const uint32_t block_bytes = 23u;
    std::vector<uint8_t> source;
    Check(wirehair_wh2_bench::MakeDeterministicSource(
              K, block_bytes, UINT64_C(0xe7d4a10362c9b85f), source),
        "arm recovery source generation failed");
    const std::vector<uint32_t> ids = ConsecutiveIds(K);

    const NativeArmSpec controls[] = {
        wirehair_wh2_bench::MakeWirehair1Arm(),
        wirehair_wh2_bench::MakeCertifiedWh2Arm(0u)
    };
    for (const NativeArmSpec& spec : controls)
    {
        const RecoveryCellResult recovery =
            wirehair_wh2_bench::RunRecoveryCell(
                spec, K, block_bytes, source, ids);
        Check(recovery.Outcome == RecoveryOutcome::Success &&
                  recovery.FirstOverhead == 0u &&
                  recovery.Result == Wirehair_Success,
            "control failed exact nested systematic recovery");
    }

    TransformState transform_state;
    const NativeArmSpec candidate =
        wirehair_wh2_bench::MakeExperimentalWh2Arm(
            1u, RuntimeCandidateTransform, &transform_state);
    NativeArm candidate_arm;
    Check(candidate_arm.Initialize(
              candidate, K, block_bytes, source) == Wirehair_Success,
        "runtime candidate exact construction failed");
    Check(candidate_arm.Kind() == NativeArmKind::Wirehair2Experiment &&
              candidate_arm.ConstructionAttempt() == 1u &&
              transform_state.Calls == 1u,
        "runtime candidate transform/identity receipt is wrong");

    NativeRecoveryFixture fixture;
    Check(fixture.Initialize(candidate_arm, ids) == Wirehair_Success,
        "runtime candidate recovery fixture failed");
    const RecoveryCellResult candidate_recovery = fixture.RunNested();
    Check(candidate_recovery.Outcome == RecoveryOutcome::Success &&
              candidate_recovery.FirstOverhead == 0u,
        "runtime candidate failed byte-verified nested recovery");

    std::vector<uint32_t> duplicate_ids = ids;
    duplicate_ids.back() = duplicate_ids.front();
    NativeRecoveryFixture duplicate_fixture;
    Check(duplicate_fixture.Initialize(candidate_arm, duplicate_ids) ==
              Wirehair_InvalidInput,
        "duplicate delivered packet id was accepted");

    NativeArm invalid_candidate;
    const WirehairResult invalid_result = invalid_candidate.Initialize(
        wirehair_wh2_bench::MakeExperimentalWh2Arm(
            0u, InvalidCandidateTransform),
        K,
        block_bytes,
        source);
    Check(invalid_result == Wirehair_InvalidInput,
        "invalid runtime equation transform was not rejected");
}

void CheckTransactionalArmAndValidation()
{
    const uint32_t K = 24u;
    const uint32_t block_bytes = 31u;
    std::vector<uint8_t> source;
    Check(wirehair_wh2_bench::MakeDeterministicSource(
              K, block_bytes, UINT64_C(0xb36f29841c05e7da), source),
        "transactional source generation failed");

    NativeArm arm;
    Check(arm.Initialize(
              wirehair_wh2_bench::MakeCertifiedWh2Arm(0u),
              K,
              block_bytes,
              source) == Wirehair_Success,
        "transactional arm initial construction failed");
    std::vector<uint8_t> before;
    Check(arm.Encode(K + 5u, before) == Wirehair_Success,
        "transactional arm initial encode failed");

    Check(arm.Initialize(
              wirehair_wh2_bench::MakeCertifiedWh2Arm(256u),
              K,
              block_bytes,
              source) == Wirehair_InvalidInput,
        "out-of-range construction attempt accepted");
    std::vector<uint8_t> after;
    Check(arm.IsInitialized() &&
              arm.ConstructionAttempt() == 0u &&
              arm.Encode(K + 5u, after) == Wirehair_Success &&
              after == before,
        "failed reinitialization modified the prior exact arm");

    NativeArmSpec bad_legacy = wirehair_wh2_bench::MakeWirehair1Arm();
    bad_legacy.ConstructionAttempt = 1u;
    NativeArm legacy;
    Check(legacy.Initialize(bad_legacy, K, block_bytes, source) ==
              Wirehair_InvalidInput,
        "Wirehair1 accepted a construction attempt");

    NativeArmSpec bad_control =
        wirehair_wh2_bench::MakeCertifiedWh2Arm(0u);
    bad_control.Transform = RuntimeCandidateTransform;
    TransformState state;
    bad_control.TransformContext = &state;
    NativeArm control;
    Check(control.Initialize(bad_control, K, block_bytes, source) ==
              Wirehair_InvalidInput && state.Calls == 0u,
        "certified control accepted an experimental transform");

    NativeArmSpec cached_candidate =
        wirehair_wh2_bench::MakeExperimentalWh2Arm(
            0u, RuntimeCandidateTransform, &state);
    cached_candidate.Wh2Options.CacheSystematicSource = true;
    NativeArm cached;
    Check(cached.Initialize(cached_candidate, K, block_bytes, source) ==
              Wirehair_InvalidInput,
        "direct adapter silently accepted a facade cache timing axis");

    NativeArmSpec context_without_transform =
        wirehair_wh2_bench::MakeExperimentalWh2Arm(0u, nullptr, &state);
    NativeArm missing_transform;
    Check(missing_transform.Initialize(
              context_without_transform, K, block_bytes, source) ==
              Wirehair_InvalidInput,
        "experimental context without a transform was accepted");

    NativeEncoderFixture failed_encoder;
    Check(failed_encoder.Initialize(
              wirehair_wh2_bench::MakeExperimentalWh2Arm(
                  0u, InvalidCandidateTransform),
              K,
              block_bytes,
              source) == Wirehair_Success,
        "invalid candidate could not reach its timed construction");
    const wirehair_wh2_bench::TimedArmResult failed_timing =
        failed_encoder.Run();
    Check(failed_timing.Result == Wirehair_InvalidInput &&
              failed_timing.ElapsedNanoseconds == 0u &&
              !failed_timing.BytesVerified,
        "failed encoder timing did not clear its duration");
}

void CheckIsolatedSolveFixture()
{
    const uint32_t K = 40u;
    const uint32_t block_bytes = 64u;
    std::vector<uint8_t> source;
    Check(wirehair_wh2_bench::MakeDeterministicSource(
              K, block_bytes, UINT64_C(0x24c7e91b65ad803f), source),
        "isolated solve source generation failed");
    const std::vector<uint32_t> ids = ConsecutiveIds(K);

    NativeArm arm;
    Check(arm.Initialize(
              wirehair_wh2_bench::MakeCertifiedWh2Arm(0u),
              K,
              block_bytes,
              source) == Wirehair_Success,
        "isolated solve exact construction failed");
    NativeSolveFixture fixture;
    Check(fixture.Initialize(arm, ids, 4u) == Wirehair_Success,
        "isolated solve fixture preparation failed");
    for (uint32_t run = 0u; run < 2u; ++run)
    {
        const IsolatedSolveResult result = fixture.Run();
        Check(result.Result == Wirehair_Success &&
                  result.BytesVerified &&
                  result.ElapsedNanoseconds > 0u &&
                  result.Stats.PacketRows == K + 4u,
            "isolated solve invocation was not timed and byte-verified");
    }

    NativeArm wh1;
    Check(wh1.Initialize(
              wirehair_wh2_bench::MakeWirehair1Arm(),
              K,
              block_bytes,
              source) == Wirehair_Success,
        "Wirehair1 solve-fixture rejection setup failed");
    NativeSolveFixture rejected;
    Check(rejected.Initialize(wh1, ids, 4u) == Wirehair_InvalidInput,
        "isolated solve incorrectly accepted Wirehair1");
    Check(rejected.Initialize(arm, ids, 5u) == Wirehair_InvalidInput,
        "isolated solve accepted overhead above the frozen cap");
}

void CheckOtherTimingScopes()
{
    const uint32_t K = 36u;
    const uint32_t block_bytes = 64u;
    std::vector<uint8_t> source;
    Check(wirehair_wh2_bench::MakeDeterministicSource(
              K, block_bytes, UINT64_C(0x943a0f6de5b728c1), source),
        "other timing-scope source generation failed");
    const std::vector<uint32_t> ids = ConsecutiveIds(K);

    TransformState transform_state;
    const NativeArmSpec specs[] = {
        wirehair_wh2_bench::MakeWirehair1Arm(),
        wirehair_wh2_bench::MakeCertifiedWh2Arm(0u),
        wirehair_wh2_bench::MakeExperimentalWh2Arm(
            1u, RuntimeCandidateTransform, &transform_state)
    };
    for (const NativeArmSpec& spec : specs)
    {
        const uint32_t calls_before_prepare = transform_state.Calls;
        NativeEncoderFixture encoder_fixture;
        Check(encoder_fixture.Initialize(
                  spec, K, block_bytes, source) == Wirehair_Success,
            "encoder timing fixture preflight failed");
        Check(transform_state.Calls == calls_before_prepare,
            "encoder fixture preparation constructed the timed candidate");
        const wirehair_wh2_bench::TimedArmResult encoder_result =
            encoder_fixture.Run();
        Check(encoder_result.Result == Wirehair_Success &&
                  encoder_result.BytesVerified &&
                  encoder_result.ElapsedNanoseconds > 0u &&
                  encoder_result.DecodedOverhead == UINT32_MAX,
            "fresh encoder timing scope failed");
        const uint32_t expected_transform_calls =
            calls_before_prepare +
            (spec.Kind == NativeArmKind::Wirehair2Experiment ? 1u : 0u);
        Check(transform_state.Calls == expected_transform_calls,
            "fresh encoder timing did not construct its arm exactly once");

        NativeArm arm;
        Check(arm.Initialize(spec, K, block_bytes, source) == Wirehair_Success,
            "receive timing arm construction failed");
        NativeReceiveFixture receive_fixture;
        Check(receive_fixture.Initialize(arm, ids) == Wirehair_Success,
            "receive timing fixture preparation failed");
        const wirehair_wh2_bench::TimedArmResult receive_result =
            receive_fixture.Run();
        Check(receive_result.Result == Wirehair_Success &&
                  receive_result.BytesVerified &&
                  receive_result.ElapsedNanoseconds > 0u &&
                  receive_result.DecodedOverhead == 0u,
            "fresh receive-to-success timing scope failed");
    }
}

} // namespace

int main()
{
    Check(wirehair_init() == Wirehair_Success,
        "wirehair initialization failed");
    CheckDeterministicSource();
    CheckCertifiedMatchesPublicProfile();
    CheckSourceIndependentRankClassification();
    CheckArmsAndNestedRecovery();
    CheckTransactionalArmAndValidation();
    CheckIsolatedSolveFixture();
    CheckOtherTimingScopes();
    if (Failures != 0)
    {
        std::fprintf(stderr, "%d native codec test(s) failed\n", Failures);
        return 1;
    }
    std::puts("WH2 native codec tests passed");
    return 0;
}
