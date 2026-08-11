#include "Wh2NativeCodec.h"

#include "../codec/WirehairV2PrecodeDecode.h"

#include <wirehair/wirehair.h>

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <map>
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
using wirehair_wh2_bench::NativeTimingControlProbe;
using wirehair_wh2_bench::NativeTimingControlQualification;
using wirehair_wh2_bench::NativeTimingControlQualificationResult;
using wirehair_wh2_bench::NativeWh2BaseKind;
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

std::vector<uint32_t> ConsecutiveIds(
    uint32_t K,
    uint32_t first = 0u,
    uint32_t overhead = 4u)
{
    std::vector<uint32_t> ids((size_t)K + overhead);
    std::iota(ids.begin(), ids.end(), first);
    return ids;
}

bool FindFirstPacketRowCollision(
    const wirehair_v2::PrecodeSystem& system,
    const wirehair_v2::PacketRowConfig& config,
    uint32_t& first_id,
    uint32_t& second_id)
{
    const uint32_t K = system.Params.BlockCount;
    const uint32_t P = system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    std::map<std::vector<uint32_t>, uint32_t> seen;
    for (uint32_t id = K; id < K + 1000000u; ++id)
    {
        std::vector<uint32_t> row =
            wirehair_v2::GeneratePacketMatrixRow(K, P, id, config);
        std::sort(row.begin(), row.end());
        const std::pair<
            std::map<std::vector<uint32_t>, uint32_t>::iterator,
            bool> inserted = seen.insert(std::make_pair(row, id));
        if (!inserted.second)
        {
            first_id = inserted.first->second;
            second_id = id;
            return true;
        }
    }
    return false;
}

struct IidTrace
{
    std::vector<uint32_t> Ids;
    uint64_t AttemptedCandidates = 0u;
};

uint64_t NextSplitMix64(uint64_t& state)
{
    state += UINT64_C(0x9e3779b97f4a7c15);
    uint64_t value = state;
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

IidTrace FrozenIidIds(
    uint32_t K,
    uint32_t block_bytes,
    uint64_t loss_seed,
    uint32_t delivered_overhead)
{
    IidTrace trace;
    trace.Ids.reserve(static_cast<size_t>(K) + delivered_overhead);
    uint64_t state = loss_seed ^
        (static_cast<uint64_t>(K) * UINT64_C(0x9e3779b97f4a7c15)) ^
        (static_cast<uint64_t>(block_bytes) *
            UINT64_C(0xbf58476d1ce4e5b9));
    while (trace.Ids.size() <
        static_cast<size_t>(K) + delivered_overhead)
    {
        const uint32_t id =
            static_cast<uint32_t>(trace.AttemptedCandidates++);
        const double unit = static_cast<double>(NextSplitMix64(state) >> 11) /
            9007199254740992.0;
        if (unit >= 0.1) {
            trace.Ids.push_back(id);
        }
    }
    return trace;
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

struct PublicReceiveResult
{
    WirehairV2Result Result = WirehairV2_Error;
    uint32_t DecodedOverhead = UINT32_MAX;
    bool BytesVerified = false;
};

bool RunPublicWh2Receive(
    const NativeArm& arm,
    const std::vector<uint8_t>& source,
    uint32_t attempt,
    const std::vector<uint32_t>& packet_ids,
    uint32_t receive_overhead_cap,
    PublicReceiveResult& result_out)
{
    const uint64_t required =
        (uint64_t)arm.BlockCount() + receive_overhead_cap;
    if (!arm.IsInitialized() ||
        arm.Kind() != NativeArmKind::Wirehair2Certified ||
        attempt != arm.ConstructionAttempt() ||
        source.size() !=
            (size_t)arm.BlockCount() * arm.BlockBytes() ||
        required > packet_ids.size())
    {
        return false;
    }

    WirehairV2Profile profile = {};
    profile.struct_bytes = (uint32_t)sizeof(profile);
    profile.profile_version = WIREHAIR_V2_PROFILE_VERSION;
    profile.profile_id = WIREHAIR_V2_PROFILE_CERTIFIED_2026_07;
    profile.message_bytes = source.size();
    profile.block_bytes = arm.BlockBytes();
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

    WirehairV2Codec decoder = nullptr;
    if (wirehair_v2_decoder_create(
            serialized, sizeof(serialized), &decoder) !=
        WirehairV2_Success)
    {
        return false;
    }

    PublicReceiveResult result;
    std::vector<uint8_t> packet(arm.BlockBytes(), 0u);
    std::vector<uint8_t> recovered(source.size(), 0u);
    for (uint32_t packet_index = 0u;
         packet_index < required;
         ++packet_index)
    {
        if (arm.EncodeInto(
                packet_ids[packet_index],
                packet.data(),
                arm.BlockBytes()) != Wirehair_Success)
        {
            wirehair_v2_free(decoder);
            return false;
        }
        result.Result = wirehair_v2_decode(
            decoder,
            packet_ids[packet_index],
            packet.data(),
            arm.BlockBytes());
        if (result.Result == WirehairV2_Success)
        {
            result.DecodedOverhead =
                packet_index + 1u - arm.BlockCount();
            uint64_t recovered_bytes = 0u;
            result.Result = wirehair_v2_recover(
                decoder,
                recovered.data(),
                recovered.size(),
                &recovered_bytes);
            result.BytesVerified =
                result.Result == WirehairV2_Success &&
                recovered_bytes == source.size() && recovered == source;
            break;
        }
        if (result.Result != WirehairV2_NeedMore) {
            break;
        }
    }
    wirehair_v2_free(decoder);
    result_out = result;
    return true;
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

bool SamePrecodeParams(
    const wirehair_v2::PrecodeParams& left,
    const wirehair_v2::PrecodeParams& right)
{
    return left.BlockCount == right.BlockCount &&
        left.Staircase == right.Staircase &&
        left.DenseRows == right.DenseRows &&
        left.HeavyRows == right.HeavyRows &&
        left.SourceHits == right.SourceHits &&
        left.DenseIdentityCorner == right.DenseIdentityCorner &&
        left.DenseAnchors == right.DenseAnchors &&
        left.HeavyFamily == right.HeavyFamily &&
        left.Seed == right.Seed;
}

struct CanonicalBaseState
{
    uint32_t ExpectedBlockCount = 0u;
    uint32_t ExpectedBlockBytes = 0u;
    uint32_t Calls = 0u;
    bool SawCanonicalBase = false;
};

bool CanonicalBaseTransform(
    uint32_t block_count,
    uint32_t block_bytes,
    wirehair_v2::PrecodeParams& params,
    wirehair_v2::PacketRowConfig& packet_config,
    void* context)
{
    CanonicalBaseState* const state =
        static_cast<CanonicalBaseState*>(context);
    if (!state) {
        return false;
    }
    ++state->Calls;
    const wirehair_v2::PrecodeParams expected =
        wirehair_v2::MakeCertifiedParams(block_count, 0u);
    state->SawCanonicalBase =
        block_count == state->ExpectedBlockCount &&
        block_bytes == state->ExpectedBlockBytes &&
        SamePrecodeParams(params, expected) &&
        packet_config.PeelSeed == 0u &&
        packet_config.MixCount == wirehair_v2::kCertifiedPacketMixCount;
    if (!state->SawCanonicalBase) {
        return false;
    }

    // Install a fixed diagnostic schedule only after observing the closed,
    // seed-free base.  These values are intentionally independent of K and
    // block width and are known to admit the representative test shapes.
    params.Seed = UINT64_C(0x487468302aad7105);
    packet_config.PeelSeed = UINT32_C(0x4ec72102);
    return true;
}

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

bool InvalidDenseAnchorTransform(
    uint32_t,
    uint32_t,
    wirehair_v2::PrecodeParams& params,
    wirehair_v2::PacketRowConfig&,
    void*)
{
    params.DenseAnchors =
        static_cast<wirehair_v2::DenseAnchorLayout>(UINT32_MAX);
    return true;
}

enum class InvalidDenseAnchorCombination
{
    DenseRows,
    HeavyRows,
    IdentityCorner,
    HeavyFamily
};

bool InvalidDenseAnchorCombinationTransform(
    uint32_t,
    uint32_t,
    wirehair_v2::PrecodeParams& params,
    wirehair_v2::PacketRowConfig&,
    void* context)
{
    const InvalidDenseAnchorCombination* const combination =
        static_cast<const InvalidDenseAnchorCombination*>(context);
    if (!combination) {
        return false;
    }
    params.DenseAnchors = wirehair_v2::DenseAnchorLayout::Two07;
    switch (*combination)
    {
    case InvalidDenseAnchorCombination::DenseRows:
        params.DenseRows = 13u;
        break;
    case InvalidDenseAnchorCombination::HeavyRows:
        params.HeavyRows = 11u;
        break;
    case InvalidDenseAnchorCombination::IdentityCorner:
        params.DenseIdentityCorner = true;
        break;
    case InvalidDenseAnchorCombination::HeavyFamily:
        params.HeavyFamily =
            wirehair_v2::HeavyCoefficientFamily::HashedNonzero;
        break;
    }
    return true;
}

bool TwoAnchorTransform(
    uint32_t block_count,
    uint32_t block_bytes,
    wirehair_v2::PrecodeParams& params,
    wirehair_v2::PacketRowConfig&,
    void*)
{
    if (block_count != params.BlockCount || block_bytes == 0u ||
        params.DenseAnchors != wirehair_v2::DenseAnchorLayout::Disabled)
    {
        return false;
    }
    params.DenseAnchors = wirehair_v2::DenseAnchorLayout::Two07;
    return true;
}

void CheckCanonicalCertifiedBase()
{
    struct Shape
    {
        uint32_t K;
        uint32_t BlockBytes;
    };
    const Shape shapes[] = {
        { 4u, 2u },
        { 32u, 17u },
        { 1001u, 1280u }
    };
    for (const Shape& shape : shapes)
    {
        std::vector<uint8_t> source;
        Check(wirehair_wh2_bench::MakeDeterministicSource(
                  shape.K,
                  shape.BlockBytes,
                  UINT64_C(0xa12f05e46379bdc8) ^ shape.K ^
                      shape.BlockBytes,
                  source),
            "canonical-base source generation failed");
        CanonicalBaseState state;
        state.ExpectedBlockCount = shape.K;
        state.ExpectedBlockBytes = shape.BlockBytes;
        const NativeArmSpec spec =
            wirehair_wh2_bench::MakeExperimentalWh2Arm(
                0u,
                CanonicalBaseTransform,
                &state,
                nullptr,
                NativeWh2BaseKind::CanonicalCertifiedStructure);
        Check(spec.Kind == NativeArmKind::Wirehair2Experiment &&
                  spec.BaseKind ==
                      NativeWh2BaseKind::CanonicalCertifiedStructure,
            "canonical-base factory produced the wrong arm identity");

        NativeArm arm;
        Check(arm.Initialize(
                  spec, shape.K, shape.BlockBytes, source) == Wirehair_Success,
            "canonical-base exact construction failed");
        Check(state.Calls == 1u && state.SawCanonicalBase,
            "canonical transform did not observe the zero-seed certified base");
        std::vector<uint8_t> packet;
        Check(arm.Encode(shape.K + 7u, packet) == Wirehair_Success &&
                  packet.size() == shape.BlockBytes,
            "canonical-base arm did not remain usable after its transform");
    }

    // The timed encoder copies a spec and reconstructs the arm later through
    // the private initialization seam.  Exercise that path and prove that a
    // rejected reinitialization cannot replace the previously frozen spec.
    const Shape timed_shape = { 4u, 2u };
    std::vector<uint8_t> timed_source;
    Check(wirehair_wh2_bench::MakeDeterministicSource(
              timed_shape.K,
              timed_shape.BlockBytes,
              UINT64_C(0x6ad54e1c830fb972),
              timed_source),
        "canonical timed source generation failed");
    CanonicalBaseState timed_state;
    timed_state.ExpectedBlockCount = timed_shape.K;
    timed_state.ExpectedBlockBytes = timed_shape.BlockBytes;
    const NativeArmSpec timed_spec =
        wirehair_wh2_bench::MakeExperimentalWh2Arm(
            0u,
            CanonicalBaseTransform,
            &timed_state,
            nullptr,
            NativeWh2BaseKind::CanonicalCertifiedStructure);
    NativeEncoderFixture timed_fixture;
    Check(timed_fixture.Initialize(
              timed_spec,
              timed_shape.K,
              timed_shape.BlockBytes,
              timed_source) == Wirehair_Success &&
              timed_state.Calls == 0u,
        "canonical timed fixture constructed before its measured run");
    NativeArmSpec invalid_replacement = timed_spec;
    invalid_replacement.Wh2Options.PrecodeSeedSalt = 1u;
    Check(timed_fixture.Initialize(
              invalid_replacement,
              timed_shape.K,
              timed_shape.BlockBytes,
              timed_source) == Wirehair_InvalidInput &&
              timed_fixture.IsInitialized() &&
              timed_state.Calls == 0u,
        "invalid canonical timed reinitialization replaced the prior spec");
    const wirehair_wh2_bench::TimedArmResult timed_result =
        timed_fixture.Run();
    Check(timed_result.Result == Wirehair_Success &&
              timed_result.BytesVerified &&
              timed_result.ElapsedNanoseconds > 0u &&
              timed_state.Calls == 1u &&
              timed_state.SawCanonicalBase,
        "canonical timed path did not use its preserved zero-seed base");
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
    Check(candidate.BaseKind == NativeWh2BaseKind::ProductionProfile,
        "experimental factory no longer defaults to the production profile");
    NativeArm candidate_arm;
    Check(candidate_arm.Initialize(
              candidate, K, block_bytes, source) == Wirehair_Success,
        "runtime candidate exact construction failed");
    Check(candidate_arm.Kind() == NativeArmKind::Wirehair2Experiment &&
              candidate_arm.ConstructionAttempt() == 1u &&
              transform_state.Calls == 1u,
        "runtime candidate transform/identity receipt is wrong");

    TransformState explicit_production_state;
    const NativeArmSpec explicit_production =
        wirehair_wh2_bench::MakeExperimentalWh2Arm(
            1u,
            RuntimeCandidateTransform,
            &explicit_production_state,
            nullptr,
            NativeWh2BaseKind::ProductionProfile);
    NativeArm explicit_production_arm;
    Check(explicit_production_arm.Initialize(
              explicit_production, K, block_bytes, source) == Wirehair_Success &&
              explicit_production_state.Calls == 1u,
        "explicit production-profile experiment changed behavior");
    std::vector<uint8_t> default_packet;
    std::vector<uint8_t> explicit_packet;
    Check(candidate_arm.Encode(K + 9u, default_packet) == Wirehair_Success &&
              explicit_production_arm.Encode(K + 9u, explicit_packet) ==
                  Wirehair_Success &&
              default_packet == explicit_packet,
        "default and explicit production experiments diverged");

    NativeRecoveryFixture fixture;
    Check(fixture.Initialize(candidate_arm, ids) == Wirehair_Success,
        "runtime candidate recovery fixture failed");
    const RecoveryCellResult candidate_recovery = fixture.RunNested();
    Check(candidate_recovery.Outcome == RecoveryOutcome::Success &&
              candidate_recovery.FirstOverhead == 0u,
        "runtime candidate failed byte-verified nested recovery");

    NativeArm two_anchor_arm;
    Check(two_anchor_arm.Initialize(
              wirehair_wh2_bench::MakeExperimentalWh2Arm(
                  0u, TwoAnchorTransform),
              K, block_bytes, source) == Wirehair_Success,
        "enabled two-anchor codec failed exact construction");
    NativeRecoveryFixture two_anchor_fixture;
    Check(two_anchor_fixture.Initialize(two_anchor_arm, ids) ==
              Wirehair_Success,
        "enabled two-anchor recovery fixture failed");
    const RecoveryCellResult two_anchor_recovery =
        two_anchor_fixture.RunNested();
    Check(two_anchor_recovery.Outcome == RecoveryOutcome::Success &&
              two_anchor_recovery.FirstOverhead == 0u,
        "enabled two-anchor codec failed byte-verified recovery");

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
    Check(invalid_candidate.Initialize(
              wirehair_wh2_bench::MakeExperimentalWh2Arm(
                  0u, InvalidDenseAnchorTransform),
              K, block_bytes, source) == Wirehair_InvalidInput,
        "invalid dense-anchor layout was not rejected by native resolution");
    const InvalidDenseAnchorCombination invalid_combinations[] = {
        InvalidDenseAnchorCombination::DenseRows,
        InvalidDenseAnchorCombination::HeavyRows,
        InvalidDenseAnchorCombination::IdentityCorner,
        InvalidDenseAnchorCombination::HeavyFamily
    };
    for (InvalidDenseAnchorCombination combination : invalid_combinations)
    {
        NativeArm invalid_combination;
        Check(invalid_combination.Initialize(
                  wirehair_wh2_bench::MakeExperimentalWh2Arm(
                      0u, InvalidDenseAnchorCombinationTransform,
                      &combination),
                  K, block_bytes, source) == Wirehair_InvalidInput,
            "invalid dense-anchor parameter combination was accepted");
    }
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

    TransformState canonical_rejection_state;
    const NativeArmSpec canonical =
        wirehair_wh2_bench::MakeExperimentalWh2Arm(
            0u,
            RuntimeCandidateTransform,
            &canonical_rejection_state,
            nullptr,
            NativeWh2BaseKind::CanonicalCertifiedStructure);
    std::vector<NativeArmSpec> invalid_canonical_specs;
    NativeArmSpec bad_base_legacy = wirehair_wh2_bench::MakeWirehair1Arm();
    bad_base_legacy.BaseKind =
        NativeWh2BaseKind::CanonicalCertifiedStructure;
    invalid_canonical_specs.push_back(bad_base_legacy);
    NativeArmSpec bad_base_control =
        wirehair_wh2_bench::MakeCertifiedWh2Arm(0u);
    bad_base_control.BaseKind =
        NativeWh2BaseKind::CanonicalCertifiedStructure;
    invalid_canonical_specs.push_back(bad_base_control);
    invalid_canonical_specs.push_back(
        wirehair_wh2_bench::MakeExperimentalWh2Arm(
            0u,
            nullptr,
            nullptr,
            nullptr,
            NativeWh2BaseKind::CanonicalCertifiedStructure));
    NativeArmSpec unknown_base = canonical;
    unknown_base.BaseKind = static_cast<NativeWh2BaseKind>(UINT32_MAX);
    invalid_canonical_specs.push_back(unknown_base);
    for (uint32_t option = 0u; option < 6u; ++option)
    {
        NativeArmSpec nondefault = canonical;
        if (option == 0u) nondefault.Wh2Options.RecoveryMixCount = 2u;
        else if (option == 1u) nondefault.Wh2Options.DenseIdentityCorner = true;
        else if (option == 2u) nondefault.Wh2Options.PrecodeSeedSalt = 1u;
        else if (option == 3u) nondefault.Wh2Options.RecoveryRowSeedSalt = 1u;
        else if (option == 4u)
            nondefault.Wh2Options.CacheSystematicSource = true;
        else
            nondefault.Wh2Options.CacheReceivedSystematicPackets = true;
        invalid_canonical_specs.push_back(nondefault);
    }
    for (const NativeArmSpec& invalid : invalid_canonical_specs)
    {
        Check(arm.Initialize(invalid, K, block_bytes, source) ==
                  Wirehair_InvalidInput,
            "invalid canonical-base combination was accepted");
        after.clear();
        Check(arm.IsInitialized() &&
                  arm.Kind() == NativeArmKind::Wirehair2Certified &&
                  arm.ConstructionAttempt() == 0u &&
                  arm.Encode(K + 5u, after) == Wirehair_Success &&
                  after == before,
            "invalid canonical-base reinitialization modified the prior arm");
    }
    Check(canonical_rejection_state.Calls == 0u,
        "invalid canonical-base spec invoked its transform");

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
    const std::vector<uint32_t> ids = ConsecutiveIds(K, 0u, 64u);

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
        const uint64_t expected_arena_bytes =
            ((uint64_t)result.Stats.PeeledColumns +
                result.Stats.InactivatedColumns) * block_bytes;
        Check(result.Result == Wirehair_Success &&
                  result.BytesVerified &&
                  result.ElapsedNanoseconds > 0u &&
                  result.Stats.PacketRows == K + 4u &&
                  result.Stats.SolveValueArenaBytes ==
                      expected_arena_bytes &&
                  result.Stats.SolveValueArenaEagerZeroBytes == 0u &&
                  result.Stats.SolveValueArenaCommitCopyBytes == 0u,
            "isolated solve did not use the owned no-init arena exactly");
    }

    // A solve fixture consumes only the frozen K+4 prefix of the shared K+64
    // timing trace.  Deliberately corrupting the unused tail must not change
    // the isolated-solve preparation or row count.
    std::vector<uint32_t> ignored_tail = ids;
    std::fill(
        ignored_tail.begin() + K + 4u,
        ignored_tail.end(),
        ignored_tail.front());
    NativeSolveFixture prefix_only;
    Check(prefix_only.Initialize(arm, ignored_tail, 4u) == Wirehair_Success,
        "isolated solve inspected the receive-only trace tail");
    const IsolatedSolveResult prefix_only_result = prefix_only.Run();
    Check(prefix_only_result.Result == Wirehair_Success &&
              prefix_only_result.BytesVerified &&
              prefix_only_result.Stats.PacketRows == K + 4u,
        "isolated solve did not retain exactly the K+4 trace prefix");

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
    Check(rejected.Initialize(arm, ids, 3u) == Wirehair_InvalidInput,
        "isolated solve accepted a prefix other than exactly K+4");
}

void CheckOtherTimingScopes()
{
    const uint32_t K = 36u;
    const uint32_t block_bytes = 64u;
    std::vector<uint8_t> source;
    Check(wirehair_wh2_bench::MakeDeterministicSource(
              K, block_bytes, UINT64_C(0x943a0f6de5b728c1), source),
        "other timing-scope source generation failed");
    const std::vector<uint32_t> ids = ConsecutiveIds(K, 0u, 64u);
    const std::vector<uint32_t> short_ids = ConsecutiveIds(K);
    std::vector<uint32_t> duplicate_ids = ids;
    duplicate_ids.back() = duplicate_ids.front();

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
        Check(receive_fixture.Initialize(arm, ids, 64u) == Wirehair_Success,
            "receive timing fixture preparation failed");
        Check(receive_fixture.Initialize(arm, short_ids, 64u) ==
                  Wirehair_InvalidInput &&
                  receive_fixture.IsInitialized(),
            "receive timing accepted a short trace or lost prior state");
        Check(receive_fixture.Initialize(arm, ids, 63u) ==
                  Wirehair_InvalidInput &&
                  receive_fixture.IsInitialized(),
            "receive timing accepted a trace longer than its exact cap");
        Check(receive_fixture.Initialize(arm, duplicate_ids, 64u) ==
                  Wirehair_InvalidInput &&
                  receive_fixture.IsInitialized(),
            "receive timing accepted duplicate IDs in its extended trace");
        Check(receive_fixture.Initialize(arm, ids, UINT32_MAX) ==
                  Wirehair_InvalidInput &&
                  receive_fixture.IsInitialized(),
            "receive timing accepted an overflowing K+cap domain");
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
        wirehair_v2::ResetDecoderReceivePathCountersForTesting();
#endif
        const wirehair_wh2_bench::TimedArmResult receive_result =
            receive_fixture.Run();
        Check(receive_result.Result == Wirehair_Success &&
                  receive_result.BytesVerified &&
                  receive_result.ElapsedNanoseconds > 0u &&
                  receive_result.DecodedOverhead == 0u,
            "fresh receive-to-success timing scope failed");
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
        const wirehair_v2::DecoderReceivePathCounters counters =
            wirehair_v2::DecoderReceivePathCountersForTesting();
        if (spec.Kind == NativeArmKind::Wirehair1)
        {
            Check(counters.ValidatedSystemInitializations == 0u &&
                      counters.ColdSolveAttempts == 0u &&
                      counters.RecoveryPreflights == 0u,
                "Wirehair1 timing unexpectedly entered the WH2 decoder");
        }
        else
        {
            Check(counters.ValidatedSystemInitializations == 1u &&
                      counters.LastValidatedPacketSeedAttempt ==
                          spec.ConstructionAttempt,
                "WH2 timing did not initialize its exact native attempt");
            Check(counters.DirectSystematicCompletions == 1u &&
                      counters.DirectSystematicCanonicalizations == 1u &&
                      counters.DirectSystematicMaterializationAttempts == 0u &&
                      counters.ColdSolveAttempts == 0u &&
                      counters.ColdSolvePacketAssemblies == 0u &&
                      counters.ColdSolveSlotEntries == 0u &&
                      counters.ColdSolvePayloadBytes == 0u,
                "WH2 systematic timing did not use direct completion");
            Check(counters.CheckpointAdoptions == 0u &&
                      counters.PendingPacketCopies == 0u &&
                      counters.ResumeAttempts == 0u,
                "full-rank systematic timing unexpectedly used a checkpoint");
            Check(counters.RecoveryPreflights == 1u &&
                      counters.RecoveryPacketEvaluations == 0u &&
                      counters.RecoveryColdPacketCopies == K &&
                      counters.RecoveryColdPacketCopyBytes ==
                          (uint64_t)K * block_bytes &&
                      counters.RecoveryColdStorageReleases == 0u,
                "WH2 timing bypassed repeatable direct recovery reuse");

            if (spec.Kind == NativeArmKind::Wirehair2Experiment)
            {
                // A systematic-only trace cannot distinguish packet-row
                // configurations.  Require the transformed, nonzero-attempt
                // arm to recover from its own repair-only symbols as well.
                const std::vector<uint32_t> repair_ids =
                    ConsecutiveIds(K, K, 64u);
                NativeReceiveFixture repair_receive;
                Check(repair_receive.Initialize(
                          arm, repair_ids, 64u) == Wirehair_Success,
                    "experimental repair-only receive setup failed");
                wirehair_v2::ResetDecoderReceivePathCountersForTesting();
                const wirehair_wh2_bench::TimedArmResult repair_result =
                    repair_receive.Run();
                const wirehair_v2::DecoderReceivePathCounters
                    repair_counters =
                        wirehair_v2::
                            DecoderReceivePathCountersForTesting();
                Check(repair_result.Result == Wirehair_Success &&
                          repair_result.BytesVerified &&
                          repair_result.DecodedOverhead <= 64u,
                    "experimental exact configuration failed repair-only receive");
                Check(repair_counters.ValidatedSystemInitializations == 1u &&
                          repair_counters.LastValidatedPacketSeedAttempt ==
                              spec.ConstructionAttempt &&
                          repair_counters.RecoveryPreflights == 1u,
                    "experimental repair-only receive lost exact attempt state");
            }
        }
#endif
    }
}

void CheckNativeReceiveTrustedResume()
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    const uint32_t K = 1000u;
    const uint32_t block_bytes = 257u;
    const NativeArmSpec spec =
        wirehair_wh2_bench::MakeCertifiedWh2Arm(0u);
    wirehair_wh2_bench::ResolvedNativeWh2Configuration resolved;
    wirehair_v2::PrecodeSystem system;
    uint32_t collision_a = 0u;
    uint32_t collision_b = 0u;
    Check(wirehair_wh2_bench::ResolveNativeWh2Configuration(
              spec, K, block_bytes, resolved) &&
              wirehair_v2::BuildPrecodeSystem(resolved.Params, system) &&
              FindFirstPacketRowCollision(
                  system,
                  resolved.PacketConfig,
                  collision_a,
                  collision_b),
        "trusted receive collision fixture construction failed");
    if (collision_a == collision_b) {
        return;
    }

    std::vector<uint8_t> source;
    Check(wirehair_wh2_bench::MakeDeterministicSource(
              K,
              block_bytes,
              UINT64_C(0xe9d348672bc150af),
              source),
        "trusted receive source generation failed");
    NativeArm arm;
    Check(arm.Initialize(spec, K, block_bytes, source) == Wirehair_Success,
        "trusted receive arm construction failed");

    std::vector<uint32_t> ids;
    ids.reserve((size_t)K + 2u);
    for (uint32_t id = 0u; id + 2u < K; ++id) {
        ids.push_back(id);
    }
    ids.push_back(collision_a);
    ids.push_back(collision_b);
    ids.push_back(K - 2u);
    ids.push_back(K - 1u);

    NativeReceiveFixture receive;
    Check(receive.Initialize(arm, ids, 2u) == Wirehair_Success,
        "trusted receive fixture initialization failed");
    wirehair_v2::ResetResumeSystemFingerprintChecksForTesting();
    wirehair_v2::ResetPackedBinaryResidualUsesForTesting();
    wirehair_v2::ResetDecoderReceivePathCountersForTesting();
    const wirehair_wh2_bench::TimedArmResult result = receive.Run();
    Check(result.Result == Wirehair_Success &&
              result.BytesVerified &&
              result.DecodedOverhead >= 1u &&
              result.DecodedOverhead <= 2u,
        "trusted receive fixture did not exercise rank-deficient resume");
    Check(wirehair_v2::ResumeSystemFingerprintChecksForTesting() == 0u,
        "native receive recomputed the immutable system fingerprint");
    Check(wirehair_v2::PackedBinaryResidualUsesForTesting() == 0u,
        "trusted receive unexpectedly cold-fell back under its budget");
    const wirehair_v2::DecoderReceivePathCounters counters =
        wirehair_v2::DecoderReceivePathCountersForTesting();
    Check(counters.ValidatedSystemInitializations == 1u &&
              counters.LastValidatedPacketSeedAttempt == 0u &&
              counters.ColdSolveAttempts == 1u &&
              counters.ColdSolvePacketAssemblies == K &&
              counters.ColdSolveSlotEntries == K &&
              counters.ColdSolvePayloadBytes ==
                  (uint64_t)K * block_bytes,
        "trusted receive bypassed exact decoder cold-receive work");
    Check(counters.CheckpointPendingAllocationAttempts == 1u &&
              counters.CheckpointAdoptions == 1u &&
              counters.PendingPacketCopies == result.DecodedOverhead &&
              counters.PendingPacketCopyBytes ==
                  (uint64_t)result.DecodedOverhead * block_bytes &&
              counters.ResumeAttempts == result.DecodedOverhead,
        "trusted receive bypassed checkpoint pending-copy/resume work");
    Check(counters.RecoveryPreflights == 1u &&
              counters.RecoveryPacketEvaluations == K &&
              counters.RecoveryColdPacketCopies == 0u &&
              counters.RecoveryColdStorageReleases == 0u,
        "trusted receive bypassed decoder recovery preflight");

    PublicReceiveResult public_result;
    Check(RunPublicWh2Receive(
              arm, source, 0u, ids, 2u, public_result),
        "trusted receive public differential setup failed");
    Check(public_result.Result == WirehairV2_Success &&
              public_result.BytesVerified &&
              public_result.DecodedOverhead == result.DecodedOverhead,
        "native exact-system receive differs from the public decoder");
#endif
}

void CheckNativeReceivePackedDispatchSeam()
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    const uint32_t K = 1000u;
    const uint32_t block_sizes[] = {
        64u,
        wirehair_v2::kDecoderPackedBinaryResidualMinBlockBytes - 1u,
        wirehair_v2::kDecoderPackedBinaryResidualMinBlockBytes
    };
    const NativeArmSpec spec =
        wirehair_wh2_bench::MakeCertifiedWh2Arm(0u);
    for (uint32_t block_bytes : block_sizes)
    {
        std::vector<uint8_t> source;
        if (!Check(wirehair_wh2_bench::MakeDeterministicSource(
                      K,
                      block_bytes,
                      UINT64_C(0x7f49d8e12ab630c5) ^ block_bytes,
                      source),
                "native receive seam source generation failed"))
        {
            return;
        }
        NativeArm arm;
        if (!Check(
                arm.Initialize(spec, K, block_bytes, source) ==
                    Wirehair_Success,
                "native receive seam arm initialization failed"))
        {
            return;
        }
        NativeReceiveFixture receive;
        std::vector<uint32_t> ids;
        ids.reserve((size_t)K + 64u);
        for (uint32_t id = 0u; id + 1u < K; ++id) {
            ids.push_back(id);
        }
        for (uint32_t repair = 0u; repair < 65u; ++repair) {
            ids.push_back(K + repair);
        }
        if (!Check(
                receive.Initialize(arm, ids, 64u) == Wirehair_Success,
                "native receive seam fixture initialization failed"))
        {
            return;
        }
        wirehair_v2::ResetPackedBinaryResidualUsesForTesting();
        const wirehair_wh2_bench::TimedArmResult result = receive.Run();
        const uint64_t expected_uses =
            block_bytes == 64u || block_bytes >=
                    wirehair_v2::
                        kDecoderPackedBinaryResidualMinBlockBytes ?
                1u : 0u;
        Check(result.Result == Wirehair_Success && result.BytesVerified,
            "native receive seam did not decode exact bytes");
        Check(wirehair_v2::PackedBinaryResidualUsesForTesting() ==
                  expected_uses,
            "native receive seam selected wrong packed policy");
    }
#endif
}

void CheckNativeReceiveBudgetColdFallback()
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    const uint32_t K = 1000u;
    const uint32_t block_bytes = 64u;
    const NativeArmSpec spec =
        wirehair_wh2_bench::MakeCertifiedWh2Arm(0u);
    wirehair_wh2_bench::ResolvedNativeWh2Configuration resolved;
    wirehair_v2::PrecodeSystem system;
    uint32_t collision_a = 0u;
    uint32_t collision_b = 0u;
    if (!Check(wirehair_wh2_bench::ResolveNativeWh2Configuration(
                  spec, K, block_bytes, resolved) &&
                  wirehair_v2::BuildPrecodeSystem(
                      resolved.Params, system) &&
                  FindFirstPacketRowCollision(
                      system,
                      resolved.PacketConfig,
                      collision_a,
                      collision_b),
            "native budget fallback collision fixture failed"))
    {
        return;
    }
    std::vector<uint8_t> source;
    if (!Check(wirehair_wh2_bench::MakeDeterministicSource(
                  K,
                  block_bytes,
                  UINT64_C(0x3d98b12fa7406ce5),
                  source),
            "native budget fallback source generation failed"))
    {
        return;
    }
    NativeArm arm;
    if (!Check(
            arm.Initialize(spec, K, block_bytes, source) == Wirehair_Success,
            "native budget fallback arm initialization failed"))
    {
        return;
    }
    std::vector<uint32_t> ids;
    ids.reserve((size_t)K + 2u);
    for (uint32_t id = 0u; id + 2u < K; ++id) {
        ids.push_back(id);
    }
    ids.push_back(collision_a);
    ids.push_back(collision_b);
    ids.push_back(K - 2u);
    ids.push_back(K - 1u);
    NativeReceiveFixture receive;
    if (!Check(
            receive.Initialize(arm, ids, 2u) == Wirehair_Success,
            "native budget fallback receive initialization failed"))
    {
        return;
    }
    wirehair_v2::ResetPackedBinaryResidualUsesForTesting();
    const wirehair_wh2_bench::TimedArmResult result = receive.Run();
    if (!Check(
            result.Result == Wirehair_Success && result.BytesVerified &&
                result.DecodedOverhead >= 1u &&
                result.DecodedOverhead <= 2u,
            "native budget fallback did not recover after deficient K") ||
        !Check(
            wirehair_v2::PackedBinaryResidualUsesForTesting() ==
                (uint64_t)result.DecodedOverhead + 1u,
            "native budget fallback retained a rejected checkpoint"))
    {
        return;
    }
#endif
}

void CheckValidatedDecoderInitializationTransactionality()
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    const uint32_t K = 36u;
    const uint32_t block_bytes = 64u;
    const NativeArmSpec spec =
        wirehair_wh2_bench::MakeCertifiedWh2Arm(0u);
    wirehair_wh2_bench::ResolvedNativeWh2Configuration resolved;
    wirehair_v2::PrecodeSystem system;
    std::vector<uint8_t> source;
    NativeArm arm;
    if (!Check(wirehair_wh2_bench::ResolveNativeWh2Configuration(
                  spec, K, block_bytes, resolved) &&
                  wirehair_v2::BuildPrecodeSystem(
                      resolved.Params, system) &&
                  wirehair_wh2_bench::MakeDeterministicSource(
                      K,
                      block_bytes,
                      UINT64_C(0xd63b8a79f142c05e),
                      source) &&
                  arm.Initialize(
                      spec, K, block_bytes, source) == Wirehair_Success,
            "validated decoder transaction fixture setup failed"))
    {
        return;
    }

    wirehair_v2::MessagePrecodeDecoder decoder;
    if (!Check(decoder.InitializeForValidatedSystemForTesting(
                  source.size(),
                  block_bytes,
                  system,
                  resolved.PacketConfig,
                  resolved.PacketAttempt) == Wirehair_Success,
            "validated decoder exact initialization failed"))
    {
        return;
    }
    const uint64_t preserved_seed = decoder.System().Params.Seed;
    const uint32_t preserved_attempt = decoder.PacketSeedAttempt();

    wirehair_v2::PrecodeSystem malformed(system);
    malformed.StaircaseRows.pop_back();
    Check(decoder.InitializeForValidatedSystemForTesting(
              source.size(),
              block_bytes,
              malformed,
              resolved.PacketConfig,
              resolved.PacketAttempt) == Wirehair_InvalidInput &&
              decoder.IsInitialized() && !decoder.IsDecoded() &&
              decoder.System().Params.Seed == preserved_seed &&
              decoder.PacketSeedAttempt() == preserved_attempt,
        "invalid exact-system initialization modified prior decoder state");

    wirehair_v2::PacketRowConfig invalid_config = resolved.PacketConfig;
    invalid_config.MixCount = UINT32_MAX;
    Check(decoder.InitializeForValidatedSystemForTesting(
              source.size(),
              block_bytes,
              system,
              invalid_config,
              resolved.PacketAttempt) == Wirehair_InvalidInput &&
              decoder.IsInitialized() && !decoder.IsDecoded() &&
              decoder.System().Params.Seed == preserved_seed &&
              decoder.PacketSeedAttempt() == preserved_attempt,
        "invalid exact packet configuration modified prior decoder state");

    for (int64_t failure = 0; failure < 4; ++failure)
    {
        wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(failure);
        const WirehairResult allocation_result =
            decoder.InitializeForValidatedSystemForTesting(
                source.size(),
                block_bytes,
                system,
                resolved.PacketConfig,
                resolved.PacketAttempt);
        wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(-1);
        Check(allocation_result == Wirehair_OOM &&
                  decoder.IsInitialized() && !decoder.IsDecoded() &&
                  decoder.System().Params.Seed == preserved_seed &&
                  decoder.PacketSeedAttempt() == preserved_attempt,
            "exact-system initialization OOM modified prior decoder state");
    }

    std::vector<uint8_t> packet(block_bytes, 0u);
    WirehairResult decode_result = Wirehair_NeedMore;
    for (uint32_t block_id = 0u; block_id < K; ++block_id)
    {
        if (!Check(arm.EncodeInto(
                      block_id, packet.data(), block_bytes) ==
                      Wirehair_Success,
                "preserved decoder packet encode failed"))
        {
            return;
        }
        decode_result = decoder.DecodeResult(
            block_id, packet.data(), block_bytes);
    }
    std::vector<uint8_t> recovered(source.size(), 0u);
    Check(decode_result == Wirehair_Success &&
              decoder.RecoverResult(
                  recovered.data(), recovered.size()) == Wirehair_Success &&
              recovered == source,
        "failed exact-system reinitialization corrupted prior decoder");
#endif
}

void CheckNativeReceiveAllocationCoverage()
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    const uint32_t K = 1000u;
    const uint32_t block_bytes = 257u;
    const NativeArmSpec spec =
        wirehair_wh2_bench::MakeCertifiedWh2Arm(0u);
    wirehair_wh2_bench::ResolvedNativeWh2Configuration resolved;
    wirehair_v2::PrecodeSystem system;
    uint32_t collision_a = 0u;
    uint32_t collision_b = 0u;
    std::vector<uint8_t> source;
    NativeArm arm;
    if (!Check(wirehair_wh2_bench::ResolveNativeWh2Configuration(
                  spec, K, block_bytes, resolved) &&
                  wirehair_v2::BuildPrecodeSystem(
                      resolved.Params, system) &&
                  FindFirstPacketRowCollision(
                      system,
                      resolved.PacketConfig,
                      collision_a,
                      collision_b) &&
                  wirehair_wh2_bench::MakeDeterministicSource(
                      K,
                      block_bytes,
                      UINT64_C(0xf067e3a29bc45d18),
                      source) &&
                  arm.Initialize(
                      spec, K, block_bytes, source) == Wirehair_Success,
            "native receive allocation fixture setup failed"))
    {
        return;
    }

    std::vector<uint32_t> ids;
    ids.reserve((size_t)K + 2u);
    for (uint32_t id = 0u; id + 2u < K; ++id) {
        ids.push_back(id);
    }
    ids.push_back(collision_a);
    ids.push_back(collision_b);
    ids.push_back(K - 2u);
    ids.push_back(K - 1u);
    NativeReceiveFixture receive;
    if (!Check(receive.Initialize(arm, ids, 2u) == Wirehair_Success,
            "native receive allocation fixture initialization failed"))
    {
        return;
    }

    wirehair_v2::ResetDecoderReceivePathCountersForTesting();
    wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(0);
    const wirehair_wh2_bench::TimedArmResult setup_oom = receive.Run();
    wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(-1);
    const wirehair_v2::DecoderReceivePathCounters setup_counters =
        wirehair_v2::DecoderReceivePathCountersForTesting();
    Check(setup_oom.Result == Wirehair_OOM &&
              setup_oom.ElapsedNanoseconds == 0u &&
              setup_counters.ValidatedSystemInitializations == 0u &&
              setup_counters.ColdSolveAttempts == 0u,
        "native receive did not expose exact-decoder setup allocation failure");

    // Four guarded allocations initialize the exact decoder (owned system,
    // id reserve, payload reserve, slot table).  The fifth is the first cold
    // solve boundary and must occur after K timed payload/slot insertions and
    // SolvePacket assembly.
    wirehair_v2::ResetDecoderReceivePathCountersForTesting();
    wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(4);
    const wirehair_wh2_bench::TimedArmResult solve_oom = receive.Run();
    wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(-1);
    const wirehair_v2::DecoderReceivePathCounters solve_counters =
        wirehair_v2::DecoderReceivePathCountersForTesting();
    Check(solve_oom.Result == Wirehair_OOM &&
              solve_oom.ElapsedNanoseconds == 0u &&
              solve_counters.ValidatedSystemInitializations == 1u &&
              solve_counters.ColdSolveAttempts == 1u &&
              solve_counters.ColdSolvePacketAssemblies == K &&
              solve_counters.ColdSolveSlotEntries == K &&
              solve_counters.ColdSolvePayloadBytes ==
                  (uint64_t)K * block_bytes &&
              solve_counters.RecoveryPreflights == 0u,
        "native receive solve OOM bypassed timed receive bookkeeping");

    // Let the deficient solve finish, then fail the production pending-slot
    // allocation immediately before checkpoint adoption.
    wirehair_v2::ResetDecoderReceivePathCountersForTesting();
    wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(5);
    const wirehair_wh2_bench::TimedArmResult pending_oom = receive.Run();
    wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(-1);
    const wirehair_v2::DecoderReceivePathCounters pending_counters =
        wirehair_v2::DecoderReceivePathCountersForTesting();
    Check(pending_oom.Result == Wirehair_OOM &&
              pending_oom.ElapsedNanoseconds == 0u &&
              pending_counters.ValidatedSystemInitializations == 1u &&
              pending_counters.ColdSolveAttempts == 1u &&
              pending_counters.CheckpointPendingAllocationAttempts == 1u &&
              pending_counters.CheckpointAdoptions == 0u &&
              pending_counters.PendingPacketCopies == 0u,
        "native receive bypassed checkpoint pending-slot allocation");

    wirehair_v2::ResetDecoderReceivePathCountersForTesting();
    const wirehair_wh2_bench::TimedArmResult retry = receive.Run();
    const wirehair_v2::DecoderReceivePathCounters retry_counters =
        wirehair_v2::DecoderReceivePathCountersForTesting();
    Check(retry.Result == Wirehair_Success && retry.BytesVerified &&
              retry.DecodedOverhead >= 1u &&
              retry.DecodedOverhead <= 2u &&
              retry_counters.CheckpointAdoptions == 1u &&
              retry_counters.PendingPacketCopies == retry.DecodedOverhead &&
              retry_counters.RecoveryPreflights == 1u,
        "native receive fixture was not reusable after allocation failures");
#endif
}

void CheckTimingControlProbe()
{
    struct Shape
    {
        uint32_t K;
        uint32_t BlockBytes;
        uint64_t SourceSeed;
    };
    const Shape shapes[] = {
        { 8u, 64u, UINT64_C(0x3c79ac492ba7d113) },
        { 32u, 1280u, UINT64_C(0x79d8e31b5a40c6f2) },
        { 101u, 4096u, UINT64_C(0xa44e91b0367fc825) }
    };

    NativeTimingControlProbe uninitialized;
    const NativeTimingControlQualificationResult uninitialized_result =
        uninitialized.Run(ConsecutiveIds(8u, 1u, 256u));
    Check(uninitialized_result.Qualification ==
              NativeTimingControlQualification::Fatal &&
              uninitialized_result.Wirehair1Result ==
                  Wirehair_InvalidInput &&
              uninitialized_result.Wirehair2HeadResult ==
                  Wirehair_InvalidInput,
        "uninitialized timing-control probe did not fail closed");

    for (const Shape& shape : shapes)
    {
        const NativeArmSpec head_spec =
            wirehair_wh2_bench::MakeCertifiedWh2Arm(0u);
        wirehair_wh2_bench::ResolvedNativeWh2Configuration actual_width;
        wirehair_wh2_bench::ResolvedNativeWh2Configuration unit_width;
        Check(wirehair_wh2_bench::ResolveNativeWh2Configuration(
                  head_spec, shape.K, shape.BlockBytes, actual_width) &&
              wirehair_wh2_bench::ResolveNativeWh2Configuration(
                  head_spec, shape.K, 1u, unit_width) &&
              (actual_width.Params.Seed != unit_width.Params.Seed ||
               actual_width.PacketConfig.PeelSeed !=
                   unit_width.PacketConfig.PeelSeed),
            "representative WH2 width did not affect resolved equations");

        std::vector<uint8_t> source;
        Check(wirehair_wh2_bench::MakeDeterministicSource(
                  shape.K, shape.BlockBytes, shape.SourceSeed, source),
            "timing-control differential source generation failed");
        const std::vector<uint32_t> ids =
            ConsecutiveIds(shape.K, 1u, 256u);

        NativeTimingControlProbe probe;
        Check(probe.Initialize(
                  head_spec, shape.K, shape.BlockBytes) == Wirehair_Success &&
              probe.IsInitialized(),
            "timing-control probe initialization failed");
        const NativeTimingControlQualificationResult structural =
            probe.Run(ids);

        NativeArm wh1;
        Check(wh1.Initialize(
                  wirehair_wh2_bench::MakeWirehair1Arm(),
                  shape.K, shape.BlockBytes, source) == Wirehair_Success,
            "full-data Wirehair1 differential construction failed");
        NativeReceiveFixture receive;
        Check(receive.Initialize(wh1, ids, 256u) == Wirehair_Success,
            "full-data Wirehair1 differential fixture failed");
        const wirehair_wh2_bench::TimedArmResult full_receive = receive.Run();

        NativeArm wh2;
        Check(wh2.Initialize(
                  head_spec, shape.K, shape.BlockBytes, source) ==
                  Wirehair_Success,
            "full-data WH2-head differential construction failed");
        NativeSolveFixture solve;
        Check(solve.Initialize(wh2, ids, 4u) == Wirehair_Success,
            "full-data WH2-head differential fixture failed");
        const IsolatedSolveResult full_solve = solve.Run();

        Check(structural.Wirehair1Result == full_receive.Result &&
                  structural.Wirehair2HeadResult == full_solve.Result,
            "minimal-payload timing-control outcome differs from full data");
        if (full_receive.Result == Wirehair_Success)
        {
            Check(full_receive.BytesVerified &&
                      structural.Wirehair1DecodedOverhead ==
                          full_receive.DecodedOverhead,
                "minimal-payload Wirehair1 overhead differs from full data");
        }
        if (full_solve.Result == Wirehair_Success) {
            Check(full_solve.BytesVerified,
                "full-data WH2-head differential did not recover bytes");
        }

        const bool both_success =
            full_receive.Result == Wirehair_Success &&
            full_solve.Result == Wirehair_Success;
        const bool both_retryable =
            (full_receive.Result == Wirehair_Success ||
             full_receive.Result == Wirehair_NeedMore) &&
            (full_solve.Result == Wirehair_Success ||
             full_solve.Result == Wirehair_NeedMore);
        const NativeTimingControlQualification expected = both_success ?
            NativeTimingControlQualification::Success :
            (both_retryable ? NativeTimingControlQualification::NeedMore :
                NativeTimingControlQualification::Fatal);
        Check(structural.Qualification == expected,
            "timing-control qualification misclassified exact outcomes");

        std::vector<uint32_t> short_ids(ids.begin(), ids.end() - 1u);
        const NativeTimingControlQualificationResult short_result =
            probe.Run(short_ids);
        Check(short_result.Qualification ==
                  NativeTimingControlQualification::Fatal &&
                  short_result.Wirehair1Result == Wirehair_InvalidInput &&
                  short_result.Wirehair2HeadResult == Wirehair_InvalidInput &&
                  probe.IsInitialized(),
            "timing-control probe accepted a short K+256 trace");
        std::vector<uint32_t> duplicate_ids(ids);
        duplicate_ids.back() = duplicate_ids.front();
        const NativeTimingControlQualificationResult duplicate_result =
            probe.Run(duplicate_ids);
        Check(duplicate_result.Qualification ==
                  NativeTimingControlQualification::Fatal &&
                  duplicate_result.Wirehair1Result == Wirehair_InvalidInput &&
                  duplicate_result.Wirehair2HeadResult ==
                      Wirehair_InvalidInput &&
                  probe.IsInitialized(),
            "timing-control probe accepted duplicate IDs");

        Check(probe.Initialize(
                  wirehair_wh2_bench::MakeWirehair1Arm(),
                  shape.K, shape.BlockBytes) == Wirehair_InvalidInput &&
                  probe.IsInitialized(),
            "timing-control probe accepted WH1 as its head or lost state");
        const NativeTimingControlQualificationResult repeated = probe.Run(ids);
        Check(repeated.Qualification == structural.Qualification &&
                  repeated.Wirehair1Result == structural.Wirehair1Result &&
                  repeated.Wirehair1DecodedOverhead ==
                      structural.Wirehair1DecodedOverhead &&
                  repeated.Wirehair2HeadResult ==
                      structural.Wirehair2HeadResult,
            "reused timing-control probe changed deterministic outcome");
    }

    // This is the frozen K=64000, width=64, replicate=84 loss stream.  It is
    // the known LegacyCurrent identity that remains NeedMore through +256;
    // the WH2-head K+4 system succeeds.  Pin the retryable classification and
    // compare both minimal-payload outcomes with the real-width full-data
    // fixtures so a future probe cannot silently treat rank as byte recovery.
    const uint32_t K = 64000u;
    const uint32_t block_bytes = 64u;
    const IidTrace retry_trace = FrozenIidIds(
        K,
        block_bytes,
        UINT64_C(0x1b52cf340f1a5879),
        256u);
    Check(retry_trace.Ids.size() == static_cast<size_t>(K) + 256u &&
              retry_trace.AttemptedCandidates == UINT64_C(71472),
        "known retryable timing trace changed");
    const NativeArmSpec head_spec =
        wirehair_wh2_bench::MakeCertifiedWh2Arm(0u);
    NativeTimingControlProbe retry_probe;
    Check(retry_probe.Initialize(head_spec, K, block_bytes) ==
              Wirehair_Success,
        "known retryable timing-control probe initialization failed");
    const NativeTimingControlQualificationResult retry_result =
        retry_probe.Run(retry_trace.Ids);
    Check(retry_result.Qualification ==
              NativeTimingControlQualification::NeedMore &&
              retry_result.Wirehair1Result == Wirehair_NeedMore &&
              retry_result.Wirehair1DecodedOverhead == UINT32_MAX &&
              retry_result.Wirehair2HeadResult == Wirehair_Success,
        "known cap-256 timing trace was not classified retryable");

    std::vector<uint8_t> source;
    Check(wirehair_wh2_bench::MakeDeterministicSource(
              K, block_bytes, UINT64_C(0x68e29714c35ba0df), source),
        "known retryable differential source generation failed");
    NativeArm wh1;
    Check(wh1.Initialize(
              wirehair_wh2_bench::MakeWirehair1Arm(),
              K, block_bytes, source) == Wirehair_Success,
        "known retryable full-data WH1 construction failed");
    NativeReceiveFixture receive;
    Check(receive.Initialize(wh1, retry_trace.Ids, 256u) ==
              Wirehair_Success,
        "known retryable full-data WH1 fixture failed");
    const wirehair_wh2_bench::TimedArmResult full_receive = receive.Run();
    Check(full_receive.Result == retry_result.Wirehair1Result &&
              !full_receive.BytesVerified &&
              full_receive.DecodedOverhead == UINT32_MAX,
        "known retryable WH1 minimal/full-data outcomes differ");

    NativeArm wh2;
    Check(wh2.Initialize(head_spec, K, block_bytes, source) ==
              Wirehair_Success,
        "known retryable full-data WH2 construction failed");
    NativeSolveFixture solve;
    Check(solve.Initialize(wh2, retry_trace.Ids, 4u) == Wirehair_Success,
        "known retryable full-data WH2 solve fixture failed");
    const IsolatedSolveResult full_solve = solve.Run();
    Check(full_solve.Result == retry_result.Wirehair2HeadResult &&
              full_solve.BytesVerified,
        "known retryable WH2 minimal/full-data outcomes differ");
}

} // namespace

int main()
{
    Check(wirehair_init() == Wirehair_Success,
        "wirehair initialization failed");
    CheckDeterministicSource();
    CheckCertifiedMatchesPublicProfile();
    CheckSourceIndependentRankClassification();
    CheckCanonicalCertifiedBase();
    CheckArmsAndNestedRecovery();
    CheckTransactionalArmAndValidation();
    CheckIsolatedSolveFixture();
    CheckOtherTimingScopes();
    CheckNativeReceiveTrustedResume();
    CheckNativeReceivePackedDispatchSeam();
    CheckNativeReceiveBudgetColdFallback();
    CheckValidatedDecoderInitializationTransactionality();
    CheckNativeReceiveAllocationCoverage();
    CheckTimingControlProbe();
    if (Failures != 0)
    {
        std::fprintf(stderr, "%d native codec test(s) failed\n", Failures);
        return 1;
    }
    std::puts("WH2 native codec tests passed");
    return 0;
}
