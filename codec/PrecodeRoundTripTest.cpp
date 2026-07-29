#include "WirehairV2Codec.h"
#include "WirehairV2PrecodeDecode.h"

#include "../WirehairTools.h"

#include <algorithm>
#include <atomic>
#include <cstdio>
#include <chrono>
#include <cstring>
#include <limits>
#include <thread>
#include <vector>

namespace {

struct TrialRng
{
    uint64_t State;
    explicit TrialRng(uint64_t seed) : State(seed) {}
    uint64_t Next()
    {
        uint64_t z = (State += UINT64_C(0x9e3779b97f4a7c15));
        z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
        return z ^ (z >> 31);
    }
    bool Drop()
    {
        const double unit = static_cast<double>(Next() >> 11) *
            (1.0 / 9007199254740992.0);
        return unit < 0.10;
    }
};

struct EncodedPacket
{
    uint32_t Id = 0;
    uint32_t Bytes = 0;
    std::vector<uint8_t> Data;
};

bool EncodePacket(
    wirehair_v2::Codec& encoder,
    uint32_t id,
    uint32_t block_bytes,
    EncodedPacket& packet)
{
    packet.Id = id;
    packet.Data.assign(block_bytes, 0xa5u);
    packet.Bytes = UINT32_MAX;
    return encoder.Encode(
        id, packet.Data.data(), block_bytes, &packet.Bytes) ==
        Wirehair_Success;
}

bool RunFacadeLossCase(
    uint32_t K,
    uint32_t block_bytes,
    uint32_t tail_bytes,
    uint32_t loss_stride,
    bool repair_only,
    bool mixed_profile = false)
{
    const uint64_t message_bytes =
        (uint64_t)(K - 1u) * block_bytes + tail_bytes;
    std::vector<uint8_t> message((size_t)message_bytes);
    for (size_t i = 0; i < message.size(); ++i) {
        message[i] = (uint8_t)(i * 67u + K * 13u + 5u);
    }

    wirehair_v2::Codec encoder;
    wirehair_v2::Codec decoder;
    wirehair_v2::MessagePrecodeEncoderOptions options;
    if (mixed_profile) {
        options.Completion =
            wirehair_v2::CompletionField::MixedGF256GF16;
    }
    if (encoder.InitializePrecodeEncoder(
            message.data(), message_bytes, block_bytes,
            nullptr, &options) != Wirehair_Success ||
        decoder.InitializePrecodeDecoder(
            message_bytes, block_bytes, &encoder.Profile()) != Wirehair_Success)
    {
        std::fprintf(stderr,
            "roundtrip: default init failed K=%u bb=%u\n", K, block_bytes);
        return false;
    }

    std::vector<uint8_t> recovered((size_t)message_bytes, 0xccu);
    if (decoder.Recover(recovered.data(), message_bytes) != Wirehair_NeedMore) {
        std::fprintf(stderr, "roundtrip: early recover was not NeedMore\n");
        return false;
    }

    WirehairResult decode_result = Wirehair_NeedMore;
    uint32_t delivered = 0u;
    if (!repair_only)
    {
        for (uint32_t id = K; id-- > 0u;)
        {
            if (id % loss_stride == 0u) {
                continue;
            }
            EncodedPacket packet;
            if (!EncodePacket(encoder, id, block_bytes, packet)) {
                return false;
            }
            decode_result = decoder.Decode(
                packet.Id, packet.Data.data(), packet.Bytes);
            ++delivered;
            if (decode_result != Wirehair_NeedMore &&
                decode_result != Wirehair_Success)
            {
                std::fprintf(stderr,
                    "roundtrip: source feed failed id=%u result=%d\n",
                    id, (int)decode_result);
                return false;
            }
        }
    }

    uint32_t repair_id = K;
    while (decode_result != Wirehair_Success && repair_id < K + K + 64u)
    {
        EncodedPacket packet;
        if (!EncodePacket(encoder, repair_id++, block_bytes, packet)) {
            return false;
        }
        decode_result = decoder.Decode(
            packet.Id, packet.Data.data(), packet.Bytes);
        ++delivered;
        if (decode_result != Wirehair_NeedMore &&
            decode_result != Wirehair_Success)
        {
            std::fprintf(stderr,
                "roundtrip: repair feed failed id=%u result=%d\n",
                packet.Id, (int)decode_result);
            return false;
        }
    }
    if (decode_result != Wirehair_Success)
    {
        std::fprintf(stderr,
            "roundtrip: recover failed K=%u bb=%u delivered=%u result=%d\n",
            K, block_bytes, delivered, (int)decode_result);
        return false;
    }
    const std::vector<uint8_t> before_oom = recovered;
    wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(0);
    const WirehairResult recover_oom =
        decoder.Recover(recovered.data(), message_bytes);
    wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(-1);
    if (tail_bytes < block_bytes &&
        (recover_oom != Wirehair_OOM || recovered != before_oom))
    {
        std::fprintf(stderr,
            "roundtrip: Recover OOM partially modified output\n");
        return false;
    }
    if (tail_bytes == block_bytes &&
        (recover_oom != Wirehair_Success || recovered != message))
    {
        std::fprintf(stderr,
            "roundtrip: exact Recover hit guarded decoder allocation\n");
        return false;
    }
    wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(1);
    const WirehairResult recover_result =
        decoder.Recover(recovered.data(), message_bytes);
    wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(-1);
    if (recover_result != Wirehair_Success || recovered != message)
    {
        std::fprintf(stderr,
            "roundtrip: allocation-free row recovery failed\n");
        return false;
    }

    // Completed decoders remain idempotently successful on duplicates.
    EncodedPacket duplicate;
    if (!EncodePacket(encoder, K, block_bytes, duplicate) ||
        decoder.Decode(
            duplicate.Id, duplicate.Data.data(), duplicate.Bytes) !=
            Wirehair_Success)
    {
        std::fprintf(stderr, "roundtrip: completed duplicate failed\n");
        return false;
    }
    duplicate.Data[0] ^= 1u;
    if (decoder.Decode(
            duplicate.Id, duplicate.Data.data(), duplicate.Bytes) !=
            Wirehair_Error)
    {
        std::fprintf(stderr,
            "roundtrip: completed corruption was not rejected\n");
        return false;
    }

    std::printf(
        "precode E2E K=%u bb=%u mode=%s delivered=%u overhead=%d: PASS\n",
        K, block_bytes,
        repair_only ?
            (mixed_profile ? "mixed-repair-only" : "repair-only") :
            (mixed_profile ? "mixed-reverse-loss" : "mixed"),
        delivered, (int)delivered - (int)K);
    return true;
}

bool RunDirectLifecycleCase()
{
    const uint32_t K = 64u;
    const uint32_t block_bytes = 29u;
    const uint64_t message_bytes = (uint64_t)K * block_bytes;
    std::vector<uint8_t> message((size_t)message_bytes, 0x5au);
    wirehair_v2::Codec encoder;
    if (encoder.InitializePrecodeEncoder(
            message.data(), message_bytes, block_bytes) != Wirehair_Success)
    {
        return false;
    }

    wirehair_v2::MessagePrecodeDecoder decoder;
    if (decoder.InitializeResult(
            message_bytes, block_bytes, &encoder.Profile()) !=
            Wirehair_Success)
    {
        return false;
    }
    EncodedPacket packet;
    if (!EncodePacket(encoder, 0u, block_bytes, packet) ||
        decoder.DecodeResult(
            packet.Id, packet.Data.data(), packet.Bytes) != Wirehair_NeedMore ||
        decoder.ReceivedCount() != 1u ||
        decoder.DecodeResult(
            packet.Id, packet.Data.data(), packet.Bytes) != Wirehair_NeedMore ||
        decoder.ReceivedCount() != 1u)
    {
        std::fprintf(stderr, "roundtrip: duplicate accounting failed\n");
        return false;
    }
    packet.Data[0] ^= 1u;
    if (decoder.DecodeResult(
            packet.Id, packet.Data.data(), packet.Bytes) !=
            Wirehair_InvalidInput || decoder.ReceivedCount() != 1u)
    {
        std::fprintf(stderr,
            "roundtrip: conflicting duplicate was not rejected\n");
        return false;
    }
    packet.Data[0] ^= 1u;
    if (decoder.DecodeResult(1u, packet.Data.data(), block_bytes - 1u) !=
            Wirehair_InvalidInput || decoder.ReceivedCount() != 1u)
    {
        std::fprintf(stderr, "roundtrip: invalid packet mutated state\n");
        return false;
    }
    std::vector<uint8_t> output((size_t)message_bytes, 0u);
    if (decoder.RecoverResult(output.data(), message_bytes - 1u) !=
            Wirehair_InvalidInput ||
        decoder.RecoverResult(output.data(), message_bytes) !=
            Wirehair_NeedMore)
    {
        std::fprintf(stderr, "roundtrip: recover lifecycle failed\n");
        return false;
    }
    return true;
}

bool RunForcedNeedMoreResumeCase()
{
    const uint32_t K = 1000u;
    const uint32_t block_bytes = 1280u;
    const uint64_t message_bytes = (uint64_t)K * block_bytes;
    std::vector<uint8_t> message((size_t)message_bytes);
    for (size_t i = 0; i < message.size(); ++i) {
        message[i] = (uint8_t)(i * 31u + 9u);
    }
    wirehair_v2::Codec encoder;
    wirehair_v2::MessagePrecodeDecoder decoder;
    if (encoder.InitializePrecodeEncoder(
            message.data(), message_bytes, block_bytes) != Wirehair_Success ||
        decoder.InitializeResult(
            message_bytes, block_bytes, &encoder.Profile()) != Wirehair_Success)
    {
        return false;
    }

    const uint64_t trial_seed =
        UINT64_C(0x5eedf411) ^
        ((uint64_t)K * UINT64_C(0x9e3779b97f4a7c15)) ^
        ((uint64_t)block_bytes * UINT64_C(0xbf58476d1ce4e5b9)) ^
        (UINT64_C(7) * UINT64_C(0xd6e8feb86659fd93));
    TrialRng rng(trial_seed);
    std::vector<uint8_t> block(block_bytes);
    uint32_t block_id = 0u;
    uint32_t delivered = 0u;
    uint32_t last_delivered_id = 0u;
    WirehairResult result = Wirehair_NeedMore;
    double first_solve_ms = 0.0;
    while (delivered < K)
    {
        const uint32_t id = block_id++;
        if (rng.Drop()) {
            continue;
        }
        uint32_t bytes = 0u;
        if (encoder.Encode(id, block.data(), block_bytes, &bytes) !=
            Wirehair_Success)
        {
            return false;
        }
        const std::chrono::steady_clock::time_point begin =
            std::chrono::steady_clock::now();
        result = decoder.DecodeResult(id, block.data(), bytes);
        last_delivered_id = id;
        const double elapsed = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - begin).count();
        ++delivered;
        if (delivered == K) {
            first_solve_ms = elapsed;
        }
    }
    if (result != Wirehair_NeedMore) {
        std::fprintf(stderr,
            "resume: deterministic K-row set result=%d, expected NeedMore\n",
            (int)result);
        return false;
    }
    const uint32_t attempts_before_duplicate = decoder.SolveAttemptCount();
    if (decoder.DecodeResult(
            last_delivered_id, block.data(), block_bytes) !=
            Wirehair_NeedMore ||
        decoder.SolveAttemptCount() != attempts_before_duplicate)
    {
        std::fprintf(stderr,
            "resume: identical duplicate repeated a no-progress solve\n");
        return false;
    }

    double resume_ms = 0.0;
    while (result == Wirehair_NeedMore && delivered < K + 32u)
    {
        const uint32_t id = block_id++;
        if (rng.Drop()) {
            continue;
        }
        uint32_t bytes = 0u;
        if (encoder.Encode(id, block.data(), block_bytes, &bytes) !=
            Wirehair_Success)
        {
            return false;
        }
        const std::chrono::steady_clock::time_point begin =
            std::chrono::steady_clock::now();
        result = decoder.DecodeResult(id, block.data(), bytes);
        resume_ms += std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - begin).count();
        ++delivered;
    }
    std::vector<uint8_t> recovered((size_t)message_bytes, 0u);
    if (result != Wirehair_Success ||
        decoder.RecoverResult(recovered.data(), message_bytes) != Wirehair_Success ||
        recovered != message)
    {
        std::fprintf(stderr,
            "resume: failed result=%d delivered=%u\n",
            (int)result, delivered);
        return false;
    }
    std::printf(
        "precode NeedMore resume K=%u delivered=%u first=%.3fms "
        "resume=%.3fms: PASS\n",
        K, delivered, first_solve_ms, resume_ms);
    return true;
}

bool RunMixedColdRetryCase()
{
    const uint32_t K = 64u;
    const uint32_t block_bytes = 16u;
    const uint64_t message_bytes = (uint64_t)K * block_bytes;
    std::vector<uint8_t> message((size_t)message_bytes);
    for (size_t i = 0; i < message.size(); ++i) {
        message[i] = (uint8_t)(i * 43u + 17u);
    }
    wirehair_v2::MessagePrecodeEncoderOptions options;
    options.Completion =
        wirehair_v2::CompletionField::MixedGF256GF16;
    wirehair_v2::Codec encoder;
    if (encoder.InitializePrecodeEncoder(
            message.data(), message_bytes, block_bytes,
            nullptr, &options) != Wirehair_Success)
    {
        return false;
    }

    std::vector<uint8_t> block(block_bytes);
    // Pinned distinct packet row that leaves the cold solve one rank short
    // for this canonical mixed profile/message fixture.
    for (uint32_t candidate = 1936u; candidate < 1937u; ++candidate)
    {
        wirehair_v2::MessagePrecodeDecoder decoder;
        if (decoder.InitializeResult(
                message_bytes, block_bytes, &encoder.Profile()) !=
                Wirehair_Success)
        {
            return false;
        }
        WirehairResult result = Wirehair_NeedMore;
        for (uint32_t id = 0; id + 1u < K; ++id)
        {
            uint32_t bytes = 0u;
            if (encoder.Encode(id, block.data(), block_bytes, &bytes) !=
                    Wirehair_Success)
            {
                return false;
            }
            result = decoder.DecodeResult(id, block.data(), bytes);
        }
        uint32_t bytes = 0u;
        if (encoder.Encode(
                candidate, block.data(), block_bytes, &bytes) !=
                Wirehair_Success)
        {
            return false;
        }
        result = decoder.DecodeResult(candidate, block.data(), bytes);
        if (result == Wirehair_Success) continue;
        if (result != Wirehair_NeedMore ||
            decoder.HasIncrementalResumeStateForTesting() ||
            decoder.ColdReceiveCapacityBytesForTesting() < message_bytes)
        {
            std::fprintf(stderr,
                "mixed cold retry: deficient state contract failed\n");
            return false;
        }
        const uint32_t attempts_before = decoder.SolveAttemptCount();
        if (decoder.DecodeResult(
                candidate, block.data(), bytes) != Wirehair_NeedMore ||
            decoder.SolveAttemptCount() != attempts_before)
        {
            std::fprintf(stderr,
                "mixed cold retry: duplicate triggered another solve\n");
            return false;
        }
        for (uint32_t repair = K;
             result == Wirehair_NeedMore && repair < K + 128u; ++repair)
        {
            if (repair == candidate) continue;
            if (encoder.Encode(
                    repair, block.data(), block_bytes, &bytes) !=
                    Wirehair_Success)
            {
                return false;
            }
            result = decoder.DecodeResult(repair, block.data(), bytes);
        }
        std::vector<uint8_t> recovered(message.size());
        if (result != Wirehair_Success ||
            decoder.RecoverResult(recovered.data(), recovered.size()) !=
                Wirehair_Success || recovered != message)
        {
            std::fprintf(stderr,
                "mixed cold retry: recovery failed candidate=%u\n",
                candidate);
            return false;
        }
        std::printf(
            "mixed deficient cold retry candidate=%u attempts=%u: PASS\n",
            candidate, decoder.SolveAttemptCount());
        return true;
    }
    std::fprintf(stderr, "mixed cold retry: no deficient fixture found\n");
    return false;
}

bool RunUnauthenticatedCorruptionBoundary()
{
    const uint32_t K = 64u;
    const uint32_t block_bytes = 17u;
    const uint64_t message_bytes = (uint64_t)K * block_bytes;
    std::vector<uint8_t> message((size_t)message_bytes);
    for (size_t i = 0; i < message.size(); ++i) {
        message[i] = (uint8_t)(i * 73u + 11u);
    }

    wirehair_v2::Codec encoder;
    wirehair_v2::MessagePrecodeDecoder decoder;
    if (encoder.InitializePrecodeEncoder(
            message.data(), message_bytes, block_bytes) != Wirehair_Success ||
        decoder.InitializeResult(
            message_bytes, block_bytes, &encoder.Profile()) != Wirehair_Success)
    {
        return false;
    }

    std::vector<uint8_t> block(block_bytes, 0u);
    WirehairResult result = Wirehair_NeedMore;
    for (uint32_t id = 0; id < K; ++id)
    {
        uint32_t bytes = 0u;
        if (encoder.Encode(id, block.data(), block_bytes, &bytes) !=
            Wirehair_Success)
        {
            return false;
        }
        if (id == 3u) {
            block[0] ^= 1u;
        }
        result = decoder.DecodeResult(id, block.data(), bytes);
    }
    std::vector<uint8_t> recovered(message.size(), 0u);
    if (result != Wirehair_Success ||
        decoder.RecoverResult(recovered.data(), recovered.size()) !=
            Wirehair_Success ||
        recovered == message)
    {
        std::fprintf(stderr,
            "integrity boundary: K-row corruption contract changed\n");
        return false;
    }

    uint32_t bytes = 0u;
    if (encoder.Encode(K, block.data(), block_bytes, &bytes) !=
            Wirehair_Success ||
        decoder.DecodeResult(K, block.data(), bytes) != Wirehair_Error)
    {
        std::fprintf(stderr,
            "integrity boundary: overdetermined conflict was not rejected\n");
        return false;
    }
    return true;
}

WirehairResult EncodeDirectPacket(
    wirehair_v2::Codec& encoder,
    uint32_t id,
    uint8_t* block,
    uint32_t block_bytes,
    uint32_t* bytes)
{
    return encoder.Encode(id, block, block_bytes, bytes);
}

WirehairResult EncodeDirectPacket(
    wirehair_v2::MessagePrecodeEncoder& encoder,
    uint32_t id,
    uint8_t* block,
    uint32_t block_bytes,
    uint32_t* bytes)
{
    return encoder.EncodeResult(id, block, block_bytes, bytes);
}

template <typename Encoder>
bool CompleteDirectRepairOnly(
    Encoder& encoder,
    wirehair_v2::MessagePrecodeDecoder& decoder,
    const std::vector<uint8_t>& message,
    uint32_t K,
    uint32_t block_bytes,
    const char* label)
{
    std::vector<uint8_t> block(block_bytes, 0u);
    WirehairResult result = Wirehair_NeedMore;
    uint32_t id = K;
    for (; result == Wirehair_NeedMore && id < 2u * K + 128u; ++id)
    {
        uint32_t bytes = 0u;
        if (EncodeDirectPacket(
                encoder, id, block.data(), block_bytes, &bytes) !=
                Wirehair_Success ||
            bytes != block_bytes)
        {
            std::fprintf(stderr, "contract %s: repair encode failed\n", label);
            return false;
        }
        result = decoder.DecodeResult(id, block.data(), bytes);
    }
    std::vector<uint8_t> recovered(message.size(), 0u);
    if (result != Wirehair_Success ||
        decoder.RecoverResult(recovered.data(), recovered.size()) !=
            Wirehair_Success ||
        recovered != message)
    {
        std::fprintf(stderr,
            "contract %s: repair-only roundtrip failed result=%d id=%u\n",
            label, (int)result, id);
        return false;
    }
    return true;
}

bool SameOptions(
    const wirehair_v2::MessagePrecodeEncoderOptions& a,
    const wirehair_v2::MessagePrecodeEncoderOptions& b)
{
    return a.Architecture == b.Architecture &&
        a.RecoveryMixCount == b.RecoveryMixCount &&
        a.DenseIdentityCorner == b.DenseIdentityCorner &&
        a.AdaptiveDenseTwoAnchor == b.AdaptiveDenseTwoAnchor &&
        a.PrecodeSeedSalt == b.PrecodeSeedSalt &&
        a.RecoveryRowSeedSalt == b.RecoveryRowSeedSalt &&
        a.Completion == b.Completion &&
        a.CacheSystematicSource == b.CacheSystematicSource &&
        a.CacheReceivedSystematicPackets ==
            b.CacheReceivedSystematicPackets;
}

bool SamePrecodeParams(
    const wirehair_v2::PrecodeParams& a,
    const wirehair_v2::PrecodeParams& b)
{
    return a.BlockCount == b.BlockCount &&
        a.Staircase == b.Staircase &&
        a.DenseRows == b.DenseRows &&
        a.HeavyRows == b.HeavyRows &&
        a.SourceHits == b.SourceHits &&
        a.Field == b.Field &&
        a.DegreeBalancedStaircase == b.DegreeBalancedStaircase &&
        a.StaircaseDegreeScale == b.StaircaseDegreeScale &&
        a.DenseIdentityCorner == b.DenseIdentityCorner &&
        a.DenseTwoAnchor == b.DenseTwoAnchor &&
        a.DenseTwoAnchorPhase == b.DenseTwoAnchorPhase &&
        a.SegmentedDenseAnchors == b.SegmentedDenseAnchors &&
        a.HeavyFamily == b.HeavyFamily &&
        a.Seed == b.Seed;
}

bool CheckExplicitDiagnostics(
    const wirehair_v2::MessagePrecodeEncoder& encoder,
    const wirehair_v2::MessagePrecodeDecoder& decoder,
    const wirehair_v2::ExplicitMessagePrecodeConfigForTesting& config,
    uint32_t K,
    uint32_t block_bytes,
    const char* label)
{
    const auto profile_matches =
        [&](const wirehair_v2::SeedProfile& profile) {
            return profile.BlockCount == K &&
                profile.BlockBytes == block_bytes &&
                profile.V2SeedSelected &&
                profile.V2SeedAttempt == 0u &&
                profile.V2PrecodeContractVersion == 0u &&
                profile.V2PacketRowContractVersion == 0u &&
                profile.V2SeedPolicy ==
                    wirehair_v2::V2SeedDerivation::RawUniform &&
                profile.V2StaircaseCount == config.Params.Staircase &&
                profile.V2DenseRowCount == config.Params.DenseRows &&
                profile.V2HeavyRowCount == config.Params.HeavyRows &&
                profile.V2CompletionField == config.Params.Field &&
                profile.V2SourceHits == config.Params.SourceHits &&
                profile.V2PrecodeSeed == config.Params.Seed &&
                profile.V2PacketPeelSeed == config.Packet.PeelSeed &&
                profile.V2RecoveryMixCount == config.Packet.MixCount &&
                profile.V2DenseIdentityCorner ==
                    config.Params.DenseIdentityCorner &&
                profile.V2DenseTwoAnchor == config.Params.DenseTwoAnchor &&
                !profile.V2AdaptiveDenseTwoAnchor &&
                profile.V2PrecodeSeedSalt == 0u &&
                profile.V2RecoveryRowSeedSalt == 0u;
        };
    const auto options_match =
        [&](const wirehair_v2::MessagePrecodeEncoderOptions& options) {
            return options.RecoveryMixCount == config.Packet.MixCount &&
                options.DenseIdentityCorner ==
                    config.Params.DenseIdentityCorner &&
                !options.AdaptiveDenseTwoAnchor &&
                options.PrecodeSeedSalt == 0u &&
                options.RecoveryRowSeedSalt == 0u &&
                options.Completion == config.Params.Field &&
                options.CacheSystematicSource ==
                    config.CacheSystematicSource &&
                options.CacheReceivedSystematicPackets ==
                    config.CacheReceivedSystematicPackets;
        };
    if (!encoder.IsInitialized() || !decoder.IsInitialized() ||
        !profile_matches(encoder.Profile()) ||
        !profile_matches(decoder.Profile()) ||
        !options_match(encoder.Options()) ||
        !options_match(decoder.Options()) ||
        encoder.DiagnosticIdentityForTesting() !=
            wirehair_v2::MessagePrecodeDiagnosticIdentityForTesting::
                ExplicitUnknownArchitecture ||
        decoder.DiagnosticIdentityForTesting() !=
            wirehair_v2::MessagePrecodeDiagnosticIdentityForTesting::
                ExplicitUnknownArchitecture ||
        encoder.PinnedEquationStateForTesting().Precode !=
            config.EquationState.Precode ||
        encoder.PinnedEquationStateForTesting().PacketRows !=
            config.EquationState.PacketRows ||
        encoder.PinnedEquationStateForTesting().MixedBandTrackingX !=
            config.EquationState.MixedBandTrackingX ||
        decoder.PinnedEquationStateForTesting().Precode !=
            config.EquationState.Precode ||
        decoder.PinnedEquationStateForTesting().PacketRows !=
            config.EquationState.PacketRows ||
        decoder.PinnedEquationStateForTesting().MixedBandTrackingX !=
            config.EquationState.MixedBandTrackingX ||
        !SamePrecodeParams(
            encoder.BlockEncoder().System().Params, config.Params) ||
        !SamePrecodeParams(decoder.System().Params, config.Params) ||
        encoder.SolveStats().PacketSeedAttempt != 0u ||
        decoder.PacketSeedAttempt() != 0u ||
        decoder.PacketPeelSeed() != config.Packet.PeelSeed ||
        encoder.HasSystematicSourceCache() !=
            config.CacheSystematicSource ||
        decoder.HasSystematicPacketCache() !=
            config.CacheReceivedSystematicPackets)
    {
        std::fprintf(stderr,
            "explicit %s: diagnostic/configuration state mismatch\n", label);
        return false;
    }
    return true;
}

bool RunOneExplicitConfiguration(
    wirehair_v2::ExplicitMessagePrecodeConfigForTesting config,
    uint32_t block_bytes,
    uint32_t tail_bytes,
    const char* label)
{
    const uint32_t K = config.Params.BlockCount;
    const uint64_t message_bytes =
        (uint64_t)(K - 1u) * block_bytes + tail_bytes;
    std::vector<uint8_t> message((size_t)message_bytes);
    for (size_t i = 0u; i < message.size(); ++i) {
        message[i] = (uint8_t)(i * 73u + K * 11u + block_bytes);
    }
    wirehair_v2::MessagePrecodeEncoder encoder;
    wirehair_v2::MessagePrecodeDecoder decoder;
    const WirehairResult encode_init =
        encoder.InitializeExplicitResultForTesting(
            message.data(), message_bytes, block_bytes, config);
    const WirehairResult decode_init =
        decoder.InitializeExplicitResultForTesting(
            message_bytes, block_bytes, config);
    if (encode_init != Wirehair_Success ||
        decode_init != Wirehair_Success ||
        !CheckExplicitDiagnostics(
            encoder, decoder, config, K, block_bytes, label) ||
        !CompleteDirectRepairOnly(
            encoder, decoder, message, K, block_bytes, label))
    {
        std::fprintf(stderr,
            "explicit %s: roundtrip failed enc=%d dec=%d\n",
            label, (int)encode_init, (int)decode_init);
        return false;
    }
    return true;
}

bool RunExplicitEquationStatePinningCases()
{
    struct RestoreHooks
    {
        ~RestoreHooks()
        {
            wirehair_v2::ClearMixedBandTrackingXForTesting();
            (void)wirehair_v2::SetMixedCoefficientGeometryForTesting(
                wirehair_v2::MixedCoefficientGeometry::FrozenPowerX);
            (void)wirehair_v2::SetMixedGF16RowsForTesting(
                wirehair_v2::kMixedGF16Rows);
            (void)wirehair_v2::SetMixedGF256RowsForTesting(
                wirehair_v2::kMixedGF256Rows);
            (void)wirehair_v2::SetMixedCoefficientPeriodForTesting(
                wirehair_v2::kMixedCoefficientPeriod);
            (void)wirehair_v2::SetPacketRowSeedMultiplierForTesting(1u);
            wirehair_v2::SetPacketRowSeedAvalancheForTesting(false);
            wirehair_v2::SetOddPacketPeelSeedXorForTesting(0u);
            wirehair_v2::ClearPeelDegreesForTesting();
            (void)wirehair_v2::SetMixedResidueBucketModeForTesting(
                wirehair_v2::MixedResidueBucketMode::Automatic);
        }
    } restore;

    static const uint32_t K = 43u;
    static const uint32_t block_bytes = 18u;
    static const uint32_t tail_bytes = 7u;
    const uint64_t message_bytes =
        (uint64_t)(K - 1u) * block_bytes + tail_bytes;
    std::vector<uint8_t> message((size_t)message_bytes);
    for (size_t i = 0u; i < message.size(); ++i) {
        message[i] = (uint8_t)(i * 29u + 0x53u);
    }

    if (!wirehair_v2::SetMixedGF16RowsForTesting(0u) ||
        !wirehair_v2::SetMixedGF256RowsForTesting(8u))
    {
        return false;
    }
    wirehair_v2::SetMixedBandTrackingXForTesting(true);

    wirehair_v2::ExplicitMessagePrecodeConfigForTesting config;
    config.Params = wirehair_v2::MakeMixedParams(
        K, UINT64_C(0xb1716d94ac2058e3));
    config.Params.Staircase = 9u;
    config.Params.DenseRows = 4u;
    config.Params.SourceHits = 2u;
    config.Params.DegreeBalancedStaircase = true;
    config.Params.StaircaseDegreeScale = 8.5;
    config.Packet.PeelSeed = UINT32_C(0x41f06c33);
    config.Packet.MixCount = 3u;
    if (config.Params.HeavyRows != 8u ||
        !config.PinActiveEquationStateForTesting())
    {
        return false;
    }

    wirehair_v2::MessagePrecodeEncoder encoder;
    wirehair_v2::MessagePrecodeDecoder decoder;
    if (encoder.InitializeExplicitResultForTesting(
            message.data(), message_bytes, block_bytes, config) !=
            Wirehair_Success ||
        decoder.InitializeExplicitResultForTesting(
            message_bytes, block_bytes, config) != Wirehair_Success)
    {
        std::fprintf(stderr,
            "explicit pinning: initial 8+0 construction failed\n");
        return false;
    }

    const uint32_t probe_id = 3u * K + 17u;
    std::vector<uint8_t> probe(block_bytes, 0u);
    uint32_t probe_bytes = 0u;
    if (encoder.EncodeResult(
            probe_id, probe.data(), block_bytes, &probe_bytes) !=
            Wirehair_Success ||
        !CompleteDirectRepairOnly(
            encoder, decoder, message, K, block_bytes, "pinning"))
    {
        return false;
    }

    // The same S/D/H and total H=8 used to accept both 8+0 and 6+2 ambient
    // row splits even though they generate different mixed equations.
    if (!wirehair_v2::SetMixedGF256RowsForTesting(6u) ||
        !wirehair_v2::SetMixedGF16RowsForTesting(2u))
    {
        return false;
    }
    wirehair_v2::MessagePrecodeEncoder stale_encoder;
    wirehair_v2::MessagePrecodeDecoder stale_decoder;
    if (config.EquationState.MatchesActive(config.Params) ||
        stale_encoder.InitializeExplicitResultForTesting(
            message.data(), message_bytes, block_bytes, config) !=
            Wirehair_InvalidInput ||
        stale_decoder.InitializeExplicitResultForTesting(
            message_bytes, block_bytes, config) != Wirehair_InvalidInput ||
        encoder.InitializeExplicitResultForTesting(
            message.data(), message_bytes, block_bytes, config) !=
            Wirehair_InvalidInput ||
        decoder.InitializeExplicitResultForTesting(
            message_bytes, block_bytes, config) != Wirehair_InvalidInput ||
        !encoder.IsInitialized() || !decoder.IsDecoded())
    {
        std::fprintf(stderr,
            "explicit pinning: stale 8+0 descriptor accepted at 6+2\n");
        return false;
    }

    wirehair_v2::ExplicitMessagePrecodeConfigForTesting split = config;
    if (!split.PinActiveEquationStateForTesting() ||
        split.EquationState.Precode == config.EquationState.Precode)
    {
        std::fprintf(stderr,
            "explicit pinning: same-total split identity did not change\n");
        return false;
    }
    wirehair_v2::MessagePrecodeEncoder split_encoder;
    wirehair_v2::MessagePrecodeDecoder split_decoder;
    if (split_encoder.InitializeExplicitResultForTesting(
            message.data(), message_bytes, block_bytes, split) !=
            Wirehair_Success ||
        split_decoder.InitializeExplicitResultForTesting(
            message_bytes, block_bytes, split) != Wirehair_Success)
    {
        std::fprintf(stderr,
            "explicit pinning: repinned 6+2 construction failed\n");
        return false;
    }
    const uint64_t columns = (uint64_t)K + config.Params.Staircase +
        config.Params.DenseRows + config.Params.HeavyRows;
    if (std::memcmp(
            encoder.IntermediateBlocks(),
            split_encoder.IntermediateBlocks(),
            (size_t)columns * block_bytes) == 0)
    {
        std::fprintf(stderr,
            "explicit pinning: 8+0 and 6+2 intermediates unexpectedly match\n");
        return false;
    }

    if (!wirehair_v2::SetMixedGF16RowsForTesting(0u) ||
        !wirehair_v2::SetMixedGF256RowsForTesting(8u) ||
        !config.EquationState.MatchesActive(config.Params))
    {
        return false;
    }

    const auto reject_drift = [&](const char* label) {
        std::vector<uint8_t> encoded(block_bytes, 0xa5u);
        uint32_t encoded_bytes = UINT32_C(0xdeadbeef);
        uint64_t encoded_ops = UINT64_MAX;
        const uint32_t received = decoder.ReceivedCount();
        std::vector<uint8_t> recovered(message.size(), 0x5au);
        if (encoder.EncodeResult(
                probe_id,
                encoded.data(),
                block_bytes,
                &encoded_bytes,
                &encoded_ops) != Wirehair_Error ||
            encoded_bytes != UINT32_C(0xdeadbeef) ||
            encoded_ops != UINT64_MAX ||
            !std::all_of(encoded.begin(), encoded.end(),
                [](uint8_t value) { return value == 0xa5u; }) ||
            decoder.DecodeResult(
                probe_id, probe.data(), probe_bytes) != Wirehair_Error ||
            decoder.ReceivedCount() != received ||
            decoder.RecoverResult(
                recovered.data(), recovered.size()) != Wirehair_Error ||
            !std::all_of(recovered.begin(), recovered.end(),
                [](uint8_t value) { return value == 0x5au; }))
        {
            std::fprintf(stderr,
                "explicit pinning: %s drift was not fail-closed\n", label);
            return false;
        }
        return true;
    };

    wirehair_v2::SetMixedBandTrackingXForTesting(false);
    if (config.EquationState.MatchesActive(config.Params) ||
        !reject_drift("tracking-X"))
    {
        return false;
    }
    wirehair_v2::SetMixedBandTrackingXForTesting(true);
    if (!config.EquationState.MatchesActive(config.Params) ||
        !wirehair_v2::SetPacketRowSeedMultiplierForTesting(3u) ||
        config.EquationState.MatchesActive(config.Params) ||
        !reject_drift("packet-row"))
    {
        return false;
    }
    if (!wirehair_v2::SetPacketRowSeedMultiplierForTesting(1u) ||
        !config.EquationState.MatchesActive(config.Params))
    {
        return false;
    }

    // Bucket selection changes only how the same mixed RHS is accumulated.
    // It must remain outside the equation identity and preserve the solve.
    if (!wirehair_v2::SetMixedResidueBucketModeForTesting(
            wirehair_v2::MixedResidueBucketMode::Separate) ||
        !config.EquationState.MatchesActive(config.Params))
    {
        return false;
    }
    wirehair_v2::MessagePrecodeEncoder bucket_encoder;
    if (bucket_encoder.InitializeExplicitResultForTesting(
            message.data(), message_bytes, block_bytes, config) !=
            Wirehair_Success ||
        std::memcmp(
            encoder.IntermediateBlocks(),
            bucket_encoder.IntermediateBlocks(),
            (size_t)columns * block_bytes) != 0 ||
        !wirehair_v2::SetMixedResidueBucketModeForTesting(
            wirehair_v2::MixedResidueBucketMode::Automatic) ||
        !config.EquationState.MatchesActive(config.Params))
    {
        std::fprintf(stderr,
            "explicit pinning: execution-only bucket mode changed identity "
            "or equations\n");
        return false;
    }

    std::vector<uint8_t> recovered(message.size(), 0u);
    if (decoder.RecoverResult(recovered.data(), recovered.size()) !=
            Wirehair_Success ||
        recovered != message)
    {
        std::fprintf(stderr,
            "explicit pinning: restored hook state did not resume\n");
        return false;
    }
    std::printf(
        "explicit equation-state pinning 8+0->6+2, tracking-X, and packet "
        "drift: PASS\n");
    return true;
}

bool RunExplicitTrackingXConcurrencyCase()
{
    struct RestoreHooks
    {
        ~RestoreHooks()
        {
            wirehair_v2::ClearMixedBandTrackingXForTesting();
            (void)wirehair_v2::SetMixedGF16RowsForTesting(
                wirehair_v2::kMixedGF16Rows);
            (void)wirehair_v2::SetMixedGF256RowsForTesting(
                wirehair_v2::kMixedGF256Rows);
        }
    } restore;

    static const uint32_t K = 12000u;
    static const uint32_t block_bytes = 2u;
    std::vector<uint8_t> message((size_t)K * block_bytes);
    for (size_t i = 0u; i < message.size(); ++i) {
        message[i] = (uint8_t)(i * 73u + 19u);
    }
    if (!wirehair_v2::SetMixedGF16RowsForTesting(0u) ||
        !wirehair_v2::SetMixedGF256RowsForTesting(8u))
    {
        return false;
    }

    wirehair_v2::ExplicitMessagePrecodeConfigForTesting base;
    base.Params = wirehair_v2::MakeMixedParams(
        K, UINT64_C(0xb1716d94ac2058e3));
    base.Params.DenseRows = 4u;
    base.Packet.PeelSeed = UINT32_C(0x41f06c33);
    base.Packet.MixCount = 3u;

    // First prove the scoped primitive itself.  Algebra must follow the
    // innermost snapshot, while the equation identity continues to describe
    // the caller-visible ambient state.  The thrown marker verifies that
    // unwinding an inner operation restores the outer operation's snapshot.
    wirehair_v2::SetMixedBandTrackingXForTesting(false);
    const wirehair_v2::MixedCoefficientRows* const frozen =
        wirehair_v2::GetMixedCoefficientRows();
    if (!frozen) {
        return false;
    }
    const wirehair_v2::MixedCoefficientRows frozen_copy = *frozen;
    wirehair_v2::SetMixedBandTrackingXForTesting(true);
    const wirehair_v2::MixedCoefficientRows* const tracked =
        wirehair_v2::GetMixedCoefficientRows();
    if (!tracked ||
        std::memcmp(
            &frozen_copy, tracked,
            sizeof(wirehair_v2::MixedCoefficientRows)) == 0)
    {
        std::fprintf(stderr,
            "explicit tracking-X scope: stable witnesses match\n");
        return false;
    }
    const wirehair_v2::MixedCoefficientRows tracked_copy = *tracked;
    const uint64_t ambient_fingerprint =
        wirehair_v2::ActivePrecodeEquationStateFingerprintForTesting(
            base.Params);
    wirehair_v2::SetMixedBandTrackingXForTesting(false);
    wirehair_v2::ExplicitMessagePrecodeConfigForTesting stable_identity_off =
        base;
    if (!stable_identity_off.PinActiveEquationStateForTesting()) {
        return false;
    }
    wirehair_v2::SetMixedBandTrackingXForTesting(true);
    wirehair_v2::ExplicitMessagePrecodeConfigForTesting stable_identity_on =
        base;
    if (!stable_identity_on.PinActiveEquationStateForTesting() ||
        stable_identity_off.EquationState.Precode !=
            stable_identity_on.EquationState.Precode ||
        stable_identity_off.EquationState.PacketRows !=
            stable_identity_on.EquationState.PacketRows ||
        stable_identity_off.EquationState.MixedBandTrackingX ||
        !stable_identity_on.EquationState.MixedBandTrackingX)
    {
        std::fprintf(stderr,
            "explicit tracking-X identity: stable split is inconsistent\n");
        return false;
    }
    const auto identity_matches = [](
        const wirehair_v2::ExplicitEquationStateIdentityForTesting& left,
        const wirehair_v2::ExplicitEquationStateIdentityForTesting& right)
    {
        return left.Precode == right.Precode &&
            left.PacketRows == right.PacketRows &&
            left.MixedBandTrackingX == right.MixedBandTrackingX;
    };
    std::atomic<bool> stop_identity_toggle{false};
    std::thread identity_toggler([&] {
        bool enabled = false;
        while (!stop_identity_toggle.load(std::memory_order_acquire))
        {
            wirehair_v2::SetMixedBandTrackingXForTesting(enabled);
            enabled = !enabled;
        }
    });
    bool identity_ok = true;
    for (uint32_t trial = 0u; trial < 200000u; ++trial)
    {
        wirehair_v2::ExplicitMessagePrecodeConfigForTesting sampled = base;
        if (!sampled.PinActiveEquationStateForTesting() ||
            (
                !identity_matches(
                    sampled.EquationState,
                    stable_identity_off.EquationState) &&
                !identity_matches(
                    sampled.EquationState,
                    stable_identity_on.EquationState)
            ))
        {
            identity_ok = false;
            break;
        }
    }
    stop_identity_toggle.store(true, std::memory_order_release);
    identity_toggler.join();
    if (!identity_ok)
    {
        std::fprintf(stderr,
            "explicit tracking-X identity: concurrent pin was torn\n");
        return false;
    }
    wirehair_v2::SetMixedBandTrackingXForTesting(true);
    {
        wirehair_v2::ScopedMixedBandTrackingXForTesting outer(true, false);
        wirehair_v2::ScopedMixedBandTrackingXForTesting inactive(false, true);
        if (!wirehair_v2::AmbientMixedBandTrackingXForTesting() ||
            wirehair_v2::
                ActivePrecodeEquationStateFingerprintForTesting(base.Params) !=
                ambient_fingerprint ||
            std::memcmp(
                &frozen_copy, wirehair_v2::GetMixedCoefficientRows(),
                sizeof(frozen_copy)) != 0)
        {
            std::fprintf(stderr,
                "explicit tracking-X scope: outer snapshot mismatch\n");
            return false;
        }
        bool caught = false;
        try
        {
            wirehair_v2::ScopedMixedBandTrackingXForTesting inner(true, true);
            if (std::memcmp(
                    &tracked_copy, wirehair_v2::GetMixedCoefficientRows(),
                    sizeof(tracked_copy)) != 0)
            {
                std::fprintf(stderr,
                    "explicit tracking-X scope: inner snapshot mismatch\n");
                return false;
            }
            throw 17;
        }
        catch (int marker)
        {
            caught = marker == 17;
        }
        if (!caught ||
            std::memcmp(
                &frozen_copy, wirehair_v2::GetMixedCoefficientRows(),
                sizeof(frozen_copy)) != 0)
        {
            std::fprintf(stderr,
                "explicit tracking-X scope: unwind did not restore outer\n");
            return false;
        }
    }
    if (!wirehair_v2::AmbientMixedBandTrackingXForTesting() ||
        std::memcmp(
            &tracked_copy, wirehair_v2::GetMixedCoefficientRows(),
            sizeof(tracked_copy)) != 0)
    {
        std::fprintf(stderr,
            "explicit tracking-X scope: outer teardown did not restore "
            "ambient state\n");
        return false;
    }

    const auto initialize_at = [&](
        bool tracking,
        wirehair_v2::MessagePrecodeEncoder& encoder)
    {
        wirehair_v2::SetMixedBandTrackingXForTesting(tracking);
        wirehair_v2::ExplicitMessagePrecodeConfigForTesting config = base;
        return config.PinActiveEquationStateForTesting() &&
            encoder.InitializeExplicitResultForTesting(
                message.data(), message.size(), block_bytes, config) ==
                Wirehair_Success;
    };

    wirehair_v2::MessagePrecodeEncoder stable_on;
    wirehair_v2::MessagePrecodeEncoder stable_off;
    if (!initialize_at(true, stable_on) ||
        !initialize_at(false, stable_off))
    {
        std::fprintf(stderr,
            "explicit tracking-X concurrency: stable construction failed\n");
        return false;
    }
    const size_t intermediate_bytes = (size_t)(
        K + base.Params.Staircase + base.Params.DenseRows +
        base.Params.HeavyRows) * block_bytes;
    if (std::memcmp(
            stable_on.IntermediateBlocks(),
            stable_off.IntermediateBlocks(),
            intermediate_bytes) == 0)
    {
        std::fprintf(stderr,
            "explicit tracking-X concurrency: stable witnesses match\n");
        return false;
    }

    std::vector<uint8_t> transaction_off_expected(block_bytes, 0u);
    uint32_t transaction_off_expected_bytes = 0u;
    if (stable_off.EncodeResult(
            K,
            transaction_off_expected.data(),
            block_bytes,
            &transaction_off_expected_bytes) != Wirehair_Success)
    {
        std::fprintf(stderr,
            "explicit tracking-X transaction: frozen setup failed\n");
        return false;
    }

    wirehair_v2::SetMixedBandTrackingXForTesting(true);
    wirehair_v2::ExplicitMessagePrecodeConfigForTesting
        transaction_config = base;
    std::vector<uint8_t> transaction_expected(block_bytes, 0u);
    uint32_t transaction_expected_bytes = 0u;
    if (!transaction_config.PinActiveEquationStateForTesting() ||
        stable_on.EncodeResult(
            K,
            transaction_expected.data(),
            block_bytes,
            &transaction_expected_bytes) != Wirehair_Success)
    {
        std::fprintf(stderr,
            "explicit tracking-X transaction: setup failed\n");
        return false;
    }
    std::vector<uint8_t> transaction_actual(block_bytes, 0xa5u);
    uint32_t transaction_actual_bytes = UINT32_MAX;
    {
        wirehair_v2::ScopedExplicitEquationStateTransactionForTesting
            transaction(transaction_config);
        if (!transaction.IsValid()) {
            return false;
        }
        wirehair_v2::SetMixedBandTrackingXForTesting(false);
        if (stable_on.EncodeResult(
                K,
                transaction_actual.data(),
                block_bytes,
                &transaction_actual_bytes) != Wirehair_Success ||
            transaction_actual_bytes != transaction_expected_bytes ||
            transaction_actual != transaction_expected)
        {
            std::fprintf(stderr,
                "explicit tracking-X transaction: pinned encode drifted\n");
            return false;
        }

        // A transaction is an authorization for the descriptor context it
        // validated, not for an arbitrary Params value carrying the same
        // copied ambient identity.  In particular, StaircaseDegreeScale is a
        // fallback input to the precode fingerprint.
        wirehair_v2::SetMixedBandTrackingXForTesting(true);
        wirehair_v2::ExplicitMessagePrecodeConfigForTesting stale_config =
            transaction_config;
        stale_config.Params.StaircaseDegreeScale += 1.0;
        wirehair_v2::MessagePrecodeDecoder stale_decoder;
        if (stale_config.EquationState.MatchesActive(stale_config.Params) ||
            stale_decoder.InitializeExplicitResultForTesting(
                message.size(), block_bytes, stale_config) !=
                    Wirehair_InvalidInput ||
            stale_decoder.IsInitialized())
        {
            std::fprintf(stderr,
                "explicit tracking-X transaction: stale Params inherited "
                "authorization\n");
            return false;
        }

        // An invalid inner scope must suppress, rather than inherit, the
        // outer transaction.  Its teardown must then restore the outer scope.
        {
            wirehair_v2::ScopedExplicitEquationStateTransactionForTesting
                invalid_nested(stale_config);
            if (invalid_nested.IsValid()) {
                return false;
            }
            wirehair_v2::SetMixedBandTrackingXForTesting(false);
            std::fill(
                transaction_actual.begin(),
                transaction_actual.end(),
                uint8_t{0x3c});
            transaction_actual_bytes = UINT32_MAX;
            if (stable_on.EncodeResult(
                    K,
                    transaction_actual.data(),
                    block_bytes,
                    &transaction_actual_bytes) != Wirehair_Error ||
                transaction_actual_bytes != UINT32_MAX ||
                transaction_actual != std::vector<uint8_t>(
                    block_bytes, uint8_t{0x3c}))
            {
                std::fprintf(stderr,
                    "explicit tracking-X transaction: invalid nested scope "
                    "inherited outer authorization\n");
                return false;
            }
        }
        std::fill(
            transaction_actual.begin(),
            transaction_actual.end(),
            uint8_t{0xa5});
        transaction_actual_bytes = UINT32_MAX;
        if (stable_on.EncodeResult(
                K,
                transaction_actual.data(),
                block_bytes,
                &transaction_actual_bytes) != Wirehair_Success ||
            transaction_actual_bytes != transaction_expected_bytes ||
            transaction_actual != transaction_expected)
        {
            std::fprintf(stderr,
                "explicit tracking-X transaction: invalid nested teardown "
                "did not restore outer authorization\n");
            return false;
        }

        // A valid inner transaction authorizes only its own identity.  Throw
        // through it to exercise both the authorization-stack and tracking-X
        // snapshot destructors, then prove that the outer transaction resumes.
        wirehair_v2::ExplicitMessagePrecodeConfigForTesting nested_config =
            base;
        if (!nested_config.PinActiveEquationStateForTesting()) {
            return false;
        }
        bool caught = false;
        try
        {
            wirehair_v2::ScopedExplicitEquationStateTransactionForTesting
                nested(nested_config);
            if (!nested.IsValid()) {
                return false;
            }
            std::fill(
                transaction_actual.begin(),
                transaction_actual.end(),
                uint8_t{0x69});
            transaction_actual_bytes = UINT32_MAX;
            if (stable_on.EncodeResult(
                    K,
                    transaction_actual.data(),
                    block_bytes,
                    &transaction_actual_bytes) != Wirehair_Error ||
                transaction_actual_bytes != UINT32_MAX ||
                transaction_actual != std::vector<uint8_t>(
                    block_bytes, uint8_t{0x69}))
            {
                std::fprintf(stderr,
                    "explicit tracking-X transaction: nested identity "
                    "authorized the outer endpoint\n");
                return false;
            }
            if (stable_off.EncodeResult(
                    K,
                    transaction_actual.data(),
                    block_bytes,
                    &transaction_actual_bytes) != Wirehair_Success ||
                transaction_actual_bytes !=
                    transaction_off_expected_bytes ||
                transaction_actual != transaction_off_expected)
            {
                std::fprintf(stderr,
                    "explicit tracking-X transaction: nested endpoint "
                    "authorization failed\n");
                return false;
            }
            throw 23;
        }
        catch (int marker)
        {
            caught = marker == 23;
        }
        std::fill(
            transaction_actual.begin(),
            transaction_actual.end(),
            uint8_t{0x96});
        transaction_actual_bytes = UINT32_MAX;
        if (!caught ||
            stable_on.EncodeResult(
                K,
                transaction_actual.data(),
                block_bytes,
                &transaction_actual_bytes) != Wirehair_Success ||
            transaction_actual_bytes != transaction_expected_bytes ||
            transaction_actual != transaction_expected)
        {
            std::fprintf(stderr,
                "explicit tracking-X transaction: nested unwind did not "
                "restore outer authorization\n");
            return false;
        }
    }
    std::fill(
        transaction_actual.begin(), transaction_actual.end(), uint8_t{0x5a});
    transaction_actual_bytes = UINT32_MAX;
    if (stable_on.EncodeResult(
            K,
            transaction_actual.data(),
            block_bytes,
            &transaction_actual_bytes) != Wirehair_Error ||
        transaction_actual_bytes != UINT32_MAX ||
        transaction_actual !=
            std::vector<uint8_t>(block_bytes, uint8_t{0x5a}))
    {
        std::fprintf(stderr,
            "explicit tracking-X transaction: teardown retained authority\n");
        return false;
    }
    wirehair_v2::SetMixedBandTrackingXForTesting(true);
    {
        // A rejected transaction is not an equation-state override.  In
        // particular, its untrusted tracking bit must not perturb ordinary
        // hook consumers for the lifetime of the invalid scope.
        wirehair_v2::ExplicitMessagePrecodeConfigForTesting invalid =
            transaction_config;
        invalid.EquationState.MixedBandTrackingX = false;
        wirehair_v2::ScopedExplicitEquationStateTransactionForTesting
            invalid_transaction(invalid);
        if (invalid_transaction.IsValid() ||
            std::memcmp(
                &tracked_copy, wirehair_v2::GetMixedCoefficientRows(),
                sizeof(tracked_copy)) != 0)
        {
            std::fprintf(stderr,
                "explicit tracking-X transaction: invalid scope changed "
                "effective hook state\n");
            return false;
        }
    }
    if (std::memcmp(
            &tracked_copy, wirehair_v2::GetMixedCoefficientRows(),
            sizeof(tracked_copy)) != 0)
    {
        std::fprintf(stderr,
            "explicit tracking-X transaction: invalid scope teardown "
            "changed ambient hook state\n");
        return false;
    }

    bool observed_success = false;
    for (uint32_t delay_us = 0u; delay_us <= 5000u; delay_us += 100u)
    {
        wirehair_v2::SetMixedBandTrackingXForTesting(true);
        wirehair_v2::ExplicitMessagePrecodeConfigForTesting config = base;
        if (!config.PinActiveEquationStateForTesting()) {
            return false;
        }
        std::atomic<bool> entered{false};
        std::thread toggler([&] {
            entered.store(true, std::memory_order_release);
            std::this_thread::sleep_for(
                std::chrono::microseconds(delay_us));
            wirehair_v2::SetMixedBandTrackingXForTesting(false);
            std::this_thread::sleep_for(
                std::chrono::microseconds(500u));
            wirehair_v2::SetMixedBandTrackingXForTesting(true);
        });
        while (!entered.load(std::memory_order_acquire)) {
        }

        wirehair_v2::MessagePrecodeEncoder raced;
        const WirehairResult result =
            raced.InitializeExplicitResultForTesting(
                message.data(), message.size(), block_bytes, config);
        toggler.join();
        if (result == Wirehair_InvalidInput) {
            continue;
        }
        if (result != Wirehair_Success ||
            std::memcmp(
                raced.IntermediateBlocks(),
                stable_on.IntermediateBlocks(),
                intermediate_bytes) != 0)
        {
            std::fprintf(stderr,
                "explicit tracking-X concurrency: mixed-state publication "
                "at delay=%u result=%d\n",
                delay_us, (int)result);
            return false;
        }
        observed_success = true;
    }
    if (!observed_success) {
        std::fprintf(stderr,
            "explicit tracking-X concurrency: no successful race witness\n");
        return false;
    }
    std::printf(
        "explicit tracking-X concurrent ABA snapshot: PASS\n");
    return true;
}

bool RunExplicitConfigurationCases()
{
    static const uint32_t K = 64u;
    static const uint32_t block_bytes = 16u;
    static const uint32_t tail_bytes = 7u;
    const uint64_t message_bytes =
        (uint64_t)(K - 1u) * block_bytes + tail_bytes;
    std::vector<uint8_t> message((size_t)message_bytes);
    for (size_t i = 0u; i < message.size(); ++i) {
        message[i] = (uint8_t)(i * 47u + 19u);
    }

    wirehair_v2::MessagePrecodeEncoderOptions normal_options;
    normal_options.Architecture =
        wirehair_v2::V2PrecodeArchitecture::SmallBandD4;
    normal_options.Completion =
        wirehair_v2::CompletionField::MixedGF256GF16;
    wirehair_v2::MessagePrecodeEncoder normal;
    if (normal.InitializeResult(
            message.data(), message_bytes, block_bytes,
            nullptr, &normal_options) != Wirehair_Success)
    {
        std::fprintf(stderr, "explicit canonical: normal init failed\n");
        return false;
    }

    wirehair_v2::ExplicitMessagePrecodeConfigForTesting config;
    config.Params = normal.BlockEncoder().System().Params;
    config.Packet.PeelSeed = normal.Profile().V2PacketPeelSeed;
    config.Packet.MixCount = normal.Profile().V2RecoveryMixCount;
    config.CacheSystematicSource = false;
    config.CacheReceivedSystematicPackets = false;
    if (!config.PinActiveEquationStateForTesting() ||
        normal.DiagnosticIdentityForTesting() !=
            wirehair_v2::MessagePrecodeDiagnosticIdentityForTesting::
                NamedContract)
    {
        std::fprintf(stderr,
            "explicit canonical: equation-state pin failed\n");
        return false;
    }

    wirehair_v2::MessagePrecodeEncoder explicit_encoder;
    wirehair_v2::MessagePrecodeDecoder explicit_decoder;
    wirehair_v2::ExplicitMessagePrecodeConfigForTesting unpinned = config;
    unpinned.EquationState = {};
    wirehair_v2::MessagePrecodeEncoder unpinned_encoder;
    wirehair_v2::MessagePrecodeDecoder unpinned_decoder;
    if (unpinned_encoder.InitializeExplicitResultForTesting(
            message.data(), message_bytes, block_bytes, unpinned) !=
            Wirehair_InvalidInput ||
        unpinned_decoder.InitializeExplicitResultForTesting(
            message_bytes, block_bytes, unpinned) != Wirehair_InvalidInput)
    {
        std::fprintf(stderr,
            "explicit canonical: unpinned descriptor was accepted\n");
        return false;
    }
    if (explicit_encoder.InitializeExplicitResultForTesting(
            message.data(), message_bytes, block_bytes, config) !=
            Wirehair_Success ||
        explicit_decoder.InitializeExplicitResultForTesting(
            message_bytes, block_bytes, config) != Wirehair_Success ||
        !CheckExplicitDiagnostics(
            explicit_encoder, explicit_decoder, config,
            K, block_bytes, "canonical"))
    {
        return false;
    }
    const uint64_t columns = (uint64_t)K + config.Params.Staircase +
        config.Params.DenseRows + config.Params.HeavyRows;
    const size_t intermediate_bytes = (size_t)(columns * block_bytes);
    if (std::memcmp(
            normal.IntermediateBlocks(),
            explicit_encoder.IntermediateBlocks(),
            intermediate_bytes) != 0)
    {
        std::fprintf(stderr,
            "explicit canonical: intermediate vector differs from normal\n");
        return false;
    }
    for (uint32_t id = 0u; id < K + 24u; ++id)
    {
        std::vector<uint8_t> expected(block_bytes, 0xa5u);
        std::vector<uint8_t> actual(block_bytes, 0x5au);
        uint32_t expected_bytes = UINT32_MAX;
        uint32_t actual_bytes = UINT32_MAX;
        uint64_t expected_ops = UINT64_MAX;
        uint64_t actual_ops = UINT64_MAX;
        if (normal.EncodeResult(
                id, expected.data(), block_bytes,
                &expected_bytes, &expected_ops) != Wirehair_Success ||
            explicit_encoder.EncodeResult(
                id, actual.data(), block_bytes,
                &actual_bytes, &actual_ops) != Wirehair_Success ||
            expected_bytes != actual_bytes ||
            expected_ops != actual_ops ||
            std::memcmp(
                expected.data(), actual.data(), expected_bytes) != 0)
        {
            std::fprintf(stderr,
                "explicit canonical: packet %u differs from normal\n", id);
            return false;
        }
    }

    wirehair_v2::MessagePrecodeEncoder replay_encoder;
    wirehair_v2::MessagePrecodeDecoder replay_decoder;
    if (replay_encoder.InitializeResult(
            message.data(), message_bytes, block_bytes,
            &explicit_encoder.Profile(), &explicit_encoder.Options()) !=
            Wirehair_InvalidInput ||
        replay_decoder.InitializeResult(
            message_bytes, block_bytes,
            &explicit_decoder.Profile(), &explicit_decoder.Options()) !=
            Wirehair_InvalidInput)
    {
        std::fprintf(stderr,
            "explicit canonical: diagnostic profile replay did not fail closed\n");
        return false;
    }

    // Exercise all four symmetric cache receipts and replacement transitions
    // on the same live objects.
    wirehair_v2::MessagePrecodeEncoder cache_encoder;
    wirehair_v2::MessagePrecodeDecoder cache_decoder;
    const uint32_t cache_masks[] = {0u, 1u, 2u, 3u, 0u};
    for (uint32_t mask : cache_masks)
    {
        config.CacheSystematicSource = (mask & 1u) != 0u;
        config.CacheReceivedSystematicPackets = (mask & 2u) != 0u;
        if (cache_encoder.InitializeExplicitResultForTesting(
                message.data(), message_bytes, block_bytes, config) !=
                Wirehair_Success ||
            cache_decoder.InitializeExplicitResultForTesting(
                message_bytes, block_bytes, config) != Wirehair_Success ||
            !CheckExplicitDiagnostics(
                cache_encoder, cache_decoder, config,
                K, block_bytes, "cache replacement"))
        {
            return false;
        }
        std::vector<uint8_t> source(block_bytes, 0u);
        uint32_t source_bytes = 0u;
        if (cache_encoder.EncodeResult(
                K - 1u, source.data(), block_bytes, &source_bytes) !=
                Wirehair_Success ||
            source_bytes != tail_bytes ||
            std::memcmp(
                source.data(),
                message.data() + (size_t)(K - 1u) * block_bytes,
                tail_bytes) != 0)
        {
            std::fprintf(stderr,
                "explicit cache replacement: systematic source mismatch\n");
            return false;
        }
        if (cache_decoder.DecodeResult(
                K - 1u, source.data(), source_bytes) != Wirehair_NeedMore ||
            cache_decoder.CachedSystematicPacketCount() !=
                (config.CacheReceivedSystematicPackets ? 1u : 0u) ||
            !CompleteDirectRepairOnly(
                cache_encoder,
                cache_decoder,
                message,
                K,
                block_bytes,
                "cache replacement"))
        {
            std::fprintf(stderr,
                "explicit cache replacement: roundtrip/cache count failed\n");
            return false;
        }
    }

    // Invalid dimensions, mixed odd widths, and equation-active overrides are
    // rejected before replacing either endpoint.
    const wirehair_v2::SeedProfile preserved_encoder_profile =
        cache_encoder.Profile();
    const wirehair_v2::SeedProfile preserved_decoder_profile =
        cache_decoder.Profile();
    wirehair_v2::ExplicitMessagePrecodeConfigForTesting invalid = config;
    ++invalid.Params.BlockCount;
    if (cache_encoder.InitializeExplicitResultForTesting(
            message.data(), message_bytes, block_bytes, invalid) !=
            Wirehair_InvalidInput ||
        cache_decoder.InitializeExplicitResultForTesting(
            message_bytes, block_bytes, invalid) != Wirehair_InvalidInput ||
        cache_encoder.Profile().V2PrecodeSeed !=
            preserved_encoder_profile.V2PrecodeSeed ||
        cache_decoder.Profile().V2PrecodeSeed !=
            preserved_decoder_profile.V2PrecodeSeed)
    {
        std::fprintf(stderr,
            "explicit invalid dimensions were accepted or changed state\n");
        return false;
    }
    invalid = config;
    invalid.Packet.MixCount = 0u;
    if (cache_encoder.InitializeExplicitResultForTesting(
            message.data(), message_bytes, block_bytes, invalid) !=
            Wirehair_InvalidInput ||
        cache_decoder.InitializeExplicitResultForTesting(
            message_bytes, block_bytes, invalid) != Wirehair_InvalidInput)
    {
        std::fprintf(stderr,
            "explicit invalid packet-row configuration was accepted\n");
        return false;
    }
    const uint32_t odd_block_bytes = block_bytes - 1u;
    const uint64_t odd_message_bytes =
        (uint64_t)(K - 1u) * odd_block_bytes + tail_bytes;
    if (cache_encoder.InitializeExplicitResultForTesting(
            message.data(), odd_message_bytes, odd_block_bytes, config) !=
            Wirehair_InvalidInput ||
        cache_decoder.InitializeExplicitResultForTesting(
            odd_message_bytes, odd_block_bytes, config) !=
            Wirehair_InvalidInput)
    {
        std::fprintf(stderr, "explicit mixed odd width was accepted\n");
        return false;
    }
    wirehair_v2::SetStaircaseDegreeScaleForTesting(3.25);
    const WirehairResult overridden_encoder =
        cache_encoder.InitializeExplicitResultForTesting(
            message.data(), message_bytes, block_bytes, config);
    const WirehairResult overridden_decoder =
        cache_decoder.InitializeExplicitResultForTesting(
            message_bytes, block_bytes, config);
    wirehair_v2::ClearStaircaseDegreeScaleForTesting();
    if (overridden_encoder != Wirehair_InvalidInput ||
        overridden_decoder != Wirehair_InvalidInput)
    {
        std::fprintf(stderr,
            "explicit effective-parameter override was not rejected\n");
        return false;
    }

    // Exhaust every encoder-guarded allocation, checking that the complete
    // live object remains byte-for-byte operational after each OOM.
    config.CacheSystematicSource = true;
    config.CacheReceivedSystematicPackets = true;
    wirehair_v2::MessagePrecodeEncoder oom_encoder;
    if (oom_encoder.InitializeExplicitResultForTesting(
            message.data(), message_bytes, block_bytes, config) !=
        Wirehair_Success)
    {
        return false;
    }
    const std::vector<uint8_t> encoder_intermediate(
        oom_encoder.IntermediateBlocks(),
        oom_encoder.IntermediateBlocks() + intermediate_bytes);
    const wirehair_v2::SeedProfile encoder_profile = oom_encoder.Profile();
    const wirehair_v2::MessagePrecodeEncoderOptions encoder_options =
        oom_encoder.Options();
    unsigned encoder_ooms = 0u;
    bool encoder_succeeded = false;
    for (int64_t countdown = 0; countdown < 16; ++countdown)
    {
        wirehair_v2::SetAllocationFailureCountdownForTesting(countdown);
        const WirehairResult result =
            oom_encoder.InitializeExplicitResultForTesting(
                message.data(), message_bytes, block_bytes, config);
        wirehair_v2::SetAllocationFailureCountdownForTesting(-1);
        if (result == Wirehair_OOM)
        {
            ++encoder_ooms;
            std::vector<uint8_t> preserved_source(block_bytes, 0u);
            uint32_t preserved_source_bytes = 0u;
            if (oom_encoder.Profile().V2PrecodeSeed !=
                    encoder_profile.V2PrecodeSeed ||
                oom_encoder.Profile().BlockCount !=
                    encoder_profile.BlockCount ||
                oom_encoder.Profile().BlockBytes !=
                    encoder_profile.BlockBytes ||
                oom_encoder.Profile().V2PrecodeContractVersion !=
                    encoder_profile.V2PrecodeContractVersion ||
                oom_encoder.Profile().V2PacketRowContractVersion !=
                    encoder_profile.V2PacketRowContractVersion ||
                oom_encoder.Profile().V2StaircaseCount !=
                    encoder_profile.V2StaircaseCount ||
                oom_encoder.Profile().V2DenseRowCount !=
                    encoder_profile.V2DenseRowCount ||
                oom_encoder.Profile().V2HeavyRowCount !=
                    encoder_profile.V2HeavyRowCount ||
                !SameOptions(oom_encoder.Options(), encoder_options) ||
                !oom_encoder.HasSystematicSourceCache() ||
                oom_encoder.SystematicSourceCacheBytes() != message.size() ||
                std::memcmp(
                    oom_encoder.IntermediateBlocks(),
                    encoder_intermediate.data(),
                    intermediate_bytes) != 0 ||
                oom_encoder.EncodeResult(
                    K - 1u,
                    preserved_source.data(),
                    block_bytes,
                    &preserved_source_bytes) != Wirehair_Success ||
                preserved_source_bytes != tail_bytes ||
                std::memcmp(
                    preserved_source.data(),
                    message.data() + (size_t)(K - 1u) * block_bytes,
                    tail_bytes) != 0)
            {
                std::fprintf(stderr,
                    "explicit encoder OOM did not preserve live state\n");
                return false;
            }
            continue;
        }
        if (result != Wirehair_Success || encoder_ooms != 5u) {
            std::fprintf(stderr,
                "explicit encoder OOM sweep ended result=%d failures=%u\n",
                (int)result, encoder_ooms);
            return false;
        }
        encoder_succeeded = true;
        break;
    }
    if (!encoder_succeeded) {
        std::fprintf(stderr,
            "explicit encoder OOM sweep never reached success\n");
        return false;
    }

    // Decoder OOM sweeps start from a partially fed cached decoder.  Duplicate
    // validation after every failed reinit proves its packet table, payload,
    // progress counters, and cache remained coherent.
    wirehair_v2::MessagePrecodeEncoder feed_encoder;
    if (feed_encoder.InitializeExplicitResultForTesting(
            message.data(), message_bytes, block_bytes, config) !=
        Wirehair_Success)
    {
        return false;
    }
    std::vector<uint8_t> packet(block_bytes, 0u);
    uint32_t packet_bytes = 0u;
    if (feed_encoder.EncodeResult(
            0u, packet.data(), block_bytes, &packet_bytes) !=
        Wirehair_Success)
    {
        return false;
    }
    unsigned decoder_ooms = 0u;
    bool decoder_succeeded = false;
    for (int64_t countdown = 0; countdown < 16; ++countdown)
    {
        wirehair_v2::MessagePrecodeDecoder decoder;
        if (decoder.InitializeExplicitResultForTesting(
                message_bytes, block_bytes, config) != Wirehair_Success ||
            decoder.DecodeResult(0u, packet.data(), packet_bytes) !=
                Wirehair_NeedMore)
        {
            return false;
        }
        const wirehair_v2::SeedProfile decoder_profile = decoder.Profile();
        const wirehair_v2::MessagePrecodeEncoderOptions decoder_options =
            decoder.Options();
        wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(countdown);
        const WirehairResult result =
            decoder.InitializeExplicitResultForTesting(
                message_bytes, block_bytes, config);
        wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(-1);
        if (result == Wirehair_OOM)
        {
            ++decoder_ooms;
            if (decoder.Profile().V2PrecodeSeed !=
                    decoder_profile.V2PrecodeSeed ||
                decoder.Profile().BlockCount !=
                    decoder_profile.BlockCount ||
                decoder.Profile().BlockBytes !=
                    decoder_profile.BlockBytes ||
                decoder.Profile().V2PrecodeContractVersion !=
                    decoder_profile.V2PrecodeContractVersion ||
                decoder.Profile().V2PacketRowContractVersion !=
                    decoder_profile.V2PacketRowContractVersion ||
                decoder.Profile().V2StaircaseCount !=
                    decoder_profile.V2StaircaseCount ||
                decoder.Profile().V2DenseRowCount !=
                    decoder_profile.V2DenseRowCount ||
                decoder.Profile().V2HeavyRowCount !=
                    decoder_profile.V2HeavyRowCount ||
                !SameOptions(decoder.Options(), decoder_options) ||
                decoder.ReceivedCount() != 1u ||
                decoder.SolveAttemptCount() != 0u ||
                decoder.CachedSystematicPacketCount() != 1u ||
                !decoder.HasSystematicPacketCache() ||
                !SamePrecodeParams(decoder.System().Params, config.Params) ||
                decoder.DecodeResult(
                    0u, packet.data(), packet_bytes) != Wirehair_NeedMore ||
                decoder.ReceivedCount() != 1u)
            {
                std::fprintf(stderr,
                    "explicit decoder OOM did not preserve partial state\n");
                return false;
            }
            continue;
        }
        if (result != Wirehair_Success || decoder_ooms != 6u) {
            std::fprintf(stderr,
                "explicit decoder OOM sweep ended result=%d failures=%u\n",
                (int)result, decoder_ooms);
            return false;
        }
        decoder_succeeded = true;
        break;
    }
    if (!decoder_succeeded) {
        std::fprintf(stderr,
            "explicit decoder OOM sweep never reached success\n");
        return false;
    }

    wirehair_v2::ExplicitMessagePrecodeConfigForTesting pure;
    pure.Params = wirehair_v2::MakeCertifiedParams(
        37u, UINT64_C(0x91f4b3c28765a10d));
    pure.Params.Staircase = 11u;
    pure.Params.DenseRows = 6u;
    pure.Params.HeavyRows = 9u;
    pure.Params.SourceHits = 3u;
    pure.Params.DegreeBalancedStaircase = true;
    pure.Params.StaircaseDegreeScale = 7.25;
    pure.Params.DenseIdentityCorner = true;
    pure.Packet.PeelSeed = UINT32_C(0x6bd073a1);
    pure.Packet.MixCount = 3u;
    if (!pure.PinActiveEquationStateForTesting() ||
        !RunOneExplicitConfiguration(
            pure, 15u, 4u, "noncanonical GF256"))
    {
        return false;
    }

    wirehair_v2::ExplicitMessagePrecodeConfigForTesting tracking;
    const bool tracking_rows_configured =
        wirehair_v2::SetMixedGF16RowsForTesting(0u) &&
        wirehair_v2::SetMixedGF256RowsForTesting(8u);
    tracking.Params = wirehair_v2::MakeMixedParams(
        43u, UINT64_C(0xb1716d94ac2058e3));
    tracking.Params.Staircase = 9u;
    tracking.Params.DenseRows = 4u;
    tracking.Params.SourceHits = 2u;
    tracking.Params.DegreeBalancedStaircase = true;
    tracking.Params.StaircaseDegreeScale = 8.5;
    tracking.Packet.PeelSeed = UINT32_C(0x41f06c33);
    tracking.Packet.MixCount = 3u;
    wirehair_v2::SetMixedBandTrackingXForTesting(true);
    const bool tracking_ok = tracking_rows_configured &&
        tracking.Params.HeavyRows == 8u &&
        wirehair_v2::ActiveMixedCoefficientPeriod() == 244u &&
        tracking.PinActiveEquationStateForTesting() &&
        RunOneExplicitConfiguration(
            tracking, 18u, 7u, "noncanonical mixed tracking-X");
    wirehair_v2::ClearMixedBandTrackingXForTesting();
    const bool tracking_rows_restored =
        wirehair_v2::SetMixedGF256RowsForTesting(
            wirehair_v2::kMixedGF256Rows) &&
        wirehair_v2::SetMixedGF16RowsForTesting(
            wirehair_v2::kMixedGF16Rows);
    if (!tracking_ok || !tracking_rows_restored ||
        !RunExplicitTrackingXConcurrencyCase() ||
        !RunExplicitEquationStatePinningCases())
    {
        return false;
    }

    std::printf(
        "explicit exact configuration, pinned equation state, cache "
        "replacement, and OOM: PASS\n");
    return true;
}

bool RunArchitectureBoundaryCases()
{
    const uint32_t block_bytes = 2u;
    std::vector<uint32_t> block_counts;
    for (uint32_t K = 2u; K <= 100u; ++K) {
        block_counts.push_back(K);
    }
    const uint32_t boundary_counts[] = {
        101u, 1024u, 4095u, 4096u, 9999u, 10000u, 64000u
    };
    block_counts.insert(
        block_counts.end(),
        boundary_counts,
        boundary_counts +
            sizeof(boundary_counts) / sizeof(boundary_counts[0]));

    wirehair_v2::MessagePrecodeEncoderOptions legacy_options;
    legacy_options.Completion =
        wirehair_v2::CompletionField::MixedGF256GF16;
    wirehair_v2::MessagePrecodeEncoderOptions small_band_options =
        legacy_options;
    small_band_options.Architecture =
        wirehair_v2::V2PrecodeArchitecture::SmallBandD4;

    for (uint32_t K : block_counts)
    {
        const uint64_t message_bytes = (uint64_t)K * block_bytes;
        std::vector<uint8_t> message((size_t)message_bytes);
        for (size_t i = 0; i < message.size(); ++i) {
            message[i] = (uint8_t)(i * 53u + K * 11u + 7u);
        }
        const wirehair_v2::SeedProfile base =
            wirehair_v2::SelectSeedProfile(K, block_bytes);
        const uint32_t inherited_staircase = wirehair::GetDenseCount(K);
        const uint32_t small_band_staircase =
            wirehair_v2::SmallBandStaircaseCount(K);
        if (base.DenseCount != inherited_staircase)
        {
            std::fprintf(stderr,
                "architecture boundary K=%u: base DenseCount=%u expected=%u\n",
                K, base.DenseCount, inherited_staircase);
            return false;
        }

        wirehair_v2::MessagePrecodeEncoder legacy_encoder;
        wirehair_v2::MessagePrecodeEncoder small_band_encoder;
        if (legacy_encoder.InitializeResult(
                message.data(), message_bytes, block_bytes,
                &base, &legacy_options) != Wirehair_Success ||
            small_band_encoder.InitializeResult(
                message.data(), message_bytes, block_bytes,
                &base, &small_band_options) != Wirehair_Success)
        {
            std::fprintf(stderr,
                "architecture boundary K=%u: encoder init failed\n", K);
            return false;
        }

        const wirehair_v2::SeedProfile& legacy =
            legacy_encoder.Profile();
        const wirehair_v2::SeedProfile& small_band =
            small_band_encoder.Profile();
        if (legacy.DenseCount != inherited_staircase ||
            small_band.DenseCount != inherited_staircase ||
            legacy.DenseCount != base.DenseCount ||
            small_band.DenseCount != base.DenseCount ||
            legacy.DenseSeed != base.DenseSeed ||
            small_band.DenseSeed != base.DenseSeed ||
            legacy.V2Architecture !=
                wirehair_v2::V2PrecodeArchitecture::LegacyD12 ||
            small_band.V2Architecture !=
                wirehair_v2::V2PrecodeArchitecture::SmallBandD4 ||
            legacy.V2SeedPolicy !=
                wirehair_v2::V2SeedDerivation::ProfileDerived ||
            small_band.V2SeedPolicy !=
                wirehair_v2::V2SeedDerivation::ProfileDerived ||
            legacy.V2StaircaseCount != inherited_staircase ||
            small_band.V2StaircaseCount != small_band_staircase ||
            legacy.V2DenseRowCount != 12u ||
            small_band.V2DenseRowCount != 4u ||
            legacy.V2PrecodeContractVersion !=
                wirehair_v2::PrecodeContractVersion(
                    wirehair_v2::CompletionField::MixedGF256GF16,
                    false,
                    wirehair_v2::V2PrecodeArchitecture::LegacyD12) ||
            small_band.V2PrecodeContractVersion !=
                wirehair_v2::PrecodeContractVersion(
                    wirehair_v2::CompletionField::MixedGF256GF16,
                    false,
                    wirehair_v2::V2PrecodeArchitecture::SmallBandD4))
        {
            std::fprintf(stderr,
                "architecture boundary K=%u: selected profile mismatch "
                "legacy(S=%u,D2=%u) small-band(S=%u,D2=%u)\n",
                K,
                legacy.V2StaircaseCount,
                legacy.V2DenseRowCount,
                small_band.V2StaircaseCount,
                small_band.V2DenseRowCount);
            return false;
        }

        wirehair_v2::MessagePrecodeDecoder legacy_decoder;
        wirehair_v2::MessagePrecodeDecoder small_band_decoder;
        if (legacy_decoder.InitializeResult(
                message_bytes, block_bytes, &legacy) != Wirehair_Success ||
            small_band_decoder.InitializeResult(
                message_bytes, block_bytes, &small_band) != Wirehair_Success ||
            !SameOptions(legacy_decoder.Options(), legacy_options) ||
            !SameOptions(small_band_decoder.Options(), small_band_options) ||
            !SamePrecodeParams(
                legacy_encoder.BlockEncoder().System().Params,
                legacy_decoder.System().Params) ||
            !SamePrecodeParams(
                small_band_encoder.BlockEncoder().System().Params,
                small_band_decoder.System().Params) ||
            legacy_decoder.PacketSeedAttempt() != legacy.V2SeedAttempt ||
            small_band_decoder.PacketSeedAttempt() !=
                small_band.V2SeedAttempt ||
            legacy_decoder.PacketPeelSeed() != legacy.V2PacketPeelSeed ||
            small_band_decoder.PacketPeelSeed() !=
                small_band.V2PacketPeelSeed)
        {
            std::fprintf(stderr,
                "architecture boundary K=%u: decoder reconstruction "
                "mismatch\n", K);
            return false;
        }

        if (K == 64u)
        {
            const auto rejects = [&](wirehair_v2::SeedProfile malformed,
                                     const char* field) {
                wirehair_v2::MessagePrecodeDecoder decoder;
                if (decoder.InitializeResult(
                        message_bytes, block_bytes, &malformed) ==
                    Wirehair_InvalidInput)
                {
                    return true;
                }
                std::fprintf(stderr,
                    "architecture boundary: malformed %s was accepted\n",
                    field);
                return false;
            };

            wirehair_v2::SeedProfile malformed = small_band;
            malformed.V2Architecture =
                wirehair_v2::V2PrecodeArchitecture::LegacyD12;
            if (!rejects(malformed, "architecture/version pair")) {
                return false;
            }
            malformed = small_band;
            malformed.V2PrecodeContractVersion =
                wirehair_v2::kMixedPrecodeContractVersion;
            if (!rejects(malformed, "small-band contract version")) {
                return false;
            }
            malformed = small_band;
            ++malformed.V2StaircaseCount;
            if (!rejects(malformed, "small-band staircase count")) {
                return false;
            }
            malformed = small_band;
            malformed.V2DenseRowCount = 12u;
            if (!rejects(malformed, "small-band dense-row count")) {
                return false;
            }
            malformed = legacy;
            malformed.V2SeedPolicy =
                wirehair_v2::V2SeedDerivation::RawUniform;
            if (!rejects(malformed, "legacy raw-seed policy")) {
                return false;
            }
            malformed = small_band;
            malformed.V2SeedPolicy =
                wirehair_v2::V2SeedDerivation::RawUniform;
            malformed.V2SeedAttempt = 1u;
            if (!rejects(malformed, "raw-seed retry attempt")) {
                return false;
            }
            malformed = small_band;
            malformed.V2SeedPolicy =
                static_cast<wirehair_v2::V2SeedDerivation>(UINT32_MAX);
            if (!rejects(malformed, "unknown seed policy")) {
                return false;
            }

            // RawUniform is a complete table-free overlay: its direct seeds,
            // architecture, and attempt zero are sufficient even when every
            // legacy per-K seed-table field is absent.
            wirehair_v2::SeedProfile raw;
            raw.BlockCount = K;
            raw.BlockBytes = block_bytes;
            raw.V2SeedPolicy =
                wirehair_v2::V2SeedDerivation::RawUniform;
            wirehair_v2::PrecodeSystem raw_system;
            raw_system.Params =
                wirehair_v2::MakeMixedParams(K, UINT64_C(0x1234));
            raw_system.Params.Staircase =
                wirehair_v2::SmallBandStaircaseCount(K);
            raw_system.Params.DenseRows = 4u;
            wirehair_v2::PacketRowConfig raw_packet;
            raw_packet.PeelSeed = 0x1234u;
            raw_packet.MixCount =
                wirehair_v2::kCertifiedPacketMixCount;
            wirehair_v2::MessagePrecodeEncoderOptions raw_options =
                small_band_options;
            raw_options.PrecodeSeedSalt = 0u;
            raw_options.RecoveryRowSeedSalt = 0u;
            wirehair_v2::BindMessagePrecodeProfile(
                raw, raw_options, raw_system, raw_packet, 0u);
            wirehair_v2::MessagePrecodeEncoder raw_encoder;
            wirehair_v2::MessagePrecodeDecoder raw_decoder;
            if (raw.DenseCount != 0u || raw.PeelSeed != 0u ||
                raw.DenseSeed != 0u || raw.UsedPeelFixup ||
                raw.UsedDenseFixup ||
                raw_encoder.InitializeResult(
                    message.data(), message_bytes, block_bytes,
                    &raw, &raw_options) != Wirehair_Success ||
                raw_decoder.InitializeResult(
                    message_bytes, block_bytes, &raw_encoder.Profile(),
                    &raw_options) != Wirehair_Success ||
                !SamePrecodeParams(
                    raw_encoder.BlockEncoder().System().Params,
                    raw_decoder.System().Params) ||
                raw_encoder.Profile().V2SeedAttempt != 0u ||
                raw_decoder.PacketSeedAttempt() != 0u ||
                !CompleteDirectRepairOnly(
                    raw_encoder, raw_decoder, message, K, block_bytes,
                    "raw dispatch-v1"))
            {
                std::fprintf(stderr,
                    "architecture boundary: table-free raw profile "
                    "did not replay and recover exactly\n");
                return false;
            }

            malformed = raw;
            malformed.V2PacketPeelSeed ^= 1u;
            if (!rejects(malformed, "raw packet-seed fold")) {
                return false;
            }
            malformed = raw;
            malformed.DenseCount = 1u;
            if (!rejects(malformed, "raw dense-count table state")) {
                return false;
            }
            malformed = raw;
            malformed.PeelSeed = 1u;
            if (!rejects(malformed, "raw peel-seed table state")) {
                return false;
            }
            malformed = raw;
            malformed.DenseSeed = 1u;
            if (!rejects(malformed, "raw dense-seed table state")) {
                return false;
            }
            malformed = raw;
            malformed.PeelSeedBucket = 1u;
            if (!rejects(malformed, "raw peel-seed bucket state")) {
                return false;
            }
            malformed = raw;
            malformed.UsedPeelFixup = true;
            if (!rejects(malformed, "raw peel fixup")) {
                return false;
            }
            malformed = raw;
            malformed.UsedDenseFixup = true;
            if (!rejects(malformed, "raw dense fixup")) {
                return false;
            }
            malformed = raw;
            malformed.Tuned = true;
            if (!rejects(malformed, "raw tuned seed state")) {
                return false;
            }
            malformed = raw;
            malformed.Policy.Codec.MinDegree = 37u;
            if (!rejects(malformed, "raw legacy peel policy")) {
                return false;
            }
            malformed = raw;
            malformed.TuningResidualMean =
                std::numeric_limits<double>::quiet_NaN();
            if (!rejects(malformed, "raw tuning residual")) {
                return false;
            }
            malformed = raw;
            malformed.TuningCandidatesRequested = UINT16_MAX;
            if (!rejects(malformed, "raw tuning candidate diagnostics")) {
                return false;
            }
            malformed = raw;
            malformed.V2PrecodeSeedSalt = UINT64_C(1);
            if (!rejects(malformed, "raw precode salt")) {
                return false;
            }
            malformed = raw;
            malformed.V2RecoveryRowSeedSalt = UINT64_C(1);
            if (!rejects(malformed, "raw packet salt")) {
                return false;
            }
        }
    }
    return true;
}

bool RunOptionContractCase(
    const char* label,
    const wirehair_v2::MessagePrecodeEncoderOptions& options,
    const wirehair_v2::MessagePrecodeEncoderOptions& mismatched,
    uint32_t K = 64u)
{
    const uint32_t block_bytes =
        options.Completion ==
            wirehair_v2::CompletionField::MixedGF256GF16 ? 16u : 17u;
    const uint64_t message_bytes = (uint64_t)K * block_bytes;
    std::vector<uint8_t> message((size_t)message_bytes);
    for (size_t i = 0; i < message.size(); ++i) {
        message[i] = (uint8_t)(i * 41u + 23u);
    }

    wirehair_v2::Codec encoder;
    if (encoder.InitializePrecodeEncoder(
            message.data(), message_bytes, block_bytes, nullptr, &options) !=
        Wirehair_Success)
    {
        std::fprintf(stderr, "contract %s: encoder init failed\n", label);
        return false;
    }
    const wirehair_v2::SeedProfile profile = encoder.Profile();
    if (wirehair_v2::kPacketRowContractVersion != 4u ||
        !profile.V2SeedSelected ||
        profile.V2PrecodeContractVersion !=
            wirehair_v2::PrecodeContractVersion(
                options.Completion, options.AdaptiveDenseTwoAnchor,
                options.Architecture) ||
        profile.V2PacketRowContractVersion !=
            wirehair_v2::kPacketRowContractVersion ||
        profile.V2Architecture != options.Architecture ||
        profile.V2SeedPolicy !=
            wirehair_v2::V2SeedDerivation::ProfileDerived ||
        profile.V2StaircaseCount !=
            (options.Architecture ==
                    wirehair_v2::V2PrecodeArchitecture::SmallBandD4 ?
                wirehair_v2::SmallBandStaircaseCount(K) :
                profile.DenseCount) ||
        profile.V2DenseRowCount !=
            (options.Architecture ==
                    wirehair_v2::V2PrecodeArchitecture::SmallBandD4 ?
                4u : 12u) ||
        profile.V2RecoveryMixCount != options.RecoveryMixCount ||
        profile.V2DenseIdentityCorner != options.DenseIdentityCorner ||
        profile.V2DenseTwoAnchor !=
            (options.AdaptiveDenseTwoAnchor &&
             K >= wirehair_v2::kDenseTwoAnchorMinBlockCount) ||
        profile.V2AdaptiveDenseTwoAnchor !=
            options.AdaptiveDenseTwoAnchor ||
        profile.V2CompletionField != options.Completion ||
        profile.V2PrecodeSeedSalt != options.PrecodeSeedSalt ||
        profile.V2RecoveryRowSeedSalt != options.RecoveryRowSeedSalt)
    {
        std::fprintf(stderr, "contract %s: published profile mismatch\n", label);
        return false;
    }

    wirehair_v2::MessagePrecodeDecoder inherited;
    if (inherited.InitializeResult(message_bytes, block_bytes, &profile) !=
            Wirehair_Success ||
        !SameOptions(inherited.Options(), options) ||
        inherited.System().Params.Staircase != profile.V2StaircaseCount ||
        inherited.System().Params.DenseRows != profile.V2DenseRowCount ||
        inherited.System().Params.HeavyRows != profile.V2HeavyRowCount ||
        inherited.System().Params.Field != profile.V2CompletionField ||
        inherited.System().Params.SourceHits != profile.V2SourceHits ||
        inherited.System().Params.Seed != profile.V2PrecodeSeed ||
        inherited.System().Params.DenseTwoAnchor !=
            profile.V2DenseTwoAnchor ||
        inherited.PacketPeelSeed() != profile.V2PacketPeelSeed ||
        !CompleteDirectRepairOnly(
            encoder, inherited, message, K, block_bytes, label))
    {
        std::fprintf(stderr, "contract %s: inherited control failed\n", label);
        return false;
    }

    wirehair_v2::MessagePrecodeDecoder explicit_match;
    if (explicit_match.InitializeResult(
            message_bytes, block_bytes, &profile, &options) !=
            Wirehair_Success ||
        !CompleteDirectRepairOnly(
            encoder, explicit_match, message, K, block_bytes, label))
    {
        std::fprintf(stderr, "contract %s: explicit control failed\n", label);
        return false;
    }

    wirehair_v2::MessagePrecodeDecoder rejected;
    if (rejected.InitializeResult(
            message_bytes, block_bytes, &profile, &mismatched) !=
        Wirehair_InvalidInput)
    {
        std::fprintf(stderr,
            "contract %s: explicit mismatch was not rejected\n", label);
        return false;
    }
    return true;
}

bool RunBoundContractCases()
{
    wirehair_v2::MessagePrecodeEncoderOptions defaults;
    wirehair_v2::MessagePrecodeEncoderOptions variant = defaults;
    variant.RecoveryMixCount = 2u;
    const uint32_t invalid_K = 64u;
    const uint32_t invalid_block_bytes = 17u;
    const uint64_t invalid_message_bytes =
        (uint64_t)invalid_K * invalid_block_bytes;
    std::vector<uint8_t> invalid_message(
        (size_t)invalid_message_bytes, 0x6bu);
    wirehair_v2::Codec invalid_encoder;
    wirehair_v2::MessagePrecodeDecoder invalid_decoder;
    if (invalid_encoder.InitializePrecodeEncoder(
            invalid_message.data(),
            invalid_message_bytes,
            invalid_block_bytes,
            nullptr,
            &variant) != Wirehair_InvalidInput ||
        invalid_decoder.InitializeResult(
            invalid_message_bytes,
            invalid_block_bytes,
            nullptr,
            &variant) != Wirehair_InvalidInput)
    {
        std::fprintf(stderr,
            "contract: GF256/mix2 option pair was accepted\n");
        return false;
    }

    const uint32_t mixed_invalid_block_bytes = 16u;
    const uint64_t mixed_invalid_message_bytes =
        (uint64_t)invalid_K * mixed_invalid_block_bytes;
    std::vector<uint8_t> mixed_invalid_message(
        (size_t)mixed_invalid_message_bytes, 0x5du);
    const auto reject_options = [&](const char* label,
        const wirehair_v2::MessagePrecodeEncoderOptions& rejected_options)
    {
        wirehair_v2::Codec rejected_encoder;
        wirehair_v2::MessagePrecodeDecoder rejected_decoder;
        if (rejected_encoder.InitializePrecodeEncoder(
                mixed_invalid_message.data(), mixed_invalid_message_bytes,
                mixed_invalid_block_bytes, nullptr, &rejected_options) !=
                    Wirehair_InvalidInput ||
            rejected_decoder.InitializeResult(
                mixed_invalid_message_bytes, mixed_invalid_block_bytes,
                nullptr, &rejected_options) != Wirehair_InvalidInput)
        {
            std::fprintf(stderr,
                "contract: invalid option pair %s was accepted\n",
                label);
            return false;
        }
        return true;
    };
    wirehair_v2::MessagePrecodeEncoderOptions invalid_two_anchor = defaults;
    invalid_two_anchor.AdaptiveDenseTwoAnchor = true;
    if (!reject_options("GF256/mix3/two-anchor", invalid_two_anchor)) {
        return false;
    }
    invalid_two_anchor.Completion =
        wirehair_v2::CompletionField::MixedGF256GF16;
    if (!reject_options("mixed/mix3/two-anchor", invalid_two_anchor)) {
        return false;
    }
    invalid_two_anchor.RecoveryMixCount = 2u;
    invalid_two_anchor.DenseIdentityCorner = true;
    if (!reject_options(
            "mixed/mix2/identity-corner", invalid_two_anchor))
    {
        return false;
    }
    variant = defaults;
    variant.DenseIdentityCorner = true;
    if (!RunOptionContractCase("dense-corner", variant, defaults)) {
        return false;
    }
    variant = defaults;
    variant.PrecodeSeedSalt ^= UINT64_C(0x123456789abcdef0);
    if (!RunOptionContractCase("precode-salt", variant, defaults)) {
        return false;
    }
    variant = defaults;
    variant.RecoveryRowSeedSalt ^= UINT64_C(0xfedcba9876543210);
    if (!RunOptionContractCase("packet-salt", variant, defaults)) {
        return false;
    }
    variant = defaults;
    variant.Completion =
        wirehair_v2::CompletionField::MixedGF256GF16;
    if (!RunOptionContractCase("mixed-completion", variant, defaults)) {
        return false;
    }
    wirehair_v2::MessagePrecodeEncoderOptions mixed_mix3 = variant;
    variant.RecoveryMixCount = 2u;
    if (!RunOptionContractCase("mixed-mix2", variant, mixed_mix3)) {
        return false;
    }
    mixed_mix3 = variant;
    variant.AdaptiveDenseTwoAnchor = true;
    if (!RunOptionContractCase(
            "mixed-mix2-two-anchor", variant, mixed_mix3))
    {
        return false;
    }
    if (!RunOptionContractCase(
            "mixed-mix2-two-anchor-active", variant, mixed_mix3,
            wirehair_v2::kDenseTwoAnchorMinBlockCount))
    {
        return false;
    }

    wirehair_v2::MessagePrecodeEncoderOptions legacy_mixed = defaults;
    legacy_mixed.Completion =
        wirehair_v2::CompletionField::MixedGF256GF16;
    wirehair_v2::MessagePrecodeEncoderOptions small_band = legacy_mixed;
    small_band.Architecture =
        wirehair_v2::V2PrecodeArchitecture::SmallBandD4;
    if (!RunOptionContractCase(
            "small-band-d4", small_band, legacy_mixed, 64u) ||
        !RunOptionContractCase(
            "small-band-d4-above-band", small_band, legacy_mixed, 101u))
    {
        return false;
    }
    wirehair_v2::MessagePrecodeEncoderOptions invalid_architecture = defaults;
    invalid_architecture.Architecture =
        static_cast<wirehair_v2::V2PrecodeArchitecture>(UINT32_MAX);
    if (!reject_options("unknown architecture", invalid_architecture)) {
        return false;
    }
    wirehair_v2::MessagePrecodeEncoderOptions invalid_small_band = small_band;
    invalid_small_band.Completion = wirehair_v2::CompletionField::GF256;
    if (!reject_options("small-band/GF256", invalid_small_band)) {
        return false;
    }
    invalid_small_band = small_band;
    invalid_small_band.RecoveryMixCount = 2u;
    if (!reject_options("small-band/mix2", invalid_small_band)) {
        return false;
    }
    invalid_small_band = small_band;
    invalid_small_band.AdaptiveDenseTwoAnchor = true;
    if (!reject_options("small-band/two-anchor", invalid_small_band)) {
        return false;
    }
    invalid_small_band = small_band;
    invalid_small_band.DenseIdentityCorner = true;
    if (!reject_options("small-band/identity-corner", invalid_small_band)) {
        return false;
    }

    const uint32_t K = 320u;
    const uint32_t block_bytes = 7u;
    const uint64_t message_bytes = (uint64_t)K * block_bytes;
    std::vector<uint8_t> message((size_t)message_bytes, 0x39u);
    wirehair_v2::SeedProfile alternate =
        wirehair_v2::SelectSeedProfile(K, block_bytes);
    alternate.DenseCount = (uint16_t)(alternate.DenseCount + 4u);
    wirehair_v2::Codec alternate_encoder;
    if (alternate_encoder.InitializePrecodeEncoder(
            message.data(), message_bytes, block_bytes, &alternate) !=
            Wirehair_Success ||
        alternate_encoder.Profile().DenseCount != alternate.DenseCount ||
        alternate_encoder.Profile().V2StaircaseCount != alternate.DenseCount)
    {
        std::fprintf(stderr,
            "contract: alternate valid DenseCount did not drive S\n");
        return false;
    }
    wirehair_v2::MessagePrecodeDecoder alternate_decoder;
    if (alternate_decoder.InitializeResult(
            message_bytes, block_bytes, &alternate_encoder.Profile()) !=
            Wirehair_Success ||
        alternate_decoder.System().Params.Staircase != alternate.DenseCount ||
        !CompleteDirectRepairOnly(
            alternate_encoder,
            alternate_decoder,
            message,
            K,
            block_bytes,
            "alternate-dense-count"))
    {
        std::fprintf(stderr,
            "contract: alternate DenseCount profile reuse failed\n");
        return false;
    }

    wirehair_v2::SeedProfile malformed = alternate_encoder.Profile();
    wirehair_v2::MessagePrecodeDecoder rejected;
    const auto reject_profile = [&](const wirehair_v2::SeedProfile& candidate,
                                    const char* field) {
        if (rejected.InitializeResult(
                message_bytes, block_bytes, &candidate) ==
            Wirehair_InvalidInput)
        {
            return true;
        }
        std::fprintf(stderr,
            "contract: forged %s profile was accepted\n", field);
        return false;
    };
    malformed.DenseCount = 0u;
    if (!reject_profile(malformed, "zero DenseCount")) {
        return false;
    }
    malformed = alternate_encoder.Profile();
    ++malformed.V2StaircaseCount;
    if (!reject_profile(malformed, "inconsistent bound S")) {
        return false;
    }
    malformed = alternate_encoder.Profile();
    ++malformed.V2PrecodeContractVersion;
    if (!reject_profile(malformed, "precode version")) {
        return false;
    }
    malformed = alternate_encoder.Profile();
    ++malformed.V2PacketRowContractVersion;
    if (!reject_profile(malformed, "packet version")) {
        return false;
    }
    malformed = alternate_encoder.Profile();
    ++malformed.V2DenseRowCount;
    if (!reject_profile(malformed, "dense-row count")) {
        return false;
    }
    malformed = alternate_encoder.Profile();
    ++malformed.V2HeavyRowCount;
    if (!reject_profile(malformed, "heavy-row count")) {
        return false;
    }
    malformed = alternate_encoder.Profile();
    ++malformed.V2SourceHits;
    if (!reject_profile(malformed, "source-hit count")) {
        return false;
    }
    malformed = alternate_encoder.Profile();
    malformed.V2PrecodeSeed ^= UINT64_C(1);
    if (!reject_profile(malformed, "precode seed")) {
        return false;
    }
    malformed = alternate_encoder.Profile();
    malformed.V2PacketPeelSeed ^= UINT32_C(1);
    if (!reject_profile(malformed, "packet seed")) {
        return false;
    }
    malformed = alternate_encoder.Profile();
    malformed.DenseCount = (uint16_t)(malformed.DenseCount + 4u);
    malformed.V2StaircaseCount += 4u;
    if (!reject_profile(malformed, "paired staircase count")) {
        return false;
    }
    malformed = alternate_encoder.Profile();
    ++malformed.V2SeedAttempt;
    if (!reject_profile(malformed, "seed attempt")) {
        return false;
    }
    malformed = alternate_encoder.Profile();
    malformed.V2PrecodeSeedSalt ^= UINT64_C(1);
    if (!reject_profile(malformed, "precode salt")) {
        return false;
    }
    malformed = alternate_encoder.Profile();
    malformed.V2RecoveryRowSeedSalt ^= UINT64_C(1);
    if (!reject_profile(malformed, "packet salt")) {
        return false;
    }
    malformed = alternate_encoder.Profile();
    malformed.V2DenseIdentityCorner =
        !malformed.V2DenseIdentityCorner;
    if (!reject_profile(malformed, "dense-corner option")) {
        return false;
    }
    malformed = alternate_encoder.Profile();
    malformed.V2DenseTwoAnchor = !malformed.V2DenseTwoAnchor;
    if (!reject_profile(malformed, "active dense-two-anchor state")) {
        return false;
    }
    malformed = alternate_encoder.Profile();
    malformed.V2AdaptiveDenseTwoAnchor =
        !malformed.V2AdaptiveDenseTwoAnchor;
    if (!reject_profile(malformed, "adaptive dense-two-anchor policy")) {
        return false;
    }
    malformed = alternate_encoder.Profile();
    malformed.V2RecoveryMixCount = 2u;
    if (!reject_profile(malformed, "recovery-mix option")) {
        return false;
    }
    malformed = alternate_encoder.Profile();
    malformed.V2SeedSelected = false;
    if (!reject_profile(malformed, "selected-flag downgrade")) {
        return false;
    }

    malformed = wirehair_v2::SelectSeedProfile(K, block_bytes);
    malformed.V2DenseRowCount = 12u;
    if (!reject_profile(malformed, "mixed unselected state")) {
        return false;
    }

    wirehair_v2::MessagePrecodeDecoder preserved;
    if (preserved.InitializeResult(
            message_bytes,
            block_bytes,
            &alternate_encoder.Profile()) != Wirehair_Success)
    {
        return false;
    }
    const wirehair_v2::SeedProfile preserved_profile = preserved.Profile();
    malformed = alternate_encoder.Profile();
    malformed.V2SeedSelected = false;
    if (preserved.InitializeResult(
            message_bytes, block_bytes, &malformed) !=
            Wirehair_InvalidInput ||
        preserved.Profile().V2PrecodeSeed !=
            preserved_profile.V2PrecodeSeed ||
        preserved.System().Params.Staircase !=
            preserved_profile.V2StaircaseCount ||
        alternate_encoder.InitializePrecodeDecoder(
            message_bytes, block_bytes, &malformed) != Wirehair_InvalidInput ||
        !CompleteDirectRepairOnly(
            alternate_encoder,
            preserved,
            message,
            K,
            block_bytes,
            "failed-init-preservation"))
    {
        std::fprintf(stderr,
            "contract: failed init did not preserve decoder/encoder state\n");
        return false;
    }

    wirehair_v2::SeedProfile invalid_unselected =
        wirehair_v2::SelectSeedProfile(K, block_bytes);
    invalid_unselected.DenseCount = 1u;
    wirehair_v2::MessagePrecodeEncoder direct_encoder;
    if (direct_encoder.InitializeResult(
            message.data(), message_bytes, block_bytes, &invalid_unselected) !=
            Wirehair_InvalidInput ||
        rejected.InitializeResult(
            message_bytes, block_bytes, &invalid_unselected) !=
            Wirehair_InvalidInput)
    {
        std::fprintf(stderr, "contract: wrong-congruence S was accepted\n");
        return false;
    }
    invalid_unselected.DenseCount =
        (uint16_t)(wirehair::kMaxDenseCount + 1u);
    if (direct_encoder.InitializeResult(
            message.data(), message_bytes, block_bytes, &invalid_unselected) !=
            Wirehair_InvalidInput ||
        rejected.InitializeResult(
            message_bytes, block_bytes, &invalid_unselected) !=
            Wirehair_InvalidInput)
    {
        std::fprintf(stderr, "contract: oversized S was accepted\n");
        return false;
    }
    return true;
}

bool RunMixedProfileBenchmark(bool mixed)
{
    struct BenchmarkCase {
        uint32_t K;
        uint32_t BlockBytes;
        uint32_t Trials;
    };
    const BenchmarkCase cases[] = {
        {1000u, 1280u, 5u},
        {320u, 102400u, 3u},
        {64u, 1024u * 1024u, 2u}
    };
    for (const BenchmarkCase& c : cases)
    {
        const uint64_t message_bytes = (uint64_t)c.K * c.BlockBytes;
        std::vector<uint8_t> message((size_t)message_bytes);
        for (size_t i = 0; i < message.size(); ++i) {
            message[i] = (uint8_t)(i * 29u + c.K + c.BlockBytes);
        }
        uint64_t create_ns = 0u;
        uint64_t encode_ns = 0u;
        uint64_t decode_ns = 0u;
        uint64_t packets = 0u;
        uint64_t cold_capacity = 0u;
        uint64_t intermediate_bytes = 0u;
        uint32_t selected_attempt = UINT32_MAX;
        for (uint32_t trial = 0; trial < c.Trials; ++trial)
        {
            wirehair_v2::MessagePrecodeEncoderOptions options;
            if (mixed) {
                options.Completion =
                    wirehair_v2::CompletionField::MixedGF256GF16;
            }
            wirehair_v2::MessagePrecodeEncoder encoder;
            std::chrono::steady_clock::time_point begin =
                std::chrono::steady_clock::now();
            if (encoder.InitializeResult(
                    message.data(), message_bytes, c.BlockBytes,
                    nullptr, &options) != Wirehair_Success)
            {
                return false;
            }
            create_ns += (uint64_t)
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now() - begin).count();
            if (selected_attempt == UINT32_MAX) {
                selected_attempt = encoder.Profile().V2SeedAttempt;
            }
            else if (selected_attempt != encoder.Profile().V2SeedAttempt) {
                return false;
            }
            intermediate_bytes = (uint64_t)(
                encoder.Profile().BlockCount +
                encoder.Profile().V2StaircaseCount +
                encoder.Profile().V2DenseRowCount +
                encoder.Profile().V2HeavyRowCount) * c.BlockBytes;

            wirehair_v2::MessagePrecodeDecoder decoder;
            if (decoder.InitializeResult(
                    message_bytes, c.BlockBytes, &encoder.Profile()) !=
                    Wirehair_Success)
            {
                return false;
            }
            std::vector<uint8_t> block(c.BlockBytes);
            WirehairResult result = Wirehair_NeedMore;
            for (uint32_t id = c.K;
                 result == Wirehair_NeedMore && id < 2u * c.K + 128u; ++id)
            {
                uint32_t bytes = 0u;
                begin = std::chrono::steady_clock::now();
                if (encoder.EncodeResult(
                        id, block.data(), c.BlockBytes, &bytes) !=
                        Wirehair_Success)
                {
                    return false;
                }
                encode_ns += (uint64_t)
                    std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::steady_clock::now() - begin).count();
                begin = std::chrono::steady_clock::now();
                result = decoder.DecodeResult(id, block.data(), bytes);
                decode_ns += (uint64_t)
                    std::chrono::duration_cast<std::chrono::nanoseconds>(
                        std::chrono::steady_clock::now() - begin).count();
                cold_capacity = std::max<uint64_t>(
                    cold_capacity,
                    decoder.ColdReceiveCapacityBytesForTesting());
                ++packets;
            }
            std::vector<uint8_t> recovered(message.size());
            if (result != Wirehair_Success ||
                decoder.RecoverResult(
                    recovered.data(), recovered.size()) != Wirehair_Success ||
                recovered != message)
            {
                return false;
            }
        }
        std::printf(
            "mixed_profile_bench,profile=%s,K=%u,bb=%u,trials=%u,"
            "create_ms=%.3f,encode_ns_per_packet=%.1f,decode_ms=%.3f,"
            "packets_per_trial=%.2f,intermediate_mib=%.3f,"
            "max_cold_capacity_mib=%.3f,attempt=%u\n",
            mixed ? "mixed" : "certified", c.K, c.BlockBytes, c.Trials,
            create_ns / 1000000.0 / c.Trials,
            packets == 0u ? 0.0 : (double)encode_ns / packets,
            decode_ns / 1000000.0 / c.Trials,
            (double)packets / c.Trials,
            intermediate_bytes / 1048576.0,
            cold_capacity / 1048576.0,
            selected_attempt);
    }
    return true;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc == 2 &&
        (std::strcmp(argv[1], "--benchmark-certified") == 0 ||
         std::strcmp(argv[1], "--benchmark-mixed") == 0))
    {
        return RunMixedProfileBenchmark(
            std::strcmp(argv[1], "--benchmark-mixed") == 0) ? 0 : 1;
    }
    if (argc != 1) return 2;
    bool ok = true;
    ok = RunFacadeLossCase(64u, 37u, 13u, 7u, false) && ok;
    ok = RunFacadeLossCase(96u, 128u, 91u, 9u, true) && ok;
    ok = RunFacadeLossCase(320u, 17u, 5u, 11u, false) && ok;
    ok = RunFacadeLossCase(64u, 16u, 7u, 7u, false, true) && ok;
    ok = RunFacadeLossCase(64u, 16u, 16u, 7u, true, true) && ok;
    ok = RunDirectLifecycleCase() && ok;
    ok = RunForcedNeedMoreResumeCase() && ok;
    ok = RunMixedColdRetryCase() && ok;
    ok = RunUnauthenticatedCorruptionBoundary() && ok;
    ok = RunExplicitConfigurationCases() && ok;
    ok = RunArchitectureBoundaryCases() && ok;
    ok = RunBoundContractCases() && ok;
    return ok ? 0 : 1;
}
