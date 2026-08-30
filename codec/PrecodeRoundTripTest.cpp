#include "WirehairV2Codec.h"
#include "WirehairV2Plan.h"
#include "WirehairV2PrecodeDecode.h"

#include "../WirehairTools.h"

#include <algorithm>
#include <cstdio>
#include <chrono>
#include <cstring>
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
    bool repair_only)
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
        repair_only ? "repair-only" : "reverse-loss",
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
        recovered == message || decoder.IntermediateBlocks() != nullptr)
    {
        std::fprintf(stderr,
            "integrity boundary: K-row corruption contract changed\n");
        return false;
    }

    uint32_t bytes = 0u;
    if (encoder.Encode(K, block.data(), block_bytes, &bytes) !=
            Wirehair_Success ||
        decoder.DecodeResult(K, block.data(), bytes) != Wirehair_Error ||
        decoder.IntermediateBlocks() == nullptr)
    {
        std::fprintf(stderr,
            "integrity boundary: overdetermined conflict was not rejected\n");
        return false;
    }
    return true;
}

bool CompleteDirectRepairOnly(
    wirehair_v2::Codec& encoder,
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
        if (encoder.Encode(id, block.data(), block_bytes, &bytes) !=
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
    return a.RecoveryMixCount == b.RecoveryMixCount &&
        a.DenseAnchors == b.DenseAnchors &&
        a.DenseIdentityCorner == b.DenseIdentityCorner &&
        a.PrecodeSeedSalt == b.PrecodeSeedSalt &&
        a.RecoveryRowSeedSalt == b.RecoveryRowSeedSalt;
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
        a.DenseIdentityCorner == b.DenseIdentityCorner &&
        a.HeavyFamily == b.HeavyFamily &&
        a.DenseAnchors == b.DenseAnchors &&
        a.Seed == b.Seed;
}

bool SamePacketConfig(
    const wirehair_v2::PacketRowConfig& a,
    const wirehair_v2::PacketRowConfig& b)
{
    return a.PeelSeed == b.PeelSeed && a.MixCount == b.MixCount;
}

bool RunGraphSeedWidthContractCases()
{
    const uint32_t K = 64u;
    const uint32_t narrow_block_bytes =
        wirehair_v2::kTwo07Mix2GraphSeedBlockBytes;
    // Exercise a different payload-width policy class, not just a different
    // numeric width inside the small-block class.
    const uint32_t wide_block_bytes = 320u * 1024u;
    const wirehair_v2::SeedProfile narrow =
        wirehair_v2::SelectSeedProfile(K, narrow_block_bytes);
    const wirehair_v2::SeedProfile wide =
        wirehair_v2::SelectSeedProfile(K, wide_block_bytes);

    wirehair_v2::MessagePrecodeEncoderOptions certified_options;
    wirehair_v2::PrecodeParams certified_narrow;
    wirehair_v2::PrecodeParams certified_wide;
    wirehair_v2::PacketRowConfig certified_narrow_packet;
    wirehair_v2::PacketRowConfig certified_wide_packet;
    wirehair_v2::MakeMessagePrecodeConfiguration(
        narrow, certified_options, certified_narrow,
        certified_narrow_packet);
    wirehair_v2::MakeMessagePrecodeConfiguration(
        wide, certified_options, certified_wide,
        certified_wide_packet);

    wirehair_v2::PrecodeParams expected_certified_narrow =
        wirehair_v2::MakeCertifiedParams(
            K,
            wirehair_v2::MatrixSeedFromProfile(
                narrow, 0u, wirehair_v2::kMessagePrecodeSeedSalt));
    wirehair_v2::PrecodeParams expected_certified_wide =
        wirehair_v2::MakeCertifiedParams(
            K,
            wirehair_v2::MatrixSeedFromProfile(
                wide, 0u, wirehair_v2::kMessagePrecodeSeedSalt));
    wirehair_v2::PacketRowConfig expected_certified_narrow_packet;
    expected_certified_narrow_packet.PeelSeed =
        wirehair_v2::PacketPeelSeedFromProfile(
            narrow, wirehair_v2::kMessageRecoveryRowSeedSalt);
    wirehair_v2::PacketRowConfig expected_certified_wide_packet;
    expected_certified_wide_packet.PeelSeed =
        wirehair_v2::PacketPeelSeedFromProfile(
            wide, wirehair_v2::kMessageRecoveryRowSeedSalt);

    if (narrow.BlockBytes != narrow_block_bytes ||
        wide.BlockBytes != wide_block_bytes ||
        narrow.Policy.ByteClass == wide.Policy.ByteClass ||
        !SamePrecodeParams(
            certified_narrow, expected_certified_narrow) ||
        !SamePrecodeParams(certified_wide, expected_certified_wide) ||
        !SamePacketConfig(
            certified_narrow_packet,
            expected_certified_narrow_packet) ||
        !SamePacketConfig(
            certified_wide_packet,
            expected_certified_wide_packet) ||
        certified_narrow.Seed != UINT64_C(0x3a880cf93d36784a) ||
        certified_narrow_packet.PeelSeed != UINT32_C(0x4df9ea57) ||
        certified_wide.Seed != UINT64_C(0x0e6034e9935ef8e3) ||
        certified_wide_packet.PeelSeed != UINT32_C(0x7902a554) ||
        certified_narrow.Seed == certified_wide.Seed ||
        certified_narrow_packet.PeelSeed == certified_wide_packet.PeelSeed)
    {
        std::fprintf(stderr,
            "contract: certified profile cross-width seed behavior changed\n");
        return false;
    }

    wirehair_v2::MessagePrecodeEncoderOptions two07_mix2;
    two07_mix2.DenseAnchors = wirehair_v2::DenseAnchorLayout::Two07;
    two07_mix2.RecoveryMixCount = 2u;
    wirehair_v2::PrecodeParams two07_narrow;
    wirehair_v2::PrecodeParams two07_wide;
    wirehair_v2::PacketRowConfig two07_narrow_packet;
    wirehair_v2::PacketRowConfig two07_wide_packet;
    wirehair_v2::MakeMessagePrecodeConfiguration(
        narrow, two07_mix2, two07_narrow, two07_narrow_packet);
    wirehair_v2::MakeMessagePrecodeConfiguration(
        wide, two07_mix2, two07_wide, two07_wide_packet);

    const wirehair_v2::SeedProfile normalized =
        wirehair_v2::SelectSeedProfile(
            K, wirehair_v2::kTwo07Mix2GraphSeedBlockBytes);
    wirehair_v2::PrecodeParams expected_two07 =
        wirehair_v2::MakeCertifiedParams(
            K,
            wirehair_v2::MatrixSeedFromProfile(
                normalized, 0u, wirehair_v2::kMessagePrecodeSeedSalt));
    expected_two07.DenseAnchors = wirehair_v2::DenseAnchorLayout::Two07;
    wirehair_v2::PacketRowConfig expected_two07_packet;
    expected_two07_packet.PeelSeed =
        wirehair_v2::PacketPeelSeedFromProfile(
            normalized, wirehair_v2::kMessageRecoveryRowSeedSalt);
    expected_two07_packet.MixCount = 2u;

    if (!SamePrecodeParams(two07_narrow, two07_wide) ||
        !SamePacketConfig(two07_narrow_packet, two07_wide_packet) ||
        !SamePrecodeParams(two07_narrow, expected_two07) ||
        !SamePacketConfig(two07_narrow_packet, expected_two07_packet) ||
        two07_narrow.Seed != UINT64_C(0x3a880cf93d36784a) ||
        two07_narrow_packet.PeelSeed != UINT32_C(0x4df9ea57) ||
        two07_narrow.DenseAnchors !=
            wirehair_v2::DenseAnchorLayout::Two07 ||
        two07_narrow_packet.MixCount != 2u)
    {
        std::fprintf(stderr,
            "contract: Two07/mix2 graph seeds depend on payload width\n");
        return false;
    }

    const uint32_t attempts[] = { 0u, 1u, 2u, 255u };
    for (uint32_t attempt : attempts)
    {
        const wirehair_v2::PrecodeParams narrow_attempt =
            wirehair_v2::PrecodeParamsForAttempt(two07_narrow, attempt);
        const wirehair_v2::PrecodeParams wide_attempt =
            wirehair_v2::PrecodeParamsForAttempt(two07_wide, attempt);
        const wirehair_v2::PacketRowConfig narrow_packet_attempt =
            wirehair_v2::PacketConfigForAttempt(
                two07_narrow_packet, attempt);
        const wirehair_v2::PacketRowConfig wide_packet_attempt =
            wirehair_v2::PacketConfigForAttempt(two07_wide_packet, attempt);
        if (!SamePrecodeParams(narrow_attempt, wide_attempt) ||
            !SamePacketConfig(narrow_packet_attempt, wide_packet_attempt))
        {
            std::fprintf(stderr,
                "contract: Two07/mix2 attempt %u depends on payload width\n",
                attempt);
            return false;
        }
    }
    const wirehair_v2::PrecodeParams attempt_one =
        wirehair_v2::PrecodeParamsForAttempt(two07_narrow, 1u);
    const wirehair_v2::PacketRowConfig packet_attempt_one =
        wirehair_v2::PacketConfigForAttempt(two07_narrow_packet, 1u);
    if (attempt_one.Seed != UINT64_C(0xd8bf86b2bc80f45f) ||
        packet_attempt_one.PeelSeed != UINT32_C(0xec316410))
    {
        std::fprintf(stderr,
            "contract: Two07/mix2 normalized attempt stepping changed\n");
        return false;
    }
    return true;
}

bool RunOptionContractCase(
    const char* label,
    const wirehair_v2::MessagePrecodeEncoderOptions& options,
    const wirehair_v2::MessagePrecodeEncoderOptions& mismatched)
{
    const uint32_t K = 64u;
    const uint32_t block_bytes = 17u;
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
            wirehair_v2::PrecodeContractVersion() ||
        profile.V2PacketRowContractVersion !=
            wirehair_v2::kPacketRowContractVersion ||
        profile.V2StaircaseCount != profile.DenseCount ||
        profile.V2RecoveryMixCount != options.RecoveryMixCount ||
        profile.V2DenseAnchorLayout != (uint32_t)options.DenseAnchors ||
        profile.V2DenseIdentityCorner != options.DenseIdentityCorner ||
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
        inherited.System().Params.SourceHits != profile.V2SourceHits ||
        inherited.System().Params.DenseAnchors != options.DenseAnchors ||
        inherited.System().Params.Seed != profile.V2PrecodeSeed ||
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
    if (!RunGraphSeedWidthContractCases()) {
        return false;
    }
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
            "contract: unsupported recovery mix count was accepted\n");
        return false;
    }
    variant = defaults;
    variant.DenseAnchors = wirehair_v2::DenseAnchorLayout::Two07;
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
            "contract: Two07/mix3 pair was accepted\n");
        return false;
    }
    variant.RecoveryMixCount = 2u;
    variant.DenseIdentityCorner = true;
    const wirehair_v2::SeedProfile invalid_profile =
        wirehair_v2::SelectSeedProfile(invalid_K, invalid_block_bytes);
    wirehair_v2::MessagePrecodeEncoderOptions ignored_options;
    wirehair_v2::PrecodeParams ignored_params;
    wirehair_v2::PacketRowConfig ignored_packet;
    if (wirehair_v2::ResolveMessagePrecodeOptions(
            invalid_profile, &variant, ignored_options) ||
        wirehair_v2::ResolveMessagePrecodeConfiguration(
            invalid_profile, variant, ignored_params, ignored_packet))
    {
        std::fprintf(stderr,
            "contract: Two07/mix2 identity corner was accepted\n");
        return false;
    }
    variant.DenseIdentityCorner = false;
    wirehair_v2::SeedProfile noncanonical_two07 = invalid_profile;
    noncanonical_two07.DenseCount =
        (uint16_t)(noncanonical_two07.DenseCount + 4u);
    if (!wirehair_v2::ResolveMessagePrecodeOptions(
            noncanonical_two07, &variant, ignored_options) ||
        wirehair_v2::ResolveMessagePrecodeConfiguration(
            noncanonical_two07,
            variant,
            ignored_params,
            ignored_packet) ||
        invalid_encoder.InitializePrecodeEncoder(
            invalid_message.data(),
            invalid_message_bytes,
            invalid_block_bytes,
            &noncanonical_two07,
            &variant) != Wirehair_InvalidInput ||
        invalid_decoder.InitializeResult(
            invalid_message_bytes,
            invalid_block_bytes,
            &noncanonical_two07,
            &variant) != Wirehair_InvalidInput)
    {
        std::fprintf(stderr,
            "contract: noncanonical Two07 staircase was accepted\n");
        return false;
    }
    if (!RunOptionContractCase("two07-mix2", variant, defaults)) {
        return false;
    }
    variant.DenseAnchors = wirehair_v2::DenseAnchorLayout::Four0369;
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
            "contract: Four0369/mix2 pair was accepted\n");
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
    malformed.V2DenseAnchorLayout =
        (uint32_t)wirehair_v2::DenseAnchorLayout::Two07;
    if (!reject_profile(malformed, "dense-anchor option")) {
        return false;
    }
    malformed = alternate_encoder.Profile();
    malformed.V2DenseIdentityCorner =
        !malformed.V2DenseIdentityCorner;
    if (!reject_profile(malformed, "dense-corner option")) {
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
    if (!reject_profile(malformed, "unselected contract state")) {
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

bool RunCertifiedProfileBenchmark()
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
            "precode_profile_bench,profile=certified,K=%u,bb=%u,trials=%u,"
            "create_ms=%.3f,encode_ns_per_packet=%.1f,decode_ms=%.3f,"
            "packets_per_trial=%.2f,intermediate_mib=%.3f,"
            "max_cold_capacity_mib=%.3f,attempt=%u\n",
            c.K, c.BlockBytes, c.Trials,
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
        std::strcmp(argv[1], "--benchmark-certified") == 0)
    {
        return RunCertifiedProfileBenchmark() ? 0 : 1;
    }
    if (argc != 1) return 2;
    bool ok = true;
    ok = RunFacadeLossCase(64u, 37u, 13u, 7u, false) && ok;
    ok = RunFacadeLossCase(96u, 128u, 91u, 9u, true) && ok;
    ok = RunFacadeLossCase(320u, 17u, 5u, 11u, false) && ok;
    ok = RunDirectLifecycleCase() && ok;
    ok = RunForcedNeedMoreResumeCase() && ok;
    ok = RunUnauthenticatedCorruptionBoundary() && ok;
    ok = RunBoundContractCases() && ok;
    return ok ? 0 : 1;
}
