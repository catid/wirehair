#include "WirehairV2Codec.h"

#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace {

bool Check(bool condition, const char* message)
{
    if (!condition) {
        std::fprintf(stderr, "%s\n", message);
        return false;
    }
    return true;
}

std::vector<uint8_t> MakeMessage(size_t bytes)
{
    std::vector<uint8_t> message(bytes);
    for (size_t i = 0; i < bytes; ++i) {
        message[i] = (uint8_t)(i * 37u + i / 11u + 19u);
    }
    return message;
}

bool RunDecodeCase(
    uint32_t block_count,
    uint32_t block_bytes,
    bool recovery_only,
    bool partial_final,
    bool keep_all_sources,
    uint32_t loss_modulus)
{
    const uint64_t full_bytes = (uint64_t)block_count * block_bytes;
    const uint64_t message_bytes = partial_final ?
        full_bytes - (block_bytes > 1u ? block_bytes / 2u : 0u) :
        full_bytes;
    std::vector<uint8_t> message = MakeMessage((size_t)message_bytes);
    std::vector<uint8_t> recovered((size_t)message_bytes, 0u);
    std::vector<uint8_t> block(block_bytes, 0u);

    wirehair_v2::Codec encoder;
    wirehair_v2::Codec decoder;
    if (!Check(
            encoder.InitializePrecodeEncoder(
                message.data(), message_bytes, block_bytes) ==
                Wirehair_Success,
            "precode encoder initialization failed") ||
        !Check(
            decoder.InitializePrecodeDecoder(
                message_bytes, block_bytes) ==
                Wirehair_Success,
            "precode decoder initialization failed"))
    {
        return false;
    }

    WirehairResult decode_result = Wirehair_NeedMore;
    if (recovery_only)
    {
        const uint32_t shuffled_ids[] = {
            UINT32_MAX, block_count + 17u, block_count + 3u
        };
        for (const uint32_t id : shuffled_ids)
        {
            uint32_t written = 0u;
            if (!Check(
                    encoder.Encode(
                        id, block.data(), block_bytes, &written) ==
                        Wirehair_Success && written == block_bytes,
                    "shuffled recovery encode failed"))
            {
                return false;
            }
            decode_result = decoder.Decode(id, block.data(), written);
            if (!Check(
                    decode_result == Wirehair_NeedMore,
                    "shuffled recovery decode ended early"))
            {
                return false;
            }
        }
        uint32_t written = 0u;
        if (!Check(
                encoder.Encode(
                    shuffled_ids[0], block.data(), block_bytes, &written) ==
                    Wirehair_Success,
                "duplicate high-id recovery encode failed") ||
            !Check(
                decoder.Decode(shuffled_ids[0], block.data(), written) ==
                    Wirehair_NeedMore,
                "duplicate high-id recovery packet was not ignored"))
        {
            return false;
        }
    }
    if (!recovery_only)
    {
        for (uint32_t id = 0; id < block_count; ++id)
        {
            if (!keep_all_sources && id % loss_modulus == 0u) {
                continue;
            }
            uint32_t written = 0u;
            if (!Check(
                    encoder.Encode(
                        id, block.data(), block_bytes, &written) ==
                        Wirehair_Success,
                    "systematic encode failed"))
            {
                return false;
            }
            decode_result = decoder.Decode(id, block.data(), written);
            const WirehairResult expected =
                keep_all_sources && id + 1u == block_count ?
                    Wirehair_Success : Wirehair_NeedMore;
            if (!Check(
                    decode_result == expected,
                    "systematic decode completed at the wrong point"))
            {
                return false;
            }
        }

        uint32_t written = 0u;
        if (!keep_all_sources &&
            (!Check(
                 encoder.Encode(
                     1u, block.data(), block_bytes, &written) ==
                     Wirehair_Success,
                 "duplicate source encode failed") ||
             !Check(
                 decoder.Decode(1u, block.data(), written) ==
                     Wirehair_NeedMore,
                 "duplicate packet was not ignored")))
        {
            return false;
        }
    }

    for (uint32_t recovery = 0;
        recovery < block_count + 512u &&
            decode_result == Wirehair_NeedMore;
        ++recovery)
    {
        const uint32_t id = block_count + recovery;
        uint32_t written = 0u;
        if (!Check(
                encoder.Encode(
                    id, block.data(), block_bytes, &written) ==
                    Wirehair_Success &&
                    written == block_bytes,
                "recovery encode failed"))
        {
            return false;
        }
        decode_result = decoder.Decode(id, block.data(), written);
    }

    if (!Check(
            decode_result == Wirehair_Success,
            "precode decode did not reach full rank") ||
        !Check(
            decoder.Recover(recovered.data(), message_bytes) ==
                Wirehair_Success,
            "precode recovery failed") ||
        !Check(
            recovered == message,
            "precode recovered message mismatch") ||
        !Check(
            decoder.Encode(0u, block.data(), block_bytes, &block_count) ==
                Wirehair_InvalidInput,
            "decoder mode unexpectedly allowed Encode"))
    {
        return false;
    }
    return true;
}

bool CheckInvalidAndOom()
{
    wirehair_v2::MessagePrecodeEncoderOptions options;
    wirehair_v2::MessagePrecodeDecoder decoder;
    if (!Check(
            decoder.InitializeResult(1024u, 16u, nullptr, &options) ==
                Wirehair_Success,
            "direct decoder initialization failed") ||
        !Check(
            decoder.InitializeResult(UINT64_MAX, 16u, nullptr, &options) ==
                Wirehair_InvalidInput && decoder.IsInitialized(),
            "decoder overflow input was not rejected transactionally") ||
        !Check(
            decoder.InitializeResult(
                UINT64_C(0x100000000), UINT32_C(0x80000000),
                nullptr, &options) == Wirehair_InvalidInput &&
                decoder.IsInitialized(),
            "decoder oversized block input was not rejected transactionally") ||
        !Check(
            decoder.RecoverResult(nullptr, 1024u) == Wirehair_InvalidInput,
            "null recovery output was accepted"))
    {
        return false;
    }
    uint8_t block[16] = {};
    if (!Check(
            decoder.DecodeResult(0u, nullptr, 16u) == Wirehair_InvalidInput,
            "null decode input was accepted") ||
        !Check(
            decoder.DecodeResult(0u, block, 15u) == Wirehair_InvalidInput,
            "short non-final source packet was accepted"))
    {
        return false;
    }

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    const uint32_t K = 64u;
    const uint32_t solve_block_bytes = 37u;
    std::vector<uint8_t> source((size_t)K * solve_block_bytes, 0u);
    wirehair_v2::MessagePrecodeDecoder solve_oom_decoder;
    if (!Check(
            solve_oom_decoder.InitializeResult(
                (uint64_t)K * solve_block_bytes,
                solve_block_bytes,
                nullptr,
                &options) == Wirehair_Success,
            "solve OOM decoder initialization failed"))
    {
        return false;
    }
    for (uint32_t id = 0; id + 1u < K; ++id) {
        uint8_t* const packet =
            source.data() + (size_t)id * solve_block_bytes;
        for (uint32_t b = 0; b < solve_block_bytes; ++b) {
            packet[b] = (uint8_t)(id * 7u + b * 11u + 3u);
        }
        if (!Check(
                solve_oom_decoder.DecodeResult(
                    id, packet, solve_block_bytes) == Wirehair_NeedMore,
                "solve OOM decoder ended early"))
        {
            return false;
        }
    }
    uint8_t* const final_packet =
        source.data() + (size_t)(K - 1u) * solve_block_bytes;
    for (uint32_t b = 0; b < solve_block_bytes; ++b) {
        final_packet[b] = (uint8_t)((K - 1u) * 7u + b * 11u + 3u);
    }
    wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(0);
    const WirehairResult solve_oom = solve_oom_decoder.DecodeResult(
        K - 1u, final_packet, solve_block_bytes);
    wirehair_v2::SetDecoderAllocationFailureCountdownForTesting(-1);
    const bool retryable_state =
        solve_oom == Wirehair_OOM &&
        solve_oom_decoder.ReceivedCount() == K &&
        !solve_oom_decoder.IsDecoded();
    const WirehairResult retry_result = solve_oom_decoder.DecodeResult(
        K - 1u, final_packet, solve_block_bytes);
    std::vector<uint8_t> recovered(source.size(), 0u);
    if (!Check(
            retryable_state,
            "solve-path OOM did not preserve retryable packet state") ||
        !Check(
            retry_result == Wirehair_Success,
            "duplicate retry after solve-path OOM did not decode") ||
        !Check(
            solve_oom_decoder.RecoverResult(
                recovered.data(), recovered.size()) == Wirehair_Success &&
                recovered == source,
            "duplicate retry after solve-path OOM recovered wrong data"))
    {
        return false;
    }
#endif
    return true;
}

bool CheckFacadeModeTransitions()
{
    const uint32_t K = 16u;
    const uint32_t block_bytes = 11u;
    const uint64_t message_bytes = (uint64_t)K * block_bytes;
    std::vector<uint8_t> message = MakeMessage((size_t)message_bytes);
    std::vector<uint8_t> block(block_bytes, 0u);
    uint32_t written = 0u;

    wirehair_v2::Codec codec;
    if (!Check(
            codec.InitializeEncoder(
                message.data(), message_bytes, block_bytes) ==
                Wirehair_Success,
            "initial V1 encoder mode failed"))
    {
        return false;
    }
    wirehair_v2::SeedProfile mismatch = codec.Profile();
    ++mismatch.BlockCount;
    if (!Check(
            codec.InitializePrecodeDecoder(
                message_bytes, block_bytes, &mismatch) ==
                Wirehair_InvalidInput,
            "mismatched V2 decoder profile was accepted") ||
        !Check(
            codec.Encode(0u, block.data(), block_bytes, &written) ==
                Wirehair_Success,
            "failed V2 decoder init did not preserve V1 encoder mode"))
    {
        return false;
    }

    if (!Check(
            codec.InitializePrecodeDecoder(message_bytes, block_bytes) ==
                Wirehair_Success,
            "V1 encoder to V2 decoder transition failed") ||
        !Check(
            codec.Encode(0u, block.data(), block_bytes, &written) ==
                Wirehair_InvalidInput,
            "V2 decoder mode allowed Encode") ||
        !Check(
            codec.Recover(message.data(), message_bytes) == Wirehair_NeedMore,
            "fresh V2 decoder did not request packets"))
    {
        return false;
    }

    const wirehair_v2::SeedProfile decoder_profile = codec.Profile();
    if (!Check(
            codec.InitializePrecodeDecoder(0u, block_bytes) ==
                Wirehair_InvalidInput,
            "invalid V2 decoder reinitialize was accepted") ||
        !Check(
            codec.Profile().BlockCount == decoder_profile.BlockCount &&
                codec.Encode(0u, block.data(), block_bytes, &written) ==
                    Wirehair_InvalidInput,
            "failed V2 decoder reinitialize changed active mode"))
    {
        return false;
    }

    if (!Check(
            codec.InitializePrecodeEncoder(
                message.data(), message_bytes, block_bytes) ==
                Wirehair_Success,
            "V2 decoder to V2 encoder transition failed") ||
        !Check(
            codec.Encode(K, block.data(), block_bytes, &written) ==
                Wirehair_Success,
            "V2 encoder unusable after decoder transition") ||
        !Check(
            codec.Decode(K, block.data(), written) == Wirehair_InvalidInput,
            "V2 encoder mode allowed Decode"))
    {
        return false;
    }

    if (!Check(
            codec.InitializeDecoder(message_bytes, block_bytes) ==
                Wirehair_Success,
            "V2 encoder to V1 decoder transition failed") ||
        !Check(
            codec.InitializeEncoder(
                message.data(), message_bytes, block_bytes) ==
                Wirehair_Success,
            "V1 decoder to V1 encoder transition failed") ||
        !Check(
            codec.Encode(0u, block.data(), block_bytes, &written) ==
                Wirehair_Success,
            "V1 encoder unusable after mode transitions"))
    {
        return false;
    }
    return true;
}

} // namespace

int main(int argc, char** argv)
{
    const bool large = argc == 2 && std::string(argv[1]) == "--large";
    const bool large_recovery =
        argc == 2 && std::string(argv[1]) == "--large-recovery";
    if (argc > 2 || (argc == 2 && !large && !large_recovery)) {
        std::fprintf(
            stderr, "usage: %s [--large|--large-recovery]\n", argv[0]);
        return 2;
    }
    if (!CheckInvalidAndOom() || !CheckFacadeModeTransitions()) {
        return 1;
    }
    if (large) {
        return RunDecodeCase(64000u, 1u, false, false, true, 10u) ? 0 : 1;
    }
    if (large_recovery) {
        return RunDecodeCase(64000u, 1u, true, false, false, 10u) ? 0 : 1;
    }
    if (!RunDecodeCase(64u, 37u, false, true, false, 5u) ||
        !RunDecodeCase(1000u, 3u, false, false, false, 5u) ||
        !RunDecodeCase(32u, 19u, true, true, false, 5u))
    {
        return 1;
    }
    return 0;
}
