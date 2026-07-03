#include "WirehairV2Codec.h"

#include <cstdio>
#include <cstring>
#include <vector>

int main()
{
    const uint32_t N = 4u;
    const uint32_t block_bytes = 16u;
    const uint64_t message_bytes = (uint64_t)N * block_bytes;

    std::vector<uint8_t> message((size_t)message_bytes);
    std::vector<uint8_t> recovered((size_t)message_bytes, 0);
    std::vector<uint8_t> block(block_bytes);
    for (size_t i = 0; i < message.size(); ++i) {
        message[i] = (uint8_t)(i * 13u + 7u);
    }

    wirehair_v2::Codec encoder;
    wirehair_v2::Codec decoder;
    WirehairResult result = encoder.InitializeEncoder(
        &message[0], message_bytes, block_bytes);
    if (result != Wirehair_Success) {
        std::fprintf(stderr, "encoder init failed: %d\n", result);
        return 1;
    }
    result = decoder.InitializeDecoder(message_bytes, block_bytes);
    if (result != Wirehair_Success) {
        std::fprintf(stderr, "decoder init failed: %d\n", result);
        return 1;
    }

    for (uint32_t block_id = 0; block_id < N; ++block_id)
    {
        uint32_t written = 0;
        result = encoder.Encode(block_id, &block[0], block_bytes, &written);
        if (result != Wirehair_Success) {
            std::fprintf(stderr, "encode failed: %d\n", result);
            return 1;
        }
        result = decoder.Decode(block_id, &block[0], written);
        if (block_id + 1u < N && result != Wirehair_NeedMore) {
            std::fprintf(stderr, "decode ended early: %d\n", result);
            return 1;
        }
    }
    if (result != Wirehair_Success) {
        std::fprintf(stderr, "decode did not complete: %d\n", result);
        return 1;
    }
    result = decoder.Recover(&recovered[0], message_bytes);
    if (result != Wirehair_Success) {
        std::fprintf(stderr, "recover failed: %d\n", result);
        return 1;
    }
    if (std::memcmp(&recovered[0], &message[0], (size_t)message_bytes) != 0) {
        std::fprintf(stderr, "recovered message mismatch\n");
        return 1;
    }
    return 0;
}
