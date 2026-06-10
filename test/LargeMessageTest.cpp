#include <wirehair/wirehair.h>

#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <memory>
#include <new>

namespace {

uint64_t Mix64(uint64_t x)
{
    x ^= x >> 30;
    x *= UINT64_C(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x *= UINT64_C(0x94d049bb133111eb);
    x ^= x >> 31;
    return x;
}

bool ReadEnvU64(const char* name, uint64_t& out)
{
    const char* value = std::getenv(name);
    if (!value || !value[0]) {
        return true;
    }
    if (value[0] == '-') {
        std::cerr << "Invalid " << name << " = " << value << std::endl;
        return false;
    }

    errno = 0;
    char* end = nullptr;
    const unsigned long long parsed = std::strtoull(value, &end, 0);
    if (errno != 0 || !end || *end != '\0') {
        std::cerr << "Invalid " << name << " = " << value << std::endl;
        return false;
    }

    out = static_cast<uint64_t>(parsed);
    return true;
}

void FillBlock(uint32_t block_id, uint8_t* block, uint32_t bytes)
{
    uint64_t state = Mix64(UINT64_C(0x99f0b41d5b6a71e3) ^ block_id);

    uint32_t offset = 0;
    while (offset + sizeof(uint64_t) <= bytes) {
        state = Mix64(state + UINT64_C(0x9e3779b97f4a7c15));
        std::memcpy(block + offset, &state, sizeof(state));
        offset += sizeof(state);
    }

    if (offset < bytes) {
        state = Mix64(state + UINT64_C(0x9e3779b97f4a7c15));
        std::memcpy(block + offset, &state, bytes - offset);
    }
}

bool FillMessage(uint8_t* message, uint32_t block_count, uint32_t block_bytes)
{
    for (uint32_t block_id = 0; block_id < block_count; ++block_id) {
        FillBlock(
            block_id,
            message + static_cast<size_t>(block_id) * block_bytes,
            block_bytes);

        if ((block_id & 1023u) == 1023u) {
            std::cout << "filled blocks: " << (block_id + 1u) << std::endl;
        }
    }

    return true;
}

bool VerifyMessage(const uint8_t* recovered, uint32_t block_count, uint32_t block_bytes)
{
    std::unique_ptr<uint8_t[]> expected(new (std::nothrow) uint8_t[block_bytes]);
    if (!expected) {
        std::cerr << "Failed to allocate verification block" << std::endl;
        return false;
    }

    for (uint32_t block_id = 0; block_id < block_count; ++block_id) {
        FillBlock(block_id, expected.get(), block_bytes);
        const uint8_t* actual =
            recovered + static_cast<size_t>(block_id) * block_bytes;
        if (std::memcmp(actual, expected.get(), block_bytes) != 0) {
            std::cerr << "Recovered block mismatch at block " << block_id
                      << std::endl;
            return false;
        }

        if ((block_id & 1023u) == 1023u) {
            std::cout << "verified blocks: " << (block_id + 1u) << std::endl;
        }
    }

    return true;
}

} // namespace

int main()
{
    uint64_t block_count64 = 8192;
    uint64_t block_bytes64 = 1024 * 1024;
    uint64_t drop_block64 = block_count64 / 2;
    uint64_t max_repair64 = 256;

    if (!ReadEnvU64("WIREHAIR_LARGE_N", block_count64) ||
        !ReadEnvU64("WIREHAIR_LARGE_BLOCK_BYTES", block_bytes64) ||
        !ReadEnvU64("WIREHAIR_LARGE_DROP_BLOCK", drop_block64) ||
        !ReadEnvU64("WIREHAIR_LARGE_MAX_REPAIR", max_repair64))
    {
        return 2;
    }

    if (block_count64 < 2 ||
        block_count64 > std::numeric_limits<uint32_t>::max() ||
        block_bytes64 == 0 ||
        block_bytes64 > 0x7fffffffu ||
        drop_block64 >= block_count64 ||
        max_repair64 == 0 ||
        max_repair64 > std::numeric_limits<uint32_t>::max() ||
        block_count64 > std::numeric_limits<uint64_t>::max() / block_bytes64 ||
        block_count64 - 1 >
            std::numeric_limits<uint32_t>::max() - max_repair64)
    {
        std::cerr << "Invalid large-message test parameters" << std::endl;
        return 2;
    }

    const uint32_t block_count = static_cast<uint32_t>(block_count64);
    const uint32_t block_bytes = static_cast<uint32_t>(block_bytes64);
    const uint32_t drop_block = static_cast<uint32_t>(drop_block64);
    const uint32_t max_repair = static_cast<uint32_t>(max_repair64);
    const uint64_t message_bytes = block_count64 * block_bytes64;

    if (message_bytes / block_bytes != block_count ||
        message_bytes > static_cast<uint64_t>(std::numeric_limits<size_t>::max()))
    {
        std::cerr << "Message size overflows this platform" << std::endl;
        return 2;
    }

    std::cout << "large_message_test: N=" << block_count
              << " block_bytes=" << block_bytes
              << " message_bytes=" << message_bytes
              << " drop_block=" << drop_block << std::endl;

    const WirehairResult init_result = wirehair_init();
    if (init_result != Wirehair_Success) {
        std::cerr << "wirehair_init failed: " << init_result << std::endl;
        return 1;
    }

    std::unique_ptr<uint8_t[]> message(
        new (std::nothrow) uint8_t[static_cast<size_t>(message_bytes)]);
    std::unique_ptr<uint8_t[]> recovered(
        new (std::nothrow) uint8_t[static_cast<size_t>(message_bytes)]);
    std::unique_ptr<uint8_t[]> repair(new (std::nothrow) uint8_t[block_bytes]);
    if (!message || !recovered || !repair) {
        std::cerr << "Failed to allocate large-message buffers" << std::endl;
        return 1;
    }

    FillMessage(message.get(), block_count, block_bytes);

    WirehairCodec encoder = wirehair_encoder_create(
        nullptr, message.get(), message_bytes, block_bytes);
    if (!encoder) {
        std::cerr << "wirehair_encoder_create failed" << std::endl;
        return 1;
    }

    WirehairCodec decoder = wirehair_decoder_create(
        nullptr, message_bytes, block_bytes);
    if (!decoder) {
        std::cerr << "wirehair_decoder_create failed" << std::endl;
        wirehair_free(encoder);
        return 1;
    }

    WirehairResult decode_result = Wirehair_NeedMore;
    for (uint32_t block_id = 0; block_id < block_count; ++block_id) {
        if (block_id == drop_block) {
            continue;
        }

        const uint8_t* block =
            message.get() + static_cast<size_t>(block_id) * block_bytes;
        decode_result = wirehair_decode(decoder, block_id, block, block_bytes);
        if (decode_result != Wirehair_NeedMore &&
            decode_result != Wirehair_Success)
        {
            std::cerr << "wirehair_decode original failed at block "
                      << block_id << ": " << decode_result << std::endl;
            wirehair_free(decoder);
            wirehair_free(encoder);
            return 1;
        }
    }

    uint32_t repair_count = 0;
    for (uint32_t repair_id = block_count;
         decode_result != Wirehair_Success && repair_count < max_repair;
         ++repair_id, ++repair_count)
    {
        uint32_t bytes_out = 0;
        const WirehairResult encode_result = wirehair_encode(
            encoder, repair_id, repair.get(), block_bytes, &bytes_out);
        if (encode_result != Wirehair_Success) {
            std::cerr << "wirehair_encode repair failed at id " << repair_id
                      << ": " << encode_result << std::endl;
            wirehair_free(decoder);
            wirehair_free(encoder);
            return 1;
        }

        decode_result = wirehair_decode(
            decoder, repair_id, repair.get(), bytes_out);
        if (decode_result != Wirehair_NeedMore &&
            decode_result != Wirehair_Success)
        {
            std::cerr << "wirehair_decode repair failed at id " << repair_id
                      << ": " << decode_result << std::endl;
            wirehair_free(decoder);
            wirehair_free(encoder);
            return 1;
        }
    }

    if (decode_result != Wirehair_Success) {
        std::cerr << "Decoder still needs more after " << repair_count
                  << " repair blocks" << std::endl;
        wirehair_free(decoder);
        wirehair_free(encoder);
        return 1;
    }

    std::cout << "decode succeeded after repair blocks: " << repair_count
              << std::endl;

    const WirehairResult recover_result =
        wirehair_recover(decoder, recovered.get(), message_bytes);
    if (recover_result != Wirehair_Success) {
        std::cerr << "wirehair_recover failed: " << recover_result
                  << std::endl;
        wirehair_free(decoder);
        wirehair_free(encoder);
        return 1;
    }

    wirehair_free(decoder);
    wirehair_free(encoder);

    if (!VerifyMessage(recovered.get(), block_count, block_bytes)) {
        return 1;
    }

    std::cout << "large_message_test: PASS" << std::endl;
    return 0;
}
