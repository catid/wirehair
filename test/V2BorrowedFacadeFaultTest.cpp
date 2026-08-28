#include <wirehair/wirehair.h>
#include <wirehair/wirehair.hpp>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace {

int Failures = 0;

#define CHECK(expression) do { \
    if (!(expression)) { \
        std::fprintf(stderr, "CHECK failed at %s:%d: %s\n", \
            __FILE__, __LINE__, #expression); \
        ++Failures; \
    } \
} while (0)

bool SetEnv(const char* name, const char* value)
{
#if defined(_WIN32)
    return _putenv_s(name, value ? value : "") == 0;
#else
    if (value) {
        return setenv(name, value, 1) == 0;
    }
    return unsetenv(name) == 0;
#endif
}

void FillMessage(std::vector<uint8_t>& message, uint64_t state)
{
    for (size_t i = 0; i < message.size(); ++i)
    {
        state ^= state >> 12;
        state ^= state << 25;
        state ^= state >> 27;
        message[i] = static_cast<uint8_t>(
            state * UINT64_C(2685821657736338717));
    }
}

bool EncodePacket(
    wirehair::v2::Encoder& encoder,
    uint32_t id,
    uint32_t block_bytes,
    std::vector<uint8_t>& packet)
{
    packet.assign(block_bytes, 0xa5u);
    uint32_t written = 0u;
    const WirehairV2Result result = encoder.Encode(
        id, packet.data(), block_bytes, written);
    if (result != WirehairV2_Success || written == 0u ||
        written > block_bytes)
    {
        return false;
    }
    packet.resize(written);
    return true;
}

void CheckSelectingRouteOom(
    bool explicit_profile,
    const std::vector<uint8_t>& message,
    uint32_t block_bytes)
{
    WirehairV2EncoderOptions options = WIREHAIR_V2_ENCODER_OPTIONS_INIT;
    options.source_policy = WirehairV2EncoderSource_BorrowedImmutable;
    std::array<uint8_t, WIREHAIR_V2_PROFILE_SERIALIZED_BYTES> profile;
    profile.fill(0x5au);
    const std::array<uint8_t, WIREHAIR_V2_PROFILE_SERIALIZED_BYTES> before =
        profile;
    uint32_t written = UINT32_MAX;
    WirehairV2Codec codec = reinterpret_cast<WirehairV2Codec>(
        static_cast<uintptr_t>(1u));
    const WirehairV2Result result = explicit_profile ?
        wirehair_v2_encoder_create_profile_id_with_options(
            WIREHAIR_V2_PROFILE_CURRENT,
            message.data(), message.size(), block_bytes, &options,
            profile.data(), static_cast<uint32_t>(profile.size()),
            &written, &codec) :
        wirehair_v2_encoder_create_with_options(
            message.data(), message.size(), block_bytes, &options,
            profile.data(), static_cast<uint32_t>(profile.size()),
            &written, &codec);
    CHECK(result == WirehairV2_OOM);
    CHECK(codec == nullptr);
    CHECK(written == WIREHAIR_V2_PROFILE_SERIALIZED_BYTES);
    CHECK(profile == before);
}

} // namespace

int main()
{
    static const char* const OomVariable =
        "WIREHAIR_TEST_FORCE_V2_FACADE_OOM";
    enum : uint32_t { BlockBytes = 16u, BlockCount = 8u };

    CHECK(SetEnv(OomVariable, nullptr));
    std::vector<uint8_t> original(BlockBytes * BlockCount - 3u);
    FillMessage(original, UINT64_C(0x243f6a8885a308d3));

    wirehair::v2::SerializedProfile profile;
    wirehair::v2::Encoder encoder;
    CHECK(encoder.CreateBorrowed(
        original.data(), original.size(), BlockBytes, profile) ==
        WirehairV2_Success);
    std::vector<uint8_t> expected_systematic;
    std::vector<uint8_t> expected_repair;
    CHECK(EncodePacket(encoder, 0u, BlockBytes, expected_systematic));
    CHECK(EncodePacket(
        encoder, BlockCount + 7u, BlockBytes, expected_repair));

    CHECK(SetEnv(OomVariable, "1"));
    CheckSelectingRouteOom(false, original, BlockBytes);
    CheckSelectingRouteOom(true, original, BlockBytes);

    WirehairV2EncoderOptions options = WIREHAIR_V2_ENCODER_OPTIONS_INIT;
    options.source_policy = WirehairV2EncoderSource_BorrowedImmutable;
    WirehairV2Codec raw = reinterpret_cast<WirehairV2Codec>(
        static_cast<uintptr_t>(1u));
    CHECK(wirehair_v2_encoder_create_profile_with_options(
        original.data(), profile.data(), profile.size(), &options, &raw) ==
        WirehairV2_OOM);
    CHECK(raw == nullptr);

    std::vector<uint8_t> replacement(original.size());
    FillMessage(replacement, UINT64_C(0x13198a2e03707344));
    wirehair::v2::SerializedProfile replacement_profile;
    CHECK(encoder.CreateBorrowed(
        replacement.data(), replacement.size(), BlockBytes,
        replacement_profile) == WirehairV2_OOM);
    CHECK(encoder);

    std::vector<uint8_t> actual;
    CHECK(EncodePacket(encoder, 0u, BlockBytes, actual));
    CHECK(actual == expected_systematic);
    CHECK(EncodePacket(
        encoder, BlockCount + 7u, BlockBytes, actual));
    CHECK(actual == expected_repair);

    CHECK(encoder.DetachInput() == WirehairV2_Success);
    std::fill(original.begin(), original.end(), 0xccu);
    std::vector<uint8_t>().swap(original);
    CHECK(SetEnv(OomVariable, nullptr));
    CHECK(EncodePacket(encoder, 0u, BlockBytes, actual));
    CHECK(actual == expected_systematic);
    CHECK(EncodePacket(
        encoder, BlockCount + 7u, BlockBytes, actual));
    CHECK(actual == expected_repair);

    if (Failures != 0) {
        std::fprintf(stderr, "V2 borrowed facade fault test: %d failures\n",
            Failures);
        return 1;
    }
    std::puts("V2 borrowed facade fault test passed");
    return 0;
}
