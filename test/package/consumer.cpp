#include "roundtrip.h"

#include <wirehair/wirehair.hpp>

#include <array>
#include <cstdint>
#include <cstring>

namespace {

int CppV2RoundTrip()
{
    static_assert(sizeof(wirehair::v2::Profile) == 32,
        "installed V2 profile ABI size");
    std::array<std::uint8_t, 64> message{};
    for (std::size_t i = 0; i < message.size(); ++i) {
        message[i] = static_cast<std::uint8_t>(i * 17u + 3u);
    }

    wirehair::v2::SerializedProfile profile;
    wirehair::v2::Encoder encoder;
    wirehair::v2::Decoder decoder;
    if (encoder.Create(
            WIREHAIR_V2_PROFILE_MIXED_2026_07,
            message.data(), message.size(), 16u, profile) !=
                WirehairV2_Success ||
        decoder.Create(profile) != WirehairV2_Success)
    {
        return 1;
    }
    WirehairV2Profile parsed{};
    if (wirehair_v2_profile_deserialize(
            profile.data(), profile.size(), &parsed) != WirehairV2_Success ||
        parsed.profile_id != WIREHAIR_V2_PROFILE_MIXED_2026_07)
    {
        return 4;
    }

    std::array<std::uint8_t, 16> block{};
    WirehairV2Result result = WirehairV2_NeedMore;
    for (std::uint32_t id = 4u;
         id < 68u && result == WirehairV2_NeedMore;
         ++id)
    {
        std::uint32_t bytes = 0;
        if (encoder.Encode(
                id, block.data(), static_cast<std::uint32_t>(block.size()),
                bytes) !=
                WirehairV2_Success)
        {
            return 2;
        }
        result = decoder.Decode(id, block.data(), bytes);
    }
    std::array<std::uint8_t, 64> recovered{};
    std::uint64_t recovered_bytes = 0;
    if (result != WirehairV2_Success ||
        decoder.Recover(
            recovered.data(), recovered.size(), recovered_bytes) !=
                WirehairV2_Success ||
        recovered_bytes != message.size() || recovered != message)
    {
        return 3;
    }
    return 0;
}

int CppTinyMdsRoundTrip()
{
    std::array<std::uint8_t, 31> message{};
    for (std::size_t i = 0; i < message.size(); ++i) {
        message[i] = static_cast<std::uint8_t>(i * 29u + 11u);
    }

    wirehair::v2::SerializedProfile profile;
    wirehair::v2::Encoder encoder;
    wirehair::v2::Decoder decoder;
    if (encoder.Create(
            WIREHAIR_V2_PROFILE_TINY_MDS_2026_07,
            message.data(), message.size(), 16u, profile) !=
                WirehairV2_Success ||
        decoder.Create(profile) != WirehairV2_Success)
    {
        return 5;
    }

    std::array<std::uint8_t, 16> block{};
    std::uint32_t bytes = 0u;
    if (encoder.Encode(
            2u, block.data(), static_cast<std::uint32_t>(block.size()),
            bytes) != WirehairV2_Success ||
        decoder.Decode(2u, block.data(), bytes) != WirehairV2_NeedMore ||
        encoder.Encode(
            256u, block.data(), static_cast<std::uint32_t>(block.size()),
            bytes) != WirehairV2_Success ||
        decoder.Decode(256u, block.data(), bytes) != WirehairV2_Success)
    {
        return 6;
    }

    std::array<std::uint8_t, 31> recovered{};
    std::uint64_t recovered_bytes = 0u;
    if (decoder.Recover(
            recovered.data(), recovered.size(), recovered_bytes) !=
                WirehairV2_Success ||
        recovered_bytes != message.size() || recovered != message)
    {
        return 7;
    }
    return 0;
}

} // namespace

int main()
{
    const int c_result = wirehair_package_round_trip();
    if (c_result != 0) {
        return c_result;
    }
    const int v2_result = CppV2RoundTrip();
    return v2_result == 0 ? CppTinyMdsRoundTrip() : v2_result;
}
