#include "roundtrip.h"

#include <wirehair/wirehair.hpp>

#include <array>
#include <cstdint>
#include <cstring>
#include <utility>

namespace {

int CppV2RoundTrip()
{
    static_assert(sizeof(wirehair::v2::Profile) == 32,
        "installed V2 profile ABI size");
    static_assert(sizeof(WirehairV2EncoderOptions) == 16,
        "installed V2 encoder options ABI size");
    static_assert(offsetof(WirehairV2EncoderOptions, struct_bytes) == 0 &&
            offsetof(WirehairV2EncoderOptions, options_version) == 4 &&
            offsetof(WirehairV2EncoderOptions, source_policy) == 8 &&
            offsetof(WirehairV2EncoderOptions, reserved) == 12,
        "installed V2 encoder options ABI offsets");
    static_assert(sizeof(WirehairV2EncoderSourcePolicy) == 4 &&
            WirehairV2EncoderSource_Invalid == 0 &&
            WirehairV2EncoderSource_Independent == 1 &&
            WirehairV2EncoderSource_BorrowedImmutable == 2 &&
            WirehairV2EncoderSource_Count == 3 &&
            WirehairV2EncoderSource_Padding == 0x7fffffff,
        "installed V2 encoder source policy ABI");
    std::array<std::uint8_t, 64> message{};
    for (std::size_t i = 0; i < message.size(); ++i) {
        message[i] = static_cast<std::uint8_t>(i * 17u + 3u);
    }

    wirehair::v2::SerializedProfile profile;
    wirehair::v2::SerializedProfile borrowed_profile;
    wirehair::v2::SerializedProfile borrowed_explicit_profile;
    wirehair::v2::Encoder encoder;
    wirehair::v2::Encoder borrowed;
    wirehair::v2::Encoder borrowed_explicit;
    wirehair::v2::Encoder borrowed_from_profile;
    wirehair::v2::Decoder decoder;
    if (encoder.Create(
            WIREHAIR_V2_PROFILE_CERTIFIED_2026_07,
            message.data(), message.size(), 16u, profile) !=
                WirehairV2_Success ||
        decoder.Create(profile) != WirehairV2_Success)
    {
        return 1;
    }
    if (borrowed.CreateBorrowed(
            message.data(), message.size(), 16u, borrowed_profile) !=
                WirehairV2_Success ||
        borrowed_explicit.CreateBorrowed(
            WIREHAIR_V2_PROFILE_CERTIFIED_2026_07,
            message.data(), message.size(), 16u,
            borrowed_explicit_profile) != WirehairV2_Success ||
        borrowed_from_profile.CreateBorrowed(message.data(), profile) !=
            WirehairV2_Success ||
        std::memcmp(
            borrowed_profile.data(), profile.data(), profile.size()) != 0 ||
        std::memcmp(
            borrowed_explicit_profile.data(),
            profile.data(), profile.size()) != 0)
    {
        return 5;
    }
    {
        std::array<std::uint8_t, 16> systematic{};
        std::uint32_t written = 0;
        if (borrowed.Encode(
                0u, systematic.data(),
                static_cast<std::uint32_t>(systematic.size()), written) !=
                    WirehairV2_Success ||
            written != static_cast<std::uint32_t>(systematic.size()) ||
            std::memcmp(
                systematic.data(), message.data(), systematic.size()) != 0)
        {
            return 6;
        }
    }
    wirehair::v2::Encoder moved(std::move(borrowed));
    if (borrowed || !moved || moved.DetachInput() != WirehairV2_Success ||
        moved.DetachInput() != WirehairV2_Success ||
        borrowed_explicit.DetachInput() != WirehairV2_Success ||
        borrowed_from_profile.DetachInput() != WirehairV2_Success)
    {
        return 7;
    }
    WirehairV2Profile parsed{};
    if (wirehair_v2_profile_deserialize(
            profile.data(), profile.size(), &parsed) != WirehairV2_Success ||
        parsed.profile_id != WIREHAIR_V2_PROFILE_CERTIFIED_2026_07)
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

} // namespace

int main()
{
    const int c_result = wirehair_package_round_trip();
    return c_result == 0 ? CppV2RoundTrip() : c_result;
}
