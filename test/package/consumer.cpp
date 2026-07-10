#include "roundtrip.h"

#include <wirehair/wirehair.hpp>

#include <array>
#include <cstdint>
#include <cstring>

#if __cplusplus < 201103L
#error "wirehair::wirehair must publish the C++11 requirement of wirehair.hpp"
#endif

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
            message.data(), message.size(), 16u, profile) !=
                WirehairV2_Success ||
        decoder.Create(profile) != WirehairV2_Success)
    {
        return 1;
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
