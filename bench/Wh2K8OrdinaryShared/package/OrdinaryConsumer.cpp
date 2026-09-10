// Installed headers and actual ordinary APIs only. No private codec linkage.
#include <wirehair/wirehair.hpp>
#include <array>
#include <cstdint>
#include <cstring>
#include <utility>

#ifndef WIREHAIR_DLL
#error "The installed shared target must publish WIREHAIR_DLL"
#endif
#ifndef WH2_EXPECT_SMALL_K8
#error "The expected ordinary selector must be explicit"
#endif

static bool RoundTrip(unsigned block, bool partial, bool borrowed)
{
    const unsigned message_bytes = 8 * block - (partial ? block - 1 : 0);
    std::array<std::uint8_t, 512> original{};
    for (unsigned i = 0; i < message_bytes; ++i) original[i] = static_cast<std::uint8_t>(i * 17 + 3);
    auto source = original;
    std::array<std::array<std::uint8_t, 64>, 40> packets{};
    wirehair::v2::SerializedProfile profile;
    {
        wirehair::v2::Encoder encoder;
        const auto result = borrowed ? encoder.CreateBorrowed(source.data(), message_bytes, block, profile) :
            encoder.Create(source.data(), message_bytes, block, profile);
        if (result != WirehairV2_Success) return false;
        WirehairV2Profile host{};
        if (wirehair_v2_profile_deserialize(profile.data(), profile.size(), &host) != WirehairV2_Success ||
            host.profile_id != (WH2_EXPECT_SMALL_K8 ? WIREHAIR_V2_PROFILE_SMALL_K8_2026_09 :
                               WIREHAIR_V2_PROFILE_CERTIFIED_2026_07) ||
            host.message_bytes != message_bytes || host.block_bytes != block ||
            (WH2_EXPECT_SMALL_K8 && host.seed_attempt != 0)) return false;
        wirehair::v2::Encoder moved(std::move(encoder));
        if (encoder || !moved || moved.DetachInput() != WirehairV2_Success ||
            moved.DetachInput() != WirehairV2_Success) return false;
        source.fill(0xcc);
        for (unsigned index = 0; index < packets.size(); ++index) {
            std::uint32_t written = 0;
            if (moved.Encode(8 + index, packets[index].data(), block, written) != WirehairV2_Success ||
                written != block) return false;
        }
    } // No sender survives receiver construction.
    wirehair::v2::Decoder decoder;
    if (decoder.Create(profile) != WirehairV2_Success) return false;
    bool succeeded = false;
    for (unsigned index = 0; index < packets.size(); ++index) {
        const auto result = decoder.Decode(8 + index, packets[index].data(), block);
        if (result == WirehairV2_Success) { succeeded = true; break; }
        if (result != WirehairV2_NeedMore) return false;
    }
    if (!succeeded) return false;
    for (unsigned repeat = 0; repeat < 2; ++repeat) {
        std::array<std::uint8_t, 514> output;
        output.fill(0xa5);
        std::uint64_t written = 0;
        if (decoder.Recover(output.data() + 1, message_bytes, written) != WirehairV2_Success ||
            written != message_bytes || std::memcmp(output.data() + 1, original.data(), message_bytes) ||
            output[0] != 0xa5 || output[message_bytes + 1] != 0xa5) return false;
    }
    return true;
}

int main()
{
    if (wirehair_init() != Wirehair_Success) return 1;
    for (unsigned block : {2u, 64u}) for (bool partial : {false, true})
        for (bool borrowed : {false, true}) if (!RoundTrip(block, partial, borrowed)) return 2;
    return 0;
}
