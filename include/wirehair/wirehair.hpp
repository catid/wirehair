#ifndef WIREHAIR_HPP
#define WIREHAIR_HPP

#include <wirehair/wirehair.h>

#include <array>
#include <cstddef>
#include <cstdint>

namespace wirehair {
namespace v2 {

using Result = WirehairV2Result;
using Profile = WirehairV2Profile;

/** Fixed-size owner for one canonical serialized V2 profile. */
class SerializedProfile
{
public:
    std::uint8_t* data() noexcept { return Bytes.data(); }
    const std::uint8_t* data() const noexcept { return Bytes.data(); }
    std::uint32_t size() const noexcept
    {
        return static_cast<std::uint32_t>(Bytes.size());
    }

    Result Deserialize(Profile& profile) const noexcept
    {
        return wirehair_v2_profile_deserialize(
            Bytes.data(), size(), &profile);
    }

    Result Validate() const noexcept
    {
        return wirehair_v2_profile_validate(Bytes.data(), size());
    }

    Result Serialize(const Profile& profile) noexcept
    {
        std::uint32_t written = 0;
        return wirehair_v2_profile_serialize(
            &profile, Bytes.data(), size(), &written);
    }

private:
    std::array<std::uint8_t,
        WIREHAIR_V2_PROFILE_SERIALIZED_BYTES> Bytes{};
};

/** Move-only RAII wrapper for a public V2 encoder. */
class Encoder
{
public:
    Encoder() noexcept = default;
    ~Encoder() { wirehair_v2_free(Handle); }

    Encoder(const Encoder&) = delete;
    Encoder& operator=(const Encoder&) = delete;

    Encoder(Encoder&& other) noexcept : Handle(other.Handle)
    {
        other.Handle = nullptr;
    }

    Encoder& operator=(Encoder&& other) noexcept
    {
        if (this != &other) {
            wirehair_v2_free(Handle);
            Handle = other.Handle;
            other.Handle = nullptr;
        }
        return *this;
    }

    Result Create(
        const void* message,
        std::uint64_t messageBytes,
        std::uint32_t blockBytes,
        SerializedProfile& profile) noexcept
    {
        WirehairV2Codec next = nullptr;
        std::uint32_t written = 0;
        const Result result = wirehair_v2_encoder_create(
            message, messageBytes, blockBytes,
            profile.data(), profile.size(), &written, &next);
        if (result == WirehairV2_Success) {
            Reset(next);
        }
        return result;
    }

    Result Create(
        std::uint64_t profileId,
        const void* message,
        std::uint64_t messageBytes,
        std::uint32_t blockBytes,
        SerializedProfile& profile) noexcept
    {
        WirehairV2Codec next = nullptr;
        std::uint32_t written = 0;
        const Result result = wirehair_v2_encoder_create_profile_id(
            profileId, message, messageBytes, blockBytes,
            profile.data(), profile.size(), &written, &next);
        if (result == WirehairV2_Success) {
            Reset(next);
        }
        return result;
    }

    Result Create(
        const void* message,
        const SerializedProfile& profile) noexcept
    {
        WirehairV2Codec next = nullptr;
        const Result result = wirehair_v2_encoder_create_profile(
            message, profile.data(), profile.size(), &next);
        if (result == WirehairV2_Success) {
            Reset(next);
        }
        return result;
    }

    Result Encode(
        std::uint32_t blockId,
        void* blockOut,
        std::uint32_t outputCapacity,
        std::uint32_t& dataBytesOut) noexcept
    {
        return wirehair_v2_encode(
            Handle, blockId, blockOut, outputCapacity, &dataBytesOut);
    }

    bool valid() const noexcept { return Handle != nullptr; }
    explicit operator bool() const noexcept { return valid(); }

private:
    void Reset(WirehairV2Codec next) noexcept
    {
        wirehair_v2_free(Handle);
        Handle = next;
    }

    WirehairV2Codec Handle = nullptr;
};

/** Move-only RAII wrapper for a decoder created only from a descriptor. */
class Decoder
{
public:
    Decoder() noexcept = default;
    ~Decoder() { wirehair_v2_free(Handle); }

    Decoder(const Decoder&) = delete;
    Decoder& operator=(const Decoder&) = delete;

    Decoder(Decoder&& other) noexcept : Handle(other.Handle)
    {
        other.Handle = nullptr;
    }

    Decoder& operator=(Decoder&& other) noexcept
    {
        if (this != &other) {
            wirehair_v2_free(Handle);
            Handle = other.Handle;
            other.Handle = nullptr;
        }
        return *this;
    }

    Result Create(const SerializedProfile& profile) noexcept
    {
        WirehairV2Codec next = nullptr;
        const Result result = wirehair_v2_decoder_create(
            profile.data(), profile.size(), &next);
        if (result == WirehairV2_Success) {
            Reset(next);
        }
        return result;
    }

    Result Decode(
        std::uint32_t blockId,
        const void* blockData,
        std::uint32_t dataBytes) noexcept
    {
        return wirehair_v2_decode(
            Handle, blockId, blockData, dataBytes);
    }

    Result Recover(
        void* messageOut,
        std::uint64_t outputCapacity,
        std::uint64_t& bytesOut) noexcept
    {
        return wirehair_v2_recover(
            Handle, messageOut, outputCapacity, &bytesOut);
    }

    bool valid() const noexcept { return Handle != nullptr; }
    explicit operator bool() const noexcept { return valid(); }

private:
    void Reset(WirehairV2Codec next) noexcept
    {
        wirehair_v2_free(Handle);
        Handle = next;
    }

    WirehairV2Codec Handle = nullptr;
};

} // namespace v2
} // namespace wirehair

#endif // WIREHAIR_HPP
