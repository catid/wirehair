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

static_assert(sizeof(WirehairV2EncoderOptions) == 16,
    "installed V2 encoder-options ABI size");
static_assert(WIREHAIR_V2_ENCODER_OPTIONS_VERSION == 1u,
    "installed V2 encoder-options version");
static_assert(offsetof(WirehairV2EncoderOptions, struct_bytes) == 0,
    "installed V2 encoder-options struct_bytes offset");
static_assert(offsetof(WirehairV2EncoderOptions, options_version) == 4,
    "installed V2 encoder-options options_version offset");
static_assert(offsetof(WirehairV2EncoderOptions, source_policy) == 8,
    "installed V2 encoder-options source_policy offset");
static_assert(offsetof(WirehairV2EncoderOptions, reserved) == 12,
    "installed V2 encoder-options reserved offset");
static_assert(sizeof(WirehairV2EncoderSourcePolicy) == 4,
    "installed V2 encoder source-policy enum size");
static_assert(WirehairV2EncoderSource_Invalid == 0 &&
        WirehairV2EncoderSource_Independent == 1 &&
        WirehairV2EncoderSource_BorrowedImmutable == 2 &&
        WirehairV2EncoderSource_Count == 3 &&
        WirehairV2EncoderSource_Padding == 0x7fffffff,
    "installed V2 encoder source-policy values");

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

    /**
        Create an encoder that borrows immutable caller-owned message storage.

        The candidate source must be readable, allocated, and byte-for-byte
        immutable on entry and throughout the call.  Failure retains none of
        that candidate source; success extends its obligation until successful
        DetachInput(), successful replacement, or destruction.  A failed
        replacement preserves the current handle and its existing
        source-lifetime obligation.
    */
    Result CreateBorrowed(
        const void* message,
        std::uint64_t messageBytes,
        std::uint32_t blockBytes,
        SerializedProfile& profile) noexcept
    {
        WirehairV2EncoderOptions options =
            WIREHAIR_V2_ENCODER_OPTIONS_INIT;
        options.source_policy = WirehairV2EncoderSource_BorrowedImmutable;
        WirehairV2Codec next = nullptr;
        std::uint32_t written = 0;
        const Result result = wirehair_v2_encoder_create_with_options(
            message, messageBytes, blockBytes, &options,
            profile.data(), profile.size(), &written, &next);
        if (result == WirehairV2_Success) {
            Reset(next);
        }
        return result;
    }

    /** Explicit-profile form of CreateBorrowed(). */
    Result CreateBorrowed(
        std::uint64_t profileId,
        const void* message,
        std::uint64_t messageBytes,
        std::uint32_t blockBytes,
        SerializedProfile& profile) noexcept
    {
        WirehairV2EncoderOptions options =
            WIREHAIR_V2_ENCODER_OPTIONS_INIT;
        options.source_policy = WirehairV2EncoderSource_BorrowedImmutable;
        WirehairV2Codec next = nullptr;
        std::uint32_t written = 0;
        const Result result =
            wirehair_v2_encoder_create_profile_id_with_options(
                profileId, message, messageBytes, blockBytes, &options,
                profile.data(), profile.size(), &written, &next);
        if (result == WirehairV2_Success) {
            Reset(next);
        }
        return result;
    }

    /** Serialized-profile form of CreateBorrowed(). */
    Result CreateBorrowed(
        const void* message,
        const SerializedProfile& profile) noexcept
    {
        WirehairV2EncoderOptions options =
            WIREHAIR_V2_ENCODER_OPTIONS_INIT;
        options.source_policy = WirehairV2EncoderSource_BorrowedImmutable;
        WirehairV2Codec next = nullptr;
        const Result result =
            wirehair_v2_encoder_create_profile_with_options(
                message, profile.data(), profile.size(), &options, &next);
        if (result == WirehairV2_Success) {
            Reset(next);
        }
        return result;
    }

    /** Release a borrowed input while retaining the solved encoder. */
    Result DetachInput() noexcept
    {
        return wirehair_v2_encoder_detach_input(Handle);
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
