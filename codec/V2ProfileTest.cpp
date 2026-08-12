#include <wirehair/wirehair.h>
#include <wirehair/wirehair.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <utility>
#include <vector>

namespace {

static_assert(sizeof(WirehairV2Profile) == 32,
    "public V2 profile ABI size changed");
static_assert(offsetof(WirehairV2Profile, profile_id) == 8,
    "public V2 profile ABI layout changed");
static_assert(offsetof(WirehairV2Profile, message_bytes) == 16,
    "public V2 profile ABI layout changed");
static_assert(offsetof(WirehairV2Profile, block_bytes) == 24,
    "public V2 profile ABI layout changed");
static_assert(offsetof(WirehairV2Profile, seed_attempt) == 28,
    "public V2 profile ABI layout changed");
static_assert(WirehairV2_UnsupportedPlatform == 14,
    "public V2 result values are stable ABI");

bool Check(bool condition, const char* what)
{
    if (condition) {
        return true;
    }
    std::fprintf(stderr, "V2 profile test failed: %s\n", what);
    return false;
}

bool IsZeroProfile(const WirehairV2Profile& profile)
{
    const uint8_t zero[sizeof(profile)] = {};
    return std::memcmp(&profile, zero, sizeof(profile)) == 0;
}

void FillMessage(std::vector<uint8_t>& message)
{
    for (size_t i = 0; i < message.size(); ++i) {
        message[i] = (uint8_t)(i * 73u + 19u);
    }
}

WirehairV2Result CallSelectingEncoderConstructor(
    bool explicit_profile,
    const void* message,
    uint64_t message_bytes,
    uint32_t block_bytes,
    void* serialized_profile_out,
    uint32_t serialized_profile_capacity,
    uint32_t* serialized_profile_bytes_out,
    WirehairV2Codec* codec_out)
{
    if (explicit_profile)
    {
        return wirehair_v2_encoder_create_profile_id(
            WIREHAIR_V2_PROFILE_CURRENT,
            message, message_bytes, block_bytes,
            serialized_profile_out, serialized_profile_capacity,
            serialized_profile_bytes_out, codec_out);
    }
    return wirehair_v2_encoder_create(
        message, message_bytes, block_bytes,
        serialized_profile_out, serialized_profile_capacity,
        serialized_profile_bytes_out, codec_out);
}

bool CheckSelectingConstructorOverlapGuards(
    const std::vector<uint8_t>& disjoint_message)
{
    static const uint32_t kProfileBytes =
        WIREHAIR_V2_PROFILE_SERIALIZED_BYTES;
    enum OverlapCase
    {
        DescriptorSize,
        DescriptorCodec,
        SizeCodec,
        SizeMessage,
        CodecMessage,
        OverlapCaseCount
    };
    enum RequestCase
    {
        FullRequest,
        ShortRequest,
        InvalidRequest,
        RequestCaseCount
    };

    const std::vector<uint8_t> message_before = disjoint_message;
    for (unsigned explicit_profile = 0u;
         explicit_profile < 2u;
         ++explicit_profile)
    {
        for (unsigned request_case = 0u;
             request_case < RequestCaseCount;
             ++request_case)
        {
            for (unsigned overlap_case = 0u;
                 overlap_case < OverlapCaseCount;
                 ++overlap_case)
            {
                for (unsigned partial = 0u; partial < 2u; ++partial)
                {
                    if (overlap_case == SizeCodec && partial != 0u &&
                        sizeof(WirehairV2Codec) == sizeof(uint32_t))
                    {
                        // Two naturally aligned four-byte output objects can
                        // only have identical or disjoint ranges.
                        continue;
                    }
                    alignas(std::max_align_t) uint8_t storage[512];
                    std::memset(storage, 0xa5, sizeof(storage));
                    uint8_t storage_before[sizeof(storage)];
                    std::memcpy(
                        storage_before, storage, sizeof(storage_before));

                    const void* message = disjoint_message.data();
                    void* descriptor = storage + 256u;
                    uint32_t* size_out =
                        reinterpret_cast<uint32_t*>(storage + 320u);
                    WirehairV2Codec* codec_out =
                        reinterpret_cast<WirehairV2Codec*>(storage + 336u);

                    switch ((OverlapCase)overlap_case)
                    {
                    case DescriptorSize:
                        descriptor = storage + (partial ? 66u : 64u);
                        size_out = reinterpret_cast<uint32_t*>(storage + 64u);
                        break;
                    case DescriptorCodec:
                        descriptor = storage + (partial ? 66u : 64u);
                        codec_out = reinterpret_cast<WirehairV2Codec*>(
                            storage + 64u);
                        break;
                    case SizeCodec:
                        size_out = reinterpret_cast<uint32_t*>(
                            storage + (partial ? 68u : 64u));
                        codec_out = reinterpret_cast<WirehairV2Codec*>(
                            storage + 64u);
                        break;
                    case SizeMessage:
                        message = storage + (partial ? 66u : 64u);
                        size_out = reinterpret_cast<uint32_t*>(storage + 64u);
                        break;
                    case CodecMessage:
                        message = storage + (partial ? 66u : 64u);
                        codec_out = reinterpret_cast<WirehairV2Codec*>(
                            storage + 64u);
                        break;
                    default:
                        return false;
                    }

                    const uint32_t capacity =
                        request_case == ShortRequest ?
                            kProfileBytes - 1u : kProfileBytes;
                    const uint32_t block_bytes =
                        request_case == InvalidRequest ? 0u : 16u;
                    const WirehairV2Result result =
                        CallSelectingEncoderConstructor(
                            explicit_profile != 0u,
                            message, disjoint_message.size(), block_bytes,
                            descriptor, capacity, size_out, codec_out);
                    if (result != WirehairV2_InvalidInput ||
                        std::memcmp(
                            storage, storage_before, sizeof(storage)) != 0 ||
                        disjoint_message != message_before)
                    {
                        std::fprintf(stderr,
                            "V2 constructor overlap changed storage: "
                            "explicit=%u request=%u overlap=%u partial=%u "
                            "result=%d\n",
                            explicit_profile, request_case, overlap_case,
                            partial, (int)result);
                        return false;
                    }
                }
            }
        }
    }

    for (unsigned explicit_profile = 0u;
         explicit_profile < 2u;
         ++explicit_profile)
    {
        uint8_t descriptor[kProfileBytes];
        std::memset(descriptor, 0x5a, sizeof(descriptor));
        uint8_t descriptor_before[sizeof(descriptor)];
        std::memcpy(
            descriptor_before, descriptor, sizeof(descriptor_before));
        uint32_t profile_bytes = UINT32_C(0xa5a5a5a5);
        const WirehairV2Codec codec_sentinel =
            reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
        WirehairV2Codec codec = codec_sentinel;
        const WirehairV2Result result = CallSelectingEncoderConstructor(
            explicit_profile != 0u,
            disjoint_message.data(), UINT64_MAX, 16u,
            descriptor, sizeof(descriptor), &profile_bytes, &codec);
        if (result != WirehairV2_InvalidDimensions ||
            profile_bytes != kProfileBytes || codec != nullptr ||
            std::memcmp(
                descriptor, descriptor_before, sizeof(descriptor)) != 0 ||
            disjoint_message != message_before)
        {
            std::fprintf(stderr,
                "V2 constructor UINT64_MAX precedence changed: "
                "explicit=%u result=%d\n",
                explicit_profile, (int)result);
            if (result == WirehairV2_Success &&
                codec != codec_sentinel)
            {
                wirehair_v2_free(codec);
            }
            return false;
        }

        std::memset(descriptor, 0x5a, sizeof(descriptor));
        std::memcpy(
            descriptor_before, descriptor, sizeof(descriptor_before));
        codec = codec_sentinel;
        const WirehairV2Result missing_size =
            CallSelectingEncoderConstructor(
                explicit_profile != 0u,
                disjoint_message.data(), disjoint_message.size(), 16u,
                descriptor, sizeof(descriptor), nullptr, &codec);
        if (missing_size != WirehairV2_InvalidInput || codec != nullptr ||
            std::memcmp(
                descriptor, descriptor_before, sizeof(descriptor)) != 0 ||
            disjoint_message != message_before)
        {
            std::fprintf(stderr,
                "V2 constructor missing size-output behavior changed: "
                "explicit=%u result=%d\n",
                explicit_profile, (int)missing_size);
            if (missing_size == WirehairV2_Success &&
                codec != codec_sentinel)
            {
                wirehair_v2_free(codec);
            }
            return false;
        }

        std::memset(descriptor, 0x5a, sizeof(descriptor));
        std::memcpy(
            descriptor_before, descriptor, sizeof(descriptor_before));
        profile_bytes = UINT32_C(0xa5a5a5a5);
        const WirehairV2Result missing_codec =
            CallSelectingEncoderConstructor(
                explicit_profile != 0u,
                disjoint_message.data(), disjoint_message.size(), 16u,
                descriptor, sizeof(descriptor), &profile_bytes, nullptr);
        if (missing_codec != WirehairV2_InvalidInput ||
            profile_bytes != kProfileBytes ||
            std::memcmp(
                descriptor, descriptor_before, sizeof(descriptor)) != 0 ||
            disjoint_message != message_before)
        {
            std::fprintf(stderr,
                "V2 constructor missing codec-output behavior changed: "
                "explicit=%u result=%d\n",
                explicit_profile, (int)missing_codec);
            return false;
        }
    }
    return true;
}

bool CheckStagedDescriptorMessageAliases(
    const std::vector<uint8_t>& message,
    const uint8_t* expected_profile)
{
    static const uint32_t kProfileBytes =
        WIREHAIR_V2_PROFILE_SERIALIZED_BYTES;
    for (unsigned explicit_profile = 0u;
         explicit_profile < 2u;
         ++explicit_profile)
    {
        for (unsigned boundary = 0u; boundary < 2u; ++boundary)
        {
            for (unsigned short_output = 0u;
                 short_output < 2u;
                 ++short_output)
            {
                alignas(std::max_align_t) uint8_t storage[256];
                std::memset(storage, 0x5a, sizeof(storage));
                uint8_t* const message_in = storage + 64u;
                std::memcpy(message_in, message.data(), message.size());
                void* const descriptor_out = boundary ?
                    message_in + message.size() - 2u : message_in;
                uint8_t storage_before[sizeof(storage)];
                std::memcpy(storage_before, storage, sizeof(storage_before));
                uint32_t profile_bytes = UINT32_C(0xa5a5a5a5);
                const WirehairV2Codec codec_sentinel =
                    reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
                WirehairV2Codec codec = codec_sentinel;

                const uint32_t capacity = short_output ?
                    kProfileBytes - 1u : kProfileBytes;
                const WirehairV2Result result =
                    CallSelectingEncoderConstructor(
                        explicit_profile != 0u,
                        message_in, message.size(), 16u,
                        descriptor_out, capacity, &profile_bytes, &codec);
                if (short_output)
                {
                    if (result != WirehairV2_BufferTooSmall ||
                        profile_bytes != kProfileBytes || codec != nullptr ||
                        std::memcmp(
                            storage, storage_before, sizeof(storage)) != 0)
                    {
                        std::fprintf(stderr,
                            "V2 staged short descriptor/message alias failed: "
                            "explicit=%u boundary=%u result=%d\n",
                            explicit_profile, boundary, (int)result);
                        if (result == WirehairV2_Success &&
                            codec != codec_sentinel)
                        {
                            wirehair_v2_free(codec);
                        }
                        return false;
                    }
                    continue;
                }

                WirehairV2Result encode_result = WirehairV2_Error;
                bool payload_ok =
                    result == WirehairV2_Success && codec &&
                    codec != codec_sentinel;
                for (uint32_t block_id = 0u;
                     payload_ok && block_id < 8u;
                     ++block_id)
                {
                    const uint32_t expected_bytes =
                        block_id == 7u ? 5u : 16u;
                    uint8_t encoded[16] = {};
                    uint32_t encoded_bytes = 0u;
                    encode_result = wirehair_v2_encode(
                        codec, block_id, encoded, sizeof(encoded),
                        &encoded_bytes);
                    payload_ok = encode_result == WirehairV2_Success &&
                        encoded_bytes == expected_bytes &&
                        std::memcmp(
                            encoded,
                            message.data() + (size_t)block_id * 16u,
                            expected_bytes) == 0;
                }
                const bool ok =
                    result == WirehairV2_Success &&
                    profile_bytes == kProfileBytes && codec != nullptr &&
                    std::memcmp(
                        descriptor_out, expected_profile, kProfileBytes) == 0 &&
                    payload_ok;
                if (codec != codec_sentinel) {
                    wirehair_v2_free(codec);
                }
                if (!ok)
                {
                    std::fprintf(stderr,
                        "V2 staged descriptor/message alias failed: "
                        "explicit=%u boundary=%u result=%d encode=%d\n",
                        explicit_profile, boundary,
                        (int)result, (int)encode_result);
                    return false;
                }
            }
        }
    }
    return true;
}

bool CheckDescriptorConstructorOverlapGuards(
    const std::vector<uint8_t>& message,
    const uint8_t* canonical_profile)
{
    static const uint32_t kProfileBytes =
        WIREHAIR_V2_PROFILE_SERIALIZED_BYTES;
    enum DescriptorRequest
    {
        FullDescriptor,
        ShortDescriptor,
        InvalidDescriptor,
        DescriptorRequestCount
    };

    const std::vector<uint8_t> message_before = message;
    for (unsigned decoder = 0u; decoder < 2u; ++decoder)
    {
        for (unsigned request = 0u;
             request < DescriptorRequestCount;
             ++request)
        {
            for (unsigned partial = 0u; partial < 2u; ++partial)
            {
                alignas(std::max_align_t) uint8_t storage[128];
                std::memset(storage, 0xa5, sizeof(storage));
                uint8_t* const descriptor =
                    storage + (partial ? 66u : 64u);
                std::memcpy(
                    descriptor, canonical_profile, kProfileBytes);
                if (request == InvalidDescriptor) {
                    descriptor[kProfileBytes - 1u] = 1u;
                }
                WirehairV2Codec* const codec_out =
                    reinterpret_cast<WirehairV2Codec*>(storage + 64u);
                uint8_t storage_before[sizeof(storage)];
                std::memcpy(
                    storage_before, storage, sizeof(storage_before));
                const uint32_t descriptor_bytes =
                    request == ShortDescriptor ?
                        kProfileBytes - 1u : kProfileBytes;
                const WirehairV2Result result = decoder ?
                    wirehair_v2_decoder_create(
                        descriptor, descriptor_bytes, codec_out) :
                    wirehair_v2_encoder_create_profile(
                        message.data(), descriptor, descriptor_bytes,
                        codec_out);
                if (result != WirehairV2_InvalidInput ||
                    std::memcmp(
                        storage, storage_before, sizeof(storage)) != 0 ||
                    message != message_before)
                {
                    std::fprintf(stderr,
                        "V2 descriptor/codec alias changed storage: "
                        "decoder=%u request=%u partial=%u result=%d\n",
                        decoder, request, partial, (int)result);
                    return false;
                }
            }
        }
    }

    for (unsigned partial = 0u; partial < 2u; ++partial)
    {
        alignas(std::max_align_t) uint8_t storage[256];
        std::memset(storage, 0xa5, sizeof(storage));
        uint8_t* const message_in = storage + (partial ? 66u : 64u);
        std::memcpy(message_in, message.data(), message.size());
        WirehairV2Codec* const codec_out =
            reinterpret_cast<WirehairV2Codec*>(storage + 64u);
        uint8_t storage_before[sizeof(storage)];
        std::memcpy(storage_before, storage, sizeof(storage_before));
        const WirehairV2Result result = wirehair_v2_encoder_create_profile(
            message_in, canonical_profile, kProfileBytes, codec_out);
        if (result != WirehairV2_InvalidInput ||
            std::memcmp(storage, storage_before, sizeof(storage)) != 0)
        {
            std::fprintf(stderr,
                "V2 profile encoder message/codec alias changed storage: "
                "partial=%u result=%d\n",
                partial, (int)result);
            return false;
        }
    }

    for (unsigned decoder = 0u; decoder < 2u; ++decoder)
    {
        for (unsigned request = ShortDescriptor;
             request <= InvalidDescriptor;
             ++request)
        {
            uint8_t descriptor[kProfileBytes];
            std::memcpy(
                descriptor, canonical_profile, sizeof(descriptor));
            if (request == InvalidDescriptor) {
                descriptor[kProfileBytes - 1u] = 1u;
            }
            uint8_t descriptor_before[sizeof(descriptor)];
            std::memcpy(
                descriptor_before, descriptor, sizeof(descriptor_before));
            const WirehairV2Codec codec_sentinel =
                reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
            WirehairV2Codec codec = codec_sentinel;
            const uint32_t descriptor_bytes =
                request == ShortDescriptor ?
                    kProfileBytes - 1u : kProfileBytes;
            const WirehairV2Result expected =
                request == ShortDescriptor ?
                    WirehairV2_InvalidSize : WirehairV2_ReservedNonzero;
            const WirehairV2Result result = decoder ?
                wirehair_v2_decoder_create(
                    descriptor, descriptor_bytes, &codec) :
                wirehair_v2_encoder_create_profile(
                    message.data(), descriptor, descriptor_bytes, &codec);
            if (result != expected || codec != nullptr ||
                std::memcmp(
                    descriptor, descriptor_before,
                    sizeof(descriptor)) != 0 ||
                message != message_before)
            {
                std::fprintf(stderr,
                    "V2 disjoint descriptor failure changed: "
                    "decoder=%u request=%u result=%d\n",
                    decoder, request, (int)result);
                if (result == WirehairV2_Success &&
                    codec != codec_sentinel)
                {
                    wirehair_v2_free(codec);
                }
                return false;
            }
        }

        const WirehairV2Codec codec_sentinel =
            reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
        WirehairV2Codec codec = codec_sentinel;
        const WirehairV2Result result = decoder ?
            wirehair_v2_decoder_create(
                nullptr, kProfileBytes, &codec) :
            wirehair_v2_encoder_create_profile(
                message.data(), nullptr, kProfileBytes, &codec);
        if (result != WirehairV2_InvalidInput || codec != nullptr ||
            message != message_before)
        {
            std::fprintf(stderr,
                "V2 null descriptor behavior changed: "
                "decoder=%u result=%d\n",
                decoder, (int)result);
            if (result == WirehairV2_Success &&
                codec != codec_sentinel)
            {
                wirehair_v2_free(codec);
            }
            return false;
        }
    }
    return true;
}

bool CheckMalformedProfiles(const uint8_t* golden)
{
    WirehairV2Profile parsed;
    std::memset(&parsed, 0xa5, sizeof(parsed));
    if (!Check(wirehair_v2_profile_deserialize(
            nullptr, WIREHAIR_V2_PROFILE_SERIALIZED_BYTES, &parsed) ==
            WirehairV2_InvalidInput && IsZeroProfile(parsed),
            "null serialized input") ||
        !Check(wirehair_v2_profile_deserialize(
            golden, WIREHAIR_V2_PROFILE_SERIALIZED_BYTES, nullptr) ==
            WirehairV2_InvalidInput,
            "null profile output"))
    {
        return false;
    }

    for (uint32_t bytes = 0;
        bytes < WIREHAIR_V2_PROFILE_SERIALIZED_BYTES;
        ++bytes)
    {
        std::memset(&parsed, 0xa5, sizeof(parsed));
        if (!Check(wirehair_v2_profile_deserialize(
                golden, bytes, &parsed) == WirehairV2_InvalidSize &&
                IsZeroProfile(parsed),
                "truncated profile"))
        {
            return false;
        }
    }
    uint8_t overlong[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES + 1u] = {};
    std::memcpy(
        overlong, golden, WIREHAIR_V2_PROFILE_SERIALIZED_BYTES);
    if (!Check(wirehair_v2_profile_deserialize(
            overlong, sizeof(overlong), &parsed) ==
            WirehairV2_InvalidSize,
            "overlong profile"))
    {
        return false;
    }

    struct Mutation
    {
        uint32_t Offset;
        uint8_t Value;
        WirehairV2Result Expected;
        const char* Name;
    };
    const Mutation mutations[] = {
        {0u, 0u, WirehairV2_InvalidMagic, "magic"},
        {4u, 2u, WirehairV2_UnsupportedVersion, "version"},
        {6u, 31u, WirehairV2_InvalidSize, "declared size"},
        {8u, 0u, WirehairV2_UnsupportedProfile, "profile id"},
        {16u, 0u, WirehairV2_InvalidDimensions, "message bytes"},
        {24u, 0u, WirehairV2_InvalidDimensions, "block bytes"},
        {29u, 1u, WirehairV2_ReservedNonzero, "reserved byte 0"},
        {30u, 1u, WirehairV2_ReservedNonzero, "reserved byte 1"},
        {31u, 1u, WirehairV2_ReservedNonzero, "reserved byte 2"}
    };
    for (const Mutation& mutation : mutations)
    {
        uint8_t malformed[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES];
        std::memcpy(malformed, golden, sizeof(malformed));
        malformed[mutation.Offset] = mutation.Value;
        std::memset(&parsed, 0xa5, sizeof(parsed));
        if (!Check(wirehair_v2_profile_deserialize(
                malformed, sizeof(malformed), &parsed) ==
                    mutation.Expected && IsZeroProfile(parsed),
                mutation.Name))
        {
            return false;
        }
    }

    uint8_t too_many_blocks[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES];
    std::memcpy(too_many_blocks, golden, sizeof(too_many_blocks));
    // 64001 little-endian message bytes with block_bytes == 1.
    too_many_blocks[16] = 0x01u;
    too_many_blocks[17] = 0xfau;
    too_many_blocks[18] = 0x00u;
    too_many_blocks[19] = 0x00u;
    too_many_blocks[24] = 1u;
    if (!Check(wirehair_v2_profile_deserialize(
            too_many_blocks, sizeof(too_many_blocks), &parsed) ==
            WirehairV2_InvalidDimensions,
            "too many blocks"))
    {
        return false;
    }
    uint8_t one_block[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES];
    std::memcpy(one_block, golden, sizeof(one_block));
    one_block[16] = 1u;
    if (!Check(wirehair_v2_profile_deserialize(
            one_block, sizeof(one_block), &parsed) ==
            WirehairV2_InvalidDimensions,
            "too few blocks"))
    {
        return false;
    }
    uint8_t oversized_block[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES];
    std::memcpy(oversized_block, golden, sizeof(oversized_block));
    oversized_block[24] = 0u;
    oversized_block[25] = 0u;
    oversized_block[26] = 0u;
    oversized_block[27] = 0x80u;
    if (!Check(wirehair_v2_profile_deserialize(
            oversized_block, sizeof(oversized_block), &parsed) ==
            WirehairV2_InvalidDimensions,
            "block size exceeds signed GF256 domain"))
    {
        return false;
    }

    WirehairV2Profile profile = {};
    if (!Check(wirehair_v2_profile_deserialize(
            golden, WIREHAIR_V2_PROFILE_SERIALIZED_BYTES, &profile) ==
            WirehairV2_Success,
            "decode profile for host validation"))
    {
        return false;
    }
    uint8_t untouched[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES];
    std::memset(untouched, 0x5a, sizeof(untouched));
    uint32_t required = 0;
    if (!Check(wirehair_v2_profile_serialize(
            &profile, nullptr, 0u, &required) ==
                WirehairV2_BufferTooSmall &&
            required == WIREHAIR_V2_PROFILE_SERIALIZED_BYTES,
            "serialize size negotiation") ||
        !Check(wirehair_v2_profile_serialize(
            &profile, untouched,
            WIREHAIR_V2_PROFILE_SERIALIZED_BYTES - 1u, &required) ==
                WirehairV2_BufferTooSmall && untouched[0] == 0x5au,
            "serialize short buffer is transactional"))
    {
        return false;
    }

    profile.struct_bytes = 0u;
    if (!Check(wirehair_v2_profile_serialize(
            &profile, untouched, sizeof(untouched), &required) ==
            WirehairV2_InvalidSize,
            "host struct size"))
    {
        return false;
    }
    profile.struct_bytes = (uint32_t)sizeof(profile);
    profile.profile_version = 2u;
    if (!Check(wirehair_v2_profile_serialize(
            &profile, untouched, sizeof(untouched), &required) ==
            WirehairV2_UnsupportedVersion,
            "host struct version"))
    {
        return false;
    }
    profile.profile_version = WIREHAIR_V2_PROFILE_VERSION;
    profile.reserved[2] = 1u;
    if (!Check(wirehair_v2_profile_serialize(
            &profile, untouched, sizeof(untouched), &required) ==
            WirehairV2_ReservedNonzero,
            "host reserved field"))
    {
        return false;
    }

    WirehairV2Codec codec = reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
    uint8_t malformed[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES];
    std::memcpy(malformed, golden, sizeof(malformed));
    malformed[4] = 2u;
    if (!Check(wirehair_v2_decoder_create(
            malformed, sizeof(malformed), &codec) ==
                WirehairV2_UnsupportedVersion && codec == nullptr,
            "decoder rejects unknown version transactionally"))
    {
        return false;
    }
    return true;
}

bool CheckCppApi(
    const std::vector<uint8_t>& message,
    const uint8_t* golden_profile)
{
    wirehair::v2::SerializedProfile profile;
    std::memcpy(profile.data(), golden_profile, profile.size());
    WirehairV2Profile native = {};
    if (!Check(profile.Deserialize(native) == WirehairV2_Success,
            "C++ profile deserialize") ||
        !Check(profile.Validate() == WirehairV2_Success,
            "C++ profile validate"))
    {
        return false;
    }

    wirehair::v2::Encoder encoder;
    if (!Check(encoder.Create(message.data(), profile) ==
            WirehairV2_Success && encoder,
            "C++ encoder create from profile"))
    {
        return false;
    }
    wirehair::v2::Encoder moved(std::move(encoder));
    if (!Check(!encoder && moved, "C++ encoder move ownership")) {
        return false;
    }

    wirehair::v2::Decoder decoder;
    if (!Check(decoder.Create(profile) == WirehairV2_Success && decoder,
            "C++ decoder create"))
    {
        return false;
    }
    uint8_t block[16];
    WirehairV2Result result = WirehairV2_NeedMore;
    for (uint32_t id = 8u; id < 80u && result == WirehairV2_NeedMore; ++id)
    {
        uint32_t bytes = 0u;
        if (!Check(moved.Encode(id, block, sizeof(block), bytes) ==
                WirehairV2_Success && bytes == sizeof(block),
                "C++ repair encode"))
        {
            return false;
        }
        result = decoder.Decode(id, block, bytes);
    }
    std::vector<uint8_t> recovered(message.size(), 0u);
    uint64_t recovered_bytes = 0u;
    return Check(result == WirehairV2_Success,
            "C++ repair-only decode") &&
        Check(decoder.Recover(
            recovered.data(), recovered.size(), recovered_bytes) ==
                WirehairV2_Success,
            "C++ recover") &&
        Check(recovered_bytes == message.size() && recovered == message,
            "C++ recovered bytes");
}

bool CheckDescriptorOnlyDecoder(
    const std::vector<uint8_t>& message,
    const uint8_t* serialized_profile)
{
    const uint32_t block_bytes = 16u;
    const uint32_t block_count = 8u;
    WirehairV2Codec decoder = nullptr;
    if (!Check(wirehair_v2_decoder_create(
            serialized_profile, WIREHAIR_V2_PROFILE_SERIALIZED_BYTES,
            &decoder) == WirehairV2_Success,
            "golden descriptor-only decoder creation"))
    {
        return false;
    }

    WirehairV2Result result = WirehairV2_NeedMore;
    for (uint32_t id = 0u; id < block_count; ++id)
    {
        const uint32_t bytes = id + 1u == block_count ? 5u : block_bytes;
        result = wirehair_v2_decode(
            decoder,
            id,
            message.data() + (size_t)id * block_bytes,
            bytes);
        const WirehairV2Result expected = id + 1u == block_count ?
            WirehairV2_Success : WirehairV2_NeedMore;
        if (!Check(result == expected,
                "golden descriptor-only systematic feed"))
        {
            wirehair_v2_free(decoder);
            return false;
        }
    }

    std::vector<uint8_t> recovered(message.size(), 0u);
    uint64_t recovered_bytes = 0u;
    const bool ok = Check(wirehair_v2_recover(
            decoder, recovered.data(), recovered.size(), &recovered_bytes) ==
                WirehairV2_Success,
            "golden descriptor-only recover") &&
        Check(recovered_bytes == message.size() && recovered == message,
            "golden descriptor-only message bytes");
    wirehair_v2_free(decoder);
    return ok;
}

bool CheckNonzeroAttemptProfile()
{
    static const uint8_t ExpectedProfile[
        WIREHAIR_V2_PROFILE_SERIALIZED_BYTES] = {
        0x57, 0x48, 0x56, 0x32, 0x01, 0x00, 0x20, 0x00,
        0xc9, 0xf9, 0xf4, 0x47, 0xbb, 0x5b, 0x29, 0x4b,
        0x02, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x01, 0x00, 0x00, 0x00, 0x02, 0x00, 0x00, 0x00
    };
    const uint8_t message[2] = {0x31u, 0xa7u};
    uint8_t serialized[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES] = {};
    uint32_t serialized_bytes = 0u;
    WirehairV2Codec encoder = nullptr;
    WirehairV2Codec decoder = nullptr;
    if (!Check(wirehair_v2_encoder_create(
            message, sizeof(message), 1u,
            serialized, sizeof(serialized), &serialized_bytes, &encoder) ==
                WirehairV2_Success,
            "nonzero-attempt encoder") ||
        !Check(std::memcmp(
            serialized, ExpectedProfile, sizeof(serialized)) == 0,
            "nonzero-attempt profile golden") ||
        !Check(wirehair_v2_decoder_create(
            serialized, serialized_bytes, &decoder) == WirehairV2_Success,
            "nonzero-attempt decoder"))
    {
        wirehair_v2_free(encoder);
        wirehair_v2_free(decoder);
        return false;
    }

    WirehairV2Result result = WirehairV2_NeedMore;
    for (uint32_t id = 0u; id < 2u; ++id)
    {
        uint8_t block = 0u;
        uint32_t bytes = 0u;
        if (!Check(wirehair_v2_encode(
                encoder, id, &block, sizeof(block), &bytes) ==
                    WirehairV2_Success && bytes == 1u,
                "nonzero-attempt packet encode"))
        {
            wirehair_v2_free(encoder);
            wirehair_v2_free(decoder);
            return false;
        }
        result = wirehair_v2_decode(decoder, id, &block, bytes);
    }
    uint8_t repair = 0u;
    uint32_t repair_bytes = 0u;
    if (!Check(result == WirehairV2_Success,
            "nonzero-attempt systematic decode") ||
        !Check(wirehair_v2_encode(
                encoder, 2u, &repair, sizeof(repair), &repair_bytes) ==
                    WirehairV2_Success && repair_bytes == 1u,
            "nonzero-attempt repair encode"))
    {
        wirehair_v2_free(encoder);
        wirehair_v2_free(decoder);
        return false;
    }
    repair ^= 1u;
    if (!Check(wirehair_v2_decode(
                decoder, 2u, &repair, repair_bytes) == WirehairV2_Error,
            "nonzero-attempt corrupt repair rejection"))
    {
        wirehair_v2_free(encoder);
        wirehair_v2_free(decoder);
        return false;
    }
    repair ^= 1u;
    if (!Check(wirehair_v2_decode(
                decoder, 2u, &repair, repair_bytes) == WirehairV2_Success,
            "nonzero-attempt repair validation"))
    {
        wirehair_v2_free(encoder);
        wirehair_v2_free(decoder);
        return false;
    }
    uint8_t recovered[2] = {};
    uint64_t recovered_bytes = 0u;
    const bool ok = Check(wirehair_v2_recover(
            decoder, recovered, sizeof(recovered), &recovered_bytes) ==
                WirehairV2_Success,
            "nonzero-attempt recovery") &&
        Check(recovered_bytes == sizeof(message) &&
            std::memcmp(recovered, message, sizeof(message)) == 0,
            "nonzero-attempt recovered bytes");
    wirehair_v2_free(encoder);
    wirehair_v2_free(decoder);
    return ok;
}

bool CheckExplicitCertifiedProfile(
    const std::vector<uint8_t>& message,
    const uint8_t* golden,
    WirehairV2Codec default_encoder)
{
    uint8_t selected[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES] = {};
    uint32_t serialized_bytes = 0u;
    WirehairV2Codec explicit_encoder = nullptr;
    const bool created = Check(wirehair_v2_encoder_create_profile_id(
            WIREHAIR_V2_PROFILE_CERTIFIED_2026_07,
            message.data(), message.size(), 16u,
            selected, sizeof(selected), &serialized_bytes,
            &explicit_encoder) == WirehairV2_Success,
            "explicit certified-profile selection") &&
        Check(serialized_bytes == sizeof(selected) &&
            std::memcmp(selected, golden, sizeof(selected)) == 0,
            "explicit certified selector preserves descriptor golden");

    bool payload_ok = created;
    const uint32_t compare_ids[] = {0u, 12345u};
    for (uint32_t id : compare_ids)
    {
        uint8_t explicit_block[16] = {};
        uint8_t default_block[16] = {};
        uint32_t explicit_bytes = 0u;
        uint32_t default_bytes = 0u;
        payload_ok = payload_ok && Check(
            wirehair_v2_encode(
                explicit_encoder, id, explicit_block, sizeof(explicit_block),
                &explicit_bytes) == WirehairV2_Success &&
            wirehair_v2_encode(
                default_encoder, id, default_block, sizeof(default_block),
                &default_bytes) == WirehairV2_Success &&
            explicit_bytes == default_bytes &&
            std::memcmp(explicit_block, default_block, explicit_bytes) == 0,
            id == 0u ?
                "explicit certified systematic equation identity" :
                "explicit certified repair equation identity");
    }
    wirehair_v2_free(explicit_encoder);

    wirehair::v2::SerializedProfile cpp_profile;
    wirehair::v2::Encoder cpp_encoder;
    return payload_ok && Check(cpp_encoder.Create(
            WIREHAIR_V2_PROFILE_CERTIFIED_2026_07,
            message.data(), message.size(), 16u, cpp_profile) ==
                WirehairV2_Success,
            "C++ explicit certified-profile selection") &&
        Check(std::memcmp(
            cpp_profile.data(), golden, cpp_profile.size()) == 0,
            "C++ explicit certified selector preserves golden");
}

bool CheckRetiredProfileRejection(
    const std::vector<uint8_t>& message,
    const uint8_t* certified_golden)
{
    // These identifiers were once assigned to experimental equation systems.
    // Keep them as private tombstones: they must never be accepted or reused.
    static const uint64_t RetiredProfileIds[] = {
        UINT64_C(0xe161ce5d456f9bb7),
        UINT64_C(0x20a4f27a870612a2)
    };

    for (uint64_t retired_id : RetiredProfileIds)
    {
        WirehairV2Profile native = {};
        native.struct_bytes = (uint32_t)sizeof(native);
        native.profile_version = WIREHAIR_V2_PROFILE_VERSION;
        native.profile_id = retired_id;
        native.message_bytes = message.size();
        native.block_bytes = 16u;

        uint8_t output[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES];
        const auto reset_output = [&](uint8_t value) {
            std::memset(output, value, sizeof(output));
        };
        const auto output_is = [&](uint8_t value) {
            return std::all_of(
                output, output + sizeof(output),
                [=](uint8_t byte) { return byte == value; });
        };

        uint32_t output_bytes = 0u;
        reset_output(0x5au);
        if (!Check(wirehair_v2_profile_serialize(
                &native, output, sizeof(output), &output_bytes) ==
                    WirehairV2_UnsupportedProfile &&
                output_bytes == sizeof(output) && output_is(0x5au),
                "retired profile serialization rejected transactionally"))
        {
            return false;
        }

        uint8_t serialized[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES];
        std::memcpy(serialized, certified_golden, sizeof(serialized));
        for (unsigned byte = 0; byte < 8u; ++byte) {
            serialized[8u + byte] =
                (uint8_t)(retired_id >> (8u * byte));
        }

        if (!Check(wirehair_v2_profile_validate(
                serialized, sizeof(serialized)) ==
                    WirehairV2_UnsupportedProfile,
                "retired serialized profile validation rejected"))
        {
            return false;
        }

        WirehairV2Profile parsed;
        std::memset(&parsed, 0xa5, sizeof(parsed));
        if (!Check(wirehair_v2_profile_deserialize(
                serialized, sizeof(serialized), &parsed) ==
                    WirehairV2_UnsupportedProfile &&
                IsZeroProfile(parsed),
                "retired serialized profile parse rejected transactionally"))
        {
            return false;
        }

        WirehairV2Codec codec =
            reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
        if (!Check(wirehair_v2_encoder_create_profile(
                message.data(), serialized, sizeof(serialized), &codec) ==
                    WirehairV2_UnsupportedProfile &&
                codec == nullptr,
                "retired descriptor encoder rejected transactionally"))
        {
            return false;
        }

        codec = reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
        if (!Check(wirehair_v2_decoder_create(
                serialized, sizeof(serialized), &codec) ==
                    WirehairV2_UnsupportedProfile &&
                codec == nullptr,
                "retired descriptor decoder rejected transactionally"))
        {
            return false;
        }

        output_bytes = 0u;
        reset_output(0xa5u);
        codec = reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
        if (!Check(wirehair_v2_encoder_create_profile_id(
                retired_id,
                message.data(), message.size(), 16u,
                output, sizeof(output), &output_bytes, &codec) ==
                    WirehairV2_UnsupportedProfile &&
                output_bytes == sizeof(output) && output_is(0xa5u) &&
                codec == nullptr,
                "retired selector rejected transactionally"))
        {
            return false;
        }
    }

    return true;
}
} // namespace

int main()
{
    enum : uint32_t { BlockBytes = 16u, BlockCount = 8u };
    const uint64_t message_bytes = 117u;
    std::vector<uint8_t> message((size_t)message_bytes);
    FillMessage(message);

    static const uint8_t ExpectedProfile[
        WIREHAIR_V2_PROFILE_SERIALIZED_BYTES] = {
        0x57, 0x48, 0x56, 0x32, 0x01, 0x00, 0x20, 0x00,
        0xc9, 0xf9, 0xf4, 0x47, 0xbb, 0x5b, 0x29, 0x4b,
        0x75, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x10, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
    };

    uint8_t serialized[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES];
    std::memset(serialized, 0xcc, sizeof(serialized));
    uint32_t serialized_bytes = 0u;
    WirehairV2Codec encoder = nullptr;
    {
        WirehairV2Codec rejected_missing_size =
            reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
        if (!Check(wirehair_v2_encoder_create(
                message.data(), message.size(), BlockBytes,
                serialized, sizeof(serialized), nullptr,
                &rejected_missing_size) == WirehairV2_InvalidInput &&
                rejected_missing_size == nullptr,
                "encoder requires descriptor-size output"))
        {
            return 1;
        }
        uint8_t short_profile[
            WIREHAIR_V2_PROFILE_SERIALIZED_BYTES - 1u];
        std::memset(short_profile, 0x5a, sizeof(short_profile));
        uint32_t required = 0u;
        WirehairV2Codec rejected =
            reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
        if (!Check(wirehair_v2_encoder_create(
                message.data(), message.size(), BlockBytes,
                short_profile, sizeof(short_profile), &required, &rejected) ==
                    WirehairV2_BufferTooSmall && rejected == nullptr &&
                required == WIREHAIR_V2_PROFILE_SERIALIZED_BYTES &&
                short_profile[0] == 0x5au,
                "encoder profile size negotiation is transactional"))
        {
            return 1;
        }
    }
    if (!Check(wirehair_v2_encoder_create(
            message.data(), message.size(), BlockBytes,
            serialized, sizeof(serialized), &serialized_bytes, &encoder) ==
            WirehairV2_Success,
            "encoder selection") ||
        !Check(encoder != nullptr, "encoder handle") ||
        !Check(serialized_bytes == sizeof(serialized),
            "serialized profile size") ||
        !Check(std::memcmp(
            serialized, ExpectedProfile, sizeof(serialized)) == 0,
            "cross-endian golden serialized profile"))
    {
        wirehair_v2_free(encoder);
        return 1;
    }

    if (!CheckSelectingConstructorOverlapGuards(message) ||
        !CheckStagedDescriptorMessageAliases(message, ExpectedProfile) ||
        !CheckDescriptorConstructorOverlapGuards(message, serialized) ||
        !Check(wirehair_v2_encoder_create_profile(
            message.data(), serialized, sizeof(serialized), nullptr) ==
                WirehairV2_InvalidInput,
            "profile encoder still rejects a null codec output") ||
        !Check(wirehair_v2_decoder_create(
            serialized, sizeof(serialized), nullptr) ==
                WirehairV2_InvalidInput,
            "profile decoder still rejects a null codec output"))
    {
        wirehair_v2_free(encoder);
        return 1;
    }

    WirehairV2Profile in_place_profile;
    std::memcpy(&in_place_profile, serialized, sizeof(in_place_profile));
    if (!Check(wirehair_v2_profile_deserialize(
            &in_place_profile,
            sizeof(in_place_profile),
            &in_place_profile) == WirehairV2_Success &&
            in_place_profile.message_bytes == message.size() &&
            in_place_profile.block_bytes == BlockBytes,
            "in-place profile deserialize"))
    {
        wirehair_v2_free(encoder);
        return 1;
    }

    WirehairV2Profile profile = {};
    uint8_t reserialized[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES] = {};
    uint32_t reserialized_bytes = 0u;
    if (!Check(wirehair_v2_profile_deserialize(
            serialized, sizeof(serialized), &profile) == WirehairV2_Success,
            "deserialize selected profile") ||
        !Check(wirehair_v2_profile_validate(
            serialized, sizeof(serialized)) == WirehairV2_Success,
            "validate selected profile") ||
        !Check(profile.struct_bytes == sizeof(profile) &&
            profile.profile_version == WIREHAIR_V2_PROFILE_VERSION &&
            profile.profile_id == WIREHAIR_V2_PROFILE_CURRENT &&
            profile.message_bytes == message.size() &&
            profile.block_bytes == BlockBytes && profile.seed_attempt == 0u,
            "deserialized selected fields") ||
        !Check(wirehair_v2_profile_serialize(
            &profile, reserialized, sizeof(reserialized),
            &reserialized_bytes) == WirehairV2_Success &&
            reserialized_bytes == sizeof(reserialized) &&
            std::memcmp(reserialized, serialized, sizeof(serialized)) == 0,
            "profile canonical round trip"))
    {
        wirehair_v2_free(encoder);
        return 1;
    }

    const auto check_profile_output_alias = [
        &profile](uint32_t capacity, bool partial, const char* what)
    {
        std::vector<uint32_t> alias_words(
            (WIREHAIR_V2_PROFILE_SERIALIZED_BYTES + 2u +
                sizeof(uint32_t) - 1u) / sizeof(uint32_t),
            UINT32_C(0xa5a5a5a5));
        const std::vector<uint32_t> before = alias_words;
        uint8_t* const base =
            reinterpret_cast<uint8_t*>(alias_words.data());
        void* const output = base + (partial ? 2u : 0u);
        return Check(wirehair_v2_profile_serialize(
                &profile, output, capacity, alias_words.data()) ==
                    WirehairV2_InvalidInput &&
                alias_words == before,
            what);
    };
    if (!check_profile_output_alias(
            sizeof(reserialized), false,
            "exact profile size/output alias is rejected transactionally") ||
        !check_profile_output_alias(
            sizeof(reserialized) - 1u, false,
            "short exact profile size/output alias is transactional") ||
        !check_profile_output_alias(
            sizeof(reserialized), true,
            "partial profile size/output alias is rejected transactionally") ||
        !check_profile_output_alias(
            sizeof(reserialized) - 1u, true,
            "short partial profile size/output alias is transactional"))
    {
        wirehair_v2_free(encoder);
        return 1;
    }

    WirehairV2Profile profile_size_alias = profile;
    const WirehairV2Profile profile_size_alias_before = profile_size_alias;
    uint8_t profile_alias_output[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES];
    std::memset(profile_alias_output, 0xa5, sizeof(profile_alias_output));
    uint8_t profile_alias_output_before[sizeof(profile_alias_output)];
    std::memcpy(
        profile_alias_output_before,
        profile_alias_output,
        sizeof(profile_alias_output));
    if (!Check(wirehair_v2_profile_serialize(
            &profile_size_alias,
            profile_alias_output,
            sizeof(profile_alias_output),
            &profile_size_alias.profile_version) ==
                WirehairV2_InvalidInput &&
            std::memcmp(
                &profile_size_alias,
                &profile_size_alias_before,
                sizeof(profile_size_alias)) == 0 &&
            std::memcmp(
                profile_alias_output,
                profile_alias_output_before,
                sizeof(profile_alias_output)) == 0,
            "profile size/input alias is rejected transactionally"))
    {
        wirehair_v2_free(encoder);
        return 1;
    }

    WirehairV2Profile staged_profile_output = profile;
    uint32_t staged_profile_bytes = 0u;
    if (!Check(wirehair_v2_profile_serialize(
            &staged_profile_output,
            &staged_profile_output,
            sizeof(staged_profile_output),
            &staged_profile_bytes) == WirehairV2_Success &&
            staged_profile_bytes == sizeof(staged_profile_output) &&
            std::memcmp(
                &staged_profile_output,
                serialized,
                sizeof(staged_profile_output)) == 0,
            "staged profile input/output alias remains supported"))
    {
        wirehair_v2_free(encoder);
        return 1;
    }

    if (!CheckMalformedProfiles(ExpectedProfile)) {
        wirehair_v2_free(encoder);
        return 1;
    }
    if (!CheckDescriptorOnlyDecoder(message, ExpectedProfile)) {
        wirehair_v2_free(encoder);
        return 1;
    }

    WirehairV2Codec recreated = nullptr;
    WirehairV2Codec decoder = nullptr;
    if (!Check(wirehair_v2_encoder_create_profile(
            message.data(), serialized, sizeof(serialized), &recreated) ==
                WirehairV2_Success,
            "encoder recreation from descriptor") ||
        !Check(wirehair_v2_decoder_create(
            serialized, sizeof(serialized), &decoder) == WirehairV2_Success,
            "decoder creation from descriptor only"))
    {
        wirehair_v2_free(encoder);
        wirehair_v2_free(recreated);
        wirehair_v2_free(decoder);
        return 1;
    }
    const auto fail_after_create = [&]() {
        wirehair_v2_free(encoder);
        wirehair_v2_free(recreated);
        wirehair_v2_free(decoder);
        return 1;
    };

    const auto check_encode_output_alias = [
        encoder](uint32_t capacity, bool partial, const char* what)
    {
        std::vector<uint32_t> alias_words(
            (BlockBytes + 2u + sizeof(uint32_t) - 1u) /
                sizeof(uint32_t),
            UINT32_C(0xa5a5a5a5));
        const std::vector<uint32_t> before = alias_words;
        uint8_t* const base =
            reinterpret_cast<uint8_t*>(alias_words.data());
        void* const output = base + (partial ? 2u : 0u);
        return Check(wirehair_v2_encode(
                encoder, BlockCount, output, capacity,
                alias_words.data()) == WirehairV2_InvalidInput &&
                alias_words == before,
            what);
    };
    if (!check_encode_output_alias(
            BlockBytes, false,
            "exact encode size/output alias is rejected transactionally") ||
        !check_encode_output_alias(
            BlockBytes - 1u, false,
            "short exact encode size/output alias is transactional") ||
        !check_encode_output_alias(
            BlockBytes, true,
            "partial encode size/output alias is rejected transactionally") ||
        !check_encode_output_alias(
            BlockBytes - 1u, true,
            "short partial encode size/output alias is transactional"))
    {
        return fail_after_create();
    }

    uint8_t disjoint_output[BlockBytes];
    std::memset(disjoint_output, 0xa5, sizeof(disjoint_output));
    uint32_t disjoint_bytes = UINT32_C(0xa5a5a5a5);
    if (!Check(wirehair_v2_encode(
            encoder, BlockCount,
            disjoint_output, sizeof(disjoint_output), &disjoint_bytes) ==
                WirehairV2_Success && disjoint_bytes == BlockBytes,
            "disjoint encode size/output remains supported"))
    {
        return fail_after_create();
    }

    uint8_t short_output[BlockBytes];
    uint8_t short_before[BlockBytes];
    std::memset(short_output, 0xa5, sizeof(short_output));
    std::memcpy(short_before, short_output, sizeof(short_before));
    uint32_t required_packet_bytes = 0u;
    if (!Check(wirehair_v2_encode(
            encoder, BlockCount, short_output, BlockBytes - 1u,
            &required_packet_bytes) == WirehairV2_BufferTooSmall &&
            required_packet_bytes == BlockBytes &&
            std::memcmp(short_output, short_before, sizeof(short_output)) == 0,
            "short repair encode is transactional") ||
        !Check(wirehair_v2_encode(
            encoder, BlockCount - 1u, short_output, 4u,
            &required_packet_bytes) == WirehairV2_BufferTooSmall &&
            required_packet_bytes == 5u &&
            std::memcmp(short_output, short_before, sizeof(short_output)) == 0,
            "short final systematic encode is transactional"))
    {
        return fail_after_create();
    }
    uint8_t exact_final[5] = {};
    if (!Check(wirehair_v2_encode(
            encoder, BlockCount - 1u, exact_final, sizeof(exact_final),
            &required_packet_bytes) == WirehairV2_Success &&
            required_packet_bytes == sizeof(exact_final) &&
            std::memcmp(
                exact_final,
                message.data() + (BlockCount - 1u) * BlockBytes,
                sizeof(exact_final)) == 0,
            "exact-size final systematic encode"))
    {
        return fail_after_create();
    }

    const uint32_t packet_ids[] = {0u, 7u, 8u, 19u, 12345u};
    uint8_t block[BlockBytes];
    uint8_t recreated_block[BlockBytes];
    for (uint32_t id : packet_ids)
    {
        uint32_t bytes = 0u;
        uint32_t recreated_bytes = 0u;
        if (!Check(wirehair_v2_encode(
                encoder, id, block, sizeof(block), &bytes) ==
                    WirehairV2_Success,
                "golden packet encode") ||
            !Check(wirehair_v2_encode(
                recreated, id, recreated_block, sizeof(recreated_block),
                &recreated_bytes) == WirehairV2_Success,
                "recreated packet encode") ||
            !Check(bytes == recreated_bytes &&
                std::memcmp(block, recreated_block, bytes) == 0,
                "packet rows reproduced from descriptor"))
        {
            return fail_after_create();
        }
    }

    static const uint8_t ExpectedRepair12345[BlockBytes] = {
        0xb7, 0xe1, 0x38, 0xf6, 0x21, 0xc8, 0xf3, 0x68,
        0x21, 0x7f, 0x7e, 0xc4, 0xbd, 0x14, 0x0c, 0xc6
    };
    uint32_t golden_bytes = 0u;
    if (!Check(wirehair_v2_encode(
            encoder, 12345u, block, sizeof(block), &golden_bytes) ==
                WirehairV2_Success && golden_bytes == sizeof(block),
            "repair golden generation"))
    {
        return fail_after_create();
    }
    if (std::memcmp(block, ExpectedRepair12345, sizeof(block)) != 0)
    {
        std::fprintf(stderr, "V2 profile repair golden actual:");
        for (uint8_t byte : block) {
            std::fprintf(stderr, " 0x%02x", byte);
        }
        std::fprintf(stderr, "\n");
        return fail_after_create();
    }

    std::vector<uint8_t> recovered(message.size(), 0xa5u);
    const std::vector<uint8_t> before = recovered;
    uint64_t recovered_bytes = 0u;
    if (!Check(wirehair_v2_recover(
            decoder, recovered.data(), recovered.size(), &recovered_bytes) ==
                WirehairV2_NeedMore && recovered == before &&
            recovered_bytes == message.size(),
            "early recover is transactional") ||
        !Check(wirehair_v2_recover(
            decoder, recovered.data(), recovered.size() - 1u,
            &recovered_bytes) == WirehairV2_BufferTooSmall &&
            recovered == before,
            "short recover buffer"))
    {
        return fail_after_create();
    }

    WirehairV2Result decode_result = WirehairV2_NeedMore;
    for (uint32_t id = BlockCount;
        id < BlockCount + 72u && decode_result == WirehairV2_NeedMore;
        ++id)
    {
        uint32_t bytes = 0u;
        if (!Check(wirehair_v2_encode(
                encoder, id, block, sizeof(block), &bytes) ==
                    WirehairV2_Success && bytes == sizeof(block),
                "repair-only encode"))
        {
            return fail_after_create();
        }
        decode_result = wirehair_v2_decode(decoder, id, block, bytes);
    }
    if (!Check(decode_result == WirehairV2_Success,
            "repair-only descriptor decode") ||
        !Check(wirehair_v2_recover(
            decoder, recovered.data(), recovered.size(), &recovered_bytes) ==
                WirehairV2_Success,
            "descriptor message recover") ||
        !Check(recovered_bytes == message.size() && recovered == message,
            "descriptor recovered message"))
    {
        return fail_after_create();
    }

    std::vector<uint64_t> aliased_words(
        (message.size() + sizeof(uint64_t) - 1u) / sizeof(uint64_t),
        UINT64_C(0xa5a5a5a5a5a5a5a5));
    const std::vector<uint64_t> aliased_before = aliased_words;
    if (!Check(wirehair_v2_recover(
            decoder, aliased_words.data(), message.size(),
            aliased_words.data()) == WirehairV2_InvalidInput,
            "aliased recovery size output rejection") ||
        !Check(aliased_words == aliased_before,
            "aliased recovery size output was not no-write"))
    {
        return fail_after_create();
    }

    std::vector<uint8_t> unaligned(message.size() + 2u, 0x5au);
    recovered_bytes = 0u;
    if (!Check(wirehair_v2_recover(
            decoder, unaligned.data() + 1u, message.size(),
            &recovered_bytes) == WirehairV2_Success,
            "unaligned descriptor message recover") ||
        !Check(recovered_bytes == message.size() &&
            unaligned.front() == 0x5au && unaligned.back() == 0x5au &&
            std::equal(
                message.begin(), message.end(), unaligned.begin() + 1u),
            "unaligned recovery payload/bounds"))
    {
        return fail_after_create();
    }

    uint32_t ignored = 0u;
    if (!Check(wirehair_v2_encode(
            decoder, 0u, block, sizeof(block), &ignored) ==
                WirehairV2_InvalidInput,
            "decoder cannot encode") ||
        !Check(wirehair_v2_decode(
            encoder, 0u, block, sizeof(block)) == WirehairV2_InvalidInput,
            "encoder cannot decode") ||
        !Check(wirehair_v2_recover(
            encoder, recovered.data(), recovered.size(), &recovered_bytes) ==
                WirehairV2_InvalidInput,
            "encoder cannot recover") ||
        !Check(std::strcmp(wirehair_v2_result_string(
            WirehairV2_ReservedNonzero), "WirehairV2_ReservedNonzero") == 0,
            "stable result string") ||
        !Check(std::strcmp(wirehair_v2_result_string(
            (WirehairV2Result)999), "Unknown") == 0,
            "unknown result string"))
    {
        return fail_after_create();
    }

    if (!CheckCppApi(message, ExpectedProfile)) {
        return fail_after_create();
    }
    if (!CheckNonzeroAttemptProfile()) {
        return fail_after_create();
    }
    if (!CheckExplicitCertifiedProfile(message, ExpectedProfile, encoder)) {
        return fail_after_create();
    }
    if (!CheckRetiredProfileRejection(message, ExpectedProfile)) {
        return fail_after_create();
    }

    wirehair_v2_free(encoder);
    wirehair_v2_free(recreated);
    wirehair_v2_free(decoder);
    std::puts("V2 serialized profile golden/malformed/roundtrip: PASS");
    return 0;
}
