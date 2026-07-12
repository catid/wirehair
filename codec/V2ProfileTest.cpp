#include <wirehair/wirehair.h>
#include <wirehair/wirehair.hpp>

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
    uint8_t recovered[2] = {};
    uint64_t recovered_bytes = 0u;
    const bool ok = Check(result == WirehairV2_Success,
            "nonzero-attempt decode") &&
        Check(wirehair_v2_recover(
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

bool CheckMixedDescriptorContract(
    const std::vector<uint8_t>& message,
    const uint8_t* old_golden,
    WirehairV2Codec default_encoder)
{
    static const uint8_t MixedGolden[
        WIREHAIR_V2_PROFILE_SERIALIZED_BYTES] = {
        0x57, 0x48, 0x56, 0x32, 0x01, 0x00, 0x20, 0x00,
        0xb7, 0x9b, 0x6f, 0x45, 0x5d, 0xce, 0x61, 0xe1,
        0x75, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x10, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00
    };
    WirehairV2Profile mixed = {};
    mixed.struct_bytes = (uint32_t)sizeof(mixed);
    mixed.profile_version = WIREHAIR_V2_PROFILE_VERSION;
    mixed.profile_id = WIREHAIR_V2_PROFILE_MIXED_2026_07;
    mixed.message_bytes = message.size();
    mixed.block_bytes = 16u;
    uint8_t serialized[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES] = {};
    uint32_t serialized_bytes = 0u;
    WirehairV2Profile parsed = {};
    if (!Check(wirehair_v2_profile_serialize(
            &mixed, serialized, sizeof(serialized), &serialized_bytes) ==
                WirehairV2_Success &&
            serialized_bytes == sizeof(serialized) &&
            std::memcmp(serialized, MixedGolden, sizeof(serialized)) == 0,
            "mixed descriptor golden") ||
        !Check(wirehair_v2_profile_deserialize(
            serialized, sizeof(serialized), &parsed) ==
                WirehairV2_Success &&
            parsed.profile_id == WIREHAIR_V2_PROFILE_MIXED_2026_07 &&
            parsed.block_bytes == 16u,
            "mixed descriptor roundtrip"))
    {
        return false;
    }

    uint8_t untouched[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES];
    const auto reset_untouched = [&]() {
        std::memset(untouched, 0x5a, sizeof(untouched));
    };
    const auto all_untouched = [&]() {
        for (uint8_t value : untouched) {
            if (value != 0x5au) return false;
        }
        return true;
    };
    reset_untouched();
    mixed.block_bytes = 15u;
    if (!Check(wirehair_v2_profile_serialize(
            &mixed, untouched, sizeof(untouched), &serialized_bytes) ==
                WirehairV2_InvalidDimensions && all_untouched(),
            "mixed odd block rejected transactionally"))
    {
        return false;
    }

    WirehairV2Codec rejected =
        reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
    serialized_bytes = 0u;
    reset_untouched();
    if (!Check(wirehair_v2_encoder_create_profile_id(
            WIREHAIR_V2_PROFILE_MIXED_2026_07,
            message.data(), message.size(), 15u,
            untouched, sizeof(untouched), &serialized_bytes, &rejected) ==
                WirehairV2_InvalidDimensions &&
            rejected == nullptr && all_untouched(),
            "mixed selector rejects odd block before output write"))
    {
        return false;
    }
    rejected = reinterpret_cast<WirehairV2Codec>(uintptr_t(1));
    reset_untouched();
    if (!Check(wirehair_v2_encoder_create_profile_id(
            UINT64_C(0x0123456789abcdef),
            message.data(), message.size(), 16u,
            untouched, sizeof(untouched), &serialized_bytes, &rejected) ==
                WirehairV2_UnsupportedProfile &&
            rejected == nullptr && all_untouched(),
            "selector rejects unknown profile transactionally"))
    {
        return false;
    }

    uint8_t selected_mixed[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES] = {};
    WirehairV2Codec mixed_encoder = nullptr;
    WirehairV2Codec mixed_decoder = nullptr;
    serialized_bytes = 0u;
    if (!Check(wirehair_v2_encoder_create_profile_id(
            WIREHAIR_V2_PROFILE_MIXED_2026_07,
            message.data(), message.size(), 16u,
            selected_mixed, sizeof(selected_mixed), &serialized_bytes,
            &mixed_encoder) == WirehairV2_Success &&
            serialized_bytes == sizeof(selected_mixed) &&
            std::memcmp(
                selected_mixed, MixedGolden, sizeof(selected_mixed)) == 0,
            "mixed selector creates canonical encoder") ||
        !Check(wirehair_v2_decoder_create(
            selected_mixed, sizeof(selected_mixed), &mixed_decoder) ==
                WirehairV2_Success,
            "mixed descriptor creates decoder"))
    {
        wirehair_v2_free(mixed_encoder);
        wirehair_v2_free(mixed_decoder);
        return false;
    }
    const uint32_t mixed_block_count =
        (uint32_t)((message.size() + 15u) / 16u);
    WirehairV2Result decode_result = WirehairV2_NeedMore;
    uint8_t mixed_block[16] = {};
    for (uint32_t id = 0; id < mixed_block_count; ++id)
    {
        if (id == 3u) continue;
        uint32_t data_bytes = 0u;
        if (!Check(wirehair_v2_encode(
                mixed_encoder, id, mixed_block, sizeof(mixed_block),
                &data_bytes) == WirehairV2_Success,
                "mixed systematic encode"))
        {
            wirehair_v2_free(mixed_encoder);
            wirehair_v2_free(mixed_decoder);
            return false;
        }
        decode_result = wirehair_v2_decode(
            mixed_decoder, id, mixed_block, data_bytes);
    }
    for (uint32_t id = mixed_block_count;
         decode_result == WirehairV2_NeedMore &&
         id < mixed_block_count + 64u; ++id)
    {
        uint32_t data_bytes = 0u;
        if (!Check(wirehair_v2_encode(
                mixed_encoder, id, mixed_block, sizeof(mixed_block),
                &data_bytes) == WirehairV2_Success,
                "mixed repair encode"))
        {
            wirehair_v2_free(mixed_encoder);
            wirehair_v2_free(mixed_decoder);
            return false;
        }
        decode_result = wirehair_v2_decode(
            mixed_decoder, id, mixed_block, data_bytes);
    }
    std::vector<uint8_t> mixed_recovered(message.size(), 0u);
    uint64_t mixed_recovered_bytes = 0u;
    const bool mixed_e2e_ok = Check(
        decode_result == WirehairV2_Success &&
        wirehair_v2_recover(
            mixed_decoder, mixed_recovered.data(), mixed_recovered.size(),
            &mixed_recovered_bytes) == WirehairV2_Success &&
        mixed_recovered_bytes == message.size() &&
        mixed_recovered == message,
        "mixed public loss/repair recovery");
    wirehair_v2_free(mixed_encoder);
    wirehair_v2_free(mixed_decoder);
    if (!mixed_e2e_ok) return false;

    uint8_t selected_old[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES] = {};
    WirehairV2Codec old_encoder = nullptr;
    const bool old_ok = Check(wirehair_v2_encoder_create_profile_id(
            WIREHAIR_V2_PROFILE_CERTIFIED_2026_07,
            message.data(), message.size(), 16u,
            selected_old, sizeof(selected_old), &serialized_bytes,
            &old_encoder) == WirehairV2_Success,
            "explicit old-profile selection") &&
        Check(std::memcmp(
            selected_old, old_golden, sizeof(selected_old)) == 0,
            "explicit old selector preserves golden");
    bool old_payload_ok = old_ok;
    const uint32_t compare_ids[] = {0u, 12345u};
    for (uint32_t id : compare_ids)
    {
        uint8_t selected_block[16] = {};
        uint8_t default_block[16] = {};
        uint32_t selected_bytes = 0u;
        uint32_t default_bytes = 0u;
        old_payload_ok = old_payload_ok && Check(
            wirehair_v2_encode(
                old_encoder, id, selected_block, sizeof(selected_block),
                &selected_bytes) == WirehairV2_Success &&
            wirehair_v2_encode(
                default_encoder, id, default_block, sizeof(default_block),
                &default_bytes) == WirehairV2_Success &&
            selected_bytes == default_bytes &&
            std::memcmp(selected_block, default_block, selected_bytes) == 0,
            id == 0u ?
                "explicit old selector systematic equation identity" :
                "explicit old selector repair equation identity");
    }
    wirehair_v2_free(old_encoder);
    wirehair::v2::SerializedProfile cpp_profile;
    wirehair::v2::Encoder cpp_encoder;
    return old_payload_ok && Check(cpp_encoder.Create(
            WIREHAIR_V2_PROFILE_CERTIFIED_2026_07,
            message.data(), message.size(), 16u, cpp_profile) ==
                WirehairV2_Success,
            "C++ explicit old-profile selection") &&
        Check(std::memcmp(
            cpp_profile.data(), old_golden, cpp_profile.size()) == 0,
            "C++ explicit selector preserves golden");
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
    if (!CheckMixedDescriptorContract(message, ExpectedProfile, encoder)) {
        return fail_after_create();
    }

    wirehair_v2_free(encoder);
    wirehair_v2_free(recreated);
    wirehair_v2_free(decoder);
    std::puts("V2 serialized profile golden/malformed/roundtrip: PASS");
    return 0;
}
