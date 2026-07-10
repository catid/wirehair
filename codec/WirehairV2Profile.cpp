#include <wirehair/wirehair.h>

#include "WirehairV2Codec.h"
#include "WirehairV2Plan.h"
#include "WirehairV2PrecodeEncode.h"
#include "WirehairV2Solve.h"

#include "../WirehairTools.h"

#include <cstring>
#include <limits>
#include <new>
#include <stdexcept>
#include <vector>

namespace {

static_assert(sizeof(WirehairV2Profile) == 32,
    "WirehairV2Profile is a fixed C ABI descriptor");
static_assert(WIREHAIR_V2_PROFILE_SERIALIZED_BYTES == 32u,
    "the canonical V2 profile encoding changed unexpectedly");
static_assert(wirehair_v2::kMaxPacketSeedAttempts == 256u,
    "the one-byte serialized seed attempt no longer covers the domain");

const uint8_t kProfileMagic[4] = { 'W', 'H', 'V', '2' };

enum class CodecMode
{
    Encoder,
    Decoder
};

struct PublicCodec
{
    wirehair_v2::Codec Impl;
    uint64_t MessageBytes = 0u;
    uint32_t BlockBytes = 0u;
    CodecMode Mode = CodecMode::Encoder;
    bool Decoded = false;
};

bool NamedV2ProfileAvailable()
{
    // WH_SEED_KNOBS can change the base seeds and packet degree distribution
    // at runtime.  Such experiment binaries must not publish or consume the
    // immutable certified profile ID.
#if defined(WH_SEED_KNOBS)
    return false;
#else
    return true;
#endif
}

uint16_t Load16LE(const uint8_t* data)
{
    return (uint16_t)((uint16_t)data[0] |
        ((uint16_t)data[1] << 8));
}

uint32_t Load32LE(const uint8_t* data)
{
    return (uint32_t)data[0] |
        ((uint32_t)data[1] << 8) |
        ((uint32_t)data[2] << 16) |
        ((uint32_t)data[3] << 24);
}

uint64_t Load64LE(const uint8_t* data)
{
    return (uint64_t)Load32LE(data) |
        ((uint64_t)Load32LE(data + 4) << 32);
}

void Store16LE(uint8_t* data, uint16_t value)
{
    data[0] = (uint8_t)value;
    data[1] = (uint8_t)(value >> 8);
}

void Store32LE(uint8_t* data, uint32_t value)
{
    data[0] = (uint8_t)value;
    data[1] = (uint8_t)(value >> 8);
    data[2] = (uint8_t)(value >> 16);
    data[3] = (uint8_t)(value >> 24);
}

void Store64LE(uint8_t* data, uint64_t value)
{
    Store32LE(data, (uint32_t)value);
    Store32LE(data + 4, (uint32_t)(value >> 32));
}

uint64_t BlockCountWide(uint64_t message_bytes, uint32_t block_bytes)
{
    if (block_bytes == 0u) {
        return 0u;
    }
    return message_bytes / block_bytes +
        (message_bytes % block_bytes != 0u ? 1u : 0u);
}

WirehairV2Result ValidateDimensions(
    uint64_t message_bytes,
    uint32_t block_bytes,
    uint32_t* block_count_out = nullptr)
{
    if (block_count_out) {
        *block_count_out = 0u;
    }
    if (message_bytes == 0u || block_bytes == 0u ||
        block_bytes > UINT32_C(0x7fffffff))
    {
        return WirehairV2_InvalidDimensions;
    }
    const uint64_t block_count = BlockCountWide(message_bytes, block_bytes);
    if (block_count < CAT_WIREHAIR_MIN_N ||
        block_count > CAT_WIREHAIR_MAX_N)
    {
        return WirehairV2_InvalidDimensions;
    }
    if (block_count_out) {
        *block_count_out = (uint32_t)block_count;
    }
    return WirehairV2_Success;
}

WirehairV2Result ValidateHostProfile(const WirehairV2Profile& profile)
{
    if (profile.struct_bytes != sizeof(WirehairV2Profile)) {
        return WirehairV2_InvalidSize;
    }
    if (profile.profile_version != WIREHAIR_V2_PROFILE_VERSION) {
        return WirehairV2_UnsupportedVersion;
    }
    if (profile.reserved[0] != 0u ||
        profile.reserved[1] != 0u ||
        profile.reserved[2] != 0u)
    {
        return WirehairV2_ReservedNonzero;
    }
    if (profile.profile_id != WIREHAIR_V2_PROFILE_CURRENT) {
        return WirehairV2_UnsupportedProfile;
    }
    if (!NamedV2ProfileAvailable()) {
        return WirehairV2_UnsupportedPlatform;
    }
    return ValidateDimensions(profile.message_bytes, profile.block_bytes);
}

WirehairV2Result MapResult(WirehairResult result)
{
    switch (result)
    {
    case Wirehair_Success:             return WirehairV2_Success;
    case Wirehair_NeedMore:            return WirehairV2_NeedMore;
    case Wirehair_InvalidInput:        return WirehairV2_InvalidInput;
    case Wirehair_BadDenseSeed:
    case Wirehair_BadPeelSeed:         return WirehairV2_BadSeed;
    case Wirehair_BadInput_SmallN:
    case Wirehair_BadInput_LargeN:     return WirehairV2_InvalidDimensions;
    case Wirehair_ExtraInsufficient:   return WirehairV2_ExtraInsufficient;
    case Wirehair_OOM:                 return WirehairV2_OOM;
    case Wirehair_UnsupportedPlatform: return WirehairV2_UnsupportedPlatform;
    case Wirehair_Error:
    default:                           return WirehairV2_Error;
    }
}

void ClearProfile(WirehairV2Profile* profile)
{
    if (profile) {
        std::memset(profile, 0, sizeof(*profile));
    }
}

PublicCodec* FromHandle(WirehairV2Codec codec)
{
    return reinterpret_cast<PublicCodec*>(codec);
}

WirehairV2Codec ToHandle(PublicCodec* codec)
{
    return reinterpret_cast<WirehairV2Codec>(codec);
}

wirehair_v2::SeedProfile ExpandProfile(const WirehairV2Profile& profile)
{
    const uint32_t block_count =
        (uint32_t)BlockCountWide(profile.message_bytes, profile.block_bytes);
    wirehair_v2::SeedProfile expanded =
        wirehair_v2::SelectSeedProfile(block_count, profile.block_bytes);
    const wirehair_v2::MessagePrecodeEncoderOptions options;

    wirehair_v2::PrecodeParams params = wirehair_v2::MakeCertifiedParams(
        block_count,
        wirehair_v2::MatrixSeedFromProfile(
            expanded, 0u, options.PrecodeSeedSalt));
    params.Staircase = expanded.DenseCount;
    params.DenseIdentityCorner = options.DenseIdentityCorner;
    params = wirehair_v2::PrecodeParamsForAttempt(
        params, profile.seed_attempt);

    wirehair_v2::PacketRowConfig packet_config;
    packet_config.PeelSeed = wirehair_v2::PacketPeelSeedFromProfile(
        expanded, options.RecoveryRowSeedSalt);
    packet_config.MixCount = options.RecoveryMixCount;
    packet_config = wirehair_v2::PacketConfigForAttempt(
        packet_config, profile.seed_attempt);

    wirehair_v2::PrecodeSystem system;
    system.Params = params;
    wirehair_v2::BindMessagePrecodeProfile(
        expanded, options, system, packet_config, profile.seed_attempt);
    return expanded;
}

WirehairV2Result MakePublicProfile(
    const wirehair_v2::SeedProfile& expanded,
    uint64_t message_bytes,
    WirehairV2Profile& profile)
{
    if (!expanded.V2SeedSelected ||
        expanded.V2SeedAttempt >= wirehair_v2::kMaxPacketSeedAttempts ||
        expanded.V2PrecodeContractVersion !=
            wirehair_v2::kPrecodeContractVersion ||
        expanded.V2PacketRowContractVersion !=
            wirehair_v2::kPacketRowContractVersion)
    {
        return WirehairV2_Error;
    }
    profile = WirehairV2Profile{};
    profile.struct_bytes = (uint32_t)sizeof(WirehairV2Profile);
    profile.profile_version = WIREHAIR_V2_PROFILE_VERSION;
    profile.profile_id = WIREHAIR_V2_PROFILE_CURRENT;
    profile.message_bytes = message_bytes;
    profile.block_bytes = expanded.BlockBytes;
    profile.seed_attempt = (uint8_t)expanded.V2SeedAttempt;

    const WirehairV2Result validation = ValidateHostProfile(profile);
    if (validation != WirehairV2_Success) {
        return validation;
    }

    // The compact descriptor deliberately omits every derived value.  Rebuild
    // them and compare here so a future implementation cannot publish the
    // current profile ID after changing one of its frozen equation constants.
    const wirehair_v2::SeedProfile canonical = ExpandProfile(profile);
    if (canonical.BlockCount != expanded.BlockCount ||
        canonical.BlockBytes != expanded.BlockBytes ||
        canonical.DenseCount != expanded.DenseCount ||
        canonical.PeelSeed != expanded.PeelSeed ||
        canonical.DenseSeed != expanded.DenseSeed ||
        canonical.V2SeedAttempt != expanded.V2SeedAttempt ||
        canonical.V2PrecodeContractVersion !=
            expanded.V2PrecodeContractVersion ||
        canonical.V2PacketRowContractVersion !=
            expanded.V2PacketRowContractVersion ||
        canonical.V2StaircaseCount != expanded.V2StaircaseCount ||
        canonical.V2DenseRowCount != expanded.V2DenseRowCount ||
        canonical.V2HeavyRowCount != expanded.V2HeavyRowCount ||
        canonical.V2SourceHits != expanded.V2SourceHits ||
        canonical.V2PrecodeSeed != expanded.V2PrecodeSeed ||
        canonical.V2PacketPeelSeed != expanded.V2PacketPeelSeed ||
        canonical.V2RecoveryMixCount != expanded.V2RecoveryMixCount ||
        canonical.V2DenseIdentityCorner !=
            expanded.V2DenseIdentityCorner ||
        canonical.V2PrecodeSeedSalt != expanded.V2PrecodeSeedSalt ||
        canonical.V2RecoveryRowSeedSalt !=
            expanded.V2RecoveryRowSeedSalt)
    {
        return WirehairV2_Error;
    }
    return WirehairV2_Success;
}

WirehairV2Result CreateEncoderForProfile(
    const void* message,
    const WirehairV2Profile& profile,
    PublicCodec*& codec_out)
{
    codec_out = nullptr;
    if (!message) {
        return WirehairV2_InvalidInput;
    }
    if (profile.message_bytes >
        (uint64_t)std::numeric_limits<size_t>::max())
    {
        return WirehairV2_UnsupportedPlatform;
    }

    PublicCodec* codec = new (std::nothrow) PublicCodec;
    if (!codec) {
        return WirehairV2_OOM;
    }
    const wirehair_v2::SeedProfile expanded = ExpandProfile(profile);
    const WirehairResult result = codec->Impl.InitializePrecodeEncoder(
        message,
        profile.message_bytes,
        profile.block_bytes,
        &expanded);
    if (result != Wirehair_Success) {
        delete codec;
        return MapResult(result);
    }
    codec->MessageBytes = profile.message_bytes;
    codec->BlockBytes = profile.block_bytes;
    codec->Mode = CodecMode::Encoder;
    codec->Decoded = false;
    codec_out = codec;
    return WirehairV2_Success;
}

} // namespace

WIREHAIR_EXPORT const char* wirehair_v2_result_string(
    WirehairV2Result result)
{
    switch (result)
    {
    case WirehairV2_Success:             return "WirehairV2_Success";
    case WirehairV2_NeedMore:            return "WirehairV2_NeedMore";
    case WirehairV2_InvalidInput:        return "WirehairV2_InvalidInput";
    case WirehairV2_BufferTooSmall:      return "WirehairV2_BufferTooSmall";
    case WirehairV2_InvalidMagic:        return "WirehairV2_InvalidMagic";
    case WirehairV2_UnsupportedVersion:  return "WirehairV2_UnsupportedVersion";
    case WirehairV2_InvalidSize:         return "WirehairV2_InvalidSize";
    case WirehairV2_ReservedNonzero:     return "WirehairV2_ReservedNonzero";
    case WirehairV2_UnsupportedProfile:  return "WirehairV2_UnsupportedProfile";
    case WirehairV2_InvalidDimensions:   return "WirehairV2_InvalidDimensions";
    case WirehairV2_BadSeed:             return "WirehairV2_BadSeed";
    case WirehairV2_ExtraInsufficient:   return "WirehairV2_ExtraInsufficient";
    case WirehairV2_Error:               return "WirehairV2_Error";
    case WirehairV2_OOM:                 return "WirehairV2_OOM";
    case WirehairV2_UnsupportedPlatform: return "WirehairV2_UnsupportedPlatform";
    default:                             return "Unknown";
    }
}

WIREHAIR_EXPORT WirehairV2Result wirehair_v2_profile_serialize(
    const WirehairV2Profile* profile,
    void* output,
    uint32_t outputCapacity,
    uint32_t* bytesOut)
{
    if (bytesOut) {
        *bytesOut = WIREHAIR_V2_PROFILE_SERIALIZED_BYTES;
    }
    if (!profile) {
        return WirehairV2_InvalidInput;
    }
    const WirehairV2Result validation = ValidateHostProfile(*profile);
    if (validation != WirehairV2_Success) {
        return validation;
    }
    if (outputCapacity < WIREHAIR_V2_PROFILE_SERIALIZED_BYTES || !output) {
        return WirehairV2_BufferTooSmall;
    }

    uint8_t encoded[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES] = {};
    std::memcpy(encoded, kProfileMagic, sizeof(kProfileMagic));
    Store16LE(encoded + 4,
        (uint16_t)WIREHAIR_V2_PROFILE_ENCODING_VERSION);
    Store16LE(encoded + 6,
        (uint16_t)WIREHAIR_V2_PROFILE_SERIALIZED_BYTES);
    Store64LE(encoded + 8, profile->profile_id);
    Store64LE(encoded + 16, profile->message_bytes);
    Store32LE(encoded + 24, profile->block_bytes);
    encoded[28] = profile->seed_attempt;
    std::memcpy(output, encoded, sizeof(encoded));
    return WirehairV2_Success;
}

WIREHAIR_EXPORT WirehairV2Result wirehair_v2_profile_deserialize(
    const void* serializedProfile,
    uint32_t serializedBytes,
    WirehairV2Profile* profileOut)
{
    if (!profileOut) {
        return WirehairV2_InvalidInput;
    }
    if (!serializedProfile) {
        ClearProfile(profileOut);
        return WirehairV2_InvalidInput;
    }
    // Snapshot before clearing profileOut so exact or partial input/output
    // aliasing cannot destroy the serialized record while it is parsed.
    uint8_t snapshot[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES] = {};
    const uint32_t snapshot_bytes =
        serializedBytes < WIREHAIR_V2_PROFILE_SERIALIZED_BYTES ?
            serializedBytes : WIREHAIR_V2_PROFILE_SERIALIZED_BYTES;
    if (snapshot_bytes != 0u) {
        std::memcpy(snapshot, serializedProfile, snapshot_bytes);
    }
    ClearProfile(profileOut);
    const uint8_t* encoded = snapshot;
    if (serializedBytes < sizeof(kProfileMagic)) {
        return WirehairV2_InvalidSize;
    }
    if (std::memcmp(encoded, kProfileMagic, sizeof(kProfileMagic)) != 0) {
        return WirehairV2_InvalidMagic;
    }
    if (serializedBytes < 8u) {
        return WirehairV2_InvalidSize;
    }
    const uint16_t version = Load16LE(encoded + 4);
    if (version != WIREHAIR_V2_PROFILE_ENCODING_VERSION) {
        return WirehairV2_UnsupportedVersion;
    }
    const uint16_t declared_bytes = Load16LE(encoded + 6);
    if (declared_bytes != WIREHAIR_V2_PROFILE_SERIALIZED_BYTES ||
        serializedBytes != declared_bytes)
    {
        return WirehairV2_InvalidSize;
    }

    WirehairV2Profile parsed = {};
    parsed.struct_bytes = (uint32_t)sizeof(WirehairV2Profile);
    parsed.profile_version = WIREHAIR_V2_PROFILE_VERSION;
    parsed.profile_id = Load64LE(encoded + 8);
    parsed.message_bytes = Load64LE(encoded + 16);
    parsed.block_bytes = Load32LE(encoded + 24);
    parsed.seed_attempt = encoded[28];
    parsed.reserved[0] = encoded[29];
    parsed.reserved[1] = encoded[30];
    parsed.reserved[2] = encoded[31];
    const WirehairV2Result validation = ValidateHostProfile(parsed);
    if (validation != WirehairV2_Success) {
        return validation;
    }
    *profileOut = parsed;
    return WirehairV2_Success;
}

WIREHAIR_EXPORT WirehairV2Result wirehair_v2_profile_validate(
    const void* serializedProfile,
    uint32_t serializedBytes)
{
    WirehairV2Profile ignored = {};
    return wirehair_v2_profile_deserialize(
        serializedProfile, serializedBytes, &ignored);
}

WIREHAIR_EXPORT WirehairV2Result wirehair_v2_encoder_create(
    const void* message,
    uint64_t messageBytes,
    uint32_t blockBytes,
    void* serializedProfileOut,
    uint32_t serializedProfileCapacity,
    uint32_t* serializedProfileBytesOut,
    WirehairV2Codec* codecOut)
{
    if (serializedProfileBytesOut) {
        *serializedProfileBytesOut = WIREHAIR_V2_PROFILE_SERIALIZED_BYTES;
    }
    if (!codecOut) {
        return WirehairV2_InvalidInput;
    }
    *codecOut = nullptr;
    if (!serializedProfileBytesOut) {
        return WirehairV2_InvalidInput;
    }
    if (serializedProfileCapacity < WIREHAIR_V2_PROFILE_SERIALIZED_BYTES ||
        !serializedProfileOut)
    {
        return WirehairV2_BufferTooSmall;
    }
    if (!NamedV2ProfileAvailable()) {
        return WirehairV2_UnsupportedPlatform;
    }
    const WirehairV2Result dimensions =
        ValidateDimensions(messageBytes, blockBytes);
    if (dimensions != WirehairV2_Success) {
        return dimensions;
    }
    if (!message) {
        return WirehairV2_InvalidInput;
    }
    if (messageBytes > (uint64_t)std::numeric_limits<size_t>::max()) {
        return WirehairV2_UnsupportedPlatform;
    }

    PublicCodec* codec = new (std::nothrow) PublicCodec;
    if (!codec) {
        return WirehairV2_OOM;
    }
    const WirehairResult create_result = codec->Impl.InitializePrecodeEncoder(
        message, messageBytes, blockBytes);
    if (create_result != Wirehair_Success) {
        delete codec;
        return MapResult(create_result);
    }

    WirehairV2Profile profile = {};
    WirehairV2Result result = MakePublicProfile(
        codec->Impl.Profile(), messageBytes, profile);
    uint8_t encoded[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES];
    uint32_t encoded_bytes = 0u;
    if (result == WirehairV2_Success) {
        result = wirehair_v2_profile_serialize(
            &profile, encoded, sizeof(encoded), &encoded_bytes);
    }
    if (result != WirehairV2_Success ||
        encoded_bytes != WIREHAIR_V2_PROFILE_SERIALIZED_BYTES)
    {
        delete codec;
        return result == WirehairV2_Success ? WirehairV2_Error : result;
    }

    codec->MessageBytes = messageBytes;
    codec->BlockBytes = blockBytes;
    codec->Mode = CodecMode::Encoder;
    codec->Decoded = false;
    std::memcpy(serializedProfileOut, encoded, sizeof(encoded));
    *codecOut = ToHandle(codec);
    return WirehairV2_Success;
}

WIREHAIR_EXPORT WirehairV2Result wirehair_v2_encoder_create_profile(
    const void* message,
    const void* serializedProfile,
    uint32_t serializedProfileBytes,
    WirehairV2Codec* codecOut)
{
    if (!codecOut) {
        return WirehairV2_InvalidInput;
    }
    *codecOut = nullptr;
    WirehairV2Profile profile = {};
    const WirehairV2Result parse_result = wirehair_v2_profile_deserialize(
        serializedProfile, serializedProfileBytes, &profile);
    if (parse_result != WirehairV2_Success) {
        return parse_result;
    }
    PublicCodec* codec = nullptr;
    const WirehairV2Result create_result =
        CreateEncoderForProfile(message, profile, codec);
    if (create_result == WirehairV2_Success) {
        *codecOut = ToHandle(codec);
    }
    return create_result;
}

WIREHAIR_EXPORT WirehairV2Result wirehair_v2_decoder_create(
    const void* serializedProfile,
    uint32_t serializedProfileBytes,
    WirehairV2Codec* codecOut)
{
    if (!codecOut) {
        return WirehairV2_InvalidInput;
    }
    *codecOut = nullptr;
    WirehairV2Profile profile = {};
    const WirehairV2Result parse_result = wirehair_v2_profile_deserialize(
        serializedProfile, serializedProfileBytes, &profile);
    if (parse_result != WirehairV2_Success) {
        return parse_result;
    }
    if (profile.message_bytes >
        (uint64_t)std::numeric_limits<size_t>::max())
    {
        return WirehairV2_UnsupportedPlatform;
    }

    PublicCodec* codec = new (std::nothrow) PublicCodec;
    if (!codec) {
        return WirehairV2_OOM;
    }
    const wirehair_v2::SeedProfile expanded = ExpandProfile(profile);
    const WirehairResult result = codec->Impl.InitializePrecodeDecoder(
        profile.message_bytes, profile.block_bytes, &expanded);
    if (result != Wirehair_Success) {
        delete codec;
        return MapResult(result);
    }
    codec->MessageBytes = profile.message_bytes;
    codec->BlockBytes = profile.block_bytes;
    codec->Mode = CodecMode::Decoder;
    codec->Decoded = false;
    *codecOut = ToHandle(codec);
    return WirehairV2_Success;
}

WIREHAIR_EXPORT WirehairV2Result wirehair_v2_encode(
    WirehairV2Codec codec,
    uint32_t blockId,
    void* blockDataOut,
    uint32_t outputCapacity,
    uint32_t* dataBytesOut)
{
    if (dataBytesOut) {
        *dataBytesOut = 0u;
    }
    PublicCodec* impl = FromHandle(codec);
    if (!impl || impl->Mode != CodecMode::Encoder ||
        !blockDataOut || !dataBytesOut)
    {
        return WirehairV2_InvalidInput;
    }
    const uint64_t block_count = BlockCountWide(
        impl->MessageBytes, impl->BlockBytes);
    uint32_t required = impl->BlockBytes;
    if (block_count != 0u && blockId == block_count - 1u) {
        required = (uint32_t)(impl->MessageBytes -
            (block_count - 1u) * impl->BlockBytes);
    }
    *dataBytesOut = required;
    if (outputCapacity < required) {
        return WirehairV2_BufferTooSmall;
    }
    return MapResult(impl->Impl.Encode(
        blockId, blockDataOut, outputCapacity, dataBytesOut));
}

WIREHAIR_EXPORT WirehairV2Result wirehair_v2_decode(
    WirehairV2Codec codec,
    uint32_t blockId,
    const void* blockData,
    uint32_t dataBytes)
{
    PublicCodec* impl = FromHandle(codec);
    if (!impl || impl->Mode != CodecMode::Decoder || !blockData) {
        return WirehairV2_InvalidInput;
    }
    const WirehairV2Result result = MapResult(
        impl->Impl.Decode(blockId, blockData, dataBytes));
    if (result == WirehairV2_Success) {
        impl->Decoded = true;
    }
    return result;
}

WIREHAIR_EXPORT WirehairV2Result wirehair_v2_recover(
    WirehairV2Codec codec,
    void* messageOut,
    uint64_t outputCapacity,
    uint64_t* bytesOut)
{
    PublicCodec* impl = FromHandle(codec);
    if (!impl || impl->Mode != CodecMode::Decoder) {
        if (bytesOut) {
            *bytesOut = 0u;
        }
        return WirehairV2_InvalidInput;
    }
    if (bytesOut) {
        *bytesOut = impl->MessageBytes;
    }
    if (outputCapacity < impl->MessageBytes) {
        return WirehairV2_BufferTooSmall;
    }
    if (!messageOut) {
        return WirehairV2_InvalidInput;
    }
    if (!impl->Decoded) {
        return WirehairV2_NeedMore;
    }
    if (impl->MessageBytes >
        (uint64_t)std::numeric_limits<size_t>::max())
    {
        return WirehairV2_UnsupportedPlatform;
    }

    try
    {
        // Recover transactionally: the internal decoder may allocate while
        // evaluating later rows, so stage the full message before publishing.
        std::vector<uint8_t> recovered((size_t)impl->MessageBytes);
        const WirehairV2Result result = MapResult(impl->Impl.Recover(
            recovered.data(), impl->MessageBytes));
        if (result != WirehairV2_Success) {
            return result;
        }
        std::memcpy(messageOut, recovered.data(), recovered.size());
        return WirehairV2_Success;
    }
    catch (const std::bad_alloc&) {
        return WirehairV2_OOM;
    }
    catch (const std::length_error&) {
        return WirehairV2_OOM;
    }
}

WIREHAIR_EXPORT void wirehair_v2_free(WirehairV2Codec codec)
{
    delete FromHandle(codec);
}
