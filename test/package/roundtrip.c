#include "roundtrip.h"

#include <wirehair/wirehair.h>

#include <stdint.h>
#include <string.h>

typedef char WirehairV2ProfileAbiSizeCheck[
    sizeof(WirehairV2Profile) == 32 ? 1 : -1];

static int wirehair_package_v2_round_trip(
    const uint8_t* message,
    uint32_t message_bytes,
    uint32_t block_bytes,
    uint32_t block_count,
    uint64_t profile_id)
{
    uint8_t profile[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES];
    uint8_t block[32];
    uint8_t recreated_block[32];
    uint8_t recovered[257];
    uint32_t profile_bytes = 0;
    uint64_t recovered_bytes = 0;
    WirehairV2Codec encoder = 0;
    WirehairV2Codec recreated = 0;
    WirehairV2Codec decoder = 0;
    WirehairV2Result decode_result = WirehairV2_NeedMore;
    uint32_t id;
    WirehairV2Profile parsed;

    if (block_bytes > sizeof(block) || message_bytes > sizeof(recovered) ||
        wirehair_v2_encoder_create_profile_id(
            profile_id, message, message_bytes, block_bytes,
            profile, sizeof(profile), &profile_bytes, &encoder) !=
                WirehairV2_Success ||
        profile_bytes != sizeof(profile) ||
        wirehair_v2_profile_deserialize(
            profile, profile_bytes, &parsed) != WirehairV2_Success ||
        parsed.profile_id != profile_id ||
        wirehair_v2_encoder_create_profile(
            message, profile, profile_bytes, &recreated) !=
                WirehairV2_Success ||
        wirehair_v2_decoder_create(
            profile, profile_bytes, &decoder) != WirehairV2_Success)
    {
        wirehair_v2_free(encoder);
        wirehair_v2_free(recreated);
        wirehair_v2_free(decoder);
        return 1;
    }

    for (id = block_count;
         id < block_count + 64u && decode_result == WirehairV2_NeedMore;
         ++id)
    {
        uint32_t written = 0;
        uint32_t recreated_written = 0;
        if (wirehair_v2_encode(
                encoder, id, block, sizeof(block), &written) !=
                    WirehairV2_Success ||
            wirehair_v2_encode(
                recreated, id, recreated_block, sizeof(recreated_block),
                &recreated_written) != WirehairV2_Success ||
            written != recreated_written ||
            memcmp(block, recreated_block, written) != 0)
        {
            wirehair_v2_free(encoder);
            wirehair_v2_free(recreated);
            wirehair_v2_free(decoder);
            return 2;
        }
        decode_result = wirehair_v2_decode(decoder, id, block, written);
    }
    if (decode_result != WirehairV2_Success ||
        wirehair_v2_recover(
            decoder, recovered, sizeof(recovered), &recovered_bytes) !=
                WirehairV2_Success ||
        recovered_bytes != message_bytes ||
        memcmp(recovered, message, message_bytes) != 0)
    {
        wirehair_v2_free(encoder);
        wirehair_v2_free(recreated);
        wirehair_v2_free(decoder);
        return 3;
    }

    wirehair_v2_free(encoder);
    wirehair_v2_free(recreated);
    wirehair_v2_free(decoder);
    return 0;
}

static int wirehair_package_v2_selector_failures(
    const uint8_t* message,
    uint32_t message_bytes)
{
    uint8_t profile[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES];
    uint32_t profile_bytes = 0;
    WirehairV2Codec codec;
    uint32_t i;

    memset(profile, 0xa5, sizeof(profile));
    codec = (WirehairV2Codec)(uintptr_t)1u;
    if (wirehair_v2_encoder_create_profile_id(
            WIREHAIR_V2_PROFILE_MIXED_2026_07,
            message, message_bytes, 31u,
            profile, sizeof(profile), &profile_bytes, &codec) !=
                WirehairV2_InvalidDimensions || codec != 0)
    {
        return 1;
    }
    for (i = 0; i < sizeof(profile); ++i) {
        if (profile[i] != 0xa5u) return 2;
    }

    memset(profile, 0x5a, sizeof(profile));
    codec = (WirehairV2Codec)(uintptr_t)1u;
    if (wirehair_v2_encoder_create_profile_id(
            UINT64_C(0x0123456789abcdef),
            message, message_bytes, 32u,
            profile, sizeof(profile), &profile_bytes, &codec) !=
                WirehairV2_UnsupportedProfile || codec != 0)
    {
        return 3;
    }
    for (i = 0; i < sizeof(profile); ++i) {
        if (profile[i] != 0x5au) return 4;
    }
    return 0;
}

int wirehair_package_round_trip(void)
{
    enum { MessageBytes = 257, BlockBytes = 32, BlockCount = 9 };
    static const unsigned reordered_ids[] = {
        5, 0, 12, 1, 7, 3, 20, 4, 8, 6, 9, 10, 11, 13, 14, 15
    };
    uint8_t message[MessageBytes];
    uint8_t recovered[MessageBytes];
    uint8_t block[BlockBytes];
    WirehairCodec encoder = 0;
    WirehairCodec decoder = 0;
    WirehairWireProfile profile;
    WirehairResult decode_result = Wirehair_NeedMore;
    uint32_t i;
    size_t order_i;

    if (wirehair_init() != Wirehair_Success) {
        return 1;
    }
    for (i = 0; i < (uint32_t)MessageBytes; ++i) {
        message[i] = (uint8_t)(i * 29u + 11u);
    }

    if (wirehair_wire_profile_init(
            WIREHAIR_LEGACY_PROFILE_CURRENT, &profile) != Wirehair_Success ||
        wirehair_encoder_create_profile_ex(
            0, message, MessageBytes, BlockBytes, &profile, 0, &encoder) !=
            Wirehair_Success ||
        wirehair_decoder_create_profile_ex(
            0, MessageBytes, BlockBytes, &profile, &decoder) !=
            Wirehair_Success)
    {
        wirehair_free(encoder);
        wirehair_free(decoder);
        return 2;
    }
    if (wirehair_encoder_detach_input(encoder) != Wirehair_Success ||
        wirehair_encoder_detach_input(encoder) != Wirehair_Success)
    {
        wirehair_free(encoder);
        wirehair_free(decoder);
        return 8;
    }

    for (order_i = 0;
         order_i < sizeof(reordered_ids) / sizeof(reordered_ids[0]);
         ++order_i)
    {
        uint32_t written = 0;
        const uint32_t block_id = reordered_ids[order_i];
        if (block_id == 2) {
            continue;
        }
        if (wirehair_encode(
                encoder, block_id, block, (uint32_t)BlockBytes, &written) !=
            Wirehair_Success)
        {
            wirehair_free(encoder);
            wirehair_free(decoder);
            return 3;
        }
        decode_result = wirehair_decode(decoder, block_id, block, written);
        if (decode_result == Wirehair_Success) {
            break;
        }
        if (decode_result != Wirehair_NeedMore) {
            wirehair_free(encoder);
            wirehair_free(decoder);
            return 4;
        }
    }
    if (decode_result != Wirehair_Success ||
        wirehair_recover(decoder, recovered, sizeof(recovered)) !=
            Wirehair_Success ||
        memcmp(recovered, message, sizeof(message)) != 0)
    {
        wirehair_free(encoder);
        wirehair_free(decoder);
        return 5;
    }

    for (i = 0; i < (uint32_t)BlockCount; ++i)
    {
        const uint32_t expected_bytes = i + 1u == (uint32_t)BlockCount ?
            1u : (uint32_t)BlockBytes;
        uint32_t written = 0;
        memset(block, 0, sizeof(block));
        if (wirehair_recover_block_ex(
                decoder, i, block, (uint32_t)BlockBytes, &written) !=
                Wirehair_Success ||
            written != expected_bytes ||
            memcmp(block,
                message + i * (uint32_t)BlockBytes,
                expected_bytes) != 0)
        {
            wirehair_free(encoder);
            wirehair_free(decoder);
            return 6;
        }
    }

    if (wirehair_decoder_becomes_encoder(decoder) != Wirehair_Success)
    {
        wirehair_free(encoder);
        wirehair_free(decoder);
        return 7;
    }
    {
        uint32_t written = 0;
        if (wirehair_encode(
                decoder, 2, block, (uint32_t)BlockBytes, &written) !=
                Wirehair_Success ||
            written != (uint32_t)BlockBytes ||
            memcmp(block, message + 2 * BlockBytes, BlockBytes) != 0)
        {
            wirehair_free(encoder);
            wirehair_free(decoder);
            return 8;
        }
    }

    wirehair_free(encoder);
    wirehair_free(decoder);
    if (wirehair_package_v2_round_trip(
            message, MessageBytes, BlockBytes, BlockCount,
            WIREHAIR_V2_PROFILE_CERTIFIED_2026_07) != 0)
    {
        return 9;
    }
    if (wirehair_package_v2_round_trip(
            message, MessageBytes, BlockBytes, BlockCount,
            WIREHAIR_V2_PROFILE_MIXED_2026_07) != 0)
    {
        return 10;
    }
    if (wirehair_package_v2_selector_failures(
            message, MessageBytes) != 0)
    {
        return 11;
    }
    return 0;
}
