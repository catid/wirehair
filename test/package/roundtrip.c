#include "roundtrip.h"

#include <wirehair/wirehair.h>

#include <stdint.h>
#include <string.h>

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
    return 0;
}
