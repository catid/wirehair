#include <wirehair/wirehair.h>

#include <stdint.h>
#include <string.h>

static int check_round_trip(int use_new_api)
{
    enum { MessageBytes = 257, BlockBytes = 32, BlockCount = 9 };
    uint8_t message[MessageBytes];
    uint8_t expected[MessageBytes];
    uint8_t block[BlockBytes];
    uint8_t recovered[MessageBytes];
    WirehairCodec encoder = 0;
    WirehairCodec decoder = 0;
    WirehairResult result = Wirehair_NeedMore;
    unsigned id;

    for (id = 0; id < MessageBytes; ++id) {
        message[id] = (uint8_t)(id * 19u + 7u);
    }
    memcpy(expected, message, sizeof(expected));

    if (use_new_api)
    {
        if (wirehair_encoder_create_owned_ex(
                0, message, MessageBytes, BlockBytes, &encoder) !=
            Wirehair_Success)
        {
            return 1;
        }
        if (wirehair_decoder_create_ex(
                0, MessageBytes, BlockBytes, &decoder) != Wirehair_Success)
        {
            wirehair_free(encoder);
            return 2;
        }
        memset(message, 0, sizeof(message));
    }
    else
    {
        encoder = wirehair_encoder_create(
            0, message, MessageBytes, BlockBytes);
        decoder = wirehair_decoder_create(0, MessageBytes, BlockBytes);
        if (!encoder || !decoder) {
            wirehair_free(encoder);
            wirehair_free(decoder);
            return 3;
        }
    }

    for (id = 0; id < BlockCount + 40u; ++id)
    {
        uint32_t written = 0;
        unsigned block_id = id;
        if (id == 2u) {
            continue;
        }
        if (wirehair_encode(
                encoder, block_id, block, sizeof(block), &written) !=
            Wirehair_Success)
        {
            wirehair_free(encoder);
            wirehair_free(decoder);
            return 4;
        }
        result = wirehair_decode(decoder, block_id, block, written);
        if (result == Wirehair_Success) {
            break;
        }
        if (result != Wirehair_NeedMore) {
            wirehair_free(encoder);
            wirehair_free(decoder);
            return 5;
        }
    }
    if (result != Wirehair_Success ||
        wirehair_recover(decoder, recovered, sizeof(recovered)) !=
            Wirehair_Success ||
        memcmp(recovered, expected, sizeof(expected)) != 0)
    {
        wirehair_free(encoder);
        wirehair_free(decoder);
        return 6;
    }

    for (id = 0; id < BlockCount; ++id)
    {
        const uint32_t expected_bytes = id + 1u == BlockCount ?
            1u : (uint32_t)BlockBytes;
        uint32_t written = 0;
        if (wirehair_recover_block_ex(
                decoder, id, block, expected_bytes, &written) !=
                Wirehair_Success ||
            written != expected_bytes ||
            memcmp(block, expected + id * BlockBytes, expected_bytes) != 0)
        {
            wirehair_free(encoder);
            wirehair_free(decoder);
            return 7;
        }
    }

    wirehair_free(encoder);
    wirehair_free(decoder);
    return 0;
}

int main(void)
{
    if (wirehair_init() != Wirehair_Success) {
        return 10;
    }
    if (check_round_trip(0) != 0) {
        return 11;
    }
    if (check_round_trip(1) != 0) {
        return 12;
    }
    return 0;
}
