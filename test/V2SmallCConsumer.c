#include <wirehair/wirehair.h>
#include <string.h>

/* Literal WHV2 descriptor; a receiver works before any encoder or explicit
 * wirehair_init call. This also checks implicit GF256 initialization. */
int main(void)
{
    const unsigned char literal[32] = {
        'W','H','V','2',1,0,32,0, 0x84,0xe1,0xa9,0xca,0x3e,0x04,0xc1,0x67,
        5,0,0,0,0,0,0,0, 2,0,0,0,0,0,0,0
    };
    const unsigned char expected[5] = {0,1,2,3,4};
    unsigned char source[5] = {0,1,2,3,4};
    unsigned char profile[32], packet[2], output[5];
    WirehairV2Codec decoder = NULL, encoder = NULL;
    WirehairV2EncoderOptions options = WIREHAIR_V2_ENCODER_OPTIONS_INIT;
    unsigned id;
    uint32_t bytes = 0;
    uint64_t recovered = 0;
    int failed = 1;
    if (wirehair_v2_decoder_create(literal, sizeof(literal), &decoder) != WirehairV2_Success)
        goto cleanup;
    for (id = 0; id < 3; ++id) {
        if (wirehair_v2_decode(decoder, id, expected + 2 * id, id == 2 ? 1 : 2) !=
            (id == 2 ? WirehairV2_Success : WirehairV2_NeedMore)) goto cleanup;
    }
    if (wirehair_v2_recover(decoder, output, sizeof(output), &recovered) != WirehairV2_Success ||
        recovered != sizeof(output) || memcmp(output, expected, sizeof(output))) goto cleanup;
    wirehair_v2_free(decoder);
    decoder = NULL;
    options.source_policy = WirehairV2EncoderSource_BorrowedImmutable;
    if (wirehair_v2_encoder_create_with_options(source, sizeof(source), 2, &options,
        profile, sizeof(profile), &bytes, &encoder) != WirehairV2_Success ||
        bytes != sizeof(profile) || memcmp(profile, literal, sizeof(profile))) goto cleanup;
    if (wirehair_v2_encoder_detach_input(encoder) != WirehairV2_Success) goto cleanup;
    memset(source, 0xcc, sizeof(source));
    if (wirehair_v2_decoder_create(profile, sizeof(profile), &decoder) != WirehairV2_Success)
        goto cleanup;
    for (id = 3; id < 6; ++id) {
        if (wirehair_v2_encode(encoder, id, packet, sizeof(packet), &bytes) != WirehairV2_Success ||
            bytes != sizeof(packet) || wirehair_v2_decode(decoder, id, packet, bytes) !=
                (id == 5 ? WirehairV2_Success : WirehairV2_NeedMore)) goto cleanup;
    }
    if (wirehair_v2_recover(decoder, output, sizeof(output), &recovered) != WirehairV2_Success ||
        recovered != sizeof(output) || memcmp(output, expected, sizeof(output))) goto cleanup;
    failed = 0;
cleanup:
    wirehair_v2_free(encoder);
    wirehair_v2_free(decoder);
    return failed;
}
