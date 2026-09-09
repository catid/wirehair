#include <wirehair/wirehair.h>
#include <string.h>

/* Literal WHV2 K5 descriptor and repair packets, derived with carryless
 * multiplication modulo 0x14d from the sealed pair. No live encoder or explicit
 * initialization is required by the first receiver. */
int main(void)
{
    const unsigned char literal[32] = {
        'W','H','V','2',1,0,32,0, 0xf1,0x75,0xe3,0xbf,0x81,0x0c,0x07,0x80,
        9,0,0,0,0,0,0,0, 2,0,0,0,0,0,0,0
    };
    const unsigned char repairs[5][2] = {
        {193,39},{67,142},{27,167},{133,190},{80,164}
    };
    const unsigned char expected[9] = {0,1,2,3,4,5,6,7,8};
    unsigned char source[9] = {0,1,2,3,4,5,6,7,8};
    unsigned char profile[32], packet[2], output[9];
    WirehairV2Codec decoder = NULL, encoder = NULL;
    WirehairV2EncoderOptions options = WIREHAIR_V2_ENCODER_OPTIONS_INIT;
    unsigned id;
    uint32_t bytes = 0;
    uint64_t recovered = 0;
    int failed = 1;
    if (wirehair_v2_decoder_create(literal, sizeof(literal), &decoder) != WirehairV2_Success)
        goto cleanup;
    for (id = 5; id < 10; ++id) {
        if (wirehair_v2_decode(decoder, id, repairs[id - 5], 2) !=
            (id == 9 ? WirehairV2_Success : WirehairV2_NeedMore)) goto cleanup;
    }
    if (wirehair_v2_recover(decoder, output, sizeof(output), &recovered) != WirehairV2_Success ||
        recovered != sizeof(output) || memcmp(output, expected, sizeof(output))) goto cleanup;
    wirehair_v2_free(decoder);
    decoder = NULL;
    options.source_policy = WirehairV2EncoderSource_BorrowedImmutable;
    if (wirehair_v2_encoder_create_profile_id_with_options(WIREHAIR_V2_PROFILE_SMALL_K5_2026_09,
        source, sizeof(source), 2, &options, profile, sizeof(profile), &bytes, &encoder) !=
        WirehairV2_Success || bytes != sizeof(profile) || memcmp(profile, literal, sizeof(profile)))
        goto cleanup;
    if (wirehair_v2_encoder_detach_input(encoder) != WirehairV2_Success) goto cleanup;
    memset(source, 0xcc, sizeof(source));
    if (wirehair_v2_decoder_create(profile, sizeof(profile), &decoder) != WirehairV2_Success)
        goto cleanup;
    for (id = 5; id < 10; ++id) {
        if (wirehair_v2_encode(encoder, id, packet, sizeof(packet), &bytes) != WirehairV2_Success ||
            bytes != sizeof(packet) || memcmp(packet, repairs[id - 5], sizeof(packet)) ||
            wirehair_v2_decode(decoder, id, packet, bytes) !=
                (id == 9 ? WirehairV2_Success : WirehairV2_NeedMore)) goto cleanup;
    }
    if (wirehair_v2_recover(decoder, output, sizeof(output), &recovered) != WirehairV2_Success ||
        recovered != sizeof(output) || memcmp(output, expected, sizeof(output))) goto cleanup;
    failed = 0;
cleanup:
    wirehair_v2_free(encoder);
    wirehair_v2_free(decoder);
    return failed;
}
