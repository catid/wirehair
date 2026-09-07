#include "wirehair/wirehair_k6.h"
#include "wirehair/wirehair.h"
#include <string.h>

/* Runs before any encoder has existed in this process. Literal descriptor and
 * systematic packets suffice; no fixture builder can hide private state. */
static int standalone(void)
{
    const unsigned char profile[32] = {
        'W','H','K','6',1,0,32,0, 0x31,0x4d,0x54,0x36,0x4b,0x32,0x48,0x57,
        11,0,0,0,0,0,0,0, 2,0,0,0,0,0,0,0
    };
    const unsigned char source[11] = {0,1,2,3,4,5,6,7,8,9,10};
    unsigned char recovered[11];
    WirehairK6CreateResult d = wirehair_k6_decoder_create(profile, sizeof(profile));
    unsigned id;
    if (d.status != WirehairK6_Success) return 1;
    for (id = 0; id < 6; ++id) {
        WirehairK6Status s = wirehair_k6_decode(d.codec, id, source + 2 * id, id == 5 ? 1 : 2);
        if (s != (id == 5 ? WirehairK6_Success : WirehairK6_NeedMore)) {
            wirehair_k6_free(d.codec); return 1;
        }
    }
    {
        WirehairK6Result r = wirehair_k6_recover(d.codec, recovered, sizeof(recovered));
        int failed = r.status != WirehairK6_Success || r.bytes_written != sizeof(source) ||
            memcmp(source, recovered, sizeof(source));
        wirehair_k6_free(d.codec);
        return failed;
    }
}

/* A separate C translation unit exercises the actual C ABI, not C++ casts. */
int main(void)
{
    unsigned char source[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    const unsigned char expected[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    unsigned char profile[WIREHAIR_K6_PROFILE_BYTES], packet[2], recovered[11];
    WirehairK6CreateResult encoder, decoder;
    unsigned id;
    if (wirehair_init() != Wirehair_Success) return 1;
    if (standalone()) return 6;
    encoder = wirehair_k6_encoder_create(source, sizeof(source), 2, WirehairK6_BorrowedImmutable,
                                    profile, sizeof(profile));
    if (encoder.status != WirehairK6_Success) return 2;
    if (wirehair_k6_profile_validate(profile, sizeof(profile)) != WirehairK6_Success) {
        wirehair_k6_free(encoder.codec); return 7;
    }
    wirehair_k6_free(encoder.codec);
    encoder = wirehair_k6_encoder_create_profile(source, profile, sizeof(profile),
                                                  WirehairK6_BorrowedImmutable);
    if (encoder.status != WirehairK6_Success) return 8;
    if (wirehair_k6_encoder_detach_input(encoder.codec) != WirehairK6_Success) {
        wirehair_k6_free(encoder.codec); return 9;
    }
    memset(source, 0, sizeof(source));
    decoder = wirehair_k6_decoder_create(profile, sizeof(profile));
    if (decoder.status != WirehairK6_Success) { wirehair_k6_free(encoder.codec); return 3; }
    for (id = 6; id < 12; ++id) {
        WirehairK6Result r = wirehair_k6_encode(encoder.codec, id, packet, sizeof(packet));
        WirehairK6Status s = wirehair_k6_decode(decoder.codec, id, packet, (size_t)r.bytes_written);
        if (r.status != WirehairK6_Success || s != (id == 11 ? WirehairK6_Success : WirehairK6_NeedMore)) {
            wirehair_k6_free(encoder.codec); wirehair_k6_free(decoder.codec); return 4;
        }
    }
    wirehair_k6_free(encoder.codec);
    {
        WirehairK6Result r = wirehair_k6_recover(decoder.codec, recovered, sizeof(recovered));
        int failed = r.status != WirehairK6_Success || r.bytes_written != sizeof(source) ||
            memcmp(expected, recovered, sizeof(expected));
        wirehair_k6_free(decoder.codec);
        return failed ? 5 : 0;
    }
}
