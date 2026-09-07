#include "Wh2K6Serialized.h"
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
    Wh2K6CreateResult d = wh2_k6_decoder_create(profile, sizeof(profile));
    unsigned id;
    if (d.status != Wh2K6_Success) return 1;
    for (id = 0; id < 6; ++id) {
        Wh2K6Status s = wh2_k6_decode(d.codec, id, source + 2 * id, id == 5 ? 1 : 2);
        if (s != (id == 5 ? Wh2K6_Success : Wh2K6_NeedMore)) {
            wh2_k6_free(d.codec); return 1;
        }
    }
    {
        Wh2K6Result r = wh2_k6_recover(d.codec, recovered, sizeof(recovered));
        int failed = r.status != Wh2K6_Success || r.bytes_written != sizeof(source) ||
            memcmp(source, recovered, sizeof(source));
        wh2_k6_free(d.codec);
        return failed;
    }
}

/* A separate C translation unit exercises the actual C ABI, not C++ casts. */
int main(void)
{
    unsigned char source[11] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    unsigned char profile[WH2_K6_PROFILE_BYTES], packet[2], recovered[11];
    Wh2K6CreateResult encoder, decoder;
    unsigned id;
    if (wirehair_init() != Wirehair_Success) return 1;
    if (standalone()) return 6;
    encoder = wh2_k6_encoder_create(source, sizeof(source), 2, Wh2K6_Independent,
                                    profile, sizeof(profile));
    if (encoder.status != Wh2K6_Success) return 2;
    decoder = wh2_k6_decoder_create(profile, sizeof(profile));
    if (decoder.status != Wh2K6_Success) { wh2_k6_free(encoder.codec); return 3; }
    for (id = 0; id < 6; ++id) {
        Wh2K6Result r = wh2_k6_encode(encoder.codec, id, packet, sizeof(packet));
        Wh2K6Status s = wh2_k6_decode(decoder.codec, id, packet, (size_t)r.bytes_written);
        if (r.status != Wh2K6_Success || s != (id == 5 ? Wh2K6_Success : Wh2K6_NeedMore)) {
            wh2_k6_free(encoder.codec); wh2_k6_free(decoder.codec); return 4;
        }
    }
    wh2_k6_free(encoder.codec);
    {
        Wh2K6Result r = wh2_k6_recover(decoder.codec, recovered, sizeof(recovered));
        int failed = r.status != Wh2K6_Success || r.bytes_written != sizeof(source) ||
            memcmp(source, recovered, sizeof(source));
        wh2_k6_free(decoder.codec);
        return failed ? 5 : 0;
    }
}
