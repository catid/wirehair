#include "Wh2SmallSerialized.h"
#include "wirehair/wirehair.h"
#include <string.h>

/* Runs before any encoder has existed in this process. Literal descriptor and
 * systematic packets suffice; no fixture builder can hide private state. */
static int standalone(void)
{
    const unsigned char profile[32] = {
        'W','H','K','3',1,0,32,0, 0x31,0x4d,0x54,0x33,0x4b,0x32,0x48,0x57,
        5,0,0,0,0,0,0,0, 2,0,0,0,0,0,0,0
    };
    const unsigned char source[5] = {0,1,2,3,4};
    unsigned char recovered[5];
    Wh2SmallCreateResult d = wh2_small_decoder_create(profile, sizeof(profile));
    unsigned id;
    if (d.status != Wh2Small_Success) return 1;
    for (id = 0; id < 3; ++id) {
        Wh2SmallStatus s = wh2_small_decode(d.codec, id, source + 2 * id, id == 2 ? 1 : 2);
        if (s != (id == 2 ? Wh2Small_Success : Wh2Small_NeedMore)) {
            wh2_small_free(d.codec); return 1;
        }
    }
    {
        Wh2SmallResult r = wh2_small_recover(d.codec, recovered, sizeof(recovered));
        int failed = r.status != Wh2Small_Success || r.bytes_written != sizeof(source) ||
            memcmp(source, recovered, sizeof(source));
        wh2_small_free(d.codec);
        return failed;
    }
}
/* A separate C translation unit exercises the actual C ABI, not C++ casts. */
int main(void)
{
    unsigned char source[5] = {0, 1, 2, 3, 4};
    unsigned char profile[WH2_SMALL_PROFILE_BYTES], packet[2], recovered[5];
    Wh2SmallCreateResult encoder, decoder;
    unsigned id;
    if (wirehair_init() != Wirehair_Success) return 1;
    if (standalone()) return 6;
    encoder = wh2_small_encoder_create(source, sizeof(source), 2, Wh2Small_Independent,
                                    profile, sizeof(profile));
    if (encoder.status != Wh2Small_Success) return 2;
    decoder = wh2_small_decoder_create(profile, sizeof(profile));
    if (decoder.status != Wh2Small_Success) { wh2_small_free(encoder.codec); return 3; }
    for (id = 0; id < 3; ++id) {
        Wh2SmallResult r = wh2_small_encode(encoder.codec, id, packet, sizeof(packet));
        Wh2SmallStatus s = wh2_small_decode(decoder.codec, id, packet, (size_t)r.bytes_written);
        if (r.status != Wh2Small_Success || s != (id == 2 ? Wh2Small_Success : Wh2Small_NeedMore)) {
            wh2_small_free(encoder.codec); wh2_small_free(decoder.codec); return 4;
        }
    }
    wh2_small_free(encoder.codec);
    {
        Wh2SmallResult r = wh2_small_recover(decoder.codec, recovered, sizeof(recovered));
        int failed = r.status != Wh2Small_Success || r.bytes_written != sizeof(source) ||
            memcmp(source, recovered, sizeof(source));
        wh2_small_free(decoder.codec);
        return failed ? 5 : 0;
    }
}
