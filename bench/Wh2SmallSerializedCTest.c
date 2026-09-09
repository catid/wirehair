#include "Wh2SmallSerialized.h"
#include "wirehair/wirehair.h"
#include <string.h>
#define K WH2_SMALL_CODEC_K
#define MESSAGE_BYTES (2 * K - 1)

/* Runs before any encoder has existed in this process. Literal descriptor and
 * systematic packets suffice; no fixture builder can hide private state. */
static int standalone(void)
{
    const unsigned char profile[32] = {
        'W','H','K','0'+K,1,0,32,0, 0x31,0x4d,0x54,0x30+K,0x4b,0x32,0x48,0x57,
        MESSAGE_BYTES,0,0,0,0,0,0,0, 2,0,0,0,0,0,0,0
    };
    unsigned char source[MESSAGE_BYTES], recovered[MESSAGE_BYTES];
    Wh2SmallCreateResult d = wh2_small_decoder_create(profile, sizeof(profile));
    unsigned id;
    for (id = 0; id < sizeof(source); ++id) source[id] = (unsigned char)id;
    if (d.status != Wh2Small_Success) return 1;
    for (id = 0; id < K; ++id) {
        Wh2SmallStatus s = wh2_small_decode(d.codec, id, source + 2 * id, id == K - 1 ? 1 : 2);
        if (s != (id == K - 1 ? Wh2Small_Success : Wh2Small_NeedMore)) {
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
    unsigned char source[MESSAGE_BYTES], recovered[MESSAGE_BYTES];
    unsigned char profile[WH2_SMALL_PROFILE_BYTES], packets[K][2];
    Wh2SmallCreateResult encoder, decoder;
    unsigned id;
    if (wirehair_init() != Wirehair_Success) return 1;
    if (standalone()) return 6;
    for (id = 0; id < sizeof(source); ++id) source[id] = (unsigned char)id;
    encoder = wh2_small_encoder_create(source, sizeof(source), 2, Wh2Small_Independent,
                                    profile, sizeof(profile));
    if (encoder.status != Wh2Small_Success) return 2;
    for (id = 0; id < K; ++id) {
        Wh2SmallResult r = wh2_small_encode(encoder.codec, id, packets[id], sizeof(packets[id]));
        if (r.status != Wh2Small_Success || r.bytes_written != (id == K - 1 ? 1u : 2u)) {
            wh2_small_free(encoder.codec); return 4;
        }
    }
    wh2_small_free(encoder.codec);
    decoder = wh2_small_decoder_create(profile, sizeof(profile));
    if (decoder.status != Wh2Small_Success) return 3;
    for (id = 0; id < K; ++id) {
        Wh2SmallStatus s = wh2_small_decode(decoder.codec, id, packets[id], id == K - 1 ? 1 : 2);
        if (s != (id == K - 1 ? Wh2Small_Success : Wh2Small_NeedMore)) {
            wh2_small_free(decoder.codec); return 4;
        }
    }
    {
        Wh2SmallResult r = wh2_small_recover(decoder.codec, recovered, sizeof(recovered));
        int failed = r.status != Wh2Small_Success || r.bytes_written != sizeof(source) ||
            memcmp(source, recovered, sizeof(source));
        wh2_small_free(decoder.codec);
        return failed ? 5 : 0;
    }
}
