#include "wirehair/wirehair_small.h"
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
    WirehairSmallCreateResult d = wirehair_small_decoder_create(profile, sizeof(profile));
    unsigned id;
    if (d.status != WirehairSmall_Success) return 1;
    for (id = 0; id < 3; ++id) {
        WirehairSmallStatus s = wirehair_small_decode(d.codec, id, source + 2 * id, id == 2 ? 1 : 2);
        if (s != (id == 2 ? WirehairSmall_Success : WirehairSmall_NeedMore)) {
            wirehair_small_free(d.codec); return 1;
        }
    }
    {
        WirehairSmallResult r = wirehair_small_recover(d.codec, recovered, sizeof(recovered));
        int failed = r.status != WirehairSmall_Success || r.bytes_written != sizeof(source) ||
            memcmp(source, recovered, sizeof(source));
        wirehair_small_free(d.codec);
        return failed;
    }
}

/* A separate C translation unit exercises the actual C ABI, not C++ casts. */
int main(void)
{
    unsigned char source[5] = {0, 1, 2, 3, 4};
    const unsigned char expected[5] = {0, 1, 2, 3, 4};
    unsigned char profile[WIREHAIR_SMALL_PROFILE_BYTES], packet[2], recovered[5];
    WirehairSmallCreateResult encoder, decoder;
    unsigned id;
    if (wirehair_init() != Wirehair_Success) return 1;
    if (standalone()) return 6;
    encoder = wirehair_small_encoder_create(source, sizeof(source), 2, WirehairSmall_BorrowedImmutable,
                                    profile, sizeof(profile));
    if (encoder.status != WirehairSmall_Success) return 2;
    if (wirehair_small_profile_validate(profile, sizeof(profile)) != WirehairSmall_Success) {
        wirehair_small_free(encoder.codec); return 7;
    }
    wirehair_small_free(encoder.codec);
    encoder = wirehair_small_encoder_create_profile(source, profile, sizeof(profile),
                                                  WirehairSmall_BorrowedImmutable);
    if (encoder.status != WirehairSmall_Success) return 8;
    if (wirehair_small_encoder_detach_input(encoder.codec) != WirehairSmall_Success) {
        wirehair_small_free(encoder.codec); return 9;
    }
    memset(source, 0, sizeof(source));
    decoder = wirehair_small_decoder_create(profile, sizeof(profile));
    if (decoder.status != WirehairSmall_Success) { wirehair_small_free(encoder.codec); return 3; }
    for (id = 3; id < 6; ++id) {
        WirehairSmallResult r = wirehair_small_encode(encoder.codec, id, packet, sizeof(packet));
        WirehairSmallStatus s = wirehair_small_decode(decoder.codec, id, packet, (size_t)r.bytes_written);
        if (r.status != WirehairSmall_Success || s != (id == 5 ? WirehairSmall_Success : WirehairSmall_NeedMore)) {
            wirehair_small_free(encoder.codec); wirehair_small_free(decoder.codec); return 4;
        }
    }
    wirehair_small_free(encoder.codec);
    {
        WirehairSmallResult r = wirehair_small_recover(decoder.codec, recovered, sizeof(recovered));
        int failed = r.status != WirehairSmall_Success || r.bytes_written != sizeof(source) ||
            memcmp(expected, recovered, sizeof(expected));
        wirehair_small_free(decoder.codec);
        return failed ? 5 : 0;
    }
}
