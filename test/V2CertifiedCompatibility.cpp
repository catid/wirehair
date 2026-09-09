// Read-only compatibility replay: compile this same source against the pinned
// pre-admission archive and candidate archive, then compare complete stdout.
// It is deterministic byte evidence, not a timing or recovery-rate workload.
#include <wirehair/wirehair.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {
void Check(bool condition)
{
    if (!condition) { std::fputs("Certified compatibility replay failed\n", stderr); std::exit(1); }
}
void Emit(const void* p, size_t n) { Check(std::fwrite(p, 1, n, stdout) == n); }
void Number(uint32_t n)
{
    unsigned char bytes[4];
    for (unsigned i = 0; i < 4; ++i) bytes[i] = static_cast<unsigned char>(n >> (8 * i));
    Emit(bytes, sizeof(bytes));
}
}
int main()
{
    Check(wirehair_init() == Wirehair_Success);
    for (uint32_t k : {2u,3u,4u,6u,8u,16u,64u,128u})
        for (uint32_t b : {2u,64u,1280u}) for (uint32_t tail : {1u,b}) {
            const size_t message = size_t(k - 1) * b + tail;
            std::vector<unsigned char> source(message);
            for (size_t i = 0; i < message; ++i) source[i] = static_cast<unsigned char>(37 * i + i / 11);
            const auto original = source;
            unsigned char profile[32];
            uint32_t bytes = 0;
            WirehairV2Codec encoder = nullptr, decoder = nullptr;
            Check(wirehair_v2_encoder_create_profile_id(WIREHAIR_V2_PROFILE_CERTIFIED_2026_07,
                source.data(), message, b, profile, sizeof(profile), &bytes, &encoder) == WirehairV2_Success);
            Check(bytes == sizeof(profile));
            Emit(profile, sizeof(profile));
            std::fill(source.begin(), source.end(), 0xcc);
            Check(wirehair_v2_decoder_create(profile, sizeof(profile), &decoder) == WirehairV2_Success);
            std::vector<unsigned char> packet(b);
            for (uint32_t id = 0; id < k + 8; ++id) {
                Check(wirehair_v2_encode(encoder, id, packet.data(), b, &bytes) == WirehairV2_Success);
                Number(id); Number(bytes); Emit(packet.data(), bytes);
                if (id < k) {
                    auto result = wirehair_v2_decode(decoder, id, packet.data(), bytes);
                    Check(result == (id + 1 == k ? WirehairV2_Success : WirehairV2_NeedMore));
                    Number(static_cast<uint32_t>(result));
                }
            }
            for (uint32_t id : {1023u,1024u,131071u,131072u,16777215u,16777216u,UINT32_MAX}) {
                Check(wirehair_v2_encode(encoder, id, packet.data(), b, &bytes) == WirehairV2_Success);
                Number(id); Number(bytes); Emit(packet.data(), bytes);
            }
            uint64_t recovered_bytes = 0;
            for (unsigned repeat = 0; repeat < 2; ++repeat) {
                Check(wirehair_v2_recover(decoder, source.data(), source.size(), &recovered_bytes) ==
                    WirehairV2_Success && recovered_bytes == source.size() && source == original);
                Emit(source.data(), source.size());
            }
            wirehair_v2_free(encoder);
            wirehair_v2_free(decoder);
        }
    Check(std::fflush(stdout) == 0);
}
