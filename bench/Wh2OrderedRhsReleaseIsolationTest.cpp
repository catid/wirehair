// Untimed parity of the two public APIs in one process. No scientific mode.
#include "OrderedReleaseApi.h"
#include <algorithm>
#include <array>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <vector>

namespace {
void Check(bool value, const char* why) {
    if (!value) throw std::runtime_error(why);
}
struct Api {
    decltype(&wirehair_v2_encoder_create_with_options) create;
    decltype(&wirehair_v2_encode) encode;
    decltype(&wirehair_v2_decoder_create) decoder;
    decltype(&wirehair_v2_decode) feed;
    decltype(&wirehair_v2_recover) recover;
    decltype(&wirehair_v2_encoder_detach_input) detach;
    decltype(&wirehair_v2_free) free;
};
const Api apis[] = {
    {wirehair_v2_encoder_create_with_options, wirehair_v2_encode, wirehair_v2_decoder_create,
     wirehair_v2_decode, wirehair_v2_recover, wirehair_v2_encoder_detach_input, wirehair_v2_free},
    {wh2_ordered_release_wirehair_v2_encoder_create_with_options,
     wh2_ordered_release_wirehair_v2_encode, wh2_ordered_release_wirehair_v2_decoder_create,
     wh2_ordered_release_wirehair_v2_decode, wh2_ordered_release_wirehair_v2_recover,
     wh2_ordered_release_wirehair_v2_encoder_detach_input, wh2_ordered_release_wirehair_v2_free}
};
struct Owner {
    const Api& api;
    WirehairV2Codec handle = nullptr;
    explicit Owner(const Api& a) : api(a) {}
    Owner(const Owner&) = delete;
    Owner& operator=(const Owner&) = delete;
    ~Owner() { api.free(handle); }
};
using Profile = std::array<uint8_t, 32>;
void Guards(const std::vector<uint8_t>& data, size_t bytes) {
    for (size_t i=0; i<data.size(); ++i)
        if (i < 16 || i >= 16+bytes) Check(data[i] == 0xa5, "output guards");
}
void Shape(unsigned k, unsigned width, unsigned tail, unsigned policy) {
    const size_t bytes = size_t(k)*width-tail;
    std::vector<uint8_t> source(bytes);
    for (size_t i=0; i<bytes; ++i) source[i] = uint8_t(37*i+i/11);
    const auto original = source;
    Owner encoder0(apis[0]), encoder1(apis[1]);
    Owner* encoders[] = {&encoder0, &encoder1};
    Profile profiles[2] = {};
    WirehairV2EncoderOptions options = WIREHAIR_V2_ENCODER_OPTIONS_INIT;
    options.source_policy = policy;
    for (unsigned a=0; a<2; ++a) {
        uint32_t written = 0;
        Check(apis[a].create(source.data(), bytes, width, &options, profiles[a].data(),
              32, &written, &encoders[a]->handle) == WirehairV2_Success &&
              written == 32 && encoders[a]->handle, "encoder create");
    }
    Check(profiles[0] == profiles[1], "profile parity");
    if (policy == WirehairV2EncoderSource_Independent) std::fill(source.begin(), source.end(), 0);
    std::vector<uint32_t> ids;
    for (unsigned id=0; id<2*k; ++id) ids.push_back(id);
    for (unsigned id=0; id<k; ++id) ids.push_back(UINT32_MAX-2*id);
    std::vector<std::vector<uint8_t>> packets;
    for (uint32_t id : ids) {
        std::vector<uint8_t> out[2] = {std::vector<uint8_t>(width+32, 0xa5),
                                     std::vector<uint8_t>(width+32, 0xa5)};
        const unsigned expected = id == k-1 ? width-tail : width;
        for (unsigned a=0; a<2; ++a) {
            uint32_t written = 0;
            Check(apis[a].encode(encoders[a]->handle, id, out[a].data()+16, expected-1,
                  &written) == WirehairV2_BufferTooSmall && written == expected, "short encode");
            Guards(out[a], 0);
            Check(apis[a].encode(encoders[a]->handle, id, out[a].data()+16, width, &written) ==
                  WirehairV2_Success && written == expected, "encode");
            Guards(out[a], expected);
        }
        Check(out[0] == out[1], "packet parity");
        packets.emplace_back(out[0].begin()+16, out[0].begin()+16+expected);
    }
    // Detach makes both ownership policies independent, without changing bytes.
    for (unsigned a=0; a<2; ++a)
        Check(apis[a].detach(encoders[a]->handle) == WirehairV2_Success, "detach");
    std::fill(source.begin(), source.end(), 0);
    for (unsigned a=0; a<2; ++a) {
        std::vector<uint8_t> packet(width); uint32_t written = 0;
        Check(apis[a].encode(encoders[a]->handle, 0, packet.data(), width, &written) ==
              WirehairV2_Success && written == width && packet == packets[0], "detached packet");
    }
    for (unsigned family=0; family<2; ++family) {
        Owner decoder0(apis[0]), decoder1(apis[1]);
        Owner* decoders[] = {&decoder0, &decoder1};
        for (unsigned a=0; a<2; ++a)
            Check(apis[a].decoder(profiles[a].data(), 32, &decoders[a]->handle) ==
                  WirehairV2_Success && decoders[a]->handle, "decoder create");
        bool success = false;
        for (unsigned step=0; step<2*k; ++step) {
            const unsigned slot = step < k ? (family+1)*k+step : step-k;
            WirehairV2Result statuses[2];
            for (unsigned a=0; a<2; ++a)
                statuses[a] = apis[a].feed(decoders[a]->handle, ids[slot], packets[slot].data(),
                                          uint32_t(packets[slot].size()));
            Check(statuses[0] == statuses[1] && (statuses[0] == WirehairV2_NeedMore ||
                  statuses[0] == WirehairV2_Success), "feed parity");
            if (step == 0) {
                for (unsigned a=0; a<2; ++a)
                    Check(apis[a].feed(decoders[a]->handle, ids[slot], packets[slot].data(),
                          uint32_t(packets[slot].size())) == statuses[a], "duplicate feed");
            }
            if (statuses[0] == WirehairV2_Success) { success = true; break; }
        }
        Check(success, "bounded systematic completion");
        for (unsigned a=0; a<2; ++a) {
            std::vector<uint8_t> out(bytes+32, 0xa5); uint64_t written = 0;
            Check(apis[a].recover(decoders[a]->handle, out.data()+16, bytes-1, &written) ==
                  WirehairV2_BufferTooSmall && written == bytes, "short recovery");
            Guards(out, 0);
            Check(apis[a].recover(decoders[a]->handle, out.data()+16, bytes, &written) ==
                  WirehairV2_Success && written == bytes &&
                  std::equal(original.begin(), original.end(), out.begin()+16), "recovered bytes");
            Guards(out, bytes);
        }
    }
}
} // namespace

int main(int argc, char** argv) {
    try {
        Check(argc == 2 && std::strcmp(argv[1], "--neutral") == 0, "neutral mode required");
        Check(apis[0].create != apis[1].create && apis[0].feed != apis[1].feed &&
              apis[0].free != apis[1].free, "separate public APIs");
        Check(wirehair_init() == Wirehair_Success, "single shared initialization");
        unsigned shapes = 0;
        for (unsigned k : {3u, 6u, 17u, 128u}) for (unsigned width : {2u, 64u, 1280u})
            for (unsigned tail : {0u, 1u}) for (unsigned policy : {1u, 2u}) {
                Shape(k, width, tail, policy); ++shapes;
            }
        std::printf("PASS %u neutral public B/C shapes; no timing or recovery-rate claim\n", shapes);
        return 0;
    } catch (const std::exception& error) {
        std::fprintf(stderr, "FAIL: %s\n", error.what()); return 1;
    }
}
