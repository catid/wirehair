// README diagnostic only: actual public APIs, identical packet IDs and payloads.
// Standalone build instructions and frozen workload are in README.md.
#include <wirehair/wirehair.h>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
uint64_t Next(uint64_t& state) {
    uint64_t z = (state += UINT64_C(0x9e3779b97f4a7c15));
    z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31);
}
void Require(bool ok, const char* what) {
    if (!ok) throw std::runtime_error(what);
}
unsigned Number(const char* s) {
    Require(s && *s, "empty number");
    uint64_t n = 0;
    for (; *s; ++s) {
        Require(*s >= '0' && *s <= '9', "nondecimal number");
        n = n * 10 + unsigned(*s - '0');
        Require(n <= 1000000, "number too large");
    }
    return static_cast<unsigned>(n);
}
struct Handles {
    WirehairCodec e1 = nullptr, d1 = nullptr;
    WirehairV2Codec e2 = nullptr, d2 = nullptr;
    void Free() {
        wirehair_free(e1); wirehair_free(d1);
        wirehair_v2_free(e2); wirehair_v2_free(d2);
        e1 = d1 = nullptr; e2 = d2 = nullptr;
    }
    ~Handles() { Free(); }
};
struct Result { unsigned first = 0, attempt = 0; uint64_t ns = 0, profile = 0; };

Result Run(bool v2, const std::vector<uint8_t>& message, unsigned b,
           const std::vector<uint32_t>& ids, bool timing) {
    // Application buffers are outside the clock for both arms.
    std::vector<uint8_t> packet(b + 2, 0xa5), recovered(message.size() + 2, 0xa5);
    uint8_t descriptor[WIREHAIR_V2_PROFILE_SERIALIZED_BYTES] = {};
    uint32_t descriptor_bytes = 0;
    Handles h;
    Result r;
    const auto start = std::chrono::steady_clock::now();
    if (v2) {
        Require(wirehair_v2_encoder_create(message.data(), message.size(), b,
            descriptor, sizeof(descriptor), &descriptor_bytes, &h.e2) == WirehairV2_Success,
            "WH2 encoder creation failed");
        Require(descriptor_bytes == sizeof(descriptor), "WH2 descriptor length");
        Require(wirehair_v2_decoder_create(descriptor, descriptor_bytes, &h.d2) ==
            WirehairV2_Success, "WH2 decoder creation failed");
    } else {
        Require(wirehair_encoder_create_ex(nullptr, message.data(), message.size(), b,
            &h.e1) == Wirehair_Success, "WH1 encoder creation failed");
        h.d1 = wirehair_decoder_create(nullptr, message.size(), b);
        Require(h.d1 != nullptr, "WH1 decoder creation failed");
    }
    for (unsigned i = 0; i < ids.size(); ++i) {
        uint32_t written = 0;
        bool success;
        if (v2) {
            Require(wirehair_v2_encode(h.e2, ids[i], packet.data() + 1, b, &written) ==
                WirehairV2_Success && written == b, "WH2 encode failed");
            const auto status = wirehair_v2_decode(h.d2, ids[i], packet.data() + 1, written);
            Require(status == WirehairV2_Success || status == WirehairV2_NeedMore,
                "WH2 decode error");
            success = status == WirehairV2_Success;
        } else {
            Require(wirehair_encode(h.e1, ids[i], packet.data() + 1, b, &written) ==
                Wirehair_Success && written == b, "WH1 encode failed");
            const auto status = wirehair_decode(h.d1, ids[i], packet.data() + 1, written);
            Require(status == Wirehair_Success || status == Wirehair_NeedMore,
                "WH1 decode error");
            success = status == Wirehair_Success;
        }
        if (success) { r.first = i + 1; break; }
    }
    if (r.first) {
        if (v2) {
            uint64_t written = 0;
            Require(wirehair_v2_recover(h.d2, recovered.data() + 1, message.size(), &written) ==
                WirehairV2_Success && written == message.size(), "WH2 recover failed");
        } else {
            Require(wirehair_recover(h.d1, recovered.data() + 1, message.size()) ==
                Wirehair_Success, "WH1 recover failed");
        }
    }
    h.Free();
    const auto stop = std::chrono::steady_clock::now();
    if (timing) r.ns = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start).count());
    Require(packet.front() == 0xa5 && packet.back() == 0xa5 && recovered.front() == 0xa5 &&
        recovered.back() == 0xa5, "output guard overwritten");
    if (r.first) Require(std::memcmp(recovered.data() + 1, message.data(), message.size()) == 0,
        "recovered bytes differ");
    if (v2) {
        WirehairV2Profile p = {};
        Require(wirehair_v2_profile_deserialize(descriptor, descriptor_bytes, &p) ==
            WirehairV2_Success && p.profile_id == WIREHAIR_V2_PROFILE_CERTIFIED_2026_07 &&
            p.message_bytes == message.size() && p.block_bytes == b, "unexpected default WH2 profile");
        r.attempt = p.seed_attempt; r.profile = p.profile_id;
    }
    return r;
}
}
int main(int argc, char** argv) {
    try {
        Require(argc == 5, "usage: public_api_sweep timing|recovery K B trials");
        const std::string mode(argv[1]);
        const bool timing = mode == "timing";
        Require(timing || mode == "recovery", "unknown mode");
        const unsigned k = Number(argv[2]), b = Number(argv[3]), trials = Number(argv[4]);
        Require(k >= 8 && k <= 64000 && (b == 64 || b == 1280) && trials >= 1 &&
            trials <= 4096, "dimensions out of range");
        Require(wirehair_init() == Wirehair_Success, "library init failed");
        std::puts("mode,blocks,block_bytes,trial,seed,wh1_first,wh2_first,wh1_ns,wh2_ns,wh2_profile,wh2_attempt");
        // Two unreported warm-up pairs for timing only; no recovery samples discarded.
        for (unsigned round = 0; round < trials + (timing ? 2u : 0u); ++round) {
            const uint64_t seed = UINT64_C(0xa736f91700000000) ^ (uint64_t(k) << 24) ^
                (uint64_t(b) << 8) ^ (uint64_t(round) * UINT64_C(0x9e3779b97f4a7c15));
            std::vector<uint8_t> message(size_t(k) * b);
            uint64_t payload_state = seed;
            for (size_t i = 0; i < message.size(); i += 8) {
                const uint64_t v = Next(payload_state);
                for (unsigned j = 0; j < 8 && i + j < message.size(); ++j)
                    message[i + j] = static_cast<uint8_t>(v >> (j * 8));
            }
            uint64_t loss_state = seed ^ UINT64_C(0x10fade);
            std::vector<uint32_t> ids;
            const unsigned delivered = k + (timing ? 0u : 4u);
            for (uint32_t id = 0; ids.size() < delivered; ++id) {
                Require(id < 4 * delivered + 4096, "loss schedule exceeded bound");
                // Uniform 53-bit variates, exactly specified ten-percent IID threshold.
                if (timing || (Next(loss_state) >> 11) >= UINT64_C(900719925474099))
                    ids.push_back(id);
            }
            Result wh1, wh2;
            if (round % 2 == 0) { wh1 = Run(false, message, b, ids, timing); wh2 = Run(true, message, b, ids, timing); }
            else { wh2 = Run(true, message, b, ids, timing); wh1 = Run(false, message, b, ids, timing); }
            Require((!wh1.first || wh1.first >= k) && (!wh2.first || wh2.first >= k), "early success");
            if (timing) Require(wh1.first == k && wh2.first == k, "no-loss recovery failed");
            if (timing && round < 2) continue;
            std::printf("%s,%u,%u,%u,%llu,%u,%u,%llu,%llu,%llu,%u\n", mode.c_str(), k, b,
                round - (timing ? 2u : 0u), static_cast<unsigned long long>(seed), wh1.first, wh2.first,
                static_cast<unsigned long long>(wh1.ns), static_cast<unsigned long long>(wh2.ns),
                static_cast<unsigned long long>(wh2.profile), wh2.attempt);
        }
        Require(std::fflush(stdout) == 0 && !std::ferror(stdout), "output failure");
        return 0;
    } catch (const std::exception& e) { std::fprintf(stderr, "%s\n", e.what()); return 1; }
}
