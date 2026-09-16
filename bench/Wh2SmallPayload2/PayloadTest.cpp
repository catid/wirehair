#include "Payload.h"
#include <algorithm>
#include <array>
#include <climits>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <vector>

static void Check(bool ok) { if (!ok) std::abort(); }
static uint8_t Multiply(unsigned a, unsigned b)
{
    unsigned p = 0;
    for (unsigned i = 0; i < 8; ++i) if (b & (1u << i)) p ^= a << i;
    for (int i = 14; i >= 8; --i) if (p & (1u << i)) p ^= 0x14du << (i - 8);
    return static_cast<uint8_t>(p);
}

template<unsigned K> static void Test()
{
    using wh2_small_payload2::Payload;
    std::array<std::array<uint8_t, 2>, K> data = {};
    const void* sources[K];
    uint8_t scales[K] = {};
    for (unsigned i = 0; i < K; ++i) sources[i] = data[i].data();
    for (unsigned pos = 0; pos < K; ++pos) {
        for (unsigned s = 0; s < 256; ++s) {
            scales[pos] = static_cast<uint8_t>(s);
            for (unsigned x = 0; x < 256; ++x) {
                data[pos] = {{static_cast<uint8_t>(x), static_cast<uint8_t>(x ^ 255)}};
                uint8_t out[4] = {173, 0, 0, 173};
                Payload<K>(out + 1, sources, scales, 2);
                Check(out[0] == 173 && out[3] == 173);
                Check(out[1] == Multiply(x, s) && out[2] == Multiply(x ^ 255, s));
            }
        }
        scales[pos] = 0;
    }
    for (unsigned bytes : {1u,2u,3u,15u,16u,17u,31u,32u,33u,63u,64u,65u,1280u})
        for (unsigned trial = 0; trial < 16; ++trial) {
            std::vector<std::vector<uint8_t>> input(K, std::vector<uint8_t>(bytes + 2, 173));
            for (unsigned i = 0; i < K; ++i) {
                scales[i] = trial < 2 ? static_cast<uint8_t>(trial) :
                    trial == 2 ? static_cast<uint8_t>(i % 3) : static_cast<uint8_t>(i * 43 + trial * 17);
                for (unsigned j = 0; j < bytes; ++j) input[i][j+1] = static_cast<uint8_t>(i*71+j*29+trial*13);
                sources[i] = input[i].data() + 1;
            }
            if (trial == 15) std::fill(sources, sources + K, sources[0]);
            const auto saved = input;
            std::vector<uint8_t> out(bytes + 2, 173), reference(bytes + 2, 173);
            Payload<K>(out.data() + 1, sources, scales, bytes);
            wirehair_k6_payload::Payload(reference.data() + 1, sources, scales, K, bytes);
            Check(out == reference && input == saved && out.front() == 173 && out.back() == 173);
            for (unsigned j = 0; j < bytes; ++j) {
                uint8_t expected = 0;
                for (unsigned i = 0; i < K; ++i)
                    expected ^= Multiply(static_cast<const uint8_t*>(sources[i])[j], scales[i]);
                Check(out[j+1] == expected);
            }
        }
    for (unsigned offset : {0u, 1u}) {
        // ASan exact-sized allocations detect reads as well as writes past B=2.
        std::unique_ptr<uint8_t[]> input[K], output(new uint8_t[2+offset]);
        for (unsigned i = 0; i < K; ++i) {
            input[i].reset(new uint8_t[2+offset]);
            std::fill(input[i].get(), input[i].get()+2+offset, static_cast<uint8_t>(i+1));
            sources[i] = input[i].get() + offset;
        }
        Payload<K>(output.get() + offset, sources, scales, 2);
    }
    Payload<K>(nullptr, nullptr, nullptr, 0);
    Payload<K>(nullptr, nullptr, nullptr, -1);
    Payload<K>(nullptr, nullptr, nullptr, INT_MIN);
    std::cout << "PASS K=" << K << " exhaustive=" << K*65536u << " fallback=208\n";
}
int main()
{
    Check(gf256_init() == 0);
    Test<2>(); Test<3>(); Test<4>(); Test<5>(); Test<6>(); Test<8>();
}
