#include "../WirehairHeavy.h"

#include <cstdint>
#include <cstdio>
#include <cstring>

namespace {

uint64_t NextRandom(uint64_t& state)
{
    state ^= state >> 12;
    state ^= state << 25;
    state ^= state >> 27;
    return state * UINT64_C(2685821657736338717);
}

bool CheckWindow(
    unsigned offset,
    unsigned bits,
    uint8_t coefficient,
    uint64_t& state)
{
    uint8_t actual[20];
    uint8_t expected[20];
    for (unsigned i = 0; i < sizeof(actual); ++i) {
        actual[i] = static_cast<uint8_t>(NextRandom(state));
    }
    std::memcpy(expected, actual, sizeof(actual));

    for (unsigned i = 0; i < 4; ++i) {
        if ((bits & (1u << i)) != 0) {
            expected[offset + i] ^= coefficient;
        }
    }
    wirehair::AddHeavyWindow(actual + offset, bits, coefficient);

    if (std::memcmp(actual, expected, sizeof(actual)) == 0) {
        return true;
    }
    std::fprintf(stderr,
        "heavy window mismatch: offset=%u bits=%u coefficient=%u\n",
        offset, bits, coefficient);
    return false;
}

} // namespace

int main()
{
    uint64_t state = UINT64_C(0x9e3779b97f4a7c15);

    // Exhaust every four-bit window and every possible GF(256) coefficient.
    // Varying the offset exercises every native alignment plus guard bytes.
    for (unsigned offset = 0; offset < 8; ++offset) {
        for (unsigned bits = 0; bits < 16; ++bits) {
            for (unsigned coefficient = 0; coefficient < 256; ++coefficient) {
                if (!CheckWindow(
                        offset, bits, static_cast<uint8_t>(coefficient), state)) {
                    return 1;
                }
            }
        }
    }

    // Supplement the exhaustive primitive domain with randomized initial data
    // and deliberately out-of-range high bits (which the helper must ignore).
    for (unsigned trial = 0; trial < 1000000; ++trial) {
        const uint64_t random = NextRandom(state);
        const unsigned offset = static_cast<unsigned>(random % 13u);
        const unsigned bits = static_cast<unsigned>(random >> 8);
        const uint8_t coefficient = static_cast<uint8_t>(random >> 40);
        if (!CheckWindow(offset, bits & 15u, coefficient, state)) {
            return 1;
        }

        uint8_t low[4] = {};
        uint8_t high[4] = {};
        wirehair::AddHeavyWindow(low, bits & 15u, coefficient);
        wirehair::AddHeavyWindow(high, bits, coefficient);
        if (std::memcmp(low, high, sizeof(low)) != 0) {
            std::fprintf(stderr,
                "heavy window high-bit mismatch: bits=%u coefficient=%u\n",
                bits, static_cast<unsigned>(coefficient));
            return 1;
        }
    }

    std::puts("heavy window exhaustive/property tests passed");
    return 0;
}
