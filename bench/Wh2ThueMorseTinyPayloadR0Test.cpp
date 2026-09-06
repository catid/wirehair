#include "Wh2ThueMorseTinyPayloadR0.h"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <vector>

namespace {
using Byte = uint8_t;
using wh2_thue_tiny_payload_r0::Payload;

void Check(bool condition, const char* message)
{
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        std::exit(1);
    }
}

// Polynomial product/reduction, independent of the shared lookup tables.
Byte Multiply(Byte a, Byte b)
{
    unsigned product = 0;
    for (unsigned bit = 0; bit < 8; ++bit)
        if (b & (1u << bit)) product ^= static_cast<unsigned>(a) << bit;
    for (int bit = 14; bit >= 8; --bit)
        if (product & (1u << bit)) product ^= 0x14du << (bit - 8);
    return static_cast<Byte>(product);
}

void Exhaustive()
{
    std::array<std::array<Byte, 4>, 6> data;
    for (auto& source : data) source = {{0xa7, 0, 0, 0xa7}};
    const void* sources[6];
    Byte scales[6] = {};
    for (unsigned i = 0; i < 6; ++i) sources[i] = data[i].data() + 1;
    unsigned cases = 0;
    for (unsigned position = 0; position < 6; ++position) {
        for (unsigned coefficient = 0; coefficient < 256; ++coefficient) {
            scales[position] = static_cast<Byte>(coefficient);
            for (unsigned value = 0; value < 256; ++value) {
                data[position][1] = static_cast<Byte>(value);
                data[position][2] = static_cast<Byte>(255 - value);
                const auto saved = data;
                std::array<Byte, 4> output = {{0xa7, 0x35, 0x8c, 0xa7}};
                Payload(output.data() + 1, sources, scales, 6, 2);
                Check(output[1] == Multiply(data[position][1], scales[position]) &&
                      output[2] == Multiply(data[position][2], scales[position]),
                      "exhaustive coefficient/input at every source position");
                Check(output.front() == 0xa7 && output.back() == 0xa7 && data == saved,
                      "exhaustive exact output span and immutable input");
                ++cases;
            }
        }
        scales[position] = 0;
    }
    Check(cases == 393216, "complete exhaustive roster");
}

void MixedAndFallback()
{
    const unsigned widths[] = {1, 2, 3, 15, 16, 17, 31, 32, 33, 63, 64, 65, 1280};
    const unsigned counts[] = {1, 5, 6, 7, 16, 17};
    for (unsigned count : counts) for (unsigned width : widths) {
        for (unsigned trial = 0; trial < 16; ++trial) {
            std::vector<std::vector<Byte>> data(count, std::vector<Byte>(width + 34, 0xa7));
            std::vector<const void*> sources(count);
            std::vector<Byte> scales(count);
            for (unsigned source = 0; source < count; ++source) {
                for (unsigned byte = 0; byte < width; ++byte)
                    data[source][17 + byte] = static_cast<Byte>(
                        source * 71 + byte * 29 + trial * 13 + (byte >> 2));
                // Include all-zero, all-unit, mixed-zero/unit and dense scales.
                scales[source] = trial < 2 ? static_cast<Byte>(trial) :
                    (trial == 2 ? static_cast<Byte>(source % 3) :
                     static_cast<Byte>(source * 43 + trial * 17));
                sources[source] = data[source].data() + 17;
            }
            // Read-only source spans may coincide; the destination stays disjoint.
            if (trial == 15)
                std::fill(sources.begin(), sources.end(), sources.front());
            const auto saved = data;
            std::vector<Byte> baseline(width + 34, 0xa7), candidate = baseline;
            gf256_mulset_multi_mem(baseline.data() + 17, sources.data(), scales.data(),
                                   static_cast<int>(count), static_cast<int>(width));
            Payload(candidate.data() + 17, sources.data(), scales.data(),
                    static_cast<int>(count), static_cast<int>(width));
            Check(candidate == baseline, "generic fallback/fast path exact parity");
            for (unsigned byte = 0; byte < width; ++byte) {
                Byte expected = 0;
                for (unsigned source = 0; source < count; ++source)
                    expected ^= Multiply(static_cast<const Byte*>(sources[source])[byte], scales[source]);
                Check(candidate[17 + byte] == expected, "mixed independent polynomial oracle");
            }
            Check(std::count(candidate.begin(), candidate.begin() + 17, 0xa7) == 17 &&
                  std::count(candidate.end() - 17, candidate.end(), 0xa7) == 17 && data == saved,
                  "mixed guards and immutable input");
        }
    }
    Byte untouched[2] = {0x35, 0x8c};
    // Literal no-op shapes also keep the compiler from diagnosing the excluded
    // positive/positive case as a possible null read after loop unrolling.
    Payload(untouched, nullptr, nullptr, -1, -1);
    Payload(untouched, nullptr, nullptr, -1, 0);
    Payload(untouched, nullptr, nullptr, -1, 2);
    Payload(untouched, nullptr, nullptr, 0, -1);
    Payload(untouched, nullptr, nullptr, 0, 0);
    Payload(untouched, nullptr, nullptr, 0, 2);
    Payload(untouched, nullptr, nullptr, 6, -1);
    Payload(untouched, nullptr, nullptr, 6, 0);
    Check(untouched[0] == 0x35 && untouched[1] == 0x8c, "nonpositive work is a no-op");
}

void ExactSpans()
{
    // Allocator redzones in the sanitizer build catch reads/writes beyond the
    // two-byte promise; readable canaries alone cannot establish that bound.
    for (unsigned offset : {0u, 1u}) {
        std::unique_ptr<Byte[]> data[6];
        const void* sources[6];
        Byte scales[6] = {0, 1, 37, 89, 171, 255};
        Byte expected[2] = {};
        for (unsigned source = 0; source < 6; ++source) {
            data[source].reset(new Byte[2 + offset]);
            if (offset) data[source][0] = 0xa7;
            sources[source] = data[source].get() + offset;
            for (unsigned byte = 0; byte < 2; ++byte) {
                const Byte value = static_cast<Byte>(source * 71 + byte * 37 + 13);
                data[source][offset + byte] = value;
                expected[byte] ^= Multiply(value, scales[source]);
            }
        }
        std::unique_ptr<Byte[]> output(new Byte[2 + offset]);
        if (offset) output[0] = 0xa7;
        Payload(output.get() + offset, sources, scales, 6, 2);
        Check(output[offset] == expected[0] && output[offset + 1] == expected[1],
              "exact-size input/output oracle");
        if (offset) Check(output[0] == 0xa7, "unaligned exact-size output prefix");
    }
}
} // namespace

int main()
{
    Check(gf256_init() == 0 && GF256Ctx.Polynomial == 0x14d, "shared ordinary GF256 runtime");
    Exhaustive();
    MixedAndFallback();
    ExactSpans();
    std::cout << "PASS 393216 exhaustive, 1248 mixed/fallback and 2 exact-span cases; no timing or candidate data\n";
}
