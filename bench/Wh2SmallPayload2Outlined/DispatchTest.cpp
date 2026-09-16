#include "Payload.h"
#include <algorithm>
#include <array>
#include <climits>
#include <cstdlib>
#include <iostream>

static void Check(bool ok) { if (!ok) std::abort(); }

int main()
{
    Check(gf256_init() == 0);
    using wh2_small_payload2_outlined::Apply;
    for (int count : {0, -1, INT_MIN})
        for (int bytes : {INT_MIN, -1, 0, 1, 2, 64})
            Apply(nullptr, nullptr, nullptr, count, bytes);
    for (int count : {1, 3, 5, 6, 8, 17})
        for (int bytes : {INT_MIN, -1, 0})
            Apply(nullptr, nullptr, nullptr, count, bytes);

    std::array<std::array<std::uint8_t, 64>, 17> input;
    const void* sources[17];
    std::uint8_t scales[17];
    unsigned cases = 0;
#ifdef WH_COUNT
    std::uint64_t observed_calls = 0;
#endif
    for (int count = 1; count <= 17; ++count)
        for (int bytes : {1, 2, 3, 64})
            for (unsigned pattern = 0; pattern < 4; ++pattern) {
                for (unsigned i = 0; i < input.size(); ++i) {
                    for (unsigned j = 0; j < input[i].size(); ++j)
                        input[i][j] = static_cast<std::uint8_t>(i * 29 + j * 71 + pattern);
                    sources[i] = input[i].data();
                    scales[i] = pattern < 2 ? static_cast<std::uint8_t>(pattern) :
                        pattern == 2 ? static_cast<std::uint8_t>(i % 3) :
                        static_cast<std::uint8_t>(2 + i * 13);
                }
                const auto saved = input;
                std::array<std::uint8_t, 66> candidate, baseline;
                candidate.fill(173); baseline.fill(173);
#ifdef WH_COUNT
                gf256_count_reset();
#endif
                Apply(candidate.data() + 1, sources, scales, count, bytes);
#ifdef WH_COUNT
                std::uint64_t calls[6], work[6];
                for (int op = 0; op < 6; ++op) {
                    calls[op] = gf256_count_calls(op);
                    work[op] = gf256_count_bytes(op);
                    observed_calls += calls[op];
                }
                gf256_count_reset();
#endif
                wirehair_k6_payload::Payload(baseline.data() + 1, sources, scales, count, bytes);
#ifdef WH_COUNT
                for (int op = 0; op < 6; ++op) {
                    Check(calls[op] == gf256_count_calls(op));
                    Check(work[op] == gf256_count_bytes(op));
                }
#endif
                Check(candidate == baseline && input == saved && candidate.front() == 173);
                Check(std::all_of(candidate.begin() + bytes + 1, candidate.end(),
                                  [](std::uint8_t value) { return value == 173; }));
                ++cases;
            }
#ifdef WH_COUNT
    Check(observed_calls > 0); // A disabled counter backend cannot pass vacuously.
#endif
    std::cout << "PASS runtime count dispatch, 36 null no-ops, " << cases << " fallbacks";
#ifdef WH_COUNT
    std::cout << ", six call/byte counters match";
#endif
    std::cout << '\n';
}
