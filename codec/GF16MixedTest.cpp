#include "WirehairV2GF16.h"

#include <chrono>
#include <cstdio>
#include <cstring>
#include <thread>
#include <vector>

namespace {

bool ConcurrentInitializationTest()
{
    std::vector<uint8_t> results(64u, 0u);
    std::vector<std::thread> workers;
    workers.reserve(results.size());
    for (size_t i = 0; i < results.size(); ++i) {
        workers.push_back(std::thread([i, &results]() {
            results[i] = (uint8_t)(
                wirehair_v2::InitializeGF16() &&
                wirehair_v2::MixedGF16Coefficient(0u, 0u) == 34916u);
        }));
    }
    for (std::thread& worker : workers) worker.join();
    for (uint8_t result : results) if (result != 1u) return false;
    return true;
}

uint64_t Mix64(uint64_t x)
{
    x += UINT64_C(0x9e3779b97f4a7c15);
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

bool FieldTest()
{
    using namespace wirehair_v2;
    if (!InitializeGF16() ||
        MixedGF16Coefficient(0u, 0u) != 34916u ||
        MixedGF16Coefficient(1u, 0u) != 2472u ||
        MixedGF16Coefficient(0u, 243u) != 59155u ||
        MixedGF16Coefficient(0u, 244u) != 34916u ||
        MixedGF16Coefficient(2u, 0u) != 0u)
    {
        return false;
    }
    if (GF16MultiplyInitialized(0x1234u, 0x5678u) !=
            GF16Multiply(0x1234u, 0x5678u) ||
        GF16InverseInitialized(0x1234u) != GF16Inverse(0x1234u))
    {
        return false;
    }
    for (uint32_t x = 1u; x < 65536u; ++x) {
        if (GF16Multiply((uint16_t)x, GF16Inverse((uint16_t)x)) != 1u)
            return false;
    }
    for (uint32_t trial = 0; trial < 100000u; ++trial)
    {
        const uint16_t a = (uint16_t)Mix64(3u * trial);
        const uint16_t b = (uint16_t)Mix64(3u * trial + 1u);
        const uint16_t c = (uint16_t)Mix64(3u * trial + 2u);
        if (GF16Multiply(a, b) != GF16Multiply(b, a) ||
            GF16Multiply(GF16Multiply(a, b), c) !=
                GF16Multiply(a, GF16Multiply(b, c)) ||
            GF16Multiply(a, (uint16_t)(b ^ c)) !=
                (uint16_t)(GF16Multiply(a, b) ^ GF16Multiply(a, c)))
        {
            return false;
        }
    }
    return true;
}

bool KernelTest()
{
    using namespace wirehair_v2;
    const uint32_t lengths[] = {2u, 4u, 30u, 128u, 1280u};
    const uint16_t scales[] = {
        0u, 1u, 0x00d3u, 0x0100u, 0x0101u, 0xc753u, 0xffffu
    };
    for (uint32_t bytes : lengths) for (uint16_t scale : scales)
    {
        std::vector<uint8_t> source(bytes + 2u, 0xa5u);
        std::vector<uint8_t> direct(bytes + 2u, 0x5au);
        std::vector<uint8_t> planar = direct;
        for (uint32_t i = 0; i < bytes; ++i) {
            source[i + 1u] = (uint8_t)Mix64(i + bytes);
            direct[i + 1u] = planar[i + 1u] =
                (uint8_t)Mix64(i + scale);
        }
        std::vector<uint8_t> source_low(bytes / 2u);
        std::vector<uint8_t> source_high(bytes / 2u);
        std::vector<uint8_t> destination_low(bytes / 2u);
        std::vector<uint8_t> destination_high(bytes / 2u);
        if (!GF16Deinterleave(
                source.data() + 1u,
                source_low.data(), source_high.data(), bytes) ||
            !GF16Deinterleave(
                planar.data() + 1u,
                destination_low.data(), destination_high.data(), bytes) ||
            !GF16MulAddMem(
                direct.data() + 1u, scale, source.data() + 1u, bytes) ||
            !GF16MulAddPlanar(
                destination_low.data(), destination_high.data(), scale,
                source_low.data(), source_high.data(), bytes / 2u) ||
            !GF16Interleave(
                destination_low.data(), destination_high.data(),
                planar.data() + 1u, bytes) ||
            direct != planar ||
            direct.front() != 0x5au || direct.back() != 0x5au)
        {
            return false;
        }

        std::vector<uint8_t> scaled = source;
        std::vector<uint8_t> oracle(bytes + 2u, 0u);
        std::vector<uint8_t> scale_low(bytes / 2u);
        std::vector<uint8_t> scale_high(bytes / 2u);
        std::vector<uint8_t> scale_scratch(bytes / 2u, 0xa5u);
        if (!GF16ScaleMem(scaled.data() + 1u, scale, bytes) ||
            !GF16MulAddMem(
                oracle.data() + 1u, scale, source.data() + 1u, bytes) ||
            !GF16Deinterleave(
                source.data() + 1u,
                scale_low.data(), scale_high.data(), bytes) ||
            !GF16ScalePlanar(
                scale_low.data(), scale_high.data(), scale,
                scale_scratch.data(), bytes / 2u) ||
            !GF16Interleave(
                scale_low.data(), scale_high.data(),
                planar.data() + 1u, bytes) ||
            std::memcmp(
                scaled.data() + 1u, oracle.data() + 1u, bytes) != 0 ||
            std::memcmp(
                scaled.data() + 1u, planar.data() + 1u, bytes) != 0)
        {
            return false;
        }
    }

    for (uint32_t bytes : lengths) for (uint16_t scale0 : scales)
        for (uint16_t scale1 : scales)
    {
        const uint32_t elements = bytes / 2u;
        std::vector<uint8_t> source_low(elements);
        std::vector<uint8_t> source_high(elements);
        std::vector<uint8_t> separate0_low(elements);
        std::vector<uint8_t> separate0_high(elements);
        std::vector<uint8_t> separate1_low(elements);
        std::vector<uint8_t> separate1_high(elements);
        for (uint32_t i = 0; i < elements; ++i) {
            source_low[i] = (uint8_t)Mix64(i + bytes);
            source_high[i] = (uint8_t)Mix64(i + bytes + 1u);
            separate0_low[i] = (uint8_t)Mix64(i + scale0);
            separate0_high[i] = (uint8_t)Mix64(i + scale0 + 1u);
            separate1_low[i] = (uint8_t)Mix64(i + scale1 + 2u);
            separate1_high[i] = (uint8_t)Mix64(i + scale1 + 3u);
        }
        std::vector<uint8_t> paired0_low = separate0_low;
        std::vector<uint8_t> paired0_high = separate0_high;
        std::vector<uint8_t> paired1_low = separate1_low;
        std::vector<uint8_t> paired1_high = separate1_high;
        if (!GF16MulAddPlanar(
                separate0_low.data(), separate0_high.data(), scale0,
                source_low.data(), source_high.data(), elements) ||
            !GF16MulAddPlanar(
                separate1_low.data(), separate1_high.data(), scale1,
                source_low.data(), source_high.data(), elements) ||
            !GF16MulAddPlanar2(
                paired0_low.data(), paired0_high.data(), scale0,
                paired1_low.data(), paired1_high.data(), scale1,
                source_low.data(), source_high.data(), elements) ||
            separate0_low != paired0_low ||
            separate0_high != paired0_high ||
            separate1_low != paired1_low ||
            separate1_high != paired1_high)
        {
            return false;
        }
    }

    uint8_t source[4] = {1u, 2u, 3u, 4u};
    uint8_t destination[4] = {9u, 8u, 7u, 6u};
    const uint8_t snapshot[4] = {9u, 8u, 7u, 6u};
    if (GF16MulAddMem(destination, 7u, source, 3u) ||
        GF16MulAddMem(nullptr, 7u, source, 4u) ||
        GF16MulAddMem(destination, 7u, nullptr, 4u) ||
        GF16MulAddMem(destination, 7u, source, 0u) ||
        GF16ScaleMem(destination, 7u, 3u) ||
        GF16Deinterleave(source, destination, destination + 2u, 3u) ||
        GF16Interleave(source, source + 2u, destination, 3u) ||
        std::memcmp(destination, snapshot, sizeof(snapshot)) != 0)
    {
        return false;
    }

    uint8_t alias[8] = {1u, 2u, 3u, 4u, 5u, 6u, 7u, 8u};
    uint8_t alias_oracle[8];
    uint8_t alias_source[8];
    std::memcpy(alias_oracle, alias, sizeof(alias));
    std::memcpy(alias_source, alias, sizeof(alias));
    if (!GF16MulAddMem(
            alias_oracle, 0xc753u, alias_source, sizeof(alias)) ||
        !GF16MulAddMem(alias, 0xc753u, alias, sizeof(alias)) ||
        std::memcmp(alias, alias_oracle, sizeof(alias)) != 0)
    {
        return false;
    }
    uint8_t overlap[16];
    for (uint32_t i = 0; i < sizeof(overlap); ++i) overlap[i] = (uint8_t)i;
    uint8_t overlap_snapshot[16];
    std::memcpy(overlap_snapshot, overlap, sizeof(overlap));
    if (GF16MulAddMem(overlap + 1u, 7u, overlap, 8u) ||
        GF16Deinterleave(overlap, overlap, overlap + 8u, 8u) ||
        GF16Interleave(overlap, overlap + 4u, overlap, 8u) ||
        GF16MulAddPlanar(
            overlap, overlap + 4u, 7u,
            overlap + 2u, overlap + 12u, 4u) ||
        GF16MulAddPlanar2(
            overlap, overlap + 4u, 7u,
            overlap + 8u, overlap + 12u, 11u,
            overlap + 2u, overlap + 6u, 4u) ||
        GF16ScalePlanar(
            overlap, overlap + 4u, 7u, overlap + 2u, 4u) ||
        std::memcmp(overlap, overlap_snapshot, sizeof(overlap)) != 0)
    {
        return false;
    }

    // Endpoint-adjacent planes are valid, while every invalid pair-kernel
    // shape below must be rejected before a zero scale can short-circuit.
    uint8_t adjacent[24];
    for (uint32_t i = 0; i < sizeof(adjacent); ++i) {
        adjacent[i] = (uint8_t)(0x80u + i);
    }
    uint8_t adjacent_snapshot[24];
    std::memcpy(adjacent_snapshot, adjacent, sizeof(adjacent));
    if (!GF16MulAddPlanar2(
            adjacent, adjacent + 4u, 0u,
            adjacent + 8u, adjacent + 12u, 0u,
            adjacent + 16u, adjacent + 20u, 4u) ||
        GF16MulAddPlanar2(
            nullptr, adjacent + 4u, 0u,
            adjacent + 8u, adjacent + 12u, 0u,
            adjacent + 16u, adjacent + 20u, 4u) ||
        GF16MulAddPlanar2(
            adjacent, adjacent + 4u, 0u,
            adjacent + 8u, nullptr, 0u,
            adjacent + 16u, adjacent + 20u, 4u) ||
        GF16MulAddPlanar2(
            adjacent, adjacent + 4u, 0u,
            adjacent + 8u, adjacent + 12u, 0u,
            nullptr, adjacent + 20u, 4u) ||
        GF16MulAddPlanar2(
            adjacent, adjacent + 4u, 0u,
            adjacent + 8u, adjacent + 12u, 0u,
            adjacent + 16u, nullptr, 4u) ||
        GF16MulAddPlanar2(
            adjacent, adjacent + 4u, 0u,
            adjacent + 8u, adjacent + 12u, 0u,
            adjacent + 16u, adjacent + 20u, 0u) ||
        GF16MulAddPlanar2(
            adjacent, adjacent + 4u, 0u,
            adjacent + 8u, adjacent + 12u, 0u,
            adjacent + 16u, adjacent + 20u, UINT32_C(0x80000000)) ||
        GF16MulAddPlanar2(
            adjacent, adjacent + 2u, 0u,
            adjacent + 8u, adjacent + 12u, 0u,
            adjacent + 16u, adjacent + 20u, 4u) ||
        GF16MulAddPlanar2(
            adjacent, adjacent + 4u, 0u,
            adjacent + 8u, adjacent + 16u, 0u,
            adjacent + 16u, adjacent + 20u, 4u) ||
        GF16MulAddPlanar2(
            adjacent, adjacent + 4u, 0u,
            adjacent + 8u, adjacent + 12u, 0u,
            adjacent + 16u, adjacent + 16u, 4u) ||
        std::memcmp(adjacent, adjacent_snapshot, sizeof(adjacent)) != 0)
    {
        return false;
    }
    return true;
}

} // namespace

int main(int argc, char** argv)
{
    const bool benchmark = argc == 2 &&
        std::strcmp(argv[1], "--benchmark") == 0;
    if (argc > 2 || (argc == 2 && !benchmark)) {
        std::fprintf(stderr, "usage: GF16MixedTest [--benchmark]\n");
        return 2;
    }
    if (!ConcurrentInitializationTest() || !FieldTest() || !KernelTest()) {
        std::fprintf(stderr, "GF(2^16) mixed-field test failed\n");
        return 1;
    }
    std::printf(
        "GF(2^16) mixed-field test passed; coefficients=34916/2472/59155\n");
    if (benchmark)
    {
        const uint32_t iterations = 10000001u;
        volatile uint16_t safe_sink = 0u;
        volatile uint16_t ready_sink = 0u;
        const auto safe_start = std::chrono::steady_clock::now();
        for (uint32_t i = 0; i < iterations; ++i)
            safe_sink ^= wirehair_v2::GF16Multiply((uint16_t)i, 0xc753u);
        const auto safe_stop = std::chrono::steady_clock::now();
        const auto ready_start = std::chrono::steady_clock::now();
        for (uint32_t i = 0; i < iterations; ++i)
            ready_sink ^= wirehair_v2::GF16MultiplyInitialized(
                (uint16_t)i, 0xc753u);
        const auto ready_stop = std::chrono::steady_clock::now();
        const double safe_seconds = std::chrono::duration<double>(
            safe_stop - safe_start).count();
        const double ready_seconds = std::chrono::duration<double>(
            ready_stop - ready_start).count();
        std::printf(
            "GF16 scalar benchmark safe=%.6fs initialized=%.6fs ratio=%.3f "
            "checksums=%u/%u\n",
            safe_seconds, ready_seconds, safe_seconds / ready_seconds,
            (unsigned)safe_sink, (unsigned)ready_sink);
    }
    return 0;
}
