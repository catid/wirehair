// Reproducible A/B microbenchmark for wirehair-xqb5 and wirehair-x8rs.12.
//
// This benchmarks the cold-solve packet-row build stage and the hot packet-
// evaluation loop.  Run under taskset on an otherwise idle physical core;
// storage and output are reused and warmed, and A/B order rotates between
// samples.
//
//   cmake -S . -B build -DWIREHAIR_BUILD_BENCHMARKS=ON
//   cmake --build build --target wirehair_v2_packet_eval_bench
//   taskset -c 2 build/codec/wirehair_v2_packet_eval_bench 40

#include "../codec/WirehairV2Precode.h"
#include "../codec/WirehairV2Solve.h"
#include "../WirehairTools.h"
#include "../gf256.h"

#include <algorithm>
#include <chrono>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

struct AlignedStorage
{
    explicit AlignedStorage(size_t bytes)
        : Raw(bytes + 63u)
    {
        const uintptr_t address = reinterpret_cast<uintptr_t>(Raw.data());
        Data = reinterpret_cast<uint8_t*>((address + 63u) & ~uintptr_t(63u));
    }

    std::vector<uint8_t> Raw;
    uint8_t* Data = nullptr;
};

bool EvaluateBaseline(
    const wirehair_v2::PrecodeSystem& system,
    const wirehair_v2::PacketRowConfig& config,
    const uint8_t* intermediate,
    uint32_t block_bytes,
    uint32_t block_id,
    uint8_t* output)
{
    const uint32_t K = system.Params.BlockCount;
    const uint32_t P = system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    wirehair::PeelRowParameters params;
    params.Initialize(block_id, config.PeelSeed, (uint16_t)K, (uint16_t)P);
    wirehair::PeelRowIterator source(
        params, (uint16_t)K, wirehair::NextPrime16((uint16_t)K));
    const wirehair::RowMixIterator mix(
        params, (uint16_t)P, wirehair::NextPrime16((uint16_t)P));

    std::memcpy(
        output,
        intermediate + (size_t)source.GetColumn() * block_bytes,
        block_bytes);
    while (source.Iterate()) {
        gf256_add_mem(
            output,
            intermediate + (size_t)source.GetColumn() * block_bytes,
            (int)block_bytes);
    }
    for (uint32_t i = 0; i < config.MixCount; ++i) {
        gf256_add_mem(
            output,
            intermediate + (size_t)(K + mix.Columns[i]) * block_bytes,
            (int)block_bytes);
    }
    return true;
}

bool EvaluateFused(
    const wirehair_v2::PrecodeSystem& system,
    const wirehair_v2::PacketRowConfig& config,
    const uint8_t* intermediate,
    uint32_t block_bytes,
    uint32_t block_id,
    uint8_t* output)
{
    return wirehair_v2::EvaluatePacketBlockForValidatedSystem(
        system, config, intermediate, block_bytes, block_id, output);
}

bool EvaluateCached(
    const wirehair_v2::PrecodeSystem& system,
    const wirehair_v2::PacketRowConfig& config,
    const wirehair_v2::PacketRowRuntime& runtime,
    const uint8_t* intermediate,
    uint32_t block_bytes,
    uint32_t block_id,
    uint8_t* output)
{
    return wirehair_v2::EvaluatePacketBlockForValidatedSystemWithRuntime(
        system, config, runtime, intermediate, block_bytes, block_id, output);
}

double Median(std::vector<double> values)
{
    std::sort(values.begin(), values.end());
    const size_t middle = values.size() / 2u;
    return values.size() % 2u ? values[middle] :
        (values[middle - 1u] + values[middle]) * 0.5;
}

std::vector<uint32_t> FixedRows(uint32_t K)
{
    std::vector<uint32_t> rows;
    rows.reserve(96u);
    for (uint32_t i = 0; i < 32u; ++i) {
        rows.push_back(i);
        rows.push_back(K + i * 7919u);
        rows.push_back(UINT32_MAX - i * 104729u);
    }
    return rows;
}

bool RunColdSolveRowBuild(unsigned samples)
{
    const uint32_t K = 10000u;
    const uint32_t P = 211u;
    wirehair_v2::PacketRowConfig config;
    config.PeelSeed = UINT32_C(0x4d241359);
    config.MixCount = wirehair_v2::kCertifiedPacketMixCount;
    wirehair_v2::PacketRowRuntime runtime;
    if (!runtime.Initialize(K, P, config.MixCount)) {
        return false;
    }

    volatile uint64_t sink = 0u;
    const auto measure = [&](bool cached) {
        const Clock::time_point start = Clock::now();
        for (uint32_t id = 0u; id < K; ++id)
        {
            const std::vector<uint32_t> row = cached ?
                wirehair_v2::GeneratePacketMatrixRowWithRuntime(
                    K, P, id, config, runtime) :
                wirehair_v2::GeneratePacketMatrixRow(K, P, id, config);
            sink += row.size() + row.front();
        }
        return std::chrono::duration<double>(Clock::now() - start).count();
    };
    (void)measure(false);
    (void)measure(true);

    std::vector<double> wrapper_samples;
    std::vector<double> cached_samples;
    std::vector<double> ratios;
    wrapper_samples.reserve(samples);
    cached_samples.reserve(samples);
    ratios.reserve(samples);
    for (unsigned sample = 0u; sample < samples; ++sample)
    {
        double wrapper_seconds;
        double cached_seconds;
        if ((sample & 1u) == 0u) {
            wrapper_seconds = measure(false);
            cached_seconds = measure(true);
        }
        else {
            cached_seconds = measure(true);
            wrapper_seconds = measure(false);
        }
        wrapper_samples.push_back(wrapper_seconds);
        cached_samples.push_back(cached_seconds);
        ratios.push_back(cached_seconds / wrapper_seconds);
    }
    std::printf(
        "cold_solve_row_build K=%u P=%u rows=%u samples=%u "
        "wrapper_median_ms=%.6f cached_median_ms=%.6f "
        "paired_ratio_median=%.6f cache_improvement=%.3f%% "
        "sink=%" PRIu64 "\n",
        K, P, K, samples,
        Median(wrapper_samples) * 1000.0,
        Median(cached_samples) * 1000.0,
        Median(ratios),
        (1.0 - Median(ratios)) * 100.0,
        (uint64_t)sink);
    return true;
}

bool RunSize(
    const wirehair_v2::PrecodeSystem& system,
    const wirehair_v2::PacketRowConfig& config,
    uint32_t block_bytes,
    unsigned repetitions,
    unsigned samples)
{
    const uint32_t K = system.Params.BlockCount;
    const uint32_t P = system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    const size_t total_bytes = (size_t)(K + P) * block_bytes;
    AlignedStorage intermediate(total_bytes);
    AlignedStorage baseline_output(block_bytes);
    AlignedStorage fused_output(block_bytes);
    AlignedStorage cached_output(block_bytes);
    for (size_t i = 0; i < total_bytes; ++i) {
        intermediate.Data[i] = (uint8_t)(i * 131u + (i >> 13) + 17u);
    }
    const std::vector<uint32_t> rows = FixedRows(K);
    wirehair_v2::PacketRowRuntime runtime;
    if (!runtime.Initialize(K, P, config.MixCount)) {
        return false;
    }

    uint32_t degree_histogram[65] = {};
    uint64_t degree_sum = 0u;
    for (uint32_t id : rows) {
        const std::vector<uint32_t> row =
            wirehair_v2::GeneratePacketMatrixRow(K, P, id, config);
        const uint32_t source_degree =
            (uint32_t)row.size() - config.MixCount;
        ++degree_histogram[std::min<uint32_t>(source_degree, 64u)];
        degree_sum += source_degree;
        EvaluateBaseline(
            system, config, intermediate.Data, block_bytes, id,
            baseline_output.Data);
        EvaluateFused(
            system, config, intermediate.Data, block_bytes, id,
            fused_output.Data);
        EvaluateCached(
            system, config, runtime, intermediate.Data, block_bytes, id,
            cached_output.Data);
        if (std::memcmp(
                baseline_output.Data, fused_output.Data, block_bytes) != 0 ||
            std::memcmp(
                fused_output.Data, cached_output.Data, block_bytes) != 0)
        {
            std::fprintf(stderr,
                "packet evaluation mismatch bb=%u id=%u degree=%u\n",
                block_bytes, id, source_degree);
            return false;
        }
    }

    volatile uint64_t sink = 0u;
    const auto measure = [&](unsigned mode) {
        const Clock::time_point start = Clock::now();
        for (unsigned rep = 0; rep < repetitions; ++rep) {
            for (size_t i = 0; i < rows.size(); ++i) {
                uint8_t* output = mode == 0u ? baseline_output.Data :
                    (mode == 1u ? fused_output.Data : cached_output.Data);
                if (mode == 2u) {
                    EvaluateCached(
                        system, config, runtime, intermediate.Data,
                        block_bytes, rows[i], output);
                }
                else if (mode == 1u) {
                    EvaluateFused(
                        system, config, intermediate.Data, block_bytes,
                        rows[i], output);
                }
                else {
                    EvaluateBaseline(
                        system, config, intermediate.Data, block_bytes,
                        rows[i], output);
                }
                sink ^= output[(i * 257u + rep) % block_bytes];
            }
        }
        return std::chrono::duration<double>(Clock::now() - start).count();
    };

    // Warm all code paths and every reusable source/output page.
    (void)measure(0u);
    (void)measure(1u);
    (void)measure(2u);

    std::vector<double> baseline_samples;
    std::vector<double> fused_samples;
    std::vector<double> cached_samples;
    std::vector<double> paired_ratios;
    std::vector<double> cached_ratios;
    baseline_samples.reserve(samples);
    fused_samples.reserve(samples);
    cached_samples.reserve(samples);
    paired_ratios.reserve(samples);
    cached_ratios.reserve(samples);
    for (unsigned sample = 0; sample < samples; ++sample) {
        double timings[3] = {};
        if (sample % 3u == 0u) {
            timings[0] = measure(0u);
            timings[1] = measure(1u);
            timings[2] = measure(2u);
        }
        else if (sample % 3u == 1u) {
            timings[1] = measure(1u);
            timings[2] = measure(2u);
            timings[0] = measure(0u);
        }
        else {
            timings[2] = measure(2u);
            timings[0] = measure(0u);
            timings[1] = measure(1u);
        }
        baseline_samples.push_back(timings[0]);
        fused_samples.push_back(timings[1]);
        cached_samples.push_back(timings[2]);
        paired_ratios.push_back(timings[1] / timings[0]);
        cached_ratios.push_back(timings[2] / timings[1]);
    }

    std::printf(
        "bb=%u rows=%zu reps=%u samples=%u mean_source_degree=%.3f "
        "baseline_median_ms=%.6f fused_median_ms=%.6f cached_median_ms=%.6f "
        "paired_ratio_median=%.6f improvement=%.3f%% "
        "cached_vs_wrapper_ratio=%.6f cache_improvement=%.3f%% "
        "destination_traffic_reduction=%.3f%% "
        "total_traffic_reduction=%.3f%% sink=%" PRIu64 "\n",
        block_bytes, rows.size(), repetitions, samples,
        (double)degree_sum / rows.size(),
        Median(baseline_samples) * 1000.0,
        Median(fused_samples) * 1000.0,
        Median(cached_samples) * 1000.0,
        Median(paired_ratios),
        (1.0 - Median(paired_ratios)) * 100.0,
        Median(cached_ratios),
        (1.0 - Median(cached_ratios)) * 100.0,
        100.0 * (double)(4u * rows.size()) /
            (double)(2u * degree_sum + 5u * rows.size()),
        100.0 * (double)(4u * rows.size()) /
            (double)(3u * degree_sum + 8u * rows.size()),
        (uint64_t)sink);
    std::printf("degree_histogram:");
    for (unsigned degree = 0; degree < 65u; ++degree) {
        if (degree_histogram[degree] != 0u) {
            std::printf(" %u:%u", degree, degree_histogram[degree]);
        }
    }
    std::printf("\n");
    return true;
}

} // namespace

int main(int argc, char** argv)
{
    unsigned samples = 40u;
    if (argc == 2) {
        char* end = nullptr;
        const unsigned long parsed = std::strtoul(argv[1], &end, 10);
        if (!end || *end != '\0' || parsed < 30u || parsed > 1000u) {
            std::fprintf(stderr, "usage: %s [samples>=30]\n", argv[0]);
            return 2;
        }
        samples = (unsigned)parsed;
    }
    else if (argc != 1) {
        std::fprintf(stderr, "usage: %s [samples>=30]\n", argv[0]);
        return 2;
    }
    if (gf256_init() != 0) {
        return 1;
    }

    const uint32_t K = 128u;
    wirehair_v2::PrecodeSystem system;
    if (!wirehair_v2::BuildPrecodeSystem(
            wirehair_v2::MakeCertifiedParams(
                K, UINT64_C(0x786f72667573696f)),
            system))
    {
        return 1;
    }
    wirehair_v2::PacketRowConfig config;
    config.PeelSeed = UINT32_C(0x4d241359);
    config.MixCount = wirehair_v2::kCertifiedPacketMixCount;
    if (!wirehair_v2::ValidatePrecodeSystem(system) ||
        !wirehair_v2::IsPacketRowDomainValid(
            K,
            system.Params.Staircase + system.Params.DenseRows +
                system.Params.HeavyRows,
            config.MixCount))
    {
        return 1;
    }

    bool ok = true;
    ok = RunColdSolveRowBuild(samples) && ok;
    ok = RunSize(system, config, 1u, 4096u, samples) && ok;
    ok = RunSize(system, config, 1280u, 2048u, samples) && ok;
    ok = RunSize(system, config, 100u * 1024u, 32u, samples) && ok;
    ok = RunSize(system, config, 1024u * 1024u, 3u, samples) && ok;
    return ok ? 0 : 1;
}
