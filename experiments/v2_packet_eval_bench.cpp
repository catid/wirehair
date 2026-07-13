// Reproducible A/B microbenchmark for wirehair-xqb5 and wirehair-x8rs.12.
//
// This benchmarks the cold-solve packet-row build stage and the hot packet-
// evaluation loop.  Run under taskset on an otherwise idle physical core;
// storage and output are reused and warmed, and A/B order rotates between
// samples.
//
//   cmake -S . -B build -DWIREHAIR_BUILD_BENCHMARKS=ON
//   cmake --build build --target wirehair_v2_packet_eval_bench
//   taskset -c 2 build/codec/wirehair_v2_packet_eval_bench 40 2

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

bool RunColdSolveRowBuild(unsigned samples, uint32_t mix_count)
{
    const uint32_t K = 10000u;
    const uint32_t P = 211u;
    wirehair_v2::PacketRowConfig config;
    config.PeelSeed = UINT32_C(0x4d241359);
    config.MixCount = mix_count;
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
        "cold_solve_row_build K=%u P=%u mix_count=%u rows=%u samples=%u "
        "wrapper_median_ms=%.6f cached_median_ms=%.6f "
        "paired_ratio_median=%.6f cache_improvement=%.3f%% "
        "sink=%" PRIu64 "\n",
        K, P, mix_count, K, samples,
        Median(wrapper_samples) * 1000.0,
        Median(cached_samples) * 1000.0,
        Median(ratios),
        (1.0 - Median(ratios)) * 100.0,
        (uint64_t)sink);
    return true;
}

void AddSourcesPaired(
    uint8_t* destination,
    const void* const* sources,
    uint32_t source_count,
    uint32_t block_bytes)
{
    uint32_t source = 0u;
    for (; source + 1u < source_count; source += 2u) {
        gf256_add2_mem(
            destination, sources[source], sources[source + 1u],
            (int)block_bytes);
    }
    if (source < source_count) {
        gf256_add_mem(destination, sources[source], (int)block_bytes);
    }
}

bool RunGatherXorSize(
    uint32_t block_bytes,
    uint32_t source_count,
    unsigned samples)
{
    if (block_bytes == 0u || block_bytes > 0x7fffffffu ||
        source_count < 3u || source_count > 8u ||
        block_bytes > std::numeric_limits<size_t>::max() / source_count)
    {
        return false;
    }
    const size_t source_bytes = (size_t)source_count * block_bytes;
    AlignedStorage source_storage(source_bytes);
    AlignedStorage initial(block_bytes);
    AlignedStorage reference_output(block_bytes);
    AlignedStorage paired_output(block_bytes);
    AlignedStorage gather_output(block_bytes);
    const void* sources[8] = {};
    for (uint32_t source = 0u; source < source_count; ++source) {
        sources[source] = source_storage.Data + (size_t)source * block_bytes;
    }
    for (size_t i = 0u; i < source_bytes; ++i) {
        source_storage.Data[i] = (uint8_t)(i * 193u + (i >> 9) + 0x6du);
    }
    for (uint32_t i = 0u; i < block_bytes; ++i) {
        initial.Data[i] = (uint8_t)(i * 131u + (i >> 3) + 0x27u);
    }

    // Use an independent byte oracle through both repetition parities.
    // Checking only the result after the adaptive repetition count is
    // insufficient: an even count cancels every XOR term and could conceal a
    // broken kernel.
    static const unsigned kOracleRepetitions[] = { 1u, 2u, 3u };
    for (unsigned oracle_repetitions : kOracleRepetitions)
    {
        std::memcpy(reference_output.Data, initial.Data, block_bytes);
        std::memcpy(paired_output.Data, initial.Data, block_bytes);
        std::memcpy(gather_output.Data, initial.Data, block_bytes);
        for (unsigned repetition = 0u;
             repetition < oracle_repetitions;
             ++repetition)
        {
            for (uint32_t source = 0u; source < source_count; ++source) {
                const uint8_t* const input = static_cast<const uint8_t*>(
                    sources[source]);
                for (uint32_t i = 0u; i < block_bytes; ++i) {
                    reference_output.Data[i] ^= input[i];
                }
            }
            AddSourcesPaired(
                paired_output.Data, sources, source_count, block_bytes);
            gf256_add_multi_mem(
                gather_output.Data, sources,
                (int)source_count, (int)block_bytes);
        }
        if (std::memcmp(
                paired_output.Data, reference_output.Data, block_bytes) != 0 ||
            std::memcmp(
                gather_output.Data, reference_output.Data, block_bytes) != 0)
        {
            std::fprintf(stderr,
                "XOR oracle mismatch bb=%u sources=%u repetitions=%u\n",
                block_bytes, source_count, oracle_repetitions);
            return false;
        }
    }

    const uint64_t target_source_bytes = UINT64_C(32) << 20;
    uint64_t repetitions_wide = target_source_bytes /
        ((uint64_t)source_count * block_bytes);
    repetitions_wide = std::max<uint64_t>(1u, repetitions_wide);
    repetitions_wide = std::min<uint64_t>(UINT64_C(262144), repetitions_wide);
    const unsigned repetitions = (unsigned)repetitions_wide;
    volatile uint64_t sink = 0u;
    const auto measure = [&](bool gather) {
        uint8_t* const output = gather ?
            gather_output.Data : paired_output.Data;
        std::memcpy(output, initial.Data, block_bytes);
        const Clock::time_point start = Clock::now();
        for (unsigned repetition = 0u;
             repetition < repetitions; ++repetition)
        {
            if (gather) {
                gf256_add_multi_mem(
                    output, sources, (int)source_count, (int)block_bytes);
            }
            else {
                AddSourcesPaired(
                    output, sources, source_count, block_bytes);
            }
        }
        const double elapsed =
            std::chrono::duration<double>(Clock::now() - start).count();
        sink ^= output[(uint32_t)(
            ((uint64_t)repetitions * 257u + source_count) % block_bytes)];
        return elapsed;
    };

    // Warm both real kernels and verify the algebra before timing.
    (void)measure(false);
    (void)measure(true);
    if (std::memcmp(
            paired_output.Data, gather_output.Data, block_bytes) != 0)
    {
        std::fprintf(stderr,
            "paired/gather XOR mismatch bb=%u sources=%u\n",
            block_bytes, source_count);
        return false;
    }

    std::vector<double> paired_seconds;
    std::vector<double> gather_seconds;
    std::vector<double> ratios;
    paired_seconds.reserve(samples);
    gather_seconds.reserve(samples);
    ratios.reserve(samples);
    for (unsigned sample = 0u; sample < samples; ++sample)
    {
        double paired;
        double gather;
        if ((sample & 1u) == 0u) {
            paired = measure(false);
            gather = measure(true);
        }
        else {
            gather = measure(true);
            paired = measure(false);
        }
        if (std::memcmp(
                paired_output.Data, gather_output.Data, block_bytes) != 0)
        {
            std::fprintf(stderr,
                "timed paired/gather XOR mismatch bb=%u sources=%u "
                "sample=%u\n",
                block_bytes, source_count, sample);
            return false;
        }
        paired_seconds.push_back(paired);
        gather_seconds.push_back(gather);
        ratios.push_back(gather / paired);
    }
    const double ratio = Median(ratios);
    std::printf(
        "paired_gather bb=%u sources=%u reps=%u samples=%u "
        "paired_median_ms=%.6f gather_median_ms=%.6f "
        "gather_over_paired_ratio=%.6f gather_improvement=%.3f%% "
        "sink=%" PRIu64 "\n",
        block_bytes, source_count, repetitions, samples,
        Median(paired_seconds) * 1000.0,
        Median(gather_seconds) * 1000.0,
        ratio, (1.0 - ratio) * 100.0, (uint64_t)sink);
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
    std::vector<uint32_t> rows = FixedRows(K);
    if (config.MixCount == 1u)
    {
        // FixedRows represents the common distribution but happens not to
        // contain the singleton peel shape.  Add one deterministic singleton
        // so the production addset(source, mix) branch is covered by both the
        // byte oracle and the timed workload.
        bool found_singleton = false;
        for (uint32_t id = 0u; id < 1000000u; ++id)
        {
            wirehair::PeelRowParameters params;
            params.Initialize(
                id, config.PeelSeed, (uint16_t)K, (uint16_t)P);
            if (params.PeelCount == 1u)
            {
                if (std::find(rows.begin(), rows.end(), id) == rows.end()) {
                    rows.push_back(id);
                }
                found_singleton = true;
                break;
            }
        }
        if (!found_singleton) {
            std::fprintf(stderr, "mix1 singleton fixture not found\n");
            return false;
        }
    }
    wirehair_v2::PacketRowRuntime runtime;
    if (!runtime.Initialize(K, P, config.MixCount)) {
        return false;
    }

    uint32_t degree_histogram[65] = {};
    uint64_t degree_sum = 0u;
    uint64_t baseline_destination_traffic = 0u;
    uint64_t baseline_total_traffic = 0u;
    uint64_t fused_traffic_reduction = 0u;
    for (uint32_t id : rows) {
        const std::vector<uint32_t> row =
            wirehair_v2::GeneratePacketMatrixRow(K, P, id, config);
        const uint32_t source_degree =
            (uint32_t)row.size() - config.MixCount;
        ++degree_histogram[std::min<uint32_t>(source_degree, 64u)];
        degree_sum += source_degree;
        baseline_destination_traffic +=
            2u * source_degree + 2u * config.MixCount - 1u;
        baseline_total_traffic +=
            3u * source_degree + 3u * config.MixCount - 1u;
        if (config.MixCount == 3u) {
            fused_traffic_reduction += 4u;
        }
        else if (config.MixCount == 2u) {
            fused_traffic_reduction += source_degree == 1u ? 2u : 4u;
        }
        else if (config.MixCount == 1u) {
            fused_traffic_reduction += source_degree <= 2u ? 2u : 4u;
        }
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
    std::vector<double> production_ratios;
    baseline_samples.reserve(samples);
    fused_samples.reserve(samples);
    cached_samples.reserve(samples);
    paired_ratios.reserve(samples);
    cached_ratios.reserve(samples);
    production_ratios.reserve(samples);
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
        production_ratios.push_back(timings[2] / timings[0]);
    }

    std::printf(
        "bb=%u mix_count=%u rows=%zu reps=%u samples=%u "
        "mean_source_degree=%.3f "
        "baseline_median_ms=%.6f fused_median_ms=%.6f cached_median_ms=%.6f "
        "paired_ratio_median=%.6f improvement=%.3f%% "
        "cached_vs_wrapper_ratio=%.6f cache_improvement=%.3f%% "
        "production_ratio_median=%.6f production_improvement=%.3f%% "
        "destination_traffic_reduction=%.3f%% "
        "total_traffic_reduction=%.3f%% sink=%" PRIu64 "\n",
        block_bytes, config.MixCount, rows.size(), repetitions, samples,
        (double)degree_sum / rows.size(),
        Median(baseline_samples) * 1000.0,
        Median(fused_samples) * 1000.0,
        Median(cached_samples) * 1000.0,
        Median(paired_ratios),
        (1.0 - Median(paired_ratios)) * 100.0,
        Median(cached_ratios),
        (1.0 - Median(cached_ratios)) * 100.0,
        Median(production_ratios),
        (1.0 - Median(production_ratios)) * 100.0,
        100.0 * fused_traffic_reduction /
            (double)baseline_destination_traffic,
        100.0 * fused_traffic_reduction /
            (double)baseline_total_traffic,
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
    uint32_t mix_count = wirehair_v2::kCertifiedPacketMixCount;
    bool run_packet = true;
    bool run_gather = false;
    if (argc >= 2) {
        char* end = nullptr;
        const unsigned long parsed = std::strtoul(argv[1], &end, 10);
        if (!end || *end != '\0' || parsed < 30u || parsed > 1000u) {
            std::fprintf(stderr,
                "usage: %s [samples=30..1000] [mix_count=1|2|3] "
                "[mode=packet|gather|all]\n",
                argv[0]);
            return 2;
        }
        samples = (unsigned)parsed;
    }
    if (argc >= 3)
    {
        char* end = nullptr;
        const unsigned long parsed = std::strtoul(argv[2], &end, 10);
        if (!end || *end != '\0' || parsed < 1u || parsed > 3u) {
            std::fprintf(stderr,
                "usage: %s [samples=30..1000] [mix_count=1|2|3] "
                "[mode=packet|gather|all]\n",
                argv[0]);
            return 2;
        }
        mix_count = (uint32_t)parsed;
    }
    if (argc == 4)
    {
        if (std::strcmp(argv[3], "packet") == 0) {
            run_packet = true;
            run_gather = false;
        }
        else if (std::strcmp(argv[3], "gather") == 0) {
            run_packet = false;
            run_gather = true;
        }
        else if (std::strcmp(argv[3], "all") == 0) {
            run_packet = true;
            run_gather = true;
        }
        else {
            std::fprintf(stderr,
                "usage: %s [samples=30..1000] [mix_count=1|2|3] "
                "[mode=packet|gather|all]\n",
                argv[0]);
            return 2;
        }
    }
    else if (argc > 4) {
        std::fprintf(stderr,
            "usage: %s [samples=30..1000] [mix_count=1|2|3] "
            "[mode=packet|gather|all]\n",
            argv[0]);
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
    config.MixCount = mix_count;
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
    if (run_packet) {
        ok = RunColdSolveRowBuild(samples, mix_count) && ok;
    }
    if (run_gather)
    {
        static const uint32_t kGatherSizes[] = {
            1u, 2u, 7u, 8u, 15u, 16u, 17u, 31u, 32u, 33u,
            63u, 64u, 65u, 127u, 128u, 129u,
            255u, 256u, 257u, 511u, 512u, 513u,
            1023u, 1024u, 1025u, 1279u, 1280u, 1281u, 4096u,
            16u * 1024u, 64u * 1024u, 100u * 1024u,
            128u * 1024u - 1u, 128u * 1024u, 128u * 1024u + 1u,
            256u * 1024u, 1024u * 1024u
        };
        for (uint32_t source_count = 3u; source_count <= 8u; ++source_count) {
            for (uint32_t block_bytes : kGatherSizes) {
                ok = RunGatherXorSize(
                    block_bytes, source_count, samples) && ok;
            }
        }
    }
    if (run_packet)
    {
        ok = RunSize(system, config, 1u, 4096u, samples) && ok;
        // The mixed completion candidate requires even block bytes.  Keep a
        // short sub-cacheline sweep so wrapper overhead at bb=1 cannot conceal
        // a regression (or win) in the smallest reachable production packets.
        ok = RunSize(system, config, 2u, 4096u, samples) && ok;
        ok = RunSize(system, config, 8u, 4096u, samples) && ok;
        ok = RunSize(system, config, 32u, 4096u, samples) && ok;
        ok = RunSize(system, config, 1280u, 2048u, samples) && ok;
        ok = RunSize(system, config, 100u * 1024u, 32u, samples) && ok;
        ok = RunSize(system, config, 1024u * 1024u, 3u, samples) && ok;
    }
    return ok ? 0 : 1;
}
