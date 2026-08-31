// Reproducible A/B microbenchmark for wirehair-xqb5 and wirehair-x8rs.12.
//
// This benchmarks the cold-solve packet-row build stage and the hot packet-
// evaluation loop.  Run under taskset on an otherwise idle physical core;
// storage and output are reused and warmed, and A/B order rotates between
// samples.
//
//   cmake -S . -B build -DBUILD_TESTS=ON -DCMAKE_BUILD_TYPE=Release -DWIREHAIR_BUILD_BENCHMARKS=ON
//   cmake --build build --target wirehair_v2_packet_eval_bench
//   taskset -c 2 build/codec/wirehair_v2_packet_eval_bench 40 2
//
// Frozen five-term MIX3 same-binary screens (one K/B cell per process):
//
//   taskset -c 2 <bench> 32 3 count5-packet 10000 1280
//   taskset -c 2 <bench> 12 3 count5-recover 10000 1280

#include "../codec/WirehairV2Precode.h"
#include "../codec/WirehairV2PrecodeDecode.h"
#include "../codec/WirehairV2Solve.h"
#include "../WirehairTools.h"
#include "../gf256.h"

#include <algorithm>
#include <chrono>
#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>
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
    uint8_t* output,
    uint64_t* block_ops_out = nullptr)
{
    return wirehair_v2::EvaluatePacketBlockForValidatedSystemWithRuntime(
        system, config, runtime, intermediate, block_bytes, block_id, output,
        block_ops_out);
}

#if defined(_MSC_VER)
#define WH2_PACKET_BENCH_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define WH2_PACKET_BENCH_NOINLINE __attribute__((noinline))
#else
#define WH2_PACKET_BENCH_NOINLINE
#endif

bool PacketBenchRangesOverlap(
    const void* first,
    size_t first_bytes,
    const void* second,
    size_t second_bytes)
{
    const uintptr_t first_begin = reinterpret_cast<uintptr_t>(first);
    const uintptr_t second_begin = reinterpret_cast<uintptr_t>(second);
    const uintptr_t limit = std::numeric_limits<uintptr_t>::max();
    if (first_bytes > limit - first_begin ||
        second_bytes > limit - second_begin)
    {
        return true;
    }
    return first_begin < second_begin + second_bytes &&
        second_begin < first_begin + first_bytes;
}

// Narrow packet-only candidate: initialize from two terms exactly as the
// production evaluator does, then consume the remaining source/mix tail two
// terms per destination pass.  This deliberately uses only the existing
// addset/add2/add kernels; it is not the rejected general gather integration.
WH2_PACKET_BENCH_NOINLINE bool EvaluateTailPaired(
    const wirehair_v2::PrecodeSystem& system,
    const wirehair_v2::PacketRowConfig& config,
    const wirehair_v2::PacketRowRuntime& runtime,
    const uint8_t* intermediate,
    uint32_t block_bytes,
    uint32_t block_id,
    uint8_t* output,
    uint64_t* block_ops_out = nullptr)
{
    const uint32_t K = system.Params.BlockCount;
    const uint64_t P_wide = (uint64_t)system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    if (!intermediate || !output || block_bytes == 0u ||
        block_bytes > 0x7fffffffu || P_wide > UINT32_MAX ||
        !runtime.IsValidFor(K, (uint32_t)P_wide, config.MixCount))
    {
        return false;
    }
    const uint64_t intermediate_bytes_wide =
        ((uint64_t)K + P_wide) * block_bytes;
    if (intermediate_bytes_wide >
            (uint64_t)std::numeric_limits<size_t>::max() ||
        PacketBenchRangesOverlap(
            output,
            block_bytes,
            intermediate,
            (size_t)intermediate_bytes_wide))
    {
        return false;
    }
    const uint32_t P = (uint32_t)P_wide;
    wirehair::PeelRowParameters params;
    params.Initialize(block_id, config.PeelSeed, (uint16_t)K, (uint16_t)P);
    wirehair::PeelRowIterator source(
        params, (uint16_t)K, runtime.SourcePrime());

    uint16_t mix_columns[wirehair::RowMixIterator::kColumnCount] = {};
    if (config.MixCount == 1u) {
        mix_columns[0] = params.MixFirst;
    }
    else
    {
        const wirehair::RowMixIterator mix(
            params, (uint16_t)P, runtime.PrecodePrime());
        for (uint32_t i = 0u; i < config.MixCount; ++i) {
            mix_columns[i] = mix.Columns[i];
        }
    }

    const uint8_t* const first_source =
        intermediate + (size_t)source.GetColumn() * block_bytes;
    uint32_t first_pending_mix = 0u;
    if (source.Iterate())
    {
        gf256_addset_mem(
            output,
            first_source,
            intermediate + (size_t)source.GetColumn() * block_bytes,
            (int)block_bytes);
    }
    else
    {
        gf256_addset_mem(
            output,
            first_source,
            intermediate + (size_t)(K + mix_columns[0]) * block_bytes,
            (int)block_bytes);
        first_pending_mix = 1u;
    }

    const uint8_t* pending = nullptr;
    const auto consume = [&](const uint8_t* term) {
        if (pending)
        {
            gf256_add2_mem(output, pending, term, (int)block_bytes);
            pending = nullptr;
        }
        else {
            pending = term;
        }
    };
    while (source.Iterate()) {
        consume(intermediate + (size_t)source.GetColumn() * block_bytes);
    }
    for (uint32_t i = first_pending_mix; i < config.MixCount; ++i) {
        consume(intermediate + (size_t)(K + mix_columns[i]) * block_bytes);
    }
    if (pending) {
        gf256_add_mem(output, pending, (int)block_bytes);
    }
    if (block_ops_out) {
        *block_ops_out = (uint64_t)params.PeelCount + config.MixCount;
    }
    return true;
}

#undef WH2_PACKET_BENCH_NOINLINE

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
    unsigned samples,
    uint32_t target_source_degree = 0u,
    bool cached_only = false,
    size_t row_limit = 0u)
{
    const uint32_t K = system.Params.BlockCount;
    const uint32_t P = system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    const uint64_t total_bytes_wide =
        ((uint64_t)K + P) * block_bytes;
    if (total_bytes_wide >
        (uint64_t)std::numeric_limits<size_t>::max())
    {
        return false;
    }
    const size_t total_bytes = (size_t)total_bytes_wide;
    AlignedStorage intermediate(total_bytes);
    AlignedStorage baseline_output(block_bytes);
    AlignedStorage fused_output(block_bytes);
    AlignedStorage cached_output(block_bytes);
    AlignedStorage tail_paired_output(block_bytes);
    for (size_t i = 0; i < total_bytes; ++i) {
        intermediate.Data[i] = (uint8_t)(i * 131u + (i >> 13) + 17u);
    }
    std::vector<uint32_t> rows;
    if (target_source_degree != 0u)
    {
        for (uint32_t id = 0u; id < 1000000u && rows.size() < 256u; ++id)
        {
            wirehair::PeelRowParameters params;
            params.Initialize(
                id, config.PeelSeed, (uint16_t)K, (uint16_t)P);
            if (params.PeelCount == target_source_degree) {
                rows.push_back(id);
            }
        }
        if (rows.size() != 256u) {
            std::fprintf(stderr,
                "degree fixture shortfall K=%u degree=%u rows=%zu\n",
                K, target_source_degree, rows.size());
            return false;
        }
    }
    else {
        rows = FixedRows(K);
    }
    if (config.MixCount == 1u && target_source_degree == 0u)
    {
        // FixedRows represents the common distribution but happens not to
        // contain the singleton peel shape.  Add one deterministic singleton
        // so the production addset(source, mix) branch is covered by both the
        // byte oracle and the timed workload.
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
                break;
            }
        }
        // Production deliberately removes weight-one rows above its tuned K
        // threshold.  The K=128/1000 screens cover the singleton branch;
        // larger-K screens correctly retain only their reachable degrees.
    }
    if (row_limit != 0u && rows.size() > row_limit) {
        rows.resize(row_limit);
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
        uint64_t cached_ops = 0u;
        uint64_t tail_paired_ops = 0u;
        EvaluateCached(
            system, config, runtime, intermediate.Data, block_bytes, id,
            cached_output.Data, &cached_ops);
        EvaluateTailPaired(
            system, config, runtime, intermediate.Data, block_bytes, id,
            tail_paired_output.Data, &tail_paired_ops);
        if (std::memcmp(
                baseline_output.Data, fused_output.Data, block_bytes) != 0 ||
            std::memcmp(
                fused_output.Data, cached_output.Data, block_bytes) != 0 ||
            std::memcmp(
                cached_output.Data, tail_paired_output.Data, block_bytes) != 0 ||
            cached_ops != tail_paired_ops)
        {
            std::fprintf(stderr,
                "packet evaluation mismatch bb=%u id=%u degree=%u "
                "cached_ops=%" PRIu64 " tail_paired_ops=%" PRIu64 "\n",
                block_bytes, id, source_degree,
                cached_ops, tail_paired_ops);
            return false;
        }
    }

    volatile uint64_t sink = 0u;
    const auto measure = [&](unsigned mode) {
        const Clock::time_point start = Clock::now();
        for (unsigned rep = 0; rep < repetitions; ++rep) {
            for (size_t i = 0; i < rows.size(); ++i) {
                uint8_t* output = mode == 0u ? baseline_output.Data :
                    (mode == 1u ? fused_output.Data :
                        (mode == 2u ? cached_output.Data :
                            tail_paired_output.Data));
                if (mode == 3u) {
                    EvaluateTailPaired(
                        system, config, runtime, intermediate.Data,
                        block_bytes, rows[i], output);
                }
                else if (mode == 2u) {
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

    // Cross-binary production A/B runs need only the cached evaluator after
    // the byte/work oracle above has checked all implementations.  Avoid
    // spending three quarters of a long fallback run timing reference paths.
    if (cached_only)
    {
        (void)measure(2u);
        std::vector<double> cached_samples;
        cached_samples.reserve(samples);
        for (unsigned sample = 0u; sample < samples; ++sample) {
            cached_samples.push_back(measure(2u));
        }
        std::printf(
            "cached_only bb=%u K=%u mix_count=%u "
            "target_source_degree=%u rows=%zu reps=%u samples=%u "
            "mean_source_degree=%.3f cached_median_ms=%.6f "
            "sink=%" PRIu64 "\n",
            block_bytes, K, config.MixCount, target_source_degree,
            rows.size(), repetitions, samples,
            (double)degree_sum / rows.size(),
            Median(cached_samples) * 1000.0, (uint64_t)sink);
        return true;
    }

    // Warm all code paths and every reusable source/output page.
    (void)measure(0u);
    (void)measure(1u);
    (void)measure(2u);
    (void)measure(3u);

    std::vector<double> baseline_samples;
    std::vector<double> fused_samples;
    std::vector<double> cached_samples;
    std::vector<double> tail_paired_samples;
    std::vector<double> paired_ratios;
    std::vector<double> cached_ratios;
    std::vector<double> production_ratios;
    std::vector<double> tail_paired_ratios;
    baseline_samples.reserve(samples);
    fused_samples.reserve(samples);
    cached_samples.reserve(samples);
    tail_paired_samples.reserve(samples);
    paired_ratios.reserve(samples);
    cached_ratios.reserve(samples);
    production_ratios.reserve(samples);
    tail_paired_ratios.reserve(samples);
    for (unsigned sample = 0; sample < samples; ++sample) {
        double timings[4] = {};
        if (sample % 4u == 0u) {
            timings[0] = measure(0u);
            timings[1] = measure(1u);
            timings[2] = measure(2u);
            timings[3] = measure(3u);
        }
        else if (sample % 4u == 1u) {
            timings[1] = measure(1u);
            timings[2] = measure(2u);
            timings[3] = measure(3u);
            timings[0] = measure(0u);
        }
        else if (sample % 4u == 2u) {
            timings[2] = measure(2u);
            timings[3] = measure(3u);
            timings[0] = measure(0u);
            timings[1] = measure(1u);
        }
        else {
            timings[3] = measure(3u);
            timings[0] = measure(0u);
            timings[1] = measure(1u);
            timings[2] = measure(2u);
        }
        baseline_samples.push_back(timings[0]);
        fused_samples.push_back(timings[1]);
        cached_samples.push_back(timings[2]);
        tail_paired_samples.push_back(timings[3]);
        paired_ratios.push_back(timings[1] / timings[0]);
        cached_ratios.push_back(timings[2] / timings[1]);
        production_ratios.push_back(timings[2] / timings[0]);
        tail_paired_ratios.push_back(timings[3] / timings[2]);
    }

    std::printf(
        "bb=%u K=%u mix_count=%u target_source_degree=%u "
        "rows=%zu reps=%u samples=%u "
        "mean_source_degree=%.3f "
        "baseline_median_ms=%.6f fused_median_ms=%.6f cached_median_ms=%.6f "
        "tail_paired_median_ms=%.6f "
        "paired_ratio_median=%.6f improvement=%.3f%% "
        "cached_vs_wrapper_ratio=%.6f cache_improvement=%.3f%% "
        "production_ratio_median=%.6f production_improvement=%.3f%% "
        "tail_paired_vs_cached_ratio=%.6f tail_paired_improvement=%.3f%% "
        "destination_traffic_reduction=%.3f%% "
        "total_traffic_reduction=%.3f%% sink=%" PRIu64 "\n",
        block_bytes, K, config.MixCount, target_source_degree,
        rows.size(), repetitions, samples,
        (double)degree_sum / rows.size(),
        Median(baseline_samples) * 1000.0,
        Median(fused_samples) * 1000.0,
        Median(cached_samples) * 1000.0,
        Median(tail_paired_samples) * 1000.0,
        Median(paired_ratios),
        (1.0 - Median(paired_ratios)) * 100.0,
        Median(cached_ratios),
        (1.0 - Median(cached_ratios)) * 100.0,
        Median(production_ratios),
        (1.0 - Median(production_ratios)) * 100.0,
        Median(tail_paired_ratios),
        (1.0 - Median(tail_paired_ratios)) * 100.0,
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

#if defined(WIREHAIR_V2_ENABLE_COUNT5_BENCHMARK)

enum class Count5Gate
{
    Active,
    Fallback,
    SameArm
};

struct Count5LogStats
{
    double GeometricMean = 0.0;
    double Lower95 = 0.0;
    double Upper95 = 0.0;
};

struct Count5PanelMetadata
{
    const char* Workload = nullptr;
    const char* Panel = nullptr;
    const char* Shape = nullptr;
    uint32_t K = 0u;
    uint32_t BlockBytes = 0u;
    uint32_t Units = 0u;
    uint32_t Repetitions = 0u;
    unsigned Samples = 0u;
    uint32_t EligibleUnits = 0u;
};

double Count5StudentT975(unsigned samples)
{
    if (samples == 12u) return 2.200985160;
    if (samples == 32u) return 2.039513446;
    return 0.0;
}

bool ComputeCount5LogStats(
    const std::vector<double>& log_ratios,
    Count5LogStats& stats_out)
{
    const double t_value = Count5StudentT975(
        (unsigned)log_ratios.size());
    if (log_ratios.size() < 2u || t_value <= 0.0) {
        return false;
    }
    long double sum = 0.0L;
    for (double value : log_ratios)
    {
        if (!std::isfinite(value)) return false;
        sum += (long double)value;
    }
    const long double mean = sum / log_ratios.size();
    long double squares = 0.0L;
    for (double value : log_ratios)
    {
        const long double difference = (long double)value - mean;
        squares += difference * difference;
    }
    const long double standard_error = std::sqrt(
        squares / (log_ratios.size() - 1u) / log_ratios.size());
    const long double half_width = t_value * standard_error;
    Count5LogStats stats;
    stats.GeometricMean = (double)std::exp(mean);
    stats.Lower95 = (double)std::exp(mean - half_width);
    stats.Upper95 = (double)std::exp(mean + half_width);
    if (!std::isfinite(stats.GeometricMean) ||
        !std::isfinite(stats.Lower95) ||
        !std::isfinite(stats.Upper95) || stats.Lower95 <= 0.0 ||
        stats.Lower95 > stats.GeometricMean ||
        stats.GeometricMean > stats.Upper95)
    {
        return false;
    }
    stats_out = stats;
    return true;
}

bool Count5GatePasses(Count5Gate gate, const Count5LogStats& stats)
{
    if (gate == Count5Gate::Active) {
        return stats.Upper95 < 0.99;
    }
    const bool bounded =
        stats.Lower95 >= 0.99 && stats.Upper95 <= 1.01;
    if (gate == Count5Gate::Fallback) {
        return bounded;
    }
    return bounded && stats.Lower95 <= 1.0 && stats.Upper95 >= 1.0;
}

template<class Measure>
bool RunCount5Panel(
    const Count5PanelMetadata& metadata,
    int left_mode,
    int right_mode,
    const char* left_name,
    const char* right_name,
    Count5Gate gate,
    Measure& measure)
{
    uint64_t ignored = 0u;
    if (!measure(left_mode, ignored) || !measure(right_mode, ignored)) {
        std::fprintf(stderr,
            "count5 warmup failed workload=%s panel=%s\n",
            metadata.Workload, metadata.Panel);
        return false;
    }

    std::vector<double> log_ratios;
    std::vector<uint64_t> left_timings;
    std::vector<uint64_t> right_timings;
    log_ratios.reserve(metadata.Samples);
    left_timings.reserve(metadata.Samples);
    right_timings.reserve(metadata.Samples);
    for (unsigned sample = 0u; sample < metadata.Samples; ++sample)
    {
        uint64_t left_ns = 0u;
        uint64_t right_ns = 0u;
        const bool left_first = (sample & 1u) == 0u;
        const bool measured = left_first ?
            (measure(left_mode, left_ns) &&
                measure(right_mode, right_ns)) :
            (measure(right_mode, right_ns) &&
                measure(left_mode, left_ns));
        if (!measured || left_ns == 0u || right_ns == 0u) {
            std::fprintf(stderr,
                "count5 measurement failed workload=%s panel=%s sample=%u\n",
                metadata.Workload, metadata.Panel, sample);
            return false;
        }
        const double ratio = (double)right_ns / (double)left_ns;
        const double log_ratio = std::log(ratio);
        if (!std::isfinite(log_ratio)) return false;
        log_ratios.push_back(log_ratio);
        left_timings.push_back(left_ns);
        right_timings.push_back(right_ns);
    }

    // Emit only after the complete panel so terminal or pipe I/O cannot
    // perturb the subsequent paired sample.
    for (unsigned sample = 0u; sample < metadata.Samples; ++sample)
    {
        const double ratio =
            (double)right_timings[sample] / (double)left_timings[sample];
        std::printf(
            "count5_raw,version=1,workload=%s,panel=%s,shape=%s,"
            "K=%u,block_bytes=%u,units=%u,repetitions=%u,samples=%u,"
            "eligible_units=%u,sample=%u,first=%s,left_arm=%s,"
            "right_arm=%s,left_ns=%" PRIu64 ",right_ns=%" PRIu64
            ",right_over_left=%.12f\n",
            metadata.Workload, metadata.Panel, metadata.Shape,
            metadata.K, metadata.BlockBytes, metadata.Units,
            metadata.Repetitions, metadata.Samples, metadata.EligibleUnits,
            sample, (sample & 1u) == 0u ? "left" : "right",
            left_name, right_name,
            left_timings[sample], right_timings[sample], ratio);
    }

    Count5LogStats stats;
    if (!ComputeCount5LogStats(log_ratios, stats)) {
        return false;
    }
    const bool gate_pass = Count5GatePasses(gate, stats);
    std::printf(
        "count5_summary,version=1,workload=%s,panel=%s,shape=%s,"
        "K=%u,block_bytes=%u,units=%u,repetitions=%u,samples=%u,"
        "eligible_units=%u,left_arm=%s,right_arm=%s,"
        "geometric_mean=%.12f,lower95=%.12f,upper95=%.12f,"
        "student_t_975=%.9f,gate=%s\n",
        metadata.Workload, metadata.Panel, metadata.Shape,
        metadata.K, metadata.BlockBytes, metadata.Units,
        metadata.Repetitions, metadata.Samples, metadata.EligibleUnits,
        left_name, right_name, stats.GeometricMean,
        stats.Lower95, stats.Upper95,
        Count5StudentT975(metadata.Samples), gate_pass ? "pass" : "fail");
    return gate_pass;
}

struct Count5ModeReset
{
    ~Count5ModeReset()
    {
        (void)wirehair_v2::SetPacketMix3FiveTermSetXorModeForBenchmark(0);
    }
};

bool BuildCount5System(uint32_t K, wirehair_v2::PrecodeSystem& system)
{
    return wirehair_v2::BuildPrecodeSystem(
            wirehair_v2::MakeCertifiedParams(
                K, UINT64_C(0x786f72667573696f)),
            system) &&
        wirehair_v2::ValidatePrecodeSystem(system);
}

std::vector<uint32_t> Count5RowsWithDegree(
    const wirehair_v2::PrecodeSystem& system,
    const wirehair_v2::PacketRowConfig& config,
    uint32_t target_degree,
    size_t count)
{
    const uint32_t K = system.Params.BlockCount;
    const uint32_t P = system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    std::vector<uint32_t> rows;
    rows.reserve(count);
    for (uint32_t id = 0u; id < UINT32_C(1000000) && rows.size() < count;
         ++id)
    {
        wirehair::PeelRowParameters params;
        params.Initialize(
            id, config.PeelSeed, (uint16_t)K, (uint16_t)P);
        if (params.PeelCount == target_degree) {
            rows.push_back(id);
        }
    }
    return rows;
}

uint32_t Count5EligibleSourceIds(
    const wirehair_v2::PrecodeSystem& system,
    const wirehair_v2::PacketRowConfig& config)
{
    const uint32_t K = system.Params.BlockCount;
    const uint32_t P = system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    uint32_t eligible = 0u;
    for (uint32_t id = 0u; id < K; ++id)
    {
        wirehair::PeelRowParameters params;
        params.Initialize(
            id, config.PeelSeed, (uint16_t)K, (uint16_t)P);
        eligible += params.PeelCount == 2u ? 1u : 0u;
    }
    return eligible;
}

unsigned Count5PacketRepetitions(uint32_t block_bytes)
{
    if (block_bytes == 1280u) return 256u;
    if (block_bytes == 4096u) return 80u;
    if (block_bytes == 32768u) return 10u;
    return 0u;
}

bool RunCount5PacketScreen(
    uint32_t K,
    uint32_t block_bytes,
    unsigned samples)
{
    Count5ModeReset mode_reset;
    wirehair_v2::PrecodeSystem system;
    if (!BuildCount5System(K, system)) return false;
    wirehair_v2::PacketRowConfig config;
    config.PeelSeed = UINT32_C(0x4d241359);
    config.MixCount = 3u;
    const uint32_t P = system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    wirehair_v2::PacketRowRuntime runtime;
    if (!runtime.Initialize(K, P, config.MixCount)) return false;

    const std::vector<uint32_t> active_rows =
        Count5RowsWithDegree(system, config, 2u, 256u);
    const std::vector<uint32_t> fallback_rows =
        Count5RowsWithDegree(system, config, 3u, 256u);
    if (active_rows.size() != 256u || fallback_rows.size() != 256u) {
        std::fprintf(stderr,
            "count5 row fixture shortfall K=%u active=%zu fallback=%zu\n",
            K, active_rows.size(), fallback_rows.size());
        return false;
    }
    const uint64_t intermediate_bytes_wide =
        ((uint64_t)K + P) * block_bytes;
    if (intermediate_bytes_wide >
        (uint64_t)std::numeric_limits<size_t>::max())
    {
        return false;
    }
    AlignedStorage intermediate((size_t)intermediate_bytes_wide);
    AlignedStorage left_output(block_bytes);
    AlignedStorage right_output(block_bytes);
    for (size_t i = 0u; i < (size_t)intermediate_bytes_wide; ++i) {
        intermediate.Data[i] = (uint8_t)(i * 131u + (i >> 13) + 17u);
    }

    const auto verify_rows = [&](
        const std::vector<uint32_t>& rows, uint64_t expected_ops) {
        for (uint32_t id : rows)
        {
            uint64_t left_ops = 0u;
            uint64_t right_ops = 0u;
            if (!wirehair_v2::SetPacketMix3FiveTermSetXorModeForBenchmark(-1) ||
                !EvaluateCached(
                    system, config, runtime, intermediate.Data, block_bytes,
                    id, left_output.Data, &left_ops) ||
                !wirehair_v2::SetPacketMix3FiveTermSetXorModeForBenchmark(+1) ||
                !EvaluateCached(
                    system, config, runtime, intermediate.Data, block_bytes,
                    id, right_output.Data, &right_ops) ||
                std::memcmp(
                    left_output.Data, right_output.Data, block_bytes) != 0 ||
                left_ops != expected_ops || right_ops != expected_ops)
            {
                std::fprintf(stderr,
                    "count5 packet oracle failed K=%u bb=%u id=%u "
                    "left_ops=%" PRIu64 " right_ops=%" PRIu64
                    " expected_ops=%" PRIu64 "\n",
                    K, block_bytes, id,
                    left_ops, right_ops, expected_ops);
                return false;
            }
        }
        return true;
    };
    if (!verify_rows(active_rows, 5u) ||
        !verify_rows(fallback_rows, 6u))
    {
        return false;
    }

    const unsigned repetitions = Count5PacketRepetitions(block_bytes);
    const std::vector<uint32_t>* timed_rows = &active_rows;
    volatile uint64_t sink = 0u;
    const auto measure = [&](int mode, uint64_t& elapsed_ns) {
        if (!wirehair_v2::SetPacketMix3FiveTermSetXorModeForBenchmark(mode)) {
            return false;
        }
        const Clock::time_point start = Clock::now();
        for (unsigned repetition = 0u; repetition < repetitions; ++repetition)
        {
            for (uint32_t id : *timed_rows)
            {
                if (!EvaluateCached(
                        system, config, runtime, intermediate.Data,
                        block_bytes, id, left_output.Data))
                {
                    return false;
                }
            }
        }
        elapsed_ns = (uint64_t)
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                Clock::now() - start).count();
        sink ^= left_output.Data[
            ((size_t)(mode + 1) + timed_rows->size() + repetitions) %
                block_bytes];
        return elapsed_ns != 0u;
    };

    Count5PanelMetadata metadata;
    metadata.Workload = "packet";
    metadata.K = K;
    metadata.BlockBytes = block_bytes;
    metadata.Units = 256u;
    metadata.Repetitions = repetitions;
    metadata.Samples = samples;
    metadata.EligibleUnits = 256u;
    metadata.Panel = "active";
    metadata.Shape = "degree2_mix3_terms5";
    bool ok = RunCount5Panel(
        metadata, -1, +1, "legacy", "candidate", Count5Gate::Active,
        measure);

    timed_rows = &fallback_rows;
    metadata.Panel = "fallback";
    metadata.Shape = "degree3_mix3_terms6";
    metadata.EligibleUnits = 0u;
    ok = RunCount5Panel(
        metadata, -1, +1, "legacy_mode", "candidate_mode",
        Count5Gate::Fallback, measure) && ok;

    timed_rows = &active_rows;
    metadata.Panel = "aa";
    metadata.Shape = "degree2_mix3_terms5";
    metadata.EligibleUnits = 256u;
    ok = RunCount5Panel(
        metadata, -1, -1, "legacy_a", "legacy_b", Count5Gate::SameArm,
        measure) && ok;
    std::printf(
        "count5_screen,version=1,workload=packet,K=%u,block_bytes=%u,"
        "sink=%" PRIu64 ",gate=%s\n",
        K, block_bytes, (uint64_t)sink, ok ? "pass" : "fail");
    return ok;
}

class Count5RecoverFixture
{
public:
    bool Initialize(uint32_t K, uint32_t block_bytes)
    {
        KValue = K;
        BlockBytesValue = block_bytes;
        if (!BuildCount5System(K, System)) return false;
        Config.PeelSeed = UINT32_C(0x4d241359);
        Config.MixCount = 3u;
        const uint64_t message_bytes_wide = (uint64_t)K * block_bytes;
        if (message_bytes_wide >
            (uint64_t)std::numeric_limits<size_t>::max())
        {
            return false;
        }
        const size_t message_bytes = (size_t)message_bytes_wide;
        Decoder.reset(new wirehair_v2::MessagePrecodeDecoder);
        if (Decoder->InitializeForValidatedSystemForBenchmark(
                message_bytes_wide, block_bytes, System, Config, 0u) !=
            Wirehair_Success)
        {
            return false;
        }
        std::vector<uint8_t> zero_block(block_bytes, 0u);
        WirehairResult result = Wirehair_NeedMore;
        for (uint32_t packet_index = 0u;
             packet_index < K + 256u && result == Wirehair_NeedMore;
             ++packet_index)
        {
            result = Decoder->DecodeResult(
                K + packet_index, zero_block.data(), block_bytes);
            Delivered = packet_index + 1u;
        }
        if (result != Wirehair_Success || !Decoder->IsDecoded()) {
            std::fprintf(stderr,
                "count5 repair-only fixture failed K=%u bb=%u delivered=%u "
                "result=%d\n",
                K, block_bytes, Delivered, (int)result);
            return false;
        }
        Output.assign(message_bytes, uint8_t{0xa5});
        if (!wirehair_v2::SetPacketMix3FiveTermSetXorModeForBenchmark(-1) ||
            Decoder->RecoverResult(Output.data(), Output.size()) !=
                Wirehair_Success ||
            std::find_if(
                Output.begin(), Output.end(),
                [](uint8_t value) { return value != 0u; }) != Output.end())
        {
            std::fprintf(stderr,
                "count5 recovery warm/oracle failed K=%u bb=%u\n",
                K, block_bytes);
            return false;
        }
        // The first legacy recovery may release retained cold receive state.
        // Check the candidate on the stable repeated-RecoverResult path that
        // the timed panels use, and overwrite a nonzero destination so a
        // missing write cannot pass the all-zero message oracle.
        std::fill(Output.begin(), Output.end(), uint8_t{0xa5});
        if (!wirehair_v2::SetPacketMix3FiveTermSetXorModeForBenchmark(+1) ||
            Decoder->RecoverResult(Output.data(), Output.size()) !=
                Wirehair_Success ||
            std::find_if(
                Output.begin(), Output.end(),
                [](uint8_t value) { return value != 0u; }) != Output.end())
        {
            std::fprintf(stderr,
                "count5 candidate recovery oracle failed K=%u bb=%u\n",
                K, block_bytes);
            return false;
        }
        EligibleSourceIds = Count5EligibleSourceIds(System, Config);
        return true;
    }

    bool Measure(int mode, uint64_t& elapsed_ns)
    {
        if (!Decoder ||
            !wirehair_v2::SetPacketMix3FiveTermSetXorModeForBenchmark(mode))
        {
            return false;
        }
        const Clock::time_point start = Clock::now();
        const WirehairResult result =
            Decoder->RecoverResult(Output.data(), Output.size());
        elapsed_ns = (uint64_t)
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                Clock::now() - start).count();
        if (result != Wirehair_Success || elapsed_ns == 0u) return false;
        Sink ^= Output[
            ((size_t)(mode + 1) + Delivered + EligibleSourceIds) %
                Output.size()];
        return true;
    }

    uint32_t K() const { return KValue; }
    uint32_t BlockBytes() const { return BlockBytesValue; }
    uint32_t Eligible() const { return EligibleSourceIds; }
    uint32_t DeliveredPackets() const { return Delivered; }
    uint64_t SinkValue() const { return (uint64_t)Sink; }

private:
    wirehair_v2::PrecodeSystem System;
    wirehair_v2::PacketRowConfig Config;
    std::unique_ptr<wirehair_v2::MessagePrecodeDecoder> Decoder;
    std::vector<uint8_t> Output;
    uint32_t KValue = 0u;
    uint32_t BlockBytesValue = 0u;
    uint32_t EligibleSourceIds = 0u;
    uint32_t Delivered = 0u;
    volatile uint64_t Sink = 0u;
};

bool RunCount5RecoverScreen(
    uint32_t K,
    uint32_t block_bytes,
    unsigned samples)
{
    Count5ModeReset mode_reset;
    bool ok = true;
    {
        Count5RecoverFixture active;
        if (!active.Initialize(K, block_bytes)) return false;
        const auto measure = [&](int mode, uint64_t& elapsed_ns) {
            return active.Measure(mode, elapsed_ns);
        };
        Count5PanelMetadata metadata;
        metadata.Workload = "recover";
        metadata.Panel = "active";
        metadata.Shape = "all_source_ids_mix3";
        metadata.K = K;
        metadata.BlockBytes = block_bytes;
        metadata.Units = K;
        metadata.Repetitions = 1u;
        metadata.Samples = samples;
        metadata.EligibleUnits = active.Eligible();
        if (metadata.EligibleUnits == 0u) return false;
        ok = RunCount5Panel(
            metadata, -1, +1, "legacy", "candidate",
            Count5Gate::Active, measure) && ok;

        metadata.Panel = "aa";
        ok = RunCount5Panel(
            metadata, -1, -1, "legacy_a", "legacy_b",
            Count5Gate::SameArm, measure) && ok;
        std::printf(
            "count5_fixture,version=1,workload=recover,panel=active,"
            "K=%u,block_bytes=%u,delivered=%u,eligible_units=%u,"
            "sink=%" PRIu64 "\n",
            active.K(), active.BlockBytes(), active.DeliveredPackets(),
            active.Eligible(), active.SinkValue());
    }

    // The exact lower-K boundary is the recovery fallback: it exercises the
    // same MIX3 evaluator and block width while making every five-term row
    // ineligible for the new K >= 10000 dispatch.
    if (K == 10000u)
    {
        Count5RecoverFixture fallback;
        if (!fallback.Initialize(9999u, block_bytes)) return false;
        const auto measure = [&](int mode, uint64_t& elapsed_ns) {
            return fallback.Measure(mode, elapsed_ns);
        };
        Count5PanelMetadata metadata;
        metadata.Workload = "recover";
        metadata.Panel = "fallback";
        metadata.Shape = "all_source_ids_mix3_K9999";
        metadata.K = fallback.K();
        metadata.BlockBytes = block_bytes;
        metadata.Units = fallback.K();
        metadata.Repetitions = 1u;
        metadata.Samples = samples;
        metadata.EligibleUnits = 0u;
        ok = RunCount5Panel(
            metadata, -1, +1, "legacy_mode", "candidate_mode",
            Count5Gate::Fallback, measure) && ok;
        std::printf(
            "count5_fixture,version=1,workload=recover,panel=fallback,"
            "K=%u,block_bytes=%u,delivered=%u,eligible_units=0,"
            "sink=%" PRIu64 "\n",
            fallback.K(), fallback.BlockBytes(),
            fallback.DeliveredPackets(), fallback.SinkValue());
    }
    std::printf(
        "count5_screen,version=1,workload=recover,K=%u,block_bytes=%u,"
        "gate=%s\n",
        K, block_bytes, ok ? "pass" : "fail");
    return ok;
}

#endif // WIREHAIR_V2_ENABLE_COUNT5_BENCHMARK

} // namespace

int main(int argc, char** argv)
{
    unsigned samples = 40u;
    uint32_t mix_count = wirehair_v2::kCertifiedPacketMixCount;
    uint32_t K = 128u;
    bool run_packet = true;
    bool run_gather = false;
    bool mtu_only = false;
    bool crossover_only = false;
    bool boundary_only = false;
    bool fallback_only = false;
    bool degrees_only = false;
    bool count5_packet = false;
    bool count5_recover = false;
    if (argc >= 2) {
        char* end = nullptr;
        const unsigned long parsed = std::strtoul(argv[1], &end, 10);
        if (!end || *end != '\0' || parsed < 12u || parsed > 1000u) {
            std::fprintf(stderr,
                "usage: %s [samples=30..1000] [mix_count=1|2|3] "
                "[mode=packet|mtu|crossover|boundary|fallback|degrees|"
                "gather|all] [K=2..64000]\n"
                "       %s 32 3 count5-packet "
                "<K=10000|64000> <block_bytes=1280|4096|32768>\n"
                "       %s 12 3 count5-recover "
                "<K=10000|64000> <block_bytes=1280|4096|32768>\n",
                argv[0], argv[0],
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
                "[mode=packet|mtu|crossover|boundary|fallback|degrees|"
                "gather|all] "
                "[K=2..64000]\n",
                argv[0]);
            return 2;
        }
        mix_count = (uint32_t)parsed;
    }
    if (argc >= 4)
    {
        if (std::strcmp(argv[3], "packet") == 0) {
            run_packet = true;
            run_gather = false;
        }
        else if (std::strcmp(argv[3], "mtu") == 0) {
            run_packet = true;
            run_gather = false;
            mtu_only = true;
        }
        else if (std::strcmp(argv[3], "crossover") == 0) {
            run_packet = true;
            run_gather = false;
            crossover_only = true;
        }
        else if (std::strcmp(argv[3], "boundary") == 0) {
            run_packet = true;
            run_gather = false;
            boundary_only = true;
        }
        else if (std::strcmp(argv[3], "fallback") == 0) {
            run_packet = true;
            run_gather = false;
            fallback_only = true;
        }
        else if (std::strcmp(argv[3], "degrees") == 0) {
            run_packet = true;
            run_gather = false;
            degrees_only = true;
        }
        else if (std::strcmp(argv[3], "count5-packet") == 0) {
            count5_packet = true;
        }
        else if (std::strcmp(argv[3], "count5-recover") == 0) {
            count5_recover = true;
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
                "[mode=packet|mtu|crossover|boundary|fallback|degrees|"
                "gather|all|count5-packet|count5-recover] "
                "[K=2..64000]\n",
                argv[0]);
            return 2;
        }
    }
    if (count5_packet || count5_recover)
    {
        if (argc != 6 || mix_count != 3u ||
            samples != (count5_packet ? 32u : 12u))
        {
            std::fprintf(stderr,
                "usage: %s %u 3 %s <K=10000|64000> "
                "<block_bytes=1280|4096|32768>\n",
                argv[0], count5_packet ? 32u : 12u,
                count5_packet ? "count5-packet" : "count5-recover");
            return 2;
        }
        char* K_end = nullptr;
        char* bytes_end = nullptr;
        const unsigned long parsed_K = std::strtoul(argv[4], &K_end, 10);
        const unsigned long parsed_bytes =
            std::strtoul(argv[5], &bytes_end, 10);
        if (!K_end || *K_end != '\0' ||
            (parsed_K != 10000u && parsed_K != 64000u) ||
            !bytes_end || *bytes_end != '\0' ||
            (parsed_bytes != 1280u && parsed_bytes != 4096u &&
                parsed_bytes != 32768u))
        {
            std::fprintf(stderr,
                "count5 screen is frozen to K={10000,64000}, "
                "block_bytes={1280,4096,32768}\n");
            return 2;
        }
        if (gf256_init() != 0) return 1;
#if defined(WIREHAIR_V2_ENABLE_COUNT5_BENCHMARK)
        try
        {
            const bool ok = count5_packet ?
                RunCount5PacketScreen(
                    (uint32_t)parsed_K, (uint32_t)parsed_bytes, samples) :
                RunCount5RecoverScreen(
                    (uint32_t)parsed_K, (uint32_t)parsed_bytes, samples);
            return ok ? 0 : 1;
        }
        catch (const std::bad_alloc&)
        {
            std::fprintf(stderr, "count5 screen allocation failed\n");
            return 1;
        }
        catch (const std::length_error&)
        {
            std::fprintf(stderr, "count5 screen size overflow\n");
            return 1;
        }
#else
        std::fprintf(stderr,
            "count5 screens require BUILD_TESTS=ON timing policy\n");
        return 2;
#endif
    }
    if (samples < 30u) {
        std::fprintf(stderr, "legacy modes require at least 30 samples\n");
        return 2;
    }
    if (argc == 5)
    {
        char* end = nullptr;
        const unsigned long parsed = std::strtoul(argv[4], &end, 10);
        if (!end || *end != '\0' || parsed < 2u || parsed > 64000u) {
            std::fprintf(stderr,
                "usage: %s [samples=30..1000] [mix_count=1|2|3] "
                "[mode=packet|mtu|crossover|boundary|fallback|degrees|"
                "gather|all] "
                "[K=2..64000]\n",
                argv[0]);
            return 2;
        }
        K = (uint32_t)parsed;
    }
    else if (argc > 5) {
        std::fprintf(stderr,
            "usage: %s [samples=30..1000] [mix_count=1|2|3] "
            "[mode=packet|mtu|crossover|boundary|fallback|degrees|"
            "gather|all] "
            "[K=2..64000]\n",
            argv[0]);
        return 2;
    }
    if (K != 128u && !mtu_only && !degrees_only)
    {
        std::fprintf(stderr,
            "non-default K is supported only by bounded MTU/degree screens\n");
        return 2;
    }
    if (gf256_init() != 0) {
        return 1;
    }

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
    if (run_packet && !mtu_only && !crossover_only && !boundary_only &&
        !fallback_only && !degrees_only)
    {
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
        if (mtu_only) {
            ok = RunSize(
                system, config, 1280u, 256u, samples, 0u, true) && ok;
            return ok ? 0 : 1;
        }
        if (crossover_only)
        {
            ok = RunSize(system, config, 4u * 1024u, 128u, samples) && ok;
            ok = RunSize(system, config, 16u * 1024u, 32u, samples) && ok;
            ok = RunSize(system, config, 32u * 1024u, 16u, samples) && ok;
            ok = RunSize(system, config, 64u * 1024u, 8u, samples) && ok;
            ok = RunSize(system, config, 100u * 1024u, 5u, samples) && ok;
            return ok ? 0 : 1;
        }
        if (boundary_only)
        {
            ok = RunSize(
                system, config, 32767u, 8u, samples, 0u, true) && ok;
            ok = RunSize(
                system, config, 32768u, 8u, samples, 0u, true) && ok;
            ok = RunSize(
                system, config, 32769u, 8u, samples, 0u, true) && ok;
            return ok ? 0 : 1;
        }
        if (fallback_only)
        {
            ok = RunSize(
                system, config, 32769u, 8u, samples, 0u, true, 64u) && ok;
            ok = RunSize(
                system, config, 100u * 1024u, 8u, samples, 0u, true, 32u) && ok;
            ok = RunSize(
                system, config, 1024u * 1024u, 1u, samples, 0u, true, 12u) && ok;
            return ok ? 0 : 1;
        }
        if (degrees_only)
        {
            for (uint32_t degree = 2u; degree <= 6u; ++degree) {
                ok = RunSize(
                    system, config, 1280u, 256u, samples, degree, true) && ok;
            }
            return ok ? 0 : 1;
        }
        ok = RunSize(system, config, 1u, 4096u, samples) && ok;
        // Keep a short sub-cacheline sweep so wrapper overhead at bb=1 cannot
        // conceal a regression (or win) in small production packets.
        ok = RunSize(system, config, 2u, 4096u, samples) && ok;
        ok = RunSize(system, config, 8u, 4096u, samples) && ok;
        ok = RunSize(system, config, 32u, 4096u, samples) && ok;
        ok = RunSize(system, config, 1280u, 2048u, samples) && ok;
        ok = RunSize(system, config, 100u * 1024u, 32u, samples) && ok;
        ok = RunSize(system, config, 1024u * 1024u, 3u, samples) && ok;
    }
    return ok ? 0 : 1;
}
