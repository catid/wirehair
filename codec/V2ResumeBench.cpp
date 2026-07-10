#include "WirehairV2PrecodeEncode.h"
#include "WirehairV2Solve.h"

#include <chrono>
#include <cinttypes>
#include <cstdio>
#include <vector>

#if defined(__linux__)
#include <sched.h>
#endif

namespace {

static const uint32_t kOverheadPackets = 8u;
static const unsigned kSamples = 20u;

typedef std::chrono::steady_clock BenchClock;

uint64_t ElapsedNanoseconds(
    const BenchClock::time_point& start,
    const BenchClock::time_point& end)
{
    return (uint64_t)std::chrono::duration_cast<std::chrono::nanoseconds>(
        end - start).count();
}

int PinToFirstAvailableCpu()
{
#if defined(__linux__)
    cpu_set_t available;
    CPU_ZERO(&available);
    if (sched_getaffinity(0, sizeof(available), &available) != 0) {
        return -1;
    }
    for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu)
    {
        if (!CPU_ISSET(cpu, &available)) {
            continue;
        }
        cpu_set_t selected;
        CPU_ZERO(&selected);
        CPU_SET(cpu, &selected);
        if (sched_setaffinity(0, sizeof(selected), &selected) == 0) {
            return cpu;
        }
        return -1;
    }
#endif
    return -1;
}

struct Fixture
{
    uint32_t K = 0u;
    uint32_t BlockBytes = 0u;
    wirehair_v2::PrecodeSystem System;
    wirehair_v2::PacketRowConfig Config;
    std::vector<uint8_t> PacketData;
    std::vector<wirehair_v2::SolvePacket> FullPackets;
    std::vector<wirehair_v2::SolvePacket> InitialPackets;
    std::vector<std::vector<wirehair_v2::SolvePacket> > ColdPrefixes;
    std::vector<uint8_t> ExpectedIntermediate;

    bool Initialize(uint32_t source_count, uint32_t block_bytes)
    {
        K = source_count;
        BlockBytes = block_bytes;
        wirehair_v2::PrecodeParams params =
            wirehair_v2::MakeCertifiedParams(
                K,
                UINT64_C(0x726573756d656265) ^ K);
        wirehair_v2::PacketRowConfig base_config;
        base_config.PeelSeed = UINT32_C(0x9e3779b9) ^ K;
        base_config.MixCount = wirehair_v2::kCertifiedPacketMixCount;
        if (wirehair_v2::SelectSystematicConfiguration(
                params,
                base_config,
                System,
                Config) != Wirehair_Success)
        {
            return false;
        }

        PacketData.resize((size_t)K * BlockBytes);
        for (size_t i = 0; i < PacketData.size(); ++i) {
            PacketData[i] = (uint8_t)(
                i * 197u + (i >> 9) + K * 11u + BlockBytes);
        }
        FullPackets.resize(K);
        for (uint32_t id = 0u; id < K; ++id)
        {
            FullPackets[id].BlockId = id;
            FullPackets[id].Data =
                PacketData.data() + (size_t)id * BlockBytes;
        }

        InitialPackets.reserve(K);
        for (uint32_t id = 0u; id + kOverheadPackets < K; ++id) {
            InitialPackets.push_back(FullPackets[id]);
        }
        for (uint32_t i = 0u; i < kOverheadPackets; ++i) {
            InitialPackets.push_back(FullPackets[0]);
        }
        if (InitialPackets.size() != K) {
            return false;
        }

        ColdPrefixes.resize(kOverheadPackets);
        for (uint32_t i = 0u; i < kOverheadPackets; ++i)
        {
            ColdPrefixes[i] = InitialPackets;
            for (uint32_t extra = 0u; extra <= i; ++extra) {
                ColdPrefixes[i].push_back(
                    FullPackets[K - kOverheadPackets + extra]);
            }
        }

        return wirehair_v2::SolvePrecodeSystem(
                   System,
                   Config,
                   FullPackets,
                   BlockBytes,
                   ExpectedIntermediate) == Wirehair_Success;
    }
};

struct Sample
{
    uint64_t InitialNanoseconds = 0u;
    uint64_t OverheadNanoseconds = 0u;
    size_t CheckpointBytes = 0u;
};

bool RunWarm(const Fixture& fixture, Sample& sample)
{
    std::vector<uint8_t> output;
    wirehair_v2::PrecodeSolveResumeState checkpoint;
    std::vector<uint8_t> pending_packet;
    const BenchClock::time_point initial_start = BenchClock::now();
    const WirehairResult initial_result = wirehair_v2::SolvePrecodeSystem(
        fixture.System,
        fixture.Config,
        fixture.InitialPackets,
        fixture.BlockBytes,
        output,
        nullptr,
        &checkpoint);
    if (initial_result == Wirehair_NeedMore && checkpoint.Active) {
        pending_packet.resize(fixture.BlockBytes);
    }
    const BenchClock::time_point initial_end = BenchClock::now();
    if (initial_result != Wirehair_NeedMore || !checkpoint.Active ||
        !output.empty())
    {
        return false;
    }
    sample.InitialNanoseconds =
        ElapsedNanoseconds(initial_start, initial_end);
    sample.CheckpointBytes =
        checkpoint.PersistentBytes() + pending_packet.capacity();

    const BenchClock::time_point overhead_start = BenchClock::now();
    for (uint32_t i = 0u; i < kOverheadPackets; ++i)
    {
        const uint32_t id = fixture.K - kOverheadPackets + i;
        const WirehairResult result = wirehair_v2::ResumePrecodeSystem(
            fixture.System,
            fixture.Config,
            id,
            fixture.FullPackets[id].Data,
            fixture.BlockBytes,
            checkpoint,
            output);
        const WirehairResult expected =
            i + 1u == kOverheadPackets ?
                Wirehair_Success : Wirehair_NeedMore;
        if (result != expected) {
            return false;
        }
    }
    const BenchClock::time_point overhead_end = BenchClock::now();
    sample.OverheadNanoseconds =
        ElapsedNanoseconds(overhead_start, overhead_end);
    return !checkpoint.Active && output == fixture.ExpectedIntermediate;
}

bool RunCold(const Fixture& fixture, Sample& sample)
{
    std::vector<uint8_t> initial_output;
    const BenchClock::time_point initial_start = BenchClock::now();
    const WirehairResult initial_result = wirehair_v2::SolvePrecodeSystem(
        fixture.System,
        fixture.Config,
        fixture.InitialPackets,
        fixture.BlockBytes,
        initial_output);
    const BenchClock::time_point initial_end = BenchClock::now();
    if (initial_result != Wirehair_NeedMore || !initial_output.empty()) {
        return false;
    }
    sample.InitialNanoseconds =
        ElapsedNanoseconds(initial_start, initial_end);

    std::vector<uint8_t> output;
    const BenchClock::time_point overhead_start = BenchClock::now();
    for (uint32_t i = 0u; i < kOverheadPackets; ++i)
    {
        output.clear();
        const WirehairResult result = wirehair_v2::SolvePrecodeSystem(
            fixture.System,
            fixture.Config,
            fixture.ColdPrefixes[i],
            fixture.BlockBytes,
            output);
        const WirehairResult expected =
            i + 1u == kOverheadPackets ?
                Wirehair_Success : Wirehair_NeedMore;
        if (result != expected) {
            return false;
        }
    }
    const BenchClock::time_point overhead_end = BenchClock::now();
    sample.OverheadNanoseconds =
        ElapsedNanoseconds(overhead_start, overhead_end);
    return output == fixture.ExpectedIntermediate;
}

bool BenchmarkCase(uint32_t K, uint32_t block_bytes)
{
    Fixture fixture;
    if (!fixture.Initialize(K, block_bytes)) {
        std::fprintf(stderr, "resume benchmark fixture failed for K=%u\n", K);
        return false;
    }

    // Warm instruction/data caches and allocator paths before collecting the
    // pinned, interleaved samples used by the acceptance ratios.
    for (unsigned warmup = 0u; warmup < 2u; ++warmup)
    {
        Sample ignored;
        if (!RunWarm(fixture, ignored) || !RunCold(fixture, ignored)) {
            return false;
        }
    }

    uint64_t warm_initial_sum = 0u;
    uint64_t cold_initial_sum = 0u;
    uint64_t warm_overhead_sum = 0u;
    uint64_t cold_overhead_sum = 0u;
    size_t checkpoint_bytes = 0u;
    for (unsigned sample_i = 0u; sample_i < kSamples; ++sample_i)
    {
        Sample warm;
        Sample cold;
        const bool ok = (sample_i & 1u) == 0u ?
            (RunCold(fixture, cold) && RunWarm(fixture, warm)) :
            (RunWarm(fixture, warm) && RunCold(fixture, cold));
        if (!ok) {
            std::fprintf(stderr,
                "resume benchmark sample failed K=%u sample=%u\n",
                K,
                sample_i);
            return false;
        }
        warm_initial_sum += warm.InitialNanoseconds;
        cold_initial_sum += cold.InitialNanoseconds;
        warm_overhead_sum += warm.OverheadNanoseconds;
        cold_overhead_sum += cold.OverheadNanoseconds;
        checkpoint_bytes = warm.CheckpointBytes;
    }

    const size_t receive_bytes =
        (size_t)K * block_bytes + (size_t)K * sizeof(uint32_t);
    const double initial_ratio =
        (double)warm_initial_sum / (double)cold_initial_sum;
    const double overhead_ratio =
        (double)warm_overhead_sum / (double)cold_overhead_sum;
    const double saved_percent = (1.0 - overhead_ratio) * 100.0;
    const double memory_growth =
        ((double)checkpoint_bytes / (double)receive_bytes - 1.0) * 100.0;

    std::printf(
        "K=%u bb=%u samples=%u overhead=8 "
        "initial_cold_ms=%.3f initial_warm_ms=%.3f initial_ratio=%.4f "
        "overhead_cold_ms=%.3f overhead_warm_ms=%.3f saved=%.2f%% "
        "receive_bytes=%zu checkpoint_bytes=%zu memory_delta=%.2f%%\n",
        K,
        block_bytes,
        kSamples,
        cold_initial_sum / 1000000.0,
        warm_initial_sum / 1000000.0,
        initial_ratio,
        cold_overhead_sum / 1000000.0,
        warm_overhead_sum / 1000000.0,
        saved_percent,
        receive_bytes,
        checkpoint_bytes,
        memory_growth);

    const size_t memory_limit = receive_bytes + receive_bytes / 4u;
    if (initial_ratio > 1.05 || overhead_ratio > 0.50 ||
        checkpoint_bytes > memory_limit)
    {
        std::fprintf(stderr,
            "resume benchmark acceptance failed K=%u: initial=%.4f "
            "overhead=%.4f memory=%zu/%zu\n",
            K,
            initial_ratio,
            overhead_ratio,
            checkpoint_bytes,
            memory_limit);
        return false;
    }
    return true;
}

} // namespace

int main()
{
    const int cpu = PinToFirstAvailableCpu();
    if (cpu >= 0) {
        std::printf("resume benchmark pinned_cpu=%d\n", cpu);
    }
    else {
        std::printf("resume benchmark pinned_cpu=unavailable\n");
    }
    const bool k1000 = BenchmarkCase(1000u, 1280u);
    const bool k10000 = BenchmarkCase(10000u, 1280u);
    if (!k1000 || !k10000) {
        return 1;
    }
    std::printf("resume benchmark acceptance: PASS\n");
    return 0;
}
