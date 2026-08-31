// Bounded Stage-0 falsifier for wirehair-sxvz.16.1.20.42.
//
// This executable measures only the XOR schedules over 256 packet rows whose
// exact production equation columns and source pointers are prepared before
// timing.  It does not time packet-row generation, pointer gathering, codec
// recovery, allocation, counters, or correctness checks.  A universal pass
// licenses a separately frozen whole-Recover experiment; it is not itself a
// production result.
//
// The primary helper directly calls the narrow, no-write-on-false target-wide
// fixed-five primitive.  Unsupported targets execute the exact current
// three-pass schedule.

#include "../codec/WirehairV2Precode.h"
#include "../codec/WirehairV2Solve.h"
#include "../bench/Wh2RdpruTargetIdentityV2.h"
#include "../gf256.h"

#include <algorithm>
#include <array>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <limits>
#include <new>
#include <string>
#include <vector>

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS) || \
    defined(WIREHAIR_V2_ENABLE_BENCHMARK_EQUATIONS)
#error "the fixed5 falsifier must link production-equivalent V2 equations"
#endif
#if defined(WH_COUNT)
#error "the fixed5 falsifier forbids timed GF256 instrumentation"
#endif

#ifndef WIREHAIR_FIXED5_SOURCE_GIT_COMMIT
#define WIREHAIR_FIXED5_SOURCE_GIT_COMMIT \
    "0000000000000000000000000000000000000000"
#endif

namespace {

using Clock = std::chrono::steady_clock;

static const uint32_t kMixCount = 3u;
static const uint32_t kRowsPerDegree = 256u;
static const uint32_t kSamples = 32u;
static const uint32_t kMaximumBlockId = 1000000u;
static const uint64_t kPrecodeSeed = UINT64_C(0x786f72667573696f);
static const uint32_t kPeelSeed = UINT32_C(0x4d241359);
static const int32_t kTargetCpu = 50;
static const uint32_t kTargetFullApicId = UINT32_C(0x00000064);
static const uint64_t kTargetIdentityCanonicalBytes = UINT64_C(617);
static const char kTargetIdentitySha256[] =
    "3288e0ef61cf3e628dcd827f9cf003c9d6ec6b5a12169e7a8bfc796baacddba7";
static const double kStudentT95Df31 = 2.039513446;
static const double kControlLower = 0.99;
static const double kControlUpper = 1.01;
static const uint8_t kPrefixGuard = 0xa7u;
static const uint8_t kSuffixGuard = 0x5du;
static const size_t kGuardBytes = 64u;

#if defined(GF256_TRY_WIDE_XOR)
static const int kWideXorBuild = 1;
#else
static const int kWideXorBuild = 0;
#endif

#if defined(_MSC_VER)
#define WH2_FIXED5_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define WH2_FIXED5_NOINLINE __attribute__((noinline))
#else
#define WH2_FIXED5_NOINLINE
#endif

struct FrozenCell
{
    uint32_t Ordinal;
    uint32_t K;
    uint32_t BlockBytes;
    uint32_t Repetitions;
};

static const FrozenCell kFrozenCells[] = {
    { 0u, 10000u, 1280u, 256u },
    { 1u, 64000u, 32768u, 10u },
    { 2u, 64000u, 1280u, 256u },
    { 3u, 10000u, 32768u, 10u },
    { 4u, 10000u, 4096u, 80u },
    { 5u, 64000u, 4096u, 80u }
};

class GuardedAlignedStorage
{
public:
    GuardedAlignedStorage() = default;

    ~GuardedAlignedStorage()
    {
        std::free(Raw);
    }

    GuardedAlignedStorage(const GuardedAlignedStorage&) = delete;
    GuardedAlignedStorage& operator=(const GuardedAlignedStorage&) = delete;

    bool Allocate(size_t bytes)
    {
        if (Raw || bytes == 0u ||
            bytes > std::numeric_limits<size_t>::max() -
                (2u * kGuardBytes + 63u))
        {
            return false;
        }
        const size_t allocation_bytes = bytes + 2u * kGuardBytes + 63u;
        Raw = std::malloc(allocation_bytes);
        if (!Raw) {
            return false;
        }
        const uintptr_t after_prefix =
            reinterpret_cast<uintptr_t>(Raw) + kGuardBytes;
        const uintptr_t aligned = (after_prefix + 63u) & ~uintptr_t(63u);
        DataValue = reinterpret_cast<uint8_t*>(aligned);
        BytesValue = bytes;
        FillGuards();
        return true;
    }

    void FillGuards()
    {
        std::memset(DataValue - kGuardBytes, kPrefixGuard, kGuardBytes);
        std::memset(DataValue + BytesValue, kSuffixGuard, kGuardBytes);
    }

    bool GuardsIntact() const
    {
        for (size_t i = 0u; i < kGuardBytes; ++i)
        {
            if (DataValue[-static_cast<std::ptrdiff_t>(i + 1u)] !=
                    kPrefixGuard ||
                DataValue[BytesValue + i] != kSuffixGuard)
            {
                return false;
            }
        }
        return true;
    }

    uint8_t* Data() { return DataValue; }
    const uint8_t* Data() const { return DataValue; }
    size_t Bytes() const { return BytesValue; }

private:
    void* Raw = nullptr;
    uint8_t* DataValue = nullptr;
    size_t BytesValue = 0u;
};

struct RowFixture
{
    uint32_t BlockId = 0u;
    uint32_t PeelCount = 0u;
    uint32_t TermCount = 0u;
    std::array<uint32_t, 6u> Columns = {{ 0u, 0u, 0u, 0u, 0u, 0u }};
    std::array<const void*, 6u> Sources = {{
        nullptr, nullptr, nullptr, nullptr, nullptr, nullptr
    }};
};

typedef void (*Schedule)(
    const RowFixture&, uint8_t*, uint32_t, uint64_t*);

struct RawPair
{
    uint64_t LegacyNs = 0u;
    uint64_t CandidateNs = 0u;
    bool LegacyFirst = true;
};

struct PanelStats
{
    double GeometricMean = 0.0;
    double Lower95 = 0.0;
    double Upper95 = 0.0;
};

struct PanelResult
{
    const char* Name = nullptr;
    const char* Candidate = nullptr;
    std::vector<RawPair> Pairs;
    PanelStats Stats;
};

bool RangesOverlap(
    const void* first,
    size_t first_bytes,
    const void* second,
    size_t second_bytes)
{
    const uintptr_t a = reinterpret_cast<uintptr_t>(first);
    const uintptr_t b = reinterpret_cast<uintptr_t>(second);
    const uintptr_t limit = std::numeric_limits<uintptr_t>::max();
    if (first_bytes > limit - a || second_bytes > limit - b) {
        return true;
    }
    return a < b + second_bytes && b < a + first_bytes;
}

uint64_t HashBytes(const uint8_t* data, size_t bytes)
{
    uint64_t hash = UINT64_C(1469598103934665603);
    for (size_t i = 0u; i < bytes; ++i) {
        hash = (hash ^ data[i]) * UINT64_C(1099511628211);
    }
    return hash;
}

void HashU32(uint64_t& hash, uint32_t value)
{
    for (unsigned shift = 0u; shift != 32u; shift += 8u) {
        hash = (hash ^ static_cast<uint8_t>(value >> shift)) *
            UINT64_C(1099511628211);
    }
}

uint64_t HashRows(const std::vector<RowFixture>& rows)
{
    uint64_t hash = UINT64_C(1469598103934665603);
    for (const RowFixture& row : rows)
    {
        HashU32(hash, row.BlockId);
        HashU32(hash, row.PeelCount);
        HashU32(hash, row.TermCount);
        for (uint32_t term = 0u; term < row.TermCount; ++term) {
            HashU32(hash, row.Columns[term]);
        }
    }
    return hash;
}

void FillDeterministic(uint8_t* data, size_t bytes)
{
    for (size_t i = 0u; i < bytes; ++i) {
        data[i] = static_cast<uint8_t>(
            i * size_t(131u) + (i >> 13u) + size_t(17u));
    }
}

class ScopedWideXor
{
public:
    explicit ScopedWideXor(bool enable)
        : Previous(gf256_set_thread_wide_xor(enable ? 1 : 0))
    {
    }

    ~ScopedWideXor()
    {
        (void)gf256_set_thread_wide_xor(Previous);
    }

    ScopedWideXor(const ScopedWideXor&) = delete;
    ScopedWideXor& operator=(const ScopedWideXor&) = delete;

private:
    int Previous;
};

WH2_FIXED5_NOINLINE void EvaluateLegacyCurrent(
    const RowFixture& row,
    uint8_t* output,
    uint32_t block_bytes,
    uint64_t* block_ops_out)
{
    if (row.TermCount == 5u)
    {
        // Exact current MIX3, PeelCount=2 production schedule: three full
        // destination passes and five logical source operations.
        gf256_addset_mem(
            output, row.Sources[0], row.Sources[1], (int)block_bytes);
        gf256_add_mem(output, row.Sources[2], (int)block_bytes);
        gf256_add2_mem(
            output, row.Sources[3], row.Sources[4], (int)block_bytes);
    }
    else
    {
        // Exact current K>=10000, B<=32768, six-term set-XOR route.  Frozen
        // fallback panels run with the wide selector disabled.
        gf256_addset_multi_mem(
            output, row.Sources.data(), (int)row.TermCount,
            (int)block_bytes);
    }
    if (block_ops_out) {
        *block_ops_out = row.TermCount;
    }
}

// Sole call site for the direct target-wide primitive.
WH2_FIXED5_NOINLINE bool TryFixed5Wide(
    const RowFixture& row,
    uint8_t* output,
    uint32_t block_bytes)
{
    return gf256_try_addset5_wide_mem(
        output, row.Sources.data(), (int)block_bytes) != 0;
}

WH2_FIXED5_NOINLINE void EvaluateFixed5Wide(
    const RowFixture& row,
    uint8_t* output,
    uint32_t block_bytes,
    uint64_t* block_ops_out)
{
    if (row.TermCount == 5u)
    {
        if (!TryFixed5Wide(row, output, block_bytes))
        {
            gf256_addset_mem(
                output, row.Sources[0], row.Sources[1],
                (int)block_bytes);
            gf256_add_mem(output, row.Sources[2], (int)block_bytes);
            gf256_add2_mem(
                output, row.Sources[3], row.Sources[4],
                (int)block_bytes);
        }
    }
    else
    {
        gf256_addset_multi_mem(
            output, row.Sources.data(), (int)row.TermCount,
            (int)block_bytes);
    }
    if (block_ops_out) {
        *block_ops_out = row.TermCount;
    }
}

WH2_FIXED5_NOINLINE void EvaluateSplit23Wide(
    const RowFixture& row,
    uint8_t* output,
    uint32_t block_bytes,
    uint64_t* block_ops_out)
{
    if (row.TermCount == 5u)
    {
        gf256_addset_mem(
            output, row.Sources[0], row.Sources[1], (int)block_bytes);
        gf256_add_multi_mem(
            output, row.Sources.data() + 2u, 3, (int)block_bytes);
    }
    else
    {
        gf256_addset_multi_mem(
            output, row.Sources.data(), (int)row.TermCount,
            (int)block_bytes);
    }
    if (block_ops_out) {
        *block_ops_out = row.TermCount;
    }
}

void ScalarOracle(
    const RowFixture& row,
    uint8_t* output,
    uint32_t block_bytes)
{
    for (uint32_t byte = 0u; byte < block_bytes; ++byte)
    {
        uint8_t value = 0u;
        for (uint32_t term = 0u; term < row.TermCount; ++term) {
            value ^= static_cast<const uint8_t*>(row.Sources[term])[byte];
        }
        output[byte] = value;
    }
}

bool ParseU32(const char* text, uint32_t& value)
{
    if (!text || !*text || *text == '-') {
        return false;
    }
    errno = 0;
    char* end = nullptr;
    const unsigned long parsed = std::strtoul(text, &end, 10);
    if (errno != 0 || !end || *end != '\0' ||
        parsed > std::numeric_limits<uint32_t>::max())
    {
        return false;
    }
    value = static_cast<uint32_t>(parsed);
    return true;
}

bool CaptureExactTargetIdentity(
    wirehair_wh2_bench::TargetIdentityReceiptV2& receipt,
    std::string& diagnostic)
{
    if (!wirehair_wh2_bench::CaptureNativeTargetIdentityV2(
            receipt, diagnostic))
    {
        return false;
    }
    if (!wirehair_wh2_bench::ValidateFrozenCpu50IdentityV2(
            receipt, diagnostic))
    {
        return false;
    }
    if (receipt.RequestedCpu != kTargetCpu ||
        receipt.Derived.FullApicId != kTargetFullApicId ||
        receipt.CanonicalSha256 != kTargetIdentitySha256 ||
        wirehair_wh2_bench::FrozenCpu50TargetIdentityCanonicalBytesV2() !=
            kTargetIdentityCanonicalBytes ||
        receipt.Before.Cpu != kTargetCpu ||
        receipt.After.Cpu != kTargetCpu ||
        receipt.Before.Affinity.size() != 1u ||
        receipt.After.Affinity.size() != 1u ||
        receipt.Before.Affinity[0] != kTargetCpu ||
        receipt.After.Affinity[0] != kTargetCpu)
    {
        diagnostic = "frozen target receipt identity mismatch";
        return false;
    }
    diagnostic.clear();
    return true;
}

const FrozenCell* FindFrozenCell(
    uint32_t ordinal,
    uint32_t K,
    uint32_t block_bytes)
{
    for (const FrozenCell& cell : kFrozenCells) {
        if (cell.Ordinal == ordinal && cell.K == K &&
            cell.BlockBytes == block_bytes)
        {
            return &cell;
        }
    }
    return nullptr;
}

bool PrepareRows(
    uint32_t K,
    uint32_t P,
    uint32_t block_bytes,
    const wirehair_v2::PacketRowConfig& config,
    const wirehair_v2::PacketRowRuntime& runtime,
    const GuardedAlignedStorage& intermediate,
    std::vector<RowFixture>& active,
    std::vector<RowFixture>& fallback)
{
    active.clear();
    fallback.clear();
    active.reserve(kRowsPerDegree);
    fallback.reserve(kRowsPerDegree);
    for (uint32_t block_id = 0u;
         block_id < kMaximumBlockId &&
            (active.size() != kRowsPerDegree ||
             fallback.size() != kRowsPerDegree);
         ++block_id)
    {
        const std::vector<uint32_t> columns =
            wirehair_v2::GeneratePacketMatrixRowWithRuntime(
                K, P, block_id, config, runtime);
        if (columns.size() != 5u && columns.size() != 6u) {
            continue;
        }
        const uint32_t peel_count =
            static_cast<uint32_t>(columns.size()) - kMixCount;
        std::vector<RowFixture>* destination = nullptr;
        if (peel_count == 2u && active.size() != kRowsPerDegree) {
            destination = &active;
        }
        else if (peel_count == 3u && fallback.size() != kRowsPerDegree) {
            destination = &fallback;
        }
        if (!destination) {
            continue;
        }

        RowFixture row;
        row.BlockId = block_id;
        row.PeelCount = peel_count;
        row.TermCount = static_cast<uint32_t>(columns.size());
        for (uint32_t term = 0u; term < row.TermCount; ++term)
        {
            const uint32_t column = columns[term];
            if (column >= K + P) {
                return false;
            }
            for (uint32_t previous = 0u; previous < term; ++previous) {
                if (columns[previous] == column) {
                    return false;
                }
            }
            row.Columns[term] = column;
            const size_t offset = (size_t)column * block_bytes;
            if (offset > intermediate.Bytes() ||
                block_bytes > intermediate.Bytes() - offset)
            {
                return false;
            }
            row.Sources[term] = intermediate.Data() + offset;
        }
        destination->push_back(row);
    }
    return active.size() == kRowsPerDegree &&
        fallback.size() == kRowsPerDegree;
}

bool VerifyStorageState(
    const GuardedAlignedStorage& intermediate,
    uint64_t expected_input_hash,
    const GuardedAlignedStorage& first_output,
    const GuardedAlignedStorage& second_output,
    const GuardedAlignedStorage& oracle_output)
{
    return intermediate.GuardsIntact() && first_output.GuardsIntact() &&
        second_output.GuardsIntact() && oracle_output.GuardsIntact() &&
        HashBytes(intermediate.Data(), intermediate.Bytes()) ==
            expected_input_hash;
}

bool VerifySchedule(
    const std::vector<RowFixture>& rows,
    Schedule schedule,
    uint32_t block_bytes,
    GuardedAlignedStorage& output,
    GuardedAlignedStorage& oracle)
{
    static const uint32_t kOracleRepetitions[] = { 1u, 2u };
    for (uint32_t repetitions : kOracleRepetitions)
    {
        for (const RowFixture& row : rows)
        {
            uint64_t operations = UINT64_MAX;
            ScalarOracle(row, oracle.Data(), block_bytes);
            for (uint32_t repetition = 0u;
                 repetition < repetitions;
                 ++repetition)
            {
                operations = UINT64_MAX;
                schedule(row, output.Data(), block_bytes, &operations);
                if (operations != row.TermCount) {
                    return false;
                }
            }
            if (std::memcmp(
                    output.Data(), oracle.Data(), block_bytes) != 0)
            {
                return false;
            }
        }
    }
    return output.GuardsIntact() && oracle.GuardsIntact();
}

bool VerifyAllRows(
    const std::vector<RowFixture>& active,
    const std::vector<RowFixture>& fallback,
    uint32_t block_bytes,
    GuardedAlignedStorage& first_output,
    GuardedAlignedStorage& second_output,
    GuardedAlignedStorage& oracle_output)
{
    // Exercise portable/current semantics with the wide selector disabled.
    {
        ScopedWideXor wide_off(false);
        if (!VerifySchedule(
                active, EvaluateLegacyCurrent, block_bytes,
                first_output, oracle_output) ||
            !VerifySchedule(
                active, EvaluateFixed5Wide, block_bytes,
                second_output, oracle_output) ||
            !VerifySchedule(
                active, EvaluateSplit23Wide, block_bytes,
                second_output, oracle_output) ||
            !VerifySchedule(
                fallback, EvaluateLegacyCurrent, block_bytes,
                first_output, oracle_output) ||
            !VerifySchedule(
                fallback, EvaluateFixed5Wide, block_bytes,
                second_output, oracle_output))
        {
            return false;
        }
    }

    // Exercise the selected target-wide path when present.  Unsupported
    // builds intentionally accept the selector as a no-op and must still be
    // byte/work identical in this correctness-only self-check.
    {
        ScopedWideXor wide_on(true);
        if (!VerifySchedule(
                active, EvaluateFixed5Wide, block_bytes,
                first_output, oracle_output) ||
            !VerifySchedule(
                active, EvaluateSplit23Wide, block_bytes,
                second_output, oracle_output))
        {
            return false;
        }
    }
    return true;
}

enum class DirectWideCheck
{
    Available,
    Unavailable,
    Incorrect
};

DirectWideCheck VerifyDirectWideRows(
    const std::vector<RowFixture>& active,
    uint32_t block_bytes,
    GuardedAlignedStorage& output,
    GuardedAlignedStorage& oracle)
{
    for (const RowFixture& row : active)
    {
        std::memset(output.Data(), 0x6b, block_bytes);
        const uint64_t before = HashBytes(output.Data(), block_bytes);
        if (!TryFixed5Wide(row, output.Data(), block_bytes))
        {
            return output.GuardsIntact() &&
                    HashBytes(output.Data(), block_bytes) == before ?
                DirectWideCheck::Unavailable : DirectWideCheck::Incorrect;
        }
        ScalarOracle(row, oracle.Data(), block_bytes);
        if (std::memcmp(output.Data(), oracle.Data(), block_bytes) != 0 ||
            !output.GuardsIntact() || !oracle.GuardsIntact())
        {
            return DirectWideCheck::Incorrect;
        }
    }
    return DirectWideCheck::Available;
}

bool VerifyDirectTryContract(
    const RowFixture& row,
    uint32_t block_bytes,
    GuardedAlignedStorage& output,
    GuardedAlignedStorage& oracle)
{
    std::memset(output.Data(), 0x93, block_bytes);
    const uint64_t untouched_hash = HashBytes(output.Data(), block_bytes);
    const bool available = TryFixed5Wide(row, output.Data(), block_bytes);
    if (!available) {
        return output.GuardsIntact() &&
            HashBytes(output.Data(), block_bytes) == untouched_hash;
    }
    ScalarOracle(row, oracle.Data(), block_bytes);
    return output.GuardsIntact() && oracle.GuardsIntact() &&
        std::memcmp(output.Data(), oracle.Data(), block_bytes) == 0;
}

uint64_t MeasureArm(
    const std::vector<RowFixture>& rows,
    Schedule schedule,
    uint32_t block_bytes,
    uint32_t repetitions,
    GuardedAlignedStorage& output)
{
    const Clock::time_point start = Clock::now();
    for (uint32_t repetition = 0u; repetition < repetitions; ++repetition) {
        for (const RowFixture& row : rows) {
            // A null work pointer guarantees that the measured path performs
            // no diagnostic counter write.
            schedule(row, output.Data(), block_bytes, nullptr);
        }
    }
    const Clock::time_point finish = Clock::now();
    return static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            finish - start).count());
}

bool VerifyLastOutput(
    const std::vector<RowFixture>& rows,
    uint32_t block_bytes,
    const GuardedAlignedStorage& output,
    GuardedAlignedStorage& oracle)
{
    ScalarOracle(rows.back(), oracle.Data(), block_bytes);
    return std::memcmp(output.Data(), oracle.Data(), block_bytes) == 0 &&
        output.GuardsIntact() && oracle.GuardsIntact();
}

bool ComputeStats(const std::vector<RawPair>& pairs, PanelStats& stats)
{
    if (pairs.size() != kSamples) {
        return false;
    }
    std::array<double, kSamples> logs = {{ 0.0 }};
    double mean = 0.0;
    for (size_t i = 0u; i < pairs.size(); ++i)
    {
        if (pairs[i].LegacyNs == 0u || pairs[i].CandidateNs == 0u) {
            return false;
        }
        const double ratio = static_cast<double>(pairs[i].CandidateNs) /
            static_cast<double>(pairs[i].LegacyNs);
        if (!(ratio > 0.0) || !std::isfinite(ratio)) {
            return false;
        }
        logs[i] = std::log(ratio);
        mean += logs[i];
    }
    mean /= static_cast<double>(pairs.size());
    double squared_deviation = 0.0;
    for (double value : logs) {
        const double delta = value - mean;
        squared_deviation += delta * delta;
    }
    const double sample_variance = squared_deviation /
        static_cast<double>(pairs.size() - 1u);
    const double half_width = kStudentT95Df31 *
        std::sqrt(sample_variance / static_cast<double>(pairs.size()));
    stats.GeometricMean = std::exp(mean);
    stats.Lower95 = std::exp(mean - half_width);
    stats.Upper95 = std::exp(mean + half_width);
    return std::isfinite(stats.GeometricMean) &&
        std::isfinite(stats.Lower95) && std::isfinite(stats.Upper95);
}

bool RunPanel(
    const char* name,
    const char* candidate,
    const std::vector<RowFixture>& rows,
    Schedule legacy_schedule,
    Schedule candidate_schedule,
    uint32_t block_bytes,
    uint32_t repetitions,
    uint64_t expected_input_hash,
    const GuardedAlignedStorage& intermediate,
    GuardedAlignedStorage& legacy_output,
    GuardedAlignedStorage& candidate_output,
    GuardedAlignedStorage& oracle_output,
    PanelResult& result)
{
    result.Name = name;
    result.Candidate = candidate;
    result.Pairs.clear();
    result.Pairs.reserve(kSamples);

    // Exactly one full-workload warmup per arm, outside all recorded samples.
    (void)MeasureArm(
        rows, legacy_schedule, block_bytes, repetitions, legacy_output);
    (void)MeasureArm(
        rows, candidate_schedule, block_bytes, repetitions, candidate_output);
    if (!VerifyLastOutput(
            rows, block_bytes, legacy_output, oracle_output) ||
        !VerifyLastOutput(
            rows, block_bytes, candidate_output, oracle_output) ||
        !VerifyStorageState(
            intermediate, expected_input_hash, legacy_output,
            candidate_output, oracle_output))
    {
        return false;
    }

    for (uint32_t sample = 0u; sample < kSamples; ++sample)
    {
        RawPair pair;
        pair.LegacyFirst = (sample & 1u) == 0u;
        if (pair.LegacyFirst)
        {
            pair.LegacyNs = MeasureArm(
                rows, legacy_schedule, block_bytes, repetitions,
                legacy_output);
            pair.CandidateNs = MeasureArm(
                rows, candidate_schedule, block_bytes, repetitions,
                candidate_output);
        }
        else
        {
            pair.CandidateNs = MeasureArm(
                rows, candidate_schedule, block_bytes, repetitions,
                candidate_output);
            pair.LegacyNs = MeasureArm(
                rows, legacy_schedule, block_bytes, repetitions,
                legacy_output);
        }
        if (pair.LegacyNs == 0u || pair.CandidateNs == 0u) {
            return false;
        }
        result.Pairs.push_back(pair);
    }

    return VerifyLastOutput(
               rows, block_bytes, legacy_output, oracle_output) &&
        VerifyLastOutput(
               rows, block_bytes, candidate_output, oracle_output) &&
        VerifyStorageState(
               intermediate, expected_input_hash, legacy_output,
               candidate_output, oracle_output) &&
        ComputeStats(result.Pairs, result.Stats);
}

void EmitPanel(
    const FrozenCell& cell,
    const PanelResult& result,
    uint32_t degree)
{
    for (size_t sample = 0u; sample < result.Pairs.size(); ++sample)
    {
        const RawPair& pair = result.Pairs[sample];
        const double ratio = static_cast<double>(pair.CandidateNs) /
            static_cast<double>(pair.LegacyNs);
        std::printf(
            "wh2_fixed5_raw,version=1,cell=%s,candidate=%s,ordinal=%u,"
            "K=%u,block_bytes=%u,degree=%u,rows=%u,repetitions=%u,"
            "sample=%zu,order=%s,legacy_ns=%" PRIu64
            ",candidate_ns=%" PRIu64 ",ratio=%.12f\n",
            result.Name, result.Candidate, cell.Ordinal, cell.K,
            cell.BlockBytes, degree, kRowsPerDegree, cell.Repetitions,
            sample, pair.LegacyFirst ? "AB" : "BA", pair.LegacyNs,
            pair.CandidateNs, ratio);
    }
    std::printf(
        "wh2_fixed5_summary,version=1,cell=%s,candidate=%s,ordinal=%u,"
        "K=%u,block_bytes=%u,degree=%u,samples=%u,geomean=%.12f,"
        "lower95=%.12f,upper95=%.12f\n",
        result.Name, result.Candidate, cell.Ordinal, cell.K,
        cell.BlockBytes, degree, kSamples, result.Stats.GeometricMean,
        result.Stats.Lower95, result.Stats.Upper95);
}

bool ControlPasses(const PanelResult& panel, bool require_one)
{
    return panel.Stats.Lower95 >= kControlLower &&
        panel.Stats.Upper95 <= kControlUpper &&
        (!require_one ||
         (panel.Stats.Lower95 <= 1.0 && panel.Stats.Upper95 >= 1.0));
}

void EmitFrozenDescription()
{
    std::printf(
        "wh2_fixed5_design,version=1,source_commit=%s,rows=%u,samples=%u,"
        "mix_count=%u,peel_seed=%08" PRIx32
        ",precode_seed=%016" PRIx64
        ",target_cpu=%d,target_full_apic_id=%08" PRIx32
        ",target_identity_sha256=%s"
        ",candidate=direct_no_write_fixed5_try,"
        "primary_gate=upper95_lt_0.99,control_band=0.99_to_1.01,"
        "student_t_df31=%.9f,timed_scope=xor_schedules_only,"
        "timed_counters=0,binary_sha256=controller_required\n",
        WIREHAIR_FIXED5_SOURCE_GIT_COMMIT, kRowsPerDegree, kSamples,
        kMixCount, kPeelSeed, kPrecodeSeed, kTargetCpu,
        kTargetFullApicId, kTargetIdentitySha256, kStudentT95Df31);
    for (const FrozenCell& cell : kFrozenCells) {
        std::printf(
            "wh2_fixed5_cell,version=1,ordinal=%u,K=%u,block_bytes=%u,"
            "repetitions=%u\n",
            cell.Ordinal, cell.K, cell.BlockBytes, cell.Repetitions);
    }
}

bool RunSelfTest()
{
    static const uint32_t kSizes[] = {
        1u, 15u, 16u, 31u, 32u, 63u, 64u, 127u, 128u, 257u
    };
    for (uint32_t block_bytes : kSizes)
    {
        GuardedAlignedStorage inputs;
        GuardedAlignedStorage first_output;
        GuardedAlignedStorage second_output;
        GuardedAlignedStorage oracle;
        if (!inputs.Allocate((size_t)6u * block_bytes) ||
            !first_output.Allocate(block_bytes) ||
            !second_output.Allocate(block_bytes) ||
            !oracle.Allocate(block_bytes))
        {
            return false;
        }
        FillDeterministic(inputs.Data(), inputs.Bytes());
        const uint64_t input_hash = HashBytes(inputs.Data(), inputs.Bytes());
        RowFixture five;
        five.PeelCount = 2u;
        five.TermCount = 5u;
        RowFixture six;
        six.PeelCount = 3u;
        six.TermCount = 6u;
        for (uint32_t term = 0u; term < 6u; ++term)
        {
            six.Columns[term] = term;
            six.Sources[term] = inputs.Data() + (size_t)term * block_bytes;
            if (term < 5u)
            {
                five.Columns[term] = term;
                five.Sources[term] = six.Sources[term];
            }
        }
        const std::vector<RowFixture> active(1u, five);
        const std::vector<RowFixture> fallback(1u, six);
        if (!VerifyDirectTryContract(
                active.front(), block_bytes, first_output, oracle) ||
            !VerifyAllRows(
                active, fallback, block_bytes, first_output,
                second_output, oracle) ||
            !inputs.GuardsIntact() ||
            HashBytes(inputs.Data(), inputs.Bytes()) != input_hash)
        {
            return false;
        }
    }
    const wirehair_wh2_bench::TargetIdentityReceiptV2 frozen_target =
        wirehair_wh2_bench::FrozenCpu50TargetIdentityReceiptV2();
    return frozen_target.RequestedCpu == kTargetCpu &&
        frozen_target.Derived.FullApicId == kTargetFullApicId &&
        frozen_target.CanonicalSha256 == kTargetIdentitySha256 &&
        wirehair_wh2_bench::FrozenCpu50TargetIdentitySha256V2() ==
            kTargetIdentitySha256;
}

int RunFrozenCell(
    const FrozenCell& cell,
    const wirehair_wh2_bench::TargetIdentityReceiptV2& pre_target)
{
    gf256_x86_cpu_features features = {};
    gf256_get_active_x86_cpu_features(&features);
    if (kWideXorBuild != 1 || features.AVX2 != 1 || features.AVX512 != 1)
    {
        std::fprintf(stderr,
            "frozen AVX-512 stratum unavailable: wide_build=%d avx2=%d "
            "avx512=%d\n",
            kWideXorBuild, features.AVX2, features.AVX512);
        return 3;
    }

    const wirehair_v2::PrecodeParams params =
        wirehair_v2::MakeCertifiedParams(cell.K, kPrecodeSeed);
    const uint64_t P_wide = (uint64_t)params.Staircase + params.DenseRows +
        params.HeavyRows;
    if (P_wide > UINT32_MAX) {
        std::fprintf(stderr, "invalid certified precode width\n");
        return 1;
    }
    const uint32_t P = static_cast<uint32_t>(P_wide);
    wirehair_v2::PacketRowConfig config;
    config.PeelSeed = kPeelSeed;
    config.MixCount = kMixCount;
    wirehair_v2::PacketRowRuntime runtime;
    if (!runtime.Initialize(cell.K, P, config.MixCount)) {
        std::fprintf(stderr, "packet runtime initialization failed\n");
        return 1;
    }

    const uint64_t total_bytes_wide =
        ((uint64_t)cell.K + P) * cell.BlockBytes;
    if (total_bytes_wide >
        static_cast<uint64_t>(std::numeric_limits<size_t>::max()))
    {
        std::fprintf(stderr, "intermediate arena size overflow\n");
        return 1;
    }
    GuardedAlignedStorage intermediate;
    GuardedAlignedStorage legacy_output;
    GuardedAlignedStorage candidate_output;
    GuardedAlignedStorage oracle_output;
    if (!intermediate.Allocate(static_cast<size_t>(total_bytes_wide)) ||
        !legacy_output.Allocate(cell.BlockBytes) ||
        !candidate_output.Allocate(cell.BlockBytes) ||
        !oracle_output.Allocate(cell.BlockBytes))
    {
        std::fprintf(stderr, "bounded arena allocation failed\n");
        return 1;
    }
    if (RangesOverlap(
            intermediate.Data(), intermediate.Bytes(), legacy_output.Data(),
            legacy_output.Bytes()) ||
        RangesOverlap(
            intermediate.Data(), intermediate.Bytes(),
            candidate_output.Data(), candidate_output.Bytes()) ||
        RangesOverlap(
            intermediate.Data(), intermediate.Bytes(), oracle_output.Data(),
            oracle_output.Bytes()) ||
        RangesOverlap(
            legacy_output.Data(), legacy_output.Bytes(),
            candidate_output.Data(), candidate_output.Bytes()) ||
        RangesOverlap(
            legacy_output.Data(), legacy_output.Bytes(),
            oracle_output.Data(), oracle_output.Bytes()) ||
        RangesOverlap(
            candidate_output.Data(), candidate_output.Bytes(),
            oracle_output.Data(), oracle_output.Bytes()))
    {
        std::fprintf(stderr, "benchmark storage overlap\n");
        return 1;
    }
    FillDeterministic(intermediate.Data(), intermediate.Bytes());
    std::memset(legacy_output.Data(), 0x39, legacy_output.Bytes());
    std::memset(candidate_output.Data(), 0xc6, candidate_output.Bytes());
    std::memset(oracle_output.Data(), 0x71, oracle_output.Bytes());
    const uint64_t input_hash =
        HashBytes(intermediate.Data(), intermediate.Bytes());

    std::vector<RowFixture> active;
    std::vector<RowFixture> fallback;
    if (!PrepareRows(
            cell.K, P, cell.BlockBytes, config, runtime, intermediate,
            active, fallback))
    {
        std::fprintf(stderr, "exact-degree fixture generation failed\n");
        return 1;
    }
    if (!VerifyAllRows(
            active, fallback, cell.BlockBytes, legacy_output,
            candidate_output, oracle_output) ||
        !VerifyStorageState(
            intermediate, input_hash, legacy_output, candidate_output,
            oracle_output))
    {
        std::fprintf(stderr, "pre-timing byte/work/guard/input oracle failed\n");
        return 1;
    }

    const DirectWideCheck direct_check = VerifyDirectWideRows(
        active, cell.BlockBytes, candidate_output, oracle_output);
    if (direct_check != DirectWideCheck::Available)
    {
        if (direct_check == DirectWideCheck::Unavailable)
        {
            std::fprintf(stderr,
                "primary direct target-wide route unavailable: "
                "wide_build=%d avx2=%d\n",
                kWideXorBuild, features.AVX2);
            return 3;
        }
        std::fprintf(stderr, "direct target-wide row oracle failed\n");
        return 1;
    }
    if (!VerifyStorageState(
            intermediate, input_hash, legacy_output, candidate_output,
            oracle_output))
    {
        std::fprintf(stderr, "direct target-wide availability probe failed\n");
        return 1;
    }

    std::printf(
        "wh2_fixed5_meta,version=1,source_commit=%s,ordinal=%u,K=%u,"
        "P=%u,block_bytes=%u,repetitions=%u,rows=%u,samples=%u,"
        "mix_count=%u,peel_seed=%08" PRIx32
        ",precode_seed=%016" PRIx64
        ",active_degree=2,fallback_degree=3,active_terms=5,"
        "fallback_terms=6,active_rows_hash=%016" PRIx64
        ",fallback_rows_hash=%016" PRIx64
        ",input_hash=%016" PRIx64
        ",active_first_id=%u,active_last_id=%u,"
        "fallback_first_id=%u,fallback_last_id=%u,wide_build=%d,"
        "avx2=%d,avx512=%d,candidate=direct_no_write_fixed5_try,"
        "direct_wide_available=1,"
        "timed_row_generation=0,timed_pointer_gather=0,timed_counters=0,"
        "panel_order=active_fallback_aa_split\n",
        WIREHAIR_FIXED5_SOURCE_GIT_COMMIT, cell.Ordinal, cell.K, P,
        cell.BlockBytes, cell.Repetitions, kRowsPerDegree, kSamples,
        kMixCount, kPeelSeed, kPrecodeSeed,
        HashRows(active), HashRows(fallback),
        input_hash, active.front().BlockId, active.back().BlockId,
        fallback.front().BlockId, fallback.back().BlockId, kWideXorBuild,
        features.AVX2, features.AVX512);

    PanelResult primary;
    {
        if (!RunPanel(
                "active", "direct_fixed5", active,
                EvaluateLegacyCurrent, EvaluateFixed5Wide,
                cell.BlockBytes, cell.Repetitions, input_hash, intermediate,
                legacy_output, candidate_output, oracle_output, primary))
        {
            std::fprintf(stderr, "active panel validation failed\n");
            return 1;
        }
    }
    PanelResult fallback_panel;
    {
        // Production Recover leaves the wide selector disabled.  Both arms
        // therefore execute the exact existing six-term set-XOR path.
        ScopedWideXor wide_off(false);
        if (!RunPanel(
                "fallback", "direct_fixed5_dispatch", fallback,
                EvaluateLegacyCurrent, EvaluateFixed5Wide,
                cell.BlockBytes, cell.Repetitions, input_hash, intermediate,
                legacy_output, candidate_output, oracle_output,
                fallback_panel))
        {
            std::fprintf(stderr, "fallback panel validation failed\n");
            return 1;
        }
    }
    PanelResult aa;
    {
        ScopedWideXor wide_off(false);
        if (!RunPanel(
                "aa", "legacy_copy", active,
                EvaluateLegacyCurrent, EvaluateLegacyCurrent,
                cell.BlockBytes, cell.Repetitions, input_hash, intermediate,
                legacy_output, candidate_output, oracle_output, aa))
        {
            std::fprintf(stderr, "A/A panel validation failed\n");
            return 1;
        }
    }
    PanelResult split;
    {
        ScopedWideXor wide_on(true);
        if (!RunPanel(
                "split", "wide_split2_plus3", active,
                EvaluateLegacyCurrent, EvaluateSplit23Wide,
                cell.BlockBytes, cell.Repetitions, input_hash, intermediate,
                legacy_output, candidate_output, oracle_output, split))
        {
            std::fprintf(stderr, "split diagnostic validation failed\n");
            return 1;
        }
    }
    wirehair_wh2_bench::TargetIdentityReceiptV2 post_target;
    std::string target_diagnostic;
    if (!CaptureExactTargetIdentity(post_target, target_diagnostic))
    {
        std::fprintf(stderr, "post-timing target identity failed: %s\n",
            target_diagnostic.c_str());
        return 1;
    }
    // Publish no scientific timing or CI output until the post-timing target
    // receipt has authenticated the same singleton CPU/APIC identity.
    EmitPanel(cell, primary, 2u);
    EmitPanel(cell, fallback_panel, 3u);
    EmitPanel(cell, aa, 2u);
    EmitPanel(cell, split, 2u);
    std::printf(
        "wh2_fixed5_target,version=1,source_commit=%s,ordinal=%u,K=%u,"
        "block_bytes=%u,"
        "target_cpu=%d,full_apic_id=%08" PRIx32
        ",pre_identity_sha256=%s,post_identity_sha256=%s,"
        "canonical_bytes=%" PRIu64 ","
        "pre_before_cpu=%d,pre_after_cpu=%d,post_before_cpu=%d,"
        "post_after_cpu=%d,pre_before_affinity_count=%zu,"
        "pre_after_affinity_count=%zu,post_before_affinity_count=%zu,"
        "post_after_affinity_count=%zu,pre_voluntary_delta=%" PRId64
        ",pre_involuntary_delta=%" PRId64
        ",post_voluntary_delta=%" PRId64
        ",post_involuntary_delta=%" PRId64 ",gate=pass\n",
        WIREHAIR_FIXED5_SOURCE_GIT_COMMIT, cell.Ordinal, cell.K,
        cell.BlockBytes, kTargetCpu,
        post_target.Derived.FullApicId, pre_target.CanonicalSha256.c_str(),
        post_target.CanonicalSha256.c_str(), kTargetIdentityCanonicalBytes,
        pre_target.Before.Cpu, pre_target.After.Cpu, post_target.Before.Cpu,
        post_target.After.Cpu,
        pre_target.Before.Affinity.size(), pre_target.After.Affinity.size(),
        post_target.Before.Affinity.size(), post_target.After.Affinity.size(),
        pre_target.After.VoluntaryContextSwitches -
            pre_target.Before.VoluntaryContextSwitches,
        pre_target.After.InvoluntaryContextSwitches -
            pre_target.Before.InvoluntaryContextSwitches,
        post_target.After.VoluntaryContextSwitches -
            post_target.Before.VoluntaryContextSwitches,
        post_target.After.InvoluntaryContextSwitches -
            post_target.Before.InvoluntaryContextSwitches);

    const bool primary_pass = primary.Stats.Upper95 < 0.99;
    const bool fallback_pass = ControlPasses(fallback_panel, false);
    const bool aa_pass = ControlPasses(aa, true);
    const bool controls_pass = fallback_pass && aa_pass;
    std::printf(
        "wh2_fixed5_result,version=1,ordinal=%u,K=%u,block_bytes=%u,"
        "primary=%s,fallback_control=%s,aa_control=%s,controls=%s,"
        "split_promotional=0,mismatch_sink=0,status=%s\n",
        cell.Ordinal, cell.K, cell.BlockBytes,
        primary_pass ? "pass" : "reject",
        fallback_pass ? "pass" : "invalid",
        aa_pass ? "pass" : "invalid",
        controls_pass ? "pass" : "invalid",
        !controls_pass ? "invalid" :
            (primary_pass ? "pass" : "reject"));
    if (!controls_pass) {
        return 11;
    }
    return primary_pass ? 0 : 10;
}

void PrintUsage(const char* program)
{
    std::fprintf(stderr,
        "usage: %s --self-test | --describe | "
        "--run-cell ORDINAL K BLOCK_BYTES EXPECTED_SOURCE_COMMIT TARGET_CPU\n",
        program);
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        if (argc == 2 && std::strcmp(argv[1], "--describe") == 0)
        {
            EmitFrozenDescription();
            return 0;
        }
        if (argc == 2 && std::strcmp(argv[1], "--self-test") == 0)
        {
            if (gf256_init() != 0 || !RunSelfTest()) {
                std::fprintf(stderr, "fixed5 falsifier self-test failed\n");
                return 1;
            }
            gf256_x86_cpu_features features = {};
            gf256_get_active_x86_cpu_features(&features);
            std::printf(
                "wh2_fixed5_selftest,version=1,status=pass,wide_build=%d,"
                "avx2=%d,avx512=%d,timing_workload=0\n",
                kWideXorBuild, features.AVX2, features.AVX512);
            return 0;
        }
        if (argc != 7 || std::strcmp(argv[1], "--run-cell") != 0)
        {
            PrintUsage(argv[0]);
            return 2;
        }
        uint32_t ordinal = 0u;
        uint32_t K = 0u;
        uint32_t block_bytes = 0u;
        uint32_t target_cpu = 0u;
        if (!ParseU32(argv[2], ordinal) || !ParseU32(argv[3], K) ||
            !ParseU32(argv[4], block_bytes) ||
            !ParseU32(argv[6], target_cpu))
        {
            PrintUsage(argv[0]);
            return 2;
        }
        const FrozenCell* cell = FindFrozenCell(ordinal, K, block_bytes);
        if (!cell ||
            std::strcmp(argv[5], WIREHAIR_FIXED5_SOURCE_GIT_COMMIT) != 0 ||
            target_cpu != static_cast<uint32_t>(kTargetCpu))
        {
            std::fprintf(stderr,
                "frozen command/source identity mismatch\n");
            return 2;
        }
        wirehair_wh2_bench::TargetIdentityReceiptV2 pre_target;
        std::string target_diagnostic;
        if (!CaptureExactTargetIdentity(pre_target, target_diagnostic))
        {
            std::fprintf(stderr, "pre-timing target identity failed: %s\n",
                target_diagnostic.c_str());
            return 1;
        }
        if (gf256_init() != 0) {
            std::fprintf(stderr, "gf256 initialization failed\n");
            return 1;
        }
        // Establish the production default once.  Candidate panels use a
        // nested scope and restore this state before fallback/A/A timing.
        ScopedWideXor production_wide_off(false);
        return RunFrozenCell(*cell, pre_target);
    }
    catch (const std::bad_alloc&) {
        std::fprintf(stderr, "fixed5 falsifier allocation exception\n");
        return 1;
    }
    catch (const std::exception& exception) {
        std::fprintf(stderr, "fixed5 falsifier exception: %s\n",
            exception.what());
        return 1;
    }
    catch (...) {
        std::fprintf(stderr, "fixed5 falsifier unknown exception\n");
        return 1;
    }
}
