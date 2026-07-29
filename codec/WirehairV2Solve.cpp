#include "WirehairV2Solve.h"

#include "../WirehairTools.h"
#include "../gf256.h"
#include "WirehairV2MixedBuckets.h"
#include "WirehairV2Plan.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>

#if defined(__linux__)
#include <sys/mman.h>
#include <unistd.h>
#endif

namespace wirehair_v2 {
namespace {

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
thread_local uint32_t OddPacketPeelSeedXor = 0u;
thread_local uint32_t PacketRowSeedMultiplier = 1u;
thread_local bool PacketRowSeedAvalanche = false;
thread_local MixedNullWitnessDiagnostic* MixedNullWitnessSink = nullptr;

struct PeelDegreeCdfState
{
    bool IsSet = false;
    bool Valid = true;
    std::vector<double> Cdf;
};

class PacketHookStateHasher
{
public:
    void U8(uint8_t value)
    {
        Value ^= value;
        Value *= UINT64_C(1099511628211);
    }

    void U32(uint32_t value)
    {
        for (unsigned i = 0; i < 4; ++i) {
            U8((uint8_t)(value >> (8u * i)));
        }
    }

    void U64(uint64_t value)
    {
        for (unsigned i = 0; i < 8; ++i) {
            U8((uint8_t)(value >> (8u * i)));
        }
    }

    void Double(double value)
    {
        uint64_t bits = 0u;
        static_assert(sizeof(bits) == sizeof(value),
            "double fingerprint representation");
        std::memcpy(&bits, &value, sizeof(bits));
        U64(bits);
    }

    uint64_t Finish() const
    {
        return Value != 0u ? Value : UINT64_C(0x9e3779b97f4a7c15);
    }

private:
    uint64_t Value = UINT64_C(14695981039346656037);
};

uint64_t FingerprintPeelDegreeCdf(const std::vector<double>& cdf)
{
    PacketHookStateHasher hash;
    hash.U64(UINT64_C(0x5748325043444631)); // "WH2PCDF1"
    hash.U64((uint64_t)cdf.size());
    for (double value : cdf) {
        hash.Double(value);
    }
    return hash.Finish();
}

// PeelRowParameters clamps every sampled degree to 64.  Canonicalizing every
// override onto exactly 64 bins is therefore semantics-preserving: missing
// bins have zero mass and all input mass at degrees >= 64 belongs to bin 64.
// A power-of-two size also permits a genuinely fixed six-comparison sampler.
constexpr size_t kPeelOverrideCdfSize = 64u;
static_assert(
    (kPeelOverrideCdfSize & (kPeelOverrideCdfSize - 1u)) == 0u,
    "peel override CDF must have a power-of-two size");

/// Peeling row-degree override; see SetPeelDegreesForTesting.  Invalid input
/// remains an active-invalid state so it cannot silently fall through to an
/// environment override or the shipped law.
thread_local PeelDegreeCdfState PeelDegreeOverride;
thread_local uint64_t PeelDegreeOverrideFingerprint =
    FingerprintPeelDegreeCdf(PeelDegreeOverride.Cdf);

bool IsAsciiSpace(char c)
{
    return c == ' ' || c == '\t' || c == '\r' ||
        c == '\n' || c == '\f' || c == '\v';
}

void SkipAsciiSpace(const char*& p)
{
    while (IsAsciiSpace(*p)) {
        ++p;
    }
}

bool MakePeelDegreeCdf(
    const std::vector<double>& weights,
    std::vector<double>& cdf)
{
    if (weights.empty()) {
        return false;
    }
    double total = 0.0;
    for (double weight : weights)
    {
        if (!std::isfinite(weight) || weight < 0.0) {
            return false;
        }
        total += weight;
        if (!std::isfinite(total)) {
            return false;
        }
    }
    if (!(total > 0.0)) {
        return false;
    }

    std::vector<double> parsed(kPeelOverrideCdfSize, 1.0);
    double running = 0.0;
    const size_t prefix_count = std::min(
        weights.size(), kPeelOverrideCdfSize - 1u);
    for (size_t i = 0; i < prefix_count; ++i)
    {
        // Divide the cumulative prefix once instead of summing independently
        // rounded fractions.  The latter can drift above one before a
        // trailing zero-weight bin; pinning only the final entry would then
        // make the CDF decrease and violate the binary decision tree's
        // sorted-CDF invariant.
        running += weights[i];
        const double probability = running / total;
        if (!std::isfinite(probability) ||
            probability < 0.0 || probability > 1.0)
        {
            return false;
        }
        parsed[i] = probability;
    }
    // The final interval owns every input degree at or above 64, plus all
    // remaining uniforms.  It is an exclusive upper bound because u < 1.
    parsed.back() = 1.0;
    cdf.swap(parsed);
    return true;
}

bool ParsePeelDegreeCdf(
    const char* text,
    std::vector<double>& cdf)
{
    if (!text) {
        return false;
    }
    const char* p = text;
    SkipAsciiSpace(p);
    if (*p == '\0') {
        return false;
    }

    std::vector<double> weights;
    for (;;)
    {
        errno = 0;
        char* end = nullptr;
        const double weight = std::strtod(p, &end);
        if (end == p || errno == ERANGE ||
            !std::isfinite(weight) || weight < 0.0)
        {
            return false;
        }
        weights.push_back(weight);
        p = end;
        SkipAsciiSpace(p);
        if (*p == '\0') {
            break;
        }
        if (*p != ',') {
            return false;
        }
        ++p;
        SkipAsciiSpace(p);
        if (*p == '\0') {
            return false;
        }
    }
    return MakePeelDegreeCdf(weights, cdf);
}

/**
    WIREHAIR_V2_PEEL_DEGREES, parsed once per process.

    Deliberately a fallback rather than the primary control: it is a process
    global, so every worker in a concurrent search would share one
    distribution, which is exactly what makes it useless to a search.  It
    exists for single-config harnesses (the timing protocol) that run one
    distribution per process.
*/
const PeelDegreeCdfState& EnvironmentPeelDegreeCdf()
{
    static const PeelDegreeCdfState parsed = [] {
        PeelDegreeCdfState result;
        const char* raw = std::getenv("WIREHAIR_V2_PEEL_DEGREES");
        if (raw == nullptr) {
            return result;
        }
        result.IsSet = true;
        result.Valid = ParsePeelDegreeCdf(raw, result.Cdf);
        return result;
    }();
    return parsed;
}

uint64_t EnvironmentPeelDegreeCdfFingerprint()
{
    static const uint64_t fingerprint =
        FingerprintPeelDegreeCdf(EnvironmentPeelDegreeCdf().Cdf);
    return fingerprint;
}

/// Resolve thread-local override first, environment second, stock otherwise.
/// False means an active source was malformed and the row must be refused.
bool ActivePeelDegreeCdf(const std::vector<double>*& cdf)
{
    cdf = nullptr;
    if (PeelDegreeOverride.IsSet)
    {
        if (!PeelDegreeOverride.Valid) {
            return false;
        }
        cdf = &PeelDegreeOverride.Cdf;
        return true;
    }
    const PeelDegreeCdfState& environment = EnvironmentPeelDegreeCdf();
    if (environment.IsSet)
    {
        if (!environment.Valid) {
            return false;
        }
        cdf = &environment.Cdf;
    }
    return true;
}

/**
    Replace params.PeelCount from an overridden degree distribution.

    STREAM PRESERVING.  PeelRowParameters::Initialize consumes its first PRNG
    draw for the weight and its later draws for the column-selection
    parameters, so re-seeding an identical PRNG here and taking ONE value
    recovers exactly the uniform the stock generator used.  Mapping that same
    uniform through the supplied CDF is what makes an ideal-soliton override
    reproduce the stock construction bit-identically instead of merely
    reproducing its distribution; everything downstream of the weight -- the
    column stride, the first column, the mix columns -- is untouched.

    The N/2 and kMaxPeelCount clamps are applied exactly as the stock path
    applies them, so an override cannot change the row-weight invariant.
*/
template<bool CountComparisons>
uint32_t SamplePeelDegreeOverrideImpl(
    const std::vector<double>& cdf,
    uint32_t rv,
    uint16_t peel_column_count,
    uint32_t fallback,
    uint32_t* comparison_count)
{
    if (cdf.empty() || peel_column_count == 0u) {
        return fallback;
    }
    // The stock generator compares against uint32 thresholds, so the uniform
    // is taken on the same [0,1) scale here.
    const double u = (double)rv * (1.0 / 4294967296.0);
    /*
        FIXED-DEPTH binary decision tree, not std::upper_bound or a linear
        scan, and the reason is measurement bias rather than speed.

        A linear scan exits at the SAMPLED index, so its cost depends on the
        distribution being sampled: a flat law scans further per row than a
        soliton does.  The override's overhead then differs between the arms of
        a comparison, and that difference does NOT cancel when the arms are
        re-baselined against a common reference, because it is not common.
        std::upper_bound is not fixed-depth either: libstdc++ takes seven
        comparisons near one end of a 64-entry range and six near the other.
        The canonical 64-bin tree below performs exactly six comparisons for
        every row, so what remains is a common fixed per-row offset.

        The comparison is STRICT.  A CDF entry is the first excluded uniform
        value for that degree, so an exact boundary belongs to the following
        degree.  lower_bound used >= and was off by one at every nonempty exact
        boundary; uniformly spaced self-test probes almost never landed on
        those 64 values and therefore missed it.
    */
    if (cdf.size() != kPeelOverrideCdfSize) {
        return fallback;
    }
    size_t index = 0u;
    for (size_t step = kPeelOverrideCdfSize >> 1u;
         step != 0u;
         step >>= 1u)
    {
        if (CountComparisons) {
            ++*comparison_count;
        }
        const size_t boundary = index + step - 1u;
        // Equality advances to the next interval: CDF entries are exclusive
        // upper bounds.  Written as an integer add so optimized builds use a
        // compare/select rather than a data-dependent loop exit.
        index += (size_t)(u >= cdf[boundary]) * step;
    }
    uint32_t weight = (uint32_t)index + 1u;
    const uint32_t max_weight = peel_column_count / 2u;
    if (max_weight != 0u && weight > max_weight) {
        weight = max_weight;
    }
    // Same ceiling PeelRowParameters::Initialize asserts against.
    if (weight > 64u) {
        weight = 64u;
    }
    if (weight == 0u) {
        weight = 1u;
    }
    return weight;
}

uint32_t SamplePeelDegreeOverride(
    const std::vector<double>& cdf,
    uint32_t rv,
    uint16_t peel_column_count,
    uint32_t fallback)
{
    // The false specialization compiles out the counter entirely, so the
    // test oracle cannot perturb the timing path it verifies.
    return SamplePeelDegreeOverrideImpl<false>(
        cdf, rv, peel_column_count, fallback, nullptr);
}

bool ApplyPeelDegreeOverride(
    uint32_t row_seed,
    uint32_t p_seed,
    uint16_t peel_column_count,
    wirehair::PeelRowParameters& params)
{
    const std::vector<double>* cdf = nullptr;
    if (!ActivePeelDegreeCdf(cdf)) {
        return false;
    }
    if (cdf == nullptr || peel_column_count == 0u) {
        return true;
    }
    wirehair::PCGRandom prng;
    prng.Seed(row_seed, p_seed);
    params.PeelCount = (uint16_t)SamplePeelDegreeOverride(
        *cdf, prng.Next(), peel_column_count, params.PeelCount);
    return true;
}
#endif

constexpr uint32_t PackedWordCount(uint32_t bit_count)
{
    return bit_count / 64u + ((bit_count & 63u) != 0u ? 1u : 0u);
}

static_assert(
    PackedWordCount(UINT32_MAX) == UINT32_C(67108864),
    "packed word count must not wrap at the uint32 boundary");

#if defined(GF256_TRY_TARGET_AVX2) && !defined(__AVX2__)
// DispatchPrecodeSolve enables this only after the shared x86 capability
// check.  Native AVX2 builds compile the same projection loop directly.
thread_local bool TargetWideProjectionXor = false;
#endif

struct ColumnSpan
{
    const uint32_t* First = nullptr;
    const uint32_t* Last = nullptr;

    const uint32_t* begin() const { return First; }
    const uint32_t* end() const { return Last; }
    size_t size() const { return (size_t)(Last - First); }
};

struct BinaryEquationView
{
    ColumnSpan Columns;
    const uint8_t* Data = nullptr;
};

class BinaryEquationArena
{
public:
    void Initialize(size_t row_count, size_t reference_count)
    {
        RowOffsets.resize(row_count + 1u);
        RowData.resize(row_count);
        Columns.reserve(reference_count);
        RowOffsets[0] = 0u;
    }

    void BeginRow(const uint8_t* data)
    {
        RowData[NextRow] = data;
    }

    void AppendColumn(uint32_t column)
    {
        Columns.push_back(column);
    }

    void AppendRow(
        const std::vector<uint32_t>& columns,
        const uint8_t* data)
    {
        BeginRow(data);
        Columns.insert(Columns.end(), columns.begin(), columns.end());
        EndRow();
    }

    void EndRow()
    {
        ++NextRow;
        RowOffsets[NextRow] = Columns.size();
    }

    bool IsComplete(size_t row_count, size_t reference_count) const
    {
        return NextRow == row_count && Columns.size() == reference_count;
    }

    size_t size() const { return RowData.size(); }

    uint64_t StorageBytes() const
    {
        return (uint64_t)RowOffsets.capacity() * sizeof(size_t) +
            (uint64_t)Columns.capacity() * sizeof(uint32_t) +
            (uint64_t)RowData.capacity() * sizeof(const uint8_t*);
    }

    uint32_t StorageAllocations() const
    {
        return (RowOffsets.capacity() != 0u ? 1u : 0u) +
            (Columns.capacity() != 0u ? 1u : 0u) +
            (RowData.capacity() != 0u ? 1u : 0u);
    }

    BinaryEquationView operator[](size_t row) const
    {
        BinaryEquationView view;
        view.Columns.First = Columns.data() + RowOffsets[row];
        view.Columns.Last = Columns.data() + RowOffsets[row + 1u];
        view.Data = RowData[row];
        return view;
    }

    // Projection discovers the solve column while this row is cache-hot.  Keep
    // it last so reconstruction can traverse only the dependencies without an
    // unpredictable self-column test in every sparse equation.
    void MoveSolveColumnToEnd(size_t row, size_t column_offset)
    {
        const size_t first = RowOffsets[row];
        const size_t last = RowOffsets[row + 1u];
        CAT_DEBUG_ASSERT(first + column_offset < last);
        std::swap(Columns[first + column_offset], Columns[last - 1u]);
    }

    BinaryEquationView SolveDependencies(size_t row) const
    {
        BinaryEquationView view = (*this)[row];
        CAT_DEBUG_ASSERT(view.Columns.First != view.Columns.Last);
        --view.Columns.Last;
        return view;
    }

private:
    std::vector<size_t> RowOffsets;
    std::vector<uint32_t> Columns;
    std::vector<const uint8_t*> RowData;
    size_t NextRow = 0u;
};

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
void CaptureMixedNullWitness(
    const PrecodeSystem& system,
    uint32_t column_count,
    uint32_t inactive_count,
    uint32_t projection_words,
    const BinaryEquationArena& rows,
    const std::vector<uint32_t>& inactive_index,
    const std::vector<uint32_t>& inactive_columns,
    const std::vector<uint64_t>& projection,
    const std::vector<uint64_t>& binary_pivot_coeff,
    const std::vector<uint8_t>& binary_have_pivot,
    uint32_t binary_rank,
    uint32_t expected_quotient_rank,
    MixedNullWitnessDiagnostic* sink) noexcept;
#endif

inline uint32_t PacketRowSeedForBlockId(uint32_t block_id)
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    uint32_t seed = block_id * PacketRowSeedMultiplier;
    if (PacketRowSeedAvalanche)
    {
        seed = (seed ^ (seed >> 16)) * UINT32_C(0x7feb352d);
        seed = (seed ^ (seed >> 15)) * UINT32_C(0x846ca68b);
        seed ^= seed >> 16;
    }
    return seed;
#else
    return block_id;
#endif
}

inline uint32_t PacketPeelSeedForBlockId(
    uint32_t block_id,
    const PacketRowConfig& config)
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    if ((block_id & 1u) != 0u) {
        return config.PeelSeed ^ OddPacketPeelSeedXor;
    }
#else
    (void)block_id;
#endif
    return config.PeelSeed;
}

bool InitializePacketRowParameters(
    uint32_t source_count,
    uint32_t precode_count,
    uint32_t block_id,
    const PacketRowConfig& config,
    const PacketRowRuntime& runtime,
    wirehair::PeelRowParameters& params)
{
    if (!runtime.IsValidFor(
            source_count, precode_count, config.MixCount))
    {
        return false;
    }
    const uint32_t row_seed = PacketRowSeedForBlockId(block_id);
    const uint32_t peel_seed = PacketPeelSeedForBlockId(block_id, config);
    params.Initialize(
        row_seed,
        peel_seed,
        (uint16_t)source_count,
        (uint16_t)precode_count);
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    // V2 overrides the weight AFTER V1 has produced the row, so WirehairTools
    // and the wire-format row layout stay untouched; only PeelCount moves, and
    // only when an override is installed.
    if (!ApplyPeelDegreeOverride(
            row_seed, peel_seed, (uint16_t)source_count, params))
    {
        return false;
    }
#endif
    return true;
}

template<class Prepare, class Append>
bool ForEachPacketMatrixColumn(
    uint32_t source_count,
    uint32_t precode_count,
    uint32_t block_id,
    const PacketRowConfig& config,
    const PacketRowRuntime& runtime,
    const Prepare& prepare,
    const Append& append)
{
    wirehair::PeelRowParameters params;
    if (!InitializePacketRowParameters(
            source_count, precode_count, block_id, config, runtime, params))
    {
        return false;
    }
    if (!prepare((size_t)params.PeelCount + config.MixCount)) {
        return false;
    }
    wirehair::PeelRowIterator source(
        params, (uint16_t)source_count, runtime.SourcePrime());
    do {
        append(source.GetColumn());
    } while (source.Iterate());
    if (config.MixCount == 1u)
    {
        // RowMixIterator's first output is exactly MixFirst.  Avoid producing
        // its unused second and third columns for the one-mix experiment.
        append(source_count + params.MixFirst);
        return true;
    }
    const wirehair::RowMixIterator mix(
        params, (uint16_t)precode_count, runtime.PrecodePrime());
    for (uint32_t i = 0; i < config.MixCount; ++i) {
        append(source_count + mix.Columns[i]);
    }
    return true;
}

struct PeelResult
{
    std::vector<uint32_t> SolveRow;
    std::vector<uint32_t> PeelOrder;
    std::vector<uint32_t> InactiveOrder;
    std::vector<uint8_t> UsedRows;
    uint64_t AdjacencyStorageBytes = 0u;
    uint32_t AdjacencyStorageAllocations = 0u;
    /// Deterministic peel work receipts, see PrecodeSolveStats.
    uint64_t AdjacencyVisits = 0u;
    uint64_t RowScanSteps = 0u;
    uint64_t HeapOperations = 0u;
};

struct PeelRowState
{
    // Validated WH2 systems have at most UINT16_MAX columns, so the live
    // degree and XOR of the live column ids at degree one or two fit in the
    // same four bytes previously used by the live-degree vector.  The XOR is
    // first recorded by the existing degree-two scan, then reduced to the
    // sole live column when the degree falls to one.
    uint16_t Live = 0u;
    uint16_t LowDegreeXor = 0u;
};

static_assert(
    sizeof(PeelRowState) == sizeof(uint32_t),
    "peel row state must not increase scratch storage");

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
std::atomic<uint32_t> MixedProjectionOracleUsers(0u);
std::atomic<uint64_t> MixedProjectionOracleComparisons(0u);
std::atomic<uint32_t> BinaryPeelOracleUsers(0u);
std::atomic<uint64_t> BinaryPeelOracleComparisons(0u);
#endif

bool CheckedBlockStorage(
    uint32_t block_count,
    uint32_t block_bytes,
    size_t& bytes_out)
{
    const uint64_t bytes = (uint64_t)block_count * block_bytes;
    if (block_count == 0u || block_bytes == 0u ||
        block_bytes > 0x7fffffffu ||
        bytes > (uint64_t)std::numeric_limits<size_t>::max())
    {
        return false;
    }
    bytes_out = (size_t)bytes;
    return true;
}

/**
    Backing-store capacity for a block buffer of `bytes`.

    A zero-width solve -- "perform every block operation, carry no payload" --
    leaves every block buffer empty, and an empty std::vector has a null
    data().  The solver still forms block base addresses from that pointer and
    still hands them to memcpy and the GF kernels with a zero length; both are
    undefined on a null pointer even though they move no bytes.  One byte of
    capacity keeps every base address real.  Nothing ever reads or writes that
    byte, because every block is zero bytes wide.
*/
size_t BlockStorageCapacity(size_t bytes)
{
    return bytes == 0u ? 1u : bytes;
}

/// Zero-filled block storage of `bytes`, whose data() is never null.
std::vector<uint8_t> MakeBlockStorage(size_t bytes)
{
    std::vector<uint8_t> storage;
    storage.reserve(BlockStorageCapacity(bytes));
    storage.resize(bytes, uint8_t{0});
    return storage;
}

/**
    Solve width that a payload-free solve DISPATCHES as.

    Several routines in this file pick between algorithm variants by payload
    width: batched versus direct RHS reduction, the binary quotient, the
    projected back-substitution, the fused block initializer, the residue
    bucket routes.  Every one of those thresholds is a performance crossover,
    and a zero-byte solve sits below all of them, so it used to select the
    NARROW variant of each -- a different algorithm from the one the shipping
    decoder runs at any real payload width.

    That matters because zero width is a MEASUREMENT mode, not a payload size:
    it exists so the block-operation counters can be harvested without moving
    bytes, and those counters ARE the deliverable.  Counting a variant the
    codec never executes makes the counts describe a phantom algorithm.
    Measured at K=128, mixed completion, on one fixed loss stream, varying
    only the width:

        width 0     xors 1819   muladds 1493   copies 266   zerofill 336
        width >=64  xors 1645   muladds 1493   copies 254   zerofill 331

    Mixed completion is constant from 64 up, but that is not the crossover for
    every configuration -- generic GF(256) completion still moves at 2048
    (the binary quotient), and the residue bucket routes at 1280 and 4096 --
    so "large enough" has to mean large enough for ALL of them at once.

    The dispatch width must be the width the cost model was calibrated at, not
    merely "large."  WirehairV2EsCostModel.inc is a 1280-byte model: choosing
    65536 here used different residual/binary-quotient algorithms and made the
    model multiply coefficients fitted to one operation mix by counters from
    another.  The parity harness therefore compares zero against 1280 exactly.
*/
constexpr uint32_t kCountingDispatchBlockBytes = 1280u;

/**
    Width used for DISPATCH decisions only.

    A zero-byte solve moves no payload, but must choose the same code paths
    the shipping decoder chooses, or the harvested counters describe an
    algorithm that never runs.

    Use this ONLY where the comparison selects an algorithm variant.  Never
    use it for allocation sizing, buffer capacity, memory caps, or address
    arithmetic: there the real (zero) byte count is the correct answer, and
    substituting a fictional width would size or address memory that does not
    exist.
*/
constexpr uint32_t DispatchBlockBytes(uint32_t block_bytes)
{
    return block_bytes != 0u ? block_bytes : kCountingDispatchBlockBytes;
}

/*
    Payload-plane wrappers for the GF(2^16) kernels.

    The GF(16) block API rejects a zero-length operation as invalid input --
    a deliberate, tested contract, since every caller that means to move bytes
    has bytes to move.  A zero-width solve is the one caller that legitimately
    has none: it performs the same sequence of operations on blocks that are
    zero bytes wide, so each conversion below has nothing to convert.

    These wrappers make that the explicit local rule -- skip the kernel, keep
    the caller's control flow and its operation counters -- rather than
    weakening the kernel contract for every other caller.  At any nonzero
    width they are the kernel call and nothing else.
*/
bool DeinterleavePayloadPlanes(
    const void* source, void* low, void* high, uint32_t bytes)
{
    return bytes == 0u || GF16Deinterleave(source, low, high, bytes);
}

bool InterleavePayloadPlanes(
    const void* low, const void* high, void* destination, uint32_t bytes)
{
    return bytes == 0u || GF16Interleave(low, high, destination, bytes);
}

bool MulAddPayloadPlanes(
    void* destination_low,
    void* destination_high,
    uint16_t scale,
    const void* source_low,
    const void* source_high,
    uint32_t elements)
{
    return elements == 0u ||
        GF16MulAddPlanar(
            destination_low, destination_high, scale,
            source_low, source_high, elements);
}

bool MulAddPayloadPlanes2(
    void* destination0_low,
    void* destination0_high,
    uint16_t scale0,
    void* destination1_low,
    void* destination1_high,
    uint16_t scale1,
    const void* source_low,
    const void* source_high,
    uint32_t elements)
{
    return elements == 0u ||
        GF16MulAddPlanar2(
            destination0_low, destination0_high, scale0,
            destination1_low, destination1_high, scale1,
            source_low, source_high, elements);
}

bool ScalePayloadPlanes(
    void* low,
    void* high,
    uint16_t scale,
    void* scratch,
    uint32_t elements)
{
    return elements == 0u ||
        GF16ScalePlanar(low, high, scale, scratch, elements);
}

bool MemoryRangesOverlap(
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
        // A real object cannot wrap the address space.  Treat artificial
        // wrapping pointer/length pairs as overlapping so callers fail closed
        // before any memory access.
        return true;
    }
    const uintptr_t first_end = first_begin + first_bytes;
    const uintptr_t second_end = second_begin + second_bytes;
    return first_begin < second_end && second_begin < first_end;
}

void AddScaledBlock(
    uint8_t* dst,
    uint8_t scale,
    const uint8_t* src,
    uint32_t block_bytes,
    PrecodeSolveStats& stats)
{
    if (scale == 0u) {
        return;
    }
    if (scale == 1u)
    {
        gf256_add_mem(dst, src, (int)block_bytes);
        ++stats.BlockXors;
    }
    else
    {
        gf256_muladd_mem(dst, scale, src, (int)block_bytes);
        ++stats.BlockMulAdds;
    }
}

void AddScaledBlocks(
    void* const* destinations,
    const uint8_t* scales,
    uint32_t destination_count,
    const uint8_t* source,
    uint32_t block_bytes,
    PrecodeSolveStats& stats)
{
    for (uint32_t i = 0; i < destination_count; ++i)
    {
        if (scales[i] == 1u) {
            ++stats.BlockXors;
        }
        else if (scales[i] > 1u) {
            ++stats.BlockMulAdds;
        }
    }
    gf256_muladd_multi_mem(
        destinations,
        scales,
        (int)destination_count,
        source,
        (int)block_bytes);
}

// Add one source block to several independent destinations while it remains
// hot in registers.  The packed binary residual back-elimination below often
// clears one new pivot from several existing RHS rows.  Calling gf256_add_mem
// separately rereads the same source for each destination; this loop loads it
// once per SIMD lane and fans that value out to the destination rows.
#if (defined(__AVX2__) || defined(GF256_TRY_TARGET_AVX2)) && \
    !defined(WH_COUNT)
static GF256_AVX2_TARGET void XorBlockIntoDestinationsAVX2(
    uint8_t* const* destinations,
    uint32_t destination_count,
    const uint8_t* source,
    uint32_t block_bytes)
{
    uint32_t offset = 0u;
    for (; block_bytes - offset >= 32u; offset += 32u)
    {
        const __m256i value = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(source + offset));
        for (uint32_t destination = 0u;
             destination < destination_count;
             ++destination)
        {
            __m256i* const output = reinterpret_cast<__m256i*>(
                destinations[destination] + offset);
            _mm256_storeu_si256(
                output,
                _mm256_xor_si256(_mm256_loadu_si256(output), value));
        }
    }
    if (offset < block_bytes)
    {
        const uint32_t tail = block_bytes - offset;
        for (uint32_t destination = 0u;
             destination < destination_count;
             ++destination)
        {
            gf256_add_mem(
                destinations[destination] + offset,
                source + offset,
                (int)tail);
        }
    }
}
#endif

#if defined(GF256_TRY_TARGET_AVX2) && !defined(WH_COUNT)
static const gf256_x86_cpu_features& ActiveXorCpuFeatures()
{
    static const gf256_x86_cpu_features features = []() {
        gf256_x86_cpu_features result = {};
        gf256_get_active_x86_cpu_features(&result);
        return result;
    }();
    return features;
}
#endif

static bool CanXorBlockIntoDestinationsVectorized()
{
#if defined(WH_COUNT)
    // Keep gf256's optional thread-local byte/call instrumentation exact.
    // Instrumented builds retain the individually accounted baseline calls.
    return false;
#elif defined(__AVX2__)
    return true;
#elif defined(GF256_TRY_TARGET_AVX2)
    return ActiveXorCpuFeatures().AVX2 != 0;
#else
    return false;
#endif
}

static void XorBlockIntoDestinations(
    uint8_t* const* destinations,
    uint32_t destination_count,
    const uint8_t* source,
    uint32_t block_bytes)
{
    if (destination_count == 0u) {
        return;
    }
    if (destination_count == 1u)
    {
        gf256_add_mem(destinations[0], source, (int)block_bytes);
        return;
    }
#if defined(WH_COUNT)
    for (uint32_t destination = 0u;
         destination < destination_count;
         ++destination)
    {
        gf256_add_mem(
            destinations[destination], source, (int)block_bytes);
    }
    return;
#elif defined(__AVX2__)
    XorBlockIntoDestinationsAVX2(
        destinations, destination_count, source, block_bytes);
    return;
#elif defined(GF256_TRY_TARGET_AVX2)
    // The only caller reaches this helper after the runtime capability gate.
    // Avoid repeating CPUID-derived dispatch for every four destinations.
    CAT_DEBUG_ASSERT(ActiveXorCpuFeatures().AVX2 != 0);
    XorBlockIntoDestinationsAVX2(
        destinations, destination_count, source, block_bytes);
    return;
#else
    for (uint32_t destination = 0u;
         destination < destination_count;
         ++destination)
    {
        gf256_add_mem(
            destinations[destination], source, (int)block_bytes);
    }
#endif
}

class BatchedBlockXorAccumulator
{
public:
    BatchedBlockXorAccumulator(uint8_t* destination, uint32_t block_bytes)
        : Destination(destination)
        , BlockBytes(block_bytes)
    {
    }

    void Add(const uint8_t* source)
    {
        // All callers traverse validated equation columns or residue buckets,
        // whose source blocks are distinct within one accumulator.
        PendingSources[PendingCount++] = source;
        if (PendingCount == kBatchSize) {
            Flush();
        }
    }

    void Flush()
    {
        if (PendingCount == 0u) return;
        gf256_add_multi_mem(
            Destination, PendingSources, (int)PendingCount,
            (int)BlockBytes);
        PendingCount = 0u;
    }

private:
    static const uint32_t kBatchSize = 16u;
    uint8_t* Destination;
    uint32_t BlockBytes;
    const void* PendingSources[kBatchSize];
    uint32_t PendingCount = 0u;
};

// Initializes the destination from the XOR of its sources.  The first batch
// uses the set-form SIMD kernel so a packet payload and its dependent values
// are consumed in one pass rather than memcpy followed by a read/modify/write.
class BatchedBlockXorInitializer
{
public:
    BatchedBlockXorInitializer(
        uint8_t* destination,
        uint32_t block_bytes,
        const uint8_t* first_source,
        bool destination_initially_zero = false)
        : Destination(destination)
        , BlockBytes(block_bytes)
        , DestinationInitiallyZero(destination_initially_zero)
        , BatchCapacity(kFusedBatchSize)
    {
        if (first_source) {
            PendingSources[PendingCount++] = first_source;
        }
    }

    void Add(const uint8_t* source)
    {
        PendingSources[PendingCount++] = source;
        if (PendingCount == BatchCapacity) {
            Flush();
        }
    }

    void Flush()
    {
        if (!Initialized)
        {
            if (PendingCount == 0u) {
                if (!DestinationInitiallyZero) {
                    std::memset(Destination, 0, BlockBytes);
                }
            }
            else {
                gf256_addset_multi_mem(
                    Destination, PendingSources, (int)PendingCount,
                    (int)BlockBytes);
            }
            Initialized = true;
            BatchCapacity = kRegularBatchSize;
        }
        else if (PendingCount != 0u)
        {
            gf256_add_multi_mem(
                Destination, PendingSources, (int)PendingCount,
                (int)BlockBytes);
        }
        PendingCount = 0u;
    }

private:
    static const uint32_t kRegularBatchSize = 8u;
    static const uint32_t kFusedBatchSize = 16u;
    uint8_t* Destination;
    uint32_t BlockBytes;
    bool DestinationInitiallyZero;
    const void* PendingSources[kFusedBatchSize];
    uint32_t PendingCount = 0u;
    uint32_t BatchCapacity;
    bool Initialized = false;
};

// Reconstruct one binary pivot directly from its affine RHS and the solved
// free variables.  Shipping mixed completion has H == 12, so the complete
// source set (pivot RHS plus every selected free value) fits in the
// initializer's 16-source set-form batch; larger test-hook geometries safely
// spill into its additive batch.  BlockXors remains a logical equation-work
// receipt: the affine RHS initializes the value and only the selected relation
// terms count as XORs.
void InitializeMixedBinaryPivotValue(
    uint8_t* value,
    uint32_t block_bytes,
    const uint8_t* pivot_rhs,
    const uint64_t* relation,
    const std::vector<uint32_t>& free_columns,
    const std::vector<uint32_t>& inactive_columns,
    const std::vector<uint8_t>& values,
    PrecodeSolveStats& stats)
{
    BatchedBlockXorInitializer value_xor(
        value, block_bytes, pivot_rhs);
    for (uint32_t free_column : free_columns)
    {
        if ((relation[free_column >> 6] &
                (UINT64_C(1) << (free_column & 63u))) == 0u)
        {
            continue;
        }
        value_xor.Add(
            values.data() +
                (size_t)inactive_columns[free_column] * block_bytes);
        ++stats.BlockXors;
    }
    value_xor.Flush();
}

template<uint32_t Count>
struct ProjectionSourceBatch
{
    static GF256_FORCE_INLINE uint64_t Xor64(
        uint64_t value,
        const uint64_t* const* GF256_RESTRICT sources,
        uint32_t word)
    {
        return ProjectionSourceBatch<Count - 1u>::Xor64(
            value, sources, word) ^ sources[Count - 1u][word];
    }

#if defined(GF256_TARGET_X86_SIMD)
    static GF256_FORCE_INLINE __m128i Xor128(
        __m128i value,
        const uint64_t* const* GF256_RESTRICT sources,
        uint32_t word)
    {
        return _mm_xor_si128(
            ProjectionSourceBatch<Count - 1u>::Xor128(
                value, sources, word),
            _mm_loadu_si128(reinterpret_cast<const __m128i*>(
                sources[Count - 1u] + word)));
    }
#endif

#if defined(__AVX2__) || defined(GF256_TRY_TARGET_AVX2)
    // Projection rows are short (typically a handful of packed words), so a
    // single unrolled YMM chain avoids half of the load/XOR instructions
    // without the frequency cost measured for an AVX-512 version.
    static GF256_AVX2_TARGET GF256_FORCE_INLINE __m256i Xor256(
        __m256i value,
        const uint64_t* const* GF256_RESTRICT sources,
        uint32_t word)
    {
        return _mm256_xor_si256(
            ProjectionSourceBatch<Count - 1u>::Xor256(
                value, sources, word),
            _mm256_loadu_si256(reinterpret_cast<const __m256i*>(
                sources[Count - 1u] + word)));
    }
#endif

};

template<>
struct ProjectionSourceBatch<0u>
{
    static GF256_FORCE_INLINE uint64_t Xor64(
        uint64_t value,
        const uint64_t* const* GF256_RESTRICT,
        uint32_t)
    {
        return value;
    }

#if defined(GF256_TARGET_X86_SIMD)
    static GF256_FORCE_INLINE __m128i Xor128(
        __m128i value,
        const uint64_t* const* GF256_RESTRICT,
        uint32_t)
    {
        return value;
    }
#endif

#if defined(__AVX2__) || defined(GF256_TRY_TARGET_AVX2)
    static GF256_AVX2_TARGET GF256_FORCE_INLINE __m256i Xor256(
        __m256i value,
        const uint64_t* const* GF256_RESTRICT,
        uint32_t)
    {
        return value;
    }
#endif

};

#if defined(GF256_TRY_TARGET_AVX2) && !defined(__AVX2__)
template<uint32_t Count>
static GF256_AVX2_TARGET uint32_t XorProjectionSourceBatchAVX2(
    uint64_t* GF256_RESTRICT destination,
    const uint64_t* const* GF256_RESTRICT sources,
    uint32_t words)
{
    uint32_t word = 0u;
    for (; words - word >= 4u; word += 4u)
    {
        const __m256i destination0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(destination + word));
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(destination + word),
            ProjectionSourceBatch<Count>::Xor256(
                destination0, sources, word));
    }
    return word;
}
#endif

template<uint32_t Count>
static GF256_FORCE_INLINE void XorProjectionSourceBatch(
    uint64_t* GF256_RESTRICT destination,
    const uint64_t* const* GF256_RESTRICT sources,
    uint32_t words)
{
    uint32_t word = 0u;
#if defined(__AVX2__)
    for (; words - word >= 4u; word += 4u)
    {
        const __m256i destination0 = _mm256_loadu_si256(
            reinterpret_cast<const __m256i*>(destination + word));
        _mm256_storeu_si256(
            reinterpret_cast<__m256i*>(destination + word),
            ProjectionSourceBatch<Count>::Xor256(
                destination0, sources, word));
    }
#elif defined(GF256_TRY_TARGET_AVX2)
    if (TargetWideProjectionXor) {
        word = XorProjectionSourceBatchAVX2<Count>(
            destination, sources, words);
    }
#endif
#if defined(GF256_TARGET_X86_SIMD)
    for (; words - word >= 4u; word += 4u)
    {
        const __m128i destination0 = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(destination + word));
        const __m128i destination1 = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(destination + word + 2u));
        _mm_storeu_si128(
            reinterpret_cast<__m128i*>(destination + word),
            ProjectionSourceBatch<Count>::Xor128(
                destination0, sources, word));
        _mm_storeu_si128(
            reinterpret_cast<__m128i*>(destination + word + 2u),
            ProjectionSourceBatch<Count>::Xor128(
                destination1, sources, word + 2u));
    }
    if (words - word >= 2u)
    {
        const __m128i destination0 = _mm_loadu_si128(
            reinterpret_cast<const __m128i*>(destination + word));
        _mm_storeu_si128(
            reinterpret_cast<__m128i*>(destination + word),
            ProjectionSourceBatch<Count>::Xor128(
                destination0, sources, word));
        word += 2u;
    }
#endif
    for (; word < words; ++word) {
        destination[word] = ProjectionSourceBatch<Count>::Xor64(
            destination[word], sources, word);
    }
}

class BatchedProjectionXorAccumulator
{
public:
    BatchedProjectionXorAccumulator(
        uint64_t* destination,
        const uint64_t* source_base,
        uint32_t words)
        : Destination(destination), SourceBase(source_base), Words(words)
    {
    }

    GF256_FORCE_INLINE void Add(uint32_t source_index)
    {
        if (Words == 0u) {
            return;
        }
        Sources[Count++] =
            SourceBase + (size_t)source_index * Words;
        if (Count == kBatchSize) {
            Flush();
        }
    }

    GF256_FORCE_INLINE void Flush()
    {
        switch (Count)
        {
#define WIREHAIR_PROJECTION_BATCH_CASE(n) \
        case n: XorProjectionSourceBatch<n>( \
            Destination, Sources, Words); break
        WIREHAIR_PROJECTION_BATCH_CASE(1u);
        WIREHAIR_PROJECTION_BATCH_CASE(2u);
        WIREHAIR_PROJECTION_BATCH_CASE(3u);
        WIREHAIR_PROJECTION_BATCH_CASE(4u);
        WIREHAIR_PROJECTION_BATCH_CASE(5u);
        WIREHAIR_PROJECTION_BATCH_CASE(6u);
        WIREHAIR_PROJECTION_BATCH_CASE(7u);
        WIREHAIR_PROJECTION_BATCH_CASE(8u);
        WIREHAIR_PROJECTION_BATCH_CASE(9u);
        WIREHAIR_PROJECTION_BATCH_CASE(10u);
        WIREHAIR_PROJECTION_BATCH_CASE(11u);
        WIREHAIR_PROJECTION_BATCH_CASE(12u);
        WIREHAIR_PROJECTION_BATCH_CASE(13u);
        WIREHAIR_PROJECTION_BATCH_CASE(14u);
        WIREHAIR_PROJECTION_BATCH_CASE(15u);
        WIREHAIR_PROJECTION_BATCH_CASE(16u);
#undef WIREHAIR_PROJECTION_BATCH_CASE
        default: break;
        }
        Count = 0u;
    }

private:
    static const uint32_t kBatchSize = 16u;
    uint64_t* Destination;
    const uint64_t* SourceBase;
    uint32_t Words;
    const uint64_t* Sources[kBatchSize];
    uint32_t Count = 0u;
};

template<class XorAccumulator>
static GF256_FORCE_INLINE uint32_t AccumulatePeeledProjectionConstant(
    uint32_t column,
    const BinaryEquationView& equation,
    const std::vector<uint32_t>& inactive_index,
    uint32_t words,
    const std::vector<uint64_t>& projection,
    const std::vector<uint8_t>& values,
    uint32_t block_bytes,
    std::vector<uint64_t>& accumulator,
    XorAccumulator& constant_xor,
    PrecodeSolveStats& stats)
{
    BatchedProjectionXorAccumulator projection_xor(
        accumulator.data(), projection.data(), words);
    uint32_t solve_column_offset = UINT32_MAX;
    for (const uint32_t* current = equation.Columns.begin();
         current != equation.Columns.end();
         ++current)
    {
        const uint32_t other = *current;
        if (other == column) {
            solve_column_offset = (uint32_t)(
                current - equation.Columns.begin());
            continue;
        }
        const uint32_t index = inactive_index[other];
        if (index != UINT32_MAX) {
            accumulator[index >> 6] ^=
                UINT64_C(1) << (index & 63u);
        }
        else
        {
            projection_xor.Add(other);
            // Inactive value slots are still the zero constant at this
            // stage.  Only peeled columns can contribute to the affine RHS;
            // XORing an inactive slot would have no algebraic effect.
            constant_xor.Add(
                values.data() + (size_t)other * block_bytes);
            ++stats.BlockXors;
            stats.ProjectionWordXors += words;
        }
    }
    projection_xor.Flush();
    return solve_column_offset;
}

bool RowIsZero(const uint8_t* data, uint32_t bytes);

enum class ResidualInsertResult
{
    Dependent,
    Inserted,
    Inconsistent,
    Independent
};

constexpr uint32_t kResidualCoefficientBulkThreshold = 16u;
// Mixed-solver measurements put the multi-source RHS crossover at a 4-KiB
// payload; the extra wide-kernel setup is neutral or slower at MTU sizes.
constexpr uint32_t kBatchedResidualRhsMinBlockBytes = 4096u;
static_assert(
    kCountingDispatchBlockBytes < kBatchedResidualRhsMinBlockBytes,
    "count-only dispatch must match the 1280-byte unbatched-RHS model");
// CheckedBlockStorage caps valid payloads below this sentinel.
constexpr uint32_t kNeverBatchResidualRhs = UINT32_MAX;
// The sentinel must stay above the dispatch width so "never batch" keeps
// meaning never, including for a payload-free solve.
static_assert(
    kNeverBatchResidualRhs > kCountingDispatchBlockBytes,
    "the never-batch sentinel must outrank the counting dispatch width");

// Keep the 4-KiB path out of line and in the compiler's cold section.
// Placing it beside the literal loop reproducibly regressed 1280-byte solves,
// while 4-KiB payloads amortize the size-optimized helper and wide XOR setup.
#if defined(_MSC_VER)
#define WH2_RESIDUAL_COLD_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define WH2_RESIDUAL_COLD_NOINLINE __attribute__((noinline, cold))
#else
#define WH2_RESIDUAL_COLD_NOINLINE
#endif

static WH2_RESIDUAL_COLD_NOINLINE void ReduceResidualRowWithBatchedRhs(
    std::vector<uint8_t>& coeff,
    std::vector<uint8_t>& rhs,
    uint32_t R,
    uint32_t block_bytes,
    const std::vector<uint8_t>& pivot_coeff,
    const std::vector<uint8_t>& pivot_rhs,
    const std::vector<uint8_t>& have_pivot,
    PrecodeSolveStats& stats);

#undef WH2_RESIDUAL_COLD_NOINLINE

#if defined(_MSC_VER)
#define WH2_RESIDUAL_NOINLINE __declspec(noinline)
// GNU ld and lld collect .text.* in their normal executable text segment.
// Bespoke ELF scripts that list only an exact .text input section can disable
// the placement while retaining noinline; see codec/README.md.
#elif defined(__ELF__) && (defined(__GNUC__) || defined(__clang__)) && \
    !defined(WIREHAIR_V2_DISABLE_PACKED_RESIDUAL_TEXT_SECTION)
#define WH2_RESIDUAL_NOINLINE \
    __attribute__((noinline, section(".text.wh2_packed_residual")))
#elif defined(__GNUC__) || defined(__clang__)
#define WH2_RESIDUAL_NOINLINE __attribute__((noinline))
#else
#define WH2_RESIDUAL_NOINLINE
#endif

static WH2_RESIDUAL_NOINLINE ResidualInsertResult
InsertPackedBinaryResidualRow(
    std::vector<uint64_t>& coeff,
    std::vector<uint8_t>& rhs,
    uint32_t R,
    uint32_t words,
    uint32_t block_bytes,
    std::vector<uint64_t>& pivot_coeff,
    std::vector<uint8_t>& pivot_rhs,
    std::vector<uint8_t>& have_pivot,
    uint32_t& rank,
    PrecodeSolveStats& stats,
    uint32_t batched_rhs_min_block_bytes);

#undef WH2_RESIDUAL_NOINLINE
constexpr uint32_t kProjectedBackSubMinBlockBytes = 64u;
static_assert(
    kCountingDispatchBlockBytes >= kProjectedBackSubMinBlockBytes,
    "a payload-free solve must dispatch into projected back-substitution");
// Paired whole-solver runs show the fused path loses below these scales even
// though the isolated payload kernel is faster.
constexpr uint32_t kFusedBlockXorInitMinBlockBytes = 1280u;
static_assert(
    kCountingDispatchBlockBytes >= kFusedBlockXorInitMinBlockBytes,
    "a payload-free solve must dispatch into the fused block initializer");
constexpr uint32_t kFusedBlockXorInitMinBlockCount = 10000u;
// Remaining payload-width dispatch thresholds in this file, asserted here so
// the counting dispatch stays in the exact 1280-byte regime.  It is above the
// projected/fused/dual-bucket gates and below the binary-quotient/batched-RHS
// gates; changing any of those relationships requires a new cost calibration.
static_assert(
    kCountingDispatchBlockBytes < kBinaryQuotientMinBlockBytes,
    "count-only dispatch must match the 1280-byte pre-quotient model");
static_assert(
    kCountingDispatchBlockBytes >= 512u,
    "a payload-free solve must dispatch into the wide/batched block routes");
static_assert(
    kCountingDispatchBlockBytes >= 1024u,
    "a payload-free solve must dispatch into the dual residue bucket route");
// The tiny mixed fast path is a MAXIMUM gate: the dispatch width must stay
// clear of it so a payload-free solve keeps running the production algorithm
// while genuinely tiny payloads keep the fast path.
static_assert(
    kCountingDispatchBlockBytes >
        WIREHAIR_V2_TINY_MIXED_FASTPATH_MAX_BLOCK_BYTES,
    "the counting dispatch width must not select the tiny fast path");
static_assert(
    kCountingDispatchBlockBytes == 1280u,
    "changing the count-only dispatch width requires refitting the cost model");

// The production GF(256) profile fixes H=12, so its periodic coefficient
// table is immutable across every message.  Keeping that small table in
// process-local read-only storage avoids rebuilding or allocating it for each
// solve.  Other heavy-row profiles retain the on-demand baseline path below.
constexpr uint32_t kCachedPeriodicHeavyRows = 12u;
constexpr uint32_t kCachedPeriodicWindow =
    256u - kCachedPeriodicHeavyRows;
constexpr uint32_t kCachedPeriodicWords =
    (kCachedPeriodicHeavyRows + 7u) / 8u;

const std::array<
    uint64_t,
    kCachedPeriodicWindow * kCachedPeriodicWords>&
CachedPeriodicHeavyTable()
{
    static const std::array<
        uint64_t,
        kCachedPeriodicWindow * kCachedPeriodicWords> table = []() {
            std::array<
                uint64_t,
                kCachedPeriodicWindow * kCachedPeriodicWords> result{};
            for (uint32_t residue = 0;
                 residue < kCachedPeriodicWindow;
                 ++residue)
            {
                uint64_t* packed = result.data() +
                    (size_t)residue * kCachedPeriodicWords;
                for (uint32_t heavy = 0;
                     heavy < kCachedPeriodicHeavyRows;
                     ++heavy)
                {
                    packed[heavy >> 3] |=
                        (uint64_t)HeavyCoefficient(
                            heavy, residue, kCachedPeriodicHeavyRows) <<
                        ((heavy & 7u) * 8u);
                }
            }
            return result;
        }();
    return table;
}

void AddScaledResidualCoefficients(
    uint8_t* destination,
    uint8_t scale,
    const uint8_t* source,
    uint32_t count)
{
    if (count >= kResidualCoefficientBulkThreshold) {
        gf256_muladd_mem(destination, scale, source, (int)count);
        return;
    }
    for (uint32_t i = 0; i < count; ++i) {
        destination[i] ^= gf256_mul(source[i], scale);
    }
}

void ScaleResidualCoefficients(
    uint8_t* coefficients,
    uint8_t scale,
    uint32_t count)
{
    if (count >= kResidualCoefficientBulkThreshold) {
        gf256_mul_mem(coefficients, coefficients, scale, (int)count);
        return;
    }
    for (uint32_t i = 0; i < count; ++i) {
        coefficients[i] = gf256_mul(coefficients[i], scale);
    }
}

ResidualInsertResult InsertResidualRow(
    std::vector<uint8_t>& coeff,
    std::vector<uint8_t>& rhs,
    uint32_t R,
    uint32_t block_bytes,
    std::vector<uint8_t>& pivot_coeff,
    std::vector<uint8_t>& pivot_rhs,
    std::vector<uint8_t>& have_pivot,
    uint32_t& rank,
    PrecodeSolveStats& stats,
    bool allow_insert,
    uint32_t batched_rhs_min_block_bytes)
{
    // DISPATCH, not sizing: the branch picks a reduction algorithm.  Every
    // address and length below still uses the real block_bytes.
    if (DispatchBlockBytes(block_bytes) < batched_rhs_min_block_bytes)
    {
        for (uint32_t j = 0; j < R; ++j)
        {
            const uint8_t scale = coeff[j];
            if (scale == 0u || !have_pivot[j]) {
                continue;
            }
            const uint8_t* pivot =
                pivot_coeff.data() + (size_t)j * R;
            if (scale == 1u) {
                gf256_add_mem(
                    coeff.data() + j, pivot + j, (int)(R - j));
            }
            else {
                AddScaledResidualCoefficients(
                    coeff.data() + j, scale, pivot + j, R - j);
            }
            stats.ResidualCoeffByteOps += R - j;
            AddScaledBlock(
                rhs.data(), scale,
                pivot_rhs.data() + (size_t)j * block_bytes,
                block_bytes, stats);
        }
    }
    else
    {
        ReduceResidualRowWithBatchedRhs(
            coeff, rhs, R, block_bytes,
            pivot_coeff, pivot_rhs, have_pivot, stats);
    }

    uint32_t pivot_column = R;
    for (uint32_t j = 0; j < R; ++j) {
        if (coeff[j] != 0u) {
            pivot_column = j;
            break;
        }
    }
    if (pivot_column == R) {
        return RowIsZero(rhs.data(), block_bytes) ?
            ResidualInsertResult::Dependent :
            ResidualInsertResult::Inconsistent;
    }
    if (!allow_insert) {
        return ResidualInsertResult::Independent;
    }
    if (have_pivot[pivot_column]) {
        return ResidualInsertResult::Inconsistent;
    }

    const uint8_t pivot_value = coeff[pivot_column];
    if (pivot_value != 1u)
    {
        const uint8_t inverse = gf256_inv(pivot_value);
        ScaleResidualCoefficients(
            coeff.data() + pivot_column, inverse, R - pivot_column);
        stats.ResidualCoeffByteOps += R - pivot_column;
        gf256_div_mem(
            rhs.data(), rhs.data(), pivot_value, (int)block_bytes);
        ++stats.BlockMulAdds;
    }

    for (uint32_t existing = 0; existing < R; ++existing)
    {
        if (!have_pivot[existing]) {
            continue;
        }
        uint8_t* existing_coeff =
            pivot_coeff.data() + (size_t)existing * R;
        const uint8_t scale = existing_coeff[pivot_column];
        if (scale == 0u) {
            continue;
        }
        if (scale == 1u) {
            gf256_add_mem(
                existing_coeff + pivot_column,
                coeff.data() + pivot_column,
                (int)(R - pivot_column));
        }
        else {
            AddScaledResidualCoefficients(
                existing_coeff + pivot_column, scale,
                coeff.data() + pivot_column, R - pivot_column);
        }
        stats.ResidualCoeffByteOps += R - pivot_column;
        AddScaledBlock(
            pivot_rhs.data() + (size_t)existing * block_bytes,
            scale, rhs.data(), block_bytes, stats);
    }

    std::memcpy(
        pivot_coeff.data() + (size_t)pivot_column * R,
        coeff.data(), R);
    std::memcpy(
        pivot_rhs.data() + (size_t)pivot_column * block_bytes,
        rhs.data(), block_bytes);
    ++stats.BlockCopies;
    have_pivot[pivot_column] = 1u;
    ++rank;
    return ResidualInsertResult::Inserted;
}

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
static bool CacheMixedGroupedGF256ResidueBlockShifts(
    uint32_t column_limit,
    uint32_t coefficient_period,
    std::vector<uint8_t>& block_shifts)
{
    static_assert(
        kMixedCoefficientPeriod <= UINT8_MAX,
        "grouped residue shifts must fit their compact cache");
    if (coefficient_period == 0u ||
        coefficient_period > kMixedCoefficientPeriod)
    {
        return false;
    }
    const uint32_t block_count = column_limit / coefficient_period +
        (column_limit % coefficient_period != 0u ? 1u : 0u);
    static const uint32_t kMixedCompletionRows =
        kMixedGF256Rows + kMixedGF16Rows;
    static const uint32_t kMaxNonCornerColumns =
        UINT16_MAX - kMixedCompletionRows;
    static const uint32_t kMinGroupedPeriod = kMixedCompletionRows + 1u;
    static const uint32_t kMaxGroupedShiftBlocks =
        (kMaxNonCornerColumns + kMinGroupedPeriod - 1u) /
            kMinGroupedPeriod;
    if (block_count > kMaxGroupedShiftBlocks ||
        block_count > block_shifts.max_size())
    {
        return false;
    }
    block_shifts.resize(block_count);
    for (uint32_t block = 0u; block < block_count; ++block)
    {
        const uint32_t shift =
            ActiveMixedGroupedGF256ResidueBlockShift(block);
        if (shift >= coefficient_period) {
            return false;
        }
        block_shifts[block] = (uint8_t)shift;
    }
    return true;
}
#endif

bool ProjectMixedCompletionCoefficientsByResidueBuckets(
    uint32_t column_count,
    uint32_t inactive_count,
    uint32_t projection_words,
    const std::vector<uint32_t>& inactive_index,
    const std::vector<uint64_t>& projection,
    const MixedPackedCoefficients& cached_packed,
    const uint8_t* grouped_block_shifts,
    uint32_t grouped_block_shift_count,
    std::vector<uint64_t>& projected)
{
    const uint32_t packed_words = ActiveMixedPackedCoefficientWords();
    static_assert(
        kMixedPackedCoefficientWords <= 4u,
        "mixed completion packing must fit the unrolled projection");
    const uint32_t expected_projection_words =
        PackedWordCount(inactive_count);
    const uint32_t coefficient_period = ActiveMixedCoefficientPeriod();
    const uint32_t subfield_rows = ActiveMixedGF256Rows();
    const uint32_t extension_rows = ActiveMixedGF16Rows();
    const uint32_t grouped_gf256_rows = ActiveMixedGroupedGF256Rows();
    const uint32_t H = subfield_rows + extension_rows;
    const uint32_t populated_residues =
        std::min(coefficient_period, column_count);
    const uint64_t projection_elements =
        (uint64_t)column_count * projection_words;
    const uint64_t projected_elements =
        (uint64_t)inactive_count * packed_words;
    if (coefficient_period < H || column_count < H ||
        coefficient_period > kMixedCoefficientPeriod ||
        grouped_gf256_rows > subfield_rows ||
        inactive_index.size() != column_count ||
        projection_words != expected_projection_words ||
        projection_words == 0u ||
        projection_elements > projection.max_size() ||
        projection.size() != (size_t)projection_elements ||
        projected_elements > projected.max_size())
    {
        return false;
    }
    const uint32_t first_grouped_gf256_row =
        subfield_rows - grouped_gf256_rows;
    const uint32_t first_heavy_column = column_count - H;
    if (grouped_gf256_rows != 0u)
    {
        const uint32_t grouped_block_count =
            first_heavy_column / coefficient_period +
            (first_heavy_column % coefficient_period != 0u ? 1u : 0u);
        if (grouped_block_shift_count != grouped_block_count ||
            (grouped_block_count != 0u && !grouped_block_shifts))
        {
            return false;
        }
    }

    // C is cached once per complete or partial non-corner P-column block.
    // The final H columns deliberately keep their A residue even when they
    // share a physical block with the last non-corner columns.
    const auto grouped_residue_at = [&](
        uint32_t column,
        uint32_t primary_residue)
    {
        if (column >= first_heavy_column) {
            return primary_residue;
        }
        uint32_t residue = column % coefficient_period;
        residue += grouped_block_shifts[column / coefficient_period];
        return residue < coefficient_period ?
            residue : residue - coefficient_period;
    };

    projected.assign((size_t)projected_elements, uint64_t{0});
    const auto xor_projected = [&](uint32_t index,
                                   const uint64_t* coefficients) {
        uint64_t* destination =
            projected.data() + (size_t)index * packed_words;
        destination[0] ^= coefficients[0];
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
        // Production is exactly three words.  Test configurations span one
        // through four as H8/H16 are swept, so guard every optional lane.
        if (packed_words > 1u) destination[1] ^= coefficients[1];
        if (packed_words > 2u) destination[2] ^= coefficients[2];
        if (packed_words == 4u) destination[3] ^= coefficients[3];
#else
        destination[1] ^= coefficients[1];
        destination[2] ^= coefficients[2];
#endif
    };

    const bool independent_extension_residues =
        ActiveMixedIndependentExtensionResidues();
    const bool secondary_schedule =
        independent_extension_residues || grouped_gf256_rows != 0u;
    if (independent_extension_residues && grouped_gf256_rows != 0u) {
        return false;
    }
    uint64_t primary_masks[kMixedPackedCoefficientWords] = {};
    uint64_t secondary_masks[kMixedPackedCoefficientWords] = {};
    if (secondary_schedule)
    {
        for (uint32_t row = 0u; row < subfield_rows; ++row)
        {
            uint64_t* const masks =
                grouped_gf256_rows != 0u &&
                row >= first_grouped_gf256_row ?
                    secondary_masks : primary_masks;
            masks[row >> 2] |=
                UINT64_C(0xffff) << ((row & 3u) * 16u);
        }
        for (uint32_t er = 0u; er < extension_rows; ++er)
        {
            const uint32_t row = subfield_rows + er;
            uint64_t* const masks = independent_extension_residues ?
                secondary_masks : primary_masks;
            masks[row >> 2] |=
                UINT64_C(0xffff) << ((row & 3u) * 16u);
        }
    }

    // At or below one complete coefficient period, no two columns share a
    // coefficient vector.  Retain the direct expansion so tiny systems pay no
    // bucket setup cost.
    if (column_count <= coefficient_period)
    {
        for (uint32_t column = 0; column < column_count; ++column)
        {
            uint64_t combined_coefficients[
                kMixedPackedCoefficientWords] = {};
            const uint32_t primary_residue =
                ActiveMixedCoefficientResidue(column);
            const uint64_t* coefficients = cached_packed.ByResidue[
                primary_residue];
            // Preserve the established one-period independent-extension
            // solve path exactly.  Grouped GF(256) still composes A/C masks
            // explicitly even though both schedules currently start at zero;
            // that keeps the row partition and corner override invariant
            // visible instead of relying on today's block-zero coincidence.
            if (grouped_gf256_rows != 0u)
            {
                const uint32_t secondary_residue =
                    grouped_residue_at(column, primary_residue);
                const uint64_t* secondary_coefficients =
                    cached_packed.ByResidue[secondary_residue];
                for (uint32_t word = 0u; word < packed_words; ++word) {
                    combined_coefficients[word] =
                        (coefficients[word] & primary_masks[word]) |
                        (secondary_coefficients[word] &
                            secondary_masks[word]);
                }
                coefficients = combined_coefficients;
            }
            const uint32_t inactive = inactive_index[column];
            if (inactive != UINT32_MAX)
            {
                if (inactive >= inactive_count) {
                    return false;
                }
                xor_projected(inactive, coefficients);
                continue;
            }
            const uint64_t* bits =
                projection.data() + (size_t)column * projection_words;
            for (uint32_t word_index = 0;
                 word_index < projection_words; ++word_index)
            {
                uint64_t word = bits[word_index];
                while (word != 0u)
                {
                    const uint32_t bit =
                        wirehair::NonzeroLowestBitIndex64(word);
                    const uint32_t index = (word_index << 6) + bit;
                    if (index < inactive_count) {
                        xor_projected(index, coefficients);
                    }
                    word &= word - 1u;
                }
            }
        }
        return true;
    }

    // Mixed completion coefficients are indexed by a period-sized residue
    // bucket.  The experiment-only skew rotates the labels of each complete
    // period block without changing its balance.  Transpose the projection by
    // the active bucket before expanding any bits: each of the first
    // `coefficient_period` scratch rows become the parity buckets for each
    // residue.  Keep the original affine projections intact: back
    // substitution can then choose between each sparse equation and its
    // projected inactive-variable relation.  This is algebraically identical
    // to dense per-column expansion, while the scratch is only P*ceil(R/64)
    // words rather than L*ceil(R/64).
    std::vector<uint64_t> residue_projection(
        (size_t)coefficient_period * projection_words, uint64_t{0});
    std::vector<uint64_t> secondary_residue_projection(
        secondary_schedule ?
            (size_t)coefficient_period * projection_words : 0u,
        uint64_t{0});

    const bool rotate_residues = ActiveMixedResiduesRotated();
    const auto accumulate_column = [&](uint32_t column, uint32_t residue) {
        uint64_t* bucket = residue_projection.data() +
            (size_t)residue * projection_words;
        uint64_t* secondary_bucket = nullptr;
        if (secondary_schedule) {
            const uint32_t secondary_residue =
                independent_extension_residues ?
                    ActiveMixedExtensionCoefficientResidue(column) :
                    grouped_residue_at(column, residue);
            secondary_bucket = secondary_residue_projection.data() +
                (size_t)secondary_residue * projection_words;
        }
        const uint32_t inactive = inactive_index[column];
        if (inactive != UINT32_MAX)
        {
            if (inactive >= inactive_count) {
                return false;
            }
            bucket[inactive >> 6] ^=
                UINT64_C(1) << (inactive & 63u);
            if (secondary_bucket) {
                secondary_bucket[inactive >> 6] ^=
                    UINT64_C(1) << (inactive & 63u);
            }
        }
        else
        {
            const uint64_t* bits =
                projection.data() + (size_t)column * projection_words;
            for (uint32_t word = 0; word < projection_words; ++word) {
                bucket[word] ^= bits[word];
                if (secondary_bucket) {
                    secondary_bucket[word] ^= bits[word];
                }
            }
        }
        return true;
    };
    for (uint32_t column = 0u;
         column < populated_residues; ++column)
    {
        if (!accumulate_column(column, column)) {
            return false;
        }
    }
    if (!rotate_residues)
    {
        uint32_t residue = 0u;
        for (uint32_t column = coefficient_period;
             column < column_count; ++column)
        {
            if (!accumulate_column(column, residue)) {
                return false;
            }
            if (++residue == coefficient_period) {
                residue = 0u;
            }
        }
    }
    else
    {
        uint32_t block_index = 1u;
        uint32_t residue = ActiveMixedResidueBlockShift(block_index);
        uint32_t block_column = 0u;
        for (uint32_t column = coefficient_period;
             column < column_count; ++column)
        {
            if (!accumulate_column(column, residue)) {
                return false;
            }
            if (++block_column == coefficient_period)
            {
                block_column = 0u;
                residue = ActiveMixedResidueBlockShift(++block_index);
            }
            else if (++residue == coefficient_period) {
                residue = 0u;
            }
        }
    }

    for (uint32_t residue = 0u; residue < populated_residues; ++residue)
    {
        const uint64_t* coefficients = cached_packed.ByResidue[residue];
        const uint64_t* bucket = residue_projection.data() +
            (size_t)residue * projection_words;
        for (uint32_t word_index = 0;
             word_index < projection_words; ++word_index)
        {
            uint64_t word = bucket[word_index];
            while (word != 0u)
            {
                const uint32_t bit =
                    wirehair::NonzeroLowestBitIndex64(word);
                const uint32_t index = (word_index << 6) + bit;
                if (index < inactive_count)
                {
                    if (!secondary_schedule) {
                        xor_projected(index, coefficients);
                    }
                    else
                    {
                        uint64_t* destination = projected.data() +
                            (size_t)index * packed_words;
                        for (uint32_t packed = 0u;
                             packed < packed_words; ++packed)
                        {
                            destination[packed] ^=
                                coefficients[packed] & primary_masks[packed];
                        }
                    }
                }
                word &= word - 1u;
            }
        }
    }
    if (secondary_schedule)
    {
        for (uint32_t residue = 0u;
             residue < populated_residues; ++residue)
        {
            const uint64_t* coefficients =
                cached_packed.ByResidue[residue];
            const uint64_t* bucket = secondary_residue_projection.data() +
                (size_t)residue * projection_words;
            for (uint32_t word_index = 0u;
                 word_index < projection_words; ++word_index)
            {
                uint64_t word = bucket[word_index];
                while (word != 0u)
                {
                    const uint32_t bit =
                        wirehair::NonzeroLowestBitIndex64(word);
                    const uint32_t index = (word_index << 6) + bit;
                    if (index < inactive_count)
                    {
                        uint64_t* destination = projected.data() +
                            (size_t)index * packed_words;
                        for (uint32_t packed = 0u;
                             packed < packed_words; ++packed)
                        {
                            destination[packed] ^=
                                coefficients[packed] & secondary_masks[packed];
                        }
                    }
                    word &= word - 1u;
                }
            }
        }
    }
    return true;
}

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
bool ProjectMixedCompletionCoefficientsByDenseExpansion(
    uint32_t column_count,
    uint32_t inactive_count,
    uint32_t projection_words,
    const std::vector<uint32_t>& inactive_index,
    const std::vector<uint64_t>& projection,
    const MixedPackedCoefficients& cached_packed,
    std::vector<uint64_t>& projected)
{
    const uint32_t packed_words = ActiveMixedPackedCoefficientWords();
    const uint64_t projection_elements =
        (uint64_t)column_count * projection_words;
    const uint64_t projected_elements =
        (uint64_t)inactive_count * packed_words;
    if (inactive_index.size() != column_count ||
        projection_elements > projection.max_size() ||
        projection.size() != (size_t)projection_elements ||
        projected_elements > projected.max_size())
    {
        return false;
    }
    projected.assign((size_t)projected_elements, uint64_t{0});
    const uint32_t coefficient_period = ActiveMixedCoefficientPeriod();
    const uint32_t subfield_rows = ActiveMixedGF256Rows();
    const uint32_t extension_rows = ActiveMixedGF16Rows();
    const uint32_t grouped_gf256_rows = ActiveMixedGroupedGF256Rows();
    const uint32_t H = subfield_rows + extension_rows;
    if (coefficient_period < H || column_count < H ||
        grouped_gf256_rows > subfield_rows ||
        coefficient_period > kMixedCoefficientPeriod)
    {
        return false;
    }
    const uint32_t first_heavy_column = column_count - H;
    uint32_t block_index = 0u;
    uint32_t block_column = 0u;
    uint32_t residue = 0u;
    const bool independent_extension_residues =
        ActiveMixedIndependentExtensionResidues();
    const bool secondary_schedule =
        independent_extension_residues || grouped_gf256_rows != 0u;
    if (independent_extension_residues && grouped_gf256_rows != 0u) {
        return false;
    }
    const uint32_t first_grouped_gf256_row =
        subfield_rows - grouped_gf256_rows;
    uint64_t primary_masks[kMixedPackedCoefficientWords] = {};
    uint64_t secondary_masks[kMixedPackedCoefficientWords] = {};
    if (secondary_schedule)
    {
        for (uint32_t row = 0u; row < subfield_rows; ++row)
        {
            uint64_t* const masks =
                grouped_gf256_rows != 0u &&
                row >= first_grouped_gf256_row ?
                    secondary_masks : primary_masks;
            masks[row >> 2] |=
                UINT64_C(0xffff) << ((row & 3u) * 16u);
        }
        for (uint32_t er = 0u; er < extension_rows; ++er)
        {
            const uint32_t row = subfield_rows + er;
            uint64_t* const masks = independent_extension_residues ?
                secondary_masks : primary_masks;
            masks[row >> 2] |=
                UINT64_C(0xffff) << ((row & 3u) * 16u);
        }
    }
    for (uint32_t column = 0; column < column_count; ++column)
    {
        uint64_t combined_coefficients[kMixedPackedCoefficientWords] = {};
        const uint64_t* column_coefficients =
            cached_packed.ByResidue[residue];
        if (secondary_schedule)
        {
            const uint32_t secondary_residue =
                independent_extension_residues ?
                    ActiveMixedExtensionCoefficientResidue(column) :
                    ActiveMixedGroupedGF256CoefficientResidue(
                        column, first_heavy_column);
            const uint64_t* secondary_coefficients =
                cached_packed.ByResidue[
                    secondary_residue];
            for (uint32_t word = 0u; word < packed_words; ++word) {
                combined_coefficients[word] =
                    (column_coefficients[word] & primary_masks[word]) |
                    (secondary_coefficients[word] & secondary_masks[word]);
            }
            column_coefficients = combined_coefficients;
        }
        if (++block_column == coefficient_period)
        {
            block_column = 0u;
            residue = ActiveMixedResidueBlockShift(++block_index);
        }
        else if (++residue == coefficient_period) {
            residue = 0u;
        }
        const auto xor_projected = [&](uint32_t index) {
            uint64_t* destination =
                projected.data() + (size_t)index * packed_words;
            for (uint32_t word = 0; word < packed_words; ++word) {
                destination[word] ^= column_coefficients[word];
            }
        };
        const uint32_t inactive = inactive_index[column];
        if (inactive != UINT32_MAX)
        {
            if (inactive >= inactive_count) {
                return false;
            }
            xor_projected(inactive);
            continue;
        }
        const uint64_t* bits =
            projection.data() + (size_t)column * projection_words;
        for (uint32_t word_index = 0;
             word_index < projection_words; ++word_index)
        {
            uint64_t word = bits[word_index];
            while (word != 0u)
            {
                const uint32_t bit =
                    wirehair::NonzeroLowestBitIndex64(word);
                const uint32_t index = (word_index << 6) + bit;
                if (index < inactive_count) {
                    xor_projected(index);
                }
                word &= word - 1u;
            }
        }
    }
    return true;
}
#endif

#if defined(_MSC_VER)
#define WH2_MIXED_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define WH2_MIXED_NOINLINE __attribute__((noinline))
#else
#define WH2_MIXED_NOINLINE
#endif

static WH2_MIXED_NOINLINE bool AccumulateDualMixedCompletionRhs(
    uint32_t column_count,
    uint32_t coefficient_period,
    uint32_t populated_residues,
    uint32_t subfield_rows,
    uint32_t block_bytes,
    uint32_t elements,
    uint32_t extension_rows,
    const std::vector<uint32_t>& inactive_index,
    const std::vector<uint8_t>& values,
    const MixedCoefficientRows& cached_rows,
    void* const* subfield_destinations,
    std::vector<uint8_t>& rhs_low,
    std::vector<uint8_t>& rhs_high,
    std::vector<uint8_t>& source_low,
    std::vector<uint8_t>& source_high,
    PrecodeSolveStats& stats)
{
    std::vector<uint8_t> subfield_buckets(
        (size_t)coefficient_period * block_bytes, uint8_t{0});
    std::vector<uint8_t> extension_buckets(
        (size_t)coefficient_period * block_bytes, uint8_t{0});
    std::vector<BatchedBlockXorAccumulator> subfield_accumulators;
    std::vector<BatchedBlockXorAccumulator> extension_accumulators;
    subfield_accumulators.reserve(coefficient_period);
    extension_accumulators.reserve(coefficient_period);
    for (uint32_t residue = 0u; residue < coefficient_period; ++residue)
    {
        subfield_accumulators.emplace_back(
            subfield_buckets.data() + (size_t)residue * block_bytes,
            block_bytes);
        extension_accumulators.emplace_back(
            extension_buckets.data() + (size_t)residue * block_bytes,
            block_bytes);
    }
    // Walk complete coefficient blocks so the comparatively expensive hashed
    // schedule shift and integer division are paid once per P columns rather
    // than twice per column.  Both sums are below 2*P, so one subtraction is
    // the exact modulo operation used by ActiveMixed*CoefficientResidue().
    uint32_t block_start = 0u;
    uint32_t block_index = 0u;
    while (block_start < column_count)
    {
        const uint32_t subfield_shift =
            ActiveMixedResidueBlockShift(block_index);
        const uint32_t extension_shift =
            ActiveMixedExtensionResidueBlockShift(block_index);
        const uint32_t block_columns = std::min(
            coefficient_period, column_count - block_start);
        for (uint32_t offset = 0u; offset < block_columns; ++offset)
        {
            const uint32_t column = block_start + offset;
            if (inactive_index[column] != UINT32_MAX) continue;
            uint32_t subfield_residue = offset + subfield_shift;
            if (subfield_residue >= coefficient_period) {
                subfield_residue -= coefficient_period;
            }
            uint32_t extension_residue = offset + extension_shift;
            if (extension_residue >= coefficient_period) {
                extension_residue -= coefficient_period;
            }
            const uint8_t* value =
                values.data() + (size_t)column * block_bytes;
            subfield_accumulators[subfield_residue].Add(value);
            extension_accumulators[extension_residue].Add(value);
            stats.BlockXors += 2u;
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
            ++stats.MixedDualSourceColumns;
#endif
        }
        block_start += block_columns;
        ++block_index;
    }
    for (uint32_t residue = 0u; residue < coefficient_period; ++residue)
    {
        subfield_accumulators[residue].Flush();
        extension_accumulators[residue].Flush();
    }

    uint8_t subfield_scales[kMixedGF256RowsMax];
    for (uint32_t residue = 0u; residue < populated_residues; ++residue)
    {
        for (uint32_t row = 0u; row < subfield_rows; ++row) {
            subfield_scales[row] = cached_rows.Subfield[row][residue];
        }
        AddScaledBlocks(
            subfield_destinations, subfield_scales,
            subfield_rows,
            subfield_buckets.data() + (size_t)residue * block_bytes,
            block_bytes, stats);
        const uint8_t* extension_bucket =
            extension_buckets.data() + (size_t)residue * block_bytes;
        if (!DeinterleavePayloadPlanes(
                extension_bucket,
                source_low.data(), source_high.data(), block_bytes))
        {
            return false;
        }
        const uint32_t row0 = subfield_rows;
        const uint32_t row1 = row0 + 1u;
        if (!MulAddPayloadPlanes2(
                rhs_low.data() + (size_t)row0 * elements,
                rhs_high.data() + (size_t)row0 * elements,
                cached_rows.Extension[0][residue],
                rhs_low.data() + (size_t)row1 * elements,
                rhs_high.data() + (size_t)row1 * elements,
                cached_rows.Extension[1][residue],
                source_low.data(), source_high.data(), elements))
        {
            return false;
        }
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
        for (uint32_t er = 2u; er < extension_rows; ++er)
        {
            const uint32_t row = subfield_rows + er;
            if (!MulAddPayloadPlanes(
                    rhs_low.data() + (size_t)row * elements,
                    rhs_high.data() + (size_t)row * elements,
                    cached_rows.Extension[er][residue],
                    source_low.data(), source_high.data(), elements))
            {
                return false;
            }
        }
#endif
        stats.BlockMulAdds += extension_rows;
    }
    return true;
}

static WH2_MIXED_NOINLINE bool AccumulateJointMixedCompletionRhs(
    uint32_t column_count,
    uint32_t coefficient_period,
    uint32_t populated_residues,
    uint32_t subfield_rows,
    uint32_t block_bytes,
    uint32_t elements,
    uint32_t extension_rows,
    const std::vector<uint32_t>& inactive_index,
    const std::vector<uint8_t>& values,
    const MixedCoefficientRows& cached_rows,
    void* const* subfield_destinations,
    std::vector<uint8_t>& rhs_low,
    std::vector<uint8_t>& rhs_high,
    std::vector<uint8_t>& source_low,
    std::vector<uint8_t>& source_high,
    PrecodeSolveStats& stats)
{
    MixedJointResidueBuckets buckets;
    if (!AccumulateMixedJointResidueBuckets(
            column_count,
            coefficient_period,
            block_bytes,
            [&](uint32_t column) {
                return values.data() + (size_t)column * block_bytes;
            },
            [&](uint32_t column) {
                return inactive_index[column] == UINT32_MAX;
            },
            true,
            buckets,
            // The solve is the one caller whose blocks can legitimately be
            // empty; see the parameter's declaration.
            true))
    {
        return false;
    }
    stats.BlockXors += buckets.SourceXors + buckets.MarginalXors;
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    stats.MixedJointSourceXors = buckets.SourceXors;
    stats.MixedJointMarginalXors = buckets.MarginalXors;
    stats.MixedJointMarginalCopies = buckets.MarginalCopies;
    stats.MixedJointScratchBytes = buckets.ScratchBytes;
    stats.MixedJointActiveDeltas = buckets.ActiveDeltas;
#endif

    uint8_t subfield_scales[kMixedGF256RowsMax];
    for (uint32_t residue = 0u; residue < populated_residues; ++residue)
    {
        for (uint32_t row = 0u; row < subfield_rows; ++row) {
            subfield_scales[row] = cached_rows.Subfield[row][residue];
        }
        AddScaledBlocks(
            subfield_destinations, subfield_scales,
            subfield_rows,
            buckets.Subfield.get() + (size_t)residue * block_bytes,
            block_bytes, stats);
        const uint8_t* extension_bucket =
            buckets.Extension.get() + (size_t)residue * block_bytes;
        if (!DeinterleavePayloadPlanes(
                extension_bucket,
                source_low.data(), source_high.data(), block_bytes))
        {
            return false;
        }
        const uint32_t row0 = subfield_rows;
        const uint32_t row1 = row0 + 1u;
        if (!MulAddPayloadPlanes2(
                rhs_low.data() + (size_t)row0 * elements,
                rhs_high.data() + (size_t)row0 * elements,
                cached_rows.Extension[0][residue],
                rhs_low.data() + (size_t)row1 * elements,
                rhs_high.data() + (size_t)row1 * elements,
                cached_rows.Extension[1][residue],
                source_low.data(), source_high.data(), elements))
        {
            return false;
        }
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
        for (uint32_t er = 2u; er < extension_rows; ++er)
        {
            const uint32_t row = subfield_rows + er;
            if (!MulAddPayloadPlanes(
                    rhs_low.data() + (size_t)row * elements,
                    rhs_high.data() + (size_t)row * elements,
                    cached_rows.Extension[er][residue],
                    source_low.data(), source_high.data(), elements))
            {
                return false;
            }
        }
#endif
        stats.BlockMulAdds += extension_rows;
    }
    return true;
}

static WH2_MIXED_NOINLINE bool ApplyGroupedMixedCompletionBuckets(
    uint32_t coefficient_period,
    uint32_t populated_residues,
    uint32_t subfield_rows,
    uint32_t first_grouped_gf256_row,
    uint32_t block_bytes,
    uint32_t elements,
    uint32_t extension_rows,
    const MixedCoefficientRows& cached_rows,
    void* const* subfield_destinations,
    std::vector<uint8_t>& rhs_low,
    std::vector<uint8_t>& rhs_high,
    std::vector<uint8_t>& source_low,
    std::vector<uint8_t>& source_high,
    const uint8_t* primary_buckets,
    const uint8_t* grouped_buckets,
    PrecodeSolveStats& stats)
{
    if (first_grouped_gf256_row >= subfield_rows ||
        extension_rows < 2u ||
        populated_residues > coefficient_period)
    {
        return false;
    }
    const uint32_t grouped_gf256_rows =
        subfield_rows - first_grouped_gf256_row;
    uint8_t scales[kMixedGF256RowsMax];
    for (uint32_t residue = 0u; residue < populated_residues; ++residue)
    {
        const uint8_t* const primary_bucket =
            primary_buckets + (size_t)residue * block_bytes;
        const uint8_t* const grouped_bucket =
            grouped_buckets + (size_t)residue * block_bytes;
        for (uint32_t row = 0u;
             row < first_grouped_gf256_row; ++row)
        {
            scales[row] = cached_rows.Subfield[row][residue];
        }
        AddScaledBlocks(
            subfield_destinations, scales,
            first_grouped_gf256_row, primary_bucket,
            block_bytes, stats);
        for (uint32_t row = 0u; row < grouped_gf256_rows; ++row) {
            scales[row] = cached_rows.Subfield[
                first_grouped_gf256_row + row][residue];
        }
        AddScaledBlocks(
            subfield_destinations + first_grouped_gf256_row,
            scales, grouped_gf256_rows, grouped_bucket,
            block_bytes, stats);

        if (!DeinterleavePayloadPlanes(
                primary_bucket,
                source_low.data(), source_high.data(), block_bytes))
        {
            return false;
        }
        const uint32_t row0 = subfield_rows;
        const uint32_t row1 = row0 + 1u;
        if (!MulAddPayloadPlanes2(
                rhs_low.data() + (size_t)row0 * elements,
                rhs_high.data() + (size_t)row0 * elements,
                cached_rows.Extension[0][residue],
                rhs_low.data() + (size_t)row1 * elements,
                rhs_high.data() + (size_t)row1 * elements,
                cached_rows.Extension[1][residue],
                source_low.data(), source_high.data(), elements))
        {
            return false;
        }
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
        for (uint32_t er = 2u; er < extension_rows; ++er)
        {
            const uint32_t row = subfield_rows + er;
            if (!MulAddPayloadPlanes(
                    rhs_low.data() + (size_t)row * elements,
                    rhs_high.data() + (size_t)row * elements,
                    cached_rows.Extension[er][residue],
                    source_low.data(), source_high.data(), elements))
            {
                return false;
            }
        }
#endif
        stats.BlockMulAdds += extension_rows;
    }
    return true;
}

static WH2_MIXED_NOINLINE bool AccumulateDualGroupedMixedCompletionRhs(
    uint32_t column_count,
    uint32_t first_heavy_column,
    uint32_t coefficient_period,
    uint32_t populated_residues,
    uint32_t subfield_rows,
    uint32_t first_grouped_gf256_row,
    uint32_t block_bytes,
    uint32_t elements,
    uint32_t extension_rows,
    const uint8_t* grouped_block_shifts,
    uint32_t grouped_block_shift_count,
    const std::vector<uint32_t>& inactive_index,
    const std::vector<uint8_t>& values,
    const MixedCoefficientRows& cached_rows,
    void* const* subfield_destinations,
    std::vector<uint8_t>& rhs_low,
    std::vector<uint8_t>& rhs_high,
    std::vector<uint8_t>& source_low,
    std::vector<uint8_t>& source_high,
    PrecodeSolveStats& stats)
{
    if (coefficient_period == 0u || first_heavy_column > column_count) {
        return false;
    }
    const uint32_t expected_grouped_blocks =
        first_heavy_column / coefficient_period +
        (first_heavy_column % coefficient_period != 0u ? 1u : 0u);
    if (grouped_block_shift_count != expected_grouped_blocks ||
        (expected_grouped_blocks != 0u && !grouped_block_shifts))
    {
        return false;
    }
    std::vector<uint8_t> primary_buckets(
        (size_t)coefficient_period * block_bytes, uint8_t{0});
    std::vector<uint8_t> grouped_buckets(
        (size_t)coefficient_period * block_bytes, uint8_t{0});
    std::vector<BatchedBlockXorAccumulator> primary_accumulators;
    std::vector<BatchedBlockXorAccumulator> grouped_accumulators;
    primary_accumulators.reserve(coefficient_period);
    grouped_accumulators.reserve(coefficient_period);
    for (uint32_t residue = 0u; residue < coefficient_period; ++residue)
    {
        primary_accumulators.emplace_back(
            primary_buckets.data() + (size_t)residue * block_bytes,
            block_bytes);
        grouped_accumulators.emplace_back(
            grouped_buckets.data() + (size_t)residue * block_bytes,
            block_bytes);
    }
    uint32_t block_start = 0u;
    uint32_t block_index = 0u;
    while (block_start < column_count)
    {
        const uint32_t primary_shift =
            ActiveMixedResidueBlockShift(block_index);
        uint32_t grouped_shift = primary_shift;
        if (block_start < first_heavy_column)
        {
            CAT_DEBUG_ASSERT(block_index < grouped_block_shift_count);
            grouped_shift = grouped_block_shifts[block_index];
        }
        const uint32_t block_columns = std::min(
            coefficient_period, column_count - block_start);
        for (uint32_t offset = 0u; offset < block_columns; ++offset)
        {
            const uint32_t column = block_start + offset;
            if (inactive_index[column] != UINT32_MAX) continue;
            uint32_t primary_residue = offset + primary_shift;
            if (primary_residue >= coefficient_period) {
                primary_residue -= coefficient_period;
            }
            uint32_t grouped_residue = primary_residue;
            if (column < first_heavy_column)
            {
                grouped_residue = offset + grouped_shift;
                if (grouped_residue >= coefficient_period) {
                    grouped_residue -= coefficient_period;
                }
            }
            const uint8_t* value =
                values.data() + (size_t)column * block_bytes;
            primary_accumulators[primary_residue].Add(value);
            grouped_accumulators[grouped_residue].Add(value);
            stats.BlockXors += 2u;
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
            ++stats.MixedDualSourceColumns;
#endif
        }
        block_start += block_columns;
        ++block_index;
    }
    for (uint32_t residue = 0u; residue < coefficient_period; ++residue)
    {
        primary_accumulators[residue].Flush();
        grouped_accumulators[residue].Flush();
    }
    return ApplyGroupedMixedCompletionBuckets(
        coefficient_period, populated_residues, subfield_rows,
        first_grouped_gf256_row, block_bytes, elements, extension_rows,
        cached_rows, subfield_destinations,
        rhs_low, rhs_high, source_low, source_high,
        primary_buckets.data(), grouped_buckets.data(), stats);
}

static WH2_MIXED_NOINLINE bool AccumulateJointGroupedMixedCompletionRhs(
    uint32_t column_count,
    uint32_t first_heavy_column,
    uint32_t coefficient_period,
    uint32_t populated_residues,
    uint32_t subfield_rows,
    uint32_t first_grouped_gf256_row,
    uint32_t block_bytes,
    uint32_t elements,
    uint32_t extension_rows,
    const uint8_t* grouped_block_shifts,
    uint32_t grouped_block_shift_count,
    const std::vector<uint32_t>& inactive_index,
    const std::vector<uint8_t>& values,
    const MixedCoefficientRows& cached_rows,
    void* const* subfield_destinations,
    std::vector<uint8_t>& rhs_low,
    std::vector<uint8_t>& rhs_high,
    std::vector<uint8_t>& source_low,
    std::vector<uint8_t>& source_high,
    PrecodeSolveStats& stats)
{
    if (coefficient_period == 0u || first_heavy_column > column_count) {
        return false;
    }
    const uint32_t expected_grouped_blocks =
        first_heavy_column / coefficient_period +
        (first_heavy_column % coefficient_period != 0u ? 1u : 0u);
    if (grouped_block_shift_count != expected_grouped_blocks ||
        (expected_grouped_blocks != 0u && !grouped_block_shifts))
    {
        return false;
    }
    MixedJointResidueBuckets buckets;
    if (!AccumulateMixedJointResidueBucketsWithShifts(
            first_heavy_column,
            coefficient_period,
            block_bytes,
            [](uint32_t block) {
                return ActiveMixedResidueBlockShift(block);
            },
            [&](uint32_t block) {
                CAT_DEBUG_ASSERT(block < grouped_block_shift_count);
                return grouped_block_shifts[block];
            },
            [&](uint32_t column) {
                return values.data() + (size_t)column * block_bytes;
            },
            [&](uint32_t column) {
                return inactive_index[column] == UINT32_MAX;
            },
            true,
            buckets,
            // The solve is the one caller whose blocks can legitimately be
            // empty; see the parameter's declaration.
            true))
    {
        return false;
    }

    uint64_t tail_xors = 0u;
    for (uint32_t column = first_heavy_column;
         column < column_count; ++column)
    {
        if (inactive_index[column] != UINT32_MAX) continue;
        const uint32_t residue = ActiveMixedCoefficientResidue(column);
        const uint8_t* const source =
            values.data() + (size_t)column * block_bytes;
        gf256_add_mem(
            buckets.Subfield.get() + (size_t)residue * block_bytes,
            source, (int)block_bytes);
        gf256_add_mem(
            buckets.Extension.get() + (size_t)residue * block_bytes,
            source, (int)block_bytes);
        tail_xors += 2u;
    }
    stats.BlockXors +=
        buckets.SourceXors + buckets.MarginalXors + tail_xors;
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    stats.MixedJointSourceXors = buckets.SourceXors + tail_xors;
    stats.MixedJointMarginalXors = buckets.MarginalXors;
    stats.MixedJointMarginalCopies = buckets.MarginalCopies;
    stats.MixedJointScratchBytes = buckets.ScratchBytes;
    stats.MixedJointActiveDeltas = buckets.ActiveDeltas;
#endif
    return ApplyGroupedMixedCompletionBuckets(
        coefficient_period, populated_residues, subfield_rows,
        first_grouped_gf256_row, block_bytes, elements, extension_rows,
        cached_rows, subfield_destinations,
        rhs_low, rhs_high, source_low, source_high,
        buckets.Subfield.get(), buckets.Extension.get(), stats);
}

#undef WH2_MIXED_NOINLINE

static const uint32_t kMixedCompletionRowsMax =
    kMixedGF256RowsMax + kMixedGF16RowsMax;
static_assert(
    kMixedCompletionRowsMax < UINT8_MAX,
    "mixed quotient row indices require an out-of-band byte sentinel");

struct MixedQuotientRowPlan
{
    std::array<uint16_t, kMixedCompletionRowsMax> ForwardScales = {};
    std::array<uint16_t, kMixedCompletionRowsMax> BackScales = {};
    uint16_t NormalizeScale = 1u;
    uint8_t PivotColumn = UINT8_MAX;
};

struct MixedQuotientFactorization
{
    std::array<MixedQuotientRowPlan, kMixedCompletionRowsMax> Rows = {};
    std::array<uint8_t, kMixedCompletionRowsMax> PivotRows = {};
    std::array<
        uint16_t,
        kMixedCompletionRowsMax * kMixedCompletionRowsMax> Transform = {};
    uint32_t Rank = 0u;
};

// Factor the tiny mixed-field quotient before touching any payload bytes.
// The row plans replay the exact historical Gauss-Jordan order on RHS blocks,
// without allocating or copying separate pivot payloads.  On the uncommon
// deficient path BuildMixedQuotientTransform replays the same plans on an
// identity matrix: each dependent row is then an exact syndrome equation,
// allowing NeedMore to remain distinguishable from an inconsistent RHS
// without constructing all H completion payload rows.
bool FactorMixedCompletionQuotient(
    std::vector<uint16_t>& quotient_rows,
    uint32_t row_count,
    uint32_t column_count,
    MixedQuotientFactorization& factor)
{
    if (row_count > kMixedCompletionRowsMax ||
        column_count > row_count ||
        quotient_rows.size() != (size_t)row_count * column_count)
    {
        return false;
    }

    factor = MixedQuotientFactorization{};
    factor.PivotRows.fill(UINT8_MAX);

    for (uint32_t row = 0u; row < row_count; ++row)
    {
        MixedQuotientRowPlan& plan = factor.Rows[row];
        uint16_t* coeff = column_count == 0u ? nullptr :
            quotient_rows.data() + (size_t)row * column_count;
        for (uint32_t column = 0u; column < column_count; ++column)
        {
            const uint16_t scale = coeff[column];
            const uint8_t pivot_row = factor.PivotRows[column];
            if (scale == 0u || pivot_row == UINT8_MAX) continue;
            plan.ForwardScales[column] = scale;
            const uint16_t* pivot = quotient_rows.data() +
                (size_t)pivot_row * column_count;
            for (uint32_t k = column; k < column_count; ++k) {
                coeff[k] ^= GF16MultiplyInitialized(scale, pivot[k]);
            }
        }

        uint32_t pivot_column = column_count;
        for (uint32_t column = 0u; column < column_count; ++column) {
            if (coeff[column] != 0u) {
                pivot_column = column;
                break;
            }
        }
        if (pivot_column == column_count) {
            continue;
        }

        const uint16_t pivot_value = coeff[pivot_column];
        if (pivot_value != 1u)
        {
            const uint16_t inverse = GF16InverseInitialized(pivot_value);
            if (inverse == 0u) return false;
            plan.NormalizeScale = inverse;
            for (uint32_t k = pivot_column; k < column_count; ++k) {
                coeff[k] = GF16MultiplyInitialized(coeff[k], inverse);
            }
        }

        for (uint32_t existing = 0u;
             existing < column_count; ++existing)
        {
            const uint8_t existing_row = factor.PivotRows[existing];
            if (existing_row == UINT8_MAX) continue;
            uint16_t* existing_coeff = quotient_rows.data() +
                (size_t)existing_row * column_count;
            const uint16_t scale = existing_coeff[pivot_column];
            if (scale == 0u) continue;
            plan.BackScales[existing] = scale;
            for (uint32_t k = pivot_column; k < column_count; ++k) {
                existing_coeff[k] ^=
                    GF16MultiplyInitialized(scale, coeff[k]);
            }
        }
        plan.PivotColumn = (uint8_t)pivot_column;
        factor.PivotRows[pivot_column] = (uint8_t)row;
        ++factor.Rank;
    }
    return true;
}

bool BuildMixedQuotientTransform(
    uint32_t row_count,
    uint32_t column_count,
    MixedQuotientFactorization& factor)
{
    if (row_count > kMixedCompletionRowsMax ||
        column_count > row_count)
    {
        return false;
    }
    factor.Transform.fill(0u);
    for (uint32_t row = 0u; row < row_count; ++row) {
        factor.Transform[(size_t)row * kMixedCompletionRowsMax + row] = 1u;
    }
    const auto muladd = [&factor, row_count](
        uint32_t destination,
        uint32_t source,
        uint16_t scale)
    {
        uint16_t* dst = factor.Transform.data() +
            (size_t)destination * kMixedCompletionRowsMax;
        const uint16_t* src = factor.Transform.data() +
            (size_t)source * kMixedCompletionRowsMax;
        for (uint32_t i = 0u; i < row_count; ++i) {
            dst[i] ^= GF16MultiplyInitialized(scale, src[i]);
        }
    };
    const auto scale_row = [&factor, row_count](
        uint32_t row,
        uint16_t scale)
    {
        uint16_t* values = factor.Transform.data() +
            (size_t)row * kMixedCompletionRowsMax;
        for (uint32_t i = 0u; i < row_count; ++i) {
            values[i] = GF16MultiplyInitialized(values[i], scale);
        }
    };

    for (uint32_t row = 0u; row < row_count; ++row)
    {
        const MixedQuotientRowPlan& plan = factor.Rows[row];
        for (uint32_t column = 0u; column < column_count; ++column)
        {
            const uint16_t value = plan.ForwardScales[column];
            if (value == 0u) continue;
            const uint8_t pivot_row = factor.PivotRows[column];
            if (pivot_row == UINT8_MAX || pivot_row == row) return false;
            muladd(row, pivot_row, value);
        }
        if (plan.PivotColumn == UINT8_MAX) continue;
        if (plan.PivotColumn >= column_count ||
            factor.PivotRows[plan.PivotColumn] != row)
        {
            return false;
        }
        if (plan.NormalizeScale != 1u) {
            scale_row(row, plan.NormalizeScale);
        }
        for (uint32_t existing = 0u;
             existing < column_count; ++existing)
        {
            const uint16_t value = plan.BackScales[existing];
            if (value == 0u) continue;
            const uint8_t existing_row = factor.PivotRows[existing];
            if (existing_row == UINT8_MAX || existing_row == row) {
                return false;
            }
            muladd(existing_row, row, value);
        }
    }
    return true;
}

// ---------------------------------------------------------------------------
// Tiny-K mixed completion fast path.
//
// At tiny system sizes the general mixed quotient pipeline (packed
// coefficient projection by residue buckets, bulk-kernel residue-bucket RHS
// accumulation, GF(2^16) factorization plans and planar Gauss-Jordan replay)
// is dominated by per-call and setup overhead: the actual algebra is a dense
// solve over a handful of columns.  The helpers below evaluate the SAME
// equations with scalar GF(2^16) arithmetic:
//
//   * per-column completion coefficients follow the exact dense-expansion
//     rule (the oracle the residue-bucket projection is verified against),
//   * the heavy RHS is the same XOR/muladd sum accumulated column by column,
//   * the quotient is solved by plain Gauss-Jordan elimination.
//
// Every algebraic outcome is identical to the general path: full-rank
// systems have a unique solution in exact field arithmetic, rank is
// pivot-order independent, and the dependent completion rows' reduced RHS
// blocks are a basis of the same left-nullspace syndromes the factorization
// transform checks, so Success/NeedMore/Error classification and every
// output byte match.  The fast path never writes to `values` before success
// is fully decided, preserving the transactional failure contract.
// ---------------------------------------------------------------------------

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
thread_local int TinyMixedFastPathTestModeValue = 0;
#endif

// Scalar GF(2^16) row helpers.  Elements are little-endian 16-bit values in
// the interleaved block layout.  Multiplication by scale (a | b << 8)
// decomposes over the GF(256) subfield exactly as GF16MulAddPlanar():
//   low'  ^= a * low ^ (lambda * b) * high
//   high' ^= b * low ^ (a ^ b) * high
struct TinyGF16RowTables
{
    const uint8_t* A;
    const uint8_t* B;
    const uint8_t* LB;
    const uint8_t* AB;
};

inline TinyGF16RowTables TinyGF16Rows(uint16_t scale)
{
    const uint8_t a = (uint8_t)scale;
    const uint8_t b = (uint8_t)(scale >> 8);
    const uint8_t lb = gf256_mul(kGF16Lambda, b);
    const uint8_t ab = (uint8_t)(a ^ b);
    TinyGF16RowTables tables;
    tables.A = &GF256Ctx.GF256_MUL_TABLE[(size_t)a << 8];
    tables.B = &GF256Ctx.GF256_MUL_TABLE[(size_t)b << 8];
    tables.LB = &GF256Ctx.GF256_MUL_TABLE[(size_t)lb << 8];
    tables.AB = &GF256Ctx.GF256_MUL_TABLE[(size_t)ab << 8];
    return tables;
}

inline void TinyGF16MulAddRow(
    uint8_t* GF256_RESTRICT destination,
    const uint8_t* GF256_RESTRICT source,
    uint16_t scale,
    uint32_t block_bytes)
{
    if (scale == 0u) {
        return;
    }
    if (scale == 1u) {
        gf256_add_mem(destination, source, (int)block_bytes);
        return;
    }
    const TinyGF16RowTables tables = TinyGF16Rows(scale);
    for (uint32_t offset = 0; offset < block_bytes; offset += 2u)
    {
        const uint8_t low = source[offset];
        const uint8_t high = source[offset + 1u];
        destination[offset] ^= (uint8_t)(tables.A[low] ^ tables.LB[high]);
        destination[offset + 1u] ^=
            (uint8_t)(tables.B[low] ^ tables.AB[high]);
    }
}

inline void TinyGF16ScaleRow(
    uint8_t* block,
    uint16_t scale,
    uint32_t block_bytes)
{
    if (scale == 1u) {
        return;
    }
    if (scale == 0u) {
        std::memset(block, 0, block_bytes);
        return;
    }
    const TinyGF16RowTables tables = TinyGF16Rows(scale);
    for (uint32_t offset = 0; offset < block_bytes; offset += 2u)
    {
        const uint8_t low = block[offset];
        const uint8_t high = block[offset + 1u];
        block[offset] = (uint8_t)(tables.A[low] ^ tables.LB[high]);
        block[offset + 1u] = (uint8_t)(tables.B[low] ^ tables.AB[high]);
    }
}

// Coefficient sub-rows are uint16 arrays.  On little-endian hosts their
// in-memory layout matches the little-endian block element format, so the
// fused byte kernels apply directly; other hosts use the equivalent scalar
// loop.  Both compute the identical field products.
#if defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__) && \
    (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)
#define WH2_TINY_LITTLE_ENDIAN 1
#elif defined(_MSC_VER) && (defined(_M_IX86) || defined(_M_X64) || \
    defined(_M_ARM) || defined(_M_ARM64))
#define WH2_TINY_LITTLE_ENDIAN 1
#else
#define WH2_TINY_LITTLE_ENDIAN 0
#endif

inline void TinyGF16MulAddCoefficients(
    uint16_t* destination,
    const uint16_t* source,
    uint16_t scale,
    uint32_t count);

inline void TinyGF16ScaleCoefficients(
    uint16_t* coefficients,
    uint16_t scale,
    uint32_t count);

// Exact GF(2^16) inverse through the GF(256) subfield norm.  With
// x^2 = x + lambda the conjugate of z = a + b x is (a ^ b) + b x, so
// z * conj(z) = a^2 ^ a b ^ lambda b^2 lies in GF(256) and
// z^-1 = conj(z) * norm^-1.  This equals GF16InverseInitialized() for every
// input (the field inverse is unique); the A/B unit test proves it
// exhaustively.
inline uint16_t TinyGF16Inverse(uint16_t value)
{
    const uint8_t a = (uint8_t)value;
    const uint8_t b = (uint8_t)(value >> 8);
    const uint8_t norm = (uint8_t)(
        gf256_mul(a, (uint8_t)(a ^ b)) ^
        gf256_mul(kGF16Lambda, gf256_mul(b, b)));
    if (norm == 0u) {
        return 0u; // norm(z) == 0 only for z == 0
    }
    const uint8_t norm_inverse = gf256_inv(norm);
    return (uint16_t)gf256_mul((uint8_t)(a ^ b), norm_inverse) |
        (uint16_t)((uint16_t)gf256_mul(b, norm_inverse) << 8);
}

inline void TinyGF16MulAddCoefficients(
    uint16_t* destination,
    const uint16_t* source,
    uint16_t scale,
    uint32_t count)
{
    if (count == 0u || scale == 0u) {
        return;
    }
#if WH2_TINY_LITTLE_ENDIAN
    TinyGF16MulAddRow(
        reinterpret_cast<uint8_t*>(destination),
        reinterpret_cast<const uint8_t*>(source),
        scale, count * 2u);
#else
    for (uint32_t i = 0u; i < count; ++i) {
        destination[i] = (uint16_t)(destination[i] ^
            GF16MultiplyInitialized(scale, source[i]));
    }
#endif
}

inline void TinyGF16ScaleCoefficients(
    uint16_t* coefficients,
    uint16_t scale,
    uint32_t count)
{
    if (count == 0u) {
        return;
    }
#if WH2_TINY_LITTLE_ENDIAN
    TinyGF16ScaleRow(
        reinterpret_cast<uint8_t*>(coefficients), scale, count * 2u);
#else
    for (uint32_t i = 0u; i < count; ++i) {
        coefficients[i] = GF16MultiplyInitialized(scale, coefficients[i]);
    }
#endif
}

bool UseTinyMixedCompletionFastPath(
    uint32_t block_count,
    uint32_t block_bytes)
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    if (TinyMixedFastPathTestModeValue < 0) {
        return false;
    }
    if (TinyMixedFastPathTestModeValue > 0) {
        return true;
    }
    // An explicitly requested residue-bucket mode is a demand for that
    // projection machinery (experiments assert its use); the automatic
    // tiny-path gate must not silently bypass it.
    if (ActiveMixedResidueBucketModeForTesting() !=
        MixedResidueBucketMode::Automatic)
    {
        return false;
    }

#endif
    // Zero width is a MEASUREMENT mode, not a tiny payload, and must run the
    // production algorithm.
    //
    // A zero-byte solve exists so op counters can be harvested without moving
    // payload -- the counts are the point.  Routing it here defeated exactly
    // that: the tiny path is a different algorithm with a different operation
    // mix, so the counters described work the codec never performs.  Measured
    // at K=128 on one fixed loss stream, varying only the width:
    //
    //     width 0,2,4,8       xors 1696   muladds   80   <- this path
    //     width 64 and above  xors 1645   muladds 1493   <- production path
    //
    // Mul-adds were understated 18.7x, and a per-site tally showed a disjoint
    // set of call sites rather than a shortfall on shared ones.  Because the
    // search's objective is a linear cost model over these counters, it ranked
    // candidates by that phantom mix: Spearman against real 64 KiB timing was
    // 0.527 and predicted cost sat at ~1.6M ns where the stopwatch read ~3.1M.
    //
    // Genuinely tiny payloads (1..8 bytes) still take the fast path, which is
    // what it is for.  Zero is not a payload size.
    //
    // This is now the same DISPATCH rule the rest of the file uses: zero asks
    // the gate at kCountingDispatchBlockBytes, which a static_assert keeps
    // strictly above TinyMixedFastPathMaxBlockBytes, so the answer is "no
    // fast path" for exactly the reason it should be -- the width it
    // dispatches as is far too wide -- and 1..8 still answer "yes".
    const uint32_t dispatch_block_bytes = DispatchBlockBytes(block_bytes);
    return block_count <= TinyMixedFastPathMaxSourceBlocks() &&
        dispatch_block_bytes <= TinyMixedFastPathMaxBlockBytes() &&
        (uint64_t)block_count * dispatch_block_bytes <=
            TinyMixedFastPathMaxProductBytes();
}

/**
    Direct dense solve of the mixed completion quotient.

    Returns false only for structurally unsupported shapes (the caller then
    falls through to the general machinery having built nothing).  When it
    returns true, `result_out` and all outputs are exactly what the general
    path would have produced for the same inputs.
*/
bool TrySolveMixedCompletionQuotientTiny(
    uint32_t column_count,
    uint32_t inactive_count,
    uint32_t quotient_columns,
    uint32_t subfield_rows,
    uint32_t extension_rows,
    uint32_t grouped_gf256_rows,
    bool independent_extension_residues,
    uint32_t projection_words,
    uint32_t block_bytes,
    const std::vector<uint32_t>& inactive_index,
    const std::vector<uint32_t>& inactive_columns,
    const std::vector<uint64_t>& projection,
    const std::vector<uint64_t>& binary_pivot_coeff,
    const std::vector<uint8_t>& binary_pivot_rhs,
    const std::vector<uint8_t>& binary_have_pivot,
    uint32_t binary_rank,
    const std::vector<uint32_t>& free_columns,
    std::vector<uint8_t>& values,
    PrecodeSolveStats& stats,
    WirehairResult& result_out)
{
    const uint32_t H = subfield_rows + extension_rows;
    const uint32_t packed_words = ActiveMixedPackedCoefficientWords();
    const uint32_t coefficient_period = ActiveMixedCoefficientPeriod();
    if (H == 0u || H > kMixedCompletionRowsMax ||
        quotient_columns > H ||
        packed_words == 0u ||
        packed_words > kMixedPackedCoefficientWords ||
        packed_words * 4u < H ||
        coefficient_period == 0u ||
        coefficient_period > kMixedCoefficientPeriod ||
        column_count < H ||
        grouped_gf256_rows > subfield_rows ||
        (grouped_gf256_rows != 0u && independent_extension_residues) ||
        free_columns.size() != quotient_columns)
    {
        return false;
    }
    const MixedPackedCoefficients* cached_packed =
        GetMixedPackedCoefficients();
    if (!cached_packed) {
        return false;
    }

    const uint32_t first_heavy_column = column_count - H;
    const bool secondary_schedule =
        independent_extension_residues || grouped_gf256_rows != 0u;
    const uint32_t first_grouped_gf256_row =
        subfield_rows - grouped_gf256_rows;
    uint64_t primary_masks[kMixedPackedCoefficientWords] = {};
    uint64_t secondary_masks[kMixedPackedCoefficientWords] = {};
    if (secondary_schedule)
    {
        for (uint32_t row = 0u; row < subfield_rows; ++row)
        {
            uint64_t* const masks =
                grouped_gf256_rows != 0u &&
                row >= first_grouped_gf256_row ?
                    secondary_masks : primary_masks;
            masks[row >> 2] |=
                UINT64_C(0xffff) << ((row & 3u) * 16u);
        }
        for (uint32_t er = 0u; er < extension_rows; ++er)
        {
            const uint32_t row = subfield_rows + er;
            uint64_t* const masks = independent_extension_residues ?
                secondary_masks : primary_masks;
            masks[row >> 2] |=
                UINT64_C(0xffff) << ((row & 3u) * 16u);
        }
    }

    std::vector<uint64_t> projected(
        (size_t)inactive_count * packed_words, uint64_t{0});
    std::vector<uint8_t> rhs = MakeBlockStorage((size_t)H * block_bytes);
    // The general path counts the same H-block RHS allocation zero-fill,
    // so counting it here keeps the receipt independent of which path the
    // block width selected.
    stats.BlockZeroFills += H;

    // One dense pass over every solver column: completion coefficients
    // project onto the inactivated variables while each peeled column's
    // affine constant accumulates directly into the heavy RHS blocks.
    {
        uint32_t block_index = 0u;
        uint32_t block_column = 0u;
        uint32_t residue = 0u;
        for (uint32_t column = 0; column < column_count; ++column)
        {
            uint64_t combined[kMixedPackedCoefficientWords] = {};
            const uint64_t* column_coefficients =
                cached_packed->ByResidue[residue];
            if (secondary_schedule)
            {
                const uint32_t secondary_residue =
                    independent_extension_residues ?
                        ActiveMixedExtensionCoefficientResidue(column) :
                        ActiveMixedGroupedGF256CoefficientResidue(
                            column, first_heavy_column);
                if (secondary_residue >= coefficient_period) {
                    return false;
                }
                const uint64_t* secondary_coefficients =
                    cached_packed->ByResidue[secondary_residue];
                for (uint32_t word = 0u; word < packed_words; ++word) {
                    combined[word] =
                        (column_coefficients[word] &
                            primary_masks[word]) |
                        (secondary_coefficients[word] &
                            secondary_masks[word]);
                }
                column_coefficients = combined;
            }
            if (++block_column == coefficient_period)
            {
                block_column = 0u;
                residue = ActiveMixedResidueBlockShift(++block_index);
                if (residue >= coefficient_period) {
                    return false;
                }
            }
            else if (++residue == coefficient_period) {
                residue = 0u;
            }
            const uint32_t inactive = inactive_index[column];
            if (inactive != UINT32_MAX)
            {
                if (inactive >= inactive_count) {
                    return false;
                }
                uint64_t* destination =
                    projected.data() + (size_t)inactive * packed_words;
                for (uint32_t word = 0u; word < packed_words; ++word) {
                    destination[word] ^= column_coefficients[word];
                }
                continue;
            }
            const uint64_t* bits =
                projection.data() + (size_t)column * projection_words;
            for (uint32_t word_index = 0;
                 word_index < projection_words; ++word_index)
            {
                uint64_t word = bits[word_index];
                while (word != 0u)
                {
                    const uint32_t bit =
                        wirehair::NonzeroLowestBitIndex64(word);
                    const uint32_t index = (word_index << 6) + bit;
                    if (index < inactive_count)
                    {
                        uint64_t* destination = projected.data() +
                            (size_t)index * packed_words;
                        for (uint32_t w = 0u; w < packed_words; ++w) {
                            destination[w] ^= column_coefficients[w];
                        }
                    }
                    word &= word - 1u;
                }
            }
            const uint8_t* value =
                values.data() + (size_t)column * block_bytes;
            if (RowIsZero(value, block_bytes)) {
                continue;
            }
            for (uint32_t row = 0u; row < H; ++row)
            {
                const uint16_t scale = (uint16_t)(
                    column_coefficients[row >> 2] >>
                        ((row & 3u) * 16u));
                if (scale == 0u) {
                    continue;
                }
                TinyGF16MulAddRow(
                    rhs.data() + (size_t)row * block_bytes,
                    value, scale, block_bytes);
                if (scale == 1u) {
                    ++stats.BlockXors;
                }
                else {
                    ++stats.BlockMulAdds;
                }
            }
        }
    }

    // Reduce every binary pivot RHS through its projected completion
    // scales, exactly as the general path does after its bucket pass.
    for (uint32_t pivot = 0u; pivot < inactive_count; ++pivot)
    {
        if (!binary_have_pivot[pivot]) {
            continue;
        }
        const uint64_t* packed_scales =
            projected.data() + (size_t)pivot * packed_words;
        bool have_scale = false;
        for (uint32_t row = 0u; row < H && !have_scale; ++row) {
            have_scale = (uint16_t)(
                packed_scales[row >> 2] >> ((row & 3u) * 16u)) != 0u;
        }
        if (!have_scale) {
            continue;
        }
        const uint8_t* source =
            binary_pivot_rhs.data() + (size_t)pivot * block_bytes;
        if (RowIsZero(source, block_bytes)) {
            continue;
        }
        for (uint32_t row = 0u; row < H; ++row)
        {
            const uint16_t scale = (uint16_t)(
                packed_scales[row >> 2] >> ((row & 3u) * 16u));
            if (scale == 0u) {
                continue;
            }
            TinyGF16MulAddRow(
                rhs.data() + (size_t)row * block_bytes,
                source, scale, block_bytes);
            if (scale == 1u) {
                ++stats.BlockXors;
            }
            else {
                ++stats.BlockMulAdds;
            }
        }
    }

    // Quotient coefficients over the free columns via the same binary
    // pivot substitution the general path performs.
    uint16_t quotient[kMixedCompletionRowsMax][kMixedCompletionRowsMax];
    for (uint32_t i = 0u; i < quotient_columns; ++i)
    {
        const uint32_t free_column = free_columns[i];
        if (free_column >= inactive_count) {
            return false;
        }
        uint64_t packed[kMixedPackedCoefficientWords] = {};
        const uint64_t* free_coefficients =
            projected.data() + (size_t)free_column * packed_words;
        for (uint32_t word = 0u; word < packed_words; ++word) {
            packed[word] = free_coefficients[word];
        }
        for (uint32_t pivot = 0u; pivot < inactive_count; ++pivot)
        {
            const uint64_t* relation =
                binary_pivot_coeff.data() +
                    (size_t)pivot * projection_words;
            if (!binary_have_pivot[pivot] ||
                (relation[free_column >> 6] &
                    (UINT64_C(1) << (free_column & 63u))) == 0u)
            {
                continue;
            }
            const uint64_t* pivot_coefficients =
                projected.data() + (size_t)pivot * packed_words;
            for (uint32_t word = 0u; word < packed_words; ++word) {
                packed[word] ^= pivot_coefficients[word];
            }
        }
        for (uint32_t row = 0u; row < H; ++row) {
            quotient[row][i] = (uint16_t)(
                packed[row >> 2] >> ((row & 3u) * 16u));
        }
    }

    // Dense scalar Gauss-Jordan over GF(2^16) carrying the RHS blocks.
    // Any pivot choice reaches the same rank and, at full rank, the same
    // unique solution; skipped columns only occur on the NeedMore path.
    uint8_t row_used[kMixedCompletionRowsMax] = {};
    uint8_t pivot_row_of[kMixedCompletionRowsMax];
    std::memset(pivot_row_of, 0xff, sizeof(pivot_row_of));
    uint32_t rank = 0u;
    for (uint32_t column = 0u; column < quotient_columns; ++column)
    {
        uint32_t pivot_row = H;
        for (uint32_t row = 0u; row < H; ++row)
        {
            if (!row_used[row] && quotient[row][column] != 0u) {
                pivot_row = row;
                break;
            }
        }
        if (pivot_row == H) {
            continue;
        }
        row_used[pivot_row] = 1u;
        pivot_row_of[column] = (uint8_t)pivot_row;
        ++rank;
        // Earlier pivot columns are already zero in every row, so only the
        // tail columns need updating; the fused coefficient kernels avoid
        // one out-of-line scalar multiply per element.
        const uint32_t tail_columns = quotient_columns - column - 1u;
        uint16_t* const pivot_tail = &quotient[pivot_row][column + 1u];
        const uint16_t pivot_value = quotient[pivot_row][column];
        if (pivot_value != 1u)
        {
            const uint16_t inverse = TinyGF16Inverse(pivot_value);
            if (inverse == 0u) {
                return false;
            }
            quotient[pivot_row][column] = 1u;
            TinyGF16ScaleCoefficients(pivot_tail, inverse, tail_columns);
            TinyGF16ScaleRow(
                rhs.data() + (size_t)pivot_row * block_bytes,
                inverse, block_bytes);
            ++stats.BlockMulAdds;
        }
        for (uint32_t row = 0u; row < H; ++row)
        {
            if (row == pivot_row) {
                continue;
            }
            const uint16_t scale = quotient[row][column];
            if (scale == 0u) {
                continue;
            }
            quotient[row][column] = 0u;
            TinyGF16MulAddCoefficients(
                &quotient[row][column + 1u], pivot_tail, scale,
                tail_columns);
            TinyGF16MulAddRow(
                rhs.data() + (size_t)row * block_bytes,
                rhs.data() + (size_t)pivot_row * block_bytes,
                scale, block_bytes);
            if (scale == 1u) {
                ++stats.BlockXors;
            }
            else {
                ++stats.BlockMulAdds;
            }
        }
    }

    // The unused completion rows carry exact left-nullspace syndromes: any
    // nonzero RHS is an inconsistent payload, matching the general path's
    // dependency checks in both its full-rank and deficient branches.
    for (uint32_t row = 0u; row < H; ++row)
    {
        if (row_used[row]) {
            continue;
        }
        if (!RowIsZero(
                rhs.data() + (size_t)row * block_bytes, block_bytes))
        {
            stats.ResidualRows += row + 1u;
            result_out = Wirehair_Error;
            return true;
        }
    }
    stats.ResidualRows += H;
    stats.ResidualRank = binary_rank + rank;
    if (rank < quotient_columns)
    {
        result_out = Wirehair_NeedMore;
        return true;
    }

    // Unique solution: publish the free values, then rebuild every binary
    // pivot from its preserved GF(2) relation.
    for (uint32_t i = 0u; i < quotient_columns; ++i)
    {
        const uint8_t pivot_row = pivot_row_of[i];
        if (pivot_row == UINT8_MAX) {
            return false;
        }
        std::memcpy(
            values.data() +
                (size_t)inactive_columns[free_columns[i]] * block_bytes,
            rhs.data() + (size_t)pivot_row * block_bytes,
            block_bytes);
        ++stats.BlockCopies;
    }
    for (uint32_t pivot = 0u; pivot < inactive_count; ++pivot)
    {
        if (!binary_have_pivot[pivot]) {
            continue;
        }
        InitializeMixedBinaryPivotValue(
            values.data() +
                (size_t)inactive_columns[pivot] * block_bytes,
            block_bytes,
            binary_pivot_rhs.data() + (size_t)pivot * block_bytes,
            binary_pivot_coeff.data() + (size_t)pivot * projection_words,
            free_columns, inactive_columns, values, stats);
    }
    result_out = Wirehair_Success;
    return true;
}

WirehairResult SolveMixedCompletionQuotient(
    const PrecodeSystem& system,
    uint32_t column_count,
    uint32_t inactive_count,
    uint32_t projection_words,
    uint32_t block_bytes,
    const std::vector<uint32_t>& inactive_index,
    const std::vector<uint32_t>& inactive_columns,
    const std::vector<uint64_t>& projection,
    const std::vector<uint64_t>& binary_pivot_coeff,
    const std::vector<uint8_t>& binary_pivot_rhs,
    const std::vector<uint8_t>& binary_have_pivot,
    uint32_t binary_rank,
    std::vector<uint8_t>& values,
    PrecodeSolveStats& stats)
{
    const uint32_t subfield_rows = ActiveMixedGF256Rows();
    const uint32_t extension_rows = ActiveMixedGF16Rows();
    const uint32_t H = subfield_rows + extension_rows;
    const uint32_t grouped_gf256_rows = ActiveMixedGroupedGF256Rows();
    const bool independent_extension_residues =
        ActiveMixedIndependentExtensionResidues();
    const uint64_t projection_elements =
        (uint64_t)column_count * projection_words;
    const uint64_t binary_coefficient_words =
        (uint64_t)inactive_count * projection_words;
    const uint64_t binary_rhs_bytes =
        (uint64_t)inactive_count * block_bytes;
    if (system.Params.Field != CompletionField::MixedGF256GF16 ||
        system.Params.HeavyRows != H || column_count < H ||
        grouped_gf256_rows > subfield_rows ||
        (grouped_gf256_rows != 0u &&
         independent_extension_residues) ||
        (block_bytes & 1u) != 0u ||
        binary_rank > inactive_count ||
        projection_words != PackedWordCount(inactive_count) ||
        inactive_index.size() != column_count ||
        inactive_columns.size() != inactive_count ||
        binary_coefficient_words > binary_pivot_coeff.max_size() ||
        binary_pivot_coeff.size() != (size_t)binary_coefficient_words ||
        binary_rhs_bytes > binary_pivot_rhs.max_size() ||
        binary_pivot_rhs.size() != (size_t)binary_rhs_bytes ||
        binary_have_pivot.size() != inactive_count ||
        projection_elements > projection.max_size() ||
        projection.size() != (size_t)projection_elements)
    {
        return Wirehair_InvalidInput;
    }
    const uint32_t first_grouped_gf256_row =
        subfield_rows - grouped_gf256_rows;
    const uint32_t first_heavy_column = column_count - H;
    for (uint32_t index = 0; index < inactive_count; ++index)
    {
        const uint32_t column = inactive_columns[index];
        if (column >= column_count || inactive_index[column] != index) {
            return Wirehair_InvalidInput;
        }
    }
    stats.BinaryResidualRank = binary_rank;
    stats.ResidualRank = binary_rank;
    const uint32_t quotient_columns = inactive_count - binary_rank;
    if (quotient_columns > H) {
        return Wirehair_NeedMore;
    }
    std::vector<uint32_t> free_columns;
    free_columns.reserve(quotient_columns);
    for (uint32_t column = 0; column < inactive_count; ++column) {
        if (!binary_have_pivot[column]) free_columns.push_back(column);
    }
    if (free_columns.size() != quotient_columns) return Wirehair_Error;

    // Tiny-K fast path: identical equations solved by one dense scalar
    // elimination, skipping the residue-bucket/factorization machinery
    // (which is then never constructed).  The projection-oracle test hook
    // compares that machinery against its dense oracle, so it pins the
    // general path.
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    const bool tiny_fast_path_blocked =
        MixedProjectionOracleUsers.load(std::memory_order_relaxed) != 0u;
#else
    const bool tiny_fast_path_blocked = false;
#endif
    if (!tiny_fast_path_blocked &&
        UseTinyMixedCompletionFastPath(
            system.Params.BlockCount, block_bytes))
    {
        WirehairResult tiny_result = Wirehair_Error;
        if (TrySolveMixedCompletionQuotientTiny(
                column_count, inactive_count, quotient_columns,
                subfield_rows, extension_rows, grouped_gf256_rows,
                independent_extension_residues,
                projection_words, block_bytes,
                inactive_index, inactive_columns, projection,
                binary_pivot_coeff, binary_pivot_rhs, binary_have_pivot,
                binary_rank, free_columns, values, stats, tiny_result))
        {
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
            ++stats.TinyMixedFastPathAcceptances;
#endif
            return tiny_result;
        }
        // Structurally unsupported shapes fall through with no state built.
    }

    // Every mixed solve previously built the same complete coefficient period.
    // Share immutable row-major and packed representations across all sizes.
    const MixedCoefficientRows* cached_rows = GetMixedCoefficientRows();
    const MixedPackedCoefficients* cached_packed =
        GetMixedPackedCoefficients();
    if (!cached_rows || !cached_packed) {
        return Wirehair_Error;
    }

    const uint32_t coefficient_period = ActiveMixedCoefficientPeriod();
    const uint8_t* grouped_block_shifts = nullptr;
    uint32_t grouped_block_shift_count = 0u;
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    if (grouped_gf256_rows != 0u)
    {
        // Reuse this thread's bounded (at most 5,041-byte) allocation across
        // decoder instances.  Grouped schedules are test-hook experiments,
        // so production and r=0 solves do not touch this TLS cache.
        static thread_local std::vector<uint8_t> grouped_block_shift_cache;
        if (!CacheMixedGroupedGF256ResidueBlockShifts(
                first_heavy_column, coefficient_period,
                grouped_block_shift_cache))
        {
            return Wirehair_Error;
        }
        grouped_block_shifts = grouped_block_shift_cache.data();
        grouped_block_shift_count =
            (uint32_t)grouped_block_shift_cache.size();
    }
#endif

    const uint32_t packed_words = ActiveMixedPackedCoefficientWords();
    static_assert(
        kMixedPackedCoefficientWords <= 4u,
        "mixed completion packing changed unexpectedly");
    if (packed_words == 0u ||
        packed_words > kMixedPackedCoefficientWords)
    {
        return Wirehair_Error;
    }
    if ((uint64_t)inactive_count * packed_words >
            (uint64_t)std::numeric_limits<size_t>::max() /
                sizeof(uint64_t))
    {
        return Wirehair_OOM;
    }
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    const bool check_projection_oracle =
        MixedProjectionOracleUsers.load(std::memory_order_relaxed) != 0u;
    std::vector<uint64_t> dense_projected;
    if (check_projection_oracle &&
        !ProjectMixedCompletionCoefficientsByDenseExpansion(
            column_count, inactive_count, projection_words,
            inactive_index, projection, *cached_packed, dense_projected))
    {
        return Wirehair_Error;
    }
#endif
    std::vector<uint64_t> projected;
    if (!ProjectMixedCompletionCoefficientsByResidueBuckets(
            column_count, inactive_count, projection_words,
            inactive_index, projection, *cached_packed,
            grouped_block_shifts, grouped_block_shift_count, projected) ||
        projected.size() != (size_t)inactive_count * packed_words)
    {
        return Wirehair_Error;
    }
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    if (check_projection_oracle)
    {
        if (dense_projected != projected) {
            return Wirehair_Error;
        }
        MixedProjectionOracleComparisons.fetch_add(
            1u, std::memory_order_relaxed);
    }
#endif

    if ((uint64_t)H * quotient_columns >
            (uint64_t)std::numeric_limits<size_t>::max() /
                sizeof(uint16_t))
    {
        return Wirehair_OOM;
    }
    std::vector<uint16_t> quotient_rows(
        (size_t)H * quotient_columns, 0u);
    // Only the free columns survive binary elimination.  Computing those
    // coefficients directly avoids materializing and Gauss-Jordan updating
    // H full inactive-width rows.  Substituting
    //
    //   x[p] = rhs[p] + sum(relation[p,f] * x[f])
    //
    // gives C'[f] = C[f] + sum(C[p] * relation[p,f]).  Each relation bit is
    // binary, so the packed GF(2^16) coefficient vectors combine by XOR.
    for (uint32_t i = 0u; i < quotient_columns; ++i)
    {
        const uint32_t free_column = free_columns[i];
        uint64_t packed[kMixedPackedCoefficientWords] = {};
        const uint64_t* free_coefficients =
            projected.data() + (size_t)free_column * packed_words;
        for (uint32_t word = 0u; word < packed_words; ++word) {
            packed[word] = free_coefficients[word];
        }
        for (uint32_t pivot = 0u; pivot < inactive_count; ++pivot)
        {
            const uint64_t* relation =
                binary_pivot_coeff.data() +
                    (size_t)pivot * projection_words;
            if (!binary_have_pivot[pivot] ||
                (relation[free_column >> 6] &
                    (UINT64_C(1) << (free_column & 63u))) == 0u)
            {
                continue;
            }
            const uint64_t* pivot_coefficients =
                projected.data() + (size_t)pivot * packed_words;
            for (uint32_t word = 0u; word < packed_words; ++word) {
                packed[word] ^= pivot_coefficients[word];
            }
        }
        for (uint32_t row = 0u; row < H; ++row) {
            quotient_rows[(size_t)row * quotient_columns + i] =
                (uint16_t)(
                    packed[row >> 2] >> ((row & 3u) * 16u));
        }
    }

    MixedQuotientFactorization factor;
    if (!FactorMixedCompletionQuotient(
            quotient_rows, H, quotient_columns, factor))
    {
        return Wirehair_Error;
    }

    const uint32_t elements = block_bytes / 2u;
    const bool rotate_residues = ActiveMixedResiduesRotated();
    const uint32_t populated_residues =
        std::min(coefficient_period, column_count);

    if (factor.Rank < quotient_columns)
    {
        if (!BuildMixedQuotientTransform(
                H, quotient_columns, factor))
        {
            return Wirehair_Error;
        }
        // Preserve the historical distinction between an underdetermined,
        // consistent system (NeedMore) and corrupted/inconsistent payloads
        // (Error).  The dependent factor rows are exact left-nullspace
        // syndromes.  Evaluate only those rows instead of constructing and
        // reducing every completion RHS block.
        std::array<uint8_t, kMixedCompletionRowsMax> dependency_rows = {};
        uint32_t dependency_count = 0u;
        for (uint32_t row = 0u; row < H; ++row) {
            if (factor.Rows[row].PivotColumn == UINT8_MAX) {
                dependency_rows[dependency_count++] = (uint8_t)row;
            }
        }
        if (dependency_count != H - factor.Rank) {
            return Wirehair_Error;
        }
        if ((uint64_t)dependency_count * elements >
                (uint64_t)std::numeric_limits<size_t>::max())
        {
            return Wirehair_OOM;
        }
        std::vector<uint8_t> syndrome_low =
            MakeBlockStorage((size_t)dependency_count * elements);
        std::vector<uint8_t> syndrome_high =
            MakeBlockStorage((size_t)dependency_count * elements);
        std::vector<uint8_t> source_low = MakeBlockStorage(elements);
        std::vector<uint8_t> source_high = MakeBlockStorage(elements);
        std::vector<uint8_t> residue_bucket = MakeBlockStorage(block_bytes);
        ++stats.BlockZeroFills;

        const auto add_syndrome_source = [&](
            uint32_t residue,
            const uint8_t* source,
            uint32_t first_subfield_row,
            uint32_t last_subfield_row,
            bool include_extension) -> bool
        {
            if (!DeinterleavePayloadPlanes(
                    source, source_low.data(),
                    source_high.data(), block_bytes))
            {
                return false;
            }
            for (uint32_t d = 0u; d < dependency_count; ++d)
            {
                const uint32_t dependency_row = dependency_rows[d];
                const uint16_t* transform = factor.Transform.data() +
                    (size_t)dependency_row * kMixedCompletionRowsMax;
                uint16_t scale = 0u;
                for (uint32_t row = first_subfield_row;
                     row < last_subfield_row; ++row)
                {
                    scale ^= GF16MultiplyInitialized(
                        transform[row],
                        cached_rows->Subfield[row][residue]);
                }
                if (include_extension)
                {
                    for (uint32_t er = 0u; er < extension_rows; ++er) {
                        scale ^= GF16MultiplyInitialized(
                            transform[subfield_rows + er],
                            cached_rows->Extension[er][residue]);
                    }
                }
                if (scale == 0u) continue;
                uint8_t* low = syndrome_low.data() + (size_t)d * elements;
                uint8_t* high = syndrome_high.data() + (size_t)d * elements;
                if (scale == 1u)
                {
                    gf256_add_mem(low, source_low.data(), (int)elements);
                    gf256_add_mem(high, source_high.data(), (int)elements);
                    ++stats.BlockXors;
                }
                else if (!MulAddPayloadPlanes(
                        low, high, scale,
                        source_low.data(), source_high.data(), elements))
                {
                    return false;
                }
                else {
                    ++stats.BlockMulAdds;
                }
            }
            return true;
        };

        const auto accumulate_syndrome_schedule = [&](
            uint32_t column_limit,
            uint32_t schedule,
            uint32_t first_subfield_row,
            uint32_t last_subfield_row,
            bool include_extension) -> bool
        {
            for (uint32_t residue = 0u;
                 residue < populated_residues; ++residue)
            {
                std::fill(
                    residue_bucket.begin(), residue_bucket.end(), uint8_t{0});
                ++stats.BlockZeroFills;
                BatchedBlockXorAccumulator bucket_xor(
                    residue_bucket.data(), block_bytes);
                if (schedule == 0u && !rotate_residues)
                {
                    for (uint32_t column = residue;
                         column < column_limit;
                         column += coefficient_period)
                    {
                        if (inactive_index[column] != UINT32_MAX) continue;
                        bucket_xor.Add(
                            values.data() + (size_t)column * block_bytes);
                        ++stats.BlockXors;
                    }
                }
                else
                {
                    uint32_t block_index = 0u;
                    for (uint32_t block_base = 0u;
                         block_base < column_limit;
                         block_base += coefficient_period)
                    {
                        uint32_t block_shift =
                            ActiveMixedResidueBlockShift(block_index);
                        if (schedule == 1u) {
                            block_shift =
                                ActiveMixedExtensionResidueBlockShift(
                                    block_index);
                        }
                        else if (schedule == 2u) {
                            CAT_DEBUG_ASSERT(
                                block_index < grouped_block_shift_count);
                            block_shift = grouped_block_shifts[block_index];
                        }
                        ++block_index;
                        const uint32_t unshifted = residue >= block_shift ?
                            residue - block_shift :
                            residue + coefficient_period - block_shift;
                        const uint32_t column = block_base + unshifted;
                        if (column >= column_limit ||
                            inactive_index[column] != UINT32_MAX)
                        {
                            continue;
                        }
                        bucket_xor.Add(
                            values.data() + (size_t)column * block_bytes);
                        ++stats.BlockXors;
                    }
                }
                if (schedule == 2u)
                {
                    for (uint32_t column = first_heavy_column;
                         column < column_count; ++column)
                    {
                        if (inactive_index[column] != UINT32_MAX ||
                            ActiveMixedCoefficientResidue(column) != residue)
                        {
                            continue;
                        }
                        bucket_xor.Add(
                            values.data() +
                                (size_t)column * block_bytes);
                        ++stats.BlockXors;
                    }
                }
                bucket_xor.Flush();
                if (!add_syndrome_source(
                        residue, residue_bucket.data(),
                        first_subfield_row, last_subfield_row,
                        include_extension))
                {
                    return false;
                }
            }
            return true;
        };

        bool syndromes_accumulated = false;
        if (grouped_gf256_rows != 0u)
        {
            syndromes_accumulated =
                accumulate_syndrome_schedule(
                    column_count, 0u, 0u,
                    first_grouped_gf256_row, true) &&
                accumulate_syndrome_schedule(
                    first_heavy_column, 2u,
                    first_grouped_gf256_row, subfield_rows, false);
        }
        else
        {
            syndromes_accumulated = accumulate_syndrome_schedule(
                column_count, 0u, 0u, subfield_rows,
                !independent_extension_residues);
            if (syndromes_accumulated && independent_extension_residues) {
                syndromes_accumulated = accumulate_syndrome_schedule(
                    column_count, 1u, 0u, 0u, true);
            }
        }
        if (!syndromes_accumulated)
        {
            return Wirehair_Error;
        }

        for (uint32_t pivot = 0u; pivot < inactive_count; ++pivot)
        {
            if (!binary_have_pivot[pivot]) continue;
            const uint64_t* packed_scales =
                projected.data() + (size_t)pivot * packed_words;
            std::array<uint16_t, kMixedCompletionRowsMax> scales = {};
            bool have_scale = false;
            for (uint32_t d = 0u; d < dependency_count; ++d)
            {
                const uint16_t* transform = factor.Transform.data() +
                    (size_t)dependency_rows[d] *
                        kMixedCompletionRowsMax;
                uint16_t scale = 0u;
                for (uint32_t row = 0u; row < H; ++row)
                {
                    const uint16_t coefficient = (uint16_t)(
                        packed_scales[row >> 2] >>
                            ((row & 3u) * 16u));
                    scale ^= GF16MultiplyInitialized(
                        transform[row], coefficient);
                }
                scales[d] = scale;
                have_scale = have_scale || scale != 0u;
            }
            if (!have_scale) continue;
            if (!DeinterleavePayloadPlanes(
                    binary_pivot_rhs.data() +
                        (size_t)pivot * block_bytes,
                    source_low.data(), source_high.data(), block_bytes))
            {
                return Wirehair_Error;
            }
            for (uint32_t d = 0u; d < dependency_count; ++d)
            {
                const uint16_t scale = scales[d];
                if (scale == 0u) continue;
                uint8_t* low = syndrome_low.data() + (size_t)d * elements;
                uint8_t* high = syndrome_high.data() + (size_t)d * elements;
                if (scale == 1u)
                {
                    gf256_add_mem(low, source_low.data(), (int)elements);
                    gf256_add_mem(high, source_high.data(), (int)elements);
                    ++stats.BlockXors;
                }
                else if (!MulAddPayloadPlanes(
                        low, high, scale,
                        source_low.data(), source_high.data(), elements))
                {
                    return Wirehair_Error;
                }
                else {
                    ++stats.BlockMulAdds;
                }
            }
        }

        for (uint32_t d = 0u; d < dependency_count; ++d)
        {
            if (!RowIsZero(
                    syndrome_low.data() + (size_t)d * elements, elements) ||
                !RowIsZero(
                    syndrome_high.data() + (size_t)d * elements, elements))
            {
                stats.ResidualRows += (uint32_t)dependency_rows[d] + 1u;
                return Wirehair_Error;
            }
        }
        stats.ResidualRows += H;
        stats.ResidualRank = binary_rank + factor.Rank;
        return Wirehair_NeedMore;
    }

    std::vector<uint8_t> subfield_rhs =
        MakeBlockStorage((size_t)subfield_rows * block_bytes);
    // MakeBlockStorage, not a bare vector: at solve width 0, elements is 0, so a
    // plain vector is EMPTY and data() returns NULL.  The back-substitution
    // batcher below tracks "is a row pending?" as (pending_existing_low !=
    // nullptr), and every row address is rhs_low.data() + row*elements == NULL,
    // so the pending pointer never becomes non-null and the batched
    // MulAddPayloadPlanes2 branch NEVER RUNS.  That silently dropped 448 of 3956
    // BlockMulAdds (1.7%) at width 0 while every other counter matched -- exactly
    // the class of count-only error that would corrupt a calibrated cost model.
    // Verified: width 0 agrees with its declared dispatch width 1280; other
    // real widths may deliberately choose different algorithms and counters.
    std::vector<uint8_t> rhs_low = MakeBlockStorage((size_t)H * elements);
    std::vector<uint8_t> rhs_high = MakeBlockStorage((size_t)H * elements);
    std::vector<uint8_t> source_low = MakeBlockStorage(elements);
    std::vector<uint8_t> source_high = MakeBlockStorage(elements);
    // Block-sized scratch, zero-filled by the allocation itself.  The
    // planar buffers hold half a block per plane, so rhs_low with rhs_high
    // is one block per completion row and source_low with source_high is
    // one more block.
    stats.BlockZeroFills += subfield_rows + H + 1u;
    void* subfield_destinations[kMixedGF256RowsMax];
    uint8_t subfield_scales[kMixedGF256RowsMax];
    for (uint32_t row = 0; row < subfield_rows; ++row) {
        subfield_destinations[row] =
            subfield_rhs.data() + (size_t)row * block_bytes;
    }
    const bool secondary_schedule =
        independent_extension_residues || grouped_gf256_rows != 0u;
    // Every residue-bucket route below is selected by payload width, so all
    // of them read the DISPATCH width.  The bucket planes themselves are
    // still allocated and addressed at the real width -- a payload-free solve
    // walks the same bucket schedule on planes that are zero bytes wide.
    const uint32_t dispatch_block_bytes = DispatchBlockBytes(block_bytes);
    const bool automatic_joint_residue_buckets =
        secondary_schedule &&
        UseAutomaticMixedJointResidueBuckets(
            system.Params.BlockCount, dispatch_block_bytes,
            coefficient_period);
    bool request_joint_residue_buckets =
        automatic_joint_residue_buckets;
    bool use_dual_residue_buckets = false;
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    // Dispatch-only quantity: the cache-residency estimate the dual route is
    // gated on, never an allocation size.
    const uint64_t one_bucket_plane_bytes =
        (uint64_t)coefficient_period * dispatch_block_bytes;
    const bool dual_buckets_fit =
        one_bucket_plane_bytes <=
            kMixedJointResidueBucketDataByteCap / 2u;
    const MixedResidueBucketMode bucket_mode =
        ActiveMixedResidueBucketModeForTesting();
    request_joint_residue_buckets =
        bucket_mode == MixedResidueBucketMode::JointDelta ||
        (bucket_mode == MixedResidueBucketMode::Automatic &&
         automatic_joint_residue_buckets);
    // A sequential dual-bucket scan improves bandwidth-bound large decoders,
    // but its accumulator setup hurts small systems and its random writes
    // hurt once the two bucket sets no longer fit comfortably in cache.
    use_dual_residue_buckets =
        secondary_schedule &&
        // The bucketed completion routines mul-add Extension[0] and
        // Extension[1] unconditionally and reject extension_rows < 2 outright,
        // and the caller turns that rejection into Wirehair_Error with no
        // fallback.  With the GF(2^16) rows removed the streamed path below is
        // the only one that handles a pure GF(256) band, so keep buckets off
        // there instead of failing every solve wide enough to select them.
        extension_rows >= 2u &&
        dual_buckets_fit &&
        (bucket_mode == MixedResidueBucketMode::Dual ||
         (bucket_mode == MixedResidueBucketMode::Automatic &&
          !automatic_joint_residue_buckets &&
          column_count >= 30000u &&
          dispatch_block_bytes >= 1024u &&
          one_bucket_plane_bytes * 2u <= (UINT64_C(128) << 10)));
#endif
    const bool use_joint_residue_buckets =
        secondary_schedule &&
        request_joint_residue_buckets &&
        // Same two-extension-row assumption as the dual route above.
        extension_rows >= 2u &&
        // Scratch-budget gate asked at the DISPATCH width.  A payload-free
        // solve carries no bucket bytes at all, so the real byte cap could
        // neither admit nor exclude it: the predicate rejects a zero width
        // outright, which would silently drop the joint route and fall back
        // to the streamed path, reporting a different BlockXors count and
        // zero for every MixedJoint* counter.  Asking the question at the
        // dispatch width answers it the way a real payload of that size
        // answers it, instead of exempting zero from the budget entirely.
        // The predicate itself keeps rejecting zero for the encoder, which
        // always has bytes to move.
        MixedJointResidueBucketStorageFits(
            coefficient_period, dispatch_block_bytes,
            kMixedJointResidueBucketDataByteCap);
    // Keep the streamed fallback unallocated on joint/dual paths.  Besides
    // avoiding one block-sized zero-fill, this makes their dominant scratch
    // exactly the selected bucket data planes plus scheduling metadata.
    std::vector<uint8_t> residue_bucket;
    const auto accumulate_extension_rhs = [&](
        uint32_t residue,
        const uint8_t* bucket) -> bool
    {
        // Removing the GF(2^16) rows degenerates the completion to a pure
        // GF(256) band: there is no extension RHS to accumulate at all.  This
        // early-out is load-bearing, not an optimization -- rhs_low/rhs_high
        // are sized H = subfield_rows + extension_rows, so with zero extension
        // rows the pair kernel below addresses rows subfield_rows and
        // subfield_rows+1 past the end of both.  The static_assert only pins
        // the compile-time production constant, so it never catches a runtime
        // row count of 0 or 1.
        if (extension_rows < 2u)
        {
            if (extension_rows == 0u)
            {
                return true;
            }
            if (!DeinterleavePayloadPlanes(
                    bucket, source_low.data(), source_high.data(),
                    block_bytes))
            {
                return false;
            }
            if (!MulAddPayloadPlanes(
                    rhs_low.data() + (size_t)subfield_rows * elements,
                    rhs_high.data() + (size_t)subfield_rows * elements,
                    cached_rows->Extension[0][residue],
                    source_low.data(), source_high.data(), elements))
            {
                return false;
            }
            stats.BlockMulAdds += extension_rows;
            return true;
        }
        if (!DeinterleavePayloadPlanes(
                bucket, source_low.data(), source_high.data(), block_bytes))
        {
            return false;
        }
        static_assert(
            kMixedGF16Rows >= 2u,
            "mixed completion pair kernel requires two GF16 rows");
        const uint32_t row0 = subfield_rows;
        const uint32_t row1 = row0 + 1u;
        if (!MulAddPayloadPlanes2(
                rhs_low.data() + (size_t)row0 * elements,
                rhs_high.data() + (size_t)row0 * elements,
                cached_rows->Extension[0][residue],
                rhs_low.data() + (size_t)row1 * elements,
                rhs_high.data() + (size_t)row1 * elements,
                cached_rows->Extension[1][residue],
                source_low.data(), source_high.data(), elements))
        {
            return false;
        }
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
        for (uint32_t er = 2u; er < extension_rows; ++er)
        {
            const uint32_t row = subfield_rows + er;
            if (!MulAddPayloadPlanes(
                    rhs_low.data() + (size_t)row * elements,
                    rhs_high.data() + (size_t)row * elements,
                    cached_rows->Extension[er][residue],
                    source_low.data(), source_high.data(), elements))
            {
                return false;
            }
        }
#endif
        stats.BlockMulAdds += extension_rows;
        return true;
    };
#if defined(__GNUC__) || defined(__clang__)
    if (__builtin_expect(
            use_dual_residue_buckets || use_joint_residue_buckets,
            false))
#else
    if (use_dual_residue_buckets || use_joint_residue_buckets)
#endif
    {
        bool accumulated = false;
        if (grouped_gf256_rows != 0u)
        {
            accumulated = use_joint_residue_buckets ?
                AccumulateJointGroupedMixedCompletionRhs(
                    column_count, first_heavy_column,
                    coefficient_period, populated_residues,
                    subfield_rows, first_grouped_gf256_row,
                    block_bytes, elements, extension_rows,
                    grouped_block_shifts, grouped_block_shift_count,
                    inactive_index, values, *cached_rows,
                    subfield_destinations,
                    rhs_low, rhs_high, source_low, source_high, stats) :
                AccumulateDualGroupedMixedCompletionRhs(
                    column_count, first_heavy_column,
                    coefficient_period, populated_residues,
                    subfield_rows, first_grouped_gf256_row,
                    block_bytes, elements, extension_rows,
                    grouped_block_shifts, grouped_block_shift_count,
                    inactive_index, values, *cached_rows,
                    subfield_destinations,
                    rhs_low, rhs_high, source_low, source_high, stats);
        }
        else
        {
            accumulated = use_joint_residue_buckets ?
                AccumulateJointMixedCompletionRhs(
                    column_count, coefficient_period, populated_residues,
                    subfield_rows, block_bytes, elements, extension_rows,
                    inactive_index, values, *cached_rows,
                    subfield_destinations,
                    rhs_low, rhs_high, source_low, source_high, stats) :
                AccumulateDualMixedCompletionRhs(
                    column_count, coefficient_period, populated_residues,
                    subfield_rows, block_bytes, elements, extension_rows,
                    inactive_index, values, *cached_rows,
                    subfield_destinations,
                    rhs_low, rhs_high, source_low, source_high, stats);
        }
        if (!accumulated)
        {
            return Wirehair_Error;
        }
        goto mixed_rhs_accumulated;
    }
    // MakeBlockStorage() keeps data() non-null at zero solve width, where
    // this bucket is still handed to the XOR kernels with a zero length.
    residue_bucket = MakeBlockStorage(block_bytes);
    ++stats.BlockZeroFills;
    for (uint32_t residue = 0u; residue < populated_residues; ++residue)
    {
        std::fill(
            residue_bucket.begin(), residue_bucket.end(), uint8_t{0});
        ++stats.BlockZeroFills;
        BatchedBlockXorAccumulator bucket_xor(
            residue_bucket.data(), block_bytes);
        if (!rotate_residues)
        {
            for (uint32_t column = residue;
                 column < column_count;
                 column += coefficient_period)
            {
                if (inactive_index[column] != UINT32_MAX) continue;
                bucket_xor.Add(
                    values.data() + (size_t)column * block_bytes);
                ++stats.BlockXors;
            }
        }
        else
        {
            uint32_t block_index = 0u;
            for (uint32_t block_base = 0u;
                 block_base < column_count;
                 block_base += coefficient_period)
            {
                const uint32_t block_shift =
                    ActiveMixedResidueBlockShift(block_index++);
                const uint32_t unshifted = residue >= block_shift ?
                    residue - block_shift :
                    residue + coefficient_period - block_shift;
                const uint32_t column = block_base + unshifted;
                if (column >= column_count) continue;
                if (inactive_index[column] != UINT32_MAX) continue;
                bucket_xor.Add(
                    values.data() + (size_t)column * block_bytes);
                ++stats.BlockXors;
            }
        }
        bucket_xor.Flush();
        const uint32_t primary_subfield_rows =
            grouped_gf256_rows != 0u ?
                first_grouped_gf256_row : subfield_rows;
        for (uint32_t row = 0u; row < primary_subfield_rows; ++row) {
            subfield_scales[row] = cached_rows->Subfield[row][residue];
        }
        AddScaledBlocks(
            subfield_destinations, subfield_scales,
            primary_subfield_rows,
            residue_bucket.data(), block_bytes, stats);
        if (!independent_extension_residues &&
            !accumulate_extension_rhs(residue, residue_bucket.data()))
        {
            return Wirehair_Error;
        }
    }
    if (independent_extension_residues)
    {
        for (uint32_t residue = 0u;
             residue < populated_residues; ++residue)
        {
            std::fill(
                residue_bucket.begin(), residue_bucket.end(), uint8_t{0});
            ++stats.BlockZeroFills;
            BatchedBlockXorAccumulator bucket_xor(
                residue_bucket.data(), block_bytes);
            uint32_t block_index = 0u;
            for (uint32_t block_base = 0u;
                 block_base < column_count;
                 block_base += coefficient_period)
            {
                const uint32_t block_shift =
                    ActiveMixedExtensionResidueBlockShift(block_index++);
                const uint32_t unshifted = residue >= block_shift ?
                    residue - block_shift :
                    residue + coefficient_period - block_shift;
                const uint32_t column = block_base + unshifted;
                if (column >= column_count) continue;
                if (inactive_index[column] != UINT32_MAX) continue;
                bucket_xor.Add(
                    values.data() + (size_t)column * block_bytes);
                ++stats.BlockXors;
            }
            bucket_xor.Flush();
            if (!accumulate_extension_rhs(
                    residue, residue_bucket.data()))
            {
                return Wirehair_Error;
            }
        }
    }
    if (grouped_gf256_rows != 0u)
    {
        for (uint32_t residue = 0u;
             residue < populated_residues; ++residue)
        {
            std::fill(
                residue_bucket.begin(), residue_bucket.end(), uint8_t{0});
            ++stats.BlockZeroFills;
            BatchedBlockXorAccumulator bucket_xor(
                residue_bucket.data(), block_bytes);
            uint32_t block_index = 0u;
            for (uint32_t block_base = 0u;
                 block_base < first_heavy_column;
                 block_base += coefficient_period)
            {
                CAT_DEBUG_ASSERT(block_index < grouped_block_shift_count);
                const uint32_t block_shift =
                    grouped_block_shifts[block_index++];
                const uint32_t unshifted = residue >= block_shift ?
                    residue - block_shift :
                    residue + coefficient_period - block_shift;
                const uint32_t column = block_base + unshifted;
                if (column >= first_heavy_column ||
                    inactive_index[column] != UINT32_MAX)
                {
                    continue;
                }
                bucket_xor.Add(
                    values.data() + (size_t)column * block_bytes);
                ++stats.BlockXors;
            }
            for (uint32_t column = first_heavy_column;
                 column < column_count; ++column)
            {
                if (inactive_index[column] != UINT32_MAX ||
                    ActiveMixedCoefficientResidue(column) != residue)
                {
                    continue;
                }
                bucket_xor.Add(
                    values.data() + (size_t)column * block_bytes);
                ++stats.BlockXors;
            }
            bucket_xor.Flush();
            for (uint32_t row = 0u; row < grouped_gf256_rows; ++row) {
                subfield_scales[row] = cached_rows->Subfield[
                    first_grouped_gf256_row + row][residue];
            }
            AddScaledBlocks(
                subfield_destinations + first_grouped_gf256_row,
                subfield_scales, grouped_gf256_rows,
                residue_bucket.data(), block_bytes, stats);
        }
    }
mixed_rhs_accumulated:

    // Reduce all completion RHS blocks through each binary pivot together.
    // Gauss-Jordan left the binary pivot matrix in reduced form, so every
    // other pivot column is zero in this relation.  A completion row's scale
    // at `pivot` is therefore still its original projected coefficient.  The
    // subfield rows stay in GF(256); the extension rows share a
    // single pivot-RHS deinterleave.
    for (uint32_t pivot = 0; pivot < inactive_count; ++pivot)
    {
        if (!binary_have_pivot[pivot]) continue;
        const uint64_t* packed_scales =
            projected.data() + (size_t)pivot * packed_words;
        for (uint32_t row = 0; row < subfield_rows; ++row)
        {
            subfield_scales[row] = (uint8_t)(
                packed_scales[row >> 2] >> ((row & 3u) * 16u));
        }
        bool have_subfield_scale = false;
        for (uint32_t row = 0; row < subfield_rows; ++row) {
            have_subfield_scale =
                have_subfield_scale || subfield_scales[row] != 0u;
        }
        if (have_subfield_scale) {
            AddScaledBlocks(
                subfield_destinations, subfield_scales,
                subfield_rows,
                binary_pivot_rhs.data() + (size_t)pivot * block_bytes,
                block_bytes, stats);
        }

        uint16_t extension_scales[kMixedGF16RowsMax] = {};
        bool have_extension_scale = false;
        for (uint32_t er = 0; er < extension_rows; ++er)
        {
            const uint32_t row = subfield_rows + er;
            const uint16_t scale = (uint16_t)(
                packed_scales[row >> 2] >> ((row & 3u) * 16u));
            extension_scales[er] = scale;
            have_extension_scale = have_extension_scale || scale != 0u;
        }
        if (have_extension_scale)
        {
            if (!DeinterleavePayloadPlanes(
                    binary_pivot_rhs.data() +
                        (size_t)pivot * block_bytes,
                    source_low.data(), source_high.data(), block_bytes))
            {
                return Wirehair_Error;
            }
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
            if (extension_rows == 1u)
            {
                if (!MulAddPayloadPlanes(
                        rhs_low.data() +
                            (size_t)subfield_rows * elements,
                        rhs_high.data() +
                            (size_t)subfield_rows * elements,
                        extension_scales[0],
                        source_low.data(), source_high.data(), elements))
                {
                    return Wirehair_Error;
                }
                ++stats.BlockMulAdds;
                continue;
            }
#endif
            if (!MulAddPayloadPlanes2(
                    rhs_low.data() + (size_t)subfield_rows * elements,
                    rhs_high.data() + (size_t)subfield_rows * elements,
                    extension_scales[0],
                    rhs_low.data() +
                        (size_t)(subfield_rows + 1u) * elements,
                    rhs_high.data() +
                        (size_t)(subfield_rows + 1u) * elements,
                    extension_scales[1],
                    source_low.data(), source_high.data(), elements))
            {
                return Wirehair_Error;
            }
            stats.BlockMulAdds +=
                (extension_scales[0] != 0u ? 1u : 0u) +
                (extension_scales[1] != 0u ? 1u : 0u);
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
            for (uint32_t er = 2u; er < extension_rows; ++er)
            {
                const uint16_t scale = extension_scales[er];
                if (scale == 0u) continue;
                const uint32_t row = subfield_rows + er;
                if (!MulAddPayloadPlanes(
                        rhs_low.data() + (size_t)row * elements,
                        rhs_high.data() + (size_t)row * elements,
                        scale, source_low.data(), source_high.data(),
                        elements))
                {
                    return Wirehair_Error;
                }
                ++stats.BlockMulAdds;
            }
#endif
        }
    }

    // One block-sized pass per subfield row that no multiply-add absorbs:
    // it reads the planar RHS and writes the interleaved block, the same
    // traffic as a block copy.
    for (uint32_t row = 0; row < subfield_rows; ++row) {
        if (!DeinterleavePayloadPlanes(
                subfield_rhs.data() + (size_t)row * block_bytes,
                rhs_low.data() + (size_t)row * elements,
                rhs_high.data() + (size_t)row * elements,
                block_bytes))
        {
            return Wirehair_Error;
        }
        ++stats.BlockCopies;
    }

    std::vector<uint8_t> scale_scratch = MakeBlockStorage(elements);
    for (uint32_t row = 0u; row < H; ++row)
    {
        const MixedQuotientRowPlan& plan = factor.Rows[row];
        uint8_t* low = rhs_low.data() + (size_t)row * elements;
        uint8_t* high = rhs_high.data() + (size_t)row * elements;
        for (uint32_t column = 0u;
             column < quotient_columns; ++column)
        {
            const uint16_t scale = plan.ForwardScales[column];
            if (scale == 0u) continue;
            const uint8_t pivot_row = factor.PivotRows[column];
            if (pivot_row == UINT8_MAX || pivot_row == row) {
                return Wirehair_Error;
            }
            const uint8_t* pivot_low =
                rhs_low.data() + (size_t)pivot_row * elements;
            const uint8_t* pivot_high =
                rhs_high.data() + (size_t)pivot_row * elements;
            if (scale == 1u)
            {
                gf256_add_mem(
                    low, pivot_low, (int)elements);
                gf256_add_mem(high, pivot_high, (int)elements);
                ++stats.BlockXors;
            }
            else if (!MulAddPayloadPlanes(
                    low, high, scale,
                    pivot_low, pivot_high, elements))
            {
                return Wirehair_Error;
            }
            else {
                ++stats.BlockMulAdds;
            }
        }

        if (plan.PivotColumn == UINT8_MAX)
        {
            if (!RowIsZero(low, elements) || !RowIsZero(high, elements))
            {
                stats.ResidualRows += row + 1u;
                return Wirehair_Error;
            }
            continue;
        }
        if (plan.PivotColumn >= quotient_columns ||
            factor.PivotRows[plan.PivotColumn] != row)
        {
            return Wirehair_Error;
        }
        if (plan.NormalizeScale != 1u)
        {
            if (!ScalePayloadPlanes(
                    low, high, plan.NormalizeScale,
                    scale_scratch.data(), elements))
            {
                return Wirehair_Error;
            }
            ++stats.BlockMulAdds;
        }
        uint8_t* pending_existing_low = nullptr;
        uint8_t* pending_existing_high = nullptr;
        uint16_t pending_existing_scale = 0u;
        for (uint32_t existing = 0u;
             existing < quotient_columns; ++existing)
        {
            const uint16_t scale = plan.BackScales[existing];
            if (scale == 0u) continue;
            const uint8_t existing_row = factor.PivotRows[existing];
            if (existing_row == UINT8_MAX || existing_row == row) {
                return Wirehair_Error;
            }
            uint8_t* existing_low =
                rhs_low.data() + (size_t)existing_row * elements;
            uint8_t* existing_high =
                rhs_high.data() + (size_t)existing_row * elements;
            if (!pending_existing_low)
            {
                pending_existing_low = existing_low;
                pending_existing_high = existing_high;
                pending_existing_scale = scale;
            }
            else
            {
                if (!MulAddPayloadPlanes2(
                        pending_existing_low, pending_existing_high,
                        pending_existing_scale,
                        existing_low, existing_high, scale,
                        low, high, elements))
                {
                    return Wirehair_Error;
                }
                stats.BlockXors +=
                    (pending_existing_scale == 1u ? 1u : 0u) +
                    (scale == 1u ? 1u : 0u);
                stats.BlockMulAdds +=
                    (pending_existing_scale != 1u ? 1u : 0u) +
                    (scale != 1u ? 1u : 0u);
                pending_existing_low = nullptr;
            }
        }
        if (pending_existing_low)
        {
            if (pending_existing_scale == 1u) {
                gf256_add_mem(
                    pending_existing_low, low, (int)elements);
                gf256_add_mem(
                    pending_existing_high, high, (int)elements);
                ++stats.BlockXors;
            }
            else if (!MulAddPayloadPlanes(
                    pending_existing_low, pending_existing_high,
                    pending_existing_scale, low, high, elements))
            {
                return Wirehair_Error;
            }
            else {
                ++stats.BlockMulAdds;
            }
        }
    }
    stats.ResidualRows += H;
    stats.ResidualRank = binary_rank + factor.Rank;

    // Publishing a free value writes one block, exactly what the tiny
    // path's memcpy of the same value counts as a block copy.
    for (uint32_t i = 0u; i < quotient_columns; ++i)
    {
        const uint8_t pivot_row = factor.PivotRows[i];
        if (pivot_row == UINT8_MAX || !InterleavePayloadPlanes(
                rhs_low.data() + (size_t)pivot_row * elements,
                rhs_high.data() + (size_t)pivot_row * elements,
                values.data() +
                    (size_t)inactive_columns[free_columns[i]] * block_bytes,
                block_bytes))
        {
            return Wirehair_Error;
        }
        ++stats.BlockCopies;
    }
    for (uint32_t pivot = 0; pivot < inactive_count; ++pivot)
    {
        if (!binary_have_pivot[pivot]) continue;
        uint8_t* value = values.data() +
            (size_t)inactive_columns[pivot] * block_bytes;
        const uint64_t* relation =
            binary_pivot_coeff.data() +
                (size_t)pivot * projection_words;
        InitializeMixedBinaryPivotValue(
            value, block_bytes,
            binary_pivot_rhs.data() + (size_t)pivot * block_bytes,
            relation, free_columns, inactive_columns, values, stats);
    }
    return Wirehair_Success;
}

/**
    Reusable per-thread scratch for PeelBinaryRowsImplementation().

    Every buffer here is sized from the system and fully rewritten at the top
    of each peel, so retaining capacity across calls is observationally
    identical to fresh vectors while removing ten small allocations per solve.
    Measured on its own the saving is small -- glibc already serves these
    same-sized blocks from tcache, and the assign() that refills each buffer
    costs what the fresh allocation's zero-fill did -- so this is bookkeeping
    hygiene rather than the reason the peel got faster.

    The two adjacency buffers (column_offsets / column_rows) deliberately stay
    per-call vectors: PeelResult reports their allocation count and capacity
    bytes as the solver's binary-adjacency footprint, and pooling them would
    make that report a fiction.
*/
struct PeelScratch
{
    std::vector<PeelRowState> RowState;
    std::vector<uint8_t> Resolved;
    std::vector<uint32_t> Queue;
    std::vector<uint32_t> DegreeTwoRefs;
    std::vector<uint64_t> DegreeTwoHeap;
    std::vector<size_t> ReferenceBucketOffsets;
    std::vector<size_t> ReferenceBucketCursor;
    std::vector<uint32_t> ReferenceColumns;
    std::vector<uint32_t> DegreeTwoTieRank;
    std::vector<uint32_t> DegreeTwoRankColumn;
};

template<bool UseLowDegreeXor>
PeelScratch& GetPeelScratch()
{
    static thread_local PeelScratch scratch;
    return scratch;
}

template<bool UseLowDegreeXor>
PeelResult PeelBinaryRowsImplementation(
    uint32_t column_count,
    const BinaryEquationArena& rows)
{
    PeelResult out;
    out.SolveRow.assign(column_count, UINT32_MAX);
    out.UsedRows.assign(rows.size(), 0u);
    out.PeelOrder.reserve(column_count);

    PeelScratch& scratch = GetPeelScratch<UseLowDegreeXor>();
    std::vector<PeelRowState>& row_state = scratch.RowState;
    row_state.assign(rows.size(), PeelRowState());
    std::vector<size_t> column_offsets((size_t)column_count + 1u, 0u);
    std::vector<uint8_t>& resolved = scratch.Resolved;
    resolved.assign(column_count, 0u);
    std::vector<uint32_t>& queue = scratch.Queue;
    queue.clear();
    queue.reserve(rows.size());
    std::vector<uint32_t>& degree_two_refs = scratch.DegreeTwoRefs;
    degree_two_refs.assign(column_count, 0u);
    // A hand-rolled binary max-heap over the same 64-bit keys as the previous
    // std::priority_queue, so the pop/push order (and therefore the selection
    // sequence) is unchanged, but the backing store can be pooled.
    std::vector<uint64_t>& degree_two_heap = scratch.DegreeTwoHeap;
    degree_two_heap.clear();
    degree_two_heap.reserve(column_count);

    // LowDegreeXor is now the XOR of *every* live column in the row, not just
    // of the one or two live columns at low degree.  The two agree exactly at
    // degree one and degree two (the only degrees that read it), and keeping
    // the invariant globally costs one XOR against a word resolve() already
    // has in a register.
    for (uint32_t r = 0; r < (uint32_t)rows.size(); ++r)
    {
        CAT_DEBUG_ASSERT(rows[r].Columns.size() <= UINT16_MAX);
        row_state[r].Live = (uint16_t)rows[r].Columns.size();
        if (row_state[r].Live == 1u) {
            queue.push_back(r);
        }
        uint32_t live_xor = 0u;
        for (uint32_t column : rows[r].Columns) {
            CAT_DEBUG_ASSERT(column < column_count);
            ++column_offsets[(size_t)column + 1u];
            live_xor ^= column;
        }
        row_state[r].LowDegreeXor = (uint16_t)live_xor;
    }

    for (uint32_t column = 0; column < column_count; ++column) {
        column_offsets[(size_t)column + 1u] += column_offsets[column];
    }
    std::vector<uint32_t> column_rows(column_offsets[column_count]);
    out.AdjacencyStorageBytes =
        (uint64_t)column_offsets.capacity() * sizeof(size_t) +
        (uint64_t)column_rows.capacity() * sizeof(uint32_t);
    out.AdjacencyStorageAllocations =
        (column_offsets.capacity() != 0u ? 1u : 0u) +
        (column_rows.capacity() != 0u ? 1u : 0u);
    for (uint32_t r = 0; r < (uint32_t)rows.size(); ++r)
    {
        for (uint32_t column : rows[r].Columns)
        {
            const size_t destination =
                column_offsets[column] + degree_two_refs[column]++;
            column_rows[destination] = r;
        }
    }
    std::fill(degree_two_refs.begin(), degree_two_refs.end(), 0u);

    // The fallback priority never changes: Prefer the largest original
    // reference count, then the lowest column id.  Stable counting-sort
    // buckets preserve that exact order without a 12-byte heap node and
    // logarithmic insertion for every column.
    uint32_t max_reference_count = 0u;
    for (uint32_t column = 0; column < column_count; ++column) {
        max_reference_count = std::max(
            max_reference_count,
            (uint32_t)(column_offsets[(size_t)column + 1u] -
                column_offsets[column]));
    }
    std::vector<size_t>& reference_bucket_offsets =
        scratch.ReferenceBucketOffsets;
    reference_bucket_offsets.assign((size_t)max_reference_count + 2u, 0u);
    for (uint32_t column = 0; column < column_count; ++column)
    {
        const uint32_t references = (uint32_t)(
            column_offsets[(size_t)column + 1u] - column_offsets[column]);
        ++reference_bucket_offsets[(size_t)references + 1u];
    }
    for (uint32_t references = 0;
         references <= max_reference_count;
         ++references)
    {
        reference_bucket_offsets[(size_t)references + 1u] +=
            reference_bucket_offsets[references];
    }
    std::vector<size_t>& reference_bucket_cursor =
        scratch.ReferenceBucketCursor;
    reference_bucket_cursor.assign(
        reference_bucket_offsets.begin(), reference_bucket_offsets.end());
    std::vector<uint32_t>& reference_columns = scratch.ReferenceColumns;
    reference_columns.resize(column_count);
    for (uint32_t column = 0; column < column_count; ++column)
    {
        const uint32_t references = (uint32_t)(
            column_offsets[(size_t)column + 1u] - column_offsets[column]);
        reference_columns[reference_bucket_cursor[references]++] = column;
    }
    reference_bucket_cursor.assign(
        reference_bucket_offsets.begin(), reference_bucket_offsets.end());

    // Degree-two priorities compare a changing reference count followed by
    // the immutable total-reference count and reverse column id.  Replace
    // the three-field heap node with one 64-bit key.  The low word is a
    // precomputed rank that preserves the exact original tie order: larger
    // total-reference counts first, then lower column ids.
    std::vector<uint32_t>& degree_two_tie_rank = scratch.DegreeTwoTieRank;
    degree_two_tie_rank.resize(column_count);
    std::vector<uint32_t>& degree_two_rank_column =
        scratch.DegreeTwoRankColumn;
    degree_two_rank_column.resize(column_count);
    uint32_t next_tie_rank = 0u;
    for (uint32_t references = 0u;
         references <= max_reference_count;
         ++references)
    {
        const size_t begin = reference_bucket_offsets[references];
        const size_t end =
            reference_bucket_offsets[(size_t)references + 1u];
        for (size_t index = end; index > begin; --index)
        {
            const uint32_t column = reference_columns[index - 1u];
            degree_two_tie_rank[column] = next_tie_rank;
            degree_two_rank_column[next_tie_rank++] = column;
        }
    }
    CAT_DEBUG_ASSERT(next_tie_rank == column_count);

    const auto degree_two_key = [&](uint32_t column) {
        return (uint64_t)degree_two_refs[column] << 32 |
            degree_two_tie_rank[column];
    };

    // Used rows are selected only at live degree one.  Resolving their sole
    // column necessarily visits the selected row through this adjacency and
    // reduces Live to zero before resolve() returns.  Live therefore already
    // excludes them from every later queue/degree-two operation; UsedRows is
    // retained solely to classify the unused residual equations afterward.
    const auto add_degree_two = [&](uint32_t row) {
        PeelRowState& state = row_state[row];
        if (state.Live != 2u) {
            return;
        }
        // state.LowDegreeXor already holds the XOR of the two live columns, so
        // only the *first* of them has to be located.  The scan stops there
        // and recovers the second by XOR instead of walking the rest of the
        // row, and because the second live column is by construction later in
        // the row's column order the two heap keys are still pushed in the
        // original row order.  Every row reaches degree two exactly once, so
        // the old full-row scan cost one pass over the whole system per solve:
        // at K=100 it examined 1184 columns, this examines 339.
        const uint32_t pair_xor = state.LowDegreeXor;
        uint32_t first = UINT32_MAX;
        for (uint32_t column : rows[row].Columns)
        {
            ++out.RowScanSteps;
            if (!resolved[column]) {
                first = column;
                break;
            }
        }
        if (first == UINT32_MAX) {
            return;
        }
        const uint32_t second = pair_xor ^ first;
        CAT_DEBUG_ASSERT(
            second < column_count && second != first && !resolved[second]);
        ++degree_two_refs[first];
        degree_two_heap.push_back(degree_two_key(first));
        std::push_heap(degree_two_heap.begin(), degree_two_heap.end());
        ++degree_two_refs[second];
        degree_two_heap.push_back(degree_two_key(second));
        std::push_heap(degree_two_heap.begin(), degree_two_heap.end());
        out.HeapOperations += 2u;
    };
    const auto remove_degree_two = [&](
        uint32_t row,
        uint32_t resolved_column)
    {
        PeelRowState& state = row_state[row];
        if (state.Live != 2u) {
            return;
        }
        const uint32_t other =
            (uint32_t)state.LowDegreeXor ^ resolved_column;
        CAT_DEBUG_ASSERT(other < column_count && !resolved[other]);
        if (degree_two_refs[other] > 0u) {
            --degree_two_refs[other];
        }
        // If the lower degree remains nonzero, its matching key is still in
        // the lazy heap from the earlier increment.  It cannot have reached
        // the top while a higher key for this unresolved column existed, so
        // pushing it again here would only create a duplicate.  Degree zero
        // needs no live heap key.  resolve() reduces LowDegreeXor to `other`
        // for every degree, so this path no longer writes it.
    };
    for (uint32_t r = 0; r < (uint32_t)rows.size(); ++r) {
        add_degree_two(r);
    }

    const auto resolve = [&](uint32_t column) {
        resolved[column] = 1u;
        out.AdjacencyVisits += column_offsets[(size_t)column + 1u] -
            column_offsets[column];
        for (size_t ref = column_offsets[column];
             ref < column_offsets[(size_t)column + 1u];
             ++ref)
        {
            const uint32_t row = column_rows[ref];
            if (row_state[row].Live == 0u) {
                continue;
            }
            remove_degree_two(row, column);
            row_state[row].LowDegreeXor ^= (uint16_t)column;
            --row_state[row].Live;
            add_degree_two(row);
            if (row_state[row].Live == 1u) {
                queue.push_back(row);
            }
        }
    };

    uint32_t remaining = column_count;
    size_t queue_head = 0u;
    while (remaining > 0u)
    {
        while (queue_head < queue.size())
        {
            const uint32_t row = queue[queue_head++];
            if (row_state[row].Live != 1u) {
                continue;
            }
            uint32_t column = UINT32_MAX;
            if (UseLowDegreeXor) {
                column = row_state[row].LowDegreeXor;
                CAT_DEBUG_ASSERT(
                    column < column_count && !resolved[column]);
            }
            else
            {
                for (uint32_t candidate : rows[row].Columns)
                {
                    if (!resolved[candidate]) {
                        column = candidate;
                        break;
                    }
                }
                if (column == UINT32_MAX) {
                    continue;
                }
            }
            out.UsedRows[row] = 1u;
            out.SolveRow[column] = row;
            out.PeelOrder.push_back(column);
            resolve(column);
            CAT_DEBUG_ASSERT(row_state[row].Live == 0u);
            --remaining;
        }
        if (remaining == 0u) {
            break;
        }

        uint32_t best = UINT32_MAX;
        while (!degree_two_heap.empty())
        {
            const uint64_t candidate = degree_two_heap.front();
            const uint32_t column =
                degree_two_rank_column[(uint32_t)candidate];
            if (resolved[column] ||
                degree_two_refs[column] != (uint32_t)(candidate >> 32))
            {
                std::pop_heap(
                    degree_two_heap.begin(), degree_two_heap.end());
                degree_two_heap.pop_back();
                ++out.HeapOperations;
                continue;
            }
            best = column;
            break;
        }
        while (best == UINT32_MAX)
        {
            size_t& cursor = reference_bucket_cursor[max_reference_count];
            const size_t end =
                reference_bucket_offsets[(size_t)max_reference_count + 1u];
            while (cursor < end && resolved[reference_columns[cursor]]) {
                ++cursor;
                ++out.RowScanSteps;
            }
            if (cursor < end) {
                best = reference_columns[cursor];
                break;
            }
            if (max_reference_count == 0u) {
                break;
            }
            --max_reference_count;
        }
        if (best == UINT32_MAX) {
            break;
        }
        out.InactiveOrder.push_back(best);
        resolve(best);
        --remaining;
    }
    return out;
}

PeelResult PeelBinaryRows(
    uint32_t column_count,
    const BinaryEquationArena& rows)
{
    PeelResult out = PeelBinaryRowsImplementation<true>(column_count, rows);
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    if (BinaryPeelOracleUsers.load(std::memory_order_relaxed) != 0u)
    {
        const PeelResult reference =
            PeelBinaryRowsImplementation<false>(column_count, rows);
        std::vector<uint32_t> last_row(column_count, UINT32_MAX);
        bool duplicate_free = true;
        for (uint32_t row = 0u; row < (uint32_t)rows.size(); ++row)
        {
            for (uint32_t column : rows[row].Columns)
            {
                if (column >= column_count || last_row[column] == row) {
                    duplicate_free = false;
                    break;
                }
                last_row[column] = row;
            }
            if (!duplicate_free) {
                break;
            }
        }
        if (!duplicate_free || out.SolveRow != reference.SolveRow ||
            out.PeelOrder != reference.PeelOrder ||
            out.InactiveOrder != reference.InactiveOrder ||
            out.UsedRows != reference.UsedRows ||
            out.AdjacencyStorageBytes !=
                reference.AdjacencyStorageBytes ||
            out.AdjacencyStorageAllocations !=
                reference.AdjacencyStorageAllocations)
        {
            // Valid systems have at least two columns, so clearing both order
            // vectors turns an oracle disagreement into a terminal solve
            // error at the caller's existing completeness check.
            out.PeelOrder.clear();
            out.InactiveOrder.clear();
        }
        else {
            BinaryPeelOracleComparisons.fetch_add(
                1u, std::memory_order_relaxed);
        }
    }
#endif
    return out;
}

bool RowIsZero(const uint8_t* data, uint32_t bytes)
{
    for (uint32_t i = 0; i < bytes; ++i) {
        if (data[i] != 0u) {
            return false;
        }
    }
    return true;
}

} // namespace

uint32_t TinyMixedFastPathMaxSourceBlocks()
{
    return WIREHAIR_V2_TINY_MIXED_FASTPATH_MAX_SOURCE_BLOCKS;
}

uint32_t TinyMixedFastPathMaxBlockBytes()
{
    return WIREHAIR_V2_TINY_MIXED_FASTPATH_MAX_BLOCK_BYTES;
}

uint32_t TinyMixedFastPathMaxProductBytes()
{
    return WIREHAIR_V2_TINY_MIXED_FASTPATH_MAX_PRODUCT_BYTES;
}

bool IsCanonicalStableTargetPacketRowState()
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    return PacketRowSeedMultiplier == 1u &&
        !PacketRowSeedAvalanche &&
        OddPacketPeelSeedXor == 0u;
#else
    return true;
#endif
}

bool IsCanonicalPacketDegreeState()
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    const std::vector<double>* cdf = nullptr;
    return ActivePeelDegreeCdf(cdf) && cdf == nullptr;
#else
    return true;
#endif
}

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
uint64_t ActivePacketRowEquationStateFingerprintForTesting()
{
    PacketHookStateHasher hash;
    hash.U64(UINT64_C(0x574832504b465031)); // "WH2PKFP1"
    hash.U32(PacketRowSeedMultiplier);
    hash.U8(PacketRowSeedAvalanche ? 1u : 0u);
    hash.U32(OddPacketPeelSeedXor);

    const std::vector<double>* cdf = nullptr;
    const bool valid = ActivePeelDegreeCdf(cdf);
    hash.U8(valid ? 1u : 0u);
    hash.U8(cdf ? 1u : 0u);
    if (cdf) {
        hash.U64(cdf == &PeelDegreeOverride.Cdf ?
            PeelDegreeOverrideFingerprint :
            EnvironmentPeelDegreeCdfFingerprint());
    }
    return hash.Finish();
}

bool SetTinyMixedFastPathModeForTesting(int mode)
{
    if (mode < -1 || mode > 1) {
        return false;
    }
    TinyMixedFastPathTestModeValue = mode;
    return true;
}

int TinyMixedFastPathModeForTesting()
{
    return TinyMixedFastPathTestModeValue;
}

bool CheckTinyMixedScalarHelpersForTesting()
{
    if (!InitializeGF16()) {
        return false;
    }
    // Exhaustive inverse identity against the power-chain reference.
    if (TinyGF16Inverse(0u) != 0u) {
        return false;
    }
    for (uint32_t value = 1u; value < 65536u; ++value)
    {
        const uint16_t inverse = TinyGF16Inverse((uint16_t)value);
        if (inverse != GF16InverseInitialized((uint16_t)value) ||
            GF16MultiplyInitialized((uint16_t)value, inverse) != 1u)
        {
            return false;
        }
    }
    // Randomized row muladd/scale against the scalar reference kernel,
    // covering the zero/one scale specials and odd element counts.
    uint64_t state = UINT64_C(0x243F6A8885A308D3);
    const auto next = [&state]() {
        state += UINT64_C(0x9E3779B97F4A7C15);
        uint64_t z = state;
        z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
        return z ^ (z >> 31);
    };
    static const uint32_t kWidths[] = {2u, 4u, 6u, 30u, 62u, 64u, 258u};
    for (uint32_t width : kWidths)
    {
        for (uint32_t trial = 0u; trial < 64u; ++trial)
        {
            std::vector<uint8_t> source(width);
            std::vector<uint8_t> destination(width);
            for (uint32_t i = 0u; i < width; ++i) {
                source[i] = (uint8_t)next();
                destination[i] = (uint8_t)next();
            }
            uint16_t scale = (uint16_t)next();
            if (trial == 0u) scale = 0u;
            if (trial == 1u) scale = 1u;
            std::vector<uint8_t> expected = destination;
            if (scale != 0u &&
                !GF16MulAddMem(
                    expected.data(), scale, source.data(), width))
            {
                return false;
            }
            std::vector<uint8_t> actual = destination;
            TinyGF16MulAddRow(
                actual.data(), source.data(), scale, width);
            if (actual != expected) {
                return false;
            }
            std::vector<uint8_t> scaled = destination;
            if (!GF16ScaleMem(scaled.data(), scale, width)) {
                return false;
            }
            std::vector<uint8_t> tiny_scaled = destination;
            TinyGF16ScaleRow(tiny_scaled.data(), scale, width);
            if (tiny_scaled != scaled) {
                return false;
            }
        }
    }
    return true;
}
#endif

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
bool SetPacketRowSeedMultiplierForTesting(uint32_t multiplier)
{
    if (multiplier == 0u || (multiplier & 1u) == 0u) {
        return false;
    }
    PacketRowSeedMultiplier = multiplier;
    return true;
}

bool SetPeelDegreesForTesting(const std::vector<double>& weights)
{
    if (weights.empty())
    {
        // Empty has always been the explicit clear spelling.
        PeelDegreeOverride = PeelDegreeCdfState();
        PeelDegreeOverrideFingerprint =
            FingerprintPeelDegreeCdf(PeelDegreeOverride.Cdf);
        return true;
    }
    PeelDegreeCdfState candidate;
    candidate.IsSet = true;
    candidate.Valid = MakePeelDegreeCdf(weights, candidate.Cdf);
    PeelDegreeOverride = std::move(candidate);
    PeelDegreeOverrideFingerprint =
        FingerprintPeelDegreeCdf(PeelDegreeOverride.Cdf);
    return PeelDegreeOverride.Valid;
}

void ClearPeelDegreesForTesting()
{
    PeelDegreeOverride = PeelDegreeCdfState();
    PeelDegreeOverrideFingerprint =
        FingerprintPeelDegreeCdf(PeelDegreeOverride.Cdf);
}

bool PeelDegreeConfigurationValidForTesting()
{
    const std::vector<double>* cdf = nullptr;
    return ActivePeelDegreeCdf(cdf);
}

uint32_t SamplePeelDegreeForTesting(
    uint32_t random_value,
    uint16_t peel_column_count)
{
    const std::vector<double>* cdf = nullptr;
    return !ActivePeelDegreeCdf(cdf) || cdf == nullptr ?
        0u : SamplePeelDegreeOverride(
        *cdf, random_value, peel_column_count, 0u);
}

uint32_t PeelDegreeSampleComparisonCountForTesting()
{
    uint32_t comparisons = 0u;
    const std::vector<double> canonical(kPeelOverrideCdfSize, 1.0);
    (void)SamplePeelDegreeOverrideImpl<true>(
        canonical, 0u, 128u, 0u, &comparisons);
    return comparisons;
}

void SetPacketRowSeedAvalancheForTesting(bool enabled)
{
    PacketRowSeedAvalanche = enabled;
}

void SetOddPacketPeelSeedXorForTesting(uint32_t seed_xor)
{
    OddPacketPeelSeedXor = seed_xor;
}

void SetMixedNullWitnessDiagnosticForTesting(
    MixedNullWitnessDiagnostic* diagnostic)
{
    if (diagnostic) {
        *diagnostic = MixedNullWitnessDiagnostic{};
    }
    MixedNullWitnessSink = diagnostic;
}

void SetMixedProjectionOracleForTesting(bool enabled)
{
    if (enabled) {
        MixedProjectionOracleUsers.fetch_add(1u, std::memory_order_relaxed);
        return;
    }
    uint32_t users =
        MixedProjectionOracleUsers.load(std::memory_order_relaxed);
    while (users != 0u &&
           !MixedProjectionOracleUsers.compare_exchange_weak(
               users, users - 1u,
               std::memory_order_relaxed,
               std::memory_order_relaxed))
    {
    }
    CAT_DEBUG_ASSERT(users != 0u);
}

void ResetMixedProjectionOracleComparisonsForTesting()
{
    MixedProjectionOracleComparisons.store(0u, std::memory_order_relaxed);
}

uint64_t MixedProjectionOracleComparisonsForTesting()
{
    return MixedProjectionOracleComparisons.load(std::memory_order_relaxed);
}

void SetBinaryPeelOracleForTesting(bool enabled)
{
    if (enabled) {
        BinaryPeelOracleUsers.fetch_add(1u, std::memory_order_relaxed);
        return;
    }
    uint32_t users = BinaryPeelOracleUsers.load(std::memory_order_relaxed);
    while (users != 0u &&
           !BinaryPeelOracleUsers.compare_exchange_weak(
               users, users - 1u,
               std::memory_order_relaxed,
               std::memory_order_relaxed))
    {
    }
    CAT_DEBUG_ASSERT(users != 0u);
}

void ResetBinaryPeelOracleComparisonsForTesting()
{
    BinaryPeelOracleComparisons.store(0u, std::memory_order_relaxed);
}

uint64_t BinaryPeelOracleComparisonsForTesting()
{
    return BinaryPeelOracleComparisons.load(std::memory_order_relaxed);
}
#endif

void PrecodeSolveResumeState::Clear()
{
    PrecodeSolveResumeState empty;
    Swap(empty);
}

void PrecodeSolveResumeState::Swap(
    PrecodeSolveResumeState& other) noexcept
{
    using std::swap;
    swap(SourceCount, other.SourceCount);
    swap(PrecodeCount, other.PrecodeCount);
    swap(ColumnCount, other.ColumnCount);
    swap(BlockBytes, other.BlockBytes);
    swap(InactiveCount, other.InactiveCount);
    swap(ProjectionWords, other.ProjectionWords);
    swap(Rank, other.Rank);
    swap(Config, other.Config);
    swap(Runtime, other.Runtime);
    swap(Stats, other.Stats);
    InactiveIndex.swap(other.InactiveIndex);
    InactiveColumns.swap(other.InactiveColumns);
    Projection.swap(other.Projection);
    Values.swap(other.Values);
    PivotCoefficients.swap(other.PivotCoefficients);
    PivotRhs.swap(other.PivotRhs);
    HavePivot.swap(other.HavePivot);
    CoefficientScratch.swap(other.CoefficientScratch);
    RhsScratch.swap(other.RhsScratch);
    swap(Active, other.Active);
}

size_t PrecodeSolveResumeState::PersistentBytes() const
{
    size_t bytes = 0u;
    const auto add = [&bytes](size_t count, size_t width) {
        if (width != 0u && count >
            (std::numeric_limits<size_t>::max() - bytes) / width)
        {
            bytes = std::numeric_limits<size_t>::max();
        }
        else {
            bytes += count * width;
        }
    };
    add(InactiveIndex.capacity(), sizeof(uint32_t));
    add(InactiveColumns.capacity(), sizeof(uint32_t));
    add(Projection.capacity(), sizeof(uint64_t));
    add(Values.capacity(), sizeof(uint8_t));
    add(PivotCoefficients.capacity(), sizeof(uint8_t));
    add(PivotRhs.capacity(), sizeof(uint8_t));
    add(HavePivot.capacity(), sizeof(uint8_t));
    add(CoefficientScratch.capacity(), sizeof(uint8_t));
    add(RhsScratch.capacity(), sizeof(uint8_t));
    return bytes;
}

bool IsPacketRowDomainValid(
    uint32_t source_count,
    uint32_t precode_count,
    uint32_t mix_count)
{
    return source_count >= 2u && source_count <= 64000u &&
        precode_count >= kMinPacketPrecodeCount &&
        precode_count <= kMaxPacketPrecodeCount &&
        (uint64_t)source_count + precode_count <= UINT16_MAX &&
        mix_count >= 1u && mix_count <= kCertifiedPacketMixCount &&
        mix_count <= precode_count;
}

bool PacketRowRuntime::Initialize(
    uint32_t source_count,
    uint32_t precode_count,
    uint32_t mix_count)
{
    if (!IsPacketRowDomainValid(source_count, precode_count, mix_count))
    {
        *this = PacketRowRuntime();
        return false;
    }
    const uint16_t source_prime =
        wirehair::NextPrime16((uint16_t)source_count);
    const uint16_t precode_prime =
        wirehair::NextPrime16((uint16_t)precode_count);
    if (source_prime == 0u || precode_prime == 0u)
    {
        *this = PacketRowRuntime();
        return false;
    }
    SourceCount = source_count;
    PrecodeCount = precode_count;
    MixCount = mix_count;
    SourcePrimeValue = source_prime;
    PrecodePrimeValue = precode_prime;
    return true;
}

bool PacketRowRuntime::IsValidFor(
    uint32_t source_count,
    uint32_t precode_count,
    uint32_t mix_count) const
{
    return SourcePrimeValue != 0u && PrecodePrimeValue != 0u &&
        SourceCount == source_count && PrecodeCount == precode_count &&
        MixCount == mix_count;
}

std::vector<uint32_t> GeneratePacketMatrixRowWithRuntime(
    uint32_t source_count,
    uint32_t precode_count,
    uint32_t block_id,
    const PacketRowConfig& config,
    const PacketRowRuntime& runtime)
{
    std::vector<uint32_t> row;
    const bool generated = ForEachPacketMatrixColumn(
        source_count,
        precode_count,
        block_id,
        config,
        runtime,
        [&row](size_t count) {
            row.reserve(count);
            return true;
        },
        [&row](uint32_t column) { row.push_back(column); });
    if (!generated) {
        row.clear();
        return row;
    }
    return row;
}

bool GeneratePacketMatrixRowIntoWithRuntime(
    uint32_t source_count,
    uint32_t precode_count,
    uint32_t block_id,
    const PacketRowConfig& config,
    const PacketRowRuntime& runtime,
    uint32_t* columns_out,
    size_t capacity,
    size_t& count_out)
{
    size_t produced_count = 0u;
    const bool generated = ForEachPacketMatrixColumn(
        source_count,
        precode_count,
        block_id,
        config,
        runtime,
        [&](size_t count) {
            if (columns_out == nullptr || capacity < count) {
                return false;
            }
            const uintptr_t output_begin =
                reinterpret_cast<uintptr_t>(columns_out);
            if (count > (UINTPTR_MAX - output_begin) / sizeof(uint32_t)) {
                return false;
            }
            const uintptr_t output_end =
                output_begin + count * sizeof(uint32_t);
            const uintptr_t count_begin =
                reinterpret_cast<uintptr_t>(&count_out);
            if (count_begin > UINTPTR_MAX - sizeof(count_out)) {
                return false;
            }
            const uintptr_t count_end = count_begin + sizeof(count_out);
            return output_begin >= count_end || count_begin >= output_end;
        },
        [&](uint32_t column) {
            columns_out[produced_count++] = column;
        });
    if (!generated) return false;
    count_out = produced_count;
    return true;
}

std::vector<uint32_t> GeneratePacketMatrixRow(
    uint32_t source_count,
    uint32_t precode_count,
    uint32_t block_id,
    const PacketRowConfig& config)
{
    PacketRowRuntime runtime;
    if (!runtime.Initialize(
            source_count, precode_count, config.MixCount))
    {
        return std::vector<uint32_t>();
    }
    return GeneratePacketMatrixRowWithRuntime(
        source_count, precode_count, block_id, config, runtime);
}

uint32_t PacketPeelSeedFromProfile(
    const SeedProfile& profile,
    uint64_t salt)
{
    const uint64_t seed = MatrixSeedFromProfile(profile, 0u, salt);
    return (uint32_t)seed ^ (uint32_t)(seed >> 32);
}

PacketRowConfig PacketConfigForAttempt(
    const PacketRowConfig& base,
    uint32_t attempt)
{
    PacketRowConfig candidate = base;
    candidate.PeelSeed += attempt * UINT32_C(0x9e3779b9);
    return candidate;
}

PrecodeParams PrecodeParamsForAttempt(
    const PrecodeParams& base,
    uint32_t attempt)
{
    PrecodeParams candidate = base;
    candidate.Seed +=
        (uint64_t)attempt * UINT64_C(0x9e3779b97f4a7c15);
    return candidate;
}

static_assert(
    kCertifiedPacketMixCount == 3u &&
        kCertifiedPacketMixCount == wirehair::RowMixIterator::kColumnCount,
    "the fused packet evaluator requires the certified three-mix contract");

static const uint32_t kPacketTailPairMaxBlockBytes = 32u * 1024u;
static const uint32_t kPacketTailPairMinTerms = 6u;
static const uint32_t kPacketSetXorMaxTerms = 16u;

#if defined(_MSC_VER)
#define WH2_PACKET_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define WH2_PACKET_NOINLINE __attribute__((noinline))
#else
#define WH2_PACKET_NOINLINE
#endif

// Packet rows contain distinct source columns, and their precode columns live
// in a disjoint suffix.  Pair the tail directly without the generic
// accumulator's duplicate-source check.  Keeping this out of line preserves
// the compact common evaluator when the size/degree gate does not select it.
static WH2_PACKET_NOINLINE void EvaluatePacketTailPaired(
    const wirehair::PeelRowParameters& params,
    uint32_t source_count,
    uint32_t precode_count,
    const PacketRowConfig& config,
    const PacketRowRuntime& runtime,
    const uint8_t* intermediate_blocks,
    uint32_t block_bytes,
    uint8_t* block_out)
{
    wirehair::PeelRowIterator source(
        params, (uint16_t)source_count, runtime.SourcePrime());
    const uint8_t* const first_source = intermediate_blocks +
        (size_t)source.GetColumn() * block_bytes;
    (void)source.Iterate();
    gf256_addset_mem(
        block_out,
        first_source,
        intermediate_blocks + (size_t)source.GetColumn() * block_bytes,
        (int)block_bytes);

    const uint8_t* pending = nullptr;
    const auto consume = [&](const uint8_t* term) {
        if (pending)
        {
            gf256_add2_mem(block_out, pending, term, (int)block_bytes);
            pending = nullptr;
        }
        else {
            pending = term;
        }
    };
    while (source.Iterate()) {
        consume(intermediate_blocks +
            (size_t)source.GetColumn() * block_bytes);
    }
    if (config.MixCount == 1u)
    {
        consume(intermediate_blocks +
            (size_t)(source_count + params.MixFirst) * block_bytes);
    }
    else
    {
        const wirehair::RowMixIterator mix(
            params, (uint16_t)precode_count, runtime.PrecodePrime());
        for (uint32_t i = 0u; i < config.MixCount; ++i) {
            consume(intermediate_blocks +
                (size_t)(source_count + mix.Columns[i]) * block_bytes);
        }
    }
    if (pending) {
        gf256_add_mem(block_out, pending, (int)block_bytes);
    }
}

// Complete small packet rows fit the fixed-count set-XOR family.  Gather the
// already-distinct source terms and disjoint precode suffix once, then consume
// every input in a single destination pass.
static WH2_PACKET_NOINLINE void EvaluatePacketSetXor(
    const wirehair::PeelRowParameters& params,
    uint32_t source_count,
    uint32_t precode_count,
    const PacketRowConfig& config,
    const PacketRowRuntime& runtime,
    const uint8_t* intermediate_blocks,
    uint32_t block_bytes,
    uint8_t* block_out)
{
    wirehair::PeelRowIterator source(
        params, (uint16_t)source_count, runtime.SourcePrime());
    const void* sources[kPacketSetXorMaxTerms];
    uint32_t count = 0u;
    do {
        sources[count++] = intermediate_blocks +
            (size_t)source.GetColumn() * block_bytes;
    } while (source.Iterate());

    if (config.MixCount == 1u)
    {
        sources[count++] = intermediate_blocks +
            (size_t)(source_count + params.MixFirst) * block_bytes;
    }
    else
    {
        const wirehair::RowMixIterator mix(
            params, (uint16_t)precode_count, runtime.PrecodePrime());
        for (uint32_t i = 0u; i < config.MixCount; ++i) {
            sources[count++] = intermediate_blocks +
                (size_t)(source_count + mix.Columns[i]) * block_bytes;
        }
    }
    gf256_addset_multi_mem(block_out, sources, (int)count, (int)block_bytes);
}

#undef WH2_PACKET_NOINLINE

static bool EvaluatePacketBlockImpl(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const PacketRowRuntime& runtime,
    const uint8_t* intermediate_blocks,
    uint32_t block_bytes,
    uint32_t block_id,
    uint8_t* block_out,
    uint64_t* block_ops_out,
    bool validate_system)
{
    const uint32_t K = system.Params.BlockCount;
    const uint64_t P_wide = (uint64_t)system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    if (!intermediate_blocks ||
        !block_out || block_bytes == 0u || block_bytes > 0x7fffffffu ||
        P_wide > UINT32_MAX ||
        !runtime.IsValidFor(K, (uint32_t)P_wide, config.MixCount))
    {
        return false;
    }
    // The block-size cap above makes this product fit uint64_t even when K
    // and P are both UINT32_MAX.
    const uint64_t intermediate_bytes_wide =
        ((uint64_t)K + P_wide) * block_bytes;
    if (intermediate_bytes_wide >
        (uint64_t)std::numeric_limits<size_t>::max())
    {
        return false;
    }
    if (MemoryRangesOverlap(
            block_out,
            block_bytes,
            intermediate_blocks,
            (size_t)intermediate_bytes_wide))
    {
        return false;
    }
    if (validate_system && !ValidatePrecodeSystem(system)) {
        return false;
    }
    const uint32_t P = (uint32_t)P_wide;

    wirehair::PeelRowParameters params;
    const uint32_t packet_row_seed = PacketRowSeedForBlockId(block_id);
    const uint32_t packet_peel_seed = PacketPeelSeedForBlockId(block_id, config);
    params.Initialize(
        packet_row_seed,
        packet_peel_seed,
        (uint16_t)K,
        (uint16_t)P);
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    /*
        The peel-degree override MUST be applied here as well as in
        InitializePacketRowParameters.

        This is the ENCODER's row builder and that one is the DECODER's, and an
        override applied to only one of them makes the two disagree about the
        degree of every block id whose weight the override moves.  The encoder
        then emits a packet built from one equation while the solver builds a
        different one, the solve returns Wirehair_Success over a system that
        does not describe the data, and the payload comes back corrupt.

        That is not hypothetical -- it was the state of this file, and it made
        every non-stock distribution fail 100% end to end while looking like a
        large win in the precode paths, which build and solve through this same
        overridden routine and so never disagree with themselves.

        THE IDENTITY CONTROL CANNOT CATCH THIS, which is why it survived a full
        campaign.  Feeding the stock PMF maps to the stock degrees on BOTH
        sides by construction, so the two paths agree and the control passes;
        the control is a faithful test of one hook and is structurally blind to
        an asymmetry between two, because its only input is the fixed point of
        both.  Any future hook that alters row construction has to be applied
        at every site that builds a row, and the check for that is a NON-stock
        distribution decoding end to end, not a stock one.

        WHAT THE ASYMMETRY LOOKED LIKE FROM THE OUTSIDE, recorded because it was
        mistaken for physics for most of a day: every candidate distribution
        "failed 100% of trials", always all-or-nothing and never partial, because
        the corruption is deterministic in block id.  That produced a published
        "decodable region" table whose tolerances were really the L1 distance at
        which the first IN-RANGE block id's sampled degree flips -- so p1 x1.5 at
        L1 0.00772 passed (first flip at id 159, out of range) while P(2) -1% at
        the SMALLER L1 0.00503 failed (first flip at id 26).  Distance did not
        predict; flip position did.  With both sides agreeing, distributions as
        far as L1 0.676 from stock decode every trial at zero overhead.

        AND NOTE WHAT THE END-TO-END GATE MUST READ.  After this fix, `compare`'s
        FAIL column does not discriminate peel distributions at all: an
        all-mass-on-degree-1 peel matrix still reports 0 failures, because the
        codec spends OVERHEAD rather than failing.  The discriminating column is
        OH_mean (with OH95) -- all-degree-1 costs 5.02 extra packets where the
        shipped law costs 0.0000.  A gate reading FAIL is vacuous.
    */
    if (!ApplyPeelDegreeOverride(
            packet_row_seed, packet_peel_seed, (uint16_t)K, params))
    {
        return false;
    }
#endif
    uint64_t operations = 1u;
    // The existing schedules are already optimal until the row contains six
    // total terms.  Above that crossover, pairing the complete tail removes
    // at least one destination read/write pass.
    const uint32_t packet_terms =
        (uint32_t)params.PeelCount + config.MixCount;
    // These two width gates deliberately read the REAL width rather than
    // DispatchBlockBytes: this routine re-encodes a packet into block_out and
    // its entry check above rejects a zero width outright, so a payload-free
    // solve never reaches them and there is no counting mode to preserve.
    if (block_bytes <= kPacketTailPairMaxBlockBytes &&
        packet_terms >= kPacketTailPairMinTerms)
    {
        if (K >= 10000u && block_bytes >= 1280u &&
            packet_terms <= kPacketSetXorMaxTerms)
        {
            EvaluatePacketSetXor(
                params, K, P, config, runtime, intermediate_blocks,
                block_bytes, block_out);
        }
        else {
            EvaluatePacketTailPaired(
                params, K, P, config, runtime, intermediate_blocks,
                block_bytes, block_out);
        }
        if (block_ops_out) {
            *block_ops_out = packet_terms;
        }
        return true;
    }
    wirehair::PeelRowIterator source(
        params, (uint16_t)K, runtime.SourcePrime());
    const uint8_t* first_source =
        intermediate_blocks + (size_t)source.GetColumn() * block_bytes;
    if (config.MixCount == kCertifiedPacketMixCount)
    {
        const wirehair::RowMixIterator mix(
            params, (uint16_t)P, runtime.PrecodePrime());
        // The certified three-mix contract mirrors the production codec's
        // fused evaluation schedule: initialize from two sources with addset
        // (or source 0 + mix 0 for a weight-one row), then consume the final
        // two mix terms together with add2.  This removes two full destination
        // read/write passes without changing the equation or work accounting.
        if (source.Iterate())
        {
            gf256_addset_mem(
                block_out,
                first_source,
                intermediate_blocks +
                    (size_t)source.GetColumn() * block_bytes,
                (int)block_bytes);
            ++operations;
            while (source.Iterate())
            {
                gf256_add_mem(
                    block_out,
                    intermediate_blocks +
                        (size_t)source.GetColumn() * block_bytes,
                    (int)block_bytes);
                ++operations;
            }
            gf256_add_mem(
                block_out,
                intermediate_blocks +
                    (size_t)(K + mix.Columns[0]) * block_bytes,
                (int)block_bytes);
            ++operations;
        }
        else
        {
            gf256_addset_mem(
                block_out,
                first_source,
                intermediate_blocks +
                    (size_t)(K + mix.Columns[0]) * block_bytes,
                (int)block_bytes);
            ++operations;
        }
        gf256_add2_mem(
            block_out,
            intermediate_blocks +
                (size_t)(K + mix.Columns[1]) * block_bytes,
            intermediate_blocks +
                (size_t)(K + mix.Columns[2]) * block_bytes,
            (int)block_bytes);
        operations += 2u;
    }
    else if (config.MixCount == 2u)
    {
        const wirehair::RowMixIterator mix(
            params, (uint16_t)P, runtime.PrecodePrime());
        // The two-mix packet contract has the same fused opportunities:
        // initialize from two sources when possible, then consume both mix
        // terms in one destination pass.  A weight-one source row instead
        // initializes from source 0 + mix 0 before adding mix 1.
        if (source.Iterate())
        {
            gf256_addset_mem(
                block_out,
                first_source,
                intermediate_blocks +
                    (size_t)source.GetColumn() * block_bytes,
                (int)block_bytes);
            ++operations;
            while (source.Iterate())
            {
                gf256_add_mem(
                    block_out,
                    intermediate_blocks +
                        (size_t)source.GetColumn() * block_bytes,
                    (int)block_bytes);
                ++operations;
            }
            gf256_add2_mem(
                block_out,
                intermediate_blocks +
                    (size_t)(K + mix.Columns[0]) * block_bytes,
                intermediate_blocks +
                    (size_t)(K + mix.Columns[1]) * block_bytes,
                (int)block_bytes);
            operations += 2u;
        }
        else
        {
            gf256_addset_mem(
                block_out,
                first_source,
                intermediate_blocks +
                    (size_t)(K + mix.Columns[0]) * block_bytes,
                (int)block_bytes);
            ++operations;
            gf256_add_mem(
                block_out,
                intermediate_blocks +
                    (size_t)(K + mix.Columns[1]) * block_bytes,
                (int)block_bytes);
            ++operations;
        }
    }
    else
    {
        // Runtime validation leaves exactly the one-mix experiment here.
        // RowMixIterator's first output is MixFirst, so avoid producing the
        // unused second and third columns.
        const uint8_t* const mix_source = intermediate_blocks +
            (size_t)(K + params.MixFirst) * block_bytes;
        if (params.PeelCount == 1u)
        {
            gf256_addset_mem(
                block_out, first_source, mix_source, (int)block_bytes);
            ++operations;
        }
        else
        {
            // Initialize from the first two source terms.  For three or more
            // sources, leave the final source pending so it can be consumed
            // with the mix term in one destination pass.
            (void)source.Iterate();
            gf256_addset_mem(
                block_out,
                first_source,
                intermediate_blocks +
                    (size_t)source.GetColumn() * block_bytes,
                (int)block_bytes);
            ++operations;
            if (params.PeelCount == 2u)
            {
                gf256_add_mem(
                    block_out, mix_source, (int)block_bytes);
                ++operations;
            }
            else
            {
                for (uint32_t source_index = 2u;
                     source_index + 1u < params.PeelCount;
                     ++source_index)
                {
                    (void)source.Iterate();
                    gf256_add_mem(
                        block_out,
                        intermediate_blocks +
                            (size_t)source.GetColumn() * block_bytes,
                        (int)block_bytes);
                    ++operations;
                }
                (void)source.Iterate();
                gf256_add2_mem(
                    block_out,
                    intermediate_blocks +
                        (size_t)source.GetColumn() * block_bytes,
                    mix_source,
                    (int)block_bytes);
                operations += 2u;
            }
        }
    }
    if (block_ops_out) {
        *block_ops_out = operations;
    }
    return true;
}

bool EvaluatePacketBlock(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const uint8_t* intermediate_blocks,
    uint32_t block_bytes,
    uint32_t block_id,
    uint8_t* block_out,
    uint64_t* block_ops_out)
{
    if (gf256_init() != 0) {
        return false;
    }
    const uint64_t P_wide = (uint64_t)system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    PacketRowRuntime runtime;
    if (P_wide > UINT32_MAX ||
        !runtime.Initialize(
            system.Params.BlockCount,
            (uint32_t)P_wide,
            config.MixCount))
    {
        return false;
    }
    return EvaluatePacketBlockImpl(
        system, config, runtime, intermediate_blocks, block_bytes, block_id,
        block_out, block_ops_out, true);
}

bool EvaluatePacketBlockForValidatedSystem(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const uint8_t* intermediate_blocks,
    uint32_t block_bytes,
    uint32_t block_id,
    uint8_t* block_out,
    uint64_t* block_ops_out)
{
    const uint64_t P_wide = (uint64_t)system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    PacketRowRuntime runtime;
    if (P_wide > UINT32_MAX ||
        !runtime.Initialize(
            system.Params.BlockCount,
            (uint32_t)P_wide,
            config.MixCount))
    {
        return false;
    }
    return EvaluatePacketBlockForValidatedSystemWithRuntime(
        system, config, runtime, intermediate_blocks, block_bytes, block_id,
        block_out, block_ops_out);
}

bool EvaluatePacketBlockForValidatedSystemWithRuntime(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const PacketRowRuntime& runtime,
    const uint8_t* intermediate_blocks,
    uint32_t block_bytes,
    uint32_t block_id,
    uint8_t* block_out,
    uint64_t* block_ops_out)
{
    return EvaluatePacketBlockImpl(
        system, config, runtime, intermediate_blocks, block_bytes, block_id,
        block_out, block_ops_out, false);
}

WirehairResult SolvePrecodeSystem(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const std::vector<SolvePacket>& packets,
    uint32_t block_bytes,
    std::vector<uint8_t>& intermediate_blocks_out,
    PrecodeSolveStats* stats,
    PrecodeSolveResumeState* resume_state)
{
    const uint64_t P_wide = (uint64_t)system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    PacketRowRuntime runtime;
    if (P_wide > UINT32_MAX ||
        !runtime.Initialize(
            system.Params.BlockCount,
            (uint32_t)P_wide,
            config.MixCount))
    {
        return Wirehair_InvalidInput;
    }
    return SolvePrecodeSystemWithRuntime(
        system, config, runtime, packets, block_bytes,
        intermediate_blocks_out, stats, resume_state);
}

static WirehairResult SolvePrecodeSystemImpl(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const PacketRowRuntime& runtime,
    const std::vector<SolvePacket>& packets,
    uint32_t block_bytes,
    std::vector<uint8_t>& intermediate_blocks_out,
    PrecodeSolveStats* stats,
    PrecodeSolveResumeState* resume_state,
    bool validate_system)
{
    PrecodeSolveStats st = {};
    const uint32_t K = system.Params.BlockCount;
    const uint32_t S = system.Params.Staircase;
    const uint32_t D2 = system.Params.DenseRows;
    const uint32_t H = system.Params.HeavyRows;
    const uint64_t P_wide = (uint64_t)S + D2 + H;
    if (P_wide > UINT32_MAX ||
        !runtime.IsValidFor(K, (uint32_t)P_wide, config.MixCount) ||
        (validate_system && !ValidatePrecodeSystem(system)) ||
        (system.Params.Field == CompletionField::MixedGF256GF16 &&
         (block_bytes & 1u) != 0u))
    {
        return Wirehair_InvalidInput;
    }
    const uint32_t P = (uint32_t)P_wide;
    const uint32_t L = K + P;
    // A zero block size is a payload-free solve: the caller wants the
    // sequence of block operations, not the bytes they would move.  The
    // storage arithmetic is then exact rather than rejected -- count blocks
    // of zero bytes really do need zero bytes -- and every block buffer below
    // still holds real backing store, so no block pointer is null.  The block
    // count is still required to be nonzero, exactly as CheckedBlockStorage
    // requires it at every other width.
    const auto checked_block_storage =
        [block_bytes](uint32_t block_count, size_t& bytes_out) -> bool
    {
        if (block_bytes == 0u) {
            bytes_out = 0u;
            return block_count != 0u;
        }
        return CheckedBlockStorage(block_count, block_bytes, bytes_out);
    };
    size_t value_bytes = 0u;
    if (!checked_block_storage(L, value_bytes))
    {
        return Wirehair_InvalidInput;
    }
    for (const SolvePacket& packet : packets) {
        if (!packet.Data) {
            return Wirehair_InvalidInput;
        }
    }
    if (packets.size() < K) {
        return Wirehair_NeedMore;
    }
    if (gf256_init() != 0) {
        return Wirehair_UnsupportedPlatform;
    }
    if (system.Params.Field == CompletionField::MixedGF256GF16 &&
        !InitializeGF16())
    {
        return Wirehair_UnsupportedPlatform;
    }

    const auto terminal_error = [&]() -> WirehairResult {
        if (stats) {
            *stats = st;
        }
        return Wirehair_Error;
    };

    try
    {
        typedef std::chrono::steady_clock SolveClock;
        SolveClock::time_point phase_start = SolveClock::now();
        const uint64_t row_count_wide =
            (uint64_t)S + D2 + packets.size();
        if (row_count_wide > UINT32_MAX) {
            return Wirehair_OOM;
        }
        const size_t row_count = (size_t)row_count_wide;
        size_t reference_count = 0u;
        const auto add_references = [&reference_count](size_t count) {
            if (count > std::numeric_limits<size_t>::max() - reference_count) {
                return false;
            }
            reference_count += count;
            return true;
        };
        for (const std::vector<uint32_t>& columns : system.StaircaseRows) {
            if (!add_references(columns.size())) {
                return Wirehair_OOM;
            }
        }
        for (const std::vector<uint32_t>& columns : system.DenseRowColumns) {
            if (!add_references(columns.size())) {
                return Wirehair_OOM;
            }
        }
        // The peel distribution averages fewer than seven source references.
        // Reserve one modestly padded estimate and generate each packet row
        // only once.  The previous exact reserve initialized every packet PRNG
        // twice; vector growth remains a safe fallback for an unusually heavy
        // finite sample.
        static const size_t kReservedPeelReferencesPerPacket = 7u;
        const size_t reserve_references_per_packet =
            kReservedPeelReferencesPerPacket + config.MixCount;
        if (packets.size() >
            (std::numeric_limits<size_t>::max() - reference_count) /
                reserve_references_per_packet)
        {
            return Wirehair_OOM;
        }
        const size_t reserved_reference_count = reference_count +
            packets.size() * reserve_references_per_packet;

        BinaryEquationArena rows;
        rows.Initialize(row_count, reserved_reference_count);
        for (const std::vector<uint32_t>& columns : system.StaircaseRows)
        {
            rows.AppendRow(columns, nullptr);
            st.BinaryRowReferences += columns.size();
        }
        for (const std::vector<uint32_t>& columns : system.DenseRowColumns)
        {
            rows.AppendRow(columns, nullptr);
            st.BinaryRowReferences += columns.size();
        }
        for (const SolvePacket& packet : packets)
        {
            size_t packet_references = 0u;
            rows.BeginRow(packet.Data);
            const bool generated = ForEachPacketMatrixColumn(
                K,
                P,
                packet.BlockId,
                config,
                runtime,
                [&packet_references](size_t count) {
                    packet_references = count;
                    return true;
                },
                [&rows](uint32_t column) { rows.AppendColumn(column); });
            if (!generated || packet_references == 0u) {
                return Wirehair_InvalidInput;
            }
            if (!add_references(packet_references)) {
                return Wirehair_OOM;
            }
            rows.EndRow();
            st.BinaryRowReferences += packet_references;
        }
        if (!rows.IsComplete(row_count, reference_count)) {
            return Wirehair_Error;
        }
        st.BinaryRowStorageBytes = rows.StorageBytes();
        st.BinaryRowStorageAllocations = rows.StorageAllocations();
        st.PacketRows = (uint32_t)packets.size();
        SolveClock::time_point phase_end = SolveClock::now();
        st.BuildNanoseconds = (uint64_t)
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                phase_end - phase_start).count();

        phase_start = phase_end;
        PeelResult peel = PeelBinaryRows(L, rows);
        st.BinaryAdjacencyStorageBytes = peel.AdjacencyStorageBytes;
        st.BinaryAdjacencyStorageAllocations =
            peel.AdjacencyStorageAllocations;
        st.PeelAdjacencyVisits = peel.AdjacencyVisits;
        st.PeelRowScanSteps = peel.RowScanSteps;
        st.PeelHeapOperations = peel.HeapOperations;
        if (peel.PeelOrder.size() + peel.InactiveOrder.size() != L) {
            return terminal_error();
        }
        const uint32_t R = (uint32_t)peel.InactiveOrder.size();
        st.PeeledColumns = (uint32_t)peel.PeelOrder.size();
        st.InactivatedColumns = R;
        phase_end = SolveClock::now();
        st.PeelNanoseconds = (uint64_t)
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                phase_end - phase_start).count();
        if (R > kMaxInactiveColumns)
        {
            if (stats) {
                *stats = st;
            }
            return Wirehair_NeedMore;
        }
        phase_start = phase_end;
        const uint32_t words = PackedWordCount(R);

        std::vector<uint32_t> inactive_index(L, UINT32_MAX);
        for (uint32_t i = 0; i < R; ++i) {
            inactive_index[peel.InactiveOrder[i]] = i;
        }

        size_t projection_words = 0u;
        if (words != 0u &&
            (uint64_t)L * words >
                (uint64_t)std::numeric_limits<size_t>::max() /
                    sizeof(uint64_t))
        {
            return Wirehair_OOM;
        }
        projection_words = (size_t)L * words;
        std::vector<uint64_t> projection(projection_words, 0u);
        std::vector<uint8_t> values;
        values.reserve(BlockStorageCapacity(value_bytes));
#if defined(__linux__) && defined(MADV_HUGEPAGE)
        // Request transparent huge pages before the first touch.  The value
        // workspace is both large and repeatedly scanned by payload XORs;
        // normal 4-KiB first-touch faults are a material part of cold solve
        // latency.  madvise is advisory, so unsupported/disabled THP retains
        // the standard vector behavior.
        if (value_bytes >= (2u << 20))
        {
            const long system_page_bytes = sysconf(_SC_PAGESIZE);
            const uintptr_t page_bytes = system_page_bytes > 0 ?
                (uintptr_t)system_page_bytes : 0u;
            const uintptr_t begin =
                reinterpret_cast<uintptr_t>(values.data());
            if (page_bytes != 0u &&
                (page_bytes & (page_bytes - 1u)) == 0u &&
                begin <= std::numeric_limits<uintptr_t>::max() -
                    (page_bytes - 1u))
            {
                const uintptr_t aligned =
                    (begin + page_bytes - 1u) & ~(page_bytes - 1u);
                const size_t skipped = (size_t)(aligned - begin);
                const size_t advised_bytes = skipped < value_bytes ?
                    (value_bytes - skipped) &
                        ~(size_t)(page_bytes - 1u) : 0u;
                if (advised_bytes >= (2u << 20)) {
                    (void)madvise(
                        reinterpret_cast<void*>(aligned),
                        advised_bytes,
                        MADV_HUGEPAGE);
#if defined(MADV_POPULATE_WRITE)
                    // Bulk prefaulting has a fixed syscall/page-table cost.
                    // Paired cold-solve measurements show that it pays for
                    // K >= 32000 systems with at least 32 MiB of full pages;
                    // smaller K with unusually wide blocks does not reliably
                    // amortize it.
                    if (K >= 32000u && advised_bytes >= (32u << 20)) {
                        (void)madvise(
                            reinterpret_cast<void*>(aligned),
                            advised_bytes,
                            MADV_POPULATE_WRITE);
                    }
#endif
                }
            }
        }
#endif
        values.resize(value_bytes, 0u);
        st.BlockZeroFills += L;
        std::vector<uint64_t> accumulator(words, 0u);
        // DISPATCH, not sizing: the fused initializer changes how the peel
        // projection consumes its sources.
        const bool enable_fused_block_initialization =
            K >= kFusedBlockXorInitMinBlockCount &&
            DispatchBlockBytes(block_bytes) >=
                kFusedBlockXorInitMinBlockBytes;

        // Affine projection of peeled columns onto inactive variables.  The
        // block stored in values[column] is the constant term.
        for (uint32_t column : peel.PeelOrder)
        {
            std::fill(accumulator.begin(), accumulator.end(), uint64_t{0});
            uint8_t* constant =
                values.data() + (size_t)column * block_bytes;
            const uint32_t solve_row = peel.SolveRow[column];
            const BinaryEquationView equation =
                rows[solve_row];
            if (equation.Columns.size() == 0u) {
                return terminal_error();
            }
            const size_t initialization_sources =
                equation.Columns.size() - 1u +
                (equation.Data ? 1u : 0u);
            uint32_t solve_column_offset = UINT32_MAX;
            if (enable_fused_block_initialization &&
                initialization_sources <= 16u)
            {
                BatchedBlockXorInitializer constant_xor(
                    constant, block_bytes, equation.Data, true);
                // The payload is consumed as the initializer's first
                // source instead of being copied in.  It is the same
                // logical block copy the unfused branch below counts, and
                // it is never Add()ed, so this cannot double count.
                if (equation.Data) {
                    ++st.BlockCopies;
                }
                solve_column_offset = AccumulatePeeledProjectionConstant(
                    column, equation, inactive_index, words, projection,
                    values, block_bytes, accumulator, constant_xor, st);
                constant_xor.Flush();
            }
            else
            {
                if (equation.Data) {
                    std::memcpy(constant, equation.Data, block_bytes);
                    ++st.BlockCopies;
                }
                BatchedBlockXorAccumulator constant_xor(
                    constant, block_bytes);
                solve_column_offset = AccumulatePeeledProjectionConstant(
                    column, equation, inactive_index, words, projection,
                    values, block_bytes, accumulator, constant_xor, st);
                constant_xor.Flush();
            }
            CAT_DEBUG_ASSERT(solve_column_offset != UINT32_MAX);
            if (solve_column_offset == UINT32_MAX) {
                return terminal_error();
            }
            rows.MoveSolveColumnToEnd(solve_row, solve_column_offset);
            for (uint32_t w = 0; w < words; ++w) {
                projection[(size_t)column * words + w] = accumulator[w];
            }
        }
        phase_end = SolveClock::now();
        st.ProjectNanoseconds = (uint64_t)
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                phase_end - phase_start).count();
        phase_start = phase_end;

        if (R == 0u)
        {
            // No residual variables: every unused binary row and every actual
            // heavy equation must reduce to zero.
            std::vector<uint8_t> rhs;
            for (uint32_t r = 0; r < (uint32_t)rows.size(); ++r)
            {
                if (peel.UsedRows[r]) {
                    continue;
                }
                // Allocated on first use.  The test is on capacity, not
                // size, because a zero-width solve leaves the buffer
                // permanently empty and would otherwise keep reallocating a
                // null data() into the XOR kernels below.
                if (rhs.capacity() == 0u) {
                    rhs = MakeBlockStorage(block_bytes);
                    ++st.BlockZeroFills;
                }
                BatchedBlockXorInitializer rhs_xor(
                    rhs.data(), block_bytes, rows[r].Data);
                if (rows[r].Data) {
                    ++st.BlockCopies;
                }
                for (uint32_t column : rows[r].Columns) {
                    rhs_xor.Add(
                        values.data() + (size_t)column * block_bytes);
                    ++st.BlockXors;
                }
                rhs_xor.Flush();
                if (!RowIsZero(rhs.data(), block_bytes)) {
                    return terminal_error();
                }
            }
            // A zero-width solve has no payload to re-derive, so every row
            // it could compare is empty and the check is vacuous.
            // VerifyPrecodeSolution() keeps rejecting a zero width as invalid
            // input -- it exists to compare bytes -- so it is skipped here
            // rather than called with nothing to compare.
            if (block_bytes != 0u &&
                !VerifyPrecodeSolution(
                    system,
                    config,
                    packets,
                    values.data(),
                    block_bytes))
            {
                return terminal_error();
            }
            std::vector<uint8_t> committed;
            committed.swap(values);
            intermediate_blocks_out.swap(committed);
            phase_end = SolveClock::now();
            st.ResidualNanoseconds = (uint64_t)
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    phase_end - phase_start).count();
            if (stats) {
                *stats = st;
            }
            return Wirehair_Success;
        }

        const bool mixed_completion =
            system.Params.Field == CompletionField::MixedGF256GF16;
        // The mixed path's residual rows are binary until its small extension
        // quotient.  Other completion fields retain their GF(256) byte rows.
        if ((!mixed_completion &&
                (uint64_t)R * R >
                    (uint64_t)std::numeric_limits<size_t>::max()) ||
            (mixed_completion &&
                (uint64_t)R * words >
                    (uint64_t)std::numeric_limits<size_t>::max() /
                        sizeof(uint64_t)))
        {
            return Wirehair_OOM;
        }
        std::vector<uint8_t> pivot_coeff(
            mixed_completion ? 0u : (size_t)R * R, 0u);
        std::vector<uint64_t> binary_pivot_coeff(
            mixed_completion ? (size_t)R * words : 0u, uint64_t{0});
        size_t residual_value_bytes = 0u;
        if (!checked_block_storage(R, residual_value_bytes)) {
            return Wirehair_OOM;
        }
        std::vector<uint8_t> pivot_rhs =
            MakeBlockStorage(residual_value_bytes);
        st.BlockZeroFills += R + 1u;
        std::vector<uint8_t> have_pivot(R, 0u);
        std::vector<uint8_t> coeff(mixed_completion ? 0u : R, 0u);
        std::vector<uint8_t> rhs = MakeBlockStorage(block_bytes);
        uint32_t rank = 0u;
        const uint32_t batched_rhs_min_block_bytes =
            system.Params.Field == CompletionField::MixedGF256GF16 ?
                kBatchedResidualRhsMinBlockBytes :
                kNeverBatchResidualRhs;

        // Project every unused binary row.
        for (uint32_t r = 0; r < (uint32_t)rows.size(); ++r)
        {
            if (peel.UsedRows[r]) {
                continue;
            }
            ++st.ResidualRows;
            std::fill(
                accumulator.begin(), accumulator.end(), uint64_t{0});
            BatchedBlockXorInitializer rhs_xor(
                rhs.data(), block_bytes, rows[r].Data);
            BatchedProjectionXorAccumulator projection_xor(
                accumulator.data(), projection.data(), words);
            for (uint32_t column : rows[r].Columns)
            {
                const uint32_t index = inactive_index[column];
                if (index != UINT32_MAX) {
                    accumulator[index >> 6] ^=
                        UINT64_C(1) << (index & 63u);
                }
                else
                {
                    projection_xor.Add(column);
                    rhs_xor.Add(
                        values.data() + (size_t)column * block_bytes);
                    ++st.BlockXors;
                    st.ProjectionWordXors += words;
                }
            }
            projection_xor.Flush();
            rhs_xor.Flush();
            ResidualInsertResult insertion;
            if (mixed_completion)
            {
                insertion = InsertPackedBinaryResidualRow(
                    accumulator, rhs, R, words, block_bytes,
                    binary_pivot_coeff, pivot_rhs, have_pivot,
                    rank, st, batched_rhs_min_block_bytes);
            }
            else
            {
                // Generic GF(256) completion consumes one byte per projected
                // coefficient.  Expand each final parity bit exactly once.
                std::fill(coeff.begin(), coeff.end(), uint8_t{0});
                for (uint32_t w = 0; w < words; ++w)
                {
                    uint64_t word = accumulator[w];
                    while (word != 0u)
                    {
                        const uint32_t bit =
                            wirehair::NonzeroLowestBitIndex64(word);
                        const uint32_t projected = (w << 6) + bit;
                        if (projected < R) {
                            coeff[projected] = 1u;
                        }
                        word &= word - 1u;
                    }
                }
                insertion = InsertResidualRow(
                    coeff, rhs, R, block_bytes,
                    pivot_coeff, pivot_rhs, have_pivot,
                    rank, st, true, batched_rhs_min_block_bytes);
            }
            if (insertion == ResidualInsertResult::Inconsistent)
            {
                return terminal_error();
            }
        }

        if (mixed_completion)
        {
            const WirehairResult mixed_result =
                SolveMixedCompletionQuotient(
                    system, L, R, words, block_bytes,
                    inactive_index, peel.InactiveOrder, projection,
                    binary_pivot_coeff, pivot_rhs, have_pivot, rank,
                    values, st);
            phase_end = SolveClock::now();
            st.ResidualNanoseconds = (uint64_t)
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    phase_end - phase_start).count();
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
            if (mixed_result == Wirehair_NeedMore &&
                MixedNullWitnessSink && R >= rank)
            {
                const uint32_t quotient_columns = R - rank;
                if (quotient_columns != 0u &&
                    quotient_columns <=
                        kMaxMixedNullWitnessQuotientColumns &&
                    quotient_columns <= system.Params.HeavyRows &&
                    st.ResidualRank >= rank && st.ResidualRank < R)
                {
                    CaptureMixedNullWitness(
                        system, L, R, words, rows, inactive_index,
                        peel.InactiveOrder, projection,
                        binary_pivot_coeff, have_pivot, rank,
                        st.ResidualRank - rank, MixedNullWitnessSink);
                }
            }
#endif
            if (mixed_result != Wirehair_Success)
            {
                if (stats) *stats = st;
                return mixed_result;
            }
            phase_start = phase_end;
        }
        else
        {
        // Heavy RHS values are bucketed by the coefficient period, avoiding
        // H*L full-block multiplications.  Heavy coefficient vectors are
        // packed eight rows per word so each projection bit is visited once.
        st.BinaryResidualRank = rank;
        const uint32_t window = 256u - H;
        const uint32_t heavy_words = (H + 7u) / 8u;
        // Building the process-local table on first use costs one complete
        // coefficient period.  Tiny systems retain the on-demand path so
        // their cold first solve cannot pay more coefficient work merely to
        // populate entries they will not visit.
        const bool cached_periodic =
            system.Params.HeavyFamily ==
                HeavyCoefficientFamily::PeriodicCauchy &&
            H == kCachedPeriodicHeavyRows &&
            L >= kCachedPeriodicWindow;
        if (heavy_words != 0u &&
            (uint64_t)R * heavy_words >
                (uint64_t)std::numeric_limits<size_t>::max() /
                    sizeof(uint64_t))
        {
            return Wirehair_OOM;
        }
        std::vector<uint64_t> projected_heavy(
            (size_t)R * heavy_words, 0u);
        std::vector<uint64_t> packed_heavy(
            cached_periodic ? 0u : heavy_words, uint64_t{0});
        const uint64_t* periodic_packed = cached_periodic ?
            CachedPeriodicHeavyTable().data() :
            nullptr;
        if (H > 0u)
        {
            uint32_t residue = 0u;
            for (uint32_t column = 0; column < L; ++column)
            {
                const uint64_t* column_heavy = nullptr;
                if (cached_periodic)
                {
                    column_heavy = periodic_packed +
                        (size_t)residue * heavy_words;
                    if (++residue == window) {
                        residue = 0u;
                    }
                }
                else
                {
                    std::fill(
                        packed_heavy.begin(), packed_heavy.end(),
                        uint64_t{0});
                    for (uint32_t heavy = 0; heavy < H; ++heavy) {
                        packed_heavy[heavy >> 3] |=
                            (uint64_t)HeavyCoefficientForParams(
                                system.Params, heavy, column) <<
                            ((heavy & 7u) * 8u);
                    }
                    column_heavy = packed_heavy.data();
                }
                const auto xor_packed = [&](uint32_t index) {
                    uint64_t* destination = projected_heavy.data() +
                        (size_t)index * heavy_words;
                    for (uint32_t w = 0; w < heavy_words; ++w) {
                        destination[w] ^= column_heavy[w];
                    }
                };
                const uint32_t inactive = inactive_index[column];
                if (inactive != UINT32_MAX) {
                    xor_packed(inactive);
                    continue;
                }
                const uint64_t* bits =
                    projection.data() + (size_t)column * words;
                for (uint32_t w = 0; w < words; ++w)
                {
                    uint64_t word = bits[w];
                    while (word != 0u)
                    {
                        const uint32_t bit =
                            wirehair::NonzeroLowestBitIndex64(word);
                        const uint32_t projected = (w << 6) + bit;
                        if (projected < R) {
                            xor_packed(projected);
                        }
                        word &= word - 1u;
                    }
                }
            }
        }

        std::vector<uint8_t> heavy_rhs =
            MakeBlockStorage((size_t)H * block_bytes);
        st.BlockZeroFills += H;
        void* heavy_destinations[128];
        uint8_t heavy_scales[128];
        for (uint32_t heavy = 0; heavy < H; ++heavy) {
            heavy_destinations[heavy] =
                heavy_rhs.data() + (size_t)heavy * block_bytes;
        }
        if (system.Params.HeavyFamily ==
            HeavyCoefficientFamily::PeriodicCauchy)
        {
            std::vector<uint8_t> residue_bucket =
                MakeBlockStorage(block_bytes);
            ++st.BlockZeroFills;
            // Residues at or above L cannot contain a column when L is smaller
            // than the coefficient period.  Processing those empty buckets
            // would issue H full-block muladds from an all-zero source, which
            // is especially expensive for tiny messages and cannot affect the
            // RHS.
            const uint32_t populated_residues = std::min(window, L);
            for (uint32_t residue = 0;
                 residue < populated_residues;
                 ++residue)
            {
                std::fill(
                    residue_bucket.begin(), residue_bucket.end(), uint8_t{0});
                ++st.BlockZeroFills;
                BatchedBlockXorAccumulator bucket_xor(
                    residue_bucket.data(), block_bytes);
                for (uint32_t column = residue; column < L; column += window)
                {
                    if (inactive_index[column] != UINT32_MAX) {
                        continue;
                    }
                    bucket_xor.Add(
                        values.data() + (size_t)column * block_bytes);
                    ++st.BlockXors;
                }
                bucket_xor.Flush();
                if (cached_periodic)
                {
                    const uint64_t* packed = periodic_packed +
                        (size_t)residue * heavy_words;
                    for (uint32_t heavy = 0; heavy < H; ++heavy) {
                        heavy_scales[heavy] = (uint8_t)(
                            packed[heavy >> 3] >>
                            ((heavy & 7u) * 8u));
                    }
                }
                else
                {
                    for (uint32_t heavy = 0; heavy < H; ++heavy) {
                        heavy_scales[heavy] = HeavyCoefficientForParams(
                            system.Params, heavy, residue);
                    }
                }
                AddScaledBlocks(
                    heavy_destinations,
                    heavy_scales,
                    H,
                    residue_bucket.data(),
                    block_bytes,
                    st);
            }
        }
        else
        {
            // Experiment families may depend on the complete column id and
            // therefore cannot use the periodic residue-bucket optimization.
            for (uint32_t column = 0; column < L; ++column)
            {
                if (inactive_index[column] != UINT32_MAX) {
                    continue;
                }
                for (uint32_t heavy = 0; heavy < H; ++heavy) {
                    heavy_scales[heavy] = HeavyCoefficientForParams(
                        system.Params, heavy, column);
                }
                AddScaledBlocks(
                    heavy_destinations,
                    heavy_scales,
                    H,
                    values.data() + (size_t)column * block_bytes,
                    block_bytes,
                    st);
            }
        }
        // The rows inserted so far are binary.  Keep that GF(2) factorization
        // intact and solve the heavy equations only on its free-variable
        // quotient.  Updating all binary pivots with each heavy pivot would
        // turn cheap XOR relationships into an R-wide GF(256) Gauss-Jordan
        // solve.  The quotient has at most H columns in the useful regime.
        const uint32_t binary_rank = rank;
        // The quotient replaces GF(256) block operations with a few more
        // scalar coefficient passes.  Keep the original insertion strategy
        // for MTU-sized blocks where measurements show that trade is neutral
        // or slightly negative; large blocks receive the material win.
        // DISPATCH, not sizing: the quotient trades block operations for
        // scalar coefficient passes, so it changes the operation mix the
        // counters report and a payload-free solve must select it exactly as
        // a large payload does.
        const bool use_binary_quotient =
            DispatchBlockBytes(block_bytes) >= kBinaryQuotientMinBlockBytes;
        std::vector<uint32_t> free_columns;
        if (use_binary_quotient)
        {
            free_columns.reserve(R - binary_rank);
            for (uint32_t column = 0; column < R; ++column) {
                if (!have_pivot[column]) {
                    free_columns.push_back(column);
                }
            }
        }
        const uint32_t quotient_columns =
            (uint32_t)free_columns.size();
        if (use_binary_quotient &&
            binary_rank + quotient_columns != R)
        {
            return terminal_error();
        }

        std::vector<uint8_t> quotient_pivot_coeff(
            (size_t)quotient_columns * quotient_columns, 0u);
        size_t quotient_value_bytes = 0u;
        if (quotient_columns != 0u &&
            !checked_block_storage(quotient_columns, quotient_value_bytes))
        {
            return Wirehair_OOM;
        }
        std::vector<uint8_t> quotient_pivot_rhs =
            MakeBlockStorage(quotient_value_bytes);
        st.BlockZeroFills += quotient_columns +
            (use_binary_quotient ? 1u : 0u);
        std::vector<uint8_t> quotient_have_pivot(
            quotient_columns, 0u);
        std::vector<uint8_t> quotient_coeff(quotient_columns, 0u);
        std::vector<uint8_t> quotient_rhs =
            MakeBlockStorage(use_binary_quotient ? block_bytes : 0u);
        uint32_t quotient_rank = 0u;

        if (use_binary_quotient)
        {
            for (uint32_t heavy = 0; heavy < H; ++heavy)
            {
                ++st.ResidualRows;
                std::fill(coeff.begin(), coeff.end(), uint8_t{0});
                std::memcpy(
                    rhs.data(),
                    heavy_rhs.data() + (size_t)heavy * block_bytes,
                    block_bytes);
                ++st.BlockCopies;
                for (uint32_t index = 0; index < R; ++index) {
                    coeff[index] = (uint8_t)(
                        projected_heavy[
                            (size_t)index * heavy_words + (heavy >> 3)] >>
                        ((heavy & 7u) * 8u));
                }

                // Reduce by the binary basis without inserting into it.  This
                // leaves coefficients and RHS for the free-variable system.
                const ResidualInsertResult reduced = InsertResidualRow(
                    coeff, rhs, R, block_bytes,
                    pivot_coeff, pivot_rhs, have_pivot,
                    rank, st, false,
                    batched_rhs_min_block_bytes);
                if (reduced == ResidualInsertResult::Inconsistent) {
                    return terminal_error();
                }
                for (uint32_t i = 0; i < quotient_columns; ++i) {
                    quotient_coeff[i] = coeff[free_columns[i]];
                }
                std::memcpy(
                    quotient_rhs.data(), rhs.data(), block_bytes);
                ++st.BlockCopies;
                if (InsertResidualRow(
                        quotient_coeff,
                        quotient_rhs,
                        quotient_columns,
                        block_bytes,
                        quotient_pivot_coeff,
                        quotient_pivot_rhs,
                        quotient_have_pivot,
                        quotient_rank,
                        st,
                        true,
                        batched_rhs_min_block_bytes) ==
                    ResidualInsertResult::Inconsistent)
                {
                    return terminal_error();
                }
            }
            rank = binary_rank + quotient_rank;
        }
        else
        {
            for (uint32_t heavy = 0; heavy < H; ++heavy)
            {
                ++st.ResidualRows;
                std::fill(coeff.begin(), coeff.end(), uint8_t{0});
                std::memcpy(
                    rhs.data(),
                    heavy_rhs.data() + (size_t)heavy * block_bytes,
                    block_bytes);
                ++st.BlockCopies;
                for (uint32_t index = 0; index < R; ++index) {
                    coeff[index] = (uint8_t)(
                        projected_heavy[
                            (size_t)index * heavy_words + (heavy >> 3)] >>
                        ((heavy & 7u) * 8u));
                }
                if (InsertResidualRow(
                        coeff, rhs, R, block_bytes,
                        pivot_coeff, pivot_rhs, have_pivot,
                        rank, st, true,
                        batched_rhs_min_block_bytes) ==
                    ResidualInsertResult::Inconsistent)
                {
                    return terminal_error();
                }
            }
        }
        st.ResidualRank = rank;

        // ResumePrecodeSystem accepts arbitrary new binary packet rows in the
        // original R-column coordinates.  A rare deficient cold solve must
        // therefore materialize the legacy combined pivot form before it can
        // publish a checkpoint.  Successful solves stay on the fast quotient
        // path, and callers that do not request resume avoid this fallback.
        if (rank < R && resume_state && use_binary_quotient)
        {
            rank = binary_rank;
            for (uint32_t heavy = 0; heavy < H; ++heavy)
            {
                std::fill(coeff.begin(), coeff.end(), uint8_t{0});
                std::memcpy(
                    rhs.data(),
                    heavy_rhs.data() + (size_t)heavy * block_bytes,
                    block_bytes);
                ++st.BlockCopies;
                for (uint32_t index = 0; index < R; ++index) {
                    coeff[index] = (uint8_t)(
                        projected_heavy[
                            (size_t)index * heavy_words + (heavy >> 3)] >>
                        ((heavy & 7u) * 8u));
                }
                if (InsertResidualRow(
                        coeff, rhs, R, block_bytes,
                        pivot_coeff, pivot_rhs, have_pivot,
                        rank, st, true,
                        batched_rhs_min_block_bytes) ==
                    ResidualInsertResult::Inconsistent)
                {
                    return terminal_error();
                }
            }
            st.ResidualRank = rank;
        }

        phase_end = SolveClock::now();
        st.ResidualNanoseconds = (uint64_t)
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                phase_end - phase_start).count();
        if (rank < R) {
            if (resume_state)
            {
                PrecodeSolveResumeState checkpoint;
                checkpoint.SourceCount = K;
                checkpoint.PrecodeCount = P;
                checkpoint.ColumnCount = L;
                checkpoint.BlockBytes = block_bytes;
                checkpoint.InactiveCount = R;
                checkpoint.ProjectionWords = words;
                checkpoint.Rank = rank;
                checkpoint.Config = config;
                checkpoint.Runtime = runtime;
                checkpoint.Stats = st;
                checkpoint.InactiveIndex.swap(inactive_index);
                checkpoint.InactiveColumns.swap(peel.InactiveOrder);
                checkpoint.Projection.swap(projection);
                checkpoint.Values.swap(values);
                checkpoint.PivotCoefficients.swap(pivot_coeff);
                checkpoint.PivotRhs.swap(pivot_rhs);
                checkpoint.HavePivot.swap(have_pivot);
                checkpoint.CoefficientScratch.swap(coeff);
                checkpoint.RhsScratch.swap(rhs);
                checkpoint.Active = true;
                resume_state->Swap(checkpoint);
            }
            if (stats) {
                *stats = st;
            }
            return Wirehair_NeedMore;
        }

        // Full quotient rank gives each free variable directly.  Reconstruct
        // the binary pivots with XORs from their preserved GF(2) relations.
        if (!use_binary_quotient)
        {
            for (uint32_t i = 0; i < R; ++i)
            {
                if (!have_pivot[i]) {
                    return Wirehair_NeedMore;
                }
                std::memcpy(
                    values.data() +
                        (size_t)peel.InactiveOrder[i] * block_bytes,
                    pivot_rhs.data() + (size_t)i * block_bytes,
                    block_bytes);
                ++st.BlockCopies;
            }
        }
        else
        {
            for (uint32_t i = 0; i < quotient_columns; ++i)
            {
                if (!quotient_have_pivot[i]) {
                    return Wirehair_NeedMore;
                }
                std::memcpy(
                    values.data() + (size_t)peel.InactiveOrder[
                        free_columns[i]] * block_bytes,
                    quotient_pivot_rhs.data() + (size_t)i * block_bytes,
                    block_bytes);
                ++st.BlockCopies;
            }
            for (uint32_t pivot = 0; pivot < R; ++pivot)
            {
                if (!have_pivot[pivot]) {
                    continue;
                }
                uint8_t* value = values.data() +
                    (size_t)peel.InactiveOrder[pivot] * block_bytes;
                std::memcpy(
                    value,
                    pivot_rhs.data() + (size_t)pivot * block_bytes,
                    block_bytes);
                ++st.BlockCopies;
                const uint8_t* relation =
                    pivot_coeff.data() + (size_t)pivot * R;
                for (uint32_t i = 0; i < quotient_columns; ++i)
                {
                    const uint8_t scale = relation[free_columns[i]];
                    AddScaledBlock(
                        value,
                        scale,
                        values.data() + (size_t)peel.InactiveOrder[
                            free_columns[i]] * block_bytes,
                        block_bytes,
                        st);
                }
            }
        }
        phase_start = phase_end;
        }

        // Dependencies of a peeled column were resolved earlier, so forward
        // chronological substitution reconstructs every remaining value.
        for (uint32_t column : peel.PeelOrder)
        {
            uint8_t* value =
                values.data() + (size_t)column * block_bytes;
            const BinaryEquationView equation =
                rows.SolveDependencies(peel.SolveRow[column]);
            const uint32_t sparse_xors =
                (uint32_t)equation.Columns.size();

            // Projection left the affine constant in this value slot.  Once
            // the inactive variables are solved, that relation and the
            // original sparse equation are equivalent reconstructions.  A
            // dense affine relation is usually worse, so count only until it
            // cannot beat the sparse row and use it solely when it removes
            // full payload-block XORs.  Tiny rank proxies avoid this scalar
            // selection work entirely.
            uint32_t projected_xors = 0u;
            // DISPATCH, not sizing: this gate is the one that moves BlockXors
            // between the sparse row and the affine relation, so a
            // payload-free solve has to evaluate it in the production regime.
            if (DispatchBlockBytes(block_bytes) >=
                    kProjectedBackSubMinBlockBytes &&
                words != 0u && sparse_xors != 0u)
            {
                const uint64_t* relation = projection.data() +
                    (size_t)column * words;
                for (uint32_t word_i = 0;
                     word_i < words && projected_xors < sparse_xors;
                     ++word_i)
                {
                    uint64_t word = relation[word_i];
                    while (word != 0u && projected_xors < sparse_xors) {
                        ++projected_xors;
                        word &= word - 1u;
                    }
                }
                if (projected_xors < sparse_xors)
                {
                    BatchedBlockXorAccumulator value_xor(
                        value, block_bytes);
                    for (uint32_t word_i = 0; word_i < words; ++word_i)
                    {
                        uint64_t word = relation[word_i];
                        while (word != 0u)
                        {
                            const uint32_t bit =
                                wirehair::NonzeroLowestBitIndex64(word);
                            const uint32_t index = (word_i << 6) + bit;
                            if (index < R) {
                                value_xor.Add(
                                    values.data() +
                                    (size_t)peel.InactiveOrder[index] *
                                        block_bytes);
                                ++st.BlockXors;
                            }
                            word &= word - 1u;
                        }
                    }
                    value_xor.Flush();
                    continue;
                }
            }
            const size_t initialization_sources =
                equation.Columns.size() +
                (equation.Data ? 1u : 0u);
            if (enable_fused_block_initialization &&
                initialization_sources <= 16u)
            {
                BatchedBlockXorInitializer value_xor(
                    value, block_bytes, equation.Data);
                // Same logical initialization the unfused branch below
                // counts: the payload is a block copy, and its absence a
                // block zero fill.  Neither is ever Add()ed.
                if (equation.Data) {
                    ++st.BlockCopies;
                }
                else {
                    ++st.BlockZeroFills;
                }
                for (uint32_t other : equation.Columns)
                {
                    value_xor.Add(
                        values.data() + (size_t)other * block_bytes);
                    ++st.BlockXors;
                }
                value_xor.Flush();
            }
            else
            {
                if (equation.Data) {
                    std::memcpy(value, equation.Data, block_bytes);
                    ++st.BlockCopies;
                }
                else {
                    std::memset(value, 0, block_bytes);
                    ++st.BlockZeroFills;
                }
                BatchedBlockXorAccumulator value_xor(value, block_bytes);
                for (uint32_t other : equation.Columns)
                {
                    value_xor.Add(
                        values.data() + (size_t)other * block_bytes);
                    ++st.BlockXors;
                }
                value_xor.Flush();
            }
        }
        phase_end = SolveClock::now();
        st.BackSubNanoseconds = (uint64_t)
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                phase_end - phase_start).count();

        intermediate_blocks_out.swap(values);
        if (stats) {
            *stats = st;
        }
        return Wirehair_Success;
    }
    catch (const std::bad_alloc&) {
        return Wirehair_OOM;
    }
    catch (const std::length_error&) {
        return Wirehair_OOM;
    }
}

static bool ShouldUseWideBlockXor(
    const PrecodeSystem& system,
    uint32_t block_bytes)
{
#if !defined(GF256_TARGET_MOBILE) && defined(GF256_TRY_WIDE_XOR)
    // WH2's 3-16 source batches amortize the wide AVX2/AVX-512 kernels at
    // this width.  Smaller payloads and non-mixed solvers retain the compact
    // kernel.
    static const uint32_t kMaxWideBlockCount = 64000u;
    static const uint32_t kMinWideBlockBytes = 512u;
    // DISPATCH, not sizing: this selects the wide GF kernel family for the
    // whole solve.  A payload-free solve selects it too -- the kernels then
    // run on zero-byte blocks -- so its control flow matches production.
    if (system.Params.Field != CompletionField::MixedGF256GF16 ||
        system.Params.BlockCount > kMaxWideBlockCount ||
        DispatchBlockBytes(block_bytes) < kMinWideBlockBytes ||
        gf256_init() != 0)
    {
        return false;
    }
    gf256_x86_cpu_features features = {};
    gf256_get_active_x86_cpu_features(&features);
    return features.AVX2 != 0;
#else
    (void)system;
    (void)block_bytes;
    return false;
#endif
}

class ScopedThreadWideXor
{
public:
    ScopedThreadWideXor()
        : Previous(gf256_set_thread_wide_xor(1))
#if defined(GF256_TRY_TARGET_AVX2) && !defined(__AVX2__)
        , PreviousProjection(TargetWideProjectionXor)
#endif
    {
#if defined(GF256_TRY_TARGET_AVX2) && !defined(__AVX2__)
        TargetWideProjectionXor = true;
#endif
    }

    ~ScopedThreadWideXor()
    {
#if defined(GF256_TRY_TARGET_AVX2) && !defined(__AVX2__)
        TargetWideProjectionXor = PreviousProjection;
#endif
        (void)gf256_set_thread_wide_xor(Previous);
    }

    ScopedThreadWideXor(const ScopedThreadWideXor&) = delete;
    ScopedThreadWideXor& operator=(const ScopedThreadWideXor&) = delete;

private:
    int Previous;
#if defined(GF256_TRY_TARGET_AVX2) && !defined(__AVX2__)
    bool PreviousProjection;
#endif
};

static WirehairResult DispatchPrecodeSolve(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const PacketRowRuntime& runtime,
    const std::vector<SolvePacket>& packets,
    uint32_t block_bytes,
    std::vector<uint8_t>& intermediate_blocks_out,
    PrecodeSolveStats* stats,
    PrecodeSolveResumeState* resume_state,
    bool validate_system)
{
    if (ShouldUseWideBlockXor(system, block_bytes))
    {
        ScopedThreadWideXor wide_xor;
        return SolvePrecodeSystemImpl(
            system, config, runtime, packets, block_bytes,
            intermediate_blocks_out, stats, resume_state, validate_system);
    }
    return SolvePrecodeSystemImpl(
        system, config, runtime, packets, block_bytes,
        intermediate_blocks_out, stats, resume_state, validate_system);
}

namespace {

// Mixed completion starts from binary packet equations.  Keep that residual
// factorization packed in GF(2) instead of expanding each bit into a GF(256)
// byte: all pivots and elimination scales remain exactly zero or one until
// the small GF(2^16) quotient is formed.  This definition intentionally lives
// after the hot solver so its size cannot relocate unrelated hot functions.
#if defined(_MSC_VER)
#define WH2_RESIDUAL_NOINLINE __declspec(noinline)
#elif defined(__ELF__) && (defined(__GNUC__) || defined(__clang__)) && \
    !defined(WIREHAIR_V2_DISABLE_PACKED_RESIDUAL_TEXT_SECTION)
#define WH2_RESIDUAL_NOINLINE \
    __attribute__((noinline, section(".text.wh2_packed_residual")))
#elif defined(__GNUC__) || defined(__clang__)
#define WH2_RESIDUAL_NOINLINE __attribute__((noinline))
#else
#define WH2_RESIDUAL_NOINLINE
#endif

static WH2_RESIDUAL_NOINLINE ResidualInsertResult
InsertPackedBinaryResidualRow(
    std::vector<uint64_t>& coeff,
    std::vector<uint8_t>& rhs,
    uint32_t R,
    uint32_t words,
    uint32_t block_bytes,
    std::vector<uint64_t>& pivot_coeff,
    std::vector<uint8_t>& pivot_rhs,
    std::vector<uint8_t>& have_pivot,
    uint32_t& rank,
    PrecodeSolveStats& stats,
    uint32_t batched_rhs_min_block_bytes)
{
    CAT_DEBUG_ASSERT(
        R > 0u && words == PackedWordCount(R) &&
        coeff.size() == words &&
        pivot_coeff.size() == (size_t)R * words &&
        pivot_rhs.size() == (size_t)R * block_bytes &&
        have_pivot.size() == R);

    if ((R & 63u) != 0u) {
        coeff[words - 1u] &=
            (UINT64_C(1) << (R & 63u)) - UINT64_C(1);
    }

    // DISPATCH, not sizing: the branch picks a reduction algorithm.  Every
    // address and length below still uses the real block_bytes.
    if (DispatchBlockBytes(block_bytes) < batched_rhs_min_block_bytes)
    {
        for (uint32_t column = 0; column < R; ++column)
        {
            const uint64_t bit = UINT64_C(1) << (column & 63u);
            if ((coeff[column >> 6] & bit) == 0u ||
                !have_pivot[column])
            {
                continue;
            }
            const uint64_t* pivot =
                pivot_coeff.data() + (size_t)column * words;
            for (uint32_t word = column >> 6; word < words; ++word) {
                coeff[word] ^= pivot[word];
            }
            stats.ResidualCoeffWordXors += words - (column >> 6);
            AddScaledBlock(
                rhs.data(), 1u,
                pivot_rhs.data() + (size_t)column * block_bytes,
                block_bytes, stats);
        }
    }
    else
    {
        BatchedBlockXorAccumulator rhs_xor(rhs.data(), block_bytes);
        for (uint32_t column = 0; column < R; ++column)
        {
            const uint64_t bit = UINT64_C(1) << (column & 63u);
            if ((coeff[column >> 6] & bit) == 0u ||
                !have_pivot[column])
            {
                continue;
            }
            const uint64_t* pivot =
                pivot_coeff.data() + (size_t)column * words;
            for (uint32_t word = column >> 6; word < words; ++word) {
                coeff[word] ^= pivot[word];
            }
            stats.ResidualCoeffWordXors += words - (column >> 6);
            rhs_xor.Add(
                pivot_rhs.data() + (size_t)column * block_bytes);
            ++stats.BlockXors;
        }
        rhs_xor.Flush();
    }

    uint32_t pivot_column = R;
    for (uint32_t word = 0; word < words; ++word)
    {
        if (coeff[word] != 0u) {
            pivot_column = (word << 6) +
                wirehair::NonzeroLowestBitIndex64(coeff[word]);
            break;
        }
    }
    if (pivot_column >= R) {
        return RowIsZero(rhs.data(), block_bytes) ?
            ResidualInsertResult::Dependent :
            ResidualInsertResult::Inconsistent;
    }
    if (have_pivot[pivot_column]) {
        return ResidualInsertResult::Inconsistent;
    }

    const uint32_t pivot_word = pivot_column >> 6;
    const uint64_t pivot_bit =
        UINT64_C(1) << (pivot_column & 63u);
    // Small payloads retain the compact direct loop.  At wider payloads,
    // batch RHS destinations so the newly inserted source row is loaded once
    // per group while the coefficient elimination remains in exact row order.
    // DISPATCH, not sizing: a payload-free solve takes the batched route the
    // shipping decoder takes, on rows that are zero bytes wide.
    if (DispatchBlockBytes(block_bytes) < 512u ||
        !CanXorBlockIntoDestinationsVectorized())
    {
        for (uint32_t existing = 0; existing < R; ++existing)
        {
            if (!have_pivot[existing]) {
                continue;
            }
            uint64_t* existing_coeff =
                pivot_coeff.data() + (size_t)existing * words;
            if ((existing_coeff[pivot_word] & pivot_bit) == 0u) {
                continue;
            }
            for (uint32_t word = pivot_word; word < words; ++word) {
                existing_coeff[word] ^= coeff[word];
            }
            stats.ResidualCoeffWordXors += words - pivot_word;
            AddScaledBlock(
                pivot_rhs.data() + (size_t)existing * block_bytes,
                1u, rhs.data(), block_bytes, stats);
        }
    }
    else
    {
        // Four destinations plus the shared source stay within the practical
        // stream budget on current cores.  Two destinations waste source
        // loads, while eight or sixteen interleaved wide rows create enough
        // memory stalls to erase the saved instructions.
        static const uint32_t kDestinationBatch = 4u;
        uint8_t* destinations[kDestinationBatch];
        uint32_t destination_count = 0u;
        const auto flush_destinations = [&]() {
            XorBlockIntoDestinations(
                destinations, destination_count,
                rhs.data(), block_bytes);
            destination_count = 0u;
        };
        for (uint32_t existing = 0; existing < R; ++existing)
        {
            if (!have_pivot[existing]) {
                continue;
            }
            uint64_t* existing_coeff =
                pivot_coeff.data() + (size_t)existing * words;
            if ((existing_coeff[pivot_word] & pivot_bit) == 0u) {
                continue;
            }
            for (uint32_t word = pivot_word; word < words; ++word) {
                existing_coeff[word] ^= coeff[word];
            }
            stats.ResidualCoeffWordXors += words - pivot_word;
            destinations[destination_count++] =
                pivot_rhs.data() + (size_t)existing * block_bytes;
            ++stats.BlockXors;
            if (destination_count == kDestinationBatch) {
                flush_destinations();
            }
        }
        flush_destinations();
    }

    std::memcpy(
        pivot_coeff.data() + (size_t)pivot_column * words,
        coeff.data(), (size_t)words * sizeof(uint64_t));
    std::memcpy(
        pivot_rhs.data() + (size_t)pivot_column * block_bytes,
        rhs.data(), block_bytes);
    ++stats.BlockCopies;
    have_pivot[pivot_column] = 1u;
    ++rank;
    return ResidualInsertResult::Inserted;
}

#undef WH2_RESIDUAL_NOINLINE

#if defined(_MSC_VER)
#define WH2_RESIDUAL_COLD_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define WH2_RESIDUAL_COLD_NOINLINE __attribute__((noinline, cold))
#else
#define WH2_RESIDUAL_COLD_NOINLINE
#endif

static WH2_RESIDUAL_COLD_NOINLINE void ReduceResidualRowWithBatchedRhs(
    std::vector<uint8_t>& coeff,
    std::vector<uint8_t>& rhs,
    uint32_t R,
    uint32_t block_bytes,
    const std::vector<uint8_t>& pivot_coeff,
    const std::vector<uint8_t>& pivot_rhs,
    const std::vector<uint8_t>& have_pivot,
    PrecodeSolveStats& stats)
{
    BatchedBlockXorAccumulator rhs_xor(rhs.data(), block_bytes);
    for (uint32_t j = 0; j < R; ++j)
    {
        const uint8_t scale = coeff[j];
        if (scale == 0u || !have_pivot[j]) {
            continue;
        }
        const uint8_t* pivot =
            pivot_coeff.data() + (size_t)j * R;
        if (scale == 1u) {
            gf256_add_mem(
                coeff.data() + j, pivot + j, (int)(R - j));
        }
        else {
            AddScaledResidualCoefficients(
                coeff.data() + j, scale, pivot + j, R - j);
        }
        stats.ResidualCoeffByteOps += R - j;
        const uint8_t* pivot_value =
            pivot_rhs.data() + (size_t)j * block_bytes;
        if (scale == 1u)
        {
            rhs_xor.Add(pivot_value);
            ++stats.BlockXors;
        }
        else {
            AddScaledBlock(
                rhs.data(), scale, pivot_value,
                block_bytes, stats);
        }
    }
    rhs_xor.Flush();
}

#undef WH2_RESIDUAL_COLD_NOINLINE

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)

#if defined(_MSC_VER)
#define WH2_WITNESS_NOINLINE __declspec(noinline)
#elif defined(__ELF__) && (defined(__GNUC__) || defined(__clang__)) && \
    !defined(WIREHAIR_V2_DISABLE_PACKED_RESIDUAL_TEXT_SECTION)
#define WH2_WITNESS_NOINLINE \
    __attribute__((noinline, section(".text.wh2_test_oracle")))
#elif defined(__GNUC__) || defined(__clang__)
#define WH2_WITNESS_NOINLINE __attribute__((noinline))
#else
#define WH2_WITNESS_NOINLINE
#endif

struct CanonicalMixedKernel
{
    uint32_t QuotientRank = 0u;
    uint32_t KernelDimension = 0u;
    uint64_t HashLow = 0u;
    uint64_t HashHigh = 0u;
    std::vector<uint16_t> Basis;
};

WH2_WITNESS_NOINLINE bool CheckedMatrixElements(
    uint32_t rows,
    uint32_t columns,
    size_t& elements)
{
    if (rows != 0u && columns >
            std::numeric_limits<size_t>::max() / rows)
    {
        return false;
    }
    elements = (size_t)rows * columns;
    return true;
}

WH2_WITNESS_NOINLINE bool RrefGF16(
    std::vector<uint16_t>& matrix,
    uint32_t rows,
    uint32_t columns,
    std::vector<uint32_t>& pivot_columns)
{
    size_t elements = 0u;
    if (!CheckedMatrixElements(rows, columns, elements) ||
        matrix.size() != elements)
    {
        return false;
    }
    pivot_columns.clear();
    pivot_columns.reserve(std::min(rows, columns));
    uint32_t rank = 0u;
    for (uint32_t column = 0u;
         column < columns && rank < rows;
         ++column)
    {
        uint32_t selected = rank;
        while (selected < rows &&
               matrix[(size_t)selected * columns + column] == 0u)
        {
            ++selected;
        }
        if (selected == rows) continue;
        if (selected != rank)
        {
            for (uint32_t c = 0u; c < columns; ++c) {
                std::swap(
                    matrix[(size_t)selected * columns + c],
                    matrix[(size_t)rank * columns + c]);
            }
        }
        uint16_t* pivot = matrix.data() + (size_t)rank * columns;
        const uint16_t inverse = GF16InverseInitialized(pivot[column]);
        if (inverse == 0u) return false;
        if (pivot[column] != 1u)
        {
            for (uint32_t c = column; c < columns; ++c) {
                pivot[c] = GF16MultiplyInitialized(pivot[c], inverse);
            }
        }
        for (uint32_t row = 0u; row < rows; ++row)
        {
            if (row == rank) continue;
            uint16_t* destination =
                matrix.data() + (size_t)row * columns;
            const uint16_t scale = destination[column];
            if (scale == 0u) continue;
            for (uint32_t c = column; c < columns; ++c) {
                destination[c] ^=
                    GF16MultiplyInitialized(scale, pivot[c]);
            }
        }
        pivot_columns.push_back(column);
        ++rank;
    }
    return true;
}

WH2_WITNESS_NOINLINE void HashMixedKernelByte(
    uint8_t value,
    uint64_t& low,
    uint64_t& high)
{
    low ^= value;
    low *= UINT64_C(0x100000001b3);
    high ^= (uint8_t)(value ^ UINT8_C(0xa5));
    high *= UINT64_C(0x100000001b3);
    high ^= high >> 29;
}

WH2_WITNESS_NOINLINE void HashMixedKernelU32(
    uint32_t value,
    uint64_t& low,
    uint64_t& high)
{
    for (uint32_t byte = 0u; byte < 4u; ++byte) {
        HashMixedKernelByte(
            (uint8_t)(value >> (8u * byte)), low, high);
    }
}

WH2_WITNESS_NOINLINE void HashCanonicalMixedKernel(
    uint32_t column_count,
    uint32_t quotient_columns,
    uint32_t quotient_rank,
    uint32_t kernel_dimension,
    const std::vector<uint16_t>& basis,
    uint64_t& low,
    uint64_t& high)
{
    low = UINT64_C(0xcbf29ce484222325);
    high = UINT64_C(0x84222325cbf29ce4);
    static const char domain[] = "wh2-mixed-null-v1";
    for (size_t i = 0u; i + 1u < sizeof(domain); ++i) {
        HashMixedKernelByte((uint8_t)domain[i], low, high);
    }
    HashMixedKernelU32(column_count, low, high);
    HashMixedKernelU32(quotient_columns, low, high);
    HashMixedKernelU32(quotient_rank, low, high);
    HashMixedKernelU32(kernel_dimension, low, high);
    for (uint16_t coefficient : basis)
    {
        HashMixedKernelByte((uint8_t)coefficient, low, high);
        HashMixedKernelByte((uint8_t)(coefficient >> 8), low, high);
    }
}

WH2_WITNESS_NOINLINE bool BuildCanonicalMixedKernel(
    const std::vector<uint16_t>& quotient,
    uint32_t quotient_rows,
    uint32_t quotient_columns,
    const std::vector<uint16_t>& full_binary_basis,
    uint32_t column_count,
    CanonicalMixedKernel& output)
{
    size_t quotient_elements = 0u;
    size_t basis_elements = 0u;
    if (!CheckedMatrixElements(
            quotient_rows, quotient_columns, quotient_elements) ||
        !CheckedMatrixElements(
            quotient_columns, column_count, basis_elements) ||
        quotient.size() != quotient_elements ||
        full_binary_basis.size() != basis_elements)
    {
        return false;
    }

    std::vector<uint16_t> reduced = quotient;
    std::vector<uint32_t> quotient_pivots;
    if (!RrefGF16(
            reduced, quotient_rows, quotient_columns, quotient_pivots))
    {
        return false;
    }
    const uint32_t quotient_rank = (uint32_t)quotient_pivots.size();
    const uint32_t kernel_dimension = quotient_columns - quotient_rank;
    size_t null_elements = 0u;
    size_t mapped_elements = 0u;
    if (!CheckedMatrixElements(
            kernel_dimension, quotient_columns, null_elements) ||
        !CheckedMatrixElements(
            kernel_dimension, column_count, mapped_elements))
    {
        return false;
    }
    std::vector<uint16_t> nullspace(null_elements, 0u);
    std::vector<uint8_t> is_pivot(quotient_columns, uint8_t{0});
    for (uint32_t pivot : quotient_pivots) is_pivot[pivot] = 1u;
    uint32_t null_row = 0u;
    for (uint32_t free_column = 0u;
         free_column < quotient_columns;
         ++free_column)
    {
        if (is_pivot[free_column]) continue;
        uint16_t* vector =
            nullspace.data() + (size_t)null_row * quotient_columns;
        vector[free_column] = 1u;
        for (uint32_t row = 0u; row < quotient_rank; ++row) {
            vector[quotient_pivots[row]] =
                reduced[(size_t)row * quotient_columns + free_column];
        }
        ++null_row;
    }
    if (null_row != kernel_dimension) return false;

    for (uint32_t row = 0u; row < quotient_rows; ++row)
    {
        for (uint32_t vector_index = 0u;
             vector_index < kernel_dimension;
             ++vector_index)
        {
            uint16_t sum = 0u;
            const uint16_t* vector = nullspace.data() +
                (size_t)vector_index * quotient_columns;
            for (uint32_t column = 0u;
                 column < quotient_columns;
                 ++column)
            {
                sum ^= GF16MultiplyInitialized(
                    quotient[(size_t)row * quotient_columns + column],
                    vector[column]);
            }
            if (sum != 0u) return false;
        }
    }

    std::vector<uint16_t> mapped(mapped_elements, 0u);
    for (uint32_t vector_index = 0u;
         vector_index < kernel_dimension;
         ++vector_index)
    {
        uint16_t* destination =
            mapped.data() + (size_t)vector_index * column_count;
        const uint16_t* vector = nullspace.data() +
            (size_t)vector_index * quotient_columns;
        for (uint32_t basis_index = 0u;
             basis_index < quotient_columns;
             ++basis_index)
        {
            const uint16_t scale = vector[basis_index];
            if (scale == 0u) continue;
            const uint16_t* source = full_binary_basis.data() +
                (size_t)basis_index * column_count;
            for (uint32_t column = 0u; column < column_count; ++column)
            {
                if (source[column] != 0u) {
                    destination[column] ^= GF16MultiplyInitialized(
                        scale, source[column]);
                }
            }
        }
    }
    std::vector<uint32_t> canonical_pivots;
    if (!RrefGF16(
            mapped, kernel_dimension, column_count, canonical_pivots) ||
        canonical_pivots.size() != kernel_dimension)
    {
        return false;
    }

    CanonicalMixedKernel result;
    result.QuotientRank = quotient_rank;
    result.KernelDimension = kernel_dimension;
    result.Basis.swap(mapped);
    HashCanonicalMixedKernel(
        column_count, quotient_columns, quotient_rank, kernel_dimension,
        result.Basis, result.HashLow, result.HashHigh);
    output = std::move(result);
    return true;
}

WH2_WITNESS_NOINLINE uint32_t Parity64(uint64_t value)
{
    value ^= value >> 32;
    value ^= value >> 16;
    value ^= value >> 8;
    value ^= value >> 4;
    return (UINT16_C(0x6996) >> (value & 15u)) & 1u;
}

WH2_WITNESS_NOINLINE bool VerifyBinaryNullVectors(
    const BinaryEquationArena& rows,
    uint32_t column_count,
    const std::vector<uint16_t>& vectors,
    uint32_t vector_count)
{
    size_t elements = 0u;
    if (!CheckedMatrixElements(
            vector_count, column_count, elements) ||
        vectors.size() != elements)
    {
        return false;
    }
    for (uint32_t vector_index = 0u;
         vector_index < vector_count;
         ++vector_index)
    {
        const uint16_t* vector =
            vectors.data() + (size_t)vector_index * column_count;
        for (uint32_t row = 0u; row < (uint32_t)rows.size(); ++row)
        {
            uint16_t sum = 0u;
            for (uint32_t column : rows[row].Columns)
            {
                if (column >= column_count) return false;
                sum ^= vector[column];
            }
            if (sum != 0u) return false;
        }
    }
    return true;
}

WH2_WITNESS_NOINLINE uint16_t DirectMixedCoefficient(
    const MixedCoefficientRows& coefficients,
    uint32_t subfield_rows,
    uint32_t heavy_row,
    uint32_t column,
    uint32_t first_heavy_column)
{
    if (heavy_row < subfield_rows)
    {
        const uint32_t grouped_gf256_rows =
            ActiveMixedGroupedGF256Rows();
        const uint32_t first_grouped_gf256_row =
            grouped_gf256_rows <= subfield_rows ?
                subfield_rows - grouped_gf256_rows : subfield_rows;
        const uint32_t residue = grouped_gf256_rows != 0u &&
                heavy_row >= first_grouped_gf256_row ?
            ActiveMixedGroupedGF256CoefficientResidue(
                column, first_heavy_column) :
            ActiveMixedCoefficientResidue(column);
        return coefficients.Subfield[heavy_row][
            residue];
    }
    return coefficients.Extension[heavy_row - subfield_rows][
        ActiveMixedExtensionCoefficientResidue(column)];
}

WH2_WITNESS_NOINLINE bool EvaluateDirectMixedSyndromes(
    const PrecodeSystem& system,
    uint32_t column_count,
    const std::vector<uint16_t>& vectors,
    uint32_t vector_count,
    std::vector<uint16_t>* syndromes)
{
    const uint32_t subfield_rows = ActiveMixedGF256Rows();
    const uint32_t extension_rows = ActiveMixedGF16Rows();
    const uint32_t heavy_rows = subfield_rows + extension_rows;
    const MixedCoefficientRows* coefficients = GetMixedCoefficientRows();
    size_t vector_elements = 0u;
    size_t syndrome_elements = 0u;
    if (system.Params.Field != CompletionField::MixedGF256GF16 ||
        system.Params.HeavyRows != heavy_rows ||
        column_count < heavy_rows ||
        ActiveMixedGroupedGF256Rows() > subfield_rows || !coefficients ||
        !CheckedMatrixElements(
            vector_count, column_count, vector_elements) ||
        vectors.size() != vector_elements ||
        !CheckedMatrixElements(
            heavy_rows, vector_count, syndrome_elements))
    {
        return false;
    }
    const uint32_t first_heavy_column = column_count - heavy_rows;
    if (syndromes) syndromes->assign(syndrome_elements, 0u);
    for (uint32_t heavy = 0u; heavy < heavy_rows; ++heavy)
    {
        for (uint32_t vector_index = 0u;
             vector_index < vector_count;
             ++vector_index)
        {
            const uint16_t* vector = vectors.data() +
                (size_t)vector_index * column_count;
            uint16_t sum = 0u;
            for (uint32_t column = 0u;
                 column < column_count;
                 ++column)
            {
                const uint16_t value = vector[column];
                if (value == 0u) continue;
                sum ^= GF16MultiplyInitialized(
                    DirectMixedCoefficient(
                        *coefficients, subfield_rows, heavy, column,
                        first_heavy_column),
                    value);
            }
            if (syndromes) {
                (*syndromes)[(size_t)heavy * vector_count + vector_index] =
                    sum;
            }
            else if (sum != 0u) {
                return false;
            }
        }
    }
    return true;
}

WH2_WITNESS_NOINLINE bool VerifyCompletionNullVectors(
    const PrecodeSystem& system,
    uint32_t column_count,
    const std::vector<uint16_t>& vectors,
    uint32_t vector_count)
{
    return EvaluateDirectMixedSyndromes(
        system, column_count, vectors, vector_count, nullptr);
}

WH2_WITNESS_NOINLINE bool BuildDirectMixedQuotient(
    const PrecodeSystem& system,
    uint32_t column_count,
    const std::vector<uint16_t>& full_binary_basis,
    uint32_t quotient_columns,
    std::vector<uint16_t>& quotient)
{
    return EvaluateDirectMixedSyndromes(
        system, column_count, full_binary_basis, quotient_columns,
        &quotient);
}

WH2_WITNESS_NOINLINE bool BuildFullBinaryNullspaceBasis(
    uint32_t column_count,
    uint32_t inactive_count,
    uint32_t projection_words,
    const std::vector<uint32_t>& inactive_index,
    const std::vector<uint32_t>& inactive_columns,
    const std::vector<uint64_t>& projection,
    const std::vector<uint64_t>& binary_pivot_coeff,
    const std::vector<uint8_t>& binary_have_pivot,
    uint32_t binary_rank,
    std::vector<uint16_t>& full_basis)
{
    if (inactive_count == 0u || binary_rank > inactive_count ||
        projection_words != PackedWordCount(inactive_count) ||
        inactive_index.size() != column_count ||
        inactive_columns.size() != inactive_count ||
        binary_have_pivot.size() != inactive_count)
    {
        return false;
    }
    size_t projection_elements = 0u;
    size_t pivot_elements = 0u;
    if (!CheckedMatrixElements(
            column_count, projection_words, projection_elements) ||
        !CheckedMatrixElements(
            inactive_count, projection_words, pivot_elements) ||
        projection.size() != projection_elements ||
        binary_pivot_coeff.size() != pivot_elements)
    {
        return false;
    }
    const uint32_t quotient_columns = inactive_count - binary_rank;
    if (quotient_columns == 0u ||
        quotient_columns > kMaxMixedNullWitnessQuotientColumns)
    {
        return false;
    }
    std::vector<uint32_t> free_columns;
    free_columns.reserve(quotient_columns);
    uint32_t counted_rank = 0u;
    for (uint32_t column = 0u; column < inactive_count; ++column)
    {
        if (binary_have_pivot[column]) ++counted_rank;
        else free_columns.push_back(column);
        const uint32_t absolute = inactive_columns[column];
        if (absolute >= column_count || inactive_index[absolute] != column) {
            return false;
        }
    }
    if (counted_rank != binary_rank ||
        free_columns.size() != quotient_columns)
    {
        return false;
    }
    const uint64_t tail_mask = (inactive_count & 63u) == 0u ?
        UINT64_MAX :
        (UINT64_C(1) << (inactive_count & 63u)) - UINT64_C(1);
    for (uint32_t column = 0u; column < column_count; ++column) {
        if ((projection[(size_t)column * projection_words +
                        projection_words - 1u] & ~tail_mask) != 0u)
        {
            return false;
        }
    }
    for (uint32_t pivot = 0u; pivot < inactive_count; ++pivot)
    {
        if (!binary_have_pivot[pivot]) continue;
        const uint64_t* relation = binary_pivot_coeff.data() +
            (size_t)pivot * projection_words;
        if ((relation[projection_words - 1u] & ~tail_mask) != 0u ||
            (relation[pivot >> 6] &
                (UINT64_C(1) << (pivot & 63u))) == 0u)
        {
            return false;
        }
        for (uint32_t other = 0u; other < inactive_count; ++other)
        {
            if (other != pivot && binary_have_pivot[other] &&
                (relation[other >> 6] &
                    (UINT64_C(1) << (other & 63u))) != 0u)
            {
                return false;
            }
        }
    }

    size_t inactive_basis_elements = 0u;
    size_t full_basis_elements = 0u;
    if (!CheckedMatrixElements(
            quotient_columns, projection_words,
            inactive_basis_elements) ||
        !CheckedMatrixElements(
            quotient_columns, column_count, full_basis_elements))
    {
        return false;
    }
    std::vector<uint64_t> inactive_basis(
        inactive_basis_elements, uint64_t{0});
    full_basis.assign(full_basis_elements, 0u);
    for (uint32_t basis_index = 0u;
         basis_index < quotient_columns;
         ++basis_index)
    {
        const uint32_t free_column = free_columns[basis_index];
        uint64_t* inactive_vector = inactive_basis.data() +
            (size_t)basis_index * projection_words;
        inactive_vector[free_column >> 6] |=
            UINT64_C(1) << (free_column & 63u);
        for (uint32_t pivot = 0u; pivot < inactive_count; ++pivot)
        {
            if (!binary_have_pivot[pivot]) continue;
            const uint64_t* relation = binary_pivot_coeff.data() +
                (size_t)pivot * projection_words;
            if ((relation[free_column >> 6] &
                    (UINT64_C(1) << (free_column & 63u))) != 0u)
            {
                inactive_vector[pivot >> 6] |=
                    UINT64_C(1) << (pivot & 63u);
            }
        }
        uint16_t* full_vector = full_basis.data() +
            (size_t)basis_index * column_count;
        for (uint32_t inactive = 0u;
             inactive < inactive_count;
             ++inactive)
        {
            if ((inactive_vector[inactive >> 6] &
                    (UINT64_C(1) << (inactive & 63u))) != 0u)
            {
                full_vector[inactive_columns[inactive]] = 1u;
            }
        }
        for (uint32_t column = 0u; column < column_count; ++column)
        {
            if (inactive_index[column] != UINT32_MAX) continue;
            uint32_t parity = 0u;
            const uint64_t* projected = projection.data() +
                (size_t)column * projection_words;
            for (uint32_t word = 0u; word < projection_words; ++word) {
                parity ^= Parity64(projected[word] & inactive_vector[word]);
            }
            full_vector[column] = (uint16_t)(parity & 1u);
        }
    }
    return true;
}

WH2_WITNESS_NOINLINE void CaptureMixedNullWitness(
    const PrecodeSystem& system,
    uint32_t column_count,
    uint32_t inactive_count,
    uint32_t projection_words,
    const BinaryEquationArena& rows,
    const std::vector<uint32_t>& inactive_index,
    const std::vector<uint32_t>& inactive_columns,
    const std::vector<uint64_t>& projection,
    const std::vector<uint64_t>& binary_pivot_coeff,
    const std::vector<uint8_t>& binary_have_pivot,
    uint32_t binary_rank,
    uint32_t expected_quotient_rank,
    MixedNullWitnessDiagnostic* sink) noexcept
{
    if (!sink) return;
    MixedNullWitnessDiagnostic diagnostic;
    diagnostic.ColumnCount = column_count;
    diagnostic.InactiveCount = inactive_count;
    diagnostic.BinaryRank = binary_rank;
    diagnostic.QuotientColumns = inactive_count >= binary_rank ?
        inactive_count - binary_rank : 0u;
    try
    {
        std::vector<uint16_t> full_binary_basis;
        if (!BuildFullBinaryNullspaceBasis(
                column_count, inactive_count, projection_words,
                inactive_index, inactive_columns, projection,
                binary_pivot_coeff, binary_have_pivot, binary_rank,
                full_binary_basis) ||
            !VerifyBinaryNullVectors(
                rows, column_count, full_binary_basis,
                diagnostic.QuotientColumns))
        {
            diagnostic.Status =
                MixedNullWitnessStatus::AlgebraMismatch;
            *sink = std::move(diagnostic);
            return;
        }
        std::vector<uint16_t> quotient;
        if (!BuildDirectMixedQuotient(
                system, column_count, full_binary_basis,
                diagnostic.QuotientColumns, quotient))
        {
            diagnostic.Status =
                MixedNullWitnessStatus::AlgebraMismatch;
            *sink = std::move(diagnostic);
            return;
        }
        CanonicalMixedKernel canonical;
        if (!BuildCanonicalMixedKernel(
                quotient, system.Params.HeavyRows,
                diagnostic.QuotientColumns, full_binary_basis,
                column_count, canonical))
        {
            diagnostic.Status =
                MixedNullWitnessStatus::AlgebraMismatch;
            *sink = std::move(diagnostic);
            return;
        }
        diagnostic.QuotientRank = canonical.QuotientRank;
        diagnostic.KernelDimension = canonical.KernelDimension;
        diagnostic.BasisHashLow = canonical.HashLow;
        diagnostic.BasisHashHigh = canonical.HashHigh;
        diagnostic.CanonicalBasis = std::move(canonical.Basis);
        diagnostic.InactiveMask.resize(column_count, uint8_t{0});
        for (uint32_t column = 0u; column < column_count; ++column) {
            diagnostic.InactiveMask[column] =
                inactive_index[column] != UINT32_MAX ? 1u : 0u;
        }
        diagnostic.QuotientVerified =
            diagnostic.QuotientRank == expected_quotient_rank &&
            diagnostic.KernelDimension ==
                diagnostic.QuotientColumns - diagnostic.QuotientRank &&
            diagnostic.KernelDimension != 0u;
        diagnostic.BinaryVerified = VerifyBinaryNullVectors(
            rows, column_count, diagnostic.CanonicalBasis,
            diagnostic.KernelDimension);
        diagnostic.CompletionVerified = VerifyCompletionNullVectors(
            system, column_count, diagnostic.CanonicalBasis,
            diagnostic.KernelDimension);

        std::vector<uint16_t> canonical_copy = diagnostic.CanonicalBasis;
        std::vector<uint32_t> canonical_pivots;
        uint64_t hash_low = 0u;
        uint64_t hash_high = 0u;
        diagnostic.CanonicalVerified =
            RrefGF16(
                canonical_copy, diagnostic.KernelDimension,
                column_count, canonical_pivots) &&
            canonical_pivots.size() == diagnostic.KernelDimension &&
            canonical_copy == diagnostic.CanonicalBasis;
        HashCanonicalMixedKernel(
            column_count, diagnostic.QuotientColumns,
            diagnostic.QuotientRank, diagnostic.KernelDimension,
            diagnostic.CanonicalBasis, hash_low, hash_high);
        diagnostic.CanonicalVerified =
            diagnostic.CanonicalVerified &&
            hash_low == diagnostic.BasisHashLow &&
            hash_high == diagnostic.BasisHashHigh;
        diagnostic.Status =
            diagnostic.QuotientVerified && diagnostic.BinaryVerified &&
                diagnostic.CompletionVerified &&
                diagnostic.CanonicalVerified ?
            MixedNullWitnessStatus::Captured :
            MixedNullWitnessStatus::AlgebraMismatch;
    }
    catch (const std::bad_alloc&) {
        diagnostic.Status = MixedNullWitnessStatus::AllocationFailure;
        diagnostic.CanonicalBasis.clear();
        diagnostic.InactiveMask.clear();
    }
    catch (...) {
        diagnostic.Status = MixedNullWitnessStatus::InternalError;
        diagnostic.CanonicalBasis.clear();
        diagnostic.InactiveMask.clear();
    }
    *sink = std::move(diagnostic);
}

#undef WH2_WITNESS_NOINLINE

#endif

} // namespace

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
#if defined(_MSC_VER)
#define WH2_TEST_NOINLINE __declspec(noinline)
#elif defined(__ELF__) && (defined(__GNUC__) || defined(__clang__)) && \
    !defined(WIREHAIR_V2_DISABLE_PACKED_RESIDUAL_TEXT_SECTION)
// BUILD_TESTS also links this oracle into the benchmark.  Keep its sizeable
// test-only body out of .text.unlikely, which GNU linkers place before hot
// text and would otherwise perturb cross-revision solver layout.
#define WH2_TEST_NOINLINE \
    __attribute__((noinline, section(".text.wh2_test_oracle")))
#define WH2_TEST_LAMBDA_NOINLINE \
    __attribute__((noinline, section(".text.wh2_test_oracle")))
#elif defined(__GNUC__) || defined(__clang__)
#define WH2_TEST_NOINLINE __attribute__((noinline))
#define WH2_TEST_LAMBDA_NOINLINE __attribute__((noinline))
#else
#define WH2_TEST_NOINLINE
#define WH2_TEST_LAMBDA_NOINLINE
#endif

#if defined(_MSC_VER)
#define WH2_TEST_LAMBDA_NOINLINE
#endif

WH2_TEST_NOINLINE bool CheckMixedQuotientFactorReplayForTesting()
{
    if (!InitializeGF16()) return false;
    try
    {
        uint64_t random_state = UINT64_C(0x243f6a8885a308d3);
        const auto next_random = [&random_state]() -> uint16_t {
            random_state ^= random_state << 13;
            random_state ^= random_state >> 7;
            random_state ^= random_state << 17;
            return (uint16_t)random_state;
        };
        const auto direct_factor = [](
            std::vector<uint16_t>& coefficients,
            std::vector<uint16_t>& rhs,
            uint32_t row_count,
            uint32_t column_count,
            uint32_t rhs_lanes,
            std::array<uint8_t, kMixedCompletionRowsMax>& pivot_rows,
            std::array<uint8_t, kMixedCompletionRowsMax>& row_pivots,
            uint32_t& rank) -> bool
        {
            pivot_rows.fill(UINT8_MAX);
            row_pivots.fill(UINT8_MAX);
            rank = 0u;
            for (uint32_t row = 0u; row < row_count; ++row)
            {
                uint16_t* coeff = column_count == 0u ? nullptr :
                    coefficients.data() + (size_t)row * column_count;
                uint16_t* row_rhs =
                    rhs.data() + (size_t)row * rhs_lanes;
                for (uint32_t column = 0u;
                     column < column_count; ++column)
                {
                    const uint16_t scale = coeff[column];
                    const uint8_t pivot_row = pivot_rows[column];
                    if (scale == 0u || pivot_row == UINT8_MAX) continue;
                    const uint16_t* pivot = coefficients.data() +
                        (size_t)pivot_row * column_count;
                    for (uint32_t k = column; k < column_count; ++k) {
                        coeff[k] ^= GF16MultiplyInitialized(scale, pivot[k]);
                    }
                    const uint16_t* pivot_rhs = rhs.data() +
                        (size_t)pivot_row * rhs_lanes;
                    for (uint32_t lane = 0u; lane < rhs_lanes; ++lane) {
                        row_rhs[lane] ^= GF16MultiplyInitialized(
                            scale, pivot_rhs[lane]);
                    }
                }

                uint32_t pivot_column = column_count;
                for (uint32_t column = 0u;
                     column < column_count; ++column)
                {
                    if (coeff[column] != 0u) {
                        pivot_column = column;
                        break;
                    }
                }
                if (pivot_column == column_count) continue;

                const uint16_t inverse =
                    GF16InverseInitialized(coeff[pivot_column]);
                if (inverse == 0u) return false;
                if (inverse != 1u)
                {
                    for (uint32_t k = pivot_column;
                         k < column_count; ++k)
                    {
                        coeff[k] = GF16MultiplyInitialized(coeff[k], inverse);
                    }
                    for (uint32_t lane = 0u; lane < rhs_lanes; ++lane) {
                        row_rhs[lane] = GF16MultiplyInitialized(
                            row_rhs[lane], inverse);
                    }
                }
                for (uint32_t existing = 0u;
                     existing < column_count; ++existing)
                {
                    const uint8_t existing_row = pivot_rows[existing];
                    if (existing_row == UINT8_MAX) continue;
                    uint16_t* existing_coeff = coefficients.data() +
                        (size_t)existing_row * column_count;
                    const uint16_t scale =
                        existing_coeff[pivot_column];
                    if (scale == 0u) continue;
                    for (uint32_t k = pivot_column;
                         k < column_count; ++k)
                    {
                        existing_coeff[k] ^= GF16MultiplyInitialized(
                            scale, coeff[k]);
                    }
                    uint16_t* existing_rhs = rhs.data() +
                        (size_t)existing_row * rhs_lanes;
                    for (uint32_t lane = 0u; lane < rhs_lanes; ++lane) {
                        existing_rhs[lane] ^= GF16MultiplyInitialized(
                            scale, row_rhs[lane]);
                    }
                }
                pivot_rows[pivot_column] = (uint8_t)row;
                row_pivots[row] = (uint8_t)pivot_column;
                ++rank;
            }
            return true;
        };

        const uint32_t row_counts[] = {
            kMixedGF256Rows + kMixedGF16Rows,
            kMixedCompletionRowsMax
        };
        static const uint32_t rhs_lanes = 3u;
        for (uint32_t row_count : row_counts)
        {
            const uint32_t column_counts[] = {
                0u, 1u, 2u, row_count - 2u, row_count
            };
            for (uint32_t column_count : column_counts)
            {
                for (uint32_t shape = 0u; shape < 4u; ++shape)
                {
                    std::vector<uint16_t> original(
                        (size_t)row_count * column_count, 0u);
                    for (uint16_t& value : original) {
                        value = next_random();
                    }
                    if (shape == 0u)
                    {
                        std::fill(original.begin(), original.end(), 0u);
                        for (uint32_t i = 0u; i < column_count; ++i) {
                            original[(size_t)i * column_count + i] = 1u;
                        }
                    }
                    else if (shape == 1u && column_count >= 2u)
                    {
                        uint16_t scale = next_random();
                        if (scale == 0u) scale = 1u;
                        for (uint32_t row = 0u; row < row_count; ++row) {
                            original[(size_t)row * column_count +
                                column_count - 1u] =
                                GF16MultiplyInitialized(
                                    original[(size_t)row * column_count],
                                    scale);
                        }
                    }
                    else if (shape == 2u && column_count != 0u)
                    {
                        std::fill(original.begin(), original.end(), 0u);
                        original[0] = 1u;
                        if (row_count > 1u) {
                            original[column_count] = UINT16_C(0x1234);
                        }
                        for (uint32_t column = 1u;
                             column < column_count &&
                                 column + 1u < row_count;
                             ++column)
                        {
                            original[(size_t)(column + 1u) * column_count +
                                column] = 1u;
                        }
                    }

                    std::vector<uint16_t> original_rhs(
                        (size_t)row_count * rhs_lanes);
                    for (uint16_t& value : original_rhs) {
                        value = next_random();
                    }
                    std::vector<uint16_t> reference_coeff = original;
                    std::vector<uint16_t> reference_rhs = original_rhs;
                    std::array<
                        uint8_t, kMixedCompletionRowsMax> reference_pivots;
                    std::array<
                        uint8_t, kMixedCompletionRowsMax> reference_row_pivots;
                    uint32_t reference_rank = 0u;
                    if (!direct_factor(
                            reference_coeff, reference_rhs,
                            row_count, column_count, rhs_lanes,
                            reference_pivots, reference_row_pivots,
                            reference_rank))
                    {
                        return false;
                    }

                    std::vector<uint16_t> planned_coeff = original;
                    MixedQuotientFactorization factor;
                    if (!FactorMixedCompletionQuotient(
                            planned_coeff, row_count, column_count, factor) ||
                        planned_coeff != reference_coeff ||
                        factor.Rank != reference_rank ||
                        factor.PivotRows != reference_pivots)
                    {
                        return false;
                    }
                    std::vector<uint16_t> planned_rhs = original_rhs;
                    for (uint32_t row = 0u; row < row_count; ++row)
                    {
                        const MixedQuotientRowPlan& plan = factor.Rows[row];
                        if (plan.PivotColumn != reference_row_pivots[row]) {
                            return false;
                        }
                        uint16_t* row_rhs = planned_rhs.data() +
                            (size_t)row * rhs_lanes;
                        for (uint32_t column = 0u;
                             column < column_count; ++column)
                        {
                            const uint16_t scale =
                                plan.ForwardScales[column];
                            if (scale == 0u) continue;
                            const uint8_t pivot_row =
                                factor.PivotRows[column];
                            if (pivot_row == UINT8_MAX || pivot_row == row) {
                                return false;
                            }
                            const uint16_t* pivot_rhs = planned_rhs.data() +
                                (size_t)pivot_row * rhs_lanes;
                            for (uint32_t lane = 0u;
                                 lane < rhs_lanes; ++lane)
                            {
                                row_rhs[lane] ^= GF16MultiplyInitialized(
                                    scale, pivot_rhs[lane]);
                            }
                        }
                        if (plan.PivotColumn == UINT8_MAX) continue;
                        if (plan.NormalizeScale != 1u) {
                            for (uint32_t lane = 0u;
                                 lane < rhs_lanes; ++lane)
                            {
                                row_rhs[lane] = GF16MultiplyInitialized(
                                    row_rhs[lane], plan.NormalizeScale);
                            }
                        }
                        for (uint32_t existing = 0u;
                             existing < column_count; ++existing)
                        {
                            const uint16_t scale = plan.BackScales[existing];
                            if (scale == 0u) continue;
                            const uint8_t existing_row =
                                factor.PivotRows[existing];
                            if (existing_row == UINT8_MAX ||
                                existing_row == row)
                            {
                                return false;
                            }
                            uint16_t* existing_rhs = planned_rhs.data() +
                                (size_t)existing_row * rhs_lanes;
                            for (uint32_t lane = 0u;
                                 lane < rhs_lanes; ++lane)
                            {
                                existing_rhs[lane] ^=
                                    GF16MultiplyInitialized(
                                        scale, row_rhs[lane]);
                            }
                        }
                    }
                    if (planned_rhs != reference_rhs ||
                        !BuildMixedQuotientTransform(
                            row_count, column_count, factor))
                    {
                        return false;
                    }
                    for (uint32_t row = 0u; row < row_count; ++row)
                    {
                        const uint16_t* transform = factor.Transform.data() +
                            (size_t)row * kMixedCompletionRowsMax;
                        for (uint32_t column = 0u;
                             column < column_count; ++column)
                        {
                            uint16_t value = 0u;
                            for (uint32_t source = 0u;
                                 source < row_count; ++source)
                            {
                                value ^= GF16MultiplyInitialized(
                                    transform[source],
                                    original[(size_t)source * column_count +
                                        column]);
                            }
                            if (value != planned_coeff[
                                    (size_t)row * column_count + column])
                            {
                                return false;
                            }
                        }
                        for (uint32_t lane = 0u; lane < rhs_lanes; ++lane)
                        {
                            uint16_t value = 0u;
                            for (uint32_t source = 0u;
                                 source < row_count; ++source)
                            {
                                value ^= GF16MultiplyInitialized(
                                    transform[source],
                                    original_rhs[
                                        (size_t)source * rhs_lanes + lane]);
                            }
                            if (value != planned_rhs[
                                    (size_t)row * rhs_lanes + lane])
                            {
                                return false;
                            }
                        }
                    }
                }
            }
        }
        return true;
    }
    catch (...) {
        return false;
    }
}

WH2_TEST_NOINLINE bool CheckMixedQuotientDeficientSyndromeForTesting()
{
    if (!InitializeGF16()) return false;
    try
    {
        const uint32_t period = ActiveMixedCoefficientPeriod();
        const uint32_t subfield_rows = ActiveMixedGF256Rows();
        const uint32_t extension_rows = ActiveMixedGF16Rows();
        const uint32_t heavy_rows = subfield_rows + extension_rows;
        const bool independent = ActiveMixedIndependentExtensionResidues();
        const uint32_t grouped_gf256_rows =
            ActiveMixedGroupedGF256Rows();
        const MixedCoefficientRows* rows = GetMixedCoefficientRows();
        const MixedPackedCoefficients* packed_rows =
            GetMixedPackedCoefficients();
        if (!rows || !packed_rows || period == 0u ||
            period > kMixedCoefficientPeriod ||
            heavy_rows > kMixedCompletionRowsMax ||
            grouped_gf256_rows > subfield_rows ||
            (grouped_gf256_rows != 0u && independent))
        {
            return false;
        }
        if (grouped_gf256_rows != 0u)
        {
            // With no non-corner prefix there are no C shifts to cache.  The
            // complete H-column suffix must still project on canonical A, and
            // an empty vector is permitted to expose a null data() pointer.
            const uint32_t inactive_count = 1u;
            const uint32_t projection_words = PackedWordCount(inactive_count);
            std::vector<uint32_t> inactive_index(
                heavy_rows, UINT32_MAX);
            inactive_index[0] = 0u;
            const std::vector<uint64_t> projection(
                (size_t)heavy_rows * projection_words, UINT64_C(0));
            std::vector<uint64_t> dense_projected;
            std::vector<uint64_t> bucket_projected;
            if (!ProjectMixedCompletionCoefficientsByDenseExpansion(
                    heavy_rows, inactive_count, projection_words,
                    inactive_index, projection, *packed_rows,
                    dense_projected) ||
                !ProjectMixedCompletionCoefficientsByResidueBuckets(
                    heavy_rows, inactive_count, projection_words,
                    inactive_index, projection, *packed_rows,
                    nullptr, 0u, bucket_projected) ||
                bucket_projected != dense_projected)
            {
                return false;
            }
        }
        const uint32_t first_grouped_gf256_row =
            subfield_rows - grouped_gf256_rows;
        const uint64_t pair_count_wide = (uint64_t)period * period;
        if (pair_count_wide >= UINT32_MAX - heavy_rows) return false;
        const uint32_t pair_count = (uint32_t)pair_count_wide;
        // Search one more than P^2 non-corner columns, guaranteeing a repeated
        // A/secondary-residue pair.  Reserve a separate final-H suffix so the
        // selected columns remain on C during the actual grouped solve.
        const uint32_t first_heavy_column = pair_count + 1u;
        const uint32_t column_count = first_heavy_column + heavy_rows;
        const auto coefficient_vector = [&](uint32_t column) {
            std::array<uint16_t, kMixedCompletionRowsMax> result = {};
            const uint32_t subfield_residue =
                ActiveMixedCoefficientResidue(column);
            const uint32_t extension_residue = independent ?
                ActiveMixedExtensionCoefficientResidue(column) :
                subfield_residue;
            const uint32_t grouped_residue =
                ActiveMixedGroupedGF256CoefficientResidue(
                    column, first_heavy_column);
            for (uint32_t row = 0u; row < subfield_rows; ++row) {
                const uint32_t residue = grouped_gf256_rows != 0u &&
                        row >= first_grouped_gf256_row ?
                    grouped_residue : subfield_residue;
                result[row] = rows->Subfield[row][residue];
            }
            for (uint32_t row = 0u; row < extension_rows; ++row) {
                result[subfield_rows + row] =
                    rows->Extension[row][extension_residue];
            }
            return result;
        };
        std::vector<uint32_t> first(pair_count, UINT32_MAX);
        uint32_t duplicate0 = UINT32_MAX;
        uint32_t duplicate1 = UINT32_MAX;
        for (uint32_t column = 0u;
             column < first_heavy_column; ++column)
        {
            const uint32_t a = ActiveMixedCoefficientResidue(column);
            const uint32_t b = grouped_gf256_rows != 0u ?
                ActiveMixedGroupedGF256CoefficientResidue(
                    column, first_heavy_column) :
                (independent ?
                    ActiveMixedExtensionCoefficientResidue(column) : a);
            if (a >= period || b >= period) return false;
            const uint32_t key = a * period + b;
            if (first[key] != UINT32_MAX) {
                duplicate0 = first[key];
                duplicate1 = column;
                break;
            }
            first[key] = column;
        }
        if (duplicate0 == UINT32_MAX) return false;

        const std::array<uint16_t, kMixedCompletionRowsMax> base =
            coefficient_vector(duplicate0);
        uint32_t pivot = heavy_rows;
        for (uint32_t row = 0u; row < heavy_rows; ++row) {
            if (base[row] != 0u) {
                pivot = row;
                break;
            }
        }
        if (pivot == heavy_rows) return false;
        const uint16_t inverse = GF16InverseInitialized(base[pivot]);
        if (inverse == 0u) return false;
        uint32_t inconsistent_column = UINT32_MAX;
        for (uint32_t column = 0u;
             column < first_heavy_column; ++column)
        {
            if (column == duplicate0 || column == duplicate1) continue;
            const std::array<uint16_t, kMixedCompletionRowsMax> candidate =
                coefficient_vector(column);
            const uint16_t scale = GF16MultiplyInitialized(
                candidate[pivot], inverse);
            bool multiple = true;
            for (uint32_t row = 0u; row < heavy_rows; ++row) {
                if (candidate[row] != GF16MultiplyInitialized(
                        scale, base[row]))
                {
                    multiple = false;
                    break;
                }
            }
            if (!multiple) {
                inconsistent_column = column;
                break;
            }
        }
        if (inconsistent_column == UINT32_MAX) return false;

        PrecodeSystem system;
        system.Params = MakeMixedParams(
            column_count, UINT64_C(0x13198a2e03707344));
        const uint32_t block_bytes = 2u;
        const int previous_fast_path_mode =
            TinyMixedFastPathTestModeValue;
        struct RestoreTinyFastPathMode
        {
            int Previous;
            ~RestoreTinyFastPathMode()
            {
                TinyMixedFastPathTestModeValue = Previous;
            }
        } restore_fast_path_mode = {previous_fast_path_mode};

        // The duplicate coefficient pair also gives a one-variable,
        // full-rank system with a nonzero solution.  Exercise it under every
        // schedule used by this oracle so secondary-residue coefficient
        // assembly and the tiny solver's successful writeback cannot escape
        // the forced on/off comparison.
        {
            const uint32_t success_inactive_count = 1u;
            const uint32_t success_projection_words =
                PackedWordCount(success_inactive_count);
            std::vector<uint32_t> success_inactive_index(
                column_count, UINT32_MAX);
            success_inactive_index[duplicate0] = 0u;
            const std::vector<uint32_t> success_inactive_columns = {
                duplicate0
            };
            const std::vector<uint64_t> success_projection(
                (size_t)column_count * success_projection_words,
                UINT64_C(0));
            const std::vector<uint64_t> success_binary_coeff(
                (size_t)success_inactive_count *
                    success_projection_words,
                UINT64_C(0));
            const std::vector<uint8_t> success_binary_rhs(
                (size_t)success_inactive_count * block_bytes, uint8_t{0});
            const std::vector<uint8_t> success_have_pivot(
                success_inactive_count, uint8_t{0});
            std::vector<uint8_t> success_initial_values(
                (size_t)column_count * block_bytes, uint8_t{0});
            success_initial_values[
                (size_t)duplicate1 * block_bytes] = 0x5au;
            success_initial_values[
                (size_t)duplicate1 * block_bytes + 1u] = 0xa5u;
            const auto run_success = [&](
                int mode,
                PrecodeSolveStats& stats,
                std::vector<uint8_t>& values) -> WirehairResult
            {
                TinyMixedFastPathTestModeValue = mode;
                values = success_initial_values;
                return SolveMixedCompletionQuotient(
                    system, column_count, success_inactive_count,
                    success_projection_words, block_bytes,
                    success_inactive_index, success_inactive_columns,
                    success_projection, success_binary_coeff,
                    success_binary_rhs, success_have_pivot, 0u,
                    values, stats);
            };
            PrecodeSolveStats success_general_stats;
            PrecodeSolveStats success_tiny_stats;
            std::vector<uint8_t> success_general_values;
            std::vector<uint8_t> success_tiny_values;
            if (run_success(
                    -1, success_general_stats,
                    success_general_values) != Wirehair_Success ||
                run_success(
                    1, success_tiny_stats,
                    success_tiny_values) != Wirehair_Success ||
                success_general_stats.TinyMixedFastPathAcceptances != 0u ||
                success_tiny_stats.TinyMixedFastPathAcceptances != 1u ||
                success_tiny_values != success_general_values ||
                success_tiny_stats.ResidualRows !=
                    success_general_stats.ResidualRows ||
                success_tiny_stats.BinaryResidualRank !=
                    success_general_stats.BinaryResidualRank ||
                success_tiny_stats.ResidualRank !=
                    success_general_stats.ResidualRank ||
                std::memcmp(
                    success_tiny_values.data() +
                        (size_t)duplicate0 * block_bytes,
                    success_initial_values.data() +
                        (size_t)duplicate1 * block_bytes,
                    block_bytes) != 0)
            {
                return false;
            }
        }

        const uint32_t inactive_count = 2u;
        const uint32_t projection_words = PackedWordCount(inactive_count);
        std::vector<uint32_t> inactive_index(column_count, UINT32_MAX);
        inactive_index[duplicate0] = 0u;
        inactive_index[duplicate1] = 1u;
        const std::vector<uint32_t> inactive_columns = {
            duplicate0, duplicate1
        };
        const std::vector<uint64_t> projection(
            (size_t)column_count * projection_words, UINT64_C(0));
        const std::vector<uint64_t> binary_coeff(
            (size_t)inactive_count * projection_words, UINT64_C(0));
        const std::vector<uint8_t> binary_rhs(
            (size_t)inactive_count * block_bytes, uint8_t{0});
        const std::vector<uint8_t> have_pivot(inactive_count, uint8_t{0});
        const std::vector<uint8_t> zero_values(
            (size_t)column_count * block_bytes, uint8_t{0});
        const auto run_quotient = [&](
            int mode,
            const std::vector<uint8_t>& initial_values,
            PrecodeSolveStats& stats,
            std::vector<uint8_t>& values) -> WirehairResult
        {
            TinyMixedFastPathTestModeValue = mode;
            values = initial_values;
            return SolveMixedCompletionQuotient(
                system, column_count, inactive_count, projection_words,
                block_bytes, inactive_index, inactive_columns, projection,
                binary_coeff, binary_rhs, have_pivot, 0u, values, stats);
        };

        PrecodeSolveStats deficient_general;
        PrecodeSolveStats deficient_tiny;
        std::vector<uint8_t> deficient_general_values;
        std::vector<uint8_t> deficient_tiny_values;
        if (run_quotient(
                -1, zero_values, deficient_general,
                deficient_general_values) != Wirehair_NeedMore ||
            run_quotient(
                1, zero_values, deficient_tiny,
                deficient_tiny_values) != Wirehair_NeedMore ||
            deficient_general.TinyMixedFastPathAcceptances != 0u ||
            deficient_tiny.TinyMixedFastPathAcceptances != 1u ||
            deficient_general_values != zero_values ||
            deficient_tiny_values != zero_values ||
            deficient_general.BinaryResidualRank != 0u ||
            deficient_general.ResidualRank != 1u ||
            deficient_tiny.ResidualRows !=
                deficient_general.ResidualRows ||
            deficient_tiny.BinaryResidualRank !=
                deficient_general.BinaryResidualRank ||
            deficient_tiny.ResidualRank != deficient_general.ResidualRank)
        {
            return false;
        }
        std::vector<uint8_t> inconsistent_values = zero_values;
        inconsistent_values[
            (size_t)inconsistent_column * block_bytes] = 1u;
        PrecodeSolveStats inconsistent_general;
        PrecodeSolveStats inconsistent_tiny;
        std::vector<uint8_t> inconsistent_general_values;
        std::vector<uint8_t> inconsistent_tiny_values;
        return
            run_quotient(
                -1, inconsistent_values,
                inconsistent_general,
                inconsistent_general_values) == Wirehair_Error &&
            run_quotient(
                1, inconsistent_values,
                inconsistent_tiny,
                inconsistent_tiny_values) == Wirehair_Error &&
            inconsistent_general.TinyMixedFastPathAcceptances == 0u &&
            inconsistent_tiny.TinyMixedFastPathAcceptances == 1u &&
            inconsistent_general_values == inconsistent_values &&
            inconsistent_tiny_values == inconsistent_values &&
            inconsistent_tiny.ResidualRows ==
                inconsistent_general.ResidualRows &&
            inconsistent_tiny.BinaryResidualRank ==
                inconsistent_general.BinaryResidualRank &&
            inconsistent_tiny.ResidualRank ==
                inconsistent_general.ResidualRank;
    }
    catch (...) {
        return false;
    }
}

WH2_TEST_NOINLINE bool CheckMixedNullWitnessCanonicalizationForTesting()
{
    if (!InitializeGF16()) return false;
    try
    {
        const uint32_t q = 3u;
        const uint32_t L = 5u;
        const uint16_t a = UINT16_C(0x0102);
        const uint16_t b = UINT16_C(0x3456);
        const std::vector<uint16_t> binary_basis = {
            1u, 0u, 1u, 0u, 0u,
            0u, 1u, 1u, 0u, 1u,
            0u, 0u, 0u, 1u, 1u
        };
        const std::vector<uint16_t> quotient = {
            1u, 0u, a,
            0u, 1u, b
        };
        CanonicalMixedKernel reference;
        if (!BuildCanonicalMixedKernel(
                quotient, 2u, q, binary_basis, L, reference) ||
            reference.QuotientRank != 2u ||
            reference.KernelDimension != 1u ||
            reference.Basis.size() != L ||
            reference.Basis[0] != 1u)
        {
            return false;
        }
        const uint16_t inverse_a = GF16InverseInitialized(a);
        const uint16_t expected[] = {
            1u,
            GF16MultiplyInitialized(b, inverse_a),
            GF16MultiplyInitialized((uint16_t)(a ^ b), inverse_a),
            inverse_a,
            GF16MultiplyInitialized((uint16_t)(b ^ 1u), inverse_a)
        };
        if (!std::equal(
                reference.Basis.begin(), reference.Basis.end(), expected))
        {
            return false;
        }

        // B' = B*T and Q' = Q*T for one fixed invertible binary T.
        std::vector<uint16_t> transformed_basis((size_t)q * L, 0u);
        std::vector<uint16_t> transformed_quotient(2u * q, 0u);
        const uint32_t combinations[q][2] = {
            {0u, UINT32_MAX}, {0u, 1u}, {1u, 2u}
        };
        for (uint32_t column = 0u; column < q; ++column)
        {
            for (uint32_t term_index = 0u; term_index < 2u; ++term_index)
            {
                const uint32_t term = combinations[column][term_index];
                if (term == UINT32_MAX) continue;
                for (uint32_t c = 0u; c < L; ++c) {
                    transformed_basis[(size_t)column * L + c] ^=
                        binary_basis[(size_t)term * L + c];
                }
                for (uint32_t row = 0u; row < 2u; ++row) {
                    transformed_quotient[(size_t)row * q + column] ^=
                        quotient[(size_t)row * q + term];
                }
            }
        }
        std::swap_ranges(
            transformed_quotient.begin(),
            transformed_quotient.begin() + q,
            transformed_quotient.begin() + q);
        CanonicalMixedKernel transformed;
        if (!BuildCanonicalMixedKernel(
                transformed_quotient, 2u, q, transformed_basis, L,
                transformed) ||
            transformed.Basis != reference.Basis ||
            transformed.HashLow != reference.HashLow ||
            transformed.HashHigh != reference.HashHigh)
        {
            return false;
        }

        CanonicalMixedKernel dimension_two;
        const std::vector<uint16_t> one_row = {1u, a, b};
        if (!BuildCanonicalMixedKernel(
                one_row, 1u, q, binary_basis, L, dimension_two) ||
            dimension_two.QuotientRank != 1u ||
            dimension_two.KernelDimension != 2u ||
            dimension_two.Basis.size() != 2u * L)
        {
            return false;
        }
        std::vector<uint16_t> transformed_one_row(q, 0u);
        for (uint32_t column = 0u; column < q; ++column)
        {
            uint16_t value = 0u;
            for (uint32_t term_index = 0u; term_index < 2u; ++term_index)
            {
                const uint32_t term = combinations[column][term_index];
                if (term != UINT32_MAX) value ^= one_row[term];
            }
            transformed_one_row[column] =
                GF16MultiplyInitialized(value, b);
        }
        CanonicalMixedKernel transformed_dimension_two;
        if (!BuildCanonicalMixedKernel(
                transformed_one_row, 1u, q, transformed_basis, L,
                transformed_dimension_two) ||
            transformed_dimension_two.Basis != dimension_two.Basis ||
            transformed_dimension_two.HashLow != dimension_two.HashLow ||
            transformed_dimension_two.HashHigh != dimension_two.HashHigh)
        {
            return false;
        }
        CanonicalMixedKernel full_rank;
        const std::vector<uint16_t> identity = {
            1u, 0u, 0u,
            0u, 1u, 0u,
            0u, 0u, 1u
        };
        if (!BuildCanonicalMixedKernel(
                identity, q, q, binary_basis, L, full_rank) ||
            full_rank.QuotientRank != q ||
            full_rank.KernelDimension != 0u ||
            !full_rank.Basis.empty())
        {
            return false;
        }
        std::vector<uint16_t> malformed_basis = binary_basis;
        malformed_basis.pop_back();
        return !BuildCanonicalMixedKernel(
            quotient, 2u, q, malformed_basis, L, full_rank);
    }
    catch (...) {
        return false;
    }
}

WH2_TEST_NOINLINE bool CheckPackedBinaryResidualOracleForTesting()
{
    if (gf256_init() != 0 ||
        PackedWordCount(0u) != 0u ||
        PackedWordCount(63u) != 1u ||
        PackedWordCount(64u) != 1u ||
        PackedWordCount(65u) != 2u ||
        PackedWordCount(UINT32_MAX) != UINT32_C(67108864))
    {
        return false;
    }

    try
    {
        const uint32_t widths[] = {
            63u, 64u, 65u,
            127u, 128u, 129u,
            191u, 192u, 193u
        };
        const uint32_t block_sizes[] = {
            2u,
            511u, 512u, 513u,
            543u,
            4096u
        };
        for (uint32_t R : widths)
        {
            const uint32_t words = PackedWordCount(R);
            const uint64_t tail_mask = (R & 63u) == 0u ?
                UINT64_MAX :
                (UINT64_C(1) << (R & 63u)) - UINT64_C(1);
            for (uint32_t block_bytes : block_sizes)
            {
                // One complete set around the first word boundary covers the
                // direct/batched threshold, exact SIMD widths, short tails,
                // and a 4-KiB payload.  Wider cases focus on coefficient
                // packing with the compact payload.
                if (R > 65u && block_bytes != 2u) continue;

                uint64_t random_state =
                    UINT64_C(0x6a09e667f3bcc909) ^
                    ((uint64_t)R << 32) ^ block_bytes;
                const auto next_random = [&random_state]()
                    WH2_TEST_LAMBDA_NOINLINE -> uint64_t {
                    random_state += UINT64_C(0x9e3779b97f4a7c15);
                    uint64_t value = random_state;
                    value = (value ^ (value >> 30)) *
                        UINT64_C(0xbf58476d1ce4e5b9);
                    value = (value ^ (value >> 27)) *
                        UINT64_C(0x94d049bb133111eb);
                    return value ^ (value >> 31);
                };

                std::vector<uint8_t> solution((size_t)R * block_bytes);
                for (uint8_t& value : solution) {
                    value = (uint8_t)next_random();
                }

                // A shuffled random upper-triangular basis is guaranteed full
                // rank but still exercises cross-word pivots and dense RREF
                // elimination.  Additional random rows below are then known
                // dependencies of this complete basis.
                std::vector<std::vector<uint64_t> > basis(
                    R, std::vector<uint64_t>(words, uint64_t{0}));
                for (uint32_t row = 0u; row < R; ++row)
                {
                    basis[row][row >> 6] |=
                        UINT64_C(1) << (row & 63u);
                    for (uint32_t column = row + 1u;
                         column < R; ++column)
                    {
                        if ((next_random() & 1u) != 0u) {
                            basis[row][column >> 6] |=
                                UINT64_C(1) << (column & 63u);
                        }
                    }
                }
                std::vector<uint32_t> order(R);
                for (uint32_t i = 0u; i < R; ++i) order[i] = i;
                for (uint32_t count = R; count > 1u; --count) {
                    std::swap(
                        order[count - 1u],
                        order[(uint32_t)(next_random() % count)]);
                }

                std::vector<uint8_t> byte_pivots((size_t)R * R, 0u);
                std::vector<uint64_t> packed_pivots(
                    (size_t)R * words, uint64_t{0});
                std::vector<uint8_t> byte_pivot_rhs(
                    (size_t)R * block_bytes, 0u);
                std::vector<uint8_t> packed_pivot_rhs(
                    (size_t)R * block_bytes, 0u);
                std::vector<uint8_t> byte_have(R, 0u);
                std::vector<uint8_t> packed_have(R, 0u);
                uint32_t byte_rank = 0u;
                uint32_t packed_rank = 0u;
                PrecodeSolveStats byte_stats = {};
                PrecodeSolveStats packed_stats = {};

                const auto states_match = [&]()
                    WH2_TEST_LAMBDA_NOINLINE -> bool {
                    if (byte_rank != packed_rank ||
                        byte_have != packed_have ||
                        byte_pivot_rhs != packed_pivot_rhs ||
                        byte_stats.BlockXors != packed_stats.BlockXors ||
                        byte_stats.BlockMulAdds !=
                            packed_stats.BlockMulAdds)
                    {
                        return false;
                    }
                    for (uint32_t row = 0u; row < R; ++row)
                    {
                        const uint64_t* packed = packed_pivots.data() +
                            (size_t)row * words;
                        if ((packed[words - 1u] & ~tail_mask) != 0u) {
                            return false;
                        }
                        for (uint32_t column = 0u;
                             column < R; ++column)
                        {
                            const uint8_t bit = (uint8_t)(
                                (packed[column >> 6] >>
                                    (column & 63u)) & 1u);
                            if (byte_pivots[(size_t)row * R + column] != bit) {
                                return false;
                            }
                        }
                    }
                    return true;
                };

                const auto insert_row = [&, tail_mask](
                    const std::vector<uint64_t>& row,
                    bool corrupt_rhs,
                    ResidualInsertResult expected)
                    WH2_TEST_LAMBDA_NOINLINE -> bool
                {
                    std::vector<uint8_t> byte_coeff(R, 0u);
                    std::vector<uint64_t> packed_coeff = row;
                    std::vector<uint8_t> rhs(block_bytes, 0u);
                    for (uint32_t word_index = 0u;
                         word_index < words; ++word_index)
                    {
                        uint64_t word = row[word_index];
                        while (word != 0u)
                        {
                            const uint32_t bit =
                                wirehair::NonzeroLowestBitIndex64(word);
                            const uint32_t column =
                                (word_index << 6) + bit;
                            if (column < R)
                            {
                                byte_coeff[column] = 1u;
                                const uint8_t* source = solution.data() +
                                    (size_t)column * block_bytes;
                                for (uint32_t i = 0u;
                                     i < block_bytes; ++i)
                                {
                                    rhs[i] ^= source[i];
                                }
                            }
                            word &= word - 1u;
                        }
                    }
                    if (corrupt_rhs) rhs[0] ^= 1u;
                    std::vector<uint8_t> byte_rhs = rhs;
                    std::vector<uint8_t> packed_rhs = rhs;
                    if (tail_mask != UINT64_MAX) {
                        packed_coeff[words - 1u] |= ~tail_mask;
                    }

                    const ResidualInsertResult byte_result =
                        InsertResidualRow(
                            byte_coeff, byte_rhs, R, block_bytes,
                            byte_pivots, byte_pivot_rhs, byte_have,
                            byte_rank, byte_stats, true,
                            kBatchedResidualRhsMinBlockBytes);
                    const ResidualInsertResult packed_result =
                        InsertPackedBinaryResidualRow(
                            packed_coeff, packed_rhs, R, words, block_bytes,
                            packed_pivots, packed_pivot_rhs, packed_have,
                            packed_rank, packed_stats,
                            kBatchedResidualRhsMinBlockBytes);
                    return byte_result == expected &&
                        packed_result == expected &&
                        (packed_coeff[words - 1u] & ~tail_mask) == 0u &&
                        states_match();
                };

                for (uint32_t index : order) {
                    if (!insert_row(
                            basis[index], false,
                            ResidualInsertResult::Inserted))
                    {
                        return false;
                    }
                }
                if (byte_rank != R || packed_rank != R) return false;
                for (uint32_t row = 0u; row < R; ++row)
                {
                    if (!byte_have[row] ||
                        std::memcmp(
                            byte_pivot_rhs.data() +
                                (size_t)row * block_bytes,
                            solution.data() + (size_t)row * block_bytes,
                            block_bytes) != 0)
                    {
                        return false;
                    }
                    for (uint32_t column = 0u;
                         column < R; ++column)
                    {
                        if (byte_pivots[(size_t)row * R + column] !=
                            (row == column ? 1u : 0u))
                        {
                            return false;
                        }
                    }
                }

                for (uint32_t extra = 0u; extra < 8u; ++extra)
                {
                    std::vector<uint64_t> dependent(words);
                    for (uint64_t& word : dependent) word = next_random();
                    dependent[words - 1u] &= tail_mask;
                    if (!insert_row(
                            dependent, false,
                            ResidualInsertResult::Dependent))
                    {
                        return false;
                    }
                }
                if (!insert_row(
                        basis[order[0]], true,
                        ResidualInsertResult::Inconsistent))
                {
                    return false;
                }
            }
        }
    }
    catch (...) {
        return false;
    }
    return true;
}

WH2_TEST_NOINLINE bool CheckMixedRhsFusionOracleForTesting()
{
    if (gf256_init() != 0) return false;

    try
    {
        // Include scalar tails and both sides of the major SIMD/fusion policy
        // boundaries.  The primitive is byte-oriented even though production
        // mixed completion admits only even payload widths.
        const uint32_t block_sizes[] = {
            1u, 2u, 3u,
            15u, 16u, 17u,
            31u, 32u, 33u,
            63u, 64u, 65u,
            127u, 128u, 129u,
            255u, 256u, 257u,
            511u, 512u, 513u,
            1279u, 1280u, 1281u,
            4095u, 4096u, 4097u
        };
        uint64_t random_state = UINT64_C(0x243f6a8885a308d3);
        const auto next_random = [&random_state]() -> uint64_t {
            random_state += UINT64_C(0x9e3779b97f4a7c15);
            uint64_t value = random_state;
            value = (value ^ (value >> 30)) *
                UINT64_C(0xbf58476d1ce4e5b9);
            value = (value ^ (value >> 27)) *
                UINT64_C(0x94d049bb133111eb);
            return value ^ (value >> 31);
        };

        for (uint32_t block_bytes : block_sizes)
        {
            // q=0..12 covers shipping H.  q=13..16 exercises the test-hook
            // extension and the 16-source set-form/additive spill boundary.
            for (uint32_t q = 0u; q <= 16u; ++q)
            {
                std::vector<uint8_t> pivot_rhs(block_bytes);
                std::vector<uint8_t> values(
                    (size_t)(q + 1u) * block_bytes);
                for (uint8_t& byte : pivot_rhs) {
                    byte = (uint8_t)next_random();
                }
                for (uint8_t& byte : values) {
                    byte = (uint8_t)next_random();
                }
                std::vector<uint32_t> free_columns(q);
                std::vector<uint32_t> inactive_columns(q + 1u);
                for (uint32_t i = 0u; i < q; ++i) {
                    free_columns[i] = i;
                    inactive_columns[i] = q - i - 1u;
                }
                inactive_columns[q] = q;

                const uint64_t all = q == 0u ? UINT64_C(0) :
                    (UINT64_C(1) << q) - UINT64_C(1);
                const uint64_t patterns[] = {
                    UINT64_C(0), all,
                    all & UINT64_C(0x5555555555555555),
                    all & UINT64_C(0xaaaaaaaaaaaaaaaa)
                };
                for (uint64_t relation_word : patterns)
                {
                    std::vector<uint8_t> expected = pivot_rhs;
                    PrecodeSolveStats expected_stats = {};
                    expected_stats.BlockXors = 19u;
                    expected_stats.BlockMulAdds = 23u;
                    for (uint32_t free_column : free_columns)
                    {
                        if ((relation_word &
                                (UINT64_C(1) << free_column)) == 0u)
                        {
                            continue;
                        }
                        gf256_add_mem(
                            expected.data(),
                            values.data() +
                                (size_t)inactive_columns[free_column] *
                                    block_bytes,
                            (int)block_bytes);
                        ++expected_stats.BlockXors;
                    }

                    std::vector<uint8_t> actual_values = values;
                    std::memset(
                        actual_values.data() + (size_t)q * block_bytes,
                        0xa5, block_bytes);
                    PrecodeSolveStats actual_stats = {};
                    actual_stats.BlockXors = 19u;
                    actual_stats.BlockMulAdds = 23u;
                    InitializeMixedBinaryPivotValue(
                        actual_values.data() + (size_t)q * block_bytes,
                        block_bytes, pivot_rhs.data(), &relation_word,
                        free_columns, inactive_columns, actual_values,
                        actual_stats);
                    if (std::memcmp(
                            expected.data(),
                            actual_values.data() + (size_t)q * block_bytes,
                            block_bytes) != 0 ||
                        actual_stats.BlockXors !=
                            expected_stats.BlockXors ||
                        actual_stats.BlockMulAdds !=
                            expected_stats.BlockMulAdds)
                    {
                        return false;
                    }
                }

                // Exercise packet-data and null-data initialization with the
                // same number of logical known-constant XORs as production.
                std::vector<uint8_t> packet_data(block_bytes);
                std::vector<uint8_t> sources((size_t)q * block_bytes);
                for (uint8_t& byte : packet_data) {
                    byte = (uint8_t)next_random();
                }
                for (uint8_t& byte : sources) {
                    byte = (uint8_t)next_random();
                }
                for (uint32_t data_mode = 0u; data_mode < 2u; ++data_mode)
                {
                    const uint8_t* data = data_mode != 0u ?
                        packet_data.data() : nullptr;
                    std::vector<uint8_t> expected(block_bytes, 0u);
                    if (data) {
                        std::memcpy(expected.data(), data, block_bytes);
                    }
                    BatchedBlockXorAccumulator expected_xor(
                        expected.data(), block_bytes);
                    for (uint32_t i = 0u; i < q; ++i) {
                        expected_xor.Add(
                            sources.data() + (size_t)i * block_bytes);
                    }
                    expected_xor.Flush();

                    std::vector<uint8_t> actual(block_bytes, 0xa5u);
                    BatchedBlockXorInitializer actual_xor(
                        actual.data(), block_bytes, data);
                    for (uint32_t i = 0u; i < q; ++i) {
                        actual_xor.Add(
                            sources.data() + (size_t)i * block_bytes);
                    }
                    actual_xor.Flush();
                    if (actual != expected) return false;
                }
            }

            // The first set-form batch also has to transition correctly to
            // additive batches for unusually dense unused equations.
            const uint32_t source_count = 33u;
            std::vector<uint8_t> sources(
                (size_t)source_count * block_bytes);
            std::vector<uint8_t> data(block_bytes);
            for (uint8_t& byte : sources) {
                byte = (uint8_t)next_random();
            }
            for (uint8_t& byte : data) {
                byte = (uint8_t)next_random();
            }
            for (uint32_t data_mode = 0u; data_mode < 2u; ++data_mode)
            {
                const uint8_t* first = data_mode != 0u ?
                    data.data() : nullptr;
                std::vector<uint8_t> expected(block_bytes, 0u);
                if (first) {
                    std::memcpy(expected.data(), first, block_bytes);
                }
                BatchedBlockXorAccumulator expected_xor(
                    expected.data(), block_bytes);
                std::vector<uint8_t> actual(block_bytes, 0x5au);
                BatchedBlockXorInitializer actual_xor(
                    actual.data(), block_bytes, first);
                for (uint32_t i = 0u; i < source_count; ++i)
                {
                    const uint8_t* source =
                        sources.data() + (size_t)i * block_bytes;
                    expected_xor.Add(source);
                    actual_xor.Add(source);
                }
                expected_xor.Flush();
                actual_xor.Flush();
                if (actual != expected) return false;
            }
        }

        // Real inactive-coordinate layouts are not generally contiguous.
        // Cross all packed relation-word boundaries used by larger test-hook
        // systems and verify that the helper indexes the inverse map, not the
        // ordinal position in free_columns.
        {
            const uint32_t block_bytes = 257u;
            const uint32_t inactive_count = 130u;
            const uint32_t free_count = 7u;
            const uint32_t destination_index = free_count;
            const uint32_t free_data[free_count] = {
                0u, 1u, 63u, 64u, 65u, 127u, 129u
            };
            const std::vector<uint32_t> free_columns(
                free_data, free_data + free_count);
            std::vector<uint32_t> inactive_columns(inactive_count, 0u);
            for (uint32_t i = 0u; i < free_count; ++i) {
                inactive_columns[free_columns[i]] = free_count - i - 1u;
            }
            std::vector<uint8_t> values(
                (size_t)(free_count + 1u) * block_bytes);
            std::vector<uint8_t> pivot_rhs(block_bytes);
            for (uint8_t& byte : values) {
                byte = (uint8_t)next_random();
            }
            for (uint8_t& byte : pivot_rhs) {
                byte = (uint8_t)next_random();
            }
            std::vector<uint64_t> relation(
                PackedWordCount(inactive_count), UINT64_C(0));
            const uint32_t selected_data[] = {
                0u, 63u, 64u, 127u, 129u
            };
            for (uint32_t selected : selected_data) {
                relation[selected >> 6] |=
                    UINT64_C(1) << (selected & 63u);
            }
            std::vector<uint8_t> expected = pivot_rhs;
            PrecodeSolveStats expected_stats = {};
            expected_stats.BlockXors = 11u;
            for (uint32_t free_column : free_columns)
            {
                if ((relation[free_column >> 6] &
                        (UINT64_C(1) << (free_column & 63u))) == 0u)
                {
                    continue;
                }
                gf256_add_mem(
                    expected.data(),
                    values.data() +
                        (size_t)inactive_columns[free_column] * block_bytes,
                    (int)block_bytes);
                ++expected_stats.BlockXors;
            }
            std::vector<uint8_t> actual_values = values;
            std::memset(
                actual_values.data() +
                    (size_t)destination_index * block_bytes,
                0xa5, block_bytes);
            PrecodeSolveStats actual_stats = {};
            actual_stats.BlockXors = 11u;
            InitializeMixedBinaryPivotValue(
                actual_values.data() +
                    (size_t)destination_index * block_bytes,
                block_bytes, pivot_rhs.data(), relation.data(),
                free_columns, inactive_columns, actual_values,
                actual_stats);
            if (std::memcmp(
                    expected.data(),
                    actual_values.data() +
                        (size_t)destination_index * block_bytes,
                    block_bytes) != 0 ||
                actual_stats.BlockXors != expected_stats.BlockXors ||
                actual_stats.BlockMulAdds != expected_stats.BlockMulAdds)
            {
                return false;
            }
        }

        // Feed fused and legacy initializations through packed residual
        // insertion.  This locks down independent/dependent/inconsistent
        // classification as well as exact pivot bytes and logical receipts.
        const uint32_t R = 13u;
        const uint32_t words = PackedWordCount(R);
        const uint32_t block_bytes = 513u;
        for (uint32_t data_mode = 0u; data_mode < 2u; ++data_mode)
        {
            std::vector<uint8_t> data(block_bytes);
            std::vector<uint8_t> source(block_bytes);
            for (uint8_t& byte : data) {
                byte = (uint8_t)next_random();
            }
            for (uint8_t& byte : source) {
                byte = (uint8_t)next_random();
            }
            const uint8_t* first = data_mode != 0u ? data.data() : nullptr;
            std::vector<uint8_t> legacy_rhs(block_bytes, 0u);
            if (first) {
                std::memcpy(legacy_rhs.data(), first, block_bytes);
            }
            BatchedBlockXorAccumulator legacy_xor(
                legacy_rhs.data(), block_bytes);
            legacy_xor.Add(source.data());
            legacy_xor.Flush();
            std::vector<uint8_t> fused_rhs(block_bytes, 0xa5u);
            BatchedBlockXorInitializer fused_xor(
                fused_rhs.data(), block_bytes, first);
            fused_xor.Add(source.data());
            fused_xor.Flush();
            if (legacy_rhs != fused_rhs) return false;

            std::vector<uint64_t> legacy_pivots((size_t)R * words, 0u);
            std::vector<uint64_t> fused_pivots((size_t)R * words, 0u);
            std::vector<uint8_t> legacy_pivot_rhs(
                (size_t)R * block_bytes, 0u);
            std::vector<uint8_t> fused_pivot_rhs(
                (size_t)R * block_bytes, 0u);
            std::vector<uint8_t> legacy_have(R, 0u);
            std::vector<uint8_t> fused_have(R, 0u);
            uint32_t legacy_rank = 0u;
            uint32_t fused_rank = 0u;
            PrecodeSolveStats legacy_stats = {};
            PrecodeSolveStats fused_stats = {};
            const auto insert_pair = [&](
                const std::vector<uint8_t>& left_rhs,
                const std::vector<uint8_t>& right_rhs,
                ResidualInsertResult expected) -> bool
            {
                std::vector<uint64_t> left_coeff(words, 0u);
                left_coeff[0] =
                    (UINT64_C(1) << 0u) |
                    (UINT64_C(1) << 6u) |
                    (UINT64_C(1) << 12u);
                std::vector<uint64_t> right_coeff = left_coeff;
                std::vector<uint8_t> left = left_rhs;
                std::vector<uint8_t> right = right_rhs;
                const ResidualInsertResult left_result =
                    InsertPackedBinaryResidualRow(
                        left_coeff, left, R, words, block_bytes,
                        legacy_pivots, legacy_pivot_rhs, legacy_have,
                        legacy_rank, legacy_stats,
                        kBatchedResidualRhsMinBlockBytes);
                const ResidualInsertResult right_result =
                    InsertPackedBinaryResidualRow(
                        right_coeff, right, R, words, block_bytes,
                        fused_pivots, fused_pivot_rhs, fused_have,
                        fused_rank, fused_stats,
                        kBatchedResidualRhsMinBlockBytes);
                return left_result == expected && right_result == expected &&
                    left_coeff == right_coeff && left == right &&
                    legacy_pivots == fused_pivots &&
                    legacy_pivot_rhs == fused_pivot_rhs &&
                    legacy_have == fused_have &&
                    legacy_rank == fused_rank &&
                    legacy_stats.BlockXors == fused_stats.BlockXors &&
                    legacy_stats.BlockMulAdds == fused_stats.BlockMulAdds;
            };
            if (!insert_pair(
                    legacy_rhs, fused_rhs,
                    ResidualInsertResult::Inserted) ||
                !insert_pair(
                    legacy_rhs, fused_rhs,
                    ResidualInsertResult::Dependent))
            {
                return false;
            }
            legacy_rhs[block_bytes / 2u] ^= 1u;
            fused_rhs[block_bytes / 2u] ^= 1u;
            if (!insert_pair(
                    legacy_rhs, fused_rhs,
                    ResidualInsertResult::Inconsistent))
            {
                return false;
            }
        }
    }
    catch (...) {
        return false;
    }
    return true;
}

#undef WH2_TEST_NOINLINE
#undef WH2_TEST_LAMBDA_NOINLINE
#endif

WirehairResult SolvePrecodeSystemWithRuntime(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const PacketRowRuntime& runtime,
    const std::vector<SolvePacket>& packets,
    uint32_t block_bytes,
    std::vector<uint8_t>& intermediate_blocks_out,
    PrecodeSolveStats* stats,
    PrecodeSolveResumeState* resume_state)
{
    return DispatchPrecodeSolve(
        system, config, runtime, packets, block_bytes,
        intermediate_blocks_out, stats, resume_state, true);
}

WirehairResult SolvePrecodeSystemForValidatedSystemWithRuntime(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const PacketRowRuntime& runtime,
    const std::vector<SolvePacket>& packets,
    uint32_t block_bytes,
    std::vector<uint8_t>& intermediate_blocks_out,
    PrecodeSolveStats* stats,
    PrecodeSolveResumeState* resume_state)
{
    return DispatchPrecodeSolve(
        system, config, runtime, packets, block_bytes,
        intermediate_blocks_out, stats, resume_state, false);
}

WirehairResult ResumePrecodeSystem(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    uint32_t block_id,
    const uint8_t* block_data,
    uint32_t block_bytes,
    PrecodeSolveResumeState& state,
    std::vector<uint8_t>& intermediate_blocks_out,
    PrecodeSolveStats* stats,
    bool allow_insert)
{
    const uint32_t K = system.Params.BlockCount;
    const uint64_t P_wide = (uint64_t)system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    if (system.Params.Field == CompletionField::MixedGF256GF16 ||
        !block_data || !state.Active || P_wide > UINT32_MAX ||
        state.SourceCount != K || state.PrecodeCount != (uint32_t)P_wide ||
        state.ColumnCount != K + (uint32_t)P_wide ||
        state.BlockBytes != block_bytes || block_bytes == 0u ||
        state.Config.PeelSeed != config.PeelSeed ||
        state.Config.MixCount != config.MixCount ||
        !state.Runtime.IsValidFor(
            K, (uint32_t)P_wide, config.MixCount) ||
        state.InactiveCount == 0u ||
        state.Rank >= state.InactiveCount ||
        state.ProjectionWords != PackedWordCount(state.InactiveCount) ||
        state.InactiveIndex.size() != state.ColumnCount ||
        state.InactiveColumns.size() != state.InactiveCount ||
        state.Projection.size() !=
            (size_t)state.ColumnCount * state.ProjectionWords ||
        state.Values.size() != (size_t)state.ColumnCount * block_bytes ||
        state.PivotCoefficients.size() !=
            (size_t)state.InactiveCount * state.InactiveCount ||
        state.PivotRhs.size() !=
            (size_t)state.InactiveCount * block_bytes ||
        state.HavePivot.size() != state.InactiveCount ||
        state.CoefficientScratch.size() != state.InactiveCount ||
        state.RhsScratch.size() != block_bytes)
    {
        return Wirehair_InvalidInput;
    }

    try
    {
        typedef std::chrono::steady_clock SolveClock;
        const SolveClock::time_point build_start = SolveClock::now();
        const std::vector<uint32_t> columns =
            GeneratePacketMatrixRowWithRuntime(
                K, (uint32_t)P_wide, block_id, config, state.Runtime);
        if (columns.empty()) {
            return Wirehair_InvalidInput;
        }

        // Duplicate validation is contractually non-mutating.  In particular,
        // do not use the checkpoint's reusable scratch vectors for that path:
        // callers may retry an allocation failure and tests may compare the
        // complete checkpoint byte-for-byte.  The inserting path stays
        // allocation-free after row generation by reusing persistent scratch.
        std::vector<uint8_t> checked_coeff;
        std::vector<uint8_t> checked_rhs;
        if (!allow_insert)
        {
            checked_coeff.resize(state.InactiveCount);
            checked_rhs.resize(block_bytes);
        }
        std::vector<uint8_t>& coeff = allow_insert ?
            state.CoefficientScratch : checked_coeff;
        std::vector<uint8_t>& rhs = allow_insert ?
            state.RhsScratch : checked_rhs;
        std::fill(coeff.begin(), coeff.end(), uint8_t{0});
        std::memcpy(rhs.data(), block_data, block_bytes);
        for (uint32_t column : columns)
        {
            if (column >= state.ColumnCount) {
                return Wirehair_InvalidInput;
            }
            const uint32_t inactive = state.InactiveIndex[column];
            if (inactive != UINT32_MAX) {
                coeff[inactive] ^= 1u;
            }
            else
            {
                const uint64_t* bits = state.Projection.data() +
                    (size_t)column * state.ProjectionWords;
                for (uint32_t word_i = 0;
                     word_i < state.ProjectionWords;
                     ++word_i)
                {
                    uint64_t word = bits[word_i];
                    while (word != 0u)
                    {
                        const uint32_t bit =
                            wirehair::NonzeroLowestBitIndex64(word);
                        const uint32_t index = (word_i << 6) + bit;
                        if (index < state.InactiveCount) {
                            coeff[index] ^= 1u;
                        }
                        word &= word - 1u;
                    }
                }
            }
            gf256_add_mem(
                rhs.data(),
                state.Values.data() + (size_t)column * block_bytes,
                (int)block_bytes);
        }
        const SolveClock::time_point residual_start = SolveClock::now();

        PrecodeSolveStats local_stats = state.Stats;
        PrecodeSolveStats& insertion_stats = allow_insert ?
            state.Stats : local_stats;
        uint32_t checked_rank = state.Rank;
        const ResidualInsertResult insertion = InsertResidualRow(
            coeff,
            rhs,
            state.InactiveCount,
            block_bytes,
            state.PivotCoefficients,
            state.PivotRhs,
            state.HavePivot,
            checked_rank,
            insertion_stats,
            allow_insert,
            kNeverBatchResidualRhs);

        if (!allow_insert)
        {
            return insertion == ResidualInsertResult::Dependent ?
                Wirehair_NeedMore : Wirehair_Error;
        }
        state.Rank = checked_rank;
        ++state.Stats.PacketRows;
        state.Stats.BinaryRowReferences += columns.size();
        ++state.Stats.ResidualRows;
        state.Stats.BuildNanoseconds += (uint64_t)
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                residual_start - build_start).count();
        const SolveClock::time_point residual_end = SolveClock::now();
        state.Stats.ResidualNanoseconds += (uint64_t)
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                residual_end - residual_start).count();
        state.Stats.ResidualRank = state.Rank;

        if (insertion == ResidualInsertResult::Inconsistent ||
            insertion == ResidualInsertResult::Independent)
        {
            if (stats) {
                *stats = state.Stats;
            }
            return Wirehair_Error;
        }
        if (state.Rank < state.InactiveCount)
        {
            if (stats) {
                *stats = state.Stats;
            }
            return Wirehair_NeedMore;
        }

        const SolveClock::time_point backsub_start = SolveClock::now();
        for (uint32_t i = 0; i < state.InactiveCount; ++i)
        {
            if (!state.HavePivot[i]) {
                return Wirehair_Error;
            }
            std::memcpy(
                state.Values.data() +
                    (size_t)state.InactiveColumns[i] * block_bytes,
                state.PivotRhs.data() + (size_t)i * block_bytes,
                block_bytes);
        }
        for (uint32_t column = 0; column < state.ColumnCount; ++column)
        {
            if (state.InactiveIndex[column] != UINT32_MAX) {
                continue;
            }
            uint8_t* value =
                state.Values.data() + (size_t)column * block_bytes;
            const uint64_t* bits = state.Projection.data() +
                (size_t)column * state.ProjectionWords;
            for (uint32_t word_i = 0;
                 word_i < state.ProjectionWords;
                 ++word_i)
            {
                uint64_t word = bits[word_i];
                while (word != 0u)
                {
                    const uint32_t bit =
                        wirehair::NonzeroLowestBitIndex64(word);
                    const uint32_t index = (word_i << 6) + bit;
                    if (index < state.InactiveCount) {
                        gf256_add_mem(
                            value,
                            state.Values.data() +
                                (size_t)state.InactiveColumns[index] *
                                    block_bytes,
                            (int)block_bytes);
                        ++state.Stats.BlockXors;
                    }
                    word &= word - 1u;
                }
            }
        }
        state.Stats.BackSubNanoseconds += (uint64_t)
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                SolveClock::now() - backsub_start).count();
        const PrecodeSolveStats completed_stats = state.Stats;
        intermediate_blocks_out.swap(state.Values);
        state.Clear();
        if (stats) {
            *stats = completed_stats;
        }
        return Wirehair_Success;
    }
    catch (const std::bad_alloc&) {
        if (stats) {
            *stats = state.Stats;
        }
        return Wirehair_OOM;
    }
    catch (const std::length_error&) {
        if (stats) {
            *stats = state.Stats;
        }
        return Wirehair_OOM;
    }
}

WirehairResult SelectSystematicPacketConfig(
    const PrecodeSystem& system,
    const PacketRowConfig& base_config,
    PacketRowConfig& selected_config,
    uint32_t* attempt_out)
{
    const uint32_t K = system.Params.BlockCount;
    const uint64_t P_wide = (uint64_t)system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    if (P_wide > UINT32_MAX ||
        !IsPacketRowDomainValid(
            K, (uint32_t)P_wide, base_config.MixCount) ||
        !ValidatePrecodeSystem(system))
    {
        return Wirehair_InvalidInput;
    }

    try
    {
        const uint8_t zero[2] = {0u, 0u};
        const uint32_t probe_bytes =
            system.Params.Field == CompletionField::MixedGF256GF16 ? 2u : 1u;
        std::vector<SolvePacket> packets(K);
        for (uint32_t block_id = 0; block_id < K; ++block_id)
        {
            packets[block_id].BlockId = block_id;
            packets[block_id].Data = zero;
        }
        for (uint32_t attempt = 0;
             attempt < kMaxPacketSeedAttempts;
             ++attempt)
        {
            const PacketRowConfig candidate =
                PacketConfigForAttempt(base_config, attempt);
            std::vector<uint8_t> intermediate;
            const WirehairResult result = SolvePrecodeSystem(
                system, candidate, packets, probe_bytes, intermediate);
            if (result == Wirehair_Success)
            {
                selected_config = candidate;
                if (attempt_out) {
                    *attempt_out = attempt;
                }
                return Wirehair_Success;
            }
            if (result != Wirehair_NeedMore) {
                return result;
            }
        }
        return Wirehair_BadPeelSeed;
    }
    catch (const std::bad_alloc&) {
        return Wirehair_OOM;
    }
    catch (const std::length_error&) {
        return Wirehair_OOM;
    }
}

WirehairResult SelectSystematicConfiguration(
    const PrecodeParams& base_params,
    const PacketRowConfig& base_config,
    PrecodeSystem& selected_system,
    PacketRowConfig& selected_config,
    uint32_t* attempt_out)
{
    try
    {
        const uint32_t K = base_params.BlockCount;
        const uint64_t P_wide = (uint64_t)base_params.Staircase +
            base_params.DenseRows + base_params.HeavyRows;
        if (P_wide > UINT32_MAX ||
            !IsPacketRowDomainValid(
                K, (uint32_t)P_wide, base_config.MixCount))
        {
            return Wirehair_InvalidInput;
        }
        const uint8_t zero[2] = {0u, 0u};
        const uint32_t probe_bytes =
            base_params.Field == CompletionField::MixedGF256GF16 ? 2u : 1u;
        std::vector<SolvePacket> packets(K);
        for (uint32_t block_id = 0; block_id < K; ++block_id)
        {
            packets[block_id].BlockId = block_id;
            packets[block_id].Data = zero;
        }

        for (uint32_t attempt = 0;
             attempt < kMaxPacketSeedAttempts;
             ++attempt)
        {
            PrecodeSystem system;
            if (!BuildPrecodeSystem(
                    PrecodeParamsForAttempt(base_params, attempt), system))
            {
                return Wirehair_InvalidInput;
            }
            const PacketRowConfig config =
                PacketConfigForAttempt(base_config, attempt);
            std::vector<uint8_t> intermediate;
            const WirehairResult result = SolvePrecodeSystem(
                system, config, packets, probe_bytes, intermediate);
            if (result == Wirehair_Success)
            {
                selected_system = std::move(system);
                selected_config = config;
                if (attempt_out) {
                    *attempt_out = attempt;
                }
                return Wirehair_Success;
            }
            if (result != Wirehair_NeedMore) {
                return result;
            }
        }
        return Wirehair_BadPeelSeed;
    }
    catch (const std::bad_alloc&) {
        return Wirehair_OOM;
    }
    catch (const std::length_error&) {
        return Wirehair_OOM;
    }
}

bool VerifyPrecodeSolution(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const std::vector<SolvePacket>& packets,
    const uint8_t* intermediate_blocks,
    uint32_t block_bytes)
{
    if (!intermediate_blocks ||
        block_bytes == 0u || block_bytes > 0x7fffffffu ||
        (system.Params.Field == CompletionField::MixedGF256GF16 &&
         (block_bytes & 1u) != 0u))
    {
        return false;
    }
    const uint32_t K = system.Params.BlockCount;
    const uint32_t S = system.Params.Staircase;
    const uint32_t D2 = system.Params.DenseRows;
    const uint32_t H = system.Params.HeavyRows;
    const uint64_t P_wide = (uint64_t)S + D2 + H;
    if (P_wide > UINT32_MAX ||
        !IsPacketRowDomainValid(K, (uint32_t)P_wide, config.MixCount) ||
        !ValidatePrecodeSystem(system))
    {
        return false;
    }
    const uint32_t P = (uint32_t)P_wide;
    const uint32_t L = K + P;
    if ((uint64_t)L * block_bytes >
            (uint64_t)std::numeric_limits<size_t>::max())
    {
        return false;
    }
    if (gf256_init() != 0) {
        return false;
    }
    if (system.Params.Field == CompletionField::MixedGF256GF16 &&
        !InitializeGF16())
    {
        return false;
    }
    const bool cached_mixed_coefficients =
        system.Params.Field == CompletionField::MixedGF256GF16;
    const uint32_t mixed_coefficient_period = cached_mixed_coefficients ?
        ActiveMixedCoefficientPeriod() : 0u;
    const bool mixed_residues_rotated = cached_mixed_coefficients &&
        ActiveMixedResiduesRotated();
    const bool independent_extension_residues =
        cached_mixed_coefficients &&
        ActiveMixedIndependentExtensionResidues();
    const uint32_t mixed_subfield_rows = cached_mixed_coefficients ?
        ActiveMixedGF256Rows() : 0u;
    const uint32_t mixed_grouped_gf256_rows = cached_mixed_coefficients ?
        ActiveMixedGroupedGF256Rows() : 0u;
    const uint32_t mixed_first_grouped_gf256_row =
        mixed_grouped_gf256_rows <= mixed_subfield_rows ?
            mixed_subfield_rows - mixed_grouped_gf256_rows : 0u;
    const uint32_t first_heavy_column = L - H;
    const MixedCoefficientRows* mixed_rows = cached_mixed_coefficients ?
        GetMixedCoefficientRows() : nullptr;
    if (cached_mixed_coefficients &&
        (!mixed_rows ||
         mixed_grouped_gf256_rows > mixed_subfield_rows ||
         (mixed_grouped_gf256_rows != 0u &&
          independent_extension_residues)))
    {
        return false;
    }
    std::vector<uint8_t> value(block_bytes, 0u);

    const auto verify_binary = [&](const std::vector<uint32_t>& columns,
                                   const uint8_t* expected) {
        std::fill(value.begin(), value.end(), uint8_t{0});
        for (uint32_t column : columns) {
            gf256_add_mem(
                value.data(),
                intermediate_blocks + (size_t)column * block_bytes,
                (int)block_bytes);
        }
        return expected ?
            std::memcmp(value.data(), expected, block_bytes) == 0 :
            RowIsZero(value.data(), block_bytes);
    };

    for (const std::vector<uint32_t>& row : system.StaircaseRows) {
        if (!verify_binary(row, nullptr)) {
            return false;
        }
    }
    for (const std::vector<uint32_t>& row : system.DenseRowColumns) {
        if (!verify_binary(row, nullptr)) {
            return false;
        }
    }
    for (const SolvePacket& packet : packets)
    {
        if (!packet.Data) {
            return false;
        }
        const std::vector<uint32_t> row =
            GeneratePacketMatrixRow(K, P, packet.BlockId, config);
        if (row.empty() || !verify_binary(row, packet.Data))
        {
            return false;
        }
    }
    for (uint32_t heavy = 0; heavy < H; ++heavy)
    {
        std::fill(value.begin(), value.end(), uint8_t{0});
        uint32_t coefficient_residue = 0u;
        uint32_t coefficient_block_index = 0u;
        uint32_t coefficient_block_column = 0u;
        for (uint32_t column = 0; column < L; ++column)
        {
            const uint32_t active_coefficient_residue = coefficient_residue;
            if (cached_mixed_coefficients)
            {
                if (!mixed_residues_rotated)
                {
                    if (++coefficient_residue == mixed_coefficient_period) {
                        coefficient_residue = 0u;
                    }
                }
                else
                {
                    if (++coefficient_block_column == mixed_coefficient_period)
                    {
                        coefficient_block_column = 0u;
                        coefficient_residue =
                            ActiveMixedResidueBlockShift(
                                ++coefficient_block_index);
                    }
                    else if (++coefficient_residue == mixed_coefficient_period) {
                        coefficient_residue = 0u;
                    }
                }
            }
            if (cached_mixed_coefficients &&
                heavy >= mixed_subfield_rows)
            {
                const uint32_t extension_residue =
                    independent_extension_residues ?
                        ActiveMixedExtensionCoefficientResidue(column) :
                        active_coefficient_residue;
                if (!GF16MulAddMem(
                        value.data(),
                        mixed_rows->Extension[
                            heavy - mixed_subfield_rows][
                                extension_residue],
                        intermediate_blocks +
                            (size_t)column * block_bytes,
                        block_bytes))
                {
                    return false;
                }
            }
            else if (system.Params.Field ==
                    CompletionField::MixedGF256GF16 &&
                heavy >= mixed_subfield_rows)
            {
                if (!GF16MulAddMem(
                        value.data(),
                        MixedGF16Coefficient(
                            heavy - mixed_subfield_rows, column),
                        intermediate_blocks +
                            (size_t)column * block_bytes,
                        block_bytes))
                {
                    return false;
                }
            }
            else
            {
                const uint32_t subfield_residue =
                    cached_mixed_coefficients &&
                    mixed_grouped_gf256_rows != 0u &&
                    heavy >= mixed_first_grouped_gf256_row ?
                        ActiveMixedGroupedGF256CoefficientResidue(
                            column, first_heavy_column) :
                        active_coefficient_residue;
                const uint8_t scale = cached_mixed_coefficients ?
                    mixed_rows->Subfield[heavy][
                        subfield_residue] :
                    HeavyCoefficientForParams(
                        system.Params, heavy, column);
                if (scale == 1u) {
                    gf256_add_mem(
                        value.data(),
                        intermediate_blocks +
                            (size_t)column * block_bytes,
                        (int)block_bytes);
                }
                else {
                    gf256_muladd_mem(
                        value.data(), scale,
                        intermediate_blocks +
                            (size_t)column * block_bytes,
                        (int)block_bytes);
                }
            }
        }
        if (!RowIsZero(value.data(), block_bytes)) {
            return false;
        }
    }
    return true;
}

} // namespace wirehair_v2
