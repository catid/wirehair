#include "WirehairV2Solve.h"

#include "../WirehairTools.h"
#include "../gf256.h"
#include "WirehairV2Plan.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstring>
#include <limits>
#include <new>
#include <queue>
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
thread_local bool TinyDirectSolveEnabled = false;
thread_local int64_t TinyDirectAllocationFailureCountdown = -1;
thread_local bool TinyDirectAllocationThrowsLengthError = false;
#endif

constexpr uint32_t PackedWordCount(uint32_t bit_count)
{
    return bit_count / 64u + ((bit_count & 63u) != 0u ? 1u : 0u);
}

static_assert(
    PackedWordCount(UINT32_MAX) == UINT32_C(67108864),
    "packed word count must not wrap at the uint32 boundary");

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
    params.Initialize(
        PacketRowSeedForBlockId(block_id),
        PacketPeelSeedForBlockId(block_id, config),
        (uint16_t)source_count,
        (uint16_t)precode_count);
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
    prepare((size_t)params.PeelCount + config.MixCount);
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
std::atomic<uint32_t> BinaryPeelOracleUsers(0u);
std::atomic<uint64_t> BinaryPeelOracleComparisons(0u);
thread_local uint64_t TinyDirectSolveAttempts = 0u;
thread_local uint64_t TinyDirectSolveCompletions = 0u;
thread_local uint64_t TinyDirectSolveFallbacks = 0u;
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
// The multi-source RHS crossover is at a 4-KiB payload; the extra kernel
// setup is neutral or slower at MTU sizes.
// CheckedBlockStorage caps valid payloads below this sentinel.
constexpr uint32_t kNeverBatchResidualRhs = UINT32_MAX;

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

constexpr uint32_t kProjectedBackSubMinBlockBytes = 64u;
// Paired whole-solver runs show the fused path loses below these scales even
// though the isolated payload kernel is faster.
constexpr uint32_t kFusedBlockXorInitMinBlockBytes = 1280u;
constexpr uint32_t kFusedBlockXorInitMinBlockCount = 10000u;

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
    if (block_bytes < batched_rhs_min_block_bytes)
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
    have_pivot[pivot_column] = 1u;
    ++rank;
    return ResidualInsertResult::Inserted;
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

    std::vector<PeelRowState> row_state(rows.size());
    std::vector<size_t> column_offsets((size_t)column_count + 1u, 0u);
    std::vector<uint8_t> resolved(column_count, 0u);
    std::vector<uint32_t> queue;
    queue.reserve(rows.size());
    std::vector<uint32_t> degree_two_refs(column_count, 0u);
    class DegreeTwoQueue : public std::priority_queue<uint64_t>
    {
    public:
        void Reserve(size_t count) { this->c.reserve(count); }
    };
    DegreeTwoQueue degree_two_queue;
    degree_two_queue.Reserve(column_count);

    for (uint32_t r = 0; r < (uint32_t)rows.size(); ++r)
    {
        CAT_DEBUG_ASSERT(rows[r].Columns.size() <= UINT16_MAX);
        row_state[r].Live = (uint16_t)rows[r].Columns.size();
        if (row_state[r].Live == 1u) {
            queue.push_back(r);
            row_state[r].LowDegreeXor =
                (uint16_t)*rows[r].Columns.begin();
        }
        for (uint32_t column : rows[r].Columns) {
            CAT_DEBUG_ASSERT(column < column_count);
            ++column_offsets[(size_t)column + 1u];
        }
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
    std::vector<size_t> reference_bucket_offsets(
        (size_t)max_reference_count + 2u, 0u);
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
    std::vector<size_t> reference_bucket_cursor = reference_bucket_offsets;
    std::vector<uint32_t> reference_columns(column_count);
    for (uint32_t column = 0; column < column_count; ++column)
    {
        const uint32_t references = (uint32_t)(
            column_offsets[(size_t)column + 1u] - column_offsets[column]);
        reference_columns[reference_bucket_cursor[references]++] = column;
    }
    reference_bucket_cursor = reference_bucket_offsets;

    // Degree-two priorities compare a changing reference count followed by
    // the immutable total-reference count and reverse column id.  Replace
    // the three-field heap node with one 64-bit key.  The low word is a
    // precomputed rank that preserves the exact original tie order: larger
    // total-reference counts first, then lower column ids.
    std::vector<uint32_t> degree_two_tie_rank(column_count);
    std::vector<uint32_t> degree_two_rank_column(column_count);
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
        uint32_t pair_xor = 0u;
        uint32_t pair_count = 0u;
        for (uint32_t column : rows[row].Columns)
        {
            if (resolved[column]) {
                continue;
            }
            pair_xor ^= column;
            ++pair_count;
            ++degree_two_refs[column];
            if (degree_two_refs[column] > 0u) {
                degree_two_queue.push(degree_two_key(column));
            }
        }
        (void)pair_count;
        CAT_DEBUG_ASSERT(pair_count == 2u && pair_xor <= UINT16_MAX);
        state.LowDegreeXor = (uint16_t)pair_xor;
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
        // needs no live heap key.
        state.LowDegreeXor = (uint16_t)other;
    };
    for (uint32_t r = 0; r < (uint32_t)rows.size(); ++r) {
        add_degree_two(r);
    }

    const auto resolve = [&](uint32_t column) {
        resolved[column] = 1u;
        for (size_t ref = column_offsets[column];
             ref < column_offsets[(size_t)column + 1u];
             ++ref)
        {
            const uint32_t row = column_rows[ref];
            if (row_state[row].Live == 0u) {
                continue;
            }
            remove_degree_two(row, column);
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
        while (!degree_two_queue.empty())
        {
            const uint64_t candidate = degree_two_queue.top();
            const uint32_t column =
                degree_two_rank_column[(uint32_t)candidate];
            if (resolved[column] ||
                degree_two_refs[column] != (uint32_t)(candidate >> 32))
            {
                degree_two_queue.pop();
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

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
bool SetPacketRowSeedMultiplierForTesting(uint32_t multiplier)
{
    if (multiplier == 0u || (multiplier & 1u) == 0u) {
        return false;
    }
    PacketRowSeedMultiplier = multiplier;
    return true;
}

void SetPacketRowSeedAvalancheForTesting(bool enabled)
{
    PacketRowSeedAvalanche = enabled;
}

void SetOddPacketPeelSeedXorForTesting(uint32_t seed_xor)
{
    OddPacketPeelSeedXor = seed_xor;
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

void SetTinyDirectSolveForTesting(bool enabled)
{
    TinyDirectSolveEnabled = enabled;
}

void SetTinyDirectSolveAllocationFailureForTesting(
    int64_t countdown,
    bool length_error)
{
    TinyDirectAllocationFailureCountdown = countdown;
    TinyDirectAllocationThrowsLengthError = length_error;
}

void ResetTinyDirectSolveCountersForTesting()
{
    TinyDirectSolveAttempts = 0u;
    TinyDirectSolveCompletions = 0u;
    TinyDirectSolveFallbacks = 0u;
}

uint64_t TinyDirectSolveAttemptsForTesting()
{
    return TinyDirectSolveAttempts;
}

uint64_t TinyDirectSolveCompletionsForTesting()
{
    return TinyDirectSolveCompletions;
}

uint64_t TinyDirectSolveFallbacksForTesting()
{
    return TinyDirectSolveFallbacks;
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
        [&row](size_t count) { row.reserve(count); },
        [&row](uint32_t column) { row.push_back(column); });
    if (!generated) {
        row.clear();
        return row;
    }
    return row;
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
    params.Initialize(
        PacketRowSeedForBlockId(block_id),
        PacketPeelSeedForBlockId(block_id, config),
        (uint16_t)K,
        (uint16_t)P);
    uint64_t operations = 1u;
    // The existing schedules are already optimal until the row contains six
    // total terms.  Above that crossover, pairing the complete tail removes
    // at least one destination read/write pass.
    const uint32_t packet_terms =
        (uint32_t)params.PeelCount + config.MixCount;
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

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
namespace {

constexpr uint32_t kTinyDirectMaxSourceBlocks = 100u;
constexpr uint32_t kTinyDirectMaxColumns = 160u;
constexpr uint32_t kTinyDirectMaxRows = 160u;
constexpr uint32_t kTinyDirectMaxBlockBytes = 4096u;

enum class TinyDirectDisposition
{
    Ineligible,
    Fallback,
    Complete
};

struct TinyDirectOutcome
{
    TinyDirectOutcome(
        TinyDirectDisposition disposition = TinyDirectDisposition::Ineligible,
        WirehairResult result = Wirehair_Error)
        : Disposition(disposition)
        , Result(result)
    {
    }

    TinyDirectDisposition Disposition;
    WirehairResult Result;
};

void TinyDirectAllocationCheckpoint()
{
    if (TinyDirectAllocationFailureCountdown < 0) {
        return;
    }
    if (TinyDirectAllocationFailureCountdown > 0)
    {
        --TinyDirectAllocationFailureCountdown;
        return;
    }
    TinyDirectAllocationFailureCountdown = -1;
    const bool throw_length_error = TinyDirectAllocationThrowsLengthError;
    TinyDirectAllocationThrowsLengthError = false;
    if (throw_length_error) {
        throw std::length_error("injected tiny-direct allocation failure");
    }
    throw std::bad_alloc();
}

bool IsTinyDirectDenseAnchor(
    DenseAnchorLayout layout,
    uint32_t dense_row)
{
    return dense_row == 0u ||
        (layout == DenseAnchorLayout::Two07 && dense_row == 7u);
}

#if defined(_MSC_VER)
#define WH2_TINY_DIRECT_NOINLINE __declspec(noinline)
#elif defined(__GNUC__) || defined(__clang__)
#define WH2_TINY_DIRECT_NOINLINE __attribute__((noinline))
#else
#define WH2_TINY_DIRECT_NOINLINE
#endif

WH2_TINY_DIRECT_NOINLINE TinyDirectOutcome TryTinyDirectSolve(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const PacketRowRuntime& runtime,
    const std::vector<SolvePacket>& packets,
    uint32_t block_bytes,
    size_t value_bytes,
    std::vector<uint8_t>& output,
    PrecodeSolveStats& st)
{
    const PrecodeParams& params = system.Params;
    const uint32_t K = params.BlockCount;
    const uint32_t S = params.Staircase;
    const uint32_t D2 = params.DenseRows;
    const uint32_t H = params.HeavyRows;
    const bool supported_layout =
        params.DenseAnchors == DenseAnchorLayout::Disabled ||
        params.DenseAnchors == DenseAnchorLayout::Two07;
    if (!TinyDirectSolveEnabled || K < 2u ||
        K > kTinyDirectMaxSourceBlocks ||
        S != wirehair::GetDenseCount(K) || D2 != 12u || H != 12u ||
        params.SourceHits != 2u || params.DenseIdentityCorner ||
        params.HeavyFamily != HeavyCoefficientFamily::PeriodicCauchy ||
        !supported_layout || config.MixCount != kCertifiedPacketMixCount ||
        block_bytes > kTinyDirectMaxBlockBytes)
    {
        return TinyDirectOutcome();
    }

    const uint64_t P_wide = (uint64_t)S + D2 + H;
    const uint64_t L_wide = (uint64_t)K + P_wide;
    if (packets.size() > (size_t)K + 4u) {
        return TinyDirectOutcome();
    }
    const uint64_t row_count_wide = P_wide + packets.size();
    if (L_wide == 0u || L_wide > kTinyDirectMaxColumns ||
        row_count_wide == 0u || row_count_wide > kTinyDirectMaxRows ||
        row_count_wide > UINT16_MAX ||
        !runtime.IsValidFor(K, (uint32_t)P_wide, config.MixCount))
    {
        return TinyDirectOutcome();
    }
    const uint32_t P = (uint32_t)P_wide;
    const uint32_t L = (uint32_t)L_wide;
    const uint32_t row_count = (uint32_t)row_count_wide;
    const uint64_t stride_wide = L_wide + block_bytes;
    if (stride_wide > (uint64_t)std::numeric_limits<size_t>::max() ||
        row_count_wide >
            (uint64_t)std::numeric_limits<size_t>::max() / stride_wide)
    {
        return TinyDirectOutcome();
    }
    const size_t stride = (size_t)stride_wide;
    const size_t matrix_bytes = (size_t)(row_count_wide * stride_wide);

    ++TinyDirectSolveAttempts;
    try
    {
        typedef std::chrono::steady_clock SolveClock;
        const SolveClock::time_point build_start = SolveClock::now();

        // Keep coefficients and RHS bytes for each equation in one allocation.
        // The fixed row permutation changes pivot order without either another
        // allocation or copying a payload-sized augmented row.
        TinyDirectAllocationCheckpoint();
        std::vector<uint8_t> matrix(matrix_bytes, 0u);
        std::array<uint16_t, kTinyDirectMaxRows> logical_rows{};

        size_t physical_row = 0u;
        const auto add_binary_row = [&](
            const std::vector<uint32_t>& columns,
            const uint8_t* rhs) -> bool {
            if (physical_row >= row_count) {
                return false;
            }
            if (columns.size() > L) {
                return false;
            }
            uint8_t* const row = matrix.data() + physical_row * stride;
            uint32_t previous = 0u;
            bool have_previous = false;
            for (uint32_t column : columns)
            {
                // Recheck even on the trusted-system entry point: this test
                // optimization must not turn a malformed cached system into an
                // out-of-bounds write.
                if (column >= L) {
                    return false;
                }
                if (have_previous && column <= previous) {
                    return false;
                }
                row[column] = 1u;
                previous = column;
                have_previous = true;
            }
            if (rhs) {
                std::memcpy(row + L, rhs, block_bytes);
            }
            st.BinaryRowReferences += columns.size();
            ++physical_row;
            return true;
        };

        if (system.StaircaseRows.size() != S ||
            system.DenseBasisRowColumns.size() != D2)
        {
            return {TinyDirectDisposition::Complete, Wirehair_InvalidInput};
        }
        for (const std::vector<uint32_t>& columns : system.StaircaseRows) {
            if (!add_binary_row(columns, nullptr)) {
                return {TinyDirectDisposition::Complete, Wirehair_InvalidInput};
            }
        }
        for (const std::vector<uint32_t>& columns :
                system.DenseBasisRowColumns) {
            if (!add_binary_row(columns, nullptr)) {
                return {TinyDirectDisposition::Complete, Wirehair_InvalidInput};
            }
        }
        for (const SolvePacket& packet : packets)
        {
            if (physical_row >= row_count) {
                return {TinyDirectDisposition::Complete, Wirehair_Error};
            }
            uint8_t* const row = matrix.data() + physical_row * stride;
            size_t references = 0u;
            bool valid_column = true;
            const bool generated = ForEachPacketMatrixColumn(
                K,
                P,
                packet.BlockId,
                config,
                runtime,
                [&references](size_t count) { references = count; },
                [&](uint32_t column) {
                    if (column >= L) {
                        valid_column = false;
                    }
                    else {
                        row[column] ^= 1u;
                    }
                });
            if (!generated || !valid_column || references == 0u) {
                return {TinyDirectDisposition::Complete, Wirehair_InvalidInput};
            }
            std::memcpy(row + L, packet.Data, block_bytes);
            if (references >
                std::numeric_limits<uint64_t>::max() -
                    st.BinaryRowReferences)
            {
                return {TinyDirectDisposition::Complete, Wirehair_Error};
            }
            st.BinaryRowReferences += references;
            ++physical_row;
        }
        for (uint32_t heavy = 0u; heavy < H; ++heavy)
        {
            if (physical_row >= row_count) {
                return {TinyDirectDisposition::Complete, Wirehair_Error};
            }
            uint8_t* const row = matrix.data() + physical_row * stride;
            for (uint32_t column = 0u; column < L; ++column) {
                row[column] = HeavyCoefficientForParams(
                    params, heavy, column);
            }
            ++physical_row;
        }
        if (physical_row != row_count) {
            return {TinyDirectDisposition::Complete, Wirehair_Error};
        }

        const uint32_t staircase_base = 0u;
        const uint32_t dense_base = S;
        const uint32_t packet_base = S + D2;
        const uint32_t heavy_base = packet_base + (uint32_t)packets.size();
        uint32_t logical = 0u;
        const auto append_logical = [&](uint32_t row) {
            logical_rows[logical++] = (uint16_t)row;
        };
        for (uint32_t row = 0u; row < S; ++row) {
            append_logical(staircase_base + row);
        }
        // Consume the sparse equation-preserving deltas before packet rows,
        // and leave each half-dense anchor until after the sparse rows.
        for (uint32_t row = 0u; row < D2; ++row) {
            if (!IsTinyDirectDenseAnchor(params.DenseAnchors, row)) {
                append_logical(dense_base + row);
            }
        }
        for (uint32_t row = 0u; row < packets.size(); ++row) {
            append_logical(packet_base + row);
        }
        for (uint32_t row = 0u; row < D2; ++row) {
            if (IsTinyDirectDenseAnchor(params.DenseAnchors, row)) {
                append_logical(dense_base + row);
            }
        }
        for (uint32_t row = 0u; row < H; ++row) {
            append_logical(heavy_base + row);
        }
        if (logical != row_count) {
            return {TinyDirectDisposition::Complete, Wirehair_Error};
        }

        st.PacketRows = (uint32_t)packets.size();
        st.InactivatedColumns = L;
        st.BinaryRowStorageBytes = matrix.capacity() * sizeof(uint8_t);
        st.BinaryRowStorageAllocations = 1u;
        st.ResidualRows = row_count;
        const SolveClock::time_point build_end = SolveClock::now();
        st.BuildNanoseconds = (uint64_t)
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                build_end - build_start).count();

        std::array<uint16_t, kTinyDirectMaxColumns> pivot_columns{};
        uint32_t rank = 0u;
        for (uint32_t column = 0u; column < L && rank < row_count; ++column)
        {
            uint32_t pivot_logical = rank;
            while (pivot_logical < row_count)
            {
                const uint8_t* const candidate = matrix.data() +
                    (size_t)logical_rows[pivot_logical] * stride;
                if (candidate[column] != 0u) {
                    break;
                }
                ++pivot_logical;
            }
            if (pivot_logical == row_count) {
                continue;
            }
            std::swap(logical_rows[rank], logical_rows[pivot_logical]);
            uint8_t* const pivot = matrix.data() +
                (size_t)logical_rows[rank] * stride;
            const uint8_t pivot_value = pivot[column];
            if (pivot_value != 1u)
            {
                gf256_div_mem(
                    pivot + column,
                    pivot + column,
                    pivot_value,
                    (int)(stride - column));
                ++st.BlockMulAdds;
            }
            for (uint32_t target_logical = rank + 1u;
                 target_logical < row_count;
                 ++target_logical)
            {
                uint8_t* const target = matrix.data() +
                    (size_t)logical_rows[target_logical] * stride;
                const uint8_t scale = target[column];
                if (scale == 0u) {
                    continue;
                }
                if (scale == 1u) {
                    gf256_add_mem(
                        target + column,
                        pivot + column,
                        (int)(stride - column));
                    ++st.BlockXors;
                }
                else
                {
                    gf256_muladd_mem(
                        target + column,
                        scale,
                        pivot + column,
                        (int)(stride - column));
                    ++st.BlockMulAdds;
                }
            }
            pivot_columns[rank] = (uint16_t)column;
            ++rank;
        }
        st.ResidualRank = rank;

        // Every non-pivot equation is reduced by every chosen pivot.  Inspect
        // all of them, including surplus rows after full rank, before deciding
        // between Success and NeedMore.
        bool inconsistent = false;
        for (uint32_t logical_row = rank;
             logical_row < row_count;
             ++logical_row)
        {
            const uint8_t* const row = matrix.data() +
                (size_t)logical_rows[logical_row] * stride;
            if (!RowIsZero(row, L) ||
                !RowIsZero(row + L, block_bytes))
            {
                inconsistent = true;
                break;
            }
        }
        const SolveClock::time_point residual_end = SolveClock::now();
        st.ResidualNanoseconds = (uint64_t)
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                residual_end - build_end).count();
        if (inconsistent) {
            return {TinyDirectDisposition::Complete, Wirehair_Error};
        }
        if (rank < L) {
            return {TinyDirectDisposition::Complete, Wirehair_NeedMore};
        }

        TinyDirectAllocationCheckpoint();
        std::vector<uint8_t> solved(value_bytes, 0u);
        for (uint32_t reverse = rank; reverse-- > 0u;)
        {
            const uint8_t* const row = matrix.data() +
                (size_t)logical_rows[reverse] * stride;
            const uint32_t pivot_column = pivot_columns[reverse];
            uint8_t* const value = solved.data() +
                (size_t)pivot_column * block_bytes;
            std::memcpy(value, row + L, block_bytes);
            for (uint32_t column = pivot_column + 1u;
                 column < L;
                 ++column)
            {
                AddScaledBlock(
                    value,
                    row[column],
                    solved.data() + (size_t)column * block_bytes,
                    block_bytes,
                    st);
            }
        }
        const SolveClock::time_point backsub_end = SolveClock::now();
        st.BackSubNanoseconds = (uint64_t)
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                backsub_end - residual_end).count();
        output.swap(solved);
        return {TinyDirectDisposition::Complete, Wirehair_Success};
    }
    catch (const std::bad_alloc&) {
        return {TinyDirectDisposition::Fallback, Wirehair_OOM};
    }
    catch (const std::length_error&) {
        return {TinyDirectDisposition::Fallback, Wirehair_OOM};
    }
}

#undef WH2_TINY_DIRECT_NOINLINE

} // namespace
#endif

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
        (validate_system && !ValidatePrecodeSystem(system)))
    {
        return Wirehair_InvalidInput;
    }
    const uint32_t P = (uint32_t)P_wide;
    const uint32_t L = K + P;
    size_t value_bytes = 0u;
    if (!CheckedBlockStorage(L, block_bytes, value_bytes))
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
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    if (TinyDirectSolveEnabled)
    {
        std::vector<uint8_t> direct_output;
        PrecodeSolveStats direct_stats = {};
        const TinyDirectOutcome direct = TryTinyDirectSolve(
            system,
            config,
            runtime,
            packets,
            block_bytes,
            value_bytes,
            direct_output,
            direct_stats);
        if (direct.Disposition == TinyDirectDisposition::Fallback ||
            (direct.Disposition == TinyDirectDisposition::Complete &&
             direct.Result == Wirehair_NeedMore && resume_state))
        {
            ++TinyDirectSolveFallbacks;
        }
        else if (direct.Disposition == TinyDirectDisposition::Complete)
        {
            ++TinyDirectSolveCompletions;
            if (direct.Result == Wirehair_Success) {
                intermediate_blocks_out.swap(direct_output);
            }
            if (stats && direct.Result != Wirehair_InvalidInput) {
                *stats = direct_stats;
            }
            return direct.Result;
        }
    }
#endif
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
        const std::vector<std::vector<uint32_t>>& dense_basis_rows =
            system.DenseBasisRowColumns;
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
        for (const std::vector<uint32_t>& columns : dense_basis_rows) {
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
        for (const std::vector<uint32_t>& columns : dense_basis_rows)
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
        values.reserve(value_bytes);
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
        std::vector<uint64_t> accumulator(words, 0u);
        const bool enable_fused_block_initialization =
            K >= kFusedBlockXorInitMinBlockCount &&
            block_bytes >= kFusedBlockXorInitMinBlockBytes;

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
                solve_column_offset = AccumulatePeeledProjectionConstant(
                    column, equation, inactive_index, words, projection,
                    values, block_bytes, accumulator, constant_xor, st);
                constant_xor.Flush();
            }
            else
            {
                if (equation.Data) {
                    std::memcpy(constant, equation.Data, block_bytes);
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
            for (uint32_t r = 0; r < (uint32_t)rows.size(); ++r)
            {
                if (peel.UsedRows[r]) {
                    continue;
                }
                std::vector<uint8_t> rhs(block_bytes, 0u);
                if (rows[r].Data) {
                    std::memcpy(rhs.data(), rows[r].Data, block_bytes);
                }
                BatchedBlockXorAccumulator rhs_xor(
                    rhs.data(), block_bytes);
                for (uint32_t column : rows[r].Columns) {
                    rhs_xor.Add(
                        values.data() + (size_t)column * block_bytes);
                }
                rhs_xor.Flush();
                if (!RowIsZero(rhs.data(), block_bytes)) {
                    return terminal_error();
                }
            }
            if (!VerifyPrecodeSolution(
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

        if ((uint64_t)R * R >
            (uint64_t)std::numeric_limits<size_t>::max())
        {
            return Wirehair_OOM;
        }
        std::vector<uint8_t> pivot_coeff((size_t)R * R, 0u);
        size_t residual_value_bytes = 0u;
        if (!CheckedBlockStorage(R, block_bytes, residual_value_bytes)) {
            return Wirehair_OOM;
        }
        std::vector<uint8_t> pivot_rhs(residual_value_bytes, 0u);
        std::vector<uint8_t> have_pivot(R, 0u);
        std::vector<uint8_t> coeff(R, 0u);
        std::vector<uint8_t> rhs(block_bytes, 0u);
        uint32_t rank = 0u;
        const uint32_t batched_rhs_min_block_bytes =
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
            std::fill(rhs.begin(), rhs.end(), uint8_t{0});
            if (rows[r].Data) {
                std::memcpy(rhs.data(), rows[r].Data, block_bytes);
            }
            BatchedBlockXorAccumulator rhs_xor(rhs.data(), block_bytes);
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
                }
            }
            projection_xor.Flush();
            rhs_xor.Flush();
            // GF(256) completion consumes one byte per projected
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
            const ResidualInsertResult insertion = InsertResidualRow(
                coeff, rhs, R, block_bytes,
                pivot_coeff, pivot_rhs, have_pivot,
                rank, st, true, batched_rhs_min_block_bytes);
            if (insertion == ResidualInsertResult::Inconsistent)
            {
                return terminal_error();
            }
        }

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

        std::vector<uint8_t> heavy_rhs((size_t)H * block_bytes, 0u);
        void* heavy_destinations[128];
        uint8_t heavy_scales[128];
        for (uint32_t heavy = 0; heavy < H; ++heavy) {
            heavy_destinations[heavy] =
                heavy_rhs.data() + (size_t)heavy * block_bytes;
        }
        if (system.Params.HeavyFamily ==
            HeavyCoefficientFamily::PeriodicCauchy)
        {
            std::vector<uint8_t> residue_bucket(block_bytes, 0u);
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
        const bool use_binary_quotient =
            block_bytes >= kBinaryQuotientMinBlockBytes;
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
            !CheckedBlockStorage(
                quotient_columns, block_bytes, quotient_value_bytes))
        {
            return Wirehair_OOM;
        }
        std::vector<uint8_t> quotient_pivot_rhs(
            quotient_value_bytes, 0u);
        std::vector<uint8_t> quotient_have_pivot(
            quotient_columns, 0u);
        std::vector<uint8_t> quotient_coeff(quotient_columns, 0u);
        std::vector<uint8_t> quotient_rhs(
            use_binary_quotient ? block_bytes : 0u, 0u);
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
            if (block_bytes >= kProjectedBackSubMinBlockBytes &&
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
                }
                else {
                    std::memset(value, 0, block_bytes);
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

namespace {

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


} // namespace


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
    return SolvePrecodeSystemImpl(
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
    return SolvePrecodeSystemImpl(
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
    if (!block_data || !state.Active || P_wide > UINT32_MAX ||
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
        const uint8_t zero = 0u;
        std::vector<SolvePacket> packets(K);
        for (uint32_t block_id = 0; block_id < K; ++block_id)
        {
            packets[block_id].BlockId = block_id;
            packets[block_id].Data = &zero;
        }
        for (uint32_t attempt = 0;
             attempt < kMaxPacketSeedAttempts;
             ++attempt)
        {
            const PacketRowConfig candidate =
                PacketConfigForAttempt(base_config, attempt);
            std::vector<uint8_t> intermediate;
            const WirehairResult result = SolvePrecodeSystem(
                system, candidate, packets, 1u, intermediate);
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
        const uint8_t zero = 0u;
        std::vector<SolvePacket> packets(K);
        for (uint32_t block_id = 0; block_id < K; ++block_id)
        {
            packets[block_id].BlockId = block_id;
            packets[block_id].Data = &zero;
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
                system, config, packets, 1u, intermediate);
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

WirehairResult ClassifyExactSystematicConstructionFailure(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const PacketRowRuntime& runtime,
    WirehairResult solve_result)
{
    if (solve_result == Wirehair_NeedMore) {
        return Wirehair_BadPeelSeed;
    }
    if (solve_result != Wirehair_Error) {
        return solve_result;
    }

    try
    {
        const uint32_t K = system.Params.BlockCount;
        const uint8_t zero = 0u;
        std::vector<SolvePacket> packets(K);
        for (uint32_t block_id = 0u; block_id < K; ++block_id)
        {
            packets[block_id].BlockId = block_id;
            packets[block_id].Data = &zero;
        }
        std::vector<uint8_t> intermediate;
        const WirehairResult probe =
            SolvePrecodeSystemForValidatedSystemWithRuntime(
                system,
                config,
                runtime,
                packets,
                1u,
                intermediate);
        if (probe == Wirehair_NeedMore) {
            return Wirehair_BadPeelSeed;
        }
        // A full-rank zero-RHS probe proves that the original Error was not a
        // rank failure.  Preserve it rather than hiding an internal defect.
        return probe == Wirehair_Success ? Wirehair_Error : probe;
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
        block_bytes == 0u || block_bytes > 0x7fffffffu)
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
    for (const std::vector<uint32_t>& row :
            system.DenseBasisRowColumns) {
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
        for (uint32_t column = 0; column < L; ++column)
        {
            const uint8_t scale = HeavyCoefficientForParams(
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
        if (!RowIsZero(value.data(), block_bytes)) {
            return false;
        }
    }
    return true;
}

} // namespace wirehair_v2
