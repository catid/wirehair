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

namespace wirehair_v2 {
namespace {

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

private:
    std::vector<size_t> RowOffsets;
    std::vector<uint32_t> Columns;
    std::vector<const uint8_t*> RowData;
    size_t NextRow = 0u;
};

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
        block_id,
        config.PeelSeed,
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

struct ColumnCandidate
{
    uint32_t Primary = 0;
    uint32_t References = 0;
    uint32_t ReverseColumn = 0;

    bool operator<(const ColumnCandidate& other) const
    {
        if (Primary != other.Primary) {
            return Primary < other.Primary;
        }
        if (References != other.References) {
            return References < other.References;
        }
        return ReverseColumn < other.ReverseColumn;
    }
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

class PairedBlockXorAccumulator
{
public:
    PairedBlockXorAccumulator(uint8_t* destination, uint32_t block_bytes)
        : Destination(destination)
        , BlockBytes(block_bytes)
    {
    }

    void Add(const uint8_t* source)
    {
        if (!PendingSource) {
            PendingSource = source;
            return;
        }
        // Equal sources cancel.  Avoid passing aliased restrict-qualified
        // source pointers to gf256_add2_mem in that unusual case.
        if (PendingSource != source) {
            gf256_add2_mem(
                Destination, PendingSource, source, (int)BlockBytes);
        }
        PendingSource = nullptr;
    }

    void Flush()
    {
        if (!PendingSource) {
            return;
        }
        gf256_add_mem(Destination, PendingSource, (int)BlockBytes);
        PendingSource = nullptr;
    }

private:
    uint8_t* Destination;
    uint32_t BlockBytes;
    const uint8_t* PendingSource = nullptr;
};

bool RowIsZero(const uint8_t* data, uint32_t bytes);

enum class ResidualInsertResult
{
    Dependent,
    Inserted,
    Inconsistent,
    Independent
};

constexpr uint32_t kResidualCoefficientBulkThreshold = 16u;

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
    bool allow_insert)
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

bool ProjectMixedCompletionCoefficientsByResidueBuckets(
    uint32_t column_count,
    uint32_t inactive_count,
    uint32_t projection_words,
    const std::vector<uint32_t>& inactive_index,
    std::vector<uint64_t>& projection,
    const MixedPackedCoefficients& cached_packed,
    std::vector<uint64_t>& projected)
{
    const uint32_t packed_words = ActiveMixedPackedCoefficientWords();
    static_assert(
        kMixedPackedCoefficientWords >= 3u &&
            kMixedPackedCoefficientWords <= 4u,
        "mixed completion packing must fit the unrolled projection");
    const uint32_t expected_projection_words =
        inactive_count / 64u + ((inactive_count & 63u) != 0u ? 1u : 0u);
    const uint32_t coefficient_period = ActiveMixedCoefficientPeriod();
    const uint32_t populated_residues =
        std::min(coefficient_period, column_count);
    const uint64_t projection_elements =
        (uint64_t)column_count * projection_words;
    const uint64_t projected_elements =
        (uint64_t)inactive_count * packed_words;
    if (coefficient_period < kMixedGF256Rows + ActiveMixedGF16Rows() ||
        coefficient_period > kMixedCoefficientPeriod ||
        inactive_index.size() != column_count ||
        projection_words != expected_projection_words ||
        projection_words == 0u ||
        projection_elements > projection.max_size() ||
        projection.size() != (size_t)projection_elements ||
        projected_elements > projected.max_size())
    {
        return false;
    }

    projected.assign((size_t)projected_elements, uint64_t{0});
    const auto xor_projected = [&](uint32_t index,
                                   const uint64_t* coefficients) {
        uint64_t* destination =
            projected.data() + (size_t)index * packed_words;
        destination[0] ^= coefficients[0];
        destination[1] ^= coefficients[1];
        destination[2] ^= coefficients[2];
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
        if (packed_words == 4u) destination[3] ^= coefficients[3];
#endif
    };

    // At or below one complete coefficient period, no two columns share a
    // coefficient vector.  Retain the direct expansion so tiny systems pay no
    // bucket setup cost.
    if (column_count <= coefficient_period)
    {
        for (uint32_t column = 0; column < column_count; ++column)
        {
            const uint64_t* coefficients =
                cached_packed.ByResidue[column];
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

    // Mixed completion coefficients repeat every `coefficient_period`
    // columns.  Transpose the projection by coefficient residue before
    // expanding any bits: each of the first `coefficient_period` projection
    // rows becomes the parity bucket for that residue.  Those rows already
    // contain the first affine vector, and all remaining source rows start at
    // the period, so no source can alias a destination.  The caller has
    // finished all other uses of the projection.
    // This is algebraically identical to dense per-column expansion and
    // requires no additional bucket allocation.
    for (uint32_t initial = 0;
         initial < coefficient_period; ++initial)
    {
        const uint32_t inactive = inactive_index[initial];
        if (inactive == UINT32_MAX) {
            continue;
        }
        if (inactive >= inactive_count) {
            return false;
        }
        uint64_t* bucket = projection.data() +
            (size_t)initial * projection_words;
        std::fill(bucket, bucket + projection_words, uint64_t{0});
        bucket[inactive >> 6] =
            UINT64_C(1) << (inactive & 63u);
    }

    uint32_t residue = 0u;
    for (uint32_t column = coefficient_period;
         column < column_count; ++column)
    {
        uint64_t* bucket = projection.data() +
            (size_t)residue * projection_words;
        const uint32_t inactive = inactive_index[column];
        if (inactive != UINT32_MAX)
        {
            if (inactive >= inactive_count) {
                return false;
            }
            bucket[inactive >> 6] ^=
                UINT64_C(1) << (inactive & 63u);
        }
        else
        {
            const uint64_t* bits =
                projection.data() + (size_t)column * projection_words;
            for (uint32_t word = 0; word < projection_words; ++word) {
                bucket[word] ^= bits[word];
            }
        }
        if (++residue == coefficient_period) {
            residue = 0u;
        }
    }

    for (residue = 0u; residue < populated_residues; ++residue)
    {
        const uint64_t* coefficients = cached_packed.ByResidue[residue];
        const uint64_t* bucket = projection.data() +
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
                    xor_projected(index, coefficients);
                }
                word &= word - 1u;
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
    if (coefficient_period < kMixedGF256Rows + ActiveMixedGF16Rows() ||
        coefficient_period > kMixedCoefficientPeriod)
    {
        return false;
    }
    uint32_t residue = 0u;
    for (uint32_t column = 0; column < column_count; ++column)
    {
        const uint64_t* column_coefficients = cached_packed.ByResidue[residue];
        if (++residue == coefficient_period) {
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

WirehairResult SolveMixedCompletionQuotient(
    const PrecodeSystem& system,
    uint32_t column_count,
    uint32_t inactive_count,
    uint32_t projection_words,
    uint32_t block_bytes,
    const std::vector<uint32_t>& inactive_index,
    const std::vector<uint32_t>& inactive_columns,
    std::vector<uint64_t>& projection,
    const std::vector<uint8_t>& binary_pivot_coeff,
    const std::vector<uint8_t>& binary_pivot_rhs,
    const std::vector<uint8_t>& binary_have_pivot,
    uint32_t binary_rank,
    std::vector<uint8_t>& values,
    PrecodeSolveStats& stats)
{
    const uint32_t extension_rows = ActiveMixedGF16Rows();
    const uint32_t H = kMixedGF256Rows + extension_rows;
    const uint64_t projection_elements =
        (uint64_t)column_count * projection_words;
    if (system.Params.Field != CompletionField::MixedGF256GF16 ||
        system.Params.HeavyRows != H || (block_bytes & 1u) != 0u ||
        binary_rank > inactive_count ||
        inactive_index.size() != column_count ||
        inactive_columns.size() != inactive_count ||
        projection_elements > projection.max_size() ||
        projection.size() != (size_t)projection_elements)
    {
        return Wirehair_InvalidInput;
    }
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

    // Every mixed solve previously built the same complete coefficient period.
    // Share immutable row-major and packed representations across all sizes.
    const MixedCoefficientRows* cached_rows = GetMixedCoefficientRows();
    const MixedPackedCoefficients* cached_packed =
        GetMixedPackedCoefficients();
    if (!cached_rows || !cached_packed) {
        return Wirehair_Error;
    }

    const uint32_t packed_words = ActiveMixedPackedCoefficientWords();
    static_assert(
        kMixedPackedCoefficientWords >= 3u &&
            kMixedPackedCoefficientWords <= 4u,
        "mixed completion packing changed unexpectedly");
    if ((uint64_t)inactive_count * packed_words >
            (uint64_t)std::numeric_limits<size_t>::max() /
                sizeof(uint64_t) ||
        (uint64_t)H * inactive_count >
            (uint64_t)std::numeric_limits<size_t>::max() /
                sizeof(uint16_t))
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
            inactive_index, projection, *cached_packed, projected) ||
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

    std::vector<uint8_t> subfield_coeff(
        (size_t)kMixedGF256Rows * inactive_count, 0u);
    std::vector<uint16_t> extension_coeff(
        (size_t)extension_rows * inactive_count, 0u);
    for (uint32_t index = 0; index < inactive_count; ++index) {
        const uint64_t* source =
            projected.data() + (size_t)index * packed_words;
        for (uint32_t row = 0; row < kMixedGF256Rows; ++row) {
            subfield_coeff[(size_t)row * inactive_count + index] =
                (uint8_t)(source[row >> 2] >> ((row & 3u) * 16u));
        }
        for (uint32_t er = 0; er < extension_rows; ++er) {
            const uint32_t row = kMixedGF256Rows + er;
            extension_coeff[(size_t)er * inactive_count + index] =
                (uint16_t)(source[row >> 2] >> ((row & 3u) * 16u));
        }
    }

    const uint32_t elements = block_bytes / 2u;
    std::vector<uint8_t> subfield_rhs(
        (size_t)kMixedGF256Rows * block_bytes, 0u);
    std::vector<uint8_t> rhs_low((size_t)H * elements, 0u);
    std::vector<uint8_t> rhs_high((size_t)H * elements, 0u);
    std::vector<uint8_t> source_low(elements);
    std::vector<uint8_t> source_high(elements);
    std::vector<uint8_t> residue_bucket(block_bytes, 0u);
    void* subfield_destinations[kMixedGF256Rows];
    void* subfield_coefficient_destinations[kMixedGF256Rows];
    uint8_t subfield_scales[kMixedGF256Rows];
    for (uint32_t row = 0; row < kMixedGF256Rows; ++row) {
        subfield_destinations[row] =
            subfield_rhs.data() + (size_t)row * block_bytes;
        subfield_coefficient_destinations[row] =
            subfield_coeff.data() + (size_t)row * inactive_count;
    }
#if defined(GF256_TRY_AVX2) || defined(GF256_TRY_SSSE3) || \
    defined(GF256_TRY_GFNI) || defined(GF256_TRY_NEON)
    const bool use_bulk_subfield_coefficients = inactive_count >= 64u;
#else
    const bool use_bulk_subfield_coefficients = false;
#endif
    const uint32_t coefficient_period = ActiveMixedCoefficientPeriod();
    const uint32_t populated_residues =
        std::min(coefficient_period, column_count);
    for (uint32_t residue = 0; residue < populated_residues; ++residue)
    {
        std::fill(
            residue_bucket.begin(), residue_bucket.end(), uint8_t{0});
        PairedBlockXorAccumulator bucket_xor(
            residue_bucket.data(), block_bytes);
        for (uint32_t column = residue;
             column < column_count; column += coefficient_period)
        {
            if (inactive_index[column] != UINT32_MAX) continue;
            bucket_xor.Add(
                values.data() + (size_t)column * block_bytes);
            ++stats.BlockXors;
        }
        bucket_xor.Flush();
        for (uint32_t row = 0; row < kMixedGF256Rows; ++row) {
            subfield_scales[row] = cached_rows->Subfield[row][residue];
        }
        AddScaledBlocks(
            subfield_destinations, subfield_scales,
            kMixedGF256Rows, residue_bucket.data(), block_bytes, stats);
        if (!GF16Deinterleave(
                residue_bucket.data(), source_low.data(), source_high.data(),
                block_bytes))
        {
            return Wirehair_Error;
        }
        static_assert(
            kMixedGF16Rows >= 2u,
            "mixed completion pair kernel requires two GF16 rows");
        const uint32_t row0 = kMixedGF256Rows;
        const uint32_t row1 = row0 + 1u;
        if (!GF16MulAddPlanar2(
                rhs_low.data() + (size_t)row0 * elements,
                rhs_high.data() + (size_t)row0 * elements,
                cached_rows->Extension[0][residue],
                rhs_low.data() + (size_t)row1 * elements,
                rhs_high.data() + (size_t)row1 * elements,
                cached_rows->Extension[1][residue],
                source_low.data(), source_high.data(), elements))
        {
            return Wirehair_Error;
        }
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
        for (uint32_t er = 2u; er < extension_rows; ++er)
        {
            const uint32_t row = kMixedGF256Rows + er;
            if (!GF16MulAddPlanar(
                    rhs_low.data() + (size_t)row * elements,
                    rhs_high.data() + (size_t)row * elements,
                    cached_rows->Extension[er][residue],
                    source_low.data(), source_high.data(), elements))
            {
                return Wirehair_Error;
            }
        }
#endif
        stats.BlockMulAdds += extension_rows;
    }

    // Reduce all completion rows through each binary pivot together.  The
    // first ten stay in the GF(256) subfield; the extension rows share a
    // single pivot-RHS deinterleave.
    for (uint32_t pivot = 0; pivot < inactive_count; ++pivot)
    {
        if (!binary_have_pivot[pivot]) continue;
        const uint8_t* relation =
            binary_pivot_coeff.data() + (size_t)pivot * inactive_count;
        for (uint32_t k = 0; k < inactive_count; ++k) {
            if (relation[k] > 1u) return Wirehair_Error;
        }
        for (uint32_t row = 0; row < kMixedGF256Rows; ++row)
        {
            subfield_scales[row] =
                subfield_coeff[(size_t)row * inactive_count + pivot];
        }
        if (use_bulk_subfield_coefficients)
        {
            // The relation is binary, so multiplying it by each pivot scale
            // exactly reproduces the conditional XOR below while sharing the
            // relation scan across all ten coefficient rows where supported.
            gf256_muladd_multi_mem(
                subfield_coefficient_destinations,
                subfield_scales,
                (int)kMixedGF256Rows,
                relation,
                (int)inactive_count);
        }
        else
        {
            for (uint32_t row = 0; row < kMixedGF256Rows; ++row)
            {
                const uint8_t scale = subfield_scales[row];
                if (scale == 0u) continue;
                uint8_t* coefficients =
                    subfield_coeff.data() + (size_t)row * inactive_count;
                for (uint32_t k = 0; k < inactive_count; ++k) {
                    if (relation[k]) {
                        coefficients[k] ^= scale;
                    }
                }
            }
        }
        bool have_subfield_scale = false;
        for (uint32_t row = 0; row < kMixedGF256Rows; ++row) {
            have_subfield_scale =
                have_subfield_scale || subfield_scales[row] != 0u;
        }
        if (have_subfield_scale) {
            AddScaledBlocks(
                subfield_destinations, subfield_scales,
                kMixedGF256Rows,
                binary_pivot_rhs.data() + (size_t)pivot * block_bytes,
                block_bytes, stats);
        }

        uint16_t extension_scales[kMixedGF16RowsMax] = {};
        bool have_extension_scale = false;
        for (uint32_t er = 0; er < extension_rows; ++er)
        {
            const uint16_t scale =
                extension_coeff[(size_t)er * inactive_count + pivot];
            extension_scales[er] = scale;
            have_extension_scale = have_extension_scale || scale != 0u;
            if (scale == 0u) continue;
            for (uint32_t k = 0; k < inactive_count; ++k) {
                if (relation[k]) {
                    extension_coeff[(size_t)er * inactive_count + k] ^=
                        scale;
                }
            }
        }
        if (have_extension_scale)
        {
            if (!GF16Deinterleave(
                    binary_pivot_rhs.data() +
                        (size_t)pivot * block_bytes,
                    source_low.data(), source_high.data(), block_bytes))
            {
                return Wirehair_Error;
            }
            if (!GF16MulAddPlanar2(
                    rhs_low.data() + (size_t)kMixedGF256Rows * elements,
                    rhs_high.data() + (size_t)kMixedGF256Rows * elements,
                    extension_scales[0],
                    rhs_low.data() +
                        (size_t)(kMixedGF256Rows + 1u) * elements,
                    rhs_high.data() +
                        (size_t)(kMixedGF256Rows + 1u) * elements,
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
                const uint32_t row = kMixedGF256Rows + er;
                if (!GF16MulAddPlanar(
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

    for (uint32_t row = 0; row < kMixedGF256Rows; ++row) {
        if (!GF16Deinterleave(
                subfield_rhs.data() + (size_t)row * block_bytes,
                rhs_low.data() + (size_t)row * elements,
                rhs_high.data() + (size_t)row * elements,
                block_bytes))
        {
            return Wirehair_Error;
        }
    }

    if ((uint64_t)H * quotient_columns >
            (uint64_t)std::numeric_limits<size_t>::max() /
                sizeof(uint16_t) ||
        (uint64_t)quotient_columns * quotient_columns >
            (uint64_t)std::numeric_limits<size_t>::max() /
                sizeof(uint16_t))
    {
        return Wirehair_OOM;
    }
    std::vector<uint16_t> quotient_rows(
        (size_t)H * quotient_columns, 0u);
    for (uint32_t row = 0; row < H; ++row) {
        for (uint32_t i = 0; i < quotient_columns; ++i) {
            const uint32_t column = free_columns[i];
            quotient_rows[(size_t)row * quotient_columns + i] =
                row < kMixedGF256Rows ?
                subfield_coeff[(size_t)row * inactive_count + column] :
                extension_coeff[
                    (size_t)(row - kMixedGF256Rows) * inactive_count +
                    column];
        }
    }

    std::vector<uint16_t> pivot_coeff(
        (size_t)quotient_columns * quotient_columns, 0u);
    std::vector<uint8_t> pivot_low(
        (size_t)quotient_columns * elements, 0u);
    std::vector<uint8_t> pivot_high(
        (size_t)quotient_columns * elements, 0u);
    std::vector<uint8_t> have_pivot(quotient_columns, 0u);
    std::vector<uint8_t> scale_scratch(elements);
    uint32_t quotient_rank = 0u;
    for (uint32_t row = 0; row < H; ++row)
    {
        ++stats.ResidualRows;
        uint16_t* coeff = quotient_columns == 0u ? nullptr :
            quotient_rows.data() + (size_t)row * quotient_columns;
        uint8_t* low = rhs_low.data() + (size_t)row * elements;
        uint8_t* high = rhs_high.data() + (size_t)row * elements;
        for (uint32_t column = 0; column < quotient_columns; ++column)
        {
            const uint16_t scale = coeff[column];
            if (scale == 0u || !have_pivot[column]) continue;
            const uint16_t* pivot =
                pivot_coeff.data() + (size_t)column * quotient_columns;
            for (uint32_t k = column; k < quotient_columns; ++k) {
                coeff[k] ^= GF16MultiplyInitialized(scale, pivot[k]);
            }
            if (scale == 1u) {
                gf256_add_mem(
                    low, pivot_low.data() + (size_t)column * elements,
                    (int)elements);
                gf256_add_mem(
                    high, pivot_high.data() + (size_t)column * elements,
                    (int)elements);
                ++stats.BlockXors;
            }
            else if (!GF16MulAddPlanar(
                    low, high, scale,
                    pivot_low.data() + (size_t)column * elements,
                    pivot_high.data() + (size_t)column * elements,
                    elements))
            {
                return Wirehair_Error;
            }
            else {
                ++stats.BlockMulAdds;
            }
        }

        uint32_t pivot_column = quotient_columns;
        for (uint32_t column = 0; column < quotient_columns; ++column) {
            if (coeff[column] != 0u) {
                pivot_column = column;
                break;
            }
        }
        if (pivot_column == quotient_columns)
        {
            if (!RowIsZero(low, elements) || !RowIsZero(high, elements)) {
                return Wirehair_Error;
            }
            continue;
        }
        const uint16_t pivot_value = coeff[pivot_column];
        if (pivot_value != 1u)
        {
            const uint16_t inverse = GF16InverseInitialized(pivot_value);
            if (inverse == 0u) return Wirehair_Error;
            for (uint32_t k = pivot_column;
                 k < quotient_columns; ++k)
            {
                coeff[k] = GF16MultiplyInitialized(coeff[k], inverse);
            }
            if (!GF16ScalePlanar(
                    low, high, inverse, scale_scratch.data(), elements))
            {
                return Wirehair_Error;
            }
            ++stats.BlockMulAdds;
        }
        uint8_t* pending_existing_low = nullptr;
        uint8_t* pending_existing_high = nullptr;
        uint16_t pending_existing_scale = 0u;
        for (uint32_t existing = 0;
             existing < quotient_columns; ++existing)
        {
            if (!have_pivot[existing]) continue;
            uint16_t* existing_coeff = pivot_coeff.data() +
                (size_t)existing * quotient_columns;
            const uint16_t scale = existing_coeff[pivot_column];
            if (scale == 0u) continue;
            for (uint32_t k = pivot_column;
                 k < quotient_columns; ++k)
            {
                existing_coeff[k] ^=
                    GF16MultiplyInitialized(scale, coeff[k]);
            }
            uint8_t* existing_low =
                pivot_low.data() + (size_t)existing * elements;
            uint8_t* existing_high =
                pivot_high.data() + (size_t)existing * elements;
            if (!pending_existing_low)
            {
                pending_existing_low = existing_low;
                pending_existing_high = existing_high;
                pending_existing_scale = scale;
            }
            else
            {
                if (!GF16MulAddPlanar2(
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
            else if (!GF16MulAddPlanar(
                    pending_existing_low, pending_existing_high,
                    pending_existing_scale, low, high, elements))
            {
                return Wirehair_Error;
            }
            else {
                ++stats.BlockMulAdds;
            }
        }
        std::memcpy(
            pivot_coeff.data() +
                (size_t)pivot_column * quotient_columns,
            coeff, (size_t)quotient_columns * sizeof(uint16_t));
        std::memcpy(
            pivot_low.data() + (size_t)pivot_column * elements,
            low, elements);
        std::memcpy(
            pivot_high.data() + (size_t)pivot_column * elements,
            high, elements);
        have_pivot[pivot_column] = 1u;
        ++quotient_rank;
    }
    stats.ResidualRank = binary_rank + quotient_rank;
    if (quotient_rank < quotient_columns) return Wirehair_NeedMore;

    for (uint32_t i = 0; i < quotient_columns; ++i)
    {
        if (!have_pivot[i] || !GF16Interleave(
                pivot_low.data() + (size_t)i * elements,
                pivot_high.data() + (size_t)i * elements,
                values.data() +
                    (size_t)inactive_columns[free_columns[i]] * block_bytes,
                block_bytes))
        {
            return Wirehair_Error;
        }
    }
    for (uint32_t pivot = 0; pivot < inactive_count; ++pivot)
    {
        if (!binary_have_pivot[pivot]) continue;
        uint8_t* value = values.data() +
            (size_t)inactive_columns[pivot] * block_bytes;
        std::memcpy(
            value,
            binary_pivot_rhs.data() + (size_t)pivot * block_bytes,
            block_bytes);
        const uint8_t* relation =
            binary_pivot_coeff.data() + (size_t)pivot * inactive_count;
        for (uint32_t i = 0; i < quotient_columns; ++i) {
            const uint8_t scale = relation[free_columns[i]];
            if (scale > 1u) return Wirehair_Error;
            if (scale == 0u) continue;
            gf256_add_mem(
                value,
                values.data() +
                    (size_t)inactive_columns[free_columns[i]] * block_bytes,
                (int)block_bytes);
            ++stats.BlockXors;
        }
    }
    return Wirehair_Success;
}

template<bool UseLowDegreeXor>
PeelResult PeelBinaryRowsImplementation(
    uint32_t column_count,
    const BinaryEquationArena& rows)
{
    PeelResult out;
    out.SolveRow.assign(column_count, UINT32_MAX);
    out.UsedRows.assign(rows.size(), 0u);

    std::vector<PeelRowState> row_state(rows.size());
    std::vector<size_t> column_offsets((size_t)column_count + 1u, 0u);
    std::vector<uint8_t> resolved(column_count, 0u);
    std::vector<uint32_t> queue;
    std::vector<uint32_t> degree_two_refs(column_count, 0u);
    std::priority_queue<ColumnCandidate> degree_two_queue;
    std::priority_queue<ColumnCandidate> reference_queue;

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

    for (uint32_t column = 0; column < column_count; ++column)
    {
        ColumnCandidate candidate;
        candidate.Primary = (uint32_t)(
            column_offsets[(size_t)column + 1u] - column_offsets[column]);
        candidate.References = candidate.Primary;
        candidate.ReverseColumn = UINT32_MAX - column;
        reference_queue.push(candidate);
    }

    const auto add_degree_two = [&](uint32_t row) {
        PeelRowState& state = row_state[row];
        if (state.Live != 2u || out.UsedRows[row]) {
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
                ColumnCandidate candidate;
                candidate.Primary = degree_two_refs[column];
                candidate.References = (uint32_t)(
                    column_offsets[(size_t)column + 1u] -
                    column_offsets[column]);
                candidate.ReverseColumn = UINT32_MAX - column;
                degree_two_queue.push(candidate);
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
        if (state.Live != 2u || out.UsedRows[row]) {
            return;
        }
        const uint32_t other =
            (uint32_t)state.LowDegreeXor ^ resolved_column;
        CAT_DEBUG_ASSERT(other < column_count && !resolved[other]);
        if (degree_two_refs[other] > 0u) {
            --degree_two_refs[other];
        }
        if (degree_two_refs[other] > 0u) {
            ColumnCandidate candidate;
            candidate.Primary = degree_two_refs[other];
            candidate.References = (uint32_t)(
                column_offsets[(size_t)other + 1u] -
                column_offsets[other]);
            candidate.ReverseColumn = UINT32_MAX - other;
            degree_two_queue.push(candidate);
        }
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
            if (row_state[row].Live == 1u && !out.UsedRows[row]) {
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
            if (row_state[row].Live != 1u || out.UsedRows[row]) {
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
            --remaining;
        }
        if (remaining == 0u) {
            break;
        }

        uint32_t best = UINT32_MAX;
        while (!degree_two_queue.empty())
        {
            const ColumnCandidate candidate = degree_two_queue.top();
            const uint32_t column = UINT32_MAX - candidate.ReverseColumn;
            if (resolved[column] ||
                degree_two_refs[column] != candidate.Primary)
            {
                degree_two_queue.pop();
                continue;
            }
            best = column;
            break;
        }
        while (best == UINT32_MAX && !reference_queue.empty())
        {
            const uint32_t column = UINT32_MAX -
                reference_queue.top().ReverseColumn;
            if (resolved[column]) {
                reference_queue.pop();
                continue;
            }
            best = column;
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
    wirehair::PeelRowIterator source,
    uint32_t source_count,
    uint32_t precode_count,
    const PacketRowConfig& config,
    const PacketRowRuntime& runtime,
    const uint8_t* intermediate_blocks,
    uint32_t block_bytes,
    uint8_t* block_out)
{
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
        block_id, config.PeelSeed, (uint16_t)K, (uint16_t)P);
    wirehair::PeelRowIterator source(
        params, (uint16_t)K, runtime.SourcePrime());
    uint64_t operations = 1u;
    const uint8_t* first_source =
        intermediate_blocks + (size_t)source.GetColumn() * block_bytes;
    // The existing schedules are already optimal until the row contains six
    // total terms.  Above that crossover, pairing the complete tail removes
    // at least one destination read/write pass.
    if (block_bytes <= kPacketTailPairMaxBlockBytes &&
        (uint32_t)params.PeelCount + config.MixCount >=
            kPacketTailPairMinTerms)
    {
        EvaluatePacketTailPaired(
            params, source, K, P, config, runtime, intermediate_blocks,
            block_bytes, block_out);
        if (block_ops_out) {
            *block_ops_out =
                (uint64_t)params.PeelCount + config.MixCount;
        }
        return true;
    }
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
    PrecodeSolveStats st = {};
    const uint32_t K = system.Params.BlockCount;
    const uint32_t S = system.Params.Staircase;
    const uint32_t D2 = system.Params.DenseRows;
    const uint32_t H = system.Params.HeavyRows;
    const uint64_t P_wide = (uint64_t)S + D2 + H;
    if (P_wide > UINT32_MAX ||
        !runtime.IsValidFor(K, (uint32_t)P_wide, config.MixCount) ||
        !ValidatePrecodeSystem(system) ||
        (system.Params.Field == CompletionField::MixedGF256GF16 &&
         (block_bytes & 1u) != 0u))
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
        for (const SolvePacket& packet : packets)
        {
            wirehair::PeelRowParameters params;
            if (!InitializePacketRowParameters(
                    K, P, packet.BlockId, config, runtime, params) ||
                !add_references(
                    (size_t)params.PeelCount + config.MixCount))
            {
                return Wirehair_OOM;
            }
        }

        BinaryEquationArena rows;
        rows.Initialize(row_count, reference_count);
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
                },
                [&rows](uint32_t column) { rows.AppendColumn(column); });
            if (!generated || packet_references == 0u) {
                return Wirehair_InvalidInput;
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
        const uint32_t words = (R + 63u) / 64u;

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
        std::vector<uint8_t> values(value_bytes, 0u);
        std::vector<uint64_t> accumulator(words, 0u);

        // Affine projection of peeled columns onto inactive variables.  The
        // block stored in values[column] is the constant term.
        for (uint32_t column : peel.PeelOrder)
        {
            std::fill(accumulator.begin(), accumulator.end(), uint64_t{0});
            uint8_t* constant =
                values.data() + (size_t)column * block_bytes;
            const BinaryEquationView equation =
                rows[peel.SolveRow[column]];
            if (equation.Data) {
                std::memcpy(constant, equation.Data, block_bytes);
            }
            PairedBlockXorAccumulator constant_xor(constant, block_bytes);
            for (uint32_t other : equation.Columns)
            {
                if (other == column) {
                    continue;
                }
                const uint32_t index = inactive_index[other];
                if (index != UINT32_MAX) {
                    accumulator[index >> 6] ^=
                        UINT64_C(1) << (index & 63u);
                }
                else {
                    for (uint32_t w = 0; w < words; ++w) {
                        accumulator[w] ^=
                            projection[(size_t)other * words + w];
                    }
                    // Inactive value slots are still the zero constant at
                    // this stage.  Only peeled columns can contribute to the
                    // affine RHS; XORing an inactive slot would be a full-
                    // block read/write pass with no algebraic effect.
                    constant_xor.Add(
                        values.data() + (size_t)other * block_bytes);
                    ++st.BlockXors;
                }
            }
            constant_xor.Flush();
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
                PairedBlockXorAccumulator rhs_xor(
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

        // Reduced GF(256) pivot rows, indexed by pivot column.  Gauss-Jordan
        // insertion leaves an identity matrix when rank reaches R.
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
            PairedBlockXorAccumulator rhs_xor(rhs.data(), block_bytes);
            for (uint32_t column : rows[r].Columns)
            {
                const uint32_t index = inactive_index[column];
                if (index != UINT32_MAX) {
                    accumulator[index >> 6] ^=
                        UINT64_C(1) << (index & 63u);
                }
                else
                {
                    const uint64_t* bits =
                        projection.data() + (size_t)column * words;
                    for (uint32_t w = 0; w < words; ++w)
                    {
                        accumulator[w] ^= bits[w];
                    }
                    rhs_xor.Add(
                        values.data() + (size_t)column * block_bytes);
                    ++st.BlockXors;
                }
            }
            rhs_xor.Flush();
            // Accumulate the complete GF(2) row in packed form first.  A bit
            // that appears through several peeled projections is expanded
            // only once after its final parity is known.
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
            if (InsertResidualRow(
                    coeff, rhs, R, block_bytes,
                    pivot_coeff, pivot_rhs, have_pivot,
                    rank, st, true) ==
                ResidualInsertResult::Inconsistent)
            {
                return terminal_error();
            }
        }

        if (system.Params.Field == CompletionField::MixedGF256GF16)
        {
            const WirehairResult mixed_result =
                SolveMixedCompletionQuotient(
                    system, L, R, words, block_bytes,
                    inactive_index, peel.InactiveOrder, projection,
                    pivot_coeff, pivot_rhs, have_pivot, rank,
                    values, st);
            phase_end = SolveClock::now();
            st.ResidualNanoseconds = (uint64_t)
                std::chrono::duration_cast<std::chrono::nanoseconds>(
                    phase_end - phase_start).count();
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
                PairedBlockXorAccumulator bucket_xor(
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
                    rank, st, false);
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
                        true) == ResidualInsertResult::Inconsistent)
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
                        rank, st, true) ==
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
                        rank, st, true) ==
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
        }

        // Dependencies of a peeled column were resolved earlier, so forward
        // chronological substitution reconstructs every remaining value.
        for (uint32_t column : peel.PeelOrder)
        {
            uint8_t* value =
                values.data() + (size_t)column * block_bytes;
            const BinaryEquationView equation =
                rows[peel.SolveRow[column]];
            if (equation.Data) {
                std::memcpy(value, equation.Data, block_bytes);
            }
            else {
                std::memset(value, 0, block_bytes);
            }
            PairedBlockXorAccumulator value_xor(value, block_bytes);
            for (uint32_t other : equation.Columns)
            {
                if (other == column) {
                    continue;
                }
                value_xor.Add(
                    values.data() + (size_t)other * block_bytes);
                ++st.BlockXors;
            }
            value_xor.Flush();
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
        state.ProjectionWords != (state.InactiveCount + 63u) / 64u ||
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
            allow_insert);

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
    const uint32_t mixed_coefficient_period =
        ActiveMixedCoefficientPeriod();
    const MixedCoefficientRows* mixed_rows = cached_mixed_coefficients ?
        GetMixedCoefficientRows() : nullptr;
    if (cached_mixed_coefficients && !mixed_rows) {
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
        uint32_t cached_residue = 0u;
        for (uint32_t column = 0; column < L; ++column)
        {
            const uint32_t coefficient_residue = cached_residue;
            if (cached_mixed_coefficients &&
                ++cached_residue == mixed_coefficient_period)
            {
                cached_residue = 0u;
            }
            if (cached_mixed_coefficients &&
                heavy >= kMixedGF256Rows)
            {
                if (!GF16MulAddMem(
                        value.data(),
                        mixed_rows->Extension[
                            heavy - kMixedGF256Rows][coefficient_residue],
                        intermediate_blocks +
                            (size_t)column * block_bytes,
                        block_bytes))
                {
                    return false;
                }
            }
            else if (system.Params.Field ==
                    CompletionField::MixedGF256GF16 &&
                heavy >= kMixedGF256Rows)
            {
                if (!GF16MulAddMem(
                        value.data(),
                        MixedGF16Coefficient(
                            heavy - kMixedGF256Rows, column),
                        intermediate_blocks +
                            (size_t)column * block_bytes,
                        block_bytes))
                {
                    return false;
                }
            }
            else
            {
                const uint8_t scale = cached_mixed_coefficients ?
                    mixed_rows->Subfield[heavy][coefficient_residue] :
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
