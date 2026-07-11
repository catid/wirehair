#include "WirehairV2Solve.h"

#include "../WirehairTools.h"
#include "../gf256.h"
#include "WirehairV2Plan.h"

#include <algorithm>
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
    const wirehair::RowMixIterator mix(
        params, (uint16_t)precode_count, runtime.PrecodePrime());
    do {
        append(source.GetColumn());
    } while (source.Iterate());
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

bool RowIsZero(const uint8_t* data, uint32_t bytes);

enum class ResidualInsertResult
{
    Dependent,
    Inserted,
    Inconsistent,
    Independent
};

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
        for (uint32_t k = j; k < R; ++k) {
            coeff[k] ^= gf256_mul(scale, pivot[k]);
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
        for (uint32_t k = pivot_column; k < R; ++k) {
            coeff[k] = gf256_mul(coeff[k], inverse);
        }
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
        for (uint32_t k = pivot_column; k < R; ++k) {
            existing_coeff[k] ^= gf256_mul(scale, coeff[k]);
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

PeelResult PeelBinaryRows(
    uint32_t column_count,
    const BinaryEquationArena& rows)
{
    PeelResult out;
    out.SolveRow.assign(column_count, UINT32_MAX);
    out.UsedRows.assign(rows.size(), 0u);

    std::vector<uint32_t> live(rows.size(), 0u);
    std::vector<size_t> column_offsets((size_t)column_count + 1u, 0u);
    std::vector<uint8_t> resolved(column_count, 0u);
    std::vector<uint32_t> queue;
    std::vector<uint32_t> degree_two_refs(column_count, 0u);
    std::priority_queue<ColumnCandidate> degree_two_queue;
    std::priority_queue<ColumnCandidate> reference_queue;

    for (uint32_t r = 0; r < (uint32_t)rows.size(); ++r)
    {
        live[r] = (uint32_t)rows[r].Columns.size();
        if (live[r] == 1u) {
            queue.push_back(r);
        }
        for (uint32_t column : rows[r].Columns) {
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

    const auto adjust_degree_two = [&](uint32_t row, bool add) {
        if (live[row] != 2u || out.UsedRows[row]) {
            return;
        }
        for (uint32_t column : rows[row].Columns)
        {
            if (resolved[column]) {
                continue;
            }
            if (add) {
                ++degree_two_refs[column];
            }
            else if (degree_two_refs[column] > 0u) {
                --degree_two_refs[column];
            }
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
    };
    for (uint32_t r = 0; r < (uint32_t)rows.size(); ++r) {
        adjust_degree_two(r, true);
    }

    const auto resolve = [&](uint32_t column) {
        resolved[column] = 1u;
        for (size_t ref = column_offsets[column];
             ref < column_offsets[(size_t)column + 1u];
             ++ref)
        {
            const uint32_t row = column_rows[ref];
            if (live[row] == 0u) {
                continue;
            }
            adjust_degree_two(row, false);
            --live[row];
            adjust_degree_two(row, true);
            if (live[row] == 1u && !out.UsedRows[row]) {
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
            if (live[row] != 1u || out.UsedRows[row]) {
                continue;
            }
            uint32_t column = UINT32_MAX;
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
    const wirehair::RowMixIterator mix(
        params, (uint16_t)P, runtime.PrecodePrime());
    uint64_t operations = 1u;
    const uint8_t* first_source =
        intermediate_blocks + (size_t)source.GetColumn() * block_bytes;
    if (config.MixCount == kCertifiedPacketMixCount)
    {
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
    else
    {
        // Non-certified experiment configurations retain the generic loop.
        std::memcpy(
            block_out,
            first_source,
            block_bytes);
        while (source.Iterate())
        {
            gf256_add_mem(
                block_out,
                intermediate_blocks +
                    (size_t)source.GetColumn() * block_bytes,
                (int)block_bytes);
            ++operations;
        }
        for (uint32_t i = 0; i < config.MixCount; ++i)
        {
            gf256_add_mem(
                block_out,
                intermediate_blocks +
                    (size_t)(K + mix.Columns[i]) * block_bytes,
                (int)block_bytes);
            ++operations;
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
        !ValidatePrecodeSystem(system))
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
                }
                gf256_add_mem(
                    constant,
                    values.data() + (size_t)other * block_bytes,
                    (int)block_bytes);
                ++st.BlockXors;
            }
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
                for (uint32_t column : rows[r].Columns) {
                    gf256_add_mem(
                        rhs.data(),
                        values.data() + (size_t)column * block_bytes,
                        (int)block_bytes);
                }
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
            std::fill(coeff.begin(), coeff.end(), uint8_t{0});
            std::fill(rhs.begin(), rhs.end(), uint8_t{0});
            if (rows[r].Data) {
                std::memcpy(rhs.data(), rows[r].Data, block_bytes);
            }
            for (uint32_t column : rows[r].Columns)
            {
                const uint32_t index = inactive_index[column];
                if (index != UINT32_MAX) {
                    coeff[index] ^= 1u;
                }
                else
                {
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
                                coeff[projected] ^= 1u;
                            }
                            word &= word - 1u;
                        }
                    }
                }
                gf256_add_mem(
                    rhs.data(),
                    values.data() + (size_t)column * block_bytes,
                    (int)block_bytes);
                ++st.BlockXors;
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

        // Heavy RHS values are bucketed by the coefficient period, avoiding
        // H*L full-block multiplications.  Heavy coefficient vectors are
        // packed eight rows per word so each projection bit is visited once.
        st.BinaryResidualRank = rank;
        const uint32_t window = 256u - H;
        const uint32_t heavy_words = (H + 7u) / 8u;
        if (heavy_words != 0u &&
            (uint64_t)R * heavy_words >
                (uint64_t)std::numeric_limits<size_t>::max() /
                    sizeof(uint64_t))
        {
            return Wirehair_OOM;
        }
        std::vector<uint64_t> projected_heavy(
            (size_t)R * heavy_words, 0u);
        std::vector<uint64_t> packed_heavy(heavy_words, 0u);
        if (H > 0u)
        {
            for (uint32_t column = 0; column < L; ++column)
            {
                std::fill(
                    packed_heavy.begin(), packed_heavy.end(), uint64_t{0});
                for (uint32_t heavy = 0; heavy < H; ++heavy) {
                    packed_heavy[heavy >> 3] |=
                        (uint64_t)HeavyCoefficientForParams(
                            system.Params, heavy, column) <<
                        ((heavy & 7u) * 8u);
                }
                const auto xor_packed = [&](uint32_t index) {
                    uint64_t* destination = projected_heavy.data() +
                        (size_t)index * heavy_words;
                    for (uint32_t w = 0; w < heavy_words; ++w) {
                        destination[w] ^= packed_heavy[w];
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
                for (uint32_t column = residue; column < L; column += window)
                {
                    gf256_add_mem(
                        residue_bucket.data(),
                        values.data() + (size_t)column * block_bytes,
                        (int)block_bytes);
                    ++st.BlockXors;
                }
                for (uint32_t heavy = 0; heavy < H; ++heavy) {
                    heavy_scales[heavy] = HeavyCoefficientForParams(
                        system.Params, heavy, residue);
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
            for (uint32_t other : equation.Columns)
            {
                if (other == column) {
                    continue;
                }
                gf256_add_mem(
                    value,
                    values.data() + (size_t)other * block_bytes,
                    (int)block_bytes);
                ++st.BlockXors;
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
        for (uint32_t column = 0; column < L; ++column)
        {
            const uint8_t scale = HeavyCoefficientForParams(
                system.Params, heavy, column);
            if (scale == 1u) {
                gf256_add_mem(
                    value.data(),
                    intermediate_blocks + (size_t)column * block_bytes,
                    (int)block_bytes);
            }
            else {
                gf256_muladd_mem(
                    value.data(), scale,
                    intermediate_blocks + (size_t)column * block_bytes,
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
