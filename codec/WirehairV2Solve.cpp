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

struct BinaryEquation
{
    std::vector<uint32_t> Columns;
    const uint8_t* Data = nullptr;
};

struct PeelResult
{
    std::vector<uint32_t> SolveRow;
    std::vector<uint32_t> PeelOrder;
    std::vector<uint32_t> InactiveOrder;
    std::vector<uint8_t> UsedRows;
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

PeelResult PeelBinaryRows(
    uint32_t column_count,
    const std::vector<BinaryEquation>& rows)
{
    PeelResult out;
    out.SolveRow.assign(column_count, UINT32_MAX);
    out.UsedRows.assign(rows.size(), 0u);

    std::vector<uint32_t> live(rows.size(), 0u);
    std::vector<std::vector<uint32_t> > column_rows(column_count);
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
            column_rows[column].push_back(r);
        }
    }

    for (uint32_t column = 0; column < column_count; ++column) {
        ColumnCandidate candidate;
        candidate.Primary = (uint32_t)column_rows[column].size();
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
                candidate.References =
                    (uint32_t)column_rows[column].size();
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
        for (uint32_t row : column_rows[column])
        {
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

uint32_t LowestSetBitIndex(uint64_t word)
{
    uint32_t index = 0u;
    while ((word & 1u) == 0u)
    {
        word >>= 1;
        ++index;
    }
    return index;
}

} // namespace

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

std::vector<uint32_t> GeneratePacketMatrixRow(
    uint32_t source_count,
    uint32_t precode_count,
    uint32_t block_id,
    const PacketRowConfig& config)
{
    std::vector<uint32_t> row;
    if (!IsPacketRowDomainValid(
            source_count, precode_count, config.MixCount))
    {
        return row;
    }

    wirehair::PeelRowParameters params;
    params.Initialize(
        block_id,
        config.PeelSeed,
        (uint16_t)source_count,
        (uint16_t)precode_count);
    const uint16_t source_prime =
        wirehair::NextPrime16((uint16_t)source_count);
    const uint16_t precode_prime =
        wirehair::NextPrime16((uint16_t)precode_count);
    wirehair::PeelRowIterator source(
        params, (uint16_t)source_count, source_prime);
    const wirehair::RowMixIterator mix(
        params, (uint16_t)precode_count, precode_prime);

    row.reserve((size_t)params.PeelCount + config.MixCount);
    do {
        row.push_back(source.GetColumn());
    } while (source.Iterate());
    for (uint32_t i = 0; i < config.MixCount; ++i) {
        row.push_back(source_count + mix.Columns[i]);
    }
    return row;
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

static bool EvaluatePacketBlockImpl(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
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
        !IsPacketRowDomainValid(K, (uint32_t)P_wide, config.MixCount) ||
        ((uint64_t)K + P_wide) * block_bytes >
            (uint64_t)std::numeric_limits<size_t>::max())
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
        params, (uint16_t)K, wirehair::NextPrime16((uint16_t)K));
    const wirehair::RowMixIterator mix(
        params, (uint16_t)P, wirehair::NextPrime16((uint16_t)P));
    uint64_t operations = 1u;
    std::memcpy(
        block_out,
        intermediate_blocks + (size_t)source.GetColumn() * block_bytes,
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
    return EvaluatePacketBlockImpl(
        system, config, intermediate_blocks, block_bytes, block_id,
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
    return EvaluatePacketBlockImpl(
        system, config, intermediate_blocks, block_bytes, block_id,
        block_out, block_ops_out, false);
}

WirehairResult SolvePrecodeSystem(
    const PrecodeSystem& system,
    const PacketRowConfig& config,
    const std::vector<SolvePacket>& packets,
    uint32_t block_bytes,
    std::vector<uint8_t>& intermediate_blocks_out,
    PrecodeSolveStats* stats)
{
    PrecodeSolveStats st = {};
    const uint32_t K = system.Params.BlockCount;
    const uint32_t S = system.Params.Staircase;
    const uint32_t D2 = system.Params.DenseRows;
    const uint32_t H = system.Params.HeavyRows;
    const uint64_t P_wide = (uint64_t)S + D2 + H;
    if (P_wide > UINT32_MAX ||
        !IsPacketRowDomainValid(K, (uint32_t)P_wide, config.MixCount) ||
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

    try
    {
        typedef std::chrono::steady_clock SolveClock;
        SolveClock::time_point phase_start = SolveClock::now();
        std::vector<BinaryEquation> rows;
        rows.reserve((size_t)S + D2 + packets.size());
        for (const std::vector<uint32_t>& columns : system.StaircaseRows)
        {
            BinaryEquation equation;
            equation.Columns = columns;
            st.BinaryRowReferences += columns.size();
            rows.push_back(std::move(equation));
        }
        for (const std::vector<uint32_t>& columns : system.DenseRowColumns)
        {
            BinaryEquation equation;
            equation.Columns = columns;
            st.BinaryRowReferences += columns.size();
            rows.push_back(std::move(equation));
        }
        for (const SolvePacket& packet : packets)
        {
            BinaryEquation equation;
            equation.Columns = GeneratePacketMatrixRow(
                K, P, packet.BlockId, config);
            if (equation.Columns.empty()) {
                return Wirehair_InvalidInput;
            }
            equation.Data = packet.Data;
            st.BinaryRowReferences += equation.Columns.size();
            rows.push_back(std::move(equation));
        }
        st.PacketRows = (uint32_t)packets.size();
        SolveClock::time_point phase_end = SolveClock::now();
        st.BuildNanoseconds = (uint64_t)
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                phase_end - phase_start).count();

        phase_start = phase_end;
        const PeelResult peel = PeelBinaryRows(L, rows);
        if (peel.PeelOrder.size() + peel.InactiveOrder.size() != L) {
            return Wirehair_Error;
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
            const BinaryEquation& equation = rows[peel.SolveRow[column]];
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
                    return Wirehair_Error;
                }
            }
            if (!VerifyPrecodeSolution(
                    system,
                    config,
                    packets,
                    values.data(),
                    block_bytes))
            {
                return Wirehair_Error;
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

        const auto insert_residual = [&]() -> bool {
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
                    block_bytes, st);
            }

            uint32_t pivot_column = R;
            for (uint32_t j = 0; j < R; ++j) {
                if (coeff[j] != 0u)
                {
                    pivot_column = j;
                    break;
                }
            }
            if (pivot_column == R) {
                return RowIsZero(rhs.data(), block_bytes);
            }
            if (have_pivot[pivot_column]) {
                return false;
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
                ++st.BlockMulAdds;
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
                    existing_coeff[k] ^=
                        gf256_mul(scale, coeff[k]);
                }
                AddScaledBlock(
                    pivot_rhs.data() + (size_t)existing * block_bytes,
                    scale, rhs.data(), block_bytes, st);
            }

            std::memcpy(
                pivot_coeff.data() + (size_t)pivot_column * R,
                coeff.data(), R);
            std::memcpy(
                pivot_rhs.data() + (size_t)pivot_column * block_bytes,
                rhs.data(), block_bytes);
            have_pivot[pivot_column] = 1u;
            ++rank;
            return true;
        };

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
                            const uint32_t bit = LowestSetBitIndex(word);
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
            if (!insert_residual()) {
                return Wirehair_Error;
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
                        (uint64_t)HeavyCoefficient(heavy, column, H) <<
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
                        const uint32_t bit = LowestSetBitIndex(word);
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
        std::vector<uint8_t> residue_bucket(block_bytes, 0u);
        for (uint32_t residue = 0; residue < window; ++residue)
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
            for (uint32_t heavy = 0; heavy < H; ++heavy)
            {
                AddScaledBlock(
                    heavy_rhs.data() + (size_t)heavy * block_bytes,
                    HeavyCoefficient(heavy, residue, H),
                    residue_bucket.data(),
                    block_bytes, st);
            }
        }
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
            if (!insert_residual()) {
                return Wirehair_Error;
            }
        }

        st.ResidualRank = rank;
        phase_end = SolveClock::now();
        st.ResidualNanoseconds = (uint64_t)
            std::chrono::duration_cast<std::chrono::nanoseconds>(
                phase_end - phase_start).count();
        if (rank < R) {
            if (stats) {
                *stats = st;
            }
            return Wirehair_NeedMore;
        }

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
        phase_start = phase_end;

        // Dependencies of a peeled column were resolved earlier, so forward
        // chronological substitution reconstructs every remaining value.
        for (uint32_t column : peel.PeelOrder)
        {
            uint8_t* value =
                values.data() + (size_t)column * block_bytes;
            const BinaryEquation& equation = rows[peel.SolveRow[column]];
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
            const uint8_t scale = HeavyCoefficient(heavy, column, H);
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
