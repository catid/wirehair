#include "WirehairV2PrecodeEncode.h"

#include "../WirehairTools.h"
#include "../gf256.h"

#include <algorithm>
#include <cstring>
#include <vector>

namespace wirehair_v2 {
namespace {

/**
    Block accumulator: first term is a copy, later terms XOR, an empty sum
    zero-fills.  Copies count as ops (Phase B convention), zero-fill does not.
*/
class BlockAccumulator
{
public:
    BlockAccumulator(uint8_t* dst, uint32_t bytes, uint64_t& ops)
        : Dst(dst), Bytes(bytes), Ops(ops)
    {
    }

    void Add(const uint8_t* src)
    {
        if (Any) {
            gf256_add_mem(Dst, src, (int)Bytes);
        }
        else {
            std::memcpy(Dst, src, Bytes);
            Any = true;
        }
        ++Ops;
    }

    void Finish()
    {
        if (!Any) {
            std::memset(Dst, 0, Bytes);
        }
    }

private:
    uint8_t* const Dst;
    const uint32_t Bytes;
    uint64_t& Ops;
    bool Any = false;
};

/// GF(2) rank of up to 64-wide bit-mask rows (destructive)
uint32_t MaskRank(std::vector<uint64_t>& masks, uint32_t width)
{
    uint32_t rank = 0;
    for (uint32_t col = 0; col < width; ++col)
    {
        const uint64_t bit = UINT64_C(1) << col;
        uint32_t pivot = rank;
        while (pivot < masks.size() && 0u == (masks[pivot] & bit)) {
            ++pivot;
        }
        if (pivot >= masks.size()) {
            continue;
        }
        std::swap(masks[rank], masks[pivot]);
        for (uint32_t r = 0; r < masks.size(); ++r)
        {
            if (r != rank && 0u != (masks[r] & bit)) {
                masks[r] ^= masks[rank];
            }
        }
        ++rank;
    }
    return rank;
}

/// Dense-column bit mask of one dense row (bits over [K + S, K + S + D2))
bool DenseColumnMask(
    const std::vector<uint32_t>& row_columns,
    uint32_t dense_base,
    uint32_t dense_rows,
    uint64_t& mask_out)
{
    uint64_t mask = 0;
    for (uint32_t col : row_columns)
    {
        if (col < dense_base) {
            continue;
        }
        if (col - dense_base >= dense_rows) {
            return false;
        }
        mask |= UINT64_C(1) << (col - dense_base);
    }
    mask_out = mask;
    return true;
}

} // namespace

bool DenseCornerInvertible(const PrecodeSystem& system)
{
    const uint32_t D2 = system.Params.DenseRows;
    if (D2 == 0u) {
        return true;
    }
    if (D2 > 64u || system.DenseRowColumns.size() != D2) {
        return false;
    }

    const uint32_t dense_base =
        system.Params.BlockCount + system.Params.Staircase;
    std::vector<uint64_t> masks(D2);
    for (uint32_t r = 0; r < D2; ++r) {
        if (!DenseColumnMask(
                system.DenseRowColumns[r], dense_base, D2, masks[r]))
        {
            return false;
        }
    }
    return MaskRank(masks, D2) == D2;
}

bool ComputePrecodeValues(
    const PrecodeSystem& system,
    const uint8_t* source_blocks,
    uint32_t block_bytes,
    uint8_t* parity_blocks,
    PrecodeEncodeStats* stats)
{
    const uint32_t K = system.Params.BlockCount;
    const uint32_t S = system.Params.Staircase;
    const uint32_t D2 = system.Params.DenseRows;
    const uint32_t H = system.Params.HeavyRows;
    const uint32_t dense_base = K + S;
    const uint32_t heavy_base = K + S + D2;

    PrecodeEncodeStats local_stats = {};
    PrecodeEncodeStats& st = stats ? *stats : local_stats;
    st = PrecodeEncodeStats{};

    if (!source_blocks || !parity_blocks ||
        block_bytes == 0u || block_bytes > 0x7fffffffu)
    {
        return false;
    }
    // Dense masks live in uint64_t; HeavyCoefficient requires H <= 128
    if (D2 > 64u || H > 128u ||
        system.StaircaseRows.size() != S ||
        system.DenseRowColumns.size() != D2)
    {
        return false;
    }

    const int bytes = (int)block_bytes;

    // Value of global column c: source for c < K, else already-computed
    // parity K + i.  Phase order (staircase -> dense -> heavy) guarantees
    // that every referenced parity is filled in before it is read.
    const auto column_value = [&](uint32_t c) -> const uint8_t* {
        return c < K ?
            source_blocks + (size_t)c * block_bytes :
            parity_blocks + (size_t)(c - K) * block_bytes;
    };

    // --- Staircase forward pass ---
    // Row j: sources + link K+j-1 + own column K+j sum to zero, so parity j
    // is the XOR of everything but the own column.  Rows are sorted, so the
    // link (already computed) precedes the own column.
    for (uint32_t j = 0; j < S; ++j)
    {
        const uint32_t own = K + j;
        BlockAccumulator acc(
            parity_blocks + (size_t)j * block_bytes,
            block_bytes, st.StaircaseBlockOps);
        for (const uint32_t col : system.StaircaseRows[j])
        {
            if (col == own) {
                continue;
            }
            if (col > own || (col >= K && col + 1u != own)) {
                return false; // not a staircase row (stray parity column)
            }
            acc.Add(column_value(col));
        }
        acc.Finish();
    }

    // --- Shuffle-2 dense rows ---
    // Transform to consecutive-row differences: diff row 0 = row 0, diff
    // row r = row r XOR row r-1 (unit bidiagonal row transform: same rank,
    // same solution).  Each difference is the flip pair, so the known part
    // of diff row r > 0 costs at most 2 block ops; dense-column flips only
    // toggle the row's corner mask.  Then Gauss-Jordan solves the D2 x D2
    // GF(2) corner with the known sums as block-valued RHS.
    if (D2 > 0u)
    {
        std::vector<uint8_t> rhs((size_t)D2 * block_bytes);
        std::vector<uint64_t> masks(D2);
        std::vector<uint32_t> diff;

        for (uint32_t r = 0; r < D2; ++r)
        {
            const std::vector<uint32_t>& cur = system.DenseRowColumns[r];
            const uint32_t* cols = cur.data();
            size_t count = cur.size();
            uint64_t mask = 0; // dense bits of THIS transformed row only
            if (r > 0u)
            {
                const std::vector<uint32_t>& prev =
                    system.DenseRowColumns[r - 1u];
                diff.clear();
                std::set_symmetric_difference(
                    prev.begin(), prev.end(),
                    cur.begin(), cur.end(),
                    std::back_inserter(diff));
                cols = diff.data();
                count = diff.size();
            }

            BlockAccumulator acc(
                rhs.data() + (size_t)r * block_bytes,
                block_bytes, st.DenseKnownBlockOps);
            for (size_t i = 0; i < count; ++i)
            {
                const uint32_t col = cols[i];
                if (col >= dense_base + D2) {
                    return false; // out-of-span dense column
                }
                if (col >= dense_base) {
                    mask ^= UINT64_C(1) << (col - dense_base);
                }
                else {
                    acc.Add(column_value(col));
                }
            }
            acc.Finish();
            masks[r] = mask;
        }

        // GF(2) Gauss-Jordan: singular corner => no (unique) solution for
        // generic source data => encoder-infeasible seed
        std::vector<uint32_t> pivot_row(D2);
        std::vector<uint8_t> used(D2, 0);
        for (uint32_t j = 0; j < D2; ++j)
        {
            const uint64_t bit = UINT64_C(1) << j;
            uint32_t p = D2;
            for (uint32_t r = 0; r < D2; ++r)
            {
                if (!used[r] && 0u != (masks[r] & bit)) {
                    p = r;
                    break;
                }
            }
            if (p >= D2) {
                return false; // singular dense corner
            }
            used[p] = 1;
            pivot_row[j] = p;
            for (uint32_t r = 0; r < D2; ++r)
            {
                if (r != p && 0u != (masks[r] & bit))
                {
                    masks[r] ^= masks[p];
                    gf256_add_mem(
                        rhs.data() + (size_t)r * block_bytes,
                        rhs.data() + (size_t)p * block_bytes,
                        bytes);
                    ++st.DenseSolveBlockOps;
                }
            }
        }
        for (uint32_t j = 0; j < D2; ++j)
        {
            std::memcpy(
                parity_blocks + (size_t)(S + j) * block_bytes,
                rhs.data() + (size_t)pivot_row[j] * block_bytes,
                block_bytes);
            ++st.DenseSolveBlockOps;
        }
    }

    // --- Cauchy heavy rows ---
    // Row r: sum over ALL L columns of HeavyCoefficient(r, c, H) * value(c)
    // is zero, so the H heavy parities solve the H x H Cauchy corner against
    // the known-column accumulation.  Coefficients depend only on
    // c mod (256 - H), so one H x window table serves every column.
    if (H > 0u)
    {
        const uint32_t window = 256u - H;
        std::vector<uint8_t> coef((size_t)H * window);
        for (uint32_t r = 0; r < H; ++r) {
            for (uint32_t m = 0; m < window; ++m) {
                coef[(size_t)r * window + m] = HeavyCoefficient(r, m, H);
            }
        }

        // Known parts.  Coefficients repeat with period `window`, so when
        // there are many more known columns than residues it is much
        // cheaper to XOR same-residue columns into buckets first and
        // muladd each bucket once per row: L XORs + H*window muladds
        // instead of H*L muladds (~6x fewer block ops at K=3200).  The
        // bucket buffer costs window*block_bytes, so huge blocks fall back
        // to the direct path rather than allocating hundreds of MB.
        std::vector<uint8_t> rhs((size_t)H * block_bytes, 0);
        const bool bucketed =
            heavy_base >= 2u * window &&
            (uint64_t)window * block_bytes <= (UINT64_C(64) << 20);
        if (bucketed)
        {
            std::vector<uint8_t> bucket((size_t)window * block_bytes, 0);
            uint32_t m = 0;
            for (uint32_t c = 0; c < heavy_base; ++c)
            {
                gf256_add_mem(
                    bucket.data() + (size_t)m * block_bytes,
                    column_value(c),
                    bytes);
                ++st.HeavyBucketXors;
                if (++m >= window) {
                    m = 0;
                }
            }
            for (uint32_t r = 0; r < H; ++r)
            {
                uint8_t* dst = rhs.data() + (size_t)r * block_bytes;
                for (uint32_t mm = 0; mm < window; ++mm)
                {
                    const uint8_t y = coef[(size_t)r * window + mm];
                    const uint8_t* v =
                        bucket.data() + (size_t)mm * block_bytes;
                    if (y == 1u) {
                        gf256_add_mem(dst, v, bytes);
                    }
                    else {
                        gf256_muladd_mem(dst, y, v, bytes);
                    }
                    ++st.HeavyMulAdds;
                }
            }
        }
        else
        {
            uint32_t m = 0;
            for (uint32_t c = 0; c < heavy_base; ++c)
            {
                const uint8_t* v = column_value(c);
                for (uint32_t r = 0; r < H; ++r)
                {
                    const uint8_t y = coef[(size_t)r * window + m];
                    uint8_t* dst = rhs.data() + (size_t)r * block_bytes;
                    if (y == 1u) {
                        gf256_add_mem(dst, v, bytes);
                    }
                    else {
                        gf256_muladd_mem(dst, y, v, bytes);
                    }
                    ++st.HeavyMulAdds;
                }
                if (++m >= window) {
                    m = 0;
                }
            }
        }

        // H x H corner over heavy columns heavy_base + j: consecutive GE
        // columns, distinct within one 244-window (H <= 128 <= window), so
        // the Cauchy submatrix is invertible; the pivot search below is
        // defensive only.
        std::vector<uint8_t> corner((size_t)H * H);
        for (uint32_t r = 0; r < H; ++r) {
            for (uint32_t j = 0; j < H; ++j) {
                corner[(size_t)r * H + j] =
                    coef[(size_t)r * window + (heavy_base + j) % window];
            }
        }

        std::vector<uint32_t> pivot_row(H);
        std::vector<uint8_t> used(H, 0);
        for (uint32_t j = 0; j < H; ++j)
        {
            uint32_t p = H;
            for (uint32_t r = 0; r < H; ++r)
            {
                if (!used[r] && corner[(size_t)r * H + j] != 0u) {
                    p = r;
                    break;
                }
            }
            if (p >= H) {
                return false; // cannot happen for a true Cauchy corner
            }
            used[p] = 1;
            pivot_row[j] = p;

            // Normalize the pivot row
            const uint8_t pivot_value = corner[(size_t)p * H + j];
            if (pivot_value != 1u)
            {
                const uint8_t inv = gf256_inv(pivot_value);
                for (uint32_t k = 0; k < H; ++k) {
                    corner[(size_t)p * H + k] =
                        gf256_mul(corner[(size_t)p * H + k], inv);
                }
                gf256_div_mem(
                    rhs.data() + (size_t)p * block_bytes,
                    rhs.data() + (size_t)p * block_bytes,
                    pivot_value, bytes);
                ++st.HeavySolveBlockOps;
            }

            // Eliminate column j everywhere else
            for (uint32_t r = 0; r < H; ++r)
            {
                const uint8_t scale = corner[(size_t)r * H + j];
                if (r == p || scale == 0u) {
                    continue;
                }
                for (uint32_t k = 0; k < H; ++k)
                {
                    corner[(size_t)r * H + k] = (uint8_t)(
                        corner[(size_t)r * H + k] ^
                        gf256_mul(scale, corner[(size_t)p * H + k]));
                }
                if (scale == 1u) {
                    gf256_add_mem(
                        rhs.data() + (size_t)r * block_bytes,
                        rhs.data() + (size_t)p * block_bytes,
                        bytes);
                }
                else {
                    gf256_muladd_mem(
                        rhs.data() + (size_t)r * block_bytes,
                        scale,
                        rhs.data() + (size_t)p * block_bytes,
                        bytes);
                }
                ++st.HeavySolveBlockOps;
            }
        }
        for (uint32_t j = 0; j < H; ++j)
        {
            std::memcpy(
                parity_blocks + (size_t)(S + D2 + j) * block_bytes,
                rhs.data() + (size_t)pivot_row[j] * block_bytes,
                block_bytes);
            ++st.HeavySolveBlockOps;
        }
    }

    return true;
}

} // namespace wirehair_v2
