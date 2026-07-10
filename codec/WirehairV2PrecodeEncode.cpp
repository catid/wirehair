#include "WirehairV2PrecodeEncode.h"

#include "../WirehairTools.h"
#include "../gf256.h"

#include "WirehairV2Plan.h"

#include <algorithm>
#include <cstring>
#include <iterator>
#include <limits>
#include <new>
#include <stdexcept>
#include <utility>
#include <vector>

namespace wirehair_v2 {
namespace {

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
thread_local int64_t AllocationFailureCountdown = -1;

void GuardedAllocation()
{
    if (AllocationFailureCountdown == 0) {
        throw std::bad_alloc();
    }
    if (AllocationFailureCountdown > 0) {
        --AllocationFailureCountdown;
    }
}
#else
void GuardedAllocation() {}
#endif

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

bool MessageBlockCount(
    uint64_t message_bytes,
    uint32_t block_bytes,
    uint32_t& block_count)
{
    block_count = 0u;
    if (message_bytes == 0u ||
        block_bytes == 0u ||
        block_bytes > 0x7fffffffu)
    {
        return false;
    }
    const uint64_t count =
        message_bytes / block_bytes +
        ((message_bytes % block_bytes) != 0u ? 1u : 0u);
    if (count < CAT_WIREHAIR_MIN_N || count > CAT_WIREHAIR_MAX_N) {
        return false;
    }
    block_count = (uint32_t)count;
    return true;
}

uint32_t SourceBlockBytes(
    uint64_t message_bytes,
    uint32_t block_bytes,
    uint32_t block_count,
    uint32_t block_id)
{
    if (block_id + 1u < block_count) {
        return block_bytes;
    }
    const uint64_t offset = (uint64_t)block_id * block_bytes;
    const uint64_t remaining = message_bytes - offset;
    return remaining < block_bytes ? (uint32_t)remaining : block_bytes;
}

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
        mask ^= UINT64_C(1) << (col - dense_base);
    }
    mask_out = mask;
    return true;
}

} // namespace

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
void SetAllocationFailureCountdownForTesting(int64_t countdown)
{
    AllocationFailureCountdown = countdown;
}
#endif

bool DenseCornerInvertible(const PrecodeSystem& system)
{
    const uint32_t D2 = system.Params.DenseRows;
    if (!ValidatePrecodeSystem(system)) {
        return false;
    }
    if (D2 == 0u) {
        return true;
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
    PrecodeEncodeStats local_stats = {};
    PrecodeEncodeStats& st = stats ? *stats : local_stats;
    st = PrecodeEncodeStats{};

    if (gf256_init() != 0) {
        return false;
    }
    if (!source_blocks || !parity_blocks ||
        block_bytes == 0u || block_bytes > 0x7fffffffu ||
        !ValidatePrecodeSystem(system))
    {
        return false;
    }

    const uint32_t K = system.Params.BlockCount;
    const uint32_t S = system.Params.Staircase;
    const uint32_t D2 = system.Params.DenseRows;
    const uint32_t H = system.Params.HeavyRows;
    const uint32_t dense_base = K + S;
    const uint32_t heavy_base = dense_base + D2;
    const uint64_t size_limit =
        (uint64_t)std::numeric_limits<size_t>::max();
    if ((uint64_t)K * block_bytes > size_limit ||
        ((uint64_t)S + D2 + H) * block_bytes > size_limit)
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

static bool ComputeRecoveryBlockImpl(
    const PrecodeSystem& system,
    const uint8_t* source_blocks,
    const uint8_t* parity_blocks,
    uint32_t block_bytes,
    const std::vector<uint32_t>& row_columns,
    uint8_t* block_out,
    uint64_t* block_ops_out,
    bool validate_system)
{
    const uint32_t K = system.Params.BlockCount;
    const uint32_t S = system.Params.Staircase;
    const uint32_t D2 = system.Params.DenseRows;
    const uint32_t H = system.Params.HeavyRows;
    const uint64_t precode_count = (uint64_t)S + D2 + H;
    const uint64_t total_columns = (uint64_t)K + precode_count;

    uint64_t local_ops = 0;
    uint64_t& ops = block_ops_out ? *block_ops_out : local_ops;
    ops = 0;

    if (!source_blocks || !parity_blocks || !block_out ||
        block_bytes == 0u || block_bytes > 0x7fffffffu ||
        (uint64_t)K * block_bytes >
            (uint64_t)std::numeric_limits<size_t>::max() ||
        precode_count * block_bytes >
            (uint64_t)std::numeric_limits<size_t>::max() ||
        (validate_system && !ValidatePrecodeSystem(system)))
    {
        return false;
    }

    const auto column_value = [&](uint32_t c) -> const uint8_t* {
        return c < K ?
            source_blocks + (size_t)c * block_bytes :
            parity_blocks + (size_t)(c - K) * block_bytes;
    };

    for (const uint32_t col : row_columns)
    {
        if (col >= total_columns) {
            return false;
        }
    }

    BlockAccumulator acc(block_out, block_bytes, ops);
    for (const uint32_t col : row_columns) {
        acc.Add(column_value(col));
    }
    acc.Finish();
    return true;
}

bool ComputeRecoveryBlock(
    const PrecodeSystem& system,
    const uint8_t* source_blocks,
    const uint8_t* parity_blocks,
    uint32_t block_bytes,
    const std::vector<uint32_t>& row_columns,
    uint8_t* block_out,
    uint64_t* block_ops_out)
{
    if (gf256_init() != 0) {
        return false;
    }
    return ComputeRecoveryBlockImpl(
        system, source_blocks, parity_blocks, block_bytes, row_columns,
        block_out, block_ops_out, true);
}

static bool ComputeEncodedBlockImpl(
    const PrecodeSystem& system,
    const PeelingCodec& codec,
    uint64_t row_seed,
    uint32_t mix_count,
    const uint8_t* source_blocks,
    const uint8_t* parity_blocks,
    uint32_t block_bytes,
    uint32_t block_id,
    uint8_t* block_out,
    uint64_t* block_ops_out,
    bool validate_system)
{
    const uint32_t K = system.Params.BlockCount;
    const uint32_t S = system.Params.Staircase;
    const uint32_t D2 = system.Params.DenseRows;
    const uint32_t H = system.Params.HeavyRows;
    const uint64_t precode_count = (uint64_t)S + D2 + H;

    uint64_t local_ops = 0;
    uint64_t& ops = block_ops_out ? *block_ops_out : local_ops;
    ops = 0;

    if (!source_blocks || !block_out ||
        block_bytes == 0u || block_bytes > 0x7fffffffu ||
        (uint64_t)K * block_bytes >
            (uint64_t)std::numeric_limits<size_t>::max() ||
        precode_count * block_bytes >
            (uint64_t)std::numeric_limits<size_t>::max() ||
        (validate_system && !ValidatePrecodeSystem(system)))
    {
        return false;
    }

    if (block_id < K)
    {
        std::memcpy(
            block_out,
            source_blocks + (size_t)block_id * block_bytes,
            block_bytes);
        ops = 1u;
        return true;
    }

    if (!parity_blocks) {
        return false;
    }
    const uint32_t recovery_index = block_id - K;

    const std::vector<uint32_t> row = GenerateRecoveryMatrixRow(
        codec, K, (uint32_t)precode_count, recovery_index, mix_count,
        row_seed);
    if (row.empty()) {
        return false;
    }
    return ComputeRecoveryBlockImpl(
        system, source_blocks, parity_blocks, block_bytes, row, block_out,
        &ops, false);
}

bool ComputeEncodedBlock(
    const PrecodeSystem& system,
    const PeelingCodec& codec,
    uint64_t row_seed,
    uint32_t mix_count,
    const uint8_t* source_blocks,
    const uint8_t* parity_blocks,
    uint32_t block_bytes,
    uint32_t block_id,
    uint8_t* block_out,
    uint64_t* block_ops_out)
{
    if (gf256_init() != 0) {
        return false;
    }
    return ComputeEncodedBlockImpl(
        system, codec, row_seed, mix_count, source_blocks, parity_blocks,
        block_bytes, block_id, block_out, block_ops_out, true);
}

PrecodeEncoder::PrecodeEncoder()
{
}

void PrecodeEncoder::Swap(PrecodeEncoder& other) noexcept
{
    using std::swap;
    swap(SystemValue.Params, other.SystemValue.Params);
    SystemValue.StaircaseRows.swap(other.SystemValue.StaircaseRows);
    SystemValue.DenseRowColumns.swap(other.SystemValue.DenseRowColumns);
    swap(CodecValue, other.CodecValue);
    swap(RowSeed, other.RowSeed);
    swap(MixCount, other.MixCount);
    swap(SourceBlocks, other.SourceBlocks);
    swap(BlockBytesValue, other.BlockBytesValue);
    ParityBlockStorage.swap(other.ParityBlockStorage);
    SolvedIntermediateStorage.swap(other.SolvedIntermediateStorage);
    swap(PacketConfigValue, other.PacketConfigValue);
    swap(StatsValue, other.StatsValue);
    swap(UsesPacketContract, other.UsesPacketContract);
    swap(Initialized, other.Initialized);
}

WirehairResult PrecodeEncoder::InitializeSolvedSystem(
    const PrecodeSystem& system,
    const PacketRowConfig& packet_config,
    std::vector<uint8_t>& intermediate_blocks,
    uint32_t block_bytes)
{
    const uint64_t precode_count_wide =
        (uint64_t)system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    if (precode_count_wide > UINT32_MAX ||
        !IsPacketRowDomainValid(
            system.Params.BlockCount,
            (uint32_t)precode_count_wide,
            packet_config.MixCount) ||
        !ValidatePrecodeSystem(system) ||
        block_bytes == 0u || block_bytes > 0x7fffffffu)
    {
        return Wirehair_InvalidInput;
    }
    const uint32_t L = system.Params.BlockCount + system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    const uint64_t expected = (uint64_t)L * block_bytes;
    if (expected > (uint64_t)((size_t)-1) ||
        intermediate_blocks.size() != (size_t)expected)
    {
        return Wirehair_InvalidInput;
    }

    try
    {
        PrecodeEncoder next;
        GuardedAllocation();
        next.SystemValue = system;
        next.RowSeed = packet_config.PeelSeed;
        next.MixCount = packet_config.MixCount;
        next.BlockBytesValue = block_bytes;
        next.PacketConfigValue = packet_config;
        next.SolvedIntermediateStorage.swap(intermediate_blocks);
        next.SourceBlocks = next.SolvedIntermediateStorage.data();
        next.UsesPacketContract = true;
        next.Initialized = true;
        Swap(next);
        return Wirehair_Success;
    }
    catch (const std::bad_alloc&) {
        return Wirehair_OOM;
    }
}

bool PrecodeEncoder::Initialize(
    const PrecodeSystem& system,
    const PeelingCodec& codec,
    uint64_t row_seed,
    uint32_t mix_count,
    const uint8_t* source_blocks,
    uint32_t block_bytes)
{
    return InitializeResult(
        system, codec, row_seed, mix_count, source_blocks, block_bytes) ==
        Wirehair_Success;
}

WirehairResult PrecodeEncoder::InitializeResult(
    const PrecodeSystem& system,
    const PeelingCodec& codec,
    uint64_t row_seed,
    uint32_t mix_count,
    const uint8_t* source_blocks,
    uint32_t block_bytes)
{
    if (gf256_init() != 0) {
        return Wirehair_UnsupportedPlatform;
    }
    if (!source_blocks || block_bytes == 0u || block_bytes > 0x7fffffffu)
    {
        return Wirehair_InvalidInput;
    }

    try
    {
        if (!ValidatePrecodeSystem(system)) {
            return Wirehair_InvalidInput;
        }
        const uint64_t parity_count_wide =
            (uint64_t)system.Params.Staircase +
            system.Params.DenseRows + system.Params.HeavyRows;
        if ((uint64_t)system.Params.BlockCount * block_bytes >
                (uint64_t)std::numeric_limits<size_t>::max() ||
            parity_count_wide == 0u ||
            parity_count_wide > (uint64_t)(size_t)-1 / block_bytes)
        {
            return Wirehair_InvalidInput;
        }

        GuardedAllocation();
        std::vector<uint8_t> parity(
            (size_t)parity_count_wide * (size_t)block_bytes);
        PrecodeEncodeStats stats = {};
        GuardedAllocation();
        if (!ComputePrecodeValues(
                system, source_blocks, block_bytes, parity.data(), &stats))
        {
            return Wirehair_BadDenseSeed;
        }

        GuardedAllocation();
        PrecodeSystem next_system = system;
        PrecodeEncoder next;
        next.SystemValue = std::move(next_system);
        next.CodecValue = codec;
        next.RowSeed = row_seed;
        next.MixCount = mix_count;
        next.SourceBlocks = source_blocks;
        next.BlockBytesValue = block_bytes;
        next.ParityBlockStorage.swap(parity);
        next.StatsValue = stats;
        next.Initialized = true;
        Swap(next);
        return Wirehair_Success;
    }
    catch (const std::bad_alloc&) {
        return Wirehair_OOM;
    }
    catch (const std::length_error&) {
        return Wirehair_OOM;
    }
}

bool PrecodeEncoder::Encode(
    uint32_t block_id,
    uint8_t* block_out,
    uint64_t* block_ops_out) const
{
    return EncodeResult(block_id, block_out, block_ops_out) == Wirehair_Success;
}

WirehairResult PrecodeEncoder::EncodeResult(
    uint32_t block_id,
    uint8_t* block_out,
    uint64_t* block_ops_out) const
{
    if (!Initialized || !block_out) {
        return Wirehair_InvalidInput;
    }
    try
    {
        if (UsesPacketContract)
        {
            GuardedAllocation();
            uint64_t local_ops = 0u;
            if (!EvaluatePacketBlockForValidatedSystem(
                    SystemValue,
                    PacketConfigValue,
                    SolvedIntermediateStorage.data(),
                    BlockBytesValue,
                    block_id,
                    block_out,
                    &local_ops))
            {
                return Wirehair_InvalidInput;
            }
            if (block_ops_out) {
                *block_ops_out = local_ops;
            }
            return Wirehair_Success;
        }
        if (block_id >= SystemValue.Params.BlockCount) {
            GuardedAllocation();
        }
        uint64_t local_ops = 0;
        if (!ComputeEncodedBlockImpl(
                SystemValue,
                CodecValue,
                RowSeed,
                MixCount,
                SourceBlocks,
                ParityBlockStorage.data(),
                BlockBytesValue,
                block_id,
                block_out,
                &local_ops,
                false))
        {
            return Wirehair_InvalidInput;
        }
        if (block_ops_out) {
            *block_ops_out = local_ops;
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

bool PrecodeEncoder::IsInitialized() const
{
    return Initialized;
}

uint32_t PrecodeEncoder::SourceBlockCount() const
{
    return Initialized ? SystemValue.Params.BlockCount : 0u;
}

uint32_t PrecodeEncoder::ParityBlockCount() const
{
    return Initialized ?
        SystemValue.Params.Staircase + SystemValue.Params.DenseRows +
            SystemValue.Params.HeavyRows :
        0u;
}

uint32_t PrecodeEncoder::BlockBytes() const
{
    return Initialized ? BlockBytesValue : 0u;
}

uint64_t PrecodeEncoder::RecoveryRowSeed() const
{
    return Initialized ? RowSeed : 0u;
}

uint32_t PrecodeEncoder::RecoveryMixCount() const
{
    return Initialized ? MixCount : 0u;
}

const PrecodeEncodeStats& PrecodeEncoder::EncodeStats() const
{
    return StatsValue;
}

const uint8_t* PrecodeEncoder::ParityBlocks() const
{
    if (!Initialized) {
        return nullptr;
    }
    if (UsesPacketContract) {
        return SolvedIntermediateStorage.data() +
            (size_t)SystemValue.Params.BlockCount * BlockBytesValue;
    }
    return ParityBlockStorage.data();
}

const uint8_t* PrecodeEncoder::IntermediateBlocks() const
{
    return Initialized && UsesPacketContract ?
        SolvedIntermediateStorage.data() : nullptr;
}

const PrecodeSystem& PrecodeEncoder::System() const
{
    return SystemValue;
}

MessagePrecodeEncoder::MessagePrecodeEncoder()
{
}

namespace {

uint64_t MessagePrecodeMatrixSeed(
    const SeedProfile& profile,
    const MessagePrecodeEncoderOptions& options)
{
    uint64_t salt = options.PrecodeSeedSalt;
    if (options.DenseIdentityCorner) {
        salt ^= UINT64_C(0x4f6eecb28d4a9137);
    }
    return MatrixSeedFromProfile(profile, 0u, salt);
}

uint32_t MessagePacketPeelSeed(
    const SeedProfile& profile,
    const MessagePrecodeEncoderOptions& options)
{
    return PacketPeelSeedFromProfile(profile, options.RecoveryRowSeedSalt);
}

} // namespace

bool HasMessagePrecodeContractState(const SeedProfile& profile)
{
    return profile.V2SeedSelected || profile.V2SeedAttempt != 0u ||
        profile.V2PrecodeContractVersion != 0u ||
        profile.V2PacketRowContractVersion != 0u ||
        profile.V2StaircaseCount != 0u ||
        profile.V2DenseRowCount != 0u ||
        profile.V2HeavyRowCount != 0u ||
        profile.V2SourceHits != 0u ||
        profile.V2PrecodeSeed != 0u ||
        profile.V2PacketPeelSeed != 0u ||
        profile.V2RecoveryMixCount != 0u ||
        profile.V2DenseIdentityCorner ||
        profile.V2PrecodeSeedSalt != 0u ||
        profile.V2RecoveryRowSeedSalt != 0u;
}

bool ResolveMessagePrecodeOptions(
    const SeedProfile& profile,
    const MessagePrecodeEncoderOptions* requested_options,
    MessagePrecodeEncoderOptions& resolved_options)
{
    if (!profile.V2SeedSelected)
    {
        if (HasMessagePrecodeContractState(profile)) {
            return false;
        }
        resolved_options = requested_options ? *requested_options :
            MessagePrecodeEncoderOptions();
        return resolved_options.RecoveryMixCount ==
            kCertifiedPacketMixCount;
    }

    if (profile.V2SeedAttempt >= kMaxPacketSeedAttempts ||
        profile.V2PrecodeContractVersion != kPrecodeContractVersion ||
        profile.V2PacketRowContractVersion != kPacketRowContractVersion ||
        profile.V2StaircaseCount == 0u ||
        profile.V2StaircaseCount != profile.DenseCount ||
        profile.V2DenseRowCount == 0u ||
        profile.V2HeavyRowCount == 0u ||
        profile.V2SourceHits == 0u ||
        profile.V2RecoveryMixCount != kCertifiedPacketMixCount)
    {
        return false;
    }

    MessagePrecodeEncoderOptions bound;
    bound.RecoveryMixCount = profile.V2RecoveryMixCount;
    bound.DenseIdentityCorner = profile.V2DenseIdentityCorner;
    bound.PrecodeSeedSalt = profile.V2PrecodeSeedSalt;
    bound.RecoveryRowSeedSalt = profile.V2RecoveryRowSeedSalt;
    if (requested_options &&
        (requested_options->RecoveryMixCount != bound.RecoveryMixCount ||
         requested_options->DenseIdentityCorner != bound.DenseIdentityCorner ||
         requested_options->PrecodeSeedSalt != bound.PrecodeSeedSalt ||
         requested_options->RecoveryRowSeedSalt != bound.RecoveryRowSeedSalt))
    {
        return false;
    }
    resolved_options = bound;
    return true;
}

bool ResolveMessagePrecodeConfiguration(
    const SeedProfile& profile,
    const MessagePrecodeEncoderOptions& options,
    PrecodeParams& params,
    PacketRowConfig& packet_config)
{
    MessagePrecodeEncoderOptions validated_options;
    if (!ResolveMessagePrecodeOptions(
            profile, &options, validated_options))
    {
        return false;
    }
    if (profile.DenseCount == 0u ||
        profile.DenseCount > wirehair::kMaxDenseCount ||
        (profile.BlockCount >= wirehair::kTinyTableCount &&
         profile.DenseCount % 4u != 2u))
    {
        return false;
    }
    if (profile.V2SeedSelected)
    {
        PrecodeParams expected = MakeCertifiedParams(
            profile.BlockCount,
            MessagePrecodeMatrixSeed(profile, validated_options));
        expected.Staircase = profile.DenseCount;
        expected.DenseIdentityCorner = validated_options.DenseIdentityCorner;
        expected = PrecodeParamsForAttempt(
            expected, profile.V2SeedAttempt);

        PacketRowConfig expected_packet;
        expected_packet.PeelSeed = MessagePacketPeelSeed(
            profile, validated_options);
        expected_packet.MixCount = validated_options.RecoveryMixCount;
        expected_packet = PacketConfigForAttempt(
            expected_packet, profile.V2SeedAttempt);

        if (profile.V2StaircaseCount != expected.Staircase ||
            profile.V2DenseRowCount != expected.DenseRows ||
            profile.V2HeavyRowCount != expected.HeavyRows ||
            profile.V2SourceHits != expected.SourceHits ||
            profile.V2PrecodeSeed != expected.Seed ||
            profile.V2PacketPeelSeed != expected_packet.PeelSeed ||
            profile.V2RecoveryMixCount != expected_packet.MixCount ||
            profile.V2DenseIdentityCorner !=
                expected.DenseIdentityCorner)
        {
            return false;
        }
        params = expected;
        packet_config = expected_packet;
    }
    else
    {
        params = MakeCertifiedParams(
            profile.BlockCount,
            MessagePrecodeMatrixSeed(profile, validated_options));
        params.Staircase = profile.DenseCount;
        params.DenseIdentityCorner = validated_options.DenseIdentityCorner;
        packet_config.PeelSeed = MessagePacketPeelSeed(
            profile, validated_options);
        packet_config.MixCount = validated_options.RecoveryMixCount;
    }

    const uint64_t precode_count_wide = (uint64_t)params.Staircase +
        params.DenseRows + params.HeavyRows;
    return precode_count_wide <= UINT32_MAX &&
        IsPacketRowDomainValid(
            params.BlockCount,
            (uint32_t)precode_count_wide,
            packet_config.MixCount);
}

void BindMessagePrecodeProfile(
    SeedProfile& profile,
    const MessagePrecodeEncoderOptions& options,
    const PrecodeSystem& system,
    const PacketRowConfig& packet_config,
    uint32_t packet_seed_attempt)
{
    profile.DenseCount = (uint16_t)system.Params.Staircase;
    profile.V2SeedSelected = true;
    profile.V2SeedAttempt = packet_seed_attempt;
    profile.V2PrecodeContractVersion = kPrecodeContractVersion;
    profile.V2PacketRowContractVersion = kPacketRowContractVersion;
    profile.V2StaircaseCount = system.Params.Staircase;
    profile.V2DenseRowCount = system.Params.DenseRows;
    profile.V2HeavyRowCount = system.Params.HeavyRows;
    profile.V2SourceHits = system.Params.SourceHits;
    profile.V2PrecodeSeed = system.Params.Seed;
    profile.V2PacketPeelSeed = packet_config.PeelSeed;
    profile.V2RecoveryMixCount = options.RecoveryMixCount;
    profile.V2DenseIdentityCorner = options.DenseIdentityCorner;
    profile.V2PrecodeSeedSalt = options.PrecodeSeedSalt;
    profile.V2RecoveryRowSeedSalt = options.RecoveryRowSeedSalt;
}

bool MessagePrecodeEncoder::Initialize(
    const void* message,
    uint64_t message_bytes,
    uint32_t block_bytes,
    const SeedProfile* seed_override,
    const MessagePrecodeEncoderOptions* options)
{
    return InitializeResult(
        message, message_bytes, block_bytes, seed_override, options) ==
        Wirehair_Success;
}

WirehairResult MessagePrecodeEncoder::InitializeResult(
    const void* message,
    uint64_t message_bytes,
    uint32_t block_bytes,
    const SeedProfile* seed_override,
    const MessagePrecodeEncoderOptions* options)
{
    uint32_t block_count = 0u;
    if (!message ||
        !MessageBlockCount(message_bytes, block_bytes, block_count))
    {
        return Wirehair_InvalidInput;
    }

    const SeedProfile profile = seed_override ? *seed_override :
        SelectSeedProfile(block_count, block_bytes);
    if (profile.BlockCount != block_count ||
        profile.BlockBytes != block_bytes)
    {
        return Wirehair_InvalidInput;
    }
    MessagePrecodeEncoderOptions opts;
    if (!ResolveMessagePrecodeOptions(profile, options, opts)) {
        return Wirehair_InvalidInput;
    }

    const uint64_t source_bytes =
        (uint64_t)block_count * (uint64_t)block_bytes;
    if (source_bytes > (uint64_t)((size_t)-1)) {
        return Wirehair_InvalidInput;
    }
    try
    {
        GuardedAllocation();
        std::vector<uint8_t> source_storage((size_t)source_bytes, 0u);
        std::memcpy(source_storage.data(), message, (size_t)message_bytes);

        PrecodeParams base_params;
        PacketRowConfig base_config;
        if (!ResolveMessagePrecodeConfiguration(
                profile, opts, base_params, base_config))
        {
            return Wirehair_InvalidInput;
        }
        uint32_t packet_seed_attempt = profile.V2SeedSelected ?
            profile.V2SeedAttempt : 0u;
        PrecodeParams params = base_params;
        PacketRowConfig packet_config = base_config;

        GuardedAllocation();
        PrecodeSystem system;
        if (!BuildPrecodeSystem(params, system)) {
            return Wirehair_InvalidInput;
        }

        GuardedAllocation();
        std::vector<SolvePacket> packets;
        packets.reserve(block_count);
        for (uint32_t block_id = 0; block_id < block_count; ++block_id)
        {
            SolvePacket packet;
            packet.BlockId = block_id;
            packet.Data = source_storage.data() +
                (size_t)block_id * block_bytes;
            packets.push_back(packet);
        }
        GuardedAllocation();
        std::vector<uint8_t> intermediate_blocks;
        PrecodeSolveStats solve_stats;
        WirehairResult solve_result = SolvePrecodeSystem(
            system,
            packet_config,
            packets,
            block_bytes,
            intermediate_blocks,
            &solve_stats);
        if (solve_result == Wirehair_NeedMore ||
            solve_result == Wirehair_Error)
        {
            if (profile.V2SeedSelected) {
                return solve_result == Wirehair_NeedMore ?
                    Wirehair_BadPeelSeed : solve_result;
            }
            PrecodeSystem selected_system;
            PacketRowConfig selected_config;
            const WirehairResult select_result = SelectSystematicConfiguration(
                base_params,
                base_config,
                selected_system,
                selected_config,
                &packet_seed_attempt);
            if (select_result != Wirehair_Success) {
                return select_result;
            }
            if (packet_seed_attempt == 0u) {
                return Wirehair_Error;
            }
            intermediate_blocks.clear();
            system = std::move(selected_system);
            solve_result = SolvePrecodeSystem(
                system,
                selected_config,
                packets,
                block_bytes,
                intermediate_blocks,
                &solve_stats);
            packet_config = selected_config;
        }
        if (solve_result != Wirehair_Success) {
            return solve_result;
        }
        solve_stats.PacketSeedAttempt = packet_seed_attempt;

        PrecodeEncoder next_encoder;
        const WirehairResult encoder_result = next_encoder.InitializeSolvedSystem(
            system,
            packet_config,
            intermediate_blocks,
            block_bytes);
        if (encoder_result != Wirehair_Success) {
            return encoder_result;
        }

        EncoderValue.Swap(next_encoder);
        ProfileValue = profile;
        BindMessagePrecodeProfile(
            ProfileValue,
            opts,
            system,
            packet_config,
            packet_seed_attempt);
        OptionsValue = opts;
        SolveStatsValue = solve_stats;
        MessageBytesValue = message_bytes;
        BlockBytesValue = block_bytes;
        Initialized = true;
        return Wirehair_Success;
    }
    catch (const std::bad_alloc&) {
        return Wirehair_OOM;
    }
    catch (const std::length_error&) {
        return Wirehair_OOM;
    }
}

bool MessagePrecodeEncoder::Encode(
    uint32_t block_id,
    uint8_t* block_out,
    uint32_t out_bytes,
    uint32_t* data_bytes_out,
    uint64_t* block_ops_out) const
{
    return EncodeResult(
        block_id, block_out, out_bytes, data_bytes_out, block_ops_out) ==
        Wirehair_Success;
}

WirehairResult MessagePrecodeEncoder::EncodeResult(
    uint32_t block_id,
    uint8_t* block_out,
    uint32_t out_bytes,
    uint32_t* data_bytes_out,
    uint64_t* block_ops_out) const
{
    if (!Initialized || !block_out || !data_bytes_out) {
        return Wirehair_InvalidInput;
    }

    const uint32_t K = ProfileValue.BlockCount;
    if (block_id < K)
    {
        const uint32_t bytes = SourceBlockBytes(
            MessageBytesValue, BlockBytesValue, K, block_id);
        if (out_bytes < bytes) {
            return Wirehair_InvalidInput;
        }
        if (bytes == BlockBytesValue)
        {
            const WirehairResult result = EncoderValue.EncodeResult(
                block_id, static_cast<uint8_t*>(block_out), block_ops_out);
            if (result != Wirehair_Success) {
                return result;
            }
        }
        else
        {
            try
            {
                std::vector<uint8_t> full_block(BlockBytesValue, 0u);
                uint64_t ops = 0u;
                const WirehairResult result = EncoderValue.EncodeResult(
                    block_id, full_block.data(), &ops);
                if (result != Wirehair_Success) {
                    return result;
                }
                std::memcpy(block_out, full_block.data(), bytes);
                if (block_ops_out) {
                    *block_ops_out = ops;
                }
            }
            catch (const std::bad_alloc&) {
                return Wirehair_OOM;
            }
        }
        *data_bytes_out = bytes;
        return Wirehair_Success;
    }

    if (out_bytes < BlockBytesValue) {
        return Wirehair_InvalidInput;
    }
    const WirehairResult result =
        EncoderValue.EncodeResult(block_id, block_out, block_ops_out);
    if (result != Wirehair_Success) {
        return result;
    }
    *data_bytes_out = BlockBytesValue;
    return Wirehair_Success;
}

bool MessagePrecodeEncoder::IsInitialized() const
{
    return Initialized;
}

uint64_t MessagePrecodeEncoder::MessageBytes() const
{
    return Initialized ? MessageBytesValue : 0u;
}

uint32_t MessagePrecodeEncoder::SourceBlockCount() const
{
    return Initialized ? ProfileValue.BlockCount : 0u;
}

uint32_t MessagePrecodeEncoder::BlockBytes() const
{
    return Initialized ? BlockBytesValue : 0u;
}

const SeedProfile& MessagePrecodeEncoder::Profile() const
{
    return ProfileValue;
}

const MessagePrecodeEncoderOptions& MessagePrecodeEncoder::Options() const
{
    return OptionsValue;
}

const PrecodeEncodeStats& MessagePrecodeEncoder::EncodeStats() const
{
    return EncoderValue.EncodeStats();
}

const PrecodeSolveStats& MessagePrecodeEncoder::SolveStats() const
{
    return SolveStatsValue;
}

const uint8_t* MessagePrecodeEncoder::IntermediateBlocks() const
{
    return Initialized ? EncoderValue.IntermediateBlocks() : nullptr;
}

const PrecodeEncoder& MessagePrecodeEncoder::BlockEncoder() const
{
    return EncoderValue;
}

} // namespace wirehair_v2
