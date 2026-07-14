#include "WirehairV2PrecodeEncode.h"
#include "WirehairV2Codec.h"
#include "WirehairV2Peel.h"

#include "../WirehairTools.h"
#include "../gf256.h"

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <utility>
#include <vector>

// Encoder value phase tests for the certified precode (wirehair-axd
// Phase 3, first slice):
//
// 1. Correctness on as-built systems: after ComputePrecodeValues succeeds,
//    every staircase and dense constraint row XORs to zero over the full
//    column values and every heavy row's GF(256) weighted sum is zero.
//    Seeds with a singular dense corner are skipped and counted, and the
//    success/failure verdict is cross-checked against the structure-only
//    DenseCornerInvertible on every seed.
// 2. Correctness on doctored systems (dense corner forced to identity) so
//    the dense-solve and heavy phases are exercised at every K even where
//    as-built corners are singular.
// 3. Feasibility measurement: dense-corner invertibility rate and rank
//    distribution over many seeds -- the Phase 4 seed-gating input.
// 4. Cost check against the certified model 2K + S - 1 + ceil(span/2) +
//    2*(D2 - 1) block ops for the staircase + dense known-part path.

namespace {

class MixedCoefficientPeriodScope
{
public:
    explicit MixedCoefficientPeriodScope(uint32_t period)
        : Previous(wirehair_v2::ActiveMixedCoefficientPeriod())
        , Valid(wirehair_v2::SetMixedCoefficientPeriodForTesting(period))
    {
    }

    ~MixedCoefficientPeriodScope()
    {
        (void)wirehair_v2::SetMixedCoefficientPeriodForTesting(Previous);
    }

    bool IsValid() const { return Valid; }

private:
    uint32_t Previous;
    bool Valid;
};

bool ParsePositiveU32(const char* text, uint32_t& out)
{
    if (!text || !*text || *text < '0' || *text > '9') {
        return false;
    }
    errno = 0;
    char* end = nullptr;
    const unsigned long value = std::strtoul(text, &end, 10);
    if (errno != 0 || !end || *end != '\0' ||
        value == 0 || value > UINT32_MAX)
    {
        return false;
    }
    out = (uint32_t)value;
    return true;
}

const uint8_t* ColumnValue(
    const wirehair_v2::PrecodeSystem& system,
    const uint8_t* source, const uint8_t* parity,
    uint32_t block_bytes, uint32_t c)
{
    const uint32_t K = system.Params.BlockCount;
    return c < K ?
        source + (size_t)c * block_bytes :
        parity + (size_t)(c - K) * block_bytes;
}

bool ManualPacketBlock(
    const wirehair_v2::PrecodeSystem& system,
    const wirehair_v2::PacketRowConfig& config,
    const uint8_t* intermediate,
    uint32_t block_bytes,
    uint32_t block_id,
    uint8_t* block_out)
{
    const uint32_t K = system.Params.BlockCount;
    const uint32_t P = system.Params.Staircase +
        system.Params.DenseRows + system.Params.HeavyRows;
    const uint32_t L = K + P;
    const std::vector<uint32_t> row =
        wirehair_v2::GeneratePacketMatrixRow(K, P, block_id, config);
    if (!intermediate || !block_out || row.empty()) {
        return false;
    }
    std::memset(block_out, 0, block_bytes);
    for (const uint32_t column : row)
    {
        if (column >= L) {
            return false;
        }
        const uint8_t* const value =
            intermediate + (size_t)column * block_bytes;
        for (uint32_t b = 0; b < block_bytes; ++b) {
            block_out[b] ^= value[b];
        }
    }
    return true;
}

/// Verify every constraint row of the system sums to zero over the values
bool VerifyValues(
    const wirehair_v2::PrecodeSystem& system,
    const uint8_t* source, const uint8_t* parity,
    uint32_t block_bytes, const char* tag)
{
    const uint32_t K = system.Params.BlockCount;
    const uint32_t S = system.Params.Staircase;
    const uint32_t D2 = system.Params.DenseRows;
    const uint32_t H = system.Params.HeavyRows;
    const uint32_t L = K + S + D2 + H;
    std::vector<uint8_t> acc(block_bytes);

    const auto is_zero = [&]() {
        for (uint32_t i = 0; i < block_bytes; ++i) {
            if (acc[i] != 0u) {
                return false;
            }
        }
        return true;
    };

    // Binary rows: staircase then dense
    for (uint32_t family = 0; family < 2u; ++family)
    {
        const std::vector<std::vector<uint32_t>>& rows =
            family == 0u ? system.StaircaseRows : system.DenseRowColumns;
        for (size_t r = 0; r < rows.size(); ++r)
        {
            std::memset(acc.data(), 0, block_bytes);
            for (const uint32_t col : rows[r])
            {
                gf256_add_mem(acc.data(),
                    ColumnValue(system, source, parity, block_bytes, col),
                    (int)block_bytes);
            }
            if (!is_zero())
            {
                std::fprintf(stderr,
                    "%s: K=%u bb=%u: %s row %zu does not sum to zero\n",
                    tag, K, block_bytes,
                    family == 0u ? "staircase" : "dense", r);
                return false;
            }
        }
    }

    // Completion rows: legacy uses GF(256) throughout; the mixed profile
    // embeds its first ten rows in the subfield and uses GF(2^16) for two.
    for (uint32_t r = 0; r < H; ++r)
    {
        std::memset(acc.data(), 0, block_bytes);
        for (uint32_t c = 0; c < L; ++c)
        {
            const uint32_t coefficient_column =
                system.Params.Field ==
                    wirehair_v2::CompletionField::MixedGF256GF16 ?
                c % wirehair_v2::ActiveMixedCoefficientPeriod() : c;
            const uint8_t* v =
                ColumnValue(system, source, parity, block_bytes, c);
            if (system.Params.Field ==
                    wirehair_v2::CompletionField::MixedGF256GF16 &&
                r >= wirehair_v2::kMixedGF256Rows)
            {
                const uint16_t y = wirehair_v2::MixedGF16Coefficient(
                    r - wirehair_v2::kMixedGF256Rows,
                    coefficient_column);
                if (!wirehair_v2::GF16MulAddMem(
                        acc.data(), y, v, block_bytes))
                {
                    return false;
                }
            }
            else
            {
                const uint8_t y = wirehair_v2::HeavyCoefficient(
                    r, coefficient_column, H);
                if (y == 1u) {
                    gf256_add_mem(acc.data(), v, (int)block_bytes);
                }
                else {
                    gf256_muladd_mem(acc.data(), y, v, (int)block_bytes);
                }
            }
        }
        if (!is_zero())
        {
            std::fprintf(stderr,
                "%s: K=%u bb=%u: heavy row %u weighted sum nonzero\n",
                tag, K, block_bytes, r);
            return false;
        }
    }
    return true;
}

void FillRandomBlocks(
    uint8_t* blocks, size_t bytes, uint64_t seed)
{
    wirehair::PCGRandom prng;
    prng.Seed(seed, 0xb10cda7au);
    for (size_t i = 0; i < bytes; ++i) {
        blocks[i] = (uint8_t)prng.Next();
    }
}

uint32_t MixedCornerRank(const std::vector<uint32_t>& columns)
{
    const uint32_t H = wirehair_v2::kMixedGF256Rows +
        wirehair_v2::kMixedGF16Rows;
    if (columns.size() != H) return 0u;
    std::vector<uint16_t> matrix((size_t)H * H);
    for (uint32_t r = 0; r < wirehair_v2::kMixedGF256Rows; ++r) {
        for (uint32_t j = 0; j < H; ++j) {
            matrix[(size_t)r * H + j] = wirehair_v2::HeavyCoefficient(
                r, columns[j], H);
        }
    }
    for (uint32_t er = 0; er < wirehair_v2::kMixedGF16Rows; ++er) {
        const uint32_t r = wirehair_v2::kMixedGF256Rows + er;
        for (uint32_t j = 0; j < H; ++j) {
            matrix[(size_t)r * H + j] =
                wirehair_v2::MixedGF16Coefficient(er, columns[j]);
        }
    }

    uint32_t rank = 0u;
    for (uint32_t col = 0; col < H; ++col)
    {
        uint32_t pivot = rank;
        while (pivot < H && matrix[(size_t)pivot * H + col] == 0u) {
            ++pivot;
        }
        if (pivot >= H) continue;
        for (uint32_t k = 0; k < H; ++k) {
            std::swap(
                matrix[(size_t)rank * H + k],
                matrix[(size_t)pivot * H + k]);
        }
        const uint16_t inverse = wirehair_v2::GF16Inverse(
            matrix[(size_t)rank * H + col]);
        for (uint32_t k = 0; k < H; ++k) {
            matrix[(size_t)rank * H + k] = wirehair_v2::GF16Multiply(
                matrix[(size_t)rank * H + k], inverse);
        }
        for (uint32_t r = 0; r < H; ++r)
        {
            if (r == rank) continue;
            const uint16_t scale = matrix[(size_t)r * H + col];
            for (uint32_t k = 0; k < H; ++k) {
                matrix[(size_t)r * H + k] ^=
                    wirehair_v2::GF16Multiply(
                        scale, matrix[(size_t)rank * H + k]);
            }
        }
        ++rank;
    }
    return rank;
}

bool TestMixedCornerRank()
{
    const uint32_t H = wirehair_v2::kMixedGF256Rows +
        wirehair_v2::kMixedGF16Rows;
    // The smaller periods exercised by the implementation oracles below are
    // experimental sweep points, not certified coefficient sets.  Only keep
    // candidates here after they pass this stronger arbitrary-column screen;
    // period 64 was rejected after producing a singular sampled corner.
    const uint32_t periods[] = {
        wirehair_v2::kMixedCoefficientPeriod, 96u
    };
    std::vector<uint32_t> columns(H);
    wirehair::PCGRandom prng;
    prng.Seed(UINT64_C(0x16c0a4e7), UINT64_C(0x244));
    for (const uint32_t period : periods)
    {
        for (uint32_t start = 0; start < period; ++start)
        {
            for (uint32_t j = 0; j < H; ++j) {
                columns[j] = (start + j) % period;
            }
            if (MixedCornerRank(columns) != H) {
                std::fprintf(stderr,
                    "mixed corner: period %u consecutive start %u is "
                    "singular\n",
                    period, start);
                return false;
            }
        }

        std::vector<uint32_t> deck(period);
        for (uint32_t trial = 0; trial < 10000u; ++trial)
        {
            for (uint32_t i = 0; i < period; ++i) deck[i] = i;
            for (uint32_t i = 0; i < H; ++i)
            {
                const uint32_t pick = i + prng.Next() % (period - i);
                std::swap(deck[i], deck[pick]);
                columns[i] = deck[i];
            }
            if (MixedCornerRank(columns) != H) {
                std::fprintf(stderr,
                    "mixed corner: period %u nonconsecutive trial %u is "
                    "singular\n",
                    period, trial);
                return false;
            }
        }
    }
    const uint32_t original_period =
        wirehair_v2::ActiveMixedCoefficientPeriod();
    if (wirehair_v2::SetMixedCoefficientPeriodForTesting(H - 1u) ||
        wirehair_v2::SetMixedCoefficientPeriodForTesting(
            wirehair_v2::kMixedCoefficientPeriod + 1u) ||
        wirehair_v2::ActiveMixedCoefficientPeriod() != original_period)
    {
        std::fprintf(stderr,
            "mixed corner: invalid period override was accepted\n");
        return false;
    }
    std::printf(
        "mixed 12x12 corner rank (periods 244/96, all starts + "
        "10000 samples each): PASS\n");
    return true;
}

/// GF(2) rank of the D2 x D2 dense corner (structure only)
uint32_t DenseCornerRank(const wirehair_v2::PrecodeSystem& system)
{
    const uint32_t D2 = system.Params.DenseRows;
    const uint32_t dense_base =
        system.Params.BlockCount + system.Params.Staircase;
    std::vector<uint64_t> masks(D2, 0);
    for (uint32_t r = 0; r < D2; ++r)
    {
        for (const uint32_t col : system.DenseRowColumns[r])
        {
            if (col >= dense_base) {
                masks[r] ^= UINT64_C(1) << (col - dense_base);
            }
        }
    }
    uint32_t rank = 0;
    for (uint32_t col = 0; col < D2; ++col)
    {
        const uint64_t bit = UINT64_C(1) << col;
        uint32_t pivot = rank;
        while (pivot < D2 && 0u == (masks[pivot] & bit)) {
            ++pivot;
        }
        if (pivot >= D2) {
            continue;
        }
        std::swap(masks[rank], masks[pivot]);
        for (uint32_t r = 0; r < D2; ++r) {
            if (r != rank && 0u != (masks[r] & bit)) {
                masks[r] ^= masks[rank];
            }
        }
        ++rank;
    }
    return rank;
}

bool TestCorrectnessAsBuilt()
{
    const uint32_t Ks[] = {2u, 16u, 1000u, 3200u};
    const uint32_t bbs[] = {1u, 13u, 1280u};
    const uint32_t seeds = 30u;
    bool ok = true;

    for (const uint32_t K : Ks)
    {
        for (const uint32_t bb : bbs)
        {
            uint32_t succeeded = 0, singular = 0;
            for (uint32_t seed = 0; seed < seeds; ++seed)
            {
                wirehair_v2::PrecodeSystem system;
                if (!BuildPrecodeSystem(
                        wirehair_v2::MakeCertifiedParams(K, seed), system))
                {
                    std::fprintf(stderr, "K=%u: build failed\n", K);
                    return false;
                }
                const uint32_t parity_count = system.Params.Staircase +
                    system.Params.DenseRows + system.Params.HeavyRows;
                std::vector<uint8_t> source((size_t)K * bb);
                std::vector<uint8_t> parity((size_t)parity_count * bb);
                FillRandomBlocks(source.data(), source.size(),
                    seed * 1000003ull + K * 17ull + bb);

                const bool feasible =
                    wirehair_v2::DenseCornerInvertible(system);
                const bool computed = wirehair_v2::ComputePrecodeValues(
                    system, source.data(), bb, parity.data());
                if (computed != feasible)
                {
                    std::fprintf(stderr,
                        "K=%u bb=%u seed=%u: ComputePrecodeValues=%d but "
                        "DenseCornerInvertible=%d\n",
                        K, bb, seed, (int)computed, (int)feasible);
                    return false;
                }
                if (!computed)
                {
                    ++singular;
                    continue;
                }
                ++succeeded;
                if (!VerifyValues(system, source.data(), parity.data(),
                        bb, "as-built"))
                {
                    ok = false;
                }
            }
            std::printf(
                "as-built K=%5u bb=%4u: %2u/%u seeds feasible "
                "(%u singular dense corners)\n",
                K, bb, succeeded, seeds, singular);
        }
    }
    return ok;
}

bool TestCorrectnessDoctored()
{
    const uint32_t Ks[] = {16u, 1000u, 3200u};
    const uint32_t bbs[] = {1u, 13u, 1280u};

    for (const uint32_t K : Ks)
    {
        for (const uint32_t bb : bbs)
        {
            for (uint32_t seed = 100; seed < 103u; ++seed)
            {
                wirehair_v2::PrecodeParams params =
                    wirehair_v2::MakeCertifiedParams(K, seed);
                params.DenseIdentityCorner = true;
                wirehair_v2::PrecodeSystem system;
                if (!BuildPrecodeSystem(params, system))
                {
                    std::fprintf(stderr, "K=%u: build failed\n", K);
                    return false;
                }
                if (!wirehair_v2::DenseCornerInvertible(system))
                {
                    std::fprintf(stderr,
                        "K=%u: doctored corner reported singular\n", K);
                    return false;
                }

                const uint32_t parity_count = system.Params.Staircase +
                    system.Params.DenseRows + system.Params.HeavyRows;
                std::vector<uint8_t> source((size_t)K * bb);
                std::vector<uint8_t> parity((size_t)parity_count * bb);
                FillRandomBlocks(source.data(), source.size(),
                    seed * 999983ull + K * 29ull + bb);

                if (!wirehair_v2::ComputePrecodeValues(
                        system, source.data(), bb, parity.data()))
                {
                    std::fprintf(stderr,
                        "K=%u bb=%u seed=%u: doctored encode failed\n",
                        K, bb, seed);
                    return false;
                }
                if (!VerifyValues(system, source.data(), parity.data(),
                        bb, "doctored"))
                {
                    return false;
                }
            }
        }
    }
    std::printf("doctored-corner correctness: PASS\n");
    return true;
}

bool TestFeasibility(uint32_t trials)
{
    const uint32_t Ks[] = {1000u, 3200u, 32000u};

    std::printf(
        "\ndense-corner feasibility (%u seeds per K; Phase 4 gating "
        "input):\n", trials);
    for (const uint32_t K : Ks)
    {
        uint32_t invertible = 0;
        uint32_t max_rank = 0;
        uint64_t rank_sum = 0;
        for (uint32_t seed = 0; seed < trials; ++seed)
        {
            wirehair_v2::PrecodeSystem system;
            if (!BuildPrecodeSystem(
                    wirehair_v2::MakeCertifiedParams(K, seed), system))
            {
                std::fprintf(stderr, "K=%u: build failed\n", K);
                return false;
            }
            const uint32_t rank = DenseCornerRank(system);
            rank_sum += rank;
            max_rank = std::max(max_rank, rank);
            const bool inv = rank == system.Params.DenseRows;
            if (inv != wirehair_v2::DenseCornerInvertible(system))
            {
                std::fprintf(stderr,
                    "K=%u seed=%u: DenseCornerInvertible disagrees with "
                    "rank\n", K, seed);
                return false;
            }
            if (inv) {
                ++invertible;
            }
        }
        std::printf(
            "  K=%5u: invertible %u/%u (%.4f%%), corner rank mean %.2f "
            "max %u of %u\n",
            K, invertible, trials, 100.0 * invertible / trials,
            (double)rank_sum / trials, max_rank, 12u);
    }
    return true;
}

bool TestCostModel()
{
    const uint32_t Ks[] = {16u, 1000u, 3200u, 10000u};
    bool ok = true;

    std::printf("\ncost check vs model min(N1,S)*K + S - 1 + ceil(span/2) + "
        "2*(D2-1):\n");
    for (const uint32_t K : Ks)
    {
        wirehair_v2::PrecodeSystem system;
        if (!BuildPrecodeSystem(
                wirehair_v2::MakeCertifiedParams(K, 0x5eedu), system))
        {
            std::fprintf(stderr, "K=%u: build failed\n", K);
            return false;
        }
        const uint32_t S = system.Params.Staircase;
        const uint32_t D2 = system.Params.DenseRows;
        const uint32_t span = K + S + D2;
        const uint32_t bb = 64u;

        const uint32_t parity_count =
            S + D2 + system.Params.HeavyRows;
        std::vector<uint8_t> source((size_t)K * bb);
        std::vector<uint8_t> parity((size_t)parity_count * bb);
        FillRandomBlocks(source.data(), source.size(), K);

        // Stats for the staircase + dense known-part path are recorded even
        // when the dense corner is singular (the common case at these K)
        wirehair_v2::PrecodeEncodeStats stats;
        (void)wirehair_v2::ComputePrecodeValues(
            system, source.data(), bb, parity.data(), &stats);

        const uint32_t hits =
            system.Params.SourceHits < S ? system.Params.SourceHits : S;
        const uint64_t staircase_model = (uint64_t)hits * K + S - 1u;
        const uint64_t dense_model =
            ((span + 1u) >> 1) + 2ull * (D2 - 1u);
        const uint64_t measured =
            stats.StaircaseBlockOps + stats.DenseKnownBlockOps;
        const uint64_t model = staircase_model + dense_model;

        std::printf(
            "  K=%5u: staircase %llu (model %llu), dense known %llu "
            "(model %llu), total %llu vs %llu (%+lld)\n",
            K,
            (unsigned long long)stats.StaircaseBlockOps,
            (unsigned long long)staircase_model,
            (unsigned long long)stats.DenseKnownBlockOps,
            (unsigned long long)dense_model,
            (unsigned long long)measured,
            (unsigned long long)model,
            (long long)(measured - model));

        if (stats.StaircaseBlockOps != staircase_model)
        {
            std::fprintf(stderr,
                "K=%u: staircase ops off model\n", K);
            ok = false;
        }
        // The dense known part undershoots the model by exactly the number
        // of dense-column touches (row-0 dense columns + dense-hitting
        // flips), each of which costs a mask toggle instead of a block op
        if (stats.DenseKnownBlockOps > dense_model ||
            dense_model - stats.DenseKnownBlockOps > D2 + 2u * (D2 - 1u))
        {
            std::fprintf(stderr,
                "K=%u: dense known-part ops outside model window\n", K);
            ok = false;
        }
    }
    return ok;
}

bool TestHeavyResidueDispatch()
{
    const uint32_t K = 500u;
    const uint32_t bb = 37u;
    wirehair_v2::PrecodeParams params =
        wirehair_v2::MakeCertifiedParams(K, UINT64_C(0x57ea4));
    params.DenseIdentityCorner = true;
    wirehair_v2::PrecodeSystem system;
    if (!BuildPrecodeSystem(params, system)) {
        std::fprintf(stderr, "heavy residue dispatch: build failed\n");
        return false;
    }
    const uint32_t parity_count = params.Staircase + params.DenseRows +
        params.HeavyRows;
    std::vector<uint8_t> source((size_t)K * bb);
    std::vector<uint8_t> full((size_t)parity_count * bb);
    std::vector<uint8_t> streamed(full.size());
    FillRandomBlocks(source.data(), source.size(), UINT64_C(0x57ea45eed));

    wirehair_v2::PrecodeEncodeStats full_stats;
    wirehair_v2::PrecodeEncodeStats streamed_stats;
    wirehair_v2::SetHeavyBucketStorageLimitForTesting(UINT64_MAX);
    const bool full_ok = wirehair_v2::ComputePrecodeValues(
        system, source.data(), bb, full.data(), &full_stats);
    wirehair_v2::SetHeavyBucketStorageLimitForTesting(0u);
    const bool streamed_ok = wirehair_v2::ComputePrecodeValues(
        system, source.data(), bb, streamed.data(), &streamed_stats);
    wirehair_v2::SetHeavyBucketStorageLimitForTesting(UINT64_C(64) << 20);

    const bool stats_equal =
        full_stats.StaircaseBlockOps == streamed_stats.StaircaseBlockOps &&
        full_stats.DenseKnownBlockOps == streamed_stats.DenseKnownBlockOps &&
        full_stats.DenseSolveBlockOps == streamed_stats.DenseSolveBlockOps &&
        full_stats.HeavyBucketXors == streamed_stats.HeavyBucketXors &&
        full_stats.HeavyMulAdds == streamed_stats.HeavyMulAdds &&
        full_stats.HeavySolveBlockOps == streamed_stats.HeavySolveBlockOps;
    const uint32_t heavy_base = K + params.Staircase + params.DenseRows;
    const uint32_t window = 256u - params.HeavyRows;
    if (!full_ok || !streamed_ok || full != streamed || !stats_equal ||
        streamed_stats.HeavyBucketXors != heavy_base ||
        streamed_stats.HeavyMulAdds != (uint64_t)params.HeavyRows * window ||
        !VerifyValues(
            system, source.data(), streamed.data(), bb,
            "heavy residue streaming"))
    {
        std::fprintf(stderr,
            "heavy residue dispatch: full/streaming mismatch\n");
        return false;
    }
    std::printf("heavy residue full/streaming differential: PASS\n");
    return true;
}

bool StatsAreZero(const wirehair_v2::PrecodeEncodeStats& stats);

bool TestMixedCompletionForPeriod(uint32_t period)
{
    MixedCoefficientPeriodScope period_scope(period);
    if (!period_scope.IsValid()) {
        return false;
    }
    // Exercise the exact production geometry (default full-span dense
    // corner).  Most such direct-source corners are singular by design, so
    // select a deterministically feasible tiny seed; the public packet solve
    // covers larger systems independently.
    bool production_geometry_checked = false;
    for (uint64_t seed = 0; seed < 1000u && !production_geometry_checked;
         ++seed)
    {
        wirehair_v2::PrecodeParams params =
            wirehair_v2::MakeMixedParams(2u, seed);
        wirehair_v2::PrecodeSystem system;
        if (!BuildPrecodeSystem(params, system) ||
            !wirehair_v2::DenseCornerInvertible(system))
        {
            continue;
        }
        const uint32_t bb = 16u;
        const uint32_t parity_count = params.Staircase +
            params.DenseRows + params.HeavyRows;
        std::vector<uint8_t> source((size_t)params.BlockCount * bb);
        std::vector<uint8_t> parity((size_t)parity_count * bb);
        FillRandomBlocks(source.data(), source.size(), seed ^ UINT64_C(0xcafe));
        if (!wirehair_v2::ComputePrecodeValues(
                system, source.data(), bb, parity.data()) ||
            !VerifyValues(
                system, source.data(), parity.data(), bb,
                "mixed production geometry"))
        {
            std::fprintf(stderr,
                "mixed completion: production geometry oracle failed\n");
            return false;
        }
        production_geometry_checked = true;
    }
    if (!production_geometry_checked) {
        std::fprintf(stderr,
            "mixed completion: no feasible production geometry seed\n");
        return false;
    }

    const uint32_t Ks[] = {16u, 500u, 1000u};
    const uint32_t bbs[] = {2u, 16u, 1280u};
    for (uint32_t K : Ks)
    {
        for (uint32_t bb : bbs)
        {
            wirehair_v2::PrecodeParams params =
                wirehair_v2::MakeMixedParams(K, UINT64_C(0x16c0de) + K + bb);
            params.DenseIdentityCorner = true;
            wirehair_v2::PrecodeSystem system;
            if (!BuildPrecodeSystem(params, system)) {
                std::fprintf(stderr, "mixed completion: build failed\n");
                return false;
            }
            const uint32_t parity_count = params.Staircase +
                params.DenseRows + params.HeavyRows;
            std::vector<uint8_t> source((size_t)K * bb);
            std::vector<uint8_t> full((size_t)parity_count * bb, 0xa5u);
            std::vector<uint8_t> streamed(full.size(), 0x5au);
            FillRandomBlocks(
                source.data(), source.size(), UINT64_C(0x16f1e1d) ^ K ^ bb);

            wirehair_v2::PrecodeEncodeStats full_stats;
            wirehair_v2::PrecodeEncodeStats streamed_stats;
            wirehair_v2::SetHeavyBucketStorageLimitForTesting(UINT64_MAX);
            const bool full_ok = wirehair_v2::ComputePrecodeValues(
                system, source.data(), bb, full.data(), &full_stats);
            wirehair_v2::SetHeavyBucketStorageLimitForTesting(0u);
            const bool streamed_ok = wirehair_v2::ComputePrecodeValues(
                system, source.data(), bb, streamed.data(), &streamed_stats);
            wirehair_v2::SetHeavyBucketStorageLimitForTesting(
                UINT64_C(64) << 20);

            const uint32_t heavy_base =
                K + params.Staircase + params.DenseRows;
            const uint32_t coefficient_period =
                wirehair_v2::ActiveMixedCoefficientPeriod();
            const bool bucketed =
                heavy_base >= 2u * coefficient_period;
            const uint64_t residues = bucketed ?
                coefficient_period : heavy_base;
            if (!full_ok || !streamed_ok || full != streamed ||
                !VerifyValues(
                    system, source.data(), streamed.data(), bb, "mixed") ||
                full_stats.HeavyBucketXors !=
                    (bucketed ? heavy_base : 0u) ||
                full_stats.HeavyMulAdds != params.HeavyRows * residues ||
                full_stats.MixedGF16MulAdds !=
                    wirehair_v2::kMixedGF16Rows * residues ||
                full_stats.MixedPlaneConversions !=
                    residues + wirehair_v2::kMixedGF256Rows +
                        params.HeavyRows ||
                full_stats.HeavyBucketXors !=
                    streamed_stats.HeavyBucketXors ||
                full_stats.HeavyMulAdds != streamed_stats.HeavyMulAdds ||
                full_stats.HeavySolveBlockOps !=
                    streamed_stats.HeavySolveBlockOps ||
                full_stats.MixedGF16MulAdds !=
                    streamed_stats.MixedGF16MulAdds ||
                full_stats.MixedGF16SolveBlockOps !=
                    streamed_stats.MixedGF16SolveBlockOps ||
                full_stats.MixedPlaneConversions !=
                    streamed_stats.MixedPlaneConversions)
            {
                std::fprintf(stderr,
                    "mixed completion: K=%u bb=%u differential/oracle failed\n",
                    K, bb);
                return false;
            }

            if (K == 16u && bb == 16u)
            {
                const wirehair_v2::PeelingCodec codec =
                    wirehair_v2::MakePeelingCodec(
                        wirehair_v2::PeelStructure::LtM1C32,
                        wirehair_v2::PeelSolver::KsBmaxTop16);
                wirehair_v2::PrecodeEncoder encoder;
                if (encoder.InitializeResult(
                        system, codec, UINT64_C(0x51eed), 5u,
                        source.data(), bb) != Wirehair_Success ||
                    encoder.System().Params.Field !=
                        wirehair_v2::CompletionField::MixedGF256GF16 ||
                    std::memcmp(
                        encoder.ParityBlocks(), full.data(), full.size()) != 0 ||
                    encoder.EncodeStats().MixedGF16MulAdds !=
                        full_stats.MixedGF16MulAdds ||
                    encoder.EncodeStats().MixedGF16SolveBlockOps !=
                        full_stats.MixedGF16SolveBlockOps)
                {
                    std::fprintf(stderr,
                        "mixed completion: stateful encoder mismatch\n");
                    return false;
                }
            }
        }
    }

    wirehair_v2::PrecodeParams params =
        wirehair_v2::MakeMixedParams(16u, UINT64_C(0x0dd));
    params.DenseIdentityCorner = true;
    wirehair_v2::PrecodeSystem system;
    if (!BuildPrecodeSystem(params, system)) return false;
    const uint32_t parity_count =
        params.Staircase + params.DenseRows + params.HeavyRows;
    std::vector<uint8_t> source((size_t)params.BlockCount * 3u, 0x11u);
    std::vector<uint8_t> parity((size_t)parity_count * 3u, 0x5au);
    const std::vector<uint8_t> before = parity;
    wirehair_v2::PrecodeEncodeStats stats;
    if (wirehair_v2::ComputePrecodeValues(
            system, source.data(), 3u, parity.data(), &stats) ||
        parity != before || !StatsAreZero(stats))
    {
        std::fprintf(stderr,
            "mixed completion: odd block rejection was not transactional\n");
        return false;
    }
    const wirehair_v2::PeelingCodec codec = wirehair_v2::MakePeelingCodec(
        wirehair_v2::PeelStructure::LtM1C32,
        wirehair_v2::PeelSolver::KsBmaxTop16);
    wirehair_v2::PrecodeEncoder encoder;
    if (encoder.InitializeResult(
            system, codec, UINT64_C(0x51eed), 5u,
            source.data(), 3u) != Wirehair_InvalidInput ||
        encoder.IsInitialized())
    {
        std::fprintf(stderr,
            "mixed completion: odd stateful encoder classification failed\n");
        return false;
    }

    std::printf(
        "mixed completion scalar/full/streamed oracle period=%u: PASS\n",
        period);
    return true;
}

bool TestMixedCompletion()
{
    const uint32_t periods[] = {
        wirehair_v2::kMixedCoefficientPeriod, 96u, 64u, 32u
    };
    for (const uint32_t period : periods) {
        if (!TestMixedCompletionForPeriod(period)) {
            return false;
        }
    }
    return true;
}

bool TestMalformedDenseCorner()
{
    wirehair_v2::PrecodeSystem system;
    system.Params.BlockCount = 2u;
    system.Params.Staircase = 1u;
    system.Params.DenseRows = 1u;
    system.Params.HeavyRows = 0u;
    system.Params.SourceHits = 1u;
    system.Params.DenseIdentityCorner = false;
    system.Params.Seed = 0u;
    system.StaircaseRows.push_back(std::vector<uint32_t>{0u, 2u});

    // dense_base is K + S = 3.  Column 67 is dense_base + 64, which used
    // to shift by 64 in DenseCornerInvertible before row validation.  Keep
    // a valid dense column at the tail so the test also catches validators
    // that only scan the sorted dense tail and miss malformed unsorted rows.
    system.DenseRowColumns.push_back(std::vector<uint32_t>{67u, 0u, 3u});
    if (wirehair_v2::DenseCornerInvertible(system)) {
        std::fprintf(stderr,
            "malformed dense corner should not be reported invertible\n");
        return false;
    }

    system.DenseRowColumns[0] = std::vector<uint32_t>{3u, 3u};
    if (wirehair_v2::DenseCornerInvertible(system)) {
        std::fprintf(stderr,
            "duplicate dense columns should cancel in GF(2)\n");
        return false;
    }
    return true;
}

bool StatsAreZero(const wirehair_v2::PrecodeEncodeStats& stats)
{
    return stats.StaircaseBlockOps == 0u &&
        stats.DenseKnownBlockOps == 0u &&
        stats.DenseSolveBlockOps == 0u &&
        stats.HeavyBucketXors == 0u &&
        stats.HeavyMulAdds == 0u &&
        stats.HeavySolveBlockOps == 0u &&
        stats.MixedGF16MulAdds == 0u &&
        stats.MixedGF16SolveBlockOps == 0u &&
        stats.MixedPlaneConversions == 0u;
}

bool TestStrictSystemValidation()
{
    wirehair_v2::PrecodeParams params =
        wirehair_v2::MakeCertifiedParams(16u, UINT64_C(0x51a1d));
    params.DenseIdentityCorner = true;
    wirehair_v2::PrecodeSystem valid;
    if (!BuildPrecodeSystem(params, valid) ||
        !wirehair_v2::ValidatePrecodeSystem(valid))
    {
        std::fprintf(stderr, "strict validation: valid builder system failed\n");
        return false;
    }

    std::vector<wirehair_v2::PrecodeSystem> invalid;
    wirehair_v2::PrecodeSystem bad;

    bad = valid;
    bad.Params.BlockCount = 0u;
    invalid.push_back(bad);
    bad = valid;
    bad.Params.Staircase = UINT32_MAX;
    invalid.push_back(bad);
    bad = valid;
    bad.Params.DenseRows = UINT32_MAX;
    invalid.push_back(bad);
    bad = valid;
    bad.Params.HeavyRows = 129u;
    invalid.push_back(bad);
    bad = valid;
    bad.StaircaseRows.pop_back();
    invalid.push_back(bad);
    bad = valid;
    bad.DenseRowColumns.pop_back();
    invalid.push_back(bad);

    const uint32_t K = valid.Params.BlockCount;
    const uint32_t S = valid.Params.Staircase;
    const uint32_t dense_base = K + S;
    const uint32_t binary_span = dense_base + valid.Params.DenseRows;

    bad = valid;
    bad.StaircaseRows[0].erase(
        std::find(bad.StaircaseRows[0].begin(),
            bad.StaircaseRows[0].end(), K));
    invalid.push_back(bad);
    bad = valid;
    bad.StaircaseRows[1].erase(
        std::find(bad.StaircaseRows[1].begin(),
            bad.StaircaseRows[1].end(), K));
    invalid.push_back(bad);
    bad = valid;
    bad.StaircaseRows[1].push_back(K + 2u);
    std::sort(bad.StaircaseRows[1].begin(), bad.StaircaseRows[1].end());
    invalid.push_back(bad);
    bad = valid;
    bad.StaircaseRows[0].insert(
        bad.StaircaseRows[0].begin(), bad.StaircaseRows[0][0]);
    invalid.push_back(bad);
    bad = valid;
    std::swap(bad.StaircaseRows[0][0], bad.StaircaseRows[0][1]);
    invalid.push_back(bad);

    bad = valid;
    bad.DenseRowColumns[0].insert(
        bad.DenseRowColumns[0].begin(), bad.DenseRowColumns[0][0]);
    invalid.push_back(bad);
    bad = valid;
    std::swap(bad.DenseRowColumns[0][0], bad.DenseRowColumns[0][1]);
    invalid.push_back(bad);
    bad = valid;
    bad.DenseRowColumns[0].push_back(binary_span);
    invalid.push_back(bad);
    bad = valid;
    bad.DenseRowColumns[0].back() = dense_base + 1u;
    std::sort(
        bad.DenseRowColumns[0].begin(), bad.DenseRowColumns[0].end());
    invalid.push_back(bad);
    bad = valid;
    bad.DenseRowColumns[1].erase(bad.DenseRowColumns[1].begin());
    invalid.push_back(bad);

    wirehair::PCGRandom mutation_rng;
    mutation_rng.Seed(UINT64_C(0xf0225eed), UINT64_C(0x51a1d));
    for (uint32_t trial = 0; trial < 256u; ++trial)
    {
        bad = valid;
        const uint32_t row = mutation_rng.Next() % S;
        switch (mutation_rng.Next() % 8u)
        {
        case 0:
            bad.Params.BlockCount = mutation_rng.Next() & 1u;
            break;
        case 1:
            bad.StaircaseRows[row].erase(std::find(
                bad.StaircaseRows[row].begin(),
                bad.StaircaseRows[row].end(), K + row));
            break;
        case 2:
            bad.StaircaseRows[row].insert(
                bad.StaircaseRows[row].begin(),
                bad.StaircaseRows[row][0]);
            break;
        case 3:
            bad.DenseRowColumns[mutation_rng.Next() % valid.Params.DenseRows]
                .push_back(binary_span);
            break;
        case 4:
        {
            std::vector<uint32_t>& dense = bad.DenseRowColumns[
                mutation_rng.Next() % valid.Params.DenseRows];
            std::swap(dense[0], dense[1]);
            break;
        }
        case 5:
        {
            const uint32_t dense_row =
                mutation_rng.Next() % valid.Params.DenseRows;
            std::vector<uint32_t>& dense = bad.DenseRowColumns[dense_row];
            dense.erase(std::find(
                dense.begin(), dense.end(), dense_base + dense_row));
            break;
        }
        case 6:
            bad.Params.HeavyRows = 129u + mutation_rng.Next() % 100u;
            break;
        default:
            bad.DenseRowColumns.pop_back();
            break;
        }
        invalid.push_back(bad);
    }

    const uint32_t bb = 7u;
    const uint32_t parity_count = S + valid.Params.DenseRows +
        valid.Params.HeavyRows;
    std::vector<uint8_t> source((size_t)K * bb, 0x3cu);
    for (size_t i = 0; i < invalid.size(); ++i)
    {
        std::vector<uint8_t> parity((size_t)parity_count * bb, 0xa5u);
        const std::vector<uint8_t> before = parity;
        wirehair_v2::PrecodeEncodeStats stats;
        stats.StaircaseBlockOps = UINT64_MAX;
        if (wirehair_v2::ValidatePrecodeSystem(invalid[i]) ||
            wirehair_v2::DenseCornerInvertible(invalid[i]) ||
            wirehair_v2::ComputePrecodeValues(
                invalid[i], source.data(), bb, parity.data(), &stats) ||
            parity != before || !StatsAreZero(stats))
        {
            std::fprintf(stderr,
                "strict validation: malformed case %zu was accepted/written\n",
                i);
            return false;
        }
    }

    wirehair_v2::PrecodeSystem singular;
    if (!BuildPrecodeSystem(
            wirehair_v2::MakeCertifiedParams(16u, 0u), singular) ||
        !wirehair_v2::ValidatePrecodeSystem(singular) ||
        wirehair_v2::DenseCornerInvertible(singular))
    {
        std::fprintf(stderr,
            "strict validation: expected structurally valid singular system\n");
        return false;
    }
    std::vector<uint8_t> singular_parity((size_t)parity_count * bb, 0xa5u);
    wirehair_v2::PrecodeEncodeStats singular_stats;
    if (wirehair_v2::ComputePrecodeValues(
            singular, source.data(), bb, singular_parity.data(),
            &singular_stats) ||
        singular_stats.StaircaseBlockOps == 0u ||
        singular_stats.DenseKnownBlockOps == 0u)
    {
        std::fprintf(stderr,
            "strict validation: singular system lost partial cost behavior\n");
        return false;
    }

    wirehair::PCGRandom property_rng;
    property_rng.Seed(UINT64_C(0xc0de51a1d), UINT64_C(0x5eed));
    for (uint32_t trial = 0; trial < 32u; ++trial)
    {
        const uint32_t random_k = 16u + property_rng.Next() % 1009u;
        wirehair_v2::PrecodeParams random_params =
            wirehair_v2::MakeCertifiedParams(random_k, property_rng.Next());
        random_params.DenseIdentityCorner = true;
        wirehair_v2::PrecodeSystem random_system;
        if (!BuildPrecodeSystem(random_params, random_system) ||
            !wirehair_v2::ValidatePrecodeSystem(random_system))
        {
            std::fprintf(stderr,
                "strict validation: random feasible build failed\n");
            return false;
        }
        const uint32_t random_parity_count = random_params.Staircase +
            random_params.DenseRows + random_params.HeavyRows;
        std::vector<uint8_t> random_source((size_t)random_k * 3u);
        std::vector<uint8_t> random_parity(
            (size_t)random_parity_count * 3u);
        FillRandomBlocks(
            random_source.data(), random_source.size(), property_rng.Next());
        if (!wirehair_v2::ComputePrecodeValues(
                random_system, random_source.data(), 3u,
                random_parity.data()) ||
            !VerifyValues(
                random_system, random_source.data(), random_parity.data(),
                3u, "random feasible property"))
        {
            std::fprintf(stderr,
                "strict validation: random values violated constraints\n");
            return false;
        }
    }

    std::printf("strict precode validation: PASS\n");
    return true;
}

void ManualRecoveryBlock(
    const wirehair_v2::PrecodeSystem& system,
    const uint8_t* source, const uint8_t* parity,
    uint32_t block_bytes,
    const std::vector<uint32_t>& row_columns,
    uint8_t* block_out)
{
    std::memset(block_out, 0, block_bytes);
    for (const uint32_t col : row_columns)
    {
        gf256_add_mem(block_out,
            ColumnValue(system, source, parity, block_bytes, col),
            (int)block_bytes);
    }
}

bool EqualBlock(
    const uint8_t* a, const uint8_t* b, uint32_t block_bytes)
{
    return std::memcmp(a, b, block_bytes) == 0;
}

bool IsZeroBlock(const uint8_t* a, uint32_t block_bytes)
{
    for (uint32_t i = 0; i < block_bytes; ++i) {
        if (a[i] != 0u) {
            return false;
        }
    }
    return true;
}

bool TestRecoveryBlockEncoding()
{
    const uint32_t K = 1000u;
    const uint32_t bb = 37u;
    wirehair_v2::PrecodeParams params =
        wirehair_v2::MakeCertifiedParams(K, UINT64_C(0x5eed1234));
    params.DenseIdentityCorner = true;

    wirehair_v2::PrecodeSystem system;
    if (!BuildPrecodeSystem(params, system)) {
        std::fprintf(stderr, "recovery encode: build failed\n");
        return false;
    }

    const uint32_t parity_count =
        system.Params.Staircase + system.Params.DenseRows +
        system.Params.HeavyRows;
    std::vector<uint8_t> source((size_t)K * bb);
    std::vector<uint8_t> parity((size_t)parity_count * bb);
    FillRandomBlocks(source.data(), source.size(), UINT64_C(0xbeadfeed));
    wirehair_v2::PrecodeEncodeStats precode_stats;
    if (!wirehair_v2::ComputePrecodeValues(
            system, source.data(), bb, parity.data(), &precode_stats))
    {
        std::fprintf(stderr, "recovery encode: precode values failed\n");
        return false;
    }
    if (!VerifyValues(
            system, source.data(), parity.data(), bb, "recovery encode"))
    {
        return false;
    }

    const wirehair_v2::PeelingCodec codec =
        wirehair_v2::MakePeelingCodec(
            wirehair_v2::PeelStructure::LtM1C32,
            wirehair_v2::PeelSolver::KsBmaxTop16);
    const uint64_t row_seed = UINT64_C(0xdec0de);
    const uint32_t recovery_mix = 5u;
    const std::vector<std::vector<uint32_t> > rows =
        wirehair_v2::GenerateRecoveryMatrixRows(
            codec, K, parity_count, 8u, recovery_mix, row_seed);
    if (rows.size() != 8u) {
        std::fprintf(stderr, "recovery encode: row generation failed\n");
        return false;
    }

    std::vector<uint8_t> got(bb), want(bb);
    for (size_t r = 0; r < rows.size(); ++r)
    {
        uint64_t ops = UINT64_MAX;
        if (!wirehair_v2::ComputeRecoveryBlock(
                system, source.data(), parity.data(), bb, rows[r],
                got.data(), &ops))
        {
            std::fprintf(stderr,
                "recovery encode: ComputeRecoveryBlock failed for row %zu\n",
                r);
            return false;
        }
        ManualRecoveryBlock(
            system, source.data(), parity.data(), bb, rows[r], want.data());
        if (!EqualBlock(got.data(), want.data(), bb)) {
            std::fprintf(stderr,
                "recovery encode: row %zu did not match manual XOR\n", r);
            return false;
        }
        if (ops != rows[r].size()) {
            std::fprintf(stderr,
                "recovery encode: row %zu ops %llu != columns %zu\n",
                r, (unsigned long long)ops, rows[r].size());
            return false;
        }
    }

    wirehair_v2::PrecodeEncoder uninitialized_encoder;
    uint64_t ops = UINT64_MAX;
    if (uninitialized_encoder.HasCompleteSystem() ||
        uninitialized_encoder.Encode(0u, got.data(), &ops) ||
        ops != UINT64_MAX)
    {
        std::fprintf(stderr,
            "recovery encode: uninitialized encoder should fail\n");
        return false;
    }

    wirehair_v2::PrecodeEncoder encoder_state;
    if (!encoder_state.Initialize(
            system, codec, row_seed, recovery_mix, source.data(), bb))
    {
        std::fprintf(stderr,
            "recovery encode: PrecodeEncoder initialize failed\n");
        return false;
    }
    if (!encoder_state.IsInitialized() ||
        !encoder_state.HasCompleteSystem() ||
        !wirehair_v2::ValidatePrecodeSystem(encoder_state.System()) ||
        encoder_state.System().StaircaseRows != system.StaircaseRows ||
        encoder_state.System().DenseRowColumns != system.DenseRowColumns ||
        encoder_state.SourceBlockCount() != K ||
        encoder_state.ParityBlockCount() != parity_count ||
        encoder_state.BlockBytes() != bb ||
        !encoder_state.ParityBlocks())
    {
        std::fprintf(stderr,
            "recovery encode: PrecodeEncoder accessors failed\n");
        return false;
    }
    if (encoder_state.EncodeStats().StaircaseBlockOps !=
            precode_stats.StaircaseBlockOps ||
        encoder_state.EncodeStats().DenseKnownBlockOps !=
            precode_stats.DenseKnownBlockOps ||
        encoder_state.EncodeStats().DenseSolveBlockOps !=
            precode_stats.DenseSolveBlockOps ||
        encoder_state.EncodeStats().HeavyBucketXors !=
            precode_stats.HeavyBucketXors ||
        encoder_state.EncodeStats().HeavyMulAdds !=
            precode_stats.HeavyMulAdds ||
        encoder_state.EncodeStats().HeavySolveBlockOps !=
            precode_stats.HeavySolveBlockOps)
    {
        std::fprintf(stderr,
            "recovery encode: PrecodeEncoder stats mismatch\n");
        return false;
    }
    if (std::memcmp(
            encoder_state.ParityBlocks(), parity.data(),
            (size_t)parity_count * bb) != 0)
    {
        std::fprintf(stderr,
            "recovery encode: PrecodeEncoder parity mismatch\n");
        return false;
    }
    if (encoder_state.Initialize(
            system, codec, row_seed, recovery_mix, nullptr, bb) ||
        !encoder_state.IsInitialized() ||
        !encoder_state.HasCompleteSystem() ||
        !wirehair_v2::ValidatePrecodeSystem(encoder_state.System()) ||
        !encoder_state.ParityBlocks() ||
        encoder_state.BlockBytes() != bb)
    {
        std::fprintf(stderr,
            "recovery encode: failed initialize should preserve state\n");
        return false;
    }
    if (!encoder_state.Initialize(
            system, codec, row_seed, recovery_mix, source.data(), bb))
    {
        std::fprintf(stderr,
            "recovery encode: PrecodeEncoder reinitialize failed\n");
        return false;
    }

    const uint32_t source_id = 17u;
    if (!wirehair_v2::ComputeEncodedBlock(
            system, codec, row_seed, recovery_mix,
            source.data(), nullptr, bb, source_id, got.data(), &ops) ||
        !EqualBlock(
            got.data(), source.data() + (size_t)source_id * bb, bb) ||
        ops != 1u)
    {
        std::fprintf(stderr, "recovery encode: encoded source block failed\n");
        return false;
    }
    if (!encoder_state.Encode(source_id, got.data(), &ops) ||
        !EqualBlock(
            got.data(), source.data() + (size_t)source_id * bb, bb) ||
        ops != 1u)
    {
        std::fprintf(stderr,
            "recovery encode: state encoded source block failed\n");
        return false;
    }

    const uint32_t recovery_index = 4u;
    if (!wirehair_v2::ComputeEncodedBlock(
            system, codec, row_seed, recovery_mix,
            source.data(), parity.data(), bb, K + recovery_index,
            got.data(), &ops))
    {
        std::fprintf(stderr,
            "recovery encode: encoded recovery block failed\n");
        return false;
    }
    ManualRecoveryBlock(
        system, source.data(), parity.data(), bb,
        rows[recovery_index], want.data());
    if (!EqualBlock(got.data(), want.data(), bb) ||
        ops != rows[recovery_index].size())
    {
        std::fprintf(stderr,
            "recovery encode: encoded recovery block mismatch\n");
        return false;
    }
    if (!encoder_state.Encode(K + recovery_index, got.data(), &ops) ||
        !EqualBlock(got.data(), want.data(), bb) ||
        ops != rows[recovery_index].size())
    {
        std::fprintf(stderr,
            "recovery encode: state encoded recovery block mismatch\n");
        return false;
    }

    std::fill(got.begin(), got.end(), uint8_t{0xac});
    if (wirehair_v2::ComputeEncodedBlock(
            system, codec, row_seed, recovery_mix,
            source.data(), nullptr, bb, K + recovery_index,
            got.data(), &ops) ||
        !std::all_of(got.begin(), got.end(),
            [](uint8_t x) { return x == 0xacu; }) ||
        ops != 0u)
    {
        std::fprintf(stderr,
            "recovery encode: recovery block should require parity values\n");
        return false;
    }

    const uint32_t high_block_id = UINT32_MAX;
    const uint32_t high_row_id = high_block_id - K;
    const std::vector<uint32_t> high_row =
        wirehair_v2::GenerateRecoveryMatrixRow(
            codec, K, parity_count, high_row_id, recovery_mix, row_seed);
    ManualRecoveryBlock(
        system, source.data(), parity.data(), bb, high_row, want.data());
    if (!wirehair_v2::ComputeEncodedBlock(
            system, codec, row_seed, recovery_mix,
            source.data(), parity.data(), bb,
            high_block_id, got.data(), &ops) ||
        !EqualBlock(got.data(), want.data(), bb) ||
        ops != high_row.size())
    {
        std::fprintf(stderr,
            "recovery encode: full uint32 block id failed\n");
        return false;
    }
    if (!encoder_state.Encode(high_block_id, got.data(), &ops) ||
        !EqualBlock(got.data(), want.data(), bb) ||
        ops != high_row.size())
    {
        std::fprintf(stderr,
            "recovery encode: state full uint32 block id failed\n");
        return false;
    }

    const std::vector<uint32_t> source_only{0u};
    if (!wirehair_v2::ComputeRecoveryBlock(
            system, source.data(), parity.data(), bb, source_only,
            got.data(), &ops) ||
        !EqualBlock(got.data(), source.data(), bb) || ops != 1u)
    {
        std::fprintf(stderr, "recovery encode: source-only row failed\n");
        return false;
    }

    const std::vector<uint32_t> parity_only{K + system.Params.Staircase};
    if (!wirehair_v2::ComputeRecoveryBlock(
            system, source.data(), parity.data(), bb, parity_only,
            got.data(), &ops) ||
        !EqualBlock(
            got.data(), parity.data() + (size_t)system.Params.Staircase * bb,
            bb) ||
        ops != 1u)
    {
        std::fprintf(stderr, "recovery encode: parity-only row failed\n");
        return false;
    }

    std::fill(got.begin(), got.end(), uint8_t{0xac});
    const std::vector<uint32_t> empty_row;
    if (!wirehair_v2::ComputeRecoveryBlock(
            system, source.data(), parity.data(), bb, empty_row,
            got.data(), &ops) ||
        !IsZeroBlock(got.data(), bb) || ops != 0u)
    {
        std::fprintf(stderr, "recovery encode: empty row failed\n");
        return false;
    }

    const std::vector<uint32_t> bad_row{
        0u,
        K + parity_count
    };
    std::fill(got.begin(), got.end(), uint8_t{0xac});
    if (wirehair_v2::ComputeRecoveryBlock(
            system, source.data(), parity.data(), bb, bad_row,
            got.data(), &ops) ||
        !std::all_of(got.begin(), got.end(),
            [](uint8_t x) { return x == 0xacu; }) ||
        ops != 0u)
    {
        std::fprintf(stderr,
            "recovery encode: invalid row should fail without writes\n");
        return false;
    }

    if (wirehair_v2::ComputeRecoveryBlock(
            system, source.data(), parity.data(), 0u, source_only,
            got.data()))
    {
        std::fprintf(stderr,
            "recovery encode: zero block_bytes should fail\n");
        return false;
    }

    std::printf("recovery block encoding: PASS\n");
    return true;
}

bool TestMessagePrecodeEncoder()
{
    const uint32_t K = 1000u;
    const uint32_t bb = 37u;
    const uint32_t tail = 13u;
    const uint64_t message_bytes = (uint64_t)(K - 1u) * bb + tail;
    std::vector<uint8_t> message((size_t)message_bytes);
    FillRandomBlocks(message.data(), message.size(), UINT64_C(0x1234beef));

    wirehair_v2::MessagePrecodeEncoderOptions options;
    options.DenseIdentityCorner = true;
    options.RecoveryMixCount = wirehair_v2::kCertifiedPacketMixCount;

    wirehair_v2::MessagePrecodeEncoder encoder;
    if (!encoder.Initialize(
            message.data(), message_bytes, bb, nullptr, &options))
    {
        std::fprintf(stderr,
            "message encoder: identity-corner initialize failed\n");
        return false;
    }
    if (!encoder.IsInitialized() ||
        encoder.MessageBytes() != message_bytes ||
        encoder.SourceBlockCount() != K ||
        encoder.BlockBytes() != bb ||
        encoder.Profile().BlockCount != K ||
        encoder.Profile().BlockBytes != bb ||
        encoder.Options().RecoveryMixCount != options.RecoveryMixCount ||
        !encoder.Options().DenseIdentityCorner ||
        !encoder.IntermediateBlocks() ||
        !encoder.BlockEncoder().IsInitialized())
    {
        std::fprintf(stderr, "message encoder: accessors failed\n");
        return false;
    }

    std::vector<uint8_t> got(bb, 0xac), want(bb, 0xbd);
    uint32_t data_bytes = UINT32_MAX;
    uint64_t ops = UINT64_MAX;
    if (!encoder.Encode(0u, got.data(), bb, &data_bytes, &ops) ||
        data_bytes != bb ||
        ops == 0u ||
        std::memcmp(got.data(), message.data(), bb) != 0)
    {
        std::fprintf(stderr, "message encoder: source block failed\n");
        return false;
    }

    std::fill(got.begin(), got.end(), uint8_t{0xac});
    const uint8_t* tail_src =
        message.data() + (size_t)(K - 1u) * bb;
    if (!encoder.Encode(K - 1u, got.data(), tail, &data_bytes, &ops) ||
        data_bytes != tail ||
        ops == 0u ||
        std::memcmp(got.data(), tail_src, tail) != 0 ||
        !std::all_of(got.begin() + tail, got.end(),
            [](uint8_t x) { return x == 0xacu; }))
    {
        std::fprintf(stderr, "message encoder: final partial source failed\n");
        return false;
    }

    std::fill(got.begin(), got.end(), uint8_t{0xac});
    if (!encoder.Encode(K - 1u, got.data(), bb, &data_bytes, &ops) ||
        data_bytes != tail ||
        std::memcmp(got.data(), tail_src, tail) != 0 ||
        !std::all_of(got.begin() + tail, got.end(),
            [](uint8_t x) { return x == 0xacu; }))
    {
        std::fprintf(stderr,
            "message encoder: full-capacity partial source over-wrote tail\n");
        return false;
    }

    std::fill(got.begin(), got.end(), uint8_t{0xac});
    data_bytes = UINT32_MAX;
    ops = UINT64_MAX;
    if (encoder.Encode(K - 1u, got.data(), tail - 1u, &data_bytes, &ops) ||
        data_bytes != UINT32_MAX ||
        ops != UINT64_MAX ||
        !std::all_of(got.begin(), got.end(),
            [](uint8_t x) { return x == 0xacu; }))
    {
        std::fprintf(stderr,
            "message encoder: short partial-source buffer should fail\n");
        return false;
    }

    const uint32_t recovery_id = K + 3u;
    std::fill(got.begin(), got.end(), uint8_t{0xac});
    if (!encoder.Encode(recovery_id, got.data(), bb, &data_bytes, &ops) ||
        data_bytes != bb ||
        ops == 0u)
    {
        std::fprintf(stderr, "message encoder: recovery block failed\n");
        return false;
    }
    uint64_t want_ops = UINT64_MAX;
    if (!encoder.BlockEncoder().Encode(
            recovery_id, want.data(), &want_ops) ||
        want_ops != ops ||
        !EqualBlock(got.data(), want.data(), bb))
    {
        std::fprintf(stderr,
            "message encoder: recovery block mismatch\n");
        return false;
    }

    const wirehair_v2::PrecodeEncoder& blocks = encoder.BlockEncoder();
    const wirehair_v2::PrecodeSystem& encoder_system = blocks.System();
    wirehair_v2::PacketRowConfig packet_config;
    packet_config.PeelSeed = (uint32_t)blocks.RecoveryRowSeed();
    packet_config.MixCount = blocks.RecoveryMixCount();
    if (!blocks.IntermediateBlocks() || blocks.HasCompleteSystem() ||
        !encoder_system.StaircaseRows.empty() ||
        !encoder_system.DenseRowColumns.empty() ||
        encoder_system.Params.BlockCount != K ||
        encoder_system.Params.Staircase != encoder.Profile().V2StaircaseCount ||
        encoder_system.Params.DenseRows != encoder.Profile().V2DenseRowCount ||
        encoder_system.Params.HeavyRows != encoder.Profile().V2HeavyRowCount)
    {
        std::fprintf(stderr,
            "message encoder: solved state retained or lost system data\n");
        return false;
    }
    const uint32_t reordered_ids[] = {
        K + 9u, 2u, K + 1u, K - 1u, K + 15u, 0u, K + 4u
    };
    for (uint32_t packet_id : reordered_ids)
    {
        std::fill(got.begin(), got.end(), uint8_t{0xac});
        std::fill(want.begin(), want.end(), uint8_t{0xbd});
        data_bytes = UINT32_MAX;
        if (!encoder.Encode(
                packet_id, got.data(), bb, &data_bytes, &ops))
        {
            std::fprintf(stderr,
                "message encoder: reordered packet %u failed\n", packet_id);
            return false;
        }
        if (packet_id < K)
        {
            const uint32_t expected_bytes =
                packet_id + 1u < K ? bb : tail;
            if (data_bytes != expected_bytes ||
                std::memcmp(
                    got.data(),
                    message.data() + (size_t)packet_id * bb,
                    expected_bytes) != 0)
            {
                std::fprintf(stderr,
                    "message encoder: reordered source oracle mismatch\n");
                return false;
            }
        }
        else
        {
            if (!ManualPacketBlock(
                encoder_system,
                packet_config,
                blocks.IntermediateBlocks(),
                bb,
                packet_id,
                want.data()))
            {
                std::fprintf(stderr,
                    "message encoder: packet oracle evaluation failed\n");
                return false;
            }
            if (data_bytes != bb || !EqualBlock(got.data(), want.data(), bb))
            {
                std::fprintf(stderr,
                    "message encoder: reordered repair oracle mismatch\n");
                return false;
            }
        }
    }

    const uint32_t max_packet_id = UINT32_MAX;
    if (!ManualPacketBlock(
            encoder_system, packet_config, blocks.IntermediateBlocks(), bb,
            max_packet_id, want.data()))
    {
        std::fprintf(stderr,
            "message encoder: maximum packet-id oracle failed\n");
        return false;
    }
    if (!encoder.Encode(
            max_packet_id, got.data(), bb, &data_bytes, &ops) ||
        data_bytes != bb || !EqualBlock(got.data(), want.data(), bb))
    {
        std::fprintf(stderr,
            "message encoder: maximum packet-id oracle mismatch\n");
        return false;
    }

    std::fill(got.begin(), got.end(), uint8_t{0xac});
    data_bytes = UINT32_MAX;
    ops = UINT64_MAX;
    if (encoder.Encode(recovery_id, got.data(), bb - 1u, &data_bytes, &ops) ||
        data_bytes != UINT32_MAX ||
        ops != UINT64_MAX ||
        !std::all_of(got.begin(), got.end(),
            [](uint8_t x) { return x == 0xacu; }))
    {
        std::fprintf(stderr,
            "message encoder: short recovery buffer should fail\n");
        return false;
    }

    wirehair_v2::SeedProfile mismatch =
        wirehair_v2::SelectSeedProfile(K, bb);
    ++mismatch.BlockCount;
    if (encoder.Initialize(
            message.data(), message_bytes, bb, &mismatch, &options) ||
        !encoder.IsInitialized() ||
        !encoder.IntermediateBlocks() ||
        encoder.BlockEncoder().HasCompleteSystem() ||
        !encoder.BlockEncoder().System().StaircaseRows.empty() ||
        !encoder.BlockEncoder().System().DenseRowColumns.empty() ||
        encoder.SourceBlockCount() != K ||
        !encoder.BlockEncoder().IsInitialized())
    {
        std::fprintf(stderr,
            "message encoder: failed reinitialize should preserve state\n");
        return false;
    }

    if (encoder.Initialize(nullptr, message_bytes, bb, nullptr, &options) ||
        !encoder.IsInitialized())
    {
        std::fprintf(stderr,
            "message encoder: null message should fail\n");
        return false;
    }
    const uint8_t dummy = 0u;
    if (encoder.Initialize(
            &dummy,
            UINT64_C(0x100000000),
            UINT32_C(0x80000000),
            nullptr,
            &options) ||
        !encoder.IsInitialized())
    {
        std::fprintf(stderr,
            "message encoder: oversized block should fail before allocation\n");
        return false;
    }

    std::printf("message precode encoder: PASS\n");
    return true;
}

bool TestBorrowedMessageLifetime()
{
    const uint32_t K = 64u;
    const uint32_t bb = 67u;
    const uint32_t final_sizes[] = { bb, 13u };
    const uint32_t packet_ids[] = { 0u, K - 1u, K + 7u };

    for (uint32_t final_bytes : final_sizes)
    {
        const uint64_t message_bytes =
            (uint64_t)(K - 1u) * bb + final_bytes;
        std::vector<uint8_t> message((size_t)message_bytes);
        FillRandomBlocks(
            message.data(), message.size(),
            UINT64_C(0xb0770ed) ^ final_bytes);

        wirehair_v2::MessagePrecodeEncoder encoder;
        if (encoder.InitializeResult(
                message.data(), message_bytes, bb) != Wirehair_Success)
        {
            std::fprintf(stderr,
                "borrowed lifetime: initialize failed for final=%u\n",
                final_bytes);
            return false;
        }

        std::vector<std::vector<uint8_t> > expected;
        std::vector<uint32_t> expected_bytes;
        for (uint32_t packet_id : packet_ids)
        {
            expected.push_back(std::vector<uint8_t>(bb, 0xa5u));
            uint32_t data_bytes = UINT32_MAX;
            if (encoder.EncodeResult(
                    packet_id, expected.back().data(), bb, &data_bytes) !=
                    Wirehair_Success)
            {
                std::fprintf(stderr,
                    "borrowed lifetime: initial encode failed id=%u\n",
                    packet_id);
                return false;
            }
            expected_bytes.push_back(data_bytes);
        }

        if (final_bytes != bb)
        {
            std::vector<uint8_t> padded(bb, 0xa5u);
            uint64_t ops = UINT64_MAX;
            if (!encoder.BlockEncoder().Encode(
                    K - 1u, padded.data(), &ops) ||
                std::memcmp(
                    padded.data(),
                    message.data() + (size_t)(K - 1u) * bb,
                    final_bytes) != 0 ||
                !std::all_of(
                    padded.begin() + final_bytes, padded.end(),
                    [](uint8_t value) { return value == 0u; }))
            {
                std::fprintf(stderr,
                    "borrowed lifetime: partial final padding was not zero\n");
                return false;
            }
        }

        // Release the caller allocation completely.  ASan will catch any
        // post-return read through a stale borrowed pointer; byte comparisons
        // also cover builds without sanitizers.
        std::vector<uint8_t>().swap(message);
        for (size_t i = 0; i <
                sizeof(packet_ids) / sizeof(packet_ids[0]); ++i)
        {
            std::vector<uint8_t> actual(bb, 0x5au);
            uint32_t data_bytes = UINT32_MAX;
            if (encoder.EncodeResult(
                    packet_ids[i], actual.data(), bb, &data_bytes) !=
                    Wirehair_Success ||
                data_bytes != expected_bytes[i] ||
                std::memcmp(
                    actual.data(), expected[i].data(), expected_bytes[i]) != 0 ||
                !std::all_of(
                    actual.begin() + expected_bytes[i], actual.end(),
                    [](uint8_t value) { return value == 0x5au; }))
            {
                std::fprintf(stderr,
                    "borrowed lifetime: encoder retained caller storage "
                    "for final=%u id=%u\n",
                    final_bytes, packet_ids[i]);
                return false;
            }
        }
    }

    std::printf("borrowed message lifetime: PASS\n");
    return true;
}

bool TestSystematicSourceCache()
{
    const uint32_t K = 64u;
    const uint32_t bb = 67u;
    const uint32_t tail = 13u;
    const uint64_t message_bytes = (uint64_t)(K - 1u) * bb + tail;
    std::vector<uint8_t> message((size_t)message_bytes);
    FillRandomBlocks(message.data(), message.size(), UINT64_C(0xca5ced));
    const std::vector<uint8_t> original = message;

    wirehair_v2::MessagePrecodeEncoder uncached;
    if (uncached.InitializeResult(
            message.data(), message_bytes, bb) != Wirehair_Success ||
        uncached.HasSystematicSourceCache() ||
        uncached.SystematicSourceCacheBytes() != 0u)
    {
        std::fprintf(stderr, "systematic cache: default path retained source\n");
        return false;
    }

    wirehair_v2::MessagePrecodeEncoderOptions cached_options;
    cached_options.CacheSystematicSource = true;
    wirehair_v2::MessagePrecodeEncoder cached;
    const wirehair_v2::SeedProfile selected_profile = uncached.Profile();
    if (cached.InitializeResult(
            message.data(), message_bytes, bb,
            &selected_profile, &cached_options) != Wirehair_Success ||
        !cached.Options().CacheSystematicSource ||
        !cached.HasSystematicSourceCache() ||
        cached.SystematicSourceCacheBytes() != message_bytes)
    {
        std::fprintf(stderr,
            "systematic cache: selected-profile opt-in failed\n");
        return false;
    }

    // A local storage policy must not alter any serialized matrix or packet
    // contract field selected by the uncached encoder.
    const wirehair_v2::SeedProfile& cached_profile = cached.Profile();
    if (cached_profile.BlockCount != selected_profile.BlockCount ||
        cached_profile.BlockBytes != selected_profile.BlockBytes ||
        cached_profile.DenseCount != selected_profile.DenseCount ||
        cached_profile.V2SeedAttempt != selected_profile.V2SeedAttempt ||
        cached_profile.V2PrecodeContractVersion !=
            selected_profile.V2PrecodeContractVersion ||
        cached_profile.V2PacketRowContractVersion !=
            selected_profile.V2PacketRowContractVersion ||
        cached_profile.V2StaircaseCount !=
            selected_profile.V2StaircaseCount ||
        cached_profile.V2DenseRowCount != selected_profile.V2DenseRowCount ||
        cached_profile.V2HeavyRowCount != selected_profile.V2HeavyRowCount ||
        cached_profile.V2CompletionField !=
            selected_profile.V2CompletionField ||
        cached_profile.V2SourceHits != selected_profile.V2SourceHits ||
        cached_profile.V2PrecodeSeed != selected_profile.V2PrecodeSeed ||
        cached_profile.V2PacketPeelSeed !=
            selected_profile.V2PacketPeelSeed ||
        cached_profile.V2RecoveryMixCount !=
            selected_profile.V2RecoveryMixCount ||
        cached_profile.V2DenseIdentityCorner !=
            selected_profile.V2DenseIdentityCorner ||
        cached_profile.V2PrecodeSeedSalt !=
            selected_profile.V2PrecodeSeedSalt ||
        cached_profile.V2RecoveryRowSeedSalt !=
            selected_profile.V2RecoveryRowSeedSalt)
    {
        std::fprintf(stderr,
            "systematic cache: local option changed selected profile\n");
        return false;
    }

    const uint32_t packet_ids[] = { 0u, K - 1u, K + 7u };
    std::vector<std::vector<uint8_t> > expected;
    std::vector<uint32_t> expected_bytes;
    for (const uint32_t packet_id : packet_ids)
    {
        expected.push_back(std::vector<uint8_t>(bb, 0xa5u));
        uint32_t bytes = UINT32_MAX;
        uint64_t uncached_ops = UINT64_MAX;
        if (uncached.EncodeResult(
                packet_id, expected.back().data(), bb,
                &bytes, &uncached_ops) != Wirehair_Success ||
            uncached_ops == 0u)
        {
            std::fprintf(stderr,
                "systematic cache: uncached reference encode failed\n");
            return false;
        }
        expected_bytes.push_back(bytes);

        std::vector<uint8_t> actual(bb, 0xa5u);
        uint32_t actual_bytes = UINT32_MAX;
        uint64_t cached_ops = UINT64_MAX;
        if (cached.EncodeResult(
                packet_id, actual.data(), bb,
                &actual_bytes, &cached_ops) != Wirehair_Success ||
            actual_bytes != bytes || actual != expected.back() ||
            (packet_id < K ? cached_ops != 0u : cached_ops != uncached_ops))
        {
            std::fprintf(stderr,
                "systematic cache: attached wire mismatch id=%u\n",
                packet_id);
            return false;
        }
    }

    // The cache is an owned byte-for-byte snapshot, not another borrowed
    // view of the caller allocation.
    std::fill(message.begin(), message.end(), uint8_t{0x5a});
    for (size_t i = 0u; i < 2u; ++i)
    {
        std::vector<uint8_t> actual(bb, 0xa5u);
        uint32_t bytes = UINT32_MAX;
        uint64_t ops = UINT64_MAX;
        if (cached.EncodeResult(
                packet_ids[i], actual.data(), bb, &bytes, &ops) !=
                Wirehair_Success ||
            bytes != expected_bytes[i] || actual != expected[i] || ops != 0u)
        {
            std::fprintf(stderr,
                "systematic cache: caller mutation changed cached bytes\n");
            return false;
        }
    }

    cached.ReleaseSystematicSourceCache();
    cached.ReleaseSystematicSourceCache();
    if (cached.HasSystematicSourceCache() ||
        cached.SystematicSourceCacheBytes() != 0u ||
        cached.Options().CacheSystematicSource)
    {
        std::fprintf(stderr, "systematic cache: release was not idempotent\n");
        return false;
    }
    for (size_t i = 0u; i < 3u; ++i)
    {
        std::vector<uint8_t> actual(bb, 0xa5u);
        uint32_t bytes = UINT32_MAX;
        uint64_t ops = UINT64_MAX;
        if (cached.EncodeResult(
                packet_ids[i], actual.data(), bb, &bytes, &ops) !=
                Wirehair_Success ||
            bytes != expected_bytes[i] || actual != expected[i] || ops == 0u)
        {
            std::fprintf(stderr,
                "systematic cache: detached fallback mismatch id=%u\n",
                packet_ids[i]);
            return false;
        }
    }

    // Failed cached initialization is transactional, including cache state.
    message = original;
    if (cached.InitializeResult(
            message.data(), message_bytes, bb,
            &selected_profile, &cached_options) != Wirehair_Success)
    {
        std::fprintf(stderr, "systematic cache: reattach initialize failed\n");
        return false;
    }
    const uint8_t* const working_intermediate = cached.IntermediateBlocks();
    wirehair_v2::SetAllocationFailureCountdownForTesting(0);
    const WirehairResult cache_oom = cached.InitializeResult(
        message.data(), message_bytes, bb,
        &selected_profile, &cached_options);
    wirehair_v2::SetAllocationFailureCountdownForTesting(-1);
    if (cache_oom != Wirehair_OOM ||
        cached.IntermediateBlocks() != working_intermediate ||
        !cached.HasSystematicSourceCache() ||
        cached.SystematicSourceCacheBytes() != message_bytes)
    {
        std::fprintf(stderr,
            "systematic cache: OOM did not preserve active encoder\n");
        return false;
    }

    // A successful default reinitialization returns to x8rs.4's zero-copy
    // source-lifetime policy and releases the prior opt-in cache.
    if (cached.InitializeResult(
            message.data(), message_bytes, bb,
            &selected_profile, nullptr) != Wirehair_Success ||
        cached.HasSystematicSourceCache() ||
        cached.Options().CacheSystematicSource)
    {
        std::fprintf(stderr,
            "systematic cache: default reinitialize retained cache\n");
        return false;
    }

    wirehair_v2::Codec facade;
    if (facade.ReleasePrecodeEncoderSystematicCache() !=
            Wirehair_InvalidInput ||
        facade.InitializePrecodeEncoder(
            message.data(), message_bytes, bb,
            &selected_profile, &cached_options) != Wirehair_Success ||
        facade.ReleasePrecodeEncoderSystematicCache() != Wirehair_Success ||
        facade.ReleasePrecodeEncoderSystematicCache() != Wirehair_Success)
    {
        std::fprintf(stderr, "systematic cache: facade detach failed\n");
        return false;
    }
    std::vector<uint8_t> facade_block(bb, 0xa5u);
    uint32_t facade_bytes = UINT32_MAX;
    if (facade.Encode(
            K - 1u, facade_block.data(), bb, &facade_bytes) !=
            Wirehair_Success ||
        facade_bytes != tail ||
        std::memcmp(
            facade_block.data(), original.data() + (size_t)(K - 1u) * bb,
            tail) != 0)
    {
        std::fprintf(stderr,
            "systematic cache: facade fallback changed systematic bytes\n");
        return false;
    }

    std::printf("systematic source cache: PASS\n");
    return true;
}

bool TestTypedFailuresAndAllocationContainment()
{
    const uint32_t K = 16u;
    const uint32_t bb = 19u;
    const uint64_t message_bytes = (uint64_t)K * bb;
    std::vector<uint8_t> message((size_t)message_bytes);
    FillRandomBlocks(message.data(), message.size(), UINT64_C(0xa110ca7e));

    wirehair_v2::MessagePrecodeEncoderOptions identity;
    identity.DenseIdentityCorner = true;
    wirehair_v2::MessagePrecodeEncoder encoder;
    if (encoder.InitializeResult(
            message.data(), message_bytes, bb, nullptr, &identity) !=
            Wirehair_Success)
    {
        std::fprintf(stderr, "typed failures: initial encoder failed\n");
        return false;
    }

    const wirehair_v2::SeedProfile working_profile = encoder.Profile();
    const uint8_t* working_source = encoder.IntermediateBlocks();
    const uint8_t dummy = 0u;
    if (encoder.InitializeResult(
            &dummy, UINT64_C(0x100000000), UINT32_C(0x80000000),
            nullptr, &identity) != Wirehair_InvalidInput ||
        !encoder.IsInitialized() ||
        encoder.IntermediateBlocks() != working_source ||
        encoder.Profile().BlockCount != working_profile.BlockCount)
    {
        std::fprintf(stderr,
            "typed failures: invalid input classification/state failed\n");
        return false;
    }

    wirehair_v2::MessagePrecodeEncoderOptions certified;
    if (encoder.InitializeResult(
            message.data(), message_bytes, bb, nullptr, &certified) !=
            Wirehair_Success || !encoder.IsInitialized() ||
        encoder.Options().DenseIdentityCorner)
    {
        std::fprintf(stderr,
            "typed failures: certified global solve failed\n");
        return false;
    }
    working_source = encoder.IntermediateBlocks();

    // The solved encoder no longer copies the nested precode row graph, so
    // only the three real guarded allocations before ownership transfer remain.
    for (int64_t failure = 0; failure < 3; ++failure)
    {
        wirehair_v2::SetAllocationFailureCountdownForTesting(failure);
        const WirehairResult init_oom = encoder.InitializeResult(
            message.data(), message_bytes, bb, nullptr, &identity);
        wirehair_v2::SetAllocationFailureCountdownForTesting(-1);
        if (init_oom != Wirehair_OOM || !encoder.IsInitialized() ||
            encoder.IntermediateBlocks() != working_source)
        {
            std::fprintf(stderr,
                "typed failures: init OOM containment failed at %lld\n",
                (long long)failure);
            return false;
        }
    }

    // A partial input owns exactly one padded tail during the synchronous
    // solve.  Failure of that allocation must preserve the prior encoder.
    wirehair_v2::SetAllocationFailureCountdownForTesting(0);
    const WirehairResult tail_oom = encoder.InitializeResult(
        message.data(), message_bytes - 1u, bb, nullptr, &identity);
    wirehair_v2::SetAllocationFailureCountdownForTesting(-1);
    if (tail_oom != Wirehair_OOM || !encoder.IsInitialized() ||
        encoder.IntermediateBlocks() != working_source)
    {
        std::fprintf(stderr,
            "typed failures: padded-tail OOM containment failed\n");
        return false;
    }

    std::vector<uint8_t> output(bb, 0xa5u);
    const std::vector<uint8_t> before = output;
    uint32_t data_bytes = UINT32_MAX;
    uint64_t ops = UINT64_MAX;
    wirehair_v2::SetAllocationFailureCountdownForTesting(0);
    const WirehairResult encode_oom = encoder.EncodeResult(
        K, output.data(), bb, &data_bytes, &ops);
    wirehair_v2::SetAllocationFailureCountdownForTesting(-1);
    if (encode_oom != Wirehair_OOM || output != before ||
        data_bytes != UINT32_MAX || ops != UINT64_MAX)
    {
        std::fprintf(stderr,
            "typed failures: recovery Encode OOM modified output\n");
        return false;
    }
    if (encoder.EncodeResult(
            K, output.data(), bb, &data_bytes, &ops) != Wirehair_Success ||
        data_bytes != bb || ops == 0u)
    {
        std::fprintf(stderr, "typed failures: post-OOM encode failed\n");
        return false;
    }

    wirehair_v2::Codec facade;
    if (facade.InitializePrecodeEncoder(
            message.data(), message_bytes, bb, nullptr, &identity) !=
            Wirehair_Success)
    {
        std::fprintf(stderr, "typed failures: facade initial encode failed\n");
        return false;
    }
    if (facade.InitializePrecodeEncoder(
            &dummy, UINT64_C(0x100000000), UINT32_C(0x80000000),
            nullptr, &identity) != Wirehair_InvalidInput ||
        facade.Profile().BlockCount != K ||
        facade.InitializePrecodeEncoder(
            message.data(), message_bytes, bb, nullptr, &certified) !=
            Wirehair_Success ||
        facade.Profile().BlockCount != K)
    {
        std::fprintf(stderr,
            "typed failures: facade invalid/singular classification failed\n");
        return false;
    }
    wirehair_v2::SetAllocationFailureCountdownForTesting(0);
    const WirehairResult facade_oom = facade.InitializePrecodeEncoder(
        message.data(), message_bytes, bb, nullptr, &identity);
    wirehair_v2::SetAllocationFailureCountdownForTesting(-1);
    if (facade_oom != Wirehair_OOM ||
        facade.Profile().BlockCount != K)
    {
        std::fprintf(stderr,
            "typed failures: facade OOM classification/state failed\n");
        return false;
    }

    std::fill(output.begin(), output.end(), uint8_t{0x5a});
    const std::vector<uint8_t> facade_before = output;
    data_bytes = UINT32_MAX;
    wirehair_v2::SetAllocationFailureCountdownForTesting(0);
    const WirehairResult facade_encode_oom = facade.Encode(
        K, output.data(), bb, &data_bytes);
    wirehair_v2::SetAllocationFailureCountdownForTesting(-1);
    if (facade_encode_oom != Wirehair_OOM || output != facade_before ||
        data_bytes != UINT32_MAX)
    {
        std::fprintf(stderr,
            "typed failures: facade Encode OOM modified output\n");
        return false;
    }

    std::printf("typed failures and OOM containment: PASS\n");
    return true;
}

} // namespace

int main(int argc, char** argv)
{
    gf256_init();

    if (argc == 2 && std::strcmp(argv[1], "--oom-probe") == 0)
    {
        const uint8_t dummy = 0u;
        wirehair_v2::MessagePrecodeEncoderOptions options;
        options.DenseIdentityCorner = true;
        wirehair_v2::Codec codec;
        const WirehairResult result = codec.InitializePrecodeEncoder(
            &dummy,
            UINT64_C(512) * 1024u * 1024u,
            UINT32_C(8) * 1024u * 1024u,
            nullptr,
            &options);
        if (result != Wirehair_OOM) {
            std::fprintf(stderr,
                "constrained precode init returned %d, expected OOM\n",
                (int)result);
            return 1;
        }
        std::printf("constrained precode init: Wirehair_OOM\n");
        return 0;
    }

    // Feasibility trials tunable so the sanitizer build stays quick
    uint32_t trials = 2000u;
    if (argc > 2 || (argc == 2 && !ParsePositiveU32(argv[1], trials))) {
        std::fprintf(stderr,
            "usage: %s [positive feasibility-trials]\n", argv[0]);
        return 1;
    }

    bool ok = true;
    ok = TestCorrectnessAsBuilt() && ok;
    ok = TestCorrectnessDoctored() && ok;
    ok = TestMalformedDenseCorner() && ok;
    ok = TestStrictSystemValidation() && ok;
    ok = TestCostModel() && ok;
    ok = TestHeavyResidueDispatch() && ok;
    ok = TestMixedCornerRank() && ok;
    ok = TestMixedCompletion() && ok;
    ok = TestRecoveryBlockEncoding() && ok;
    ok = TestMessagePrecodeEncoder() && ok;
    ok = TestBorrowedMessageLifetime() && ok;
    ok = TestSystematicSourceCache() && ok;
    ok = TestTypedFailuresAndAllocationContainment() && ok;
    ok = TestFeasibility(trials) && ok;

    if (!ok) {
        std::fprintf(stderr, "precode_encode_test: FAIL\n");
        return 1;
    }
    std::printf("\nprecode_encode_test: PASS\n");
    return 0;
}
