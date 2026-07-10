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

    // Heavy rows: GF(256) weighted sum over ALL L columns
    for (uint32_t r = 0; r < H; ++r)
    {
        std::memset(acc.data(), 0, block_bytes);
        for (uint32_t c = 0; c < L; ++c)
        {
            const uint8_t y = wirehair_v2::HeavyCoefficient(r, c, H);
            const uint8_t* v =
                ColumnValue(system, source, parity, block_bytes, c);
            if (y == 1u) {
                gf256_add_mem(acc.data(), v, (int)block_bytes);
            }
            else {
                gf256_muladd_mem(acc.data(), y, v, (int)block_bytes);
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
        stats.HeavySolveBlockOps == 0u;
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
    if (uninitialized_encoder.Encode(0u, got.data(), &ops) ||
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
    options.RecoveryMixCount = 5u;

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
        !encoder.Options().UseWirehairRowDistribution ||
        !encoder.SourceBlocks() ||
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
        ops != 1u ||
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
        ops != 1u ||
        std::memcmp(got.data(), tail_src, tail) != 0 ||
        !std::all_of(got.begin() + tail, got.end(),
            [](uint8_t x) { return x == 0xacu; }))
    {
        std::fprintf(stderr, "message encoder: final partial source failed\n");
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
    const uint32_t parity_count = blocks.ParityBlockCount();
    wirehair_v2::PeelingCodec recovery_codec =
        encoder.Profile().Policy.Codec;
    recovery_codec.UseWirehairRowDistribution =
        encoder.Options().UseWirehairRowDistribution;
    const std::vector<std::vector<uint32_t> > oracle_rows =
        wirehair_v2::GenerateRecoveryMatrixRows(
            recovery_codec,
            K,
            parity_count,
            16u,
            blocks.RecoveryMixCount(),
            blocks.RecoveryRowSeed());
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
            const uint32_t row_index = packet_id - K;
            ManualRecoveryBlock(
                encoder_system,
                encoder.SourceBlocks(),
                blocks.ParityBlocks(),
                bb,
                oracle_rows[row_index],
                want.data());
            if (data_bytes != bb || !EqualBlock(got.data(), want.data(), bb))
            {
                std::fprintf(stderr,
                    "message encoder: reordered repair oracle mismatch\n");
                return false;
            }
        }
    }

    const uint32_t max_packet_id = UINT32_MAX;
    const std::vector<uint32_t> max_row =
        wirehair_v2::GenerateRecoveryMatrixRow(
            recovery_codec,
            K,
            parity_count,
            max_packet_id - K,
            blocks.RecoveryMixCount(),
            blocks.RecoveryRowSeed());
    ManualRecoveryBlock(
        encoder_system, encoder.SourceBlocks(), blocks.ParityBlocks(), bb,
        max_row, want.data());
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
        !encoder.SourceBlocks() ||
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
    const uint8_t* const working_source = encoder.SourceBlocks();
    const uint8_t dummy = 0u;
    if (encoder.InitializeResult(
            &dummy, UINT64_C(0x100000000), UINT32_C(0x80000000),
            nullptr, &identity) != Wirehair_InvalidInput ||
        !encoder.IsInitialized() || encoder.SourceBlocks() != working_source ||
        encoder.Profile().BlockCount != working_profile.BlockCount)
    {
        std::fprintf(stderr,
            "typed failures: invalid input classification/state failed\n");
        return false;
    }

    wirehair_v2::MessagePrecodeEncoderOptions certified;
    certified.DenseIdentityCorner = false;
    if (encoder.InitializeResult(
            message.data(), message_bytes, bb, nullptr, &certified) !=
            Wirehair_BadDenseSeed || !encoder.IsInitialized() ||
        encoder.SourceBlocks() != working_source)
    {
        std::fprintf(stderr,
            "typed failures: singular dense classification/state failed\n");
        return false;
    }

    for (int64_t failure = 0; failure < 5; ++failure)
    {
        wirehair_v2::SetAllocationFailureCountdownForTesting(failure);
        const WirehairResult init_oom = encoder.InitializeResult(
            message.data(), message_bytes, bb, nullptr, &identity);
        wirehair_v2::SetAllocationFailureCountdownForTesting(-1);
        if (init_oom != Wirehair_OOM || !encoder.IsInitialized() ||
            encoder.SourceBlocks() != working_source)
        {
            std::fprintf(stderr,
                "typed failures: init OOM containment failed at %lld\n",
                (long long)failure);
            return false;
        }
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
            Wirehair_BadDenseSeed ||
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
    ok = TestRecoveryBlockEncoding() && ok;
    ok = TestMessagePrecodeEncoder() && ok;
    ok = TestTypedFailuresAndAllocationContainment() && ok;
    ok = TestFeasibility(trials) && ok;

    if (!ok) {
        std::fprintf(stderr, "precode_encode_test: FAIL\n");
        return 1;
    }
    std::printf("\nprecode_encode_test: PASS\n");
    return 0;
}
