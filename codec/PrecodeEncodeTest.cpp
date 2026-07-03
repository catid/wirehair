#include "WirehairV2PrecodeEncode.h"

#include "../WirehairTools.h"
#include "../gf256.h"

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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
                masks[r] |= UINT64_C(1) << (col - dense_base);
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

/// Rewrite the dense rows' corner: strip all dense columns, then set row
/// r's dense-column membership from masks[r].  The known (source +
/// staircase) parts are untouched; consecutive rows no longer differ in
/// exactly 2 columns, which also exercises the encoder's generic
/// difference path.
void ForceCorner(
    wirehair_v2::PrecodeSystem& system, const uint64_t* masks)
{
    const uint32_t dense_base =
        system.Params.BlockCount + system.Params.Staircase;
    for (uint32_t r = 0; r < system.Params.DenseRows; ++r)
    {
        std::vector<uint32_t>& row = system.DenseRowColumns[r];
        row.erase(
            std::remove_if(row.begin(), row.end(),
                [&](uint32_t col) { return col >= dense_base; }),
            row.end());
        for (uint32_t j = 0; j < system.Params.DenseRows; ++j)
        {
            if (0u != (masks[r] >> j & 1u)) {
                row.push_back(dense_base + j); // ascending j keeps it sorted
            }
        }
    }
}

/// Doctor the dense corner invertible: identity for the trivial-solve
/// case, otherwise a random invertible corner so the GF(2) Gauss-Jordan
/// block eliminations run at every K
void DoctorCorner(wirehair_v2::PrecodeSystem& system, uint32_t seed)
{
    const uint32_t D2 = system.Params.DenseRows;
    std::vector<uint64_t> masks(D2);
    if (seed % 3u == 0u)
    {
        for (uint32_t r = 0; r < D2; ++r) {
            masks[r] = UINT64_C(1) << r;
        }
        ForceCorner(system, masks.data());
        return;
    }
    wirehair::PCGRandom prng;
    prng.Seed(seed, 0xc042e2u);
    do
    {
        for (uint32_t r = 0; r < D2; ++r) {
            masks[r] = prng.Next() & ((UINT64_C(1) << D2) - 1u);
        }
        ForceCorner(system, masks.data());
    } while (DenseCornerRank(system) != D2);
}

bool TestCorrectnessDoctored()
{
    const uint32_t Ks[] = {2u, 16u, 1000u, 3200u};
    const uint32_t bbs[] = {1u, 13u, 1280u};

    for (const uint32_t K : Ks)
    {
        for (const uint32_t bb : bbs)
        {
            for (uint32_t seed = 100; seed < 103u; ++seed)
            {
                wirehair_v2::PrecodeSystem system;
                if (!BuildPrecodeSystem(
                        wirehair_v2::MakeCertifiedParams(K, seed), system))
                {
                    std::fprintf(stderr, "K=%u: build failed\n", K);
                    return false;
                }
                DoctorCorner(system, seed);
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
    const uint32_t Ks[] = {16u, 1000u, 3200u};
    bool ok = true;

    std::printf("\ncost check vs model 2K + S - 1 + ceil(span/2) + "
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

        const uint64_t staircase_model = 2ull * K + S - 1u;
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
    return true;
}

} // namespace

int main(int argc, char** argv)
{
    gf256_init();

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
    ok = TestCostModel() && ok;
    ok = TestFeasibility(trials) && ok;

    if (!ok) {
        std::fprintf(stderr, "precode_encode_test: FAIL\n");
        return 1;
    }
    std::printf("\nprecode_encode_test: PASS\n");
    return 0;
}
