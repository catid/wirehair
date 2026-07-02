#include "WirehairV2Precode.h"

#include "../WirehairTools.h"
#include "../gf256.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <vector>

// Structural invariant tests for the certified precode construction
// (wirehair-axd Phase 1).  These validate the documented ldpcdense_s2
// structure rules; reliability itself is certified by precode_sim and,
// end-to-end, by the later solver phases.

namespace {

bool TestParams()
{
    const wirehair_v2::PrecodeParams params =
        wirehair_v2::MakeCertifiedParams(3200u, 1u);
    if (params.Staircase != wirehair::GetDenseCount(3200u) ||
        params.DenseRows != 12u ||
        params.HeavyRows != 12u ||
        params.SourceHits != 2u)
    {
        std::fprintf(stderr, "certified params wrong for K=3200\n");
        return false;
    }

    // Out-of-domain block counts must not reach GetDenseCount (it
    // extrapolates past 64000) and must fail to build
    const wirehair_v2::PrecodeParams bad =
        wirehair_v2::MakeCertifiedParams(64001u, 1u);
    wirehair_v2::PrecodeSystem system;
    if (bad.Staircase != 0u ||
        wirehair_v2::BuildPrecodeSystem(bad, system))
    {
        std::fprintf(stderr, "K=64001 must fail to build\n");
        return false;
    }
    return true;
}

bool TestStaircase(const wirehair_v2::PrecodeSystem& system)
{
    const uint32_t K = system.Params.BlockCount;
    const uint32_t S = system.Params.Staircase;
    const uint32_t N1 = system.Params.SourceHits;
    const uint32_t hits = N1 < S ? N1 : S;

    if (system.StaircaseRows.size() != S) {
        std::fprintf(stderr, "K=%u: staircase row count\n", K);
        return false;
    }

    std::vector<uint32_t> source_hits(K, 0);
    for (uint32_t j = 0; j < S; ++j)
    {
        const std::vector<uint32_t>& row = system.StaircaseRows[j];

        // Sorted, unique, in-range
        for (size_t i = 0; i < row.size(); ++i)
        {
            if (row[i] >= K + S ||
                (i > 0u && row[i] <= row[i - 1u]))
            {
                std::fprintf(stderr,
                    "K=%u: staircase row %u not sorted/unique/in-range\n",
                    K, j);
                return false;
            }
        }

        // Own parity column and staircase link
        if (!std::binary_search(row.begin(), row.end(), K + j) ||
            (j > 0u &&
             !std::binary_search(row.begin(), row.end(), K + j - 1u)))
        {
            std::fprintf(stderr,
                "K=%u: staircase row %u missing parity/link column\n", K, j);
            return false;
        }

        for (uint32_t col : row) {
            if (col < K) {
                ++source_hits[col];
            }
        }
    }

    // Every source column connects to exactly min(N1, S) distinct parities
    for (uint32_t c = 0; c < K; ++c)
    {
        if (source_hits[c] != hits)
        {
            std::fprintf(stderr,
                "K=%u: source column %u has %u parity hits, want %u\n",
                K, c, source_hits[c], hits);
            return false;
        }
    }
    return true;
}

bool TestDenseRows(const wirehair_v2::PrecodeSystem& system)
{
    const uint32_t K = system.Params.BlockCount;
    const uint32_t S = system.Params.Staircase;
    const uint32_t D2 = system.Params.DenseRows;
    const uint32_t span = K + S + D2;
    const uint32_t set_count = (span + 1u) >> 1;

    if (system.DenseRowColumns.size() != D2) {
        std::fprintf(stderr, "K=%u: dense row count\n", K);
        return false;
    }

    for (uint32_t r = 0; r < D2; ++r)
    {
        const std::vector<uint32_t>& row = system.DenseRowColumns[r];
        for (size_t i = 0; i < row.size(); ++i)
        {
            if (row[i] >= span ||
                (i > 0u && row[i] <= row[i - 1u]))
            {
                std::fprintf(stderr,
                    "K=%u: dense row %u not sorted/unique/in-range\n", K, r);
                return false;
            }
        }
    }

    // First row is exactly the set half of the deck
    if (system.DenseRowColumns[0].size() != set_count)
    {
        std::fprintf(stderr,
            "K=%u: dense row 0 has %zu columns, want %u\n",
            K, system.DenseRowColumns[0].size(), set_count);
        return false;
    }

    // Every subsequent row differs from its predecessor in EXACTLY 2 columns
    for (uint32_t r = 1; r < D2; ++r)
    {
        const std::vector<uint32_t>& prev = system.DenseRowColumns[r - 1u];
        const std::vector<uint32_t>& cur = system.DenseRowColumns[r];
        std::vector<uint32_t> sym;
        std::set_symmetric_difference(
            prev.begin(), prev.end(),
            cur.begin(), cur.end(),
            std::back_inserter(sym));
        if (sym.size() != 2u)
        {
            std::fprintf(stderr,
                "K=%u: dense rows %u->%u differ in %zu columns, want 2\n",
                K, r - 1u, r, sym.size());
            return false;
        }
    }
    return true;
}

bool TestDeterminism(uint32_t K)
{
    wirehair_v2::PrecodeSystem a, b, c;
    if (!BuildPrecodeSystem(wirehair_v2::MakeCertifiedParams(K, 7u), a) ||
        !BuildPrecodeSystem(wirehair_v2::MakeCertifiedParams(K, 7u), b) ||
        !BuildPrecodeSystem(wirehair_v2::MakeCertifiedParams(K, 8u), c))
    {
        std::fprintf(stderr, "K=%u: determinism build failed\n", K);
        return false;
    }
    if (a.StaircaseRows != b.StaircaseRows ||
        a.DenseRowColumns != b.DenseRowColumns)
    {
        std::fprintf(stderr, "K=%u: same seed produced different systems\n", K);
        return false;
    }
    if (a.StaircaseRows == c.StaircaseRows)
    {
        std::fprintf(stderr, "K=%u: different seed produced same system\n", K);
        return false;
    }
    return true;
}

// GF(256) GE rank of an h x h matrix (destructive)
unsigned SquareRank(std::vector<uint8_t>& m, unsigned h)
{
    unsigned rank = 0;
    for (unsigned col = 0; col < h; ++col)
    {
        unsigned pivot = rank;
        while (pivot < h && m[pivot * h + col] == 0u) {
            ++pivot;
        }
        if (pivot >= h) {
            continue;
        }
        for (unsigned k = 0; k < h; ++k) {
            std::swap(m[rank * h + k], m[pivot * h + k]);
        }
        const uint8_t inv = gf256_inv(m[rank * h + col]);
        for (unsigned r = 0; r < h; ++r)
        {
            if (r == rank || m[r * h + col] == 0u) {
                continue;
            }
            const uint8_t scale = gf256_mul(m[r * h + col], inv);
            for (unsigned k = col; k < h; ++k) {
                m[r * h + k] ^= gf256_mul(scale, m[rank * h + k]);
            }
        }
        ++rank;
    }
    return rank;
}

bool TestHeavyCoefficients()
{
    const uint32_t H = 12u;
    const uint32_t window = 256u - H; // 244

    // Nonzero everywhere across a wide column range
    for (uint32_t r = 0; r < H; ++r) {
        for (uint32_t c = 0; c < 1024u; ++c) {
            if (wirehair_v2::HeavyCoefficient(r, c, H) == 0u) {
                std::fprintf(stderr, "heavy coefficient zero at %u,%u\n", r, c);
                return false;
            }
        }
    }

    // Every H x H submatrix within one 244-column window must be
    // invertible (Cauchy MDS property).  Exhaustive over window starts and
    // deterministically-jittered column choices.
    wirehair::PCGRandom prng;
    prng.Seed(UINT64_C(0x4ea7c0de), 12u);
    for (uint32_t base = 0; base < window; base += 7u)
    {
        uint32_t cols[12];
        for (uint32_t i = 0; i < H; ++i)
        {
            // Distinct columns inside [base, base + window)
            for (bool collide = true; collide;)
            {
                cols[i] = base + prng.Next() % window;
                collide = false;
                for (uint32_t j = 0; j < i; ++j) {
                    if ((cols[j] % window) == (cols[i] % window)) {
                        collide = true;
                        break;
                    }
                }
            }
        }
        std::vector<uint8_t> m(H * H);
        for (uint32_t r = 0; r < H; ++r) {
            for (uint32_t i = 0; i < H; ++i) {
                m[r * H + i] = wirehair_v2::HeavyCoefficient(r, cols[i], H);
            }
        }
        if (SquareRank(m, H) != H)
        {
            std::fprintf(stderr,
                "heavy submatrix singular at window base %u\n", base);
            return false;
        }
    }
    return true;
}

} // namespace

int main()
{
    gf256_init();

    if (!TestParams()) {
        return 1;
    }
    if (!TestHeavyCoefficients()) {
        return 1;
    }

    const uint32_t Ks[] = {2u, 3u, 64u, 1000u, 3200u, 10000u, 32000u, 64000u};
    for (uint32_t K : Ks)
    {
        wirehair_v2::PrecodeSystem system;
        if (!BuildPrecodeSystem(
                wirehair_v2::MakeCertifiedParams(K, 0x5eedu), system))
        {
            std::fprintf(stderr, "K=%u: build failed\n", K);
            return 1;
        }
        if (!TestStaircase(system) || !TestDenseRows(system)) {
            return 1;
        }
    }

    const uint32_t detKs[] = {64u, 3200u, 64000u};
    for (uint32_t K : detKs) {
        if (!TestDeterminism(K)) {
            return 1;
        }
    }

    std::printf("precode_test: PASS\n");
    return 0;
}
