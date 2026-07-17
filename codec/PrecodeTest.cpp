#include "WirehairV2Precode.h"

#include "../WirehairTools.h"
#include "../gf256.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <iterator>
#include <utility>
#include <vector>

// Structural invariant tests for the certified precode construction
// (wirehair-axd Phase 1).  These validate the documented ldpcdense_s2
// structure rules; reliability itself is certified by precode_sim and,
// end-to-end, by the later solver phases.

namespace {

bool TestParams()
{
    struct ParamCase
    {
        uint32_t K;
        uint32_t SourceHits;
    };
    const ParamCase cases[] = {
        { 3200u, 2u },
        { 9999u, 2u },
        { 10000u, 3u },
        { 32000u, 3u },
        { 64000u, 3u },
    };
    for (const ParamCase& c : cases)
    {
        const wirehair_v2::PrecodeParams params =
            wirehair_v2::MakeCertifiedParams(c.K, 1u);
        if (params.Staircase != wirehair::GetDenseCount(c.K) ||
            params.DenseRows != 12u ||
            params.HeavyRows != 12u ||
            params.Field != wirehair_v2::CompletionField::GF256 ||
            params.HeavyFamily !=
                wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy ||
            params.SourceHits != c.SourceHits ||
            params.DenseTwoAnchor || params.DenseTwoAnchorPhase != 0u)
        {
            std::fprintf(stderr,
                "certified params wrong for K=%u (N1=%u, want %u)\n",
                c.K, params.SourceHits, c.SourceHits);
            return false;
        }
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

    wirehair_v2::PrecodeSystem sentinel;
    sentinel.Params.BlockCount = 7u;
    sentinel.StaircaseRows.push_back(std::vector<uint32_t>{1u, 2u});
    sentinel.DenseRowColumns.push_back(std::vector<uint32_t>{3u});

    std::vector<wirehair_v2::PrecodeParams> invalid_params;
    wirehair_v2::PrecodeParams invalid =
        wirehair_v2::MakeCertifiedParams(16u, 1u);
    invalid.DenseRows = 65u;
    invalid_params.push_back(invalid);
    invalid = wirehair_v2::MakeCertifiedParams(16u, 1u);
    invalid.HeavyRows = 129u;
    invalid_params.push_back(invalid);
    invalid = wirehair_v2::MakeCertifiedParams(16u, 1u);
    invalid.HeavyFamily =
        static_cast<wirehair_v2::HeavyCoefficientFamily>(UINT32_MAX);
    invalid_params.push_back(invalid);
    invalid = wirehair_v2::MakeCertifiedParams(16u, 1u);
    invalid.Field = static_cast<wirehair_v2::CompletionField>(UINT32_MAX);
    invalid_params.push_back(invalid);
    invalid = wirehair_v2::MakeMixedParams(16u, 1u);
    invalid.HeavyRows = 11u;
    invalid_params.push_back(invalid);
    invalid = wirehair_v2::MakeCertifiedParams(64000u, 1u);
    invalid.Staircase = 1500u;
    invalid.DenseRows = 36u;
    invalid.HeavyRows = 0u;
    invalid_params.push_back(invalid);
    invalid = wirehair_v2::MakeCertifiedParams(64000u, 1u);
    invalid.Staircase = 1500u;
    invalid.DenseRows = 0u;
    invalid.HeavyRows = 36u;
    invalid_params.push_back(invalid);
    invalid = wirehair_v2::MakeCertifiedParams(2u, 1u);
    invalid.Staircase = 1u;
    invalid.DenseRows = 4u;
    invalid.HeavyRows = 0u;
    invalid.DenseIdentityCorner = true;
    invalid_params.push_back(invalid);
    invalid = wirehair_v2::MakeCertifiedParams(16u, 1u);
    invalid.DenseRows = 13u;
    invalid.DenseTwoAnchor = true;
    invalid_params.push_back(invalid);
    invalid = wirehair_v2::MakeCertifiedParams(16u, 1u);
    invalid.DenseIdentityCorner = true;
    invalid.DenseTwoAnchor = true;
    invalid_params.push_back(invalid);
    invalid = wirehair_v2::MakeCertifiedParams(16u, 1u);
    invalid.DenseTwoAnchorPhase = 1u;
    invalid_params.push_back(invalid);
    invalid = wirehair_v2::MakeCertifiedParams(16u, 1u);
    invalid.DenseTwoAnchor = true;
    invalid.DenseTwoAnchorPhase = 3u;
    invalid_params.push_back(invalid);

    for (size_t i = 0; i < invalid_params.size(); ++i)
    {
        wirehair_v2::PrecodeSystem out = sentinel;
        if (wirehair_v2::BuildPrecodeSystem(invalid_params[i], out) ||
            out.Params.BlockCount != sentinel.Params.BlockCount ||
            out.StaircaseRows != sentinel.StaircaseRows ||
            out.DenseRowColumns != sentinel.DenseRowColumns)
        {
            std::fprintf(stderr,
                "invalid parameter case %zu must fail before modifying output\n",
                i);
            return false;
        }
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

        for (uint32_t col : row)
        {
            if (col < K) {
                ++source_hits[col];
            }
            // Parity-range membership must be EXACTLY own + link: any other
            // parity column means a wrong staircase direction or link bug
            else if (col != K + j &&
                     !(j > 0u && col == K + j - 1u))
            {
                std::fprintf(stderr,
                    "K=%u: staircase row %u has stray parity column %u\n",
                    K, j, col);
                return false;
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

    // Each anchor is exactly the set half of its deck.
    const uint32_t second_anchor = 1u + (D2 >> 1);
    for (uint32_t r = 0; r < D2; ++r)
    {
        if ((r == 0u ||
             (system.Params.DenseTwoAnchor && r == second_anchor)) &&
            system.DenseRowColumns[r].size() != set_count)
        {
            std::fprintf(stderr,
                "K=%u: dense anchor row %u has %zu columns, want %u\n",
                K, r, system.DenseRowColumns[r].size(), set_count);
            return false;
        }
    }

    // Every non-anchor row differs from its predecessor in EXACTLY 2 columns,
    // and within each reshuffle half the flip pairs come from distinct deck
    // positions, so they must be pairwise disjoint.  A missing reshuffle
    // (half 2 reusing half 1's flips) creates exact linear dependences among
    // the D2 rows — the failure class these rows exist to prevent — and would
    // only be caught here.
    std::vector<std::vector<uint32_t>> flips(D2);
    for (uint32_t r = 1; r < D2; ++r)
    {
        if (system.Params.DenseTwoAnchor && r == second_anchor) {
            continue;
        }
        const std::vector<uint32_t>& prev = system.DenseRowColumns[r - 1u];
        const std::vector<uint32_t>& cur = system.DenseRowColumns[r];
        std::vector<uint32_t>& sym = flips[r];
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
    const uint32_t half1_end = 1u + (D2 >> 1); // flips[1 .. half1_end)
    for (uint32_t half = 0; half < 2u; ++half)
    {
        const uint32_t begin = half == 0u ? 1u :
            half1_end + (system.Params.DenseTwoAnchor ? 1u : 0u);
        const uint32_t end = half == 0u ? half1_end : D2;
        std::vector<uint32_t> seen;
        for (uint32_t r = begin; r < end; ++r) {
            seen.insert(seen.end(), flips[r].begin(), flips[r].end());
        }
        std::sort(seen.begin(), seen.end());
        if (std::adjacent_find(seen.begin(), seen.end()) != seen.end())
        {
            std::fprintf(stderr,
                "K=%u: dense flip columns repeat within half %u "
                "(reshuffle cadence broken)\n", K, half);
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
    if (a.StaircaseRows == c.StaircaseRows ||
        a.DenseRowColumns == c.DenseRowColumns)
    {
        std::fprintf(stderr, "K=%u: different seed produced same system\n", K);
        return false;
    }
    return true;
}

uint64_t DenseRowsFingerprint(const wirehair_v2::PrecodeSystem& system)
{
    uint64_t hash = UINT64_C(14695981039346656037);
    const auto mix_u32 = [&](uint32_t value) {
        for (uint32_t shift = 0u; shift < 32u; shift += 8u) {
            hash ^= (uint8_t)(value >> shift);
            hash *= UINT64_C(1099511628211);
        }
    };
    mix_u32(system.Params.BlockCount);
    mix_u32((uint32_t)system.DenseRowColumns.size());
    for (size_t r = 0u; r < system.DenseRowColumns.size(); ++r)
    {
        const std::vector<uint32_t>& row = system.DenseRowColumns[r];
        mix_u32((uint32_t)r);
        mix_u32((uint32_t)row.size());
        for (uint32_t col : row) {
            mix_u32(col);
        }
    }
    return hash;
}

bool TestTwoAnchorDeterminism(uint32_t K)
{
    wirehair_v2::PrecodeParams params =
        wirehair_v2::MakeCertifiedParams(K, UINT64_C(0x5eed7a12));
    wirehair_v2::PrecodeSystem baseline;
    if (!BuildPrecodeSystem(params, baseline)) {
        std::fprintf(stderr, "K=%u: two-anchor baseline build failed\n", K);
        return false;
    }

    params.DenseTwoAnchor = true;
    wirehair_v2::PrecodeSystem phased[3];
    for (uint32_t phase = 0u; phase <= 2u; ++phase)
    {
        params.DenseTwoAnchorPhase = phase;
        wirehair_v2::PrecodeSystem repeat;
        if (!BuildPrecodeSystem(params, phased[phase]) ||
            !BuildPrecodeSystem(params, repeat) ||
            !TestStaircase(phased[phase]) ||
            !TestDenseRows(phased[phase]))
        {
            std::fprintf(stderr,
                "K=%u: two-anchor phase %u build failed\n", K, phase);
            return false;
        }
        if (phased[phase].StaircaseRows != repeat.StaircaseRows ||
            phased[phase].DenseRowColumns != repeat.DenseRowColumns)
        {
            std::fprintf(stderr,
                "K=%u: same seed produced different two-anchor phase %u "
                "systems\n", K, phase);
            return false;
        }
        if (phased[phase].StaircaseRows != baseline.StaircaseRows)
        {
            std::fprintf(stderr,
                "K=%u: two-anchor phase %u changed staircase "
                "construction\n", K, phase);
            return false;
        }
        for (uint32_t r = 0u; r <= 6u; ++r)
        {
            if (phased[phase].DenseRowColumns[r] !=
                baseline.DenseRowColumns[r])
            {
                std::fprintf(stderr,
                    "K=%u: two-anchor phase %u changed preserved dense "
                    "row %u\n", K, phase, r);
                return false;
            }
        }
    }
    const wirehair_v2::PrecodeSystem& a = phased[0];
    // q1/q2 are exact forward windows through the same independently
    // shuffled second-deck cadence.  Pinning these overlaps proves q0 remains
    // byte-identical and the experiment changes phase rather than RNG use.
    for (uint32_t phase = 1u; phase <= 2u; ++phase)
    {
        for (uint32_t row = 7u; row + phase < 12u; ++row)
        {
            if (phased[phase].DenseRowColumns[row] !=
                phased[0].DenseRowColumns[row + phase])
            {
                std::fprintf(stderr,
                    "K=%u: two-anchor phase %u row %u is not q0 row %u\n",
                    K, phase, row, row + phase);
                return false;
            }
        }
        if (phased[phase].DenseRowColumns[11] ==
            phased[0].DenseRowColumns[11])
        {
            std::fprintf(stderr,
                "K=%u: two-anchor phase %u did not advance final row\n",
                K, phase);
            return false;
        }
    }
    if (a.DenseRowColumns[7] == baseline.DenseRowColumns[7])
    {
        std::fprintf(stderr,
            "K=%u: two-anchor failed to replace row 7\n", K);
        return false;
    }
    std::vector<uint32_t> anchor_transition;
    std::set_symmetric_difference(
        a.DenseRowColumns[6].begin(), a.DenseRowColumns[6].end(),
        a.DenseRowColumns[7].begin(), a.DenseRowColumns[7].end(),
        std::back_inserter(anchor_transition));
    if (anchor_transition.size() <= 2u)
    {
        std::fprintf(stderr,
            "K=%u: two-anchor row 6->7 transition is not dense\n", K);
        return false;
    }

    for (uint32_t phase = 0u; phase <= 2u; ++phase)
    {
        wirehair_v2::PrecodeSystem bad = phased[phase];
        bad.DenseRowColumns[7].pop_back();
        if (wirehair_v2::ValidatePrecodeSystem(bad))
        {
            std::fprintf(stderr,
                "K=%u: validator accepted an unbalanced phase %u second "
                "anchor\n", K, phase);
            return false;
        }
    }

    if (K == 945u)
    {
        const uint64_t baseline_hash = DenseRowsFingerprint(baseline);
        const uint64_t candidate_hash = DenseRowsFingerprint(a);
        const uint64_t expected_baseline =
            UINT64_C(0xcde1e21be25b9081);
        const uint64_t expected_candidate =
            UINT64_C(0x7da0674ba8931e64);
        if (baseline_hash != expected_baseline ||
            candidate_hash != expected_candidate)
        {
            std::fprintf(stderr,
                "K=945: dense construction golden mismatch "
                "baseline=%016llx candidate=%016llx\n",
                (unsigned long long)baseline_hash,
                (unsigned long long)candidate_hash);
            return false;
        }
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

    // H x H submatrices within a 244-column window must be invertible
    // (Cauchy MDS property).  The coefficient depends only on c mod 244, so
    // sampled distinct-column subsets of one window cover the structure;
    // the MDS guarantee itself is analytic (Cauchy determinant).
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

    const uint32_t Ks[] = {
        2u, 3u, 64u, 945u, 1000u, 3200u, 10000u, 32000u, 64000u
    };
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

    for (uint32_t K : Ks) {
        if (!TestTwoAnchorDeterminism(K)) {
            return 1;
        }
    }

    // Identity-corner dense variant: deck spans only K + S, and dense row
    // r references exactly its own dense column K + S + r
    for (uint32_t K : Ks)
    {
        wirehair_v2::PrecodeParams params =
            wirehair_v2::MakeCertifiedParams(K, 0x5eedu);
        params.DenseIdentityCorner = true;
        wirehair_v2::PrecodeSystem system;
        const bool feasible =
            params.BlockCount + params.Staircase >=
            2u * (params.DenseRows >> 1);
        if (!BuildPrecodeSystem(params, system))
        {
            if (feasible) {
                std::fprintf(stderr, "K=%u: ic build failed\n", K);
                return 1;
            }
            continue; // tiny K: rejection is the contract
        }
        if (!feasible) {
            std::fprintf(stderr, "K=%u: infeasible ic build accepted\n", K);
            return 1;
        }
        if (!TestStaircase(system)) {
            return 1;
        }
        const uint32_t S = params.Staircase;
        const uint32_t D2 = params.DenseRows;
        const uint32_t deck_span = K + S;
        for (uint32_t r = 0; r < D2; ++r)
        {
            const std::vector<uint32_t>& row = system.DenseRowColumns[r];
            uint32_t own = 0;
            for (uint32_t col : row)
            {
                if (col >= deck_span)
                {
                    if (col != deck_span + r) {
                        std::fprintf(stderr,
                            "K=%u: ic row %u has foreign dense column %u\n",
                            K, r, col);
                        return 1;
                    }
                    ++own;
                }
            }
            if (own != 1u) {
                std::fprintf(stderr,
                    "K=%u: ic row %u own-column count %u\n", K, r, own);
                return 1;
            }
        }
        if (system.DenseRowColumns[0].size() != ((deck_span + 1u) >> 1) + 1u)
        {
            std::fprintf(stderr, "K=%u: ic row 0 size wrong\n", K);
            return 1;
        }
        // Consecutive rows: 2 deck flips + the two distinct own columns
        for (uint32_t r = 1; r < D2; ++r)
        {
            std::vector<uint32_t> sym;
            std::set_symmetric_difference(
                system.DenseRowColumns[r - 1u].begin(),
                system.DenseRowColumns[r - 1u].end(),
                system.DenseRowColumns[r].begin(),
                system.DenseRowColumns[r].end(),
                std::back_inserter(sym));
            if (sym.size() != 4u) {
                std::fprintf(stderr,
                    "K=%u: ic rows %u->%u differ in %zu columns, want 4\n",
                    K, r - 1u, r, sym.size());
                return 1;
            }
        }
    }

    wirehair::PCGRandom random;
    random.Seed(UINT64_C(0x51a1d5eed), UINT64_C(0xb00d));
    for (uint32_t trial = 0; trial < 96u; ++trial)
    {
        const uint32_t K = 2u + random.Next() % 63999u;
        wirehair_v2::PrecodeParams params =
            wirehair_v2::MakeCertifiedParams(K, random.Next());
        wirehair_v2::PrecodeSystem system;
        if (!BuildPrecodeSystem(params, system) ||
            !wirehair_v2::ValidatePrecodeSystem(system))
        {
            std::fprintf(stderr,
                "random builder validation failed at K=%u trial=%u\n",
                K, trial);
            return 1;
        }
        if (K + params.Staircase >= 2u * (params.DenseRows >> 1))
        {
            params.DenseIdentityCorner = true;
            if (!BuildPrecodeSystem(params, system) ||
                !wirehair_v2::ValidatePrecodeSystem(system))
            {
                std::fprintf(stderr,
                    "random identity builder validation failed at K=%u\n", K);
                return 1;
            }
        }

        params.DenseIdentityCorner = false;
        params.DenseRows = 12u;
        params.DenseTwoAnchor = true;
        params.DenseTwoAnchorPhase = trial % 3u;
        if (!BuildPrecodeSystem(params, system) ||
            !wirehair_v2::ValidatePrecodeSystem(system) ||
            !TestDenseRows(system))
        {
            std::fprintf(stderr,
                "random two-anchor builder validation failed at K=%u\n", K);
            return 1;
        }
    }

    std::printf("precode_test: PASS\n");
    return 0;
}
