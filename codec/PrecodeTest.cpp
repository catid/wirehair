#include "WirehairV2Precode.h"

#include "../WirehairTools.h"
#include "../gf256.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <iterator>
#include <new>
#include <stdexcept>
#include <utility>
#include <vector>

// Structural invariant tests for the certified precode construction
// (wirehair-axd Phase 1).  These validate the documented ldpcdense_s2
// structure rules; reliability itself is certified by precode_sim and,
// end-to-end, by the later solver phases.

namespace {

bool HaveSameBinaryRows(
    const wirehair_v2::PrecodeSystem& left,
    const wirehair_v2::PrecodeSystem& right)
{
    return left.BinaryRowOffsets == right.BinaryRowOffsets &&
        left.BinaryRowColumns == right.BinaryRowColumns;
}

bool HaveSameRowRange(
    const wirehair_v2::PrecodeSystem& left,
    const wirehair_v2::PrecodeSystem& right,
    uint32_t first,
    uint32_t count)
{
    for (uint32_t row = 0u; row < count; ++row) {
        if (left.BinaryRow(first + row) != right.BinaryRow(first + row)) {
            return false;
        }
    }
    return true;
}

bool HaveSameParams(
    const wirehair_v2::PrecodeParams& left,
    const wirehair_v2::PrecodeParams& right)
{
    return left.BlockCount == right.BlockCount &&
        left.Staircase == right.Staircase &&
        left.DenseRows == right.DenseRows &&
        left.HeavyRows == right.HeavyRows &&
        left.SourceHits == right.SourceHits &&
        left.DenseIdentityCorner == right.DenseIdentityCorner &&
        left.HeavyFamily == right.HeavyFamily &&
        left.DenseAnchors == right.DenseAnchors &&
        left.Seed == right.Seed;
}

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
class PrecodeBuildAllocationFailureScope
{
public:
    PrecodeBuildAllocationFailureScope(
        wirehair_v2::test::PrecodeBuildAllocationFailureException exception)
    {
        wirehair_v2::test::
            SetPrecodeBuildAllocationFailurePointForTesting(
                wirehair_v2::test::
                    PrecodeBuildAllocationFailurePoint::FinalValidation,
                exception);
    }

    ~PrecodeBuildAllocationFailureScope()
    {
        wirehair_v2::test::
            SetPrecodeBuildAllocationFailurePointForTesting(
                wirehair_v2::test::
                    PrecodeBuildAllocationFailurePoint::None,
                wirehair_v2::test::
                    PrecodeBuildAllocationFailureException::BadAlloc);
    }

private:
    PrecodeBuildAllocationFailureScope(
        const PrecodeBuildAllocationFailureScope&);
    PrecodeBuildAllocationFailureScope& operator=(
        const PrecodeBuildAllocationFailureScope&);
};
#endif

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
            params.HeavyFamily !=
                wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy ||
            params.DenseAnchors !=
                wirehair_v2::DenseAnchorLayout::Disabled ||
            params.SourceHits != c.SourceHits)
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
    sentinel.BinaryRowOffsets = { 0u, 2u, 3u };
    sentinel.BinaryRowColumns = { 1u, 2u, 4u };

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
    invalid.DenseAnchors =
        static_cast<wirehair_v2::DenseAnchorLayout>(UINT32_MAX);
    invalid_params.push_back(invalid);
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    invalid = wirehair_v2::MakeCertifiedParams(16u, 1u);
    invalid.DenseAnchors = wirehair_v2::DenseAnchorLayout::Two07;
    invalid.DenseRows = 11u;
    invalid_params.push_back(invalid);
    invalid = wirehair_v2::MakeCertifiedParams(16u, 1u);
    invalid.DenseAnchors = wirehair_v2::DenseAnchorLayout::Four0369;
    invalid.HeavyRows = 11u;
    invalid_params.push_back(invalid);
    invalid = wirehair_v2::MakeCertifiedParams(16u, 1u);
    invalid.DenseAnchors = wirehair_v2::DenseAnchorLayout::Two07;
    invalid.DenseIdentityCorner = true;
    invalid_params.push_back(invalid);
    invalid = wirehair_v2::MakeCertifiedParams(16u, 1u);
    invalid.DenseAnchors = wirehair_v2::DenseAnchorLayout::Four0369;
    invalid.HeavyFamily =
        wirehair_v2::HeavyCoefficientFamily::HashedNonzero;
    invalid_params.push_back(invalid);
#endif

    for (size_t i = 0; i < invalid_params.size(); ++i)
    {
        wirehair_v2::PrecodeSystem out = sentinel;
        if (wirehair_v2::BuildPrecodeSystem(invalid_params[i], out) ||
            out.Params.BlockCount != sentinel.Params.BlockCount ||
            !HaveSameBinaryRows(out, sentinel))
        {
            std::fprintf(stderr,
                "invalid parameter case %zu must fail before modifying output\n",
                i);
            return false;
        }
    }
    return true;
}

bool TestMalformedFlatOffsets()
{
    wirehair_v2::PrecodeSystem valid;
    if (!wirehair_v2::BuildPrecodeSystem(
            wirehair_v2::MakeCertifiedParams(
                64u, UINT64_C(0x4353524f46465345)),
            valid) ||
        !wirehair_v2::HasValidPrecodeSystemShape(valid) ||
        !wirehair_v2::ValidatePrecodeSystem(valid))
    {
        std::fprintf(stderr, "flat-offset fixture did not validate\n");
        return false;
    }

    const auto rejected = [](const char* name,
                             const wirehair_v2::PrecodeSystem& system) {
        if (wirehair_v2::HasValidPrecodeSystemShape(system) ||
            wirehair_v2::ValidatePrecodeSystem(system))
        {
            std::fprintf(stderr,
                "flat-offset validator accepted %s\n", name);
            return false;
        }
        return true;
    };

    wirehair_v2::PrecodeSystem malformed = valid;
    malformed.BinaryRowOffsets.clear();
    if (!rejected("empty offsets", malformed)) {
        return false;
    }
    malformed = valid;
    malformed.BinaryRowOffsets.pop_back();
    if (!rejected("missing terminal offset", malformed)) {
        return false;
    }
    malformed = valid;
    malformed.BinaryRowOffsets.push_back(
        malformed.BinaryRowOffsets.back());
    if (!rejected("extra terminal offset", malformed)) {
        return false;
    }
    malformed = valid;
    malformed.BinaryRowOffsets[0] = 1u;
    if (!rejected("nonzero origin", malformed)) {
        return false;
    }
    malformed = valid;
    malformed.BinaryRowOffsets[2] =
        malformed.BinaryRowOffsets[1] - 1u;
    if (!rejected("decreasing offsets", malformed)) {
        return false;
    }
    malformed = valid;
    malformed.BinaryRowOffsets[1] =
        (uint32_t)malformed.BinaryRowColumns.size() + 1u;
    if (!rejected("out-of-range offset", malformed)) {
        return false;
    }
    malformed = valid;
    --malformed.BinaryRowOffsets.back();
    if (!rejected("short terminal offset", malformed)) {
        return false;
    }
    malformed = valid;
    ++malformed.BinaryRowOffsets.back();
    if (!rejected("long terminal offset", malformed)) {
        return false;
    }
    return true;
}

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
bool TestBuildAllocationTransaction()
{
    wirehair_v2::PrecodeSystem retained;
    if (!wirehair_v2::BuildPrecodeSystem(
            wirehair_v2::MakeCertifiedParams(
                64u, UINT64_C(0x52455441494e4544)),
            retained) ||
        !wirehair_v2::ValidatePrecodeSystem(retained))
    {
        std::fprintf(stderr, "precode OOM retained fixture failed\n");
        return false;
    }
    const wirehair_v2::PrecodeSystem retained_before = retained;
    const uint32_t* const offsets_data = retained.BinaryRowOffsets.data();
    const uint32_t* const columns_data = retained.BinaryRowColumns.data();
    const size_t offsets_capacity = retained.BinaryRowOffsets.capacity();
    const size_t columns_capacity = retained.BinaryRowColumns.capacity();
    const wirehair_v2::PrecodeParams replacement_params =
        wirehair_v2::MakeCertifiedParams(
            2048u, UINT64_C(0x4f4f4d43414e4449));

    using wirehair_v2::test::PrecodeBuildAllocationFailureException;
    const PrecodeBuildAllocationFailureException exceptions[] = {
        PrecodeBuildAllocationFailureException::BadAlloc,
        PrecodeBuildAllocationFailureException::LengthError
    };
    for (PrecodeBuildAllocationFailureException exception : exceptions)
    {
        bool caught_expected = false;
        uint32_t hits = 0u;
        {
            PrecodeBuildAllocationFailureScope failure(exception);
            try {
                (void)wirehair_v2::BuildPrecodeSystem(
                    replacement_params, retained);
            }
            catch (const std::bad_alloc&) {
                caught_expected = exception ==
                    PrecodeBuildAllocationFailureException::BadAlloc;
            }
            catch (const std::length_error&) {
                caught_expected = exception ==
                    PrecodeBuildAllocationFailureException::LengthError;
            }
            catch (...) {
                caught_expected = false;
            }
            hits = wirehair_v2::test::
                PrecodeBuildAllocationFailureHitsForTesting();
        }
        if (!caught_expected || hits != 1u ||
            !HaveSameParams(retained.Params, retained_before.Params) ||
            !HaveSameBinaryRows(retained, retained_before) ||
            retained.BinaryRowOffsets.data() != offsets_data ||
            retained.BinaryRowColumns.data() != columns_data ||
            retained.BinaryRowOffsets.capacity() != offsets_capacity ||
            retained.BinaryRowColumns.capacity() != columns_capacity ||
            !wirehair_v2::ValidatePrecodeSystem(retained))
        {
            std::fprintf(stderr,
                "precode OOM build escaped or modified retained graph\n");
            return false;
        }
    }

    wirehair_v2::PrecodeParams invalid = replacement_params;
    invalid.BlockCount = 1u;
    uint32_t invalid_hits = UINT32_MAX;
    bool invalid_built = true;
    {
        PrecodeBuildAllocationFailureScope failure(
            PrecodeBuildAllocationFailureException::BadAlloc);
        invalid_built = wirehair_v2::BuildPrecodeSystem(invalid, retained);
        invalid_hits = wirehair_v2::test::
            PrecodeBuildAllocationFailureHitsForTesting();
    }
    if (invalid_built || invalid_hits != 0u ||
        !HaveSameParams(retained.Params, retained_before.Params) ||
        !HaveSameBinaryRows(retained, retained_before) ||
        retained.BinaryRowOffsets.data() != offsets_data ||
        retained.BinaryRowColumns.data() != columns_data ||
        retained.BinaryRowOffsets.capacity() != offsets_capacity ||
        retained.BinaryRowColumns.capacity() != columns_capacity)
    {
        std::fprintf(stderr,
            "invalid precode build reached hook or modified output\n");
        return false;
    }

    if (!wirehair_v2::BuildPrecodeSystem(replacement_params, retained) ||
        retained.Params.BlockCount != replacement_params.BlockCount ||
        !wirehair_v2::ValidatePrecodeSystem(retained))
    {
        std::fprintf(stderr,
            "precode build did not recover after allocation injection\n");
        return false;
    }
    return true;
}
#endif

bool TestStaircase(const wirehair_v2::PrecodeSystem& system)
{
    const uint32_t K = system.Params.BlockCount;
    const uint32_t S = system.Params.Staircase;
    const uint32_t N1 = system.Params.SourceHits;
    const uint32_t hits = N1 < S ? N1 : S;

    if (system.BinaryRowCount() != (size_t)S + system.Params.DenseRows) {
        std::fprintf(stderr, "K=%u: staircase row count\n", K);
        return false;
    }

    std::vector<uint32_t> source_hits(K, 0);
    for (uint32_t j = 0; j < S; ++j)
    {
        const wirehair_v2::PrecodeRowView row = system.StaircaseRow(j);

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

bool IsAnchorRow(wirehair_v2::DenseAnchorLayout layout, uint32_t row)
{
    if (row == 0u) {
        return true;
    }
    if (layout == wirehair_v2::DenseAnchorLayout::Two07) {
        return row == 7u;
    }
    if (layout == wirehair_v2::DenseAnchorLayout::Four0369) {
        return row == 3u || row == 6u || row == 9u;
    }
    return false;
}

bool ReconstructDenseRows(
    const wirehair_v2::PrecodeSystem& system,
    std::vector<std::vector<uint32_t>>& rows_out)
{
    const uint32_t D2 = system.Params.DenseRows;
    if (system.BinaryRowCount() !=
        (size_t)system.Params.Staircase + D2)
    {
        return false;
    }
    std::vector<std::vector<uint32_t>> rows(D2);
    for (uint32_t row = 0; row < D2; ++row)
    {
        if (IsAnchorRow(system.Params.DenseAnchors, row)) {
            const wirehair_v2::PrecodeRowView basis =
                system.DenseBasisRow(row);
            rows[row].assign(basis.begin(), basis.end());
            continue;
        }
        const wirehair_v2::PrecodeRowView basis =
            system.DenseBasisRow(row);
        std::set_symmetric_difference(
            rows[row - 1u].begin(), rows[row - 1u].end(),
            basis.begin(), basis.end(),
            std::back_inserter(rows[row]));
    }
    rows_out.swap(rows);
    return true;
}

bool TestDenseRows(const wirehair_v2::PrecodeSystem& system)
{
    const uint32_t K = system.Params.BlockCount;
    const uint32_t S = system.Params.Staircase;
    const uint32_t D2 = system.Params.DenseRows;
    const uint32_t span = K + S + D2;
    const uint32_t set_count = (span + 1u) >> 1;

    std::vector<std::vector<uint32_t>> dense_rows;
    if (!ReconstructDenseRows(system, dense_rows) ||
        dense_rows.size() != D2)
    {
        std::fprintf(stderr, "K=%u: dense row count\n", K);
        return false;
    }

    for (uint32_t r = 0; r < D2; ++r)
    {
        const std::vector<uint32_t>& row = dense_rows[r];
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
    if (dense_rows[0].size() != set_count)
    {
        std::fprintf(stderr,
            "K=%u: dense row 0 has %zu columns, want %u\n",
            K, dense_rows[0].size(), set_count);
        return false;
    }

    // Every subsequent row differs from its predecessor in EXACTLY 2
    // columns, and within each reshuffle half the flip pairs come from
    // distinct deck positions, so they must be pairwise disjoint.  A missing
    // reshuffle (half 2 reusing half 1's flips) creates exact linear
    // dependences among the D2 rows — the failure class these rows exist to
    // prevent — and would only be caught here.
    std::vector<std::vector<uint32_t>> flips(D2);
    for (uint32_t r = 1; r < D2; ++r)
    {
        const std::vector<uint32_t>& prev = dense_rows[r - 1u];
        const std::vector<uint32_t>& cur = dense_rows[r];
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
        const uint32_t begin = half == 0u ? 1u : half1_end;
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
    if (!HaveSameBinaryRows(a, b))
    {
        std::fprintf(stderr, "K=%u: same seed produced different systems\n", K);
        return false;
    }
    if (HaveSameRowRange(a, c, 0u, a.Params.Staircase) ||
        HaveSameRowRange(
            a, c, a.Params.Staircase, a.Params.DenseRows))
    {
        std::fprintf(stderr, "K=%u: different seed produced same system\n", K);
        return false;
    }
    return true;
}

void MixGraphFingerprintU32(uint64_t& hash, uint32_t value)
{
    for (uint32_t shift = 0u; shift < 32u; shift += 8u) {
        hash ^= (uint8_t)(value >> shift);
        hash *= UINT64_C(1099511628211);
    }
}

void MixGraphFingerprintU64(uint64_t& hash, uint64_t value)
{
    MixGraphFingerprintU32(hash, (uint32_t)value);
    MixGraphFingerprintU32(hash, (uint32_t)(value >> 32));
}

void MixGraphFingerprint(
    uint64_t& hash,
    const wirehair_v2::PrecodeSystem& system)
{
    const wirehair_v2::PrecodeParams& params = system.Params;
    MixGraphFingerprintU32(hash, params.BlockCount);
    MixGraphFingerprintU32(hash, params.Staircase);
    MixGraphFingerprintU32(hash, params.DenseRows);
    MixGraphFingerprintU32(hash, params.HeavyRows);
    MixGraphFingerprintU32(hash, params.SourceHits);
    MixGraphFingerprintU32(hash, params.DenseIdentityCorner ? 1u : 0u);
    MixGraphFingerprintU32(hash, (uint32_t)params.DenseAnchors);
    MixGraphFingerprintU32(hash, (uint32_t)params.HeavyFamily);
    MixGraphFingerprintU64(hash, params.Seed);
    MixGraphFingerprintU32(hash, params.Staircase);
    for (uint32_t row_index = 0;
         row_index < params.Staircase;
         ++row_index)
    {
        const wirehair_v2::PrecodeRowView row =
            system.StaircaseRow(row_index);
        MixGraphFingerprintU32(hash, (uint32_t)row.size());
        for (uint32_t column : row) {
            MixGraphFingerprintU32(hash, column);
        }
    }
    MixGraphFingerprintU32(hash, params.DenseRows);
    for (uint32_t row_index = 0;
         row_index < params.DenseRows;
         ++row_index)
    {
        const wirehair_v2::PrecodeRowView row =
            system.DenseBasisRow(row_index);
        MixGraphFingerprintU32(hash, (uint32_t)row.size());
        for (uint32_t column : row) {
            MixGraphFingerprintU32(hash, column);
        }
    }
}

bool TestExactGraphGolden()
{
    // This digest was frozen from 95e6fd0 before the build-only allocation
    // cleanup.  It covers both sides of the inline scratch boundaries, the
    // N1 transition, the complete production range endpoints, three seeds,
    // every feasible identity-corner variant, and a valid maximum-span
    // zero-mean staircase geometry.  Allocation-only changes must reproduce
    // the exact staircase and dense-basis column streams.
    static const uint64_t kGolden = UINT64_C(0xc7924ce359571489);
    static const uint32_t kBlockCounts[] = {
        2u, 3u, 4u, 8u, 64u, 100u, 450u, 500u, 1000u, 1024u,
        1025u, 2048u, 10000u, 32000u, 64000u
    };
    uint64_t hash = UINT64_C(14695981039346656037);
    for (uint32_t K : kBlockCounts)
    {
        for (uint32_t seed = 0u; seed < 3u; ++seed)
        {
            wirehair_v2::PrecodeParams params =
                wirehair_v2::MakeCertifiedParams(
                    K, UINT64_C(0x50524f4245530000) ^ seed ^ K);
            wirehair_v2::PrecodeSystem system;
            if (!BuildPrecodeSystem(params, system)) {
                return false;
            }
            MixGraphFingerprint(hash, system);

            const bool identity_feasible =
                K + params.Staircase >= 2u * (params.DenseRows >> 1);
            if (!identity_feasible) {
                continue;
            }
            params.Seed = UINT64_C(0x4944454e54495459) ^ seed ^ K;
            params.DenseIdentityCorner = true;
            if (!BuildPrecodeSystem(params, system)) {
                return false;
            }
            MixGraphFingerprint(hash, system);
        }
    }

    wirehair_v2::PrecodeParams exotic;
    exotic.BlockCount = 2u;
    exotic.Staircase = 65533u;
    exotic.DenseRows = 0u;
    exotic.HeavyRows = 0u;
    exotic.SourceHits = 1u;
    exotic.Seed = UINT64_C(0x5a45524f4d45414e);
    wirehair_v2::PrecodeSystem exotic_system;
    if (!BuildPrecodeSystem(exotic, exotic_system)) {
        return false;
    }
    const size_t expected_exotic_references =
        (size_t)exotic.Staircase * 2u - 1u + exotic.BlockCount;
    if (exotic_system.BinaryRowOffsets.size() !=
            (size_t)exotic.Staircase + 1u ||
        exotic_system.BinaryRowColumns.size() !=
            expected_exotic_references)
    {
        std::fprintf(stderr,
            "zero-mean staircase flat shape offsets=%zu refs=%zu, "
            "want %u/%zu\n",
            exotic_system.BinaryRowOffsets.size(),
            exotic_system.BinaryRowColumns.size(),
            exotic.Staircase + 1u,
            expected_exotic_references);
        return false;
    }
    MixGraphFingerprint(hash, exotic_system);
    if (hash != kGolden)
    {
        std::fprintf(stderr,
            "precode graph fingerprint=%016llx want=%016llx\n",
            (unsigned long long)hash,
            (unsigned long long)kGolden);
        return false;
    }
    return true;
}

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
uint64_t DenseRowsFingerprint(const wirehair_v2::PrecodeSystem& system)
{
    std::vector<std::vector<uint32_t>> dense_rows;
    if (!ReconstructDenseRows(system, dense_rows)) {
        return 0u;
    }
    uint64_t hash = UINT64_C(14695981039346656037);
    const auto mix_u32 = [&](uint32_t value) {
        for (uint32_t shift = 0u; shift < 32u; shift += 8u) {
            hash ^= (uint8_t)(value >> shift);
            hash *= UINT64_C(1099511628211);
        }
    };
    mix_u32(system.Params.BlockCount);
    mix_u32((uint32_t)dense_rows.size());
    for (size_t row = 0u; row < dense_rows.size(); ++row)
    {
        mix_u32((uint32_t)row);
        mix_u32((uint32_t)dense_rows[row].size());
        for (uint32_t column : dense_rows[row]) {
            mix_u32(column);
        }
    }
    return hash;
}

bool TestDenseAnchorLayout(
    uint32_t K,
    uint64_t seed,
    wirehair_v2::DenseAnchorLayout layout)
{
    wirehair_v2::PrecodeParams params =
        wirehair_v2::MakeCertifiedParams(K, seed);
    params.DenseAnchors = layout;
    wirehair_v2::PrecodeSystem system, repeat;
    if (!BuildPrecodeSystem(params, system) ||
        !BuildPrecodeSystem(params, repeat) ||
        !wirehair_v2::ValidatePrecodeSystem(system) ||
        !HaveSameBinaryRows(system, repeat))
    {
        std::fprintf(stderr,
            "K=%u layout=%u: anchor build/validation/determinism failed\n",
            K, (unsigned)layout);
        return false;
    }

    const uint32_t span = K + params.Staircase + params.DenseRows;
    const size_t anchor_size = (span + 1u) / 2u;
    std::vector<std::vector<uint32_t>> dense_rows;
    if (system.BinaryRowCount() !=
            (size_t)params.Staircase + params.DenseRows ||
        !ReconstructDenseRows(system, dense_rows))
    {
        std::fprintf(stderr,
            "K=%u layout=%u: dense basis build failed\n",
            K, (unsigned)layout);
        return false;
    }

    size_t expected_references = 0u;
    size_t actual_references = 0u;
    std::vector<uint32_t> reconstructed;
    std::vector<uint32_t> segment_columns;
    for (uint32_t row = 0; row < params.DenseRows; ++row)
    {
        const bool anchor = IsAnchorRow(layout, row);
        const size_t expected_size = anchor ? anchor_size : 2u;
        const wirehair_v2::PrecodeRowView basis =
            system.DenseBasisRow(row);
        if (basis.size() != expected_size)
        {
            std::fprintf(stderr,
                "K=%u layout=%u row=%u: basis refs=%zu want=%zu\n",
                K, (unsigned)layout, row,
                basis.size(), expected_size);
            return false;
        }
        expected_references += expected_size;
        actual_references += basis.size();

        if (anchor) {
            reconstructed.assign(basis.begin(), basis.end());
            segment_columns.clear();
        }
        else
        {
            if (layout == wirehair_v2::DenseAnchorLayout::Disabled &&
                row == 1u + (params.DenseRows >> 1))
            {
                segment_columns.clear();
            }
            segment_columns.insert(
                segment_columns.end(), basis.begin(), basis.end());
            std::sort(segment_columns.begin(), segment_columns.end());
            if (std::adjacent_find(
                    segment_columns.begin(), segment_columns.end()) !=
                segment_columns.end())
            {
                std::fprintf(stderr,
                    "K=%u layout=%u row=%u: repeated segment delta column\n",
                    K, (unsigned)layout, row);
                return false;
            }
            std::vector<uint32_t> next;
            std::set_symmetric_difference(
                reconstructed.begin(), reconstructed.end(),
                basis.begin(), basis.end(),
                std::back_inserter(next));
            reconstructed.swap(next);
        }
        if (reconstructed != dense_rows[row])
        {
            std::fprintf(stderr,
                "K=%u layout=%u row=%u: basis reconstruction mismatch\n",
                K, (unsigned)layout, row);
            return false;
        }
    }
    if (actual_references != expected_references) {
        return false;
    }
    if (K == 64u)
    {
        wirehair_v2::PrecodeSystem bad = system;
        wirehair_v2::test::ReplacePrecodeRowForTesting(
            bad, (size_t)params.Staircase + 1u, {});
        if (wirehair_v2::ValidatePrecodeSystem(bad)) {
            std::fprintf(stderr,
                "layout=%u: validator accepted missing direct delta\n",
                (unsigned)layout);
            return false;
        }
        bad = system;
        std::vector<uint32_t> damaged_anchor =
            wirehair_v2::test::CopyPrecodeRowForTesting(
                bad.DenseBasisRow(0u));
        damaged_anchor.pop_back();
        wirehair_v2::test::ReplacePrecodeRowForTesting(
            bad, params.Staircase, damaged_anchor);
        if (wirehair_v2::ValidatePrecodeSystem(bad)) {
            std::fprintf(stderr,
                "layout=%u: validator accepted damaged direct anchor\n",
                (unsigned)layout);
            return false;
        }
    }
    return true;
}

bool TestDenseAnchorGoldens()
{
    struct GoldenCase
    {
        wirehair_v2::DenseAnchorLayout Layout;
        uint64_t Fingerprint;
    };
    // These are the pure-binary dense-equation fingerprints first measured
    // by historical commit 57d7990.  Pinning them proves this GF(256)-only
    // port did not inherit or alter any retired completion-field behavior.
    static const GoldenCase kCases[] = {
        { wirehair_v2::DenseAnchorLayout::Disabled,
          UINT64_C(0xcde1e21be25b9081) },
        { wirehair_v2::DenseAnchorLayout::Two07,
          UINT64_C(0x7da0674ba8931e64) },
        { wirehair_v2::DenseAnchorLayout::Four0369,
          UINT64_C(0x89f570b207ae565d) }
    };
    for (const GoldenCase& golden : kCases)
    {
        wirehair_v2::PrecodeParams params =
            wirehair_v2::MakeCertifiedParams(
                945u, UINT64_C(0x5eed7a12));
        params.DenseAnchors = golden.Layout;
        wirehair_v2::PrecodeSystem system;
        if (!BuildPrecodeSystem(params, system) ||
            DenseRowsFingerprint(system) != golden.Fingerprint)
        {
            std::fprintf(stderr,
                "K=945 layout=%u: dense equation golden mismatch\n",
                (unsigned)golden.Layout);
            return false;
        }
    }
    return true;
}
#endif

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
    static const uint32_t kHeavyRows[] = {11u, 12u, 13u};
    for (uint32_t H : kHeavyRows)
    {
        const uint32_t window = 256u - H;

        // Nonzero everywhere across a wide column range
        for (uint32_t r = 0; r < H; ++r) {
            for (uint32_t c = 0; c < 1024u; ++c) {
                if (wirehair_v2::HeavyCoefficient(r, c, H) == 0u) {
                    std::fprintf(stderr,
                        "H=%u: heavy coefficient zero at %u,%u\n",
                        H, r, c);
                    return false;
                }
            }
        }

        // H x H submatrices within one coefficient period must be invertible
        // (Cauchy MDS property).  The coefficient depends only on c modulo
        // 256-H, so sampled distinct-column subsets of one period cover the
        // structure; the MDS guarantee itself is analytic (Cauchy
        // determinant).
        wirehair::PCGRandom prng;
        prng.Seed(UINT64_C(0x4ea7c0de), H);
        for (uint32_t base = 0; base < window; base += 7u)
        {
            std::vector<uint32_t> cols(H);
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
                    m[r * H + i] =
                        wirehair_v2::HeavyCoefficient(r, cols[i], H);
                }
            }
            if (SquareRank(m, H) != H)
            {
                std::fprintf(stderr,
                    "H=%u: heavy submatrix singular at window base %u\n",
                    H, base);
                return false;
            }
        }
    }
    return true;
}

bool TestCandidateGeometries()
{
    struct Geometry
    {
        uint32_t DenseRows;
        uint32_t HeavyRows;
    };
    static const uint32_t kBlockCounts[] = {3u, 4u, 1001u};
    static const Geometry kGeometries[] = {
        {12u, 11u},
        {12u, 13u},
        {13u, 12u}
    };
    static const wirehair_v2::HeavyCoefficientFamily kFamilies[] = {
        wirehair_v2::HeavyCoefficientFamily::PeriodicCauchy,
        wirehair_v2::HeavyCoefficientFamily::HashedNonzero
    };

    for (uint32_t K : kBlockCounts)
    {
        for (const Geometry& geometry : kGeometries)
        {
            for (wirehair_v2::HeavyCoefficientFamily family : kFamilies)
            {
                wirehair_v2::PrecodeParams params =
                    wirehair_v2::MakeCertifiedParams(
                        K,
                        UINT64_C(0x4448434f56455241) ^
                            ((uint64_t)K << 24) ^
                            ((uint64_t)geometry.DenseRows << 16) ^
                            ((uint64_t)geometry.HeavyRows << 8) ^
                            (uint32_t)family);
                params.DenseRows = geometry.DenseRows;
                params.HeavyRows = geometry.HeavyRows;
                params.HeavyFamily = family;

                wirehair_v2::PrecodeSystem system;
                if (!BuildPrecodeSystem(params, system) ||
                    !wirehair_v2::ValidatePrecodeSystem(system) ||
                    system.Params.DenseRows != geometry.DenseRows ||
                    system.Params.HeavyRows != geometry.HeavyRows ||
                    system.Params.HeavyFamily != family ||
                    !TestStaircase(system) || !TestDenseRows(system))
                {
                    std::fprintf(stderr,
                        "candidate geometry failed K=%u D=%u H=%u family=%u\n",
                        K, geometry.DenseRows, geometry.HeavyRows,
                        (unsigned)family);
                    return false;
                }
            }
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
    if (!TestMalformedFlatOffsets()) {
        return 1;
    }
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    if (!TestBuildAllocationTransaction()) {
        return 1;
    }
#endif
    if (!TestHeavyCoefficients()) {
        return 1;
    }
    if (!TestCandidateGeometries()) {
        return 1;
    }
    if (!TestExactGraphGolden()) {
        return 1;
    }

    const uint32_t Ks[] = {
        2u, 3u, 4u, 64u, 1000u, 1001u, 3200u, 10000u, 32000u, 64000u
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

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    static const wirehair_v2::DenseAnchorLayout kAnchorLayouts[] = {
        wirehair_v2::DenseAnchorLayout::Disabled,
        wirehair_v2::DenseAnchorLayout::Two07,
        wirehair_v2::DenseAnchorLayout::Four0369
    };
    static const uint32_t kAnchorKs[] = {2u, 64u, 1001u, 64000u};
    static const uint64_t kAnchorSeeds[] = {
        UINT64_C(0),
        UINT64_C(1),
        UINT64_C(0x616e63686f727331)
    };
    for (uint32_t K : kAnchorKs) {
        for (uint64_t seed : kAnchorSeeds) {
            for (wirehair_v2::DenseAnchorLayout layout : kAnchorLayouts) {
                if (!TestDenseAnchorLayout(K, seed ^ K, layout))
                {
                    return 1;
                }
            }
        }
    }
    if (!TestDenseAnchorGoldens()) {
        return 1;
    }
#endif

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
        std::vector<std::vector<uint32_t>> dense_rows;
        if (!ReconstructDenseRows(system, dense_rows)) {
            std::fprintf(stderr, "K=%u: ic reconstruction failed\n", K);
            return 1;
        }
        for (uint32_t r = 0; r < D2; ++r)
        {
            const std::vector<uint32_t>& row = dense_rows[r];
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
        if (dense_rows[0].size() != ((deck_span + 1u) >> 1) + 1u)
        {
            std::fprintf(stderr, "K=%u: ic row 0 size wrong\n", K);
            return 1;
        }
        // Consecutive rows: 2 deck flips + the two distinct own columns
        for (uint32_t r = 1; r < D2; ++r)
        {
            std::vector<uint32_t> sym;
            std::set_symmetric_difference(
                dense_rows[r - 1u].begin(),
                dense_rows[r - 1u].end(),
                dense_rows[r].begin(),
                dense_rows[r].end(),
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
    }

    std::printf("precode_test: PASS\n");
    return 0;
}
