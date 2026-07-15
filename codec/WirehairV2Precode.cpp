#include "WirehairV2Precode.h"

#include "../WirehairTools.h"
#include "../gf256.h"

#include <algorithm>
#include <limits>

namespace wirehair_v2 {
namespace {

uint32_t CertifiedSourceHits(uint32_t block_count)
{
    // K=32000 codec-port certification rejects N1=2 at large K, while
    // K=10000 is the transition point tracked in the Phase 4 notes.
    return block_count >= 10000u ? 3u : 2u;
}

// Unbiased uniform in [0, bound) via Lemire multiply-shift rejection,
// matching the simulator's corrected Rng::Below used for the codecport
// certification (the earlier Phase B runs predate the Below threshold fix;
// its bias was <= bound/2^32, statistically immaterial at those rates)
uint32_t UniformBelow(wirehair::PCGRandom& prng, uint32_t bound)
{
    if (bound <= 1u) {
        return 0;
    }
    // 2^32 mod bound, computed in 32-bit arithmetic
    const uint32_t threshold = (0u - bound) % bound;
    for (;;)
    {
        const uint32_t x = prng.Next();
        const uint64_t m = (uint64_t)x * bound;
        if ((uint32_t)m >= threshold) {
            return (uint32_t)(m >> 32);
        }
    }
}

// Sattolo-style inside-out deck shuffle with UNBIASED position draws.
// Production ShuffleDeck16 draws 16-bit chunks modulo ii, which is fine at
// its <= ~500-entry dense-count sizes but measurably biased once the deck
// spans the whole precode (set-half membership 0.455-0.577 instead of 0.500
// at span ~64000).  The Phase B certification assumed unbiased decks, so
// this keeps the same deck structure with a rejection-sampled draw.
void UnbiasedShuffleDeck(
    wirehair::PCGRandom& prng,
    uint16_t* deck,
    uint32_t count)
{
    deck[0] = 0;
    for (uint32_t ii = 1; ii < count; ++ii)
    {
        const uint32_t jj = UniformBelow(prng, ii);
        deck[ii] = deck[jj];
        deck[jj] = (uint16_t)ii;
    }
}

bool IsStrictlyIncreasingBelow(
    const std::vector<uint32_t>& row,
    uint64_t limit)
{
    uint32_t previous = 0;
    bool have_previous = false;
    for (uint32_t column : row)
    {
        if ((uint64_t)column >= limit ||
            (have_previous && column <= previous))
        {
            return false;
        }
        previous = column;
        have_previous = true;
    }
    return true;
}

size_t SymmetricDifferenceCountBelow(
    const std::vector<uint32_t>& a,
    const std::vector<uint32_t>& b,
    uint32_t limit)
{
    size_t ai = 0, bi = 0, count = 0;
    while (ai < a.size() || bi < b.size())
    {
        while (ai < a.size() && a[ai] >= limit) {
            ++ai;
        }
        while (bi < b.size() && b[bi] >= limit) {
            ++bi;
        }
        if (ai >= a.size() && bi >= b.size()) {
            break;
        }
        if (bi >= b.size() || (ai < a.size() && a[ai] < b[bi])) {
            ++ai;
            ++count;
        }
        else if (ai >= a.size() || b[bi] < a[ai]) {
            ++bi;
            ++count;
        }
        else {
            ++ai;
            ++bi;
        }
    }
    return count;
}

bool ValidatePrecodeParams(const PrecodeParams& params)
{
    const uint64_t binary_span =
        (uint64_t)params.BlockCount + params.Staircase + params.DenseRows;
    const uint64_t total_span = binary_span + params.HeavyRows;

    if (params.BlockCount < 2u || params.BlockCount > 64000u ||
        params.Staircase == 0u ||
        params.SourceHits == 0u || params.SourceHits > 8u ||
        params.DenseRows > 64u || params.HeavyRows > 128u ||
        (params.Field != CompletionField::GF256 &&
         params.Field != CompletionField::MixedGF256GF16) ||
        (params.HeavyFamily != HeavyCoefficientFamily::PeriodicCauchy &&
         params.HeavyFamily != HeavyCoefficientFamily::HashedNonzero) ||
        binary_span > UINT16_MAX || total_span > UINT16_MAX)
    {
        return false;
    }

    // Identity-corner flips address both halves of the K + S deck.
    const uint64_t known_span =
        (uint64_t)params.BlockCount + params.Staircase;
    if (params.Field == CompletionField::MixedGF256GF16 &&
        (params.HeavyRows != ActiveMixedGF256Rows() +
                ActiveMixedGF16Rows() ||
         params.HeavyFamily != HeavyCoefficientFamily::PeriodicCauchy))
    {
        return false;
    }
    return !params.DenseIdentityCorner ||
        known_span >= 2u * (uint64_t)(params.DenseRows >> 1);
}

} // namespace

PrecodeParams MakeCertifiedParams(uint32_t block_count, uint64_t seed)
{
    PrecodeParams params;
    params.BlockCount = block_count;
    params.Staircase =
        (block_count >= 2u && block_count <= 64000u) ?
        wirehair::GetDenseCount(block_count) : 0u;
    params.DenseRows = 12u;
    params.HeavyRows = 12u;
    params.SourceHits = CertifiedSourceHits(block_count);
    params.DenseIdentityCorner = false;
    params.Seed = seed;
    return params;
}

PrecodeParams MakeMixedParams(uint32_t block_count, uint64_t seed)
{
    PrecodeParams params = MakeCertifiedParams(block_count, seed);
    params.Field = CompletionField::MixedGF256GF16;
    params.HeavyRows = ActiveMixedGF256Rows() + ActiveMixedGF16Rows();
    return params;
}

bool BuildPrecodeSystem(const PrecodeParams& params, PrecodeSystem& out)
{
    const uint32_t K = params.BlockCount;
    const uint32_t S = params.Staircase;
    const uint32_t D2 = params.DenseRows;
    const uint32_t N1 = params.SourceHits;
    // Reject the complete parameter domain before modifying `out` or
    // allocating row/deck storage.
    if (!ValidatePrecodeParams(params)) {
        return false;
    }
    const uint64_t span_wide = (uint64_t)K + S + D2;
    const uint32_t span = (uint32_t)span_wide;

    out.Params = params;
    out.StaircaseRows.assign(S, std::vector<uint32_t>());
    out.DenseRowColumns.assign(D2, std::vector<uint32_t>());

    wirehair::PCGRandom prng;
    prng.Seed(params.Seed, K);

    // --- Staircase rows ---
    // Reserve for the expected load: K*N1/S source hits + parity + link
    {
        const uint32_t expected = (K * N1) / S + 3u;
        for (uint32_t j = 0; j < S; ++j) {
            out.StaircaseRows[j].reserve(expected);
        }
    }

    // Each source column hits min(N1, S) distinct parities.  Distinctness
    // by retry: the tiny-K dense-count table goes as low as S=2 (K=2), but
    // hits <= S always leaves a free residue, so the loop terminates.
    {
        uint32_t picks[8];
        const uint32_t hits = N1 < S ? N1 : S;
        for (uint32_t c = 0; c < K; ++c)
        {
            for (uint32_t hit = 0; hit < hits; ++hit)
            {
                uint32_t p;
                bool collide;
                do {
                    p = UniformBelow(prng, S);
                    collide = false;
                    for (uint32_t j = 0; j < hit; ++j) {
                        if (picks[j] == p) {
                            collide = true;
                            break;
                        }
                    }
                } while (collide);
                picks[hit] = p;
                out.StaircaseRows[p].push_back(c);
            }
        }
    }

    // Own parity column and staircase link.  Source entries were appended
    // in ascending c, so appending K+j-1 then K+j keeps rows sorted, and
    // no column can appear twice (per-source hits are distinct parities).
    for (uint32_t j = 0; j < S; ++j)
    {
        if (j > 0u) {
            out.StaircaseRows[j].push_back(K + j - 1u);
        }
        out.StaircaseRows[j].push_back(K + j);
    }

    // --- Shuffle-2 dense rows ---
    // Documented rule (experiments/precode/README.md, D2 section):
    //   1. deck = shuffle over span columns; set_count = ceil(span/2);
    //      first row = deck[0 .. set_count).
    //   2. Reshuffle; next floor(D2/2) rows each flip {deck[ii],
    //      deck[set_count + ii]}, ii = 0, 1, ...
    //   3. Reshuffle again; remaining floor(D2/2) - 1 + (D2 & 1) rows by
    //      the same flip rule, ii restarting at 0.
    if (D2 > 0u)
    {
        // Certified construction decks over all binary columns; the
        // identity-corner variant excludes the D2 dense columns and gives
        // each row its own dense column instead (see the header)
        const uint32_t deck_span =
            params.DenseIdentityCorner ? (K + S) : span;
        const uint32_t set_count = (deck_span + 1u) >> 1;
        std::vector<uint16_t> deck(deck_span);
        std::vector<uint8_t> bitmap(deck_span, 0);

        UnbiasedShuffleDeck(prng, deck.data(), deck_span);
        for (uint32_t i = 0; i < set_count; ++i) {
            bitmap[deck[i]] = 1;
        }

        uint32_t row_i = 0;
        const auto emit_row = [&]() {
            std::vector<uint32_t>& columns = out.DenseRowColumns[row_i];
            columns.reserve(set_count + 8u);
            for (uint32_t col = 0; col < deck_span; ++col) {
                if (bitmap[col]) {
                    columns.push_back(col);
                }
            }
            if (params.DenseIdentityCorner) {
                // Own dense column is above every deck column: stays sorted
                columns.push_back(K + S + row_i);
            }
            ++row_i;
        };
        emit_row();

        const uint32_t halves[2] = {
            D2 >> 1,
            (D2 >> 1) + (D2 & 1u) - 1u
        };
        for (uint32_t half = 0; half < 2u; ++half)
        {
            UnbiasedShuffleDeck(prng, deck.data(), deck_span);
            for (uint32_t ii = 0; ii < halves[half]; ++ii)
            {
                // Deck entries at distinct positions are distinct columns,
                // so each flip pair changes exactly two columns
                bitmap[deck[ii]] ^= 1u;
                bitmap[deck[set_count + ii]] ^= 1u;
                emit_row();
            }
        }
    }

    return ValidatePrecodeSystem(out);
}

bool ValidatePrecodeSystem(const PrecodeSystem& system)
{
    const PrecodeParams& params = system.Params;
    const uint32_t K = params.BlockCount;
    const uint32_t S = params.Staircase;
    const uint32_t D2 = params.DenseRows;
    const uint64_t binary_span = (uint64_t)K + S + D2;

    if (!ValidatePrecodeParams(params) ||
        system.StaircaseRows.size() != S ||
        system.DenseRowColumns.size() != D2)
    {
        return false;
    }

    const uint32_t staircase_end = K + S;
    std::vector<uint8_t> source_hits(K, 0u);
    for (uint32_t row_index = 0; row_index < S; ++row_index)
    {
        const std::vector<uint32_t>& row = system.StaircaseRows[row_index];
        if (!IsStrictlyIncreasingBelow(row, staircase_end)) {
            return false;
        }
        const uint32_t own = K + row_index;
        const uint32_t link = own - 1u;
        bool have_own = false;
        bool have_link = row_index == 0u;
        for (uint32_t column : row)
        {
            if (column < K)
            {
                if (source_hits[column] == UINT8_MAX) {
                    return false;
                }
                ++source_hits[column];
            }
            else if (column == own) {
                have_own = true;
            }
            else if (row_index > 0u && column == link) {
                have_link = true;
            }
            else {
                return false;
            }
        }
        if (!have_own || !have_link) {
            return false;
        }
    }
    const uint8_t expected_hits = (uint8_t)std::min(params.SourceHits, S);
    for (uint8_t hits : source_hits) {
        if (hits != expected_hits) {
            return false;
        }
    }

    if (D2 == 0u) {
        return true;
    }

    const uint32_t known_span = staircase_end;
    const size_t first_expected = params.DenseIdentityCorner ?
        (known_span + 1u) / 2u : ((size_t)binary_span + 1u) / 2u;
    for (uint32_t row_index = 0; row_index < D2; ++row_index)
    {
        const std::vector<uint32_t>& row = system.DenseRowColumns[row_index];
        if (!IsStrictlyIncreasingBelow(row, binary_span)) {
            return false;
        }

        size_t known_count = 0;
        size_t dense_count = 0;
        bool have_own = false;
        for (uint32_t column : row)
        {
            if (column < known_span) {
                ++known_count;
            }
            else {
                ++dense_count;
                have_own = have_own || column == known_span + row_index;
            }
        }
        if (params.DenseIdentityCorner &&
            (dense_count != 1u || !have_own))
        {
            return false;
        }
        if (row_index == 0u)
        {
            const size_t first_count = params.DenseIdentityCorner ?
                known_count : row.size();
            if (first_count != first_expected) {
                return false;
            }
        }
        else
        {
            const uint32_t difference_limit = params.DenseIdentityCorner ?
                known_span : (uint32_t)binary_span;
            if (SymmetricDifferenceCountBelow(
                    system.DenseRowColumns[row_index - 1u], row,
                    difference_limit) != 2u)
            {
                return false;
            }
        }
    }
    return true;
}

uint8_t HeavyCoefficient(
    uint32_t heavy_row,
    uint32_t ge_column,
    uint32_t heavy_rows)
{
    // Past 128 rows the X window shrinks below H and the per-window MDS
    // claim silently degrades
    CAT_DEBUG_ASSERT(heavy_rows >= 1u && heavy_rows <= 128u);
    CAT_DEBUG_ASSERT(heavy_row < heavy_rows);

    // Y values occupy [0, H); X values occupy [H, 256).  X ^ Y is never
    // zero because the sets are disjoint, and within one 244-column window
    // all X are distinct, giving the Cauchy MDS property there.
    const uint32_t x = heavy_rows + (ge_column % (256u - heavy_rows));
    const uint32_t y = heavy_row;
    return gf256_inv((uint8_t)(x ^ y));
}

const MixedCoefficientRows* GetMixedCoefficientRows()
{
    if (!InitializeGF16()) {
        return nullptr;
    }
    const auto build_rows = [](MixedCoefficientGeometry geometry) {
        MixedCoefficientRows result = {};
        // The first ten rows, optional Y=11 test row, and shared-X extension
        // rows deliberately keep
        // the frozen H12 X coordinates [12, 256).  H13/H14 append Y values;
        // moving X with H would rewrite every existing coefficient.
        const uint32_t H = kMixedGF256Rows + kMixedGF16Rows;
        for (uint32_t residue = 0;
             residue < kMixedCoefficientPeriod;
             ++residue)
        {
            for (uint32_t row = 0; row < kMixedGF256RowsMax; ++row) {
                // Y=10 develops a singular H15 independently scheduled
                // corner at K=23092.  Y=11 is the other unused subfield
                // coordinate below X=12 and is exhaustively nonsingular for
                // every K=2..64000 under the active keyed P32 construction.
                const uint32_t cauchy_y =
                    row == kMixedGF256Rows ? row + 1u : row;
                result.Subfield[row][residue] =
                    HeavyCoefficient(cauchy_y, residue, H);
            }
            for (uint32_t row = 0; row < kMixedGF16RowsMax; ++row)
            {
                uint16_t coefficient = MixedGF16Coefficient(row, residue);
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
                if (geometry == MixedCoefficientGeometry::SharedCauchyX) {
                    coefficient =
                        MixedGF16SharedXCoefficient(row, residue);
                }
#else
                (void)geometry;
#endif
                result.Extension[row][residue] = coefficient;
            }
        }
        return result;
    };
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    if (ActiveMixedCoefficientGeometry() ==
        MixedCoefficientGeometry::SharedCauchyX)
    {
        static const MixedCoefficientRows shared_rows =
            build_rows(MixedCoefficientGeometry::SharedCauchyX);
        return &shared_rows;
    }
#endif
    static const MixedCoefficientRows frozen_rows =
        build_rows(MixedCoefficientGeometry::FrozenPowerX);
    return &frozen_rows;
}

const MixedPackedCoefficients* GetMixedPackedCoefficients()
{
    const MixedCoefficientRows* rows = GetMixedCoefficientRows();
    if (!rows) {
        return nullptr;
    }
    const auto pack_rows = [](
        const MixedCoefficientRows* source_rows,
        uint32_t subfield_rows) {
        MixedPackedCoefficients result = {};
        for (uint32_t residue = 0;
             residue < kMixedCoefficientPeriod;
             ++residue)
        {
            for (uint32_t row = 0; row < subfield_rows; ++row) {
                result.ByResidue[residue][row >> 2] |=
                    (uint64_t)source_rows->Subfield[row][residue] <<
                    ((row & 3u) * 16u);
            }
            for (uint32_t er = 0; er < kMixedGF16RowsMax; ++er)
            {
                const uint32_t row = subfield_rows + er;
                result.ByResidue[residue][row >> 2] |=
                    (uint64_t)source_rows->Extension[er][residue] <<
                    ((row & 3u) * 16u);
            }
        }
        return result;
    };
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    if (ActiveMixedCoefficientGeometry() ==
        MixedCoefficientGeometry::SharedCauchyX)
    {
        static const MixedPackedCoefficients shared_packed_base =
            pack_rows(rows, kMixedGF256Rows);
        static const MixedPackedCoefficients shared_packed_extra =
            pack_rows(rows, kMixedGF256RowsMax);
        return ActiveMixedGF256Rows() == kMixedGF256Rows ?
            &shared_packed_base : &shared_packed_extra;
    }
#endif
    static const MixedPackedCoefficients frozen_packed =
        pack_rows(rows, kMixedGF256Rows);
    return &frozen_packed;
}

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
static thread_local uint32_t MixedCoefficientPeriodForTesting =
    kMixedCoefficientPeriod;
static thread_local uint32_t MixedResidueSkewForTesting = 0u;
static thread_local MixedResidueSchedule MixedResidueScheduleForTesting =
    MixedResidueSchedule::Constant;
static thread_local uint32_t MixedResidueHashSeedForTesting = 0u;
static thread_local MixedCoefficientGeometry MixedGeometryForTesting =
    MixedCoefficientGeometry::FrozenPowerX;
static thread_local uint32_t MixedGF256RowsForTesting = kMixedGF256Rows;
static thread_local uint32_t MixedGF16RowsForTesting = kMixedGF16Rows;
static thread_local bool MixedIndependentExtensionResiduesForTesting = false;
static thread_local uint32_t MixedExtensionResidueHashSeedForTesting = 0u;
static thread_local uint32_t MixedIndependentExtensionSeedXorForTesting = 78u;

static bool HasFullCycleMixedResidueSeed(
    uint32_t period,
    uint32_t step_count,
    uint32_t seed)
{
    static const uint32_t kStepCycle = 127u;
    uint64_t cycle_sum = 0u;
    for (uint32_t i = 0u; i < kStepCycle; ++i)
    {
        uint32_t x = i + UINT32_C(0x9e3779b9) +
            seed * UINT32_C(0x85ebca6b);
        x = (x ^ (x >> 16)) * UINT32_C(0x85ebca6b);
        x = (x ^ (x >> 13)) * UINT32_C(0xc2b2ae35);
        x ^= x >> 16;
        cycle_sum += 1u + x % step_count;
    }
    uint32_t a = (uint32_t)(cycle_sum % period);
    uint32_t b = period;
    while (b != 0u) {
        const uint32_t remainder = a % b;
        a = b;
        b = remainder;
    }
    return a == 1u;
}

static bool SelectIndependentExtensionResidueSeed(
    uint32_t period,
    uint32_t step_count,
    uint32_t base_seed,
    uint32_t& selected_seed)
{
    uint32_t candidate =
        base_seed ^ MixedIndependentExtensionSeedXorForTesting;
    for (uint32_t attempt = 0u; attempt < 1024u; ++attempt, ++candidate)
    {
        if (candidate != base_seed &&
            HasFullCycleMixedResidueSeed(period, step_count, candidate))
        {
            selected_seed = candidate;
            return true;
        }
    }
    return false;
}

bool SetMixedCoefficientPeriodForTesting(uint32_t period)
{
    const uint32_t H =
        MixedGF256RowsForTesting + MixedGF16RowsForTesting;
    if (period < H || period > kMixedCoefficientPeriod) {
        return false;
    }
    MixedCoefficientPeriodForTesting = period;
    MixedResidueSkewForTesting = 0u;
    MixedResidueScheduleForTesting = MixedResidueSchedule::Constant;
    MixedIndependentExtensionResiduesForTesting = false;
    return true;
}

bool SetMixedResidueSkewForTesting(uint32_t skew)
{
    const uint32_t period = MixedCoefficientPeriodForTesting;
    const uint32_t H =
        MixedGF256RowsForTesting + MixedGF16RowsForTesting;
    if (skew >= period ||
        (skew != 0u &&
         (MixedGeometryForTesting !=
                MixedCoefficientGeometry::SharedCauchyX ||
          skew > period - H ||
          MixedResidueScheduleForTesting !=
                MixedResidueSchedule::Constant)))
    {
        return false;
    }
    MixedResidueSkewForTesting = skew;
    MixedIndependentExtensionResiduesForTesting = false;
    return true;
}

bool SetMixedResidueScheduleForTesting(MixedResidueSchedule schedule)
{
    if (schedule != MixedResidueSchedule::Constant &&
        schedule != MixedResidueSchedule::Ramp &&
        schedule != MixedResidueSchedule::Hashed)
    {
        return false;
    }
    const uint32_t H =
        MixedGF256RowsForTesting + MixedGF16RowsForTesting;
    if (schedule != MixedResidueSchedule::Constant &&
        (MixedGeometryForTesting !=
                MixedCoefficientGeometry::SharedCauchyX ||
         MixedCoefficientPeriodForTesting <= H ||
         MixedResidueSkewForTesting != 0u))
    {
        return false;
    }
    MixedResidueScheduleForTesting = schedule;
    MixedIndependentExtensionResiduesForTesting = false;
    return true;
}

void SetMixedResidueHashSeedForTesting(uint32_t seed)
{
    MixedResidueHashSeedForTesting = seed;
    MixedIndependentExtensionResiduesForTesting = false;
}

bool SelectFullCycleMixedResidueKeyedSeedForTesting(
    uint32_t base_seed,
    uint32_t block_count,
    uint32_t& selected_seed)
{
    uint32_t candidate = base_seed ^
        (block_count + UINT32_C(0x9e3779b9));
    candidate = (candidate ^ (candidate >> 16)) * UINT32_C(0x85ebca6b);
    candidate = (candidate ^ (candidate >> 13)) * UINT32_C(0xc2b2ae35);
    candidate ^= candidate >> 16;

    const uint32_t period = ActiveMixedCoefficientPeriod();
    for (uint32_t attempt = 0u; attempt < 1024u; ++attempt, ++candidate)
    {
        SetMixedResidueHashSeedForTesting(candidate);
        uint32_t a = ActiveMixedResidueBlockShift(127u);
        uint32_t b = period;
        while (b != 0u) {
            const uint32_t remainder = a % b;
            a = b;
            b = remainder;
        }
        if (a == 1u) {
            selected_seed = candidate;
            return true;
        }
    }
    return false;
}

bool SetMixedIndependentExtensionResiduesForTesting(bool enabled)
{
    if (enabled &&
        (MixedGeometryForTesting !=
                MixedCoefficientGeometry::SharedCauchyX ||
         MixedResidueScheduleForTesting != MixedResidueSchedule::Hashed ||
         MixedCoefficientPeriodForTesting <=
            MixedGF256RowsForTesting + MixedGF16RowsForTesting))
    {
        return false;
    }
    if (enabled &&
        !SelectIndependentExtensionResidueSeed(
            MixedCoefficientPeriodForTesting,
            MixedCoefficientPeriodForTesting -
                MixedGF256RowsForTesting - MixedGF16RowsForTesting,
            MixedResidueHashSeedForTesting,
            MixedExtensionResidueHashSeedForTesting))
    {
        return false;
    }
    MixedIndependentExtensionResiduesForTesting = enabled;
    return true;
}

void SetMixedIndependentExtensionSeedXorForTesting(uint32_t seed_xor)
{
    MixedIndependentExtensionSeedXorForTesting = seed_xor;
    MixedIndependentExtensionResiduesForTesting = false;
}

bool SetMixedCoefficientGeometryForTesting(
    MixedCoefficientGeometry geometry)
{
    if (geometry != MixedCoefficientGeometry::FrozenPowerX &&
        geometry != MixedCoefficientGeometry::SharedCauchyX)
    {
        return false;
    }
    if (geometry != MixedCoefficientGeometry::SharedCauchyX &&
        MixedGF256RowsForTesting != kMixedGF256Rows)
    {
        return false;
    }
    MixedGeometryForTesting = geometry;
    if (geometry != MixedCoefficientGeometry::SharedCauchyX) {
        MixedResidueSkewForTesting = 0u;
        MixedResidueScheduleForTesting = MixedResidueSchedule::Constant;
    }
    MixedIndependentExtensionResiduesForTesting = false;
    return true;
}

bool SetMixedGF16RowsForTesting(uint32_t rows)
{
    if (rows < kMixedGF16Rows || rows > kMixedGF16RowsMax ||
        MixedCoefficientPeriodForTesting <
            MixedGF256RowsForTesting + rows)
    {
        return false;
    }
    MixedGF16RowsForTesting = rows;
    MixedResidueSkewForTesting = 0u;
    MixedResidueScheduleForTesting = MixedResidueSchedule::Constant;
    MixedIndependentExtensionResiduesForTesting = false;
    return true;
}

bool SetMixedGF256RowsForTesting(uint32_t rows)
{
    if (rows < kMixedGF256Rows || rows > kMixedGF256RowsMax ||
        MixedCoefficientPeriodForTesting < rows + MixedGF16RowsForTesting ||
        (rows != kMixedGF256Rows &&
         MixedGeometryForTesting != MixedCoefficientGeometry::SharedCauchyX))
    {
        return false;
    }
    MixedGF256RowsForTesting = rows;
    MixedResidueSkewForTesting = 0u;
    MixedResidueScheduleForTesting = MixedResidueSchedule::Constant;
    MixedIndependentExtensionResiduesForTesting = false;
    return true;
}
#endif

uint32_t ActiveMixedCoefficientPeriod()
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    return MixedCoefficientPeriodForTesting;
#else
    return kMixedCoefficientPeriod;
#endif
}

uint32_t ActiveMixedCoefficientResidue(uint32_t column)
{
    const uint32_t period = ActiveMixedCoefficientPeriod();
    return (column % period +
        ActiveMixedResidueBlockShift(column / period)) % period;
}

uint32_t ActiveMixedExtensionCoefficientResidue(uint32_t column)
{
    const uint32_t period = ActiveMixedCoefficientPeriod();
    return (column % period +
        ActiveMixedExtensionResidueBlockShift(column / period)) % period;
}

uint32_t ActiveMixedResidueBlockShift(uint32_t block_index)
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    const uint32_t period = ActiveMixedCoefficientPeriod();
    if (MixedResidueScheduleForTesting == MixedResidueSchedule::Ramp)
    {
        const uint32_t H =
            MixedGF256RowsForTesting + MixedGF16RowsForTesting;
        const uint32_t step_count = period - H;
        if (step_count == 0u) return 0u;
        const uint64_t complete_cycles = block_index / step_count;
        const uint64_t remainder = block_index % step_count;
        // Boundary steps cycle through 2,3,...,P-H,1.  Every step is at most
        // P-H, preserving distinct residues across an H-column boundary,
        // while the cumulative block labels have a much longer period than
        // a single constant rotation.
        const uint64_t cycle_sum =
            (uint64_t)step_count * (step_count + 1u) / 2u;
        const uint64_t prefix_sum = remainder * (remainder + 3u) / 2u;
        return (uint32_t)(
            (complete_cycles * cycle_sum + prefix_sum) % period);
    }
    if (MixedResidueScheduleForTesting == MixedResidueSchedule::Hashed)
    {
        const uint32_t H =
            MixedGF256RowsForTesting + MixedGF16RowsForTesting;
        const uint32_t step_count = period - H;
        if (step_count == 0u) return 0u;
        static const uint32_t kStepCycle = 127u;
        static thread_local uint32_t cached_step_count = 0u;
        static thread_local uint32_t cached_seed = UINT32_MAX;
        static thread_local uint64_t prefix[kStepCycle + 1u] = {};
        if (cached_step_count != step_count ||
            cached_seed != MixedResidueHashSeedForTesting)
        {
            prefix[0] = 0u;
            for (uint32_t i = 0u; i < kStepCycle; ++i)
            {
                uint32_t x = i + UINT32_C(0x9e3779b9) +
                    MixedResidueHashSeedForTesting *
                        UINT32_C(0x85ebca6b);
                x = (x ^ (x >> 16)) * UINT32_C(0x85ebca6b);
                x = (x ^ (x >> 13)) * UINT32_C(0xc2b2ae35);
                x ^= x >> 16;
                prefix[i + 1u] = prefix[i] + 1u + x % step_count;
            }
            cached_step_count = step_count;
            cached_seed = MixedResidueHashSeedForTesting;
        }
        const uint64_t complete_cycles = block_index / kStepCycle;
        const uint32_t remainder = block_index % kStepCycle;
        return (uint32_t)(
            (complete_cycles * prefix[kStepCycle] + prefix[remainder]) %
            period);
    }
    return (uint32_t)(
        (uint64_t)block_index * MixedResidueSkewForTesting % period);
#else
    (void)block_index;
    return 0u;
#endif
}

uint32_t ActiveMixedExtensionResidueBlockShift(uint32_t block_index)
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    if (!MixedIndependentExtensionResiduesForTesting) {
        return ActiveMixedResidueBlockShift(block_index);
    }
    const uint32_t period = MixedCoefficientPeriodForTesting;
    const uint32_t H =
        MixedGF256RowsForTesting + MixedGF16RowsForTesting;
    const uint32_t step_count = period - H;
    if (step_count == 0u) return 0u;
    static const uint32_t kStepCycle = 127u;
    static thread_local uint32_t cached_period = 0u;
    static thread_local uint32_t cached_step_count = 0u;
    static thread_local uint32_t cached_seed = UINT32_MAX;
    static thread_local uint64_t prefix[kStepCycle + 1u] = {};
    if (cached_period != period ||
        cached_step_count != step_count ||
        cached_seed != MixedExtensionResidueHashSeedForTesting)
    {
        const uint32_t seed = MixedExtensionResidueHashSeedForTesting;
        prefix[0] = 0u;
        for (uint32_t i = 0u; i < kStepCycle; ++i)
        {
            uint32_t x = i + UINT32_C(0x9e3779b9) +
                seed * UINT32_C(0x85ebca6b);
            x = (x ^ (x >> 16)) * UINT32_C(0x85ebca6b);
            x = (x ^ (x >> 13)) * UINT32_C(0xc2b2ae35);
            x ^= x >> 16;
            prefix[i + 1u] = prefix[i] + 1u + x % step_count;
        }
        cached_period = period;
        cached_step_count = step_count;
        cached_seed = seed;
    }
    const uint64_t complete_cycles = block_index / kStepCycle;
    const uint32_t remainder = block_index % kStepCycle;
    return (uint32_t)(
        (complete_cycles * prefix[kStepCycle] + prefix[remainder]) %
        period);
#else
    return ActiveMixedResidueBlockShift(block_index);
#endif
}

uint32_t ActiveMixedResidueSkew()
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    return MixedResidueSkewForTesting;
#else
    return 0u;
#endif
}

MixedResidueSchedule ActiveMixedResidueSchedule()
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    return MixedResidueScheduleForTesting;
#else
    return MixedResidueSchedule::Constant;
#endif
}

uint32_t ActiveMixedResidueHashSeed()
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    return MixedResidueHashSeedForTesting;
#else
    return 0u;
#endif
}

bool ActiveMixedResiduesRotated()
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    return MixedResidueSkewForTesting != 0u ||
        MixedResidueScheduleForTesting != MixedResidueSchedule::Constant;
#else
    return false;
#endif
}

bool ActiveMixedIndependentExtensionResidues()
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    return MixedIndependentExtensionResiduesForTesting;
#else
    return false;
#endif
}

MixedCoefficientGeometry ActiveMixedCoefficientGeometry()
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    return MixedGeometryForTesting;
#else
    return MixedCoefficientGeometry::FrozenPowerX;
#endif
}

uint32_t ActiveMixedGF16Rows()
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    return MixedGF16RowsForTesting;
#else
    return kMixedGF16Rows;
#endif
}

uint32_t ActiveMixedGF256Rows()
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    return MixedGF256RowsForTesting;
#else
    return kMixedGF256Rows;
#endif
}

uint32_t ActiveMixedPackedCoefficientWords()
{
    return (ActiveMixedGF256Rows() + ActiveMixedGF16Rows() + 3u) / 4u;
}

uint8_t HeavyCoefficientForParams(
    const PrecodeParams& params,
    uint32_t heavy_row,
    uint32_t ge_column)
{
    if (params.HeavyFamily == HeavyCoefficientFamily::PeriodicCauchy) {
        return HeavyCoefficient(heavy_row, ge_column, params.HeavyRows);
    }

    // SplitMix64 finalizer over the complete column and row ids.  Mapping zero
    // to one keeps this comparable to the nonzero Cauchy coefficients while
    // removing their 244-column periodicity entirely.
    uint64_t x = params.Seed ^
        ((uint64_t)ge_column * UINT64_C(0x9e3779b97f4a7c15)) ^
        ((uint64_t)heavy_row * UINT64_C(0xd6e8feb86659fd93)) ^
        UINT64_C(0x6865617679686173);
    x ^= x >> 30;
    x *= UINT64_C(0xbf58476d1ce4e5b9);
    x ^= x >> 27;
    x *= UINT64_C(0x94d049bb133111eb);
    x ^= x >> 31;
    const uint8_t coefficient = (uint8_t)x;
    return coefficient != 0u ? coefficient : 1u;
}

} // namespace wirehair_v2
