#include "WirehairV2Precode.h"

#include "../WirehairTools.h"
#include "../gf256.h"

#include <algorithm>
#include <new>
#include <stdexcept>

namespace wirehair_v2 {
namespace {

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
thread_local test::PrecodeBuildAllocationFailurePoint
    ActivePrecodeBuildAllocationFailurePoint =
        test::PrecodeBuildAllocationFailurePoint::None;
thread_local test::PrecodeBuildAllocationFailureException
    ActivePrecodeBuildAllocationFailureException =
        test::PrecodeBuildAllocationFailureException::BadAlloc;
thread_local uint32_t ActivePrecodeBuildAllocationFailureHits = 0u;
#endif

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
    const PrecodeRowView& row,
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

bool IsKnownDenseAnchorLayout(DenseAnchorLayout layout)
{
    return layout == DenseAnchorLayout::Disabled ||
        layout == DenseAnchorLayout::Two07 ||
        layout == DenseAnchorLayout::Four0369;
}

bool IsAdditionalDenseAnchor(
    DenseAnchorLayout layout,
    uint32_t row_index)
{
    switch (layout)
    {
    case DenseAnchorLayout::Disabled:
        return false;
    case DenseAnchorLayout::Two07:
        return row_index == 7u;
    case DenseAnchorLayout::Four0369:
        return row_index == 3u || row_index == 6u || row_index == 9u;
    }
    return false;
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
        (params.HeavyFamily != HeavyCoefficientFamily::PeriodicCauchy &&
         params.HeavyFamily != HeavyCoefficientFamily::HashedNonzero) ||
        !IsKnownDenseAnchorLayout(params.DenseAnchors) ||
        binary_span > UINT16_MAX || total_span > UINT16_MAX)
    {
        return false;
    }

    // Identity-corner flips address both halves of the K + S deck.
    const uint64_t known_span =
        (uint64_t)params.BlockCount + params.Staircase;
    if (params.DenseAnchors != DenseAnchorLayout::Disabled)
    {
#if !defined(WIREHAIR_V2_ENABLE_TEST_HOOKS) && \
    !defined(WIREHAIR_V2_ENABLE_BENCHMARK_EQUATIONS)
        return false;
#else
        if (params.DenseRows != 12u || params.HeavyRows != 12u ||
            params.DenseIdentityCorner ||
            params.HeavyFamily != HeavyCoefficientFamily::PeriodicCauchy)
        {
            return false;
        }
#endif
    }
    return !params.DenseIdentityCorner ||
        known_span >= 2u * (uint64_t)(params.DenseRows >> 1);
}

} // namespace

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
namespace test {

void SetPrecodeBuildAllocationFailurePointForTesting(
    PrecodeBuildAllocationFailurePoint point,
    PrecodeBuildAllocationFailureException exception)
{
    ActivePrecodeBuildAllocationFailurePoint = point;
    ActivePrecodeBuildAllocationFailureException = exception;
    ActivePrecodeBuildAllocationFailureHits = 0u;
}

uint32_t PrecodeBuildAllocationFailureHitsForTesting()
{
    return ActivePrecodeBuildAllocationFailureHits;
}

void TriggerPrecodeBuildAllocationFailureForTesting(
    PrecodeBuildAllocationFailurePoint point)
{
    if (ActivePrecodeBuildAllocationFailurePoint != point) {
        return;
    }
    ++ActivePrecodeBuildAllocationFailureHits;
    if (ActivePrecodeBuildAllocationFailureException ==
        PrecodeBuildAllocationFailureException::LengthError)
    {
        throw std::length_error(
            "injected WH2 precode-build allocation failure");
    }
    throw std::bad_alloc();
}

} // namespace test
#endif

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
    params.DenseAnchors = DenseAnchorLayout::Disabled;
    params.Seed = seed;
    return params;
}

bool HasValidPrecodeSystemShape(const PrecodeSystem& system)
{
    if (!ValidatePrecodeParams(system.Params)) {
        return false;
    }
    const uint64_t row_count = (uint64_t)system.Params.Staircase +
        system.Params.DenseRows;
    if (row_count + 1u != system.BinaryRowOffsets.size() ||
        system.BinaryRowOffsets.empty() ||
        system.BinaryRowOffsets[0] != 0u ||
        system.BinaryRowColumns.size() > UINT32_MAX)
    {
        return false;
    }
    uint32_t previous = 0u;
    for (uint32_t offset : system.BinaryRowOffsets)
    {
        if (offset < previous || offset > system.BinaryRowColumns.size()) {
            return false;
        }
        previous = offset;
    }
    return previous == system.BinaryRowColumns.size();
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

    // Construct transactionally.  Status-bearing callers translate allocation
    // failures to Wirehair_OOM, and a direct caller that observes an exception
    // must not receive a partially converted CSR graph.
    PrecodeSystem candidate;
    candidate.Params = params;
    const uint32_t row_count = S + D2;
    candidate.BinaryRowOffsets.assign((size_t)row_count + 1u, 0u);

    wirehair::PCGRandom prng;
    prng.Seed(params.Seed, K);

    // --- Staircase rows ---
    // Each source column hits min(N1, S) distinct parities.  Distinctness
    // by retry: the tiny-K dense-count table goes as low as S=2 (K=2), but
    // hits <= S always leaves a free residue, so the loop terminates.  Save
    // only the selected row indices: source columns are implicit in the flat
    // position and are scattered in ascending order after prefixing counts.
    const uint32_t hits = N1 < S ? N1 : S;
    const size_t hit_assignment_count = (size_t)K * hits;
    std::vector<uint16_t> hit_rows(hit_assignment_count);
    {
        uint32_t picks[8];
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
                hit_rows[(size_t)c * hits + hit] = (uint16_t)p;
                ++candidate.BinaryRowOffsets[p + 1u];
            }
        }
    }

    // Account for own parity and staircase-link columns.  Source entries will
    // be scattered in ascending c, so writing K+j-1 then K+j keeps each row
    // sorted and unique.
    for (uint32_t j = 0; j < S; ++j)
    {
        candidate.BinaryRowOffsets[j + 1u] += j == 0u ? 1u : 2u;
    }

    // Dense row lengths are known before their contents.  Reserve their exact
    // slices in the same CSR allocation as the staircase rows.
    const uint32_t deck_span = params.DenseIdentityCorner ? (K + S) : span;
    const uint32_t set_count = D2 == 0u ? 0u : (deck_span + 1u) >> 1;
    for (uint32_t row = 0; row < D2; ++row)
    {
        const bool is_anchor = row == 0u ||
            IsAdditionalDenseAnchor(params.DenseAnchors, row);
        const uint32_t length = is_anchor ?
            set_count + (params.DenseIdentityCorner ? 1u : 0u) :
            (params.DenseIdentityCorner ? 4u : 2u);
        candidate.BinaryRowOffsets[S + row + 1u] = length;
    }

    uint64_t reference_count = 0u;
    for (uint32_t row = 0; row < row_count; ++row)
    {
        reference_count += candidate.BinaryRowOffsets[row + 1u];
        if (reference_count > UINT32_MAX) {
            return false;
        }
        candidate.BinaryRowOffsets[row + 1u] = (uint32_t)reference_count;
    }
    candidate.BinaryRowColumns.resize((size_t)reference_count);

    // Temporarily use each staircase start offset as its fill cursor.  Once all
    // sources are scattered, reconstruct the starts while appending the fixed
    // link/own suffixes.  This avoids a fourth heap allocation for cursors.
    for (uint32_t c = 0; c < K; ++c)
    {
        for (uint32_t hit = 0; hit < hits; ++hit)
        {
            const uint32_t row = hit_rows[(size_t)c * hits + hit];
            candidate.BinaryRowColumns[
                candidate.BinaryRowOffsets[row]++] = c;
        }
    }
    uint32_t next_start = 0u;
    for (uint32_t row = 0; row < S; ++row)
    {
        uint32_t write = candidate.BinaryRowOffsets[row];
        candidate.BinaryRowOffsets[row] = next_start;
        if (row > 0u) {
            candidate.BinaryRowColumns[write++] = K + row - 1u;
        }
        candidate.BinaryRowColumns[write++] = K + row;
        next_start = write;
    }
    if (candidate.BinaryRowOffsets[S] != next_start) {
        return false;
    }

    // --- Shuffle-2 dense rows ---
    // Documented rule (experiments/precode/README.md, D2 section):
    //   1. deck = shuffle over span columns; set_count = ceil(span/2);
    //      first row = deck[0 .. set_count).
    //   2. Reshuffle; next floor(D2/2) rows each flip {deck[ii],
    //      deck[set_count + ii]}, ii = 0, 1, ...
    //   3. Reshuffle again; remaining floor(D2/2) - 1 + (D2 & 1) rows by
    //      the same flip rule, ii restarting at 0.
    // Materialize that system directly in its row-equivalent sparse basis:
    // anchors retain the balanced set-half, and every flip becomes its
    // two-column adjacent-row delta.  The shuffle calls and deck indices stay
    // identical to the documented construction; only redundant original-row
    // materialization is skipped.
    if (D2 > 0u)
    {
        // Keep the small-K shuffle scratch inline.  Both arrays are fully
        // initialized below, so the large-K vectors remain simple fallbacks
        // and their initialization cannot affect the generated equations.
        static const uint32_t kInlineDeckSpan = 512u;
        uint16_t deck_inline[kInlineDeckSpan];
        uint8_t anchor_flags_inline[kInlineDeckSpan];
        std::vector<uint16_t> deck_heap;
        std::vector<uint8_t> anchor_flags_heap;
        uint16_t* deck = deck_inline;
        uint8_t* anchor_flags = anchor_flags_inline;
        if (deck_span > kInlineDeckSpan)
        {
            deck_heap.resize(deck_span);
            anchor_flags_heap.resize(deck_span);
            deck = deck_heap.data();
            anchor_flags = anchor_flags_heap.data();
        }
        std::fill(anchor_flags, anchor_flags + deck_span, uint8_t{0});
        // One bit per anchor lets every balanced row be populated by the same
        // ascending column scan.  Four0369 is the largest accepted layout.
        uint32_t anchor_rows[4] = {};
        uint8_t anchor_bits[4] = {};
        uint32_t anchor_write[4] = {};
        uint32_t anchor_count = 0u;
        uint32_t row_i = 0;
        const auto mark_anchor = [&]() {
            CAT_DEBUG_ASSERT(anchor_count < 4u);
            const uint8_t anchor_bit = (uint8_t)(1u << anchor_count);
            for (uint32_t i = 0; i < set_count; ++i) {
                // Anchors overlap heavily; assignment would erase membership
                // recorded for an earlier anchor at the same column.
                anchor_flags[deck[i]] |= anchor_bit;
            }
            anchor_rows[anchor_count] = row_i;
            anchor_bits[anchor_count] = anchor_bit;
            anchor_write[anchor_count] =
                candidate.BinaryRowOffsets[S + row_i];
            ++anchor_count;
            ++row_i;
        };
        const auto emit_delta = [&](uint32_t first, uint32_t second) {
            uint32_t columns[4];
            uint32_t count = 0u;
            columns[count++] = first;
            columns[count++] = second;
            if (params.DenseIdentityCorner)
            {
                columns[count++] = K + S + row_i - 1u;
                columns[count++] = K + S + row_i;
            }
            std::sort(columns, columns + count);
            const uint32_t write = candidate.BinaryRowOffsets[S + row_i];
            std::copy(
                columns,
                columns + count,
                candidate.BinaryRowColumns.begin() + write);
            CAT_DEBUG_ASSERT(
                write + count ==
                candidate.BinaryRowOffsets[S + row_i + 1u]);
            ++row_i;
        };

        UnbiasedShuffleDeck(prng, deck, deck_span);
        mark_anchor();

        if (params.DenseAnchors != DenseAnchorLayout::Disabled)
        {
            // Row zero uses the initial balanced deck.  Each segment gets a
            // fresh deck: an additional anchor emits its balanced set-half,
            // while other rows consume one distinct set/clear flip pair.
            UnbiasedShuffleDeck(prng, deck, deck_span);
            uint32_t flip_index = 0u;
            while (row_i < D2)
            {
                if (IsAdditionalDenseAnchor(params.DenseAnchors, row_i))
                {
                    UnbiasedShuffleDeck(prng, deck, deck_span);
                    flip_index = 0u;
                    mark_anchor();
                    continue;
                }
                emit_delta(
                    deck[flip_index], deck[set_count + flip_index]);
                ++flip_index;
            }
        }
        else
        {
            const uint32_t halves[2] = {
                D2 >> 1,
                (D2 >> 1) + (D2 & 1u) - 1u
            };
            for (uint32_t half = 0; half < 2u; ++half)
            {
                UnbiasedShuffleDeck(prng, deck, deck_span);
                for (uint32_t ii = 0; ii < halves[half]; ++ii)
                {
                    // Deck entries at distinct positions are distinct
                    // columns, so each flip pair changes exactly two columns.
                    emit_delta(deck[ii], deck[set_count + ii]);
                }
            }
        }

        // Populate every anchor in one ordered pass.  Delta rows were already
        // emitted above and consume no membership scan.
        for (uint32_t col = 0; col < deck_span; ++col)
        {
            const uint8_t flags = anchor_flags[col];
            for (uint32_t anchor = 0; anchor < anchor_count; ++anchor) {
                if ((flags & anchor_bits[anchor]) != 0u) {
                    candidate.BinaryRowColumns[anchor_write[anchor]++] = col;
                }
            }
        }
        if (params.DenseIdentityCorner)
        {
            // Owned dense columns follow the known-column deck and preserve
            // sorted order.  Accepted identity-corner layouts have one
            // anchor, but retaining the loop keeps the basis rule explicit.
            for (uint32_t anchor = 0; anchor < anchor_count; ++anchor) {
                const uint32_t row = anchor_rows[anchor];
                candidate.BinaryRowColumns[anchor_write[anchor]++] =
                    K + S + row;
            }
        }
        for (uint32_t anchor = 0; anchor < anchor_count; ++anchor)
        {
            if (anchor_write[anchor] !=
                candidate.BinaryRowOffsets[
                    S + anchor_rows[anchor] + 1u])
            {
                return false;
            }
        }
        if (row_i != D2) {
            return false;
        }
    }

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    test::TriggerPrecodeBuildAllocationFailureForTesting(
        test::PrecodeBuildAllocationFailurePoint::FinalValidation);
#endif
    if (!ValidatePrecodeSystem(candidate)) {
        return false;
    }
    out.Params = candidate.Params;
    out.BinaryRowOffsets.swap(candidate.BinaryRowOffsets);
    out.BinaryRowColumns.swap(candidate.BinaryRowColumns);
    return true;
}

bool ValidatePrecodeSystem(const PrecodeSystem& system)
{
    const PrecodeParams& params = system.Params;
    const uint32_t K = params.BlockCount;
    const uint32_t S = params.Staircase;
    const uint32_t D2 = params.DenseRows;
    const uint64_t binary_span = (uint64_t)K + S + D2;

    if (!HasValidPrecodeSystemShape(system))
    {
        return false;
    }

    const uint32_t staircase_end = K + S;
    // Validation is part of every generated system build.  Avoid a heap
    // allocation for the small-K range while retaining the same bounded
    // heap fallback for larger block counts.
    static const uint32_t kInlineSourceHits = 1024u;
    uint8_t source_hits_inline[kInlineSourceHits];
    std::vector<uint8_t> source_hits_heap;
    uint8_t* source_hits = source_hits_inline;
    if (K > kInlineSourceHits) {
        source_hits_heap.assign(K, uint8_t{0});
        source_hits = source_hits_heap.data();
    }
    else {
        std::fill(source_hits, source_hits + K, uint8_t{0});
    }
    for (uint32_t row_index = 0; row_index < S; ++row_index)
    {
        const PrecodeRowView row = system.StaircaseRow(row_index);
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
    for (uint32_t column = 0; column < K; ++column) {
        if (source_hits[column] != expected_hits) {
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
        const PrecodeRowView basis = system.DenseBasisRow(row_index);
        if (!IsStrictlyIncreasingBelow(basis, binary_span))
        {
            return false;
        }

        size_t known_count = 0;
        size_t dense_count = 0;
        bool have_previous = false;
        bool have_current = false;
        for (uint32_t column : basis)
        {
            if (column < known_span) {
                ++known_count;
            }
            else {
                ++dense_count;
                have_previous = have_previous ||
                    (row_index > 0u &&
                     column == known_span + row_index - 1u);
                have_current = have_current ||
                    column == known_span + row_index;
            }
        }
        const bool is_anchor = row_index == 0u ||
            IsAdditionalDenseAnchor(params.DenseAnchors, row_index);
        if (is_anchor)
        {
            if (params.DenseIdentityCorner) {
                if (known_count != first_expected || dense_count != 1u ||
                    !have_current)
                {
                    return false;
                }
            }
            else if (basis.size() != first_expected) {
                return false;
            }
        }
        else
        {
            const size_t expected_basis_size =
                params.DenseIdentityCorner ? 4u : 2u;
            if (basis.size() != expected_basis_size)
            {
                return false;
            }
            if (params.DenseIdentityCorner &&
                (known_count != 2u || dense_count != 2u ||
                 !have_previous || !have_current))
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
