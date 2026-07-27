#include "WirehairV2Precode.h"

#include "../WirehairEnvironment.h"
#include "../WirehairTools.h"
#include "../gf256.h"

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <vector>

namespace wirehair_v2 {

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
namespace {

bool IsAsciiSpace(char c)
{
    return c == ' ' || c == '\t' || c == '\r' ||
        c == '\n' || c == '\f' || c == '\v';
}

void SkipAsciiSpace(const char*& p)
{
    while (IsAsciiSpace(*p)) {
        ++p;
    }
}

} // namespace

bool ParseStaircaseDegreesForTesting(
    const char* text,
    std::vector<double>& weights)
{
    std::vector<double> parsed;
    if (!text) {
        return false;
    }
    const char* p = text;
    SkipAsciiSpace(p);
    if (*p == '\0') {
        return false;
    }
    for (;;)
    {
        errno = 0;
        char* end = nullptr;
        const double weight = std::strtod(p, &end);
        if (end == p || errno == ERANGE || !std::isfinite(weight) ||
            weight < 0.0)
        {
            return false;
        }
        parsed.push_back(weight);
        p = end;
        SkipAsciiSpace(p);
        if (*p == '\0') {
            break;
        }
        if (*p != ',') {
            return false;
        }
        ++p;
        SkipAsciiSpace(p);
        if (*p == '\0') {
            return false;
        }
    }
    weights.swap(parsed);
    return true;
}

bool ParseStaircaseRowDegreesForTesting(
    const char* text,
    std::vector<uint32_t>& degrees)
{
    std::vector<uint32_t> parsed;
    if (!text) {
        return false;
    }
    const char* p = text;
    SkipAsciiSpace(p);
    if (*p == '\0') {
        return false;
    }
    for (;;)
    {
        // strtoull accepts signs and silently wraps a leading minus.  Row
        // degrees are unsigned decimal values, so require the first digit.
        if (*p < '0' || *p > '9') {
            return false;
        }
        errno = 0;
        char* end = nullptr;
        const unsigned long long degree = std::strtoull(p, &end, 10);
        if (end == p || errno == ERANGE || degree > UINT32_MAX) {
            return false;
        }
        parsed.push_back((uint32_t)degree);
        p = end;
        SkipAsciiSpace(p);
        if (*p == '\0') {
            break;
        }
        if (*p != ',') {
            return false;
        }
        ++p;
        SkipAsciiSpace(p);
        if (*p == '\0') {
            return false;
        }
    }
    degrees.swap(parsed);
    return true;
}

bool ParseStaircaseDegreeScaleForTesting(
    const char* text,
    double& scale)
{
    if (!text) {
        return false;
    }
    const char* p = text;
    SkipAsciiSpace(p);
    errno = 0;
    char* end = nullptr;
    const double parsed = std::strtod(p, &end);
    if (end == p || errno == ERANGE || !std::isfinite(parsed)) {
        return false;
    }
    const char* tail = end;
    SkipAsciiSpace(tail);
    if (*tail != '\0' ||
        !(parsed >= kStaircaseDegreeScaleMin) ||
        !(parsed <= kStaircaseDegreeScaleMax))
    {
        return false;
    }
    scale = parsed;
    return true;
}

/// Staircase row-degree distribution for the shaped-construction experiment,
/// read once from WIREHAIR_V2_STAIRCASE_DEGREES as comma-separated weights for
/// degrees 1, 2, 3, ...  Unset (the default) keeps the stock construction.
/// Per-thread override of the staircase degree distribution.
///
/// The environment variable below is parsed ONCE per process, which makes it
/// useless to a search: every candidate in a run would share one distribution.
/// This is a real optimization target because redistributing the same edge
/// budget changes both inactivation and XOR work.  No quantitative speed claim
/// is attached here: the inherited operation-cost model is not reproducible
/// and remains explicitly opt-in until it is recalibrated.
/// Thread-local, matching the convention of the other *ForTesting knobs, so
/// worker threads can each evaluate a different shape concurrently.
///
/// These two definitions must stay at wirehair_v2 scope, NOT in the anonymous
/// namespace below: the header declares them at wirehair_v2 scope, so an
/// anonymous-namespace definition is a different, internal-linkage function
/// that leaves the declared symbol undefined (it was silently discarded as an
/// unused static, and essearch could not link against it).
static thread_local std::vector<double> g_staircase_degrees_override;

void SetStaircaseDegreesForTesting(const std::vector<double>& weights)
{
    g_staircase_degrees_override = weights;
}

void ClearStaircaseDegreesForTesting()
{
    g_staircase_degrees_override.clear();
}

/// Per-thread override of the target mean staircase row degree.  The exact
/// kStaircaseDegreeScaleUnset value means unset; other negative values are
/// invalid rather than aliases for the sentinel.  Zero cannot be the sentinel
/// because zero is a legal target mean row degree (an empty edge budget).
///
/// Same scope rule as the two definitions above: these MUST stay at
/// wirehair_v2 scope, NOT in the anonymous namespace below, or the symbol the
/// header declares is left undefined and essearch fails to link.
static thread_local double g_staircase_degree_scale_override =
    kStaircaseDegreeScaleUnset;

void SetStaircaseDegreeScaleForTesting(double scale)
{
    // Stored VERBATIM.  A NaN is not negative, so it counts as supplied and
    // MakeStaircaseDegreeMixture refuses it; silently reinterpreting it as
    // "unset" would turn a broken caller into a legacy-budget measurement.
    g_staircase_degree_scale_override = scale;
}

void ClearStaircaseDegreeScaleForTesting()
{
    g_staircase_degree_scale_override = kStaircaseDegreeScaleUnset;
}

/// Per-thread override of the REALIZED row-degree sequence -- see the header.
/// Same scope rule as the two overrides above: wirehair_v2 scope, not the
/// anonymous namespace, or the declared symbol is left undefined.
static thread_local std::vector<uint32_t> g_staircase_row_degrees_override;

void SetStaircaseRowDegreesForTesting(const std::vector<uint32_t>& degrees)
{
    g_staircase_row_degrees_override = degrees;
}

void ClearStaircaseRowDegreesForTesting()
{
    g_staircase_row_degrees_override.clear();
}
#endif

namespace {

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
template <typename T>
struct ParsedListEnvironment
{
    bool IsSet = false;
    bool Valid = true;
    std::vector<T> Values;
};

struct ParsedScaleEnvironment
{
    bool IsSet = false;
    bool Valid = true;
    double Value = kStaircaseDegreeScaleUnset;
};

static const ParsedListEnvironment<double>& StaircaseDegreesEnvironment()
{
    static const ParsedListEnvironment<double> parsed = [] {
        ParsedListEnvironment<double> out;
        const wirehair::EnvironmentValue value(
            "WIREHAIR_V2_STAIRCASE_DEGREES");
        if (!value.IsSet()) {
            return out;
        }
        out.IsSet = true;
        out.Valid = ParseStaircaseDegreesForTesting(
            value.Get(), out.Values);
        return out;
    }();
    return parsed;
}

static const ParsedListEnvironment<uint32_t>&
StaircaseRowDegreesEnvironment()
{
    static const ParsedListEnvironment<uint32_t> parsed = [] {
        ParsedListEnvironment<uint32_t> out;
        const wirehair::EnvironmentValue value(
            "WIREHAIR_V2_STAIRCASE_ROW_DEGREES");
        if (!value.IsSet()) {
            return out;
        }
        out.IsSet = true;
        out.Valid = ParseStaircaseRowDegreesForTesting(
            value.Get(), out.Values);
        return out;
    }();
    return parsed;
}

static const ParsedScaleEnvironment& StaircaseDegreeScaleEnvironment()
{
    static const ParsedScaleEnvironment parsed = [] {
        ParsedScaleEnvironment out;
        const wirehair::EnvironmentValue value(
            "WIREHAIR_V2_STAIRCASE_DEGREE_SCALE");
        if (!value.IsSet()) {
            return out;
        }
        out.IsSet = true;
        out.Valid = ParseStaircaseDegreeScaleForTesting(
            value.Get(), out.Value);
        return out;
    }();
    return parsed;
}

static const std::vector<double>* ActiveStaircaseDegreeDistribution()
{
    if (!g_staircase_degrees_override.empty()) {
        return &g_staircase_degrees_override;
    }
    const ParsedListEnvironment<double>& parsed =
        StaircaseDegreesEnvironment();
    return parsed.IsSet && parsed.Valid ? &parsed.Values : nullptr;
}

/// Realized row-degree sequence in force, or null when none is set.
///
/// Precedence matches the shape knob: the thread-local the search uses, then
/// WIREHAIR_V2_STAIRCASE_ROW_DEGREES parsed once per process.  A sequence wins
/// over a shape, because a shape only SAMPLES a sequence and this names one.
static const std::vector<uint32_t>* ActiveStaircaseRowDegrees()
{
    if (!g_staircase_row_degrees_override.empty()) {
        return &g_staircase_row_degrees_override;
    }
    const ParsedListEnvironment<uint32_t>& parsed =
        StaircaseRowDegreesEnvironment();
    return parsed.IsSet && parsed.Valid ? &parsed.Values : nullptr;
}

/// Validate only the environment arm that wins the documented precedence.
/// A malformed active value is a refused construction, never a partial list
/// or a silent fallback to a different experiment.
static bool StaircaseHookConfigurationValid()
{
    bool pinned = !g_staircase_row_degrees_override.empty();
    if (!pinned)
    {
        const ParsedListEnvironment<uint32_t>& row_degrees =
            StaircaseRowDegreesEnvironment();
        if (row_degrees.IsSet)
        {
            if (!row_degrees.Valid) {
                return false;
            }
            pinned = true;
        }
    }
    if (!pinned && g_staircase_degrees_override.empty())
    {
        const ParsedListEnvironment<double>& degrees =
            StaircaseDegreesEnvironment();
        if (degrees.IsSet && !degrees.Valid) {
            return false;
        }
    }
    if (g_staircase_degree_scale_override ==
        kStaircaseDegreeScaleUnset)
    {
        const ParsedScaleEnvironment& scale =
            StaircaseDegreeScaleEnvironment();
        if (scale.IsSet && !scale.Valid) {
            return false;
        }
    }
    return true;
}

/// True when the shaped-degree experiment is active.
static bool ShapedStaircaseActive()
{
    return ActiveStaircaseDegreeDistribution() != nullptr ||
        ActiveStaircaseRowDegrees() != nullptr;
}

static bool ActiveStaircaseHookValuesValid(
    const PrecodeParams& params,
    const StaircaseDegreeMixture& mixture)
{
    const std::vector<uint32_t>* pinned = ActiveStaircaseRowDegrees();
    if (pinned)
    {
        if (pinned->size() != (size_t)params.Staircase) {
            return false;
        }
        uint64_t total = 0u;
        for (uint32_t degree : *pinned)
        {
            if (degree > params.BlockCount) {
                return false;
            }
            total += degree;
        }
        return total == mixture.EdgeCount;
    }

    const std::vector<double>* shape =
        ActiveStaircaseDegreeDistribution();
    if (!shape) {
        return true;
    }
    double total = 0.0;
    for (double weight : *shape)
    {
        if (!std::isfinite(weight) || weight < 0.0) {
            return false;
        }
        total += weight;
        if (!std::isfinite(total)) {
            return false;
        }
    }
    return total > 0.0;
}

/// Effective target mean staircase row degree for a build.
///
/// Precedence, matching the shaped-degree knob: the thread-local override the
/// search uses, then WIREHAIR_V2_STAIRCASE_DEGREE_SCALE (parsed once per
/// process, which is why it is useless to a search), then the params field.
/// An unparseable or out-of-domain active environment value makes the build
/// fail rather than silently selecting the params field instead.
///
/// Every "not supplied" state uses the exact
/// kStaircaseDegreeScaleUnset sentinel, because zero is a legal target mean
/// row degree and other negative values must remain invalid.
static double ActiveStaircaseDegreeScale(double param_scale)
{
    // Compare only with the exact sentinel.  Thus -2, -infinity and NaN all
    // reach the mixture validator and fail rather than falling through to an
    // environment value or the legacy budget.
    if (g_staircase_degree_scale_override !=
        kStaircaseDegreeScaleUnset)
    {
        return g_staircase_degree_scale_override;
    }
    const ParsedScaleEnvironment& parsed =
        StaircaseDegreeScaleEnvironment();
    return parsed.IsSet && parsed.Valid ? parsed.Value : param_scale;
}
#else
/// Production builds carry no shaped-degree experiment.
static bool ShapedStaircaseActive() { return false; }
static bool StaircaseHookConfigurationValid() { return true; }
static bool ActiveStaircaseHookValuesValid(
    const PrecodeParams&,
    const StaircaseDegreeMixture&)
{
    return true;
}
/// ... and no scale override: the params field is the whole story.
static double ActiveStaircaseDegreeScale(double param_scale)
{
    return param_scale;
}
#define WIREHAIR_V2_SHAPED_FALLBACK_DEFINED 1
#endif


uint32_t CertifiedSourceHits(uint32_t block_count)
{
    // N1 = 1 over 8 <= K <= 100 measured ~5% faster at equal-or-better
    // untuned fail@+0 on top of the small-band staircase rule (K=8
    // 2.22% -> 1.94%, K=12 1.74% -> 1.46%, noise through K=100), but it
    // moves the construction-seed attempt counts that
    // PrecodeSeedSelectionTest pins and removes the deficient fixture that
    // PrecodeRoundTripTest's cold-retry case searches for.  Regenerating
    // those two fixtures is a separate, deliberate change.
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

// Uniform Fisher-Yates permutation for matching-only random choices.  The
// certified dense-row construction above deliberately uses a Sattolo cycle;
// a prefix of such a cycle is not a uniform subset (for example deck[0] can
// never be zero), so it must not select the degree-remainder rows.
void UnbiasedShufflePermutation(
    wirehair::PCGRandom& prng,
    uint16_t* deck,
    uint32_t count)
{
    for (uint32_t i = 0; i < count; ++i) {
        deck[i] = (uint16_t)i;
    }
    for (uint32_t remaining = count; remaining > 1u; --remaining)
    {
        const uint32_t pick = UniformBelow(prng, remaining);
        std::swap(deck[pick], deck[remaining - 1u]);
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

/// Per-source-column staircase degree under a resolved mixture.
///
/// `high_column` is empty in the exactly-regular case (HighColumns == 0), so
/// the certified construction pays no per-column lookup at all.
uint32_t StaircaseColumnHits(
    const StaircaseDegreeMixture& mixture,
    const std::vector<uint8_t>& high_column,
    uint32_t column)
{
    return mixture.LowHits +
        (high_column.empty() ? 0u : (uint32_t)high_column[column]);
}

bool AppendDegreeBalancedStaircaseEdges(
    const PrecodeParams& params,
    const StaircaseDegreeMixture& mixture,
    const std::vector<uint8_t>& high_column,
    std::vector<std::vector<uint32_t>>& rows)
{
    const uint32_t K = params.BlockCount;
    const uint32_t S = params.Staircase;
    const uint32_t edge_count = mixture.EdgeCount;
    const uint32_t low_degree = edge_count / S;
    const uint32_t extra_rows = edge_count % S;
    uint32_t high_degree = low_degree + (extra_rows != 0u);

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    // EXPERIMENT: shaped row-degree target.
    //
    // The column side is an exact two-valued mixture (ValidatePrecodeSystem
    // enforces it, counts included), so the edge budget is fixed once
    // StaircaseDegreeScale is chosen and only its DISTRIBUTION across rows is
    // free.  The peeling literature says that shape matters: a few
    // low-degree rows start the ripple, a few high-degree rows finish coverage.
    // The bucket greedy below realizes ANY feasible target sequence, not just the
    // balanced one, so shaping is a matter of choosing `remaining` differently.
    const std::vector<double>* shape = ActiveStaircaseDegreeDistribution();
    // A PINNED sequence names the object the shape only samples, so it takes
    // precedence over one and the shape path below is skipped entirely.
    const std::vector<uint32_t>* pinned = ActiveStaircaseRowDegrees();
    if (pinned) { shape = nullptr; }
#endif

    // This stream is deliberately independent of the certified construction
    // stream.  BuildPrecodeSystem still consumes the certified staircase draws
    // before calling here, so dense rows remain bit-identical between arms.
    //
    // The key carries K and S because a construction of a different SIZE must
    // draw a different stream; that is stream selection, and it names no
    // parameter.  The third word used to be params.SourceHits -- that is,
    // CertifiedSourceHits(K), the K-THRESHOLDED TABLE `K >= 10000 ? 3 : 2`.
    // It set no parameter either, but it was still a built-in K table feeding
    // the construction, and it demonstrably changed the produced system (7 of
    // 8 SourceHits values gave a different system at every S and scale tested
    // on this path), so it is replaced by the frozen constant below.  The
    // constant is the value that table returned for every K < 10000, which
    // leaves every existing balanced/shaped fixture and the whole K = 128
    // campaign bit-identical; only balanced/shaped experiments at K >= 10000
    // re-roll onto a different, equally valid stream.  Nothing here may be
    // computed from K, from a table, or from another parameter again.
    static const uint64_t kMatchingKeySalt = 2u;
    wirehair::PCGRandom matching_prng;
    matching_prng.Seed(
        params.Seed ^ UINT64_C(0x8b8b8b8bd3c4a56f),
        ((uint64_t)K << 32) ^ ((uint64_t)S << 16) ^
            kMatchingKeySalt ^ UINT64_C(0x6465677265656d61));

    // Randomly choose which rows receive the remainder socket, then match each
    // source's sockets to distinct rows with the greatest residual capacity.
    // Holding a source's selected rows out of the buckets until all its sockets
    // are assigned prevents duplicate edges.  Highest-capacity selection keeps
    // the residual sequence graphical through the final source.
    std::vector<uint16_t> row_deck(S);
    UnbiasedShufflePermutation(matching_prng, row_deck.data(), S);
    std::vector<uint32_t> remaining(S, low_degree);
    for (uint32_t i = 0; i < extra_rows; ++i) {
        ++remaining[row_deck[i]];
    }

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    if (pinned)
    {
        // Consume exactly the draws the sampler would have made, so a pinned
        // sequence lands on the SAME matching stream the shape path reaches
        // when it happens to sample this sequence.  Without this a pinned run
        // and the shaped run it reproduces would differ in the bucket picks
        // and could not be compared cell for cell.
        for (uint32_t row = 0; row < S; ++row) { (void)matching_prng.Next(); }
        if (pinned->size() != (size_t)S) { return false; }
        uint64_t assigned = 0;
        for (uint32_t row = 0; row < S; ++row)
        {
            const uint32_t d = (*pinned)[row];
            if (d > K) { return false; }
            remaining[row] = d;
            assigned += d;
        }
        // The budget is fixed by StaircaseDegreeScale, so a sequence that does
        // not sum to it is a different construction than the one requested and
        // is refused rather than repaired -- a repair would silently evaluate a
        // NEIGHBOUR of the sequence under test, which is the whole failure mode
        // this hook exists to remove.
        if (assigned != (uint64_t)edge_count) { return false; }
        high_degree = 0u;
        for (uint32_t row = 0; row < S; ++row) {
            high_degree = std::max(high_degree, remaining[row]);
        }
    }
    else if (shape)
    {
        // Draw each row's target degree from the shape, capped at K (a row
        // cannot touch a source column twice), then repair the total back to
        // the fixed budget: hand out or reclaim single sockets in shuffled row
        // order so the correction itself does not bias any particular row.
        double total = 0.0;
        for (double w : *shape) { total += w > 0.0 ? w : 0.0; }
        if (total <= 0.0) { return false; }
        uint64_t assigned = 0;
        for (uint32_t row = 0; row < S; ++row)
        {
            double u = (double)(matching_prng.Next() >> 8) / 16777216.0 * total;
            uint32_t degree = 1u;
            for (uint32_t d = 0; d < (uint32_t)shape->size(); ++d) {
                const double w = (*shape)[d] > 0.0 ? (*shape)[d] : 0.0;
                degree = d + 1u;
                if (u < w) { break; }
                u -= w;
            }
            remaining[row] = degree > K ? K : degree;
            assigned += remaining[row];
        }
        // The edge budget is already fixed by StaircaseDegreeScale (or, with
        // no scale supplied, by the legacy rule), so only the SHAPE across
        // rows is free here.  Rescale the sampled degrees multiplicatively
        // onto that budget -- an additive +1/-1 repair would flatten the very
        // shape being searched, which makes most of the parameter space
        // degenerate.
        if (assigned > 0u && assigned != edge_count)
        {
            const double scale = (double)edge_count / (double)assigned;
            assigned = 0;
            for (uint32_t row = 0; row < S; ++row) {
                double scaled = (double)remaining[row] * scale + 0.5;
                uint32_t d = (uint32_t)scaled;
                if (d > K) { d = K; }
                remaining[row] = d;
                assigned += d;
            }
        }
        // Settle rounding and saturation in one deterministic pass.  A
        // single-edge loop used to stop after 4096*S attempts, which made
        // feasible shapes fail merely because a capped sample left a large
        // correction.  Each row takes its even share of the remaining delta,
        // plus any amount that MUST be taken because the later rows lack the
        // capacity.  The suffix-capacity term proves completion in O(S).
        if (assigned < edge_count)
        {
            uint64_t needed = (uint64_t)edge_count - assigned;
            uint64_t suffix_capacity = 0u;
            for (uint32_t row = 0; row < S; ++row) {
                suffix_capacity += K - remaining[row];
            }
            if (needed > suffix_capacity) {
                return false;
            }
            for (uint32_t i = 0; i < S && needed != 0u; ++i)
            {
                const uint32_t row = row_deck[i];
                const uint64_t capacity = K - remaining[row];
                suffix_capacity -= capacity;
                const uint64_t rows_left = S - i;
                const uint64_t even_share =
                    (needed + rows_left - 1u) / rows_left;
                const uint64_t required =
                    needed > suffix_capacity ? needed - suffix_capacity : 0u;
                const uint64_t take = std::min(
                    capacity, std::max(even_share, required));
                remaining[row] += (uint32_t)take;
                needed -= take;
            }
            assigned = (uint64_t)edge_count - needed;
        }
        else if (assigned > edge_count)
        {
            uint64_t excess = assigned - (uint64_t)edge_count;
            uint64_t suffix_capacity = assigned;
            if (excess > suffix_capacity) {
                return false;
            }
            for (uint32_t i = 0; i < S && excess != 0u; ++i)
            {
                const uint32_t row = row_deck[i];
                const uint64_t capacity = remaining[row];
                suffix_capacity -= capacity;
                const uint64_t rows_left = S - i;
                const uint64_t even_share =
                    (excess + rows_left - 1u) / rows_left;
                const uint64_t required =
                    excess > suffix_capacity ? excess - suffix_capacity : 0u;
                const uint64_t take = std::min(
                    capacity, std::max(even_share, required));
                remaining[row] -= (uint32_t)take;
                excess -= take;
            }
            assigned = (uint64_t)edge_count + excess;
        }
        if (assigned != edge_count) { return false; }
        high_degree = 0u;
        for (uint32_t row = 0; row < S; ++row) {
            high_degree = std::max(high_degree, remaining[row]);
        }
    }
#endif

    std::vector<std::vector<uint16_t>> buckets(high_degree + 1u);
    for (uint32_t row = 0; row < S; ++row) {
        buckets[remaining[row]].push_back((uint16_t)row);
    }

    uint32_t active_degree = high_degree;
    // Per-column degrees are bounded by min(S, kStaircaseColumnHitsMax), not
    // by the [1,8] SourceHits domain, once the budget is scaled.  The inline
    // arrays still cover every unscaled configuration.
    uint16_t selected_rows_inline[8];
    uint32_t selected_degrees_inline[8];
    std::vector<uint16_t> selected_rows_heap;
    std::vector<uint32_t> selected_degrees_heap;
    uint16_t* selected_rows = selected_rows_inline;
    uint32_t* selected_degrees = selected_degrees_inline;
    if (mixture.MaxHits > 8u)
    {
        selected_rows_heap.resize(mixture.MaxHits);
        selected_degrees_heap.resize(mixture.MaxHits);
        selected_rows = selected_rows_heap.data();
        selected_degrees = selected_degrees_heap.data();
    }
    for (uint32_t source = 0; source < K; ++source)
    {
        // Bipartite Havel-Hakimi: connecting ANY column to the columns' worth
        // of rows with the largest residual capacity preserves realizability,
        // so a mixture of column degrees needs no reordering of the sources --
        // which matters because appending in ascending source order is what
        // keeps every row's column list sorted.
        const uint32_t hits =
            StaircaseColumnHits(mixture, high_column, source);
        for (uint32_t hit = 0; hit < hits; ++hit)
        {
            while (active_degree > 0u && buckets[active_degree].empty()) {
                --active_degree;
            }
            if (active_degree == 0u) {
                return false;
            }
            std::vector<uint16_t>& bucket = buckets[active_degree];
            const uint32_t pick = UniformBelow(
                matching_prng, (uint32_t)bucket.size());
            const uint16_t row = bucket[pick];
            bucket[pick] = bucket.back();
            bucket.pop_back();
            selected_rows[hit] = row;
            selected_degrees[hit] = active_degree;
        }
        for (uint32_t hit = 0; hit < hits; ++hit)
        {
            const uint32_t next_degree = selected_degrees[hit] - 1u;
            buckets[next_degree].push_back(selected_rows[hit]);
            active_degree = std::max(active_degree, next_degree);
            rows[selected_rows[hit]].push_back(source);
        }
    }

    if (buckets[0].size() != S) {
        return false;
    }
    for (size_t degree = 1; degree < buckets.size(); ++degree) {
        if (!buckets[degree].empty()) {
            return false;
        }
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

bool IsSegmentedDenseAnchor(
    DenseAnchorLayout layout,
    uint32_t row_index)
{
    switch (layout)
    {
    case DenseAnchorLayout::Disabled:
        return false;
    case DenseAnchorLayout::Three048:
        return row_index == 4u || row_index == 8u;
    case DenseAnchorLayout::Three059:
        return row_index == 5u || row_index == 9u;
    case DenseAnchorLayout::Four0369:
        return row_index == 3u || row_index == 6u || row_index == 9u;
    }
    return false;
}

bool IsValidDenseAnchorLayout(DenseAnchorLayout layout)
{
    return layout == DenseAnchorLayout::Disabled ||
        layout == DenseAnchorLayout::Three048 ||
        layout == DenseAnchorLayout::Three059 ||
        layout == DenseAnchorLayout::Four0369;
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
        !IsValidDenseAnchorLayout(params.SegmentedDenseAnchors) ||
        binary_span > UINT16_MAX || total_span > UINT16_MAX)
    {
        return false;
    }

    // The scaled edge budget has to resolve to a realizable column mixture.
    // Rejecting here rather than asserting keeps an out-of-domain scale the
    // same kind of outcome as any other bad parameter: a refused build.
    StaircaseDegreeMixture mixture;
    if (!MakeStaircaseDegreeMixture(params, mixture)) {
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
    if ((params.DenseTwoAnchor &&
         (params.DenseRows != 12u || params.DenseIdentityCorner ||
          params.DenseTwoAnchorPhase > 2u ||
          params.SegmentedDenseAnchors != DenseAnchorLayout::Disabled)) ||
        (!params.DenseTwoAnchor && params.DenseTwoAnchorPhase != 0u))
    {
        return false;
    }
    if (params.SegmentedDenseAnchors != DenseAnchorLayout::Disabled &&
        (params.DenseRows != 12u || params.DenseIdentityCorner))
    {
        return false;
    }
    return !params.DenseIdentityCorner ||
        known_span >= 2u * (uint64_t)(params.DenseRows >> 1);
}

} // namespace

bool MakeStaircaseDegreeMixture(
    const PrecodeParams& params,
    StaircaseDegreeMixture& out)
{
    out = StaircaseDegreeMixture();

    const uint32_t K = params.BlockCount;
    const uint32_t S = params.Staircase;
    // Deliberately re-checked here instead of calling ValidatePrecodeParams:
    // that function calls this one, and these are the only inputs the mixture
    // actually reads.  SourceHits keeps its domain check even when a scale is
    // supplied and it feeds nothing: it is still a recorded, fingerprinted
    // field of the certified profile, so an out-of-domain value stays a
    // refused parameter set rather than a silently ignored one.
    if (K < 2u || K > 64000u || S == 0u || S > UINT16_MAX ||
        params.SourceHits == 0u || params.SourceHits > 8u)
    {
        return false;
    }
    const double scale = params.StaircaseDegreeScale;
    const bool scale_supplied = scale != kStaircaseDegreeScaleUnset;
    // Only the declared sentinel means "not supplied."  Treating every
    // negative value (including -infinity) as the sentinel would let a broken
    // direct caller silently select the legacy budget.  Written as positive
    // domain tests so NaN and both infinities are rejected too.
    if (scale_supplied &&
        (!(scale >= kStaircaseDegreeScaleMin) ||
         !(scale <= kStaircaseDegreeScaleMax)))
    {
        return false;
    }

    uint64_t edge_count;
    if (scale_supplied)
    {
        // total = mean x count.  The scale NAMES the mean staircase row
        // degree and there are S rows, so this product is a definitional
        // identity, not a model: no K, no SourceHits, no table, no rule.
        const long long rounded = std::llround(scale * (double)S);
        edge_count = rounded <= 0 ? 0u : (uint64_t)rounded;
    }
    else
    {
        // Legacy rule, and the ONLY place it survives: when no scale is
        // supplied the certified construction's K * min(N1,S) budget applies
        // unchanged, so production and every un-scaled test is bit-identical.
        edge_count = (uint64_t)K * std::min(params.SourceHits, S);
    }
    // A budget above K*S cannot be realized at all -- a row may not touch a
    // column twice, so K*S is every edge the bipartite graph can hold.  That
    // ceiling CLAMPS rather than rejects, so a search walking off the end of
    // its box lands on the saturated construction instead of on a hole.
    //
    // There is deliberately NO floor.  A budget below S leaves some rows with
    // no source column and a budget below K leaves lo = 0, i.e. some columns
    // with no parity at all; both are legal here.  Flooring either one up
    // would put a K-derived (or S-derived) rule back into the construction,
    // and the search's own failure constraint is what rejects a budget that
    // is too small.
    const uint64_t max_edges = (uint64_t)K * S;
    if (edge_count > max_edges) {
        edge_count = max_edges;
    }
    if (edge_count > UINT32_MAX) {
        return false;
    }

    const uint32_t low_hits = (uint32_t)(edge_count / K);
    const uint32_t high_columns =
        (uint32_t)(edge_count - (uint64_t)low_hits * K);
    const uint32_t max_hits = low_hits + (high_columns != 0u ? 1u : 0u);
    // A column cannot take the same staircase row twice, and the validator
    // tallies per-column hits in a byte.
    if (max_hits > S || max_hits > kStaircaseColumnHitsMax) {
        return false;
    }

    out.EdgeCount = (uint32_t)edge_count;
    out.LowHits = low_hits;
    out.HighHits = low_hits + 1u;
    out.HighColumns = high_columns;
    out.LowColumns = K - high_columns;
    out.MaxHits = max_hits;
    return true;
}

uint32_t SmallBandStaircaseCount(uint32_t block_count)
{
    if (block_count < 2u || block_count > 64000u) {
        return 0u;
    }
    const uint32_t inherited = wirehair::GetDenseCount(block_count);
    if (block_count > kSmallBandStaircaseMaxBlockCount) {
        return inherited;
    }
    // wirehair::GetDenseCount() sizes Wirehair 1's *dense GE row* count.  V2
    // reuses that table for a structurally different quantity -- the LDPC
    // staircase parity count S -- and in the K <= 100 band the inherited
    // value is far larger than the staircase needs: S is 26 at K=100 and 13
    // at K=32.  Each surplus parity adds one binary row *and* one binary
    // column to every solve phase.  A measured sweep of S against untuned
    // fail@+0 (3000-4000 paired constructions per point) is flat from the
    // inherited value down to roughly 1.25*sqrt(K), and only starts to
    // degrade below it -- at K=100, S=26 fails 0.93% and S=12 fails 0.87%,
    // while S=8 fails 1.13% and S=4 fails 1.77%.  Sizing S at the knee keeps
    // reliability and drops ~11% of encode time.
    uint32_t root = 1u;
    while ((root + 1u) * (root + 1u) <= block_count) {
        ++root;
    }
    // floor(1.25 * (root + K/root) / 2): one Newton refinement of the integer
    // square root, scaled by 5/4, all in exact integer arithmetic so the rule
    // is reproducible on every platform.  Only reached for K <= 100, so the
    // 5*(K + root^2) product cannot overflow.  Gives S = 3 at K=8, 5 at
    // K=16, 7 at K=32, 10 at K=64 and 12 at K=100.
    const uint32_t scaled = (5u * (block_count + root * root)) /
        (8u * root);
    uint32_t staircase = scaled < 2u ? 2u : scaled;
    if (staircase > inherited) {
        staircase = inherited;
    }
    return staircase;
}

PrecodeParams MakeCertifiedParams(uint32_t block_count, uint64_t seed)
{
    PrecodeParams params;
    params.BlockCount = block_count;
    params.Staircase = SmallBandStaircaseCount(block_count);
    params.DenseRows = 12u;
    params.HeavyRows = 12u;
    params.SourceHits = CertifiedSourceHits(block_count);
    // The certified path supplies NO target mean row degree, so the legacy
    // K * min(N1,S) budget applies and this function is unchanged.
    params.StaircaseDegreeScale = kStaircaseDegreeScaleUnset;
    params.DegreeBalancedStaircase = false;
    params.DenseIdentityCorner = false;
    params.DenseTwoAnchor = false;
    params.DenseTwoAnchorPhase = 0u;
    params.SegmentedDenseAnchors = DenseAnchorLayout::Disabled;
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
    // N1 no longer appears anywhere below: the per-column parity count comes
    // from the resolved degree mixture, which is a function of the target mean
    // row degree and S alone whenever a scale is supplied.
    // The target mean row degree is the one parameter a test hook may
    // override per thread, so the EFFECTIVE parameter set is resolved before
    // anything else and is what gets validated, built and recorded in
    // out.Params -- otherwise ValidatePrecodeSystem(out) would check the
    // system against a scale it was not built at.
    PrecodeParams effective = params;
    effective.StaircaseDegreeScale =
        ActiveStaircaseDegreeScale(params.StaircaseDegreeScale);
    // Hooks that fail to parse are not silently replaced by the stock graph,
    // and a shaped or pinned row sequence cannot truthfully be recorded as
    // DegreeBalancedStaircase.  Reject both conflicts before touching `out`.
    if (!StaircaseHookConfigurationValid()) {
        return false;
    }
    const bool shaped_staircase = ShapedStaircaseActive();
    if (effective.DegreeBalancedStaircase && shaped_staircase) {
        return false;
    }
    // Reject the complete parameter domain before modifying `out` or
    // allocating row/deck storage.
    if (!ValidatePrecodeParams(effective)) {
        return false;
    }
    StaircaseDegreeMixture mixture;
    if (!MakeStaircaseDegreeMixture(effective, mixture) ||
        !ActiveStaircaseHookValuesValid(effective, mixture))
    {
        return false;
    }
    const uint64_t span_wide = (uint64_t)K + S + D2;
    const uint32_t span = (uint32_t)span_wide;

    out.Params = effective;
    out.StaircaseRows.assign(S, std::vector<uint32_t>());
    out.DenseRowColumns.assign(D2, std::vector<uint32_t>());

    wirehair::PCGRandom prng;
    prng.Seed(params.Seed, K);

    // --- Staircase rows ---
    // Reserve for the heaviest plausible load, not the mean.  Source hits are
    // EdgeCount balls thrown into S bins, so a row's source degree is
    // concentrated around the mean with a spread of order sqrt(mean); a
    // mean-sized reserve leaves a large fraction of the rows to reallocate
    // mid-build, which costs a second allocation and a copy on each of them.
    {
        const uint32_t mean = mixture.EdgeCount / S;
        uint32_t slack = 4u;
        while (slack * slack < 16u * mean) {
            ++slack;
        }
        const uint32_t expected = mean + slack + 3u;
        for (uint32_t j = 0; j < S; ++j) {
            out.StaircaseRows[j].reserve(expected);
        }
    }

    // Which source columns carry the extra parity, when the scaled budget is
    // not a whole multiple of K.  Drawn from the construction stream so the
    // choice is reproducible, and drawn ONLY when the mixture is non-trivial
    // so that the exactly-regular case (every unscaled construction) consumes
    // the identical PRNG sequence it always has.
    std::vector<uint8_t> high_column;
    if (mixture.HighColumns != 0u)
    {
        // K <= 64000 fits the uint16 deck the shuffle helpers use.
        std::vector<uint16_t> column_deck(K);
        UnbiasedShufflePermutation(prng, column_deck.data(), K);
        high_column.assign(K, uint8_t{0});
        for (uint32_t i = 0; i < mixture.HighColumns; ++i) {
            high_column[column_deck[i]] = 1u;
        }
    }

    // Consume the certified independent-placement stream in both modes so the
    // later dense construction remains bit-identical.  In the experimental
    // balanced mode these placements are discarded and replaced below.
    {
        // The [1,8] SourceHits domain bounds the unscaled per-column degree,
        // but a scaled budget can push it up to min(S, 255).
        uint32_t picks_inline[8];
        std::vector<uint32_t> picks_heap;
        uint32_t* picks = picks_inline;
        if (mixture.MaxHits > 8u)
        {
            picks_heap.resize(mixture.MaxHits);
            picks = picks_heap.data();
        }
        for (uint32_t c = 0; c < K; ++c)
        {
            const uint32_t hits = StaircaseColumnHits(mixture, high_column, c);
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
                if (!effective.DegreeBalancedStaircase &&
                    !shaped_staircase) {
                    out.StaircaseRows[p].push_back(c);
                }
            }
        }
    }
    if ((effective.DegreeBalancedStaircase || shaped_staircase) &&
        !AppendDegreeBalancedStaircaseEdges(
            effective, mixture, high_column, out.StaircaseRows))
    {
        return false;
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
        // The deck is overwritten in full by every shuffle and the bitmap is
        // cleared before use, so neither needs a zero-initialised heap block.
        // Spans up to the inline bound -- which covers the whole small-K band
        // this codec is tuned for -- avoid the allocation entirely.
        static const uint32_t kInlineDeckSpan = 512u;
        uint16_t deck_inline[kInlineDeckSpan];
        uint8_t bitmap_inline[kInlineDeckSpan];
        std::vector<uint16_t> deck_heap;
        std::vector<uint8_t> bitmap_heap;
        uint16_t* deck = deck_inline;
        uint8_t* bitmap = bitmap_inline;
        if (deck_span > kInlineDeckSpan)
        {
            deck_heap.resize(deck_span);
            bitmap_heap.resize(deck_span);
            deck = deck_heap.data();
            bitmap = bitmap_heap.data();
        }
        std::memset(bitmap, 0, deck_span);

        UnbiasedShuffleDeck(prng, deck, deck_span);
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

        if (params.SegmentedDenseAnchors != DenseAnchorLayout::Disabled)
        {
            // Match the certified first segment's RNG convention: row 0 is
            // the balanced half of the initial deck, while its cheap flips
            // use one fresh shuffle.  Every later segment starts with the
            // balanced half of a fresh deck and continues through distinct
            // set/clear swap pairs from that same deck.
            UnbiasedShuffleDeck(prng, deck, deck_span);
            uint32_t flip_index = 0u;
            while (row_i < D2)
            {
                if (IsSegmentedDenseAnchor(
                        params.SegmentedDenseAnchors, row_i))
                {
                    UnbiasedShuffleDeck(prng, deck, deck_span);
                    std::memset(bitmap, 0, deck_span);
                    for (uint32_t i = 0; i < set_count; ++i) {
                        bitmap[deck[i]] = 1u;
                    }
                    flip_index = 0u;
                    emit_row();
                    continue;
                }
                bitmap[deck[flip_index]] ^= 1u;
                bitmap[deck[set_count + flip_index]] ^= 1u;
                ++flip_index;
                emit_row();
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
                // A half's shuffled deck is observable only through the rows
                // that half emits and through the PRNG state a later half
                // consumes.  When neither exists -- the D2 == 1 configuration
                // makes both halves empty, and the generator is not consulted
                // after this loop -- the shuffle is dead work, and skipping it
                // leaves the produced system bit-identical.
                bool deck_observable =
                    halves[half] != 0u ||
                    (params.DenseTwoAnchor && half == 1u);
                for (uint32_t rest = half + 1u; rest < 2u; ++rest) {
                    if (halves[rest] != 0u ||
                        (params.DenseTwoAnchor && rest == 1u))
                    {
                        deck_observable = true;
                    }
                }
                if (!deck_observable) {
                    continue;
                }
                UnbiasedShuffleDeck(prng, deck, deck_span);
                uint32_t flip_count = halves[half];
                uint32_t flip_offset = 0u;
                if (params.DenseTwoAnchor && half == 1u)
                {
                    // Reuse the certified second-half shuffle, but turn its
                    // balanced set-half into a new dense equation.  One of
                    // the five baseline flips becomes the anchor emission,
                    // keeping D2, RNG consumption, and rows 0..6 unchanged.
                    std::memset(bitmap, 0, deck_span);
                    for (uint32_t i = 0; i < set_count; ++i) {
                        bitmap[deck[i]] = 1u;
                    }
                    // Experimental q1/q2 phases rotate the five-row anchor
                    // window forward within this same deck.
                    flip_offset = params.DenseTwoAnchorPhase;
                    for (uint32_t ii = 0; ii < flip_offset; ++ii)
                    {
                        bitmap[deck[ii]] ^= 1u;
                        bitmap[deck[set_count + ii]] ^= 1u;
                    }
                    emit_row();
                    --flip_count;
                }
                for (uint32_t ii = 0; ii < flip_count; ++ii)
                {
                    // Deck entries at distinct positions are distinct
                    // columns, so each flip pair changes exactly two columns.
                    const uint32_t deck_index = flip_offset + ii;
                    bitmap[deck[deck_index]] ^= 1u;
                    bitmap[deck[set_count + deck_index]] ^= 1u;
                    emit_row();
                }
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

    StaircaseDegreeMixture mixture;
    if (!ValidatePrecodeParams(params) ||
        !MakeStaircaseDegreeMixture(params, mixture) ||
        system.StaircaseRows.size() != S ||
        system.DenseRowColumns.size() != D2)
    {
        return false;
    }

    const uint32_t staircase_end = K + S;
    // Validation runs on every generated system, so keep the per-source hit
    // tally off the heap for the block counts this codec is tuned for.
    static const uint32_t kInlineSourceHits = 1024u;
    uint8_t source_hits_inline[kInlineSourceHits];
    std::vector<uint8_t> source_hits_heap;
    uint8_t* source_hits = source_hits_inline;
    if (K > kInlineSourceHits) {
        source_hits_heap.assign(K, uint8_t{0});
        source_hits = source_hits_heap.data();
    }
    else {
        std::memset(source_hits, 0, K);
    }
    uint32_t balanced_high_rows = 0u;
    const uint32_t balanced_edges = mixture.EdgeCount;
    const uint32_t balanced_low = balanced_edges / S;
    const uint32_t balanced_extra = balanced_edges % S;
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
        uint32_t source_degree = 0u;
        for (uint32_t column : row)
        {
            if (column < K)
            {
                if (source_hits[column] == UINT8_MAX) {
                    return false;
                }
                ++source_hits[column];
                ++source_degree;
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
        if (params.DegreeBalancedStaircase)
        {
            if (source_degree == balanced_low + 1u && balanced_extra != 0u) {
                ++balanced_high_rows;
            }
            else if (source_degree != balanced_low) {
                return false;
            }
        }
    }
    if (params.DegreeBalancedStaircase &&
        balanced_high_rows != balanced_extra)
    {
        return false;
    }
    // Exact two-valued column mixture.  This is the generalization of the old
    // "every column has exactly min(N1,S) parities" invariant and it is just
    // as strict: a column outside {lo, hi} is rejected, AND the counts on each
    // of the two allowed values must match the mixture exactly.  A range check
    // would accept any split between lo and hi, which is precisely the edge
    // budget this parameter exists to control.
    uint32_t low_columns = 0u;
    uint32_t high_columns = 0u;
    for (uint32_t column = 0; column < K; ++column)
    {
        const uint32_t degree = source_hits[column];
        if (degree == mixture.LowHits) {
            ++low_columns;
        }
        else if (mixture.HighColumns != 0u && degree == mixture.HighHits) {
            ++high_columns;
        }
        else {
            return false;
        }
    }
    if (low_columns != mixture.LowColumns ||
        high_columns != mixture.HighColumns)
    {
        return false;
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
        const uint32_t second_anchor = 1u + (D2 >> 1);
        const bool is_anchor = row_index == 0u ||
            (params.DenseTwoAnchor && row_index == second_anchor) ||
            IsSegmentedDenseAnchor(
                params.SegmentedDenseAnchors, row_index);
        if (is_anchor)
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

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
/// When set, a trimmed completion band re-derives its Cauchy X coordinates
/// from its own row count instead of inheriting the frozen H12 ones.  Read
/// once from WIREHAIR_V2_BAND_TRACKING_X so a sweep can toggle it per process
/// without threading a new flag through every bench mode.
static std::atomic<int> MixedBandTrackingXOverride{-1};

static bool MixedBandTrackingXEnabled()
{
    const int override_value =
        MixedBandTrackingXOverride.load(std::memory_order_relaxed);
    if (override_value >= 0) {
        return override_value != 0;
    }
    static const bool enabled = [] {
        const wirehair::EnvironmentValue value(
            "WIREHAIR_V2_BAND_TRACKING_X");
        return value.IsSet() && value.Get()[0] == '1';
    }();
    return enabled;
}
#define MixedBandTrackingXForTesting MixedBandTrackingXEnabled()

#endif

const MixedCoefficientRows* GetMixedCoefficientRows()
{
    if (!InitializeGF16()) {
        return nullptr;
    }
    const auto build_rows = [](MixedCoefficientGeometry geometry,
                               uint32_t band_h) {
        MixedCoefficientRows result = {};
        // The first ten rows, optional Y=11/Y=46 test rows, and shared-X
        // extension rows deliberately keep
        // the frozen H12 X coordinates [12, 256).  Later test geometries
        // append Y values; moving X with H would rewrite every existing
        // coefficient.
        //
        // EXPERIMENT (small-K speed): band_h lets a trimmed band re-derive its
        // own X coordinates instead of inheriting the H12 ones.  A Cauchy
        // matrix needs X and Y drawn from disjoint sets; trimming Y to
        // [0, H') while X stays at 12 + residue keeps them disjoint but no
        // longer matches the geometry the coefficients were chosen for, which
        // is the suspected source of the per-K rank spikes below H=12.
        const uint32_t H = band_h;
        for (uint32_t residue = 0;
             residue < kMixedCoefficientPeriod;
             ++residue)
        {
            for (uint32_t row = 0; row < kMixedGF256RowsMax; ++row) {
                // Y=10 develops a singular H15 independently scheduled
                // corner at K=23092.  Y=11 is the other unused subfield
                // coordinate below X=12 and is exhaustively nonsingular for
                // every K=2..64000 under the active keyed P31/P32
                // constructions.
                // With Y=11 retained, Y=46 also keeps the H16 corner
                // nonsingular for every K=2..64000 under the corresponding
                // independently scheduled P31/P32 constructions.
                uint32_t cauchy_y = row;
                if (row == kMixedGF256Rows) cauchy_y = 11u;
                else if (row == kMixedGF256Rows + 1u) cauchy_y = 46u;
                const uint32_t cauchy_x = H + residue;
                result.Subfield[row][residue] = cauchy_y < H ?
                    HeavyCoefficient(cauchy_y, residue, H) :
                    gf256_inv((uint8_t)(cauchy_x ^ cauchy_y));
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
            build_rows(MixedCoefficientGeometry::SharedCauchyX,
                       kMixedGF256Rows + kMixedGF16Rows);
        const uint32_t grouped_rows = ActiveMixedGroupedGF256Rows();
        const uint32_t row_mask = ActiveMixedGroupedGF256RowMask();
        const uint32_t suffix_mask = grouped_rows == 0u ? 0u :
            ((UINT32_C(1) << grouped_rows) - 1u) <<
                (kMixedGF256Rows - grouped_rows);
        if (row_mask != suffix_mask)
        {
            // The grouped kernels deliberately retain one compact A prefix
            // and C suffix.  Reorder only the equation rows so the requested
            // Cauchy Y coordinates occupy that suffix.  A row permutation
            // changes neither the solution set nor the final all-A corner,
            // and avoids adding membership branches to any payload loop.
            static thread_local MixedCoefficientRows permuted_rows = {};
            static thread_local uint32_t cached_row_mask = UINT32_MAX;
            if (cached_row_mask != row_mask)
            {
                permuted_rows = shared_rows;
                uint32_t destination = 0u;
                for (uint32_t source = 0u;
                     source < kMixedGF256Rows; ++source)
                {
                    if ((row_mask & (UINT32_C(1) << source)) != 0u) {
                        continue;
                    }
                    std::copy(
                        shared_rows.Subfield[source],
                        shared_rows.Subfield[source] +
                            kMixedCoefficientPeriod,
                        permuted_rows.Subfield[destination++]);
                }
                for (uint32_t source = 0u;
                     source < kMixedGF256Rows; ++source)
                {
                    if ((row_mask & (UINT32_C(1) << source)) == 0u) {
                        continue;
                    }
                    std::copy(
                        shared_rows.Subfield[source],
                        shared_rows.Subfield[source] +
                            kMixedCoefficientPeriod,
                        permuted_rows.Subfield[destination++]);
                }
                CAT_DEBUG_ASSERT(destination == kMixedGF256Rows);
                cached_row_mask = row_mask;
            }
            return &permuted_rows;
        }
        return &shared_rows;
    }
    if (MixedBandTrackingXForTesting)
    {
        // Rebuilt whenever the active band size changes: the table is a pure
        // function of (geometry, band H), so caching on H is exact.
        const uint32_t band_h =
            ActiveMixedGF256Rows() + ActiveMixedGF16Rows();
        static thread_local MixedCoefficientRows tracked_rows = {};
        static thread_local uint32_t cached_band_h = UINT32_MAX;
        if (cached_band_h != band_h) {
            tracked_rows =
                build_rows(MixedCoefficientGeometry::FrozenPowerX, band_h);
            cached_band_h = band_h;
        }
        return &tracked_rows;
    }
#endif
    static const MixedCoefficientRows frozen_rows =
        build_rows(MixedCoefficientGeometry::FrozenPowerX,
                   kMixedGF256Rows + kMixedGF16Rows);
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
        if (ActiveMixedGroupedGF256Rows() != 0u)
        {
            // `rows` may be the thread-local row-permuted view.  Keep its
            // packed lanes in a matching thread-local cache; static shared
            // caches below remain canonical for every ungrouped caller.
            static thread_local MixedPackedCoefficients grouped_packed = {};
            static thread_local uint32_t cached_row_mask = UINT32_MAX;
            const uint32_t row_mask = ActiveMixedGroupedGF256RowMask();
            if (cached_row_mask != row_mask) {
                grouped_packed = pack_rows(rows, kMixedGF256Rows);
                cached_row_mask = row_mask;
            }
            return &grouped_packed;
        }
        const uint32_t subfield_rows = ActiveMixedGF256Rows();
        // Retain the original shared immutable caches for the historically
        // benchmarked 8..12-row arms.  Besides avoiding first-use TLS packing
        // in their timings, the production-sized H8=10 table is deliberately
        // published at one process-wide address and has a concurrency oracle.
        static const MixedPackedCoefficients shared_packed_trim_two =
            pack_rows(rows, kMixedGF256Rows - 2u);
        static const MixedPackedCoefficients shared_packed_trim_one =
            pack_rows(rows, kMixedGF256Rows - 1u);
        static const MixedPackedCoefficients shared_packed_base =
            pack_rows(rows, kMixedGF256Rows);
        static const MixedPackedCoefficients shared_packed_one_extra =
            pack_rows(rows, kMixedGF256Rows + 1u);
        static const MixedPackedCoefficients shared_packed_two_extra =
            pack_rows(rows, kMixedGF256RowsMax);
        if (subfield_rows == kMixedGF256Rows - 2u)
            return &shared_packed_trim_two;
        if (subfield_rows == kMixedGF256Rows - 1u)
            return &shared_packed_trim_one;
        if (subfield_rows == kMixedGF256Rows) return &shared_packed_base;
        if (subfield_rows == kMixedGF256Rows + 1u)
            return &shared_packed_one_extra;
        if (subfield_rows == kMixedGF256RowsMax)
            return &shared_packed_two_extra;
        // Test builds expose every leading GF(256) prefix from 1 through 12.
        // The GF(2^16) lanes must begin immediately after that ACTIVE prefix:
        // reusing an H8=9/12 packed table for a smaller prefix silently shifts
        // every extension lane and makes the payload equations disagree with
        // the row-major encoder.  One TLS cache avoids twelve large static
        // tables while keeping repeated solves at one configuration free of
        // repacking work.
        static thread_local MixedPackedCoefficients shared_packed = {};
        static thread_local uint32_t cached_subfield_rows = UINT32_MAX;
        if (cached_subfield_rows != subfield_rows) {
            shared_packed = pack_rows(rows, subfield_rows);
            cached_subfield_rows = subfield_rows;
        }
        return &shared_packed;
    }
    // Frozen geometry supports every leading prefix through H8=10.  With
    // band-tracking enabled the row coefficients also change with H8+H16, so
    // both quantities belong in the cache key even though the row-table
    // pointer itself remains the same thread-local object.
    const uint32_t subfield_rows = ActiveMixedGF256Rows();
    const bool track_band = MixedBandTrackingXForTesting;
    if (!track_band)
    {
        static const MixedPackedCoefficients frozen_packed_trim_two =
            pack_rows(rows, kMixedGF256Rows - 2u);
        static const MixedPackedCoefficients frozen_packed_trim_one =
            pack_rows(rows, kMixedGF256Rows - 1u);
        static const MixedPackedCoefficients frozen_packed_base =
            pack_rows(rows, kMixedGF256Rows);
        if (subfield_rows == kMixedGF256Rows - 2u)
            return &frozen_packed_trim_two;
        if (subfield_rows == kMixedGF256Rows - 1u)
            return &frozen_packed_trim_one;
        if (subfield_rows == kMixedGF256Rows)
            return &frozen_packed_base;
    }
    const uint32_t band_h = track_band ?
        subfield_rows + ActiveMixedGF16Rows() : 0u;
    const uint32_t key = (band_h << 8) | subfield_rows;
    static thread_local MixedPackedCoefficients frozen_test_packed = {};
    static thread_local uint32_t cached_frozen_key = UINT32_MAX;
    if (cached_frozen_key != key) {
        frozen_test_packed = pack_rows(rows, subfield_rows);
        cached_frozen_key = key;
    }
    return &frozen_test_packed;
#endif
    static const MixedPackedCoefficients frozen_packed =
        pack_rows(rows, kMixedGF256Rows);
    return &frozen_packed;
}

bool MixedJointResidueBucketStorageFits(
    uint32_t coefficient_period,
    uint32_t block_bytes,
    uint64_t data_byte_limit)
{
    if (coefficient_period == 0u ||
        coefficient_period > kMixedCoefficientPeriod ||
        block_bytes == 0u || block_bytes > 0x7fffffffu)
    {
        return false;
    }
    const uint64_t plane_bytes =
        (uint64_t)coefficient_period * block_bytes;
    return plane_bytes <= std::numeric_limits<size_t>::max() &&
        plane_bytes <= data_byte_limit / 3u;
}

bool UseAutomaticMixedJointResidueBuckets(
    uint32_t block_count,
    uint32_t block_bytes,
    uint32_t coefficient_period)
{
    // Pinned ABBA measurements cover P=32 only.  Larger periods multiply the
    // marginal work and scratch, so they remain explicitly opt-in in test
    // builds until separately benchmarked.
    return coefficient_period == 32u &&
        ((block_bytes >= 4096u && block_count >= 3200u) ||
         (block_bytes >= 1280u && block_count >= 10000u));
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
static thread_local uint32_t MixedGroupedGF256RowsForTesting = 0u;
static thread_local uint32_t MixedGroupedGF256RowMaskForTesting = 0u;
static thread_local uint32_t MixedGroupedGF256ResidueHashSeedForTesting = 0u;
static const uint32_t kMixedGroupedGF256SeedBase = UINT32_C(0xb7e15162);
static thread_local MixedResidueBucketMode MixedResidueBucketModeForTesting =
    MixedResidueBucketMode::Automatic;

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

static bool SelectGroupedGF256ResidueSeed(
    uint32_t period,
    uint32_t step_count,
    uint32_t& selected_seed)
{
    uint32_t candidate = kMixedGroupedGF256SeedBase;
    for (uint32_t attempt = 0u; attempt < 1024u; ++attempt, ++candidate)
    {
        if (HasFullCycleMixedResidueSeed(period, step_count, candidate))
        {
            selected_seed = candidate;
            return true;
        }
    }
    return false;
}

static void DisableMixedScheduleExperiments()
{
    MixedIndependentExtensionResiduesForTesting = false;
    MixedGroupedGF256RowsForTesting = 0u;
    MixedGroupedGF256RowMaskForTesting = 0u;
}

static uint32_t CountSetBits(uint32_t value)
{
    uint32_t count = 0u;
    while (value != 0u) {
        value &= value - 1u;
        ++count;
    }
    return count;
}

static bool IsValidatedH16Period(uint32_t period)
{
    return period == 31u || period == 32u;
}

bool SetMixedCoefficientPeriodForTesting(uint32_t period)
{
    const uint32_t H =
        MixedGF256RowsForTesting + MixedGF16RowsForTesting;
    if (period < H || period > kMixedCoefficientPeriod ||
        (MixedGF256RowsForTesting >= kMixedGF256Rows + 2u &&
         !IsValidatedH16Period(period)))
    {
        return false;
    }
    MixedCoefficientPeriodForTesting = period;
    MixedResidueSkewForTesting = 0u;
    MixedResidueScheduleForTesting = MixedResidueSchedule::Constant;
    DisableMixedScheduleExperiments();
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
    DisableMixedScheduleExperiments();
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
    DisableMixedScheduleExperiments();
    return true;
}

void SetMixedResidueHashSeedForTesting(uint32_t seed)
{
    MixedResidueHashSeedForTesting = seed;
    DisableMixedScheduleExperiments();
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
        (MixedGroupedGF256RowsForTesting != 0u ||
         MixedGeometryForTesting !=
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

bool SetMixedGroupedGF256RowsForTesting(uint32_t rows)
{
    if (rows > 9u ||
        (rows != 0u &&
         (MixedGeometryForTesting != MixedCoefficientGeometry::SharedCauchyX ||
          MixedGF256RowsForTesting != kMixedGF256Rows ||
          MixedGF16RowsForTesting != kMixedGF16Rows ||
          MixedCoefficientPeriodForTesting <=
              kMixedGF256Rows + kMixedGF16Rows ||
          MixedResidueSkewForTesting != 0u ||
          MixedResidueScheduleForTesting != MixedResidueSchedule::Constant ||
          MixedIndependentExtensionResiduesForTesting)))
    {
        return false;
    }
    if (rows != 0u &&
        !SelectGroupedGF256ResidueSeed(
            MixedCoefficientPeriodForTesting,
            MixedCoefficientPeriodForTesting -
                kMixedGF256Rows - kMixedGF16Rows,
            MixedGroupedGF256ResidueHashSeedForTesting))
    {
        return false;
    }
    MixedGroupedGF256RowsForTesting = rows;
    MixedGroupedGF256RowMaskForTesting = rows == 0u ? 0u :
        ((UINT32_C(1) << rows) - 1u) << (kMixedGF256Rows - rows);
    return true;
}

bool SetMixedGroupedGF256RowMaskForTesting(uint32_t row_mask)
{
    static const uint32_t kValidRowMask =
        (UINT32_C(1) << kMixedGF256Rows) - 1u;
    if ((row_mask & ~kValidRowMask) != 0u ||
        CountSetBits(row_mask) != MixedGroupedGF256RowsForTesting)
    {
        return false;
    }
    MixedGroupedGF256RowMaskForTesting = row_mask;
    return true;
}

bool SetMixedResidueBucketModeForTesting(MixedResidueBucketMode mode)
{
    if (mode != MixedResidueBucketMode::Automatic &&
        mode != MixedResidueBucketMode::Separate &&
        mode != MixedResidueBucketMode::Dual &&
        mode != MixedResidueBucketMode::JointDelta)
    {
        return false;
    }
    MixedResidueBucketModeForTesting = mode;
    return true;
}

MixedResidueBucketMode ActiveMixedResidueBucketModeForTesting()
{
    return MixedResidueBucketModeForTesting;
}

bool UseAutomaticMixedJointResidueBucketsForTesting(
    uint32_t block_count,
    uint32_t block_bytes,
    uint32_t coefficient_period)
{
    return UseAutomaticMixedJointResidueBuckets(
        block_count, block_bytes, coefficient_period);
}

bool SetMixedCoefficientGeometryForTesting(
    MixedCoefficientGeometry geometry)
{
    if (geometry != MixedCoefficientGeometry::FrozenPowerX &&
        geometry != MixedCoefficientGeometry::SharedCauchyX)
    {
        return false;
    }
    // Frozen geometry supports the base row count and the trimmed leading
    // subsets; only EXTRA rows (11/12) require shared-X coordinates.
    if (geometry != MixedCoefficientGeometry::SharedCauchyX &&
        MixedGF256RowsForTesting > kMixedGF256Rows)
    {
        return false;
    }
    MixedGeometryForTesting = geometry;
    if (geometry != MixedCoefficientGeometry::SharedCauchyX) {
        MixedResidueSkewForTesting = 0u;
        MixedResidueScheduleForTesting = MixedResidueSchedule::Constant;
    }
    DisableMixedScheduleExperiments();
    return true;
}

bool SetMixedGF16RowsForTesting(uint32_t rows)
{
    if (rows < kMixedGF16RowsMin || rows > kMixedGF16RowsMax ||
        (MixedGF256RowsForTesting >= kMixedGF256Rows + 2u &&
         rows != kMixedGF16RowsMax) ||
        MixedCoefficientPeriodForTesting <
            MixedGF256RowsForTesting + rows)
    {
        return false;
    }
    MixedGF16RowsForTesting = rows;
    MixedResidueSkewForTesting = 0u;
    MixedResidueScheduleForTesting = MixedResidueSchedule::Constant;
    DisableMixedScheduleExperiments();
    return true;
}

void SetMixedBandTrackingXForTesting(bool enabled)
{
    // Process-wide by design: the cost-model regime is fixed for a whole
    // benchmark invocation, and evaluator workers must all inherit the same
    // choice.  Atomic storage also makes an accidental late test change
    // race-free, although callers should configure it before launching work.
    MixedBandTrackingXOverride.store(
        enabled ? 1 : 0, std::memory_order_relaxed);
}

void ClearMixedBandTrackingXForTesting()
{
    MixedBandTrackingXOverride.store(-1, std::memory_order_relaxed);
}

bool SetMixedGF256RowsForTesting(uint32_t rows)
{
    // Row trims (one through nine rows) take the leading subset of the frozen
    // ten-row Cauchy table: Y stays in [0, rows) against the unchanged
    // X = 12 + residue coordinates, so the trimmed system remains one Cauchy
    // matrix in either geometry.  Extra rows (11/12) still require shared-X
    // coordinates as before.
    if (rows < kMixedGF256RowsMin || rows > kMixedGF256RowsMax ||
        MixedCoefficientPeriodForTesting < rows + MixedGF16RowsForTesting ||
        (rows >= kMixedGF256Rows + 2u &&
         (!IsValidatedH16Period(MixedCoefficientPeriodForTesting) ||
          MixedGF16RowsForTesting != kMixedGF16RowsMax)) ||
        (rows > kMixedGF256Rows &&
         MixedGeometryForTesting != MixedCoefficientGeometry::SharedCauchyX))
    {
        return false;
    }
    MixedGF256RowsForTesting = rows;
    MixedResidueSkewForTesting = 0u;
    MixedResidueScheduleForTesting = MixedResidueSchedule::Constant;
    DisableMixedScheduleExperiments();
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
    uint32_t residue = column % period;
    const uint32_t shift =
        ActiveMixedResidueBlockShift(column / period);
    CAT_DEBUG_ASSERT(shift < period);
    residue += shift;
    // Every block-shift schedule returns a reduced residue, so the sum is
    // below 2 * period and one subtraction avoids a second integer divide.
    return residue < period ? residue : residue - period;
}

uint32_t ActiveMixedExtensionCoefficientResidue(uint32_t column)
{
    const uint32_t period = ActiveMixedCoefficientPeriod();
    uint32_t residue = column % period;
    const uint32_t shift =
        ActiveMixedExtensionResidueBlockShift(column / period);
    CAT_DEBUG_ASSERT(shift < period);
    residue += shift;
    return residue < period ? residue : residue - period;
}

uint32_t ActiveMixedGroupedGF256CoefficientResidue(
    uint32_t column,
    uint32_t first_heavy_column)
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    if (MixedGroupedGF256RowsForTesting != 0u &&
        column < first_heavy_column)
    {
        const uint32_t period = ActiveMixedCoefficientPeriod();
        uint32_t residue = column % period;
        residue += ActiveMixedGroupedGF256ResidueBlockShift(column / period);
        return residue < period ? residue : residue - period;
    }
#else
    (void)first_heavy_column;
#endif
    return ActiveMixedCoefficientResidue(column);
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

uint32_t ActiveMixedGroupedGF256ResidueBlockShift(uint32_t block_index)
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    if (MixedGroupedGF256RowsForTesting == 0u) {
        return ActiveMixedResidueBlockShift(block_index);
    }
    const uint32_t period = MixedCoefficientPeriodForTesting;
    const uint32_t step_count =
        period - kMixedGF256Rows - kMixedGF16Rows;
    if (step_count == 0u) return 0u;
    static const uint32_t kStepCycle = 127u;
    static thread_local uint32_t cached_period = 0u;
    static thread_local uint32_t cached_seed = UINT32_MAX;
    static thread_local uint64_t prefix[kStepCycle + 1u] = {};
    if (cached_period != period ||
        cached_seed != MixedGroupedGF256ResidueHashSeedForTesting)
    {
        const uint32_t seed = MixedGroupedGF256ResidueHashSeedForTesting;
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
        cached_seed = seed;
    }
    const uint64_t complete_cycles = block_index / kStepCycle;
    const uint32_t remainder = block_index % kStepCycle;
    return (uint32_t)(
        (complete_cycles * prefix[kStepCycle] + prefix[remainder]) % period);
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

uint32_t ActiveMixedGroupedGF256Rows()
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    return MixedGroupedGF256RowsForTesting;
#else
    return 0u;
#endif
}

uint32_t ActiveMixedGroupedGF256RowMask()
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    return MixedGroupedGF256RowsForTesting != 0u ?
        MixedGroupedGF256RowMaskForTesting : 0u;
#else
    return 0u;
#endif
}

uint32_t ActiveMixedGroupedGF256HashSeed()
{
#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
    return MixedGroupedGF256RowsForTesting != 0u ?
        MixedGroupedGF256ResidueHashSeedForTesting : 0u;
#else
    return 0u;
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
