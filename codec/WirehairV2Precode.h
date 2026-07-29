#pragma once

#include "WirehairV2GF16.h"

#include <stdint.h>

#include <vector>

/*
    Certified precode construction for the V2 codec (wirehair-axd).

    Implements the Phase-B-certified `ldpcdense_s<S>_d12_s2_h12` legacy
    structure and the versioned SmallBandD4 architecture used by the stable
    dispatch-v1 benchmark target:

    - LegacyD12 uses S = GetDenseCount(K).  SmallBandD4 instead uses
      SmallBandStaircaseCount(K), which differs only at K <= 100.  Each
      SOURCE column connects to N1 distinct staircase parities (N1 = 2 below
      K=10000, N1 = 3 from K=10000 upward); parity row j additionally
      references its own parity column K+j and the staircase link K+j-1.
    - LegacyD12 uses D2 = 12 dense binary rows; SmallBandD4 uses D2 = 4.
      The rows span all K + S + D2 binary columns,
      generated Shuffle-2 style: first row = the ceil(span/2) set-half of a
      shuffled deck, every subsequent row = previous row XOR two deck-driven
      flips (one set-half entry, one clear-half entry), with a reshuffle at
      the halves boundary.  Consecutive rows differ in exactly two columns.
      In the experiment-only two-anchor variant, row 7 is instead reset to
      the balanced half of the second shuffled deck; rows 8..11 resume the
      two-flip cadence.  Segmented experiments instead reset anchors at
      {0,4,8}, {0,5,9}, or {0,3,6,9} and retain that cadence within segments.
    - H = 12 explicit Cauchy heavy rows.  Coefficients come from a Cauchy
      construction that is exactly MDS within any window of up to
      256 - H = 244 GE columns (W99 band requirement is >= 20); columns wrap
      modulo 244, so rank guarantees hold per-window, not globally.

    This module builds STRUCTURE only (column lists / coefficients); the
    solver and block-data phases consume it.  RNG is the production
    PCGRandom / ShuffleDeck16 idiom, NOT the simulator's splitmix stream:
    validation is structural invariants plus measured decode failure, not
    bit-identity with precode_sim.
*/

namespace wirehair_v2 {

// The opt-in two-anchor profile keeps the certified dense construction below
// this boundary and replaces only the second D12 shuffle half at larger K.
static const uint32_t kDenseTwoAnchorMinBlockCount = 4096u;

enum class HeavyCoefficientFamily : uint32_t
{
    PeriodicCauchy = 0,
    /// Experiment-only full-column hash family.  Named/public profiles never
    /// select this; it distinguishes a periodic-coefficient artifact from the
    /// generic GF(256) square-completion floor in actual-solver trials.
    HashedNonzero = 1
};

/**
    Experiment-only segmented D=12 dense-row layouts.

    Each nonzero layout replaces sparse inter-segment transitions with fresh
    balanced anchors at the named row indices.  Rows inside a segment retain
    the Shuffle-2 one-set/one-clear flip cadence.  Named/public profiles never
    select these values.
*/
enum class DenseAnchorLayout : uint32_t
{
    Disabled = 0,
    Three048 = 1, ///< anchors at rows {0, 4, 8}
    Three059 = 2, ///< anchors at rows {0, 5, 9}
    Four0369 = 3 ///< anchors at rows {0, 3, 6, 9}
};

/**
    Sentinel for PrecodeParams::StaircaseDegreeScale meaning "not supplied".

    Negative rather than zero because ZERO IS A LEGAL SCALE: a learned target
    mean row degree of 0.0 asks for an empty staircase edge budget, and the
    search is allowed to walk there and be rejected by its own failure
    constraint rather than by a floor built into the codec.

    When the scale is unset the certified legacy rule supplies the budget
    instead (K * min(SourceHits, S)).  That rule is the one thing this
    parameter replaces, so it must not run when a scale IS supplied.
*/
static const double kStaircaseDegreeScaleUnset = -1.0;

/**
    Accepted domain of a SUPPLIED PrecodeParams::StaircaseDegreeScale.

    The scale is an absolute mean row degree, so its domain is exactly the
    realizable range of a row degree: zero through the largest block count the
    codec accepts.  This is a validity domain, not a model -- nothing here
    picks a value, and the construction clamps the resolved budget to the
    K*S edges a simple bipartite graph can hold.
*/
static const double kStaircaseDegreeScaleMin = 0.0;
static const double kStaircaseDegreeScaleMax = 64000.0;

/// Largest per-source-column staircase degree the system can carry.
/// ValidatePrecodeSystem tallies per-column hits in a byte; mixture
/// construction rejects a larger degree before building, and validation also
/// guards the increment so malformed external systems fail rather than wrap.
static const uint32_t kStaircaseColumnHitsMax = 255u;

struct PrecodeParams
{
    uint32_t BlockCount = 0;   ///< K: source blocks
    uint32_t Staircase = 0;    ///< S: staircase parity columns
    uint32_t DenseRows = 0;    ///< D2: Shuffle-2 dense binary rows
    uint32_t HeavyRows = 0;    ///< H: Cauchy heavy rows
    uint32_t SourceHits = 0;   ///< N1: legacy hits when degree scale is unset
    CompletionField Field = CompletionField::GF256;

    /**
        Experiment-only degree-balanced staircase matching.

        The certified construction independently places each source column's
        distinct hits, so staircase-row source degrees fluctuate.  This option
        instead gives every row floor(EdgeCount/S) or ceil(EdgeCount/S) source
        columns while keeping every source degree, edge count, parity/link
        column, and downstream dense RNG stream unchanged.  A separately keyed
        randomized capacity matching avoids affine/cyclic placement.  Named
        profiles never enable this flag.
    */
    bool DegreeBalancedStaircase = false;

    /**
        Experiment-only TARGET MEAN STAIRCASE ROW DEGREE.

        This is an absolute number of source columns per staircase row.  It is
        not a multiple of anything, not a fraction of K, and not derived from
        S, H or D2.  With S staircase rows the total source->staircase edge
        budget is the definitional identity total = mean x count:

            edge_count = llround(StaircaseDegreeScale * S)

        clamped only at the K*S edges a simple bipartite graph can hold.  K
        does not appear, and neither does SourceHits.

        WHY IT EXISTS.  The certified construction pins the edge total at
        K * min(N1,S), because every source column carries EXACTLY min(N1,S)
        parities.  N1 is an integer, so that budget moves only in steps of K
        edges -- 128 at K=128 -- which is why N1 is near-frozen in any search
        over this structure (87.7% of mirrored pairs tied) and why a shaped
        row-degree distribution can only ever redistribute a budget it cannot
        resize.  A learned mean row degree moves the same budget by ONE edge
        at its finest step.

        COLUMN SIDE.  The edge total counted row-wise equals the total counted
        column-wise, so a freely chosen budget cannot leave every column at a
        constant hit count.  The exact-regularity invariant generalizes to an
        exact TWO-VALUED MIXTURE (see StaircaseDegreeMixture): every source
        column carries either lo or lo+1 parities, with exactly n_hi columns
        at lo+1.  lo may legitimately be ZERO -- the budget is not floored up
        to cover every column, because that would be a K-derived rule.  This
        is still a strict invariant and ValidatePrecodeSystem enforces it
        exactly, counts included.

        kStaircaseDegreeScaleUnset (the default, and what MakeCertifiedParams
        sets) means NOT SUPPLIED: the legacy K * min(N1,S) rule provides the
        budget instead and the construction is the certified one, bit for bit.
        Setting the scale to K*min(N1,S)/S reproduces that construction
        exactly -- same budget, n_hi = 0, every column at min(N1,S), and the
        identical PRNG draw order, because the mixture is exactly regular and
        no extra draw is made.  Named and public profiles never set it.
    */
    double StaircaseDegreeScale = kStaircaseDegreeScaleUnset;

    /**
        Identity-corner dense variant: the Shuffle-2 deck spans only the
        K + S source/staircase columns and dense row r additionally
        references exactly its own dense column K + S + r.

        The certified construction (false) decks over ALL K + S + D2
        binary columns, which makes the D2 x D2 dense-column corner rank
        ~1 (consecutive rows differ in just 2 columns), so a phased
        encoder can essentially never solve for the dense parity values
        directly — measured 0/2000 feasible seeds at K >= 1000.  With the
        identity corner each dense parity is simply its row's known-column
        sum (encoder-feasible by construction, same 2-XOR incremental
        generation).  The variant cleared paired 20k reliability comparisons
        through K=64000, but remains experimental: it changes the system, is
        unavailable at K=2..5 under this cadence, and the version-4 joint
        packet/precode solver already makes the certified full-span system
        encoder-feasible.
    */
    bool DenseIdentityCorner = false;

    /**
        D=12 Shuffle-2 variant with two balanced anchors.

        The certified construction keeps one balanced row and derives all
        eleven later rows through two-column flips.  This variant keeps rows
        0..6 unchanged, resets row 7 to the balanced half of the already
        scheduled second deck, and derives rows 8..11 through two-column
        flips.  It therefore trades one sparse difference direction for a
        second independently shuffled dense equation without adding a row.

        The flag records the active equation construction for one K.  Named
        profiles bind any adaptive K cutoff separately.  It is deliberately
        incompatible with DenseIdentityCorner and with DenseRows != 12.
    */
    bool DenseTwoAnchor = false;

    /**
        Experiment-only phase within the independently shuffled second deck.

        Zero emits the second balanced anchor before consuming any swap pair,
        preserving the versioned two-anchor construction exactly.  Phases one
        and two consume that many balanced set/clear swap pairs before row 7,
        then continue with four fresh pairs.  This changes no row count, RNG
        consumption, or production default.  Nonzero values are valid only
        when DenseTwoAnchor is enabled.
    */
    uint32_t DenseTwoAnchorPhase = 0;

    /**
        Experiment-only segmented three/four-anchor layout.

        This is mutually exclusive with DenseTwoAnchor and
        DenseIdentityCorner and is valid only for DenseRows == 12.  It does
        not change the staircase, dense-row, or heavy-row counts.
    */
    DenseAnchorLayout SegmentedDenseAnchors = DenseAnchorLayout::Disabled;

    HeavyCoefficientFamily HeavyFamily =
        HeavyCoefficientFamily::PeriodicCauchy;

    uint64_t Seed = 0;         ///< constraint-generation seed
};

/// Largest block count that takes the measured small-band staircase rule.
/// Above it the inherited wirehair::GetDenseCount() value is used unchanged.
/// The reliability sweep behind the rule covers K = 2..100 only, so the
/// boundary sits at the top of that measured band; S is discontinuous across
/// it (S = 12 at K = 100, S = 26 at K = 101), as the inherited table already
/// is at its own K = 64/65 seam.
static const uint32_t kSmallBandStaircaseMaxBlockCount = 100u;

/**
    Staircase parity count S for a block count.

    Equals wirehair::GetDenseCount(block_count) above
    kSmallBandStaircaseMaxBlockCount.  At or below it, S is the smaller of the
    inherited value and the exact integer value floor(5 * sqrt(K) / 4), never
    below 2.  See the definition for the reliability sweep behind the
    coefficient.  Returns 0 for block counts outside [2, 64000].

    The frozen public profiles continue to select LegacyD12.  The stable
    benchmark-only dispatch-v1 contract selects SmallBandD4 explicitly, and
    its all-K fingerprint pins that equation system.  SeedProfile::DenseCount
    intentionally remains the legacy table output: MatrixSeedFromProfile(),
    packet seed derivation, and legacy seed/fixup tables consume it as a seed
    input.  The selected staircase count and dense-row count are separate
    versioned architecture fields and must not be inferred from DenseCount.
*/
uint32_t SmallBandStaircaseCount(uint32_t block_count);

/// Construction helper used by experiments: S = SmallBandStaircaseCount(K),
/// D2 = 12, H = 12, N1 = 2 below K=10000 and N1 = 3 from K=10000 upward.
/// Message-codec contracts resolve their architecture after this helper.
PrecodeParams MakeCertifiedParams(uint32_t block_count, uint64_t seed);

/// Versioned mixed 10-row GF(256) + 2-row GF(2^16) completion rule.
PrecodeParams MakeMixedParams(uint32_t block_count, uint64_t seed);

/**
    Exact two-valued source-column degree mixture implied by a parameter set.

    The staircase edge total counted row-wise must equal the total counted
    column-wise, so a freely chosen edge budget necessarily moves the column
    side off a single integer.  It moves it by the smallest possible amount:
    to a mixture of exactly two ADJACENT degrees.

        EdgeCount   = min(llround(scale * S), K*S)   [scale supplied]
        EdgeCount   = K * min(SourceHits, S)         [scale unset: legacy]
        LowHits     = EdgeCount / K                (floor; may be 0)
        HighColumns = EdgeCount - LowHits * K      (columns at LowHits + 1)
        LowColumns  = K - HighColumns

    HighColumns == 0 is the exactly-regular case, which is what the legacy
    budget always produces -- the certified construction unchanged.
*/
struct StaircaseDegreeMixture
{
    uint32_t EdgeCount = 0;   ///< total source->staircase edges
    uint32_t LowHits = 0;     ///< parities carried by a LowColumns column
    uint32_t HighHits = 0;    ///< LowHits + 1; meaningful when HighColumns > 0
    uint32_t LowColumns = 0;  ///< columns carrying LowHits parities
    uint32_t HighColumns = 0; ///< columns carrying HighHits parities
    uint32_t MaxHits = 0;     ///< largest realized column degree
};

/**
    Resolve the column-degree mixture for a parameter set.

    Returns false when the block count, staircase count, source hits or scale
    are outside their domains, or when the budget cannot be realized as a
    two-valued mixture (a column cannot touch a staircase row twice, so the
    top degree must stay within [0, min(S, kStaircaseColumnHitsMax)]).  It is
    a pure function of `params` -- no PRNG, no thread-local state -- so the
    builder and the validator agree by construction.
*/
bool MakeStaircaseDegreeMixture(
    const PrecodeParams& params,
    StaircaseDegreeMixture& out);

/**
    True when no non-receipted staircase shape or realized-row override is
    active.

    The target-mean StaircaseDegreeScale override is deliberately excluded:
    stable benchmark targets receipt that overlay explicitly.  Both the
    thread-local and environment-backed shape/row-sequence controls are
    equation-active but are not part of that receipt, so stable target
    identities must reject them.
*/
bool IsCanonicalStableTargetStaircaseState();

/**
    True only when no target-mean staircase scale overlay is active.

    Stable benchmark targets may receipt that overlay separately, which is why
    IsCanonicalStableTargetStaircaseState() deliberately permits it.  Frozen
    equation contracts that do not name such an overlay use this stricter
    predicate as well.
*/
bool IsCanonicalStaircaseDegreeScaleState();

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
/**
    Compact semantic fingerprint of every ambient precode hook that can alter
    the supplied parameter set's equations.

    This includes the effective staircase scale/shape/realized degrees and,
    for mixed completion, the complete active row split, coefficient geometry,
    and residue schedules.  The process-wide tracking-X regime has a dedicated
    single-read field in ExplicitEquationStateIdentityForTesting, so it is
    deliberately excluded here.  Execution-only controls such as the
    residue-bucket implementation are also excluded.
*/
uint64_t ActivePrecodeEquationStateFingerprintForTesting(
    const PrecodeParams& params);

/**
    Return the process-wide/environment tracking-X selection, ignoring any
    scoped explicit-operation snapshot on the current thread.
*/
bool AmbientMixedBandTrackingXForTesting();

/**
    Pin the effective tracking-X selection on this thread for one complete
    explicit operation.

    The benchmark control remains process-wide, but an already-pinned explicit
    endpoint must not assemble equations from two states if another worker
    changes that control concurrently.  Nesting is supported so decoder retry
    recursion retains the outer snapshot.
*/
class ScopedMixedBandTrackingXForTesting
{
public:
    ScopedMixedBandTrackingXForTesting(
        bool active,
        bool enabled) noexcept;
    ~ScopedMixedBandTrackingXForTesting() noexcept;

    ScopedMixedBandTrackingXForTesting(
        const ScopedMixedBandTrackingXForTesting&) = delete;
    ScopedMixedBandTrackingXForTesting& operator=(
        const ScopedMixedBandTrackingXForTesting&) = delete;

private:
    int Previous = -1;
};
#endif

struct PrecodeSystem
{
    PrecodeParams Params = {};

    /**
        Staircase parity rows.

        Row j (j in [0, S)) is a GF(2) constraint over binary columns
        [0, K + S): the source columns whose staircase hits landed on parity
        j, plus the own-parity column K + j, plus the staircase link column
        K + j - 1 for j > 0.  Column lists are sorted and deduplicated
        (every column appears exactly once; construction cannot produce
        duplicates because per-source hits are distinct).  The optional
        degree-balanced experiment additionally makes every row's source
        degree floor(EdgeCount/S) or ceil(EdgeCount/S), where EdgeCount is
        the resolved StaircaseDegreeMixture budget.
    */
    std::vector<std::vector<uint32_t>> StaircaseRows;

    /**
        Shuffle-2 dense rows.

        Row r (r in [0, D2)) is a GF(2) constraint over binary columns
        [0, K + S + D2), stored as a sorted column list.  In the certified
        construction, row 0 has exactly ceil((K+S+D2)/2) columns and every
        row r > 0 differs from row r - 1 in exactly two columns.  In the
        identity-corner variant, the deck spans only K+S known columns and
        each row additionally owns dense column K+S+r; consecutive full rows
        differ in four columns, but their known-column part still changes by
        exactly two deck flips.  In the two-anchor variant, row 7 is another
        balanced row and only the row 6 -> 7 transition is exempt from the
        two-column-difference rule.  In a segmented-anchor experiment, every
        listed anchor transition is exempt and all within-segment transitions
        still differ in exactly two columns.

        Known limitation (inherited from the certified reference
        construction): at tiny EVEN spans (K=2 and K=4 with the certified
        table) a later row's weight can walk down to exactly zero — a dead
        constraint — at roughly 1e-4 systems.  Guarding would break the
        exact-2-difference invariant.  Version-4 message initialization
        instead rejects rank-deficient attempts and deterministically selects
        a full-rank joint packet/precode seed; exhaustive tiny-K tests pin the
        selected attempts.
    */
    std::vector<std::vector<uint32_t>> DenseRowColumns;
};

/**
    Build the staircase + Shuffle-2 constraint structure.

    Returns false for invalid parameters (BlockCount outside [2, 64000],
    Staircase == 0, SourceHits outside [1, 8], DenseRows > 64,
    HeavyRows > 128, an invalid dense experiment combination, or a full
    symbol domain that does not fit uint16) or if the generated structure
    fails ValidatePrecodeSystem().
*/
bool BuildPrecodeSystem(const PrecodeParams& params, PrecodeSystem& out);

/**
    Validate every structural invariant consumed by the encoder.  Validation
    uses widened arithmetic and performs no writes to block data.
*/
bool ValidatePrecodeSystem(const PrecodeSystem& system);

/**
    GF(256) coefficient of heavy row r at GE column c.

    Cauchy element 1 / (X ^ Y) with Y = r and X = H + (c mod (256 - H)):
    nonzero everywhere, and every square submatrix drawn from one
    244-column window is invertible.  heavy_rows must be <= 128.  Direct
    callers must initialize GF256 first (`wirehair_init()` or `gf256_init()`);
    all encoder/decoder/solver entry points do this internally.
*/
uint8_t HeavyCoefficient(
    uint32_t heavy_row,
    uint32_t ge_column,
    uint32_t heavy_rows);

static const uint32_t kMixedPackedCoefficientWords =
    (kMixedGF256RowsMax + kMixedGF16RowsMax + 3u) / 4u;

enum class MixedCoefficientGeometry : uint32_t
{
    FrozenPowerX = 0,
    SharedCauchyX = 1
};

enum class MixedResidueSchedule : uint32_t
{
    Constant = 0,
    Ramp = 1,
    Hashed = 2
};

/// Immutable row-major coefficient period shared by mixed encode, solve, and
/// verification paths.  The accessor initializes both fields and returns null
/// only if their arithmetic tables could not be initialized.
struct MixedCoefficientRows
{
    uint8_t Subfield[kMixedGF256RowsMax][kMixedCoefficientPeriod];
    uint16_t Extension[kMixedGF16RowsMax][kMixedCoefficientPeriod];
};

const MixedCoefficientRows* GetMixedCoefficientRows();

/// The same mixed coefficient period packed as four uint16 lanes per word,
/// residue-major.  This separate lazy cache prevents encode/verify callers
/// from paying the solve-only packing cost on first use.
struct MixedPackedCoefficients
{
    uint64_t ByResidue[
        kMixedCoefficientPeriod][kMixedPackedCoefficientWords];
};

const MixedPackedCoefficients* GetMixedPackedCoefficients();

/**
    Active mixed coefficient period.

    Production builds always return kMixedCoefficientPeriod.  Test builds may
    select a smaller period on the current thread to measure the speed/rank
    tradeoff without changing any named or serialized wire profile.
*/
uint32_t ActiveMixedCoefficientPeriod();
/// Coefficient residue for a GE column under the active balanced block skew.
uint32_t ActiveMixedCoefficientResidue(uint32_t column);
/// Residue used by extension-field rows.  Production aliases the shared
/// residue; test builds may select an independently keyed schedule.
uint32_t ActiveMixedExtensionCoefficientResidue(uint32_t column);
/**
    Residue used by the experiment-only grouped GF(256) suffix rows.

    Non-corner columns use one fixed full-cycle schedule C.  The final H
    completion columns deliberately remain on canonical schedule A so the
    encoder corner is byte-for-byte unchanged.  Production and a zero active
    grouped-row count always alias ActiveMixedCoefficientResidue().
*/
uint32_t ActiveMixedGroupedGF256CoefficientResidue(
    uint32_t column,
    uint32_t first_heavy_column);
/// Rotation applied to one complete period-sized block.
uint32_t ActiveMixedResidueBlockShift(uint32_t block_index);
uint32_t ActiveMixedExtensionResidueBlockShift(uint32_t block_index);
uint32_t ActiveMixedGroupedGF256ResidueBlockShift(uint32_t block_index);
uint32_t ActiveMixedResidueSkew();
MixedResidueSchedule ActiveMixedResidueSchedule();
uint32_t ActiveMixedResidueHashSeed();
bool ActiveMixedResiduesRotated();
bool ActiveMixedIndependentExtensionResidues();
uint32_t ActiveMixedGroupedGF256Rows();
/// Cauchy-Y rows assigned to grouped schedule C (bit r selects Y=r).
uint32_t ActiveMixedGroupedGF256RowMask();
uint32_t ActiveMixedGroupedGF256HashSeed();
MixedCoefficientGeometry ActiveMixedCoefficientGeometry();
uint32_t ActiveMixedGF256Rows();
uint32_t ActiveMixedGF16Rows();
uint32_t ActiveMixedPackedCoefficientWords();

/**
    True when every equation-active mixed-completion control matches the
    frozen production 10+2/P244 geometry.

    Inactive keyed seeds and execution-only bucket policy are intentionally
    excluded.  Stable raw contracts call this both when binding and replaying
    a profile so thread-local test-hook changes cannot inherit their identity.
*/
bool IsCanonicalMixedCompletionState();

// Production data-plane cap for joint A/B residue accumulation.  The small
// scheduling vectors are bounded by K and P separately; this limit covers the
// three P * block_bytes planes that dominate scratch use.
static const uint64_t kMixedJointResidueBucketDataByteCap =
    UINT64_C(64) << 20;

/// True when three joint-delta data planes fit the supplied scratch budget.
bool MixedJointResidueBucketStorageFits(
    uint32_t coefficient_period,
    uint32_t block_bytes,
    uint64_t data_byte_limit);

/**
    Production implementation policy for independently scheduled A/B mixed
    residues.  This is an execution choice only: it changes neither equations
    nor wire/profile bytes.  Pinned ABBA measurements currently justify the
    joint-delta helper only at P=32 and the two conservative K/payload
    crossovers below; callers must additionally enforce the scratch cap.
*/
bool UseAutomaticMixedJointResidueBuckets(
    uint32_t block_count,
    uint32_t block_bytes,
    uint32_t coefficient_period);

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
enum class MixedResidueBucketMode : uint32_t
{
    Automatic = 0,
    Separate = 1,
    Dual = 2,
    JointDelta = 3
};

/// Strict parsers shared by the environment hooks and their focused tests.
/// Lists are comma-separated with optional surrounding ASCII whitespace.
/// Invalid, non-finite, overflowing, partially parsed, or empty inputs return
/// false without modifying the output.
bool ParseStaircaseDegreesForTesting(
    const char* text,
    std::vector<double>& weights);
bool ParseStaircaseRowDegreesForTesting(
    const char* text,
    std::vector<uint32_t>& degrees);
bool ParseStaircaseDegreeScaleForTesting(
    const char* text,
    double& scale);

/// Per-thread staircase row-degree distribution (weights for degrees 1,2,3...).
/// Empty clears the override and restores the stock construction.  Thread-local
/// so concurrent workers can evaluate different shapes.
void SetStaircaseDegreesForTesting(const std::vector<double>& weights);
void ClearStaircaseDegreesForTesting();

/// Per-thread override of PrecodeParams::StaircaseDegreeScale -- the target
/// mean staircase row degree -- in [kStaircaseDegreeScaleMin,
/// kStaircaseDegreeScaleMax].  Set wins over the params field and over
/// WIREHAIR_V2_STAIRCASE_DEGREE_SCALE; Clear restores the params field (or the
/// environment value, if one is set).  Thread-local, matching
/// SetStaircaseDegreesForTesting, so concurrent search workers can each
/// evaluate a different scale.  Passing the exact
/// kStaircaseDegreeScaleUnset sentinel is equivalent to Clear; every other
/// negative value remains an active invalid override and fails construction.
void SetStaircaseDegreeScaleForTesting(double scale);
void ClearStaircaseDegreeScaleForTesting();

/// Per-thread override of the REALIZED staircase row-degree sequence: exactly
/// S integers, one per staircase row in row order, summing to the resolved edge
/// budget.  This names the object the shape knob only samples -- the sampler
/// draws S degrees from a distribution and rescales them onto the budget, so
/// one shape produces a different sequence on every construction seed -- and
/// lets a sequence be evaluated as itself.  Empty clears the override.
///
/// The matching stream still consumes the S draws the sampler would have made,
/// so a pinned sequence reproduces the shaped path CONDITIONED on that sequence
/// bit for bit rather than landing on a different random stream.
///
/// Wins over SetStaircaseDegreesForTesting and over
/// WIREHAIR_V2_STAIRCASE_ROW_DEGREES; a sequence whose length is not S, whose
/// entries are not in [0, K], or whose sum is not the resolved edge budget is a
/// REFUSED construction, not a silently repaired one.
void SetStaircaseRowDegreesForTesting(const std::vector<uint32_t>& degrees);
void ClearStaircaseRowDegreesForTesting();

/// Set the current thread's experiment-only period in [H, 244].
bool SetMixedCoefficientPeriodForTesting(uint32_t period);
/// Rotate each period block by a corner-preserving skew in [0, P-H].
bool SetMixedResidueSkewForTesting(uint32_t skew);
/// Select a constant, ramp, or hashed corner-safe residue schedule.
bool SetMixedResidueScheduleForTesting(MixedResidueSchedule schedule);
/// Select the deterministic hashed-schedule sequence (benchmark/test only).
void SetMixedResidueHashSeedForTesting(uint32_t seed);
/// Derive and activate the first K-keyed hashed schedule with a full cycle.
bool SelectFullCycleMixedResidueKeyedSeedForTesting(
    uint32_t base_seed,
    uint32_t block_count,
    uint32_t& selected_seed);
/// Give GF(2^16) rows an independently keyed full-cycle hashed schedule.
bool SetMixedIndependentExtensionResiduesForTesting(bool enabled);
/// Select the XOR used to derive the independent extension schedule seed.
void SetMixedIndependentExtensionSeedXorForTesting(uint32_t seed_xor);
///
/// Assign one shared schedule C to the final 0..9 rows of the fixed H12
/// GF(256) prefix.  This raw-architecture hook is intentionally restricted to
/// shared-X, constant-A, 10 GF(256) + 2 GF(2^16), P>H configurations and is
/// mutually exclusive with independent GF(2^16) scheduling.
bool SetMixedGroupedGF256RowsForTesting(uint32_t rows);
///
/// Select which of the ten fixed Cauchy-Y rows consume schedule C.  The mask
/// must contain exactly ActiveMixedGroupedGF256Rows() bits in the low ten
/// positions.  Internally this is only an equation-row permutation: the
/// existing compact A/C destination groups, source scans, and arithmetic work
/// remain identical to the suffix selected by SetMixedGroupedGF256RowsForTesting().
bool SetMixedGroupedGF256RowMaskForTesting(uint32_t row_mask);
/// Select the secondary-schedule RHS accumulation implementation.
bool SetMixedResidueBucketModeForTesting(MixedResidueBucketMode mode);
MixedResidueBucketMode ActiveMixedResidueBucketModeForTesting();
/// Compatibility test accessor for the production implementation policy.
bool UseAutomaticMixedJointResidueBucketsForTesting(
    uint32_t block_count,
    uint32_t block_bytes,
    uint32_t coefficient_period);
/// Select frozen or shared-X mixed coefficients on the current test thread.
bool SetMixedCoefficientGeometryForTesting(MixedCoefficientGeometry geometry);
/// Select a leading 1..10-row subset of the frozen GF(256) table, or guarded
/// shared-X 11/12-row experiments (12 requires the validated four-row
/// extension geometry).
bool SetMixedGF256RowsForTesting(uint32_t rows);
/// Select zero through four extension rows; twelve GF(256) rows require four.
bool SetMixedGF16RowsForTesting(uint32_t rows);
/// Force every worker onto the cost model's trimmed-band X-coordinate regime.
/// Process-wide: configure before launching workers.  Clear restores the
/// WIREHAIR_V2_BAND_TRACKING_X environment-controlled experiment.
void SetMixedBandTrackingXForTesting(bool enabled);
void ClearMixedBandTrackingXForTesting();
#endif

/// Coefficient dispatch for actual encoder/decoder equations.  Alternate
/// families are confined to explicitly constructed experiment systems.
uint8_t HeavyCoefficientForParams(
    const PrecodeParams& params,
    uint32_t heavy_row,
    uint32_t ge_column);

} // namespace wirehair_v2
