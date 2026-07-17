#pragma once

#include "WirehairV2GF16.h"

#include <stdint.h>

#include <vector>

/*
    Certified precode construction for the V2 codec (wirehair-axd).

    Implements the Phase-B-certified `ldpcdense_s<S>_d12_s2_h12` structure
    from experiments/precode (see experiments/precode/README.md and
    results/SUMMARY.md "Ship Gate Decision"):

    - S = GetDenseCount(K) LDPC-staircase parity columns.  Each SOURCE column
      connects to N1 distinct staircase parities (N1 = 2 below K=10000,
      N1 = 3 from K=10000 upward); parity row j additionally references its
      own parity column K+j and the staircase link K+j-1.
    - D2 = 12 dense binary rows over all span = K + S + D2 binary columns,
      generated Shuffle-2 style: first row = the ceil(span/2) set-half of a
      shuffled deck, every subsequent row = previous row XOR two deck-driven
      flips (one set-half entry, one clear-half entry), with a reshuffle at
      the halves boundary.  Consecutive rows differ in exactly two columns.
      In the experiment-only two-anchor variant, row 7 is instead reset to
      the balanced half of the second shuffled deck; rows 8..11 resume the
      two-flip cadence.
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

enum class HeavyCoefficientFamily : uint32_t
{
    PeriodicCauchy = 0,
    /// Experiment-only full-column hash family.  Named/public profiles never
    /// select this; it distinguishes a periodic-coefficient artifact from the
    /// generic GF(256) square-completion floor in actual-solver trials.
    HashedNonzero = 1
};

struct PrecodeParams
{
    uint32_t BlockCount = 0;   ///< K: source blocks
    uint32_t Staircase = 0;    ///< S: staircase parity columns
    uint32_t DenseRows = 0;    ///< D2: Shuffle-2 dense binary rows
    uint32_t HeavyRows = 0;    ///< H: Cauchy heavy rows
    uint32_t SourceHits = 0;   ///< N1: staircase parities per source column
    CompletionField Field = CompletionField::GF256;

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
        Experiment-only D=12 Shuffle-2 variant with two balanced anchors.

        The certified construction keeps one balanced row and derives all
        eleven later rows through two-column flips.  This variant keeps rows
        0..6 unchanged, resets row 7 to the balanced half of the already
        scheduled second deck, and derives rows 8..11 through two-column
        flips.  It therefore trades one sparse difference direction for a
        second independently shuffled dense equation without adding a row.

        Named/public profiles never select this flag.  It is deliberately
        incompatible with DenseIdentityCorner and with DenseRows != 12.
    */
    bool DenseTwoAnchor = false;

    HeavyCoefficientFamily HeavyFamily =
        HeavyCoefficientFamily::PeriodicCauchy;

    uint64_t Seed = 0;         ///< constraint-generation seed
};

/// Certified rule: S = GetDenseCount(K), D2 = 12, H = 12,
/// N1 = 2 below K=10000 and N1 = 3 from K=10000 upward
PrecodeParams MakeCertifiedParams(uint32_t block_count, uint64_t seed);

/// Versioned mixed 10-row GF(256) + 2-row GF(2^16) completion rule.
PrecodeParams MakeMixedParams(uint32_t block_count, uint64_t seed);

struct PrecodeSystem
{
    PrecodeParams Params = {};

    /**
        Staircase parity rows.

        Row j (j in [0, S)) is a GF(2) constraint over binary columns
        [0, K + S): the source columns whose N1 hits landed on parity j,
        plus the own-parity column K + j, plus the staircase link column
        K + j - 1 for j > 0.  Column lists are sorted and deduplicated
        (every column appears exactly once; construction cannot produce
        duplicates because per-source hits are distinct).
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
        two-column-difference rule.

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
/// Rotation applied to one complete period-sized block.
uint32_t ActiveMixedResidueBlockShift(uint32_t block_index);
uint32_t ActiveMixedExtensionResidueBlockShift(uint32_t block_index);
uint32_t ActiveMixedResidueSkew();
MixedResidueSchedule ActiveMixedResidueSchedule();
uint32_t ActiveMixedResidueHashSeed();
bool ActiveMixedResiduesRotated();
bool ActiveMixedIndependentExtensionResidues();
MixedCoefficientGeometry ActiveMixedCoefficientGeometry();
uint32_t ActiveMixedGF256Rows();
uint32_t ActiveMixedGF16Rows();
uint32_t ActiveMixedPackedCoefficientWords();

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
/// Select the independent-schedule RHS accumulation implementation.
bool SetMixedResidueBucketModeForTesting(MixedResidueBucketMode mode);
MixedResidueBucketMode ActiveMixedResidueBucketModeForTesting();
/// Compatibility test accessor for the production implementation policy.
bool UseAutomaticMixedJointResidueBucketsForTesting(
    uint32_t block_count,
    uint32_t block_bytes,
    uint32_t coefficient_period);
/// Select frozen or shared-X mixed coefficients on the current test thread.
bool SetMixedCoefficientGeometryForTesting(MixedCoefficientGeometry geometry);
/// Select 10/11 rows generally, or a validated 12+4-row test geometry.
bool SetMixedGF256RowsForTesting(uint32_t rows);
/// Select two, three, or four extension rows; twelve GF(256) rows require four.
bool SetMixedGF16RowsForTesting(uint32_t rows);
#endif

/// Coefficient dispatch for actual encoder/decoder equations.  Alternate
/// families are confined to explicitly constructed experiment systems.
uint8_t HeavyCoefficientForParams(
    const PrecodeParams& params,
    uint32_t heavy_row,
    uint32_t ge_column);

} // namespace wirehair_v2
