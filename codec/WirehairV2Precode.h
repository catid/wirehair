#pragma once

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS) && \
    defined(WIREHAIR_V2_ENABLE_BENCHMARK_EQUATIONS)
#error "WH2 test hooks and benchmark equations are mutually exclusive"
#endif

#include <stdint.h>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <stdexcept>
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
      the halves boundary.  Consecutive rows differ in exactly two columns,
      so encoder parity generation costs 2 block-XORs per row after the
      first (the measured 54% precode-gen cut at K=3200 for n12+s2).
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

/**
    Experiment-only segmentation of the D=12 Shuffle-2 dense equations.

    An anchor is a freshly shuffled balanced half-row.  Rows between anchors
    retain the certified one-set/one-clear two-column delta cadence.  The
    layouts are pure binary/GF(256) architecture arms; named/public profiles
    always use Disabled.
*/
enum class DenseAnchorLayout : uint32_t
{
    Disabled = 0,
    Two07 = 1,    ///< independent anchors at rows {0, 7}
    Four0369 = 2 ///< independent anchors at rows {0, 3, 6, 9}
};

struct PrecodeParams
{
    uint32_t BlockCount = 0;   ///< K: source blocks
    uint32_t Staircase = 0;    ///< S: staircase parity columns
    uint32_t DenseRows = 0;    ///< D2: Shuffle-2 dense binary rows
    uint32_t HeavyRows = 0;    ///< H: Cauchy heavy rows
    uint32_t SourceHits = 0;   ///< N1: staircase parities per source column
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

    HeavyCoefficientFamily HeavyFamily =
        HeavyCoefficientFamily::PeriodicCauchy;

    DenseAnchorLayout DenseAnchors = DenseAnchorLayout::Disabled;

    uint64_t Seed = 0;         ///< constraint-generation seed
};

/// Certified rule: S = GetDenseCount(K), D2 = 12, H = 12,
/// N1 = 2 below K=10000 and N1 = 3 from K=10000 upward
PrecodeParams MakeCertifiedParams(uint32_t block_count, uint64_t seed);

struct PrecodeRowView
{
    const uint32_t* Data = nullptr;
    size_t Count = 0u;

    const uint32_t* begin() const noexcept { return Data; }
    const uint32_t* end() const noexcept
    {
        return Count == 0u ? Data : Data + Count;
    }
    size_t size() const noexcept { return Count; }
    bool empty() const noexcept { return Count == 0u; }
    const uint32_t& operator[](size_t index) const noexcept
    {
        return Data[index];
    }
    const uint32_t& front() const noexcept { return Data[0]; }
    const uint32_t& back() const noexcept { return Data[Count - 1u]; }
};

struct MutablePrecodeRowView
{
    uint32_t* Data = nullptr;
    size_t Count = 0u;

    uint32_t* begin() const noexcept { return Data; }
    uint32_t* end() const noexcept
    {
        return Count == 0u ? Data : Data + Count;
    }
    size_t size() const noexcept { return Count; }
    bool empty() const noexcept { return Count == 0u; }
    uint32_t& operator[](size_t index) const noexcept { return Data[index]; }
    uint32_t& front() const noexcept { return Data[0]; }
    uint32_t& back() const noexcept { return Data[Count - 1u]; }
};

inline bool operator==(
    const PrecodeRowView& left,
    const PrecodeRowView& right) noexcept
{
    return left.size() == right.size() &&
        (left.empty() ||
         std::equal(left.begin(), left.end(), right.begin()));
}

inline bool operator!=(
    const PrecodeRowView& left,
    const PrecodeRowView& right) noexcept
{
    return !(left == right);
}

struct PrecodeSystem
{
    PrecodeParams Params = {};

    /**
        Flat binary-row storage shared by staircase rows followed by dense
        basis rows.  BinaryRowOffsets always uses uint32_t because the fully
        validated construction envelope contains far fewer than UINT32_MAX
        references.  A valid non-empty system has S + D2 + 1 offsets, starts
        at zero, and ends at BinaryRowColumns.size().

        Keeping these two vectors unconditional is important: V2 sources are
        compiled into archives with different test-hook macros, so conditional
        row-storage members would create an ODR/ABI mismatch.

        Staircase row j (j in [0, S)) contains the source columns whose N1
        hits landed on parity j, the own-parity column K + j, and for j > 0
        the staircase link K + j - 1.  Each later dense-basis segment starts
        with a balanced half-dense anchor; subsequent rows store sorted
        two-column deltas (four with the identity corner).  All slices are
        sorted and deduplicated.  These elementary row additions preserve the
        original equation span without retaining duplicate half-dense rows.

        The certified tiny-even-span limitation remains unchanged: a later
        reconstructed dense row can become zero for rare K=2/K=4 seeds.
        Message initialization rejects rank-deficient attempts and selects a
        full-rank packet/precode seed.
    */
    std::vector<uint32_t> BinaryRowOffsets;
    std::vector<uint32_t> BinaryRowColumns;

    /// Row accessors require HasValidPrecodeSystemShape() to have succeeded.
    size_t BinaryRowCount() const noexcept
    {
        return BinaryRowOffsets.empty() ? 0u :
            BinaryRowOffsets.size() - 1u;
    }

    PrecodeRowView BinaryRow(size_t row) const noexcept
    {
        const uint32_t first = BinaryRowOffsets[row];
        const uint32_t last = BinaryRowOffsets[row + 1u];
        PrecodeRowView view;
        view.Data = BinaryRowColumns.empty() ? nullptr :
            BinaryRowColumns.data() + first;
        view.Count = (size_t)(last - first);
        return view;
    }

    MutablePrecodeRowView MutableBinaryRow(size_t row) noexcept
    {
        const uint32_t first = BinaryRowOffsets[row];
        const uint32_t last = BinaryRowOffsets[row + 1u];
        MutablePrecodeRowView view;
        view.Data = BinaryRowColumns.empty() ? nullptr :
            BinaryRowColumns.data() + first;
        view.Count = (size_t)(last - first);
        return view;
    }

    PrecodeRowView StaircaseRow(size_t row) const noexcept
    {
        return BinaryRow(row);
    }

    MutablePrecodeRowView MutableStaircaseRow(size_t row) noexcept
    {
        return MutableBinaryRow(row);
    }

    PrecodeRowView DenseBasisRow(size_t row) const noexcept
    {
        return BinaryRow((size_t)Params.Staircase + row);
    }

    MutablePrecodeRowView MutableDenseBasisRow(size_t row) noexcept
    {
        return MutableBinaryRow((size_t)Params.Staircase + row);
    }
};

#if defined(WIREHAIR_V2_ENABLE_TEST_HOOKS)
namespace test {

inline std::vector<uint32_t> CopyPrecodeRowForTesting(
    const PrecodeRowView& row)
{
    if (row.empty()) {
        return {};
    }
    return std::vector<uint32_t>(row.begin(), row.end());
}

inline void ReplacePrecodeRowForTesting(
    PrecodeSystem& system,
    size_t row,
    const std::vector<uint32_t>& replacement)
{
    if (system.BinaryRowOffsets.empty() ||
        row >= system.BinaryRowOffsets.size() - 1u)
    {
        throw std::out_of_range("WH2 test precode row is out of range");
    }
    const uint32_t first = system.BinaryRowOffsets[row];
    const uint32_t last = system.BinaryRowOffsets[row + 1u];
    if (first > last || last > system.BinaryRowColumns.size()) {
        throw std::out_of_range("WH2 test precode offsets are malformed");
    }
    const size_t old_size = (size_t)(last - first);
    const size_t retained_size = system.BinaryRowColumns.size() - old_size;
    const size_t max_size = std::numeric_limits<uint32_t>::max();
    if (retained_size > max_size ||
        replacement.size() > max_size - retained_size)
    {
        throw std::length_error("WH2 test precode row is too large");
    }
    system.BinaryRowColumns.erase(
        system.BinaryRowColumns.begin() + first,
        system.BinaryRowColumns.begin() + last);
    system.BinaryRowColumns.insert(
        system.BinaryRowColumns.begin() + first,
        replacement.begin(),
        replacement.end());
    for (size_t offset = row + 1u;
         offset < system.BinaryRowOffsets.size();
         ++offset)
    {
        if (replacement.size() >= old_size) {
            system.BinaryRowOffsets[offset] +=
                (uint32_t)(replacement.size() - old_size);
        }
        else {
            system.BinaryRowOffsets[offset] -=
                (uint32_t)(old_size - replacement.size());
        }
    }
}

inline void ErasePrecodeRowForTesting(PrecodeSystem& system, size_t row)
{
    if (system.BinaryRowOffsets.empty() ||
        row >= system.BinaryRowOffsets.size() - 1u)
    {
        throw std::out_of_range("WH2 test precode row is out of range");
    }
    const uint32_t first = system.BinaryRowOffsets[row];
    const uint32_t last = system.BinaryRowOffsets[row + 1u];
    if (first > last || last > system.BinaryRowColumns.size()) {
        throw std::out_of_range("WH2 test precode offsets are malformed");
    }
    const uint32_t removed = last - first;
    system.BinaryRowColumns.erase(
        system.BinaryRowColumns.begin() + first,
        system.BinaryRowColumns.begin() + last);
    for (size_t offset = row + 1u;
         offset < system.BinaryRowOffsets.size();
         ++offset)
    {
        system.BinaryRowOffsets[offset] -= removed;
    }
    system.BinaryRowOffsets.erase(
        system.BinaryRowOffsets.begin() + row + 1u);
}

inline void AppendPrecodeRowForTesting(
    PrecodeSystem& system,
    const std::vector<uint32_t>& row)
{
    if (system.BinaryRowOffsets.empty()) {
        system.BinaryRowOffsets.push_back(0u);
    }
    const size_t max_size = std::numeric_limits<uint32_t>::max();
    if (system.BinaryRowColumns.size() > max_size ||
        row.size() > max_size - system.BinaryRowColumns.size())
    {
        throw std::length_error("WH2 test precode graph is too large");
    }
    system.BinaryRowColumns.insert(
        system.BinaryRowColumns.end(), row.begin(), row.end());
    system.BinaryRowOffsets.push_back(
        (uint32_t)system.BinaryRowColumns.size());
}

enum class PrecodeBuildAllocationFailurePoint : uint8_t
{
    None,
    FinalValidation
};

enum class PrecodeBuildAllocationFailureException : uint8_t
{
    BadAlloc,
    LengthError
};

void SetPrecodeBuildAllocationFailurePointForTesting(
    PrecodeBuildAllocationFailurePoint point,
    PrecodeBuildAllocationFailureException exception);
uint32_t PrecodeBuildAllocationFailureHitsForTesting();
void TriggerPrecodeBuildAllocationFailureForTesting(
    PrecodeBuildAllocationFailurePoint point);

} // namespace test
#endif

/**
    Allocation-free validation of Params and the complete flat CSR shape.
    This is the exact prefix checked by ValidatePrecodeSystem before row views
    or source-hit scratch are used.
*/
bool HasValidPrecodeSystemShape(const PrecodeSystem& system);

/**
    Build the staircase + Shuffle-2 constraint structure.

    Returns false for invalid parameters (BlockCount outside [2, 64000],
    Staircase == 0, SourceHits outside [1, 8], DenseRows > 64,
    HeavyRows > 128, or a full symbol domain that does not fit uint16) or if
    the generated structure fails ValidatePrecodeSystem().  This low-level
    construction primitive may throw std::bad_alloc or std::length_error;
    status-bearing API boundaries translate those exceptions to Wirehair_OOM.
*/
bool BuildPrecodeSystem(const PrecodeParams& params, PrecodeSystem& out);

/**
    Validate every structural invariant consumed by the encoder.  Validation
    uses widened arithmetic and performs no writes to block data.  Validation
    allocates bounded source-hit scratch and may throw std::bad_alloc or
    std::length_error; bool and WirehairResult API boundaries that invoke it
    define their own exception translation.
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

/// Coefficient dispatch for actual encoder/decoder equations.  Alternate
/// families are confined to explicitly constructed experiment systems.
uint8_t HeavyCoefficientForParams(
    const PrecodeParams& params,
    uint32_t heavy_row,
    uint32_t ge_column);

} // namespace wirehair_v2
