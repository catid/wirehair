#pragma once

#include <stdint.h>

#include <vector>

/*
    Certified precode construction for the V2 codec (wirehair-axd).

    Implements the Phase-B-certified `ldpcdense_s<S>_d12_s2_h12` structure
    from experiments/precode (see experiments/precode/README.md and
    results/SUMMARY.md "Ship Gate Decision"):

    - S = GetDenseCount(K) LDPC-staircase parity columns.  Each SOURCE column
      connects to N1 distinct staircase parities; parity row j additionally
      references its own parity column K+j and the staircase link K+j-1.
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

struct PrecodeParams
{
    uint32_t BlockCount;   ///< K: source blocks
    uint32_t Staircase;    ///< S: staircase parity columns
    uint32_t DenseRows;    ///< D2: Shuffle-2 dense binary rows
    uint32_t HeavyRows;    ///< H: Cauchy heavy rows
    uint32_t SourceHits;   ///< N1: staircase parities per source column
    uint64_t Seed;         ///< constraint-generation seed
};

/// Certified rule: S = GetDenseCount(K), D2 = 12, H = 12, N1 = 2
PrecodeParams MakeCertifiedParams(uint32_t block_count, uint64_t seed);

struct PrecodeSystem
{
    PrecodeParams Params;

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
        [0, K + S + D2), stored as a sorted column list.  Row 0 has exactly
        ceil(span/2) columns; every row r > 0 differs from row r - 1 in
        exactly two columns.

        Known limitation (inherited from the certified reference
        construction): at tiny EVEN spans (K=2 and K=4 with the certified
        table) a later row's weight can walk down to exactly zero — a dead
        constraint — at roughly 1e-4 systems.  Guarding would break the
        exact-2-difference invariant, and the production policy does not
        rely on the precode at such tiny K, so it is documented rather
        than patched.
    */
    std::vector<std::vector<uint32_t>> DenseRowColumns;
};

/**
    Build the staircase + Shuffle-2 constraint structure.

    Returns false only for invalid parameters (BlockCount outside
    [2, 64000], Staircase == 0, SourceHits outside [1, 8], or a span that
    does not fit the 16-bit deck domain).
*/
bool BuildPrecodeSystem(const PrecodeParams& params, PrecodeSystem& out);

/**
    GF(256) coefficient of heavy row r at GE column c.

    Cauchy element 1 / (X ^ Y) with Y = r and X = H + (c mod (256 - H)):
    nonzero everywhere, and every square submatrix drawn from one
    244-column window is invertible.  heavy_rows must be <= 128.
*/
uint8_t HeavyCoefficient(
    uint32_t heavy_row,
    uint32_t ge_column,
    uint32_t heavy_rows);

} // namespace wirehair_v2
