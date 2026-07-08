#pragma once

#include "WirehairV2Peel.h"
#include "WirehairV2Precode.h"

#include <stdint.h>
#include <vector>

/*
    Encoder value phase for the certified V2 precode (wirehair-axd Phase 3).

    Given the K source block values, computes the S + D2 + H intermediate
    parity block values so that every constraint row of the PrecodeSystem
    sums to zero over the full column value vector
    [source | staircase | dense | heavy]:

    - STAIRCASE rows solve by one forward pass: parity j is the XOR of its
      row's other columns (source hits + the j-1 link), exactly
      min(N1,S)*K + S - 1 block ops.
    - Shuffle-2 DENSE rows couple all D2 dense parity values, because every
      row spans source, staircase AND dense columns.  The rows are
      transformed to consecutive-row differences (an invertible unit
      bidiagonal row transform, so solvability and the solution are
      unchanged): the difference of two consecutive rows is exactly the two
      flipped columns, so each transformed row's known part costs at most 2
      block ops -- the incremental structure behind the certified
      ceil(span/2) + 2*(D2-1) cost model.  What remains is a D2 x D2 GF(2)
      system over the dense columns, solved by Gauss-Jordan with block RHS.
      That corner is NOT guaranteed invertible -- when singular the
      constraints have no (unique) solution for these source values and the
      function returns false.  Use DenseCornerInvertible to measure or gate
      on this without touching block data.
    - Cauchy HEAVY rows: the known part is a GF(256) muladd accumulation
      over the K + S + D2 known columns; the H x H corner over the heavy
      columns (consecutive GE columns, always distinct within one
      256 - H = 244 window) is Cauchy, hence invertible, and is solved by
      GF(256) Gauss-Jordan over blocks.

    Cost accounting (PrecodeEncodeStats) counts every block-sized memory
    operation: a first-term copy into an accumulator counts like an XOR,
    matching the Phase B `precode_gen_xors_mu` convention; zero-filling an
    empty accumulator counts as nothing.
*/

namespace wirehair_v2 {

struct PrecodeEncodeStats
{
    /// Staircase forward pass block ops
    /// (certified model: min(N1,S)*K + S - 1)
    uint64_t StaircaseBlockOps;

    /// Dense difference-row known-part block ops
    /// (certified model: ceil(span/2) + 2*(D2 - 1), counted even when the
    /// dense corner turns out singular)
    uint64_t DenseKnownBlockOps;

    /// GF(2) Gauss-Jordan block XORs + solution copies for the dense corner
    uint64_t DenseSolveBlockOps;

    /// Plain block XORs folding same-residue known columns into the
    /// mod-(256 - H) residue buckets before the heavy accumulation
    /// (model: K + S + D2 when the bucketed path is taken, else 0)
    uint64_t HeavyBucketXors;

    /// GF(256) muladd-class block ops accumulating the heavy known parts
    /// (model: H * min(K + S + D2, 256 - H) bucketed, H * (K + S + D2)
    /// direct)
    uint64_t HeavyMulAdds;

    /// GF(256) block ops (muladd/div/copy) solving the H x H Cauchy corner
    uint64_t HeavySolveBlockOps;
};

/**
    Structure-only feasibility check: GF(2) rank of the D2 x D2 submatrix of
    the dense rows over the dense columns [K + S, K + S + D2).  When this is
    rank-deficient, ComputePrecodeValues must fail for generic source data
    (the staircase part is always solvable and the heavy corner is always
    invertible, so this corner is the ONLY feasibility gate).  Requires
    DenseRows <= 64.  D2 == 0 counts as invertible.
*/
bool DenseCornerInvertible(const PrecodeSystem& system);

/**
    Compute the S + D2 + H intermediate parity block values.

    source_blocks: K contiguous blocks of block_bytes each.
    parity_blocks: out, S + D2 + H contiguous blocks -- staircase parities
    [0, S), dense parities [S, S + D2), heavy parities [S + D2, S + D2 + H),
    i.e. parity i is global column K + i.

    Returns false on invalid arguments (block_bytes 0 or > 2^31 - 1,
    DenseRows > 64, HeavyRows > 128, row vectors inconsistent with Params)
    or when the D2 x D2 dense corner is singular.  On the singular-corner
    path stats (if given) still hold the staircase and dense known-part
    counts, so cost measurements do not require a feasible seed.
*/
bool ComputePrecodeValues(
    const PrecodeSystem& system,
    const uint8_t* source_blocks,
    uint32_t block_bytes,
    uint8_t* parity_blocks,
    PrecodeEncodeStats* stats = nullptr);

/**
    Compute one V2 recovery symbol from a row over the intermediate domain
    [source | staircase | dense | heavy].

    `row_columns` uses the same global column numbering as
    GenerateRecoveryMatrixRows(): columns < K reference source_blocks, and
    columns >= K reference parity_blocks at offset column - K.  The output is
    the GF(2) XOR of all referenced block values.  A first-term copy counts as
    one block op, matching the precode cost convention; an empty row produces
    the zero block and zero ops.

    `block_out` must not overlap the source/parity block arrays.
*/
bool ComputeRecoveryBlock(
    const PrecodeSystem& system,
    const uint8_t* source_blocks,
    const uint8_t* parity_blocks,
    uint32_t block_bytes,
    const std::vector<uint32_t>& row_columns,
    uint8_t* block_out,
    uint64_t* block_ops_out = nullptr);

/**
    Compute one encoder output block.

    block_id < K copies the corresponding source block.  block_id >= K treats
    block_id - K as the zero-based recovery row index, generates that row with
    GenerateRecoveryMatrixRow(), then evaluates it with ComputeRecoveryBlock().
    `row_seed` is the seed for the recovery row stream, and `mix_count` is the
    number of precode columns appended to each recovery row.

    `parity_blocks` may be null only for source block ids; recovery ids require
    the S + D2 + H parity blocks already computed by ComputePrecodeValues().
    `block_out` must not overlap the source/parity block arrays.
*/
bool ComputeEncodedBlock(
    const PrecodeSystem& system,
    const PeelingCodec& codec,
    uint64_t row_seed,
    uint32_t mix_count,
    const uint8_t* source_blocks,
    const uint8_t* parity_blocks,
    uint32_t block_bytes,
    uint32_t block_id,
    uint8_t* block_out,
    uint64_t* block_ops_out = nullptr);

} // namespace wirehair_v2
