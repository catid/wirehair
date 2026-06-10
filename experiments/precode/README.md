# Precode / Dense-Matrix-Replacement Experiments

This directory holds the standalone precode simulator used to evaluate
replacements for the production dense binary rows (+ heavy rows).  Unlike
`experiments/peeling/peel_sweep`, which models only the peel graph and reports
residual sizes, `precode_sim` models the FULL intermediate-block linear system
and reports real decode-failure rates:

- `L = K` source columns `+ D` binary precode columns `+ H` heavy columns.
- Constraint rows: `D` binary precode rows (scheme-dependent) and `H` heavy
  GF(256) rows.
- Received rows: `K + OH` packets, each with a peel degree over the `K` source
  columns (production `GeneratePeelRowWeight` by default) plus `--mix` (default
  3, production-like) distinct columns over the `D + H` precode columns.

After peeling, the simulator replays peeled-column substitution in solve order
to compute the exact GF(2) projection of every unused row onto the inactivated
column set, then computes the binary residual rank.  Decode success requires
`def = inactivated - rank <= H`.

Heavy rows are modeled as a full-coverage MDS rank patch (explicit Cauchy rows
over all GE columns).  The production heavy matrix only covers the last 18 GE
columns, so production failure rates are an upper bound of the model's for the
same `def` distribution.  Heavy GF(256) costs are reported separately as
block-unit muladd counts (`heavy_muladds_mu = H*inact + H^2`, `heavy_divs_mu =
H`); `rank_total.py` folds them into a single XOR-equivalent total using the
measured muladd:xor ratio per block size.

`N` in other harnesses is `K` here; payload bytes scale as `K * block_bytes`.
All cost columns are counted in block units (block-size independent); convert
with `experiments/peeling/xor_bench` timings.

## Correctness

The self-test proves the projection replay correct via rank invariance: for
thousands of random instances across all scheme kinds, `peeled + residual_rank`
must equal the brute-force GE rank of the entire binary system, because peeling
is just a solving strategy and cannot change rank.  It additionally checks that
the `--ge-replay` explicit elimination (both row orders) reproduces the rank
accumulator's rank on every one of those instances, and that `--paired`
received-row generation is scheme-independent (see below).

```bash
bash experiments/precode/build.sh
./experiments/precode/precode_sim --self-test
python3 experiments/precode/rank_total.py --self-test
```

## `--ge-replay`: real GE elimination replay

`--ge-replay` builds the residual binary system explicitly after the rank
analysis: rows are the unused-row projections onto the inactivated columns
(the exact bitsets the rank accumulator consumes), columns ordered in
INACTIVATION order.  Rows are sorted once, up front, by ascending initial
popcount (cheapest-pivot-first, a production-like heuristic;
`--ge-replay-reverse` flips to descending order as a sensitivity check, and
implies `--ge-replay`).  Plain GE then runs column by column: the first
remaining row with the column bit set becomes the pivot and is XORed into
every remaining row holding that bit.  Recorded per trial:

- word XORs: each elimination costs `words - col/64` 64-bit word XORs, the
  words a triangular implementation actually touches (`ge_real_bitops_mu`);
- row eliminations, i.e. real GE block XORs (`ge_real_rowops_mu`);
- fill-in: the popcount growth of each pivot row above its initial weight,
  measured at the moment it is selected as pivot (`fill_in_mu`);
- deficient (pivotless) column positions, counted from the END of the
  inactivation order (last inactivated column = position 0).

Per trial the replay rank is asserted equal to the rank accumulator's rank.
The deficiency-position census feeds the heavy-band question (production's
heavy matrix covers only the last 18 GE columns):

- `def_outside_w18_rate`: fraction of trials with at least one deficient
  column at from-the-end position >= 18 (i.e. outside a width-18 band);
- `def_band_w95` / `def_band_w99`: 95th/99th percentile (nearest-rank) over
  trials of the band width needed to cover every deficiency, i.e. deepest
  deficient from-the-end position + 1; trials with no deficiency count as 0.

## `--paired`: common-random-number trials across schemes

`--paired` generates received rows from a second RNG whose per-trial seed
omits the scheme token (constraint rows keep the scheme-dependent stream), so
cross-scheme fail-rate deltas at the same `(K, oh, trial, --seed)` are
CRN-paired.  Exactly what is and is not paired:

- PAIRED for any two schemes at the same `(K, oh, trial, --seed, --rowdist,
  --mix)`: the per-row peel degree sequence and the source-column sets
  (columns `< K`), including their draw order.
- PAIRED only when the schemes have equal precode width `D + H`: the mix
  columns.  Each row draws exactly `mix` raw 64-bit values from the paired
  stream (consumed even when `D + H = 0`, keeping streams in lockstep) and
  maps them into the scheme-local `D + H` precode space by modulo plus
  deterministic linear probing for distinctness.  Equal widths give identical
  mix columns; different widths necessarily give different mix columns, but
  the RNG stream itself never diverges.
- NOT PAIRED, by design: the binary precode constraint rows (LDPC parity
  assignments, dense row bits), which remain a function of the scheme token.

## `rank_total.py`: heavy-aware total-cost ranking

Folds the per-scheme XOR ledger and the GF(256) heavy ledger into one
calibrated `total_block_ops(bb)` per block size and prints a ranked table per
`(K, oh, bb)`, re-scoring failure for alternative heavy counts
`H' in {6,8,12,16}` from `def_pdf` (asserting per row that re-scoring at the
modeled H reproduces `fail_rate`), plus the muladd:xor crossover ratio at
which the ldpcdense H=12 scheme stops beating the dense H=6 baseline.
`--pessimistic` doubles the heavy term.  Old CSVs without the heavy/replay
columns are accepted (heavy is derived as `H*inact_mu + H^2`, GE falls back to
the `R^2/2` estimate).

```bash
python3 experiments/precode/rank_total.py results/hybrid_K*.csv \
    --xor-csv experiments/peeling/results/cold_xor_calibration.csv \
    --muladd-csv experiments/peeling/results/muladd_calibration.csv \
    --bb 1280,102400,1048576
```

Missing/empty calibration files are skipped gracefully (`--assume-ratio`,
default 4.0, fills in).

## Schemes

| token | meaning |
| --- | --- |
| `none` | no precode (pure LT control) |
| `dense` | `D = GetDenseCount(K)` iid p=0.5 binary rows (production idealization), H=6 |
| `dense_d<D>` | explicit dense row count |
| `densesparse_w<W>[_d<D>]` | fixed-weight-`W` random binary rows |
| `ldpc[_s<S>]` | LDPC-staircase: each source column in 3 random parities, double-diagonal parity chain |
| `ldpc2x` | staircase with `S = 2 * GetDenseCount(K)` |
| `ldpctri[_s<S>]` | circulant triple `{a, a+b, a+2b} mod S` variant |
| `heavyonly` | no binary precode, heavy only (default H=16) |
| any token + `_h<H>` | override heavy row count |

## Example

```bash
./experiments/precode/precode_sim --K 3200 --schemes dense,ldpc,dense_d31_h12 \
    --oh 0,1,2,5 --trials 4000 --threads 120
```

Sweep driver used for `results/sweep_K*.csv`: see git history (`/tmp/precode_sweep.sh`
pattern): K in {1000, 3200, 10000, 32000, 64000}, 14 schemes, OH in {0,1,2,5}.

## Output columns

- `fail_rate`: P(def > H) — decode failure under the MDS heavy-patch model.
- `fail_rate_noheavy`: P(def > 0). With H heavy columns and OH < H this is 1 by
  construction (binary rows minus unknowns = OH - H), so it is only meaningful
  for `none`.
- `def_pdf`: empirical def distribution (`value:probability`), enough to
  re-evaluate any alternative heavy-row count after the fact.
- `inact_*`: inactivated column count (GE width).
- `backsub_xors_mu`: sum of popcount(projection) over peeled columns — the
  block-XOR cost of pushing inactivated solutions back into peeled columns;
  this is the dominant block-op term and scales with `inact`.
- `sparse_solve_xors_mu`: peel-side block-XOR proxy (degree-1 per used row).
- `precode_gen_xors_mu`: encoder-side precode value generation block XORs
  (dense uses the Shuffle-2 incremental estimate `2.5 * (K + D)`).
- `ge_block_xors_mu`, `ge_bitops_mu`: `R^2/2` and `rows * R^2 / 64` estimates.
- `heavy_muladds_mu`: mean `H*inact + H^2` block-unit GF(256) muladds (heavy
  value substitution + heavy GE elimination); `heavy_divs_mu`: `H` block-unit
  GF(256) divides.
- With `--ge-replay` only (columns appended after the above):
  `ge_real_bitops_mu` (measured 64-bit word XORs), `ge_real_rowops_mu`
  (measured row eliminations = real GE block XORs), `fill_in_mu`,
  `def_outside_w18_rate`, `def_band_w95`, `def_band_w99` (see above).

Existing column order is unchanged; new columns are appended at the end so
old CSV parsers keep working.
