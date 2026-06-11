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
thousands of random instances across all scheme kinds (including the `_s2`
and `_n1<X>` variants), `peeled + residual_rank` must equal the brute-force GE
rank of the entire binary system, because peeling is just a solving strategy
and cannot change rank.  It additionally checks that the `--ge-replay`
explicit elimination (both row orders) reproduces the rank accumulator's rank
on every one of those instances, and that `--paired` received-row generation
is scheme-independent (see below).  Further clauses:

- `_s2` structure: the first D2 row has exactly `ceil(span/2)` columns and
  every subsequent D2 row differs from its predecessor in EXACTLY 2 columns;
- `_n1<X>`: every source column hits exactly X distinct staircase parities
  (and tokens without `_n1` keep the historical X = 3);
- degree-law chi-square: for every non-Wirehair `--rowdist`, ~200k sampled
  degrees at K = 3200 are tested against an independently transcribed copy of
  the analytic weight law (bins pooled to expected count >= 10; bound
  `df + 8*sqrt(2*df) + 30` — loose against noise, far below the chi-square a
  mis-ported spike/tau/cap produces);
- `--max-inact`: capping below a trial's natural inactivation count must
  flag it as a runaway failure, capping at exactly the natural count must
  reproduce the unlimited result bit-for-bit;
- determinism of the new kinds (Shuffle-2 D2 rows, n1 variants, robust
  soliton rowdist).

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
every remaining row holding that bit.  `--ge-pivot-window N` implies
`--ge-replay` and changes only the replay pivot choice: it scans the next `N`
remaining rows for the current column and chooses the candidate with the lowest
remaining tail popcount, falling back to the first later candidate if the
bounded window has none.  Recorded per trial:

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
| `ldpcdense_s<S>_d<D2>` | S staircase parity columns + D2 iid p=0.5 dense rows over all `K + S + D2` binary columns |
| `ldpcdense_n1<X>_s<S>_d<D2>` | same, but each SOURCE column connects to X (2, 3 or 4) distinct staircase parities; `_n1` absent = 3 (token and RNG-stream compatible with the historical behavior) |
| `ldpcdense...`+`_s2` | D2 rows generated Shuffle-2-style instead of iid (exact rule below); combines with `_n1<X>`, e.g. `ldpcdense_n12_s50_d12_s2` |
| `heavyonly` | no binary precode, heavy only (default H=16) |
| any token + `_h<H>` | override heavy row count (after `_s2` when both present: `..._s2_h<H>`) |

### `_s2`: Shuffle-2 structured D2 rows (exact generation rule)

This is the production-implementable form of the D2 dense rows, mirroring the
window mechanics of `MultiplyDenseRows` (WirehairCodec.cpp:786) /
`ShuffleDeck16` (WirehairTools.cpp:398).  The rule below is what the simulator
certifies and what a codec implementation must reproduce bit-exactly over the
matrix bits (`span = K + S + D2`, the full binary column space the D2 rows
cover: source + staircase parity + dense columns):

1. **Deck construction.** `deck` = permutation of the `span` column ids,
   produced by the Sattolo-style inside-out shuffle used by `ShuffleDeck16`:
   `deck[0] = 0; for ii in [1, span): jj = rand() % ii; deck[ii] = deck[jj];
   deck[jj] = ii;`.  The simulator draws `jj` as an unbiased full-width
   uniform in `[0, ii)` from the trial's constraint-row splitmix64 RNG
   (immediately after the staircase rows are generated, so the deck stream is
   seeded by the scheme/trial seed).  Production will instead consume a
   `PCGRandom` seeded with the dense seed in 8-/16-bit chunks exactly as
   `ShuffleDeck16` does — the *structure* of the rule is what the simulator
   certifies, not the raw bit source.
2. **First-row window.** `set_count = ceil(span/2) = (span+1) >> 1`.  The
   first D2 row has exactly the columns `deck[0 .. set_count)` set.
3. **Per-row flip selection.** Reshuffle the deck (same shuffle, continuing
   the same RNG stream).  Interpret `set_half = deck[0 .. set_count)` and
   `clear_half = deck[set_count .. span)`.  The next `floor(D2/2)` rows are
   each `row[i+1] = row[i] XOR {set_half[ii], clear_half[ii]}` for
   `ii = 0, 1, 2, ...` — exactly two bit flips per row; the two flipped
   columns are always distinct because deck entries at distinct positions are
   distinct.  Then reshuffle the deck once more and generate the remaining
   `floor(D2/2) - 1 + (D2 & 1)` rows by the same flip rule with `ii`
   restarting at 0.  Row count check: `1 + floor(D2/2) + floor(D2/2) - 1 +
   (D2 & 1) = D2`, mirroring production's two half-loops around the second
   reshuffle.

Documented deviations from production `MultiplyDenseRows` (both deliberate):
(a) ONE window spanning all `span` columns instead of production's successive
blocks of `dense_count` columns — production XOR-accumulates every block's
pattern into the destination rows, which destroys the global 2-flip chain;
with a single window the "each row = previous row + 2 flips" property holds
over the whole row, which is the structure E5 is testing; (b) no `rows`
destination-permutation deck — rows are emitted in generation order, which
cannot change the linear system (a codec implementation MAY reintroduce the
rows deck for stream parity with `MultiplyDenseRows`).

Value-generation cost charged to `precode_gen_xors_mu`: `set_count +
2*(D2-1)` block XORs (first-row accumulation plus two incremental flips per
subsequent row), vs the `2.5 * span` window estimate charged to the iid
model.

## `--rowdist`: received-row degree distributions

| token | law |
| --- | --- |
| `wirehair` (default) | production `GeneratePeelRowWeight` |
| `lt_m1_c64` | truncated LT (soliton-like, `w(1) = 1/K`, `w(d) = 1/(d(d-1))`), min degree 1, cap 64, renormalized |
| `lt_m1_c16` | same LT family, min degree 1, cap 16 |
| `lt_m2_c1024` | LT family, min degree 2, cap 1024 |
| `rs_c001_d50_c128` | robust soliton `c = 0.01, delta = 0.50`, min degree 1, cap 128 |

`rs_c001_d50_c128` mirrors `peel_sweep`'s `robust_soliton_weight` exactly:
`R = max(1, c*ln(K/delta)*sqrt(K))`, `spike = clamp(floor(K/R), 1, K)`,
`w(d) = lt(d) + tau(d)` with `tau(d) = R/(d*K)` for `d < spike`,
`R*ln(R/delta)/K` at `d == spike`, 0 above; the cap truncates the law (at
K = 3200 the spike sits at d = 645, above the 128 cap, so only the `R/(d*K)`
shoulder survives) and the cumulative sum renormalizes.  All non-Wirehair
laws are chi-square-validated against an independent transcription in the
self-test.

## `--max-inact <count>`: runaway guard

Retuned low-cap rowdists (especially `lt_m1_c16` at large K) can drive the
peeler into pathological inactivation counts, making the residual-rank step
quadratically expensive.  `--max-inact N` aborts a trial as soon as the
inactivated column count EXCEEDS N (0 = unlimited, the default).  Aborted
trials are recorded as *runaways*: they count as decode failures in
`fail_rate` / `fail_rate_noheavy`, are reported in the `runaway_rate` column
(appended at the end of the CSV), and are EXCLUDED from every def/inact/cost
mean (those statistics do not exist for an abandoned solve, and a partial
inactivation count would bias the means of the surviving trials).  In
`def_pdf` runaways appear as a sentinel `999999:<rate>` bucket so the pdf
stays normalized over all trials and re-scoring it at any H still reproduces
the fail rate exactly (`rank_total.py` asserts this).  A nonzero
`runaway_rate` is signal, not an error — it certifies that the
(rowdist, scheme, K) cell blows past the bound.

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
  re-evaluate any alternative heavy-row count after the fact.  Normalized
  over ALL trials; `--max-inact` runaways show up as a `999999:<rate>`
  sentinel bucket (fails at every H), keeping `P(def > H) == fail_rate`.
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
- `runaway_rate` (always present, last column): fraction of trials aborted by
  `--max-inact` (0 when the guard is off).  Runaway trials are inside
  `fail_rate`/`fail_rate_noheavy` but outside every mean column, whose
  denominators are the completed trials only.

Existing column order is unchanged; new columns are appended at the end so
old CSV parsers keep working.
