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
same `def` distribution.  Heavy GF(256) costs (muladd vs XOR) are not folded
into the XOR cost columns; interpret `H`-heavy schemes with that in mind.

`N` in other harnesses is `K` here; payload bytes scale as `K * block_bytes`.
All cost columns are counted in block units (block-size independent); convert
with `experiments/peeling/xor_bench` timings.

## Correctness

The self-test proves the projection replay correct via rank invariance: for
thousands of random instances across all scheme kinds, `peeled + residual_rank`
must equal the brute-force GE rank of the entire binary system, because peeling
is just a solving strategy and cannot change rank.

```bash
bash experiments/precode/build.sh
./experiments/precode/precode_sim --self-test
```

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
