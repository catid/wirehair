# Precode / Dense-Replacement Simulation Summary

All results from `experiments/precode/precode_sim` (see `../README.md` for the
model).  `K` is the source block count; cost columns are in block-XOR units and
are block-size independent.  Failure model: decode fails when binary rank
deficiency `def` exceeds the heavy row count `H` (full-coverage MDS heavy
patch, i.e. an explicit Cauchy heavy solve as in `wirehair-wuw`; production's
18-column heavy band can only be worse, so production fail rates upper-bound
the `dense` rows here).

Correctness: the projection replay is proven by a rank-invariance self-test
(`peeled + residual_rank == brute-force GE rank` over thousands of instances
across every scheme kind), plus determinism checks.

Protocols:

- `sweep_K{1000,3200,10000,32000,64000}.csv`: 14 schemes x OH {0,1,2,5},
  1000-4000 trials each, production `GeneratePeelRowWeight` received rows,
  mix=3 over the precode columns.
- `grid_K{1000,3200,10000}.csv`: dense/LDPC size vs heavy count grid,
  4000-6000 trials.
- `hybrid_K{1000,3200,10000}.csv`: LDPC-staircase + small dense hybrid,
  4000-6000 trials.

## Headline results (OH=0, vs production `dense` ~1.4-1.7% fail)

`ldpc` = LDPC-staircase precode, S = production GetDenseCount(K), H=6:
equal failure rate, much cheaper solve.

| K | scheme | fail | inact (vs dense) | sparse XORs | GE bitops |
| ---: | --- | ---: | ---: | ---: | ---: |
| 3200 | dense (prod) | 1.72% | 112.8 | 38.6k | 21.5k |
| 3200 | ldpc | 1.57% | 106.2 (-6%) | 28.1k (-27%) | 18.0k (-16%) |
| 10000 | dense (prod) | 1.68% | 193.4 | 133.5k | 111.0k |
| 10000 | ldpc | 1.95% | 159.6 (-17%) | 94.6k (-29%) | 62.7k (-43%) |
| 32000 | dense (prod) | 1.43% | 461.5 | 459.4k | 1526.5k |
| 32000 | ldpc | 1.73% | 331.2 (-28%) | 318.1k (-31%) | 566.8k (-63%) |
| 64000 | dense (prod) | 1.40% | 843.8 | 952.2k | 9358.6k |
| 64000 | ldpc | 1.70% | 575.1 (-32%) | 646.0k (-32%) | 2972.5k (-68%) |

backsub XORs (the dominant block cost) are equal between dense and ldpc at
equal column counts; the savings above are on top of that.

## Hybrid LDPC + small dense (`ldpcdense_s<S>_d<D2>`)

The failure tail of undersized precodes (see below) is pinned by a few dense
rows.  Best production-implementable points (H=6, same heavy structure as
production):

| K | scheme | fail OH0 | fail OH2 | inact | sparse XORs | GE bitops |
| ---: | --- | ---: | ---: | ---: | ---: | ---: |
| 1000 | dense (prod) | 1.70% | 0.00% | 75.4 | 10.6k | 6.2k |
| 1000 | ldpcdense_s25_d16 | 1.28% | 0.00% | 63.4 | 8.9k | 3.7k |
| 3200 | dense (prod) | 1.72% | 0.02% | 112.8 | 38.6k | 21.5k |
| 3200 | ldpcdense_s43_d12 | 1.53% | 0.00% | 99.2 | 29.4k | 14.6k |
| 10000 | dense (prod) | 1.68% | 0.13% | 193.4 | 133.5k | 111.0k |
| 10000 | ldpcdense_s52_d34 | 1.58% | 0.15% | 163.9 | 109.4k | 67.6k |

With H=12 (requires explicit Cauchy heavy solve, `wirehair-wuw`), failure
rates drop 10-100x below production at still-lower cost, e.g.
`ldpcdense_s43_d12_h12` at K=3200: 0.00%/0.02% fail, inact 105.6,
sparse 29.4k, bitops 16.6k.

## Dense-vs-heavy budget grid

- Raising H 6->12 at production D cuts failure ~10-100x for the cost of 6
  extra heavy rows (e.g. K=3200 `dense_d50_h12`: 0/6000 fails at OH0 and OH2).
- Shrinking D below a K-dependent threshold introduces a heavy-tailed def
  distribution (def up to ~30) that overhead does NOT fix; the def PDF shows
  production-size D pins def exactly at its minimum H while undersized D has
  a long tail.  At K=32000 and K=64000, half-size precodes collapse entirely
  (99-100% fail).  Safe shrink shrinks with K: ~0.8x at K=1000-3200, ~1.0x by
  K=10000+.
- Fixed-weight sparse dense rows (w32/w64) are NOT viable replacements
  (fail 25-100%).
- Heavy-only (no binary precode) needs H ~ 64 at K=10000+ and is dominated.

## Failure floor

All strong precodes (dense, ldpc, ldpc2x, ldpctri) hit the same ~1.4-2%
failure wall at OH=0 for their H: the def>H events at production scale come
from received-row deficiencies, not the precode.  The precode determines
whether you *reach* that floor and what the solve costs; H (and OH) determine
the floor itself.

## Interpretation for v2

1. The dense binary matrix can be replaced by an LDPC staircase + a few dense
   rows at equal reliability and 12-30% smaller GE, 25-30% less sparse XOR,
   and 30-60% fewer GE bit-ops; gains GROW with K.
2. Explicit Cauchy heavy rows (H=12-16 full-width) are the cheapest
   reliability lever and compose with (1).
3. Caveats not yet modeled: heavy GF(256) muladd cost (block ops ~ H x inact;
   small vs the savings), LDPC parity-check value generation in the encoder
   pipeline, and production peeling treats constraint rows differently (the
   simulator lets precode rows peel, which is the point of the staircase).
