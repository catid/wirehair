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

## Pivot-window GE replay (`wirehair-ht0`)

`ht0_ge_pivot_window.csv` adds a replay-only `--ge-pivot-window` screen over
K={1000,3200,10000}, schemes={dense,ldpc}, OH=0, 300 paired trials.  Windows
16/64/256 choose the candidate pivot with the lowest remaining tail popcount
within the next N replay rows.

Readout: pivot windows cut the fill-in counter by 80-98%, but they do not cut
the actual replayed GE work.  `ge_real_bitops_mu` and `ge_real_rowops_mu` stay
flat or drift slightly upward (roughly -0.1% to +0.3% across this screen).  This
does not justify a production pivot-selection change; keep production pivoting
simple unless a later RHS/value-schedule metric shows a real block-XOR win.

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

## June 2026 Big-K Session

New artifacts:

- `e3a_K16000.csv`, `e3a_K20000.csv`, `e3a_K32000.csv`,
  `e3a_K48000.csv`, `e3a_K64000.csv`: paired big-K hybrid grids.
- `e1_hgrid_K32000.csv`, `e1_hgrid_K64000.csv`: GE replay and heavy-band
  census for the chosen frontier.
- `e5_K32000.csv`, `e5_K64000.csv`: shuffle2 and n1/n12 implementability
  checks at the full-size D2=12 rule point.
- `phaseb_cert_20260626.csv`: paired Phase-B certification run with seed
  `0x5eed0002`, OH0/OH1, and 200k/100k/50k trials by size band.

### Rule

Use the existing production dense-count table for S, add D2=12 extra dense
binary rows, and use H=12 full-width heavy rows for the conservative prototype
and Phase-B certification candidate:

`ldpcdense_s<GetDenseCount(K)>_d12_h12`

This rule is intentionally conservative.  The H=6 `ldpcdense_sDprod_d12`
rows are already strong at K=32000 and K=64000, but H=12 provides the desired
tail suppression without giving up the LDPC sparse-solve savings.

| K | dense H6 fail OH0 | selected H12 rule | fail OH0/OH1/OH2 | sparse XORs vs dense | GE vs dense |
| ---: | ---: | --- | ---: | ---: | ---: |
| 32000 | 2.20% | `ldpcdense_s190_d12_h12` | 0.033% / 0.033% / 0.000% | -29.6% | -54.4% real |
| 48000 | 1.77% | `ldpcdense_s370_d12_h12` | 0.000% / 0.000% / 0.000% | -35.0% | -55.2% estimated |
| 64000 | 1.93% | `ldpcdense_s346_d12_h12` | 0.067% / 0.000% / 0.000% | -31.4% | -64.2% real |

K=48000 keeps the production table value S=370 even though the table is
non-monotone relative to K=64000 S=346.  The measured K=48000 point is strong,
so do not smooth the table before certification.

### Phase B Certification

`phaseb_cert_20260626.csv` certifies the primary
`ldpcdense_s<GetDenseCount(K)>_d12_h12` rule and secondary
`ldpcdense_s<GetDenseCount(K)>_d12_s2_h12` rule against dense with paired
received-row randomness.  Trial counts were 200k for K<=20000, 100k at
K=32000, and 50k at K=48000/64000.

The primary rule clears the failure gate at every K and OH0/OH1; every measured
failure rate is below dense, with the widest primary 95% binomial half-width
0.020pp at K=64000 OH0.  For K>=10000 it also clears the sparse-solve gate by
27.9-34.9%.  The inactivation gate needs the same interpretation as the H=12
rule selection: raw `inact_mu` is only -6.7% at K=10000 because the candidate
has six additional heavy rows, but `rank_mu` (the non-heavy residual width) is
-10.1% there and improves to -30.1% by K=64000.  If the ship gate is interpreted
as literal raw `inact_mu`, K=10000 needs a small retune; under residual-width
interpretation the primary table is certified.

| K | dense fail OH0/OH1 | primary fail OH0/OH1 | primary rank / raw inact / sparse | s2 fail OH0/OH1 | s2 gen vs primary |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 10000 | 1.754%/0.179% | 0.042%/0.018% | -10.1% / -6.7% / -27.9% | 0.085%/0.046% | -36.5% |
| 16000 | 1.654%/0.150% | 0.053%/0.017% | -17.0% / -14.3% / -28.4% | 0.093%/0.046% | -36.4% |
| 20000 | 1.618%/0.101% | 0.049%/0.009% | -19.7% / -17.5% / -28.8% | 0.071%/0.026% | -36.4% |
| 32000 | 1.724%/0.112% | 0.038%/0.002% | -24.9% / -23.3% / -29.6% | 0.079%/0.030% | -36.4% |
| 48000 | 1.610%/0.010% | 0.034%/0.000% | -24.2% / -23.3% / -34.9% | 0.030%/0.000% | -36.5% |
| 64000 | 1.710%/0.186% | 0.052%/0.032% | -30.1% / -29.2% / -31.4% | 0.156%/0.108% | -36.4% |

The `_s2_h12` variant is production-attractive because it cuts D2 generation by
36.4-36.5% and clears the dense-relative failure and sparse-solve gates.  It is
not the primary certified rule: it leaks relative to iid D2 at max K
(K=64000 OH1 is 0.108% vs 0.032%), and its K=10000 residual-width reduction is
just under 10%.

### Cost Model Correction

The old GE-bitop columns overstate absolute GE replay work by about 5x at big
K.  GE replay measured:

- K=32000 dense H6: estimated 1.52M bitops, real 297.8k bitops.
- K=64000 dense H6: estimated 9.42M bitops, real 1.71M bitops.

The relative LDPC and hybrid savings survive replay.  At K=64000, pure LDPC
H6 has -32.3% sparse XORs and -67.7% real GE bitops vs dense.  The conservative
`ldpcdense_s346_d12_h12` point has -31.4% sparse XORs and -64.2% real GE
bitops vs dense.

`rank_total.py` folds receive XORs, precode generation, sparse solve,
backsubstitution, replayed GE row ops where available, and calibrated heavy
GF(256) muladds into a single block-op score.  At 1280-byte blocks and OH0,
pure `ldpc_h12` is about 1.1-1.2% cheaper than dense at K=32000/64000, but it
has no dense binary safety rows.  The conservative full hybrid H12 point is
about +4.9% total block ops at K=32000 and +2.6% at K=64000 vs dense because
of LDPC parity generation, while buying much lower failure tails.

### Heavy Band

The big-K GE replay census gives the production band target:

- H=6 rows: W99 is 13-14, with `def_outside_w18_rate` <= 0.0013.
- H=12 rows: W99 is 17-20, with OH0 `def_outside_w18_rate` around
  0.012-0.021 for the tested dense/LDPC/hybrid rows.

An H=12 production path should not rely on the current 18-column heavy band
being sufficient for every row.  Either support a wider/full heavy solve in the
v2 path or certify a precise band/table regeneration rule.

### E5 Implementability

At K=64000, the baseline H6 point is:

`ldpcdense_s346_d12`: OH0/OH1/OH2 fail = 1.133% / 0.100% / 0.133%.

Shuffle2 alone reduces precode generation work from 353.2k to 224.5k XORs, but
it is not reliability-neutral at max K under H=6:

`ldpcdense_s346_d12_s2`: 1.733% / 0.167% / 0.200%.

Re-scoring the same deficit PDFs at H=12 brings shuffle2 OH0 back to the same
observed 0.067% tail as the baseline H12 rule, so shuffle2 is a certification
candidate only if the implementation commits to H=12.  It should not be the
default H=6-compatible spec.

n12 is not viable at max K:

- `ldpcdense_n12_s346_d12`: 5.033% / 2.533% / 2.400%.
- `ldpcdense_n12_s346_d12_s2`: 35.067% / 26.767% / 24.900%.

The high-D2 control `ldpcdense_s346_d104` gives 1.267% / 0.000% / 0.000%, but
GE grows to about 5.52M estimated bitops vs 3.22M for D2=12, so it is not the
frontier unless certification requires OH1/OH2 zero-failure margins.

Recommendation for Phase B: certify the iid D2=12 H12 rule as the primary
candidate, optionally include shuffle2 H12 as a secondary generation-saving
candidate, and do not include n12 or n12+shuffle2 except as rejected controls.
