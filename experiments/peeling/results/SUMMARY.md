# Peeling Sweep Summary

These experiments measure only the sparse peel graph and inactivation schedule.
They do not model dense solver cost or dense residual success probability.
`N` is the number of blocks; payload bytes and memory pressure scale as
`N * block_bytes`.

## Fresh Random-Only Sweep

Protocol:

- Binary: current `experiments/peeling/peel_sweep`.
- Structures: all 111 random-compatible structures.
- Methods: all 115 peeling/inactivation schedules.
- Sizes: anchor `N=320,3200,32000`, with actual `N` sampled from
  `anchor +/- 10` per trial.
- Trials: 100 per aggregate row.
- Matrix seeds: paired, so every method sees the same generated matrices for a
  fixed structure, N anchor, row count, and trial.
- PDFs: enabled for residual columns and rows.
- Production deterministic Wirehair row IDs are excluded. `wirehair_rand` is
  the fair Wirehair row-weight distribution with fully random columns.

Result directories:

- `fresh_fixedOH_N320_3200_32000_t100`: fixed packet overhead 0, 1, and 2.
- `fresh_pct_N320_3200_32000_t100`: percent overhead 0, 1, 2, and 5.

Validation:

- `bash experiments/peeling/build.sh`
- `./experiments/peeling/peel_sweep --self-test`
- ASan/UBSan self-test build and run.
- Post-patch smoke over representative new structures and methods.
- CSV integrity checks: expected file/row counts, 100 trials, paired-random
  source, all PDFs present, and no production Wirehair structure.
- Duplicate OH0 protocol check after fixing multistart tie-breaking:
  percent-overhead 0 and fixed-overhead 0 match on all non-timing fields.

The bug-fix passes found one issue: multistart schedules used elapsed greedy
time as a tie-breaker after equal residual quality. That could change secondary
component metrics for tied child schedules. The tie-breaker is now deterministic
and structural: residual columns, residual rows, component square sum, max
component, component count, then choice count. Existing fresh-run residual
means/PDFs are unaffected; do not treat the old pre-patch multistart component
fields as deterministic evidence.

## Pure Random-Row Winners

This table excludes `raptorq_ldpc_*` because those variants add LDPC precode
rows on top of the random repair rows and therefore use a different row budget.
`Thr/s` is `1e6 / mean_total_us` from the harness timing; it is useful for
relative comparisons, not as a stable benchmark.

| N | fixed OH | structure | method | mean cols | var | median | p95 | thr/s | cost |
| ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| 320 | 0 | lt_m2_c128 | rqcc_lowref | 11.81 | 5.15 | 12 | 16 | 2827.3 | O(rows*degree + cols) largest degree-2 component + low refs |
| 320 | 1 | d1mix_lt_p2 | rqcc_lowref | 11.24 | 8.01 | 11 | 16 | 3786.4 | O(rows*degree + cols) largest degree-2 component + low refs |
| 320 | 2 | lt_m2_c128 | hyb_rqbeam | 10.70 | 6.60 | 11 | 16 | 1104.9 | O(rqcc_lowref then beam near tail) |
| 3200 | 0 | lt_m2_c256_fold | rqd2_default | 36.42 | 49.84 | 35 | 49 | 164.3 | O(rows*degree + cols) largest degree-2 component/default tie |
| 3200 | 1 | lt_m2_c128_fold | rqcc_lowref | 35.58 | 37.33 | 35 | 47 | 218.8 | O(rows*degree + cols) largest degree-2 component + low refs |
| 3200 | 2 | lt_m2_c256_fold | rqd2_default | 34.45 | 46.24 | 33 | 46 | 157.7 | O(rows*degree + cols) largest degree-2 component/default tie |
| 32000 | 0 | lt_m2_c256_fold | rqd2_default | 123.26 | 481.80 | 123 | 160 | 6.0 | O(rows*degree + cols) largest degree-2 component/default tie |
| 32000 | 1 | lt_m2_c256_fold | rqd2_default | 123.25 | 629.51 | 118 | 168 | 7.9 | O(rows*degree + cols) largest degree-2 component/default tie |
| 32000 | 2 | lt_m2_c256_fold | rqcc_lowref | 125.52 | 365.96 | 125 | 161 | 6.9 | O(rows*degree + cols) largest degree-2 component + low refs |

At large N, `lt_m2_c256_fold` is the strongest pure-random structure found.
The fast RaptorQ-style degree-2 component schedules are effectively tied on
residuals, with `rqd2_default` and `rqcc_lowref` the most practical winners.

## Percent Overhead At N=32000

Pure random-row results, excluding LDPC-precode variants:

| percent OH | structure | method | mean cols | var | median | p95 | thr/s |
| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| 0 | lt_m2_c256_fold | rqcc_lowref | 123.26 | 481.80 | 123 | 160 | 7.2 |
| 1 | lt_m2_c256_fold | rqcc_lowref | 42.20 | 182.52 | 40 | 67 | 12.8 |
| 2 | lt_m2_c256_fold | rqcc_lowref | 31.60 | 52.27 | 32 | 44 | 11.5 |
| 5 | lt_m2_c256_fold | rqcc_lowref | 23.54 | 34.11 | 23 | 35 | 11.7 |

The earlier surprising small residuals at `N=32000` were due to percent
overhead. With zero percent overhead, the pure-random best is about 123
residual columns, close to the expected `sqrt(N)` scale.

## Including RaptorQ LDPC Precode Rows

These are not row-budget-apples-to-apples with the pure random structures:
`raptorq_ldpc_struct` and `raptorq_ldpc_rand` add S LDPC check rows before the
random rows. They are useful to compare structured LDPC checks against a random
same-degree LDPC control.

| N | fixed OH | structure | method | mean cols | var | median | p95 | matrix rows mean |
| ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| 320 | 0 | raptorq_ldpc_rand | raptorq_d2cc | 6.44 | 8.01 | 6 | 11 | 351.01 |
| 320 | 1 | raptorq_ldpc_rand | rqd2_default | 6.95 | 11.70 | 7 | 13 | 352.01 |
| 320 | 2 | raptorq_ldpc_struct | rqcc_lowref | 5.85 | 7.45 | 6 | 11 | 353.01 |
| 3200 | 0 | raptorq_ldpc_struct | raptorq_d2cc | 23.54 | 73.79 | 23 | 40 | 3312.17 |
| 3200 | 1 | raptorq_ldpc_struct | rqd2_default | 24.51 | 72.76 | 24 | 41 | 3313.17 |
| 3200 | 2 | raptorq_ldpc_struct | rqcc_lowref | 23.38 | 78.85 | 23 | 38 | 3314.17 |
| 32000 | 0 | raptorq_ldpc_struct | hyb_rqbeam | 61.74 | 819.68 | 57 | 114 | 32577.50 |
| 32000 | 1 | raptorq_ldpc_struct | rqd2_default | 63.60 | 913.25 | 62 | 122 | 32578.50 |
| 32000 | 2 | raptorq_ldpc_struct | rqcc_lowref | 63.87 | 904.20 | 60 | 117 | 32579.50 |

At medium and large N, the actual structured RaptorQ LDPC rows beat the random
same-degree LDPC control. At N around 320, the random control is slightly better
for fixed overhead 0 and 1.
