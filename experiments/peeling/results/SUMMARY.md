# Peeling Sweep Summary

These experiments measure only the sparse peel graph and inactivation schedule.
They do not model dense solver cost or dense residual success probability.
`N` is the number of blocks; payload bytes and memory pressure scale as
`N * block_bytes`.

## Quiet-Window Block-Operation Calibration (2026-06-10)

Calibration files:

- `cold_xor_calibration.csv`: cold-pool block-XOR baseline.
- `muladd_calibration.csv`: `gf256_muladd_mem` vs XOR on the same cold-pool
  protocol.
- `fanin_gather.csv`: chained pairwise XORs vs fused k-ary gather XOR.

Protocol: `xor_bench`, 4 GiB cold pool, `--target-gib 32`, 5 repeats, seed 1,
single-core `taskset -c 8`.  The cold XOR baseline reports 15.2 GiB/s at
1280-byte blocks, 33.7 GiB/s at 100 KiB, and 35.1 GiB/s at 1 MiB.

`gf256_muladd_mem` is much slower than XOR only for 1280-byte blocks
(`muladd_xor_ratio=1.83`).  At 100 KiB it is close to XOR (`1.10`), and at
1 MiB it is effectively bandwidth-equivalent (`1.00`).  Heavy-row cost models
should therefore not apply a large muladd penalty at large block sizes.

The fused gather-XOR experiment is not strong enough to justify a production
gather-kernel phase.  The session-plan kill gate was add8 >= 1.15x at 1 MiB;
the measured add8 ratio is 1.127x at 1 MiB.  Gather remains useful as evidence
for the cost model (`k=8` is 1.34x at 100 KiB and 2.22x at 1280 bytes), but
the production implementation branch should stop unless a real schedule replay
shows a stronger bottleneck-specific win.

## Fresh Random-Only Sweep

Original broad-sweep protocol:

- Binary: current `experiments/peeling/peel_sweep`.
- Structures: the then-current 111 random-compatible structures.
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

In the original broad sweep, `lt_m2_c256_fold` was the strongest pure-random
structure found. The follow-up high-cap sweep below supersedes that large-N
recommendation.

## High-Cap LT Follow-Up

After `lt_m2_c256_fold` won the original zero-overhead pure-random comparison,
the harness added `lt_m2_c512`, `lt_m2_c1024`, `lt_m2_c512_fold`, and
`lt_m2_c1024_fold`.

Focused protocol:

- Structures: `lt_m2_c128_fold`, `lt_m2_c256`, `lt_m2_c256_fold`,
  `lt_m2_c512`, `lt_m2_c512_fold`, `lt_m2_c1024`,
  `lt_m2_c1024_fold`, and `wirehair_rand`.
- Methods: current fast/frontier degree-2-component schedules
  `10,18,40,41,42,43,44,78,79,99`.
- Sizes: anchor `N=320,3200,32000`, N-jitter 10.
- Fixed overhead: 0.
- Trials: 100, paired matrix seeds, PDFs enabled.

Best method per structure, breaking equal mean residual count by lower
`total_us`:

| N | rank | structure | method | mean cols | sd | median | p95 | thr/s |
| ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| 320 | 1 | lt_m2_c1024_fold | hyb_rqbeam | 11.75 | 2.45 | 12 | 15 | 898.1 |
| 320 | 2 | wirehair_rand | hyb_rqbeam | 11.84 | 2.63 | 12 | 17 | 1309.1 |
| 320 | 3 | lt_m2_c512_fold | hyb_rqbeam | 12.02 | 2.30 | 12 | 16 | 891.3 |
| 3200 | 1 | lt_m2_c1024 | hyb_rqbeam | 35.72 | 6.16 | 35 | 47 | 47.9 |
| 3200 | 2 | lt_m2_c512_fold | rqcc_lowref | 35.95 | 6.22 | 35 | 46 | 106.2 |
| 3200 | 3 | lt_m2_c512 | raptorq_d2cc | 35.97 | 7.04 | 35 | 51 | 120.7 |
| 32000 | 1 | lt_m2_c1024_fold | rqd2_default | 106.11 | 18.60 | 104 | 140 | 4.8 |
| 32000 | 2 | lt_m2_c512_fold | rqcc_lowref | 108.45 | 18.31 | 104 | 141 | 7.0 |
| 32000 | 3 | lt_m2_c1024 | rqcc_lowref | 109.38 | 20.47 | 105 | 147 | 5.5 |
| 32000 | 4 | lt_m2_c256_fold | rqcc_lowref | 123.26 | 21.95 | 123 | 160 | 9.3 |

At `N=32000`, `lt_m2_c1024_fold` improves the previous best
`lt_m2_c256_fold` mean residual count by about 14 percent. `lt_m2_c512_fold`
is only about 2 percent worse than `lt_m2_c1024_fold` and is meaningfully
faster in this focused run, so it is the more practical current candidate if
peel-side runtime matters.

Audit caveats:

- This is an equal-row comparison, not an equal sparse-edge or symbol-XOR
  comparison. At `N=32000`, expected row degree is about 7.11 for
  `lt_m2_c256_fold`, 7.80 for `lt_m2_c512_fold`, and 8.48 for
  `lt_m2_c1024_fold`. Some quality gain is therefore purchased with more
  sparse references per row.
- A naive equal-edge rerun is not a valid replacement for the zero-overhead
  comparison because reducing high-cap row counts below `N` makes the system
  underdetermined and leaves thousands of residuals by construction.
- Matrix seeds are paired across methods within a structure. Different
  structures still use structure-name-derived seed streams, so structure
  comparisons are not common-random-number paired. Four separate seed-family
  reruns at `N=32000` kept the same ranking: `lt_m2_c1024_fold` averaged
  105.98 residual columns, `lt_m2_c512_fold` 109.86, `lt_m2_c256_fold`
  121.69, and `wirehair_rand` 207.57.
- For `N <= cap`, folded and non-folded high-cap LT structures clamp to the
  same distribution. Small-N differences among `c512`/`c1024` variants should
  be treated as Monte Carlo noise unless a common-random-number structure test
  is added.
- The high-cap run was focused on the fast/frontier method set. An all-method
  audit completed for `lt_m2_c256_fold` and `lt_m2_c512_fold` and found the
  same best methods. It was stopped in known-slow methods while processing
  `lt_m2_c1024_fold`; among completed methods, the best residual remained
  106.11.

At `N=320`, caps above `N` clamp to `N`. A same-matrix check with `N-jitter=0`
confirmed that `lt_m2_c512`, `lt_m2_c1024`, `lt_m2_c512_fold`, and
`lt_m2_c1024_fold` produce identical non-timing results. The differences in
the original `N=320` focused table were from structure-name-derived independent
matrix streams, not from a real distribution difference.

## N=320 Recipe Ablation Sweep

Follow-up protocol for the small-N recipe search:

- Binary: current `experiments/peeling/peel_sweep`.
- Structures: 52 N=320-focused structures, including cap sweeps, fold-scale
  sweeps, min-degree controls, tiny degree-1 mixtures, tail-alpha controls,
  explicit high-degree rows, degree-2/3/4 mass controls, `wirehair_rand`, and
  `binary_p50`.
- Methods: all 115 peeling/inactivation schedules.
- Size: fixed `N=320`, `N-jitter=0`.
- Fixed row overheads: 0, 1, 2, 4, and 8 packets.
- Trials: 4 seed families x 100 trials, reported below as combined 400-trial
  empirical distributions for each exact `structure + method + overhead` row.
- Matrix seeds: paired across methods inside each structure. Different
  structures still use structure-name-derived streams, so exact ties should be
  treated cautiously unless a common-random-number structure comparison is
  added.
- Production deterministic Wirehair rows remain excluded. `wirehair_rand` is
  the Wirehair row-weight distribution with fully random column choices.

Validation:

- `bash experiments/peeling/build.sh`
- `./experiments/peeling/peel_sweep --self-test`
- ASan/UBSan self-test build and run.
- ASan/UBSan smoke sweep over representative new structures and top methods.
- CSV integrity checks over all 119600 result rows: expected row counts, 100
  trials, fixed `N=320`, paired-random source, full 115-method coverage, PDFs
  present, expected row overheads, and no production Wirehair structure.

Top exact `structure + method` combinations by fixed row overhead:

| OH | rank | structure | method | mean cols | var | median | p95 | total_us |
| ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| 0 | 1 | lt_m2_c320_d2_003 | hyb_rqbeam | 11.918 | 5.20 | 12 | 16 | 987.0 |
| 0 | 2 | lt_c320_d2_003_d3_003_d4_003 | rqd2_default | 11.950 | 6.50 | 12 | 16 | 516.4 |
| 0 | 3 | lt_c320_d2_003_d3_003_d4_003 | raptorq_d2cc | 11.950 | 6.50 | 12 | 16 | 532.2 |
| 0 | 4 | lt_c320_d2_003_d3_003_d4_003 | ks_random | 11.950 | 6.50 | 12 | 16 | 847.8 |
| 0 | 5 | lt_c320_d2_003_d3_003_d4_003 | rqd2_minfill | 11.950 | 6.50 | 12 | 16 | 859.2 |
| 1 | 1 | lt_m2_c320_d3_003 | hyb_rqbeam | 11.375 | 6.16 | 11 | 16 | 920.4 |
| 1 | 2 | lt_m1_c320 | hyb_rqbeam | 11.398 | 6.32 | 11 | 16 | 879.0 |
| 1 | 3 | lt_m2_c320_d3_003 | rqd2_default | 11.408 | 6.20 | 11 | 16 | 522.0 |
| 1 | 4 | lt_m2_c320_d3_003 | raptorq_d2cc | 11.408 | 6.20 | 11 | 16 | 531.1 |
| 1 | 5 | lt_m2_c320_d3_003 | ks_boundary_min | 11.408 | 6.20 | 11 | 16 | 899.8 |
| 2 | 1 | lt_m2_c256_fold25 | hyb_rqbeam | 10.870 | 6.10 | 11 | 15 | 799.7 |
| 2 | 2 | lt_m2_c256_fold25 | raptorq_d2cc | 10.893 | 6.13 | 11 | 15 | 451.4 |
| 2 | 3 | lt_m2_c256_fold25 | rqd2_default | 10.893 | 6.13 | 11 | 15 | 466.9 |
| 2 | 4 | lt_m2_c256_fold25 | ks_boundary_min | 10.893 | 6.13 | 11 | 15 | 770.0 |
| 2 | 5 | lt_m2_c256_fold25 | ks_random | 10.893 | 6.13 | 11 | 15 | 772.9 |
| 4 | 1 | wirehair_rand | rqd2_default | 9.762 | 5.51 | 10 | 14 | 147.2 |
| 4 | 2 | wirehair_rand | raptorq_d2cc | 9.762 | 5.51 | 10 | 14 | 151.2 |
| 4 | 3 | wirehair_rand | rqcc_lowref | 9.762 | 5.51 | 10 | 14 | 151.2 |
| 4 | 4 | wirehair_rand | rqcc_minfill | 9.762 | 5.51 | 10 | 14 | 314.3 |
| 4 | 5 | wirehair_rand | rqcc_livemin | 9.762 | 5.51 | 10 | 14 | 317.0 |
| 8 | 1 | lt_c320_d2_003 | hyb_rqbeam | 7.997 | 5.97 | 8 | 12 | 839.4 |
| 8 | 2 | lt_c320_d2_003 | rqcc_dup | 8.015 | 5.97 | 8 | 12 | 987.9 |
| 8 | 3 | lt_c320_d2_003 | rqd2_default | 8.018 | 5.98 | 8 | 12 | 510.8 |
| 8 | 4 | lt_c320_d2_003 | raptorq_d2cc | 8.018 | 5.98 | 8 | 12 | 529.7 |
| 8 | 5 | lt_c320_d2_003 | rqd2_minfill | 8.018 | 5.98 | 8 | 12 | 871.2 |

Top distinct structures by their best exact method:

| OH | rank | structure | method | mean cols | var | median | p95 | total_us |
| ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| 0 | 1 | lt_m2_c320_d2_003 | hyb_rqbeam | 11.918 | 5.20 | 12 | 16 | 987.0 |
| 0 | 2 | lt_c320_d2_003_d3_003_d4_003 | rqd2_default | 11.950 | 6.50 | 12 | 16 | 516.4 |
| 0 | 3 | lt_c320_d3_003 | hyb_rqbeam | 12.033 | 6.59 | 12 | 17 | 923.5 |
| 0 | 4 | lt_m2_c256_fold | hyb_rqbeam | 12.037 | 6.04 | 12 | 16 | 868.8 |
| 0 | 5 | lt_m2_c128 | rqcc_maxrow | 12.050 | 6.60 | 12 | 17 | 974.0 |
| 1 | 1 | lt_m2_c320_d3_003 | hyb_rqbeam | 11.375 | 6.16 | 11 | 16 | 920.4 |
| 1 | 2 | lt_m1_c320 | hyb_rqbeam | 11.398 | 6.32 | 11 | 16 | 879.0 |
| 1 | 3 | wirehair_rand | hyb_rqbeam | 11.463 | 5.70 | 11 | 16 | 592.6 |
| 1 | 4 | lt_c320_d3_003 | hyb_rqbeam | 11.495 | 6.63 | 11 | 16 | 955.1 |
| 1 | 5 | d1mix_lt_p2 | hyb_rqbeam | 11.535 | 6.92 | 11 | 16 | 991.0 |
| 2 | 1 | lt_m2_c256_fold25 | hyb_rqbeam | 10.870 | 6.10 | 11 | 15 | 799.7 |
| 2 | 2 | lt_c320_d4_003 | hyb_rqbeam | 10.900 | 6.19 | 11 | 16 | 894.9 |
| 2 | 3 | lt_m2_c256 | hyb_rqbeam | 10.947 | 6.29 | 11 | 15 | 775.0 |
| 2 | 4 | lt_c320_d3_003 | hyb_rqbeam | 10.975 | 5.79 | 11 | 15 | 904.9 |
| 2 | 5 | lt_c320_d2_003 | hyb_rqbeam | 10.980 | 6.15 | 11 | 16 | 908.9 |
| 4 | 1 | wirehair_rand | rqd2_default | 9.762 | 5.51 | 10 | 14 | 147.2 |
| 4 | 2 | lt_m2_c128_fold | hyb_rqbeam | 9.910 | 5.82 | 10 | 14 | 607.0 |
| 4 | 3 | lt_c320_d2_003 | hyb_rqbeam | 9.915 | 5.66 | 10 | 14 | 934.9 |
| 4 | 4 | lt_c320_d3_003 | hyb_rqbeam | 10.002 | 6.48 | 10 | 14 | 919.3 |
| 4 | 5 | lt_c320_d2_008_d3_003_d4_003 | hyb_rqbeam | 10.055 | 6.90 | 10 | 15 | 1100.0 |
| 8 | 1 | lt_c320_d2_003 | hyb_rqbeam | 7.997 | 5.97 | 8 | 12 | 839.4 |
| 8 | 2 | lt_m1_c320 | rqd2_default | 8.172 | 5.36 | 8 | 12 | 532.4 |
| 8 | 3 | lt_c320_d3_003 | hyb_rqbeam | 8.190 | 5.86 | 8 | 13 | 850.0 |
| 8 | 4 | wirehair_rand | rqcc_lowref | 8.262 | 5.92 | 8 | 13 | 146.8 |
| 8 | 5 | lt_c320_d2_003_d3_003_d4_003 | hyb_rqbeam | 8.325 | 5.97 | 8 | 13 | 890.8 |

Method costs for the top methods:

- `hyb_rqbeam`: `O(rqcc_lowref then beam near tail)`.
- `rqd2_default`: `O(rows*degree + cols)`, largest degree-2 component/default
  tie.
- `raptorq_d2cc`: `O(rows*degree + cols)`, largest degree-2 component.
- `rqcc_lowref`: `O(rows*degree + cols)`, largest degree-2 component plus low
  refs.
- `rqcc_dup`: `O(rows*degree + cols + cc*rowrefs)`, largest degree-2
  component plus duplicate partners.

Ablation readout:

- At zero overhead, the cleanest new signal is a very small explicit
  degree-2 mass on the min-2 cap-320 LT family:
  `lt_m2_c320_d2_003 + hyb_rqbeam` averaged 11.918 residual columns.
- The fastest near-tie at zero overhead is
  `lt_c320_d2_003_d3_003_d4_003 + rqd2_default` at 11.950 residual columns and
  about half the measured harness time of the best-mean row.
- The winner changes with row overhead. `lt_m2_c320_d3_003` is best at +1,
  `lt_m2_c256_fold25` at +2, `wirehair_rand` at +4, and `lt_c320_d2_003` at
  +8. The differences among good sparse candidates are small: usually a few
  tenths of a residual column over 400 trials.
- Full cap-fold is not consistently best at N=320. A 25 percent folded-tail
  scale on cap 256 wins the +2 row case, while folded-tail scale 0 and exact
  caps are competitive elsewhere.
- Forced high-degree rows are a clear failure mode for this peel-only metric:
  the best explicit-high variant still leaves about 140 residual columns across
  all overheads tested.
- Tail-alpha/temperature variants are also poor at N=320, with best means
  around 21 to 24 residual columns.
- Dense `binary_p50` remains only a smoke baseline. Peeling has to inactivate
  nearly the whole matrix before rows become peelable.

## N=3200 Recipe Ablation Sweep

Follow-up protocol for the medium-N recipe search:

- Binary: current `experiments/peeling/peel_sweep`.
- Normal sparse structures: 76 structures, adding N=3200 cap sweeps,
  cap-fold sweeps, fold-scale controls, min-degree controls, tail-alpha
  controls, degree-2/3/4 explicit mass controls, and `wirehair_rand`.
- Methods: all 115 peeling/inactivation schedules for the normal sparse set.
- Size: fixed `N=3200`, `N-jitter=0`.
- Fixed row overheads: 0, 1, 2, 4, and 8 packets.
- Trials: 4 seed families x 100 trials, reported below as combined 400-trial
  empirical distributions for each exact `structure + method + overhead` row.
- Matrix seeds: paired across methods inside each structure. Different
  structures still use structure-name-derived streams, so exact ties should be
  treated cautiously unless a common-random-number structure comparison is
  added.
- Production deterministic Wirehair rows remain excluded. `wirehair_rand` is
  the Wirehair row-weight distribution with fully random column choices.
- Dense and explicit-high-row controls were run as limited-method smoke tests,
  not mixed into the normal sparse rankings.

Validation:

- `bash experiments/peeling/build.sh`
- `./experiments/peeling/peel_sweep --self-test`
- CSV integrity checks over all completed shards: normal `174800/174800`
  rows, high-row smoke `3520/3520` rows, binary smoke `140/140` rows, fixed
  `N=3200`, paired-random source, expected overheads, PDFs present, and no
  production Wirehair structure.
- A bug-fix pass found that explicit-high-row structures were using default
  `1/d` weights for their ordinary base rows instead of the intended LT
  weights. That was fixed and covered by self-tests. The normal sparse sweep
  did not include these high-row structures.

Top exact `structure + method` combinations by fixed row overhead:

| OH | rank | structure | method | mean cols | var | median | p95 | total_us |
| ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| 0 | 1 | lt_m2_c1920 | rqcc_dup | 35.023 | 39.44 | 34 | 47 | 57014.6 |
| 0 | 2 | lt_m2_c1920 | raptorq_d2cc | 35.030 | 39.39 | 34 | 47 | 26023.8 |
| 0 | 3 | lt_m2_c1920 | rqd2_default | 35.030 | 39.39 | 34 | 47 | 26234.0 |
| 0 | 4 | lt_m2_c1920 | ks_boundary_min | 35.030 | 39.39 | 34 | 47 | 56435.7 |
| 0 | 5 | lt_m2_c1920 | ks_lowref | 35.030 | 39.39 | 34 | 47 | 56813.2 |
| 1 | 1 | lt_m2_c3200 | rqd2_default | 34.752 | 37.07 | 34 | 45 | 40931.3 |
| 1 | 2 | lt_m2_c3200 | raptorq_d2cc | 34.752 | 37.07 | 34 | 45 | 40936.3 |
| 1 | 3 | lt_m2_c3200 | rqd2_minfill | 34.752 | 37.07 | 34 | 45 | 84698.8 |
| 1 | 4 | lt_m2_c3200 | ks_boundary_max | 34.752 | 37.07 | 34 | 45 | 85813.8 |
| 1 | 5 | lt_m2_c3200 | ks_lowref | 34.752 | 37.07 | 34 | 45 | 86184.1 |
| 2 | 1 | lt_m2_c1024_fold25 | rqd2_default | 34.380 | 38.82 | 34 | 45 | 15823.7 |
| 2 | 2 | lt_m2_c1024_fold25 | raptorq_d2cc | 34.380 | 38.82 | 34 | 45 | 15899.9 |
| 2 | 3 | lt_m2_c1024_fold25 | rqd2_minfill | 34.380 | 38.82 | 34 | 45 | 33698.1 |
| 2 | 4 | lt_m2_c1024_fold25 | ks_random | 34.380 | 38.82 | 34 | 45 | 33754.5 |
| 2 | 5 | lt_m2_c1024_fold25 | ks_boundary_min | 34.380 | 38.82 | 34 | 45 | 34468.3 |
| 4 | 1 | lt_m2_c1024 | rqd2_default | 33.240 | 38.93 | 32 | 44 | 15573.1 |
| 4 | 2 | lt_m2_c1024 | raptorq_d2cc | 33.240 | 38.93 | 32 | 44 | 15739.9 |
| 4 | 3 | lt_m2_c1024 | rqd2_minfill | 33.240 | 38.93 | 32 | 44 | 32085.5 |
| 4 | 4 | lt_m2_c1024 | rqcc_dup | 33.240 | 38.93 | 32 | 44 | 33139.9 |
| 4 | 5 | lt_m2_c1024 | ks_boundary_min | 33.240 | 38.93 | 32 | 44 | 33199.3 |
| 8 | 1 | lt_m2_c1280 | rqd2_default | 31.165 | 38.93 | 30 | 43 | 18528.8 |
| 8 | 2 | lt_m2_c1280 | raptorq_d2cc | 31.165 | 38.93 | 30 | 43 | 18694.7 |
| 8 | 3 | lt_m2_c1280 | rqd2_minfill | 31.165 | 38.93 | 30 | 43 | 39229.9 |
| 8 | 4 | lt_m2_c1280 | ks_random | 31.165 | 38.93 | 30 | 43 | 39421.3 |
| 8 | 5 | lt_m2_c1280 | ks_lowref | 31.165 | 38.93 | 30 | 43 | 39517.4 |

Top distinct structures by their best exact method:

| OH | rank | structure | method | mean cols | var | median | p95 | total_us |
| ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| 0 | 1 | lt_m2_c1920 | rqcc_dup | 35.023 | 39.44 | 34 | 47 | 57014.6 |
| 0 | 2 | lt_m1_c3200 | raptorq_d2cc | 35.365 | 40.57 | 35 | 47 | 41306.2 |
| 0 | 3 | lt_m2_c3200 | raptorq_d2cc | 35.468 | 37.35 | 35 | 47 | 40354.7 |
| 0 | 4 | lt_m2_c512_fold | rqd2_default | 35.485 | 41.92 | 34 | 48 | 10675.4 |
| 0 | 5 | lt_m2_c1024 | raptorq_d2cc | 35.528 | 38.65 | 35 | 46 | 15868.5 |
| 1 | 1 | lt_m2_c3200 | rqd2_default | 34.752 | 37.07 | 34 | 45 | 40931.3 |
| 1 | 2 | lt_m2_c1920_fold | rqd2_default | 34.825 | 45.06 | 34 | 48 | 26662.2 |
| 1 | 3 | lt_m2_c1920 | raptorq_d2cc | 34.858 | 33.63 | 35 | 45 | 25930.2 |
| 1 | 4 | lt_m2_c960 | rqcc_lowref | 34.917 | 42.28 | 34 | 48 | 15218.9 |
| 1 | 5 | lt_m2_c1024_fold0 | raptorq_d2cc | 34.990 | 33.82 | 34 | 46 | 15643.3 |
| 2 | 1 | lt_m2_c1024_fold25 | rqd2_default | 34.380 | 38.82 | 34 | 45 | 15823.7 |
| 2 | 2 | lt_m2_c1920 | rqcc_dup | 34.523 | 43.17 | 34 | 46 | 55372.6 |
| 2 | 3 | lt_m2_c1024_fold0 | raptorq_d2cc | 34.550 | 38.89 | 34 | 46 | 15443.4 |
| 2 | 4 | lt_m2_c640_fold | rqcc_dup | 34.583 | 42.71 | 34 | 47 | 32275.9 |
| 2 | 5 | lt_m2_c1024_fold50 | rqd2_default | 34.587 | 36.71 | 34 | 46 | 16069.5 |
| 4 | 1 | lt_m2_c1024 | rqd2_default | 33.240 | 38.93 | 32 | 44 | 15573.1 |
| 4 | 2 | lt_m2_c3200 | rqcc_dup | 33.375 | 44.30 | 32 | 46 | 76660.1 |
| 4 | 3 | lt_m2_c3200_fold | rqcc_dup | 33.405 | 39.59 | 33 | 45 | 80733.4 |
| 4 | 4 | lt_m2_c1024_fold25 | rqcc_dup | 33.455 | 45.53 | 33 | 46 | 34273.1 |
| 4 | 5 | lt_m2_c960_fold | raptorq_d2cc | 33.490 | 44.56 | 33 | 45 | 15746.1 |
| 8 | 1 | lt_m2_c1280 | rqd2_default | 31.165 | 38.93 | 30 | 43 | 18528.8 |
| 8 | 2 | lt_m2_c512_fold | raptorq_d2cc | 31.250 | 40.96 | 30 | 44 | 10576.1 |
| 8 | 3 | lt_m2_c1280_fold | raptorq_d2cc | 31.317 | 42.74 | 31 | 44 | 19467.0 |
| 8 | 4 | lt_m2_c1920_fold | rqd2_default | 31.440 | 42.72 | 30 | 43 | 26844.9 |
| 8 | 5 | lt_m2_c1920 | rqd2_default | 31.468 | 43.08 | 31 | 44 | 25143.4 |

Method costs for the top methods:

- `rqd2_default`: `O(rows*degree + cols)`, largest degree-2 component/default
  tie.
- `raptorq_d2cc`: `O(rows*degree + cols)`, largest degree-2 component.
- `rqcc_lowref`: `O(rows*degree + cols)`, largest degree-2 component plus low
  refs.
- `rqcc_dup`: `O(rows*degree + cols + cc*rowrefs)`, largest degree-2
  component plus duplicate partners.
- `rqd2_minfill`: `O(rows*degree + cols + rowrefs)`, largest degree-2
  component plus min-fill tie.
- `ks_boundary_min` and `ks_boundary_max`:
  `O(rows*degree + cols + cc*rowrefs)`, largest degree-2 component plus
  component-boundary tie.

Ablation readout:

- At N=3200, the strongest normal sparse recipes are high-cap min-2 LT
  variants. The strict zero-overhead winner is `lt_m2_c1920 + rqcc_dup` at
  35.023 residual columns, but `raptorq_d2cc`/`rqd2_default` are within 0.007
  residual columns and are about 2x faster in this harness.
- The preferred cap shifts with overhead: `c1920` at +0, `c3200` at +1,
  `c1024_fold25` at +2, `c1024` at +4, and `c1280` at +8.
- `wirehair_rand` is not competitive in this N=3200 cohort: at zero overhead
  its best exact method is `raptorq_d2cc` at 41.395 mean residual columns,
  versus 35.023 for the best high-cap LT recipe.
- The explicit degree-2/3/4 mass and tail-alpha controls do not improve the
  high-cap LT frontier at N=3200. Tail-alpha variants are especially poor:
  `sol_m2_alpha075_c3200` leaves about 530 residual columns at zero overhead.

Limited high-row smoke top rows:

| OH | rank | structure | method | mean cols | var | median | p95 | total_us |
| ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| 0 | 1 | lt_m2_c320_hi2_h160 | rqd2_default | 35.843 | 37.09 | 35 | 47 | 8137.9 |
| 0 | 2 | lt_m2_c320_hi2_h160 | raptorq_d2cc | 35.843 | 37.09 | 35 | 47 | 8359.7 |
| 0 | 3 | lt_m2_c320_hi2_h160 | hyb_rqbeam | 35.928 | 37.24 | 35 | 47 | 24816.0 |
| 1 | 1 | lt_m2_c3200_hi1_h1600 | rqd2_default | 34.972 | 41.06 | 34 | 47 | 40503.7 |
| 2 | 1 | lt_m2_c320_hi2_h160 | raptorq_d2cc | 34.885 | 44.40 | 34 | 47 | 8306.2 |
| 4 | 1 | lt_m2_c3200_hi1_h3200 | raptorq_d2cc | 33.398 | 41.94 | 32 | 45 | 44375.3 |
| 8 | 1 | lt_m2_c3200_hi2_h1600 | rqd2_default | 31.663 | 37.90 | 31 | 43 | 42482.4 |

Dense `binary_p50` smoke top rows:

| OH | rank | structure | method | mean cols | var | median | p95 | total_us |
| ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| 0 | 1 | binary_p50 | rqd2_default | 3171.613 | 0.61 | 3172 | 3173 | 1442434.0 |
| 0 | 2 | binary_p50 | raptorq_d2cc | 3171.613 | 0.61 | 3172 | 3173 | 1525862.6 |
| 0 | 3 | binary_p50 | allmin_default | 3171.620 | 0.55 | 3172 | 3173 | 1368721.9 |
| 1 | 1 | binary_p50 | allmin_default | 3171.622 | 0.60 | 3172 | 3173 | 1406812.2 |
| 2 | 1 | binary_p50 | allmin_default | 3171.610 | 0.53 | 3172 | 3173 | 1386085.7 |
| 4 | 1 | binary_p50 | raptorq_d2cc | 3171.645 | 0.59 | 3172 | 3173 | 1247087.7 |
| 8 | 1 | binary_p50 | allmin_default | 3171.605 | 0.60 | 3172 | 3173 | 1381012.8 |

The dense baseline behaves as a smoke test should: residuals stay near `N`
because almost every column must be inactivated before dense random rows become
peelable. It should not be used as a normal sparse-code candidate.

## N=3200 XOR-Cost Follow-Up

Residual count alone is not enough because larger caps buy lower dense
residuals with more sparse references per generated row. A follow-up run added
the following structural XOR metrics to `peel_sweep`:

- `row_xors_mu`: average `degree - 1` block XORs needed to generate one row.
- `matrix_xors_mu`: total row-generation block XORs for the received matrix.
- `solve_sparse_xors_mu`: block XORs from propagating peeled symbols into live
  rows. This includes rows that arrive after a referenced column has already
  peeled.
- `solve_dense_xors_est_mu`: dense residual block-XOR estimate,
  `residual_cols * max(residual_cols, residual_rows)`.
- `solve_total_xors_est_mu`: sparse solve XORs plus dense residual estimate.
- `combined_xors`: `matrix_xors_mu + solve_total_xors_est_mu`. This mixes
  encoder-side row generation and decoder-side solve work, so treat it as an
  end-to-end block-XOR proxy, not as a single-machine benchmark.

Validation:

- `bash experiments/peeling/build.sh`
- `./experiments/peeling/peel_sweep --self-test`
- Added regression coverage for sparse solve XOR accounting when a row consumes
  a column that peeled before the row was added.
- Frontier rerun: 27 structures x 13 methods x overheads 0/1/2/4/8 x 4 seed
  families x 100 trials. Residual means matched the original full sweep to
  floating-point roundoff.

Best residual rows with XOR metrics:

| OH | structure | method | mean cols | row XOR/packet | encode XORs | solve total | combined |
| ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| 0 | lt_m2_c1920 | rqcc_dup | 35.023 | 7.174 | 22956 | 23872 | 46828 |
| 1 | lt_m2_c3200 | ks_boundary_max | 34.752 | 7.597 | 24317 | 25170 | 49487 |
| 2 | lt_m2_c1024_fold25 | ks_boundary_max | 34.380 | 6.649 | 21290 | 22200 | 43490 |
| 4 | lt_m2_c1024 | ks_boundary_max | 33.240 | 6.511 | 20861 | 21784 | 42645 |
| 8 | lt_m2_c1280 | ks_boundary_max | 31.165 | 6.713 | 21537 | 22454 | 43991 |

`ks_boundary_max` often has the lowest data-XOR count among methods tied on
residual quality, but it is not the fastest schedule by CPU time because it
uses component-boundary scoring. If CPU time is the tie-breaker,
`rqd2_default` or `raptorq_d2cc` usually remain the practical choices.

Pareto readout, residual quality versus combined block-XOR proxy:

| OH | structure | method | mean cols | row XOR/packet | solve total | combined |
| ---: | --- | --- | ---: | ---: | ---: | ---: |
| 0 | lt_m2_c1920 | rqcc_dup | 35.023 | 7.174 | 23872 | 46828 |
| 0 | lt_m2_c1024 | ks_boundary_max | 35.528 | 6.507 | 21740 | 42561 |
| 0 | lt_m2_c512 | ks_boundary_max | 36.225 | 5.839 | 19699 | 38385 |
| 0 | wirehair_rand | ks_boundary_max | 41.395 | 4.731 | 16612 | 31750 |
| 1 | lt_m2_c3200 | ks_boundary_max | 34.752 | 7.597 | 25170 | 49487 |
| 1 | lt_m2_c960 | ks_boundary_max | 34.917 | 6.467 | 21618 | 42317 |
| 1 | lt_m2_c640 | ks_boundary_max | 35.325 | 6.068 | 20396 | 39821 |
| 1 | wirehair_rand | ks_boundary_max | 40.657 | 4.729 | 16597 | 31733 |
| 2 | lt_m2_c1024_fold25 | ks_boundary_max | 34.380 | 6.649 | 22200 | 43490 |
| 2 | lt_m2_c1024_fold0 | ks_boundary_max | 34.550 | 6.525 | 21825 | 42719 |
| 2 | lt_m2_c512 | ks_boundary_max | 34.835 | 5.858 | 19756 | 38512 |
| 2 | wirehair_rand | ks_boundary_max | 40.303 | 4.721 | 16603 | 31721 |
| 4 | lt_m2_c1024 | ks_boundary_max | 33.240 | 6.511 | 21784 | 42645 |
| 4 | lt_m2_c512 | ks_boundary_max | 34.025 | 5.816 | 19659 | 38294 |
| 4 | wirehair_rand | ks_boundary_max | 38.752 | 4.713 | 16544 | 31645 |
| 8 | lt_m2_c1280 | ks_boundary_max | 31.165 | 6.713 | 22454 | 43991 |
| 8 | lt_m2_c512_fold | ks_boundary_max | 31.250 | 6.691 | 22379 | 43845 |
| 8 | lt_m2_c1024_fold0 | ks_boundary_max | 31.582 | 6.498 | 21811 | 42656 |
| 8 | lt_m2_c192_fold | ks_boundary_max | 32.303 | 5.792 | 19613 | 38194 |

This changes the interpretation of the N=3200 sweep. If the objective is only
to minimize residual dense size, use the high-cap winners above. If block XOR
budget matters, `lt_m2_c512`/`lt_m2_c512_fold` and `lt_m2_c1024` variants are
more attractive tradeoffs: they usually give up less than one residual column
while saving roughly 10 to 18 percent of the combined block-XOR proxy versus
the highest-cap residual winner.

## Dense Binary Smoke Baseline

The harness includes `binary_p50`, a fully random dense binary matrix where
each row includes each column independently with probability 0.5. This is a
smoke-test baseline rather than a practical sparse-code candidate: rows have
average degree `N/2`, so storage, row updates, and real symbol XOR work scale
quadratically.

Input: `--structures binary_p50 --N 320 --N-jitter 10 --methods all
--trials 100 --overhead 0 --pdf`

All 115 methods completed. Because the matrix is dense, peeling has to
inactivate almost every column before rows become peelable; the residual dense
count is therefore around 300 columns for `N ~= 320`.

| rank | method | mean cols | sd | median | p95 | total_us | cost |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | --- |
| 1 | bm_minrand | 300.05 | 6.23 | 300 | 310 | 3256.8 | O(rows*degree) min row, random all-but-one |
| 2 | minrow_ovmax | 300.08 | 6.05 | 301 | 309 | 61091.3 | O(rows*degree + candidates) all min rows/max overlap |
| 3 | allmin_default | 300.11 | 5.98 | 300 | 309 | 3628.0 | O(rows*degree + candidates) all minimum rows + default score |
| 4 | allmin_ratio | 300.11 | 6.16 | 301 | 309 | 60569.6 | O(rows*degree + candidates*rowrefs) all minimum rows + d2/fill ratio |
| 5 | minrow_best | 300.11 | 5.98 | 300 | 309 | 61242.3 | O(rows*degree) best column among all minimum rows |
| 10 | raptorq_d2cc | 300.18 | 6.09 | 301 | 309 | 3636.7 | O(rows*degree + cols) largest degree-2 component |
| 11 | rqd2_default | 300.18 | 6.09 | 301 | 309 | 3712.7 | O(rows*degree + cols) largest degree-2 component/default tie |
| 19 | rq_minrow | 300.20 | 6.09 | 301 | 309 | 3518.6 | O(rows*degree) Raptor-style minimum row degree |
| 54 | default | 305.38 | 6.26 | 306 | 315 | 4822.2 | O(cols) cached degree-2 refs |
| 101 | rqcc_lowref | 309.09 | 6.26 | 310 | 319 | 3587.2 | O(rows*degree + cols) largest degree-2 component + low refs |
| 113 | global_lowref | 309.99 | 6.21 | 311 | 320 | 3496.9 | O(cols) all columns + low reference count |

The dense baseline is useful because it behaves very differently from the
sparse LT-style matrices: algorithm choices only move the residual by about 10
columns, while the absolute residual is close to `N`. That is expected for a
dense random matrix under a peel/inactivation schedule.

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

## Reduced Frontier Sweep

To cut down the experiment space, the reduced frontier sweep used the union of
structures and methods that had appeared on the best residual, runtime, and XOR
frontiers in the earlier runs. It then ran the full Cartesian product of that
reduced list at `N=1000`, `N=6400`, and `N=32000`.

Protocol: `N-jitter=10`, fixed overhead rows `0,1,2`, paired matrix seeds, 4
seed families, and 100 trials per seed family. Each aggregate row is therefore
400 trials. The compact aggregate CSV is
`experiments/peeling/results/reduced_frontier_N1000_6400_32000_j10.csv`.

Selected structures:
`wirehair_rand`, `rs_c001_d50_c128`, `lt_m2_c96_fold`, `lt_m2_c128_fold`,
`lt_m2_c192`, `lt_m2_c256`, `lt_m2_c256_fold`, `lt_m2_c256_fold0`,
`lt_m2_c256_fold25`, `lt_m2_c320`, `lt_m2_c320_d2_003`,
`lt_m2_c320_d3_003`, `lt_c320_d3_003`, `lt_m2_c512`, `lt_m2_c512_fold`,
`lt_m2_c640`, `lt_m2_c960`, `lt_m2_c1024`, `lt_m2_c1024_fold0`,
`lt_m2_c1024_fold25`, `lt_m2_c1024_fold50`, `lt_m2_c1280`,
`lt_m2_c1920`, `lt_m2_c2560_fold`, `lt_m2_c3200`.

Selected methods:
`default`, `random_tie`, `raptorq_d2cc`, `rqd2_default`, `rqd2_minfill`,
`top16_lowref`, `rqcc_lowref`, `rqcc_dup`, `rqcc_ratio`,
`ks_boundary_min`, `ks_boundary_max`, `ks_lowref`, `ks_random`,
`hyb_lowref10`, `hyb_rqbeam`.

Validation:

- 48 shard files checked: 12 each for N=1000, N=6400, N=32000 prefix, and
  N=32000 remainder.
- stderr logs were empty.
- All aggregate rows have 400 trials.
- N samples stayed inside the requested +/-10 band.
- Final aggregate row count is 3375:
  3 N values x 3 overheads x 25 structures x 15 methods.

Residual winners:

| N | OH | structure | method | mean cols | var | median | p95 | row XOR/packet | combined XORs | total_us |
| ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1000 | 0 | lt_m2_c256_fold | raptorq_d2cc | 20.63 | 15.42 | 20 | 28 | 5.81 | 12014 | 2384 |
| 1000 | 1 | lt_m2_c256_fold | raptorq_d2cc | 20.05 | 14.58 | 20 | 26 | 5.77 | 11946 | 2209 |
| 1000 | 2 | lt_m2_c640 | raptorq_d2cc | 19.57 | 14.37 | 19 | 26 | 6.04 | 12480 | 4437 |
| 6400 | 0 | lt_m2_c1920 | raptorq_d2cc | 48.62 | 71.01 | 48 | 64 | 7.06 | 92653 | 51520 |
| 6400 | 1 | lt_m2_c3200 | raptorq_d2cc | 48.48 | 82.26 | 47 | 64 | 7.68 | 100649 | 79920 |
| 6400 | 2 | lt_m2_c1024_fold0 | rqd2_default | 47.54 | 78.49 | 46 | 63 | 6.52 | 85510 | 32325 |
| 32000 | 0 | lt_m2_c2560_fold | rqd2_default | 104.85 | 453.32 | 101 | 143 | 8.37 | 545568 | 453523 |
| 32000 | 1 | lt_m2_c3200 | raptorq_d2cc | 104.42 | 384.11 | 101 | 142 | 7.64 | 500153 | 479702 |
| 32000 | 2 | lt_m2_c3200 | raptorq_d2cc | 103.12 | 414.41 | 100 | 142 | 7.66 | 501280 | 476483 |

Best combined-XOR rows within residual slack:

| N | OH | structure | method | mean cols | row XOR/packet | solve XORs | combined XORs | total_us |
| ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |
| 1000 | 0 | wirehair_rand | ks_boundary_max | 21.83 | 4.20 | 4583 | 8786 | 2119 |
| 1000 | 1 | wirehair_rand | ks_boundary_max | 21.27 | 4.24 | 4624 | 8866 | 1764 |
| 1000 | 2 | wirehair_rand | ks_boundary_max | 21.03 | 4.23 | 4636 | 8877 | 2039 |
| 6400 | 0 | lt_m2_c320 | ks_boundary_max | 53.43 | 5.36 | 36894 | 71206 | 28264 |
| 6400 | 1 | lt_m2_c320 | ks_boundary_max | 53.11 | 5.37 | 36957 | 71329 | 29488 |
| 6400 | 2 | lt_m2_c320 | ks_boundary_max | 51.93 | 5.37 | 36861 | 71213 | 29468 |
| 32000 | 0 | lt_m2_c960 | ks_boundary_max | 112.09 | 6.44 | 218228 | 424463 | 563118 |
| 32000 | 1 | lt_m2_c960 | ks_boundary_max | 110.39 | 6.46 | 218354 | 425031 | 551895 |
| 32000 | 2 | lt_m2_c960 | ks_boundary_max | 110.86 | 6.44 | 217880 | 423855 | 572682 |

Readout:

- Pure residual quality continues to favor higher-cap min-2 LT variants.
- The XOR frontier is different. Within a 10 percent or +1 residual-column
  slack window, lower-cap structures save substantial block XORs:
  `wirehair_rand` at N=1000, `lt_m2_c320` at N=6400, and `lt_m2_c960` at
  N=32000.
- For ties in residual count, `rqd2_default`/`raptorq_d2cc` remain the
  practical CPU choices. `ks_boundary_max` often trims solve/combined XORs
  slightly, but costs much more wall time.

## Boundary Approximation And XOR Calibration Sweep

Follow-up protocol: same reduced frontier style, but focused on
`wirehair_rand`, promising min-2 LT cap variants, and the boundary-method
families.  Each row is 100 paired trials with `N-jitter=10` and fixed packet
overhead `0,1,2`.  Aggregate CSV:
`experiments/peeling/results/boundary_xor_focus_N1000_6400_32000_j10.csv`.

XOR calibration used a solve-like random block-XOR benchmark:
`./experiments/peeling/xor_bench --sizes 1280,102400,1048576 --target-gib 32 --repeats 5`.

| block bytes | median us per block XOR | median GiB/s |
| ---: | ---: | ---: |
| 1280 | 0.022577343 | 105.601 |
| 102400 | 4.309740442 | 44.257 |
| 1048576 | 45.927485992 | 42.526 |

Note: these calibration values predate the `xor_bench` default working-set
change to 256 MiB. Regenerate before using the 1280-byte row for block-XOR
conversion; the old small-block run may have been cache-resident.

The derived `greedy_xor_eq_*` columns convert greedy peel scheduler time into
equivalent block XORs.  `matrix_solve_plus_greedy_xors_*` adds row generation
XORs, solve XOR estimates, and greedy scheduler time in the same units.

Boundary approximation check versus full `ks_boundary_max` over 81 paired
structure/N/OH comparisons:

| method | mean residual delta | mean solve XOR delta | mean greedy-time ratio | same mean residual rows |
| --- | ---: | ---: | ---: | ---: |
| ks_bmax_top16 | 0.0000 | 11.24 | 0.1817 | 81 / 81 |
| ks_bmax_top64 | 0.0000 | 1.11 | 0.1947 | 81 / 81 |
| hyb_bmax5 | 0.0000 | 89.27 | 0.1523 | 81 / 81 |
| hyb_bmax10 | 0.0000 | 81.54 | 0.1665 | 81 / 81 |
| rqd2_default | 0.0000 | 97.09 | 0.1431 | 81 / 81 |
| rqcc_lowref | 0.0101 | 490.44 | 0.1444 | 51 / 81 |
| raptorq_d2cc | 0.0000 | 489.86 | 0.1423 | 81 / 81 |

The first top-K implementation accidentally built the whole metric cache and
was not faster.  After changing limited boundary scoring to local
per-candidate metrics, top16/top64 preserved residual quality and ran in about
18-20 percent of full boundary-max greedy time.

Pure residual winners:

| N | OH | structure | method | mean cols | solve XORs | greedy us |
| ---: | ---: | --- | --- | ---: | ---: | ---: |
| 1000 | 0 | lt_m2_c640 | ks_bmax_top64 | 20.45 | 6388 | 622 |
| 1000 | 1 | lt_m2_c256_fold | ks_boundary_max | 19.79 | 5978 | 2796 |
| 1000 | 2 | lt_m2_c1024_fold0 | ks_boundary_max | 19.59 | 6748 | 3553 |
| 6400 | 0 | lt_m2_c1920 | ks_boundary_max | 48.80 | 47629 | 54274 |
| 6400 | 1 | lt_m2_c3200 | ks_boundary_max | 48.07 | 50465 | 68573 |
| 6400 | 2 | lt_m2_c640 | ks_boundary_max | 47.47 | 40774 | 25240 |
| 32000 | 0 | lt_m2_c3200 | ks_boundary_max | 103.26 | 254134 | 604143 |
| 32000 | 1 | lt_m2_c3200 | ks_boundary_max | 103.10 | 253904 | 612870 |
| 32000 | 2 | lt_m2_c2560_fold | ks_boundary_max | 103.79 | 275441 | 922097 |

Best total cost at 1280-byte pieces:

| N | OH | structure | method | mean cols | matrix+solve+greedy eq XORs | greedy eq XORs |
| ---: | ---: | --- | --- | ---: | ---: | ---: |
| 1000 | 0 | wirehair_rand | raptorq_d2cc | 21.68 | 18039 | 9226 |
| 1000 | 1 | wirehair_rand | rqd2_default | 21.51 | 20183 | 11317 |
| 1000 | 2 | wirehair_rand | raptorq_d2cc | 21.08 | 17631 | 8690 |
| 6400 | 0 | lt_m2_c320 | rqcc_lowref | 54.14 | 201592 | 130117 |
| 6400 | 1 | lt_m2_c320 | rqd2_default | 55.12 | 190384 | 118606 |
| 6400 | 2 | lt_m2_c320 | raptorq_d2cc | 51.78 | 188028 | 116728 |
| 32000 | 0 | lt_m2_c320 | rqd2_default | 138.28 | 1506881 | 1143735 |
| 32000 | 1 | lt_m2_c320 | raptorq_d2cc | 142.43 | 1476476 | 1112894 |
| 32000 | 2 | lt_m2_c320 | rqd2_default | 142.30 | 1529122 | 1164783 |

Best total cost at 100 KiB pieces:

| N | OH | structure | method | mean cols | matrix+solve+greedy eq XORs | greedy eq XORs |
| ---: | ---: | --- | --- | ---: | ---: | ---: |
| 1000 | 0 | wirehair_rand | ks_bmax_top16 | 21.68 | 8791 | 64.9 |
| 1000 | 1 | wirehair_rand | ks_bmax_top16 | 21.51 | 8923 | 70.0 |
| 1000 | 2 | wirehair_rand | ks_bmax_top16 | 21.08 | 8921 | 66.9 |
| 6400 | 0 | wirehair_rand | rqd2_default | 62.91 | 65045 | 962.5 |
| 6400 | 1 | wirehair_rand | rqd2_default | 61.31 | 64748 | 995.2 |
| 6400 | 2 | wirehair_rand | rqd2_default | 62.08 | 65281 | 1022.8 |
| 32000 | 0 | wirehair_rand | hyb_bmax10 | 204.43 | 353855 | 9489.7 |
| 32000 | 1 | wirehair_rand | hyb_bmax5 | 205.47 | 354056 | 9467.2 |
| 32000 | 2 | wirehair_rand | hyb_bmax10 | 206.89 | 356526 | 10475.0 |

Best total cost at 1 MiB pieces:

| N | OH | structure | method | mean cols | matrix+solve+greedy eq XORs | greedy eq XORs |
| ---: | ---: | --- | --- | ---: | ---: | ---: |
| 1000 | 0 | wirehair_rand | ks_bmax_top16 | 21.68 | 8732 | 6.1 |
| 1000 | 1 | wirehair_rand | ks_bmax_top16 | 21.51 | 8859 | 6.6 |
| 1000 | 2 | wirehair_rand | ks_bmax_top16 | 21.08 | 8860 | 6.3 |
| 6400 | 0 | wirehair_rand | ks_bmax_top16 | 62.91 | 64138 | 104.2 |
| 6400 | 1 | wirehair_rand | ks_bmax_top16 | 61.31 | 63812 | 107.5 |
| 6400 | 2 | wirehair_rand | ks_bmax_top16 | 62.08 | 64309 | 103.4 |
| 32000 | 0 | wirehair_rand | ks_bmax_top64 | 204.43 | 345250 | 1006.5 |
| 32000 | 1 | wirehair_rand | hyb_bmax5 | 205.47 | 345478 | 888.4 |
| 32000 | 2 | wirehair_rand | hyb_bmax10 | 206.89 | 347034 | 983.0 |

Readout:

- For 1280-byte pieces, scheduler CPU dominates quickly; cheap schedules on
  lower-XOR structures win the total-cost metric even with worse residuals.
- At 100 KiB, scheduler CPU is still visible at large N, but cheap Wirehair-like
  rows dominate total cost and the best schedule shifts between default and
  boundary hybrids.
- For 1 MiB pieces, scheduler CPU is nearly free in block-XOR units; the best
  rows favor lower row-generation and solve XOR counts, and top-K boundary
  scoring becomes attractive.
- Pure residual minimization still favors higher-cap min-2 LT structures, but
  those can be a poor total-cost choice once row-generation XORs and scheduler
  CPU are included.

## Degree Retune Cost Frontier

This pass targets `wirehair-vbk`: retune the random row degree distribution for
lower symbol-XOR traffic while keeping overhead and residual size visible.  The
broad screen covered all 191 pure random-row structures, excluding dense binary
and RaptorQ-LDPC controls, with methods `10,18,40,115,116,117,118`.
`N=1000` and `N=6400` completed for all structures at 50 trials; `N=32000` was
stopped after the low-degree, d1mix, uniform-mix, and robust-LT families because
the remaining speculative heavy-tail families were not likely cost-frontier
contenders and were much slower.

The confirmation run used the union of broad-screen frontier structures and
prior high-quality structures: 45 structures x 7 methods x 3 N anchors x 3
fixed-overhead values, 100 paired trials per row, `N-jitter=10`, with residual
PDFs.  Aggregate CSV:
`experiments/peeling/results/degree_retune_frontier_N1000_6400_32000_j10.csv`.

Pure residual winners:

| N | OH | structure | method | mean cols | var | median | p95 | row XOR/packet | solve XORs |
| ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1000 | 0 | lt_m2_c192_fold | raptorq_d2cc | 20.08 | 13.69 | 20 | 27 | 5.5788 | 5948 |
| 1000 | 1 | lt_m1_c256 | rqcc_lowref | 19.68 | 19.36 | 19 | 27 | 5.1770 | 5568 |
| 1000 | 2 | lt_m2_c1024_fold25 | rqd2_default | 19.40 | 10.76 | 19 | 25 | 6.3961 | 6653 |
| 6400 | 0 | lt_m2_c1920 | raptorq_d2cc | 48.80 | 53.29 | 48 | 62 | 7.1535 | 48091 |
| 6400 | 1 | lt_m2_c3200 | raptorq_d2cc | 48.07 | 83.17 | 46 | 65 | 7.5992 | 50944 |
| 6400 | 2 | lt_m1_c3200 | raptorq_d2cc | 47.05 | 78.32 | 47 | 61 | 7.6538 | 51244 |
| 32000 | 0 | lt_m2_c3200 | raptorq_d2cc | 103.26 | 348.94 | 101 | 137 | 7.6383 | 255225 |
| 32000 | 1 | lt_m2_c2560 | raptorq_d2cc | 101.41 | 396.81 | 98 | 135 | 7.4401 | 248672 |
| 32000 | 2 | lt_m2_c2560 | rqd2_default | 99.07 | 172.92 | 97 | 123 | 7.4076 | 246244 |

Best total cost at 1280-byte pieces:

| N | OH | structure | method | mean cols | var | p95 | row XOR/packet | solve XORs | total eq XORs |
| ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1000 | 0 | robust_d1_001_d2_020 | rqcc_lowref | 38.48 | 64.96 | 52 | 3.3089 | 4846 | 14722 |
| 1000 | 1 | lt_robust64 | rqcc_lowref | 32.49 | 56.25 | 45 | 3.4854 | 4627 | 14918 |
| 1000 | 2 | robust_d1_001_d2_012 | rqd2_default | 31.00 | 61.94 | 43 | 3.4320 | 4465 | 14480 |
| 6400 | 0 | robust_d1_001_d2_003 | rqcc_lowref | 114.85 | 359.48 | 146 | 3.6795 | 37116 | 105548 |
| 6400 | 1 | robust_d1_001_d2_003 | rqcc_lowref | 116.32 | 368.64 | 150 | 3.6827 | 37605 | 107940 |
| 6400 | 2 | robust_d1_001_d2_003 | rqcc_lowref | 117.11 | 351.19 | 149 | 3.6768 | 37862 | 109207 |
| 32000 | 0 | rs_c001_d50_c128 | rqcc_lowref | 220.45 | 355.32 | 252 | 4.5285 | 193785 | 879720 |
| 32000 | 1 | lt_m2_c96 | rqcc_lowref | 300.25 | 600.25 | 338 | 4.1840 | 224819 | 884498 |
| 32000 | 2 | lt_m2_c96 | rqcc_lowref | 302.08 | 823.69 | 351 | 4.1876 | 226564 | 897877 |

Best total cost at 100 KiB pieces:

| N | OH | structure | method | mean cols | var | p95 | row XOR/packet | solve XORs | total eq XORs |
| ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1000 | 0 | lt_m2_c16 | rqd2_default | 44.30 | 35.05 | 55 | 2.5402 | 4494 | 7079 |
| 1000 | 1 | lt_m1_c16 | rqd2_default | 43.81 | 37.45 | 53 | 2.5234 | 4484 | 7056 |
| 1000 | 2 | lt_m2_c16 | rqd2_default | 43.38 | 46.10 | 55 | 2.5581 | 4532 | 7146 |
| 6400 | 0 | d1mix_lt_p005 | rqd2_default | 91.58 | 180.90 | 115 | 3.8032 | 32808 | 57543 |
| 6400 | 1 | lt_m1_c64 | rqd2_default | 91.25 | 175.30 | 109 | 3.8095 | 32837 | 57658 |
| 6400 | 2 | lt_no1_64 | rqd2_default | 91.52 | 151.54 | 112 | 3.8080 | 32947 | 57765 |
| 32000 | 0 | rs_c001_d50_c128 | rqcc_lowref | 220.45 | 355.32 | 252 | 4.5285 | 193785 | 341534 |
| 32000 | 1 | rs_c003_d10_c128 | rqd2_default | 196.54 | 465.70 | 236 | 4.6734 | 188480 | 342320 |
| 32000 | 2 | rs_c001_d50_c128 | rqcc_lowref | 218.09 | 470.46 | 252 | 4.5374 | 193592 | 342037 |

Best total cost at 1 MiB pieces:

| N | OH | structure | method | mean cols | var | p95 | row XOR/packet | solve XORs | total eq XORs |
| ---: | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1000 | 0 | lt_m1_c16 | ks_bmax_top16 | 44.20 | 37.82 | 55 | 2.5412 | 4489 | 7034 |
| 1000 | 1 | lt_m1_c16 | ks_bmax_top16 | 43.81 | 37.45 | 53 | 2.5234 | 4483 | 7013 |
| 1000 | 2 | lt_m2_c16 | ks_bmax_top16 | 43.38 | 46.10 | 55 | 2.5581 | 4531 | 7099 |
| 6400 | 0 | d1mix_lt_p005 | ks_bmax_top64 | 91.58 | 180.90 | 115 | 3.8032 | 32804 | 57182 |
| 6400 | 1 | lt_m1_c64 | ks_bmax_top64 | 91.25 | 175.30 | 109 | 3.8095 | 32830 | 57256 |
| 6400 | 2 | lt_no1_64 | ks_bmax_top64 | 91.52 | 151.54 | 112 | 3.8080 | 32938 | 57361 |
| 32000 | 0 | rs_c001_d50_c128 | ks_bmax_top64 | 220.45 | 355.32 | 252 | 4.5285 | 193562 | 338789 |
| 32000 | 1 | rs_c003_d10_c128 | ks_bmax_top16 | 196.54 | 465.70 | 236 | 4.6734 | 188434 | 338412 |
| 32000 | 2 | rs_c001_d50_c128 | ks_bmax_top16 | 218.09 | 470.46 | 252 | 4.5374 | 193346 | 338892 |

Readout:

- The cost-optimal degree distribution moves upward with N: cap16 around
  N=1000, cap64/d1mix around N=6400, and robust-soliton-ish cap128 around
  N=32000 for 100 KiB and 1 MiB pieces.
- For 1280-byte pieces, scheduler CPU dominates enough that `rqcc_lowref` on
  low-degree robust LT rows wins most total-cost rows, even though residuals are
  much larger than the pure-residual optimum.
- For 100 KiB and 1 MiB pieces, scheduler time is much less important and
  `ks_bmax_top16/top64` becomes useful again when solve XORs tie closely.
- These are experiment-harness results only. Shipping a degree-table retune
  still needs dense/heavy precode co-optimization and real codec failure-rate
  validation; the cap-only production prototype was already rejected in
  `wirehair-vbk`.

## Bounded Lookahead Cost Pass

This pass tested cheaper alternatives to exact bounded path lookahead in the
standalone peeling harness.  New methods:

- `av_top16`, `av_top64`: top-K/default candidates scored by a local exact
  avalanche-size estimate, without copying the full graph.
- `av_rqcc16`, `av_rqcc64`: largest degree-2 component candidates, capped by
  top-K/default score, then local avalanche estimate.
- `sim1_top16`, `sim1_rqcc16`: exact one-step copied simulation controls.
- `hyb_rqav5`, `hyb_rqav10`: `rqd2_default` until the residual tail, then
  `av_rqcc16`.

Validation:

- Normal `peel_sweep --self-test`.
- ASan/UBSan `peel_sweep --self-test`.
- Added a self-test that compares the local avalanche removed-column count
  against an exact one-step copied simulation on top candidates and during a
  short inactivation walk.
- N=32000 smoke showed exact one-step copy controls take multi-second trials,
  so they were kept out of the large-N main sweep.

Main aggregate:
`experiments/peeling/results/lookahead_cost_N1000_6400_32000_j10.csv`.
This is 10 frontier/prior structures x 13 methods x N=1000/6400/32000 x
fixed overhead 0/1/2, 50 paired trials per row, `N-jitter=10`, with PDFs.
Exact-copy control aggregate:
`experiments/peeling/results/lookahead_exact_controls_N1000_6400_j10.csv`.

Best pure residual winners in the main sweep:

| N | OH | structure | method | mean cols | p95 | greedy us |
| ---: | ---: | --- | --- | ---: | ---: | ---: |
| 1000 | 0 | lt_m2_c3200 | av_top64 | 18.30 | 23 | 5175 |
| 1000 | 1 | lt_m2_c2560 | av_top64 | 18.40 | 24 | 8072 |
| 1000 | 2 | lt_m2_c2560 | av_top64 | 17.90 | 26 | 9224 |
| 6400 | 0 | lt_m2_c3200 | av_top64 | 44.98 | 60 | 97458 |
| 6400 | 1 | lt_m2_c3200 | av_top64 | 42.98 | 57 | 142522 |
| 6400 | 2 | lt_m2_c2560 | av_top64 | 43.32 | 60 | 124273 |
| 32000 | 0 | lt_m2_c3200 | av_top64 | 93.88 | 127 | 809841 |
| 32000 | 1 | lt_m2_c2560 | av_top64 | 93.70 | 113 | 847280 |
| 32000 | 2 | lt_m2_c2560 | av_top64 | 91.54 | 119 | 846878 |

Old-vs-new residual frontier:

| N | OH | best old method | old cols | best new method | new cols |
| ---: | ---: | --- | ---: | --- | ---: |
| 1000 | 0 | rqd2_default | 20.44 | av_top64 | 18.30 |
| 1000 | 1 | rqd2_default | 20.36 | av_top64 | 18.40 |
| 1000 | 2 | rqd2_default | 19.52 | av_top64 | 17.90 |
| 6400 | 0 | rqd2_default | 49.34 | av_top64 | 44.98 |
| 6400 | 1 | rqd2_default | 47.68 | av_top64 | 42.98 |
| 6400 | 2 | rqd2_default | 47.58 | av_top64 | 43.32 |
| 32000 | 0 | rqd2_default | 100.66 | av_top64 | 93.88 |
| 32000 | 1 | ks_bmax_top16 | 99.80 | av_top64 | 93.70 |
| 32000 | 2 | rqd2_default | 97.94 | av_top64 | 91.54 |

Best 1 MiB total-cost winners:

| N | OH | structure | method | mean cols | total eq XORs |
| ---: | ---: | --- | --- | ---: | ---: |
| 1000 | 0 | lt_m2_c16 | ks_bmax_top16 | 43.56 | 6987 |
| 1000 | 1 | lt_m1_c16 | ks_bmax_top16 | 44.04 | 7020 |
| 1000 | 2 | lt_m2_c16 | ks_bmax_top16 | 43.40 | 7112 |
| 6400 | 0 | lt_m1_c64 | ks_bmax_top64 | 90.88 | 56945 |
| 6400 | 1 | lt_m1_c64 | ks_bmax_top64 | 91.08 | 57142 |
| 6400 | 2 | lt_no1_64 | ks_bmax_top16 | 91.00 | 57266 |
| 32000 | 0 | rs_c001_d50_c128 | ks_bmax_top16 | 224.48 | 340286 |
| 32000 | 1 | rs_c003_d10_c128 | ks_bmax_top16 | 193.80 | 337449 |
| 32000 | 2 | rs_c003_d10_c128 | hyb_rqav5 | 194.92 | 338435 |

Readout:

- `av_top64` is the best residual-quality schedule in this focused set,
  reducing mean residual columns by about 2 at N=1000, 4-5 at N=6400, and 6-7
  at N=32000 versus the best prior methods on the same selected structures.
- That residual gain is not cost-free.  Across paired rows, `av_top64` had a
  median greedy-time ratio of about 9.9x versus `rqd2_default`; `av_top16` was
  about 5.3x.  The capped rqcc avalanche variants did not improve residuals
  over the rqcc baseline enough to justify their cost.
- Exact one-step copy controls did not buy meaningful extra quality.  On the
  N=1000/6400 control set, `sim1_top16` was usually equal to or worse than
  `av_top16/av_top64`, while often 4x-220x slower.  `sim1_rqcc16` typically
  matched `rqcc_lowref` residuals while being tens to hundreds of times slower.
- Cost frontiers still favor `ks_bmax_top16/top64` for 1 MiB pieces and mostly
  `rqd2_default`/`rqcc_lowref` at smaller pieces.  `hyb_rqav5` only reached the
  1 MiB top row once, at N=32000/OH=2, and the margin was tiny.
- Conclusion: local avalanche lookahead is useful as a residual-quality
  diagnostic, but not a clear shipping scheduler improvement.  If we pursue it
  further, the next step should be a truly incremental avalanche/bucket
  implementation or a tail-only trigger tuned specifically for large block
  sizes.

## Cached Degree Sampler Optimization

Follow-up implementation pass: weighted random-row structures now build their
degree CDF once per generated matrix instead of recomputing degree weights for
every row.  This is a speed-only change; deterministic non-timing columns were
identical before and after on a paired `lt_m2_c320`, `rs_c001_d50_c128`, and
`wirehair_rand` smoke sweep.

Measured on `N=1000`, 5 trials, `ks_bmax_top16`, overhead 0:

| structure | build_us before | build_us after | change |
| --- | ---: | ---: | ---: |
| `lt_m2_c320` | 1193.2 | 381.6 | -68.0% |
| `rs_c001_d50_c128` | 1143.2 | 275.6 | -75.9% |
| `wirehair_rand` | 270.6 | 259.4 | -4.1% |

The same cache was added to the smaller codec-side v2 peel evaluator.  The
focused seed-table case
`N=32000 --peel-candidates 128 --trials 4 --bb-list 1280` stayed byte-identical
and improved from 121.65s to 116.83s (-4.0%).  Validation:

- `cmake --build build --target wirehair_v2_policy_test wirehair_v2_bench -j 8`
- `./build/codec/wirehair_v2_policy_test`
- `./experiments/peeling/build.sh`
- `./experiments/peeling/peel_sweep --self-test`
- ASan/UBSan `peel_sweep --self-test`

## Codec V2 Peel Scorer Shortcuts

Follow-up codec-side pass: the v2 peel evaluator no longer rescans incident
rows to count live references for todo columns.  For a todo column, every
incident row necessarily still has at least that column live, so the live-row
count is the stored column incidence count.  Boundary scoring for the top-K
degree-2 component candidates now uses the degree-2 endpoint counts already
collected during candidate construction instead of rescanning incident rows.

The focused seed-table case
`N=32000 --peel-candidates 128 --trials 4 --bb-list 1280` stayed byte-identical
and improved further from 116.83s to 51.35s.  Relative to the pre-cache
baseline from the previous pass, the same case moved from 121.65s to 51.35s
(-57.8%).  Validation:

- `cmake --build build --target wirehair_v2_policy_test wirehair_v2_bench -j 8`
- `./build/codec/wirehair_v2_policy_test`
- `wirehair_v2_bench seedtable --N 1000,6400,32000 --bb-list 1280,102400 --peel-candidates 8 --trials 2`
