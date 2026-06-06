# Peeling Experiments

This directory is reserved for standalone peeling-schedule experiments.
Experimental heuristics should live here instead of behind compile-time branches
in `WirehairCodec.cpp`.

Keep production codec changes separate from sweep harnesses:
- reuse public headers or copied experiment-only model code here;
- report residual unpeeled rows/cols and peel-side runtime separately from dense
  matrix behavior;
- remember that `N` is the number of blocks, while bytes are `N * block_bytes`.

Build:

```bash
bash experiments/peeling/build.sh
```

Self-test the sweep model and method invariants:

```bash
./experiments/peeling/peel_sweep --self-test
```

Calibrate block-XOR solve cost for 1280-byte, 100 KiB, and 1 MiB pieces:

```bash
./experiments/peeling/xor_bench --target-gib 16 --repeats 3
```

The benchmark reports median microseconds per `dst_block ^= src_block`, which
can convert `peel_sweep` timing columns into block-XOR-equivalent cost.

Example sweep:

```bash
./experiments/peeling/peel_sweep --N 128,2048,32000 --trials 40
```

`--N` values are anchor block counts.  By default each trial samples the actual
block count uniformly from `anchor +/- 10` and reports `actual_N_*` summary
columns.  Use `--N-jitter 0` for fixed-size trials.

Useful method selectors:

```bash
./experiments/peeling/peel_sweep --list-methods
./experiments/peeling/peel_sweep --list-structures
./experiments/peeling/peel_sweep --list-dimensions
./experiments/peeling/peel_sweep --methods legacy
./experiments/peeling/peel_sweep --methods combo
./experiments/peeling/peel_sweep --methods fast
```

`fast` is a convenience subset for broad scans; it skips full live-ref scans
and the slower row-action modes.

Matrix structure selectors:

```bash
./experiments/peeling/peel_sweep --structures random
./experiments/peeling/peel_sweep --structures all
./experiments/peeling/peel_sweep --structures uniform2_8,lt_trunc64,one10_u2_8
```

The harness intentionally excludes production Wirehair row-id generation.  The
`wirehair_rand` structure uses the Wirehair row-weight distribution with fully
random duplicate-free columns, making it comparable to the other random
structures.

Matrix seeds are paired by default: each method sees the same generated matrix
for the same structure, `N`, row count, trial, and base seed.  Use
`--matrix-seeds independent` only when intentionally measuring independent
Monte Carlo samples.

Rows default to `N + --overhead + ceil(N * overhead_pct / 100)`.  To sweep small
positive overheads in one command:

```bash
./experiments/peeling/peel_sweep --overhead-pct 0,1,2,5
```

The CSV reports `overhead`, `overhead_pct`, and sampled `matrix_rows_*`
columns for each aggregate row.

For empirical residual PDFs, add `--pdf`. The PDF columns are compact
`residual:probability` maps built from the trial samples.

Current broad comparison target:

```bash
./experiments/peeling/peel_sweep --N 320,3200,32000 --methods all \
    --structures all --trials 100 --threads 100 --overhead-pct 0,1,2,5 --pdf
```

Fixed packet overhead should be run as a separate suite:

```bash
./experiments/peeling/peel_sweep --N 320,3200,32000 --methods all \
    --structures all --trials 100 --threads 100 --overhead 0 --pdf
./experiments/peeling/peel_sweep --N 320,3200,32000 --methods all \
    --structures all --trials 100 --threads 100 --overhead 1 --pdf
./experiments/peeling/peel_sweep --N 320,3200,32000 --methods all \
    --structures all --trials 100 --threads 100 --overhead 2 --pdf
```

The harness models only the peel graph. It feeds generated random sparse rows,
allows singleton rows to avalanche immediately, then applies a selected greedy
inactivation schedule. Dense matrix success/failure is outside the measurement
target.

The self-test covers every method over random generated structures and over
small synthetic structures for min-1, min-2, min-3, degree-2 components,
LT-like, uniform, gaussian-like, degree-1 fraction, and max-degree-capped rows.
It checks row/column back references, live counts, final residual counts,
candidate pool membership, all-but-one row-action shape at every greedy step,
raw random-structure degree bounds before simulator feed, and deterministic
paired row generation.  It also checks that the N-jitter protocol is
deterministic for a fixed seed, varies across trials, stays within the requested
window, and keeps paired matrices identical across methods.
