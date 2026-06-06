# Wirehair v2 codec scaffold

This folder is isolated from the shipping Wirehair codec.  It is for the next
codec version and contains the peel policy, seed selector, H=0 peel evaluator,
and a codec wrapper used for seed-tuned experiments.

Initial solver policy:

| Block bytes | N <= 1000 | N <= 3200 | N <= 12000 | N > 12000 |
| --- | --- | --- | --- | --- |
| small, around 1280 bytes | `robust_d1_001_d2_012` | `lt_m1_c32` | `robust_d1_001_d2_003` | `lt_m2_c96` |
| medium, around 100 KiB | `lt_m1_c16` | `lt_m1_c32` | `lt_m1_c64` | `rs_c001_d50_c128` |
| large, around 1 MiB | `lt_m1_c16` | `lt_m1_c32` | `lt_m1_c64` | `rs_c003_d10_c128` |

All selected rows currently use `ks_bmax_top16`.

`SelectPeelingCodec()` returns a `PeelingCodec` struct that can be passed to
row-generation and peeling code.  The struct includes:

- `Solver`: `ks_bmax_top16` for the selected policy.
- `Structure`: the named experimental structure.
- `Family`: LT, robust degree-1/degree-2 mixture, or robust soliton.
- `MinDegree`/`MaxDegree`: row degree limits.
- `Degree1Mass`/`Degree2Mass`: explicit low-degree mass for robust mixtures.
- `RobustC`/`RobustDelta`: robust soliton parameters.
- `SolverCandidateLimit`: 16 for `ks_bmax_top16`.
- `FullyRandomRows`: currently true for all selected v2 rows.

The byte-class boundaries are geometric midpoints between the measured block
sizes: small below 12 KiB, medium from 12 KiB to below 320 KiB, and large at
320 KiB and above.  These thresholds should be retuned once the actual v2 codec
can replay full encode/decode schedules instead of the standalone peel model.

## Seeds

`SelectSeedProfile()` mirrors the original Wirehair seed-table shape:

- `PeelSeedBucket = N % 2048`
- table seed first
- exact-N peel/dense fixups take precedence
- block bytes still select the v2 peel policy

`TuneSeedProfile()` is the offline tuning hook.  It evaluates candidate peel
seeds for the same modulo bucket at H=0 (`rows = N`) and picks a non-outlier
candidate with low estimated XOR cost.  This is intended to prevent degenerate
N where the residual or solve XOR count spikes.

`BuildPeelSolvePlan()` is the shared policy-driven row-generation path.  It
combines `PeelingCodec` + `SeedProfile` into a deterministic matrix seed,
generates rows, and runs the selected peel solver/evaluator.  The seed-table
and benchmark commands use this path so matrix statistics stay consistent.

Dense-seed checks are handled by the benchmark's `densecheck` and `densetune`
modes, which run real encode/decode trials with candidate dense seeds and
report overhead.  Both modes report `-1` for mean overhead when every trial
failed, and `densetune` compares candidates for explicit N/block-byte lists.

## Codec Wrapper

`wirehair_v2::Codec` is a full encoder/decoder facade.  It selects a
`SeedProfile`, injects the selected dense count, peel seed, and dense seed into
the underlying solver, then exposes encode/decode/recover methods.  Production
codec files are not modified by this wrapper.

Validation:

```bash
cmake --build build --target wirehair_v2_policy_test
./build/codec/wirehair_v2_policy_test
```

Benchmark smoke checks:

```bash
cmake --build build --target wirehair_v2_bench
./build/codec/wirehair_v2_bench compare --nlo 64 --nhi 256 --trials 2 --bb-list 1280,102400,1048576 --max-message-mib 96
./build/codec/wirehair_v2_bench seedtable --N 320,1000,3200 --bb-list 1280,102400 --peel-candidates 8 --trials 2
./build/codec/wirehair_v2_bench densecheck --N 7533 --bb 1280 --candidates 4 --trials 1
./build/codec/wirehair_v2_bench densetune --N 320,1000 --bb-list 1280 --candidates 4 --trials 2
```
