# Wirehair v2 codec scaffold

This folder is isolated from the shipping Wirehair codec.  It is for the next
codec version and contains the peel policy, seed selector, H=0 peel evaluator,
and a codec wrapper used for v2 experiments.

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

**Model-transfer caveat (wirehair-b6h):** TuneSeedProfile scores candidates in
the v2 synthetic row model, but the winning seed cannot be injected into the
PRODUCTION solver as `_p_seed` without changing the seed-to-matrix mapping.
The bench `compare --v2-profile tuned` arm is therefore retired; use
`--v2-profile auto`, which A/B-validates the synthetic candidate against the
base profile with real trials, or `seedtable` for offline diagnostics.

`BuildPeelSolvePlan()` is the shared policy-driven row-generation path.  It
combines `PeelingCodec` + `SeedProfile` into a deterministic matrix seed,
generates rows, and runs the selected peel solver/evaluator.  The seed-table
and benchmark commands use this path so matrix statistics stay consistent.
`GenerateRecoveryMatrixRows()` is the corresponding V2 recovery-row helper for
the real precode path: it keeps the source-column prefix identical to
`GeneratePeelMatrixRows()` for the same seed, then appends distinct columns from
the intermediate precode range `[K, K + S + D2 + H)`.
`ComputePrecodeValues()` computes the concrete staircase, dense, and heavy
intermediate block values for encoder-feasible precode systems, and
`ComputeRecoveryBlock()` evaluates one generated recovery row over the full
`[source | precode]` intermediate block vector.

Dense-seed checks are handled by the benchmark's `densecheck` and `densetune`
modes, which run real encode/decode trials with candidate dense seeds and
report overhead.  Both modes report `-1` for mean overhead when every trial
failed, and `densetune` compares candidates for explicit N/block-byte lists.
`densecount` sweeps dense row-count deltas with paired trial seeds, and
`densegrid` sweeps dense row-count deltas plus dense seed candidates.

`peelcost` is an offline scoring helper for degree-distribution retuning.  It
keeps the production codec untouched, generates policy or named experimental
peel distributions, and reports a block-XOR cost proxy that combines sparse
matrix XORs with a dense/LDPC-dense/certified-codec-port precode width estimate:

```text
matrix_xors + precode_gen + N * solve_width + solve_width^2 / 2 + H * solve_width
```

This is not a decode-failure model, and it does not replay LDPC constraint-row
rank effects.  The certified precode proxy names are `codecport` and
`codecport_ic`; both use `MakeCertifiedParams()` for S/D2/H/N1 and differ only
in dense identity-corner accounting, so the `--heavy` override only applies to
the older dense/LDPC proxy models.  Use `experiments/precode/precode_sim` for
rank/failure validation before promoting any degree/precode candidate.

## Codec Wrapper

`wirehair_v2::Codec` is a full encoder/decoder facade.  It selects a
`SeedProfile`, injects the selected dense count, peel seed, and dense seed into
the underlying solver, then exposes encode/decode/recover methods.  Production
codec files are not modified by this wrapper.  The facade still emits and
decodes V1-compatible recovery packets; the certified V2 precode block-data
path is exposed by the helper APIs above until the matching V2 decoder is
ported.

Validation:

```bash
cmake --build build --target wirehair_v2_policy_test
./build/codec/wirehair_v2_policy_test
```

Benchmark smoke checks:

```bash
cmake --build build --target wirehair_v2_bench
./build/codec/wirehair_v2_bench compare --nlo 64 --nhi 256 --trials 2 --bb-list 1280,102400,1048576 --max-message-mib 96
./build/codec/wirehair_v2_bench compare --nlo 320 --nhi 320 --trials 20 --bb-list 1280 --v2-profile auto --auto-trials 8 --auto-min-delta 0.10
./build/codec/wirehair_v2_bench seedtable --N 320,1000,3200 --bb-list 1280,102400 --peel-candidates 8 --trials 2
./build/codec/wirehair_v2_bench compare --nlo 320 --nhi 320 --trials 20 --bb-list 102400 --dense-delta 4 --dense-candidate 6
./build/codec/wirehair_v2_bench peelcost --N 320,3200,32000 --bb-list 1280,102400 --structures policy,lt_m2_c96,lt_m2_c256,lt_m2_c512 --precode dense,ldpcdense,codecport --overhead 0,2 --trials 8
./build/codec/wirehair_v2_bench densecheck --N 7533 --bb 1280 --candidates 4 --trials 1
./build/codec/wirehair_v2_bench densetune --N 320,1000 --bb-list 1280 --candidates 4 --trials 2
./build/codec/wirehair_v2_bench densecount --N 320,1000 --bb-list 1280 --deltas -8,0,8 --trials 20
./build/codec/wirehair_v2_bench densegrid --N 320 --bb-list 102400 --deltas -16,0,4 --candidates 8 --trials 20
```
