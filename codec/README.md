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
`GeneratePeelMatrixRow()` uses the version-2 row-index-addressable stream to
return any uint32 zero-based row without replaying earlier rows.
`GenerateRecoveryMatrixRows()` is the corresponding V2 recovery-row helper for
the real precode path: it keeps the source-column prefix identical to
`GeneratePeelMatrixRows()` for the same seed, then appends distinct columns from
the intermediate precode range `[K, K + S + D2 + H)`.
`GenerateRecoveryMatrixRow()` provides the matching constant-time seek for one
recovery block index.  Batch and single generation are bit-identical, the
source and precode-mix streams are independently keyed, and the source prefix
continues to match the corresponding peel row.
The real packet contract is version 3 and is exposed by
`GeneratePacketMatrixRow()`.  It uses production Wirehair's integer
`GeneratePeelRowWeight()` distribution, three precode mix columns in the
certified default (the bound experimental option permits one through three),
and the public packet ID directly.  IDs below K are therefore distinct
systematic equations rather than fixed intermediate columns, and repair IDs
cannot collide with them.
`ComputePrecodeValues()` computes the concrete staircase, dense, and heavy
intermediate block values for encoder-feasible precode systems, and
`ComputeRecoveryBlock()` evaluates one generated recovery row over the full
`[source | precode]` intermediate block vector.  `ComputeEncodedBlock()` wraps
that mapping for encoder block IDs: `block_id < K` copies a source block, and
`block_id >= K` evaluates recovery row `block_id - K`.  `PrecodeEncoder`
caches the computed precode parity blocks and exposes the same block-id encode
mapping without requiring callers to manage the parity buffer directly.
`MessagePrecodeEncoder` is the message-level adapter: it zero-pads the message
during initialization, solves the K systematic packet equations plus every certified
precode constraint for all intermediate blocks, and reports the partial byte
count for the final systematic packet while emitting full-size repairs.
`MessagePrecodeDecoder` incrementally collects the same versioned packet
equations, peels the binary system, solves the projected residual together with
the actual Cauchy GF(256) rows, and reconstructs the original systematic
equations.  The identity-corner option remains an explicit experiment; the
default full-span Shuffle-2 construction no longer depends on its singular
D2-only parity corner.
Its `IntermediateBlocks()` accessor exposes the solved full intermediate
vector; it intentionally does not retain or mislabel a second copy of the
original systematic message.

The encoder expands the public seed profile into a deterministic joint
precode/packet seed sequence and selects the first attempt whose K systematic
equations plus precode constraints have full rank.  Attempt zero is the base
profile; attempt `i` adds fixed full-period constants to both seeds.  The
published profile transports that attempt together with every derived matrix
field.  A decoder recomputes those fields from the base profile, bound options,
and attempt and rejects any mismatch before building the system; it does not
independently select or blindly trust the derived values.  Exhaustive
K=2..2048 tests lock the current fixup attempts (13 K values, maximum attempt
2), and representative large K values exercise the same contract.

Dense-seed checks are handled by the benchmark's `densecheck` and `densetune`
modes, which run real encode/decode trials with candidate dense seeds and
report overhead.  Both modes report `-1` for mean overhead when every trial
failed, and `densetune` compares candidates for explicit N/block-byte lists.
`densecount` sweeps dense row-count deltas with paired trial seeds, and
`densegrid` sweeps dense row-count deltas plus dense seed candidates.
All loss-driven modes accept loss probabilities only in `[0, 0.99]` and print
the exact parsed value they execute.  `seedtable` accepts trial counts in
`[1, 1000000]` and records requested and completed counts in each output row.
Allocating modes validate widened message sizes and configured caps before
emitting result headers.

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
codec files are not modified by this wrapper.  `InitializeEncoder()` and
`InitializeDecoder()` retain the V1-compatible path and reject profiles that
carry any selected or mixed V2 contract state; V2 packets must use the explicit
precode initializers.
`InitializePrecodeEncoder()` and `InitializePrecodeDecoder()` select the
version-3 certified-precode path.  A successfully selected profile binds the
packet/precode contract versions, exact precode dimensions and seeds, and all
matrix-affecting options.  A decoder given that profile inherits the bound
options when its options argument is null; explicit mismatches are rejected
before matrix construction.  `DenseCount` is validated and directly determines
the staircase width, so a serialized selected profile does not silently depend
on the local dense-count table.  The precode decoder supports reordered
systematic and repair packets, duplicate-ID validation, resumable
`Wirehair_NeedMore`, and exact message recovery.

The codec provides erasure recovery, not packet authentication.  Any K
independent packet equations, including altered right-hand sides, can define a
different valid message, so callers must authenticate packets or verify the
recovered message with a cryptographic digest/MAC.  Identical/conflicting
duplicates and packets checked against a completed or overdetermined solution
are consistency-checked, but those checks are not a substitute for integrity
metadata.

Validation:

```bash
cmake --build build --target wirehair_v2_policy_test
./build/codec/wirehair_v2_policy_test
cmake --build build --target wirehair_v2_precode_roundtrip_test
./build/codec/wirehair_v2_precode_roundtrip_test
```

Benchmark smoke checks:

```bash
cmake --build build --target wirehair_v2_bench
./build/codec/wirehair_v2_bench compare --nlo 64 --nhi 256 --trials 2 --bb-list 1280,102400,1048576 --max-message-mib 96
./build/codec/wirehair_v2_bench compare --nlo 64 --nhi 3200 --trials 10 --bb-list 17,1280,102400 --max-message-mib 128 --loss 0.10 --precode
./build/codec/wirehair_v2_bench precodefail --N 1000,3200,10000,32000,64000 --bb-list 1280 --overhead 0,1 --trials 100 --threads 16 --loss 0.10
./build/codec/wirehair_v2_bench compare --nlo 320 --nhi 320 --trials 20 --bb-list 1280 --v2-profile auto --auto-trials 8 --auto-min-delta 0.10
./build/codec/wirehair_v2_bench seedtable --N 320,1000,3200 --bb-list 1280,102400 --peel-candidates 8 --trials 2
./build/codec/wirehair_v2_bench compare --nlo 320 --nhi 320 --trials 20 --bb-list 102400 --dense-delta 4 --dense-candidate 6
./build/codec/wirehair_v2_bench peelcost --N 320,3200,32000 --bb-list 1280,102400 --structures policy,lt_m2_c96,lt_m2_c256,lt_m2_c512 --precode dense,ldpcdense,codecport --overhead 0,2 --trials 8
./build/codec/wirehair_v2_bench densecheck --N 7533 --bb 1280 --candidates 4 --trials 1
./build/codec/wirehair_v2_bench densetune --N 320,1000 --bb-list 1280 --candidates 4 --trials 2
./build/codec/wirehair_v2_bench densecount --N 320,1000 --bb-list 1280 --deltas -8,0,8 --trials 20
./build/codec/wirehair_v2_bench densegrid --N 320 --bb-list 102400 --deltas -16,0,4 --candidates 8 --trials 20
```
