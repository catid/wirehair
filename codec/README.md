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
continues to match the corresponding peel row.  These helpers retain the
experimental recovery-row model; the shipping message packet contract is the
separate version-4 path below.

`GeneratePacketMatrixRow()` defines that version-4 packet equation.  It uses
production Wirehair's integer `GeneratePeelRowWeight()` and row iterators,
adds exactly three precode columns, and addresses the generator with the
public packet ID.  Encoder, decoder, selected profile, golden tests, and
benchmarks all share this one mapping.  Its evaluation domain is intentionally
narrower than the structure builder's domain: `BuildPrecodeSystem()` can
represent any valid 16-bit total span (including one precode column and
precode spans above 65521), while exact packet generation/evaluation/solving
requires `2 <= P <= 65521`, `K + P <= 65535`, and a mix count no larger than
both three and `P`.  This upper bound is the last prime supported by the
version-4 `NextPrime16()` distinct-column iterator.  Structure-only tooling
may retain systems outside it, but message facades and packet APIs reject them
before evaluating a row or writing packet/solve output.
`ComputePrecodeValues()` computes the concrete staircase, dense, and heavy
intermediate block values for encoder-feasible precode systems, and
`ComputeRecoveryBlock()` evaluates one generated recovery row over the full
`[source | precode]` intermediate block vector.  `ComputeEncodedBlock()` wraps
that mapping for encoder block IDs: `block_id < K` copies a source block, and
`block_id >= K` evaluates recovery row `block_id - K`.  `PrecodeEncoder`
caches the computed precode parity blocks and exposes the same block-id encode
mapping without requiring callers to manage the parity buffer directly.
`MessagePrecodeEncoder` zero-pads an arbitrary byte-length message, then solves
the K systematic packet equations together with every certified precode
constraint for the complete intermediate vector.  This joint solve makes the
certified full-span Shuffle-2 construction encoder-feasible without changing
it to an identity corner.  The final systematic packet reports its original
short byte count; repair packets remain full blocks.

`MessagePrecodeDecoder` incrementally collects the same equations, peels the
binary system, projects the unresolved residual, and solves it exactly over
GF(256), including the Cauchy heavy coefficients and block right-hand sides.
After a bounded rank failure it retains the peeled-column affine projection and
reduced GF(256) pivots, then projects each new equation into that fixed basis.
This is equivalent to a cold solve even when the added row would change the
peeling order, because the checkpoint represents the complete original row
space rather than assuming that order remains optimal.  The receive payload
buffer is released only when the checkpoint plus one retry packet is no more
than 25% larger; residuals over `kMaxInactiveColumns` and checkpoints outside
that policy use the original cold re-solve fallback.  Both paths remain bounded
by K+1024 accepted packet ids.  Conflicting duplicates are rejected, and an
identical packet retries a transient solve OOM without allowing a later packet
to overwrite the pending equation.  The selected profile binds the precode and
packet contract versions, exact dimensions, seeds, seed-attempt index, salts,
dense-corner choice, and mix count; V1 modes reject profiles carrying any V2
contract state.  Deterministic seed selection supports K=2..64000; tests
exhaust K=2..2048 and pin representative large-K attempts through K=64000.

This is an erasure codec, not packet authentication.  Altered independent
equations can define a different valid message, so applications must verify a
cryptographic digest/MAC obtained from authenticated or otherwise trusted
metadata, or authenticate packets directly.  The expected value must not come
from the unauthenticated FEC packet stream.  Duplicate and
overdetermined consistency checks are useful diagnostics, not an integrity
substitute.

The benchmark's `precodecheck` mode exercises the explicit V2 encoder and
decoder and reports terminal categories, successful-trial overhead, and stage
throughput over a K/block-byte grid.  `compare --precode` adds the same V2 path
beside production and the V1-compatible wrapper arms.  `precodefail` runs a
threaded fixed-overhead V2 rank/failure grid and reports inactivation and solve
cost rather than relying on a peel-only proxy.
`wirehair_v2_resume_bench` pins itself to the first available CPU and runs 20
alternating cold/warm samples at K=1000 and K=10000.  Its fixed deficient
K-packet stream becomes full rank across exactly eight appended equations; the
executable fails unless cumulative resume time saves at least 50%, initial
K-packet time stays within 5%, and checkpoint memory stays within 25% of the
replaced receive payload/id buffers.

Dense-seed checks are handled by the benchmark's `densecheck` and `densetune`
modes, which run real encode/decode trials with candidate dense seeds and
report overhead.  Both modes report `-1` for mean overhead when every trial
failed, and `densetune` compares candidates for explicit N/block-byte lists.
`densecount` sweeps dense row-count deltas with paired trial seeds, and
`densegrid` sweeps dense row-count deltas plus dense seed candidates.
All loss-driven modes accept loss probabilities only in `[0, 0.99]` and print
the exact parsed value they execute.  `seedtable` accepts trial counts in
`[1, 1000000]` and records requested and completed counts in each output row.
Peel and dense candidate indices traverse keyed odd-step permutations of the
8-bit seed domain, so prefixes contain no repeated work and indices 0..255
enumerate every byte seed exactly once.  Tuning CSV rows report requested,
unique, and completed candidate counts explicitly; tiny-N dense candidate zero
retains its full 16-bit shipped seed before later candidates enter the byte
domain.  Requests above the 256-value domain remain visible as requested while
unique/completed counts report the bounded 256 candidates actually evaluated.
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

`wirehair_v2::Codec` is a full encoder/decoder facade.  Its existing
`InitializeEncoder()` and `InitializeDecoder()` modes select a
`SeedProfile`, inject the selected dense count, peel seed, and dense seed into
the underlying solver, then expose encode/decode/recover methods.  Production
codec files are not modified by this wrapper.  The facade still emits and
decodes V1-compatible recovery packets in those modes.  The explicit
`InitializePrecodeEncoder()` and `InitializePrecodeDecoder()` modes route the
facade through the matching experimental V2 packet format.  V1 and V2 modes
remain deliberately separate, successful mode changes discard the prior
mode's state, and failed changes preserve the last valid mode.

The installed API exposes only the precode/packet V2 mode through the separate
opaque `WirehairV2Codec` handle. `wirehair_v2_encoder_create()` publishes the
selected canonical 32-byte profile; `wirehair_v2_decoder_create()` reconstructs
all internal `SeedProfile` state from those bytes. The C++ RAII facade is in
`<wirehair/wirehair.hpp>`. The byte layout, typed parse errors, frozen profile
rules, and migration policy are specified in `V2_WIRE_PROFILE.md`.

Validation:

```bash
cmake --build build --target wirehair_v2_policy_test
./build/codec/wirehair_v2_policy_test
cmake --build build --target wirehair_v2_precode_decode_test
./build/codec/wirehair_v2_precode_decode_test
./build/codec/wirehair_v2_precode_decode_test --large
./build/codec/wirehair_v2_precode_decode_test --large-recovery
./build/codec/wirehair_v2_precode_solve_test
./build/codec/wirehair_v2_precode_roundtrip_test
./build/codec/wirehair_v2_precode_seed_selection_test
./build/codec/wirehair_v2_profile_test
```

Benchmark smoke checks:

```bash
cmake --build build --target wirehair_v2_bench
cmake --build build --target wirehair_v2_resume_bench
./build/codec/wirehair_v2_resume_bench
./build/codec/wirehair_v2_bench compare --nlo 64 --nhi 256 --trials 2 --bb-list 1280,102400,1048576 --max-message-mib 96
./build/codec/wirehair_v2_bench precodecheck --N 64,320,1000 --bb-list 16,1280 --trials 10 --loss 0.10
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

`precodefail` reports the actual production solver's binary-deficiency and
heavy-rank-gain histograms.  Its optional
`--heavy-family periodic,hashed` comparison replays identical packet schedules
against the frozen periodic Cauchy coefficients and an experiment-only
non-periodic full-column hash family.  `hashed` is diagnostic only: named and
serialized V2 profiles always use `periodic`.

`compare` and `precodecheck` accept deterministic common packet schedules via
`--schedule iid|burst|permutation|systematic-first|repair-only|adversarial`.
The reported schedule seed reproduces the exact candidate-compatible ID
stream; `--trial-details` prints each arm's success, overhead, and paired
overhead delta.  Schedule construction is bounded even at high requested loss.
