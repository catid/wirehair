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
adds the named profile's two or three distinct precode columns, and addresses
the generator with the public packet ID.  Encoder, decoder, selected profile,
golden tests, and benchmarks all share this one mapping.  Its evaluation
domain is intentionally narrower than the structure builder's domain:
`BuildPrecodeSystem()` can
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

`compare --precode-profile certified|mixed|both` selects the WH2 equation
profile for the precode and cached-precode arms; the default is `certified`.
`both` replays the same message seed and packet-ID schedule through both WH2
profiles, alternates their execution order by trial, and emits separate
`v2_precode` and `v2_mixed` rows.  This option is deliberately distinct from
`--v2-profile`, which controls the older
V1-compatible wrapper comparison rather than WH2.  The mixed profile represents
payloads as byte pairs, so `mixed` and `both` reject any odd block-byte entry
before printing benchmark results.
Test builds accept `compare --mixed-mix-count 2|3` with a mixed precode
profile.  The encoder binds that recovery-row fanout into its serialized
profile and the decoder consumes the bound value, allowing end-to-end speed
and recovery comparisons without changing a named profile.

The `compare` cached arm has separate default-off local storage controls.  Use
`--precode-encoder-cache` to retain one exact message copy for direct
systematic encoding, or `--precode-decoder-cache` to retain accepted
systematic payloads and evaluate only missing source IDs during `Recover()`.
The decoder option reserves one message-sized address range and touches the
pages for received systematic packets, so common 10% loss costs roughly 90%
of a message in resident memory and adds one receive-side copy in exchange for
substantially faster recovery.  `--precode-cache` enables both controls.  All
three flags affect only the benchmark's cached-precode arm(s), `v2_cached`
and/or `v2_mixed_cached`: they do not alter packet bytes or serialized
profiles, and the installed public C API remains default-off pending a
versioned options surface.
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
./build/codec/wirehair_v2_bench compare --nlo 1000 --nhi 1000 --trials 30 --bb-list 1280,102400 --loss 0.10 --precode --precode-profile both --trial-details
./build/codec/wirehair_v2_bench compare --nlo 1000 --nhi 1000 --trials 10 --bb-list 1280 --loss 0.10 --precode --precode-profile mixed
./build/codec/wirehair_v2_bench compare --nlo 1000 --nhi 1000 --trials 10 --bb-list 1280 --loss 0.10 --precode-encoder-cache
./build/codec/wirehair_v2_bench compare --nlo 1000 --nhi 1000 --trials 10 --bb-list 1280 --loss 0.10 --precode-decoder-cache
./build/codec/wirehair_v2_bench compare --nlo 1000 --nhi 1000 --trials 10 --bb-list 1280 --loss 0.10 --precode-cache
./build/codec/wirehair_v2_bench precodefail --N 1000,3200,10000,32000,64000 --bb-list 1280 --overhead 0,1 --trials 100 --threads 16 --loss 0.10
./build/codec/wirehair_v2_bench precodefail --N 1000,3200,10000,32000,64000 --bb-list 1280 --overhead 0,1 --trials 100 --threads 16 --loss 0.10 --completion mixed --mix-count 2,3 --payload-e2e
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

The default `precodefail` completion is `certified`.  The `--completion mixed`
arm exercises the mixed GF(256)/GF(2^16) H12 solver and requires even block
bytes and the `periodic` heavy family.  Mix count 3 corresponds to
`WIREHAIR_V2_PROFILE_MIXED_2026_07`; mix count 2 corresponds to the distinct
opt-in `WIREHAIR_V2_PROFILE_MIXED_MIX2_2026_07`.  Neither changes the default
`WIREHAIR_V2_PROFILE_CURRENT`.  With `--mix-count 2,3`, both selected candidate
profiles receive the same packet-ID trace for each trial; optional payload E2E
checks also share their message and loss stream.  The
`precodefail_paired` comment reports all four paired outcome cells, the exact
two-sided McNemar p-value, marginal Wilson 95% intervals, and each arm's
systematic seed-attempt index.  Systematic selection is intentionally performed
per candidate, so differing attempt indices mean the result compares two valid
selected profiles rather than isolating only the mix-count variable.  The CSV
`failure_trials` field retains the corresponding per-arm trial indices.
By default each overhead value has an independently salted loss stream for
broad robustness sampling.  Test builds accept `--paired-overhead-stream` to
remove only that overhead salt: trial `t` then receives the same delivered-ID
prefix at every requested overhead, with each larger arm appending packets to
the smaller one.  This makes `failure_trials` sets directly comparable and
should produce a monotone failure curve for rank-only solves.

Test builds also expose `--source-hits N` and
`--packet-peel-seed-xor U32` on `precodefail`.  The latter perturbs the
profile-derived sparse-row seed without replacing its K-dependent derivation,
allowing solve-graph cost and rank to be screened independently of row degree.
Both `compare` and `precodefail` accept
`--packet-row-seed-multiplier U32` plus the optional
`--packet-row-seed-avalanche` test hooks.  These apply the same odd permutation
and avalanche to recovery-row seeds, allowing rank candidates found by
`precodefail` to be replayed through the end-to-end encoder/decoder benchmark.
The mixed
`compare` and `precodefail` arms additionally expose `--mixed-period P`,
`--mixed-gf16-rows 2|3|4`,
`--mixed-geometry frozen|shared-x`, `--mixed-residue-skew S`, and
`--mixed-residue-schedule constant|ramp|hashed`.  The thread-local period
override accepts H through 244.  A nonzero residue skew is restricted to the
`1..P-H`: each successive P-column block rotates its coefficient buckets by S,
changing cross-block collision groups while keeping the bucket count,
per-block balance, operation count, and every H-column encoder corner
invertible.  `--mixed-residue-schedule ramp` instead cycles the safe boundary
steps through `2,3,...,P-H,1`; its longer block-label sequence targets the
sharp K-dependent resonances left by every constant skew without changing the
same balance, operation-count, or corner-rank invariants.  The ramp requires
shared-X geometry, `P>H`, and zero constant skew.  The `hashed` schedule uses
the same safe step range in a deterministic 127-step sequence chosen so its
cumulative shift is coprime to periods 28 through 32; the complete block-label
cycle therefore extends beyond the supported K=64000 domain.  It has the same
geometry requirements as the ramp.  Test builds additionally accept
`--mixed-residue-hash-seed U32` with the hashed schedule to search alternate
safe step sequences; seed zero selects the long-cycle sequence described
above.  `--mixed-residue-hash-keyed` instead derives a distinct sequence from
the base hash seed and K, then advances to the first sequence whose 127-step
cumulative shift is coprime to P.  This preserves the full `127*P` block-label
cycle while breaking the fixed schedule/K coupling that produces sharp
cross-K resonances.  With shared-X hashed scheduling and `P>H`,
`--mixed-independent-extension-residues` gives the GF(2^16) rows a second,
distinct full-cycle sequence derived from the active keyed seed.  This breaks
rank collisions shared by the GF(256) and GF(2^16) rows at the cost of one
additional active-value XOR pass; the encoder, decoder, projection oracle,
and verifier all use the same extension schedule.  Test-hook `compare` and
`precodefail` builds accept `--mixed-gf256-rows 10|11|12`; the extra subfield
rows use the exhaustively corner-safe Y=11 and Y=46 coordinates, preserve the
frozen rows' X=12 coordinate window, and are only available with `shared-x`.
This provides H15 11+4 and H16 12+4 experiments without adding another
GF(2^16) row or changing the default 10+2 production profile.  The H16 pair
is restricted to P31/P32, whose independently scheduled corners are
nonsingular for every K=2..64000.  Neither period is a uniform promotion:
independent normalized-graph holdouts found severe K-dependent recovery
resonances (including a P32 failure concentration at K=11), while
cross-payload holdouts found additional graph-dependent resonances at
different payload widths.  The two periods remain available to research a
versioned payload-independent graph seed or K-dependent selection; other
periods are rejected because they are singular, allow the Y=46 coordinate to
overlap the shared-X window, or produced worse recovery resonances despite a
nonsingular corner.
`precodefail` also accepts
`--mixed-extension-residue-seed-xor U32` with this mode to screen alternate
full-cycle extension derivations without changing the base GF(256) schedule;
the default XOR is 78.
The test-only
`--packet-peel-seed-table normalized-h15-v1` selects an offline-tuned packet
seed XOR at the 23 hard block counts in `[4,41]` and leaves every other K at
zero.  It is intentionally restricted to the exact normalized seed-profile
1280, H15 11+4/P32, keyed-hash-68, independent-extension-XOR-78, mix-2
experiment used to derive it.  Each CSV row appends the active XOR so a
multi-K run remains reproducible.  Three independent hard-loss stages totaling
5.796 million selected-arm trials on those K reduced exact-K failures by
98.5% versus salt zero and by 87.0% versus the best K-dependent H16/D13/D14
row-count comparator.  A separate 6.9-million-trial-per-arm cross-payload
holdout at block sizes 2/64/256/1280/4096 reproduced 98.48--98.52% reductions
at every width.  A ten-pair core-pinned 1280-byte full-payload solve check at
overhead four found a 0.992 aggregate time ratio (0.991 median pair ratio),
with 0.8% fewer block XORs and 1.0% fewer muladds.  This remains a benchmark
hook rather than a named wire profile; the v2 table below supersedes it for
the final all-K validation.

`--packet-peel-seed-table normalized-h15-v2` keeps all v1 entries immutable
and adds six large-K graph fixes: `1683:19`, `15182:98`, `21394:26`,
`24432:75`, `34207:213`, and `62039:2`.  The selection is deliberately sparse.
An exhaustive K=2..64000 hard-loss sweep (burst, adversarial, and repair-only
at losses .35 and .50, 1,919,970 trials per arm) found 1,284 failures for salt
zero, 1,344 for global salt one, and 1,372 for global salt two.  The shifted
graphs repaired nearly every base failure but created a similarly sized,
mostly disjoint set, so neither global shift is an improvement.  H16 12+4/P32
reduced that base total to 1,199, but introduced 196 failures that H15 avoided
and required about 7% more block muladds; it is therefore not selected either.

The final v2 salts were chosen jointly for recovery and payload solve speed.
A fresh two-seed holdout covered block sizes 2/64/256/1280/4096, IID, burst,
permutation, systematic-first, repair-only, and adversarial schedules, and
losses .05/.10/.25/.35/.50/.75.  Across 72,000 trials at each of the six K
values, the selected graphs reduced 23,341 base failures to 40 (99.83%) and
heavy shortfalls from 146 to 10.  Alternating, physical-core-isolated
1280-byte full-payload solves at overhead four measured candidate/base cycle
ratios of 0.991, 0.971, 0.968, 1.002, 1.000, and 0.994 respectively; one
decode at each K aggregates to 0.993.  Thus the sparse table is about 0.7%
faster over those six cases while removing the recovery hotspots.  Like v1,
v2 remains a reproducible benchmark hook rather than a named or production
wire profile.

An independent-seed direct base-versus-v2 K=2..64000 sweep then repeated the
six hard-loss cells above with 1,919,970 trials per arm.  V2 reduced total
failures from 1,316 to 1,239 (5.85%).  At the 29 deliberately remapped K
values it repaired all 79 base failures and introduced two disjoint failures,
a 97.5% reduction.  All 383,820 paired rows at unlisted K were identical in
every deterministic result field, confirming that the sparse table leaves
the rest of the graph space unchanged.  Both arms reported zero internal
errors and 30 heavy shortfalls.

`--packet-peel-seed-table normalized-h15-v3` keeps all v1/v2 entries immutable
and adds nine residual all-K hotspot repairs: `10:139`, `20:140`, `11414:86`,
`48567:209`, `49312:52`, `49842:188`, `50281:121`, `51375:192`, and
`53503:238`.  A 256-salt discovery screen followed by two fresh-seed holdouts
tested burst, adversarial, repair-only, and IID schedules at losses .10, .35,
and .50.  Across 383,400 selected-arm trials on these nine K values, the final
salts reduced 13,074 salt-zero failures to 388 (33.7x); every schedule/loss
cell improved in the deepest holdout.  On that 324,000-trial slice they also
used 0.48% fewer block XORs and 1.98% fewer field multiply-adds.  D13 and D14
binary-row controls were both less reliable and did more aggregate work.  Like
v1/v2, v3 is a reproducible benchmark hook rather than a named or production
wire profile.  A subsequent fresh-seed K=2..64000 sweep repeated the six hard
loss cells with 1,919,970 trials per table.  V3 reduced v2's 1,254 failures to
1,249 by repairing all five failures at the nine new K values and introducing
none.  All 383,940 rows at unlisted K matched in every deterministic result
field; both tables reported 32 heavy shortfalls and zero internal errors.
Eight alternating, physical-core-isolated phases then measured 1280- and
4096-byte full-payload solves at overhead four.  V3 won every phase and had a
0.981 aggregate solve-time ratio versus v2, with 2.16% fewer block XORs and
6.66% fewer field multiply-adds.  V2 had one rank failure in the timing set;
v3 had none.

For decoders with at least
30000 solver columns and 1024-byte-or-larger blocks, the two residue families
are accumulated in one sequential scan when their combined scratch is at most
128 KiB; other shapes retain the lower-setup streamed passes.  With ten
GF(256) rows, the three- and four-row settings are H13/H14 experiments that
append one or two GF(2^16) completion rows;
`shared-x` builds all active extension rows from the same Cauchy column
coordinates as the ten subfield rows, while `frozen` preserves the named H12
profile's coefficients.  Together these hooks allow coefficient-bucket speed
and rank to be swept without changing a named or serialized profile.
The source-hit override accepts `[1,8]`; an omitted override is reported as
zero.  It changes only the experiment's staircase fanout, allowing its rank
and block-XOR tradeoff to be measured without changing the certified
block-count policy.
`precodefail --full-payload-solve`
uses each requested `--bb-list` value in the solver instead of the default
one- or two-byte rank proxy, making `solve_ms_mu` include the real RHS cost.
The adjacent `build_ms_mu`, `peel_ms_mu`, `project_ms_mu`,
`residual_ms_mu`, and `backsub_ms_mu` columns split that solve time into its
major phases so payload-size bottlenecks can be distinguished from matrix
construction and peeling costs.
These flags are benchmarking hooks, not wire-format selection controls.

`compare`, `precodecheck`, and `precodefail` accept deterministic common packet
schedules via
`--schedule iid|burst|permutation|systematic-first|repair-only|adversarial`.
The reported schedule seed reproduces the exact candidate-compatible ID
stream; `precodefail` uses the requested trace to select its exact K+overhead
delivered rows, while `--trial-details` prints each end-to-end arm's success,
overhead, and paired overhead delta.  Schedule construction is bounded even at
high requested loss.

Test-hook builds of `precodefail` also accept
`--binary-dense-rows 1..64` for an absolute Shuffle-2 binary-row count and
`--gf256-heavy-rows 1..128` with `--completion certified` for an absolute
all-GF(256) completion-row count.  These controls compare denser binary and
larger subfield alternatives on identical packet traces; they do not change
named or production profiles.
