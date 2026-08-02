# Wirehair2 benchmark contract v3

The machine-readable source of truth is
[`wh2_benchmark_contract_v3.json`](wh2_benchmark_contract_v3.json).  The
companion checker, [`wh2_benchmark_contract.py`](wh2_benchmark_contract.py),
validates its exact domains and rejects result ledgers that omit or replace a
cell.  This is intentionally a small contract and checker, not another
campaign framework.

## Target

Wirehair2 uses only GF(256).  The primary recovery metric is failure to decode
after exactly `K` packets have been received (overhead zero).  A selectable
final codec must have at most 1% failures both overall and in every declared K
band × loss/schedule stratum.  Overhead 1, 2, and 4 are retained as nested tail
metrics rather than being averaged away.

Decoder solve time is the primary speed metric.  Within Wirehair2, candidates
are compared using isolated solve time.  Wirehair1 peels while packets are fed,
so the fair cross-codec metric is full receive-to-success time; timing only its
last API call would hide Wirehair1 work.  Encoder initialization plus the first
K symbols is measured separately.

## Bounded development screen

Every candidate first receives the exact 30-K `short` cohort.  Recovery pairs
the three explicit public WH2 base seed attempts 0, 1, and 2 with three frozen
training-loss roots and uses these four strata:

- IID loss at 10%;
- burst loss at 50%;
- adversarial loss at 50%;
- repair-only loss at 50%.

That is 360 recovery cells per arm at the two-byte structure width.  The
separate speed screen is 384 common cells: 12 K anchors, widths 64 and 1280,
and sixteen paired repetitions.  Unlike raw recovery, development timing uses
the production base attempt zero with the three training roots.  Timing is a
speed measurement of one predeclared production construction, not another raw
seed census; raw attempts 0, 1, and 2 remain fully represented by the recovery
screen.  A candidate gets at most 7,200 seconds of wall time.  A timeout is
unresolved and non-selectable, not permission to shrink the domain.  Fatal
errors, byte mismatches, semantic/equation mismatches, and domain mismatches
stop immediately.

Each raw cell uses its declared base attempt exactly once: no retry, rescue, or
per-K fix.  A weak raw seed unit is one
`(K, base_seed_attempt, block_bytes)` with any overhead-zero failure.  Weak units,
multiplicity, repairs, and introductions are reported, but an individual weak
unit is not an architecture veto.  This preserves the distinction between
choosing the best architecture and repairing its seeds afterward.

For every overhead tail and every aggregate/band/width/stratum scope, a candidate
may exceed each control by at most 0.1 percentage point, rounded up to at least
one cell.  Find the eligible minimum aggregate overhead-zero count, then form a
recovery-equivalent set containing every eligible candidate within the same
0.1-percentage-point allowance of that minimum.  Order that set by the
equal-cell aggregate paired mean log ratio for isolated decoder solve, then
failures at overhead 1, 2, and 4, and finally stable arm identifier.  The
checker implements this ordering in `select_development_architecture`.  Raw
weak-seed count is descriptive and is deliberately absent from the tie-break
list.  The resulting canonical selection receipt binds both domain hashes,
both freeze hashes, the shared architecture-artifact hash, exact cardinalities,
the complete candidate roster, every eligible overhead-zero count, and the
ordered result.  It separately seals the winner's semantic identity from its
arm name, codec kind, and architecture descriptor.  This descriptor remains
fixed across later equation-preserving implementation optimization even when
the source commit or binary hash changes.

### Raw seed-policy evidence boundary

Closed uniform-raw architecture recovery uses
`wirehair.wh2.benchmark-freeze.v2`. Every frozen arm separately binds its
`construction_seed_basis` and `seed_schedule_sha256`: the certified head uses
`production-profile-v1`, Wirehair1 uses `not-applicable`, and raw experimental
arms use `uniform-raw-v1` plus schedule SHA-256
`90a98a3db207852dabdf5fb27573ef48bce52e0228cee4e291d96fa44ed509a7`.
The seed policy is deliberately alongside, not inside, the structure-only arm
descriptor, so structure and seed-schedule changes cannot masquerade as one
another.

Raw native records, logical ledgers, and execution receipts use their v2
schemas. Each raw row binds the paired precode/packet attempts, effective
seeds, and realized staircase, binary-dense, and GF(256)-heavy geometry; the
realized-construction hash covers all of those fields. The work/rank sidecar
independently emits the same identity, and combination requires an exact
native-to-sidecar trace and construction join for every raw cell. Missing or
forged policy fields and legacy v1 freezes that claim a closed raw arm name or
descriptor are rejected.

The four current raw descriptors form a closed descriptor-to-arm/geometry
map: D12/H11, D12/H12, D12/H13, and D13/H12, all periodic-Cauchy with mix
count three. Validation also reproduces the certified `GetDenseCount(K)`
staircase exactly for every K in 2..64,000, requires the certified source-hit
boundary (two below K=10,000 and three thereafter), and requires the dense
identity corner to remain disabled. Recomputing a realized-construction hash
therefore cannot relabel one descriptor with another arm or geometry.

Production and timing evidence remains on `benchmark-freeze.v1` and the
existing v1 native receipt schemas. Selecting a raw-named architecture does
not rewrite later repaired/validation production receipts; the schema boundary
distinguishes raw construction evidence from the arm's human-readable name.

## Final recovery phases

Only the selected architecture expands beyond the short screen.

1. `final_raw` is every K from 2 through 64,000, three training seed pairs,
   and the three 50% hard schedules.  It has exactly **575,991 cells per arm**.
   It must satisfy the 1% overall and per-band/stratum rules and remain
   noninferior to both frozen controls at all four overhead thresholds.
2. After the architecture is frozen, a separate deterministic map-derivation
   run uses the frozen training traces.  For each K it tries retry offsets 0
   through 255 in ascending order and keeps the first offset that decodes all
   three hard schedules crossed with all three training loss roots at overhead
   zero.  With production base attempt zero, public `seed_attempt` is
   `(base + retry_offset) mod 256`.  The 63,999-byte offset vector binds the
   exact training trace manifest, binary, arm descriptor, and source commit.
   `final_repaired` is a scored replay of that already-frozen map, not the map
   derivation run.  It repeats 575,991 cells and requires zero weak production
   K units.  Manual exceptions are forbidden.
3. `final_validation` repeats the every-K hard census on three disjoint loss
   roots while replaying exactly the frozen attempt for the single production
   construction at each K.  It requires zero overhead-zero weak K units and no
   unsupported or construction errors; the 1% rule is then implied.  Validation roots may not alter the
   repair procedure or map.
4. `cross_width_validation` adds 1,440 bounded cells across widths 64, 256,
   1280, and 4096 on the 30-K cohort, including IID and all hard strata.  A
   payload-dependent construction must report the realized construction at
   every width; a requested parameter that the codec ignores is not an
   experiment axis.

The raw architecture census pairs base attempts 0, 1, and 2 with three training
loss roots rather than crossing them.  The repaired and validation censuses use
the one mapped production attempt per K crossed with three loss roots.  Both
shapes remain 63,999 × 3 × 3 = 575,991 cells per arm.

Final acceptance is a joint decision, never a collection of unrelated
per-phase passes.  `final_raw`, `final_repaired`, `final_validation`,
`cross_width_validation`, and final timing contain the two controls plus the
same one selected candidate and reuse the exact source commit, host, codec
kinds, binaries, descriptors, repair-training trace, and production map
hashes.  `validate-final-continuity` also requires the exact development
selection receipt and rejects a different sole candidate.  It checks this
semantic descriptor identity—not merely the arm name—before the phase
summaries may be combined.  Every individual phase still freezes its exact
source commit and binary.

The public attempt has its existing codec meaning: it adds
`attempt*0x9e3779b97f4a7c15` modulo 2^64 to the certified precode seed and
`attempt*0x9e3779b9` modulo 2^32 to the packet peel seed.  For the certified
control, the checker reconstructs the canonical 32-byte `WHV2` profile using
`message_bytes=K*block_bytes` and authenticates it together with the frozen arm
descriptor.  Experimental candidates use a benchmark-only construction
receipt until the winning equations receive a new public profile ID; they must
not claim the existing certified profile ID.  Codec identities are reserved in
both directions: only the named Wirehair1 control may declare `wirehair1`, only
the named Wirehair2 head may declare `wirehair2_certified`, and a candidate is
either `wirehair2_experiment` or an explicitly measured `routed_composite`.

For each cell, convert `loss_ppm` to `double(loss_ppm) / 1,000,000` and derive
the trace seed as `loss_seed xor K*0x9e3779b97f4a7c15 xor
block_bytes*0xbf58476d1ce4e5b9`, modulo 2^64.  The JSON contract freezes the
complete SplitMix64/drop and hard-schedule algorithms.  Generate one prefix of
K+4 delivered packet IDs and replay its first K+h IDs at threshold h; rerunning
four independently randomized thresholds is invalid.  The trace receipt is
SHA-256 over all K+4 packet IDs encoded as consecutive little-endian uint32s.

## Common-cell and failure rules

The recovery cell key is:

```text
(phase, band, K, block_bytes, loss_ppm, schedule, trial,
 base_seed_attempt, loss_seed, overhead_cap)
```

The arm name is intentionally excluded.  Before any result is read, one
canonical freeze manifest becomes the sole authority for the exact contract,
source commit, ordered arm roster, per-arm binary and descriptor hashes, domain
hash, trace-manifest hash, repair-map hashes, commands, affinity, and host.
The command-line caller cannot substitute a smaller roster.  Every arm must
then provide exactly one result for every frozen key.  Missing,
duplicate, extra, seed-swapped, partial, wrong-band, wrong-loss, or roster-drift
evidence invalidates the arm.  Observed-key intersection is never a valid
denominator.

The ordered trace manifest has one canonical JSONL row containing
`ordinal`, `cell_sha256`, and `trace_sha256` for every arm-free cell.  Duplicate
trace hashes are valid because distinct small-K/loss cells can legitimately
deliver the same packet-ID prefix.  The freeze hashes the logical canonical
stream with schema, contract, domain, and phase separation.  The checker proves
that every arm replayed the pre-frozen cell-to-trace mapping; the native trace
generator remains responsible for checking the packet IDs against the frozen
SplitMix64 algorithm before it seals the manifest.

`success` carries the first decoded overhead in `[0,4]`.  `need_more_at_cap`,
`construct_failed`, and `unsupported` remain explicit rows, count as failures
at all retained thresholds, and contribute four to capped-overhead sums.
Fatal or internal errors abort instead of becoming scoreable rows.  An
architecture with a limited K range is selectable only as an explicitly
measured routed composite whose fallback participates in all common cells.
An `unsupported` row in either mandatory control invalidates architecture
selection rather than gifting the candidate an easier baseline; a measured
control construction failure remains a scored recovery outcome.
Every result also binds the realized public construction-attempt index, the
SHA-256 of the realized construction descriptor, and the arm's repair-map
SHA-256.  Raw WH2 phases require the cell's base attempt and an all-zero no-map
marker.  Repaired WH2 phases load the content-addressed map and require exact
`map[K]` replay; Wirehair1 declares construction `not_applicable`, attempt zero,
and the no-map marker.  The checker recomputes every realized construction
receipt rather than accepting a free-form digest.

Single-object evidence is canonical JSON.  JSONL may physically use LF or
CRLF, but hashes are over parsed, compact, key-sorted logical records separated
by LF.  Duplicate keys, noncanonical records, booleans masquerading as
integers, missing rows, and extra rows fail closed.

## Timing decision

Timing uses warm-cache, pinned four-slot panels.  With one candidate, each base
cell has seven A/A panels (one for every arm/scope used) and four A/B panels:
candidate versus Wirehair2 solve, candidate versus Wirehair1
receive-to-success, and candidate versus each control for encoder work.  Thus
development has exactly 384 × 11 = 4,224 timing rows and final timing has
3,600 × 11 = 39,600.  Missing, duplicate, extra, side-swapped, or omitted A/A
rows invalidate the receipt.

For the canonical panel key `(kind, scope, left_arm, right_arm)`, take the low
bit of its SHA-256.  Use ABBA exactly when that bit equals replicate parity and
BAAB otherwise.  Run one untimed fresh warm-up per logical side immediately
before the four measured slots.

The frozen protocol replaces each underpowered one-shot measured slot with a
deterministic batch.  The machine-readable block target is 65,536 and every
timing cell receipts

```text
invocations_per_slot = max(1, ceil(65536 / K))
```

Each slot constructs and executes exactly that many fresh invocations of its
declared scope, checks every invocation's stable identity and outcome, and
sums only the exact-scope duration reported by each invocation.  Setup outside
the declared scope is not added to the sum.  The count is arm-independent,
fixed before results, part of the timing cell key and domain hash, and may not
be calibrated or shortened from observed timings.  A sum that cannot fit in a
positive signed 63-bit nanosecond value is fatal rather than wrapped.  The
arm-free timing cell key is:

```text
(phase, band, K, block_bytes, loss_ppm, schedule, replicate,
 base_seed_attempt, loss_seed, fixed_received_overhead,
 invocations_per_slot)
```

Changing the batch count therefore cannot preserve either the cell identity
or its frozen trace/domain binding.

Contract v3 additionally freezes the execution geometry that removed the
heterogeneous-load artifact.  Timing uses exactly eight singleton-pinned
workers on distinct physical cores.  For one candidate, the 12 K values, two
widths, and eleven frozen panels form 264 homogeneous cohorts ordered by width,
K-set position, then panel ordinal.  A cohort contains only one K, width, and
panel.  Its sixteen development replicates run as two waves of eight jobs:
wave `w` contains replicates with `floor(replicate / 8) == w`, and the next
wave or cohort may not start until all eight jobs in the current wave finish.
The one-candidate final domain analogously has 1,650 cohorts and three waves;
both phase bounds reject out-of-range cohort indices.
For zero-based cohort and wave indices, replicate `r` runs on

```text
sorted_worker_cpus[(r % 8 + cohort_index + wave_index) % 8]
```

This Latin rotation balances replicate identity over worker slots without
mixing K, payload width, or panel load inside a wave.  Before freezing results,
the runner requires the eight cores to cover as many L3/LLC groups as the
eligible CPU affinity permits, then fills them deterministically.  Unique L3
groups are preferred, not required: this host exposes only seven eligible L3
groups, so an honest eight-worker run necessarily reuses one group.  Terminal
validation independently rechecks that workers, controller, and sampler occupy
distinct physical cores, then checks every record's CPU, timestamp, placement,
and wave barrier.  A complete row set with the wrong placement or overlap is
not selectable.  Maximum-LLC coverage itself is a runner pre-freeze check;
reconstructing that decision offline would also require freezing the complete
pre-pin eligible CPU roster.

For the four aggregate elapsed nanoseconds `e0..e3`, one panel observation is:

```text
ABBA: ((ln e0 - ln e1) + (ln e3 - ln e2)) / 2
BAAB: ((ln e1 - ln e0) + (ln e2 - ln e3)) / 2
```

The independent statistical unit is the replicate, not an individual K or one
of the two adjacent pairs.  Within each band × width, average all frozen K
observations for a replicate, then use the 16 or 24 replicate means with unbiased
sample variance and the exact two-sided Student-t critical value (2.131449545559323
or 2.0686576104190477).  Aggregate encoder results give every frozen K × width
cell equal weight within its replicate.  No rounding occurs before a decision.

Timing may not form an observed common-success intersection: every compared
arm must succeed on every predeclared cell.  A non-success row carries four
null durations, remains explicit, and makes the affected panel non-selectable;
it is never dropped.  Successful timings require the exact derived
`invocations_per_slot` and four integer aggregate durations in `[1,2^63-1]`.
Final WH2 timing replays and receipts the same frozen repair map and
construction attempt used by recovery validation.

Isolated solve uses a fresh decoder after the identical K+4 received prefix is
loaded.  Encoder timing covers initialization and exactly IDs 0 through K-1.
Receive-to-success replays the nested trace and includes all feed and final
recovery work through each arm's first successful byte recovery.

The fixed effective timing floor is `log1p(0.02)`.  Before an A/B decision,
both corresponding A/A 95% intervals must lie strictly inside the symmetric
range from `-log1p(0.02)` to `+log1p(0.02)`; otherwise the comparison is noisy
and non-selectable.  A/A drift can therefore reject evidence but can never
expand the allowed regression.  A speed win requires the A/B 95% upper
log-ratio bound to be strictly below the negative floor.  Noninferiority
requires it to be strictly below the positive floor.  A resolved slowdown has
its lower bound strictly above the positive floor and rejects a development
candidate; a successful result spanning a decision boundary is unresolved.

Final promotion requires:

- Wirehair2 isolated decoder solve faster than the frozen Wirehair2 control in
  every band × payload width;
- receive-to-success faster than Wirehair1 in every band × payload width;
- encoder timing noninferior to both controls in every band × payload width and
  faster than both controls under the equal-cell aggregate.

Run the frozen homogeneous eight-worker waves on pinned cores.  No worker may
share a physical core with another worker, the controller, or the sampler.  The
existing sole CPU/DIMM/EDAC sampler remains active for the complete load
interval and its thermal/error policy is authoritative.

## Checker commands

Validate the frozen contract and print its hashes/cardinalities:

```bash
PYTHONWARNINGS=error python3 bench/wh2_benchmark_contract.py describe
```

Validate a development JSONL ledger containing both frozen controls:

```bash
PYTHONWARNINGS=error python3 bench/wh2_benchmark_contract.py \
  validate-ledger --phase development \
  --freeze-manifest freeze.json --trace-manifest traces.jsonl results.jsonl
```

Validate the complete paired timing receipt:

```bash
PYTHONWARNINGS=error python3 bench/wh2_benchmark_contract.py \
  validate-timing --phase development \
  --freeze-manifest timing-freeze.json \
  --trace-manifest timing-traces.jsonl timing.jsonl
```

Repaired recovery and final timing additionally repeat
`--repair-map ARM=PATH` for each mapped WH2/composite arm.  The checked-in
native worker and short-screen runner emit and validate the development
raw-cell and timing streams.  Separate final-phase native emitters remain
required before promotion.  Aggregate `compare` or `precodefail` output remains
diagnostic and cannot satisfy this contract.

Before final promotion, run the continuity check with one argument for each
required freeze:

```bash
PYTHONWARNINGS=error python3 bench/wh2_benchmark_contract.py \
  validate-final-continuity \
  --selection-receipt development-selection.json \
  --freeze recovery:final_raw=final-raw-freeze.json \
  --freeze recovery:final_repaired=final-repaired-freeze.json \
  --freeze recovery:final_validation=final-validation-freeze.json \
  --freeze recovery:cross_width_validation=cross-width-freeze.json \
  --freeze timing:final=final-timing-freeze.json
```

The checker authenticates and replays a frozen retry vector; proving that its
entries came from the lowest-success search, validating generated packet IDs
against the specified SplitMix64 algorithm, and binding worker-affinity plus
thermal/EDAC terminal receipts are responsibilities of the narrow native
emitters.  Their raw derivation/trace evidence is required before any campaign
is promotable; a hand-authored map or trace manifest is only a checker fixture,
not benchmark evidence.
