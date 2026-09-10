# Serialized V2 wire profile

The V2 packet API carries its complete equation identity in a canonical
32-byte descriptor. An encoder selects and serializes this record; a decoder
needs only the record and packet `(id, payload)` pairs. No in-process
`wirehair_v2::SeedProfile` value is exchanged.

The descriptor identifies equations. It does not authenticate the descriptor,
packets, or recovered message. Applications must authenticate framing and
verify a trusted digest or MAC before accepting recovered bytes.

## Canonical encoding version 1

Every multibyte integer is unsigned little-endian. The record is exactly 32
bytes; shorter and longer inputs are rejected. This exact-length rule prevents
a version-1 parser from silently ignoring data that a newer sender considers
meaningful.

| Offset | Bytes | Field | Version-1 rule |
| ---: | ---: | --- | --- |
| 0 | 4 | magic | ASCII `WHV2` |
| 4 | 2 | encoding version | `WIREHAIR_V2_PROFILE_ENCODING_VERSION` (`1`) |
| 6 | 2 | encoded bytes | `32` |
| 8 | 8 | equation profile ID | a supported `WIREHAIR_V2_PROFILE_*` value |
| 16 | 8 | message bytes | exact nonzero original length |
| 24 | 4 | block bytes | `1..2^31-1` |
| 28 | 1 | seed attempt | certified: `0..255`; small K3: `0` |
| 29 | 3 | reserved | all zero |

The derived block count is `ceil(message_bytes / block_bytes)` and must be in
`2..64000`. It is intentionally not serialized a second time.

The golden record for a 117-byte message, 16-byte blocks, attempt zero, and the
current profile is:

```text
57 48 56 32 01 00 20 00 c9 f9 f4 47 bb 5b 29 4b
75 00 00 00 00 00 00 00 10 00 00 00 00 00 00 00
```

`wirehair_v2_profile_serialize()` and
`wirehair_v2_profile_deserialize()` are the only supported conversion between
this byte encoding and the host-native `WirehairV2Profile` ABI structure.
`wirehair_v2_profile_validate()` performs the same byte validation when the
host representation is not needed. Copying that C structure directly to the
wire is not portable.

## Supported certified equation profile

`WIREHAIR_V2_PROFILE_CERTIFIED_2026_07` has numeric ID
`4b295bbb47f4f9c9`. The ID is the first 64 bits of SHA-256 over this exact
canonical UTF-8/ASCII input, with no newline:

```text
wirehair:v2:precode-v2:packet-v4:certified-2026-07
```

The full digest is
`4b295bbb47f4f9c91ebf12ba77afc33cf1c7c36131d154aea03c320f0f13dcf4`.
As with the legacy profile IDs, this hash is a stable compatibility name, not a
security primitive.

The profile freezes all equation-affecting rules used by the public V2 API:

- the current legacy dense-count, dense-seed, and peel-seed selection tables;
- V2 precode contract 2 and packet-row contract 4;
- the certified staircase, Shuffle-2 dense, and 12-row Cauchy-heavy geometry;
- source-hit selection, packet degree sampling, packet-ID mapping, and exactly
  three distinct precode mix columns;
- the precode salt `0x763263707265636f` and recovery-row salt
  `0x76327265636f7665`;
- the full-span dense corner (the experimental identity corner is disabled);
- deterministic attempt stepping and the first selected attempt published by
  the encoder; and
- GF(256) arithmetic and all integer PRNG/shuffle rules consumed by those
  equations.

Builds compiled with the private `WH_SEED_KNOBS` experiment switch can change
base seeds or packet degree weights at runtime. They deliberately return
`WirehairV2_UnsupportedPlatform` for the named profile instead of publishing
or consuming incompatible equations under its ID.

The descriptor serializes only canonical inputs: profile ID, message length,
block size, and the selected attempt. The block count, base seeds, peel policy,
precode dimensions, expanded precode/packet seeds, salts, mix count, fixup
diagnostics, and tuning statistics are derived under the named profile and are
not duplicated on the wire. Publishing the current ID after changing a frozen
rule is a compatibility bug.

## Small K3 equation profile

`WIREHAIR_V2_PROFILE_SMALL_K3_2026_09` has ID `67c1043ecaa9e184`, the first
64 bits of SHA-256 over the exact name (no newline):

```text
wirehair:v2:small-k3:thue-morse-8-14-7:gf256-14d:2026-09
```

The full digest is
`67c1043ecaa9e1847ffb885eaf845b3af93ce52df8493bb1c6cf5428e863bdb9`.
It freezes the same GF(256), polynomial `0x14d`, Thue-Morse companion pair
`(8,14,7)` / `(9,14,7)`, packet-ID map, systematic identity rows and zero
tail padding as the qualified [WHK3 equations](SMALL_WIRE_PROFILES.md).
The two facades share one immutable 13,056-byte lookup; WHV2 framing and
source-ownership rules remain distinct from WHK3.

This profile requires exactly three source blocks, block bytes at most
67,108,864 and seed attempt zero. Other shapes return `InvalidDimensions`;
a nonzero attempt on a valid shape returns `BadSeed`. There is no seed search.
The prepared encoder basis is the original three source blocks, with private
padding for a partial final block. It is owned for both storage policies:
borrowing adds no second full-message cache, no extra construction allocation,
and no extra copy. Repair packets never read borrowed input; detach remains
allocation-free. Dependent contradictory packets report `Error` without
changing the retained decoder basis; recovery does not authenticate data.

On this development branch, ordinary constructors select this profile for its
supported K3 shapes. Other shapes, including K3 above the small-profile block
bound, continue to select the certified profile. Correctness and compatibility
checks, the scoped ordinary-path speed gate and retained paired recovery
qualification below pass. Old-path/opt-in performance regression checks remain
open.
The published WHK3 opt-in timings do **not** qualify this owned-basis path.
This is not an all-K performance or construction-seed claim.

### Ordinary K3 full-lifecycle speed

At source `d77a9db`, `wirehair.wh2.k3-ordinary-cost-r0` measured the actual
ordinary source-independent constructor and ordinary borrowed constructor
against actual WH1 and explicitly selected certified WH2. All 72 same-code
timing controls and all 72 candidate comparisons passed. Every candidate
comparison's upper 95% time-ratio bound was strictly below one; there was no 5% minimum
gain, treatment regression tolerance, sample filtering or order pooling.

Observed elapsed-time reductions versus WH1, spanning both measurement orders:

| Source policy | Block bytes | Encoder | Low-ID decoder | Distant-ID decoder |
|---|---:|---:|---:|---:|
| Independent | 2 | 62.55–62.59% | 94.35% | 89.85–89.92% |
| Independent | 64 | 70.11–70.16% | 94.07–94.24% | 90.53–90.57% |
| Independent | 1280 | 71.38–72.14% | 93.29–93.40% | 91.90–92.09% |
| Borrowed | 2 | 61.58–61.68% | 94.31–94.33% | 89.88–89.90% |
| Borrowed | 64 | 69.11–69.19% | 94.14–94.22% | 90.58% |
| Borrowed | 1280 | 70.71–71.84% | 93.29–93.36% | 91.90–92.03% |

The run used full three-block messages in the static library on one
GFNI-capable host. Encoder time includes create, descriptor output, 18 packets
and free. Decoder time includes create, feed through its own first success,
recover and free; both frozen streams needed three packets for every arm.
The two source policies share the decoder: its repeated observations are not
extra recovery samples. Both WH2 source policies were compared with borrowed-input
WH1 and borrowed certified-WH2 controls. The source-independent WH2 rows are
therefore not ownership-matched comparisons; owned-input WH1 was not timed in
this gate. The historical measurements and decisions are unchanged. New
ownership-matched qualification must use the corresponding owned WH1 control.
This screen does not establish partial-tail,
additional-width, shared-call, cold-start, non-GFNI or all-K performance, nor
does it isolate regressions in the preserved old profile or opt-in facade.

All 31,104 callbacks, 3,981,312 fresh codec lifecycles and 42,467,328 API calls
were retained, with every output, descriptor, endpoint and buffer guard
checked. Native, portable-arithmetic and sanitizer neutral fixtures agree;
portable and sanitizer builds were not scientific timing arms. Python 3.8
replayed the complete result and all 582 receipt pins exactly. A separately
written standard-library chronology, ledger and statistical audit reproduced
all 144 decisions before source HEAD or pinned documentation changed.

The immutable outcome is `/var/tmp/wh2-k3-ordinary-cost-r0`, raw SHA-256
`03ad310ea6d44e7a98ece8ba27c96ab8284741df38c43b14715be4f7c66503d2`.
This namespace is spent. No recovery-rate sample or full promotion claim is
made by this timing result.

### Ordinary K3 retained recovery comparison

At source `734150a`, `wirehair.wh2.k3-ordinary-recovery-r0` replayed the
retained cohort through actual ordinary WH2, separately for independent and
borrowed creation, alongside actual WH1 and explicitly selected certified
WH2. Failures with zero extra packets on the 6,144 retained loss traces:

| Codec | Failures / traces | Observed failure rate |
|---|---:|---:|
| Ordinary WH2 K3, each source policy | 0 / 6144 | 0% |
| WH1 | 14 / 6144 | 0.23% |
| Certified WH2 | 269 / 6144 | 4.38% |

Both ordinary policies fixed all 14 WH1 and 269 certified-WH2 failures,
introducing none against either control in this cohort. Each of the twelve
width/loss-schedule cells retained 512 traces and had no K3 failures. WH1
recovered all failures with one extra packet; certified WH2 still had five
failures after one extra packet, one after two, and none after three.

The separate 72 hard cases also had no K3 failures at zero overhead, versus
one WH1 and five certified-WH2 failures. K3 and WH1 recovered all 53 historical
prefixes at their original widths; certified WH2 remained unresolved on 43.
Those prefixes were never extended to rescue an unresolved decode. Hard and
historical cases are not included in the 6,144-trace denominator.

Every encoder was freed before receiver creation. Both ordinary policies
produced identical descriptors, packet hashes, feed statuses and first-success
counts. Each K3 packet matched an independent polynomial-GF(256) payload
oracle; first success matched independently computed rank, and each successful
decode recovered the original message twice with guarded buffers. Native,
portable-arithmetic and ASan/UBSan builds agreed on every one of the 6,269
cases and all four API routes, with 50,152 fresh codec handles per backend.

Python 3.8 independently replayed all 647 receipt pins and the complete raw
results. A separately written standard-library audit rebuilt the GF(256)
payload/rank oracle and reproduced the full chronology, all group/cell counts
and every paired fix/introduction list before HEAD or pinned files changed.
The immutable bundle is `/var/tmp/wh2-k3-ordinary-recovery-r0`; native raw
SHA-256 is `1beaa881fb2a8ce1f13b58fd463271dc0f31b6dc7179cb1bddbce0bb6a370f86`.

This namespace is spent. The cohort is retained, not a fresh independent
holdout; source-policy and backend replays do not multiply its sample size.
The result is not a universal recovery guarantee, all-K construction-seed
validation or an old-path performance regression check.

### Preserved-path regression and small-state isolation

The pre/post shared-library screen at source `26304b3`,
`wirehair.wh2.admission-regression-cost-r0`, finished **REGRESSION**.
All 320 same-code controls passed across the two DSO load orders, but
57 and 56 treatment cells respectively showed resolved slowdowns. Certified
K2/K3/K4/K6 encoder and decoder lifecycles took 4.8-9.1% more time after the
ordinary K3 admission. The screen also retained order-sensitive K128 results
and changes in WH1 and opt-in paths whose producing objects were unchanged;
it does not establish allocation size as the sole cause. This is a native,
co-resident shared-library result, not a static-link or all-K claim.

The exact Python 3.8 replay and a separately written raw/statistical audit
reproduced the outcome before source advancement: 643 pins, 103,680 records
and 480 statistical decisions. The spent immutable bundle is
`/var/tmp/wh2-admission-regression-cost-r0`, with analysis SHA-256
`4bd1c590a36e677655524e04215959ef26fe3342ca01d4453b089d965ee15ada`.
These regressions remain an unmet ordinary-admission qualification requirement.

Subsequently, a `ctest -N` listing overwrote the pre-isolation build's
`Testing/Temporary/LastTest.log`, a non-producing diagnostic pinned by an older
receipt. The library objects and sealed measurement bundles are unchanged,
but that historical receipt's full artifact closure is no longer intact.
The follow-on screen records the lost log explicitly and proves all producing
inputs separately; the original strict verifier still rejects the mismatch.

The follow-on implementation moves K3-only owning pointers to a private
derived handle and keeps its operation bodies out of the common dispatch
functions. The immutable type tag uses existing header padding on the tested
ABI; no virtual dispatch, new allocation, equation or source-policy change
is introduced. Certified handle allocation returns from 296 to 272 bytes on
the native host. Allocation/OOM tests cover correct concrete-type destruction
and unchanged small-profile allocation counts. These are structural and
correctness checks, **not evidence that the timing regressions are fixed**;
performance qualification remains unmet, including preservation of the
ordinary K3 speed gain.

At source `6245943`, the follow-on
`wirehair.wh2.small-isolation-preserved-cost-r1` completed **CONTROL_FAIL**.
Both native workers completed normally with empty stderr and no codec or
observer errors. Two of 320 same-code controls failed the frozen equivalence
bounds, one in each library load order. Both were old-library opt-in K3
two-byte decoders; both confidence intervals included one. The complete run
therefore provides no qualified speed-restoration or regression result, even
though individual treatment comparisons have directional estimates.

The exact Python 3.8 replay and a separate standard-library audit reproduced
all 103,680 records and 480 statistical decisions and verified all 721 current
receipt pins before source advancement. The historical diagnostic-log loss
above remained explicitly disclosed and checked, not silently repaired.
The spent immutable bundle is `/var/tmp/wh2-small-isolation-preserved-cost-r1`,
analysis SHA-256
`a0970624f8b624498048f1af7c9c5647d3d58d96a9f5c015509970aa68b04933`.
No samples were discarded, pooled to rescue a failed control, or remeasured.
The preceding isolation R0 namespace is separately spent **INVALID**: its
worker read a stale claim-file path and stopped before any codec work. R1
fixed and positively tested that launch binding without changing the workload.

### Uncovered-K recovery inventory

At source `450eb82`, `wirehair.wh2.uncovered-k-recovery-inventory-r0`
compared actual ordinary source-independent WH2 with owned-input WH1 in the
unchanged native shared library produced at `5e902df`. This baseline-only
diagnostic completed; it tested no new candidate and measured no performance.
The frozen roster has 64 roots across three block widths (2, 64, 1280 bytes)
and four loss schedules (IID, burst, adversarial, repair-only) for each K.

Observed failures after receiving exactly K packets, with 768 traces per K:

| K | Ordinary WH2 | WH1 | WH2 minus WH1 |
|---|---:|---:|---:|
| 2 | 49 (6.38%) | 5 (0.65%) | 44 |
| 4 | 13 (1.69%) | 16 (2.08%) | -3 |
| 5 | 61 (7.94%) | 13 (1.69%) | 48 |
| 8 | 38 (4.95%) | 4 (0.52%) | 34 |

The prospectively fixed priority rule selects K5, then K2, K8 and K4. This
ordering is a diagnostic work priority, not a statistical superiority claim.
K5's 64-byte burst cell had 24/64 WH2 failures versus 3/64 for WH1; its
64-byte repair-only cell had 21/64 versus zero. After one additional packet,
K5 still had five WH2 failures and K8 had two; both had none after two.
Across all four K values WH1 had seven failures after one extra packet and
two after two, all at K4, with none after three.

The 24 separate full/partial low-repair hard cases are not in the rate
denominator. WH2 needed one extra packet in two K2 cases; WH1 needed none.
In every case, the first actual decoder success matched the first full-rank
source-equation prefix: no decoder lag after full rank was observed. Thus
the observed recovery gaps call for equation-structure work, not merely
earlier recognition of sufficient rank. This does not identify the internal
construction mechanism responsible for each deficiency.

Native basis-message probes observed each generator's coefficient rows.
Independent GF(256) arithmetic reconstructed every real packet and checked
rank; these are observed native rows, not an independent derivation of the
native generator. The worker checked source independence, packet guards and
two exact recoveries, freeing each encoder before creating its receiver.
Python 3.8 and 3.12 replayed all 683 receipt pins and the complete outcome.
A separately written audit with no project imports reproduced all 3,096
paired cases, 54,180 real-message packets, cell counts, priority decisions
and API ledgers before source advancement. The explicit historical diagnostic
log loss described above remained disclosed and checked.

The immutable bundle is `/var/tmp/wh2-uncovered-k-recovery-inventory-r0`;
raw SHA-256 is
`dc4bbbac93cad98669f325755c42f47339eed349d3064b0cca0908819c4c5729`.
Its namespace is spent. This inventory is not a candidate holdout, speed gate,
all-K construction validation or promotion. K5 is the next structural-screen
target; a survivor still needs independent recovery and actual full encoder
and decoder WH1 speed qualification.

### K5 structural-screen survivor (not an installed profile)

At source `13988b9`, the separately frozen
`wirehair.wh2.k5-thue-morse-r0` mathematical screen passed. The first local
algebra parameter, lambda 1, selected the GF(256) companion feedbacks
`(121,110,207,198,31)` and `(120,110,207,198,31)`. Selection did not use loss
traces or retry after observing recovery results. The dimension-five packed
lookup has 29,440 bytes, SHA-256
`4ac8059aba3b5797c8789c4258a1bda52e5fdb005592d601591705940cfe76c9`.

The selected pair passed all 1,260 local five-column minors, all 3,780
dyadic-seam/window minors, all 72 main-contract hard traces at zero overhead,
and all 54 distinct retained failure prefixes from 79 WH2/WH1 inventory
origins. Prefixes retained their original lengths; no extra packets were
appended to rescue a historical failure.

On 6,144 fresh loss traces, source-equation rank was deficient after exactly
five packets in 11 cases (0.1790%). Every case reached full rank with one
additional packet. All twelve 512-trace width/schedule cells met the frozen
1% bound; the largest count was four failures (0.78125%). These are rank
results on a different cohort from the baseline inventory, not a paired WH1
recovery comparison, native payload result or universal guarantee.

Independent Python 3.8 and 3.12 audits reproduced the selected parameter,
all 2,270 recorded packet rows, the entire lookup, all local/seam minors,
31,080 hard/fresh prefix ranks, history origins, and every cell decision
before source advancement. The audit used separate field arithmetic,
fraction-free elimination, recursive prefix products and trace generation.
All eleven source pins, interpreter bytes and immutable historical inputs
were checked; the exact receipt also replayed without rerunning the worker.

The spent immutable bundle is `/var/tmp/wh2-k5-thue-morse-r0`, raw SHA-256
`acddcd2d6c8dbfa7284aa900980b2845a03130b5f2f721515ac8d45bb580af76`.
This licenses the next bounded native-core/payload validation and actual
full encoder/decoder timing comparison. It does not add an installed K5
profile, change ordinary WH2 selection, repair the existing admission timing
regressions, or establish all-K performance or construction-seed coverage.

### K5 native payload qualification (benchmark boundary only)

The sealed K5 equations now instantiate the existing compile-time small core
and benchmark C boundary, using `WHK5` / `0x5748324b35544d31`. The decoder
algorithm and shared arithmetic are unchanged; K5 does not use the six-source
tiny-payload specialization. The generated lookup exactly matches the sealed
29,440-byte hash above. No installed K5 profile or ordinary selector is added.

Native, portable-arithmetic, and ASan/UBSan builds each passed all seven
standalone correctness tests. The retained corpus has 10,053 cases and 75,134
packet/feed/prefix-recovery checks, plus 2,270 recorded rows. It includes all
6,144 fresh and 72 hard traces, 57 original-width historical cases from the
54 unique prefixes, and every five-of-nine subset of the 30 seam windows.
Both independent and borrowed serialized constructors reproduce the same
packets and recovery; every encoder is freed before any receiver is created.
The independent polynomial oracle checks packet bytes, rank, no-write errors,
guards and repeated recovery. The observed prefix statuses match the original
eleven five-packet deficiencies and full recovery after one additional packet.

Separate neutral cases exercise partial tails, every packed selector, ownership
and detach, allocation failures, aliases, malformed and cross-dimension
descriptors, permanent public conflict poison, and original K6 parity.
Recompiling both affected native production translation units with their
recorded flags produces byte-identical `WirehairSmall.cpp.o` and
`WirehairV2Profile.cpp.o`. Existing production selection is unchanged.

This is an engineering replay of retained evidence, not additional independent
recovery samples or a paired WH1 recovery comparison. Qualified full-lifecycle
WH1 timing, installed admission, preserved-path regression qualification, and
the full all-K objective remain outstanding. See the standalone
[build instructions](bench/Wh2SmallNative/README.md).

### K5 serialized lifecycle screen: timing not qualified

At source `14d14e4`, `wirehair.wh2.k5-serialized-cost-r0` completed
**CONTROL_FAIL**. The exact qualified K5 benchmark boundary and actual WH1
and explicitly selected certified-WH2 libraries were compared separately for
independent and borrowed input. Full five-block messages used 2-, 64-, and
1,280-byte blocks. Encoder work includes create, descriptor output, 18 packets
and free; decoder work includes create, feed through its own first success,
recover and free, for separate low-ID and distant-ID streams.

All 38,880 callbacks completed with correct packets, descriptors, recovery,
guards and API ledgers: 4,976,640 fresh codec lifecycles and 59,719,680 API
calls. Every arm needed five packets in these particular decoder streams.
Native, portable-arithmetic and ASan/UBSan preflights reproduced all six arms'
fixtures, with every packet checked using independent GF(256) arithmetic and
first success checked against independently computed source-equation rank.

Three of 108 same-code controls failed the frozen equivalence bounds. All
three compare the borrowed K5 path with itself: two-byte encoder order 0 and
1,280-byte distant decoder orders 0 and 1. Their 95% ratio intervals were
`[0.96493, 1.01463]`, `[0.98785, 1.03696]` and `[0.98482, 1.03627]`.
Each includes one: these are failed precision/equivalence checks, not resolved
directional bias. Each affected cell contains one approximately threefold
CPU-and-wall-time excursion, with no recorded fault or context switch and
the same handle address as surrounding callbacks. The cause is unproven.
All 72 treatment constraints nominally passed, but cannot override the global
control failure or establish a K5 speed win.

Exact Python 3.8 and 3.12 replays and a separate standard-library audit agree
on all 180 statistical decisions, the complete chronology and ledgers, and
all 691 receipt pins before source advancement. The controller completed in
69.263 seconds, including 22.452 seconds of timed work, with empty stderr.
The immutable bundle is `/var/tmp/wh2-k5-serialized-cost-r0`, raw SHA-256
`59c5e626f8c6dffd6176e932d0753d2e5ee67d9ccd3daf54be6c8d53ae527019`.
This namespace is spent: no rerun, sample removal, resizing or subset rescue.
No installed/default, recovery-rate or all-K claim follows from this timing
screen. The separate paired recovery comparison is reported below; another
speed screen needs new evidence addressing the measurement failure, not a blind retry.

### K5 timing-excursion diagnosis: millisecond-boundary association

A read-only audit at source `b3037aa` inspected all 38,880 retained timing
records, including warm positions, without executing a codec or changing the
spent gate. Descriptive medians keep width, metric, comparison, order, logical
side and warm/measured position separate: 720 classes. Multiples of these
medians label observations for inspection; they are not new acceptance bounds
or rules for excluding samples.

Thirty callbacks exceed 1.5 times their respective wall-time median. Twenty-two
record an involuntary context switch. The remaining eight exceed 1.5 times
both their wall and thread-CPU-envelope medians with all four recorded counter
deltas zero:

| Path | Retained callback indices |
|---|---|
| K5 two-byte borrowed encoder | 4321 (warm), 5402 |
| WH1 two-byte borrowed encoder | 5431 |
| K5 1280-byte low-ID decoder | 2600, 13933 (warm), 14960 |
| K5 1280-byte distant-ID decoder | 9272, 19068 |

All eight start at absolute monotonic-clock phases 950,192–997,008 ns modulo
one millisecond and end 93,348–237,082 ns into the next millisecond. This
post-hoc temporal association uses the gate's existing one-millisecond phase
period; it is not a causal finding, a significance test, or permission to
avoid those phases in later speed measurements.

There are also eighteen codec-free waiting intervals exceeding their target
by more than 50 microseconds with no intervening recorded counter change.
Their overshoots span 56,964–272,985 ns, while wait wall time exceeds its
thread-CPU interval by only 130–432 ns. The waiting loop allocates no codec
and runs no codec operation. Thus similar charged-CPU delays do not require
K5 work, although these observations do not prove that every WORK excursion
has the same cause. Clock-read latency, interrupted computation, and internal
allocation effects are not separately observed by the existing envelope.

Source and exact measured-executable inspection also confirm that the recorded
outer handle omits a separate 88-byte K5 Encoder allocation, or a 136-byte
Decoder plus its six-block slab. The full-message borrowed encoder has two
allocations; a decoder has three. Equality of outer addresses cannot establish
equality of internal storage. Both decoder source-policy labels call the same
decoder. The thread-CPU envelope includes monotonic reads excluded from the
wall WORK interval, so a small positive CPU-minus-wall difference is expected
and is not evidence of extra codec instructions.

Python 3.8 and 3.12 produce byte-identical complete descriptive reports. A
separate raw traversal, nested-roster reconstruction and median calculation
reproduce every event index, all thirty wall excursions, all eighteen waiting
events, complete counter totals and the 22.452193143-second WORK ledger. Original
bundle hashes are unchanged before and after analysis. The audit is retained
under `/tmp/wh2-k5-excursion-audit.RlGLbm`; `report.json` SHA-256 is
`383750e4487173d0039dd5ceb7d79c4947330b5231449deef3a4524a6db8af15`.

The resulting next question is clock-read stalls versus interrupted fixed-count
computation, first in a separately frozen codec-free diagnostic with all phases
and observations retained. No allocator rewrite, speed qualification, or
production change follows from this audit. The unprivileged `perf` probe is
denied on this host; no host settings were changed. The latest Fable invocation
returned a usage-limit error and no report, so none of these findings is
attributed to Fable. K5 cost R0 remains **CONTROL_FAIL**.

### Codec-free diagnostic: slow computation after output boundaries

At source `3a2e099`, the separately frozen
`wirehair.wh2.clock-boundary-r0` completed its single codec-free run. All
524,288 observations and 34,359,738,368 fixed-count xorshift iterations passed
checksum, timestamp, CPU identity, affinity, ordinal and counter checks. All
48 passive absolute-millisecond bins are populated (9,741–11,063 observations
each). No phase targeting, warm exclusion, filtering or codec operation occurs.
See the [frozen protocol](bench/Wh2ClockBoundaryR0.md).

Four ordered TSC/AUX stamps separate first-clock, register-only computation
and second-clock windows. Descriptive times use the frozen full-run endpoint
conversion of approximately 3.195199 TSC ticks per nanosecond:

| Window | Median (microseconds) | Maximum (microseconds) |
|---|---:|---:|
| First monotonic read and capture bookkeeping | 0.0501 | 12.509 |
| Fixed-count register computation | 78.167 | 165.809 |
| Second monotonic read and capture bookkeeping | 0.0501 | 8.533 |

None of the 1,048,576 clock windows reached the frozen 50-microsecond label.
Fourteen computation windows exceeded 1.5 times their overall median, all with
zero recorded in-window fault or context-switch changes. They form just two
clusters: records 252672–252676 and 505344–505352, immediately after output
chunks. The preceding gaps were 9.694065 and 9.513142 milliseconds, but only
0.264537 and 0.165660 milliseconds of thread CPU time; each gap records one
voluntary context switch. These are the only voluntary switches and the only
greater-than-one-millisecond inter-observation gaps in the run. Computation
returns to its usual duration after five and nine observations respectively.

Thus long delays occur in fixed-count computation, not just in monotonic-clock
reads, and timed-window counters alone miss the preceding descheduling event.
The association with synchronous result output supplies a concrete reason to
remove interleaved file I/O from a future timing instrument. It does not prove
the precise blocking syscall, a CPU-frequency cause, or the cause of every
earlier K5 excursion. This diagnostic did not reproduce a 50-microsecond clock
window; that is not a universal bound on clock-read latency. The last four
absolute-phase bins also retain 2.08–2.24% higher computation medians, so output
changes cannot be assumed to remove all phase-dependent measurement effects.

Strict native and ASan/UBSan neutral checks, final-record failure retention,
claim authentication, sanitizer science-mode rejection and eleven Python
reader tests passed before the sole run. Both Python 3.8 and 3.12 replay the
full result exactly. A separate audit derives the checksum by binary-operator
exponentiation, traverses every raw record, and reproduces all windows, events,
phase summaries and receipt pins. Worker wall/CPU times were 41.718906 and
41.699199 seconds; controller plus analysis took 45.610457 seconds. Stderr is
empty. The immutable bundle is `/var/tmp/wh2-clock-boundary-r0`, raw SHA-256
`2584d2ef6eb3e0fa7702cd5adc0d038155acc71885fa9f64d55a1c94aff4727e`.

This namespace is spent. The next instrumentation change is deferred result
output using preallocated record storage, with complete failure retention and
all declared phases preserved. It requires neutral qualification before a new
ownership-matched K5 lifecycle gate. No old speed result is rescored, no clock
or wake-up cost is subtracted, and no K5 speed/admission or production change
is claimed here. The full all-K objective remains open.

### K5 paired retained recovery comparison

At source `3b82f35`, `wirehair.wh2.k5-serialized-recovery-r0` passed its
actual-codec retained-cohort comparison. The benchmark K5 boundary, actual WH1
and explicitly selected certified WH2 were tested separately with independent
and borrowed input. Failures after receiving exactly five packets:

| Codec, either source policy | Failures / retained traces | Observed failure rate |
|---|---:|---:|
| K5 candidate | 11 / 6144 | 0.179% |
| WH1 | 122 / 6144 | 1.986% |
| Certified WH2 | 439 / 6144 | 7.145% |

K5 has about 91.0% fewer zero-extra-packet failures than WH1 on this cohort.
It fixes all 122 WH1 failures but introduces 11 different failures. Against
certified WH2 it fixes 438 and introduces ten; one failure is shared. Thus this
is an aggregate recovery improvement, not per-trace dominance. All eleven K5
failures recover with one extra packet. WH1 still has two failures after one
extra packet and certified WH2 has ten; both have none after two.

All twelve width/schedule cells retain 512 traces and meet the candidate's
frozen 1% bound; the worst K5 cell has four failures (0.78125%). K5 is not
better in every individual cell: the 64-byte adversarial cell has two failures
versus WH1's one, and the 1,280-byte adversarial cell has four versus three.
These counts do not establish population-level or universal guarantees.

All 72 separate hard cases recover without extra packets for K5 and WH1,
versus nine certified-WH2 failures. K5 recovers all 57 historical original-width
cases from 54 distinct prefixes; WH1 remains unresolved on twelve and certified
WH2 on 45. Each historical prefix retains its original length; none is extended
to rescue an unresolved decode. Hard and historical cases are outside the
6,144-trace rate denominator.

Native, portable-arithmetic and ASan/UBSan runs agree on all 6,273 cases and
six API routes. Each backend uses 263,466 fresh handles, checks 337,404 real
packets plus their five basis probes, and performs 75,048 guarded successful
recoveries. Independent polynomial arithmetic verifies every arm's packet
bytes and first-success source-equation rank; candidate rows are also derived
from the fixed companion pair. All real and basis encoders are freed before
any receiver, and independent inputs are destroyed before encoding.

Exact Python 3.8 and 3.12 replays and a separate standard-library arithmetic,
chronology, ledger and paired-count audit agree, including all 657 receipt pins
before source advancement. The original qualification document is retained
at its exact Git version and the current report is pinned separately; producing
code and archives are unchanged. The immutable bundle is
`/var/tmp/wh2-k5-serialized-recovery-r0`, native raw SHA-256
`d6951f2594434c739776d0c37e1e4ed3b44317d91a49819dc832c28fee837303`.

This is a retained cohort, not a new independent holdout. Source-policy and
backend replays do not multiply its sample size. The namespace is spent. The
original K5 timing result remains CONTROL_FAIL; a separately qualified new
instrument provides the scoped speed result below. Installed admission,
preserved-path regressions and the all-K objective remain unresolved.

### K5 deferred-output lifecycle speed: qualified benchmark boundary

At source `4046600`, `wirehair.wh2.k5-deferred-cost-r0` completed **PASS**.
All 108 same-code equivalence controls and all 72 candidate comparisons pass,
separately for each width, lifecycle, source policy and measurement order.
Every candidate comparison's upper 95% time-ratio bound is below one against
its ownership-matched actual WH1 or explicitly selected certified-WH2 control.

Elapsed-time reductions versus WH1, spanning both source policies and orders:

| Block bytes | Full encoder | Low-ID decoder | Distant-ID decoder |
|---|---:|---:|---:|
| 2 | 67.31–67.65% | 91.85–91.95% | 85.45–85.50% |
| 64 | 76.12–76.55% | 91.39–91.42% | 86.03–86.20% |
| 1280 | 73.57–75.24% | 90.18–90.24% | 88.04–88.24% |

This is the actual benchmark K5 C boundary, not an installed ordinary-WH2 K5
profile. Full five-block messages use 2-, 64-, and 1,280-byte blocks on the
frozen native GFNI-capable CPU50 host. Encoder time includes create, descriptor
output, 18 packets and free; decoder time includes create, feed through its own
first success, recover and free. Both decoder streams require five packets
for every arm. These measurements do not establish partial-tail, other-width,
shared-library, cold-start, non-GFNI or all-K speed.

The new instrument defers all result formatting/output until measurement ends,
using the existing full static record array. The sealed K5 core, qualified
archives, WORK body, 38,880-callback roster, batch 128, 48 phase offsets and all
acceptance thresholds are unchanged. No new logging queue, thread or allocation
is introduced. No observations are dropped or rescored, no clock/wake-up cost
is subtracted, and the original spent CONTROL_FAIL remains unchanged. The new
PASS does not prove a unique cause for all old timing excursions.

All 4,976,640 fresh lifecycles and 59,719,680 attempted API calls passed their
packet, descriptor, recovery, buffer-guard and ownership checks. Timed WORK
totaled 18.361103009 seconds; the complete controller took 68.680927 seconds
and worker stderr was empty. All 249 involuntary context switches and the one
minor fault observed in/between captures remain in the data; no voluntary
switch was recorded there.

Native, portable-arithmetic and ASan/UBSan neutral qualification covered 108
WORK cells per backend, the full roster, late errors/exceptions, partial final
clock capture, source-corruption detection and output-device failures.
Twenty-eight reader tests pass under each Python 3.8 and 3.12. Exact replays
and an independently written standard-library packet/rank, chronology, ledger,
phase-coverage and statistical audit agree on all 38,880 records, 180 decisions
and 783 receipt pins before HEAD or pinned documentation advances.

The immutable bundle is `/var/tmp/wh2-k5-deferred-cost-r0`; raw SHA-256 is
`bfd61a30307f11bd28101dabdcd0249b5d751ad2ba35a3dd78d57a79a74dc18b`,
COMPLETE SHA-256
`2c10634555ad447b1c2264c854fd170b57550fe66426b4b6874a6d910b3e804e`.
See the [instrument contract](bench/Wh2K5DeferredCostR0.md). This namespace is
spent. The next requirement is installed/ordinary K5 admission with its own
actual-path speed, retained recovery and preserved-path regression gates;
benchmark timings cannot be inherited by a different owning public facade.

## Small K5 equation profile: explicit admission candidate

`WIREHAIR_V2_PROFILE_SMALL_K5_2026_09` has ID `80070c81bfe375f1`, the first
64 bits of SHA-256 over this exact name (no newline):

```text
wirehair:v2:small-k5:thue-morse-121-110-207-198-31:gf256-14d:2026-09
```

The full digest is
`80070c81bfe375f17b654786a037aa26fa1a3630417e0a3d6a5562199d8fc355`.
It freezes the sealed K5 GF(256) equations described above: polynomial `0x14d`,
companion feedbacks `(121,110,207,198,31)` / `(120,110,207,198,31)`,
right-product Thue-Morse packet-ID map, identity systematic rows, and zero tail
padding. The exact 29,440-byte lookup has SHA-256
`4ac8059aba3b5797c8789c4258a1bda52e5fdb005592d601591705940cfe76c9`.
Its packed layout and generic evaluator are the same as the qualified benchmark
core; a separate production translation unit owns the table, with no build-time
dependency on benchmark artifacts or Python.

The profile requires `4 * block_bytes < message_bytes <= 5 * block_bytes`,
`1 <= block_bytes <= 44739242`, and seed attempt zero. The bound keeps the six
decoder slab blocks at or below 256 MiB. Invalid shapes return
`InvalidDimensions`; nonzero attempts on valid shapes return `BadSeed`.
There is no construction seed search.

Use `wirehair_v2_encoder_create_profile_id()` (or its options/C++ counterpart)
to select it, or recreate from its serialized WHV2 descriptor. The ordinary
constructors still select the certified profile at K5; **defaults are unchanged**.
The benchmark WHK5 descriptor/ID is not accepted as a WHV2 descriptor and is
not reinterpreted as this profile.

K3 and K5 reuse the typed small-core ownership and dispatch code. Both storage
policies prepare one owned source basis; repair packets do not read borrowed
input and detach is allocation-free. An immutable byte-sized dimension tag
occupies existing common-header padding. Only small handles contain the three
small-state pointers; certified handles do not grow. Dependent contradictions
return `Error` without poisoning the retained decoder basis, matching ordinary
K3 semantics, not the separate benchmark facade's permanent-poison policy.

This is explicit public-path integration, **not performance promotion**.
Benchmark speed/recovery results above cannot qualify the different owning
facade. The installed retained-recovery gate below now passes; valid
actual-public-path speed and preserved-path regression/K3-retention gates
are still required before default selection.

### Public K5 integration: neutral correctness checks

Strict GCC 13.3 checks pass with ordinary compiler defaults (13 selected tests)
and host tuning enabled (29 selected tests). Portable arithmetic, ASan/UBSan
with leak/fake-stack checks, and Clang 18.1.3 each pass five selected tests.
Static Debug, shared Release, and dual RelWithDebInfo relocated-package checks
all pass, including installed C and C++ K5 consumers. The literal C receiver
uses independently derived repair packets before any encoder or explicit
runtime initialization exists.

The shared K3/K5 API test covers 69 width/tail shapes per dimension, nine K3
and six explicit K5 constructor/policy routes, standalone systematic/low/distant
decoders, and 4,864 K3 / 6,400 K5 packed-selector and mixed-chunk coefficient
oracle packets. Every constructor/decoder allocation is fault-injected;
packet operations and detach allocate nothing. Checks include protected
borrowed-source repair/detach, source release, short buffers, aliases, invalid
packets, conflicts before/after full rank and recovery, repeated recovery,
C++ transactional moves, all 255 nonzero seed-attempt rejections, and the exact
block-size bounds without large allocations.

With matching default compiler flags, 16 of the 17 existing archive members
are byte-identical to the preserved small-isolation library; only
`WirehairV2Profile.cpp.o` changes, plus the new K5 table object. Certified
handles remain 272 bytes and small outer handles remain 296 bytes on this host.
The four common encode/decode/recover/free function sizes remain unchanged;
that is layout evidence, not a performance result. The certified compatibility
replay matches all 2,180,292 output bytes across 48 cases, SHA-256
`2e6536dcd86a7c2892399ddf1f14c3ff2290c2ef9e270aaa0c5ed87d0928907b`.
These are engineering checks, not a new loss cohort, comparative recovery
qualification, or evidence that the earlier performance regression is fixed.

### Installed K5 lifecycle speed: controls remain unqualified

At source `1cc38ba`, `wirehair.wh2.k5-public-cost-r0` completed
**CONTROL_FAIL**. The explicit installed WHV2 K5 profile, certified WH2 and
WH1 used one qualified default-flags static archive, both ownership policies,
full K5 messages, three block widths and full encoder/low-ID/distant-ID
decoder lifecycles. The worker completed normally with empty stderr and no
codec, oracle, guard, ledger or publication error.

Three of 108 same-code controls failed the frozen equivalence bounds. All
three concern 1,280-byte K5 decoders, and every interval includes one:

| Same-code K5 control | Order | Estimated time ratio | 95% interval |
|---|---:|---:|---:|
| Independent, low-ID decoder | 1 | 1.00500 | 0.98600–1.02438 |
| Independent, distant-ID decoder | 0 | 0.99273 | 0.97753–1.00816 |
| Borrowed, distant-ID decoder | 1 | 1.00640 | 0.98860–1.02452 |

These are failed precision/equivalence checks, not established directional
bias. All 72 candidate comparisons nominally meet their speed threshold,
but cannot override the failed controls. This run supplies **no qualified
installed K5 speed or default-promotion result**. No samples were discarded,
orders pooled, pauses subtracted or measurements repeated.

All 38,880 records, 4,976,640 fresh lifecycles and 59,719,680 attempted API
calls remain retained. Timed WORK was 22.120847526 seconds; the complete
controller took 69.612015 seconds. Native, portable-arithmetic and ASan/UBSan
neutral qualification passed all 108 WORK cells, the full roster, independent
packet/rank fixtures, five deferred-publication streams and output-device
failures. All 32 reader tests pass under Python 3.8 and 3.12. Both exact
replays and a separately written standard-library matrix/payload/rank,
chronology, ledger, phase, statistical and provenance audit reproduced all
180 decisions and 735 receipt pins before HEAD or pinned documentation
advanced. The independent audit also agrees under both Python versions.

Retained-raw inspection finds an involuntary context switch in the largest
pair excursion of two failed replicates; the third has increased thread CPU
work without a recorded switch. Across all captures there are 1,699
involuntary switches and 2,247,239 minor faults, with no major faults or
voluntary switches. All minor faults occur in wide certified-WH2/WH1 decoder
callbacks, not K5 callbacks. This does not establish why K5 controls failed.

A separate bounded syscall trace of the existing neutral worker passes all
108 WORK cells and records 8,228 `brk` calls, including repeated heap growth
and shrinkage. Its 22 `mmap` calls and single `munmap` belong to startup;
there is no repeated mapping/unmapping or `madvise` stream. Thus heap
growth/trimming is an observed allocation mechanism in this neutral workload.
A subsequent bounded stack trace attributes WH1 growth to `AllocateWorkspace`
and shrinkage to `FreeWorkspace`/`FreeInput`; certified WH2 grows at
`quotient_pivot_rhs` allocation and trims during RHS-buffer destruction and
outer teardown. This does not establish the scientific control failures' cause.
The earlier five-buffer release-order experiment already measured small encoder
regressions; a duplicate variant was stopped before implementation or testing.

The immutable, permanently spent bundle is
`/var/tmp/wh2-k5-public-cost-r0`, raw SHA-256
`7537c8aa1610c8d6df47875a4d5455a18553b84a9c496fd1c8e5fb39823768ac`,
COMPLETE SHA-256
`113d5d9630c6a4f8bf7548ed2b0f15fc41ef9e81fd858a711c22839ebbd1f11a`.
Qualification and independent audit artifacts are retained at
`/tmp/wh2-k5-public-cost-neutral.IR0xrEBi`.
See the [installed gate contract](bench/Wh2K5PublicCostR0.md).
Valid installed speed, preserved-path regressions and ordinary K3 speed
retention remain required before default promotion.

### Installed K5 retained recovery: qualified

At source `2a873f9`, `wirehair.wh2.k5-public-recovery-r0` passes through the
actual installed WHV2 APIs, using the explicit K5 and certified profiles plus
ownership-matched WH1. Native, portable arithmetic and ASan/UBSan agree on all
6,273 cases and both source policies. The full records reproduce the earlier
prototype results with only the candidate's serialized descriptor changed.

On the same 6,144 retained loss traces:

| Codec | Zero-overhead failures | Rate | Failures with one extra packet |
|---|---:|---:|---:|
| Installed K5 | 11 | 0.1790% | 0 |
| WH1 | 122 | 1.9857% | 2 |
| Certified WH2 | 439 | 7.1452% | 10 |

K5 has about 91% fewer aggregate zero-overhead failures than WH1. It fixes all
122 WH1 failures but introduces 11 different failures; there is no per-trace
or per-cell dominance. Every 512-trace candidate cell is at most 1%, with a
worst cell of 4/512 (0.78125%). All 72 hard cases recover at zero overhead for
K5 and WH1; certified WH2 needs one extra packet in nine. All 57 original-width
historical cases recover for K5, versus 12 unresolved WH1 and 45 unresolved
certified cases at their original, unextended horizons.

All real and basis encoders are destroyed before receivers. Independent
polynomial coefficients, every-arm basis-derived packet hashes, rank/first
success, twice-recovered messages, guards, immutable inputs and complete API
ledgers pass. Each backend checks 337,404 real packets and 1,687,020 basis
packets. These are repeated engineering checks of retained traces, not new
independent samples, a population failure bound, speed evidence, or an all-K
guarantee. Defaults remain unchanged.

Final neutral qualification covers 36 cases and four late-call cleanup failures
per backend; all 54 recovery/shared-build tests pass under Python 3.8 and 3.12.
An initial sanitizer-driver context-layout mismatch was caught and fixed before
the retained run; compiler ISA flags and linked GF context sizes now match.
The sole controller completed in 32.703651 seconds with all workers successful
and empty stderr. Exact replays and a separately written arithmetic/raw/build
provenance audit agree under both Python versions, checking all 739 receipt
pins before HEAD or pinned documentation advanced.

The permanently spent bundle is `/var/tmp/wh2-k5-public-recovery-r0`, COMPLETE
SHA-256 `66986b694e8d9a75eb19d913f5c18a580aaa38aaf341dbc85660376b7f798644`,
native raw SHA-256
`4fdc142519395b79de905be46c0fc336c1039b48fffe754c999412f476306fda`.
Qualification and audits are retained at
`/tmp/wh2-k5-public-recovery-qualified.bksSezzu`.
See the [frozen recovery contract](bench/Wh2K5PublicRecoveryR0.md).

### Current preserved paths: valid regression result

At source `bdf5c83`, `wirehair.wh2.current-preserved-deferred-cost-r0`
compares the exact pre-admission shared library against the current
K5-integrated native-default library. It retains the original 20-case
preserved-path roster, workloads, batches, 48 delay phases and decision rules,
but publishes the header and every retained record after measurement ends.
No codec implementation or equation changes for this screen.

All 320 same-code controls pass. Both load orders nevertheless report
**REGRESSION**: 42/80 and 49/80 treatment cells have resolved slowdowns;
15/80 and 28/80 have upper confidence bounds reaching 2% (these sets overlap).
Certified K2/3/4/6 full-encoder point estimates are 0.30–3.09% slower than
pre-admission WH2 across the declared widths and orders; not every individual
estimate is resolved. Opt-in K6 two-byte encoders improve by 4.39–4.91%, but
that does not satisfy the complete preserved-path gate. No pooling, trimming,
rerun, threshold change or promotion follows from this result.

The mechanism remains unresolved. WH1 two-byte decoder timings also regress
despite its unchanged codec object, and K128 results change direction with
measurement order. Fifteen of the original 17 archive objects are byte-identical;
only `WirehairSmall.cpp.o` and `WirehairV2Profile.cpp.o` differ, with the K5 table
object added. This motivates distinguishing facade overhead from linked-code
placement and physical-context effects before another candidate is measured.
It does not establish any one of those as the cause.

Both workers finish normally with empty stderr. All 103,680 callbacks, payload
and ledger checks are retained; WORK totals are 36.446777265 and 36.446111810
seconds, and the controller takes 230.228120 seconds. Native and ASan/UBSan
driver qualification pass, including late failures, incomplete clocks and
output errors; the loaded native DSOs are not sanitizer-instrumented. All 32
reader tests pass under Python 3.8 and 3.12. Exact replays and a separately
written full raw/statistical/library-provenance audit agree under both versions
before HEAD or this document advances: 770 receipt pins, 480 decisions and a
maximum independent arithmetic difference below `1e-15`.

The permanently spent bundle is
`/var/tmp/wh2-current-preserved-deferred-cost-r0`, COMPLETE SHA-256
`4b0e106c180aee5460d467c93171a903735b82654159ff8a462ddf2430384125`.
Qualification and audits are retained at
`/tmp/wh2-current-preserved-deferred-qualified.B6qTTW13`.
See the [frozen gate contract](bench/Wh2CurrentPreservedDeferredCostR0.md).
Earlier admission REGRESSION, isolation CONTROL_FAIL and installed K5
CONTROL_FAIL outcomes remain unchanged. Preserved-path restoration, current
ordinary K3 retention and installed K5 speed are still unmet; this screen
does not qualify strict WH1 speed, recovery, static speed, defaults or all K.

## Small K8 equation profile: explicit admission candidate

`WIREHAIR_V2_PROFILE_SMALL_K8_2026_09` has ID `7a9276b85c730ae0`, the first
64 bits of SHA-256 over the exact name (no newline):

```text
wirehair:v2:small-k8:thue-morse-96-19-186-153-85-252-7-255:lambda-2:gf256-14d:2026-09
```

The full digest is
`7a9276b85c730ae04a8ac1c2648e556d1586f771d5f35b5db8ad1bbeead2ad95`.
This freezes GF(256), polynomial `0x14d`, and the selected Thue-Morse
companion pair `(96,19,186,153,85,252,7,255)` /
`(98,19,186,153,85,252,7,255)`. The phase perturbation is two, not the
one used by K3/K5. The right-product packet-ID mapping, systematic identity
rows and zero tail padding equal the qualified benchmark K8 equations.
Its 65,536-byte immutable lookup has byte SHA-256
`512c6646e44517964e7e6a7cd0ffa41057802182ddc500a818af541c05770817`.
The table is compiled into the ordinary library; installation and consumers
need neither Python, benchmark artifacts, a live encoder nor a seed search.

Require exactly eight source blocks, `7*block_bytes < message_bytes <=
8*block_bytes`, `block_bytes <= 29826161`, and attempt zero. The bound limits
the decoder's eight pivot blocks plus scratch to at most 256 MiB. Invalid
dimensions precede `BadSeed` for a nonzero attempt. ID7 carries only the
meaningful final-block bytes; all other IDs carry one full block.

This is an explicit WHV2 profile, not a new default. Ordinary K8 constructors
still select certified WH2. K3 selection, K5 explicit selection, the CURRENT
alias, all existing profile IDs and the opt-in WHK3/K6 APIs are unchanged.
WHK8 prototype descriptors and its distinct profile ID are rejected, never
silently reinterpreted as WHV2.

K8 reuses ordinary WHV2's small-profile owning state and shared core. Both
source policies prepare one owned source basis, with private zero padding for
a partial tail. Borrowing adds no full-message copy or allocation; repairs
do not read borrowed input, and detach is allocation-free. Contradictory
dependent packets report `Error` without poisoning the retained basis, even
after recovery; successful recovery can be repeated. This differs from the
prototype facade's borrowed-source and permanent-poison contracts.

The benchmark boundary has independently audited
[paired retained recovery](bench/Wh2K8SerializedRecoveryR0.md) and
[full-lifecycle speed](bench/Wh2K8SerializedCostR0.md) passes. Those results
do not qualify this installed owned-basis path. Its own installed recovery
and full-lifecycle speed gates below now pass; preserved-path restoration and
current K3 speed retention remain separate admission requirements.
No default, all-K, non-GFNI-host speed or all-construction-seed qualification
is implied by adding the explicit profile.

### Installed K8 neutral qualification

The explicit profile passes correctness and package qualification at
`/tmp/wh2-v2-k8-admission.YNN2XmAx`. GNU 13.3 native static/shared testing
passes 34 selected gates, including all three relocated package configurations
(static Debug, shared Release, and dual RelWithDebInfo). Portable arithmetic
passes six selected gates; ASan/UBSan with leak and fake-stack detection passes
six; Clang 18 shared-primary passes eight. Its C++ codec tests link
`libwirehair_test_support.a`; its C consumer and export checks exercise the
actual shared DSO. The separately retained `LAUNCH_ENVIRONMENT.json` records
the explicit sanitizer launch options; it is execution metadata, not part of
the earlier neutral inventory. Portable arithmetic here is a forced backend
on this host, not a separate non-GFNI-machine performance result.

The public K8 correctness test covers 69 width/tail shapes through six explicit
constructor/policy routes, independent lambda-2 packet oracles, every constructor
allocation failure, aliases, source release, protected borrowed-source repairs,
allocation-free detach, conflict/basis preservation, repeated recovery and C++
transactional ownership. An additional 96 rotated/reversed systematic-order
cases exercise pivot bit 7, full/one-byte tails and conflicts before and after
recovery. Bounds include the exact maximum block size, both message boundaries,
all 255 invalid nonzero attempts and prototype/retired descriptor rejection.
K3/K5 tests and ordinary K8 default-preservation checks pass too.

`bench/Wh2V2SmallProductionParity.cpp` additionally replays the unchanged sealed
native corpus through both the original benchmark boundary archive and the
actual WHV2 library, translating descriptors explicitly in test code. Per
native/portable/sanitizer backend, all 21,110 K8 corpus cases and 193,746 packet
oracles pass. The 2,347 retained coefficient-row checks call the reusable
core's `S::Row`, not the installed WHV2 API. Both source policies produce
42,220 paired external lifecycles and 387,492 packet comparisons. Every prefix
status and recovery agrees, including repeated recovery, original historical
widths/tails/horizons and all 14,850 seam subsets; encoders are freed before
receivers are created. K3's 7,774-case and K5's 10,053-case retained corpora also
pass through this same adapter on all three backends. These are semantic
replays, not new recovery samples, timing, or equality of the two APIs' distinct
ownership/detach/conflict contracts. Core-only subtests in the reused harness
are not counted as installed contract coverage.

The same certified-compatibility program linked against the pre-K8 and current
native libraries emits exactly 2,180,292 identical bytes over 48 cases; SHA-256
is `2e6536dcd86a7c2892399ddf1f14c3ff2290c2ef9e270aaa0c5ed87d0928907b`.
Seventeen of eighteen existing archive members are byte-identical: only
`WirehairV2Profile.cpp.o` changes, plus the new K8 lookup object. This is byte
compatibility and change isolation, not evidence of preserved-path speed.

Clang testing exposed an export-checker portability bug: LLVM nm decorates
the zero-valued absolute ABI version marker while GNU nm leaves it bare.
Both describe the same 53 function exports. The checker now normalizes only
either exact zero-valued/zero-sized absolute marker spelling, at most once;
19 synthetic parser fixtures retain strict rejection of malformed metadata,
duplicate markers, unexpected functions and wrong function versions. The
original failed Clang log is retained. Repeated main and independent complete
source reviews of the integration, parity adapter and checker fix found no
remaining defects.

Native archive SHA-256:
`0d4a0a7ec3cfa4b8e0f08df75e2c12935b51f28653c8c953a32bb49a5658eff5`;
native DSO:
`ef684abac606667e30e3b0de1204b0897aa0fa93a4f5ec2268d1ba951dce03b2`.
Actual installed paired recovery and full-lifecycle speed are tracked separately
in `.82.5.3.2` and `.82.5.3.3`. The existing preserved-path regressions and
current K3 speed retention remain unqualified; no default promotion follows
from these neutral results.

### Installed K8 retained recovery: qualified

At source `5a30050`, `wirehair.wh2.k8-public-recovery-r0` passes through the
actual installed WHV2 APIs, comparing explicit K8, ordinary certified WH2 and
WH1 under matched independent/owned and borrowed source policies. Native,
portable arithmetic and ASan/UBSan agree on all 6,260 retained cases and all
six arms. No benchmark facade is linked and candidate descriptors are the
actual WHV2 records, not translated prototype identities.

On the 6,144 retained loss traces, each source policy gives:

| Codec | Zero-overhead failures | Observed rate | Failures with one extra packet |
|---|---:|---:|---:|
| Installed explicit K8 | 12 | 0.1953125% | 0 |
| WH1 | 37 | 0.6022% | 5 |
| Ordinary WH2 | 246 | 4.0039% | 3 |

K8 has about 67.6% fewer aggregate zero-overhead failures than WH1 and 95.1%
fewer than ordinary WH2. It fixes 37 WH1 failures but introduces 12 different
ones; against ordinary WH2 it fixes 245 and introduces 11, with one shared
failure. There is no per-trace or per-cell dominance. All twelve 512-trace
candidate cells meet the frozen 1% threshold; the worst is 3/512 (0.5859375%).
All twelve candidate failures recover with one additional packet.

All 72 separate hard cases recover without extra packets for K8 and WH1,
versus six ordinary-WH2 failures. All 44 historical origins recover for K8 at
their original widths/tails and unextended horizons, versus five unresolved
WH1 and 41 unresolved ordinary-WH2 origins. These cases are outside the loss
denominator. The twelve cells share 512 roots; policy/backend replays do not
increase the independent sample size.

Independent coefficients, every-arm packet arithmetic and prefix ranks,
twice-recovered guarded messages and complete attempted API ledgers agree.
All senders are destroyed before receivers. Complete six-arm prototype parity
permits only separately validated candidate descriptor differences, without
rewriting either record stream. The different API ownership/detach/conflict
contracts retain their separate neutral qualification.

All nineteen production objects and full archives reproduce byte-for-byte on
all three backends. Forty-eight neutral cases and twelve late-call cleanup
checks pass per backend; all 27 controller/build tests pass on Python 3.12 and
3.8. Exact replays and a separately written full arithmetic, summary,
archive/ELF and 1,287-pin provenance audit pass on both interpreters before
source advancement. All sixteen auditor selftests and repeated source reviews
are clean. Worker stderr is empty; the complete controller took 66.611779
seconds. These are retained-data checks, not new recovery samples or timing.

The spent bundle is `/var/tmp/wh2-k8-public-recovery-r0`, COMPLETE SHA256
`699415eef17fcc486aeb3e5487feef0a9dcd4812110497df1370014c47b7af02`;
native raw SHA256
`bab96ca60d697f9f5a751fdbc74dce306c6c0d12a98e651f7e9735beb6af50e1`.
See the [installed recovery contract](bench/Wh2K8PublicRecoveryR0.md).
The separately qualified installed speed result follows below. Preserved-path
restoration and current K3 speed retention remain unmet; defaults are unchanged
and no all-K or default-promotion claim follows.

### Installed K8 full-lifecycle speed: qualified

At source `6f557dd`, `wirehair.wh2.k8-public-cost-r0` completed **PASS**
through the actual installed WHV2 APIs. All 216 same-code controls and all
144 candidate comparisons pass, separately for each shape, lifecycle, source
policy and measurement order. Every candidate comparison's upper 95%
time-ratio bound is below one against both ownership-matched actual WH1 and
ordinary certified WH2. There is no 5% minimum, filtering, timing subtraction,
order pooling or candidate regression tolerance.

Observed elapsed-time reductions versus WH1, spanning both source policies
and measurement orders (ranges of separate estimates, not pooled results):

| Block bytes | Tail bytes | Full encoder | Low-ID decoder | Distant-ID decoder |
|---|---|---|---|---|
| 2 | 2 | 53.61–54.72% | 86.29–86.44% | 66.15–66.42% |
| 2 | 1 | 52.62–53.56% | 86.35–86.51% | 66.36–66.60% |
| 64 | 64 | 64.72–65.24% | 85.71–85.81% | 68.82–69.01% |
| 64 | 1 | 63.48–63.78% | 85.67–85.73% | 68.86–69.01% |
| 1280 | 1280 | 63.90–65.27% | 82.73–82.92% | 77.31–77.42% |
| 1280 | 1 | 63.99–64.85% | 82.88–82.90% | 77.22–77.37% |

The candidate uses explicit profile `7a9276b85c730ae0` with independent and
borrowed options, matching the installed recovery routes. Both policies'
owned-basis allocation and preparation stay inside timed creation. Encoder
time includes create, descriptor output, 24 encodes and free. Decoder time
includes create, feeds through its own first success, recover and free; all
arms require eight packets on both frozen decoder streams. These are actual
installed static-library calls, not the prototype facade or translated
descriptors. Source policies share the receiver implementation.

All 77,760 callbacks, 9,953,280 fresh codec lifecycles and 159,252,480 attempted
API calls are retained with packet, descriptor, recovery, ownership and buffer
checks. All 160 involuntary switches and the one minor fault within WORK are
retained. WORK totals 38.577614860 seconds; the controller completes in
149.622703180 seconds with empty worker stderr.

All nineteen original production objects and complete archives reproduce
byte-for-byte in native, portable-arithmetic and sanitizer qualification.
The three backends agree on neutral packet/rank fixtures; only native is timed.
All 38 reader/build tests and all 16 independent-auditor selftests pass under
Python 3.12 and 3.8. Both exact retained replays and a separately written full
arithmetic, chronology, API-ledger, statistical and build-provenance audit
pass before HEAD or pinned documentation advances. Independent decisions
agree exactly across interpreters; the maximum numeric difference from the
controller is below `9e-16`. All 915 receipt pins are verified.

The permanently spent bundle is `/var/tmp/wh2-k8-public-cost-r0`, COMPLETE
SHA-256 `787021740f10034a3c7acba4ef033b990c4b3fe215f7692de240178f0953f942`;
native raw SHA-256
`c28f29059652ccc3cfce7d3b29615b4e15224018c3411818ad478e7947fe65fa`.
See the [installed timing contract](bench/Wh2K8PublicCostR0.md).
This qualifies only the six explicit-profile shapes and frozen lifecycles
after process initialization on this GFNI-capable host. It is not shared-call,
cold-start, non-GFNI-host or all-K performance, a new recovery sample, restored
preserved-path speed or current K3 speed retention. Defaults remain unchanged.

### Ordinary K8 selector: candidate lifecycle speed qualified

The isolated six-line ordinary-K8 selector now has its own independently
audited full-lifecycle speed result at source `7088ddf`, separate from the
explicit-profile timings above. All 216 same-code controls and all 144 strict
comparisons pass against ownership-matched WH1 and explicitly selected
certified WH2. Encoder time is 51.32–66.04% lower than WH1; low-ID decoder time
82.81–86.23% lower, and distant-ID decoder time 66.21–77.42% lower. These are
ranges of separate point estimates across 2/64/1280-byte blocks, full/one-byte
tails, both source policies and both measurement orders, not pooled results.

The actual ordinary entry points include owned-basis preparation, descriptor
output, 24 encoder packets or decoding through each arm's own first success,
recovery and teardown. All 77,760 observations are retained. Exact replays and
independent full raw/statistical/build audits pass under Python 3.12 and 3.8
before source advancement. See the [complete result and scope](bench/Wh2K8OrdinaryCost/README.md).

The selector remains benchmark-only: **ordinary production K8 still selects
certified WH2**. Its separate retained-recovery comparison below now passes;
pre-admission preserved-path restoration and current K3 speed retention remain
outstanding.
This is not shared-call, cold-start, non-GFNI-host or all-K speed, a new
recovery sample or zero-bad-construction-seed qualification.

### Ordinary K8 retained recovery: qualified

At source `adeb96c`, `wirehair.wh2.k8-ordinary-recovery-r0` passes through
the actual ordinary independent-bare and borrowed-options constructors in the
same six-line selector candidate. Each source policy gives these results on
the 6,144 retained loss traces:

| Codec | OH0 failures | Observed rate | Failures with one extra packet |
| --- | ---: | ---: | ---: |
| Ordinary K8 candidate | 12 | 0.1953125% | 0 |
| WH1 | 37 | 0.6022% | 5 |
| Explicit certified WH2 | 246 | 4.0039% | 3 |

K8 fixes 37 WH1 failures and introduces twelve different ones; it fixes 245
certified-WH2 failures and introduces eleven, sharing one failure. Every
candidate cell meets the frozen 1% bound; the worst is 3/512. All candidate
failures recover with one extra packet, all 72 hard cases recover at OH0 and
all 44 historical origins recover within their original, unextended horizons.
This is aggregate retained-cohort improvement, not per-trace dominance or a
population guarantee. Policy/backend replays do not multiply the shared roots.

Native, portable and full ASan/UBSan agree on all 6,260 records and 74,946
packet IDs. Every record exactly matches installed explicit-K8 evidence,
including the actual WHV2 descriptor; the separate prototype parity remains
unchanged. Independent payload/rank, first-success and guarded repeated
recovery checks pass. All nineteen original archive objects and all three
selector objects reproduce exactly; link maps exclude the old facade and
select one matching GF runtime.

Exact replays and independent full arithmetic/raw/build audits pass under
Python 3.12 and 3.8 before source or pinned documentation advances. All 1,346
receipt pins and every summary cell/fix/introduction list agree. The spent
bundle is `/var/tmp/wh2-k8-ordinary-recovery-r0`, COMPLETE SHA256
`22c461bb67fe7458d6aa36e5272fe96ec3deeb7f3cf8228485305336ae4ac992`.
See the [complete ordinary recovery result](bench/Wh2K8OrdinaryRecovery/README.md).

Both scoped ordinary K8 gates now pass. Production defaults remain unchanged:
pre-admission preserved-path restoration and current K3 speed retention still
precede promotion. Neither this retained replay nor the separate static speed
screen establishes all-K, shared-call, cold-start or non-GFNI-host performance.
The separate [native shared and relocated-package qualification](bench/Wh2K8OrdinaryShared/README.md)
now passes actual shared C/C++ ownership/provider and installed consumer/plugin
checks. Its [prospective shared gates](bench/Wh2K8OrdinaryShared/Gates.md)
separate current-path retention, ordinary K8 superiority and ownership-matched
K3 speed/retention. The independently qualified gate-A instrument has now
completed **CONTROL_FAIL**: current-then-candidate has two failed same-code
controls, while the reverse load passes every control and finds 17 resolved
slowdown cells. The exact replay and independent full audits agree. Current
shared-path retention is not qualified; gates B/C did not launch and the
selector remains unpromoted. See the [complete retained gate-A result](bench/Wh2K8SharedRetention/README.md#independently-audited-gate-a-result-not-retained).

The [read-only retained-data diagnosis](bench/Wh2K8RetentionDiagnosis/README.md)
and independent source/ELF review now cover all 196,992 observations and 912
cells. For valid full K3 B2/B64, the bare ordinary constructor skips the added
K8 checks and executes 48 non-NOP instructions in each library, counting its
call/return but excluding callees. The counts, not the instructions, are equal.
All 108 paired public-handle
vectors match in each of 32 investigated encoder cells; 31 cells have zero
capture fault/switch deltas, with one retained WH1 B1280 prelude switch in the
exception. Failed-control excursions resemble the previously observed
millisecond-boundary association. These findings do not establish a cause or
license sample exclusion, a gate rerun or another selector/layout variant.
The original verdict, stopped B/C gates and unpromoted selector are unchanged.

## K4 structural-screen survivor: not an installed profile

The separately frozen `wirehair.wh2.k4-thue-morse-r0` mathematical screen
passes at source `86acfff`. First local lambda 1 selects the GF(256)
companion feedbacks `(64,120,54,15)` and `(65,120,54,15)`, polynomial `0x14d`,
with the same dimension-parametric mapping approach, not inherited K3/K5/K8
qualification. The packed lookup has 20,480 bytes, SHA256
`f4c3bc7d9fffac727acefc217d0481206b4c2d608af051993f34fb6a19f47c69`.

All 700 local minors, 2,170 seam minors, 72 hard OH0 traces and 37 original-
length historical prefixes pass, preserving all 38 WH2/WH1 failure origins.
On 6,144 fresh loss traces, nine have deficient rank after exactly four
packets (0.146484375%); all become full rank with one extra packet. Every
512-trace width/schedule cell meets the frozen 1% bound; the worst has four
failures (0.78125%, B2 repair-only).

Exact and independent full arithmetic, lookup, trace/rank and provenance
audits pass under Python 3.12 and 3.8 before source advancement; the independent
reports are byte-identical. The unchanged native FrozenTrace generator agrees
on all 6,216 hard/fresh traces, without executing a codec. Thirty-four neutral
screen/controller tests and five independent-auditor tests pass per interpreter.

See the [complete K4 contract and result](bench/Wh2K4ThueMorseR0.md).
This is structural feasibility, not a paired WH1 recovery comparison or native
speed evidence. Ordinary K4 still selects certified WH2. Native payload and
ownership-matched full encoder/decoder qualification remain required, alongside
the unresolved existing-path admission and all-K objectives. The scientific
namespace is spent and must not be rerun or retuned.

## Retired equation profile identifiers

The identifiers `e161ce5d456f9bb7` and `20a4f27a870612a2` are permanently
reserved tombstones for rejected experimental contracts. They are unsupported,
have no public selector constants, and must never be assigned to another
equation system. For otherwise well-formed descriptors and API arguments,
serialization, parsing, validation, explicit encoder selection, and codec
creation return `WirehairV2_UnsupportedProfile` for either identifier without
writing a descriptor or publishing a codec handle.

## APIs and errors

`wirehair_v2_encoder_create()` prepares owned equation state, chooses the
profile (and deterministic seed attempt where applicable), and returns the
serialized descriptor. A null or short descriptor
buffer returns `WirehairV2_BufferTooSmall`, reports the required 32 bytes, and
does not create a codec. Its descriptor-size output pointer is required.
`wirehair_v2_encoder_create_profile()` recreates an encoder under an existing
descriptor. In both forms, the message pointer must provide at least the exact
message byte count supplied directly or recorded in that descriptor; the
implementation copies those bytes before returning.

`wirehair_v2_encoder_create_profile_id()` performs the same operation for an
explicit supported profile ID. `WIREHAIR_V2_PROFILE_CURRENT` deliberately
remains an alias for `WIREHAIR_V2_PROFILE_CERTIFIED_2026_07`; explicitly selecting
either name continues to emit the original GF(256)-only equations byte-for-byte,
including at K3. The ordinary selector is a dispatch policy, not this constant.
The corresponding C++
`Encoder::Create(profileId, ...)` overload provides the same explicit
selection. Unknown or retired IDs return `WirehairV2_UnsupportedProfile`
without falling back to the current profile.

The three additive `*_with_options()` encoder constructors accept an exact
versioned `WirehairV2EncoderOptions` record. The default initializer selects
`WirehairV2EncoderSource_Independent`, which has the same source-independent
lifetime contract as the original constructors: after successful return the
caller may change or release its message. A zero-initialized record selects no
policy and is rejected, as are unknown policy values.

`WirehairV2EncoderSource_BorrowedImmutable` is an explicit opt-in local storage
policy. Construction still completes the same eager full solve before success,
but the encoder then retains the caller's exact message range without owning it
or adding a copy beyond the independent constructor's prepared state. The
caller must keep that range readable,
allocated, and byte-for-byte immutable from constructor entry throughout the
call. Failure retains nothing and ends that obligation; after success it
continues until successful `wirehair_v2_encoder_detach_input()` or encoder
destruction. While attached, systematic packet IDs copy only their meaningful
source bytes; repair IDs use the existing solved-intermediate evaluator and do
not read the source. Detach is allocation-free and idempotent, and later
systematic packets use equation evaluation while remaining byte-identical.
Operations on one codec are externally serialized: encode, detach, free, and
source mutation must not race, including with another codec operation.

Source-storage policy is not serialized, does not change the profile ID or any
descriptor byte, and is not visible to a decoder. Independent and borrowed
encoders therefore use the same selected GF(256)-only equations, selected seed
attempt, solved intermediate state, systematic bytes, and repair bytes. Unknown
or retired equation profiles remain rejected for every storage policy.

`wirehair_v2_encode()` reports `WirehairV2_BufferTooSmall` and the exact
required packet size without modifying a short non-null output buffer.

`wirehair_v2_decoder_create()` takes only the serialized descriptor. It derives
all internal profile state from the record before accepting packets. Encode,
decode, recover, and free use the separate opaque `WirehairV2Codec` type, so V2
handles cannot be confused with legacy `WirehairCodec` handles.

Parsing distinguishes invalid magic, unsupported encoding version, invalid
size, nonzero reserved fields, unknown equation profile, and invalid dimensions
through stable `WirehairV2Result` values. Codec failures additionally report
need-more, bad seed, resource exhaustion, OOM, and unsupported platform.

The installed C++ header `<wirehair/wirehair.hpp>` provides move-only RAII
`wirehair::v2::Encoder` and `wirehair::v2::Decoder` wrappers plus a fixed-size
`SerializedProfile`. Its three explicit `Encoder::CreateBorrowed()` overloads
retain no hidden source owner, and `Encoder::DetachInput()` releases the same
native lifetime obligation. Failed replacement preserves the prior encoder and
its borrow; moves transfer the handle and its obligation. The wrappers use the
same C ABI and byte contract.

## Version migration

Encoding version 1 accepts only a declared and supplied size of 32 and requires
all reserved bytes to be zero. A future record that changes semantics or adds
fields must use a new encoding version and receive explicit parser support;
older libraries return `WirehairV2_UnsupportedVersion`. Reserved bytes cannot
be repurposed while retaining version 1.

Any change to a frozen equation rule requires a new equation profile ID and new
golden packet vectors. A new equation profile may continue using the 32-byte
encoding when the existing fields suffice. Decoders must reject unknown IDs
instead of guessing a compatible profile. Migrations recover and authenticate
under the old supported profile, then re-encode under the new profile.

## Relationship to legacy profiles

V2 and legacy share the same profile-ID vocabulary: the first 64 bits of
SHA-256 over a documented canonical name, exposed as a `uint64_t` compatibility
identifier. Both require an application integrity/authentication layer and
reject unknown named contracts.

The representations intentionally differ. Legacy packets retain their
historical unframed API, so `WirehairWireProfile` is a 16-byte in-process
selector and the application chooses framing and byte order. V2 was not yet a
public wire contract, so it starts with a canonical endian-stable record that
also carries message dimensions and the selected attempt. Legacy and V2 profile
IDs are different equation namespaces and are never interchangeable.
