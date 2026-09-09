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
extra recovery samples. This screen does not establish partial-tail,
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
backend replays do not multiply its sample size. The namespace is spent; the
K5 speed-control failure, installed admission, preserved-path regressions and
full all-K speed/recovery/construction objective remain unresolved.

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
