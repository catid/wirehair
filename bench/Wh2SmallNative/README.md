# K2/K3/K4/K5/K6/K8 native correctness harness

This standalone build tests the private compile-time small-block core. It does
not change the public library, select a wire profile, or measure performance.
The fixture generators authenticate the retained recovery bundles and
reconstruct their sealed tables without rerunning a selector or loss campaign.

```sh
cmake -S bench/Wh2SmallNative -B /tmp/wh2-small-native-build \
  -DWH2_SMALL_LIBRARY=/absolute/path/to/qualified/libwirehair.a
cmake --build /tmp/wh2-small-native-build
ctest --test-dir /tmp/wh2-small-native-build --output-on-failure
```

Use a fresh, task-owned build directory. This harness needs the original local
`/var/tmp/wh2-k3-thue-morse-r0` evidence; missing or changed evidence fails closed.
The ordinary library and installed package do not depend on this directory.
Pass `-DWH2_SMALL_TEST_DIMENSION=5` to select the separately sealed K5 evidence
at `/var/tmp/wh2-k5-thue-morse-r0` instead. Dimension 2 selects the separately
sealed `/var/tmp/wh2-k2-thue-morse-r0` evidence. Dimension 8 selects
`/var/tmp/wh2-k8-thue-morse-r0`. Dimension 4 selects
`/var/tmp/wh2-k4-thue-morse-r0`. All five dimensions build two direct-core
tests and five serialized-boundary tests. The K2/K4/K8 boundaries
are benchmark-only; they do not introduce installed profiles.
Other dimensions are rejected. K3 remains the default build; this option never
changes the installed library.

The private `gf256.h` compilation settings must match the linked library's
settings. For an `ANDROID`-selected portable-GF library, pass
`-DWH2_SMALL_PORTABLE=ON`. This checks the portable backend, **not** an Android
platform or real non-GFNI-host performance. Sanitizer compile and link flags
must likewise match the selected instrumented library.

`--neutral` checks the selected dimension and K6 widths and partial tails,
field products, allocation failures, aliases, dependent/conflicting packets, and K6 parity with the
untouched original core. `--corpus` checks all 7,774 retained K3 payload cases,
including the 6,144 fresh and 72 hard traces, historical original-width
prefixes, and every triple in the development/seam windows. Neither mode
compares speed with WH1 or establishes a new failure-rate sample.
At K5 the corpus has 10,053 cases: 6,216 retained fresh/hard traces, 57
original-width historical prefixes (54 unique ID sequences), and all 3,780
five-of-nine seam subsets. It checks 75,134 packets and every corresponding
prefix rank and recovery, plus 2,270 recorded coefficient rows. Every encoder
is freed before creating any receiver; repeated recovery uses reset output
buffers and every actual packet is checked against a polynomial-field oracle.

At K2 the corpus has 15,956 cases and 56,776 packets: all 6,216 retained
hard/fresh traces; all 56 historical origins without deduplication or changes
to their original widths, tails or prefix lengths; 450 seam subsets; and the
three legacy plus 1,536 stride pairs at all three widths with full and one-byte
tails. It also checks all 5,274 recorded coefficient rows. Five retained stride
pairs are deficient, giving 30 shape replays that must stay at rank one and
leave recovery buffers unchanged. These cases are not extended with rescue
packets. All other cases must recover exactly. Neutral coverage additionally
includes 19 selected widths spanning 1 to 4,096, partial tails, every packed selector, aliases,
OOM, conflict preservation and allocation-free packet operations.

At K4 the corpus has 8,424 cases and 58,571 packets: all 6,216 hard/fresh
eight-ID traces with all five recorded prefix ranks; all 38 historical origins
at their original widths, tails and horizons (163 packets, 37 unique ID
prefixes); and all 2,170 four-of-eight seam subsets, replayed at B2/full tail.
All 2,252 recorded coefficient rows are checked against native mapping and an
independent polynomial oracle. The lookup is exactly 20,480 bytes. No history
is deduplicated or extended and no new loss traces are generated.

At K8 the corpus preserves 21,110 cases and 193,746 packets: 6,216 retained
hard/fresh twelve-ID traces with all five recorded ranks, all 44 historical
origins with original widths/tails/prefix lengths (354 packets), and all
14,850 eight-of-twelve seam subsets. It checks all 2,347 recorded coefficient
rows against both native mapping and an independent polynomial oracle. The
selected companion parameter is 2, not the parameter 1 used by smaller
dimensions. The native core admits K8 only through compile-time assertions;
no decoder algorithm, production profile, installed lookup or default changes.
Additional highest-pivot tests feed ID7 first and last over 108 width/tail/order
shapes, including duplicates, contradictions, repeated recovery and allocation
checks. Corpus replay uses recorded terminal ranks and never appends rescue
packets. This is retained-data correctness validation, not new recovery or
timing evidence.

## K4 native qualification

The sealed `(64,120,54,15)` / `(65,120,54,15)` pair now instantiates the
existing native core. Only three compile-time dimension assertions and a
comment change in the core; there is no new runtime algorithm, payload
kernel, installed lookup, public API, wire profile or default. K4 uses the
generic payload operation, not the K6-only two-byte specialization.

All 90 tests pass: two K4 direct-core tests plus all seven existing K2/K3/K5/K8
direct/serialized/K6-parity tests, in each of native, portable-arithmetic and
ASan/UBSan builds. Sanitizers instrument the harness and linked library;
leak and fake-stack detection are enabled. Every backend checks the complete
8,424-case K4 corpus, 58,571 packet/feed/prefix-recovery checks and 2,252 rows.
All nine retained OH0 deficiencies agree with the original ranks and recover
with the recorded extra packet. These replays add no recovery samples.

K4 neutral testing checks 162 width/tail/order cases alongside 162 K6 parity
cases, 5,632 selector rows, ownership lifetimes, every constructor allocation
failure, aliases, error outputs, dependent conflicts and repeated recovery.
Another 108 K4 highest-pivot shapes exercise 432 systematic feeds. Every
payload replay destroys its sender before receiver creation. This is the
private core's borrowed-source contract, not qualification of a public ownership
facade.

All 31 fixture-generator tests pass on Python 3.12 and 3.8. Independent
source review caught a copied negative test that left the first K4 seam ID
unchanged; it was corrected before the first test/build run. Repeated main
and independent complete source reviews find no remaining confirmed bugs.

The four affected native production translation units (`Small`, `SmallK5`,
`SmallK8`, `V2Profile`) were each recompiled with baseline-header, candidate-
header and normal inclusion. All twelve outputs exactly match the original
archive members. The archive and shared library consume those same original
objects. This proves native code-generation identity, not speed qualification.

Builds and complete logs are retained at `/tmp/wh2-k4-native.j0iJxQDC`.
`QUALIFIED.json` SHA256:
`1317c222c7b9a32e02a3831a76b817e44edbfa4bb72d6ff0ae0d2947546443a2`.
The generated fixture has 534,629 bytes, SHA256
`608bafbe37ac0ba3aa94f5d030390af6cc623f9cf0791e86c5bcb8d72f277763`.
The independent audit at `/tmp/wh2-k4-native-independent.f65PiBZz` verifies
all 895 input identities, every initializer, all test logs and dependency/ABI
settings, and the object proof. Its Python 3.12/3.8 reports are byte-identical,
SHA256 `ca7d5cf1e3ef38758d7044c9bd836800b6d53225cc8eb14872f32bff327ece83`.
Five auditor neutral tests pass per interpreter. Four auditor-only development
exits (formatting, access-time comparison, interpreter hardlinks and Ninja
link evidence) are disclosed in its `DEVELOPMENT.json`; they changed no codec
or retained test result.

Do not reconfigure or run CTest discovery in these archived builds. Source-
sensitive audits completed before this documentation advanced. At this native
milestone, the next step was the serialized ownership boundary; its separate
qualification is recorded below. Actual ownership-matched WH1 full-lifecycle
speed and retained paired recovery gates remain required. Ordinary K4 still
selects certified WH2; no promotion or all-K claim follows from this milestone.

## K8 native qualification

All 69 tests passed: K8 direct neutral/corpus tests and all seven K2/K3/K5
direct/serialized/K6-parity checks, each in native, portable arithmetic and
ASan+UBSan builds with leak and fake-stack detection. Each K8 backend checked
21,110 corpus cases and 193,746 packets, 324 neutral cases, and 108 additional
highest-pivot shapes with 864 packet feeds. All 25 fixture-generator tests
passed under both Python 3.12 and 3.8. Backend replays are correctness checks
of retained data, not extra recovery samples or performance measurements.

The separate source reviews and independent artifact audit found no confirmed
bugs. Both interpreter audits checked every K8 fixture value, all 69 test
results, build flags, library/private-GF ABI and dependency identities. Eighteen
production compilations (`Small`, `SmallK5`, `V2Profile`, static/shared,
baseline-forced/candidate-forced/normal) produced six byte-identical three-way
groups. The object-check script initially stopped because shared CMake also
exports internal-support compile entries; selecting the exact public target
resolved this audit-only ambiguity. The completed final objects and original
partial attempt are both retained. No codec result was discarded or retimed.

Artifacts: `/tmp/wh2-k8-native.ncbVY3vR`, with `QUALIFIED.json` SHA256
`dbb31bb68bcee09934537de91453a094c695df3d3815d92d98b280b7a5c4997f`.
The independent Python 3.12/3.8 reports at
`/tmp/wh2-k8-native-independent.pdP8n3ii` both have SHA256
`aef362ab0f718fa5e91e4e7f45fb4f15cfca9d16b13f437738eae0e529d33ab8`.
The generated fixture SHA256 is
`99e2b20da1405755136311d58bd7391474866bdee1a10ccf22a2de7482b4aaea`.
Evidence was audited and sealed before this result paragraph was appended;
the qualification reports pin the pre-result README, not this later text.

This milestone accepted native K8 only. Its external ownership boundary, actual
paired WH1 recovery and full lifecycle timing require separate qualification. No
installed lookup, public profile, default, speed or all-K claim is made.

## K2 native qualification

The fixed `(2,3)` / `(3,3)` companion pair now works through the existing
native core without a new decoder, payload kernel or runtime branch. The only
core edits admit dimension two in its three compile-time assertions.
The lookup remains a benchmark-generated, 7,168-byte fixture; no installed
lookup, public profile, defaults or exported API change.

Qualification artifacts are retained at `/tmp/wh2-k2-native.jtQp2mfD`.
The final `QUALIFIED.json` SHA256 is
`e5df706707e92698368aa464cbe318d2d21eda528521d641369feb50a958ed0c`.
All 48 tests passed: K2 direct neutral/corpus checks and all seven existing
K3 and K5 direct/serialized/K6-parity checks, each under native, portable
arithmetic and ASan+UBSan+leak/fake-stack configurations. The fixture and
source identities, complete generated values, reported counts, build flags
and mode agreement passed a separate artifact audit on Python 3.8 and 3.12.
All 19 K2/K3/K5 fixture-generator tests also passed on both interpreters.
Repeated source-reading passes found no core or harness bugs. The artifact
auditor initially needed dependency paths containing `..` normalized; this
was corrected before both complete audits passed, without changing codec
or test outcomes.

The generated K2 fixture SHA256 is
`3c7000ed80af9f7f15cdb51dabb1a016cf24afe12988aa355605387c6f855cba`.
For every production translation unit that includes the core (`Small`,
`SmallK5`, `V2Profile`), baseline-forced, candidate-forced and normal builds
produced identical object bytes under both static and shared Release flags:
18 compilations, six identical three-way comparisons. The retained
`object-identity.json` SHA256 is
`539349f52389f7d30521a9b2d8f6d353938b4a4816038715fa281239ad899115`.
This is code-generation identity, not a new timing result or restoration
of previously measured regressions.

Native replay confirms the already retained recovery cases; it does not add
a fresh sample or a paired WH1 recovery-rate comparison. At this direct-core
milestone, the external boundary and full lifecycle speed were still unqualified.

## Serialized correctness for K2, K3, K4, K5 and K8

The build also produces a separate, non-LTO serialized C boundary and five
additional tests for K2, K3, K4, K5 and K8. Its shared `Facade<K,Traits>` implementation retains K6's
independent/borrowed input, transactional allocating detach, descriptor
validation and permanent conflict-poison behavior. Only the selected benchmark
wrapper instantiates that boundary for external callers; the installed K6
implementation remains untouched. K3 uses `WHK3` and profile ID
`0x5748324b33544d31`, not an existing or retired WH2/K6 identity.
K5 uses `WHK5` and `0x5748324b35544d31`, a 29,440-byte lookup, and the same
bounded allocation policy (`block_bytes <= floor(268435456 / 6)`). There is no
K5-specific decoder, duplicated API, or inherited six-source tiny-payload kernel.
K2 uses `WHK2` and `0x5748324b32544d31`, the sealed 7,168-byte lookup,
and `block_bytes <= 89478485` (`floor(268435456 / 3)`). Its admission changes
only compile-time dimension checks and fixture selection, not the shared
runtime algorithm or the ordinary WHV2 conflict contract.
K4 uses `WHK4` and `0x5748324b34544d31`, the sealed 20,480-byte lookup,
and `block_bytes <= 53687091` (`floor(268435456 / 5)`). This bounds the
five-block decoder slab, not aggregate or transient process memory. K4 uses
the existing generic runtime and lambda 1, not a new payload kernel.
K8 uses `WHK8` and `0x5748324b38544d31`, the same sealed 65,536-byte lookup,
and `block_bytes <= 29826161` (`floor(268435456 / 9)`). This limit bounds the
nine-block decoder slab, not aggregate or transient process memory. No K8
production facade or exported installed API is introduced. The independent
serialized oracle uses K8's lambda 2; smaller dimensions retain lambda 1.
K8's retained serialized corpus has 42,220 cases and 387,492 packets across
both policies, without changing any trace, width, tail or terminal rank.
Its additional external highest-pivot gate checks 108 ID7-first/last shapes,
864 primary feeds, and 216 separate permanent-poison receivers before and
after full recovery. The private-core conflict-and-resume checks remain
separate because the external boundary intentionally poisons permanently.

`small_serialized` tests the selected dimension's ownership/error/byte behavior;
`small_serialized_k6`
instantiates the same template at K6 and checks installed K6 descriptor,
packet, feed, recovery and detach parity. `small_serialized_c` verifies the C
ABI, including decoder creation before any encoder exists. The remaining two
tests replay the neutral and retained selected-dimension corpus through the
external serialized functions, separately for independent and borrowed input. These two policy
replays do not increase the number of recovery samples. All seven tests are
correctness checks, not speed measurements.
When enabling sanitizers, instrument both C and C++ compile flags and the link
flags; the supplied archive must use the same sanitizers/backend configuration.

### K2 serialized qualification

All 63 tests passed: all seven checks at K2/K3/K5 under native, portable
arithmetic, and ASan+UBSan with leak and fake-stack detection. K2's external
boundary passed 78 lifecycle shapes, including both input policies, partial
tails, transactional detach, every allocation-failure site, aliases, malformed
descriptors and permanent conflict poison. The separate C consumer created a
decoder from a literal descriptor before any encoder existed.

The retained K2 corpus replayed 31,912 serialized cases and 113,552 packets
across both policies, with independent packet/prefix-rank/recovery checks.
All 60 deficient shape/policy replays remained `NeedMore` with untouched output;
no extra packets were appended. These are repeated engineering checks of the
same retained cases, not additional recovery samples.

Artifacts are retained at `/tmp/wh2-k2-serialized.rnxDJH4G`; `QUALIFIED.json`
SHA256 is `616ecfb9727a4ad48990d3fa00cf8b559ca2f6fd5e767c83c906269b2fbbdf6c`.
A separate artifact audit passed on Python 3.8 and 3.12, checking source and
library hashes, build flags, every generated K2 fixture value, K3/K5 fixture
identity, test counts and backend agreement. All 19 fixture-generator tests
also passed on both interpreters. Repeated source-reading passes, including a
separate local review, found no boundary bugs. An initial artifact-auditor
syntax error was corrected before either recorded audit; codec tests were
unaffected. No Fable report was obtained for this milestone.

Full lifecycle speed against ownership-matched WH1 and ordinary WH2 remains
unmeasured for K2. No installed API, profile, default or production source
changed in this serialized-boundary milestone.

### K4 serialized qualification

At base source `84eeea9`, five benchmark files admit the exact sealed K4
fixture and extend dimension-specific oracle, descriptor and highest-pivot
tests. The facade runtime bodies, private core, production sources, installed
profiles and existing library archives are unchanged. There is no installed
K4 API or ordinary-WH2 ownership/poison semantic change.

All 105 tests pass: all seven checks at K4/K2/K3/K5/K8 under native, portable
arithmetic and full ASan/UBSan with leak and fake-stack detection. K4 passes
78 lifecycle shapes, 108 highest-pivot shapes with 432 primary feeds, and
216 separate permanent-poison receivers. The literal-descriptor C consumer
creates its receiver before any encoder exists. Independent source copying,
borrowed-source lifetime, transactional allocating detach, every allocation
failure, aliases, partial tails, malformed descriptors and repeated recovery
are checked. Packet operations allocate no memory.

Each backend preserves all 8,424 direct-core cases, 58,571 packets and 2,252
coefficient rows. The serialized corpus checks 16,848 cases and 117,142 packets
across both source policies; serialized neutral replay adds 324 policy cases
and 4,212 packets. Every original history keeps its width, tail and horizon.
The nine retained OH0 deficiencies remain deficient until the recorded extra
packet. These engineering replays do not add recovery samples.

Builds and complete logs are at `/tmp/wh2-k4-serialized.yegpP9w9`.
The native, scalar and sanitizer `RESULT.json` SHA256 values are respectively:

- `90f743af46b9d7d2ae4e7869e2ed148900df130289d28504ebdf500c98147a8c`
- `55be9432eaca5fb9399dc0ef18ef56688f459c03cfec70c784a86de35501427e`
- `c7ef4675b33335dbc5f6ffaddbd864c4cbc716b82192173bcfff1c95e52d077f`

Independent source/artifact audits at
`/tmp/wh2-k4-serialized-independent.VQMU5xiC` pass on Python 3.12 and 3.8.
The byte-identical 3,611-byte reports have SHA256
`c0a5afa2e2250da819bde1b9eeefb9e13e895ce7a81666564fb143ff766054d9`.
They verify all 924 pre/post input identities, every fixture initializer,
all 105 complete test outputs, nine actual C exports, non-LTO translation-unit
boundaries, compiler dependencies and private GF ABI/sanitizer settings.
All 398 other source files match the native qualification pins; its four
production object-identity groups are reauthenticated without recompiling.
Six auditor selftests pass per interpreter; no diagnostic development exits
occurred. Repeated main and independent source-reading passes find no
remaining confirmed bug.

These audits completed before this result documentation advanced. Do not
rerun source-sensitive historical auditors, reconfigure archived builds or
run CTest discovery there. Separate ownership-matched actual WH1/certified-WH2
retained recovery and full encoder/decoder lifecycle timing remain required.
This qualification changes no production default and establishes no speed,
paired recovery superiority or all-K result.

### K8 serialized qualification

All 84 tests passed: all seven checks at K8/K2/K3/K5 under native, portable
arithmetic, and ASan+UBSan with leak and fake-stack detection. K8's external
boundary passed 78 ownership/lifecycle shapes, the separate literal-descriptor
C receiver, 108 highest-pivot shapes and 216 permanent-poison receivers.
The complete retained corpus passed 42,220 serialized cases and 387,492 packet
checks across both ownership policies per backend; no history was extended,
deduplicated or given different widths/tails. Native prefix-rank and guarded
payload recovery match the same sealed mathematical records.

Two source reviewers and repeated local reading found no boundary bugs. An
independent artifact auditor verified all tests, every generated fixture value,
222 pinned artifacts, compilation/link settings, private GF ABI and complete
repository dependency coverage. Python 3.12 and 3.8 audit reports are identical.
All 25 fixture-generator tests also passed under both interpreters. The core,
production sources, installed profiles and existing library archives remained
unchanged; only six benchmark files implement/document this boundary extension.

Artifacts: `/tmp/wh2-k8-serialized.qQHyZqaa`, with `QUALIFIED.json` SHA256
`991c5a69b5fa22d34ad0d51922e81dd64a678a867f38add7d121c05b342bd1a4`.
The independent reports at `/tmp/wh2-k8-serialized-independent.sbcgjBaH`
both have SHA256
`d59235955675ebfda2eded46e25054f96d42752b8ba4dcb9ad6fecd5f11f8f33`.
The source/fixture hashes were audited and sealed before this result text.
`README-pre-result.md` in each K8 native/serialized artifact directory preserves
the exact non-producing README bytes pinned by its qualification report.

This is benchmark-boundary correctness qualification, not an installed profile,
new recovery sample, paired WH1 recovery result or speed promotion. Actual
ownership-matched WH1 recovery and full encoder/decoder lifecycle comparisons
remain required in `wirehair-sxvz.16.1.20.82.5.2`.

## Separate K3 serialized performance qualification

At source `5caa4a9`, the separate `Wh2K3SerializedCostR1.py` gate passed all
54 same-code timing controls and all 36 comparisons against actual WH1 and
current public WH2. It retained 2,488,320 fresh codec lifecycles, checking every
output. Time reductions versus WH1, spanning both measurement orders:

| Block bytes | Encoder | Low-ID decoder | Distant-ID decoder |
|---|---|---|---|
| 2 | 65.8-65.9% | 95.4% | 91.0% |
| 64 | 72.7% | 95.1% | 91.3% |
| 1280 | 75.2-75.5% | 93.9% | 92.2% |

These are serialized-prototype, static-call results on one GFNI-capable host,
using borrowed immutable input and three full source blocks. Encoder work
includes creation, descriptor output, 18 packets and free; decoder work
includes creation, feed through first success, recovery and free. Both decoder
streams reached success after three packets for every arm. This is not an
integrated-library, cold-start, shared-call, partial-tail, non-GFNI or all-K
speed result. This prototype gate did not change the normal library or defaults.

The immutable outcome is `/var/tmp/wh2-k3-serialized-cost-r1`. Its raw stream
SHA-256 is `35d3ce567924ae21418f0b16986e6cf89a02a162e0b2c64bab6d89036a250604`.
A Python 3.8 full replay verified all 601 receipt pins; a separately written
raw chronology, API ledger and confidence-interval audit reproduced all 90
decisions. The R0 run remains invalid because its neutral diagnostic truncated
the binary CPU identity at an embedded NUL. R1 corrected that transport and
added prelaunch validation without changing the candidate or performance gate.
Neither namespace may be rerun or its observations filtered or rescored.

The retained K3 recovery screen had no zero-overhead failures in 6,144 traces
and passed all 72 hard traces. Native and serialized payload replay confirmed
those same cases; this timing run adds no recovery-rate sample or paired WH1
recovery-rate comparison. Production integration must preserve the sealed
equations and explicit descriptor, then qualify the actual library separately.

The later opt-in integration is documented in
[`SMALL_WIRE_PROFILES.md`](../../SMALL_WIRE_PROFILES.md). Its separately frozen
`Wh2K3ProductionCostR0.py` actual-library gate at `ffa7739` also passed all
54 same-code controls and 36 WH1/current WH2 comparisons; see that document
for its distinct measurements and scope. `Wh2SmallProductionParity.cpp` reuses
this correctness corpus
through both the old serialized prototype archive and the new library,
comparing descriptors, packet bytes, prefix status and recovered bytes. It is
an engineering replay, not a timing or fresh recovery experiment, and does
not alter the frozen source harness or old outcome bundles.
