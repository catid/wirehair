# K2/K3/K5/K6/K8 native correctness harness

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
sealed `/var/tmp/wh2-k2-thue-morse-r0` evidence. All three selected dimensions
build two direct-core tests and five serialized-boundary tests. Dimension 8
selects `/var/tmp/wh2-k8-thue-morse-r0` and builds only the two direct-core tests;
its external ownership boundary is not yet implemented. The K2 boundary
is benchmark-only; it does not introduce an installed profile.
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

This accepts native K8 only. Its external ownership boundary, actual paired
WH1 recovery and full lifecycle timing require separate qualification. No
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

## Serialized correctness for K2, K3 and K5

The build also produces a separate, non-LTO serialized C boundary and five
additional tests for K2, K3 and K5. Its shared `Facade<K,Traits>` implementation retains K6's
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
