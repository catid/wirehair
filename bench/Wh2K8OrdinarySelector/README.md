# Ordinary K8 selector: isolated correctness qualification

This benchmark-only candidate adds one six-line `SmallShape<8>` dispatch after
the existing K3 dispatch in `wirehair_v2_encoder_create`. It calls the existing
`WIREHAIR_V2_PROFILE_SMALL_K8_2026_09` constructor. The ordinary options
constructor already delegates to that entry point. Production source and
defaults are unchanged; neither serializer nor validator shortcuts are included.

The candidate preserves the sealed K8 equations, table, owned prepared basis,
source-independent repairs, allocation-free detach and nonpoisoning conflict
semantics. Explicit certified/CURRENT, K3/K5 and retired profile identities are
not changed. Larger K8 shapes still fall back to the certified profile.

## Retained result

At production source HEAD `7c4feb92f3538f7a37635d514931de6668f2b639`, the final
external build `/tmp/wh2-k8-ordinary-qualified.97E2f2VQ` passed **40/40** CTest
checks on 2026-09-10. These are correctness tests, not timing observations.

The build was configured with `CMAKE_BUILD_TYPE=Release` and
`WH2_K8_ORDINARY_NEUTRAL=ON`. The completed build/test commands were:

```sh
cmake --build /tmp/wh2-k8-ordinary-qualified.97E2f2VQ -j 12
ctest --test-dir /tmp/wh2-k8-ordinary-qualified.97E2f2VQ --output-on-failure -j 6
```

These commands document the completed qualification. For further changes use
a fresh external build; do not reconfigure or run CTest discovery in retained
qualification directories, where it can overwrite diagnostic logs.

The build uses GNU 13.3/Linux, strict C11/C++11 and no LTO. Native and forced
portable-arithmetic targets link the corresponding hash-checked installed K8
archives. ASan/UBSan instruments the candidate Profile TU, tests and retained
library; `-march=native` matches the archive's private GF ABI. Its fresh Release
Profile TU retains optimization and is not object/optimization-flag identical
to the retained debug sanitizer archive. Leak and fake-stack checks are enabled.
Separate `WIREHAIR_TESTING=1` facade objects serve only the borrowed-fault suites.

## Coverage and limits

- Derived full small-profile tests: K3/K8 each 69 width/tail shapes across nine
  constructor/policy routes; K5 remains explicit-only across six routes. K8
  includes 8,704 packed-lookup packet checks and 96 pivot-order cases.
- 198 neighboring shape/constructor/policy cases at K2/3/4/5/6/7/8/9/16/64/128,
  block widths 1/2/64 and full/one-byte tails. The two B1 tail cases coincide;
  these are not 198 independently distinct message shapes.
- 48 first-allocation probes around `7B`, `7B+1`, `8B`, `8B+1`, with
  `B=1,2,29826161,29826162`. Inaccessible virtual source memory and immediate
  allocation failure check the actual small/certified route without reading
  the payload or allocating a huge physical source. Explicit K8 just above its
  bound rejects without allocating. This is not maximum-size payload recovery.
- Source ownership, protected-source repair/detach, allocation failures,
  capacity/alias errors, duplicates/conflicts, C++ replacement across codec
  types and failed-replacement preservation. Partial descriptor/source overlap
  checks every systematic packet, the partial tail and low/distant repairs
  after source poisoning, as well as exact publication-span guards.
- Full unchanged K9 borrowed-source suite and the hash-derived K8 counterpart;
  full borrowed-fault suites against baseline and candidate. The C11 consumer
  calls the ordinary options constructor with unchanged literal descriptor and
  independent repair oracle.
- Independent metadata oracle on baseline/candidate in all three modes:
  66,928 host and 66,820 wire cases each. The full legacy Profile suite runs
  on baseline only: it intentionally expects ordinary K8 to remain certified.
  The candidate runs selected equation-independent overlap matrices plus the
  new ordinary-K8 assertions, not that full legacy suite.
- The entire retained K8 corpus runs through actual ordinary construction in
  all three modes: 21,110 cases, 193,746 packet/feed/recovery checks and 2,347
  core row checks. Both-policy serialized parity covers 42,220 cases and
  387,492 packet pairs per backend. Receivers are created after encoders are
  destroyed. Descriptor translation is explicit test code; valid operations
  are compared without equating the differing prototype ownership/conflict
  contracts. Row checks target the core, not a nonexistent installed row API.
  These replays reuse the original corpus and do not multiply recovery samples.
- Baseline/candidate explicit certified compatibility emits the same 2,180,292
  bytes over 48 cases, matching the pre-admission certified stream below.

The first neutral build, `/tmp/wh2-k8-ordinary-neutral.ZMnQ8KLr`, passed 28/31
tests. Its three K8 failures were a test mistake: the neighbor probe used a
small-codec helper that required allocation-free packet encoding, a guarantee
not provided by certified neighbors. The corrected probe checks packet bytes,
lengths and guards without imposing that unrelated contract. Candidate code
did not change. The initial failed log is retained. The final build also keeps
`initial-34-tests.log` and `intermediate-37-tests.log` from the successive test
additions. Independent source review prompted the stronger partial-overlap
payload assertions before the final test run.

## Exact artifact identities

| Artifact | SHA-256 |
| --- | --- |
| Production Profile source | `975da8d892363d05de3ad4535ec79b4a386bf3f6fab377b114af393f82d96d15` |
| Generated candidate source | `172324fd048daf154267f60a6e76df1e1b9a2545b9eac49059d0be21cfe99a3f` |
| Reproduced native baseline Profile object | `229da2af98169c140774fff3bc970a23b6688bba1586ae49fae862554aa64071` |
| Native candidate Profile object | `4c5e748f1fa618c88d8d9b8cd58045a1d4bf28211bceea63db4d726362ebf784` |
| Portable candidate Profile object | `78e0eb810723ed89cd7f18cbaa61984265578f5e8f242e06359a7d5f44bf493b` |
| ASan/UBSan candidate Profile object | `59404f5d4f1572df89ed7d32425bf2c741cdea12759fd097f988deb6853b750a` |
| Final `Testing/Temporary/LastTest.log` | `0665ec2927728dc0f890acdbca5a9ffdd3d762a8de99e53a2b5295a547cb4eef` |
| Both certified compatibility streams | `2e6536dcd86a7c2892399ddf1f14c3ff2290c2ef9e270aaa0c5ed87d0928907b` |
| Original retained K8 fixture | `99e2b20da1405755136311d58bd7391474866bdee1a10ccf22a2de7482b4aaea` |

The baseline object is byte-identical to the installed qualification object.
Source diff inspection confirms that the complete candidate differs only by
the six-line selector insertion. Final compile databases, link maps, generated
test derivations and executable artifacts remain in the build directory.
Repeated main-agent and independent bugworker source/artifact reviews are
clean. All 41 linker maps select the intended facade; allocation hooks occur
only in the six dedicated fault objects. Corpus links select only the
serialized facade from the prototype archive and use the matching retained
library's core/GF backend, not a second arithmetic context.

## Admission is still open

Issue `wirehair-sxvz.16.1.20.82.5.3.4` remains active. Actual ordinary-path
full-lifecycle speed, its separately recorded retained recovery comparison,
preserved-path restoration and current K3 speed retention are required before
default promotion. Explicit K8 speed evidence cannot be inherited by this
selector. This qualification claims no new recovery rate, shared-call or cold
speed, real non-GFNI-host speed, all-K success or zero bad construction seeds.
