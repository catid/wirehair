# Ordinary K8: native shared-library preparation

Tracked in `wirehair-sxvz.16.1.20.82.5.3.4.3`. This directory contains neutral
build and API checks, not a timing worker, loss campaign, or default change.

The candidate is the already qualified six-line ordinary K8 selector. Its
Profile source and object must reproduce exactly. All nineteen original native
producer objects are rebuilt with the original flags; the full baseline DSO
must reproduce byte-for-byte. Candidate linking replaces only Profile, in its
original position, preserving the export map, SONAME and link options. It does
not copy the static benchmark's candidate-before-archive arrangement.

Use a fresh task-owned external parent directory, and separate new children:

```sh
python3 bench/Wh2K8OrdinaryShared/Prepare.py /absolute/fresh/parent/prepared
python3 bench/Wh2K8OrdinaryShared/Qualify.py /absolute/fresh/parent/prepared /absolute/fresh/parent/neutral
python3 -m unittest discover -s bench/Wh2K8OrdinaryShared -p 'test_*.py' -v
```

Never point either builder at a retained qualification directory. The original
K8 qualification, static selector, prior scientific evidence, and historical
CTest logs are read-only inputs. Neutral fixes can be checked in new directories;
this does not authorize repeating any spent scientific namespace.

## Actual shared-call checks

All thirty direct consumer executables link only the selected DSO, with no
static codec support. A provider check runs before `main`, resolves five public
entry points with `dladdr`, and requires the exact intended library. Each
consumer is also deliberately launched against the opposite arm and must fail
before running its test. This matters because both libraries preserve the same
`libwirehair.so.2` SONAME.

The 31 required CTests include baseline/candidate K3/K5/K8 C++ ownership,
allocation-failure, alias, conflict and partial-tail suites; borrowed-source
contracts; C/C++ consumers; independent metadata validation; actual dynamic
exports; and complete certified-byte parity. The candidate adds ordinary K8
route, boundary and partial-overlap checks. The exact roster and every passed
JUnit record are checked; a missing, disabled, skipped or duplicated test fails.

The allocation suites interpose C++ allocation functions and retain their exact
failure/count assertions. They do not substitute the root build's static
`wirehair_test_support` library. The separate environment-controlled facade
fault hook is absent from production-style DSOs and is not claimed by this
qualification. Its prior test-only/static evidence remains separate.

In two fresh processes, `Bindings.py` loads baseline/candidate in opposite
orders with `RTLD_NOW | RTLD_LOCAL`, with no globally linked Wirehair. It checks
all 53 public exports, six internal public GOT targets, 37 runtime GOT targets,
and distinct matching GF contexts. It then exercises 48 deterministic ordinary
API cases per library, including K3/K5/K8/K16, both widths/tails and bare or
independent/borrowed options, with sender destruction before receiver creation
and guarded repeated recovery. It checks bindings again afterward. These are
engineering cases, not independent loss samples or timing observations.

This explicitly isolated loading model is not a claim that two globally loaded
libraries with the same public symbol names are immune to ELF interposition.
No `-Bsymbolic`, renamed SONAME, or other binding change is used to hide that
distinction.

## Native qualification result

The final native preparation and direct-shared qualification both completed
successfully on 2026-09-10 at
`/tmp/wh2-k8-ordinary-shared-final.Iy6tBm1H`. All 31 required CTests passed;
all thirty opposite-provider launches rejected the wrong library. Each of the
two load orders passed 48 actual API cases per DSO and all public/internal/runtime
binding checks. The complete certified streams match the retained pre-admission
2,180,292-byte reference. Eleven read-only/synthetic tests pass under both Python
3.12 and 3.8; the complete new preparation verifier also passes under Python 3.8.

The nineteen original objects and complete baseline DSO reproduce exactly.
The candidate uses the exact previously qualified selector object, in the
original Profile slot, and keeps the original 53 exported functions and one
matching 141,328-byte native GF context per DSO. No equation or production
source changes accompany these tests.

| Retained artifact | SHA-256 |
| --- | --- |
| `prepared/PREPARED.json` | `d98b53db95464544d73b4f1cc5f68eaf5b00d8c6062b23e5a91112e7a55b7379` |
| Baseline shared library | `ef684abac606667e30e3b0de1204b0897aa0fa93a4f5ec2268d1ba951dce03b2` |
| Candidate shared library | `fe7527a7470e7ea8e1790eb84b3761689e2bb61c883d6df905a1b3889699df3b` |
| `neutral-closed/QUALIFIED.json` | `db047c9c3d7eb1f526cb99596630a255bc8717b042498317af0784c24decce82` |
| `neutral-closed/ctest-results.xml` | `52420e5d58ea24bff8c5d45553179c54c33b6e702b64df4d537dfc4fe7bc4a2e` |

Preparation records 627 input/output pins, 62 successful commands, and 268
captured command/dependency artifacts. The separate consumer qualification
records 1,646 file pins and 186 successful commands, in addition to all thirty
wrong-provider outcomes. Exact command graphs, twenty producer dependency files,
compiler/tool/runtime inputs, consumer source/header closure, generated test
roster, checker scripts, binding reports and negative captures are authenticated.

Earlier neutral builds at `/tmp/wh2-k8-ordinary-shared-prep.NdDkOaWn` and
`/tmp/wh2-k8-ordinary-shared-closed.LWi1csHD` are retained. Their actual codec
checks passed, but source reviews exposed incomplete verifier input/capture
closure. Those omissions were fixed and synthetic deletion/alteration checks
added before the final fresh qualification. A Python 3.8 test initially assumed
`_ctypes` was an external module; that interpreter builds it in, so the test now
accepts the separately pinned interpreter representation. No codec candidate,
equation, measurement threshold, loss trace or scientific result was changed.

The additional `/neutral` run also passed every codec check. Its independent
records audit found that the invoked `ctest` executable itself was not pinned.
The final `/neutral-closed` run adds that tool/runtime pin and requires every
invoked tool and consumer to be frozen before execution. Independent final
source and record review verifies 1,917 unique combined evidence pins and finds
no remaining issue. The earlier result and its omission remain retained.

This milestone is native shared correctness and provenance qualification only.
Shared timing, all-K claims and default promotion do not follow. Relocated
package checks are qualified separately below.

## Relocated installed C/C++ and plugin qualification

`Package.py` copies the retained native shared installation into two fresh
staging prefixes, replaces only the candidate DSO, and moves both copies to
new relocated prefixes before configuring consumers. The original installation
is never modified. All nineteen regular files and both relative SONAME links
are preserved, except the exact candidate DSO bytes/length. This is an overlay
of the installed package metadata, not a new production install/default change.

```sh
python3 bench/Wh2K8OrdinaryShared/Package.py /absolute/prepared /absolute/fresh/package-output
```

The unchanged `test/package` C/C++/K3/K5/K6/K8/plugin consumers build against
the actual imported CMake target. Its shared-library type, relocated location
and include directory, `WIREHAIR_DLL=1`, and C++11 requirement are checked.
No repository public headers or static codec implementation are used. Two
ordinary C++ consumers start with a C++98 requirement and rely on the installed
target to raise it to C++11. Each checks both policies, B2/B64, full/one-byte
tails, the expected descriptor, move/detach behavior, repair packets after
source release, sender destruction before receiver creation, and twice-guarded
recovery. The candidate also runs the ordinary C constructor against the
unchanged independently derived literal K8 descriptor and repair bytes.

The final run at `/tmp/wh2-k8-ordinary-package.8vNek285/neutral-helpers`
passes all nineteen required tests (baseline nine, candidate ten), including
both actual dynamically loaded plugins. Every one of the nineteen codec-linked
executables/plugins also rejects the opposite-arm library before application
work; the two plugin hosts themselves have no Wirehair dependency. There are
46 compiled TUs, 21 link maps including the hosts, 159 captured commands and
1,317 frozen input/output pins. Before/after compiler dependencies, actual
link/runtime providers, installed metadata and unchanged original package are
checked. All seventeen synthetic/read-only tests pass on Python 3.12 and 3.8.
Repeated main-agent and independent complete source passes are clean. The
independent standard-library-only records audit verifies 1,935 unique combined
evidence pins, all command/capture/environment cross-links, both imported
targets, every consumer/provider and both unchanged package inventories except
the candidate DSO. It executes no codec or historical verifier.

| Package artifact | SHA-256 |
| --- | --- |
| `PACKAGE_QUALIFIED.json` | `625e5236976fea2d315320a67f926572b8f136badc909cd9d4644200b68eaeaf` |
| Baseline JUnit | `567f4fb1aafcc28d960d017c67e095bf2b4864aff8f0cb4beffdaf933d71d5d6` |
| Candidate JUnit | `567d424fe5bbb31bc7b8634cd0aba57dc676dfabf3558eabcc038f947909f00d` |

The initial `/neutral` attempt stopped before configuration: the package checker
incorrectly assumed equal DSO sizes (baseline 686,992 bytes, candidate 687,008).
It now uses the exact prepared size and hash. The `/neutral-sized` attempt
built the baseline but rejected the package's own `wirehair_package_round_trip`
test helper as a codec implementation. Only the exact upstream test-helper
names are now allowed on their expected targets; actual codec/GF symbols still
reject. Source review also corrected a per-TU header expectation: thin upstream
wrapper files legitimately include only `roundtrip.h`, so installed public
headers are checked per target, while DLL definitions and forbidden repository
header fallback are checked per TU. All failed artifacts are retained. These
are harness corrections, not codec changes or repeated scientific experiments.

This qualifies native relocated C/C++/plugin use only. It is not Python-package,
non-native-platform, sanitizer-codec, shared-speed, recovery-rate, or cold-start
qualification. Existing backend correctness evidence remains separate.

## Remaining admission work

Native shared and relocated-package correctness are preparation milestones.
The [prospective gate specification](Gates.md) separates current-path retention,
ordinary K8 superiority, and ownership-matched K3 speed/retention. The actual
shared timing instruments and their independent qualification remain required
before any new scientific timing launch.

The next gates must distinguish three questions:

1. Does the exact selector candidate preserve current shared-path performance,
   including current ordinary K3 and the explicit certified/K5/K8 paths?
2. Does actual ordinary K8 retain its speed advantage over ownership-matched
   WH1 through the shared-call boundary? Static timing cannot answer this.
3. What genuinely reduced-work change restores the previously measured
   preserved-path regressions against the pre-admission library?

The selector alone is not a mechanism for the third question. The retained
certified cases explicitly select their profile and never execute the new K8
branch; ordinary K3 returns before it. A change in code placement is not itself
evidence of restoration. Do not launch a blind selector-only rerun of the old
restoration workload, or revive rejected branch-hint, validator, serializer,
division, alignment or release-order variants.

For new K3 timing, match source ownership explicitly. The earlier ordinary K3
worker used borrowed WH1 and borrowed certified controls against both WH2
source policies. Its historical measurements remain unchanged, but independent
WH2 was not compared with owned-input WH1. Independent controls in a new gate
must use `wirehair_encoder_create_owned_ex`, not the borrowed constructor.

Keep candidate-versus-current retention, pre-admission restoration, and actual
WH1 superiority as separate decisions. Freeze the workload, ownership, first
success endpoints, measurement/load orders, controls and resource bounds only
after the final actual shared boundary is independently reviewed. No timing
namespace or default promotion is established by these preparation scripts.
