# Dormant-core lifecycle screen

Issue `wirehair-fcsn.1`. This is a new screen of conditional unused-core lifetime,
not a retry of the rejected alignment or serialization candidates. No production
change or recovery-rate claim. Sole namespace:
`/var/tmp/wh2-small-dormant-core-screen-r0`.
**Terminal CONTROL_FAIL; candidate not qualified or promoted.**

The observer uses the exact qualified libraries from
`/tmp/wh2-small-dormant-core.qARBuvHF/`; it never rebuilds them. Native baseline
SHA256 `bd3c353847ec2838d95d5b73804664d3a2df8ac76b22ce4f5cf8c27fd58c1bbe`, candidate
`f9c2e226f180fb885561e65b0556df5b39c43e49bb53141ba9f4724cc5a8f9be`. The qualified
manifest, generated sources and original test logs must remain unchanged.

## Frozen workload and decisions

Six routes: ordinary K3; explicit small K3, K5 and K8; explicit certified K3
and K16. Each has B2-full, B2-tail1, B64-full, B1280-full and B1280-tail1, under
independent and borrowed-immutable source policies. Ordinary K3 is explicitly
distinguished from both explicit K3 routes. The actual selected profile is
deserialized and checked, not inferred from K. This is 60 fixtures / 180 cells.

Three complete lifecycle metrics: encoder construction + 3K packets + free;
low-repair-only decoder to its own first success + recover + free; distant-repair
decoder likewise. The 3K encoder packets are K systematic, K consecutive low
repairs and K distant repairs. WH1 uses corresponding owned/borrowed encoder
policies and its own repair packets and first-success count. WH2 arms must match
every descriptor, packet, repair input, first-success count and recovered byte.
No replay truncates WH1 to WH2's success count. Packet/recovery bytes and
padding/guards are checked outside timed WORK; descriptor equality and status
checks remain inside WORK for every arm.

All arms share actual source/output, decoder descriptor input, encoder descriptor
output and staged repair addresses within each replicate/cell. Raw records carry
these six addresses and the first-success count. The descriptor-output pointer
is the actual constructor destination, unlike the old screen's expected-profile
pointer. These are public workspaces, not private allocation addresses or every
scalar output/options pointer on helper stacks. No allocator interposition is
used. Arm-specific fixture staging can differ in cache state; public-address
matching does not prove identical private heap/cache state.

Each load order has five comparisons: B/B, C/C, WH1/WH1, C/B, C/WH1. CPU50
singleton; **64 lifecycles per observation**, 12 replicates, two observation
orders, 18 observations per panel: two retained warmups followed by eight
adjacent paired contrasts. Batch64 is chosen prospectively, before timing,
to reduce fixed timing noise relative to previous Batch32 lifecycle screens.
A Batch128 draft passed neutral correctness but its full-workload cost projection
from CTest elapsed time threatened the bounded screen cap; it was reduced before
freezing or claiming any timing namespace. No candidate timing ratios were seen.
The cell/pair rotation and side sequence are fixed in source.

For each cell/comparison/order, average the eight paired log ratios within a
replicate; use the resulting 12 values for a t11 95% interval. Strict gates:

- All 1,080 A/A intervals per load order must lie inside `[1/1.02,1.02]`.
- All 360 C/B cell intervals per load order must have upper bounds below 1.02,
  including all certified K3/K16 and all small routes.
- Six primary C/B intervals per load order (three metrics × two orders) must
  have upper bounds below 1. Each first averages all 40 small-fixture log ratios
  **within each replicate**, then applies t11 to the 12 replicate aggregates.
  This preserves within-replicate covariance; it never treats 480 samples as
  independent. Equal fixture weighting gives K3 half the small-roster weight
  (ordinary and explicit routes), and K5/K8 one-quarter each.
- Every C/WH1 interval is separately tested for upper bound below 1 and reported;
  it neither rescues failures nor turns an existing certified-route WH1 deficit
  into a claim that this change caused it.

Both DSO load orders must pass; no cross-order pooling. Any A/A failure overrides
all apparent gains; otherwise any failed retention or primary constraint fails
the screen. A geomean win is not every-cell improvement. The 2% retention envelope
permits small regressions; it does not prove zero regression. These are per-
interval screening bounds under replicate assumptions, not simultaneous 95%
confidence across the matrix. A precision-limited rejection remains a rejection.

Per load order: 388,800 observations (43,200 retained warmups), 1,800 cell
intervals and six primary intervals. Across both: 777,600 observations, 2,160
A/A constraints, 720 retention constraints and 12 primary constraints.

## Execution and preservation

Build observers in a fresh external directory with
`cmake -S bench/Wh2SmallDormantCoreScreen -B <build> -G Ninja -DDSO_DIR=<qualified-mode>`.
Use native, portable-DSO, and matching ASan+UBSan observer/DSO neutral checks;
the latter adds `-DSANITIZE=ON`. No neutral mode emits timing samples.
Run synthetic parser/statistics tests and repeated source review before launch.
Commit the frozen protocol, then run once:
`python3 -B bench/Wh2SmallDormantCoreScreen/Screen.py --run <native-build>`.

Each worker has a 240-second elapsed cap including preparation and checks, plus
a separate controller timeout of 270 seconds **per load-order worker**. WORK
sums including warmups must also stay below 240 seconds. These are bounded
screens, not all-K campaigns. Complete both orders after statistical failure;
abort on worker/input/format failure. Preserve all successful raw observations,
warmups and any emitted failure prefix/stderr. A hard termination can prevent
the in-memory raw buffer from publishing; it can never yield a successful receipt.
Every claimed namespace is permanently spent regardless of outcome: no retries,
trimming, resizing, retuning, subset rescue or repooling.

The controller checks the unchanged qualified codec manifest and a separately
hash-pinned observer qualification manifest. The latter binds final CMake/main/
test sources, all three generated workers, observer binaries/build definitions
and neutral test logs. Screen.py is committed/pinned separately to avoid a
circular seal. Only that exact native observer directory is accepted; a Ninja
dry-run must report no pending rebuild. Current C/C++ sources, screen source,
and all members of both manifests are pinned before launch and after each worker.
It rejects loader/allocator override
environment variables. This is not a complete transitive toolchain/loader seal.
Offline `--analyze` validates the frozen bundle's exact members/hashes and
recomputes the recorded decision; independent reconstruction and current-pin
checks must occur before advancing pinned source/HEAD after timing.

Even PASS only allows considering further adoption gates. It does not establish
all-K retention, universal speed against WH1, or better recovery. Conditional core
lifetime changes no recovery equation. Broader objective remains active.

## Pre-launch qualification

Final Batch64 observers under `/tmp/wh2-small-dormant-core-screen.mTN6m7uY/`
pass both neutral load orders in native, portable-DSO and ASan+UBSan builds:
six CTests total. Each verifies 60 fixtures × three lifecycle metrics × three
APIs. Sanitizer tests enable leaks, fake-stack detection and immediate UB failure.
The native observer defines no new/delete replacements. Native and portable-DSO
observers are byte-identical; only their loaded codec libraries differ.

Twelve synthetic tests pass on Python 3.8 and 3.12, including full chronology,
warmup exclusion/charge, first-success parity, descriptor-output address drift,
certified retention, strict gates, primary covariance and manifest replacement.
Review found and fixed missing binding of the tested observer to the final
protocol; a separate observer seal now prevents a stale/replaced worker from
being accepted by newly hashing it. Source reading also caught the original
codec manifest's repository-relative paths; these now have an explicit root.

`OBSERVER.sha256` pins 29 members; SHA256
`1d35223f3ca357c4b565c2d44d08950a5b8fbabad76914a61c3b8295cc567693`.
Generated worker SHA256
`0bcc44642106424932ade7c5beed4fc23c5d8dc6566d0f92c251fdbe843ae9af`;
native executable SHA256
`825bc1ae5167b798956340fe9e0525f5460ee3f7e2442f83b6d548bf4b5b46ef`.
Do not rebuild these observers or rerun CTest into their now-pinned test logs.

## Terminal result (2026-09-16)

The sole run at frozen commit `a2d669d` completed both workers with exit zero and
empty stderr. Every raw observation and warmup was retained. No retry, trimming,
subset selection, resizing or pooling. The namespace is permanently spent.

| DSO load order | Outcome | Failed A/A | Failed C/B retention | Nominal primary wins |
| --- | --- | ---: | ---: | ---: |
| Baseline first | CONTROL_FAIL | 3/1,080 | 0/360 | 6/6 |
| Candidate first | CONTROL_FAIL | 6/1,080 | 0/360 | 6/6 |

The equal-fixture primary point estimates below are **descriptive, not qualified
speed improvements**. Each range spans both observation orders and both library
load orders; it is not a pooled estimate or confidence interval.

| Small-profile lifecycle | Nominal time reduction versus current WH2 |
| --- | ---: |
| Full encoder | 3.32–3.65% |
| Low-repair decoder | 6.10–6.46% |
| Distant-repair decoder | 3.76–4.25% |

All 720 C/B point estimates lie below one, and all retention upper bounds lie
below 1.02. Neither observation overrides the failed controls. All three normal
A/A failure intervals and five of six reversed intervals include one. The other
failure is reversed explicit-small K3 / B1280-full / independent low decoder,
B/B observation-order1: ratio `0.9899870738`, CI
`[0.9803769758,0.9996913743]`. Do not describe every control failure as merely an
inconclusive interval containing equality, or infer its cause from this sample.

The C/WH1 roster has 104/360 intervals per load order not proven faster, all
certified: every certified K3 cell, every certified K16 encoder, and K16 low/
distant decoders at B2-full, B2-tail1 and B64. All 240 small-route comparisons per
load order nominally have upper bounds below one, but control failure prevents
qualification. This still does not address every K, default admission, or the
separate existing certified-path restoration requirement.

Exact Python 3.8/3.12 replay and independent reconstruction passed **before**
advancing HEAD or pinned documents: 777,600 observations, 86,400 retained warmups,
3,600 cell intervals plus 12 primary intervals, all 389 current source/artifact
pins, all 48 codec-manifest and 29 observer-manifest members. Maximum numerical
discrepancy `1.11e-16`. Recorded WORK totals `183.499771298` seconds across both
workers. This is validation of the rejection and its evidence, not a speed gate
pass. No recovery equations or production sources changed.

Bundle `complete.json` SHA256:
`19d7f29303de3f35ee7412b3fa1f64876183d9f7e82e762f67a39f724ab111fa`.
Independent auditor and report:
`/tmp/wh2-dormant-core-independent-audit.PxMjCIlT/{audit.py,report.json}`;
report SHA256:
`fded79c3b13a588ae6726768b5466bdb673246b11aa91fd40bd6a9c050d2b8f6`.
Do not rerun a historical current-pin audit after these report-only edits or
rebind its receipt to a later HEAD. Existing admission/uncovered-K work remains
open; this candidate and cohort are not a new retry opportunity.
