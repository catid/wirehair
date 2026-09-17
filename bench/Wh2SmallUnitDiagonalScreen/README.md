# Unit-diagonal decoder lifecycle screen

Prospectively frozen in `wirehair-sxvz.16.1.20.86.1`. This measures only the
two arithmetic substitutions in `bench/Wh2SmallUnitDiagonal`. Production,
equations, allocation counts and layouts are unchanged. No recovery-rate or
all-K qualification follows from this bounded screen.

The sole screen is now terminal **CONTROL_FAIL; candidate not retained**.
The independently audited result is below. Do not repeat the launch command.

The workload, statistics, launcher, generator and synthetic tests are unchanged
from the reviewed `Wh2SmallDecoderStorageScreen` observer. This is a distinct
arithmetic candidate, not a retry or combination of that rejected allocation
candidate. Historical namespaces and results remain untouched.

## Immutable inputs

Qualified codec root: `/tmp/wh2-small-unit-diagonal.MfPjFJaS/`.
Its 1,039-member `NEUTRAL.sha256` records build/generated artifacts, logs,
compiler-reported source/header dependencies and relevant build/test inputs.
All entries were independently rehashed after sealing. It is not a complete
transitive compiler/linker/runtime provenance seal.

| Input | SHA-256 |
| --- | --- |
| Codec neutral manifest | `93bfff5491f451c15ff8be50ec4453bef8770c0dbbc24e290ad48caa1822f532` |
| Native baseline DSO | `bd3c353847ec2838d95d5b73804664d3a2df8ac76b22ce4f5cf8c27fd58c1bbe` |
| Native candidate DSO | `af37e1286af35fe7695e9d423aafd3d4cba2ea7127ef432ce89f807ef990d0ba` |

Observer-only builds use
`/tmp/wh2-small-unit-diagonal-screen.tSZm35BW/{native,portable,asan}`.
They load the corresponding existing codec libraries; no codec is rebuilt.
Portable selects non-GFNI arithmetic on this host, not a different machine.
The sanitizer build instruments observer and libraries and rejects timing.
Do not rebuild or test into sealed directories.

All six observer neutral CTests passed (both loading orders in each of the
three builds), as did all 24 synthetic tests under both Python 3.12 and 3.8.
The 270-member observer manifest SHA-256 is
`ca5fa7ff12828b0c52f03470b95cf4a9914950b03ac0923fcb934848511fc3e3`.
It excludes the controller and this document, which the launch claim pins
separately. Generated worker SHA-256 is
`d395c62ecac4819adfc6125aa50f0d7efb997d5a3d0e43e7533e07a3b0608565`;
native/portable executable SHA-256 is
`b63fcbc9cece6c072b76a3f06e2994d2a6e4ea72714a807b12e06e9a01482ff5`.
The worker source and native binary exactly match the reviewed observer;
only its input libraries and separate scientific namespace change.

## Frozen workload and gates

Seven routes: ordinary WHV2 K3, explicit small WHV2 K3/K5/K8, explicit
certified WHV2 K3/K16, and the separate WHK3 API. Each uses B2/full, B2/tail1,
B64/full, B1280/full and B1280/tail1, with independent and borrowed ownership:
70 fixtures and 210 metric cells.

The three metrics are complete encoder construction + 3K packets + destruction,
complete low-repair decoding, and complete distant-repair decoding. Each
decoder includes construction, feeds to its verified first-success endpoint,
recovery and destruction. WH1 has its own packets and first-success endpoint.
All packet and recovered bytes, padding and guards are checked. Public source,
output, descriptor and repair-array addresses match within every cell/replicate.

Each of two DSO loading orders uses CPU50, batch64, 12 replicates, two
observation orders, five comparisons (BB/CC/WH1WH1/CB/CWH1), and 18 positions
per panel: two retained warmups and eight adjacent pairs. All 453,600 rows per
loading order are retained (907,200 total, including 100,800 warmups).

Paired log ratios are averaged within replicates. Two-sided t11 95% intervals
use the 12 replicate means. Candidate/reference time ratios below one favor
the candidate. Each loading order must independently pass:

- All 1,260 same-code intervals strictly inside `[1/1.02, 1.02]`.
- All 420 candidate/baseline upper bounds strictly below 1.02.
- All 300 affected-small-route candidate/WH1 upper bounds strictly below one.
  The other 120 certified WH1 bounds are reported; their existing deficits
  remain unresolved. Certified routes still must pass retention.
- Four decoder primary upper bounds strictly below one: low/distant by
  observation order, averaging all 50 affected fixtures within each replicate
  before constructing intervals to preserve covariance. Encoders must retain
  performance but are not decoder primaries.

Failed same-code gates yield CONTROL_FAIL, overriding nominal gains. Otherwise
any failed retention, required WH1 or primary gate yields FAIL. Both loading
orders must pass; no pooling, subsets, reweighting or simultaneous-confidence
claim. A passing aggregate is not a claim of improvement in every cell.

## One-shot execution and audit

The sole namespace is `/var/tmp/wh2-small-unit-diagonal-screen-r0`.
Existence permanently spends it. No retry, renaming, tuning, trimming or rescue.
Finish both loading orders after statistical failure, but stop for invalid
infrastructure, identities or payloads. Before timing: neutral qualification,
independent source review, artifact sealing and exact source commit/push.

`Launch.py` supplies a clean environment, a fresh process group, CPU250s,
address-space384MiB, core0 and file-size128MiB limits. Capture limits are
stdout128MiB, stderr64KiB and wall270s. Failure kills the process group and
preserves bounded prefixes. The worker's 240s cap includes preparation,
checking and publication; summed measured WORK including warmups is also
bounded by 240s. Nonzero exit, stderr, missing rows or mismatched inputs fail.

`Screen.py` verifies the codec and observer seals, pins committed source and
build/test inputs in `claim.json`, and rechecks all pins before/after each
worker. Exact replay and a separate complete raw/statistical/provenance audit
must finish before HEAD or result documentation advances. Only a passing
screen permits further consideration of production adoption.

Synthetic checks: `python3 -B -m unittest discover -s bench/Wh2SmallUnitDiagonalScreen -v`
and the equivalent Python3.8 command. Fresh observer qualification builds use
this CMake directory, the matching `DSO_DIR`, and `SANITIZE=ON` only for ASan.
The sole native launch, once all prerequisites are sealed and reviewed, is:

```sh
python3 -B bench/Wh2SmallUnitDiagonalScreen/Screen.py --run /tmp/wh2-small-unit-diagonal-screen.tSZm35BW/native
```

## Terminal result, 2026-09-17

Producing commit `f4723a150a54c4dff8cae9cebfad8a2eae9caa54` was pushed before
the sole launch. Both workers exited zero with empty stderr, complete captures
and unchanged inputs. Wall times were 177.315/186.703 seconds; measured WORK
including warmups was 95.914/96.322 seconds. No codec rebuild, timing retry,
trimming, retuning or subset selection occurred.

| DSO loading order | Same-code failures / 1,260 | Retention failures / 420 | Required WH1 failures / 300 | Decoder primary passes / 4 | Decision |
| --- | ---: | ---: | ---: | ---: | --- |
| Baseline first | 14 | 95 | 0 | 4 | CONTROL_FAIL |
| Candidate first | 15 | 87 | 0 | 4 | CONTROL_FAIL |

The eight decoder primary point estimates show nominally 1.91–2.51% less
time, but failed controls override those gains. All 29 failed same-code
intervals include one while extending beyond the equivalence bounds. The
retention counts are failed upper-bound gates, not 95/87 individually proven
slowdowns greater than 2%. Both encoder and unchanged certified routes have
nominal regressions; their cause is not established by this experiment.
All required small-route WH1 bounds passed, but 104 of the other 120 certified
WH1 bounds failed in each loading order. Those deficits remain unresolved.

Before any producing HEAD or pinned documentation advanced, exact replay
passed under Python 3.12 and 3.8 and the complete current input map matched all
1,383 claim pins and the producing HEAD. A separate implementation importing
neither controller nor launcher reconstructed all 907,200 rows, 100,800
warmups, 58,060,800 lifecycle invocations, 4,200 cell intervals and eight
decoder primary intervals. It checked chronology, all fixtures and endpoints,
public addresses, covariance, strict gates, captures, seals and unchanged
before/after source pins. Independent reports under both Python versions are
byte-identical; maximum numerical discrepancy from the controller's different
log-ratio arithmetic is `1.0658141036401503e-14`, with identical decisions.
Repeated root and independent source/artifact review found no outstanding bugs.

Bundle: `/var/tmp/wh2-small-unit-diagonal-screen-r0/`.
Exact replay/input-map logs:
`/tmp/wh2-small-unit-diagonal-root-audit.Eo3lToxJ/`.
Independent script, reports and logs:
`/tmp/wh2-small-unit-diagonal-independent-audit.wS4UtJzx/`.

| Artifact | SHA-256 |
| --- | --- |
| `complete.json` | `875508969c502d7f1cc373b5ca8dd7bc084e29afe0058e0b4c05883d3e836ac3` |
| `normal.csv` | `dd2f6d27bb3b65480dd63f839cc80fc6ebfd78f85f9d984de7891a969e3e5d4a` |
| `reverse.csv` | `4068bb1fc147df4cf27de071768bcc65aaa70aeb627db30ce4574ebe6d5d653d` |
| Independent audit script | `235c44e0e365545e96fea3155d9d09a8e805b5f40c7e2c5b316ece248668a19e` |
| Either independent report | `aac3660d77318fa3a29402dd6128b7c6cd50fa1071d5053c0eeedd16bb8f7981` |

This result section intentionally changes a document pinned by the producing
claim. Strict current-tree provenance checks are now historical: use the
recorded revision and completed audits, not a fresh current-tree revalidation
after documentation advancement. Mathematical replay remains available.

The experiment is rejected, not promoted. Its proven arithmetic reduction
does not establish an acceptable lifecycle speedup. No retry, alignment or
layout rescue, subset promotion or combination with rejected candidates is
authorized by these results. The broader speed/recovery objective remains open.
