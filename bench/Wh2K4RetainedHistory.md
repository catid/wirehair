# K4 retained timing history: no fresh qualification

Issue `wirehair-sxvz.16.1.20.82.6.1.1.2.2`, audited 2026-09-16.
The previously closed fresh-qualification issue's R22 speed claim was unsupported.
No production/default admission is justified by these records. This audit ran
no codecs, tests of the historical libraries, new timing or old current-source
verifiers. It preserved every original bundle and outcome.

## Complete retained cohort

| Namespace | Producing HEAD | Recorded outcome | Raw observations | Failed A/A |
| --- | --- | --- | ---: | ---: |
| R0 | `6621dbb` | PASS, already rejected for provenance | 77,760 | 0/216 |
| R17 | `1a38c4d` | INVALID: claim bytes, empty raw output | 0 | Not measured |
| R18 | `3cd347e` | PASS | 77,760 | 0/216 |
| R22 | `a802b66` | PASS | 77,760 | 0/216 |
| R23 | `6b2851e` | CONTROL_FAIL | 77,760 | 1/216 |
| R24 | `68196dd` | PASS | 77,760 | 0/216 |

All six bundles still use protocol `wirehair.wh2.k4-serialized-cost-r0`.
Their exact COMPLETE hashes are fixed in `Wh2K4RetainedHistory.py`.
R0 is under `/var/tmp/wh2-k4-serialized-cost-r0`; the others are under
`/var/tmp/wh2-k4-serialized-cost-r1.<namespace>/science`.

The read-only reconstruction checks all bundle member sizes/hashes before and
after traversal, every raw coordinate, ordered WORK clocks, declared own-endpoint
API counts, handle-count records, warmup accounting, per-replicate paired logs,
t11 intervals and unchanged decisions. Total: 388,800 observations, 43,200
warmups and 1,800 separate statistical cells. All 720 treatment constraints
across the five measured runs nominally pass; this does not make any run
promotional. R23's failed control is small-K4 borrowed low decoder, B2-tail1,
observation order0. No control was dropped or reinterpreted.

Every measured run has identical retained fixture bytes, canonical SHA256
`b873efca75d701eef4ca0057a016bccc215b43d48c53b3f79068c4a8b6670955`.
Actual R22/R23/R24 native production and K4 boundary archives also match:
`0d4a0a7ec3cfa4b8e0f08df75e2c12935b51f28653c8c953a32bb49a5658eff5` and
`c475f75d9e183222cf55a9b75a9f672d683d426261b32a9dffb14fd4149f1790`.

The two later code commits change only claim-directory strings in the fresh
adapter and six-line worker. Between R22 and R23 there is also a Beads closure
commit; it changes no codec or measurement code. In particular, the failed
R23 run followed by a path-only R24 rerun violates the frozen no-retry rule.
Renaming a directory is not a new experiment. R24 cannot rescue R23, and an
earlier or later favorable run must not be selected as a replacement.

## R22 provenance findings

Independent source/artifact review, locally reproduced, found:

- The committed `Wh2K4SerializedCostR0.current()` calls the unconditional R0
  producer rejection in `Wh2K4CostBuildR0.verify_qualified_library()`. Both run
  and replay call `current()`. These files match producing commit `a802b66`
  byte-for-byte. The fresh adapter patches its build-time modules, not the
  launch reader. A runtime override or different launcher must therefore have
  been involved, but no authenticated R22 override/launcher or independently
  retained R22 audit report was located. The old located audit concerns R0.
- All three R22 backend receipts omit six actual test dependencies:
  `/usr/include/c++/13/{iostream,thread,cinttypes,bits/std_thread.h,bits/this_thread_sleep.h}`
  and `/usr/include/inttypes.h`. The qualifier covers the 19 production
  translations and serialized boundary translation, not test producer closure.
  For example, native `GF256InPlaceTest.cpp.o.d` lines188–189,235–236,248 and
  `Wh2SmallSerializedTest.cpp.o.d` line255 name these missing inputs.
- The native DSO uses `abi/wirehair.map` (retained native `build.ninja`
  lines292/294), excluded by the qualifier's suffix-based source roster.
  Actual PIE startup inputs `Scrt1.o`, `crtbeginS.o`, `crtendS.o` are also absent;
  the qualifier records their non-PIE counterparts. The retained configure log
  records PIE linking, and native `gf256_inplace_test` is ELF type DYN.

This is not an allegation of modified codec bytes or failed recovery. Positive
evidence remains: all 2,239 R22 receipt files exist; all 546 repository pins
match producing Git; three observer manifests exactly merge into that receipt;
all proof pins survive the handoff; all 57 production archive members match
their objects. All 68 retained neutral tests passed and all 27 recorded commands
exited successfully with empty stderr. All 185 built-object recipes retain
the required backend flags and no LTO, and all 55 retained tool-runtime paths
are pinned by the cost builder. These facts do not fill the omitted input or
launch provenance gaps.

The history auditor intentionally does **not** verify wait/target/phase clocks,
thread-CPU chronology, rusage or CPU identity, certify full producer closure,
rederive the polynomial payload oracle, or establish newly measured performance.
It authenticates the retained numerical outcomes and the path-only retry history.
The separate retained K4 recovery assessment is not changed.

## Preventing accidental reuse

R1's cost build entrypoint now fails before proof reads, imported-module changes
or outputs. Its qualifier likewise refuses to publish another incomplete R1
closure assertion. The R1 C++ wrapper is a compile-time error. Historical source
is still available at its producing commits; no retained artifacts were edited.
These are safeguards, not repairs to historical receipts or permission to retry.

Read-only reproduction, writing only a new external report:
`python3 -B bench/Wh2K4RetainedHistory.py --report <new-external-report.json>`.
Final retained report:
`/tmp/wh2-k4-retained-history.E4krwc2x/report-final.json`, SHA256
`09b36913d0731596f74c08d1a4108890bcb3397c76342a50be96c3ab474d521f`.
Separate Python 3.8 traversal produced identical bytes at
`/tmp/wh2-k4-retained-review.VScFL0bR/report-python3.8.json`.
Review corrected an initially too-strict `<180s` WORK check to the original
inclusive bound and clarified the deliberately limited audit scope. An earlier
history traversal also caught the separate Beads metadata commit between R22
and R23; that metadata is now explicitly distinguished from the two path-only
code changes. These were read-only audit corrections, not timing retries.

The 62-test history/qualifier/old-controller synthetic suite passed on Python
3.12 and 3.8; after adding the compiler-rejection test, all 14 focused history/
qualifier tests passed again on both. No codec or historical worker was run.
Repeated independent source and documentation review is clean.

Remaining work is explicit: any future neutral provenance utility must cover
actual test/link/tool inputs and its real launcher. The spent K4 timing family
must not be renamed or remeasured. Installed/default admission and all-K speed,
recovery and construction-seed objectives remain open.
