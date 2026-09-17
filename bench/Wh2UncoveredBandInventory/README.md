# Uncovered-band recovery inventory

Baseline-only diagnostic, tracked by Beads `wirehair-sxvz.16.1.20.82.3.2`.
No production change, candidate, timing, holdout or all-K claim. The previous
disabled checkpoint is preserved in commit `58db84d`; this controller requires
actual neutral evidence and committed source before admitting its sole launch.

## Frozen scope

Protocol `wirehair.wh2.uncovered-band-inventory-r0`, sole scientific namespace
`/var/tmp/wh2-uncovered-band-inventory-r0`. No retry or overwrite.

K7/9/12/16, widths 2/64/1280, 64 roots (first 16 SHA256 hex digits of the
protocol plus `:inventory/` plus decimal indices 0 through 63). Each trace
delivers exactly K+4 distinct IDs with a 65,536-candidate bound. Schedules are
IID 10%, burst 50% with length 8, adversarial 50% at `UINT32_MAX-2*candidate`,
and repair-only 50% at `K+candidate`. The frozen SplitMix64 trace implementation
is covered by an independently written full-roster reconstruction test.

There are 3,072 full-tail loss cases, plus 24 separate full/tail1 low-repair
cases at IDs K through 2K+3. Widths share roots, not independent seed samples.
The low-repair cases are diagnostics, not an exhaustive hard-pattern gate.

Both actual ordinary WH2 and owned WH1 use the unchanged baseline DSOs under
`/tmp/wh2-small-unit-diagonal.MfPjFJaS`, never the rejected candidate library.
Basis messages observe each generator row; independent polynomial GF(256)
arithmetic checks every packet and information rank. Actual first-success
endpoints and two exact guarded recoveries define decode results; full-rank
decoder lag is reported separately. Decoder creation follows encoder teardown.

Targets are ordered descending by `(WH2_OH0_failures - WH1_OH0_failures,
WH2_OH4_failures, WH2_OH0_failures, -K)`. Recommend the first target only if
its excess is positive and its WH2 OH0 failure rate exceeds 1%. This selects
a dimension for later work; it cannot establish candidate superiority.

## Admission and audit

Synthetic checks run without native codec calls:

```bash
python3 -B -m unittest discover -s bench/Wh2UncoveredBandInventory -p 'test_*.py'
```

After source review, `Control.py neutral native` and `neutral portable` each
run 16 systematic-prefix cases at separate widths 17/65, full/tail1, under
the same bounded capture mechanism. Both decoders must succeed at K. Neutral
claims bind all four runtime modules, the producing interpreter, library,
limits and environment before/after execution. Every normalized descriptor,
coefficient, packet, endpoint and API-ledger record must match across backends.
The scientific producing interpreter must match both neutral producers.
Failed or partial neutral outputs remain retained and cannot be overwritten.

Limits: CPU 120s, wall 150s, address space 512 MiB, file/stdout 128 MiB,
stderr 64 KiB, core 0. Subprocesses receive an explicit clean environment;
timeouts/output overflow kill the process group without retry.

`Control.py preflight` verifies both neutral proofs, the retained 1,039-member
engineering manifest, and that this harness matches committed source.
Review/freeze/push before `Control.py run`. The run exclusively creates its
namespace and seals either complete captures/analysis or a spent failure.
`Control.py replay` checks complete retained evidence without codec execution;
an alternate reader Python is allowed while the producing interpreter remains
pinned. Perform exact and independent full audits before advancing pinned
HEAD or documents. The DSO checks authenticate bytes and public-function
ownership, not full transitive linker/toolchain/runtime provenance.

## Completed baseline diagnostic — 2026-09-17

The sole run at `cf5f416e44f109f22c5ea1c5b332b6646fc5236a` completed with exit 0,
empty stderr and 86,272,424 captured bytes in 3.973 seconds. That wall time
describes the diagnostic process, **not codec speed**. All 3,096 cases were
retained and audited; the namespace is permanently spent.

Zero-overhead failures among 768 loss cases per K:

| K | WH2 | WH1 | WH2 failures at +1 | WH1 failures at +1 |
|---|---:|---:|---:|---:|
| 7 | 9 (1.17%) | 33 (4.30%) | 0 | 0 |
| 9 | 7 (0.91%) | 23 (2.99%) | 0 | 2 |
| 12 | 49 (6.38%) | 11 (1.43%) | 3 | 0 |
| 16 | 1 (0.13%) | 24 (3.13%) | 0 | 0 |

Every loss case recovered by +2. Both arms also passed all 24 separate low-repair
cases at zero overhead. Actual decoder endpoints matched information-rank
endpoints throughout; there was no full-rank decoder lag in this inventory.

The frozen priority is **K12, K9, K16, K7**, with **K12 recommended**. Its worst
cell is B1280/adversarial: WH2 33/64 zero-overhead failures versus WH1 0/64.
All 49 observed high-ID coefficient rows in that cell are distinct and their
union has rank 12; the 64 zero-overhead prefixes have ranks 11 (33 cases) or
12 (31 cases), with no duplicate coefficient rows within a prefix. This points
to linear dependence, not exact row collisions. These are retained-sample
diagnostics, not universal recovery-rate or construction-seed certification.

The both-arm K12 history projection retains 63 failed-prefix origins (60 at
length 12, three at length 13), 60 distinct ID prefixes and 64 inventory roots.
Do not shorten failed prefixes or treat their later successful suffixes as new
tests. Follow-up `wirehair-sxvz.16.1.20.82.3.3` will freeze and screen a finite
K12 GF(256) construction before any native integration or lifecycle timing.

Evidence anchors:

- Bundle: `/var/tmp/wh2-uncovered-band-inventory-r0`;
  `complete.json` SHA256
  `16fcf13214cd25362fe66ee35f62c1c9616ed082e343fa5e2009d81d61359ea0`.
- Neutral: `/var/tmp/wh2-uncovered-band-qualification.C0lOL3VL`;
  16 cases per backend, identical normalized parity
  `92de12b58a8731bb50a938fffce23d0c0feef3dbabf331f5f3905c5fe033315b`;
  31 synthetic tests pass under each Python version.
- Exact Python 3.12/3.8 replays:
  `/tmp/wh2-uncovered-band-root-audit.6uwksoXP`.
- Independent no-harness-import audits:
  `/tmp/wh2-uncovered-band-independent-audit.kw63F6qe`;
  both report SHA256
  `21cd4ed9c81d94d7aa1e9a912a2ec3130139aa0624ba605dd81791260b16012a`.
  All 92,880 packet payloads (41,672,160 bytes), rank/actual decode endpoints,
  API ledgers, 48 cells, four totals, neutral evidence and 1,062 input pins were
  independently verified before HEAD or pinned documentation advanced.
- K12 derived history/row checks in the same independent-audit directory:
  `k12-projection-python312.json` and `k12-projection-python38.json`, identical
  SHA256 `0798dc50291a1b5db9d6cd3538fea90611b6a610bd21a0853f2d32017638f4a7`.

Strict current-tree replay is now historical: use the producing commit and
unchanged pinned inputs, not this advanced README, for that check. No candidate,
speed improvement, production change or promotion resulted from this inventory.
