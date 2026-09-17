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
