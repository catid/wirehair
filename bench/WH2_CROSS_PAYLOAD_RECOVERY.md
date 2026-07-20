# WH2 cross-payload recovery protocol

`wh2_cross_payload_recovery.py` compares two raw WH2 architectures after the
separate solve-timing experiment has chosen which pair deserves the expensive
recovery campaign.  It does not select an architecture, apply seed fixes, make
a speed claim, or start benchmark work during preparation.

## Why the extra campaign is required

`MatrixSeedFromProfile` includes `BlockBytes`, so changing the payload width
changes the solve graph even when every explicit precode option is identical.
The tracked degree-balanced all-K result covered 575,991 cells per arm and
reduced failures from 733 to 711, but its graph domain was only
`BlockBytes=64`.  Its `full_payload_solve=0` setting makes that excellent rank
evidence, not evidence that the same weak K or recovery delta persists at
other widths.  In particular, the 1280- and 4096-byte timing rows cannot close
this recovery gap.

This protocol therefore uses the natural graph seeds at block sizes 64, 256,
1280, and 4096.  `--seed-block-bytes` is forbidden.  A K-indexed
`--packet-peel-seed-table` (and any architecture option named as a fixup or K
fix) is also forbidden: architecture selection must compare unpatched designs.
Weak-K counts, multiplicity, repairs, introductions, and work are reported as
raw metrics.  Seed repair is a later stage only after the architecture is
selected.

## Leakage boundary

Preparation freezes two disjoint seed namespaces.  Discovery scans every
`K=2..64000` at all four widths with three development roots, one exact-K
trial, 50% loss, and `iid`, `burst`, `adversarial`, and `repair-only`
schedules.  The candidate cohort is exactly:

1. every K with at least one rank failure in either arm, at any discovery
   width, root, or schedule; plus
2. controls declared in the controller source before any discovery output is
   read, including boundary/scale K and neighborhoods representing 100 KiB and
   1 MiB messages at every width.

`seal-cohort` reads only the complete hash-bound discovery receipt cube.  It
stages the deterministic evaluation ledger, publishes the cohort, and only
then atomically publishes the evaluation directory.  An interrupted or
repeated invocation may finish that exact transaction, but it refuses any
cohort or ledger that differs from a fresh discovery reduction.  Evaluation
roots were already committed in `config.json` but are never inputs to cohort
selection.

Evaluation uses four fresh roots, the same four schedules, losses 0.10, 0.35,
and 0.50, and 64 trials per root (256 trials per aggregate cell).  Each
invocation requests OH0/OH1/OH2 with
`--paired-overhead-stream`, so the larger overheads append packets to the
exact smaller prefix.  Control and candidate are adjacent jobs with identical
K, width, root, schedule, loss, and trial axes.  Reduction reconstructs the
per-trial paired quadrants as well as aggregate, width, schedule, loss,
overhead, root-to-root, weak-K, work, repair, and introduction totals.  XOR and
field-operation work sums are retained in each one-dimensional axis and each
width/schedule/loss/overhead stratum so reliability changes can be compared
with their decode-work cost.

The recovery gate requires zero codec errors and no aggregate failure increase
at any of OH0, OH1, or OH2.  Per-width/schedule/loss deltas and weak-seed
introductions remain report-only during raw architecture selection, matching
the project methodology.  A passing recovery result still records
`promotion_ready=false`: independent full-payload timing evidence is required.

## Frozen input and commands

The canonical input JSON has schema
`wirehair.wh2.cross_payload_recovery.v1.input`, a campaign label, and exactly
two ordered arms named `control` and `candidate`.  Each arm supplies an
executable, its full source commit/tree, a hash-bound build receipt, a label,
and one explicit mixed-completion `architecture_argv` list.  The controller
copies those artifacts, probes each copied binary once, and freezes the full
effective static preamble before generating any task.

```text
python3 bench/wh2_cross_payload_recovery.py prepare \
  --root /tmp/wh2-cross-payload-v1 --design /path/to/canonical-input.json

/tmp/wh2-cross-payload-v1/frozen/wh2_cross_payload_recovery.py run \
  --root /tmp/wh2-cross-payload-v1 --phase discovery --workers 128

/tmp/wh2-cross-payload-v1/frozen/wh2_cross_payload_recovery.py run \
  --root /tmp/wh2-cross-payload-v1 --phase discovery --workers 128 --execute

/tmp/wh2-cross-payload-v1/frozen/wh2_cross_payload_recovery.py seal-cohort \
  --root /tmp/wh2-cross-payload-v1

/tmp/wh2-cross-payload-v1/frozen/wh2_cross_payload_recovery.py run \
  --root /tmp/wh2-cross-payload-v1 --phase evaluation --workers 128 --execute

/tmp/wh2-cross-payload-v1/frozen/wh2_cross_payload_recovery.py reduce \
  --root /tmp/wh2-cross-payload-v1
```

The first `run` command is a receipt preflight; only `--execute` starts codec
work.  The controller intentionally does not own or signal the machine's sole
CPU/DIMM/EDAC sampler.  Run execution under the existing externally managed
thermal monitor.  Timing fields remain in hash-bound raw CSV but are excluded
from the semantic recovery hash and cannot support a speed claim.
