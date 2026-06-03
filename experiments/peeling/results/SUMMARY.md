# Peeling Sweep Summary

These runs measure only the peel graph and inactivation schedule. They do not
model dense solver cost or success probability after the residual matrix is
formed.

All rows report statistical residual unpeeled columns/rows, component shape,
and peel-side runtime. `N` is the number of blocks; total payload bytes are
`N * block_bytes`.

## Main Loss Sweep

Input: `--source loss --loss 0.10 --N 128,2048,32000 --trials 40`

At `N=32000`, the best residual column counts were:

| method | mean residual cols | mean residual rows | mean total us | cost |
| --- | ---: | ---: | ---: | --- |
| rqd2_minfill | 193.30 | 193.30 | 49279.9 | O(rows*degree + cols + rowrefs) |
| raptorq_d2cc | 196.55 | 196.55 | 43559.3 | O(rows*degree + cols) |
| rqd2_default | 196.68 | 196.68 | 44287.8 | O(rows*degree + cols) |
| topk8_lookfill | 205.80 | 205.80 | 19885.7 | O(cols + 8*rowrefs) |
| default | 210.57 | 210.57 | 21628.8 | O(cols) |

Avoid the Markowitz-style min-fill family for this peel-only objective: in this
sweep it was both slower and left more residuals than default.

## Focused Loss Sweep

Input: `--source loss --loss 0.10 --N 128,2048,32000 --trials 60`

The focused rerun favored the plain RaptorQ degree-2 component variants over
the min-fill tie-break at large `N`:

| method | N=32000 mean residual cols | mean total us |
| --- | ---: | ---: |
| raptorq_d2cc | 194.77 | 43510.0 |
| rqd2_default | 195.63 | 43132.3 |
| rqd2_minfill | 196.97 | 51595.2 |
| topk16_lookfill | 209.07 | 24545.1 |
| default | 213.08 | 22773.9 |

## Source Shape Checks

Input: `--N 128,2048,32000 --trials 40` with `systematic` and `repair` sources.

At `N=32000`, Raptor/RaptorQ-style degree-2 component selection improved the
residual versus default for both source shapes:

| source | method | mean residual cols | mean total us |
| --- | --- | ---: | ---: |
| systematic | default | 207.00 | 22885.2 |
| systematic | rqd2_default | 191.00 | 41686.7 |
| systematic | rqd2_minfill | 191.00 | 49223.8 |
| repair | default | 233.50 | 26754.6 |
| repair | rqd2_default | 211.53 | 61686.1 |
| repair | rqd2_minfill | 212.10 | 67135.1 |

## Literature Cross-Check

RFC 6330's RaptorQ example decoder first chooses a row of minimum reduced
degree. When that minimum degree is two, it builds the graph whose nodes are
columns and whose edges are degree-2 rows, then chooses a row in a maximum-size
component. The `raptorq_d2cc`, `rqd2_default`, and `rqd2_minfill` methods are
peel-only column-selection approximations of that schedule:

https://www.rfc-editor.org/rfc/rfc6330#section-5.4.2.2

Recent LT/Raptor inactivation papers mainly analyze or optimize the output
degree distribution and often use random inactivation as the schedule model.
That suggests future work on Wirehair's row distribution or permanent
inactivation/precode structure, but it did not reveal a clearly better fast
local greedy schedule for this peel-only harness:

- https://arxiv.org/abs/1408.2660
- https://arxiv.org/abs/1510.08364
- https://arxiv.org/abs/1706.05814

## Current Recommendation

`rqd2_default` is the strongest practical global candidate from these runs. It
gives most of the residual reduction seen from the Raptor/RaptorQ degree-2
component schedule without the extra min-fill tie-break cost. `raptorq_d2cc`
is effectively tied if optimizing only residual count, while `topk8_lookfill`
is the cheaper low-risk ablation if the component scan is considered too
expensive.
