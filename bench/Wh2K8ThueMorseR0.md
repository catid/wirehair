# K8 noncommuting GF256 short screen

Protocol `wirehair.wh2.k8-thue-morse-r0`; sole namespace
`/var/tmp/wh2-k8-thue-morse-r0`. Mathematical feasibility only: no native
codec, timing, comparative WH1 rate, new profile, defaults or all-K claim.
The prospective freeze is recorded in `wirehair-sxvz.16.1.20.82.5` before
candidate arithmetic or implementation edits.

Use GF(256), polynomial `0x14d`, and dimension-eight companion matrices.
Base feedback is the coefficient vector of `product(j=0..7)(x + 2^j)`
in that field, excluding its monic leading term. The second feedback differs only
in its constant coefficient, XOR lambda. In ascending lambda 1..255,
excluding the base constant, select the FIRST pair passing all 495
eight-of-twelve minors for each of the ten fixed length-four Thue-Morse
factors: 4,950 local tests. Selection sees only local algebra, never
historical or loss results. No survivor means EXHAUSTED; any subsequent
failure ends the experiment without reselection or retuning.

Require invertible, noncommuting matrices and the first eight rows
systematic. The existing dimension-parametric mapper stores exactly
65,536 lookup bytes: `2*1024*8 + 4*128*64 + 256*64`. Verify each accessed
row against full dyadic prefix products, and IDs 0..2048 against literal
sequential multiplication. No production core is modified or presumed
to support K8.

Before fresh testing, require all 14,850 minors in 30 twelve-ID seam
windows: starts `2^e-4`, e=3..31, plus `UINT32_MAX-11`. Require OH0 rank
eight in all 72 frozen hard traces (six main training/validation roots,
three widths, four schedules), and full rank in all 42 original-length,
width-preserving historical failure prefixes.

Authenticate the complete original 3,096-case inventory, including all
774 K8 cases and all five sealed manifest members. Preserve 44 failure
origins from both old WH2 and WH1, at every overhead 0..4, with their
original widths/tails, and all 64 inventory roots.

| Input | SHA256 |
| --- | --- |
| Inventory COMPLETE | `93a9068986dc24e1caa2931f5be897022e34f13b8f5ca1170a294d87365c5f2c` |
| 44 origins | `392b45ed4aba311d43fe4246db2c85053719186f10d0a7718c113d14558407a5` |
| 42 width-preserving prefixes | `aa6e36367301f687f6606e7d58bbb83c20fbdedb9e571c84688e21917d9a4fba` |
| 64 roots | `6a03f629cae261c18128ba4707b2656dda7833a9f836948442f9fb92e8c20bbf` |

Only after those gates pass, evaluate 512 fresh roots from the first
16 SHA256 hex digits of protocol + `:fresh/` + decimal(i), i=0..511.
Exclude inventory/main roots and prior K2/K3/K5/K6 fresh roots. A collision
is INVALID, never replaced. The unchanged main loss law uses K8,
B2/64/1280, IID10/burst50/adversarial50/repair-only50, twelve delivered
IDs and a 68,608 candidate-ID cap. Preserve every OH0..4 rank. Require
at most 1% OH0 failures overall and in each of twelve 512-trace cells.

Reuse the existing one-shot controller unchanged: worker wall60 seconds,
AS512 MiB/core0, outer70 seconds, raw4 MiB/stderr1 MiB/aggregate8 MiB.
Neutral tests use unrelated matrices and synthetic outcomes, never score
the frozen pair. Commit/push before receipt and the sole scientific launch.
Before advancing HEAD/report, audit exact source plus independently
implemented field arithmetic, first-local selection, minors, lookup, every
trace/rank, historical projection and complete provenance. A passing
screen permits native feasibility work, not speed or production promotion.

## Retained result

The sole screen is **PASS**, produced at pushed source
`8ac0848066f21c32da845e4c3ca440037d406628`. Ascending local-only selection
rejected lambda 1 at its 311th minor and selected lambda 2. The feedbacks
are `(96,19,186,153,85,252,7,255)` and `(98,19,186,153,85,252,7,255)`.
No history or loss outcome influenced that selection.

| Frozen check | Result |
| --- | --- |
| Selected local minors | 4,950/4,950 full rank |
| Seam minors | 14,850/14,850 full rank |
| Hard traces, OH0 | 72/72 full rank |
| Original-length historical prefixes | 42/42 full rank, all 44 origins preserved |
| Fresh traces, OH0 | 12 failures / 6,144 = 0.1953125% |
| Worst width/schedule cell, OH0 | 3/512 = 0.5859375%, B2/adversarial |
| Fresh traces, OH1 through OH4 | 0 failures at every overhead |
| Packed lookup | 65,536 bytes |

This is not a paired recovery-rate comparison against WH1 or old WH2:
the baseline inventory and fresh screen use different root cohorts.
The result qualifies the structural candidate for native work only.

Exact-helper retained replays, an independent auditor importing no experiment
code, and retained-trace comparison against the native frozen generator all
passed under Python 3.12 and 3.8 before producing HEAD advanced. The
independent audit regenerated the complete 3,096-case inventory, checked
both candidate selection records, all local/seam minors, every lookup byte,
2,049 literal sequential rows, all 2,347 accessed rows, and all 31,080 hard/fresh
prefix ranks. Its two interpreter reports are byte-identical. Both source
review passes were clean; all 34 neutral/shared tests passed per interpreter.

Worker wall time was 2.606742502 seconds; observed whole-controller time was
2.936050238 seconds. These are mathematical experiment execution times,
**not encoder or decoder timing**.

Retained artifacts:

- Bundle `/var/tmp/wh2-k8-thue-morse-r0`; COMPLETE SHA256
  `4911bedbb288c8a7e39577c35dbe6c31a0a2990e0d4dfab1ee4466a32f99fdb0`.
- Raw SHA256
  `f8b976f155fa755d01d8184df3ee8c911f56bcd96f00808bf41fc12fda5f9345`.
- Lookup SHA256
  `512c6646e44517964e7e6a7cd0ffa41057802182ddc500a818af541c05770817`.
- Receipt and exact/native audits `/tmp/wh2-k8-thue-receipt.so3PSauw`;
  receipt SHA256 `db7daa5b8f7bed6e6117e3da617b58fb26f4e10dc0307f3e28e568891fff313f`.
- Audit manifest `AUDITS.json` SHA256
  `a6b889fa739dd84ad3c7d0bc209d1263fc89834b9fabbac7bb42127776c6cdaf`.
- Independent auditor `/tmp/wh2-k8-independent-audit.PACGdwEN/audit.py`, SHA256
  `334370d00be177593ee463652ff02282a4c539f9fc59e5a4356d54e71d10e6e0`.
- Independent `report-python312.json` and `report-python38.json` both SHA256
  `20792372c1d71901f3531a0d2711ca7a392d52bb8dc014e7515f1cc7a14778ab`.

Native qualification is tracked in `wirehair-sxvz.16.1.20.82.5.1`.
No production/profile/default code was changed, and no native speed,
actual-codec comparative recovery, zero-bad-seed all-K validation or
promotion is claimed.
