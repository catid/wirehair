# Two-byte small-profile payload screen

Benchmark-only candidate for `wirehair-0wio`, based on production `68196dd`.
There is no production, lookup, equation, default-selection or recovery-policy
change. For K3/K5/K8 with B2, inline the two GF(256) dot products instead of
entering the generic bulk kernel. Other dimensions, wider payloads and WH_COUNT
delegate unchanged. A single helper replaces three repeated specializations.

Early exploratory programs from this session are **not speed evidence**:
the direct kernel allowed invariant hoisting, public loops were unpinned with
unmatched build configurations, and their checks/rosters were incomplete.
Do not use their approximately 30–40% estimates to promote this candidate.

## Neutral qualification

Use a fresh external build with `cmake -S bench/Wh2SmallPayload2 -B <build>
-G Ninja -DCMAKE_BUILD_TYPE=Release`, then build and run CTest. Separate builds
use `-DPORTABLE=ON`, `-DSANITIZE=ON`, or `-DCOUNT=ON`. Sanitizer tests use
`ASAN_OPTIONS=detect_leaks=1:detect_stack_use_after_return=1` and
`UBSAN_OPTIONS=halt_on_error=1`. The two shared libraries reuse identical common
objects; only the four small-profile translation units are rebuilt against the
external modified header. Default project builds never include this experiment.

The polynomial oracle exhausts 256 input bytes × 256 coefficients at each
source position for K2/3/4/5/6/8. It also checks 208 mixed/fallback shapes per K,
unaligned exact-sized allocations, zero/unit coefficients, repeated source
pointers and nonpositive no-ops. Existing WHV2 tests cover both ownership
policies, partial tails, failure paths and allocation injection. Public shared
fixtures check matching baseline/candidate packets and first-success recovery.

## Prospective diagnostic protocol

After neutral tests pass, run `python3 bench/Wh2SmallPayload2/Screen.py --run
<native-build>`. Sole namespace: `/var/tmp/wh2-small-payload2-screen-r0`. No
retry, trimming, resizing, load-order pooling or promotion from failed controls.

Roster: K3/K5/K8 × (B2 full, B2 one-byte tail, B64 full, B1280 full) ×
independent/borrowed × prebuilt repairs/full encoder/low-ID decoder/distant-ID
decoder. K5/K8 use explicit installed profiles; **not ordinary defaults**.
WH1 independent matches the owning constructor and borrowed matches zero-copy.
Each decoder pays through its own first success. Encoder lifecycles include
creation, K systematic, K low repairs, K distant repairs and destruction.
Prebuilt-repair cost excludes creation and emits K low plus K distant repairs.

Two separate processes reverse DSO load order while maintaining arm identity.
CPU50, baseline/candidate/WH1 A/A controls plus C/B and C/WH1 comparisons,
12 replicates, both side orders, 18 observations/panel, 32 cycles/observation.
First two observations per panel are retained warmups; the following eight
paired log-ratios form each replicate. Every output is checked outside WORK;
rows are retained in memory until publication. Compiler cannot hoist DSO calls.
An explicit caught failure publishes its completed prefix and cannot pass.
Worker cap 180 seconds per process; controller timeout 240 seconds per process.

All A/A 95% t11 intervals must lie strictly inside [1/1.02,1.02]. All B2
prebuilt-repair and encoder C/B upper bounds must be <1. Other C/B cells must
be <1.02 (a **screening envelope**, not a zero-regression claim). C/WH1 upper
bounds are separately reported and never repair an A/A failure. Both load
orders must pass to justify further qualification. No 5% minimum gain.

This is an early mechanism screen, not a formal installed-package/production
promotion gate. It records current source, generated headers, binaries, build
recipes and raw hashes but does not implement the older campaigns' complete
transitive toolchain/runtime receipts. No all-K, recovery-rate, universal speed
or preserved-certified-path claim follows even if it passes.

## Terminal result (2026-09-16): not qualified

The sole namespace is spent and sealed. Both workers completed with empty
stderr and 207,360 rows each. Do not rerun it or use descriptive treatment
gains to override failed controls.

| DSO load order | Decision | Failed A/A controls | Failed C/B constraints | C/WH1 cells not proven faster |
| --- | --- | ---: | ---: | ---: |
| Baseline first | CONTROL_FAIL | 47/576 | 11/192 | 40/192 |
| Candidate first | CONTROL_FAIL | 55/576 | 15/192 | 40/192 |

All B2 encoder/repair treatment constraints nominally passed. Descriptive B2
full-encoder time reductions were approximately 14–17% for K3, 22–26% for K5,
and 15–16% for K8. However, K5/B64 independent full encoding took approximately
3–4% longer in all four load/side-order combinations; some decoder retention
constraints also failed. None of these estimates is qualified speed evidence.
Production remains unchanged; this did not measure improved recovery.

Native, portable, ASan/UBSan, and WH_COUNT builds each passed 10/10 neutral
CTest gates. The controller's seven tests and exact bundle replay passed on
Python 3.8 and 3.12. An independent audit (without importing `Screen.py`)
verified both chronologies, 11,520 panels and all 960 statistical decisions per
load order, all seven manifest hashes, and 286 claim pins at `7f27c92`.
Maximum ratio/interval discrepancy was `6.66e-16`. Failed A/A intervals excluded
1 in 22 normal and 24 reverse cells, so failure was not solely imprecision.

Bundle `/var/tmp/wh2-small-payload2-screen-r0`:

- `complete.json`: `8d0a34ee7971740b135bbee1af988d458de8108b5a9b81248083f147890b7a33`
- `run.csv`: `d3627c4a21e74c09211bc0de924e88fb444463ac428f81c68e23847d5127dc13`
- `run-reverse.csv`: `fb577d9276680b9bfde2df73ac90f7852c791f22e3c424d5525baafcfa9eff9c`

Static inspection of the recorded native libraries found `EncodeSmall` grew
from 5,468 to 6,298 bytes; `DecodeSmall` retained its size but shifted by 832
bytes. An extra B2 branch now precedes each wider repair's generic-kernel call.
These are leads, not established causes of the regression. Follow-up
`wirehair-6llt` inspects a separate out-of-line dispatch without rerunning or
reinterpreting this failed candidate.
