# Alignment-only candidate: prospective two-layer screen

Issue `wirehair-hgiu`. Reuses the exact neutral-qualified libraries in
`Wh2AlignedSmallBasis`; this project only builds observers. No production
change, new recovery sample, or rescue of the rejected two-byte candidates.
Sole namespace: `/var/tmp/wh2-aligned-small-basis-screen-r0`.

The three executables are intentionally separate. `trace` interposes ordinary
malloc/free-backed new/delete during non-timing constructors. `controlled`
redirects the constructor allocations to fixed storage for a hot repair test.
`natural` has **no allocator interposition** and times complete public API
lifecycles, including real allocation and deallocation. Traced addresses are
never treated as addresses observed in natural timed samples.

## Controlled repair layer

K3/K5/K8, B1280 full messages, both source policies. Two 16-KiB page-aligned
carriers and raw offsets 0/16/32/48. Handle, evaluator, public source and output
addresses are fixed. The candidate owns M+63 bytes but copies/reads its aligned
interior view; all raw/effective addresses are recorded. Original-pointer frees,
exact copied ranges, guards, packet parity and allocation-free WORK are checked.
Manual ASan poisoning covers live allocation extents and lifetime in neutral
tests. Both carriers are prepared/scanned in fixed order; control gates must
catch preparation asymmetry. This is hot repair work, not allocator cost.

16 pairs per cell: four baseline A/A raw offsets across carriers, four candidate
A/A offsets across carriers, and four C/B offsets within each of two carriers.
All 96 A/A t11 intervals must fit strictly inside `[1/1.02,1.02]`. All 72
nonzero-offset C/B upper bounds must be below 1; the 24 aligned-offset bounds
must be below 1.02 (screening tolerance, not proof of zero regression).

## Natural lifecycle layer

K3 uses its ordinary default constructor; K5/K8 use the explicit small profiles.
Each tests B2-full, B2-tail1, B64, **B256**, B257, B1280-full and B1280-tail1,
under both policies. Explicit certified K16 adds B64/B1280 controls under both
policies. It is only a certified-path sentinel: it does not resolve the still-
open preserved-certified-K3 regression or establish all-K retention.

Three metrics: full constructor/3K-packet encoder/free, low-repair decoder to
its own first success plus recover/free, and distant-repair decoder likewise.
No prebuilt-handle metric is used here; A/A includes fresh construction too.
All arms share public source/output and descriptor/repair-staging addresses within
each cell. Output verification and staging copies are outside WORK; staging
still reads arm-specific fixtures, so identical cache preparation is not claimed.
Raw records carry public workspace pointers, not private allocation addresses.
The `profile` column is the staged descriptor: decoder input or the encoder's
expected descriptor. Encoder constructors write a separate Work-local descriptor;
that output pointer is not recorded or claimed as part of the address checks.

138 cells; comparisons B/B, C/C, WH1/WH1, C/B and C/WH1. All 828 A/A intervals
must fit the reciprocal 2% envelope. The four primary K8/B1280-full encoder
C/B intervals (both policies/observation orders) must have upper bounds below
1; all other C/B upper bounds must be below 1.02. B256 is a retention falsifier
at the first affected width, not assumed to benefit. C/WH1 is reported separately
and does not override a failed control or retention constraint.

## Shared protocol and interpretation

Both timing layers use CPU50 singleton, 12 replicates, two observation orders,
18 observations per panel: first two retained warmups, then eight adjacent
matched pairs. No trimming, pooling, resizing or retries. Each replicate's
mean paired log ratio enters a separate t11 95% interval. Controlled WORK uses
64 batches of K low plus K distant repairs; natural WORK uses 32 lifecycles.
Per DSO load order: 41,472 controlled observations and 298,080 natural
observations, with 192 and 1,380 intervals respectively. Both DSO load orders
must pass. Every control failure overrides gains in either layer.

Before timing, two separate ordinary-malloc traces observe K3/K5/K8,
B64/B1280, both policies and creation orders in both DSO load orders. Each
has 144 allocation records, recording raw ownership and effective copy views.
Partial-tail allocation safety is covered by the earlier neutral qualification,
not these full-message traces or the controlled repair layer.

Build externally with `-DDSO_DIR=<qualified-library-directory>`; use native,
portable-DSO and matching fully sanitized observer/DSO neutral runs. Commit
final code after independent review, then launch once with
`python3 -B bench/Wh2AlignedSmallBasisScreen/Screen.py --run <native-observer>`.
Each timed worker has a 120-second cap; controller timeout is 150 seconds per
worker. Worker failures preserve a spent raw prefix and cannot produce a
successful completion manifest. The controller pins source/binary/build inputs,
not the full transitive compiler/loader environment.

Even an overall PASS is only a bounded screening result. Production adoption
still needs the applicable ordinary/default, preserved-path, package and wider
workload gates. Neither alignment nor these tests change recovery equations.

## Pre-launch neutral checks

Observer builds under `/tmp/wh2-aligned-small-basis-screen.tdqLWHon/` passed
6/6 CTests each for `native`, `scalar` (portable qualified DSOs), and `asan`.
Sanitizer runs enabled leak/fake-stack checks and immediate UB failure. Each
load order checked 384 controlled constructions and 46 natural fixtures ×
three lifecycle metrics × three APIs. All six 144-row non-timing traces also
passed the Python allocation/address parser. Eight synthetic parser/decision
tests passed on Python 3.8 and 3.12. Symbol inspection confirmed the native
natural observer defines no new/delete replacements. The qualified codec
libraries and their previously recorded test logs were not modified.

Initial build-only checks caught and fixed a quoted CMake variable collision,
a pointer comparison requiring `void*`, and the certified constant's versioned
name. These were compile errors before any timing launch, not timing retries.
