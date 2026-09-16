# Small-profile private-basis placement diagnosis

Issue `wirehair-m9f5`. Pure GF(256); no production changes. This is a new
non-timing allocation observation, not a replay of the rejected payload screens
and not a reconstruction of their unrecorded heap addresses.

`Trace.cpp` loads the exact baseline/outlined DSOs without rebuilding them.
Exported C++ new/delete replacements preserve ordinary malloc/free semantics
and record the three allocations made by each full-message small constructor.
It requires the handle, private byte-array basis, and evaluator allocation order;
it checks the copied bytes, descriptor and repair parity while handles are live.
It never casts an opaque handle or assumes its layout. No timer is present.
Output is deferred until all 144 allocation records are complete.

Build externally with CMake/Ninja, then run `trace <baseline.so> <candidate.so>
normal` and `reverse` for the two DSO load orders. Both constructor orders,
K3/K5/K8, B64/B1280, and both ownership policies are observed. Sanitizer builds
must load fully instrumented DSOs and use leak/fake-stack/UB checking.

The production constructor copies into its own `SmallBasis` under both source
policies; repair packets read that basis. Thus matching the caller's source
pointer does not match actual repair-source addresses. The native GFNI wide
kernel performs 64-byte unaligned vector loads with 64-byte stride. An unaligned
basis at B64/B1280 makes those vector loads cross cache-line boundaries. These
facts motivate checking placement; they do not establish the cause of a prior
timing regression or justify a production alignment change by themselves.

## Prospective baseline-only intervention

`Controlled.cpp` uses the **unchanged baseline DSO only**. During each
constructor, its exact three allocations are directed into fixed, adequately
aligned storage: one handle buffer, one evaluator buffer, and a basis view
at offset 0/16/32/48 in one of two page-aligned 16-KiB carriers. All three
deletions are intercepted and checked before the next observation. No opaque
handle is inspected or mutated. Construction/destruction is outside WORK.
Manual ASan poisoning enforces exact live-object bounds and lifetime in
sanitizer neutral tests; byte guards and source/packet checks run in all modes.

Handle, evaluator, public source and output addresses are fixed throughout
the process. Every raw record carries them plus the active basis address.
Both carriers are reset and checked in fixed order before WORK. Those reads
and constructor copies warm memory; this is a **hot repair-kernel** experiment,
not a lifecycle, allocator-performance or cold-cache benchmark. Cross-carrier
controls also test preparation-order/cache effects and must all pass.

Frozen roster: K5/K8 × B64/B1280 × both source policies; 12 replicates;
64 batches of K low plus K distant repairs per observation; CPU50 singleton;
18 observations per panel, first two retained warmups, eight matched pairs;
both observation orders. Four control pairs compare equal offsets across the
two carriers. Six treatment pairs compare each nonzero offset against offset
zero within each carrier. Handle/evaluator placement is identical in every
pair. Total 34,560 observations and 160 separate t11 95% intervals.

Every one of the 64 cross-carrier A/A intervals must lie strictly inside
`[1/1.02,1.02]`. The primary effect requires all 24 K8/B1280 treatment lower
bounds (both policies/carriers/orders, offsets 16/32/48) to exceed 1. Any
control failure overrides apparent effects. Other effects are descriptive;
no pooling, trimming, resizing, retry or candidate rescue. Even a detected
effect would not prove the prior slowdown's cause, because its actual private
addresses were not retained, or establish this candidate's pure overhead.

Sole namespace `/var/tmp/wh2-small-basis-alignment-r0`. After neutral tests and
independent review, `Screen.py --run <observer-build> --baseline <baseline.so>
--candidate <outlined.so>` records the two non-timing allocation traces, then
the baseline-only intervention. The outlined DSO is used only by the trace,
never by the timed worker. Both DSO hashes are fixed to the previous screen's
recorded native libraries. Worker wall cap 120 seconds, controller cap 150.
A created namespace is spent even if a worker fails; raw prefixes are kept.
The receipt pins sources and binaries, not a full transitive toolchain closure.

## Integration checkpoint

The diagnostic is committed as tooling, not as a production optimization or
a completed timing result. The one-shot timing namespace has not been claimed.
Neutral checks passed under Python 3.8 and 3.12 (seven controller tests each).
Native, portable-arithmetic, and ASan/UBSan baseline checks each passed all
320 controlled constructions, including packet parity, placement, guards and
deallocation. Native and ASan/UBSan traces passed in both DSO load orders,
with all 144 records per trace accepted by the controller. Sanitizer checks
enabled leak detection, fake-stack checks and immediate UB failure.

Issue `wirehair-m9f5` remains open for the separately launched timing diagnostic
and its independent result audit. Neither previous payload candidate is
promoted by these neutral checks.
