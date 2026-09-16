# Small-decoder lifecycle attribution (non-timing)

Issue `wirehair-6juk`. This closes the bounded source/object/allocation inquiry
after the rejected aligned-basis screen. It does **not** identify the cause of
the old timing shifts, rescue that candidate, or qualify a speed improvement.
Production and recovery equations are unchanged.

## Static evidence

The retained native baseline/candidate libraries are the exact qualified pair
under `/tmp/wh2-aligned-small-basis.Zu0KNnhY/native/`:

- Baseline SHA256: `bd3c353847ec2838d95d5b73804664d3a2df8ac76b22ce4f5cf8c27fd58c1bbe`.
- Candidate SHA256: `59a98ab34920cf40c2cfba6aa26b63d8ff929e9ec0823050ab07d21bb2ee05c7`.

Of 580 defined, sized function symbols (including aliases), 559 have identical
addresses, sizes and bytes. `CreateEncoderForProfile` grows 127 bytes, moving
later facade functions and custom recovery sections 128 bytes. The decoder
constructor's 770 normalized instructions remain identical, including six saved
registers and its `0x148` stack reservation. Its two switch tables retain their
addresses and bytes. K3/K8 exception-cleanup helpers exchange addresses
`0x5ef50` and `0x5efe0`; corresponding normalized instruction sequences match.
Normal small-codec destruction is unchanged. Certified decode/recover bodies
have changed call displacements to moved recovery sections, not new work.

The complete `.rodata` (162,824 bytes), initialized data, PLT/GOT and dynamic
relocation **payloads** match. Data/BSS offsets match. Relocation section file/
virtual positions themselves can move (for example `.rela.dyn` moves eight
bytes); equal payloads are not a claim of identical ELF layout. Resolved runtime
global contents/cache state are not established by these static comparisons.

The source confirms three small-decoder allocations: 296-byte facade,
104/136/200-byte evaluator for K3/K5/K8, and `(K+1)*B`-byte slab. Feed and recover
allocate nothing. The encoder's basis allocation is not read by the decoder.
Source policy is not a decoder argument: the two policy-labelled cells for a
shape have the same descriptor/packets but different fixture/heap context.

## Observer and limits

This project builds only an observer, never a codec library. It hash-pins the
retained generated `Natural.cpp` at
`d97b7a4d289764c8729e752dea534a97f1fa6ffdf51f4f5f9965e8a976e1fc6c`, reuses its
fixtures/Work/verification, removes the clock and original main, and records
actual arguments immediately before each WH2 decoder API call. It interposes
C++ new/delete with ordinary malloc/free, without choosing addresses. The fixed
event buffers allocate nothing while recording and fail closed on overflow.

The instrumented worker is newly compiled: its stack addresses are **not**
historical timing addresses. `stack` on allocation events is the address of an
observer-helper local; `arg2` on constructor/recovery calls is the actual API
codec-pointer/count output address. Argument recording itself can change frames.
Caller offsets are relative to the indicated caller owner: zero is the selected
DSO, one the observer. A tail-called final handle delete returns into Work, so
its observer-relative offset must not be added to the recorded codec DSO base.

K3/K5/K8 each use B2-full, B2-tail1, B64, B256, B257, B1280-full and B1280-tail1,
under both source policies: 42 fixtures. Each has low/distant repairs, ABBA/BAAB
orders, and three immediate histories: no additional encoder, baseline encoder,
or candidate encoder. History zero is **not** a pristine heap. Each observation
runs 32 fresh decoder lifecycles using shared descriptor/packet/output buffers.
Printing, fixture preparation and earlier observations still affect heap history.
This is deliberately not a replay of the old natural timing process.

Runtime and parser verify the complete roster, three allocations and reverse
original-pointer frees, no feed/recover allocation, phases, live/private/public
range disjointness, descriptor/output overlap checks, exact packet and recovery
strides. Work checks every status, first-success count, output byte and guard.
Initial RHS offsets are **not** treated as post-feed row identities: the core
rotates RHS/scratch pointers within the unchanged slab.

K16 and WH1 are excluded. K16 has a different allocation contract, while a
new/delete-only observer would miss WH1's direct calloc allocation. No all-K,
new recovery-rate or internal post-feed pointer-trace claim is made.

## Result, 2026-09-16

Final artifacts and observer builds:
`/tmp/wh2-decoder-allocation-trace.CHVlJPAs/`.
For each backend/order, `final-<backend>-<order>.csv.gz` and `.json` contain the
raw trace and strict reconstruction. `.py38.json` matches the Python 3.12 report
exactly. Each run checks 2,016 observations, 64,512 decoder lifecycles, 387,072
allocation/free events and 537,600 API invocations.
`manifest.sha256` pins the five observer/parser source files, retained worker,
two observer executables, six codec DSOs and all 18 final trace/report artifacts;
its SHA256 is `97bba9e5619cb486bd508535309b5bdbf00aa8bebdf7fe466fe7f83370cc9770`.

| Backend / load orders | Lifecycle checks | C/B identical private / helper-stack / public-argument sequences |
| --- | --- | --- |
| Native / normal and reverse | Both pass | 1,008/1,008 for all three comparisons in each order |
| Portable codec / normal and reverse | Both pass | 1,008/1,008 for all three comparisons in each order |
| ASan+UBSan / normal and reverse | Both pass | 0/1,008; allocator quarantine and fake-stack addresses differ |

Native/portable B/B and C/C comparisons likewise match all 504/504 sequences
per order. Sanitizer B/B and C/C addresses also differ in all comparisons; those
runs check safety, not natural address identity. Sanitizers enable leak detection,
fake stacks and immediate UB failure. Eight synthetic tests pass in Python 3.8
and 3.12, including malformed frees, phases, arguments, overlap and stride cases.
Independent reconstruction without importing this analyzer verified both final
native traces, including every event/call, all comparisons, public overlaps and
strides. Both native traces reach first success at K for these fixed fixtures;
this is not a new recovery sample or general recovery-rate claim. Allocation
caller offsets differ by 128 bytes, while the free call offsets are unchanged.

An initial observer assertion incorrectly required the tail-called handle delete
to report a DSO caller; that diagnostic failed before completion and is retained
as `native-normal.csv.gz`. The caller ownership fix is explicit above. Independent
source review also caught missing public-output overlap/stride reconstruction;
those checks and synthetic corruptions were added before the six final traces.
These were non-timing observer corrections, not timing retries. Qualified codec
libraries, generated sources, prior test logs and all spent screens were untouched.

The native/portable observations do not support an arm-specific decoder allocation
or API-stack geometry explanation **in this diagnostic**. They do not establish
the heap/cache state of the rejected screen or prove code-placement causation.
No padding search, unchanged-candidate timing retry or decoder alignment change
is justified by this result.

## Separate optimization lead

Issue `wirehair-fcsn` targets real work visible in the baseline: every small
handle constructs an unused `wirehair_v2::Codec`. Its native constructor at
`0x2d6a0` zeroes 24 bytes of pointer owners plus 216 bytes of `CurrentProfile`.
Free then calls the general destructor's null-state path. These are extra
initialization/destruction operations, not extra allocations.

All nine direct `Impl` expressions are on fresh certified handles or dominated by `SmallK`
guards. A same-layout anonymous union with lifetime conditional on immutable
`SmallK` is therefore a candidate for removing that work without moving fields,
changing allocation sizes, or touching equations. It needs explicit noncopyability,
complete error/cleanup/lifetime review, layout and packet parity, and neutral
native/portable/sanitizer gates before a new prospectively frozen short speed
screen. This is not a cause assigned to the previous regression and not a speed
claim. The rejected alignment modification is not part of this proposed change.

## Reproduction

Use a fresh external 64-bit Linux build. Pass the retained worker source with
`-DNATURAL_SOURCE=/tmp/wh2-aligned-small-basis-screen.tdqLWHon/native/Natural.cpp`.
Configure `-DSANITIZE=ON` only with matching qualified sanitizer DSOs. Invoke
`trace <baseline.so> <candidate.so> normal` or `reverse`; pipe stdout to a fresh
artifact with shell `pipefail` enabled, then run `Analyze.py <trace.csv.gz>`.
Never build in the qualified-library or old screening directories. The observer
contains no timer and does not create a new timing namespace.
