# WH2 borrowed-source encoder design

Status: reviewed design and validation plan.  This document does not add a
public symbol, change the default V2 API, or authorize production promotion.

The retained-source systematic-emission screen recorded by
`wirehair-sxvz.16.1.20.23` established one narrow result: when an encoder
already has stable source storage, direct emission of packet IDs below `K` is
materially faster than reevaluating those packet equations, while repair bytes
and the tested repair timing remain unchanged.  The result does not license an
extra source copy, deferred or lazy solve, a receive-path change, or a claim
against Wirehair1.

## Compatibility boundary

The existing `wirehair_v2_encoder_create()`,
`wirehair_v2_encoder_create_profile_id()`, and
`wirehair_v2_encoder_create_profile()` entry points remain source-independent
after a successful return.  Their signatures, result priority, supported
message/profile aliasing, serialized profile bytes, and implementation policy
remain unchanged.

Source lifetime is a local encoder-storage policy.  It is not part of the
32-byte serialized equation profile, does not change the profile ID, and is
not visible to a decoder.

The future API is additive under the existing `WIREHAIR_2.0` ABI namespace.
It consists of a fixed-layout options record, three constructor variants that
cover every existing V2 creation route, and one detach operation:

```c
#define WIREHAIR_V2_ENCODER_OPTIONS_VERSION 1u

typedef enum WirehairV2EncoderSourcePolicy_t
{
    WirehairV2EncoderSource_Invalid = 0,
    WirehairV2EncoderSource_Independent = 1,
    WirehairV2EncoderSource_BorrowedImmutable = 2,

    WirehairV2EncoderSource_Count,
    WirehairV2EncoderSource_Padding = 0x7fffffff
} WirehairV2EncoderSourcePolicy;

typedef struct WirehairV2EncoderOptions_t
{
    uint32_t struct_bytes;
    uint32_t options_version;
    uint32_t source_policy;
    uint32_t reserved;
} WirehairV2EncoderOptions;

#define WIREHAIR_V2_ENCODER_OPTIONS_INIT \
    { sizeof(WirehairV2EncoderOptions), \
      WIREHAIR_V2_ENCODER_OPTIONS_VERSION, \
      WirehairV2EncoderSource_Independent, 0u }

WIREHAIR_EXPORT WirehairV2Result wirehair_v2_encoder_create_with_options(
    const void* message,
    uint64_t messageBytes,
    uint32_t blockBytes,
    const WirehairV2EncoderOptions* options,
    void* serializedProfileOut,
    uint32_t serializedProfileCapacity,
    uint32_t* serializedProfileBytesOut,
    WirehairV2Codec* codecOut);

WIREHAIR_EXPORT WirehairV2Result
wirehair_v2_encoder_create_profile_id_with_options(
    uint64_t profileId,
    const void* message,
    uint64_t messageBytes,
    uint32_t blockBytes,
    const WirehairV2EncoderOptions* options,
    void* serializedProfileOut,
    uint32_t serializedProfileCapacity,
    uint32_t* serializedProfileBytesOut,
    WirehairV2Codec* codecOut);

WIREHAIR_EXPORT WirehairV2Result
wirehair_v2_encoder_create_profile_with_options(
    const void* message,
    const void* serializedProfile,
    uint32_t serializedProfileBytes,
    const WirehairV2EncoderOptions* options,
    WirehairV2Codec* codecOut);

WIREHAIR_EXPORT WirehairV2Result wirehair_v2_encoder_detach_input(
    WirehairV2Codec codec);
```

The options structure is exactly 16 bytes in version 1.  Its stored policy is
an explicit `uint32_t`; the enum supplies stable named values without relying
on an implementation-defined enum layout.  Static assertions in the library,
C package consumer, and C++ package consumer must pin the structure size,
every member offset, enum size, and numeric values.  The new constructors
require a non-null options pointer that remains readable only for the duration
of the call.  Validation is fail-closed:

| Condition | Result |
| --- | --- |
| `struct_bytes != sizeof(WirehairV2EncoderOptions)` | `WirehairV2_InvalidSize` |
| unknown `options_version` | `WirehairV2_UnsupportedVersion` |
| nonzero `reserved` | `WirehairV2_ReservedNonzero` |
| null options, zero/invalid source policy, or unknown source policy | `WirehairV2_InvalidInput` |

A zero-initialized record therefore never silently selects a policy.  No new
result code is needed.  A null options pointer returns `InvalidInput`; for a
non-null record, simultaneous validation defects are resolved in the table's
top-to-bottom order.

## Source lifetime

`WirehairV2EncoderSource_Independent` is exactly the current public behavior.
The message is needed only during synchronous construction and may be changed
or released after success.  The word “independent” describes the lifetime
guarantee; it does not require the implementation to retain a second full
message copy.

`WirehairV2EncoderSource_BorrowedImmutable` has this contract:

- The exact range `[message, message + messageBytes)` must be readable on
  entry and remain readable, allocated, and byte-for-byte immutable throughout
  the constructor call.  A failed call retains nothing and ends the obligation;
  after success all three requirements continue through every encode until
  successful detach or `wirehair_v2_free()`.
- Construction still completes the existing eager full solve before it
  succeeds.  It adds no full-message allocation or copy.  The existing private
  final-block padding remains permitted.
- Packet IDs below `K` copy only their meaningful source bytes.  A partial
  final packet never reads beyond `messageBytes`.
- Every ID at or above `K`, including `UINT32_MAX`, uses the existing solved-
  intermediate repair evaluator and does not read the borrowed source.
- The library never frees caller storage.  Mutation or early release while
  attached violates the caller contract and is not dynamically diagnosed.
- Operations on one codec are externally serialized.  In particular, encode,
  detach, free, and source mutation must not race.

There is no V2 raw-C reuse operation.  A future C++ `CreateBorrowed()` wrapper
retains the current transactional replacement rule: success frees the old
handle and ends its old borrow; failure preserves both the old handle and its
lifetime obligation.  Move construction transfers the borrow and leaves the
moved-from object with no obligation.  Move assignment first ends the
destination's prior borrow, then transfers the source borrow and likewise
leaves the source empty.  A failed replacement preserves the old handle, so
its source must remain readable, allocated, and immutable throughout the
attempted replacement call.

`wirehair_v2_encoder_detach_input()` is allocation-free and idempotent.  On a
borrowed encoder it clears the borrowed pointer and makes later systematic
packets use equation evaluation.  On an already independent or detached
encoder it succeeds without changing state.  A null handle or decoder handle
returns `WirehairV2_InvalidInput`.  After success the caller may immediately
change or release the former source, and every later successful packet remains
byte-identical.  Equation evaluation of a partial final packet may still
report an existing allocation failure.  The function may not race encode,
free, or source mutation.

An attach-after-construction operation is deliberately excluded: cheaply
proving that a new byte range matches the already solved intermediate state is
not possible.

## Overlap and result priority

Borrowing makes the complete message an immutable retained input.  Rejecting
overlap with only the packet currently being emitted is insufficient because a
write could corrupt another source packet needed later.

The two selecting constructors (`create_with_options` and
`create_profile_id_with_options`) use this exact order:

1. Perform the old constructor's fixed-range output/output and message/output
   alias checks, plus reject every writable output range that overlaps the
   readable 16-byte options object.  A prohibited overlap returns
   `WirehairV2_InvalidInput` without writing any output and outranks every
   options, capacity, profile, dimension, and construction error.
2. Stage the options into a local value and determine its exact validation
   result without publishing anything.  “Structurally selects
   `BorrowedImmutable`” below means that the non-null record has the exact v1
   size and version, zero reserved field, and the
   `BorrowedImmutable` policy value.  A record with any options error does not
   structurally select a policy and cannot publish descriptor bytes or retain
   source storage.
3. If and only if the staged options structurally select
   `BorrowedImmutable`, reject an unrepresentable or wrapping complete message
   range and reject exact or partial overlap between that complete range and
   `serializedProfileOut`.  These are also no-write `InvalidInput` failures and
   outrank short capacity and all later ordinary errors.
4. After alias preflight, preserve the old bookkeeping order: publish the
   required descriptor byte count when its pointer is valid, clear a non-null
   `codecOut`, and validate the required output pointers.
5. Return the staged options error, if any.  In particular, an options error
   outranks `WirehairV2_BufferTooSmall`.
6. Check descriptor capacity, profile selection, dimensions, and construction
   in the existing order.  A malformed options record never allocates a codec
   or publishes descriptor bytes.

The borrowed-only message/profile-output rule must not be added to the old
constructors or to a new constructor selecting `Independent`: those paths
intentionally stage the descriptor and preserve the existing exception that
allows descriptor output to overlap a message that is not retained after
success.

The serialized-profile constructor has a different staged order because the
message length is available only after successfully parsing the descriptor:

1. A null `codecOut` returns `WirehairV2_InvalidInput` immediately, without
   reading or writing any other argument, and therefore outranks descriptor
   and options errors exactly as in the old constructor.
2. Reject fixed-range `codecOut` overlap with the descriptor bytes or readable
   options object.  This is a no-write `WirehairV2_InvalidInput` failure and
   outranks every other result.
3. Stage the options and parse the complete descriptor into local values,
   without publishing a facade.
4. If descriptor parsing succeeds, reject `codecOut` overlap with the
   descriptor-defined complete message range.  If the staged policy is
   structurally `BorrowedImmutable`, also reject a wrapping or otherwise
   unrepresentable retained range.  Either failure is no-write
   `WirehairV2_InvalidInput` and outranks options and construction errors.
5. Only after those checks, clear `codecOut`.  Return the exact
   existing descriptor parse error before any staged options error; otherwise
   return the options error before attempting construction.
6. With both inputs valid, retain the existing message/profile/dimension and
   construction result order.

Thus the serialized-profile cross-product is pinned as follows:

| Condition | Result priority and writes |
| --- | --- |
| null `codecOut` | `InvalidInput`, no writes or other input reads; wins over malformed descriptor and options |
| fixed `codecOut`/descriptor or `codecOut`/options overlap | `InvalidInput`, no writes; wins over malformed descriptor and options |
| malformed or short descriptor, no fixed overlap | existing descriptor error; `codecOut` is cleared when non-null; wins over options errors |
| valid descriptor and `codecOut` overlaps its message | `InvalidInput`, no writes; wins over options and construction errors |
| valid descriptor/ranges and malformed options | exact options error after clearing `codecOut` |
| valid descriptor, ranges, and options | existing construction result order |

Message, descriptor, and options are read-only inputs and may overlap each
other after staging.  The tests must cover every row above with valid and
malformed options, a short and malformed descriptor, and exact and partial
overlaps.  The three new `Independent` constructors must also reproduce the
old descriptor/message alias exception byte-for-byte.

For a valid attached borrowed encoder, `wirehair_v2_encode()` first computes
the exact required packet length without writing an output.  It then performs
these checks in order:

1. If `dataBytesOut` is non-null, reject its overlap with any byte of the
   complete borrowed-message range.  This check precedes the existing invalid-
   output path that may zero the counter.  In particular, a null
   `blockDataOut` plus a counter inside the source must reject without changing
   that counter.
2. Keep the existing packet-output/counter overlap and address-range checks,
   using the exact required packet range even when the supplied capacity is
   short.
3. Reject that same exact required packet range overlapping any byte of the
   complete borrowed message, for systematic and repair IDs alike.
4. Only after all no-write alias checks may the function publish the required
   byte count, report `BufferTooSmall`, or write packet storage.

Such a call returns `WirehairV2_InvalidInput` without changing the source,
packet output, or counter.  Once detach succeeds, only the existing ordinary
overlap rules apply.

## Implementation plan

The preferred implementation seam is the private `PublicCodec` facade in
`codec/WirehairV2Profile.cpp`, not the serialized profile and not the internal
`CacheSystematicSource` option.  The latter means “own a second full copy” and
has different memory and timing semantics.

The implementation steps are:

1. Run the exact fixed-alias preflight and stage the public options into a
   local value in the result-priority order above.
2. Run the unchanged `Codec::InitializePrecodeEncoder()` eager solve.
3. Only after complete success, bind the original message pointer and the
   borrowed policy in the unpublished `PublicCodec`.
4. In the public encode facade, use the retained source only for IDs below
   `K`; otherwise call the existing `Impl.Encode()` repair path.
5. Detach only clears facade state.  Destruction never dereferences or frees
   the borrowed pointer.

An unpublished or failed facade stays in the invalid policy state.  Every
successful old constructor and every new constructor selecting `Independent`
publishes the independent state explicitly; only a successful borrowed
constructor publishes a non-null source pointer.  This makes detach on an old
encoder a defined successful no-op rather than relying on a zero/default enum.

This keeps `MessagePrecodeEncoder`, its owning cache, the equation profile,
precode construction, and the evaluator source unchanged.  Because the facade
does add dispatch around an otherwise unchanged core call, promotion requires
a separately sealed pre-change/default timing closure; policy-specific A/B
neutrality alone is not enough.

The C++ wrapper adds three explicit `CreateBorrowed()` overloads and
`DetachInput()`.  It must retain no hidden source owner: the C++ caller has the
same lifetime obligation as the C caller.  A later Python wrapper would have
to pin the exporting object until detach/close and release it only after a
successful native detach.

## Required implementation gates

Before any production promotion:

- Prove the three old constructors retain their existing source-independent
  behavior and descriptor/message alias contract.
- Compare all three new independent constructors with their old counterparts:
  descriptor, configuration, intermediate state, every systematic packet, and
  representative/high-ID repairs must match.
- Test borrowed exact and partial-tail lengths `(K-1)B+1`, `(K-1)B+B-1`, and
  `KB`, including guard-page or sanitizer coverage for over-read.
- Test borrowed detach twice, then poison or unmap the source and reproduce all
  tested systematic and repair packets through equation evaluation.  Each of
  the three old encoder routes and each new `Independent` route must also
  return `Success` for repeated no-op detach while retaining packet identity;
  null and decoder handles must retain their pinned `InvalidInput` result.
- Exhaust exact/partial overlap boundaries for options, all constructor
  outputs, packet output, and the counter, including short-buffer calls.
- Exercise no-dereference address-range failures near `UINTPTR_MAX` for the
  message, packet, and counter on supported flat-pointer test builds, and a
  `messageBytes > SIZE_MAX` case on narrower-size platforms; each must return
  its pinned result without wrapping an address or modifying an output.
- Test null/zero/short/new-version/reserved/unknown policy records, allocation
  failures, null/decoder detach, and C++ failed replacement transactionality.
- For each new constructor route, include at least one options record with
  simultaneous size/version/reserved/policy defects and cross it with another
  ordinary failure so the pinned top-to-bottom result priority is observed,
  not merely inferred from single-defect tests.
- Test borrowed C++ move construction and move assignment: moved-from wrappers
  are invalid and retain no obligation, moved-to wrappers retain direct packet
  identity and the original source obligation, assignment ends the
  destination's prior borrow, and wrapper `DetachInput()` is idempotent and
  releases the transferred obligation only after success.
- Update the public header, C++ wrapper, `V2_WIRE_PROFILE.md`, ELF export map,
  C and C++ package consumers, fuzz operations, shared/static packaging, and
  sanitizer/E2E gates.
- Prove the default and repair core objects are unchanged or provide a sealed
  pre-change/default A/B closure with no regression.

## Broader timing license

The broader timing campaign remains encoder-only.  It may compare the current
equation and direct modes in one binary, but that binary must expose, schedule,
and serialize only those two WH2 arms.  Dormant Wirehair1 implementation text
from the shared `Wh2NativeCodec.cpp` translation unit may remain linked; it is
not a selectable axis and supplies no comparison evidence.  Any actual
Wirehair1 or literal pre-change comparison uses separately sealed role
binaries, because even equation-neutral hooks have produced percent-level
layout bias historically.

The screen must retain:

- full constructor plus IDs `[0,K)` as the primary measured scope;
- true equation/equation and direct/direct A/A panels;
- exact source, configuration, descriptor, systematic-byte, and repair-byte
  receipts.  This NativeArm screen uses exact `K*B` sources; partial-final-
  block evidence belongs to the later public-API implementation gates, not to
  this timing falsifier.  In the same-binary falsifier, a binary-bound
  `construction_equivalent:true` may attest the field-by-field comparison of
  system rows, runtime primes, intermediate bytes, and normalized non-time
  statistics; the later cross-binary campaign must serialize a stable state
  fingerprint for those values rather than relying on that Boolean alone;
- exactly `K` qualified direct packets in the primary scope;
- a branchless/pre-change repair control in addition to direct/equation policy
  neutrality;
- canonical ABBA/BAAB records, a hard external deadline, and a single newly
  created immutable output directory.  Before the one-shot, an independent
  fresh-build audit must bind compiler identity, normalized compile/link
  commands, options, object roster, source manifest, and the executable hash;
  the runner revalidates that source and executable seal before and after the
  child.  The later cross-role campaign embeds the complete per-role build
  receipts in its own sealed provenance;
- preregistered cell and pooled confidence gates.

The development timing domain is the existing 24-cell Cartesian product:

```text
K = {8, 32, 100, 128, 512, 1000, 2048, 5000, 8192,
     20000, 32768, 64000}
B = {64, 1280}
```

It is deliberately smaller than the later production domain.  Work for one
cell is distributed over 12 independent replicates as:

```text
N(K) = max(24, ceil(65536 / K))
n[r] = floor(N / 12) + (r < N mod 12)
```

Thus every replicate executes at least two invocations per slot and contains
both the primary and opposite four-slot counterbalanced blocks.

The first gate is a new, never-before-timed 10-cell falsifier:

```text
K = {8, 128, 512, 5000, 64000}
B = {64, 1280}
```

It uses only the current same-binary equation (`E`) and direct (`C`) WH2 arms,
never Wirehair1, with a 115-second internal and 120-second external deadline.
For every cell it measures baseline A/A, candidate A/A, and `C/E` in these
three scopes:

1. fresh eager initialization plus systematic IDs `[0,K)`;
2. fresh eager initialization plus repair IDs `[K,2K)`, retaining both the
   nested initialization-plus-first-repair timestamp and the full total;
3. a prebuilt encoder evaluating repair IDs `[K,2K)`.

Every A/A 95% interval must lie strictly inside
`[-log1p(0.02), +log1p(0.02)]`.  For the systematic `C/E` contrast, no cell
point estimate may exceed `log1p(0.02)`, and the equal-cell aggregate and each
width and size-band aggregate must have upper 95% bounds below `log(0.99)`.
The preregistered bands are `small={8,128}`, `medium={512,5000}`, and
`large={64000}`, each pooled equally across its two payload widths.  For every
repair metric, no cell point may exceed `log1p(0.02)` and all cell, width,
size-band, and equal-cell aggregate upper 95% bounds must be below
`log1p(0.02)`.  Statistics use paired log contrasts, `math.fsum`, unbiased
variance, and
`t11 = 2.200985160082949`.  A valid pass licenses only the full role-separated
campaign.

The full campaign uses two independently built and sealed workers:

- current source exposes `C` (direct) and `E` (equation evaluation);
- exact pre-change source commit
  `87a3f1cb1d6562a6652c1ec0b99f6131eaeade09` exposes `P` (pre-change WH2)
  and `L` (pre-change Wirehair1).

The pre-change worker must never compile the current `Wh2NativeCodec.cpp`.
Each role receipt binds its Git commit, ordered source/blob manifest, compiler,
normalized compile and link command, CMake options, object/archive roster,
executable hash and metadata, and dynamic dependency closure.  A fresh
reproducibility build must match before the one-shot campaign.  Workers run
sequentially on one pinned CPU while the controller and health sampler stay on
other physical cores.

The full campaign covers all 24 development cells and the same three timing
scopes.  It interprets `C/E` as the policy effect, `E/P` as the contamination
closure, `C/P` as the historical-WH2 comparison, and `C/L` only as isolated
encoder evidence against historical Wirehair1.  `P`, `E`, and `C` must agree
on the exact configuration and attempts, source, system rows, runtime primes,
intermediate bytes, normalized non-time statistics, all systematic bytes,
repairs `[K,2K)`, and high-ID repair controls.  `C` records exactly `K` direct
systematic hits and zero repair hits.

`E/P` must be equivalent within the strict two-sided 2% interval for every
cell and scope.  `C/E`, `C/P`, and the repair gates retain the falsifier's
noninferiority and aggregate-improvement directions; `C/L` is descriptive
unless a separate stricter WH1 gate is preregistered.  The campaign has a
1770-second internal and 1800-second external deadline, canonical exact-roster
JSONL, empty stderr, full CPU/thermal/EDAC coverage, and immutable raw and
summary artifacts.  Missing data, failed outcomes, deadline/health failures,
or any A/A failure invalidate the whole campaign; no cell may be dropped.

No result from that campaign establishes a receive-path improvement, changes
the default API, proves Wirehair1 superiority, or authorizes the public symbols
above without the implementation gates in this document.
