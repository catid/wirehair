# WH2 borrowed-source encoder design

Status: implemented and functionally validated by
`wirehair-sxvz.16.1.20.23.1.3`.  The API remains explicit opt-in and does not
change the default V2 behavior or authorize production performance promotion.

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

The API is additive under the existing `WIREHAIR_2.0` ABI namespace.
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

There is no V2 raw-C reuse operation.  The C++ `CreateBorrowed()` wrapper
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

## Public-facade timing closure

Performance promotion is a separate, encoder-only evidence task tracked by
`wirehair-sxvz.16.1.20.23.1.3.1`.  Earlier retained-source screens used an
internal `NativeArm` and cannot establish the timing of the exported facade,
default-policy neutrality, an exact pre-feature comparison, or superiority to
Wirehair1.  Their observations are not combined with this task's statistics.

The required first gate is exactly one sealed public-C-API falsifier using two
independently built role binaries.  Both binaries compile the identical worker
translation unit from the later harness commit, but each includes the public
header and links the library built from its own exact implementation commit:

- current commit `d6ddb35e046956174202584fb5a26c9a79679ea8` exposes `C`,
  `E`, and `D`;
- exact parent commit `101584b7e5c30326b1791429221c331c82a00807` exposes
  `P` and `L`.

The arms use only exported constructors and encode entry points:

- `C` calls `wirehair_v2_encoder_create_with_options()` with
  `WirehairV2EncoderSource_BorrowedImmutable`;
- `E` calls the same constructor with
  `WirehairV2EncoderSource_Independent`;
- `D` calls the existing `wirehair_v2_encoder_create()` at the current
  commit;
- `P` calls `wirehair_v2_encoder_create()` at the exact parent commit;
- `L` calls the parent commit's public Wirehair1 constructor and encoder.

The role libraries are never built from the harness worktree.  Each build
receipt binds the implementation and harness Git trees, ordered source/blob
manifest, compiler and linker identities and commands, CMake configuration,
object and archive closure, executable and library hashes, and resolved
dynamic-dependency closure.  Before the sole run, strict builds, selftests,
sanitizers, and repeated local and independent adversarial reviews must all
complete without a P0 or P1 finding.

The falsifier covers this exact 10-cell Cartesian product:

```text
K = {8, 128, 512, 5000, 64000}
B = {64, 1280}
```

For each cell, work is distributed across 12 matched replicates as:

```text
N(K) = max(24, ceil(65536 / K))
n[r] = floor(N / 12) + (r < N mod 12)
```

Every replicate measures these four scopes:

1. a prebuilt encoder emitting systematic IDs `[0,K)`;
2. fresh eager construction followed by systematic IDs `[0,K)`;
3. fresh eager construction followed by repair IDs `[K,2K)`, including a
   nested construction-plus-first-repair clock;
4. a prebuilt encoder emitting repair IDs `[K,2K)`.

Each scope contains the five true A/A controls `D/D`, `E/E`, `C/C`, `P/P`,
and `L/L`, plus `E/D`, `C/E`, `D/P`, and `C/L`.  The controller schedules the
cross-role invocations as matched ABBA/BAAB blocks on one pinned worker CPU.
The complete roster is therefore exactly
`12 replicates * 10 cells * 4 scopes * 9 comparisons = 4,320` panels.

The public evidence must bind the exact source and work allocation, canonical
per-arm descriptors, serialized WH2 configuration, public-state fingerprint,
systematic bytes, repairs `[K,2K)`, first repair, representative high IDs, and
public decode roundtrip.  `C`, `E`, `D`, and `P` must have identical WH2
source, configuration, public state, systematic, repair, and high-ID receipts.
`L` is a different equation system: it must reproduce the source systematics,
retain exact stable role-local repair and high-ID receipts, and pass its public
decode roundtrip, but it is never required to equal WH2 repairs or
configuration.

The timing worker must not link an internal codec or expose a timing-only path
counter.  Instead, it authenticates the attached `BorrowedImmutable` policy
and records exactly `K` systematic invocations eligible for the borrowed route
and zero eligible repair invocations.  The already-passed production tests and
source review, not this timing receipt, remain authoritative for direct-route
execution.

Statistics use matched log contrasts, `math.fsum`, unbiased variance, and
`t11 = 2.200985160082949`.  The frozen gates are:

- every A/A interval, and every `E/D` and `D/P` interval, lies strictly inside
  the two-sided `[-log1p(0.02), +log1p(0.02)]` range for every cell and scope;
- in both systematic scopes, neither `C/E` nor `C/L` has a cell point estimate
  above `1.02`, and each equal-cell, width, and preregistered size-band upper
  95% bound is strictly below `0.99`;
- the bands are `small={8,128}`, `medium={512,5000}`, and `large={64000}`;
- in both repair scopes, `C/E` has no cell point estimate above `1.02`, and
  every cell and group upper 95% bound is at most `1.02`.  Wirehair1 repair
  timing is descriptive because `L` uses different equations.

The controller's exact internal deadline is 840 seconds inside the 900-second
external scientific guardian.  The controller stops worker timing at T+820
and reserves its final 20 seconds for health finalization; the systemd
activation/stop backstop is separate from the scientific deadline.  Evidence
must be canonical, complete, immutable,
and independently replayable, including raw records, summary, role build
receipts, source and executable closure, affinity, containment, thermal,
DIMM/EDAC health, durable attempt accounting, and `COMPLETE` published last.

A valid pass licenses only creation and execution of a separately tracked and
preregistered role campaign over the full 24-cell development domain:

```text
K = {8, 32, 100, 128, 512, 1000, 2048, 5000, 8192,
     20000, 32768, 64000}
B = {64, 1280}
```

Any A/A gate failure is invalid evidence; a candidate/equivalence gate failure
with otherwise complete infrastructure is a valid scientific reject.  The
screen does not itself establish receive or decoder performance, all-K recovery,
default-policy promotion, production promotion, or final Wirehair1
superiority.  A valid reject is final for this claim.  Reject or invalid
evidence licenses neither a rerun nor post-hoc cell or scope narrowing, and in
all outcomes the functional opt-in API remains unchanged.
