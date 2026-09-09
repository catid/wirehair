# Serialized V2 wire profile

The V2 packet API carries its complete equation identity in a canonical
32-byte descriptor. An encoder selects and serializes this record; a decoder
needs only the record and packet `(id, payload)` pairs. No in-process
`wirehair_v2::SeedProfile` value is exchanged.

The descriptor identifies equations. It does not authenticate the descriptor,
packets, or recovered message. Applications must authenticate framing and
verify a trusted digest or MAC before accepting recovered bytes.

## Canonical encoding version 1

Every multibyte integer is unsigned little-endian. The record is exactly 32
bytes; shorter and longer inputs are rejected. This exact-length rule prevents
a version-1 parser from silently ignoring data that a newer sender considers
meaningful.

| Offset | Bytes | Field | Version-1 rule |
| ---: | ---: | --- | --- |
| 0 | 4 | magic | ASCII `WHV2` |
| 4 | 2 | encoding version | `WIREHAIR_V2_PROFILE_ENCODING_VERSION` (`1`) |
| 6 | 2 | encoded bytes | `32` |
| 8 | 8 | equation profile ID | a supported `WIREHAIR_V2_PROFILE_*` value |
| 16 | 8 | message bytes | exact nonzero original length |
| 24 | 4 | block bytes | `1..2^31-1` |
| 28 | 1 | seed attempt | certified: `0..255`; small K3: `0` |
| 29 | 3 | reserved | all zero |

The derived block count is `ceil(message_bytes / block_bytes)` and must be in
`2..64000`. It is intentionally not serialized a second time.

The golden record for a 117-byte message, 16-byte blocks, attempt zero, and the
current profile is:

```text
57 48 56 32 01 00 20 00 c9 f9 f4 47 bb 5b 29 4b
75 00 00 00 00 00 00 00 10 00 00 00 00 00 00 00
```

`wirehair_v2_profile_serialize()` and
`wirehair_v2_profile_deserialize()` are the only supported conversion between
this byte encoding and the host-native `WirehairV2Profile` ABI structure.
`wirehair_v2_profile_validate()` performs the same byte validation when the
host representation is not needed. Copying that C structure directly to the
wire is not portable.

## Supported certified equation profile

`WIREHAIR_V2_PROFILE_CERTIFIED_2026_07` has numeric ID
`4b295bbb47f4f9c9`. The ID is the first 64 bits of SHA-256 over this exact
canonical UTF-8/ASCII input, with no newline:

```text
wirehair:v2:precode-v2:packet-v4:certified-2026-07
```

The full digest is
`4b295bbb47f4f9c91ebf12ba77afc33cf1c7c36131d154aea03c320f0f13dcf4`.
As with the legacy profile IDs, this hash is a stable compatibility name, not a
security primitive.

The profile freezes all equation-affecting rules used by the public V2 API:

- the current legacy dense-count, dense-seed, and peel-seed selection tables;
- V2 precode contract 2 and packet-row contract 4;
- the certified staircase, Shuffle-2 dense, and 12-row Cauchy-heavy geometry;
- source-hit selection, packet degree sampling, packet-ID mapping, and exactly
  three distinct precode mix columns;
- the precode salt `0x763263707265636f` and recovery-row salt
  `0x76327265636f7665`;
- the full-span dense corner (the experimental identity corner is disabled);
- deterministic attempt stepping and the first selected attempt published by
  the encoder; and
- GF(256) arithmetic and all integer PRNG/shuffle rules consumed by those
  equations.

Builds compiled with the private `WH_SEED_KNOBS` experiment switch can change
base seeds or packet degree weights at runtime. They deliberately return
`WirehairV2_UnsupportedPlatform` for the named profile instead of publishing
or consuming incompatible equations under its ID.

The descriptor serializes only canonical inputs: profile ID, message length,
block size, and the selected attempt. The block count, base seeds, peel policy,
precode dimensions, expanded precode/packet seeds, salts, mix count, fixup
diagnostics, and tuning statistics are derived under the named profile and are
not duplicated on the wire. Publishing the current ID after changing a frozen
rule is a compatibility bug.

## Small K3 equation profile

`WIREHAIR_V2_PROFILE_SMALL_K3_2026_09` has ID `67c1043ecaa9e184`, the first
64 bits of SHA-256 over the exact name (no newline):

```text
wirehair:v2:small-k3:thue-morse-8-14-7:gf256-14d:2026-09
```

The full digest is
`67c1043ecaa9e1847ffb885eaf845b3af93ce52df8493bb1c6cf5428e863bdb9`.
It freezes the same GF(256), polynomial `0x14d`, Thue-Morse companion pair
`(8,14,7)` / `(9,14,7)`, packet-ID map, systematic identity rows and zero
tail padding as the qualified [WHK3 equations](SMALL_WIRE_PROFILES.md).
The two facades share one immutable 13,056-byte lookup; WHV2 framing and
source-ownership rules remain distinct from WHK3.

This profile requires exactly three source blocks, block bytes at most
67,108,864 and seed attempt zero. Other shapes return `InvalidDimensions`;
a nonzero attempt on a valid shape returns `BadSeed`. There is no seed search.
The prepared encoder basis is the original three source blocks, with private
padding for a partial final block. It is owned for both storage policies:
borrowing adds no second full-message cache, no extra construction allocation,
and no extra copy. Repair packets never read borrowed input; detach remains
allocation-free. Dependent contradictory packets report `Error` without
changing the retained decoder basis; recovery does not authenticate data.

On this development branch, ordinary constructors select this profile for its
supported K3 shapes. Other shapes, including K3 above the small-profile block
bound, continue to select the certified profile. Correctness and compatibility
checks, the scoped ordinary-path speed gate and retained paired recovery
qualification below pass. Old-path/opt-in performance regression checks remain
open.
The published WHK3 opt-in timings do **not** qualify this owned-basis path.
This is not an all-K performance or construction-seed claim.

### Ordinary K3 full-lifecycle speed

At source `d77a9db`, `wirehair.wh2.k3-ordinary-cost-r0` measured the actual
ordinary source-independent constructor and ordinary borrowed constructor
against actual WH1 and explicitly selected certified WH2. All 72 same-code
timing controls and all 72 candidate comparisons passed. Every candidate
comparison's upper 95% time-ratio bound was strictly below one; there was no 5% minimum
gain, treatment regression tolerance, sample filtering or order pooling.

Observed elapsed-time reductions versus WH1, spanning both measurement orders:

| Source policy | Block bytes | Encoder | Low-ID decoder | Distant-ID decoder |
|---|---:|---:|---:|---:|
| Independent | 2 | 62.55–62.59% | 94.35% | 89.85–89.92% |
| Independent | 64 | 70.11–70.16% | 94.07–94.24% | 90.53–90.57% |
| Independent | 1280 | 71.38–72.14% | 93.29–93.40% | 91.90–92.09% |
| Borrowed | 2 | 61.58–61.68% | 94.31–94.33% | 89.88–89.90% |
| Borrowed | 64 | 69.11–69.19% | 94.14–94.22% | 90.58% |
| Borrowed | 1280 | 70.71–71.84% | 93.29–93.36% | 91.90–92.03% |

The run used full three-block messages in the static library on one
GFNI-capable host. Encoder time includes create, descriptor output, 18 packets
and free. Decoder time includes create, feed through its own first success,
recover and free; both frozen streams needed three packets for every arm.
The two source policies share the decoder: its repeated observations are not
extra recovery samples. This screen does not establish partial-tail,
additional-width, shared-call, cold-start, non-GFNI or all-K performance, nor
does it isolate regressions in the preserved old profile or opt-in facade.

All 31,104 callbacks, 3,981,312 fresh codec lifecycles and 42,467,328 API calls
were retained, with every output, descriptor, endpoint and buffer guard
checked. Native, portable-arithmetic and sanitizer neutral fixtures agree;
portable and sanitizer builds were not scientific timing arms. Python 3.8
replayed the complete result and all 582 receipt pins exactly. A separately
written standard-library chronology, ledger and statistical audit reproduced
all 144 decisions before source HEAD or pinned documentation changed.

The immutable outcome is `/var/tmp/wh2-k3-ordinary-cost-r0`, raw SHA-256
`03ad310ea6d44e7a98ece8ba27c96ab8284741df38c43b14715be4f7c66503d2`.
This namespace is spent. No recovery-rate sample or full promotion claim is
made by this timing result.

### Ordinary K3 retained recovery comparison

At source `734150a`, `wirehair.wh2.k3-ordinary-recovery-r0` replayed the
retained cohort through actual ordinary WH2, separately for independent and
borrowed creation, alongside actual WH1 and explicitly selected certified
WH2. Failures with zero extra packets on the 6,144 retained loss traces:

| Codec | Failures / traces | Observed failure rate |
|---|---:|---:|
| Ordinary WH2 K3, each source policy | 0 / 6144 | 0% |
| WH1 | 14 / 6144 | 0.23% |
| Certified WH2 | 269 / 6144 | 4.38% |

Both ordinary policies fixed all 14 WH1 and 269 certified-WH2 failures,
introducing none against either control in this cohort. Each of the twelve
width/loss-schedule cells retained 512 traces and had no K3 failures. WH1
recovered all failures with one extra packet; certified WH2 still had five
failures after one extra packet, one after two, and none after three.

The separate 72 hard cases also had no K3 failures at zero overhead, versus
one WH1 and five certified-WH2 failures. K3 and WH1 recovered all 53 historical
prefixes at their original widths; certified WH2 remained unresolved on 43.
Those prefixes were never extended to rescue an unresolved decode. Hard and
historical cases are not included in the 6,144-trace denominator.

Every encoder was freed before receiver creation. Both ordinary policies
produced identical descriptors, packet hashes, feed statuses and first-success
counts. Each K3 packet matched an independent polynomial-GF(256) payload
oracle; first success matched independently computed rank, and each successful
decode recovered the original message twice with guarded buffers. Native,
portable-arithmetic and ASan/UBSan builds agreed on every one of the 6,269
cases and all four API routes, with 50,152 fresh codec handles per backend.

Python 3.8 independently replayed all 647 receipt pins and the complete raw
results. A separately written standard-library audit rebuilt the GF(256)
payload/rank oracle and reproduced the full chronology, all group/cell counts
and every paired fix/introduction list before HEAD or pinned files changed.
The immutable bundle is `/var/tmp/wh2-k3-ordinary-recovery-r0`; native raw
SHA-256 is `1beaa881fb2a8ce1f13b58fd463271dc0f31b6dc7179cb1bddbce0bb6a370f86`.

This namespace is spent. The cohort is retained, not a fresh independent
holdout; source-policy and backend replays do not multiply its sample size.
The result is not a universal recovery guarantee, all-K construction-seed
validation or an old-path performance regression check.

### Preserved-path regression and small-state isolation

The pre/post shared-library screen at source `26304b3`,
`wirehair.wh2.admission-regression-cost-r0`, finished **REGRESSION**.
All 320 same-code controls passed across the two DSO load orders, but
57 and 56 treatment cells respectively showed resolved slowdowns. Certified
K2/K3/K4/K6 encoder and decoder lifecycles took 4.8-9.1% more time after the
ordinary K3 admission. The screen also retained order-sensitive K128 results
and changes in WH1 and opt-in paths whose producing objects were unchanged;
it does not establish allocation size as the sole cause. This is a native,
co-resident shared-library result, not a static-link or all-K claim.

The exact Python 3.8 replay and a separately written raw/statistical audit
reproduced the outcome before source advancement: 643 pins, 103,680 records
and 480 statistical decisions. The spent immutable bundle is
`/var/tmp/wh2-admission-regression-cost-r0`, with analysis SHA-256
`4bd1c590a36e677655524e04215959ef26fe3342ca01d4453b089d965ee15ada`.
These regressions remain an unmet ordinary-admission qualification requirement.

Subsequently, a `ctest -N` listing overwrote the pre-isolation build's
`Testing/Temporary/LastTest.log`, a non-producing diagnostic pinned by an older
receipt. The library objects and sealed measurement bundles are unchanged,
but that historical receipt's full artifact closure is no longer intact.
The follow-on screen records the lost log explicitly and proves all producing
inputs separately; the original strict verifier still rejects the mismatch.

The follow-on implementation moves K3-only owning pointers to a private
derived handle and keeps its operation bodies out of the common dispatch
functions. The immutable type tag uses existing header padding on the tested
ABI; no virtual dispatch, new allocation, equation or source-policy change
is introduced. Certified handle allocation returns from 296 to 272 bytes on
the native host. Allocation/OOM tests cover correct concrete-type destruction
and unchanged small-profile allocation counts. These are structural and
correctness checks, **not evidence that the timing regressions are fixed**;
performance qualification remains unmet, including preservation of the
ordinary K3 speed gain.

At source `6245943`, the follow-on
`wirehair.wh2.small-isolation-preserved-cost-r1` completed **CONTROL_FAIL**.
Both native workers completed normally with empty stderr and no codec or
observer errors. Two of 320 same-code controls failed the frozen equivalence
bounds, one in each library load order. Both were old-library opt-in K3
two-byte decoders; both confidence intervals included one. The complete run
therefore provides no qualified speed-restoration or regression result, even
though individual treatment comparisons have directional estimates.

The exact Python 3.8 replay and a separate standard-library audit reproduced
all 103,680 records and 480 statistical decisions and verified all 721 current
receipt pins before source advancement. The historical diagnostic-log loss
above remained explicitly disclosed and checked, not silently repaired.
The spent immutable bundle is `/var/tmp/wh2-small-isolation-preserved-cost-r1`,
analysis SHA-256
`a0970624f8b624498048f1af7c9c5647d3d58d96a9f5c015509970aa68b04933`.
No samples were discarded, pooled to rescue a failed control, or remeasured.
The preceding isolation R0 namespace is separately spent **INVALID**: its
worker read a stale claim-file path and stopped before any codec work. R1
fixed and positively tested that launch binding without changing the workload.

## Retired equation profile identifiers

The identifiers `e161ce5d456f9bb7` and `20a4f27a870612a2` are permanently
reserved tombstones for rejected experimental contracts. They are unsupported,
have no public selector constants, and must never be assigned to another
equation system. For otherwise well-formed descriptors and API arguments,
serialization, parsing, validation, explicit encoder selection, and codec
creation return `WirehairV2_UnsupportedProfile` for either identifier without
writing a descriptor or publishing a codec handle.

## APIs and errors

`wirehair_v2_encoder_create()` prepares owned equation state, chooses the
profile (and deterministic seed attempt where applicable), and returns the
serialized descriptor. A null or short descriptor
buffer returns `WirehairV2_BufferTooSmall`, reports the required 32 bytes, and
does not create a codec. Its descriptor-size output pointer is required.
`wirehair_v2_encoder_create_profile()` recreates an encoder under an existing
descriptor. In both forms, the message pointer must provide at least the exact
message byte count supplied directly or recorded in that descriptor; the
implementation copies those bytes before returning.

`wirehair_v2_encoder_create_profile_id()` performs the same operation for an
explicit supported profile ID. `WIREHAIR_V2_PROFILE_CURRENT` deliberately
remains an alias for `WIREHAIR_V2_PROFILE_CERTIFIED_2026_07`; explicitly selecting
either name continues to emit the original GF(256)-only equations byte-for-byte,
including at K3. The ordinary selector is a dispatch policy, not this constant.
The corresponding C++
`Encoder::Create(profileId, ...)` overload provides the same explicit
selection. Unknown or retired IDs return `WirehairV2_UnsupportedProfile`
without falling back to the current profile.

The three additive `*_with_options()` encoder constructors accept an exact
versioned `WirehairV2EncoderOptions` record. The default initializer selects
`WirehairV2EncoderSource_Independent`, which has the same source-independent
lifetime contract as the original constructors: after successful return the
caller may change or release its message. A zero-initialized record selects no
policy and is rejected, as are unknown policy values.

`WirehairV2EncoderSource_BorrowedImmutable` is an explicit opt-in local storage
policy. Construction still completes the same eager full solve before success,
but the encoder then retains the caller's exact message range without owning it
or adding a copy beyond the independent constructor's prepared state. The
caller must keep that range readable,
allocated, and byte-for-byte immutable from constructor entry throughout the
call. Failure retains nothing and ends that obligation; after success it
continues until successful `wirehair_v2_encoder_detach_input()` or encoder
destruction. While attached, systematic packet IDs copy only their meaningful
source bytes; repair IDs use the existing solved-intermediate evaluator and do
not read the source. Detach is allocation-free and idempotent, and later
systematic packets use equation evaluation while remaining byte-identical.
Operations on one codec are externally serialized: encode, detach, free, and
source mutation must not race, including with another codec operation.

Source-storage policy is not serialized, does not change the profile ID or any
descriptor byte, and is not visible to a decoder. Independent and borrowed
encoders therefore use the same selected GF(256)-only equations, selected seed
attempt, solved intermediate state, systematic bytes, and repair bytes. Unknown
or retired equation profiles remain rejected for every storage policy.

`wirehair_v2_encode()` reports `WirehairV2_BufferTooSmall` and the exact
required packet size without modifying a short non-null output buffer.

`wirehair_v2_decoder_create()` takes only the serialized descriptor. It derives
all internal profile state from the record before accepting packets. Encode,
decode, recover, and free use the separate opaque `WirehairV2Codec` type, so V2
handles cannot be confused with legacy `WirehairCodec` handles.

Parsing distinguishes invalid magic, unsupported encoding version, invalid
size, nonzero reserved fields, unknown equation profile, and invalid dimensions
through stable `WirehairV2Result` values. Codec failures additionally report
need-more, bad seed, resource exhaustion, OOM, and unsupported platform.

The installed C++ header `<wirehair/wirehair.hpp>` provides move-only RAII
`wirehair::v2::Encoder` and `wirehair::v2::Decoder` wrappers plus a fixed-size
`SerializedProfile`. Its three explicit `Encoder::CreateBorrowed()` overloads
retain no hidden source owner, and `Encoder::DetachInput()` releases the same
native lifetime obligation. Failed replacement preserves the prior encoder and
its borrow; moves transfer the handle and its obligation. The wrappers use the
same C ABI and byte contract.

## Version migration

Encoding version 1 accepts only a declared and supplied size of 32 and requires
all reserved bytes to be zero. A future record that changes semantics or adds
fields must use a new encoding version and receive explicit parser support;
older libraries return `WirehairV2_UnsupportedVersion`. Reserved bytes cannot
be repurposed while retaining version 1.

Any change to a frozen equation rule requires a new equation profile ID and new
golden packet vectors. A new equation profile may continue using the 32-byte
encoding when the existing fields suffice. Decoders must reject unknown IDs
instead of guessing a compatible profile. Migrations recover and authenticate
under the old supported profile, then re-encode under the new profile.

## Relationship to legacy profiles

V2 and legacy share the same profile-ID vocabulary: the first 64 bits of
SHA-256 over a documented canonical name, exposed as a `uint64_t` compatibility
identifier. Both require an application integrity/authentication layer and
reject unknown named contracts.

The representations intentionally differ. Legacy packets retain their
historical unframed API, so `WirehairWireProfile` is a 16-byte in-process
selector and the application chooses framing and byte order. V2 was not yet a
public wire contract, so it starts with a canonical endian-stable record that
also carries message dimensions and the selected attempt. Legacy and V2 profile
IDs are different equation namespaces and are never interchangeable.
