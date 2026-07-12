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
| 28 | 1 | seed attempt | deterministic attempt `0..255` |
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

## Supported equation profiles

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

`WIREHAIR_V2_PROFILE_MIXED_2026_07` has numeric ID
`e161ce5d456f9bb7`. The ID is the first 64 bits of SHA-256 over this exact
canonical UTF-8/ASCII input, with no newline:

```text
wirehair:v2:precode-v3-mixed10gf256-2gf65536-even:packet-v4:certified-2026-07
```

The full digest is
`e161ce5d456f9bb748d7d221fa9bcb702c8dd64130edaafdcadc083178eb96e5`.
This opt-in profile retains packet-row contract 4 and the certified binary
staircase/dense geometry, but replaces the 12-row completion contract with:

- ten periodic Cauchy rows over the embedded GF(256) subfield;
- two periodic rows over `GF(256)[u]/(u^2 + u + 32)`, with each adjacent
  low/high payload-byte pair representing one extension-field element; and
- precode contract 3, which binds this mixed coefficient family and field
  representation.

The mixed profile therefore requires a positive even block size. Serialization,
encoder creation, and decoder creation reject odd block sizes as
`WirehairV2_InvalidDimensions` before writing descriptor output or publishing a
codec handle. The coefficient period is 244 and the exact generator,
coefficient exponents, arithmetic, seed stepping, and all other equation rules
are part of the named contract.

In exact coefficient notation, a `uint16_t` field element stores `a + b*u`
with `a` in the low byte and `b` in the high byte. The extension generator is
`g = 266`. For residue `m = c mod 244`, completion rows `r=0..9` retain the
embedded subfield coefficient `inv((12 + m) XOR r)`. Extension rows `r=0,1`
use `1 / (g^m XOR g^(1000+r))`. The pinned coefficient goldens are decimal
`34916` (extension row 0, residue 0), `2472` (row 1, residue 0), and `59155`
(row 0, residue 243).

The mixed completion rows improve the rank floor that dominates rare recovery
failures in the GF(256)-only profile. They also add two planar extension-field
rows and low/high conversion work during residual solving. Implementations keep
the first ten rows on batched GF(256) kernels and solve only the at-most-12-wide
binary quotient over GF(65536), bounding the extra CPU and scratch cost. Exact
paired recovery, throughput, and memory measurements belong to the release
certification record rather than the wire contract; applications should
benchmark their own block-size and loss distribution.

### Non-normative July 2026 certification snapshot

The production solver was screened with 100,000 common deterministic packet
schedules per profile at two-byte blocks. At zero packet overhead, certified
versus mixed failures were 389 versus 20 (19.45x reduction) for K=1,000 with
10% loss scheduling, 434 versus 29 (14.97x) for K=1,000 repair-only, 470 versus
68 (6.91x) for K=10,000 with loss, and 434 versus 34 (12.76x) for K=10,000
repair-only. Both profiles selected attempt zero in these cases. One extra
packet reduced mixed failures to zero except for 2/100,000 K=10,000
repair-only schedules; those remaining binary quotients exceeded the 12-row
completion width rather than failing in the extension field.

Controlled alternating whole-codec medians on the reference x86-64 machine
showed the mixed profile effectively neutral at K=1,000 and 1,280-byte blocks
(create, repair encode, and repair-only decode within 1%). At 100 KiB and 1
MiB blocks, repair encoding remained within 0.2%, while create/decode were
about 6.5-7.3% slower. Persistent intermediate and cold-receive capacities
were identical; separate-process peak RSS for the three-case benchmark was
1.6% higher for mixed (433,628 versus 426,908 KiB). These figures are
reproducible with `wirehair_v2_precode_solve_test --recovery-benchmark 100000`
and
`wirehair_v2_precode_roundtrip_test --benchmark-{certified,mixed}` and are not
part of the compatibility contract.

## APIs and errors

`wirehair_v2_encoder_create()` copies the message, chooses the deterministic
seed attempt, and returns the serialized descriptor. A null or short descriptor
buffer returns `WirehairV2_BufferTooSmall`, reports the required 32 bytes, and
does not create a codec. Its descriptor-size output pointer is required.
`wirehair_v2_encoder_create_profile()` recreates an encoder under an existing
descriptor. In both forms, the message pointer must provide at least the exact
message byte count supplied directly or recorded in that descriptor; the
implementation copies those bytes before returning.

`wirehair_v2_encoder_create_profile_id()` performs the same operation for an
explicit supported profile ID. `WIREHAIR_V2_PROFILE_CURRENT` deliberately
remains an alias for `WIREHAIR_V2_PROFILE_CERTIFIED_2026_07`; existing callers
and `wirehair_v2_encoder_create()` continue to emit the original GF(256)-only
equations byte-for-byte. Mixed encoding is an explicit opt-in through the new
selector (or the corresponding C++ `Encoder::Create(profileId, ...)` overload).
Unknown IDs return `WirehairV2_UnsupportedProfile` without falling back to the
current profile.

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
`SerializedProfile`. It uses the same C ABI and byte contract.

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
