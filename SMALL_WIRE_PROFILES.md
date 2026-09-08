# Explicit small-block codecs

`wirehair/wirehair_small.h` provides the opt-in `wirehair_small_*` C API in
the normal static and shared libraries. Initialize with `wirehair_init()`.
Currently only exactly three source blocks (K3) are supported. Existing WH1,
WH2 and K6 APIs never select this path; their profiles, defaults and source
guarantees are unchanged. Handles from different APIs are not interchangeable.

## K3 identity and qualification

The GF(256) equations use the fixed Thue-Morse companion pair `(8,14,7)` and
`(9,14,7)`. The immutable 13,056-byte lookup is bundled in the library.
Library-only builds (`BUILD_TESTS=OFF`, `BUILD_CODEC_V2=OFF`) and consumers
need neither Python nor local benchmark evidence. Its
byte SHA-256 is
`c78d6f350767bc5336f36eae30424914347f42b4314465b590fee6c1612e9d15`.
No run-time seed selection or construction search occurs.

The sealed equation family passed 6,144 frozen recovery traces with no
zero-overhead failures and all 72 hard traces. Native and serialized prototype
tests recovered the actual payloads on the same cohort. The integrated library
also matches the prototype on all 7,774 retained payload cases, including
48,187 packet/feed/recovery checks in native, portable-backend and sanitizer
builds. Those replays are not fresh failure-rate samples or paired WH1
recovery-rate comparisons. Installed C consumers, ownership/error tests and
package relocation are checked separately from performance.

The separate serialized prototype passed its encoder and low/distant decoder
timing gate against WH1/current WH2 at 2-, 64- and 1,280-byte blocks on one
GFNI-capable host. **Those speed results do not qualify this integrated
library.** Actual-library performance qualification is still pending. Neither
the retained recovery screen nor prototype timing establishes all-K,
shared-call, cold-start or non-GFNI performance or universal recovery.

## Ownership and operations

`WirehairSmall_Independent` copies the source before encoder creation returns.
`WirehairSmall_BorrowedImmutable` retains it: all packets, including repairs,
may read every message byte until successful detach or destruction. All calls
on a handle, detach/free and source mutation must be externally serialized.

Unlike existing WH2, the first borrowed detach allocates and copies. Failure
preserves the old encoder and source-lifetime obligation; success permits
immediate source mutation or release. Independent and repeated detach calls
allocate nothing. Encode, decode and recover also allocate nothing.

Descriptors are exactly 32 bytes, little endian:

| Offset | Bytes | K3 value |
|---|---|---|
| 0 | 4 | ASCII `WHK3` |
| 4 | 2 | Version 1 |
| 6 | 2 | Descriptor size 32 |
| 8 | 8 | Profile ID `0x5748324b33544d31` |
| 16 | 8 | Message bytes |
| 24 | 4 | Block bytes |
| 28 | 4 | Reserved, zero |

Require `2 * block_bytes < message_bytes <= 3 * block_bytes` and
`block_bytes <= 67,108,864`, bounding the decoder slab to 256 MiB. The descriptor
is compatible with the qualified K3 prototype, not with current WH2 or K6.
Unknown and retired IDs fail safely and are never reinterpreted.

IDs 0, 1 and 2 are systematic. ID 2 carries only the meaningful tail bytes;
every other ID carries a full block. Repair IDs span the remaining 32-bit
space. Matching dependent packets are idempotent; contradictory payloads
permanently poison the decoder. NeedMore preserves equations, and successful
recovery is repeatable. See the installed header for alias and capacity rules.

Recovery is not authentication. Protect descriptors, IDs and payloads with
trusted framing and suitable integrity/authentication checks.
