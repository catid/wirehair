# Explicit small-block codecs

`wirehair/wirehair_small.h` provides the opt-in `wirehair_small_*` C API in
the normal static and shared libraries. Initialize with `wirehair_init()`.
Currently only exactly three source blocks (K3) are supported. WH1 and K6 never
select this facade. The development-branch ordinary WH2 K3 admission reuses
these equations through a distinct WHV2 descriptor and owned prepared basis;
see `V2_WIRE_PROFILE.md`. It preserves WH2's different ownership/detach rules
and requires its own speed/recovery gates. The results below qualify this
opt-in WHK3 facade only. Handles from different APIs are not interchangeable.

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

The integrated static library passed its own full-lifecycle timing gate,
`wirehair.wh2.k3-production-cost-r0`, at source `ffa7739`. All 54 same-code
timing controls and all 36 comparisons with actual WH1/current public WH2
passed. The run retained 2,488,320 fresh codec lifecycles and checked every
output. Observed time reductions versus WH1, spanning both measurement orders:

| Block bytes | Encoder | Low-ID decoder | Distant-ID decoder |
|---|---|---|---|
| 2 | 64.4-64.7% | 95.2% | 90.6-90.7% |
| 64 | 72.4% | 94.9% | 91.0-91.1% |
| 1280 | 75.3-75.6% | 93.8% | 92.2% |

These are borrowed-immutable, full-three-block measurements on one GFNI-capable
host. Encoder time includes create, descriptor output, 18 packets and free;
decoder time includes create, feed through first success, recover and free.
Both decoder streams needed three packets for every arm. This is separate
evidence from the earlier serialized-prototype timing result, not an inherited
qualification. It does not establish partial-tail, all-K, shared-call,
cold-start or non-GFNI speed, and does not change the default profile.

The immutable outcome is `/var/tmp/wh2-k3-production-cost-r0`; its raw stream
SHA-256 is `5bed150f33d02e7b8cf0db08225d695637284002f4adb4ad3ec64c961201b6c6`.
A Python 3.8 full replay verified all 602 receipt pins; a separately written
raw chronology, API-ledger and confidence-interval audit reproduced all 90
decisions. The namespace is spent: no rerun, filtering or rescoring.
The timing run adds no recovery-rate sample.

## Retained-cohort recovery comparison

At source `1861689`, `wirehair.wh2.k3-recovery-controls-r0` replayed the same
6,144 frozen loss traces through the integrated K3, WH1 and current public WH2
APIs. Observed failures with zero extra packets:

| Codec | Failures / traces | Failure rate |
|---|---|---|
| Integrated K3 | 0 / 6144 | 0% |
| WH1 | 14 / 6144 | 0.23% |
| Current public WH2 | 269 / 6144 | 4.38% |

K3 introduced no failures against either control in this cohort. All 14 WH1
failures recovered with one extra packet. Current WH2 had 5 failures after
one extra packet, 1 after two, and none after three. Each of the twelve
block-width/loss-schedule cells retained all 512 traces; K3 had no failures
in any cell. These are observed results on the retained sample, not a new
independent holdout or a universal recovery guarantee.

The 72 hard traces and 53 historical cases at their original block widths
are separate from that denominator. K3 passed all 72 hard traces at zero
overhead, versus 1 WH1 and 5 current-WH2 failures. K3 and WH1 recovered all
53 historical prefixes; current WH2 remained unresolved on 43. No packets
were appended to rescue those prefixes.

Every encoder was freed before creating its receiver. Every K3 packet matched
an independent polynomial payload oracle; every successful decode recovered
the original message twice with buffer guards. Native, portable-backend and
ASan/UBSan builds agreed on all 6,269 cases, including every descriptor,
packet hash, feed status and first-success count. Backend replays do not
increase the sample size or establish portable-backend speed.

The immutable bundle is `/var/tmp/wh2-k3-recovery-controls-r0`; native raw
SHA-256 is `60a5f0d502d2fe4b49917d47e5ad75d7be80580d878ebdf588bff0fb1fcb613c`.
A Python 3.8 replay verified all 636 input pins and all backend records;
a separately written raw-count and paired-outcome audit reproduced the
results. This namespace is also spent. All-K construction-seed validation
and broader performance/recovery coverage remain outstanding.

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
