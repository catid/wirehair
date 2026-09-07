# Explicit K6 codec

`wirehair/wirehair_k6.h` provides the opt-in `wirehair_k6_*` C API in the
normal Wirehair static and shared libraries. Initialize with `wirehair_init()`.
Existing WH1 and WH2 functions never select this codec; their wire formats,
source policies, and detach guarantees are unchanged. Handles from these APIs
are not interchangeable.

## Qualification and scope

This is a GF(256)-only codec for exactly six source blocks, including a partial
last block. The equations are the fixed K6 Thue-Morse pair selected in the
recovery experiment, with feedback vectors `(124,127,152,84,241,63)` and
`(125,127,152,84,241,63)`. There is no run-time seed selection or construction
search. The immutable 39,936-byte lookup is included in the library; building
or consuming the package requires neither Python nor a live encoder.

The integrated static library passed its own full-lifecycle timing screen,
`wirehair.wh2.k6-production-cost-r0`, at source `0487537`. All 54 same-code
timing controls and all 36 comparisons with WH1/current public WH2 passed.
The run retained 2,488,320 fresh codec lifecycles and checked every output.
Observed time reductions versus WH1, spanning the two measurement orders:

| Block bytes | Encoder | Low-ID decoder | Distant-ID decoder |
|---|---|---|---|
| 2 | 74.3% | 90.0% | 81.5% |
| 64 | 74.8-75.0% | 89.3-89.4% | 82.4% |
| 1280 | 77.2% | 88.5% | 85.5% |

These measurements use borrowed immutable input on one GFNI-capable host.
Encoder time includes create, descriptor output, 18 packets and free; decoder
time includes create, feed through first success, recover and free. They are
fresh-handle measurements after process initialization, not process cold-start,
shared-library call-overhead, non-GFNI, or all-K performance claims.

The unchanged equation family had 8 zero-overhead failures in 6,144 frozen
recovery traces; all eight recovered with one additional packet. The selected
lookup and receiver body are preserved in this integration, with independent
packet/selector and ownership/recovery tests. The timing screen is not a new
failure-rate sample. Neither result is a universal recovery guarantee, and
this codec is not a new default profile.

## Input ownership

`WirehairK6_Independent` copies the input before encoder creation returns.
`WirehairK6_BorrowedImmutable` retains it: **all packets, including repairs,
may read the source** until successful detach or encoder destruction. All
operations on a handle, source mutation, detach and free must be externally
serialized.

Unlike existing WH2, the first borrowed `wirehair_k6_encoder_detach_input()`
allocates and copies. Out-of-memory failure preserves the original encoder
and its source-lifetime obligation. Success permits immediate source mutation
or release. Independent and already detached encoders detach without work.
Encode, decode and recover allocate no memory.

## Descriptor and packets

Descriptors are exactly 32 bytes, little endian:

| Offset | Size | Value |
|---|---|---|
| 0 | 4 | ASCII `WHK6` |
| 4 | 2 | Encoding version 1 |
| 6 | 2 | Descriptor size 32 |
| 8 | 8 | Profile ID `0x5748324b36544d31` |
| 16 | 8 | Message bytes |
| 24 | 4 | Block bytes |
| 28 | 4 | Reserved, zero |

The identity and descriptor are byte-compatible with the qualified serialized
K6 prototype, not with existing `WHV2` profiles. Unknown, current WH2 and
retired profile records are rejected, never reinterpreted. Dimensions require
`5 * block_bytes < message_bytes <= 6 * block_bytes` and
`block_bytes <= floor(268435456 / 7)`.

IDs 0 through 5 are systematic. ID 5 carries only its meaningful tail bytes;
every other ID carries one full block. Repair IDs span the remainder of the
32-bit ID space. Matching dependent packets are idempotent; contradictory
dependent payloads permanently poison that decoder. A `NeedMore` response
preserves the accepted equations so the caller can feed another packet.
Successful recovery can be repeated.

The lookup byte SHA-256 is
`27b105e1449bec190bd3c83f07feefa639cd32bc356baebfb03828ea7cbccb6d`.
The tiny GFNI kernel changes the byte basis within GF(256); the portable
fallback computes the same products. It shares Wirehair's existing field
initialization and CPU/OS feature checks. Counter builds use the existing
counted generic payload operation.

Recovery is not authentication. Protect descriptors, packet IDs and payloads
with trusted framing and an appropriate integrity/authentication mechanism.
See the installed header for exact capacities, error results and alias rules.
