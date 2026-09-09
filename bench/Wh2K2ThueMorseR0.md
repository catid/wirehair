# K2 noncommuting GF256 short screen

Protocol `wirehair.wh2.k2-thue-morse-r0`; sole namespace
`/var/tmp/wh2-k2-thue-morse-r0`. Mathematical feasibility only: no native
codec, timing, WH1 comparative recovery, new profile, defaults or all-K claim.

This is not the retired piecewise low-ID/avalanche high-ID projective mapping
or literal modulo mapping. It uses two noncommuting 2×2 matrices over the
existing byte field, ordered by Thue-Morse bits. Its base feedback `(2,3)`
comes from `(x+1)(x+2)` in GF(256), polynomial `0x14d`. The second feedback
is `(2 XOR lambda,3)`. Select the first lambda in 1..255, excluding 2, that
passes every two-column minor of the six local columns for all ten valid
length-four binary factors: 150 minors. Selection never sees loss/history.
Exhaustion or any later failure stops this finite experiment without retuning.

The common companion first column makes the first two packets systematic.
Invertibility plus all local minors provides a six-packet-window independence
certificate. The two matrices must actually fail to commute. The existing
generic dyadic mapper provides a 7,168-byte lookup, independently checked
against full prefix products and literal sequential products at IDs 0..2048.
No production core is modified or assumed to support K2.

Before fresh testing, require all 450 minors in 30 six-ID seam windows:
starts `2^e-4`, e=3..31, plus `UINT32_MAX-5`. Also require rank two at OH0 in
all 72 main-contract hard traces (six training/validation roots × three
widths × four schedules), all 13 original-length inventory failure prefixes,
and fixed legacy witness pairs `(265,270)`, `(2,966)`, `(1056,1313)`.

The authenticated inventory supplies 56 failure origins, preserving original
widths/tails and all 64 inventory roots. Its COMPLETE SHA256 is
`93a9068986dc24e1caa2931f5be897022e34f13b8f5ca1170a294d87365c5f2c`.
Origin projection SHA256:
`74361b5ac2588afd50a58e9a0ccf4c6b10f5f464e577f7c3d3c9daaeaa354516`.
Width-preserving prefix projection SHA256:
`2e0968c8ffb53546547cb4034c6a9c4c0d247981f6a18825f44b258a2d926a8b`.

Only after those gates pass, evaluate 512 fresh roots from the first 16 hex
digits of SHA256(protocol + `:fresh/` + decimal(i)), i=0..511. Exclude all
inventory/main roots and the prior K3/K5/K6 fresh roots. A collision is INVALID,
not replaced. The unchanged main loss law uses K2, B2/64/1280,
IID10/burst50/adversarial50/repair-only50, six delivered IDs, and a 67,072
candidate-ID cap. Retain every OH0..4 rank. Require ≤1% OH0 failures overall
and separately in each of the twelve 512-trace cells.

In the same fresh phase, strides 255, 257 and 65537 each receive 512 distinct
ordered pairs. Anchor i is the first 16 SHA256 hex digits of protocol +
`:stride/` + decimal(stride) + `/` + decimal(i), reduced modulo
`2^32-stride`; the pair is `(anchor, anchor+stride)`. Record every determinant
and rank; require ≤1% deficient pairs per stride. This is a finite nonperiodicity check, not global
injectivity: a GF256-linear K2 codec has only 257 projective directions.

Use the existing one-shot controller unchanged: worker wall60 seconds,
AS512 MiB/core0, outer70 seconds, raw4 MiB/stderr1 MiB/aggregate8 MiB.
Tests use unrelated matrices and synthetic outcomes, never the frozen
candidate. Commit/push before receipt and the sole launch. Complete exact and
independently implemented arithmetic, trace, rank and provenance audits before
HEAD/report advancement. Never rerun, filter, resize or reselect after results.
