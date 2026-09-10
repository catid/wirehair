# K4 noncommuting GF256 short screen

Protocol `wirehair.wh2.k4-thue-morse-r0`; sole namespace
`/var/tmp/wh2-k4-thue-morse-r0`. Prospectively frozen before code or candidate
arithmetic in `wirehair-sxvz.16.1.20.82.6`. This is mathematical feasibility,
not native codec timing, comparative WH1 recovery, an installed profile or
default promotion. Existing K2/K3/K5/K6/K8 evidence does not qualify K4.

## Motivation and finite candidate

K4 is the next untested structural dimension in the sealed inventory priority
K5, K2, K8, K4. Ordinary WH2 had 13/768 zero-extra-packet failures versus
WH1's 16/768. At overheads 1/2/3/4, WH2 had 0/0/0/0 failures and WH1 had
7/2/0/0. All six separate hard cases recovered at four packets in both arms.
These baseline observations do not prove candidate superiority. The historical
`tri-bmq4-schur-v0` experiment is K64000, not a dimension-four companion-pair
screen. Earlier generic-mapper dimension-four neutral tests used unrelated
matrices; they did not select or score this equation family.

Use GF(256), polynomial `0x14d`, primitive element 2, dimension four. The base
feedback is the low-to-high coefficient vector below the monic leading term
of `product(j=0..3)(x + 2^j)` in this field. A0 has subdiagonal ones and that
last column. A1 changes only its constant feedback by XOR lambda.

Try ascending lambda 1..255 except the base constant (which would make A1
singular). Select the FIRST pair passing every four-of-eight minor of the local
columns for each of the ten fixed Thue-Morse factors:
`0010,0011,0100,0101,0110,1001,1010,1011,1100,1101`.
Local columns are identity4 followed by `P * feedback(bit)`, initially P=I,
then P=P*A_bit for each factor bit. There are 70 minors per factor, 700 per
candidate, at most 254 candidates and 177,800 selection checks. Preserve each
rejection's parameter, count and first witness. No survivor is EXHAUSTED.
Selection never consults loss/history; no later failure permits reselection.

Packet row(id) is the first column of the right product of matrices selected
by popcount(j) parity for j=0..id-1. IDs are uint32, never wrapped. Reuse the
dimension-parametric 10/7/7/8 mapper: exactly 20,480 immutable lookup bytes,
`2*1024*4 + 4*128*16 + 256*16`. Its three dense matrix-vector products have
48 GF products and 36 nontrivial algebraic XORs. The Python helper actually
uses 48 XOR assignments because each dot starts at zero; reference-oracle and
caching work are additional. None of these counts is a timing result.

## Gates after the one local selection

Require invertible noncommuting matrices, systematic identity rows 0..3, an
independent polynomial-arithmetic/reverse-pivot rank check for all 700 selected
minors, and agreement of every accessed packed row with full dyadic products.
Additionally compare IDs 0..2048 with literal sequential matrix products.

Require all 70 minors in each actual eight-ID window starting `2^e-4` for
e=2..31, plus `UINT32_MAX-7`: 31 windows and 2,170 seam minors. Require rank
four at zero overhead on the 72 main-contract hard traces (six training and
validation roots, three widths, four schedules), and full rank on every one
of the 37 original-length historical prefixes. Preserve every seam, hard and
historical result. Any failure is terminal FAIL without entering fresh testing.
No global injectivity or universal MDS claim is made.

The six-member inventory bundle is authenticated in full, including all 3,096
cases, 774 K4 cases and all WH2/WH1 failures at overheads 0..4. Its immutable
COMPLETE hash is
`93a9068986dc24e1caa2931f5be897022e34f13b8f5ca1170a294d87365c5f2c`;
raw hash is
`dc4bbbac93cad98669f325755c42f47339eed349d3064b0cca0908819c4c5729`.
Independent pre-candidate projections agree, using canonical JSON without a
trailing newline and the same origin/width-preserving-prefix schema as K8:

| Projection | Count | SHA256 |
| --- | ---: | --- |
| All failure origins, including original widths/tails | 38 | `2efbce06e3db3928635764c51c74dd0f79966b7050b4751f592c60f14d75b578` |
| Sorted prefixes with original widths | 37 | `96cb80f4458adca9f5cc6ed773f9cca27ba9e0b5a048dc0cd8c9f1ca11be42ff` |
| Inventory roots | 64 | `6a03f629cae261c18128ba4707b2656dda7833a9f836948442f9fb92e8c20bbf` |

No prefix is extended or origin dropped. Combined projection hash is
`44621346e5d5f969ed890d76d29019a00554553ba77e57a821d3d3860d5f6e61`;
with the complete source-inventory provenance it is
`cd6a1a5edc96d717b0dacfa44d24edab2e819565dfa44130608bf5c5710fdaef`.

Only after every structural gate passes, test 512 fresh roots, the first
16 SHA256 hex digits of protocol + `:fresh/` + decimal(i), i=0..511. Exclude
all inventory/main roots and the 512 fresh roots of each prior K2/K3/K5/K6/K8
structural recovery domain. A collision is INVALID, never replaced.

Keep the original FrozenTrace law: K4, B2/64/1280, IID10%, burst50%,
adversarial50%, repair-only50%, eight delivered IDs, candidate cap 67,584.
Bursts have length eight; adversarial IDs descend from UINT32_MAX at stride
two; repair-only IDs start at four. Record all ranks at overheads 0..4 for
all 6,144 traces. Require at most 1% OH0 failures overall AND in each of
twelve 512-trace width/schedule cells. Preserve first-success histograms and
all observed rows. Width cells share roots; do not claim 6,144 independent
random seeds or treat later native/backend replays as fresh rate samples.

## Execution and audit boundaries

Reuse the unchanged one-shot capture/controller: worker wall 60 seconds,
address space 512 MiB, core zero; outer 70 seconds, raw 4 MiB, stderr 1 MiB,
aggregate 8 MiB, historical input 64 MiB. No production source is changed.
Neutral tests use unrelated matrices or synthetic outcomes, never the frozen
feedback/candidate, and exercise selection/exhaustion, all failure stops,
trace generation, per-cell thresholds, full accounting, positive/negative
claims, source closure and inherited capture/publication/deadline checks.

Before the sole launch, require clean repeated main/independent source
reviews, neutral tests on Python 3.12 and 3.8, committed/pushed exact sources
and an authenticated receipt. Worker claim checks precede selection; complete
source/interpreter/history identities are checked again after it ends. Every
launched namespace is permanently spent, including failures. No rerun,
filtering, pooling, resizing, retuning or subset rescue.

Exact retained replay and a separately implemented full arithmetic, lookup,
row, local/seam/rank, trace, history, decision and provenance audit must finish
before producing HEAD or pinned documentation advances. A survivor licenses
separate native payload and ownership-matched actual-WH1 encoder/decoder
qualification only. K8 retention, K3 preservation, installed K5 speed and the
full all-K speed/recovery/zero-bad-construction objective remain required.

## Independently audited terminal result

The sole screen at pushed source
`86acfff0ff4a7a831a5a4ba9651a208dbc1cea19` completed **PASS** on 2026-09-10.
The first local parameter, lambda 1, selected feedbacks `(64,120,54,15)` and
`(65,120,54,15)`. No second parameter was tried and no loss/history outcome
influenced selection.

| Frozen check | Result |
| --- | --- |
| Selected local minors | 700/700 full rank |
| Seam minors | 2170/2170 full rank |
| Main hard traces, OH0 | 72/72 full rank |
| Original-length historical prefixes | 37/37 full rank; all 38 origins retained |
| Fresh traces, OH0 | 9/6144 failures = 0.146484375% |
| Worst width/schedule cell, OH0 | 4/512 = 0.78125%, B2 repair-only |
| Fresh traces, OH1 through OH4 | Zero failures at every overhead |
| Packed lookup | 20480 bytes |

All nine deficient fresh prefixes become full rank with one additional
packet. All twelve width/schedule cells meet the frozen 1% threshold. The
baseline inventory used different roots, so its 13 WH2 / 16 WH1 failures
must not be treated as a paired comparison with this candidate. This is not
per-trace dominance, an installed-codec recovery rate or a universal guarantee.

The worker exited zero with empty stderr in 0.823798982 seconds. The complete
observed controller took 1.139050763 seconds and retained every output with
matching pre/post source, interpreter and inventory identities. These are
experiment execution times, not encoder/decoder performance.

Both exact-helper retained replays passed on Python 3.12 and 3.8 without
invoking a worker or extending the search. The separately written auditor
imports no experiment module. It independently checked field arithmetic,
the finite selection, all local/seam minors, the complete lookup, 2049 literal
sequential rows, all 2252 accessed rows, every one of 31080 hard/fresh prefix
ranks, full history projection, all cell decisions and fourteen source pins.
Its two complete reports are byte-identical; all inputs matched again after
each audit. A trace-only C++ program linked to the unchanged FrozenTrace
generator also agrees on every one of 6216 hard/fresh ID lists and their
packet-ID hashes under both Python readers. No native codec was executed.

All 34 neutral screen/controller tests pass under each Python version. The
independent auditor's five neutral tests and standalone unrelated-matrix
checks pass under both as well. Repeated main/independent source reviews are
clean. A replay command initially named a nonexistent `/usr/bin/python3.8`;
it executed nothing. Using the installed `python3.8` completed the read-only
check, preserving the preceding Python 3.12 report and the sole scientific run.

Retained evidence:

- Bundle `/var/tmp/wh2-k4-thue-morse-r0`; COMPLETE SHA256
  `8f1c9357b125edffc282df4cde29592130e27566f86946a7a85b0e92f175d281`.
- Raw SHA256
  `1e6422eeb1a62631fb9f0e8df86dac286b39ec7021cc5bbcb1d3f2335a351ce2`.
- Lookup SHA256
  `f4c3bc7d9fffac727acefc217d0481206b4c2d608af051993f34fb6a19f47c69`.
- Receipt/exact/native-trace artifacts `/tmp/wh2-k4-thue-receipt.kV7XJ04c`;
  receipt SHA256
  `96253f1b876b170af9bee75be946822d1240c5355a8d3b37428833c744a12b70`;
  `AUDITS.json` SHA256
  `e2f79b46a2242b2bf7ad5b8966af2289adcf1a34444fca3c57f11848851b08b6`.
- Independent auditor directory `/tmp/wh2-k4-independent-audit.y5uXMxNj`;
  `audit.py` SHA256
  `17d71ad9ed6075b4f66329d9d741c1bb6444cef0095bcc14fd7d71cbbb8fc6f0`;
  both final report SHA256
  `e3ede776af6bb83eabac1728184075e75cb411840970f7144765946e493153ae`;
  `POSTRUN.json` SHA256
  `9a0199c760ecaea524edd2e46f131351d39b41bdc7a20afacef280fcea6aece9`.

All exact and independent audits and the sealed crosslinks completed before
producing HEAD or this pinned contract advanced. The namespace is permanently
spent: do not rerun its strict current-HEAD checks after this documentation
change. Native qualification now proceeds separately in
`wirehair-sxvz.16.1.20.82.6.1`, preserving the exact selected equations. No
production/profile/default change, WH1 speed win, comparative recovery win or
all-K qualification follows from this mathematical result alone.
