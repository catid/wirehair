# WH2 small-band timing screens — 2026-07-28

These are the immutable raw receipts and exact runners for the initial
`wirehair-za5v` throughput screen and its local S/D2 neighbor screen.  They are
screening evidence only: neither archive is the all-K, six-schedule selection
or the independent holdout required for promotion.

These archives use historical receipt protocol
`wirehair-v2-bench:bandtiming:dispatch-v1:v1` and native row schema
`wirehair.wh2.bandtiming.v1`.  Later source review found that v1 direct-solve
timing was repeatable but did not independently commit the direct intermediate
to the encoder intermediate.  The 600 receipts remain exactly replayable for
historical screening, but are explicitly non-promotable.  Promotion-grade
reruns use the witnessed v2 protocol and must repeat the complete protected
selection before any holdout.

## Frozen completion screen

- Raw archive: `frozen_screen_raw.tar.gz`
- Archive SHA-256:
  `4d03f4af545113cc8c8c8e8bcaff73886cf025270a33b4ccb0b2d7de03a7e85f42`
- Exact runner: `run_frozen_screen.py`
- Runner SHA-256:
  `67295019b6eaeffe8a7820a21c763987f05f3099123dc3cd7a934148c846ff94`
- Complete jobs: 220/220, with 79,200 authenticated native rows and 660
  independently replayed contrasts.
- Frozen manifest, summary, and completion SHA-256:
  `8c5608fc2dae429df5a072c5c8e7bd1d19ac13919a58a9a5cb08edf1b830feac`,
  `5fca0907b9dd414909bd9a7cb1f92a26ffea7355b493fc2625febf39d765e7d7`,
  and
  `729454bb5cb5fd6f0b0159bf86ad3a6da24d6050293925b863e7c84ed9a4a5cc`.

Pure eight-row GF(256) completion was faster than pure nine-row GF(256) in all
60 paired comparisons.  Their only weak construction was the same
control-shared K=5 seed; there was no candidate-only weakness in this screen.
The cell-balanced pure8/dispatch ratios were approximately 0.811 encoder,
0.778 decoder, and 0.766 direct solve.  Tracking-X did not improve the result.

## Local S/D2 neighbors

- Raw archive: `local_neighbors_raw.tar.gz`
- Archive SHA-256:
  `939252436007f0d1c6119a535a05e5558c0e68eba61fb6ffb2168bbb665bd063`
- Exact runner: `run_local_neighbors.py`
- Runner SHA-256:
  `3140359a60e70a9abc2aaa7eb932f30cb28d00e9c8c56e7b46592ea8edb6a297`
- Complete jobs: 380/380, with 4,171 authenticated thermal rows, no read or
  EDAC errors, and 100% CPU busy throughout.
- Frozen manifest, summary, and completion SHA-256:
  `e6b8795e99bbe9e78fb078ec0bd144deb68c61106eaffe7a1d08aea6dd53999e`,
  `0a77146c5a450c4556650daaa62fc9a70570494ad914556ecffc8899821d5129`,
  and
  `9f40d230e07eeb122bbfcf651e65321f84a366e189f9414438621452498b69e9`.

The leading cell was pure8 with `S = SmallBandStaircaseCount(K) - 1` and
`D2 = 3`.  Across the 20 K/cache cells its geometric mean ratios versus
dispatch were exactly 0.754434441 encoder, 0.728543396 decoder, and
0.711185307 direct solve, with zero candidate-weak units.  The only failure
was the same control-only K=5 construction duplicated by cache state.

The predeclared large-selection roster is:

1. pure8, S-1, D2=3 — timing leader;
2. pure8, S+0, D2=3 — extra staircase hedge;
3. pure8, S-1, D2=5 — denser-binary hedge;
4. pure9, S-1, D2=3 — extra-GF(256) hedge.

No per-K or per-seed repair was applied.  Architectures are compared on their
raw weak-seed populations; any seed repair remains deliberately deferred until
after the architecture is selected.

## Historical superseded v1 all-K checkpoint

A first protected execution completed its 2,376 native jobs but was rejected
by the final authenticated replay before any summary, completion record, or
decision was published.  It is not selection evidence.  Result-blind
forensics found that the runner had incorrectly treated the direct dispatch
overhead as candidate-independent even though the native direct panel uses the
pair-local maximum candidate/dispatch recovery prefix.  The corrected runner
excludes only that overhead, retains the invariant direct result class, and
requires a completely fresh promotion-grade execution.

The following hashes identify only the final historical v1 runner checkpoint.
They are retained for archive provenance and are not the current
promotion-grade runner:

- campaign runner SHA-256:
  `72ee3a1f1576304fa66ebcd084045f54cff45f0ca32dd9f88e27144772739871`;
- campaign runner tests SHA-256:
  `af096019e98d13fca634dcf5f887e9b2e36024376517dec8c4f8a221bc813d6a`;
- strict band-timing parser SHA-256:
  `66a1c7e83149914b8063d03eafac771d9ff9edbaf3be084f30557d1938fd92c2`;
- strict parser tests SHA-256:
  `27f3ca2761126630b985846948fefd1a7740417a6fc298ab8f1ca61ab92002ae`;

## Current witnessed v2 prelaunch checkpoint

The promotion receipt protocol is
`wirehair-v2-bench:bandtiming:dispatch-v1:v2`, the native row schema is
`wirehair.wh2.bandtiming.v2`, and the campaign schema is
`wirehair.wh2.za5v.all-k-campaign.v2`.  Successful WH2 encoder, decoder, and
direct rows bind the same solved-intermediate byte count and SHA-256 witness;
historical v1 receipts are replay-only and cannot satisfy this contract.

The current promotion artifacts are:

- benchmark binary SHA-256:
  `84c6d4fd1ab410ffee55db68ab4decd65fd498ee37e8d1de3d1d0362697e5c4d`;
- native benchmark source SHA-256:
  `db9e1b11bf809b958a9b9869fc55d6a62addb067df0779c5401eac5fd36138b0`;
- native CLI test SHA-256:
  `e41501899e6b82c2b9d45262ce439f77e72401f4dda6e05e361e18303e799775`;
- strict band-timing parser SHA-256:
  `ceeebe3cab0507fdfe86db9067d2e1aea1741aaeb1629cdcaafb589bb568e519`;
- strict parser tests SHA-256:
  `c9921d227d4ce535ae27cbcd696ca8d80fe35590c05a03e5f31e7d0318f7f02b`;
- paired-context tool SHA-256:
  `07cd2dd509e0847caaaab813ba2d7fe6f4cf57a4fe2d78e2c9400c2545f87c82`;
- campaign runner SHA-256:
  `bae4336508d9d5fd9dd514950da709cb01cc129dcf9e0b36f3e8fe96e17346ef`;
- campaign runner tests SHA-256:
  `1a5d2378c11bf2f69998dedf43204f0493fdf47f8082c61feae9d266ec25ec95`;
- exact result-free selection plan SHA-256:
  `3a0413150506aa777c995774e80e7ca89e2dc34f97801d8089c8e085697807c0`;
- frozen-roster SHA-256:
  `d7f032731d3f5b1a52d37cb27cc4236e92ef8be46efa5445934c8356288421fc`;
- selection job-list SHA-256:
  `021dcb0607d83e5ce5da1d4ebbafbab43f1dfcef34209b64326a877fe97e2b63`;
- matched-cell-set SHA-256:
  `6cc049811a0826e1be9fa750497ba03028fd263a2dc3328f9409b8cbaf9a8cdf`;
- selection-policy SHA-256:
  `bd7f3bb0497e141af20e3c3189006702b77e832ddfd3581475a0303bd589780c`.

Selection contains 2,376 jobs: four candidates, every K=2..100, and all six
loss schedules.  Each candidate receives 768 matched raw replicates per
K/schedule, or 456,192 byte-recovery cells.  The four possible holdout job-list
hashes were also fixed in advance; each holdout has 594 jobs and 256 disjoint
replicates per K/schedule:

- pure8, S-1, D2=3:
  `621b29e4c897bc47dd4497dc2e990fedbf2236333906d36df77b4ace6c4560b5`;
- pure8, S+0, D2=3:
  `bb15989bcca31a956100cdae64e2ee281d36ca341eca814f0b49b624dd293e4f`;
- pure8, S-1, D2=5:
  `04b32e54498252315c74a1d64270d1669cc89fa308ebbbff230353e6bd4466e9`;
- pure9, S-1, D2=3:
  `26b0b72a4868aafac50b0f1071e5fe856631d5aa69fe4736bdd41a796511a672`.

The result-free policy computes the raw six-dimensional Pareto surface from
final unrecovered count, the complete overhead-tail area, unique weak
constructions, direct-solve cost versus dispatch, and encoder/decoder cost
versus WH1.  Weak seeds are descriptive and are neither repaired nor used as
an eligibility veto.  Timing points are arithmetic means of equal-weight
K/schedule job means.  Confirmation requires at least half of the raw
replicates in every cross and both corresponding A/A arms, a Student-t 95%
interval across job means, and an aggregate upper bound that also covers the
mean of the per-job interval highs and the mean effective timing floor.

The runner derives the sole eligible survivor itself and cryptographically
binds that decision into the completed selection evidence; the holdout path
rejects a caller-supplied alternative or a no-survivor selection.  Selection
and holdout reliability blockers are cumulative.  Holdout reports recovery
and direct-speed confirmation separately from the stricter production
requirement to beat WH1 encoder and decoder throughput.  Any raw recovery or
throughput gap retains dispatch; nonrepairable candidate-only regressions
require investigation, while only constructor-weak blockers are eligible for
the deferred seed-repair step.

All 600 archived screen receipts were independently replayed against the
enriched parser before this freeze.  All 1,800 legacy replicate records
re-derived their exact candidate, dispatch, and WH1 construction result/class;
all 600 native stream hashes were unique.  The canonical-JSON SHA-256 of the
sorted stream-hash list is
`179367f473a787b15141100a29f6df525061f9fa18c5f2abf23cbb75bc085dd9`;
the equivalent newline-delimited list with a terminal newline hashes to
`f51ed68d80706ad1b1a0c65fe24f262ca06f5e1d5efb888fc6c8bf431e7b5c35`.
