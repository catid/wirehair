# Resume notes — WH2 peel-training audit

Audit date: 2026-07-27
Branch: `feat/wh2-opcount-cost-model`
Primary Beads: `wirehair-d45d`, `wirehair-68o1`, `wirehair-tz5b`

## Stop: the inherited peel tables are not production evidence

Do not deploy, interpolate, benchmark as “trained,” or seed-tune from any of
these legacy artifacts:

- `tools/peel_table.json`
- `tools/peel_table_small.json`
- `tools/peel_table_small_validated.json`
- `tools/peel_table_bigk.json`
- `tools/peel_table_final.json`

The previous `resume.md` described `peel_table_final.json` as a validated
137-K deliverable with median `+17.4%` goodput. That conclusion does not survive
replay. The table is retained only as invalidated historical evidence; the
hardened tools reject its legacy schema.

The audit found:

1. The Python proxy used an analytic degree-one probability of `1/K`.
   The codec actually ships an approximately `1/128` degree-one prefix through
   K=2048 and no degree-one prefix above that range.
2. Proxy and real-codec stages did not interpret all coordinates identically.
   In particular, the funnel and validator ignored the staircase degree scale,
   while 131 of 137 final entries did not record one at all.
3. Search and validation did not carry enough provenance to replay a result:
   no independent seed domains, no complete tail/failure metrics, no executable
   identity, and subprocess failures could be mistaken for partial success.
4. The checked-in cost-model table mixed calibration regimes. Its K=128 row
   had been replaced with a 64-KiB/GF(2^16)-enabled measurement inside a table
   whose other rows and receipt describe 1280-byte, zero-GF(2^16) solves.
5. The pushed branch was not self-contained. Several source/header changes and
   `WirehairV2EsCostModel.inc` existed only as uncommitted handoff files, so a
   clean checkout could not build the benchmark.

No gain percentage from the inherited tables should be quoted as a current
WH2 result. There is presently no production peel-parameter table.

## Repairs made by this audit

The training pipeline now uses one versioned, fail-closed schema and a shared
codec helper:

- `wirehair_v2_bench peelpmf --N K` exports the native PMF and structural
  metadata instead of asking Python to transcribe the codec.
- Native metadata, exact PMF identity, parameter semantics, executable
  identity, independent search/validation seeds, and complete metrics are
  recorded and checked.
- The true shipped control is now a run with no peel or staircase override.
  Passing an identity-looking vector through a test-hook path is not treated
  as the shipped arm.
- Search receipts store the exact binary64 PMF that was measured.  Exact
  anchors replay that stored PMF; shorthand coordinates are descriptive and
  are not claimed to be an exact cross-language reconstruction protocol.
- Writes are atomic. Partial subprocess output and nonzero exits are failures,
  not candidate records. The source table, executable, and source tree are
  rechecked as stable snapshots before publication.
- Legacy/unversioned tables and incomplete entries are rejected.
- The unverified counter-cost proxy is disabled by default. It requires an
  explicit experimental opt-in until its raw calibration data and generator
  are reproducible.
- Sweep warm starts are bidirectional from an explicit pivot rather than
  accidentally depending on traversal direction.
- The shipped arm is the safe default at unmeasured K. Parameter interpolation
  or clamping requires a separate explicit experimental opt-in.
- A selected arm with any recovery failure is refused, strict JSON rejects
  duplicate keys and every non-finite spelling (including exponent overflow),
  and funnel deduplication keys both PMF and staircase scale.

The codec/benchmark bug hunt also repaired:

- strict transactional parsing for staircase shape, pinned row degree, and
  degree-scale hooks;
- legal zero-degree shaped/pinned rows and a deterministic O(S) bulk repair
  with no arbitrary `4096*S` failure cap;
- rejection of contradictory balanced-plus-shaped configurations;
- missing degree-scale equality/default/fuzz coverage;
- the peel-CDF exact-boundary off-by-one (`lower_bound` versus `upper_bound`);
- variable-depth override sampling, replaced by a canonical 64-bin,
  six-comparison decision tree;
- malformed peel PMFs silently falling through to another arm instead of
  remaining active-invalid and failing closed;
- cost-only dispatch so width zero follows the table's 1280-byte algorithm
  paths rather than unrelated 64-KiB paths;
- silent aliasing of ES cell identifiers at or above 2^32;
- ignored `precodecost --threads` (parallel counter harvesting now works only
  at payload width zero; timed runs remain serial);
- thread-order-dependent mixed-core configuration. Reused workers now pass
  through a canonical valid bridge before applying each target, and exact
  receipts are invariant under forward/reverse task order and 1/4 threads;
- accepted-but-unusable `precodecost` cells, unchecked pool-size arithmetic,
  invalid overhead/mix ranges, inverted ES boxes, and excessive ES thread
  counts;
- successful ES runs with no usable cost, non-finite objectives from extreme
  finite scales, and NaN fail-model predictions that `std::max` had converted
  into an apparent zero failure rate;
- ignored `--emit-noise` close/write failures;
- ordinary CTests inheriting `WIREHAIR_V2_*` experiment hooks from the
  developer shell. Every normal V2 test now starts from a hermetic hook
  environment, and the strict `BUILD_TESTS=OFF` gate verifies that production
  benchmark compilation contains no test-hook macro;
- packed mixed-field coefficients for experimental GF(256) row counts below
  eight. The old cache silently substituted a 9- or 12-row layout, shifting the
  following GF(2^16) lanes and causing deterministic payload-recovery failures.
  Historic 8--12-row caches remain immutable and low-row experiments now use
  an exact active-layout cache; the solver test path also handles every
  one-through-four packed-word count;
- phased encode and pivot back-substitution paths that entered two-row
  GF(2^16) kernels when an experiment selected zero or one extension row.
  Those legal test configurations now use no extension operation or the
  one-row kernel respectively, with direct encoder and lossy-solve coverage;
- stale GF(256), GF(2^16), tiny-fast-path, and null-witness test fixtures.

The K=128 cost row has been replaced with the handoff row claimed to belong to
the 1280-byte regime, and a runtime receipt checks row count, ordering, ranges,
and total fit rows. The regime and coefficients cannot be independently
verified: the original raw observations and exact generator invocation are
still absent.

## Current validation commands

From the repository root:

```bash
cmake -S . -B build-audit \
  -DBUILD_TESTS=ON \
  -DWIREHAIR_BUILD_BENCHMARKS=ON \
  -DCMAKE_BUILD_TYPE=Release
cmake --build build-audit -j"$(nproc)"
ctest --test-dir build-audit --output-on-failure -j"$(nproc)"

PYTHONWARNINGS=error python3 -m unittest tools.test_peel_tools
cmake -E env \
  --unset=WIREHAIR_V2_PEEL_DEGREES \
  --unset=WIREHAIR_V2_STAIRCASE_DEGREES \
  --unset=WIREHAIR_V2_STAIRCASE_ROW_DEGREES \
  --unset=WIREHAIR_V2_STAIRCASE_DEGREE_SCALE \
  --unset=WIREHAIR_V2_BAND_TRACKING_X \
  build-audit/codec/wirehair_v2_bench selftest
cmake -E env \
  --unset=WIREHAIR_V2_PEEL_DEGREES \
  --unset=WIREHAIR_V2_STAIRCASE_DEGREES \
  --unset=WIREHAIR_V2_STAIRCASE_ROW_DEGREES \
  --unset=WIREHAIR_V2_STAIRCASE_DEGREE_SCALE \
  --unset=WIREHAIR_V2_BAND_TRACKING_X \
  build-audit/codec/wirehair_v2_bench peelpmf --N 128
```

The 2026-07-27 audit receipt is:

- strict GCC release: 43/43 CTests passed, including the 1,113-second all-K
  fingerprint; after the final low-H16 test-hook repair, the exact final tree
  again passed all 42 other CTests and the default H16=2 fingerprint path was
  unchanged;
- Clang 18 ASan+UBSan: the final encode/solve regressions passed for H16=0/1,
  both coefficient geometries, and constant/independent residue schedules;
- duration fuzz: 128 parallel workers, 463,477,941 final-tree mutations, all
  passed with no crash or sanitizer-style diagnostic;
- Python: 43/43 hardened peel-tool tests passed with warnings promoted to
  errors, followed by bytecode compilation and the benchmark self-test;
- two independent final source-reading passes found no additional defect after
  the low-H8 packed-layout and low-H16 pair-kernel repairs.

## What remains

1. Recreate the operation-cost model from checked-in raw observations and a
   deterministic generator, with separate receipts for every payload width and
   completion regime used. Until then, prefer direct real-codec ranking
   (`wirehair-znwn`).
2. Train a new table from scratch with the hardened schema. Never import or
   “fill in” fields from the invalid legacy JSON.
3. Validate candidates on independent seeds, hard/random loss families, the
   intended payload widths, and all relevant K. Record raw weak-seed counts and
   overhead distributions for architecture comparisons before seed fixups.
   Counterbalance/interleave trained and shipped timing arms before trusting
   small throughput deltas (`wirehair-lnfk`).
4. Only after the architecture and parameter table are frozen, search targeted
   seed fixups and compare WH2 against mainline for solve speed and recovery.
5. Define either a versioned native coordinate-to-PMF reconstruction or consume
   the stored exact PMF directly, then gate a future production handoff on an
   exact cross-language grid (`wirehair-2qlw`).
6. Finish the separate ES-trainer repair (`wirehair-d45d`).  The historical
   `meta_es.py` and `es_degree.py` named by that Bead are not present in this
   checkout, so their invalid meta-grid and unconstrained degree-shape results
   have not been repaired and must not be reused.  A replacement trainer must
   score its final candidate, constrain recovery on independent held-out seeds,
   handle ties and both-fail pairs, and either condition inactivation cost on
   successful solves or assign failed solves an explicit penalty.  The current
   `precodefail` report stores `InactivatedColumns` for every result at
   `WirehairV2Bench.cpp:9538-9539`, accumulates it at line 9801, and divides by
   all trials at line 9867; that all-trial statistic is not by itself a safe ES
   objective.
7. Repeat full-codebase bug hunts after each substantial codec change; a clean
   test run is necessary but is not a source review.
