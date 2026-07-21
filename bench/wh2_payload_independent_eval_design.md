# WH2 Payload-Independent Graph Seeding: Cross-Width Evaluation Design

Bead: wirehair-sxvz.16.1.5.2
Branch: exp/wh2-payload-independent-profile
Status: DESIGN ONLY. Implementation is committed; evaluation runs are
deferred until the saturated campaign on this host completes.

## Background and hypothesis

WH2 equation algebra never consumes the payload width, but
`MatrixSeedFromProfile()` hashes `SeedProfile::BlockBytes` plus the
block-byte-classified peel policy, so precode and packet-graph seeds -- and
therefore recovery hotspots -- are payload-dependent at matched K.  A prior
holdout that normalized graph seeding to bb=2 (via the test-only precodefail
flag `--seed-block-bytes 2`) made K-indexed fixups portable across
bb={2,64,256,1280,4096} with 88-91% failure reductions.

This experiment promotes that normalization to a first-class experiment-only
profile: the test-hook `--payload-independent-seeding` (codec knob
`SetPayloadIndependentEquationSeedingForTesting`, canonical constant
`wirehair_v2::kPayloadIndependentEquationSeedBlockBytes` = 2) normalizes every
equation-seed derivation inside `MatrixSeedFromProfile()` while memory layout
and kernels keep the true BlockBytes.  It is NOT a public V2 profile ID and is
unreachable in production builds.

Hypothesis: with equation seeding normalized to the canonical width, (H1)
recovery-failure behavior at matched K becomes identical across all payload
widths, so per-K fixups and seed-attempt selections transfer 1:1 from the
cheap bb=2 stratum to every other width; and (H2) the normalization has no
measurable full-payload throughput cost, because it changes only which 64-bit
seeds are hashed, never the data plane.

## Arms

- Arm A (baseline, inherited policy): production seeding.  Seeds inherit the
  true-bb `SelectPeelPolicy(K, bb)` and the direct BlockBytes hash.  No
  experiment flags.
- Arm B (payload-independent profile): identical invocation plus
  `--payload-independent-seeding`.  Receipts must echo
  `payload_independent_seeding=1 payload_independent_seed_bb=2`.

Sanity gate: at bb=2 the two arms are bitwise-identical by construction
(the unit test `wirehair_v2_payload_independent_seeding_test` pins this).
Any bb=2 divergence between arms in the campaign output invalidates the run.

## Strata (matched across both arms)

- Payload widths: bb IN {2, 64, 256, 1280, 4096} -- all five, including the
  bb=2 identity stratum.
- Block counts: K IN {1000, 3200, 10000, 32000, 64000} (the precodefail
  default N list, covering every BlockCountBand).
- Overheads: {0, 1}.
- Completion: certified (mix 3) primary; mixed mix-2 secondary (all five
  widths are even, so mixed accepts them).
- Schedules: iid primary; adversarial and burst as secondary robustness
  strata.
- Seeds: one shared `--seed 0x5eedf411` for the primary run plus two
  replication seeds (0x5eedf412, 0x5eedf413).  The precodefail loss stream
  for a trial depends only on (seed, K, bb, overhead, trial), never on the
  arm, so A/B pairs at a cell share identical loss realizations and support
  paired (McNemar-style) analysis.

## Recovery-count arm (required)

Per (seed, schedule, completion) combination, run both arms:

    wirehair_v2_bench precodefail \
      --N 1000,3200,10000,32000,64000 \
      --bb-list 2,64,256,1280,4096 \
      --overhead 0,1 --trials 20000 --threads 16 \
      --loss 0.10 --seed 0x5eedf411 --schedule iid
      [Arm B adds: --payload-independent-seeding]

Collected metrics (existing CSV columns): `fail_rate`, `rank_fail`,
`first_rank_fail`, `failure_trials`, `inact_mu/max`, `binary_def_mu/max`,
`heavy_gain_mu/min`, `seed_attempt`.

Analyses:

1. Cross-width invariance (H1): within Arm B, for each K the per-width
   failure sets over the SAME loss stream seeds must coincide only where the
   loss streams coincide; the operative check is that `seed_attempt` and the
   K-indexed failure *rates* match across widths within binomial noise, and
   that any K flagged as a hotspot at bb=2 is flagged at every width.
2. Paired A-vs-B deltas per (K, bb) cell with McNemar tests (the
   `precodefail_paired` machinery already exists for mix pairing; A/B pairing
   is done offline from `failure_trials` lists, which record failing trial
   indices).
3. Fixup portability replay: take the K-indexed fixups derived on bb=2 under
   Arm B and confirm the same fixups reproduce the 88-91% reduction at
   bb={64,256,1280,4096}, replacing the prior ad-hoc `--seed-block-bytes 2`
   runs.

## Full-payload timing arm (required)

Both arms, all five widths including bb=2, matched trials and seeds.

1. Solver-path timing with real payload buffers (per-phase receipts:
   `solve_ms_mu`, `build_ms_mu`, `peel_ms_mu`, `project_ms_mu`,
   `residual_ms_mu`, `backsub_ms_mu`, `block_xors_mu`, `block_muladds_mu`):

       wirehair_v2_bench precodefail \
         --N 1000,3200,10000,32000,64000 \
         --bb-list 2,64,256,1280,4096 \
         --overhead 0,1 --trials 200 --threads 1 \
         --loss 0.10 --seed 0x5eedf411 --schedule iid \
         --full-payload-solve
         [Arm B adds: --payload-independent-seeding]

2. End-to-end codec throughput (encode/decode MB/s through the public facade
   paths, including precode arms):

       wirehair_v2_bench compare \
         --nlo 1000 --nhi 1000 --trials 25 \
         --bb-list 2,64,256,1280,4096 --max-message-mib 512 \
         --loss 0.10 --seed 0xc0decafe --precode
         [Arm B adds: --payload-independent-seeding]

       (repeat with --nlo/--nhi at 3200, 10000, 32000, 64000, capped by
       --max-message-mib exactly as compare already enforces)

   Timing runs are single-threaded and pinned per the existing thermal-run
   protocol (sole I2C thermal reader, ABBA interleave of Arm A / Arm B cells
   to cancel drift).

Acceptance for the timing arm: Arm B throughput within noise of Arm A at
every width (expected: identical data plane; only seed hashing differs), and
no per-phase regression beyond run-to-run variance.

## Receipts and bookkeeping

- Every Arm B output must carry
  `payload_independent_seeding=1 payload_independent_seed_bb=2` in its `#`
  header; every Arm A output must carry
  `payload_independent_seeding=0 payload_independent_seed_bb=0`.
- `--payload-independent-seeding` is mutually exclusive with
  `--seed-block-bytes` (the older ad-hoc normalization) and the CLI enforces
  this.
- Store raw CSVs under bench/results/ with the arm, seed, and schedule in the
  filename, as previous campaigns do.

## What must be built and run later

1. Build with tests + benchmarks (blocked now by the saturated campaign):
   `cmake -DBUILD_TESTS=ON -DWIREHAIR_BUILD_BENCHMARKS=ON` then build
   `wirehair_v2_bench` and the test targets.
2. Run `ctest` including the new `wirehair_v2_payload_independent_seeding_test`
   and the extended `wirehair_v2_bench_cli_test` receipts checks.
3. Execute the recovery-count and timing matrices above (both arms, all five
   widths, three seeds, iid + adversarial + burst).
4. Analyze per the three H1/H2 analyses; file the ship/no-ship decision on
   the bead.
