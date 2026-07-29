# Resume notes — WH2 codec performance and reliability

Audit date: 2026-07-29
Branch: `feat/wh2-opcount-cost-model`
Primary Beads: `wirehair-rv4a`, `wirehair-rv4a.1`,
`wirehair-rv4a.2`, `wirehair-rv4a.3`, `wirehair-rv4a.4`

## 2026-07-29 current state and next work

The authenticated `wirehair-za5v` architecture campaign is closed with the
policy-required `no-survivor` decision recorded below; dispatch-v1 remains
unchanged.  Its source-pinned parallel verifier reproduces the complete
2,376-job result in 6 minutes 44 seconds rather than roughly 3 hours
21 minutes.  Stored-summary replay is bounded to 32 MiB compressed and
256 MiB decompressed by `c2d18fa`.

`wirehair-4azp` is also closed.  Commits `d029397` and `64d884e` add and fully
freeze the explicit `WIREHAIR_V2_PROFILE_TINY_MDS_2026_07` profile without
changing `WIREHAIR_V2_PROFILE_CURRENT`.  K=1 is exact repetition for every
uint32 packet ID.  K=2 maps IDs 0..256 to all points of
`P^1(GF(256))`, rejects larger IDs, and recovers from every two distinct
supported IDs with no weak seeds.  The complete public-API mapping and
GF(256)/0x14d arithmetic witness is 72,094 bytes with SHA-256
`03bced599bfe6bb916cc26d15987ee86ed680e129fed45c3b9b897f767eb67ca`.
Release, strict ASan+UBSan, package, exhaustive ordered-pair/tail, OOM/fault,
and all-K fingerprint gates pass; the final all-K run took 2,673.07 seconds.

The non-retroactive `wirehair-rv4a` seed-repair experiment has completed with
the authenticated no-survivor result below.  Its selector/contracts child
`wirehair-rv4a.1` is verified and closed by commits `7b43515` and `c120c92`.
The frozen repair-v1 policy SHA-256 is
`5e67a150d1f909d6ed80468185fa2dd0e82eb2fc3486c0fa662e213cf3100b42`.
The provisional pure8 ID/SHA is `19cccf775ce0bf09` /
`19cccf775ce0bf098c9a425cb349714c4c4a880e7cf136c3bc365e13c05089a5`;
pure9 is `a530f9105beaa450` /
`a530f9105beaa450dee70ad9b2a5cc54c944d3cd47f0aa6534630b8971608541`.
Healthy attempt-zero roots pay no probe; weak roots use the frozen
payload-independent zero-RHS first-full-rank selector over attempts 0..7, and
the decoder consumes exactly the serialized selected attempt without search.
No per-K tables or seed hand fixes are allowed.

Selector verification covers every K=2..100, all eight attempts, seven even
widths with full/partial tails and distinct payloads, exact semantic bridges,
OOM/error transactions and telemetry, endpoint state restoration, and
concurrent pure8/pure9/real dispatch solves.  Release, strict ASan+UBSan,
Clang ThreadSanitizer, production-symbol isolation, and the 575,991-cell
public fingerprint gate pass.  Repeated primary and independent source audits
ended clean.

`wirehair-rv4a.2` is verified and closed at pushed commit `7383e6b`.  The
additive v3 native/parser/campaign protocol leaves bandtiming v2 and public
dispatch unchanged.  Campaign v3 uses block_bytes=2, every K=2..100, all six
hard schedules, 768 training roots, and prebound joint sealed random-root plus
production-root lanes.  Representative wider even block sizes remain
prelaunch correctness gates rather than roster axes.

The retained canonical result-free plan is
`/home/catid/wh2-rv4a2-result-free-plan-568a9673.json`, with exact file
SHA-256
`762a3c01984cb57b42f8fb771ba13a14e81f9fedb17e5451c16d93743a094283`.
Its semantic plan, roster, policy, and seed-disjointness SHA-256 values are,
respectively,
`568a967393ab1b5a23f9cc518b2e87f91752972c0f5c94d4c00099f4d9087291`,
`c4d29e2cacfe8fdf762dc449e803a8d5bfcdd7817b948c59da61580f219d0d7a`,
`81d04c117d8fa414971057f8bd69a39d320139a872eac9c785f5ff8c73c02b20`,
and
`ca29960e58a231f2d0bb7212843ccdee84e7ff6e6c122e18fd29873e2bd8e7f2`.
Training contains exactly 1,188 jobs, 912,384 recovery cells, 152,064 unique
selectors, 7,299,072 attempt rows, and 110,398,464 native rows.  Its job/cell
hashes are
`70e0e4cc184fe478e53d8b877b35b434134d1d38a10b7d24c8eb52ecde259c0e`
and
`9dc8406a564df1f1d2bf5572e0580688a0dc0032453d60ecc31715e576e7d3cb`.
Each possible winner has one prebound sealed template of 1,188 jobs, 304,128
recovery cells, 25,443 unique selectors, 2,433,024 attempt rows, and
36,799,488 native rows.  Pure8's sealed job/cell hashes are
`dc7713a1ad4a4dcd4866c1617fc5b80b7e9b8011fc2e817b86f09f919ed31199`
and
`b54169a3a539087b96f58934d6ea5009e6d7a6d382bb8bb3c0677492d43f9d46`;
pure9's are
`4a9b72575f7a005cd839ab0c705c56a40b445e11e0122e46dbdfebb7f65ce2e9`
and
`c163f8f37a89dc9437304be3e4ea3ff757077049e8c55d45d04fe80e901e76b2`.

The external runtime receipt is
`/home/catid/wh2-rv4a2-preflight-pin-568a9673.json`, exact SHA-256
`22ec26d578fcba5bf71547ff5fe626cca8bc151561292426070e7310c00d7132`.
The source-forced campaign and independent source-forced verifier generated
byte-identical canonical receipts.  It binds benchmark SHA-256
`7b767372ba090e42e8120d58dc5fabad5f42e55811802a78edeae1deb5b9022`,
all 17 native/parser/runner/test/build inputs, and the live thermal source
identity.  Final frozen-source gates pass: 123/123 warning-strict v3 tests,
181/181 complete tool tests, focused Release and strict ASan+UBSan 5/5 each,
broad Release 43/43, packages 3/3, build-policy 1/1, TSan selector 3/3 plus
repairtiming CLI, and an exact 3/3 live receipt replay.  Primary and three
independent code-reading audits ended with a full clean pass.

`wirehair-rv4a.3` consumed those exact frozen inputs and completed all 1,188
training jobs.  The independently authenticated manifest, completion,
decision, and verification-receipt SHA-256 values are, respectively,
`c11d15c1c0c17e226ea45930cf11fb9d8cd9b43c05a4a3ff8d95c7b9f23a534a`,
`ab2a269488e64a5fba03716df9c284084cc3a71f6cdf69cb07b61d969456059c`,
`ee676d07ce82f71e9252fa2a26b031630470cebc650d9910b098fe83655a7e52`,
and
`32ff54529883e5c10d53cdd1340c22ff61dbc01e7fa366ba673595b7d16bffcd`.
The compressed summary and ledger SHA-256 values are
`21a23a5d940350dfd60e72276c82b929ce3738ce878d8adcd4fa510c3ed3b816`
and
`df659fbde0965ea38f03bc44f20ccebd32368240da456df4184d0943b56adae1`.
Independent hashing found zero mismatches across all 3,564 manifest-listed
job/result/log files, all worker logs are empty, and all 912,384 recovery
cells, 152,064 unique selectors, 7,299,072 attempt rows, 102,187,008 timing
rows, and 110,398,464 total native rows match the result-free freeze.

Repair was completely effective on the authenticated training population:

| Arm | Raw unrecovered @ +64 | Raw tail AUC | Repaired unrecovered @ +64 | Repaired tail AUC | Selected attempt > 0 |
|---|---:|---:|---:|---:|---:|
| `pure8_s0m1_d3` | 7,716 | 505,762 | 0 | 4,326 | 1,286 / 76,032 |
| `pure9_s0m1_d3` | 6,762 | 442,667 | 0 | 3,212 | 1,127 / 76,032 |
| dispatch-v1 control | 4,302 | 280,333 | — | — | — |
| WH1 control | 11,607 | 765,140 | — | — | — |

Pure8's selected-attempt histogram is `{0:74746, 1:1223, 2:58, 3:5}`;
pure9's is `{0:74905, 1:1066, 2:51, 3:9, 4:1}`.  Mean selector calls are
1.05164 and 1.04542; overall p99/max attempts are 2/5 and p99/max calls are
4/7.  Repaired recovery succeeds in all 456,192 cells per arm, with zero cap
exhaustion, final weakness, mismatch, uncommitted result, nonrepairable
regression, or OOM.

Neither arm passes the frozen WH1 throughput gates:

| Arm | Decoder feed / WH1 | Full decoder / WH1 | Full encoder / WH1 | Selected direct / dispatch |
|---|---:|---:|---:|---:|
| `pure8_s0m1_d3` | 2.13635 | 2.31527 | 1.78164 | 0.63822 |
| `pure9_s0m1_d3` | 2.22392 | 2.40291 | 1.84330 | 0.68231 |

All 594 jobs per arm have complete cross-arm and A/A timing panels.  Direct
solve passes decisively, but both arms independently fail decoder feed, full
decoder, and full encoder.  The authenticated decision is therefore
`status=no-survivor`, `candidate_order=[]`, `selected_survivor=null`, and
`sealed_authorization=forbidden`.  No sealed job was launched and dispatch-v1
remains the production profile.  The non-gating full-selector/forced-selected
encoder ratios are 1.04615 for pure8 and 1.04252 for pure9.

Thermal evidence covers all 1,188 native job intervals with 7,906 valid rows:
maximum CPU Tctl 62.5 C, maximum DIMM temperature 53.0 C, and EDAC CE/UE 0/0.
The runtime monitor saw mean occupancy 30.434/32 and all groups exited.

The frozen v3 summary has one reporting erratum, tracked by
`wirehair-rv4a.4`.  It blanket-counts the expected repair-triggering
attempt-zero real-payload `Wirehair_Error` as a fatal runtime error and omits
it from raw structural weakness.  The six schedules repeat each selector, so
the 5,412 pure8 and 4,422 pure9 physical Error rows correspond to 902 and 737
unique selectors; adding the 384 and 390 unique NeedMore selectors gives the
correct deduplicated retry-needed counts of 1,286 and 1,127.  Every one is
corroborated by the attempt-zero zero-RHS structural probe and later successful
repair, so genuine runtime errors are zero.  This additive erratum does not
change no-survivor: the three independent timing failures remain decisive.
The immutable v3 summary, completion, decision, and verification bytes must
not be rewritten.  Before any future outcome, fix the classifier, version the
post-v3 schema/policy, rerun mutation/replay tests, and establish a new
source/plan/pin freeze.

## 2026-07-29 authenticated all-K za5v result

The fresh witnessed-v2 `wirehair-za5v` selection completed and independently
replayed all 2,376 jobs and 1,824,768 matched raw cells.  The immutable
selection directory is
`/home/catid/wh2-za5v-selection-51fd529-v2`.  Its manifest SHA-256 is
`aa962f62b50b310987b62c4fc70abc7314f82ee9027ea1d6f180bbf38d918786`;
completion SHA-256 is
`692a0f58fc3e717b9fe6d25dd78ec9eb77c3fb29308612b9005b81fdd18fe704`;
ledger SHA-256 is
`fa63a95b233cbae26c30648f2f6554e173be622f31f55e4b547874ad7a79e8b1`;
and summary SHA-256 is
`250b6cfa3195637079a96f89a6fd4497eab9deaa4bd4f8280306ee0034f0c280`.
Every worker log is empty.  A source-pinned 32-process independent verifier
reproduced the exact completion hash, job count, and raw-cell count in
6 minutes 44 seconds after its serial authenticated finalizer completed.

The frozen decision is `no-survivor`; no holdout is permitted.  Dispatch-v1
remains the production profile.  Three candidates are on the descriptive raw
Pareto surface, but every candidate misses the recovery and cross-codec
throughput gates:

| Candidate | Weak constructions | Unrecovered @ +64 | Tail AUC | Direct / dispatch | Decoder / WH1 | Encoder / WH1 |
|---|---:|---:|---:|---:|---:|---:|
| `pure8_s0_d3` | 1,166 | 6,996 | 458,921 | 0.8161 | 1.4521 | 1.6314 |
| `pure8_s0m1_d3` | 1,208 | 7,248 | 475,309 | 0.8024 | 1.4294 | 1.6081 |
| `pure8_s0m1_d5` | 1,786 | 10,716 | 700,757 | 0.8526 | 1.5097 | 1.7277 |
| `pure9_s0m1_d3` | 1,034 | 6,204 | 406,561 | 0.8549 | 1.5116 | 1.6825 |
| dispatch-v1 control | — | 4,404 | 286,958 | 1.0000 | — | — |
| WH1 control | — | 12,491 | 822,794 | — | 1.0000 | 1.0000 |

Thus the best direct-solve arm is `pure8_s0m1_d3` at about 19.8% faster than
dispatch, while the most reliable arm is `pure9_s0m1_d3`: it halves the WH1
unrecovered count and tail area, but is still about 41% worse than dispatch on
both recovery measures and about 51% slower than WH1 end-to-end decode.  All
candidate-only recovery regressions are constructor-weak/seed-repairable:
there are zero nonrepairable regressions and zero constructor, encoder,
decoder, or direct runtime errors.  Preserve these raw counts.  A later
seed-repair experiment may retest the reliability-leading and speed-leading
architectures, but it is a new experiment and cannot retroactively promote
this selection.

Authenticated campaign thermal evidence covers 13,102/13,102 valid rows on
CPUs 0-31, with maximum CPU Tctl 66.875 C, maximum DIMM temperature 57.0 C,
and zero corrected or uncorrected EDAC errors.

## 2026-07-29 direct-solve evidence witness repair

The second protected `wirehair-za5v` selection attempt was stopped after
1,472/2,376 result files, before any summary, completion marker, decision, or
candidate outcome was read.  Source review found that direct timing proved
repeatability only against its own preflight solve; it did not independently
bind the solved intermediate to the encoder intermediate.  Decoder byte
recovery does not close that gap for every pair-local direct prefix.  Because
direct cost participates in architecture ranking, neither protected attempt
nor the earlier v1 screens can certify promotion.

The replacement receipt protocol is
`wirehair-v2-bench:bandtiming:dispatch-v1:v2`; its native row schema is
`wirehair.wh2.bandtiming.v2`.  Successful WH2 encoder, decoder, and direct
rows carry an `intermediate_sha256` witness; WH1 and unsuccessful rows carry
the exact `not_applicable` sentinel.  The semantic bridge separately carries
both direct intermediate hashes.  Native code fails closed when encoder and
scheduled payload intermediates differ, when direct preflight differs from
its encoder payload, or when any timed direct solve differs from both.  The
independent Python consumer selects that exact
protocol/schema/semantic/header format and requires one common witness per
WH2 arm and replicate across panels and scopes.

Historical v1 receipts retain their exact 52-column and semantic layouts and
remain replayable, but `valid_for_promotion` is false.  The all-K campaign is
now schema v2, binds the v2 bandtiming protocol in its result-free plan, and
rejects a legacy receipt both immediately after a worker and during final
result replay.  All 600 archived v1 screen receipts replay exactly as
non-promotable, while a live 360-row native v2 probe and replay pass with
authenticated thermal context.

Current prelaunch gates pass: 73/73 warning-strict parser/campaign tests,
39/39 broad Release CTests excluding the separately frozen long fingerprint
and `ci-once` tests, the focused Release and Clang ASan+UBSan native CLI
gates, and production `build_policy_e2e` in 132.38 seconds.  Independent
current-diff native and integrated Python source reviews are clean.  The
result-free selection plan SHA-256 is
`3a0413150506aa777c995774e80e7ca89e2dc34f97801d8089c8e085697807c0`.
Its logical roster, job-list, cell-set, policy, and holdout hashes remain
unchanged.

Current promotion artifacts are:

- benchmark binary:
  `84c6d4fd1ab410ffee55db68ab4decd65fd498ee37e8d1de3d1d0362697e5c4d`;
- native benchmark source:
  `db9e1b11bf809b958a9b9869fc55d6a62addb067df0779c5401eac5fd36138b0`;
- native CLI test:
  `e41501899e6b82c2b9d45262ce439f77e72401f4dda6e05e361e18303e799775`;
- strict parser:
  `ceeebe3cab0507fdfe86db9067d2e1aea1741aaeb1629cdcaafb589bb568e519`;
- parser tests:
  `c9921d227d4ce535ae27cbcd696ca8d80fe35590c05a03e5f31e7d0318f7f02b`;
- paired-context tool:
  `07cd2dd509e0847caaaab813ba2d7fe6f4cf57a4fe2d78e2c9400c2545f87c82`;
- campaign runner:
  `bae4336508d9d5fd9dd514950da709cb01cc129dcf9e0b36f3e8fe96e17346ef`;
- campaign tests:
  `1a5d2378c11bf2f69998dedf43204f0493fdf47f8082c61feae9d266ec25ec95`.

Historical next step (completed): push this witnessed-v2 checkpoint and run a
fresh 2,376-job selection without inspecting outcomes before authenticated
completion.  The completed result is recorded above.

## 2026-07-29 protected selection replay fix

The first protected `wirehair-za5v` selection attempt is invalid and produced
no selection decision.  All 2,376 native jobs finished, but the final
authenticated replay failed closed before publishing the summary, completion
record, or completion checksum.  The failed directory must not be used for
promotion.

The failure exposed a campaign-runner bug rather than a native-codec result.
The direct candidate/dispatch panel intentionally solves both arms at the
pair-local maximum first-success prefix when both recover, otherwise K+64.
Its recorded dispatch overhead can therefore differ when the candidate
changes, but the recovery aggregator incorrectly included that overhead in
the candidate-independent duplicate signature.  Result-blind forensics across
all 594 K/schedule groups found this expected direct-overhead drift in every
group and zero drift in the other shared controls.  No ranking, survivor, or
aggregate candidate performance or recovery result was derived or inspected.

The corrected runner excludes only the pair-local direct overhead.  It retains
the direct dispatch result class as an invariant, because a dispatch that
already recovered must remain solvable at a larger consistent prefix, while an
unrecovered dispatch uses the same K+64 prefix in every copy.  The native
parser independently enforces the same distinction.  Regression coverage now
accepts a realistic candidate-dependent prefix change and rejects drift in
every retained signature field.

Corrected protected artifacts:

- campaign runner:
  `72ee3a1f1576304fa66ebcd084045f54cff45f0ca32dd9f88e27144772739871`;
- campaign tests:
  `af096019e98d13fca634dcf5f887e9b2e36024376517dec8c4f8a221bc813d6a`;
- strict parser:
  `66a1c7e83149914b8063d03eafac771d9ff9edbaf3be084f30557d1938fd92c2`;
- parser tests:
  `27f3ca2761126630b985846948fefd1a7740417a6fc298ab8f1ca61ab92002ae`.

Warning-strict parser plus campaign tests pass 65/65; campaign-only tests pass
44/44.  Two independent source reviews and an exhaustive result-blind
2,376-job control audit found no analogous candidate-coupled field.  The
logical experiment remains frozen: roster, selection job list, matched cell
set, policy, and all four disjoint holdout hashes are byte-for-byte unchanged
from the prelaunch checkpoint.

Historical next step (superseded by the witnessed-v2 repair above): commit and
push the corrected v1 runner, then launch a fresh selection directory.  No v1
execution is promotion evidence.

## 2026-07-28 protected all-K campaign prelaunch checkpoint

`wirehair-za5v` remains open.  No protected all-K selection or holdout outcome
has been launched or read at this checkpoint.

The native test-only `bandtiming` path is implemented end to end.  Candidate,
dispatch-v1, and WH1 use shared nonzero payloads, matched raw construction and
loss seeds, all six loss schedules, counterbalanced candidate/control and A/A
panels, and authenticated CPU/cache/thermal context.  Encoder timing includes
initialization plus the first K symbols; cross-codec decoder timing covers the
first feed through first success; candidate-versus-dispatch direct timing
isolates the common solver.  Explicit candidate parameters are transactional
per encoder/decoder instance, and the untimed semantic bridge proves the
canonical explicit descriptor identical to real dispatch.

The immutable screens contain 220/220 frozen jobs and 380/380 local-neighbor
jobs.  Their raw archives hash to
`4d03f4af545113cc8c8c8e8bcaff73886cf025270a33b4ccb0b2d7de03a7e85f42`
and
`939252436007f0d1c6119a535a05e5558c0e68eba61fb6ffb2168bbb665bd063`.
Pure8 with S-1 and D2=3 led the local screen at exact geometric
candidate/dispatch ratios 0.754434441 encoder, 0.728543396 decoder, and
0.711185307 direct solve.  It had no candidate-only weak construction; the
single K=5 weakness was control-only and duplicated by cache state.  The
predeclared all-K roster keeps that leader plus the S+0, D2=5, and pure9
hedges, with no retries or seed repair before architecture selection.

The stable protected runner artifacts are:

- strict parser:
  `66a1c7e83149914b8063d03eafac771d9ff9edbaf3be084f30557d1938fd92c2`;
- parser tests:
  `27f3ca2761126630b985846948fefd1a7740417a6fc298ab8f1ca61ab92002ae`;
- campaign runner:
  `daa12df8587b0072d97c68aa86902fcc843fc05497887b9b56a0deb0fbe04899`;
- campaign tests:
  `92a24becb41918e20d2faf2d8bf754a3a2d133d49fb1626e624a8995a691c568`.

The exact result-free plan has roster SHA-256
`d7f032731d3f5b1a52d37cb27cc4236e92ef8be46efa5445934c8356288421fc`,
selection job-list SHA-256
`021dcb0607d83e5ce5da1d4ebbafbab43f1dfcef34209b64326a877fe97e2b63`,
matched-cell-set SHA-256
`6cc049811a0826e1be9fa750497ba03028fd263a2dc3328f9409b8cbaf9a8cdf`,
and policy SHA-256
`bd7f3bb0497e141af20e3c3189006702b77e832ddfd3581475a0303bd589780c`.
Selection is 2,376 jobs and 456,192 raw recovery cells per candidate across
every K=2..100 and every schedule.  Each possible holdout is 594 jobs with
disjoint construction/loss seed sets and 256 replicates per K/schedule.

All 600 archived screen receipts replayed exactly; all 1,800 legacy replicate
records enriched to the same candidate, dispatch, and WH1 construction
result/class, and all 600 stream hashes are unique.  Warning-strict Python
tests pass 64/64; the formerly racy real-signal probe also passes 50/50 author
stress and 20/20 independent root stress.  Release benchmark CLI and seed
selection tests pass, as do the prior broad Release 43/43 non-long tests,
strict Clang ASan+UBSan benchmark CLI gate, production build-policy gate, and
the exact all-K fingerprint run in 2,681.15 seconds.

Two independent exact-hash reviews certified the four frozen Python artifacts
without findings.  The broader review independently passed 120 signal tests,
200 exited-leader repetitions, three direct cleanup attacks, 199 malformed
schema cases, 4,624 paired recovery states/300,560 threshold outcomes, an
independent Student-t oracle, 200 randomized Pareto/ranking cases, every
result-free plan, and all 600 archived receipt replays.  Start and end hashes
were identical.

Repeated source-reading and adversarial passes repaired aggregate timing
weighting/coverage, raw holdout action routing, candidate-only failure
classification, malformed nested evidence, selection-manifest TOCTOU,
exited-leader descendant cleanup, asynchronous TERM/HUP cleanup races, and
caller-controlled unhashable schema fields.  A later independent gate exposed
and repaired the ready-file race in the signal regression itself.  The final
runner records signals and raises only at cleanup-safe points; a traced TERM at
the former exception-to-cleanup transition preserves the original exception
and leaves no live process group.

Five consecutive full-core fuzz waves completed 640/640 workers and
45,378,600,820 mutations with zero artifacts.  The hardened thermal sampler
remains live; under 100 percent CPU load the post-EXPO system is typically
about 62-64 C CPU and 50-52.25 C hottest DIMM, with zero sensor or EDAC
errors.

Historical next step (superseded by the witnessed-v2 repair above): this was
the prelaunch instruction for the v1 campaign.  Only the runner-derived sole
survivor from a fresh authenticated v2 selection may enter its predeclared
disjoint holdout.  Promote a new per-instance contract/profile and fingerprint
only if
the cumulative recovery, direct-solve, WH1 encoder/decoder throughput, and
reliability gates all pass; otherwise retain dispatch-v1 and record the exact
remaining gap.

## 2026-07-28 canonical tiny-completion checkpoint

The first `wirehair-za5v` checkpoint is complete on the exact current source
diff, SHA-256
`54365c4f6c734b97786bce21b022c0d5c48abc4189928bd3248f98cee978e5b1`.
The bead remains open: this checkpoint proves the existing canonical
dispatch-v1 completion path before the separate native band-timing experiment;
it does not select or promote a replacement band.

`PrecodeSolveStats` now has a test-only
`TinyMixedFastPathAcceptances` counter incremented only after
`TrySolveMixedCompletionQuotientTiny` accepts.  The forced on/off/automatic
test covers the exact canonical P244 frozen 10 GF(256) + 2 GF(2^16),
SmallBandD4, raw mix3 descriptor at widths two and four for every K=2..100.
It requires exactly one acceptance for every successful forced-on mixed solve
and zero for forced-off, nonmixed, and ineligible solves.  Explicit bb=10,
K=385, and separate-bucket seams prove the exclusion boundary.  A direct
deficient-syndrome oracle compares forced general and tiny paths for Success,
NeedMore, and Error under constant, independent, and grouped schedules,
including nonzero successful writeback and immutable failure buffers.
Deterministic solve-stat comparisons were expanded to cover every field in
both the decoder and benchmark helpers.

Exact final receipts were:

- Release and strict Clang ASan+UBSan default grids: 2,998 encoder and 6,342
  decoder cases; 3,970 tiny acceptances, all 3,970 successful forced-on mixed
  solves accepted, and every forced-off/nonmixed/ineligible count was zero.
- Release and two independent strict sanitizer `--full` runs: 47,139 encoder
  and 138,927 decoder cases, including 73,706 tiny acceptances, all 73,706
  successful forced-on mixed solves accepted, canonical cells=12,473, and
  K coverage=99/99.  All 120 Release shards emitted the same exact receipt.
- The direct solve oracle passed all mixed quotient factor/replay and
  constant/independent/grouped deficient-syndrome cases in Release and strict
  sanitizer builds.  Focused decode, full build, and production build-policy
  gates also passed.
- The final exact-tree Release CTest run passed 44/44, including the all-K
  fingerprint in 2,525.73 seconds and `build_policy_e2e` in 107.75 seconds.
  `LastTest.log` SHA-256 is
  `5790015edeca8b5d8e9fa35c0a5d9ea158b96cc7172453f68b0864ea69dfe2ac`.
- Three final full-core fuzz waves passed 125/125, 125/125, and 128/128
  workers with 3,052,662,567, 1,540,135,489, and 943,646,060 mutations,
  respectively, and zero artifacts.  Two dedicated replacement workers added
  6,679,964 stateful and 49,600,647 profile mutations, also with zero
  artifacts.  Together with the two earlier ten-minute waves, this checkpoint
  accumulated more than 11.7 billion clean fuzz mutations.

Three source-reading passes found and repaired proof gaps: three deterministic
stats were initially omitted, automatic-mode Success did not prove engagement,
and failure-buffer immutability was not checked.  A fresh adversarial pass on
the frozen final diff found no further bug.

The hardened thermal worker was stopped gracefully only after the complete
CTest and fuzz interval.  It captured 6,559 valid one-second rows over
6,558.000074 seconds with no gaps, invalid samples, read errors, or EDAC
corrected/uncorrected errors.  Whole-session maxima were 95.25 C CPU and
82.25 C DIMM.  Across 4,842 samples at or above 90 percent CPU busy, average
busy was 99.435801 percent and CPU/DIMM maxima were 91.25/82.25 C.  The
immutable CSV SHA-256 is
`0ba1b4d80635c3683f336292a0ac77b03511b3cb42db6a841447fd8d2243b301`;
the hardened sampler SHA-256 is
`2b84efa91375a96a4a64e09ce5bfd7cba0b85b75028f5a93470cd4ae58aadb01`.

Historical next step (completed and superseded): add a transactional,
test-build-only explicit `PrecodeParams` plus
`PacketRowConfig` encoder/decoder initialization seam, then build
`wirehair.wh2.bandtiming.v1`.  Candidate and dispatch must use that same
explicit path; an untimed semantic bridge must prove the canonical descriptor
identical to real dispatch.  The native receipt must retain all raw weak seeds,
use panel-local censoring, and measure encoder init plus first K symbols,
decoder first feed through first success, and candidate-vs-dispatch direct
solve as distinct estimands.  Do not close `wirehair-za5v` until the full
all-K, six-schedule speed/recovery acceptance contract is verified.

## 2026-07-28 verified same-path timing closure

`wirehair-lnfk` is closed after clause-by-clause acceptance verification on the
current tree.  The test-only
native `peeltiming` command now emits
`wirehair.wh2.peeltiming.v2`: candidate and identity controls both use explicit
thread-local PMF and staircase-scale hooks through the exact nonpublic
`dispatch-v1` target.  The no-hook path is untimed and must independently prove
identical equation rows, packets, intermediate solve bytes, counters, overhead,
full recovery, and payload.  Candidate/control and identity-A/A panels use
ABBABAAB slots with odd-replicate label swaps, one exact shared loss trace,
common recovery prefix, startup and cache preparation outside the timer, and
decoder-solve-only timing.

The native protocol supports all six bound schedules: IID, burst, permutation,
systematic-first, repair-only, and adversarial.  Every row receipts raw
construction/loss seeds, trace identity, full recovery, weak-seed result class,
phase and operation counters, CPU migration/fault/saturation state, and the
external cache/CPU/thermal context hash.  Python independently authenticates
the complete native stream and exact cardinality/order/balance, replays
semantics, aggregates four adjacent paired log ratios per replicate, forms the
95 percent Student-t interval across replicates, and applies an independent
identity-A/A timing floor.  Direct, Funnel, Sweep, and Validate use this paired
measurement for ranking or promotion; `compare` remains only the cheap recovery
screen.

Final bug-hunt repairs include:

- zero-RHS structural probing, so RHS-dependent solve failures are not
  mislabeled as weak equation seeds;
- ratio-consistent asymmetric log bounds for exact practical margins and
  half-width-only cross-calibration, preventing equal stock-panel bias from
  legitimizing itself;
- complete source/build provenance for `.in`, `.map`, `CMakeLists.txt`,
  tracked, untracked, and self-ignored `.gitignore` files, repository/global
  excludes, exact environment/configuration and include candidates,
  linked-worktree `commondir`, local/worktree config, split-index backing,
  index/HEAD/multi-hop symbolic refs, legacy `GIT_CONFIG`, relative HOME/XDG
  semantics, and lexical, resolved, symlink, and hardlink aliases;
- a sanitizer-specific 7200-second all-K fingerprint timeout while preserving
  the 3600-second release timeout; and
- exact production `BUILD_TESTS=OFF` rejection of test-only modes and usage.

Verified current-tree gates are: Release 44/44 CTests, including the exact
all-K fingerprint in 2397.25 seconds; the other 42/43 Clang
ASan+UBSan CTests plus the separately rerun monolithic fingerprint in 3823.74
seconds; five independently sharded sanitizer fingerprints with the exact
certified/mixed/mix2/two-anchor/dispatch digests below; production
`build_policy_e2e` in 145.20 seconds; and 113/113 warning-strict Python tests,
bytecode compilation, and `git diff --check`.  A final post-review focused
replay of seed selection plus both benchmark CLI suites passed 3/3 in both
Release and ASan+UBSan builds.  Independent final C++ and Python source-reading
passes are clean.  Through fuzz seed 2130, 1,758/1,758 sanitizer workers passed
9,119,320,266 mutations with zero artifacts.

The concurrent hardened thermal worker captured 21,130 complete one-second
rows over 21,129.000251 seconds.  Whole-session maxima were 85.125 C CPU and
76.000 C DIMM.  Across the 8,655 samples at or above 90 percent CPU busy,
average busy was 99.019950 percent, CPU peaked at 75.625 C, and DIMM peaked at
76.000 C.  Sensor read errors and EDAC corrected/uncorrected counts remained
zero.  The retained CSV SHA-256 is
`3b602510fc7b52d5b72360e4e64b8db6b4287704d48909cc2513401a2ece639a`.

Historical next step (completed and superseded): `wirehair-za5v` first added
the canonical mixed P244/D4 10+2 profile to the forced tiny-path grid and then
the separate v1 timing experiment.  The current promotion path is the
witnessed-v2 campaign described at the top of this file.

## 2026-07-28 architecture reconciliation

The four published V2 profiles remain byte-frozen on the
`LegacyD12/ProfileDerived` equation system.  They still use the legacy
`GetDenseCount(K)` staircase rule and up to 256 profile-derived seed attempts.
They were not rewritten to match an experiment.

Fresh peel work instead names the nonpublic `dispatch-v1` target explicitly:
contract ID `a98c37c23ee7feae`, `SmallBandD4`, mixed 10 GF(256) + 2 GF(2^16)
completion rows, P244 frozen geometry, mix3, and one table-free raw construction
seed.  `compare`, `precodefail`, `peelpmf`, and every Python training/validation
receipt must agree on that complete identity.  The target is not a public
serialized profile; production promotion remains `wirehair-2qlw`.

All earlier small-K/proxy evidence that did not bind this exact target is
invalid for architecture or parameter selection, independently of its reported
gain.  The target/profile plumbing (`wirehair-umrh`, `wirehair-uvkr`) and the
same-process, same-hook-path timing gate (`wirehair-lnfk`) are now complete.
Fresh table training must use those exact contracts and still requires
independent validation.

## 2026-07-28 verified target-plumbing closure

`wirehair-umrh` and `wirehair-uvkr` are closed. Their acceptance clauses were
mapped independently and replayed against the final tree before either Bead
was closed. Two final source-reading passes were clean.

The five equation contracts produced their exact frozen SHA-256 over every
K=2..64000:

- certified: `e6146e7fea89089689a819c72e7e82f799344b451f5e4653b125f622dee3de0b`
- mixed: `47bca161ff7b51684f39d19db3b5b0d11137f21335ec5c74d818d06756d93627`
- mixed mix2: `858321c2c0a07103b2615bb6586ce310105d37dd9b2120eb5b63527a6fcb5404`
- mixed mix2 two-anchor:
  `272b24d707befdeb08afdba83472ccdaeed9ebfac17b05bc6c251793b6b17a7b`
- dispatch-v1:
  `66902904804e0efb31971855da6d36dc74e06ed1653b28d6d993d271817a69b1`

Every all-K process also checked the grouped-GF(256) K=4096 seam against
`4c2af25b58d403d6df117a1bb872d25a27b4225ba2c927a57821511b8c42c7be`.
The remaining final gates passed: 43/43 audit CTests and 42/42 ASan+UBSan
CTests (the separately executed full all-K test excluded from those counts),
5/5 focused architecture/lifecycle seams, 71/71 warning-strict Python tests,
Python bytecode compilation, and `git diff --check`.

The raw target CLI is now singleton and fail-closed. An independent matrix
rejected all 30 duplicate/cross-alias cases and all 26 partial, omitted, or
multi-selector cases. Candidate/control gates require an actual successful
decode/solve, not merely exit zero and a receipt. Receipts serialize parsed
staircase scales canonically, so raw whitespace cannot split a line and
negative zero is exactly `0` in both native and Python paths.

During the final all-K run, 1,200 profile and 1,200 stateful fuzz workers all
passed, covering 2,980,603,266 mutations. Two additional captured load batches
brought the final-session total to 2,656/2,656 passing workers and
3,449,849,883 mutations. The hardened thermal sampler recorded 4,953
one-second samples over 4,952 seconds with a monotonic, complete 19-column
stream. Across the 602-sample all-K/load window CPU busy averaged 96.697%, CPU
peaked at 63.125 C, and the hottest DIMM peaked at 50.500 C; the complete
session CPU peak was 80.250 C. Sensor read errors and EDAC corrected and
uncorrected counts remained zero.

The next ready correctness gate is `wirehair-lnfk`: implement same-process,
same-hook-path paired peel timing before training any replacement table.

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

- `wirehair_v2_bench peelpmf --N K --target-profile dispatch-v1` exports the
  native PMF and structural
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
  build-audit/codec/wirehair_v2_bench peelpmf --N 128 \
    --target-profile dispatch-v1
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
   Use the closed `wirehair-lnfk` same-process paired protocol for peel
   rankings; use the separately versioned `wirehair.wh2.bandtiming.v2` scope
   specified by `wirehair-za5v` for small-K completion-band work.  Archived
   v1 receipts are replay-only and cannot support promotion.
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
