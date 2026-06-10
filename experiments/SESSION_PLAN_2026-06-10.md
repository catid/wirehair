# Wirehair Research Session Plan — Synthesis (2026-06-10)

Eight experiments (E1–E8) selected from 24 proposals. Coverage: 4 dense-replacement (E1, E3, E4, E5), 2 structure/scheduling (E2, E8), 2 implementation (E6, E7). Two merges: `dense-heavy-aware-total-rank` + `dense-staircase-ge-replay` → E1 (both are precode_sim instrumentation feeding one ranking script); `impl-lazy-window-dest-batch` → rider inside E7 (it is a 2-core-hour add-on once E7's schedules and E6's kernels exist).

**Dependency spine:** E1's ranking script (`rank_total.py`) scores E3/E4/E5. E3 Phase A's rule point supplies scheme tokens for E4/E5 (fallback if delayed: `S=round(0.8*Dprod)`, `D2=max(12, round(0.2*Dprod))`, Dprod = production `GetDenseCount(K)`). E6 Phase 1 gates E6 Phase 2 and the E7 rider. E7 Phase 2 gates the Phase 3 production refactor (next session).

**Machine discipline:** ≤120 threads across concurrent runs, always. Bandwidth-sensitive microbenchmarks (xor_bench, rowop_tiling replay, codec A/B) run only in **quiet windows** with no concurrent sweeps — co-running a 96-thread sweep pollutes DRAM bandwidth and invalidates GiB/s numbers.

---

## Wave schedule (≤120 threads per wave)

| Wave | When | Jobs (threads) | Wall |
|---|---|---|---|
| W0 quiet | 0:00 | Builds + self-tests; E1a muladd calibration (solo); E6 Phase 1 fanin microbench (solo, sequential) | ~0.7 h |
| W1 | 0:45 | E2 peel_sweep grid (64) + E1b precode_sim H-grid w/ --ge-replay (56) | ~3 h |
| W2 | ~4:00 | E3 Phase A grid (104) + E8 deck screening (8) + E7 Phase 1 WH_OPLOG capture (8) | ~7 h |
| W3 quiet | ~11:00 | E7 Phase 2 tiled replay + dest-batch rider (single core, perf counters); E3 fit script; E1 rank_total.py on all CSVs | ~2 h |
| W4 | ~13:00 | E3 Phase B certification (120), then E5 shuffle2/N1 matrix (120) | ~3.5 h |
| W5 | ~16:30 | E4 rowdist×precode cross, 5 invocations (120) | ~3 h |
| W6 | ~19:30 | E5 certification (80) + E4 best-pair certification (40) | ~3 h |
| W7 quiet | ~22:30 | E6 Phase 2 codec A/B + Phase 3 parity (single core); E7 rider re-runs if needed | ~2 h |

Total ≈ 24 h wall. W4–W7 are gated; if upstream kills fire, the session ends early with negatives recorded.

---

## E1 — Cost-model backbone: heavy-aware total ranking + real GE replay census
*(merged: dense-heavy-aware-total-rank + dense-staircase-ge-replay; runs first because every dense verdict downstream consumes its outputs)*

**Code changes**
- `experiments/peeling/xor_bench.cpp` (+~80 lines): `gf256_muladd_mem` kernel, same 32 GiB-target 5-repeat-median protocol as the XOR kernel; link gf256.cpp via `experiments/peeling/build.sh` (+2 lines). Output `results/muladd_calibration.csv` (us/op + muladd:xor ratio per block size).
- `experiments/precode/precode_sim.cpp` (+~200 lines):
  - TrialResult (:231) + RunTrial (:921-926) + Aggregate (:1085) + CSV (:1357/:1426): `heavy_value_muladds = H*Inactivated + H*H`, `heavy_divs = H`.
  - `--ge-replay`: new `ReplayGe(...)` (~140 lines) after `AnalyzeResidual` (:762) — build residual_rows × R bit matrix from the existing projections (rowbits loop :846-859), columns in inactivation order, rows by ascending projection weight; count word-XORs, row-eliminations, fill-in, pivotless-column positions. New columns: `ge_real_bitops_mu, ge_real_rowops_mu, fill_in_mu, def_outside_w18_rate, def_band_w95, def_band_w99`. Add `--ge-replay-reverse` sensitivity flag. Per-trial assert: `R - pivots_found == def`. Self-test clause: replay rank == RankAccumulator on the existing K=8..48 scheme loop (~20 lines).
- New `experiments/precode/rank_total.py` (~150 lines): `total_block_ops(bb) = recv_gen + precode_gen + sparse + backsub + ge_block_xors(→ge_real_rowops where present) + ratio(bb)*heavy_value_muladds`; separate bit-op ledger at measured word-op rate; re-score all existing CSVs for H∈{6,8,12,16} from def_pdf (assert def is H-independent in-script); pessimistic 2× heavy-term variant; report the crossover muladd:xor ratio.

**Run commands**
```bash
# W0, quiet machine:
cd /home/catid/wirehair/experiments/peeling && ./build.sh
./xor_bench --sizes 1280,102400,1048576 --target-gib 32 --repeats 5 --seed 1   # ~10 min

# W1 (56 threads, alongside E2):
cd /home/catid/wirehair/experiments/precode && ./build.sh && ./precode_sim --self-test
./precode_sim --K 1000,3200,10000,32000,64000 --oh 0,1,2 --threads 56 --ge-replay \
  --schemes dense,dense_h12,ldpc,ldpc_h12,ldpcdense_s40_d12,ldpcdense_s40_d12_h8,ldpcdense_s40_d12_h12,ldpcdense_s40_d12_h16 \
  --trials 20000   # use per-K trials 20000/20000/20000/10000/4000 (split invocations per K if the harness applies one count globally)
```
(Scheme sizes: rule point per K, i.e. `s=round(0.8*Dprod)`, `d2=max(12,round(0.2*Dprod))` — Dprod 50/62/86/190/346.)

**Wall:** 10 min calibration + ~3 h at 56 threads (ge-replay adds ≤30%).

**Decision rules**
- **Promote H=12:** ldpcdense_h12 beats dense_h6 on `total_block_ops` at all three block sizes for K≥3200 by ≥10% with fail ≤ dense_h6's → H=12 hybrids become the default tokens in E3/E4/E5; wirehair-wuw+edj quantitatively justified. If muladd:xor at 1280B exceeds the crossover ratio → record the H=8 fallback ranking instead (still a deliverable).
- **GE claim:** measured `ge_real_bitops` ratio (staircase/dense) ≤ closed-form ratio +20% at K=32000/64000 → confirm the −32..−68% GE headline; if fill-in pushes it above dense → formally correct digest 3, and rank_total.py uses ge_real_* from here on.
- **Heavy band:** if `def_outside_w18_rate` > 5% of def≤6 trials at K≥10000 → production's 18-column heavy band is materially under-spec; publish W99 as the band-width recommendation for wirehair-wuw.

**On success:** rank_total.py becomes the single objective for E3/E4/E5; W99 spec filed to wirehair-wuw; corrected GE ratios fed into the E4 25%-combined bar.

---

## E2 — Fold-scale × cap fine grid at percent overhead (struct-fold-cap-pct-grid)

**Code changes:** `experiments/peeling/peel_sweep.cpp` only, ~40 lines: append ~34 `kStructures[]` entries (L563) — caps {384,512,768,1024,1536,2048} × folds {0.10,0.25,0.40,0.60} as `kStructureLTFoldScale` (existing kind, entries precedent L678-685), `lt_m2_cC_d2_003`/`_d3_003` for C∈{512,1024,2048} (`kStructureLTExtra34`), plus comparators `lt_m2_c512_fold, lt_m2_c1024_fold, lt_m2_c2560, lt_m2_c3200, lt_m2_c1024_fold25, lt_m2_c256_fold`. Optionally extend the fold-scale self-test clause (L5782) with one new (cap, fold) point.

**Run commands**
```bash
cd /home/catid/wirehair/experiments/peeling && ./build.sh && ./peel_sweep --self-test
# W1, 64 threads; seeds sequential:
for S in 1 2 3 4; do
  ./peel_sweep --N 320,3200,32000 --N-jitter 10 --overhead-pct 0,1,2,5 --trials 200 \
    --threads 64 --seed $S --matrix-seeds paired --methods 10,18,40,115 \
    --structures <the ~40-name csv>
done
# Confirmation: top-3 per (N,pct) + method 120 at N=32000, --trials 400
```

**Wall:** minutes per seed at 64 threads; ~1 h total including confirmation.

**Decision rules**
- **Regression gate (validity):** OH0 cells must reproduce the 106.11 / 108.45 / 123.26 ordering, else stop and debug.
- **Promote:** a new structure beats `lt_m2_c256_fold` + rqcc_lowref at N=32000 pct1 (42.20 mean) by ≥10% (≤38.0) in ≥3/4 seed families → winner enters precode_sim cross-check (next session, E4-style) and feeds wirehair-jwg seed regeneration + wirehair-vbk.
- **Secondary:** interior fold optimum (some fold* ∈ {0.10..0.60} beats both endpoints by ≥2% at fixed cap at N=3200/32000) → fold becomes a standing dimension in future sweeps. Exclude N=320 from high-cap conclusions (cap clamping).
- **Kill:** if no structure beats the stale frontier by ≥5% at pct1/pct2, record that the high-cap advantage washes out with overhead — closes the gap with evidence.

---

## E3 — (S, D2, H) hybrid scaling rule across all K (dense-sd2h-scaling-rule)

**Code changes**
- `precode_sim.cpp` (+~20 lines): `--paired` flag — second Rng for received-row generation whose seed omits `Mix64(HashString64(token))` (main loop :1387-1392; thread into GenerateSystem's received-row loop :499-514) so cross-scheme fail deltas are CRN-paired.
- New `experiments/precode/run_sd2h_grid.sh` driver + fit script (~120 lines Python): least-squares S(K), D2(K) on grid points matching dense fail within CI.

**Run commands**
```bash
# W2, 104 threads. Phase A — per K in {1000,2000,3200,6400,10000,16000,20000,32000,48000,64000}
# (Dprod: 50,62,62,74,86,114,134,190,370,346):
./precode_sim --K <K> --oh 0,1,2,5 --trials <20000 if K<=20000 else 10000> --threads 104 \
  --paired --seed 0x5eed0001 \
  --schemes dense,ldpc,<ldpcdense_s{round(f*Dprod)}_d{D2} for f in 0.6,0.7,0.8,0.9,1.0 x D2 in {12,round(0.2*Dprod),round(0.3*Dprod)}, dedup>,<rule point>_h12 \
  > results/hybrid_K<K>.csv

# W4, 120 threads. Phase B — fitted rule + dense control, DISJOINT seed:
./precode_sim --K <each K> --oh 0,1 --paired --seed 0x5eed0002 --threads 120 \
  --schemes dense,<fitted rule point> \
  --trials <200000 K<=20000 | 100000 K=32000 | 50000 K=48000,64000>
```

**Wall:** Phase A ~7 h at 104 threads; Phase B ~2.5 h at 120.

**Decision rules**
- **Promote (ship table):** fitted rule fail ≤ dense +0.3pp (Phase A) and ≤ dense +0.15pp (Phase B) at every K, oh∈{0,1}, AND ≥10% inact + ≥20% sparse-XOR reduction vs dense at K≥10000 → deliver S(K)/D2(K) table to wirehair-edj.
- **Fallback:** any K failing the constant-fraction rule → fit a K-dependent fraction table (like GetDenseCount) — still the deliverable. f=0.9/1.0 columns are the hedge against the razor-sharp cliff at 32k/64k.
- **Bonus check:** K=48000, where production GetDenseCount is non-monotone (370 vs 346 at 64000) — flag if the fitted rule beats the table there.

**On success:** rule point becomes the canonical hybrid token for E4/E5; table filed to wirehair-edj; codec prototype scoped next session.

---

## E4 — Rowdist × precode interaction (dense-rowdist-precode-cross)

**Code changes:** `precode_sim.cpp` ~50 lines: two new RowDists mirroring peel_sweep frontier winners — `rs_c001_d50_c128` (robust-soliton tau-spike, weight law per peel_sweep `robust_soliton_weight` L4584) and `lt_m1_c16`. Registration: enum :119, DegreeSampler ctor :126-171, `--rowdist` parse :1315-1330 (LtM1C64/LtM2C1024 already exist unswept). Add chi-square degree-histogram self-test vs the analytic law (mirror peel_sweep L5782) — guards against silent weight-law port mismatch. Add `--max-inact` abort guard (runaway protection for lt_m1_c16 at K=64000).

**Run commands**
```bash
# W5, 120 threads; one invocation per rowdist (flag is global):
for RD in wirehair lt_m1_c16 lt_m1_c64 lt_m2_c1024 rs_c001_d50_c128; do
  ./precode_sim --rowdist $RD --paired --threads 120 --max-inact <2x dense inact at K> \
    --schemes dense,ldpc,ldpcdense_<E3 rule>,ldpcdense_<E3 rule>_h12 \
    --K 1000,3200,10000,32000,64000 --oh 0,1,2 \
    --trials <20000 K<=10000 | 10000 K=32000 | 4000 K=64000>
done
python3 rank_total.py --bb 1280,102400,1048576 results/rowdist_*.csv
# W6, 40 threads: certify best non-production (rowdist, scheme) pair per K: --oh 0 --trials 100000
```

**Wall:** ~3 h grid + ~1 h certification.

**Decision rules**
- **Promote:** some (rowdist ≠ wirehair, ldpc-family scheme) cell has fail ≤ dense+wirehair baseline +0.3pp at oh0 AND oh1, with combined recv+sparse+backsub+GE block-op total ≥25% below baseline at K∈{10000, 32000} (use E1's ge_real where available) → certify at 100k trials, flag pair for codec prototyping.
- **Interaction finding:** hybrid-vs-dense relative gain shifts >5pp under retuned rowdists → record that wirehair-vbk and wirehair-edj cannot be planned independently (decision-grade either way).
- **Kill branch:** every retuned rowdist breaks the oh0 floor (fail >3%) even with `_h12` → peel-only cost frontier certified non-shippable; close that branch of wirehair-vbk with failure-rate evidence.
- Report per-term deltas — backsub grows with inact width and may cap the 25% bar at 1MiB; partial wins (sparse/GE-only) must stay visible.

---

## E5 — Production-implementable D2: Shuffle-2 structure + N1 staircase degree (dense-prod-implementable-d2)

**Code changes:** `precode_sim.cpp` ~155 lines, no production files: (1) token suffix `_s2` → `SchemeKind::LdpcDenseShuffle2` — D2 rows from a Shuffle-2 deck (first row ~half window from seeded deck shuffle, each next row = previous + 2 deterministic flips; mirror ShuffleDeck16 logic from WirehairTools.cpp / MultiplyDenseRows WirehairCodec.cpp:786 — **document the exact deck/flip rule**, it must carry into the codec unchanged); (2) `ldpcdense_n1<X>_s<S>_d<D2>` parameterizing source-column parity hits X∈{2,3,4} (3-hit loop :394-415). Registration: GenerateSystem :338, MakeScheme :959, self-test kinds list +5 lines.

**Run commands**
```bash
# W4 (after E3 Phase B), 120 threads — matrix:
./precode_sim --paired --threads 120 --K 1000,3200,10000,32000,64000 --oh 0,1,2 \
  --schemes <{iid,_s2} x n1{2,3} x E3-rule sizes, plus dense control> \
  --trials <20000/10000/4000 by K>
# W6, 80 threads — certification (shuffle2 best-N1 vs iid control):
./precode_sim --paired --oh 0 --threads 80 --schemes <iid rule>,<shuffle2 rule> \
  --trials <200000 K<=10000 | 100000 K=32000 | 50000 K=64000>
# If N1=2 fail rises: D2-compensation arm, D2 in {0.2,0.3,0.4}*Dprod; report net encoder+decoder XOR delta.
```

**Wall:** ~1 h matrix + ~2.5 h certification (+~1 h 10^6-trial escalation at K≤10000 if the leak is marginal 0.1–0.2pp).

**Decision rules**
- **(a) Green-light codec prototype:** shuffle2-D2 fail within +0.2pp of iid-D2 at every K at certification counts, no def_pdf tail divergence beyond def H+3. If it leaks → deliverable is the minimal extra D2 (or per-row reseeding rate) restoring parity — still a production spec.
- **(b) N1 verdict:** N1=2 cuts precode_gen ≥25% at K≥32000 with fail within +0.3pp of N1=3 at equal-or-compensated D2 → adopt N1=2 (flips the encoder-side regression: ~130k vs 161k gen XORs at K=64000); else certify N1=3 as necessary with cost quantified. Check oh=0 inact **tail**, not just means.

**On success:** full implementable hybrid spec (S, D2, N1, D2-generation rule, seeds) → wirehair-edj codec prototype in `codec/` next session, carrying the documented bit-exact rule.

---

## E6 — k-ary gather-XOR kernels (impl-multisrc-gather-xor)

**Code changes**
- Phase 1: `experiments/peeling/xor_bench.cpp` (+~80 lines, parser L277/usage L257): `--fanin k` mode, dst ^= s1^..^sk for k∈{2,4,8}, sources from a 4 GiB pool (force DRAM), vs chained-pairwise baseline.
- Phase 2 (gated): `gf256.h/gf256.cpp` `gf256_add4_mem`/`gf256_add8_mem` SSE/AVX2/AVX-512 (~150 lines); `WirehairCodec.cpp` Substitute (:2916) + Encode (:4350) source-pointer gathering (~80 lines) behind `WH_GATHER`. Verify generated asm vectorizes (restrict pitfall).

**Run commands**
```bash
# W0, QUIET machine, right after E1a:
./xor_bench --fanin 2,4,8 --sizes 1280,102400,1048576 --target-gib 32 --repeats 5 --seed 1
# W7, QUIET, single core, performance governor — paired A/B (WH_GATHER whx vs master whx):
for N in 1000 4000 12000; do for BB in 1280 102400 1048576; do <whx A/B, 200 paired trials>; done; done
# Parity:
codec/wirehair_v2_bench compare --nlo 2 --nhi 64000 --trials 1000 --bb-list 1280,102400,1048576 --loss 0.1
```

**Wall:** Phase 1 <1 h (quiet); Phase 2/3 ~2 h quiet + coding time during W2–W5.

**Decision rules**
- **Kill gate (W0!):** add8 logical throughput <1.15× chained-pairwise at 1MiB → kill Phase 2 entirely; also lowers the prior on the E7 dest-batch rider (note in log).
- **Pass:** add8 ≥1.4× at 1MiB and ≥1.25× at 100KiB → build Phase 2.
- **Promote:** decode throughput +10% at N=1000/100KiB and N=1000/1MiB vs master whx, 1280B regression bounded at −2%, zero parity delta + bit-identical recovery → land behind WH_GATHER; **mandatory composed gather+tiling cell with E7** before claiming additive gains (overlapping traffic component).

---

## E7 — Real-schedule capture + byte-stripe tiled replay (impl-tiled-real-schedule-replay, with dest-batch rider)

**Code changes**
- Phase 1: `WirehairCodec.cpp` `WH_OPLOG` (~120 lines, WH_COUNT precedent :66/:463): CSV tuples (stage, op_type∈{xor,addset,add2,muladd,div,memcpy}, dst_block, src_block(s), scalar) for the four value stages.
- Phase 2: `experiments/tiling/rowop_tiling.cpp` (+~250 lines): schedule-file loader replacing the synthetic Op generator (struct Op L42, replay_untiled L110, replay_tiled L123), gf256 muladd support, tile sizes {4,8,16,32,64,128}KiB (tiles array L247), per-schedule checksum tiled==untiled.
- Rider (merged impl-lazy-window-dest-batch, +~120 lines): dest-batched replay mode reordering backsub window-apply ops into per-destination groups B∈{2,4,8} using E6's add4/add8; composed cell with 16KiB stripes.

**Run commands**
```bash
# W2 (8 threads): capture 50 decode schedules per (N,bb) in {1000,12000,40000} x {102400,1048576}, 10% loss, via bench/whx WH_OPLOG build.
# W3, QUIET, single core:
./rowop_tiling --schedules <dir> --tiles 4096,8192,16384,32768,65536,131072 --repeats 5
perf stat -e LLC-load-misses,l2_rqsts.miss <100KiB cells, 16KiB vs 32KiB tiles>   # the cliff
./rowop_tiling --dest-batch 2,4,8 --composed-with-stripe 16384 <N=12000,40000 schedules>  # rider
```

**Wall:** capture ~1 h background; replay ~2 h quiet.

**Decision rules**
- **Promote to Phase 3:** median real-schedule speedup ≥1.5× at 1MiB and ≥1.15× at 100KiB (best tile, baseline = untiled replay of the same real schedule; synthetic 2.09×/1.285× is the ceiling) → schedule the production record-then-replay refactor (~400 lines behind WH_STRIPED, gated bb≥64KiB, windowed-backsub table-aliasing trap covered by checksums) — **next session**, target end-to-end decode ≥1.3× at N=12000/1MiB.
- **Kill 100KiB targeting** if real-schedule gain there <1.08×; keep 1MiB-only path if it clears 1.5×.
- **Rider kill:** dest-batch incremental gain over 16KiB striping <8% → striping subsumes it, drop; ≥1.15× incremental → add to the Phase 3 design.
- Cliff deliverable regardless: perf-counter explanation of the 16→32KiB drop → tile-size auto-selection rule.

---

## E8 — Deck-dealt column-degree regularization (struct-deck-regular-columns)

**Code changes:** `peel_sweep.cpp` ~150–250 lines: `kStructureDeckColumns` (enum L496); deck dealing in `generate_random_structure_rows` (L5221, balanced/binary_p50 precedent) — per-matrix column-id deck replicated to the edge budget, shuffled, dealt with skip-and-requeue; `deckloss10/30` variants (deal across 1.11×/1.43× virtual rows, keep configured count, ~20 lines); `p2c` soft variant = BalanceMode 3 in `balanced_row_columns` (L5002, sample 2d candidates keep lowest-load d, ~30 lines); ~10 `kStructures[]` entries; self-test clause asserting exact column-degree balance (L5782).

**Run commands**
```bash
# W2, 8 threads (background; slow is fine), seeds 1-4 sequential:
./peel_sweep --N 320,3200,32000 --N-jitter 10 --overhead-pct 0,1,2,5 --trials 400 \
  --threads 8 --seed $S --matrix-seeds paired --methods 10,18,40,115 \
  --structures lt_m2_c256_fold,lt_m2_c512_fold,lt_m2_c1024_fold,deck_lt_m2_c256_fold,deck_lt_m2_c512_fold,deck_lt_m2_c128,p2c_lt_m2_c256_fold,deckloss10_lt_m2_c256_fold,deckloss30_lt_m2_c256_fold
# Cost arm: --N 1000,6400,32000 --structures deck_lt_m2_c16,deck_lt_m1_c64,deck_rs_c001_d50_c128,<iid baselines>
# Confirmation: --methods 120 on top-3 deck structures at N=32000.
```

**Wall:** ~2 h per seed family at 8 threads, fully overlapped with W2.

**Decision rules**
- **Promote:** `deck_lt_m2_c256_fold` mean residual ≤108.45 (c512_fold quality) at N=32000 OH0 with row_xors ≤7.3 (edge parity: matrix_xors within 2% of iid baseline) in ≥3/4 seed families → c1024_fold quality at ~16% fewer edges; send winner to precode_sim cross-check next session (peel-only wins are never shipping claims).
- **Loss-robustness gate:** deckloss30 retains ≥50% of the deck-vs-iid improvement, else only the p2c soft variant survives as the shippable form.
- **Cost arm:** deck_rs_c001_d50_c128 cuts solve_total_xors_est ≥5% vs iid at N=32000 → feeds wirehair-vbk.

---

## Session-close deliverables
1. `rank_total.py` verdict table (H=12 vs H=8, crossover ratio) + W99 heavy-band spec → wirehair-wuw.
2. S(K)/D2(K) rule or table with certification CIs → wirehair-edj.
3. Rowdist×precode interaction matrix + promote/kill of the retune branch → wirehair-vbk.
4. Shuffle-2/N1 implementability verdict with the documented bit-exact deck/flip rule.
5. pct-overhead frontier refresh + fold curve → wirehair-jwg.
6. Tiling real-schedule numbers + 100KiB cliff explanation; go/no-go for the WH_STRIPED Phase 3 refactor.
7. Gather-kernel verdict (and composed cell if both E6/E7 pass).
8. Deck-regularization verdict with loss-robustness bound.

File new beads issues for every gated next step (E7 Phase 3, E5 codec prototype, E4-winner certification follow-ups, impl-solve-order revival); close/update wirehair-vbk/edj/wuw/jwg/n34 per outcomes; commit results CSVs and `git push` per CLAUDE.md session protocol.

---

## Rejected but noteworthy
- **impl-solve-order-recovery-layout (25.3):** strong MTU lever, but a pervasive ~30-site indexing change with encoder/GE consistency traps; queue next session reusing E7's WH_OPLOG captures for its Phase 1 microbench.
- **dense-alt-precode-sweep (24.3):** the Gray-row arm is substantively the same family as E5's Shuffle-2 (2-bit-flip chain); run the IRA/2stair/band arms only after E5's verdict to avoid double-spending on the same question.
- **sched-lazy-avalanche (23.7):** clean engineering discharging wirehair-n34, but av_top64's ~7 cols is ~1.5% of GE width; run only if E2/E8 produce winners that need av-class scheduling online.
- **impl-shuffle2-window-seam-mst (23.5):** cheap exact counter, but E3/E5 may structurally shrink the same MultiplyDenseValues stage — wait for the edj verdict.
- **sched-d2-incremental-uf (23.3):** payoff conditional on v2 adopting rqd2-class online; do its 10-minute Step-0 perf gate next session.
- **struct-sc-window-locality (23.3):** a gate whose value is entirely contingent on E7 Phase 3 landing; run after the tiling verdict.
- **struct-d2-chain-backbone (22.7):** rank-suspect degree-2 chains invisible to peel-only metric; revisit via precode_sim directly if E8 succeeds.
- **sched-merge-scorekind / sched-av-frontier-credit (22.7/22.3):** sub-1% quality-frontier movement; cheap but not worth core-hours this session.
- **struct-ripple-shaped-degree-law (22):** highest implementation risk (fragile numeric inversion) for a speculative 3-8 col gain; needs a recursion-vs-simulated ripple self-check before it's trustworthy.
- **sched-tail-portfolio (22):** prior says the gain accrues before the 10% switch; the informative negative isn't worth the 150-250 lines now.
- **impl-m4rm-batched-triangle (19.7):** GE bitops lever redundant with the −68% structural cut E3/E4 are chasing; revisit only if the precode replacement dies.
- **sched-endgame-exact (17.7):** noise-level effect; pure diagnostic.

---

## Full proposal ranking

- **dense-rowdist-precode-cross** (score 31.0): Peel degree-distribution x precode interaction: do hybrid-precode gains survive (and compound with) the retuned low-cap row distributions?
- **dense-sd2h-scaling-rule** (score 30.3): Fit and validate a closed-form (S, D2, H) hybrid-precode rule across K, including all untested K points and the missing K>=32000 hybrids
- **impl-tiled-real-schedule-replay** (score 29.7): Make wirehair-2sc concrete: capture real row-op schedules, byte-stripe tiled replay, then production stripe loop
- **dense-prod-implementable-d2** (score 28.3): Production-implementability of the hybrid: Shuffle-2-structured D2 rows, seeded generation, and staircase N1=2 encoder-cost reduction
- **struct-fold-cap-pct-grid** (score 27.7): Fold-scale x cap fine grid plus d2/d3 mass, run at percent overhead: close the two largest named frontier gaps with zero new structure code
- **dense-staircase-ge-replay** (score 27.3): Real GE elimination replay on the residual system: validate the staircase GE-cost claim under a production-like elimination order, and census deficient-column positions for the heavy-band gap
- **dense-heavy-aware-total-rank** (score 27.0): Heavy-cost-aware single-objective ranking: fold GF(256) muladd, heavy GE, and encoder parity generation into one calibrated block-op total per block-size regime
- **impl-multisrc-gather-xor** (score 26.3): k-ary gather-XOR kernels (gf256_add4/add8_mem) for Substitute, Encode, and value stages
- **impl-solve-order-recovery-layout** (score 25.3): Make wirehair-g4c concrete: solve-order physical layout of recovery_blocks for the MTU scatter stream
- **struct-deck-regular-columns** (score 24.7): Deck-dealt column-degree regularization: remove the low-coverage column tail at zero extra edge cost
- **dense-alt-precode-sweep** (score 24.3): One controlled sweep each for four structured alternatives: IRA-accumulator (N1-profile), double staircase, banded dense rows, and Gray-sequence dense rows
- **sched-lazy-avalanche** (score 23.7): Lazy-greedy cached avalanche scoring (incremental av_top64 at shipping cost)
- **impl-shuffle2-window-seam-mst** (score 23.5): Make wirehair-dmx concrete: offline op-count study of MST/seam-optimized derivation order for MultiplyDenseValues
- **struct-sc-window-locality** (score 23.3): Spatially-coupled / banded column windows: test whether column locality is residual-free, enabling cache-local codec layouts
- **sched-d2-incremental-uf** (score 23.3): Incremental degree-2 component maintenance for rqd2-class schedulers (speed at equal quality)
- **struct-d2-chain-backbone** (score 22.7): Permutation-chain degree-2 backbone: equal-row, row-budget-fair transfer of the LDPC-staircase win into the received rows
- **sched-merge-scorekind** (score 22.7): Component-bridge ScoreKind: pick the big-component column whose d3 shell bridges to other components
- **sched-av-frontier-credit** (score 22.3): Frontier-credit avalanche key: score cascades by the degree-2 mass and component merges they create
- **struct-ripple-shaped-degree-law** (score 22.0): Analytic decreasing-ripple degree family with a designed stall point: replace black-box cap/fold search with a 2-parameter trajectory family
- **sched-tail-portfolio** (score 22.0): Deterministic tail portfolio: snapshot at 10% todo, finish with 3 complementary schedules, keep the best
- **impl-lazy-window-dest-batch** (score 20.0): Deferred window-table application in BackSubstituteAboveDiagonal: batch B windows per destination pass
- **impl-m4rm-batched-triangle** (score 19.7): Make wirehair-akm concrete: Method-of-Four-Russians batched elimination for the binary GE Triangle at K>=10k, with vpopcntdq pivot scoring in the same harness
- **sched-endgame-exact** (score 17.7): Exact component-decomposed endgame solver (replaces depth-3 beam in hybrids)
