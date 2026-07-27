# Resume notes — WH2 peel degree distribution training

Branch `feat/wh2-opcount-cost-model`, pushed through `bbf34ab`.
Beads: `wirehair-tz5b` (stages 3–4) is blocked by `wirehair-mggg` (large-K regression).
Durable summary is in `bd remember` under `wh2-peel-training-state`.
Stage 1–2 report: https://claude.ai/code/artifact/406d25e1-6c08-4f02-b0be-1bd3948ef287

---

## Where things actually stand

`tools/peel_table_final.json` is the deliverable: **137 block counts, K=2…8192,
115 trained and 22 held at shipped, median +17.4% goodput, max +35.3%.** Every
entry was measured end to end against a codec with no overrides.

The gain is concentrated, and this matters more than the median:

| K band | median gain |
|---|---|
| 2–15 | +4.9% |
| 16–31 | +5.2% |
| 32–63 | +16.4% |
| 64–127 | **+20.4%** |
| 128–255 | +20.9% |
| 256–511 | +12.5% |
| 512–2048 | +4.2% |
| above 2048 | **negative — ships stock** |

If the deployment target is large block counts, this work buys close to nothing.
Worth confirming the target before investing in stages 3–4.

## The tools, in the order they run

| file | what it does |
|---|---|
| `tools/peel_funnel.py` | one K: LHS screen on the cost proxy → local refine → real-codec rank |
| `tools/peel_direct.py` | one K, **real codec only**. Used for every K ≤ 128 |
| `tools/peel_sweep.py` | walks the calibrated ladder, warm-starting each K from its neighbour |
| `tools/peel_params.py` | linear interpolation between anchors; `--check` does hold-out |
| `tools/peel_validate.py` | **authoritative.** trained vs unmodified codec, reverts losers |

Nothing enters a table without passing `peel_validate.py`. That is not ceremony —
see below.

---

## Blocker: the search regresses above K=2048

Measured, 2000 trials per arm, same payload:

```
K=1536   trained 2444.4   shipped 2417.4    +1.1%   keep
K=2048   trained 2520.8   shipped 2532.1    -0.4%   revert
K=4096   trained 1023.3   shipped 1222.9   -16.3%   revert
K=8192   trained  742.5   shipped  993.5   -25.3%   revert
K=16384  trained  429.4   shipped  945.2   -54.6%   revert
```

Monotone, and large. **Do not re-run the existing sweep at large K expecting a
different answer** — it costs 461 s per block count and the validator reverts it.

Cause: stages 1–2 rank on the counter-based cost proxy and only stage 3 sees the
real codec, so at large K the pool handed to stage 3 is already bad. Best-of-a-bad-pool.

### Two ways forward — this is the decision to make first

**(a) Accept shipped above K=2048.** Defensible: headroom up there was already
+1.1% at 1536 and −0.4% at 2048 before things went wrong, because the cache cliff
dominates (at K=4096 the intermediates are 16 MB and *both* arms halve). Cost: zero.
Unblocks stages 3–4 immediately.

**(b) Search large K directly on the real codec** — the `peel_direct.py` approach
that worked below K=128, with no proxy in the loop at all. Expect roughly 1–3
minutes per candidate at K=16384 with a realistic payload, so a 40-point screen is
already an hour per block count. Before spending that, run the cheap experiment
below to find out whether there is anything to win.

### Cheap experiment worth running before choosing (b)

Take the shipped law at K=4096 and sweep tilt alone over, say, {−100, −50, −20, 0,
+20, +50} at 400 trials, everything else at stock. Two outcomes:

- If nothing beats stock by more than noise, the answer is (a) and it is settled
  in about ten minutes.
- If something does, (b) is worth funding and you already know which direction the
  optimum lies in.

This is the same one-axis probe that settled the interpolation question at K=256,
and it is far cheaper than a search.

---

## Stage 3 — production seeds (not started)

The WH1 mechanism is confirmed and is exactly the "fallback table" pattern:

- base tables `kDenseSeeds[100]` (N/4 → seed, for N ≥ 2048) and `kPeelSeeds[2048]`
  (N % 2048 → seed)
- a sorted **per-N fixup table** consulted *first* by binary search:
  `kDenseSeedFixups` (164 entries, `WirehairDenseFixups.inc`) and
  `WirehairPeelFixups.inc` (769 entries)
- `GetDenseSeedPreFixup` / `GetPeelSeedPreFixup` preserve the originals
- `WirehairTools.cpp:1109` — "Correction table first (binary search)"

Weakness is judged on **overhead**, and the two tables are searched *jointly* per N
— the peel fixup header notes cases where a peel-seed change alone cannot fix the
defect.

WH2 already has the harness: **`wirehair_v2_bench seedtable`**, with
`--N`, `--bb-list`, `--peel-candidates`, `--trials`. It emits
`N,bb,solver,structure,bucket,base_peel,tuned_peel,dense,…`. Verified working:

```
build-fast/codec/wirehair_v2_bench seedtable --N 320 --bb-list 4096 \
  --peel-candidates 8 --trials 3
```

**Stage 3 must run after the parameter table is frozen** — seeds are searched
against a fixed parameter set, so any later parameter change invalidates them.

Design sketch: for each N, generate the peel PMF from `peel_params.py`, sweep
`--peel-candidates`, keep N where the best candidate's overhead is materially
better than the base seed's, and emit only those as fixups. Sort ascending, binary
search at runtime — mirror the WH1 file format so the loader is unchanged.

## Stage 4 — WH2 vs WH1 mainline (not started)

10% of K values chosen at random, speed and overhead, graphs, report.
`tools/peel_validate.py` already emits exactly this shape (two arms measured back
to back at identical trials and payload) — extend it to take WH1 as the baseline
arm instead of shipped-WH2, rather than writing something new.

---

## Traps that cost real time this session — do not re-learn these

**`pgrep -f` lies.** It matches your own shell wrapper because the command line
contains the repo path or script name. I reported a job as running for ~40 minutes
that had never been created. Verify with `ps -eo pid=,args= | grep '[b]igk.py'`
(bracket trick) **and** confirm the output file exists.

**Never merge from a file that is still being written.** I built a table from the
validator's JSON mid-run and silently dropped K=8192. Check the process has exited
first.

**The proxy is not trustworthy.** Six of 11 small-K anchors and *every* large-K
anchor were regressions the proxy scored as wins. Nothing reaches a table without
`peel_validate.py`.

**An in-search incumbent guard cannot work here.** The shipped staircase degree
scale is a sentinel (`kStaircaseDegreeScaleUnset = -1.0`) that makes the codec
compute the scale internally, so any pool entry must name *some* scale and is a
hybrid — shipped law at a searched scale, which is neither arm. That hybrid read
1145.9 MB/s at K=4096 where true shipped reads 1239.2. The only honest baseline is
a run with **no overrides at all**.

**Solve rate does not need a payload.** At bb=64 a passing candidate costs 0.25 s
against 2.44 s at bb=4096, identical verdict; a non-decoder fails 10 of 10. Gate
cheap, rank expensive, survivors only. (`--bb-list 0` is *rejected* — the "bb 0"
that skips work is `--solve-block-bytes 0`, and only on the proxy harnesses.)

**Coordinate boxes must span what the codec builds.** `SmallBandStaircaseCount(2)=1`
and `GetDenseCount(16384)=114` rising to 390 near K=52224. A `[2,80]` box silently
loses both ends and reports only "FAILED to produce a winner."

**essearch has a calibrated cost model at only 26 block counts** — 2 3 4 6 8 12 16
24 32 48 64 96 128 192 256 384 512 768 1024 1536 2048 4096 8192 16384 32768 64000 —
and rejects every other K outright. An every-1024 grid produces nothing above 4096.

---

## Things I would do next, in order

1. **Run the cheap tilt probe at K=4096** (above). Ten minutes, settles the
   large-K question either way.
2. **Confirm the deployment K range** before funding anything else. The whole
   result lives in K=32–512.
3. Unblock `wirehair-mggg` with the answer from (1), then start stage 3.
4. **Scale the validator's trial count with K.** It uses a flat 2000 everywhere; at
   K=8192 with a 4096-byte payload that is ~67 GB per arm and 8 minutes for one
   measurement. Throughput converges fast and failures are ~0 up there.
5. **Re-check `p1`.** The search pins it at the box ceiling of 400 almost
   everywhere, but a direct sweep at K=128 showed the objective is *flat* above it
   (1549.7 at p1=400 vs 1558.2 at p1=3200). It is a near-free parameter and its
   tabulated values carry no meaning — consider fixing it rather than storing it.
6. **Revisit `scale`.** The funnel optimises it alongside the peel parameters and
   picks values well away from the shipped default. That is a second shipped
   default being changed, and it interacts with stage 3. Worth an ablation: how
   much of the +20% is peel shape and how much is scale?

## Open questions I did not resolve

- Is the K=32–255 peak a property of the distribution, or of where the completion
  band and the cache both happen to be comfortable? An ablation against
  `wirehair-za5v` (the fixed 12-row completion band at K≤100) would tell you.
- Small-K anchors are erratic (K=27 → tilt +349, K=28 → +668, K=29 → +1218) yet all
  three validate as wins. That is consistent with a very flat optimum, but it means
  those exact values are arbitrary within a wide band — do not read them as precise.
- Nothing here has been tested with gf16 rows enabled; the whole campaign ran
  `--mixed-gf16-rows 0` per the standing constraint.
