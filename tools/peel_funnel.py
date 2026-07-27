#!/usr/bin/env python3
"""Train the peeling degree distribution for a given block count K.

REPLACES the evolution-strategy search (essearch --peel ...) for this problem.
The ES was measured against uniform random, Latin hypercube and grid search at
matched evaluation budgets: random search DOMINATED it for the first ~200
evaluations and the ES only edged ahead at 240 by 0.76%, a margin not resolved
at that sample size. The funnel below reaches the same objective value in about
seven seconds instead of tens of minutes, and is roughly a hundred lines.

WHY IT IS FAST -- the objective has two terms whose measurements need cell
counts three orders of magnitude apart, and every previous search paid the
expensive one's price for both:

  COST (pred_ns, from operation counters) is converged at 1,000 cells. Over 40
  random configs, ranking against a 500,000-cell reference gives Spearman
  +1.0000, 0.05% median error, and an identical top ten -- at 2 ms per candidate
  against 898 ms. The counters average a near-deterministic quantity, since the
  structure is fixed and only the matrix varies.

  FAILURE is a binomial estimate of a ~1% event, so it improves only as sqrt(n).
  Three disjoint 1,000-cell windows on one config gave 0.500 / 0.600 / 0.800%;
  at 100,000 they gave 0.793 / 0.768 / 0.784%. It cannot be read as a number at
  1,000 cells.

  But it CAN be read as a VERDICT there. Judging "failure <= 1.5%" at 1,000
  cells against 500,000-cell truth was 30/30 correct with no false passes,
  because the candidate population is bimodal: configs are either fine
  (0.70-1.38%) or badly broken (5.6-99.9%), and almost nothing lands near the
  threshold.

STAGES
  1  screen a Latin hypercube at 1,000 cells; rank on cost, gate on <= 1.5%
  2  refine locally, still at 1,000 cells (axis-aligned, halving radius)
  3  re-measure the finalists at 100,000 cells -- the only stage whose failure
     numbers may be quoted

STAGE 3 EXISTS BECAUSE STAGE 2 LIES. The minimum of thousands of noisy estimates
is biased low: a refined best of 4,684,435 at 1,000 cells re-measured at
4,806,796 at 100,000, a 2.6% winner's-curse gap. Treat 1,000-cell values as an
ordering and nothing else.

WHAT THIS DOES NOT FIX. The objective's failure term is a single-shot precode
solve that ignores the encoder's seed escalation, and it is NON-MONOTONE against
real codec overhead -- it reports 22.68% where the codec measures fail=0 with
OH_mean 0.0053 over 3,000 trials. So the funnel finds the best point under a
proxy known to be wrong in places. Finalists must clear a real-codec gate
(--check-real) before any of them means anything, and the ranking among
near-optimal finalists is not trustworthy until that term is replaced.

PER-K USE. Everything here is K-specific and none of it transfers:
  * the structure the codec builds  (S = GetDenseCount(K) for K > 100)
  * the shipped PMF                 (P(1) = 1/K)
  * the optimum itself              (tilt about 329 at K=128, 9 at K=1024;
                                     the K=128 winners FAIL EVERY TRIAL at K=512)
"""
import argparse
import random
import subprocess
import sys
import time

BENCH_DEFAULT = "build-fast/codec/wirehair_v2_bench"

# S = SmallBandStaircaseCount(K): floor(1.25*(root + K/root)/2) for K <= 100,
# else wirehair::GetDenseCount(K). D2 is hardcoded 12 in MakeCertifiedParams and
# H is ActiveMixedGF256Rows + ActiveMixedGF16Rows = 10 with gf16 rows 0. S and D2
# cannot be moved by any flag on the mixed/certified path, so they are pinned,
# not searched -- leaving them free makes the search WANDER OFF-STRUCTURE, where
# the cost proxy inverts (Spearman +0.10 at 6,9,1 against -0.80 at 30,10,12).
DENSE_COUNT = {
    2: 2, 4: 3, 8: 6, 16: 13, 32: 13, 64: 19, 128: 30, 192: 30, 256: 30, 384: 34, 512: 38,
    768: 42, 1024: 50, 1536: 58, 2048: 54, 3072: 62, 4096: 70, 5120: 70, 6144: 74, 7168: 74,
    8192: 78, 9216: 78, 10240: 86, 11264: 94, 12288: 98, 13312: 106, 14336: 114, 15360: 114,
    16384: 114, 17408: 118, 18432: 122, 19456: 130, 20480: 138, 21504: 138, 22528: 146,
    23552: 154, 24576: 158, 25600: 162, 26624: 170, 27648: 174, 28672: 178, 29696: 178,
    30720: 182, 31744: 190, 32768: 194, 33792: 202, 34816: 206, 35840: 214, 36864: 218,
    37888: 226, 38912: 234, 39936: 250, 40960: 266, 41984: 286, 43008: 302, 44032: 318,
    45056: 330, 46080: 346, 47104: 358, 48128: 370, 49152: 378, 50176: 378, 51200: 386,
    52224: 390, 53248: 390, 54272: 386, 55296: 378, 56320: 370, 57344: 366, 58368: 358,
    59392: 354, 60416: 346, 61440: 338, 62464: 338, 63488: 342, 64000: 346
}

BOX = [("scale", 100, 4000),   # centi; mean staircase row degree
       ("p1", 0, 400),         # hundredths of the incumbent's degree-1 mass
       ("tilt", 0, 1700),      # offset lattice; signed = value - 100
       ("dmax", 2, 64),
       ("absorb", 0, 100)]


def small_band_s(k):
    if k > 100:
        if k in DENSE_COUNT:
            return DENSE_COUNT[k]
        raise SystemExit(
            f"no GetDenseCount entry for K={k}; add it to DENSE_COUNT rather "
            f"than guessing -- pinning the wrong S puts the search off-structure")
    root = 1
    while (root + 1) * (root + 1) <= k:
        root += 1
    return (5 * (root + k // root)) // 8


class Funnel:
    def __init__(self, args):
        self.a = args
        self.k = args.K
        self.struct = f"{small_band_s(args.K)}:10:12"
        self.calls = 0

    def token(self, v):
        s, h, d2 = self.struct.split(":")
        return (f"{s},{h},{d2},{v[0]/100.0:.2f},0,20,20,20,100,100,100,100,100,16,1,"
                f"{v[1]},{v[2]-100},{v[3]},{v[4]}")

    def measure(self, vecs, cells):
        """[(objective, fail_percent)] -- columns read BY HEADER NAME.

        Positional indexing produced three separate wrong analyses during this
        work, including a table of failure rates in the tens of millions of
        percent, because the column order differs between subcommands.
        """
        out = []
        for i in range(0, len(vecs), self.a.batch):
            chunk = vecs[i:i + self.a.batch]
            self.calls += 1
            p = subprocess.run(
                [self.a.bench, "essearch", "--N", str(self.k),
                 # --s-lo 1 because SmallBandStaircaseCount(2) = 1 and the ES
                 # coordinate box otherwise floors S at 2, which rejects the
                 # structure the codec actually builds at the smallest K.
                 # --s-lo 1 / --s-hi 400 must SPAN the S the codec builds across the
                 # whole K range: SmallBandStaircaseCount(2) = 1 and
                 # GetDenseCount(16384) = 114, rising to 390 near K=52224. A
                 # default box of [2,80] silently rejects both ends -- it cost the
                 # K=2 anchor and the entire K >= 16384 tail of the first sweep.
                 "--fix-structure", self.struct, "--s-lo", "1", "--s-hi", "400",
                 "--target", str(self.a.target), "--lam", str(int(self.a.lam)),
                 "--measure", "--threads", str(self.a.threads),
                 "--cells", f"{self.a.cell_base}:{self.a.cell_base + cells}",
                 "--configs", ";".join(self.token(v) for v in chunk)],
                capture_output=True, text=True)
            hdr, got = None, []
            for line in p.stdout.splitlines():
                if "fail_rate" in line:
                    hdr = [c.strip() for c in line.split(",")]
                    continue
                if hdr and line[:1].isdigit() and line.count(",") > 10:
                    f = line.split(",")
                    try:
                        fail = float(f[hdr.index("fail_rate")])
                        pred = float(f[hdr.index("pred_ns")])
                    except (ValueError, IndexError):
                        got.append((None, None)); continue
                    pen = self.a.lam * (fail - self.a.target)
                    got.append((pred + (pen if fail > self.a.target else 0.0),
                                100.0 * fail))
            while len(got) < len(chunk):      # a rejected config emits no row
                got.append((None, None))
            out.extend(got[:len(chunk)])
        return out

    def lhs(self, rng, n, centre=None, frac=0.20):
        """Latin hypercube over the box, or over a local window around `centre`.

        With a warm start, HALF the screen is drawn from a window of +/- frac
        around the neighbouring K's answer and half from the full box. The local
        half exploits the fact that the optimum moves smoothly in K (tilt is
        straight in log K); the global half is insurance, because that smoothness
        is an empirical observation and not a guarantee -- and at K=64 the
        objective is known to be wrong, so a warm start there would inherit a
        wrong answer with no way to escape it.
        """
        def block(m, lows, highs):
            cols = []
            for lo, hi in zip(lows, highs):
                if hi <= lo:
                    cols.append([int(round(lo))] * m); continue
                cut = [lo + (hi - lo) * (i + rng.random()) / m for i in range(m)]
                rng.shuffle(cut)
                cols.append([int(round(c)) for c in cut])
            return [[cols[j][i] for j in range(len(BOX))] for i in range(m)]

        glo = [lo for _, lo, hi in BOX]
        ghi = [hi for _, lo, hi in BOX]
        if centre is None:
            return block(n, glo, ghi)
        half = n // 2
        lows, highs = [], []
        for c, (_, lo, hi) in zip(centre, BOX):
            w = frac * (hi - lo)
            lows.append(max(lo, c - w)); highs.append(min(hi, c + w))
        return block(half, lows, highs) + block(n - half, glo, ghi)

    @staticmethod
    def clamp(v):
        return [max(lo, min(hi, int(round(x)))) for x, (_, lo, hi) in zip(v, BOX)]

    def incumbent_fail(self):
        """Failure of the SHIPPED law at this K, at screen resolution.

        The feasibility gate cannot be an absolute percentage across the whole K
        range: at K=2 the shipped law itself fails most of the time (measured
        83.2% for one structure at K=2, 3.5% at K=32), so a flat 1.5% gate
        rejects every candidate including the incumbent and the search reports
        no winner. Gate RELATIVE to the incumbent instead -- a candidate has to
        be no worse than the law we are trying to beat, not better than a
        constant chosen at K=128.
        """
        neutral = [1400, 100, 100, 64, 100]      # scale is irrelevant to the probe
        (_, f), = self.measure([neutral], self.a.screen_cells)
        return f

    def run(self, seed):
        rng = random.Random(seed)
        t0 = time.time()

        # NO PROXY FEASIBILITY GATE. The proxy failure term is not a usable
        # gate at any K and is actively misleading at small K -- it reports 100%
        # at K=2, 24.4% at K=8 and 3.2% at K=32 where the real codec decodes
        # 0/2000 with OH_mean <= 0.0065, because it ignores seed escalation. A
        # relative gate merely rescales a meaningless number. So stages 1-2 rank
        # on COST ALONE, which IS sound (Spearman -1.000 against real throughput
        # on the tilt axis at the shipped structure), and reality decides in
        # stage 3. The risk this takes on is that cost alone prefers degenerate
        # laws, so the finals pool is deliberately wide and every member is put
        # to the real codec.
        self.gate = float("inf")

        pts = self.lhs(rng, self.a.screen, self.a.init)
        res = self.measure(pts, self.a.screen_cells)
        alive = [(o, f, p) for p, (o, f) in zip(pts, res)
                 if o is not None and f is not None and f <= self.gate]
        print(f"  screen  {self.a.screen} pts @{self.a.screen_cells} cells  "
              f"{time.time()-t0:5.1f}s  {len(alive)} passed <= {self.gate:.2f}%")
        if not alive:
            print("  nothing passed the feasibility gate")
            return None
        alive.sort(key=lambda r: r[0])
        best, bestv = alive[0][0], alive[0][2]

        t1, r, used = time.time(), 0.25, 0
        while r > 1e-3 and used < self.a.refine:
            cand = []
            for j, (_, lo, hi) in enumerate(BOX):
                step = max(1, int(round(r * (hi - lo))))
                for s in (-step, step):
                    w = self.clamp([x + (s if kk == j else 0)
                                    for kk, x in enumerate(bestv)])
                    if w != bestv:
                        cand.append(w)
            cand = cand[:self.a.refine - used]
            if not cand:
                break
            res = self.measure(cand, self.a.screen_cells)
            used += len(cand)
            moved = False
            for p, (o, f) in zip(cand, res):
                if o is not None and f is not None and f <= self.gate and o < best:
                    best, bestv, moved = o, p, True
            if not moved:
                r *= 0.5
        print(f"  refine  {used} pts @{self.a.screen_cells} cells  "
              f"{time.time()-t1:5.1f}s  {alive[0][0]:,.0f} -> {best:,.0f}")

        # THE INCUMBENT IS ALWAYS A FINALIST. Without it the search can return a
        # REGRESSION and not notice: at K=4096 the funnel picked 1034.9 MB/s
        # while the shipped law measures 1239.2, because LHS never sampled the
        # shipped point and nothing else in the pool beat it. Including it costs
        # one probe and makes "no worse than shipping" a property of the output
        # rather than a hope.
        # The incumbent is the SHIPPED point in full: shipped scale 14.00 AND the
        # shipped law (p1 100 = unchanged degree-1 mass, tilt 0 on the offset
        # lattice, dmax 64, absorb 100). Carrying the found scale instead makes
        # this a hybrid that is neither arm -- it measured 1145.9 at K=4096 where
        # the true shipped point measures 1239.2, so the "guarantee" guaranteed
        # nothing.
        incumbent = [1400, 100, 100, 64, 100]
        pool = [bestv] + [v for _, _, v in alive if v != bestv]
        pool = pool[:max(1, self.a.finals - 1)] + [incumbent]
        t2 = time.time()
        ok, dead = self.real_select(pool)
        print(f"  finals  {len(pool)} gated @bb{self.a.gate_bb}/"
              f"{self.a.gate_trials}tr, top {self.a.rank_top} timed @bb"
              f"{self.a.rank_bb}/{self.a.real_trials}tr  {time.time()-t2:5.1f}s  "
              f"{len(dead)} rejected as non-decoding")
        if not ok:
            print("  no candidate decoded -- the incumbent stands at this K")
            return None
        ranked = [(g, oh, mb, fail, v) for g, fail, oh, mb, v in ok]
        print(f"\n  K={self.k}  structure {self.struct}  "
              f"total {time.time()-t0:.1f}s, {self.calls} bench calls\n")
        print(f"  {'rank':>4} {'goodput':>10} {'MB/s':>9} {'OH_mean':>9}  parameters")
        for i, (g, oh, mb, fail, v) in enumerate(ranked[:self.a.show], 1):
            print(f"  {i:>4} {g:>10.1f} {mb:>9.1f} {oh:>9.4f}  scale={v[0]/100:.2f} "
                  f"p1={v[1]} tilt={v[2]-100} dmax={v[3]} absorb={v[4]}")
        # The 1,000-cell stages are an ordering; only this stage is a value.
        print(f"\n  NOTE goodput = MB/s * K/(K+OH); failure is a hard gate, not a "
              f"term.\n       Top finalists are typically within noise of each other.")
        return ranked

    def real_probe(self, v, trials, bb):
        """(fail, OH_mean, decode_MBps) from the REAL codec for one candidate.

        This is the only trustworthy recovery measurement. The proxy failure
        term is a single-shot precode solve at zero overhead and ignores the
        encoder's seed escalation (up to kMaxPacketSeedAttempts = 256), which
        absorbs essentially all of it at small K: the proxy reports 100% failure
        at K=2, 24.4% at K=8 and 3.2% at K=32, while the real codec decodes
        0/2000 at every one of those with OH_mean <= 0.0065. It is also
        NON-MONOTONE against real overhead at K=64, so no calibration rescues it.
        """
        import os
        w = self.family(self.stock_pmf(), v[1], v[2] - 100, v[3], v[4])
        if w is None:
            return None
        env = dict(os.environ,
                   WIREHAIR_V2_PEEL_DEGREES=",".join(f"{x:.9f}" for x in w))
        p = subprocess.run(
            [self.a.bench, "compare", "--nlo", str(self.k), "--nhi", str(self.k),
             "--bb-list", str(bb), "--trials", str(trials),
             "--loss", "0.10", "--precode", "--precode-profile", "mixed",
             "--mixed-gf16-rows", "0"],
            capture_output=True, text=True, env=env)
        hdr = None
        for line in p.stdout.splitlines():
            if line.startswith("codec"):
                hdr = line.split()
            elif line.startswith("v2_mixed") and hdr:
                g = line.split()
                try:
                    return (int(g[hdr.index("fail")]),
                            float(g[hdr.index("OH_mean")]),
                            float(g[hdr.index("decode_MBps")]))
                except (ValueError, IndexError):
                    return None
        return None

    def real_select(self, pool):
        """Two-tier: a cheap solve-rate GATE, then a real-payload RANK.

        The two jobs need completely different measurements, and conflating them
        was costing about 50x:

          GATE (does it decode?) is a solve-rate question, so it needs neither a
          real payload nor many trials. A non-decoder fails 10 of 10 -- measured
          at K=512 -- so 25 trials is ample, and dropping the payload from 4096
          to 64 bytes takes a passing candidate from 2.44 s to 0.25 s with an
          identical verdict and overhead. Moving bytes to answer "is the matrix
          solvable" is pure waste.

          RANK (how fast is it?) is a timing question and DOES need a realistic
          payload, or decode_MBps measures loop overhead instead of memory
          traffic. But only the survivors need it, not the whole pool.

        J = decode_MBps * K/(K+OH): overhead is a bandwidth term, since
        delivering K blocks costs K+OH packets sent, and failure is a HARD GATE
        rather than a large number -- a distribution that does not decode has no
        goodput at all. Both facts are measured: at K=512, seven of nineteen
        arms failed EVERY trial, and they were exactly the arms that win at
        K=128.
        """
        passed, dead = [], []
        for v in pool:
            r = self.real_probe(v, self.a.gate_trials, self.a.gate_bb)
            if r is None:
                continue
            fail, oh, mb = r
            (dead if fail > 0 else passed).append((fail, oh, mb, v))
        out = []
        for _, _, _, v in passed[:self.a.rank_top]:
            r = self.real_probe(v, self.a.real_trials, self.a.rank_bb)
            if r is None:
                continue
            fail, oh, mb = r
            if fail > 0:                      # slipped through the cheap gate
                dead.append((fail, oh, mb, v)); continue
            out.append((mb * self.k / (self.k + oh), fail, oh, mb, v))
        return sorted(out, key=lambda x: -x[0]), dead

    def stock_pmf(self):
        w = [1.0 / self.k] + [1.0 / (d * (d - 1)) for d in range(2, 64)]
        w.append(max(0.0, 1.0 - sum(w)))
        return w

    @staticmethod
    def family(stock, p1, tilt, dmax, absorb):
        dmax = max(2, min(int(dmax), len(stock)))
        w = [0.0] * dmax
        w[0] = (p1 / 100.0) * stock[0]
        for d in range(2, dmax):
            w[d - 1] = stock[d - 1] * (d ** (-tilt / 100.0))
        a = min(absorb, 100) / 100.0
        w[dmax - 1] = (a * sum(stock[dmax - 1:]) + (1 - a) * stock[dmax - 1]) \
            * (dmax ** (-tilt / 100.0))
        tot = sum(x for x in w if x > 0)
        return [max(x, 0.0) / tot for x in w] if tot > 0 else None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--K", type=int, required=True)
    ap.add_argument("--bench", default=BENCH_DEFAULT)
    ap.add_argument("--screen", type=int, default=3000)
    ap.add_argument("--refine", type=int, default=400)
    ap.add_argument("--finals", type=int, default=20)
    ap.add_argument("--screen-cells", type=int, default=1000)
    ap.add_argument("--final-cells", type=int, default=100000)
    ap.add_argument("--gate-mult", type=float, default=1.25,
                    help="gate is max(--gate, this x the incumbent's failure at "
                         "this K); a flat gate rejects everything at small K")
    ap.add_argument("--gate", type=float, default=1.5,
                    help="coarse feasibility threshold in percent at screen-cells")
    ap.add_argument("--target", type=float, default=0.0095)
    ap.add_argument("--lam", type=float, default=85e6)
    ap.add_argument("--threads", type=int, default=64)
    ap.add_argument("--batch", type=int, default=60)
    ap.add_argument("--cell-base", type=int, default=900_000_000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--show", type=int, default=10)
    ap.add_argument("--init", type=str, default=None,
                    help="warm start 'scale,p1,tilt,dmax,absorb' on the OFFSET "
                         "lattice (tilt = signed + 100), normally the answer "
                         "from a neighbouring K")
    ap.add_argument("--init-frac", type=float, default=0.20)
    ap.add_argument("--real-trials", type=int, default=2000,
                    help="trials for the RANK tier (real payload)")
    ap.add_argument("--gate-trials", type=int, default=25,
                    help="trials for the cheap solve-rate gate")
    ap.add_argument("--gate-bb", type=int, default=64,
                    help="payload for the gate; solve rate does not need bytes")
    ap.add_argument("--rank-bb", type=int, default=4096,
                    help="payload for throughput ranking; must be realistic")
    ap.add_argument("--rank-top", type=int, default=3,
                    help="how many gate survivors get a real timing measurement")
    a = ap.parse_args()
    if a.init:
        a.init = [int(x) for x in a.init.split(',')]
    f = Funnel(a)
    ranked = f.run(a.seed)
    return 0 if ranked else 1


if __name__ == "__main__":
    sys.exit(main())
