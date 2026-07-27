#!/usr/bin/env python3
"""Solve every small K individually, measuring only on the real codec.

For K <= 128 there is no reason to interpolate: each point is cheap to measure
directly, the anchors from the calibrated ladder are visibly noisy down there
(K=6 picks tilt -61 between neighbours at +112 and +53), and the optimum moves
faster with K than anywhere else. So this trains all of them, one at a time.

WHY THIS DOES NOT USE THE FUNNEL. The funnel's cheap stages rank on a predicted
cost that only exists at the 26 block counts essearch was calibrated for -- it
refuses everything else ("no calibrated cost model for N=5120"), which rules out
K = 5, 7, 9, 10, 11 and most of the range below 128. The real codec has no such
restriction, and at small K it is fast enough to search with directly: a 2,000
trial run at K=8 costs well under a second.

TWO TIERS, the same split the funnel uses and for the same reason:
  GATE  solve-rate only, so no payload and few trials -- a non-decoder fails
        10 of 10, and at bb=64 a passing candidate costs 0.25 s against 2.44 s
        at bb=4096, with an identical verdict.
  RANK  throughput, which DOES need a realistic payload or decode_MBps measures
        loop overhead instead of memory traffic. Only survivors pay for it.

Selection is goodput = decode_MBps * K/(K+OH) with failure a hard reject, not a
penalty term: a distribution that does not decode has no goodput at all.
"""
import argparse
import json
import os
import random
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from peel_funnel import Funnel                     # noqa: E402  (family/pmf math)

BOX = [("p1", 0, 400), ("tilt", -100, 1600), ("dmax", 2, 64), ("absorb", 0, 100)]


def stock_pmf(k):
    w = [1.0 / k] + [1.0 / (d * (d - 1)) for d in range(2, 64)]
    w.append(max(0.0, 1.0 - sum(w)))
    return w


class Direct:
    def __init__(self, a):
        self.a = a
        self.probes = 0

    def probe(self, k, v, trials, bb):
        """(fail, OH_mean, decode_MBps) or None."""
        w = Funnel.family(stock_pmf(k), v[0], v[1], v[2], v[3])
        if w is None:
            return None
        self.probes += 1
        env = dict(os.environ,
                   WIREHAIR_V2_PEEL_DEGREES=",".join(f"{x:.9f}" for x in w))
        p = subprocess.run(
            [self.a.bench, "compare", "--nlo", str(k), "--nhi", str(k),
             "--bb-list", str(bb), "--trials", str(trials), "--loss", "0.10",
             "--precode", "--precode-profile", "mixed", "--mixed-gf16-rows", "0"],
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

    def lhs(self, rng, n, centre=None, frac=0.25):
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
        lows = [max(lo, c - frac * (hi - lo)) for c, (_, lo, hi) in zip(centre, BOX)]
        highs = [min(hi, c + frac * (hi - lo)) for c, (_, lo, hi) in zip(centre, BOX)]
        return block(half, lows, highs) + block(n - half, glo, ghi)

    @staticmethod
    def clamp(v):
        return [max(lo, min(hi, int(round(x)))) for x, (_, lo, hi) in zip(v, BOX)]

    def solve(self, k, seed, centre):
        rng = random.Random(seed)
        t0 = time.time()
        pool = self.lhs(rng, self.a.screen, centre)
        # gate everything cheaply, keep the decoders
        alive = []
        for v in pool:
            r = self.probe(k, v, self.a.gate_trials, self.a.gate_bb)
            if r and r[0] == 0:
                alive.append((r[2] * k / (k + r[1]), v))
        if not alive:
            return None
        alive.sort(key=lambda x: -x[0])
        best_g, best_v = alive[0]

        # local refine, still on the cheap tier
        r_frac, used = 0.25, 0
        while r_frac > 0.02 and used < self.a.refine:
            moved = False
            for j, (_, lo, hi) in enumerate(BOX):
                step = max(1, int(round(r_frac * (hi - lo))))
                for s in (-step, step):
                    w = self.clamp([x + (s if jj == j else 0)
                                    for jj, x in enumerate(best_v)])
                    if w == best_v:
                        continue
                    used += 1
                    rr = self.probe(k, w, self.a.gate_trials, self.a.gate_bb)
                    if rr and rr[0] == 0:
                        g = rr[2] * k / (k + rr[1])
                        if g > best_g:
                            best_g, best_v, moved = g, w, True
            if not moved:
                r_frac *= 0.5

        # rank tier: real payload, more trials, on the top few
        # Same guarantee as the funnel: the shipped law is always ranked, so the
        # table can never contain a point worse than what already ships.
        incumbent = [100, 0, 64, 100]        # shipped law; scale is not searched here
        cands = [best_v] + [v for _, v in alive[:self.a.rank_top] if v != best_v]
        cands = cands[:max(1, self.a.rank_top - 1)] + [incumbent]
        out = []
        for v in cands:
            r = self.probe(k, v, self.a.rank_trials, self.a.rank_bb)
            if r and r[0] == 0:
                out.append((r[2] * k / (k + r[1]), r[1], r[2], v))
        if not out:
            return None
        out.sort(key=lambda x: -x[0])
        g, oh, mb, v = out[0]
        return {"K": k, "p1": v[0], "tilt": v[1], "dmax": v[2], "absorb": v[3],
                "goodput": round(g, 1), "oh_mean": oh, "mbps": mb,
                "seconds": round(time.time() - t0, 1), "probes": self.probes}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", default="build-fast/codec/wirehair_v2_bench")
    ap.add_argument("--kmin", type=int, default=2)
    ap.add_argument("--kmax", type=int, default=128)
    ap.add_argument("--out", default="tools/peel_table_small.json")
    ap.add_argument("--screen", type=int, default=40)
    ap.add_argument("--refine", type=int, default=48)
    ap.add_argument("--gate-trials", type=int, default=400)
    ap.add_argument("--gate-bb", type=int, default=64)
    ap.add_argument("--rank-trials", type=int, default=2000)
    ap.add_argument("--rank-bb", type=int, default=4096)
    ap.add_argument("--rank-top", type=int, default=3)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()

    table, centre = {}, None
    # walk upward so each K warm-starts from the one below it
    for k in range(a.kmin, a.kmax + 1):
        d = Direct(a)
        r = d.solve(k, a.seed, centre)
        if r is None:
            print(f"  K={k:<5} no candidate decoded", flush=True)
            continue
        table[k] = r
        centre = [r["p1"], r["tilt"], r["dmax"], r["absorb"]]
        print(f"  K={k:<5} p1={r['p1']:>3} tilt={r['tilt']:>+5} dmax={r['dmax']:>2} "
              f"absorb={r['absorb']:>3}  {r['mbps']:>7.1f}MB/s OH={r['oh_mean']:.4f} "
              f"({r['seconds']}s, {r['probes']} probes)", flush=True)
        with open(a.out, "w") as fh:
            json.dump({str(kk): vv for kk, vv in sorted(table.items())}, fh, indent=1)
    print(f"\n  {len(table)} K values written to {a.out}")


if __name__ == "__main__":
    main()
