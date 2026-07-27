#!/usr/bin/env python3
"""Train peel parameters across the full K range, warm-starting from neighbours.

ORDER. The sweep starts at K=128 -- the best-characterised point, and the only
one whose answer has been validated end to end -- and walks OUTWARD in both
directions, handing each K the previous one's answer as a warm start. Half the
screen is then drawn from a window around that start and half from the full box,
so a wrong neighbour cannot trap the search.

BUDGETS SCALE WITH K, because a cell at K=64000 solves a 64000-column system.
Measured per-evaluation cost at 1,000 cells: 2.2 ms at K=128, 4.5 ms at K=512,
15.8 ms at K=2048, 79.5 ms at K=8192. A fixed 3,000-point screen would cost about
25 minutes at K=64000 on its own, and the 100,000-cell finals stage would cost far
more, so both shrink as K grows. That trades resolution for reach: at large K the
failure gate is coarser and the finals are noisier, and the recorded per-K entry
says which budget produced it.

THE OBJECTIVE IS KNOWN WRONG IN PLACES. Its failure term is a single-shot precode
solve that ignores the encoder's seed escalation and is non-monotone against real
codec overhead. At K=64 that costs the search most of the available gain -- it
returns tilt ~2 (real decode 1324 MB/s) where a direct sweep found tilt ~800
(~1500 MB/s). So a table built from this sweep is a table of "best under a proxy
known to understate steepening", not of optima. --check-real gates each K's winner
on the real codec so the table at least never contains a distribution that fails
to decode; it cannot make the proxy correct.
"""
import argparse
import json
import subprocess
import sys
import time

HERE = __import__("os").path.dirname(__import__("os").path.abspath(__file__))
sys.path.insert(0, HERE)
from peel_funnel import DENSE_COUNT, small_band_s          # noqa: E402


# The block counts essearch has a CALIBRATED COST MODEL for. Any other K is
# rejected outright ("essearch has no calibrated cost model for N=5120"), which
# is why an every-1024 grid silently lost 5120, 6144, 7168 and so on. This ladder
# is also the better training grid: it is geometric, so it spends points where
# the optimum moves fastest (tilt falls roughly 3x per doubling of K) and few
# where it is flat, and stage 2 interpolates the gaps.
CALIBRATED_K = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384,
                512, 768, 1024, 1536, 2048, 4096, 8192, 16384, 32768, 64000]


def k_list():
    return list(CALIBRATED_K)


def budget(k):
    """(screen, refine, finals, screen_cells, final_cells) sized for K."""
    if k <= 1024:
        return 3000, 400, 16, 1000, 100000
    if k <= 8192:
        return 1200, 250, 12, 1000, 20000
    if k <= 32768:
        return 500, 150, 8, 600, 6000
    return 250, 100, 6, 400, 2000


def gate_cfg(k):
    """(gate_trials, gate_bb, rank_top) -- the cheap solve-rate tier.

    A non-decoder fails 10 of 10, measured at K=512, so 25 trials is ample; and
    at bb=64 a passing candidate costs 0.25 s against 2.44 s at bb=4096, with an
    identical verdict. Only the top few survivors get a real-payload timing run,
    because THAT is the measurement that needs bytes moved.
    """
    return (25, 64, 3)


def real_trials_for(k):
    """Trials for the real-codec stage. This stage is now the whole cost of the
    sweep, and it scales with K, so it is budgeted rather than fixed: enough to
    reject a non-decoder and rank goodput to about a percent, not enough to
    resolve ties. Shipping candidates get re-measured properly."""
    if k <= 512:   return 800
    if k <= 4096:  return 400
    if k <= 16384: return 200
    return 100


def run_one(bench, k, init, seed, real_trials):
    s, r, f, sc, fc = budget(k)
    rt = real_trials or real_trials_for(k)
    cmd = [sys.executable, f"{HERE}/peel_funnel.py", "--K", str(k), "--bench", bench,
           "--screen", str(s), "--refine", str(r), "--finals", str(f),
           "--screen-cells", str(sc), "--final-cells", str(fc),
           "--seed", str(seed), "--show", "1"]
    if init:
        cmd += ["--init", ",".join(str(x) for x in init)]
    gt, gbb, rtop = gate_cfg(k)
    cmd += ["--real-trials", str(rt), "--gate-trials", str(gt),
            "--gate-bb", str(gbb), "--rank-top", str(rtop)]
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    el = time.time() - t0
    best, rejected = None, None
    for line in p.stdout.splitlines():
        t = line.split()
        if "rejected as non-decoding" in line:
            try: rejected = int(t[t.index("rejected") - 1])
            except (ValueError, IndexError): pass
        if len(t) > 6 and t[0] == "1" and "scale=" in line:
            try:
                kv = dict(x.split("=") for x in t if "=" in x)
                best = {"goodput": float(t[1]), "mbps": float(t[2]),
                        "oh_mean": float(t[3]),
                        "scale": float(kv["scale"]), "p1": int(kv["p1"]),
                        "tilt": int(kv["tilt"]), "dmax": int(kv["dmax"]),
                        "absorb": int(kv["absorb"])}
            except (ValueError, KeyError, IndexError):
                pass
    if best:
        best.update(K=k, S=small_band_s(k), seconds=round(el, 1),
                    screen=s, screen_cells=sc, finals=f, rejected=rejected,
                    real_trials=rt)
    return best, p.stdout, p.stderr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", default="build-fast/codec/wirehair_v2_bench")
    ap.add_argument("--out", default="tools/peel_table.json")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--real-trials", type=int, default=0,
                    help="0 = scale with K via real_trials_for()")
    ap.add_argument("--kmax", type=int, default=64000)
    ap.add_argument("--log", default=None)
    a = ap.parse_args()

    ks = [k for k in k_list() if k <= a.kmax]
    pivot = 128 if 128 in ks else ks[len(ks) // 2]
    up = [k for k in ks if k >= pivot]
    down = [k for k in ks if k < pivot][::-1]

    table, logf = {}, (open(a.log, "w") if a.log else None)

    def walk(seq, label):
        init = None
        for k in seq:
            best, out, err = run_one(a.bench, k, init, a.seed, a.real_trials)
            if best is None:
                print(f"  K={k:<6} FAILED to produce a winner", flush=True)
                if logf: logf.write(f"=== K={k} FAILED ===\n{out}\n{err}\n")
                continue
            table[k] = best
            # warm start the NEXT K from this one, on the offset lattice
            init = [int(round(best["scale"] * 100)), best["p1"],
                    best["tilt"] + 100, best["dmax"], best["absorb"]]
            print(f"  K={k:<6} tilt={best['tilt']:>+5} p1={best['p1']:>3} "
                  f"dmax={best['dmax']:>2} absorb={best['absorb']:>3} "
                  f"scale={best['scale']:>7.2f} {best['mbps']:>7.1f}MB/s "
                  f"OH={best['oh_mean']:.4f} ({best['seconds']}s)", flush=True)
            if logf: logf.write(f"=== K={k} ({label}) ===\n{out}\n")
            with open(a.out, "w") as fh:
                json.dump({str(kk): vv for kk, vv in sorted(table.items())}, fh, indent=1)
    print(f"  walking UP from K={pivot} ({len(up)} values)", flush=True)
    walk(up, "up")
    print(f"  walking DOWN from K={pivot} ({len(down)} values)", flush=True)
    walk(down, "down")
    print(f"\n  {len(table)} K values written to {a.out}")
    if logf: logf.close()


if __name__ == "__main__":
    main()
