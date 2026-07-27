#!/usr/bin/env python3
"""Per-K peel parameters: measured anchors plus LINEAR interpolation between them.

Stage 2 of the training plan. The anchors come from tools/peel_table.json, which
tools/peel_sweep.py produces by running the funnel at every block count essearch
has a calibrated cost model for:

    2 3 4 6 8 12 16 24 32 48 64 96 128 192 256 384 512 768 1024 1536 2048
    4096 8192 16384 32768 64000

That ladder is not a design choice so much as a hard constraint -- essearch
rejects any other K outright ("no calibrated cost model for N=5120"), which is
why an every-1024 grid silently produced nothing above 4096. It happens to be a
good grid anyway: it is geometric, and the optimum moves fastest at small K
(tilt falls roughly 3x per doubling) and is nearly flat at large K.

INTERPOLATION IS LINEAR IN K, as specified. Note this differs from an earlier
log-space fit: ln(tilt) is close to straight in log2(K), so log-space
interpolation is more accurate in the middle of a wide gap. Linear is used here
because it is what the codec will implement and because the anchors are dense
enough at small K -- where tilt moves fastest -- that the difference is small.
Between the widest gaps (16384 -> 32768 -> 64000) tilt is nearly flat and close
to zero, so the choice barely matters there either.

WHICH COORDINATES ARE WORTH INTERPOLATING
  tilt   yes -- it carries the objective and moves smoothly and steeply with K
  scale  yes -- moves smoothly, real effect on throughput
  dmax   weakly -- moves smoothly but its effect is within noise at fixed tilt
  p1     NO REAL MEANING. The search pins it at the box ceiling of 400 for almost
         every K, and a direct sweep at K=128 shows the objective is FLAT above
         it (decode 1549.7 at p1=400 against 1558.2 at p1=3200, inside noise).
         It is interpolated for completeness, but do not read significance into
         the value.
  absorb likewise weak; it wanders across the anchors with no trend.
"""
import argparse
import bisect
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_TABLE = os.path.join(HERE, "peel_table.json")

FIELDS = ("scale", "p1", "tilt", "dmax", "absorb")


def load(path=DEFAULT_TABLE):
    with open(path) as fh:
        raw = json.load(fh)
    return {int(k): v for k, v in raw.items()}


def params_for_k(k, table):
    """(scale, p1, tilt, dmax, absorb) for any K, linear between anchors.

    Outside the anchor range this returns the nearest anchor UNCHANGED rather
    than extrapolating. Extrapolation here is not merely inaccurate, it is
    dangerous: the parameters that win at one K can fail EVERY trial at another
    (the K=128 winners fail 4000/4000 at K=512), so inventing a point beyond the
    measured range risks producing a distribution that does not decode.
    """
    ks = sorted(table)
    if not ks:
        raise SystemExit("empty anchor table")
    if k <= ks[0]:
        return {f: table[ks[0]][f] for f in FIELDS}, "clamped-low"
    if k >= ks[-1]:
        return {f: table[ks[-1]][f] for f in FIELDS}, "clamped-high"
    if k in table:
        return {f: table[k][f] for f in FIELDS}, "anchor"
    i = bisect.bisect_left(ks, k)
    lo, hi = ks[i - 1], ks[i]
    t = (k - lo) / (hi - lo)
    out = {}
    for f in FIELDS:
        a, b = table[lo][f], table[hi][f]
        v = a + t * (b - a)
        out[f] = v if f == "scale" else int(round(v))
    return out, f"interp {lo}..{hi}"


def pmf_for_k(k, table):
    """The peel PMF this K should use, in the shipped family."""
    p, _ = params_for_k(k, table)
    stock = [1.0 / k] + [1.0 / (d * (d - 1)) for d in range(2, 64)]
    stock.append(max(0.0, 1.0 - sum(stock)))
    dmax = max(2, min(int(p["dmax"]), len(stock)))
    w = [0.0] * dmax
    w[0] = (p["p1"] / 100.0) * stock[0]
    for d in range(2, dmax):
        w[d - 1] = stock[d - 1] * (d ** (-p["tilt"] / 100.0))
    a = min(p["absorb"], 100) / 100.0
    w[dmax - 1] = (a * sum(stock[dmax - 1:]) + (1 - a) * stock[dmax - 1]) \
        * (dmax ** (-p["tilt"] / 100.0))
    tot = sum(x for x in w if x > 0)
    return [max(x, 0.0) / tot for x in w] if tot > 0 else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", default=DEFAULT_TABLE)
    ap.add_argument("--K", type=int, help="print parameters for one K")
    ap.add_argument("--pmf", action="store_true", help="also print the PMF")
    ap.add_argument("--check", action="store_true",
                    help="hold out each interior anchor and interpolate it")
    a = ap.parse_args()
    table = load(a.table)

    if a.check:
        ks = sorted(table)
        print(f"  HOLD-OUT CHECK: drop each interior anchor, interpolate from "
              f"its neighbours\n")
        print(f"  {'K':>7} {'tilt true':>10} {'tilt interp':>12} {'err':>8} "
              f"{'scale true':>11} {'scale interp':>13}")
        for k in ks[1:-1]:
            held = {kk: v for kk, v in table.items() if kk != k}
            got, _ = params_for_k(k, held)
            tt, gt = table[k]["tilt"], got["tilt"]
            err = abs(gt - tt)
            print(f"  {k:>7} {tt:>+10} {gt:>+12} {err:>8} "
                  f"{table[k]['scale']:>11.2f} {got['scale']:>13.2f}")
        return

    if a.K:
        p, how = params_for_k(a.K, table)
        print(f"  K={a.K} ({how}): scale={p['scale']:.2f} p1={p['p1']} "
              f"tilt={p['tilt']} dmax={p['dmax']} absorb={p['absorb']}")
        if a.pmf:
            w = pmf_for_k(a.K, table)
            print("  " + ",".join(f"{x:.9f}" for x in w))
        return

    print(f"  {len(table)} anchors from {a.table}\n")
    print(f"  {'K':>7} {'scale':>8} {'p1':>5} {'tilt':>6} {'dmax':>5} "
          f"{'absorb':>7} {'MB/s':>9} {'OH':>8}")
    for k in sorted(table):
        v = table[k]
        print(f"  {k:>7} {v['scale']:>8.2f} {v['p1']:>5} {v['tilt']:>+6} "
              f"{v['dmax']:>5} {v['absorb']:>7} {v.get('mbps', 0):>9.1f} "
              f"{v.get('oh_mean', 0):>8.4f}")


if __name__ == "__main__":
    sys.exit(main())
