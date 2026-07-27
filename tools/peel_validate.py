#!/usr/bin/env python3
"""Verify every trained K against the TRUE shipped codec, and reject regressions.

WHY THIS EXISTS RATHER THAN AN IN-SEARCH GUARD. The obvious fix for "the search
returned something worse than shipping" is to put the incumbent in the candidate
pool. That cannot express the incumbent here: the shipped staircase degree scale
is a sentinel (kStaircaseDegreeScaleUnset = -1.0) that makes the codec compute
the scale internally, so any pool entry must name SOME scale and is therefore a
hybrid -- shipped law at a searched scale, which is neither arm. That hybrid
measured 1145.9 MB/s at K=4096 while the true shipped point measures 1239.2, so
the guard silently guaranteed nothing.

The only honest shipped measurement is a run with NO overrides at all: no
WIREHAIR_V2_PEEL_DEGREES, no scale flag. That is what this does.

WHAT IT PRODUCES. For each K: the trained point and the shipped point, measured
back to back at the same trial count and payload, and a verdict. Any K where
training did not beat shipping is rewritten to shipped, so the table cannot ship
a regression. The output doubles as the stage-4 comparison data.

Regressions are EXPECTED at some K and are not a bug in the codec: above about
K=2048 the working set leaves cache (K=4096 x 4096 B is 16 MB of intermediates)
and throughput falls by more than half for shipped and trained alike, which
compresses the margin the search is trying to find.
"""
import argparse
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from peel_funnel import Funnel                              # noqa: E402


def stock_pmf(k):
    w = [1.0 / k] + [1.0 / (d * (d - 1)) for d in range(2, 64)]
    w.append(max(0.0, 1.0 - sum(w)))
    return w


def run(bench, k, trials, bb, env_extra=None):
    """(fail, OH_mean, decode_MBps) for one arm."""
    env = dict(os.environ)
    env.pop("WIREHAIR_V2_PEEL_DEGREES", None)
    if env_extra:
        env.update(env_extra)
    p = subprocess.run(
        [bench, "compare", "--nlo", str(k), "--nhi", str(k), "--bb-list", str(bb),
         "--trials", str(trials), "--loss", "0.10", "--precode",
         "--precode-profile", "mixed", "--mixed-gf16-rows", "0"],
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", default="build-fast/codec/wirehair_v2_bench")
    ap.add_argument("--table", required=True)
    ap.add_argument("--out", default=None, help="write the corrected table here")
    ap.add_argument("--trials", type=int, default=3000)
    ap.add_argument("--bb", type=int, default=4096)
    ap.add_argument("--kmax", type=int, default=10 ** 9)
    ap.add_argument("--margin", type=float, default=2.0,
                    help="percent goodput a trained point must beat shipped by")
    a = ap.parse_args()

    table = {int(k): v for k, v in json.load(open(a.table)).items()}
    ks = [k for k in sorted(table) if k <= a.kmax]

    print(f"  {len(ks)} K values, {a.trials} trials each arm, bb={a.bb}\n")
    print(f"  {'K':>6} {'trained':>9} {'shipped':>9} {'delta':>8} "
          f"{'tr OH':>7} {'sh OH':>7} {'verdict':>10}")
    fixed, wins, losses = {}, 0, 0
    for k in ks:
        v = table[k]
        w = Funnel.family(stock_pmf(k), v["p1"], v["tilt"], v["dmax"], v["absorb"])
        tr = run(a.bench, k, a.trials, a.bb,
                 {"WIREHAIR_V2_PEEL_DEGREES": ",".join(f"{x:.9f}" for x in w)})
        sh = run(a.bench, k, a.trials, a.bb)          # NO overrides: true shipped
        if tr is None or sh is None:
            print(f"  {k:>6}   measurement failed"); continue
        tf, toh, tmb = tr
        sf, soh, smb = sh
        tg = tmb * k / (k + toh) if tf == 0 else 0.0
        sg = smb * k / (k + soh) if sf == 0 else 0.0
        # A MARGIN, not a bare comparison. Single-shot goodput carries a few
        # percent of noise, so "trained > shipped by 0.2%" is a coin flip -- and
        # it is not a free coin flip, because the trained point can carry real
        # extra overhead (K=11 won by +0.2% while needing 1.3% of K in extra
        # packets). Shipped is the safer default, so training has to clear a
        # margin to displace it.
        keep = tg > sg * (1.0 + a.margin / 100.0)
        d = 100.0 * (tg - sg) / sg if sg > 0 else 0.0
        if keep:
            wins += 1; fixed[k] = dict(v, verified_mbps=tmb, verified_oh=toh,
                                       shipped_mbps=smb, gain_pct=round(d, 2))
        else:
            losses += 1
            fixed[k] = {"K": k, "p1": 100, "tilt": 0, "dmax": 64, "absorb": 100,
                        "reverted_to_shipped": True, "verified_mbps": smb,
                        "verified_oh": soh, "shipped_mbps": smb,
                        "gain_pct": 0.0, "search_would_have_lost_pct": round(d, 2)}
        print(f"  {k:>6} {tmb:>9.1f} {smb:>9.1f} {d:>+7.1f}% {toh:>7.4f} "
              f"{soh:>7.4f} {'keep' if keep else 'REVERT':>10}"
              + ("" if tf == 0 and sf == 0 else f"  fail tr={tf} sh={sf}"), flush=True)
        if a.out:
            with open(a.out, "w") as fh:
                json.dump({str(kk): vv for kk, vv in sorted(fixed.items())}, fh, indent=1)

    g = [v["gain_pct"] for v in fixed.values() if not v.get("reverted_to_shipped")]
    print(f"\n  {wins} kept, {losses} reverted to shipped")
    if g:
        g.sort()
        print(f"  gains where kept: median {g[len(g)//2]:+.1f}%  "
              f"min {g[0]:+.1f}%  max {g[-1]:+.1f}%")
    if a.out:
        print(f"  corrected table -> {a.out}")


if __name__ == "__main__":
    main()
