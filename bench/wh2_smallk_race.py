#!/usr/bin/env python3
"""Race Wirehair 2 against Wirehair 1 over one K range, untuned on both sides.

Protocol (user-set):
  * sample 10% of the K values in the range, 3 construction seeds each,
  * NO per-K seed tuning on either codec -- one uniform construction seed,
  * each K range is measured on its own, because the ranges are structurally
    different regimes and are expected to need different code.

Throughput is encoder init plus emitting the first K systematic symbols, at the
minimal block size, single threaded, median of --reps. Overhead is measured on
the same sampled cells: the fraction of (K, seed) cells that fail to decode at
zero extra symbols, plus the mean extra symbols needed.

The point of the harness is the RATIO: WH2 nanoseconds over WH1 nanoseconds per
K. Below 1.0 means WH2 is faster, which is the goal.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import statistics
import subprocess
import sys


def sample_k(lo: int, hi: int, pct: float, seed: int) -> list[int]:
    """Deterministic ~pct% sample of [lo, hi], always including both ends."""
    span = list(range(lo, hi + 1))
    want = max(2, round(len(span) * pct / 100.0))
    if want >= len(span):
        return span
    # SplitMix64 keyed by the sample seed, so the sample is reproducible and
    # does not depend on Python's hash seed.
    state = seed & 0xFFFFFFFFFFFFFFFF
    def nxt() -> int:
        nonlocal state
        state = (state + 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
        z = state
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
        return z ^ (z >> 31)
    pool = span[1:-1]
    picked = {span[0], span[-1]}
    while len(picked) < want and pool:
        picked.add(pool.pop(nxt() % len(pool)))
    return sorted(picked)


def run(cmd: list[str]) -> str:
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise SystemExit(f"command failed ({p.returncode}): {' '.join(cmd)}\n{p.stderr[:2000]}")
    return p.stdout


def wh2_timing(bench: str, ks: list[int], bb: int, reps: int, flags: list[str]) -> dict[int, int]:
    out = run([bench, "precodefail", "--N", ",".join(map(str, ks)), "--bb-list", str(bb),
               "--overhead", "0", "--trials", "1", "--threads", "1", "--loss", "0",
               "--seed", "1", "--encode-timing", str(reps)] + flags)
    got = {}
    for line in out.splitlines():
        if not line.startswith("# encode_timing,N="):
            continue
        f = dict(kv.split("=", 1) for kv in line[2:].split(",") if "=" in kv)
        got[int(f["N"])] = int(f["median_ns"])
    return got


def wh1_timing(whx: str, ks: list[int], bb: int, reps: int) -> dict[int, int]:
    out = run([whx, "enctime", "--nlist", ",".join(map(str, ks)),
               "--reps", str(reps), "--bb", str(bb)])
    got = {}
    for line in out.splitlines():
        if line.startswith("#") or line.startswith("K,"):
            continue
        c = line.split(",")
        if len(c) > 3 and c[0].isdigit():
            got[int(c[0])] = int(c[3])
    return got


def wh2_overhead(bench: str, ks: list[int], bb: int, seeds: list[int], loss: float,
                 flags: list[str]) -> dict:
    """Untuned single-attempt cells: one construction seed per pass."""
    fails, extras = 0, []
    for cs in seeds:
        out = run([bench, "precodefail", "--N", ",".join(map(str, ks)), "--bb-list", str(bb),
                   "--overhead", "0,1,2,3,4", "--trials", "1", "--threads", "1",
                   "--loss", repr(loss), "--seed", "1", "--schedule", "iid",
                   "--paired-overhead-stream", "--overhead-early-stop",
                   "--construction-seed", "0x%x" % cs] + flags)
        per_k: dict[int, int | None] = {}
        for line in out.splitlines():
            if line.startswith("#") or line.startswith("N,"):
                continue
            c = line.split(",")
            if len(c) < 10 or not c[0].isdigit():
                continue
            k, oh, success = int(c[0]), int(c[4]), int(c[6])
            if success and per_k.get(k) is None:
                per_k[k] = oh
        for k in ks:
            v = per_k.get(k)
            if v is None:
                fails += 1
                extras.append(4)   # censored at the cap, same rule as the sealed campaign
            else:
                extras.append(v)
    n = len(ks) * len(seeds)
    return {"cells": n, "censored": fails,
            "fail_at_0": sum(1 for e in extras if e > 0) / n,
            "mean_extra": sum(extras) / n}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bench", required=True, help="wirehair_v2_bench binary")
    ap.add_argument("--whx", required=True, help="whx binary (Wirehair 1 arm)")
    ap.add_argument("--k-min", type=int, default=1)
    ap.add_argument("--k-max", type=int, default=100)
    ap.add_argument("--pct", type=float, default=10.0, help="percent of K values to sample")
    ap.add_argument("--seeds", type=int, default=3, help="construction seeds per K")
    ap.add_argument("--sample-seed", type=lambda s: int(s, 0), default=0x6789CAFE)
    ap.add_argument("--bb", type=int, default=2)
    ap.add_argument("--reps", type=int, default=9)
    ap.add_argument("--loss", type=float, default=0.1)
    ap.add_argument("--skip-overhead", action="store_true")
    ap.add_argument("--label", default="wh2")
    ap.add_argument("--json-out", type=pathlib.Path)
    ap.add_argument("flags", nargs="*", help="WH2 config flags, after --")
    a = ap.parse_args()

    # Both benches take N >= 2; K=1 is degenerate (the single source block is
    # the only systematic symbol) and neither arm measures it.
    k_min = max(2, a.k_min)
    if k_min != a.k_min:
        print(f"# note: K range clamped to {k_min}..{a.k_max} (benches require N >= 2)")
    ks = sample_k(k_min, a.k_max, a.pct, a.sample_seed)
    seeds = [(a.sample_seed + i * 0x9E3779B9) & 0xFFFFFFFF for i in range(a.seeds)]
    flags = [f for f in a.flags if f != "--"]

    t2 = wh2_timing(a.bench, ks, a.bb, a.reps, flags)
    t1 = wh1_timing(a.whx, ks, a.bb, a.reps)
    common = [k for k in ks if k in t1 and k in t2]
    if not common:
        raise SystemExit("no K measured on both arms")

    ratios = [t2[k] / t1[k] for k in common]
    print(f"# range K={a.k_min}..{a.k_max}  sampled {len(ks)} K ({a.pct}%)  "
          f"x {a.seeds} seeds  bb={a.bb}  reps={a.reps}")
    print(f"# config: {' '.join(flags) or '(default)'}")
    print(f"{'K':>6} {'WH1 ns':>9} {'WH2 ns':>9} {'WH2/WH1':>8}")
    for k in common:
        print(f"{k:>6} {t1[k]:>9} {t2[k]:>9} {t2[k]/t1[k]:>8.2f}")
    med = statistics.median(ratios)
    print(f"\n{a.label}: median WH2/WH1 = {med:.3f}   "
          f"({'FASTER' if med < 1 else 'slower'} than WH1)   "
          f"best {min(ratios):.2f} worst {max(ratios):.2f}")
    # aggregate: total time to build every sampled generation
    agg = sum(t2[k] for k in common) / sum(t1[k] for k in common)
    print(f"{a.label}: aggregate WH2/WH1 = {agg:.3f}")

    result = {"k_min": a.k_min, "k_max": a.k_max, "ks": common, "flags": flags,
              "wh1_ns": {str(k): t1[k] for k in common},
              "wh2_ns": {str(k): t2[k] for k in common},
              "median_ratio": med, "aggregate_ratio": agg}

    if not a.skip_overhead:
        oh = wh2_overhead(a.bench, ks, a.bb, seeds, a.loss, flags)
        print(f"{a.label}: overhead over {oh['cells']} untuned cells -- "
              f"fail@+0 {oh['fail_at_0']:.4f}, mean extra {oh['mean_extra']:.4f}, "
              f"censored {oh['censored']}")
        result["overhead"] = oh

    if a.json_out:
        a.json_out.write_text(json.dumps(result, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
