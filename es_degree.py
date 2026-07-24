#!/usr/bin/env python3
"""Learn the staircase row-degree distribution by evolution strategies.

Objective: mean inactivated-column count at one K.  It is the right target
because R sets the dimension of three phases at once -- projection, back
substitution, and the (R-H)*H term of the GF(256) muladd model -- and unlike a
timing measurement it is structural, so it is immune to CPU contention.

Two variance-reduction choices matter more than population size here:

  * Common random numbers.  Every candidate in a generation is evaluated on the
    SAME construction seeds and the SAME loss-pattern seed, so the comparison is
    paired and the ES gradient measures shape differences rather than sampling
    noise.
  * Antithetic sampling.  Perturbations are evaluated in +/- pairs, which
    cancels the first-order noise term in the gradient estimate.

The search space is the SHAPE only.  The precode is column-regular (every source
column carries exactly SourceHits edges), so the edge budget is fixed and the
codec rescales any proposed degree vector onto it; a scale direction in the
parameters would therefore be a flat direction in the objective.  Parameters are
control points on a log-degree grid, expanded to dense per-degree weights.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import math
import os
import pathlib
import subprocess
import numpy as np

BENCH = os.environ.get("WH_BENCH", "./build-fast/codec/wirehair_v2_bench")
H8 = ["--completion", "mixed", "--mixed-geometry", "frozen", "--mixed-period", "244",
      "--mixed-gf256-rows", "8", "--mixed-gf16-rows", "0",
      "--binary-dense-rows", "1", "--source-hits", "1"]


def shape_to_weights(theta: np.ndarray, K: int, max_degree: int) -> str:
    """Control points on a log-degree grid -> dense weights for degrees 1..D."""
    D = min(max_degree, max(2, K))
    knots = np.geomspace(1.0, max(2.0, float(D)), num=len(theta))
    grid = np.arange(1, D + 1, dtype=float)
    # softmax over interpolated logits keeps weights positive and scale-free
    logits = np.interp(np.log(grid), np.log(knots), theta)
    logits -= logits.max()
    w = np.exp(logits)
    w /= w.sum()
    return ",".join(f"{x:.6g}" for x in w)


def evaluate(args) -> float:
    """Mean inactivation for one shape at one K, averaged over the CRN seed set."""
    theta, K, seeds, trials, max_degree, loss = args
    weights = shape_to_weights(np.asarray(theta), K, max_degree)
    env = dict(os.environ)
    env["WIREHAIR_V2_BAND_TRACKING_X"] = "1"
    env["WIREHAIR_V2_STAIRCASE_DEGREES"] = weights
    total, n = 0.0, 0
    for cs in seeds:
        p = subprocess.run(
            [BENCH, "precodefail", "--N", str(K), "--bb-list", "2", "--overhead", "0",
             "--trials", str(trials), "--threads", "1", "--loss", repr(loss),
             "--seed", "9", "--construction-seed", hex(cs)] + H8,
            capture_output=True, text=True, env=env)
        for line in p.stdout.splitlines():
            if line.startswith("#") or line.startswith("N,"):
                continue
            f = line.split(",")
            if len(f) > 11 and f[0].isdigit():
                total += float(f[10])   # inact_mu
                n += 1
    return total / n if n else float("inf")


def run_es(K: int, dims: int, gens: int, pop: int, sigma: float, lr: float,
           seeds: list[int], trials: int, max_degree: int, loss: float,
           pool: cf.ProcessPoolExecutor, log) -> tuple[np.ndarray, float, float]:
    rng = np.random.default_rng(12345 + K)
    theta = np.zeros(dims)
    base = evaluate((theta.copy(), K, seeds, trials, max_degree, loss))
    best_theta, best_val = theta.copy(), base
    half = pop // 2
    for g in range(gens):
        eps = rng.normal(size=(half, dims))
        cand = [theta + sigma * e for e in eps] + [theta - sigma * e for e in eps]
        vals = np.array(list(pool.map(
            evaluate, [(c, K, seeds, trials, max_degree, loss) for c in cand])))
        # centered-rank fitness: scale-free and robust to outliers
        order = vals.argsort()                       # lower inactivation is better
        ranks = np.empty(pop); ranks[order] = np.arange(pop)
        adv = (ranks / (pop - 1) - 0.5)              # best -> -0.5, worst -> +0.5
        grad = (adv[:half, None] * eps).sum(0) - (adv[half:, None] * eps).sum(0)
        theta -= lr / (half * sigma) * grad          # descend: minimize inactivation
        cur = evaluate((theta.copy(), K, seeds, trials, max_degree, loss))
        if cur < best_val:
            best_theta, best_val = theta.copy(), cur
        if g % 10 == 0 or g == gens - 1:
            log(f"  K={K:<4} gen {g:>3}  current {cur:.3f}  best {best_val:.3f}  "
                f"(flat baseline {base:.3f})")
    return best_theta, best_val, base


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ks", default="2,4,8,16,32,64,128")
    ap.add_argument("--dims", type=int, default=10)
    ap.add_argument("--gens", type=int, default=60)
    ap.add_argument("--pop", type=int, default=32)
    ap.add_argument("--sigma", type=float, default=0.25)
    ap.add_argument("--lr", type=float, default=0.30)
    ap.add_argument("--seeds", type=int, default=6, help="CRN construction seeds")
    ap.add_argument("--trials", type=int, default=300)
    ap.add_argument("--max-degree", type=int, default=48)
    ap.add_argument("--loss", type=float, default=0.1)
    ap.add_argument("--workers", type=int, default=max(1, (os.cpu_count() or 8) - 2))
    ap.add_argument("--out", type=pathlib.Path, required=True)
    a = ap.parse_args()

    seeds = [int((0x9E3779B9 * (i + 1)) & 0xFFFFFFFF) for i in range(a.seeds)]
    out = {"config": vars(a) | {"out": str(a.out)}, "results": {}}
    print(f"ES over {a.dims} shape dims, pop {a.pop}, {a.gens} gens, "
          f"{a.seeds} CRN seeds x {a.trials} trials, {a.workers} workers")
    with cf.ProcessPoolExecutor(max_workers=a.workers) as pool:
        for K in [int(x) for x in a.ks.split(",")]:
            theta, val, base = run_es(K, a.dims, a.gens, a.pop, a.sigma, a.lr,
                                      seeds, a.trials, a.max_degree, a.loss,
                                      pool, print)
            weights = shape_to_weights(theta, K, a.max_degree)
            stock = evaluate((None,) and (np.zeros(a.dims), K, seeds, a.trials,
                                          a.max_degree, a.loss))
            out["results"][str(K)] = {
                "theta": theta.tolist(), "weights": weights,
                "best_inact": val, "flat_baseline": base, "flat_recheck": stock,
            }
            print(f"K={K}: {base:.3f} -> {val:.3f}  ({100*(base-val)/base:+.1f}%)")
            a.out.write_text(json.dumps(out, indent=1))
    print(f"wrote {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
