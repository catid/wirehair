#!/usr/bin/env python3
"""Search individual small K values using measurements from the real codec.

Every compare arm receives independently derived construction and loss seeds
paired within its tier.  The stock target PMF is always included in the final
ranking, and output is published
atomically in a versioned table only after every requested K has a result.
"""
import argparse
import math
import os
import random
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from peel_codec import (                                # noqa: E402
    MeasurementError,
    capture_artifact_identity,
    compare_probe,
    derive_seed,
    family,
    make_table_document,
    make_search_receipt,
    stock_profile,
    valid_loss_rate,
    write_json_atomic,
)

BOX = [("p1", 0, 400), ("tilt", -100, 1600), ("dmax", 2, 64), ("absorb", 0, 100)]

class Direct:
    def __init__(self, a):
        self.a = a
        self.probes = 0
        self.profile = None

    def probe(self, k, v, trials, bb, tier):
        """Return full metrics under the tier's explicit paired seeds."""
        construction_seed = derive_seed(
            self.a.construction_seed, "direct-search", k, tier, trials, bb,
            "construction")
        loss_seed = derive_seed(
            self.a.loss_seed, "direct-search", k, tier, trials, bb, "loss")
        exact = {
            "native_profile": self.profile,
            "target_profile": self.a.target_profile,
            "seed_policy": self.a.seed_policy,
            "loss": self.a.loss,
            "schedule": self.a.schedule,
            "construction_seed": construction_seed,
            "loss_seed": loss_seed,
        }
        if v is None:
            self.probes += 1
            return compare_probe(
                self.a.bench, k, trials, bb, **exact)
        w = family(
            self.profile.pmf, v[0], v[1], v[2], v[3])
        if w is None:
            return None
        self.probes += 1
        # This search intentionally leaves the staircase scale unset.  Clearing
        # every inherited WIREHAIR_V2_* hook in compare_probe makes that an
        # exact arm rather than whatever the launching shell happened to carry.
        return compare_probe(
            self.a.bench, k, trials, bb, peel_weights=w, **exact)

    def candidate_pmf(self, vector):
        return family(
            self.profile.pmf,
            vector[0], vector[1], vector[2], vector[3])

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

    def solve(self, k, centre):
        if centre is not None and (
                not isinstance(centre, list) or len(centre) != len(BOX) or
                any(isinstance(value, bool) or not isinstance(value, int)
                    for value in centre)):
            raise ValueError("warm start must contain four integers")
        self.profile = stock_profile(
            self.a.bench, k, target_profile=self.a.target_profile)
        sampling_seed = derive_seed(
            self.a.construction_seed, "direct-search", k, "sampling")
        rng = random.Random(sampling_seed)
        t0 = time.time()
        raw_pool = self.lhs(rng, self.a.screen, centre)
        seen_pmfs = {tuple(self.profile.pmf)}
        pool = []
        for vector in raw_pool:
            candidate_pmf = self.candidate_pmf(vector)
            if candidate_pmf is None:
                continue
            key = tuple(candidate_pmf)
            if key in seen_pmfs:
                continue
            seen_pmfs.add(key)
            pool.append(vector)
        # gate everything cheaply, keep the decoders
        alive = []
        for v in pool:
            r = self.probe(
                k, v, self.a.gate_trials, self.a.gate_bb, "gate")
            if r and r.fail == 0:
                alive.append((r.goodput(k), v))
        alive.sort(key=lambda x: -x[0])
        best_g, best_v = alive[0] if alive else (0.0, None)

        # local refine, still on the cheap tier
        r_frac, used = 0.25, 0
        while best_v is not None and r_frac > 0.02 and used < self.a.refine:
            moved = False
            for j, (_, lo, hi) in enumerate(BOX):
                step = max(1, int(round(r_frac * (hi - lo))))
                for s in (-step, step):
                    if used >= self.a.refine:
                        break
                    w = self.clamp([x + (s if jj == j else 0)
                                    for jj, x in enumerate(best_v)])
                    if w == best_v:
                        continue
                    candidate_pmf = self.candidate_pmf(w)
                    if candidate_pmf is None:
                        continue
                    key = tuple(candidate_pmf)
                    if key in seen_pmfs:
                        continue
                    seen_pmfs.add(key)
                    used += 1
                    rr = self.probe(
                        k, w, self.a.gate_trials, self.a.gate_bb, "gate")
                    if rr and rr.fail == 0:
                        g = rr.goodput(k)
                        if g > best_g:
                            best_g, best_v, moved = g, w, True
                if used >= self.a.refine:
                    break
            if not moved:
                r_frac *= 0.5

        # rank tier: real payload, more trials, on the top few
        # The shipped law is a mandatory control in the final ranking.
        # None is the true shipped codec with no overrides.  Feeding the
        # identity weights through the environment hook is equation-identical
        # but measurably biased by the alternate code path.
        incumbent = None
        cands = (
            ([best_v] if best_v is not None else []) +
            [v for _, v in alive[:self.a.rank_top] if v != best_v]
        )
        cands = cands[:self.a.rank_top] + [incumbent]
        out, rank_measurements = [], []
        for v in cands:
            r = self.probe(
                k, v, self.a.rank_trials, self.a.rank_bb, "rank")
            if r is not None:
                rank_measurements.append((r, v))
            if r and r.fail == 0:
                out.append((r.goodput(k), r, v))
        if not out:
            return None
        # The incumbent wins an exact tie; a trained arm must demonstrate a
        # strict goodput improvement to become the recorded search winner.
        out.sort(key=lambda x: (-x[0], x[2] is not None))
        g, measurement, v = out[0]
        shipped_control = next(
            (rank_metrics for rank_metrics, vector in rank_measurements
             if vector is None),
            None)
        if shipped_control is None:
            raise MeasurementError("mandatory shipped rank control is missing")
        if v is None:
            result = {
                "K": k, "scale": -1.0, "p1": 100, "tilt": 0,
                "dmax": 64, "absorb": 100, "reverted_to_shipped": True,
            }
        else:
            result = {
                "K": k, "scale": -1.0, "p1": v[0], "tilt": v[1],
                "dmax": v[2], "absorb": v[3],
            }
        mode = "shipped" if v is None else "trained"
        profile = self.profile
        selected_pmf = (
            list(profile.pmf) if v is None else
            family(profile.pmf, v[0], v[1], v[2], v[3])
        )
        if selected_pmf is None:
            return None
        result.update(
            goodput=g,
            **measurement.as_dict(),
            native=profile.as_dict(),
            peel_pmf=selected_pmf,
            search_receipt=make_search_receipt(
                measurement,
                mode=mode,
                goodput=g,
                trials=self.a.rank_trials,
                block_bytes=self.a.rank_bb,
                search_kind="direct-real-codec",
                construction_seed_base=self.a.construction_seed,
                loss_seed_base=self.a.loss_seed,
                seed_domain="direct-search",
                coordinates={
                    name: result[name]
                    for name in ("scale", "p1", "tilt", "dmax", "absorb")
                },
                peel_pmf=selected_pmf,
                shipped_control=shipped_control,
                shipped_goodput=shipped_control.goodput(k),
                context={
                    "warm_start": list(centre) if centre is not None else None,
                    "sampling_seed": sampling_seed,
                    "screen": self.a.screen,
                    "refine": self.a.refine,
                    "gate_trials": self.a.gate_trials,
                    "gate_block_bytes": self.a.gate_bb,
                    "rank_top": self.a.rank_top,
                },
            ),
            seconds=round(time.time() - t0, 1), probes=self.probes)
        return result


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", default="build-fast/codec/wirehair_v2_bench")
    ap.add_argument("--kmin", type=int, default=2)
    ap.add_argument("--kmax", type=int, default=128)
    ap.add_argument("--out", default="tools/peel_table_small.json")
    ap.add_argument("--screen", type=int, default=40)
    ap.add_argument("--refine", type=int, default=48)
    ap.add_argument("--gate-trials", type=int, default=400)
    ap.add_argument("--gate-bb", type=int, required=True)
    ap.add_argument("--rank-trials", type=int, default=2000)
    ap.add_argument("--rank-bb", type=int, required=True)
    ap.add_argument("--rank-top", type=int, default=3)
    ap.add_argument("--target-profile", required=True, choices=["dispatch-v1"])
    ap.add_argument("--seed-policy", required=True, choices=["raw"])
    ap.add_argument("--schedule", required=True, choices=["iid"])
    ap.add_argument("--loss", type=float, required=True)
    ap.add_argument("--construction-seed", type=int, required=True)
    ap.add_argument("--loss-seed", type=int, required=True)
    a = ap.parse_args(argv)
    if (a.kmin < 2 or a.kmax > 64000 or a.kmax < a.kmin or
            a.screen < 1 or a.refine < 0 or
            a.gate_trials < 1 or a.rank_trials < 1 or a.gate_bb < 1 or
            a.rank_bb < 1 or a.rank_top < 1 or
            not valid_loss_rate(a.loss) or
            not 0 <= a.construction_seed <= 0xffffffffffffffff or
            not 0 <= a.loss_seed <= 0xffffffffffffffff):
        ap.error(
            "invalid range, budget, payload, rank-top, loss, or uint64 seed")

    try:
        identity = capture_artifact_identity(
            a.bench, "tools/peel_direct.py")
    except MeasurementError as error:
        print(f"  REFUSED measurement: {error}", file=sys.stderr)
        return 2
    table, centre = {}, None
    failed = []
    # walk upward so each K warm-starts from the one below it
    for k in range(a.kmin, a.kmax + 1):
        d = Direct(a)
        try:
            r = d.solve(k, centre)
        except (MeasurementError, OSError, ValueError) as error:
            print(
                f"  REFUSED publication: measurement failed for K={k}: "
                f"{error}", file=sys.stderr)
            return 1
        if r is None:
            print(f"  K={k:<5} no candidate decoded", flush=True)
            failed.append(k)
            continue
        table[k] = r
        centre = None if r.get("reverted_to_shipped") else [
            r["p1"], r["tilt"], r["dmax"], r["absorb"]]
        print(f"  K={k:<5} p1={r['p1']:>3} tilt={r['tilt']:>+5} dmax={r['dmax']:>2} "
              f"absorb={r['absorb']:>3}  {r['decode_mbps']:>7.1f}MB/s "
              f"OH={r['oh_mean']:.4f} "
              f"({r['seconds']}s, {r['probes']} probes)", flush=True)
    if failed:
        print(
            f"\n  REFUSED publication: no result for K values {failed}",
            file=sys.stderr)
        return 1
    try:
        document = make_table_document(
            table,
            generator="tools/peel_direct.py",
            bench=a.bench,
            construction_seed_base=a.construction_seed,
            loss_seed_base=a.loss_seed,
            target_profile=a.target_profile,
            seed_policy=a.seed_policy,
            loss=a.loss,
            schedule=a.schedule,
            settings={
                "kmin": a.kmin,
                "kmax": a.kmax,
                "screen": a.screen,
                "refine": a.refine,
                "gate_trials": a.gate_trials,
                "gate_block_bytes": a.gate_bb,
                "rank_trials": a.rank_trials,
                "rank_block_bytes": a.rank_bb,
                "rank_top": a.rank_top,
                "target_profile": a.target_profile,
                "seed_policy": a.seed_policy,
                "loss": a.loss,
                "schedule": a.schedule,
            },
            artifact_identity=identity,
        )
        write_json_atomic(a.out, document)
    except (MeasurementError, OSError, ValueError) as error:
        print(f"  REFUSED publication: {error}", file=sys.stderr)
        return 1
    print(f"\n  {len(table)} K values written to {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
