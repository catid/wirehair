#!/usr/bin/env python3
"""Screen peel candidates with native counters, then rank on the real codec.

The counter stage is only a candidate ordering. Recovery and throughput come
from explicit common-random-number ``compare`` pairs within each tier. The
native ``peelpmf``
receipt supplies both the exact shipped distribution and the staircase size;
this module deliberately does not duplicate either construction in Python.
The shipped arm is always measured in both the cheap gate and final ranking.
"""
import argparse
import json
import math
import random
import subprocess
import sys
import time

from peel_codec import (
    MeasurementError,
    compare_probe,
    derive_seed,
    family as peel_family,
    isolated_codec_env,
    stock_profile,
)

BENCH_DEFAULT = "build-fast/codec/wirehair_v2_bench"
PROXY_COST_MODEL = "embedded-es-cost-model/raw-calibration-unavailable"
SEARCH_BOX_PROTOCOL = "native-mean:[0,min(K,max(40,4*mean))]/v1"
PROXY_MEASURE_REGIME = {
    "solve_block_bytes": 0,
    "cost_model_block_bytes": 1280,
    "cost_model_verified": 0,
    "band_tracking_x": 1,
    "loss": "0.100000",
    "seed_base": 55,
    "completion": "mixed",
    "geometry": "frozen",
    "period": 244,
    "gf16_rows": 0,
}
PROXY_ORDERING_PROTOCOL = "fail-rate-then-pred-ns/v1"
FUNNEL_RESULT_SCHEMA = "wirehair-v2-peel-funnel-result/v1"
FUNNEL_RESULT_PREFIX = "# peel_funnel_result="
MEASURE_COLUMNS = (
    "S", "H", "D2", "scale", "shape", "p1", "p2", "p3", "c1", "c2",
    "c3", "c4", "c5", "dmax", "peel", "peel_p1", "peel_tilt",
    "peel_dmax", "peel_absorb", "cells", "failures", "fail_rate",
    "solved", "pred_ns", "xors", "muladds", "copies", "zerofills", "status",
)

NON_SCALE_BOX = [
    ("p1", 0, 400),         # hundredths of the incumbent's degree-1 mass
    ("tilt", 0, 1700),      # offset lattice; signed = value - 100
    ("dmax", 2, 64),
    ("absorb", 0, 100),
]


def search_box(profile):
    """Build a broad scale range that contains the native-density region."""
    scale_max = min(
        float(profile.block_count),
        max(40.0, 4.0 * profile.shipped_mean))
    return [("scale", 0, int(math.ceil(100.0 * scale_max)))] + NON_SCALE_BOX


def proxy_order(result):
    """Order proxy rows by recovery surrogate first, predicted work second."""
    return result[1], result[0]


class Funnel:
    def __init__(self, args):
        if not getattr(args, "allow_unverified_cost_model", False):
            raise MeasurementError(
                f"proxy cost model {PROXY_COST_MODEL!r} requires explicit "
                "allow_unverified_cost_model opt-in")
        self.a = args
        self.k = args.K
        self.native_profile = stock_profile(args.bench, args.K)
        self.box = search_box(self.native_profile)
        self.struct = f"{self.native_profile.staircase}:10:12"
        self.calls = 0

    def token(self, v):
        s, h, d2 = self.struct.split(":")
        return (f"{s},{h},{d2},{v[0]/100.0:.2f},0,20,20,20,100,100,100,100,100,16,1,"
                f"{v[1]},{v[2]-100},{v[3]},{v[4]}")

    def measure(self, vecs, cells):
        """[(predicted_cost, fail_percent)] -- columns read BY HEADER NAME.

        The proxy failure surrogate is the primary finite-screen ordering and
        predicted work breaks ties. The later real-codec gate remains the
        authority for actual recovery.
        """
        out = []
        if (cells < 1 or self.a.cell_base < 0 or
                self.a.cell_base + cells > 0x100000000):
            raise MeasurementError(
                "essearch cell range must fit the uint32 cell domain")
        for i in range(0, len(vecs), self.a.batch):
            chunk = vecs[i:i + self.a.batch]
            self.calls += 1
            command = [
                self.a.bench, "essearch", "--N", str(self.k),
                 # The native receipt pins the actual structure. Keep the broad
                 # parser domain so essearch accepts every production S.
                 "--fix-structure", self.struct, "--s-lo", "1", "--s-hi", "400",
                 # CmdEsSearch requires these fields even in --measure mode.
                 # They do not enter the rows, and this caller ranks pred_ns
                 # alone, so spell an explicitly inert objective.
                 "--target", "0.5", "--lam", "0",
                 "--measure", "--threads", str(self.a.threads),
                 "--cells", f"{self.a.cell_base}:{self.a.cell_base + cells}",
                 "--configs", ";".join(self.token(v) for v in chunk),
                 "--allow-unverified-cost-model"]
            try:
                p = subprocess.run(
                    command, capture_output=True, text=True,
                    env=isolated_codec_env())
            except OSError as error:
                raise MeasurementError(
                    f"could not execute essearch benchmark: {error}")
            if p.returncode != 0:
                detail = p.stderr.strip() or p.stdout.strip() or "no output"
                raise MeasurementError(
                    f"essearch --measure exited {p.returncode}: {detail}")
            hdr, got, banner = None, [], None
            expected_banner = (
                f"# essearch measure,N={self.k},"
                f"cells=[{self.a.cell_base},{self.a.cell_base + cells}),"
                "solve_bb=0,cost_model_bb=1280,cost_model_verified=0,"
                "band_tracking_x=1,loss=0.100000,seed_base=55,"
                "completion=mixed,geometry=frozen,period=244,gf16_rows=0,"
                f"threads={self.a.threads}"
            )
            for line in p.stdout.splitlines():
                if not line.strip():
                    continue
                if line.startswith("# essearch measure,"):
                    if banner is not None or line != expected_banner:
                        raise MeasurementError(
                            f"essearch emitted an unexpected measure regime: "
                            f"{line}")
                    banner = line
                    continue
                if line.startswith("S,H,D2,scale,") and "fail_rate" in line:
                    parsed_header = tuple(c.strip() for c in line.split(","))
                    if (banner is None or hdr is not None or got or
                            parsed_header != MEASURE_COLUMNS):
                        raise MeasurementError(
                            "essearch emitted an unexpected measure header")
                    hdr = list(parsed_header)
                    continue
                if hdr and line[:1].isdigit() and line.count(",") > 10:
                    fields = [field.strip() for field in line.split(",")]
                    try:
                        if len(fields) != len(hdr):
                            raise ValueError("wrong field count")
                        expected = self.token(chunk[len(got)]).split(",")
                        if fields[:len(expected)] != expected:
                            raise ValueError("echoed configuration mismatch")
                        measured_cells = int(fields[hdr.index("cells")])
                        failures = int(fields[hdr.index("failures")])
                        solved = int(fields[hdr.index("solved")])
                        status = fields[hdr.index("status")]
                        fail = float(fields[hdr.index("fail_rate")])
                        pred = float(fields[hdr.index("pred_ns")])
                        counters = [
                            float(fields[hdr.index(name)])
                            for name in (
                                "xors", "muladds", "copies", "zerofills")
                        ]
                    except (ValueError, IndexError) as error:
                        raise MeasurementError(
                            f"malformed essearch measure row ({error}): {line}")
                    expected_fail = (
                        failures / measured_cells
                        if measured_cells > 0 else math.nan)
                    if (measured_cells != cells or measured_cells < 1 or
                            not 0 <= failures <= measured_cells or
                            not 0 <= solved <= measured_cells or
                            failures + solved != measured_cells or
                            status not in ("ok", "rejected") or
                            not math.isfinite(fail) or
                            not math.isclose(
                                fail, expected_fail,
                                rel_tol=0.0, abs_tol=5.1e-9) or
                            not math.isfinite(pred) or pred < 0.0 or
                            any(not math.isfinite(value) or value < 0.0
                                for value in counters)):
                        raise MeasurementError(
                            f"invalid essearch measure row: {line}")
                    if status != "ok" or solved == 0:
                        got.append((None, None))
                        continue
                    got.append((pred, 100.0 * fail))
                    continue
                raise MeasurementError(
                    f"unexpected essearch measure output: {line}")
            if banner is None or hdr is None:
                raise MeasurementError(
                    "essearch emitted an incomplete measure-regime receipt")
            if len(got) != len(chunk):
                raise MeasurementError(
                    f"essearch returned {len(got)} rows for "
                    f"{len(chunk)} configurations")
            out.extend(got)
        return out

    def lhs(self, rng, n, centre=None, frac=0.20):
        """Latin hypercube over the box, or over a local window around `centre`.

        With a warm start, HALF the screen is drawn from a window of +/- frac
        around the neighbouring K's answer and half from the full box.
        """
        def block(m, lows, highs):
            cols = []
            for lo, hi in zip(lows, highs):
                if hi <= lo:
                    cols.append([int(round(lo))] * m); continue
                cut = [lo + (hi - lo) * (i + rng.random()) / m for i in range(m)]
                rng.shuffle(cut)
                cols.append([int(round(c)) for c in cut])
            return [[cols[j][i] for j in range(len(self.box))]
                    for i in range(m)]

        glo = [lo for _, lo, hi in self.box]
        ghi = [hi for _, lo, hi in self.box]
        if centre is None:
            return block(n, glo, ghi)
        half = n // 2
        lows, highs = [], []
        for c, (_, lo, hi) in zip(centre, self.box):
            w = frac * (hi - lo)
            lows.append(max(lo, c - w)); highs.append(min(hi, c + w))
        return block(half, lows, highs) + block(n - half, glo, ghi)

    def clamp(self, v):
        return [
            max(lo, min(hi, int(round(x))))
            for x, (_, lo, hi) in zip(v, self.box)
        ]

    def run(self, seed):
        if seed != self.a.seed:
            raise ValueError("sampling seed must match the recorded base seed")
        sampling_seed = derive_seed(
            seed, "funnel-search", self.k, "sampling")
        rng = random.Random(sampling_seed)
        t0 = time.time()
        print(
            f"# proxy_cost_model={PROXY_COST_MODEL} "
            "allow_unverified_cost_model=1 "
            f"proxy_ordering={PROXY_ORDERING_PROTOCOL} "
            f"search_box={SEARCH_BOX_PROTOCOL} "
            f"scale_centi={self.box[0][1]}:{self.box[0][2]}")

        # The counter stage uses failure first and predicted work second. It
        # only orders a finite candidate screen; the real codec still decides.
        pts = self.lhs(
            rng, self.a.screen, self.a.init, frac=self.a.init_frac)
        res = self.measure(pts, self.a.screen_cells)
        alive = [(o, f, p) for p, (o, f) in zip(pts, res)
                 if o is not None and f is not None]
        print(f"  screen  {self.a.screen} pts @{self.a.screen_cells} cells  "
              f"{time.time()-t0:5.1f}s  {len(alive)} measured")
        if not alive:
            print("  counter stage returned no candidates")
            return None
        alive.sort(key=proxy_order)
        best, best_fail, bestv = alive[0]

        t1, r, used = time.time(), 0.25, 0
        while r > 1e-3 and used < self.a.refine:
            cand = []
            for j, (_, lo, hi) in enumerate(self.box):
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
                if (o is not None and f is not None and
                        proxy_order((o, f, p)) <
                        proxy_order((best, best_fail, bestv))):
                    best, best_fail, bestv, moved = o, f, p, True
            if not moved:
                r *= 0.5
        print(f"  refine  {used} pts @{self.a.screen_cells} cells  "
              f"{time.time()-t1:5.1f}s  {alive[0][0]:,.0f} -> {best:,.0f}")

        # None means the true shipped arm: no peel override and no staircase
        # scale override. It is a mandatory finalist and control.
        incumbent = None
        pool = [bestv] + [v for _, _, v in alive if v != bestv]
        pool = pool[:self.a.finals - 1] + [incumbent]
        t2 = time.time()
        ok, dead = self.real_select(pool)
        print(f"  finals  {len(pool)} gated @bb{self.a.gate_bb}/"
              f"{self.a.gate_trials}tr, top {self.a.rank_top} timed @bb"
              f"{self.a.rank_bb}/{self.a.real_trials}tr  {time.time()-t2:5.1f}s  "
              f"{len(dead)} rejected as non-decoding")
        if not ok:
            print("  no candidate decoded -- the incumbent stands at this K")
            return None
        winner_goodput, winner_metrics, winner_vector = ok[0]
        shipped_control = next(
            (metrics for metrics, vector in self.rank_measurements
             if vector is None),
            None)
        if shipped_control is None:
            raise MeasurementError("mandatory shipped rank control is missing")
        if winner_vector is None:
            winner_coordinates = {
                "scale": -1.0, "p1": 100, "tilt": 0,
                "dmax": 64, "absorb": 100,
            }
            winner_pmf = list(self.native_profile.pmf)
            winner_mode = "shipped"
        else:
            winner_coordinates = {
                "scale": winner_vector[0] / 100.0,
                "p1": winner_vector[1],
                "tilt": winner_vector[2] - 100,
                "dmax": winner_vector[3],
                "absorb": winner_vector[4],
            }
            winner_pmf = peel_family(
                self.native_profile.pmf,
                winner_coordinates["p1"], winner_coordinates["tilt"],
                winner_coordinates["dmax"], winner_coordinates["absorb"])
            winner_mode = "trained"
        if winner_pmf is None:
            raise MeasurementError("winner has an invalid peel PMF")
        winner_receipt = {
            "schema": FUNNEL_RESULT_SCHEMA,
            "K": self.k,
            "mode": winner_mode,
            "coordinates": winner_coordinates,
            "peel_pmf": winner_pmf,
            "goodput": winner_goodput,
            "trials": self.a.real_trials,
            "block_bytes": self.a.rank_bb,
            "rejected": len(dead),
            "shipped_control": shipped_control.as_dict(),
            "shipped_goodput": shipped_control.goodput(self.k),
            **winner_metrics.as_dict(),
        }
        print(
            FUNNEL_RESULT_PREFIX +
            json.dumps(
                winner_receipt, sort_keys=True, separators=(",", ":"),
                allow_nan=False))
        print(f"\n  K={self.k}  structure {self.struct}  "
              f"total {time.time()-t0:.1f}s, {self.calls} bench calls\n")
        print(f"  {'rank':>4} {'goodput':>10} {'MB/s':>9} {'OH_mean':>9}  parameters")
        for i, (g, measurement, v) in enumerate(ok[:self.a.show], 1):
            metrics = (
                f"fail={measurement.fail} OH_sd={measurement.oh_sd:.17g} "
                f"OH50={measurement.oh50:.17g} "
                f"OH95={measurement.oh95:.17g} "
                f"OH99={measurement.oh99:.17g} "
                f"OH_max={measurement.oh_max:.17g} "
                f"seed={measurement.seed}")
            if v is None:
                print(f"  {i:>4} {g:>10.1f} "
                      f"{measurement.decode_mbps:>9.1f} "
                      f"{measurement.oh_mean:>9.4f}  "
                      "mode=shipped scale=-1.00 p1=100 tilt=0 "
                      f"dmax=64 absorb=100 {metrics}")
            else:
                print(f"  {i:>4} {g:>10.1f} "
                      f"{measurement.decode_mbps:>9.1f} "
                      f"{measurement.oh_mean:>9.4f}  "
                      f"mode=trained scale={v[0]/100:.2f} "
                      f"p1={v[1]} tilt={v[2]-100} dmax={v[3]} "
                      f"absorb={v[4]} {metrics}")
        print(f"\n  NOTE goodput = MB/s * K/(K+OH); failure is a hard gate, not a "
              f"term.\n       Top finalists are typically within noise of each other.")
        return ok

    def real_probe(self, v, trials, bb, tier):
        """Return full real-codec metrics under the tier's paired seed."""
        seed = derive_seed(
            self.a.seed, "funnel-search", self.k, tier, trials, bb)
        if v is None:
            return compare_probe(
                self.a.bench, self.k, trials, bb, seed=seed)
        w = peel_family(
            self.native_profile.pmf,
            v[1], v[2] - 100, v[3], v[4])
        if w is None:
            return None
        return compare_probe(
            self.a.bench, self.k, trials, bb,
            peel_weights=w, degree_scale=v[0] / 100.0, seed=seed)

    def real_select(self, pool):
        """Use a cheap recovery gate, then rank survivors at real payload."""
        self.rank_measurements = []
        # Deduplicate by the complete active configuration. Scale is an
        # independent staircase override, so two equal PMFs at different
        # scales are distinct arms. The true shipped arm has no override at all
        # and is represented only by None.
        unique_pool = []
        seen_configs = set()
        have_shipped = False
        for vector in pool:
            if vector is None:
                if not have_shipped:
                    unique_pool.append(None)
                    have_shipped = True
                continue
            candidate_pmf = peel_family(
                self.native_profile.pmf,
                vector[1], vector[2] - 100, vector[3], vector[4])
            if candidate_pmf is None:
                continue
            key = (vector[0], tuple(candidate_pmf))
            if key in seen_configs:
                continue
            seen_configs.add(key)
            unique_pool.append(vector)
        if not have_shipped:
            unique_pool.append(None)
        pool = unique_pool
        passed, dead = [], []
        for v in pool:
            r = self.real_probe(
                v, self.a.gate_trials, self.a.gate_bb, "gate")
            if r is None:
                continue
            (dead if r.fail > 0 else passed).append((r, v))
        out = []
        # The shipped arm is a mandatory control, not merely the last member of
        # a cost-sorted candidate list.  The old slice silently omitted it
        # whenever at least rank_top trained candidates passed the cheap gate.
        rank_vectors = [
            v for _, v in passed if v is not None
        ][:self.a.rank_top]
        if any(v is None for v in pool):
            rank_vectors.append(None)
        for v in rank_vectors:
            r = self.real_probe(
                v, self.a.real_trials, self.a.rank_bb, "rank")
            if r is None:
                continue
            self.rank_measurements.append((r, v))
            if r.fail > 0:                    # slipped through the cheap gate
                if not any(dead_v is v for _, dead_v in dead):
                    dead.append((r, v))
                continue
            out.append((r.goodput(self.k), r, v))
        return sorted(out, key=lambda x: (-x[0], x[2] is not None)), dead

    def stock_pmf(self):
        return list(self.native_profile.pmf)

    @staticmethod
    def family(stock, p1, tilt, dmax, absorb):
        return peel_family(stock, p1, tilt, dmax, absorb)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--K", type=int, required=True)
    ap.add_argument("--bench", default=BENCH_DEFAULT)
    ap.add_argument("--screen", type=int, default=3000)
    ap.add_argument("--refine", type=int, default=400)
    ap.add_argument("--finals", type=int, default=20)
    ap.add_argument("--screen-cells", type=int, default=1000)
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
    ap.add_argument(
        "--allow-unverified-cost-model", action="store_true",
        help="opt in to the embedded ES cost model whose raw calibration "
             "provenance is unavailable")
    a = ap.parse_args(argv)
    if (not 2 <= a.K <= 64000 or a.screen < 1 or
            a.refine < 0 or a.finals < 2 or a.show < 1 or
            a.screen_cells < 1 or a.threads < 1 or a.batch < 1 or
            a.cell_base < 0 or a.gate_trials < 1 or a.gate_bb < 1 or
            a.real_trials < 1 or a.rank_bb < 1 or a.rank_top < 1 or
            not 0.0 < a.init_frac <= 1.0 or
            not 0 <= a.seed <= 0xffffffffffffffff):
        ap.error("invalid K, budget, payload, fraction, or uint64 seed")
    if not a.allow_unverified_cost_model:
        print(
            f"REFUSED: proxy cost model {PROXY_COST_MODEL!r} has no replayable "
            "raw calibration; use --allow-unverified-cost-model only for an "
            "explicit experiment",
            file=sys.stderr)
        return 2
    if a.init:
        try:
            a.init = [int(x) for x in a.init.split(',')]
        except ValueError:
            ap.error("--init must contain five comma-separated integers")
        if len(a.init) != 1 + len(NON_SCALE_BOX):
            ap.error("--init must contain five comma-separated integers")
    try:
        f = Funnel(a)
        if a.init:
            a.init = f.clamp(a.init)
        ranked = f.run(a.seed)
    except (MeasurementError, OSError, ValueError) as error:
        print(f"REFUSED funnel measurement: {error}", file=sys.stderr)
        return 1
    return 0 if ranked else 1


if __name__ == "__main__":
    sys.exit(main())
