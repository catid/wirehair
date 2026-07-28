#!/usr/bin/env python3
"""Run the per-K funnel over the embedded proxy ladder.

The pivot result warm-starts both directions.  Child failures abort the entire
sweep, and a versioned table is atomically published only when every requested
K produced a complete real-codec receipt.
"""
import argparse
import math
import subprocess
import sys
import time

HERE = __import__("os").path.dirname(__import__("os").path.abspath(__file__))
sys.path.insert(0, HERE)
from peel_codec import (                                  # noqa: E402
    MeasurementError,
    PROXY_K_LADDER,
    RecoveryMetrics,
    capture_artifact_identity,
    derive_seed,
    family,
    make_table_document,
    make_search_receipt,
    stock_profile,
    strict_json_loads,
    valid_loss_rate,
    write_json_atomic,
)
from peel_funnel import (                                 # noqa: E402
    FUNNEL_RESULT_PREFIX,
    FUNNEL_RESULT_SCHEMA,
    PROXY_COST_MODEL,
    PROXY_MEASURE_REGIME,
    PROXY_ORDERING_PROTOCOL,
    SEARCH_BOX_PROTOCOL,
    search_box,
)


# The block counts represented by the embedded proxy table.  Its raw
# calibration provenance is unavailable, so use still requires an explicit
# opt-in and every finalist is measured on the real codec.
PROXY_K = list(PROXY_K_LADDER)


def k_list():
    return list(PROXY_K)


def budget(k):
    """(screen, refine, finals, screen_cells) sized for K."""
    if k <= 1024:
        return 3000, 400, 16, 1000
    if k <= 8192:
        return 1200, 250, 12, 1000
    if k <= 32768:
        return 500, 150, 8, 600
    return 250, 100, 6, 400


def gate_cfg(k):
    """(gate_trials, gate_bb, rank_top) for the cheap solve-rate tier."""
    return (25, 64, 3)


def real_trials_for(k):
    """Trials for the real-codec ranking stage."""
    if k <= 512:   return 800
    if k <= 4096:  return 400
    if k <= 16384: return 200
    return 100


def warm_start(best):
    """Return an offset-lattice neighbour seed, or none for shipped control."""
    if best.get("reverted_to_shipped") or best["scale"] < 0.0:
        return None
    return [
        int(round(best["scale"] * 100)), best["p1"],
        best["tilt"] + 100, best["dmax"], best["absorb"],
    ]


def run_one(
        bench, k, init, construction_seed_base, loss_seed_base, real_trials,
        allow_unverified_cost_model, *, target_profile, seed_policy, loss,
        schedule):
    s, r, f, sc = budget(k)
    rt = real_trials or real_trials_for(k)
    rank_bb = 4096
    threads = 64
    batch = 60
    cell_base = 900_000_000
    profile = stock_profile(bench, k, target_profile=target_profile)
    box = search_box(profile)
    if init is not None:
        if (not isinstance(init, list) or len(init) != len(box) or
                any(isinstance(value, bool) or
                    not isinstance(value, (int, float)) or
                    not math.isfinite(value)
                    for value in init)):
            raise MeasurementError(
                f"invalid warm start for K={k}")
        init = [
            max(lo, min(hi, int(round(value))))
            for value, (_, lo, hi) in zip(init, box)
        ]
    cmd = [sys.executable, f"{HERE}/peel_funnel.py", "--K", str(k), "--bench", bench,
           "--screen", str(s), "--refine", str(r), "--finals", str(f),
           "--screen-cells", str(sc),
           "--threads", str(threads), "--batch", str(batch),
           "--cell-base", str(cell_base),
           "--target-profile", target_profile,
           "--seed-policy", seed_policy,
           "--loss", f"{loss:.17g}", "--schedule", schedule,
           "--construction-seed", str(construction_seed_base),
           "--loss-seed", str(loss_seed_base), "--show", "1"]
    if allow_unverified_cost_model:
        cmd.append("--allow-unverified-cost-model")
    if init is not None:
        cmd += ["--init", ",".join(str(x) for x in init)]
    gt, gbb, rtop = gate_cfg(k)
    cmd += ["--real-trials", str(rt), "--gate-trials", str(gt),
            "--gate-bb", str(gbb), "--rank-bb", str(rank_bb),
            "--rank-top", str(rtop)]
    t0 = time.time()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True)
    except OSError as error:
        raise MeasurementError(
            f"could not run peel_funnel for K={k}: {error}")
    el = time.time() - t0
    if p.returncode != 0:
        detail = p.stderr.strip() or p.stdout.strip() or "no output"
        raise MeasurementError(
            f"peel_funnel failed for K={k} with exit {p.returncode}: {detail}")
    receipt = None
    for line in p.stdout.splitlines():
        if line.startswith(FUNNEL_RESULT_PREFIX):
            if receipt is not None:
                raise MeasurementError(
                    f"peel_funnel emitted duplicate winner receipts for K={k}")
            try:
                receipt = strict_json_loads(
                    line[len(FUNNEL_RESULT_PREFIX):])
            except (UnicodeDecodeError, TypeError, ValueError) as error:
                raise MeasurementError(
                    f"peel_funnel emitted malformed JSON for K={k}: {error}")
    if not isinstance(receipt, dict):
        raise MeasurementError(
            f"peel_funnel emitted no complete winner receipt for K={k}")
    expected_receipt_fields = {
        "schema", "K", "mode", "coordinates", "peel_pmf", "goodput",
        "trials", "block_bytes", "rejected", "shipped_control",
        "shipped_goodput", "construction_seed", "loss_seed",
        "target_receipt", "fail", "oh_mean", "OH_sd", "OH50", "OH95",
        "OH99", "OH_max", "decode_mbps",
    }
    if (set(receipt) != expected_receipt_fields or
            receipt.get("schema") != FUNNEL_RESULT_SCHEMA or
            receipt.get("K") != k or
            type(receipt.get("K")) is not int or
            receipt.get("mode") not in
            ("trained", "scale-only", "shipped") or
            receipt.get("trials") != rt or
            type(receipt.get("trials")) is not int or
            receipt.get("block_bytes") != rank_bb or
            type(receipt.get("block_bytes")) is not int or
            type(receipt.get("rejected")) is not int or
            # Each proxy finalist can add one canonical scale-only gate arm.
            not 0 <= receipt["rejected"] <= 2 * f or
            not isinstance(receipt.get("coordinates"), dict) or
            not isinstance(receipt.get("peel_pmf"), list) or
            not isinstance(receipt.get("shipped_control"), dict)):
        raise MeasurementError(
            f"peel_funnel emitted invalid winner metadata for K={k}")
    coordinates = receipt["coordinates"]
    coordinate_names = ("scale", "p1", "tilt", "dmax", "absorb")
    if (set(coordinates) != set(coordinate_names) or
            type(coordinates["scale"]) is not float or
            any(type(coordinates[name]) is not int
                for name in ("p1", "tilt", "dmax", "absorb")) or
            any(
                type(receipt.get(name)) is not int
                for name in ("construction_seed", "loss_seed")) or
            not isinstance(receipt.get("target_receipt"), dict) or
            type(receipt.get("fail")) is not int or
            any(type(receipt.get(name)) is not float
                for name in (
                    "goodput", "decode_mbps", "oh_mean", "OH_sd", "OH50",
                    "OH95", "OH99", "OH_max", "shipped_goodput"))):
        raise MeasurementError(
            f"peel_funnel emitted invalid typed metrics for K={k}")
    try:
        best = {
            "goodput": float(receipt["goodput"]),
            "decode_mbps": float(receipt["decode_mbps"]),
            "oh_mean": float(receipt["oh_mean"]),
            "scale": float(coordinates["scale"]),
            "p1": int(coordinates["p1"]),
            "tilt": int(coordinates["tilt"]),
            "dmax": int(coordinates["dmax"]),
            "absorb": int(coordinates["absorb"]),
            "fail": int(receipt["fail"]),
            "OH_sd": float(receipt["OH_sd"]),
            "OH50": float(receipt["OH50"]),
            "OH95": float(receipt["OH95"]),
            "OH99": float(receipt["OH99"]),
            "OH_max": float(receipt["OH_max"]),
            "construction_seed": int(receipt["construction_seed"]),
            "loss_seed": int(receipt["loss_seed"]),
            "target_receipt": dict(receipt["target_receipt"]),
            **({"reverted_to_shipped": True}
               if receipt["mode"] == "shipped" else {}),
            "peel_pmf": receipt["peel_pmf"],
        }
    except (KeyError, TypeError, ValueError, OverflowError) as error:
        raise MeasurementError(
            f"peel_funnel emitted incomplete winner metrics for K={k}: "
            f"{error}")
    if (
            any(best[name] != coordinates[name] for name in coordinate_names)):
        raise MeasurementError(
            f"peel_funnel emitted invalid winner coordinates for K={k}")
    if (receipt["mode"] == "shipped" and
            coordinates != {
                "scale": -1.0, "p1": 100, "tilt": 0,
                "dmax": 64, "absorb": 100,
            }):
        raise MeasurementError(
            f"peel_funnel emitted non-production shipped coordinates for K={k}")
    if receipt["mode"] in ("trained", "scale-only") and best["scale"] < 0.0:
        raise MeasurementError(
            f"peel_funnel emitted an unset scale for a trained arm at K={k}")
    finite_fields = (
        "goodput", "decode_mbps", "oh_mean", "OH_sd", "OH50",
        "OH95", "OH99", "OH_max")
    if (any(not math.isfinite(best[name]) or best[name] < 0.0
            for name in finite_fields) or
            not math.isfinite(best["scale"]) or
            (best["scale"] != -1.0 and
             not 0.0 <= best["scale"] <= 64000.0) or
            not 0 <= best["p1"] <= 400 or
            not -100 <= best["tilt"] <= 1600 or
            not 2 <= best["dmax"] <= 64 or
            not 0 <= best["absorb"] <= 100 or
            best["fail"] != 0 or
            not 0 <= best["construction_seed"] <= 0xffffffffffffffff or
            not 0 <= best["loss_seed"] <= 0xffffffffffffffff or
            best["oh_mean"] > best["OH_max"] or
            not (best["OH50"] <= best["OH95"] <= best["OH99"] <=
                 best["OH_max"])):
        raise MeasurementError(
            f"peel_funnel emitted invalid recovery metrics for K={k}")
    if receipt["mode"] in ("trained", "scale-only"):
        scale_centi = int(round(best["scale"] * 100.0))
        if (not box[0][1] <= scale_centi <= box[0][2] or
                best["scale"] != scale_centi / 100.0):
            raise MeasurementError(
                f"peel_funnel emitted an off-lattice scale for K={k}")
    expected_construction_seed = derive_seed(
        construction_seed_base, "funnel-search", k, "rank", rt, rank_bb,
        "construction")
    expected_loss_seed = derive_seed(
        loss_seed_base, "funnel-search", k, "rank", rt, rank_bb, "loss")
    if (best["construction_seed"] != expected_construction_seed or
            best["loss_seed"] != expected_loss_seed):
        raise MeasurementError(
            f"peel_funnel emitted invalid rank construction/loss seeds "
            f"for K={k}")
    expected_pmf = (
        list(profile.pmf)
        if receipt["mode"] in ("shipped", "scale-only") else
        family(
            profile.pmf, best["p1"], best["tilt"],
            best["dmax"], best["absorb"])
    )
    if (any(type(probability) is not float
            for probability in best["peel_pmf"]) or
            expected_pmf is None or best["peel_pmf"] != expected_pmf):
        raise MeasurementError(
            f"peel_funnel emitted a PMF/coordinate mismatch for K={k}")
    measurement = RecoveryMetrics(
        construction_seed=best["construction_seed"],
        loss_seed=best["loss_seed"],
        target_receipt=best["target_receipt"],
        fail=best["fail"],
        oh_mean=best["oh_mean"],
        oh_sd=best["OH_sd"],
        oh50=best["OH50"],
        oh95=best["OH95"],
        oh99=best["OH99"],
        oh_max=best["OH_max"],
        decode_mbps=best["decode_mbps"],
    )
    candidate_goodput = measurement.goodput(k)
    if best["goodput"] != candidate_goodput:
        raise MeasurementError(
            f"peel_funnel emitted inconsistent goodput for K={k}")
    shipped_raw = receipt["shipped_control"]
    expected_metric_fields = {
        "construction_seed", "loss_seed", "target_receipt", "fail",
        "oh_mean", "OH_sd", "OH50", "OH95", "OH99", "OH_max",
        "decode_mbps",
    }
    if (set(shipped_raw) != expected_metric_fields or
            any(
                type(shipped_raw.get(name)) is not int
                for name in ("construction_seed", "loss_seed")) or
            not isinstance(shipped_raw.get("target_receipt"), dict) or
            type(shipped_raw.get("fail")) is not int or
            any(type(shipped_raw.get(name)) is not float
                for name in (
                    "oh_mean", "OH_sd", "OH50", "OH95", "OH99",
                    "OH_max", "decode_mbps"))):
        raise MeasurementError(
            f"peel_funnel emitted invalid shipped-control types for K={k}")
    try:
        shipped_control = RecoveryMetrics(
            construction_seed=int(shipped_raw["construction_seed"]),
            loss_seed=int(shipped_raw["loss_seed"]),
            target_receipt=dict(shipped_raw["target_receipt"]),
            fail=int(shipped_raw["fail"]),
            oh_mean=float(shipped_raw["oh_mean"]),
            oh_sd=float(shipped_raw["OH_sd"]),
            oh50=float(shipped_raw["OH50"]),
            oh95=float(shipped_raw["OH95"]),
            oh99=float(shipped_raw["OH99"]),
            oh_max=float(shipped_raw["OH_max"]),
            decode_mbps=float(shipped_raw["decode_mbps"]),
        )
        shipped_goodput = receipt["shipped_goodput"]
    except (KeyError, TypeError, ValueError, OverflowError) as error:
        raise MeasurementError(
            f"peel_funnel emitted incomplete shipped-control metrics "
            f"for K={k}: {error}")
    shipped_values = (
        shipped_control.oh_mean, shipped_control.oh_sd,
        shipped_control.oh50, shipped_control.oh95,
        shipped_control.oh99, shipped_control.oh_max,
        shipped_control.decode_mbps,
    )
    canonical_shipped_goodput = shipped_control.goodput(k)
    if (shipped_control.construction_seed != expected_construction_seed or
            shipped_control.loss_seed != expected_loss_seed or
            not 0 <= shipped_control.fail <= rt or
            any(not math.isfinite(value) or value < 0.0
                for value in shipped_values) or
            shipped_control.oh_mean > shipped_control.oh_max or
            not (shipped_control.oh50 <= shipped_control.oh95 <=
                 shipped_control.oh99 <= shipped_control.oh_max) or
            not math.isfinite(shipped_goodput) or
            shipped_goodput != canonical_shipped_goodput):
        raise MeasurementError(
            f"peel_funnel emitted invalid shipped-control metrics for K={k}")
    if (receipt["mode"] in ("trained", "scale-only") and
            not candidate_goodput > canonical_shipped_goodput):
        raise MeasurementError(
            f"peel_funnel selected a candidate that did not beat shipped "
            f"for K={k}")
    if (receipt["mode"] == "shipped" and
            (measurement.as_dict() != shipped_control.as_dict() or
             candidate_goodput != canonical_shipped_goodput)):
        raise MeasurementError(
            f"peel_funnel shipped winner contradicts its control for K={k}")
    best["goodput"] = candidate_goodput
    shipped_goodput = canonical_shipped_goodput
    mode = receipt["mode"]
    best.update(K=k, S=profile.staircase,
                source_hits=profile.source_hits,
                target_mean=profile.target_mean,
                native=profile.as_dict(),
                search_receipt=make_search_receipt(
                    measurement,
                    mode=mode,
                    goodput=best["goodput"],
                    trials=rt,
                    block_bytes=rank_bb,
                    search_kind="unverified-proxy-funnel",
                    construction_seed_base=construction_seed_base,
                    loss_seed_base=loss_seed_base,
                    seed_domain="funnel-search",
                    coordinates={
                        name: best[name]
                        for name in (
                            "scale", "p1", "tilt", "dmax", "absorb")
                    },
                    peel_pmf=best["peel_pmf"],
                    shipped_control=shipped_control,
                    shipped_goodput=shipped_goodput,
                    context={
                        "proxy_cost_model": PROXY_COST_MODEL,
                        "proxy_measure_regime": dict(PROXY_MEASURE_REGIME),
                        "proxy_ordering": PROXY_ORDERING_PROTOCOL,
                        "search_box": SEARCH_BOX_PROTOCOL,
                        "scale_centi": [box[0][1], box[0][2]],
                        "warm_start": list(init) if init is not None else None,
                        "sampling_seed": derive_seed(
                            construction_seed_base,
                            "funnel-search", k, "sampling"),
                        "screen": s,
                        "refine": r,
                        "finals": f,
                        "screen_cells": sc,
                        "gate_trials": gt,
                        "gate_block_bytes": gbb,
                        "rank_top": rtop,
                        "threads": threads,
                        "batch": batch,
                        "cell_base": cell_base,
                    },
                ),
                seconds=round(el, 1),
                screen=s, screen_cells=sc, finals=f,
                rejected=receipt["rejected"],
                real_trials=rt)
    return best, p.stdout, p.stderr


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", default="build-fast/codec/wirehair_v2_bench")
    ap.add_argument("--out", default="tools/peel_table.json")
    ap.add_argument("--target-profile", required=True, choices=["dispatch-v1"])
    ap.add_argument("--seed-policy", required=True, choices=["raw"])
    ap.add_argument("--schedule", required=True, choices=["iid"])
    ap.add_argument("--loss", type=float, required=True)
    ap.add_argument("--construction-seed", type=int, required=True)
    ap.add_argument("--loss-seed", type=int, required=True)
    ap.add_argument("--real-trials", type=int, default=0,
                    help="0 = scale with K via real_trials_for()")
    ap.add_argument("--kmax", type=int, default=64000)
    ap.add_argument("--log", default=None)
    ap.add_argument(
        "--allow-unverified-cost-model", action="store_true",
        help="explicitly opt in to the proxy whose raw calibration is missing")
    a = ap.parse_args(argv)
    if (a.real_trials < 0 or not 2 <= a.kmax <= 64000 or
            not valid_loss_rate(a.loss) or
            not 0 <= a.construction_seed <= 0xffffffffffffffff or
            not 0 <= a.loss_seed <= 0xffffffffffffffff):
        ap.error("invalid trial count, K, loss, or uint64 seed")
    if not a.allow_unverified_cost_model:
        print(
            f"REFUSED: proxy cost model {PROXY_COST_MODEL!r} is unverified; "
            "pass --allow-unverified-cost-model only for an explicit "
            "experiment",
            file=sys.stderr)
        return 2

    try:
        identity = capture_artifact_identity(
            a.bench, "tools/peel_sweep.py")
    except MeasurementError as error:
        print(f"REFUSED measurement: {error}", file=sys.stderr)
        return 2
    ks = [k for k in k_list() if k <= a.kmax]
    if not ks:
        print("REFUSED: --kmax selects no proxy-table K values", file=sys.stderr)
        return 2
    pivot = 128 if 128 in ks else ks[len(ks) // 2]
    up = [k for k in ks if k >= pivot]
    down = [k for k in ks if k < pivot][::-1]

    try:
        logf = open(a.log, "w") if a.log else None
    except OSError as error:
        print(f"REFUSED log output: {error}", file=sys.stderr)
        return 2
    table = {}

    def walk(seq, label, init=None):
        for k in seq:
            best, out, err = run_one(
                a.bench, k, init, a.construction_seed, a.loss_seed,
                a.real_trials, a.allow_unverified_cost_model,
                target_profile=a.target_profile,
                seed_policy=a.seed_policy,
                loss=a.loss,
                schedule=a.schedule)
            table[k] = best
            # warm start the NEXT K from this one, on the offset lattice
            init = warm_start(best)
            print(f"  K={k:<6} tilt={best['tilt']:>+5} p1={best['p1']:>3} "
                  f"dmax={best['dmax']:>2} absorb={best['absorb']:>3} "
                  f"scale={best['scale']:>7.2f} "
                  f"{best['decode_mbps']:>7.1f}MB/s "
                  f"OH={best['oh_mean']:.4f} ({best['seconds']}s)", flush=True)
            if logf:
                logf.write(f"=== K={k} ({label}) ===\n{out}\n{err}\n")
        return init
    try:
        print(f"  walking UP from K={pivot} ({len(up)} values)", flush=True)
        walk(up, "up")
        pivot_entry = table[pivot]
        pivot_init = warm_start(pivot_entry)
        print(f"  walking DOWN from K={pivot} ({len(down)} values)", flush=True)
        walk(down, "down", pivot_init)
    except (MeasurementError, OSError, ValueError) as error:
        print(f"  REFUSED publication: {error}", file=sys.stderr)
        if logf:
            logf.write(f"=== FAILED ===\n{error}\n")
            logf.close()
        return 1
    try:
        document = make_table_document(
            table,
            generator="tools/peel_sweep.py",
            bench=a.bench,
            construction_seed_base=a.construction_seed,
            loss_seed_base=a.loss_seed,
            target_profile=a.target_profile,
            seed_policy=a.seed_policy,
            loss=a.loss,
            schedule=a.schedule,
            settings={
                "proxy_k_ladder": ks,
                "real_trials_override": a.real_trials,
                "target_profile": a.target_profile,
                "seed_policy": a.seed_policy,
                "loss": a.loss,
                "schedule": a.schedule,
                "proxy_cost_model": PROXY_COST_MODEL,
                "proxy_measure_regime": dict(PROXY_MEASURE_REGIME),
                "proxy_ordering": PROXY_ORDERING_PROTOCOL,
                "allow_unverified_cost_model": True,
                "search_box": SEARCH_BOX_PROTOCOL,
            },
            artifact_identity=identity,
        )
        write_json_atomic(a.out, document)
    except (MeasurementError, OSError, ValueError) as error:
        print(f"  REFUSED publication: {error}", file=sys.stderr)
        if logf:
            logf.close()
        return 1
    print(f"\n  {len(table)} K values written to {a.out}")
    if logf: logf.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
