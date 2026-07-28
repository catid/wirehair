#!/usr/bin/env python3
"""Run the per-K funnel over the embedded proxy ladder.

The pivot result warm-starts both directions.  Child failures abort the entire
sweep, and a versioned table is atomically published only when every requested
K produced a complete real-codec receipt.
"""
import argparse
import math
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from peel_codec import (                                  # noqa: E402
    MeasurementError,
    PEELTIMING_EVICT_BYTES_MAX,
    PROXY_K_LADDER,
    _native_peel_cdf_signature,
    _validate_paired_measurement_receipt,
    capture_artifact_identity,
    derive_seed,
    family,
    make_paired_search_receipt,
    make_table_document,
    paired_context_thermal_source,
    paired_selected_metrics,
    require_distinct_from_source_provenance,
    require_distinct_paths,
    stock_profile,
    strict_json_loads,
    validate_paired_stock_control_receipt,
    validate_peeltiming_dimensions,
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


def _close_optional_log(logf, failure=None):
    """Best-effort failure note and close for a non-authoritative text log."""
    if logf is None:
        return []
    errors = []
    if failure is not None:
        try:
            logf.write(f"=== FAILED ===\n{failure}\n")
        except OSError as error:
            errors.append(f"could not append failure note: {error}")
    try:
        logf.close()
    except OSError as error:
        errors.append(f"could not close log: {error}")
    return errors


def _report_log_errors(errors):
    for error in errors:
        print(f"  WARNING: optional log {error}", file=sys.stderr)


def warm_start(best):
    """Return an offset-lattice neighbour seed, or none for shipped control."""
    if best.get("reverted_to_shipped") or best["scale"] < 0.0:
        return None
    return [
        int(round(best["scale"] * 100)), best["p1"],
        best["tilt"] + 100, best["dmax"], best["absorb"],
    ]


def run_one(
        bench, k, init, construction_seed_base, loss_seed_base,
        paired_replicates, allow_unverified_cost_model, *,
        target_profile, seed_policy, loss, schedule, paired_context,
        paired_warmups, paired_inner_reps, max_overhead, cache_state,
        evict_bytes, rank_margin):
    s, r, f, sc = budget(k)
    rt = paired_replicates
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
    cmd += ["--paired-replicates", str(rt), "--gate-trials", str(gt),
            "--gate-bb", str(gbb), "--rank-bb", str(rank_bb),
            "--rank-top", str(rtop),
            "--paired-context", paired_context,
            "--paired-warmups", str(paired_warmups),
            "--paired-inner-reps", str(paired_inner_reps),
            "--max-overhead", str(max_overhead),
            "--cache-state", cache_state,
            "--evict-bytes", str(evict_bytes),
            "--rank-margin", f"{rank_margin:.17g}"]
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
        "trials", "block_bytes", "rejected", "evaluated_coordinates",
        "evaluated_pmf", "paired_measurement", "stock_control_source",
        "stock_control_measurement", "shipped_goodput",
        "fail", "oh_mean", "OH_sd", "OH50", "OH95", "OH99", "OH_max",
        "solve_mbps",
    }
    coordinate_names = ("scale", "p1", "tilt", "dmax", "absorb")
    stock_coordinates = {
        "scale": -1.0, "p1": 100, "tilt": 0,
        "dmax": 64, "absorb": 100,
    }
    if (set(receipt) != expected_receipt_fields or
            receipt.get("schema") != FUNNEL_RESULT_SCHEMA or
            type(receipt.get("K")) is not int or receipt["K"] != k or
            receipt.get("mode") not in ("trained", "scale-only", "shipped") or
            type(receipt.get("trials")) is not int or
            receipt["trials"] != rt or
            type(receipt.get("block_bytes")) is not int or
            receipt["block_bytes"] != rank_bb or
            type(receipt.get("rejected")) is not int or
            not 0 <= receipt["rejected"] <= 2 * f or
            not isinstance(receipt.get("paired_measurement"), dict) or
            any(type(receipt.get(name)) is not float for name in (
                "goodput", "shipped_goodput", "oh_mean", "OH_sd", "OH50",
                "OH95", "OH99", "OH_max", "solve_mbps")) or
            type(receipt.get("fail")) is not int):
        raise MeasurementError(
            f"peel_funnel emitted invalid paired winner metadata for K={k}")

    def checked_coordinates(value, label):
        if (not isinstance(value, dict) or
                set(value) != set(coordinate_names) or
                type(value.get("scale")) is not float or
                any(type(value.get(name)) is not int
                    for name in ("p1", "tilt", "dmax", "absorb"))):
            raise MeasurementError(
                f"peel_funnel emitted invalid {label} for K={k}")
        scale = value["scale"]
        if (not math.isfinite(scale) or
                (scale != -1.0 and not 0.0 <= scale <= 64000.0) or
                not 0 <= value["p1"] <= 400 or
                not -100 <= value["tilt"] <= 1600 or
                not 2 <= value["dmax"] <= 64 or
                not 0 <= value["absorb"] <= 100):
            raise MeasurementError(
                f"peel_funnel emitted out-of-domain {label} for K={k}")
        if scale != -1.0:
            scale_centi = int(round(scale * 100.0))
            if (not box[0][1] <= scale_centi <= box[0][2] or
                    scale != scale_centi / 100.0):
                raise MeasurementError(
                    f"peel_funnel emitted off-lattice {label} for K={k}")
        return dict(value)

    coordinates = checked_coordinates(receipt["coordinates"], "coordinates")
    evaluated_coordinates = checked_coordinates(
        receipt["evaluated_coordinates"], "evaluated coordinates")
    for name in ("peel_pmf", "evaluated_pmf"):
        pmf = receipt.get(name)
        if (not isinstance(pmf, list) or
                any(type(probability) is not float for probability in pmf)):
            raise MeasurementError(
                f"peel_funnel emitted invalid {name} for K={k}")

    evaluated_shape_is_stock = (
        evaluated_coordinates["p1"] == 100 and
        evaluated_coordinates["tilt"] == 0 and
        evaluated_coordinates["dmax"] == 64 and
        evaluated_coordinates["absorb"] == 100
    )
    expected_evaluated_pmf = (
        list(profile.pmf) if evaluated_shape_is_stock else
        family(
            profile.pmf,
            evaluated_coordinates["p1"],
            evaluated_coordinates["tilt"],
            evaluated_coordinates["dmax"],
            evaluated_coordinates["absorb"],
        )
    )
    if (expected_evaluated_pmf is None or
            receipt["evaluated_pmf"] != expected_evaluated_pmf):
        raise MeasurementError(
            f"peel_funnel emitted an evaluated PMF mismatch for K={k}")
    evaluated_semantics_are_stock = (
        _native_peel_cdf_signature(receipt["evaluated_pmf"], k) ==
        _native_peel_cdf_signature(profile.pmf, k)
    )
    expected_construction_seed = derive_seed(
        construction_seed_base, "funnel-search", k, "rank", rt, rank_bb,
        "construction")
    expected_loss_seed = derive_seed(
        loss_seed_base, "funnel-search", k, "rank", rt, rank_bb, "loss")
    measurement = _validate_paired_measurement_receipt(
        receipt["paired_measurement"],
        block_count=k,
        block_bytes=rank_bb,
        candidate_pmf=receipt["evaluated_pmf"],
        degree_scale=(
            None if evaluated_coordinates["scale"] == -1.0
            else evaluated_coordinates["scale"]),
        native=profile.as_dict(),
        measurement_policy={
            "target_profile": target_profile,
            "seed_policy": seed_policy,
            "loss": float(loss),
            "schedule": schedule,
        },
        construction_seed=expected_construction_seed,
        loss_seed=expected_loss_seed,
        label=f"peel_funnel K={k} paired_measurement",
    )
    manifest = measurement.manifest
    if (manifest["warmup_replicates"] != paired_warmups or
            manifest["replicates"] != paired_replicates or
            manifest["inner_reps"] != paired_inner_reps or
            manifest["max_overhead"] != max_overhead or
            manifest["cache_state"] != cache_state or
            manifest["evict_bytes"] != evict_bytes or
            measurement.required_margin != rank_margin):
        raise MeasurementError(
            f"peel_funnel paired receipt contradicts invocation for K={k}")

    mode = receipt["mode"]
    control_measurement = validate_paired_stock_control_receipt(
        mode=mode,
        source=receipt["stock_control_source"],
        receipt=receipt["stock_control_measurement"],
        selected_measurement=measurement,
        block_count=k,
        block_bytes=rank_bb,
        native=profile.as_dict(),
        measurement_policy={
            "target_profile": target_profile,
            "seed_policy": seed_policy,
            "loss": float(loss),
            "schedule": schedule,
        },
        construction_seed=expected_construction_seed,
        loss_seed=expected_loss_seed,
        label=f"peel_funnel K={k} stock control",
    )
    if mode == "shipped":
        if (coordinates != stock_coordinates or
                evaluated_coordinates != stock_coordinates or
                receipt["peel_pmf"] != list(profile.pmf) or
                receipt["evaluated_pmf"] != list(profile.pmf)):
            raise MeasurementError(
                f"peel_funnel shipped fallback is not its predeclared "
                f"stock control for K={k}")
        selected = measurement.identity
    else:
        expected_mode = (
            "scale-only"
            if evaluated_semantics_are_stock
            else "trained")
        if (mode != expected_mode or coordinates != evaluated_coordinates or
                receipt["peel_pmf"] != receipt["evaluated_pmf"] or
                (expected_mode == "scale-only" and
                 not evaluated_shape_is_stock) or
                not measurement.valid_for_promotion):
            raise MeasurementError(
                f"peel_funnel selected an ineligible candidate for K={k}")
        selected = measurement.candidate
    selected_metrics = selected.as_dict()
    if (any(receipt[name] != selected_metrics[name]
            for name in selected_metrics) or
            receipt["goodput"] != selected.goodput(k) or
            receipt["shipped_goodput"] != measurement.identity.goodput(k)):
        raise MeasurementError(
            f"peel_funnel emitted stale paired summaries for K={k}")
    best = {
        **coordinates,
        **selected_metrics,
        "goodput": selected.goodput(k),
        "peel_pmf": list(receipt["peel_pmf"]),
        **({"reverted_to_shipped": True} if mode == "shipped" else {}),
    }
    best.update(K=k, S=profile.staircase,
                source_hits=profile.source_hits,
                target_mean=profile.target_mean,
                native=profile.as_dict(),
                search_receipt=make_paired_search_receipt(
                    measurement,
                    mode=mode,
                    block_count=k,
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
                    evaluated_coordinates=evaluated_coordinates,
                    evaluated_pmf=receipt["evaluated_pmf"],
                    stock_control_measurement=(
                        None if mode == "shipped" else
                        control_measurement),
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
                paired_replicates=rt)
    return best, p.stdout, p.stderr


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", default="build-fast/codec/wirehair_v2_bench")
    ap.add_argument("--out", default="tools/peel_table.json")
    ap.add_argument("--target-profile", required=True, choices=["dispatch-v1"])
    ap.add_argument("--seed-policy", required=True, choices=["raw"])
    ap.add_argument(
        "--schedule", required=True,
        choices=[
            "iid", "burst", "permutation", "systematic-first",
            "repair-only", "adversarial",
        ])
    ap.add_argument("--loss", type=float, required=True)
    ap.add_argument("--construction-seed", type=int, required=True)
    ap.add_argument("--loss-seed", type=int, required=True)
    ap.add_argument(
        "--paired-replicates", type=int, default=16,
        help="even measured replicate count for each final paired timing run")
    ap.add_argument("--paired-context", required=True)
    ap.add_argument("--paired-warmups", type=int, default=2)
    ap.add_argument("--paired-inner-reps", type=int, default=1)
    ap.add_argument("--max-overhead", type=int, default=512)
    ap.add_argument("--cache-state", choices=["warm", "cold"], default="warm")
    ap.add_argument("--evict-bytes", type=int, default=64 * 1024 * 1024)
    ap.add_argument("--rank-margin", type=float, default=0.0)
    ap.add_argument("--kmax", type=int, default=64000)
    ap.add_argument("--log", default=None)
    ap.add_argument(
        "--allow-unverified-cost-model", action="store_true",
        help="explicitly opt in to the proxy whose raw calibration is missing")
    a = ap.parse_args(argv)
    if (a.paired_replicates < 4 or a.paired_replicates % 2 != 0 or
            a.paired_warmups < 0 or a.paired_warmups % 2 != 0 or
            not 1 <= a.paired_inner_reps <= 1024 or
            a.max_overhead < 0 or
            not 4096 <= a.evict_bytes <= PEELTIMING_EVICT_BYTES_MAX or
            (a.cache_state == "cold" and a.paired_inner_reps != 1) or
            not 0.0 <= a.rank_margin <= 1.0 or
            not 2 <= a.kmax <= 64000 or
            not valid_loss_rate(a.loss) or
            not 0 <= a.construction_seed <= 0xffffffffffffffff or
            not 0 <= a.loss_seed <= 0xffffffffffffffff):
        ap.error("invalid paired settings, K, loss, or uint64 seed")
    if not a.allow_unverified_cost_model:
        print(
            f"REFUSED: proxy cost model {PROXY_COST_MODEL!r} is unverified; "
            "pass --allow-unverified-cost-model only for an explicit "
            "experiment",
            file=sys.stderr)
        return 2

    ks = [k for k in k_list() if k <= a.kmax]
    if not ks:
        print("REFUSED: --kmax selects no proxy-table K values", file=sys.stderr)
        return 2
    try:
        for k in ks:
            validate_peeltiming_dimensions(
                block_count=k,
                block_bytes=4096,
                target_profile=a.target_profile,
                seed_policy=a.seed_policy,
                construction_seed=a.construction_seed,
                loss=a.loss,
                loss_seed=a.loss_seed,
                schedule=a.schedule,
                warmup_replicates=a.paired_warmups,
                replicates=a.paired_replicates,
                inner_reps=a.paired_inner_reps,
                max_overhead=a.max_overhead,
                cache_state=a.cache_state,
                evict_bytes=a.evict_bytes,
                required_margin=a.rank_margin,
            )
    except ValueError as error:
        ap.error(str(error))
    try:
        require_distinct_paths(
            a.out, a.bench, "--out", "the benchmark executable")
        require_distinct_paths(
            a.out, a.paired_context, "--out", "--paired-context")
        require_distinct_from_source_provenance(
            a.out, "--out", "tools/peel_sweep.py")
        if a.log is not None:
            require_distinct_paths(
                a.log, a.bench, "--log", "the benchmark executable")
            require_distinct_paths(
                a.log, a.paired_context, "--log", "--paired-context")
            require_distinct_paths(a.log, a.out, "--log", "--out")
            require_distinct_from_source_provenance(
                a.log, "--log", "tools/peel_sweep.py")
        thermal_source = paired_context_thermal_source(a.paired_context)
        require_distinct_paths(
            a.out, thermal_source, "--out", "the live thermal CSV")
        if a.log is not None:
            require_distinct_paths(
                a.log, thermal_source, "--log", "the live thermal CSV")
    except MeasurementError as error:
        print(f"REFUSED output: {error}", file=sys.stderr)
        return 2
    try:
        identity = capture_artifact_identity(
            a.bench, "tools/peel_sweep.py")
    except MeasurementError as error:
        print(f"REFUSED measurement: {error}", file=sys.stderr)
        return 2
    pivot = 128 if 128 in ks else ks[len(ks) // 2]
    up = [k for k in ks if k >= pivot]
    down = [k for k in ks if k < pivot][::-1]

    try:
        if a.log is not None:
            require_distinct_from_source_provenance(
                a.log, "--log", "tools/peel_sweep.py")
        logf = open(a.log, "w") if a.log else None
    except (MeasurementError, OSError) as error:
        print(f"REFUSED log output: {error}", file=sys.stderr)
        return 2
    table = {}

    def walk(seq, label, init=None):
        for k in seq:
            best, out, err = run_one(
                a.bench, k, init, a.construction_seed, a.loss_seed,
                a.paired_replicates, a.allow_unverified_cost_model,
                target_profile=a.target_profile,
                seed_policy=a.seed_policy,
                loss=a.loss,
                schedule=a.schedule,
                paired_context=a.paired_context,
                paired_warmups=a.paired_warmups,
                paired_inner_reps=a.paired_inner_reps,
                max_overhead=a.max_overhead,
                cache_state=a.cache_state,
                evict_bytes=a.evict_bytes,
                rank_margin=a.rank_margin)
            table[k] = best
            # warm start the NEXT K from this one, on the offset lattice
            init = warm_start(best)
            print(f"  K={k:<6} tilt={best['tilt']:>+5} p1={best['p1']:>3} "
                  f"dmax={best['dmax']:>2} absorb={best['absorb']:>3} "
                  f"scale={best['scale']:>7.2f} "
                  f"{best['solve_mbps']:>7.1f}MB/s "
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
        _report_log_errors(_close_optional_log(logf, error))
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
                "paired_replicates": a.paired_replicates,
                "paired_context": os.path.realpath(a.paired_context),
                "paired_warmups": a.paired_warmups,
                "paired_inner_reps": a.paired_inner_reps,
                "max_overhead": a.max_overhead,
                "cache_state": a.cache_state,
                "evict_bytes": a.evict_bytes,
                "rank_margin": a.rank_margin,
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
        require_distinct_paths(
            a.out, a.bench, "--out", "the benchmark executable")
        require_distinct_paths(
            a.out, a.paired_context, "--out", "--paired-context")
        if paired_context_thermal_source(a.paired_context) != thermal_source:
            raise MeasurementError(
                "--paired-context changed its thermal source during sweep")
        require_distinct_paths(
            a.out, thermal_source, "--out", "the live thermal CSV")
        if a.log is not None:
            require_distinct_paths(
                a.log, thermal_source, "--log", "the live thermal CSV")
            require_distinct_paths(a.log, a.out, "--log", "--out")
        require_distinct_from_source_provenance(
            a.out, "--out", "tools/peel_sweep.py")
        write_json_atomic(a.out, document)
    except (MeasurementError, OSError, ValueError) as error:
        print(f"  REFUSED publication: {error}", file=sys.stderr)
        _report_log_errors(_close_optional_log(logf))
        return 1
    print(f"\n  {len(table)} K values written to {a.out}")
    _report_log_errors(_close_optional_log(logf))
    return 0


if __name__ == "__main__":
    sys.exit(main())
