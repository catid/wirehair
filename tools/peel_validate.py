#!/usr/bin/env python3
"""Validate a versioned peel table with one native paired experiment per K.

The selected search configuration and an explicit identity-PMF/identity-scale
control run through the same hook path in one native process with common
traces, counterbalanced order, and an identity A/A panel.  A separate untimed
semantic proof establishes equivalence to the true no-hook shipped path.
Candidate promotion depends only on the paired solve-cost confidence bound
plus recovery non-regression.  Overhead is retained as a diagnostic and never
becomes an implicit selection rule.
"""
import argparse
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from peel_codec import (                                 # noqa: E402
    MeasurementError,
    NATIVE_PAIRED_PROTOCOL,
    PEELTIMING_EVICT_BYTES_MAX,
    capture_artifact_identity,
    derive_seed,
    make_table_document,
    paired_context_thermal_source,
    paired_probe,
    read_table_document_snapshot,
    require_distinct_from_source_provenance,
    require_distinct_paths,
    require_paired_stock_control,
    stock_profile,
    validate_peeltiming_dimensions,
    valid_loss_rate,
    verify_benchmark_identity,
    write_json_atomic,
)


def require_distinct_output(source_path, output_path):
    """Fail closed unless an output path is demonstrably distinct."""
    require_distinct_paths(
        output_path, source_path, "--out", "the source table")


def require_safe_output(
        source_path, output_path, bench, paired_context,
        expected_thermal_source=None):
    """Protect every persistent input that validation must not overwrite."""
    require_distinct_output(source_path, output_path)
    require_distinct_paths(
        output_path, bench, "--out", "the benchmark executable")
    require_distinct_paths(
        output_path, paired_context, "--out", "--paired-context")
    require_distinct_from_source_provenance(
        output_path, "--out", "tools/peel_validate.py")
    thermal_source = paired_context_thermal_source(paired_context)
    if (expected_thermal_source is not None and
            thermal_source != expected_thermal_source):
        raise MeasurementError(
            "--paired-context changed its thermal source during validation")
    require_distinct_paths(
        output_path, thermal_source, "--out", "the live thermal CSV")
    return thermal_source


def canonical_percent(value):
    """Round a displayed percentage without emitting a signed-zero alias."""
    rounded = float(round(value, 2))
    return 0.0 if rounded == 0.0 else rounded


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", default="build-fast/codec/wirehair_v2_bench")
    ap.add_argument("--table", required=True)
    ap.add_argument("--out", default=None, help="write the corrected table here")
    ap.add_argument(
        "--paired-replicates", type=int, default=32,
        help="even measured replicate count for each validation experiment")
    ap.add_argument("--bb", type=int, required=True)
    ap.add_argument("--kmax", type=int, default=10 ** 9)
    ap.add_argument(
        "--margin", type=float, default=2.0,
        help="required paired solve-rate margin in percent")
    ap.add_argument("--paired-context", required=True)
    ap.add_argument("--paired-warmups", type=int, default=2)
    ap.add_argument("--paired-inner-reps", type=int, default=1)
    ap.add_argument("--max-overhead", type=int, default=512)
    ap.add_argument("--cache-state", choices=["warm", "cold"], default="warm")
    ap.add_argument("--evict-bytes", type=int, default=64 * 1024 * 1024)
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
    a = ap.parse_args(argv)
    if (a.paired_replicates < 4 or a.paired_replicates % 2 != 0 or
            a.paired_warmups < 0 or a.paired_warmups % 2 != 0 or
            not 1 <= a.paired_inner_reps <= 1024 or
            a.max_overhead < 0 or
            not 4096 <= a.evict_bytes <= PEELTIMING_EVICT_BYTES_MAX or
            (a.cache_state == "cold" and a.paired_inner_reps != 1) or
            a.bb < 2 or a.bb % 2 != 0 or a.kmax < 2 or
            not math.isfinite(a.margin) or not 0.0 <= a.margin <= 100.0 or
            not valid_loss_rate(a.loss) or
            not 0 <= a.construction_seed <= 0xffffffffffffffff or
            not 0 <= a.loss_seed <= 0xffffffffffffffff):
        ap.error(
            "invalid paired settings, payload, K, margin, loss, or uint64 seed")
    required_margin = a.margin / 100.0
    if a.out is not None:
        try:
            thermal_source = require_safe_output(
                a.table, a.out, a.bench, a.paired_context)
        except MeasurementError as error:
            print(f"  REFUSED output: {error}", file=sys.stderr)
            return 2

    try:
        identity = capture_artifact_identity(
            a.bench, "tools/peel_validate.py")
        source_document, table, source_table_sha256 = (
            read_table_document_snapshot(a.table))
        verify_benchmark_identity(source_document, a.bench)
        if source_document["provenance"]["generator"] == "tools/peel_validate.py":
            raise MeasurementError("input table is already a validation result")
    except (MeasurementError, OSError, ValueError) as error:
        print(f"  REFUSED input table: {error}", file=sys.stderr)
        return 2
    ks = [k for k in sorted(table) if k <= a.kmax]
    if not ks:
        print("  REFUSED input table: selected K set is empty", file=sys.stderr)
        return 2
    try:
        for k in ks:
            validate_peeltiming_dimensions(
                block_count=k,
                block_bytes=a.bb,
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
                required_margin=required_margin,
            )
    except ValueError as error:
        ap.error(str(error))

    print(
        f"  {len(ks)} K values, {a.paired_replicates} paired replicates, "
        f"bb={a.bb}\n")
    print(
        f"  {'K':>6} {'candidate':>9} {'shipped':>9} {'solve':>8} "
        f"{'ca OH':>7} {'sh OH':>7} {'verdict':>10}")
    fixed, wins, losses, failed = {}, 0, 0, []
    for k in ks:
        source = table[k]
        try:
            native = stock_profile(
                a.bench, k, target_profile=a.target_profile)
            if native.as_dict() != source["native"]:
                raise MeasurementError(
                    "native peel profile does not match the search receipt")
            source_mode = source["search_receipt"]["mode"]
            source_is_shipped = source_mode == "shipped"
            scale = source["search_receipt"]["coordinates"]["scale"]
            construction_seed = derive_seed(
                a.construction_seed, "validation", k,
                a.paired_replicates, a.bb, "construction")
            loss_seed = derive_seed(
                a.loss_seed, "validation", k,
                a.paired_replicates, a.bb, "loss")
            measurement = paired_probe(
                a.bench, k, a.bb, source["search_receipt"]["peel_pmf"],
                degree_scale=None if scale == -1.0 else scale,
                native_profile=native,
                target_profile=a.target_profile,
                seed_policy=a.seed_policy,
                loss=a.loss,
                schedule=a.schedule,
                construction_seed=construction_seed,
                loss_seed=loss_seed,
                warmup_replicates=a.paired_warmups,
                replicates=a.paired_replicates,
                inner_reps=a.paired_inner_reps,
                max_overhead=a.max_overhead,
                cache_state=a.cache_state,
                evict_bytes=a.evict_bytes,
                context=a.paired_context,
                required_margin=required_margin,
            )
            if source_is_shipped:
                require_paired_stock_control(
                    measurement, "paired validation shipped-source control")
            keep = (
                not source_is_shipped and
                measurement.valid_for_promotion
            )
            selected = (
                measurement.candidate if keep else measurement.identity)
        except (KeyError, TypeError, ValueError, OSError,
                MeasurementError) as error:
            failed.append(k)
            print(f"  {k:>6}   measurement failed: {error}", flush=True)
            continue

        candidate = measurement.candidate
        shipped = measurement.identity
        candidate_goodput = candidate.goodput(k)
        shipped_goodput = shipped.goodput(k)
        solve_gain = (
            100.0 * (candidate.solve_mbps - shipped.solve_mbps) /
            shipped.solve_mbps if shipped.solve_mbps > 0.0 else None)
        validation_receipt = {
            "protocol": NATIVE_PAIRED_PROTOCOL,
            "verdict": "keep" if keep else "control",
            "selected_arm": "candidate" if keep else "identity",
            "source_mode": source_mode,
            "trials": a.paired_replicates,
            "block_bytes": a.bb,
            "scale": scale,
            "trained_pmf_sha256":
                source["search_receipt"]["peel_pmf_sha256"],
            "trained_goodput": candidate_goodput,
            "shipped_goodput": shipped_goodput,
            "paired_measurement": measurement.as_dict(),
        }
        shared = {
            **selected.as_dict(),
            "goodput": selected.goodput(k),
            "verified_solve_mbps": selected.solve_mbps,
            "verified_oh": selected.oh_mean,
            "shipped_solve_mbps": shipped.solve_mbps,
            "solve_gain_pct": (
                canonical_percent(solve_gain)
                if keep and solve_gain is not None
                else 0.0),
            "validation_receipt": validation_receipt,
        }
        if keep:
            wins += 1
            fixed[k] = dict(source, **shared)
        else:
            losses += 1
            fixed[k] = {
                "K": k,
                "scale": -1.0,
                "p1": 100,
                "tilt": 0,
                "dmax": 64,
                "absorb": 100,
                "reverted_to_shipped": True,
                **shared,
                "search_would_have_solve_delta_pct": (
                    canonical_percent(solve_gain)
                    if solve_gain is not None else None),
                "native": source["native"],
                "peel_pmf": list(native.pmf),
                "search_receipt": source["search_receipt"],
            }
        delta_text = (
            f"{solve_gain:+7.1f}%" if solve_gain is not None else "     n/a")
        print(
            f"  {k:>6} {candidate.solve_mbps:>9.1f} "
            f"{shipped.solve_mbps:>9.1f} {delta_text:>8} "
            f"{candidate.oh_mean:>7.4f} {shipped.oh_mean:>7.4f} "
            f"{'keep' if keep else 'REVERT':>10}"
            + ("" if candidate.fail == 0 and shipped.fail == 0 else
               f"  fail ca={candidate.fail} sh={shipped.fail}"),
            flush=True)
    if failed:
        print(
            f"\n  REFUSED publication: {len(failed)} of {len(ks)} "
            f"measurements failed: {failed}", file=sys.stderr)
        return 1

    gains = [
        value["solve_gain_pct"] for value in fixed.values()
        if not value.get("reverted_to_shipped")
    ]
    print(f"\n  {wins} kept, {losses} selected true shipped")
    if gains:
        gains.sort()
        print(
            f"  solve gains where kept: median "
            f"{gains[len(gains) // 2]:+.1f}%  min {gains[0]:+.1f}%  "
            f"max {gains[-1]:+.1f}%")
    if a.out:
        try:
            _, _, current_source_sha256 = read_table_document_snapshot(a.table)
            if current_source_sha256 != source_table_sha256:
                raise MeasurementError(
                    "source table changed during validation")
            require_safe_output(
                a.table, a.out, a.bench, a.paired_context,
                expected_thermal_source=thermal_source)
        except MeasurementError as error:
            print(f"  REFUSED publication: {error}", file=sys.stderr)
            return 1
        try:
            output_document = make_table_document(
                fixed,
                generator="tools/peel_validate.py",
                bench=a.bench,
                construction_seed_base=a.construction_seed,
                loss_seed_base=a.loss_seed,
                target_profile=a.target_profile,
                seed_policy=a.seed_policy,
                loss=a.loss,
                schedule=a.schedule,
                settings={
                    "source_table": os.path.realpath(a.table),
                    "source_table_sha256": source_table_sha256,
                    "source_entry_count": len(table),
                    "selected_entry_count": len(ks),
                    "selected_K": list(ks),
                    "paired_replicates": a.paired_replicates,
                    "block_bytes": a.bb,
                    "kmax": a.kmax,
                    "margin_percent": a.margin,
                    "paired_context": os.path.realpath(a.paired_context),
                    "paired_warmups": a.paired_warmups,
                    "paired_inner_reps": a.paired_inner_reps,
                    "max_overhead": a.max_overhead,
                    "cache_state": a.cache_state,
                    "evict_bytes": a.evict_bytes,
                    "target_profile": a.target_profile,
                    "seed_policy": a.seed_policy,
                    "loss": a.loss,
                    "schedule": a.schedule,
                },
                source_provenance={
                    "schema": source_document["schema"],
                    "document_sha256": source_table_sha256,
                    "provenance": source_document["provenance"],
                    "entry_count": len(table),
                    "selected_entry_count": len(ks),
                    "selected_K": list(ks),
                },
                artifact_identity=identity,
            )
            require_distinct_from_source_provenance(
                a.out, "--out", "tools/peel_validate.py")
            write_json_atomic(a.out, output_document)
        except (MeasurementError, OSError, ValueError) as error:
            print(f"  REFUSED publication: {error}", file=sys.stderr)
            return 1
        print(f"  corrected table -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
