#!/usr/bin/env python3
"""Validate a versioned peel table against the true shipped codec.

Candidate and true shipped arms use explicit paired construction/loss seeds.
The candidate replays its measured PMF and scale, while the shipped arm uses
neither hook. Every recovery metric emitted by the native benchmark is
retained. Legacy unversioned tables are refused: their PMF source and
measurement conditions cannot be replayed. The output is replaced atomically
only after all requested K values succeed.
"""
import argparse
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from peel_codec import (                                 # noqa: E402
    MeasurementError,
    capture_artifact_identity,
    compare_probe,
    derive_seed,
    make_table_document,
    read_table_document_snapshot,
    stock_profile,
    valid_loss_rate,
    verify_benchmark_identity,
    write_json_atomic,
)


def require_distinct_output(source_path, output_path):
    """Fail closed unless an output path is demonstrably distinct."""
    try:
        if os.path.realpath(source_path) == os.path.realpath(output_path):
            raise MeasurementError(
                "--out must not replace the source table")
        if os.path.lexists(output_path) and os.path.samefile(
                source_path, output_path):
            raise MeasurementError(
                "--out must not replace the source table")
    except OSError as error:
        raise MeasurementError(
            f"could not establish that --out is distinct from the "
            f"source table: {error}")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", default="build-fast/codec/wirehair_v2_bench")
    ap.add_argument("--table", required=True)
    ap.add_argument("--out", default=None, help="write the corrected table here")
    ap.add_argument("--trials", type=int, default=3000)
    ap.add_argument("--bb", type=int, required=True)
    ap.add_argument("--kmax", type=int, default=10 ** 9)
    ap.add_argument("--margin", type=float, default=2.0,
                    help="percent goodput a trained point must beat shipped by")
    ap.add_argument("--target-profile", required=True, choices=["dispatch-v1"])
    ap.add_argument("--seed-policy", required=True, choices=["raw"])
    ap.add_argument("--schedule", required=True, choices=["iid"])
    ap.add_argument("--loss", type=float, required=True)
    ap.add_argument("--construction-seed", type=int, required=True)
    ap.add_argument("--loss-seed", type=int, required=True)
    a = ap.parse_args(argv)
    if (a.trials < 1 or a.bb < 1 or a.kmax < 2 or
            not math.isfinite(a.margin) or a.margin < 0.0 or
            not valid_loss_rate(a.loss) or
            not 0 <= a.construction_seed <= 0xffffffffffffffff or
            not 0 <= a.loss_seed <= 0xffffffffffffffff):
        ap.error(
            "invalid trials, payload, K, margin, loss, or uint64 seed")
    if a.out is not None:
        try:
            require_distinct_output(a.table, a.out)
        except MeasurementError as error:
            print(
                f"  REFUSED output: {error}",
                file=sys.stderr)
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

    print(f"  {len(ks)} K values, {a.trials} trials each arm, bb={a.bb}\n")
    print(f"  {'K':>6} {'trained':>9} {'shipped':>9} {'delta':>8} "
          f"{'tr OH':>7} {'sh OH':>7} {'verdict':>10}")
    fixed, wins, losses, failed = {}, 0, 0, []
    for k in ks:
        v = table[k]
        try:
            scale = float(v["scale"])
            native = stock_profile(
                a.bench, k, target_profile=a.target_profile)
            if native.as_dict() != v["native"]:
                raise MeasurementError(
                    "native peel profile does not match the search receipt")
            construction_seed = derive_seed(
                a.construction_seed, "validation", k, a.trials, a.bb,
                "construction")
            loss_seed = derive_seed(
                a.loss_seed, "validation", k, a.trials, a.bb, "loss")
            exact = {
                "native_profile": native,
                "target_profile": a.target_profile,
                "seed_policy": a.seed_policy,
                "loss": a.loss,
                "schedule": a.schedule,
                "construction_seed": construction_seed,
                "loss_seed": loss_seed,
            }
            source_mode = v["search_receipt"]["mode"]
            source_is_shipped = source_mode == "shipped"
            source_is_scale_only = source_mode == "scale-only"
            if source_is_shipped:
                tr = compare_probe(
                    a.bench, k, a.trials, a.bb, **exact)
            elif source_is_scale_only:
                # A stock PMF with an active scale is a real candidate, not
                # the shipped codec. Leave the peel hook unset and replay only
                # its staircase scale.
                tr = compare_probe(
                    a.bench, k, a.trials, a.bb,
                    degree_scale=scale, **exact)
            else:
                tr = compare_probe(
                    a.bench, k, a.trials, a.bb,
                    peel_weights=v["peel_pmf"],
                    degree_scale=scale if scale >= 0.0 else None,
                    **exact)
            sh = compare_probe(
                a.bench, k, a.trials, a.bb, **exact)
            if tr.fail > 0 and sh.fail > 0:
                raise MeasurementError(
                    "both trained and shipped validation arms failed recovery")
        except (KeyError, TypeError, ValueError, MeasurementError) as error:
            failed.append(k)
            print(f"  {k:>6}   measurement failed: {error}", flush=True)
            continue
        tg = tr.goodput(k)
        sg = sh.goodput(k)
        # Shipped is the safer default, so training must clear an explicit
        # margin rather than win a single noisy timing comparison.
        keep = (
            not source_is_shipped and
            tg > sg * (1.0 + a.margin / 100.0)
        )
        selected = tr if keep else sh
        if selected.fail != 0:
            failed.append(k)
            print(
                f"  {k:>6}   measurement failed: selected "
                f"{'trained' if keep else 'shipped'} arm did not recover",
                flush=True)
            continue
        d = 100.0 * (tg - sg) / sg if sg > 0.0 else None
        validation_receipt = {
            "verdict": "keep" if keep else "control",
            "margin_percent": a.margin,
            "trials": a.trials,
            "block_bytes": a.bb,
            "scale": scale,
            "trained_pmf_sha256": tr.target_receipt["pmf_sha256"],
            "trained_goodput": tg,
            "shipped_goodput": sg,
            "trained": tr.as_dict(),
            "shipped": sh.as_dict(),
        }
        if keep:
            wins += 1
            fixed[k] = dict(
                v, **tr.as_dict(), goodput=tg, scale=scale,
                peel_pmf=list(v["peel_pmf"]),
                verified_mbps=tr.decode_mbps,
                verified_oh=tr.oh_mean, shipped_mbps=sh.decode_mbps,
                gain_pct=round(d, 2) if d is not None else None,
                validation_receipt=validation_receipt)
        else:
            losses += 1
            fixed[k] = {
                "K": k, "scale": -1.0, "p1": 100, "tilt": 0,
                "dmax": 64, "absorb": 100,
                "reverted_to_shipped": True,
                **sh.as_dict(),
                "goodput": sg,
                "verified_mbps": sh.decode_mbps,
                "verified_oh": sh.oh_mean,
                "shipped_mbps": sh.decode_mbps,
                "gain_pct": 0.0,
                "search_would_have_lost_pct": (
                    round(d, 2) if d is not None else None),
                "native": v["native"],
                "peel_pmf": list(native.pmf),
                "search_receipt": v["search_receipt"],
                "validation_receipt": validation_receipt,
            }
        delta_text = f"{d:+7.1f}%" if d is not None else "     n/a"
        print(f"  {k:>6} {tr.decode_mbps:>9.1f} "
              f"{sh.decode_mbps:>9.1f} {delta_text:>8} "
              f"{tr.oh_mean:>7.4f} {sh.oh_mean:>7.4f} "
              f"{'keep' if keep else 'REVERT':>10}"
              + ("" if tr.fail == 0 and sh.fail == 0 else
                 f"  fail tr={tr.fail} sh={sh.fail}"), flush=True)
    if failed:
        print(
            f"\n  REFUSED publication: {len(failed)} of {len(ks)} "
            f"measurements failed: {failed}", file=sys.stderr)
        return 1

    g = [
        v["gain_pct"] for v in fixed.values()
        if not v.get("reverted_to_shipped") and v["gain_pct"] is not None
    ]
    print(f"\n  {wins} kept, {losses} selected true shipped")
    if g:
        g.sort()
        print(f"  gains where kept: median {g[len(g)//2]:+.1f}%  "
              f"min {g[0]:+.1f}%  max {g[-1]:+.1f}%")
    if a.out:
        try:
            _, _, current_source_sha256 = read_table_document_snapshot(a.table)
            if current_source_sha256 != source_table_sha256:
                raise MeasurementError(
                    "source table changed during validation")
            require_distinct_output(a.table, a.out)
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
                    "trials": a.trials,
                    "block_bytes": a.bb,
                    "kmax": a.kmax,
                    "margin_percent": a.margin,
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
            write_json_atomic(a.out, output_document)
        except (MeasurementError, OSError, ValueError) as error:
            print(f"  REFUSED publication: {error}", file=sys.stderr)
            return 1
        print(f"  corrected table -> {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
