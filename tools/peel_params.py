#!/usr/bin/env python3
"""Inspect an independently validated, provenance-bearing peel table.

Legacy flat JSON and raw search tables are refused by default. Exact validated
anchors are replayable; every unmeasured K resolves to shipped by default.
Raw-search inspection, interpolation, and clamping each require an explicit
experimental opt-in because search evidence is not independent validation and
an arm measured at one K has no recovery guarantee at another.
"""
import argparse
import bisect
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_TABLE = os.path.join(HERE, "peel_table.json")
sys.path.insert(0, HERE)
from peel_codec import (                                  # noqa: E402
    MeasurementError,
    family,
    read_table_document,
    stock_pmf,
    verify_benchmark_identity,
)

FIELDS = ("scale", "p1", "tilt", "dmax", "absorb")


def load(
        path=DEFAULT_TABLE, allow_unverified_cost_model=False, bench=None,
        allow_unvalidated_search=False):
    document, entries = read_table_document(path)
    if bench is not None:
        verify_benchmark_identity(document, bench)
    generator = document["provenance"]["generator"]
    if (generator != "tools/peel_validate.py" and
            not allow_unvalidated_search):
        raise MeasurementError(
            "peel table is an unvalidated search result; pass "
            "--allow-unvalidated-search only for an explicit experiment, "
            "or independently validate it into a new table first")
    settings = document["provenance"]["settings"]
    if (settings.get("allow_unverified_cost_model") and
            not allow_unverified_cost_model):
        raise MeasurementError(
            "peel table was selected with an unverified proxy cost model; "
            "pass --allow-unverified-cost-model only for an explicit "
            "experiment, or validate it into a new table first")
    return entries


def shipped_params():
    return {
        "scale": -1.0,
        "p1": 100,
        "tilt": 0,
        "dmax": 64,
        "absorb": 100,
    }


def params_for_k(k, table, allow_unverified_interpolation=False):
    """Resolve one K, using shipped unless unverified synthesis is opted in.

    With the explicit experimental opt-in, interior trained anchors are
    interpolated and values outside the measured range are clamped to the
    nearest anchor. A shipped-control boundary still propagates shipped rather
    than creating a meaningless hybrid with the ``scale=-1`` sentinel.
    """
    if (isinstance(k, bool) or not isinstance(k, int) or
            not 2 <= k <= 64000):
        raise MeasurementError("K must be an integer in [2,64000]")
    ks = sorted(table)
    if not ks:
        raise MeasurementError("empty anchor table")
    if k in table:
        return {f: table[k][f] for f in FIELDS}, "anchor"
    if not allow_unverified_interpolation:
        return shipped_params(), "shipped-unmeasured"
    if k <= ks[0]:
        if table[ks[0]].get("reverted_to_shipped"):
            return shipped_params(), "shipped-control"
        return (
            {f: table[ks[0]][f] for f in FIELDS},
            f"EXPERIMENTAL clamped-low from {ks[0]}",
        )
    if k >= ks[-1]:
        if table[ks[-1]].get("reverted_to_shipped"):
            return shipped_params(), "shipped-control"
        return (
            {f: table[ks[-1]][f] for f in FIELDS},
            f"EXPERIMENTAL clamped-high from {ks[-1]}",
        )
    i = bisect.bisect_left(ks, k)
    lo, hi = ks[i - 1], ks[i]
    if (table[lo].get("reverted_to_shipped") or
            table[hi].get("reverted_to_shipped") or
            ((table[lo]["scale"] < 0.0) != (table[hi]["scale"] < 0.0))):
        return shipped_params(), f"shipped-control {lo}..{hi}"
    t = (k - lo) / (hi - lo)
    out = {}
    for f in FIELDS:
        a, b = table[lo][f], table[hi][f]
        v = a + t * (b - a)
        out[f] = v if f == "scale" else int(round(v))
    return out, f"EXPERIMENTAL interp {lo}..{hi}"


def uses_shipped_arm(k, table, allow_unverified_interpolation=False):
    """Return whether K resolves to the no-override shipped control."""
    if (isinstance(k, bool) or not isinstance(k, int) or
            not 2 <= k <= 64000):
        raise MeasurementError("K must be an integer in [2,64000]")
    ks = sorted(table)
    if not ks:
        raise MeasurementError("empty anchor table")
    if k in table:
        return bool(table[k].get("reverted_to_shipped"))
    if not allow_unverified_interpolation:
        return True
    if k <= ks[0]:
        return bool(table[ks[0]].get("reverted_to_shipped"))
    if k >= ks[-1]:
        return bool(table[ks[-1]].get("reverted_to_shipped"))
    i = bisect.bisect_left(ks, k)
    return bool(
        table[ks[i - 1]].get("reverted_to_shipped") or
        table[ks[i]].get("reverted_to_shipped") or
        ((table[ks[i - 1]]["scale"] < 0.0) !=
         (table[ks[i]]["scale"] < 0.0)))


def pmf_for_k(
        k, table, bench, allow_unverified_interpolation=False, *,
        target_profile):
    """The peel PMF this K should use, based on the exact native PMF."""
    p, _ = params_for_k(k, table, allow_unverified_interpolation)
    stock = stock_pmf(bench, k, target_profile=target_profile)
    if uses_shipped_arm(k, table, allow_unverified_interpolation):
        return stock
    if k in table:
        # Search receipts preserve the exact binary64 law measured by compare;
        # do not rebuild it from shorthand coordinates on an anchor replay.
        return list(table[k]["peel_pmf"])
    return family(stock, p["p1"], p["tilt"], p["dmax"], p["absorb"])


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", default=DEFAULT_TABLE)
    ap.add_argument("--bench", default="build-fast/codec/wirehair_v2_bench")
    ap.add_argument("--target-profile", required=True, choices=["dispatch-v1"])
    ap.add_argument("--K", type=int, help="print parameters for one K")
    ap.add_argument("--pmf", action="store_true", help="also print the PMF")
    ap.add_argument("--check", action="store_true",
                    help="hold out each interior anchor and interpolate it")
    ap.add_argument(
        "--allow-unverified-cost-model", action="store_true",
        help="accept a sweep table that records the unverified proxy opt-in")
    ap.add_argument(
        "--allow-unvalidated-search", action="store_true",
        help="EXPERIMENTAL: inspect a raw search result without independent "
             "validation")
    ap.add_argument(
        "--allow-unverified-interpolation", action="store_true",
        help="EXPERIMENTAL: interpolate or clamp trained coordinates for an "
             "unmeasured K instead of using shipped")
    a = ap.parse_args(argv)
    if a.K is not None and not 2 <= a.K <= 64000:
        ap.error("--K must be in [2,64000]")
    if a.pmf and a.K is None:
        ap.error("--pmf requires --K")
    try:
        table = load(
            a.table, a.allow_unverified_cost_model, a.bench,
            a.allow_unvalidated_search)
    except (MeasurementError, OSError, ValueError) as error:
        print(f"REFUSED input table: {error}", file=sys.stderr)
        return 2

    if a.check:
        ks = sorted(table)
        print(f"  EXPERIMENTAL HOLD-OUT CHECK: drop each interior anchor, "
              f"interpolate from "
              f"its neighbours\n")
        print(f"  {'K':>7} {'tilt true':>10} {'tilt interp':>12} {'err':>8} "
              f"{'scale true':>11} {'scale interp':>13}")
        for k in ks[1:-1]:
            held = {kk: v for kk, v in table.items() if kk != k}
            got, _ = params_for_k(
                k, held, allow_unverified_interpolation=True)
            tt, gt = table[k]["tilt"], got["tilt"]
            err = abs(gt - tt)
            print(f"  {k:>7} {tt:>+10} {gt:>+12} {err:>8} "
                  f"{table[k]['scale']:>11.2f} {got['scale']:>13.2f}")
        return 0

    if a.K is not None:
        p, how = params_for_k(
            a.K, table, a.allow_unverified_interpolation)
        if a.allow_unverified_interpolation and a.K not in table:
            print(
                "WARNING: using unverified trained coordinates at an "
                "unmeasured K",
                file=sys.stderr)
        print(f"  K={a.K} ({how}): scale={p['scale']:.2f} p1={p['p1']} "
              f"tilt={p['tilt']} dmax={p['dmax']} absorb={p['absorb']}")
        if a.pmf:
            w = pmf_for_k(
                a.K, table, a.bench, a.allow_unverified_interpolation,
                target_profile=a.target_profile)
            print("  " + ",".join(f"{x:.17g}" for x in w))
        return 0

    print(f"  {len(table)} anchors from {a.table}\n")
    print(f"  {'K':>7} {'scale':>8} {'p1':>5} {'tilt':>6} {'dmax':>5} "
          f"{'absorb':>7} {'MB/s':>9} {'OH':>8}")
    for k in sorted(table):
        v = table[k]
        print(f"  {k:>7} {v['scale']:>8.2f} {v['p1']:>5} {v['tilt']:>+6} "
                  f"{v['dmax']:>5} {v['absorb']:>7} "
                  f"{v.get('decode_mbps', v.get('verified_mbps', 0)):>9.1f} "
                  f"{v.get('oh_mean', 0):>8.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
