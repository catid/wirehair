#!/usr/bin/env python3
"""Render the WH2 untuned-Pareto screen to a self-contained HTML page.

The screen plots expected decoding overhead against encoder throughput for every
precode tuning point in the sealed campaign, faceted across the six K regions.

Two inputs, kept deliberately separate:

  * the sealed per-method result JSONs, which carry the overhead statistics, and
  * one timing JSON, which carries the throughput used for the x-axis.

The split matters. A result JSON also contains the throughput measured during the
run that produced it, and that number is *not* what belongs on the axis: round-1
JSONs carry contended 60-worker timing and round-2 JSONs carry pre-fast-path loki
timing. Only the quiet single-host re-time is comparable across the whole grid, so
the timing file is always supplied explicitly and always wins.

Usage
-----
  # extract from the sealed campaign, cache the distilled inputs, and render
  render_untuned_pareto_html.py \
      --v1 <sealed_v1_dir> --v2 <sealed_v2_dir> --timing <quiet_timing.json> \
      --desc bench/results/pareto/descriptions.json \
      --prose bench/results/pareto/prose.json \
      --rows-out bench/results/pareto/rows.json \
      --alt-out bench/results/pareto/r1_common.json \
      --out results/untuned_pareto_screen.html

  # re-render from the cached inputs -- no access to the sealed campaign needed
  render_untuned_pareto_html.py \
      --rows bench/results/pareto/rows.json --alt bench/results/pareto/r1_common.json \
      --desc bench/results/pareto/descriptions.json \
      --prose bench/results/pareto/prose.json \
      --out results/untuned_pareto_screen.html

Add --fragment to emit body-only HTML for a host that supplies its own <head>.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
import pathlib
import sys

REGIONS = [
    ("r1_2_100", "K 2 – 100"),
    ("r2_100_1000", "K 100 – 1,000"),
    ("r3_1000_5000", "K 1,000 – 5,000"),
    ("r4_5000_10000", "K 5,000 – 10,000"),
    ("r5_10000_20000", "K 10,000 – 20,000"),
    ("r6_20000_64000", "K 20,000 – 64,000"),
]

WH1 = "wh1-legacy"
SELECTED = "wh2-dispatched-v1"

# This alternate r1 view predates the paired-cell ledgers.  The sealed
# cells_detail arrays retain only [K, outcome], not the construction/loss seed
# identities needed to prove that two arms saw the same cell.  Keep the cache
# self-describing and fail closed: this is useful for diagnosing the historical
# plot, but it must never be promoted into architecture-selection evidence.
HISTORICAL_DIAGNOSTIC_SCHEMA = "wirehair.untuned-pareto.r1-historical-diagnostic.v2"
ROWS_CACHE_SCHEMA = "wirehair.untuned-pareto.rows.v1"
RESULT_SOURCE_SCHEMA = "wirehair.untuned-pareto.result-sources.v1"
HISTORICAL_DIAGNOSTIC_REGION = "r1_2_100"
HISTORICAL_DIAGNOSTIC_K_MIN = 6
HISTORICAL_DIAGNOSTIC_K_MAX = 100
HISTORICAL_DIAGNOSTIC_CELLS_PER_K = 10
HISTORICAL_DIAGNOSTIC_CELLS_PER_ARM = 950
HISTORICAL_DIAGNOSTIC_CLASSIFICATION = {
    "historical": True,
    "diagnostic_only": True,
    "selection_eligible": False,
    "architecture_evidence_eligible": False,
    "seed_ids_present": False,
}
HISTORICAL_DIAGNOSTIC_GRID = {
    "region": HISTORICAL_DIAGNOSTIC_REGION,
    "k_min": HISTORICAL_DIAGNOSTIC_K_MIN,
    "k_max": HISTORICAL_DIAGNOSTIC_K_MAX,
    "k_count": HISTORICAL_DIAGNOSTIC_K_MAX - HISTORICAL_DIAGNOSTIC_K_MIN + 1,
    "cells_per_k": HISTORICAL_DIAGNOSTIC_CELLS_PER_K,
    "cells_per_arm": HISTORICAL_DIAGNOSTIC_CELLS_PER_ARM,
}
HISTORICAL_DIAGNOSTIC_NOTICE = (
    "Historical diagnostic only; non-selectable and never architecture evidence. "
    "cells_detail has no seed IDs, so cross-arm cell identity cannot be verified."
)

# The historical inputs are composite: round two's dispatch arm was appended
# after its original index/provenance were sealed.  We therefore authenticate
# the exact result-file sets directly instead of trusting the stale index, and
# preserve their Git tree IDs as provenance.  Any data refresh requires an
# explicit schema/identity update rather than silently accepting a same-shaped
# campaign.
HISTORICAL_RESULT_SOURCES = {
    "schema": RESULT_SOURCE_SCHEMA,
    "metadata_status": (
        "exact-result-files-authenticated-directly; round2 index/provenance "
        "predate the dispatch append and are not authority"
    ),
    "round1": {
        "git_tree_oid": "d5dd90d18cb1d442414997f690e3ebec6a74ca52",
        "result_file_count": 63,
        "result_file_manifest_sha256":
            "8e12eae35125c0dfaf55e1c517f9ecb5a8987b4303d0da75760a95558897df34",
    },
    "round2": {
        "git_tree_oid": "cfd3d0b86ece0d3ffc005a4f4c4363999f7f3625",
        "result_file_count": 21,
        "result_file_manifest_sha256":
            "becc3ec542ec9645e3d22edec89d84008cffade86fdda1fbc319ed00c46eb404",
    },
}
CANONICAL_TIMING_SOURCE = {
    "git_blob_oid": "713c9374a58000144241265218d7028d57cbd249",
    "sha256": "10ee9918b4a732bb14abcbba828c40ae44197508a665d0b946bc02bde0878068",
}
ROWS_SOURCE_IDENTITIES = {
    "result_sources": HISTORICAL_RESULT_SOURCES,
    "timing": CANONICAL_TIMING_SOURCE,
}
HISTORICAL_DIAGNOSTIC_ARMS_SHA256 = (
    "6c9666703e526b9006a68403ac5395eefdd9ad47c019875f67bbdc9f322d2eb2"
)
ROWS_PAYLOAD_SHA256 = (
    "ed1c2a1779eb87d3a63338567ff1d9c5e41c019f1a9ddf8cc61727465c82643f"
)

# Hue encodes the precode architecture, shape the sub-family within it.
#
# Only four categorical hues are used. That set -- blue / green / magenta / yellow --
# is one of exactly two 4-subsets of the documented palette that clear the all-pairs
# CVD and normal-vision floors in BOTH light and dark mode, which is the gate a
# scatter has to pass because any two marks can end up adjacent. Two validator
# warnings survive, each discharged by secondary encoding:
#
#   light, contrast < 3:1 for magenta (2.62) and yellow (2.11)
#       -> relief is the table view plus direct labels, both shipped.
#   dark, protanopia dE 6.9 for yellow vs green
#       -> yellow is Wirehair 1 alone, drawn as a square and direct-labelled in
#          every facet, so hue never carries that distinction by itself.
GROUPS = {
    "g1": ("Mixed — period & dense-row count", "--s1", "c"),
    "g2": ("Mixed — row shape & H15", "--s1", "d"),
    "g3": ("Certified (all-GF(256))", "--s2", "c"),
    "g4": ("Certified + identity corner", "--s2", "d"),
    "g5": ("Identity corner (mixed)", "--s3", "c"),
    "g6": ("Anchors, staircase & masks", "--s3", "d"),
    "g7": ("Wirehair 1 (legacy reference)", "--s4", "s"),
    "g8": ("WH2 dispatch profile v1 (selected)", "--s1", "r"),
}

# Sealed-data family -> display group. The family string comes from the result JSON,
# so the grouping follows what was actually run rather than the method's name.
FAMILY_GROUP = {
    "wh2-mixed-frozen": "g1",
    "wh2-mixed-sharedx": "g1",
    "wh2-r2-denseknee": "g1",
    "wh2-mixed-h15": "g2",
    "wh2-mixed-11rows": "g2",
    "wh2-mixed-11x4": "g2",
    "wh2-mixed-12rows": "g2",
    "wh2-mixed-ramp": "g2",
    "wh2-mixed-skew": "g2",
    "wh2-r2-rowtrim": "g2",
    "wh2-certified": "g3",
    "wh2-r2-cert-ic": "g4",
    "wh2-r2-identcorner": "g5",
    "wh2-identity-corner": "g5",
    "wh2-r2-kitchensink": "g5",
    "wh2-two-anchor": "g6",
    "wh2-r2-twoanchor": "g6",
    "wh2-segmented-anchors": "g6",
    "wh2-degree-balanced": "g6",
    "wh2-grouped-masks": "g6",
    "wh1": "g7",
    "wh2-dispatched": "g8",
}

ROW_FIELDS = ("m", "fam", "reg", "eoh", "sps", "f0", "f2", "dead", "cf", "hk", "n", "exc", "dl", "sk")


def extract_rows(
        v1: pathlib.Path,
        v2: pathlib.Path,
        timing_path: pathlib.Path,
        *,
        expected_result_sources: dict | None = None,
        expected_timing_source: dict | None = None) -> list[dict]:
    """Join the sealed per-method results to the canonical timing baseline."""
    if expected_result_sources is None:
        expected_result_sources = HISTORICAL_RESULT_SOURCES
    if expected_timing_source is None:
        expected_timing_source = CANONICAL_TIMING_SOURCE
    _authenticated_sources, result_documents = _authenticated_result_documents(
        v1, v2, expected_result_sources,
    )
    _authenticated_timing, timing = _authenticated_timing_document(
        timing_path, expected_timing_source,
    )
    if type(timing) is not dict:
        raise SystemExit("canonical timing input must be a JSON object")
    results: dict[str, dict] = {}
    for path, doc in result_documents:
        if type(doc) is not dict:
            raise SystemExit(f"{path}: historical result must be a JSON object")
        name = doc["method"]
        if name in results:
            raise SystemExit(f"method {name} appears in both result sets; refusing to guess")
        results[name] = doc

    missing_timing = sorted(set(results) - set(timing))
    missing_results = sorted(set(timing) - set(results))
    if missing_timing or missing_results:
        raise SystemExit(
            f"result/timing method sets disagree; only in results: {missing_timing}; "
            f"only in timing: {missing_results}"
        )

    rows = []
    for name in sorted(results):
        doc = results[name]
        for reg, _ in REGIONS:
            r = doc["regions"][reg]
            work = r.get("work") or {}
            rows.append(
                {
                    "m": name,
                    "fam": doc["family"],
                    "reg": reg,
                    "eoh": r["expected_overhead"],
                    "sps": timing[name][reg]["median_symbols_per_sec"],
                    "f0": r["fail_at_0_fraction"],
                    "f2": r["fail_at_2_fraction"],
                    "dead": r["dead_fraction"],
                    "cf": r["construct_failure_fraction"],
                    "hk": r["hard_k_fraction"],
                    "n": r["cells"],
                    "exc": r["excluded_cells"],
                    "dl": r["decoder_limit_fraction"],
                    "sk": r["sampled_k"],
                    "xors": work.get("block_xors_mu"),
                    "mula": work.get("block_muladds_mu"),
                }
            )
    return rows


def _require_exact_keys(value: object, expected: set[str], label: str) -> dict:
    if type(value) is not dict:
        raise SystemExit(f"{label} must be a JSON object")
    actual = set(value)
    if actual != expected:
        raise SystemExit(
            f"{label} fields disagree; missing: {sorted(expected - actual)}; "
            f"unexpected: {sorted(actual - expected)}"
        )
    return value


def _require_exact_value(value: object, expected: object, label: str) -> None:
    """Compare JSON values without Python's bool/int or int/float coercion."""
    if type(value) is not type(expected):
        raise SystemExit(
            f"{label} type mismatch: expected {type(expected).__name__}, "
            f"got {type(value).__name__}"
        )
    if type(expected) is dict:
        value_dict = _require_exact_keys(value, set(expected), label)
        for key in expected:
            _require_exact_value(value_dict[key], expected[key], f"{label}.{key}")
        return
    if type(expected) is list:
        if len(value) != len(expected):
            raise SystemExit(
                f"{label} length mismatch: expected {len(expected)}, got {len(value)}"
            )
        for index, (item, expected_item) in enumerate(zip(value, expected)):
            _require_exact_value(item, expected_item, f"{label}[{index}]")
        return
    if value != expected:
        raise SystemExit(f"{label} mismatch: expected {expected!r}, got {value!r}")


def _canonical_json_sha256(value: object) -> str:
    payload = json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _read_json_snapshot(path: pathlib.Path, label: str) -> tuple[bytes, object]:
    """Read once so authentication and parsing cover the identical bytes."""
    if path.is_symlink() or not path.is_file():
        raise SystemExit(f"{path}: {label} must be a regular non-symlink file")
    try:
        payload = path.read_bytes()
    except OSError as exc:
        raise SystemExit(f"{path}: cannot read {label}: {exc}") from exc
    try:
        document = json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise SystemExit(f"{path}: invalid JSON in {label}: {exc}") from exc
    return payload, document


def _result_root_snapshot(
        root: pathlib.Path,
        git_tree_oid: str) -> tuple[dict, list[tuple[pathlib.Path, object]]]:
    result_dir = root / "results"
    files = sorted(result_dir.glob("*.json"))
    data_files = [path for path in files if path.name != "index.json"]
    if not data_files:
        raise SystemExit(f"{result_dir}: no historical result JSON files")
    digest = hashlib.sha256()
    documents = []
    for path in data_files:
        payload, document = _read_json_snapshot(path, "historical result")
        file_sha256 = hashlib.sha256(payload).hexdigest()
        digest.update(f"{file_sha256}  results/{path.name}\n".encode("ascii"))
        documents.append((path, document))
    return (
        {
            "git_tree_oid": git_tree_oid,
            "result_file_count": len(data_files),
            "result_file_manifest_sha256": digest.hexdigest(),
        },
        documents,
    )


def _result_source_snapshots(
        v1: pathlib.Path,
        v2: pathlib.Path,
        template: dict) -> tuple[dict, list[tuple[pathlib.Path, object]]]:
    template = _require_exact_keys(
        template, {"schema", "metadata_status", "round1", "round2"},
        "historical result source template",
    )
    for name in ("round1", "round2"):
        source = _require_exact_keys(
            template[name],
            {"git_tree_oid", "result_file_count", "result_file_manifest_sha256"},
            f"historical result source template.{name}",
        )
        if (type(source["git_tree_oid"]) is not str or
                len(source["git_tree_oid"]) != 40 or
                any(c not in "0123456789abcdef" for c in source["git_tree_oid"])):
            raise SystemExit(f"historical result source template.{name} has invalid Git tree OID")
    round1, round1_documents = _result_root_snapshot(
        v1, template["round1"]["git_tree_oid"],
    )
    round2, round2_documents = _result_root_snapshot(
        v2, template["round2"]["git_tree_oid"],
    )
    return (
        {
            "schema": template["schema"],
            "metadata_status": template["metadata_status"],
            "round1": round1,
            "round2": round2,
        },
        round1_documents + round2_documents,
    )


def _result_source_identities(
        v1: pathlib.Path,
        v2: pathlib.Path,
        template: dict) -> dict:
    identities, _documents = _result_source_snapshots(v1, v2, template)
    return identities


def _authenticated_result_documents(
        v1: pathlib.Path,
        v2: pathlib.Path,
        expected_sources: dict) -> tuple[dict, list[tuple[pathlib.Path, object]]]:
    actual, documents = _result_source_snapshots(v1, v2, expected_sources)
    _require_exact_value(actual, expected_sources, "historical result sources")
    return json.loads(json.dumps(expected_sources)), documents


def _authenticated_timing_document(path: pathlib.Path, expected: dict) -> tuple[dict, object]:
    expected = _require_exact_keys(
        expected, {"git_blob_oid", "sha256"}, "canonical timing source",
    )
    payload, document = _read_json_snapshot(path, "canonical timing input")
    actual = {
        "git_blob_oid": expected["git_blob_oid"],
        "sha256": hashlib.sha256(payload).hexdigest(),
    }
    _require_exact_value(actual, expected, "canonical timing source")
    return json.loads(json.dumps(expected)), document


def validate_historical_diagnostic_cache(
        cache: object,
        expected_methods: set[str] | None = None,
        *,
        expected_sources: dict | None = None) -> dict[str, dict]:
    """Validate and unwrap the deliberately non-selectable historical cache."""
    if expected_sources is None:
        expected_sources = HISTORICAL_RESULT_SOURCES
    top = _require_exact_keys(
        cache,
        {"schema", "classification", "grid", "notice", "sources", "arms_sha256", "arms"},
        "historical diagnostic cache",
    )
    _require_exact_value(top["schema"], HISTORICAL_DIAGNOSTIC_SCHEMA, "cache schema")
    _require_exact_value(
        top["classification"], HISTORICAL_DIAGNOSTIC_CLASSIFICATION,
        "historical diagnostic cache classification",
    )
    _require_exact_value(
        top["grid"], HISTORICAL_DIAGNOSTIC_GRID,
        "historical diagnostic cache grid",
    )
    _require_exact_value(top["notice"], HISTORICAL_DIAGNOSTIC_NOTICE, "cache notice")
    _require_exact_value(top["sources"], expected_sources, "historical diagnostic cache sources")

    arms = top["arms"]
    if type(arms) is not dict or not arms:
        raise SystemExit("historical diagnostic cache arms must be a non-empty JSON object")
    for method, arm in arms.items():
        if type(method) is not str or not method:
            raise SystemExit("historical diagnostic cache has an invalid method name")
        arm = _require_exact_keys(arm, {"eohc", "nc"}, f"historical diagnostic arm {method!r}")
        eohc = arm["eohc"]
        if type(eohc) not in (int, float) or not math.isfinite(eohc) or not 0 <= eohc <= 4:
            raise SystemExit(f"historical diagnostic arm {method!r} has invalid eohc")
        if type(arm["nc"]) is not int or arm["nc"] != HISTORICAL_DIAGNOSTIC_CELLS_PER_ARM:
            raise SystemExit(
                f"historical diagnostic arm {method!r} must have exactly "
                f"{HISTORICAL_DIAGNOSTIC_CELLS_PER_ARM} cells"
            )
    arms_sha256 = _canonical_json_sha256(arms)
    if type(top["arms_sha256"]) is not str or top["arms_sha256"] != arms_sha256:
        raise SystemExit("historical diagnostic cache arms SHA-256 mismatch")
    if expected_sources == HISTORICAL_RESULT_SOURCES and \
            arms_sha256 != HISTORICAL_DIAGNOSTIC_ARMS_SHA256:
        raise SystemExit("historical diagnostic cache payload is not the frozen artifact")
    if expected_methods is not None and set(arms) != expected_methods:
        raise SystemExit(
            "historical diagnostic/result method sets disagree; "
            f"only in results: {sorted(expected_methods - set(arms))}; "
            f"only in diagnostic: {sorted(set(arms) - expected_methods)}"
        )
    return arms


def validate_rows_cache(cache: object) -> list[dict]:
    top = _require_exact_keys(
        cache, {"schema", "sources", "rows_sha256", "rows"}, "rows cache",
    )
    _require_exact_value(top["schema"], ROWS_CACHE_SCHEMA, "rows cache schema")
    _require_exact_value(top["sources"], ROWS_SOURCE_IDENTITIES, "rows cache sources")
    rows = top["rows"]
    if type(rows) is not list or len(rows) != len(REGIONS) * 84:
        raise SystemExit("rows cache must contain exactly 504 rows")
    rows_sha256 = _canonical_json_sha256(rows)
    if type(top["rows_sha256"]) is not str or top["rows_sha256"] != rows_sha256:
        raise SystemExit("rows cache payload SHA-256 mismatch")
    if rows_sha256 != ROWS_PAYLOAD_SHA256:
        raise SystemExit("rows cache payload is not the frozen artifact")
    methods = {row.get("m") for row in rows if type(row) is dict}
    if len(methods) != 84:
        raise SystemExit("rows cache must contain exactly 84 methods")
    for region, _label in REGIONS:
        if sum(type(row) is dict and row.get("reg") == region for row in rows) != 84:
            raise SystemExit(f"rows cache region {region!r} must contain exactly 84 methods")
    return rows


def make_rows_cache(rows: list[dict]) -> dict:
    cache = {
        "schema": ROWS_CACHE_SCHEMA,
        "sources": json.loads(json.dumps(ROWS_SOURCE_IDENTITIES)),
        "rows_sha256": _canonical_json_sha256(rows),
        "rows": rows,
    }
    validate_rows_cache(cache)
    return cache


def historical_diagnostic_overhead(
        v1: pathlib.Path,
        v2: pathlib.Path,
        reg: str,
        *,
        expected_sources: dict | None = None) -> dict:
    """Re-score the historical r1 files on the fixed K=6..100 diagnostic grid.

    This deliberately does not infer a grid by intersecting whatever happened to
    survive in the inputs.  Each arm must contribute exactly ten observations at
    each of the 95 diagnostic K.  Because cells_detail omits seed IDs, satisfying
    that shape does *not* establish paired cell identity; the returned artifact is
    explicitly classified as diagnostic-only and ineligible for selection.
    """
    if reg != HISTORICAL_DIAGNOSTIC_REGION:
        raise SystemExit(
            f"historical diagnostic supports only {HISTORICAL_DIAGNOSTIC_REGION!r}, got {reg!r}"
        )

    if expected_sources is None:
        expected_sources = HISTORICAL_RESULT_SOURCES
    authenticated_sources, result_documents = _authenticated_result_documents(
        v1, v2, expected_sources,
    )

    docs: dict[str, list] = {}
    for path, doc in result_documents:
        if type(doc) is not dict:
            raise SystemExit(f"{path}: historical result must be a JSON object")
        method = doc["method"]
        if type(method) is not str or not method:
            raise SystemExit(f"{path}: invalid method name")
        if method in docs:
            raise SystemExit(f"method {method!r} appears in both result sets; refusing to guess")
        docs[method] = doc["regions"][reg]["cells_detail"]
    if not docs:
        raise SystemExit("historical diagnostic found no result documents")

    arms = {}
    wanted = set(range(HISTORICAL_DIAGNOSTIC_K_MIN, HISTORICAL_DIAGNOSTIC_K_MAX + 1))
    for method, cells_detail in docs.items():
        if type(cells_detail) is not list:
            raise SystemExit(f"method {method!r} cells_detail must be a JSON array")
        counts: Counter[int] = Counter()
        selected = []
        for index, entry in enumerate(cells_detail):
            if type(entry) is not list or len(entry) != 2:
                raise SystemExit(f"method {method!r} cells_detail[{index}] must be [K, outcome]")
            k, outcome = entry
            if type(k) is not int:
                raise SystemExit(f"method {method!r} cells_detail[{index}] has a non-integer K")
            if k < 2 or k > HISTORICAL_DIAGNOSTIC_K_MAX:
                raise SystemExit(f"method {method!r} has unexpected K={k}")
            if type(outcome) is not int or outcome not in (-1, 0, 1, 2, 3, 4):
                raise SystemExit(
                    f"method {method!r} cells_detail[{index}] has invalid outcome {outcome!r}"
                )
            counts[k] += 1
            if k in wanted:
                selected.append(outcome)

        for k in range(HISTORICAL_DIAGNOSTIC_K_MIN, HISTORICAL_DIAGNOSTIC_K_MAX + 1):
            if counts[k] == 0:
                raise SystemExit(f"method {method!r} historical diagnostic is missing K={k}")
            if counts[k] != HISTORICAL_DIAGNOSTIC_CELLS_PER_K:
                raise SystemExit(
                    f"method {method!r} K={k} has multiplicity {counts[k]}, expected "
                    f"{HISTORICAL_DIAGNOSTIC_CELLS_PER_K}"
                )
        # The historical source arms legitimately start at different K below 6.
        # If those optional prefix K are present, their multiplicity must still
        # match the campaign's ten-cell construction-seed count.
        for k in range(2, HISTORICAL_DIAGNOSTIC_K_MIN):
            if counts[k] not in (0, HISTORICAL_DIAGNOSTIC_CELLS_PER_K):
                raise SystemExit(
                    f"method {method!r} K={k} has multiplicity {counts[k]}, expected 0 or "
                    f"{HISTORICAL_DIAGNOSTIC_CELLS_PER_K}"
                )
        if len(selected) != HISTORICAL_DIAGNOSTIC_CELLS_PER_ARM:
            raise SystemExit(
                f"method {method!r} historical diagnostic has {len(selected)} cells, expected "
                f"{HISTORICAL_DIAGNOSTIC_CELLS_PER_ARM}"
            )
        # -1 is the censor sentinel (dead / construct failure); the metric charges
        # it the absolute cap of 4, exactly as the published aggregate does.
        arms[method] = {
            "eohc": sum(4 if outcome == -1 else outcome for outcome in selected) / len(selected),
            "nc": len(selected),
        }

    cache = {
        "schema": HISTORICAL_DIAGNOSTIC_SCHEMA,
        "classification": dict(HISTORICAL_DIAGNOSTIC_CLASSIFICATION),
        "grid": dict(HISTORICAL_DIAGNOSTIC_GRID),
        "notice": HISTORICAL_DIAGNOSTIC_NOTICE,
        "sources": authenticated_sources,
        "arms_sha256": _canonical_json_sha256(arms),
        "arms": arms,
    }
    validate_historical_diagnostic_cache(
        cache, set(docs), expected_sources=expected_sources,
    )
    return cache


def compact(rows: list[dict], desc: dict[str, str], alt: dict | None = None) -> list[dict]:
    """Trim rows to what the page draws, and attach the display group."""
    alt_arms = None
    if alt is not None:
        alt_arms = validate_historical_diagnostic_cache(alt, {r["m"] for r in rows})
    out = []
    for r in rows:
        g = FAMILY_GROUP.get(r["fam"])
        if g is None:
            raise SystemExit(f"no display group for family {r['fam']!r} (method {r['m']})")
        row = {
            "m": r["m"],
            "g": g,
            "reg": r["reg"],
            "eoh": r["eoh"],
            "sps": round(r["sps"]),
            "f0": r["f0"],
            "f2": r["f2"],
            "dead": r["dead"],
            "cf": r["cf"],
            "hk": r["hk"],
            "n": r["n"],
            "exc": r["exc"],
        }
        if alt_arms is not None and r["reg"] == HISTORICAL_DIAGNOSTIC_REGION:
            row["eohc"] = alt_arms[r["m"]]["eohc"]
            row["nc"] = alt_arms[r["m"]]["nc"]
        out.append(row)
    unknown = sorted({r["m"] for r in out} - set(desc))
    if unknown:
        print(f"warning: no description for {len(unknown)} method(s): {unknown[:6]}", file=sys.stderr)
    return out


def build_summary(rows: list[dict]) -> list[dict]:
    """Per region: the fastest WH2 arm measured against Wirehair 1."""
    tiles = []
    for reg, label in REGIONS:
        rr = [r for r in rows if r["reg"] == reg]
        wh1 = next(r for r in rr if r["m"] == WH1)
        best = max((r for r in rr if r["m"] != WH1), key=lambda r: r["sps"])
        tiles.append(
            {
                "reg": reg,
                "label": label,
                "pct": 100.0 * best["sps"] / wh1["sps"],
                "best": best["m"],
                "bestSps": best["sps"],
                "wh1Sps": wh1["sps"],
            }
        )
    return tiles


# Standalone pages get a real document shell -- without a viewport meta a phone
# renders the whole thing at desktop width. The --fragment form omits it because the
# Artifact host supplies its own <head>, and a nested doctype there would be invalid.
SHELL = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
{head}
</head>
<body>
{body}
</body>
</html>
"""


def render(rows: list[dict], desc: dict[str, str], prose: dict[str, str]) -> str:
    data = json.dumps(rows, separators=(",", ":"))
    groups = json.dumps({k: {"label": v[0], "hue": v[1], "shape": v[2]} for k, v in GROUPS.items()},
                        separators=(",", ":"))
    regions = json.dumps(REGIONS, separators=(",", ":"))
    descs = json.dumps(desc, separators=(",", ":"), ensure_ascii=False)
    tiles = json.dumps(build_summary(rows), separators=(",", ":"))
    return TEMPLATE.format(
        DATA=data, GROUPS=groups, REGIONS=regions, DESC=descs, TILES=tiles,
        WH1=json.dumps(WH1), SELECTED=json.dumps(SELECTED), **prose,
    )


TEMPLATE = r"""<title>{TITLE}</title>
<style>
:root{{
  color-scheme:light;
  --plane:#f9f9f7; --surface:#fcfcfb; --sunk:#f2f1ec;
  --ink:#0b0b0b; --ink2:#52514e; --ink3:#898781;
  --grid:#e1e0d9; --axis:#c3c2b7; --hair:rgba(11,11,11,.10);
  --s1:#2a78d6; --s2:#008300; --s3:#e87ba4; --s4:#eda100;
  --acc:#2a78d6; --ring:#fcfcfb;
}}
@media (prefers-color-scheme:dark){{
  :root:where(:not([data-theme="light"])){{
    color-scheme:dark;
    --plane:#0d0d0d; --surface:#1a1a19; --sunk:#232321;
    --ink:#fff; --ink2:#c3c2b7; --ink3:#898781;
    --grid:#2c2c2a; --axis:#383835; --hair:rgba(255,255,255,.10);
    --s1:#3987e5; --s2:#008300; --s3:#d55181; --s4:#c98500;
    --acc:#3987e5; --ring:#1a1a19;
  }}
}}
:root[data-theme="dark"]{{
  color-scheme:dark;
  --plane:#0d0d0d; --surface:#1a1a19; --sunk:#232321;
  --ink:#fff; --ink2:#c3c2b7; --ink3:#898781;
  --grid:#2c2c2a; --axis:#383835; --hair:rgba(255,255,255,.10);
  --s1:#3987e5; --s2:#008300; --s3:#d55181; --s4:#c98500;
  --acc:#3987e5; --ring:#1a1a19;
}}
*{{box-sizing:border-box;}}
body{{background:var(--plane);color:var(--ink);margin:0;padding:0;
  font:15px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif;
  -webkit-font-smoothing:antialiased;}}
.wrap{{max-width:1320px;margin:0 auto;padding:34px 26px 60px;display:flex;flex-direction:column;gap:26px;}}
.mono,code{{font-family:ui-monospace,"SF Mono",Menlo,Consolas,monospace;}}

header{{display:flex;flex-direction:column;gap:9px;}}
.eyebrow{{margin:0;font-size:11.5px;letter-spacing:.1em;text-transform:uppercase;
  color:var(--ink3);font-weight:600;}}
h1{{margin:0;font-size:27px;line-height:1.2;font-weight:640;letter-spacing:-.015em;text-wrap:balance;}}
.sub{{margin:0;color:var(--ink2);font-size:14.5px;max-width:76ch;}}
.readnote{{margin:0;color:var(--ink3);font-size:13px;max-width:82ch;}}

.tiles{{display:grid;grid-template-columns:repeat(auto-fit,minmax(168px,1fr));gap:10px;}}
.tile{{background:var(--surface);border:1px solid var(--hair);border-radius:9px;
  padding:12px 14px 13px;display:flex;flex-direction:column;gap:3px;}}
.tile .tl{{font-size:11.5px;color:var(--ink3);font-weight:600;letter-spacing:.04em;}}
.tile .tv{{font-size:26px;font-weight:640;letter-spacing:-.02em;line-height:1.1;}}
.tile .td{{font-size:11.5px;color:var(--ink2);}}
.tile .tm{{font-size:11px;color:var(--ink3);font-family:ui-monospace,Menlo,monospace;
  overflow:hidden;text-overflow:ellipsis;white-space:nowrap;}}
.tilehead{{display:flex;align-items:baseline;gap:10px;flex-wrap:wrap;}}
.tilehead h2{{margin:0;font-size:15px;font-weight:620;}}
.tilehead p{{margin:0;font-size:12.5px;color:var(--ink3);}}

.controls{{display:flex;gap:8px;align-items:center;flex-wrap:wrap;}}
button{{background:var(--surface);color:var(--ink2);border:1px solid var(--hair);
  border-radius:7px;padding:6px 13px;font:13px system-ui;cursor:pointer;}}
button:hover{{color:var(--ink);}}
button[aria-pressed="true"]{{border-color:var(--acc);color:var(--acc);font-weight:600;}}
button:focus-visible,summary:focus-visible,a:focus-visible{{outline:2px solid var(--acc);outline-offset:2px;}}
.hint{{color:var(--ink3);font-size:12.5px;margin-left:2px;}}

.legend{{display:flex;flex-wrap:wrap;gap:7px;}}
.lg{{display:inline-flex;gap:7px;align-items:center;font-size:12.5px;color:var(--ink2);
  background:var(--surface);border:1px solid var(--hair);border-radius:999px;
  padding:4px 12px 4px 8px;cursor:pointer;}}
.lg[aria-pressed="false"]{{opacity:.4;}}
.lg svg{{flex:none;}}

.facets{{display:grid;grid-template-columns:repeat(auto-fit,minmax(392px,1fr));gap:14px;}}
.facet{{background:var(--surface);border:1px solid var(--hair);border-radius:10px;padding:13px 12px 8px;}}
.fhead{{font-size:13px;font-weight:640;margin:0 0 1px 4px;letter-spacing:-.01em;}}
.fsub{{font-size:11.5px;color:var(--ink3);margin:0 0 2px 4px;font-variant-numeric:tabular-nums;}}
.facet svg{{display:block;width:100%;height:auto;touch-action:none;}}
svg text{{font-family:system-ui;fill:var(--ink3);font-size:10.5px;font-variant-numeric:tabular-nums;}}
svg text.dl{{fill:var(--ink2);font-size:10.5px;font-weight:650;}}
svg text.axt{{fill:var(--ink3);font-size:10.5px;}}
.dot{{cursor:pointer;}}

.tt{{position:fixed;pointer-events:none;background:var(--surface);border:1px solid var(--hair);
  border-radius:8px;padding:9px 12px;font-size:12.5px;box-shadow:0 4px 18px rgba(0,0,0,.19);
  display:none;z-index:20;max-width:310px;}}
.tt .ttm{{font-family:ui-monospace,Menlo,monospace;font-size:12px;font-weight:650;
  display:flex;align-items:center;gap:6px;}}
.tt .ttg{{color:var(--ink3);font-size:11.5px;margin-top:1px;}}
.tt .ttd{{color:var(--ink2);font-size:11.5px;line-height:1.45;margin-top:5px;}}
.tt table{{border-collapse:collapse;margin-top:6px;width:100%;}}
.tt td{{padding:1.5px 0;font-variant-numeric:tabular-nums;}}
.tt td:first-child{{color:var(--ink3);padding-right:12px;}}
.tt td:last-child{{text-align:right;font-weight:600;}}

.panel{{background:var(--surface);border:1px solid var(--hair);border-radius:10px;
  padding:6px 16px;overflow-x:auto;}}
#tblwrap{{display:none;}}
table.data{{border-collapse:collapse;font-size:12.5px;min-width:940px;width:100%;}}
table.data th{{text-align:left;color:var(--ink2);font-weight:620;padding:9px 14px 8px 0;
  border-bottom:1px solid var(--grid);white-space:nowrap;position:sticky;top:0;background:var(--surface);}}
table.data th[data-k]{{cursor:pointer;}}
table.data td{{padding:5px 14px 5px 0;border-bottom:1px solid var(--grid);
  font-variant-numeric:tabular-nums;white-space:nowrap;}}
table.data td.mono{{font-family:ui-monospace,Menlo,monospace;font-size:11.5px;}}
table.data tr:last-child td{{border-bottom:none;}}
.swatch{{display:inline-block;width:9px;height:9px;border-radius:2px;margin-right:6px;
  vertical-align:baseline;}}

details.gd{{background:var(--surface);border:1px solid var(--hair);border-radius:10px;padding:12px 16px;}}
details.gd+details.gd{{margin-top:9px;}}
details.gd summary{{cursor:pointer;font-size:14px;font-weight:620;}}
details.gd[open] summary{{margin-bottom:8px;}}
.gd p{{color:var(--ink2);font-size:13.5px;max-width:92ch;margin:9px 0 0;}}
.gd ul{{color:var(--ink2);font-size:13.5px;max-width:92ch;margin:8px 0 2px;padding-left:19px;}}
.gd li{{margin:6px 0;}}
.gd dl{{margin:6px 0 2px;}}
.gd dt{{margin-top:11px;font-size:12.5px;font-family:ui-monospace,Menlo,monospace;color:var(--ink);}}
.gd dd{{margin:2px 0 0;color:var(--ink2);font-size:12.5px;max-width:92ch;}}
.sec{{display:flex;flex-direction:column;gap:9px;}}
.sech{{margin:0;font-size:15px;font-weight:620;}}
footer{{color:var(--ink3);font-size:11.5px;max-width:104ch;border-top:1px solid var(--hair);
  padding-top:14px;line-height:1.6;}}
@media (max-width:560px){{ .wrap{{padding:22px 15px 44px;}} h1{{font-size:22px;}} }}
</style>

<div class="wrap">
<header>
  <p class="eyebrow">Untuned architecture comparison · round 3</p>
  <h1>{TITLE}</h1>
  <p class="sub">{SUBTITLE}</p>
  <p class="readnote">{READNOTE}</p>
</header>

<section class="sec">
  <div class="tilehead">
    <h2>Fastest WH2 arm as a share of Wirehair 1</h2>
    <p>encoder throughput, same host and build for both</p>
  </div>
  <div class="tiles" id="tiles"></div>
  <p class="readnote">{TILENOTE}</p>
</section>

<div class="controls">
  <button id="btnChart" aria-pressed="true">Chart</button>
  <button id="btnTable" aria-pressed="false">Table</button>
  <button id="btnHistorical" aria-pressed="false" aria-describedby="historicalNote" title="Historical K=6–100 diagnostic only; non-selectable and never architecture evidence">r1: historical diagnostic</button>
  <span class="hint">The table carries every plotted value — it is the keyboard and screen-reader path.</span>
</div>
<p class="readnote" id="historicalNote">{HISTORICALNOTE}</p>

<div class="legend" id="legend"></div>
<div class="facets" id="facets"></div>
<div class="panel" id="tblwrap"></div>

<section class="sec">
  <h2 class="sech">What the run shows</h2>
  <div id="findings"></div>
</section>

<section class="sec">
  <h2 class="sech">Configurations</h2>
  <div id="glossary"></div>
</section>

<footer>{FOOTER}</footer>
</div>
<div class="tt" id="tt" role="tooltip"></div>

<script>
const DATA={DATA};
const GROUPS={GROUPS};
const REGIONS={REGIONS};
const DESC={DESC};
const TILES={TILES};
const WH1={WH1}, SELECTED={SELECTED};
const FINDINGS={FINDINGS};
const PANELS={PANELS};

const W=460,H=344,ML=52,MR=13,MT=11,MB=40,BAND=27;
const plotW=W-ML-MR, plotH=H-MT-MB-BAND;
const XMIN=1.0e6,XMAX=8.4e6,YMIN=1e-4,YMAX=0.35;
const xs=v=>ML+(v-XMIN)/(XMAX-XMIN)*plotW;
const ys=v=>MT+(Math.log10(YMAX)-Math.log10(Math.max(v,YMIN)))/(Math.log10(YMAX)-Math.log10(YMIN))*plotH;
const bandY=MT+plotH+11, bandMid=bandY+(BAND-8)/2;
/* OH() can display the fixed historical K=6..100 diagnostic.  The old cells_detail
   schema has no seed IDs, so this view is non-selectable and never architecture
   evidence even though every arm passed the 10-cells-per-K shape check. */
let historicalDiagnostic=false;
const OH=d=>(historicalDiagnostic&&d.eohc!=null)?d.eohc:d.eoh;
const PY=d=>OH(d)===0?bandMid:ys(OH(d));
const fmt=(v,d=4)=>v==null?"–":(+v).toFixed(d).replace(/0+$/,"").replace(/\.$/,"");
const M=v=>(v/1e6).toFixed(2)+"M";
const esc=s=>String(s).replace(/[&<>"]/g,c=>({{"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}})[c]);

/* Marks are painted with `var(--token)` rather than a resolved hex so the viewer's
   theme toggle repaints them; resolving at build time freezes the first theme. */
function markSVG(g,x,y,r,cls,extra){{
  x=+x; y=+y;  /* coerce: the diamond/square arms are arithmetic, not concatenation */
  const f=`var(${{GROUPS[g].hue}})`, sh=GROUPS[g].shape;
  const a=`fill="${{f}}" stroke="var(--ring)" stroke-width="2" ${{cls?`class="${{cls}}"`:""}} ${{extra||""}}`;
  const n=v=>(+v).toFixed(1);
  if(sh==="c") return `<circle cx="${{n(x)}}" cy="${{n(y)}}" r="${{r}}" ${{a}}></circle>`;
  if(sh==="d"){{const k=r*1.28;
    return `<path d="M ${{n(x)}} ${{n(y-k)}} L ${{n(x+k)}} ${{n(y)}} L ${{n(x)}} ${{n(y+k)}} L ${{n(x-k)}} ${{n(y)}} Z" ${{a}}></path>`;}}
  if(sh==="s"){{const h=r*1.08;
    return `<rect x="${{n(x-h)}}" y="${{n(y-h)}}" width="${{n(h*2)}}" height="${{n(h*2)}}" rx="1.5" ${{a}}></rect>`;}}
  /* "r" — the selected profile: a ring, so it reads as chosen rather than as one more dot */
  return `<circle cx="${{n(x)}}" cy="${{n(y)}}" r="${{r}}" fill="var(--surface)" stroke="var(${{GROUPS[g].hue}})" stroke-width="3.4" ${{cls?`class="${{cls}}"`:""}} ${{extra||""}}></circle>`;
}}

function facetSVG(reg,rows){{
  let s=`<svg viewBox="0 0 ${{W}} ${{H}}" data-reg="${{reg}}" role="img" aria-label="Expected overhead versus encoder throughput, ${{reg}}">`;
  for(const e of [-4,-3,-2,-1]){{const y=ys(Math.pow(10,e));
    s+=`<line x1="${{ML}}" x2="${{W-MR}}" y1="${{y}}" y2="${{y}}" stroke="var(--grid)" stroke-width="1"></line>`;
    s+=`<text x="${{ML-7}}" y="${{y+3.5}}" text-anchor="end">1e${{e}}</text>`;}}
  for(let m=2;m<=8;m+=2){{const x=xs(m*1e6);
    s+=`<line x1="${{x}}" x2="${{x}}" y1="${{MT}}" y2="${{MT+plotH}}" stroke="var(--grid)" stroke-width="1"></line>`;
    s+=`<text class="axt" x="${{x}}" y="${{H-16}}" text-anchor="middle">${{m}}M</text>`;}}
  s+=`<line x1="${{ML}}" x2="${{W-MR}}" y1="${{MT+plotH}}" y2="${{MT+plotH}}" stroke="var(--axis)" stroke-width="1"></line>`;
  /* Exact-zero arms cannot sit on a log axis, so they get their own band below it.
     The band is labelled in the y gutter like any other tick -- an in-band caption
     collides with the marks that land there. */
  s+=`<rect x="${{ML}}" y="${{bandY}}" width="${{plotW}}" height="${{BAND-8}}" fill="var(--sunk)" rx="4"></rect>`;
  s+=`<text x="${{ML-7}}" y="${{bandMid+3.5}}" text-anchor="end">0</text>`;
  s+=`<text class="axt" x="${{ML+plotW/2}}" y="${{H-3}}" text-anchor="middle">encoder throughput (symbols / s)</text>`;
  s+=`<text class="axt" x="${{-(MT+plotH/2)}}" y="12" transform="rotate(-90)" text-anchor="middle">expected overhead</text>`;

  const pts=rows.map(d=>({{d,x:xs(d.sps),y:PY(d)}}));
  /* draw the crowd first, the two labelled entities last so they are never buried */
  const rank=d=>d.m===SELECTED?2:d.m===WH1?1:0;
  pts.sort((a,b)=>rank(a.d)-rank(b.d)||OH(b.d)-OH(a.d));
  for(const p of pts){{
    const big=p.d.m===SELECTED||p.d.m===WH1;
    s+=markSVG(p.d.g,p.x.toFixed(1),p.y.toFixed(1),big?6.4:4.6,"dot",
      `data-m="${{esc(p.d.m)}}" data-reg="${{reg}}"`);
  }}
  for(const nm of [WH1,SELECTED]){{
    const p=pts.find(q=>q.d.m===nm); if(!p) continue;
    const label=nm===WH1?"Wirehair 1":"dispatch v1";
    let tx=p.x+10, anchor="start";
    if(tx>W-MR-62){{tx=p.x-10;anchor="end";}}
    s+=`<text class="dl" x="${{tx.toFixed(1)}}" y="${{(p.y-9).toFixed(1)}}" text-anchor="${{anchor}}">${{label}}</text>`;
  }}
  return s+`</svg>`;
}}

const facets=document.getElementById("facets");
function renderFacets(){{
  facets.textContent="";
  for(const [reg,name] of REGIONS){{
    const rows=DATA.filter(d=>d.reg===reg);
    const zero=rows.filter(d=>OH(d)===0).length;
    const alt=historicalDiagnostic&&rows.some(d=>d.eohc!=null);
    const cells=alt?rows.find(d=>d.nc!=null).nc:rows[0].n;
    const div=document.createElement("div"); div.className="facet";
    div.innerHTML=`<p class="fhead">${{name}}</p><p class="fsub">${{cells.toLocaleString()}} cells per point`+
      `${{alt?" (historical K=6–100 diagnostic)":""}} · ${{zero}} of ${{rows.length}} arms never needed a repair symbol`+
      `${{(!alt&&rows.some(d=>d.eohc!=null))?" · 12 arms omit the smallest K":""}}</p>${{facetSVG(reg,rows)}}`;
    facets.appendChild(div);
  }}
}}
renderFacets();

/* legend — real buttons, so groups can be isolated from the keyboard too */
const legend=document.getElementById("legend");
const off=new Set();
/* re-applied after every facet re-render so an isolated group survives a redraw */
function applyGroupFilter(){{
  const byName={{}};
  for(const r of DATA) byName[r.m]=r.g;
  document.querySelectorAll(".dot").forEach(el=>{{
    el.style.display=off.has(byName[el.dataset.m])?"none":"";
  }});
}}
for(const [g,spec] of Object.entries(GROUPS)){{
  const b=document.createElement("button");
  b.className="lg"; b.type="button"; b.setAttribute("aria-pressed","true"); b.dataset.g=g;
  const c=`var(${{spec.hue}})`;
  const sh=spec.shape==="c"?`<circle cx="7" cy="7" r="4.6" fill="${{c}}"></circle>`
    :spec.shape==="d"?`<path d="M 7 1.9 L 12.1 7 L 7 12.1 L 1.9 7 Z" fill="${{c}}"></path>`
    :spec.shape==="s"?`<rect x="2.4" y="2.4" width="9.2" height="9.2" rx="1.5" fill="${{c}}"></rect>`
    :`<circle cx="7" cy="7" r="4.4" fill="var(--surface)" stroke="${{c}}" stroke-width="3"></circle>`;
  b.innerHTML=`<svg width="14" height="14" aria-hidden="true">${{sh}}</svg><span></span>`;
  b.querySelector("span").textContent=spec.label;
  b.addEventListener("click",()=>{{
    off.has(g)?off.delete(g):off.add(g);
    b.setAttribute("aria-pressed",off.has(g)?"false":"true");
    applyGroupFilter();
  }});
  legend.appendChild(b);
}}

/* summary tiles */
const tiles=document.getElementById("tiles");
for(const t of TILES){{
  const d=document.createElement("div"); d.className="tile";
  d.innerHTML=`<span class="tl"></span><span class="tv"></span><span class="td"></span><span class="tm"></span>`;
  d.querySelector(".tl").textContent=t.label;
  d.querySelector(".tv").textContent=t.pct.toFixed(0)+"%";
  d.querySelector(".td").textContent=M(t.bestSps)+" vs "+M(t.wh1Sps)+" sym/s";
  d.querySelector(".tm").textContent=t.best;
  tiles.appendChild(d);
}}

/* hover — nearest point within the facet, so the pointer only has to be closest */
const tt=document.getElementById("tt");
let active=null;
function show(row,cx,cy){{
  const g=GROUPS[row.g];
  const historicalCell=historicalDiagnostic&&row.eohc!=null;
  tt.innerHTML=`<div class="ttm"><svg width="11" height="11" aria-hidden="true"><circle cx="5.5" cy="5.5" r="4.4" fill="var(${{g.hue}})"></circle></svg><span></span></div>`+
    `<div class="ttg"></div><div class="ttd"></div>`+
    `<table><tr><td>expected overhead</td><td>${{OH(row)===0?"0 in sample":fmt(OH(row),5)}}</td></tr>`+
    (row.eohc!=null?`<tr><td>${{historicalDiagnostic?"as published":"historical diagnostic"}}</td><td>`+
      `${{fmt(historicalDiagnostic?row.eoh:row.eohc,5)}}</td></tr>`:"")+
    `<tr><td>throughput</td><td>${{M(row.sps)}} sym/s</td></tr>`+
    (historicalCell
      ? `<tr><td colspan="2">Failure, dead, hard-K and construction metrics hidden: they exist only for the published r1 domain.</td></tr>`
      : `<tr><td>fail @ +0 / +2</td><td>${{fmt(row.f0,5)}} / ${{fmt(row.f2,5)}}</td></tr>`+
        `<tr><td>dead @ cap / hard-K</td><td>${{fmt(row.dead,5)}} / ${{fmt(row.hk,5)}}</td></tr>`+
        `<tr><td>construct failures</td><td>${{fmt(row.cf,5)}}</td></tr>`)+`</table>`;
  tt.querySelector(".ttm span").textContent=row.m;
  tt.querySelector(".ttg").textContent=g.label;
  tt.querySelector(".ttd").textContent=DESC[row.m]||"";
  tt.style.display="block";
  const r=tt.getBoundingClientRect();
  tt.style.left=Math.min(cx+15,innerWidth-r.width-10)+"px";
  tt.style.top=Math.min(cy+13,innerHeight-r.height-10)+"px";
  document.querySelectorAll(".dot").forEach(el=>{{el.style.opacity=el.dataset.m===row.m?1:.22;}});
}}
function clear(){{
  tt.style.display="none"; active=null;
  document.querySelectorAll(".dot").forEach(el=>{{el.style.opacity=1;}});
}}
facets.addEventListener("pointermove",e=>{{
  const svg=e.target.closest("svg"); if(!svg){{clear();return;}}
  const reg=svg.dataset.reg, box=svg.getBoundingClientRect();
  const sx=(e.clientX-box.left)/box.width*W, sy=(e.clientY-box.top)/box.height*H;
  let best=null,bd=1e9;
  for(const d of DATA){{
    if(d.reg!==reg||off.has(d.g)) continue;
    const px=xs(d.sps), py=PY(d);
    const dist=(px-sx)**2+(py-sy)**2;
    if(dist<bd){{bd=dist;best=d;}}
  }}
  /* generous catch radius in viewBox units — an 9px mark is not a hit target */
  if(best&&bd<=18*18){{ if(!active||active.m!==best.m||active.reg!==best.reg){{active=best;}} show(best,e.clientX,e.clientY); }}
  else clear();
}});
facets.addEventListener("pointerleave",clear);

/* table view — the complete, non-hover path to every plotted value */
const tw=document.getElementById("tblwrap");
let sortKey="eoh",sortDir=1;
function buildTable(){{
  const historicalCell=r=>historicalDiagnostic&&r.eohc!=null;
  const VAL=(r,k)=>{{
    if(k==="eoh") return OH(r);
    if(historicalCell(r)&&(k==="f0"||k==="f2"||k==="dead"||k==="cf")) return Infinity;
    return r[k];
  }};
  const rows=[...DATA].sort((a,b)=>
    a.reg.localeCompare(b.reg)||(VAL(a,sortKey)-VAL(b,sortKey))*sortDir||a.m.localeCompare(b.m));
  let h=`<table class="data"><thead><tr><th>Region</th><th>Configuration</th><th>Group</th>`+
    `<th data-k="eoh">Expected overhead</th><th data-k="sps">sym/s</th>`+
    `<th data-k="f0">fail @ +0 (published)</th><th data-k="f2">fail @ +2 (published)</th>`+
    `<th data-k="dead">dead (published)</th><th data-k="cf">construct fail (published)</th><th data-k="n">cells</th></tr></thead><tbody>`;
  for(const r of rows){{
    const hidden=historicalCell(r);
    h+=`<tr><td>${{r.reg.split("_")[0]}}</td><td class="mono">${{esc(r.m)}}</td>`+
      `<td><span class="swatch" style="background:var(${{GROUPS[r.g].hue}})"></span>${{esc(GROUPS[r.g].label)}}</td>`+
      `<td>${{OH(r)===0?"0":fmt(OH(r),5)}}</td><td>${{M(r.sps)}}</td><td>${{hidden?"—":fmt(r.f0,5)}}</td>`+
      `<td>${{hidden?"—":fmt(r.f2,5)}}</td><td>${{hidden?"—":fmt(r.dead,5)}}</td><td>${{hidden?"—":fmt(r.cf,5)}}</td>`+
      `<td>${{(historicalDiagnostic&&r.nc!=null?r.nc:r.n).toLocaleString()}}</td></tr>`;
  }}
  tw.innerHTML=h+"</tbody></table>";
  tw.querySelectorAll("th[data-k]").forEach(th=>th.addEventListener("click",()=>{{
    const k=th.dataset.k; sortDir=(k===sortKey)?-sortDir:1; sortKey=k; buildTable();
  }}));
}}
buildTable();

/* findings + protocol/caveat panels */
const fw=document.getElementById("findings");
fw.innerHTML=`<details class="gd" open><summary>Findings</summary><ul>`+
  FINDINGS.map(t=>`<li>${{t}}</li>`).join("")+`</ul></details>`;
for(const [title,body] of PANELS)
  fw.insertAdjacentHTML("beforeend",`<details class="gd"><summary>${{title}}</summary>${{body}}</details>`);

/* per-group configuration glossary */
(function(){{
  const g=document.getElementById("glossary"), by={{}};
  for(const r of DATA) (by[r.g]=by[r.g]||new Set()).add(r.m);
  for(const [grp,spec] of Object.entries(GROUPS)){{
    const ms=[...(by[grp]||[])].sort(); if(!ms.length) continue;
    let inner="";
    for(const m of ms) inner+=`<dt>${{esc(m)}}</dt><dd>${{esc(DESC[m]||"")}}</dd>`;
    g.insertAdjacentHTML("beforeend",
      `<details class="gd"><summary><span class="swatch" style="background:var(${{spec.hue}})"></span>`+
      `${{esc(spec.label)}} <span style="color:var(--ink3);font-weight:400">(${{ms.length}})</span></summary>`+
      `<dl>${{inner}}</dl></details>`);
  }}
}})();

const bC=document.getElementById("btnChart"),bT=document.getElementById("btnTable");
const view=t=>{{
  facets.style.display=t?"none":"grid";
  legend.style.display=t?"none":"flex";
  tw.style.display=t?"block":"none";
  bT.setAttribute("aria-pressed",t?"true":"false");
  bC.setAttribute("aria-pressed",t?"false":"true");
  if(t) clear();
}};
bC.onclick=()=>view(false); bT.onclick=()=>view(true);

const bHistorical=document.getElementById("btnHistorical");
bHistorical.onclick=()=>{{
  historicalDiagnostic=!historicalDiagnostic;
  bHistorical.setAttribute("aria-pressed",historicalDiagnostic?"true":"false");
  clear(); renderFacets(); applyGroupFilter(); buildTable();
}};
</script>
"""


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--v1", type=pathlib.Path, help="round-1 sealed result directory")
    ap.add_argument("--v2", type=pathlib.Path, help="round-2 sealed result directory")
    ap.add_argument("--timing", type=pathlib.Path, help="canonical quiet timing JSON (x-axis)")
    ap.add_argument("--rows", type=pathlib.Path, help="previously extracted rows JSON")
    ap.add_argument("--rows-out", type=pathlib.Path, help="write the extracted rows here")
    ap.add_argument(
        "--alt", type=pathlib.Path,
        help="previously computed non-selectable historical r1 diagnostic cache",
    )
    ap.add_argument(
        "--alt-out", type=pathlib.Path,
        help="write the non-selectable historical r1 diagnostic cache here",
    )
    ap.add_argument("--desc", type=pathlib.Path, required=True, help="method description JSON")
    ap.add_argument("--prose", type=pathlib.Path, required=True, help="page prose JSON")
    ap.add_argument("--out", type=pathlib.Path, required=True, help="HTML to write")
    ap.add_argument("--fragment", action="store_true",
                    help="emit body-only HTML for a host that supplies its own document shell")
    a = ap.parse_args()

    if a.rows:
        rows_cache = json.loads(a.rows.read_text())
        rows = validate_rows_cache(rows_cache)
    else:
        if not (a.v1 and a.v2 and a.timing):
            ap.error("either --rows, or all of --v1/--v2/--timing")
        rows = extract_rows(a.v1, a.v2, a.timing)
        rows_cache = make_rows_cache(rows)
    if a.rows_out:
        a.rows_out.write_text(json.dumps(rows_cache, indent=0, sort_keys=True))

    desc = json.loads(a.desc.read_text())
    prose = json.loads(a.prose.read_text())
    prose["FINDINGS"] = json.dumps(prose["FINDINGS"], ensure_ascii=False)
    prose["PANELS"] = json.dumps(prose["PANELS"], ensure_ascii=False)

    alt = None
    if a.v1 and a.v2:
        alt = historical_diagnostic_overhead(a.v1, a.v2, HISTORICAL_DIAGNOSTIC_REGION)
        if a.alt_out:
            a.alt_out.write_text(json.dumps(alt, indent=0, sort_keys=True))
    elif a.alt:
        alt = json.loads(a.alt.read_text())

    compacted = compact(rows, desc, alt)
    page = render(compacted, desc, prose)
    if not a.fragment:
        # split the <title>/<style> preamble from the body so each lands in the
        # right half of the document shell
        cut = page.index("</style>") + len("</style>")
        page = SHELL.format(head=page[:cut].strip(), body=page[cut:].strip())
    a.out.write_text(page)
    n_m = len({r["m"] for r in compacted})
    print(f"wrote {a.out} — {len(compacted)} rows, {n_m} configurations, {len(REGIONS)} regions")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
