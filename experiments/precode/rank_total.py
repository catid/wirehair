#!/usr/bin/env python3
"""Validate and rank precode_sim results by calibrated cost and reliability.

The cost ledger is expressed in block-XOR equivalents:

  base = recv_xors_per_packet * (K + oh)
       + precode_gen_xors_mu + sparse_solve_xors_mu + backsub_xors_mu
       + (ge_real_rowops_mu when recorded, otherwise ge_block_xors_mu)
  total = base + muladd:xor_ratio * heavy_muladds_mu

Reliability at a selected heavy-row count is recovered from ``def_pdf``.
The one-sided confidence gate uses the Wilson score upper bound.  By default
cost remains the complete ledger measured at each row's modeled H.  With
``--cost-h reliability``, only the analytic heavy term is recomputed as
``H * inact_mu + H^2``; the remaining ledger is explicitly held at modeled H.

Parsing is atomic: every result and calibration source is strictly validated
before either the human report or optional deterministic JSON is emitted.
"""

import argparse
import csv
import io
import json
import math
import os
import re
import statistics
import subprocess
import sys
import tempfile


DEFAULT_BB = "1280,102400,1048576"
H_PRIMES = (6, 8, 12, 16)
MAX_UINT32 = (1 << 32) - 1
MAX_TRIALS = 1000000
DISPLAY_TOLERANCE = 0.500001e-6

BASE_FIELDS = (
    "K", "scheme", "D", "H", "oh", "trials", "fail_rate",
    "fail_rate_noheavy", "def_mu", "def_max", "def_pdf", "inact_mu",
    "inact_sd", "inact_max", "rank_mu", "peeled_mu", "residual_rows_mu",
    "recv_xors_per_packet", "precode_gen_xors_mu",
    "sparse_solve_xors_mu", "backsub_xors_mu", "ge_block_xors_mu",
    "ge_bitops_mu",
)
HEAVY_FIELDS = ("heavy_muladds_mu", "heavy_divs_mu")
REPLAY_FIELDS = (
    "ge_real_bitops_mu", "ge_real_rowops_mu", "fill_in_mu",
    "def_outside_w18_rate", "def_band_w95", "def_band_w99",
)
RESULT_SCHEMAS = {
    BASE_FIELDS,
    BASE_FIELDS + HEAVY_FIELDS,
    BASE_FIELDS + HEAVY_FIELDS + ("runaway_rate",),
    BASE_FIELDS + HEAVY_FIELDS + REPLAY_FIELDS,
    BASE_FIELDS + HEAVY_FIELDS + REPLAY_FIELDS + ("runaway_rate",),
    ("pivot_window",) + BASE_FIELDS + HEAVY_FIELDS + REPLAY_FIELDS
    + ("runaway_rate",),
}

INTEGER_FIELDS = {
    "pivot_window", "K", "D", "H", "oh", "trials", "def_max",
    "inact_max", "def_band_w95", "def_band_w99",
}
RATE_FIELDS = {
    "fail_rate", "fail_rate_noheavy", "def_outside_w18_rate",
    "runaway_rate",
}

XOR_CALIBRATION_FIELDS = (
    "block_bytes", "working_blocks", "ops_per_repeat", "repeats",
    "median_total_us", "median_us_per_xor", "median_gib_per_s", "checksum",
)
MULADD_CALIBRATION_FIELDS = (
    "block_bytes", "pool_blocks", "ops_per_repeat", "repeats",
    "xor_median_us_per_op", "xor_median_gib_per_s",
    "muladd_median_us_per_op", "muladd_median_gib_per_s",
    "muladd_xor_ratio", "checksum",
)


class RankError(ValueError):
    """A concise user-facing validation error."""


def location(source, line=None):
    return "%s%s" % (source, ":%d" % line if line is not None else "")


def fail(source, line, message):
    raise RankError("%s: %s" % (location(source, line), message))


def parse_uint(text, source, line, name, minimum=0, maximum=MAX_UINT32):
    if not isinstance(text, str) or not re.fullmatch(r"[0-9]+", text):
        fail(source, line, "%s is not an unsigned integer: %r" % (name, text))
    try:
        value = int(text, 10)
    except (ValueError, OverflowError):
        fail(source, line, "%s is outside the integer parser range" % name)
    if value < minimum or value > maximum:
        fail(source, line, "%s=%d outside [%d,%d]" %
             (name, value, minimum, maximum))
    return value


def parse_finite(text, source, line, name, minimum=0.0, maximum=None,
                 strictly_positive=False):
    try:
        value = float(text)
    except (TypeError, ValueError, OverflowError):
        fail(source, line, "%s is not numeric: %r" % (name, text))
    if not math.isfinite(value):
        fail(source, line, "%s must be finite: %r" % (name, text))
    if strictly_positive and value <= 0.0:
        fail(source, line, "%s must be positive: %r" % (name, text))
    if value < minimum or (maximum is not None and value > maximum):
        upper = "inf" if maximum is None else "%.9g" % maximum
        fail(source, line, "%s=%.9g outside [%.9g,%s]" %
             (name, value, minimum, upper))
    return value


def parse_bb_list(text):
    if not text:
        raise RankError("--bb must be a nonempty comma-separated list")
    values = []
    for token in text.split(","):
        values.append(parse_uint(token, "--bb", None, "block_bytes", 1))
    if len(values) != len(set(values)):
        raise RankError("--bb contains duplicate block sizes")
    return values


def parse_cli_float(text, name, minimum, maximum=None, positive=False):
    return parse_finite(text, name, None, name, minimum, maximum, positive)


def parse_def_pdf(text, source="def_pdf", line=None, trials=None):
    """Parse a normalized empirical deficit PDF and recover exact counts."""
    if not text:
        fail(source, line, "empty def_pdf")
    probabilities = {}
    for item in text.split("|"):
        parts = item.split(":")
        if len(parts) != 2 or not parts[0] or not parts[1]:
            fail(source, line, "malformed def_pdf item %r" % item)
        deficit = parse_uint(parts[0], source, line, "def_pdf deficit", 0, 999999)
        if deficit > 65535 and deficit != 999999:
            fail(source, line, "def_pdf deficit %d outside model domain" % deficit)
        probability = parse_finite(
            parts[1], source, line, "def_pdf probability", 0.0, 1.0)
        if deficit in probabilities:
            fail(source, line, "duplicate def_pdf deficit %d" % deficit)
        probabilities[deficit] = probability

    tolerance = 2e-5 + len(probabilities) * 1e-6
    mass = math.fsum(probabilities.values())
    if abs(mass - 1.0) > tolerance:
        fail(source, line, "def_pdf mass %.9g is not 1 (tolerance %.3g)" %
             (mass, tolerance))

    counts = None
    if trials is not None:
        counts = {}
        for deficit, probability in probabilities.items():
            count = int(round(probability * trials))
            if abs(probability - count / float(trials)) > DISPLAY_TOLERANCE:
                fail(source, line,
                     "def_pdf probability %.9g is not a six-decimal empirical "
                     "count for trials=%d" % (probability, trials))
            counts[deficit] = count
        if sum(counts.values()) != trials:
            fail(source, line, "def_pdf rounded counts sum to %d, expected %d" %
                 (sum(counts.values()), trials))
    return probabilities, counts, tolerance


def fail_count_at(counts, heavy_rows):
    return sum(count for deficit, count in counts.items() if deficit > heavy_rows)


def fail_at(pdf, heavy_rows):
    return math.fsum(probability for deficit, probability in pdf.items()
                     if deficit > heavy_rows)


def _read_csv(path, allow_preamble=False):
    source = str(path)
    try:
        with open(path, "r", encoding="utf-8", newline="") as handle:
            lines = handle.readlines()
    except (OSError, UnicodeError) as exc:
        raise RankError("%s: cannot read: %s" % (source, exc)) from exc
    if not lines:
        fail(source, None, "empty CSV")

    start = 0
    if allow_preamble:
        while start < len(lines) and lines[start].startswith("#"):
            start += 1
    if start >= len(lines):
        fail(source, None, "CSV contains no header")
    if not lines[start].strip():
        fail(source, start + 1, "blank line where CSV header was expected")

    reader = csv.reader(io.StringIO("".join(lines[start:])), strict=True)
    try:
        header = tuple(next(reader))
        rows = []
        for values in reader:
            physical_line = start + reader.line_num
            if not values:
                fail(source, physical_line, "blank CSV row")
            rows.append((physical_line, values))
    except (csv.Error, StopIteration) as exc:
        line = start + max(1, reader.line_num)
        fail(source, line, "malformed CSV: %s" % exc)
    return source, start + 1, header, rows


def load_calibration(path, label):
    """Strictly load one documented calibration schema."""
    if not path:
        return {}
    source, header_line, header, rows = _read_csv(path, allow_preamble=True)
    if label == "xor":
        supported = {XOR_CALIBRATION_FIELDS, MULADD_CALIBRATION_FIELDS}
    else:
        supported = {MULADD_CALIBRATION_FIELDS}
    if header not in supported:
        fail(source, header_line,
             "unsupported %s calibration header" % label)
    if not rows:
        fail(source, header_line, "calibration CSV has no data rows")

    table = {}
    for line, values in rows:
        if len(values) != len(header):
            fail(source, line, "wrong field count: got %d, expected %d" %
                 (len(values), len(header)))
        raw = dict(zip(header, values))
        block_bytes = parse_uint(raw["block_bytes"], source, line,
                                 "block_bytes", 1)
        if block_bytes in table:
            fail(source, line, "duplicate block_bytes %d" % block_bytes)
        integer_names = {"working_blocks", "pool_blocks", "ops_per_repeat",
                         "repeats"}
        for name in integer_names & set(header):
            parse_uint(raw[name], source, line, name, 1)
        for name in set(header) - integer_names - {"block_bytes", "checksum"}:
            parse_finite(raw[name], source, line, name, 0.0,
                         strictly_positive=True)
        if not re.fullmatch(r"0x[0-9A-Fa-f]+", raw["checksum"]):
            fail(source, line, "invalid checksum %r" % raw["checksum"])

        if header == MULADD_CALIBRATION_FIELDS:
            internal = (float(raw["muladd_median_us_per_op"])
                        / float(raw["xor_median_us_per_op"]))
            recorded = float(raw["muladd_xor_ratio"])
            if abs(internal - recorded) > 0.00055:
                fail(source, line,
                     "muladd_xor_ratio %.9g inconsistent with timings %.9g" %
                     (recorded, internal))
            column = ("xor_median_us_per_op" if label == "xor"
                      else "muladd_median_us_per_op")
        else:
            column = "median_us_per_xor"
        table[block_bytes] = float(raw[column])
    return table


def build_ratios(bb_list, xor_us, muladd_us, assume_ratio):
    ratios = {}
    for block_bytes in bb_list:
        if block_bytes in xor_us and block_bytes in muladd_us:
            ratio = muladd_us[block_bytes] / xor_us[block_bytes]
            if not math.isfinite(ratio) or ratio <= 0.0:
                raise RankError("calibration ratio for block_bytes=%d is invalid" %
                                block_bytes)
            ratios[block_bytes] = {"value": ratio, "source": "calibrated"}
        else:
            ratios[block_bytes] = {"value": assume_ratio,
                                   "source": "assumed"}
    return ratios


def _validate_result_domains(values, source, line):
    if values["K"] < 2 or values["K"] > 64000:
        fail(source, line, "K outside [2,64000]")
    if values["D"] > 4096 or values["H"] > 128:
        fail(source, line, "D/H outside supported model domain")
    if values["K"] + values["D"] + values["H"] > 65535:
        fail(source, line, "K+D+H exceeds 65535")
    if values["K"] + values["oh"] > 65535:
        fail(source, line, "K+oh exceeds 65535")
    if values["trials"] > MAX_TRIALS:
        fail(source, line,
             "trials exceeds %d; six-decimal PDFs cannot recover exact counts" %
             MAX_TRIALS)


def score_result_row(raw, header, source, line):
    if (len(raw["scheme"]) > 128
            or not re.fullmatch(r"[A-Za-z0-9_]+", raw["scheme"])):
        fail(source, line, "invalid scheme %r" % raw["scheme"])
    values = {}
    for name in header:
        if name in {"scheme", "def_pdf"}:
            continue
        if name in INTEGER_FIELDS:
            minimum = 1 if name == "trials" else 0
            maximum = MAX_TRIALS if name == "trials" else MAX_UINT32
            values[name] = parse_uint(raw[name], source, line, name,
                                      minimum, maximum)
        elif name in RATE_FIELDS:
            values[name] = parse_finite(raw[name], source, line, name, 0.0, 1.0)
        else:
            values[name] = parse_finite(raw[name], source, line, name, 0.0)
    _validate_result_domains(values, source, line)

    pdf, counts, _tolerance = parse_def_pdf(
        raw["def_pdf"], source, line, values["trials"])
    completed = [(deficit, count) for deficit, count in counts.items()
                 if deficit != 999999]
    observed_max = max((deficit for deficit, count in completed if count > 0),
                       default=0)
    if values["def_max"] != observed_max:
        fail(source, line, "def_max %d != def_pdf maximum %d" %
             (values["def_max"], observed_max))

    completed_trials = sum(count for _deficit, count in completed)
    expected_def_mu = (sum(deficit * count for deficit, count in completed)
                       / float(completed_trials) if completed_trials else 0.0)
    if abs(values["def_mu"] - expected_def_mu) > 0.0000501:
        fail(source, line, "def_mu %.9g != def_pdf empirical mean %.9g" %
             (values["def_mu"], expected_def_mu))

    runaway_count = counts.get(999999, 0)
    runaway_empirical = runaway_count / float(values["trials"])
    runaway_rate = values.get("runaway_rate", 0.0)
    if "runaway_rate" not in values and runaway_count != 0:
        fail(source, line, "def_pdf runaway sentinel requires runaway_rate")
    if abs(runaway_rate - runaway_empirical) > DISPLAY_TOLERANCE:
        fail(source, line, "runaway_rate %.9g != def_pdf empirical rate %.9g" %
             (runaway_rate, runaway_empirical))

    failures = fail_count_at(counts, values["H"])
    noheavy_failures = fail_count_at(counts, 0)
    empirical_fail = failures / float(values["trials"])
    empirical_noheavy = noheavy_failures / float(values["trials"])
    if abs(values["fail_rate"] - empirical_fail) > DISPLAY_TOLERANCE:
        fail(source, line, "fail_rate %.9g != def_pdf empirical rate %.9g" %
             (values["fail_rate"], empirical_fail))
    if abs(values["fail_rate_noheavy"] - empirical_noheavy) > DISPLAY_TOLERANCE:
        fail(source, line,
             "fail_rate_noheavy %.9g != def_pdf empirical rate %.9g" %
             (values["fail_rate_noheavy"], empirical_noheavy))

    if values["inact_mu"] > values["inact_max"] + 0.011:
        fail(source, line, "inact_mu exceeds inact_max")
    if values["rank_mu"] > values["inact_mu"] + 0.011:
        fail(source, line, "rank_mu exceeds inact_mu")
    if values["rank_mu"] > values["residual_rows_mu"] + 0.011:
        fail(source, line, "rank_mu exceeds residual_rows_mu")
    if ("def_band_w95" in values
            and values["def_band_w95"] > values["def_band_w99"]):
        fail(source, line, "def_band_w95 exceeds def_band_w99")
    expected_heavy_divs = values["H"] if completed_trials else 0
    if ("heavy_divs_mu" in values
            and abs(values["heavy_divs_mu"] - expected_heavy_divs) > 0.050001):
        fail(source, line, "heavy_divs_mu does not match completed-trial H")

    recv = values["recv_xors_per_packet"] * (values["K"] + values["oh"])
    if "ge_real_rowops_mu" in values:
        ge = values["ge_real_rowops_mu"]
        ge_source = "real"
    else:
        ge = values["ge_block_xors_mu"]
        ge_source = "estimated"
    base = math.fsum((recv, values["precode_gen_xors_mu"],
                      values["sparse_solve_xors_mu"],
                      values["backsub_xors_mu"], ge))
    if not math.isfinite(base):
        fail(source, line, "base cost overflows")
    if "heavy_muladds_mu" in values:
        heavy = values["heavy_muladds_mu"]
        heavy_source = "recorded"
    else:
        heavy = values["H"] * values["inact_mu"] + values["H"] ** 2
        heavy_source = "derived"
    if not math.isfinite(heavy):
        fail(source, line, "heavy cost overflows")

    return {
        "row_id": "%s:%d" % (source, line),
        "source": source,
        "line": line,
        "K": values["K"],
        "D": values["D"],
        "modeled_h": values["H"],
        "oh": values["oh"],
        "trials": values["trials"],
        "scheme": raw["scheme"],
        "pivot_window": values.get("pivot_window"),
        "base_ops": base,
        "heavy_muladds": heavy,
        "heavy_source": heavy_source,
        "inact_mu": values["inact_mu"],
        "ge_source": ge_source,
        "runaway_rate": runaway_empirical,
        "pdf": pdf,
        "def_counts": counts,
    }


def load_rows(paths):
    if not paths:
        raise RankError("no precode result CSVs")
    normalized = [os.path.normpath(path) for path in paths]
    identities = [os.path.realpath(path) for path in normalized]
    if len(identities) != len(set(identities)):
        raise RankError("duplicate precode result path")
    rows = []
    for path in sorted(normalized):
        source, header_line, header, csv_rows = _read_csv(path)
        if header not in RESULT_SCHEMAS:
            fail(source, header_line, "unsupported precode result header")
        if not csv_rows:
            fail(source, header_line, "precode result CSV has no data rows")
        seen = set()
        for line, fields in csv_rows:
            if len(fields) != len(header):
                fail(source, line, "wrong field count: got %d, expected %d" %
                     (len(fields), len(header)))
            raw = dict(zip(header, fields))
            row = score_result_row(raw, header, source, line)
            key = (row["pivot_window"], row["K"], row["scheme"], row["D"],
                   row["modeled_h"], row["oh"])
            if key in seen:
                fail(source, line, "duplicate experiment row %r" % (key,))
            seen.add(key)
            rows.append(row)
    return rows


def wilson_upper(failures, trials, confidence):
    """One-sided Wilson score upper confidence bound for a binomial rate."""
    if trials <= 0 or failures < 0 or failures > trials:
        raise ValueError("invalid binomial observation")
    if not (0.5 < confidence < 1.0):
        raise ValueError("confidence must be between 0.5 and 1")
    z = statistics.NormalDist().inv_cdf(confidence)
    p = failures / float(trials)
    z2 = z * z
    denominator = 1.0 + z2 / trials
    center = p + z2 / (2.0 * trials)
    radius = z * math.sqrt((p * (1.0 - p) + z2 / (4.0 * trials)) / trials)
    return min(1.0, (center + radius) / denominator)


def pareto_ids(rows, metric="upper_failure"):
    """Return IDs on the minimize(total_ops, failure metric) Pareto front."""
    result = []
    for candidate in rows:
        dominated = False
        for other in rows:
            if other is candidate:
                continue
            no_worse = (other["total_ops"] <= candidate["total_ops"]
                        and other[metric] <= candidate[metric])
            strictly_better = (other["total_ops"] < candidate["total_ops"]
                               or other[metric] < candidate[metric])
            if no_worse and strictly_better:
                dominated = True
                break
        if not dominated:
            result.append(candidate["row_id"])
    return sorted(result)


def crossover_ratio(dense_row, alternate_row, pessimistic=False):
    factor = 2.0 if pessimistic else 1.0
    heavy_delta = factor * (alternate_row["heavy_muladds"]
                            - dense_row["heavy_muladds"])
    base_advantage = dense_row["base_ops"] - alternate_row["base_ops"]
    if heavy_delta > 0.0:
        ratio = base_advantage / heavy_delta
        if ratio <= 0.0:
            return None, "never beats (base cost already higher)"
        return ratio, "beats below this ratio"
    if base_advantage > 0.0:
        return None, "beats at every ratio (cheaper base and heavy)"
    if heavy_delta < 0.0:
        return base_advantage / heavy_delta, "beats above this ratio"
    return None, "never beats (heavier base, no heavy advantage)"


def _selected_h(row, reliability_h):
    return row["modeled_h"] if reliability_h == "modeled" else reliability_h


def evaluate_row(row, block_bytes, ratio, options):
    reliability_h = _selected_h(row, options["reliability_h"])
    failures = fail_count_at(row["def_counts"], reliability_h)
    observed = failures / float(row["trials"])
    upper = wilson_upper(failures, row["trials"], options["confidence"])

    if options["cost_h"] == "modeled":
        cost_h = row["modeled_h"]
        heavy = row["heavy_muladds"]
        heavy_source = row["heavy_source"]
        cost_semantics = "modeled_h_ledger"
    else:
        cost_h = reliability_h
        heavy = cost_h * row["inact_mu"] + cost_h ** 2
        heavy_source = "recomputed"
        cost_semantics = "heavy_recomputed_using_modeled_inact_base_held_at_modeled_h"
    factor = 2.0 if options["pessimistic"] else 1.0
    total = row["base_ops"] + ratio * factor * heavy
    if not math.isfinite(total):
        fail(row["source"], row["line"], "total cost overflows")

    reasons = []
    if row["runaway_rate"] > 0.0:
        reasons.append("runaway_rate>0")
    if row["trials"] < options["min_trials"]:
        reasons.append("trials<%d" % options["min_trials"])
    if observed > options["max_observed_fail"]:
        reasons.append("observed_failure>%.9g" % options["max_observed_fail"])
    if upper > options["max_upper_fail"]:
        reasons.append("wilson_upper>%.9g" % options["max_upper_fail"])

    failure_by_h = {}
    for heavy_rows in sorted(set(H_PRIMES + (reliability_h,))):
        count = fail_count_at(row["def_counts"], heavy_rows)
        failure_by_h[str(heavy_rows)] = {
            "failures": count,
            "rate": count / float(row["trials"]),
        }
    return {
        "row_id": row["row_id"],
        "source": row["source"],
        "line": row["line"],
        "scheme": row["scheme"],
        "D": row["D"],
        "pivot_window": row["pivot_window"],
        "modeled_h": row["modeled_h"],
        "reliability_h": reliability_h,
        "cost_h": cost_h,
        "cost_semantics": cost_semantics,
        "trials": row["trials"],
        "failures": failures,
        "observed_failure": observed,
        "wilson_upper": upper,
        "failure_by_h": failure_by_h,
        "base_ops": row["base_ops"],
        "heavy_muladds": heavy,
        "total_ops": None if row["runaway_rate"] > 0.0 else total,
        "ge_source": row["ge_source"],
        "heavy_source": heavy_source,
        "runaway_rate": row["runaway_rate"],
        "eligible": not reasons,
        "rejection_reasons": reasons,
        "eligible_rank": None,
        "pareto": False,
        "winner": False,
        "block_bytes": block_bytes,
    }


def build_report(rows, bb_list, ratios, options, inputs):
    groups = {}
    modeled_groups = {}
    for row in rows:
        modeled_groups.setdefault((row["K"], row["oh"]), []).append(row)
        for block_bytes in bb_list:
            key = (row["K"], row["oh"], block_bytes)
            evaluated = evaluate_row(
                row, block_bytes, ratios[block_bytes]["value"], options)
            groups.setdefault(key, []).append(evaluated)

    output_groups = []
    pareto_metric = options["pareto_metric"]
    for key in sorted(groups):
        candidates = groups[key]
        eligible = [row for row in candidates if row["eligible"]]
        eligible.sort(key=lambda row: (
            row["total_ops"], row[pareto_metric], row["scheme"],
            row["source"], row["line"]))
        prior_cost = None
        prior_rank = None
        for index, row in enumerate(eligible):
            if prior_cost is None or row["total_ops"] != prior_cost:
                prior_rank = index + 1
                prior_cost = row["total_ops"]
            row["eligible_rank"] = prior_rank
        front = pareto_ids(eligible, pareto_metric)
        for row in eligible:
            row["pareto"] = row["row_id"] in front
        winners = []
        if eligible:
            minimum = eligible[0]["total_ops"]
            minimum_failure = eligible[0][pareto_metric]
            winners = sorted(row["row_id"] for row in eligible
                             if row["total_ops"] == minimum
                             and row[pareto_metric] == minimum_failure)
            for row in eligible:
                row["winner"] = row["row_id"] in winners

        rejected = sorted(
            (row for row in candidates if not row["eligible"]),
            key=lambda row: (row["scheme"], row["source"], row["line"]))
        ordered = eligible + rejected
        k, oh, block_bytes = key
        output_groups.append({
            "K": k,
            "oh": oh,
            "block_bytes": block_bytes,
            "muladd_xor_ratio": ratios[block_bytes]["value"],
            "ratio_source": ratios[block_bytes]["source"],
            "winner_ids": winners,
            "pareto_ids": front,
            "rows": ordered,
        })

    crossovers = []
    for (k, oh), modeled in sorted(modeled_groups.items()):
        dense = [
            row for row in modeled
            if row["runaway_rate"] == 0.0 and row["modeled_h"] == 6
            and (row["scheme"] == "dense"
                 or row["scheme"].startswith("dense_d")
                 or row["scheme"] == "dense_h6")
        ]
        alternate = [
            row for row in modeled
            if row["runaway_rate"] == 0.0 and row["modeled_h"] == 12
            and row["scheme"].startswith("ldpcdense")
        ]
        for dense_row in sorted(dense, key=lambda row: row["row_id"]):
            for alternate_row in sorted(alternate,
                                        key=lambda row: row["row_id"]):
                ratio, verdict = crossover_ratio(
                    dense_row, alternate_row, options["pessimistic"])
                crossovers.append({
                    "K": k,
                    "oh": oh,
                    "dense_id": dense_row["row_id"],
                    "alternate_id": alternate_row["row_id"],
                    "muladd_xor_ratio": ratio,
                    "verdict": verdict,
                    "cost_semantics": "modeled_h_ledger",
                })

    reliability_h = options["reliability_h"]
    return {
        "schema": "wirehair.precode-ranking.v1",
        "inputs": sorted(os.path.normpath(path) for path in inputs),
        "assumptions": {
            "confidence": options["confidence"],
            "confidence_method": "one-sided Wilson score upper bound",
            "cost_h": options["cost_h"],
            "cost_semantics": (
                "complete ledger held at each row's modeled H"
                if options["cost_h"] == "modeled" else
                "heavy H*inact_mu+H^2 recomputed at reliability H using "
                "modeled-H inact_mu; base ledger held at modeled H"),
            "max_observed_fail": options["max_observed_fail"],
            "max_upper_fail": options["max_upper_fail"],
            "min_trials": options["min_trials"],
            "pareto_metric": pareto_metric,
            "pessimistic_heavy_factor": 2 if options["pessimistic"] else 1,
            "reliability_h": reliability_h,
        },
        "crossovers": crossovers,
        "groups": output_groups,
    }


def render_human(report):
    assumptions = report["assumptions"]
    out = io.StringIO()
    out.write("rank_total report v1\n")
    out.write("confidence: %.6g one-sided Wilson score upper bound\n" %
              assumptions["confidence"])
    out.write("reliability: H=%s min_trials=%d observed<=%.9g upper<=%.9g\n" %
              (assumptions["reliability_h"], assumptions["min_trials"],
               assumptions["max_observed_fail"],
               assumptions["max_upper_fail"]))
    out.write("cost: %s; heavy_factor=%d\n" %
              (assumptions["cost_semantics"],
               assumptions["pessimistic_heavy_factor"]))
    out.write("pareto: minimize total_ops and %s among eligible rows\n\n" %
              assumptions["pareto_metric"])

    for group in report["groups"]:
        out.write("== K=%d oh=%d bb=%d (muladd:xor=%.6g, %s) ==\n" %
                  (group["K"], group["oh"], group["block_bytes"],
                   group["muladd_xor_ratio"], group["ratio_source"]))
        out.write("rank status  P W scheme                 source:line"
                  "                         mH rH cH trials fail upper total_ops\n")
        for row in group["rows"]:
            rank = "-" if row["eligible_rank"] is None else str(row["eligible_rank"])
            status = "ok" if row["eligible"] else "reject"
            source = "%s:%d" % (os.path.basename(row["source"]), row["line"])
            total = ("unavailable" if row["total_ops"] is None
                     else "%.1f" % row["total_ops"])
            out.write("%4s %-6s  %s %s %-22s %-35s %2d %2d %2d %6d "
                      "%.6f %.6f %s" %
                      (rank, status, "Y" if row["pareto"] else "-",
                       "Y" if row["winner"] else "-", row["scheme"], source,
                       row["modeled_h"], row["reliability_h"], row["cost_h"],
                       row["trials"], row["observed_failure"],
                       row["wilson_upper"], total))
            if row["rejection_reasons"]:
                out.write(" [%s]" % ",".join(row["rejection_reasons"]))
            out.write("\n")
        if group["winner_ids"]:
            out.write("winner: %s\n" % ", ".join(group["winner_ids"]))
        else:
            out.write("winner: no feasible scheme\n")
        out.write("pareto: %s\n\n" %
                  (", ".join(group["pareto_ids"]) if group["pareto_ids"]
                   else "empty"))
    if report["crossovers"]:
        out.write("modeled-H ldpcdense/dense crossovers:\n")
        for item in report["crossovers"]:
            ratio = ("n/a" if item["muladd_xor_ratio"] is None
                     else "%.6g" % item["muladd_xor_ratio"])
            out.write("K=%d oh=%d alternate=%s dense=%s ratio=%s (%s)\n" %
                      (item["K"], item["oh"], item["alternate_id"],
                       item["dense_id"], ratio, item["verdict"]))
    return out.getvalue()


def deterministic_json(report):
    return json.dumps(report, sort_keys=True, indent=2, allow_nan=False) + "\n"


def atomic_write(path, content):
    directory = os.path.dirname(os.path.abspath(path))
    basename = os.path.basename(path)
    fd, temporary = tempfile.mkstemp(prefix=".%s." % basename, dir=directory,
                                     text=True)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except OSError:
            pass
        raise


def self_test():
    pdf, counts, _ = parse_def_pdf("0:0.900000|7:0.100000", trials=100)
    assert fail_at(pdf, 6) == 0.1
    assert fail_count_at(counts, 6) == 10
    assert abs(wilson_upper(0, 100, 0.95) - 0.026342720) < 1e-8
    candidates = [
        {"row_id": "a", "total_ops": 1.0, "upper_failure": 0.2},
        {"row_id": "b", "total_ops": 2.0, "upper_failure": 0.1},
        {"row_id": "c", "total_ops": 3.0, "upper_failure": 0.3},
    ]
    assert pareto_ids(candidates) == ["a", "b"]
    for bad in ("", "0:nan|1:1", "-1:1", "0:0.5|0:0.5", "0:0.2"):
        try:
            parse_def_pdf(bad, trials=10)
        except RankError:
            pass
        else:
            raise AssertionError("accepted malformed def_pdf %r" % bad)
    suite = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "test_rank_total.py")
    completed = subprocess.run([sys.executable, suite, "-q"], check=False)
    if completed.returncode != 0:
        return completed.returncode
    print("rank_total self-test passed")
    return 0


def parse_options(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("csvs", nargs="*", help="precode_sim result CSVs")
    parser.add_argument("--xor-csv", help="XOR calibration CSV")
    parser.add_argument("--muladd-csv", help="muladd calibration CSV")
    parser.add_argument("--bb", default=DEFAULT_BB,
                        help="comma-separated block sizes (default %(default)s)")
    parser.add_argument("--assume-ratio", default="4.0",
                        help="positive fallback muladd:xor ratio")
    parser.add_argument("--reliability-h", default="modeled",
                        help="modeled or an alternative H in [0,128]")
    parser.add_argument("--cost-h", choices=("modeled", "reliability"),
                        default="modeled",
                        help="modeled ledger or recompute its heavy term at reliability H")
    parser.add_argument("--min-trials", default="100",
                        help="minimum trials for eligibility (default %(default)s)")
    parser.add_argument("--max-observed-fail", default="0.005",
                        help="observed failure gate (default %(default)s)")
    parser.add_argument("--max-upper-fail", default="0.05",
                        help="Wilson upper-bound gate (default %(default)s)")
    parser.add_argument("--confidence", default="0.95",
                        help="one-sided Wilson confidence in (0.5,1)")
    parser.add_argument("--pareto-metric", choices=("observed_failure",
                                                     "wilson_upper"),
                        default="wilson_upper")
    parser.add_argument("--pessimistic", action="store_true",
                        help="double the GF(256) heavy term")
    parser.add_argument("--json-out", help="write deterministic JSON report")
    parser.add_argument("--self-test", action="store_true")
    return parser, parser.parse_args(argv)


def main(argv=None):
    parser, args = parse_options(argv)
    if args.self_test:
        return self_test()
    if not args.csvs:
        parser.error("no precode result CSVs (or use --self-test)")

    try:
        bb_list = parse_bb_list(args.bb)
        assume_ratio = parse_cli_float(args.assume_ratio, "--assume-ratio", 0.0,
                                       positive=True)
        confidence = parse_cli_float(args.confidence, "--confidence", 0.0, 1.0)
        if not 0.5 < confidence < 1.0:
            raise RankError("--confidence must be strictly between 0.5 and 1")
        min_trials = parse_uint(args.min_trials, "--min-trials", None,
                                "--min-trials", 1, MAX_TRIALS)
        max_observed = parse_cli_float(args.max_observed_fail,
                                       "--max-observed-fail", 0.0, 1.0)
        max_upper = parse_cli_float(args.max_upper_fail,
                                    "--max-upper-fail", 0.0, 1.0)
        if args.reliability_h == "modeled":
            reliability_h = "modeled"
        else:
            reliability_h = parse_uint(args.reliability_h, "--reliability-h",
                                       None, "--reliability-h", 0, 128)

        output_real = os.path.realpath(args.json_out) if args.json_out else None
        protected = args.csvs + [p for p in (args.xor_csv, args.muladd_csv) if p]
        if output_real and any(output_real == os.path.realpath(path)
                               for path in protected):
            raise RankError("--json-out must not overwrite an input")

        xor_us = load_calibration(args.xor_csv, "xor")
        muladd_us = load_calibration(args.muladd_csv, "muladd")
        ratios = build_ratios(bb_list, xor_us, muladd_us, assume_ratio)
        rows = load_rows(args.csvs)
        options = {
            "confidence": confidence,
            "cost_h": args.cost_h,
            "max_observed_fail": max_observed,
            "max_upper_fail": max_upper,
            "min_trials": min_trials,
            "pareto_metric": args.pareto_metric,
            "pessimistic": args.pessimistic,
            "reliability_h": reliability_h,
        }
        report = build_report(rows, bb_list, ratios, options, args.csvs)
        human = render_human(report)
        machine = deterministic_json(report)
        if args.json_out:
            atomic_write(args.json_out, machine)
    except (RankError, OSError, UnicodeError, csv.Error) as exc:
        sys.stderr.write("rank_total: %s\n" % exc)
        return 1

    sys.stdout.write(human)
    return 0


if __name__ == "__main__":
    sys.exit(main())
