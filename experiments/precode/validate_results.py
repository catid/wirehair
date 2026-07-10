#!/usr/bin/env python3
"""Strict validator for precode_sim result CSV files."""

import argparse
import csv
import itertools
import json
import math
import re
import sys
from pathlib import Path


BASE = [
    "K", "scheme", "D", "H", "oh", "trials", "fail_rate",
    "fail_rate_noheavy", "def_mu", "def_max", "def_pdf", "inact_mu",
    "inact_sd", "inact_max", "rank_mu", "peeled_mu", "residual_rows_mu",
    "recv_xors_per_packet", "precode_gen_xors_mu", "sparse_solve_xors_mu",
    "backsub_xors_mu", "ge_block_xors_mu", "ge_bitops_mu",
]
HEAVY = ["heavy_muladds_mu", "heavy_divs_mu"]
REPLAY = [
    "ge_real_bitops_mu", "ge_real_rowops_mu", "fill_in_mu",
    "def_outside_w18_rate", "def_band_w95", "def_band_w99",
]
SCHEMAS = {
    tuple(BASE),
    tuple(BASE + HEAVY + REPLAY),
    tuple(BASE + HEAVY + REPLAY + ["runaway_rate"]),
    tuple(BASE + HEAVY + ["runaway_rate"]),
    tuple(["pivot_window"] + BASE + HEAVY + REPLAY + ["runaway_rate"]),
}
INTEGER_FIELDS = {
    "pivot_window", "K", "D", "H", "oh", "trials", "def_max",
    "inact_max", "def_band_w95", "def_band_w99",
}
RATE_FIELDS = {
    "fail_rate", "fail_rate_noheavy", "def_outside_w18_rate", "runaway_rate",
}


class ValidationError(ValueError):
    pass


def uint(text, source, line, name, minimum=0, maximum=None):
    if not re.fullmatch(r"[0-9]+", text):
        raise ValidationError(f"{source}:{line}: {name} is not an unsigned integer: {text!r}")
    value = int(text, 10)
    if value < minimum or (maximum is not None and value > maximum):
        raise ValidationError(f"{source}:{line}: {name}={value} outside supported range")
    return value


def finite(text, source, line, name, minimum=0.0, maximum=None):
    try:
        value = float(text)
    except ValueError as exc:
        raise ValidationError(f"{source}:{line}: {name} is not numeric: {text!r}") from exc
    if not math.isfinite(value) or value < minimum or (maximum is not None and value > maximum):
        raise ValidationError(f"{source}:{line}: {name}={text!r} outside supported finite range")
    return value


def deficit_pdf(text, source, line):
    if not text:
        raise ValidationError(f"{source}:{line}: empty def_pdf")
    result = {}
    for item in text.split("|"):
        fields = item.split(":")
        if len(fields) != 2:
            raise ValidationError(f"{source}:{line}: malformed def_pdf item {item!r}")
        deficit = uint(fields[0], source, line, "def_pdf deficit")
        probability = finite(fields[1], source, line, "def_pdf probability", 0.0, 1.0)
        if deficit in result:
            raise ValidationError(f"{source}:{line}: duplicate def_pdf deficit {deficit}")
        result[deficit] = probability
    total = sum(result.values())
    tolerance = 2e-5 + len(result) * 1e-6
    if abs(total - 1.0) > tolerance:
        raise ValidationError(f"{source}:{line}: def_pdf mass {total:.9g} is not 1")
    return result, tolerance


def validate_file(path, seen=None):
    source = str(path)
    rows = []
    seen = {} if seen is None else seen
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, strict=True)
        if reader.fieldnames is None:
            raise ValidationError(f"{source}: empty CSV")
        if tuple(reader.fieldnames) not in SCHEMAS:
            raise ValidationError(f"{source}: unsupported precode CSV header")
        for line, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise ValidationError(f"{source}:{line}: wrong field count")
            if not re.fullmatch(r"[A-Za-z0-9_]+", row["scheme"]):
                raise ValidationError(f"{source}:{line}: invalid scheme {row['scheme']!r}")
            values = {}
            for name, text in row.items():
                if name == "scheme" or name == "def_pdf":
                    continue
                if name in INTEGER_FIELDS:
                    minimum = 1 if name == "trials" else 0
                    maximum = 64000 if name == "K" else None
                    values[name] = uint(text, source, line, name, minimum, maximum)
                elif name in RATE_FIELDS:
                    values[name] = finite(text, source, line, name, 0.0, 1.0)
                else:
                    values[name] = finite(text, source, line, name)
            if (values["K"] < 2 or values["D"] > 4096 or values["H"] > 128 or
                    values["K"] + values["D"] + values["H"] > 65535 or
                    values["K"] + values["oh"] > 65535):
                raise ValidationError(f"{source}:{line}: K/D/H/oh outside model domain")
            pdf, tolerance = deficit_pdf(row["def_pdf"], source, line)
            if any(deficit > 65535 and deficit != 999999 for deficit in pdf):
                raise ValidationError(f"{source}:{line}: def_pdf deficit outside model domain")
            completed_deficits = [deficit for deficit in pdf if deficit != 999999]
            observed_max = max(completed_deficits, default=0)
            if values["def_max"] != observed_max:
                raise ValidationError(
                    f"{source}:{line}: def_max {values['def_max']} != def_pdf maximum "
                    f"{observed_max}"
                )
            runaway_probability = pdf.get(999999, 0.0)
            if "runaway_rate" not in values and runaway_probability > 0.0:
                raise ValidationError(
                    f"{source}:{line}: sentinel def_pdf bin requires runaway_rate"
                )
            if ("runaway_rate" in values and
                    abs(values["runaway_rate"] - runaway_probability) > tolerance):
                raise ValidationError(
                    f"{source}:{line}: runaway_rate {values['runaway_rate']} != "
                    f"def_pdf sentinel probability {runaway_probability}"
                )
            fail = sum(probability for deficit, probability in pdf.items() if deficit > values["H"])
            noheavy = sum(probability for deficit, probability in pdf.items() if deficit > 0)
            if abs(fail - values["fail_rate"]) > tolerance:
                raise ValidationError(
                    f"{source}:{line}: fail_rate {values['fail_rate']} != def_pdf rescore {fail}"
                )
            if abs(noheavy - values["fail_rate_noheavy"]) > tolerance:
                raise ValidationError(
                    f"{source}:{line}: fail_rate_noheavy {values['fail_rate_noheavy']} "
                    f"!= def_pdf rescore {noheavy}"
                )
            key = (
                values.get("pivot_window"), values["K"], row["scheme"],
                values["D"], values["H"], values["oh"],
            )
            previous = seen.get(key)
            if previous is not None:
                previous_source, previous_line = previous
                raise ValidationError(
                    f"{source}:{line}: duplicate experiment row {key}; "
                    f"first seen at {previous_source}:{previous_line}"
                )
            seen[key] = (source, line)
            rows.append({**values, "scheme": row["scheme"]})
    if not rows:
        raise ValidationError(f"{source}: CSV has no data rows")
    return rows


def comma_values(text, converter=str):
    if text is None:
        return None
    values = [converter(item) for item in text.split(",") if item != ""]
    if not values:
        raise argparse.ArgumentTypeError("expected a nonempty comma-separated list")
    return values


def load_manifest(path):
    source = str(path)
    try:
        with path.open("r", encoding="utf-8") as handle:
            manifest = json.load(handle)
    except (OSError, json.JSONDecodeError) as exc:
        raise ValidationError(f"cannot read experiment manifest {source}: {exc}") from exc
    if not isinstance(manifest, dict):
        raise ValidationError(f"{source}: manifest must be a JSON object")
    expected_keys = {"K", "schemes", "oh", "trials"}
    if set(manifest) != expected_keys:
        raise ValidationError(
            f"{source}: manifest keys must be exactly {sorted(expected_keys)}"
        )
    for key in ("K", "schemes", "oh"):
        if not isinstance(manifest[key], list) or not manifest[key]:
            raise ValidationError(f"{source}: manifest {key} must be a nonempty list")
    if (not all(isinstance(value, int) and not isinstance(value, bool)
                for value in manifest["K"] + manifest["oh"]) or
            not all(isinstance(value, str) and value for value in manifest["schemes"]) or
            not isinstance(manifest["trials"], int) or
            isinstance(manifest["trials"], bool) or manifest["trials"] < 1):
        raise ValidationError(f"{source}: manifest contains invalid value types or domains")
    if (any(value < 2 or value > 64000 for value in manifest["K"]) or
            any(value < 0 for value in manifest["oh"])):
        raise ValidationError(f"{source}: manifest contains values outside supported ranges")
    if any(len(values) != len(set(values))
           for values in (manifest["K"], manifest["schemes"], manifest["oh"])):
        raise ValidationError(f"{source}: manifest lists must not contain duplicates")
    return manifest


def merge_expectation(name, command_value, manifest_value):
    if command_value is not None and manifest_value is not None:
        if command_value != manifest_value:
            raise ValidationError(
                f"command-line {name} expectation is incompatible with manifest"
            )
    return command_value if command_value is not None else manifest_value


def validate_grid(rows, ks, schemes, overheads, trials):
    if trials is not None:
        for row in rows:
            if row["trials"] != trials:
                raise ValidationError(
                    f"expected {trials} trials, got {row['trials']} for "
                    f"K={row['K']} scheme={row['scheme']} oh={row['oh']}"
                )
    if ks is None and schemes is None and overheads is None:
        return
    ks = ks if ks is not None else sorted({row["K"] for row in rows})
    schemes = schemes if schemes is not None else sorted({row["scheme"] for row in rows})
    overheads = overheads if overheads is not None else sorted({row["oh"] for row in rows})
    actual = {(row["K"], row["scheme"], row["oh"]) for row in rows}
    expected = set(itertools.product(ks, schemes, overheads))
    missing = sorted(expected - actual)
    extra = sorted(actual - expected)
    if missing or extra:
        raise ValidationError(f"experiment grid mismatch: missing={missing[:8]} extra={extra[:8]}")


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+")
    parser.add_argument("--expect-k")
    parser.add_argument("--expect-schemes")
    parser.add_argument("--expect-oh")
    parser.add_argument("--expect-trials", type=int)
    parser.add_argument("--manifest")
    args = parser.parse_args(argv)
    try:
        ks = comma_values(args.expect_k, int)
        schemes = comma_values(args.expect_schemes)
        overheads = comma_values(args.expect_oh, int)
        trials = args.expect_trials
        if args.manifest:
            manifest = load_manifest(Path(args.manifest))
            ks = merge_expectation("K", ks, manifest["K"])
            schemes = merge_expectation("schemes", schemes, manifest["schemes"])
            overheads = merge_expectation("oh", overheads, manifest["oh"])
            trials = merge_expectation("trials", trials, manifest["trials"])
        rows = []
        seen = {}
        for name in args.files:
            path = Path(name)
            if not path.is_file():
                raise ValidationError(f"missing precode result: {path}")
            rows.extend(validate_file(path, seen))
        validate_grid(rows, ks, schemes, overheads, trials)
    except (OSError, csv.Error, ValidationError, ValueError) as exc:
        print(exc, file=sys.stderr)
        return 1
    print(f"validated {len(args.files)} precode result file(s), {len(rows)} row(s)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
