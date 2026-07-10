#!/usr/bin/env python3
"""Strictly validate and canonically merge table-regeneration result shards."""

import argparse
import csv
import os
import sys
from pathlib import Path


MASK64 = (1 << 64) - 1

SCHEMAS = {
    "dense-count": [
        "N", "DenseCount", "SelectedFailures", "Trials", "MaxFailures",
        "LowCountRun", "BaseSeed", "TrialSeed", "QualifyingCount",
        "SourceTablesHash", "MethodVersion", "Status",
    ],
    "small": [
        "N", "DenseCount", "DenseSeed", "PeelSeed", "DenseFailures",
        "PeelScore", "CurrentDenseCount", "CurrentDenseSeed",
        "CurrentPeelSeed", "Trials", "BaseSeed", "TrialSeed",
        "SmallTablesHash", "MethodVersion",
    ],
    "dense-seeds": [
        "DenseIndex", "Used", "N", "DenseCount", "DenseSeed", "Failures",
        "CurrentDenseSeed", "Trials", "BaseSeed", "TrialSeed",
        "SourceTablesHash", "MethodVersion",
    ],
    "peel": [
        "Subdivision", "FirstN", "LastN", "PeelSeed", "Score",
        "CurrentPeelSeed", "Trials", "MaxTries", "EvaluatedCandidates",
        "BaseSeed", "RequestedNLo", "RequestedNHi", "BaseTablesHash",
        "MethodVersion", "Status",
    ],
}

RESULT_GLOBS = {
    "dense-count": "dense_count.results.tsv",
    "small": "small_dense_seeds.results.tsv",
    "dense-seeds": "most_dense_seeds.results.tsv",
    "peel": "peel_seeds.results.tsv",
}


class ValidationError(ValueError):
    pass


DOMAINS = {
    "dense-count": 0x44434F554E54,
    "small": 0x534D414C4C,
    "dense-seeds": 0x4D4F5354,
    "peel": 0x5045454C,
}


def derive_trial_seed(base_seed, key, kind):
    value = ((base_seed ^ DOMAINS[kind]) + 0x9E3779B97F4A7C15 * key) & MASK64
    value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & MASK64
    return (value ^ (value >> 31)) & MASK64


def _uint(text, source, line, name, minimum=0, maximum=None):
    if not text or not text.isascii() or not text.isdecimal():
        raise ValidationError(
            f"{source}:{line}: {name} is not an unsigned decimal: {text!r}"
        )
    value = int(text, 10)
    if value < minimum or (maximum is not None and value > maximum):
        upper = "unbounded" if maximum is None else maximum
        raise ValidationError(
            f"{source}:{line}: {name}={value} outside [{minimum}, {upper}]"
        )
    return value


def _resolve_input(kind, value):
    path = Path(value)
    if path.is_dir():
        path = path / RESULT_GLOBS[kind]
    if not path.is_file():
        raise ValidationError(f"missing {kind} results input: {path}")
    return path


def _read_rows(kind, path):
    schema = SCHEMAS[kind]
    rows = []
    try:
        with path.open("r", encoding="ascii", newline="") as handle:
            reader = csv.reader(handle, delimiter="\t", strict=True)
            try:
                header = next(reader)
            except StopIteration as exc:
                raise ValidationError(f"{path}: empty results file") from exc
            if header != schema:
                raise ValidationError(
                    f"{path}:1: expected tab-separated header {schema!r}, got {header!r}"
                )
            for line, fields in enumerate(reader, 2):
                if not fields or fields == [""]:
                    raise ValidationError(f"{path}:{line}: blank rows are not allowed")
                if len(fields) != len(schema):
                    raise ValidationError(
                        f"{path}:{line}: expected {len(schema)} fields, got {len(fields)}"
                    )
                rows.append((dict(zip(schema, fields)), str(path), line))
    except (UnicodeError, csv.Error) as exc:
        raise ValidationError(f"{path}: malformed ASCII TSV: {exc}") from exc
    if not rows:
        raise ValidationError(f"{path}: no result rows")
    return rows


def _validate_small(row, source, line):
    n = _uint(row["N"], source, line, "N", 2, 2047)
    dense_count = _uint(row["DenseCount"], source, line, "DenseCount", 1, n)
    dense_seed_max = 65535 if n < 65 else 255
    _uint(row["DenseSeed"], source, line, "DenseSeed", 0, dense_seed_max)
    _uint(row["PeelSeed"], source, line, "PeelSeed", 0, 255)
    dense_failures = _uint(row["DenseFailures"], source, line, "DenseFailures")
    _uint(row["PeelScore"], source, line, "PeelScore")
    _uint(row["CurrentDenseCount"], source, line, "CurrentDenseCount", 1, n)
    _uint(
        row["CurrentDenseSeed"], source, line, "CurrentDenseSeed", 0,
        dense_seed_max,
    )
    _uint(row["CurrentPeelSeed"], source, line, "CurrentPeelSeed", 0, 255)
    trials = _uint(row["Trials"], source, line, "Trials", 1, (1 << 31) - 1)
    if dense_failures > trials:
        raise ValidationError(f"{source}:{line}: DenseFailures exceeds Trials")
    base_seed = _uint(row["BaseSeed"], source, line, "BaseSeed", 0, MASK64)
    trial_seed = _uint(row["TrialSeed"], source, line, "TrialSeed", 0, MASK64)
    if trial_seed != derive_trial_seed(base_seed, n, "small"):
        raise ValidationError(f"{source}:{line}: TrialSeed does not match BaseSeed/N")
    tables_hash = _uint(
        row["SmallTablesHash"], source, line, "SmallTablesHash", 0, MASK64
    )
    if row["MethodVersion"] != "small-v2":
        raise ValidationError(f"{source}:{line}: unsupported MethodVersion")
    return n, (trials, base_seed, tables_hash, row["MethodVersion"])


def _validate_dense_count(row, source, line):
    n = _uint(row["N"], source, line, "N", 2, 64000)
    dense_count = _uint(row["DenseCount"], source, line, "DenseCount", 1, n)
    failures = _uint(row["SelectedFailures"], source, line, "SelectedFailures")
    trials = _uint(row["Trials"], source, line, "Trials", 1, (1 << 31) - 1)
    max_failures = _uint(row["MaxFailures"], source, line, "MaxFailures")
    low_count_run = _uint(row["LowCountRun"], source, line, "LowCountRun", 1)
    base_seed = _uint(row["BaseSeed"], source, line, "BaseSeed", 0, MASK64)
    trial_seed = _uint(row["TrialSeed"], source, line, "TrialSeed", 0, MASK64)
    qualifying = _uint(row["QualifyingCount"], source, line, "QualifyingCount")
    source_hash = _uint(
        row["SourceTablesHash"], source, line, "SourceTablesHash", 0, MASK64
    )
    if failures > trials or max_failures > trials:
        raise ValidationError(
            f"{source}:{line}: failure threshold/count exceeds Trials"
        )
    if trial_seed != derive_trial_seed(base_seed, n, "dense-count"):
        raise ValidationError(f"{source}:{line}: TrialSeed does not match BaseSeed/N")
    if row["MethodVersion"] != "dense-count-v2":
        raise ValidationError(f"{source}:{line}: unsupported MethodVersion")
    status = row["Status"]
    if status == "selected-knee":
        if qualifying < low_count_run or failures > max_failures:
            raise ValidationError(f"{source}:{line}: inconsistent selected-knee row")
    elif status == "fallback-minimum":
        if qualifying >= low_count_run:
            raise ValidationError(f"{source}:{line}: inconsistent fallback row")
    else:
        raise ValidationError(f"{source}:{line}: unsupported Status={status!r}")
    return (base_seed, n), (
        trials, max_failures, low_count_run, source_hash, row["MethodVersion"]
    )


def _validate_dense_seed(row, source, line):
    index = _uint(row["DenseIndex"], source, line, "DenseIndex", 2, 99)
    used = _uint(row["Used"], source, line, "Used", 0, 1)
    n = _uint(row["N"], source, line, "N", 0, 64000)
    dense_count = _uint(row["DenseCount"], source, line, "DenseCount", 10, 398)
    if dense_count != index * 4 + 2:
        raise ValidationError(
            f"{source}:{line}: DenseCount does not map to DenseIndex"
        )
    _uint(row["DenseSeed"], source, line, "DenseSeed", 0, 255)
    failures = _uint(row["Failures"], source, line, "Failures")
    _uint(row["CurrentDenseSeed"], source, line, "CurrentDenseSeed", 0, 255)
    trials = _uint(row["Trials"], source, line, "Trials", 1, (1 << 31) - 1)
    if failures > trials:
        raise ValidationError(f"{source}:{line}: Failures exceeds Trials")
    base_seed = _uint(row["BaseSeed"], source, line, "BaseSeed", 0, MASK64)
    trial_seed = _uint(row["TrialSeed"], source, line, "TrialSeed", 0, MASK64)
    if trial_seed != derive_trial_seed(base_seed, index, "dense-seeds"):
        raise ValidationError(
            f"{source}:{line}: TrialSeed does not match BaseSeed/DenseIndex"
        )
    graph_hash = _uint(
        row["SourceTablesHash"], source, line, "SourceTablesHash", 0, MASK64
    )
    if row["MethodVersion"] != "most-dense-v2":
        raise ValidationError(f"{source}:{line}: unsupported MethodVersion")
    if used and n < 2048:
        raise ValidationError(f"{source}:{line}: used dense index has invalid N={n}")
    if not used and n != 0:
        raise ValidationError(f"{source}:{line}: unused dense index must have N=0")
    if not used and failures != 0:
        raise ValidationError(f"{source}:{line}: unused dense index has failures")
    return index, (trials, base_seed, graph_hash, row["MethodVersion"])


def _validate_peel(row, source, line):
    subdivision = _uint(
        row["Subdivision"], source, line, "Subdivision", 0, 2047
    )
    first_n = _uint(row["FirstN"], source, line, "FirstN", 0, 64000)
    last_n = _uint(row["LastN"], source, line, "LastN", 0, 64000)
    _uint(row["PeelSeed"], source, line, "PeelSeed", 0, 255)
    _uint(row["Score"], source, line, "Score", 0, (1 << 32) - 1)
    _uint(row["CurrentPeelSeed"], source, line, "CurrentPeelSeed", 0, 255)
    trials = _uint(row["Trials"], source, line, "Trials", 1, (1 << 31) - 1)
    max_tries = _uint(row["MaxTries"], source, line, "MaxTries", 1)
    evaluated = _uint(
        row["EvaluatedCandidates"], source, line, "EvaluatedCandidates", 0, 256
    )
    base_seed = _uint(row["BaseSeed"], source, line, "BaseSeed", 0, MASK64)
    requested_nlo = _uint(
        row["RequestedNLo"], source, line, "RequestedNLo", 2048, 64000
    )
    requested_nhi = _uint(
        row["RequestedNHi"], source, line, "RequestedNHi", 2048, 64000
    )
    if requested_nlo > requested_nhi:
        raise ValidationError(f"{source}:{line}: requested N range is reversed")
    base_tables_hash = _uint(
        row["BaseTablesHash"], source, line, "BaseTablesHash", 0, MASK64
    )
    if row["MethodVersion"] != "peel-v2":
        raise ValidationError(f"{source}:{line}: unsupported MethodVersion")
    status = row["Status"]
    statuses = {"searched", "retained-no-candidate", "retained-no-N"}
    if status not in statuses:
        raise ValidationError(f"{source}:{line}: unsupported Status={status!r}")
    if status == "retained-no-N":
        if first_n != 0 or last_n != 0 or evaluated != 0:
            raise ValidationError(
                f"{source}:{line}: retained-no-N row has evaluated N/candidates"
            )
    else:
        if first_n < requested_nlo or last_n > requested_nhi or last_n < first_n:
            raise ValidationError(
                f"{source}:{line}: invalid first/last N for subdivision"
            )
        if first_n % 2048 != subdivision or last_n % 2048 != subdivision:
            raise ValidationError(
                f"{source}:{line}: first/last N is not on the subdivision lattice"
            )
    return subdivision, (
        trials, max_tries, base_seed, requested_nlo, requested_nhi,
        base_tables_hash, row["MethodVersion"],
    )


VALIDATORS = {
    "dense-count": _validate_dense_count,
    "small": _validate_small,
    "dense-seeds": _validate_dense_seed,
    "peel": _validate_peel,
}

FULL_KEYS = {
    "small": set(range(2, 2048)),
    "dense-seeds": set(range(2, 100)),
    "peel": set(range(0, 2048)),
}


def canonical_dense_count_ns():
    values = []
    n = 2
    while n <= 64000:
        values.append(n)
        if n == 64000:
            break
        n = min(64000, n + max(1, n // 64))
    return values


def validate_and_merge(
    kind, inputs, require_full=False, allow_retained=False, min_runs=1
):
    records = {}
    methodology = None
    for value in inputs:
        path = _resolve_input(kind, value)
        for row, source, line in _read_rows(kind, path):
            key, method = VALIDATORS[kind](row, source, line)
            if key in records:
                previous = records[key][1]
                raise ValidationError(
                    f"duplicate {kind} key {key}: {previous} and {source}:{line}"
                )
            if methodology is None:
                methodology = method
            elif method != methodology:
                raise ValidationError(
                    f"{source}:{line}: methodology {method!r} does not match "
                    f"{methodology!r}"
                )
            records[key] = (row, f"{source}:{line}")
    if not records:
        raise ValidationError(f"no {kind} records")
    if require_full:
        if kind == "dense-count":
            expected_ns = set(canonical_dense_count_ns())
            seeds = sorted({key[0] for key in records})
            if len(seeds) < min_runs:
                raise ValidationError(
                    f"dense-count has {len(seeds)} independent runs; need {min_runs}"
                )
            for seed in seeds:
                actual_ns = {key[1] for key in records if key[0] == seed}
                if actual_ns != expected_ns:
                    missing = sorted(expected_ns - actual_ns)
                    extra = sorted(actual_ns - expected_ns)
                    raise ValidationError(
                        f"incomplete dense-count seed {seed}: missing={missing[:12]} "
                        f"({len(missing)} total) extra={extra[:12]} ({len(extra)} total)"
                    )
            fallback = [
                key for key, (row, _) in records.items()
                if row["Status"] == "fallback-minimum" and key[1] >= 2048
            ]
            if fallback and not allow_retained:
                raise ValidationError(
                    f"dense-count has large-N fallback rows: {fallback[:12]}"
                )
        else:
            expected = FULL_KEYS[kind]
            actual = set(records)
            if actual != expected:
                missing = sorted(expected - actual)
                extra = sorted(actual - expected)
                raise ValidationError(
                    f"incomplete {kind} coverage: missing={missing[:12]} "
                    f"({len(missing)} total) extra={extra[:12]} ({len(extra)} total)"
                )
        if kind == "peel":
            retained_no_n = [
                key for key, (row, _) in records.items()
                if row["Status"] == "retained-no-N"
            ]
            if retained_no_n:
                raise ValidationError(
                    f"full peel coverage contains retained-no-N rows: {retained_no_n[:12]}"
                )
            expected_last = lambda subdivision: (
                2048 + subdivision +
                ((64000 - (2048 + subdivision)) // 2048) * 2048
            )
            incomplete = [
                key for key, (row, _) in records.items()
                if int(row["RequestedNLo"]) != 2048
                or int(row["RequestedNHi"]) != 64000
                or int(row["FirstN"]) != 2048 + key
                or int(row["LastN"]) != expected_last(key)
            ]
            if incomplete:
                raise ValidationError(
                    f"full peel rows do not cover the complete lattice: {incomplete[:12]}"
                )
            retained = [
                key for key, (row, _) in records.items()
                if row["Status"] != "searched"
            ]
            if retained and not allow_retained:
                raise ValidationError(
                    f"full peel coverage has unsearched retained rows: {retained[:12]}"
                )
    return [records[key][0] for key in sorted(records)]


def write_rows(kind, rows, output=None):
    schema = SCHEMAS[kind]
    if output is None:
        writer = csv.DictWriter(
            sys.stdout, fieldnames=schema, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
        return
    path = Path(output)
    temporary = Path(str(path) + ".tmp")
    if path.exists() or temporary.exists():
        raise ValidationError(f"output already exists: {path}")
    try:
        with temporary.open("x", encoding="ascii", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=schema, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
            handle.flush()
            os.fsync(handle.fileno())
        os.rename(temporary, path)
    except Exception:
        # Preserve a partial temporary file as interruption evidence.
        raise


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("kind", choices=sorted(SCHEMAS))
    parser.add_argument("inputs", nargs="+")
    parser.add_argument("--full", action="store_true", help="require full domain coverage")
    parser.add_argument(
        "--allow-retained", action="store_true",
        help="permit explicit retained peel rows during full validation",
    )
    parser.add_argument("--output", help="write canonical merged TSV atomically")
    parser.add_argument(
        "--min-runs", type=int, default=1,
        help="minimum independent dense-count base seeds for full validation",
    )
    args = parser.parse_args(argv)
    try:
        rows = validate_and_merge(
            args.kind, args.inputs, args.full, args.allow_retained, args.min_runs
        )
        write_rows(args.kind, rows, args.output)
    except (OSError, ValidationError) as exc:
        print(exc, file=sys.stderr)
        return 1
    print(f"validated {len(rows)} {args.kind} rows", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
