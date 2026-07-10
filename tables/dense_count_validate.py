#!/usr/bin/env python3
"""Validate dense-count shard and aggregate records before consumption."""

import argparse
import csv
import math
import re
import sys
from pathlib import Path


SHARD_HEADER = ["N", "DenseCount", "LowestFailures"]
AGGREGATE_HEADER = [
    "N",
    "Observations",
    "MinGeneratedDense",
    "MaxGeneratedDense",
    "LowestFailuresAtMax",
    "SourcesAtMax",
]
MIN_N = 2
MAX_N = 64000
MAX_DENSE = 400


class ValidationError(ValueError):
    pass


def _uint(text, source, line, name, minimum=0, maximum=None):
    if not re.fullmatch(r"[0-9]+", text):
        raise ValidationError(f"{source}:{line}: {name} is not an unsigned integer: {text!r}")
    value = int(text, 10)
    if value < minimum or (maximum is not None and value > maximum):
        upper = "unbounded" if maximum is None else str(maximum)
        raise ValidationError(
            f"{source}:{line}: {name}={value} outside [{minimum}, {upper}]"
        )
    return value


def _rate(text, source, line, name):
    try:
        value = float(text)
    except ValueError as exc:
        raise ValidationError(f"{source}:{line}: {name} is not numeric: {text!r}") from exc
    if not math.isfinite(value) or value < 0.0 or value > 1.0:
        raise ValidationError(f"{source}:{line}: {name}={text!r} outside finite [0, 1]")
    return value


def _manifest(path):
    values = {}
    manifest = path.parent / "manifest.txt"
    if manifest.is_file():
        for raw in manifest.read_text(encoding="utf-8").splitlines():
            if "=" not in raw or raw.startswith("#"):
                continue
            key, value = raw.split("=", 1)
            values[key.strip()] = value.strip()

    name = path.parent.name
    trial_match = re.search(r"(?:^|_)t([0-9]+)(?:_|$)", name)
    seed_match = re.search(r"(?:^|_)s([0-9]+)(?:_|$)", name)
    trials = values.get("dcount_trials")
    if trials is None and trial_match:
        trials = trial_match.group(1)
    seed = values.get("dcount_seed", values.get("seed"))
    if seed is None and seed_match:
        seed = seed_match.group(1)
    return {
        "trials": trials,
        "max_failures": values.get("dcount_max_failures"),
        "low_count_run": values.get("dcount_low_count_run"),
        "seed": seed,
    }


def _open_tsv(path_text):
    if path_text == "-":
        return sys.stdin, "<stdin>", None
    path = Path(path_text)
    if path.is_dir():
        path = path / "dense_count.out"
    if not path.is_file():
        raise ValidationError(f"missing dense-count input: {path}")
    return path.open("r", encoding="utf-8", newline=""), str(path), path


def read_shard(path_text):
    handle, source, path = _open_tsv(path_text)
    close = handle is not sys.stdin
    try:
        reader = csv.reader(handle, delimiter="\t", strict=True)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValidationError(f"{source}: empty dense-count shard") from exc
        if header != SHARD_HEADER:
            raise ValidationError(
                f"{source}:1: expected tab-separated header {SHARD_HEADER!r}, got {header!r}"
            )
        rows = []
        seen = set()
        for line, fields in enumerate(reader, 2):
            if fields == [] or fields == [""]:
                continue
            if len(fields) != 3:
                raise ValidationError(f"{source}:{line}: expected 3 fields, got {len(fields)}")
            n = _uint(fields[0], source, line, "N", MIN_N, MAX_N)
            dense = _uint(fields[1], source, line, "DenseCount", 1, MAX_DENSE)
            _rate(fields[2], source, line, "LowestFailures")
            if n in seen:
                raise ValidationError(f"{source}:{line}: duplicate N={n} within shard")
            seen.add(n)
            rows.append((n, dense, fields[2], source))
        if not rows:
            raise ValidationError(f"{source}: no dense-count rows found")
        return rows, (_manifest(path) if path else {})
    except csv.Error as exc:
        raise ValidationError(f"{source}: malformed TSV: {exc}") from exc
    finally:
        if close:
            handle.close()


def validate_shards(paths):
    all_rows = []
    methods = {}
    seen_seed_n = {}
    for path_text in paths:
        rows, metadata = read_shard(path_text)
        source = rows[0][3]
        for key in ("trials", "max_failures", "low_count_run"):
            value = metadata.get(key)
            if value is None:
                continue
            previous = methods.setdefault(key, value)
            if previous != value:
                raise ValidationError(
                    f"incompatible dense-count methodology: {key}={previous} and {value} ({source})"
                )
        seed = metadata.get("seed")
        if seed is not None:
            for n, _, _, _ in rows:
                key = (seed, n)
                if key in seen_seed_n:
                    raise ValidationError(
                        f"reused seed {seed} overlaps N={n}: {seen_seed_n[key]} and {source}"
                    )
                seen_seed_n[key] = source
        all_rows.extend(rows)
    return all_rows


def read_aggregate(path_text):
    handle, source, _ = _open_tsv(path_text)
    close = handle is not sys.stdin
    try:
        reader = csv.reader(handle, delimiter="\t", strict=True)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise ValidationError(f"{source}: empty dense-count aggregate") from exc
        if header != AGGREGATE_HEADER:
            raise ValidationError(
                f"{source}:1: expected tab-separated header {AGGREGATE_HEADER!r}, got {header!r}"
            )
        rows = []
        seen = set()
        for line, fields in enumerate(reader, 2):
            if fields == [] or fields == [""]:
                continue
            if len(fields) != 6:
                raise ValidationError(f"{source}:{line}: expected 6 fields, got {len(fields)}")
            n = _uint(fields[0], source, line, "N", MIN_N, MAX_N)
            observations = _uint(fields[1], source, line, "Observations", 1)
            min_dense = _uint(fields[2], source, line, "MinGeneratedDense", 1, MAX_DENSE)
            max_dense = _uint(fields[3], source, line, "MaxGeneratedDense", 1, MAX_DENSE)
            _rate(fields[4], source, line, "LowestFailuresAtMax")
            if min_dense > max_dense:
                raise ValidationError(f"{source}:{line}: min dense exceeds max dense")
            if not fields[5]:
                raise ValidationError(f"{source}:{line}: SourcesAtMax is empty")
            if n in seen:
                raise ValidationError(f"{source}:{line}: duplicate N={n} within aggregate")
            seen.add(n)
            rows.append((n, observations, min_dense, max_dense, fields[4], fields[5]))
        if not rows:
            raise ValidationError(f"{source}: no dense-count aggregate rows found")
        return rows
    except csv.Error as exc:
        raise ValidationError(f"{source}: malformed TSV: {exc}") from exc
    finally:
        if close:
            handle.close()


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("kind", choices=("shards", "aggregate"))
    parser.add_argument("inputs", nargs="+")
    args = parser.parse_args(argv)
    try:
        if args.kind == "shards":
            for row in validate_shards(args.inputs):
                print("\t".join(map(str, row)))
        else:
            for input_name in args.inputs:
                for row in read_aggregate(input_name):
                    print("\t".join(map(str, row)))
    except (OSError, ValidationError) as exc:
        print(exc, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
