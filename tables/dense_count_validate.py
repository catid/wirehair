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
DCOUNT_SCHEMA_VERSION = 1
DCOUNT_METHODOLOGY = "wirehair-generate-dense-count-v1"
UINT_MAX = (1 << 32) - 1
INT_MAX = (1 << 31) - 1
UINT64_MAX = (1 << 64) - 1
REQUIRED_MULTI_SHARD_METADATA = (
    "schema",
    "methodology",
    "seed",
    "trials",
    "max_failures",
    "low_count_run",
    "nlo",
    "nhi",
)


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
        for line, raw in enumerate(
                manifest.read_text(encoding="utf-8").splitlines(), 1):
            stripped = raw.strip()
            if not stripped or stripped.startswith("#"):
                continue
            # run_regeneration_slice.sh appends sha256sum records after
            # the key/value campaign metadata.  Check these first because
            # an output path may itself contain '='.
            if re.fullmatch(r"\\?[0-9a-fA-F]{64} [ *].+", raw):
                continue
            if "=" not in raw:
                raise ValidationError(
                    f"{manifest}:{line}: malformed manifest entry: {raw!r}"
                )
            key, value = raw.split("=", 1)
            key = key.strip()
            if not key:
                raise ValidationError(f"{manifest}:{line}: empty manifest key")
            if key in values:
                raise ValidationError(
                    f"{manifest}:{line}: duplicate manifest key {key!r}"
                )
            values[key] = (value.strip(), line)

    name = path.parent.name
    trial_match = re.search(r"(?:^|_)t([0-9]+)(?:_|$)", name)
    seed_match = re.search(r"(?:^|_)s([0-9]+)(?:_|$)", name)

    def text(key):
        entry = values.get(key)
        return entry[0] if entry is not None else None

    def number(key, name, minimum=0, maximum=None):
        entry = values.get(key)
        if entry is None:
            return None
        value, line = entry
        return _uint(value, str(manifest), line, name, minimum, maximum)

    trials = number("dcount_trials", "dcount_trials", 1, INT_MAX)
    if trials is None and trial_match:
        trials = _uint(
            trial_match.group(1), str(path.parent), 1,
            "directory-name dcount_trials", 1, INT_MAX,
        )
    seed = number("dcount_seed", "dcount_seed", 0, UINT64_MAX)
    if seed is None:
        seed = number("seed", "seed", 0, UINT64_MAX)
    if seed is None and seed_match:
        seed = _uint(
            seed_match.group(1), str(path.parent), 1,
            "directory-name dcount_seed", 0, UINT64_MAX,
        )

    schema = number("dcount_schema", "dcount_schema", 1, UINT_MAX)
    if schema is not None and schema != DCOUNT_SCHEMA_VERSION:
        line = values["dcount_schema"][1]
        raise ValidationError(
            f"{manifest}:{line}: unsupported dcount_schema={schema}; "
            f"expected {DCOUNT_SCHEMA_VERSION}"
        )
    methodology = text("dcount_methodology")
    if methodology is not None and methodology != DCOUNT_METHODOLOGY:
        line = values["dcount_methodology"][1]
        raise ValidationError(
            f"{manifest}:{line}: unsupported dcount_methodology={methodology!r}; "
            f"expected {DCOUNT_METHODOLOGY!r}"
        )

    return {
        "schema": schema,
        "methodology": methodology,
        "trials": trials,
        "max_failures": number(
            "dcount_max_failures", "dcount_max_failures", 0, UINT_MAX
        ),
        "low_count_run": number(
            "dcount_low_count_run", "dcount_low_count_run", 1, UINT_MAX
        ),
        "seed": seed,
        "nlo": number("dcount_nlo", "dcount_nlo", MIN_N, MAX_N),
        "nhi": number("dcount_nhi", "dcount_nhi", MIN_N, MAX_N),
        "provided": {
            field for field, key in (
                ("schema", "dcount_schema"),
                ("methodology", "dcount_methodology"),
                ("seed", "dcount_seed"),
                ("trials", "dcount_trials"),
                ("max_failures", "dcount_max_failures"),
                ("low_count_run", "dcount_low_count_run"),
                ("nlo", "dcount_nlo"),
                ("nhi", "dcount_nhi"),
            ) if key in values
        },
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


def _input_identity(path_text):
    if path_text == "-":
        return None, "<stdin>"
    path = Path(path_text)
    if path.is_dir():
        path = path / "dense_count.out"
    if not path.is_file():
        raise ValidationError(f"missing dense-count input: {path}")
    resolved = path.resolve()
    stat = resolved.stat()
    return (stat.st_dev, stat.st_ino), str(resolved)


def _reject_repeated_inputs(paths, kind):
    seen = {}
    for path_text in paths:
        identity, source = _input_identity(path_text)
        if identity is None:
            if len(paths) > 1:
                raise ValidationError(
                    f"standard input cannot be combined with other {kind} inputs"
                )
            continue
        if identity in seen:
            raise ValidationError(
                f"repeated {kind} input: {seen[identity]} and {source}"
            )
        seen[identity] = source


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
    _reject_repeated_inputs(paths, "dense-count shard")
    all_rows = []
    methods = {}
    seen_seed_n = {}
    seed_ranges = {}
    require_complete_metadata = len(paths) > 1
    for path_text in paths:
        rows, metadata = read_shard(path_text)
        source = rows[0][3]
        if require_complete_metadata:
            missing = [
                key for key in REQUIRED_MULTI_SHARD_METADATA
                if key not in metadata["provided"]
            ]
            if missing:
                raise ValidationError(
                    f"{source}: multi-shard input requires complete provenance; "
                    f"missing {', '.join(missing)}"
                )

        nlo = metadata.get("nlo")
        nhi = metadata.get("nhi")
        if (nlo is None) != (nhi is None):
            raise ValidationError(
                f"{source}: dcount_nlo and dcount_nhi must be provided together"
            )
        if nlo is not None:
            if nlo > nhi:
                raise ValidationError(
                    f"{source}: dense-count provenance range {nlo}..{nhi} is reversed"
                )
            for n, _, _, _ in rows:
                if n < nlo or n > nhi:
                    raise ValidationError(
                        f"{source}: N={n} outside declared provenance range {nlo}..{nhi}"
                    )

        for key in (
                "schema", "methodology", "trials", "max_failures",
                "low_count_run"):
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
            if nlo is not None:
                for previous_lo, previous_hi, previous_source in seed_ranges.setdefault(seed, []):
                    if nlo <= previous_hi and previous_lo <= nhi:
                        raise ValidationError(
                            f"reused seed {seed} has overlapping provenance ranges "
                            f"{previous_lo}..{previous_hi} ({previous_source}) and "
                            f"{nlo}..{nhi} ({source})"
                        )
                seed_ranges[seed].append((nlo, nhi, source))
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


def validate_aggregates(paths):
    _reject_repeated_inputs(paths, "dense-count aggregate")
    all_rows = []
    seen_n = {}
    for path_text in paths:
        rows = read_aggregate(path_text)
        source = "<stdin>" if path_text == "-" else str(Path(path_text))
        for row in rows:
            n = row[0]
            if n in seen_n:
                raise ValidationError(
                    f"duplicate aggregate N={n} across {seen_n[n]} and {source}; "
                    "combine the provenance-bearing shards before aggregation"
                )
            seen_n[n] = source
        all_rows.extend(rows)
    return all_rows


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
            for row in validate_aggregates(args.inputs):
                print("\t".join(map(str, row)))
    except (OSError, ValidationError) as exc:
        print(exc, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
