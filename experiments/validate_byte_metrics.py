#!/usr/bin/env python3
"""Validate legacy and schema-v2 experiment byte-accounting CSV files."""

import argparse
import csv
import math
import re
import sys
from pathlib import Path


XOR_V1 = (
    "block_bytes", "working_blocks", "ops_per_repeat", "repeats",
    "median_total_us", "median_us_per_xor", "median_gib_per_s", "checksum",
)
MULADD_V1 = (
    "block_bytes", "pool_blocks", "ops_per_repeat", "repeats",
    "xor_median_us_per_op", "xor_median_gib_per_s",
    "muladd_median_us_per_op", "muladd_median_gib_per_s",
    "muladd_xor_ratio", "checksum",
)
FANIN_V1 = (
    "k", "block_bytes", "pool_blocks", "ops_per_repeat", "repeats",
    "chained_median_us_per_op", "chained_logical_gib_per_s",
    "gather_median_us_per_op", "gather_logical_gib_per_s",
    "gather_chained_ratio", "checksum",
)
TILING_SYNTH_V1 = (
    "case", "block_bytes", "blocks", "ops", "mode", "tile_bytes",
    "median_ms", "logical_gib", "memory_gib_3stream", "gib_per_s", "checksum",
)
TILING_TRACE_V1 = (
    "case", "block_bytes", "recovery_blocks", "input_blocks", "ops", "mode",
    "tile_bytes", "median_ms", "logical_gib", "memory_gib_3stream",
    "gib_per_s", "checksum",
)

XOR_V2 = XOR_V1 + (
    "schema_version", "logical_bytes_per_op", "estimated_read_bytes_per_op",
    "estimated_write_bytes_per_op", "logical_gib_per_s",
    "estimated_traffic_gib_per_s",
)
MULADD_V2 = MULADD_V1 + (
    "schema_version", "logical_bytes_per_op", "estimated_read_bytes_per_op",
    "estimated_write_bytes_per_op", "xor_logical_gib_per_s",
    "xor_estimated_traffic_gib_per_s", "muladd_logical_gib_per_s",
    "muladd_estimated_traffic_gib_per_s",
)
FANIN_V2 = (
    "k", "block_bytes", "pool_blocks", "ops_per_repeat", "repeats",
    "chained_median_us_per_op", "chained_legacy_normalized_gib_per_s",
    "gather_median_us_per_op", "gather_legacy_normalized_gib_per_s",
    "gather_chained_ratio", "checksum", "schema_version", "logical_bytes_per_op",
    "chained_estimated_read_bytes_per_op",
    "chained_estimated_write_bytes_per_op", "gather_estimated_read_bytes_per_op",
    "gather_estimated_write_bytes_per_op", "chained_logical_gib_per_s",
    "chained_estimated_traffic_gib_per_s", "gather_logical_gib_per_s",
    "gather_estimated_traffic_gib_per_s",
)
TILING_SUFFIX_V2 = (
    "schema_version", "logical_work_gib", "estimated_read_gib",
    "estimated_write_gib", "estimated_traffic_gib", "logical_work_gib_per_s",
    "estimated_traffic_gib_per_s",
)
TILING_SYNTH_V2 = TILING_SYNTH_V1 + TILING_SUFFIX_V2
TILING_TRACE_V2 = TILING_TRACE_V1 + TILING_SUFFIX_V2

SCHEMAS = {
    XOR_V1: ("xor", 1),
    MULADD_V1: ("muladd", 1),
    FANIN_V1: ("fanin", 1),
    TILING_SYNTH_V1: ("tiling-synthetic", 1),
    TILING_TRACE_V1: ("tiling-trace", 1),
    XOR_V2: ("xor", 2),
    MULADD_V2: ("muladd", 2),
    FANIN_V2: ("fanin", 2),
    TILING_SYNTH_V2: ("tiling-synthetic", 2),
    TILING_TRACE_V2: ("tiling-trace", 2),
}
TEXT_FIELDS = {"case", "mode", "checksum"}
INTEGER_FIELDS = {
    "k", "block_bytes", "working_blocks", "pool_blocks", "ops_per_repeat",
    "repeats", "blocks", "recovery_blocks", "input_blocks", "ops", "tile_bytes",
    "schema_version", "logical_bytes_per_op", "estimated_read_bytes_per_op",
    "estimated_write_bytes_per_op", "chained_estimated_read_bytes_per_op",
    "chained_estimated_write_bytes_per_op", "gather_estimated_read_bytes_per_op",
    "gather_estimated_write_bytes_per_op",
}


class ValidationError(ValueError):
    pass


def integer(text, source, name, minimum=0):
    if not re.fullmatch(r"[0-9]+", text or ""):
        raise ValidationError("%s: %s is not an unsigned integer" % (source, name))
    value = int(text)
    if value < minimum or value > (1 << 63) - 1:
        raise ValidationError("%s: %s outside supported range" % (source, name))
    return value


def number(text, source, name):
    try:
        value = float(text)
    except ValueError as exc:
        raise ValidationError("%s: %s is not numeric" % (source, name)) from exc
    if not math.isfinite(value) or value < 0.0:
        raise ValidationError("%s: %s is not finite nonnegative" % (source, name))
    return value


def close(actual, expected, relative=0.006, absolute=0.002):
    return abs(actual - expected) <= max(absolute, relative * max(abs(expected), 1.0))


def rate(bytes_count, usec):
    return (bytes_count / (1024.0 ** 3)) / (usec / 1000000.0)


def validate_v2(kind, values, source):
    if values["schema_version"] != 2:
        raise ValidationError("%s: schema_version must be 2" % source)
    block = values["block_bytes"]
    if kind == "xor":
        expected = (block, 2 * block, block)
        actual = tuple(values[name] for name in (
            "logical_bytes_per_op", "estimated_read_bytes_per_op",
            "estimated_write_bytes_per_op"))
        if actual != expected:
            raise ValidationError("%s: XOR byte ledger mismatch" % source)
        usec = values["median_us_per_xor"]
        if (usec <= 0.0 or
                not close(values["logical_gib_per_s"], rate(block, usec)) or
                not close(values["estimated_traffic_gib_per_s"], rate(3 * block, usec))):
            raise ValidationError("%s: XOR throughput identity mismatch" % source)
    elif kind == "muladd":
        actual = tuple(values[name] for name in (
            "logical_bytes_per_op", "estimated_read_bytes_per_op",
            "estimated_write_bytes_per_op"))
        if actual != (block, 2 * block, block):
            raise ValidationError("%s: muladd byte ledger mismatch" % source)
        for prefix in ("xor", "muladd"):
            usec = values[prefix + "_median_us_per_op"]
            if (usec <= 0.0 or
                    not close(values[prefix + "_logical_gib_per_s"], rate(block, usec)) or
                    not close(values[prefix + "_estimated_traffic_gib_per_s"],
                              rate(3 * block, usec))):
                raise ValidationError("%s: %s throughput identity mismatch" %
                                      (source, prefix))
    elif kind == "fanin":
        fanin = values["k"]
        expected = (
            block, 2 * fanin * block, fanin * block,
            (fanin + 1) * block, block,
        )
        names = (
            "logical_bytes_per_op", "chained_estimated_read_bytes_per_op",
            "chained_estimated_write_bytes_per_op",
            "gather_estimated_read_bytes_per_op",
            "gather_estimated_write_bytes_per_op",
        )
        if tuple(values[name] for name in names) != expected:
            raise ValidationError("%s: fan-in byte ledger mismatch" % source)
        for prefix, traffic in (("chained", 3 * fanin * block),
                                ("gather", (fanin + 2) * block)):
            usec = values[prefix + "_median_us_per_op"]
            if (usec <= 0.0 or
                    not close(values[prefix + "_logical_gib_per_s"], rate(block, usec)) or
                    not close(values[prefix + "_estimated_traffic_gib_per_s"],
                              rate(traffic, usec))):
                raise ValidationError("%s: %s throughput identity mismatch" %
                                      (source, prefix))
    else:
        if not close(values["logical_work_gib"], values["logical_gib"], 1e-6, 2e-6):
            raise ValidationError("%s: tiling logical-work mismatch" % source)
        traffic = values["estimated_read_gib"] + values["estimated_write_gib"]
        if not close(values["estimated_traffic_gib"], traffic, 1e-6, 2e-6):
            raise ValidationError("%s: tiling read/write traffic mismatch" % source)
        seconds = values["median_ms"] / 1000.0
        if (seconds <= 0.0 or
                not close(values["logical_work_gib_per_s"],
                          values["logical_work_gib"] / seconds) or
                not close(values["estimated_traffic_gib_per_s"],
                          values["estimated_traffic_gib"] / seconds)):
            raise ValidationError("%s: tiling throughput identity mismatch" % source)
        if kind == "tiling-synthetic":
            if (not close(values["estimated_read_gib"], 2 * values["logical_gib"],
                          1e-6, 2e-6) or
                    not close(values["estimated_write_gib"], values["logical_gib"],
                              1e-6, 2e-6)):
                raise ValidationError("%s: synthetic XOR ledger mismatch" % source)


def validate_file(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        lines = [line for line in handle if not line.startswith("#")]
    reader = csv.DictReader(lines, strict=True)
    header = tuple(reader.fieldnames or ())
    if header not in SCHEMAS:
        raise ValidationError("%s: unsupported byte-metric CSV schema" % path)
    kind, version = SCHEMAS[header]
    rows = 0
    for line, row in enumerate(reader, 2):
        source = "%s:%d" % (path, line)
        if None in row or any(value is None for value in row.values()):
            raise ValidationError("%s: wrong field count" % source)
        values = {}
        for name, text in row.items():
            if name == "checksum":
                if not re.fullmatch(r"0x[0-9a-fA-F]{16}", text):
                    raise ValidationError("%s: malformed checksum" % source)
            elif name in TEXT_FIELDS:
                if not text:
                    raise ValidationError("%s: empty %s" % (source, name))
            elif name in INTEGER_FIELDS:
                values[name] = integer(text, source, name)
            else:
                values[name] = number(text, source, name)
        if version == 2:
            validate_v2(kind, values, source)
        rows += 1
    if rows == 0:
        raise ValidationError("%s: no data rows" % path)
    return kind, version, rows


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="+")
    args = parser.parse_args(argv)
    summaries = []
    try:
        for name in args.files:
            path = Path(name)
            if not path.is_file():
                raise ValidationError("missing byte-metric CSV: %s" % path)
            summaries.append((path,) + validate_file(path))
    except (OSError, csv.Error, ValidationError) as exc:
        print(exc, file=sys.stderr)
        return 1
    for path, kind, version, rows in summaries:
        print("%s: kind=%s schema=v%d rows=%d" % (path, kind, version, rows))
    return 0


if __name__ == "__main__":
    sys.exit(main())
