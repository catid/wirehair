#!/usr/bin/env python3
"""Verify regeneration provenance and render candidate Wirehair base tables."""

import argparse
import csv
import re
import sys
from pathlib import Path

try:
    from . import atomic_publish
    from . import regeneration_results as results
except ImportError:
    import atomic_publish
    import regeneration_results as results


MASK64 = (1 << 64) - 1
FNV_OFFSET = 1469598103934665603
FNV_PRIME = 1099511628211

ARRAY_LENGTHS = {
    "kTinyDenseCounts": 65,
    "kTinyDenseSeeds": 65,
    "kSmallDenseSeeds": 1983,
    "kSmallPeelSeeds": 2048,
    "kDenseSeeds": 100,
    "kPeelSeeds": 2048,
}


class ApplyError(ValueError):
    pass


def fnv_add(hash_value, value):
    return ((hash_value ^ value) * FNV_PRIME) & MASK64


def _unique_match(pattern, text, description):
    matches = list(re.finditer(pattern, text, re.DOTALL))
    if len(matches) != 1:
        raise ApplyError(f"expected one {description}, found {len(matches)}")
    return matches[0]


def parse_numeric_array(text, name):
    pattern = (
        rf"(?P<prefix>(?:static\s+)?const\s+uint(?:8|16)_t\s+{name}"
        rf"\s*\[[^\]]+\]\s*=\s*\{{)(?P<body>.*?)(?P<suffix>\n\}};)"
    )
    match = _unique_match(pattern, text, name)
    body = re.sub(r"//[^\n]*", "", match.group("body"))
    values = [int(value, 0) for value in re.findall(r"0x[0-9a-fA-F]+|[0-9]+", body)]
    cleaned = re.sub(r"0x[0-9a-fA-F]+|[0-9]+|[\s,]", "", body)
    if cleaned:
        raise ApplyError(f"unsupported syntax in {name}: {cleaned[:40]!r}")
    expected = ARRAY_LENGTHS[name]
    if len(values) != expected:
        raise ApplyError(f"{name} has {len(values)} entries; expected {expected}")
    return values, match


def replace_numeric_array(text, name, values, columns=32):
    current, match = parse_numeric_array(text, name)
    if len(values) != len(current):
        raise ApplyError(f"replacement length mismatch for {name}")
    lines = []
    for start in range(0, len(values), columns):
        lines.append("    " + ",".join(str(value) for value in values[start:start + columns]) + ",")
    body = "\n" + "\n".join(lines)
    return text[:match.start()] + match.group("prefix") + body + match.group("suffix") + text[match.end():]


def parse_dense_points(text):
    count_match = _unique_match(
        r"static\s+const\s+unsigned\s+kDensePointCount\s*=\s*([0-9]+)\s*;",
        text,
        "kDensePointCount",
    )
    declared = int(count_match.group(1), 10)
    pattern = (
        r"(?P<prefix>static\s+const\s+DensePoint\s+kDensePoints"
        r"\s*\[[^\]]+\]\s*=\s*\{)(?P<body>.*?)(?P<suffix>\n\};)"
    )
    match = _unique_match(pattern, text, "kDensePoints")
    body = re.sub(r"//[^\n]*", "", match.group("body"))
    points = [
        (int(n, 10), int(count, 10))
        for n, count in re.findall(r"\{\s*([0-9]+)\s*,\s*([0-9]+)\s*\}\s*,", body)
    ]
    cleaned = re.sub(r"\{\s*[0-9]+\s*,\s*[0-9]+\s*\}\s*,|\s", "", body)
    if cleaned:
        raise ApplyError(f"unsupported syntax in kDensePoints: {cleaned[:40]!r}")
    if len(points) != declared:
        raise ApplyError(
            f"kDensePointCount={declared}, but parsed {len(points)} points"
        )
    return points, count_match, match


def replace_dense_points(text, points):
    if len(points) < 2 or points[0][0] != 2048 or points[-1][0] != 64000:
        raise ApplyError("dense points must span exactly N=2048..64000")
    if any(a[0] >= b[0] for a, b in zip(points, points[1:])):
        raise ApplyError("replacement dense points must be strictly increasing")
    _, _, array_match = parse_dense_points(text)
    body = "\n" + "\n".join(
        f"    {{ {n}, {count} }}," for n, count in points
    )
    text = (
        text[:array_match.start()] + array_match.group("prefix") + body +
        array_match.group("suffix") + text[array_match.end():]
    )
    count_match = _unique_match(
        r"static\s+const\s+unsigned\s+kDensePointCount\s*=\s*([0-9]+)\s*;",
        text,
        "kDensePointCount",
    )
    text = (
        text[:count_match.start(1)] + str(len(points)) +
        text[count_match.end(1):]
    )
    if len(parse_dense_points(text)[0]) != len(points):
        raise ApplyError("dense point replacement did not round-trip")
    return text


class SourceTables:
    def __init__(self, text):
        self.text = text
        self.tiny_counts, _ = parse_numeric_array(text, "kTinyDenseCounts")
        self.tiny_dense, _ = parse_numeric_array(text, "kTinyDenseSeeds")
        self.small_dense, _ = parse_numeric_array(text, "kSmallDenseSeeds")
        self.small_peel, _ = parse_numeric_array(text, "kSmallPeelSeeds")
        self.dense_seeds, _ = parse_numeric_array(text, "kDenseSeeds")
        self.peel_seeds, _ = parse_numeric_array(text, "kPeelSeeds")
        self.dense_points, _, _ = parse_dense_points(text)

    @staticmethod
    def _interpolate(n0, n1, c0, c1, n):
        # C++ integer division truncates toward zero.
        delta = int(((n - n0) * (c1 - c0)) / (n1 - n0))
        return c0 + delta

    def dense_count(self, n):
        if n < 65:
            return self.tiny_counts[n]
        if n < 2048:
            if n <= 500:
                raw = self._interpolate(64, 500, 26, 35, n)
            elif n <= 1000:
                raw = self._interpolate(500, 1000, 35, 48, n)
            else:
                raw = self._interpolate(1000, 2048, 48, 62, n)
        else:
            low = 0
            high = len(self.dense_points) - 1
            while True:
                mid = (high + low) // 2
                if mid == low:
                    break
                if n > self.dense_points[mid][0]:
                    low = mid
                else:
                    high = mid
            n0, c0 = self.dense_points[low]
            n1, c1 = self.dense_points[low + 1]
            if not n0 <= n <= n1:
                raise ApplyError(f"dense point graph does not cover N={n}")
            raw = self._interpolate(n0, n1, c0, c1, n)
        return raw + ((2 - raw) % 4)

    def dense_seed(self, n, dense_count=None):
        count = self.dense_count(n) if dense_count is None else dense_count
        if n < 65:
            return self.tiny_dense[n]
        if n < 2048:
            return self.small_dense[n - 65]
        return self.dense_seeds[count // 4]

    def peel_seed(self, n):
        if n < 2048:
            return self.small_peel[n]
        return self.peel_seeds[n % 2048]

    def small_hash(self):
        value = FNV_OFFSET
        for n in range(2, 2048):
            count = self.dense_count(n)
            for item in (n, count, self.dense_seed(n, count), self.peel_seed(n)):
                value = fnv_add(value, item)
        return value

    def dense_count_hash(self):
        value = FNV_OFFSET
        for n in range(2, 64001):
            value = fnv_add(value, n)
            value = fnv_add(value, self.dense_count(n))
        return value

    def dense_seed_hash(self):
        value = FNV_OFFSET
        for n in range(2048, 64001):
            value = fnv_add(value, n)
            value = fnv_add(value, self.dense_count(n))
        for seed in self.dense_seeds:
            value = fnv_add(value, seed)
        return value

    def dense_graph_hash(self):
        value = FNV_OFFSET
        for n in range(2048, 64001):
            value = fnv_add(value, n)
            value = fnv_add(value, self.dense_count(n))
        return value

    def peel_hash(self):
        value = FNV_OFFSET
        for n in range(2048, 64001):
            count = self.dense_count(n)
            for item in (n, count, self.dense_seed(n, count)):
                value = fnv_add(value, item)
        for seed in self.peel_seeds:
            value = fnv_add(value, seed)
        return value


def _read_validated(phase, paths, allow_retained=False):
    return results.validate_and_merge(
        phase,
        [str(path) for path in paths],
        require_full=True,
        allow_retained=allow_retained,
        min_runs=4 if phase == "dense-count" else 1,
    )


def _assert_single_hash(rows, field, expected):
    values = {int(row[field], 10) for row in rows}
    if values != {expected}:
        raise ApplyError(
            f"{field} does not match source: results={sorted(values)} source={expected}"
        )


def apply_dense_count(source, rows):
    tables = SourceTables(source)
    _assert_single_hash(rows, "SourceTablesHash", tables.dense_count_hash())
    by_n = {}
    ledger = []
    for row in rows:
        n = int(row["N"])
        by_n.setdefault(n, []).append(row)
    points = []
    for n in results.canonical_dense_count_ns():
        group = by_n[n]
        generated = [int(row["DenseCount"]) for row in group]
        chosen = max(generated)
        chosen += (2 - chosen) % 4
        if n >= 2048:
            if chosen > 400:
                raise ApplyError(f"normalized dense count exceeds 400 at N={n}")
            points.append((n, chosen))
        ledger.append({
            "Key": str(n),
            "Current": str(tables.dense_count(n)),
            "Candidate": str(chosen),
            "Decision": (
                "defer-small-generator" if n < 65 else
                "retain-fixed-interpolation" if n < 2048 else
                "replace" if chosen != tables.dense_count(n) else "retain-same"
            ),
            "Evidence": (
                f"max-over-{len(group)}-independent-runs;"
                + ("rendered-by-small-phase" if n < 65 else
                   "runtime-fixed-interpolation" if n < 2048 else
                   "rendered-dense-point")
            ),
        })
    return replace_dense_points(source, points), ledger


def apply_small(source, rows):
    tables = SourceTables(source)
    _assert_single_hash(rows, "SmallTablesHash", tables.small_hash())
    tiny_counts = list(tables.tiny_counts)
    tiny_dense = list(tables.tiny_dense)
    small_dense = list(tables.small_dense)
    small_peel = list(tables.small_peel)
    ledger = []
    for row in rows:
        n = int(row["N"])
        current_count = tables.dense_count(n)
        current_dense = tables.dense_seed(n, current_count)
        current_peel = tables.peel_seed(n)
        if (
            int(row["CurrentDenseCount"]) != current_count
            or int(row["CurrentDenseSeed"]) != current_dense
            or int(row["CurrentPeelSeed"]) != current_peel
        ):
            raise ApplyError(f"small current-value mismatch at N={n}")
        candidate_count = int(row["DenseCount"])
        if n < 65:
            tiny_counts[n] = candidate_count
            tiny_dense[n] = int(row["DenseSeed"])
        else:
            if candidate_count != current_count:
                raise ApplyError(
                    f"small N={n} selected count {candidate_count}, but runtime graph is {current_count}"
                )
            small_dense[n - 65] = int(row["DenseSeed"])
        small_peel[n] = int(row["PeelSeed"])
        candidate = f"{candidate_count}/{row['DenseSeed']}/{row['PeelSeed']}"
        current = f"{current_count}/{current_dense}/{current_peel}"
        ledger.append({
            "Key": str(n), "Current": current, "Candidate": candidate,
            "Decision": "replace" if candidate != current else "retain-same",
            "Evidence": f"dense_fail={row['DenseFailures']};peel_score={row['PeelScore']}",
        })
    rendered = replace_numeric_array(source, "kTinyDenseCounts", tiny_counts)
    rendered = replace_numeric_array(rendered, "kTinyDenseSeeds", tiny_dense)
    rendered = replace_numeric_array(rendered, "kSmallDenseSeeds", small_dense)
    rendered = replace_numeric_array(rendered, "kSmallPeelSeeds", small_peel)
    return rendered, ledger


def apply_dense_seeds(source, rows):
    tables = SourceTables(source)
    _assert_single_hash(rows, "SourceTablesHash", tables.dense_seed_hash())
    expected = {}
    for n in range(2048, 64001):
        expected[tables.dense_count(n) // 4] = n
    seeds = list(tables.dense_seeds)
    ledger = []
    for row in rows:
        index = int(row["DenseIndex"])
        used = int(row["Used"])
        expected_n = expected.get(index, 0)
        if used != int(expected_n != 0) or int(row["N"]) != expected_n:
            raise ApplyError(
                f"dense index {index} usage/N does not match candidate graph"
            )
        if int(row["CurrentDenseSeed"]) != seeds[index]:
            raise ApplyError(f"dense index {index} current seed mismatch")
        candidate = int(row["DenseSeed"])
        if used:
            seeds[index] = candidate
        elif candidate != seeds[index]:
            raise ApplyError(f"unused dense index {index} changed unexpectedly")
        ledger.append({
            "Key": str(index), "Current": str(tables.dense_seeds[index]),
            "Candidate": str(candidate),
            "Decision": "replace" if used and candidate != tables.dense_seeds[index] else "retain",
            "Evidence": f"used={used};N={expected_n};failures={row['Failures']}",
        })
    return replace_numeric_array(source, "kDenseSeeds", seeds), ledger


def apply_peel(source, rows, allow_retained=False):
    tables = SourceTables(source)
    _assert_single_hash(rows, "BaseTablesHash", tables.peel_hash())
    seeds = list(tables.peel_seeds)
    ledger = []
    for row in rows:
        subdivision = int(row["Subdivision"])
        current = tables.peel_seeds[subdivision]
        if int(row["CurrentPeelSeed"]) != current:
            raise ApplyError(f"peel subdivision {subdivision} current seed mismatch")
        status = row["Status"]
        candidate = int(row["PeelSeed"])
        if status == "searched":
            seeds[subdivision] = candidate
            decision = "replace" if candidate != current else "retain-same"
        elif allow_retained:
            if candidate != current:
                raise ApplyError(
                    f"retained peel subdivision {subdivision} changed seed"
                )
            decision = "retain-unsearched"
        else:
            raise ApplyError(f"unsearched peel subdivision {subdivision}")
        ledger.append({
            "Key": str(subdivision), "Current": str(current),
            "Candidate": str(candidate), "Decision": decision,
            "Evidence": f"status={status};score={row['Score']}",
        })
    return replace_numeric_array(source, "kPeelSeeds", seeds), ledger


APPLIERS = {
    "dense-count": apply_dense_count,
    "small": apply_small,
    "dense-seeds": apply_dense_seeds,
    "peel": apply_peel,
}


def write_atomic(path, content):
    try:
        atomic_publish.publish_files_no_replace([(path, content)])
    except atomic_publish.AtomicPublishError as exc:
        raise ApplyError(str(exc)) from exc


def render_ledger(rows):
    output = []
    stream = _ListWriter(output)
    writer = csv.DictWriter(
        stream,
        fieldnames=["Key", "Current", "Candidate", "Decision", "Evidence"],
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(rows)
    return "".join(output)


def write_ledger(path, rows):
    write_atomic(path, render_ledger(rows))


def publish_outputs(entries):
    try:
        atomic_publish.publish_files_no_replace(entries)
    except atomic_publish.AtomicPublishError as exc:
        raise ApplyError(str(exc)) from exc


class _ListWriter:
    def __init__(self, target):
        self.target = target

    def write(self, value):
        self.target.append(value)
        return len(value)


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] == "fingerprint":
        fingerprint_parser = argparse.ArgumentParser()
        fingerprint_parser.add_argument("fingerprint")
        fingerprint_parser.add_argument("--source", required=True)
        fingerprint_args = fingerprint_parser.parse_args(argv)
        try:
            source = Path(fingerprint_args.source).read_text(encoding="ascii")
            tables = SourceTables(source)
        except (OSError, ApplyError) as exc:
            print(exc, file=sys.stderr)
            return 1
        print(f"small_hash={tables.small_hash()}")
        print(f"dense_count_hash={tables.dense_count_hash()}")
        print(f"dense_graph_hash={tables.dense_graph_hash()}")
        print(f"dense_seed_hash={tables.dense_seed_hash()}")
        print(f"peel_hash={tables.peel_hash()}")
        return 0

    parser = argparse.ArgumentParser()
    parser.add_argument("phase", choices=sorted(APPLIERS))
    parser.add_argument("--source", required=True)
    parser.add_argument("--results", required=True, nargs="+")
    parser.add_argument("--output", required=True)
    parser.add_argument("--ledger", required=True)
    parser.add_argument("--table-generator")
    parser.add_argument("--table-generator-output")
    parser.add_argument("--allow-retained", action="store_true")
    args = parser.parse_args(argv)
    try:
        source = Path(args.source).read_text(encoding="ascii")
        rows = _read_validated(
            args.phase, args.results, allow_retained=args.allow_retained
        )
        if args.phase == "peel":
            rendered, ledger = apply_peel(
                source, rows, allow_retained=args.allow_retained
            )
        else:
            rendered, ledger = APPLIERS[args.phase](source, rows)
        table_generator_rendered = None
        if args.phase == "dense-count":
            if not args.table_generator or not args.table_generator_output:
                raise ApplyError(
                    "dense-count requires --table-generator and "
                    "--table-generator-output"
                )
            table_generator = Path(args.table_generator).read_text(encoding="ascii")
            candidate_points, _, _ = parse_dense_points(rendered)
            table_generator_rendered = replace_dense_points(
                table_generator, candidate_points
            )
            if parse_dense_points(table_generator_rendered)[0] != candidate_points:
                raise ApplyError("TableGenerator dense points differ after rendering")
        outputs = [(args.output, rendered)]
        if table_generator_rendered is not None:
            outputs.append((args.table_generator_output, table_generator_rendered))
        outputs.append((args.ledger, render_ledger(ledger)))
        publish_outputs(outputs)
    except (OSError, ApplyError, results.ValidationError) as exc:
        print(exc, file=sys.stderr)
        return 1
    print(
        f"rendered {args.phase}: {len(ledger)} decisions -> {args.output}",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
