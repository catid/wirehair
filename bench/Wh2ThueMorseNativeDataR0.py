#!/usr/bin/env python3
"""Authenticated, deterministic .63 recorded evidence -> native build data.

Imports are inert. Generation only reconstructs the already-fixed lookup
tables and copies recorded IDs/ranks; it never derives loss roots, selects a
parameter, evaluates packet rows, scores ranks, or runs a native executable.
"""
import argparse
import hashlib
import importlib.util
import os
from pathlib import Path
import stat
import struct
import sys


HERE = Path(__file__).resolve().parent
PROTOCOL = "wirehair.wh2.thue-morse-native-r0"
RECOVERY = Path("/var/tmp/wh2-thue-morse-recovery-r0")
MANIFEST_SHA = "5e4416e6701e9d725d431629ebfec6c918f2f70d769be51a73ee4dbc3e3b623f"
RAW_SHA = "475af51aa7c312e31298d3744f53c12b2283c642a2e6b2875f298bd823a42486"
PACKED_SHA = "27b105e1449bec190bd3c83f07feefa639cd32bc356baebfb03828ea7cbccb6d"
VISIT_SHA = "b9187d801c8c92e31672dc53fe0b9ccc1c186858196c1f9c00a2ad37028a7d61"
MEMBERS = ("CLAIM.json", "raw.json", "stderr.txt", "summary.json")
WIDTHS = (2, 64, 1280)
SCHEDULES = ("iid", "burst", "adversarial", "repair-only")
BYTE_CAP = 4 * 1024 * 1024


def sibling(name, filename):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


H = sibling("_thue_native_history", "Wh2ThueMorseRecoveryHistoryR0.py")
C = H.C
canonical, require = C.canonical, H.require


def digest(data):
    return hashlib.sha256(data).hexdigest()


def exact(actual, expected, label):
    require(canonical(actual) == canonical(expected), label)


def integer(value, low, high, label):
    require(type(value) is int and low <= value <= high, label)
    return value


def recovery_bundle(deadline=None):
    raw = C.read_regular(RECOVERY / "COMPLETE.json", BYTE_CAP, deadline=deadline)
    require(digest(raw) == MANIFEST_SHA, "recovery manifest pin")
    manifest = C.strict_json(raw)
    require(raw == canonical(manifest), "recovery canonical manifest")
    require(type(manifest) is dict and set(manifest) == {"protocol", "outcome", "files"} and
            manifest["protocol"] == "wirehair.wh2.thue-morse-recovery-r0" and
            manifest["outcome"] == "PASS" and set(manifest["files"]) == set(MEMBERS),
            "recovery manifest schema/verdict")
    seen = set()
    with os.scandir(RECOVERY) as entries:
        for entry in entries:
            C.time_left(deadline)
            require(entry.name in set(MEMBERS) | {"COMPLETE.json"}, "unexpected recovery member")
            seen.add(entry.name)
    require(seen == set(MEMBERS) | {"COMPLETE.json"}, "incomplete recovery directory")
    content, pins, total = {}, {}, len(raw)
    for name in MEMBERS:
        C.time_left(deadline)
        data = C.read_regular(RECOVERY / name, BYTE_CAP, deadline=deadline)
        total += len(data)
        require(total <= 8 * 1024 * 1024, "recovery aggregate input cap")
        pins[name] = dict(bytes=len(data), sha256=digest(data))
        exact(manifest["files"][name], pins[name], "recovery member pin: " + name)
        content[name] = data
    require(not content["stderr.txt"] and digest(content["raw.json"]) == RAW_SHA,
            "recovery raw/stderr pin")
    report = C.strict_json(content["raw.json"])
    require(content["raw.json"] == canonical(report) + b"\n", "recovery raw canonical encoding")
    C.time_left(deadline)
    return report, dict(manifest_sha256=MANIFEST_SHA, manifest_bytes=len(raw), files=pins)


def authenticated_inputs(deadline=None):
    """Existing evidence only; safe for the controller/reducer to authenticate."""
    report, _ = recovery_bundle(deadline=deadline)
    inputs = H.load_inputs(deadline=deadline)
    exact(report["inputs"], inputs, "historical inputs differ from sealed recovery")
    require(report["inputs_sha256"] == digest(canonical(inputs)), "recovery input seal")
    require(report["protocol"] == "wirehair.wh2.thue-morse-recovery-r0" and report["outcome"] == "PASS",
            "recovery result protocol/verdict")
    require(report["evidence"]["lookup"]["tables_sha256"] == PACKED_SHA and
            report["evidence"]["coefficient_visit_sha256"] == VISIT_SHA,
            "fixed lookup/visit pins")
    C.time_left(deadline)
    return report


def checked_ids(ids, count):
    require(type(ids) is list and len(ids) == count, "packet ID count")
    for packet_id in ids:
        integer(packet_id, 0, (1 << 32) - 1, "packet ID type/range")
    require(len(set(ids)) == count, "duplicate packet ID")
    return list(ids)


def checked_ranks(ranks):
    require(type(ranks) is list and len(ranks) == 5, "prefix rank count")
    for rank in ranks:
        integer(rank, 1, 6, "rank type/range")
    require(all(a <= b <= a + 1 for a, b in zip(ranks, ranks[1:])), "nonmonotone ranks")
    return list(ranks)


def extract_data(report):
    """Pure typed projection of already-recorded evidence, usable with neutral fixtures."""
    inputs = report["inputs"]
    for key, count in (("fresh", 6144), ("hard", 18), ("history", 82), ("fixtures", 33)):
        require(type(report[key]) is list and len(report[key]) == count, "recorded roster: " + key)
    require(len(inputs["prefixes"]) == 82 and len(inputs["fixtures"]) == 33, "input prefix roster")
    pair = []
    for key in ("A0", "A1"):
        matrix = inputs["pair"][key]
        require(type(matrix) in (tuple, list) and len(matrix) == 6, "matrix rows")
        copied = []
        for row in matrix:
            require(type(row) in (tuple, list) and len(row) == 6, "matrix columns")
            copied.append([integer(x, 0, 255, "matrix byte") for x in row])
        pair.append(copied)
    traces = []
    for ordinal, row in enumerate(report["fresh"] + report["hard"]):
        if ordinal < 6144:
            index, B, schedule = ordinal // 12, WIDTHS[(ordinal % 12) // 4], SCHEDULES[ordinal % 4]
            exact([row["root_index"], row["B"], row["schedule"]], [index, B, schedule], "fresh ordering")
        else:
            offset = ordinal - 6144
            B = 2
            exact([row["phase"], row["root_index"], row["B"], row["schedule"]],
                  ["training" if offset < 9 else "validation", (offset % 9) // 3, B,
                   SCHEDULES[1 + offset % 3]], "hard ordering")
        ids = checked_ids(row["ids"], 10)
        integer(row["attempted_candidates"], 10, 68096, "candidate count")
        require(digest(b"".join(struct.pack("<I", x) for x in ids)) == row["trace_sha256"], "ID hash")
        traces.append(dict(B=B, ids=ids, ranks=checked_ranks(row["ranks"])))
    history, fixtures = [], []
    for index, (prefix, result) in enumerate(zip(inputs["prefixes"], report["history"])):
        exact(prefix["K"], 6, "historical K")
        overhead = integer(prefix["overhead"], 0, 1, "historical overhead")
        rank = integer(result["rank"], 1, 6, "historical rank")
        exact(result, dict(index=index, rank=rank, repaired=rank == 6), "historical result")
        widths = prefix["original_widths"]
        require(type(widths) is list and widths and
                all(type(B) is int and B in WIDTHS for B in widths) and widths == sorted(set(widths)),
                "historical original widths")
        history.append(dict(ids=checked_ids(prefix["ids"], 6 + overhead), rank=rank, widths=list(widths)))
    for index, (fixture, result) in enumerate(zip(inputs["fixtures"], report["fixtures"])):
        rank = integer(fixture["expected_rank"], 5, 6, "fixture rank")
        exact(result, dict(index=index, rank=rank), "fixture result")
        count = integer(len(fixture["ids"]), 6, 8, "fixture length")
        fixtures.append(dict(ids=checked_ids(fixture["ids"], count), rank=rank, widths=list(WIDTHS)))
    history_cases = [[i, B, B] for i, row in enumerate(history) for B in row["widths"]]
    fixture_cases = [[i, B, B] for i in range(33) for B in WIDTHS]
    partial_cases = [[6144 + i, B, tail] for i in range(18) for B in WIDTHS for tail in sorted({1, B - 1})]
    require(len(history_cases) == 105 and len(fixture_cases) == 99 and len(partial_cases) == 90,
            "payload projection counts")
    require(sum(len(row["ids"]) for row in history) == 496 and
            sum(len(row["ids"]) for row in fixtures) == 231, "original coefficient visit lengths")
    return dict(protocol=PROTOCOL, pair=pair, traces=traces, history=history, fixtures=fixtures,
                history_cases=history_cases, fixture_cases=fixture_cases, partial_cases=partial_cases,
                recovery_raw_sha256=digest(canonical(report) + b"\n"), inputs_sha256=report["inputs_sha256"],
                lookup_sha256=report["evidence"]["lookup"]["tables_sha256"],
                coefficient_visit_sha256=report["evidence"]["coefficient_visit_sha256"])


def build_lookup(pair):
    """Build fixed tables only: no .row(), reference_row(), rank, or selector."""
    math = sibling("_thue_native_fixed_tables", "Wh2ThueMorseR0.py")
    mapper = math.RowMapper(tuple(tuple(tuple(row) for row in matrix) for matrix in pair))
    packed = b"".join(mapper.low + mapper.middle10 + mapper.middle17 + (mapper.high,))
    return packed, mapper.evidence()


def render_header(data, packed):
    require(type(packed) is bytes and len(packed) == 39936 and digest(packed) == data["lookup_sha256"],
            "render lookup bytes/hash")
    lines = ["// Generated from sealed .62/.63 evidence; include in the worker TU only.",
             "#ifndef WH2_THUE_MORSE_NATIVE_DATA_R0_H", "#define WH2_THUE_MORSE_NATIVE_DATA_R0_H",
             "#include <stdint.h>", "namespace wh2_thue_native_data {",
             "enum { kTraceCount=6162, kHistoryCount=82, kFixtureCount=33,",
             "       kHistoryCaseCount=105, kFixtureCaseCount=99, kPartialCaseCount=90 };",
             "struct Trace { uint16_t B; uint32_t ids[10]; uint8_t ranks[5]; };",
             "struct Prefix { uint8_t count, rank, width_mask; uint32_t ids[8]; };",
             "struct PayloadCase { uint16_t index, B, tail; };"]
    for name, value in (("kLookupSha256", data["lookup_sha256"]),
                        ("kCoefficientVisitSha256", data["coefficient_visit_sha256"]),
                        ("kDataSha256", digest(canonical(data))),
                        ("kRecoveryRawSha256", data["recovery_raw_sha256"])):
        require(type(value) is str and len(value) == 64 and all(x in "0123456789abcdef" for x in value),
                "header hash encoding")
        lines.append('extern const char {}[] = "{}";'.format(name, value))
    lines.append("alignas(64) extern const uint8_t kLookup[39936] = {")
    for offset in range(0, len(packed), 24):
        lines.append("  " + ",".join(str(x) for x in packed[offset:offset + 24]) + ",")
    lines.append("};\nextern const uint8_t kPair[2][6][6] = {")
    for matrix in data["pair"]:
        lines.append("  {" + ",".join("{" + ",".join(map(str, row)) + "}" for row in matrix) + "},")
    lines.append("};\nextern const Trace kTraces[kTraceCount] = {")
    for row in data["traces"]:
        lines.append("  {" + str(row["B"]) + ",{" + ",".join(str(x) + "u" for x in row["ids"]) +
                     "},{" + ",".join(map(str, row["ranks"])) + "}},")
    lines.append("};")
    for name, rows in (("kHistory", data["history"]), ("kFixtures", data["fixtures"])):
        lines.append("extern const Prefix {}[{}] = {{".format(name, len(rows)))
        for row in rows:
            width_mask = sum(1 << WIDTHS.index(B) for B in row["widths"])
            ids = row["ids"] + [0] * (8 - len(row["ids"]))
            lines.append("  {" + ",".join(map(str, (len(row["ids"]), row["rank"], width_mask))) +
                         ",{" + ",".join(str(x) + "u" for x in ids) + "}},")
        lines.append("};")
    for name, key in (("kHistoryCases", "history_cases"), ("kFixtureCases", "fixture_cases"),
                      ("kPartialCases", "partial_cases")):
        lines.append("extern const PayloadCase {}[{}] = {{".format(name, len(data[key])))
        lines.extend("  {" + ",".join(map(str, row)) + "}," for row in data[key])
        lines.append("};")
    lines.extend(("} // namespace wh2_thue_native_data", "#endif", ""))
    result = "\n".join(lines).encode("ascii")
    require(len(result) <= BYTE_CAP, "generated header cap")
    return result


def write_create_only(name, data, directory):
    require(type(name) is str and Path(name).name == name and name not in ("", ".", ".."),
            "artifact basename")
    descriptor = os.open(name, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                         0o400, dir_fd=directory)
    try:
        with os.fdopen(descriptor, "wb", closefd=False) as stream:
            stream.write(data)
            stream.flush()
            os.fsync(descriptor)
        info = os.fstat(descriptor)
        require(stat.S_ISREG(info.st_mode) and stat.S_IMODE(info.st_mode) == 0o400 and
                info.st_nlink == 1 and info.st_size == len(data), "created artifact identity")
    finally:
        os.close(descriptor)


def generate(header_path, receipt_path):
    header, receipt = Path(header_path), Path(receipt_path)
    require(header.is_absolute() and receipt.is_absolute() and header.parent == receipt.parent and header != receipt,
            "generated paths must be distinct absolute siblings")
    require(header.parent.resolve(strict=True) == header.parent, "generated directory must not alias a symlink")
    directory = os.open(str(header.parent), os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        for name in (header.name, receipt.name):
            try:
                os.stat(name, dir_fd=directory, follow_symlinks=False)
            except FileNotFoundError:
                pass
            else:
                raise ValueError("generated output already exists")
        report = authenticated_inputs()
        data = extract_data(report)
        packed, lookup = build_lookup(data["pair"])
        exact(lookup, report["evidence"]["lookup"], "all seven frozen table hashes/evidence")
        output = render_header(data, packed)
        repeated, provenance = recovery_bundle()
        exact(repeated, report, "recovery input changed during generation")
        exact(H.load_inputs(), report["inputs"], "historical inputs changed during generation")
        result = dict(protocol=PROTOCOL, recovery_provenance=provenance,
                      historical_provenance=report["inputs"]["provenance"], lookup=lookup,
                      data_sha256=digest(canonical(data)), header_sha256=digest(output), header_bytes=len(output),
                      counts=dict(traces=6162, history=82, fixtures=33, history_cases=105,
                                  fixture_cases=99, partial_cases=90, coefficient_visits=62347))
        write_create_only(header.name, output, directory)
        write_create_only(receipt.name, canonical(result), directory)
        os.fsync(directory)
        return result
    finally:
        os.close(directory)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generate", metavar="HEADER", required=True)
    parser.add_argument("--receipt", required=True)
    args = parser.parse_args()
    try:
        result = generate(args.generate, args.receipt)
        print(canonical(result).decode("ascii"))
    except (ValueError, OSError, TypeError, KeyError, IndexError) as error:
        print("native data: " + str(error)[:1000], file=sys.stderr)
        sys.exit(1)
