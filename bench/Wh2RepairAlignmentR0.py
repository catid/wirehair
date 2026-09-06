#!/usr/bin/env python3
"""Frozen, one-shot public-repair placement diagnostic (Beads .64).

Build/selftest modes do not construct a codec or run a timed fixture. Replay
reduces recorded evidence only. No old timing namespace or outcome is reused.
"""
import argparse
import math
import os
from pathlib import Path
import struct
import sys
import time
import importlib.util


def sibling(name, filename):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(filename))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


D = sibling("_repair_alignment_reviewed_tools", "Wh2DirectRowObserveR0.py")
M, N, C, R = D.M, D.N, D.C, D.R
ROOT = D.ROOT
require, exact, integer = M.require, M.exact, M.integer
PROTOCOL = "wirehair.wh2.repair-alignment-r0"
TARGET = "wh2_repair_alignment_r0"
OUTPUT = Path("/var/tmp/wh2-repair-alignment-r0")
HISTORICAL_RECEIPT = Path("/var/tmp/wh2-direct-row-observe-r0/CLAIM.json")
HISTORICAL_SHA = "b171445f0093f15db38c85a895c8356020547635845a9f52c32fa159989789bc"
SOURCES = tuple("bench/" + name for name in (
    "Wh2RepairAlignmentR0Bridge.h", "Wh2RepairAlignmentR0Bridge.cpp",
    "Wh2RepairAlignmentR0.cpp", "Wh2RepairAlignmentR0.py", "test_Wh2RepairAlignmentR0.py"))
RAW_CAP, ERR_CAP, TOTAL_CAP = 2 * 1024**2, 65536, 8 * 1024**2
COUNTERS = ("minflt", "majflt", "nvcsw", "nivcsw")
METRICS = ("encode_wall_ns", "thread_span_ns", "outer_wall_ns", "instrumentation_envelope_ns")
U64 = (1 << 64) - 1
PRELUDE_SEED, PRELUDE_FINAL = 0x9e3779b97f4a7c15, 0x43935dad1647741b
T11 = 2.200985160082949


def endpoints(phase, comparison):
    if comparison < 4:
        return comparison, comparison
    return (((0, 2), (1, 3), (0, 1), (2, 3)) if phase == 0 else
            ((0, 1), (2, 3)))[comparison - 4]


def roster():
    """All callbacks in chronological order, including both warm positions."""
    index = 0
    sides = (0, 1, 0, 1, 1, 0, 1, 0, 0, 1)
    for rep in range(12):
        for step in range(2):
            epoch, order = 2 * rep + step, (rep + step) % 2
            for phase_step in range(2):
                phase = (rep + step + phase_step) % 2
                count = 8 if phase == 0 else 6
                for class_step in range(count):
                    comparison = (epoch + class_step) % count
                    pair = endpoints(phase, comparison)
                    for position, side in enumerate(sides):
                        yield (index, rep, order, phase, comparison, position, pair[side ^ order])
                        index += 1


def fnv(data, value=1469598103934665603):
    for byte in data:
        value = ((value * 1099511628211) ^ byte) & U64
    return value


def expected_checksum(packets, deadline=None):
    # Every actual callback checks all packet bytes and folds the little-endian
    # FNV of the six concatenated packets into the retained running checksum.
    C.time_left(deadline)
    require(type(packets) is bytes and len(packets) == 7680, "checksum packet shape")
    block = struct.pack("<Q", fnv(packets))
    value = 1469598103934665603
    for _ in range(3360):
        C.time_left(deadline)
        value = fnv(block, value)
    return value


def hex_bytes(value, size, label):
    require(type(value) is str and len(value) == size * 2, label + " size")
    try:
        data = bytes.fromhex(value)
    except ValueError as error:
        raise ValueError(label + " encoding") from error
    require(len(data) == size and data.hex() == value, label + " canonical encoding")
    return data


def confidence(values):
    require(len(values) == 12 and all(math.isfinite(x) for x in values), "12 independent replicate contrasts")
    mean = math.fsum(values) / 12
    variance = math.fsum((x - mean)**2 for x in values) / 11
    radius = T11 * math.sqrt(variance / 12)
    return {"replicates": 12, "mean_log": mean, "lower95_log": mean - radius,
            "upper95_log": mean + radius, "ratio": math.exp(mean),
            "lower95": math.exp(mean - radius), "upper95": math.exp(mean + radius)}


def clock_record(values, previous, first, elapsed):
    require(type(values) is list and len(values) == 8, "clock envelope shape")
    m0, c0, m1, m2, c1, m3 = [integer(x, 0, (1 << 63) - 1, "timestamp") for x in values[:6]]
    require(m0 <= m1 < m2 <= m3 and c0 <= c1, "within-clock ordering")
    require(previous[0] <= m0 and previous[1] <= c0, "cross-record clock ordering")
    usage = []
    for item in values[6:]:
        require(type(item) is list and len(item) == 4, "usage shape")
        usage.append([integer(x, 0, (1 << 63) - 1, "usage counter") for x in item])
    require(all(previous[2][i] <= usage[0][i] <= usage[1][i] for i in range(4)), "usage monotonicity")
    require(m3 - first[0] <= elapsed and c1 - first[1] <= elapsed, "observation exceeds worker lifetime")
    derived = {"encode_wall_ns": m2 - m1, "thread_span_ns": c1 - c0,
               "outer_wall_ns": m3 - m0, "instrumentation_envelope_ns": m3 - m0 - (m2 - m1),
               "counter_deltas": {name: usage[1][i] - usage[0][i] for i, name in enumerate(COUNTERS)}}
    return derived, (m3, c1, usage[1])


def verify_metadata(metadata, deadline=None):
    C.time_left(deadline)
    R.keys(metadata, ("source_hex", "source_sha256", "intermediate_hex", "intermediate_sha256",
                      "expected_packets_hex", "expected_packets_sha256", "columns",
                      "expected_operations", "handles"), "metadata")
    source = hex_bytes(metadata["source_hex"], 7680, "source")
    exact(source.hex(), bytes((37 * i + i // 11) & 255 for i in range(7680)).hex(), "frozen source law")
    master = hex_bytes(metadata["intermediate_hex"], 46080, "intermediate")
    for name, data in (("source", source), ("intermediate", master)):
        exact(metadata[name + "_sha256"], C.sha(data), name + " hash")
    require(type(metadata["columns"]) is list and len(metadata["columns"]) == 6, "six packet rows")
    require(type(metadata["expected_operations"]) is list and len(metadata["expected_operations"]) == 6,
            "six operation counts")
    packets = bytearray()
    for index, columns in enumerate(metadata["columns"]):
        C.time_left(deadline)
        require(type(columns) is list and 4 <= len(columns) <= 9, "packet term count")
        for column in columns:
            integer(column, 0, 35, "packet column")
        require(len(set(columns)) == len(columns) and sum(column >= 6 for column in columns) == 3,
                "distinct source plus three mix columns")
        exact(metadata["expected_operations"][index], len(columns), "operation accounting")
        result = bytearray(1280)
        for column in columns:
            block = master[column * 1280:(column + 1) * 1280]
            for offset, byte in enumerate(block):
                result[offset] ^= byte
        packets.extend(result)
    packets = bytes(packets)
    exact(metadata["expected_packets_hex"], packets.hex(), "scalar XOR packet oracle")
    exact(metadata["expected_packets_sha256"], C.sha(packets), "packet hash")
    require(type(metadata["handles"]) is list and len(metadata["handles"]) == 2, "two live encoders")
    exact(metadata["handles"][0], metadata["handles"][1], "encoder equation/data mismatch")
    descriptor = metadata["handles"][0]
    R.keys(descriptor, ("profile_hex", "profile_sha256", "source_sha256", "intermediate_sha256",
                        "message_bytes", "block_bytes", "source_count", "precode_count", "intermediate_bytes",
                        "source_policy", "seed_attempt", "params", "config", "runtime"), "encoder descriptor")
    profile = hex_bytes(descriptor["profile_hex"], 32, "serialized profile")
    attempt = integer(descriptor["seed_attempt"], 0, 255, "seed attempt")
    expected_profile = b"WHV2" + struct.pack("<HHQQI", 1, 32, 0x4b295bbb47f4f9c9, 7680, 1280) + bytes([attempt, 0, 0, 0])
    exact(profile.hex(), expected_profile.hex(), "certified profile binding")
    for key, value in (("profile_sha256", C.sha(profile)), ("source_sha256", C.sha(source)),
                       ("intermediate_sha256", C.sha(master)), ("message_bytes", 7680), ("block_bytes", 1280),
                       ("source_count", 6), ("precode_count", 30), ("intermediate_bytes", 46080), ("source_policy", 2)):
        exact(descriptor[key], value, key)
    params = descriptor["params"]
    R.keys(params, ("block_count", "staircase", "dense_rows", "heavy_rows", "source_hits",
                    "dense_identity_corner", "heavy_family", "dense_anchors", "seed"), "precode parameters")
    expected = dict(block_count=6, staircase=6, dense_rows=12, heavy_rows=12, source_hits=2,
                    dense_identity_corner=False, heavy_family=0, dense_anchors=0)
    exact({key: value for key, value in params.items() if key != "seed"}, expected, "certified dimensions")
    integer(params["seed"], 0, U64, "precode seed")
    R.keys(descriptor["config"], ("peel_seed", "mix_count"), "packet config")
    integer(descriptor["config"]["peel_seed"], 0, (1 << 32) - 1, "peel seed")
    exact(descriptor["config"]["mix_count"], 3, "certified mix count")
    exact(descriptor["runtime"], dict(source_prime=7, precode_prime=31), "reconstructed row runtime")
    return packets


def verify_addresses(addresses):
    R.keys(addresses, ("source", "master", "handles", "intermediates", "carriers", "outputs", "scratch",
                       "public_function", "shadow_function", "runner_function", "prelude_function"), "addresses")
    spans = []
    for name, count, size, aligned in (("source", 1, 7680, False), ("master", 1, 46080, False),
            ("handles", 2, 1, False), ("intermediates", 2, 46080, False),
            ("carriers", 2, 65536, True), ("outputs", 2, 12288, True), ("scratch", 1, 492544, False)):
        items = addresses[name] if count == 2 else [addresses[name]]
        require(type(items) is list and len(items) == count, "address count: " + name)
        for index, item in enumerate(items):
            R.keys(item, ("address", "mod32", "mod64", "mod4096"), "address descriptor")
            address = integer(item["address"], 1, U64 - size, "address range")
            for modulus in (32, 64, 4096):
                exact(item["mod" + str(modulus)], address % modulus, "address residue")
            require(not aligned or address % 4096 == 0, "page-aligned carrier/output")
            spans.append((address, address + size, name + str(index)))
    for i, left in enumerate(spans):
        for right in spans[i + 1:]:
            require(left[1] <= right[0] or right[1] <= left[0], "overlapping live storage: " + left[2] + "/" + right[2])
    functions = []
    for name in ("public_function", "shadow_function", "runner_function", "prelude_function"):
        functions.append(integer(addresses[name], 1, U64, "function address"))
    require(len(set(functions)) == 4, "distinct public/shadow/runner/prelude entries")


def reduce_raw(raw, deadline=None):
    C.time_left(deadline)
    R.keys(raw, ("protocol", "schema", "outcome", "failure", "target_cpu", "identity_before", "identity_after",
                 "runtime_before", "runtime_after", "worker_start_ns", "monotonic_resolution_ns", "thread_resolution_ns",
                 "metadata", "handles_after", "addresses", "addresses_after", "preflight", "prelude", "callbacks", "encode_calls",
                 "checked_packets", "checksum", "sum_encode_wall_ns", "elapsed_ns", "records"), "raw")
    for key, value in (("protocol", PROTOCOL), ("schema", PROTOCOL + ".raw.v1"), ("outcome", "COMPLETE"),
                       ("failure", None), ("target_cpu", 50), ("callbacks", 3360), ("encode_calls", 1290240),
                       ("checked_packets", 20160)):
        exact(raw[key], value, key)
    elapsed = integer(raw["elapsed_ns"], 1, 10**10 - 1, "worker elapsed")
    start = integer(raw["worker_start_ns"], 0, (1 << 63) - 1, "worker start")
    for key in ("monotonic_resolution_ns", "thread_resolution_ns"):
        integer(raw[key], 1, 10**9, key)
    D.identity_and_runtime(raw)
    packets = verify_metadata(raw["metadata"], deadline)
    exact(raw["metadata"]["handles"], raw["handles_after"], "encoder state changed")
    exact(raw["checksum"], expected_checksum(packets, deadline), "consumed packet checksum")
    verify_addresses(raw["addresses"])
    exact(raw["addresses"], raw["addresses_after"], "live address drift")
    exact(raw["preflight"], dict(public=12, shadow=24, original_view=6, scalar=6), "preflight call roster")
    prelude = raw["prelude"]
    R.keys(prelude, ("iterations", "seed", "final_state", "m0", "c0", "m1", "m2", "c1", "m3", "ru0", "ru1"), "prelude")
    exact([prelude[key] for key in ("iterations", "seed", "final_state")],
          [1 << 20, PRELUDE_SEED, PRELUDE_FINAL], "fixed prelude identity")
    envelope = [prelude[key] for key in ("m0", "c0", "m1", "m2", "c1", "m3", "ru0", "ru1")]
    cpu_start = integer(prelude["c0"], 0, (1 << 63) - 1, "prelude CPU start")
    conditioning, previous = clock_record(envelope, (start, 0, [0] * 4), (start, cpu_start), elapsed)
    conditioning["prelude_wall_ns"] = conditioning.pop("encode_wall_ns")
    conditioning["iterations"] = 1 << 20
    require(type(raw["records"]) is list and len(raw["records"]) == 3360, "complete callback roster")
    derived, panel_samples, cells = [], {}, {}
    total = 0
    for record, coordinate in zip(raw["records"], roster()):
        C.time_left(deadline)
        require(type(record) is list and len(record) == 19, "callback shape")
        exact(record[:7], list(coordinate), "callback chronology/roster")
        metrics, previous = clock_record(record[7:15], previous, (start, cpu_start), elapsed)
        index, rep, order, phase, comparison, position, cell = coordinate
        exact(record[15], 384, "calls per callback")
        exact(record[16], [0] * 6, "all packet statuses")
        exact(record[17], [1280] * 6, "all packet lengths")
        expected_ops = [None] * 6 if phase == 0 else [64 * n for n in raw["metadata"]["expected_operations"]]
        exact(record[18], expected_ops, "shadow operations/public unavailable")
        metrics.update(index=index, rep=rep, order=order, phase=phase, comparison=comparison,
                       position=position, cell=cell, warmup=position < 2)
        derived.append(metrics)
        total += metrics["encode_wall_ns"]
        panel_samples.setdefault((phase, comparison, order, rep), []).append(metrics)
        cells.setdefault((phase, cell), []).append(metrics)
    integer(raw["sum_encode_wall_ns"], 1, 200000000, "cumulative encode cap")
    exact(raw["sum_encode_wall_ns"], total, "cumulative encode sum")
    contrasts = {}
    panels = []
    sides = (0, 1, 1, 0, 1, 0, 0, 1)
    for key, samples in panel_samples.items():
        phase, comparison, order, rep = key
        require(len(samples) == 10, "whole panel")
        measured = samples[2:]
        log_sides = [math.fsum(math.log(row["encode_wall_ns"]) for row, side in zip(measured, sides)
                               if side ^ order == label) / 4 for label in (0, 1)]
        contrast = log_sides[1] - log_sides[0]
        contrasts.setdefault((phase, comparison, order), {})[rep] = contrast
        panels.append({"phase": phase, "comparison": comparison, "order": order, "rep": rep,
                       "first_record": samples[0]["index"], "right_over_left_log": contrast})
    stats, failures, alignment = [], [], []
    limit = math.log1p(.02)
    for key in sorted(contrasts):
        phase, comparison, order = key
        exact(sorted(contrasts[key]), list(range(12)), "independent replicate closure")
        estimate = confidence([contrasts[key][rep] for rep in range(12)])
        item = {"phase": phase, "comparison": comparison, "order": order, "estimate": estimate}
        if comparison < 4:
            item["AA_pass"] = estimate["lower95_log"] > -limit and estimate["upper95_log"] < limit
            if not item["AA_pass"]:
                failures.append(list(key))
        elif phase == 1:
            item["placement_supported"] = estimate["lower95_log"] > 0
            alignment.append(item["placement_supported"])
        stats.append(item)
    require(len(stats) == 28 and len(alignment) == 4, "complete contrast classes")
    interaction = []
    for order in (0, 1):
        values = [contrasts[(0, 5, order)][rep] - contrasts[(0, 4, order)][rep] for rep in range(12)]
        interaction.append({"order": order, "definition": "log(H1/H0 at O1)-log(H1/H0 at O0)",
                            "estimate": confidence(values)})
    absolute = []
    for key in sorted(cells):
        samples = cells[key]
        # The prospective early/late split is the first/last replicate, never
        # selected using observed durations. Warm callbacks remain included.
        groups = {"all": samples, "early_rep0": [row for row in samples if row["rep"] == 0],
                  "late_rep11": [row for row in samples if row["rep"] == 11]}
        absolute.append({"phase": key[0], "cell": key[1], "summaries": {
            name: {"callbacks": len(rows), "metrics": {metric: D.summarize([row[metric] for row in rows])
                                                       for metric in METRICS}}
            for name, rows in groups.items()}})
    C.time_left(deadline)
    return {"protocol": PROTOCOL, "outcome": "CONTROL_FAIL" if failures else
            "PLACEMENT_SUPPORTED" if all(alignment) else "PLACEMENT_INCONCLUSIVE",
            "failed_controls": failures, "statistics": stats, "panels": panels, "records": derived,
            "handle_output_interaction": interaction, "absolute_cell_summaries": absolute,
            "conditioning": conditioning, "startup_to_first_panel_ns": raw["records"][0][7] - start,
            "sum_encode_wall_ns": total, "all_samples_retained": True, "WH1_compared": False,
            "production_promotion_claimed": False, "speed_claimed": False,
            "scope": "fixed_conditioned_K6_B1280_sequential_repair_placement",
            "interpretation": "Global AA failure overrides every treatment. Placement support is not unique "
                              "split-load causality, steady-state proof, historical rescue, or end-to-end speed evidence."}


def build_commands(build):
    flags = [N.COMPILER, "-DWIREHAIR_STATIC=1", "-I" + str(ROOT), "-I" + str(ROOT / "include"),
             "-I" + str(M.BASE), "-O3", "-DNDEBUG", "-std=gnu++11", "-Wall", "-Wextra",
             "-Wpedantic", "-Werror", "-march=native"]
    names = ("Wh2RepairAlignmentR0Bridge", "Wh2RepairAlignmentR0")
    commands = []
    for index, name in enumerate(names):
        own_flags = flags + (["-DWIREHAIR_BUILDING=1", "-fPIC"] if index == 0 else [])
        commands.append(own_flags + ["-MD", "-MF", str(build / (name + ".d")), "-c",
                                    str(ROOT / "bench" / (name + ".cpp")), "-o", str(build / (name + ".o"))])
    commands.append([N.COMPILER, "-O3", "-DNDEBUG"] + [str(build / (name + ".o")) for name in names] +
                    [str(path) for path in M.HELPERS] + [str(M.BASE / "libwirehair.a"), "-lm",
                    "-Wl,-Map," + str(build / "link.map"), "-o", str(build / TARGET)])
    return commands


def historical_receipt(deadline=None):
    # Authenticate .66's imported Python source as well as its transitive .65
    # and original native build closure. Historical commits/interpreters are
    # provenance, not required to equal this new launch's source commit.
    data = C.read_regular(HISTORICAL_RECEIPT, RAW_CAP, deadline=deadline)
    exact(C.sha(data), HISTORICAL_SHA, "historical .66 receipt identity")
    old = C.strict_json(data)
    D.historical_receipt(deadline)
    exact(sorted(old["sources"]), sorted(D.SOURCES), "historical .66 source roster")
    for name, digest in old["sources"].items():
        exact(C.sha(C.read_regular(ROOT / name, RAW_CAP, deadline=deadline)), digest, "historical .66 source pin")
    for item in old["dependencies"] + old["build_files"] + old["helpers"] + [old["interpreter"]]:
        exact(N.pin(item["path"], deadline, installed=True), item, "historical .66 dependency pin")
    return old


def build_inputs(dependencies):
    result = set(dependencies) | {ROOT / name for name in SOURCES} | set(M.HELPERS)
    result.update((M.BASE / "libwirehair.a", Path(N.COMPILER).resolve(strict=True)))
    return result


def build_outputs(build):
    commands = build_commands(build)
    return {Path(command[-1]) for command in commands} | {build / "link.map"} | {
        Path(command[command.index("-MF") + 1]) for command in commands[:2]}


def validate_link_map(content, build):
    text = content.decode("utf-8")
    require("libwirehair.a(WirehairV2Profile.cpp.o)" not in text, "archive facade linked")
    require(str(build / "Wh2RepairAlignmentR0Bridge.o") in text and
            str(M.BASE / "libwirehair.a") in text and "wirehair_v2_encode" in text,
            "incomplete facade link map")
    require("Wh2ThueMorseNativeCodec.o" not in text and "Wh2ThueMorseMulRowsR0.o" not in text,
            "candidate object linked")


def build_only(build):
    require(build.is_absolute() and build.parent.is_dir() and not build.exists(), "fresh explicit build path")
    build.mkdir(mode=0o700)
    historical_receipt()
    commands = build_commands(build)
    dependencies = set()
    for command in commands[:2]:
        preprocess = command[:command.index("-MD")] + ["-M", "-MT", command[-1],
                                                       command[command.index("-c") + 1]]
        dependencies.update(M.dependency_names(N.checked(preprocess, time.monotonic() + 30), command[-1]))
    inputs = build_inputs(dependencies)
    before = [N.pin(path, installed=True) for path in sorted(inputs)]
    for command in commands:
        N.checked(command, time.monotonic() + 60)
    exact([N.pin(path, installed=True) for path in sorted(inputs)], before, "build input changed")
    historical_receipt()
    validate_link_map(C.read_regular(build / "link.map", 4 * 1024**2), build)
    C.write_new(build / "commands.json", C.canonical(commands))
    C.write_new(build / "build_inputs.json", C.canonical({"inputs": before,
                "outputs": [N.pin(path) for path in sorted(build_outputs(build))]}))
    print(N.checked([str(build / TARGET), "--selftest"], time.monotonic() + 30).decode("ascii"), end="")


def current_receipt(build, deadline=None):
    require(sys.version_info[:2] in ((3, 8), (3, 12)), "unvalidated interpreter")
    build = build.resolve(strict=True)
    historical_receipt(deadline)
    exact(C.strict_json(C.read_regular(build / "commands.json", 65536, deadline=deadline)),
          build_commands(build), "build commands")
    commit = C.git("rev-parse", "HEAD", deadline=deadline)
    require(commit == C.git("rev-parse", "@{upstream}", deadline=deadline), "sources not pushed")
    require(not C.git("status", "--porcelain", "--", *SOURCES, deadline=deadline), "dirty new sources")
    exact(sorted(C.git("ls-files", "--", *SOURCES, deadline=deadline).splitlines()), sorted(SOURCES), "source roster")
    sources = {}
    for name in SOURCES:
        data = C.read_regular(ROOT / name, RAW_CAP, deadline=deadline)
        require(data == C.git_bytes("show", commit + ":" + name, deadline=deadline), "source blob differs")
        sources[name] = C.sha(data)
    dependencies = set()
    for command in build_commands(build)[:2]:
        dep = Path(command[command.index("-MF") + 1])
        dependencies.update(M.dependency_names(C.read_regular(dep, 256 * 1024, deadline=deadline), command[-1]))
    require({ROOT / "gf256.h", ROOT / "codec/WirehairV2Profile.cpp",
             ROOT / "codec/WirehairV2Codec.h"}.issubset(dependencies), "incomplete compiler dependencies")
    require(not any(path.name in ("Wh2ThueMorseNativeCodec.h", "Wh2ThueMorseMulRowsR0.h",
                                  "Wh2ThueMorseNativeDataR0.h") for path in dependencies), "candidate dependency")
    outputs = build_outputs(build)
    M.validate_build_manifest(C.strict_json(C.read_regular(build / "build_inputs.json", RAW_CAP, deadline=deadline)),
                              build_inputs(dependencies), outputs, deadline)
    validate_link_map(C.read_regular(build / "link.map", 4 * 1024**2, deadline=deadline), build)
    symbols = N.checked(["/usr/bin/nm", "--defined-only", str(build / TARGET)], deadline)
    names = [line.split()[-1] for line in symbols.decode("ascii").splitlines() if line.split()]
    for name in ("GF256Ctx", "gf256_init_", "gf256_add_mem", "gf256_addset_mem", "gf256_add2_mem",
                 "wirehair_v2_encode", "wirehair_v2_encoder_create_with_options", "wirehair_v2_free"):
        require(names.count(name) == 1, "single runtime/facade symbol count: " + name)
    runtimes = N.runtime_paths(N.checked(["/usr/bin/ldd", str(build / TARGET)], deadline))
    old = M.old_receipt(deadline)
    exact(runtimes, sorted(item["path"] for item in old["runtime_libraries"]), "runtime closure")
    return {"protocol": PROTOCOL, "source_commit": commit, "sources": sources,
            "historical_receipt_sha256": HISTORICAL_SHA, "build_directory": str(build),
            "build_files": [N.pin(path, deadline) for path in sorted(outputs | {build / "commands.json", build / "build_inputs.json"})],
            "dependencies": [N.pin(path, deadline, installed=True) for path in sorted(dependencies)],
            "helpers": [N.pin(path, deadline) for path in M.HELPERS], "symbols_sha256": C.sha(symbols),
            "interpreter": N.pin(Path(sys.executable).resolve(strict=True), deadline, installed=True),
            "python_version": list(sys.version_info[:3]), "worker_argv": [str(build / TARGET), "--worker"],
            "environment": C.ENV, "worker_seconds": 10, "outer_seconds": 20,
            "max_stdout_bytes": RAW_CAP, "max_stderr_bytes": ERR_CAP, "max_total_bytes": TOTAL_CAP}


def run(receipt_path):
    started = time.monotonic()
    receipt_bytes = C.read_regular(receipt_path, RAW_CAP, deadline=started + 4)
    receipt = C.strict_json(receipt_bytes)
    require(C.canonical(current_receipt(Path(receipt["build_directory"]), started + 4)) == receipt_bytes,
            "launch receipt differs")
    OUTPUT.mkdir(mode=0o700, exist_ok=False)
    members, total = {}, 0

    def publish(name, data):
        nonlocal total
        total += len(data)
        require(total <= TOTAL_CAP, "aggregate output cap")
        C.write_new(OUTPUT / name, data)
        members[name] = {"bytes": len(data), "sha256": C.sha(data)}

    publish("CLAIM.json", receipt_bytes)
    outcome, failure, analysis = "INVALID", None, None
    stdout, stderr, observation = b"", b"", {}
    try:
        observation = C.capture_worker(receipt["worker_argv"], min(time.monotonic() + 10, started + 14),
                                       stdout_cap=RAW_CAP, stderr_cap=ERR_CAP)
        stdout, stderr = observation.pop("stdout"), observation.pop("stderr")
        require(not observation["failure"] and observation["returncode"] == 0 and not stderr, "worker failed")
        duration = observation["elapsed_seconds"]
        require(type(duration) in (int, float) and math.isfinite(duration) and 0 < duration < 10, "observed worker bound")
        raw = C.strict_json(stdout)
        require(type(raw) is dict and type(raw.get("elapsed_ns")) is int and
                0 < raw["elapsed_ns"] <= duration * 10**9, "internal interval exceeds worker lifetime")
        require(stdout.endswith(b"\n"), "raw framing")
        analysis = reduce_raw(raw, deadline=started + 17)
        outcome = analysis["outcome"]
    except Exception as error:
        failure = type(error).__name__ + ": " + str(error)[:1000]
    try:
        require(C.canonical(current_receipt(Path(receipt["build_directory"]), started + 19)) == receipt_bytes,
                "post-run receipt differs")
        observation["post_pins"] = "match"
    except Exception as error:
        observation["post_pins"] = "failed_or_unobserved"
        outcome = "INVALID"
        failure = (failure + "; " if failure else "") + str(error)[:1000]
    publish("raw.json", stdout)
    publish("stderr.txt", stderr)
    if analysis is not None:
        publish("analysis.json", C.canonical(analysis))
    if time.monotonic() - started >= 19.5:
        outcome, failure = "INVALID", "controller deadline"
    summary = {"protocol": PROTOCOL, "outcome": outcome, "failure": failure, "worker": observation,
               "elapsed_seconds": time.monotonic() - started, "speed_claimed": False,
               "production_promotion_claimed": False, "WH1_compared": False}
    publish("summary.json", C.canonical(summary))
    publish("COMPLETE.json", C.canonical({"protocol": PROTOCOL, "outcome": outcome, "files": members}))
    fd = os.open(str(OUTPUT), os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)
    require(time.monotonic() - started < 20, "whole CLI deadline")
    print(C.canonical(summary).decode("ascii"))
    return int(outcome == "INVALID")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--build", action="store_true")
    mode.add_argument("--make-receipt", type=Path)
    mode.add_argument("--run", type=Path)
    mode.add_argument("--replay", type=Path)
    parser.add_argument("--build-dir", type=Path)
    args = parser.parse_args()
    if args.build or args.make_receipt:
        if args.build_dir is None:
            parser.error("build/receipt requires --build-dir")
        if args.build:
            build_only(args.build_dir)
        else:
            C.write_new(args.make_receipt, C.canonical(current_receipt(args.build_dir)))
    elif args.build_dir is not None:
        parser.error("--build-dir only belongs to build/receipt")
    elif args.run:
        return run(args.run)
    else:
        raw = C.strict_json(C.read_regular(args.replay, RAW_CAP))
        print(C.canonical(reduce_raw(raw)).decode("ascii"))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print(type(error).__name__ + ": " + str(error)[:1000], file=sys.stderr)
        sys.exit(1)
