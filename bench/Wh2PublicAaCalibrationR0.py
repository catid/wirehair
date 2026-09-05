#!/usr/bin/env python3
"""One-shot A/A calibration (.50); no WH2/WH1 comparison or r1 rerun.

The private r1 helper instance supplies already-reviewed probes and bounded
process supervision. Its on-disk sources and the old evidence remain immutable.
Replay checks trusted-producer consistency, not offline cryptographic authority.
"""
import argparse
import hashlib
import importlib.util
import math
import os
from pathlib import Path
import shlex
import stat
import sys
import time


_spec = importlib.util.spec_from_file_location(
    "_aa_r0_helpers", Path(__file__).with_name("Wh2PublicBorrowedCurrentVsWh1R1.py"))
h = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(h)
CAMPAIGN = "wh2-public-aa-calibration-r0"
PREFIX = "wirehair.wh2.public-aa-calibration-r0"
FIXED_OUTPUT_DIR = Path("/var/tmp/" + CAMPAIGN)
TARGET = "wirehair_wh2_public_aa_calibration_r0"
TUPLES = (
    (8, 64, "prebuilt-systematic"), (8, 1280, "prebuilt-systematic"),
    (8, 1280, "prebuilt-repair"), (128, 1280, "prebuilt-systematic"),
    (128, 1280, "fresh-repair"), (8192, 1280, "prebuilt-repair"),
    (8192, 1280, "prebuilt-systematic"), (128, 1280, "fresh-systematic"),
)
CELLS = tuple(sorted({(k, b) for k, b, _ in TUPLES}))
CONDITION_SEQUENCES = ((0, 1, 3, 2), (1, 2, 0, 3), (2, 3, 1, 0), (3, 0, 2, 1))
SOURCES = ("bench/Wh2PublicAaCalibrationR0.cpp",) + h.R1_TARGET_SOURCES[1:]
INPUTS = h.R1_TRACKED_INPUTS + (
    SOURCES[0], "bench/Wh2PublicAaCalibrationR0.py",
    "bench/test_Wh2PublicAaCalibrationR0.py",
)
CORE_FILES = ("CLAIM", "raw.jsonl", "worker.stderr", "provenance.json", "summary.json", "COMPLETE")
COUNTS = {
    "expected_panels": 768, "expected_records": 770,
    "expected_measured_batches": 6144, "expected_warmup_batches": 1536,
    "expected_measured_invocations": 2411520, "expected_warmup_invocations": 602880,
    "expected_measured_encode_calls": 37945344,
    "expected_warmup_encode_calls": 9486336, "expected_encode_calls": 47431680,
}
# This module instance is never shared with an r1 controller. The legacy
# supervisor's config-identity argument carries our immutable describe hash.
h.CAMPAIGN = CAMPAIGN
h.FIXED_OUTPUT_DIR = FIXED_OUTPUT_DIR
h.WORKER_STARTED_SCHEMA = PREFIX + ".worker-started.v1"
h.R0_ROOTS += (Path("/var/tmp/wh2-public-borrowed-current-vs-wh1-r1"),)
h.R1_TRACKED_INPUTS = INPUTS
fail = h.fail
canonical = h.canonical_bytes
digest = h.sha256_bytes


def equal(actual, expected, where):
    # Canonical bytes distinguish bool from int, unlike Python equality.
    if canonical(actual) != canonical(expected):
        fail(where + " differs")


def obj(value, keys, where):
    if type(value) is not dict:
        fail(where + " is not an object")
    h.exact_keys(value, keys, where)
    return value


def unit_order(replicate):
    return sorted(range(16), key=lambda unit: digest(
        "aa-cal-r0:rep={}:unit={}".format(replicate, unit).encode("ascii")))


def roster():
    for replicate in range(12):
        for unit in unit_order(replicate):
            for position, condition in enumerate(CONDITION_SEQUENCES[replicate % 4]):
                yield replicate, unit, position, condition


def batch_invocations(k, scope):
    return 1 if scope.startswith("fresh-") else (8192 + k - 1) // k


def environment(library):
    return {"LANG": "C", "LC_ALL": "C", "PATH": "/usr/bin:/bin", "TZ": "UTC",
            "LD_LIBRARY_PATH": str(Path(library).parent)}


def expected_build_tokens(provenance, invocation):
    source = Path(provenance["source_root"])
    build = Path(provenance["build_root"])
    compiler = provenance["compiler"]
    definitions = ["-DWIREHAIR_DLL=1", '-DWIREHAIR_WH2_CMAKE_CXX_COMPILER_ID="GNU"',
                   '-DWIREHAIR_WH2_CMAKE_CXX_COMPILER_VERSION="{}"'.format(compiler["version"]),
                   '-DWIREHAIR_WH2_CXX_COMPILER_PATH="{}"'.format(compiler["path"]),
                   '-DWIREHAIR_WH2_CXX_COMPILER_SHA256="{}"'.format(compiler["sha256"]),
                   '-DWIREHAIR_WH2_SOURCE_GIT_COMMIT="{}"'.format(provenance["source_commit"])]
    result, objects = [], []
    for name, sources, defs, flags in (
        (TARGET, SOURCES, definitions, []),
        ("wirehair", h.R1_LIBRARY_SOURCES,
         ["-DWIREHAIR_BUILDING=1", "-DWIREHAIR_DLL=1", "-Dwirehair_EXPORTS"], ["-fPIC"]),
    ):
        group = []
        for relative in sources:
            target = "CMakeFiles/{}.dir/{}.o".format(name, relative)
            group.append(target)
            result.append([invocation] + defs + ["-I" + str(source / "include"),
                          "-O3", "-DNDEBUG", "-std=gnu++11"] + flags +
                          ["-Wall", "-Wextra", "-march=native", "-MD", "-MT", target,
                           "-MF", target + ".d", "-o", target, "-c", str(source / relative)])
        objects.append(group)
    library = "libwirehair.so.2.0.0"
    result.append([":", "&&", invocation, "-O3", "-DNDEBUG"] + objects[0] +
                  ["-o", TARGET, "-Wl,-rpath," + str(build), library, "-ldl", "&&", ":"])
    result.append([":", "&&", invocation, "-fPIC", "-O3", "-DNDEBUG",
                   "-Wl,--version-script=" + str(source / "abi/wirehair.map"),
                   "-shared", "-Wl,-soname,libwirehair.so.2", "-o", library] +
                  objects[1] + ["-lm", "&&", ":"])
    return result


BUILD_CACHE = {"BUILD_SHARED_LIBS": "ON", "CMAKE_BUILD_TYPE": "Release",
               "CMAKE_GENERATOR": "Ninja", "MARCH_NATIVE": "ON", "WH_LTO": "OFF",
               "WH_PGO_MODE": "OFF"}


def validate_build(value, provenance):
    obj(value, {"cache", "compiler_invocation", "commands", "commands_sha256", "dry_run"}, "build")
    equal(value["cache"], BUILD_CACHE, "build.cache")
    h.replay_path(value["compiler_invocation"], "build.compiler_invocation")
    lines = value["dry_run"].splitlines()
    if not (len(lines) in (1, 2) and lines[-1] == "ninja: no work to do." and
            (len(lines) == 1 or lines[0].startswith("ninja: Entering directory `"))):
        fail("build is not up to date")
    h.replay_text_hash(value, "commands", "commands_sha256", "build")
    rows = [shlex.split(line) for line in value["commands"].splitlines()]
    expected = expected_build_tokens(provenance, value["compiler_invocation"])
    # Ninja also emits the standard DSO symlink helper. All actual compiler
    # commands (including unexpected injected TUs) must match the exact roster.
    actual = [row for row in rows if "-c" in row or "-o" in row]
    equal(sorted(actual), sorted(expected), "build compiler command roster")
    extras = [row for row in rows if "-c" not in row and "-o" not in row]
    expected_extra = [["/usr/bin/cmake", "-E", "cmake_symlink_library",
                       "libwirehair.so.2.0.0", "libwirehair.so.2", "libwirehair.so", "&&", ":"]]
    equal(extras, expected_extra, "build helper command roster")


def build_receipt(provenance):
    build = Path(provenance["build_root"])
    cache = h.parse_cache(build / "CMakeCache.txt")
    for key, expected in (("CMAKE_HOME_DIRECTORY", provenance["source_root"]),
                          ("CMAKE_CACHEFILE_DIR", str(build)),
                          ("CMAKE_CXX_COMPILER", provenance["compiler"]["path"]),
                          ("_Python3_EXECUTABLE", provenance["interpreter"]["path"])):
        equal(str(Path(cache.get(key, "")).resolve()), expected, "build." + key)
    commands = h.run_checked(["/usr/bin/ninja", "-C", str(build), "-t", "commands", TARGET]).decode("utf-8")
    rows = [shlex.split(line) for line in commands.splitlines()]
    invocation = next(row[0] for row in rows if "-c" in row)
    equal(str(Path(invocation).resolve(strict=True)), provenance["compiler"]["path"], "build compiler")
    result = {"cache": {key: cache.get(key) for key in BUILD_CACHE},
              "compiler_invocation": invocation, "commands": commands,
              "commands_sha256": digest(commands.encode("utf-8")),
              "dry_run": h.run_checked(["/usr/bin/ninja", "-C", str(build), "-n", TARGET]).decode("utf-8")}
    validate_build(result, provenance)
    return result


def statistics(samples):
    result, failed = {"conditions": [], "matched_replicates": []}, []
    boundary = math.log1p(0.02)
    for unit in range(16):
        for condition in range(4):
            values = [samples[(rep, unit, condition)] for rep in range(12)]
            summary = h.sample_summary(values)
            result["conditions"].append({"unit": unit, "condition": condition, "summary": summary})
            if not (-boundary < summary["lower95_log"] and summary["upper95_log"] < boundary):
                failed.append("unit={}:condition={}".format(unit, condition))
        paired = [math.fsum(samples[(rep, unit, cond)] for cond in range(4)) / 4
                  for rep in range(12)]
        summary = h.sample_summary(paired)
        result["matched_replicates"].append({"unit": unit, "summary": summary})
        if not (-boundary < summary["lower95_log"] and summary["upper95_log"] < boundary):
            failed.append("unit={}:matched".format(unit))
    return result, failed


def description():
    value = dict(COUNTS)
    value.update({
        "campaign": CAMPAIGN, "schema": PREFIX + ".describe.v1",
        "condition_sequences": [list(row) for row in CONDITION_SEQUENCES],
        "conditions": [{"condition": cond, "mapping": "direct" if cond < 2 else "swapped",
                        "order": "ABBA" if cond % 2 == 0 else "BAAB",
                        "warmup": ["left", "right"] if cond % 2 == 0 else ["right", "left"]}
                       for cond in range(4)],
        "internal_deadline_seconds": 110, "panel_replicates": 12,
        "roles": ["C", "L"], "sibling_cpu": 56, "target_cpu": 120,
        "source_seed_base": h.SOURCE_SEED_BASE_TEXT, "timing_granularity": "whole-batch",
        "tuples": [{"K": k, "block_bytes": b, "scope": scope,
                    "scope_invocations_per_batch": batch_invocations(k, scope)}
                   for k, b, scope in TUPLES],
        "unit_order": "sha256:aa-cal-r0:rep=R:unit=U",
    })
    return value


def description_hash():
    return digest(canonical(description()) + b"\n")


def validate_config(value, expected, described):
    obj(value, {"campaign", "cells", "compile", "description_sha256",
                "runtime_library_maps_sha256", "runtime_library_path", "schema",
                "target_identity"}, "config")
    equal(described, description(), "description")
    for key, wanted in (("campaign", CAMPAIGN), ("schema", PREFIX + ".config.v1"),
                        ("description_sha256", description_hash()),
                        ("runtime_library_path", expected["library_path"])):
        equal(value[key], wanted, "config." + key)
    h.lower_hash(value["runtime_library_maps_sha256"], "config maps hash")
    h.validate_target_identity(value["target_identity"], "config.target_identity")
    equal(value["compile"], {
        "compiler_path": expected["compiler_path"], "compiler_sha256": expected["compiler_sha256"],
        "compiler_version": expected["compiler_version"],
        "harness_git_commit": expected["source_commit"],
        "implementation_git_commit": expected["source_commit"],
    }, "config.compile")
    cells = value["cells"]
    if type(cells) is not list or len(cells) != len(CELLS):
        fail("config cell roster differs")
    oracles = {}
    for cell, (k, b) in zip(cells, CELLS):
        obj(cell, h.CELL_KEYS, "cell")
        source_hash = h.SOURCE_SHA256[k, b, b]
        for key, wanted in (("K", k), ("block_bytes", b), ("final_bytes", b),
                            ("message_bytes", k * b), ("source_sha256", source_hash),
                            ("source_seed", h.source_seed_text(k, b, b)),
                            ("invocations_by_replicate", h.invocations_by_replicate(k))):
            equal(cell[key], wanted, "cell." + key)
        h.validate_partial_oracles(cell["partial_final_oracles"], k, b, "cell partial oracles")
        if type(cell["oracles"]) is not list or len(cell["oracles"]) != 2:
            fail("cell oracle roster differs")
        for arm, oracle in zip(("C", "L"), cell["oracles"]):
            oracles[k, b, arm] = h.validate_oracle(oracle, arm, k, source_hash, "cell oracle")
        if oracles[k, b, "C"]["full_repair_sha256"] == oracles[k, b, "L"]["full_repair_sha256"]:
            fail("distinct public equations have matching repair streams")
    return oracles


PANEL_KEYS = {"K", "block_bytes", "condition", "description_sha256", "left_output_sha256",
              "mapping", "order", "physical_scratch_addresses", "replicate", "right_output_sha256",
              "role", "runtime_library_maps_sha256", "schema", "scope", "scope_invocations_per_batch",
              "sequence", "slots", "source_immutable_verified", "target_cpu", "tuple_index",
              "unit_index", "unit_key_sha256"}


def validate_panel(value, key, sequence, config, oracles):
    obj(value, PANEL_KEYS, "panel")
    replicate, unit, _, condition = key
    k, b, scope = TUPLES[unit // 2]
    role = ("C", "L")[unit % 2]
    order = "ABBA" if condition % 2 == 0 else "BAAB"
    wanted = {"K": k, "block_bytes": b, "scope": scope, "condition": condition,
              "description_sha256": description_hash(), "replicate": replicate,
              "role": role, "mapping": "direct" if condition < 2 else "swapped", "order": order,
              "runtime_library_maps_sha256": config["runtime_library_maps_sha256"],
              "schema": PREFIX + ".panel.v1", "sequence": sequence,
              "scope_invocations_per_batch": batch_invocations(k, scope),
              "source_immutable_verified": True, "target_cpu": 120, "tuple_index": unit // 2,
              "unit_index": unit,
              "unit_key_sha256": digest("aa-cal-r0:rep={}:unit={}".format(replicate, unit).encode("ascii"))}
    output_hash = oracles[k, b, role]["systematic_sha256" if scope.endswith("systematic") else "full_repair_sha256"]
    wanted.update(left_output_sha256=output_hash, right_output_sha256=output_hash)
    for field, expected in wanted.items():
        equal(value[field], expected, "panel." + field)
    addresses = value["physical_scratch_addresses"]
    if type(addresses) is not list or len(addresses) != 2 or any(
        type(address) is not str or h.re.fullmatch(r"0x[0-9a-f]{1,16}", address) is None
        for address in addresses
    ) or any(int(address, 16) == 0 for address in addresses) or addresses[0] == addresses[1]:
        fail("panel physical scratch addresses differ")
    if abs(int(addresses[0], 16) - int(addresses[1], 16)) < k * b:
        fail("physical scratch buffers overlap")
    slots = value["slots"]
    if type(slots) is not list or len(slots) != 8:
        fail("panel slot roster differs")
    elapsed = []
    for slot, side in zip(slots, h.sides_for(order)):
        obj(slot, {"elapsed_ns", "logical_lane", "physical_lane"}, "slot")
        equal(slot["logical_lane"], side, "slot logical lane")
        equal(slot["physical_lane"], (0 if side == "left" else 1) ^ (condition // 2), "slot physical lane")
        elapsed.append(h.exact_int(slot["elapsed_ns"], 1, h.MAX_INT63, "slot elapsed"))
    return h.lane_contrast(elapsed, order)


def validate_terminal(value, config):
    equal(value, {"campaign": CAMPAIGN, "encode_call_count": 47431680,
                  "measured_batch_count": 6144, "measured_invocation_count": 2411520,
                  "panel_count": 768, "record_count": 770, "schema": PREFIX + ".terminal.v1",
                  "status": "complete", "warmup_batch_count": 1536,
                  "warmup_invocation_count": 602880}, "terminal")


def parse_transcript(raw, expected, description):
    if len(raw) > h.MAX_RAW_BYTES:
        fail("raw exceeds byte limit")
    lines = raw.splitlines(keepends=True)
    if len(lines) != 770:
        fail("transcript record roster differs")
    records = [h.parse_canonical_line(line, "record {}".format(index))
               for index, line in enumerate(lines)]
    config, terminal = records[0], records[-1]
    oracles = validate_config(config, expected, description)
    samples, pair_addresses = {}, {}
    for index, (value, key) in enumerate(zip(records[1:-1], roster())):
        contrast = validate_panel(value, key, index, config, oracles)
        pair = key[:2]
        addresses = value["physical_scratch_addresses"]
        if pair in pair_addresses:
            equal(addresses, pair_addresses[pair], "unit retained physical pair")
        pair_addresses[pair] = addresses
        samples[(key[0], key[1], key[3])] = contrast
    validate_terminal(terminal, config)
    return statistics(samples)


def artifacts(provenance):
    return {key: h.artifact_receipt(Path(provenance[key + "_path"]))
            for key in ("worker", "library", "controller", "compiler")}


def preflight(args):
    source = h.exact_absolute_directory(args.source_root, "source root")
    build = h.exact_absolute_directory(args.build_root, "build root")
    paths = {key + "_path": str(h.exact_absolute_file(getattr(args, key), key))
             for key in ("worker", "library", "compiler")}
    paths["controller_path"] = str(Path(__file__).resolve())
    equal(paths["controller_path"], str(source / "bench/Wh2PublicAaCalibrationR0.py"), "controller entrypoint")
    equal(paths["worker_path"], str(build / TARGET), "worker build path")
    equal(paths["library_path"], str(build / "libwirehair.so.2.0.0"), "library build path")
    h.lower_commit(args.expected_source_commit, "expected source commit")
    p = dict(paths, source_root=str(source), build_root=str(build), source_commit=args.expected_source_commit)
    p["compiler"] = h.compiler_receipt(Path(p["compiler_path"]))
    p["interpreter"] = h.interpreter_receipt()
    h.validate_interpreter_receipt(p["interpreter"], "interpreter")
    p["artifacts_before"] = artifacts(p)
    for key in ("worker", "library", "controller"):
        expected_hash = getattr(args, "expected_" + key + "_sha256")
        h.lower_hash(expected_hash, "expected " + key + " hash")
        equal(p["artifacts_before"][key]["sha256"], expected_hash, key + " hash")
    p["expected"] = {"compiler_path": p["compiler_path"], "compiler_sha256": p["compiler"]["sha256"],
                     "compiler_version": p["compiler"]["version"], "source_commit": p["source_commit"],
                     "library_path": p["library_path"]}
    p["git_before"] = h.git_receipt(source, p["source_commit"])
    p["build_before"] = build_receipt(p)
    p["dynamic_before"] = h.dynamic_receipt(Path(p["worker_path"]), Path(p["library_path"]))
    described = h.run_checked([p["worker_path"], "--describe"], maximum=h.MAX_LINE_BYTES,
                             env=environment(p["library_path"]), require_empty_stderr=True)
    p["description"] = h.parse_canonical_line(described, "description")
    equal(p["description"], description(), "worker description")
    equal(digest(described), description_hash(), "worker description hash")
    p["preserved_before"] = h.snapshot_r0()
    p["health_before"] = h.cpu_receipt()
    equal(artifacts(p), p["artifacts_before"], "preflight stable artifacts")
    return p


def marker_value(claim, claim_hash):
    return {"campaign": CAMPAIGN, "claim_sha256": claim_hash,
            "description_sha256": description_hash(), "schema": PREFIX + ".worker-started.v1",
            "source_commit": claim["source_commit"], "status": "started",
            "worker_sha256": claim["worker_sha256"]}


def parse_marker(data, claim):
    try:
        marker = h.parse_canonical_line(data, "WORKER_STARTED")
        equal(marker, marker_value(claim, digest(canonical(claim) + b"\n")), "WORKER_STARTED")
        return marker, None
    except (h.ValidationError, ValueError, TypeError, KeyError) as exc:
        return None, "worker marker:" + str(exc)


def validate_provenance(p, claim):
    for key in ("source_root", "build_root", "worker_path", "library_path", "compiler_path", "controller_path"):
        h.replay_path(p[key], "provenance." + key)
    equal(p["source_commit"], claim["source_commit"], "provenance commit")
    equal(p["worker_path"], str(Path(p["build_root"]) / TARGET), "worker build path")
    equal(p["library_path"], str(Path(p["build_root"]) / "libwirehair.so.2.0.0"), "library build path")
    equal(p["controller_path"], str(Path(p["source_root"]) / "bench/Wh2PublicAaCalibrationR0.py"), "controller entrypoint")
    equal(p["description"], description(), "provenance description")
    equal(p["expected"], {"compiler_path": p["compiler_path"], "compiler_sha256": p["compiler"]["sha256"],
                           "compiler_version": p["compiler"]["version"], "source_commit": p["source_commit"],
                           "library_path": p["library_path"]}, "expected native identity")
    h.validate_interpreter_receipt(p["interpreter"], "interpreter")
    obj(p["compiler"], h.COMPILER_RECEIPT_KEYS, "compiler")
    h.replay_text_hash(p["compiler"], "version_text", "version_text_sha256", "compiler")
    if h.re.fullmatch(r"13(?:\.[0-9]+){1,2}", p["compiler"]["version"]) is None:
        fail("compiler is not GCC 13")
    equal(p["compiler"]["path"], p["compiler_path"], "compiler path")
    equal(p["compiler"]["version_sha256"], digest((p["compiler"]["version"] + "\n").encode("ascii")), "compiler short version")
    equal(p["compiler"]["sha256"], p["artifacts_before"]["compiler"]["sha256"], "compiler hash")
    for key in ("worker", "library", "controller", "compiler"):
        h.validate_stat_receipt(p["artifacts_before"][key], "artifact " + key, p[key + "_path"], artifact=True)
        if key != "compiler":
            equal(p["artifacts_before"][key]["sha256"], claim[key + "_sha256"], "claim " + key)
    equal(p["interpreter"]["artifact"]["sha256"], claim["controller_interpreter_sha256"], "claim interpreter")
    h.validate_git_receipt(p["git_before"], p["source_commit"], "git")
    validate_build(p["build_before"], p)
    h.validate_dynamic_receipt(p["dynamic_before"], p, "dynamic")
    h.validate_r0_receipt(p["preserved_before"], "preserved campaigns")
    h.validate_health_receipt(p["health_before"], "health before")
    affinity = p["controller_affinity"]
    obj(affinity, {"affinity_after", "affinity_before", "controller_cpu", "sibling_available_before",
                   "singleton_verified", "target_available_before"}, "controller affinity")
    inherited = affinity["affinity_before"]
    if type(inherited) is not list or any(type(cpu) is not int for cpu in inherited) or inherited != sorted(set(inherited)) or not {0, 56, 120}.issubset(inherited):
        fail("inherited controller affinity differs")
    equal(affinity, dict(affinity_before=inherited, affinity_after=[0], controller_cpu=0,
                         sibling_available_before=True, singleton_verified=True, target_available_before=True), "controller affinity")


def adjudicate_bundle(raw, stderr, p, marker, claim, elapsed):
    failures, stats, failed = [], None, []
    if type(elapsed) not in (int, float) or not math.isfinite(elapsed) or elapsed < 0:
        fail("elapsed receipt is malformed")
    if elapsed > 120:
        failures.append("outer_deadline")
    process = obj(p["process"], {"exit_code", "elapsed_seconds", "timed_out", "error"}, "process")
    if type(process["timed_out"]) is not bool:
        fail("process timeout receipt is malformed")
    if process["exit_code"] is not None:
        h.exact_int(process["exit_code"], -255, 255, "process exit")
    if type(process["elapsed_seconds"]) not in (int, float) or not math.isfinite(process["elapsed_seconds"]) or process["elapsed_seconds"] < 0:
        fail("process elapsed is malformed")
    if process["timed_out"]:
        failures.append("worker_deadline")
    if process["error"] is not None:
        failures.append("execution:" + h.bounded_string(process["error"], "execution error"))
    if process["exit_code"] != 0:
        failures.append("worker_exit:" + str(process["exit_code"]))
    if stderr:
        failures.append("worker_stderr")
    if marker is None:
        failures.append("worker_started_missing")
    else:
        equal(marker, marker_value(claim, digest(canonical(claim) + b"\n")), "WORKER_STARTED")
    obj(p["postflight_errors"], (key for key in POSTFLIGHT), "postflight errors")
    for key in POSTFLIGHT:
        error = p["postflight_errors"][key]
        if error is not None:
            h.bounded_string(error, "postflight error")
            equal(p[key + "_after"], None, "failed postflight receipt")
            failures.append("postflight_" + key + ":" + error)
            continue
        try:
            if key == "health":
                h.validate_health_pair(p["health_before"], p["health_after"])
            else:
                equal(p[key + "_after"], p[key + "_before"] if key != "interpreter" else p["interpreter"],
                      "stable " + key)
        except h.ValidationError as exc:
            failures.append("postflight_" + key + ":" + str(exc))
    try:
        stats, failed = parse_transcript(raw, p["expected"], p["description"])
    except (h.ValidationError, KeyError, TypeError, ValueError, IndexError, OverflowError) as exc:
        failures.append("transcript:" + str(exc))
    return {"campaign": CAMPAIGN, "schema": PREFIX + ".summary.v1",
            "outcome": "invalid" if failures else "fail" if failed else "pass",
            "elapsed_seconds": elapsed, "infrastructure_failures": failures,
            "failed_aa_gates": failed, "statistics": stats, "raw_sha256": digest(raw)}


POSTFLIGHT = {
    "health": lambda p: h.cpu_receipt(), "preserved": lambda p: h.snapshot_r0(),
    "artifacts": artifacts, "interpreter": lambda p: h.interpreter_receipt(),
    "git": lambda p: h.git_receipt(Path(p["source_root"]), p["source_commit"]),
    "build": build_receipt,
    "dynamic": lambda p: h.dynamic_receipt(Path(p["worker_path"]), Path(p["library_path"])),
}


def run_once(args):
    affinity = h.pin_controller()
    p = preflight(args)
    p["controller_affinity"] = affinity
    claim = {"campaign": CAMPAIGN, "schema": PREFIX + ".claim.v1",
             "created_unix_ns": time.time_ns(), "description_sha256": description_hash(),
             "source_commit": p["source_commit"],
             "controller_interpreter_sha256": p["interpreter"]["artifact"]["sha256"]}
    claim.update({key + "_sha256": p["artifacts_before"][key]["sha256"]
                  for key in ("worker", "library", "controller")})
    validate_provenance(p, claim)
    directory = h.claim_namespace(claim)
    started = time.monotonic()
    raw, stderr, marker = b"", b"", None
    try:
        p["process"] = {"exit_code": None, "elapsed_seconds": 0.0, "timed_out": False, "error": None}
        try:
            raw, stderr, code, elapsed, timeout = h.run_worker(
                p["worker_path"], p["library_path"], started + 120,
                digest(canonical(claim) + b"\n"), claim["worker_sha256"], description_hash())
            p["process"].update(exit_code=code, elapsed_seconds=elapsed, timed_out=timeout)
        except Exception as exc:
            p["process"].update(error=str(exc), elapsed_seconds=time.monotonic() - started)
        p["postflight_errors"] = {}
        for key, probe in POSTFLIGHT.items():
            p[key + "_after"], p["postflight_errors"][key] = None, None
            try:
                if time.monotonic() >= started + 120:
                    fail("outer deadline prevented postflight probe")
                p[key + "_after"] = probe(p)
            except Exception as exc:
                p["postflight_errors"][key] = str(exc)
        marker_path = FIXED_OUTPUT_DIR / "WORKER_STARTED"
        p["marker_sha256"] = None
        if marker_path.exists():
            try:
                with marker_path.open("rb") as source:
                    marker_data = source.read(h.MAX_PROBE_BYTES + 1)
                p["marker_sha256"] = digest(marker_data)
                marker, marker_error = parse_marker(marker_data, claim)
                if marker_error:
                    p["process"]["error"] = marker_error
            except Exception as exc:
                p["process"]["error"] = "worker marker:" + str(exc)
                marker = None
        h.write_exclusive(directory, "raw.jsonl", raw, 0o400)
        h.write_exclusive(directory, "worker.stderr", stderr, 0o400)
        h.publish_json(directory, "provenance.json", p)
        elapsed = time.monotonic() - started
        summary = adjudicate_bundle(raw, stderr, p, marker, claim, elapsed)
        # Adjudication itself is part of the bound, including JSON parsing/CIs.
        summary = dict(summary, elapsed_seconds=time.monotonic() - started)
        if summary["elapsed_seconds"] > 120 and "outer_deadline" not in summary["infrastructure_failures"]:
            summary["infrastructure_failures"].insert(0, "outer_deadline")
            summary["outcome"] = "invalid"
        h.publish_json(directory, "summary.json", summary)
        complete = {"campaign": CAMPAIGN, "schema": PREFIX + ".complete.v1", "status": "complete",
                    "outcome": summary["outcome"],
                    "files": {name: h.file_sha256(FIXED_OUTPUT_DIR / name, h.MAX_RAW_BYTES)
                              for name in CORE_FILES if name != "COMPLETE"}}
        if marker_path.exists():
            complete["files"]["WORKER_STARTED"] = p["marker_sha256"]
        h.publish_json(directory, "COMPLETE", complete)
        os.fsync(directory)
        print(canonical(summary).decode("ascii"))
        return {"pass": 0, "fail": 2, "invalid": 1}[summary["outcome"]]
    finally:
        os.close(directory)


def replay(bundle):
    bundle = h.exact_absolute_directory(str(bundle), "bundle")
    names = sorted(path.name for path in bundle.iterdir())
    if names not in (sorted(CORE_FILES), sorted(CORE_FILES + ("WORKER_STARTED",))):
        fail("bundle roster differs")
    data = {}
    for name in names:
        path = h.exact_absolute_file(str(bundle / name), name)
        info = path.stat()
        if stat.S_IMODE(info.st_mode) != 0o400 or info.st_nlink != 1:
            fail("bundle metadata differs for " + name)
        cap = h.MAX_RAW_BYTES if name == "raw.jsonl" else h.MAX_STDERR_BYTES if name == "worker.stderr" else h.MAX_PROBE_BYTES
        with path.open("rb") as source:
            data[name] = source.read(cap + 1)
        if len(data[name]) > cap:
            fail("bundle member exceeds byte cap: " + name)
    objects = {name: h.parse_canonical_line(value, name) for name, value in data.items()
               if name not in ("raw.jsonl", "worker.stderr", "WORKER_STARTED")}
    claim, p = objects["CLAIM"], objects["provenance.json"]
    obj(claim, {"campaign", "schema", "created_unix_ns", "description_sha256", "source_commit",
                "controller_interpreter_sha256", "worker_sha256", "library_sha256", "controller_sha256"}, "CLAIM")
    equal(claim["campaign"], CAMPAIGN, "claim campaign")
    equal(claim["schema"], PREFIX + ".claim.v1", "claim schema")
    equal(claim["description_sha256"], description_hash(), "claim description")
    h.exact_int(claim["created_unix_ns"], 1, h.MAX_INT63, "claim timestamp")
    h.lower_commit(claim["source_commit"], "claim commit")
    validate_provenance(p, claim)
    marker, marker_error = parse_marker(data["WORKER_STARTED"], claim) if "WORKER_STARTED" in data else (None, None)
    if marker_error:
        equal(p["process"]["error"], marker_error, "recorded marker failure")
    equal(p["marker_sha256"], digest(data["WORKER_STARTED"]) if "WORKER_STARTED" in data else None, "marker hash")
    summary = adjudicate_bundle(data["raw.jsonl"], data["worker.stderr"], p, marker, claim,
                                objects["summary.json"]["elapsed_seconds"])
    equal(objects["summary.json"], summary, "recomputed summary")
    equal(objects["COMPLETE"], {"campaign": CAMPAIGN, "schema": PREFIX + ".complete.v1",
                               "status": "complete", "outcome": summary["outcome"],
                               "files": {name: digest(value) for name, value in data.items() if name != "COMPLETE"}},
          "COMPLETE")
    print(canonical(summary).decode("ascii"))
    return {"pass": 0, "fail": 2, "invalid": 1}[summary["outcome"]]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--run", action="store_true")
    mode.add_argument("--replay", type=Path)
    mode.add_argument("--selftest", action="store_true")
    for name in ("source-root", "build-root", "worker", "library", "compiler", "expected-source-commit",
                 "expected-worker-sha256", "expected-library-sha256", "expected-controller-sha256"):
        parser.add_argument("--" + name)
    args = parser.parse_args(argv)
    try:
        if args.replay:
            return replay(args.replay)
        if args.selftest:
            equal(len(list(roster())), 768, "roster count")
            print("A/A calibration controller selftest passed")
            return 0
        required = vars(args).copy()
        for key in ("run", "replay", "selftest"):
            required.pop(key)
        if any(value is None for value in required.values()):
            fail("--run requires all exact provenance arguments")
        return run_once(args)
    except (h.ValidationError, OSError, ValueError, TypeError, KeyError, IndexError, OverflowError) as exc:
        print("INVALID: " + str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
