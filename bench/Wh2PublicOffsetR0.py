#!/usr/bin/env python3
"""One bounded candidate/current public-encode screen; never a WH1 speed claim.

Reuse the immutable r1 process supervisor, clock statistics and host checks in
a private module. Both dynamically loaded codecs share the exact I/O addresses.
Replay checks trusted-producer consistency, not cryptographic attestation.
"""
import argparse
import importlib.util
import math
import os
from pathlib import Path
import shlex
import struct
import sys
import time

ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location("_offset_r0_helpers",
    ROOT / "bench/Wh2PublicBorrowedCurrentVsWh1R1.py")
h = importlib.util.module_from_spec(spec)
spec.loader.exec_module(h)
canonical, digest, fail = h.canonical_bytes, h.sha256_bytes, h.fail
CAMPAIGN = "wh2-public-offset-r0"
PREFIX = "wirehair.wh2.public-offset-r0"
OUTPUT = Path("/var/tmp/" + CAMPAIGN)
BASELINE = Path("/tmp/wh2-page-aa.uWaCGu/libwirehair.so.2.0.0")
BASELINE_HASH = "c0bf0666e3cc51523e18847b8a07384f2eac312518d0b4f9ac48cd14e63fa038"
BASIS_PATH = Path("/var/tmp/wh2-production-mix3-k3k6-r0/freeze.json")
BASIS_HASH = "559050987e6b683c85f8b81dc69d33c1d4e2b1e54f88b364a703ae8c17384de0"
BASIS_BUILD_HASH = "beae9fbb37d15e5412810840c13eb1a29ac12d792627a1170eca1b94f16267f9"
LIBC = Path("/usr/lib/x86_64-linux-gnu/libc.so.6")
LIBC_HASH = "8db37cf3f2169f59a0f07ef1fea308c35656668c64c8ff294e1860f4121eb161"
CANDIDATE_HASH = "a5d5f51d47fe7def0862e93e8c8b526fa8d8989c38e904525d18d3b145154ff8"
CELLS = tuple((k, b) for k in (8, 128, 8192) for b in (64, 1280))
CONDITIONS = ((0, 1, 3, 2), (1, 2, 0, 3), (2, 3, 1, 0), (3, 0, 2, 1))
SOURCES = ("bench/Wh2PublicOffsetR0.cpp", "bench/Wh2FrozenTrace.cpp",
           "bench/Wh2PublicBorrowedTargetIdentity.cpp", "bench/Wh2RdpruTargetIdentityV2.cpp")
INPUTS = SOURCES + ("bench/Wh2PublicOffsetR0.py", "bench/test_Wh2PublicOffsetR0.py",
    "bench/Wh2FrozenTrace.h", "bench/Wh2PublicBorrowedTargetIdentity.h",
    "bench/Wh2RdpruTargetIdentityV2.h", "bench/Wh2PublicBorrowedCurrentVsWh1R1.py")
h.CAMPAIGN, h.FIXED_OUTPUT_DIR = CAMPAIGN, OUTPUT
h.R1_TRACKED_INPUTS = INPUTS
h.R0_ROOTS += tuple(Path("/var/tmp/" + name) for name in (
    "wh2-public-borrowed-current-vs-wh1-r1", "wh2-public-aa-calibration-r0",
    "wh2-public-page-aa-r0", "wh2-production-mix3-k3k6-r0",
    "wh2-production-mix3-k3k6-broad-r0"))
FILES = {"CLAIM", "raw.jsonl", "worker.stderr", "provenance.json", "summary.json"}


def equal(actual, expected, where):
    if canonical(actual) != canonical(expected):
        fail(where + " differs")


def obj(value, keys, where):
    if type(value) is not dict:
        fail(where + " is not an object")
    h.exact_keys(value, keys, where)
    return value


def file_hash(path):
    return h.file_sha256(Path(path), h.MAX_RAW_BYTES)


def command(argv):
    return h.run_checked(argv, require_empty_stderr=True).decode("utf-8")


def source_receipt():
    files = command(["git", "-C", str(ROOT), "ls-files", "--", "CMakeLists.txt",
        "cmake", "codec", "include", "tables", "abi/wirehair.map", "wirehair.cpp",
        "Wirehair*", "gf256.cpp", "gf256.h"]).splitlines()
    return {path: file_hash(ROOT / path) for path in sorted(set(files) | set(INPUTS))}


def validate_basis(basis):
    equal(digest(canonical(basis)), BASIS_BUILD_HASH, "sealed baseline build")
    equal(basis["library"], str(BASELINE), "baseline path")
    equal(basis["library_sha256"], BASELINE_HASH, "baseline hash")


def candidate_build(path):
    cache = h.parse_cache(path / "CMakeCache.txt")
    expected = {"BUILD_SHARED_LIBS": "ON", "CMAKE_BUILD_TYPE": "Release",
        "CMAKE_GENERATOR": "Ninja", "MARCH_NATIVE": "ON", "WH_LTO": "OFF", "WH_PGO_MODE": "OFF"}
    equal({key: cache.get(key) for key in expected}, expected, "candidate build flags")
    equal(cache["CMAKE_HOME_DIRECTORY"], str(ROOT), "candidate source directory")
    dry = command(["/usr/bin/ninja", "-C", str(path), "-n", "wirehair"])
    if dry.splitlines()[-1:] != ["ninja: no work to do."]:
        fail("candidate build is stale")
    return {"directory": str(path), "cache_sha256": file_hash(path / "CMakeCache.txt"),
        "flags": expected, "dry_run": dry,
        "commands": command(["/usr/bin/ninja", "-C", str(path), "-t", "commands", "wirehair"])}


def check_build(receipt, live=True):
    equal(receipt["campaign"], CAMPAIGN, "build campaign")
    validate_basis(receipt["baseline_build"])
    equal(receipt["baseline_path"], str(BASELINE), "build baseline")
    equal(receipt["baseline_sha256"], BASELINE_HASH, "build baseline hash")
    equal(receipt["candidate_sha256"], CANDIDATE_HASH, "frozen arithmetic candidate")
    equal(receipt["libc_sha256"], LIBC_HASH, "build libc hash")
    equal([shlex.split(line) for line in receipt["candidate_build"]["commands"].splitlines()],
          [shlex.split(line) for line in receipt["baseline_build"]["library_commands"].splitlines()],
          "identical production build command roster")
    h.lower_commit(receipt["source_commit"], "build commit")
    for name in ("worker", "candidate", "compiler"):
        h.replay_path(receipt[name + "_path"], "build " + name)
        h.lower_hash(receipt[name + "_sha256"], "build " + name + " hash")
    expected_inputs = {path for path in receipt["baseline_build"]["source_files"] if not path.startswith("bench/")}
    expected_inputs.update(INPUTS + ("abi/wirehair.map",))
    obj(receipt["source_files"], expected_inputs, "source closure")
    for path, sha in receipt["source_files"].items():
        h.lower_hash(sha, "source hash " + path)
        if path in receipt["baseline_build"]["source_files"] and path not in (
                "codec/WirehairV2Profile.cpp", "codec/V2BorrowedSourceTest.cpp"):
            equal(sha, receipt["baseline_build"]["source_files"][path], "unchanged input " + path)
    expected_command = [receipt["compiler_path"], "-std=c++11", "-O3", "-DNDEBUG", "-Wall", "-Wextra", "-Wpedantic", "-Werror",
        '-DWIREHAIR_WH2_SOURCE_GIT_COMMIT="{}"'.format(receipt["source_commit"]), "-I" + str(ROOT / "include")]
    expected_command += [str(ROOT / path) for path in SOURCES] + ["-ldl", "-o", receipt["worker_path"]]
    equal(receipt["command"], expected_command, "standalone worker build command")
    if "libwirehair" in receipt["dynamic"] or "libwirehair" in receipt["ldd"]:
        fail("worker globally links wirehair")
    if live:
        equal(command(["git", "-C", str(ROOT), "rev-parse", "HEAD"]).strip(),
              receipt["source_commit"], "build HEAD")
        equal(source_receipt(), receipt["source_files"], "build input bytes")
        for name in ("worker", "baseline", "candidate", "compiler"):
            equal(file_hash(receipt[name + "_path"]), receipt[name + "_sha256"], name + " bytes")
        equal(file_hash(LIBC), LIBC_HASH, "libc bytes")
        equal(candidate_build(Path(receipt["candidate_build"]["directory"])),
              receipt["candidate_build"], "candidate build stability")


def build(output, candidate):
    output, candidate = Path(output).resolve(), Path(candidate).resolve(strict=True)
    if output.exists():
        fail("build output already exists")
    equal(file_hash(BASIS_PATH), BASIS_HASH, "baseline freeze")
    basis = h.parse_canonical_line(BASIS_PATH.read_bytes(), "baseline freeze")["build"]
    validate_basis(basis)
    compiler = Path("/usr/bin/g++-13").resolve(strict=True)
    receipt = {"campaign": CAMPAIGN, "source_files": source_receipt(),
        "source_commit": command(["git", "-C", str(ROOT), "rev-parse", "HEAD"]).strip(),
        "baseline_build": basis, "baseline_path": str(BASELINE), "baseline_sha256": BASELINE_HASH,
        "candidate_path": str(candidate / "libwirehair.so.2.0.0"),
        "candidate_sha256": file_hash(candidate / "libwirehair.so.2.0.0"),
        "candidate_build": candidate_build(candidate), "compiler_path": str(compiler),
        "compiler_sha256": file_hash(compiler), "compiler_version": command([str(compiler), "--version"]),
        "libc_sha256": file_hash(LIBC)}
    output.mkdir()
    worker = output / "worker"
    argv = [str(compiler), "-std=c++11", "-O3", "-DNDEBUG", "-Wall", "-Wextra", "-Wpedantic", "-Werror",
        '-DWIREHAIR_WH2_SOURCE_GIT_COMMIT="{}"'.format(receipt["source_commit"]),
        "-I" + str(ROOT / "include")] + [str(ROOT / path) for path in SOURCES] + ["-ldl", "-o", str(worker)]
    receipt["command"], receipt["build_output"] = argv, command(argv)
    receipt["worker_path"], receipt["worker_sha256"] = str(worker), file_hash(worker)
    receipt["dynamic"] = command(["/usr/bin/readelf", "-d", str(worker)])
    receipt["ldd"] = command(["/usr/bin/ldd", str(worker)])
    if "libwirehair" in receipt["dynamic"] or "libwirehair" in receipt["ldd"]:
        fail("worker must not globally link either codec")
    equal(command([str(worker), "--describe"]), canonical(description()).decode("ascii") + "\n", "native description")
    check_build(receipt)
    with (output / "build.json").open("xb") as stream:
        stream.write(canonical(receipt) + b"\n")
    return receipt


def statistics(samples):
    stats, failed = {"controls": [], "cells": [], "widths": []}, []
    boundary = math.log1p(0.02)
    cell_values = {}
    for cell, (k, b) in enumerate(CELLS):
        for comparison in range(2):
            for condition in range(4):
                summary = h.sample_summary([samples[rep, cell, comparison, condition] for rep in range(12)])
                stats["controls"].append(dict(cell=cell, comparison=comparison, condition=condition, summary=summary))
                if not (-boundary < summary["lower95_log"] and summary["upper95_log"] < boundary):
                    failed.append("AA:{}:{}:{}".format(cell, comparison, condition))
        values = [math.fsum(samples[rep, cell, comparison, condition] * (1 if condition < 2 else -1)
                           for comparison in (2, 3) for condition in range(4)) / 8 for rep in range(12)]
        cell_values[cell] = values
        summary = h.sample_summary(values)
        stats["cells"].append(dict(K=k, block_bytes=b, summary=summary))
        if summary["upper95_log"] > boundary:
            failed.append("cell-regression:{}:{}".format(k, b))
    for b in (64, 1280):
        values = [math.fsum(cell_values[cell][rep] for cell, (_, width) in enumerate(CELLS) if width == b) / 3
                  for rep in range(12)]
        summary = h.sample_summary(values)
        stats["widths"].append(dict(block_bytes=b, summary=summary))
        if not (summary["upper95_log"] < 0 if b == 64 else summary["upper95_log"] <= 0):
            failed.append("width-gain:" + str(b))
    return stats, failed


def cell_key(rep, cell):
    return digest("offset-r0:rep={}:cell={}".format(rep, cell).encode("ascii"))


def roster():
    for rep in range(12):
        for cell in sorted(range(6), key=lambda cell: cell_key(rep, cell)):
            for index in range(4):
                comparison = (rep + cell + index) % 4
                for condition in CONDITIONS[rep % 4]:
                    yield rep, cell, comparison, condition


def description():
    return {"campaign": CAMPAIGN, "schema": PREFIX + ".describe.v1",
        "cell_order": "sha256:offset-r0:rep=R:cell=I", "comparison_order": "(replicate+cell+index)%4",
        "comparisons": ["B0/B1", "C0/C1", "C0/B0", "C1/B1"],
        "condition_sequences": [list(row) for row in CONDITIONS],
        "conditions": [{"condition": c, "mapping": "direct" if c < 2 else "swapped",
            "order": "ABBA" if c % 2 == 0 else "BAAB",
            "warmup": ["left", "right"] if c % 2 == 0 else ["right", "left"]} for c in range(4)],
        "constructor_orders": [["B0", "C0", "C1", "B1"], ["C0", "B0", "B1", "C1"]],
        "counter_phase": 1024, "expected_encode_calls": 94371840, "expected_measured_batches": 9216,
        "expected_measured_encode_calls": 75497472, "expected_measured_invocations": 3345408,
        "expected_panels": 1152, "expected_records": 1154, "expected_warmup_batches": 2304,
        "expected_warmup_encode_calls": 18874368, "expected_warmup_invocations": 836352,
        "internal_deadline_seconds": 60, "output_phase": 2336, "page_bytes": 4096,
        "panel_replicates": 12, "sibling_cpu": 56, "source_phase": 2048,
        "source_seed_base": h.SOURCE_SEED_BASE_TEXT, "target_cpu": 120,
        "timing_granularity": "whole-batch", "tuples": [{"K": k, "block_bytes": b,
            "scope": "prebuilt-systematic", "scope_invocations_per_batch": 8192 // k} for k, b in CELLS]}


def description_hash():
    return digest(canonical(description()) + b"\n")


def address(value):
    if type(value) is not str or h.re.fullmatch(r"0x[1-9a-f][0-9a-f]{0,15}", value) is None:
        fail("malformed or null address")
    return int(value, 16)


def profile(value, k, b, tail):
    text = value["profile_hex"]
    if type(text) is not str or h.re.fullmatch("[0-9a-f]{64}", text) is None:
        fail("profile hex")
    raw = bytes.fromhex(text)
    equal(raw[:4].hex(), "57485632", "profile magic")
    version, size, identity, message, block, attempt, reserved = struct.unpack("<HHQQIB3s", raw[4:])
    equal((version, size, identity, message, block, reserved.hex()),
          (1, 32, 0x4b295bbb47f4f9c9, (k - 1) * b + tail, b, "000000"), "public descriptor")
    equal(value["profile_sha256"], digest(raw), "profile digest")
    equal(value["source_sha256"], h.SOURCE_SHA256[k, b, tail], "source digest")
    for key in ("repair_sha256", "high_id_sha256"):
        h.lower_hash(value[key], key)


def validate_config(value, receipt):
    obj(value, {"campaign", "schema", "description_sha256", "target_identity", "bindings", "cells"}, "config")
    equal(value["campaign"], CAMPAIGN, "config campaign")
    equal(value["schema"], PREFIX + ".config.v1", "config schema")
    equal(value["description_sha256"], description_hash(), "config description")
    h.validate_target_identity(value["target_identity"], "config target")
    bindings = obj(value["bindings"], {"B", "C"}, "bindings")
    for arm, name in (("B", "baseline"), ("C", "candidate")):
        binding = obj(bindings[arm], {"path", "sha256", "symbols", "memcpy"}, "binding")
        equal(binding["path"], receipt[name + "_path"], "loaded DSO")
        equal(binding["sha256"], receipt[name + "_sha256"], "loaded DSO hash")
        obj(binding["symbols"], {"wirehair_init_", "wirehair_v2_encoder_create_with_options",
            "wirehair_v2_encode", "wirehair_v2_free"}, "symbols")
        for symbol in binding["symbols"].values():
            obj(symbol, {"path", "address", "elf_offset"}, "symbol")
            equal(symbol["path"], binding["path"], "symbol owner")
            if address(symbol["address"]) <= address(symbol["elf_offset"]):
                fail("symbol load base")
        copy = obj(binding["memcpy"], {"path", "address", "elf_offset"}, "memcpy")
        equal(copy["path"], str(LIBC), "copy provider")
        equal(copy["elf_offset"], "0x1a14c0", "copy implementation")
        address(copy["address"])
    equal(bindings["B"]["memcpy"], bindings["C"]["memcpy"], "shared memcpy")
    if bindings["B"]["symbols"]["wirehair_v2_encode"]["address"] == bindings["C"]["symbols"]["wirehair_v2_encode"]["address"]:
        fail("candidate and baseline are the same encode entrypoint")
    if type(value["cells"]) is not list or len(value["cells"]) != 6:
        fail("cell oracle roster")
    keys = {"source_sha256", "profile_hex", "profile_sha256", "repair_sha256", "high_id_sha256"}
    for cell, (k, b) in zip(value["cells"], CELLS):
        obj(cell, keys | {"K", "block_bytes", "partial"}, "cell")
        equal((cell["K"], cell["block_bytes"]), (k, b), "cell shape")
        profile(cell, k, b, b)
        if type(cell["partial"]) is not list or len(cell["partial"]) != 2:
            fail("partial oracle roster")
        for partial, tail in zip(cell["partial"], (1, b - 1)):
            obj(partial, keys | {"tail_bytes", "systematic_sha256"}, "partial")
            equal(partial["tail_bytes"], tail, "tail")
            profile(partial, k, b, tail)
            equal(partial["systematic_sha256"], partial["source_sha256"], "partial systematic bytes")


def validate_addresses(value, k, b):
    obj(value, {"arena", "arena_bytes", "span", "source", "output", "counters", "handles"}, "addresses")
    arena = address(value["arena"])
    equal(arena % 4096, 0, "arena alignment")
    span = ((k * b + 8191) // 4096) * 4096
    equal(value["span"], span, "arena span")
    equal(value["arena_bytes"], 2 * span + 1024 + 8192 * 4 + 64, "arena bytes")
    equal(address(value["source"]), arena + 2048, "shared source")
    equal(address(value["output"]), arena + span + 2336, "shared output")
    equal(address(value["counters"]), arena + 2 * span + 1024, "shared counters")
    if type(value["handles"]) is not list or len(value["handles"]) != 4:
        fail("handle roster")
    handles = [address(item) for item in value["handles"]]
    if len(set(handles)) != 4 or any(arena <= handle < arena + value["arena_bytes"] for handle in handles):
        fail("distinct retained codec handles")


def terminal_value():
    return {"campaign": CAMPAIGN, "schema": PREFIX + ".terminal.v1", "status": "complete",
        "encode_call_count": 94371840, "measured_batch_count": 9216, "measured_invocation_count": 3345408,
        "panel_count": 1152, "record_count": 1154, "warmup_batch_count": 2304, "warmup_invocation_count": 836352}


def parse_transcript(raw, receipt):
    if len(raw) > h.MAX_RAW_BYTES:
        fail("raw byte cap")
    lines = raw.splitlines(keepends=True)
    equal(len(lines), 1154, "record count")
    records = [h.parse_canonical_line(line, "record") for line in lines]
    config = records[0]
    validate_config(config, receipt)
    samples, retained = {}, {}
    for sequence, (panel, key) in enumerate(zip(records[1:-1], roster())):
        rep, cell, comparison, condition = key
        k, b = CELLS[cell]
        order = "ABBA" if condition % 2 == 0 else "BAAB"
        wanted = {"campaign": CAMPAIGN, "schema": PREFIX + ".panel.v1", "description_sha256": description_hash(),
            "sequence": sequence, "replicate": rep, "cell_index": cell, "cell_key_sha256": cell_key(rep, cell),
            "comparison_index": comparison, "condition": condition, "order": order,
            "mapping": "direct" if condition < 2 else "swapped", "K": k, "block_bytes": b,
            "scope_invocations_per_batch": 8192 // k, "source_sha256": h.SOURCE_SHA256[k, b, b],
            "output_sha256": h.SOURCE_SHA256[k, b, b], "profile_sha256": config["cells"][cell]["profile_sha256"],
            "source_immutable_verified": True}
        obj(panel, set(wanted) | {"addresses", "slots"}, "panel")
        for name, expected in wanted.items():
            equal(panel[name], expected, "panel " + name)
        validate_addresses(panel["addresses"], k, b)
        if (rep, cell) in retained:
            equal(panel["addresses"], retained[rep, cell], "retained exact I/O and handles")
        retained[rep, cell] = panel["addresses"]
        slots = panel["slots"]
        if type(slots) is not list or len(slots) != 8:
            fail("slot roster")
        elapsed = []
        for slot, side in zip(slots, h.sides_for(order)):
            obj(slot, {"elapsed_ns", "logical_lane", "physical_lane"}, "slot")
            equal(slot["logical_lane"], side, "logical lane")
            equal(slot["physical_lane"], int(side == "right") ^ (condition // 2), "physical lane")
            elapsed.append(h.exact_int(slot["elapsed_ns"], 1, h.MAX_INT63, "elapsed"))
        samples[key] = h.lane_contrast(elapsed, order)
    equal(records[-1], terminal_value(), "terminal")
    return statistics(samples)


def marker_value(claim):
    return {"campaign": CAMPAIGN, "schema": PREFIX + ".worker-started.v1",
        "claim_sha256": digest(canonical(claim) + b"\n"), "worker_sha256": claim["worker_sha256"],
        "description_sha256": description_hash()}


def adjudicate(raw, stderr, p, marker, claim, elapsed):
    obj(p, {"build", "controller_affinity", "interpreter", "interpreter_after", "git_before", "git_after",
        "preserved_before", "preserved_after", "health_before", "health_after", "errors", "process", "build_after"}, "provenance")
    obj(claim, {"baseline_path", "baseline_sha256", "candidate_path", "candidate_sha256", "worker_sha256", "source_commit",
        "campaign", "schema", "created_unix_ns", "description_sha256", "build_receipt_sha256"}, "claim")
    check_build(p["build"], live=False)
    equal(claim["build_receipt_sha256"], digest(canonical(p["build"])), "claimed build")
    for key in ("baseline_path", "baseline_sha256", "candidate_path", "candidate_sha256", "worker_sha256", "source_commit"):
        equal(claim[key], p["build"][key], "claim " + key)
    equal(claim["campaign"], CAMPAIGN, "claim campaign")
    equal(claim["schema"], PREFIX + ".claim.v1", "claim schema")
    equal(claim["description_sha256"], description_hash(), "claimed description")
    h.exact_int(claim["created_unix_ns"], 1, h.MAX_INT63, "claim timestamp")
    h.validate_interpreter_receipt(p["interpreter"], "interpreter")
    h.validate_git_receipt(p["git_before"], claim["source_commit"], "git before")
    affinity = obj(p["controller_affinity"], {"affinity_before", "affinity_after", "controller_cpu",
        "sibling_available_before", "singleton_verified", "target_available_before"}, "controller affinity")
    inherited = affinity["affinity_before"]
    if type(inherited) is not list or any(type(x) is not int or not 0 <= x <= 1048575 for x in inherited) or inherited != sorted(set(inherited)) or not {0, 56, 120}.issubset(inherited):
        fail("inherited affinity")
    equal(affinity, dict(affinity_before=inherited, affinity_after=[0], controller_cpu=0,
        sibling_available_before=True, singleton_verified=True, target_available_before=True), "controller affinity")
    if type(elapsed) not in (int, float) or not math.isfinite(elapsed) or elapsed < 0:
        fail("elapsed is invalid")
    if type(p["errors"]) is not list or len(p["errors"]) > 16 or any(type(error) is not str or len(error) > 65536 for error in p["errors"]):
        fail("error receipt")
    failures = list(p["errors"])
    process = obj(p["process"], {"exit_code", "elapsed_seconds", "timed_out"}, "process")
    if process["exit_code"] is not None:
        h.exact_int(process["exit_code"], -255, 255, "process exit")
    if type(process["timed_out"]) is not bool or type(process["elapsed_seconds"]) not in (int, float) or not math.isfinite(process["elapsed_seconds"]) or not 0 <= process["elapsed_seconds"] <= elapsed:
        fail("process time or timeout")
    if elapsed >= 70:
        failures.append("outer-deadline")
    if p["process"]["exit_code"] != 0 or p["process"]["timed_out"] or stderr:
        failures.append("worker-failure")
    try:
        equal(marker, marker_value(claim), "worker marker")
        equal(p["build_after"], None, "successful postflight build check")
        equal(p["interpreter_after"], p["interpreter"], "interpreter stability")
        equal(p["git_after"], p["git_before"], "git stability")
        h.validate_r0_receipt(p["preserved_before"], "old campaigns before")
        h.validate_r0_receipt(p["preserved_after"], "old campaigns after")
        equal(p["preserved_before"], p["preserved_after"], "old campaigns")
        h.validate_health_pair(p["health_before"], p["health_after"])
    except (h.ValidationError, KeyError, TypeError, ValueError) as error:
        failures.append("postflight:" + str(error))
    stats, failed = None, []
    try:
        stats, failed = parse_transcript(raw, p["build"])
    except (h.ValidationError, KeyError, TypeError, ValueError, IndexError, OverflowError) as error:
        failures.append("transcript:" + str(error))
    return {"campaign": CAMPAIGN, "schema": PREFIX + ".summary.v1", "raw_sha256": digest(raw),
        "outcome": "invalid" if failures else "fail" if failed else "pass",
        "infrastructure_failures": failures, "failed_gates": failed, "statistics": stats,
        "elapsed_seconds": elapsed, "WH1_compared": False, "promotion_claimed": False}


def run_once(build_path):
    receipt = h.parse_canonical_line(Path(build_path).read_bytes(), "build receipt")
    check_build(receipt)
    git = h.git_receipt(ROOT, receipt["source_commit"])
    equal(command(["git", "-C", str(ROOT), "rev-parse", "@{upstream}"]).strip(), receipt["source_commit"], "pushed source")
    p = {"build": receipt, "controller_affinity": h.pin_controller(), "interpreter": h.interpreter_receipt(),
         "preserved_before": h.snapshot_r0(), "health_before": h.cpu_receipt(), "git_before": git, "errors": []}
    h.validate_interpreter_receipt(p["interpreter"], "interpreter")
    h.validate_health_receipt(p["health_before"], "health before")
    claim = {key: receipt[key] for key in ("baseline_path", "baseline_sha256", "candidate_path", "candidate_sha256", "worker_sha256", "source_commit")}
    claim.update(campaign=CAMPAIGN, schema=PREFIX + ".claim.v1", created_unix_ns=time.time_ns(),
                 description_sha256=description_hash(), build_receipt_sha256=digest(canonical(receipt)))
    directory = h.claim_namespace(claim)
    started, raw, stderr, marker = time.monotonic(), b"", b"", None
    try:
        p["process"] = {"exit_code": None, "elapsed_seconds": 0, "timed_out": False}
        try:
            raw, stderr, code, elapsed, timed_out = h.run_worker(receipt["worker_path"], receipt["candidate_path"],
                started + 60, digest(canonical(claim) + b"\n"), claim["worker_sha256"], description_hash())
            p["process"] = {"exit_code": code, "elapsed_seconds": elapsed, "timed_out": timed_out}
        except Exception as error:
            p["errors"].append("execution:" + str(error))
        for key in ("health_after", "preserved_after", "build_after", "git_after", "interpreter_after"):
            p[key] = None
        for key, probe in (("health_after", h.cpu_receipt), ("preserved_after", h.snapshot_r0),
                           ("build_after", lambda: check_build(receipt)),
                           ("git_after", lambda: h.git_receipt(ROOT, receipt["source_commit"])),
                           ("interpreter_after", h.interpreter_receipt)):
            try:
                if time.monotonic() - started >= 70:
                    fail("outer deadline")
                p[key] = probe()
            except Exception as error:
                p["errors"].append(key + ":" + str(error))
        marker_path = OUTPUT / "WORKER_STARTED"
        if marker_path.exists():
            try:
                marker = h.parse_canonical_line(marker_path.read_bytes(), "worker marker")
            except Exception as error:
                p["errors"].append("marker:" + str(error))
        h.write_exclusive(directory, "raw.jsonl", raw, 0o400)
        h.write_exclusive(directory, "worker.stderr", stderr, 0o400)
        h.publish_json(directory, "provenance.json", p)
        summary = adjudicate(raw, stderr, p, marker, claim, time.monotonic() - started)
        summary["elapsed_seconds"] = time.monotonic() - started
        if summary["elapsed_seconds"] >= 70 and "outer-deadline" not in summary["infrastructure_failures"]:
            summary["infrastructure_failures"].append("outer-deadline")
            summary["outcome"] = "invalid"
        h.publish_json(directory, "summary.json", summary)
        names = FILES | ({"WORKER_STARTED"} if marker_path.exists() else set())
        hashes = {name: file_hash(OUTPUT / name) for name in sorted(names)}
        if time.monotonic() - started >= 70:
            fail("deadline before COMPLETE; namespace remains spent")
        h.publish_json(directory, "COMPLETE", hashes)
        os.fsync(directory)
        return summary
    finally:
        os.close(directory)


def replay(path):
    path = h.exact_absolute_directory(str(path), "bundle")
    names = {p.name for p in path.iterdir()}
    if names not in (FILES | {"COMPLETE"}, FILES | {"COMPLETE", "WORKER_STARTED"}):
        fail("bundle roster")
    data = {}
    for name in names:
        member = h.exact_absolute_file(str(path / name), name)
        if member.stat().st_mode & 0o777 != 0o400 or member.stat().st_nlink != 1:
            fail("bundle metadata")
        cap = h.MAX_RAW_BYTES if name == "raw.jsonl" else h.MAX_PROBE_BYTES
        with member.open("rb") as stream:
            data[name] = stream.read(cap + 1)
        if len(data[name]) > cap:
            fail("bundle member too large")
    objects = {name: h.parse_canonical_line(value, name) for name, value in data.items()
               if name not in ("raw.jsonl", "worker.stderr")}
    equal(objects["COMPLETE"], {name: digest(value) for name, value in data.items() if name != "COMPLETE"}, "bundle hashes")
    summary = adjudicate(data["raw.jsonl"], data["worker.stderr"], objects["provenance.json"],
        objects.get("WORKER_STARTED"), objects["CLAIM"], objects["summary.json"]["elapsed_seconds"])
    equal(objects["summary.json"], summary, "recomputed summary")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--build")
    mode.add_argument("--run", action="store_true")
    mode.add_argument("--replay")
    parser.add_argument("--candidate-build")
    parser.add_argument("--build-receipt")
    args = parser.parse_args()
    if args.build:
        if not args.candidate_build:
            fail("build needs candidate-build")
        result = build(args.build, args.candidate_build)
    elif args.run:
        if not args.build_receipt:
            fail("run needs build-receipt")
        result = run_once(args.build_receipt)
    else:
        result = replay(args.replay)
    print(canonical(result).decode("ascii"))
    return {"pass": 0, "fail": 2, "invalid": 1}.get(result.get("outcome"), 0)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (h.ValidationError, OSError, ValueError, TypeError, KeyError, IndexError) as error:
        print("Offset r0: " + str(error), file=sys.stderr)
        sys.exit(1)
