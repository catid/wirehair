#!/usr/bin/env python3
"""One frozen baseline-only stack-placement crossover; never a WH1 speed claim.

The old Offset module is imported privately for immutable byte/receipt helpers,
not its candidate build, launcher, old workload or decision rule.
"""
import argparse
import importlib.util
import itertools
import math
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location("_stack_crossover_offset",
    ROOT / "bench/Wh2PublicOffsetR0.py")
o = importlib.util.module_from_spec(spec)
spec.loader.exec_module(o)
h = o.h
canonical, digest, fail = o.canonical, o.digest, o.fail
equal, obj, file_hash, command = o.equal, o.obj, o.file_hash, o.command
address, profile, validate_addresses = o.address, o.profile, o.validate_addresses
BASELINE, BASELINE_HASH = o.BASELINE, o.BASELINE_HASH
BASIS_PATH, BASIS_HASH, BASIS_BUILD_HASH = o.BASIS_PATH, o.BASIS_HASH, o.BASIS_BUILD_HASH
LIBC, LIBC_HASH = o.LIBC, o.LIBC_HASH
COMPILER_HASH = "1353e9bdd29a7295c7226bf6c63abccce056d8cac31f112e5cdbecc3f28c2769"
HELPER_HASHES = {
    "bench/Wh2PublicOffsetR0.cpp": "da7d33bebfba6865f8fd9b621e544c5c8dcbef1c67a4b4b9b92b2ca772cc15af",
    "bench/Wh2PublicOffsetR0.py": "4c26db58391a7fe434c05d2f4098c978ed5e62661413e3a0040792a8c09c2715",
    "bench/Wh2PublicBorrowedCurrentVsWh1R1.py": "51aa47781355f20a3024438771c84cc36c16b2afb39037bdd4c7bd0459f17ac7",
}
CAMPAIGN = "wh2-public-stack-crossover-r0"
PREFIX = "wirehair.wh2.public-stack-crossover-r0"
OUTPUT = Path("/var/tmp/" + CAMPAIGN)
CELLS, CONDITIONS = o.CELLS, o.CONDITIONS
TARGET_ORDERS = tuple(itertools.permutations(range(3)))
SOURCES = ("bench/Wh2PublicStackCrossoverR0.cpp",) + o.SOURCES[1:]
INPUTS = SOURCES + ("bench/Wh2PublicStackCrossoverR0.py",
    "bench/test_Wh2PublicStackCrossoverR0.py", "bench/Wh2PublicOffsetR0.cpp",
    "bench/Wh2PublicOffsetR0.py", "bench/test_Wh2PublicOffsetR0.py",
    "bench/test_Wh2PublicBorrowedCurrentVsWh1R1.py") + o.INPUTS[6:]
h.CAMPAIGN, h.FIXED_OUTPUT_DIR, h.R1_TRACKED_INPUTS = CAMPAIGN, OUTPUT, INPUTS
h.R0_ROOTS += tuple(Path("/var/tmp/" + name) for name in (
    "wh2-public-offset-r0", "wh2-production-mix3-k3k6-portability-r0",
    "wh2-production-mix3-k3k6-joint-r0", "wh2-production-mix3-k3k6-width-local-r0"))
FILES = o.FILES


def source_receipt():
    files = command(["git", "-C", str(ROOT), "ls-files", "--", "CMakeLists.txt",
        "cmake", "codec", "include", "tables", "abi/wirehair.map", "wirehair.cpp",
        "Wirehair*", "gf256.cpp", "gf256.h"]).splitlines()
    return {path: file_hash(ROOT / path) for path in sorted(set(files) | set(INPUTS))}


def build_command(receipt):
    return [receipt["compiler_path"], "-std=c++11", "-O3", "-DNDEBUG", "-Wall",
        "-Wextra", "-Wpedantic", "-Werror",
        '-DWIREHAIR_WH2_SOURCE_GIT_COMMIT="{}"'.format(receipt["source_commit"]),
        "-I" + str(ROOT / "include")] + [str(ROOT / path) for path in SOURCES] + [
        "-ldl", "-o", receipt["worker_path"]]


def codegen_receipt(worker):
    # Read-only inspection: never objcopy or in-place object manipulation.
    return {name: digest(command(argv).encode("utf-8")) for name, argv in (
        ("disassembly_sha256", ["/usr/bin/objdump", "-d", "-C", worker]),
        ("unwind_sha256", ["/usr/bin/readelf", "--debug-dump=frames", worker]))}


def check_build(receipt, live=True):
    obj(receipt, {"campaign", "source_files", "source_commit", "baseline_build",
        "baseline_path", "baseline_sha256", "compiler_path", "compiler_sha256",
        "compiler_version", "libc_sha256", "worker_path", "worker_sha256",
        "command", "build_output", "dynamic", "ldd", "codegen"}, "build")
    equal(receipt["campaign"], CAMPAIGN, "build campaign")
    o.validate_basis(receipt["baseline_build"])
    equal(receipt["baseline_path"], str(BASELINE), "build baseline")
    equal(receipt["baseline_sha256"], BASELINE_HASH, "baseline hash")
    equal(receipt["compiler_sha256"], COMPILER_HASH, "compiler hash")
    equal(receipt["compiler_path"], "/usr/bin/x86_64-linux-gnu-g++-13", "compiler")
    equal(receipt["libc_sha256"], LIBC_HASH, "libc hash")
    h.lower_commit(receipt["source_commit"], "build commit")
    h.replay_path(receipt["worker_path"], "worker")
    h.lower_hash(receipt["worker_sha256"], "worker hash")
    expected_inputs = {path for path in receipt["baseline_build"]["source_files"]
                       if not path.startswith("bench/")}
    expected_inputs.update(INPUTS + ("abi/wirehair.map",))
    obj(receipt["source_files"], expected_inputs, "source closure")
    for path, sha in HELPER_HASHES.items():
        equal(receipt["source_files"][path], sha, "frozen helper " + path)
    for path, sha in receipt["source_files"].items():
        h.lower_hash(sha, "source hash " + path)
        if path in receipt["baseline_build"]["source_files"] and path != "codec/V2BorrowedSourceTest.cpp":
            equal(sha, receipt["baseline_build"]["source_files"][path], "unchanged input " + path)
    equal(receipt["command"], build_command(receipt), "standalone build command")
    equal(receipt["build_output"], "", "quiet strict build")
    obj(receipt["codegen"], {"disassembly_sha256", "unwind_sha256"}, "codegen")
    for sha in receipt["codegen"].values():
        h.lower_hash(sha, "codegen hash")
    for key in ("compiler_version", "dynamic", "ldd"):
        if type(receipt[key]) is not str or not receipt[key]:
            fail("missing build probe " + key)
    if "libwirehair" in receipt["dynamic"] or "libwirehair" in receipt["ldd"]:
        fail("worker globally links wirehair")
    if live:
        equal(command(["git", "-C", str(ROOT), "rev-parse", "HEAD"]).strip(),
              receipt["source_commit"], "build HEAD")
        equal(source_receipt(), receipt["source_files"], "build input bytes")
        for name in ("worker", "baseline", "compiler"):
            equal(file_hash(receipt[name + "_path"]), receipt[name + "_sha256"], name + " bytes")
        equal(file_hash(LIBC), LIBC_HASH, "libc bytes")
        equal(codegen_receipt(receipt["worker_path"]), receipt["codegen"], "codegen stability")


def build(output):
    output = Path(output).resolve()
    if output.exists():
        fail("build output already exists")
    equal(file_hash(BASIS_PATH), BASIS_HASH, "baseline freeze")
    basis = read_json(BASIS_PATH, "baseline freeze")["build"]
    o.validate_basis(basis)
    compiler = Path("/usr/bin/g++-13").resolve(strict=True)
    receipt = {"campaign": CAMPAIGN, "source_files": source_receipt(),
        "source_commit": command(["git", "-C", str(ROOT), "rev-parse", "HEAD"]).strip(),
        "baseline_build": basis, "baseline_path": str(BASELINE), "baseline_sha256": BASELINE_HASH,
        "compiler_path": str(compiler), "compiler_sha256": file_hash(compiler),
        "compiler_version": command([str(compiler), "--version"]), "libc_sha256": file_hash(LIBC),
        "worker_path": str(output / "worker")}
    output.mkdir()
    receipt["command"] = build_command(receipt)
    receipt["build_output"] = command(receipt["command"])
    receipt["worker_sha256"] = file_hash(receipt["worker_path"])
    receipt["dynamic"] = command(["/usr/bin/readelf", "-d", receipt["worker_path"]])
    receipt["ldd"] = command(["/usr/bin/ldd", receipt["worker_path"]])
    receipt["codegen"] = codegen_receipt(receipt["worker_path"])
    equal(command([receipt["worker_path"], "--describe"]),
          canonical(description()).decode("ascii") + "\n", "native description")
    check_build(receipt)
    with (output / "build.json").open("xb") as stream:
        stream.write(canonical(receipt) + b"\n")
    return receipt


def statistics(samples):
    if set(samples) != set(roster()):
        fail("complete statistics roster differs")
    if any(type(value) not in (int, float) or not math.isfinite(value) for value in samples.values()):
        fail("invalid contrast sample")
    stats, failed = {"controls": [], "deltas": []}, []
    boundary = math.log1p(0.02)
    supported, bounded = True, True
    for cell in range(6):
        for target in range(3):
            for comparison in range(2):
                for condition in range(4):
                    summary = h.sample_summary([samples[rep, cell, target, comparison, condition]
                                                for rep in range(12)])
                    stats["controls"].append(dict(cell=cell, target=target,
                        comparison=comparison, condition=condition, summary=summary))
                    if not (-boundary < summary["lower95_log"] and summary["upper95_log"] < boundary):
                        failed.append("AA:{}:{}:{}:{}".format(cell, target, comparison, condition))
        for target in range(2):
            for condition in range(4):
                sign = 1 if condition < 2 else -1
                values = [sign * (samples[rep, cell, target, 2, condition] -
                                  samples[rep, cell, 2, 2, condition]) for rep in range(12)]
                summary = h.sample_summary(values)
                primary = cell in (3, 5)
                stats["deltas"].append(dict(cell=cell, target=target,
                    condition=condition, primary=primary, summary=summary))
                if primary:
                    supported &= summary["lower95_log"] > 0 if target == 0 else summary["upper95_log"] < 0
                    bounded &= -boundary < summary["lower95_log"] and summary["upper95_log"] < boundary
    stats["mechanistic_outcome"] = "SUPPORTED" if supported else "WITHIN_2PCT" if bounded else "INCONCLUSIVE"
    return stats, failed


def cell_key(rep, cell):
    return digest("stack-crossover-r0:rep={}:cell={}".format(rep, cell).encode("ascii"))


def roster():
    for rep in range(12):
        for cell in sorted(range(6), key=lambda cell: cell_key(rep, cell)):
            for target in TARGET_ORDERS[(rep + cell) % 6]:
                for index in range(3):
                    comparison = (rep + cell + index) % 3
                    for condition in CONDITIONS[rep % 4]:
                        yield rep, cell, target, comparison, condition


def description():
    value = o.description()
    value.update(campaign=CAMPAIGN, schema=PREFIX + ".describe.v1",
        cell_order="sha256:stack-crossover-r0:rep=R:cell=I",
        comparison_order="(replicate+cell+index)%3",
        comparisons=["h0/h0", "h1/h1", "h0/h1"],
        constructor_order=[0, 1, 2, 3],
        target_orders=[list(row) for row in TARGET_ORDERS],
        target_order="permutations[(replicate+cell)%6]",
        expected_encode_calls=212336640, expected_measured_batches=20736,
        expected_measured_encode_calls=169869312, expected_measured_invocations=7527168,
        expected_panels=2592, expected_records=2594, expected_warmup_batches=5184,
        expected_warmup_encode_calls=42467328, expected_warmup_invocations=1881792,
        fixed_envelope_bytes=256, moving_store_bytes=40, stack_limit_bytes=16384)
    del value["constructor_orders"]
    return value


def description_hash():
    return digest(canonical(description()) + b"\n")


def validate_config(value, receipt):
    obj(value, {"campaign", "schema", "description_sha256", "target_identity", "bindings", "cells"}, "config")
    equal(value["campaign"], CAMPAIGN, "config campaign")
    equal(value["schema"], PREFIX + ".config.v1", "config schema")
    equal(value["description_sha256"], description_hash(), "config description")
    h.validate_target_identity(value["target_identity"], "config target")
    bindings = obj(value["bindings"], {"baseline"}, "bindings")
    binding = obj(bindings["baseline"], {"path", "sha256", "symbols", "memcpy"}, "binding")
    equal(binding["path"], receipt["baseline_path"], "loaded DSO")
    equal(binding["sha256"], receipt["baseline_sha256"], "loaded DSO hash")
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


def modular_bytes(start, size):
    return {(start + offset) & 4095 for offset in range(size)}


def select_geometry(addresses, probe):
    obj(probe, {"pre_rsp", "hot_rsp", "final_rsp", "frame_address", "correction"}, "probe")
    pre, frame = address(probe["pre_rsp"]), address(probe["frame_address"])
    if pre % 16 or frame % 16 or pre > (1 << 64) - 256 or not pre <= frame <= pre + 240:
        fail("probe fixed frame")
    equal(probe["hot_rsp"], probe["pre_rsp"], "probe hot")
    equal(probe["final_rsp"], probe["pre_rsp"], "probe final")
    equal(probe["correction"], 0, "probe correction")
    arena = address(addresses["arena"])
    size = h.exact_int(addresses["arena_bytes"], 1, h.MAX_INT63, "arena bytes")
    if arena > (1 << 64) - 1 - size:
        fail("arena address overflows")
    if type(addresses["handles"]) is not list or len(addresses["handles"]) != 4:
        fail("geometry handles")
    handles = [address(item) for item in addresses["handles"]]
    if len(set(handles)) != 4 or any(item % 16 or item > (1 << 64) - 1 - 0x110 or
                                  arena <= item < arena + size for item in handles):
        fail("geometry handle identity")
    phases = [(item + 0x110) & 4095 for item in handles]
    windows = [modular_bytes(phase - 40, 40) for phase in phases]
    fixed = modular_bytes(pre, 256)
    for i in range(4):
        for j in range(i + 1, 4):
            if windows[i] & windows[j] or windows[i] & fixed or windows[j] & fixed:
                continue
            for null in range(0, 4096, 16):
                if not modular_bytes(null - 40, 40) & (windows[i] | windows[j] | fixed):
                    return dict(selected_pair=[i, j], target_phases=[phases[i], phases[j], null],
                                probe=probe, fixed_envelope_bytes=256)
            fail("no null geometry")
    fail("no handle pair geometry")


def validate_stack(value, probe, target):
    obj(value, {"pre_rsp", "hot_rsp", "final_rsp", "frame_address", "correction"}, "stack")
    equal(value["pre_rsp"], probe["pre_rsp"], "stable pre RSP")
    equal(value["frame_address"], probe["frame_address"], "stable frame")
    pre, hot = address(value["pre_rsp"]), address(value["hot_rsp"])
    correction = 4096 + ((pre - target) & 4095)
    equal(value["correction"], correction, "stack correction")
    if not 4096 <= correction <= 8176 or hot != pre - correction or hot % 4096 != target:
        fail("actual stack placement")
    equal(value["final_rsp"], value["hot_rsp"], "final RSP")


def terminal_value():
    return {"campaign": CAMPAIGN, "schema": PREFIX + ".terminal.v1", "status": "complete",
        "encode_call_count": 212336640, "measured_batch_count": 20736, "measured_invocation_count": 7527168,
        "panel_count": 2592, "record_count": 2594, "warmup_batch_count": 5184, "warmup_invocation_count": 1881792}


def parse_transcript(raw, receipt):
    if len(raw) > h.MAX_RAW_BYTES:
        fail("raw byte cap")
    lines = raw.splitlines(keepends=True)
    equal(len(lines), 2594, "record count")
    records = [h.parse_canonical_line(line, "record") for line in lines]
    config = records[0]
    validate_config(config, receipt)
    samples, retained = {}, {}
    for sequence, (panel, key) in enumerate(zip(records[1:-1], roster())):
        rep, cell, target, comparison, condition = key
        k, b = CELLS[cell]
        order = "ABBA" if condition % 2 == 0 else "BAAB"
        wanted = {"campaign": CAMPAIGN, "schema": PREFIX + ".panel.v1", "description_sha256": description_hash(),
            "sequence": sequence, "replicate": rep, "cell_index": cell, "cell_key_sha256": cell_key(rep, cell),
            "target_index": target, "comparison_index": comparison, "condition": condition, "order": order,
            "mapping": "direct" if condition < 2 else "swapped", "K": k, "block_bytes": b,
            "scope_invocations_per_batch": 8192 // k, "source_sha256": h.SOURCE_SHA256[k, b, b],
            "output_sha256": h.SOURCE_SHA256[k, b, b], "profile_sha256": config["cells"][cell]["profile_sha256"],
            "source_immutable_verified": True}
        obj(panel, set(wanted) | {"addresses", "slots", "warmup", "geometry", "target_phase"}, "panel")
        for name, expected in wanted.items():
            equal(panel[name], expected, "panel " + name)
        validate_addresses(panel["addresses"], k, b)
        geometry = select_geometry(panel["addresses"], panel["geometry"]["probe"])
        equal(panel["geometry"], geometry, "first eligible geometry")
        pair = geometry["selected_pair"]
        phase = geometry["target_phases"][target]
        equal(panel["target_phase"], phase, "target phase")
        stable = (panel["addresses"], geometry)
        if (rep, cell) in retained:
            equal(stable, retained[rep, cell], "retained geometry")
        retained[rep, cell] = stable
        selected = (pair[0], pair[0]) if comparison == 0 else (pair[1], pair[1]) if comparison == 1 else pair
        elapsed = []
        for name, sides in (("warmup", ["left", "right"] if condition % 2 == 0 else ["right", "left"]),
                            ("slots", h.sides_for(order))):
            if type(panel[name]) is not list or len(panel[name]) != len(sides):
                fail(name + " roster")
            for slot, side in zip(panel[name], sides):
                obj(slot, {"elapsed_ns", "logical_lane", "physical_lane", "handle_index", "stack"}, name)
                physical = int(side == "right") ^ (condition // 2)
                equal(slot["logical_lane"], side, "logical lane")
                equal(slot["physical_lane"], physical, "physical lane")
                equal(slot["handle_index"], selected[physical], "physical handle")
                validate_stack(slot["stack"], geometry["probe"], phase)
                if name == "warmup":
                    equal(slot["elapsed_ns"], 0, "untimed warmup")
                else:
                    elapsed.append(h.exact_int(slot["elapsed_ns"], 1, h.MAX_INT63, "elapsed"))
        samples[key] = h.lane_contrast(elapsed, order)
    equal(records[-1], terminal_value(), "terminal")
    return statistics(samples)


def read_json(path, where):
    with Path(path).open("rb") as stream:
        data = stream.read(h.MAX_LINE_BYTES + 1)
    return h.parse_canonical_line(data, where)

def marker_value(claim):
    return {"campaign": CAMPAIGN, "schema": PREFIX + ".worker-started.v1",
        "claim_sha256": digest(canonical(claim) + b"\n"), "worker_sha256": claim["worker_sha256"],
        "description_sha256": description_hash()}


def adjudicate(raw, stderr, p, marker, claim, elapsed):
    obj(p, {"build", "controller_affinity", "interpreter", "interpreter_after", "git_before", "git_after",
        "preserved_before", "preserved_after", "health_before", "health_after", "errors", "process", "build_after"}, "provenance")
    obj(claim, {"baseline_path", "baseline_sha256", "worker_sha256", "source_commit",
        "campaign", "schema", "created_unix_ns", "description_sha256", "build_receipt_sha256"}, "claim")
    check_build(p["build"], live=False)
    equal(claim["build_receipt_sha256"], digest(canonical(p["build"])), "claimed build")
    for key in ("baseline_path", "baseline_sha256", "worker_sha256", "source_commit"):
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
    if process["elapsed_seconds"] >= 60:
        failures.append("worker-deadline")
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
        "outcome": "INVALID" if failures else "CONTROL_FAIL" if failed else stats["mechanistic_outcome"],
        "infrastructure_failures": failures, "failed_gates": failed, "statistics": stats,
        "elapsed_seconds": elapsed, "WH1_compared": False, "promotion_claimed": False}


def run_once(build_path):
    receipt = read_json(build_path, "build receipt")
    check_build(receipt)
    git = h.git_receipt(ROOT, receipt["source_commit"])
    equal(command(["git", "-C", str(ROOT), "rev-parse", "@{upstream}"]).strip(), receipt["source_commit"], "pushed source")
    p = {"build": receipt, "controller_affinity": h.pin_controller(), "interpreter": h.interpreter_receipt(),
         "preserved_before": h.snapshot_r0(), "health_before": h.cpu_receipt(), "git_before": git, "errors": []}
    h.validate_interpreter_receipt(p["interpreter"], "interpreter")
    h.validate_health_receipt(p["health_before"], "health before")
    claim = {key: receipt[key] for key in ("baseline_path", "baseline_sha256", "worker_sha256", "source_commit")}
    claim.update(campaign=CAMPAIGN, schema=PREFIX + ".claim.v1", created_unix_ns=time.time_ns(),
                 description_sha256=description_hash(), build_receipt_sha256=digest(canonical(receipt)))
    directory = h.claim_namespace(claim)
    started, raw, stderr, marker = time.monotonic(), b"", b"", None
    try:
        p["process"] = {"exit_code": None, "elapsed_seconds": 0, "timed_out": False}
        try:
            raw, stderr, code, elapsed, timed_out = h.run_worker(receipt["worker_path"], receipt["baseline_path"],
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
                marker = read_json(marker_path, "worker marker")
            except Exception as error:
                p["errors"].append("marker:" + str(error))
        h.write_exclusive(directory, "raw.jsonl", raw, 0o400)
        h.write_exclusive(directory, "worker.stderr", stderr, 0o400)
        h.publish_json(directory, "provenance.json", p)
        summary = adjudicate(raw, stderr, p, marker, claim, time.monotonic() - started)
        summary["elapsed_seconds"] = time.monotonic() - started
        if summary["elapsed_seconds"] >= 70 and "outer-deadline" not in summary["infrastructure_failures"]:
            summary["infrastructure_failures"].append("outer-deadline")
            summary["outcome"] = "INVALID"
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
        cap = h.MAX_RAW_BYTES if name == "raw.jsonl" else h.MAX_STDERR_BYTES if name == "worker.stderr" else h.MAX_PROBE_BYTES
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
    parser.add_argument("--build-receipt")
    args = parser.parse_args()
    if args.build:
        if args.build_receipt:
            fail("build does not accept build-receipt")
        result = build(args.build)
    elif args.run:
        if not args.build_receipt:
            fail("run needs build-receipt")
        result = run_once(args.build_receipt)
    else:
        if args.build_receipt:
            fail("replay does not accept build-receipt")
        result = replay(args.replay)
    print(canonical(result).decode("ascii"), flush=True)
    return {"INVALID": 1, "CONTROL_FAIL": 2}.get(result.get("outcome"), 0)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (h.ValidationError, OSError, ValueError, TypeError, KeyError, IndexError) as error:
        print("Stack crossover r0: " + str(error), file=sys.stderr)
        sys.exit(1)
