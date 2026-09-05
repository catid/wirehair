#!/usr/bin/env python3
"""Bounded page-layout A/A experiment, not a WH2/WH1 performance comparison.

Reuse the reviewed r0 supervisor, health/build probes and immutable publication
protocol in a private module instance. Only the new roster, address checks and
statistics differ. No old campaign or helper source is modified.
"""
import argparse
import importlib.util
import math
from pathlib import Path
import sys

_spec = importlib.util.spec_from_file_location(
    "_page_aa_base", Path(__file__).with_name("Wh2PublicAaCalibrationR0.py"))
a = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(a)
h = a.h
equal, obj, fail, canonical, digest = a.equal, a.obj, a.fail, a.canonical, a.digest
artifacts, build_receipt, validate_build = a.artifacts, a.build_receipt, a.validate_build
environment = a.environment
CAMPAIGN = "wh2-public-page-aa-r0"
PREFIX = "wirehair.wh2.public-page-aa-r0"
TARGET = "wirehair_wh2_public_page_aa_r0"
FIXED_OUTPUT_DIR = Path("/var/tmp/" + CAMPAIGN)
CELLS = ((8, 64), (8, 1280), (128, 1280))
TUPLES = tuple((k, b, "prebuilt-systematic") for k, b in CELLS)
ROLES = ("C", "L", "M")
CONDITION_SEQUENCES = a.CONDITION_SEQUENCES
PHASES = ((0x920, 0x920), (0x820, 0x820), (0x920, 0x820), (0x920, 0x920))
COPY_LIBRARY = "/usr/lib/x86_64-linux-gnu/libc.so.6"
COPY_LIBRARY_SHA256 = "8db37cf3f2169f59a0f07ef1fea308c35656668c64c8ff294e1860f4121eb161"
SOURCES = ("bench/Wh2PublicPageAaR0.cpp",) + a.SOURCES[1:]
INPUTS = a.INPUTS + ("bench/Wh2PublicPageAaR0.cpp", "bench/Wh2PublicPageAaR0.py",
                    "bench/test_Wh2PublicPageAaR0.py")
COUNTS = {
    "expected_panels": 1728, "expected_records": 1730,
    "expected_measured_batches": 13824, "expected_warmup_batches": 3456,
    "expected_measured_invocations": 9732096, "expected_warmup_invocations": 2433024,
    "expected_measured_encode_calls": 113246208,
    "expected_warmup_encode_calls": 28311552, "expected_encode_calls": 141557760,
}
# The module is private: rebinding its protocol callbacks cannot affect any
# old controller, evidence bundle or concurrently imported test instance.
for _name in ("CAMPAIGN", "PREFIX", "TARGET", "FIXED_OUTPUT_DIR", "CELLS", "TUPLES",
              "SOURCES", "INPUTS", "COUNTS"):
    setattr(a, _name, globals()[_name])
h.CAMPAIGN, h.FIXED_OUTPUT_DIR = CAMPAIGN, FIXED_OUTPUT_DIR
h.WORKER_STARTED_SCHEMA = PREFIX + ".worker-started.v1"
h.R1_TRACKED_INPUTS = INPUTS
h.R0_ROOTS += (Path("/var/tmp/wh2-public-aa-calibration-r0"),)
a.POSTFLIGHT["copy_library"] = lambda p: h.artifact_receipt(Path(COPY_LIBRARY))
_base_validate_config = a.validate_config


def unit_key(replicate, unit):
    return digest("page-aa-r0:rep={}:unit={}".format(replicate, unit).encode("ascii"))


def unit_order(replicate):
    return sorted(range(9), key=lambda unit: unit_key(replicate, unit))


def roster():
    for rep in range(12):
        for unit in unit_order(rep):
            for index in range(4):
                scenario = (rep + unit + index) % 4
                for position, condition in enumerate(CONDITION_SEQUENCES[rep % 4]):
                    yield rep, unit, scenario, position, condition


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
        "roles": list(ROLES), "sibling_cpu": 56, "target_cpu": 120,
        "source_seed_base": h.SOURCE_SEED_BASE_TEXT, "timing_granularity": "whole-batch",
        "tuples": [{"K": k, "block_bytes": b, "scope": scope,
                    "scope_invocations_per_batch": 8192 // k} for k, b, scope in TUPLES],
        "unit_order": "sha256:page-aa-r0:rep=R:unit=U",
        "page_bytes": 4096, "source_phase": 2048, "counter_phase": 1024,
        "primary_roles": ["C", "L"], "primary_scenario": 0,
        "scenario_order": "(replicate+unit+index)%4",
        "scenarios": [{"scenario": n, "output_phases": list(phases), "shared_handle": n == 3}
                      for n, phases in enumerate(PHASES)],
    })
    return value


def description_hash():
    return digest(canonical(description()) + b"\n")


def address(value, where, allow_zero=False):
    if type(value) is not str or h.re.fullmatch(r"0x(?:0|[1-9a-f][0-9a-f]{0,15})", value) is None:
        fail(where + " malformed address")
    number = int(value, 16)
    if not allow_zero and not number:
        fail(where + " null address")
    return number


def validate_config(value, expected, described):
    obj(value, {"campaign", "cells", "compile", "copy_binding", "description_sha256",
                "runtime_library_maps_sha256", "runtime_library_path", "schema", "target_identity"},
        "config")
    binding = obj(value["copy_binding"], {"address", "elf_offset", "path"}, "copy binding")
    address(binding["address"], "copy binding")
    equal(binding["elf_offset"], "0x1a14c0", "resolved memcpy ELF offset")
    equal(binding["path"], COPY_LIBRARY, "resolved copy library")
    return _base_validate_config({key: item for key, item in value.items() if key != "copy_binding"},
                                 expected, described)


def validate_addresses(value, k, b, scenario, role):
    obj(value, {"arena", "source", "outputs", "counters", "handles", "span"}, "addresses")
    arena = address(value["arena"], "arena")
    equal(arena % 4096, 0, "arena alignment")
    span = ((k * b + 4096 + 4095) // 4096) * 4096
    equal(value["span"], span, "arena span")
    equal(address(value["source"], "source"), arena + 0x800, "source location")
    for field in ("outputs", "counters", "handles"):
        if type(value[field]) is not list or len(value[field]) != 2:
            fail("address pair differs: " + field)
    outputs = [address(item, "output") for item in value["outputs"]]
    counters = [address(item, "counter") for item in value["counters"]]
    handles = [address(item, "handle", allow_zero=role == "M") for item in value["handles"]]
    equal(outputs, [arena + span * (i + 1) + PHASES[scenario][i] for i in range(2)], "output locations")
    equal(counters, [arena + span * (i + 3) + 0x400 for i in range(2)], "counter locations")
    ranges = [(arena + 0x800, k * b)] + [(p, k * b) for p in outputs] + [(p, 512) for p in counters]
    if any(left + size > right for (left, size), (right, _) in zip(ranges, ranges[1:])):
        fail("arena regions overlap")
    if role == "M":
        equal(handles, [0, 0], "memcpy handles")
    elif (handles[0] == handles[1]) != (scenario == 3):
        fail("effective handle sharing differs")


PANEL_KEYS = (a.PANEL_KEYS - {"physical_scratch_addresses"}) | {
    "scenario", "addresses", "copy_binding_sha256"}


def validate_panel(value, key, sequence, config, oracles):
    obj(value, PANEL_KEYS, "panel")
    rep, unit, scenario, _, condition = key
    k, b = CELLS[unit // 3]
    role = ROLES[unit % 3]
    order = "ABBA" if condition % 2 == 0 else "BAAB"
    source_hash = h.SOURCE_SHA256[k, b, b]
    expected = {
        "K": k, "block_bytes": b, "condition": condition, "scenario": scenario,
        "description_sha256": description_hash(), "replicate": rep, "role": role,
        "mapping": "direct" if condition < 2 else "swapped", "order": order,
        "runtime_library_maps_sha256": config["runtime_library_maps_sha256"],
        "copy_binding_sha256": digest(canonical(config["copy_binding"])),
        "schema": PREFIX + ".panel.v1", "scope": "prebuilt-systematic",
        "sequence": sequence, "scope_invocations_per_batch": 8192 // k,
        "source_immutable_verified": True, "target_cpu": 120,
        "tuple_index": unit // 3, "unit_index": unit, "unit_key_sha256": unit_key(rep, unit),
        "left_output_sha256": source_hash, "right_output_sha256": source_hash,
    }
    for field, wanted in expected.items():
        equal(value[field], wanted, "panel." + field)
    validate_addresses(value["addresses"], k, b, scenario, role)
    slots = value["slots"]
    if type(slots) is not list or len(slots) != 8:
        fail("slot roster differs")
    elapsed = []
    for slot, side in zip(slots, h.sides_for(order)):
        obj(slot, {"elapsed_ns", "logical_lane", "physical_lane"}, "slot")
        equal(slot["logical_lane"], side, "logical lane")
        equal(slot["physical_lane"], (side == "right") ^ (condition // 2), "physical lane")
        elapsed.append(h.exact_int(slot["elapsed_ns"], 1, h.MAX_INT63, "elapsed"))
    return h.lane_contrast(elapsed, order)


def statistics(samples):
    result, failed = {"conditions": [], "matched_replicates": []}, []
    boundary = math.log1p(0.02)
    for scenario in range(4):
        for unit in range(9):
            primary = scenario == 0 and unit % 3 != 2
            for condition in range(4):
                summary = h.sample_summary([samples[rep, unit, scenario, condition] for rep in range(12)])
                result["conditions"].append({"unit": unit, "scenario": scenario,
                                             "condition": condition, "primary": primary, "summary": summary})
                if primary and not (-boundary < summary["lower95_log"] and summary["upper95_log"] < boundary):
                    failed.append("unit={}:condition={}".format(unit, condition))
            values = [math.fsum(samples[rep, unit, scenario, cond] for cond in range(4)) / 4
                      for rep in range(12)]
            summary = h.sample_summary(values)
            result["matched_replicates"].append({"unit": unit, "scenario": scenario,
                                                 "primary": primary, "summary": summary})
            if primary and not (-boundary < summary["lower95_log"] and summary["upper95_log"] < boundary):
                failed.append("unit={}:matched".format(unit))
    return result, failed


def terminal_value():
    return {"campaign": CAMPAIGN, "encode_call_count": 141557760,
            "measured_batch_count": 13824, "measured_invocation_count": 9732096,
            "panel_count": 1728, "record_count": 1730, "schema": PREFIX + ".terminal.v1",
            "status": "complete", "warmup_batch_count": 3456, "warmup_invocation_count": 2433024}


def parse_transcript(raw, expected, described):
    if len(raw) > h.MAX_RAW_BYTES:
        fail("raw exceeds byte limit")
    lines = raw.splitlines(keepends=True)
    equal(len(lines), 1730, "record count")
    records = [h.parse_canonical_line(line, "record") for line in lines]
    config = records[0]
    oracles = validate_config(config, expected, described)
    samples, retained = {}, {}
    for sequence, (panel, key) in enumerate(zip(records[1:-1], roster())):
        rep, unit, scenario, _, condition = key
        samples[rep, unit, scenario, condition] = validate_panel(panel, key, sequence, config, oracles)
        addresses = panel["addresses"]
        pair = rep, unit
        stable = {name: addresses[name] for name in ("arena", "source", "counters", "span")}
        if pair in retained:
            equal(stable, retained[pair], "retained arena")
        retained[pair] = stable
        effective = pair + (scenario,)
        if effective in retained:
            equal(addresses, retained[effective], "retained scenario")
        retained[effective] = addresses
    # All separate-handle scenarios retain the same two handles; shared uses 0.
    for rep in range(12):
        for unit in range(9):
            for scenario in (1, 2):
                equal(retained[rep, unit, scenario]["handles"], retained[rep, unit, 0]["handles"],
                      "retained handles across layouts")
            first = retained[rep, unit, 0]["handles"][0]
            equal(retained[rep, unit, 3]["handles"], [first, first], "shared handle")
    equal(records[-1], terminal_value(), "terminal")
    return statistics(samples)


def preflight(args):
    source = h.exact_absolute_directory(args.source_root, "source root")
    build = h.exact_absolute_directory(args.build_root, "build root")
    paths = {key + "_path": str(h.exact_absolute_file(getattr(args, key), key))
             for key in ("worker", "library", "compiler")}
    paths["controller_path"] = str(Path(__file__).resolve())
    equal(paths["controller_path"], str(source / "bench/Wh2PublicPageAaR0.py"), "controller entrypoint")
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
    p["copy_library_before"] = h.artifact_receipt(Path(COPY_LIBRARY))
    equal(artifacts(p), p["artifacts_before"], "preflight stable artifacts")
    return p


def validate_provenance(p, claim):
    for key in ("source_root", "build_root", "worker_path", "library_path", "compiler_path", "controller_path"):
        h.replay_path(p[key], "provenance." + key)
    equal(p["source_commit"], claim["source_commit"], "provenance commit")
    equal(p["worker_path"], str(Path(p["build_root"]) / TARGET), "worker build path")
    equal(p["library_path"], str(Path(p["build_root"]) / "libwirehair.so.2.0.0"), "library build path")
    equal(p["controller_path"], str(Path(p["source_root"]) / "bench/Wh2PublicPageAaR0.py"), "controller entrypoint")
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

    h.validate_stat_receipt(p["copy_library_before"], "copy library", COPY_LIBRARY, artifact=True)
    equal(p["copy_library_before"]["sha256"], COPY_LIBRARY_SHA256, "copy library hash")


for _function in (description, description_hash, roster, unit_order, validate_config,
                  parse_transcript, statistics, preflight, validate_provenance):
    setattr(a, _function.__name__, _function)


def main(argv=None):
    args = sys.argv[1:] if argv is None else argv
    if args == ["--selftest"]:
        equal(len(list(roster())), 1728, "roster")
        equal(sum(8 * 8192 for _ in roster()), COUNTS["expected_measured_encode_calls"], "calls")
        print("Page A/A controller selftest passed")
        return 0
    return a.main(args)


if __name__ == "__main__":
    sys.exit(main())
