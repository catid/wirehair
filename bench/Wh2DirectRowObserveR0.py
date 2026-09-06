#!/usr/bin/env python3
"""One-shot baseline-only timing observation; no performance acceptance claim.

Beads .66 freezes all samples and bounds. Build/selftest/receipt modes do not
execute the observation worker. Historical controllers and artifacts are read
only. The optional replay mode reduces recorded bytes, never remeasures them.
"""
import argparse
import importlib.util
import math
import os
from pathlib import Path
import statistics
import sys
import time


def sibling(name, filename):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(filename))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


M = sibling("_direct_observe_reviewed_tools", "Wh2ThueMorseMulRowsR0.py")
N, C, R = M.N, M.C, M.R
ROOT = M.ROOT
PROTOCOL = "wirehair.wh2.direct-row-observe-r0"
TARGET = "wh2_direct_row_observe_r0"
OUTPUT = Path("/var/tmp/wh2-direct-row-observe-r0")
HISTORICAL_RECEIPT = Path("/tmp/wh2-mulrows-final.e5iqJ8/receipt.json")
HISTORICAL_SHA = "487360da716332c7904b4f5eed15bdf9d5ef4065079c0aa1ad2e5197d9133d46"
SOURCES = tuple("bench/" + name for name in (
    "Wh2DirectRowObserveR0.cpp", "Wh2DirectRowObserveR0.py",
    "test_Wh2DirectRowObserveR0.py"))
METRICS = ("row_wall_ns", "thread_span_ns", "outer_wall_ns", "instrumentation_envelope_ns")
COUNTERS = ("minflt", "majflt", "nvcsw", "nivcsw")
require, exact, integer = M.require, M.exact, M.integer


def roster():
    return [(rep * 20 + step * 10 + position, rep, (rep + step) % 2, position)
            for rep in range(12) for step in range(2) for position in range(10)]


def expected_rows(lookup, deadline=None):
    C.time_left(deadline)
    require(type(lookup) is bytes and len(lookup) == 39936, "lookup shape")
    # Low IDs are the authenticated direct vectors, not another row campaign.
    return b"".join(lookup[((13 * j) & 1023) * 6:((13 * j) & 1023) * 6 + 6]
                    for j in range(64))


def output_checksum(rows, deadline=None):
    require(type(rows) is bytes and len(rows) == 384, "oracle row shape")
    value = 1469598103934665603
    for _ in range(240):
        C.time_left(deadline)
        for byte in rows:
            value = ((value * 1099511628211) ^ byte) & ((1 << 64) - 1)
    return value


def identity_and_runtime(raw):
    for suffix in ("before", "after"):
        try:
            R.H.validate_target_identity(raw["identity_" + suffix], suffix)
        except R.H.ValidationError as error:
            raise ValueError(str(error)) from error
    exact(raw["identity_before"], raw["identity_after"], "CPU identity drift")
    exact(raw["identity_before"]["derived"], {
        "ccd_id": 6, "complex_id": 6, "core_id": 50, "family": 26,
        "full_apic_id": 100, "initial_apic_id_8": 100,
        "logical_processors_per_package": 128, "model": 8, "package_id": 0,
        "stepping": 1, "thread_id": 0, "threads_per_core": 2}, "frozen target")
    runtime = raw["runtime_before"]
    R.keys(runtime, ("polynomial", "address", "ssse3", "avx2", "gfni", "avx512"), "GF runtime")
    integer(runtime["address"], 1, (1 << 64) - 1, "GF context address")
    exact({k: v for k, v in runtime.items() if k != "address"},
          {"polynomial": 333, "ssse3": 1, "avx2": 1, "gfni": 1, "avx512": 1}, "GF runtime")
    exact(runtime, raw["runtime_after"], "GF runtime drift")


def summarize(values):
    return {"min": min(values), "median": statistics.median(values), "max": max(values)}


def reduce_raw(raw, lookup=None, deadline=None):
    C.time_left(deadline)
    R.keys(raw, ("protocol", "schema", "outcome", "target_cpu", "lookup_sha256",
                 "expected_rows_hex", "expected_rows_sha256", "output_address", "output_stride",
                 "baseline_function_address", "identity_before", "identity_after", "runtime_before",
                 "runtime_after", "monotonic_resolution_ns", "thread_resolution_ns", "preflight_row_calls",
                 "callbacks", "row_calls", "checked_rows", "checksum", "sum_row_wall_ns",
                 "elapsed_ns", "records", "failure"), "raw")
    for key, value in (("protocol", PROTOCOL), ("schema", PROTOCOL + ".raw.v1"),
                       ("outcome", "DIAGNOSTIC_COMPLETE"), ("failure", None),
                       ("target_cpu", 50), ("lookup_sha256", M.LOOKUP_SHA),
                       ("output_stride", 38), ("preflight_row_calls", 64),
                       ("callbacks", 240), ("row_calls", 983040), ("checked_rows", 15360)):
        exact(raw[key], value, key)
    integer(raw["output_address"], 1, (1 << 64) - 64 * 38, "output address")
    integer(raw["baseline_function_address"], 1, (1 << 64) - 1, "baseline function address")
    elapsed = integer(raw["elapsed_ns"], 1, 10 * 10**9, "worker elapsed")
    for key in ("monotonic_resolution_ns", "thread_resolution_ns"):
        integer(raw[key], 1, 10**9, key)
    identity_and_runtime(raw)
    rows = expected_rows(M.authenticated_lookup(deadline) if lookup is None else lookup, deadline)
    exact(raw["expected_rows_hex"], rows.hex(), "independent direct-row oracle")
    exact(raw["expected_rows_sha256"], C.sha(rows), "row oracle hash")
    exact(raw["checksum"], output_checksum(rows, deadline), "consumed output checksum")
    require(type(raw["records"]) is list and len(raw["records"]) == 240, "complete record roster")
    derived = []
    previous_m3, previous_c1, previous_usage = 0, 0, [0] * 4
    total = 0
    first_m0, first_c0 = None, None
    for record, coordinate in zip(raw["records"], roster()):
        C.time_left(deadline)
        require(type(record) is list and len(record) == 12, "record shape")
        exact(record[:4], list(coordinate), "record roster/order")
        m0, c0, m1, m2, c1, m3 = [integer(value, 0, (1 << 63) - 1, "timestamp")
                                  for value in record[4:10]]
        require(m0 <= m1 < m2 <= m3 and c0 <= c1, "within-clock ordering")
        require(previous_m3 <= m0 and previous_c1 <= c0, "cross-record clock ordering")
        usage = []
        for snapshot in record[10:]:
            require(type(snapshot) is list and len(snapshot) == 4, "usage snapshot shape")
            usage.append([integer(value, 0, (1 << 63) - 1, "usage counter") for value in snapshot])
        require(all(previous_usage[i] <= usage[0][i] <= usage[1][i] for i in range(4)),
                "usage counter monotonicity")
        row_wall, thread_span, outer_wall = m2 - m1, c1 - c0, m3 - m0
        require(outer_wall <= elapsed and thread_span <= elapsed, "sample exceeds worker lifetime")
        if first_m0 is None:
            first_m0 = m0
            first_c0 = c0
        require(m3 - first_m0 <= elapsed, "record span exceeds worker lifetime")
        require(c1 - first_c0 <= elapsed, "thread record span exceeds worker lifetime")
        delta = {name: usage[1][i] - usage[0][i] for i, name in enumerate(COUNTERS)}
        derived.append({"index": coordinate[0], "rep": coordinate[1], "order": coordinate[2],
                        "position": coordinate[3], "warmup": coordinate[3] < 2,
                        "row_wall_ns": row_wall, "thread_span_ns": thread_span,
                        "outer_wall_ns": outer_wall, "instrumentation_envelope_ns": outer_wall - row_wall,
                        "counter_deltas": delta})
        total += row_wall
        previous_m3, previous_c1, previous_usage = m3, c1, usage[1]
    integer(raw["sum_row_wall_ns"], 1, 100000000, "cumulative inner wall bound")
    exact(raw["sum_row_wall_ns"], total, "cumulative inner wall sum")
    bins = []
    for start in range(0, 240, 10):
        sample = derived[start:start + 10]
        bins.append({"index": start // 10, "rep": sample[0]["rep"], "order": sample[0]["order"],
                     "first_record": start, "records": 10,
                     "metrics": {name: summarize([row[name] for row in sample]) for name in METRICS},
                     "counter_totals": {name: sum(row["counter_deltas"][name] for row in sample)
                                        for name in COUNTERS}})
    ratios = {}
    for name in METRICS:
        first, last = bins[0]["metrics"][name], bins[-1]["metrics"][name]
        ratios[name] = {key: first[key] / last[key] if last[key] else None
                        for key in ("min", "median", "max")}
    C.time_left(deadline)
    return {"protocol": PROTOCOL, "outcome": "DIAGNOSTIC_COMPLETE", "records": derived,
            "bins": bins, "first_over_last_bin_ratios": ratios,
            "sum_row_wall_ns": total, "all_samples_retained": True,
            "speed_claimed": False, "AA_equivalence_claimed": False,
            "production_promotion_claimed": False, "WH1_compared": False,
            "scope": "baseline_only_instrumented_observation",
            "interpretation": "Counters enclose clock calls. Thread span includes monotonic clock edges; "
                              "wall-minus-thread is not exact off-CPU time. Instrumentation and the shortened "
                              "roster change startup. No frequency/interrupt cause or historical rescue is established."}


def historical_receipt(deadline=None):
    data = C.read_regular(HISTORICAL_RECEIPT, 2 * 1024 * 1024, deadline=deadline)
    require(C.sha(data) == HISTORICAL_SHA, "historical receipt identity")
    old = C.strict_json(data)
    current = M.current_receipt(Path(old["build_directory"]), deadline)
    # New additive commits and an independently pinned supported interpreter
    # may differ from the historical launch provenance. Authenticate the old
    # interpreter explicitly before normalizing those provenance-only fields.
    # Every historical source/build/runtime/dependency byte remains exact.
    exact(N.pin(old["interpreter"]["path"], deadline, installed=True),
          old["interpreter"], "historical interpreter pin")
    current["source_commit"] = old["source_commit"]
    current["interpreter"] = old["interpreter"]
    current["python_version"] = old["python_version"]
    exact(current, old, "historical build/source pins changed")
    return old


def build_commands(build):
    flags = [N.COMPILER, "-DWIREHAIR_STATIC=1", "-I" + str(ROOT), "-I" + str(ROOT / "include"),
             "-I" + str(M.BASE), "-O3", "-DNDEBUG", "-std=gnu++11", "-Wall", "-Wextra",
             "-Wpedantic", "-Werror", "-march=native"]
    names = ("Wh2ThueMorseNativeCodec", "Wh2DirectRowObserveR0")
    commands = [flags + ["-MD", "-MF", str(build / (name + ".d")), "-c",
                         str(ROOT / "bench" / (name + ".cpp")), "-o", str(build / (name + ".o"))]
                for name in names]
    commands.append([N.COMPILER, "-O3", "-DNDEBUG"] + [str(build / (name + ".o")) for name in names] +
                    [str(path) for path in M.HELPERS] + [str(M.BASE / "libwirehair.a"), "-lm",
                                                       "-o", str(build / TARGET)])
    return commands


def build_inputs(dependencies):
    result = set(dependencies) | {ROOT / name for name in SOURCES} | set(M.HELPERS)
    result.update((M.BASE / "libwirehair.a", Path(N.COMPILER).resolve(strict=True)))
    return result


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
    outputs = [Path(command[-1]) for command in commands]
    outputs += [Path(command[command.index("-MF") + 1]) for command in commands[:2]]
    C.write_new(build / "commands.json", C.canonical(commands))
    C.write_new(build / "build_inputs.json", C.canonical({"inputs": before,
                "outputs": [N.pin(path) for path in sorted(outputs)]}))
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
        data = C.read_regular(ROOT / name, 2 * 1024 * 1024, deadline=deadline)
        require(data == C.git_bytes("show", commit + ":" + name, deadline=deadline), "source blob differs")
        sources[name] = C.sha(data)
    dependencies, outputs = set(), {build / TARGET}
    for command in build_commands(build)[:2]:
        dep = Path(command[command.index("-MF") + 1])
        outputs.update((dep, Path(command[-1])))
        dependencies.update(M.dependency_names(C.read_regular(dep, 256 * 1024, deadline=deadline), command[-1]))
    require(ROOT / "gf256.h" in dependencies and M.BASE / "Wh2ThueMorseNativeDataR0.h" in dependencies,
            "incomplete compiler dependencies")
    require(ROOT / "bench/Wh2ThueMorseMulRowsR0.cpp" not in dependencies, "real candidate kernel dependency")
    M.validate_build_manifest(C.strict_json(C.read_regular(build / "build_inputs.json", 2 * 1024 * 1024, deadline=deadline)),
                              build_inputs(dependencies), outputs, deadline)
    symbols = N.checked(["/usr/bin/nm", "--defined-only", str(build / TARGET)], deadline)
    names = [line.split()[-1] for line in symbols.decode("ascii").splitlines() if line.split()]
    for name in ("GF256Ctx", "gf256_init_", "gf256_mulset_multi_mem", "gf256_muladd_mem"):
        require(names.count(name) == 1, "shared GF symbol count")
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
            "max_stdout_bytes": 1048576, "max_stderr_bytes": 65536, "max_total_bytes": 4194304}


def run(receipt_path):
    started = time.monotonic()
    receipt_bytes = C.read_regular(receipt_path, 2 * 1024 * 1024, deadline=started + 4)
    receipt = C.strict_json(receipt_bytes)
    require(C.canonical(current_receipt(Path(receipt["build_directory"]), started + 4)) == receipt_bytes, "launch receipt differs")
    OUTPUT.mkdir(mode=0o700, exist_ok=False)
    members, total = {}, 0

    def publish(name, data):
        nonlocal total
        total += len(data)
        require(total <= 4194304, "aggregate output cap")
        C.write_new(OUTPUT / name, data)
        members[name] = {"bytes": len(data), "sha256": C.sha(data)}

    publish("CLAIM.json", receipt_bytes)
    outcome, failure, analysis = "INVALID", None, None
    stdout, stderr, observation = b"", b"", {}
    try:
        observation = C.capture_worker(receipt["worker_argv"], min(time.monotonic() + 10, started + 14),
                                       stdout_cap=1048576, stderr_cap=65536)
        stdout, stderr = observation.pop("stdout"), observation.pop("stderr")
        require(not observation["failure"] and observation["returncode"] == 0 and not stderr, "worker failed")
        duration = observation["elapsed_seconds"]
        require(type(duration) in (int, float) and math.isfinite(duration) and 0 < duration <= 10, "observed worker bound")
        raw = C.strict_json(stdout)
        require(type(raw) is dict and type(raw.get("elapsed_ns")) is int and
                0 < raw["elapsed_ns"] <= duration * 10**9, "internal interval exceeds worker lifetime")
        require(stdout.endswith(b"\n"), "raw framing")
        analysis = reduce_raw(raw, deadline=started + 17)
        outcome = analysis["outcome"]
    except Exception as error:
        failure = type(error).__name__ + ": " + str(error)[:1000]
    try:
        require(C.canonical(current_receipt(Path(receipt["build_directory"]), started + 19)) == receipt_bytes, "post-run receipt differs")
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
        raw = C.strict_json(C.read_regular(args.replay, 1048576))
        print(C.canonical(reduce_raw(raw)).decode("ascii"))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print(type(error).__name__ + ": " + str(error)[:1000], file=sys.stderr)
        sys.exit(1)
