#!/usr/bin/env python3
"""One-shot, equation-neutral row-mapping cost screen (Beads .65).

Build/selftest and receipt modes never run the timed worker. The old native R0
sources, build and evidence stay immutable; its timing outcome is not reused.
"""
import argparse
import importlib.util
import math
import os
from pathlib import Path
import shlex
import sys
import time


def sibling(name, filename):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(filename))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


N = sibling("_mulrows_native_tools", "Wh2ThueMorseNativeRunR0.py")
R = sibling("_mulrows_native_math", "Wh2ThueMorseNativeReduceR0.py")

C = N.C
ROOT = N.ROOT
PROTOCOL = "wirehair.wh2.thue-morse-mulrows-r0"
OUTPUT = Path("/var/tmp/wh2-thue-morse-mulrows-r0")
BASE = Path("/tmp/wh2-thue-native-compile.6CJYNz")
BASE_RECEIPT = Path("/tmp/wh2-thue-native-final.Y9drbS/build.json")
BASE_RECEIPT_SHA = "53795b83ac19026f9f19ed64657245ac9ef2348784ab1fdbd731cca8d49bce04"
LOOKUP_SHA = "27b105e1449bec190bd3c83f07feefa639cd32bc356baebfb03828ea7cbccb6d"
TARGET = "wh2_thue_morse_mulrows_r0"
SOURCES = tuple("bench/" + name for name in (
    "Wh2ThueMorseMulRowsR0.h", "Wh2ThueMorseMulRowsR0.cpp",
    "Wh2ThueMorseMulRowsR0Test.cpp", "Wh2ThueMorseMulRowsR0Bench.cpp",
    "Wh2ThueMorseMulRowsR0.py", "test_Wh2ThueMorseMulRowsR0.py"))
HELPERS = tuple(BASE / "CMakeFiles/wh2_thue_morse_native_r0.dir/bench" /
                (name + ".cpp.o") for name in (
                    "Wh2RdpruTargetIdentityV2", "Wh2PublicBorrowedTargetIdentity",
                    "Wh2FrozenTrace"))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def exact(actual, expected, label):
    require(C.canonical(actual) == C.canonical(expected), label)


def integer(value, low, high, label):
    require(type(value) is int and low <= value <= high, label)
    return value


def ids(family):
    result = []
    for j in range(64):
        value = (13 * j) & 1023
        if family >= 1:
            value |= (1 + (17 * j) % 127) << 10
        if family >= 2:
            value |= (1 + (23 * j) % 127) << 17
        if family >= 3:
            value |= (1 + (37 * j) % 255) << 24
        result.append(value)
    return result


def roster():
    return [(rep, family, comparison, (rep + family + comparison + step) % 2)
            for rep in range(12) for family in range(4)
            for comparison in range(3) for step in range(2)]


def multiply(x, y):
    result = 0
    for _ in range(8):
        if y & 1:
            result ^= x
        x <<= 1
        if x & 256:
            x ^= 333
        y >>= 1
    return result


def expected_rows(lookup, deadline=None):
    require(type(lookup) is bytes and len(lookup) == 39936, "oracle lookup shape")
    rows = bytearray()
    for family in range(4):
        C.time_left(deadline)
        for block in ids(family):
            high, m17, m10, low = block >> 24, (block >> 17) & 127, (block >> 10) & 127, block & 1023
            p17 = bin(high).count("1") & 1
            p10 = p17 ^ (bin(m17).count("1") & 1)
            p0 = p10 ^ (bin(m10).count("1") & 1)
            row = lookup[p0 * 6144 + low * 6:p0 * 6144 + low * 6 + 6]
            for index, offset in ((m10, 12288 + p10 * 4608),
                                  (m17, 21504 + p17 * 4608), (high, 30720)):
                if not index:
                    continue
                matrix = lookup[offset + index * 36:offset + index * 36 + 36]
                result = []
                for r in range(6):
                    value = 0
                    for c in range(6):
                        value ^= multiply(matrix[r * 6 + c], row[c])
                    result.append(value)
                row = bytes(result)
            rows.extend(row)
    return bytes(rows)


def output_checksum(rows, deadline=None):
    require(len(rows) == 1536, "expected row bytes")
    value = 1469598103934665603
    for _, family, _, _ in roster():
        C.time_left(deadline)
        for _ in range(10):
            for byte in rows[family * 384:(family + 1) * 384]:
                value = ((value * 1099511628211) ^ byte) & ((1 << 64) - 1)
    return value


def authenticated_lookup(deadline=None):
    D = sibling("_mulrows_fixed_data", "Wh2ThueMorseNativeDataR0.py")
    data = D.extract_data(D.authenticated_inputs(deadline=deadline))
    C.time_left(deadline)
    lookup, _ = D.build_lookup(data["pair"])
    C.time_left(deadline)
    require(C.sha(lookup) == LOOKUP_SHA, "authenticated lookup identity")
    return lookup


def reduce_raw(raw, lookup=None, deadline=None):
    C.time_left(deadline)
    R.keys(raw, ("schema", "protocol", "outcome", "lookup_sha256", "panels", "warm_callbacks",
                 "measured_callbacks", "row_calls", "checksum", "elapsed_ns",
                 "identity_before", "identity_after", "runtime_before", "runtime_after",
                 "target_cpu", "calls_per_callback", "preflight_row_calls", "checked_rows",
                 "expected_rows_hex", "expected_rows_sha256", "output_address", "output_stride",
                 "baseline_function_address", "candidate_function_address"), "raw")
    exact(raw["schema"], PROTOCOL + ".raw.v1", "schema")
    exact(raw["protocol"], PROTOCOL, "protocol")
    exact(raw["outcome"], "COMPLETE", "incomplete worker")
    exact(raw["lookup_sha256"], LOOKUP_SHA, "lookup")
    exact(raw["warm_callbacks"], 576, "warm callbacks")
    exact(raw["measured_callbacks"], 2304, "measured callbacks")
    exact(raw["row_calls"], 2880 * 4096, "row calls")
    exact(raw["target_cpu"], 50, "target CPU")
    exact(raw["calls_per_callback"], 4096, "calls per callback")
    exact(raw["preflight_row_calls"], 256, "preflight calls")
    exact(raw["checked_rows"], 184320, "checked rows")
    exact(raw["output_stride"], 38, "output stride")
    integer(raw["output_address"], 1, (1 << 64) - 64 * 38, "output address")
    for key in ("baseline_function_address", "candidate_function_address"):
        integer(raw[key], 1, (1 << 64) - 1, key)
    require(raw["baseline_function_address"] != raw["candidate_function_address"], "aliased implementations")
    rows = expected_rows(authenticated_lookup(deadline) if lookup is None else lookup, deadline)
    exact(raw["expected_rows_hex"], rows.hex(), "independent row oracle")
    exact(raw["expected_rows_sha256"], C.sha(rows), "row oracle hash")
    exact(raw["checksum"], output_checksum(rows, deadline), "consumed output checksum")
    elapsed = integer(raw["elapsed_ns"], 1, 30 * 10**9, "worker elapsed")
    try:
        R.H.validate_target_identity(raw["identity_before"], "before")
        R.H.validate_target_identity(raw["identity_after"], "after")
    except R.H.ValidationError as error:
        raise ValueError(str(error)) from error
    exact(raw["identity_before"], raw["identity_after"], "CPU identity drift")
    # This is the same native host/runtime as the frozen baseline, not an ISA
    # ablation. Target identity validates all advertised APIC sources.
    exact(raw["identity_before"]["derived"], {
        "ccd_id": 6, "complex_id": 6, "core_id": 50, "family": 26,
        "full_apic_id": 100, "initial_apic_id_8": 100,
        "logical_processors_per_package": 128, "model": 8, "package_id": 0,
        "stepping": 1, "thread_id": 0, "threads_per_core": 2}, "frozen target")
    R.keys(raw["runtime_before"], ("polynomial", "address", "ssse3", "avx2", "gfni", "avx512"), "GF runtime")
    integer(raw["runtime_before"]["address"], 1, (1 << 64) - 1, "GF context address")
    exact({k: v for k, v in raw["runtime_before"].items() if k != "address"},
          {"polynomial": 333, "ssse3": 1, "avx2": 1, "gfni": 1, "avx512": 1}, "GF runtime")
    exact(raw["runtime_before"], raw["runtime_after"], "GF runtime drift")
    require(type(raw["panels"]) is list and len(raw["panels"]) == 288, "panel count")
    groups = {}
    total_ns = 0
    for panel, expected in zip(raw["panels"], roster()):
        C.time_left(deadline)
        require(type(panel) is list and len(panel) == 5, "panel shape")
        exact(panel[:4], list(expected), "panel roster/order")
        times = panel[4]
        require(type(times) is list and len(times) == 8, "slot count")
        for value in times:
            total_ns += integer(value, 1, elapsed, "slot time")
        _, family, comparison, order = expected
        contrast = R.H.lane_contrast(times, "ABBA" if order == 0 else "BAAB")
        groups.setdefault((family, comparison, order), []).append(contrast)
    require(total_ns <= elapsed, "overlapping/impossible timing intervals")
    statistics, failed_controls, failed_candidates = [], [], []
    limit = math.log1p(.02)
    for key in sorted(groups):
        family, comparison, order = key
        values = groups[key]
        # Raw contrasts are left/right. B/C becomes C/B for readable costs.
        estimate = R.H.sample_summary([-v for v in values] if comparison == 2 else values)
        item = {"family": family, "comparison": comparison, "order": order,
                "estimate": estimate}
        if comparison != 2:
            passed = estimate["lower95_log"] > -limit and estimate["upper95_log"] < limit
            if not passed:
                failed_controls.append(list(key))
        else:
            passed = estimate["upper95"] < 1 if family == 3 else estimate["upper95"] <= 1.02
            if not passed:
                failed_candidates.append(list(key))
        item["passed"] = passed
        statistics.append(item)
    return {"protocol": PROTOCOL,
            "outcome": "CONTROL_FAIL" if failed_controls else "SHORT_FAIL" if failed_candidates else "PASS",
            "statistics": statistics, "failed_controls": failed_controls,
            "failed_candidates": failed_candidates, "WH1_compared": False,
            "production_promotion_claimed": False, "scope": "row_mapping_only"}


def build_commands(build):
    flags = [N.COMPILER, "-DWIREHAIR_STATIC=1", "-I" + str(ROOT),
             "-I" + str(ROOT / "include"), "-I" + str(BASE), "-O3", "-DNDEBUG",
             "-std=gnu++11", "-Wall", "-Wextra", "-Wpedantic", "-Werror", "-march=native"]
    names = ("Wh2ThueMorseNativeCodec", "Wh2ThueMorseMulRowsR0",
             "Wh2ThueMorseMulRowsR0Bench", "Wh2ThueMorseMulRowsR0Test")
    commands = [flags + ["-MD", "-MF", str(build / (name + ".d")), "-c",
                         str(ROOT / "bench" / (name + ".cpp")), "-o", str(build / (name + ".o"))]
                for name in names]
    common = [str(build / (name + ".o")) for name in names[:2]]
    commands.append([N.COMPILER, "-O3", "-DNDEBUG"] + common +
                    [str(build / (names[2] + ".o"))] + [str(p) for p in HELPERS] +
                    [str(BASE / "libwirehair.a"), "-lm", "-o", str(build / TARGET)])
    commands.append([N.COMPILER, "-O3", "-DNDEBUG"] + common +
                    [str(build / (names[3] + ".o")), str(BASE / "libwirehair.a"),
                     "-lm", "-o", str(build / (TARGET + "_test"))])
    return commands


def build_only(build):
    require(build.is_absolute() and build.parent.is_dir() and not build.exists(), "fresh explicit build path")
    build.mkdir(mode=0o700)
    commands = build_commands(build)
    old_receipt()
    inputs = {ROOT / name for name in SOURCES}
    inputs.update(HELPERS)
    inputs.update((BASE / "libwirehair.a", Path(N.COMPILER).resolve(strict=True)))
    for command in commands[:4]:
        preprocessor = command[:command.index("-MD")] + ["-M", "-MT", command[-1],
                         command[command.index("-c") + 1]]
        inputs.update(dependency_names(N.checked(preprocessor, time.monotonic() + 30), command[-1]))
    before = [N.pin(path, installed=True) for path in sorted(inputs)]
    for command in commands:
        N.checked(command, time.monotonic() + 60)
    exact([N.pin(path, installed=True) for path in sorted(inputs)], before, "build input changed during compilation")
    old_receipt()
    C.write_new(build / "commands.json", C.canonical(commands))
    outputs = [Path(command[-1]) for command in commands]
    outputs += [Path(command[command.index("-MF") + 1]) for command in commands[:4]]
    C.write_new(build / "build_inputs.json", C.canonical({"inputs": before,
                "outputs": [N.pin(path) for path in sorted(outputs)]}))
    for command in ([str(build / (TARGET + "_test"))], [str(build / TARGET), "--selftest"]):
        print(N.checked(command, time.monotonic() + 30).decode("ascii"), end="")


def dependency_names(content, expected_target):
    target, separator, names = content.decode("utf-8").replace("\\\n", " ").partition(":")
    require(separator and target.strip() == expected_target, "dependency target")
    return {Path(name).resolve(strict=True) for name in shlex.split(names)}


def validate_build_manifest(manifest, inputs, outputs, deadline=None):
    R.keys(manifest, ("inputs", "outputs"), "build-time manifest")
    exact(manifest["inputs"], [N.pin(path, deadline, installed=True) for path in sorted(inputs)],
          "stale build input")
    exact(manifest["outputs"], [N.pin(path, deadline) for path in sorted(outputs)], "stale build output")


def old_receipt(deadline=None):
    data = C.read_regular(BASE_RECEIPT, 2 * 1024 * 1024, deadline=deadline)
    require(C.sha(data) == BASE_RECEIPT_SHA, "baseline receipt changed")
    old = C.strict_json(data)
    for name, digest in old["sources"].items():
        require(C.sha(C.read_regular(ROOT / name, 4 * 1024 * 1024, deadline=deadline)) == digest,
                "frozen baseline source changed: " + name)
    pins = old["system_headers"] + old["build_files"] + old["runtime_libraries"] + old["tools"]
    pins += [old["compiler"], old["interpreter"], old["generated_header"]]
    for item in pins:
        exact(N.pin(item["path"], deadline, installed=True), item, "baseline dependency pin")
    return old


def current_receipt(build, deadline=None):
    require(sys.version_info[:2] in ((3, 8), (3, 12)), "unvalidated interpreter")
    build = build.resolve(strict=True)
    old = old_receipt(deadline)
    commands = C.read_regular(build / "commands.json", 65536, deadline=deadline)
    exact(C.strict_json(commands), build_commands(build), "build commands")
    commit = C.git("rev-parse", "HEAD", deadline=deadline)
    require(commit == C.git("rev-parse", "@{upstream}", deadline=deadline), "sources not pushed")
    require(not C.git("status", "--porcelain", "--", *SOURCES, deadline=deadline), "dirty new sources")
    exact(sorted(C.git("ls-files", "--", *SOURCES, deadline=deadline).splitlines()),
          sorted(SOURCES), "untracked sources")
    sources = {}
    for name in SOURCES:
        data = C.read_regular(ROOT / name, 2 * 1024 * 1024, deadline=deadline)
        require(data == C.git_bytes("show", commit + ":" + name, deadline=deadline), "source blob differs")
        sources[name] = C.sha(data)
    dependencies = set()
    build_files = {build / TARGET, build / (TARGET + "_test"), build / "commands.json", build / "build_inputs.json"}
    for command in build_commands(build)[:4]:
        dep = Path(command[command.index("-MF") + 1])
        build_files.update((dep, Path(command[-1])))
        content = C.read_regular(dep, 256 * 1024, deadline=deadline)
        dependencies.update(dependency_names(content, command[-1]))
    require(ROOT / "gf256.h" in dependencies and BASE / "Wh2ThueMorseNativeDataR0.h" in dependencies,
            "incomplete compiler dependencies")
    build_inputs = dependencies | {ROOT / name for name in SOURCES} | set(HELPERS)
    build_inputs.update((BASE / "libwirehair.a", Path(N.COMPILER).resolve(strict=True)))
    validate_build_manifest(C.strict_json(C.read_regular(build / "build_inputs.json", 2 * 1024 * 1024, deadline=deadline)),
                            build_inputs, build_files - {build / "commands.json", build / "build_inputs.json"}, deadline)
    symbols = N.checked(["/usr/bin/nm", "--defined-only", str(build / TARGET)], deadline)
    symbol_names = [line.split()[-1] for line in symbols.decode("ascii").splitlines() if line.split()]
    for symbol in ("GF256Ctx", "gf256_init_", "gf256_mulset_multi_mem", "gf256_muladd_mem"):
        require(symbol_names.count(symbol) == 1, "shared GF symbol count")
    runtimes = N.runtime_paths(N.checked(["/usr/bin/ldd", str(build / TARGET)], deadline))
    exact(runtimes, sorted(item["path"] for item in old["runtime_libraries"]), "runtime closure")
    return {"protocol": PROTOCOL, "source_commit": commit, "sources": sources,
            "baseline_receipt_sha256": BASE_RECEIPT_SHA, "build_directory": str(build),
            "build_files": [N.pin(path, deadline) for path in sorted(build_files)],
            "dependencies": [N.pin(path, deadline, installed=True) for path in sorted(dependencies)],
            "helpers": [N.pin(path, deadline) for path in HELPERS],
            "symbols_sha256": C.sha(symbols),
            "interpreter": N.pin(Path(sys.executable).resolve(strict=True), deadline, installed=True),
            "python_version": list(sys.version_info[:3]),
            "worker_argv": [str(build / TARGET), "--worker"], "environment": C.ENV}


def run(receipt_path):
    started = time.monotonic()
    receipt_bytes = C.read_regular(receipt_path, 2 * 1024 * 1024, deadline=started + 5)
    receipt = C.strict_json(receipt_bytes)
    require(C.canonical(current_receipt(Path(receipt["build_directory"]), started + 5)) == receipt_bytes,
            "launch receipt differs")
    OUTPUT.mkdir(mode=0o700, exist_ok=False)
    members, total = {}, 0

    def publish(name, data):
        nonlocal total
        total += len(data)
        require(total <= 4 * 1024 * 1024, "aggregate output cap")
        C.write_new(OUTPUT / name, data)
        members[name] = {"bytes": len(data), "sha256": C.sha(data)}

    publish("CLAIM.json", receipt_bytes)
    outcome, failure, analysis = "INVALID", None, None
    stdout, stderr, observation = b"", b"", {}
    try:
        observation = C.capture_worker(receipt["worker_argv"], min(time.monotonic() + 30, started + 35),
                                       stdout_cap=2 * 1024 * 1024, stderr_cap=65536)
        stdout, stderr = observation.pop("stdout"), observation.pop("stderr")
        require(not observation["failure"] and observation["returncode"] == 0 and not stderr, "worker failed")
        require(type(observation["elapsed_seconds"]) in (int, float) and
                math.isfinite(observation["elapsed_seconds"]) and 0 < observation["elapsed_seconds"] <= 30,
                "observed worker wall-time bound")
        raw = C.strict_json(stdout)
        require(type(raw) is dict and type(raw.get("elapsed_ns")) is int and
                0 < raw["elapsed_ns"] <= observation["elapsed_seconds"] * 10**9,
                "internal interval exceeds observed worker lifetime")
        analysis = reduce_raw(raw, deadline=started + 38)
        require(analysis and stdout.endswith(b"\n"), "raw framing")
        outcome = analysis["outcome"]
    except Exception as error:
        failure = type(error).__name__ + ": " + str(error)[:1000]
    try:
        require(C.canonical(current_receipt(Path(receipt["build_directory"]), started + 39)) == receipt_bytes,
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
    if time.monotonic() - started >= 39.5:
        outcome, failure = "INVALID", "controller deadline"
    summary = {"protocol": PROTOCOL, "outcome": outcome, "failure": failure,
               "worker": observation, "elapsed_seconds": time.monotonic() - started,
               "WH1_compared": False, "production_promotion_claimed": False}
    publish("summary.json", C.canonical(summary))
    publish("COMPLETE.json", C.canonical({"protocol": PROTOCOL, "outcome": outcome, "files": members}))
    fd = os.open(str(OUTPUT), os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)
    require(time.monotonic() - started < 40, "whole CLI deadline")
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
        print(C.canonical(reduce_raw(C.strict_json(C.read_regular(args.replay, 2 * 1024 * 1024)))).decode("ascii"))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print(type(error).__name__ + ": " + str(error)[:1000], file=sys.stderr)
        sys.exit(1)
