#!/usr/bin/env python3
"""Create-only native K6 screen controller; scientific freeze is Beads .63.

Receipt construction only inspects an already-built executable. It never runs
the candidate. The sole --run invocation is authorized by an exact, pushed
source/build/input receipt and permanently consumes a fresh output namespace.
"""
import argparse
import importlib.util
import os
from pathlib import Path
import re
import stat
import subprocess
import sys
import time


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def sibling(name, filename):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


C = sibling("_tm_native_safe_io", "Wh2NoncommutingRadixRunR0.py")
PROTOCOL = "wirehair.wh2.thue-morse-native-r0"
OUTPUT = Path("/var/tmp/wh2-thue-morse-native-r0")
TARGET = "wh2_thue_morse_native_r0"
COMPILER = "/usr/bin/x86_64-linux-gnu-g++-13"
COMPILER_SHA = "1353e9bdd29a7295c7226bf6c63abccce056d8cac31f112e5cdbecc3f28c2769"
MAX_STDOUT = 8 * 1024 * 1024
MAX_STDERR = 1024 * 1024
MAX_TOTAL = 16 * 1024 * 1024
MAX_RECEIPT = 2 * 1024 * 1024
SOURCES = (
    "CMakeLists.txt", "AGENTS.md",
    "bench/Wh2ThueMorseNativeCodec.h", "bench/Wh2ThueMorseNativeCodec.cpp",
    "bench/Wh2ThueMorseNativeCodecTest.cpp", "bench/Wh2ThueMorseNativeR0.cpp",
    "bench/Wh2ThueMorseNativeDataR0.py", "bench/test_Wh2ThueMorseNativeDataR0.py",
    "bench/Wh2ThueMorseNativeRunR0.py", "bench/test_Wh2ThueMorseNativeRunR0.py",
    "bench/Wh2ThueMorseNativeReduceR0.py", "bench/test_Wh2ThueMorseNativeReduceR0.py",
    "bench/Wh2NoncommutingRadixRunR0.py", "bench/Wh2NoncommutingRadixR0.py",
    "bench/Wh2ThueMorseR0.py", "bench/Wh2ThueMorseRecoveryR0.py",
    "bench/Wh2ThueMorseRecoveryHistoryR0.py", "bench/Wh2ThueMorseRecoveryRunR0.py",
    "bench/test_Wh2ThueMorseR0.py", "bench/test_Wh2NoncommutingRadixR0.py",
    "bench/test_Wh2ThueMorseRecoveryR0.py", "bench/Wh2FrozenTraceTest.cpp",
    "bench/wh2_benchmark_contract_v4.json", "bench/wh2_mix2_seed_repair_contract.py",
    "bench/wh2_precodefail_work_screen.py",
    "bench/Wh2PublicBorrowedCurrentVsWh1R1.py",
)
CACHE_EXPECTED = {
    "CMAKE_BUILD_TYPE": "Release", "BUILD_SHARED_LIBS": "OFF",
    "WIREHAIR_BUILD_BOTH": "OFF", "BUILD_CODEC_V2": "OFF",
    "MARCH_NATIVE": "ON", "WH_LTO": "OFF", "WH_PGO_MODE": "OFF",
    "WIREHAIR_STRICT_WARNINGS": "ON", "WIREHAIR_ENABLE_THUE_MORSE_NATIVE_R0": "ON",
    "WIREHAIR_STATIC_PIC": "ON", "CMAKE_CXX_FLAGS": "",
    "CMAKE_CXX_FLAGS_RELEASE": "-O3 -DNDEBUG", "CMAKE_GENERATOR": "Ninja",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def checked(argv, deadline=None, cap=4 * 1024 * 1024):
    result = subprocess.run(argv, cwd=str(ROOT), env=C.ENV,
                            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, timeout=C.time_left(deadline),
                            check=True)
    require(not result.stderr and len(result.stdout) <= cap, "tool output/diagnostic bound")
    return result.stdout


def pin(path, deadline=None, installed=False, cap=64 * 1024 * 1024):
    path = Path(path)
    data = C.read_regular(path, cap, single_link=not installed, deadline=deadline)
    return {"path": str(path), "bytes": len(data), "sha256": C.sha(data)}


def parse_cache(data):
    result = {}
    for line in data.decode("utf-8").splitlines():
        if not line or line.startswith(("#", "//")):
            continue
        key_type, sep, value = line.partition("=")
        key, colon, _ = key_type.partition(":")
        require(bool(sep and colon and key) and key not in result, "CMake cache syntax")
        result[key] = value
    return result


def dependency_paths(data):
    """Ninja's compiler-emitted dependencies include system headers (-MD)."""
    paths = set()
    for line in data.decode("utf-8").splitlines():
        if line.startswith("    "):
            name = line[4:]
            require(name and "\x00" not in name, "empty/invalid dependency")
            paths.add(name)
        elif line:
            require(re.fullmatch(r".+: #deps [0-9]+, deps mtime [0-9]+ \(VALID\)", line)
                    is not None, "stale or unrecognized compiler dependency entry")
    require(bool(paths), "missing compiler dependency closure")
    return sorted(paths)


def runtime_paths(data):
    paths = set()
    for line in data.decode("ascii").splitlines():
        value = line.strip()
        if value.startswith("linux-vdso.so.1 "):
            continue
        match = re.fullmatch(r"(?:[^ ]+ => )?(/[^ ]+) \(0x[0-9a-f]+\)", value)
        require(match is not None, "unresolved or unexpected runtime dependency")
        paths.add(str(Path(match.group(1)).resolve(strict=True)))
    require(len(paths) >= 3, "incomplete dynamic runtime closure")
    return sorted(paths)


def verify_symbols(data):
    definitions = []
    for line in data.decode("utf-8").splitlines():
        parts = line.split(None, 2)
        if len(parts) == 3:
            definitions.append(parts[2])
    for symbol in ("GF256Ctx", "gf256_init_", "gf256_mulset_multi_mem",
                   "gf256_muladd_mem", "wirehair_encode", "wirehair_decode",
                   "wirehair_v2_encode", "wirehair_v2_decode"):
        require(definitions.count(symbol) == 1, "missing/duplicate shared symbol: " + symbol)


def verify_compile_commands(command_text, commit):
    macros = re.findall(r'WIREHAIR_WH2_SOURCE_GIT_COMMIT=\\?"([0-9a-f]{40})\\?"', command_text)
    require(macros and set(macros) == {commit}, "compiled source commit requires fresh configure")
    compiled = [line for line in command_text.splitlines() if re.search(r"\s-c\s", line)]
    require(bool(compiled), "missing native compile commands")
    for line in compiled:
        standards = re.findall(r"(?:^|\s)-std=([^\s]+)", line)
        require(standards in (["gnu++11"], ["c++11"]), "native compile is not frozen C++11")


def current_receipt(build, deadline=None):
    build = Path(build).resolve(strict=True)
    require(build != ROOT and build.is_dir(), "explicit non-root build directory required")
    compiler = pin(Path(COMPILER).resolve(strict=True), deadline, installed=True)
    require(compiler["sha256"] == COMPILER_SHA, "frozen compiler changed")
    cache_bytes = C.read_regular(build / "CMakeCache.txt", 2 * 1024 * 1024, deadline=deadline)
    cache = parse_cache(cache_bytes)
    for key, expected in CACHE_EXPECTED.items():
        require(cache.get(key) == expected, "frozen cache mismatch: " + key)
    require(Path(cache["CMAKE_HOME_DIRECTORY"]).resolve(strict=True) == ROOT,
            "foreign source tree")
    require(Path(cache["CMAKE_CXX_COMPILER"]).resolve(strict=True) == Path(compiler["path"]),
            "wrong build compiler")
    dry = checked(["/usr/bin/ninja", "-C", str(build), "-n", TARGET], deadline)
    require(dry.decode("utf-8").splitlines()[-1:] == ["ninja: no work to do."],
            "native target needs rebuilding")
    commands = checked(["/usr/bin/ninja", "-C", str(build), "-t", "commands", TARGET], deadline)
    command_text = commands.decode("utf-8")
    require("libwirehair.a" in command_text and
            not any(name in command_text for name in
                    ("libwirehair_v2_", "WIREHAIR_V2_TEST", "WIREHAIR_V2_BENCHMARK",
                     "Wh2NativeCodec.cpp", "-flto", "-fprofile")),
            "duplicate/instrumented/changed runtime build")
    deps = checked(["/usr/bin/ninja", "-C", str(build), "-t", "deps"], deadline)
    project = set(SOURCES)
    installed = set()
    generated = build / "Wh2ThueMorseNativeDataR0.h"
    for name in dependency_paths(deps):
        path = Path(name)
        if not path.is_absolute():
            path = build / path
        path = path.resolve(strict=True)
        if path == generated:
            continue
        try:
            relative = path.relative_to(ROOT)
        except ValueError:
            require(str(path).startswith(("/usr/", "/lib/", "/lib64/")),
                    "unreceipted non-system generated dependency")
            installed.add(str(path))
        else:
            project.add(str(relative))
    require("gf256.cpp" in project and "bench/Wh2ThueMorseNativeR0.cpp" in project,
            "incomplete native compiler dependency closure")
    names = sorted(project)
    commit = C.git("rev-parse", "HEAD", deadline=deadline)
    verify_compile_commands(command_text, commit)
    require(commit == C.git("rev-parse", "@{upstream}", deadline=deadline),
            "native source must already be pushed")
    require(not C.git("status", "--porcelain", "--", *names, deadline=deadline),
            "native sources must be clean")
    require(sorted(C.git("ls-files", "--", *names, deadline=deadline).splitlines()) == names,
            "native sources must be tracked")
    source_hashes = {}
    for name in names:
        data = C.read_regular(ROOT / name, 4 * 1024 * 1024, deadline=deadline)
        require(data == C.git_bytes("show", commit + ":" + name, deadline=deadline,
                                    cap=4 * 1024 * 1024), "source blob mismatch: " + name)
        source_hashes[name] = C.sha(data)
    executable = build / TARGET
    symbols = checked(["/usr/bin/nm", "--defined-only", str(executable)], deadline)
    verify_symbols(symbols)
    runtimes = runtime_paths(checked(["/usr/bin/ldd", str(executable)], deadline))
    for path in (generated, build / "Wh2ThueMorseNativeDataR0.json"):
        info = path.lstat()
        require(stat.S_ISREG(info.st_mode) and stat.S_IMODE(info.st_mode) == 0o400
                and info.st_nlink == 1, "generated evidence must be regular immutable single-link")
    data_receipt_bytes = C.read_regular(build / "Wh2ThueMorseNativeDataR0.json",
                                        MAX_RECEIPT, deadline=deadline)
    data_receipt = C.strict_json(data_receipt_bytes)
    require(data_receipt_bytes == C.canonical(data_receipt), "noncanonical generated receipt")
    header_pin = pin(generated, deadline, cap=8 * 1024 * 1024)
    require(data_receipt["header_sha256"] == header_pin["sha256"] and
            data_receipt["header_bytes"] == header_pin["bytes"], "generated header pin")
    # Authenticate originals again, without deriving a fresh root or native row.
    D = sibling("_tm_native_receipt_data", "Wh2ThueMorseNativeDataR0.py")
    authenticated = D.authenticated_inputs(deadline=deadline)
    data = D.extract_data(authenticated)
    require(data_receipt["data_sha256"] == C.sha(C.canonical(data)),
            "generated corpus differs from authenticated inputs")
    require(C.canonical(data_receipt["lookup"]) == C.canonical(authenticated["evidence"]["lookup"])
            and C.canonical(data_receipt["historical_provenance"]) ==
            C.canonical(authenticated["inputs"]["provenance"]), "generated lookup/historical evidence")
    # Derive only the fixed build tables, never rows, ranks, roots or a codec.
    packed, lookup = D.build_lookup(data["pair"])
    require(C.canonical(lookup) == C.canonical(data_receipt["lookup"]), "all generated table pins")
    require(C.sha(D.render_header(data, packed)) == header_pin["sha256"],
            "generated header is not the exact deterministic input projection")
    _, provenance = D.recovery_bundle(deadline=deadline)
    expected_data_receipt = {
        "protocol": PROTOCOL, "recovery_provenance": provenance,
        "historical_provenance": authenticated["inputs"]["provenance"], "lookup": lookup,
        "data_sha256": C.sha(C.canonical(data)), "header_sha256": header_pin["sha256"],
        "header_bytes": header_pin["bytes"],
        "counts": {"traces": 6162, "history": 82, "fixtures": 33, "history_cases": 105,
                   "fixture_cases": 99, "partial_cases": 90, "coefficient_visits": 62347},
    }
    require(C.canonical(data_receipt) == C.canonical(expected_data_receipt),
            "generated receipt schema/provenance/counts mismatch")
    C.time_left(deadline)
    interpreter = pin(Path(sys.executable).resolve(strict=True), deadline, installed=True)
    require(sys.version_info[:2] in ((3, 8), (3, 12)), "unvalidated Python version")
    return {
        "protocol": PROTOCOL, "source_commit": commit, "sources": source_hashes,
        "build_directory": str(build), "cache_sha256": C.sha(cache_bytes),
        "cache": {key: cache[key] for key in sorted(CACHE_EXPECTED)},
        "compiler": compiler, "commands": command_text, "commands_sha256": C.sha(commands),
        "dependencies_sha256": C.sha(deps),
        "system_headers": [pin(path, deadline, installed=True, cap=8 * 1024 * 1024)
                           for path in sorted(installed)],
        "build_files": [pin(build / name, deadline) for name in
                        (TARGET, "libwirehair.a", "build.ninja", "CMakeFiles/rules.ninja")],
        "generated_header": header_pin, "generated_receipt": data_receipt,
        "symbols_sha256": C.sha(symbols),
        "runtime_libraries": [pin(path, deadline, installed=True) for path in runtimes],
        "tools": [pin(Path(path).resolve(strict=True), deadline, installed=True)
                  for path in ("/usr/bin/ninja", "/usr/bin/nm", "/usr/bin/ldd", "/usr/bin/cmake")],
        "interpreter": interpreter, "python_version": list(sys.version_info[:3]),
        "worker_argv": [str(executable), "--worker"], "environment": C.ENV,
        "target_cpu": 50, "sibling_cpu": 114, "worker_seconds": 120,
        "outer_seconds": 150, "observer_seconds": 165,
        "max_stdout_bytes": MAX_STDOUT, "max_stderr_bytes": MAX_STDERR,
        "max_total_bytes": MAX_TOTAL,
    }


def run(receipt_path, started):
    receipt_bytes = C.read_regular(receipt_path, MAX_RECEIPT, deadline=started + 12)
    receipt = C.strict_json(receipt_bytes)
    require(receipt_bytes == C.canonical(receipt), "noncanonical launch receipt")
    require(receipt_bytes == C.canonical(current_receipt(receipt["build_directory"], started + 12)),
            "launch receipt changed")
    OUTPUT.mkdir(mode=0o700, parents=False, exist_ok=False)
    members, total = {}, 0

    def publish(name, data):
        nonlocal total
        require(total + len(data) <= MAX_TOTAL, "aggregate evidence byte cap")
        C.write_new(OUTPUT / name, data)
        total += len(data)
        members[name] = {"bytes": len(data), "sha256": C.sha(data)}

    publish("CLAIM.json", C.canonical({"protocol": PROTOCOL, "receipt": receipt,
                                     "receipt_sha256": C.sha(receipt_bytes)}))
    stdout, stderr, observation = b"", b"", {}
    outcome, failure, analysis = "INVALID", None, None
    try:
        worker_start = time.monotonic()
        require(worker_start - started < 14, "preflight time budget")
        observation = C.capture_worker(receipt["worker_argv"],
                                       min(worker_start + 120, started + 134),
                                       stdout_cap=MAX_STDOUT, stderr_cap=MAX_STDERR)
        stdout, stderr = observation.pop("stdout"), observation.pop("stderr")
        require(not observation["failure"] and observation["returncode"] == 0 and not stderr,
                "native worker failed, timed out or emitted diagnostics")
        require(observation["elapsed_seconds"] <= 120, "native worker wall-time bound")
        raw = C.strict_json(stdout)
        require(stdout.endswith(b"\n") and raw["protocol"] == PROTOCOL,
                "native raw framing/protocol")
        R = sibling("_tm_native_reducer", "Wh2ThueMorseNativeReduceR0.py")
        analysis = R.reduce_raw(raw, deadline=started + 145)
        outcome = analysis["outcome"]
    except Exception as error:
        failure = type(error).__name__ + ": " + str(error)[:1000]
    try:
        require(C.canonical(current_receipt(receipt["build_directory"], started + 148)) == receipt_bytes,
                "native pre/post source/build/input/runtime identity changed")
        observation["post_pins"] = "match"
    except Exception as error:
        observation["post_pins"] = "failed_or_unobserved"
        failure = (failure + "; " if failure else "") + type(error).__name__ + ": " + str(error)[:1000]
        outcome = "INVALID"
    publish("raw.json", stdout)
    publish("stderr.txt", stderr)
    if analysis is not None:
        publish("analysis.json", C.canonical(analysis))
    if time.monotonic() - started >= 149:
        outcome, failure = "INVALID", "outer publication deadline"
    summary = {"protocol": PROTOCOL, "outcome": outcome, "failure": failure,
               "worker": observation, "raw_sha256": C.sha(stdout),
               "elapsed_before_publication": time.monotonic() - started,
               "WH1_compared": analysis is not None,
               "promotion_claimed": False, "all_K_claimed": False,
               "public_candidate_api_claimed": False}
    publish("summary.json", C.canonical(summary))
    manifest = C.canonical({"protocol": PROTOCOL, "outcome": outcome, "files": members})
    require(total + len(manifest) <= MAX_TOTAL, "terminal manifest byte cap")
    C.write_new(OUTPUT / "COMPLETE.json", manifest)
    directory_fd = os.open(str(OUTPUT), os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)
    require(time.monotonic() - started < 150, "whole CLI wall-time bound")
    print(C.canonical(summary).decode("ascii"))
    return 0 if outcome != "INVALID" else 1


def main():
    started = time.monotonic()
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--make-receipt", metavar="PATH")
    mode.add_argument("--run", action="store_true")
    parser.add_argument("--build-dir")
    parser.add_argument("--receipt")
    args = parser.parse_args()
    if args.make_receipt:
        if not args.build_dir or args.receipt:
            parser.error("--make-receipt requires --build-dir and forbids --receipt")
        C.write_new(Path(args.make_receipt), C.canonical(current_receipt(args.build_dir)))
        return 0
    if not args.receipt or args.build_dir:
        parser.error("--run requires --receipt and forbids --build-dir")
    return run(Path(args.receipt), started)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print(type(error).__name__ + ": " + str(error)[:1000], file=sys.stderr)
        sys.exit(1)
