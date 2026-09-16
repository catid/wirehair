#!/usr/bin/env python3
"""Explicit, read-only qualifier for a fresh K4 neutral build.

The R0 K4 readers deliberately reject fresh provenance.  This module is an
isolated reader for a *new* ``Wh2K4FreshNeutralR0.py`` result.  It does not
compile, run tests, launch a worker, or consult the historical R0 receipts.
Every input is checked before the proof is written.
"""
import argparse
import hashlib
import json
import os
import stat
from pathlib import Path
import re
import shlex
import subprocess


ROOT = Path(__file__).resolve().parent.parent
MODES = ("native", "scalar", "asan")
PROTOCOL = "wirehair.wh2.k4-fresh-qualifier-r1"
NEUTRAL_PROTOCOL = "wh2-k4-fresh-neutral-r0"
PRODUCERS = (
    "wirehair.cpp", "codec/WirehairSmall.cpp", "codec/WirehairSmallK5.cpp",
    "codec/WirehairSmallK8.cpp", "codec/WirehairK6.cpp",
    "codec/WirehairK6Core.cpp", "gf256.cpp", "WirehairCodec.cpp",
    "WirehairTools.cpp", "codec/WirehairV2Codec.cpp",
    "codec/WirehairV2Peel.cpp", "codec/WirehairV2Plan.cpp",
    "codec/WirehairV2Policy.cpp", "codec/WirehairV2Precode.cpp",
    "codec/WirehairV2PrecodeDecode.cpp",
    "codec/WirehairV2PrecodeEncode.cpp", "codec/WirehairV2Profile.cpp",
    "codec/WirehairV2Seeds.cpp", "codec/WirehairV2Solve.cpp")
TARGETS = (
    "small_codec_test", "v2_small_codec_test", "v2_small_k5_codec_test",
    "v2_small_k8_codec_test", "small_c_consumer", "v2_small_c_consumer",
    "v2_small_k5_c_consumer", "v2_small_k8_c_consumer", "k6_codec_test",
    "k6_c_consumer", "k6_payload_test", "v2_borrowed_facade_fault_test",
    "gf256_inplace_test", "portability_roundtrip_test")
SHARED = (
    "small_c_consumer_shared", "v2_small_c_consumer_shared",
    "v2_small_k5_c_consumer_shared", "v2_small_k8_c_consumer_shared",
    "k6_c_consumer_shared")
FIXTURE_SHA = "608bafbe37ac0ba3aa94f5d030390af6cc623f9cf0791e86c5bcb8d72f277763"
SOURCE_SUFFIXES = (".cpp", ".c", ".h", ".hpp", ".inc", ".cmake", ".py")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def exact(actual, expected, message):
    require(actual == expected, message)


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                        allow_nan=False) + "\n").encode()


def _regular(path, cap=256 * 1024 ** 2, installed=False):
    path = Path(path)
    try:
        descriptor = os.open(str(path), os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    except OSError as error:
        raise ValueError("regular file: " + str(path)) from error
    try:
        before = os.fstat(descriptor)
        require(stat.S_ISREG(before.st_mode), "regular file: " + str(path))
        require(installed or before.st_nlink == 1, "hard-linked file: " + str(path))
        require(before.st_size <= cap, "file size cap: " + str(path))
        chunks, size = [], 0
        while True:
            chunk = os.read(descriptor, min(65536, cap + 1 - size))
            if not chunk:
                break
            chunks.append(chunk)
            size += len(chunk)
            require(size <= cap, "file grew past cap: " + str(path))
        after = os.fstat(descriptor)
        exact((before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns,
               before.st_ctime_ns, before.st_nlink),
              (after.st_dev, after.st_ino, size, after.st_mtime_ns,
               after.st_ctime_ns, after.st_nlink), "file changed during read")
        return b"".join(chunks)
    finally:
        os.close(descriptor)


def pin(path, cap=256 * 1024 ** 2, installed=False):
    path = Path(path)
    raw = _regular(path, cap, installed=installed)
    return {"path": str(path), "bytes": len(raw),
            "sha256": hashlib.sha256(raw).hexdigest()}


def decode(raw, label="JSON"):
    try:
        value = json.loads(raw.decode())
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("invalid " + label) from error
    return value


def read_json(path, cap=16 * 1024 ** 2):
    return decode(_regular(path, cap), str(path))


def pin_map(records):
    require(isinstance(records, list), "pin list")
    result = {}
    for record in records:
        exact(set(record), {"path", "bytes", "sha256"}, "pin schema")
        path = Path(record["path"])
        require(path.is_absolute() and ".." not in path.parts, "absolute pin path")
        require(type(record["bytes"]) is int and 0 <= record["bytes"] <= 256 * 1024 ** 2,
                "pin byte count")
        require(type(record["sha256"]) is str and
                re.fullmatch(r"[0-9a-f]{64}", record["sha256"]) is not None,
                "pin digest")
        old = result.get(str(path))
        if old is not None:
            exact(old, record, "conflicting duplicate pin")
        result[str(path)] = record
    return result


def _resolve(path):
    path = Path(path)
    require(path.is_absolute(), "absolute path")
    return path.resolve(strict=True)


def _git(args):
    result = subprocess.run(["/usr/bin/git", *map(str, args)], cwd=ROOT,
                            env={"PATH": "/usr/bin:/bin", "LANG": "C",
                                 "LC_ALL": "C", "TZ": "UTC"},
                            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, timeout=60)
    exact(result.returncode, 0, "git inspection")
    exact(result.stderr, b"", "git inspection stderr")
    return result.stdout


def current_head():
    return _git(["rev-parse", "HEAD"]).decode().strip()


def current_source_paths():
    names = _git(["ls-files", "-z"]).decode().split("\0")
    require(names and names[-1] == "", "complete git source roster")
    return tuple(name for name in names[:-1]
                 if Path(name).suffix in SOURCE_SUFFIXES or
                 Path(name).name == "CMakeLists.txt")


def _source_pins(head):
    paths = [ROOT / name for name in current_source_paths()]
    records = [pin(path) for path in sorted(set(paths))]
    return records


def _validate_source(mode_root, result, source_doc, head):
    exact(set(source_doc), {"head", "source_pins"}, "SOURCE schema")
    exact(source_doc["head"], head, "SOURCE current HEAD")
    exact(result["head"], head, "RESULT current HEAD")
    expected = _source_pins(head)
    exact(source_doc["source_pins"], expected, "source pin ordering/roster")
    exact(result["source_pins"], source_doc["source_pins"], "source pin handoff")
    source_map = pin_map(source_doc["source_pins"])
    result_map = pin_map(result["source_pins"])
    expected_map = pin_map(expected)
    exact(source_map, expected_map, "complete current source roster")
    exact(result_map, source_map, "source pins preserved through handoff")
    for record in source_map.values():
        path = Path(record["path"])
        exact(pin(path), record, "source pin changed: " + str(path))
    return source_map


def _validate_artifacts(mode_root, result):
    """Recompute the complete neutral output roster before trusting RESULT."""
    artifact_map = pin_map(result["artifacts"])
    actual = {}
    for path in mode_root.rglob("*"):
        if not path.is_file() or path.is_symlink() or path.name == "RESULT.json":
            continue
        actual[str(path)] = pin(path)
    exact(set(artifact_map), set(actual), "complete neutral artifact roster")
    for name, record in actual.items():
        exact(artifact_map[name], record, "neutral artifact changed: " + name)
    return artifact_map


def _neutral_root(mode_root):
    mode_root = Path(mode_root)
    require(mode_root.is_absolute() and mode_root.name in MODES,
            "absolute mode root")
    require(not mode_root.is_symlink() and mode_root == mode_root.resolve(strict=True),
            "resolved mode root")
    require(mode_root.parent.name.startswith("wh2-k4-fresh-neutral-r0."),
            "fresh neutral root")
    require(mode_root.parent.parent == Path("/var/tmp"),
            "durable fresh neutral root")
    return mode_root


def _command_files(mode_root, index):
    stem = mode_root / ("command-%03d" % index)
    return stem.with_suffix(".argv.json"), stem.with_suffix(".stdout"), \
        stem.with_suffix(".stderr"), stem.with_suffix(".result.json")


def _validate_command_log(mode_root, commands):
    require(type(commands) is list and commands, "complete command log")
    records = []
    for index, invocation in enumerate(commands):
        exact(set(invocation), {"argv", "cwd", "environment"},
              "neutral command schema")
        require(type(invocation["argv"]) is list and invocation["argv"] and
                all(type(arg) is str for arg in invocation["argv"]),
                "command argv")
        cwd = Path(invocation["cwd"])
        require(cwd.is_absolute() and cwd == cwd.resolve(strict=True), "command cwd")
        exact(invocation["environment"],
              {"PATH": "/usr/bin:/bin", "LANG": "C", "LC_ALL": "C", "TZ": "UTC",
               "ASAN_OPTIONS": "detect_leaks=1:detect_stack_use_after_return=1:halt_on_error=1",
               "UBSAN_OPTIONS": "halt_on_error=1:print_stacktrace=1"},
              "neutral command environment")
        argv_path, stdout_path, stderr_path, result_path = _command_files(mode_root, index)
        exact(read_json(argv_path), invocation, "command argv binding")
        stdout, stderr = _regular(stdout_path, 32 * 1024 ** 2), _regular(stderr_path, 1024 ** 2)
        command_result = read_json(result_path, 65536)
        exact(set(command_result), {"returncode"}, "command result schema")
        exact(command_result["returncode"], 0, "successful neutral command")
        records.append({"invocation": invocation, "stdout": pin(stdout_path),
                        "stderr": pin(stderr_path), "result": pin(result_path)})
    return records


def _validate_runner_roster(mode_root, commands, mode):
    expected_tests = len(TARGETS) + (len(SHARED) if mode == "native" else 0)
    exact(len(commands), 9, "exact fresh-neutral command count")
    argv = [entry["argv"] for entry in commands]
    exact(argv[0][:2], ["/usr/bin/git", "rev-parse"], "HEAD command")
    exact(argv[0][2:], ["HEAD"], "HEAD command arguments")
    exact(argv[1], ["/usr/bin/git", "ls-files", "-z"], "source roster command")
    library = mode_root / "library"
    boundary = mode_root / "k4"
    flags = {"native": "", "scalar": "-DANDROID",
             "asan": "-fsanitize=address,undefined -fno-omit-frame-pointer -march=native"}[mode]
    options = ["-G", "Ninja", "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
               "-DCMAKE_BUILD_TYPE=" + ("Debug" if mode == "asan" else "Release"),
               "-DCMAKE_C_FLAGS=" + flags, "-DCMAKE_CXX_FLAGS=" + flags,
               "-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=OFF", "-DBUILD_TESTS=ON",
               "-DBUILD_CODEC_V2=OFF", "-DMARCH_NATIVE=OFF",
               "-DWIREHAIR_STRICT_WARNINGS=ON", "-DWH_LTO=OFF", "-DWH_PGO_MODE=OFF",
               "-DWIREHAIR_BUILD_BOTH=" + ("ON" if mode == "native" else "OFF")]
    exact(argv[2], ["/usr/bin/cmake", "-S", str(ROOT), "-B", str(library)] + options,
          "exact library configure command")
    targets = list(TARGETS) + (list(SHARED) if mode == "native" else [])
    exact(argv[3], ["/usr/bin/ninja", "-C", str(library), "-d", "keepdepfile", "-j", "8",
                    "wirehair"] + targets, "exact library build command")
    tests = "^(" + "|".join(targets) + ")$"
    exact(argv[4], ["/usr/bin/ctest", "--test-dir", str(library), "--output-on-failure",
                    "-j", "4", "-R", tests], "exact library test command")
    # Boundary CMake is a separate project and has no build-both option.
    boundary_options = ["-G", "Ninja", "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
                        "-DCMAKE_BUILD_TYPE=" + ("Debug" if mode == "asan" else "Release"),
                        "-DCMAKE_C_FLAGS=" + flags, "-DCMAKE_CXX_FLAGS=" + flags,
                        "-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=OFF",
                        "-DWH2_SMALL_TEST_DIMENSION=4",
                        "-DWH2_SMALL_LIBRARY=" + str(library / "libwirehair.a"),
                        "-DWH2_SMALL_PORTABLE=" + ("ON" if mode == "scalar" else "OFF"),
                        "-DPython3_EXECUTABLE=/usr/bin/python3"]
    exact(argv[5], ["/usr/bin/cmake", "-S", str(ROOT / "bench" / "Wh2SmallNative"),
                    "-B", str(boundary)] + boundary_options, "exact boundary configure command")
    exact(argv[6], ["/usr/bin/ninja", "-C", str(boundary), "-d", "keepdepfile", "-j", "8"],
          "exact boundary build command")
    exact(argv[7], ["/usr/bin/ctest", "--test-dir", str(boundary), "--output-on-failure",
                    "-j", "4"], "exact boundary test command")
    exact(argv[8][:2], ["/usr/bin/git", "rev-parse"], "final HEAD command")
    exact(argv[8][2:], ["HEAD"], "final HEAD command arguments")
    # Test output is bound to the command result, not merely to RESULT fields.
    lib_text = _regular(mode_root / "command-004.stdout", 32 * 1024 ** 2)
    boundary_text = _regular(mode_root / "command-007.stdout", 32 * 1024 ** 2)
    require((f"100% tests passed, 0 tests failed out of {expected_tests}").encode()
            in lib_text, "exact selected library test count")
    require(b"100% tests passed, 0 tests failed out of 7" in boundary_text,
            "exact selected K4 test count")


def _dep_paths(raw, target, base=None):
    text = raw.decode().replace("\\\n", "")
    target = Path(target)
    prefixes = [str(target) + ": "]
    if base is not None:
        base = Path(base)
        prefixes.append(str(target.relative_to(base)) + ": ")
    prefix = next((value for value in prefixes if text.startswith(value)), None)
    require(prefix is not None, "exact depfile target")
    values = shlex.split(text[len(prefix):])
    require(values and all(Path(value).is_absolute() for value in values),
            "absolute depfile paths")
    paths = {_resolve(target)} | {_resolve(value) for value in values}
    return paths


def _compiler_flags(mode):
    base = ["-DWIREHAIR_BUILDING=1"]
    if mode != "native":
        base.append("-DWIREHAIR_STATIC=1")
    base += ["-I" + str(ROOT / "include")]
    if mode == "scalar":
        base.append("-DANDROID")
    if mode == "asan":
        base += ["-fsanitize=address,undefined", "-fno-omit-frame-pointer",
                 "-march=native", "-g"]
    else:
        base += ["-O3", "-DNDEBUG"]
    base += ["-std=gnu++11", "-fPIC", "-Wall", "-Wextra", "-Wpedantic", "-Werror"]
    return base


def _validate_compile_command(entry, source, object_path, library, mode, prefix):
    exact(set(entry), {"directory", "command", "file", "output"},
          "compile database entry schema")
    exact(entry["directory"], str(library), "compiler cwd")
    exact(entry["file"], str(source), "compiler source binding")
    expected_output = prefix + source.relative_to(ROOT).as_posix() + ".o"
    exact(entry["output"], expected_output, "compiler output binding")
    tokens = shlex.split(entry["command"])
    expected = ["/usr/bin/c++"] + _compiler_flags(mode) + ["-o", expected_output,
                                                              "-c", str(source)]
    exact(tokens, expected, "exact compiler recipe")
    exact(object_path, library / expected_output, "object path")


def _inspect(argv, cap=4 * 1024 ** 2):
    result = subprocess.run(list(map(str, argv)), cwd=ROOT,
                            env={"PATH": "/usr/bin:/bin", "LANG": "C",
                                 "LC_ALL": "C", "TZ": "UTC"},
                            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, timeout=60)
    exact(result.returncode, 0, "tool inspection: " + str(argv))
    exact(result.stderr, b"", "tool inspection stderr: " + str(argv))
    require(len(result.stdout) <= cap, "tool inspection output cap")
    return result.stdout


def validate_toolchain():
    names = ("cmake", "ninja", "ctest", "python3", "c++", "ar", "ranlib", "nm",
             "ld", "as", "ldd", "git")
    tools = {}
    inspections = []
    for name in names:
        path = _resolve(Path("/usr/bin") / name)
        raw = _inspect([str(path), "--version"])
        tools[name] = dict(binary=pin(path, installed=True), version_sha256=hashlib.sha256(raw).hexdigest())
        inspections.append({"argv": [str(path), "--version"], "returncode": 0,
                            "stdout": raw.hex(), "stderr": ""})
    compiler = Path("/usr/bin/c++")
    compiler_commands = {}
    for option in ("-dumpfullversion", "-dumpmachine", "-v"):
        # -v writes the useful identity to stderr on GCC; retain both streams.
        result = subprocess.run([str(compiler), option], cwd=ROOT,
                                env={"PATH": "/usr/bin:/bin", "LANG": "C",
                                     "LC_ALL": "C", "TZ": "UTC"},
                                stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, timeout=60)
        exact(result.returncode, 0, "compiler identity")
        require(len(result.stdout) + len(result.stderr) <= 4 * 1024 ** 2,
                "compiler identity cap")
        compiler_commands[option] = {"stdout": result.stdout.hex(),
                                     "stderr": result.stderr.hex()}
        inspections.append({"argv": [str(compiler), option], "returncode": 0,
                            "stdout": result.stdout.hex(), "stderr": result.stderr.hex()})
    programs = {}
    for name in ("cc1", "cc1plus", "collect2"):
        raw_bytes = _inspect([str(compiler), "-print-prog-name=" + name])
        raw = raw_bytes.decode().strip()
        programs[name] = {"binary": pin(_resolve(Path(raw)), installed=True)}
        inspections.append({"argv": [str(compiler), "-print-prog-name=" + name],
                            "returncode": 0, "stdout": raw_bytes.hex(), "stderr": ""})
    link_inputs = {}
    for name in ("crt1.o", "crti.o", "crtbegin.o", "crtend.o", "crtn.o",
                 "libstdc++.so", "libm.so", "libmvec.so", "libgcc.a",
                 "libgcc_s.so", "libc.so", "libpthread.a", "libasan_preinit.o",
                 "libasan.so", "libubsan.so"):
        raw_bytes = _inspect([str(compiler), "-print-file-name=" + name])
        raw = raw_bytes.decode().strip()
        require(Path(raw).is_absolute(), "absolute compiler link input")
        link_inputs[name] = pin(_resolve(Path(raw)), installed=True)
        inspections.append({"argv": [str(compiler), "-print-file-name=" + name],
                            "returncode": 0, "stdout": raw_bytes.hex(), "stderr": ""})
    for name, record in programs.items():
        binary = Path(record["binary"]["path"])
        raw = _inspect(["/usr/bin/ldd", binary])
        require(b"not found" not in raw, "resolved compiler runtime: " + name)
        inspections.append({"argv": ["/usr/bin/ldd", str(binary)], "returncode": 0,
                            "stdout": raw.hex(), "stderr": ""})
    for name, record in tools.items():
        binary = Path(record["binary"]["path"])
        if _regular(binary, installed=True)[:4] != b"\x7fELF":
            continue
        raw = _inspect(["/usr/bin/ldd", binary])
        require(b"not found" not in raw, "resolved tool runtime: " + name)
        inspections.append({"argv": ["/usr/bin/ldd", str(binary)], "returncode": 0,
                            "stdout": raw.hex(), "stderr": ""})
    return {"tools": tools, "compiler": compiler_commands,
            "programs": programs, "link_inputs": link_inputs,
            "inspections": inspections}


def _production(mode_root, mode, source_map):
    library = mode_root / "library"
    db_path = library / "compile_commands.json"
    database = read_json(db_path, 8 * 1024 ** 2)
    require(type(database) is list, "compile database list")
    prefix = "CMakeFiles/wirehair_objects.dir/" if mode == "native" else "CMakeFiles/wirehair.dir/"
    database = [entry for entry in database
                if str(entry.get("output", "")).startswith(prefix)]
    exact(len(database), 19, "complete production compile database")
    expected_sources = [ROOT / name for name in PRODUCERS]
    entries = {entry.get("file"): entry for entry in database}
    exact(len(entries), 19, "unique production compile entries")
    exact(set(entries), {str(path) for path in expected_sources},
          "exact production source set")
    objects = []
    dependencies = {}
    for source in expected_sources:
        entry = entries[str(source)]
        output = library / (prefix + source.relative_to(ROOT).as_posix() + ".o")
        _validate_compile_command(entry, source, output, library, mode, prefix)
        dep = Path(str(output) + ".d")
        object_pin = pin(output, 64 * 1024 ** 2)
        dep_pin = pin(dep, 2 * 1024 ** 2)
        deps = _dep_paths(_regular(dep, 2 * 1024 ** 2), output, library)
        require(source.resolve() in deps, "producer source depfile binding")
        for path in deps:
            record = pin(path)
            if ROOT in path.parents:
                exact(source_map.get(str(path)), record,
                      "repository dependency missing from source handoff")
            old = dependencies.get(str(path))
            if old is not None:
                exact(old, record, "dependency changed")
            dependencies[str(path)] = record
        objects.append({"source": str(source), "object": object_pin,
                        "depfile": dep_pin, "dependencies":
                        [dependencies[str(path)] for path in sorted(deps)]})
    archive = library / "libwirehair.a"
    archive_pin = pin(archive, 128 * 1024 ** 2)
    listing = _inspect(["/usr/bin/ar", "t", archive]).decode().splitlines()
    expected_members = [Path(item["object"]["path"]).name for item in objects]
    exact(listing, expected_members, "exact production archive order")
    for item in objects:
        member = _inspect(["/usr/bin/ar", "p", archive, Path(item["object"]["path"]).name],
                          cap=64 * 1024 ** 2)
        exact({"path": item["object"]["path"], "bytes": len(member),
               "sha256": hashlib.sha256(member).hexdigest()}, item["object"],
              "archive member/object identity")
    return {"compile_database": pin(db_path), "archive": archive_pin,
            "objects": objects, "dependencies": list(dependencies.values())}


def _boundary(mode_root, mode, source_map):
    boundary = mode_root / "k4"
    db_path = boundary / "compile_commands.json"
    database = read_json(db_path, 2 * 1024 ** 2)
    require(type(database) is list, "boundary compile database list")
    exact(len(database), 8, "exact K4 boundary compile database")
    selected = [entry for entry in database
                if str(entry.get("file", "")) == str(ROOT / "bench/Wh2SmallSerialized.cpp")]
    exact(len(selected), 1, "one K4 boundary producer")
    entry = selected[0]
    exact(set(entry), {"directory", "command", "file", "output"},
          "boundary compile schema")
    exact(entry["directory"], str(boundary), "boundary compiler cwd")
    source = ROOT / "bench/Wh2SmallSerialized.cpp"
    expected_output = "CMakeFiles/wh2_small_serialized.dir" + str(source) + ".o"
    exact(entry["output"], expected_output, "boundary output binding")
    tokens = shlex.split(entry["command"])
    boundary_flags = (["-DANDROID"] if mode == "scalar" else []) + ["-DWH2_SMALL_CODEC_K=4"]
    if mode == "scalar":
        boundary_flags.append("-DWH2_SMALL_EXPECT_PORTABLE=1")
    boundary_flags += ["-I" + str(ROOT), "-I" + str(ROOT / "include"),
                       "-I" + str(boundary)]
    if mode == "scalar":
        boundary_flags.append("-DANDROID")
    boundary_flags += (["-fsanitize=address,undefined", "-fno-omit-frame-pointer",
                        "-march=native", "-g"] if mode == "asan" else
                       ["-O3", "-DNDEBUG"])
    boundary_flags += ["-std=c++11", "-fPIC", "-Wall", "-Wextra", "-Wpedantic",
                       "-Werror", "-fno-strict-aliasing", "-fno-lto"]
    exact(tokens, ["/usr/bin/c++"] + boundary_flags + ["-o", expected_output,
                                                         "-c", str(source)],
          "exact K4 compiler recipe")
    object_path = boundary / expected_output
    object_pin = pin(object_path, 64 * 1024 ** 2)
    dep = Path(str(object_path) + ".d")
    dep_pin = pin(dep, 2 * 1024 ** 2)
    deps = _dep_paths(_regular(dep, 2 * 1024 ** 2), object_path, boundary)
    fixture = boundary / "Wh2K4NativeData.inc"
    exact(pin(fixture)["sha256"], FIXTURE_SHA, "retained K4 fixture")
    require(source.resolve() in deps and fixture.resolve() in deps,
            "boundary source/fixture dependency")
    dep_records = [pin(path) for path in sorted(deps)]
    for path, record in zip(sorted(deps), dep_records):
        if ROOT in path.parents:
            exact(source_map.get(str(path)), record,
                  "boundary dependency missing from source handoff")
    archive = boundary / "libwh2_small_serialized.a"
    archive_pin = pin(archive, 128 * 1024 ** 2)
    listing = _inspect(["/usr/bin/ar", "t", archive]).decode().splitlines()
    exact(listing, [object_path.name], "sole K4 archive member")
    member = _inspect(["/usr/bin/ar", "p", archive, object_path.name], cap=64 * 1024 ** 2)
    exact(hashlib.sha256(member).hexdigest(), object_pin["sha256"],
          "K4 archive/member identity")
    return {"compile_database": pin(db_path), "entry": entry,
            "object": object_pin, "depfile": dep_pin,
            "dependencies": dep_records, "fixture": pin(fixture),
            "archive": archive_pin}


def qualify(mode_root, mode, proof_path=None):
    """Validate one fresh neutral mode and optionally publish its R1 proof."""
    require(mode in MODES, "backend")
    mode_root = _neutral_root(mode_root)
    result = read_json(mode_root / "RESULT.json", 8 * 1024 ** 2)
    exact(set(result), {"schema", "mode", "status", "head", "library_tests",
                       "boundary_tests", "commands", "source_pins", "artifacts",
                       "producing_source_closure", "scientific_launch", "scope"},
          "RESULT schema")
    exact(result["schema"], NEUTRAL_PROTOCOL, "neutral protocol")
    exact((result["mode"], result["status"]), (mode, "PASS"), "neutral status")
    exact((result["library_tests"], result["boundary_tests"]),
          (len(TARGETS) + (len(SHARED) if mode == "native" else 0), 7),
          "neutral test counts")
    exact((result["producing_source_closure"], result["scientific_launch"]),
          (False, False), "neutral-only result")
    exact(result["scope"],
          "fresh neutral checks; independent producing/tool/runtime qualification still required",
          "neutral scope")
    source_doc = read_json(mode_root / "SOURCE.json", 8 * 1024 ** 2)
    head = current_head()
    source_map = _validate_source(mode_root, result, source_doc, head)
    artifact_map = _validate_artifacts(mode_root, result)
    command_records = _validate_command_log(mode_root, result["commands"])
    _validate_runner_roster(mode_root, result["commands"], mode)
    production = _production(mode_root, mode, source_map)
    boundary = _boundary(mode_root, mode, source_map)
    toolchain = validate_toolchain()
    # The complete source pin set is deliberately copied verbatim into the
    # proof, rather than being reconstructed from a reduced dependency set.
    proof = {"protocol": PROTOCOL, "mode": mode, "head": head,
             "neutral_result": pin(mode_root / "RESULT.json"),
             "neutral_source": pin(mode_root / "SOURCE.json"),
             "source_pins": result["source_pins"],
             "artifacts": list(result["artifacts"]),
             "commands": command_records, "production": production,
             "boundary": boundary, "toolchain": toolchain,
             "producing_source_closure": True, "scientific_launch": False,
             "scope": "fresh current-source K4 production and boundary closure"}
    if proof_path is not None:
        proof_path = Path(proof_path)
        parent = proof_path.parent
        require(proof_path.is_absolute() and not proof_path.exists() and
                not proof_path.is_symlink() and parent != mode_root and
                parent != ROOT and ROOT not in parent.parents and
                parent.parent == Path("/var/tmp") and
                parent.name.startswith("wh2-k4-fresh-qualified-r1.") and
                not parent.is_symlink() and (not parent.exists() or parent.is_dir()),
                "fresh proof destination")
        if not parent.exists():
            parent.mkdir(mode=0o700, parents=False, exist_ok=False)
        raw = canonical(proof)
        with proof_path.open("xb") as stream:
            stream.write(raw); stream.flush(); os.fsync(stream.fileno())
        proof_path.chmod(0o400)
    return proof


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=MODES)
    parser.add_argument("mode_root", type=Path)
    parser.add_argument("--proof", type=Path)
    args = parser.parse_args(argv)
    proof = qualify(args.mode_root, args.mode, args.proof)
    print(json.dumps({"protocol": proof["protocol"], "mode": proof["mode"],
                      "head": proof["head"], "producing_source_closure": True},
                     sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
