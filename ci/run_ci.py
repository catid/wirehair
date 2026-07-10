#!/usr/bin/env python3
"""Reusable local entry points for Wirehair's CI lanes."""

import argparse
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
WIREHAIR_C_API_EXPORTS = frozenset({
    "wirehair_decode",
    "wirehair_decoder_becomes_encoder",
    "wirehair_decoder_create",
    "wirehair_decoder_create_ex",
    "wirehair_encode",
    "wirehair_encoder_create",
    "wirehair_encoder_create_ex",
    "wirehair_encoder_create_owned",
    "wirehair_encoder_create_owned_ex",
    "wirehair_free",
    "wirehair_init_",
    "wirehair_recover",
    "wirehair_recover_block",
    "wirehair_recover_block_ex",
    "wirehair_result_string",
})


def display_command(command):
    if os.name == "nt":
        return subprocess.list2cmdline([str(item) for item in command])
    return shlex.join([str(item) for item in command])


def run(command, *, cwd=ROOT, env=None, log_path=None):
    command = [str(item) for item in command]
    print("+ " + display_command(command), flush=True)
    if log_path is None:
        subprocess.run(command, cwd=str(cwd), env=env, check=True)
        return

    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8", newline="") as log:
        process = subprocess.Popen(
            command,
            cwd=str(cwd),
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            errors="replace",
        )
        assert process.stdout is not None
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            log.write(line)
        return_code = process.wait()
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, command)


def run_expect_failure(
        command, *, env, log_path, expected_returncode, required_text,
        forbidden_text=()):
    command = [str(item) for item in command]
    print("+ expect-failure: " + display_command(command), flush=True)
    completed = subprocess.run(
        command,
        cwd=str(ROOT),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        errors="replace",
    )
    sys.stdout.write(completed.stdout)
    sys.stdout.flush()
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(completed.stdout, encoding="utf-8")
    if completed.returncode != expected_returncode:
        raise RuntimeError(
            "negative-path command returned %d instead of %d" %
            (completed.returncode, expected_returncode))
    if required_text not in completed.stdout:
        raise RuntimeError(
            "negative-path command failed without expected diagnostic: " +
            required_text)
    for signature in forbidden_text:
        if signature in completed.stdout:
            raise RuntimeError(
                "negative-path output contained forbidden diagnostic: " +
                signature)


def require_log_line(log_path, expected_line):
    lines = Path(log_path).read_text(encoding="utf-8").splitlines()
    if expected_line not in lines:
        raise RuntimeError(
            "expected exact diagnostic was not found in %s: %s" %
            (log_path, expected_line))


def capture_output(command, *, log_path=None):
    command = [str(item) for item in command]
    print("+ " + display_command(command), flush=True)
    completed = subprocess.run(
        command,
        cwd=str(ROOT),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        errors="replace",
    )
    if log_path is not None:
        log_path = Path(log_path)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(completed.stdout, encoding="utf-8")
    if completed.returncode != 0:
        sys.stdout.write(completed.stdout)
        sys.stdout.flush()
        raise subprocess.CalledProcessError(completed.returncode, command)
    return completed.stdout


def cmake_generator_args(args):
    result = []
    if args.generator:
        result.extend(["-G", args.generator])
    if args.architecture:
        result.extend(["-A", args.architecture])
    return result


def is_msvc_generator(generator):
    return bool(generator and generator.startswith("Visual Studio"))


def configure(args, *, shared, extra_args=(), explicit_tools=False):
    command = [
        "cmake",
        "-S", ROOT,
        "-B", args.build_dir,
        *cmake_generator_args(args),
        "-DBUILD_TESTS=ON",
        "-DBUILD_CODEC_V2=ON",
        "-DMARCH_NATIVE=OFF",
        "-DWH_LTO=OFF",
        "-DWH_PGO_MODE=OFF",
        "-DWIREHAIR_BUILD_TOOLS=" + ("ON" if explicit_tools else "OFF"),
        "-DWIREHAIR_BUILD_BENCHMARKS=" + ("ON" if explicit_tools else "OFF"),
        "-DWIREHAIR_BUILD_BOTH=OFF",
        "-DBUILD_SHARED_LIBS=" + ("ON" if shared else "OFF"),
        "-DCMAKE_BUILD_TYPE=" + args.config,
        *strict_cmake_args(args),
        *extra_args,
        *args.cmake_arg,
    ]
    run(command)


def build(args, targets=()):
    command = [
        "cmake", "--build", args.build_dir,
        "--config", args.config,
        "--parallel", str(args.jobs),
    ]
    if targets:
        command.extend(["--target", *targets])
    run(command)


def ctest(args, *, regex=None, label=None, env=None, prefix=()):
    command = [
        *prefix,
        "ctest",
        "--test-dir", args.build_dir,
        "--build-config", args.config,
        "--output-on-failure",
        "--no-tests=error",
    ]
    if regex:
        command.extend(["--tests-regex", regex])
    if label:
        command.extend(["--label-regex", label])
    run(command, env=env)


def install(args):
    run([
        "cmake", "--install", args.build_dir,
        "--config", args.config,
        "--prefix", args.install_dir,
    ])


def runtime_env(prefix, platform=None):
    platform = platform or sys.platform
    env = os.environ.copy()
    prefix = Path(prefix).resolve()
    if platform.startswith("win"):
        key = "PATH"
        directories = [prefix / "bin", prefix / "lib"]
    elif platform == "darwin":
        key = "DYLD_LIBRARY_PATH"
        directories = [prefix / "lib", prefix / "lib64"]
    else:
        key = "LD_LIBRARY_PATH"
        directories = [prefix / "lib", prefix / "lib64"]
    old_value = env.get(key)
    values = [str(path) for path in directories]
    if old_value:
        values.append(old_value)
    env[key] = os.pathsep.join(values)
    return env


def package_consumer(args):
    source_dir = ROOT / "test" / "package"
    if not (source_dir / "CMakeLists.txt").is_file():
        raise RuntimeError("installed package fixture is missing: " + str(source_dir))

    consumer_build = Path(args.build_dir) / "package-consumer"
    command = [
        "cmake",
        "-S", source_dir,
        "-B", consumer_build,
        *cmake_generator_args(args),
        "-DCMAKE_PREFIX_PATH=" + str(Path(args.install_dir).resolve()),
        "-DCMAKE_BUILD_TYPE=" + args.config,
        *consumer_toolchain_args(args),
    ]
    run(command)
    run([
        "cmake", "--build", consumer_build,
        "--config", args.config,
        "--parallel", str(args.jobs),
    ])
    run([
        "ctest",
        "--test-dir", consumer_build,
        "--build-config", args.config,
        "--output-on-failure",
        "--no-tests=error",
    ], env=runtime_env(args.install_dir))


def strict_cmake_args(args):
    if not getattr(args, "strict", False):
        return []
    return ["-DWIREHAIR_STRICT_WARNINGS=ON"]


def consumer_toolchain_args(args):
    prefixes = ("-DCMAKE_C_COMPILER=", "-DCMAKE_CXX_COMPILER=")
    forwarded = [
        argument for argument in args.cmake_arg
        if argument.startswith(prefixes)
    ]
    return [*forwarded, *strict_cmake_args(args)]


def find_executable(build_dir, name, config="Release"):
    suffix = ".exe" if os.name == "nt" else ""
    build_dir = Path(build_dir)
    candidates = [
        build_dir / (name + suffix),
        build_dir / config / (name + suffix),
        build_dir / "codec" / (name + suffix),
        build_dir / "codec" / config / (name + suffix),
    ]
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    raise RuntimeError("cannot find built executable %r under %s" % (name, build_dir))


def find_python_module(prefix):
    prefix = Path(prefix)
    direct = prefix / "python" / "whirehair.py"
    if direct.is_file():
        return direct.resolve()
    matches = sorted(prefix.rglob("whirehair.py"))
    if not matches:
        raise RuntimeError("installed whirehair.py was not found under " + str(prefix))
    return matches[0].resolve()


def find_shared_library(prefix, platform=None):
    platform = platform or sys.platform
    prefix = Path(prefix)
    if platform.startswith("win"):
        names = ("wirehair.dll", "libwirehair.dll")
    elif platform == "darwin":
        names = ("libwirehair.dylib", "libwirehair.2.dylib")
    else:
        names = ("libwirehair.so", "libwirehair.so.2")

    for directory in (prefix / "bin", prefix / "lib", prefix / "lib64"):
        for name in names:
            candidate = directory / name
            if candidate.is_file():
                return candidate.resolve()
    for name in names:
        matches = sorted(path for path in prefix.rglob(name) if path.is_file())
        if matches:
            return matches[0].resolve()
    raise RuntimeError("installed Wirehair shared library was not found under " + str(prefix))


def _visual_studio_dumpbin_candidates(installation):
    tools_root = Path(installation) / "VC" / "Tools" / "MSVC"

    def version_key(path):
        return tuple(
            (0, int(part)) if part.isdigit() else (1, part)
            for part in re.split(r"[.-]", path.name)
        )

    versions = sorted(
        (path for path in tools_root.glob("*") if path.is_dir()),
        key=version_key,
        reverse=True,
    )
    host_targets = (
        ("Hostx64", "x64"),
        ("Hostx64", "x86"),
        ("Hostx86", "x64"),
        ("Hostx86", "x86"),
    )
    return [
        version / "bin" / host / target / "dumpbin.exe"
        for version in versions
        for host, target in host_targets
    ]


def _query_vswhere_installations(vswhere):
    completed = subprocess.run(
        [
            str(vswhere),
            "-products", "*",
            "-version", "[17.0,18.0)",
            "-requires", "Microsoft.VisualStudio.Component.VC.Tools.x86.x64",
            "-property", "installationPath",
        ],
        cwd=str(ROOT),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        errors="replace",
        check=False,
    )
    if completed.returncode != 0:
        print(
            "vswhere could not locate Visual Studio 2022:\n" + completed.stdout,
            file=sys.stderr,
            flush=True,
        )
        return []
    return [Path(line.strip()) for line in completed.stdout.splitlines() if line.strip()]


def find_dumpbin(environment=None, which=None):
    environment = os.environ if environment is None else environment
    which = shutil.which if which is None else which
    candidates = []

    path_dumpbin = which("dumpbin")
    if path_dumpbin and Path(path_dumpbin).is_file():
        return Path(path_dumpbin).resolve()

    tools_install = environment.get("VCToolsInstallDir")
    if tools_install:
        tools_root = Path(tools_install)
        for host, target in (
                ("Hostx64", "x64"), ("Hostx64", "x86"),
                ("Hostx86", "x64"), ("Hostx86", "x86")):
            candidates.append(tools_root / "bin" / host / target / "dumpbin.exe")

    installations = []
    vs_install = environment.get("VSINSTALLDIR")
    if vs_install:
        installations.append(Path(vs_install))
    vc_install = environment.get("VCINSTALLDIR")
    if vc_install:
        installations.append(Path(vc_install).parent)
    common_tools = environment.get("VS170COMNTOOLS")
    if common_tools:
        installations.append(Path(common_tools).parent.parent)

    vswhere = which("vswhere")
    if not vswhere:
        program_files_x86 = environment.get("ProgramFiles(x86)")
        if program_files_x86:
            standard_vswhere = (
                Path(program_files_x86) / "Microsoft Visual Studio" /
                "Installer" / "vswhere.exe"
            )
            if standard_vswhere.is_file():
                vswhere = str(standard_vswhere)
    if vswhere:
        installations.extend(_query_vswhere_installations(vswhere))

    for variable in ("ProgramFiles", "ProgramFiles(x86)"):
        program_files = environment.get(variable)
        if program_files:
            installations.extend(sorted(
                (Path(program_files) / "Microsoft Visual Studio" / "2022").glob("*"),
                reverse=True,
            ))

    seen_installations = set()
    for installation in installations:
        normalized = str(installation)
        if normalized in seen_installations:
            continue
        seen_installations.add(normalized)
        candidates.extend(_visual_studio_dumpbin_candidates(installation))

    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    raise RuntimeError(
        "dumpbin.exe was not found in PATH, the Visual Studio environment, "
        "vswhere results, or standard Visual Studio 2022 install roots")


def dumpbin_exports_command(dumpbin, library):
    return [str(dumpbin), "/nologo", "/exports", str(library)]


def parse_dumpbin_exports(dumpbin_output):
    exports = set()
    export_line = re.compile(
        r"^\s*\d+\s+[0-9A-Fa-f]+\s+[0-9A-Fa-f]+\s+([^\s=]+)")
    for line in dumpbin_output.splitlines():
        match = export_line.match(line)
        if match:
            exports.add(match.group(1))
    return frozenset(exports)


def validate_msvc_exports(args):
    library = find_shared_library(args.install_dir, platform="win32")
    dumpbin = find_dumpbin()
    log_path = Path(args.build_dir) / "ci-logs" / "msvc-dumpbin-exports.txt"
    output = capture_output(
        dumpbin_exports_command(dumpbin, library),
        log_path=log_path,
    )
    exports = parse_dumpbin_exports(output)
    if exports != WIREHAIR_C_API_EXPORTS:
        raise RuntimeError(
            "MSVC DLL export mismatch: missing=%s extra=%s; full dumpbin "
            "output: %s" %
            (sorted(WIREHAIR_C_API_EXPORTS - exports),
             sorted(exports - WIREHAIR_C_API_EXPORTS), log_path))
    print(
        "MSVC DLL exports verified (%d): %s" %
        (len(exports), ", ".join(sorted(exports))),
        flush=True,
    )


def should_validate_msvc_exports(args):
    return (
        args.linkage == "shared" and args.strict and
        is_msvc_generator(args.generator)
    )


def python_unit_tests():
    run([sys.executable, "-m", "unittest", "-v", "ci.test_run_ci"])
    run([sys.executable, str(ROOT / "python" / "test_whirehair.py"), "-v"])


def python_native_test(args):
    module = find_python_module(args.install_dir)
    library = find_shared_library(args.install_dir)
    env = runtime_env(args.install_dir)
    env["WIREHAIR_LIBRARY"] = str(library)
    run([
        sys.executable,
        str(ROOT / "ci" / "python_native_test.py"),
        "--module-dir", module.parent,
        "--library", library,
    ], env=env)


def large_message_profiles(args, mode):
    bash = shutil.which("bash")
    if not bash:
        raise RuntimeError("bash is required for large-message profile validation")
    executable = find_executable(args.build_dir, "large_message_test", args.config)
    executable_arg = executable.as_posix() if os.name == "nt" else str(executable)
    run([
        bash,
        ROOT / "test" / "run_large_message_profiles.sh",
        executable_arg,
        mode,
    ], env=os.environ.copy())


def explicit_tools_smoke(args):
    logs = Path(args.build_dir) / "ci-logs"
    cases = (
        ("gen_small_dseeds", ["--selection-self-test"]),
        ("gen_peel_seeds", [
            "--trials", "1", "--sublo", "0", "--subhi", "0",
            "--nlo", "2048", "--nhi", "2048",
            "--max-tries", "1", "--skip-tuning",
        ]),
        ("gen_most_dseeds", [
            "--seed", "1", "--trials", "1",
            "--dense-index-lo", "13", "--dense-index-hi", "13",
        ]),
        ("gen_dcounts", [
            "--seed", "1", "--trials", "1", "--nlo", "2", "--nhi", "2",
            "--max-failures", "1", "--low-count-run", "1",
        ]),
        ("gen_tables", ["--no-benchmarks", "--heavy-trials", "0"]),
        ("wirehair_v2_bench", [
            "compare", "--nlo", "2", "--nhi", "2", "--trials", "1",
            "--bb-list", "8", "--max-message-mib", "1", "--loss", "0",
        ]),
        ("wirehair_v2_bench", [
            "precodecheck", "--N", "2", "--bb-list", "8",
            "--trials", "1", "--loss", "0",
        ]),
    )
    for name, arguments in cases:
        executable = find_executable(args.build_dir, name, args.config)
        run([executable, *arguments], log_path=logs / (name + "-smoke.log"))


def run_matrix(args):
    shared = args.linkage == "shared"
    explicit_tools = args.strict and is_msvc_generator(args.generator)
    python_unit_tests()
    configure(args, shared=shared, explicit_tools=explicit_tools)
    build(args)
    if explicit_tools:
        explicit_tools_smoke(args)
    large_message_profiles(args, "quick")
    ctest(args)
    install(args)
    if should_validate_msvc_exports(args):
        validate_msvc_exports(args)
    package_consumer(args)
    if shared:
        python_native_test(args)


def sanitizer_flags(kind):
    sanitize = "address,undefined" if kind == "asan-ubsan" else "thread"
    compile_flags = "-fsanitize=%s -fno-omit-frame-pointer" % sanitize
    link_flags = "-fsanitize=%s" % sanitize
    return [
        "-DCMAKE_C_FLAGS=" + compile_flags,
        "-DCMAKE_CXX_FLAGS=" + compile_flags,
        "-DCMAKE_EXE_LINKER_FLAGS=" + link_flags,
        "-DCMAKE_SHARED_LINKER_FLAGS=" + link_flags,
    ]


def run_sanitizer(args):
    configure(args, shared=False, extra_args=sanitizer_flags("asan-ubsan"))
    build(args)
    large_message_profiles(args, "quick")
    ctest(args, env=os.environ.copy())


def run_tsan(args):
    configure(args, shared=False, extra_args=sanitizer_flags("tsan"))
    build(args, targets=("legacy_core_test",))
    prefix = ()
    if sys.platform.startswith("linux"):
        setarch = shutil.which("setarch")
        if not setarch:
            raise RuntimeError("Linux TSan requires setarch from util-linux")
        # High-entropy ASLR on current Linux hosts can overlap TSan's fixed
        # shadow region before the test starts. The personality is inherited by
        # CTest children, so disabling ASLR here isolates that runtime issue.
        prefix = (setarch, os.uname().machine, "-R")
    ctest(
        args,
        regex=r"^legacy_core_(test|cold_start_test)$",
        env=os.environ.copy(),
        prefix=prefix,
    )


def scenario_env(**values):
    env = os.environ.copy()
    env.update({key: str(value) for key, value in values.items()})
    return env


def run_scheduled(args):
    configure(
        args,
        shared=False,
        extra_args=("-DWIREHAIR_ENABLE_SCHEDULED_TESTS=ON",),
    )
    build(args, targets=(
        "unit_test",
        "large_message_test",
        "wirehair_v2_precode_decode_test",
    ))
    ctest(args, label="^scheduled$")
    large_message_profiles(args, "scheduled")
    logs = Path(args.build_dir) / "ci-logs"
    unit_test = find_executable(args.build_dir, "unit_test", args.config)
    large_test = find_executable(args.build_dir, "large_message_test", args.config)

    common = {
        "WIREHAIR_UNIT_SKIP_BENCH": 1,
        "WIREHAIR_UNIT_SEED": "0x4d595df4d0f33173",
        "OMP_NUM_THREADS": 2,
    }
    scenarios = [
        ("small", 2, 64, 1, 16),
        ("medium", 4096, 4096, 64, 64),
        ("maximum-n", 64000, 64000, 8, 8),
    ]
    for name, min_n, max_n, min_bytes, max_bytes in scenarios:
        env = scenario_env(
            **common,
            WIREHAIR_UNIT_MIN_N=min_n,
            WIREHAIR_UNIT_MAX_N=max_n,
            WIREHAIR_UNIT_BLOCK_BYTES_MIN=min_bytes,
            WIREHAIR_UNIT_BLOCK_BYTES_MAX=max_bytes,
        )
        run([unit_test], env=env, log_path=logs / (name + ".log"))

    terminal_env = scenario_env(
        **common,
        WIREHAIR_UNIT_MIN_N=2,
        WIREHAIR_UNIT_MAX_N=2,
        WIREHAIR_UNIT_BLOCK_BYTES_MIN=1,
        WIREHAIR_UNIT_BLOCK_BYTES_MAX=1,
        WIREHAIR_UNIT_INJECT_TERMINAL=1,
    )
    run_expect_failure(
        [unit_test],
        env=terminal_env,
        log_path=logs / "terminal-failure-sentinel.log",
        expected_returncode=4,
        required_text="Terminal decoder exhaustion",
        forbidden_text=(
            "Sanitizer",
            "runtime error:",
        ),
    )

    large_env = scenario_env(
        WIREHAIR_LARGE_N=128,
        WIREHAIR_LARGE_BLOCK_BYTES=1048576,
        WIREHAIR_LARGE_FINAL_BYTES=1048573,
        WIREHAIR_LARGE_LOSS_COUNT=13,
        WIREHAIR_LARGE_MAX_REPAIR=128,
        WIREHAIR_LARGE_SEED="0x8af09d31e7c4526b",
        WIREHAIR_LARGE_MAX_PAYLOAD_BYTES=268435456,
        WIREHAIR_LARGE_MAX_MILLISECONDS=120000,
    )
    run([large_test], env=large_env, log_path=logs / "large-block.log")

    compiler = shlex.split(os.environ.get("CXX", "c++"))
    whx = Path(args.build_dir).resolve() / ("whx.exe" if os.name == "nt" else "whx")
    run([
        *compiler,
        "-O2", "-std=c++11", "-pthread", "-Wall", "-Wextra",
        "-Wno-unused-function",
        "-I" + str(ROOT),
        "-I" + str(ROOT / "include"),
        "-I" + str(ROOT / "test"),
        ROOT / "gf256.cpp",
        ROOT / "WirehairCodec.cpp",
        ROOT / "WirehairTools.cpp",
        ROOT / "wirehair.cpp",
        ROOT / "bench" / "whx.cpp",
        "-o", whx,
    ])
    run([whx, "selftest"], log_path=logs / "whx-selftest.log")
    run([
        whx, "bench", "--threads", "1",
        "--N", "2,4096,64000", "--bb", "1,257",
        "--rounds", "1", "--memory-mib", "4096", "--loss", "0.30",
    ], log_path=logs / "whx-deterministic-loss.log")
    run([
        whx, "fuzz", "--threads", "1", "--secs", "15",
        "--nmax", "64000", "--seed", "0x71c3a5e9",
    ], log_path=logs / "whx-fuzz.log")


def run_qemu(args):
    if not sys.platform.startswith("linux"):
        raise RuntimeError("the QEMU portability lane is supported on Linux hosts")
    qemu = shutil.which("qemu-x86_64") or shutil.which("qemu-x86_64-static")
    if not qemu:
        raise RuntimeError("qemu-x86_64 is required (install qemu-user or qemu-user-static)")
    compiler = shlex.split(os.environ.get("CXX", "clang++"))
    if not shutil.which(compiler[0]):
        raise RuntimeError("C++ compiler was not found: " + compiler[0])

    output_dir = Path(args.build_dir).resolve()
    logs = output_dir / "ci-logs"
    output_dir.mkdir(parents=True, exist_ok=True)
    sources = [
        ROOT / "test" / "LegacyCoreTest.cpp",
        ROOT / "gf256.cpp",
        ROOT / "WirehairCodec.cpp",
        ROOT / "WirehairTools.cpp",
        ROOT / "wirehair.cpp",
    ]
    common = [
        *compiler,
        "-O2", "-g", "-std=c++11", "-pthread", "-Wall", "-Wextra",
        "-Wno-unused-function", "-DWIREHAIR_TESTING=1",
        "-I" + str(ROOT),
        "-I" + str(ROOT / "include"),
        "-I" + str(ROOT / "test"),
    ]
    cases = [
        (
            "conroe-portable",
            "Conroe",
            [],
            "Active x86 kernels: SSSE3=0 AVX2=0 GFNI=0 AVX512=0",
        ),
        (
            "haswell-avx2",
            "Haswell",
            ["-mavx2"],
            "Active x86 kernels: SSSE3=1 AVX2=1 GFNI=0 AVX512=0",
        ),
    ]
    for name, cpu, flags, expected_kernels in cases:
        executable = output_dir / name
        runtime_log = logs / (name + ".log")
        run(
            [*common, *flags, *sources, "-o", executable],
            log_path=logs / (name + "-compile.log"),
        )
        run(
            [qemu, "-cpu", cpu, executable],
            log_path=runtime_log,
        )
        require_log_line(runtime_log, expected_kernels)


def find_unique_file(root, name):
    matches = sorted(path for path in Path(root).rglob(name) if path.is_file())
    if len(matches) != 1:
        raise RuntimeError(
            "expected one %s under %s, found %d" %
            (name, root, len(matches)))
    return matches[0].resolve()


def parse_mingw_exports(nm_output):
    exports = set()
    for line in nm_output.splitlines():
        fields = line.split()
        if fields and fields[-1].startswith("wirehair_"):
            exports.add(fields[-1])
    return frozenset(exports)


def imports_wirehair_dll(objdump_output):
    return any(
        "dll name:" in line.lower() and "libwirehair.dll" in line.lower()
        for line in objdump_output.splitlines()
    )


def run_mingw(args):
    triplet = "x86_64-w64-mingw32"
    tools = {
        "cc": triplet + "-gcc",
        "cxx": triplet + "-g++",
        "rc": triplet + "-windres",
        "nm": triplet + "-nm",
        "objdump": triplet + "-objdump",
    }
    for name, executable in tools.items():
        if not shutil.which(executable):
            raise RuntimeError("MinGW %s tool was not found: %s" % (name, executable))

    strict_flags = "-Wall -Wextra -Wpedantic -Werror"
    for linkage in ("static", "shared"):
        shared = linkage == "shared"
        build_dir = Path(args.build_dir) / linkage
        install_dir = Path(args.install_dir) / linkage
        consumer_build = build_dir / "package-consumer"
        common = [
            "-DCMAKE_SYSTEM_NAME=Windows",
            "-DCMAKE_C_COMPILER=" + tools["cc"],
            "-DCMAKE_CXX_COMPILER=" + tools["cxx"],
            "-DCMAKE_RC_COMPILER=" + tools["rc"],
            "-DCMAKE_C_FLAGS=" + strict_flags,
            "-DCMAKE_CXX_FLAGS=" + strict_flags,
            "-DCMAKE_BUILD_TYPE=Release",
        ]
        run([
            "cmake", "-S", ROOT, "-B", build_dir,
            *cmake_generator_args(args),
            *common,
            "-DBUILD_SHARED_LIBS=" + ("ON" if shared else "OFF"),
            "-DBUILD_TESTS=OFF",
            "-DBUILD_CODEC_V2=OFF",
            "-DMARCH_NATIVE=OFF",
            "-DWIREHAIR_BUILD_BOTH=OFF",
            "-DWIREHAIR_BUILD_TOOLS=OFF",
            "-DWIREHAIR_BUILD_BENCHMARKS=OFF",
        ])
        run([
            "cmake", "--build", build_dir,
            "--config", "Release", "--parallel", str(args.jobs),
        ])
        run([
            "cmake", "--install", build_dir,
            "--config", "Release", "--prefix", install_dir,
        ])

        package_dir = install_dir / "lib" / "cmake" / "wirehair"
        run([
            "cmake", "-S", ROOT / "test" / "package", "-B", consumer_build,
            *cmake_generator_args(args),
            *common,
            "-Dwirehair_DIR=" + str(package_dir.resolve()),
        ])
        run([
            "cmake", "--build", consumer_build,
            "--config", "Release", "--parallel", str(args.jobs),
        ])

        consumer = find_unique_file(consumer_build, "package_c_consumer.exe")
        logs = build_dir / "ci-logs"
        imports = capture_output(
            [tools["objdump"], "-p", consumer],
            log_path=logs / "package-c-consumer-imports.txt",
        )
        if imports_wirehair_dll(imports) != shared:
            raise RuntimeError(
                "%s MinGW C consumer has the wrong libwirehair.dll import state" %
                linkage)

        if shared:
            find_unique_file(install_dir / "bin", "libwirehair.dll")
            import_library = find_unique_file(
                install_dir / "lib", "libwirehair.dll.a")
            nm_output = capture_output(
                [tools["nm"], "-g", "--defined-only", import_library],
                log_path=logs / "wirehair-import-library-symbols.txt",
            )
            exports = parse_mingw_exports(nm_output)
            if exports != WIREHAIR_C_API_EXPORTS:
                raise RuntimeError(
                    "MinGW DLL export mismatch: missing=%s extra=%s" %
                    (sorted(WIREHAIR_C_API_EXPORTS - exports),
                     sorted(exports - WIREHAIR_C_API_EXPORTS)))
        else:
            find_unique_file(install_dir / "lib", "libwirehair.a")
            dlls = list((install_dir / "bin").glob("*wirehair*.dll"))
            if dlls:
                raise RuntimeError("static MinGW install unexpectedly contains a DLL")


def add_common_arguments(parser, default_build):
    parser.add_argument("--build-dir", type=Path, default=ROOT / default_build)
    parser.add_argument("--install-dir", type=Path, default=ROOT / (default_build + "-install"))
    parser.add_argument("--config", default="Release")
    parser.add_argument("--generator")
    parser.add_argument("--architecture")
    parser.add_argument("--jobs", type=int, default=max(1, min(os.cpu_count() or 1, 4)))
    parser.add_argument(
        "--cmake-arg",
        action="append",
        default=[],
        help="extra CMake configure argument; use --cmake-arg=-DNAME=VALUE",
    )


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    matrix = subparsers.add_parser("matrix", help="fast build/install/E2E lane")
    add_common_arguments(matrix, "build/ci-matrix")
    matrix.add_argument("--linkage", choices=("static", "shared"), required=True)
    matrix.add_argument(
        "--strict",
        action="store_true",
        help="treat portable compiler warnings as errors",
    )
    matrix.set_defaults(function=run_matrix)

    sanitizer = subparsers.add_parser("sanitizer", help="ASan+UBSan CTest lane")
    add_common_arguments(sanitizer, "build/ci-asan-ubsan")
    sanitizer.set_defaults(function=run_sanitizer)

    tsan = subparsers.add_parser("tsan", help="TSan cold-start concurrency lane")
    add_common_arguments(tsan, "build/ci-tsan")
    tsan.set_defaults(function=run_tsan)

    scheduled = subparsers.add_parser("scheduled", help="expensive deterministic E2E lane")
    add_common_arguments(scheduled, "build/ci-scheduled")
    scheduled.set_defaults(function=run_scheduled)

    qemu = subparsers.add_parser("qemu", help="emulated lower-feature x86 E2E lane")
    add_common_arguments(qemu, "build/ci-qemu")
    qemu.set_defaults(function=run_qemu)

    mingw = subparsers.add_parser("mingw", help="MinGW static/shared package ABI lane")
    add_common_arguments(mingw, "build/ci-mingw")
    mingw.set_defaults(function=run_mingw)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.jobs < 1:
        raise SystemExit("--jobs must be positive")
    args.build_dir = args.build_dir.resolve()
    args.install_dir = args.install_dir.resolve()
    args.function(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
