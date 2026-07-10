#!/usr/bin/env python3
"""Build and execute Wirehair's architecture-specific portability gates."""

import argparse
import os
from pathlib import Path
import platform
import re
import shlex
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]
EXPORT_MAP = ROOT / "abi" / "wirehair.map"
GCC_VERSION_RE = re.compile(r"^13\.3(?:\.|$)")
QEMU_VERSION_RE = re.compile(r"\bversion 8\.2(?:\.|\s)")

CROSS_TARGETS = {
    "aarch64": {
        "triplet": "aarch64-linux-gnu",
        "processor": "aarch64",
        "emulator": "qemu-aarch64",
        "sysroot": "/usr/aarch64-linux-gnu",
        "elf_class": "ELF64",
        "elf_data": "little endian",
        "elf_machine": "AArch64",
        "backend": "neon",
        "targets": ("gf256_inplace_test", "portability_roundtrip_test"),
        "test_regex": r"^(gf256_inplace_test|portability_roundtrip_test)$",
    },
    "i686": {
        "triplet": "i686-linux-gnu",
        "processor": "i686",
        "emulator": "qemu-i386",
        "sysroot": "/usr/i686-linux-gnu",
        "elf_class": "ELF32",
        "elf_data": "little endian",
        "elf_machine": "Intel 80386",
        "backend": "scalar",
        "targets": ("allocation_overflow_test", "portability_roundtrip_test"),
        "test_regex": r"^(allocation_overflow_test|portability_roundtrip_test)$",
    },
    "s390x": {
        "triplet": "s390x-linux-gnu",
        "processor": "s390x",
        "emulator": "qemu-s390x",
        "sysroot": "/usr/s390x-linux-gnu",
        "elf_class": "ELF64",
        "elf_data": "big endian",
        "elf_machine": "IBM S/390",
        "backend": "scalar",
        "targets": ("heavy_window_test", "portability_roundtrip_test"),
        "test_regex": r"^(heavy_window_test|portability_roundtrip_test)$",
    },
    "riscv64": {
        "triplet": "riscv64-linux-gnu",
        "processor": "riscv64",
        "emulator": "qemu-riscv64",
        "sysroot": "/usr/riscv64-linux-gnu",
        "elf_class": "ELF64",
        "elf_data": "little endian",
        "elf_machine": "RISC-V",
        "backend": "scalar",
        "targets": ("portability_roundtrip_test",),
        "test_regex": r"^portability_roundtrip_test$",
    },
}


def display_command(command):
    command = [str(item) for item in command]
    if os.name == "nt":
        return subprocess.list2cmdline(command)
    return shlex.join(command)


def run(command, *, env=None, log_path=None, input_text=None):
    command = [str(item) for item in command]
    print("+ " + display_command(command), flush=True)
    if log_path is None:
        subprocess.run(
            command,
            cwd=str(ROOT),
            env=env,
            input=input_text,
            text=input_text is not None,
            check=True,
        )
        return

    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8", newline="") as log:
        process = subprocess.Popen(
            command,
            cwd=str(ROOT),
            env=env,
            stdin=subprocess.PIPE if input_text is not None else None,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            errors="replace",
        )
        if input_text is not None:
            assert process.stdin is not None
            process.stdin.write(input_text)
            process.stdin.close()
        assert process.stdout is not None
        for line in process.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            log.write(line)
        return_code = process.wait()
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, command)


def capture(command, *, input_text=None, log_path=None, echo=True):
    command = [str(item) for item in command]
    print("+ " + display_command(command), flush=True)
    completed = subprocess.run(
        command,
        cwd=str(ROOT),
        input=input_text,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        errors="replace",
        check=False,
    )
    if echo:
        sys.stdout.write(completed.stdout)
        sys.stdout.flush()
    if log_path is not None:
        log_path = Path(log_path)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(completed.stdout, encoding="utf-8")
    if completed.returncode != 0:
        raise subprocess.CalledProcessError(
            completed.returncode, command, output=completed.stdout)
    return completed.stdout


def require_executable(name):
    executable = shutil.which(name)
    if not executable:
        raise RuntimeError("required portability tool was not found: " + name)
    return executable


def require_version(executable, arguments, expected, description):
    output = capture([executable, *arguments]).strip()
    if not expected.search(output):
        raise RuntimeError(
            "%s is outside the pinned version family: %s" %
            (description, output))
    return output


def compiler_macros(cxx):
    output = capture([
        cxx,
        "-std=c++11",
        "-I", ROOT,
        "-dM", "-E", "-x", "c++",
        "-include", ROOT / "gf256.h",
        "-",
    ], input_text="", echo=False)
    macros = {}
    for line in output.splitlines():
        match = re.match(r"#define\s+(\w+)(?:\s+(.*))?$", line)
        if match:
            macros[match.group(1)] = match.group(2) or ""
    return macros


def validate_gf256_backend(macros, expected_backend, *, linux_arm):
    if "GF256_TARGET_MOBILE" not in macros:
        raise RuntimeError("non-x86 portability target did not select portable GF256")
    has_neon = "GF256_TRY_NEON" in macros
    if (expected_backend == "neon") != has_neon:
        raise RuntimeError(
            "GF256 backend mismatch: expected %s, GF256_TRY_NEON=%s" %
            (expected_backend, has_neon))
    has_linux_arm = "LINUX_ARM" in macros
    if linux_arm != has_linux_arm:
        raise RuntimeError(
            "Linux ARM auxv guard mismatch: expected LINUX_ARM=%s, got %s" %
            (linux_arm, has_linux_arm))
    forbidden = (
        "GF256_TARGET_X86_SIMD", "GF256_TRY_SSSE3", "GF256_TRY_AVX2",
        "GF256_TRY_GFNI", "GF256_TRY_AVX512",
    )
    enabled_forbidden = [name for name in forbidden if name in macros]
    if enabled_forbidden:
        raise RuntimeError(
            "non-x86 portability target enabled x86 kernels: " +
            ", ".join(enabled_forbidden))
    print(
        "GF256 compile diagnostics: backend=%s LINUX_ARM=%s" %
        (expected_backend, has_linux_arm),
        flush=True,
    )


def common_cmake_arguments(build_dir, *, shared):
    return [
        "cmake", "-S", ROOT, "-B", build_dir, "-G", "Ninja",
        "-DBUILD_TESTS=ON",
        "-DBUILD_CODEC_V2=OFF",
        "-DBUILD_SHARED_LIBS=" + ("ON" if shared else "OFF"),
        "-DWIREHAIR_BUILD_BOTH=OFF",
        "-DWIREHAIR_BUILD_TOOLS=OFF",
        "-DWIREHAIR_BUILD_BENCHMARKS=OFF",
        "-DWIREHAIR_STRICT_WARNINGS=ON",
        "-DMARCH_NATIVE=OFF",
        "-DWH_LTO=OFF",
        "-DWH_PGO_MODE=OFF",
        "-DCMAKE_BUILD_TYPE=Release",
    ]


def build_targets(build_dir, jobs, targets):
    run([
        "cmake", "--build", build_dir,
        "--config", "Release",
        "--parallel", str(jobs),
        "--target", *targets,
    ], log_path=Path(build_dir) / "ci-logs" / "build.log")


def run_ctest(build_dir, test_regex, *, env=None):
    run([
        "ctest", "--test-dir", build_dir,
        "--build-config", "Release",
        "--output-on-failure",
        "--no-tests=error",
        "--verbose",
        "--tests-regex", test_regex,
    ], env=env)


def find_executable(build_dir, name):
    candidates = (
        Path(build_dir) / name,
        Path(build_dir) / "Release" / name,
    )
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    raise RuntimeError("cannot find %s under %s" % (name, build_dir))


def validate_elf_header(output, spec):
    required = (
        "Class:" + " " + spec["elf_class"],
        spec["elf_data"],
        "Machine:" + " " + spec["elf_machine"],
    )
    normalized = re.sub(r"[ \t]+", " ", output)
    missing = [fragment for fragment in required if fragment not in normalized]
    if missing:
        raise RuntimeError("ELF architecture mismatch; missing " + repr(missing))


def validate_cross_artifact(build_dir, spec):
    readelf = require_executable(spec["triplet"] + "-readelf")
    binary = find_executable(build_dir, "portability_roundtrip_test")
    output = capture(
        [readelf, "-h", binary],
        log_path=Path(build_dir) / "ci-logs" / "architecture.txt",
    )
    validate_elf_header(output, spec)


def cross_cmake_arguments(spec, cc, cxx, emulator):
    return [
        "-DCMAKE_SYSTEM_NAME=Linux",
        "-DCMAKE_SYSTEM_PROCESSOR=" + spec["processor"],
        "-DCMAKE_C_COMPILER=" + cc,
        "-DCMAKE_CXX_COMPILER=" + cxx,
        "-DCMAKE_CROSSCOMPILING_EMULATOR=%s;-L;%s" %
        (emulator, spec["sysroot"]),
    ]


def install_and_consume(
        build_dir, install_dir, jobs, cc, cxx, *, cross_arguments=(),
        runtime_key=None):
    run([
        "cmake", "--install", build_dir,
        "--config", "Release",
        "--prefix", install_dir,
    ])
    consumer_build = Path(build_dir) / "package-consumer"
    run([
        "cmake", "-S", ROOT / "test" / "package", "-B", consumer_build,
        "-G", "Ninja",
        "-DCMAKE_PREFIX_PATH=" + str(Path(install_dir).resolve()),
        "-DCMAKE_BUILD_TYPE=Release",
        "-DCMAKE_C_COMPILER=" + cc,
        "-DCMAKE_CXX_COMPILER=" + cxx,
        "-DWIREHAIR_STRICT_WARNINGS=ON",
        *cross_arguments,
    ])
    run([
        "cmake", "--build", consumer_build,
        "--config", "Release", "--parallel", str(jobs),
    ], log_path=Path(build_dir) / "ci-logs" / "consumer-build.log")
    env = os.environ.copy()
    if runtime_key:
        library_paths = [
            str(Path(install_dir).resolve() / "lib"),
            str(Path(install_dir).resolve() / "lib64"),
        ]
        old_value = env.get(runtime_key)
        if old_value:
            library_paths.append(old_value)
        env[runtime_key] = os.pathsep.join(library_paths)
    run_ctest(consumer_build, ".*", env=env)


def run_cross(args):
    spec = CROSS_TARGETS[args.target]
    cc = require_executable(spec["triplet"] + "-gcc-13")
    cxx = require_executable(spec["triplet"] + "-g++-13")
    emulator = require_executable(spec["emulator"])
    require_version(cc, ["-dumpfullversion"], GCC_VERSION_RE, "GCC cross compiler")
    require_version(emulator, ["--version"], QEMU_VERSION_RE, "QEMU emulator")
    if not Path(spec["sysroot"]).is_dir():
        raise RuntimeError("cross sysroot was not found: " + spec["sysroot"])

    macros = compiler_macros(cxx)
    validate_gf256_backend(
        macros, spec["backend"], linux_arm=args.target == "aarch64")
    cross_arguments = cross_cmake_arguments(spec, cc, cxx, emulator)
    run([
        *common_cmake_arguments(args.build_dir, shared=False),
        *cross_arguments,
    ])
    build_targets(args.build_dir, args.jobs, spec["targets"])
    validate_cross_artifact(args.build_dir, spec)
    run_ctest(args.build_dir, spec["test_regex"])

    # The native ARM lane is authoritative for packaging.  Keep the QEMU ARM
    # path capable of validating it too, so developers can reproduce the lane
    # without access to an ARM runner.
    if args.target == "aarch64":
        install_and_consume(
            args.build_dir, args.install_dir, args.jobs, cc, cxx,
            cross_arguments=cross_arguments)


def assert_host(expected_system, expected_machines):
    actual_system = platform.system()
    actual_machine = platform.machine().lower()
    print(
        "Host diagnostics: system=%s machine=%s" %
        (actual_system, actual_machine),
        flush=True,
    )
    if actual_system != expected_system or actual_machine not in expected_machines:
        raise RuntimeError(
            "wrong hosted runner: expected %s/%s, got %s/%s" %
            (expected_system, sorted(expected_machines),
             actual_system, actual_machine))


def run_native_arm64(args):
    assert_host("Linux", {"aarch64", "arm64"})
    cc = require_executable("gcc-13")
    cxx = require_executable("g++-13")
    require_version(cc, ["-dumpfullversion"], GCC_VERSION_RE, "native ARM GCC")
    macros = compiler_macros(cxx)
    validate_gf256_backend(macros, "neon", linux_arm=True)

    run([
        *common_cmake_arguments(args.build_dir, shared=True),
        "-DCMAKE_C_COMPILER=" + cc,
        "-DCMAKE_CXX_COMPILER=" + cxx,
    ])
    targets = ("wirehair", "gf256_inplace_test", "portability_roundtrip_test")
    build_targets(args.build_dir, args.jobs, targets)
    binary = find_executable(args.build_dir, "portability_roundtrip_test")
    readelf = require_executable("readelf")
    output = capture(
        [readelf, "-h", binary],
        log_path=Path(args.build_dir) / "ci-logs" / "architecture.txt")
    validate_elf_header(output, CROSS_TARGETS["aarch64"])
    run_ctest(
        args.build_dir,
        r"^(gf256_inplace_test|portability_roundtrip_test)$")
    install_and_consume(
        args.build_dir, args.install_dir, args.jobs, cc, cxx,
        runtime_key="LD_LIBRARY_PATH")


def load_public_exports(path=EXPORT_MAP):
    text = Path(path).read_text(encoding="ascii")
    match = re.fullmatch(
        r"\s*[A-Za-z_][A-Za-z0-9_.]*\s*\{\s*global:\s*"
        r"((?:wirehair_[A-Za-z0-9_]+\s*;\s*)+)"
        r"local:\s*\*\s*;\s*\}\s*;\s*",
        text,
    )
    if not match:
        raise RuntimeError("malformed public export map: " + str(path))
    symbols = re.findall(r"(wirehair_[A-Za-z0-9_]+)\s*;", match.group(1))
    if len(symbols) != len(set(symbols)):
        raise RuntimeError("duplicate public symbol in " + str(path))
    return frozenset(symbols)


def parse_macho_exports(output):
    exports = set()
    for line in output.splitlines():
        symbol = line.strip()
        if symbol.startswith("_wirehair_"):
            exports.add(symbol[1:])
        elif symbol:
            exports.add(symbol)
    return frozenset(exports)


def find_installed_dylib(install_dir):
    matches = sorted(Path(install_dir).rglob("libwirehair*.dylib"))
    if not matches:
        raise RuntimeError("installed Wirehair dylib was not found")
    # Ignore versioned symlink aliases, but reject genuinely distinct files.
    resolved = {path.resolve() for path in matches}
    if len(resolved) != 1:
        raise RuntimeError("multiple installed Wirehair dylibs found: %r" % matches)
    return next(iter(resolved))


def validate_macho_library(install_dir, build_dir):
    library = find_installed_dylib(install_dir)
    lipo = require_executable("lipo")
    architectures = capture([lipo, "-archs", library]).strip().split()
    if architectures != ["arm64"]:
        raise RuntimeError("Mach-O library is not arm64-only: %r" % architectures)

    nm = require_executable("nm")
    output = capture(
        [nm, "-gjU", library],
        log_path=Path(build_dir) / "ci-logs" / "macho-exports.txt")
    actual = parse_macho_exports(output)
    expected = load_public_exports()
    if actual != expected:
        raise RuntimeError(
            "Mach-O public export mismatch: missing=%s extra=%s" %
            (sorted(expected - actual), sorted(actual - expected)))
    print(
        "Mach-O ABI allowlist verified (%d symbols)" % len(actual),
        flush=True,
    )


def run_macos_arm64(args):
    assert_host("Darwin", {"arm64"})
    xcodebuild = require_executable("xcodebuild")
    version = capture([xcodebuild, "-version"])
    expected_xcode = "Xcode %s\nBuild version %s" % (
        args.xcode_version, args.xcode_build)
    if version.strip() != expected_xcode:
        raise RuntimeError(
            "Xcode pin mismatch: expected %r, got %r" %
            (expected_xcode, version.strip()))
    cc = require_executable("clang")
    cxx = require_executable("clang++")
    macros = compiler_macros(cxx)
    validate_gf256_backend(macros, "neon", linux_arm=False)

    run([
        *common_cmake_arguments(args.build_dir, shared=True),
        "-DCMAKE_C_COMPILER=" + cc,
        "-DCMAKE_CXX_COMPILER=" + cxx,
        "-DCMAKE_OSX_ARCHITECTURES=arm64",
    ])
    build_targets(
        args.build_dir, args.jobs,
        ("wirehair", "portability_roundtrip_test"))
    run_ctest(args.build_dir, r"^portability_roundtrip_test$")
    install_and_consume(
        args.build_dir, args.install_dir, args.jobs, cc, cxx,
        runtime_key="DYLD_LIBRARY_PATH")
    validate_macho_library(args.install_dir, args.build_dir)


def add_common_arguments(parser, default_name):
    parser.add_argument(
        "--build-dir", type=Path, default=ROOT / "build" / default_name)
    parser.add_argument(
        "--install-dir", type=Path,
        default=ROOT / "build" / (default_name + "-install"))
    parser.add_argument(
        "--jobs", type=int,
        default=max(1, min(os.cpu_count() or 1, 4)))


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    cross = subparsers.add_parser(
        "cross", help="cross-build and run a target under pinned QEMU")
    add_common_arguments(cross, "ci-portability-cross")
    cross.add_argument("--target", choices=sorted(CROSS_TARGETS), required=True)
    cross.set_defaults(function=run_cross)

    native_arm = subparsers.add_parser(
        "native-arm64", help="native Linux ARM64 SIMD and package lane")
    add_common_arguments(native_arm, "ci-portability-arm64")
    native_arm.set_defaults(function=run_native_arm64)

    macos_arm = subparsers.add_parser(
        "macos-arm64", help="native macOS ARM64 ABI and package lane")
    add_common_arguments(macos_arm, "ci-portability-macos-arm64")
    macos_arm.add_argument("--xcode-version", default="16.4")
    macos_arm.add_argument("--xcode-build", default="16F6")
    macos_arm.set_defaults(function=run_macos_arm64)
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
