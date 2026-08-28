#!/usr/bin/python3
"""Root supervisor for the one-shot WH2 retained-direct experiment.

This file deliberately has only three public modes.  ``--selftest`` is a
scratch-only, unprivileged check.  ``--seal-build`` creates the successful-only
root authority for one fresh detached build.  ``--execute`` is intended to be
run as the root-owned copy installed at :data:`INSTALLED_LAUNCHER` inside one
exact delegated systemd service.  The scientific controller, sampler, and
retained verifier always run without privilege.

The supervisor owns *process containment and durable attempt accounting*, not
scientific adjudication.  Only the retained verifier's authenticated exit
status can classify an attempt as pass, reject, invalid, or absent.
"""

from __future__ import annotations

import argparse
import array
import ctypes
import datetime
import errno
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import select
import shlex
import signal
import socket
import stat
import struct
import subprocess
import sys
import time
from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


ATTEMPT_ID = "wh2-retained-direct-systematic-complement-v2-r0"
CAMPAIGN = "wh2-retained-direct-systematic-complement-v2-r0"
UNIT_NAME = "wirehair-wh2-retained-direct-systematic-complement-v2-r0.service"
BUILD_UNIT_NAME = (
    "wirehair-wh2-retained-direct-systematic-complement-v2-r0-build.service"
)
INSTALLED_LAUNCHER = Path(
    "/usr/local/libexec/wirehair/Wh2DirectSystematicComplementLaunch.py"
)
SOURCE_LAUNCHER_RELATIVE = "bench/Wh2DirectSystematicComplementLaunch.py"
CONTROLLER_RELATIVE = "bench/Wh2DirectSystematicComplementScreen.py"
SAMPLER_RELATIVE = "bench/wirehair_expo_thermal_sampler.py"
BINARY_NAME = "wirehair_wh2_direct_systematic_complement_screen"
STATIC_SOURCE_PATHS = (
    "CMakeLists.txt",
    SOURCE_LAUNCHER_RELATIVE,
    "bench/Wh2DirectSystematicComplementScreen.cpp",
    CONTROLLER_RELATIVE,
    "bench/Wh2FrozenTrace.cpp",
    "bench/Wh2FrozenTrace.h",
    "bench/Wh2NativeCodec.cpp",
    "bench/Wh2NativeCodec.h",
    "bench/Wh2NativePanel.cpp",
    "bench/Wh2NativePanel.h",
    "bench/wh2_benchmark_contract.py",
    "bench/wh2_benchmark_contract_v4.json",
    "bench/wh2_native_short_screen.py",
    "bench/wh2_run_native_short_screen.py",
    SAMPLER_RELATIVE,
    "cmake/Wh2DirectSystematicComplementSymbolAudit.cmake",
    "cmake/Wh2TimingPolicySymbolAudit.cmake",
)
BUILD_PATHS = (
    "CMakeCache.txt", "compile_commands.json", "build.ninja",
    "CMakeFiles/rules.ninja",
)
COMPLEMENT_BUILD_SOURCES = (
    "bench/Wh2DirectSystematicComplementScreen.cpp",
    "bench/Wh2FrozenTrace.cpp",
    "bench/Wh2NativeCodec.cpp",
    "bench/Wh2NativePanel.cpp",
)
CORE_BUILD_SOURCES = (
    "wirehair.cpp", "gf256.cpp", "WirehairCodec.cpp", "WirehairTools.cpp",
    "codec/WirehairV2Codec.cpp", "codec/WirehairV2Peel.cpp",
    "codec/WirehairV2Plan.cpp", "codec/WirehairV2Policy.cpp",
    "codec/WirehairV2Precode.cpp", "codec/WirehairV2PrecodeDecode.cpp",
    "codec/WirehairV2PrecodeEncode.cpp", "codec/WirehairV2Profile.cpp",
    "codec/WirehairV2Seeds.cpp", "codec/WirehairV2Solve.cpp",
)
TIMING_POLICY_BUILD_SOURCES = (
    "codec/WirehairV2Codec.cpp", "codec/WirehairV2Peel.cpp",
    "codec/WirehairV2Plan.cpp", "codec/WirehairV2Policy.cpp",
    "codec/WirehairV2Precode.cpp", "codec/WirehairV2PrecodeDecode.cpp",
    "codec/WirehairV2PrecodeEncode.cpp", "codec/WirehairV2Seeds.cpp",
    "codec/WirehairV2Solve.cpp",
)
C_COMPILER = Path("/usr/bin/x86_64-linux-gnu-gcc-13")
CXX_COMPILER = Path("/usr/bin/x86_64-linux-gnu-g++-13")
FROZEN_TOOLCHAIN_IMAGES = {
    Path("/usr/bin/python3"): "1643dacd9feaedc58f3cc581e4d22577dfe25c09b10282936186ccf0f2e61118",
    Path("/usr/bin/setpriv"): "96b083b79c32fd2f0c29657e88e20c7495839349fc64ad5d0503f32d26bf8733",
    Path("/usr/bin/git"): "2a8c18fbf43da9f692d75474c72bea9dfd796c260b0f3dfe456376abc3bbd668",
    Path("/usr/bin/env"): "0aefff8f912fb75716c5d4de3b6acde93edbe8fa280fc8ee895c1226d3e373ef",
    Path("/bin/sh"): "86d31f6fb799e91fa21bad341484564510ca287703a16e9e46c53338776f4f42",
    Path("/usr/bin/bash"): "bc5945feb8bd26203ebfafea5ce1878bb2e32cb8fb50ab7ae395cfb1e1aaaef1",
    Path("/usr/bin/cat"): "a63158e6e5bce20616425f5d61e5bd7374bb5bccf15bbb93ae2e40238248f179",
    Path("/usr/bin/chmod"): "f7ba478ce7ff158ccd3ba01bac1c4231e68cbd490fe49b5562a6ea008a586696",
    Path("/usr/bin/dirname"): "5ed2ac7b191a3457d47a1c718397b819b04754174176c7cdd3c92faadab08055",
    Path("/usr/bin/grep"): "20a35a4c18f2fe9e6305c1cbf866b9c0fcc3f957132ae66028c99ec353bb0a80",
    Path("/usr/bin/mktemp"): "45d28b6224635b9506be3fd58b280e8267ea1418f251ccc84f9c4c1b56ce9e45",
    Path("/usr/bin/rm"): "aedf2ef491e1ca4ea2682b88a031f7c3960ebb133677a339743e859fa3d3d8b5",
    Path("/usr/bin/time"): "3b11dec50514a8473e9f6efa7a34d584d0657538c09988f61b72d38ad4991a10",
    C_COMPILER: "1b99826121ae6682a634e5efe09bd3e3df58ce58e0b28f849114ab5b89139c26",
    CXX_COMPILER: "1353e9bdd29a7295c7226bf6c63abccce056d8cac31f112e5cdbecc3f28c2769",
    Path("/usr/libexec/gcc/x86_64-linux-gnu/13/cc1"): "5d1679131184e2de4435b426eb264bf13472fe026db8e5c6bc97445814e8e2f4",
    Path("/usr/libexec/gcc/x86_64-linux-gnu/13/cc1plus"): "840b332fb62ec6f694ac77d91fe69ef7f80b0d69512ed89374af0ee7a506255d",
    Path("/usr/libexec/gcc/x86_64-linux-gnu/13/collect2"): "4d1f341ae5b763b513258ee2812422a45e063c30a2f1924a0cf63d3699f3a158",
    Path("/usr/bin/cmake"): "1c5227af4edd22d8d689def545e18ee458260c0fd579eba2187967f38817e638",
    Path("/usr/bin/ctest"): "45c60ab2448411a4f49a0bdb6cf5800fedbf0293048e3c82710a34190c7d1ef4",
    Path("/usr/bin/ninja"): "5965527e09fe2b3787772aa4f711d6a36b393e7f2fcaa744a7a96c5a4ddf59cb",
    Path("/usr/bin/as"): "21aff249b692b5c31a44007491f922dcb49f41323e362c57d2ada3f52eddb7f0",
    Path("/usr/bin/ld.bfd"): "e9ceb054c12207970f2726dfc07e9a66b411602748628baf27399f02a9bbb31b",
    Path("/usr/bin/ar"): "6452af2eea333b8c65e1adb92964fc8f97863ab003fa13f9d12bff5345cd7dbe",
    Path("/usr/bin/ranlib"): "60254978b8ee2c1b21b41d16a18a621dbbcc72e43eb2fa9f916818256c77e9ee",
    Path("/usr/bin/nm"): "3708df3d7ce16e2c1dcdcd722784ae382d4d582276a71ab7bbe2d542f5621202",
    Path("/usr/bin/objcopy"): "9276e918e3c26d3aa7ef9dd2232cd8ac692f03d3c1f270bfa5144a2b5492b713",
    Path("/usr/bin/objdump"): "325c4205a4c658a9d1e1ebc469ae55975a2b897a3d3c1e79d9b158612d37f745",
    Path("/usr/bin/strip"): "0d980587ada7ab12193f39271f060d5663aa2f289b0e80d2a0274ce7306e4e42",
    Path("/usr/bin/gcc-ar-13"): "b980c6ffedb8ceb443f2a95bd3c935a4d6aba09fa05da66864cbcaf2ec79f5c3",
    Path("/usr/bin/gcc-ranlib-13"): "c685333bc0ce5a0b432a2263184b378c4182a897e2239091adffac5c965baad1",
}

LAUNCH_DIR = Path(
    "/var/tmp/wh2-retained-direct-systematic-complement-v2-r0.launch"
)
SAMPLER_DIR = Path(
    "/var/tmp/wh2-retained-direct-systematic-complement-v2-r0.sampler"
)
FIXED_OUTPUT_DIR = Path(
    "/var/tmp/wh2-retained-direct-systematic-complement-v2-r0"
)
FIXED_CLAIM_PATH = Path(str(FIXED_OUTPUT_DIR) + ".claim")
LOCK_PATH = Path("/run/lock/wirehair-wh2-retained-direct-systematic-complement-v2-r0.lock")
BUILD_AUTHORITY_DIR = Path("/var/lib/wirehair")
BUILD_AUTHORITY_PATH = BUILD_AUTHORITY_DIR / (
    "wh2-retained-direct-systematic-complement-v2-r0-build-authority.json"
)

SAMPLER_CSV = SAMPLER_DIR / "thermal.csv"
SAMPLER_PID_FILE = SAMPLER_DIR / "sampler.pid"
SAMPLER_VALIDATION = SAMPLER_DIR / "validation.jsonl"
SAMPLER_RECEIPT = SAMPLER_DIR / "sampler-terminal.json"

PYTHON = Path("/usr/bin/python3")
SETPRIV = Path("/usr/bin/setpriv")
GIT = Path("/usr/bin/git")
SYSTEMD_RUN = Path("/usr/bin/systemd-run")
ENV_EXECUTABLE = Path("/usr/bin/env")
CTEST = Path("/usr/bin/ctest")
NINJA = Path("/usr/bin/ninja")
CGROUP_ROOT = Path("/sys/fs/cgroup")
HOST_MOUNT_ROOT = Path("/proc/1/root")

CAMPAIGN_UID = 1000
CAMPAIGN_GID = 1000
I2C_GID = 113
WORKER_CPU = 120
CONTROLLER_CPU = 121
SAMPLER_CPU = 122
VERIFIER_CPU = 123
RUN_CPUS = (WORKER_CPU, CONTROLLER_CPU, SAMPLER_CPU, VERIFIER_CPU)
PREFLIGHT_CPUS = (WORKER_CPU, CONTROLLER_CPU, SAMPLER_CPU)
MEMORY_NODES = "0"

SEALED_ENVIRONMENT = {
    "LANG": "C.UTF-8",
    "LC_ALL": "C.UTF-8",
    "PATH": "/usr/bin:/bin",
    "TZ": "UTC",
}

BUILD_ENVIRONMENT_FIXED = {
    "LANG": "C",
    "LC_ALL": "C",
    "PATH": "/usr/bin:/bin",
    "PYTHONDONTWRITEBYTECODE": "1",
    "TZ": "UTC",
}

PREFLIGHT_SCHEMA = "wirehair.wh2.direct-systematic-complement-preflight-seal.v1"
VERIFY_SCHEMA = "wirehair.wh2.direct-systematic-complement-retained-verification.v2"
CLAIM_SCHEMA = "wirehair.wh2.direct-systematic-complement-claim.v3"
SAMPLER_SCHEMA = "wirehair.wh2.thermal_sampler.v2"
ATTEMPT_SCHEMA = "wirehair.wh2.direct-systematic-complement-launch-attempt.v1"
AUTHORITY_SCHEMA = "wirehair.wh2.direct-systematic-complement-launch-authority.v1"
START_SCHEMA = "wirehair.wh2.direct-systematic-complement-launch-start.v1"
TERMINAL_SCHEMA = "wirehair.wh2.direct-systematic-complement-launch-terminal.v1"
ROLE_SCHEMA = "wirehair.wh2.direct-systematic-complement-launch-role.v1"
BUILD_AUTHORITY_SCHEMA = (
    "wirehair.wh2.direct-systematic-complement-build-authority.v1"
)
BUILD_CHILD_SECURITY_SCHEMA = (
    "wirehair.wh2.direct-systematic-complement-build-child-security.v1"
)
BUILD_TRUST_BOUNDARY = {
    "campaign_uid": (
        "trusted-no-concurrent-hostile-uid1000-process-during-build-or-launch"
    ),
    "root_runtime": (
        "trusted-root-owned-kernel-loader-shared-libraries-and-system-headers"
    ),
    "service_creation": "trusted-root-operator-invokes-exact-recorded-systemd-argv",
}
BUILD_AUTHORITY_KEYS = {
    "binary", "build_directory", "build_environment", "build_manifest_sha256",
    "build_profile", "build_subdirectories", "commands",
    "expected_source_commit", "git_sha256_after",
    "git_sha256_before", "installed_launcher", "receipt_sha256", "schema",
    "root_boundary",
    "source_manifest_sha256_after", "source_manifest_sha256_before",
    "source_root", "status", "systemd_run_argv", "toolchain_after",
    "toolchain_before",
    "trust_boundary",
}

MAX_PREFLIGHT_BYTES = 4 * 1024 * 1024
MAX_STDERR_BYTES = 1024 * 1024
MAX_SAMPLER_ARTIFACT_BYTES = 8 * 1024 * 1024
MAX_CLAIM_BYTES = 1024 * 1024
PREFLIGHT_SECONDS = 30.0
SAMPLER_READY_SECONDS = 12.0
CLAIM_ADMISSION_SECONDS = 12.0
OUTER_SECONDS = 120.0
SAMPLER_STOP_SECONDS = 6.0
VERIFIER_SECONDS = 10.0
WHOLE_RUN_SECONDS = 140.0
BUILD_CONFIGURE_SECONDS = 120.0
BUILD_COMPILE_SECONDS = 900.0
BUILD_TEST_SECONDS = 900.0
BUILD_NO_WORK_SECONDS = 60.0
MAX_BUILD_COMMAND_OUTPUT_BYTES = 4 * 1024 * 1024
MAX_BUILD_AUTHORITY_BYTES = 64 * 1024 * 1024
MAX_NINJA_TARGET_COMMAND_BYTES = 16 * 1024 * 1024
MAX_TRACKED_SOURCE_FILES = 4096
MAX_TRACKED_SOURCE_FILE_BYTES = 16 * 1024 * 1024
MAX_TRACKED_SOURCE_BYTES = 256 * 1024 * 1024

BUILD_CHILD_SECURITY_BOOTSTRAP = r'''
import hashlib,json,os,sys,time
fd=int(sys.argv[1]); tool=sys.argv[2:]
if fd < 3 or not tool or not os.path.isabs(tool[0]): raise SystemExit(121)
fields={}
with open('/proc/self/status','r',encoding='ascii') as source:
 for line in source:
  if ':' in line:
   key,value=line.split(':',1); fields[key]=value.strip()
def cap(name):
 value=fields.get(name,'')
 if len(value)!=16: raise SystemExit(122)
 return int(value,16)
value={'affinity':sorted(os.sched_getaffinity(0)),'cap_ambient':cap('CapAmb'),
'cap_bounding':cap('CapBnd'),'cap_effective':cap('CapEff'),
'cap_inheritable':cap('CapInh'),'cap_permitted':cap('CapPrm'),
'gid':list(os.getresgid()),'groups':sorted(os.getgroups()),
'no_new_privs':int(fields.get('NoNewPrivs','-1')),
'observed_monotonic_ns':time.monotonic_ns(),'pid':os.getpid(),
'receipt_sha256':None,'schema':
'wirehair.wh2.direct-systematic-complement-build-child-security.v1',
'uid':list(os.getresuid())}
encode=lambda item:json.dumps(item,allow_nan=False,ensure_ascii=True,
 separators=(',',':'),sort_keys=True).encode('ascii')
value['receipt_sha256']=hashlib.sha256(encode(value)).hexdigest()
raw=encode(value)+b'\n'; offset=0
while offset < len(raw):
 written=os.write(fd,raw[offset:])
 if written <= 0: raise SystemExit(123)
 offset+=written
os.close(fd); os.execve(tool[0],tool,dict(os.environ))
'''.strip()

LOWER40 = re.compile(r"[0-9a-f]{40}\Z")
LOWER64 = re.compile(r"[0-9a-f]{64}\Z")
CANONICAL_UINT = re.compile(r"(?:0|[1-9][0-9]*)\Z")
UTC_MICROSECOND = re.compile(
    r"[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}"
    r"\.[0-9]{6}Z\Z"
)
CGROUP2_SUPER_MAGIC = 0x63677270
PR_SET_CHILD_SUBREAPER = 36
PR_GET_CHILD_SUBREAPER = 37

PREFLIGHT_OPTION_ORDER = (
    "--binary", "--build-dir", "--sampler-pid", "--sampler-script",
    "--sampler-csv", "--sampler-pid-file", "--sampler-validation-jsonl",
    "--sampler-receipt", "--expected-source-commit",
)
RUN_OPTION_ORDER = (
    "--binary", "--build-dir", "--cpu", "--controller-cpu",
    "--sampler-pid", "--sampler-cpu", "--sampler-script", "--sampler-csv",
    "--sampler-pid-file", "--sampler-validation-jsonl", "--sampler-receipt",
    "--expected-source-commit", "--expected-binary-sha256",
    "--expected-binary-uid", "--expected-build-manifest-sha256",
    "--expected-controller-sha256", "--expected-controller-uid",
    "--expected-controller-gid", "--expected-git-sha256",
    "--expected-python-sha256", "--expected-sampler-process-start-ticks",
    "--expected-sampler-script-sha256", "--expected-sampler-csv-device",
    "--expected-sampler-csv-inode", "--expected-sampler-pid-file-device",
    "--expected-sampler-pid-file-inode",
    "--expected-sampler-validation-device",
    "--expected-sampler-validation-inode",
    "--expected-sampler-receipt-device",
    "--expected-sampler-receipt-inode", "--expected-sampler-cmdline-sha256",
    "--expected-sampler-environ-sha256",
    "--expected-sampler-executable-sha256", "--expected-sampler-uid",
    "--expected-sampler-gid", "--expected-sampler-i2c-gid",
    "--expected-source-manifest-sha256",
)

PREFLIGHT_KEYS = {
    "binary_after", "binary_before", "build_manifest_after",
    "build_manifest_before", "expected_source_commit", "git_after",
    "git_before", "health_module_loader", "preflight_controller_after",
    "preflight_controller_before", "receipt_sha256", "run_argv",
    "run_argv_sha256", "sampler_after", "sampler_before",
    "sampler_prefix_binding", "schema", "source_manifest_after",
    "source_manifest_before", "source_root",
}

VERIFY_EXIT = {"pass": 0, "reject": 2, "invalid": 3, "absent": 4}
EXIT_VERIFY = {value: key for key, value in VERIFY_EXIT.items()}


class LaunchError(RuntimeError):
    """A fail-closed launcher error."""


class AttemptConsumedError(LaunchError):
    """Irrevocable namespace consumption, possibly without a usable journal."""

    def __init__(self, message: str,
                 journal: Optional["AttemptJournal"] = None) -> None:
        super().__init__(message)
        self.journal = journal


def fail(message: str) -> None:
    raise LaunchError(message)


def exception_text(exc: BaseException) -> str:
    try:
        value = f"{type(exc).__name__}: {exc}"
    except BaseException:
        value = "BaseException"
    return value[:4096]


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def canonical_bytes(value: Any) -> bytes:
    try:
        return json.dumps(
            value, allow_nan=False, ensure_ascii=True,
            separators=(",", ":"), sort_keys=True,
        ).encode("ascii")
    except (TypeError, ValueError, UnicodeEncodeError) as exc:
        fail("canonical JSON encoding failed: " + exception_text(exc))


def canonical_equal(left: Any, right: Any) -> bool:
    return canonical_bytes(left) == canonical_bytes(right)


def unique_object(pairs: Sequence[Tuple[str, Any]]) -> Dict[str, Any]:
    value: Dict[str, Any] = {}
    for key, item in pairs:
        if key in value:
            fail("duplicate JSON member: " + key)
        value[key] = item
    return value


def reject_constant(value: str) -> None:
    fail("non-finite JSON constant: " + value)


def parse_canonical_line(raw: bytes, limit: int, where: str) -> Dict[str, Any]:
    if (
        not raw or len(raw) > limit or not raw.endswith(b"\n")
        or raw.count(b"\n") != 1 or b"\r" in raw
    ):
        fail(where + " framing differs")
    try:
        value = json.loads(
            raw[:-1].decode("ascii"), object_pairs_hook=unique_object,
            parse_constant=reject_constant,
        )
    except LaunchError:
        raise
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        fail(where + " is malformed: " + exception_text(exc))
    if type(value) is not dict or canonical_bytes(value) + b"\n" != raw:
        fail(where + " is not one canonical ASCII JSON line")
    return value


def exact_keys(value: Mapping[str, Any], expected: Iterable[str], where: str) -> None:
    if type(value) is not dict or set(value) != set(expected):
        fail(where + " member roster differs")


def exact_uint(value: Any, lower: int, upper: int, where: str) -> int:
    if type(value) is not int or value < lower or value > upper:
        fail(where + " integer differs")
    return value


def exact_hash(value: Any, where: str) -> str:
    if type(value) is not str or LOWER64.fullmatch(value) is None:
        fail(where + " SHA-256 differs")
    return value


def exact_utc(value: Any, where: str) -> datetime.datetime:
    if type(value) is not str or UTC_MICROSECOND.fullmatch(value) is None:
        fail(where + " UTC timestamp differs")
    try:
        parsed = datetime.datetime.strptime(value, "%Y-%m-%dT%H:%M:%S.%fZ")
    except ValueError as exc:
        fail(where + " UTC timestamp is invalid: " + exception_text(exc))
    return parsed


def hash_path(path: Path, limit: int, where: str,
              require_single_link: bool = True) -> Tuple[str, os.stat_result]:
    flags = os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK | getattr(os, "O_NOFOLLOW", 0)
    fd = os.open(str(path), flags)
    try:
        before = os.fstat(fd)
        if (
            not stat.S_ISREG(before.st_mode) or before.st_nlink < 1
            or (require_single_link and before.st_nlink != 1)
            or before.st_size < 1 or before.st_size > limit
        ):
            fail(where + " regular-file policy differs")
        digest = hashlib.sha256()
        offset = 0
        while offset < before.st_size:
            block = os.pread(fd, min(1024 * 1024, before.st_size - offset), offset)
            if not block:
                fail(where + " read was short")
            digest.update(block)
            offset += len(block)
        after = os.fstat(fd)
        named = os.stat(str(path), follow_symlinks=False)
        if (
            (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns,
             before.st_ctime_ns, before.st_mode, before.st_uid, before.st_gid,
             before.st_nlink)
            != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns,
                after.st_ctime_ns, after.st_mode, after.st_uid, after.st_gid,
                after.st_nlink)
            or (after.st_dev, after.st_ino) != (named.st_dev, named.st_ino)
        ):
            fail(where + " changed while hashing")
        return digest.hexdigest(), after
    finally:
        os.close(fd)


def sealed_file_receipt(path: Path, relative: str, mode: int, limit: int,
                        where: str) -> Dict[str, Any]:
    digest, info = hash_path(path, limit, where)
    if (
        info.st_uid != CAMPAIGN_UID or info.st_gid != CAMPAIGN_GID
        or stat.S_IMODE(info.st_mode) != mode
    ):
        fail(where + " ownership/mode seal differs")
    return {
        "bytes": info.st_size, "device": info.st_dev, "gid": info.st_gid,
        "inode": info.st_ino, "mode": stat.S_IMODE(info.st_mode),
        "nlink": info.st_nlink, "path": relative, "sha256": digest,
        "uid": info.st_uid,
    }


def static_source_manifest(root: Path) -> Tuple[str, Dict[str, str]]:
    hashes: Dict[str, str] = {}
    preimage = bytearray()
    for relative in STATIC_SOURCE_PATHS:
        receipt = sealed_file_receipt(
            root / relative, relative, 0o444, 1024 * 1024,
            "static source " + relative,
        )
        digest = receipt["sha256"]
        hashes[relative] = digest
        preimage.extend((digest + "  " + relative + "\n").encode("ascii"))
    return sha256_bytes(bytes(preimage)), hashes


def parse_cmake_cache_entries(raw: bytes) -> Dict[str, Tuple[str, str]]:
    if not raw or b"\0" in raw or b"\r" in raw:
        fail("build-profile CMake cache framing differs")
    try:
        lines = raw.decode("utf-8").split("\n")
    except UnicodeDecodeError as exc:
        fail("build-profile CMake cache is not UTF-8: " + exception_text(exc))
    entries: Dict[str, Tuple[str, str]] = {}
    for line in lines:
        if not line or line.startswith("//") or line.startswith("#"):
            continue
        if ":" not in line or "=" not in line:
            fail("build-profile CMake cache line differs")
        name, remainder = line.split(":", 1)
        kind, value = remainder.split("=", 1)
        if (
            re.fullmatch(r"[A-Za-z_][A-Za-z0-9_-]*", name) is None
            or not kind or name in entries
        ):
            fail("build-profile CMake cache member differs: " + name)
        entries[name] = (kind, value)
    return entries


def compiler_image_receipt(path: Path, language: str) -> Dict[str, Any]:
    if path not in (C_COMPILER, CXX_COMPILER):
        fail("build-profile compiler path differs")
    resolved = path.resolve(strict=True)
    expected_suffix = "gcc-13" if language == "c" else "g++-13"
    if language not in ("c", "cxx") or not resolved.name.endswith(expected_suffix):
        fail("build-profile compiler family differs")
    digest, info = hash_path(
        resolved, 64 * 1024 * 1024, "build-profile " + language + " compiler",
    )
    if (
        digest != FROZEN_TOOLCHAIN_IMAGES[path]
        or
        info.st_uid != 0 or info.st_gid != 0
        or stat.S_IMODE(info.st_mode) != 0o755
    ):
        fail("build-profile compiler image policy differs")
    environment = {
        "LANG": "C", "LC_ALL": "C", "PATH": "/usr/bin:/bin", "TZ": "UTC",
    }

    def query(argument: str) -> str:
        result = subprocess.run(
            [str(path), argument], env=environment, stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            timeout=5.0,
        )
        try:
            value = result.stdout.decode("ascii").rstrip("\n")
        except UnicodeDecodeError as exc:
            fail("build-profile compiler response differs: " + exception_text(exc))
        if (
            result.returncode != 0 or result.stderr or not value
            or "\n" in value or "\r" in value
        ):
            fail("build-profile compiler query failed")
        return value

    version = query("-dumpfullversion")
    machine = query("-dumpmachine")
    if re.fullmatch(r"13(?:\.[0-9]+){1,2}", version) is None:
        fail("build-profile compiler version differs")
    if machine != "x86_64-linux-gnu":
        fail("build-profile compiler machine differs")
    return {
        "configured_path": str(path), "device": info.st_dev,
        "gid": info.st_gid, "inode": info.st_ino, "machine": machine,
        "mode": stat.S_IMODE(info.st_mode), "resolved_path": str(resolved),
        "sha256": digest, "uid": info.st_uid, "version": version,
    }


def frozen_toolchain_receipts() -> List[Dict[str, Any]]:
    receipts: List[Dict[str, Any]] = []
    for configured, expected_hash in FROZEN_TOOLCHAIN_IMAGES.items():
        resolved = configured.resolve(strict=True)
        digest, info = hash_path(
            resolved, 128 * 1024 * 1024,
            "frozen toolchain image " + str(configured),
            require_single_link=(
                configured not in (Path("/usr/bin/python3"), Path("/usr/bin/bash"))
            ),
        )
        if (
            digest != expected_hash or info.st_uid != 0 or info.st_gid != 0
            or stat.S_IMODE(info.st_mode) != 0o755
        ):
            fail("frozen toolchain image differs: " + str(configured))
        receipts.append({
            "bytes": info.st_size, "configured_path": str(configured),
            "device": info.st_dev,
            "gid": info.st_gid, "inode": info.st_ino,
            "mode": stat.S_IMODE(info.st_mode), "nlink": info.st_nlink,
            "resolved_path": str(resolved), "sha256": digest,
            "uid": info.st_uid,
        })
    return receipts


def ninja_target_commands(build: Path) -> List[str]:
    result = subprocess.run(
        [str(NINJA), "-f", "build.ninja", "-t", "commands", BINARY_NAME],
        cwd=str(build), env={
            "LANG": "C", "LC_ALL": "C", "PATH": "/usr/bin:/bin", "TZ": "UTC",
        }, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, check=False, timeout=10.0,
    )
    raw = result.stdout
    if (
        result.returncode != 0 or result.stderr or not raw
        or len(raw) > MAX_NINJA_TARGET_COMMAND_BYTES
        or not raw.endswith(b"\n") or b"\0" in raw or b"\r" in raw
    ):
        fail("build-profile effective Ninja command query differs")
    try:
        lines = raw.decode("ascii").splitlines()
    except UnicodeDecodeError as exc:
        fail("build-profile effective Ninja commands are not ASCII: "
             + exception_text(exc))
    if not lines or len(lines) > 4096 or any(not line for line in lines):
        fail("build-profile effective Ninja command roster differs")
    return lines


def validate_effective_ninja_commands(
    lines: Sequence[str], compile_tokens: Mapping[str, Sequence[str]],
    direct_outputs: Sequence[str],
) -> None:
    if type(lines) not in (list, tuple) or not lines:
        fail("build-profile effective Ninja command roster differs")
    seen_compile: Dict[str, List[str]] = {}
    archives: set = set()
    direct_link_count = 0
    expected_link = (
        ": && " + str(CXX_COMPILER) + " -O3 -DNDEBUG  "
        + " ".join(direct_outputs) + " -o " + BINARY_NAME
        + "  codec/libwirehair_v2_timing_policy.a  libwirehair.a  -lm && :"
    )
    archive_objects = {
        "libwirehair.a": [
            "CMakeFiles/wirehair.dir/" + relative + ".o"
            for relative in CORE_BUILD_SOURCES
        ],
        "codec/libwirehair_v2_timing_policy.a": [
            "codec/CMakeFiles/wirehair_v2_timing_policy.dir/"
            + Path(relative).name + ".o"
            for relative in TIMING_POLICY_BUILD_SOURCES
        ],
    }
    expected_archives = {
        archive: (
            ": && /usr/bin/cmake -E rm -f " + archive
            + " && /usr/bin/x86_64-linux-gnu-ar qc " + archive + "  "
            + " ".join(objects)
            + " && /usr/bin/x86_64-linux-gnu-ranlib " + archive + " && :"
        )
        for archive, objects in archive_objects.items()
    }
    expected_compile_outputs = set(direct_outputs)
    for values in archive_objects.values():
        expected_compile_outputs.update(values)
    for line in lines:
        if type(line) is not str or any(
            character in line for character in "\r\n\0"
        ):
            fail("build-profile effective Ninja command framing differs")
        if line == expected_link:
            direct_link_count += 1
            continue
        matched_archive = next(
            (
                archive for archive, expected in expected_archives.items()
                if line == expected
            ),
            None,
        )
        if matched_archive is not None:
            if matched_archive in archives:
                fail("build-profile effective archive command differs")
            archives.add(matched_archive)
            continue
        try:
            tokens = shlex.split(line, posix=True)
        except ValueError as exc:
            fail("build-profile effective Ninja command shell syntax differs: "
                 + exception_text(exc))
        if tokens.count("-c") == 1 and tokens.count("-o") == 1:
            output_index = tokens.index("-o")
            if output_index + 1 >= len(tokens):
                fail("build-profile effective compile output differs")
            output = tokens[output_index + 1]
            if output not in compile_tokens or output in seen_compile:
                fail("build-profile effective compile roster differs")
            expected = list(compile_tokens[output])
            expected_output_index = expected.index("-o")
            expected[expected_output_index:expected_output_index] = [
                "-MD", "-MT", output, "-MF", output + ".d",
            ]
            if (
                tokens != expected
                or any("WIREHAIR_V2_ENABLE_TEST_HOOKS" in token for token in tokens)
            ):
                fail("build-profile effective compile command differs: " + output)
            seen_compile[output] = tokens
            continue
        fail("build-profile effective Ninja command is unexpected")
    if (
        direct_link_count != 1
        or set(seen_compile) != expected_compile_outputs
        or archives != set(expected_archives)
    ):
        fail("build-profile effective direct target command roster differs")


def validate_build_profile(
    build: Path, source_root: Path, expected_commit: str,
) -> Tuple[str, Dict[str, Any]]:
    """Validate and hash one held-FD snapshot of the frozen build profile."""
    opened: Dict[str, int] = {}
    before: Dict[str, os.stat_result] = {}
    payloads: Dict[str, bytes] = {}
    receipts: Dict[str, Dict[str, Any]] = {}
    try:
        for relative in BUILD_PATHS:
            path = build / relative
            fd = os.open(
                str(path), os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK
                | getattr(os, "O_NOFOLLOW", 0),
            )
            opened[relative] = fd
            info = os.fstat(fd)
            named = os.stat(str(path), follow_symlinks=False)
            if (
                not stat.S_ISREG(info.st_mode) or info.st_nlink != 1
                or info.st_size < 1 or info.st_size > 64 * 1024 * 1024
                or info.st_uid != CAMPAIGN_UID or info.st_gid != CAMPAIGN_GID
                or stat.S_IMODE(info.st_mode) != 0o444
                or (info.st_dev, info.st_ino) != (named.st_dev, named.st_ino)
            ):
                fail("build-profile artifact policy differs: " + relative)
            before[relative] = info
            data = bytearray()
            offset = 0
            while offset < info.st_size:
                block = os.pread(
                    fd, min(1024 * 1024, info.st_size - offset), offset,
                )
                if not block:
                    fail("build-profile artifact read was short: " + relative)
                data.extend(block)
                offset += len(block)
            payloads[relative] = bytes(data)
            receipts[relative] = {
                "bytes": info.st_size, "device": info.st_dev,
                "gid": info.st_gid, "inode": info.st_ino,
                "mode": stat.S_IMODE(info.st_mode), "nlink": info.st_nlink,
                "path": relative, "sha256": sha256_bytes(payloads[relative]),
                "uid": info.st_uid,
            }

        cache = parse_cmake_cache_entries(payloads["CMakeCache.txt"])
        required_cache = {
            "BUILD_CODEC_V2": ("BOOL", "ON"),
            "BUILD_SHARED_LIBS": ("BOOL", "OFF"),
            "BUILD_TESTS": ("BOOL", "ON"),
            "CMAKE_AR": ("FILEPATH", "/usr/bin/x86_64-linux-gnu-ar"),
            "CMAKE_BUILD_TYPE": ("STRING", "Release"),
            "CMAKE_CACHEFILE_DIR": ("INTERNAL", str(build)),
            "CMAKE_CXX_COMPILER": ("STRING", str(CXX_COMPILER)),
            "CMAKE_CXX_COMPILER_AR": (
                "FILEPATH", "/usr/bin/x86_64-linux-gnu-gcc-ar-13"
            ),
            "CMAKE_CXX_COMPILER_RANLIB": (
                "FILEPATH", "/usr/bin/x86_64-linux-gnu-gcc-ranlib-13"
            ),
            "CMAKE_CXX_FLAGS": ("STRING", ""),
            "CMAKE_CXX_FLAGS_RELEASE": ("STRING", "-O3 -DNDEBUG"),
            "CMAKE_C_COMPILER": ("STRING", str(C_COMPILER)),
            "CMAKE_C_COMPILER_AR": (
                "FILEPATH", "/usr/bin/x86_64-linux-gnu-gcc-ar-13"
            ),
            "CMAKE_C_COMPILER_RANLIB": (
                "FILEPATH", "/usr/bin/x86_64-linux-gnu-gcc-ranlib-13"
            ),
            "CMAKE_C_FLAGS": ("STRING", ""),
            "CMAKE_C_FLAGS_RELEASE": ("STRING", "-O3 -DNDEBUG"),
            "CMAKE_EXE_LINKER_FLAGS": ("STRING", ""),
            "CMAKE_EXE_LINKER_FLAGS_RELEASE": ("STRING", ""),
            "CMAKE_EXPORT_COMPILE_COMMANDS": ("UNINITIALIZED", "ON"),
            "CMAKE_GENERATOR": ("INTERNAL", "Ninja"),
            "CMAKE_GENERATOR_INSTANCE": ("INTERNAL", ""),
            "CMAKE_GENERATOR_PLATFORM": ("INTERNAL", ""),
            "CMAKE_GENERATOR_TOOLSET": ("INTERNAL", ""),
            "CMAKE_HOME_DIRECTORY": ("INTERNAL", str(source_root)),
            "CMAKE_LINKER": ("FILEPATH", "/usr/bin/x86_64-linux-gnu-ld"),
            "CMAKE_MAKE_PROGRAM": ("FILEPATH", "/usr/bin/ninja"),
            "CMAKE_MODULE_LINKER_FLAGS": ("STRING", ""),
            "CMAKE_MODULE_LINKER_FLAGS_RELEASE": ("STRING", ""),
            "CMAKE_NM": ("FILEPATH", "/usr/bin/x86_64-linux-gnu-nm"),
            "CMAKE_OBJCOPY": (
                "FILEPATH", "/usr/bin/x86_64-linux-gnu-objcopy"
            ),
            "CMAKE_OBJDUMP": (
                "FILEPATH", "/usr/bin/x86_64-linux-gnu-objdump"
            ),
            "CMAKE_SHARED_LINKER_FLAGS": ("STRING", ""),
            "CMAKE_SHARED_LINKER_FLAGS_RELEASE": ("STRING", ""),
            "CMAKE_RANLIB": (
                "FILEPATH", "/usr/bin/x86_64-linux-gnu-ranlib"
            ),
            "CMAKE_STATIC_LINKER_FLAGS": ("STRING", ""),
            "CMAKE_STATIC_LINKER_FLAGS_RELEASE": ("STRING", ""),
            "CMAKE_STRIP": (
                "FILEPATH", "/usr/bin/x86_64-linux-gnu-strip"
            ),
            "CMAKE_PROJECT_NAME": ("STATIC", "wirehair"),
            "MARCH_NATIVE": ("BOOL", "OFF"),
            "Python3_EXECUTABLE": ("UNINITIALIZED", str(PYTHON)),
            "WH_PGO_DIR": ("PATH", str(build / "pgo")),
            "WH_LTO": ("STRING", "OFF"),
            "WH_PGO_MODE": ("STRING", "OFF"),
            "WIREHAIR_BUILD_BENCHMARKS": ("BOOL", "OFF"),
            "WIREHAIR_BUILD_BOTH": ("BOOL", "OFF"),
            "WIREHAIR_BUILD_TOOLS": ("BOOL", "OFF"),
            "WIREHAIR_ENABLE_LIBFUZZER": ("BOOL", "OFF"),
            "WIREHAIR_ENABLE_SCHEDULED_TESTS": ("BOOL", "OFF"),
            "WIREHAIR_ENABLE_TRI_BMQ4_SCHUR_FALSIFIER": ("BOOL", "OFF"),
            "WIREHAIR_STATIC_PIC": ("BOOL", "ON"),
            "WIREHAIR_STRICT_WARNINGS": ("BOOL", "ON"),
            "_Python3_EXECUTABLE": ("INTERNAL", str(PYTHON)),
            "wirehair_BINARY_DIR": ("STATIC", str(build)),
            "wirehair_SOURCE_DIR": ("STATIC", str(source_root)),
        }
        for name, expected in required_cache.items():
            if cache.get(name) != expected:
                fail("build-profile cache value differs: " + name)
        forbidden_cache = {
            "CMAKE_C_COMPILER_ARG1", "CMAKE_C_COMPILER_LAUNCHER",
            "CMAKE_C_LINKER_LAUNCHER", "CMAKE_C_COMPILER_TARGET",
            "CMAKE_C_COMPILER_EXTERNAL_TOOLCHAIN",
            "CMAKE_CXX_COMPILER_ARG1", "CMAKE_CXX_COMPILER_LAUNCHER",
            "CMAKE_CXX_LINKER_LAUNCHER", "CMAKE_CXX_COMPILER_TARGET",
            "CMAKE_CXX_COMPILER_EXTERNAL_TOOLCHAIN",
            "CMAKE_CXX_STANDARD", "CMAKE_CXX_STANDARD_REQUIRED",
            "CMAKE_CXX_EXTENSIONS", "CMAKE_POSITION_INDEPENDENT_CODE",
            "CMAKE_PROJECT_INCLUDE", "CMAKE_PROJECT_INCLUDE_BEFORE",
            "CMAKE_PROJECT_TOP_LEVEL_INCLUDES", "CMAKE_TOOLCHAIN_FILE",
            "CMAKE_SYSROOT", "CMAKE_SYSROOT_COMPILE", "CMAKE_SYSROOT_LINK",
            "RULE_LAUNCH_COMPILE", "RULE_LAUNCH_CUSTOM", "RULE_LAUNCH_LINK",
        }
        forbidden_prefixes = (
            "CMAKE_INTERPROCEDURAL_OPTIMIZATION", "CMAKE_UNITY_BUILD",
            "CMAKE_USER_MAKE_RULES_OVERRIDE",
        )
        forbidden_present = forbidden_cache.intersection(cache)
        forbidden_present.update(
            name for name in cache
            if any(name.startswith(prefix) for prefix in forbidden_prefixes)
            or (
                name.startswith("CMAKE_PROJECT_")
                and (name.endswith("_INCLUDE") or "_INCLUDE_" in name)
            )
        )
        if forbidden_present:
            fail("build-profile cache contains a launcher/toolchain injection")

        try:
            commands = json.loads(
                payloads["compile_commands.json"].decode("utf-8"),
                object_pairs_hook=unique_object, parse_constant=reject_constant,
            )
        except LaunchError:
            raise
        except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
            fail("build-profile compile commands malformed: " + exception_text(exc))
        if (
            type(commands) is not list or not commands or len(commands) > 4096
            or any(type(item) is not dict for item in commands)
        ):
            fail("build-profile compile command roster differs")
        compile_outputs: Dict[str, str] = {}
        compile_tokens: Dict[str, List[str]] = {}
        forbidden_compile_prefixes = (
            "-march", "-mtune", "-mcpu", "-flto", "-fprofile",
            "-fsanitize", "-fno-sanitize", "-fplugin", "-specs",
            "--sysroot", "-isysroot", "-finstrument", "--coverage",
            "-fcoverage", "-fprofile-arcs", "-ftest-coverage", "-wrapper",
        )
        for item in commands:
            exact_keys(
                item, {"command", "directory", "file", "output"},
                "build-profile compile command",
            )
            if any(type(item[name]) is not str for name in item):
                fail("build-profile compile command member differs")
            try:
                tokens = shlex.split(item["command"], posix=True)
            except ValueError as exc:
                fail("build-profile compile command shell syntax differs: "
                     + exception_text(exc))
            file_path = Path(item["file"])
            output = item["output"]
            if (
                not tokens or item["directory"] != str(build)
                or not file_path.is_absolute()
                or file_path.resolve(strict=True) != file_path
                or source_root not in file_path.parents
                or Path(output).is_absolute()
                or os.path.normpath(output) != output
                or output.startswith("../") or output in compile_outputs
                or tokens.count("-c") != 1 or tokens.count("-o") != 1
                or tokens[-1] != str(file_path)
                or tokens.count("-O3") != 1 or tokens.count("-DNDEBUG") != 1
            ):
                fail("build-profile compile command shape differs")
            output_index = tokens.index("-o")
            compile_index = tokens.index("-c")
            if (
                output_index + 1 >= len(tokens)
                or tokens[output_index + 1] != output
                or compile_index + 1 >= len(tokens)
                or tokens[compile_index + 1] != str(file_path)
            ):
                fail("build-profile compile command input/output differs")
            expected_compiler = (
                str(C_COMPILER) if file_path.suffix == ".c" else str(CXX_COMPILER)
            )
            if tokens[0] != expected_compiler:
                fail("build-profile compile command compiler differs")
            if (
                expected_compiler == str(CXX_COMPILER)
                and tokens.count("-std=gnu++11") != 1
            ):
                fail("build-profile C++ language mode differs")
            for token in tokens:
                lowered = token.lower()
                if (
                    (token.startswith("-O") and token != "-O3")
                    or token.startswith("@")
                    or token in ("-pg", "-p")
                    or any(character in token for character in ";|&<>`\r\n\0")
                    or "$(" in token or "${" in token
                    or any(
                        lowered.startswith(prefix)
                        for prefix in forbidden_compile_prefixes
                    )
                    or "WIREHAIR_V2_TEST_HOOK" in token
                ):
                    fail("build-profile compile command contains a forbidden flag")
            compile_outputs[output] = str(file_path)
            compile_tokens[output] = tokens
        object_prefix = "CMakeFiles/" + BINARY_NAME + ".dir/bench/"
        target_entries = [
            item for item in commands
            if type(item.get("output")) is str
            and item["output"].startswith(object_prefix)
        ]
        expected_entries: List[Dict[str, str]] = []
        for relative in COMPLEMENT_BUILD_SOURCES:
            source = source_root / relative
            output = object_prefix + Path(relative).name + ".o"
            command = " ".join((
                str(CXX_COMPILER), "-DWIREHAIR_STATIC=1",
                "-DWIREHAIR_V2_ENABLE_BENCHMARK_EQUATIONS=1",
                "-DWIREHAIR_WH2_SOURCE_GIT_COMMIT=\\\""
                + expected_commit + "\\\"",
                "-I" + str(source_root / "codec"),
                "-I" + str(source_root / "include"),
                "-O3", "-DNDEBUG", "-std=gnu++11", "-Wall", "-Wextra",
                "-Wpedantic", "-Werror", "-o", output, "-c", str(source),
            ))
            expected_entries.append({
                "command": command, "directory": str(build),
                "file": str(source), "output": output,
            })
        if target_entries != expected_entries:
            fail("build-profile target compile commands differ")
        expected_transitive_sources: Dict[str, str] = {}
        expected_transitive_tokens: Dict[str, List[str]] = {}
        for relative in CORE_BUILD_SOURCES:
            output = "CMakeFiles/wirehair.dir/" + relative + ".o"
            source = str(source_root / relative)
            expected_transitive_sources[output] = source
            expected_transitive_tokens[output] = [
                str(CXX_COMPILER), "-DWIREHAIR_BUILDING=1",
                "-DWIREHAIR_STATIC=1", "-I" + str(source_root / "include"),
                "-O3", "-DNDEBUG", "-std=gnu++11", "-fPIC", "-Wall",
                "-Wextra", "-Wpedantic", "-Werror", "-o", output, "-c",
                source,
            ]
        for relative in TIMING_POLICY_BUILD_SOURCES:
            output = (
                "codec/CMakeFiles/wirehair_v2_timing_policy.dir/"
                + Path(relative).name + ".o"
            )
            source = str(source_root / relative)
            expected_transitive_sources[output] = source
            expected_transitive_tokens[output] = [
                str(CXX_COMPILER), "-DWIREHAIR_STATIC=1",
                "-DWIREHAIR_V2_ENABLE_BENCHMARK_EQUATIONS=1",
                "-I" + str(source_root / "codec"),
                "-I" + str(source_root / "include"), "-O3", "-DNDEBUG",
                "-std=gnu++11", "-Wall", "-Wextra", "-Wpedantic",
                "-Werror", "-o", output, "-c", source,
            ]
        if any(
            compile_outputs.get(output) != source
            or compile_tokens.get(output) != expected_transitive_tokens[output]
            for output, source in expected_transitive_sources.items()
        ):
            fail("build-profile transitive compile command differs")

        forbidden = (
            b"-march", b"-mtune", b"-mcpu", b"-flto", b"-fprofile",
            b"-fsanitize", b"-fno-sanitize", b"-fplugin", b"-specs",
            b"--sysroot", b"-isysroot", b"-finstrument", b"--coverage",
            b"-fcoverage", b"-fprofile-arcs", b"-ftest-coverage",
            b"wirehair_v2_test_hook",
        )
        for relative in ("build.ninja", "CMakeFiles/rules.ninja"):
            lowered = payloads[relative].lower()
            if any(token in lowered for token in forbidden):
                fail("build-profile generated flags differ: " + relative)
            optimization_tokens = re.findall(
                rb"(?<![A-Za-z0-9_])(-O[^ \t\r\n$]*)", payloads[relative],
            )
            if any(token != b"-O3" for token in optimization_tokens):
                fail("build-profile generated optimization differs: " + relative)
        build_ninja = payloads["build.ninja"]
        rules_ninja = payloads["CMakeFiles/rules.ninja"]
        try:
            build_text = build_ninja.decode("utf-8")
            rules_text = rules_ninja.decode("utf-8")
        except UnicodeDecodeError as exc:
            fail("build-profile Ninja files are not UTF-8: " + exception_text(exc))
        if (
            build_text.count("# Configurations: Release\n") != 1
            or build_text.count("CONFIGURATION = Release\n") != 1
            or build_text.count(
                "cmake_ninja_workdir = " + str(build) + "/\n"
            ) != 1
        ):
            fail("build-profile Ninja configuration differs")
        object_rows = [
            line for line in build_text.splitlines()
            if line.startswith("build ") and ".o:" in line
        ]
        if len(object_rows) != len(compile_outputs):
            fail("build-profile Ninja/compile-command object roster differs")
        for output, source in compile_outputs.items():
            matches = [
                line for line in object_rows
                if line.startswith("build " + output + ":")
            ]
            if len(matches) != 1 or (" " + source) not in matches[0]:
                fail("build-profile Ninja object source binding differs")
        for entry in expected_entries:
            output = entry["output"]
            source = entry["file"]
            object_block = "\n".join((
                "build " + output + ": CXX_COMPILER__" + BINARY_NAME
                + "_unscanned_Release " + source
                + " || cmake_object_order_depends_target_" + BINARY_NAME,
                "  DEFINES = -DWIREHAIR_STATIC=1"
                + " -DWIREHAIR_V2_ENABLE_BENCHMARK_EQUATIONS=1"
                + " -DWIREHAIR_WH2_SOURCE_GIT_COMMIT=\\\""
                + expected_commit + "\\\"",
                "  DEP_FILE = " + output + ".d",
                "  FLAGS = -O3 -DNDEBUG -std=gnu++11 -Wall -Wextra"
                + " -Wpedantic -Werror",
                "  INCLUDES = -I" + str(source_root / "codec")
                + " -I" + str(source_root / "include"),
                "  OBJECT_DIR = CMakeFiles/" + BINARY_NAME + ".dir",
                "  OBJECT_FILE_DIR = CMakeFiles/" + BINARY_NAME + ".dir/bench",
                "  TARGET_COMPILE_PDB = CMakeFiles/" + BINARY_NAME + ".dir/",
                "  TARGET_PDB = " + BINARY_NAME + ".pdb",
                "",
            ))
            if build_text.count(object_block + "\n") != 1:
                fail("build-profile Ninja target object differs")
        for variable in ("LAUNCHER", "CODE_CHECK", "LINK_FLAGS", "LINK_PATH"):
            assignment = re.compile(
                r"(?m)^(?:  )?" + re.escape(variable) + r"[ \t]*=",
            )
            if assignment.search(build_text) or assignment.search(rules_text):
                fail("build-profile Ninja expansion variable differs: " + variable)
        if build_ninja.count(("build " + BINARY_NAME + ":").encode("ascii")) != 1:
            fail("build-profile Ninja target link differs")
        object_outputs = [entry["output"] for entry in expected_entries]
        expected_link = (
            "build " + BINARY_NAME + ": CXX_EXECUTABLE_LINKER__" + BINARY_NAME
            + "_Release " + " ".join(object_outputs)
            + " | codec/libwirehair_v2_timing_policy.a libwirehair.a"
            + " || codec/libwirehair_v2_timing_policy.a libwirehair.a"
        )
        if build_text.count(expected_link + "\n") != 1:
            fail("build-profile exact direct link graph differs")
        link_start = build_text.index(expected_link + "\n")
        link_end = build_text.find("\n\n", link_start)
        link_block = build_text[link_start:link_end]
        link_lines = link_block.splitlines()
        expected_link_lines = [
            expected_link,
            "  FLAGS = -O3 -DNDEBUG",
            "  LINK_LIBRARIES = codec/libwirehair_v2_timing_policy.a  libwirehair.a  -lm",
            "  OBJECT_DIR = CMakeFiles/" + BINARY_NAME + ".dir",
            "  POST_BUILD = :", "  PRE_LINK = :",
            "  TARGET_COMPILE_PDB = CMakeFiles/" + BINARY_NAME + ".dir/",
            "  TARGET_FILE = " + BINARY_NAME,
            "  TARGET_PDB = " + BINARY_NAME + ".pdb",
        ]
        if link_lines != expected_link_lines:
            fail("build-profile direct link property roster differs")
        compile_rule = "rule CXX_COMPILER__" + BINARY_NAME \
            + "_unscanned_Release"
        link_rule = "rule CXX_EXECUTABLE_LINKER__" + BINARY_NAME + "_Release"
        if rules_text.count(compile_rule + "\n") != 1 \
                or rules_text.count(link_rule + "\n") != 1:
            fail("build-profile Ninja rule roster differs")
        compile_rule_block = rules_text.split(compile_rule + "\n", 1)[1].split(
            "\n\n", 1,
        )[0]
        link_rule_block = rules_text.split(link_rule + "\n", 1)[1].split(
            "\n\n", 1,
        )[0]
        expected_compile_command = (
            "  command = ${LAUNCHER}${CODE_CHECK}" + str(CXX_COMPILER)
            + " $DEFINES $INCLUDES $FLAGS -MD -MT $out -MF $DEP_FILE"
            + " -o $out -c $in"
        )
        expected_link_command = (
            "  command = $PRE_LINK && " + str(CXX_COMPILER)
            + " $FLAGS $LINK_FLAGS $in -o $TARGET_FILE $LINK_PATH"
            + " $LINK_LIBRARIES && $POST_BUILD"
        )
        if compile_rule_block.splitlines() != [
            "  depfile = $DEP_FILE", "  deps = gcc", expected_compile_command,
            "  description = Building CXX object $out",
        ] or link_rule_block.splitlines() != [
            expected_link_command,
            "  description = Linking CXX executable $TARGET_FILE",
            "  restat = $RESTAT",
        ]:
            fail("build-profile Ninja rule command differs")

        effective_commands = ninja_target_commands(build)
        validate_effective_ninja_commands(
            effective_commands, compile_tokens, object_outputs,
        )

        c_compiler = compiler_image_receipt(C_COMPILER, "c")
        cxx_compiler = compiler_image_receipt(CXX_COMPILER, "cxx")
        toolchain_images = frozen_toolchain_receipts()

        # B-side of the snapshot: semantic parsing and compiler queries all
        # happened while every generated artifact FD remained held.
        for relative, fd in opened.items():
            after = os.fstat(fd)
            named = os.stat(str(build / relative), follow_symlinks=False)
            first = before[relative]
            if (
                (first.st_dev, first.st_ino, first.st_size, first.st_mtime_ns,
                 first.st_ctime_ns,
                 first.st_mode, first.st_uid, first.st_gid, first.st_nlink)
                != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns,
                    after.st_ctime_ns,
                    after.st_mode, after.st_uid, after.st_gid, after.st_nlink)
                or (after.st_dev, after.st_ino) != (named.st_dev, named.st_ino)
            ):
                fail("build-profile artifact changed during validation: " + relative)

        manifest: Dict[str, Any] = {
            "entries": [receipts[name] for name in BUILD_PATHS],
            "root": str(build),
            "schema": "wirehair.wh2.direct-systematic-complement-build.v2",
            "sha256": None,
        }
        manifest["sha256"] = sha256_bytes(canonical_bytes(manifest))
        profile = self_hashed(
            "wirehair.wh2.direct-systematic-complement-build-profile.v1", {
                "artifact_bindings": [receipts[name] for name in BUILD_PATHS],
                "build_environment_contract": {
                    "exact_environment_keys": [
                        "HOME", "LANG", "LC_ALL", "PATH",
                        "PYTHONDONTWRITEBYTECODE", "TMPDIR", "TZ",
                    ],
                    "fixed_values": {
                        "LANG": "C", "LC_ALL": "C",
                        "PATH": "/usr/bin:/bin",
                        "PYTHONDONTWRITEBYTECODE": "1", "TZ": "UTC",
                    },
                    "home_policy": (
                        "dedicated-empty-campaign-owned-directory-crossbound-by-"
                        "root-build-receipt"
                    ),
                    "tmpdir_policy": (
                        "dedicated-empty-campaign-owned-directory-crossbound-by-"
                        "root-build-receipt"
                    ),
                    "inheritance": "env-i-no-additional-variables",
                },
                "c_compiler": c_compiler, "cache": {
                    name: {"type": value[0], "value": value[1]}
                    for name, value in sorted(required_cache.items())
                },
                "cxx_compiler": cxx_compiler,
                "expected_source_commit": expected_commit,
                "generated_forbidden_flags_absent": True,
                "target_compile_commands": expected_entries,
                "target_effective_commands": effective_commands,
                "toolchain_images": toolchain_images,
            },
        )
        return manifest["sha256"], profile
    finally:
        for fd in opened.values():
            os.close(fd)


def require_normal_git_index_flags(raw: bytes) -> None:
    if (
        not raw or not raw.endswith(b"\n")
        or any(not line.startswith(b"H ") for line in raw.splitlines())
    ):
        fail("static Git tracked index contains non-normal flags")


def parse_git_relative_path(raw: bytes, where: str) -> str:
    try:
        relative = raw.decode("ascii")
    except UnicodeDecodeError as exc:
        fail(where + " path is not ASCII: " + exception_text(exc))
    components = relative.split("/")
    if (
        not relative or len(relative) > 4096 or relative.startswith("/")
        or re.fullmatch(r"[A-Za-z0-9_.+@=/\-]+", relative) is None
        or any(component in ("", ".", "..") for component in components)
    ):
        fail(where + " path differs")
    return relative


def parse_git_blob_roster(raw: bytes, index: bool) -> Dict[str, Tuple[str, str]]:
    where = "static Git index" if index else "static Git HEAD tree"
    if not raw or not raw.endswith(b"\0"):
        fail(where + " framing differs")
    records = raw[:-1].split(b"\0")
    if not records or len(records) > MAX_TRACKED_SOURCE_FILES:
        fail(where + " file-count bound differs")
    result: Dict[str, Tuple[str, str]] = {}
    for record in records:
        try:
            metadata, raw_path = record.split(b"\t", 1)
            fields = metadata.split(b" ")
            if len(fields) != 3:
                fail(where + " metadata shape differs")
            raw_mode, raw_oid, raw_kind_or_stage = fields
            if index:
                raw_stage = raw_kind_or_stage
                if raw_stage != b"0":
                    fail(where + " contains a nonzero stage")
            else:
                raw_kind = raw_oid
                raw_oid = raw_kind_or_stage
                if raw_kind != b"blob":
                    fail(where + " contains a non-blob entry")
            mode = raw_mode.decode("ascii")
            oid = raw_oid.decode("ascii")
        except (UnicodeDecodeError, ValueError) as exc:
            fail(where + " row is malformed: " + exception_text(exc))
        if mode not in ("100644", "100755") or LOWER40.fullmatch(oid) is None:
            fail(where + " blob identity differs")
        relative = parse_git_relative_path(raw_path, where)
        if relative in result:
            fail(where + " contains a duplicate path")
        result[relative] = (mode, oid)
    return result


def hash_tracked_git_blob(root: Path, relative: str, mode: str) -> Tuple[str, int]:
    path = root / relative
    try:
        if path.resolve(strict=True) != path:
            fail("static Git worktree path is not canonical: " + relative)
    except OSError as exc:
        fail("static Git worktree path is unavailable: " + exception_text(exc))
    flags = os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK \
        | getattr(os, "O_NOFOLLOW", 0)
    try:
        fd = os.open(str(path), flags)
    except OSError as exc:
        fail("static Git worktree open failed: " + exception_text(exc))
    try:
        before = os.fstat(fd)
        expected_permissions = (0o555, 0o755) if mode == "100755" else (0o444, 0o644)
        if (
            not stat.S_ISREG(before.st_mode) or before.st_nlink != 1
            or before.st_uid != CAMPAIGN_UID or before.st_gid != CAMPAIGN_GID
            or stat.S_IMODE(before.st_mode) not in expected_permissions
            or before.st_size < 0
            or before.st_size > MAX_TRACKED_SOURCE_FILE_BYTES
        ):
            fail("static Git worktree file policy differs: " + relative)
        digest = hashlib.sha1()
        digest.update(b"blob " + str(before.st_size).encode("ascii") + b"\0")
        offset = 0
        while offset < before.st_size:
            block = os.pread(
                fd, min(1024 * 1024, before.st_size - offset), offset,
            )
            if not block:
                fail("static Git worktree read was short: " + relative)
            digest.update(block)
            offset += len(block)
        after = os.fstat(fd)
        named = os.stat(str(path), follow_symlinks=False)
        fields = (
            "st_dev", "st_ino", "st_size", "st_mtime_ns", "st_ctime_ns",
            "st_mode", "st_uid", "st_gid", "st_nlink",
        )
        before_values = tuple(getattr(before, field) for field in fields)
        if (
            before_values != tuple(getattr(after, field) for field in fields)
            or before_values != tuple(getattr(named, field) for field in fields)
        ):
            fail("static Git worktree file changed while hashing: " + relative)
        return digest.hexdigest(), before.st_size
    finally:
        os.close(fd)


def static_git_gate(root: Path, expected_commit: str) -> str:
    digest, info = hash_path(GIT, 64 * 1024 * 1024, "static Git executable")
    if (
        info.st_uid != 0 or info.st_gid != 0
        or stat.S_IMODE(info.st_mode) != 0o755
    ):
        fail("static Git executable policy differs")
    environment = {
        "GIT_ASKPASS": "/bin/false", "GIT_CONFIG_GLOBAL": "/dev/null",
        "GIT_CONFIG_NOSYSTEM": "1", "GIT_NO_LAZY_FETCH": "1",
        "GIT_NO_REPLACE_OBJECTS": "1", "GIT_OPTIONAL_LOCKS": "0",
        "GIT_TERMINAL_PROMPT": "0", "LANG": "C", "LC_ALL": "C",
        "PATH": "/usr/bin:/bin", "SSH_ASKPASS": "/bin/false",
    }

    def git(*arguments: str) -> bytes:
        result = subprocess.run(
            [
                str(GIT), "-c", "core.fsmonitor=false",
                "-c", "safe.directory=" + str(root), *arguments,
            ],
            cwd=str(root), env=environment, stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            timeout=10.0,
        )
        if result.returncode != 0 or result.stderr:
            fail("static Git command failed: " + " ".join(arguments))
        return result.stdout

    if git("rev-parse", "--show-toplevel") != (str(root) + "\n").encode("ascii"):
        fail("static Git worktree root differs")
    if git("rev-parse", "--verify", "HEAD^{commit}") != (
        expected_commit + "\n"
    ).encode("ascii"):
        fail("static Git HEAD differs")
    if git("status", "--porcelain=v1", "--untracked-files=all"):
        fail("static Git tracked worktree is dirty")
    if git("ls-files", "--others", "--exclude-standard"):
        fail("static Git worktree contains untracked files")
    if git("ls-files", "--others", "--ignored", "--exclude-standard"):
        fail("static Git worktree contains ignored files")
    require_normal_git_index_flags(git("ls-files", "-v"))
    head_roster = parse_git_blob_roster(
        git("ls-tree", "--full-tree", "-r", "-z", expected_commit), False,
    )
    index_roster = parse_git_blob_roster(git("ls-files", "--stage", "-z"), True)
    if index_roster != head_roster:
        fail("static Git index differs from HEAD")
    total_bytes = 0
    for relative in sorted(head_roster):
        mode, oid = head_roster[relative]
        observed_oid, observed_bytes = hash_tracked_git_blob(root, relative, mode)
        total_bytes += observed_bytes
        if total_bytes > MAX_TRACKED_SOURCE_BYTES:
            fail("static Git tracked byte bound differs")
        if observed_oid != oid:
            fail("static Git worktree bytes differ from HEAD: " + relative)
    expected_paths = b"".join(
        (relative + "\n").encode("ascii") for relative in STATIC_SOURCE_PATHS
    )
    if git("ls-files", "--error-unmatch", "--", *STATIC_SOURCE_PATHS) != expected_paths:
        fail("static Git source path roster differs")
    expected_flags = b"".join(
        ("H " + relative + "\n").encode("ascii")
        for relative in STATIC_SOURCE_PATHS
    )
    if git("ls-files", "-v", "--error-unmatch", "--", *STATIC_SOURCE_PATHS) \
            != expected_flags:
        fail("static Git source index flags differ")
    tree: Dict[str, str] = {}
    try:
        for raw_line in git(
            "ls-tree", "--full-tree", "-r", "HEAD", "--", *STATIC_SOURCE_PATHS
        ).decode("ascii").splitlines():
            metadata, relative = raw_line.split("\t", 1)
            mode, kind, oid = metadata.split(" ")
            if (
                mode != "100644" or kind != "blob"
                or LOWER40.fullmatch(oid) is None or relative in tree
            ):
                fail("static Git HEAD source blob row differs")
            tree[relative] = oid
    except (UnicodeDecodeError, ValueError) as exc:
        fail("static Git HEAD source blob roster malformed: " + exception_text(exc))
    if set(tree) != set(STATIC_SOURCE_PATHS):
        fail("static Git HEAD source blob roster differs")
    try:
        worktree_oids = git(
            "hash-object", "--no-filters", "--", *STATIC_SOURCE_PATHS
        ).decode("ascii").splitlines()
    except UnicodeDecodeError as exc:
        fail("static Git worktree blob roster malformed: " + exception_text(exc))
    if (
        len(worktree_oids) != len(STATIC_SOURCE_PATHS)
        or any(LOWER40.fullmatch(oid) is None for oid in worktree_oids)
        or any(
            tree[relative] != oid
            for relative, oid in zip(STATIC_SOURCE_PATHS, worktree_oids)
        )
    ):
        fail("static Git worktree bytes differ from HEAD")
    git("diff", "--quiet", "--no-ext-diff", "--no-textconv", "HEAD", "--",
        *STATIC_SOURCE_PATHS)
    after, after_info = hash_path(GIT, 64 * 1024 * 1024, "static Git executable after")
    if (
        after != digest
        or (after_info.st_dev, after_info.st_ino)
        != (info.st_dev, info.st_ino)
    ):
        fail("static Git executable changed")
    return digest


def require_root_owned_nonwritable_chain(path: Path) -> None:
    if not path.is_absolute() or path.resolve(strict=True) != path:
        fail("installed launcher path chain is not canonical")
    current = Path("/")
    components = [current]
    for part in path.parts[1:]:
        current = current / part
        components.append(current)
    for index, component in enumerate(components):
        info = os.stat(str(component), follow_symlinks=False)
        leaf = index == len(components) - 1
        if (
            info.st_uid != 0 or info.st_gid != 0
            or stat.S_IMODE(info.st_mode) & 0o022
            or (leaf and not stat.S_ISREG(info.st_mode))
            or (not leaf and not stat.S_ISDIR(info.st_mode))
        ):
            fail("installed launcher root path policy differs: " + str(component))


class AuthorityArgumentParser(argparse.ArgumentParser):
    def error(self, message: str) -> None:
        fail("argument validation: " + message)


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    tokens = list(argv)
    if tokens == ["--selftest"]:
        return argparse.Namespace(
            selftest=True, execute=False, seal_build=False,
        )
    if (
        len(tokens) == 7 and tokens[0] == "--seal-build"
        and tokens[1] == "--source-dir" and tokens[3] == "--build-dir"
        and tokens[5] == "--expected-source-commit"
        and not any("=" in token for token in tokens)
    ):
        parser = AuthorityArgumentParser(allow_abbrev=False, add_help=False)
        parser.add_argument("--seal-build", action="store_true", required=True)
        parser.add_argument("--source-dir", type=Path, required=True)
        parser.add_argument("--build-dir", type=Path, required=True)
        parser.add_argument("--expected-source-commit", required=True)
        args = parser.parse_args(tokens)
        args.selftest = False
        args.execute = False
        source = args.source_dir
        build = args.build_dir
        if (
            not source.is_absolute() or os.path.normpath(str(source)) != str(source)
            or not build.is_absolute() or os.path.normpath(str(build)) != str(build)
            or source == build or source in build.parents or build in source.parents
            or LOWER40.fullmatch(args.expected_source_commit) is None
        ):
            fail(
                "seal-build requires distinct canonical absolute source/build "
                "paths and one lowercase commit"
            )
        return args
    if (
        len(tokens) != 5 or tokens[0] != "--execute"
        or tokens[1] != "--build-dir"
        or tokens[3] != "--expected-source-commit"
        or any("=" in token for token in tokens)
    ):
        fail("argument vector does not match one exact public mode")
    parser = AuthorityArgumentParser(allow_abbrev=False, add_help=False)
    parser.add_argument("--execute", action="store_true", required=True)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--expected-source-commit", required=True)
    args = parser.parse_args(tokens)
    args.selftest = False
    args.seal_build = False
    build = args.build_dir
    if (
        not build.is_absolute() or os.path.normpath(str(build)) != str(build)
        or LOWER40.fullmatch(args.expected_source_commit) is None
    ):
        fail("execute requires one canonical absolute build path and lowercase commit")
    return args


def systemd_run_argv(build_dir: Path, expected_commit: str) -> List[str]:
    """Return the exact transient-service command used to enter execute mode."""
    if (
        not build_dir.is_absolute()
        or os.path.normpath(str(build_dir)) != str(build_dir)
        or LOWER40.fullmatch(expected_commit) is None
    ):
        fail("systemd launch inputs differ")
    properties = (
        "User=0", "Group=0", "SupplementaryGroups=",
        "Restart=no", "ExitType=cgroup", "KillMode=control-group",
        "SendSIGKILL=yes",
        "TimeoutStopSec=1s", "RuntimeMaxSec=240s", "Delegate=cpuset pids",
        "AllowedCPUs=120-123", "AllowedMemoryNodes=0", "CPUAffinity=123",
        "TasksMax=32", "LimitCORE=0", "UMask=0077", "PrivateTmp=no",
        "PrivateDevices=no", "ProtectControlGroups=no", "DevicePolicy=closed",
        "DeviceAllow=/dev/i2c-1 rw", "DeviceAllow=/dev/i2c-2 rw",
    )
    argv = [
        str(SYSTEMD_RUN), "--quiet", "--wait", "--pipe", "--collect",
        "--service-type=exec", "--expand-environment=no", "--unit=" + UNIT_NAME,
    ]
    for item in properties:
        argv.append("--property=" + item)
    argv.extend((
        str(ENV_EXECUTABLE), "-i",
        "LANG=" + SEALED_ENVIRONMENT["LANG"],
        "LC_ALL=" + SEALED_ENVIRONMENT["LC_ALL"],
        "PATH=" + SEALED_ENVIRONMENT["PATH"],
        "TZ=" + SEALED_ENVIRONMENT["TZ"],
        str(PYTHON), "-I", "-B", str(INSTALLED_LAUNCHER),
        "--execute", "--build-dir", str(build_dir),
        "--expected-source-commit", expected_commit,
    ))
    return argv


def build_systemd_run_argv(config: BuildSealConfig) -> List[str]:
    properties = (
        "User=0", "Group=0", "SupplementaryGroups=", "Restart=no",
        "KillMode=control-group", "SendSIGKILL=yes", "TimeoutStopSec=1s",
        "RuntimeMaxSec=2400s", "AllowedCPUs=0-1", "AllowedMemoryNodes=0",
        "CPUAffinity=0-1", "TasksMax=64", "LimitCORE=0", "UMask=0077",
        "PrivateTmp=yes", "PrivateDevices=yes", "DevicePolicy=closed",
        "ProtectControlGroups=yes", "ProtectProc=invisible",
    )
    argv = [
        str(SYSTEMD_RUN), "--quiet", "--wait", "--pipe", "--collect",
        "--service-type=exec", "--expand-environment=no",
        "--unit=" + BUILD_UNIT_NAME,
    ]
    for item in properties:
        argv.append("--property=" + item)
    argv.extend((
        str(ENV_EXECUTABLE), "-i",
        "LANG=" + SEALED_ENVIRONMENT["LANG"],
        "LC_ALL=" + SEALED_ENVIRONMENT["LC_ALL"],
        "PATH=" + SEALED_ENVIRONMENT["PATH"],
        "TZ=" + SEALED_ENVIRONMENT["TZ"],
        str(PYTHON), "-I", "-B", str(INSTALLED_LAUNCHER),
        "--seal-build", "--source-dir", str(config.source_dir),
        "--build-dir", str(config.build_dir),
        "--expected-source-commit", config.expected_commit,
    ))
    return argv


def build_environment(build_dir: Path) -> Dict[str, str]:
    environment = dict(BUILD_ENVIRONMENT_FIXED)
    environment["HOME"] = str(build_dir / ".wh2-build-home")
    environment["TMPDIR"] = str(build_dir / ".wh2-build-tmp")
    return environment


def build_command_vectors(config: BuildSealConfig) -> List[Tuple[str, List[str], float]]:
    source = str(config.source_dir)
    build = str(config.build_dir)
    configure = [
        "/usr/bin/cmake", "-S", source, "-B", build, "-G", "Ninja",
        "-DCMAKE_BUILD_TYPE=Release",
        "-DCMAKE_C_COMPILER=" + str(C_COMPILER),
        "-DCMAKE_CXX_COMPILER=" + str(CXX_COMPILER),
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
        "-DPython3_EXECUTABLE=" + str(PYTHON),
        "-DBUILD_TESTS=ON", "-DBUILD_CODEC_V2=ON",
        "-DBUILD_SHARED_LIBS=OFF", "-DMARCH_NATIVE=OFF",
        "-DWIREHAIR_BUILD_BOTH=OFF", "-DWIREHAIR_STATIC_PIC=ON",
        "-DWIREHAIR_BUILD_TOOLS=OFF",
        "-DWIREHAIR_BUILD_BENCHMARKS=OFF",
        "-DWIREHAIR_ENABLE_SCHEDULED_TESTS=OFF",
        "-DWIREHAIR_ENABLE_LIBFUZZER=OFF",
        "-DWIREHAIR_ENABLE_TRI_BMQ4_SCHUR_FALSIFIER=OFF",
        "-DWIREHAIR_STRICT_WARNINGS=ON", "-DWH_LTO=OFF",
        "-DWH_PGO_MODE=OFF",
    ]
    return [
        ("configure", configure, BUILD_CONFIGURE_SECONDS),
        (
            "build", ["/usr/bin/ninja", "-C", build, "-j", "2", "all"],
            BUILD_COMPILE_SECONDS,
        ),
        (
            "test", [str(CTEST), "--test-dir", build, "--output-on-failure",
                     "--no-tests=error", "-j", "2"],
            BUILD_TEST_SECONDS,
        ),
        (
            "binary-selftest", [str(config.build_dir / BINARY_NAME), "--selftest"],
            120.0,
        ),
        (
            "no-work", ["/usr/bin/ninja", "-C", build, "-n", "all"],
            BUILD_NO_WORK_SECONDS,
        ),
    ]


def parse_cmake_source_root(cache: bytes, build_dir: Path) -> Path:
    if not cache or len(cache) > 64 * 1024 * 1024 or b"\0" in cache:
        fail("CMake cache framing differs")
    try:
        lines = cache.decode("utf-8").splitlines()
    except UnicodeDecodeError as exc:
        fail("CMake cache is not UTF-8: " + exception_text(exc))
    prefix = "CMAKE_HOME_DIRECTORY:INTERNAL="
    values = [line[len(prefix):] for line in lines if line.startswith(prefix)]
    if len(values) != 1:
        fail("CMake cache source-root entry differs")
    root = Path(values[0])
    if (
        not root.is_absolute() or os.path.normpath(str(root)) != str(root)
        or not build_dir.is_absolute()
    ):
        fail("CMake source root is not canonical absolute")
    try:
        resolved = root.resolve(strict=True)
    except OSError as exc:
        fail("CMake source root cannot be resolved: " + exception_text(exc))
    if resolved != root or not root.is_dir():
        fail("CMake source root is not one canonical directory")
    return root


@dataclass(frozen=True)
class LaunchConfig:
    build_dir: Path
    expected_commit: str


@dataclass(frozen=True)
class BuildSealConfig:
    source_dir: Path
    build_dir: Path
    expected_commit: str


@dataclass(frozen=True)
class Prepared:
    config: LaunchConfig
    source_root: Path
    controller: Path
    sampler_script: Path
    binary: Path
    launcher_sha256: str
    cgroup_path: Path
    cgroup_relative: str
    controller_fd: int = -1
    controller_sha256: str = ""
    controller_bytes: int = 0
    sampler_fd: int = -1
    sampler_sha256: str = ""
    sampler_bytes: int = 0
    python_sha256: str = ""
    setpriv_sha256: str = ""
    env_sha256: str = ""
    git_sha256: str = ""
    binary_sha256: str = ""
    build_manifest_sha256: str = ""
    build_profile: Dict[str, Any] = field(default_factory=dict)
    build_authority: Dict[str, Any] = field(default_factory=dict)
    build_authority_sha256: str = ""
    source_manifest_sha256: str = ""
    i2c_devices: Tuple[Dict[str, Any], ...] = field(default_factory=tuple)
    i2c_pre_scan: Dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class PreflightAuthority:
    raw: bytes
    receipt: Dict[str, Any]
    run_argv: List[str]
    run_argv_sha256: str
    receipt_sha256: str


def argv_sha256(argv: Sequence[str]) -> str:
    if type(argv) is not list or any(type(item) is not str for item in argv):
        fail("run argv is not a string array")
    try:
        raw = b"".join(item.encode("ascii") + b"\0" for item in argv)
    except UnicodeEncodeError:
        fail("run argv is not ASCII")
    return sha256_bytes(raw)


def option_values(argv: Sequence[str], mode: str,
                  order: Sequence[str], where: str) -> Dict[str, str]:
    if (
        type(argv) is not list or len(argv) != 2 + 2 * len(order)
        or argv[1] != mode or any(type(item) is not str for item in argv)
        or any(
            argv[2 + 2 * index] != option
            for index, option in enumerate(order)
        )
    ):
        fail(where + " option roster/order differs")
    values: Dict[str, str] = {}
    for index, option in enumerate(order):
        item = argv[3 + 2 * index]
        if not item or len(item) > 4096 or "\0" in item:
            fail(where + " value differs for " + option)
        try:
            item.encode("ascii")
        except UnicodeEncodeError:
            fail(where + " value is not ASCII for " + option)
        values[option] = item
    return values


def parse_preflight_receipt(
    raw: bytes, prepared: Prepared, sampler_pid: int,
) -> PreflightAuthority:
    value = parse_canonical_line(raw, MAX_PREFLIGHT_BYTES, "preflight receipt")
    exact_keys(value, PREFLIGHT_KEYS, "preflight receipt")
    if value["schema"] != PREFLIGHT_SCHEMA:
        fail("preflight schema differs")
    expected_source = str(prepared.source_root)
    if (
        value["source_root"] != expected_source
        or value["expected_source_commit"] != prepared.config.expected_commit
    ):
        fail("preflight source authority differs")
    receipt_hash = exact_hash(value["receipt_sha256"], "preflight receipt")
    preimage = dict(value)
    preimage["receipt_sha256"] = None
    if sha256_bytes(canonical_bytes(preimage)) != receipt_hash:
        fail("preflight receipt self-hash differs")
    run_argv = value["run_argv"]
    values = option_values(run_argv, "--run-once", RUN_OPTION_ORDER, "run argv")
    digest = argv_sha256(run_argv)
    if digest != exact_hash(value["run_argv_sha256"], "preflight run argv"):
        fail("preflight run argv hash differs")
    expected_paths = {
        "--binary": str(prepared.binary),
        "--build-dir": str(prepared.config.build_dir),
        "--sampler-pid": str(sampler_pid),
        "--sampler-script": str(prepared.sampler_script),
        "--sampler-csv": str(SAMPLER_CSV),
        "--sampler-pid-file": str(SAMPLER_PID_FILE),
        "--sampler-validation-jsonl": str(SAMPLER_VALIDATION),
        "--sampler-receipt": str(SAMPLER_RECEIPT),
        "--expected-source-commit": prepared.config.expected_commit,
        "--cpu": str(WORKER_CPU),
        "--controller-cpu": str(CONTROLLER_CPU),
        "--sampler-cpu": str(SAMPLER_CPU),
        "--expected-binary-uid": str(CAMPAIGN_UID),
        "--expected-controller-uid": str(CAMPAIGN_UID),
        "--expected-controller-gid": str(CAMPAIGN_GID),
        "--expected-sampler-uid": str(CAMPAIGN_UID),
        "--expected-sampler-gid": str(CAMPAIGN_GID),
        "--expected-sampler-i2c-gid": str(I2C_GID),
    }
    if run_argv[0] != str(prepared.controller):
        fail("preflight controller path differs")
    for option, expected in expected_paths.items():
        if values[option] != expected:
            fail("preflight run authority differs for " + option)
    for option in (
        "--expected-binary-sha256", "--expected-build-manifest-sha256",
        "--expected-controller-sha256", "--expected-git-sha256",
        "--expected-python-sha256", "--expected-sampler-script-sha256",
        "--expected-sampler-cmdline-sha256",
        "--expected-sampler-environ-sha256",
        "--expected-sampler-executable-sha256",
        "--expected-source-manifest-sha256",
    ):
        exact_hash(values[option], "preflight " + option)
    if (
        values["--expected-controller-sha256"] != prepared.controller_sha256
        or values["--expected-sampler-script-sha256"] != prepared.sampler_sha256
        or values["--expected-python-sha256"] != prepared.python_sha256
        or values["--expected-git-sha256"] != prepared.git_sha256
        or values["--expected-binary-sha256"] != prepared.binary_sha256
        or values["--expected-build-manifest-sha256"]
        != prepared.build_manifest_sha256
        or values["--expected-source-manifest-sha256"]
        != prepared.source_manifest_sha256
    ):
        fail("preflight held image binding differs")
    for name, expected in (
        ("binary_before", prepared.binary_sha256),
        ("binary_after", prepared.binary_sha256),
        ("build_manifest_before", prepared.build_manifest_sha256),
        ("build_manifest_after", prepared.build_manifest_sha256),
        ("source_manifest_before", prepared.source_manifest_sha256),
        ("source_manifest_after", prepared.source_manifest_sha256),
    ):
        component = value[name]
        if type(component) is not dict or component.get("sha256") != expected:
            fail("preflight static receipt binding differs: " + name)
    source_manifest = value["source_manifest_after"]
    entries = source_manifest.get("entries") if type(source_manifest) is dict else None
    if (
        type(entries) is not list or len(entries) != len(STATIC_SOURCE_PATHS)
        or any(type(entry) is not dict for entry in entries)
        or [entry.get("path") for entry in entries] != list(STATIC_SOURCE_PATHS)
    ):
        fail("preflight source-manifest entry roster differs")
    source_hashes = {
        entry.get("path"): entry.get("sha256")
        for entry in entries if type(entry) is dict
    }
    if (
        source_hashes.get(SOURCE_LAUNCHER_RELATIVE) != prepared.launcher_sha256
        or source_hashes.get(CONTROLLER_RELATIVE) != prepared.controller_sha256
        or source_hashes.get(SAMPLER_RELATIVE) != prepared.sampler_sha256
    ):
        fail("preflight source manifest/held launcher images differ")
    for option in (
        "--expected-sampler-process-start-ticks",
        "--expected-sampler-csv-device", "--expected-sampler-csv-inode",
        "--expected-sampler-pid-file-device",
        "--expected-sampler-pid-file-inode",
        "--expected-sampler-validation-device",
        "--expected-sampler-validation-inode",
        "--expected-sampler-receipt-device",
        "--expected-sampler-receipt-inode",
    ):
        token = values[option]
        if CANONICAL_UINT.fullmatch(token) is None or int(token) == 0:
            fail("preflight integer authority differs for " + option)
    return PreflightAuthority(
        raw=raw, receipt=value, run_argv=list(run_argv),
        run_argv_sha256=digest, receipt_sha256=receipt_hash,
    )


def preflight_argv(prepared: Prepared, sampler_pid: int) -> List[str]:
    argv = [str(prepared.controller), "--preflight-seal"]
    values = {
        "--binary": str(prepared.binary),
        "--build-dir": str(prepared.config.build_dir),
        "--sampler-pid": str(sampler_pid),
        "--sampler-script": str(prepared.sampler_script),
        "--sampler-csv": str(SAMPLER_CSV),
        "--sampler-pid-file": str(SAMPLER_PID_FILE),
        "--sampler-validation-jsonl": str(SAMPLER_VALIDATION),
        "--sampler-receipt": str(SAMPLER_RECEIPT),
        "--expected-source-commit": prepared.config.expected_commit,
    }
    for option in PREFLIGHT_OPTION_ORDER:
        argv.extend((option, values[option]))
    return argv


def sampler_argv(prepared: Prepared, sampler_hash: str) -> List[str]:
    exact_hash(sampler_hash, "sampler source")
    return [
        str(prepared.sampler_script),
        "--csv", str(SAMPLER_CSV),
        "--pid-file", str(SAMPLER_PID_FILE),
        "--validation-jsonl", str(SAMPLER_VALIDATION),
        "--receipt", str(SAMPLER_RECEIPT),
        "--expected-source-sha256", sampler_hash,
        "--expected-output-owner-uid", str(CAMPAIGN_UID),
        "--interval", "1.0", "--dimm-attempts", "5",
        "--dimm-retry-delay", "0.01",
    ]


SEALED_SOURCE_FD = 198
SEALED_MODULE_BOOTSTRAP = (
    "import hashlib,os,sys,types;"
    "f=int(sys.argv[1]);p=sys.argv[2];h=sys.argv[3];n=int(sys.argv[4]);"
    "e=sys.argv[5];a=[p]+sys.argv[6:];"
    "d=b'';"
    "\nwhile len(d)<=n:"
    "\n b=os.pread(f,min(1048576,n+1-len(d)),len(d));"
    "\n if not b:break"
    "\n d+=b"
    "\nos.close(f);"
    "\nif len(d)!=n or hashlib.sha256(d).hexdigest()!=h:"
    "\n raise SystemExit('sealed source binding differs')"
    "\nm=types.ModuleType('_wh2_sealed_'+h);m.__file__=p;m.__package__=None;"
    "m.__dict__['__builtins__']=__builtins__;sys.modules[m.__name__]=m;"
    "sys.argv=a;exec(compile(d,p,'exec'),m.__dict__);"
    "raise SystemExit(getattr(m,e)(a[1:]))"
)


def setpriv_exec_argv(
    python_argv: Sequence[str], *, sampler: bool,
    sealed_source: Optional[Tuple[int, str, str, int, str]] = None,
) -> List[str]:
    if not python_argv or any(type(item) is not str or not item for item in python_argv):
        fail("role Python argv differs")
    group_args = ["--groups", str(I2C_GID)] if sampler else ["--clear-groups"]
    python_flags = ["-I", "-S", "-B"] if sampler else ["-I", "-B"]
    invocation = list(python_argv)
    if sealed_source is not None:
        fd, path, digest, size, entrypoint = sealed_source
        if (
            sampler or fd != SEALED_SOURCE_FD or not path.startswith("/")
            or LOWER64.fullmatch(digest) is None or size < 1
            or entrypoint not in ("main",)
            or python_argv[0] != path
        ):
            fail("sealed source execution authority differs")
        invocation = [
            "-c", SEALED_MODULE_BOOTSTRAP, str(fd), path, digest, str(size),
            entrypoint, *list(python_argv[1:]),
        ]
        python_flags = ["-I", "-B"]
    return [
        str(SETPRIV), "--reuid", str(CAMPAIGN_UID),
        "--regid", str(CAMPAIGN_GID), *group_args,
        "--bounding-set=-all", "--inh-caps=-all", "--ambient-caps=-all",
        "--no-new-privs", "--pdeathsig", "SIGTERM", str(PYTHON),
        *python_flags, *invocation,
    ]


def role_document(role: str, python_argv: Sequence[str],
                  authority: PreflightAuthority,
                  prepared: Optional[Prepared] = None) -> bytes:
    if role not in ("controller", "verifier"):
        fail("role command differs")
    value: Dict[str, Any] = {
        "environment": dict(SEALED_ENVIRONMENT),
        "expected_preflight_sha256": sha256_bytes(authority.raw),
        "expected_run_argv_sha256": authority.run_argv_sha256,
        "held_controller_bytes": (
            prepared.controller_bytes if prepared is not None else 1
        ),
        "held_controller_path": (
            str(prepared.controller) if prepared is not None else python_argv[0]
        ),
        "held_controller_sha256": (
            prepared.controller_sha256 if prepared is not None else "0" * 64
        ),
        "gid": CAMPAIGN_GID, "groups": [], "no_new_privs": True,
        "python_argv": list(python_argv), "receipt_sha256": None,
        "role": role, "schema": ROLE_SCHEMA, "uid": CAMPAIGN_UID,
    }
    value["receipt_sha256"] = sha256_bytes(canonical_bytes(value))
    return canonical_bytes(value) + b"\n"


def validate_role_document(raw: bytes, role: str,
                           preflight_raw: bytes) -> Dict[str, Any]:
    value = parse_canonical_line(raw, MAX_PREFLIGHT_BYTES, "role command")
    exact_keys(value, {
        "environment", "expected_preflight_sha256",
        "expected_run_argv_sha256", "gid", "groups",
        "held_controller_bytes", "held_controller_path",
        "held_controller_sha256", "no_new_privs",
        "python_argv", "receipt_sha256", "role", "schema", "uid",
    }, "role command")
    if (
        value["schema"] != ROLE_SCHEMA or value["role"] != role
        or value["uid"] != CAMPAIGN_UID or value["gid"] != CAMPAIGN_GID
        or value["groups"] != [] or value["no_new_privs"] is not True
        or value["environment"] != SEALED_ENVIRONMENT
        or value["expected_preflight_sha256"] != sha256_bytes(preflight_raw)
    ):
        fail("role command authority differs")
    exact_hash(value["expected_run_argv_sha256"], "role run argv")
    exact_hash(value["held_controller_sha256"], "role controller image")
    exact_uint(value["held_controller_bytes"], 1, 1024 * 1024,
               "role controller bytes")
    if (
        type(value["held_controller_path"]) is not str
        or not value["held_controller_path"].startswith("/")
        or value["held_controller_path"] != value["python_argv"][0]
    ):
        fail("role held controller path differs")
    receipt_hash = exact_hash(value["receipt_sha256"], "role receipt")
    preimage = dict(value)
    preimage["receipt_sha256"] = None
    if sha256_bytes(canonical_bytes(preimage)) != receipt_hash:
        fail("role command self-hash differs")
    argv = value["python_argv"]
    if (
        type(argv) is not list or not argv
        or any(type(item) is not str or not item or len(item) > 4096 for item in argv)
    ):
        fail("role Python argv differs")
    return value


def parse_proc_cgroup(raw: bytes, expected_unit: str = UNIT_NAME) -> str:
    if not raw or len(raw) > 64 * 1024 or b"\0" in raw or b"\r" in raw:
        fail("/proc/self/cgroup framing differs")
    try:
        lines = raw.decode("ascii").splitlines()
    except UnicodeDecodeError:
        fail("/proc/self/cgroup is not ASCII")
    if len(lines) != 1 or not lines[0].startswith("0::/"):
        fail("process is not in one unified cgroup-v2 hierarchy")
    relative = lines[0][3:]
    if (
        not relative.startswith("/") or "//" in relative
        or any(part in ("", ".", "..") for part in relative.split("/")[1:])
        or Path(relative).name != expected_unit
    ):
        fail("delegated service cgroup path differs")
    return relative


def proc_start_ticks(raw: bytes) -> int:
    if raw.endswith(b"\n"):
        raw = raw[:-1]
    if not raw or len(raw) > 64 * 1024 or b"\0" in raw or b"\n" in raw:
        fail("process stat framing differs")
    right = raw.rfind(b")")
    if right < 2 or right + 2 >= len(raw) or raw[right + 1:right + 2] != b" ":
        fail("process stat comm framing differs")
    fields = raw[right + 2:].split()
    # fields[0] is kernel field 3 (state); starttime is kernel field 22.
    if len(fields) < 20:
        fail("process stat field roster differs")
    token = fields[19]
    if not token.isdigit() or token.startswith(b"0"):
        fail("process start ticks differ")
    value = int(token)
    if value <= 0 or value > (1 << 64) - 1:
        fail("process start ticks outside uint64")
    return value


def i2c_device_receipts() -> Tuple[Dict[str, Any], ...]:
    result: List[Dict[str, Any]] = []
    for path in (Path("/dev/i2c-1"), Path("/dev/i2c-2")):
        info = os.stat(str(path), follow_symlinks=False)
        if (
            not stat.S_ISCHR(info.st_mode) or info.st_nlink != 1
            or info.st_uid != 0 or info.st_gid != I2C_GID
            or stat.S_IMODE(info.st_mode) != 0o660
        ):
            fail("I2C device-node policy differs: " + str(path))
        result.append({
            "device": info.st_dev, "gid": info.st_gid,
            "inode": info.st_ino, "mode": stat.S_IMODE(info.st_mode),
            "path": str(path), "rdev": info.st_rdev, "uid": info.st_uid,
        })
    return tuple(result)


def scan_i2c_holders(devices: Sequence[Mapping[str, Any]]) -> Dict[str, Any]:
    identities = {
        (item["device"], item["inode"], item["rdev"]): item["path"]
        for item in devices
    }
    if len(identities) != 2:
        fail("I2C device identity roster differs")
    holders: Dict[str, List[Dict[str, int]]] = {
        item["path"]: [] for item in devices
    }
    try:
        processes = sorted(
            int(name) for name in os.listdir("/proc") if name.isdigit()
        )
    except OSError as exc:
        fail("I2C /proc scan failed: " + exception_text(exc))
    if len(processes) > 1_000_000:
        fail("I2C /proc process roster exceeds bound")
    observed_fds = 0
    for pid in processes:
        directory = Path("/proc") / str(pid) / "fd"
        try:
            names = os.listdir(str(directory))
        except FileNotFoundError:
            continue
        except PermissionError as exc:
            fail("I2C /proc descriptor scan is not authoritative: " + exception_text(exc))
        observed_fds += len(names)
        if observed_fds > 1_000_000:
            fail("I2C /proc descriptor roster exceeds bound")
        for name in names:
            if not name.isdigit():
                fail("I2C /proc descriptor name differs")
            try:
                info = os.stat(str(directory / name))
            except FileNotFoundError:
                continue
            identity = (info.st_dev, info.st_ino, info.st_rdev)
            path = identities.get(identity)
            if path is None:
                continue
            try:
                stat_fd = os.open(
                    f"/proc/{pid}/stat", os.O_RDONLY | os.O_CLOEXEC
                    | os.O_NONBLOCK | getattr(os, "O_NOFOLLOW", 0)
                )
                try:
                    raw_stat = os.read(stat_fd, 64 * 1024 + 1)
                finally:
                    os.close(stat_fd)
                if len(raw_stat) > 64 * 1024:
                    fail("I2C holder process stat exceeds bound")
                start = proc_start_ticks(raw_stat)
            except FileNotFoundError:
                continue
            holders[path].append({
                "fd": int(name), "pid": pid, "process_start_ticks": start,
            })
    return {
        "devices": [
            {
                **dict(device),
                "holders": sorted(
                    holders[device["path"]],
                    key=lambda item: (item["pid"], item["fd"]),
                ),
            }
            for device in devices
        ],
        "observed_monotonic_ns": time.monotonic_ns(),
        "scope": "trusted-root-point-in-time-no-kernel-exclusion",
    }


def require_no_i2c_holders(receipt: Mapping[str, Any]) -> None:
    if (
        receipt.get("scope")
        != "trusted-root-point-in-time-no-kernel-exclusion"
        or type(receipt.get("devices")) is not list
        or any(device.get("holders") != [] for device in receipt["devices"])
    ):
        fail("pre-sampler I2C holder roster is not empty")


def require_sampler_sole_i2c_holder(receipt: Mapping[str, Any],
                                    sampler: ProcessHandle) -> None:
    if (
        receipt.get("scope")
        != "trusted-root-point-in-time-no-kernel-exclusion"
        or type(receipt.get("devices")) is not list
        or len(receipt["devices"]) != 2
    ):
        fail("post-open I2C holder receipt differs")
    expected = [{
        "fd": None, "pid": sampler.pid,
        "process_start_ticks": sampler.start_ticks,
    }]
    for device in receipt["devices"]:
        holders = device.get("holders")
        if type(holders) is not list or len(holders) != 1:
            fail("sampler is not the sole I2C device holder")
        holder = holders[0]
        if (
            type(holder) is not dict or set(holder) != set(expected[0])
            or type(holder["fd"]) is not int or holder["fd"] < 0
            or holder["pid"] != sampler.pid
            or holder["process_start_ticks"] != sampler.start_ticks
        ):
            fail("sampler I2C holder identity differs")


def parse_claim(raw: bytes, authority: PreflightAuthority, controller_pid: int,
                controller_start_ticks: int, release_ns: int,
                observed_ns: int) -> Tuple[Dict[str, Any], int]:
    value = parse_canonical_line(raw, MAX_CLAIM_BYTES, "campaign claim")
    exact_keys(value, {
        "binary_sha256", "build_manifest_sha256", "campaign",
        "controller_receipt_sha256", "controller_started_monotonic_ns",
        "gid", "git_sha256", "output_path", "parent_device", "parent_inode",
        "pid", "process_start_ticks", "schema", "source_commit",
        "source_manifest_sha256", "uid",
    }, "campaign claim")
    run = option_values(authority.run_argv, "--run-once", RUN_OPTION_ORDER,
                        "claim run authority")
    if (
        value["schema"] != CLAIM_SCHEMA or value["campaign"] != CAMPAIGN
        or value["output_path"] != str(FIXED_OUTPUT_DIR)
        or value["uid"] != CAMPAIGN_UID or value["gid"] != CAMPAIGN_GID
        or value["pid"] != controller_pid
        or value["process_start_ticks"] != controller_start_ticks
        or value["source_commit"] != run["--expected-source-commit"]
        or value["binary_sha256"] != run["--expected-binary-sha256"]
        or value["build_manifest_sha256"]
        != run["--expected-build-manifest-sha256"]
        or value["git_sha256"] != run["--expected-git-sha256"]
        or value["source_manifest_sha256"]
        != run["--expected-source-manifest-sha256"]
    ):
        fail("campaign claim authority differs")
    exact_hash(value["controller_receipt_sha256"], "claim controller receipt")
    started = exact_uint(
        value["controller_started_monotonic_ns"], 1, (1 << 63) - 1,
        "claim controller start",
    )
    if started < release_ns or started > observed_ns:
        fail("claim controller epoch is outside the release/observation sandwich")
    parent = os.stat("/var/tmp", follow_symlinks=False)
    if (
        value["parent_device"] != parent.st_dev
        or value["parent_inode"] != parent.st_ino
        or not stat.S_ISDIR(parent.st_mode)
        or stat.S_IMODE(parent.st_mode) != 0o1777
        or parent.st_uid != 0 or parent.st_gid != 0
    ):
        fail("claim /var/tmp authority differs")
    return value, started


def parse_verifier_output(raw: bytes, returncode: int) -> str:
    if returncode not in EXIT_VERIFY:
        fail("retained verifier did not return an authentic classification code")
    value = parse_canonical_line(raw, 4096, "retained verifier result")
    exact_keys(value, {"outcome", "schema", "status"}, "retained verifier result")
    outcome = EXIT_VERIFY[returncode]
    if (
        value["schema"] != VERIFY_SCHEMA or value["outcome"] != outcome
        or value["status"] != ("absent" if outcome == "absent" else "verified")
    ):
        fail("retained verifier output/exit binding differs")
    return outcome


def file_binding(fd: int, limit: int, *,
                 expected_nlink: Optional[int] = 1) -> Dict[str, Any]:
    before = os.fstat(fd)
    if (
        not stat.S_ISREG(before.st_mode)
        or (expected_nlink is not None and before.st_nlink != expected_nlink)
        or before.st_size < 0 or before.st_size > limit
    ):
        fail("held artifact policy differs")
    digest = hashlib.sha256()
    offset = 0
    while offset < before.st_size:
        block = os.pread(fd, min(1024 * 1024, before.st_size - offset), offset)
        if not block:
            fail("held artifact read was short")
        digest.update(block)
        offset += len(block)
    after = os.fstat(fd)
    if (
        before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns,
        before.st_ctime_ns, before.st_mode, before.st_uid, before.st_gid,
        before.st_nlink,
    ) != (
        after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns,
        after.st_ctime_ns, after.st_mode, after.st_uid, after.st_gid,
        after.st_nlink,
    ):
        fail("held artifact changed while hashing")
    return {
        "bytes": before.st_size, "device": before.st_dev, "gid": before.st_gid,
        "inode": before.st_ino, "mode": stat.S_IMODE(before.st_mode),
        "nlink": before.st_nlink, "sha256": digest.hexdigest(),
        "uid": before.st_uid,
    }


def read_held_file(fd: int, limit: int, where: str) -> bytes:
    """Read one already-held regular file without consulting its pathname."""
    before = os.fstat(fd)
    if (
        not stat.S_ISREG(before.st_mode) or before.st_nlink != 1
        or before.st_size < 1 or before.st_size > limit
    ):
        fail(where + " held-file policy differs")
    data = bytearray()
    offset = 0
    while offset < before.st_size:
        block = os.pread(fd, min(1024 * 1024, before.st_size - offset), offset)
        if not block:
            fail(where + " held-file read was short")
        data.extend(block)
        offset += len(block)
    after = os.fstat(fd)
    if (
        before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns,
        before.st_ctime_ns, before.st_mode, before.st_uid, before.st_gid,
        before.st_nlink,
    ) != (
        after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns,
        after.st_ctime_ns, after.st_mode, after.st_uid, after.st_gid,
        after.st_nlink,
    ):
        fail(where + " held-file changed while reading")
    return bytes(data)


def open_held_source(path: Path, where: str) -> Tuple[int, Dict[str, Any]]:
    fd = os.open(
        str(path), os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK
        | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        binding = file_binding(fd, 1024 * 1024)
        named = os.stat(str(path), follow_symlinks=False)
        if (
            (binding["device"], binding["inode"])
            != (named.st_dev, named.st_ino)
            or binding["mode"] != 0o444 or binding["uid"] != CAMPAIGN_UID
            or binding["gid"] != CAMPAIGN_GID
        ):
            fail(where + " held-source policy differs")
        return fd, binding
    except BaseException:
        os.close(fd)
        raise


def validate_sampler_terminal(
    raw: bytes, sampler_pid: int, prepared: Prepared,
    held: Mapping[str, int], sampler_returncode: int,
) -> Dict[str, Any]:
    value = parse_canonical_line(raw, 1024 * 1024, "sampler terminal receipt")
    if value.get("schema") != SAMPLER_SCHEMA or value.get("pid") != sampler_pid:
        fail("sampler terminal identity differs")
    if type(value.get("self_sha256_excluding_field")) is not str:
        fail("sampler terminal self-hash is absent")
    digest = exact_hash(
        value["self_sha256_excluding_field"], "sampler terminal receipt"
    )
    preimage = dict(value)
    del preimage["self_sha256_excluding_field"]
    if sha256_bytes(canonical_bytes(preimage) + b"\n") != digest:
        fail("sampler terminal self-hash differs")
    # The controller and verifier run as the same UID that owns the worktree.
    # Authenticate against the source FD retained by root before either role
    # ran.  A later pathname replacement is a launcher diagnostic, not a way
    # to erase otherwise authentic terminal evidence.
    source_binding = file_binding(
        prepared.sampler_fd, 1024 * 1024, expected_nlink=None,
    )
    if (
        source_binding["nlink"] not in (0, 1)
        or source_binding["bytes"] != prepared.sampler_bytes
        or source_binding["sha256"] != prepared.sampler_sha256
        or source_binding["mode"] != 0o444
        or source_binding["uid"] != CAMPAIGN_UID
        or source_binding["gid"] != CAMPAIGN_GID
    ):
        fail("held sampler source terminal binding differs")
    expected_sampler_hash = prepared.sampler_sha256
    source_path_still_bound = False
    try:
        named_source = os.stat(
            str(prepared.sampler_script), follow_symlinks=False,
        )
    except (FileNotFoundError, OSError):
        pass
    else:
        source_path_still_bound = (
            source_binding["nlink"] == 1
            and (source_binding["device"], source_binding["inode"])
            == (named_source.st_dev, named_source.st_ino)
            and stat.S_ISREG(named_source.st_mode)
        )
    expected_argv = sampler_argv(prepared, expected_sampler_hash)[1:]
    if (
        value.get("argv") != expected_argv
        or value.get("expected_output_owner_uid") != CAMPAIGN_UID
        or value.get("exit_code") != sampler_returncode
    ):
        fail("sampler terminal launch/exit binding differs")
    source = value.get("sampler_source")
    if (
        type(source) is not dict or source.get("path") != str(prepared.sampler_script)
        or source.get("expected_sha256") != expected_sampler_hash
        or source.get("sha256") != expected_sampler_hash
    ):
        fail("sampler terminal source binding differs")
    artifact_values = (
        ("raw_csv", "csv", SAMPLER_CSV),
        ("validation_jsonl", "validation", SAMPLER_VALIDATION),
        ("receipt_file", "receipt", SAMPLER_RECEIPT),
    )
    final: Dict[str, Any] = {}
    for receipt_name, held_name, path in artifact_values:
        fd = held[held_name]
        observed = file_binding(fd, MAX_SAMPLER_ARTIFACT_BYTES)
        named = os.stat(str(path), follow_symlinks=False)
        if (
            (observed["device"], observed["inode"])
            != (named.st_dev, named.st_ino)
            or observed["mode"] != 0o444 or observed["uid"] != CAMPAIGN_UID
            or observed["gid"] != CAMPAIGN_GID or observed["nlink"] != 1
        ):
            fail("sampler terminal artifact policy differs: " + held_name)
        record = value.get(receipt_name)
        binding = record.get("binding") if type(record) is dict else None
        if receipt_name == "receipt_file":
            # The receipt describes its initially empty reservation rather than
            # recursively hashing its final bytes.  Held identity is authoritative.
            binding = record.get("reservation_binding") if type(record) is dict else None
        if type(record) is not dict or record.get("path") != str(path):
            fail("sampler terminal artifact path differs: " + held_name)
        if receipt_name != "receipt_file" and (
            type(binding) is not dict
            or binding.get("device") != observed["device"]
            or binding.get("inode") != observed["inode"]
            or binding.get("sha256") != observed["sha256"]
            or binding.get("size") != observed["bytes"]
        ):
            fail("sampler terminal artifact content binding differs: " + held_name)
        if receipt_name == "receipt_file" and (
            type(binding) is not dict
            or binding.get("device") != observed["device"]
            or binding.get("inode") != observed["inode"]
        ):
            fail("sampler terminal receipt reservation identity differs")
        final[held_name] = observed
    pid_record = value.get("pid_file")
    pid_binding = file_binding(held["pid"], 4096, expected_nlink=0)
    if (
        type(pid_record) is not dict or pid_record.get("path") != str(SAMPLER_PID_FILE)
        or pid_record.get("removed") is not True
        or type(pid_record.get("binding")) is not dict
        or pid_record["binding"].get("device") != pid_binding["device"]
        or pid_record["binding"].get("inode") != pid_binding["inode"]
        or pid_binding["bytes"] != len(str(sampler_pid)) + 1
        or pid_binding["sha256"]
        != sha256_bytes((str(sampler_pid) + "\n").encode("ascii"))
        or pid_binding["mode"] != 0o444
        or pid_binding["uid"] != CAMPAIGN_UID
        or pid_binding["gid"] != CAMPAIGN_GID
        or pid_binding["nlink"] != 0
    ):
        fail("sampler terminal PID binding differs")
    try:
        os.stat(str(SAMPLER_PID_FILE), follow_symlinks=False)
    except FileNotFoundError:
        pass
    else:
        fail("sampler PID pathname survived terminal receipt")
    final["pid"] = pid_binding
    return {
        "artifacts": final, "receipt_sha256": sha256_bytes(raw),
        "source_path_still_bound": source_path_still_bound, "value": value,
    }


@dataclass
class ProcessHandle:
    role: str
    pid: int
    pidfd: int
    start_ticks: int
    stdout_fd: int = -1
    stderr_fd: int = -1
    control: Optional[socket.socket] = None
    returncode: Optional[int] = None


@dataclass
class CgroupLayout:
    service: Path
    supervisor: Path
    run: Path
    sampler: Path
    experiment: Path
    verifier: Path
    experiment_kill_fd: int
    sampler_kill_fd: int
    verifier_kill_fd: int
    run_kill_fd: int


class AttemptJournal:
    """O_EXCL-only durable root journal.  Nothing in it is ever replaced."""

    def __init__(self, parent_fd: int, directory_fd: int) -> None:
        self.parent_fd = parent_fd
        self.directory_fd = directory_fd
        self.directory_binding: Optional[Tuple[int, int]] = None
        self.names: List[str] = []
        self.complete = False

    def _ensure_directory(self) -> None:
        if self.directory_fd >= 0:
            return
        flags = (
            os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC
            | getattr(os, "O_NOFOLLOW", 0)
        )
        fd = os.open(LAUNCH_DIR.name, flags, dir_fd=self.parent_fd)
        try:
            directory = os.fstat(fd)
            named = os.stat(
                LAUNCH_DIR.name, dir_fd=self.parent_fd,
                follow_symlinks=False,
            )
            binding = (directory.st_dev, directory.st_ino)
            if (
                not stat.S_ISDIR(directory.st_mode)
                or binding != (named.st_dev, named.st_ino)
                or stat.S_IMODE(directory.st_mode) != 0o700
                or directory.st_uid != 0 or directory.st_gid != 0
                or (
                    self.directory_binding is not None
                    and binding != self.directory_binding
                )
            ):
                fail("attempt directory policy differs")
            self.directory_binding = binding
            self.directory_fd = fd
        except BaseException:
            os.close(fd)
            raise

    @classmethod
    def reserve(cls, attempt: Mapping[str, Any],
                consumed_callback: Optional[Any] = None) -> "AttemptJournal":
        parent_flags = os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC | getattr(os, "O_NOFOLLOW", 0)
        parent_fd = os.open("/var/tmp", parent_flags)
        consumed = False
        journal: Optional["AttemptJournal"] = None
        try:
            parent = os.fstat(parent_fd)
            if (
                not stat.S_ISDIR(parent.st_mode) or stat.S_IMODE(parent.st_mode) != 0o1777
                or parent.st_uid != 0 or parent.st_gid != 0
            ):
                fail("/var/tmp policy differs")
            try:
                os.stat(
                    LAUNCH_DIR.name, dir_fd=parent_fd,
                    follow_symlinks=False,
                )
            except FileNotFoundError:
                pass
            else:
                fail("attempt namespace already exists at reservation")
            # Allocate the recovery authority before the irrevocable mkdir.
            # After mkdir succeeds, even callback/allocation failures must be
            # able to propagate the live parent FD for terminal accounting.
            journal = cls(parent_fd, -1)
            os.mkdir(LAUNCH_DIR.name, 0o700, dir_fd=parent_fd)
            consumed = True
            if consumed_callback is not None:
                consumed_callback()
            journal._ensure_directory()
            os.fsync(parent_fd)
            journal.write_json("ATTEMPT", attempt)
            return journal
        except BaseException as exc:
            if not consumed and journal is not None:
                # This outer check also covers a signal/trace exception after
                # mkdir returned but before the following assignment executed.
                try:
                    os.stat(
                        LAUNCH_DIR.name, dir_fd=parent_fd,
                        follow_symlinks=False,
                    )
                except FileNotFoundError:
                    pass
                except BaseException:
                    consumed = True
                else:
                    consumed = True
            if consumed:
                raise AttemptConsumedError(
                    "attempt namespace consumed during reservation: "
                    + exception_text(exc), journal,
                ) from exc
            if parent_fd >= 0:
                os.close(parent_fd)
            raise

    def _reserve(self, name: str, mode: int = 0o600) -> int:
        if self.complete or not name or "/" in name or name in self.names:
            fail("attempt journal append order differs")
        self._ensure_directory()
        fd = os.open(
            name, os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
            | getattr(os, "O_NOFOLLOW", 0), mode, dir_fd=self.directory_fd,
        )
        self.names.append(name)
        return fd

    def write_bytes(self, name: str, payload: bytes, final_mode: int = 0o400) -> Dict[str, Any]:
        fd = self._reserve(name)
        try:
            view = memoryview(payload)
            while view:
                written = os.write(fd, view)
                if written <= 0:
                    fail("attempt journal write was short: " + name)
                view = view[written:]
            os.fchmod(fd, final_mode)
            os.fsync(fd)
            os.fsync(self.directory_fd)
            info = os.fstat(fd)
            named = os.stat(name, dir_fd=self.directory_fd, follow_symlinks=False)
            if (
                (info.st_dev, info.st_ino) != (named.st_dev, named.st_ino)
                or info.st_size != len(payload) or info.st_nlink != 1
                or stat.S_IMODE(info.st_mode) != final_mode
                or info.st_uid != 0 or info.st_gid != 0
            ):
                fail("attempt journal artifact policy differs: " + name)
            return {
                "bytes": len(payload), "device": info.st_dev, "inode": info.st_ino,
                "mode": final_mode, "sha256": sha256_bytes(payload),
            }
        finally:
            os.close(fd)

    def write_json(self, name: str, value: Mapping[str, Any]) -> Dict[str, Any]:
        return self.write_bytes(name, canonical_bytes(dict(value)) + b"\n")

    def finish(self) -> None:
        if self.complete or not self.names or self.names[-1] != "terminal.json":
            fail("attempt COMPLETE is not terminal-last")
        fd = self._reserve("COMPLETE")
        try:
            os.fchmod(fd, 0o000)
            os.fsync(fd)
            os.fchmod(self.directory_fd, 0o500)
            os.fsync(self.directory_fd)
            os.fsync(self.parent_fd)
            complete_info = os.fstat(fd)
            directory = os.fstat(self.directory_fd)
            named = os.stat(
                LAUNCH_DIR.name, dir_fd=self.parent_fd, follow_symlinks=False
            )
            if (
                not stat.S_ISREG(complete_info.st_mode)
                or complete_info.st_size != 0 or complete_info.st_nlink != 1
                or stat.S_IMODE(complete_info.st_mode) != 0
                or complete_info.st_uid != 0 or complete_info.st_gid != 0
                or not stat.S_ISDIR(directory.st_mode)
                or stat.S_IMODE(directory.st_mode) != 0o500
                or directory.st_uid != 0 or directory.st_gid != 0
                or (directory.st_dev, directory.st_ino)
                != (named.st_dev, named.st_ino)
            ):
                fail("attempt COMPLETE/directory policy differs")
            self.complete = True
        finally:
            os.close(fd)

    def close(self) -> None:
        for fd in (self.directory_fd, self.parent_fd):
            if fd >= 0:
                try:
                    os.close(fd)
                except OSError:
                    pass
        self.directory_fd = -1
        self.parent_fd = -1


def utc_now() -> str:
    return datetime.datetime.now(datetime.timezone.utc).isoformat(
        timespec="microseconds"
    ).replace("+00:00", "Z")


def self_hashed(schema: str, values: Mapping[str, Any]) -> Dict[str, Any]:
    result = {"schema": schema, **dict(values), "receipt_sha256": None}
    result["receipt_sha256"] = sha256_bytes(canonical_bytes(result))
    return result


def directory_binding(path: Path, uid: int, gid: int, mode: int,
                      where: str) -> Dict[str, Any]:
    info = os.stat(str(path), follow_symlinks=False)
    if (
        not stat.S_ISDIR(info.st_mode) or info.st_uid != uid or info.st_gid != gid
        or stat.S_IMODE(info.st_mode) != mode
    ):
        fail(where + " directory policy differs")
    return {
        "device": info.st_dev, "gid": info.st_gid, "inode": info.st_ino,
        "mode": stat.S_IMODE(info.st_mode), "path": str(path),
        "uid": info.st_uid,
    }


def installed_launcher_binding() -> Dict[str, Any]:
    installed = Path(__file__).resolve(strict=True)
    if installed != INSTALLED_LAUNCHER:
        fail("root mode is not the fixed installed launcher")
    require_root_owned_nonwritable_chain(installed)
    digest, info = hash_path(
        installed, 4 * 1024 * 1024, "installed launcher authority",
    )
    if (
        info.st_uid != 0 or info.st_gid != 0
        or stat.S_IMODE(info.st_mode) != 0o555 or info.st_nlink != 1
    ):
        fail("installed launcher authority differs")
    return {
        "bytes": info.st_size, "device": info.st_dev, "gid": info.st_gid,
        "inode": info.st_ino, "mode": stat.S_IMODE(info.st_mode),
        "nlink": info.st_nlink, "path": str(installed), "sha256": digest,
        "uid": info.st_uid,
    }


def validate_root_mode_boundary(expected_unit: str,
                                expected_affinity: Iterable[int]) -> Dict[str, Any]:
    if os.getresuid() != (0, 0, 0) or os.getresgid() != (0, 0, 0):
        fail("root mode requires exact root credentials")
    if (
        not sys.flags.isolated or not sys.dont_write_bytecode
        or sys.flags.optimize != 0 or Path(sys.executable) != PYTHON
        or Path(sys.argv[0]).resolve(strict=True) != INSTALLED_LAUNCHER
        or dict(os.environ) != SEALED_ENVIRONMENT
    ):
        fail("root launcher interpreter/environment boundary differs")
    affinity = sorted(os.sched_getaffinity(0))
    if affinity != sorted(set(expected_affinity)):
        fail("root launcher affinity differs")
    groups = sorted(os.getgroups())
    if groups:
        fail("root launcher supplementary groups differ")
    fd = os.open(
        "/proc/self/cgroup", os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK
        | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        raw = os.read(fd, 64 * 1024 + 1)
    finally:
        os.close(fd)
    if len(raw) > 64 * 1024:
        fail("root launcher cgroup receipt exceeds bound")
    relative = parse_proc_cgroup(raw, expected_unit)
    return {
        "affinity": affinity, "cgroup": relative,
        "gid": 0, "groups": groups, "uid": 0,
    }


def build_setpriv_argv(tool_argv: Sequence[str], security_fd: int) -> List[str]:
    if (
        not tool_argv or any(type(item) is not str or not item for item in tool_argv)
        or not Path(tool_argv[0]).is_absolute()
        or type(security_fd) is not int or security_fd < 3
        or security_fd > (1 << 20)
    ):
        fail("build child argv differs")
    return [
        str(SETPRIV), "--reuid", str(CAMPAIGN_UID),
        "--regid", str(CAMPAIGN_GID), "--clear-groups",
        "--bounding-set=-all", "--inh-caps=-all", "--ambient-caps=-all",
        "--no-new-privs", "--pdeathsig", "SIGKILL",
        str(PYTHON), "-I", "-B", "-c", BUILD_CHILD_SECURITY_BOOTSTRAP,
        str(security_fd), *list(tool_argv),
    ]


def read_fd_bounded(fd: int, limit: int, where: str) -> bytes:
    info = os.fstat(fd)
    if not stat.S_ISREG(info.st_mode) or info.st_size < 0 or info.st_size > limit:
        fail(where + " file bound differs")
    data = bytearray()
    offset = 0
    while offset < info.st_size:
        block = os.pread(fd, min(1024 * 1024, info.st_size - offset), offset)
        if not block:
            fail(where + " read was short")
        data.extend(block)
        offset += len(block)
    return bytes(data)


def require_no_build_descendants(where: str) -> None:
    children_path = Path("/proc/self/task") / str(os.getpid()) / "children"
    deadline = time.monotonic() + 2.0
    observed: List[int] = []
    while True:
        try:
            raw = children_path.read_bytes()
        except OSError as exc:
            fail("build descendant roster failed: " + exception_text(exc))
        try:
            pids = [int(token) for token in raw.decode("ascii").split()]
        except (UnicodeDecodeError, ValueError) as exc:
            fail("build descendant roster malformed: " + exception_text(exc))
        for pid in pids:
            if pid <= 1:
                fail("build descendant PID differs")
            observed.append(pid)
            pidfd = -1
            try:
                pidfd = os.pidfd_open(pid, 0)
                signal.pidfd_send_signal(pidfd, signal.SIGKILL)
            except ProcessLookupError:
                pass
            finally:
                if pidfd >= 0:
                    os.close(pidfd)
        while True:
            try:
                pid, status = os.waitpid(-1, os.WNOHANG)
            except ChildProcessError:
                pid = -1
                break
            if pid == 0:
                break
            observed.append(pid)
            wait_status_code(status)
        if not pids and pid == -1:
            break
        if time.monotonic() >= deadline:
            fail("build descendants did not quiesce: " + where)
        time.sleep(0.005)
    if observed:
        fail("build command left descendants: " + where)


def expected_no_work_stdout(tool_argv: Sequence[str]) -> str:
    try:
        index = list(tool_argv).index("-C")
        build = tool_argv[index + 1]
    except (ValueError, IndexError):
        fail("no-work command build path differs")
    return "ninja: Entering directory `" + build + "'\nninja: no work to do.\n"


def validate_build_child_security_receipt(
    value: Mapping[str, Any], expected_pid: int, lower_ns: int, upper_ns: int,
) -> None:
    expected_keys = {
        "affinity", "cap_ambient", "cap_bounding", "cap_effective",
        "cap_inheritable", "cap_permitted", "gid", "groups",
        "no_new_privs", "observed_monotonic_ns", "pid", "receipt_sha256",
        "schema", "uid",
    }
    exact_keys(value, expected_keys, "build child security receipt")
    validate_self_hashed_receipt(
        value, BUILD_CHILD_SECURITY_SCHEMA, "build child security receipt",
    )
    observed = exact_uint(
        value["observed_monotonic_ns"], lower_ns, upper_ns,
        "build child security observation",
    )
    exact_uint(value["pid"], 2, (1 << 31) - 1, "build child security PID")
    exact_uint(value["no_new_privs"], 0, 1, "build child no-new-privs")
    for name in (
        "cap_ambient", "cap_bounding", "cap_effective", "cap_inheritable",
        "cap_permitted",
    ):
        exact_uint(value[name], 0, (1 << 64) - 1, "build child " + name)
    if (
        type(value["affinity"]) is not list
        or any(type(item) is not int for item in value["affinity"])
        or value["pid"] != expected_pid or value["uid"] != [CAMPAIGN_UID] * 3
        or value["gid"] != [CAMPAIGN_GID] * 3 or value["groups"] != []
        or value["affinity"] != [0, 1] or value["no_new_privs"] != 1
        or any(value[name] != 0 for name in (
            "cap_ambient", "cap_bounding", "cap_effective",
            "cap_inheritable", "cap_permitted",
        ))
        or observed > upper_ns
    ):
        fail("build child security policy differs")


def read_build_child_security(fd: int, pid: int, started_ns: int) \
        -> Dict[str, Any]:
    os.set_blocking(fd, False)
    deadline = time.monotonic() + 5.0
    data = bytearray()
    poller = select.poll()
    poller.register(fd, select.POLLIN | select.POLLHUP | select.POLLERR)
    eof = False
    while not eof:
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            fail("build child security receipt deadline expired")
        if not poller.poll(max(1, int(remaining * 1000))):
            fail("build child security receipt deadline expired")
        while True:
            try:
                block = os.read(fd, 16 * 1024 + 1 - len(data))
            except BlockingIOError:
                break
            if not block:
                eof = True
                break
            data.extend(block)
            if len(data) > 16 * 1024:
                fail("build child security receipt exceeds bound")
    value = parse_canonical_line(
        bytes(data), 16 * 1024, "build child security receipt",
    )
    validate_build_child_security_receipt(
        value, pid, started_ns, time.monotonic_ns(),
    )
    return value


def run_build_command(
    name: str, tool_argv: Sequence[str], timeout: float,
    environment: Mapping[str, str], source_root: Path, log_dir_fd: int,
) -> Dict[str, Any]:
    if re.fullmatch(r"[a-z][a-z-]*", name) is None:
        fail("build command name differs")
    if (
        os.getresuid() != (0, 0, 0) or os.getresgid() != (0, 0, 0)
        or os.getgroups()
    ):
        fail("build command requires the exact root launcher boundary")
    descriptors: List[int] = []
    descendants_verified = False
    security_read = -1
    security_write = -1
    process: Optional[subprocess.Popen] = None
    paths = (name + ".stdout", name + ".stderr")
    try:
        for path in paths:
            descriptors.append(os.open(
                path, os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
                | getattr(os, "O_NOFOLLOW", 0), 0o600, dir_fd=log_dir_fd,
            ))
        security_read, security_write = os.pipe2(os.O_CLOEXEC)
        launch_argv = build_setpriv_argv(tool_argv, security_write)
        started_ns = time.monotonic_ns()
        started_utc = utc_now()
        process = subprocess.Popen(
            launch_argv, cwd=str(source_root), env=dict(environment),
            stdin=subprocess.DEVNULL, stdout=descriptors[0], stderr=descriptors[1],
            close_fds=True, pass_fds=(security_write,), start_new_session=True,
        )
        os.close(security_write)
        security_write = -1
        child_security = read_build_child_security(
            security_read, process.pid, started_ns,
        )
        os.close(security_read)
        security_read = -1
        try:
            returncode = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired as exc:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait(timeout=5.0)
            require_no_build_descendants(name)
            descendants_verified = True
            fail("build command deadline expired: " + name)
        require_no_build_descendants(name)
        descendants_verified = True
        finished_ns = time.monotonic_ns()
        finished_utc = utc_now()
        for fd in descriptors:
            os.fsync(fd)
        stdout = read_fd_bounded(
            descriptors[0], MAX_BUILD_COMMAND_OUTPUT_BYTES,
            "build command stdout " + name,
        )
        stderr = read_fd_bounded(
            descriptors[1], MAX_BUILD_COMMAND_OUTPUT_BYTES,
            "build command stderr " + name,
        )
        try:
            stdout_text = stdout.decode("utf-8")
            stderr_text = stderr.decode("utf-8")
        except UnicodeDecodeError as exc:
            fail("build command output is not UTF-8: " + exception_text(exc))
        if returncode != 0:
            fail("build command failed: " + name)
        if name == "no-work" and (
            stdout_text != expected_no_work_stdout(tool_argv) or stderr
        ):
            fail("fresh build did not reach an exact no-work fixed point")
        return self_hashed(
            "wirehair.wh2.direct-systematic-complement-build-command.v1", {
                "environment": dict(environment), "finished_monotonic_ns": finished_ns,
                "finished_utc": finished_utc, "launch_argv": launch_argv,
                "child_pid": process.pid,
                "child_security": child_security,
                "child_security_fd": int(launch_argv[-len(tool_argv) - 1]),
                "name": name, "returncode": returncode,
                "started_monotonic_ns": started_ns, "started_utc": started_utc,
                "stderr_bytes": len(stderr), "stderr_sha256": sha256_bytes(stderr),
                "stderr_utf8": stderr_text, "stdout_bytes": len(stdout),
                "stdout_sha256": sha256_bytes(stdout), "stdout_utf8": stdout_text,
                "tool_argv": list(tool_argv),
            },
        )
    finally:
        if process is not None and process.poll() is None:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            try:
                process.wait(timeout=5.0)
            except (subprocess.TimeoutExpired, OSError):
                pass
        if process is not None and not descendants_verified:
            require_no_build_descendants(name)
        for fd in (security_read, security_write):
            if fd >= 0:
                try:
                    os.close(fd)
                except OSError:
                    pass
        for fd in descriptors:
            try:
                os.close(fd)
            except OSError:
                pass


def seal_campaign_file(path: Path, mode: int, where: str) -> Dict[str, Any]:
    fd = os.open(
        str(path), os.O_RDWR | os.O_CLOEXEC | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        before = os.fstat(fd)
        if (
            not stat.S_ISREG(before.st_mode) or before.st_nlink != 1
            or before.st_uid != CAMPAIGN_UID or before.st_gid != CAMPAIGN_GID
            or before.st_size < 1 or before.st_size > 64 * 1024 * 1024
        ):
            fail(where + " pre-seal policy differs")
        os.fchmod(fd, mode)
        os.fsync(fd)
        sealed = os.fstat(fd)
        digest = hashlib.sha256()
        offset = 0
        while offset < sealed.st_size:
            block = os.pread(fd, min(1024 * 1024, sealed.st_size - offset), offset)
            if not block:
                fail(where + " post-seal read was short")
            digest.update(block)
            offset += len(block)
        after = os.fstat(fd)
        named = os.stat(str(path), follow_symlinks=False)
        if (
            (sealed.st_dev, sealed.st_ino, sealed.st_size, sealed.st_mtime_ns,
             sealed.st_ctime_ns, sealed.st_mode, sealed.st_uid, sealed.st_gid,
             sealed.st_nlink)
            != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns,
                after.st_ctime_ns, after.st_mode, after.st_uid, after.st_gid,
                after.st_nlink)
            or (after.st_dev, after.st_ino) != (named.st_dev, named.st_ino)
            or after.st_uid != CAMPAIGN_UID or after.st_gid != CAMPAIGN_GID
            or after.st_nlink != 1 or stat.S_IMODE(after.st_mode) != mode
        ):
            fail(where + " final seal differs")
        return {
            "bytes": after.st_size, "device": after.st_dev, "gid": after.st_gid,
            "inode": after.st_ino, "mode": stat.S_IMODE(after.st_mode),
            "nlink": after.st_nlink, "path": str(path),
            "sha256": digest.hexdigest(), "uid": after.st_uid,
        }
    finally:
        os.close(fd)


def validate_self_hashed_receipt(value: Mapping[str, Any], schema: str,
                                 where: str) -> str:
    if type(value) is not dict or value.get("schema") != schema:
        fail(where + " schema differs")
    digest = exact_hash(value.get("receipt_sha256"), where)
    preimage = dict(value)
    preimage["receipt_sha256"] = None
    if sha256_bytes(canonical_bytes(preimage)) != digest:
        fail(where + " self-hash differs")
    return digest


def validate_build_command_receipt(
    value: Mapping[str, Any], expected_name: str, expected_argv: Sequence[str],
    expected_environment: Mapping[str, str],
) -> None:
    expected_keys = {
        "child_pid", "child_security", "child_security_fd", "environment",
        "finished_monotonic_ns", "finished_utc", "launch_argv", "name",
        "receipt_sha256", "returncode", "schema",
        "started_monotonic_ns", "started_utc", "stderr_bytes",
        "stderr_sha256", "stderr_utf8", "stdout_bytes", "stdout_sha256",
        "stdout_utf8", "tool_argv",
    }
    exact_keys(value, expected_keys, "build command receipt")
    validate_self_hashed_receipt(
        value, "wirehair.wh2.direct-systematic-complement-build-command.v1",
        "build command receipt",
    )
    security_fd = exact_uint(
        value["child_security_fd"], 3, 1 << 20,
        "build command security descriptor",
    )
    if (
        value["name"] != expected_name or value["tool_argv"] != list(expected_argv)
        or value["launch_argv"] != build_setpriv_argv(expected_argv, security_fd)
        or value["environment"] != dict(expected_environment)
    ):
        fail("build command authority differs: " + expected_name)
    exact_uint(value["returncode"], 0, 0, "build command return code")
    started = exact_uint(
        value["started_monotonic_ns"], 1, (1 << 63) - 1,
        "build command start",
    )
    finished = exact_uint(
        value["finished_monotonic_ns"], started, (1 << 63) - 1,
        "build command finish",
    )
    child_pid = exact_uint(
        value["child_pid"], 2, (1 << 31) - 1, "build child PID",
    )
    security = value["child_security"]
    if type(security) is not dict:
        fail("build child security receipt differs")
    security_pid = exact_uint(
        security.get("pid"), 2, (1 << 31) - 1, "build child PID",
    )
    if security_pid != child_pid:
        fail("build child PID binding differs")
    validate_build_child_security_receipt(
        security, security_pid, started, finished,
    )
    started_utc = exact_utc(value["started_utc"], "build command start")
    finished_utc = exact_utc(value["finished_utc"], "build command finish")
    if finished_utc < started_utc:
        fail("build command UTC chronology differs")
    for prefix in ("stdout", "stderr"):
        text_value = value[prefix + "_utf8"]
        if type(text_value) is not str:
            fail("build command text differs")
        raw = text_value.encode("utf-8")
        observed_bytes = exact_uint(
            value[prefix + "_bytes"], 0, MAX_BUILD_COMMAND_OUTPUT_BYTES,
            "build command " + prefix + " byte count",
        )
        if (
            observed_bytes != len(raw)
            or value[prefix + "_sha256"] != sha256_bytes(raw)
            or len(raw) > MAX_BUILD_COMMAND_OUTPUT_BYTES
        ):
            fail("build command output binding differs")
    if expected_name == "no-work" and (
        value["stdout_utf8"] != expected_no_work_stdout(expected_argv)
        or value["stderr_utf8"] != ""
    ):
        fail("build command no-work evidence differs")


def read_build_authority() -> Tuple[bytes, Dict[str, Any]]:
    parent = directory_binding(
        BUILD_AUTHORITY_DIR, 0, 0, 0o700, "build authority parent",
    )
    del parent
    fd = os.open(
        str(BUILD_AUTHORITY_PATH), os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK
        | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        info = os.fstat(fd)
        named = os.stat(str(BUILD_AUTHORITY_PATH), follow_symlinks=False)
        if (
            not stat.S_ISREG(info.st_mode) or info.st_nlink != 1
            or info.st_uid != 0 or info.st_gid != 0
            or stat.S_IMODE(info.st_mode) != 0o444
            or (info.st_dev, info.st_ino) != (named.st_dev, named.st_ino)
        ):
            fail("build authority file policy differs")
        raw = read_fd_bounded(fd, MAX_BUILD_AUTHORITY_BYTES, "build authority")
        after = os.fstat(fd)
        named_after = os.stat(str(BUILD_AUTHORITY_PATH), follow_symlinks=False)
        if (
            (info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns,
             info.st_ctime_ns, info.st_mode, info.st_uid, info.st_gid, info.st_nlink)
            != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns,
                after.st_ctime_ns, after.st_mode, after.st_uid, after.st_gid,
                after.st_nlink)
            or (after.st_dev, after.st_ino)
            != (named_after.st_dev, named_after.st_ino)
        ):
            fail("build authority changed while reading")
    finally:
        os.close(fd)
    value = parse_canonical_line(raw, MAX_BUILD_AUTHORITY_BYTES, "build authority")
    exact_keys(value, BUILD_AUTHORITY_KEYS, "build authority")
    validate_self_hashed_receipt(value, BUILD_AUTHORITY_SCHEMA, "build authority")
    if value["status"] != "sealed":
        fail("build authority status differs")
    return raw, value


def publish_build_authority(raw: bytes) -> None:
    if (
        not raw or len(raw) > MAX_BUILD_AUTHORITY_BYTES
        or not raw.endswith(b"\n") or raw.count(b"\n") != 1
    ):
        fail("build authority publication framing differs")
    directory_binding(BUILD_AUTHORITY_DIR, 0, 0, 0o700, "build authority parent")
    parent_fd = os.open(
        str(BUILD_AUTHORITY_DIR), os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC
        | getattr(os, "O_NOFOLLOW", 0),
    )
    temporary = "." + BUILD_AUTHORITY_PATH.name + "." + str(os.getpid()) \
        + "." + str(time.monotonic_ns())
    fd = -1
    try:
        fd = os.open(
            temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
            | getattr(os, "O_NOFOLLOW", 0), 0o600, dir_fd=parent_fd,
        )
        view = memoryview(raw)
        while view:
            written = os.write(fd, view)
            if written <= 0:
                fail("build authority write was short")
            view = view[written:]
        os.fchmod(fd, 0o444)
        os.fsync(fd)
        info = os.fstat(fd)
        if (
            info.st_size != len(raw) or info.st_nlink != 1
            or info.st_uid != 0 or info.st_gid != 0
            or stat.S_IMODE(info.st_mode) != 0o444
        ):
            fail("build authority temporary seal differs")
        libc = ctypes.CDLL(None, use_errno=True)
        if not hasattr(libc, "renameat2"):
            fail("renameat2 is unavailable for no-replace authority publication")
        libc.renameat2.argtypes = (
            ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p,
            ctypes.c_uint,
        )
        libc.renameat2.restype = ctypes.c_int
        result = libc.renameat2(
            parent_fd, os.fsencode(temporary), parent_fd,
            os.fsencode(BUILD_AUTHORITY_PATH.name), 1,
        )
        if result != 0:
            error = ctypes.get_errno()
            raise OSError(error, os.strerror(error), str(BUILD_AUTHORITY_PATH))
        os.fsync(parent_fd)
    finally:
        if fd >= 0:
            os.close(fd)
        try:
            os.unlink(temporary, dir_fd=parent_fd)
        except FileNotFoundError:
            pass
        os.close(parent_fd)


def validate_build_staging_paths(config: BuildSealConfig,
                                 build_mode: str) -> Dict[str, Any]:
    if build_mode not in ("absent", "sealed"):
        fail("build staging mode differs")
    source = config.source_dir
    build = config.build_dir
    try:
        source_text = str(source).encode("ascii").decode("ascii")
        build_text = str(build).encode("ascii").decode("ascii")
    except UnicodeError:
        fail("build staging paths are not ASCII")
    if (
        source_text != str(source) or build_text != str(build)
        or re.fullmatch(r"/(?:[A-Za-z0-9._-]+/)*[A-Za-z0-9._-]+", source_text)
        is None
        or re.fullmatch(r"/(?:[A-Za-z0-9._-]+/)*[A-Za-z0-9._-]+", build_text)
        is None
        or source.parent != build.parent
        or source.name in ("", ".", "..") or build.name in ("", ".", "..")
        or re.fullmatch(r"[A-Za-z0-9._-]+", source.name) is None
        or re.fullmatch(r"[A-Za-z0-9._-]+", build.name) is None
    ):
        fail("build staging path relationship differs")
    parent = source.parent
    if parent.resolve(strict=True) != parent:
        fail("build staging parent is not canonical")
    parent_binding = directory_binding(
        parent, 0, 0, 0o755, "build staging parent",
    )
    if source.resolve(strict=True) != source:
        fail("build source root is not canonical")
    source_binding = directory_binding(
        source, CAMPAIGN_UID, CAMPAIGN_GID, 0o555, "build source root",
    )
    if build_mode == "absent":
        try:
            os.stat(str(build), follow_symlinks=False)
        except FileNotFoundError:
            pass
        else:
            fail("fresh build directory is not absent")
        build_binding: Optional[Dict[str, Any]] = None
    else:
        if build.resolve(strict=True) != build:
            fail("sealed build root is not canonical")
        build_binding = directory_binding(
            build, CAMPAIGN_UID, CAMPAIGN_GID, 0o500, "sealed build root",
        )
    return {
        "build": build_binding, "parent": parent_binding,
        "source": source_binding,
    }


def require_scientific_namespaces_absent() -> None:
    fixed = (LAUNCH_DIR, SAMPLER_DIR, FIXED_CLAIM_PATH, FIXED_OUTPUT_DIR)
    for view, root in (("service", None), ("host", HOST_MOUNT_ROOT)):
        for path in fixed:
            observed = path if root is None else root / str(path).lstrip("/")
            try:
                os.stat(str(observed), follow_symlinks=False)
            except FileNotFoundError:
                continue
            fail(
                "scientific namespace already exists during build seal "
                + "(" + view + " view): " + str(path)
            )


def validate_root_boundary_receipt(value: Mapping[str, Any]) -> None:
    exact_keys(value, {"affinity", "cgroup", "gid", "groups", "uid"},
               "build root boundary")
    exact_uint(value["uid"], 0, 0, "build root UID")
    exact_uint(value["gid"], 0, 0, "build root GID")
    if (
        type(value["affinity"]) is not list
        or any(type(item) is not int for item in value["affinity"])
        or value["affinity"] != [0, 1]
        or value["groups"] != [] or type(value["cgroup"]) is not str
        or value["cgroup"].split("/")[-1] != BUILD_UNIT_NAME
    ):
        fail("build root boundary receipt differs")


def current_binary_receipt(build: Path) -> Dict[str, Any]:
    return sealed_file_receipt(
        build / BINARY_NAME, BINARY_NAME, 0o555, 64 * 1024 * 1024,
        "sealed scientific binary",
    )


def validate_build_authority_receipt(
    value: Mapping[str, Any], config: BuildSealConfig,
) -> str:
    exact_keys(value, BUILD_AUTHORITY_KEYS, "build authority")
    receipt_hash = validate_self_hashed_receipt(
        value, BUILD_AUTHORITY_SCHEMA, "build authority",
    )
    if (
        value["status"] != "sealed"
        or value["expected_source_commit"] != config.expected_commit
        or value["source_root"] != str(config.source_dir)
        or value["build_environment"] != build_environment(config.build_dir)
        or value["systemd_run_argv"] != build_systemd_run_argv(config)
        or not canonical_equal(value["trust_boundary"], BUILD_TRUST_BOUNDARY)
    ):
        fail("build authority top-level binding differs")
    validate_root_boundary_receipt(value["root_boundary"])
    staging = validate_build_staging_paths(config, "sealed")
    if not canonical_equal(value["build_directory"], staging["build"]):
        fail("build authority directory binding differs")
    if not canonical_equal(
        value["build_subdirectories"],
        current_build_subdirectories(config.build_dir),
    ):
        fail("build authority subdirectory binding differs")

    current_launcher = installed_launcher_binding()
    if not canonical_equal(value["installed_launcher"], current_launcher):
        fail("build authority installed launcher differs")
    source_launcher = config.source_dir / SOURCE_LAUNCHER_RELATIVE
    source_launcher_hash, source_launcher_info = hash_path(
        source_launcher, 4 * 1024 * 1024,
        "build authority source launcher",
    )
    if (
        source_launcher_hash != current_launcher["sha256"]
        or source_launcher_info.st_uid != CAMPAIGN_UID
        or source_launcher_info.st_gid != CAMPAIGN_GID
        or source_launcher_info.st_nlink != 1
        or stat.S_IMODE(source_launcher_info.st_mode) != 0o444
    ):
        fail("build authority source/installed launcher differs")

    toolchain = frozen_toolchain_receipts()
    if (
        not canonical_equal(value["toolchain_before"], toolchain)
        or not canonical_equal(value["toolchain_after"], toolchain)
    ):
        fail("build authority toolchain binding differs")
    git_hash = static_git_gate(config.source_dir, config.expected_commit)
    if (
        value["git_sha256_before"] != git_hash
        or value["git_sha256_after"] != git_hash
    ):
        fail("build authority Git binding differs")
    source_manifest, _source_hashes = static_source_manifest(config.source_dir)
    if (
        value["source_manifest_sha256_before"] != source_manifest
        or value["source_manifest_sha256_after"] != source_manifest
    ):
        fail("build authority source-manifest binding differs")

    vectors = build_command_vectors(config)
    commands = value["commands"]
    if type(commands) is not list or len(commands) != len(vectors):
        fail("build authority command roster differs")
    previous_finish = 0
    for command, (name, argv, _timeout) in zip(commands, vectors):
        validate_build_command_receipt(
            command, name, argv, build_environment(config.build_dir),
        )
        started = exact_uint(
            command["started_monotonic_ns"], 1, (1 << 63) - 1,
            "build command chronology start",
        )
        finished = exact_uint(
            command["finished_monotonic_ns"], started, (1 << 63) - 1,
            "build command chronology finish",
        )
        if started < previous_finish:
            fail("build command chronology overlaps")
        previous_finish = finished

    build_manifest, build_profile = validate_build_profile(
        config.build_dir, config.source_dir, config.expected_commit,
    )
    if (
        value["build_manifest_sha256"] != build_manifest
        or not canonical_equal(value["build_profile"], build_profile)
        or not canonical_equal(
            value["binary"], current_binary_receipt(config.build_dir)
        )
    ):
        fail("build authority final artifact binding differs")
    return receipt_hash


def _create_build_subdirectory(parent_fd: int, name: str, uid: int, gid: int,
                               mode: int) -> Dict[str, Any]:
    os.mkdir(name, mode, dir_fd=parent_fd)
    os.chown(name, uid, gid, dir_fd=parent_fd, follow_symlinks=False)
    fd = os.open(
        name, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC
        | getattr(os, "O_NOFOLLOW", 0), dir_fd=parent_fd,
    )
    try:
        info = os.fstat(fd)
        named = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
        if (
            not stat.S_ISDIR(info.st_mode)
            or (info.st_dev, info.st_ino) != (named.st_dev, named.st_ino)
            or info.st_uid != uid or info.st_gid != gid
            or stat.S_IMODE(info.st_mode) != mode
        ):
            fail("fresh build subdirectory policy differs: " + name)
        os.fsync(fd)
        return {
            "device": info.st_dev, "gid": info.st_gid,
            "inode": info.st_ino, "mode": stat.S_IMODE(info.st_mode),
            "path": name, "uid": info.st_uid,
        }
    finally:
        os.close(fd)


def seal_build_directory(path: Path, uid: int, gid: int, mode: int,
                         where: str,
                         initial: Optional[Mapping[str, Any]] = None) \
        -> Dict[str, Any]:
    fd = os.open(
        str(path), os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC
        | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        before = os.fstat(fd)
        named = os.stat(str(path), follow_symlinks=False)
        if (
            not stat.S_ISDIR(before.st_mode)
            or (before.st_dev, before.st_ino) != (named.st_dev, named.st_ino)
            or before.st_uid != uid or before.st_gid != gid
        ):
            fail(where + " directory pre-seal differs")
        if initial is not None and (
            initial.get("device") != before.st_dev
            or initial.get("inode") != before.st_ino
            or initial.get("uid") != uid or initial.get("gid") != gid
        ):
            fail(where + " directory identity changed")
        os.fchmod(fd, mode)
        os.fsync(fd)
        after = os.fstat(fd)
        named = os.stat(str(path), follow_symlinks=False)
        if (
            (before.st_dev, before.st_ino) != (after.st_dev, after.st_ino)
            or (after.st_dev, after.st_ino) != (named.st_dev, named.st_ino)
            or after.st_uid != uid or after.st_gid != gid
            or stat.S_IMODE(after.st_mode) != mode
        ):
            fail(where + " directory final seal differs")
        return {
            "device": after.st_dev, "gid": after.st_gid,
            "inode": after.st_ino, "mode": stat.S_IMODE(after.st_mode),
            "path": str(path), "uid": after.st_uid,
        }
    finally:
        os.close(fd)


def current_build_subdirectories(build: Path) -> List[Dict[str, Any]]:
    return [
        directory_binding(
            build / ".wh2-build-home", CAMPAIGN_UID, CAMPAIGN_GID, 0o500,
            "sealed build HOME",
        ),
        directory_binding(
            build / ".wh2-build-tmp", CAMPAIGN_UID, CAMPAIGN_GID, 0o500,
            "sealed build TMPDIR",
        ),
        directory_binding(
            build / ".wh2-build-logs", 0, 0, 0o500,
            "sealed build logs",
        ),
        directory_binding(
            build / "CMakeFiles", CAMPAIGN_UID, CAMPAIGN_GID, 0o555,
            "sealed CMakeFiles parent",
        ),
    ]


def seal_build_authority(config: BuildSealConfig) -> str:
    root_boundary = validate_root_mode_boundary(BUILD_UNIT_NAME, (0, 1))
    require_linux_primitives()
    os.umask(0o077)
    become_subreaper()
    directory_binding(
        BUILD_AUTHORITY_DIR, 0, 0, 0o700, "build authority parent",
    )
    try:
        os.stat(str(BUILD_AUTHORITY_PATH), follow_symlinks=False)
    except FileNotFoundError:
        pass
    else:
        fail("build authority already exists")
    require_scientific_namespaces_absent()
    validate_build_staging_paths(config, "absent")

    launcher_before = installed_launcher_binding()
    source_launcher = config.source_dir / SOURCE_LAUNCHER_RELATIVE
    source_launcher_hash, source_launcher_info = hash_path(
        source_launcher, 4 * 1024 * 1024, "fresh source launcher",
    )
    if (
        source_launcher_hash != launcher_before["sha256"]
        or source_launcher_info.st_uid != CAMPAIGN_UID
        or source_launcher_info.st_gid != CAMPAIGN_GID
        or source_launcher_info.st_nlink != 1
        or stat.S_IMODE(source_launcher_info.st_mode) != 0o444
    ):
        fail("fresh source/installed launcher seal differs")
    git_before = static_git_gate(config.source_dir, config.expected_commit)
    source_before, _source_hashes_before = static_source_manifest(config.source_dir)
    toolchain_before = frozen_toolchain_receipts()

    parent_fd = os.open(
        str(config.build_dir.parent),
        os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC
        | getattr(os, "O_NOFOLLOW", 0),
    )
    build_fd = -1
    log_fd = -1
    try:
        os.mkdir(config.build_dir.name, 0o700, dir_fd=parent_fd)
        os.chown(
            config.build_dir.name, CAMPAIGN_UID, CAMPAIGN_GID,
            dir_fd=parent_fd, follow_symlinks=False,
        )
        build_fd = os.open(
            config.build_dir.name,
            os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC
            | getattr(os, "O_NOFOLLOW", 0), dir_fd=parent_fd,
        )
        initial = os.fstat(build_fd)
        named = os.stat(
            config.build_dir.name, dir_fd=parent_fd, follow_symlinks=False,
        )
        if (
            not stat.S_ISDIR(initial.st_mode)
            or (initial.st_dev, initial.st_ino) != (named.st_dev, named.st_ino)
            or initial.st_uid != CAMPAIGN_UID
            or initial.st_gid != CAMPAIGN_GID
            or stat.S_IMODE(initial.st_mode) != 0o700
        ):
            fail("fresh build root reservation differs")
        home_initial = _create_build_subdirectory(
            build_fd, ".wh2-build-home", CAMPAIGN_UID, CAMPAIGN_GID, 0o700,
        )
        tmp_initial = _create_build_subdirectory(
            build_fd, ".wh2-build-tmp", CAMPAIGN_UID, CAMPAIGN_GID, 0o700,
        )
        logs_initial = _create_build_subdirectory(
            build_fd, ".wh2-build-logs", 0, 0, 0o700,
        )
        log_fd = os.open(
            ".wh2-build-logs",
            os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC
            | getattr(os, "O_NOFOLLOW", 0), dir_fd=build_fd,
        )
        os.fsync(build_fd)
        os.fsync(parent_fd)

        environment = build_environment(config.build_dir)
        commands = [
            run_build_command(
                name, argv, timeout, environment, config.source_dir, log_fd,
            )
            for name, argv, timeout in build_command_vectors(config)
        ]

        git_after = static_git_gate(config.source_dir, config.expected_commit)
        source_after, _source_hashes_after = static_source_manifest(config.source_dir)
        toolchain_after = frozen_toolchain_receipts()
        launcher_after = installed_launcher_binding()
        if (
            git_after != git_before or source_after != source_before
            or toolchain_after != toolchain_before
            or launcher_after != launcher_before
        ):
            fail("fresh build authority changed across command execution")

        for relative in BUILD_PATHS:
            seal_campaign_file(
                config.build_dir / relative, 0o444,
                "fresh build provenance " + relative,
            )
        seal_campaign_file(
            config.build_dir / BINARY_NAME, 0o555,
            "fresh scientific binary",
        )
        build_subdirectories = [
            seal_build_directory(
                config.build_dir / ".wh2-build-home", CAMPAIGN_UID,
                CAMPAIGN_GID, 0o500, "fresh build HOME", home_initial,
            ),
            seal_build_directory(
                config.build_dir / ".wh2-build-tmp", CAMPAIGN_UID,
                CAMPAIGN_GID, 0o500, "fresh build TMPDIR", tmp_initial,
            ),
            seal_build_directory(
                config.build_dir / ".wh2-build-logs", 0, 0, 0o500,
                "fresh build logs", logs_initial,
            ),
            seal_build_directory(
                config.build_dir / "CMakeFiles", CAMPAIGN_UID,
                CAMPAIGN_GID, 0o555, "fresh CMakeFiles parent",
            ),
        ]
        os.fchmod(build_fd, 0o500)
        os.fsync(build_fd)
        os.fsync(parent_fd)
        final_info = os.fstat(build_fd)
        final_named = os.stat(
            config.build_dir.name, dir_fd=parent_fd, follow_symlinks=False,
        )
        if (
            (initial.st_dev, initial.st_ino)
            != (final_info.st_dev, final_info.st_ino)
            or (final_info.st_dev, final_info.st_ino)
            != (final_named.st_dev, final_named.st_ino)
            or final_info.st_uid != CAMPAIGN_UID
            or final_info.st_gid != CAMPAIGN_GID
            or stat.S_IMODE(final_info.st_mode) != 0o500
        ):
            fail("fresh build root final seal differs")

        build_manifest, build_profile = validate_build_profile(
            config.build_dir, config.source_dir, config.expected_commit,
        )
        staging = validate_build_staging_paths(config, "sealed")
        receipt = self_hashed(BUILD_AUTHORITY_SCHEMA, {
            "binary": current_binary_receipt(config.build_dir),
            "build_directory": staging["build"],
            "build_environment": environment,
            "build_manifest_sha256": build_manifest,
            "build_profile": build_profile,
            "build_subdirectories": build_subdirectories,
            "commands": commands,
            "expected_source_commit": config.expected_commit,
            "git_sha256_after": git_after, "git_sha256_before": git_before,
            "installed_launcher": launcher_after,
            "root_boundary": root_boundary,
            "source_manifest_sha256_after": source_after,
            "source_manifest_sha256_before": source_before,
            "source_root": str(config.source_dir), "status": "sealed",
            "systemd_run_argv": build_systemd_run_argv(config),
            "toolchain_after": toolchain_after,
            "toolchain_before": toolchain_before,
            "trust_boundary": BUILD_TRUST_BOUNDARY,
        })
        receipt_hash = validate_build_authority_receipt(receipt, config)
        raw = canonical_bytes(receipt) + b"\n"
        require_scientific_namespaces_absent()
        publish_build_authority(raw)
        published_raw, published = read_build_authority()
        if published_raw != raw:
            fail("published build authority bytes differ")
        if validate_build_authority_receipt(published, config) != receipt_hash:
            fail("published build authority validation differs")
        return receipt_hash
    finally:
        if log_fd >= 0:
            os.close(log_fd)
        if build_fd >= 0:
            os.close(build_fd)
        os.close(parent_fd)


class _StatFs(ctypes.Structure):
    _fields_ = [
        ("f_type", ctypes.c_long), ("f_bsize", ctypes.c_long),
        ("f_blocks", ctypes.c_ulong * 7), ("f_fsid", ctypes.c_int * 2),
        ("f_namelen", ctypes.c_long), ("f_frsize", ctypes.c_long),
        ("f_flags", ctypes.c_long), ("f_spare", ctypes.c_long * 4),
    ]


def require_cgroup2(path: Path) -> None:
    libc = ctypes.CDLL(None, use_errno=True)
    result = _StatFs()
    if libc.statfs(os.fsencode(str(path)), ctypes.byref(result)) != 0:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error), str(path))
    if result.f_type != CGROUP2_SUPER_MAGIC:
        fail("/sys/fs/cgroup is not a cgroup-v2 filesystem")


def become_subreaper() -> None:
    libc = ctypes.CDLL(None, use_errno=True)
    if libc.prctl(PR_SET_CHILD_SUBREAPER, 1, 0, 0, 0) != 0:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error))
    value = ctypes.c_int(0)
    if libc.prctl(PR_GET_CHILD_SUBREAPER, ctypes.byref(value), 0, 0, 0) != 0:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error))
    if value.value != 1:
        fail("root supervisor did not become a child subreaper")


def require_linux_primitives() -> None:
    required = {
        "os.fork": hasattr(os, "fork"),
        "os.memfd_create": hasattr(os, "memfd_create"),
        "os.pidfd_open": hasattr(os, "pidfd_open"),
        "os.sched_getaffinity": hasattr(os, "sched_getaffinity"),
        "os.sched_setaffinity": hasattr(os, "sched_setaffinity"),
        "signal.pidfd_send_signal": hasattr(signal, "pidfd_send_signal"),
        "socket.SCM_RIGHTS": hasattr(socket, "SCM_RIGHTS"),
        "socket.SOCK_SEQPACKET": hasattr(socket, "SOCK_SEQPACKET"),
    }
    libc = ctypes.CDLL(None)
    required["timerfd_create"] = hasattr(libc, "timerfd_create")
    required["timerfd_settime"] = hasattr(libc, "timerfd_settime")
    required["renameat2"] = hasattr(libc, "renameat2")
    missing = sorted(name for name, available in required.items() if not available)
    if missing:
        fail("required Linux/Python primitives unavailable: " + ",".join(missing))


class CgroupTree:
    def __init__(self, service: Path) -> None:
        self.service = service

    @staticmethod
    def _write(path: Path, value: str) -> None:
        fd = os.open(str(path), os.O_WRONLY | os.O_CLOEXEC | getattr(os, "O_NOFOLLOW", 0))
        try:
            raw = value.encode("ascii")
            if os.write(fd, raw) != len(raw):
                fail("cgroup write was short: " + str(path))
        finally:
            os.close(fd)

    @staticmethod
    def _read(path: Path, limit: int = 65536) -> str:
        fd = os.open(str(path), os.O_RDONLY | os.O_CLOEXEC | getattr(os, "O_NOFOLLOW", 0))
        try:
            raw = os.read(fd, limit + 1)
        finally:
            os.close(fd)
        if len(raw) > limit:
            fail("cgroup control file exceeds bound: " + str(path))
        try:
            return raw.decode("ascii").strip()
        except UnicodeDecodeError:
            fail("cgroup control file is not ASCII: " + str(path))

    def _mkdir(self, parent: Path, name: str) -> Path:
        path = parent / name
        os.mkdir(str(path), 0o755)
        info = os.stat(str(path), follow_symlinks=False)
        if not stat.S_ISDIR(info.st_mode) or info.st_uid != 0 or info.st_gid != 0:
            fail("cgroup leaf policy differs: " + str(path))
        return path

    def _configure(self, path: Path, cpus: str, pids: int) -> None:
        self._write(path / "cpuset.mems", MEMORY_NODES)
        self._write(path / "cpuset.cpus", cpus)
        self._write(path / "pids.max", str(pids))
        if (
            self._read(path / "cpuset.mems") != MEMORY_NODES
            or self._read(path / "cpuset.mems.effective") != MEMORY_NODES
            or self._read(path / "cpuset.cpus") != cpus
            or self._read(path / "cpuset.cpus.effective") != cpus
            or self._read(path / "pids.max") != str(pids)
            or self._read(path / "cgroup.type") != "domain"
        ):
            fail("cgroup cpuset/memory/pids readback differs")

    def create(self) -> CgroupLayout:
        controllers = set(self._read(self.service / "cgroup.controllers").split())
        if not {"cpuset", "pids"}.issubset(controllers):
            fail("delegated service lacks cpuset/pids controllers")
        supervisor = self._mkdir(self.service, "supervisor")
        self._write(supervisor / "cgroup.procs", str(os.getpid()))
        os.sched_setaffinity(0, {VERIFIER_CPU})
        if os.sched_getaffinity(0) != {VERIFIER_CPU}:
            fail("root supervisor did not pin to CPU123")
        if self._read(self.service / "cgroup.procs"):
            fail("delegated service root is not empty")
        self._write(self.service / "cgroup.subtree_control", "+cpuset +pids")
        if set(self._read(self.service / "cgroup.subtree_control").split()) != {
            "cpuset", "pids",
        }:
            fail("delegated service subtree-controller readback differs")
        self._configure(supervisor, str(VERIFIER_CPU), 4)
        run = self._mkdir(self.service, "run")
        self._configure(run, "120-123", 28)
        self._write(run / "cgroup.subtree_control", "+cpuset +pids")
        if set(self._read(run / "cgroup.subtree_control").split()) != {
            "cpuset", "pids",
        }:
            fail("run subtree-controller readback differs")
        sampler = self._mkdir(run, "sampler")
        experiment = self._mkdir(run, "experiment")
        verifier = self._mkdir(run, "verifier")
        self._configure(sampler, str(SAMPLER_CPU), 3)
        self._configure(experiment, "120-122", 12)
        self._configure(verifier, str(VERIFIER_CPU), 3)
        for leaf in (sampler, experiment, verifier):
            events = self._read(leaf / "cgroup.events").splitlines()
            if "populated 0" not in events:
                fail("new cgroup leaf is unexpectedly populated")
        kill_fd = os.open(
            str(experiment / "cgroup.kill"), os.O_WRONLY | os.O_CLOEXEC
            | getattr(os, "O_NOFOLLOW", 0),
        )
        sampler_kill_fd = os.open(
            str(sampler / "cgroup.kill"), os.O_WRONLY | os.O_CLOEXEC
            | getattr(os, "O_NOFOLLOW", 0),
        )
        verifier_kill_fd = os.open(
            str(verifier / "cgroup.kill"), os.O_WRONLY | os.O_CLOEXEC
            | getattr(os, "O_NOFOLLOW", 0),
        )
        run_kill_fd = os.open(
            str(run / "cgroup.kill"), os.O_WRONLY | os.O_CLOEXEC
            | getattr(os, "O_NOFOLLOW", 0),
        )
        return CgroupLayout(
            service=self.service, supervisor=supervisor, run=run,
            sampler=sampler, experiment=experiment, verifier=verifier,
            experiment_kill_fd=kill_fd, sampler_kill_fd=sampler_kill_fd,
            verifier_kill_fd=verifier_kill_fd, run_kill_fd=run_kill_fd,
        )

    @staticmethod
    def move_pid(path: Path, pid: int, cpus: Iterable[int]) -> None:
        CgroupTree._write(path / "cgroup.procs", str(pid))
        os.sched_setaffinity(pid, set(cpus))
        expected = sorted(set(cpus))
        if sorted(os.sched_getaffinity(pid)) != expected:
            fail("child affinity differs after cgroup placement")
        actual = CgroupTree._read(Path(f"/proc/{pid}/cgroup"))
        if not actual.endswith("/" + path.name):
            fail("child cgroup placement differs")

    @staticmethod
    def pids(path: Path) -> List[int]:
        raw = CgroupTree._read(path / "cgroup.procs")
        if not raw:
            return []
        result: List[int] = []
        for token in raw.splitlines():
            if CANONICAL_UINT.fullmatch(token) is None or int(token) <= 0:
                fail("cgroup process roster is malformed")
            result.append(int(token))
        return sorted(set(result))


def make_sealed_memfd(name: str, data: bytes) -> int:
    if not hasattr(os, "memfd_create"):
        fail("sealed role transfer requires memfd_create")
    flags = os.MFD_CLOEXEC | getattr(os, "MFD_ALLOW_SEALING", 0x0002)
    fd = os.memfd_create(name, flags)
    try:
        view = memoryview(data)
        while view:
            written = os.write(fd, view)
            if written <= 0:
                fail("memfd write was short")
            view = view[written:]
        seals = (
            fcntl.F_SEAL_SEAL | fcntl.F_SEAL_SHRINK
            | fcntl.F_SEAL_GROW | fcntl.F_SEAL_WRITE
        )
        fcntl.fcntl(fd, fcntl.F_ADD_SEALS, seals)
        if fcntl.fcntl(fd, fcntl.F_GET_SEALS) != seals:
            fail("memfd seal roster differs")
        os.lseek(fd, 0, os.SEEK_SET)
        return fd
    except BaseException:
        os.close(fd)
        raise


def read_sealed_memfd(fd: int, limit: int, where: str) -> bytes:
    seals = (
        fcntl.F_SEAL_SEAL | fcntl.F_SEAL_SHRINK
        | fcntl.F_SEAL_GROW | fcntl.F_SEAL_WRITE
    )
    if fcntl.fcntl(fd, fcntl.F_GET_SEALS) != seals:
        fail(where + " memfd seals differ")
    info = os.fstat(fd)
    if not stat.S_ISREG(info.st_mode) or info.st_size < 1 or info.st_size > limit:
        fail(where + " memfd size/type differs")
    data = bytearray()
    offset = 0
    while offset < info.st_size:
        block = os.pread(fd, min(1024 * 1024, info.st_size - offset), offset)
        if not block:
            fail(where + " memfd read was short")
        data.extend(block)
        offset += len(block)
    return bytes(data)


def send_fds(sock: socket.socket, fds: Sequence[int]) -> None:
    payload = array.array("i", fds)
    sent = sock.sendmsg([b"A"], [(socket.SOL_SOCKET, socket.SCM_RIGHTS, payload)])
    if sent != 1:
        fail("role descriptor transfer was short")


def receive_fds(sock: socket.socket) -> Tuple[int, int]:
    space = socket.CMSG_SPACE(2 * array.array("i").itemsize)
    data, ancdata, flags, _ = sock.recvmsg(1, space)
    if data != b"A" or flags & (socket.MSG_TRUNC | socket.MSG_CTRUNC):
        fail("role descriptor transfer framing differs")
    received: List[int] = []
    for level, kind, payload in ancdata:
        if level != socket.SOL_SOCKET or kind != socket.SCM_RIGHTS:
            fail("role descriptor transfer ancillary type differs")
        values = array.array("i")
        usable = len(payload) - (len(payload) % values.itemsize)
        values.frombytes(payload[:usable])
        received.extend(values.tolist())
    if len(received) != 2:
        for fd in received:
            os.close(fd)
        fail("role descriptor transfer count differs")
    return received[0], received[1]


def child_exec_failure(message: str) -> None:
    try:
        os.write(2, (message[:4000] + "\n").encode("ascii", "backslashreplace"))
    finally:
        os._exit(125)


def wait_status_code(status: int) -> int:
    if os.WIFEXITED(status):
        return os.WEXITSTATUS(status)
    if os.WIFSIGNALED(status):
        return -os.WTERMSIG(status)
    fail("child wait status differs")


def pidfd_signal(pidfd: int, signum: int) -> None:
    if not hasattr(signal, "pidfd_send_signal") or pidfd < 0:
        fail("pidfd signaling is unavailable")
    signal.pidfd_send_signal(pidfd, signum)


def reap_failed_spawn(pid: int, pidfd: int, seconds: float = 2.0) -> None:
    """Kill and reap a fork that failed before its startup release.

    A pidfd is authoritative whenever it was obtained.  If pidfd_open itself
    failed, the PID still names our live, unreaped direct child and therefore
    cannot have been reused; numeric signaling is confined to this pre-release
    exception path.  Normal launched-process cleanup is pidfd-only.
    """
    failures: List[str] = []
    if pidfd >= 0:
        try:
            pidfd_signal(pidfd, signal.SIGKILL)
        except ProcessLookupError:
            pass
        except (OSError, LaunchError) as exc:
            failures.append("pidfd signal: " + exception_text(exc))
    else:
        try:
            waited, status = os.waitpid(pid, os.WNOHANG)
            if waited == pid:
                wait_status_code(status)
                return
            os.kill(pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        except (OSError, LaunchError) as exc:
            failures.append("direct-child signal: " + exception_text(exc))
    deadline = time.monotonic() + seconds
    reaped = False
    while time.monotonic() < deadline:
        try:
            waited, status = os.waitpid(pid, os.WNOHANG)
        except ChildProcessError:
            reaped = True
            break
        if waited == pid:
            wait_status_code(status)
            reaped = True
            break
        time.sleep(0.005)
    if pidfd >= 0:
        try:
            os.close(pidfd)
        except OSError as exc:
            failures.append("pidfd close: " + exception_text(exc))
    if not reaped:
        failures.append("child did not exit through its closed startup channel")
    if failures:
        fail("failed-spawn cleanup differs: " + " | ".join(failures))


class _Timespec(ctypes.Structure):
    _fields_ = [("tv_sec", ctypes.c_long), ("tv_nsec", ctypes.c_long)]


class _Itimerspec(ctypes.Structure):
    _fields_ = [("it_interval", _Timespec), ("it_value", _Timespec)]


def timerfd_create() -> int:
    libc = ctypes.CDLL(None, use_errno=True)
    fd = libc.timerfd_create(1, os.O_CLOEXEC | os.O_NONBLOCK)
    if fd < 0:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error))
    return fd


def timerfd_arm(fd: int, deadline_ns: Optional[int]) -> None:
    libc = ctypes.CDLL(None, use_errno=True)
    when = 0 if deadline_ns is None else deadline_ns
    value = _Itimerspec(
        _Timespec(0, 0),
        _Timespec(when // 1_000_000_000, when % 1_000_000_000),
    )
    if libc.timerfd_settime(fd, 1, ctypes.byref(value), None) != 0:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error))


def close_fds_except(keep: Iterable[int]) -> None:
    retained = set(keep)
    for name in os.listdir("/proc/self/fd"):
        try:
            fd = int(name)
        except ValueError:
            continue
        if fd > 2 and fd not in retained:
            try:
                os.close(fd)
            except OSError:
                pass


def cgroup_descendant_pids(path: Path) -> List[int]:
    result: List[int] = []
    for root, directories, _files in os.walk(str(path), topdown=True, followlinks=False):
        directories[:] = sorted(
            name for name in directories
            if name not in (".", "..") and "/" not in name
        )
        result.extend(CgroupTree.pids(Path(root)))
    return sorted(set(result))


def guardian_event(
    directory_fd: int, name: str, event: str, deadline_ns: int,
    observed_ns: int, roster: Sequence[int], errors: Sequence[str],
) -> None:
    value = self_hashed(
        "wirehair.wh2.direct-systematic-complement-launch-deadline-event.v1",
        {
            "deadline_monotonic_ns": deadline_ns,
            "errors": list(errors)[:16], "event": event,
            "guardian_pid": os.getpid(),
            "observed_monotonic_ns": observed_ns,
            "process_roster": list(roster),
        },
    )
    raw = canonical_bytes(value) + b"\n"
    fd = os.open(
        name, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
        | getattr(os, "O_NOFOLLOW", 0), 0o600, dir_fd=directory_fd,
    )
    try:
        view = memoryview(raw)
        while view:
            written = os.write(fd, view)
            if written <= 0:
                fail("guardian event write was short")
            view = view[written:]
        os.fchmod(fd, 0o400)
        os.fsync(fd)
        os.fsync(directory_fd)
    finally:
        os.close(fd)


def guardian_response(sequence: int, state: str, outer_fired: bool,
                      backstop_fired: bool, error: str) -> bytes:
    return canonical_bytes({
        "backstop_fired": backstop_fired, "error": error,
        "outer_fired": outer_fired, "sequence": sequence,
        "state": state,
    }) + b"\n"


def guardian_child(
    sock: socket.socket, experiment_kill_fd: int, run_kill_fd: int,
    experiment: Path, run: Path, journal_fd: int,
) -> None:
    """Independent admission/outer/backstop authority.

    It remains runnable if the Python supervisor is SIGSTOPed.  A closed parent
    socket never disarms an established deadline.
    """
    outer_timer = timerfd_create()
    backstop_timer = timerfd_create()
    keep = {
        sock.fileno(), experiment_kill_fd, run_kill_fd, journal_fd,
        outer_timer, backstop_timer,
    }
    close_fds_except(keep)
    poller = select.poll()
    poller.register(sock.fileno(), select.POLLIN | select.POLLHUP | select.POLLERR)
    poller.register(outer_timer, select.POLLIN)
    poller.register(backstop_timer, select.POLLIN)
    state = "ready"
    sequence = 0
    release_not_before_ns: Optional[int] = None
    outer_deadline_ns: Optional[int] = None
    backstop_deadline_ns: Optional[int] = None
    outer_fired = False
    backstop_fired = False
    error = "none"
    if sock.send(guardian_response(sequence, state, False, False, error)) <= 0:
        fail("guardian initial ACK was short")
    parent_connected = True
    while True:
        ready = {
            fd: events for fd, events in poller.poll()
            if events & (select.POLLIN | select.POLLHUP | select.POLLERR)
        }
        # Deadline readiness is deliberately handled before cancellation.
        if outer_timer in ready:
            try:
                os.read(outer_timer, 8)
            except BlockingIOError:
                pass
            observed = time.monotonic_ns()
            roster = CgroupTree.pids(experiment)
            event = "admission" if state == "admission" else "outer"
            event_errors: List[str] = []
            try:
                guardian_event(
                    journal_fd, "deadline-fired.json", event,
                    outer_deadline_ns or observed, observed, roster, event_errors,
                )
            except BaseException as exc:
                event_errors.append("durable event: " + exception_text(exc))
                error = event_errors[-1]
            try:
                if os.write(experiment_kill_fd, b"1") != 1:
                    fail("experiment cgroup.kill write was short")
            except BaseException as exc:
                # Primary containment failure must not destroy the independent
                # +140 authority.  Record it and remain alive for run kill.
                error = "experiment kill: " + exception_text(exc)
            outer_fired = True
            outer_deadline_ns = None
            state = "outer_fired"
        if backstop_timer in ready:
            try:
                os.read(backstop_timer, 8)
            except BlockingIOError:
                pass
            observed = time.monotonic_ns()
            roster = cgroup_descendant_pids(run)
            event_errors = []
            try:
                guardian_event(
                    journal_fd, "backstop-fired.json", "backstop",
                    backstop_deadline_ns or observed, observed, roster,
                    event_errors,
                )
            except BaseException as exc:
                event_errors.append("durable event: " + exception_text(exc))
                error = event_errors[-1]
            try:
                if os.write(run_kill_fd, b"1") != 1:
                    fail("run cgroup.kill write was short")
            except BaseException as exc:
                error = "run kill: " + exception_text(exc)
            backstop_fired = True
            os._exit(0 if error == "none" else 125)
        socket_events = ready.get(sock.fileno(), 0)
        if socket_events & select.POLLIN:
            raw = sock.recv(4097)
            if not raw:
                parent_connected = False
                poller.unregister(sock.fileno())
                continue
            command = parse_canonical_line(raw, 4096, "guardian command")
            exact_keys(command, {"command", "sequence", "values"}, "guardian command")
            sequence = exact_uint(command["sequence"], sequence + 1,
                                  sequence + 1, "guardian sequence")
            kind = command["command"]
            values = command["values"]
            if type(values) is not dict:
                fail("guardian command values differ")
            if kind == "admission":
                if state != "ready":
                    fail("guardian admission state differs")
                exact_keys(values, {"release_not_before_monotonic_ns"},
                           "guardian admission values")
                release_not_before_ns = exact_uint(
                    values["release_not_before_monotonic_ns"], 1,
                    (1 << 63) - 1, "guardian release epoch",
                )
                outer_deadline_ns = release_not_before_ns + int(
                    CLAIM_ADMISSION_SECONDS * 1e9
                )
                backstop_deadline_ns = release_not_before_ns + int(
                    (CLAIM_ADMISSION_SECONDS + WHOLE_RUN_SECONDS) * 1e9
                )
                timerfd_arm(outer_timer, outer_deadline_ns)
                timerfd_arm(backstop_timer, backstop_deadline_ns)
                state = "admission"
            elif kind == "exact":
                if state not in ("admission", "outer_fired") \
                        or release_not_before_ns is None:
                    fail("guardian exact-arm state differs")
                exact_keys(values, {"controller_started_monotonic_ns"},
                           "guardian exact values")
                started = exact_uint(
                    values["controller_started_monotonic_ns"],
                    release_not_before_ns, (1 << 63) - 1,
                    "guardian controller epoch",
                )
                exact_outer_deadline_ns = started + int(OUTER_SECONDS * 1e9)
                backstop_deadline_ns = started + int(WHOLE_RUN_SECONDS * 1e9)
                timerfd_arm(backstop_timer, backstop_deadline_ns)
                if state == "admission":
                    outer_deadline_ns = exact_outer_deadline_ns
                    timerfd_arm(outer_timer, outer_deadline_ns)
                    state = "exact"
                # At the inclusive admission edge the experiment was already
                # killed above.  Preserve that fired state while still
                # re-anchoring the independent +140 backstop to the authentic
                # controller epoch; never crash the guardian in this race.
            elif kind == "cancel_outer":
                if state not in ("admission", "exact", "outer_fired"):
                    fail("guardian cancel state differs")
                if not outer_fired:
                    timerfd_arm(outer_timer, None)
                state = "outer_cancelled"
            elif kind == "finish":
                if state not in ("outer_cancelled", "outer_fired"):
                    fail("guardian finish state differs")
                timerfd_arm(outer_timer, None)
                timerfd_arm(backstop_timer, None)
                state = "finished"
            else:
                fail("guardian command differs")
            response = guardian_response(
                sequence, state, outer_fired, backstop_fired, error
            )
            try:
                sent = sock.send(response)
            except OSError:
                sent = 0
            if sent != len(response):
                # Once a deadline is armed, parent loss must not disarm it.
                parent_connected = False
                try:
                    poller.unregister(sock.fileno())
                except KeyError:
                    pass
            if state == "finished":
                os._exit(0 if error == "none" else 125)
        elif socket_events & (select.POLLHUP | select.POLLERR):
            parent_connected = False
            poller.unregister(sock.fileno())
        if not parent_connected and backstop_deadline_ns is None:
            os._exit(125)


class GuardianHandle:
    def __init__(self, process: ProcessHandle, sock: socket.socket) -> None:
        self.process = process
        self.sock = sock
        self.sequence = 0
        self.fired = False
        self.backstop_fired = False
        self.roster: List[int] = []
        self.error = "none"
        self.finished = False

    def _response(self) -> Dict[str, Any]:
        raw = self.sock.recv(4097)
        value = parse_canonical_line(raw, 4096, "guardian ACK")
        exact_keys(value, {
            "backstop_fired", "error", "outer_fired", "sequence", "state",
        }, "guardian ACK")
        if value["sequence"] != self.sequence:
            fail("guardian ACK sequence differs")
        if type(value["outer_fired"]) is not bool or type(value["backstop_fired"]) is not bool:
            fail("guardian ACK booleans differ")
        if type(value["error"]) is not str or not value["error"]:
            fail("guardian ACK error differs")
        self.fired = value["outer_fired"]
        self.backstop_fired = value["backstop_fired"]
        self.error = value["error"]
        return value

    def initial_ack(self) -> None:
        self.sock.settimeout(5.0)
        value = self._response()
        if value["state"] != "ready":
            fail("guardian initial state differs")
        self.sock.settimeout(None)

    def _request(self, command: str, values: Mapping[str, Any],
                 expected_states: Tuple[str, ...]) -> str:
        self.sequence += 1
        raw = canonical_bytes({
            "command": command, "sequence": self.sequence,
            "values": dict(values),
        }) + b"\n"
        self.sock.settimeout(5.0)
        if self.sock.send(raw) != len(raw):
            fail("guardian command write was short")
        response = self._response()
        self.sock.settimeout(None)
        if response["state"] not in expected_states:
            fail("guardian command state differs")
        return response["state"]

    def arm_admission(self, release_not_before_ns: int) -> None:
        self._request("admission", {
            "release_not_before_monotonic_ns": release_not_before_ns,
        }, ("admission",))

    def arm_exact(self, controller_started_ns: int) -> bool:
        state = self._request("exact", {
            "controller_started_monotonic_ns": controller_started_ns,
        }, ("exact", "outer_fired"))
        return state == "exact"

    def cancel(self) -> None:
        self._request("cancel_outer", {}, ("outer_cancelled",))

    def finish(self) -> None:
        self._request("finish", {}, ("finished",))
        poller = select.poll()
        poller.register(self.process.pidfd, select.POLLIN)
        if not poller.poll(2000):
            fail("guardian process did not exit after finish")
        waited, status = os.waitpid(self.process.pid, 0)
        if waited != self.process.pid:
            fail("guardian reap identity differs")
        self.process.returncode = wait_status_code(status)
        if self.process.returncode != 0:
            fail("guardian process exit differs")
        self.finished = True

    def close(self) -> bool:
        self.sock.close()
        return True


class RealBackend:
    """Linux/systemd implementation behind the testable supervisor protocol."""

    def __init__(self) -> None:
        self.lock_fd = -1
        self.layout: Optional[CgroupLayout] = None
        self.journal: Optional[AttemptJournal] = None
        self.attempt_consumed = False
        self.held_sampler: Dict[str, int] = {}
        self.held_sources: Dict[str, int] = {}
        self.children: List[ProcessHandle] = []
        self.sampler_dir_fd = -1
        self.sampler_dir_binding: Optional[Tuple[int, int]] = None
        self.sampler_dir_state = "not_reserved"

    @staticmethod
    def _read_bounded(path: Path, limit: int) -> bytes:
        fd = os.open(
            str(path), os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK
            | getattr(os, "O_NOFOLLOW", 0),
        )
        try:
            info = os.fstat(fd)
            if not stat.S_ISREG(info.st_mode) or info.st_size < 0 or info.st_size > limit:
                fail("bounded input policy differs: " + str(path))
            data = bytearray()
            while len(data) < info.st_size:
                block = os.read(fd, min(1024 * 1024, info.st_size - len(data)))
                if not block:
                    fail("bounded input read was short: " + str(path))
                data.extend(block)
            return bytes(data)
        finally:
            os.close(fd)

    @staticmethod
    def _read_proc_bounded(path: Path, limit: int) -> bytes:
        if type(limit) is not int or limit < 1:
            fail("procfs bounded input limit differs")
        fd = os.open(
            str(path), os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK
            | getattr(os, "O_NOFOLLOW", 0),
        )
        try:
            info = os.fstat(fd)
            if not stat.S_ISREG(info.st_mode):
                fail("procfs bounded input policy differs: " + str(path))
            data = bytearray()
            while len(data) <= limit:
                block = os.read(fd, min(16384, limit + 1 - len(data)))
                if not block:
                    break
                data.extend(block)
            if len(data) > limit:
                fail("procfs bounded input exceeds limit: " + str(path))
            return bytes(data)
        finally:
            os.close(fd)

    @staticmethod
    def _read_claim_if_ready() -> Optional[bytes]:
        fd = os.open(
            str(FIXED_CLAIM_PATH), os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK
            | getattr(os, "O_NOFOLLOW", 0),
        )
        try:
            before = os.fstat(fd)
            mode = stat.S_IMODE(before.st_mode)
            if (
                not stat.S_ISREG(before.st_mode) or before.st_nlink != 1
                or before.st_uid != CAMPAIGN_UID or before.st_gid != CAMPAIGN_GID
                or before.st_size < 0 or before.st_size > MAX_CLAIM_BYTES
            ):
                fail("campaign claim file policy differs")
            # O_EXCL publication is visible before the controller completes its
            # write+fchmod+fsync sequence.  Exact 0600 is the sole transient
            # state; only 0400 licenses parsing as a published claim.
            if mode == 0o600:
                return None
            if mode != 0o400 or before.st_size < 1:
                fail("campaign claim publication mode differs")
            data = bytearray()
            offset = 0
            while offset < before.st_size:
                block = os.pread(
                    fd, min(1024 * 1024, before.st_size - offset), offset
                )
                if not block:
                    fail("campaign claim read was short")
                data.extend(block)
                offset += len(block)
            after = os.fstat(fd)
            named = os.stat(str(FIXED_CLAIM_PATH), follow_symlinks=False)
            if (
                (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns,
                 before.st_mode, before.st_uid, before.st_gid, before.st_nlink)
                != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns,
                    after.st_mode, after.st_uid, after.st_gid, after.st_nlink)
                or (after.st_dev, after.st_ino) != (named.st_dev, named.st_ino)
            ):
                fail("campaign claim changed while reading")
            return bytes(data)
        finally:
            os.close(fd)

    def prepare(self, config: LaunchConfig) -> Prepared:
        if os.getresuid() != (0, 0, 0) or os.getresgid() != (0, 0, 0):
            fail("execute mode requires exact root credentials")
        if (
            not sys.flags.isolated or not sys.dont_write_bytecode
            or sys.flags.optimize != 0
            or Path(sys.executable) != PYTHON
            or Path(sys.argv[0]).resolve(strict=True) != INSTALLED_LAUNCHER
            or dict(os.environ) != SEALED_ENVIRONMENT
        ):
            fail("root launcher interpreter/environment boundary differs")
        require_linux_primitives()
        os.umask(0o077)
        become_subreaper()
        build = config.build_dir.resolve(strict=True)
        if build != config.build_dir or not build.is_dir():
            fail("build directory is not canonical")
        build_authority_raw, build_authority = read_build_authority()
        authority_source = build_authority.get("source_root")
        if type(authority_source) is not str:
            fail("build authority source root differs")
        authority_config = BuildSealConfig(
            source_dir=Path(authority_source), build_dir=build,
            expected_commit=config.expected_commit,
        )
        validate_build_authority_receipt(build_authority, authority_config)
        cache_path = build / "CMakeCache.txt"
        cache = self._read_bounded(cache_path, 64 * 1024 * 1024)
        source_root = parse_cmake_source_root(cache, build)
        if source_root != authority_config.source_dir:
            fail("build authority/CMake source root differs")
        controller = source_root / CONTROLLER_RELATIVE
        sampler = source_root / SAMPLER_RELATIVE
        source_launcher = source_root / SOURCE_LAUNCHER_RELATIVE
        binary = build / BINARY_NAME
        for path in (controller, sampler, source_launcher, binary):
            if path.resolve(strict=True) != path:
                fail("campaign path is not canonical: " + str(path))
        installed = Path(__file__).resolve(strict=True)
        if installed != INSTALLED_LAUNCHER:
            fail("execute mode is not the fixed installed launcher")
        require_root_owned_nonwritable_chain(installed)
        installed_hash, installed_info = hash_path(installed, 4 * 1024 * 1024, "installed launcher")
        source_hash, source_info = hash_path(source_launcher, 4 * 1024 * 1024, "source launcher")
        if (
            installed_hash != source_hash or installed_info.st_uid != 0
            or installed_info.st_gid != 0 or stat.S_IMODE(installed_info.st_mode) != 0o555
            or source_info.st_uid != CAMPAIGN_UID or source_info.st_gid != CAMPAIGN_GID
            or stat.S_IMODE(source_info.st_mode) != 0o444
        ):
            fail("installed/source launcher seal differs")
        # Every deterministic provenance gate runs before the one irreversible
        # operation (LAUNCH_DIR mkdir).  Preflight repeats these observations
        # after the sampler is live and is cross-bound below.
        git_hash = static_git_gate(source_root, config.expected_commit)
        source_manifest_hash, source_hashes = static_source_manifest(source_root)
        build_manifest_hash, build_profile = validate_build_profile(
            build, source_root, config.expected_commit,
        )
        binary_hash, binary_info = hash_path(
            binary, 64 * 1024 * 1024, "static screen binary"
        )
        if (
            binary_info.st_uid != CAMPAIGN_UID
            or binary_info.st_gid != CAMPAIGN_GID
            or stat.S_IMODE(binary_info.st_mode) != 0o555
        ):
            fail("static screen binary policy differs")
        controller_fd, controller_binding = open_held_source(
            controller, "controller source"
        )
        self.held_sources["controller"] = controller_fd
        sampler_fd, sampler_binding = open_held_source(
            sampler, "sampler source"
        )
        self.held_sources["sampler"] = sampler_fd
        if (
            controller_binding["sha256"] != source_hashes[CONTROLLER_RELATIVE]
            or sampler_binding["sha256"] != source_hashes[SAMPLER_RELATIVE]
            or installed_hash != source_hashes[SOURCE_LAUNCHER_RELATIVE]
        ):
            fail("held source/static source-manifest binding differs")
        python_hash, python_info = hash_path(
            PYTHON.resolve(strict=True), 64 * 1024 * 1024, "root Python",
            require_single_link=False,
        )
        setpriv_hash, setpriv_info = hash_path(
            SETPRIV, 64 * 1024 * 1024, "setpriv"
        )
        env_hash, env_info = hash_path(
            ENV_EXECUTABLE, 64 * 1024 * 1024, "root environment scrubber"
        )
        if (
            stat.S_IMODE(python_info.st_mode) != 0o755
            or python_info.st_uid != 0 or python_info.st_gid != 0
            or stat.S_IMODE(setpriv_info.st_mode) != 0o755
            or setpriv_info.st_uid != 0 or setpriv_info.st_gid != 0
            or stat.S_IMODE(env_info.st_mode) != 0o755
            or env_info.st_uid != 0 or env_info.st_gid != 0
        ):
            fail("root Python/setpriv image policy differs")
        i2c_devices = i2c_device_receipts()
        i2c_pre_scan = scan_i2c_holders(i2c_devices)
        require_no_i2c_holders(i2c_pre_scan)
        for path in (LAUNCH_DIR, SAMPLER_DIR, FIXED_CLAIM_PATH, FIXED_OUTPUT_DIR):
            try:
                os.stat(str(path), follow_symlinks=False)
            except FileNotFoundError:
                continue
            fail("one-shot namespace already exists: " + str(path))
        self.lock_fd = os.open(
            str(LOCK_PATH), os.O_RDWR | os.O_CREAT | os.O_CLOEXEC
            | getattr(os, "O_NOFOLLOW", 0), 0o600,
        )
        fcntl.flock(self.lock_fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        lock_info = os.fstat(self.lock_fd)
        if (
            not stat.S_ISREG(lock_info.st_mode) or lock_info.st_nlink != 1
            or lock_info.st_uid != 0 or lock_info.st_gid != 0
            or stat.S_IMODE(lock_info.st_mode) != 0o600
        ):
            fail("campaign lock policy differs")
        proc_cgroup = self._read_proc_bounded(
            Path("/proc/self/cgroup"), 64 * 1024,
        )
        relative = parse_proc_cgroup(proc_cgroup)
        require_cgroup2(CGROUP_ROOT)
        service = CGROUP_ROOT / relative.lstrip("/")
        if service.resolve(strict=True) != service or service.name != UNIT_NAME:
            fail("delegated cgroup is not the exact transient service")
        layout = CgroupTree(service).create()
        self.layout = layout
        return Prepared(
            config=config, source_root=source_root, controller=controller,
            sampler_script=sampler, binary=binary, launcher_sha256=installed_hash,
            cgroup_path=service, cgroup_relative=relative,
            controller_fd=controller_fd,
            controller_sha256=controller_binding["sha256"],
            controller_bytes=controller_binding["bytes"],
            sampler_fd=sampler_fd,
            sampler_sha256=sampler_binding["sha256"],
            sampler_bytes=sampler_binding["bytes"],
            python_sha256=python_hash, setpriv_sha256=setpriv_hash,
            env_sha256=env_hash, git_sha256=git_hash,
            binary_sha256=binary_hash,
            build_manifest_sha256=build_manifest_hash,
            build_profile=build_profile,
            build_authority=build_authority,
            build_authority_sha256=sha256_bytes(build_authority_raw),
            source_manifest_sha256=source_manifest_hash,
            i2c_devices=i2c_devices, i2c_pre_scan=i2c_pre_scan,
        )

    def reserve(self, prepared: Prepared) -> AttemptJournal:
        attempt = self_hashed(ATTEMPT_SCHEMA, {
            "attempt_id": ATTEMPT_ID, "build_dir": str(prepared.config.build_dir),
            "expected_source_commit": prepared.config.expected_commit,
            "installed_launcher": str(INSTALLED_LAUNCHER),
            "launcher_sha256": prepared.launcher_sha256,
            "service_cgroup": prepared.cgroup_relative,
            "started_monotonic_ns": time.monotonic_ns(), "started_utc": utc_now(),
        })
        def consumed() -> None:
            self.attempt_consumed = True

        try:
            self.journal = AttemptJournal.reserve(attempt, consumed)
        except AttemptConsumedError as exc:
            self.attempt_consumed = True
            if exc.journal is not None:
                self.journal = exc.journal
                self.sampler_dir_state = "absent"
            raise
        parent_fd = self.journal.parent_fd
        self.sampler_dir_state = "absent"
        try:
            try:
                os.stat(
                    SAMPLER_DIR.name, dir_fd=parent_fd,
                    follow_symlinks=False,
                )
            except FileNotFoundError:
                pass
            else:
                self.sampler_dir_state = "unsafe"
                fail("sampler evidence namespace already exists")
            self.sampler_dir_state = "creating"
            os.mkdir(SAMPLER_DIR.name, 0o700, dir_fd=parent_fd)
            self.sampler_dir_state = "created"
            self.sampler_dir_fd = os.open(
                SAMPLER_DIR.name,
                os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC
                | getattr(os, "O_NOFOLLOW", 0),
                dir_fd=parent_fd,
            )
            info = os.fstat(self.sampler_dir_fd)
            named = os.stat(
                SAMPLER_DIR.name, dir_fd=parent_fd, follow_symlinks=False,
            )
            if (
                not stat.S_ISDIR(info.st_mode)
                or (info.st_dev, info.st_ino) != (named.st_dev, named.st_ino)
                or stat.S_IMODE(info.st_mode) != 0o700
                or info.st_uid != 0 or info.st_gid != 0
            ):
                self.sampler_dir_state = "unsafe"
                fail("new sampler evidence directory policy differs")
            self.sampler_dir_binding = (info.st_dev, info.st_ino)
            self.sampler_dir_state = "bound"
            os.chown(SAMPLER_DIR.name, CAMPAIGN_UID, CAMPAIGN_GID,
                     dir_fd=parent_fd, follow_symlinks=False)
            os.fsync(self.sampler_dir_fd)
            os.fsync(parent_fd)
            info = os.fstat(self.sampler_dir_fd)
            named = os.stat(
                SAMPLER_DIR.name, dir_fd=parent_fd, follow_symlinks=False,
            )
            if (
                not stat.S_ISDIR(info.st_mode)
                or (info.st_dev, info.st_ino) != self.sampler_dir_binding
                or (info.st_dev, info.st_ino) != (named.st_dev, named.st_ino)
                or stat.S_IMODE(info.st_mode) != 0o700
                or info.st_uid != CAMPAIGN_UID or info.st_gid != CAMPAIGN_GID
            ):
                fail("sampler evidence directory policy differs")
            return self.journal
        except BaseException:
            # The attempt is already irrevocably consumed.  Recover the exact
            # directory if mkdir returned before a signal/trace exception, or
            # if a later setup operation failed, so terminal accounting does
            # not depend on an in-memory assignment immediately after mkdir.
            if self.sampler_dir_state != "unsafe":
                self._recover_sampler_dir_after_reservation(parent_fd)
            raise

    def _recover_sampler_dir_after_reservation(self, parent_fd: int) -> None:
        candidate_fd = self.sampler_dir_fd
        opened_here = False
        if candidate_fd < 0:
            try:
                candidate_fd = os.open(
                    SAMPLER_DIR.name,
                    os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC
                    | getattr(os, "O_NOFOLLOW", 0),
                    dir_fd=parent_fd,
                )
                opened_here = True
            except FileNotFoundError:
                self.sampler_dir_state = "absent"
                return
            except BaseException:
                self.sampler_dir_state = "unsafe"
                return
        try:
            info = os.fstat(candidate_fd)
            named = os.stat(
                SAMPLER_DIR.name, dir_fd=parent_fd,
                follow_symlinks=False,
            )
            binding = (info.st_dev, info.st_ino)
            if (
                not stat.S_ISDIR(info.st_mode)
                or binding != (named.st_dev, named.st_ino)
                or stat.S_IMODE(info.st_mode) != 0o700
                or (info.st_uid, info.st_gid)
                not in ((0, 0), (CAMPAIGN_UID, CAMPAIGN_GID))
                or os.listdir(candidate_fd)
            ):
                self.sampler_dir_state = "unsafe"
                return
            self.sampler_dir_fd = candidate_fd
            self.sampler_dir_binding = binding
            self.sampler_dir_state = "bound"
            opened_here = False
        except BaseException:
            self.sampler_dir_state = "unsafe"
        finally:
            if opened_here:
                try:
                    os.close(candidate_fd)
                except OSError:
                    pass

    @staticmethod
    def _pipe() -> Tuple[int, int]:
        read_fd, write_fd = os.pipe2(os.O_CLOEXEC)
        os.set_blocking(read_fd, False)
        return read_fd, write_fd

    def _spawn_exec(self, role: str, cgroup: Path, cpus: Iterable[int],
                    argv: Sequence[str],
                    held_source: Optional[Tuple[int, str, str, int]] = None,
                    ) -> ProcessHandle:
        out_r, out_w = self._pipe()
        try:
            err_r, err_w = self._pipe()
            gate_r, gate_w = os.pipe2(os.O_CLOEXEC)
        except BaseException:
            os.close(out_r)
            os.close(out_w)
            raise
        try:
            pid = os.fork()
        except BaseException:
            for fd in (out_r, out_w, err_r, err_w, gate_r, gate_w):
                os.close(fd)
            raise
        if pid == 0:
            try:
                os.close(out_r)
                os.close(err_r)
                os.close(gate_w)
                os.dup2(out_w, 1, inheritable=True)
                os.dup2(err_w, 2, inheritable=True)
                if out_w not in (1, 2):
                    os.close(out_w)
                if err_w not in (1, 2):
                    os.close(err_w)
                CgroupTree.move_pid(cgroup, os.getpid(), cpus)
                # Do not execute an unprivileged role until the root parent has
                # acquired its pidfd, checked identity, and repeated containment.
                if os.read(gate_r, 1) != b"R":
                    fail(role + " startup release differs")
                os.close(gate_r)
                sealed = None
                if held_source is not None:
                    source_fd, source_path, source_hash, source_size = held_source
                    os.dup2(source_fd, SEALED_SOURCE_FD, inheritable=True)
                    sealed = (
                        SEALED_SOURCE_FD, source_path, source_hash,
                        source_size, "main",
                    )
                command = setpriv_exec_argv(
                    argv, sampler=(role == "sampler"), sealed_source=sealed,
                )
                os.execve(str(SETPRIV), command, dict(SEALED_ENVIRONMENT))
            except BaseException as exc:
                child_exec_failure(role + " exec: " + exception_text(exc))
        os.close(gate_r)
        os.close(out_w)
        os.close(err_w)
        pidfd = -1
        try:
            pidfd = os.pidfd_open(pid, 0)
            start = proc_start_ticks(self._read_proc_bounded(
                Path(f"/proc/{pid}/stat"), 64 * 1024,
            ))
            CgroupTree.move_pid(cgroup, pid, cpus)
            if os.write(gate_w, b"R") != 1:
                fail(role + " startup release was short")
            os.close(gate_w)
            gate_w = -1
        except BaseException as exc:
            if gate_w >= 0:
                os.close(gate_w)
            try:
                reap_failed_spawn(pid, pidfd)
            except BaseException as cleanup_exc:
                exc = LaunchError(
                    exception_text(exc) + " | failed-spawn cleanup: "
                    + exception_text(cleanup_exc)
                )
            os.close(out_r)
            os.close(err_r)
            raise exc
        handle = ProcessHandle(role, pid, pidfd, start, out_r, err_r)
        self.children.append(handle)
        return handle

    def start_sampler(self, prepared: Prepared) -> ProcessHandle:
        if self.layout is None:
            fail("cgroup layout is absent")
        sampler_hash = hash_path(
            prepared.sampler_script, 1024 * 1024, "sampler source launch"
        )[0]
        return self._spawn_exec(
            "sampler", self.layout.sampler, (SAMPLER_CPU,),
            sampler_argv(prepared, sampler_hash),
        )

    def _process_exited(self, process: ProcessHandle) -> bool:
        poller = select.poll()
        poller.register(process.pidfd, select.POLLIN)
        return bool(poller.poll(0))

    def wait_sampler_ready(self, process: ProcessHandle,
                           prepared: Prepared) -> Dict[str, Any]:
        deadline = time.monotonic() + SAMPLER_READY_SECONDS
        while time.monotonic() < deadline:
            if self._process_exited(process):
                fail("sampler exited before readiness")
            opened: Dict[str, int] = {}
            try:
                for name, path in (
                    ("csv", SAMPLER_CSV), ("pid", SAMPLER_PID_FILE),
                    ("validation", SAMPLER_VALIDATION), ("receipt", SAMPLER_RECEIPT),
                ):
                    opened[name] = os.open(
                        str(path), os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK
                        | getattr(os, "O_NOFOLLOW", 0),
                    )
                pid_raw = os.pread(opened["pid"], 64, 0)
                if pid_raw != (str(process.pid) + "\n").encode("ascii"):
                    fail("sampler PID file differs")
                if os.fstat(opened["csv"]).st_size == 0 or os.fstat(opened["validation"]).st_size == 0:
                    raise FileNotFoundError("sampler streams not initialized")
                # This scan is an authenticated point-in-time observation, not
                # a claim of global or kernel-enforced device exclusion.
                devices = i2c_device_receipts()
                if devices != prepared.i2c_devices:
                    fail("I2C device nodes changed between root scans")
                post_scan = scan_i2c_holders(devices)
                require_sampler_sole_i2c_holder(post_scan, process)
                self.held_sampler = opened
                return post_scan
            except FileNotFoundError:
                for fd in opened.values():
                    os.close(fd)
                time.sleep(0.02)
        fail("sampler readiness deadline expired")

    def _communicate_exact(self, process: ProcessHandle, deadline: float,
                           stdout_limit: int, stderr_limit: int) -> Tuple[bytes, bytes, int]:
        out = bytearray()
        err = bytearray()
        streams: Dict[int, Tuple[bytearray, int]] = {
            process.stdout_fd: (out, stdout_limit),
            process.stderr_fd: (err, stderr_limit),
        }
        poller = select.poll()
        for fd in streams:
            poller.register(fd, select.POLLIN | select.POLLHUP | select.POLLERR)
        while process.returncode is None or streams:
            if time.monotonic() >= deadline:
                fail(process.role + " deadline expired")
            for fd, _event in poller.poll(50):
                if fd not in streams:
                    continue
                target, limit = streams[fd]
                while True:
                    try:
                        block = os.read(fd, 65536)
                    except BlockingIOError:
                        break
                    if not block:
                        poller.unregister(fd)
                        os.close(fd)
                        streams.pop(fd)
                        if process.stdout_fd == fd:
                            process.stdout_fd = -1
                        if process.stderr_fd == fd:
                            process.stderr_fd = -1
                        break
                    target.extend(block)
                    if len(target) > limit:
                        fail(process.role + " output exceeded its bound")
            if process.returncode is None:
                waited, status = os.waitpid(process.pid, os.WNOHANG)
                if waited == process.pid:
                    process.returncode = wait_status_code(status)
        return bytes(out), bytes(err), process.returncode

    def run_preflight(self, prepared: Prepared, sampler: ProcessHandle) -> PreflightAuthority:
        if self.layout is None:
            fail("cgroup layout is absent")
        process = self._spawn_exec(
            "preflight", self.layout.experiment, PREFLIGHT_CPUS,
            preflight_argv(prepared, sampler.pid),
            (
                prepared.controller_fd, str(prepared.controller),
                prepared.controller_sha256, prepared.controller_bytes,
            ),
        )
        try:
            stdout, stderr, rc = self._communicate_exact(
                process, time.monotonic() + PREFLIGHT_SECONDS,
                MAX_PREFLIGHT_BYTES, MAX_STDERR_BYTES,
            )
        except BaseException:
            self.kill_experiment()
            self.reap(process, 2.0)
            raise
        if rc != 0 or stderr:
            fail("preflight process failed: " + stderr.decode("ascii", "replace")[:4000])
        if CgroupTree.pids(self.layout.experiment):
            fail("experiment cgroup is not empty after preflight")
        return parse_preflight_receipt(stdout, prepared, sampler.pid)

    def _spawn_stub(self, role: str, cgroup: Path, cpus: Iterable[int],
                    authority: PreflightAuthority, python_argv: Sequence[str],
                    prepared: Prepared) -> ProcessHandle:
        parent_sock, child_sock = socket.socketpair(
            socket.AF_UNIX, socket.SOCK_SEQPACKET | socket.SOCK_CLOEXEC,
        )
        try:
            out_r, out_w = self._pipe()
            err_r, err_w = self._pipe()
        except BaseException:
            parent_sock.close()
            child_sock.close()
            raise
        try:
            pid = os.fork()
        except BaseException:
            parent_sock.close()
            child_sock.close()
            for fd in (out_r, out_w, err_r, err_w):
                os.close(fd)
            raise
        if pid == 0:
            try:
                parent_sock.close()
                os.close(out_r)
                os.close(err_r)
                os.dup2(out_w, 1, inheritable=True)
                os.dup2(err_w, 2, inheritable=True)
                CgroupTree.move_pid(cgroup, os.getpid(), cpus)
                receipt_fd, command_fd = receive_fds(child_sock)
                try:
                    preflight_raw = read_sealed_memfd(
                        receipt_fd, MAX_PREFLIGHT_BYTES, role + " preflight"
                    )
                    command_raw = read_sealed_memfd(
                        command_fd, MAX_PREFLIGHT_BYTES, role + " command"
                    )
                    command = validate_role_document(command_raw, role, preflight_raw)
                finally:
                    os.close(receipt_fd)
                    os.close(command_fd)
                if child_sock.send(b"K") != 1:
                    fail("role ACK write was short")
                action = child_sock.recv(1)
                if action == b"C":
                    os._exit(124)
                if action != b"R":
                    fail("role release command differs")
                child_sock.close()
                os.dup2(prepared.controller_fd, SEALED_SOURCE_FD, inheritable=True)
                command_argv = setpriv_exec_argv(
                    command["python_argv"], sampler=False,
                    sealed_source=(
                        SEALED_SOURCE_FD, str(prepared.controller),
                        prepared.controller_sha256,
                        prepared.controller_bytes, "main",
                    ),
                )
                os.execve(str(SETPRIV), command_argv, dict(SEALED_ENVIRONMENT))
            except BaseException as exc:
                child_exec_failure(role + " stub: " + exception_text(exc))
        child_sock.close()
        os.close(out_w)
        os.close(err_w)
        pidfd = -1
        try:
            pidfd = os.pidfd_open(pid, 0)
            start = proc_start_ticks(self._read_proc_bounded(
                Path(f"/proc/{pid}/stat"), 64 * 1024,
            ))
            CgroupTree.move_pid(cgroup, pid, cpus)
            receipt_fd = -1
            command_fd = -1
            try:
                receipt_fd = make_sealed_memfd("wh2-preflight", authority.raw)
                command_raw = role_document(role, python_argv, authority, prepared)
                command_fd = make_sealed_memfd("wh2-" + role, command_raw)
                send_fds(parent_sock, (receipt_fd, command_fd))
            finally:
                if receipt_fd >= 0:
                    os.close(receipt_fd)
                if command_fd >= 0:
                    os.close(command_fd)
            parent_sock.settimeout(5.0)
            if parent_sock.recv(1) != b"K":
                fail(role + " stub did not ACK authority")
            parent_sock.settimeout(None)
        except BaseException as exc:
            parent_sock.close()
            try:
                reap_failed_spawn(pid, pidfd)
            except BaseException as cleanup_exc:
                exc = LaunchError(
                    exception_text(exc) + " | failed-spawn cleanup: "
                    + exception_text(cleanup_exc)
                )
            os.close(out_r)
            os.close(err_r)
            raise exc
        process = ProcessHandle(
            role, pid, pidfd, start, out_r, err_r, parent_sock,
        )
        self.children.append(process)
        return process

    def spawn_stubs(self, prepared: Prepared, authority: PreflightAuthority
                    ) -> Tuple[ProcessHandle, ProcessHandle]:
        if self.layout is None:
            fail("cgroup layout is absent")
        controller = self._spawn_stub(
            "controller", self.layout.experiment, PREFLIGHT_CPUS,
            authority, authority.run_argv, prepared,
        )
        verifier_argv = [
            str(prepared.controller), "--verify-retained",
            "--expected-run-argv-sha256", authority.run_argv_sha256,
        ]
        verifier = self._spawn_stub(
            "verifier", self.layout.verifier, (VERIFIER_CPU,),
            authority, verifier_argv, prepared,
        )
        return controller, verifier

    def spawn_guardian(self, journal: AttemptJournal) -> GuardianHandle:
        if self.layout is None:
            fail("cgroup layout is absent")
        parent_sock, child_sock = socket.socketpair(
            socket.AF_UNIX, socket.SOCK_SEQPACKET | socket.SOCK_CLOEXEC,
        )
        try:
            pid = os.fork()
        except BaseException:
            parent_sock.close()
            child_sock.close()
            raise
        if pid == 0:
            try:
                parent_sock.close()
                CgroupTree.move_pid(
                    self.layout.supervisor, os.getpid(), (VERIFIER_CPU,)
                )
                guardian_child(
                    child_sock, self.layout.experiment_kill_fd,
                    self.layout.run_kill_fd, self.layout.experiment,
                    self.layout.run, journal.directory_fd,
                )
            except BaseException as exc:
                child_exec_failure("guardian process: " + exception_text(exc))
        child_sock.close()
        pidfd = -1
        try:
            pidfd = os.pidfd_open(pid, 0)
            start = proc_start_ticks(
                self._read_proc_bounded(
                    Path(f"/proc/{pid}/stat"), 64 * 1024,
                )
            )
            CgroupTree.move_pid(
                self.layout.supervisor, pid, (VERIFIER_CPU,)
            )
            process = ProcessHandle("guardian", pid, pidfd, start)
            handle = GuardianHandle(process, parent_sock)
            handle.initial_ack()
        except BaseException as exc:
            parent_sock.close()
            try:
                reap_failed_spawn(pid, pidfd)
            except BaseException as cleanup_exc:
                exc = LaunchError(
                    exception_text(exc) + " | failed-spawn cleanup: "
                    + exception_text(cleanup_exc)
                )
            raise exc
        self.children.append(process)
        return handle

    @staticmethod
    def release(process: ProcessHandle) -> int:
        if process.control is None:
            fail(process.role + " has no release channel")
        released = time.monotonic_ns()
        if process.control.send(b"R") != 1:
            fail(process.role + " release was short")
        process.control.close()
        process.control = None
        return released

    @staticmethod
    def cancel(process: ProcessHandle) -> None:
        if process.control is not None:
            try:
                process.control.send(b"C")
            except OSError:
                pass
            process.control.close()
            process.control = None

    def wait_claim(self, authority: PreflightAuthority, controller: ProcessHandle,
                   release_ns: int) -> Tuple[Dict[str, Any], bytes, int]:
        deadline_ns = release_ns + int(CLAIM_ADMISSION_SECONDS * 1e9)
        while time.monotonic_ns() < deadline_ns:
            try:
                raw = self._read_claim_if_ready()
            except FileNotFoundError:
                raw = None
            if raw is not None:
                observed = time.monotonic_ns()
                claim, started = parse_claim(
                    raw, authority, controller.pid, controller.start_ticks,
                    release_ns, observed,
                )
                return claim, raw, started
            if self._process_exited(controller):
                fail("controller exited before publishing a valid claim")
            time.sleep(0.001)
        fail("controller claim admission deadline expired")

    def kill_experiment(self) -> List[int]:
        if self.layout is None:
            return []
        roster = CgroupTree.pids(self.layout.experiment)
        if os.write(self.layout.experiment_kill_fd, b"1") != 1:
            fail("experiment cgroup.kill write was short")
        return roster

    def ensure_experiment_empty(self, controller: ProcessHandle
                                ) -> Tuple[bool, List[int]]:
        if self.layout is None:
            fail("cgroup layout is absent")
        roster = CgroupTree.pids(self.layout.experiment)
        forced = bool(roster)
        if forced:
            if os.write(self.layout.experiment_kill_fd, b"1") != 1:
                fail("experiment cleanup cgroup.kill write was short")
        deadline = time.monotonic() + 2.0
        while CgroupTree.pids(self.layout.experiment):
            if time.monotonic() >= deadline:
                fail("experiment cgroup did not become empty")
            time.sleep(0.005)
        for pid in roster:
            if pid == controller.pid:
                continue
            try:
                waited, _status = os.waitpid(pid, os.WNOHANG)
                if waited == 0:
                    while time.monotonic() < deadline:
                        waited, _status = os.waitpid(pid, os.WNOHANG)
                        if waited == pid:
                            break
                        time.sleep(0.005)
            except ChildProcessError:
                pass
        return forced, roster

    def wait_process(self, process: ProcessHandle, deadline: float) -> int:
        while process.returncode is None:
            waited, status = os.waitpid(process.pid, os.WNOHANG)
            if waited == process.pid:
                process.returncode = wait_status_code(status)
                break
            if time.monotonic() >= deadline:
                fail(process.role + " reap deadline expired")
            time.sleep(0.005)
        return process.returncode

    def reap(self, process: ProcessHandle, seconds: float) -> int:
        return self.wait_process(process, time.monotonic() + seconds)

    @staticmethod
    def _drain_fd(fd: int, limit: int) -> bytes:
        if fd < 0:
            return b""
        data = bytearray()
        while True:
            try:
                block = os.read(fd, 65536)
            except BlockingIOError:
                break
            if not block:
                break
            data.extend(block)
            if len(data) > limit:
                fail("child diagnostic exceeded bound")
        return bytes(data)

    def controller_wait(self, process: ProcessHandle, deadline_ns: int) -> Tuple[int, bytes, bytes]:
        rc = self.wait_process(process, deadline_ns / 1e9 + 2.0)
        out = self._drain_fd(process.stdout_fd, 4096)
        err = self._drain_fd(process.stderr_fd, MAX_STDERR_BYTES)
        return rc, out, err

    def stop_sampler(self, process: ProcessHandle, prepared: Prepared
                     ) -> Tuple[int, Dict[str, Any], bytes, bytes]:
        if process.returncode is None:
            try:
                signal.pidfd_send_signal(process.pidfd, signal.SIGTERM)
            except ProcessLookupError:
                pass
        try:
            rc = self.wait_process(process, time.monotonic() + SAMPLER_STOP_SECONDS)
        except BaseException:
            if self.layout is not None:
                try:
                    os.write(self.layout.sampler_kill_fd, b"1")
                except OSError:
                    pass
            self.reap(process, 2.0)
            raise
        # The admission-time receipt FD is authoritative.  Reading the mutable
        # pathname here would let a later same-UID replacement hide the sealed
        # terminal bytes before the identity check below can diagnose it.
        receipt_raw = read_held_file(
            self.held_sampler["receipt"], 1024 * 1024,
            "sampler terminal receipt",
        )
        terminal = validate_sampler_terminal(
            receipt_raw, process.pid, prepared, self.held_sampler, rc,
        )
        stdout = self._drain_fd(process.stdout_fd, MAX_STDERR_BYTES)
        stderr = self._drain_fd(process.stderr_fd, MAX_STDERR_BYTES)
        if stdout or stderr:
            fail("sampler emitted unexpected stdout/stderr")
        terminal["raw"] = receipt_raw
        return rc, terminal, stdout, stderr

    def run_verifier(self, process: ProcessHandle) -> Tuple[str, int, bytes, bytes]:
        self.release(process)
        stdout, stderr, rc = self._communicate_exact(
            process, time.monotonic() + VERIFIER_SECONDS, 4096, MAX_STDERR_BYTES,
        )
        if stderr:
            fail("retained verifier emitted stderr")
        if self.layout is None or cgroup_descendant_pids(self.layout.verifier):
            fail("retained verifier cgroup remained populated")
        return parse_verifier_output(stdout, rc), rc, stdout, stderr

    def contain_verifier(self, process: ProcessHandle) -> bool:
        """Kill, reap, and close a verifier whose release/result is ambiguous."""
        if self.layout is None:
            fail("cgroup layout is absent")
        self.cancel(process)
        roster = cgroup_descendant_pids(self.layout.verifier)
        if roster:
            if os.write(self.layout.verifier_kill_fd, b"1") != 1:
                fail("verifier containment cgroup.kill write was short")
        deadline = time.monotonic() + 2.0
        while cgroup_descendant_pids(self.layout.verifier):
            if time.monotonic() >= deadline:
                fail("verifier cgroup did not become empty")
            time.sleep(0.005)
        if process.returncode is None:
            self.wait_process(process, deadline)
        failures: List[str] = []
        for member, limit in (
            ("stdout_fd", 4096), ("stderr_fd", MAX_STDERR_BYTES),
        ):
            fd = getattr(process, member)
            if fd < 0:
                continue
            try:
                self._drain_fd(fd, limit)
            except BaseException as exc:
                failures.append(member + " drain: " + exception_text(exc))
            try:
                os.close(fd)
            except OSError as exc:
                failures.append(member + " close: " + exception_text(exc))
            setattr(process, member, -1)
        if cgroup_descendant_pids(self.layout.verifier):
            failures.append("verifier leaf repopulated after containment")
        if process.returncode is None:
            failures.append("verifier process remained unreaped")
        if failures:
            fail("verifier containment failed: " + " | ".join(failures[:8]))
        return True

    def finish_sampler_dir(self) -> None:
        if self.journal is None:
            fail("sampler directory has no attempt authority")
        parent_fd = self.journal.parent_fd
        if self.sampler_dir_state == "absent":
            try:
                os.stat(
                    SAMPLER_DIR.name, dir_fd=parent_fd,
                    follow_symlinks=False,
                )
            except FileNotFoundError:
                os.fsync(parent_fd)
                return
            fail("absent sampler directory acquired an unexpected name")
        if (
            self.sampler_dir_state != "bound" or self.sampler_dir_fd < 0
            or self.sampler_dir_binding is None
        ):
            fail("sampler directory was not safely bound")
        before = os.fstat(self.sampler_dir_fd)
        named = os.stat(
            SAMPLER_DIR.name, dir_fd=parent_fd, follow_symlinks=False,
        )
        if (
            not stat.S_ISDIR(before.st_mode)
            or (before.st_dev, before.st_ino) != self.sampler_dir_binding
            or (before.st_dev, before.st_ino) != (named.st_dev, named.st_ino)
            or before.st_uid not in (0, CAMPAIGN_UID)
            or before.st_gid not in (0, CAMPAIGN_GID)
            or stat.S_IMODE(before.st_mode) not in (0o700, 0o500)
        ):
            fail("bound sampler evidence directory drifted")
        os.fchown(self.sampler_dir_fd, CAMPAIGN_UID, CAMPAIGN_GID)
        os.fchmod(self.sampler_dir_fd, 0o500)
        os.fsync(self.sampler_dir_fd)
        info = os.fstat(self.sampler_dir_fd)
        named = os.stat(
            SAMPLER_DIR.name, dir_fd=parent_fd, follow_symlinks=False,
        )
        if (
            not stat.S_ISDIR(info.st_mode) or stat.S_IMODE(info.st_mode) != 0o500
            or info.st_uid != CAMPAIGN_UID or info.st_gid != CAMPAIGN_GID
            or (info.st_dev, info.st_ino) != self.sampler_dir_binding
            or (info.st_dev, info.st_ino) != (named.st_dev, named.st_ino)
        ):
            fail("final sampler evidence directory policy differs")
        os.fsync(parent_fd)
        os.close(self.sampler_dir_fd)
        self.sampler_dir_fd = -1
        self.sampler_dir_state = "closed"

    def quiesce(self) -> bool:
        failures: List[str] = []
        if self.layout is not None and cgroup_descendant_pids(self.layout.run):
            try:
                if os.write(self.layout.run_kill_fd, b"1") != 1:
                    fail("run cleanup cgroup.kill write was short")
            except BaseException as exc:
                failures.append("run kill: " + exception_text(exc))
        for process in self.children:
            if process.role == "guardian":
                continue
            if process.control is not None:
                self.cancel(process)
            if process.returncode is None:
                try:
                    pidfd_signal(process.pidfd, signal.SIGKILL)
                except (OSError, LaunchError) as exc:
                    failures.append(
                        process.role + " pidfd kill: " + exception_text(exc)
                    )
                try:
                    self.reap(process, 1.0)
                except BaseException as exc:
                    failures.append(
                        process.role + " reap: " + exception_text(exc)
                    )
        if self.layout is not None:
            deadline = time.monotonic() + 2.0
            roster = cgroup_descendant_pids(self.layout.run)
            while roster and time.monotonic() < deadline:
                time.sleep(0.005)
                roster = cgroup_descendant_pids(self.layout.run)
            if roster:
                failures.append("run cgroup remains populated: " + repr(roster))
        if failures:
            fail("runtime quiescence failed: " + " | ".join(failures[:16]))
        return True

    def close_runtime(self) -> bool:
        failures: List[str] = []
        for process in self.children:
            if process.returncode is None:
                try:
                    pidfd_signal(process.pidfd, signal.SIGKILL)
                except (OSError, LaunchError) as exc:
                    failures.append(
                        process.role + " final pidfd kill: " + exception_text(exc)
                    )
                try:
                    self.reap(process, 1.0)
                except BaseException as exc:
                    failures.append(
                        process.role + " final reap: " + exception_text(exc)
                    )
            if process.returncode is None:
                failures.append(process.role + " remained unreaped")
            for fd in (process.pidfd, process.stdout_fd, process.stderr_fd):
                if fd >= 0:
                    try:
                        os.close(fd)
                    except OSError as exc:
                        failures.append(
                            process.role + " descriptor close: " + exception_text(exc)
                        )
            if process.control is not None:
                try:
                    process.control.close()
                except OSError as exc:
                    failures.append(
                        process.role + " control close: " + exception_text(exc)
                    )
        while True:
            try:
                pid, status = os.waitpid(-1, os.WNOHANG)
            except ChildProcessError:
                break
            except OSError as exc:
                failures.append("subreaper sweep: " + exception_text(exc))
                break
            if pid == 0:
                failures.append("one or more adopted children remained live")
                break
            # An adopted, untracked descendant is a containment fault, but its
            # successful reap is still an exact closure proof.  The earlier
            # recursive run-empty gate is what prevents publication while it
            # remains alive.
            wait_status_code(status)
        for fd in self.held_sampler.values():
            try:
                os.close(fd)
            except OSError as exc:
                failures.append("sampler evidence close: " + exception_text(exc))
        self.held_sampler.clear()
        for fd in self.held_sources.values():
            try:
                os.close(fd)
            except OSError as exc:
                failures.append("held source close: " + exception_text(exc))
        self.held_sources.clear()
        if self.sampler_dir_fd >= 0:
            try:
                os.close(self.sampler_dir_fd)
                self.sampler_dir_fd = -1
            except OSError as exc:
                failures.append(
                    "sampler directory descriptor close: " + exception_text(exc)
                )
        if self.layout is not None:
            for fd in (
                self.layout.experiment_kill_fd, self.layout.sampler_kill_fd,
                self.layout.verifier_kill_fd, self.layout.run_kill_fd,
            ):
                try:
                    os.close(fd)
                except OSError as exc:
                    failures.append("cgroup descriptor close: " + exception_text(exc))
        if failures:
            fail("runtime closure failed: " + " | ".join(failures[:16]))
        return True

    def close_authority(self) -> bool:
        if self.lock_fd >= 0:
            os.close(self.lock_fd)
            self.lock_fd = -1
        return True


@dataclass
class RunResult:
    classification: str
    launcher_status: str
    exit_code: int
    phases: List[str]


class Supervisor:
    """Dependency-injected one-shot state machine used by real and fake tests."""

    def __init__(self, backend: Any) -> None:
        self.backend = backend
        self.phases: List[str] = []

    def _phase(self, name: str) -> None:
        self.phases.append(name)

    def execute(self, config: LaunchConfig) -> RunResult:
        prepared: Optional[Prepared] = None
        journal: Optional[Any] = None
        sampler: Optional[Any] = None
        controller: Optional[Any] = None
        verifier: Optional[Any] = None
        guardian: Optional[Any] = None
        authority: Optional[PreflightAuthority] = None
        consumed = False
        journal_complete = False
        runtime_closed = False
        quiesced = False
        lifecycle_resolved = False
        sampler_dir_closed = False
        authority_closed = False
        guardian_closed = False
        guardian_finished = False
        guardian_outer_cancelled = False
        sampler_terminalized = False
        sampler_terminal_write_attempted = False
        verifier_ran = False
        verifier_attempted = False
        verifier_resolved = False
        verification_write_attempted = False
        classification = "consumed_unclassifiable"
        errors: List[str] = []
        terminal: Dict[str, Any] = {}
        controller_rc: Optional[int] = None
        controller_stdout = b""
        controller_stderr = b""
        sampler_rc: Optional[int] = None
        sampler_terminal: Dict[str, Any] = {}
        verifier_rc: Optional[int] = None
        verifier_stdout = b""
        verifier_stderr = b""
        release_ns: Optional[int] = None
        controller_t0_ns: Optional[int] = None
        deadline_ns: Optional[int] = None
        exact_outer_deadline_armed = False
        exact_backstop_reanchored = False
        containment_roster: List[int] = []
        i2c_post_scan: Dict[str, Any] = {}
        experiment_empty_for_verifier = False

        def record(message: str) -> None:
            if message not in errors and len(errors) < 64:
                errors.append(message[:4096])

        try:
            self._phase("read_only_gates")
            prepared = self.backend.prepare(config)
            try:
                journal = self.backend.reserve(prepared)
                consumed = True
            except BaseException:
                consumed = bool(getattr(self.backend, "attempt_consumed", False))
                journal = getattr(self.backend, "journal", None)
                raise
            self._phase("attempt_reserved")
            journal.write_bytes(
                "build-authority.json",
                canonical_bytes(prepared.build_authority) + b"\n",
            )
            sampler = self.backend.start_sampler(prepared)
            i2c_post_scan = self.backend.wait_sampler_ready(sampler, prepared)
            self._phase("sampler_ready")
            authority = self.backend.run_preflight(prepared, sampler)
            journal.write_bytes("preflight.json", authority.raw)
            journal.write_json("authority.json", self_hashed(AUTHORITY_SCHEMA, {
                "preflight_receipt_sha256": authority.receipt_sha256,
                "preflight_raw_sha256": sha256_bytes(authority.raw),
                "run_argv": authority.run_argv,
                "run_argv_sha256": authority.run_argv_sha256,
                "held_controller_bytes": prepared.controller_bytes,
                "held_controller_sha256": prepared.controller_sha256,
                "held_sampler_bytes": prepared.sampler_bytes,
                "held_sampler_sha256": prepared.sampler_sha256,
                "i2c_post_open_scan": i2c_post_scan,
                "i2c_pre_open_scan": prepared.i2c_pre_scan,
                "root_binary_sha256": prepared.binary_sha256,
                "root_build_manifest_sha256": prepared.build_manifest_sha256,
                "root_build_profile": prepared.build_profile,
                "root_build_authority_raw_sha256": (
                    prepared.build_authority_sha256
                ),
                "root_build_authority_receipt_sha256": (
                    prepared.build_authority["receipt_sha256"]
                ),
                "root_env_sha256": prepared.env_sha256,
                "root_git_sha256": prepared.git_sha256,
                "root_python_sha256": prepared.python_sha256,
                "setpriv_sha256": prepared.setpriv_sha256,
                "root_source_manifest_sha256": prepared.source_manifest_sha256,
            }))
            self._phase("preflight_sealed")
            controller, verifier = self.backend.spawn_stubs(prepared, authority)
            guardian = self.backend.spawn_guardian(journal)
            release_not_before_ns = time.monotonic_ns()
            guardian.arm_admission(release_not_before_ns)
            start_record = self_hashed(START_SCHEMA, {
                "controller_pid": controller.pid,
                "controller_process_start_ticks": controller.start_ticks,
                "controller_release_not_before_monotonic_ns": release_not_before_ns,
                "guardian_admission_deadline_monotonic_ns": (
                    release_not_before_ns + int(CLAIM_ADMISSION_SECONDS * 1e9)
                ),
                "guardian_pid": guardian.process.pid,
                "guardian_process_start_ticks": guardian.process.start_ticks,
                "i2c_post_open_scan": i2c_post_scan,
                "i2c_pre_open_scan": prepared.i2c_pre_scan,
                "sampler_pid": sampler.pid,
                "sampler_process_start_ticks": sampler.start_ticks,
                "verifier_pid": verifier.pid,
                "verifier_process_start_ticks": verifier.start_ticks,
            })
            journal.write_json("start.json", start_record)
            release_ns = self.backend.release(controller)
            self._phase("controller_released")
            claim_raw = b""
            try:
                _claim, claim_raw, controller_t0_ns = self.backend.wait_claim(
                    authority, controller, release_ns,
                )
                exact_outer_deadline_armed = (
                    guardian.arm_exact(controller_t0_ns) is True
                )
                exact_backstop_reanchored = True
                if exact_outer_deadline_armed:
                    deadline_ns = controller_t0_ns + int(OUTER_SECONDS * 1e9)
                else:
                    # The inclusive admission timer already killed experiment.
                    # Its ACK still proves that +140 was re-anchored to T0.
                    record("guardian admission deadline fired before exact arm")
            except BaseException as exc:
                record("claim admission: " + exception_text(exc))
                try:
                    containment_roster.extend(self.backend.kill_experiment())
                except BaseException as kill_exc:
                    record("claim containment: " + exception_text(kill_exc))
            deadline_value = self_hashed(
                "wirehair.wh2.direct-systematic-complement-launch-deadline.v1",
                {
                    "admission_deadline_monotonic_ns": (
                        release_not_before_ns
                        + int(CLAIM_ADMISSION_SECONDS * 1e9)
                    ),
                    "claim_sha256": sha256_bytes(claim_raw) if claim_raw else None,
                    "controller_release_monotonic_ns": release_ns,
                    "controller_started_monotonic_ns": controller_t0_ns,
                    "exact_backstop_reanchored": exact_backstop_reanchored,
                    "exact_deadline_armed": exact_outer_deadline_armed,
                    "guardian_anchor_monotonic_ns": release_not_before_ns,
                    "experiment_deadline_monotonic_ns": (
                        deadline_ns if deadline_ns is not None else
                        release_not_before_ns
                        + int(CLAIM_ADMISSION_SECONDS * 1e9)
                    ),
                    "whole_run_deadline_monotonic_ns": (
                        controller_t0_ns + int(WHOLE_RUN_SECONDS * 1e9)
                        if exact_backstop_reanchored else
                        release_not_before_ns + int(
                            (CLAIM_ADMISSION_SECONDS + WHOLE_RUN_SECONDS) * 1e9
                        )
                    ),
                },
            )
            journal.write_json("deadline.json", deadline_value)
            self._phase("outer_armed")
            controller_wait_ns = (
                deadline_ns if deadline_ns is not None else
                release_not_before_ns + int(CLAIM_ADMISSION_SECONDS * 1e9)
            )
            try:
                controller_rc, controller_stdout, controller_stderr = (
                    self.backend.controller_wait(controller, controller_wait_ns)
                )
            except BaseException as exc:
                record("controller wait: " + exception_text(exc))
                try:
                    containment_roster.extend(self.backend.kill_experiment())
                    controller_rc, controller_stdout, controller_stderr = (
                        self.backend.controller_wait(
                            controller, time.monotonic_ns() + 2_000_000_000
                        )
                    )
                except BaseException as cleanup_exc:
                    record("controller containment: " + exception_text(cleanup_exc))
            try:
                forced, roster = self.backend.ensure_experiment_empty(controller)
                containment_roster.extend(roster)
                if forced:
                    record("experiment descendants survived controller reap")
                experiment_empty_for_verifier = True
            except BaseException as exc:
                record("experiment empty proof: " + exception_text(exc))
                raise
            # Only an empty experiment leaf licenses cancellation of +120.
            guardian.cancel()
            guardian_outer_cancelled = True
            if guardian.error != "none":
                record("guardian outer: " + guardian.error)
            self._phase("controller_reaped")
            # Retained adjudication runs while the exact sampler process and
            # its held evidence FDs remain live.  Sampler shutdown/authentication
            # is a later launcher concern and cannot erase a verifier outcome.
            verifier_attempted = True
            classification, verifier_rc, verifier_stdout, verifier_stderr = (
                self.backend.run_verifier(verifier)
            )
            verifier_ran = True
            verifier_resolved = True
            self._phase("retained_verified")
            sampler_rc, sampler_terminal, _sampler_out, _sampler_err = (
                self.backend.stop_sampler(sampler, prepared)
            )
            sampler_terminalized = True
            if sampler_terminal.get("source_path_still_bound") is False:
                record("sampler source pathname drifted after static seal")
            sampler_terminal_write_attempted = True
            journal.write_bytes("sampler-terminal.json", sampler_terminal["raw"])
            self._phase("sampler_terminal")
            verification_write_attempted = True
            journal.write_bytes("verification.json", verifier_stdout)
            terminal = {
                "claim_controller_started_monotonic_ns": controller_t0_ns,
                "controller_returncode": controller_rc,
                "controller_stderr_sha256": sha256_bytes(controller_stderr),
                "controller_stdout_sha256": sha256_bytes(controller_stdout),
                "containment_roster": sorted(set(containment_roster)),
                "errors": list(errors),
                "guardian_deadline_monotonic_ns": deadline_ns,
                "guardian_fired": guardian.fired,
                "launcher_status": "clean" if not errors else "fault",
                "outcome": classification,
                "run_argv_sha256": authority.run_argv_sha256,
                "sampler_returncode": sampler_rc,
                "sampler_terminal_receipt_sha256": sampler_terminal["receipt_sha256"],
                "verifier_returncode": verifier_rc,
                "verifier_stderr_sha256": sha256_bytes(verifier_stderr),
                "verifier_stdout_sha256": sha256_bytes(verifier_stdout),
            }
        except BaseException as exc:
            record(exception_text(exc))
            if controller is not None:
                try:
                    containment_roster.extend(self.backend.kill_experiment())
                    if not experiment_empty_for_verifier:
                        forced, roster = self.backend.ensure_experiment_empty(controller)
                        containment_roster.extend(roster)
                        experiment_empty_for_verifier = True
                        if forced:
                            record("experiment required exception-path containment")
                except BaseException as kill_exc:
                    record("experiment kill/empty: " + exception_text(kill_exc))
            if (
                verifier is not None and not verifier_attempted
                and experiment_empty_for_verifier
            ):
                try:
                    verifier_attempted = True
                    classification, verifier_rc, verifier_stdout, verifier_stderr = (
                        self.backend.run_verifier(verifier)
                    )
                    verifier_ran = True
                    verifier_resolved = True
                except BaseException as verify_exc:
                    record("retained verifier: " + exception_text(verify_exc))
            if verifier is not None and not verifier_resolved:
                try:
                    verifier_resolved = (
                        self.backend.contain_verifier(verifier) is True
                    )
                    if not verifier_resolved:
                        fail("verifier containment lacked exact proof")
                except BaseException as verify_cleanup_exc:
                    record(
                        "verifier containment: "
                        + exception_text(verify_cleanup_exc)
                    )
            if (
                sampler is not None and prepared is not None
                and not sampler_terminalized
                and (controller is None or experiment_empty_for_verifier)
                and (verifier is None or verifier_resolved)
            ):
                try:
                    sampler_rc, sampler_terminal, _out, _err = (
                        self.backend.stop_sampler(sampler, prepared)
                    )
                    sampler_terminalized = True
                    if sampler_terminal.get("source_path_still_bound") is False:
                        record("sampler source pathname drifted after static seal")
                    if journal is not None and not sampler_terminal_write_attempted:
                        sampler_terminal_write_attempted = True
                        journal.write_bytes(
                            "sampler-terminal.json", sampler_terminal["raw"]
                        )
                except BaseException as sampler_exc:
                    record("sampler cleanup: " + exception_text(sampler_exc))
            if journal is not None and verifier_ran and not verification_write_attempted:
                try:
                    verification_write_attempted = True
                    journal.write_bytes("verification.json", verifier_stdout)
                except BaseException as verify_journal_exc:
                    record(
                        "verification journal: "
                        + exception_text(verify_journal_exc)
                    )
            terminal = {
                "claim_controller_started_monotonic_ns": controller_t0_ns,
                "controller_returncode": controller_rc,
                "controller_stderr_sha256": sha256_bytes(controller_stderr),
                "controller_stdout_sha256": sha256_bytes(controller_stdout),
                "containment_roster": sorted(set(containment_roster)),
                "errors": list(errors), "launcher_status": "fault",
                "guardian_deadline_monotonic_ns": deadline_ns,
                "guardian_fired": guardian.fired if guardian is not None else False,
                "outcome": classification,
                "run_argv_sha256": (
                    authority.run_argv_sha256 if authority is not None else None
                ),
                "sampler_returncode": sampler_rc,
                "sampler_terminal_receipt_sha256": (
                    sampler_terminal.get("receipt_sha256")
                    if sampler_terminalized else None
                ),
                "verifier_returncode": verifier_rc,
                "verifier_stderr_sha256": sha256_bytes(verifier_stderr),
                "verifier_stdout_sha256": sha256_bytes(verifier_stdout),
            }
        finally:
            # Quiescence, descriptor closure, and guardian reap all precede
            # terminal.json and COMPLETE.  COMPLETE can never promise cleanup
            # that is still pending in a finally block.
            try:
                quiesced = self.backend.quiesce() is True
                if not quiesced:
                    fail("runtime quiesce did not return an exact empty proof")
            except BaseException as cleanup_exc:
                record("runtime quiesce: " + exception_text(cleanup_exc))
                classification = "consumed_unclassifiable"
            lifecycle_resolved = (
                (controller is None or experiment_empty_for_verifier)
                and (verifier is None or verifier_resolved)
                and (sampler is None or sampler_terminalized)
            )
            if consumed and quiesced:
                try:
                    self.backend.finish_sampler_dir()
                    sampler_dir_closed = True
                except BaseException as sampler_dir_exc:
                    record(
                        "sampler directory close: "
                        + exception_text(sampler_dir_exc)
                    )
                    classification = "consumed_unclassifiable"
            if guardian is not None and quiesced and lifecycle_resolved:
                try:
                    if not guardian_outer_cancelled:
                        guardian.cancel()
                        guardian_outer_cancelled = True
                    guardian.finish()
                    guardian_finished = (
                        getattr(guardian, "finished", False) is True
                    )
                    if not guardian_finished:
                        fail("guardian finish lacked exact reap proof")
                    if guardian.error != "none" or guardian.backstop_fired:
                        record("guardian terminal: " + guardian.error)
                        classification = "consumed_unclassifiable"
                except BaseException as guardian_exc:
                    record("guardian terminal: " + exception_text(guardian_exc))
                    classification = "consumed_unclassifiable"
                try:
                    guardian_closed = guardian.close() is True
                    if not guardian_closed:
                        fail("guardian channel close lacked exact proof")
                except BaseException as guardian_close_exc:
                    record(
                        "guardian channel close: "
                        + exception_text(guardian_close_exc)
                    )
                    classification = "consumed_unclassifiable"
            if (
                quiesced and lifecycle_resolved
                and (guardian is None or guardian_finished)
            ):
                try:
                    runtime_closed = self.backend.close_runtime() is True
                    if not runtime_closed:
                        fail("runtime close did not return exact closure proof")
                except BaseException as close_exc:
                    record("runtime close: " + exception_text(close_exc))
                    classification = "consumed_unclassifiable"
            # Release the flock and every non-journal authority descriptor
            # before publishing completion.  The permanent launch namespace
            # itself still prevents another attempt after the lock is gone.
            try:
                authority_closed = self.backend.close_authority() is True
                if not authority_closed:
                    fail("authority close did not return exact closure proof")
            except BaseException as close_exc:
                record("authority close: " + exception_text(close_exc))
            if (
                consumed and journal is not None and runtime_closed
                and quiesced and lifecycle_resolved
                and sampler_dir_closed and authority_closed
                and (
                    guardian is None
                    or (guardian_finished and guardian_closed)
                )
            ):
                try:
                    terminal.update({
                        "errors": list(errors),
                        "launcher_status": "clean" if not errors else "fault",
                        "outcome": classification,
                    })
                    terminal_value = self_hashed(TERMINAL_SCHEMA, terminal)
                    journal.write_json("terminal.json", terminal_value)
                    journal.finish()
                    journal_complete = getattr(journal, "complete", False) is True
                    if not journal_complete:
                        fail("journal finish did not establish COMPLETE")
                    self._phase("complete")
                except BaseException as journal_exc:
                    record("terminal journal: " + exception_text(journal_exc))
                    classification = "consumed_unclassifiable"
                finally:
                    journal.close()
        if not consumed:
            raise LaunchError(errors[0] if errors else "read-only launch gate failed")
        if (
            not journal_complete or not quiesced
            or not lifecycle_resolved
            or not runtime_closed or not authority_closed
            or not sampler_dir_closed
            or (
                guardian is not None
                and (not guardian_finished or not guardian_closed)
            )
        ):
            raise LaunchError("attempt consumed but durable completion is incomplete")
        return RunResult(
            classification=classification,
            launcher_status="clean" if not errors else "fault",
            exit_code=(
                VERIFY_EXIT.get(classification, 1) if not errors else 1
            ),
            phases=list(self.phases),
        )


def selftest() -> int:
    commit = "1" * 40
    build = Path("/tmp/wh2-launch-selftest")
    args = parse_args([
        "--execute", "--build-dir", str(build),
        "--expected-source-commit", commit,
    ])
    if args.build_dir != build or args.expected_source_commit != commit:
        fail("launcher parser selftest differs")
    command = systemd_run_argv(build, commit)
    if command[-5:] != [
        "--execute", "--build-dir", str(build),
        "--expected-source-commit", commit,
    ]:
        fail("systemd command selftest differs")
    env_tail = [
        str(ENV_EXECUTABLE), "-i", "LANG=C.UTF-8", "LC_ALL=C.UTF-8",
        "PATH=/usr/bin:/bin", "TZ=UTC", str(PYTHON), "-I", "-B",
        str(INSTALLED_LAUNCHER),
    ]
    env_index = command.index(str(ENV_EXECUTABLE))
    if command[env_index:env_index + len(env_tail)] != env_tail:
        fail("systemd root environment selftest differs")
    for required in (
        "--property=Restart=no", "--property=ExitType=cgroup",
        "--property=RuntimeMaxSec=240s",
        "--property=Delegate=cpuset pids", "--property=AllowedCPUs=120-123",
        "--property=CPUAffinity=123", "--property=DevicePolicy=closed",
    ):
        if required not in command:
            fail("systemd property selftest differs")
    role = setpriv_exec_argv(["/sealed/controller", "--run-once"], sampler=False)
    if (
        "--clear-groups" not in role or "--no-new-privs" not in role
        or "--bounding-set=-all" not in role or "SIGTERM" not in role
        or role[-4:-2] != ["-I", "-B"]
    ):
        fail("unprivileged role selftest differs")
    compile(SEALED_MODULE_BOOTSTRAP, "<sealed-bootstrap-selftest>", "exec")
    compile(
        BUILD_CHILD_SECURITY_BOOTSTRAP,
        "<build-child-security-bootstrap-selftest>", "exec",
    )
    if parse_proc_cgroup(("0::/system.slice/" + UNIT_NAME + "\n").encode("ascii")) \
            != "/system.slice/" + UNIT_NAME:
        fail("cgroup parser selftest differs")
    print("WH2 direct systematic complement launcher selftest passed")
    return 0


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)
    if args.selftest:
        return selftest()
    if args.seal_build:
        receipt_sha256 = seal_build_authority(BuildSealConfig(
            source_dir=args.source_dir, build_dir=args.build_dir,
            expected_commit=args.expected_source_commit,
        ))
        print(canonical_bytes({
            "receipt_sha256": receipt_sha256,
            "schema": (
                "wirehair.wh2.direct-systematic-complement-build-result.v1"
            ),
            "status": "sealed",
        }).decode("ascii"))
        return 0
    backend = RealBackend()
    result = Supervisor(backend).execute(
        LaunchConfig(args.build_dir, args.expected_source_commit)
    )
    print(canonical_bytes({
        "launcher_status": result.launcher_status,
        "outcome": result.classification,
        "schema": "wirehair.wh2.direct-systematic-complement-launch-result.v1",
        "status": "complete",
    }).decode("ascii"))
    return result.exit_code


if __name__ == "__main__":
    try:
        exit_code = main(sys.argv[1:])
    except Exception as exc:
        print("WH2 launcher failed: " + exception_text(exc), file=sys.stderr)
        exit_code = 1
    raise SystemExit(exit_code)
