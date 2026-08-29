#!/usr/bin/python3
"""Root build sealer and one-shot supervisor for the V2 facade falsifier.

This is deliberately a standalone authority.  It does not import the older
direct-screen launcher: doing so would make that campaign-specific root module
part of this experiment's privileged code closure.

The three public modes are intentionally exact:

``--selftest``
    Pure and scratch-only checks.  It never opens hardware or a fixed campaign
    namespace.

``--seal-build``
    Build current and exact-parent role images in fresh directories, verify the
    complete Git trees and the effective compile/link graph, and publish the
    successful-only root build authority.

``--execute``
    Run one already-sealed attempt inside the delegated transient service.  A
    root journal, not the unprivileged controller directory, is the durable
    authority.  ``COMPLETE`` is created only after containment, retained replay,
    sampler authentication, copied-output verification, and descriptor closure.

The scientific controller remains responsible for its frozen timing protocol
and statistical adjudication.  This module is responsible only for attribution,
containment, and durable attempt accounting.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, field
import ctypes
import datetime
import errno
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import select
import shlex
import signal
import stat
import subprocess
import sys
import tempfile
import time
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


CAMPAIGN = "wh2-v2-facade-default-parent-falsifier-r0"
ATTEMPT_ID = CAMPAIGN + "-root-attempt-r0"
CURRENT_IMPLEMENTATION_COMMIT = "d6ddb35e046956174202584fb5a26c9a79679ea8"
PARENT_IMPLEMENTATION_COMMIT = "101584b7e5c30326b1791429221c331c82a00807"

UNIT_NAME = "wirehair-wh2-v2-facade-default-parent-falsifier-r0.service"
BUILD_UNIT_NAME = UNIT_NAME[:-8] + "-build.service"
INSTALLED_LAUNCHER = Path(
    "/usr/local/libexec/wirehair/Wh2V2FacadeTimingLaunch.py"
)
SOURCE_LAUNCHER_RELATIVE = "bench/Wh2V2FacadeTimingLaunch.py"
CONTROLLER_RELATIVE = "bench/Wh2V2FacadeTimingScreen.py"
WORKER_RELATIVE = "bench/Wh2V2FacadeTimingWorker.cpp"
ROLE_CMAKE_RELATIVE = "bench/public_facade_timing/CMakeLists.txt"
SAMPLER_RELATIVE = "bench/wirehair_expo_thermal_sampler.py"
HEALTH_ADAPTER_RELATIVE = "bench/Wh2DirectSystematicComplementScreen.py"

FROZEN_CONTROLLER_SHA256 = (
    "78470254363a067b4809aeed8f39ed7a212cbea8e84f0e3b5c754e79d8c6444f"
)
FROZEN_WORKER_SOURCE_SHA256 = (
    "fbafe52dbd677eee9c16bca5f8255f289f9e09bb62650a5372b12caa28d4d3c6"
)
FROZEN_ROLE_CMAKE_SHA256 = (
    "2a95010b94fa027f646106b753a5fc607b0982ee4eb9893ffe56ea6c1d97b17b"
)
FROZEN_HEALTH_ADAPTER_SHA256 = (
    "d6fa8cde8293796bdb97d05f1c7e9330a8ade4991f31644474562eee05f73464"
)

BUILD_AUTHORITY_DIR = Path("/var/lib/wirehair")
BUILD_AUTHORITY_PATH = BUILD_AUTHORITY_DIR / (CAMPAIGN + "-build-authority.json")
ATTEMPT_DIR = Path("/var/tmp") / (CAMPAIGN + ".root")
CONTROLLER_PARENT = Path("/var/tmp") / (CAMPAIGN + ".controller")
CONTROLLER_OUTPUT = CONTROLLER_PARENT / CAMPAIGN
SAMPLER_DIR = Path("/var/tmp") / (CAMPAIGN + ".sampler")
SAMPLER_CSV = SAMPLER_DIR / "thermal.csv"
SAMPLER_PID_FILE = SAMPLER_DIR / "sampler.pid"
SAMPLER_VALIDATION = SAMPLER_DIR / "validation.jsonl"
SAMPLER_RECEIPT = SAMPLER_DIR / "sampler-terminal.json"
LOCK_PATH = Path("/run/lock") / (CAMPAIGN + ".lock")

PYTHON = Path("/usr/bin/python3")
GIT = Path("/usr/bin/git")
CMAKE = Path("/usr/bin/cmake")
NINJA = Path("/usr/bin/ninja")
CTEST = Path("/usr/bin/ctest")
SETPRIV = Path("/usr/bin/setpriv")
ENV = Path("/usr/bin/env")
READELF = Path("/usr/bin/x86_64-linux-gnu-readelf")
LDD = Path("/usr/bin/ldd")
CXX = Path("/usr/bin/x86_64-linux-gnu-g++-13")
CC = Path("/usr/bin/x86_64-linux-gnu-gcc-13")
SYSTEMD_RUN = Path("/usr/bin/systemd-run")

FROZEN_TOOL_SHA256 = {
    Path("/usr/bin/python3.12"):
        "1643dacd9feaedc58f3cc581e4d22577dfe25c09b10282936186ccf0f2e61118",
    GIT: "2a8c18fbf43da9f692d75474c72bea9dfd796c260b0f3dfe456376abc3bbd668",
    CMAKE: "1c5227af4edd22d8d689def545e18ee458260c0fd579eba2187967f38817e638",
    NINJA: "5965527e09fe2b3787772aa4f711d6a36b393e7f2fcaa744a7a96c5a4ddf59cb",
    CTEST: "45c60ab2448411a4f49a0bdb6cf5800fedbf0293048e3c82710a34190c7d1ef4",
    SETPRIV: "96b083b79c32fd2f0c29657e88e20c7495839349fc64ad5d0503f32d26bf8733",
    ENV: "0aefff8f912fb75716c5d4de3b6acde93edbe8fa280fc8ee895c1226d3e373ef",
    READELF: "64c58e15274bbbb5153f31078e455e9e77ee5f51489e709bba5bb788ce9df2b0",
    LDD: "4f1d37e25f27535e3f02a5b7da63e1ce18d4982445db2c25fc8f985a3d395cc3",
    CXX: "1353e9bdd29a7295c7226bf6c63abccce056d8cac31f112e5cdbecc3f28c2769",
    CC: "1b99826121ae6682a634e5efe09bd3e3df58ce58e0b28f849114ab5b89139c26",
    SYSTEMD_RUN:
        "dbc8b988a849d5c9d7ef2de7068a6f107021bc6c11e0d7864c73f373eef726a7",
    Path("/usr/bin/x86_64-linux-gnu-ld.bfd"):
        "e9ceb054c12207970f2726dfc07e9a66b411602748628baf27399f02a9bbb31b",
    Path("/usr/bin/x86_64-linux-gnu-as"):
        "21aff249b692b5c31a44007491f922dcb49f41323e362c57d2ada3f52eddb7f0",
    Path("/usr/bin/x86_64-linux-gnu-ar"):
        "6452af2eea333b8c65e1adb92964fc8f97863ab003fa13f9d12bff5345cd7dbe",
    Path("/usr/bin/x86_64-linux-gnu-ranlib"):
        "60254978b8ee2c1b21b41d16a18a621dbbcc72e43eb2fa9f916818256c77e9ee",
    Path("/usr/lib/git-core/git"):
        "2a8c18fbf43da9f692d75474c72bea9dfd796c260b0f3dfe456376abc3bbd668",
    Path("/usr/bin/uname"):
        "8fd593fd9e1902eb67cfeba15642719b0ef6c44ba97febceac819e7e5c442abb",
    Path("/usr/bin/dash"):
        "86d31f6fb799e91fa21bad341484564510ca287703a16e9e46c53338776f4f42",
    Path("/usr/libexec/gcc/x86_64-linux-gnu/13/cc1"):
        "5d1679131184e2de4435b426eb264bf13472fe026db8e5c6bc97445814e8e2f4",
    Path("/usr/libexec/gcc/x86_64-linux-gnu/13/cc1plus"):
        "840b332fb62ec6f694ac77d91fe69ef7f80b0d69512ed89374af0ee7a506255d",
    Path("/usr/libexec/gcc/x86_64-linux-gnu/13/collect2"):
        "4d1f341ae5b763b513258ee2812422a45e063c30a2f1924a0cf63d3699f3a158",
    Path("/usr/libexec/gcc/x86_64-linux-gnu/13/liblto_plugin.so"):
        "14469ff8498fca3871c29dbf678c2330123207fffb8308151314a08395309ddf",
    Path("/usr/libexec/gcc/x86_64-linux-gnu/13/lto-wrapper"):
        "f6d39bcfb1c2798c704e275a4ed2c1fe6374ebfb6aa0c5f3f3b6c669e1ff2e13",
    Path("/usr/bin/bash"):
        "bc5945feb8bd26203ebfafea5ce1878bb2e32cb8fb50ab7ae395cfb1e1aaaef1",
    Path("/usr/lib32/ld-linux.so.2"):
        "8228542107069b0ffe686f0815581c995cb7bc17ba234f04159511a1f11b6178",
    Path("/usr/lib/x86_64-linux-gnu/ld-linux-x86-64.so.2"):
        "cd4df4f3c7b83673d61189bf2eaebd33ca4f2853ab9772b8a25e025ef99b1e81",
}

FROZEN_TOOL_NLINK = {
    Path("/usr/bin/python3.12"): 2,
    Path("/usr/bin/bash"): 2,
    Path("/usr/lib/x86_64-linux-gnu/ld-linux-x86-64.so.2"): 2,
}

FROZEN_TOOL_SYMLINKS = {
    Path("/usr/bin/python3"): "python3.12",
    Path("/usr/bin/readelf"): "x86_64-linux-gnu-readelf",
    Path("/usr/bin/ld.bfd"): "x86_64-linux-gnu-ld.bfd",
    Path("/usr/bin/as"): "x86_64-linux-gnu-as",
    Path("/usr/bin/ar"): "x86_64-linux-gnu-ar",
    Path("/usr/bin/ranlib"): "x86_64-linux-gnu-ranlib",
    Path("/bin/sh"): "dash",
    Path("/usr/bin/ld"): "x86_64-linux-gnu-ld",
    Path("/usr/bin/x86_64-linux-gnu-ld"): "x86_64-linux-gnu-ld.bfd",
    Path("/usr/lib/git-core/git-upload-pack"): "git",
    Path("/lib/ld-linux.so.2"): "../lib32/ld-linux.so.2",
    Path("/lib64/ld-linux-x86-64.so.2"):
        "../lib/x86_64-linux-gnu/ld-linux-x86-64.so.2",
}

CAMPAIGN_UID = 1000
CAMPAIGN_GID = 1000
I2C_GID = 113
SIBLING_CPU = 56
WORKER_CPU = 120
CONTROLLER_CPU = 121
SAMPLER_CPU = 122
AUTHORITY_CPU = 123
RUN_CPUS = (WORKER_CPU, CONTROLLER_CPU, SAMPLER_CPU, AUTHORITY_CPU)
HEALTH_CONTROLLER_INITIAL_AFFINITY = (
    WORKER_CPU, CONTROLLER_CPU, SAMPLER_CPU, AUTHORITY_CPU,
)

INTERNAL_DEADLINE_SECONDS = 840
EXTERNAL_DEADLINE_SECONDS = 900
# systemd's timer starts before the Python supervisor and therefore before the
# bounded pre-release interval.  It is an independent last-resort backstop,
# not the scientific deadline: the durable guardian remains authoritative at
# root T0 + 900 seconds.  Keep a separately named activation/stop cushion in
# addition to every launcher-owned interval.
SERVICE_ACTIVATION_STOP_MARGIN_SECONDS = 60
SERVICE_DEADLINE_SECONDS = 1020
# Git 2.43's local pack pipeline peaks above 64 tasks even when the build unit
# is pinned to two CPUs.  The exact sealed boundary measured 77 tasks; retain a
# bounded margin without constraining the independently frozen Ninja -j2 law.
BUILD_TASKS_MAX = 128
SAMPLER_READY_SECONDS = 12.0
CONTROLLER_ADMISSION_SECONDS = 12.0
SAMPLER_STOP_SECONDS = 8.0
REPLAY_SECONDS = 20.0
PRE_RELEASE_SECONDS = 10.0
POST_CONTROLLER_SECONDS = 45.0
BUILD_CONFIGURE_SECONDS = 120.0
BUILD_COMPILE_SECONDS = 900.0
BUILD_TEST_SECONDS = 180.0
BUILD_NO_WORK_SECONDS = 30.0

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
GIT_ENVIRONMENT = {
    "GIT_ASKPASS": "/bin/false",
    "GIT_CONFIG_GLOBAL": "/dev/null",
    "GIT_CONFIG_NOSYSTEM": "1",
    "GIT_NO_LAZY_FETCH": "1",
    "GIT_NO_REPLACE_OBJECTS": "1",
    "GIT_OPTIONAL_LOCKS": "0",
    "GIT_TERMINAL_PROMPT": "0",
    "LANG": "C",
    "LC_ALL": "C",
    "PATH": "/usr/bin:/bin",
    "SSH_ASKPASS": "/bin/false",
}

BUILD_AUTHORITY_SCHEMA = "wirehair.wh2.v2-facade-timing-build-authority.v1"
BUILD_RECEIPT_SCHEMA = "wirehair.wh2.v2-facade-timing-build-receipt.v1"
ATTEMPT_SCHEMA = "wirehair.wh2.v2-facade-timing-root-attempt.v1"
TERMINAL_SCHEMA = "wirehair.wh2.v2-facade-timing-root-terminal.v1"
COMPLETE_SCHEMA = "wirehair.wh2.v2-facade-timing-root-complete.v1"
CHILD_COMPLETE_SCHEMA = "wirehair.wh2.v2-facade-timing-screen.complete.v1"
CHILD_SUMMARY_SCHEMA = "wirehair.wh2.v2-facade-timing-screen.summary.v1"
CHILD_PROVENANCE_SCHEMA = "wirehair.wh2.v2-facade-timing-screen.provenance.v1"
REPLAY_SCHEMA = "wirehair.wh2.v2-facade-timing-screen.replay.v1"
SAMPLER_PRODUCER_SCHEMA = "wirehair.wh2.thermal_sampler.v2"
HEALTH_SAMPLER_ATTESTATION_SCHEMA = (
    "wirehair.wh2.direct-systematic-complement-sampler-attestation.v3"
)
PROCESS_OBSERVATION_LIMIT = (
    "trusted-controller-snapshot-monitor-v1: pre/post Git and worker-stub "
    "images may complete between 100ms observations; no runtime no-extra-process "
    "proof is claimed for those ephemeral phases"
)
SAMPLER_VALIDATION_STREAM_SCHEMA = "wirehair.wh2.thermal_validation_stream.v1"
SAMPLER_VALIDATION_SCHEMA = "wirehair.wh2.thermal_validation_sample.v1"
SAMPLER_COLUMNS = (
    "utc", "monotonic_s", "cpu_busy_pct", "cpu_avg_mhz", "cpu_tctl_c",
    "dimm_i2c1_50_c", "dimm_i2c1_51_c", "dimm_i2c1_52_c",
    "dimm_i2c1_53_c", "dimm_i2c2_50_c", "dimm_i2c2_51_c",
    "dimm_i2c2_52_c", "dimm_i2c2_53_c", "dimm_read_errors",
    "load1", "load5", "load15", "edac_ce", "edac_ue",
)
SAMPLER_DIMM_FIELDS = SAMPLER_COLUMNS[5:13]
SAMPLER_SAMPLING = {
    "dimm_attempts": 5, "dimm_retry_delay_s": 0.01, "interval_s": 1.0,
}
SAMPLER_THRESHOLDS = {
    "dimm_safety_c_inclusive": 90.0,
    "hot_confirm_samples": 3,
    "max_dimm_jump_c": 12.0,
    "max_dimm_rate_c_per_s": 6.0,
    "max_plausible_dimm_c_exclusive": 130.0,
    "min_plausible_dimm_c_exclusive": 0.0,
    "telemetry_fault_abort_samples": 8,
}
SAMPLER_VALIDATION_SAMPLE_KEYS = {
    "consecutive_fault_rows", "decision", "edac_ce_delta", "edac_ue_delta",
    "fault_count", "hot_sensors", "monotonic_s", "read_error_count",
    "sample_index", "schema", "sensors",
}
SAMPLER_VALIDATION_SENSOR_KEYS = {
    "attempt_errors", "hot", "hot_streak", "jump_c", "rate_c_per_s",
    "raw_c", "reason", "valid",
}
HEALTH_KEYS = {
    "admission_sibling_ticks", "child_reap_monotonic_ns",
    "child_start_monotonic_ns", "collection_failures", "controller_core",
    "controller_cpu", "controller_initial_affinity",
    "controller_singleton_affinity_end", "edac_policy", "evidence_status",
    "receipt_sha256", "sampler", "sampler_admission",
    "sampler_admission_receipt_sha256", "sampler_core", "sampler_cpu",
    "sampler_receipt_sha256", "sampler_terminal",
    "sampler_terminal_receipt_sha256", "schema", "sibling_non_idle_tick_cap",
    "sibling_tick_policy", "sibling_ticks", "target_core", "target_cpu",
    "target_threads", "thermal", "thermal_max_millic", "terminal_status",
    "violations",
}
HEALTH_SAMPLER_KEYS = {
    "cmdline_argv", "cmdline_sha256", "cpu", "csv_device", "csv_inode",
    "csv_path", "csv_bytes", "csv_sha256", "csv_stat", "evidence_parent",
    "environ_sha256", "environment", "environment_sha256",
    "executable_device", "executable_inode", "executable_path",
    "executable_sha256", "executable_stat", "pid", "pid_file",
    "process_affinity", "process_gid", "process_security",
    "process_start_ticks", "process_uid", "receipt_file", "schema",
    "script_device", "script_inode", "script_path", "script_sha256",
    "script_stat", "terminal_status", "window_end_monotonic_ns",
    "window_start_monotonic_ns", "validation_header_ascii",
    "validation_jsonl",
}
HEALTH_THERMAL_KEYS = {
    "cpu", "cpu_tctl_max_millic", "csv_device", "csv_inode", "csv_path",
    "dimm_max_millic", "dimm_read_errors", "edac_ce_max", "edac_ue_max",
    "invalid_sample_count", "parse_failures", "pid", "process_start_ticks",
    "sample_count", "script_path", "script_sha256", "terminal_status",
    "valid_sample_count", "validation_attempt_errors_total",
    "validation_device", "validation_failures", "validation_inode",
    "validation_jsonl_ascii", "validation_jsonl_bytes",
    "validation_jsonl_sha256", "validation_path", "validation_sample_count",
    "window_csv_ascii", "window_csv_bytes", "window_csv_sha256",
    "window_end_monotonic_ns", "window_start_monotonic_ns",
}
FROZEN_WORKER_STUB_CODE_SHA256 = (
    "839e4cfd2a64efa46674f10e9801e92ce02325207939b26a43dc90d19a4e8cf6"
)
FROZEN_INNER_GUARDIAN_CODE_SHA256 = (
    "cf4274f50cfd7691d0386d7dbbeb8d5fae97023541991974a8c916d7b2560836"
)

CGROUP_ROOT = Path("/sys/fs/cgroup")
CGROUP2_SUPER_MAGIC = 0x63677270
PR_SET_CHILD_SUBREAPER = 36
PR_GET_CHILD_SUBREAPER = 37
MEMORY_NODES = "0"

ROLE_SECURITY_BOOTSTRAP = r'''import json,os,sys
r=int(sys.argv[1]);g=int(sys.argv[2]);command=json.loads(sys.argv[3])
fields={}
for line in open("/proc/self/status","r",encoding="ascii"):
    if ":" in line:
        key,value=line.split(":",1)
        if key in ("CapInh","CapPrm","CapEff","CapBnd","CapAmb","NoNewPrivs"):
            fields[key]=value.strip()
value={"affinity":sorted(os.sched_getaffinity(0)),"caps":fields,
       "environment":dict(os.environ),"gid":list(os.getresgid()),
       "groups":os.getgroups(),"uid":list(os.getresuid())}
raw=(json.dumps(value,sort_keys=True,separators=(",",":"),ensure_ascii=True)+"\n").encode("ascii")
if os.write(r,raw)!=len(raw): raise SystemExit("short ready write")
os.close(r)
if os.read(g,1)!=b"R": raise SystemExit("release differs")
os.close(g)
os.execve(command[0],command,dict(os.environ))
'''

# Build-time controller runtime profile.  It loads the frozen screen under a
# non-__main__ module name, then asks the screen to load the exact sealed health
# adapter/modules.  Besides closing the interpreter map set, this executes the
# real adapter interface during build sealing, so an incompatible frozen pair
# cannot receive a successful build authority.
CONTROLLER_MAP_PROFILE_CODE = r'''import importlib.util,json,os,pathlib,sys
screen=pathlib.Path(sys.argv[1]); root=pathlib.Path(sys.argv[2])
name="_wh2_v2_facade_controller_profile"
spec=importlib.util.spec_from_file_location(name,screen)
if spec is None or spec.loader is None: raise SystemExit("screen spec")
module=importlib.util.module_from_spec(spec); sys.modules[name]=module
spec.loader.exec_module(module)
module.load_health_adapter(root,sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6])
paths=set()
for raw in open("/proc/self/maps","rb"):
    fields=raw.rstrip(b"\n").split(None,5)
    if len(fields)<5: raise SystemExit("maps shape")
    perms=fields[1]
    if len(perms)!=4 or perms[0:1] not in (b"r",b"-") or perms[1:2] not in (b"w",b"-") or perms[2:3] not in (b"x",b"-") or perms[3:4] not in (b"p",b"s"): raise SystemExit("maps perms")
    if perms[1:3]==b"wx": raise SystemExit("writable executable map")
    if len(fields)==5:
        if perms[2:3]==b"x": raise SystemExit("anonymous executable map")
        continue
    path=fields[5]
    if path.startswith(b"["):
        if perms[2:3]==b"x" and path not in (b"[vdso]",b"[vsyscall]"): raise SystemExit("pseudo executable map")
        continue
    if path.endswith(b" (deleted)") or not path.startswith(b"/"): raise SystemExit("map path")
    paths.add(os.path.realpath(os.fsdecode(path)))
print(json.dumps({"paths":sorted(paths),"schema":"wirehair.wh2.v2-facade-timing-controller-map-profile.v1"},sort_keys=True,separators=(",",":"),ensure_ascii=True))
'''

BUILD_SECURITY_BOOTSTRAP = r'''import json,os,sys
r=int(sys.argv[1]);g=int(sys.argv[2]);command=json.loads(sys.argv[3])
fields={}
for line in open("/proc/self/status","r",encoding="ascii"):
    if ":" in line:
        key,value=line.split(":",1)
        if key in ("CapInh","CapPrm","CapEff","CapBnd","CapAmb","NoNewPrivs"):
            fields[key]=value.strip()
value={"affinity":sorted(os.sched_getaffinity(0)),"caps":fields,
       "environment":dict(os.environ),"gid":list(os.getresgid()),
       "groups":os.getgroups(),"uid":list(os.getresuid())}
raw=(json.dumps(value,sort_keys=True,separators=(",",":"),ensure_ascii=True)+"\n").encode("ascii")
if os.write(r,raw)!=len(raw): raise SystemExit("short ready write")
os.close(r)
if os.read(g,1)!=b"R": raise SystemExit("release differs")
os.close(g)
os.execve(command[0],command,dict(os.environ))
'''

LOWER40 = re.compile(r"^[0-9a-f]{40}$")
LOWER64 = re.compile(r"^[0-9a-f]{64}$")
SAFE_RELATIVE = re.compile(r"^[A-Za-z0-9_.+@=/\-]+$")
MAX_TRACKED_FILES = 65536
MAX_TRACKED_FILE_BYTES = 1024 * 1024 * 1024
MAX_TRACKED_TOTAL_BYTES = 4 * 1024 * 1024 * 1024
MAX_DOCUMENT_BYTES = 16 * 1024 * 1024
MAX_RAW_BYTES = 32 * 1024 * 1024
MAX_JOURNAL_ARTIFACT_BYTES = MAX_RAW_BYTES
MAX_COMMAND_OUTPUT_BYTES = 64 * 1024 * 1024
MAX_STDERR_BYTES = 1024 * 1024

IMPLEMENTATION_SOURCES = (
    "wirehair.cpp",
    "gf256.cpp",
    "WirehairCodec.cpp",
    "WirehairTools.cpp",
    "codec/WirehairV2Codec.cpp",
    "codec/WirehairV2Peel.cpp",
    "codec/WirehairV2Plan.cpp",
    "codec/WirehairV2Policy.cpp",
    "codec/WirehairV2Precode.cpp",
    "codec/WirehairV2PrecodeDecode.cpp",
    "codec/WirehairV2PrecodeEncode.cpp",
    "codec/WirehairV2Profile.cpp",
    "codec/WirehairV2Seeds.cpp",
    "codec/WirehairV2Solve.cpp",
)
WORKER_SOURCES = (
    "bench/Wh2V2FacadeTimingWorker.cpp",
    "bench/Wh2FrozenTrace.cpp",
)
CONTROLLER_OUTPUT_NAMES = (
    "raw.jsonl",
    "summary.json",
    "provenance.json",
    "current.stderr",
    "parent.stderr",
    "COMPLETE",
)
HEALTH_SOURCE_PATHS = (
    "CMakeLists.txt",
    "bench/Wh2DirectSystematicComplementLaunch.py",
    "bench/Wh2DirectSystematicComplementScreen.cpp",
    "bench/Wh2DirectSystematicComplementScreen.py",
    "bench/Wh2FrozenTrace.cpp", "bench/Wh2FrozenTrace.h",
    "bench/Wh2NativeCodec.cpp", "bench/Wh2NativeCodec.h",
    "bench/Wh2NativePanel.cpp", "bench/Wh2NativePanel.h",
    "bench/wh2_benchmark_contract.py",
    "bench/wh2_benchmark_contract_v4.json",
    "bench/wh2_native_short_screen.py",
    "bench/wh2_run_native_short_screen.py",
    SAMPLER_RELATIVE,
    "cmake/Wh2DirectSystematicComplementSymbolAudit.cmake",
    "cmake/Wh2TimingPolicySymbolAudit.cmake",
)
HEALTH_MODULE_SOURCES = (
    ("wh2_benchmark_contract", "bench/wh2_benchmark_contract.py"),
    ("wh2_native_short_screen", "bench/wh2_native_short_screen.py"),
    ("wh2_run_native_short_screen", "bench/wh2_run_native_short_screen.py"),
)


class LaunchError(RuntimeError):
    """A fail-closed launcher or build-authority error."""


class AttemptConsumedError(LaunchError):
    """The permanent attempt namespace exists and must be terminalized."""

    def __init__(self, message: str, journal: Optional["AttemptJournal"]) -> None:
        super().__init__(message)
        self.journal = journal


def fail(message: str) -> None:
    raise LaunchError(message)


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
        allow_nan=False,
    ).encode("ascii")


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def utc_now() -> str:
    return datetime.datetime.now(datetime.timezone.utc).isoformat(
        timespec="microseconds"
    ).replace("+00:00", "Z")


def exception_text(exc: BaseException) -> str:
    value = "{}: {}".format(type(exc).__name__, exc)
    return value if len(value) <= 4096 else value[:4093] + "..."


def self_hashed(schema: str, values: Mapping[str, Any]) -> Dict[str, Any]:
    result = {"schema": schema, **dict(values), "receipt_sha256": None}
    result["receipt_sha256"] = sha256_bytes(canonical_bytes(result))
    return result


def parse_canonical_document(data: bytes, where: str,
                             maximum: int = MAX_DOCUMENT_BYTES) -> Dict[str, Any]:
    if not data or len(data) > maximum or not data.endswith(b"\n"):
        fail(where + " framing differs")

    def unique(pairs: Sequence[Tuple[str, Any]]) -> Dict[str, Any]:
        result: Dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                fail(where + " repeats key " + key)
            result[key] = value
        return result

    def reject_constant(value: str) -> None:
        fail(where + " contains non-finite number " + value)

    try:
        value = json.loads(
            data.decode("ascii"), object_pairs_hook=unique,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        fail(where + " is not canonical JSON: " + exception_text(exc))
    if type(value) is not dict or canonical_bytes(value) + b"\n" != data:
        fail(where + " is not one canonical object")
    return value


def exact_absolute(path: Path, where: str, *, must_exist: bool = False) -> Path:
    text = str(path)
    if not path.is_absolute() or os.path.normpath(text) != text or len(text) > 4096:
        fail(where + " is not a canonical absolute path")
    if must_exist:
        try:
            if path.resolve(strict=True) != path:
                fail(where + " traverses a symbolic link")
        except OSError as exc:
            fail(where + " is unavailable: " + exception_text(exc))
    return path


def parse_relative(raw: bytes, where: str) -> str:
    try:
        value = raw.decode("ascii")
    except UnicodeDecodeError as exc:
        fail(where + " path is not ASCII: " + exception_text(exc))
    components = value.split("/")
    if (
        not value or len(value) > 4096 or value.startswith("/")
        or SAFE_RELATIVE.fullmatch(value) is None
        or any(component in ("", ".", "..") for component in components)
    ):
        fail(where + " path differs")
    return value


def full_stat(value: os.stat_result) -> Tuple[int, ...]:
    return (
        value.st_dev,
        value.st_ino,
        value.st_mode,
        value.st_uid,
        value.st_gid,
        value.st_nlink,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )


def stat_receipt(value: os.stat_result) -> Dict[str, Any]:
    return {
        "ctime_ns": value.st_ctime_ns,
        "device": value.st_dev,
        "gid": value.st_gid,
        "inode": value.st_ino,
        "mode": stat.S_IMODE(value.st_mode),
        "mtime_ns": value.st_mtime_ns,
        "nlink": value.st_nlink,
        "size": value.st_size,
        "uid": value.st_uid,
    }


def hash_fd(fd: int, size: int) -> str:
    if type(size) is not int or size < 0 or size > MAX_TRACKED_FILE_BYTES:
        fail("file hash size bound differs")
    digest = hashlib.sha256()
    offset = 0
    while offset < size:
        block = os.pread(fd, min(1024 * 1024, size - offset), offset)
        if not block:
            fail("file hash read was short")
        digest.update(block)
        offset += len(block)
    return digest.hexdigest()


def git_blob_oid_fd(fd: int, size: int) -> str:
    digest = hashlib.sha1()
    digest.update(b"blob " + str(size).encode("ascii") + b"\0")
    offset = 0
    while offset < size:
        block = os.pread(fd, min(1024 * 1024, size - offset), offset)
        if not block:
            fail("Git blob read was short")
        digest.update(block)
        offset += len(block)
    return digest.hexdigest()


def hash_path(path: Path, maximum: int, where: str,
              *, nofollow: bool = True,
              expected_nlink: Optional[int] = 1) -> Tuple[str, os.stat_result]:
    flags = os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK
    if nofollow:
        flags |= getattr(os, "O_NOFOLLOW", 0)
    try:
        fd = os.open(str(path), flags)
    except OSError as exc:
        fail(where + " open failed: " + exception_text(exc))
    try:
        before = os.fstat(fd)
        if (
            not stat.S_ISREG(before.st_mode)
            or (expected_nlink is not None and before.st_nlink != expected_nlink)
            or before.st_size < 0 or before.st_size > maximum
        ):
            fail(where + " file policy differs")
        digest = hash_fd(fd, before.st_size)
        after = os.fstat(fd)
        named = os.stat(str(path), follow_symlinks=not nofollow)
        if full_stat(before) != full_stat(after) or full_stat(before) != full_stat(named):
            fail(where + " changed while hashing")
        return digest, before
    finally:
        os.close(fd)


def rename_noreplace_at(directory_fd: int, source: str, target: str) -> None:
    """Atomically rename one descriptor-relative name without replacement."""
    if (
        type(source) is not str or type(target) is not str
        or not source or not target or "/" in source or "/" in target
        or source in (".", "..") or target in (".", "..")
    ):
        fail("rename-noreplace name differs")
    libc = ctypes.CDLL(None, use_errno=True)
    function = getattr(libc, "renameat2", None)
    if function is None:
        fail("renameat2 is unavailable")
    function.argtypes = [
        ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p,
        ctypes.c_uint,
    ]
    function.restype = ctypes.c_int
    if function(
        directory_fd, os.fsencode(source), directory_fd, os.fsencode(target), 1,
    ) != 0:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error), target)


def same_inode(left: os.stat_result, right: os.stat_result) -> bool:
    return (left.st_dev, left.st_ino) == (right.st_dev, right.st_ino)


def open_fd_roster() -> List[int]:
    result: List[int] = []
    for token in os.listdir("/proc/self/fd"):
        if token.isdigit():
            fd = int(token)
            try:
                os.fstat(fd)
            except OSError:
                continue
            result.append(fd)
    return sorted(set(result))


def fd_authority_receipt(fd: int) -> Tuple[Any, ...]:
    info = os.fstat(fd)
    try:
        target = os.readlink("/proc/self/fd/{}".format(fd))
    except OSError as exc:
        fail("descriptor authority link differs: " + exception_text(exc))
    return (
        info.st_dev, info.st_ino, stat.S_IFMT(info.st_mode),
        fcntl.fcntl(fd, fcntl.F_GETFD), fcntl.fcntl(fd, fcntl.F_GETFL),
        target,
    )


def fd_authority_roster(fds: Iterable[int]) -> Dict[int, Tuple[Any, ...]]:
    return {fd: fd_authority_receipt(fd) for fd in sorted(set(fds))}


def verify_fd_authority_roster(
    expected: Mapping[int, Tuple[Any, ...]], where: str,
) -> None:
    if set(open_fd_roster()).issuperset(set(expected)) is not True:
        fail(where + " descriptor number roster is incomplete")
    for fd, receipt in expected.items():
        try:
            current = fd_authority_receipt(fd)
        except (OSError, LaunchError) as exc:
            fail(where + " descriptor {} is unavailable: {}".format(
                fd, exception_text(exc)))
        if current != receipt:
            fail(where + " descriptor {} identity changed".format(fd))


def close_untracked_inode_fds(info: os.stat_result,
                              retained: Iterable[int]) -> List[int]:
    """Close leaked descriptors for one private inode after async creation."""
    keep = set(retained)
    closed: List[int] = []
    for fd in open_fd_roster():
        if fd <= 2 or fd in keep:
            continue
        try:
            candidate = os.fstat(fd)
        except OSError:
            continue
        if same_inode(candidate, info):
            os.close(fd)
            closed.append(fd)
    return closed


def close_fd_delta(
    baseline: Iterable[int], retained: Iterable[int],
    closed_sink: Optional[List[int]] = None,
) -> List[int]:
    """Close descriptors created after a single-threaded admission snapshot.

    When supplied, ``closed_sink`` is the persistent evidence owner.  Record a
    descriptor number before the irreversible close so a return-handoff fault
    cannot erase the fact that the delta existed; exact final FD identity is
    still proved separately before COMPLETE.
    """
    original = set(baseline)
    keep = set(retained)
    closed = closed_sink if closed_sink is not None else []
    for fd in open_fd_roster():
        if fd <= 2 or fd in original or fd in keep:
            continue
        if fd not in closed:
            closed.append(fd)
        close_fd_once(fd, "descriptor-delta recovery")
    return closed


def close_fd_once(fd: int, where: str) -> None:
    """Close without ever retrying a numeric FD that may have been reused.

    Linux closes a descriptor even when close(2) reports EINTR.  The same
    return-before-assignment problem applies to injected Python exceptions:
    if the number is already EBADF, closure succeeded and retrying could close
    an unrelated subsequently-opened descriptor.  If it is still the exact
    pre-close object, ownership remains with the caller and the exception is
    propagated for a safe retry.
    """
    if fd < 0:
        return
    before = os.fstat(fd)
    try:
        os.close(fd)
        return
    except BaseException as exc:
        try:
            after = os.fstat(fd)
        except OSError as probe:
            if probe.errno == errno.EBADF:
                return
            raise LaunchError(
                where + " close probe failed: " + exception_text(probe)
            ) from exc
        if full_stat(after) != full_stat(before):
            fail(where + " descriptor number was reused during close")
        raise


def close_object_fd(owner: Any, member: str, where: str) -> None:
    states = getattr(owner, "_fd_close_states", None)
    if states is None:
        states = {}
        setattr(owner, "_fd_close_states", states)
    pending = states.get(member)
    current_member = getattr(owner, member)
    if pending is None:
        if current_member < 0:
            return
        info = os.fstat(current_member)
        pending = (current_member, info.st_dev, info.st_ino)
        states[member] = pending
    fd, device, inode = pending
    if current_member not in (fd, -1):
        fail(where + " descriptor ownership changed during close")
    try:
        info = os.fstat(fd)
    except OSError as exc:
        if exc.errno != errno.EBADF:
            raise
    else:
        if (info.st_dev, info.st_ino) != (device, inode):
            fail(where + " descriptor number was reused during close")
        try:
            os.close(fd)
        except BaseException:
            try:
                info = os.fstat(fd)
            except OSError as probe:
                if probe.errno != errno.EBADF:
                    raise
            else:
                if (info.st_dev, info.st_ino) != (device, inode):
                    fail(where + " descriptor number was reused during close")
                raise
    # The immutable identity record remains authoritative across exceptions
    # after close(2) returns.  Clear the public number only after probing the
    # exact descriptor, and leave the state resumable until the final pop.
    setattr(owner, member, -1)
    states.pop(member, None)


_DICT_FD_CLOSE_STATES: Dict[Tuple[int, str], Tuple[Dict[str, int], int, int, int]] = {}


def close_dict_fd(owner: Dict[str, int], key: str, where: str) -> None:
    state_key = (id(owner), key)
    pending = _DICT_FD_CLOSE_STATES.get(state_key)
    current = owner.get(key, -1)
    if pending is None:
        if current < 0:
            return
        info = os.fstat(current)
        pending = (owner, current, info.st_dev, info.st_ino)
        _DICT_FD_CLOSE_STATES[state_key] = pending
    pending_owner, fd, device, inode = pending
    if pending_owner is not owner or current not in (fd, -1):
        fail(where + " descriptor ownership changed during close")
    try:
        info = os.fstat(fd)
    except OSError as exc:
        if exc.errno != errno.EBADF:
            raise
    else:
        if (info.st_dev, info.st_ino) != (device, inode):
            fail(where + " descriptor number was reused during close")
        try:
            os.close(fd)
        except BaseException:
            try:
                info = os.fstat(fd)
            except OSError as probe:
                if probe.errno != errno.EBADF:
                    raise
            else:
                if (info.st_dev, info.st_ino) != (device, inode):
                    fail(where + " descriptor number was reused during close")
                raise
    owner.pop(key, None)
    _DICT_FD_CLOSE_STATES.pop(state_key, None)


def open_registered(roster: List[int], path: Any, flags: int,
                    mode: Optional[int] = None, **kwargs: Any) -> None:
    """Open into an owner roster without a return-before-assignment FD gap."""
    baseline = open_fd_roster()
    fd = -1
    try:
        if mode is None:
            fd = os.open(path, flags, **kwargs)
        else:
            fd = os.open(path, flags, mode, **kwargs)
        roster.append(fd)
    except BaseException:
        retained = list(roster)
        if fd >= 0:
            retained.append(fd)
        try:
            close_fd_delta(baseline, retained)
        finally:
            if fd >= 0:
                try:
                    os.close(fd)
                except OSError:
                    pass
            if roster and roster[-1] == fd:
                roster.pop()
        raise


def require_root_owned_nonwritable_chain(path: Path) -> None:
    exact_absolute(path, "root authority path", must_exist=True)
    current = Path("/")
    for index, part in enumerate(path.parts[1:]):
        current = current / part
        info = os.stat(str(current), follow_symlinks=False)
        leaf = index == len(path.parts[1:]) - 1
        if (
            info.st_uid != 0 or info.st_gid != 0
            or stat.S_IMODE(info.st_mode) & 0o022
            or (leaf and not stat.S_ISREG(info.st_mode))
            or (not leaf and not stat.S_ISDIR(info.st_mode))
        ):
            fail("root authority path policy differs: " + str(current))


def direct_child_pids() -> List[int]:
    raw = read_bounded_path(
        Path("/proc/self/task/{}/children".format(os.getpid())),
        1024 * 1024, "direct child roster",
    )
    try:
        values = [int(token) for token in raw.split()]
    except ValueError as exc:
        fail("direct child roster differs: " + exception_text(exc))
    if any(value <= 0 for value in values) or len(values) != len(set(values)):
        fail("direct child roster identity differs")
    return sorted(values)


@dataclass
class BoundedCommandOwnership:
    """Persistent authority from before acquisition through quiescence."""

    label: str
    fd_baseline: Tuple[int, ...]
    child_baseline: Tuple[int, ...]
    service_baseline: Tuple[int, ...]
    process_owner: OwnershipSlot = field(
        default_factory=lambda: OwnershipSlot("bounded command process"),
    )
    fds: Dict[str, int] = field(default_factory=dict)
    may_have_released: bool = False
    observed_descendants: List[int] = field(default_factory=list)

    @classmethod
    def create(cls, label: str) -> "BoundedCommandOwnership":
        # Descendants that double-fork or create a new session must return to
        # this supervisor, rather than escaping process-group-only cleanup.
        become_subreaper()
        return cls(
            label=label, fd_baseline=tuple(open_fd_roster()),
            child_baseline=tuple(direct_child_pids()),
            service_baseline=tuple(current_service_pids()),
            process_owner=OwnershipSlot(label + " process"),
        )

    @property
    def process(self) -> Optional[subprocess.Popen]:
        return self.process_owner.value

    def spawn(self, *arguments: Any, **keywords: Any) -> subprocess.Popen:
        # Adopt the Python object before __init__ can fork.  Thus an
        # initializer that creates a child and then raises still leaves its
        # pid/_child_created state reachable by quiesce().
        process: Optional[subprocess.Popen] = None
        try:
            process = subprocess.Popen.__new__(subprocess.Popen)
            # Keep Popen.__del__ safe even if an injected exception lands
            # before __init__ establishes its own initial state.
            process._child_created = False
            process.returncode = None
            self.process_owner.adopt(process)
            subprocess.Popen.__init__(process, *arguments, **keywords)
            return process
        except BaseException:
            if process is not None:
                if not hasattr(process, "_child_created"):
                    process._child_created = False
                if not hasattr(process, "returncode"):
                    process.returncode = None
            raise

    @staticmethod
    def _signal_pid(pid: int) -> bool:
        baseline = open_fd_roster()
        try:
            try:
                start_ticks = process_start_ticks(pid)
                pidfd = os.pidfd_open(pid, 0)
                if process_start_ticks(pid) != start_ticks:
                    fail("bounded descendant identity changed before kill")
            except (FileNotFoundError, ProcessLookupError):
                return False
            # Authority-roster failures are not evidence that the target
            # disappeared.  Keep them outside the absence-only catch so
            # quiescence fails closed if either held namespace cannot be read.
            if (
                pid not in current_service_pids()
                and pid not in direct_child_pids()
            ):
                return False
            try:
                if process_start_ticks(pid) != start_ticks:
                    fail("bounded descendant identity changed during kill admission")
            except (FileNotFoundError, ProcessLookupError):
                return False
            readable, _, exceptional = select.select(
                [pidfd], [], [pidfd], 0.0,
            )
            if exceptional:
                fail("bounded descendant pidfd readiness differs")
            if readable:
                return False
            try:
                signal.pidfd_send_signal(pidfd, signal.SIGKILL, None, 0)
            except (FileNotFoundError, ProcessLookupError):
                # The identity-bound pidfd was nonterminal at the readiness
                # sample.  A concurrent natural exit does not erase that
                # proven post-command liveness.
                return True
            return True
        finally:
            close_fd_delta(baseline, ())

    def quiesce(self) -> List[int]:
        failures: List[str] = []
        observed: set[int] = set(self.observed_descendants)
        process = self.process
        pid = getattr(process, "pid", -1) if process is not None else -1
        isolated_service = self.service_baseline == (os.getpid(),)
        deadline = time.monotonic() + 5.0
        while True:
            try:
                service_now = (
                    set(current_service_pids())
                    if isolated_service else set(self.service_baseline)
                )
                child_now = set(direct_child_pids())
            except BaseException as exc:
                failures.append("descendant roster: " + exception_text(exc))
                break
            new_service = service_now - set(self.service_baseline)
            new_children = child_now - set(self.child_baseline)
            candidates = sorted(new_service | new_children)
            live_candidates: List[int] = []
            for candidate in candidates:
                try:
                    waited, status_value = os.waitpid(candidate, os.WNOHANG)
                    if waited == candidate:
                        # A subreaper adopts both live descendants and
                        # already-exited zombies.  Reap the latter as ordinary
                        # command closure; they cannot execute after the
                        # command has returned and are not a live escape.
                        if candidate == pid and process is not None:
                            process.returncode = wait_status_code(status_value)
                        continue
                    if waited != 0:
                        failures.append(
                            "descendant {} initial reap identity differs".format(
                                candidate))
                    live_candidates.append(candidate)
                except ChildProcessError:
                    # A service-cgroup descendant may not yet have been
                    # reparented to this subreaper.  Let the identity-bound
                    # pidfd signal distinguish a still-live member from a PID
                    # that disappeared after the roster snapshot.
                    try:
                        if self._signal_pid(candidate):
                            observed.add(candidate)
                    except BaseException as exc:
                        failures.append(
                            "descendant {} signal: {}".format(
                                candidate, exception_text(exc)))
                except BaseException as exc:
                    failures.append(
                        "descendant {} initial reap: {}".format(
                            candidate, exception_text(exc)))
                    live_candidates.append(candidate)
            observed.update(live_candidates)
            for candidate in live_candidates:
                try:
                    self._signal_pid(candidate)
                except BaseException as exc:
                    failures.append(
                        "descendant {} signal: {}".format(
                            candidate, exception_text(exc)))
            for candidate in live_candidates:
                try:
                    waited, status_value = os.waitpid(candidate, os.WNOHANG)
                    if waited not in (0, candidate):
                        failures.append(
                            "descendant {} reap identity differs".format(
                                candidate))
                    elif (
                        waited == candidate and candidate == pid
                        and process is not None
                    ):
                        process.returncode = wait_status_code(status_value)
                except ChildProcessError:
                    if candidate == pid and process is not None:
                        # A prior interrupted Popen wait may already have
                        # reaped the exact child without publishing its
                        # Python returncode.  The process is terminal; retain
                        # a conservative killed status for destructor safety.
                        process.returncode = -signal.SIGKILL
                    pass
                except BaseException as exc:
                    failures.append(
                        "descendant {} reap: {}".format(
                            candidate, exception_text(exc)))
            if not candidates:
                break
            if time.monotonic() >= deadline:
                failures.append("descendant quiescence deadline expired")
                break
            time.sleep(0.01)
        if type(pid) is int and pid > 0:
            try:
                os.killpg(pid, 0)
            except ProcessLookupError:
                pass
            except BaseException as exc:
                failures.append("process group probe: " + exception_text(exc))
            else:
                failures.append("process group remained populated")
        for name in list(self.fds):
            try:
                close_dict_fd(self.fds, name, self.label + " " + name)
            except BaseException as exc:
                failures.append(name + " closure: " + exception_text(exc))
        try:
            close_fd_delta(self.fd_baseline, ())
        except BaseException as exc:
            failures.append("descriptor delta: " + exception_text(exc))
        try:
            final_children = set(direct_child_pids())
            if final_children - set(self.child_baseline):
                failures.append("adopted descendants remain after cleanup")
            if isolated_service:
                final_service = set(current_service_pids())
                if final_service != {os.getpid()}:
                    failures.append("isolated service baseline changed")
        except BaseException as exc:
            failures.append("final descendant proof: " + exception_text(exc))
        self.observed_descendants = sorted(observed)
        if failures:
            fail(self.label + " quiescence differs: " + " | ".join(failures))
        return list(self.observed_descendants)


def _bounded_command(
    argv: Sequence[str], *, cwd: Optional[Path], environment: Mapping[str, str],
    timeout: float, stdout_limit: int = MAX_COMMAND_OUTPUT_BYTES,
    stderr_limit: int = MAX_STDERR_BYTES,
) -> Tuple[bytes, bytes]:
    if (
        not argv or not Path(argv[0]).is_absolute()
        or any(type(item) is not str or not item or len(item) > 4096 for item in argv)
        or timeout <= 0.0
    ):
        fail("bounded command contract differs")
    ownership = BoundedCommandOwnership.create("bounded command " + argv[0])
    primary: Optional[BaseException] = None
    code = -1
    output = error = b""
    try:
        with tempfile.TemporaryFile() as stdout_file, tempfile.TemporaryFile() as stderr_file:
            process = ownership.spawn(
                list(argv), cwd=None if cwd is None else str(cwd),
                env=dict(environment), stdin=subprocess.DEVNULL,
                stdout=stdout_file, stderr=stderr_file, close_fds=True,
                start_new_session=True,
            )
            ownership.may_have_released = True
            try:
                code = process.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                fail("command deadline expired: " + argv[0])
            out_size = os.fstat(stdout_file.fileno()).st_size
            err_size = os.fstat(stderr_file.fileno()).st_size
            if out_size > stdout_limit or err_size > stderr_limit:
                fail("command output bound exceeded: " + argv[0])
            stdout_file.seek(0)
            stderr_file.seek(0)
            output = stdout_file.read()
            error = stderr_file.read()
    except BaseException as exc:
        primary = exc
    cleanup: Optional[BaseException] = None
    try:
        descendants = ownership.quiesce()
        if descendants and primary is None:
            primary = LaunchError(
                "command left one or more descendants: " + repr(descendants))
    except BaseException as exc:
        cleanup = exc
    if primary is not None:
        if cleanup is not None:
            raise LaunchError(
                exception_text(primary) + " | " + exception_text(cleanup)
            ) from primary
        raise primary
    if cleanup is not None:
        raise cleanup
    if code != 0:
        fail("command failed rc={} argv={} stderr_sha256={}".format(
            code, argv[0], sha256_bytes(error)))
    return output, error


def git_command(root: Path, arguments: Sequence[str], maximum: int) -> bytes:
    exact_absolute(root, "Git root", must_exist=True)
    argv = [
        str(GIT), "-c", "core.fsmonitor=false",
        "-c", "core.untrackedCache=false",
        "-c", "core.filemode=false",
        "-c", "safe.directory=" + str(root),
        "-C", str(root),
    ] + list(arguments)
    output, error = _bounded_command(
        argv, cwd=None, environment=GIT_ENVIRONMENT, timeout=20.0,
        stdout_limit=maximum,
    )
    if error:
        fail("Git command wrote stderr: " + sha256_bytes(error))
    return output


@dataclass(frozen=True)
class GitEntry:
    mode: str
    oid: str
    relative: str


@dataclass(frozen=True)
class GitTreeReceipt:
    root: Path
    commit: str
    tree_oid: str
    tree_bytes: int
    tree_listing_sha256: str
    manifest_sha256: str
    entries: Tuple[Dict[str, Any], ...]


def parse_git_tree(data: bytes, where: str, *, index: bool) -> List[GitEntry]:
    if not data or not data.endswith(b"\0"):
        fail(where + " framing differs")
    records = data[:-1].split(b"\0")
    if not records or len(records) > MAX_TRACKED_FILES:
        fail(where + " count differs")
    result: List[GitEntry] = []
    paths: set = set()
    for record in records:
        try:
            metadata, raw_path = record.split(b"\t", 1)
            fields = metadata.split(b" ")
            if len(fields) != 3:
                fail(where + " metadata differs")
            if index:
                raw_mode, raw_oid, stage = fields
                if stage != b"0":
                    fail(where + " contains nonzero stage")
            else:
                raw_mode, kind, raw_oid = fields
                if kind != b"blob":
                    fail(where + " contains a non-blob entry")
            mode = raw_mode.decode("ascii")
            oid = raw_oid.decode("ascii")
        except (UnicodeDecodeError, ValueError) as exc:
            fail(where + " row malformed: " + exception_text(exc))
        relative = parse_relative(raw_path, where)
        if mode not in ("100644", "100755") or LOWER40.fullmatch(oid) is None:
            fail(where + " blob identity differs")
        if relative in paths:
            fail(where + " repeats a path")
        paths.add(relative)
        result.append(GitEntry(mode, oid, relative))
    return result


def verify_git_tree(
    root: Path, commit: str, *, required_uid: Optional[int] = None,
    required_gid: Optional[int] = None, sealed: bool = False,
) -> GitTreeReceipt:
    """Verify HEAD, index, ignored/untracked absence, and every tracked blob.

    Git status is used only as an additional diagnostic.  The direct descriptor
    hash of every HEAD blob is authoritative and therefore cannot be bypassed by
    index stat-cache or fsmonitor metadata.
    """
    exact_absolute(root, "Git source root", must_exist=True)
    if LOWER40.fullmatch(commit) is None:
        fail("expected Git commit differs")
    top = git_command(root, ("rev-parse", "--show-toplevel"), 8192)
    if top != (str(root) + "\n").encode("ascii"):
        fail("Git top-level differs")
    head = git_command(root, ("rev-parse", "--verify", "HEAD^{commit}"), 8192)
    if head != (commit + "\n").encode("ascii"):
        fail("Git HEAD differs")
    if git_command(root, ("rev-parse", "--abbrev-ref", "HEAD"), 8192) != b"HEAD\n":
        fail("Git source is not detached")
    if git_command(
        root, ("status", "--porcelain=v1", "-z", "--untracked-files=all"),
        32 * 1024 * 1024,
    ):
        fail("Git worktree is dirty")
    if git_command(root, ("ls-files", "--others", "--exclude-standard", "-z"),
                   32 * 1024 * 1024):
        fail("Git worktree contains untracked files")
    if git_command(
        root, ("ls-files", "--others", "--ignored", "--exclude-standard", "-z"),
        32 * 1024 * 1024,
    ):
        fail("Git worktree contains ignored files")
    tree_oid_raw = git_command(
        root, ("rev-parse", "--verify", commit + "^{tree}"), 8192)
    try:
        tree_oid = tree_oid_raw.rstrip(b"\n").decode("ascii")
    except UnicodeDecodeError as exc:
        fail("Git tree OID malformed: " + exception_text(exc))
    if LOWER40.fullmatch(tree_oid) is None or tree_oid_raw != (
        tree_oid + "\n"
    ).encode("ascii"):
        fail("Git tree OID differs")
    raw_tree = git_command(
        root, ("ls-tree", "--full-tree", "-r", "-z", commit),
        32 * 1024 * 1024,
    )
    raw_index = git_command(
        root, ("ls-files", "--stage", "-z"), 32 * 1024 * 1024,
    )
    tree = parse_git_tree(raw_tree, "Git HEAD tree", index=False)
    index = parse_git_tree(raw_index, "Git stage-0 index", index=True)
    if tree != index:
        fail("Git index does not exactly equal HEAD")
    total = 0
    manifest = hashlib.sha256()
    entries: List[Dict[str, Any]] = []
    for item in tree:
        path = root / item.relative
        if path.resolve(strict=True) != path:
            fail("tracked path traverses a link: " + item.relative)
        fd = os.open(
            str(path), os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK
            | getattr(os, "O_NOFOLLOW", 0),
        )
        try:
            before = os.fstat(fd)
            named = os.stat(str(path), follow_symlinks=False)
            live_mode = stat.S_IMODE(before.st_mode)
            mode_matches_tree = (
                live_mode == 0o444 if sealed else
                (
                    live_mode & 0o7000 == 0 and
                    bool(live_mode & 0o111) == (item.mode == "100755")
                )
            )
            if (
                not stat.S_ISREG(before.st_mode) or before.st_nlink != 1
                or before.st_size < 0 or before.st_size > MAX_TRACKED_FILE_BYTES
                or not mode_matches_tree
                or (required_uid is not None and before.st_uid != required_uid)
                or (required_gid is not None and before.st_gid != required_gid)
                or full_stat(before) != full_stat(named)
            ):
                fail("tracked file policy differs: " + item.relative)
            total += before.st_size
            if total > MAX_TRACKED_TOTAL_BYTES:
                fail("tracked source byte bound differs")
            oid = git_blob_oid_fd(fd, before.st_size)
            digest = hash_fd(fd, before.st_size)
            after = os.fstat(fd)
            named_after = os.stat(str(path), follow_symlinks=False)
            if (
                oid != item.oid or full_stat(after) != full_stat(before)
                or full_stat(named_after) != full_stat(before)
            ):
                fail("tracked bytes differ from HEAD: " + item.relative)
            receipt = {
                "bytes": before.st_size,
                "git_blob_oid": oid,
                "git_mode": item.mode,
                "repository_relative_path": item.relative,
                "sha256": digest,
                "stat": stat_receipt(before),
            }
            entries.append(receipt)
            manifest.update(canonical_bytes(receipt) + b"\n")
        finally:
            os.close(fd)
    return GitTreeReceipt(
        root=root,
        commit=commit,
        tree_oid=tree_oid,
        tree_bytes=len(raw_tree),
        tree_listing_sha256=sha256_bytes(raw_tree),
        manifest_sha256=manifest.hexdigest(),
        entries=tuple(entries),
    )


def systemd_run_argv(expected_harness_commit: str) -> List[str]:
    if LOWER40.fullmatch(expected_harness_commit) is None:
        fail("systemd harness commit differs")
    properties = (
        "User=0", "Group=0", "SupplementaryGroups=", "Restart=no",
        "ExitType=cgroup", "KillMode=control-group", "SendSIGKILL=yes",
        "TimeoutStopSec=1s",
        "RuntimeMaxSec={}s".format(SERVICE_DEADLINE_SECONDS),
        "Delegate=cpuset pids", "AllowedCPUs=120-123",
        "AllowedMemoryNodes=0", "CPUAffinity=123", "TasksMax=40",
        "LimitCORE=0", "UMask=0077", "PrivateTmp=no",
        "PrivateDevices=no", "ProtectControlGroups=no", "DevicePolicy=closed",
        "DeviceAllow=/dev/i2c-1 rw", "DeviceAllow=/dev/i2c-2 rw",
    )
    argv = [
        "/usr/bin/systemd-run", "--quiet", "--wait", "--pipe", "--collect",
        "--service-type=exec", "--expand-environment=no", "--unit=" + UNIT_NAME,
    ]
    argv.extend("--property=" + item for item in properties)
    argv.extend((
        str(ENV), "-i",
        "LANG=" + SEALED_ENVIRONMENT["LANG"],
        "LC_ALL=" + SEALED_ENVIRONMENT["LC_ALL"],
        "PATH=" + SEALED_ENVIRONMENT["PATH"],
        "TZ=" + SEALED_ENVIRONMENT["TZ"],
        str(PYTHON), "-I", "-B", str(INSTALLED_LAUNCHER),
        "--execute", "--expected-harness-commit", expected_harness_commit,
    ))
    return argv


class AttemptJournal:
    """Descriptor-relative, O_EXCL-only root journal with COMPLETE last."""

    def __init__(self, path: Path, parent_fd: int, directory_fd: int,
                 *, expected_uid: int, expected_gid: int,
                 fault_injector: Optional[Any] = None,
                 required_prefix: Sequence[str] = ("ATTEMPT",)) -> None:
        self.path = path
        self.parent_fd = parent_fd
        self.directory_fd = directory_fd
        self.parent_fd_authority = (
            fd_authority_receipt(parent_fd) if parent_fd >= 0 else None
        )
        self.directory_fd_authority = (
            fd_authority_receipt(directory_fd) if directory_fd >= 0 else None
        )
        self.binding: Optional[Tuple[int, int]] = None
        if directory_fd >= 0:
            info = os.fstat(directory_fd)
            self.binding = (info.st_dev, info.st_ino)
        self.expected_uid = expected_uid
        self.expected_gid = expected_gid
        self.names: List[str] = []
        self.stream_fds: Dict[str, int] = {}
        # Every short-lived journal descriptor is owned here before the open
        # helper returns.  This roster deliberately outlives the individual
        # publication call: an asynchronous exception at a Python handoff
        # boundary must leave enough identity to close (or diagnose reuse of)
        # the exact descriptor before COMPLETE can be named.
        self.transient_fds: List[int] = []
        self.transient_fd_stats: Dict[int, Tuple[int, int]] = {}
        self.complete = False
        self.fault_injector = fault_injector
        self.required_prefix = tuple(required_prefix)
        if (
            not self.required_prefix
            or self.required_prefix[0] != "ATTEMPT"
            or len(set(self.required_prefix)) != len(self.required_prefix)
            or any(not name or "/" in name or name in (".", "..")
                   for name in self.required_prefix)
        ):
            fail("attempt mandatory prefix differs")

    def _fault(self, point: str, name: str) -> None:
        if self.fault_injector is not None:
            self.fault_injector(point, name)

    def _open_transient(self, path: Any, flags: int,
                        mode: Optional[int] = None, **kwargs: Any) -> int:
        open_registered(self.transient_fds, path, flags, mode, **kwargs)
        fd = self.transient_fds[-1]
        info = os.fstat(fd)
        self.transient_fd_stats[fd] = (info.st_dev, info.st_ino)
        return fd

    def _close_transient(self, fd: int) -> None:
        """Close one persistent journal FD without numeric-reuse hazards."""
        if fd not in self.transient_fds:
            return
        expected = self.transient_fd_stats.get(fd)
        try:
            live = os.fstat(fd)
        except OSError as exc:
            if exc.errno != errno.EBADF:
                raise
            # A real close followed by an asynchronous exception can leave
            # the persistent roster stale.  EBADF is the only safe recovery.
            self.transient_fds.remove(fd)
            self.transient_fd_stats.pop(fd, None)
            return
        current = (live.st_dev, live.st_ino)
        if expected is None:
            self.transient_fd_stats[fd] = current
            expected = current
        if current != expected:
            # Never close a reused numeric descriptor.  Forget the stale
            # journal ownership, but fail closed so COMPLETE is forbidden.
            self.transient_fds.remove(fd)
            self.transient_fd_stats.pop(fd, None)
            fail("journal transient descriptor number was reused")
        close_fd_once(fd, "journal transient descriptor")
        self.transient_fds.remove(fd)
        self.transient_fd_stats.pop(fd, None)

    def _drain_transients(self) -> None:
        failures: List[str] = []
        for fd in list(reversed(self.transient_fds)):
            try:
                self._close_transient(fd)
            except BaseException as exc:
                failures.append(exception_text(exc))
        if failures:
            fail("journal transient descriptor closure differs: "
                 + " | ".join(failures))

    def _transient_matches_name(self, fd: int, name: str) -> bool:
        try:
            held = os.fstat(fd)
            named = os.stat(
                name, dir_fd=self.directory_fd, follow_symlinks=False,
            )
        except (FileNotFoundError, OSError):
            return False
        return same_inode(held, named)

    @staticmethod
    def _private_name(kind: str) -> str:
        return ".{}.{}.{}.{}".format(
            kind, os.getpid(), time.monotonic_ns(), os.urandom(8).hex(),
        )

    def _ensure_directory(self) -> None:
        if self.directory_fd >= 0:
            return
        baseline = open_fd_roster()
        owned: List[int] = []
        old_binding = self.binding
        old_authority = self.directory_fd_authority
        fd = -1
        try:
            open_registered(
                owned, self.path.name,
                os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC
                | getattr(os, "O_NOFOLLOW", 0),
                dir_fd=self.parent_fd,
            )
            fd = owned[0]
            info = os.fstat(fd)
            named = os.stat(
                self.path.name, dir_fd=self.parent_fd, follow_symlinks=False,
            )
            binding = (info.st_dev, info.st_ino)
            if (
                not stat.S_ISDIR(info.st_mode)
                or binding != (named.st_dev, named.st_ino)
                or stat.S_IMODE(info.st_mode) != 0o700
                or info.st_uid != self.expected_uid
                or info.st_gid != self.expected_gid
                or (self.binding is not None and binding != self.binding)
            ):
                fail("attempt directory policy differs")
            self.binding = binding
            self.directory_fd_authority = fd_authority_receipt(fd)
            self.directory_fd = fd
        except BaseException as original:
            if fd >= 0 and self.directory_fd == fd:
                # The fully validated descriptor crossed its persistent
                # ownership handoff; a retry can use it directly.
                raise
            self.binding = old_binding
            self.directory_fd_authority = old_authority
            cleanup_failure: Optional[BaseException] = None
            try:
                close_fd_delta(baseline, ())
            except BaseException as exc:
                cleanup_failure = exc
            if cleanup_failure is not None:
                raise LaunchError(
                    exception_text(original)
                    + " | attempt directory reopen cleanup: "
                    + exception_text(cleanup_failure)
                ) from original
            raise

    @classmethod
    def reserve(
        cls, path: Path, attempt: Mapping[str, Any], *,
        expected_uid: int = 0, expected_gid: int = 0,
        test_parent: bool = False, fault_injector: Optional[Any] = None,
        owner: Optional[Any] = None,
        required_prefix: Sequence[str] = ("ATTEMPT",),
    ) -> "AttemptJournal":
        exact_absolute(path, "attempt journal path")
        parent = path.parent
        parent_handles: List[int] = []
        reservation_baseline = open_fd_roster()
        parent_fd = -1
        journal: Optional[AttemptJournal] = None
        temporary = cls._private_name("attempt")
        temporary_exists = False
        consumed = False
        try:
            open_registered(
                parent_handles, str(parent),
                os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC
                | getattr(os, "O_NOFOLLOW", 0),
            )
            parent_fd = parent_handles[0]
            parent_info = os.fstat(parent_fd)
            if not stat.S_ISDIR(parent_info.st_mode):
                fail("attempt parent is not a directory")
            if not test_parent and (
                parent_info.st_uid != 0 or parent_info.st_gid != 0
                or (
                    parent == Path("/var/tmp")
                    and stat.S_IMODE(parent_info.st_mode) != 0o1777
                )
            ):
                fail("attempt parent policy differs")
            try:
                os.stat(
                    path.name, dir_fd=parent_fd, follow_symlinks=False,
                )
            except FileNotFoundError:
                pass
            else:
                fail("attempt namespace already exists at reservation")
            try:
                os.stat(temporary, dir_fd=parent_fd, follow_symlinks=False)
            except FileNotFoundError:
                pass
            else:
                fail("attempt private reservation name already exists")
            # Creation first consumes a private root-only name.  Only a
            # successful RENAME_NOREPLACE can consume the fixed attempt name;
            # its held directory FD distinguishes a completed rename from an
            # EEXIST race without adopting somebody else's namespace.
            os.mkdir(temporary, 0o700, dir_fd=parent_fd)
            temporary_exists = True
            directory_fd = os.open(
                temporary,
                os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC
                | getattr(os, "O_NOFOLLOW", 0),
                dir_fd=parent_fd,
            )
            journal = cls(
                path, parent_fd, directory_fd,
                expected_uid=expected_uid, expected_gid=expected_gid,
                fault_injector=fault_injector,
                required_prefix=required_prefix,
            )
            journal._fault("attempt-private-created", path.name)
            rename_noreplace_at(parent_fd, temporary, path.name)
            temporary_exists = False
            journal._fault("attempt-fixed-created", path.name)
            consumed = True
            info = os.fstat(journal.directory_fd)
            named = os.stat(path.name, dir_fd=parent_fd, follow_symlinks=False)
            if (
                not stat.S_ISDIR(info.st_mode)
                or not same_inode(info, named)
                or stat.S_IMODE(info.st_mode) != 0o700
                or info.st_uid != expected_uid or info.st_gid != expected_gid
            ):
                fail("attempt directory policy differs")
            os.fsync(parent_fd)
            journal.write_json("ATTEMPT", attempt)
            if owner is not None:
                owner.adopt(journal)
            return journal
        except BaseException as exc:
            if parent_fd < 0:
                close_fd_delta(reservation_baseline, ())
                raise
            # If renameat2 completed but an injected asynchronous exception
            # arrived before the assignment, the fixed name must equal our
            # already-held private directory.  EEXIST/racing foreign names do
            # not match and are never adopted or terminalized.
            if journal is not None:
                try:
                    held = os.fstat(journal.directory_fd)
                    named = os.stat(
                        path.name, dir_fd=parent_fd, follow_symlinks=False,
                    )
                except (FileNotFoundError, OSError):
                    pass
                else:
                    if same_inode(held, named):
                        consumed = True
                        temporary_exists = False
            if consumed:
                recovery_error: Optional[BaseException] = None
                if journal is not None and "ATTEMPT" not in journal.names:
                    try:
                        journal.write_json("ATTEMPT", attempt)
                    except BaseException as recovery_exc:
                        recovery_error = recovery_exc
                message = (
                    "attempt namespace consumed during reservation: "
                    + exception_text(exc)
                )
                if recovery_error is not None:
                    message += " | ATTEMPT recovery: " + exception_text(
                        recovery_error)
                raise AttemptConsumedError(
                    message, journal,
                ) from exc
            # The private name was proven absent immediately before mkdir.
            # Inspect it even when an exception landed before the following
            # Python assignment.  Under the explicit no-hostile-root boundary,
            # any such name is the incomplete private reservation and must be
            # removed before reporting that the fixed attempt was unconsumed.
            try:
                try:
                    private_info = os.stat(
                        temporary, dir_fd=parent_fd, follow_symlinks=False,
                    )
                except FileNotFoundError:
                    private_info = None
                if private_info is not None:
                    if (
                        not stat.S_ISDIR(private_info.st_mode)
                        or private_info.st_uid != expected_uid
                        or private_info.st_gid != expected_gid
                        or stat.S_IMODE(private_info.st_mode) != 0o700
                    ):
                        fail("attempt private reservation changed type")
                    retained = [parent_fd]
                    if journal is not None and journal.directory_fd >= 0:
                        retained.append(journal.directory_fd)
                    close_untracked_inode_fds(private_info, retained)
                    os.rmdir(temporary, dir_fd=parent_fd)
                    os.fsync(parent_fd)
                    temporary_exists = False
            finally:
                close_fd_delta(reservation_baseline, ())
                if journal is not None:
                    journal.directory_fd = -1
                    journal.parent_fd = -1
            raise

    def _verify_directory(self, *, mode: int = 0o700) -> None:
        self._ensure_directory()
        info = os.fstat(self.directory_fd)
        named = os.stat(
            self.path.name, dir_fd=self.parent_fd, follow_symlinks=False,
        )
        if (
            (info.st_dev, info.st_ino) != self.binding
            or (named.st_dev, named.st_ino) != self.binding
            or not stat.S_ISDIR(info.st_mode)
            or stat.S_IMODE(info.st_mode) != mode
            or info.st_uid != self.expected_uid or info.st_gid != self.expected_gid
        ):
            fail("attempt directory binding changed")

    def verify_descriptor_authority(self) -> None:
        if (
            self.parent_fd < 0 or self.directory_fd < 0
            or self.parent_fd_authority is None
            or self.directory_fd_authority is None
            or fd_authority_receipt(self.parent_fd)
            != self.parent_fd_authority
            or fd_authority_receipt(self.directory_fd)
            != self.directory_fd_authority
        ):
            fail("attempt journal descriptor authority changed")

    def _reserve(self, name: str, mode: int = 0o600,
                 owner: Optional[Dict[str, int]] = None,
                 *, transient: bool = False) -> int:
        if owner is not None and transient:
            fail("attempt artifact has two descriptor owners")
        if (
            self.complete or not name or "/" in name or name in self.names
            or name in (".", "..")
        ):
            fail("attempt append order differs")
        self._verify_directory()
        try:
            os.stat(name, dir_fd=self.directory_fd, follow_symlinks=False)
        except FileNotFoundError:
            pass
        else:
            fail("attempt artifact name already exists: " + name)
        temporary = self._private_name("artifact")
        try:
            os.stat(temporary, dir_fd=self.directory_fd, follow_symlinks=False)
        except FileNotFoundError:
            pass
        else:
            fail("attempt private artifact name already exists")
        fd = -1
        temporary_exists = False
        renamed = False
        try:
            fd = os.open(
                temporary, os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
                | getattr(os, "O_NOFOLLOW", 0), mode, dir_fd=self.directory_fd,
            )
            temporary_exists = True
            held = os.fstat(fd)
            self._fault("artifact-private-created", name)
            rename_noreplace_at(self.directory_fd, temporary, name)
            temporary_exists = False
            renamed = True
            self._fault("artifact-fixed-created", name)
            named = os.stat(
                name, dir_fd=self.directory_fd, follow_symlinks=False,
            )
            if not same_inode(held, named):
                fail("attempt artifact rename binding differs: " + name)
            self.names.append(name)
            if transient:
                self.transient_fds.append(fd)
                info = os.fstat(fd)
                self.transient_fd_stats[fd] = (info.st_dev, info.st_ino)
            elif owner is not None:
                owner[name] = fd
            return fd
        except BaseException as original:
            cleanup_error: Optional[BaseException] = None
            try:
                # A completed fixed-name rename is ours only when it equals the
                # held private inode.  This separates EEXIST from
                # return-then-raise.
                if fd >= 0:
                    try:
                        held = os.fstat(fd)
                        named = os.stat(
                            name, dir_fd=self.directory_fd,
                            follow_symlinks=False,
                        )
                    except (FileNotFoundError, OSError):
                        pass
                    else:
                        if same_inode(held, named):
                            renamed = True
                            temporary_exists = False
                # Inspect the pre-proved-absent private name even if an
                # exception interrupted assignment after O_EXCL creation.
                try:
                    private_info = os.stat(
                        temporary, dir_fd=self.directory_fd,
                        follow_symlinks=False,
                    )
                except FileNotFoundError:
                    private_info = None
                if private_info is not None:
                    if (
                        not stat.S_ISREG(private_info.st_mode)
                        or private_info.st_uid != self.expected_uid
                        or private_info.st_gid != self.expected_gid
                        or stat.S_IMODE(private_info.st_mode) != mode
                        or private_info.st_nlink != 1
                    ):
                        fail("attempt private artifact changed type: " + name)
                    close_untracked_inode_fds(
                        private_info,
                        (self.parent_fd, self.directory_fd, fd),
                    )
                    os.unlink(temporary, dir_fd=self.directory_fd)
                    temporary_exists = False
                # `_reserve` has not transferred ownership until it returns.
                # Roll back any fixed name created by this invocation.
                if renamed and fd >= 0:
                    try:
                        held = os.fstat(fd)
                        named = os.stat(
                            name, dir_fd=self.directory_fd,
                            follow_symlinks=False,
                        )
                        if not same_inode(held, named):
                            fail("failed attempt artifact binding differs: " + name)
                        os.unlink(name, dir_fd=self.directory_fd)
                        if name in self.names:
                            self.names.remove(name)
                        os.fsync(self.directory_fd)
                    except FileNotFoundError:
                        if name in self.names:
                            self.names.remove(name)
            except BaseException as exc:
                cleanup_error = exc
            finally:
                if transient and fd in self.transient_fds:
                    try:
                        self._close_transient(fd)
                    except BaseException as exc:
                        if cleanup_error is None:
                            cleanup_error = exc
                elif fd >= 0:
                    try:
                        close_fd_once(fd, "failed attempt artifact")
                    except BaseException as exc:
                        if cleanup_error is None:
                            cleanup_error = exc
            if cleanup_error is not None:
                if owner is not None and owner.get(name) == fd:
                    owner.pop(name, None)
                raise LaunchError(
                    exception_text(original) + " | reserve cleanup: "
                    + exception_text(cleanup_error)
                ) from original
            if owner is not None and owner.get(name) == fd:
                owner.pop(name, None)
            raise

    def _discard_owned_artifact(self, name: str, fd: int) -> None:
        """Remove a failed publication while its exact inode is still held."""
        failures: List[str] = []
        try:
            held = os.fstat(fd)
            named = os.stat(
                name, dir_fd=self.directory_fd, follow_symlinks=False,
            )
            if not same_inode(held, named):
                fail("failed attempt artifact name changed: " + name)
            os.unlink(name, dir_fd=self.directory_fd)
            if name in self.names:
                self.names.remove(name)
            os.fsync(self.directory_fd)
            try:
                os.stat(name, dir_fd=self.directory_fd, follow_symlinks=False)
            except FileNotFoundError:
                pass
            else:
                fail("failed attempt artifact name survived: " + name)
            if os.fstat(fd).st_nlink != 0:
                fail("failed attempt artifact retained a link: " + name)
        except BaseException as exc:
            failures.append(exception_text(exc))
        if failures:
            fail("failed attempt artifact cleanup differs: " + " | ".join(failures))

    def write_bytes(self, name: str, payload: bytes,
                    *, final_mode: int = 0o400) -> Dict[str, Any]:
        if type(payload) is not bytes or len(payload) > MAX_JOURNAL_ARTIFACT_BYTES:
            fail("attempt artifact payload bound differs: " + name)
        fd = -1
        try:
            fd = self._reserve(name, transient=True)
            offset = 0
            while offset < len(payload):
                written = os.write(fd, payload[offset:])
                if written <= 0:
                    fail("attempt artifact write was short: " + name)
                offset += written
            os.fchmod(fd, final_mode)
            self._fault("artifact-final-mode", name)
            os.fsync(fd)
            self._fault("artifact-file-fsync", name)
            info = os.fstat(fd)
            named = os.stat(name, dir_fd=self.directory_fd, follow_symlinks=False)
            if (
                full_stat(info) != full_stat(named) or info.st_size != len(payload)
                or info.st_nlink != 1 or stat.S_IMODE(info.st_mode) != final_mode
                or info.st_uid != self.expected_uid
                or info.st_gid != self.expected_gid
                or hash_fd(fd, info.st_size) != sha256_bytes(payload)
            ):
                fail("attempt artifact publication differs: " + name)
            os.fsync(self.directory_fd)
            self._fault("artifact-directory-fsync", name)
            receipt = {
                "bytes": len(payload), "device": info.st_dev,
                "inode": info.st_ino, "mode": final_mode,
                "sha256": sha256_bytes(payload),
            }
            return receipt
        except BaseException as exc:
            if fd < 0:
                candidates = [item for item in self.transient_fds
                              if self._transient_matches_name(item, name)]
                fd = candidates[-1] if candidates else -1
            if fd >= 0:
                try:
                    self._discard_owned_artifact(name, fd)
                except BaseException as cleanup_exc:
                    raise LaunchError(
                        exception_text(exc) + " | publication cleanup: "
                        + exception_text(cleanup_exc)
                    ) from exc
            raise
        finally:
            if fd >= 0:
                self._close_transient(fd)

    def write_json(self, name: str, value: Mapping[str, Any]) -> Dict[str, Any]:
        return self.write_bytes(name, canonical_bytes(dict(value)) + b"\n")

    def write_bytes_resumable(self, name: str, payload: bytes,
                              *, final_mode: int = 0o400,
                              ) -> Dict[str, Any]:
        """Recover an exact prior write across the caller return boundary."""
        self._drain_transients()
        if name not in self.names:
            return self.write_bytes(name, payload, final_mode=final_mode)
        if self.complete or type(payload) is not bytes:
            fail("resumable attempt artifact contract differs: " + name)
        self._verify_directory()
        fd = -1
        try:
            fd = self._open_transient(
                name,
                os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK
                | getattr(os, "O_NOFOLLOW", 0),
                dir_fd=self.directory_fd,
            )
            info = os.fstat(fd)
            named = os.stat(
                name, dir_fd=self.directory_fd, follow_symlinks=False,
            )
            digest = sha256_bytes(payload)
            if (
                not stat.S_ISREG(info.st_mode)
                or full_stat(info) != full_stat(named)
                or info.st_uid != self.expected_uid
                or info.st_gid != self.expected_gid or info.st_nlink != 1
                or stat.S_IMODE(info.st_mode) != final_mode
                or info.st_size != len(payload)
                or hash_fd(fd, info.st_size) != digest
            ):
                fail("resumable attempt artifact differs: " + name)
            return {
                "bytes": info.st_size, "device": info.st_dev,
                "inode": info.st_ino, "mode": final_mode,
                "sha256": digest,
            }
        finally:
            self._drain_transients()

    def open_stream(self, name: str) -> int:
        """Pre-reserve one parent-rostered stream for a trusted fork child."""
        return self._reserve(name, owner=self.stream_fds)

    def discard_empty_stream(self, name: str) -> None:
        """Roll back a parent-owned stream before any child receives it."""
        if self.complete or name not in self.stream_fds:
            fail("discarded attempt stream ownership differs: " + name)
        fd = self.stream_fds[name]
        info = os.fstat(fd)
        if (
            not stat.S_ISREG(info.st_mode) or info.st_uid != self.expected_uid
            or info.st_gid != self.expected_gid or info.st_nlink != 1
            or stat.S_IMODE(info.st_mode) != 0o600 or info.st_size != 0
        ):
            # A previous unlink may have completed before an asynchronous
            # exception reached Python.  That exact held inode is a valid
            # resumable intermediate only when its link count is now zero.
            if not (
                stat.S_ISREG(info.st_mode)
                and info.st_uid == self.expected_uid
                and info.st_gid == self.expected_gid
                and info.st_nlink == 0
                and stat.S_IMODE(info.st_mode) == 0o600
                and info.st_size == 0
            ):
                fail("discarded attempt stream policy differs: " + name)
        linked = info.st_nlink == 1
        if linked:
            named = os.stat(
                name, dir_fd=self.directory_fd, follow_symlinks=False,
            )
            if not same_inode(info, named) or name not in self.names:
                fail("discarded attempt stream name differs: " + name)
            try:
                os.unlink(name, dir_fd=self.directory_fd)
            except BaseException:
                # Suppress only a real-unlink-then-raise boundary.  A fixed
                # name that remains linked preserves the original failure.
                try:
                    os.stat(
                        name, dir_fd=self.directory_fd,
                        follow_symlinks=False,
                    )
                except FileNotFoundError:
                    if os.fstat(fd).st_nlink != 0:
                        raise
                else:
                    raise
        try:
            os.stat(name, dir_fd=self.directory_fd, follow_symlinks=False)
        except FileNotFoundError:
            pass
        else:
            fail("discarded attempt stream name survived: " + name)
        if os.fstat(fd).st_nlink != 0:
            fail("discarded attempt stream retained a link: " + name)
        if name in self.names:
            self.names.remove(name)
        os.fsync(self.directory_fd)
        close_dict_fd(self.stream_fds, name, "discarded attempt stream " + name)

    def seal_stream(self, name: str, fd: int, maximum: int) -> Tuple[bytes, Dict[str, Any]]:
        if (
            self.complete or name not in self.names or fd < 0
            or maximum < 1 or maximum > MAX_JOURNAL_ARTIFACT_BYTES
        ):
            fail("attempt stream seal contract differs: " + name)
        self._verify_directory()
        before = os.fstat(fd)
        named_before = os.stat(
            name, dir_fd=self.directory_fd, follow_symlinks=False,
        )
        if (
            not stat.S_ISREG(before.st_mode) or not same_inode(before, named_before)
            or before.st_nlink != 1 or before.st_uid != self.expected_uid
            or before.st_gid != self.expected_gid or before.st_size < 1
            or before.st_size > maximum
            or stat.S_IMODE(before.st_mode) not in (0o600, 0o400)
        ):
            fail("attempt stream pre-seal policy differs: " + name)
        os.fchmod(fd, 0o400)
        os.fsync(fd)
        data = read_held_file(fd, maximum, "attempt stream " + name)
        after = os.fstat(fd)
        named_after = os.stat(
            name, dir_fd=self.directory_fd, follow_symlinks=False,
        )
        if (
            full_stat(after) != full_stat(named_after)
            or stat.S_IMODE(after.st_mode) != 0o400
            or sha256_bytes(data) != hash_fd(fd, after.st_size)
        ):
            fail("attempt stream named/readback seal differs: " + name)
        os.fsync(self.directory_fd)
        return data, {
            "bytes": after.st_size, "device": after.st_dev,
            "inode": after.st_ino, "mode": 0o400,
            "sha256": sha256_bytes(data),
        }

    def finish(self) -> None:
        self._drain_transients()
        if (
            self.complete or self.names[-1:] != ["terminal.json"]
            or tuple(self.names[:len(self.required_prefix)])
            != self.required_prefix
        ):
            fail("attempt COMPLETE is not terminal-last")
        self._verify_directory()
        if sorted(os.listdir(self.directory_fd)) != sorted(self.names):
            fail("attempt namespace differs before COMPLETE")
        # COMPLETE authenticates a closed journal, not merely a list of
        # names.  In particular, a pre-reserved guardian stream must have
        # been durably sealed; a zero-byte/mode-0600 stream after a failed
        # guardian admission can never be hidden by an invalid terminal.
        total = 0
        for name in self.names:
            artifact_fd = -1
            try:
                artifact_fd = self._open_transient(
                    name,
                    os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK
                    | getattr(os, "O_NOFOLLOW", 0),
                    dir_fd=self.directory_fd,
                )
                before = os.fstat(artifact_fd)
                named_artifact = os.stat(
                    name, dir_fd=self.directory_fd, follow_symlinks=False,
                )
                if (
                    not stat.S_ISREG(before.st_mode)
                    or full_stat(before) != full_stat(named_artifact)
                    or before.st_uid != self.expected_uid
                    or before.st_gid != self.expected_gid
                    or before.st_nlink != 1
                    or stat.S_IMODE(before.st_mode) != 0o400
                    or before.st_size < 0
                    or before.st_size > MAX_JOURNAL_ARTIFACT_BYTES
                ):
                    fail("attempt pre-COMPLETE artifact policy differs: " + name)
                total += before.st_size
                if total > MAX_RAW_BYTES + 32 * MAX_JOURNAL_ARTIFACT_BYTES:
                    fail("attempt pre-COMPLETE aggregate byte bound differs")
                digest = hash_fd(artifact_fd, before.st_size)
                after = os.fstat(artifact_fd)
                named_after = os.stat(
                    name, dir_fd=self.directory_fd, follow_symlinks=False,
                )
                if (
                    digest != hash_fd(artifact_fd, after.st_size)
                    or full_stat(after) != full_stat(before)
                    or full_stat(named_after) != full_stat(before)
                ):
                    fail("attempt pre-COMPLETE artifact changed: " + name)
            finally:
                self._drain_transients()
        complete = self_hashed(COMPLETE_SCHEMA, {
            "attempt_id": ATTEMPT_ID,
            "journal_entries": list(self.names),
            "terminal_name": "terminal.json",
        })
        temporary = self._private_name("complete")
        try:
            os.stat(temporary, dir_fd=self.directory_fd, follow_symlinks=False)
        except FileNotFoundError:
            pass
        else:
            fail("attempt COMPLETE private name already exists")
        fd = -1
        renamed = False
        sealed_info: Optional[os.stat_result] = None
        try:
            data = canonical_bytes(complete) + b"\n"
            fd = self._open_transient(
                temporary,
                os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
                | getattr(os, "O_NOFOLLOW", 0),
                0o600, dir_fd=self.directory_fd,
            )
            offset = 0
            while offset < len(data):
                written = os.write(fd, data[offset:])
                if written <= 0:
                    fail("attempt COMPLETE write was short")
                offset += written
            os.fchmod(fd, 0o000)
            self._fault("complete-final-mode", "COMPLETE")
            os.fsync(fd)
            self._fault("complete-file-fsync", "COMPLETE")
            private_info = os.stat(
                temporary, dir_fd=self.directory_fd, follow_symlinks=False,
            )
            held_info = os.fstat(fd)
            if (
                not same_inode(private_info, held_info)
                or held_info.st_size != len(data) or held_info.st_nlink != 1
                or held_info.st_uid != self.expected_uid
                or held_info.st_gid != self.expected_gid
                or stat.S_IMODE(held_info.st_mode) != 0
                or hash_fd(fd, held_info.st_size) != sha256_bytes(data)
            ):
                fail("attempt COMPLETE private seal differs")
            sealed_info = held_info
            # No journal artifact descriptor may survive publication of the
            # success marker.  The fixed rename below therefore cannot race a
            # later close/assignment exception into a false exact-FD claim.
            self._drain_transients()
            fd = -1
            if self.transient_fds:
                fail("attempt COMPLETE retained transient descriptors")
            production_root_seal = (
                self.expected_uid == 0 and self.expected_gid == 0
            )
            if production_root_seal:
                if os.geteuid() != 0:
                    fail("production COMPLETE seal requires root authority")
                # The production root can rename within its already-held
                # owner-read/execute directory via CAP_DAC_OVERRIDE.  Seal the
                # complete pre-marker namespace first; the fixed marker rename
                # is then the final namespace/state transition.
                os.fchmod(self.directory_fd, 0o500)
                self._fault("complete-directory-mode", "COMPLETE")
                os.fsync(self.directory_fd)
                self._fault("complete-directory-fsync", "COMPLETE")
            # COMPLETE becomes publicly named only after its exact immutable
            # bytes are durable and, in production, the directory is sealed.
            rename_noreplace_at(
                self.directory_fd, temporary, "COMPLETE",
            )
            renamed = True
            self.names.append("COMPLETE")
            self._fault("complete-fixed-created", "COMPLETE")
            os.fsync(self.directory_fd)
            self._fault("complete-rename-fsync", "COMPLETE")
            if not production_root_seal:
                # Scratch-only non-root journals cannot rename in a 0500
                # directory.  They still exercise content atomicity; the
                # production branch above is the authoritative order.
                os.fchmod(self.directory_fd, 0o500)
                self._fault("complete-directory-mode", "COMPLETE")
                os.fsync(self.directory_fd)
                self._fault("complete-directory-fsync", "COMPLETE")
            os.fsync(self.parent_fd)
            self._fault("complete-parent-fsync", "COMPLETE")
            complete_named = os.stat(
                "COMPLETE", dir_fd=self.directory_fd,
                follow_symlinks=False,
            )
            directory = os.fstat(self.directory_fd)
            named = os.stat(
                self.path.name, dir_fd=self.parent_fd, follow_symlinks=False,
            )
            if (
                sealed_info is None
                or not same_inode(sealed_info, complete_named)
                or complete_named.st_size != len(data)
                or complete_named.st_nlink != 1
                or stat.S_IMODE(complete_named.st_mode) != 0
                or complete_named.st_uid != self.expected_uid
                or complete_named.st_gid != self.expected_gid
                or complete_named.st_mtime_ns != sealed_info.st_mtime_ns
                or stat.S_IMODE(directory.st_mode) != 0o500
                or directory.st_uid != self.expected_uid
                or directory.st_gid != self.expected_gid
                or (directory.st_dev, directory.st_ino) != self.binding
                or (named.st_dev, named.st_ino) != self.binding
                or sorted(os.listdir(self.directory_fd)) != sorted(self.names)
                or self.transient_fds
            ):
                fail("attempt COMPLETE seal differs")
            self.complete = True
        except BaseException as exc:
            self.complete = False
            cleanup_failures: List[str] = []
            try:
                fixed_after_error = os.stat(
                    "COMPLETE", dir_fd=self.directory_fd,
                    follow_symlinks=False,
                )
            except (FileNotFoundError, OSError):
                fixed_after_error = None
            if fixed_after_error is not None:
                renamed = True
                if "COMPLETE" not in self.names:
                    self.names.append("COMPLETE")
            try:
                os.fchmod(self.directory_fd, 0o700)
            except BaseException as cleanup_exc:
                cleanup_failures.append("directory restore: " + exception_text(cleanup_exc))
            candidate = "COMPLETE" if renamed else temporary
            try:
                candidate_info = os.stat(
                    candidate, dir_fd=self.directory_fd,
                    follow_symlinks=False,
                )
            except FileNotFoundError:
                candidate_info = None
            except BaseException as cleanup_exc:
                cleanup_failures.append(
                    "COMPLETE candidate inspection: "
                    + exception_text(cleanup_exc))
                candidate_info = None
            if candidate_info is not None:
                try:
                    matching = [
                        item for item in self.transient_fds
                        if self._transient_matches_name(item, candidate)
                    ]
                    cleanup_fd = matching[-1] if matching else (
                        self._open_transient(
                            candidate,
                            getattr(os, "O_PATH", os.O_RDONLY)
                            | os.O_CLOEXEC
                            | getattr(os, "O_NOFOLLOW", 0),
                            dir_fd=self.directory_fd,
                        )
                    )
                    held = os.fstat(cleanup_fd)
                    if (
                        not same_inode(held, candidate_info)
                        or not stat.S_ISREG(held.st_mode)
                        or held.st_uid != self.expected_uid
                        or held.st_gid != self.expected_gid
                        or held.st_nlink != 1
                        or stat.S_IMODE(held.st_mode) not in (0o600, 0)
                    ):
                        fail("failed COMPLETE cleanup binding differs")
                    os.unlink(candidate, dir_fd=self.directory_fd)
                    if "COMPLETE" in self.names:
                        self.names.remove("COMPLETE")
                    if os.fstat(cleanup_fd).st_nlink != 0:
                        fail("failed COMPLETE retained a namespace link")
                    os.fsync(self.directory_fd)
                except BaseException as cleanup_exc:
                    cleanup_failures.append(
                        "COMPLETE removal: " + exception_text(cleanup_exc))
            elif "COMPLETE" in self.names:
                self.names.remove("COMPLETE")
            try:
                self._drain_transients()
            except BaseException as cleanup_exc:
                cleanup_failures.append(
                    "COMPLETE descriptor closure: "
                    + exception_text(cleanup_exc))
            try:
                os.fsync(self.parent_fd)
            except BaseException as cleanup_exc:
                cleanup_failures.append("parent sync: " + exception_text(cleanup_exc))
            if cleanup_failures:
                raise LaunchError(
                    exception_text(exc) + " | COMPLETE cleanup: "
                    + " | ".join(cleanup_failures)
                ) from exc
            raise

    def close(self) -> None:
        failures: List[str] = []
        try:
            self._drain_transients()
        except BaseException as exc:
            failures.append("transient descriptors: " + exception_text(exc))
        for name in list(self.stream_fds):
            try:
                close_dict_fd(
                    self.stream_fds, name, "attempt stream " + name,
                )
            except BaseException as exc:
                failures.append(name + ": " + exception_text(exc))
        for member in ("directory_fd", "parent_fd"):
            try:
                close_object_fd(self, member, "attempt journal " + member)
            except BaseException as exc:
                failures.append(member + ": " + exception_text(exc))
        if failures:
            fail("attempt journal descriptor closure differs: "
                 + " | ".join(failures))


@dataclass(frozen=True)
class BuildSealConfig:
    harness_source: Path
    current_source: Path
    parent_source: Path
    build_root: Path
    expected_harness_commit: str


@dataclass(frozen=True)
class LaunchConfig:
    expected_harness_commit: str


def implementation_configure_argv(source: Path, build: Path) -> List[str]:
    return [
        str(CMAKE), "-S", str(source), "-B", str(build), "-G", "Ninja",
        "-DCMAKE_BUILD_TYPE=Release",
        "-DCMAKE_C_COMPILER=" + str(CC),
        "-DCMAKE_CXX_COMPILER=" + str(CXX),
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
        "-DCMAKE_C_FLAGS=", "-DCMAKE_CXX_FLAGS=",
        "-DCMAKE_C_FLAGS_RELEASE=-O3 -DNDEBUG",
        "-DCMAKE_CXX_FLAGS_RELEASE=-O3 -DNDEBUG",
        "-DCMAKE_EXE_LINKER_FLAGS=", "-DCMAKE_SHARED_LINKER_FLAGS=",
        "-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=OFF",
        "-DBUILD_SHARED_LIBS=ON", "-DBUILD_TESTS=OFF", "-DBUILD_CODEC_V2=ON",
        "-DMARCH_NATIVE=OFF", "-DWH_LTO=OFF", "-DWH_PGO_MODE=OFF",
        "-DWIREHAIR_BUILD_BOTH=OFF", "-DWIREHAIR_BUILD_TOOLS=OFF",
        "-DWIREHAIR_BUILD_BENCHMARKS=OFF",
        "-DWIREHAIR_ENABLE_SCHEDULED_TESTS=OFF",
        "-DWIREHAIR_ENABLE_LIBFUZZER=OFF",
    ]


def role_configure_argv(
    harness: Path, implementation: Path, implementation_library: Path,
    build: Path, role: str, harness_commit: str, implementation_commit: str,
) -> List[str]:
    if role not in ("current", "parent"):
        fail("role configure role differs")
    return [
        str(CMAKE), "-S", str(harness / "bench/public_facade_timing"),
        "-B", str(build), "-G", "Ninja", "-DCMAKE_BUILD_TYPE=Release",
        "-DCMAKE_CXX_COMPILER=" + str(CXX),
        "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
        "-DCMAKE_CXX_FLAGS=",
        "-DCMAKE_CXX_FLAGS_RELEASE=-O3 -DNDEBUG",
        "-DCMAKE_EXE_LINKER_FLAGS=",
        "-DCMAKE_INTERPROCEDURAL_OPTIMIZATION=OFF",
        "-DWIREHAIR_ROLE=" + role,
        "-DWIREHAIR_ROLE_SOURCE_ROOT=" + str(implementation),
        "-DWIREHAIR_HARNESS_SOURCE_ROOT=" + str(harness),
        "-DWIREHAIR_ROLE_LIBRARY=" + str(implementation_library),
        "-DWIREHAIR_HARNESS_GIT_COMMIT=" + harness_commit,
        "-DWIREHAIR_IMPLEMENTATION_GIT_COMMIT=" + implementation_commit,
    ]


def build_systemd_run_argv(config: BuildSealConfig) -> List[str]:
    properties = (
        "User=0", "Group=0", "SupplementaryGroups=", "Restart=no",
        "ExitType=cgroup", "KillMode=control-group", "SendSIGKILL=yes",
        "TimeoutStopSec=1s", "RuntimeMaxSec=2400s", "AllowedCPUs=0-1",
        "AllowedMemoryNodes=0", "CPUAffinity=0-1",
        "TasksMax={}".format(BUILD_TASKS_MAX),
        "LimitCORE=0", "UMask=0077", "PrivateTmp=yes",
        "PrivateDevices=yes", "PrivateNetwork=yes", "DevicePolicy=closed",
        "ProtectControlGroups=yes", "ProtectProc=invisible",
    )
    argv = [
        "/usr/bin/systemd-run", "--quiet", "--wait", "--pipe", "--collect",
        "--service-type=exec", "--expand-environment=no",
        "--unit=" + BUILD_UNIT_NAME,
    ]
    argv.extend("--property=" + item for item in properties)
    argv.extend((
        str(ENV), "-i", "LANG=C.UTF-8", "LC_ALL=C.UTF-8",
        "PATH=/usr/bin:/bin", "TZ=UTC", str(PYTHON), "-I", "-B",
        str(INSTALLED_LAUNCHER), "--seal-build",
        "--harness-source-dir", str(config.harness_source),
        "--current-source-dir", str(config.current_source),
        "--parent-source-dir", str(config.parent_source),
        "--build-root", str(config.build_root),
        "--expected-harness-commit", config.expected_harness_commit,
    ))
    return argv


class AuthorityArgumentParser(argparse.ArgumentParser):
    def error(self, message: str) -> None:
        fail("argument validation: " + message)


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    tokens = list(argv)
    if tokens == ["--selftest"]:
        return argparse.Namespace(selftest=True, seal_build=False, execute=False)
    seal_names = (
        "--harness-source-dir", "--current-source-dir", "--parent-source-dir",
        "--build-root", "--expected-harness-commit",
    )
    if (
        len(tokens) == 11 and tokens[0] == "--seal-build"
        and tuple(tokens[index] for index in range(1, 11, 2)) == seal_names
        and not any("=" in token for token in tokens)
    ):
        parser = AuthorityArgumentParser(allow_abbrev=False, add_help=False)
        parser.add_argument("--seal-build", action="store_true", required=True)
        parser.add_argument("--harness-source-dir", type=Path, required=True)
        parser.add_argument("--current-source-dir", type=Path, required=True)
        parser.add_argument("--parent-source-dir", type=Path, required=True)
        parser.add_argument("--build-root", type=Path, required=True)
        parser.add_argument("--expected-harness-commit", required=True)
        args = parser.parse_args(tokens)
        args.selftest = False
        args.execute = False
        paths = (
            args.harness_source_dir, args.current_source_dir,
            args.parent_source_dir, args.build_root,
        )
        for index, path in enumerate(paths):
            exact_absolute(path, "seal-build path {}".format(index))
        if len(set(paths)) != len(paths) or any(
            left in right.parents or right in left.parents
            for left_index, left in enumerate(paths)
            for right in paths[left_index + 1:]
        ):
            fail("seal-build paths overlap")
        if LOWER40.fullmatch(args.expected_harness_commit) is None:
            fail("seal-build harness commit differs")
        return args
    if (
        len(tokens) == 3 and tokens[0] == "--execute"
        and tokens[1] == "--expected-harness-commit"
        and not any("=" in token for token in tokens)
    ):
        parser = AuthorityArgumentParser(allow_abbrev=False, add_help=False)
        parser.add_argument("--execute", action="store_true", required=True)
        parser.add_argument("--expected-harness-commit", required=True)
        args = parser.parse_args(tokens)
        args.selftest = False
        args.seal_build = False
        if LOWER40.fullmatch(args.expected_harness_commit) is None:
            fail("execute harness commit differs")
        return args
    fail("argument vector does not match one exact public mode")


def selftest() -> int:
    commit = "1" * 40
    args = parse_args(["--execute", "--expected-harness-commit", commit])
    if not args.execute or args.expected_harness_commit != commit:
        fail("launcher argument selftest differs")
    command = systemd_run_argv(commit)
    for required in (
        "--property=RuntimeMaxSec={}s".format(SERVICE_DEADLINE_SECONDS),
        "--property=Restart=no",
        "--property=ExitType=cgroup", "--property=KillMode=control-group",
        "--property=AllowedCPUs=120-123", "--property=CPUAffinity=123",
    ):
        if required not in command:
            fail("launcher systemd selftest lacks " + required)
    if not SERVICE_DEADLINE_SECONDS > (
        PRE_RELEASE_SECONDS + EXTERNAL_DEADLINE_SECONDS
        + POST_CONTROLLER_SECONDS + SERVICE_ACTIVATION_STOP_MARGIN_SECONDS
    ):
        fail("launcher service backstop margin is not strict")
    if not (
        CONTROLLER_ADMISSION_SECONDS + INTERNAL_DEADLINE_SECONDS
        + POST_CONTROLLER_SECONDS < EXTERNAL_DEADLINE_SECONDS
    ):
        fail("launcher scientific/post-controller deadline composition differs")
    value = self_hashed("selftest.v1", {"value": 7})
    if value["receipt_sha256"] is None:
        fail("launcher receipt selftest differs")
    print("WH2 V2 facade timing launcher selftest passed")
    return 0


def require_root_mode(expected_unit: str,
                      expected_affinity: Iterable[int], *,
                      expected_cgroup: Optional[str] = None) -> Dict[str, Any]:
    if (
        not hasattr(os, "getresuid") or os.getresuid() != (0, 0, 0)
        or not hasattr(os, "getresgid") or os.getresgid() != (0, 0, 0)
        or os.getgroups()
    ):
        fail("root mode credentials differ")
    if (
        not sys.flags.isolated or not sys.dont_write_bytecode
        or sys.flags.optimize != 0 or dict(os.environ) != SEALED_ENVIRONMENT
        or Path(sys.executable).resolve(strict=True) != PYTHON.resolve(strict=True)
        or Path(sys.argv[0]).resolve(strict=True) != INSTALLED_LAUNCHER
    ):
        fail("root launcher interpreter/environment boundary differs")
    require_root_owned_nonwritable_chain(INSTALLED_LAUNCHER)
    python_info = os.stat(str(PYTHON.resolve(strict=True)), follow_symlinks=False)
    live_python_info = os.stat("/proc/self/exe", follow_symlinks=True)
    if (
        full_stat(python_info) != full_stat(live_python_info)
        or python_info.st_nlink != FROZEN_TOOL_NLINK[PYTHON.resolve(strict=True)]
    ):
        fail("root launcher Python image binding differs")
    affinity = sorted(os.sched_getaffinity(0))
    if affinity != sorted(set(expected_affinity)):
        fail("root launcher affinity differs")
    raw = Path("/proc/self/cgroup").read_bytes()
    if len(raw) > 64 * 1024:
        fail("root launcher cgroup receipt exceeds bound")
    try:
        lines = raw.decode("ascii").splitlines()
        fields = [line.split(":", 2) for line in lines]
    except (UnicodeDecodeError, ValueError) as exc:
        fail("root launcher cgroup receipt malformed: " + exception_text(exc))
    unified = [field[2] for field in fields if len(field) == 3 and field[:2] == ["0", ""]]
    if len(unified) != 1 or (
        expected_cgroup is None and Path(unified[0]).name != expected_unit
    ) or (
        expected_cgroup is not None and unified[0] != expected_cgroup
    ):
        fail("root launcher is not in the exact delegated service")
    if (
        expected_cgroup is not None
        and Path(expected_cgroup).parent.name != expected_unit
    ):
        fail("root launcher supervisor cgroup escaped its exact service")
    return {
        "affinity": affinity,
        "cgroup": unified[0],
        "gid": 0,
        "groups": [],
        "uid": 0,
    }


def toolchain_receipt() -> Dict[str, Any]:
    receipts: List[Dict[str, Any]] = []
    for path, expected in sorted(
        FROZEN_TOOL_SHA256.items(), key=lambda item: str(item[0])
    ):
        actual_path = path
        digest, info = hash_path(
            actual_path, 256 * 1024 * 1024, "toolchain " + str(path),
            expected_nlink=FROZEN_TOOL_NLINK.get(actual_path, 1),
        )
        if (
            digest != expected or info.st_uid != 0 or info.st_gid != 0
            or stat.S_IMODE(info.st_mode) != 0o755
        ):
            fail("frozen toolchain image differs: " + str(path))
        receipts.append({
            "bytes": info.st_size,
            "path": str(actual_path),
            "sha256": digest,
            "stat": stat_receipt(info),
        })
    symlinks: List[Dict[str, Any]] = []
    for path, expected_target in sorted(
        FROZEN_TOOL_SYMLINKS.items(), key=lambda item: str(item[0])
    ):
        info = os.lstat(str(path))
        target = os.readlink(str(path))
        if (
            not stat.S_ISLNK(info.st_mode) or info.st_uid != 0
            or info.st_gid != 0 or stat.S_IMODE(info.st_mode) != 0o777
            or target != expected_target
        ):
            fail("frozen tool symlink differs: " + str(path))
        resolved = path.resolve(strict=True)
        if resolved not in FROZEN_TOOL_SHA256:
            fail("frozen tool symlink target lacks a canonical image: " + str(path))
        target_info = os.stat(str(resolved), follow_symlinks=False)
        followed_info = os.stat(str(path), follow_symlinks=True)
        if full_stat(target_info) != full_stat(followed_info):
            fail("frozen tool symlink target identity differs: " + str(path))
        symlinks.append({
            "path": str(path), "stat": stat_receipt(info), "target": target,
            "resolved_path": str(resolved),
        })
    discovery: Dict[str, Any] = {}
    for compiler, language, include_hash in (
        (CC, "c", "08d357250ffc68d3796b0a52f4ec512dea56178041b2b15b3199316280878467"),
        (CXX, "c++", "482534c810c5f241c5dd11145de9b95c1f0a3c19086eed75685e8b3f10e870bb"),
    ):
        search, search_error = _bounded_command(
            [str(compiler), "-print-search-dirs"], cwd=None,
            environment=BUILD_ENVIRONMENT_FIXED, timeout=10.0,
        )
        machine, machine_error = _bounded_command(
            [str(compiler), "-dumpmachine"], cwd=None,
            environment=BUILD_ENVIRONMENT_FIXED, timeout=10.0,
        )
        version, version_error = _bounded_command(
            [str(compiler), "-dumpfullversion"], cwd=None,
            environment=BUILD_ENVIRONMENT_FIXED, timeout=10.0,
        )
        _preprocessed, verbose = _bounded_command(
            [str(compiler), "-E", "-v", "-x", language, "/dev/null"],
            cwd=None, environment=BUILD_ENVIRONMENT_FIXED, timeout=10.0,
        )
        marker = b"#include <...> search starts here:\n"
        end_marker = b"End of search list.\n"
        try:
            start = verbose.index(marker)
            end = verbose.index(end_marker, start) + len(end_marker)
        except ValueError as exc:
            fail("compiler include search receipt differs: " + exception_text(exc))
        include_block = verbose[start:end]
        if (
            search_error or machine_error or version_error
            or sha256_bytes(search)
            != "35d3d157615681245967db80fd6d13a54e57f21d4a19cf0c7d28ffd3bc0dc309"
            or machine != b"x86_64-linux-gnu\n"
            or version != b"13.3.0\n"
            or sha256_bytes(include_block) != include_hash
        ):
            fail("compiler search/profile discovery differs: " + str(compiler))
        discovery[str(compiler)] = {
            "dumpfullversion": "13.3.0", "dumpmachine": "x86_64-linux-gnu",
            "include_search_sha256": include_hash,
            "search_dirs_sha256": sha256_bytes(search),
        }
    git_exec, git_exec_error = _bounded_command(
        [str(GIT), "--exec-path"], cwd=None, environment=GIT_ENVIRONMENT,
        timeout=10.0,
    )
    if (
        git_exec_error or git_exec != b"/usr/lib/git-core\n"
        or sha256_bytes(git_exec)
        != "626caad5dbe8e41b685b1b423515134e0694d5a073ed31026ad4414d14e88078"
    ):
        fail("Git backend search path differs")
    discovery["git_exec_path"] = "/usr/lib/git-core"
    value = {
        "discovery": discovery,
        "entries": receipts,
        "schema": "wirehair.wh2.v2-facade-timing-toolchain.v1",
        "sha256": None,
        "symlinks": symlinks,
    }
    value["sha256"] = sha256_bytes(canonical_bytes(value))
    return value


def _fresh_directory(path: Path, *, uid: int, gid: int,
                     mode: int = 0o700) -> None:
    exact_absolute(path, "fresh directory")
    parent = path.parent
    exact_absolute(parent, "fresh directory parent", must_exist=True)
    parent_info = os.stat(str(parent), follow_symlinks=False)
    if (
        not stat.S_ISDIR(parent_info.st_mode)
        or parent_info.st_uid != 0 or parent_info.st_gid != 0
        or stat.S_IMODE(parent_info.st_mode) & 0o022
    ):
        fail("fresh directory parent is not root protected")
    try:
        os.mkdir(str(path), mode)
    except FileExistsError:
        fail("fresh directory already exists: " + str(path))
    os.chown(str(path), uid, gid, follow_symlinks=False)
    # systemd intentionally launches the root sealer with UMask=0077.  The
    # requested campaign mode is an authority property, not an umask hint, so
    # establish it explicitly after the exclusive create and ownership handoff.
    os.chmod(str(path), mode, follow_symlinks=False)
    info = os.stat(str(path), follow_symlinks=False)
    if (
        not stat.S_ISDIR(info.st_mode) or info.st_uid != uid or info.st_gid != gid
        or stat.S_IMODE(info.st_mode) != mode
    ):
        fail("fresh directory policy differs: " + str(path))


def _run_git_outside(
    root: Path, arguments: Sequence[str], timeout: float,
    *, environment: Mapping[str, str] = GIT_ENVIRONMENT,
) -> bytes:
    output, error = _bounded_command(
        [str(GIT), *list(arguments)], cwd=root, environment=environment,
        timeout=timeout,
    )
    if error:
        fail("snapshot Git command wrote stderr: " + sha256_bytes(error))
    return output


def snapshot_clone_environment(source: Path) -> Dict[str, str]:
    """Authorize Git's local upload-pack for one exact UID-1000 input.

    ``git clone --no-local`` launches a separate local ``git-upload-pack``.
    Git deliberately strips command-line ``-c safe.directory=...`` state
    before that helper, but preserves ``SUDO_UID`` and documents it as the
    additional owner root may trust.  Bind both repository directories to the
    campaign identity before supplying that narrowly scoped environment.
    """
    exact_absolute(source, "snapshot clone source", must_exist=True)
    git_directory = source / ".git"
    try:
        if git_directory.resolve(strict=True) != git_directory:
            fail("snapshot clone Git directory traverses a symbolic link")
        for path in (source, git_directory):
            info = os.stat(str(path), follow_symlinks=False)
            if (
                not stat.S_ISDIR(info.st_mode)
                or info.st_uid != CAMPAIGN_UID
                or info.st_gid != CAMPAIGN_GID
            ):
                fail("snapshot clone source ownership differs: " + str(path))
    except OSError as exc:
        fail("snapshot clone source policy is unavailable: " + exception_text(exc))
    environment = dict(GIT_ENVIRONMENT)
    environment["SUDO_UID"] = str(CAMPAIGN_UID)
    return environment


def seal_directory_tree(path: Path, *, allow_symlinks: bool = False) -> None:
    """Make a completed tree root-owned and non-writable, bottom-up."""
    exact_absolute(path, "tree seal root", must_exist=True)
    directories: List[Path] = []
    for raw_root, raw_directories, raw_files in os.walk(
        str(path), topdown=True, followlinks=False,
    ):
        root = Path(raw_root)
        directories.append(root)
        for name in list(raw_directories) + list(raw_files):
            child = root / name
            info = os.lstat(str(child))
            if stat.S_ISLNK(info.st_mode):
                if not allow_symlinks:
                    fail("tree seal encountered a symbolic link: " + str(child))
                continue
            if stat.S_ISDIR(info.st_mode):
                continue
            if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
                fail("tree seal encountered a nonregular entry: " + str(child))
            executable = bool(stat.S_IMODE(info.st_mode) & 0o111)
            os.chown(str(child), 0, 0, follow_symlinks=False)
            os.chmod(str(child), 0o555 if executable else 0o444,
                     follow_symlinks=False)
    for directory in reversed(directories):
        os.chown(str(directory), 0, 0, follow_symlinks=False)
        os.chmod(str(directory), 0o555, follow_symlinks=False)


def materialize_snapshot(source: Path, destination: Path,
                         commit: str) -> GitTreeReceipt:
    """Copy one exact commit into a fresh, self-contained root-owned checkout."""
    clone_environment = snapshot_clone_environment(source)
    verify_git_tree(source, commit)
    if destination.exists():
        fail("snapshot destination already exists")
    _run_git_outside(
        destination.parent,
        (
            "clone", "--quiet", "--no-checkout", "--no-local", "--no-hardlinks",
            str(source), str(destination),
        ),
        120.0,
        environment=clone_environment,
    )
    _run_git_outside(
        destination,
        ("-c", "advice.detachedHead=false", "checkout", "--quiet", "--detach", commit),
        60.0,
    )
    copied = verify_git_tree(destination, commit)
    source_receipt = verify_git_tree(source, commit)
    def content_projection(receipt: GitTreeReceipt) -> List[Tuple[Any, ...]]:
        return [
            (
                entry["repository_relative_path"], entry["git_mode"],
                entry["git_blob_oid"], entry["bytes"], entry["sha256"],
            )
            for entry in receipt.entries
        ]

    if (
        copied.tree_oid != source_receipt.tree_oid
        or copied.tree_listing_sha256 != source_receipt.tree_listing_sha256
        or content_projection(copied) != content_projection(source_receipt)
    ):
        fail("fresh snapshot differs from its verified source")
    seal_directory_tree(destination)
    # The frozen controller executes as UID1000 and explicitly requires every
    # tracked source file to be UID1000:1000/0444 while every tracked directory
    # remains root:root/0555.  Ownership is not claimed as immutable against a
    # hostile concurrent UID1000 process; that exclusion is an explicit launch
    # trust assumption recorded by the root authority.
    for entry in copied.entries:
        tracked = destination / entry["repository_relative_path"]
        os.chown(str(tracked), CAMPAIGN_UID, CAMPAIGN_GID,
                 follow_symlinks=False)
        os.chmod(str(tracked), 0o444, follow_symlinks=False)
    sealed = verify_git_tree(
        destination, commit, required_uid=CAMPAIGN_UID,
        required_gid=CAMPAIGN_GID, sealed=True,
    )
    if content_projection(sealed) != content_projection(copied):
        fail("sealed snapshot content differs from its copied commit")
    return sealed


def build_environment(build_root: Path) -> Dict[str, str]:
    result = dict(BUILD_ENVIRONMENT_FIXED)
    result["HOME"] = str(build_root / ".build-home")
    result["TMPDIR"] = str(build_root / ".build-tmp")
    return result


def setpriv_build_argv(argv: Sequence[str]) -> List[str]:
    if not argv or not Path(argv[0]).is_absolute():
        fail("build child argv differs")
    return [
        str(SETPRIV), "--reuid", str(CAMPAIGN_UID),
        "--regid", str(CAMPAIGN_GID), "--clear-groups",
        "--bounding-set=-all", "--inh-caps=-all", "--ambient-caps=-all",
        "--no-new-privs", "--pdeathsig", "SIGKILL", *list(argv),
    ]


def run_build_command(name: str, argv: Sequence[str], cwd: Path,
                      environment: Mapping[str, str], timeout: float,
                      expected_stdout: Optional[bytes] = None,
                      expected_stderr: bytes = b"",
                      output_sink: Optional[List[bytes]] = None,
                      ) -> Dict[str, Any]:
    started_ns = time.monotonic_ns()
    output, error, security = bounded_build_child(
        argv, cwd, environment, timeout,
        stdout_limit=MAX_COMMAND_OUTPUT_BYTES,
        stderr_limit=8 * 1024 * 1024,
    )
    ended_ns = time.monotonic_ns()
    if expected_stdout is not None and output != expected_stdout:
        fail(name + " stdout differs from the exact expected transcript")
    if error != expected_stderr:
        fail(name + " stderr differs from the exact expected transcript")
    if output_sink is not None:
        if output_sink:
            fail(name + " output sink was not empty")
        output_sink.append(output)
    return self_hashed("wirehair.wh2.v2-facade-timing-build-command.v1", {
        "argv": list(argv),
        "ended_monotonic_ns": ended_ns,
        "name": name,
        "process_security": security,
        "stderr_bytes": len(error),
        "stderr_sha256": sha256_bytes(error),
        "stdout_bytes": len(output),
        "stdout_sha256": sha256_bytes(output),
        "started_monotonic_ns": started_ns,
    })


def process_start_ticks(pid: int) -> int:
    raw = Path("/proc/{}/stat".format(pid)).read_bytes()
    if len(raw) > 64 * 1024 or b"\0" in raw:
        fail("process stat framing differs")
    right = raw.rfind(b")")
    fields = raw[right + 2:].split() if right >= 2 else []
    if len(fields) < 20 or not fields[19].isdigit():
        fail("process start ticks differ")
    value = int(fields[19])
    if value <= 0:
        fail("process start ticks are not positive")
    return value


def process_topology(pid: int) -> Dict[str, int]:
    raw = read_bounded_path(
        Path("/proc/{}/stat".format(pid)), 64 * 1024, "process topology stat",
    )
    marker = raw.rfind(b") ")
    fields = raw[marker + 2:].split() if marker >= 0 else []
    if len(fields) < 4:
        fail("process topology stat is incomplete")
    try:
        values = [int(fields[index]) for index in (1, 2, 3)]
    except ValueError as exc:
        fail("process topology stat differs: " + exception_text(exc))
    if any(value <= 0 for value in values):
        fail("process topology values are not positive")
    return {
        "parent_pid": values[0], "process_group": values[1],
        "session": values[2],
    }


def current_service_pids() -> List[int]:
    raw = Path("/proc/self/cgroup").read_text(encoding="ascii")
    unified = [line[3:] for line in raw.splitlines() if line.startswith("0::/")]
    if len(unified) != 1:
        fail("service cgroup path differs")
    service = CGROUP_ROOT / unified[0].lstrip("/")
    if not service.is_dir():
        fail("service cgroup directory differs")
    # Include nested delegated leaves.  A child that moves itself below the
    # service root must not disappear from command-quiescence authority.
    return cgroup_descendant_pids(service)


def bounded_build_child(
    argv: Sequence[str], cwd: Path, environment: Mapping[str, str],
    timeout: float, *, stdout_limit: int, stderr_limit: int,
) -> Tuple[bytes, bytes, Dict[str, Any]]:
    if (
        not argv or any(type(value) is not str or not value for value in argv)
        or timeout <= 0.0 or stdout_limit < 0 or stderr_limit < 0
    ):
        fail("build child command contract differs")
    command_json = json.dumps(
        list(argv), sort_keys=False, separators=(",", ":"),
    )
    bootstrap = [
        str(PYTHON), "-I", "-S", "-c", BUILD_SECURITY_BOOTSTRAP,
        "READY_FD", "GATE_FD", command_json,
    ]
    ownership = BoundedCommandOwnership.create("build child " + argv[0])
    primary: Optional[BaseException] = None
    result: Optional[Tuple[bytes, bytes, Dict[str, Any]]] = None
    try:
        ready_r, ready_w = os.pipe2(os.O_CLOEXEC)
        ownership.fds.update({"ready_r": ready_r, "ready_w": ready_w})
        gate_r, gate_w = os.pipe2(os.O_CLOEXEC)
        ownership.fds.update({"gate_r": gate_r, "gate_w": gate_w})
        bootstrap[5] = str(ready_w)
        bootstrap[6] = str(gate_r)
        command = setpriv_build_argv(bootstrap)
        with tempfile.TemporaryFile() as stdout_file, tempfile.TemporaryFile() as stderr_file:
            process = ownership.spawn(
                command, cwd=str(cwd), env=dict(environment),
                stdin=subprocess.DEVNULL, stdout=stdout_file, stderr=stderr_file,
                close_fds=True, pass_fds=(ready_w, gate_r),
                start_new_session=True,
            )
            close_dict_fd(
                ownership.fds, "ready_w", "build child parent ready writer",
            )
            close_dict_fd(
                ownership.fds, "gate_r", "build child parent gate reader",
            )
            pidfd_baseline = open_fd_roster()
            try:
                pidfd = os.pidfd_open(process.pid, 0)
                ownership.fds["pidfd"] = pidfd
            except BaseException:
                close_fd_delta(pidfd_baseline, ownership.fds.values())
                raise
            start_ticks = process_start_ticks(process.pid)
            readable, _, _ = select.select([ready_r], [], [], 5.0)
            raw_security = os.read(ready_r, 65537) if readable else b""
            if not raw_security or len(raw_security) > 65536:
                fail("build child did not publish bounded security state")
            security = parse_canonical_document(
                raw_security, "build child security", 65536,
            )
            expected_caps = {
                "CapAmb": "0000000000000000",
                "CapBnd": "0000000000000000",
                "CapEff": "0000000000000000",
                "CapInh": "0000000000000000",
                "CapPrm": "0000000000000000",
                "NoNewPrivs": "1",
            }
            if (
                security != {
                    "affinity": [0, 1], "caps": expected_caps,
                    "environment": dict(environment),
                    "gid": [CAMPAIGN_GID] * 3,
                    "groups": [], "uid": [CAMPAIGN_UID] * 3,
                }
                or sorted(os.sched_getaffinity(process.pid)) != [0, 1]
                or process_start_ticks(process.pid) != start_ticks
                or process.poll() is not None
            ):
                fail("build child security boundary differs")
            ownership.may_have_released = True
            if os.write(gate_w, b"R") != 1:
                fail("build child release was short")
            close_dict_fd(
                ownership.fds, "gate_w", "build child parent gate writer",
            )
            try:
                code = process.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                fail("build child deadline expired")
            out_size = os.fstat(stdout_file.fileno()).st_size
            err_size = os.fstat(stderr_file.fileno()).st_size
            if out_size > stdout_limit or err_size > stderr_limit:
                fail("build child output exceeded its bound")
            stdout_file.seek(0)
            stderr_file.seek(0)
            output, error = stdout_file.read(), stderr_file.read()
            if code != 0:
                fail("build child failed rc={} stderr_sha256={}".format(
                    code, sha256_bytes(error)))
            security["command_argv_sha256"] = sha256_bytes(
                canonical_bytes(list(argv)))
            security["pid"] = process.pid
            security["process_start_ticks"] = start_ticks
            result = (output, error, security)
    except BaseException as exc:
        primary = exc
    cleanup: Optional[BaseException] = None
    try:
        descendants = ownership.quiesce()
        if descendants and primary is None:
            primary = LaunchError(
                "build command left one or more descendants: "
                + repr(descendants))
    except BaseException as exc:
        cleanup = exc
    if primary is not None:
        if cleanup is not None:
            raise LaunchError(
                exception_text(primary) + " | " + exception_text(cleanup)
            ) from primary
        raise primary
    if cleanup is not None:
        raise cleanup
    if result is None:
        fail("build child result handoff differs")
    return result


def read_json_file(path: Path, maximum: int, where: str) -> Any:
    digest, info = hash_path(path, maximum, where)
    del digest
    data = path.read_bytes()
    if len(data) != info.st_size:
        fail(where + " read size differs")
    try:
        return json.loads(data.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        fail(where + " is malformed: " + exception_text(exc))


def _resolve_command_path(directory: Path, value: str) -> Path:
    candidate = Path(value)
    if not candidate.is_absolute():
        candidate = directory / candidate
    return Path(os.path.normpath(str(candidate)))


@dataclass(frozen=True)
class CompileUnit:
    source: Path
    object_path: Path
    argv: Tuple[str, ...]
    effective_argv: Tuple[str, ...]
    dep_receipt: Dict[str, Any]


def exact_compile_argv_pair(
    build: Path, source: Path, relative_object: str,
    fixed_prefix: Sequence[str],
) -> Tuple[List[str], List[str]]:
    """Return the exact compile-database and Ninja command representations.

    CMake's ``compile_commands.json`` intentionally omits Ninja's depfile
    bookkeeping flags.  ``ninja -t commands`` includes those flags because
    Ninja expands the generated rule.  Both representations are build
    authority, but they are not byte-for-byte identical.
    """
    del build
    common = [str(CXX), *list(fixed_prefix)]
    compile_database = [
        *common, "-o", relative_object, "-c", str(source),
    ]
    effective_ninja = [
        *common,
        "-MD", "-MT", relative_object,
        "-MF", relative_object + ".d",
        "-o", relative_object, "-c", str(source),
    ]
    return compile_database, effective_ninja


def parse_compile_graph(
    build: Path, target: str, expected_sources: Sequence[Path],
    object_root: Path, fixed_prefix: Sequence[str],
) -> Tuple[List[CompileUnit], List[Dict[str, Any]]]:
    value = read_json_file(
        build / "compile_commands.json", 64 * 1024 * 1024,
        "compile commands",
    )
    if type(value) is not list or len(value) > 4096:
        fail("compile command roster shape differs")
    by_source: Dict[Path, CompileUnit] = {}
    expected = {path.resolve(strict=True) for path in expected_sources}
    effective_by_source: Dict[Path, List[str]] = {}
    for raw_command in command_lines_for_target(build, target):
        command = effective_ninja_command(raw_command)
        if not command or command[0] != str(CXX) or "-c" not in command:
            continue
        try:
            source = _resolve_command_path(
                build, command[command.index("-c") + 1],
            ).resolve(strict=True)
        except (ValueError, IndexError, OSError) as exc:
            fail("effective Ninja compile command differs: " + exception_text(exc))
        if source in expected:
            if source in effective_by_source:
                fail("effective Ninja graph repeats a compile source")
            effective_by_source[source] = command
    if set(effective_by_source) != expected:
        fail("effective Ninja graph lacks the exact compile source roster")
    dep_receipts: List[Dict[str, Any]] = []
    for index, entry in enumerate(value):
        if type(entry) is not dict or "file" not in entry or "directory" not in entry:
            fail("compile command {} shape differs".format(index))
        directory = Path(entry["directory"])
        if not directory.is_absolute() or directory.resolve(strict=True) != build:
            continue
        source = _resolve_command_path(directory, entry["file"]).resolve(strict=True)
        if source not in expected:
            continue
        if ("arguments" in entry) == ("command" in entry):
            fail("compile command must contain exactly one argv encoding")
        if "arguments" in entry:
            argv = list(entry["arguments"])
        else:
            try:
                argv = shlex.split(entry["command"], posix=True)
            except (TypeError, ValueError) as exc:
                fail("compile command parse failed: " + exception_text(exc))
        if (
            not argv or Path(argv[0]).resolve(strict=True) != CXX
            or any(token.startswith("@") for token in argv)
            or "-c" not in argv or str(source) not in argv
            or "-O3" not in argv or "-DNDEBUG" not in argv
        ):
            fail("effective compile argv differs for " + str(source))
        forbidden_prefixes = (
            "-march=", "-mtune=", "-flto", "-fprofile", "-fsanitize",
            "-fplugin", "-wrapper", "-include", "-imacros",
        )
        if any(token.startswith(forbidden_prefixes) for token in argv):
            fail("effective compile argv contains a forbidden feature")
        try:
            output_index = argv.index("-o") + 1
            object_path = _resolve_command_path(directory, argv[output_index])
        except (ValueError, IndexError) as exc:
            fail("compile object output differs: " + exception_text(exc))
        if object_path != object_root and object_root not in object_path.parents:
            continue
        relative_object = str(object_path.relative_to(build))
        compile_database_argv, effective_argv = exact_compile_argv_pair(
            build, source, relative_object, fixed_prefix,
        )
        if (
            argv != compile_database_argv
            or effective_by_source[source] != effective_argv
        ):
            fail("exact compile profile differs for " + str(source))
        if object_path in (unit.object_path for unit in by_source.values()):
            fail("compile graph repeats an object output")
        deps_out, deps_err = _bounded_command(
            [str(NINJA), "-C", str(build), "-t", "deps", relative_object],
            cwd=build, environment=BUILD_ENVIRONMENT_FIXED, timeout=10.0,
            stdout_limit=8 * 1024 * 1024,
        )
        if deps_err or not deps_out or b"deps mtime" not in deps_out.splitlines()[0]:
            fail("Ninja depfile closure is absent for " + str(object_path))
        dependency_paths: List[Path] = []
        for raw_line in deps_out.splitlines()[1:]:
            if not raw_line:
                continue
            if not raw_line.startswith(b"    "):
                fail("Ninja depfile row framing differs for " + str(object_path))
            try:
                raw_dependency = raw_line.strip().decode("utf-8")
            except UnicodeDecodeError as exc:
                fail("Ninja dependency path is not UTF-8: " + exception_text(exc))
            dependency = _resolve_command_path(build, raw_dependency).resolve(
                strict=True,
            )
            if dependency not in dependency_paths:
                dependency_paths.append(dependency)
        if not dependency_paths or source not in dependency_paths:
            fail("Ninja depfile closure lacks its exact translation unit")
        dependency_artifacts: List[Dict[str, Any]] = []
        for dependency in sorted(dependency_paths, key=str):
            digest, dependency_info = hash_path(
                dependency, MAX_TRACKED_FILE_BYTES,
                "compile dependency " + str(dependency),
                expected_nlink=None,
            )
            dependency_artifacts.append({
                "bytes": dependency_info.st_size,
                "path": str(dependency),
                "sha256": digest,
            })
        dep_receipt = {
            "bytes": len(deps_out),
            "dependencies": dependency_artifacts,
            "object_path": str(object_path),
            "sha256": sha256_bytes(deps_out),
        }
        unit = CompileUnit(
            source, object_path, tuple(argv), tuple(effective_argv), dep_receipt,
        )
        if source in by_source:
            fail("compile graph repeats a source")
        by_source[source] = unit
        dep_receipts.append(dep_receipt)
    if set(by_source) != expected:
        missing = sorted(str(path) for path in expected - set(by_source))
        fail("compile graph lacks exact sources: " + ",".join(missing))
    units = [by_source[path] for path in sorted(by_source, key=str)]
    return units, sorted(dep_receipts, key=lambda item: item["object_path"])


def implementation_compile_prefix(source: Path) -> List[str]:
    return [
        "-DWIREHAIR_BUILDING=1", "-DWIREHAIR_DLL=1", "-Dwirehair_EXPORTS",
        "-I" + str(source / "include"),
        "-O3", "-DNDEBUG", "-std=gnu++11", "-fPIC", "-Wall", "-Wextra",
    ]


def role_compile_prefix(implementation: Path, role: str,
                        harness_commit: str, implementation_commit: str) -> List[str]:
    if role not in ("current", "parent"):
        fail("role compile profile differs")
    return [
        "-DWIREHAIR_WH2_FACADE_ROLE_{}=1".format(role.upper()),
        '-DWIREHAIR_WH2_HARNESS_GIT_COMMIT="{}"'.format(harness_commit),
        '-DWIREHAIR_WH2_IMPLEMENTATION_GIT_COMMIT="{}"'.format(
            implementation_commit),
        "-I" + str(implementation / "include"),
        "-O3", "-DNDEBUG", "-std=gnu++11",
        "-O3", "-DNDEBUG", "-Wall", "-Wextra",
    ]


def command_lines_for_target(build: Path, target: str) -> List[List[str]]:
    output, error = _bounded_command(
        [str(NINJA), "-C", str(build), "-t", "commands", target],
        cwd=build, environment=BUILD_ENVIRONMENT_FIXED, timeout=10.0,
    )
    if error or not output:
        fail("Ninja target command closure differs: " + target)
    try:
        lines = output.decode("utf-8").splitlines()
        commands = [shlex.split(line, posix=True) for line in lines]
    except (UnicodeDecodeError, ValueError) as exc:
        fail("Ninja command closure malformed: " + exception_text(exc))
    if not commands or any(not command for command in commands):
        fail("Ninja target command closure is empty")
    return commands


def effective_ninja_command(command: Sequence[str]) -> List[str]:
    """Remove Ninja's exact no-op shell guard without hiding other shell code."""
    result = list(command)
    if result[:2] == [":", "&&"]:
        result = result[2:]
    if result[-2:] == ["&&", ":"]:
        result = result[:-2]
    if not result or any(token in ("&&", "||", ";", "|") for token in result):
        fail("Ninja command contains an unaudited shell composition")
    return result


def exact_link_command(build: Path, target: str, output: Path,
                       objects: Sequence[Path],
                       required_library: Optional[Path] = None,
                       expected_argv: Optional[Sequence[str]] = None) -> List[str]:
    candidates: List[List[str]] = []
    for raw_command in command_lines_for_target(build, target):
        command = effective_ninja_command(raw_command)
        for index, token in enumerate(command[:-1]):
            if token != "-o":
                continue
            candidate = _resolve_command_path(build, command[index + 1])
            if candidate == output:
                candidates.append(command)
    if len(candidates) != 1:
        fail("target has no unique exact link command: " + target)
    command = candidates[0]
    if Path(command[0]).resolve(strict=True) != CXX:
        fail("link command compiler driver differs")
    if any(token.startswith("@") for token in command):
        fail("link command uses an unaudited response file")
    linked_objects = {
        _resolve_command_path(build, token) for token in command
        if token.endswith(".o")
    }
    if linked_objects != set(objects):
        fail("link object closure differs for " + target)
    if required_library is not None:
        library_tokens = [
            token for token in command
            if not token.startswith("-")
            and _resolve_command_path(build, token) == required_library
        ]
        if len(library_tokens) != 1:
            fail("link command lacks the unique exact adjacent role library")
    forbidden = ("-flto", "-fprofile", "-fsanitize", "--wrap", "-wrapper")
    if any(any(part in token for part in forbidden) for token in command):
        fail("link command contains a forbidden feature")
    if expected_argv is not None and command != list(expected_argv):
        fail("exact link profile differs for " + target)
    return command


def implementation_link_argv(source: Path, build: Path,
                             objects: Sequence[Path]) -> List[str]:
    return [
        str(CXX), "-fPIC", "-O3", "-DNDEBUG",
        "-Wl,--version-script=" + str(source / "abi/wirehair.map"),
        "-shared", "-Wl,-soname,libwirehair.so.2",
        "-o", "libwirehair.so.2.0.0",
        *[str(path.relative_to(build)) for path in objects], "-lm",
    ]


def role_link_argv(build: Path, output: Path,
                   objects: Sequence[Path]) -> List[str]:
    # ``ninja -t commands`` preserves Ninja's shell escape for the literal
    # dollar sign.  The independently parsed ELF RUNPATH below proves that the
    # executed link semantics are the intended literal ``$ORIGIN``.
    return [
        str(CXX), "-O3", "-DNDEBUG", "-Wl,-z,now", "-Wl,-z,relro",
        *[str(path.relative_to(build)) for path in objects],
        "-o", str(output.relative_to(build)), "-Wl,-rpath,\\$ORIGIN",
        "libwirehair.so.2",
    ]


def validate_target_command_roster(build: Path, target: str,
                                   compile_count: int,
                                   ancillary: Sequence[Sequence[str]]) -> None:
    commands = [
        effective_ninja_command(command)
        for command in command_lines_for_target(build, target)
    ]
    compile_commands = [
        command for command in commands
        if command and command[0] == str(CXX) and "-c" in command
    ]
    link_commands = [
        command for command in commands
        if command and command[0] == str(CXX) and "-c" not in command
        and "-o" in command
    ]
    other = [
        command for command in commands
        if command not in compile_commands and command not in link_commands
    ]
    if (
        len(compile_commands) != compile_count or len(link_commands) != 1
        or other != [list(command) for command in ancillary]
        or len(commands) != compile_count + 1 + len(ancillary)
    ):
        fail("exact Ninja target command roster differs: " + target)


def canonical_link_argv(command: Sequence[str], build: Path,
                        required_library: Path) -> List[str]:
    """Path-normalize the one role library token for the frozen screen schema.

    The root build graph retains Ninja's raw effective argv.  The screen's
    frozen receipt adapter additionally requires the absolute role-library
    path as an argv member, so this representation normalizes only that one
    already-proven token and leaves every other argument byte-for-byte equal.
    """
    result: List[str] = []
    replacements = 0
    for token in command:
        if (
            not token.startswith("-")
            and _resolve_command_path(build, token) == required_library
        ):
            result.append(str(required_library))
            replacements += 1
        else:
            result.append(token)
    if replacements != 1 or str(required_library) not in result:
        fail("canonical link argv role-library binding differs")
    return result


def artifact_receipt(path: Path, maximum: int = 1024 * 1024 * 1024,
                     *, expected_nlink: Optional[int] = 1) -> Dict[str, Any]:
    digest, info = hash_path(
        path, maximum, "artifact " + str(path),
        expected_nlink=expected_nlink,
    )
    return {"bytes": info.st_size, "path": str(path), "sha256": digest}


def build_receipt_expected_nlink(
    roster_name: str, path: str, role_library_path: str,
) -> Optional[int]:
    """Freeze single-link campaign artifacts without rejecting system DSOs.

    The adjacent role library is root-created and must remain nlink=1.  Direct
    root-owned system dependencies can legitimately be multiply linked by the
    distribution; their full observed stat receipt is still cross-bound and
    rechecked before and after execution.
    """
    if roster_name not in (
        "source_manifest", "object_closure", "archive_closure",
        "dynamic_dependencies",
    ):
        fail("build-receipt roster name differs for link-count policy")
    return (
        None
        if roster_name == "dynamic_dependencies" and path != role_library_path
        else 1
    )


def object_receipts(units: Sequence[CompileUnit]) -> List[Dict[str, Any]]:
    result = [artifact_receipt(unit.object_path) for unit in units]
    paths = [item["path"] for item in result]
    if paths != sorted(set(paths)):
        result.sort(key=lambda item: item["path"])
    return result


def parse_direct_ldd(ldd_output: bytes, needed: Sequence[str],
                     adjacent_library: Path) -> List[Path]:
    """Resolve exactly the direct DT_NEEDED roster, excluding ldd transitives."""
    resolved_by_soname: Dict[str, Path] = {}
    needed_set = set(needed)
    if len(needed_set) != len(needed):
        fail("DT_NEEDED roster repeats a soname")
    for raw_line in ldd_output.splitlines():
        line = raw_line.strip()
        if not line or line.startswith(b"linux-vdso"):
            continue
        if b"=>" not in line:
            # The ELF interpreter and any direct absolute-name dependency are
            # outside this frozen worker's DT_NEEDED roster.
            continue
        left, right = line.split(b"=>", 1)
        try:
            soname = left.strip().decode("ascii")
        except UnicodeDecodeError as exc:
            fail("ldd soname is not ASCII: " + exception_text(exc))
        if soname not in needed_set:
            continue
        if soname in resolved_by_soname:
            fail("ldd repeats a direct dependency: " + soname)
        path_token = right.strip().split(b" ", 1)[0]
        try:
            path = Path(path_token.decode("utf-8"))
        except UnicodeDecodeError as exc:
            fail("ldd path is not UTF-8: " + exception_text(exc))
        if not path.is_absolute():
            fail("ldd direct dependency path is not absolute: " + soname)
        if soname == "libwirehair.so.2":
            if path != adjacent_library:
                fail("loader did not resolve the exact adjacent role library")
            resolved = path
        else:
            resolved = path.resolve(strict=True)
        resolved_by_soname[soname] = resolved
    if set(resolved_by_soname) != needed_set:
        missing = sorted(needed_set - set(resolved_by_soname))
        fail("ldd lacks direct dependencies: " + ",".join(missing))
    result = [resolved_by_soname[soname] for soname in needed]
    if len(result) != len(set(result)):
        fail("dynamic dependency path closure repeats a file")
    return result


def parse_runtime_ldd(ldd_output: bytes, adjacent_library: Path) -> List[Path]:
    """Return the full regular-file loader closure used for maps auditing."""
    result: List[Path] = []
    for raw_line in ldd_output.splitlines():
        line = raw_line.strip()
        if not line or line.startswith(b"linux-vdso"):
            continue
        if b"=>" in line:
            token = line.split(b"=>", 1)[1].strip().split(b" ", 1)[0]
        elif line.startswith(b"/"):
            token = line.split(b" ", 1)[0]
        else:
            fail("ldd runtime closure row differs")
        try:
            path = Path(token.decode("utf-8"))
        except UnicodeDecodeError as exc:
            fail("ldd runtime path is not UTF-8: " + exception_text(exc))
        if not path.is_absolute():
            fail("ldd runtime closure contains a nonabsolute path")
        resolved = path if path == adjacent_library else path.resolve(strict=True)
        if resolved not in result:
            result.append(resolved)
    if adjacent_library not in result:
        fail("ldd runtime closure lacks the adjacent role library")
    return result


def parse_elf(worker: Path, adjacent_library: Path) -> Dict[str, Any]:
    dynamic, dynamic_error = _bounded_command(
        [str(READELF), "-W", "-d", str(worker)], cwd=worker.parent,
        environment=BUILD_ENVIRONMENT_FIXED, timeout=10.0,
    )
    program, program_error = _bounded_command(
        [str(READELF), "-W", "-l", str(worker)], cwd=worker.parent,
        environment=BUILD_ENVIRONMENT_FIXED, timeout=10.0,
    )
    if dynamic_error or program_error:
        fail("ELF inspection wrote stderr")
    needed: List[str] = []
    runpaths: List[str] = []
    for raw_line in dynamic.splitlines():
        match = re.search(rb"\(NEEDED\).*\[([^]]+)\]", raw_line)
        if match:
            needed.append(match.group(1).decode("ascii"))
        match = re.search(rb"\(RUNPATH\).*\[([^]]+)\]", raw_line)
        if match:
            runpaths.append(match.group(1).decode("ascii"))
    expected_needed = [
        "libwirehair.so.2", "libstdc++.so.6", "libgcc_s.so.1", "libc.so.6",
    ]
    if (
        needed != expected_needed or runpaths != ["$ORIGIN"]
        or b"(RPATH)" in dynamic or b"BIND_NOW" not in dynamic
        or b"Flags: NOW PIE" not in dynamic
        or b"GNU_RELRO" not in program or b"GNU_STACK" not in program
        or b"Requesting program interpreter" not in program
    ):
        fail("worker ELF policy differs")
    ldd_out, ldd_err = _bounded_command(
        [str(LDD), str(worker)], cwd=worker.parent,
        environment=BUILD_ENVIRONMENT_FIXED, timeout=10.0,
    )
    if ldd_err or b"not found" in ldd_out:
        fail("worker dynamic closure could not be resolved")
    paths = parse_direct_ldd(ldd_out, needed, adjacent_library)
    runtime_paths = parse_runtime_ldd(ldd_out, adjacent_library)
    dependencies = sorted((
        artifact_receipt(
            path, expected_nlink=build_receipt_expected_nlink(
                "dynamic_dependencies", str(path), str(adjacent_library),
            ),
        )
        for path in paths
    ), key=lambda item: item["path"])
    return {
        "dt_needed": needed,
        "dynamic_dependencies": dependencies,
        "program_headers_sha256": sha256_bytes(program),
        "readelf_dynamic_sha256": sha256_bytes(dynamic),
        "runpath": "$ORIGIN",
        "runtime_map_dependencies": sorted(
            (artifact_receipt(path, expected_nlink=None) for path in runtime_paths),
            key=lambda item: item["path"],
        ),
    }


def source_manifest_entries(receipt: GitTreeReceipt,
                            authority: str) -> List[Dict[str, Any]]:
    result: List[Dict[str, Any]] = []
    for entry in receipt.entries:
        result.append({
            "authority": authority,
            "bytes": entry["bytes"],
            "git_blob_oid": entry["git_blob_oid"],
            "git_mode": entry["git_mode"],
            "path": str(receipt.root / entry["repository_relative_path"]),
            "repository_relative_path": entry["repository_relative_path"],
            "sha256": entry["sha256"],
        })
    return result


def write_root_file_exclusive(path: Path, data: bytes,
                              *, mode: int = 0o400) -> Dict[str, Any]:
    exact_absolute(path, "root publication path")
    parent_info = os.stat(str(path.parent), follow_symlinks=False)
    if (
        not stat.S_ISDIR(parent_info.st_mode) or parent_info.st_uid != 0
        or parent_info.st_gid != 0 or stat.S_IMODE(parent_info.st_mode) & 0o022
    ):
        fail("root publication parent differs")
    fd = os.open(
        str(path), os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
        | getattr(os, "O_NOFOLLOW", 0), mode,
    )
    try:
        offset = 0
        while offset < len(data):
            written = os.write(fd, data[offset:])
            if written <= 0:
                fail("root publication write was short")
            offset += written
        os.fchmod(fd, mode)
        os.fsync(fd)
        info = os.fstat(fd)
        if (
            info.st_uid != 0 or info.st_gid != 0 or info.st_nlink != 1
            or stat.S_IMODE(info.st_mode) != mode or info.st_size != len(data)
        ):
            fail("root publication policy differs")
    finally:
        os.close(fd)
    parent_fd = os.open(str(path.parent), os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
    try:
        os.fsync(parent_fd)
    finally:
        os.close(parent_fd)
    return artifact_receipt(path, MAX_DOCUMENT_BYTES)


def publish_root_file_noreplace(path: Path, data: bytes,
                                *, mode: int = 0o400) -> Dict[str, Any]:
    """Publish complete bytes atomically without exposing a partial success."""
    exact_absolute(path, "fixed root publication")
    if path.parent != BUILD_AUTHORITY_DIR or path != BUILD_AUTHORITY_PATH:
        fail("fixed authority publication path differs")
    current = Path("/")
    for part in path.parent.parts[1:]:
        current = current / part
        info = os.stat(str(current), follow_symlinks=False)
        if (
            not stat.S_ISDIR(info.st_mode) or info.st_uid != 0
            or info.st_gid != 0 or stat.S_IMODE(info.st_mode) & 0o022
        ):
            fail("fixed authority directory chain differs: " + str(current))
    authority_dir_info = os.stat(str(BUILD_AUTHORITY_DIR), follow_symlinks=False)
    if stat.S_IMODE(authority_dir_info.st_mode) != 0o700:
        fail("fixed authority directory mode differs")
    directory_fd = os.open(
        str(path.parent), os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC
        | getattr(os, "O_NOFOLLOW", 0),
    )
    temporary = ".{}.tmp.{}.{}".format(
        path.name, os.getpid(), time.monotonic_ns(),
    )
    fd = -1
    renamed = False
    try:
        held_parent = os.fstat(directory_fd)
        named_parent = os.stat(str(path.parent), follow_symlinks=False)
        if full_stat(held_parent) != full_stat(named_parent):
            fail("fixed authority parent binding differs")
        try:
            os.stat(path.name, dir_fd=directory_fd, follow_symlinks=False)
        except FileNotFoundError:
            pass
        else:
            fail("fixed authority target already exists")
        fd = os.open(
            temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC
            | getattr(os, "O_NOFOLLOW", 0), mode, dir_fd=directory_fd,
        )
        view = memoryview(data)
        while view:
            written = os.write(fd, view)
            if written <= 0:
                fail("fixed authority temporary write was short")
            view = view[written:]
        os.fchmod(fd, mode)
        os.fsync(fd)
        before = os.fstat(fd)
        if (
            not stat.S_ISREG(before.st_mode) or before.st_nlink != 1
            or before.st_uid != 0 or before.st_gid != 0
            or stat.S_IMODE(before.st_mode) != mode
            or before.st_size != len(data)
            or hash_fd(fd, before.st_size) != sha256_bytes(data)
        ):
            fail("fixed authority temporary seal differs")
        rename_noreplace_at(directory_fd, temporary, path.name)
        renamed = True
        os.fsync(directory_fd)
        named = os.stat(path.name, dir_fd=directory_fd, follow_symlinks=False)
        after = os.fstat(fd)
        if (
            full_stat(after) != full_stat(named) or after.st_nlink != 1
            or hash_fd(fd, after.st_size) != sha256_bytes(data)
        ):
            fail("fixed authority named/readback seal differs")
    finally:
        if fd >= 0:
            os.close(fd)
        if not renamed:
            try:
                os.unlink(temporary, dir_fd=directory_fd)
            except FileNotFoundError:
                pass
        os.close(directory_fd)
    return artifact_receipt(path, MAX_DOCUMENT_BYTES)


def _prepare_build_root(config: BuildSealConfig) -> Dict[str, Path]:
    exact_absolute(config.build_root, "build root")
    if config.build_root.exists():
        fail("build root is not fresh")
    _fresh_directory(config.build_root, uid=0, gid=0, mode=0o755)
    paths = {
        "sources": config.build_root / "sources",
        "builds": config.build_root / "builds",
        "harness": config.build_root / "sources/harness",
        "current_source": config.build_root / "sources/current",
        "parent_source": config.build_root / "sources/parent",
        "current_impl_build": config.build_root / "builds/current-implementation",
        "parent_impl_build": config.build_root / "builds/parent-implementation",
        "current_role_build": config.build_root / "builds/current-role",
        "parent_role_build": config.build_root / "builds/parent-role",
    }
    _fresh_directory(paths["sources"], uid=0, gid=0, mode=0o755)
    _fresh_directory(paths["builds"], uid=0, gid=0, mode=0o755)
    for name in (
        "current_impl_build", "parent_impl_build",
        "current_role_build", "parent_role_build",
    ):
        _fresh_directory(
            paths[name], uid=CAMPAIGN_UID, gid=CAMPAIGN_GID, mode=0o700,
        )
    for name in (".build-home", ".build-tmp"):
        _fresh_directory(
            config.build_root / name, uid=CAMPAIGN_UID,
            gid=CAMPAIGN_GID, mode=0o700,
        )
    return paths


def require_science_namespaces_absent() -> None:
    candidates = (
        ATTEMPT_DIR, CONTROLLER_PARENT, CONTROLLER_OUTPUT, SAMPLER_DIR, LOCK_PATH,
    )
    for path in candidates:
        for observed in (path, Path("/proc/1/root") / str(path).lstrip("/")):
            try:
                os.stat(str(observed), follow_symlinks=False)
            except FileNotFoundError:
                continue
            fail("science one-shot namespace already exists: " + str(observed))
    science_cgroup = Path("/sys/fs/cgroup/system.slice") / UNIT_NAME
    try:
        os.stat(str(science_cgroup), follow_symlinks=False)
    except FileNotFoundError:
        pass
    else:
        fail("science transient unit already exists during build sealing")


def _validate_required_sources(receipt: GitTreeReceipt,
                               required: Iterable[str], where: str) -> None:
    available = {entry["repository_relative_path"] for entry in receipt.entries}
    missing = sorted(set(required) - available)
    if missing:
        fail(where + " lacks frozen sources: " + ",".join(missing))


def validate_frozen_harness_sources(receipt: GitTreeReceipt) -> None:
    by_name = {
        entry["repository_relative_path"]: entry for entry in receipt.entries
    }
    expected = {
        CONTROLLER_RELATIVE: FROZEN_CONTROLLER_SHA256,
        WORKER_RELATIVE: FROZEN_WORKER_SOURCE_SHA256,
        ROLE_CMAKE_RELATIVE: FROZEN_ROLE_CMAKE_SHA256,
        HEALTH_ADAPTER_RELATIVE: FROZEN_HEALTH_ADAPTER_SHA256,
    }
    for relative, digest in expected.items():
        entry = by_name.get(relative)
        if entry is None or entry.get("sha256") != digest:
            fail("frozen harness source differs: " + relative)


def reverify_dependency_closure(depfiles: Mapping[str, Any]) -> str:
    for group_name in sorted(depfiles):
        group = depfiles[group_name]
        if type(group) is not list:
            fail("dependency closure group differs: " + group_name)
        for unit in group:
            dependencies = unit.get("dependencies") if type(unit) is dict else None
            if type(dependencies) is not list or not dependencies:
                fail("dependency closure unit differs: " + group_name)
            for entry in dependencies:
                path = Path(entry["path"])
                digest, info = hash_path(
                    path, MAX_TRACKED_FILE_BYTES,
                    "final compile dependency " + str(path),
                    expected_nlink=None,
                )
                if digest != entry["sha256"] or info.st_size != entry["bytes"]:
                    fail("compile dependency changed: " + str(path))
    return sha256_bytes(canonical_bytes(depfiles))


def content_receipt(path: Path, where: str) -> Dict[str, Any]:
    canonical = path.resolve(strict=True)
    digest, info = hash_path(
        canonical, MAX_TRACKED_FILE_BYTES, where, expected_nlink=None,
    )
    return {
        "bytes": info.st_size, "path": str(canonical), "sha256": digest,
        "stat": stat_receipt(info),
    }


def ninja_input_closure(build: Path, target: str) -> Dict[str, Any]:
    output, error = _bounded_command(
        [str(NINJA), "-C", str(build), "-t", "inputs", target],
        cwd=build, environment=BUILD_ENVIRONMENT_FIXED, timeout=20.0,
        stdout_limit=64 * 1024 * 1024,
    )
    if error or not output or b"\0" in output or b"\r" in output:
        fail("Ninja input closure differs: " + target)
    try:
        names = output.decode("utf-8").splitlines()
    except UnicodeDecodeError as exc:
        fail("Ninja input closure is not UTF-8: " + exception_text(exc))
    paths = sorted({
        _resolve_command_path(build, name).resolve(strict=True)
        for name in names if name
    }, key=str)
    if not paths:
        fail("Ninja input closure is empty: " + target)
    entries = [
        content_receipt(path, "Ninja input " + str(path)) for path in paths
    ]
    return {
        "entries": entries, "raw_sha256": sha256_bytes(output),
        "target": target,
    }


def linker_input_closure() -> List[Dict[str, Any]]:
    names = (
        "Scrt1.o", "crti.o", "crtbeginS.o", "crtendS.o", "crtn.o",
        "libstdc++.so", "libgcc_s.so", "libgcc.a", "libc.so", "libm.so",
    )
    result: List[Dict[str, Any]] = []
    for name in names:
        output, error = _bounded_command(
            [str(CXX), "-print-file-name=" + name], cwd=None,
            environment=BUILD_ENVIRONMENT_FIXED, timeout=10.0,
            stdout_limit=8192,
        )
        try:
            value = output.rstrip(b"\n").decode("utf-8")
        except UnicodeDecodeError as exc:
            fail("linker input discovery is not UTF-8: " + exception_text(exc))
        if error or output != (value + "\n").encode("utf-8") or value == name:
            fail("linker input discovery differs: " + name)
        receipt = content_receipt(Path(value), "linker input " + name)
        receipt["logical_name"] = name
        result.append(receipt)
    return result


def installed_launcher_receipt(harness: GitTreeReceipt) -> Dict[str, Any]:
    matches = [
        entry for entry in harness.entries
        if entry["repository_relative_path"] == SOURCE_LAUNCHER_RELATIVE
    ]
    if len(matches) != 1:
        fail("harness snapshot lacks one launcher source")
    digest, info = hash_path(
        INSTALLED_LAUNCHER, 8 * 1024 * 1024, "installed root launcher",
    )
    if (
        digest != matches[0]["sha256"] or info.st_size != matches[0]["bytes"]
        or info.st_uid != 0 or info.st_gid != 0
        or stat.S_IMODE(info.st_mode) != 0o555
    ):
        fail("installed launcher differs from the harness snapshot")
    return {
        "bytes": info.st_size, "path": str(INSTALLED_LAUNCHER),
        "sha256": digest, "source_git_blob_oid": matches[0]["git_blob_oid"],
        "source_relative_path": SOURCE_LAUNCHER_RELATIVE,
        "stat": stat_receipt(info),
    }


def _role_receipt(
    role: str, harness_commit: str, implementation_commit: str,
    harness: GitTreeReceipt, implementation: GitTreeReceipt,
    implementation_library: Path, role_build: Path,
    configure_argv: Sequence[str], build_argv: Sequence[str],
    worker_units: Sequence[CompileUnit], worker_link: Sequence[str],
    elf: Mapping[str, Any],
) -> Dict[str, Any]:
    worker = role_build / ("wirehair_wh2_v2_facade_timing_worker_" + role)
    adjacent = role_build / "libwirehair.so.2"
    worker_artifact = artifact_receipt(worker, 256 * 1024 * 1024)
    adjacent_artifact = artifact_receipt(adjacent)
    original_artifact = artifact_receipt(implementation_library)
    if adjacent_artifact["sha256"] != original_artifact["sha256"]:
        fail(role + " original and adjacent role libraries differ")
    source_manifest = source_manifest_entries(harness, "harness")
    source_manifest.extend(source_manifest_entries(implementation, "implementation"))
    source_manifest.sort(key=lambda entry: entry["path"])
    compiler_hash = FROZEN_TOOL_SHA256[CXX]
    screen_link_argv = canonical_link_argv(
        worker_link, role_build, adjacent,
    )
    receipt: Dict[str, Any] = {
        "archive_closure": [],
        "build_argv": list(build_argv),
        "build_type": "Release",
        "compile_argv": list(worker_units[0].argv),
        "compiler_path": str(CXX),
        "compiler_sha256": compiler_hash,
        "configure_argv": list(configure_argv),
        "dt_needed": list(elf["dt_needed"]),
        "dynamic_dependencies": list(elf["dynamic_dependencies"]),
        "generator": "Ninja",
        "git_path": str(GIT),
        "git_sha256": FROZEN_TOOL_SHA256[GIT],
        "harness_git_commit": harness_commit,
        "harness_source_root": str(harness.root),
        "harness_tree_oid": harness.tree_oid,
        "harness_tree_listing_sha256": harness.tree_listing_sha256,
        "implementation_git_commit": implementation_commit,
        "implementation_library_path": str(implementation_library),
        "implementation_library_sha256": original_artifact["sha256"],
        "implementation_source_root": str(implementation.root),
        "implementation_tree_oid": implementation.tree_oid,
        "implementation_tree_listing_sha256": implementation.tree_listing_sha256,
        "link_argv": screen_link_argv,
        "object_closure": object_receipts(worker_units),
        "receipt_preimage_sha256": None,
        "role": role,
        "role_library_path": str(adjacent),
        "role_library_sha256": adjacent_artifact["sha256"],
        "runpath": elf["runpath"],
        "schema": BUILD_RECEIPT_SCHEMA,
        "source_manifest": source_manifest,
        "worker_path": str(worker),
        "worker_sha256": worker_artifact["sha256"],
    }
    receipt["receipt_preimage_sha256"] = sha256_bytes(canonical_bytes(receipt))
    return receipt


def seal_build_authority(config: BuildSealConfig) -> str:
    root_boundary = require_root_mode(BUILD_UNIT_NAME, (0, 1))
    require_runtime_primitives()
    become_subreaper()
    if current_service_pids() != [os.getpid()]:
        fail("build service is not isolated before source/tool admission")
    if LOWER40.fullmatch(config.expected_harness_commit) is None:
        fail("build harness commit differs")
    for path in (
        config.harness_source, config.current_source, config.parent_source,
    ):
        exact_absolute(path, "build input", must_exist=True)
        # Establish the only non-root Git ownership exception before the first
        # source-derived Git command executes.  materialize_snapshot repeats
        # this check immediately before cloning to close later drift.
        snapshot_clone_environment(path)
    require_science_namespaces_absent()
    if BUILD_AUTHORITY_PATH.exists():
        fail("successful build authority already exists")
    toolchain_before = toolchain_receipt()
    linker_inputs_before = linker_input_closure()
    paths = _prepare_build_root(config)
    started_ns = time.monotonic_ns()
    input_receipts_before = {
        "harness": verify_git_tree(
            config.harness_source, config.expected_harness_commit),
        "current": verify_git_tree(
            config.current_source, CURRENT_IMPLEMENTATION_COMMIT),
        "parent": verify_git_tree(
            config.parent_source, PARENT_IMPLEMENTATION_COMMIT),
    }
    snapshots = {
        "harness": materialize_snapshot(
            config.harness_source, paths["harness"],
            config.expected_harness_commit,
        ),
        "current": materialize_snapshot(
            config.current_source, paths["current_source"],
            CURRENT_IMPLEMENTATION_COMMIT,
        ),
        "parent": materialize_snapshot(
            config.parent_source, paths["parent_source"],
            PARENT_IMPLEMENTATION_COMMIT,
        ),
    }
    input_receipts_after = {
        "harness": verify_git_tree(
            config.harness_source, config.expected_harness_commit),
        "current": verify_git_tree(
            config.current_source, CURRENT_IMPLEMENTATION_COMMIT),
        "parent": verify_git_tree(
            config.parent_source, PARENT_IMPLEMENTATION_COMMIT),
    }
    for name in ("harness", "current", "parent"):
        if input_receipts_before[name] != input_receipts_after[name]:
            fail(name + " source changed across fresh snapshot materialization")
    installed_before = installed_launcher_receipt(snapshots["harness"])
    _validate_required_sources(
        snapshots["harness"],
        (
            SOURCE_LAUNCHER_RELATIVE, CONTROLLER_RELATIVE, WORKER_RELATIVE,
            ROLE_CMAKE_RELATIVE, SAMPLER_RELATIVE, *WORKER_SOURCES,
        ),
        "harness snapshot",
    )
    validate_frozen_harness_sources(snapshots["harness"])
    for role in ("current", "parent"):
        _validate_required_sources(
            snapshots[role],
            ("CMakeLists.txt", "include/wirehair/wirehair.h", *IMPLEMENTATION_SOURCES),
            role + " implementation snapshot",
        )

    environment = build_environment(config.build_root)
    commands: List[Dict[str, Any]] = []
    impl_build_argv: Dict[str, List[str]] = {}
    role_build_argv: Dict[str, List[str]] = {}
    role_configure: Dict[str, List[str]] = {}
    implementation_libraries: Dict[str, Path] = {}
    implementation_units: Dict[str, List[CompileUnit]] = {}
    implementation_links: Dict[str, List[str]] = {}
    worker_units: Dict[str, List[CompileUnit]] = {}
    worker_links: Dict[str, List[str]] = {}
    depfiles: Dict[str, Any] = {}
    elf_receipts: Dict[str, Dict[str, Any]] = {}
    receipts: Dict[str, Dict[str, Any]] = {}
    cmake_inputs_before: Dict[str, Dict[str, Any]] = {}

    for role in ("current", "parent"):
        source = snapshots[role].root
        build = paths[role + "_impl_build"]
        configure = implementation_configure_argv(source, build)
        build_argv = [
            str(CMAKE), "--build", str(build), "--target", "wirehair",
            "--parallel", "2",
        ]
        impl_build_argv[role] = build_argv
        commands.append(run_build_command(
            role + "-implementation-configure", configure,
            config.build_root, environment, BUILD_CONFIGURE_SECONDS,
        ))
        cmake_inputs_before[role + "-implementation"] = ninja_input_closure(
            build, "build.ninja",
        )
        commands.append(run_build_command(
            role + "-implementation-build", build_argv,
            config.build_root, environment, BUILD_COMPILE_SECONDS,
        ))
        library = build / "libwirehair.so.2.0.0"
        if not library.is_file() or library.is_symlink():
            fail(role + " implementation library is absent")
        implementation_libraries[role] = library
        expected_sources = [source / relative for relative in IMPLEMENTATION_SOURCES]
        units, dependency_receipts = parse_compile_graph(
            build, "wirehair", expected_sources,
            build / "CMakeFiles/wirehair.dir",
            implementation_compile_prefix(source),
        )
        implementation_units[role] = units
        depfiles[role + "-implementation"] = dependency_receipts
        ordered_objects = [
            next(unit.object_path for unit in units
                 if unit.source == (source / relative))
            for relative in IMPLEMENTATION_SOURCES
        ]
        implementation_links[role] = exact_link_command(
            build, "wirehair", library,
            ordered_objects,
            expected_argv=implementation_link_argv(
                source, build, ordered_objects,
            ),
        )
        validate_target_command_roster(
            build, "wirehair", len(IMPLEMENTATION_SOURCES),
            ([str(CMAKE), "-E", "cmake_symlink_library",
              "libwirehair.so.2.0.0", "libwirehair.so.2", "libwirehair.so"],),
        )
        commands.append(run_build_command(
            role + "-implementation-no-work",
            [str(NINJA), "-n", "wirehair"], build, environment,
            BUILD_NO_WORK_SECONDS,
            expected_stdout=b"ninja: no work to do.\n",
        ))

    for role, commit in (
        ("current", CURRENT_IMPLEMENTATION_COMMIT),
        ("parent", PARENT_IMPLEMENTATION_COMMIT),
    ):
        build = paths[role + "_role_build"]
        configure = role_configure_argv(
            snapshots["harness"].root, snapshots[role].root,
            implementation_libraries[role], build, role,
            config.expected_harness_commit, commit,
        )
        build_argv = [
            str(CMAKE), "--build", str(build), "--target",
            "wirehair_wh2_v2_facade_timing_worker", "--parallel", "2",
        ]
        role_configure[role] = configure
        role_build_argv[role] = build_argv
        commands.append(run_build_command(
            role + "-worker-configure", configure,
            config.build_root, environment, BUILD_CONFIGURE_SECONDS,
        ))
        cmake_inputs_before[role + "-worker"] = ninja_input_closure(
            build, "build.ninja",
        )
        commands.append(run_build_command(
            role + "-worker-build", build_argv,
            config.build_root, environment, BUILD_COMPILE_SECONDS,
        ))
        worker = build / ("wirehair_wh2_v2_facade_timing_worker_" + role)
        commands.append(run_build_command(
            role + "-worker-selftest", [str(worker), "--selftest"],
            build, environment, BUILD_TEST_SECONDS,
        ))
        expected_sources = [
            snapshots["harness"].root / relative for relative in WORKER_SOURCES
        ]
        units, dependency_receipts = parse_compile_graph(
            build, "wirehair_wh2_v2_facade_timing_worker", expected_sources,
            build / "CMakeFiles/wirehair_wh2_v2_facade_timing_worker.dir",
            role_compile_prefix(
                snapshots[role].root, role,
                config.expected_harness_commit, commit,
            ),
        )
        worker_units[role] = units
        depfiles[role + "-worker"] = dependency_receipts
        adjacent = build / "libwirehair.so.2"
        ordered_objects = [
            next(unit.object_path for unit in units
                 if unit.source == (snapshots["harness"].root / relative))
            for relative in WORKER_SOURCES
        ]
        worker_links[role] = exact_link_command(
            build, "wirehair_wh2_v2_facade_timing_worker", worker,
            ordered_objects, required_library=adjacent,
            expected_argv=role_link_argv(build, worker, ordered_objects),
        )
        validate_target_command_roster(
            build, "wirehair_wh2_v2_facade_timing_worker",
            len(WORKER_SOURCES), (),
        )
        elf_receipts[role] = parse_elf(worker, adjacent)
        commands.append(run_build_command(
            role + "-no-work",
            [str(NINJA), "-n", "wirehair_wh2_v2_facade_timing_worker"],
            build, environment, BUILD_NO_WORK_SECONDS,
            expected_stdout=b"ninja: no work to do.\n",
        ))

    # Sealing after every successful command makes later controller checks a
    # meaningful TOCTOU defense: UID1000 can read but cannot rewrite any input,
    # object, worker, library, or receipt.
    for name in (
        "current_impl_build", "parent_impl_build",
        "current_role_build", "parent_role_build",
    ):
        seal_directory_tree(paths[name], allow_symlinks=True)
    for role in ("current", "parent"):
        os.chmod(str(implementation_libraries[role]), 0o444)
        os.chmod(str(paths[role + "_role_build"] / "libwirehair.so.2"), 0o444)

    snapshot_commits = {
        "harness": config.expected_harness_commit,
        "current": CURRENT_IMPLEMENTATION_COMMIT,
        "parent": PARENT_IMPLEMENTATION_COMMIT,
    }
    for name in ("harness", "current", "parent"):
        final_snapshot = verify_git_tree(
            snapshots[name].root, snapshot_commits[name],
            required_uid=CAMPAIGN_UID, required_gid=CAMPAIGN_GID,
            sealed=True,
        )
        if final_snapshot != snapshots[name]:
            fail(name + " sealed snapshot changed during UID1000 build commands")
    final_inputs = {
        "harness": verify_git_tree(
            config.harness_source, config.expected_harness_commit),
        "current": verify_git_tree(
            config.current_source, CURRENT_IMPLEMENTATION_COMMIT),
        "parent": verify_git_tree(
            config.parent_source, PARENT_IMPLEMENTATION_COMMIT),
    }
    if final_inputs != input_receipts_before:
        fail("mutable build inputs changed during fresh build sealing")
    dependency_closure_sha256 = reverify_dependency_closure(depfiles)

    for role, commit in (
        ("current", CURRENT_IMPLEMENTATION_COMMIT),
        ("parent", PARENT_IMPLEMENTATION_COMMIT),
    ):
        receipt = _role_receipt(
            role, config.expected_harness_commit, commit,
            snapshots["harness"], snapshots[role],
            implementation_libraries[role], paths[role + "_role_build"],
            role_configure[role], role_build_argv[role], worker_units[role],
            worker_links[role], elf_receipts[role],
        )
        receipt_path = config.build_root / (role + "-build-receipt.json")
        write_root_file_exclusive(
            receipt_path, canonical_bytes(receipt) + b"\n", mode=0o444,
        )
        receipts[role] = receipt

    # Exercise the frozen screen/health-adapter interface and capture the
    # exact interpreter file-backed mapping roster under the same UID and
    # sealed source tree used by the campaign.  This deliberately fails build
    # publication when the two frozen Python sources are interface-incompatible.
    profile_manifest = health_source_manifest_receipt(
        snapshots["harness"].root, receipts["current"],
    )
    profile_adapter = _receipt_entry(
        receipts["current"], "harness", HEALTH_ADAPTER_RELATIVE,
    )
    profile_output: List[bytes] = []
    profile_argv = [
        str(PYTHON), "-I", "-B", "-c",
        CONTROLLER_MAP_PROFILE_CODE,
        str(snapshots["harness"].root / CONTROLLER_RELATIVE),
        str(snapshots["harness"].root), profile_manifest["sha256"],
        profile_adapter["sha256"], config.expected_harness_commit,
        receipts["current"]["git_sha256"],
    ]
    commands.append(run_build_command(
        "controller-runtime-map-profile", profile_argv,
        snapshots["harness"].root, SEALED_ENVIRONMENT,
        BUILD_TEST_SECONDS, output_sink=profile_output,
    ))
    profile = parse_canonical_document(
        profile_output[0], "controller runtime map profile",
        MAX_COMMAND_OUTPUT_BYTES,
    )
    if (
        set(profile) != {"paths", "schema"}
        or profile.get("schema")
        != "wirehair.wh2.v2-facade-timing-controller-map-profile.v1"
        or type(profile.get("paths")) is not list
        or profile["paths"] != sorted(set(profile["paths"]))
        or any(type(item) is not str or not Path(item).is_absolute()
               for item in profile["paths"])
    ):
        fail("controller runtime map profile output differs")
    controller_runtime_maps = runtime_map_closure(
        tuple(Path(item) for item in profile["paths"]), "controller profile",
    )
    for name in ("harness", "current", "parent"):
        profiled_snapshot = verify_git_tree(
            snapshots[name].root, snapshot_commits[name],
            required_uid=CAMPAIGN_UID, required_gid=CAMPAIGN_GID,
            sealed=True,
        )
        if profiled_snapshot != snapshots[name]:
            fail(name + " sealed snapshot changed during controller profiling")

    for auxiliary in (
        config.build_root / ".build-home", config.build_root / ".build-tmp",
    ):
        seal_directory_tree(auxiliary, allow_symlinks=False)
    for container in (paths["sources"], paths["builds"]):
        os.chown(str(container), 0, 0, follow_symlinks=False)
        os.chmod(str(container), 0o555, follow_symlinks=False)
    os.chmod(str(config.build_root), 0o555)

    installed_after = installed_launcher_receipt(snapshots["harness"])
    if installed_after != installed_before:
        fail("installed launcher changed during build sealing")
    toolchain_after = toolchain_receipt()
    if toolchain_after != toolchain_before:
        fail("toolchain changed during build sealing")
    linker_inputs_after = linker_input_closure()
    if linker_inputs_after != linker_inputs_before:
        fail("linker input closure changed during build sealing")
    cmake_inputs_after: Dict[str, Dict[str, Any]] = {}
    target_inputs: Dict[str, Dict[str, Any]] = {}
    for role in ("current", "parent"):
        for kind, path_name in (
            ("implementation", role + "_impl_build"),
            ("worker", role + "_role_build"),
        ):
            key = role + "-" + kind
            build = paths[path_name]
            after = ninja_input_closure(build, "build.ninja")
            cmake_inputs_after[key] = after
            before_projection = [
                (entry["path"], entry["bytes"], entry["sha256"])
                for entry in cmake_inputs_before[key]["entries"]
            ]
            after_projection = [
                (entry["path"], entry["bytes"], entry["sha256"])
                for entry in after["entries"]
            ]
            if (
                before_projection != after_projection
                or cmake_inputs_before[key]["raw_sha256"] != after["raw_sha256"]
            ):
                fail(key + " CMake input closure changed during its build")
        target_inputs[role + "-implementation"] = ninja_input_closure(
            paths[role + "_impl_build"], "libwirehair.so.2.0.0",
        )
        target_inputs[role + "-worker"] = ninja_input_closure(
            paths[role + "_role_build"],
            "wirehair_wh2_v2_facade_timing_worker_" + role,
        )
    for role in ("current", "parent"):
        expected_raw = canonical_bytes(receipts[role]) + b"\n"
        receipt_path = config.build_root / (role + "-build-receipt.json")
        digest, info = hash_path(
            receipt_path, MAX_DOCUMENT_BYTES, role + " final build receipt",
        )
        if digest != sha256_bytes(expected_raw) or info.st_size != len(expected_raw):
            fail(role + " final build receipt changed")
    reverify_dependency_closure(depfiles)
    require_science_namespaces_absent()
    if current_service_pids() != [os.getpid()]:
        fail("build service is not quiescent before success publication")
    launcher_runtime_maps = self_runtime_map_closure()

    source_authority = {
        name: {
            "commit": receipt.commit,
            "entries": len(receipt.entries),
            "manifest_sha256": receipt.manifest_sha256,
            "root": str(receipt.root),
            "tree_listing_sha256": receipt.tree_listing_sha256,
            "tree_oid": receipt.tree_oid,
        }
        for name, receipt in snapshots.items()
    }
    graph = {
        role: {
            "implementation_compile_argv": [
                list(unit.argv) for unit in implementation_units[role]
            ],
            "implementation_effective_compile_argv": [
                list(unit.effective_argv) for unit in implementation_units[role]
            ],
            "implementation_link_argv": implementation_links[role],
            "worker_compile_argv": [list(unit.argv) for unit in worker_units[role]],
            "worker_effective_compile_argv": [
                list(unit.effective_argv) for unit in worker_units[role]
            ],
            "worker_link_argv": worker_links[role],
        }
        for role in ("current", "parent")
    }
    authority = self_hashed(BUILD_AUTHORITY_SCHEMA, {
        "build_root": str(config.build_root),
        "commands": commands,
        "cmake_input_closure": cmake_inputs_after,
        "completed_monotonic_ns": time.monotonic_ns(),
        "dependency_closure_sha256": dependency_closure_sha256,
        "depfile_closure": depfiles,
        "elf": elf_receipts,
        "expected_current_implementation_commit": CURRENT_IMPLEMENTATION_COMMIT,
        "expected_harness_commit": config.expected_harness_commit,
        "expected_parent_implementation_commit": PARENT_IMPLEMENTATION_COMMIT,
        "frozen_harness_source_sha256": {
            CONTROLLER_RELATIVE: FROZEN_CONTROLLER_SHA256,
            WORKER_RELATIVE: FROZEN_WORKER_SOURCE_SHA256,
            ROLE_CMAKE_RELATIVE: FROZEN_ROLE_CMAKE_SHA256,
            HEALTH_ADAPTER_RELATIVE: FROZEN_HEALTH_ADAPTER_SHA256,
        },
        "graph": graph,
        "installed_launcher": installed_after,
        "linker_input_closure": linker_inputs_after,
        "launcher_runtime_map_closure": launcher_runtime_maps,
        "controller_runtime_map_closure": controller_runtime_maps,
        "root_boundary": root_boundary,
        "role_receipts": {
            role: {
                "path": str(config.build_root / (role + "-build-receipt.json")),
                "sha256": sha256_bytes(canonical_bytes(receipts[role]) + b"\n"),
            }
            for role in ("current", "parent")
        },
        "source_authority": source_authority,
        "started_monotonic_ns": started_ns,
        "status": "sealed",
        "systemd_build_argv": build_systemd_run_argv(config),
        "target_input_closure": target_inputs,
        "toolchain": toolchain_after,
        "trust_boundary": {
            "concurrent_hostile_uid1000_processes_excluded": True,
            "explanation": (
                "The frozen controller requires tracked source files and its "
                "child output namespace to be UID/GID1000. Root-owned 0555 "
                "directories prevent replacement but UID1000 could chmod and "
                "rewrite those regular files. Launch therefore requires an "
                "administrative guarantee that no other hostile UID1000 "
                "process exists for the attempt lifetime."
            ),
            "root_build_graph_is_sole_build_authority": True,
        },
    })
    raw = canonical_bytes(authority) + b"\n"
    published = publish_root_file_noreplace(
        BUILD_AUTHORITY_PATH, raw, mode=0o400,
    )
    if published["sha256"] != sha256_bytes(raw):
        fail("published build authority hash differs")
    return authority["receipt_sha256"]


def validate_self_hash(value: Mapping[str, Any], schema: str,
                       field_name: str = "receipt_sha256") -> str:
    if type(value) is not dict or value.get("schema") != schema:
        fail("self-hashed document schema differs: " + schema)
    digest = value.get(field_name)
    if type(digest) is not str or LOWER64.fullmatch(digest) is None:
        fail("self-hashed document digest differs: " + schema)
    preimage = dict(value)
    preimage[field_name] = None
    if sha256_bytes(canonical_bytes(preimage)) != digest:
        fail("self-hashed document preimage differs: " + schema)
    return digest


def require_artifact(
    path: Path, expected_sha256: str, expected_bytes: int, where: str,
    *, uid: int, gid: int, mode: int, expected_nlink: int = 1,
) -> Dict[str, Any]:
    if path.resolve(strict=True) != path:
        fail(where + " path is not canonical")
    digest, info = hash_path(
        path, max(expected_bytes, 1), where,
        expected_nlink=expected_nlink,
    )
    if (
        digest != expected_sha256 or info.st_size != expected_bytes
        or info.st_uid != uid or info.st_gid != gid
        or stat.S_IMODE(info.st_mode) != mode
    ):
        fail(where + " artifact policy differs")
    return {
        "bytes": info.st_size, "path": str(path), "sha256": digest,
        "stat": stat_receipt(info),
    }


@dataclass(frozen=True)
class RoleExecution:
    role: str
    worker: Path
    worker_sha256: str
    library: Path
    library_sha256: str
    implementation_library: Path
    implementation_library_sha256: str
    build_receipt: Path
    build_receipt_sha256: str
    worker_receipt: Dict[str, Any]
    runtime_map_receipts: Tuple[Dict[str, Any], ...]


@dataclass(frozen=True)
class PreparedExecution:
    root_boundary: Dict[str, Any]
    build_authority: Dict[str, Any]
    build_authority_raw: bytes
    harness_root: Path
    controller: Path
    controller_sha256: str
    sampler_script: Path
    sampler_sha256: str
    health_adapter_sha256: str
    health_manifest_sha256: str
    health_manifest: Dict[str, Any]
    health_module_loader: Dict[str, Any]
    health_source_git_receipt: Dict[str, Any]
    closure_sha256: Dict[str, str]
    roles: Dict[str, RoleExecution]
    python_image: Path
    python_sha256: str


def _receipt_entry(receipt: Mapping[str, Any], authority: str,
                   relative: str) -> Mapping[str, Any]:
    matches = [
        entry for entry in receipt["source_manifest"]
        if entry.get("authority") == authority
        and entry.get("repository_relative_path") == relative
    ]
    if len(matches) != 1:
        fail("build receipt source member differs: {}:{}".format(
            authority, relative))
    return matches[0]


def health_stat_receipt(value: os.stat_result) -> Dict[str, Any]:
    return {
        "device": value.st_dev, "gid": value.st_gid,
        "inode": value.st_ino, "mode": stat.S_IMODE(value.st_mode),
        "mtime_ns": value.st_mtime_ns, "nlink": value.st_nlink,
        "size": value.st_size, "uid": value.st_uid,
    }


def health_source_manifest_receipt(
    harness_root: Path, receipt: Mapping[str, Any],
) -> Dict[str, Any]:
    preimage = bytearray()
    entries: List[Dict[str, Any]] = []
    for relative in HEALTH_SOURCE_PATHS:
        entry = _receipt_entry(receipt, "harness", relative)
        path = harness_root / relative
        live = require_artifact(
            path, entry["sha256"], entry["bytes"],
            "health source " + relative, uid=CAMPAIGN_UID,
            gid=CAMPAIGN_GID, mode=0o444,
        )
        entries.append({
            "bytes": entry["bytes"],
            "git_blob_oid": entry["git_blob_oid"],
            "path": relative, "sha256": entry["sha256"],
            "stat": {
                key: live["stat"][key] for key in (
                    "device", "gid", "inode", "mode", "mtime_ns",
                    "nlink", "size", "uid",
                )
            },
        })
        preimage.extend((entry["sha256"] + "  " + relative + "\n").encode("ascii"))
    return {"entries": entries, "sha256": sha256_bytes(bytes(preimage))}


def health_source_manifest(harness_root: Path,
                           receipt: Mapping[str, Any]) -> str:
    return health_source_manifest_receipt(harness_root, receipt)["sha256"]


def health_module_loader_receipt(
    manifest: Mapping[str, Any],
) -> Dict[str, Any]:
    by_path = {entry["path"]: entry for entry in manifest["entries"]}
    modules = [
        {
            "bytes": by_path[relative]["bytes"], "module": module,
            "path": relative, "sha256": by_path[relative]["sha256"],
        }
        for module, relative in HEALTH_MODULE_SOURCES
    ]
    result: Dict[str, Any] = {
        "dont_write_bytecode": True, "isolated": True,
        "modules": modules, "optimize": 0, "receipt_sha256": None,
        "schema": "wirehair.wh2.direct-systematic-complement-health-source-loader.v2",
    }
    result["receipt_sha256"] = sha256_bytes(canonical_bytes(result))
    return result


def health_source_git_receipt(
    root: Path, commit: str, manifest: Mapping[str, Any],
) -> Dict[str, Any]:
    by_path = {entry["path"]: entry for entry in manifest["entries"]}
    index_flags = git_command(root, ("ls-files", "-v"), 1024 * 1024)
    if (
        not index_flags or not index_flags.endswith(b"\n")
        or any(not line.startswith(b"H ") for line in index_flags.splitlines())
    ):
        fail("health-source Git index flags differ")
    source_blob_oids = [
        {
            "head_oid": by_path[relative]["git_blob_oid"],
            "path": relative,
            "worktree_oid": by_path[relative]["git_blob_oid"],
        }
        for relative in HEALTH_SOURCE_PATHS
    ]
    git_digest, git_info = hash_path(
        GIT, MAX_TRACKED_FILE_BYTES, "health-source Git executable",
    )
    if git_digest != FROZEN_TOOL_SHA256[GIT]:
        fail("health-source Git executable differs")
    return {
        "executable": {
            "path": str(GIT), "sha256": git_digest,
            "stat": health_stat_receipt(git_info),
        },
        "head": commit,
        "source_blob_oids": source_blob_oids,
        "source_blob_roster_sha256": sha256_bytes(
            canonical_bytes(source_blob_oids)),
        "source_roster_sha256": sha256_bytes(b"".join(
            (relative + "\n").encode("ascii")
            for relative in HEALTH_SOURCE_PATHS
        )),
        "tracked_index_flags_sha256": sha256_bytes(index_flags),
        "tracked_status_sha256": sha256_bytes(b""),
        "worktree_root": str(root),
    }


def screen_worktree_receipt(tree: GitTreeReceipt,
                            authority: str) -> Dict[str, Any]:
    """Reproduce the frozen screen's tracked-worktree receipt exactly."""
    directories = {tree.root}
    for entry in tree.entries:
        parent = (tree.root / entry["repository_relative_path"]).parent
        while parent != tree.root:
            if tree.root not in parent.parents:
                fail(authority + " tracked directory escapes its root")
            directories.add(parent)
            parent = parent.parent
    directory_digest = hashlib.sha256()
    for directory in sorted(directories):
        info = os.stat(str(directory), follow_symlinks=False)
        if (
            not stat.S_ISDIR(info.st_mode) or info.st_uid != 0
            or info.st_gid != 0 or stat.S_IMODE(info.st_mode) != 0o555
            or Path(os.path.realpath(str(directory))) != directory
        ):
            fail(authority + " tracked-directory seal differs")
        directory_digest.update(canonical_bytes({
            "path": os.path.relpath(str(directory), str(tree.root)),
            "stat": stat_receipt(info),
        }) + b"\n")
    total = 0
    worktree_digest = hashlib.sha256()
    for entry in tree.entries:
        path = tree.root / entry["repository_relative_path"]
        digest, info = hash_path(
            path, MAX_TRACKED_FILE_BYTES,
            authority + " tracked file " + entry["repository_relative_path"],
        )
        if (
            digest != entry["sha256"] or info.st_size != entry["bytes"]
            or info.st_uid != CAMPAIGN_UID or info.st_gid != CAMPAIGN_GID
            or stat.S_IMODE(info.st_mode) != 0o444
        ):
            fail(authority + " tracked-worktree file differs")
        total += info.st_size
        if total > MAX_TRACKED_TOTAL_BYTES:
            fail(authority + " tracked-worktree byte bound differs")
        worktree_digest.update(canonical_bytes({
            "bytes": info.st_size,
            "git_blob_oid": entry["git_blob_oid"],
            "git_mode": entry["git_mode"],
            "repository_relative_path": entry["repository_relative_path"],
            "sha256": digest,
        }) + b"\n")
    return {
        "bytes": total,
        "entries": len(tree.entries),
        "directory_manifest_sha256": directory_digest.hexdigest(),
        "worktree_manifest_sha256": worktree_digest.hexdigest(),
    }


def screen_repository_receipt(tree: GitTreeReceipt,
                              authority: str) -> Dict[str, Any]:
    return {
        "authority": authority,
        "commit": tree.commit,
        "root": str(tree.root),
        "status_sha256": sha256_bytes(b""),
        "tree_bytes": tree.tree_bytes,
        "tree_oid": tree.tree_oid,
        "tree_sha256": tree.tree_listing_sha256,
        "worktree": screen_worktree_receipt(tree, authority),
    }


def screen_build_closure_digest(
    receipt: Mapping[str, Any], role: str,
    trees: Mapping[str, GitTreeReceipt],
) -> str:
    """Recompute the exact frozen-screen closure digest under root authority."""
    entries: Dict[str, Mapping[str, Any]] = {}
    for roster_name in (
        "source_manifest", "object_closure", "archive_closure",
        "dynamic_dependencies",
    ):
        roster = receipt.get(roster_name)
        if type(roster) is not list:
            fail(role + " screen closure roster differs: " + roster_name)
        for entry in roster:
            if type(entry) is not dict or type(entry.get("path")) is not str:
                fail(role + " screen closure entry differs")
            previous = entries.get(entry["path"])
            if previous is not None and previous != entry:
                fail(role + " screen closure contains conflicting receipts")
            entries[entry["path"]] = entry
    compiler_path = receipt.get("compiler_path")
    if type(compiler_path) is not str:
        fail(role + " screen compiler path differs")
    compiler_info = os.stat(compiler_path, follow_symlinks=False)
    compiler_entry = {
        "bytes": compiler_info.st_size,
        "path": compiler_path,
        "sha256": receipt.get("compiler_sha256"),
    }
    previous = entries.get(compiler_path)
    if previous is not None and previous != compiler_entry:
        fail(role + " screen compiler closure conflicts")
    entries[compiler_path] = compiler_entry
    source_paths = {
        entry["path"]: entry for entry in receipt["source_manifest"]
    }
    object_paths = {entry["path"] for entry in receipt["object_closure"]}
    archive_paths = {entry["path"] for entry in receipt["archive_closure"]}
    dynamic_paths = {
        entry["path"] for entry in receipt["dynamic_dependencies"]
    }
    verified: List[Dict[str, Any]] = []
    for raw_path in sorted(entries):
        entry = entries[raw_path]
        path = Path(raw_path)
        digest, info = hash_path(
            path, 1024 * 1024 * 1024,
            role + " frozen-screen closure " + raw_path,
            expected_nlink=build_receipt_expected_nlink(
                "dynamic_dependencies" if raw_path in dynamic_paths
                else "source_manifest" if raw_path in source_paths
                else "object_closure",
                raw_path, receipt["role_library_path"],
            ),
        )
        if digest != entry.get("sha256") or info.st_size != entry.get("bytes"):
            fail(role + " frozen-screen closure content differs: " + raw_path)
        if raw_path in source_paths:
            source = source_paths[raw_path]
            if (
                info.st_uid != CAMPAIGN_UID or info.st_gid != CAMPAIGN_GID
                or stat.S_IMODE(info.st_mode) != 0o444
            ):
                fail(role + " frozen-screen source seal differs")
            fd = os.open(
                raw_path, os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK
                | getattr(os, "O_NOFOLLOW", 0),
            )
            try:
                if git_blob_oid_fd(fd, info.st_size) != source["git_blob_oid"]:
                    fail(role + " frozen-screen source Git blob differs")
            finally:
                os.close(fd)
        elif raw_path == compiler_path:
            if (
                info.st_uid != 0 or info.st_gid != 0
                or stat.S_IMODE(info.st_mode) != 0o755
            ):
                fail(role + " frozen-screen compiler seal differs")
        elif raw_path in object_paths or raw_path in archive_paths:
            if (
                info.st_uid != 0 or info.st_gid != 0
                or stat.S_IMODE(info.st_mode) != 0o444
            ):
                fail(role + " frozen-screen object/archive seal differs")
        elif raw_path in dynamic_paths:
            if (
                info.st_uid != 0 or info.st_gid != 0
                or stat.S_IMODE(info.st_mode) not in (0o444, 0o555, 0o644, 0o755)
                or raw_path == receipt["role_library_path"]
                and stat.S_IMODE(info.st_mode) != 0o444
            ):
                fail(role + " frozen-screen dynamic-input seal differs")
        verified.append({
            "bytes": info.st_size, "path": raw_path, "sha256": digest,
            "stat": stat_receipt(info),
        })
    git_digest, git_info = hash_path(
        Path(receipt["git_path"]), MAX_TRACKED_FILE_BYTES,
        role + " frozen-screen Git executable",
    )
    if git_digest != receipt["git_sha256"]:
        fail(role + " frozen-screen Git executable differs")
    repositories = [
        screen_repository_receipt(trees["harness"], role + " harness"),
        screen_repository_receipt(
            trees["implementation"], role + " implementation",
        ),
    ]
    return sha256_bytes(canonical_bytes({
        "artifacts": verified,
        "git": {
            "name": "Git executable", "path": receipt["git_path"],
            "sha256": git_digest, "stat": stat_receipt(git_info),
        },
        "repositories": repositories,
    }))


class _StatFs(ctypes.Structure):
    _fields_ = [
        ("f_type", ctypes.c_long), ("f_bsize", ctypes.c_long),
        ("f_blocks", ctypes.c_ulong * 7), ("f_fsid", ctypes.c_int * 2),
        ("f_namelen", ctypes.c_long), ("f_frsize", ctypes.c_long),
        ("f_flags", ctypes.c_long), ("f_spare", ctypes.c_long * 4),
    ]


def require_runtime_primitives() -> None:
    required = {
        "fork": hasattr(os, "fork"),
        "pidfd_open": hasattr(os, "pidfd_open"),
        "pidfd_send_signal": hasattr(signal, "pidfd_send_signal"),
        "sched_getaffinity": hasattr(os, "sched_getaffinity"),
        "sched_setaffinity": hasattr(os, "sched_setaffinity"),
        "renameat2": hasattr(ctypes.CDLL(None), "renameat2"),
    }
    missing = sorted(name for name, available in required.items() if not available)
    if missing:
        fail("required runtime primitives unavailable: " + ",".join(missing))
    result = _StatFs()
    libc = ctypes.CDLL(None, use_errno=True)
    if libc.statfs(os.fsencode(str(CGROUP_ROOT)), ctypes.byref(result)) != 0:
        error = ctypes.get_errno()
        raise OSError(error, os.strerror(error), str(CGROUP_ROOT))
    if result.f_type != CGROUP2_SUPER_MAGIC:
        fail("runtime requires a cgroup-v2 hierarchy")


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


def read_bounded_path(path: Path, maximum: int, where: str) -> bytes:
    if maximum < 0:
        fail(where + " byte bound differs")
    owned: List[int] = []
    try:
        open_registered(
            owned, str(path), os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK,
        )
        fd = owned[0]
        result = bytearray()
        while len(result) <= maximum:
            try:
                block = os.read(fd, min(65536, maximum + 1 - len(result)))
            except BlockingIOError:
                break
            if not block:
                break
            result.extend(block)
        if len(result) > maximum:
            fail(where + " exceeds its byte bound")
        return bytes(result)
    finally:
        for fd in owned:
            os.close(fd)


def process_cgroup_relative(pid: int) -> str:
    raw = read_bounded_path(Path("/proc/{}/cgroup".format(pid)), 65536,
                            "process cgroup")
    try:
        lines = raw.decode("ascii").splitlines()
    except UnicodeDecodeError as exc:
        fail("process cgroup is not ASCII: " + exception_text(exc))
    matches = [line[3:] for line in lines if line.startswith("0::/")]
    if len(matches) != 1 or not matches[0].startswith("/"):
        fail("process unified cgroup differs")
    return matches[0]


@dataclass
class OwnershipSlot:
    """A shared handoff cell that closes Python return/assignment windows."""
    label: str
    value: Any = None

    def adopt(self, value: Any) -> None:
        if value is None or self.value is not None:
            fail(self.label + " ownership slot is already occupied")
        self.value = value

    def owns(self, value: Any) -> bool:
        return self.value is value


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
    created_paths: List[Path] = field(default_factory=list)
    owned_fds: List[int] = field(default_factory=list)

    def close(self) -> None:
        failures: List[str] = []
        members = (
            "experiment_kill_fd", "sampler_kill_fd", "verifier_kill_fd",
            "run_kill_fd",
        )
        for fd in list(dict.fromkeys(self.owned_fds)):
            if fd < 0:
                continue
            aliases = [name for name in members if getattr(self, name) == fd]
            self.owned_fds = [item for item in self.owned_fds if item != fd]
            for name in aliases:
                setattr(self, name, -1)
            try:
                close_fd_once(fd, "cgroup kill authority")
            except BaseException as exc:
                self.owned_fds.append(fd)
                for name in aliases:
                    setattr(self, name, fd)
                failures.append(exception_text(exc))
        if failures:
            fail("cgroup authority FD closure differs: " + " | ".join(failures))


class CgroupTree:
    def __init__(self, service: Path) -> None:
        self.service = service

    @staticmethod
    def write(path: Path, value: str) -> None:
        owned: List[int] = []
        try:
            open_registered(
                owned, str(path), os.O_WRONLY | os.O_CLOEXEC
                | getattr(os, "O_NOFOLLOW", 0),
            )
            fd = owned[0]
            raw = value.encode("ascii")
            if os.write(fd, raw) != len(raw):
                fail("cgroup write was short: " + str(path))
        finally:
            for fd in owned:
                os.close(fd)

    @staticmethod
    def read(path: Path, maximum: int = 65536) -> str:
        raw = read_bounded_path(path, maximum, "cgroup control file")
        try:
            return raw.decode("ascii").strip()
        except UnicodeDecodeError as exc:
            fail("cgroup control file is not ASCII: " + exception_text(exc))

    @staticmethod
    def mkdir(parent: Path, name: str) -> Path:
        path = parent / name
        os.mkdir(str(path), 0o755)
        info = os.stat(str(path), follow_symlinks=False)
        if (
            not stat.S_ISDIR(info.st_mode) or info.st_uid != 0
            or info.st_gid != 0
        ):
            fail("cgroup leaf policy differs: " + str(path))
        return path

    @classmethod
    def configure(cls, path: Path, cpus: str, pids: int) -> None:
        cls.write(path / "cpuset.mems", MEMORY_NODES)
        cls.write(path / "cpuset.cpus", cpus)
        cls.write(path / "pids.max", str(pids))
        if (
            cls.read(path / "cpuset.mems") != MEMORY_NODES
            or cls.read(path / "cpuset.mems.effective") != MEMORY_NODES
            or cls.read(path / "cpuset.cpus") != cpus
            or cls.read(path / "cpuset.cpus.effective") != cpus
            or cls.read(path / "pids.max") != str(pids)
            or cls.read(path / "cgroup.type") != "domain"
        ):
            fail("cgroup configuration readback differs: " + str(path))

    def create(self, owner: Optional[OwnershipSlot] = None) -> CgroupLayout:
        supervisor = self.service / "supervisor"
        run = self.service / "run"
        sampler = run / "sampler"
        experiment = run / "experiment"
        verifier = run / "verifier"
        layout = CgroupLayout(
            service=self.service, supervisor=supervisor, run=run,
            sampler=sampler, experiment=experiment, verifier=verifier,
            experiment_kill_fd=-1, sampler_kill_fd=-1,
            verifier_kill_fd=-1, run_kill_fd=-1,
        )
        ownership_transferred = False
        if owner is not None:
            try:
                owner.adopt(layout)
            except BaseException:
                ownership_transferred = owner.owns(layout)
                if not ownership_transferred:
                    layout.close()
                raise
            ownership_transferred = True
        controllers = set(self.read(self.service / "cgroup.controllers").split())
        if not {"cpuset", "pids"}.issubset(controllers):
            fail("delegated service lacks cpuset/pids controllers")
        self.mkdir(self.service, "supervisor")
        layout.created_paths.append(supervisor)
        self.write(supervisor / "cgroup.procs", str(os.getpid()))
        os.sched_setaffinity(0, {AUTHORITY_CPU})
        if os.sched_getaffinity(0) != {AUTHORITY_CPU}:
            fail("root supervisor affinity differs")
        if self.read(self.service / "cgroup.procs"):
            fail("delegated service root remained populated")
        self.write(self.service / "cgroup.subtree_control", "+cpuset +pids")
        if set(self.read(self.service / "cgroup.subtree_control").split()) != {
            "cpuset", "pids",
        }:
            fail("service subtree controller roster differs")
        self.configure(supervisor, str(AUTHORITY_CPU), 4)
        self.mkdir(self.service, "run")
        layout.created_paths.append(run)
        self.configure(run, "120-123", 36)
        self.write(run / "cgroup.subtree_control", "+cpuset +pids")
        if set(self.read(run / "cgroup.subtree_control").split()) != {
            "cpuset", "pids",
        }:
            fail("run subtree controller roster differs")
        self.mkdir(run, "sampler")
        layout.created_paths.append(sampler)
        self.mkdir(run, "experiment")
        layout.created_paths.append(experiment)
        self.mkdir(run, "verifier")
        layout.created_paths.append(verifier)
        self.configure(sampler, str(SAMPLER_CPU), 3)
        self.configure(experiment, "120-123", 24)
        self.configure(verifier, str(AUTHORITY_CPU), 3)
        for leaf in (sampler, experiment, verifier):
            if "populated 0" not in self.read(leaf / "cgroup.events").splitlines():
                fail("new cgroup leaf is populated: " + str(leaf))
        kill_fds: List[int] = []
        try:
            for member, path in (
                ("experiment_kill_fd", experiment),
                ("sampler_kill_fd", sampler),
                ("verifier_kill_fd", verifier),
                ("run_kill_fd", run),
            ):
                open_registered(
                    layout.owned_fds, str(path / "cgroup.kill"),
                    os.O_WRONLY | os.O_CLOEXEC
                    | getattr(os, "O_NOFOLLOW", 0),
                )
                fd = layout.owned_fds[-1]
                setattr(layout, member, fd)
                kill_fds.append(fd)
            return layout
        except BaseException:
            if not ownership_transferred:
                layout.close()
            raise

    @staticmethod
    def move_pid(path: Path, pid: int, cpus: Iterable[int]) -> None:
        expected = sorted(set(cpus))
        if not expected:
            fail("child affinity roster is empty")
        CgroupTree.write(path / "cgroup.procs", str(pid))
        os.sched_setaffinity(pid, set(expected))
        if sorted(os.sched_getaffinity(pid)) != expected:
            fail("child affinity differs after cgroup placement")
        expected_relative = str(path)[len(str(CGROUP_ROOT)):]
        if process_cgroup_relative(pid) != expected_relative:
            fail("child cgroup placement differs")

    @staticmethod
    def pids(path: Path) -> List[int]:
        raw = CgroupTree.read(path / "cgroup.procs")
        if not raw:
            return []
        tokens = raw.splitlines()
        if any(not token.isdigit() or int(token) <= 0 for token in tokens):
            fail("cgroup PID roster is malformed")
        values = [int(token) for token in tokens]
        if len(values) != len(set(values)):
            fail("cgroup PID roster repeats a process")
        return sorted(values)


def cgroup_descendant_pids(path: Path) -> List[int]:
    result: List[int] = []
    for root, directories, _files in os.walk(
        str(path), topdown=True, followlinks=False,
    ):
        directories[:] = sorted(
            name for name in directories
            if name not in (".", "..") and "/" not in name
        )
        result.extend(CgroupTree.pids(Path(root)))
    return sorted(set(result))


def wait_status_code(status: int) -> int:
    if os.WIFEXITED(status):
        return os.WEXITSTATUS(status)
    if os.WIFSIGNALED(status):
        return -os.WTERMSIG(status)
    fail("child wait status differs")


def pidfd_signal(pidfd: int, signum: int) -> None:
    if pidfd < 0 or not hasattr(signal, "pidfd_send_signal"):
        fail("pidfd signaling is unavailable")
    signal.pidfd_send_signal(pidfd, signum)


def reap_failed_spawn(pid: int, pidfd: int, seconds: float = 2.0) -> None:
    failures: List[str] = []
    if pidfd >= 0:
        try:
            pidfd_signal(pidfd, signal.SIGKILL)
        except ProcessLookupError:
            pass
        except BaseException as exc:
            failures.append("pidfd signal: " + exception_text(exc))
    else:
        # Before startup release this is still our unreaped direct child, so a
        # numeric PID cannot have been reused.
        try:
            waited, _status = os.waitpid(pid, os.WNOHANG)
            if waited == 0:
                os.kill(pid, signal.SIGKILL)
            elif waited == pid:
                return
        except ProcessLookupError:
            pass
        except BaseException as exc:
            failures.append("direct-child signal: " + exception_text(exc))
    deadline = time.monotonic() + seconds
    reaped = False
    while time.monotonic() < deadline:
        try:
            waited, _status = os.waitpid(pid, os.WNOHANG)
        except ChildProcessError:
            reaped = True
            break
        if waited == pid:
            reaped = True
            break
        time.sleep(0.005)
    if pidfd >= 0:
        try:
            os.close(pidfd)
        except OSError as exc:
            failures.append("pidfd close: " + exception_text(exc))
    if not reaped:
        failures.append("failed-spawn child was not reaped")
    if failures:
        fail("failed-spawn cleanup differs: " + " | ".join(failures))


def reap_direct_child_delta(baseline: Iterable[int], where: str) -> List[int]:
    """Recover a forked child whose PID return could not be assigned.

    The launcher is single-threaded during role admission.  A child appearing
    in the direct-child roster after the immediately-pre-fork snapshot is
    therefore the only safe authority available if ``fork()`` completed and
    an asynchronous exception landed before Python stored its return value.
    """
    original = set(baseline)
    children = sorted(set(direct_child_pids()) - original)
    failures: List[str] = []
    for child in children:
        try:
            reap_failed_spawn(child, -1)
        except BaseException as exc:
            failures.append(
                "child {}: {}".format(child, exception_text(exc)))
    try:
        remaining = sorted(set(direct_child_pids()) - original)
        if remaining:
            failures.append("remaining children " + repr(remaining))
    except BaseException as exc:
        failures.append("final roster: " + exception_text(exc))
    if failures:
        fail(where + " direct-child recovery differs: " + " | ".join(failures))
    return children


@dataclass
class ProcessHandle:
    role: str
    pid: int
    pidfd: int
    start_ticks: int
    argv: Tuple[str, ...]
    stdout_fd: int
    stderr_fd: int
    returncode: Optional[int] = None
    communication_started: bool = False
    communication_stdout: bytearray = field(default_factory=bytearray)
    communication_stderr: bytearray = field(default_factory=bytearray)
    communication_inflight_fd: int = -1
    communication_wait_inflight: bool = False
    communication_result: Optional[Tuple[bytes, bytes, int]] = None

    def close(self) -> None:
        for member in ("pidfd", "stdout_fd", "stderr_fd"):
            try:
                close_object_fd(self, member, self.role + " " + member)
            except OSError:
                pass


def exact_environment_bytes(environment: Mapping[str, str]) -> bytes:
    return b"".join(
        (key + "=" + value + "\0").encode("ascii")
        for key, value in environment.items()
    )


def exact_cmdline_bytes(argv: Sequence[str]) -> bytes:
    return b"".join(item.encode("utf-8") + b"\0" for item in argv)


def process_argv(pid: int, where: str) -> List[str]:
    raw = read_bounded_path(
        Path("/proc/{}/cmdline".format(pid)), 1024 * 1024, where + " cmdline",
    )
    if not raw or not raw.endswith(b"\0") or b"\0\0" in raw:
        fail(where + " cmdline is not a canonical NUL vector")
    try:
        values = [item.decode("utf-8") for item in raw[:-1].split(b"\0")]
    except UnicodeDecodeError as exc:
        fail(where + " cmdline is not UTF-8: " + exception_text(exc))
    if not values or any(not value for value in values):
        fail(where + " cmdline has an empty member")
    return values


def process_security_receipt(
    pid: int, expected_argv: Sequence[str], expected_groups: Sequence[int],
    allowed_affinities: Sequence[Sequence[int]], python_info: os.stat_result,
    expected_environment: Mapping[str, str] = SEALED_ENVIRONMENT,
) -> Dict[str, Any]:
    status_raw = read_bounded_path(
        Path("/proc/{}/status".format(pid)), 65536, "process status",
    )
    wanted = {
        b"Uid", b"Gid", b"Groups", b"CapInh", b"CapPrm", b"CapEff",
        b"CapBnd", b"CapAmb", b"NoNewPrivs",
    }
    fields: Dict[bytes, bytes] = {}
    for line in status_raw.splitlines():
        if b":" not in line:
            continue
        key, value = line.split(b":", 1)
        if key in wanted:
            if key in fields:
                fail("process status repeats a security field")
            fields[key] = value.strip()
    if set(fields) != wanted:
        fail("process security status is incomplete")
    zero_caps = (b"0", b"0000000000000000")
    if (
        fields[b"Uid"].split() != [b"1000"] * 4
        or fields[b"Gid"].split() != [b"1000"] * 4
        or [int(item) for item in fields[b"Groups"].split()]
        != list(expected_groups)
        or fields[b"NoNewPrivs"] != b"1"
        or any(fields[key] not in zero_caps for key in (
            b"CapInh", b"CapPrm", b"CapEff", b"CapBnd", b"CapAmb",
        ))
    ):
        fail("process credential/capability boundary differs")
    cmdline = read_bounded_path(
        Path("/proc/{}/cmdline".format(pid)), 1024 * 1024,
        "process cmdline",
    )
    environment = read_bounded_path(
        Path("/proc/{}/environ".format(pid)), 65536,
        "process environment",
    )
    affinity = sorted(os.sched_getaffinity(pid))
    allowed = [sorted(set(value)) for value in allowed_affinities]
    executable = os.stat("/proc/{}/exe".format(pid), follow_symlinks=True)
    if (
        cmdline != exact_cmdline_bytes(expected_argv)
        or environment != exact_environment_bytes(expected_environment)
        or affinity not in allowed or full_stat(executable) != full_stat(python_info)
    ):
        fail("process image/argv/environment/affinity boundary differs")
    return {
        "affinity": affinity,
        "cmdline_sha256": sha256_bytes(cmdline),
        "environ_sha256": sha256_bytes(environment),
        "executable_stat": stat_receipt(executable),
        "gid": CAMPAIGN_GID, "groups": list(expected_groups),
        "no_new_privs": 1, "pid": pid,
        "process_start_ticks": process_start_ticks(pid), "uid": CAMPAIGN_UID,
        "topology": process_topology(pid),
    }


def setpriv_runtime_argv(target: Sequence[str], *, sampler: bool,
                         ready_fd: int, gate_fd: int) -> List[str]:
    if not target or any(type(item) is not str or not item for item in target):
        fail("runtime target argv differs")
    groups = ["--groups", str(I2C_GID)] if sampler else ["--clear-groups"]
    bootstrap = [
        str(PYTHON), "-I", "-S", "-B", "-c", ROLE_SECURITY_BOOTSTRAP,
        str(ready_fd), str(gate_fd),
        json.dumps(list(target), separators=(",", ":")),
    ]
    return [
        str(SETPRIV), "--reuid", str(CAMPAIGN_UID),
        "--regid", str(CAMPAIGN_GID), *groups,
        "--bounding-set=-all", "--inh-caps=-all", "--ambient-caps=-all",
        "--no-new-privs", "--pdeathsig", "SIGKILL", *bootstrap,
    ]


def child_exec_failure(message: str) -> None:
    try:
        os.write(2, (message[:4000] + "\n").encode("ascii", "backslashreplace"))
    finally:
        os._exit(125)


@dataclass
class BlockedRole:
    process: ProcessHandle
    gate_fd: int
    bootstrap_security: Dict[str, Any]
    expected_cpus: Tuple[int, ...]
    expected_groups: Tuple[int, ...]
    python_info: os.stat_result
    released: bool = False
    may_have_released: bool = False
    final_security: Optional[Dict[str, Any]] = None
    cancellation_roster: Optional[List[int]] = None
    cancellation_result: Optional[Dict[str, Any]] = None

    def release(self, poll_callback: Optional[Any] = None,
                release_callback: Optional[Any] = None,
                ) -> Tuple[ProcessHandle, Dict[str, Any]]:
        if self.released or self.gate_fd < 0:
            fail(self.process.role + " blocked role was already resolved")
        if poll_callback is not None:
            poll_callback()
        # Publish the irreversible state through the owned BlockedRole before
        # the first gate write.  The caller can therefore recover this fact
        # even if a Python line event lands after write(2) but before return.
        self.may_have_released = True
        if release_callback is not None:
            release_callback()
        try:
            if os.write(self.gate_fd, b"R") != 1:
                fail(self.process.role + " security release was short")
        finally:
            close_object_fd(
                self, "gate_fd", self.process.role + " release gate",
            )
        self.released = True
        if poll_callback is not None:
            poll_callback()
        deadline = time.monotonic() + 2.0
        last_error = "final exec was not observed"
        allowed = (self.expected_cpus,)
        while time.monotonic() < deadline:
            if poll_callback is not None:
                poll_callback()
            try:
                receipt = process_security_receipt(
                    self.process.pid, self.process.argv, self.expected_groups,
                    allowed, self.python_info,
                )
                if receipt["process_start_ticks"] != self.process.start_ticks:
                    fail(self.process.role + " start identity changed across exec")
                security = dict(self.bootstrap_security)
                security["final_process"] = receipt
                self.final_security = security
                return self.process, security
            except (LaunchError, FileNotFoundError, ProcessLookupError) as exc:
                last_error = exception_text(exc)
            if process_exited(self.process):
                break
            time.sleep(0.005)
        fail(self.process.role + " final exec admission failed: " + last_error)

    def cancel(self) -> None:
        if self.gate_fd >= 0:
            close_object_fd(
                self, "gate_fd", self.process.role + " cancellation gate",
            )


@dataclass
class BlockedSpawnPreamble:
    pid: int
    out_r: int
    err_r: int
    launch_w: int
    ready_r: int
    gate_w: int
    pidfd_baseline: List[int]
    setpriv_command: Tuple[str, ...]

    def cleanup(self) -> None:
        failures: List[str] = []
        # Closing the placement and target gates first guarantees a blocked
        # child cannot advance while the numeric PID is being contained.
        for member in ("launch_w", "gate_w", "ready_r", "out_r", "err_r"):
            try:
                close_object_fd(self, member, "blocked preamble " + member)
            except BaseException as exc:
                failures.append(member + ": " + exception_text(exc))
        try:
            reap_failed_spawn(self.pid, -1)
        except BaseException as exc:
            failures.append("child: " + exception_text(exc))
        try:
            close_fd_delta(self.pidfd_baseline, ())
        except BaseException as exc:
            failures.append("descriptor delta: " + exception_text(exc))
        if failures:
            fail("blocked preamble cleanup differs: " + " | ".join(failures))


def fork_blocked_preamble(
    role: str, cgroup: Path, expected_cpus: Tuple[int, ...],
    target: Sequence[str], *, sampler: bool,
    fault_injector: Optional[Any] = None,
    owner: Optional[OwnershipSlot] = None,
) -> BlockedSpawnPreamble:
    """Fork a blocked role with cleanup authority spanning every handoff."""
    out_r = out_w = err_r = err_w = -1
    launch_r = launch_w = ready_r = ready_w = -1
    gate_r = gate_w = -1
    pid = -1
    fork_parent_pid = os.getpid()
    fork_child_baseline: Optional[Tuple[int, ...]] = None
    baseline = open_fd_roster()
    try:
        out_r, out_w = os.pipe2(os.O_CLOEXEC)
        err_r, err_w = os.pipe2(os.O_CLOEXEC)
        launch_r, launch_w = os.pipe2(os.O_CLOEXEC)
        ready_r, ready_w = os.pipe2(os.O_CLOEXEC)
        gate_r, gate_w = os.pipe2(os.O_CLOEXEC)
        os.set_blocking(out_r, False)
        os.set_blocking(err_r, False)
        setpriv_command = setpriv_runtime_argv(
            target, sampler=sampler, ready_fd=ready_w, gate_fd=gate_r,
        )
        pidfd_baseline = open_fd_roster()
        fork_child_baseline = tuple(direct_child_pids())
        pid = os.fork()
        if pid == 0:
            try:
                for fd in (out_r, err_r, launch_w, ready_r, gate_w):
                    os.close(fd)
                os.setsid()
                os.dup2(out_w, 1, inheritable=True)
                os.dup2(err_w, 2, inheritable=True)
                if out_w not in (1, 2):
                    os.close(out_w)
                if err_w not in (1, 2):
                    os.close(err_w)
                CgroupTree.move_pid(cgroup, os.getpid(), expected_cpus)
                if os.read(launch_r, 1) != b"R":
                    fail(role + " placement release differs")
                os.close(launch_r)
                os.set_inheritable(ready_w, True)
                os.set_inheritable(gate_r, True)
                os.execve(
                    str(SETPRIV), list(setpriv_command),
                    dict(SEALED_ENVIRONMENT),
                )
            except BaseException as exc:
                child_exec_failure(role + " exec: " + exception_text(exc))
        if fault_injector is not None:
            fault_injector("preamble-post-fork", role)
        # Parent ownership is published in the returned object only after all
        # child-only ends are closed.  Any exception before return is handled
        # by this same scope, including a line event immediately after fork.
        for fd in (out_w, err_w, launch_r, ready_w, gate_r):
            close_fd_once(fd, role + " parent child-end")
        preamble = BlockedSpawnPreamble(
            pid=pid, out_r=out_r, err_r=err_r, launch_w=launch_w,
            ready_r=ready_r, gate_w=gate_w,
            pidfd_baseline=pidfd_baseline,
            setpriv_command=tuple(setpriv_command),
        )
        if owner is not None:
            try:
                owner.adopt(preamble)
            except BaseException:
                if owner.owns(preamble):
                    raise
                raise
        return preamble
    except BaseException as exc:
        if os.getpid() != fork_parent_pid or pid == 0:
            child_exec_failure(
                role + " post-fork branch: " + exception_text(exc),
            )
        if (
            owner is not None and "preamble" in locals()
            and owner.owns(preamble)
        ):
            raise
        cleanup_failures: List[str] = []
        # Gate closure precedes numeric-PID cleanup.  The child either remains
        # blocked or observes EOF and exits; reap_failed_spawn is the backstop.
        for fd in (launch_w, gate_w, ready_r, out_r, err_r,
                   out_w, err_w, launch_r, ready_w, gate_r):
            if fd < 0:
                continue
            try:
                close_fd_once(fd, role + " failed preamble")
            except BaseException as cleanup_exc:
                cleanup_failures.append(exception_text(cleanup_exc))
        if pid > 0:
            try:
                reap_failed_spawn(pid, -1)
            except BaseException as cleanup_exc:
                cleanup_failures.append(exception_text(cleanup_exc))
        elif fork_child_baseline is not None:
            try:
                reap_direct_child_delta(
                    fork_child_baseline, role + " failed preamble",
                )
            except BaseException as cleanup_exc:
                cleanup_failures.append(exception_text(cleanup_exc))
        try:
            close_fd_delta(baseline, ())
        except BaseException as cleanup_exc:
            cleanup_failures.append(exception_text(cleanup_exc))
        if cleanup_failures:
            raise LaunchError(
                exception_text(exc) + " | blocked preamble cleanup: "
                + " | ".join(cleanup_failures)
            ) from exc
        raise


def spawn_blocked_role(
    role: str, cgroup: Path, cpus: Iterable[int], target: Sequence[str],
    *, sampler: bool, python_info: os.stat_result,
    fault_injector: Optional[Any] = None,
    owner: Optional[OwnershipSlot] = None,
    poll_callback: Optional[Any] = None,
) -> BlockedRole:
    expected_cpus = tuple(sorted(set(cpus)))
    groups = (I2C_GID,) if sampler else ()
    preamble: Optional[BlockedSpawnPreamble] = None
    preamble_owner = OwnershipSlot(role + " spawn preamble")
    try:
        fork_blocked_preamble(
            role, cgroup, expected_cpus, target, sampler=sampler,
            fault_injector=fault_injector,
            owner=preamble_owner,
        )
        preamble = preamble_owner.value
        if not isinstance(preamble, BlockedSpawnPreamble):
            fail(role + " spawn preamble adoption failed")
        pid = preamble.pid
        out_r = preamble.out_r
        err_r = preamble.err_r
        launch_w = preamble.launch_w
        ready_r = preamble.ready_r
        gate_w = preamble.gate_w
        setpriv_command = list(preamble.setpriv_command)
        bootstrap_index = setpriv_command.index(str(PYTHON), 1)
        bootstrap_argv = setpriv_command[bootstrap_index:]
        # The preamble closed five child-only pipe ends.  Snapshot only now so
        # a pidfd that reuses one of those numbers is not mistaken for a
        # pre-existing descriptor after a return-then-raise pidfd_open.
        pidfd_baseline = open_fd_roster()
        pidfd = -1
        launch_open = True
        gate_open = True
        blocked: Optional[BlockedRole] = None
        ownership_transferred = False
    except BaseException:
        if preamble is None and isinstance(
            preamble_owner.value, BlockedSpawnPreamble,
        ):
            preamble = preamble_owner.value
        if preamble is not None:
            preamble.cleanup()
        raise
    try:
        if fault_injector is not None:
            fault_injector("before-pidfd", role)
        pidfd = os.pidfd_open(pid, 0)
        if fault_injector is not None:
            fault_injector("after-pidfd", role)
        start = process_start_ticks(pid)
        if fault_injector is not None:
            fault_injector("after-start-ticks", role)
        CgroupTree.move_pid(cgroup, pid, expected_cpus)
        if os.write(launch_w, b"R") != 1:
            fail(role + " placement release was short")
        os.close(launch_w)
        launch_open = False
        ready_deadline = time.monotonic() + 5.0
        readable: List[int] = []
        while time.monotonic() < ready_deadline:
            if poll_callback is not None:
                poll_callback()
            readable, _, _ = select.select(
                [ready_r], [], [], min(0.05, max(
                    0.0, ready_deadline - time.monotonic(),
                )),
            )
            if readable:
                break
        raw = os.read(ready_r, 65537) if readable else b""
        if not raw or len(raw) > 65536:
            fail(role + " did not publish bounded bootstrap security")
        security = parse_canonical_document(raw, role + " bootstrap security", 65536)
        expected_caps = {
            "CapAmb": "0000000000000000", "CapBnd": "0000000000000000",
            "CapEff": "0000000000000000", "CapInh": "0000000000000000",
            "CapPrm": "0000000000000000", "NoNewPrivs": "1",
        }
        if security != {
            "affinity": list(expected_cpus), "caps": expected_caps,
            "environment": dict(SEALED_ENVIRONMENT),
            "gid": [CAMPAIGN_GID] * 3, "groups": list(groups),
            "uid": [CAMPAIGN_UID] * 3,
        }:
            fail(role + " bootstrap security differs")
        bootstrap_receipt = process_security_receipt(
            pid, bootstrap_argv, groups, (expected_cpus,), python_info,
        )
        if (
            bootstrap_receipt["process_start_ticks"] != start
            or process_start_ticks(pid) != start
            or sorted(os.sched_getaffinity(pid)) != list(expected_cpus)
            or process_cgroup_relative(pid)
            != str(cgroup)[len(str(CGROUP_ROOT)):]
        ):
            fail(role + " bootstrap process identity differs")
        security["bootstrap_process"] = bootstrap_receipt
        if poll_callback is not None:
            poll_callback()
        if fault_injector is not None:
            fault_injector("blocked-admitted", role)
        handle = ProcessHandle(
            role=role, pid=pid, pidfd=pidfd, start_ticks=start,
            argv=tuple(target), stdout_fd=out_r, stderr_fd=err_r,
        )
        blocked = BlockedRole(
            process=handle, gate_fd=gate_w,
            bootstrap_security=security, expected_cpus=expected_cpus,
            expected_groups=groups, python_info=python_info,
        )
        if owner is not None:
            try:
                owner.adopt(blocked)
            except BaseException:
                ownership_transferred = owner.owns(blocked)
                if ownership_transferred:
                    gate_open = False
                raise
            ownership_transferred = True
        gate_open = False
        return blocked
    except BaseException as exc:
        if ownership_transferred or (
            owner is not None and blocked is not None and owner.owns(blocked)
        ):
            raise
        if launch_open:
            try:
                os.close(launch_w)
            except OSError:
                pass
        if blocked is not None and blocked.gate_fd >= 0:
            try:
                blocked.cancel()
            except OSError:
                pass
        elif gate_open:
            try:
                os.close(gate_w)
            except OSError:
                pass
        try:
            reap_failed_spawn(pid, pidfd)
            pidfd = -1
        except BaseException as cleanup_exc:
            exc = LaunchError(
                exception_text(exc) + " | failed-spawn cleanup: "
                + exception_text(cleanup_exc)
            )
        # pidfd_open can itself complete and then be interrupted before its
        # return value is assigned.  Runtime admission is single-threaded, so
        # the exact FD delta identifies and closes that otherwise-anonymous
        # leaked pidfd without touching the pre-existing pipe roster.
        try:
            close_fd_delta(pidfd_baseline, ())
        except BaseException as cleanup_exc:
            exc = LaunchError(
                exception_text(exc) + " | pidfd-delta cleanup: "
                + exception_text(cleanup_exc)
            )
        for fd in (out_r, err_r):
            try:
                os.close(fd)
            except OSError:
                pass
        raise exc
    finally:
        caller_owns = (
            owner is not None and blocked is not None and owner.owns(blocked)
        )
        for fd in (ready_r, launch_w if launch_open else -1,
                   gate_w if gate_open and not caller_owns else -1):
            if fd >= 0:
                try:
                    os.close(fd)
                except OSError:
                    pass


def spawn_gated_role(
    role: str, cgroup: Path, cpus: Iterable[int], target: Sequence[str],
    *, sampler: bool, python_info: os.stat_result,
    fault_injector: Optional[Any] = None,
) -> Tuple[ProcessHandle, Dict[str, Any]]:
    blocked = spawn_blocked_role(
        role, cgroup, cpus, target, sampler=sampler,
        python_info=python_info, fault_injector=fault_injector,
    )
    try:
        return blocked.release()
    except BaseException:
        blocked.cancel()
        try:
            pidfd_signal(blocked.process.pidfd, signal.SIGKILL)
        except ProcessLookupError:
            pass
        try:
            wait_process(blocked.process, time.monotonic() + 2.0)
        except BaseException:
            pass
        blocked.process.close()
        raise


def process_exited(process: ProcessHandle) -> bool:
    poller = select.poll()
    poller.register(process.pidfd, select.POLLIN)
    return bool(poller.poll(0))


def wait_process(process: ProcessHandle, deadline: float) -> int:
    while process.returncode is None:
        try:
            waited, status_value = os.waitpid(process.pid, os.WNOHANG)
        except ChildProcessError:
            fail(process.role + " child was reaped without accounting")
        if waited == process.pid:
            process.returncode = wait_status_code(status_value)
            break
        if time.monotonic() >= deadline:
            fail(process.role + " reap deadline expired")
        time.sleep(0.005)
    return process.returncode


def communicate_process(
    process: ProcessHandle, deadline: float, stdout_limit: int,
    stderr_limit: int, poll_callback: Optional[Any] = None,
) -> Tuple[bytes, bytes, int]:
    if process.communication_result is not None:
        output, error, code = process.communication_result
        if len(output) > stdout_limit or len(error) > stderr_limit:
            fail(process.role + " retained output exceeds the requested bound")
        return process.communication_result
    if process.communication_inflight_fd >= 0:
        fail(process.role + " output read result is ambiguous")
    if process.communication_wait_inflight:
        fail(process.role + " wait result is ambiguous")
    process.communication_started = True
    output = process.communication_stdout
    error = process.communication_stderr
    streams: Dict[int, Tuple[bytearray, int]] = {
        fd: value for fd, value in (
            (process.stdout_fd, (output, stdout_limit)),
            (process.stderr_fd, (error, stderr_limit)),
        ) if fd >= 0
    }
    poller = select.poll()
    for fd in streams:
        poller.register(fd, select.POLLIN | select.POLLHUP | select.POLLERR)
    while process.returncode is None or streams:
        if time.monotonic() >= deadline:
            fail(process.role + " deadline expired")
        if poll_callback is not None:
            poll_callback()
        for fd, _event in poller.poll(20):
            if fd not in streams:
                continue
            target, limit = streams[fd]
            while True:
                process.communication_inflight_fd = fd
                try:
                    block = os.read(fd, 65536)
                except BlockingIOError:
                    process.communication_inflight_fd = -1
                    break
                if not block:
                    process.communication_inflight_fd = -1
                    poller.unregister(fd)
                    streams.pop(fd)
                    member = (
                        "stdout_fd" if process.stdout_fd == fd
                        else "stderr_fd" if process.stderr_fd == fd
                        else ""
                    )
                    if not member:
                        fail(process.role + " output descriptor ownership differs")
                    close_object_fd(
                        process, member, process.role + " " + member,
                    )
                    break
                target.extend(block)
                if len(target) > limit:
                    fail(process.role + " output exceeded its bound")
                process.communication_inflight_fd = -1
        if process.returncode is None:
            process.communication_wait_inflight = True
            waited, status_value = os.waitpid(process.pid, os.WNOHANG)
            process.communication_wait_inflight = False
            if waited == process.pid:
                process.returncode = wait_status_code(status_value)
    result = (bytes(output), bytes(error), process.returncode)
    process.communication_result = result
    return result


def close_fds_except(retained: Iterable[int]) -> None:
    keep = set(retained)
    for fd in open_fd_roster():
        if fd > 2 and fd not in keep:
            try:
                os.close(fd)
            except OSError:
                pass


def append_guardian_record(fd: int, directory_fd: int,
                           value: Mapping[str, Any]) -> None:
    raw = canonical_bytes(dict(value)) + b"\n"
    if len(raw) > 4096:
        fail("guardian journal append bound differs")
    before = os.fstat(fd)
    if (
        not stat.S_ISREG(before.st_mode) or before.st_size < 0
        or before.st_size > 64 * 1024 - len(raw)
        or os.lseek(fd, 0, os.SEEK_CUR) != before.st_size
    ):
        fail("guardian journal append position differs")
    existing = os.pread(fd, before.st_size, 0)
    if len(existing) != before.st_size:
        fail("guardian journal pre-append read was short")
    already_appended = bool(existing.endswith(raw))
    if already_appended:
        record_offset = before.st_size - len(raw)
    else:
        line_offset = existing.rfind(b"\n") + 1
        partial = existing[line_offset:]
        if partial and not raw.startswith(partial):
            fail("guardian journal partial record differs")
        record_offset = line_offset if partial else before.st_size

        def resume_append() -> None:
            while True:
                current = os.fstat(fd)
                if (
                    current.st_size < record_offset
                    or current.st_size > record_offset + len(raw)
                ):
                    fail("guardian journal append extent differs")
                progress = current.st_size - record_offset
                observed = os.pread(fd, progress, record_offset)
                if len(observed) != progress or observed != raw[:progress]:
                    fail("guardian journal partial append differs")
                if progress == len(raw):
                    if os.lseek(fd, 0, os.SEEK_END) != current.st_size:
                        fail("guardian journal completed position differs")
                    return
                if os.lseek(fd, current.st_size, os.SEEK_SET) != current.st_size:
                    fail("guardian journal resume position differs")
                try:
                    written = os.write(fd, raw[progress:])
                except BaseException:
                    after_error = os.fstat(fd)
                    error_progress = after_error.st_size - record_offset
                    if (
                        error_progress <= progress or error_progress > len(raw)
                        or os.pread(fd, error_progress, record_offset)
                        != raw[:error_progress]
                    ):
                        raise
                    continue
                if written <= 0 or written > len(raw) - progress:
                    fail("guardian journal append made no valid progress")

        # One Python line-event or syscall-return exception after irreversible
        # progress is recovered by inspecting the exact intended prefix and
        # resuming it.  A short write is ordinary progress, never a license to
        # begin a contradictory failure record on the same physical line.
        try:
            resume_append()
        except BaseException:
            resume_append()
    after = os.fstat(fd)
    if (
        after.st_size != record_offset + len(raw)
        or record_offset < 0
        or os.pread(fd, len(raw), record_offset) != raw
    ):
        fail("guardian journal append readback differs")
    # fsync is safely repeatable.  A return-boundary exception is retried
    # here so the guardian never converts an already-written intended record
    # into a contradictory third failure record.
    for target in (fd, directory_fd):
        try:
            os.fsync(target)
        except BaseException:
            os.fsync(target)


def guardian_record_is_last(fd: int, value: Mapping[str, Any]) -> bool:
    raw = canonical_bytes(dict(value)) + b"\n"
    info = os.fstat(fd)
    if not stat.S_ISREG(info.st_mode) or info.st_size < len(raw):
        return False
    return os.pread(fd, len(raw), info.st_size - len(raw)) == raw


def root_process_status_fields(pid: int) -> Dict[str, str]:
    raw = read_bounded_path(
        Path("/proc/{}/status".format(pid)), 65536, "root process status",
    )
    wanted = {
        "Uid", "Gid", "Groups", "CapInh", "CapPrm", "CapEff", "CapBnd",
        "CapAmb", "NoNewPrivs",
    }
    result: Dict[str, str] = {}
    for line in raw.decode("ascii", "strict").splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        if key in wanted:
            if key in result:
                fail("root process status repeats " + key)
            result[key] = value.strip()
    if (
        set(result) != wanted
        or result["Uid"].split() != ["0"] * 4
        or result["Gid"].split() != ["0"] * 4
        or result["Groups"].split()
        or result["NoNewPrivs"] not in ("0", "1")
        or any(re.fullmatch(r"[0-9a-f]{16}", result[name]) is None
               for name in ("CapInh", "CapPrm", "CapEff", "CapBnd", "CapAmb"))
    ):
        fail("root process credential/capability status differs")
    return result


def root_guardian_security_receipt(
    pid: int, expected_status: Mapping[str, str], expected_cmdline: bytes,
    expected_environ: bytes, expected_executable: os.stat_result,
    expected_map_closure: Mapping[str, Any], layout: CgroupLayout,
) -> Dict[str, Any]:
    status_fields = root_process_status_fields(pid)
    cmdline = read_bounded_path(
        Path("/proc/{}/cmdline".format(pid)), 1024 * 1024,
        "root guardian cmdline",
    )
    environ = read_bounded_path(
        Path("/proc/{}/environ".format(pid)), 65536,
        "root guardian environ",
    )
    executable = os.stat("/proc/{}/exe".format(pid), follow_symlinks=True)
    _, mapped = mapped_file_paths(pid)
    expected_paths = tuple(
        Path(entry["path"]) for entry in expected_map_closure["entries"]
    )
    topology = process_topology(pid)
    if (
        status_fields != dict(expected_status)
        or cmdline != expected_cmdline or environ != expected_environ
        or full_stat(executable) != full_stat(expected_executable)
        or tuple(mapped) != expected_paths
        or topology != {
            "parent_pid": os.getpid(), "process_group": pid, "session": pid,
        }
        or sorted(os.sched_getaffinity(pid)) != [AUTHORITY_CPU]
        or process_cgroup_relative(pid)
        != str(layout.supervisor)[len(str(CGROUP_ROOT)):]
    ):
        fail("independent root guardian security boundary differs")
    return {
        "affinity": [AUTHORITY_CPU], "cmdline_sha256": sha256_bytes(cmdline),
        "environ_sha256": sha256_bytes(environ),
        "executable_stat": stat_receipt(executable),
        "mapped_path_manifest_sha256": expected_map_closure[
            "path_manifest_sha256"],
        "pid": pid, "process_start_ticks": process_start_ticks(pid),
        "status": status_fields, "topology": topology,
    }


def guardian_child(control_fd: int, ready_fd: int, run_kill_fd: int,
                   event_fd: int, journal_directory_fd: int,
                   supervisor: Path, deadline_ns: int) -> None:
    try:
        close_fds_except((
            control_fd, ready_fd, run_kill_fd, event_fd,
            journal_directory_fd,
        ))
        os.setsid()
        CgroupTree.move_pid(supervisor, os.getpid(), (AUTHORITY_CPU,))
        armed_ns = time.monotonic_ns()
        armed = {
            "armed_monotonic_ns": armed_ns,
            "deadline_monotonic_ns": deadline_ns,
            "event": "armed", "pid": os.getpid(),
            "schema": "wirehair.wh2.v2-facade-timing-root-guardian-event.v1",
        }
        append_guardian_record(event_fd, journal_directory_fd, armed)
        ready = canonical_bytes({
            "armed_monotonic_ns": armed_ns,
            "deadline_monotonic_ns": deadline_ns,
            "pid": os.getpid(),
            "schema": "wirehair.wh2.v2-facade-timing-root-guardian-ready.v1",
        }) + b"\n"
        if armed_ns >= deadline_ns or os.write(ready_fd, ready) != len(ready):
            raise LaunchError("guardian could not publish an early ready record")
        os.close(ready_fd)
        poller = select.poll()
        poller.register(control_fd, select.POLLIN | select.POLLHUP | select.POLLERR)
        terminal_event = "deadline-fired"
        terminal_observed_ns = deadline_ns
        while True:
            remaining_ms = max(
                0, (deadline_ns - time.monotonic_ns() + 999999) // 1000000,
            )
            events = poller.poll(min(remaining_ms, 1000))
            if events:
                terminal_observed_ns = time.monotonic_ns()
                token = os.read(control_fd, 2)
                if (
                    token and set(token) == {ord("D")}
                    and terminal_observed_ns < deadline_ns
                ):
                    terminal_event = "completed"
                elif token and set(token) == {ord("D")}:
                    terminal_event = "control-late"
                elif token:
                    terminal_event = "control-invalid"
                else:
                    terminal_event = "control-lost"
                break
            if time.monotonic_ns() >= deadline_ns:
                terminal_observed_ns = time.monotonic_ns()
                break
        terminal = {
            "deadline_monotonic_ns": deadline_ns,
            "event": terminal_event,
            "observed_monotonic_ns": terminal_observed_ns,
            "pid": os.getpid(),
            "schema": "wirehair.wh2.v2-facade-timing-root-guardian-event.v1",
        }
        append_guardian_record(event_fd, journal_directory_fd, terminal)
        os.fchmod(event_fd, 0o400)
        os.fsync(event_fd)
        os.fsync(journal_directory_fd)
        if terminal_event == "completed":
            os.close(control_fd)
            os._exit(0)
        if os.write(run_kill_fd, b"1") != 1:
            raise LaunchError("guardian cgroup.kill write was short")
        os._exit(124)
    except BaseException as exc:
        # An asynchronous exception can land in the caller immediately after
        # the intended terminal append has completed.  Recover that exact
        # durable record and preserve its outcome; never append a conflicting
        # third `failure` terminal.
        try:
            intended = locals().get("terminal")
            if (
                type(intended) is dict
                and guardian_record_is_last(event_fd, intended)
            ):
                os.fchmod(event_fd, 0o400)
                os.fsync(event_fd)
                os.fsync(journal_directory_fd)
                if intended.get("event") == "completed":
                    os._exit(0)
                try:
                    os.write(run_kill_fd, b"1")
                finally:
                    os._exit(124)
        except BaseException:
            pass
        try:
            if stat.S_IMODE(os.fstat(event_fd).st_mode) == 0o600:
                append_guardian_record(event_fd, journal_directory_fd, {
                    "deadline_monotonic_ns": deadline_ns,
                    "error": exception_text(exc), "event": "failure",
                    "observed_monotonic_ns": time.monotonic_ns(),
                    "pid": os.getpid(),
                    "schema":
                        "wirehair.wh2.v2-facade-timing-root-guardian-event.v1",
                })
                os.fchmod(event_fd, 0o400)
                os.fsync(event_fd)
                os.fsync(journal_directory_fd)
        except BaseException:
            pass
        try:
            os.write(run_kill_fd, b"1")
        except BaseException:
            pass
        os._exit(125)


@dataclass
class GuardianHandle:
    process: ProcessHandle
    control_fd: int
    event_fd: int
    ready: Dict[str, Any]
    completed: bool = False
    completion_may_have_been_sent: bool = False
    terminal_receipt: Optional[Dict[str, Any]] = None

    def poll(self) -> None:
        if process_exited(self.process):
            fail("independent root guardian exited before completion")
        if process_start_ticks(self.process.pid) != self.process.start_ticks:
            fail("independent root guardian identity changed")

    def complete(self) -> Dict[str, Any]:
        if self.completed:
            if self.terminal_receipt is None:
                fail("resolved independent root guardian lacks its receipt")
            return self.terminal_receipt
        if not process_exited(self.process):
            if self.control_fd < 0:
                fail("independent root guardian lacks completion authority")
            self.poll()
            # Completion tokens are deliberately idempotent: the guardian
            # accepts any nonempty all-`D` read.  Therefore a retry after a
            # write-return or Python handoff exception can safely send D
            # again.  This Boolean records irreversible intent for abort(),
            # but never suppresses the retryable write itself.
            self.completion_may_have_been_sent = True
            if os.write(self.control_fd, b"D") != 1:
                fail("independent root guardian completion was short")
        if self.control_fd >= 0:
            close_object_fd(
                self, "control_fd", "independent guardian control",
            )
        code = wait_process(self.process, time.monotonic() + 3.0)
        if code != 0:
            fail("independent root guardian completion exit differs")
        events = validate_guardian_events(
            self.event_fd, self.process.pid,
            self.ready["deadline_monotonic_ns"], "completed", code,
        )
        receipt = {
            **self.ready, "completed": True, "events": events,
            "returncode": code,
        }
        # Store the exact replayable receipt before publishing the resolved
        # Boolean.  If an asynchronous exception lands at either assignment,
        # a retry can revalidate the already-reaped child or return this
        # receipt without writing the irreversible completion token again.
        self.terminal_receipt = receipt
        self.completed = True
        return receipt

    def close(self) -> None:
        if self.control_fd >= 0:
            close_object_fd(
                self, "control_fd", "independent guardian control",
            )


def validate_guardian_events(fd: int, pid: int, deadline_ns: int,
                             terminal_event: Optional[str],
                             exit_code: Optional[int] = None,
                             ) -> List[Dict[str, Any]]:
    raw = read_held_file(fd, 64 * 1024, "root guardian journal")
    records = raw.splitlines(keepends=True)
    if len(records) != 2:
        fail("root guardian journal record count differs")
    values = [
        parse_canonical_document(record, "root guardian journal", 4096)
        for record in records
    ]
    schema = "wirehair.wh2.v2-facade-timing-root-guardian-event.v1"
    armed_keys = {
        "armed_monotonic_ns", "deadline_monotonic_ns", "event", "pid",
        "schema",
    }
    ordinary_terminal_keys = {
        "deadline_monotonic_ns", "event", "observed_monotonic_ns", "pid",
        "schema",
    }
    failure_terminal_keys = ordinary_terminal_keys | {"error"}
    event = values[1].get("event")
    permitted_events = {
        "completed", "deadline-fired", "control-late", "control-invalid",
        "control-lost", "failure",
    }
    expected_exit = {
        "completed": 0,
        "deadline-fired": 124,
        "control-late": 124,
        "control-invalid": 124,
        "control-lost": 124,
        "failure": 125,
    }.get(event)
    observed = values[1].get("observed_monotonic_ns")
    if (
        set(values[0]) != armed_keys
        or set(values[1])
        != (failure_terminal_keys if event == "failure"
            else ordinary_terminal_keys)
        or any(value.get("schema") != schema or value.get("pid") != pid
               or value.get("deadline_monotonic_ns") != deadline_ns
               for value in values)
        or values[0].get("event") != "armed"
        or type(values[0].get("armed_monotonic_ns")) is not int
        or values[0]["armed_monotonic_ns"] < 0
        or values[0]["armed_monotonic_ns"] >= deadline_ns
        or event not in permitted_events
        or type(observed) is not int or observed < 0
        or observed < values[0]["armed_monotonic_ns"]
        or (event == "completed" and observed >= deadline_ns)
        or (event in ("deadline-fired", "control-late")
            and observed < deadline_ns)
        or (event == "failure"
            and (type(values[1].get("error")) is not str
                 or not values[1]["error"]
                 or len(values[1]["error"]) > 4096))
        or terminal_event is not None and event != terminal_event
        or exit_code is not None and exit_code != expected_exit
    ):
        fail("root guardian journal contract differs")
    return values


@dataclass
class GuardianSpawnPreamble:
    pid: int
    control_w: int
    ready_r: int
    pidfd_baseline: List[int]

    def cleanup(self) -> None:
        failures: List[str] = []
        for member in ("control_w", "ready_r"):
            try:
                close_object_fd(self, member, "guardian preamble " + member)
            except BaseException as exc:
                failures.append(member + ": " + exception_text(exc))
        try:
            reap_failed_spawn(self.pid, -1)
        except BaseException as exc:
            failures.append("child: " + exception_text(exc))
        try:
            close_fd_delta(self.pidfd_baseline, ())
        except BaseException as exc:
            failures.append("descriptor delta: " + exception_text(exc))
        if failures:
            fail("guardian preamble cleanup differs: " + " | ".join(failures))


def fork_guardian_preamble(
    layout: CgroupLayout, deadline_ns: int, event_fd: int,
    journal_directory_fd: int, *, fault_injector: Optional[Any] = None,
    owner: Optional[OwnershipSlot] = None,
) -> GuardianSpawnPreamble:
    control_r = control_w = ready_r = ready_w = -1
    pid = -1
    fork_parent_pid = os.getpid()
    fork_child_baseline: Optional[Tuple[int, ...]] = None
    baseline = open_fd_roster()
    try:
        control_r, control_w = os.pipe2(os.O_CLOEXEC)
        ready_r, ready_w = os.pipe2(os.O_CLOEXEC)
        pidfd_baseline = open_fd_roster()
        fork_child_baseline = tuple(direct_child_pids())
        pid = os.fork()
        if pid == 0:
            try:
                os.close(control_w)
                os.close(ready_r)
                guardian_child(
                    control_r, ready_w, layout.run_kill_fd,
                    event_fd, journal_directory_fd, layout.supervisor,
                    deadline_ns,
                )
                os._exit(125)
            except BaseException as exc:
                child_exec_failure(
                    "root guardian child: " + exception_text(exc),
                )
        if fault_injector is not None:
            fault_injector("guardian-preamble-post-fork", "guardian")
        close_fd_once(control_r, "guardian parent control read")
        close_fd_once(ready_w, "guardian parent ready write")
        preamble = GuardianSpawnPreamble(
            pid=pid, control_w=control_w, ready_r=ready_r,
            pidfd_baseline=pidfd_baseline,
        )
        if owner is not None:
            try:
                owner.adopt(preamble)
            except BaseException:
                if owner.owns(preamble):
                    raise
                raise
        return preamble
    except BaseException as exc:
        if os.getpid() != fork_parent_pid or pid == 0:
            child_exec_failure(
                "root guardian post-fork branch: " + exception_text(exc),
            )
        if (
            owner is not None and "preamble" in locals()
            and owner.owns(preamble)
        ):
            raise
        cleanup_failures: List[str] = []
        for fd in (control_w, ready_r, control_r, ready_w):
            if fd < 0:
                continue
            try:
                close_fd_once(fd, "failed guardian preamble")
            except BaseException as cleanup_exc:
                cleanup_failures.append(exception_text(cleanup_exc))
        if pid > 0:
            try:
                reap_failed_spawn(pid, -1)
            except BaseException as cleanup_exc:
                cleanup_failures.append(exception_text(cleanup_exc))
        elif fork_child_baseline is not None:
            try:
                reap_direct_child_delta(
                    fork_child_baseline, "failed guardian preamble",
                )
            except BaseException as cleanup_exc:
                cleanup_failures.append(exception_text(cleanup_exc))
        try:
            close_fd_delta(baseline, ())
        except BaseException as cleanup_exc:
            cleanup_failures.append(exception_text(cleanup_exc))
        if cleanup_failures:
            raise LaunchError(
                exception_text(exc) + " | guardian preamble cleanup: "
                + " | ".join(cleanup_failures)
            ) from exc
        raise


def spawn_guardian(layout: CgroupLayout, deadline_ns: int,
                   event_fd: int, journal_directory_fd: int,
                   fault_injector: Optional[Any] = None,
                   owner: Optional[OwnershipSlot] = None) -> GuardianHandle:
    event_info = os.fstat(event_fd)
    if (
        not stat.S_ISREG(event_info.st_mode) or event_info.st_uid != 0
        or event_info.st_gid != 0 or event_info.st_nlink != 1
        or event_info.st_size != 0 or stat.S_IMODE(event_info.st_mode) != 0o600
    ):
        fail("root guardian journal reservation differs")
    parent_status = root_process_status_fields(os.getpid())
    parent_cmdline = read_bounded_path(
        Path("/proc/self/cmdline"), 1024 * 1024, "root supervisor cmdline",
    )
    parent_environ = read_bounded_path(
        Path("/proc/self/environ"), 65536, "root supervisor environ",
    )
    parent_executable = os.stat("/proc/self/exe", follow_symlinks=True)
    parent_map_closure = self_runtime_map_closure()
    preamble: Optional[GuardianSpawnPreamble] = None
    preamble_owner = OwnershipSlot("guardian spawn preamble")
    try:
        fork_guardian_preamble(
            layout, deadline_ns, event_fd, journal_directory_fd,
            fault_injector=fault_injector, owner=preamble_owner,
        )
        preamble = preamble_owner.value
        if not isinstance(preamble, GuardianSpawnPreamble):
            fail("guardian spawn preamble adoption failed")
        pid = preamble.pid
        control_w = preamble.control_w
        ready_r = preamble.ready_r
        # The preamble has already closed its child-only pipe ends.  Snapshot
        # now so a pidfd that reuses one of those numbers cannot hide in the
        # stale pre-fork descriptor roster after a return-then-raise fault.
        baseline = open_fd_roster()
        pidfd = -1
        ownership_transferred = False
    except BaseException:
        if preamble is None and isinstance(
            preamble_owner.value, GuardianSpawnPreamble,
        ):
            preamble = preamble_owner.value
        if preamble is not None:
            preamble.cleanup()
        raise
    try:
        if fault_injector is not None:
            fault_injector("guardian-before-pidfd", "guardian")
        pidfd = os.pidfd_open(pid, 0)
        start = process_start_ticks(pid)
        CgroupTree.move_pid(layout.supervisor, pid, (AUTHORITY_CPU,))
        readable, _, _ = select.select([ready_r], [], [], 5.0)
        raw = os.read(ready_r, 8193) if readable else b""
        if not raw or len(raw) > 8192:
            fail("independent root guardian did not publish readiness")
        ready = parse_canonical_document(raw, "root guardian ready", 8192)
        if (
            ready.get("schema")
            != "wirehair.wh2.v2-facade-timing-root-guardian-ready.v1"
            or ready.get("pid") != pid
            or ready.get("deadline_monotonic_ns") != deadline_ns
            or type(ready.get("armed_monotonic_ns")) is not int
            or ready["armed_monotonic_ns"] >= deadline_ns
            or process_start_ticks(pid) != start
            or sorted(os.sched_getaffinity(pid)) != [AUTHORITY_CPU]
            or process_cgroup_relative(pid)
            != str(layout.supervisor)[len(str(CGROUP_ROOT)):]
        ):
            fail("independent root guardian admission differs")
        security = root_guardian_security_receipt(
            pid, parent_status, parent_cmdline, parent_environ,
            parent_executable, parent_map_closure, layout,
        )
        process = ProcessHandle(
            role="guardian", pid=pid, pidfd=pidfd, start_ticks=start,
            argv=("<forked-root-guardian>",), stdout_fd=-1, stderr_fd=-1,
        )
        # The ready pipe is only an admission accelerator.  The durable first
        # journal record is independently held and revalidated here.
        current = read_held_file(
            event_fd, 64 * 1024, "root guardian armed journal",
        )
        armed_records = current.splitlines(keepends=True)
        if len(armed_records) != 1:
            fail("root guardian armed journal count differs")
        armed = parse_canonical_document(
            armed_records[0], "root guardian armed journal", 4096,
        )
        if (
            armed.get("event") != "armed" or armed.get("pid") != pid
            or armed.get("deadline_monotonic_ns") != deadline_ns
            or armed.get("armed_monotonic_ns") != ready["armed_monotonic_ns"]
        ):
            fail("root guardian pipe/durable readiness differs")
        handle = GuardianHandle(process, control_w, event_fd, ready)
        handle.ready["security"] = security
        if owner is not None:
            try:
                owner.adopt(handle)
            except BaseException:
                ownership_transferred = owner.owns(handle)
                raise
            ownership_transferred = True
        return handle
    except BaseException as exc:
        if ownership_transferred or (
            owner is not None and 'handle' in locals() and owner.owns(handle)
        ):
            raise
        try:
            os.close(control_w)
        except OSError:
            pass
        try:
            reap_failed_spawn(pid, pidfd)
            pidfd = -1
        except BaseException as cleanup_exc:
            exc = LaunchError(
                exception_text(exc) + " | guardian cleanup: "
                + exception_text(cleanup_exc)
            )
        try:
            close_fd_delta(baseline, ())
        except BaseException as cleanup_exc:
            exc = LaunchError(
                exception_text(exc) + " | guardian FD cleanup: "
                + exception_text(cleanup_exc)
            )
        raise exc
    finally:
        try:
            os.close(ready_r)
        except OSError:
            pass


@dataclass
class CampaignLock:
    parent_fd: int
    fd: int
    binding: Tuple[int, int]

    @classmethod
    def acquire(cls) -> "CampaignLock":
        exact_absolute(LOCK_PATH, "campaign lock path")
        parent_owned: List[int] = []
        lock_owned: List[int] = []
        absent_before = False
        try:
            open_registered(
                parent_owned, str(LOCK_PATH.parent),
                os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC
                | getattr(os, "O_NOFOLLOW", 0),
            )
            parent_fd = parent_owned[0]
            parent = os.fstat(parent_fd)
            named_parent = os.stat(
                str(LOCK_PATH.parent), follow_symlinks=False,
            )
            if (
                not stat.S_ISDIR(parent.st_mode)
                or full_stat(parent) != full_stat(named_parent)
                or parent.st_uid != 0 or parent.st_gid != 0
                or stat.S_IMODE(parent.st_mode) != 0o1777
            ):
                fail("campaign lock parent policy differs")
            try:
                os.stat(
                    LOCK_PATH.name, dir_fd=parent_fd, follow_symlinks=False,
                )
            except FileNotFoundError:
                absent_before = True
            open_registered(
                lock_owned, LOCK_PATH.name,
                os.O_RDWR | os.O_CREAT | os.O_CLOEXEC
                | getattr(os, "O_NOFOLLOW", 0),
                0o600, dir_fd=parent_fd,
            )
            fd = lock_owned[0]
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
            info = os.fstat(fd)
            named = os.stat(
                LOCK_PATH.name, dir_fd=parent_fd, follow_symlinks=False,
            )
            if (
                not stat.S_ISREG(info.st_mode) or full_stat(info) != full_stat(named)
                or info.st_uid != 0 or info.st_gid != 0 or info.st_nlink != 1
                or stat.S_IMODE(info.st_mode) != 0o600 or info.st_size != 0
            ):
                fail("campaign lock policy differs")
            return cls(parent_fd, fd, (info.st_dev, info.st_ino))
        except BaseException:
            for opened in reversed(lock_owned):
                try:
                    os.close(opened)
                except OSError:
                    pass
            if absent_before and parent_owned:
                try:
                    created = os.stat(
                        LOCK_PATH.name, dir_fd=parent_owned[0],
                        follow_symlinks=False,
                    )
                except FileNotFoundError:
                    pass
                else:
                    if (
                        not stat.S_ISREG(created.st_mode) or created.st_uid != 0
                        or created.st_gid != 0 or created.st_nlink != 1
                        or stat.S_IMODE(created.st_mode) != 0o600
                        or created.st_size != 0
                    ):
                        fail("failed campaign-lock reservation changed type")
                    os.unlink(LOCK_PATH.name, dir_fd=parent_owned[0])
                    os.fsync(parent_owned[0])
            for opened in reversed(parent_owned):
                try:
                    os.close(opened)
                except OSError:
                    pass
            raise

    def verify(self) -> None:
        held = os.fstat(self.fd)
        named = os.stat(
            LOCK_PATH.name, dir_fd=self.parent_fd, follow_symlinks=False,
        )
        if (
            not same_inode(held, named)
            or (held.st_dev, held.st_ino) != self.binding
            or held.st_uid != 0 or held.st_gid != 0 or held.st_nlink != 1
            or stat.S_IMODE(held.st_mode) != 0o600 or held.st_size != 0
        ):
            fail("campaign lock binding changed")

    def close(self) -> None:
        for member in ("fd", "parent_fd"):
            close_object_fd(self, member, "campaign lock " + member)


def require_execute_namespaces_absent(parent_fd: int) -> None:
    parent = os.fstat(parent_fd)
    named_parent = os.stat("/var/tmp", follow_symlinks=False)
    if (
        not stat.S_ISDIR(parent.st_mode) or full_stat(parent) != full_stat(named_parent)
        or parent.st_uid != 0 or parent.st_gid != 0
        or stat.S_IMODE(parent.st_mode) != 0o1777
    ):
        fail("/var/tmp execution authority differs")
    for path in (ATTEMPT_DIR, CONTROLLER_PARENT, SAMPLER_DIR):
        try:
            os.stat(path.name, dir_fd=parent_fd, follow_symlinks=False)
        except FileNotFoundError:
            continue
        fail("execution namespace already exists before attempt: " + str(path))


@dataclass
class BoundDirectory:
    path: Path
    parent_fd: int
    uid: int
    gid: int
    fd: int = -1
    binding: Optional[Tuple[int, int]] = None
    state: str = "absent"

    def create(self) -> None:
        if self.state == "bound" and self.fd >= 0 and self.binding is not None:
            # A prior call may have completed the fixed-name rename and then
            # lost a directory/parent fsync acknowledgement.  The object still
            # owns the exact inode, so retry only the idempotent durability and
            # binding checks; never allocate a second namespace.
            self.verify(0o700)
            os.fsync(self.fd)
            os.fsync(self.parent_fd)
            self.verify(0o700)
            return
        if self.state != "absent" or self.fd >= 0:
            fail("owned directory creation state differs: " + str(self.path))
        temporary = AttemptJournal._private_name("owned-directory")
        try:
            os.stat(self.path.name, dir_fd=self.parent_fd, follow_symlinks=False)
        except FileNotFoundError:
            pass
        else:
            fail("owned directory fixed name already exists: " + str(self.path))
        try:
            os.stat(temporary, dir_fd=self.parent_fd, follow_symlinks=False)
        except FileNotFoundError:
            pass
        else:
            fail("owned directory private name already exists")
        try:
            os.mkdir(temporary, 0o700, dir_fd=self.parent_fd)
            self.state = "private"
            self.fd = os.open(
                temporary,
                os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC
                | getattr(os, "O_NOFOLLOW", 0), dir_fd=self.parent_fd,
            )
            os.fchown(self.fd, self.uid, self.gid)
            os.fchmod(self.fd, 0o700)
            held = os.fstat(self.fd)
            self.binding = (held.st_dev, held.st_ino)
            rename_noreplace_at(self.parent_fd, temporary, self.path.name)
            self.state = "bound"
            self.verify(0o700)
            os.fsync(self.fd)
            os.fsync(self.parent_fd)
        except BaseException as original:
            if self.fd >= 0:
                try:
                    held = os.fstat(self.fd)
                    named = os.stat(
                        self.path.name, dir_fd=self.parent_fd,
                        follow_symlinks=False,
                    )
                except (FileNotFoundError, OSError):
                    pass
                else:
                    if same_inode(held, named):
                        self.binding = (held.st_dev, held.st_ino)
                        self.state = "bound"
            if self.state != "bound":
                cleanup_failures: List[str] = []
                private_removed = False
                try:
                    private = os.stat(
                        temporary, dir_fd=self.parent_fd,
                        follow_symlinks=False,
                    )
                except FileNotFoundError:
                    private = None
                    private_removed = True
                except BaseException as cleanup_exc:
                    private = None
                    cleanup_failures.append(
                        "private inspection: " + exception_text(cleanup_exc))
                if private is not None:
                    try:
                        retained = (self.parent_fd, self.fd)
                        close_untracked_inode_fds(private, retained)
                        os.rmdir(temporary, dir_fd=self.parent_fd)
                        os.fsync(self.parent_fd)
                        private_removed = True
                    except BaseException as cleanup_exc:
                        cleanup_failures.append(
                            "private removal: " + exception_text(cleanup_exc))
                if self.fd >= 0:
                    try:
                        close_object_fd(
                            self, "fd",
                            "failed owned directory " + str(self.path),
                        )
                    except BaseException as cleanup_exc:
                        cleanup_failures.append(
                            "descriptor closure: "
                            + exception_text(cleanup_exc))
                if private_removed and self.fd < 0:
                    self.binding = None
                    self.state = "absent"
                if cleanup_failures:
                    raise LaunchError(
                        exception_text(original)
                        + " | owned-directory create cleanup: "
                        + " | ".join(cleanup_failures)
                    ) from original
            raise

    def verify(self, mode: int) -> os.stat_result:
        if self.fd < 0 or self.binding is None or self.state not in ("bound", "sealed"):
            fail("owned directory lacks a held binding: " + str(self.path))
        held = os.fstat(self.fd)
        named = os.stat(
            self.path.name, dir_fd=self.parent_fd, follow_symlinks=False,
        )
        if (
            not stat.S_ISDIR(held.st_mode) or not same_inode(held, named)
            or (held.st_dev, held.st_ino) != self.binding
            or held.st_uid != self.uid or held.st_gid != self.gid
            or stat.S_IMODE(held.st_mode) != mode
        ):
            fail("owned directory binding differs: " + str(self.path))
        return held

    def seal(self, expected_names: Sequence[str]) -> Dict[str, Any]:
        if self.state == "sealed":
            info = self.verify(0o500)
            if sorted(os.listdir(self.fd)) != sorted(expected_names):
                fail("owned directory sealed roster differs: " + str(self.path))
            os.fsync(self.fd)
            os.fsync(self.parent_fd)
            info = self.verify(0o500)
            return stat_receipt(info)
        self.verify(0o700)
        if sorted(os.listdir(self.fd)) != sorted(expected_names):
            fail("owned directory roster differs before seal: " + str(self.path))
        try:
            os.fchmod(self.fd, 0o500)
            os.fsync(self.fd)
            os.fsync(self.parent_fd)
            self.state = "sealed"
            info = self.verify(0o500)
            if sorted(os.listdir(self.fd)) != sorted(expected_names):
                fail("owned directory roster differs after seal: " + str(self.path))
            return stat_receipt(info)
        except BaseException as exc:
            # A failed seal is not authority.  Restore the live owner-only mode
            # and its durable binding so cleanup/invalid publication can
            # proceed.  If restoration itself fails, preserve the original
            # exception plus the cleanup failure and never claim `sealed`.
            self.state = "bound"
            try:
                os.fchmod(self.fd, 0o700)
                os.fsync(self.fd)
                os.fsync(self.parent_fd)
                self.verify(0o700)
            except BaseException as cleanup_exc:
                raise LaunchError(
                    exception_text(exc) + " | owned-directory seal cleanup: "
                    + exception_text(cleanup_exc)
                ) from exc
            raise

    def close(self) -> None:
        if self.fd >= 0:
            close_object_fd(self, "fd", "owned directory " + str(self.path))


def file_binding(fd: int, maximum: int, where: str,
                 *, expected_nlink: Optional[int] = 1,
                 allow_empty: bool = False) -> Dict[str, Any]:
    before = os.fstat(fd)
    if (
        not stat.S_ISREG(before.st_mode)
        or (expected_nlink is not None and before.st_nlink != expected_nlink)
        or before.st_size < (0 if allow_empty else 1)
        or before.st_size > maximum
    ):
        fail(where + " held-file policy differs")
    digest = hash_fd(fd, before.st_size)
    after = os.fstat(fd)
    if full_stat(before) != full_stat(after):
        fail(where + " held file changed while hashing")
    return {
        "bytes": before.st_size, "device": before.st_dev,
        "gid": before.st_gid, "inode": before.st_ino,
        "mode": stat.S_IMODE(before.st_mode), "nlink": before.st_nlink,
        "sha256": digest, "uid": before.st_uid,
    }


def read_held_file(fd: int, maximum: int, where: str,
                   *, expected_nlink: Optional[int] = 1,
                   allow_empty: bool = False) -> bytes:
    binding = file_binding(
        fd, maximum, where, expected_nlink=expected_nlink,
        allow_empty=allow_empty,
    )
    result = bytearray()
    while len(result) < binding["bytes"]:
        block = os.pread(
            fd, min(1024 * 1024, binding["bytes"] - len(result)), len(result),
        )
        if not block:
            fail(where + " held read was short")
        result.extend(block)
    if file_binding(
        fd, maximum, where, expected_nlink=expected_nlink,
        allow_empty=allow_empty,
    ) != binding:
        fail(where + " changed across held read")
    return bytes(result)


def sampler_argv(prepared: PreparedExecution) -> List[str]:
    return [
        str(PYTHON), "-I", "-S", "-B", str(prepared.sampler_script),
        "--csv", str(SAMPLER_CSV),
        "--pid-file", str(SAMPLER_PID_FILE),
        "--validation-jsonl", str(SAMPLER_VALIDATION),
        "--receipt", str(SAMPLER_RECEIPT),
        "--expected-source-sha256", prepared.sampler_sha256,
        "--expected-output-owner-uid", str(CAMPAIGN_UID),
        "--interval", "1.0", "--dimm-attempts", "5",
        "--dimm-retry-delay", "0.01",
    ]


@dataclass
class SamplerAuthority:
    process: ProcessHandle
    held: Dict[str, int]
    admission: Dict[str, Any]

    def close(self) -> None:
        failures: List[str] = []
        for name in list(self.held):
            try:
                close_dict_fd(self.held, name, "sampler held " + name)
            except BaseException as exc:
                failures.append(name + ": " + exception_text(exc))
        if failures:
            fail("sampler descriptor closure differs: "
                 + " | ".join(failures))


def close_sampler_file_map(result: Dict[str, int], where: str) -> None:
    failures: List[str] = []
    for name in list(result):
        try:
            close_dict_fd(result, name, where + " " + name)
        except BaseException as exc:
            failures.append(name + ": " + exception_text(exc))
    if failures:
        fail(where + " descriptor closure differs: " + " | ".join(failures))


def open_sampler_files(
    process: ProcessHandle, result: Optional[Dict[str, int]] = None,
) -> Dict[str, int]:
    # A caller-supplied mapping is the persistent owner across the helper's
    # return boundary.  Without it, an exception after return but before the
    # caller's STORE_FAST loses all four held descriptors.
    if result is None:
        result = {}
    if result:
        fail("sampler admission descriptor owner was not empty")
    owned: List[int] = []
    try:
        for name, path in (
            ("csv", SAMPLER_CSV), ("pid", SAMPLER_PID_FILE),
            ("validation", SAMPLER_VALIDATION), ("receipt", SAMPLER_RECEIPT),
        ):
            open_registered(
                owned, str(path),
                os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK
                | getattr(os, "O_NOFOLLOW", 0),
            )
            result[name] = owned[-1]
            owned.pop()
        expected_modes = {
            "csv": 0o600, "pid": 0o444,
            "validation": 0o444, "receipt": 0o444,
        }
        for name, fd in result.items():
            info = os.fstat(fd)
            named = os.stat(
                str({
                    "csv": SAMPLER_CSV, "pid": SAMPLER_PID_FILE,
                    "validation": SAMPLER_VALIDATION,
                    "receipt": SAMPLER_RECEIPT,
                }[name]), follow_symlinks=False,
            )
            if (
                not stat.S_ISREG(info.st_mode) or not same_inode(info, named)
                or info.st_uid != CAMPAIGN_UID or info.st_gid != CAMPAIGN_GID
                or info.st_nlink != 1
                or stat.S_IMODE(info.st_mode) != expected_modes[name]
            ):
                fail("sampler admission artifact policy differs: " + name)
        if os.pread(result["pid"], 64, 0) != (
            str(process.pid) + "\n"
        ).encode("ascii"):
            fail("sampler PID file differs")
        if (
            os.fstat(result["csv"]).st_size == 0
            or os.fstat(result["validation"]).st_size == 0
            or os.fstat(result["receipt"]).st_size != 0
        ):
            raise FileNotFoundError("sampler streams are not initialized")
        return result
    except BaseException as exc:
        result_fds = set(result.values())
        cleanup_failures: List[str] = []
        try:
            close_sampler_file_map(result, "failed sampler admission")
        except BaseException as cleanup_exc:
            cleanup_failures.append(exception_text(cleanup_exc))
        for fd in reversed(owned):
            if fd in result_fds:
                continue
            try:
                close_fd_once(fd, "failed sampler unnamed descriptor")
            except BaseException as cleanup_exc:
                cleanup_failures.append(exception_text(cleanup_exc))
        if cleanup_failures:
            raise LaunchError(
                exception_text(exc) + " | sampler open cleanup: "
                + " | ".join(cleanup_failures)
            ) from exc
        raise


def last_canonical_jsonl(fd: int, maximum: int, where: str) -> Dict[str, Any]:
    info = os.fstat(fd)
    if info.st_size < 1 or info.st_size > maximum:
        fail(where + " stream size differs")
    raw = os.pread(fd, info.st_size, 0)
    if len(raw) != info.st_size or not raw.endswith(b"\n"):
        fail(where + " stream framing differs")
    line = raw.rstrip(b"\n").split(b"\n")[-1] + b"\n"
    return parse_canonical_document(line, where, maximum)


def live_file_identity(fd: int, where: str, *, expected_mode: int,
                       allow_empty: bool = False) -> Dict[str, Any]:
    before = os.fstat(fd)
    if (
        not stat.S_ISREG(before.st_mode) or before.st_nlink != 1
        or before.st_uid != CAMPAIGN_UID or before.st_gid != CAMPAIGN_GID
        or stat.S_IMODE(before.st_mode) != expected_mode
        or before.st_size < (0 if allow_empty else 1)
    ):
        fail(where + " live-file policy differs")
    after = os.fstat(fd)
    stable_identity = (
        before.st_dev, before.st_ino, before.st_mode, before.st_uid,
        before.st_gid, before.st_nlink,
    )
    if stable_identity != (
        after.st_dev, after.st_ino, after.st_mode, after.st_uid,
        after.st_gid, after.st_nlink,
    ):
        fail(where + " live-file identity changed")
    return {
        "bytes_observed": after.st_size, "device": after.st_dev,
        "gid": after.st_gid, "inode": after.st_ino,
        "mode": stat.S_IMODE(after.st_mode), "nlink": after.st_nlink,
        "uid": after.st_uid,
    }


def stable_live_file_snapshot(fd: int, maximum: int, where: str, *,
                              expected_mode: int,
                              attempts: int = 4) -> Tuple[bytes, Dict[str, Any]]:
    last = "changed during snapshot"
    for _ in range(attempts):
        before = os.fstat(fd)
        if (
            not stat.S_ISREG(before.st_mode) or before.st_nlink != 1
            or before.st_uid != CAMPAIGN_UID or before.st_gid != CAMPAIGN_GID
            or stat.S_IMODE(before.st_mode) != expected_mode
            or before.st_size < 1 or before.st_size > maximum
        ):
            fail(where + " live stream policy differs")
        raw = os.pread(fd, before.st_size, 0)
        after = os.fstat(fd)
        identity = (
            before.st_dev, before.st_ino, before.st_mode, before.st_uid,
            before.st_gid, before.st_nlink,
        )
        if (
            len(raw) == before.st_size
            and identity == (
                after.st_dev, after.st_ino, after.st_mode, after.st_uid,
                after.st_gid, after.st_nlink,
            )
            and before.st_size == after.st_size
            and before.st_mtime_ns == after.st_mtime_ns
        ):
            return raw, {
                "bytes": len(raw), "device": after.st_dev,
                "gid": after.st_gid, "inode": after.st_ino,
                "mode": stat.S_IMODE(after.st_mode), "nlink": after.st_nlink,
                "sha256": sha256_bytes(raw), "stat": stat_receipt(after),
                "uid": after.st_uid,
            }
        last = "grew during snapshot"
        time.sleep(0.005)
    raise FileNotFoundError(where + " " + last)


def tail_canonical_jsonl(fd: int, where: str,
                         maximum_line: int = 256 * 1024) -> Dict[str, Any]:
    before = os.fstat(fd)
    if before.st_size < 1:
        raise FileNotFoundError(where + " has no records")
    start = max(0, before.st_size - maximum_line - 1)
    raw = os.pread(fd, before.st_size - start, start)
    after = os.fstat(fd)
    if (
        (before.st_dev, before.st_ino, before.st_mode, before.st_uid,
         before.st_gid, before.st_nlink)
        != (after.st_dev, after.st_ino, after.st_mode, after.st_uid,
            after.st_gid, after.st_nlink)
    ):
        fail(where + " identity changed during tail read")
    if after.st_size != before.st_size or len(raw) != before.st_size - start:
        raise FileNotFoundError(where + " grew during tail snapshot")
    records = raw.split(b"\n")
    if raw.endswith(b"\n"):
        records = records[:-1]
    else:
        records = records[:-1]
    if start:
        records = records[1:]
    if not records or not records[-1] or len(records[-1]) > maximum_line:
        raise FileNotFoundError(where + " lacks a bounded complete tail record")
    return parse_canonical_document(records[-1] + b"\n", where, maximum_line + 1)


def wait_sampler_ready(
    process: ProcessHandle, prepared: PreparedExecution,
    python_info: os.stat_result, *, owner: Optional[OwnershipSlot] = None,
) -> SamplerAuthority:
    deadline = time.monotonic() + SAMPLER_READY_SECONDS
    last_error = "sampler files were absent"
    while time.monotonic() < deadline:
        if process_exited(process):
            fail("sampler exited before readiness")
        held: Dict[str, int] = {}
        authority: Optional[SamplerAuthority] = None
        ownership_transferred = False
        try:
            open_sampler_files(process, held)
            try:
                validation = tail_canonical_jsonl(
                    held["validation"],
                    "sampler validation admission",
                )
            except LaunchError as exc:
                raise FileNotFoundError(exception_text(exc))
            if (
                validation.get("schema") != SAMPLER_VALIDATION_SCHEMA
                or validation.get("decision") != "continue"
                or type(validation.get("monotonic_s")) not in (int, float)
                or not 0.0 <= time.monotonic() - validation["monotonic_s"] <= 6.0
            ):
                raise FileNotFoundError(
                    "sampler has not published a continuing sample")
            security = process_security_receipt(
                process.pid, process.argv, (I2C_GID,), ((SAMPLER_CPU,),),
                python_info,
            )
            if security["process_start_ticks"] != process.start_ticks:
                fail("sampler start identity changed")
            bindings = {
                "csv": live_file_identity(
                    held["csv"], "sampler admission csv", expected_mode=0o600,
                ),
                "validation": live_file_identity(
                    held["validation"], "sampler admission validation",
                    expected_mode=0o444,
                ),
                "pid": file_binding(
                    held["pid"], 4096, "sampler admission pid",
                ),
                "receipt": file_binding(
                    held["receipt"], 1024 * 1024,
                    "sampler admission receipt", allow_empty=True,
                ),
            }
            csv_prefix, csv_prefix_binding = stable_live_file_snapshot(
                held["csv"], 16 * 1024 * 1024,
                "sampler admission CSV", expected_mode=0o600,
            )
            validation_prefix, validation_prefix_binding = stable_live_file_snapshot(
                held["validation"], 16 * 1024 * 1024,
                "sampler admission validation", expected_mode=0o444,
            )
            if not csv_prefix.endswith(b"\n") or not validation_prefix.endswith(b"\n"):
                raise FileNotFoundError("sampler admission streams have partial suffixes")
            csv_lines = csv_prefix.splitlines(keepends=True)
            validation_lines = validation_prefix.splitlines(keepends=True)
            if len(validation_lines) < 2 or len(csv_lines) < len(validation_lines):
                raise FileNotFoundError("sampler validation header is absent")
            if len(csv_lines) - len(validation_lines) not in (0, 1):
                fail("sampler admission raw/validation row delta differs")
            csv_prefix = b"".join(csv_lines[:len(validation_lines)])
            expected_header = {
                "expected_output_owner_uid": CAMPAIGN_UID,
                "raw_columns": list(SAMPLER_COLUMNS),
                "sampler_source_expected_sha256": prepared.sampler_sha256,
                "sampling": SAMPLER_SAMPLING,
                "schema": SAMPLER_VALIDATION_STREAM_SCHEMA,
                "thresholds": SAMPLER_THRESHOLDS,
            }
            header = parse_canonical_document(
                validation_lines[0], "sampler admission validation header",
                256 * 1024,
            )
            if header != expected_header:
                fail("sampler admission validation header differs")
            script_digest, script_info = hash_path(
                prepared.sampler_script, MAX_DOCUMENT_BYTES,
                "sampler admission source", expected_nlink=1,
            )
            if script_digest != prepared.sampler_sha256:
                fail("sampler admission source changed")
            parent_info = os.stat(str(SAMPLER_DIR), follow_symlinks=False)
            authority = SamplerAuthority(process, held, {
                "artifacts": bindings,
                "cmdline_argv": list(process.argv),
                "cmdline_sha256": sha256_bytes(exact_cmdline_bytes(process.argv)),
                "cpu": SAMPLER_CPU,
                "csv_bytes": len(csv_prefix),
                "csv_device": csv_prefix_binding["device"],
                "csv_inode": csv_prefix_binding["inode"],
                "csv_path": str(SAMPLER_CSV),
                "csv_sha256": sha256_bytes(csv_prefix),
                "csv_stat": csv_prefix_binding["stat"],
                "evidence_parent": {
                    "path": str(SAMPLER_DIR), "stat": stat_receipt(parent_info),
                },
                "environ_sha256": sha256_bytes(
                    exact_environment_bytes(SEALED_ENVIRONMENT)),
                "environment": dict(SEALED_ENVIRONMENT),
                "environment_sha256": sha256_bytes(
                    canonical_bytes(SEALED_ENVIRONMENT)),
                "executable_sha256": prepared.python_sha256,
                "executable_device": security["executable_stat"]["device"],
                "executable_inode": security["executable_stat"]["inode"],
                "executable_path": str(prepared.python_image),
                "executable_stat": dict(security["executable_stat"]),
                "pid": process.pid,
                "pid_file": {
                    "bytes": bindings["pid"]["bytes"],
                    "path": str(SAMPLER_PID_FILE),
                    "sha256": bindings["pid"]["sha256"],
                    "stat": stat_receipt(os.fstat(held["pid"])),
                },
                "process_affinity": [SAMPLER_CPU],
                "process_gid": CAMPAIGN_GID,
                "process_security": security,
                "process_start_ticks": process.start_ticks,
                "process_uid": CAMPAIGN_UID,
                "receipt_file": {
                    "bytes": bindings["receipt"]["bytes"],
                    "path": str(SAMPLER_RECEIPT),
                    "sha256": bindings["receipt"]["sha256"],
                    "stat": stat_receipt(os.fstat(held["receipt"])),
                },
                "script_path": str(prepared.sampler_script),
                "script_device": script_info.st_dev,
                "script_inode": script_info.st_ino,
                "script_sha256": script_digest,
                "script_stat": stat_receipt(script_info),
                "schema": HEALTH_SAMPLER_ATTESTATION_SCHEMA,
                "terminal_status": "ok",
                "window_end_monotonic_ns": 0,
                "window_start_monotonic_ns": 0,
                "validation": validation,
                "validation_header_ascii": validation_lines[0].decode("ascii"),
                "validation_jsonl": {
                    "bytes": len(validation_prefix),
                    "path": str(SAMPLER_VALIDATION),
                    "sha256": sha256_bytes(validation_prefix),
                    "stat": validation_prefix_binding["stat"],
                },
            })
            if owner is not None:
                try:
                    owner.adopt(authority)
                except BaseException:
                    ownership_transferred = owner.owns(authority)
                    raise
                ownership_transferred = True
            return authority
        except FileNotFoundError as exc:
            if ownership_transferred or (
                owner is not None and authority is not None
                and owner.owns(authority)
            ):
                raise
            last_error = exception_text(exc)
            close_sampler_file_map(held, "retryable sampler admission")
        except LaunchError:
            if ownership_transferred or (
                owner is not None and authority is not None
                and owner.owns(authority)
            ):
                raise
            close_sampler_file_map(held, "failed sampler admission")
            raise
        except BaseException:
            caller_owns = (
                owner is not None and authority is not None
                and owner.owns(authority)
            )
            if not ownership_transferred and not caller_owns:
                close_sampler_file_map(held, "aborted sampler admission")
            raise
        time.sleep(0.02)
    fail("sampler readiness deadline expired: " + last_error)


@dataclass
class SamplerPoller:
    authority: SamplerAuthority
    next_full_poll: float = 0.0
    last_sample_index: int = -1
    full_poll_count: int = 0

    def poll(self) -> None:
        if process_exited(self.authority.process):
            fail("sampler exited while scientific authority was required")
        receipt = os.fstat(self.authority.held["receipt"])
        if receipt.st_size != 0:
            fail("sampler published terminal intent during controller execution")
        now = time.monotonic()
        if now < self.next_full_poll:
            return
        try:
            validation = tail_canonical_jsonl(
                self.authority.held["validation"],
                "sampler live validation",
            )
        except FileNotFoundError:
            # A 1 Hz writer can be observed in the middle of one append.  The
            # pidfd and terminal-intent checks above remain continuous; retry a
            # bounded tail snapshot on the next 100 ms cadence.
            self.next_full_poll = now + 0.1
            return
        index = validation.get("sample_index")
        sample_time = validation.get("monotonic_s")
        if (
            validation.get("schema") != SAMPLER_VALIDATION_SCHEMA
            or validation.get("decision") != "continue"
            or type(index) is not int or index < self.last_sample_index
            or type(sample_time) not in (int, float)
            or not 0.0 <= now - sample_time <= 5.0
        ):
            fail("sampler live validation authority differs")
        self.last_sample_index = index
        self.full_poll_count += 1
        self.next_full_poll = now + 1.0


def controller_run_argv(
    prepared: PreparedExecution, config: LaunchConfig,
    sampler: SamplerAuthority, output_parent: BoundDirectory,
) -> List[str]:
    parent_info = output_parent.verify(0o700)
    values: Dict[str, str] = {
        "current-worker": str(prepared.roles["current"].worker),
        "parent-worker": str(prepared.roles["parent"].worker),
        "current-library": str(prepared.roles["current"].library),
        "parent-library": str(prepared.roles["parent"].library),
        "current-implementation-library":
            str(prepared.roles["current"].implementation_library),
        "parent-implementation-library":
            str(prepared.roles["parent"].implementation_library),
        "current-build-receipt":
            str(prepared.roles["current"].build_receipt),
        "parent-build-receipt":
            str(prepared.roles["parent"].build_receipt),
        "output-dir": str(CONTROLLER_OUTPUT),
        "expected-current-worker-sha256":
            prepared.roles["current"].worker_sha256,
        "expected-parent-worker-sha256":
            prepared.roles["parent"].worker_sha256,
        "expected-current-library-sha256":
            prepared.roles["current"].library_sha256,
        "expected-parent-library-sha256":
            prepared.roles["parent"].library_sha256,
        "expected-current-implementation-library-sha256":
            prepared.roles["current"].implementation_library_sha256,
        "expected-parent-implementation-library-sha256":
            prepared.roles["parent"].implementation_library_sha256,
        "expected-current-build-receipt-sha256":
            prepared.roles["current"].build_receipt_sha256,
        "expected-parent-build-receipt-sha256":
            prepared.roles["parent"].build_receipt_sha256,
        "expected-harness-commit": config.expected_harness_commit,
        "expected-controller-sha256": prepared.controller_sha256,
        "expected-python-image-path": str(prepared.python_image),
        "expected-python-image-sha256": prepared.python_sha256,
        "expected-current-implementation-commit":
            CURRENT_IMPLEMENTATION_COMMIT,
        "expected-parent-implementation-commit": PARENT_IMPLEMENTATION_COMMIT,
        "expected-health-source-manifest-sha256":
            prepared.health_manifest_sha256,
        "expected-health-adapter-sha256": prepared.health_adapter_sha256,
        "cpu": str(WORKER_CPU),
        "controller-cpu": str(CONTROLLER_CPU),
        "sampler-cpu": str(SAMPLER_CPU),
        "sampler-pid": str(sampler.process.pid),
        "sampler-script": str(prepared.sampler_script),
        "sampler-csv": str(SAMPLER_CSV),
        "sampler-pid-file": str(SAMPLER_PID_FILE),
        "sampler-validation-jsonl": str(SAMPLER_VALIDATION),
        "sampler-receipt": str(SAMPLER_RECEIPT),
        "expected-sampler-process-start-ticks":
            str(sampler.process.start_ticks),
        "expected-sampler-script-sha256": prepared.sampler_sha256,
        "expected-sampler-csv-device":
            str(sampler.admission["artifacts"]["csv"]["device"]),
        "expected-sampler-csv-inode":
            str(sampler.admission["artifacts"]["csv"]["inode"]),
        "expected-sampler-pid-file-device":
            str(sampler.admission["artifacts"]["pid"]["device"]),
        "expected-sampler-pid-file-inode":
            str(sampler.admission["artifacts"]["pid"]["inode"]),
        "expected-sampler-validation-device":
            str(sampler.admission["artifacts"]["validation"]["device"]),
        "expected-sampler-validation-inode":
            str(sampler.admission["artifacts"]["validation"]["inode"]),
        "expected-sampler-receipt-device":
            str(sampler.admission["artifacts"]["receipt"]["device"]),
        "expected-sampler-receipt-inode":
            str(sampler.admission["artifacts"]["receipt"]["inode"]),
        "expected-sampler-cmdline-sha256":
            sampler.admission["cmdline_sha256"],
        "expected-sampler-environ-sha256":
            sampler.admission["environ_sha256"],
        "expected-sampler-executable-sha256":
            sampler.admission["executable_sha256"],
        "expected-sampler-uid": str(CAMPAIGN_UID),
        "expected-sampler-gid": str(CAMPAIGN_GID),
        "expected-sampler-i2c-gid": str(I2C_GID),
        "expected-output-parent-device": str(parent_info.st_dev),
        "expected-output-parent-inode": str(parent_info.st_ino),
        "expected-output-parent-uid": str(CAMPAIGN_UID),
        "expected-output-parent-gid": str(CAMPAIGN_GID),
        "expected-output-parent-mode": str(0o700),
    }
    order = (
        "current-worker", "parent-worker", "current-library", "parent-library",
        "current-implementation-library", "parent-implementation-library",
        "current-build-receipt", "parent-build-receipt", "output-dir",
        "expected-current-worker-sha256", "expected-parent-worker-sha256",
        "expected-current-library-sha256", "expected-parent-library-sha256",
        "expected-current-implementation-library-sha256",
        "expected-parent-implementation-library-sha256",
        "expected-current-build-receipt-sha256",
        "expected-parent-build-receipt-sha256", "expected-harness-commit",
        "expected-controller-sha256", "expected-python-image-path",
        "expected-python-image-sha256",
        "expected-current-implementation-commit",
        "expected-parent-implementation-commit",
        "expected-health-source-manifest-sha256",
        "expected-health-adapter-sha256", "cpu", "controller-cpu",
        "sampler-cpu", "sampler-pid", "sampler-script", "sampler-csv",
        "sampler-pid-file", "sampler-validation-jsonl", "sampler-receipt",
        "expected-sampler-process-start-ticks",
        "expected-sampler-script-sha256", "expected-sampler-csv-device",
        "expected-sampler-csv-inode", "expected-sampler-pid-file-device",
        "expected-sampler-pid-file-inode",
        "expected-sampler-validation-device",
        "expected-sampler-validation-inode",
        "expected-sampler-receipt-device",
        "expected-sampler-receipt-inode",
        "expected-sampler-cmdline-sha256",
        "expected-sampler-environ-sha256",
        "expected-sampler-executable-sha256", "expected-sampler-uid",
        "expected-sampler-gid", "expected-sampler-i2c-gid",
        "expected-output-parent-device", "expected-output-parent-inode",
        "expected-output-parent-uid", "expected-output-parent-gid",
        "expected-output-parent-mode",
    )
    if set(values) != set(order):
        fail("controller argv adapter roster differs")
    argv = [
        str(PYTHON), "-I", "-B", str(prepared.controller), "--run-once",
    ]
    for option in order:
        argv.extend(("--" + option, values[option]))
    return argv


def replay_argv(prepared: PreparedExecution, config: LaunchConfig,
                raw_path: Path) -> List[str]:
    exact_absolute(raw_path, "retained replay raw path")
    return [
        str(PYTHON), "-I", "-B", str(prepared.controller),
        "--replay", str(raw_path), "--cpu", str(WORKER_CPU),
        "--expected-harness-commit", config.expected_harness_commit,
        "--expected-current-implementation-commit",
        CURRENT_IMPLEMENTATION_COMMIT,
    ]


def expected_role_map_receipts(role: RoleExecution) -> Dict[Path, Dict[str, Any]]:
    entries = [role.worker_receipt, *role.runtime_map_receipts]
    result: Dict[Path, Dict[str, Any]] = {}
    for entry in entries:
        path = Path(entry["path"])
        canonical = path.resolve(strict=True)
        existing = result.get(canonical)
        if existing is not None and (
            existing["sha256"], existing["bytes"]
        ) != (entry["sha256"], entry["bytes"]):
            fail(role.role + " runtime closure aliases unequal receipts")
        result[canonical] = dict(entry)
    return result


def mapped_file_paths(pid: int) -> Tuple[bytes, Tuple[Path, ...]]:
    raw = read_bounded_path(
        Path("/proc/{}/maps".format(pid)), 4 * 1024 * 1024,
        "released worker maps",
    )
    paths: List[Path] = []
    for line in raw.splitlines():
        fields = line.split(None, 5)
        if len(fields) < 5:
            fail("released worker maps row is malformed")
        try:
            permissions = fields[1].decode("ascii")
            int(fields[2], 16)
            device_parts = fields[3].split(b":")
            if len(device_parts) != 2:
                raise ValueError("device")
            int(device_parts[0], 16)
            int(device_parts[1], 16)
            inode = int(fields[4], 10)
        except (UnicodeDecodeError, ValueError) as exc:
            fail("released worker maps metadata differs: " + exception_text(exc))
        if re.fullmatch(r"[r-][w-][x-][ps]", permissions) is None:
            fail("released worker maps permission field differs")
        if permissions[1:3] == "wx":
            fail("released worker maps contains a writable executable mapping")
        if len(fields) == 5:
            if permissions[2] == "x":
                fail("released worker maps contains anonymous executable memory")
            continue
        pathname = fields[5]
        if pathname.startswith(b"[") and pathname.endswith(b"]"):
            if permissions[2] == "x" and pathname not in (b"[vdso]", b"[vsyscall]"):
                fail("released worker maps contains an unlicensed pseudo executable")
            if permissions[1] == "w" and permissions[2] == "x":
                fail("released worker maps pseudo image is writable executable")
            continue
        if pathname.endswith(b" (deleted)") or not pathname.startswith(b"/"):
            fail("released worker maps contains a deleted/nonabsolute file")
        if inode <= 0:
            fail("released worker file mapping lacks an inode")
        try:
            path = Path(pathname.decode("utf-8")).resolve(strict=True)
        except (UnicodeDecodeError, OSError) as exc:
            fail("released worker map path differs: " + exception_text(exc))
        if path not in paths:
            paths.append(path)
    if not paths:
        fail("released worker maps has no file-backed mappings")
    return raw, tuple(sorted(paths, key=str))


def runtime_map_closure(paths: Sequence[Path], where: str,
                        *, raw_maps_bytes: Optional[int] = None,
                        ) -> Dict[str, Any]:
    if not paths or tuple(paths) != tuple(sorted(set(paths), key=str)):
        fail(where + " mapped-path roster differs")
    entries = []
    for path in paths:
        digest, info = hash_path(
            path, MAX_TRACKED_FILE_BYTES,
            where + " mapped image " + str(path), expected_nlink=None,
        )
        if info.st_uid != 0 or info.st_gid != 0:
            fail(where + " mapped image is not root-owned: " + str(path))
        entries.append({
            "bytes": info.st_size, "path": str(path), "sha256": digest,
            "stat": stat_receipt(info),
        })
    result = {
        "entries": entries,
        "path_manifest_sha256": sha256_bytes(canonical_bytes(entries)),
    }
    if raw_maps_bytes is not None:
        result["raw_maps_bytes"] = raw_maps_bytes
    return result


def self_runtime_map_closure() -> Dict[str, Any]:
    raw, paths = mapped_file_paths(os.getpid())
    return runtime_map_closure(
        paths, "root launcher", raw_maps_bytes=len(raw),
    )


def validate_self_runtime_map_closure(authority: Mapping[str, Any]) -> Dict[str, Any]:
    expected = authority.get("launcher_runtime_map_closure")
    if type(expected) is not dict or type(expected.get("entries")) is not list:
        fail("build authority root-launcher map closure differs")
    actual = self_runtime_map_closure()
    if (
        actual["entries"] != expected["entries"]
        or actual["path_manifest_sha256"]
        != expected.get("path_manifest_sha256")
    ):
        fail("root launcher actual mapped-image closure differs")
    return actual


def validate_controller_runtime_map_closure(
    authority: Mapping[str, Any],
) -> Dict[str, Any]:
    expected = authority.get("controller_runtime_map_closure")
    if (
        type(expected) is not dict
        or set(expected) != {"entries", "path_manifest_sha256"}
        or type(expected.get("entries")) is not list
    ):
        fail("build authority controller map closure differs")
    paths: List[Path] = []
    for entry in expected["entries"]:
        if (
            type(entry) is not dict
            or set(entry) != {"bytes", "path", "sha256", "stat"}
            or type(entry.get("path")) is not str
        ):
            fail("controller map authority entry differs")
        paths.append(Path(entry["path"]))
    actual = runtime_map_closure(paths, "controller authority")
    if actual != expected:
        fail("controller mapped-image inputs changed after build seal")
    return actual


def audit_worker_maps(pid: int, role: RoleExecution) -> Dict[str, Any]:
    if process_start_ticks(pid) <= 0:
        fail(role.role + " worker start identity differs")
    expected = expected_role_map_receipts(role)
    raw, paths = mapped_file_paths(pid)
    actual = set(paths)
    unexpected = sorted(actual - set(expected), key=str)
    if unexpected:
        fail(role.role + " worker mapped unexpected files: "
             + ",".join(str(path) for path in unexpected))
    missing = sorted(set(expected) - actual, key=str)
    if missing:
        raise FileNotFoundError(
            role.role + " worker loader closure is not complete yet")
    verified: List[Dict[str, Any]] = []
    for path in paths:
        entry = expected[path]
        digest, info = hash_path(
            path, max(entry["bytes"], 1),
            role.role + " mapped file " + str(path), expected_nlink=None,
        )
        if digest != entry["sha256"] or info.st_size != entry["bytes"]:
            fail(role.role + " mapped file changed: " + str(path))
        verified.append({
            "bytes": info.st_size, "path": str(path), "sha256": digest,
            "stat": stat_receipt(info),
        })
    worker_security = process_security_receipt(
        pid,
        (str(role.worker), "--serve", "--role", role.role,
         "--cpu", str(WORKER_CPU)),
        (), ((WORKER_CPU,),),
        os.stat(str(role.worker), follow_symlinks=False),
    )
    return {
        "maps_bytes": len(raw), "maps_sha256": sha256_bytes(raw),
        "mapped_files": verified, "process_security": worker_security,
        "role": role.role,
    }


def audit_controller_maps(pid: int,
                          prepared: PreparedExecution) -> Dict[str, Any]:
    expected_entries = prepared.build_authority.get(
        "controller_runtime_map_closure", {}).get("entries")
    if type(expected_entries) is not list:
        fail("controller map authority roster differs")
    expected = {Path(entry["path"]): entry for entry in expected_entries}
    raw, paths = mapped_file_paths(pid)
    if set(paths) != set(expected):
        fail("controller mapped image escaped the sealed runtime closure")
    verified: List[Dict[str, Any]] = []
    for path in paths:
        authority = expected[path]
        digest, info = hash_path(
            path, max(authority["bytes"], 1),
            "controller mapped file " + str(path), expected_nlink=None,
        )
        if (
            digest != authority["sha256"]
            or info.st_size != authority["bytes"]
            or stat_receipt(info) != authority["stat"]
        ):
            fail("controller mapped file changed: " + str(path))
        verified.append({
            "bytes": info.st_size, "path": str(path), "sha256": digest,
            "stat": stat_receipt(info),
        })
    return {
        "maps_bytes": len(raw), "maps_sha256": sha256_bytes(raw),
        "mapped_files": verified,
    }


@dataclass
class WorkerMapMonitor:
    prepared: PreparedExecution
    layout: CgroupLayout
    sampler: SamplerPoller
    guardian: GuardianHandle
    controller: ProcessHandle
    controller_argv: Tuple[str, ...]
    controller_security: Mapping[str, Any]
    root_deadline_ns: int
    observed: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    observed_pids: Dict[str, int] = field(default_factory=dict)
    stubs: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    inner_guardian: Optional[Dict[str, Any]] = None
    git_observations: List[Dict[str, Any]] = field(default_factory=list)
    git_by_pid: Dict[int, Dict[str, Any]] = field(default_factory=dict)
    audit_counts: Dict[str, int] = field(default_factory=dict)
    next_map_poll: float = 0.0
    controller_pinned: bool = False
    controller_maps: Optional[Dict[str, Any]] = None
    controller_map_audit_count: int = 0
    scientific_started: bool = False
    scientific_finished: bool = False

    def __post_init__(self) -> None:
        if self.controller_security.get("topology") != {
            "parent_pid": os.getpid(),
            "process_group": self.controller.pid,
            "session": self.controller.pid,
        }:
            fail("controller process topology differs")

    def _update_phase(self, live_pids: Sequence[int]) -> None:
        live = set(live_pids)
        if self.scientific_finished:
            if (
                self.inner_guardian is not None
                and self.inner_guardian["pid"] in live
                or any(pid in live for pid in self.observed_pids.values())
            ):
                fail("scientific processes reappeared after terminal phase")
            return
        if self.inner_guardian is None:
            return
        guardian_live = self.inner_guardian["pid"] in live
        workers_live = any(pid in live for pid in self.observed_pids.values())
        if not guardian_live:
            if set(self.observed) != {"current", "parent"} or workers_live:
                fail("inner guardian disappeared before both workers terminated")
            self.scientific_finished = True

    def _cheap_identity(self, pid: int, receipt: Mapping[str, Any],
                        expected_affinity: Sequence[int],
                        expected_executable: os.stat_result,
                        where: str) -> None:
        try:
            live = os.stat("/proc/{}/exe".format(pid), follow_symlinks=True)
            if (
                process_start_ticks(pid) != receipt["process_start_ticks"]
                or process_cgroup_relative(pid)
                != str(self.layout.experiment)[len(str(CGROUP_ROOT)):]
                or sorted(os.sched_getaffinity(pid)) != list(expected_affinity)
                or full_stat(live) != full_stat(expected_executable)
            ):
                fail(where + " cheap live identity changed")
        except (FileNotFoundError, ProcessLookupError):
            return

    def _classify_python(self, pid: int) -> str:
        if self.inner_guardian is not None and self.inner_guardian["pid"] == pid:
            self._cheap_identity(
                pid, self.inner_guardian["process_security"],
                (AUTHORITY_CPU,),
                os.stat(str(self.prepared.python_image), follow_symlinks=False),
                "inner guardian",
            )
            return "guardian"
        for role_name, prior in self.stubs.items():
            if prior["pid"] == pid:
                self._cheap_identity(
                    pid, prior["process_security"], (WORKER_CPU,),
                    os.stat(str(self.prepared.python_image), follow_symlinks=False),
                    role_name + " worker stub",
                )
                return role_name + "-stub"
        argv = process_argv(pid, "experiment Python role")
        if len(argv) < 5 or argv[:4] != [str(PYTHON), "-I", "-S", "-c"]:
            fail("experiment contains an unlicensed Python command")
        code_hash = sha256_bytes(argv[4].encode("ascii", "strict"))
        if code_hash == FROZEN_INNER_GUARDIAN_CODE_SHA256:
            deadline_ns = (
                int(argv[8])
                if len(argv) == 10 and argv[8].isdigit() else -1
            )
            if (
                len(argv) != 10 or not all(item.isdigit() for item in argv[5:])
                or int(argv[7]) != self.controller.pid
                or not self.root_deadline_ns <= deadline_ns <= (
                    self.root_deadline_ns
                    + int(CONTROLLER_ADMISSION_SECONDS * 1e9)
                )
                or int(argv[9]) != AUTHORITY_CPU
            ):
                fail("inner guardian command line differs")
            receipt = process_security_receipt(
                pid, argv, (), ((AUTHORITY_CPU,),),
                os.stat(str(self.prepared.python_image), follow_symlinks=False),
            )
            if receipt["topology"] != {
                "parent_pid": self.controller.pid,
                "process_group": pid, "session": pid,
            }:
                fail("inner guardian process topology differs")
            if self.inner_guardian is not None and self.inner_guardian["pid"] != pid:
                fail("inner guardian process identity was replaced")
            self.inner_guardian = {
                "argv": argv, "code_sha256": code_hash,
                "deadline_monotonic_ns": deadline_ns, "pid": pid,
                "process_security": receipt,
            }
            self.scientific_started = True
            return "guardian"
        if code_hash != FROZEN_WORKER_STUB_CODE_SHA256:
            fail("experiment Python code image differs")
        if (
            len(argv) != 13 or not argv[5].isdigit()
            or argv[6] != str(WORKER_CPU) or argv[8] != "--serve"
            or argv[9] != "--role" or argv[10] not in ("current", "parent")
            or argv[11:] != ["--cpu", str(WORKER_CPU)]
        ):
            fail("blocked worker-stub command line differs")
        role_name = argv[10]
        if argv[7] != str(self.prepared.roles[role_name].worker):
            fail(role_name + " blocked worker-stub target differs")
        receipt = process_security_receipt(
            pid, argv, (), ((WORKER_CPU,),),
            os.stat(str(self.prepared.python_image), follow_symlinks=False),
        )
        if receipt["topology"] != {
            "parent_pid": self.controller.pid,
            "process_group": pid, "session": pid,
        }:
            fail(role_name + " blocked worker-stub topology differs")
        prior = self.stubs.get(role_name)
        if prior is not None and prior["pid"] != pid:
            fail(role_name + " blocked worker-stub identity was replaced")
        self.stubs[role_name] = {"pid": pid, "process_security": receipt}
        self.scientific_started = True
        return role_name + "-stub"

    def _classify_git(self, pid: int) -> None:
        prior = self.git_by_pid.get(pid)
        if prior is not None:
            self._cheap_identity(
                pid, prior, tuple(prior["affinity"]),
                os.stat(str(GIT), follow_symlinks=False),
                "controller Git child",
            )
            return
        if self.scientific_started and not self.scientific_finished:
            fail("controller Git child appeared during the scientific phase")
        argv = process_argv(pid, "controller Git child")
        roots = {
            str(self.prepared.harness_root),
            str(Path(self.prepared.build_authority["source_authority"]
                     ["current"]["root"])),
            str(Path(self.prepared.build_authority["source_authority"]
                     ["parent"]["root"])),
        }
        screen_fixed = [
            str(GIT), "-c", "core.fsmonitor=false", "-c",
            "core.filemode=false", "-c", "core.checkStat=default", "-c",
            "core.trustctime=true",
        ]
        if len(argv) >= len(screen_fixed) + 5 and argv[:len(screen_fixed)] == screen_fixed:
            if (
                argv[len(screen_fixed)] != "-c"
                or not argv[len(screen_fixed) + 1].startswith("safe.directory=")
                or argv[len(screen_fixed) + 2] != "-C"
            ):
                fail("controller closure-Git root option vector differs")
            root = argv[len(screen_fixed) + 3]
            if (
                argv[len(screen_fixed) + 1] != "safe.directory=" + root
                or root not in roots
            ):
                fail("controller closure-Git source root differs")
            command = tuple(argv[len(screen_fixed) + 4:])
            permitted = {
                ("rev-parse", "--show-toplevel"),
                ("rev-parse", "--verify", "HEAD"),
                ("rev-parse", "--abbrev-ref", "HEAD"),
                ("status", "--porcelain=v1", "-z", "--untracked-files=all"),
                ("ls-files", "--others", "--exclude-standard", "-z"),
                ("ls-files", "--others", "--ignored", "--exclude-standard", "-z"),
                ("ls-files", "--stage", "-z"),
            }
            if command not in permitted and not (
                len(command) == 3
                and command[:2] == ("rev-parse", "--verify")
                and command[2].endswith("^{tree}")
            ) and not (
                len(command) == 5
                and command[:4] == ("ls-tree", "-r", "-z", "--full-tree")
            ):
                fail("controller closure-Git subcommand differs")
            git_environment = {
                "LANG": "C.UTF-8", "LC_ALL": "C.UTF-8",
                "PATH": "/usr/bin:/bin", "TZ": "UTC",
                "GIT_CONFIG_GLOBAL": "/dev/null",
                "GIT_CONFIG_NOSYSTEM": "1", "GIT_NO_LAZY_FETCH": "1",
                "GIT_NO_REPLACE_OBJECTS": "1", "GIT_OPTIONAL_LOCKS": "0",
                "GIT_TERMINAL_PROMPT": "0",
            }
        else:
            match = re.fullmatch(r"/proc/self/fd/([0-9]+)", argv[0])
            if (
                match is None
                or argv[1:5] != [
                    "-c", "core.fsmonitor=false", "-c", "core.filemode=false",
                ]
                or len(argv) < 8 or argv[5] != "-c"
            ):
                fail("controller health-Git fixed option vector differs")
            fd_number = int(match.group(1))
            try:
                passed_fd = os.stat(
                    "/proc/{}/fd/{}".format(pid, fd_number),
                    follow_symlinks=True,
                )
                git_named = os.stat(str(GIT), follow_symlinks=False)
                root = os.readlink("/proc/{}/cwd".format(pid))
            except OSError as exc:
                fail("controller health-Git held authority differs: "
                     + exception_text(exc))
            if not same_inode(passed_fd, git_named) or root != str(
                self.prepared.harness_root
            ):
                fail("controller health-Git FD/cwd authority differs")
            if argv[6] != "safe.directory=" + root:
                fail("controller health-Git safe-directory authority differs")
            command = tuple(argv[7:])
            source = tuple(HEALTH_SOURCE_PATHS)
            permitted_direct = {
                ("rev-parse", "--show-toplevel"),
                ("rev-parse", "--verify", "HEAD^{commit}"),
                ("status", "--porcelain=v1", "--untracked-files=no"),
                ("ls-files", "-v"),
                ("ls-files", "--error-unmatch", "--", *source),
                ("ls-files", "-v", "--error-unmatch", "--", *source),
                ("ls-tree", "--full-tree", "-r", "HEAD", "--", *source),
                ("hash-object", "--no-filters", "--", *source),
                ("diff", "--quiet", "--no-ext-diff", "--no-textconv",
                 "HEAD", "--", *source),
            }
            if command not in permitted_direct:
                fail("controller health-Git subcommand differs")
            git_environment = {
                "GIT_ASKPASS": "/bin/false", "GIT_NO_LAZY_FETCH": "1",
                "GIT_NO_REPLACE_OBJECTS": "1", "GIT_OPTIONAL_LOCKS": "0",
                "GIT_TERMINAL_PROMPT": "0", "LANG": "C", "LC_ALL": "C",
                "PATH": "/usr/bin:/bin", "SSH_ASKPASS": "/bin/false",
            }
        receipt = process_security_receipt(
            pid, argv, (), (RUN_CPUS, (CONTROLLER_CPU,)),
            os.stat(str(GIT), follow_symlinks=False), git_environment,
        )
        if receipt["topology"] != {
            "parent_pid": self.controller.pid,
            "process_group": pid, "session": pid,
        }:
            fail("controller Git child topology differs")
        self.git_by_pid[pid] = receipt
        self.git_observations.append(receipt)

    def poll(self) -> None:
        self.sampler.poll()
        self.guardian.poll()
        now = time.monotonic()
        if now < self.next_map_poll:
            return
        self.next_map_poll = now + 0.1
        live_pids = CgroupTree.pids(self.layout.experiment)
        self._update_phase(live_pids)
        if self.controller.returncode is None:
            if self.controller.pid not in live_pids:
                fail("controller left the exact experiment cgroup")
            affinity = sorted(os.sched_getaffinity(self.controller.pid))
            if affinity == [CONTROLLER_CPU]:
                self.controller_pinned = True
                if self.controller_maps is None:
                    self.controller_maps = audit_controller_maps(
                        self.controller.pid, self.prepared,
                    )
                    self.controller_map_audit_count += 1
            if (
                process_start_ticks(self.controller.pid) != self.controller.start_ticks
                or process_cgroup_relative(self.controller.pid)
                != str(self.layout.experiment)[len(str(CGROUP_ROOT)):]
                or affinity not in (list(RUN_CPUS), [CONTROLLER_CPU])
                or self.controller_pinned and affinity != [CONTROLLER_CPU]
                or self.inner_guardian is not None and affinity != [CONTROLLER_CPU]
                or full_stat(os.stat(
                    "/proc/{}/exe".format(self.controller.pid),
                    follow_symlinks=True,
                )) != full_stat(os.stat(
                    str(self.prepared.python_image), follow_symlinks=False,
                ))
            ):
                fail("controller live identity changed")
        for pid in live_pids:
            if pid == self.controller.pid:
                continue
            classified = False
            for role_name, role in self.prepared.roles.items():
                try:
                    live = os.stat(
                        "/proc/{}/exe".format(pid), follow_symlinks=True,
                    )
                    expected = os.stat(str(role.worker), follow_symlinks=False)
                except (FileNotFoundError, ProcessLookupError):
                    continue
                if not same_inode(live, expected):
                    continue
                classified = True
                prior = self.observed_pids.get(role_name)
                if prior is not None and prior != pid:
                    fail(role_name + " worker process identity was replaced")
                if prior == pid and role_name in self.observed:
                    if (
                        process_start_ticks(pid)
                        != self.observed[role_name]["process_security"][
                            "process_start_ticks"]
                        or process_cgroup_relative(pid)
                        != str(self.layout.experiment)[len(str(CGROUP_ROOT)):]
                        or sorted(os.sched_getaffinity(pid)) != [WORKER_CPU]
                    ):
                        fail(role_name + " worker live identity changed")
                    continue
                try:
                    receipt = audit_worker_maps(pid, role)
                except FileNotFoundError:
                    continue
                self.observed_pids[role_name] = pid
                self.scientific_started = True
                if receipt["process_security"]["topology"] != {
                    "parent_pid": self.controller.pid,
                    "process_group": pid, "session": pid,
                }:
                    fail(role_name + " released worker topology differs")
                self.observed[role_name] = receipt
                self.audit_counts[role_name] = (
                    self.audit_counts.get(role_name, 0) + 1
                )
                break
            if classified:
                continue
            try:
                live = os.stat("/proc/{}/exe".format(pid), follow_symlinks=True)
                python = os.stat(
                    str(self.prepared.python_image), follow_symlinks=False,
                )
            except (FileNotFoundError, ProcessLookupError):
                continue
            if same_inode(live, python):
                self._classify_python(pid)
                continue
            git_info = os.stat(str(GIT), follow_symlinks=False)
            if same_inode(live, git_info):
                self._classify_git(pid)
                continue
            fail("experiment contains an unlicensed process image")

    def finish(self) -> Dict[str, Any]:
        self._update_phase(())
        if set(self.observed) != {"current", "parent"}:
            fail("released worker maps were not observed for both roles")
        if self.audit_counts != {"current": 1, "parent": 1}:
            fail("released worker map content audit count differs")
        if self.inner_guardian is None:
            fail("inner controller guardian was never observed")
        if (
            not self.scientific_finished or not self.controller_pinned
            or self.controller_maps is None
            or self.controller_map_audit_count != 1
        ):
            fail("controller scientific phase transition was incomplete")
        # Re-hash the frozen map inputs after worker exit without repeatedly
        # touching them during the timed campaign.
        for role in self.prepared.roles.values():
            for path, entry in expected_role_map_receipts(role).items():
                digest, info = hash_path(
                    path, max(entry["bytes"], 1),
                    role.role + " post-run map input", expected_nlink=None,
                )
                if digest != entry["sha256"] or info.st_size != entry["bytes"]:
                    fail(role.role + " post-run map input changed")
        for entry in self.controller_maps["mapped_files"]:
            digest, info = hash_path(
                Path(entry["path"]), max(entry["bytes"], 1),
                "controller post-run mapped file", expected_nlink=None,
            )
            if (
                digest != entry["sha256"]
                or stat_receipt(info) != entry["stat"]
            ):
                fail("controller post-run mapped input changed")
        return {
            "controller_maps": dict(self.controller_maps),
            "inner_guardian": dict(self.inner_guardian),
            "git_observations": list(self.git_observations),
            "observation_limit": PROCESS_OBSERVATION_LIMIT,
            "workers": dict(self.observed),
            "worker_stubs": dict(self.stubs),
        }


def validate_null_placeholder_hash(value: Mapping[str, Any], schema: str,
                                   field_name: str) -> str:
    if type(value) is not dict or value.get("schema") != schema:
        fail("document schema differs: " + schema)
    digest = value.get(field_name)
    if type(digest) is not str or LOWER64.fullmatch(digest) is None:
        fail("document self-hash field differs: " + schema)
    preimage = dict(value)
    preimage[field_name] = None
    if sha256_bytes(canonical_bytes(preimage)) != digest:
        fail("document self-hash differs: " + schema)
    return digest


@dataclass
class ChildBundle:
    directory_fd: int
    files: Dict[str, int]
    data: Dict[str, bytes]
    receipts: Dict[str, Dict[str, Any]]
    complete: Dict[str, Any]
    summary: Optional[Dict[str, Any]]
    provenance: Optional[Dict[str, Any]]
    state: str = "open"
    file_identities: Dict[str, Tuple[int, int]] = field(default_factory=dict)
    directory_identity: Optional[Tuple[int, int]] = None

    @staticmethod
    def _close_bound_number(fd: int, identity: Tuple[int, int],
                            where: str) -> bool:
        """Close one descriptor without trusting a stale numeric FD.

        False means the exact descriptor remains owned after a reported close
        failure.  True means it is closed (including a prior close followed
        by an asynchronous state-handoff exception).  A reused number is a
        hard failure and is never closed.
        """
        try:
            current = os.fstat(fd)
        except OSError as exc:
            if exc.errno == errno.EBADF:
                return True
            raise
        if (current.st_dev, current.st_ino) != identity:
            fail(where + " descriptor number was reused")
        try:
            os.close(fd)
        except BaseException:
            try:
                current = os.fstat(fd)
            except OSError as probe:
                if probe.errno == errno.EBADF:
                    return True
                raise
            if (current.st_dev, current.st_ino) != identity:
                fail(where + " descriptor number was reused during close")
            return False
        return True

    def close(self) -> None:
        if self.state == "closed":
            if (
                self.directory_fd >= 0 or self.files
                or self.directory_identity is not None
                or self.file_identities
            ):
                fail("closed child bundle retained descriptors")
            return
        if self.state not in ("open", "closing"):
            fail("child bundle closure state differs")
        # Publish the resumable state before the first irreversible close.
        # A retry after any descriptor-close handoff skips content use and
        # drains only the still-owned descriptor members.
        self.state = "closing"
        failures: List[str] = []
        for name in list(dict.fromkeys(
            list(self.files) + list(self.file_identities)
        )):
            if name not in self.files:
                self.file_identities.pop(name, None)
                continue
            if name not in self.file_identities:
                failures.append(name + ": descriptor identity is absent")
                continue
            fd = self.files[name]
            try:
                closed = self._close_bound_number(
                    fd, self.file_identities[name], "child bundle " + name,
                )
            except BaseException as exc:
                failures.append(name + ": " + exception_text(exc))
                continue
            if not closed:
                failures.append(name + ": exact descriptor remains open")
                continue
            self.files.pop(name, None)
            self.file_identities.pop(name, None)
        if self.directory_fd >= 0:
            if self.directory_identity is None:
                failures.append("directory: descriptor identity is absent")
            else:
                try:
                    closed = self._close_bound_number(
                        self.directory_fd, self.directory_identity,
                        "child bundle directory",
                    )
                except BaseException as exc:
                    failures.append("directory: " + exception_text(exc))
                else:
                    if closed:
                        self.directory_fd = -1
                        self.directory_identity = None
                    else:
                        failures.append(
                            "directory: exact descriptor remains open"
                        )
        elif self.directory_identity is not None:
            # A prior close completed before its final identity handoff.
            self.directory_identity = None
        if failures:
            fail("child bundle descriptor closure differs: "
                 + " | ".join(failures))
        self.state = "closed"


def open_child_bundle(parent: BoundDirectory, *,
                      owner: Optional[OwnershipSlot] = None) -> ChildBundle:
    parent.verify(0o700)
    if sorted(os.listdir(parent.fd)) != [CONTROLLER_OUTPUT.name]:
        fail("controller output parent roster differs")
    baseline = open_fd_roster()
    owned: List[int] = []
    files: Dict[str, int] = {}
    bundle: Optional[ChildBundle] = None
    try:
        open_registered(
            owned, CONTROLLER_OUTPUT.name,
            os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC
            | getattr(os, "O_NOFOLLOW", 0), dir_fd=parent.fd,
        )
        directory_fd = owned[0]
        directory = os.fstat(directory_fd)
        named_directory = os.stat(
            CONTROLLER_OUTPUT.name, dir_fd=parent.fd, follow_symlinks=False,
        )
        if (
            not stat.S_ISDIR(directory.st_mode)
            or full_stat(directory) != full_stat(named_directory)
            or directory.st_uid != CAMPAIGN_UID
            or directory.st_gid != CAMPAIGN_GID
            or stat.S_IMODE(directory.st_mode) != 0o500
        ):
            fail("controller child directory seal differs")
        if sorted(os.listdir(directory_fd)) != sorted(CONTROLLER_OUTPUT_NAMES):
            fail("controller child output roster differs")
        maxima = {
            "raw.jsonl": MAX_RAW_BYTES,
            "summary.json": MAX_DOCUMENT_BYTES,
            "provenance.json": MAX_DOCUMENT_BYTES,
            "current.stderr": MAX_STDERR_BYTES,
            "parent.stderr": MAX_STDERR_BYTES,
            "COMPLETE": MAX_DOCUMENT_BYTES,
        }
        data: Dict[str, bytes] = {}
        receipts: Dict[str, Dict[str, Any]] = {}
        for name in CONTROLLER_OUTPUT_NAMES:
            open_registered(
                owned, name,
                os.O_RDONLY | os.O_CLOEXEC | os.O_NONBLOCK
                | getattr(os, "O_NOFOLLOW", 0), dir_fd=directory_fd,
            )
            fd = owned[-1]
            files[name] = fd
            before = os.fstat(fd)
            named = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
            if (
                not stat.S_ISREG(before.st_mode)
                or full_stat(before) != full_stat(named)
                or before.st_uid != CAMPAIGN_UID
                or before.st_gid != CAMPAIGN_GID or before.st_nlink != 1
                or stat.S_IMODE(before.st_mode) != 0o400
                or before.st_size > maxima[name]
                or (name not in ("current.stderr", "parent.stderr")
                    and before.st_size < 1)
            ):
                fail("controller child file policy differs: " + name)
            raw = read_held_file(
                fd, maxima[name], "controller child " + name,
                allow_empty=name in ("current.stderr", "parent.stderr"),
            )
            after = os.fstat(fd)
            named_after = os.stat(
                name, dir_fd=directory_fd, follow_symlinks=False,
            )
            if full_stat(after) != full_stat(named_after):
                fail("controller child name changed across read: " + name)
            data[name] = raw
            receipts[name] = {
                "bytes": len(raw), "name": name, "sha256": sha256_bytes(raw),
                "stat": stat_receipt(after),
            }
        complete = parse_canonical_document(
            data["COMPLETE"], "controller child COMPLETE", MAX_DOCUMENT_BYTES,
        )
        validate_null_placeholder_hash(
            complete, CHILD_COMPLETE_SCHEMA, "preimage_sha256",
        )
        if (
            set(complete) != {
                "campaign", "files", "outcome", "preimage_sha256", "schema",
            }
            or complete.get("campaign") != CAMPAIGN
            or complete.get("outcome") not in ("pass", "reject", "invalid")
            or type(complete.get("files")) is not list
            or len(complete["files"]) != 5
        ):
            fail("controller child COMPLETE contract differs")
        expected_file_names = list(CONTROLLER_OUTPUT_NAMES[:-1])
        for index, entry in enumerate(complete["files"]):
            name = expected_file_names[index]
            if (
                type(entry) is not dict
                or set(entry) != {"bytes", "name", "sha256"}
                or entry.get("name") != name
                or entry.get("bytes") != receipts[name]["bytes"]
                or entry.get("sha256") != receipts[name]["sha256"]
            ):
                fail("controller child COMPLETE file receipt differs: " + name)
        summary_value = parse_canonical_document(
            data["summary.json"], "controller child summary", MAX_DOCUMENT_BYTES,
        )
        provenance_value = parse_canonical_document(
            data["provenance.json"], "controller child provenance",
            MAX_DOCUMENT_BYTES,
        )
        summary: Optional[Dict[str, Any]] = None
        provenance: Optional[Dict[str, Any]] = None
        if complete["outcome"] in ("pass", "reject"):
            validate_null_placeholder_hash(
                summary_value, CHILD_SUMMARY_SCHEMA,
                "summary_preimage_sha256",
            )
            validate_null_placeholder_hash(
                provenance_value, CHILD_PROVENANCE_SCHEMA, "preimage_sha256",
            )
            if (
                summary_value.get("campaign") != CAMPAIGN
                or summary_value.get("outcome") != complete["outcome"]
                or summary_value.get("valid_evidence") is not True
                or summary_value.get("raw_complete") is not True
                or summary_value.get("raw_bytes") != len(data["raw.jsonl"])
                or summary_value.get("raw_sha256")
                != sha256_bytes(data["raw.jsonl"])
                or summary_value.get("provenance_sha256")
                != sha256_bytes(data["provenance.json"])
                or provenance_value.get("campaign") != CAMPAIGN
            ):
                fail("controller child summary/provenance cross-binding differs")
            summary = summary_value
            provenance = provenance_value
        elif (
            summary_value.get("schema")
            != "wirehair.wh2.v2-facade-timing-screen.emergency-invalid.v1"
            and summary_value.get("schema") != CHILD_SUMMARY_SCHEMA
        ):
            fail("invalid child summary has an unknown schema")
        directory_after = os.fstat(directory_fd)
        named_after = os.stat(
            CONTROLLER_OUTPUT.name, dir_fd=parent.fd, follow_symlinks=False,
        )
        if (
            full_stat(directory_after) != full_stat(named_after)
            or full_stat(directory_after) != full_stat(directory)
            or sorted(os.listdir(directory_fd)) != sorted(CONTROLLER_OUTPUT_NAMES)
        ):
            fail("controller child directory changed across held read")
        bundle = ChildBundle(
            directory_fd, files, data, receipts, complete, summary, provenance,
            file_identities={
                name: (
                    receipts[name]["stat"]["device"],
                    receipts[name]["stat"]["inode"],
                )
                for name in files
            },
            directory_identity=(directory.st_dev, directory.st_ino),
        )
        if owner is not None:
            owner.adopt(bundle)
        return bundle
    except BaseException as original:
        # Once the caller's persistent slot owns the complete bundle it is
        # responsible for resumable closure, including an exception after the
        # adoption callback but before this function returns.  Before that
        # handoff, close the exact single-threaded descriptor delta so an
        # open-return or dict-assignment fault cannot leak a child artifact.
        if bundle is not None and owner is not None and owner.owns(bundle):
            raise
        cleanup_failure: Optional[BaseException] = None
        try:
            close_fd_delta(baseline, ())
        except BaseException as exc:
            cleanup_failure = exc
        if cleanup_failure is not None:
            raise LaunchError(
                exception_text(original) + " | child-bundle cleanup: "
                + exception_text(cleanup_failure)
            ) from original
        raise


def validate_child_bundle_binding(parent: BoundDirectory,
                                  bundle: ChildBundle) -> None:
    """Reprove the sealed nested namespace before sealing its parent."""
    parent.verify(0o700)
    if bundle.directory_fd < 0 or set(bundle.files) != set(CONTROLLER_OUTPUT_NAMES):
        fail("controller child held descriptor roster differs")
    directory = os.fstat(bundle.directory_fd)
    named_directory = os.stat(
        CONTROLLER_OUTPUT.name, dir_fd=parent.fd, follow_symlinks=False,
    )
    if (
        not stat.S_ISDIR(directory.st_mode)
        or full_stat(directory) != full_stat(named_directory)
        or directory.st_uid != CAMPAIGN_UID
        or directory.st_gid != CAMPAIGN_GID
        or stat.S_IMODE(directory.st_mode) != 0o500
        or tuple(sorted(os.listdir(bundle.directory_fd)))
        != tuple(sorted(CONTROLLER_OUTPUT_NAMES))
    ):
        fail("controller child nested seal changed")
    for name in CONTROLLER_OUTPUT_NAMES:
        fd = bundle.files[name]
        held = os.fstat(fd)
        named = os.stat(
            name, dir_fd=bundle.directory_fd, follow_symlinks=False,
        )
        receipt = bundle.receipts[name]
        raw = bundle.data[name]
        if (
            not stat.S_ISREG(held.st_mode)
            or full_stat(held) != full_stat(named)
            or held.st_uid != CAMPAIGN_UID or held.st_gid != CAMPAIGN_GID
            or held.st_nlink != 1 or stat.S_IMODE(held.st_mode) != 0o400
            or receipt.get("bytes") != len(raw)
            or receipt.get("name") != name
            or receipt.get("sha256") != sha256_bytes(raw)
            or receipt.get("stat") != stat_receipt(held)
        ):
            fail("controller child nested file binding changed: " + name)


def copy_child_bundle(bundle: ChildBundle, journal: AttemptJournal,
                      state: Optional[Dict[str, Dict[str, Any]]] = None,
                      ) -> Dict[str, Dict[str, Any]]:
    """Idempotently copy every held child artifact into the root journal."""
    result = state if state is not None else {}
    for child_name, root_name in (
        ("raw.jsonl", "screen.raw.jsonl"),
        ("summary.json", "screen.summary.json"),
        ("provenance.json", "screen.provenance.json"),
        ("current.stderr", "screen.current.stderr"),
        ("parent.stderr", "screen.parent.stderr"),
        ("COMPLETE", "screen.child-COMPLETE.json"),
    ):
        receipt = journal.write_bytes_resumable(
            root_name, bundle.data[child_name],
        )
        if child_name in result and result[child_name] != receipt:
            fail("root child-copy receipt changed: " + child_name)
        result[child_name] = receipt
    if set(result) != set(CONTROLLER_OUTPUT_NAMES):
        fail("root child-copy receipt roster differs")
    return result


def _security_projection_matches(child: Any, root: Mapping[str, Any]) -> bool:
    if type(child) is not dict or type(root) is not dict:
        return False
    return child == {
        "affinity": root.get("affinity"),
        "caps": {
            name: "0000000000000000" for name in (
                "CapInh", "CapPrm", "CapEff", "CapBnd", "CapAmb",
            )
        },
        "cmdline_sha256": root.get("cmdline_sha256"),
        "environment_sha256": root.get("environ_sha256"),
        "executable_stat": root.get("executable_stat"),
        "gid": CAMPAIGN_GID,
        "groups": root.get("groups"),
        "no_new_privs": 1,
        "pid": root.get("pid"),
        "process_start_ticks": root.get("process_start_ticks"),
        "uid": CAMPAIGN_UID,
    }


def _artifact_projection(path: Path, digest: str, name: str,
                         expected_nlink: Optional[int] = 1) -> Dict[str, Any]:
    actual, info = hash_path(
        path, MAX_TRACKED_FILE_BYTES, "child provenance artifact " + name,
        expected_nlink=expected_nlink,
    )
    if actual != digest:
        fail("child provenance artifact changed: " + name)
    return {
        "name": name, "path": str(path), "sha256": actual,
        "stat": stat_receipt(info),
    }


def validate_child_authority(
    bundle: ChildBundle, prepared: PreparedExecution, config: LaunchConfig,
    controller: ProcessHandle, controller_security: Mapping[str, Any],
    map_receipt: Mapping[str, Any], root_t0_ns: Optional[int],
    root_deadline_ns: Optional[int], controller_reap_ns: Optional[int],
) -> None:
    if bundle.summary is None or bundle.provenance is None:
        fail("valid child evidence lacks summary/provenance")
    summary = bundle.summary
    provenance = bundle.provenance
    if set(summary) != {
        "campaign", "failures", "health", "outcome", "provenance_sha256",
        "raw_bytes", "raw_complete", "raw_sha256", "schema", "statistics",
        "summary_preimage_sha256", "valid_evidence", "worker_terminals",
    }:
        fail("valid child summary key roster differs")
    expected_provenance_keys = {
        "artifacts", "build_receipts", "campaign",
        "closure_verification_sha256", "controller_end_monotonic_ns",
        "controller_security_boundary", "controller_start_monotonic_ns",
        "error", "expected_current_implementation_commit",
        "expected_controller_sha256", "expected_harness_commit",
        "expected_parent_implementation_commit", "health",
        "health_adapter_sha256", "health_module_loader",
        "health_source_git_receipt", "health_source_manifest", "guardian",
        "internal_deadline_seconds", "launch_argv", "outer_guardian",
        "output_reservation", "preimage_sha256", "schema",
        "sealed_environment", "target_cpu", "trusted_security_assumption",
        "worker_processes",
    }
    if set(provenance) != expected_provenance_keys:
        fail("valid child provenance key roster differs")
    if (
        type(root_t0_ns) is not int or type(root_deadline_ns) is not int
        or type(controller_reap_ns) is not int
        or root_deadline_ns
        != root_t0_ns + EXTERNAL_DEADLINE_SECONDS * 1_000_000_000
        or
        provenance.get("campaign") != CAMPAIGN
        or provenance.get("expected_harness_commit")
        != config.expected_harness_commit
        or provenance.get("expected_current_implementation_commit")
        != CURRENT_IMPLEMENTATION_COMMIT
        or provenance.get("expected_parent_implementation_commit")
        != PARENT_IMPLEMENTATION_COMMIT
        or provenance.get("expected_controller_sha256")
        != prepared.controller_sha256
        or provenance.get("health_adapter_sha256")
        != prepared.health_adapter_sha256
        or provenance.get("internal_deadline_seconds")
        != INTERNAL_DEADLINE_SECONDS
        or provenance.get("sealed_environment") != SEALED_ENVIRONMENT
        or provenance.get("target_cpu") != WORKER_CPU
        or provenance.get("error") != "none"
        or type(provenance.get("controller_start_monotonic_ns")) is not int
        or type(provenance.get("controller_end_monotonic_ns")) is not int
        or provenance["controller_start_monotonic_ns"]
        > provenance["controller_end_monotonic_ns"]
        or not root_t0_ns <= provenance["controller_start_monotonic_ns"]
        <= provenance["controller_end_monotonic_ns"] <= controller_reap_ns
        < root_deadline_ns
    ):
        fail("valid child provenance campaign/deadline authority differs")
    if (
        provenance.get("health_source_manifest") != prepared.health_manifest
        or provenance.get("health_module_loader")
        != prepared.health_module_loader
        or provenance.get("health_source_git_receipt")
        != prepared.health_source_git_receipt
        or provenance.get("trusted_security_assumption")
        != (
            "root-owned-0555 source directories and no concurrent hostile "
            "UID1000 process; root launcher independently journals outputs"
        )
        or provenance.get("outer_guardian") != {
            "code_sha256": FROZEN_INNER_GUARDIAN_CODE_SHA256,
            "deadline_seconds": EXTERNAL_DEADLINE_SECONDS,
            "python": str(PYTHON),
        }
    ):
        fail("valid child health-source/trust authority differs")
    controller_root = controller_security.get("final_process")
    if (
        type(controller_root) is not dict
        or controller_root.get("pid") != controller.pid
        or controller_root.get("process_start_ticks") != controller.start_ticks
        or not _security_projection_matches(
            provenance.get("controller_security_boundary"), controller_root,
        )
    ):
        fail("valid child controller security boundary differs")

    expected_receipts: Dict[str, Dict[str, Any]] = {}
    expected_artifacts: List[Dict[str, Any]] = []
    for role in ("current", "parent"):
        item = prepared.roles[role]
        receipt_raw = item.build_receipt.read_bytes()
        receipt = parse_canonical_document(
            receipt_raw, role + " child-provenance build receipt",
            MAX_DOCUMENT_BYTES,
        )
        expected_receipts[role] = receipt
        expected_artifacts.extend((
            _artifact_projection(item.worker, item.worker_sha256,
                                 role + "_worker"),
            _artifact_projection(item.library, item.library_sha256,
                                 role + "_library"),
            _artifact_projection(
                item.implementation_library,
                item.implementation_library_sha256,
                role + "_implementation_library",
            ),
            _artifact_projection(item.build_receipt,
                                 item.build_receipt_sha256,
                                 role + "_build_receipt"),
        ))
    expected_artifacts.extend((
        _artifact_projection(
            prepared.controller, prepared.controller_sha256, "controller",
        ),
        _artifact_projection(
            prepared.python_image, prepared.python_sha256, "python_image",
            expected_nlink=FROZEN_TOOL_NLINK[prepared.python_image],
        ),
        _artifact_projection(
            Path(expected_receipts["current"]["git_path"]),
            expected_receipts["current"]["git_sha256"], "Git executable",
        ),
    ))
    if (
        provenance.get("build_receipts") != expected_receipts
        or sorted(provenance.get("artifacts", []), key=lambda value: value["name"])
        != sorted(expected_artifacts, key=lambda value: value["name"])
    ):
        fail("valid child build/artifact provenance differs")
    closure = provenance.get("closure_verification_sha256")
    if (
        type(closure) is not dict or set(closure) != {"current", "parent"}
        or closure != prepared.closure_sha256
    ):
        fail("valid child closure digest authority differs")

    expected_launch = {
        role: [str(prepared.roles[role].worker), "--serve", "--role", role,
               "--cpu", str(WORKER_CPU)]
        for role in ("current", "parent")
    }
    if provenance.get("launch_argv") != expected_launch:
        fail("valid child worker argv provenance differs")
    workers = provenance.get("worker_processes")
    observed_workers = map_receipt.get("workers")
    observed_stubs = map_receipt.get("worker_stubs")
    if (
        type(workers) is not dict or set(workers) != {"current", "parent"}
        or type(observed_workers) is not dict
        or type(observed_stubs) is not dict
        or set(observed_workers) != {"current", "parent"}
        or not set(observed_stubs).issubset({"current", "parent"})
        or map_receipt.get("observation_limit") != PROCESS_OBSERVATION_LIMIT
    ):
        fail("valid child worker process roster differs")
    expected_stub_environment_sha256 = sha256_bytes(
        exact_environment_bytes(SEALED_ENVIRONMENT),
    )
    expected_stub_executable = stat_receipt(os.stat(
        str(prepared.python_image), follow_symlinks=False,
    ))
    for role in ("current", "parent"):
        value = workers[role]
        released_root = observed_workers[role]["process_security"]
        blocked_child = (
            value.get("blocked_stub") if type(value) is dict else None
        )
        released_child = (
            value.get("released_worker") if type(value) is dict else None
        )
        terminal_child = (
            value.get("terminal_boundary") if type(value) is dict else None
        )
        if (
            type(value) is not dict
            or set(value) != {
                "artifact", "blocked_stub", "pgid", "released_worker",
                "terminal_boundary",
            }
            or value.get("pgid") != released_root["pid"]
            or type(blocked_child) is not dict
            or type(released_child) is not dict
            or type(terminal_child) is not dict
            or set(blocked_child) != {
                "affinity", "caps", "cmdline_sha256", "environment_sha256",
                "executable_stat", "gid", "groups", "no_new_privs", "pid",
                "process_start_ticks", "uid",
            }
            or blocked_child.get("pid") != value.get("pgid")
            or released_child.get("pid") != value.get("pgid")
            or terminal_child.get("pid") != value.get("pgid")
            or blocked_child.get("affinity") != [WORKER_CPU]
            or blocked_child.get("caps") != {
                name: "0000000000000000" for name in (
                    "CapInh", "CapPrm", "CapEff", "CapBnd", "CapAmb",
                )
            }
            or blocked_child.get("environment_sha256")
            != expected_stub_environment_sha256
            or blocked_child.get("executable_stat") != expected_stub_executable
            or blocked_child.get("gid") != CAMPAIGN_GID
            or blocked_child.get("groups") != []
            or blocked_child.get("no_new_privs") != 1
            or blocked_child.get("uid") != CAMPAIGN_UID
            or type(blocked_child.get("cmdline_sha256")) is not str
            or LOWER64.fullmatch(blocked_child["cmdline_sha256"]) is None
            or not _security_projection_matches(
                value.get("released_worker"), released_root)
            or not _security_projection_matches(
                value.get("terminal_boundary"), released_root)
            or value["blocked_stub"].get("process_start_ticks")
            != value["released_worker"].get("process_start_ticks")
            or value["terminal_boundary"].get("process_start_ticks")
            != value["released_worker"].get("process_start_ticks")
            or value.get("terminal_boundary") != value.get("released_worker")
            or value.get("artifact") != next(
                item for item in expected_artifacts
                if item["name"] == role + "_worker")
        ):
            fail(role + " valid child worker process authority differs")
        if role in observed_stubs and not _security_projection_matches(
            blocked_child, observed_stubs[role].get("process_security", {}),
        ):
            fail(role + " observed worker-stub diagnostic differs")

    inner = map_receipt.get("inner_guardian")
    guardian = provenance.get("guardian")
    if (
        type(inner) is not dict or type(guardian) is not dict
        or set(guardian) != {
            "affinity", "armed_monotonic_ns", "code_sha256",
            "deadline_monotonic_ns", "deadline_seconds", "pid",
            "process_start_ticks", "registered_worker_pgids", "security",
            "started_monotonic_ns", "worker_stub_sha256",
        }
        or guardian.get("affinity") != [AUTHORITY_CPU]
        or guardian.get("code_sha256")
        != FROZEN_INNER_GUARDIAN_CODE_SHA256
        or guardian.get("worker_stub_sha256")
        != FROZEN_WORKER_STUB_CODE_SHA256
        or guardian.get("deadline_seconds") != EXTERNAL_DEADLINE_SECONDS
        or guardian.get("pid") != inner.get("pid")
        or guardian.get("process_start_ticks")
        != inner.get("process_security", {}).get("process_start_ticks")
        or guardian.get("armed_monotonic_ns") is None
        or guardian.get("deadline_monotonic_ns")
        != provenance["controller_start_monotonic_ns"]
        + EXTERNAL_DEADLINE_SECONDS * 1_000_000_000
        or guardian.get("deadline_monotonic_ns")
        != inner.get("deadline_monotonic_ns")
        or type(guardian.get("started_monotonic_ns")) is not int
        or not guardian["started_monotonic_ns"]
        <= guardian["armed_monotonic_ns"]
        or not guardian["armed_monotonic_ns"]
        < guardian["deadline_monotonic_ns"]
        or not root_deadline_ns <= guardian["deadline_monotonic_ns"]
        <= root_deadline_ns + int(CONTROLLER_ADMISSION_SECONDS * 1e9)
        or guardian.get("registered_worker_pgids")
        != [workers[role]["pgid"] for role in ("current", "parent")]
        or not _security_projection_matches(
            guardian.get("security"), inner.get("process_security", {}),
        )
    ):
        fail("valid child inner guardian authority differs")

    health = provenance.get("health")
    if (
        type(health) is not dict or summary.get("health") != health
        or health.get("evidence_status") != "complete"
        or health.get("terminal_status") != "ok"
        or health.get("collection_failures") != []
        or health.get("violations") != []
        or health.get("target_cpu") != WORKER_CPU
        or health.get("controller_cpu") != CONTROLLER_CPU
        or health.get("sampler_cpu") != SAMPLER_CPU
        or type(health.get("receipt_sha256")) is not str
    ):
        fail("valid child health authority differs")
    health_preimage = dict(health)
    health_digest = health_preimage.get("receipt_sha256")
    health_preimage["receipt_sha256"] = None
    if sha256_bytes(canonical_bytes(health_preimage)) != health_digest:
        fail("valid child health self-hash differs")
    terminals = summary.get("worker_terminals")
    if type(terminals) is not dict or set(terminals) != {"current", "parent"}:
        fail("valid child worker terminal roster differs")
    for role in ("current", "parent"):
        if (
            type(terminals[role]) is not dict
            or set(terminals[role]) != {
                "invocation_count", "role", "schema", "status", "target_cpu",
            }
            or terminals[role].get("role") != role
            or terminals[role].get("status") != "complete"
            or terminals[role].get("target_cpu") != WORKER_CPU
            or terminals[role].get("schema")
            != "wirehair.wh2.v2-facade-timing-worker.terminal.v1"
        ):
            fail(role + " valid child worker terminal differs")
    reservation = provenance.get("output_reservation")
    parent_info = os.stat(str(CONTROLLER_PARENT), follow_symlinks=False)
    child_info = os.fstat(bundle.directory_fd)
    if reservation != {
        "directory": {
            "device": child_info.st_dev, "gid": child_info.st_gid,
            "inode": child_info.st_ino, "mode": 0o700,
            "uid": child_info.st_uid,
        },
        "parent": {
            "device": parent_info.st_dev, "gid": parent_info.st_gid,
            "inode": parent_info.st_ino, "mode": 0o700,
            "uid": parent_info.st_uid,
        },
        "path": str(CONTROLLER_OUTPUT),
    }:
        fail("valid child output-reservation provenance differs")


def validate_replay_result(
    raw_stdout: bytes, raw_stderr: bytes, returncode: int,
    bundle: ChildBundle,
) -> Dict[str, Any]:
    if raw_stderr or returncode not in (0, 1, 2):
        fail("retained replay process failed or emitted stderr")
    replay = parse_canonical_document(
        raw_stdout, "retained replay result", MAX_DOCUMENT_BYTES,
    )
    expected_outcome = {0: "pass", 1: "invalid", 2: "reject"}[returncode]
    raw_terminal = replay.get("raw_terminal")
    worker_terminals = (
        bundle.summary.get("worker_terminals")
        if bundle.summary is not None else None
    )
    if (
        set(replay) != {
            "campaign", "failed_gates", "outcome", "raw_bytes",
            "raw_sha256", "raw_terminal", "schema", "statistics",
        }
        or replay.get("schema") != REPLAY_SCHEMA
        or replay.get("campaign") != CAMPAIGN
        or replay.get("outcome") != expected_outcome
        or replay.get("raw_bytes") != len(bundle.data["raw.jsonl"])
        or replay.get("raw_sha256") != sha256_bytes(bundle.data["raw.jsonl"])
        or bundle.complete.get("outcome") != expected_outcome
        or bundle.summary is None
        or bundle.summary.get("outcome") != expected_outcome
        or bundle.summary.get("statistics") != replay.get("statistics")
        or bundle.summary.get("failures") != replay.get("failed_gates")
        or type(raw_terminal) is not dict
        or type(worker_terminals) is not dict
        or set(worker_terminals) != {"current", "parent"}
        or any(
            type(raw_terminal.get(role + "_invocations")) is not int
            or raw_terminal[role + "_invocations"] <= 0
            or worker_terminals[role].get("invocation_count")
            != raw_terminal[role + "_invocations"]
            for role in ("current", "parent")
        )
    ):
        fail("retained replay/controller evidence cross-check differs")
    return replay


def _sampler_binding_matches(record: Any, observed: Mapping[str, Any]) -> bool:
    return (
        type(record) is dict
        and record.get("device") == observed["device"]
        and record.get("inode") == observed["inode"]
        and record.get("size") == observed["bytes"]
        and record.get("sha256") == observed["sha256"]
    )


def sampler_binding_from_observed(observed: Mapping[str, Any]) -> Dict[str, Any]:
    return {
        "device": observed["device"], "gid": observed["gid"],
        "inode": observed["inode"], "mode": observed["mode"],
        "nlink": observed["nlink"], "sha256": observed["sha256"],
        "size": observed["bytes"], "uid": observed["uid"],
    }


def validate_sampler_destination(value: Any, expected_path: Path,
                                 expected_owner_uid: Optional[int],
                                 expected_parent: os.stat_result,
                                 where: str) -> None:
    if type(value) is not dict or set(value) != {
        "basename", "expected_owner_uid", "parent", "path",
    }:
        fail(where + " destination keys differ")
    parent = value["parent"]
    expected_parent_binding = {
        "device": expected_parent.st_dev, "gid": expected_parent.st_gid,
        "inode": expected_parent.st_ino,
        "mode": stat.S_IMODE(expected_parent.st_mode),
        "nlink": expected_parent.st_nlink, "uid": expected_parent.st_uid,
    }
    if (
        value["path"] != str(expected_path)
        or value["basename"] != expected_path.name
        or value["expected_owner_uid"] != expected_owner_uid
        or type(parent) is not dict or set(parent) != {"binding", "path"}
        or parent["path"] != str(expected_path.parent)
        or parent["binding"] != expected_parent_binding
    ):
        fail(where + " destination authority differs")


def sampler_exact_uint(value: Any, where: str) -> int:
    if type(value) is not int or value < 0 or value > (1 << 64) - 1:
        fail(where + " is not a uint64")
    return value


def sampler_finite(value: Any, where: str) -> float:
    if type(value) not in (int, float) or not math.isfinite(value):
        fail(where + " is not finite")
    return float(value)


def sampler_finite_text(value: str, where: str) -> float:
    if type(value) is not str or not value:
        fail(where + " is not a producer scalar")
    try:
        number = float(value)
    except ValueError as exc:
        fail(where + " is not numeric: " + exception_text(exc))
    if not math.isfinite(number) or str(number) != value:
        fail(where + " is not a canonical producer scalar")
    return number


def sampler_optional_finite_text(value: str, where: str) -> Optional[float]:
    if value == "":
        return None
    return sampler_finite_text(value, where)


def sampler_monotonic_text(value: str, where: str) -> float:
    if type(value) is not str or re.fullmatch(r"(?:0|[1-9][0-9]*)\.[0-9]{6}", value) is None:
        fail(where + " is not a fixed-six producer timestamp")
    number = float(value)
    if not math.isfinite(number) or number < 0 or f"{number:.6f}" != value:
        fail(where + " is not a canonical fixed-six producer timestamp")
    return number


def sampler_utc_ns(value: Any, where: str) -> int:
    if (
        type(value) is not str
        or re.fullmatch(
            r"[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:"
            r"[0-9]{2}\.[0-9]{3}Z", value,
        ) is None
    ):
        fail(where + " UTC syntax differs")
    try:
        parsed = datetime.datetime.strptime(
            value, "%Y-%m-%dT%H:%M:%S.%fZ",
        ).replace(tzinfo=datetime.timezone.utc)
    except ValueError as exc:
        fail(where + " UTC calendar value differs: " + exception_text(exc))
    canonical = parsed.isoformat(timespec="milliseconds").replace("+00:00", "Z")
    if canonical != value:
        fail(where + " UTC canonical value differs")
    return int(parsed.timestamp() * 1_000_000_000)


def sampler_counter_text(value: str, where: str) -> int:
    if not value.isdecimal() or (len(value) > 1 and value.startswith("0")):
        fail(where + " is not a canonical counter")
    return sampler_exact_uint(int(value), where)


def sampler_edac_controller_roster(value: Any, counter: str,
                                   where: str) -> Tuple[int, ...]:
    if type(value) is not list or not value or value != sorted(set(value)):
        fail(where + " EDAC path roster differs")
    result: List[int] = []
    pattern = re.compile(
        r"/sys/devices/system/edac/mc/mc(0|[1-9][0-9]*)/"
        + re.escape(counter) + r"_count"
    )
    for path in value:
        if type(path) is not str:
            fail(where + " EDAC path type differs")
        match = pattern.fullmatch(path)
        if match is None:
            fail(where + " EDAC path shape differs: " + path)
        result.append(int(match.group(1)))
    if len(result) != len(set(result)):
        fail(where + " EDAC controller roster differs")
    return tuple(result)


def parse_sampler_csv_line(raw: bytes, where: str) -> List[str]:
    if not raw.endswith(b"\n") or b"\r" in raw or len(raw) > 16384:
        fail(where + " physical line differs")
    try:
        rows = list(csv.reader([raw[:-1].decode("ascii")], strict=True))
    except (UnicodeDecodeError, csv.Error) as exc:
        fail(where + " CSV syntax differs: " + exception_text(exc))
    if (
        len(rows) != 1
        or (",".join(rows[0]) + "\n").encode("ascii") != raw
    ):
        fail(where + " is not canonical unquoted ASCII CSV")
    return rows[0]


def parse_sampler_validation_sample(value: Any, where: str) -> Dict[str, Any]:
    if type(value) is not dict or set(value) != SAMPLER_VALIDATION_SAMPLE_KEYS:
        fail(where + " validation sample keys differ")
    if value["schema"] != SAMPLER_VALIDATION_SCHEMA:
        fail(where + " validation sample schema differs")
    for name in (
        "consecutive_fault_rows", "edac_ce_delta", "edac_ue_delta",
        "fault_count", "read_error_count", "sample_index",
    ):
        sampler_exact_uint(value[name], where + "." + name)
    if value["decision"] not in (
        "continue", "thermal_abort", "telemetry_abort", "edac_abort",
    ):
        fail(where + " validation decision differs")
    if type(value["monotonic_s"]) is not float:
        fail(where + " validation timestamp type differs")
    sample_time = sampler_finite(value["monotonic_s"], where + ".monotonic_s")
    if sample_time < 0:
        fail(where + " validation timestamp is negative")
    hot = value["hot_sensors"]
    if (
        type(hot) is not list or hot != sorted(set(hot))
        or any(type(name) is not str or name not in SAMPLER_DIMM_FIELDS
               for name in hot)
    ):
        fail(where + " hot-sensor roster differs")
    sensors = value["sensors"]
    if type(sensors) is not dict or set(sensors) != set(SAMPLER_DIMM_FIELDS):
        fail(where + " validation sensor roster differs")
    for name in SAMPLER_DIMM_FIELDS:
        sensor = sensors[name]
        sensor_where = where + ".sensors." + name
        if type(sensor) is not dict or set(sensor) != SAMPLER_VALIDATION_SENSOR_KEYS:
            fail(sensor_where + " keys differ")
        sampler_exact_uint(sensor["attempt_errors"], sensor_where + ".attempt_errors")
        sampler_exact_uint(sensor["hot_streak"], sensor_where + ".hot_streak")
        if type(sensor["hot"]) is not bool or type(sensor["valid"]) is not bool:
            fail(sensor_where + " booleans differ")
        if (
            type(sensor["reason"]) is not str or not sensor["reason"]
            or len(sensor["reason"]) > 128
        ):
            fail(sensor_where + " reason differs")
        for field_name in ("jump_c", "rate_c_per_s", "raw_c"):
            scalar = sensor[field_name]
            if scalar is not None:
                sampler_finite(scalar, sensor_where + "." + field_name)
    return value


def validate_sampler_stream_replay(raw_csv: bytes, validation_raw: bytes,
                                   receipt: Mapping[str, Any]) -> Dict[str, Any]:
    if (
        not raw_csv.endswith(b"\n") or b"\r" in raw_csv
        or not validation_raw.endswith(b"\n") or b"\r" in validation_raw
    ):
        fail("sampler streams contain a partial/non-LF suffix")
    csv_lines = raw_csv.splitlines(keepends=True)
    rows = [
        parse_sampler_csv_line(line, "sampler CSV row {}".format(index))
        for index, line in enumerate(csv_lines)
    ]
    if not rows or tuple(rows[0]) != SAMPLER_COLUMNS or len(rows) < 2:
        fail("sampler CSV header/sample roster differs")
    validation_lines = validation_raw.splitlines(keepends=True)
    values = [
        parse_canonical_document(
            line, "sampler validation record", 256 * 1024,
        )
        for line in validation_lines
    ]
    if values[0] != {
        "expected_output_owner_uid": CAMPAIGN_UID,
        "raw_columns": list(SAMPLER_COLUMNS),
        "sampler_source_expected_sha256": receipt["sampler_source"][
            "expected_sha256"],
        "sampling": SAMPLER_SAMPLING,
        "schema": SAMPLER_VALIDATION_STREAM_SCHEMA,
        "thresholds": SAMPLER_THRESHOLDS,
    }:
        fail("sampler validation stream header differs")
    samples = values[1:]
    if len(samples) != len(rows) - 1:
        fail("sampler raw/validation streams are not fully paired")
    states: Dict[str, Dict[str, Any]] = {
        name: {
            "attempt_errors": 0, "hot_candidate_temp": None,
            "hot_candidate_time": None, "hot_streak": 0,
            "invalid_samples": 0, "last_valid_temp": None,
            "last_valid_time": None, "max_hot_streak": 0,
            "max_raw_c": None, "max_valid_c": None, "raw_samples": 0,
            "read_error_samples": 0, "valid_samples": 0,
        }
        for name in SAMPLER_DIMM_FIELDS
    }
    previous_time: Optional[float] = None
    previous_raw_time: Optional[float] = None
    consecutive_fault_rows = 0
    max_consecutive_fault_rows = 0
    ce_baseline: Optional[int] = None
    ue_baseline: Optional[int] = None
    last_ce = 0
    last_ue = 0
    cpu_values: List[float] = []
    started_ns = receipt.get("started_monotonic_ns")
    finished_ns = receipt.get("finished_monotonic_ns")
    started_utc_ns = sampler_utc_ns(
        receipt.get("started_utc"), "sampler receipt start",
    )
    finished_utc_ns = sampler_utc_ns(
        receipt.get("finished_utc"), "sampler receipt finish",
    )
    if (
        type(started_ns) is not int or type(finished_ns) is not int
        or started_ns < 0 or finished_ns < started_ns
        or finished_utc_ns < started_utc_ns
    ):
        fail("sampler receipt interval differs during stream replay")
    for index, (row, value) in enumerate(zip(rows[1:], samples)):
        where = "sampler paired sample {}".format(index)
        if len(row) != len(SAMPLER_COLUMNS):
            fail("sampler CSV row width differs")
        value = parse_sampler_validation_sample(value, where)
        raw_time = sampler_monotonic_text(row[1], where + " raw time")
        sample_time = value["monotonic_s"]
        cpu_tctl = sampler_finite_text(row[4], where + " CPU Tctl")
        # Every other raw scalar is retained as part of the sampler authority,
        # even though it does not drive this campaign's abort state machine.
        cpu_busy = sampler_optional_finite_text(row[2], where + " CPU busy")
        cpu_mhz = sampler_optional_finite_text(row[3], where + " CPU MHz")
        loads = [
            sampler_finite_text(row[column], where + " load " + str(column))
            for column in (14, 15, 16)
        ]
        if (
            f"{sample_time:.6f}" != row[1]
            or raw_time < 0
            or cpu_busy is not None and not 0.0 <= cpu_busy <= 100.0
            or cpu_mhz is not None and cpu_mhz < 0.0
            or any(value < 0.0 for value in loads)
            or previous_time is not None and sample_time <= previous_time
            or previous_raw_time is not None and raw_time <= previous_raw_time
            or value["sample_index"] != index
        ):
            fail("sampler paired timestamp/index differs")
        sample_ns = int(round(sample_time * 1_000_000_000.0))
        sample_utc_ns = sampler_utc_ns(row[0], where + " UTC")
        if (
            not started_ns <= sample_ns <= finished_ns
            or not started_utc_ns <= sample_utc_ns <= finished_utc_ns
        ):
            fail("sampler paired sample lies outside the producer interval")
        previous_time = sample_time
        previous_raw_time = raw_time
        raw_ce = sampler_counter_text(row[17], where + " EDAC CE")
        raw_ue = sampler_counter_text(row[18], where + " EDAC UE")
        current_ce_baseline = raw_ce - value["edac_ce_delta"]
        current_ue_baseline = raw_ue - value["edac_ue_delta"]
        if (
            current_ce_baseline < 0 or current_ue_baseline < 0
            or ce_baseline is not None
            and (current_ce_baseline != ce_baseline
                 or current_ue_baseline != ue_baseline)
            or index > 0 and (raw_ce < last_ce or raw_ue < last_ue)
        ):
            fail("sampler EDAC replay differs")
        if ce_baseline is None:
            ce_baseline = current_ce_baseline
            ue_baseline = current_ue_baseline
        last_ce, last_ue = raw_ce, raw_ue

        expected_sensors: Dict[str, Dict[str, Any]] = {}
        fault_count = 0
        read_error_count = 0
        hot_sensors: List[str] = []
        for offset, name in enumerate(SAMPLER_DIMM_FIELDS, start=5):
            reported = value["sensors"][name]
            state = states[name]
            attempts = reported["attempt_errors"]
            state["attempt_errors"] += attempts
            temperature = (
                None if row[offset] == ""
                else sampler_finite_text(row[offset], where + " " + name)
            )
            if temperature is not None:
                quarter_degrees = temperature * 4.0
                if (
                    not quarter_degrees.is_integer()
                    or not -1024 <= int(quarter_degrees) <= 1023
                ):
                    fail("sampler DIMM value is outside the SPD5118 domain: " + name)
            if (
                temperature is None and attempts != SAMPLER_SAMPLING["dimm_attempts"]
                or temperature is not None
                and attempts >= SAMPLER_SAMPLING["dimm_attempts"]
            ):
                fail("sampler DIMM retry accounting differs: " + name)
            if temperature is None:
                state["read_error_samples"] += 1
                state["hot_candidate_temp"] = None
                state["hot_candidate_time"] = None
                state["hot_streak"] = 0
                fault_count += 1
                read_error_count += 1
                expected_sensors[name] = {
                    "attempt_errors": attempts, "hot": False,
                    "hot_streak": 0, "jump_c": None,
                    "rate_c_per_s": None, "raw_c": None,
                    "reason": "read_error", "valid": False,
                }
                continue
            state["raw_samples"] += 1
            state["max_raw_c"] = (
                temperature if state["max_raw_c"] is None
                else max(state["max_raw_c"], temperature)
            )
            absolute_valid = (
                SAMPLER_THRESHOLDS["min_plausible_dimm_c_exclusive"]
                < temperature
                < SAMPLER_THRESHOLDS["max_plausible_dimm_c_exclusive"]
            )
            is_hot = temperature >= SAMPLER_THRESHOLDS[
                "dimm_safety_c_inclusive"]
            if is_hot:
                candidate = state["hot_candidate_temp"]
                candidate_time = state["hot_candidate_time"]
                coherent = False
                if candidate is not None and candidate_time is not None:
                    elapsed = sample_time - candidate_time
                    jump = abs(temperature - candidate)
                    coherent = (
                        elapsed > 0
                        and jump <= SAMPLER_THRESHOLDS["max_dimm_jump_c"]
                        and jump / elapsed
                        <= SAMPLER_THRESHOLDS["max_dimm_rate_c_per_s"]
                    )
                state["hot_streak"] = state["hot_streak"] + 1 if coherent else 1
                state["hot_candidate_temp"] = temperature
                state["hot_candidate_time"] = sample_time
                state["max_hot_streak"] = max(
                    state["max_hot_streak"], state["hot_streak"],
                )
                if state["hot_streak"] >= SAMPLER_THRESHOLDS[
                    "hot_confirm_samples"]:
                    hot_sensors.append(name)
            else:
                state["hot_candidate_temp"] = None
                state["hot_candidate_time"] = None
                state["hot_streak"] = 0
            jump_value: Optional[float] = None
            rate_value: Optional[float] = None
            reasons: List[str] = []
            if not absolute_valid:
                reasons.append("absolute_range")
            elif state["last_valid_temp"] is not None:
                elapsed = sample_time - state["last_valid_time"]
                jump_value = abs(temperature - state["last_valid_temp"])
                rate_value = jump_value / elapsed
                if jump_value > SAMPLER_THRESHOLDS["max_dimm_jump_c"]:
                    reasons.append("jump")
                if rate_value > SAMPLER_THRESHOLDS["max_dimm_rate_c_per_s"]:
                    reasons.append("rate")
            valid = not reasons
            if valid:
                state["valid_samples"] += 1
                state["last_valid_temp"] = temperature
                state["last_valid_time"] = sample_time
                state["max_valid_c"] = (
                    temperature if state["max_valid_c"] is None
                    else max(state["max_valid_c"], temperature)
                )
                reason = "ok"
            else:
                state["invalid_samples"] += 1
                fault_count += 1
                reason = "+".join(reasons)
            expected_sensors[name] = {
                "attempt_errors": attempts, "hot": is_hot,
                "hot_streak": state["hot_streak"], "jump_c": jump_value,
                "rate_c_per_s": rate_value, "raw_c": temperature,
                "reason": reason, "valid": valid,
            }
        if sampler_counter_text(row[13], where + " DIMM read errors") \
                != read_error_count:
            fail("sampler raw DIMM read-error count differs")
        consecutive_fault_rows = (
            consecutive_fault_rows + 1 if fault_count else 0
        )
        max_consecutive_fault_rows = max(
            max_consecutive_fault_rows, consecutive_fault_rows,
        )
        decision = (
            "thermal_abort" if hot_sensors else
            "edac_abort" if raw_ce != ce_baseline or raw_ue != ue_baseline else
            "telemetry_abort" if consecutive_fault_rows >= SAMPLER_THRESHOLDS[
                "telemetry_fault_abort_samples"] else
            "continue"
        )
        expected_sample = {
            "consecutive_fault_rows": consecutive_fault_rows,
            "decision": decision,
            "edac_ce_delta": raw_ce - ce_baseline,
            "edac_ue_delta": raw_ue - ue_baseline,
            "fault_count": fault_count, "hot_sensors": hot_sensors,
            "monotonic_s": value["monotonic_s"],
            "read_error_count": read_error_count, "sample_index": index,
            "schema": SAMPLER_VALIDATION_SCHEMA,
            "sensors": expected_sensors,
        }
        if canonical_bytes(value) != canonical_bytes(expected_sample):
            fail("sampler validation state-machine replay differs")
        if decision != "continue":
            fail("stopped sampler contains a terminal validation decision")
        cpu_values.append(cpu_tctl)
    assert ce_baseline is not None and ue_baseline is not None
    summary_sensors = {
        name: {
            key: states[name][key]
            for key in (
                "attempt_errors", "invalid_samples", "max_hot_streak",
                "max_raw_c", "max_valid_c", "raw_samples",
                "read_error_samples", "valid_samples",
            )
        }
        for name in SAMPLER_DIMM_FIELDS
    }
    expected_summary = {
        "consecutive_fault_rows": consecutive_fault_rows,
        "decision": "continue",
        "dimm_attempt_errors_total": sum(
            state["attempt_errors"] for state in states.values()),
        "dimm_invalid_samples_total": sum(
            state["invalid_samples"] for state in states.values()),
        "dimm_read_error_samples_total": sum(
            state["read_error_samples"] for state in states.values()),
        "edac_ce_baseline": ce_baseline, "edac_ce_delta": last_ce - ce_baseline,
        "edac_ce_last": last_ce, "edac_ue_baseline": ue_baseline,
        "edac_ue_delta": last_ue - ue_baseline, "edac_ue_last": last_ue,
        "max_consecutive_fault_rows": max_consecutive_fault_rows,
        "sample_count": len(samples), "sensors": summary_sensors,
    }
    summary = receipt.get("summary")
    if (
        type(summary) is not dict
        or canonical_bytes(summary) != canonical_bytes(expected_summary)
        or canonical_bytes(receipt.get("cpu_tctl_max_c"))
        != canonical_bytes(max(cpu_values))
    ):
        fail("sampler terminal summary does not replay from paired streams")
    return {
        "cpu_tctl_max_c": max(cpu_values), "sample_count": len(samples),
        "summary": expected_summary,
        "validation_sha256": sha256_bytes(validation_raw),
        "raw_sha256": sha256_bytes(raw_csv),
    }


def validate_sampler_signal_chronology(
    started_ns: int, finished_ns: int, exit_observed_ns: Optional[int],
    signal_not_before_ns: Optional[int], signal_returned_ns: Optional[int],
) -> None:
    if signal_not_before_ns is None and signal_returned_ns is None:
        return
    if (
        type(started_ns) is not int or type(finished_ns) is not int
        or type(exit_observed_ns) is not int
        or type(signal_not_before_ns) is not int
        or type(signal_returned_ns) is not int
        or not started_ns <= signal_not_before_ns <= finished_ns <= exit_observed_ns
        or not signal_not_before_ns <= signal_returned_ns <= exit_observed_ns
    ):
        fail("sampler root SIGTERM chronology differs")


def validate_sampler_terminal(
    authority: SamplerAuthority, prepared: PreparedExecution,
    returncode: int, *, required_start_before_ns: Optional[int] = None,
    required_finish_after_ns: Optional[int] = None,
    exit_observed_ns: Optional[int] = None,
    signal_not_before_ns: Optional[int] = None,
    signal_returned_ns: Optional[int] = None,
) -> Tuple[Dict[str, Any], Dict[str, bytes], Dict[str, Dict[str, Any]]]:
    if returncode != 0:
        fail("sampler terminal return code differs")
    raw_receipt = read_held_file(
        authority.held["receipt"], 1024 * 1024,
        "sampler terminal receipt",
    )
    receipt = parse_canonical_document(
        raw_receipt, "sampler terminal receipt", 1024 * 1024,
    )
    digest = receipt.get("self_sha256_excluding_field")
    if type(digest) is not str or LOWER64.fullmatch(digest) is None:
        fail("sampler terminal self-hash is absent")
    preimage = dict(receipt)
    del preimage["self_sha256_excluding_field"]
    if sha256_bytes(canonical_bytes(preimage) + b"\n") != digest:
        fail("sampler terminal self-hash differs")
    expected_keys = {
        "argv", "cpu_tctl_max_c", "edac_ce_paths", "edac_ue_paths",
        "errors", "expected_output_owner_uid", "exit_code",
        "finished_monotonic_ns", "finished_utc", "outcome", "pid",
        "pid_file", "raw_csv", "raw_samples_preserved", "receipt_file",
        "sampler_source", "sampling", "schema",
        "self_sha256_excluding_field", "signal", "started_monotonic_ns",
        "started_utc", "summary", "thresholds", "validation_jsonl",
    }
    started = receipt.get("started_monotonic_ns")
    finished = receipt.get("finished_monotonic_ns")
    started_utc_ns = sampler_utc_ns(
        receipt.get("started_utc"), "sampler terminal start",
    )
    finished_utc_ns = sampler_utc_ns(
        receipt.get("finished_utc"), "sampler terminal finish",
    )
    ce_controllers = sampler_edac_controller_roster(
        receipt.get("edac_ce_paths"), "ce", "sampler terminal",
    )
    ue_controllers = sampler_edac_controller_roster(
        receipt.get("edac_ue_paths"), "ue", "sampler terminal",
    )
    validate_sampler_signal_chronology(
        started, finished, exit_observed_ns,
        signal_not_before_ns, signal_returned_ns,
    )
    if (
        set(receipt) != expected_keys
        or receipt.get("schema") != SAMPLER_PRODUCER_SCHEMA
        or receipt.get("pid") != authority.process.pid
        or receipt.get("outcome") != "stopped"
        or receipt.get("exit_code") != 0
        or receipt.get("signal") != "SIGTERM"
        or receipt.get("errors") != []
        or receipt.get("expected_output_owner_uid") != CAMPAIGN_UID
        or receipt.get("argv") != sampler_argv(prepared)[5:]
        or receipt.get("raw_samples_preserved") is not True
        or receipt.get("sampling") != SAMPLER_SAMPLING
        or receipt.get("thresholds") != SAMPLER_THRESHOLDS
        or type(started) is not int or type(finished) is not int
        or started < 0 or started > finished
        or started_utc_ns > finished_utc_ns
        or ce_controllers != ue_controllers
        or (required_start_before_ns is not None
            and started > required_start_before_ns)
        or (required_finish_after_ns is not None
            and finished < required_finish_after_ns)
        or (exit_observed_ns is not None and finished > exit_observed_ns)
    ):
        fail("sampler terminal identity/outcome differs")
    source = receipt.get("sampler_source")
    source_digest, source_info = hash_path(
        prepared.sampler_script, MAX_DOCUMENT_BYTES,
        "sampler terminal source", expected_nlink=1,
    )
    source_binding = {
        "device": source_info.st_dev, "gid": source_info.st_gid,
        "inode": source_info.st_ino, "mode": stat.S_IMODE(source_info.st_mode),
        "nlink": source_info.st_nlink, "sha256": source_digest,
        "size": source_info.st_size, "uid": source_info.st_uid,
    }
    if (
        type(source) is not dict
        or set(source) != {
            "binding", "binding_finished", "destination", "expected_sha256",
            "path", "sha256",
        }
        or source.get("path") != str(prepared.sampler_script)
        or source.get("expected_sha256") != prepared.sampler_sha256
        or source.get("sha256") != prepared.sampler_sha256
        or source.get("binding") != source_binding
        or source.get("binding_finished") != source_binding
        or type(source.get("destination")) is not dict
    ):
        fail("sampler terminal source authority differs")
    validate_sampler_destination(
        source["destination"], prepared.sampler_script, None,
        os.stat(str(prepared.sampler_script.parent), follow_symlinks=False),
        "sampler source",
    )
    paths = {
        "csv": SAMPLER_CSV, "validation": SAMPLER_VALIDATION,
        "receipt": SAMPLER_RECEIPT,
    }
    observed: Dict[str, Dict[str, Any]] = {}
    raw: Dict[str, bytes] = {}
    for name, path in paths.items():
        data = read_held_file(
            authority.held[name], 16 * 1024 * 1024,
            "sampler terminal " + name,
        )
        binding = file_binding(
            authority.held[name], 16 * 1024 * 1024,
            "sampler terminal " + name,
        )
        named = os.stat(str(path), follow_symlinks=False)
        held = os.fstat(authority.held[name])
        if (
            not same_inode(named, held) or binding["mode"] != 0o444
            or binding["uid"] != CAMPAIGN_UID
            or binding["gid"] != CAMPAIGN_GID or binding["nlink"] != 1
        ):
            fail("sampler terminal artifact policy differs: " + name)
        raw[name] = data
        observed[name] = binding
    raw_csv = receipt.get("raw_csv")
    validation = receipt.get("validation_jsonl")
    receipt_file = receipt.get("receipt_file")
    if (
        type(raw_csv) is not dict
        or set(raw_csv) != {"binding", "destination", "path"}
        or raw_csv.get("path") != str(SAMPLER_CSV)
        or raw_csv.get("binding")
        != sampler_binding_from_observed(observed["csv"])
        or type(raw_csv.get("destination")) is not dict
        or type(validation) is not dict
        or set(validation) != {"binding", "destination", "path"}
        or validation.get("path") != str(SAMPLER_VALIDATION)
        or validation.get("binding")
        != sampler_binding_from_observed(observed["validation"])
        or type(validation.get("destination")) is not dict
        or type(receipt_file) is not dict
        or set(receipt_file) != {"destination", "path", "reservation_binding"}
        or receipt_file.get("path") != str(SAMPLER_RECEIPT)
        or receipt_file.get("reservation_binding")
        != sampler_binding_from_observed(
            authority.admission["artifacts"]["receipt"])
        or type(receipt_file.get("destination")) is not dict
    ):
        fail("sampler terminal artifact receipt differs")
    sampler_parent_info = os.stat(str(SAMPLER_DIR), follow_symlinks=False)
    validate_sampler_destination(
        raw_csv["destination"], SAMPLER_CSV, CAMPAIGN_UID,
        sampler_parent_info, "sampler raw CSV",
    )
    validate_sampler_destination(
        validation["destination"], SAMPLER_VALIDATION, CAMPAIGN_UID,
        sampler_parent_info, "sampler validation",
    )
    validate_sampler_destination(
        receipt_file["destination"], SAMPLER_RECEIPT, CAMPAIGN_UID,
        sampler_parent_info, "sampler receipt file",
    )
    for name in ("csv", "validation"):
        admitted = authority.admission["artifacts"][name]
        if any(observed[name][field] != admitted[field] for field in (
            "device", "gid", "inode", "nlink", "uid",
        )):
            fail("sampler terminal/admission identity differs: " + name)
    if (
        authority.admission["artifacts"]["csv"]["mode"] != 0o600
        or authority.admission["artifacts"]["validation"]["mode"] != 0o444
        or sampler_binding_from_observed(
            authority.admission["artifacts"]["receipt"])["size"] != 0
        or observed["receipt"]["device"]
        != authority.admission["artifacts"]["receipt"]["device"]
        or observed["receipt"]["inode"]
        != authority.admission["artifacts"]["receipt"]["inode"]
        or observed["receipt"]["mode"] != 0o444
        or observed["receipt"]["uid"] != CAMPAIGN_UID
        or observed["receipt"]["gid"] != CAMPAIGN_GID
        or observed["receipt"]["nlink"] != 1
        or observed["receipt"]["sha256"] != sha256_bytes(raw_receipt)
        or observed["receipt"]["bytes"] != len(raw_receipt)
        or len(raw["csv"]) < authority.admission["csv_bytes"]
        or sha256_bytes(raw["csv"][:authority.admission["csv_bytes"]])
        != authority.admission["csv_sha256"]
        or len(raw["validation"])
        < authority.admission["validation_jsonl"]["bytes"]
        or sha256_bytes(raw["validation"][:
            authority.admission["validation_jsonl"]["bytes"]])
        != authority.admission["validation_jsonl"]["sha256"]
    ):
        fail("sampler initial stream/reservation policy differs")
    pid_data = read_held_file(
        authority.held["pid"], 4096, "sampler removed PID file",
        expected_nlink=0,
    )
    pid_binding = file_binding(
        authority.held["pid"], 4096, "sampler removed PID file",
        expected_nlink=0,
    )
    pid_record = receipt.get("pid_file")
    if (
        pid_data != (str(authority.process.pid) + "\n").encode("ascii")
        or pid_binding["mode"] != 0o444
        or pid_binding["uid"] != CAMPAIGN_UID
        or pid_binding["gid"] != CAMPAIGN_GID
        or type(pid_record) is not dict
        or set(pid_record) != {"binding", "destination", "path", "removed"}
        or pid_record.get("removed") is not True
        or pid_record.get("path") != str(SAMPLER_PID_FILE)
        or pid_record.get("binding") != sampler_binding_from_observed(
            authority.admission["artifacts"]["pid"])
        or type(pid_record.get("destination")) is not dict
    ):
        fail("sampler removed PID authority differs")
    validate_sampler_destination(
        pid_record["destination"], SAMPLER_PID_FILE, CAMPAIGN_UID,
        sampler_parent_info, "sampler PID file",
    )
    expected_unlinked = sampler_binding_from_observed(
        authority.admission["artifacts"]["pid"],
    )
    expected_unlinked["nlink"] = 0
    if sampler_binding_from_observed(pid_binding) != expected_unlinked:
        fail("sampler held PID unlink transition differs")
    try:
        os.stat(str(SAMPLER_PID_FILE), follow_symlinks=False)
    except FileNotFoundError:
        pass
    else:
        fail("sampler PID pathname survived terminalization")
    raw["pid"] = pid_data
    observed["pid"] = pid_binding
    observed["stream_replay"] = validate_sampler_stream_replay(
        raw["csv"], raw["validation"], receipt,
    )
    return receipt, raw, observed


def _health_stat_identity(value: Any, expected: Mapping[str, Any],
                          where: str, *, mutable: bool) -> None:
    if type(value) is not dict or set(value) != {
        "device", "gid", "inode", "mode", "mtime_ns", "nlink", "size", "uid",
    }:
        fail(where + " stat keys differ")
    fields = ("device", "gid", "inode", "mode", "nlink", "uid")
    if any(value.get(name) != expected.get(name) for name in fields):
        fail(where + " stat identity differs")
    if not mutable and value != dict(expected):
        fail(where + " immutable stat differs")


def validate_health_sampler_prefixes(
    csv_prefix: bytes, validation_prefix: bytes, where: str,
) -> Dict[str, Any]:
    """Validate one exact paired, complete, continuing sampler prefix."""
    if (
        not csv_prefix.endswith(b"\n") or b"\r" in csv_prefix
        or not validation_prefix.endswith(b"\n") or b"\r" in validation_prefix
    ):
        fail(where + " sampler prefix contains a partial/non-LF record")
    raw_lines = csv_prefix.splitlines(keepends=True)
    validation_lines = validation_prefix.splitlines(keepends=True)
    if (
        len(raw_lines) < 2 or len(validation_lines) != len(raw_lines)
        or tuple(parse_sampler_csv_line(raw_lines[0], where + " raw header"))
        != SAMPLER_COLUMNS
    ):
        fail(where + " sampler prefix framing/pairing differs")
    header = parse_canonical_document(
        validation_lines[0], where + " validation header", 256 * 1024,
    )
    if (
        type(header) is not dict
        or header.get("schema") != SAMPLER_VALIDATION_STREAM_SCHEMA
        or header.get("raw_columns") != list(SAMPLER_COLUMNS)
    ):
        fail(where + " sampler validation header differs")
    sample_ns: List[int] = []
    previous_raw_time: Optional[float] = None
    for index, (raw_line, validation_line) in enumerate(zip(
        raw_lines[1:], validation_lines[1:],
    )):
        row = parse_sampler_csv_line(
            raw_line, "{} raw sample {}".format(where, index),
        )
        if len(row) != len(SAMPLER_COLUMNS):
            fail(where + " sampler raw-prefix width differs")
        sample = parse_sampler_validation_sample(
            parse_canonical_document(
                validation_line,
                "{} validation sample {}".format(where, index),
                256 * 1024,
            ),
            "{} validation sample {}".format(where, index),
        )
        raw_time = sampler_monotonic_text(
            row[1], "{} raw sample {} time".format(where, index),
        )
        if (
            sample["sample_index"] != index
            or sample["decision"] != "continue"
            or f"{sample['monotonic_s']:.6f}" != row[1]
            or previous_raw_time is not None and raw_time <= previous_raw_time
        ):
            fail(where + " sampler paired prefix chronology differs")
        previous_raw_time = raw_time
        sample_ns.append(int(round(raw_time * 1_000_000_000.0)))
    return {
        "raw_lines": raw_lines, "sample_ns": sample_ns,
        "validation_lines": validation_lines,
    }


def select_exact_health_window(
    terminal_csv: bytes, terminal_validation: bytes,
    start_ns: int, end_ns: int, where: str,
) -> Tuple[bytes, bytes]:
    """Select the unique contiguous paired terminal slice at two endpoints."""
    prefix = validate_health_sampler_prefixes(
        terminal_csv, terminal_validation, where,
    )
    timestamps = prefix["sample_ns"]
    start_matches = [
        index for index, value in enumerate(timestamps) if value == start_ns
    ]
    end_matches = [
        index for index, value in enumerate(timestamps) if value == end_ns
    ]
    if (
        len(start_matches) != 1 or len(end_matches) != 1
        or start_matches[0] > end_matches[0]
    ):
        fail(where + " sampler window endpoints differ")
    start = start_matches[0] + 1
    end = end_matches[0] + 2
    raw_lines = prefix["raw_lines"]
    validation_lines = prefix["validation_lines"]
    return (
        raw_lines[0] + b"".join(raw_lines[start:end]),
        validation_lines[0] + b"".join(validation_lines[start:end]),
    )


def validate_health_sampler_attestation(
    value: Any, authority: SamplerAuthority, prepared: PreparedExecution,
    terminal_csv: bytes, terminal_validation: bytes, where: str,
    expected_start_ns: Optional[int], expected_end_ns: Optional[int],
) -> Dict[str, Any]:
    if type(value) is not dict or set(value) != HEALTH_SAMPLER_KEYS:
        fail(where + " sampler attestation key roster differs")
    admission = authority.admission
    if (
        value.get("schema") != HEALTH_SAMPLER_ATTESTATION_SCHEMA
        or value.get("cmdline_argv") != list(authority.process.argv)
        or value.get("cmdline_sha256") != admission["cmdline_sha256"]
        or value.get("cpu") != SAMPLER_CPU
        or value.get("csv_device") != admission["artifacts"]["csv"]["device"]
        or value.get("csv_inode") != admission["artifacts"]["csv"]["inode"]
        or value.get("csv_path") != str(SAMPLER_CSV)
        or type(value.get("csv_bytes")) is not int
        or not 1 <= value["csv_bytes"] <= len(terminal_csv)
        or sha256_bytes(terminal_csv[:value["csv_bytes"]])
        != value.get("csv_sha256")
        or value.get("evidence_parent") != admission["evidence_parent"]
        or value.get("environ_sha256") != admission["environ_sha256"]
        or value.get("environment") != SEALED_ENVIRONMENT
        or value.get("environment_sha256")
        != sha256_bytes(canonical_bytes(SEALED_ENVIRONMENT))
        or value.get("executable_device")
        != admission["executable_device"]
        or value.get("executable_inode")
        != admission["executable_inode"]
        or value.get("executable_path") != str(prepared.python_image)
        or value.get("executable_sha256") != prepared.python_sha256
        or value.get("executable_stat") != admission["executable_stat"]
        or value.get("pid") != authority.process.pid
        or value.get("process_affinity") != [SAMPLER_CPU]
        or value.get("process_gid") != CAMPAIGN_GID
        or value.get("process_start_ticks") != authority.process.start_ticks
        or value.get("process_uid") != CAMPAIGN_UID
        or value.get("script_device") != admission["script_device"]
        or value.get("script_inode") != admission["script_inode"]
        or value.get("script_path") != str(prepared.sampler_script)
        or value.get("script_sha256") != prepared.sampler_sha256
        or value.get("script_stat") != admission["script_stat"]
        or value.get("terminal_status") != "ok"
        or type(value.get("window_start_monotonic_ns")) is not int
        or type(value.get("window_end_monotonic_ns")) is not int
        or value["window_start_monotonic_ns"] > value["window_end_monotonic_ns"]
        or expected_start_ns is not None
        and value["window_start_monotonic_ns"] != expected_start_ns
        or expected_end_ns is not None
        and value["window_end_monotonic_ns"] != expected_end_ns
        or value.get("validation_header_ascii")
        != admission["validation_header_ascii"]
    ):
        fail(where + " sampler identity/prefix authority differs")
    security = value.get("process_security")
    if (
        type(security) is not dict
        or set(security) != {
            "cap_ambient", "cap_bounding", "cap_effective",
            "cap_inheritable", "cap_permitted", "gids", "groups",
            "no_new_privs", "schema", "uids",
        }
        or security.get("uids") != [CAMPAIGN_UID] * 4
        or security.get("gids") != [CAMPAIGN_GID] * 4
        or security.get("groups") != [I2C_GID]
        or security.get("no_new_privs") != 1
        or security.get("schema")
        != "wirehair.wh2.direct-systematic-complement-process-security.v1"
        or any(security.get(name) != "0000000000000000" for name in (
            "cap_ambient", "cap_bounding", "cap_effective",
            "cap_inheritable", "cap_permitted",
        ))
    ):
        fail(where + " sampler process security differs")
    _health_stat_identity(
        value["csv_stat"], admission["csv_stat"], where + " CSV", mutable=True,
    )
    if (
        value["csv_stat"]["size"] < value["csv_bytes"]
        or value["csv_stat"]["size"] > len(terminal_csv)
    ):
        fail(where + " sampler CSV stat/prefix size differs")
    for name in ("pid_file", "receipt_file"):
        child = value.get(name)
        root = admission[name]
        if (
            type(child) is not dict or set(child) != {"bytes", "path", "sha256", "stat"}
            or child.get("bytes") != root["bytes"]
            or child.get("path") != root["path"]
            or child.get("sha256") != root["sha256"]
        ):
            fail(where + " sampler " + name + " differs")
        _health_stat_identity(
            child["stat"], root["stat"], where + " " + name, mutable=False,
        )
    validation = value.get("validation_jsonl")
    if (
        type(validation) is not dict
        or set(validation) != {"bytes", "path", "sha256", "stat"}
        or validation.get("path") != str(SAMPLER_VALIDATION)
        or type(validation.get("bytes")) is not int
        or not 1 <= validation["bytes"] <= len(terminal_validation)
        or sha256_bytes(terminal_validation[:validation["bytes"]])
        != validation.get("sha256")
    ):
        fail(where + " sampler validation prefix differs")
    _health_stat_identity(
        validation["stat"], admission["validation_jsonl"]["stat"],
        where + " validation", mutable=True,
    )
    if (
        validation["stat"]["size"] < validation["bytes"]
        or validation["stat"]["size"] > len(terminal_validation)
    ):
        fail(where + " sampler validation stat/prefix size differs")
    if expected_start_ns == 0 and expected_end_ns == 0:
        root_validation = admission["validation_jsonl"]
        if (
            value["csv_bytes"] < admission["csv_bytes"]
            or validation["bytes"] < root_validation["bytes"]
            or sha256_bytes(terminal_csv[:admission["csv_bytes"]])
            != admission["csv_sha256"]
            or sha256_bytes(terminal_validation[:root_validation["bytes"]])
            != root_validation["sha256"]
            or value["csv_stat"]["size"] < admission["csv_stat"]["size"]
            or value["csv_stat"]["mtime_ns"]
            < admission["csv_stat"]["mtime_ns"]
            or validation["stat"]["size"] < root_validation["stat"]["size"]
            or validation["stat"]["mtime_ns"]
            < root_validation["stat"]["mtime_ns"]
        ):
            fail(where + " sampler admission did not grow from root admission")
    prefix = validate_health_sampler_prefixes(
        terminal_csv[:value["csv_bytes"]],
        terminal_validation[:validation["bytes"]], where,
    )
    if prefix["validation_lines"][0] != value[
        "validation_header_ascii"
    ].encode("ascii"):
        fail(where + " sampler validation-header bytes differ")
    return prefix


def validate_health_sampler_growth(initial: Mapping[str, Any],
                                   final: Mapping[str, Any]) -> None:
    immutable = set(HEALTH_SAMPLER_KEYS) - {
        "csv_bytes", "csv_sha256", "csv_stat", "validation_jsonl",
        "window_start_monotonic_ns", "window_end_monotonic_ns",
    }
    if any(initial.get(name) != final.get(name) for name in immutable):
        fail("valid child sampler live authority changed across worker interval")
    initial_validation = initial["validation_jsonl"]
    final_validation = final["validation_jsonl"]
    if (
        final["csv_bytes"] < initial["csv_bytes"]
        or final_validation["bytes"] < initial_validation["bytes"]
        or final["csv_stat"]["size"] < initial["csv_stat"]["size"]
        or final["csv_stat"]["mtime_ns"] < initial["csv_stat"]["mtime_ns"]
        or final_validation["stat"]["size"]
        < initial_validation["stat"]["size"]
        or final_validation["stat"]["mtime_ns"]
        < initial_validation["stat"]["mtime_ns"]
    ):
        fail("valid child sampler live streams did not grow monotonically")


def validate_child_health_against_root(
    bundle: ChildBundle, authority: SamplerAuthority,
    prepared: PreparedExecution, root_sampler_terminal: Mapping[str, Any],
) -> None:
    if bundle.summary is None or bundle.provenance is None:
        fail("child health/root crosscheck lacks valid evidence")
    health = bundle.provenance["health"]
    if type(health) is not dict or set(health) != HEALTH_KEYS:
        fail("valid child health key roster differs")
    if health != bundle.summary.get("health"):
        fail("valid child summary/provenance health differs")
    health_preimage = dict(health)
    health_digest = health_preimage.get("receipt_sha256")
    health_preimage["receipt_sha256"] = None
    if (
        health.get("schema")
        != "wirehair.wh2.direct-systematic-complement-health.v3"
        or health.get("evidence_status") != "complete"
        or health.get("terminal_status") != "ok"
        or health.get("collection_failures") != []
        or health.get("violations") != []
        or health.get("target_cpu") != WORKER_CPU
        or health.get("controller_cpu") != CONTROLLER_CPU
        or health.get("sampler_cpu") != SAMPLER_CPU
        or health.get("controller_initial_affinity")
        != list(HEALTH_CONTROLLER_INITIAL_AFFINITY)
        or health.get("controller_singleton_affinity_end") != [CONTROLLER_CPU]
        or health.get("target_core") != [0, 56]
        or health.get("controller_core") != [0, 57]
        or health.get("sampler_core") != [0, 58]
        or health.get("target_threads") != [SIBLING_CPU, WORKER_CPU]
        or health.get("edac_policy") != "ce-and-ue-every-sample-zero-v1"
        or health.get("thermal_max_millic") != 85000
        or health.get("sibling_non_idle_tick_cap") != 5
        or health.get("sibling_tick_policy")
        != "linux-proc-stat-user-nice-system-irq-softirq-steal-v1"
        or health.get("sampler_terminal") != {}
        or health.get("sampler_terminal_receipt_sha256") is not None
        or sha256_bytes(canonical_bytes(health_preimage)) != health_digest
    ):
        fail("valid child health root policy differs")
    child_start = health.get("child_start_monotonic_ns")
    child_reap = health.get("child_reap_monotonic_ns")
    if type(child_start) is not int or type(child_reap) is not int or child_start > child_reap:
        fail("valid child health interval differs")
    siblings = health.get("sibling_ticks")
    admission_sibling = health.get("admission_sibling_ticks")
    if (
        type(siblings) is not list or len(siblings) != 1
        or type(siblings[0]) is not dict
        or siblings[0].get("cpu") != SIBLING_CPU
        or type(siblings[0].get("delta_non_idle_ticks")) is not int
        or not 0 <= siblings[0]["delta_non_idle_ticks"] <= 5
        or siblings[0].get("interval_start_monotonic_ns") != child_start
        or siblings[0].get("interval_end_monotonic_ns") != child_reap
        or type(admission_sibling) is not dict
        or admission_sibling.get("cpu") != SIBLING_CPU
    ):
        fail("valid child sibling-tick health differs")
    terminal_csv = read_held_file(
        authority.held["csv"], 16 * 1024 * 1024,
        "root health terminal CSV",
    )
    terminal_validation = read_held_file(
        authority.held["validation"], 16 * 1024 * 1024,
        "root health terminal validation",
    )
    admission = health.get("sampler_admission")
    sampler = health.get("sampler")
    admission_prefix = validate_health_sampler_attestation(
        admission, authority, prepared, terminal_csv, terminal_validation,
        "child sampler admission", 0, 0,
    )
    sampler_prefix = validate_health_sampler_attestation(
        sampler, authority, prepared, terminal_csv, terminal_validation,
        "child sampler terminal-live", None, None,
    )
    validate_health_sampler_growth(admission, sampler)
    window_start = sampler["window_start_monotonic_ns"]
    window_end = sampler["window_end_monotonic_ns"]
    if not window_start <= child_start <= child_reap <= window_end:
        fail("valid child sampler window does not cover worker interval")
    # Both live attestations must end at complete, paired prefixes of the
    # retained terminal streams.  The final claimed endpoints must already be
    # present in that final prefix, rather than referring to later samples.
    select_exact_health_window(
        b"".join(sampler_prefix["raw_lines"]),
        b"".join(sampler_prefix["validation_lines"]),
        window_start, window_end, "child sampler terminal-live prefix",
    )
    if (
        len(admission_prefix["raw_lines"]) > len(sampler_prefix["raw_lines"])
        or len(admission_prefix["validation_lines"])
        > len(sampler_prefix["validation_lines"])
    ):
        fail("valid child sampler paired prefixes did not grow")
    if (
        health.get("sampler_admission_receipt_sha256")
        != sha256_bytes(canonical_bytes(admission))
        or health.get("sampler_receipt_sha256")
        != sha256_bytes(canonical_bytes(sampler))
    ):
        fail("valid child sampler health self-binding differs")
    thermal = health.get("thermal")
    if type(thermal) is not dict or set(thermal) != HEALTH_THERMAL_KEYS:
        fail("valid child thermal receipt key roster differs")
    try:
        window_csv = thermal["window_csv_ascii"].encode("ascii")
        window_validation = thermal["validation_jsonl_ascii"].encode("ascii")
    except (AttributeError, UnicodeEncodeError):
        fail("valid child thermal window encoding differs")
    raw_lines = window_csv.splitlines(keepends=True)
    validation_lines = window_validation.splitlines(keepends=True)
    full_raw_lines = terminal_csv.splitlines(keepends=True)
    full_validation_lines = terminal_validation.splitlines(keepends=True)
    if (
        len(raw_lines) < 3 or len(validation_lines) != len(raw_lines)
        or tuple(parse_sampler_csv_line(raw_lines[0], "health raw header"))
        != SAMPLER_COLUMNS
        or raw_lines[0] != full_raw_lines[0]
        or validation_lines[0] != full_validation_lines[0]
    ):
        fail("valid child thermal window framing differs")
    expected_window_csv, expected_window_validation = select_exact_health_window(
        terminal_csv, terminal_validation, window_start, window_end,
        "valid child thermal root selection",
    )
    if (
        window_csv != expected_window_csv
        or window_validation != expected_window_validation
    ):
        fail("valid child thermal window is not the exact contiguous root slice")
    cpu_values: List[float] = []
    dimm_values: List[float] = []
    attempt_errors = 0
    first_ns: Optional[int] = None
    previous_ns: Optional[int] = None
    for raw_line, validation_line in zip(raw_lines[1:], validation_lines[1:]):
        row = parse_sampler_csv_line(raw_line, "health raw row")
        sample = parse_sampler_validation_sample(
            parse_canonical_document(validation_line, "health validation row", 256 * 1024),
            "health validation row",
        )
        timestamp = row[1]
        sample_ns = int(round(float(timestamp) * 1_000_000_000.0))
        if (
            f"{float(sample['monotonic_s']):.6f}" != timestamp
            or not window_start <= sample_ns <= window_end
            or previous_ns is not None
            and (sample_ns <= previous_ns or sample_ns - previous_ns > 5_000_000_000)
        ):
            fail("valid child thermal window/root stream chronology differs")
        previous_ns = sample_ns
        if first_ns is None:
            first_ns = sample_ns
        cpu_values.append(sampler_finite_text(row[4], "health CPU Tctl"))
        dimm_values.extend(
            sampler_finite_text(row[index], "health DIMM")
            for index in range(5, 13)
        )
        if (
            sampler_counter_text(row[13], "health DIMM errors") != 0
            or sampler_counter_text(row[17], "health EDAC CE") != 0
            or sampler_counter_text(row[18], "health EDAC UE") != 0
            or sample["decision"] != "continue"
            or sample["fault_count"] != 0 or sample["read_error_count"] != 0
        ):
            fail("valid child thermal window contains a policy violation")
        attempt_errors += sum(
            sensor["attempt_errors"] for sensor in sample["sensors"].values()
        )
    sample_count = len(raw_lines) - 1
    if (
        thermal.get("cpu") != SAMPLER_CPU
        or thermal.get("pid") != authority.process.pid
        or thermal.get("process_start_ticks") != authority.process.start_ticks
        or thermal.get("script_path") != str(prepared.sampler_script)
        or thermal.get("script_sha256") != prepared.sampler_sha256
        or thermal.get("csv_device") != authority.admission["artifacts"]["csv"]["device"]
        or thermal.get("csv_inode") != authority.admission["artifacts"]["csv"]["inode"]
        or thermal.get("csv_path") != str(SAMPLER_CSV)
        or thermal.get("validation_device")
        != authority.admission["artifacts"]["validation"]["device"]
        or thermal.get("validation_inode")
        != authority.admission["artifacts"]["validation"]["inode"]
        or thermal.get("validation_path") != str(SAMPLER_VALIDATION)
        or thermal.get("terminal_status") != "complete"
        or first_ns != window_start or previous_ns != window_end
        or thermal.get("window_start_monotonic_ns") != window_start
        or thermal.get("window_end_monotonic_ns") != window_end
        or thermal.get("window_csv_bytes") != len(window_csv)
        or thermal.get("window_csv_sha256") != sha256_bytes(window_csv)
        or thermal.get("validation_jsonl_bytes") != len(window_validation)
        or thermal.get("validation_jsonl_sha256") != sha256_bytes(window_validation)
        or thermal.get("sample_count") != sample_count
        or thermal.get("validation_sample_count") != sample_count
        or thermal.get("valid_sample_count") != sample_count
        or thermal.get("invalid_sample_count") != 0
        or thermal.get("parse_failures") != []
        or thermal.get("validation_failures") != []
        or thermal.get("validation_attempt_errors_total") != attempt_errors
        or thermal.get("cpu_tctl_max_millic")
        != int(round(max(cpu_values) * 1000.0))
        or thermal.get("dimm_max_millic")
        != int(round(max(dimm_values) * 1000.0))
        or thermal.get("cpu_tctl_max_millic") > 85000
        or thermal.get("dimm_max_millic") > 85000
        or thermal.get("dimm_read_errors") != 0
        or thermal.get("edac_ce_max") != 0
        or thermal.get("edac_ue_max") != 0
        or root_sampler_terminal.get("terminal", {}).get("outcome") != "stopped"
    ):
        fail("valid child thermal/root terminal cross-binding differs")


@dataclass
class SamplerTerminalState:
    signal_not_before_ns: Optional[int] = None
    signal_returned_ns: Optional[int] = None
    signal_call_entered: bool = False
    exit_observed_ns: Optional[int] = None
    communication_started: bool = False
    communication: Optional[Tuple[bytes, bytes, int]] = None
    proof: Optional[Tuple[Dict[str, Any], Dict[str, bytes],
                          Dict[str, Any]]] = None
    copied: Dict[str, Dict[str, Any]] = field(default_factory=dict)
    directory_sealed: bool = False
    result: Optional[Dict[str, Any]] = None


def stop_and_validate_sampler(
    authority: SamplerAuthority, prepared: PreparedExecution,
    sampler_directory: BoundDirectory, journal: AttemptJournal,
    *, required_start_before_ns: Optional[int] = None,
    required_finish_after_ns: Optional[int] = None,
    owner: Optional[OwnershipSlot] = None,
    fault_injector: Optional[Any] = None,
) -> Dict[str, Any]:
    state: Optional[SamplerTerminalState] = None
    if owner is not None and owner.value is not None:
        if not isinstance(owner.value, SamplerTerminalState):
            fail("sampler terminal ownership type differs")
        state = owner.value
    if state is None:
        state = SamplerTerminalState()
        if owner is not None:
            try:
                owner.adopt(state)
            except BaseException:
                if not owner.owns(state):
                    raise
                raise
    if state.result is not None:
        return state.result

    if state.proof is None:
        if state.communication_started and state.communication is None:
            # communicate_process drains nonblocking pipes incrementally.  An
            # exception can therefore consume an unauthorised output prefix
            # even while the child remains live.  With no persistent byte
            # accumulator, every such interrupted attempt is irrecoverably
            # ambiguous and must never be retried into a clean transcript.
            fail("sampler communication result is ambiguous after start")
        if authority.process.returncode is None:
            if (
                process_exited(authority.process)
                or process_start_ticks(authority.process.pid)
                != authority.process.start_ticks
            ):
                fail("sampler was not live at the root SIGTERM boundary")
            if state.signal_returned_ns is None:
                # Store the lower bound before the irreversible syscall.  If
                # the syscall result itself is interrupted before assignment,
                # fail closed rather than infer successful root signaling.
                if state.signal_not_before_ns is None:
                    state.signal_not_before_ns = time.monotonic_ns()
                state.signal_call_entered = True
                pidfd_signal(authority.process.pidfd, signal.SIGTERM)
                state.signal_returned_ns = time.monotonic_ns()
            state.communication_started = True
            communication = communicate_process(
                authority.process, time.monotonic() + SAMPLER_STOP_SECONDS,
                MAX_STDERR_BYTES, MAX_STDERR_BYTES,
            )
            if fault_injector is not None:
                fault_injector("sampler-communicate-returned")
            state.communication = communication
        if state.communication is not None:
            stdout, stderr, communication_returncode = state.communication
            if communication_returncode != authority.process.returncode:
                fail("sampler communication returncode binding differs")
            if stdout or stderr:
                fail("sampler emitted unexpected stdout/stderr")
        if (
            state.signal_not_before_ns is None
            or not state.signal_call_entered
            or state.signal_returned_ns is None
        ):
            fail("sampler root SIGTERM result is ambiguous")
        if authority.process.returncode is None:
            fail("sampler was not reaped after root SIGTERM")
        if state.exit_observed_ns is None:
            state.exit_observed_ns = time.monotonic_ns()
        state.proof = validate_sampler_terminal(
            authority, prepared, authority.process.returncode,
            required_start_before_ns=required_start_before_ns,
            required_finish_after_ns=required_finish_after_ns,
            exit_observed_ns=state.exit_observed_ns,
            signal_not_before_ns=state.signal_not_before_ns,
            signal_returned_ns=state.signal_returned_ns,
        )
    assert state.proof is not None
    receipt, raw, bindings = state.proof
    publications = (
        ("csv", "sampler.thermal.csv", raw["csv"]),
        ("validation", "sampler.validation.jsonl", raw["validation"]),
        ("receipt", "sampler.terminal.json", raw["receipt"]),
        ("pid", "sampler.removed-pid.txt", raw["pid"]),
    )
    for key, name, payload in publications:
        if key not in state.copied:
            state.copied[key] = journal.write_bytes_resumable(name, payload)
        elif journal.write_bytes_resumable(name, payload) != state.copied[key]:
            fail("sampler copied artifact receipt changed: " + key)
    if not state.directory_sealed:
        sampler_directory.seal((
            SAMPLER_CSV.name, SAMPLER_VALIDATION.name, SAMPLER_RECEIPT.name,
        ))
        state.directory_sealed = True
    elif sampler_directory.state != "sealed":
        fail("sampler terminal directory seal state differs")
    result = {
        "bindings": bindings, "copied": dict(state.copied),
        "receipt_sha256": sha256_bytes(raw["receipt"]),
        "signal_not_before_monotonic_ns": state.signal_not_before_ns,
        "signal_returned_monotonic_ns": state.signal_returned_ns,
        "signal_proof": "pidfd-send-return-observed",
        "terminal": receipt,
    }
    state.result = result
    return result


def kill_cgroup_and_wait(path: Path, kill_fd: int,
                         seconds: float = 3.0, *,
                         observed_roster: Optional[List[int]] = None,
                         ) -> List[int]:
    roster = cgroup_descendant_pids(path)
    if observed_roster is not None:
        # Publish every helper-local observation into caller-owned persistent
        # state before the irreversible cgroup.kill write.  A return-handoff
        # exception therefore cannot lose a PID that this helper observed.
        for pid_value in roster:
            if pid_value not in observed_roster:
                observed_roster.append(pid_value)
        observed_roster.sort()
    if roster and os.write(kill_fd, b"1") != 1:
        fail("cgroup.kill write was short: " + str(path))
    deadline = time.monotonic() + seconds
    while cgroup_descendant_pids(path):
        if time.monotonic() >= deadline:
            fail("cgroup remained populated after kill: " + str(path))
        time.sleep(0.005)
    return roster


def strict_process_close(process: ProcessHandle) -> None:
    failures: List[str] = []
    if process.returncode is None:
        failures.append(process.role + " is not reaped")
    for member in ("pidfd", "stdout_fd", "stderr_fd"):
        fd = getattr(process, member)
        if fd < 0:
            continue
        try:
            close_object_fd(process, member, process.role + " " + member)
        except OSError as exc:
            failures.append(member + ": " + exception_text(exc))
        except BaseException as exc:
            failures.append(member + ": " + exception_text(exc))
    if failures:
        fail(process.role + " descriptor closure differs: " + " | ".join(failures))


def drain_and_reap_cancelled(blocked: BlockedRole, leaf: Path,
                             kill_fd: int,
                             before_kill: Optional[Any] = None,
                             ) -> Dict[str, Any]:
    if blocked.cancellation_result is not None:
        if cgroup_descendant_pids(leaf):
            fail(blocked.process.role + " cancelled cgroup repopulated")
        return dict(blocked.cancellation_result)
    if blocked.cancellation_roster is None:
        # Persist the exact owned stub before cancellation/exit can make its
        # cgroup membership disappear.  Later observations may only extend it.
        blocked.cancellation_roster = [blocked.process.pid]
    blocked.cancel()
    # A blocked bootstrap normally exits on gate EOF.  Give that exact child
    # a short bounded grace period so cancellation does not needlessly invoke
    # cgroup.kill while an armed guardian still owns the durable deadline.
    grace = time.monotonic() + 0.25
    while blocked.process.returncode is None and time.monotonic() < grace:
        if process_exited(blocked.process):
            wait_process(blocked.process, grace)
            break
        time.sleep(0.005)
    roster = cgroup_descendant_pids(leaf)
    for pid_value in roster:
        if pid_value not in blocked.cancellation_roster:
            blocked.cancellation_roster.append(pid_value)
    blocked.cancellation_roster.sort()
    if before_kill is not None:
        # The callback is also the resumable result-publication authority.  It
        # must run when a prior kill already emptied the leaf but an async
        # exception interrupted the corresponding journal result handoff.
        observed = before_kill()
        if observed is not None:
            if (
                type(observed) is not list
                or any(type(value) is not int or value <= 0 for value in observed)
            ):
                fail(blocked.process.role + " cancellation roster differs")
            for pid_value in observed:
                if pid_value not in blocked.cancellation_roster:
                    blocked.cancellation_roster.append(pid_value)
            blocked.cancellation_roster.sort()
    elif roster:
        kill_cgroup_and_wait(leaf, kill_fd)
    output, error, returncode = communicate_process(
        blocked.process, time.monotonic() + 2.0, 4096, MAX_STDERR_BYTES,
    )
    result = {
        "returncode": returncode,
        "roster": list(blocked.cancellation_roster),
        "stderr_sha256": sha256_bytes(error),
        "stdout_sha256": sha256_bytes(output),
    }
    blocked.cancellation_result = result
    return dict(result)


@dataclass
class ContainmentState:
    index: int
    path: str
    reason: str
    roster: List[int]
    requested_monotonic_ns: int
    observed_roster: List[int] = field(default_factory=list)
    intent_receipt: Optional[Dict[str, Any]] = None
    kill_may_have_been_issued: bool = False
    result_receipt: Optional[Dict[str, Any]] = None
    completed_monotonic_ns: Optional[int] = None

    def __post_init__(self) -> None:
        # `roster` is the immutable observation bound into the durable intent.
        # Later members are evidence about the same exact cgroup kill scope,
        # not grounds to refuse the atomic cgroup.kill operation.
        if not self.observed_roster:
            self.observed_roster = list(self.roster)


def reap_adopted_children(values: List[Dict[str, int]]) -> List[Dict[str, int]]:
    """Persist exited-child identity before the irreversible waitpid reap."""
    while True:
        try:
            observed = os.waitid(
                os.P_ALL, 0, os.WEXITED | os.WNOHANG | os.WNOWAIT,
            )
        except ChildProcessError:
            break
        if observed is None:
            fail("one or more adopted descendants remain live")
        pid = observed.si_pid
        code = observed.si_code
        status_value = observed.si_status
        if (
            type(pid) is not int or pid <= 0 or type(code) is not int
            or type(status_value) is not int or status_value < 0
            or code not in (1, 2, 3)
        ):
            fail("adopted descendant waitid receipt differs")
        returncode = status_value if code == 1 else -status_value
        receipt = {
            "pid": pid, "returncode": returncode,
            "waitid_code": code, "waitid_status": status_value,
        }
        matches = [item for item in values if item.get("pid") == pid]
        if matches and matches != [receipt]:
            fail("adopted descendant persisted identity changed")
        if not matches:
            values.append(receipt)
        try:
            waited, raw_status = os.waitpid(pid, os.WNOHANG)
        except ChildProcessError:
            # A Python exception may have landed after waitpid reaped this
            # already-persisted child.  Its WNOWAIT receipt remains authority.
            continue
        if waited != pid or wait_status_code(raw_status) != returncode:
            fail("adopted descendant reap receipt differs")
    values.sort(key=lambda item: item["pid"])
    return list(values)


def process_start_monotonic_ns(pid: int) -> int:
    ticks = process_start_ticks(pid)
    frequency = os.sysconf("SC_CLK_TCK")
    if type(frequency) is not int or frequency <= 0:
        fail("system clock-tick frequency differs")
    return ticks * 1_000_000_000 // frequency


@dataclass(frozen=True)
class RootRunResult:
    outcome: str
    launcher_status: str
    exit_code: int
    terminal_sha256: str


class ExecutionSession:
    """One-shot root lifecycle; all evidence publication is descriptor-bound."""

    @property
    def layout(self) -> Optional[CgroupLayout]:
        return self._layout_owner.value

    @property
    def journal(self) -> Optional[AttemptJournal]:
        return self._journal_owner.value

    @property
    def sampler_blocked(self) -> Optional[BlockedRole]:
        return self._sampler_blocked_owner.value

    @property
    def sampler(self) -> Optional[SamplerAuthority]:
        return self._sampler_owner.value

    @property
    def controller_blocked(self) -> Optional[BlockedRole]:
        return self._controller_blocked_owner.value

    @property
    def verifier_blocked(self) -> Optional[BlockedRole]:
        return self._verifier_blocked_owner.value

    @property
    def guardian(self) -> Optional[GuardianHandle]:
        return self._guardian_owner.value

    @property
    def bundle(self) -> Optional[ChildBundle]:
        return self._bundle_owner.value

    @property
    def root_t0_ns(self) -> Optional[int]:
        return None if self._root_window is None else self._root_window[0]

    @property
    def root_deadline_ns(self) -> Optional[int]:
        return None if self._root_window is None else self._root_window[1]

    def guardian_stream_fd(self, *, required: bool = True) -> int:
        """Recover the journal-owned guardian stream across handoff faults."""
        registered = -1
        if self.journal is not None:
            registered = self.journal.stream_fds.get("guardian.jsonl", -1)
        if self.guardian_event_fd >= 0:
            if registered >= 0 and registered != self.guardian_event_fd:
                fail("guardian event descriptor ownership differs")
            return self.guardian_event_fd
        if registered >= 0:
            self.guardian_event_fd = registered
            return registered
        if required:
            fail("guardian event descriptor is unavailable")
        return -1

    def __init__(self, config: LaunchConfig) -> None:
        self.config = config
        self.process_started_ns = process_start_monotonic_ns(os.getpid())
        self.initial_fds = open_fd_roster()
        self.initial_fd_authority = fd_authority_roster(self.initial_fds)
        self.prepared: Optional[PreparedExecution] = None
        self.lock: Optional[CampaignLock] = None
        self.var_tmp_fd = -1
        self.service_path: Optional[Path] = None
        self._layout_owner = OwnershipSlot("cgroup layout")
        self._journal_owner = OwnershipSlot("attempt journal")
        self.controller_directory: Optional[BoundDirectory] = None
        self.sampler_directory: Optional[BoundDirectory] = None
        self._sampler_owner = OwnershipSlot("sampler authority")
        self._sampler_blocked_owner = OwnershipSlot("sampler blocked role")
        self.sampler_process: Optional[ProcessHandle] = None
        self.sampler_may_have_started = False
        self._controller_blocked_owner = OwnershipSlot("controller blocked role")
        self._verifier_blocked_owner = OwnershipSlot("verifier blocked role")
        self.controller: Optional[ProcessHandle] = None
        self.verifier: Optional[ProcessHandle] = None
        self._guardian_owner = OwnershipSlot("root guardian")
        self.guardian_event_fd = -1
        self._bundle_owner = OwnershipSlot("controller child bundle")
        self.child_bundle_validated = False
        self.child_copy_receipts: Dict[str, Dict[str, Any]] = {}
        self.child_copies_complete = False
        self.monitor: Optional[WorkerMapMonitor] = None
        self.errors: List[str] = []
        self.phases: List[str] = []
        self.controller_released = False
        self.controller_communication_published = False
        self.controller_cancellation_published = False
        self.verifier_resolved = False
        self.sampler_authenticated = False
        self.sampler_admission_published = False
        self.sampler_cancellation_published = False
        self.sampler_admission_receipt: Optional[Dict[str, Any]] = None
        self.guardian_authenticated = False
        self.quiesced = False
        self.postflight_ok = False
        self.namespaces_sealed = False
        self.layout_removed = False
        self.runtime_closed = False
        self.controller_rc: Optional[int] = None
        self.controller_reap_ns: Optional[int] = None
        self.verifier_rc: Optional[int] = None
        self.verifier_reap_ns: Optional[int] = None
        self._root_window: Optional[Tuple[int, int]] = None
        self.child_outcome: Optional[str] = None
        self.replay: Optional[Dict[str, Any]] = None
        self.map_receipt: Optional[Dict[str, Any]] = None
        self.sampler_terminal: Optional[Dict[str, Any]] = None
        self._sampler_terminal_owner = OwnershipSlot("sampler terminal state")
        self.guardian_receipt: Optional[Dict[str, Any]] = None
        self.controller_security: Optional[Dict[str, Any]] = None
        self.verifier_security: Optional[Dict[str, Any]] = None
        self.containment_roster: List[int] = []
        self.containment_states: List[ContainmentState] = []
        self.containment_receipts: List[Dict[str, Any]] = []
        self.unexpected_final_containment = False
        self.adopted_reaps: List[Dict[str, int]] = []
        self.adopted_reaps_reported = False
        self.recovered_fds: List[int] = []

    def phase(self, name: str) -> None:
        self.phases.append(name)

    def record(self, where: str, exc: BaseException) -> None:
        message = where + ": " + exception_text(exc)
        if message not in self.errors and len(self.errors) < 64:
            self.errors.append(message[:4096])

    def sampler_release_possible(self) -> bool:
        """Return the conservative irreversible sampler-release state."""
        return bool(
            self.sampler_may_have_started
            or (
                self.sampler_blocked is not None
                and self.sampler_blocked.may_have_released
            )
        )

    def contain_cgroup(self, path: Path, kill_fd: int,
                       reason: str) -> List[int]:
        """Durably journal one parent containment while guardian stays live."""
        if self.layout is None or self.journal is None:
            fail("parent containment lacks runtime authority")
        expected = {
            str(self.layout.experiment): self.layout.experiment_kill_fd,
            str(self.layout.verifier): self.layout.verifier_kill_fd,
            str(self.layout.sampler): self.layout.sampler_kill_fd,
            str(self.layout.run): self.layout.run_kill_fd,
        }
        if (
            str(path) not in expected or expected[str(path)] != kill_fd
            or type(reason) is not str or not reason or len(reason) > 256
        ):
            fail("parent containment target authority differs")
        state = next((
            item for item in reversed(self.containment_states)
            if item.path == str(path) and (
                item.result_receipt is None
                or item.index >= len(self.containment_receipts)
            )
        ), None)
        if state is None:
            roster = cgroup_descendant_pids(path)
            if not roster:
                return []
            if len(self.containment_states) >= 32:
                fail("parent containment event bound differs")
            state = ContainmentState(
                index=len(self.containment_states), path=str(path),
                reason=reason, roster=roster,
                requested_monotonic_ns=time.monotonic_ns(),
            )
            self.containment_states.append(state)
        elif state.reason != reason:
            # A retry resumes the first irreversible reason; it cannot relabel
            # an already-journaled containment as a later cleanup phase.
            reason = state.reason
        if self._root_window is not None:
            if self.guardian is None:
                self.seal_failed_guardian_admission()
                if not self.guardian_authenticated:
                    fail("parent containment lacks guardian admission terminal")
        prefix = "containment-{:03d}".format(state.index)
        intent = self_hashed(
            "wirehair.wh2.v2-facade-timing-root-containment-intent.v2", {
                "attempt_id": ATTEMPT_ID,
                "guardian_deadline_monotonic_ns": self.root_deadline_ns,
                "guardian_pid": (
                    self.guardian.process.pid if self.guardian is not None
                    else None
                ),
                "initial_roster": list(state.roster),
                "kill_scope": (
                    "all-processes-in-exact-cgroup-and-descendants-at-"
                    "atomic-cgroup.kill"
                ),
                "path": state.path, "postcondition": "exact-cgroup-empty",
                "reason": state.reason,
                "requested_monotonic_ns": state.requested_monotonic_ns,
            },
        )
        intent_raw = canonical_bytes(intent) + b"\n"
        intent_receipt = self.journal.write_bytes_resumable(
            prefix + "-intent.json", intent_raw,
        )
        if (
            state.intent_receipt is not None
            and state.intent_receipt != intent_receipt
        ):
            fail("parent containment intent receipt changed")
        state.intent_receipt = intent_receipt
        if self.guardian is not None and not self.guardian_authenticated:
            try:
                self.guardian.poll()
            except LaunchError:
                # The guardian may reach its independent deadline while the
                # already-requested containment intent is being durably
                # fsynced.  Its terminal record precedes its whole-run kill.
                # Wait for that authority to empty the run, authenticate the
                # durable terminal, then finish this pending result instead
                # of retrying an impossible live-guardian poll forever.
                terminal_deadline = time.monotonic() + 2.0
                while cgroup_descendant_pids(self.layout.run):
                    if time.monotonic() >= terminal_deadline:
                        fail("terminal guardian did not empty the run")
                    time.sleep(0.01)
                self.resolve_guardian()
                if not self.guardian_authenticated:
                    fail("terminal guardian containment authority is absent")
        remaining = cgroup_descendant_pids(path)
        for pid_value in remaining:
            if pid_value not in state.observed_roster:
                state.observed_roster.append(pid_value)
        state.observed_roster.sort()
        if remaining:
            state.kill_may_have_been_issued = True
            kill_cgroup_and_wait(
                path, kill_fd, observed_roster=state.observed_roster,
            )
        if cgroup_descendant_pids(path):
            fail("parent containment target remained populated")
        if state.completed_monotonic_ns is None:
            state.completed_monotonic_ns = time.monotonic_ns()
        result = self_hashed(
            "wirehair.wh2.v2-facade-timing-root-containment-result.v2", {
                "attempt_id": ATTEMPT_ID,
                "completed_monotonic_ns": state.completed_monotonic_ns,
                "initial_roster": list(state.roster),
                "intent_sha256": intent_receipt["sha256"],
                "kill_may_have_been_issued": state.kill_may_have_been_issued,
                "kill_scope": (
                    "all-processes-in-exact-cgroup-and-descendants-at-"
                    "atomic-cgroup.kill"
                ),
                "observed_roster": list(state.observed_roster),
                "path": state.path, "reason": state.reason,
                "remaining": [],
            },
        )
        result_receipt = self.journal.write_bytes_resumable(
            prefix + "-result.json", canonical_bytes(result) + b"\n",
        )
        if (
            state.result_receipt is not None
            and state.result_receipt != result_receipt
        ):
            fail("parent containment result receipt changed")
        state.result_receipt = result_receipt
        summary = {
            "intent": intent_receipt, "path": state.path,
            "reason": state.reason, "result": result_receipt,
            "initial_roster": list(state.roster),
            "observed_roster": list(state.observed_roster),
        }
        if state.index == len(self.containment_receipts):
            self.containment_receipts.append(summary)
        elif self.containment_receipts[state.index] != summary:
            fail("parent containment receipt roster changed")
        for pid_value in state.observed_roster:
            if pid_value not in self.containment_roster:
                self.containment_roster.append(pid_value)
        return list(state.observed_roster)

    def resume_containments(self) -> None:
        if self.layout is None:
            if self.containment_states:
                fail("parent containment state outlived its layout")
            return
        authorities = {
            str(self.layout.experiment): self.layout.experiment_kill_fd,
            str(self.layout.verifier): self.layout.verifier_kill_fd,
            str(self.layout.sampler): self.layout.sampler_kill_fd,
            str(self.layout.run): self.layout.run_kill_fd,
        }
        for state in list(self.containment_states):
            # The authoritative aggregate is derived from durable per-event
            # state, not from a caller's post-return `list.extend` handoff.
            for pid_value in state.observed_roster:
                if pid_value not in self.containment_roster:
                    self.containment_roster.append(pid_value)
            if (
                state.result_receipt is not None
                and state.index < len(self.containment_receipts)
            ):
                continue
            if state.path not in authorities:
                fail("pending containment path authority differs")
            self.contain_cgroup(
                Path(state.path), authorities[state.path], state.reason,
            )

    def _pre_release_bound(self) -> None:
        if time.monotonic_ns() >= (
            self.process_started_ns + int(PRE_RELEASE_SECONDS * 1e9)
        ):
            fail("root pre-release service budget expired")

    def _post_controller_bound(self) -> None:
        if (
            self.controller_reap_ns is not None
            and time.monotonic_ns() >= self.controller_reap_ns
            + int(POST_CONTROLLER_SECONDS * 1e9)
        ):
            fail("root post-controller publication budget expired")

    def prepare(self) -> None:
        require_runtime_primitives()
        become_subreaper()
        self.prepared = prepare_execution(self.config)
        self.lock = CampaignLock.acquire()
        self.lock.verify()
        owned: List[int] = []
        open_registered(
            owned, "/var/tmp",
            os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC
            | getattr(os, "O_NOFOLLOW", 0),
        )
        self.var_tmp_fd = owned[0]
        require_execute_namespaces_absent(self.var_tmp_fd)
        relative = self.prepared.root_boundary["cgroup"]
        service = CGROUP_ROOT / relative.lstrip("/")
        if not service.is_dir():
            fail("delegated service cgroup path differs")
        self.service_path = service
        self._pre_release_bound()
        self.phase("read_only_admission")

    def reserve(self) -> None:
        if self.prepared is None or self.service_path is None:
            fail("attempt reservation precedes runtime admission")
        attempt = self_hashed(ATTEMPT_SCHEMA, {
            "attempt_id": ATTEMPT_ID,
            "build_authority_sha256": sha256_bytes(
                self.prepared.build_authority_raw),
            "campaign": CAMPAIGN,
            "expected_current_implementation_commit":
                CURRENT_IMPLEMENTATION_COMMIT,
            "expected_harness_commit": self.config.expected_harness_commit,
            "expected_parent_implementation_commit":
                PARENT_IMPLEMENTATION_COMMIT,
            "external_deadline_seconds": EXTERNAL_DEADLINE_SECONDS,
            "internal_deadline_seconds": INTERNAL_DEADLINE_SECONDS,
            "process_start_monotonic_ns": self.process_started_ns,
            "reserved_monotonic_ns": time.monotonic_ns(),
            "root_boundary": self.prepared.root_boundary,
            "service_deadline_seconds": SERVICE_DEADLINE_SECONDS,
        })
        try:
            AttemptJournal.reserve(
                ATTEMPT_DIR, attempt, owner=self._journal_owner,
                required_prefix=("ATTEMPT", "build-authority.json"),
            )
        except AttemptConsumedError as exc:
            if self.journal is None:
                self._journal_owner.adopt(exc.journal)
            raise
        if self.journal is None:
            fail("attempt journal adoption failed")
        self.journal.write_bytes_resumable(
            "build-authority.json", self.prepared.build_authority_raw,
        )
        CgroupTree(self.service_path).create(owner=self._layout_owner)
        if self.layout is None:
            fail("cgroup layout adoption failed")
        self.controller_directory = BoundDirectory(
            CONTROLLER_PARENT, self.var_tmp_fd, CAMPAIGN_UID, CAMPAIGN_GID,
        )
        self.sampler_directory = BoundDirectory(
            SAMPLER_DIR, self.var_tmp_fd, CAMPAIGN_UID, CAMPAIGN_GID,
        )
        self.controller_directory.create()
        self.sampler_directory.create()
        self.phase("attempt_reserved")

    def start_and_admit_sampler(self) -> SamplerPoller:
        assert self.prepared is not None and self.layout is not None
        target = sampler_argv(self.prepared)
        spawn_blocked_role(
            "sampler", self.layout.sampler, (SAMPLER_CPU,), target,
            sampler=True,
            python_info=os.stat(str(self.prepared.python_image),
                                follow_symlinks=False),
            owner=self._sampler_blocked_owner,
        )
        if self.sampler_blocked is None:
            fail("sampler process adoption failed")
        def mark_sampler_may_run() -> None:
            self.sampler_may_have_started = True

        process, security = self.sampler_blocked.release(
            release_callback=mark_sampler_may_run,
        )
        self.sampler_process = process
        wait_sampler_ready(
            process, self.prepared,
            os.stat(str(self.prepared.python_image), follow_symlinks=False),
            owner=self._sampler_owner,
        )
        if self.sampler is None:
            fail("sampler authority adoption failed")
        if self.sampler_blocked.final_security != security:
            fail("sampler final security handoff differs")
        self.publish_sampler_admission()
        self.phase("sampler_ready")
        return SamplerPoller(self.sampler)

    def publish_sampler_admission(self) -> None:
        """Resumably bind bootstrap, target, files, and root journal."""
        if self.sampler is None or self.sampler_blocked is None:
            fail("sampler admission lacks owned process authority")
        assert self.journal is not None
        target_security = self.sampler.admission.get("process_security")
        if type(target_security) is not dict:
            fail("sampler admission target security differs")
        security = self.sampler_blocked.final_security
        if security is None:
            security = dict(self.sampler_blocked.bootstrap_security)
            security["final_process"] = target_security
            self.sampler_blocked.final_security = security
        if (
            security.get("final_process") != target_security
            or target_security.get("pid") != self.sampler.process.pid
            or target_security.get("process_start_ticks")
            != self.sampler.process.start_ticks
        ):
            fail("sampler admission process-security cross-binding differs")
        self.sampler.admission["spawn_security"] = security
        payload = canonical_bytes(self_hashed(
            "wirehair.wh2.v2-facade-timing-root-sampler-admission.v1", {
                "admission": self.sampler.admission,
                "argv": list(self.sampler.process.argv),
                "pid": self.sampler.process.pid,
                "process_start_ticks": self.sampler.process.start_ticks,
            },
        )) + b"\n"
        receipt = self.journal.write_bytes_resumable(
            "sampler-admission.json", payload,
        )
        if (
            self.sampler_admission_receipt is not None
            and self.sampler_admission_receipt != receipt
        ):
            fail("sampler admission root receipt changed")
        self.sampler_admission_receipt = receipt
        self.sampler_admission_published = True

    def publish_sampler_cancellation(
            self, cancelled: Mapping[str, Any]) -> None:
        """Resumably bind a sampler stub proven not to have crossed its gate."""
        if self.sampler_blocked is None or self.journal is None:
            fail("sampler cancellation lacks owned authority")
        if (
            self.sampler_release_possible()
            or self.sampler_blocked.cancellation_result is None
            or dict(cancelled) != self.sampler_blocked.cancellation_result
        ):
            fail("sampler cancellation authority differs")
        self.journal.write_bytes_resumable(
            "sampler-cancelled.json",
            canonical_bytes(dict(cancelled)) + b"\n",
        )
        self.sampler_cancellation_published = True

    def spawn_blocked_roles(self, poller: SamplerPoller) -> None:
        assert self.prepared is not None and self.layout is not None
        assert self.controller_directory is not None
        controller_command = controller_run_argv(
            self.prepared, self.config, self.sampler, self.controller_directory,
        )
        verifier_command = replay_argv(
            self.prepared, self.config, CONTROLLER_OUTPUT / "raw.jsonl",
        )
        python_info = os.stat(
            str(self.prepared.python_image), follow_symlinks=False,
        )
        spawn_blocked_role(
            "controller", self.layout.experiment, RUN_CPUS,
            controller_command, sampler=False, python_info=python_info,
            owner=self._controller_blocked_owner,
            poll_callback=poller.poll,
        )
        if self.controller_blocked is None:
            fail("controller process adoption failed")
        try:
            spawn_blocked_role(
                "verifier", self.layout.verifier, (AUTHORITY_CPU,),
                verifier_command, sampler=False, python_info=python_info,
                owner=self._verifier_blocked_owner,
                poll_callback=poller.poll,
            )
            if self.verifier_blocked is None:
                fail("verifier process adoption failed")
        except BaseException:
            cancelled = drain_and_reap_cancelled(
                self.controller_blocked, self.layout.experiment,
                self.layout.experiment_kill_fd,
            )
            self.publish_controller_cancellation(cancelled)
            raise
        self.phase("roles_blocked")

    def arm_and_release(self, poller: SamplerPoller) -> None:
        assert self.journal is not None and self.layout is not None
        assert self.controller_blocked is not None
        self._pre_release_bound()
        # Reserve durable guardian evidence before declaring the external T0.
        # If reservation itself fails, no deadline phase exists.  Once T0 is
        # visible, every failure path has a journal-owned stream that must be
        # terminalized and authenticated before COMPLETE can be licensed.
        self.journal.open_stream("guardian.jsonl")
        self.guardian_event_fd = self.guardian_stream_fd()
        root_t0_ns = time.monotonic_ns()
        self._root_window = (
            root_t0_ns,
            root_t0_ns + EXTERNAL_DEADLINE_SECONDS * 1_000_000_000,
        )
        spawn_guardian(
            self.layout, self.root_deadline_ns, self.guardian_stream_fd(),
            self.journal.directory_fd,
            owner=self._guardian_owner,
        )
        if self.guardian is None:
            fail("independent guardian adoption failed")
        self.journal.write_json("start.json", self_hashed(
            "wirehair.wh2.v2-facade-timing-root-start.v1", {
                "controller_pid": self.controller_blocked.process.pid,
                "controller_process_start_ticks":
                    self.controller_blocked.process.start_ticks,
                "guardian_pid": self.guardian.process.pid,
                "guardian_process_start_ticks": self.guardian.process.start_ticks,
                "root_deadline_monotonic_ns": self.root_deadline_ns,
                "root_t0_monotonic_ns": self.root_t0_ns,
                "sampler_pid": self.sampler.process.pid,
                "sampler_process_start_ticks": self.sampler.process.start_ticks,
                "verifier_pid": self.verifier_blocked.process.pid,
                "verifier_process_start_ticks":
                    self.verifier_blocked.process.start_ticks,
            },
        ))
        self._pre_release_bound()
        def release_poll() -> None:
            poller.poll()
            assert self.guardian is not None
            self.guardian.poll()

        def mark_controller_may_run() -> None:
            self.controller_released = True

        # BlockedRole marks its owned irreversible state and invokes this
        # callback before the first possible gate write.
        self.controller, self.controller_security = self.controller_blocked.release(
            release_poll, mark_controller_may_run,
        )
        self.monitor = WorkerMapMonitor(
            self.prepared, self.layout, poller, self.guardian,
            self.controller, tuple(self.controller.argv),
            self.controller_security, self.root_deadline_ns,
        )
        self.monitor.poll()
        self.phase("controller_released")

    def run_controller(self) -> Tuple[bytes, bytes]:
        assert self.controller is not None and self.monitor is not None
        assert self.root_deadline_ns is not None
        communicate_process(
            self.controller, self.root_deadline_ns / 1e9,
            4096, MAX_STDERR_BYTES, self.monitor.poll,
        )
        return self.publish_controller_communication()

    def publish_controller_communication(self) -> Tuple[bytes, bytes]:
        """Resumably copy the exact retained controller transcript."""
        if self.controller is None or self.journal is None:
            fail("controller communication lacks owned authority")
        result = self.controller.communication_result
        if result is None:
            fail("controller communication transcript is incomplete")
        stdout, stderr, code = result
        if self.controller_rc is not None and self.controller_rc != code:
            fail("controller communication returncode changed")
        self.controller_rc = code
        if self.controller_reap_ns is None:
            self.controller_reap_ns = time.monotonic_ns()
        self.journal.write_bytes_resumable("controller.stdout", stdout)
        self.journal.write_bytes_resumable("controller.stderr", stderr)
        self.controller_communication_published = True
        if "controller_reaped" not in self.phases:
            self.phase("controller_reaped")
        return stdout, stderr

    def publish_controller_cancellation(
            self, cancelled: Mapping[str, Any]) -> None:
        """Resumably bind a controller stub that never crossed its gate."""
        if self.controller_blocked is None or self.journal is None:
            fail("controller cancellation lacks owned authority")
        if (
            self.controller_released or self.controller_blocked.may_have_released
            or self.controller_blocked.cancellation_result is None
            or dict(cancelled) != self.controller_blocked.cancellation_result
        ):
            fail("controller cancellation authority differs")
        code = cancelled.get("returncode")
        if type(code) is not int:
            fail("controller cancellation returncode differs")
        if self.controller_rc is not None and self.controller_rc != code:
            fail("controller cancellation returncode changed")
        self.controller = self.controller_blocked.process
        self.controller_rc = code
        if self.controller_reap_ns is None:
            self.controller_reap_ns = time.monotonic_ns()
        self.journal.write_bytes_resumable(
            "controller-cancelled.json",
            canonical_bytes(dict(cancelled)) + b"\n",
        )
        self.controller_cancellation_published = True

    def collect_child(self, controller_stdout: bytes,
                      controller_stderr: bytes) -> None:
        assert self.layout is not None and self.controller_directory is not None
        assert self.journal is not None
        roster = cgroup_descendant_pids(self.layout.experiment)
        if roster:
            self.containment_roster.extend(self.contain_cgroup(
                self.layout.experiment, self.layout.experiment_kill_fd,
                "experiment-survivors-after-controller",
            ))
            fail("experiment descendants survived controller reap")
        try:
            if self.bundle is None:
                open_child_bundle(
                    self.controller_directory, owner=self._bundle_owner,
                )
        except (FileNotFoundError, LaunchError):
            # An exception after the helper's persistent adoption belongs to
            # the admitted bundle, not to the no-output controller-rc1 case.
            if self.bundle is not None:
                raise
            if self.controller_rc == 1:
                if sorted(os.listdir(self.controller_directory.fd)):
                    raise
                self.child_outcome = "invalid"
                return
            raise
        if self.bundle is None:
            fail("controller child bundle ownership handoff failed")
        self.child_outcome = self.bundle.complete["outcome"]
        validate_child_bundle_binding(self.controller_directory, self.bundle)
        self.child_bundle_validated = True
        self.copy_child_evidence()
        if self.child_outcome in ("pass", "reject"):
            expected_rc = 0 if self.child_outcome == "pass" else 2
            if (
                self.controller_rc != expected_rc or controller_stdout
                or controller_stderr
                or self.bundle.data["current.stderr"]
                or self.bundle.data["parent.stderr"]
            ):
                fail("valid controller exit/stderr authority differs")
            assert self.monitor is not None
            self.map_receipt = self.monitor.finish()
            validate_child_authority(
                self.bundle, self.prepared, self.config,
                self.controller, self.controller_security,
                self.map_receipt, self.root_t0_ns, self.root_deadline_ns,
                self.controller_reap_ns,
            )
        elif self.controller_rc != 1:
            fail("invalid controller exit binding differs")
        self.phase("child_copied")

    def copy_child_evidence(self) -> None:
        if not self.child_bundle_validated or self.bundle is None:
            fail("root child-copy lacks an admitted held bundle")
        assert self.journal is not None
        copy_child_bundle(
            self.bundle, self.journal, self.child_copy_receipts,
        )
        if set(self.child_copy_receipts) != set(CONTROLLER_OUTPUT_NAMES):
            fail("root child-copy completion roster differs")
        self.child_copies_complete = True

    def run_or_cancel_verifier(self) -> None:
        self._post_controller_bound()
        assert self.verifier_blocked is not None and self.layout is not None
        assert self.journal is not None
        if self.child_outcome not in ("pass", "reject"):
            cancelled = drain_and_reap_cancelled(
                self.verifier_blocked, self.layout.verifier,
                self.layout.verifier_kill_fd,
                lambda: self.contain_cgroup(
                    self.layout.verifier, self.layout.verifier_kill_fd,
                    "cancel-invalid-verifier-stub",
                ),
            )
            self.verifier = self.verifier_blocked.process
            self.verifier_reap_ns = time.monotonic_ns()
            self.journal.write_bytes_resumable(
                "verifier-cancelled.json",
                canonical_bytes(cancelled) + b"\n",
            )
            self.verifier_resolved = True
            return
        assert self.guardian is not None and self.sampler is not None
        poller = SamplerPoller(self.sampler)
        def verifier_release_poll() -> None:
            poller.poll()
            self.guardian.poll()
            if cgroup_descendant_pids(self.layout.experiment):
                fail("experiment repopulated before retained replay release")

        self.verifier, self.verifier_security = self.verifier_blocked.release(
            verifier_release_poll,
        )
        deadline = min(
            time.monotonic() + REPLAY_SECONDS,
            (self.root_deadline_ns or 0) / 1e9,
            ((self.controller_reap_ns or 0)
             + int(POST_CONTROLLER_SECONDS * 1e9)) / 1e9,
        )
        communicate_process(
            self.verifier, deadline, MAX_DOCUMENT_BYTES,
            MAX_STDERR_BYTES,
            lambda: (poller.poll(), self.guardian.poll()),
        )
        self.publish_verifier_communication()

    def publish_verifier_communication(self) -> None:
        """Resumably validate and copy one released retained verifier."""
        if (
            self.verifier is None or self.journal is None
            or self.layout is None or self.bundle is None
            or self.child_outcome not in ("pass", "reject")
        ):
            fail("released verifier communication lacks owned authority")
        result = self.verifier.communication_result
        if result is None:
            fail("released verifier communication transcript is incomplete")
        stdout, stderr, code = result
        if self.verifier_rc is not None and self.verifier_rc != code:
            fail("released verifier communication returncode changed")
        self.verifier_rc = code
        if self.verifier_reap_ns is None:
            self.verifier_reap_ns = time.monotonic_ns()
        validation_error: Optional[BaseException] = None
        try:
            replay = validate_replay_result(stdout, stderr, code, self.bundle)
            if self.replay is not None and self.replay != replay:
                fail("retained replay validation result changed")
            self.replay = replay
        except BaseException as exc:
            validation_error = exc
        self.journal.write_bytes_resumable("replay.json", stdout)
        self.journal.write_bytes_resumable("replay.stderr", stderr)
        roster = cgroup_descendant_pids(self.layout.verifier)
        if roster:
            self.containment_roster.extend(self.contain_cgroup(
                self.layout.verifier, self.layout.verifier_kill_fd,
                "released-verifier-survivors",
            ))
            if validation_error is None:
                validation_error = LaunchError(
                    "retained replay cgroup remained populated"
                )
        self.verifier_resolved = True
        if "retained_replay" not in self.phases:
            self.phase("retained_replay")
        self._post_controller_bound()
        if validation_error is not None:
            raise validation_error

    def stop_sampler(self) -> None:
        self._post_controller_bound()
        assert self.sampler is not None and self.prepared is not None
        assert self.sampler_directory is not None and self.journal is not None
        self.sampler_terminal = stop_and_validate_sampler(
            self.sampler, self.prepared, self.sampler_directory, self.journal,
            required_start_before_ns=self.root_t0_ns,
            required_finish_after_ns=self.verifier_reap_ns,
            owner=self._sampler_terminal_owner,
        )
        # The strict producer receipt/stream/binding proof is complete here.
        # Preserve that independent fact if the subsequent child-health
        # crosscheck rejects untrusted controller evidence.
        self.sampler_authenticated = True
        if self.child_outcome in ("pass", "reject"):
            if self.bundle is None:
                fail("valid child bundle was lost before sampler crosscheck")
            validate_child_health_against_root(
                self.bundle, self.sampler, self.prepared,
                self.sampler_terminal,
            )
        self._post_controller_bound()
        self.phase("sampler_terminal")

    def finish_guardian(self) -> None:
        assert self.guardian is not None and self.journal is not None
        assert self.layout is not None
        if cgroup_descendant_pids(self.layout.run):
            fail("run cgroup is populated before guardian completion")
        self.guardian_receipt = self.guardian.complete()
        raw, receipt = self.journal.seal_stream(
            "guardian.jsonl", self.guardian_stream_fd(), 64 * 1024,
        )
        if sha256_bytes(raw) != receipt["sha256"]:
            fail("guardian root-journal copy differs")
        self.guardian_authenticated = True
        self.phase("guardian_terminal")

    def postflight(self) -> None:
        self._post_controller_bound()
        assert self.prepared is not None
        current_relative = process_cgroup_relative(os.getpid())
        allowed = {self.prepared.root_boundary["cgroup"]}
        if self.layout is not None:
            allowed.add(
                str(self.layout.supervisor)[len(str(CGROUP_ROOT)):],
            )
        if current_relative not in allowed:
            fail("root postflight cgroup phase differs")
        after = prepare_execution(
            self.config, expected_cgroup=current_relative,
        )
        before_values = dict(self.prepared.__dict__)
        after_values = dict(after.__dict__)
        before_boundary = dict(before_values.pop("root_boundary"))
        after_boundary = dict(after_values.pop("root_boundary"))
        before_boundary.pop("cgroup", None)
        after_boundary.pop("cgroup", None)
        if after_values != before_values or after_boundary != before_boundary:
            fail("postflight build/source/runtime authority changed")
        self.postflight_ok = True
        self._post_controller_bound()
        self.phase("postflight")

    def seal_namespaces(self) -> None:
        if self.var_tmp_fd < 0:
            fail("execution namespace parent authority is unavailable")

        def require_absent(path: Path, where: str) -> None:
            try:
                os.stat(path.name, dir_fd=self.var_tmp_fd,
                        follow_symlinks=False)
            except FileNotFoundError:
                return
            fail(where + " fixed namespace appeared after absent admission")

        controller = self.controller_directory
        if controller is None or controller.state == "absent":
            if (
                self.bundle is not None or self.child_bundle_validated
            ):
                fail("controller evidence exists without its parent authority")
            require_absent(CONTROLLER_PARENT, "controller parent")
        elif controller.state in ("bound", "sealed"):
            if self.bundle is not None:
                bundle = self.bundle
                if bundle.state == "open":
                    # Recover an async helper-return handoff before the normal
                    # collect_child admission assignments.  This licenses
                    # only descriptor-bound copying for an invalid root
                    # terminal; scientific pass/reject still needs the full
                    # controller/root authority path completed above.
                    if not self.child_bundle_validated:
                        validate_child_bundle_binding(controller, bundle)
                        self.child_outcome = bundle.complete["outcome"]
                        self.child_bundle_validated = True
                    # Resume any interrupted held-FD copy before surrendering
                    # the only child descriptors.
                    self.copy_child_evidence()
                    validate_child_bundle_binding(controller, bundle)
                elif bundle.state not in ("closing", "closed"):
                    fail("controller child bundle state differs before seal")
                if not self.child_copies_complete:
                    fail("controller child descriptors close before copy completion")
                bundle.close()
                if (
                    bundle.state != "closed" or bundle.directory_fd >= 0
                    or bundle.files or bundle.directory_identity is not None
                    or bundle.file_identities
                ):
                    fail("controller child descriptor closure was incomplete")
                # The closed object contains no authority.  Clearing its slot
                # is safe even if the caller is interrupted immediately after
                # this assignment; the namespace evidence is already copied.
                self._bundle_owner.value = None
            controller_names = tuple(sorted(os.listdir(controller.fd)))
            expected_controller_names = (
                (CONTROLLER_OUTPUT.name,) if self.child_bundle_validated else ()
            )
            if controller_names != expected_controller_names:
                fail("controller parent cannot be sealed with an inexact roster")
            controller.seal(controller_names)
        else:
            fail("controller namespace creation is unresolved")

        sampler = self.sampler_directory
        if sampler is None or sampler.state == "absent":
            if self.sampler_may_have_started:
                fail("sampler may have started without an evidence namespace")
            require_absent(SAMPLER_DIR, "sampler")
        elif sampler.state == "bound":
            sampler_names = tuple(sorted(os.listdir(sampler.fd)))
            if sampler_names:
                fail("unterminated sampler directory is not sealable")
            sampler.seal(())
        elif sampler.state == "sealed":
            expected_sampler_names = (
                tuple(sorted((
                    SAMPLER_CSV.name, SAMPLER_VALIDATION.name,
                    SAMPLER_RECEIPT.name,
                ))) if self.sampler_authenticated else ()
            )
            sampler.seal(expected_sampler_names)
        else:
            fail("sampler namespace creation is unresolved")
        self.namespaces_sealed = True

    def _wait_tracked_after_kill(self) -> None:
        unique: Dict[int, ProcessHandle] = {}
        for process in (
            self.controller,
            self.controller_blocked.process if self.controller_blocked else None,
            self.verifier,
            self.verifier_blocked.process if self.verifier_blocked else None,
            self.sampler_process,
            self.sampler_blocked.process if self.sampler_blocked else None,
            self.sampler.process if self.sampler else None,
        ):
            if process is not None:
                unique[process.pid] = process
        for process in unique.values():
            if process.returncode is None:
                wait_process(process, time.monotonic() + 2.0)

    def resolve_guardian(self) -> None:
        if self.guardian_authenticated:
            return
        if self.guardian is None:
            if self.controller_released:
                fail("released controller lacks an independent root guardian")
            return
        assert self.journal is not None and self.layout is not None
        if cgroup_descendant_pids(self.layout.run):
            fail("run cgroup remains populated before guardian resolution")
        if self.guardian.process.returncode is None and not process_exited(
            self.guardian.process
        ):
            self.guardian_receipt = self.guardian.complete()
        else:
            code = wait_process(
                self.guardian.process, time.monotonic() + 2.0,
            )
            events = validate_guardian_events(
                self.guardian.event_fd, self.guardian.process.pid,
                self.guardian.ready["deadline_monotonic_ns"], None, code,
            )
            self.guardian_receipt = {
                **self.guardian.ready, "completed": False,
                "events": events, "returncode": code,
            }
            self.record(
                "root guardian",
                LaunchError("terminal event {} exit {}".format(
                    events[-1]["event"], code)),
            )
        raw, receipt = self.journal.seal_stream(
            "guardian.jsonl", self.guardian_stream_fd(), 64 * 1024,
        )
        if sha256_bytes(raw) != receipt["sha256"]:
            fail("guardian root-journal seal differs")
        self.guardian_authenticated = True

    def seal_failed_guardian_admission(self) -> None:
        if self.guardian is not None or self.journal is None:
            return
        event_fd = self.guardian_stream_fd(required=False)
        if event_fd < 0 or self.root_deadline_ns is None:
            return
        raw = read_held_file(
            event_fd, 64 * 1024,
            "failed guardian admission journal",
            allow_empty=True,
        )
        if not raw:
            value = {
                "deadline_monotonic_ns": self.root_deadline_ns,
                "error": "guardian spawn did not establish child authority",
                "event": "failure", "observed_monotonic_ns": time.monotonic_ns(),
                "pid": os.getpid(),
                "schema":
                    "wirehair.wh2.v2-facade-timing-root-guardian-event.v1",
            }
            append_guardian_record(event_fd, self.journal.directory_fd, value)
            raw = read_held_file(
                event_fd, 64 * 1024,
                "failed guardian admission journal",
            )
        records = raw.splitlines(keepends=True)
        values = [
            parse_canonical_document(
                record, "failed guardian admission event", 4096,
            )
            for record in records
        ]
        schema = "wirehair.wh2.v2-facade-timing-root-guardian-event.v1"
        if len(values) == 1 and values[0].get("event") == "armed":
            armed = values[0]
            if (
                set(armed) != {
                    "armed_monotonic_ns", "deadline_monotonic_ns",
                    "event", "pid", "schema",
                }
                or armed.get("schema") != schema
                or armed.get("deadline_monotonic_ns") != self.root_deadline_ns
                or type(armed.get("pid")) is not int or armed["pid"] <= 0
                or type(armed.get("armed_monotonic_ns")) is not int
                or not 0 <= armed["armed_monotonic_ns"] < self.root_deadline_ns
            ):
                fail("failed guardian durable arm differs")
            append_guardian_record(event_fd, self.journal.directory_fd, {
                "deadline_monotonic_ns": self.root_deadline_ns,
                "error": "parent admission failed after durable guardian arm",
                "event": "failure",
                "observed_monotonic_ns": max(
                    time.monotonic_ns(), armed["armed_monotonic_ns"],
                ),
                "pid": armed["pid"], "schema": schema,
            })
            raw = read_held_file(
                event_fd, 64 * 1024,
                "failed guardian admission journal",
            )
            records = raw.splitlines(keepends=True)
            values = [
                parse_canonical_document(
                    record, "failed guardian admission event", 4096,
                )
                for record in records
            ]
        if len(values) == 2:
            events = validate_guardian_events(
                event_fd, values[0].get("pid", -1), self.root_deadline_ns,
                None, None,
            )
            if events[-1]["event"] == "completed":
                fail("unadmitted guardian reported successful completion")
        elif len(values) == 1:
            value = values[0]
            if (
                set(value) != {
                    "deadline_monotonic_ns", "error", "event",
                    "observed_monotonic_ns", "pid", "schema",
                }
                or value.get("schema") != schema
                or value.get("event") != "failure"
                or value.get("deadline_monotonic_ns") != self.root_deadline_ns
                or type(value.get("pid")) is not int or value["pid"] <= 0
                or type(value.get("observed_monotonic_ns")) is not int
                or value["observed_monotonic_ns"] < 0
                or type(value.get("error")) is not str or not value["error"]
                or len(value["error"]) > 4096
            ):
                fail("failed guardian admission terminal differs")
        else:
            fail("failed guardian admission event count differs")
        sealed, receipt = self.journal.seal_stream(
            "guardian.jsonl", event_fd, 64 * 1024,
        )
        if sealed != raw:
            fail("failed guardian admission stream changed during seal")
        self.guardian_receipt = {
            "admission_failed": True, "events": values,
            "stream_sha256": receipt["sha256"],
        }
        self.guardian_authenticated = True

    def quiesce(self) -> None:
        if self.layout is None:
            self.quiesced = True
            return
        if not self.layout.run.exists():
            if any(path.exists() for path in (
                self.layout.sampler, self.layout.experiment,
                self.layout.verifier,
            )):
                fail("partial cgroup children exist without their run parent")
            self.quiesced = True
            return
        roster = cgroup_descendant_pids(self.layout.run)
        if roster:
            self.unexpected_final_containment = True
            if self.layout.run_kill_fd < 0:
                fail("populated partial run cgroup lacks kill authority")
            self.containment_roster.extend(self.contain_cgroup(
                self.layout.run, self.layout.run_kill_fd,
                "final-run-quiescence",
            ))
            self.record(
                "final run containment",
                LaunchError("unexpected run descendants required cgroup.kill"),
            )
        self._wait_tracked_after_kill()
        if cgroup_descendant_pids(self.layout.run):
            fail("run cgroup repopulated after exact containment")
        self.quiesced = True

    def close_runtime_authority(self) -> None:
        failures: List[str] = []
        unique: Dict[int, ProcessHandle] = {}
        for process in (
            self.controller,
            self.controller_blocked.process if self.controller_blocked else None,
            self.verifier,
            self.verifier_blocked.process if self.verifier_blocked else None,
            self.sampler_process,
            self.sampler_blocked.process if self.sampler_blocked else None,
            self.sampler.process if self.sampler else None,
            self.guardian.process if self.guardian else None,
        ):
            if process is not None:
                unique[process.pid] = process
        for process in unique.values():
            try:
                strict_process_close(process)
            except BaseException as exc:
                failures.append(exception_text(exc))
        if self.sampler is not None:
            try:
                self.sampler.close()
            except BaseException as exc:
                failures.append("sampler held FD: " + exception_text(exc))
        for blocked in (
            self.sampler_blocked, self.controller_blocked,
            self.verifier_blocked,
        ):
            if blocked is not None and blocked.gate_fd >= 0:
                try:
                    blocked.cancel()
                except BaseException as exc:
                    failures.append("blocked gate FD: " + exception_text(exc))
        if self.guardian is not None:
            if self.guardian.control_fd >= 0:
                try:
                    self.guardian.close()
                except BaseException as exc:
                    failures.append("guardian control FD: " + exception_text(exc))
        event_fd = self.guardian_stream_fd(required=False)
        if event_fd >= 0:
            try:
                # Clear both aliases before the irreversible close; restore
                # them only if the exact descriptor remains live.
                self.guardian_event_fd = -1
                if self.journal is not None:
                    registered = self.journal.stream_fds.get("guardian.jsonl")
                    if registered != event_fd:
                        fail("guardian stream registry changed before closure")
                    self.journal.stream_fds.pop("guardian.jsonl", None)
                try:
                    close_fd_once(event_fd, "guardian event stream")
                except BaseException:
                    self.guardian_event_fd = event_fd
                    if self.journal is not None:
                        self.journal.stream_fds["guardian.jsonl"] = event_fd
                    raise
            except BaseException as exc:
                failures.append("guardian event FD: " + exception_text(exc))
        for directory in (self.controller_directory, self.sampler_directory):
            if directory is not None and directory.fd >= 0:
                try:
                    directory.close()
                except BaseException as exc:
                    failures.append("owned directory FD: " + exception_text(exc))
        if self.layout is not None:
            try:
                self.layout.close()
            except BaseException as exc:
                failures.append("cgroup authority FD: " + exception_text(exc))
            if not failures:
                try:
                    for path in (
                        self.layout.sampler, self.layout.experiment,
                        self.layout.verifier, self.layout.run,
                    ):
                        if not path.exists():
                            continue
                        if cgroup_descendant_pids(path):
                            fail("cgroup populated during layout removal: " + str(path))
                        os.rmdir(str(path))
                    if self.layout.supervisor.exists():
                        CgroupTree.write(
                            self.layout.service / "cgroup.subtree_control",
                            "-cpuset -pids",
                        )
                        CgroupTree.write(
                            self.layout.service / "cgroup.procs", str(os.getpid()),
                        )
                        os.sched_setaffinity(0, {AUTHORITY_CPU})
                        if (
                            process_cgroup_relative(os.getpid())
                            != str(self.layout.service)[len(str(CGROUP_ROOT)):]
                            or os.sched_getaffinity(0) != {AUTHORITY_CPU}
                            or CgroupTree.pids(self.layout.supervisor)
                        ):
                            fail("root supervisor did not leave its temporary leaf")
                        os.rmdir(str(self.layout.supervisor))
                    child_directories = [
                        name for name in os.listdir(self.layout.service)
                        if (self.layout.service / name).is_dir()
                    ]
                    if child_directories:
                        fail("service cgroup child roster remained after removal")
                    self.layout_removed = True
                except BaseException as exc:
                    failures.append("cgroup layout removal: " + exception_text(exc))
        else:
            self.layout_removed = True
        if self.var_tmp_fd >= 0:
            try:
                close_object_fd(self, "var_tmp_fd", "/var/tmp authority")
            except BaseException as exc:
                failures.append("/var/tmp authority FD: " + exception_text(exc))
        if self.lock is not None:
            try:
                self.lock.close()
            except BaseException as exc:
                failures.append("campaign lock FD: " + exception_text(exc))
        try:
            reap_adopted_children(self.adopted_reaps)
        except BaseException as exc:
            failures.append("subreaper sweep: " + exception_text(exc))
        for item in self.adopted_reaps:
            if item["pid"] not in self.containment_roster:
                self.containment_roster.append(item["pid"])
        if self.adopted_reaps and not self.adopted_reaps_reported:
            # Publish the one-shot invalidating state before returning the
            # failure.  A retry can then close runtime authority without losing
            # the exact already-reaped PID/status evidence.
            self.adopted_reaps_reported = True
            failures.append("untracked adopted descendants were reaped")
        if self.journal is not None:
            try:
                verify_fd_authority_roster(
                    self.initial_fd_authority,
                    "root initial descriptor authority",
                )
                self.journal.verify_descriptor_authority()
                close_fd_delta(
                    self.initial_fds,
                    (self.journal.parent_fd, self.journal.directory_fd),
                    closed_sink=self.recovered_fds,
                )
                if self.recovered_fds:
                    message = (
                        "recovered untracked descriptor delta: "
                        + repr(self.recovered_fds)
                    )
                    if message not in self.errors:
                        self.errors.append(message)
            except BaseException as exc:
                failures.append("final FD-delta recovery: " + exception_text(exc))
        if failures:
            fail("runtime authority closure differs: " + " | ".join(failures[:16]))
        self.runtime_closed = True

    def completion_licensed(self) -> bool:
        guardian_phase_started = bool(
            self.root_t0_ns is not None
            or self.root_deadline_ns is not None
            or self.guardian is not None
            or (
                self.journal is not None
                and "guardian.jsonl" in self.journal.names
            )
        )
        return bool(
            self.journal is not None
            and self.quiesced and self.postflight_ok and self.namespaces_sealed
            and self.runtime_closed and self.layout_removed
            and all(
                state.result_receipt is not None
                for state in self.containment_states
            )
            and len(self.containment_receipts)
            == len(self.containment_states)
            and set(
                pid_value
                for state in self.containment_states
                for pid_value in state.observed_roster
            ).issubset(set(self.containment_roster))
            and (
                (
                    not self.sampler_release_possible()
                    and (
                        self.sampler_blocked is None
                        or self.sampler_cancellation_published
                    )
                )
                or (
                    self.sampler_release_possible()
                    and self.sampler_authenticated
                    and self.sampler_admission_published
                )
            )
            and (self.verifier_blocked is None or self.verifier_resolved)
            and (
                not self.child_bundle_validated
                or self.child_copies_complete
            )
            and (not guardian_phase_started or self.guardian_authenticated)
            and (
                (
                    not self.controller_released
                    and (
                        self.controller_blocked is None
                        or self.controller_cancellation_published
                    )
                )
                or (
                    self.controller_released
                    and self.guardian is not None
                    and self.guardian_authenticated
                    and self.controller_communication_published
                )
            )
        )

    def publish_terminal(self, outcome: str) -> RootRunResult:
        if outcome not in ("pass", "reject", "invalid"):
            fail("root terminal outcome differs")
        if outcome in ("pass", "reject") and (
            self.unexpected_final_containment or self.containment_states
        ):
            fail("scientific outcome cannot survive parent containment")
        if not self.completion_licensed() or self.journal is None:
            fail("attempt is not licensed for root COMPLETE")
        if self.prepared is None:
            fail("root terminal lacks prepared build authority")
        self.journal.write_bytes_resumable(
            "build-authority.json", self.prepared.build_authority_raw,
        )
        self._post_controller_bound()
        expected_fds = set(self.initial_fds) | {
            self.journal.parent_fd, self.journal.directory_fd,
        }
        verify_fd_authority_roster(
            self.initial_fd_authority,
            "root terminal initial descriptor authority",
        )
        self.journal.verify_descriptor_authority()
        if set(open_fd_roster()) != expected_fds:
            fail("root descriptor roster is not exact before terminal publication")
        terminal = self_hashed(TERMINAL_SCHEMA, {
            "attempt_id": ATTEMPT_ID,
            "adopted_reaps": list(self.adopted_reaps),
            "campaign": CAMPAIGN,
            "child_copies": dict(self.child_copy_receipts),
            "child_outcome": self.child_outcome,
            "containment_events": list(self.containment_receipts),
            "containment_roster": sorted(set(self.containment_roster)),
            "controller_returncode": self.controller_rc,
            "errors": list(self.errors),
            "guardian": self.guardian_receipt,
            "launcher_status": "clean" if not self.errors else "fault",
            "map_receipt": self.map_receipt,
            "outcome": outcome,
            "phases": list(self.phases),
            "postflight_authority_sha256": sha256_bytes(
                self.prepared.build_authority_raw),
            "replay": self.replay,
            "recovered_fds": list(self.recovered_fds),
            "root_deadline_monotonic_ns": self.root_deadline_ns,
            "root_t0_monotonic_ns": self.root_t0_ns,
            "sampler": self.sampler_terminal,
            "sampler_admission": self.sampler_admission_receipt,
            "verifier_returncode": self.verifier_rc,
            "unexpected_final_containment": self.unexpected_final_containment,
        })
        terminal_raw = canonical_bytes(terminal) + b"\n"
        receipt: Optional[Dict[str, Any]] = None
        publication_error: Optional[BaseException] = None
        for attempt in range(2):
            try:
                receipt = self.journal.write_bytes_resumable(
                    "terminal.json", terminal_raw,
                )
                break
            except BaseException as exc:
                publication_error = exc
                if attempt:
                    raise
        if receipt is None:
            assert publication_error is not None
            raise publication_error
        for attempt in range(2):
            try:
                if not self.journal.complete:
                    self.journal.finish()
                break
            except BaseException:
                if self.journal.complete:
                    break
                if attempt:
                    raise
        if not self.journal.complete:
            fail("root journal COMPLETE publication was not established")
        self.journal.close()
        exit_code = 0 if outcome == "pass" and not self.errors else (
            2 if outcome == "reject" and not self.errors else 1
        )
        return RootRunResult(
            outcome, "clean" if not self.errors else "fault",
            exit_code, receipt["sha256"],
        )

    def execute(self) -> RootRunResult:
        try:
            self.prepare()
            self.reserve()
            poller = self.start_and_admit_sampler()
            self.spawn_blocked_roles(poller)
            self.arm_and_release(poller)
            controller_stdout, controller_stderr = self.run_controller()
            self.collect_child(controller_stdout, controller_stderr)
            self.run_or_cancel_verifier()
            self.stop_sampler()
            self.seal_namespaces()
            self.quiesce()
            self.resolve_guardian()
            self.postflight()
            outcome = self.child_outcome or "invalid"
        except BaseException as exc:
            self.record("execution", exc)
            outcome = "invalid"
            if (
                self.guardian is None and self._root_window is None
                and self.journal is not None
                and "guardian.jsonl" in self.journal.names
            ):
                try:
                    # spawn_guardian is invoked only after the atomic window
                    # assignment.  Therefore this exact empty stream was
                    # never shared with a child and can be rolled back.
                    event_fd = self.guardian_stream_fd()
                    self.guardian_event_fd = -1
                    try:
                        self.journal.discard_empty_stream("guardian.jsonl")
                    except BaseException:
                        if self.journal.stream_fds.get(
                            "guardian.jsonl"
                        ) == event_fd:
                            self.guardian_event_fd = event_fd
                        raise
                except BaseException as cleanup_exc:
                    self.record("guardian pre-arm rollback", cleanup_exc)
            # If the guardian child was not admitted after T0, durably close
            # that failure record before any containment kill.  This keeps
            # the root journal authoritative even if the supervisor dies
            # during the subsequent cleanup sequence.
            if (
                self.guardian is None
                and self.root_deadline_ns is not None
                and not self.guardian_authenticated
            ):
                try:
                    self.seal_failed_guardian_admission()
                except BaseException as cleanup_exc:
                    self.record("guardian admission terminal", cleanup_exc)
            try:
                if self.verifier_blocked is not None and not self.verifier_resolved:
                    assert self.layout is not None
                    if self.verifier_blocked.may_have_released:
                        self.verifier = self.verifier_blocked.process
                        if self.verifier.communication_result is None:
                            communicate_process(
                                self.verifier, time.monotonic() + 2.0,
                                MAX_DOCUMENT_BYTES, MAX_STDERR_BYTES,
                            )
                        self.publish_verifier_communication()
                    else:
                        cancelled = drain_and_reap_cancelled(
                            self.verifier_blocked, self.layout.verifier,
                            self.layout.verifier_kill_fd,
                            lambda: self.contain_cgroup(
                                self.layout.verifier,
                                self.layout.verifier_kill_fd,
                                "cleanup-verifier-stub",
                            ),
                        )
                        if self.journal is None:
                            fail("cancelled verifier lacks attempt journal")
                        self.journal.write_bytes_resumable(
                            "verifier-cancelled.json",
                            canonical_bytes(cancelled) + b"\n",
                        )
                        self.verifier = self.verifier_blocked.process
                        if self.verifier_reap_ns is None:
                            self.verifier_reap_ns = time.monotonic_ns()
                        self.verifier_resolved = True
            except BaseException as cleanup_exc:
                self.record("verifier containment", cleanup_exc)
            try:
                controller_process = (
                    self.controller if self.controller is not None else
                    self.controller_blocked.process
                    if self.controller_blocked is not None else None
                )
                if (
                    self.layout is not None
                    and self.controller_blocked is not None
                    and not self.controller_released
                    and not self.controller_blocked.may_have_released
                ):
                    cancelled = drain_and_reap_cancelled(
                        self.controller_blocked, self.layout.experiment,
                        self.layout.experiment_kill_fd,
                        lambda: self.contain_cgroup(
                            self.layout.experiment,
                            self.layout.experiment_kill_fd,
                            "cleanup-controller-stub",
                        ),
                    )
                    self.publish_controller_cancellation(cancelled)
                else:
                    if self.layout is not None:
                        experiment = cgroup_descendant_pids(
                            self.layout.experiment)
                        if experiment:
                            self.containment_roster.extend(self.contain_cgroup(
                                self.layout.experiment,
                                self.layout.experiment_kill_fd,
                                "cleanup-experiment",
                            ))
                if controller_process is not None and self.controller_released:
                    if controller_process.communication_result is None:
                        communicate_process(
                            controller_process, time.monotonic() + 2.0,
                            4096, MAX_STDERR_BYTES,
                        )
                    if self.controller is None:
                        self.controller = controller_process
                    self.publish_controller_communication()
                elif (
                    controller_process is not None
                    and controller_process.returncode is None
                ):
                    wait_process(controller_process, time.monotonic() + 2.0)
            except BaseException as cleanup_exc:
                self.record("controller containment", cleanup_exc)
            try:
                # A still-held bootstrap gate proves that no sampler target
                # byte was released.  Cancel/reap that exact stub before
                # clearing the conservative may-run flag.  Once the gate has
                # closed, recover target authority instead; killing an
                # unauthenticated possibly-running sampler cannot license a
                # root COMPLETE.
                if (
                    self.sampler_blocked is not None
                    and self.sampler_blocked.gate_fd >= 0
                    and not self.sampler_blocked.may_have_released
                    and not self.sampler_may_have_started
                    and self.layout is not None
                ):
                    cancelled = drain_and_reap_cancelled(
                        self.sampler_blocked, self.layout.sampler,
                        self.layout.sampler_kill_fd,
                        lambda: self.contain_cgroup(
                            self.layout.sampler,
                            self.layout.sampler_kill_fd,
                            "cleanup-blocked-sampler",
                        ),
                    )
                    self.publish_sampler_cancellation(cancelled)
                    self.sampler_may_have_started = False
                elif (
                    self.sampler is None and self.sampler_release_possible()
                    and self.sampler_blocked is not None
                    and self.prepared is not None
                ):
                    process = self.sampler_blocked.process
                    self.sampler_process = process
                    wait_sampler_ready(
                        process, self.prepared,
                        os.stat(
                            str(self.prepared.python_image),
                            follow_symlinks=False,
                        ),
                        owner=self._sampler_owner,
                    )
                    if self.sampler is None:
                        fail("cleanup sampler authority adoption failed")
                if self.sampler is not None and self.sampler_release_possible():
                    self.publish_sampler_admission()
                experiment_safe = bool(
                    self.layout is not None
                    and not cgroup_descendant_pids(self.layout.experiment)
                    and (
                        self.controller_blocked is None
                        or self.controller_blocked.process.returncode is not None
                    )
                )
                verifier_safe = bool(
                    self.layout is not None
                    and (
                        self.verifier_blocked is None
                        or self.verifier_resolved
                    )
                    and not cgroup_descendant_pids(self.layout.verifier)
                )
                if (
                    self.sampler is not None and not self.sampler_authenticated
                    and experiment_safe and verifier_safe
                ):
                    self.stop_sampler()
            except BaseException as cleanup_exc:
                self.record("sampler terminal", cleanup_exc)
            try:
                self.resume_containments()
            except BaseException as cleanup_exc:
                self.record("containment receipt recovery", cleanup_exc)
            try:
                self.quiesce()
            except BaseException as cleanup_exc:
                self.record("run quiescence", cleanup_exc)
            try:
                # quiesce() may itself create the final run-containment
                # record.  Recover a result-publication handoff from that
                # call before resolving the guardian or surrendering cgroup
                # authority.
                self.resume_containments()
            except BaseException as cleanup_exc:
                self.record("post-quiescence containment recovery", cleanup_exc)
            try:
                if self.guardian is not None:
                    self.resolve_guardian()
                else:
                    self.seal_failed_guardian_admission()
            except BaseException as cleanup_exc:
                self.record("guardian terminal", cleanup_exc)
            try:
                if self.var_tmp_fd >= 0:
                    self.seal_namespaces()
            except BaseException as cleanup_exc:
                self.record("namespace seal", cleanup_exc)
            try:
                if self.prepared is not None:
                    self.postflight()
            except BaseException as cleanup_exc:
                self.record("postflight", cleanup_exc)
        for closure_attempt in range(2):
            try:
                self.close_runtime_authority()
                break
            except BaseException as exc:
                self.record("descriptor closure", exc)
                if closure_attempt:
                    break
        if self.journal is None:
            raise LaunchError(self.errors[0] if self.errors else "launch admission failed")
        if not self.completion_licensed():
            self.journal.close()
            raise LaunchError("attempt consumed but root COMPLETE is not licensed: "
                              + " | ".join(self.errors[:8]))
        if self.errors:
            outcome = "invalid"
        return self.publish_terminal(outcome)


def execute_once(config: LaunchConfig) -> int:
    result = ExecutionSession(config).execute()
    print(canonical_bytes({
        "launcher_status": result.launcher_status,
        "outcome": result.outcome,
        "schema": "wirehair.wh2.v2-facade-timing-launch-result.v1",
        "status": "complete",
        "terminal_sha256": result.terminal_sha256,
    }).decode("ascii"))
    return result.exit_code


def prepare_execution(config: LaunchConfig, *,
                      expected_cgroup: Optional[str] = None) -> PreparedExecution:
    root_boundary = require_root_mode(
        UNIT_NAME, (AUTHORITY_CPU,), expected_cgroup=expected_cgroup,
    )
    digest, authority_info = hash_path(
        BUILD_AUTHORITY_PATH, MAX_DOCUMENT_BYTES, "fixed build authority",
    )
    if (
        authority_info.st_uid != 0 or authority_info.st_gid != 0
        or stat.S_IMODE(authority_info.st_mode) != 0o400
    ):
        fail("fixed build authority policy differs")
    authority_raw = BUILD_AUTHORITY_PATH.read_bytes()
    if sha256_bytes(authority_raw) != digest:
        fail("fixed build authority readback differs")
    authority = parse_canonical_document(
        authority_raw, "fixed build authority", MAX_DOCUMENT_BYTES,
    )
    validate_self_hash(authority, BUILD_AUTHORITY_SCHEMA)
    if (
        authority.get("status") != "sealed"
        or authority.get("expected_harness_commit") != config.expected_harness_commit
        or authority.get("expected_current_implementation_commit")
        != CURRENT_IMPLEMENTATION_COMMIT
        or authority.get("expected_parent_implementation_commit")
        != PARENT_IMPLEMENTATION_COMMIT
        or authority.get("frozen_harness_source_sha256") != {
            CONTROLLER_RELATIVE: FROZEN_CONTROLLER_SHA256,
            WORKER_RELATIVE: FROZEN_WORKER_SOURCE_SHA256,
            ROLE_CMAKE_RELATIVE: FROZEN_ROLE_CMAKE_SHA256,
            HEALTH_ADAPTER_RELATIVE: FROZEN_HEALTH_ADAPTER_SHA256,
        }
        or authority.get("trust_boundary", {}).get(
            "concurrent_hostile_uid1000_processes_excluded") is not True
        or authority.get("trust_boundary", {}).get(
            "root_build_graph_is_sole_build_authority") is not True
    ):
        fail("fixed build authority contract differs")
    build_root = Path(authority["build_root"])
    exact_absolute(build_root, "sealed build root", must_exist=True)
    build_info = os.stat(str(build_root), follow_symlinks=False)
    if (
        not stat.S_ISDIR(build_info.st_mode) or build_info.st_uid != 0
        or build_info.st_gid != 0 or stat.S_IMODE(build_info.st_mode) != 0o555
    ):
        fail("sealed build root policy differs")
    installed = authority.get("installed_launcher")
    if type(installed) is not dict:
        fail("build authority installed launcher receipt differs")
    installed_live = content_receipt(INSTALLED_LAUNCHER, "execute installed launcher")
    installed_info = os.stat(str(INSTALLED_LAUNCHER), follow_symlinks=False)
    if (
        installed_live["sha256"] != installed["sha256"]
        or installed_live["bytes"] != installed["bytes"]
        or installed_info.st_uid != 0 or installed_info.st_gid != 0
        or installed_info.st_nlink != 1
        or stat.S_IMODE(installed_info.st_mode) != 0o555
    ):
        fail("installed launcher/build authority binding differs")
    # No external source/tool command is permitted until its exact sealed
    # image roster has been rebound to the fixed build authority.
    live_toolchain = toolchain_receipt()
    if live_toolchain != authority.get("toolchain"):
        fail("execute toolchain/build authority binding differs")
    validate_self_runtime_map_closure(authority)
    validate_controller_runtime_map_closure(authority)

    roles: Dict[str, RoleExecution] = {}
    receipts: Dict[str, Dict[str, Any]] = {}
    for role, commit in (
        ("current", CURRENT_IMPLEMENTATION_COMMIT),
        ("parent", PARENT_IMPLEMENTATION_COMMIT),
    ):
        metadata = authority.get("role_receipts", {}).get(role)
        if type(metadata) is not dict:
            fail(role + " build-authority receipt member differs")
        receipt_path = Path(metadata["path"])
        receipt_digest, receipt_info = hash_path(
            receipt_path, MAX_DOCUMENT_BYTES, role + " build receipt",
        )
        if (
            receipt_digest != metadata["sha256"]
            or receipt_info.st_uid != 0 or receipt_info.st_gid != 0
            or stat.S_IMODE(receipt_info.st_mode) != 0o444
        ):
            fail(role + " build receipt policy differs")
        receipt_raw = receipt_path.read_bytes()
        if sha256_bytes(receipt_raw) != receipt_digest:
            fail(role + " build receipt readback differs")
        receipt = parse_canonical_document(
            receipt_raw, role + " build receipt", MAX_DOCUMENT_BYTES,
        )
        validate_self_hash(
            receipt, BUILD_RECEIPT_SCHEMA, "receipt_preimage_sha256",
        )
        if (
            receipt.get("role") != role
            or receipt.get("harness_git_commit") != config.expected_harness_commit
            or receipt.get("implementation_git_commit") != commit
        ):
            fail(role + " build receipt identity differs")
        for roster_name in (
            "source_manifest", "object_closure", "archive_closure",
            "dynamic_dependencies",
        ):
            roster = receipt.get(roster_name)
            if type(roster) is not list:
                fail(role + " receipt roster differs: " + roster_name)
            for entry in roster:
                if type(entry) is not dict:
                    fail(role + " receipt artifact entry differs")
                path = Path(entry["path"])
                digest_now, info_now = hash_path(
                    path, max(entry["bytes"], 1),
                    role + " receipt artifact " + str(path),
                    expected_nlink=build_receipt_expected_nlink(
                        roster_name, str(path), receipt["role_library_path"],
                    ),
                )
                if digest_now != entry["sha256"] or info_now.st_size != entry["bytes"]:
                    fail(role + " receipt artifact changed: " + str(path))
        worker = Path(receipt["worker_path"])
        library = Path(receipt["role_library_path"])
        implementation_library = Path(receipt["implementation_library_path"])
        worker_digest, worker_info = hash_path(
            worker, 256 * 1024 * 1024, role + " worker",
        )
        if (
            worker_digest != receipt["worker_sha256"]
            or worker_info.st_uid != 0 or worker_info.st_gid != 0
            or stat.S_IMODE(worker_info.st_mode) != 0o555
        ):
            fail(role + " worker policy differs")
        for path, expected, where in (
            (library, receipt["role_library_sha256"], "role library"),
            (implementation_library,
             receipt["implementation_library_sha256"],
             "implementation library"),
        ):
            artifact_digest, artifact_info = hash_path(
                path, MAX_TRACKED_FILE_BYTES, role + " " + where,
            )
            if (
                artifact_digest != expected or artifact_info.st_uid != 0
                or artifact_info.st_gid != 0
                or stat.S_IMODE(artifact_info.st_mode) != 0o444
            ):
                fail(role + " " + where + " policy differs")
        runtime_entries = authority.get("elf", {}).get(role, {}).get(
            "runtime_map_dependencies")
        if type(runtime_entries) is not list or not runtime_entries:
            fail(role + " runtime map closure differs")
        runtime_receipts: List[Dict[str, Any]] = []
        for entry in runtime_entries:
            path = Path(entry["path"])
            digest_now, info_now = hash_path(
                path, max(entry["bytes"], 1), role + " runtime map input",
                expected_nlink=None,
            )
            if digest_now != entry["sha256"] or info_now.st_size != entry["bytes"]:
                fail(role + " runtime map input changed")
            runtime_receipts.append(dict(entry))
        roles[role] = RoleExecution(
            role=role, worker=worker, worker_sha256=receipt["worker_sha256"],
            library=library, library_sha256=receipt["role_library_sha256"],
            implementation_library=implementation_library,
            implementation_library_sha256=receipt[
                "implementation_library_sha256"],
            build_receipt=receipt_path, build_receipt_sha256=receipt_digest,
            worker_receipt={
                "bytes": worker_info.st_size, "path": str(worker),
                "sha256": worker_digest, "stat": stat_receipt(worker_info),
            },
            runtime_map_receipts=tuple(runtime_receipts),
        )
        receipts[role] = receipt
    if receipts["current"]["harness_source_root"] != receipts["parent"][
        "harness_source_root"
    ]:
        fail("role receipts disagree on harness root")
    harness_root = Path(receipts["current"]["harness_source_root"])
    source_specs = {
        "harness": (harness_root, config.expected_harness_commit),
        "current": (
            Path(receipts["current"]["implementation_source_root"]),
            CURRENT_IMPLEMENTATION_COMMIT,
        ),
        "parent": (
            Path(receipts["parent"]["implementation_source_root"]),
            PARENT_IMPLEMENTATION_COMMIT,
        ),
    }
    verified_trees: Dict[str, GitTreeReceipt] = {}
    for name, (source_root, commit) in source_specs.items():
        source_authority = authority.get("source_authority", {}).get(name)
        if (
            type(source_authority) is not dict
            or source_authority.get("root") != str(source_root)
            or source_authority.get("commit") != commit
        ):
            fail("execute " + name + " source metadata differs")
        verified = verify_git_tree(
            source_root, commit, required_uid=CAMPAIGN_UID,
            required_gid=CAMPAIGN_GID, sealed=True,
        )
        if (
            verified.tree_oid != source_authority.get("tree_oid")
            or verified.tree_listing_sha256
            != source_authority.get("tree_listing_sha256")
            or verified.manifest_sha256
            != source_authority.get("manifest_sha256")
            or len(verified.entries) != source_authority.get("entries")
        ):
            fail("execute " + name + " source authority differs")
        verified_trees[name] = verified
    validate_frozen_harness_sources(verified_trees["harness"])
    current_harness_entries = [
        entry for entry in receipts["current"]["source_manifest"]
        if entry.get("authority") == "harness"
    ]
    parent_harness_entries = [
        entry for entry in receipts["parent"]["source_manifest"]
        if entry.get("authority") == "harness"
    ]
    if current_harness_entries != parent_harness_entries:
        fail("role receipts disagree on the complete harness manifest")
    controller_entry = _receipt_entry(
        receipts["current"], "harness", CONTROLLER_RELATIVE,
    )
    sampler_entry = _receipt_entry(
        receipts["current"], "harness", SAMPLER_RELATIVE,
    )
    adapter_entry = _receipt_entry(
        receipts["current"], "harness", HEALTH_ADAPTER_RELATIVE,
    )
    if controller_entry.get("sha256") != FROZEN_CONTROLLER_SHA256:
        fail("frozen controller build receipt differs")
    controller = harness_root / CONTROLLER_RELATIVE
    sampler = harness_root / SAMPLER_RELATIVE
    require_artifact(
        controller, controller_entry["sha256"], controller_entry["bytes"],
        "frozen controller", uid=CAMPAIGN_UID, gid=CAMPAIGN_GID, mode=0o444,
    )
    require_artifact(
        sampler, sampler_entry["sha256"], sampler_entry["bytes"],
        "thermal sampler", uid=CAMPAIGN_UID, gid=CAMPAIGN_GID, mode=0o444,
    )
    health_manifest = health_source_manifest_receipt(
        harness_root, receipts["current"],
    )
    health_hash = health_manifest["sha256"]
    health_loader = health_module_loader_receipt(health_manifest)
    health_git = health_source_git_receipt(
        harness_root, config.expected_harness_commit, health_manifest,
    )
    closure_sha256 = {
        role: screen_build_closure_digest(
            receipts[role], role, {
                "harness": verified_trees["harness"],
                "implementation": verified_trees[role],
            },
        )
        for role in ("current", "parent")
    }
    python_image = PYTHON.resolve(strict=True)
    python_digest, python_info = hash_path(
        python_image, 256 * 1024 * 1024, "execute Python image",
        expected_nlink=2,
    )
    if (
        python_digest != FROZEN_TOOL_SHA256[python_image]
        or python_info.st_uid != 0 or python_info.st_gid != 0
        or stat.S_IMODE(python_info.st_mode) != 0o755
    ):
        fail("execute Python image differs")
    return PreparedExecution(
        root_boundary=root_boundary, build_authority=authority,
        build_authority_raw=authority_raw, harness_root=harness_root,
        controller=controller, controller_sha256=controller_entry["sha256"],
        sampler_script=sampler, sampler_sha256=sampler_entry["sha256"],
        health_adapter_sha256=adapter_entry["sha256"],
        health_manifest_sha256=health_hash,
        health_manifest=health_manifest,
        health_module_loader=health_loader,
        health_source_git_receipt=health_git,
        closure_sha256=closure_sha256, roles=roles,
        python_image=python_image, python_sha256=python_digest,
    )


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)
    if args.selftest:
        return selftest()
    if args.seal_build:
        receipt = seal_build_authority(BuildSealConfig(
            harness_source=args.harness_source_dir,
            current_source=args.current_source_dir,
            parent_source=args.parent_source_dir,
            build_root=args.build_root,
            expected_harness_commit=args.expected_harness_commit,
        ))
        print(canonical_bytes({
            "receipt_sha256": receipt,
            "schema": "wirehair.wh2.v2-facade-timing-build-result.v1",
            "status": "sealed",
        }).decode("ascii"))
        return 0
    return execute_once(LaunchConfig(args.expected_harness_commit))


if __name__ == "__main__":
    try:
        result_code = main(sys.argv[1:])
    except Exception as exc:
        print("WH2 V2 facade launcher failed: " + exception_text(exc), file=sys.stderr)
        result_code = 1
    raise SystemExit(result_code)
