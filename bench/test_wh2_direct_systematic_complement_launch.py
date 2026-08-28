#!/usr/bin/python3
"""Scratch-only tests for the one-shot WH2 root supervisor.

No test in this file starts the sampler or scientific worker, touches a fixed
campaign namespace, asks systemd for a unit, or opens a real cgroup control
file.  The orchestration tests use an exact fake backend.
"""

from __future__ import annotations

import ast
import copy
import hashlib
import inspect
import json
import os
from pathlib import Path
import shlex
import shutil
import stat
import subprocess
import sys
import tempfile
import types
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import Wh2DirectSystematicComplementLaunch as subject


COMMIT = "1" * 40
HASH = "a" * 64


def canonical(value):
    return json.dumps(
        value, allow_nan=False, ensure_ascii=True,
        separators=(",", ":"), sort_keys=True,
    ).encode("ascii")


def prepared_fixture() -> subject.Prepared:
    source = Path("/sealed/source")
    config = subject.LaunchConfig(Path("/sealed/build"), COMMIT)
    return subject.Prepared(
        config=config, source_root=source,
        controller=source / subject.CONTROLLER_RELATIVE,
        sampler_script=source / subject.SAMPLER_RELATIVE,
        binary=config.build_dir / subject.BINARY_NAME,
        launcher_sha256=HASH,
        cgroup_path=Path("/sys/fs/cgroup/system.slice") / subject.UNIT_NAME,
        cgroup_relative="/system.slice/" + subject.UNIT_NAME,
        controller_sha256=HASH, controller_bytes=100,
        sampler_sha256=HASH, sampler_bytes=100,
        python_sha256=HASH, setpriv_sha256=HASH, env_sha256=HASH,
        git_sha256=HASH, binary_sha256=HASH,
        build_manifest_sha256=HASH,
        build_authority={"receipt_sha256": HASH},
        build_authority_sha256=HASH, source_manifest_sha256=HASH,
    )


def run_argv_fixture(prepared, sampler_pid=301):
    values = {}
    for option in subject.RUN_OPTION_ORDER:
        values[option] = "1"
    values.update({
        "--binary": str(prepared.binary),
        "--build-dir": str(prepared.config.build_dir),
        "--cpu": str(subject.WORKER_CPU),
        "--controller-cpu": str(subject.CONTROLLER_CPU),
        "--sampler-pid": str(sampler_pid),
        "--sampler-cpu": str(subject.SAMPLER_CPU),
        "--sampler-script": str(prepared.sampler_script),
        "--sampler-csv": str(subject.SAMPLER_CSV),
        "--sampler-pid-file": str(subject.SAMPLER_PID_FILE),
        "--sampler-validation-jsonl": str(subject.SAMPLER_VALIDATION),
        "--sampler-receipt": str(subject.SAMPLER_RECEIPT),
        "--expected-source-commit": prepared.config.expected_commit,
        "--expected-binary-uid": str(subject.CAMPAIGN_UID),
        "--expected-controller-uid": str(subject.CAMPAIGN_UID),
        "--expected-controller-gid": str(subject.CAMPAIGN_GID),
        "--expected-sampler-uid": str(subject.CAMPAIGN_UID),
        "--expected-sampler-gid": str(subject.CAMPAIGN_GID),
        "--expected-sampler-i2c-gid": str(subject.I2C_GID),
    })
    for option in (
        "--expected-binary-sha256", "--expected-build-manifest-sha256",
        "--expected-controller-sha256", "--expected-git-sha256",
        "--expected-python-sha256", "--expected-sampler-script-sha256",
        "--expected-sampler-cmdline-sha256",
        "--expected-sampler-environ-sha256",
        "--expected-sampler-executable-sha256",
        "--expected-source-manifest-sha256",
    ):
        values[option] = HASH
    argv = [str(prepared.controller), "--run-once"]
    for option in subject.RUN_OPTION_ORDER:
        argv.extend((option, values[option]))
    return argv


def preflight_fixture(prepared=None, sampler_pid=301):
    if prepared is None:
        prepared = prepared_fixture()
    argv = run_argv_fixture(prepared, sampler_pid)
    value = {name: {} for name in subject.PREFLIGHT_KEYS}
    value.update({
        "expected_source_commit": prepared.config.expected_commit,
        "receipt_sha256": None,
        "run_argv": argv,
        "run_argv_sha256": subject.argv_sha256(argv),
        "schema": subject.PREFLIGHT_SCHEMA,
        "source_root": str(prepared.source_root),
    })
    value["binary_before"] = {"sha256": prepared.binary_sha256}
    value["binary_after"] = {"sha256": prepared.binary_sha256}
    value["build_manifest_before"] = {
        "sha256": prepared.build_manifest_sha256
    }
    value["build_manifest_after"] = {
        "sha256": prepared.build_manifest_sha256
    }
    value["source_manifest_before"] = {
        "sha256": prepared.source_manifest_sha256
    }
    source_hashes = {relative: HASH for relative in subject.STATIC_SOURCE_PATHS}
    source_hashes[subject.SOURCE_LAUNCHER_RELATIVE] = prepared.launcher_sha256
    source_hashes[subject.CONTROLLER_RELATIVE] = prepared.controller_sha256
    source_hashes[subject.SAMPLER_RELATIVE] = prepared.sampler_sha256
    value["source_manifest_after"] = {
        "sha256": prepared.source_manifest_sha256,
        "entries": [
            {"path": relative, "sha256": source_hashes[relative]}
            for relative in subject.STATIC_SOURCE_PATHS
        ],
    }
    value["receipt_sha256"] = hashlib.sha256(canonical(value)).hexdigest()
    raw = canonical(value) + b"\n"
    return subject.parse_preflight_receipt(raw, prepared, sampler_pid)


def build_command_fixture(name, argv, environment, started):
    stdout = (
        subject.expected_no_work_stdout(argv) if name == "no-work" else ""
    )
    security_fd = 9
    security = subject.self_hashed(subject.BUILD_CHILD_SECURITY_SCHEMA, {
        "affinity": [0, 1], "cap_ambient": 0, "cap_bounding": 0,
        "cap_effective": 0, "cap_inheritable": 0, "cap_permitted": 0,
        "gid": [subject.CAMPAIGN_GID] * 3, "groups": [],
        "no_new_privs": 1, "observed_monotonic_ns": started + 1,
        "pid": 1234, "uid": [subject.CAMPAIGN_UID] * 3,
    })
    return subject.self_hashed(
        "wirehair.wh2.direct-systematic-complement-build-command.v1", {
            "child_pid": 1234, "child_security": security,
            "child_security_fd": security_fd,
            "environment": dict(environment),
            "finished_monotonic_ns": started + 10,
            "finished_utc": "2026-08-28T00:00:01.000000Z",
            "launch_argv": subject.build_setpriv_argv(argv, security_fd),
            "name": name, "returncode": 0,
            "started_monotonic_ns": started,
            "started_utc": "2026-08-28T00:00:00.000000Z",
            "stderr_bytes": 0, "stderr_sha256": hashlib.sha256(b"").hexdigest(),
            "stderr_utf8": "", "stdout_bytes": len(stdout.encode("utf-8")),
            "stdout_sha256": hashlib.sha256(stdout.encode("utf-8")).hexdigest(),
            "stdout_utf8": stdout, "tool_argv": list(argv),
        },
    )


def build_authority_fixture(config):
    environment = subject.build_environment(config.build_dir)
    commands = []
    started = 100
    for name, argv, _timeout in subject.build_command_vectors(config):
        commands.append(build_command_fixture(name, argv, environment, started))
        started += 20
    directory = {"sealed": "build"}
    subdirectories = [{"sealed": "subdirectories"}]
    installed = {"sha256": HASH}
    profile = {"sealed": "profile"}
    binary = {"sealed": "binary"}
    tools = [{"sealed": "toolchain"}]
    value = subject.self_hashed(subject.BUILD_AUTHORITY_SCHEMA, {
        "binary": binary, "build_directory": directory,
        "build_environment": environment,
        "build_manifest_sha256": HASH, "build_profile": profile,
        "build_subdirectories": subdirectories, "commands": commands,
        "expected_source_commit": config.expected_commit,
        "git_sha256_after": HASH, "git_sha256_before": HASH,
        "installed_launcher": installed,
        "root_boundary": {
            "affinity": [0, 1],
            "cgroup": "/system.slice/" + subject.BUILD_UNIT_NAME,
            "gid": 0, "groups": [], "uid": 0,
        },
        "source_manifest_sha256_after": HASH,
        "source_manifest_sha256_before": HASH,
        "source_root": str(config.source_dir), "status": "sealed",
        "systemd_run_argv": subject.build_systemd_run_argv(config),
        "toolchain_after": tools, "toolchain_before": tools,
        "trust_boundary": subject.BUILD_TRUST_BOUNDARY,
    })
    observations = {
        "binary": binary, "directory": directory, "installed": installed,
        "profile": profile, "subdirectories": subdirectories, "tools": tools,
    }
    return value, observations


def write_build_profile_fixture(build, source_root, commit=COMMIT):
    (build / "CMakeFiles").mkdir(parents=True)
    cache = {
        "BUILD_CODEC_V2": ("BOOL", "ON"),
        "BUILD_SHARED_LIBS": ("BOOL", "OFF"),
        "BUILD_TESTS": ("BOOL", "ON"),
        "CMAKE_AR": ("FILEPATH", "/usr/bin/x86_64-linux-gnu-ar"),
        "CMAKE_BUILD_TYPE": ("STRING", "Release"),
        "CMAKE_CACHEFILE_DIR": ("INTERNAL", str(build)),
        "CMAKE_CXX_COMPILER": ("STRING", str(subject.CXX_COMPILER)),
        "CMAKE_CXX_COMPILER_AR": (
            "FILEPATH", "/usr/bin/x86_64-linux-gnu-gcc-ar-13"
        ),
        "CMAKE_CXX_COMPILER_RANLIB": (
            "FILEPATH", "/usr/bin/x86_64-linux-gnu-gcc-ranlib-13"
        ),
        "CMAKE_CXX_FLAGS": ("STRING", ""),
        "CMAKE_CXX_FLAGS_RELEASE": ("STRING", "-O3 -DNDEBUG"),
        "CMAKE_C_COMPILER": ("STRING", str(subject.C_COMPILER)),
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
        "CMAKE_PROJECT_NAME": ("STATIC", "wirehair"),
        "CMAKE_RANLIB": ("FILEPATH", "/usr/bin/x86_64-linux-gnu-ranlib"),
        "CMAKE_SHARED_LINKER_FLAGS": ("STRING", ""),
        "CMAKE_SHARED_LINKER_FLAGS_RELEASE": ("STRING", ""),
        "CMAKE_STATIC_LINKER_FLAGS": ("STRING", ""),
        "CMAKE_STATIC_LINKER_FLAGS_RELEASE": ("STRING", ""),
        "CMAKE_STRIP": ("FILEPATH", "/usr/bin/x86_64-linux-gnu-strip"),
        "MARCH_NATIVE": ("BOOL", "OFF"),
        "Python3_EXECUTABLE": ("UNINITIALIZED", str(subject.PYTHON)),
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
        "_Python3_EXECUTABLE": ("INTERNAL", str(subject.PYTHON)),
        "wirehair_BINARY_DIR": ("STATIC", str(build)),
        "wirehair_SOURCE_DIR": ("STATIC", str(source_root)),
    }
    cache_raw = "".join(
        name + ":" + kind + "=" + value + "\n"
        for name, (kind, value) in sorted(cache.items())
    ).encode("utf-8")
    object_prefix = "CMakeFiles/" + subject.BINARY_NAME + ".dir/bench/"
    commands = []
    build_lines = []
    for relative in subject.COMPLEMENT_BUILD_SOURCES:
        source = source_root / relative
        source.parent.mkdir(parents=True, exist_ok=True)
        source.write_bytes(b"// sealed fixture\n")
        output = object_prefix + Path(relative).name + ".o"
        command = " ".join((
            str(subject.CXX_COMPILER), "-DWIREHAIR_STATIC=1",
            "-DWIREHAIR_V2_ENABLE_BENCHMARK_EQUATIONS=1",
            "-DWIREHAIR_WH2_SOURCE_GIT_COMMIT=\\\"" + commit + "\\\"",
            "-I" + str(source_root / "codec"),
            "-I" + str(source_root / "include"),
            "-O3", "-DNDEBUG", "-std=gnu++11", "-Wall", "-Wextra",
            "-Wpedantic", "-Werror", "-o", output, "-c", str(source),
        ))
        commands.append({
            "directory": str(build), "command": command,
            "file": str(source), "output": output,
        })
        build_lines.append(
            "build " + output + ": CXX_COMPILER__" + subject.BINARY_NAME
            + "_unscanned_Release " + str(source)
            + " || cmake_object_order_depends_target_" + subject.BINARY_NAME + "\n"
            "  DEFINES = -DWIREHAIR_STATIC=1"
            " -DWIREHAIR_V2_ENABLE_BENCHMARK_EQUATIONS=1"
            " -DWIREHAIR_WH2_SOURCE_GIT_COMMIT=\\\"" + commit + "\\\"\n"
            "  DEP_FILE = " + output + ".d\n"
            "  FLAGS = -O3 -DNDEBUG -std=gnu++11 -Wall -Wextra"
            " -Wpedantic -Werror\n"
            "  INCLUDES = -I" + str(source_root / "codec")
            + " -I" + str(source_root / "include") + "\n"
            "  OBJECT_DIR = CMakeFiles/" + subject.BINARY_NAME + ".dir\n"
            "  OBJECT_FILE_DIR = CMakeFiles/" + subject.BINARY_NAME
            + ".dir/bench\n"
            "  TARGET_COMPILE_PDB = CMakeFiles/" + subject.BINARY_NAME
            + ".dir/\n"
            "  TARGET_PDB = " + subject.BINARY_NAME + ".pdb\n\n"
        )
    direct_outputs = [item["output"] for item in commands]
    for relative in subject.CORE_BUILD_SOURCES:
        source = source_root / relative
        source.parent.mkdir(parents=True, exist_ok=True)
        source.write_bytes(b"// sealed core fixture\n")
        output = "CMakeFiles/wirehair.dir/" + relative + ".o"
        command = " ".join((
            str(subject.CXX_COMPILER), "-DWIREHAIR_BUILDING=1",
            "-DWIREHAIR_STATIC=1", "-I" + str(source_root / "include"),
            "-O3", "-DNDEBUG", "-std=gnu++11", "-fPIC", "-Wall",
            "-Wextra", "-Wpedantic", "-Werror", "-o", output, "-c",
            str(source),
        ))
        commands.append({
            "directory": str(build), "command": command,
            "file": str(source), "output": output,
        })
        build_lines.append(
            "build " + output + ": CXX_COMPILER__fixture_unscanned_Release "
            + str(source) + "\n\n"
        )
    for relative in subject.TIMING_POLICY_BUILD_SOURCES:
        source = source_root / relative
        source.parent.mkdir(parents=True, exist_ok=True)
        source.write_bytes(b"// sealed timing fixture\n")
        output = "codec/CMakeFiles/wirehair_v2_timing_policy.dir/" \
            + Path(relative).name + ".o"
        command = " ".join((
            str(subject.CXX_COMPILER), "-DWIREHAIR_STATIC=1",
            "-DWIREHAIR_V2_ENABLE_BENCHMARK_EQUATIONS=1",
            "-I" + str(source_root / "codec"),
            "-I" + str(source_root / "include"), "-O3", "-DNDEBUG",
            "-std=gnu++11", "-Wall", "-Wextra", "-Wpedantic", "-Werror",
            "-o", output, "-c", str(source),
        ))
        commands.append({
            "directory": str(build), "command": command,
            "file": str(source), "output": output,
        })
        build_lines.append(
            "build " + output + ": CXX_COMPILER__fixture_unscanned_Release "
            + str(source) + "\n\n"
        )
    link_line = (
        "build " + subject.BINARY_NAME + ": CXX_EXECUTABLE_LINKER__"
        + subject.BINARY_NAME + "_Release " + " ".join(direct_outputs)
        + " | codec/libwirehair_v2_timing_policy.a libwirehair.a"
        + " || codec/libwirehair_v2_timing_policy.a libwirehair.a\n"
        "  FLAGS = -O3 -DNDEBUG\n"
        "  LINK_LIBRARIES = codec/libwirehair_v2_timing_policy.a  libwirehair.a  -lm\n"
        "  OBJECT_DIR = CMakeFiles/" + subject.BINARY_NAME + ".dir\n"
        "  POST_BUILD = :\n"
        "  PRE_LINK = :\n"
        "  TARGET_COMPILE_PDB = CMakeFiles/" + subject.BINARY_NAME + ".dir/\n"
        "  TARGET_FILE = " + subject.BINARY_NAME + "\n"
        "  TARGET_PDB = " + subject.BINARY_NAME + ".pdb\n\n"
    )
    build_prefix = (
        "# Configurations: Release\n"
        "CONFIGURATION = Release\n"
        "cmake_ninja_workdir = " + str(build) + "/\n"
    )
    build_lines.append(link_line)
    rules = (
        "rule CXX_COMPILER__" + subject.BINARY_NAME + "_unscanned_Release\n"
        "  depfile = $DEP_FILE\n"
        "  deps = gcc\n"
        "  command = ${LAUNCHER}${CODE_CHECK}" + str(subject.CXX_COMPILER)
        + " $DEFINES $INCLUDES $FLAGS -MD -MT $out -MF $DEP_FILE"
        " -o $out -c $in\n"
        "  description = Building CXX object $out\n\n"
        "rule CXX_EXECUTABLE_LINKER__" + subject.BINARY_NAME + "_Release\n"
        "  command = $PRE_LINK && " + str(subject.CXX_COMPILER)
        + " $FLAGS $LINK_FLAGS $in -o $TARGET_FILE $LINK_PATH"
        " $LINK_LIBRARIES && $POST_BUILD\n"
        "  description = Linking CXX executable $TARGET_FILE\n"
        "  restat = $RESTAT\n\n"
    ).encode("utf-8")
    values = {
        "CMakeCache.txt": cache_raw,
        "compile_commands.json": json.dumps(commands).encode("utf-8"),
        "build.ninja": (build_prefix + "".join(build_lines)).encode("utf-8"),
        "CMakeFiles/rules.ninja": rules,
    }
    for relative, raw in values.items():
        path = build / relative
        path.write_bytes(raw)
        path.chmod(0o444)
    return values


def effective_ninja_commands_fixture(build):
    commands = json.loads((build / "compile_commands.json").read_text("utf-8"))
    lines = []
    direct_outputs = []
    for item in commands:
        tokens = shlex.split(item["command"], posix=True)
        output = item["output"]
        output_index = tokens.index("-o")
        tokens[output_index:output_index] = [
            "-MD", "-MT", output, "-MF", output + ".d",
        ]
        lines.append(shlex.join(tokens))
        if output.startswith(
            "CMakeFiles/" + subject.BINARY_NAME + ".dir/bench/"
        ):
            direct_outputs.append(output)
    archive_objects = {
        "libwirehair.a": [
            "CMakeFiles/wirehair.dir/" + relative + ".o"
            for relative in subject.CORE_BUILD_SOURCES
        ],
        "codec/libwirehair_v2_timing_policy.a": [
            "codec/CMakeFiles/wirehair_v2_timing_policy.dir/"
            + Path(relative).name + ".o"
            for relative in subject.TIMING_POLICY_BUILD_SOURCES
        ],
    }
    for archive, objects in archive_objects.items():
        lines.append(
            ": && /usr/bin/cmake -E rm -f " + archive
            + " && /usr/bin/x86_64-linux-gnu-ar qc " + archive + "  "
            + " ".join(objects)
            + " && /usr/bin/x86_64-linux-gnu-ranlib " + archive + " && :"
        )
    lines.append(
        ": && " + str(subject.CXX_COMPILER) + " -O3 -DNDEBUG  "
        + " ".join(direct_outputs) + " -o " + subject.BINARY_NAME
        + "  codec/libwirehair_v2_timing_policy.a  libwirehair.a  -lm && :"
    )
    return lines


class FakeJournal:
    SUCCESS_ROSTER = [
        "build-authority.json", "preflight.json", "authority.json",
        "start.json", "deadline.json",
        "sampler-terminal.json", "verification.json", "terminal.json",
        "COMPLETE",
    ]

    def __init__(self, fail_name=None, delay_name=None):
        self.fail_name = fail_name
        self.delay_name = delay_name
        self.names = []
        self.payloads = {}
        self.finished = False
        self.complete = False
        self.closed = False

    def _append(self, name, payload):
        if self.fail_name == name:
            self.fail_name = None
            raise subject.LaunchError("injected journal failure " + name)
        if name in self.names:
            raise AssertionError("duplicate fake journal member")
        self.names.append(name)
        self.payloads[name] = payload

    def write_bytes(self, name, payload, final_mode=0o400):
        self._append(name, bytes(payload))
        return {"bytes": len(payload), "mode": final_mode}

    def write_json(self, name, value):
        if self.delay_name == name:
            subject.time.sleep(0.01)
        self._append(name, canonical(dict(value)) + b"\n")
        return {"bytes": len(self.payloads[name]), "mode": 0o400}

    def finish(self):
        if not self.names or self.names[-1] != "terminal.json":
            raise AssertionError("fake COMPLETE was not terminal-last")
        if self.fail_name == "COMPLETE":
            self.fail_name = None
            raise OSError("injected COMPLETE fsync failure")
        self.names.append("COMPLETE")
        self.finished = True
        self.complete = True

    def close(self):
        self.closed = True
        return True


class FakeGuardian:
    def __init__(self):
        self.error = "none"
        self.fired = False
        self.backstop_fired = False
        self.roster = []
        self.cancelled = False
        self.finished = False
        self.closed = False
        self.close_fails = False
        self.cancel_fails = False
        self.finish_fails = False
        self.process = types.SimpleNamespace(pid=304, start_ticks=3004)
        self.admission_armed = False
        self.exact_armed = False
        self.exact_arm_result = True

    def arm_admission(self, release_not_before_ns):
        self.admission_armed = True
        self.admission_anchor = release_not_before_ns

    def arm_exact(self, controller_t0_ns):
        self.exact_armed = True
        if not self.exact_arm_result:
            self.fired = True
        return self.exact_arm_result

    def cancel(self):
        if self.cancel_fails:
            raise subject.LaunchError("injected guardian cancel failure")
        self.cancelled = True

    def finish(self):
        if self.finish_fails:
            raise subject.LaunchError("injected guardian finish failure")
        self.finished = True

    def close(self):
        if self.close_fails:
            raise subject.LaunchError("injected guardian channel close failure")
        self.closed = True
        return True


class FakeBackend:
    def __init__(self, classification="pass", controller_rc=0,
                 fail_at=None, journal_fail=None, journal_delay=None):
        self.classification = classification
        self.controller_rc = controller_rc
        self.fail_at = fail_at
        self.calls = []
        self.prepared = prepared_fixture()
        self.journal = FakeJournal(journal_fail, journal_delay)
        self.guardian = FakeGuardian()
        self.cleaned = False
        self.runtime_closed = False
        self.authority_closed = False
        self.quiesce_ok = True
        self.quiesce_clears_verifier = False
        self.sampler_stop_fails = False
        self.verifier_containment_fails = False
        self.verifier_live = False
        self.attempt_consumed = False
        self.experiment_roster = []
        self.processes = {
            "sampler": types.SimpleNamespace(
                role="sampler", pid=301, start_ticks=3001),
            "controller": types.SimpleNamespace(
                role="controller", pid=302, start_ticks=3002),
            "verifier": types.SimpleNamespace(
                role="verifier", pid=303, start_ticks=3003),
        }

    def _call(self, name):
        self.calls.append(name)
        if self.fail_at == name:
            self.fail_at = None
            raise subject.LaunchError("injected " + name)

    def prepare(self, config):
        self._call("prepare")
        self.prepared = subject.Prepared(
            config=config, source_root=self.prepared.source_root,
            controller=self.prepared.controller,
            sampler_script=self.prepared.sampler_script,
            binary=config.build_dir / subject.BINARY_NAME,
            launcher_sha256=HASH, cgroup_path=self.prepared.cgroup_path,
            cgroup_relative=self.prepared.cgroup_relative,
            controller_sha256=HASH, controller_bytes=100,
            sampler_sha256=HASH, sampler_bytes=100,
            python_sha256=HASH, setpriv_sha256=HASH, env_sha256=HASH,
            git_sha256=HASH, binary_sha256=HASH,
            build_manifest_sha256=HASH,
            build_authority={"receipt_sha256": HASH},
            build_authority_sha256=HASH, source_manifest_sha256=HASH,
        )
        return self.prepared

    def reserve(self, prepared):
        self.calls.append("reserve")
        if self.fail_at == "reserve_after_consumed":
            self.attempt_consumed = True
            self.fail_at = None
            raise subject.AttemptConsumedError(
                "injected post-mkdir reservation failure", self.journal
            )
        if self.fail_at == "reserve":
            self.fail_at = None
            raise subject.LaunchError("injected reserve")
        self.attempt_consumed = True
        return self.journal

    def start_sampler(self, prepared):
        self._call("start_sampler")
        return self.processes["sampler"]

    def wait_sampler_ready(self, sampler, prepared):
        self._call("wait_sampler_ready")
        return {
            "devices": [], "observed_monotonic_ns": 1,
            "scope": "trusted-root-point-in-time-no-kernel-exclusion",
        }

    def run_preflight(self, prepared, sampler):
        self._call("run_preflight")
        argv = run_argv_fixture(prepared, sampler.pid)
        raw = b'{"synthetic":"preflight"}\n'
        return subject.PreflightAuthority(
            raw=raw, receipt={}, run_argv=argv,
            run_argv_sha256=subject.argv_sha256(argv),
            receipt_sha256=HASH,
        )

    def spawn_stubs(self, prepared, authority):
        self._call("spawn_stubs")
        self.verifier_live = True
        return self.processes["controller"], self.processes["verifier"]

    def spawn_guardian(self, journal):
        self._call("spawn_guardian")
        return self.guardian

    def release(self, process):
        self._call("release_" + process.role)
        return subject.time.monotonic_ns()

    def wait_claim(self, authority, controller, release_ns):
        self._call("wait_claim")
        return {}, b'{"synthetic":"claim"}\n', release_ns + 1

    def controller_wait(self, controller, deadline_ns):
        self._call("controller_wait")
        return self.controller_rc, b"", b""

    def ensure_experiment_empty(self, controller):
        self._call("ensure_experiment_empty")
        roster = list(self.experiment_roster)
        return bool(roster), roster

    def stop_sampler(self, sampler, prepared):
        self._call("stop_sampler")
        if self.verifier_live:
            raise AssertionError("sampler stopped while verifier remained live")
        if self.sampler_stop_fails:
            raise subject.LaunchError("injected persistent sampler terminal failure")
        return 0, {
            "raw": b'{"synthetic":"sampler"}\n',
            "receipt_sha256": HASH,
        }, b"", b""

    def run_verifier(self, verifier):
        self._call("run_verifier")
        self.verifier_live = False
        rc = subject.VERIFY_EXIT[self.classification]
        raw = canonical({
            "outcome": self.classification, "schema": subject.VERIFY_SCHEMA,
            "status": (
                "absent" if self.classification == "absent" else "verified"
            ),
        }) + b"\n"
        return self.classification, rc, raw, b""

    def contain_verifier(self, verifier):
        self._call("contain_verifier")
        if self.verifier_containment_fails:
            raise subject.LaunchError("injected persistent verifier containment failure")
        self.verifier_live = False
        return True

    def finish_sampler_dir(self):
        self._call("finish_sampler_dir")

    def kill_experiment(self):
        self.calls.append("kill_experiment")
        return []

    def cancel(self, process):
        self.calls.append("cancel_" + process.role)

    def quiesce(self):
        self.calls.append("quiesce")
        self.cleaned = True
        if self.quiesce_clears_verifier:
            self.verifier_live = False
        if not self.quiesce_ok or self.verifier_live:
            raise subject.LaunchError("injected live run descendant")
        return True

    def close_runtime(self):
        self._call("close_runtime")
        self.runtime_closed = True
        return True

    def close_authority(self):
        self._call("close_authority")
        if self.journal.complete:
            raise AssertionError("authority was closed after COMPLETE")
        self.authority_closed = True
        return True


class ParserTests(unittest.TestCase):
    def test_all_tracked_git_index_flags_must_be_normal(self):
        subject.require_normal_git_index_flags(b"H CMakeLists.txt\nH src/a.cpp\n")
        for raw in (
            b"", b"H no-final-newline", b"h assume-unchanged\n",
            b"S skip-worktree\n", b"H ok\ns skip-worktree\n",
        ):
            with self.subTest(raw=raw), self.assertRaises(subject.LaunchError):
                subject.require_normal_git_index_flags(raw)

    def test_static_git_gate_hashes_every_tracked_blob_despite_clean_stat_cache(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "source"
            root.mkdir()

            def git(*arguments, capture=False):
                result = subprocess.run(
                    [str(subject.GIT), *arguments], cwd=str(root), check=True,
                    stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
                return result.stdout if capture else b""

            git("init", "-q")
            git("config", "user.email", "test@example.invalid")
            git("config", "user.name", "WH2 test")
            static_path = root / "static.txt"
            dependency = root / "dependency.cpp"
            static_path.write_bytes(b"STATIC\n")
            dependency.write_bytes(b"GOOD\n")
            static_path.chmod(0o644)
            dependency.chmod(0o644)
            old_ns = 1_600_000_000_123_456_789
            os.utime(static_path, ns=(old_ns, old_ns))
            os.utime(dependency, ns=(old_ns, old_ns))
            git("add", "--", "static.txt", "dependency.cpp")
            git("commit", "-q", "-m", "fixture")
            commit = git("rev-parse", "HEAD", capture=True).decode("ascii").strip()
            with mock.patch.object(subject, "STATIC_SOURCE_PATHS", ("static.txt",)):
                self.assertRegex(subject.static_git_gate(root, commit), r"^[0-9a-f]{64}$")
                git("config", "core.trustctime", "false")
                git("config", "core.checkStat", "minimal")
                git("update-index", "--refresh")
                before = dependency.stat()
                dependency.write_bytes(b"EVIL\n")
                os.utime(
                    dependency,
                    ns=(before.st_atime_ns, before.st_mtime_ns),
                )
                self.assertEqual(
                    git(
                        "status", "--porcelain=v1", "--untracked-files=all",
                        capture=True,
                    ),
                    b"",
                )
                with self.assertRaisesRegex(
                    subject.LaunchError,
                    "static Git worktree bytes differ from HEAD: dependency.cpp",
                ):
                    subject.static_git_gate(root, commit)

    def test_controller_and_launcher_manifest_constants_match_exactly(self):
        def literal_assignment(path, name):
            tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
            environment = {}
            for node in tree.body:
                if not isinstance(node, ast.Assign) or len(node.targets) != 1:
                    continue
                target = node.targets[0]
                if not isinstance(target, ast.Name):
                    continue
                try:
                    value = ast.literal_eval(node.value)
                except ValueError:
                    if isinstance(node.value, (ast.Tuple, ast.List)):
                        value = tuple(
                            environment[item.id]
                            if isinstance(item, ast.Name)
                            else ast.literal_eval(item)
                            for item in node.value.elts
                        )
                    else:
                        continue
                environment[target.id] = value
            self.assertIn(name, environment)
            return environment[name]

        bench = Path(subject.__file__).resolve().parent
        controller = bench / "Wh2DirectSystematicComplementScreen.py"
        self.assertEqual(
            literal_assignment(controller, "SOURCE_PATHS"),
            subject.STATIC_SOURCE_PATHS,
        )
        self.assertEqual(
            literal_assignment(controller, "BUILD_PATHS"), subject.BUILD_PATHS,
        )
        self.assertEqual(len(subject.STATIC_SOURCE_PATHS), len(set(subject.STATIC_SOURCE_PATHS)))
        self.assertEqual(len(subject.BUILD_PATHS), len(set(subject.BUILD_PATHS)))

    def test_exact_public_modes(self):
        args = subject.parse_args([
            "--execute", "--build-dir", "/sealed/build",
            "--expected-source-commit", COMMIT,
        ])
        self.assertEqual(args.build_dir, Path("/sealed/build"))
        self.assertEqual(args.expected_source_commit, COMMIT)
        self.assertTrue(subject.parse_args(["--selftest"]).selftest)
        build_args = subject.parse_args([
            "--seal-build", "--source-dir", "/sealed/source",
            "--build-dir", "/sealed/build",
            "--expected-source-commit", COMMIT,
        ])
        self.assertTrue(build_args.seal_build)
        self.assertEqual(build_args.source_dir, Path("/sealed/source"))

    def test_public_mode_rejects_reorder_equals_missing_extra_and_bad_commit(self):
        cases = [
            ["--execute", "--expected-source-commit", COMMIT,
             "--build-dir", "/sealed/build"],
            ["--execute", "--build-dir=/sealed/build",
             "--expected-source-commit", COMMIT],
            ["--execute", "--build-dir", "/sealed/build"],
            ["--execute", "--build-dir", "/sealed/build",
             "--expected-source-commit", COMMIT, "extra"],
            ["--execute", "--build-dir", "relative",
             "--expected-source-commit", COMMIT],
            ["--execute", "--build-dir", "/sealed/build",
             "--expected-source-commit", "A" * 40],
        ]
        for argv in cases:
            with self.subTest(argv=argv), self.assertRaises(subject.LaunchError):
                subject.parse_args(argv)

    def test_systemd_command_has_exact_safety_properties(self):
        argv = subject.systemd_run_argv(Path("/sealed/build"), COMMIT)
        properties = [item for item in argv if item.startswith("--property=")]
        self.assertEqual(len(properties), len(set(properties)))
        for expected in (
            "--property=Restart=no", "--property=KillMode=control-group",
            "--property=ExitType=cgroup",
            "--property=RuntimeMaxSec=240s",
            "--property=Delegate=cpuset pids",
            "--property=AllowedCPUs=120-123",
            "--property=CPUAffinity=123",
            "--property=DevicePolicy=closed",
            "--property=DeviceAllow=/dev/i2c-1 rw",
            "--property=DeviceAllow=/dev/i2c-2 rw",
        ):
            self.assertIn(expected, argv)
        python_index = argv.index(str(subject.PYTHON))
        self.assertEqual(argv[python_index:python_index + 5], [
            str(subject.PYTHON), "-I", "-B",
            str(subject.INSTALLED_LAUNCHER), "--execute",
        ])
        env_index = argv.index(str(subject.ENV_EXECUTABLE))
        self.assertEqual(argv[env_index:env_index + 7], [
            str(subject.ENV_EXECUTABLE), "-i", "LANG=C.UTF-8",
            "LC_ALL=C.UTF-8", "PATH=/usr/bin:/bin", "TZ=UTC",
            str(subject.PYTHON),
        ])
        self.assertEqual(argv[-5:], [
            "--execute", "--build-dir", "/sealed/build",
            "--expected-source-commit", COMMIT,
        ])

    def test_build_systemd_command_is_root_scrubbed_and_private(self):
        config = subject.BuildSealConfig(
            Path("/sealed/source"), Path("/sealed/build"), COMMIT,
        )
        argv = subject.build_systemd_run_argv(config)
        for expected in (
            "--property=RuntimeMaxSec=2400s",
            "--property=AllowedCPUs=0-1",
            "--property=CPUAffinity=0-1",
            "--property=PrivateTmp=yes",
            "--property=PrivateDevices=yes",
            "--property=DevicePolicy=closed",
            "--property=ProtectProc=invisible",
        ):
            self.assertIn(expected, argv)
        self.assertEqual(argv[-7:], [
            "--seal-build", "--source-dir", "/sealed/source",
            "--build-dir", "/sealed/build",
            "--expected-source-commit", COMMIT,
        ])

    def test_main_dispatches_build_without_constructing_scientific_backend(self):
        with mock.patch.object(
            subject, "seal_build_authority", return_value=HASH,
        ) as seal, mock.patch.object(
            subject, "RealBackend", side_effect=AssertionError("scientific backend")
        ), mock.patch("builtins.print") as output:
            result = subject.main([
                "--seal-build", "--source-dir", "/sealed/source",
                "--build-dir", "/sealed/build",
                "--expected-source-commit", COMMIT,
            ])
        self.assertEqual(result, 0)
        seal.assert_called_once_with(subject.BuildSealConfig(
            Path("/sealed/source"), Path("/sealed/build"), COMMIT,
        ))
        self.assertIn(HASH, output.call_args.args[0])

    def test_cmake_source_root_parser_is_exact(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary).resolve()
            build = root / "build"
            build.mkdir()
            raw = (
                "x:y=z\nCMAKE_HOME_DIRECTORY:INTERNAL=" + str(root) + "\n"
            ).encode("utf-8")
            self.assertEqual(subject.parse_cmake_source_root(raw, build), root)
            with self.assertRaises(subject.LaunchError):
                subject.parse_cmake_source_root(raw + raw, build)

    def test_cgroup_and_proc_stat_parsers(self):
        raw = ("0::/system.slice/" + subject.UNIT_NAME + "\n").encode("ascii")
        self.assertEqual(
            subject.parse_proc_cgroup(raw),
            "/system.slice/" + subject.UNIT_NAME,
        )
        with self.assertRaises(subject.LaunchError):
            subject.parse_proc_cgroup(raw + raw)
        fields = [b"S"] + [b"1"] * 18 + [b"777"] + [b"1"] * 8
        self.assertEqual(
            subject.proc_start_ticks(
                b"123 (comm with spaces) " + b" ".join(fields) + b"\n"
            ),
            777,
        )
        with self.assertRaises(subject.LaunchError):
            subject.proc_start_ticks(
                b"123 (comm) " + b" ".join(fields) + b"\n\n"
            )


class ReceiptTests(unittest.TestCase):
    def test_i2c_point_in_time_holder_policy_is_exact(self):
        scope = "trusted-root-point-in-time-no-kernel-exclusion"
        empty = {
            "scope": scope, "observed_monotonic_ns": 1,
            "devices": [
                {"path": "/dev/i2c-1", "holders": []},
                {"path": "/dev/i2c-2", "holders": []},
            ],
        }
        subject.require_no_i2c_holders(empty)
        sampler = types.SimpleNamespace(pid=123, start_ticks=456)
        sole = copy.deepcopy(empty)
        for index, device in enumerate(sole["devices"]):
            device["holders"] = [{
                "fd": 7 + index, "pid": sampler.pid,
                "process_start_ticks": sampler.start_ticks,
            }]
        subject.require_sampler_sole_i2c_holder(sole, sampler)
        sole["devices"][0]["holders"].append({
            "fd": 99, "pid": 999, "process_start_ticks": 1,
        })
        with self.assertRaises(subject.LaunchError):
            subject.require_sampler_sole_i2c_holder(sole, sampler)

    def test_claim_reader_waits_for_publish_mode_and_then_binds_inode(self):
        with tempfile.TemporaryDirectory() as temporary:
            claim = Path(temporary) / "claim"
            claim.write_bytes(b'{"ready":true}\n')
            claim.chmod(0o600)
            with mock.patch.object(subject, "FIXED_CLAIM_PATH", claim):
                self.assertIsNone(subject.RealBackend._read_claim_if_ready())
                claim.chmod(0o400)
                self.assertEqual(
                    subject.RealBackend._read_claim_if_ready(),
                    b'{"ready":true}\n',
                )

    def test_preflight_exact_fixture_passes(self):
        authority = preflight_fixture()
        self.assertEqual(authority.receipt_sha256, authority.receipt["receipt_sha256"])
        self.assertEqual(authority.run_argv_sha256, subject.argv_sha256(authority.run_argv))

    def test_preflight_rejects_duplicate_noncanonical_and_tamper(self):
        authority = preflight_fixture()
        duplicate = b'{"schema":"duplicate",' + authority.raw[1:]
        with self.assertRaisesRegex(subject.LaunchError, "duplicate"):
            subject.parse_preflight_receipt(duplicate, prepared_fixture(), 301)
        noncanonical = authority.raw.replace(b'":', b'": ', 1)
        with self.assertRaises(subject.LaunchError):
            subject.parse_preflight_receipt(noncanonical, prepared_fixture(), 301)
        value = copy.deepcopy(authority.receipt)
        index = value["run_argv"].index("--cpu") + 1
        value["run_argv"][index] = "119"
        value["run_argv_sha256"] = subject.argv_sha256(value["run_argv"])
        value["receipt_sha256"] = None
        value["receipt_sha256"] = hashlib.sha256(canonical(value)).hexdigest()
        with self.assertRaisesRegex(subject.LaunchError, "--cpu"):
            subject.parse_preflight_receipt(
                canonical(value) + b"\n", prepared_fixture(), 301,
            )

    def test_run_argv_rejects_reordered_and_non_ascii_values(self):
        argv = run_argv_fixture(prepared_fixture())
        argv[2], argv[4] = argv[4], argv[2]
        with self.assertRaises(subject.LaunchError):
            subject.option_values(
                argv, "--run-once", subject.RUN_OPTION_ORDER, "test run",
            )
        argv = run_argv_fixture(prepared_fixture())
        argv[3] = "\N{SNOWMAN}"
        with self.assertRaises(subject.LaunchError):
            subject.option_values(
                argv, "--run-once", subject.RUN_OPTION_ORDER, "test run",
            )

    def test_role_command_cross_binds_preflight_and_self_hash(self):
        authority = preflight_fixture()
        raw = subject.role_document("controller", authority.run_argv, authority)
        parsed = subject.validate_role_document(raw, "controller", authority.raw)
        self.assertEqual(parsed["python_argv"], authority.run_argv)
        with self.assertRaises(subject.LaunchError):
            subject.validate_role_document(raw, "verifier", authority.raw)
        with self.assertRaises(subject.LaunchError):
            subject.validate_role_document(raw, "controller", authority.raw + b"x")

    def test_verifier_exit_and_output_are_one_binding(self):
        for outcome, code in subject.VERIFY_EXIT.items():
            raw = canonical({
                "outcome": outcome, "schema": subject.VERIFY_SCHEMA,
                "status": "absent" if outcome == "absent" else "verified",
            }) + b"\n"
            self.assertEqual(subject.parse_verifier_output(raw, code), outcome)
            with self.assertRaises(subject.LaunchError):
                subject.parse_verifier_output(raw, 1)

    def test_setpriv_roles_have_exact_identity_and_python_flags(self):
        controller = subject.setpriv_exec_argv(["/controller", "--run-once"], sampler=False)
        sampler = subject.setpriv_exec_argv(["/sampler", "--csv", "/x"], sampler=True)
        self.assertIn("--clear-groups", controller)
        self.assertNotIn("--groups", controller)
        self.assertEqual(controller[-4:-2], ["-I", "-B"])
        self.assertIn("--groups", sampler)
        self.assertIn(str(subject.I2C_GID), sampler)
        self.assertEqual(sampler[-6:-3], ["-I", "-S", "-B"])
        for argv in (controller, sampler):
            self.assertIn("--bounding-set=-all", argv)
            self.assertIn("--ambient-caps=-all", argv)
            self.assertIn("--no-new-privs", argv)
            self.assertIn("SIGTERM", argv)

    def test_held_controller_inode_survives_path_replacement(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "controller.py"
            source.write_text(
                "def main(argv):\n"
                " print('held:' + __file__ + ':' + ','.join(argv))\n"
                " return 0\n",
                encoding="ascii",
            )
            original = source.read_bytes()
            held_fd = os.open(source, os.O_RDONLY | os.O_CLOEXEC)
            try:
                source.unlink()
                source.write_text(
                    "def main(argv):\n print('substituted')\n return 0\n",
                    encoding="ascii",
                )
                results = [
                    subprocess.run(
                        [
                            sys.executable, "-I", "-B", "-c",
                            subject.SEALED_MODULE_BOOTSTRAP, str(held_fd),
                            str(source), hashlib.sha256(original).hexdigest(),
                            str(len(original)), "main", "--verify-retained",
                            "token",
                        ],
                        pass_fds=(held_fd,), stdin=subprocess.DEVNULL,
                        stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                        timeout=5.0, check=False,
                    )
                    for _ in range(2)
                ]
            finally:
                os.close(held_fd)
            for result in results:
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertEqual(
                    result.stdout,
                    (
                        "held:" + str(source)
                        + ":--verify-retained,token\n"
                    ).encode("ascii"),
                )
                self.assertNotIn(b"substituted", result.stdout)

    def test_authentic_sampler_terminal_accepts_unlinked_held_pid_inode(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampler_source = root / "sampler.py"
            csv_path = root / "thermal.csv"
            validation_path = root / "validation.jsonl"
            pid_path = root / "sampler.pid"
            receipt_path = root / "sampler-terminal.json"
            sampler_source.write_bytes(b"# sealed sampler fixture\n")
            csv_path.write_bytes(b"utc,value\n2026-01-01T00:00:00Z,1\n")
            validation_path.write_bytes(b'{"ok":true}\n')
            pid = 43210
            pid_path.write_bytes((str(pid) + "\n").encode("ascii"))
            receipt_path.write_bytes(b"")
            for path in (
                sampler_source, csv_path, validation_path, pid_path,
                receipt_path,
            ):
                path.chmod(0o444)

            held = {
                "csv": os.open(csv_path, os.O_RDONLY | os.O_CLOEXEC),
                "validation": os.open(
                    validation_path, os.O_RDONLY | os.O_CLOEXEC
                ),
                "pid": os.open(pid_path, os.O_RDONLY | os.O_CLOEXEC),
                "receipt": os.open(receipt_path, os.O_RDONLY | os.O_CLOEXEC),
            }
            source_fd = os.open(sampler_source, os.O_RDONLY | os.O_CLOEXEC)
            try:
                csv_binding = subject.file_binding(held["csv"], 4096)
                validation_binding = subject.file_binding(
                    held["validation"], 4096
                )
                pid_binding = subject.file_binding(held["pid"], 4096)
                receipt_binding = subject.file_binding(held["receipt"], 4096)
                sampler_hash = hashlib.sha256(sampler_source.read_bytes()).hexdigest()
                prepared = prepared_fixture()
                prepared = subject.Prepared(
                    config=prepared.config, source_root=root,
                    controller=prepared.controller,
                    sampler_script=sampler_source, binary=prepared.binary,
                    launcher_sha256=prepared.launcher_sha256,
                    cgroup_path=prepared.cgroup_path,
                    cgroup_relative=prepared.cgroup_relative,
                    sampler_fd=source_fd, sampler_sha256=sampler_hash,
                    sampler_bytes=len(b"# sealed sampler fixture\n"),
                )
                with mock.patch.multiple(
                    subject, SAMPLER_CSV=csv_path,
                    SAMPLER_VALIDATION=validation_path,
                    SAMPLER_PID_FILE=pid_path, SAMPLER_RECEIPT=receipt_path,
                ):
                    expected_sampler_argv = subject.sampler_argv(
                        prepared, sampler_hash
                    )[1:]
                value = {
                    "argv": expected_sampler_argv,
                    "exit_code": 0,
                    "expected_output_owner_uid": subject.CAMPAIGN_UID,
                    "pid": pid,
                    "pid_file": {
                        "binding": {
                            "device": pid_binding["device"],
                            "inode": pid_binding["inode"],
                        },
                        "path": str(pid_path), "removed": True,
                    },
                    "raw_csv": {
                        "binding": {
                            "device": csv_binding["device"],
                            "inode": csv_binding["inode"],
                            "sha256": csv_binding["sha256"],
                            "size": csv_binding["bytes"],
                        },
                        "path": str(csv_path),
                    },
                    "receipt_file": {
                        "path": str(receipt_path),
                        "reservation_binding": {
                            "device": receipt_binding["device"],
                            "inode": receipt_binding["inode"],
                        },
                    },
                    "sampler_source": {
                        "expected_sha256": sampler_hash,
                        "path": str(sampler_source), "sha256": sampler_hash,
                    },
                    "schema": subject.SAMPLER_SCHEMA,
                    "self_sha256_excluding_field": None,
                    "validation_jsonl": {
                        "binding": {
                            "device": validation_binding["device"],
                            "inode": validation_binding["inode"],
                            "sha256": validation_binding["sha256"],
                            "size": validation_binding["bytes"],
                        },
                        "path": str(validation_path),
                    },
                }
                preimage = dict(value)
                del preimage["self_sha256_excluding_field"]
                value["self_sha256_excluding_field"] = hashlib.sha256(
                    canonical(preimage) + b"\n"
                ).hexdigest()
                raw = canonical(value) + b"\n"
                receipt_path.chmod(0o644)
                with receipt_path.open("r+b", buffering=0) as output:
                    output.write(raw)
                    output.flush()
                    os.fsync(output.fileno())
                receipt_path.chmod(0o444)
                pid_path.unlink()
                sampler_source.unlink()
                sampler_source.write_bytes(b"# substituted sampler\n")
                sampler_source.chmod(0o444)
                with mock.patch.multiple(
                    subject, SAMPLER_CSV=csv_path,
                    SAMPLER_VALIDATION=validation_path,
                    SAMPLER_PID_FILE=pid_path, SAMPLER_RECEIPT=receipt_path,
                ):
                    terminal = subject.validate_sampler_terminal(
                        raw, pid, prepared, held, 0,
                    )
                self.assertEqual(terminal["artifacts"]["pid"]["nlink"], 0)
                self.assertEqual(
                    terminal["artifacts"]["receipt"]["inode"],
                    receipt_binding["inode"],
                )
                self.assertFalse(terminal["source_path_still_bound"])
            finally:
                for fd in held.values():
                    os.close(fd)
                os.close(source_fd)


class AttemptJournalTests(unittest.TestCase):
    @unittest.skipUnless(
        os.geteuid() == 0,
        "real preexisting reservation gate requires root-owned scratch",
    )
    def test_preexisting_attempt_name_is_not_attributed_to_this_attempt(self):
        existing = Path(tempfile.mkdtemp(
            dir="/var/tmp", prefix=".wh2-attempt-preexisting-test-",
        ))
        callback = False

        def consumed_callback():
            nonlocal callback
            callback = True

        try:
            with mock.patch.object(subject, "LAUNCH_DIR", existing):
                with self.assertRaises(subject.LaunchError) as caught:
                    subject.AttemptJournal.reserve({}, consumed_callback)
            self.assertNotIsInstance(
                caught.exception, subject.AttemptConsumedError,
            )
            self.assertFalse(callback)
            self.assertTrue(existing.is_dir())
        finally:
            existing.rmdir()

    @unittest.skipUnless(
        os.geteuid() == 0,
        "real post-mkdir journal recovery requires root-owned scratch",
    )
    def test_each_post_mkdir_failure_retains_terminal_journal_authority(self):
        attempt = subject.self_hashed(subject.ATTEMPT_SCHEMA, {
            "attempt_id": "scratch", "build_dir": "/scratch/build",
            "expected_source_commit": COMMIT,
            "installed_launcher": "/scratch/launcher",
            "launcher_sha256": HASH, "service_cgroup": "/scratch",
            "started_monotonic_ns": 1,
            "started_utc": "2026-08-28T00:00:00.000000Z",
        })
        cases = (
            "mkdir-return-ambiguity", "post-return-trace", "callback",
            "directory-open", "directory-policy",
            "parent-fsync", "attempt-write",
        )
        for case in cases:
            with self.subTest(case=case):
                created = Path(tempfile.mkdtemp(
                    dir="/var/tmp", prefix=".wh2-attempt-journal-test-",
                ))
                os.rmdir(created)
                journal = None
                try:
                    real_open = subject.os.open
                    real_mkdir = subject.os.mkdir
                    real_fstat = subject.os.fstat
                    real_fsync = subject.os.fsync
                    open_failed = False
                    policy_calls = 0
                    fsync_failed = False

                    def injected_mkdir(path, mode=0o777, *args, **kwargs):
                        if (
                            case == "mkdir-return-ambiguity"
                            and path == created.name
                            and kwargs.get("dir_fd") is not None
                        ):
                            real_mkdir(path, mode, *args, **kwargs)
                            raise OSError("injected post-mkdir return ambiguity")
                        return real_mkdir(path, mode, *args, **kwargs)

                    def injected_open(path, flags, *args, **kwargs):
                        nonlocal open_failed
                        if (
                            case == "directory-open" and not open_failed
                            and path == created.name
                            and kwargs.get("dir_fd") is not None
                        ):
                            open_failed = True
                            raise OSError("injected attempt directory open")
                        return real_open(path, flags, *args, **kwargs)

                    def injected_fstat(fd):
                        nonlocal policy_calls
                        policy_calls += 1
                        if case == "directory-policy" and policy_calls == 2:
                            raise OSError("injected attempt directory policy")
                        return real_fstat(fd)

                    def injected_fsync(fd):
                        nonlocal fsync_failed
                        if case == "parent-fsync" and not fsync_failed:
                            fsync_failed = True
                            raise OSError("injected attempt parent fsync")
                        return real_fsync(fd)

                    original_write_json = subject.AttemptJournal.write_json

                    def injected_write_json(self, name, value):
                        if case == "attempt-write" and name == "ATTEMPT":
                            raise OSError("injected ATTEMPT write")
                        return original_write_json(self, name, value)

                    def consumed_callback():
                        if case == "callback":
                            raise OSError("injected consumed callback")

                    source_lines, source_start = inspect.getsourcelines(
                        subject.AttemptJournal.reserve.__func__
                    )
                    consumed_line = next(
                        source_start + index
                        for index, line in enumerate(source_lines)
                        if line.strip() == "consumed = True"
                    )

                    def trace_after_mkdir(frame, event, _argument):
                        if (
                            case == "post-return-trace" and event == "line"
                            and frame.f_code
                            is subject.AttemptJournal.reserve.__func__.__code__
                            and frame.f_lineno == consumed_line
                        ):
                            sys.settrace(None)
                            raise OSError("injected post-mkdir trace ambiguity")
                        return trace_after_mkdir

                    with mock.patch.object(subject, "LAUNCH_DIR", created), \
                         mock.patch.object(subject.os, "open", injected_open), \
                         mock.patch.object(subject.os, "mkdir", injected_mkdir), \
                         mock.patch.object(subject.os, "fstat", injected_fstat), \
                         mock.patch.object(subject.os, "fsync", injected_fsync), \
                         mock.patch.object(
                             subject.AttemptJournal, "write_json",
                             injected_write_json,
                         ):
                        if case == "post-return-trace":
                            sys.settrace(trace_after_mkdir)
                        try:
                            with self.assertRaises(
                                subject.AttemptConsumedError,
                            ) as caught:
                                subject.AttemptJournal.reserve(
                                    attempt, consumed_callback,
                                )
                        finally:
                            sys.settrace(None)
                        journal = caught.exception.journal
                        self.assertIsNotNone(journal)
                        self.assertGreaterEqual(journal.parent_fd, 0)
                        journal.write_json(
                            "terminal.json",
                            subject.self_hashed(subject.TERMINAL_SCHEMA, {
                                "errors": ["scratch reservation failure"],
                                "launcher_status": "fault",
                                "outcome": "consumed_unclassifiable",
                            }),
                        )
                        journal.finish()
                        self.assertTrue(journal.complete)
                        self.assertTrue((created / "COMPLETE").exists())
                finally:
                    if journal is not None:
                        journal.close()
                    if created.exists():
                        created.chmod(0o700)
                        shutil.rmtree(str(created))

    def test_sampler_directory_post_mkdir_trace_is_recovered_for_terminal(self):
        with tempfile.TemporaryDirectory() as temporary:
            parent = Path(temporary)
            sampler_dir = parent / ".sampler"
            parent_fd = os.open(parent, os.O_RDONLY | os.O_DIRECTORY)
            journal = types.SimpleNamespace(parent_fd=parent_fd)
            backend = subject.RealBackend()
            prepared = types.SimpleNamespace(
                config=subject.LaunchConfig(Path("/sealed/build"), COMMIT),
                launcher_sha256=HASH,
                cgroup_relative="/scratch",
            )

            def reserve_attempt(_attempt, consumed_callback):
                consumed_callback()
                return journal

            source_lines, source_start = inspect.getsourcelines(
                subject.RealBackend.reserve
            )
            created_line = next(
                source_start + index
                for index, line in enumerate(source_lines)
                if line.strip() == 'self.sampler_dir_state = "created"'
            )

            def trace_after_sampler_mkdir(frame, event, _argument):
                if (
                    event == "line"
                    and frame.f_code is subject.RealBackend.reserve.__code__
                    and frame.f_lineno == created_line
                ):
                    sys.settrace(None)
                    raise OSError("injected sampler post-mkdir trace ambiguity")
                return trace_after_sampler_mkdir

            try:
                with mock.patch.object(subject, "SAMPLER_DIR", sampler_dir), \
                     mock.patch.object(
                         subject.AttemptJournal, "reserve",
                         side_effect=reserve_attempt,
                     ):
                    sys.settrace(trace_after_sampler_mkdir)
                    try:
                        with self.assertRaisesRegex(
                            OSError, "sampler post-mkdir trace ambiguity",
                        ):
                            backend.reserve(prepared)
                    finally:
                        sys.settrace(None)
                    self.assertTrue(backend.attempt_consumed)
                    self.assertEqual(backend.sampler_dir_state, "bound")
                    self.assertGreaterEqual(backend.sampler_dir_fd, 0)
                    self.assertEqual(
                        backend.sampler_dir_binding,
                        (
                            sampler_dir.stat().st_dev,
                            sampler_dir.stat().st_ino,
                        ),
                    )
                    backend.finish_sampler_dir()
                    self.assertEqual(backend.sampler_dir_state, "closed")
                    info = sampler_dir.stat()
                    self.assertEqual(stat.S_IMODE(info.st_mode), 0o500)
                    self.assertEqual(
                        (info.st_uid, info.st_gid),
                        (subject.CAMPAIGN_UID, subject.CAMPAIGN_GID),
                    )
            finally:
                sys.settrace(None)
                if backend.sampler_dir_fd >= 0:
                    os.close(backend.sampler_dir_fd)
                    backend.sampler_dir_fd = -1
                if sampler_dir.exists():
                    sampler_dir.chmod(0o700)
                    sampler_dir.rmdir()
                os.close(parent_fd)

    @unittest.skipUnless(
        os.geteuid() == 0,
        "real reservation-to-supervisor recovery requires root scratch",
    )
    def test_mkdir_return_ambiguity_terminalizes_through_real_backend(self):
        launch = Path(tempfile.mkdtemp(
            dir="/var/tmp", prefix=".wh2-attempt-supervisor-test-",
        ))
        os.rmdir(launch)
        sampler = Path(str(launch) + ".sampler")
        backend = subject.RealBackend()
        prepared = prepared_fixture()
        real_mkdir = subject.os.mkdir
        injected = False

        def ambiguous_mkdir(path, mode=0o777, *args, **kwargs):
            nonlocal injected
            if (
                not injected and path == launch.name
                and kwargs.get("dir_fd") is not None
            ):
                injected = True
                real_mkdir(path, mode, *args, **kwargs)
                raise OSError("injected post-mkdir return ambiguity")
            return real_mkdir(path, mode, *args, **kwargs)

        try:
            with mock.patch.multiple(
                subject, LAUNCH_DIR=launch, SAMPLER_DIR=sampler,
            ), mock.patch.object(
                backend, "prepare", return_value=prepared,
            ), mock.patch.object(subject.os, "mkdir", ambiguous_mkdir):
                result = subject.Supervisor(backend).execute(prepared.config)
            self.assertTrue(injected)
            self.assertTrue(backend.attempt_consumed)
            self.assertIsNotNone(backend.journal)
            self.assertTrue(backend.journal.complete)
            self.assertEqual(result.classification, "consumed_unclassifiable")
            self.assertEqual(result.launcher_status, "fault")
            self.assertEqual(result.exit_code, 1)
            self.assertTrue((launch / "terminal.json").is_file())
            self.assertTrue((launch / "COMPLETE").is_file())
            self.assertFalse(sampler.exists())
        finally:
            if launch.exists():
                launch.chmod(0o700)
                shutil.rmtree(str(launch))
            if sampler.exists():
                sampler.chmod(0o700)
                shutil.rmtree(str(sampler))


class BuildProfileTests(unittest.TestCase):
    @staticmethod
    def validate(build, source_root):
        with mock.patch.object(
            subject, "compiler_image_receipt", return_value={"sealed": True},
        ), mock.patch.object(
            subject, "frozen_toolchain_receipts", return_value=[{"sealed": True}],
        ), mock.patch.object(
            subject, "ninja_target_commands",
            side_effect=lambda build: effective_ninja_commands_fixture(build),
        ):
            return subject.validate_build_profile(build, source_root, COMMIT)

    def test_exact_profile_returns_crossbound_manifest_and_receipt(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            build = root / "build"
            source_root = root / "source"
            build.mkdir()
            source_root.mkdir()
            write_build_profile_fixture(build, source_root)
            manifest_hash, profile = self.validate(build, source_root)
            self.assertRegex(manifest_hash, r"^[0-9a-f]{64}$")
            self.assertEqual(
                [item["path"] for item in profile["artifact_bindings"]],
                list(subject.BUILD_PATHS),
            )
            self.assertEqual(
                profile["receipt_sha256"],
                hashlib.sha256(canonical({
                    **profile, "receipt_sha256": None,
                })).hexdigest(),
            )

    def test_effective_target_commands_reject_edge_and_link_expansion(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            build = root / "build"
            source_root = root / "source"
            build.mkdir()
            source_root.mkdir()
            write_build_profile_fixture(build, source_root)
            commands = json.loads(
                (build / "compile_commands.json").read_text("utf-8")
            )
            compile_tokens = {
                item["output"]: shlex.split(item["command"], posix=True)
                for item in commands
            }
            outputs = [
                item["output"] for item in commands
                if item["output"].startswith(
                    "CMakeFiles/" + subject.BINARY_NAME + ".dir/bench/"
                )
            ]
            lines = effective_ninja_commands_fixture(build)
            subject.validate_effective_ninja_commands(
                lines, compile_tokens, outputs,
            )
            mutations = {
                "launcher": lambda value: value.__setitem__(
                    0, "/sealed/wrapper -- " + value[0]
                ),
                "test-hooks": lambda value: value.__setitem__(
                    0, value[0].replace(
                        " -MD ", " -DWIREHAIR_V2_ENABLE_TEST_HOOKS=1 -MD ", 1,
                    ),
                ),
                "link-flags": lambda value: value.__setitem__(
                    -1, value[-1].replace(" -o ", " -Wl,--wrap=malloc -o ", 1),
                ),
                "transitive-compile-omitted": lambda value: value.pop(4),
                "archive-omitted": lambda value: value.pop(-3),
                "archive-input": lambda value: value.__setitem__(
                    -3, value[-3].replace(
                        " CMakeFiles/wirehair.dir/wirehair.cpp.o ",
                        " CMakeFiles/wirehair.dir/forged.cpp.o ", 1,
                    ),
                ),
                "quoted-archive-operator": lambda value: value.__setitem__(
                    -3, value[-3].replace(" && ", " '&&' ", 1),
                ),
                "quoted-link-operator": lambda value: value.__setitem__(
                    -1, value[-1].replace(" && ", " '&&' ", 1),
                ),
                "unexpected": lambda value: value.append("/sealed/extra-command"),
            }
            for name, mutate in mutations.items():
                altered = list(lines)
                mutate(altered)
                with self.subTest(name=name), self.assertRaises(subject.LaunchError):
                    subject.validate_effective_ninja_commands(
                        altered, compile_tokens, outputs,
                    )

    def test_cache_command_and_generated_flag_mutations_fail_preconsumption(self):
        mutations = {
            "cache": ("CMakeCache.txt", b"MARCH_NATIVE:BOOL=OFF", b"MARCH_NATIVE:BOOL=ON"),
            "command": ("compile_commands.json", b"-O3", b"-O2"),
            "command-output": (
                "compile_commands.json", b"-o CMakeFiles/", b"-o forged/",
            ),
            "ninja": ("build.ninja", b"\n", b"\n# -flto\n"),
            "ninja-o4": ("build.ninja", b"  FLAGS = -O3", b"  FLAGS = -O4"),
            "duplicate-pre-link": (
                "build.ninja", b"  PRE_LINK = :\n",
                b"  PRE_LINK = :\n  PRE_LINK = /evil\n",
            ),
            "compile-edge-launcher": (
                "build.ninja", b"  DEP_FILE = CMakeFiles/",
                b"  LAUNCHER = /sealed/wrapper --\n  DEP_FILE = CMakeFiles/",
            ),
            "compile-edge-test-hooks": (
                "build.ninja",
                b"  DEFINES = -DWIREHAIR_STATIC=1",
                b"  DEFINES = -DWIREHAIR_V2_ENABLE_TEST_HOOKS=1",
            ),
            "global-link-flags": (
                "build.ninja", ("build " + subject.BINARY_NAME + ":").encode(),
                b"LINK_FLAGS = -Wl,--wrap=malloc\nbuild "
                + subject.BINARY_NAME.encode() + b":",
            ),
            "rule-suffix": (
                "CMakeFiles/rules.ninja", b" -o $out -c $in\n",
                b" -o $out -c $in && /evil\n",
            ),
        }
        for name, (relative, old, new) in mutations.items():
            with self.subTest(name=name), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                build = root / "build"
                source_root = root / "source"
                build.mkdir()
                source_root.mkdir()
                write_build_profile_fixture(build, source_root)
                path = build / relative
                raw = path.read_bytes()
                self.assertIn(old, raw)
                path.chmod(0o644)
                path.write_bytes(raw.replace(old, new, 1))
                path.chmod(0o444)
                with self.assertRaises(subject.LaunchError):
                    self.validate(build, source_root)

    def test_held_fd_snapshot_detects_in_place_post_parse_change(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            build = root / "build"
            source_root = root / "source"
            build.mkdir()
            source_root.mkdir()
            write_build_profile_fixture(build, source_root)
            changed = False

            def mutate_after_parse(_path, _language):
                nonlocal changed
                if not changed:
                    changed = True
                    target = build / "build.ninja"
                    target.chmod(0o644)
                    target.write_bytes(target.read_bytes() + b"# changed\n")
                    target.chmod(0o444)
                return {"sealed": True}

            with mock.patch.object(
                subject, "compiler_image_receipt", side_effect=mutate_after_parse,
            ), mock.patch.object(
                subject, "frozen_toolchain_receipts", return_value=[],
            ), mock.patch.object(
                subject, "ninja_target_commands",
                side_effect=lambda value: effective_ninja_commands_fixture(value),
            ), self.assertRaisesRegex(subject.LaunchError, "changed during"):
                subject.validate_build_profile(build, source_root, COMMIT)

    def test_transitive_compile_semantics_are_exact(self):
        mutations = {
            "core-building-define-removed": lambda command: command.replace(
                " -DWIREHAIR_BUILDING=1", "", 1,
            ),
            "core-unexpected-define": lambda command: command.replace(
                " -DWIREHAIR_STATIC=1",
                " -DWH2_UNEXPECTED=1 -DWIREHAIR_STATIC=1", 1,
            ),
            "timing-benchmark-equations-disabled": lambda command: (
                command.replace(
                    "-DWIREHAIR_V2_ENABLE_BENCHMARK_EQUATIONS=1",
                    "-DWIREHAIR_V2_ENABLE_BENCHMARK_EQUATIONS=0", 1,
                )
            ),
        }
        for name, mutate in mutations.items():
            with self.subTest(name=name), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                build = root / "build"
                source_root = root / "source"
                build.mkdir()
                source_root.mkdir()
                write_build_profile_fixture(build, source_root)
                path = build / "compile_commands.json"
                commands = json.loads(path.read_text("utf-8"))
                prefix = (
                    "codec/CMakeFiles/wirehair_v2_timing_policy.dir/"
                    if name.startswith("timing-")
                    else "CMakeFiles/wirehair.dir/"
                )
                target = next(
                    item for item in commands
                    if item["output"].startswith(prefix)
                )
                changed = mutate(target["command"])
                self.assertNotEqual(changed, target["command"])
                target["command"] = changed
                path.chmod(0o644)
                path.write_text(json.dumps(commands), encoding="utf-8")
                path.chmod(0o444)
                with self.assertRaisesRegex(
                    subject.LaunchError, "transitive compile command",
                ):
                    self.validate(build, source_root)


class BuildAuthorityTests(unittest.TestCase):
    def validate(self, value, config, observations):
        source_info = types.SimpleNamespace(
            st_uid=subject.CAMPAIGN_UID, st_gid=subject.CAMPAIGN_GID,
            st_nlink=1, st_mode=stat.S_IFREG | 0o444,
        )
        with mock.patch.multiple(
            subject,
            validate_build_staging_paths=mock.DEFAULT,
            current_build_subdirectories=mock.DEFAULT,
            installed_launcher_binding=mock.DEFAULT,
            hash_path=mock.DEFAULT,
            frozen_toolchain_receipts=mock.DEFAULT,
            static_git_gate=mock.DEFAULT,
            static_source_manifest=mock.DEFAULT,
            validate_build_profile=mock.DEFAULT,
            current_binary_receipt=mock.DEFAULT,
        ) as patched:
            patched["validate_build_staging_paths"].return_value = {
                "build": observations["directory"]
            }
            patched["current_build_subdirectories"].return_value = (
                observations["subdirectories"]
            )
            patched["installed_launcher_binding"].return_value = (
                observations["installed"]
            )
            patched["hash_path"].return_value = (HASH, source_info)
            patched["frozen_toolchain_receipts"].return_value = (
                observations["tools"]
            )
            patched["static_git_gate"].return_value = HASH
            patched["static_source_manifest"].return_value = (HASH, {})
            patched["validate_build_profile"].return_value = (
                HASH, observations["profile"]
            )
            patched["current_binary_receipt"].return_value = (
                observations["binary"]
            )
            return subject.validate_build_authority_receipt(value, config)

    @staticmethod
    def reseal(value):
        value["receipt_sha256"] = None
        value["receipt_sha256"] = hashlib.sha256(canonical(value)).hexdigest()

    def test_exact_build_authority_crossbinds_every_future_input(self):
        config = subject.BuildSealConfig(
            Path("/sealed/source"), Path("/sealed/build"), COMMIT,
        )
        value, observations = build_authority_fixture(config)
        self.assertEqual(
            self.validate(value, config, observations), value["receipt_sha256"],
        )
        mutations = {
            "commit": lambda item: item.__setitem__(
                "expected_source_commit", "2" * 40
            ),
            "source": lambda item: item.__setitem__(
                "source_root", "/sealed/other-source"
            ),
            "toolchain": lambda item: item.__setitem__(
                "toolchain_after", [{"sealed": "other"}]
            ),
            "binary": lambda item: item.__setitem__(
                "binary", {"sealed": "other"}
            ),
            "directory": lambda item: item.__setitem__(
                "build_directory", {"sealed": "other"}
            ),
            "systemd": lambda item: item["systemd_run_argv"].__setitem__(
                -1, "2" * 40
            ),
        }
        for name, mutation in mutations.items():
            with self.subTest(name=name):
                changed = copy.deepcopy(value)
                mutation(changed)
                self.reseal(changed)
                with self.assertRaises(subject.LaunchError):
                    self.validate(changed, config, observations)

    def test_command_reseal_cannot_hide_overlap_or_no_work_tamper(self):
        config = subject.BuildSealConfig(
            Path("/sealed/source"), Path("/sealed/build"), COMMIT,
        )
        value, observations = build_authority_fixture(config)
        overlap = copy.deepcopy(value)
        overlap["commands"][1]["started_monotonic_ns"] = 1
        self.reseal(overlap["commands"][1])
        self.reseal(overlap)
        with self.assertRaises(subject.LaunchError):
            self.validate(overlap, config, observations)
        no_work = copy.deepcopy(value)
        command = no_work["commands"][-1]
        command["stdout_utf8"] += "forged\n"
        raw = command["stdout_utf8"].encode("utf-8")
        command["stdout_bytes"] = len(raw)
        command["stdout_sha256"] = hashlib.sha256(raw).hexdigest()
        self.reseal(command)
        self.reseal(no_work)
        with self.assertRaises(subject.LaunchError):
            self.validate(no_work, config, observations)

        utc = copy.deepcopy(value)
        command = utc["commands"][0]
        command["finished_utc"] = "2026-99-99T00:00:00.000000Z"
        self.reseal(command)
        self.reseal(utc)
        with self.assertRaises(subject.LaunchError):
            self.validate(utc, config, observations)

        security = copy.deepcopy(value)
        command = security["commands"][0]
        command["child_security"]["uid"] = [0, 0, 0]
        self.reseal(command["child_security"])
        self.reseal(command)
        self.reseal(security)
        with self.assertRaises(subject.LaunchError):
            self.validate(security, config, observations)

        child_pid = copy.deepcopy(value)
        command = child_pid["commands"][0]
        command["child_pid"] += 1
        self.reseal(command)
        self.reseal(child_pid)
        with self.assertRaises(subject.LaunchError):
            self.validate(child_pid, config, observations)

        boolean_integer = copy.deepcopy(value)
        boolean_integer["root_boundary"]["uid"] = False
        self.reseal(boolean_integer)
        with self.assertRaises(subject.LaunchError):
            self.validate(boolean_integer, config, observations)

        boolean_returncode = copy.deepcopy(value)
        command = boolean_returncode["commands"][0]
        command["returncode"] = False
        self.reseal(command)
        self.reseal(boolean_returncode)
        with self.assertRaises(subject.LaunchError):
            self.validate(boolean_returncode, config, observations)

    @unittest.skipUnless(
        os.geteuid() == 0 and {0, 1}.issubset(os.sched_getaffinity(0)),
        "exact build-child credential transition requires root and CPUs 0-1",
    )
    def test_build_child_security_handshake_is_observed_before_exec(self):
        previous_affinity = os.sched_getaffinity(0)
        previous_groups = os.getgroups()
        try:
            os.setgroups([])
            os.sched_setaffinity(0, {0, 1})
            with tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                logs = root / "logs"
                logs.mkdir()
                log_fd = os.open(
                    logs, os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC,
                )
                try:
                    receipt = subject.run_build_command(
                        "probe", ["/usr/bin/true"], 5.0,
                        subject.build_environment(root), root, log_fd,
                    )
                finally:
                    os.close(log_fd)
        finally:
            os.sched_setaffinity(0, previous_affinity)
            os.setgroups(previous_groups)
        self.assertEqual(receipt["returncode"], 0)
        self.assertEqual(receipt["child_pid"], receipt["child_security"]["pid"])
        self.assertEqual(receipt["child_security"]["affinity"], [0, 1])
        self.assertEqual(
            receipt["child_security"]["uid"], [subject.CAMPAIGN_UID] * 3,
        )
        self.assertEqual(
            receipt["child_security"]["gid"], [subject.CAMPAIGN_GID] * 3,
        )
        self.assertEqual(receipt["child_security"]["groups"], [])
        self.assertEqual(receipt["child_security"]["no_new_privs"], 1)

    def test_build_namespace_gate_checks_host_view_behind_private_tmp(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            host_root = root / "host"
            host_fixed = host_root / "var/tmp/fixed-output"
            host_fixed.parent.mkdir(parents=True)
            host_fixed.mkdir()
            fixed = Path("/var/tmp/fixed-output")
            missing = tuple(
                Path("/var/tmp/missing-" + str(index)) for index in range(3)
            )
            with mock.patch.multiple(
                subject, HOST_MOUNT_ROOT=host_root, LAUNCH_DIR=fixed,
                SAMPLER_DIR=missing[0], FIXED_CLAIM_PATH=missing[1],
                FIXED_OUTPUT_DIR=missing[2],
            ), self.assertRaisesRegex(subject.LaunchError, "host view"):
                subject.require_scientific_namespaces_absent()

    @unittest.skipUnless(
        os.geteuid() == 0, "root-owned no-replace publication requires root",
    )
    def test_build_authority_publication_is_atomic_no_replace(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / "authority.json"
            first = canonical({"receipt": "first"}) + b"\n"
            second = canonical({"receipt": "second"}) + b"\n"
            with mock.patch.multiple(
                subject, BUILD_AUTHORITY_DIR=root,
                BUILD_AUTHORITY_PATH=target,
            ):
                subject.publish_build_authority(first)
                with self.assertRaises(FileExistsError):
                    subject.publish_build_authority(second)
            info = target.stat()
            self.assertEqual(target.read_bytes(), first)
            self.assertEqual(info.st_uid, 0)
            self.assertEqual(info.st_gid, 0)
            self.assertEqual(stat.S_IMODE(info.st_mode), 0o444)
            self.assertEqual(info.st_nlink, 1)

    def test_build_descendant_observation_is_killed_and_fails(self):
        waits = [(43210, 9), ChildProcessError(), ChildProcessError()]
        with mock.patch.object(
            Path, "read_bytes", side_effect=[b"43210\n", b""]
        ), mock.patch.object(
            subject.os, "pidfd_open", return_value=77, create=True,
        ), mock.patch.object(
            subject.signal, "pidfd_send_signal", create=True,
        ) as send, mock.patch.object(
            subject.os, "waitpid", side_effect=waits,
        ), mock.patch.object(subject.os, "close"):
            with self.assertRaisesRegex(subject.LaunchError, "left descendants"):
                subject.require_no_build_descendants("fixture")
        send.assert_called_once_with(77, subject.signal.SIGKILL)


class StateMachineTests(unittest.TestCase):
    def run_subject(self, backend):
        return subject.Supervisor(backend).execute(
            subject.LaunchConfig(Path("/sealed/build"), COMMIT)
        )

    def test_success_exact_order_and_roster(self):
        backend = FakeBackend()
        result = self.run_subject(backend)
        self.assertEqual(result.classification, "pass")
        self.assertEqual(result.launcher_status, "clean")
        self.assertEqual(result.exit_code, 0)
        self.assertEqual(backend.journal.names, FakeJournal.SUCCESS_ROSTER)
        self.assertTrue(backend.journal.finished)
        self.assertTrue(backend.guardian.cancelled)
        self.assertTrue(backend.guardian.finished)
        self.assertTrue(backend.guardian.closed)
        self.assertLess(
            backend.calls.index("run_verifier"),
            backend.calls.index("stop_sampler"),
        )
        self.assertEqual(backend.calls[-1], "close_authority")
        start = json.loads(backend.journal.payloads["start.json"])
        deadline = json.loads(backend.journal.payloads["deadline.json"])
        anchor = start["controller_release_not_before_monotonic_ns"]
        self.assertEqual(deadline["guardian_anchor_monotonic_ns"], anchor)
        self.assertEqual(
            deadline["admission_deadline_monotonic_ns"],
            anchor + int(subject.CLAIM_ADMISSION_SECONDS * 1e9),
        )
        self.assertGreaterEqual(
            deadline["controller_release_monotonic_ns"], anchor,
        )

    def test_only_verifier_classifies_regardless_of_controller_exit(self):
        for classification, expected in subject.VERIFY_EXIT.items():
            for controller_rc in (0, 1, 2, -9):
                backend = FakeBackend(classification, controller_rc)
                result = self.run_subject(backend)
                self.assertEqual(result.classification, classification)
                self.assertEqual(result.exit_code, expected)

    def test_start_journal_delay_does_not_reanchor_admission_deadline(self):
        backend = FakeBackend(journal_delay="start.json")
        self.run_subject(backend)
        start = json.loads(backend.journal.payloads["start.json"])
        deadline = json.loads(backend.journal.payloads["deadline.json"])
        anchor = start["controller_release_not_before_monotonic_ns"]
        released = deadline["controller_release_monotonic_ns"]
        self.assertGreaterEqual(released - anchor, 5_000_000)
        self.assertEqual(
            deadline["admission_deadline_monotonic_ns"],
            anchor + int(subject.CLAIM_ADMISSION_SECONDS * 1e9),
        )

    def test_admission_edge_records_observed_claim_but_not_exact_outer_arm(self):
        backend = FakeBackend()
        backend.guardian.exact_arm_result = False
        result = self.run_subject(backend)
        deadline = json.loads(backend.journal.payloads["deadline.json"])
        anchor = deadline["guardian_anchor_monotonic_ns"]
        self.assertIsNotNone(deadline["controller_started_monotonic_ns"])
        self.assertFalse(deadline["exact_deadline_armed"])
        self.assertTrue(deadline["exact_backstop_reanchored"])
        self.assertEqual(
            deadline["experiment_deadline_monotonic_ns"],
            anchor + int(subject.CLAIM_ADMISSION_SECONDS * 1e9),
        )
        self.assertEqual(result.classification, "pass")
        self.assertEqual(result.launcher_status, "fault")
        self.assertEqual(result.exit_code, 1)

    def test_each_post_reservation_failure_is_terminal_and_never_retried(self):
        steps = (
            "start_sampler", "wait_sampler_ready", "run_preflight",
            "spawn_stubs", "spawn_guardian", "release_controller",
            "stop_sampler", "run_verifier",
        )
        for step in steps:
            with self.subTest(step=step):
                backend = FakeBackend(fail_at=step)
                result = self.run_subject(backend)
                if step in (
                    "start_sampler", "wait_sampler_ready", "run_preflight",
                    "spawn_stubs", "run_verifier",
                ):
                    self.assertEqual(
                        result.classification, "consumed_unclassifiable"
                    )
                self.assertEqual(backend.calls.count("reserve"), 1)
                self.assertTrue(backend.journal.finished)
                self.assertEqual(backend.journal.names[-2:], ["terminal.json", "COMPLETE"])
                self.assertTrue(backend.cleaned)

    def test_read_only_failure_does_not_consume(self):
        backend = FakeBackend(fail_at="prepare")
        with self.assertRaises(subject.LaunchError):
            self.run_subject(backend)
        self.assertNotIn("reserve", backend.calls)
        self.assertEqual(backend.journal.names, [])
        self.assertTrue(backend.cleaned)

    def test_nonempty_experiment_after_reap_fails_closed(self):
        backend = FakeBackend()
        backend.experiment_roster = [999]
        result = self.run_subject(backend)
        self.assertEqual(result.classification, "pass")
        self.assertIn("run_verifier", backend.calls)
        terminal = json.loads(backend.journal.payloads["terminal.json"])
        self.assertEqual(terminal["launcher_status"], "fault")
        self.assertTrue(terminal["errors"])
        self.assertEqual(result.launcher_status, "fault")
        self.assertEqual(result.exit_code, 1)

    def test_sampler_terminal_failure_cannot_erase_retained_classification(self):
        backend = FakeBackend(classification="reject")
        backend.sampler_stop_fails = True
        with self.assertRaisesRegex(
            subject.LaunchError, "completion is incomplete",
        ):
            self.run_subject(backend)
        self.assertLess(
            backend.calls.index("run_verifier"),
            backend.calls.index("stop_sampler"),
        )
        self.assertIn("verification.json", backend.journal.names)
        self.assertNotIn("sampler-terminal.json", backend.journal.names)
        self.assertNotIn("COMPLETE", backend.journal.names)

    def test_verifier_failure_is_contained_before_sampler_terminalization(self):
        backend = FakeBackend(fail_at="run_verifier")
        result = self.run_subject(backend)
        self.assertEqual(result.classification, "consumed_unclassifiable")
        self.assertLess(
            backend.calls.index("run_verifier"),
            backend.calls.index("contain_verifier"),
        )
        self.assertLess(
            backend.calls.index("contain_verifier"),
            backend.calls.index("stop_sampler"),
        )
        self.assertFalse(backend.verifier_live)
        self.assertIn("sampler-terminal.json", backend.journal.names)

    def test_uncontained_verifier_preserves_sampler_and_backstop(self):
        backend = FakeBackend(fail_at="run_verifier")
        backend.verifier_containment_fails = True
        backend.quiesce_clears_verifier = True
        with self.assertRaisesRegex(
            subject.LaunchError, "completion is incomplete",
        ):
            self.run_subject(backend)
        self.assertFalse(backend.verifier_live)
        self.assertTrue(backend.cleaned)
        self.assertNotIn("stop_sampler", backend.calls)
        self.assertFalse(backend.guardian.finished)
        self.assertFalse(backend.runtime_closed)
        self.assertNotIn("sampler-terminal.json", backend.journal.names)
        self.assertNotIn("COMPLETE", backend.journal.names)

    def test_journal_failure_never_retries_name_or_attempt(self):
        backend = FakeBackend(journal_fail="authority.json")
        result = self.run_subject(backend)
        self.assertEqual(result.classification, "consumed_unclassifiable")
        self.assertEqual(backend.calls.count("reserve"), 1)
        self.assertNotIn("authority.json", backend.journal.names)
        self.assertEqual(backend.journal.names[-2:], ["terminal.json", "COMPLETE"])

    def test_post_mkdir_reservation_failure_remains_consumed_and_terminal(self):
        backend = FakeBackend(fail_at="reserve_after_consumed")
        result = self.run_subject(backend)
        self.assertEqual(result.classification, "consumed_unclassifiable")
        self.assertEqual(backend.calls.count("reserve"), 1)
        self.assertEqual(
            backend.journal.names[-2:], ["terminal.json", "COMPLETE"]
        )
        self.assertTrue(backend.journal.closed)

    def test_real_backend_can_close_an_absent_post_attempt_sampler_dir(self):
        with tempfile.TemporaryDirectory() as temporary:
            parent = Path(temporary)
            parent_fd = os.open(parent, os.O_RDONLY | os.O_DIRECTORY)
            backend = subject.RealBackend()
            backend.journal = types.SimpleNamespace(parent_fd=parent_fd)
            backend.sampler_dir_state = "absent"
            try:
                with mock.patch.object(
                    subject, "SAMPLER_DIR", parent / ".sampler",
                ):
                    self.assertIsNone(backend.finish_sampler_dir())
                    (parent / ".sampler").mkdir()
                    with self.assertRaisesRegex(
                        subject.LaunchError, "unexpected name",
                    ):
                        backend.finish_sampler_dir()
            finally:
                os.close(parent_fd)

    def test_complete_failure_never_returns_or_reports_complete(self):
        backend = FakeBackend(journal_fail="COMPLETE")
        supervisor = subject.Supervisor(backend)
        with self.assertRaisesRegex(subject.LaunchError, "completion is incomplete"):
            supervisor.execute(subject.LaunchConfig(Path("/sealed/build"), COMMIT))
        self.assertNotIn("COMPLETE", backend.journal.names)
        self.assertNotIn("complete", supervisor.phases)

    def test_live_run_descendant_preserves_backstop_and_forbids_complete(self):
        backend = FakeBackend()
        backend.quiesce_ok = False
        with self.assertRaisesRegex(subject.LaunchError, "completion is incomplete"):
            self.run_subject(backend)
        self.assertFalse(backend.guardian.finished)
        self.assertNotIn("COMPLETE", backend.journal.names)
        self.assertFalse(backend.runtime_closed)

    def test_guardian_channel_close_failure_forbids_complete(self):
        backend = FakeBackend()
        backend.guardian.close_fails = True
        with self.assertRaisesRegex(subject.LaunchError, "completion is incomplete"):
            self.run_subject(backend)
        self.assertNotIn("COMPLETE", backend.journal.names)

    def test_guardian_cancel_or_finish_failure_forbids_complete(self):
        for member in ("cancel_fails", "finish_fails"):
            with self.subTest(member=member):
                backend = FakeBackend()
                setattr(backend.guardian, member, True)
                with self.assertRaisesRegex(
                    subject.LaunchError, "completion is incomplete",
                ):
                    self.run_subject(backend)
                self.assertFalse(backend.guardian.finished)
                self.assertFalse(backend.runtime_closed)
                self.assertFalse(backend.journal.complete)
                self.assertNotIn("terminal.json", backend.journal.names)
                self.assertNotIn("COMPLETE", backend.journal.names)

    def test_runtime_or_authority_close_failure_forbids_complete(self):
        for step in ("finish_sampler_dir", "close_runtime", "close_authority"):
            with self.subTest(step=step):
                backend = FakeBackend(fail_at=step)
                with self.assertRaisesRegex(
                    subject.LaunchError, "completion is incomplete"
                ):
                    self.run_subject(backend)
                self.assertNotIn("COMPLETE", backend.journal.names)

    def test_no_fork_occurs_in_supervisor_after_guard_arm(self):
        source = inspect.getsource(subject.Supervisor.execute)
        suffix = source[source.index("arm_admission"):]
        self.assertNotIn("fork", suffix)
        self.assertLess(source.index("spawn_guardian"), source.index("arm_admission"))
        self.assertLess(source.index("run_verifier"), source.index("stop_sampler"))


class FailedSpawnCleanupTests(unittest.TestCase):
    @unittest.skipUnless(
        hasattr(os, "fork") and hasattr(os, "pidfd_open"),
        "real failed-spawn test requires Python pidfd bindings",
    )
    def test_spawn_exec_pidfd_open_failure_reaps_pre_release_child(self):
        backend = subject.RealBackend()
        before = len(os.listdir("/proc/self/fd"))
        real_reap = subject.reap_failed_spawn
        with mock.patch.object(
            subject.CgroupTree, "move_pid", return_value=None,
        ), mock.patch.object(
            subject.os, "pidfd_open", side_effect=OSError("injected pidfd_open"),
        ), mock.patch.object(
            subject, "reap_failed_spawn", wraps=real_reap,
        ) as reaper:
            with self.assertRaisesRegex(OSError, "injected pidfd_open"):
                backend._spawn_exec(
                    "preflight", Path("/scratch/fake-cgroup"), (120,),
                    ["/not/executed"],
                )
        self.assertEqual(len(reaper.call_args_list), 1)
        failed_pid = reaper.call_args.args[0]
        with self.assertRaises(ChildProcessError):
            os.waitpid(failed_pid, os.WNOHANG)
        self.assertEqual(backend.children, [])
        self.assertEqual(len(os.listdir("/proc/self/fd")), before)

    def test_every_real_spawn_failure_path_uses_common_reaper(self):
        for method in (
            subject.RealBackend._spawn_exec,
            subject.RealBackend._spawn_stub,
            subject.RealBackend.spawn_guardian,
        ):
            source = inspect.getsource(method)
            self.assertIn("reap_failed_spawn", source)
            self.assertNotIn("pidfd_signal(pidfd", source)


class GuardianProcessTests(unittest.TestCase):
    @unittest.skipUnless(
        hasattr(os, "pidfd_open") and hasattr(subject.signal, "pidfd_send_signal"),
        "scratch guardian test requires Python pidfd bindings",
    )
    def test_independent_admission_guard_fires_while_launcher_is_inactive(self):
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary) / "run"
            experiment = run / "experiment"
            experiment.mkdir(parents=True)
            (run / "cgroup.procs").write_text("", encoding="ascii")
            (experiment / "cgroup.procs").write_text("", encoding="ascii")
            journal_fd = os.open(run, os.O_RDONLY | os.O_DIRECTORY)
            experiment_read, experiment_write = os.pipe2(os.O_CLOEXEC)
            run_read, run_write = os.pipe2(os.O_CLOEXEC)
            parent_sock, child_sock = subject.socket.socketpair(
                subject.socket.AF_UNIX,
                subject.socket.SOCK_SEQPACKET | subject.socket.SOCK_CLOEXEC,
            )
            process = None
            handle = None
            with mock.patch.object(subject, "CLAIM_ADMISSION_SECONDS", 0.04), \
                    mock.patch.object(subject, "WHOLE_RUN_SECONDS", 0.5):
                pid = os.fork()
                if pid == 0:
                    try:
                        parent_sock.close()
                        subject.guardian_child(
                            child_sock, experiment_write, run_write,
                            experiment, run, journal_fd,
                        )
                    except BaseException:
                        os._exit(125)
                child_sock.close()
                process = subject.ProcessHandle(
                    "guardian", pid, os.pidfd_open(pid, 0), 1,
                )
                handle = subject.GuardianHandle(process, parent_sock)
                try:
                    handle.initial_ack()
                    anchor = subject.time.monotonic_ns()
                    handle.arm_admission(anchor)
                    # No supervisor polling or Python thread drives this wait.
                    subject.time.sleep(0.09)
                    os.set_blocking(experiment_read, False)
                    self.assertEqual(os.read(experiment_read, 1), b"1")
                    self.assertTrue((run / "deadline-fired.json").is_file())
                    self.assertFalse(handle.arm_exact(anchor + 1))
                    handle.cancel()
                    self.assertTrue(handle.fired)
                    handle.finish()
                    self.assertTrue(handle.finished)
                finally:
                    if handle is not None:
                        handle.close()
                    if process is not None and process.returncode is None:
                        try:
                            subject.pidfd_signal(process.pidfd, subject.signal.SIGKILL)
                        except (OSError, subject.LaunchError):
                            pass
                        try:
                            os.waitpid(process.pid, 0)
                        except ChildProcessError:
                            pass
                    if process is not None:
                        os.close(process.pidfd)
            for fd in (
                journal_fd, experiment_read, experiment_write, run_read,
                run_write,
            ):
                try:
                    os.close(fd)
                except OSError:
                    pass

    @unittest.skipUnless(
        hasattr(os, "pidfd_open") and hasattr(subject.signal, "pidfd_send_signal"),
        "scratch guardian test requires Python pidfd bindings",
    )
    def test_outer_kill_failure_preserves_independent_backstop(self):
        with tempfile.TemporaryDirectory() as temporary:
            run = Path(temporary) / "run"
            experiment = run / "experiment"
            experiment.mkdir(parents=True)
            (run / "cgroup.procs").write_text("", encoding="ascii")
            (experiment / "cgroup.procs").write_text("", encoding="ascii")
            journal_fd = os.open(run, os.O_RDONLY | os.O_DIRECTORY)
            experiment_read, experiment_write = os.pipe2(os.O_CLOEXEC)
            run_read, run_write = os.pipe2(os.O_CLOEXEC)
            parent_sock, child_sock = subject.socket.socketpair(
                subject.socket.AF_UNIX,
                subject.socket.SOCK_SEQPACKET | subject.socket.SOCK_CLOEXEC,
            )
            process = None
            handle = None
            with mock.patch.object(subject, "CLAIM_ADMISSION_SECONDS", 0.03), \
                    mock.patch.object(subject, "WHOLE_RUN_SECONDS", 0.05):
                pid = os.fork()
                if pid == 0:
                    try:
                        parent_sock.close()
                        subject.guardian_child(
                            child_sock, experiment_write, run_write,
                            experiment, run, journal_fd,
                        )
                    except BaseException:
                        os._exit(125)
                child_sock.close()
                process = subject.ProcessHandle(
                    "guardian", pid, os.pidfd_open(pid, 0), 1,
                )
                handle = subject.GuardianHandle(process, parent_sock)
                try:
                    handle.initial_ack()
                    os.close(experiment_read)
                    experiment_read = -1
                    handle.arm_admission(subject.time.monotonic_ns())
                    subject.time.sleep(0.14)
                    os.set_blocking(run_read, False)
                    self.assertEqual(os.read(run_read, 1), b"1")
                    self.assertTrue((run / "deadline-fired.json").is_file())
                    self.assertTrue((run / "backstop-fired.json").is_file())
                    waited, status = os.waitpid(pid, 0)
                    self.assertEqual(waited, pid)
                    process.returncode = subject.wait_status_code(status)
                    self.assertEqual(process.returncode, 125)
                finally:
                    if handle is not None:
                        handle.close()
                    if process is not None and process.returncode is None:
                        try:
                            subject.pidfd_signal(
                                process.pidfd, subject.signal.SIGKILL,
                            )
                        except (OSError, subject.LaunchError):
                            pass
                        try:
                            os.waitpid(process.pid, 0)
                        except ChildProcessError:
                            pass
                    if process is not None:
                        os.close(process.pidfd)
            for fd in (
                journal_fd, experiment_read, experiment_write, run_read,
                run_write,
            ):
                if fd >= 0:
                    try:
                        os.close(fd)
                    except OSError:
                        pass


class EntrypointTests(unittest.TestCase):
    def test_result_stdout_separates_science_from_launcher_fault(self):
        result = subject.RunResult("pass", "fault", 1, ["complete"])
        with mock.patch.object(subject, "RealBackend", return_value=object()), \
                mock.patch.object(
                    subject.Supervisor, "execute", return_value=result,
                ), mock.patch("builtins.print") as output:
            rc = subject.main([
                "--execute", "--build-dir", "/sealed/build",
                "--expected-source-commit", COMMIT,
            ])
        self.assertEqual(rc, 1)
        value = json.loads(output.call_args.args[0])
        self.assertEqual(value["outcome"], "pass")
        self.assertEqual(value["launcher_status"], "fault")
        self.assertEqual(value["status"], "complete")

    def test_selftest_subprocess_is_unprivileged_and_bounded(self):
        result = subprocess.run(
            [sys.executable, "-I", "-B", str(Path(subject.__file__).resolve()),
             "--selftest"],
            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, timeout=5.0, check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr.decode("utf-8"))
        self.assertEqual(
            result.stdout,
            b"WH2 direct systematic complement launcher selftest passed\n",
        )
        self.assertEqual(result.stderr, b"")


if __name__ == "__main__":
    unittest.main()
