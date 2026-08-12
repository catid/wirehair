#!/usr/bin/env python3
"""Deterministic tests for the bounded thermal transition controller."""

from __future__ import annotations

import ast
from dataclasses import replace
from contextlib import ExitStack, contextmanager, redirect_stderr
import copy
import csv
import errno
import hashlib
import inspect
import io
import json
import os
from pathlib import Path
import shutil
import signal
import stat
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

import wh2_thermal_sampler_transition as transition
import wh2_p32_dispatch_timing as frozen_p32


ROOT_STAGE_TEST_PREFIX = (
    transition.SOURCE_TEST_TRANSITION_PREFIX + "root-stage-")


def frozen_legacy_thermal_bytes(*, cpu_busy_pct: str = "50.0") -> bytes:
    """One canonical sample accepted by the unchanged frozen P32 parser."""
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(
        stream, fieldnames=frozen_p32.THERMAL_FIELDS, lineterminator="\n")
    writer.writeheader()
    row = {field: "1" for field in frozen_p32.THERMAL_FIELDS}
    row.update({
        "utc": "2026-08-12T00:00:00.000Z",
        "monotonic_s": "123.000000",
        "cpu_busy_pct": cpu_busy_pct,
        "cpu_avg_mhz": "4500.0", "cpu_tctl_c": "60.0",
        "dimm_read_errors": "0", "load1": "1.0", "load5": "1.0",
        "load15": "1.0", "edac_ce": "0", "edac_ue": "0",
    })
    for field in frozen_p32.DIMM_FIELDS:
        row[field] = "45.0"
    writer.writerow(row)
    return stream.getvalue().encode("ascii")


def remove_root_stage_test_tree(path: Path) -> None:
    """Remove only a uniquely named synthetic v6 stage-test tree."""
    path = Path(path)
    if path.parent != Path("/dev/shm") or \
            not path.name.startswith(ROOT_STAGE_TEST_PREFIX):
        raise AssertionError("refusing non-test root-stage cleanup")
    if path.is_symlink():
        raise AssertionError("refusing symlink root-stage cleanup")
    if not path.exists():
        return
    for directory, child_directories, files in os.walk(
            path, topdown=False, followlinks=False):
        current = Path(directory)
        for name in (*child_directories, *files):
            child = current / name
            if child.is_symlink():
                raise AssertionError("refusing symlink in stage-test tree")
        for name in files:
            os.chmod(current / name, 0o600)
        for name in child_directories:
            os.chmod(current / name, 0o700)
        os.chmod(current, 0o700)
    shutil.rmtree(path)


@contextmanager
def unprivileged_root_stage_fixture(label: str):
    """Provide four exact sources and a disjoint synthetic v6 stage plan."""
    workspace = Path(tempfile.mkdtemp(
        prefix=ROOT_STAGE_TEST_PREFIX + label + "-", dir="/dev/shm"))
    transition_id = workspace.name
    parent = Path("/dev/shm") / (transition_id + "-root-code")
    output = Path("/dev/shm") / (transition_id + "-output")
    raws = {
        "candidate": b"candidate-stage-v6\n",
        "controller": b"controller-stage-v6\n",
        "legacy": b"legacy-stage-v6\n",
        "p32": b"p32-stage-v6\n",
    }
    paths = {name: workspace / (name + ".py") for name in raws}
    try:
        for name, path in paths.items():
            path.write_bytes(raws[name])
            os.chmod(path, 0o444)
        plan = replace(
            transition.TransitionPlan(), transition_id=transition_id,
            root=output, controller_uid=os.geteuid(),
            controller_gid=os.getegid(), code_owner_uid=os.geteuid(),
            code_owner_gid=os.getegid(), code_anchor=Path("/dev/shm"),
            code_seal_parent=parent, code_seal_root=parent / "seal-v6",
            candidate_sha256=hashlib.sha256(raws["candidate"]).hexdigest(),
            controller_sha256=hashlib.sha256(raws["controller"]).hexdigest(),
            old_source_sha256=hashlib.sha256(raws["legacy"]).hexdigest(),
            p32_sha256=hashlib.sha256(raws["p32"]).hexdigest(),
            stage_candidate_source=paths["candidate"],
            stage_controller_source=paths["controller"],
            stage_legacy_source=paths["legacy"],
            stage_p32_source=paths["p32"], old_source=paths["legacy"])
        yield plan, paths, raws
    finally:
        for cleanup_path in (output, parent, workspace):
            remove_root_stage_test_tree(cleanup_path)


@contextmanager
def exact_prepare_runtime(
        plan, *, environment=None, orig_argv=None, flags=None,
        executable=None):
    """Expose the exact isolated prepare process contract to direct tests."""
    observed_environment = transition.prepare_environment(plan) \
        if environment is None else dict(environment)
    observed_orig_argv = transition.expected_prepare_orig_argv(plan) \
        if orig_argv is None else tuple(orig_argv)
    observed_flags = dict(transition.PREPARE_FLAG_CONTRACT) \
        if flags is None else dict(flags)
    with ExitStack() as stack:
        stack.enter_context(mock.patch.dict(
            os.environ, observed_environment, clear=True))
        stack.enter_context(mock.patch.object(
            sys, "orig_argv", list(observed_orig_argv)))
        stack.enter_context(mock.patch.object(
            sys, "flags", SimpleNamespace(**observed_flags)))
        if executable is not None:
            stack.enter_context(mock.patch.object(
                sys, "executable", executable))
        yield


@contextmanager
def exact_privileged_helper_runtime(
        runtime, *, environment=None, orig_argv=None, flags=None):
    """Expose one exact isolated identity/signal-helper runtime."""
    observed_environment = runtime["environment"] \
        if environment is None else dict(environment)
    observed_orig_argv = runtime["sys_orig_argv"] \
        if orig_argv is None else list(orig_argv)
    observed_flags = runtime["flags"] if flags is None else dict(flags)
    with mock.patch.dict(os.environ, observed_environment, clear=True), \
            mock.patch.object(sys, "orig_argv", observed_orig_argv), \
            mock.patch.object(
                sys, "flags", SimpleNamespace(**observed_flags)):
        yield


def make_test_code_seal(base: Path, *, candidate: bytes = b"candidate\n",
                        p32: bytes = b"p32\n",
                        controller: bytes = b"controller\n",
                        legacy: bytes = b"legacy\n",
                        output_name: str = "output"):
    """Create a non-root synthetic equivalent of the production code seal."""
    anchor = base / "code-anchor"
    parent = anchor / "transition-root-code"
    code_root = parent / "seal-v6"
    anchor.mkdir()
    parent.mkdir()
    code_root.mkdir()
    os.chmod(parent, 0o700)
    uid = os.geteuid()
    gid = os.getegid()
    replacement = (
        "/usr/bin/python3.12", "-I", "-S", "-B",
        str(code_root / "legacy_wirehair_expo_thermal_sampler.py"),
        "--csv", str(transition.OLD_CSV),
        "--pid-file", str(transition.OLD_PID_FILE),
    )
    replacement_sha = hashlib.sha256(
        b"\0".join(value.encode("ascii") for value in replacement) +
        b"\0").hexdigest()
    transition_id = transition.SOURCE_TEST_TRANSITION_PREFIX + \
        hashlib.sha256(str(base).encode("ascii")).hexdigest()[:16]
    plan = replace(
        transition.TransitionPlan(), transition_id=transition_id,
        root=base / output_name,
        controller_uid=uid, controller_gid=gid,
        code_owner_uid=uid, code_owner_gid=gid,
        code_anchor=anchor, code_seal_parent=parent,
        code_seal_root=code_root,
        candidate_sha256=hashlib.sha256(candidate).hexdigest(),
        p32_sha256=hashlib.sha256(p32).hexdigest(),
        controller_sha256=hashlib.sha256(controller).hexdigest(),
        old_source_sha256=hashlib.sha256(legacy).hexdigest(),
        replacement_old_argv=replacement,
        replacement_old_cmdline_sha256=replacement_sha)
    files = {
        plan.sampler: candidate,
        plan.p32: p32,
        plan.controller: controller,
        plan.legacy_sampler: legacy,
    }
    for path, raw in files.items():
        path.write_bytes(raw)
        os.chmod(path, 0o444)
    os.chmod(code_root, 0o555)
    file_records = {
        name: {"binding": transition.file_binding(path, with_hash=True),
               "path": str(path)}
        for name, path in (
            ("candidate", plan.sampler), ("controller", plan.controller),
            ("legacy", plan.legacy_sampler), ("p32", plan.p32))
    }
    source_paths = {
        "candidate": plan.stage_candidate_source,
        "controller": plan.stage_controller_source,
        "legacy": plan.stage_legacy_source,
        "p32": plan.stage_p32_source,
    }
    source_records = {
        name: {
            "binding": dict(file_records[name]["binding"]),
            "fd": 10 + index, "path": str(source_paths[name]),
            "stability_observations": 2,
        }
        for index, name in enumerate(
            ("candidate", "controller", "legacy", "p32"))
    }
    anchor_binding = transition.directory_binding(anchor)
    anchor_binding["mode"] = 0o700
    parent_binding = transition.directory_binding(parent)
    parent_binding["mode"] = 0o555
    root_binding = transition.directory_binding(code_root)
    receipt = transition.sealed_record(
        transition.ROOT_CODE_STAGE_RECEIPT_SCHEMA, {
            "authority": transition._root_stage_program_hashes(),
            "directories": {
                "anchor": anchor_binding, "parent": parent_binding,
                "root": root_binding},
            "files": file_records,
            "no_live_state_or_workload": True,
            "partial_or_residual_stage_policy":
                "hard-stop-no-blind-retry",
            "sources": source_records,
            "transition_id": plan.transition_id,
        })
    transition.write_new(
        plan.root_code_stage_receipt, transition.canonical_json(receipt),
        owner_uid=uid)
    os.chmod(parent, 0o555)
    os.chmod(anchor, 0o700)
    return plan


@contextmanager
def real_p32_old_restart_fixture(
        base: Path, csv_observations, *, now_observations,
        classifier=None, success=False):
    """Drive the actual old-restart loop around the frozen P32 CSV parser."""
    plan = make_test_code_seal(base, legacy=b"old source\n")
    plan = replace(
        plan, old_csv=base / "legacy.csv",
        old_pid_file=base / "legacy.pid")
    source_binding = transition.file_binding(
        plan.legacy_sampler, with_hash=True)
    raw_pid = b"7702\n"
    final_raw = bytes(csv_observations[-1])
    csv_binding = {
        "device": 31, "gid": 0, "inode": 310, "mode": 0o444,
        "nlink": 1, "size": len(final_raw), "uid": 0,
    }
    pid_binding = {
        "device": 32, "gid": 0, "inode": 320, "mode": 0o444,
        "nlink": 1, "sha256": hashlib.sha256(raw_pid).hexdigest(),
        "size": len(raw_pid), "uid": 0,
    }
    classifier_call = mock.Mock(wraps=(
        frozen_p32._bootstrap_thermal_csv_state
        if classifier is None else classifier))
    parser_call = mock.Mock(wraps=frozen_p32._parse_thermal_csv)
    cleanup = mock.Mock()
    backend = object.__new__(transition.LiveBackend)
    backend.plan = plan
    backend.tools = transition.capture_tool_records()
    backend.design = {"synthetic": "real-frozen-p32-bootstrap"}
    backend.recovery_budget = transition.EmergencyRecoveryBudget(60.0)
    backend.p32 = SimpleNamespace(
        TimingError=frozen_p32.TimingError,
        _bootstrap_thermal_csv_state=classifier_call,
        _parse_thermal_csv=parser_call,
        THERMAL_FIELDS=frozen_p32.THERMAL_FIELDS,
        _kill_owned_process_session=cleanup,
        capture_process_identity=mock.Mock(side_effect=AssertionError(
            "unprivileged replacement identity capture ran")),
    )
    backend._require_old_launcher_absent = mock.Mock()
    backend._i2c_readers = mock.Mock(
        side_effect=((), (7702,)) if success else ((),))
    backend._fuser = mock.Mock(return_value=(7702,))
    backend._replacement_launch_command = mock.Mock(
        return_value=("/usr/bin/false",))
    backend._recovery_wait_deadline = mock.Mock(return_value=1.0)
    backend._recovery_now = mock.Mock(side_effect=tuple(now_observations))
    helper_identity = {
        "affinity": str(plan.old_cpu),
        "argv": list(plan.replacement_old_argv),
        "cmdline_sha256": plan.replacement_old_cmdline_sha256,
        "executable": {
            "device": backend.tools["python3"]["binding"]["device"],
            "inode": backend.tools["python3"]["binding"]["inode"],
        },
        "pid": 7702, "ppid": 7701, "process_group": 7701,
        "processor": plan.old_cpu, "session": 7701,
        "start_tick": 456, "state": "S", "uids": [0, 0, 0, 0],
    }
    backend._inspect_replacement_identities = mock.Mock(return_value={
        "boot_id": "boot",
        "targets": {
            "child": {"identity": helper_identity, "state": "present"},
            "launcher": {"state": "absent"},
        },
    })
    backend._inspect_sampler_fds = mock.Mock(return_value={
        "synthetic": "exact legacy FD receipt"})
    backend._replacement_launcher_roster = mock.Mock(return_value=[])
    process = mock.Mock(pid=7701)
    process.poll.return_value = 0
    csv_values = iter(bytes(value) for value in csv_observations)

    def stable(path, attempts=20):
        del attempts
        if path == plan.old_pid_file:
            return raw_pid
        if path == plan.old_csv:
            return next(csv_values)
        raise AssertionError("unexpected stable read: %s" % path)

    def binding(path, *, with_hash):
        del with_hash
        if path == plan.legacy_sampler:
            return source_binding
        if path == plan.old_csv:
            return csv_binding
        if path == plan.old_pid_file:
            return pid_binding
        raise AssertionError("unexpected binding path: %s" % path)

    code_seal = {"files": {"legacy": {"binding": source_binding}}}
    with mock.patch.object(
            transition, "verify_code_seal", return_value=code_seal), \
            mock.patch.object(
                transition, "execute_environment", return_value={}), \
            mock.patch.object(
                transition, "stable_file_bytes", side_effect=stable) as reads, \
            mock.patch.object(
                transition, "file_binding", side_effect=binding), \
            mock.patch.object(
                transition.subprocess, "Popen", return_value=process), \
            mock.patch.object(
                transition, "capture_owned_session_leader", return_value=123), \
            mock.patch.object(
                transition.Path, "read_text", return_value="boot\n"), \
            mock.patch.object(transition.time, "sleep") as sleep:
        yield SimpleNamespace(
            backend=backend, plan=plan, process=process,
            classifier=classifier_call, parser=parser_call, cleanup=cleanup,
            reads=reads, sleep=sleep, raw_pid=raw_pid,
            source_binding=source_binding)


def make_fd_snapshot(*, kind: str = "candidate"):
    source_path = "/seal/candidate.py" if kind == "candidate" else \
        "/seal/legacy.py"
    csv_path = "/output/thermal.csv"
    argv = (
        "/usr/bin/python3.12", "-I", "-S", "-B", source_path,
        "--csv", csv_path,
    )
    source_raw = b"candidate-source\n" if kind == "candidate" else \
        b"legacy-source\n"
    source = {
        "device": 20, "gid": 0, "inode": 200, "mode": 0o444,
        "nlink": 1, "sha256": hashlib.sha256(source_raw).hexdigest(),
        "size": len(source_raw), "uid": 0,
    }
    csv_binding = {
        "device": 21, "gid": 0, "inode": 210, "mode": 0o444,
        "nlink": 1, "size": 120, "uid": 0,
    }

    def device(bus):
        return {
            "device": 22, "gid": 0, "inode": 220 + bus,
            "kind": "char", "major": 89, "minor": bus,
            "mode": 0o660, "nlink": 1, "rdev": os.makedev(89, bus),
            "size": 0, "uid": 0,
        }

    def fd(fd_number, binding, path, access_mode, *, kind_name="regular",
           status_flags=0):
        return {
            "access_mode": access_mode, "device": binding["device"],
            "fd": fd_number, "flags": access_mode | status_flags,
            "gid": binding["gid"], "inode": binding["inode"],
            "kind": kind_name, "link_path": path,
            "major": binding.get("major", 0),
            "minor": binding.get("minor", 0), "mode": binding["mode"],
            "nlink": binding["nlink"], "rdev": binding.get("rdev", 0),
            "size": binding["size"], "uid": binding["uid"],
        }

    i2c1 = device(1)
    i2c2 = device(2)
    python = {"device": 30, "inode": 300}
    identity = {
        "affinity": "127", "argv": list(argv),
        "cmdline_sha256": hashlib.sha256(
            b"\0".join(value.encode("ascii") for value in argv) +
            b"\0").hexdigest(),
        "executable": python, "pid": 7001, "ppid": 7000,
        "process_group": 7000, "processor": 127, "session": 7000,
        "start_tick": 12345, "state": "S", "uids": [0, 0, 0, 0],
    }
    fds = [
        fd(3, source, source_path, os.O_RDONLY,
           status_flags=os.O_NONBLOCK),
        fd(4, csv_binding, csv_path,
           os.O_RDWR if kind == "candidate" else os.O_WRONLY),
        fd(5, i2c1, "/dev/i2c-1", os.O_RDWR, kind_name="char"),
        fd(6, i2c2, "/dev/i2c-2", os.O_RDWR, kind_name="char"),
    ]
    if kind == "legacy":
        fds.pop(0)
    snapshot = {
        "boot_id": "00000000-0000-0000-0000-000000000001",
        "csv": {"binding": csv_binding, "path": csv_path},
        "fds": fds,
        "i2c": {"/dev/i2c-1": i2c1, "/dev/i2c-2": i2c2},
        "identity": identity, "pid": 7001,
        "source": {"binding": source, "path": source_path},
    }
    contract = {
        "kind": kind, "expected_pid": 7001,
        "expected_start_tick": 12345,
        "expected_boot_id": snapshot["boot_id"], "expected_argv": argv,
        "expected_python": python,
        "expected_source_sha256": source["sha256"],
    }
    return snapshot, contract


def make_identity_request(plan, *, profile="original",
                          allowed_absence="none", child_pid=None,
                          child_start_tick=None, launcher_pid=None,
                          launcher_start_tick=None, controller_pid=8123):
    """Build one canonical synthetic request without inspecting live state."""
    if profile == "original":
        child_pid = plan.old_pid if child_pid is None else child_pid
        child_start_tick = plan.old_start_tick \
            if child_start_tick is None else child_start_tick
        launcher_pid = plan.old_launcher_pid \
            if launcher_pid is None else launcher_pid
        launcher_start_tick = plan.old_launcher_start_tick \
            if launcher_start_tick is None else launcher_start_tick
    else:
        child_pid = 7002 if child_pid is None else child_pid
        child_start_tick = 456 if child_start_tick is None \
            else child_start_tick
        launcher_pid = 7001 if launcher_pid is None else launcher_pid
        launcher_start_tick = 123 if launcher_start_tick is None \
            else launcher_start_tick
    return transition.identity_inspection_request(
        plan, profile=profile, child_pid=child_pid,
        child_start_tick=child_start_tick, launcher_pid=launcher_pid,
        launcher_start_tick=launcher_start_tick,
        controller_pid=controller_pid, allowed_absence=allowed_absence)


def identity_from_contract(contract, *, captured_start_tick=456):
    """Materialize the exact process record produced by `_proc_identity_at`."""
    affinity = contract["affinity"] or "16-59,81-123"
    return {
        "affinity": affinity,
        "argv": list(contract["argv"]),
        "cmdline_sha256": contract["cmdline_sha256"],
        "executable": dict(contract["executable"]),
        "pid": contract["pid"],
        "ppid": contract["allowed_ppids"][0],
        "process_group": contract["process_group"],
        "processor": int(affinity.split(",", 1)[0].split("-", 1)[0]),
        "session": contract["session"],
        "start_tick": captured_start_tick
        if contract["capture_start_tick"] else contract["start_tick"],
        "state": "S",
        "uids": list(contract["uids"]),
    }


def synthetic_code_seal_receipt(plan):
    return transition.sealed_record(
        transition.CODE_SEAL_RECEIPT_SCHEMA, {
            "synthetic_test_only": True,
            "transition_id": plan.transition_id,
        })


def verify_synthetic_identity_receipt(plan, receipt, request, tools,
                                      code_seal):
    with mock.patch.object(
            transition, "verify_code_seal", return_value=code_seal):
        return transition.verify_identity_inspection_receipt(
            plan, receipt, request, tools)


def make_identity_receipt(plan, request, tools, *, absent=(), exiting=()):
    """Create and verify one exact synthetic privileged helper receipt."""
    contracts = transition.identity_target_contracts(plan, request, tools)
    targets = {}
    for name in ("child", "launcher"):
        if name in absent:
            targets[name] = {
                "expected_pid": contracts[name]["pid"],
                "expected_start_tick": contracts[name]["start_tick"],
                "pidfd_absence_observations": 2,
                "state": "absent",
            }
        elif name in exiting:
            identity = identity_from_contract(contracts[name])
            targets[name] = transition._exiting_identity_target(
                contracts[name], snapshot_kind="pidfd-ready-after-identity",
                proc_stat=transition._identity_proc_stat(identity),
                identity=identity, observations=2)
        else:
            targets[name] = {
                "identity": identity_from_contract(contracts[name]),
                "pidfd_opened_before_target_proc": True,
                "pidfd_unready_after_replay": True,
                "proc_identity_observations": 2,
                "state": "present",
            }
    return make_identity_receipt_from_capture(
        plan, request, tools, {
            "boot_id": plan.old_boot_id,
            "targets": targets,
        })


def make_identity_receipt_from_capture(plan, request, tools, captured):
    """Seal and replay an actual synthetic identity-core capture."""
    code_seal = synthetic_code_seal_receipt(plan)
    receipt = transition.sealed_record(
        transition.IDENTITY_INSPECTION_SCHEMA, {
            **copy.deepcopy(captured),
            "code_seal": code_seal,
            "helper_runtime": transition._runtime_with_code_seal(
                transition.identity_inspection_runtime_contract(
                    plan, request, tools), code_seal),
            "pidfd_policy": transition.IDENTITY_PIDFD_POLICY,
            "request": dict(request),
            "tools": dict(tools),
            "transition_id": plan.transition_id,
        })
    return verify_synthetic_identity_receipt(
        plan, receipt, request, tools, code_seal)


def make_process_signal_receipt(plan, identity_request, tools, *, target,
                                signum=signal.SIGTERM, pidfd=101):
    request = transition.process_signal_request(
        plan, identity_request, target=target, signum=signum)
    contracts = transition.identity_target_contracts(
        plan, identity_request, tools)
    targets = {
        name: {
            "identity": identity_from_contract(contract),
            "pidfd_opened_before_target_proc": True,
            "pidfd_unready_after_replay": True,
            "proc_identity_observations": 2,
            "state": "present",
        }
        for name, contract in contracts.items()
    }
    code_seal = synthetic_code_seal_receipt(plan)
    action = {
        "boot_id": plan.old_boot_id,
        "identity_target_sha256": hashlib.sha256(
            transition.canonical_json(targets[target])).hexdigest(),
        "pid": contracts[target]["pid"],
        "pidfd": pidfd,
        "pidfd_held_from_before_target_proc_through_signal": True,
        "pidfd_send_signal_result": 0,
        "signal": int(signum),
        "signal_name": signal.Signals(signum).name,
        "syscall": "pidfd_send_signal",
        "target": target,
    }
    receipt = transition.sealed_record(
        transition.PROCESS_SIGNAL_SCHEMA, {
            "boot_id": plan.old_boot_id,
            "code_seal": code_seal,
            "helper_runtime": transition._runtime_with_code_seal(
                transition.process_signal_runtime_contract(
                    plan, request, tools), code_seal),
            "pidfd_policy": transition.PROCESS_SIGNAL_PIDFD_POLICY,
            "request": request,
            "signal_action": action,
            "targets": targets,
            "tools": dict(tools),
            "transition_id": plan.transition_id,
        })
    with mock.patch.object(
            transition, "verify_code_seal", return_value=code_seal):
        return transition.verify_process_signal_receipt(
            plan, receipt, request, tools)


class FakeClock:
    def __init__(self) -> None:
        self.value = 100.0

    def __call__(self) -> float:
        return self.value


class FakeJournal:
    def __init__(self, fail_at=None, *, clock=None, advance_at=None) -> None:
        self.fail_at = fail_at
        self.clock = clock
        self.advance_at = advance_at
        self.records = []
        self.replay_calls = 0

    def record(self, phase, status, payload):
        self.records.append((phase, status, dict(payload)))
        if self.advance_at == (phase, status) and self.clock is not None:
            self.clock.value += 541.0
        if self.fail_at == (phase, status):
            raise transition.TransitionError("injected journal failure")
        return {"phase": phase, "status": status}

    def replay(self):
        self.replay_calls += 1
        if self.advance_at == (
                "receipt_chain", "replay-%d" % self.replay_calls) and \
                self.clock is not None:
            self.clock.value += 541.0
        if self.fail_at == ("receipt_chain", "replay"):
            raise transition.TransitionError("injected journal failure")
        return {
            "count": len(self.records),
            "head_sha256": "c" * 64,
            "roster": [
                {"sequence": index, "phase": phase, "status": status}
                for index, (phase, status, _payload) in enumerate(self.records)
            ],
        }


class FakeSignalGuard:
    def __init__(self) -> None:
        self.requested = False

    def raise_if_requested(self):
        if self.requested:
            raise transition.TransitionError("deferred controller signal: SIGTERM")


class FakeBackend:
    NORMAL = (
        "hard_preflight", "arm_recovery", "stop_old", "archive_old",
        "start_candidate", "exercise_candidate", "stop_candidate",
        "accept_candidate",
    )

    def __init__(self, fail_at=None, clock=None, advance_at=None,
                 *, forced_stop=False, forced_archive=False) -> None:
        self.fail_at = fail_at
        self.clock = clock
        self.advance_at = advance_at
        self.calls = []
        self.recovery_budget = None
        self.forced_stop = forced_stop
        self.forced_archive = forced_archive
        self.signal_actions = []
        self.signal_failures = {}

    def _call(self, name, value=None):
        self.calls.append(name)
        if self.advance_at == name and self.clock is not None:
            self.clock.value += 500.0
        if self.fail_at == name:
            raise transition.TransitionError("injected " + name)
        return value if value is not None else {"phase": name}

    def hard_preflight(self):
        return self._call("hard_preflight")

    def arm_recovery(self):
        return self._call("arm_recovery")

    def stop_old(self):
        return self._call("stop_old", {"forced": self.forced_stop})

    def archive_old(self):
        return self._call("archive_old", {
            "binding": {"sha256": "a" * 64},
            "path": "/archive", "forced_stop": self.forced_archive,
        })

    def start_candidate(self):
        return self._call("start_candidate")

    def exercise_candidate(self):
        return self._call("exercise_candidate")

    def stop_candidate(self):
        return self._call("stop_candidate")

    def accept_candidate(self):
        return self._call("accept_candidate", {
            "sample_count": 6, "raw_sha256": "d" * 64})

    def cleanup_candidate(self):
        return self._call("cleanup_candidate")

    def begin_emergency_recovery(self, budget):
        self.recovery_budget = budget
        return {"budget_installed": True}

    def recovery_action_evidence(self):
        return transition.sealed_record(
            transition.RECOVERY_ACTION_EVIDENCE_SCHEMA, {
                "failed_attempts": copy.deepcopy(self.signal_failures),
                "receipt_count": len(self.signal_actions),
                "signal_receipts": copy.deepcopy(self.signal_actions),
                "transition_id": transition.TRANSITION_ID,
            })

    def restore_old(self, archive):
        if archive is not None:
            assert archive["binding"]["sha256"] == "a" * 64
        return self._call("restore_old", {
            "pid": 7001,
            "csv_live_identity": {"device": 1, "inode": 2},
            "archived_pre_dry_sha256": "a" * 64,
        })

    def publish_audit_binding(self, archive, restored, dry_accepted,
                              candidate_accept, receipt_chain_prefix):
        assert restored["pid"] == 7001
        assert receipt_chain_prefix["count"] > 0
        if dry_accepted:
            assert candidate_accept["sample_count"] == 6
        return self._call("publish_audit_binding", {
            "dry_accepted": dry_accepted, "sha256": "b" * 64,
        })

    def final_replay(self, candidate_accept, archive, restored, audit):
        assert restored["pid"] == 7001
        assert audit["sha256"] == "b" * 64
        return self._call("final_replay", {
            "candidate_accept": candidate_accept,
            "archive": archive,
            "receipt": "replayed",
        })


class TransitionStateMachineTests(unittest.TestCase):
    def make_controller(self, backend, journal=None, clock=None,
                        signal_guard=None):
        clock = clock or FakeClock()
        deadline = transition.Deadline(540.0, 60.0, clock=clock)
        return transition.TransitionController(
            transition.TransitionPlan(controller_uid=os.geteuid()), backend,
            journal or FakeJournal(), deadline, signal_guard), clock

    def actual_signal_evidence_backend(self, receipts, failures=None):
        plan = transition.TransitionPlan(controller_sha256="a" * 64)
        backend = FakeBackend()
        backend.plan = plan
        backend.original_signal_receipts = list(receipts)
        backend.original_signal_failures = dict(failures or {})
        backend.recovery_action_evidence = lambda: \
            transition.LiveBackend.recovery_action_evidence(backend)
        return backend, plan

    def assert_action_evidence(self, value, expected_receipts,
                               expected_failures):
        verified = transition.verify_sealed(
            value, transition.RECOVERY_ACTION_EVIDENCE_SCHEMA,
            "synthetic controller recovery action evidence")
        self.assertEqual(set(verified), {
            "failed_attempts", "receipt_count", "schema",
            "self_sha256_excluding_field", "signal_receipts",
            "transition_id",
        })
        self.assertEqual(verified["receipt_count"], len(expected_receipts))
        self.assertEqual(verified["signal_receipts"], expected_receipts)
        self.assertEqual(verified["failed_attempts"], expected_failures)
        self.assertEqual(verified["transition_id"], transition.TRANSITION_ID)
        self.assertEqual(
            hashlib.sha256(transition.canonical_json({
                key: item for key, item in verified.items()
                if key != "self_sha256_excluding_field"
            })).hexdigest(), verified["self_sha256_excluding_field"])

    def test_success_is_nonpromotional_and_restores_old(self):
        backend = FakeBackend()
        journal = FakeJournal()
        controller, _clock = self.make_controller(backend, journal)
        terminal = controller.run()
        self.assertTrue(terminal["success"])
        self.assertTrue(terminal["dry_accepted"])
        self.assertEqual(terminal["future_audit_binding"]["sha256"], "b" * 64)
        self.assertEqual(backend.calls, [
            *FakeBackend.NORMAL, "cleanup_candidate", "restore_old",
            "publish_audit_binding", "final_replay",
        ])
        self.assertIn(("recovery_arm", "completed", {"phase": "arm_recovery"}),
                      journal.records)
        self.assertEqual(terminal["candidate_accept"]["sample_count"], 6)
        self.assertEqual(terminal["final_replay"]["receipt"], "replayed")
        self.assertGreater(terminal["receipt_chain_prefix"]["count"], 0)
        checkpoints = [payload for phase, status, payload in journal.records
                       if (phase, status) == ("terminal", "started")]
        self.assertEqual(len(checkpoints), 1)
        self.assertIsNone(checkpoints[0]["transition_success"])
        self.assertNotIn("success", checkpoints[0])

    def test_forced_old_stop_vetoes_candidate_and_enters_recovery_immediately(self):
        backend = FakeBackend(forced_stop=True)
        controller, _clock = self.make_controller(backend)
        with self.assertRaisesRegex(
                transition.TransitionError, "forced old stop vetoes"):
            controller.run()
        self.assertNotIn("archive_old", backend.calls)
        self.assertNotIn("start_candidate", backend.calls)
        self.assertEqual(backend.calls[:3], [
            "hard_preflight", "arm_recovery", "stop_old"])
        self.assertIn("cleanup_candidate", backend.calls)
        self.assertIn("restore_old", backend.calls)

    def test_unclean_archive_vetoes_candidate_before_launch(self):
        backend = FakeBackend(forced_archive=True)
        controller, _clock = self.make_controller(backend)
        with self.assertRaisesRegex(
                transition.TransitionError, "forced old stop vetoes"):
            controller.run()
        self.assertIn("archive_old", backend.calls)
        self.assertNotIn("start_candidate", backend.calls)
        self.assertIn("restore_old", backend.calls)

    def test_preflight_failure_never_enters_recovery(self):
        backend = FakeBackend(fail_at="hard_preflight")
        controller, _clock = self.make_controller(backend)
        with self.assertRaisesRegex(transition.TransitionError, "hard_preflight"):
            controller.run()
        self.assertEqual(backend.calls, ["hard_preflight"])

    def test_recovery_arm_failure_never_stops_old(self):
        backend = FakeBackend(fail_at="arm_recovery")
        controller, _clock = self.make_controller(backend)
        with self.assertRaisesRegex(transition.TransitionError, "arm_recovery"):
            controller.run()
        self.assertEqual(backend.calls, ["hard_preflight", "arm_recovery"])

    def test_every_post_arm_boundary_failure_enters_recovery(self):
        for failure in FakeBackend.NORMAL[2:]:
            with self.subTest(failure=failure):
                backend = FakeBackend(fail_at=failure)
                controller, _clock = self.make_controller(backend)
                with self.assertRaisesRegex(transition.TransitionError, failure):
                    controller.run()
                self.assertIn("cleanup_candidate", backend.calls)
                self.assertIn("restore_old", backend.calls)
                self.assertIn("publish_audit_binding", backend.calls)
                self.assertLess(
                    backend.calls.index("cleanup_candidate"),
                    backend.calls.index("restore_old"))

    def test_complete_receipt_failure_after_stop_still_recovers(self):
        backend = FakeBackend()
        journal = FakeJournal(fail_at=("old_stop", "completed"))
        controller, _clock = self.make_controller(backend, journal)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "journal failure"):
            controller.run()
        self.assertIn("cleanup_candidate", backend.calls)
        self.assertIn("restore_old", backend.calls)

    def test_stop_completion_journal_failure_retains_canonical_action_evidence(self):
        plan = transition.TransitionPlan(controller_sha256="a" * 64)
        tools = transition.capture_tool_records()
        action = make_process_signal_receipt(
            plan, make_identity_request(
                plan, allowed_absence="launcher"),
            tools, target="child")
        backend, _plan = self.actual_signal_evidence_backend([])
        backend.fail_at = "restore_old"

        def stop_old():
            backend.calls.append("stop_old")
            backend.original_signal_receipts.append(action)
            return {"forced": False, "signal_receipts": [action]}

        backend.stop_old = stop_old
        journal = FakeJournal(fail_at=("old_stop", "completed"))
        controller, _clock = self.make_controller(backend, journal)
        with self.assertRaisesRegex(
                transition.TransitionError, "journal failure"):
            controller.run()
        completion = [
            payload for phase, status, payload in journal.records
            if (phase, status) == ("old_stop", "completed")]
        recovery_finished = [
            payload for phase, status, payload in journal.records
            if (phase, status) ==
            ("process_signal_actions", "completed") and
            payload["checkpoint"] == "recovery-finished"]
        terminal = [
            payload for phase, status, payload in journal.records
            if (phase, status) == ("terminal", "started")]
        self.assertEqual(completion[0]["signal_receipts"], [action])
        self.assert_action_evidence(
            recovery_finished[-1]["evidence"], [action], {})
        self.assert_action_evidence(
            terminal[-1]["process_signal_actions"], [action], {})
        self.assertIsNone(terminal[-1]["old_sampler_restored"])
        self.assertIn("cleanup_candidate", backend.calls)
        self.assertIn("restore_old", backend.calls)

    def test_late_accept_completion_retains_result_but_is_not_accepted(self):
        clock = FakeClock()
        backend = FakeBackend()
        journal = FakeJournal(
            clock=clock, advance_at=("candidate_accept", "completed"))
        controller, _clock = self.make_controller(
            backend, journal=journal, clock=clock)
        with self.assertRaisesRegex(
                transition.TransitionError, "normal deadline"):
            controller.run()
        self.assertEqual(controller.candidate_accept["sample_count"], 6)
        self.assertFalse(controller.dry_accepted)
        self.assertIn("publish_audit_binding", backend.calls)

    def test_stop_start_receipt_failure_leaves_old_untouched(self):
        backend = FakeBackend()
        journal = FakeJournal(fail_at=("old_stop", "started"))
        controller, _clock = self.make_controller(backend, journal)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "journal failure"):
            controller.run()
        self.assertNotIn("stop_old", backend.calls)
        self.assertNotIn("cleanup_candidate", backend.calls)
        self.assertNotIn("restore_old", backend.calls)

    def test_stop_start_receipt_overrun_leaves_old_untouched(self):
        clock = FakeClock()
        backend = FakeBackend()
        journal = FakeJournal(
            clock=clock, advance_at=("old_stop", "started"))
        controller, _clock = self.make_controller(
            backend, journal=journal, clock=clock)
        with self.assertRaisesRegex(
                transition.TransitionError, "normal deadline"):
            controller.run()
        self.assertNotIn("stop_old", backend.calls)
        self.assertNotIn("cleanup_candidate", backend.calls)
        self.assertNotIn("restore_old", backend.calls)

    def test_post_stop_phase_start_overrun_recovers_without_running_action(self):
        clock = FakeClock()
        backend = FakeBackend()
        journal = FakeJournal(
            clock=clock, advance_at=("candidate_start", "started"))
        controller, _clock = self.make_controller(
            backend, journal=journal, clock=clock)
        with self.assertRaisesRegex(
                transition.TransitionError, "normal deadline"):
            controller.run()
        self.assertNotIn("start_candidate", backend.calls)
        self.assertIn("cleanup_candidate", backend.calls)
        self.assertIn("restore_old", backend.calls)

    def test_started_receipt_failure_after_arm_still_recovers(self):
        backend = FakeBackend()
        journal = FakeJournal(fail_at=("old_archive", "started"))
        controller, _clock = self.make_controller(backend, journal)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "journal failure"):
            controller.run()
        self.assertIn("cleanup_candidate", backend.calls)
        self.assertIn("restore_old", backend.calls)

    def test_every_post_stop_phase_receipt_failure_recovers(self):
        boundaries = [("old_stop", "completed")]
        for phase in (
                "old_archive", "candidate_start", "candidate_exercise",
                "candidate_stop", "candidate_accept"):
            boundaries.extend(((phase, "started"), (phase, "completed")))
        for boundary in boundaries:
            with self.subTest(boundary=boundary):
                backend = FakeBackend()
                controller, _clock = self.make_controller(
                    backend, FakeJournal(fail_at=boundary))
                with self.assertRaisesRegex(transition.TransitionError,
                                            "journal failure"):
                    controller.run()
                self.assertIn("cleanup_candidate", backend.calls)
                self.assertIn("restore_old", backend.calls)

    def test_every_post_stop_failed_receipt_failure_still_recovers(self):
        failures = {
            "stop_old": "old_stop",
            "archive_old": "old_archive",
            "start_candidate": "candidate_start",
            "exercise_candidate": "candidate_exercise",
            "stop_candidate": "candidate_stop",
            "accept_candidate": "candidate_accept",
        }
        for action, phase in failures.items():
            with self.subTest(action=action):
                backend = FakeBackend(fail_at=action)
                controller, _clock = self.make_controller(
                    backend, FakeJournal(fail_at=(phase, "failed")))
                with self.assertRaises(transition.TransitionError):
                    controller.run()
                self.assertIn("cleanup_candidate", backend.calls)
                self.assertIn("restore_old", backend.calls)

    def test_every_recovery_receipt_failure_preserves_action_order(self):
        for boundary in (
                ("emergency_recovery", "started"),
                ("candidate_cleanup", "started"),
                ("candidate_cleanup", "completed"),
                ("old_restore", "started"),
                ("old_restore", "completed"),
                ("audit_binding", "started"),
                ("audit_binding", "completed"),
                ("final_replay", "started"),
                ("final_replay", "completed"),
                ("emergency_recovery", "completed"),
                ("terminal", "started")):
            with self.subTest(boundary=boundary):
                backend = FakeBackend()
                controller, _clock = self.make_controller(
                    backend, FakeJournal(fail_at=boundary))
                with self.assertRaisesRegex(transition.TransitionError,
                                            "journal failure"):
                    controller.run()
                self.assertIn("cleanup_candidate", backend.calls)
                self.assertIn("restore_old", backend.calls)
                self.assertLess(
                    backend.calls.index("cleanup_candidate"),
                    backend.calls.index("restore_old"))

    def test_every_recovery_failed_receipt_failure_preserves_safe_progress(self):
        cases = (
            ("cleanup_candidate", "candidate_cleanup", True, True),
            ("restore_old", "old_restore", True, False),
            ("publish_audit_binding", "audit_binding", True, True),
            ("final_replay", "final_replay", True, True),
        )
        for action, phase, restore_expected, audit_expected in cases:
            with self.subTest(action=action):
                backend = FakeBackend(fail_at=action)
                controller, _clock = self.make_controller(
                    backend, FakeJournal(fail_at=(phase, "failed")))
                with self.assertRaises(transition.TransitionError):
                    controller.run()
                self.assertEqual("restore_old" in backend.calls,
                                 restore_expected)
                self.assertEqual("publish_audit_binding" in backend.calls,
                                 audit_expected)

    def test_cleanup_failure_does_not_skip_restore_attempt(self):
        backend = FakeBackend(fail_at="cleanup_candidate")
        controller, _clock = self.make_controller(backend)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "cleanup_candidate"):
            controller.run()
        self.assertIn("restore_old", backend.calls)

    def test_restore_failure_is_terminal_and_skips_audit_publication(self):
        backend = FakeBackend(fail_at="restore_old")
        controller, _clock = self.make_controller(backend)
        with self.assertRaisesRegex(transition.TransitionError, "restore_old"):
            controller.run()
        self.assertNotIn("publish_audit_binding", backend.calls)

    def test_stop_failure_action_is_in_failed_phase_and_terminal_evidence(self):
        plan = transition.TransitionPlan(controller_sha256="a" * 64)
        tools = transition.capture_tool_records()
        request = make_identity_request(
            plan, allowed_absence="launcher")
        action = make_process_signal_receipt(
            plan, request, tools, target="child")
        backend, _plan = self.actual_signal_evidence_backend([])

        def stop_old():
            backend.calls.append("stop_old")
            backend.original_signal_receipts.append(action)
            raise transition.TransitionError(
                "injected failure after verified stop signal")

        backend.stop_old = stop_old
        journal = FakeJournal()
        controller, _clock = self.make_controller(backend, journal)
        with self.assertRaisesRegex(
                transition.TransitionError, "after verified stop signal"):
            controller.run()
        failed = [payload for phase, status, payload in journal.records
                  if (phase, status) == ("old_stop", "failed")]
        terminal = [payload for phase, status, payload in journal.records
                    if (phase, status) == ("terminal", "started")]
        self.assertEqual(len(failed), 1)
        self.assertEqual(len(terminal), 1)
        self.assert_action_evidence(
            failed[0]["process_signal_actions"], [action], {})
        self.assert_action_evidence(
            terminal[0]["process_signal_actions"], [action], {})

    def test_quiesce_action_survives_later_quiesce_failure(self):
        plan = transition.TransitionPlan(controller_sha256="a" * 64)
        tools = transition.capture_tool_records()
        request = make_identity_request(
            plan, allowed_absence="launcher")
        action = make_process_signal_receipt(
            plan, request, tools, target="child")
        signal_failure = {
            "launcher": {
                "error": {
                    "message": "synthetic launcher signal receipt failure",
                    "type": "TransitionError",
                },
                "no_retry": True,
                "signal": int(signal.SIGTERM),
                "target": "launcher",
            },
        }
        backend, _plan = self.actual_signal_evidence_backend([])

        def restore_old(_archive):
            backend.calls.append("restore_old")
            backend.original_signal_receipts.append(action)
            backend.original_signal_failures.update(signal_failure)
            raise transition.TransitionError(
                "injected quiesce failure after verified signal")

        backend.restore_old = restore_old
        journal = FakeJournal()
        controller, _clock = self.make_controller(backend, journal)
        with self.assertRaisesRegex(
                transition.TransitionError, "quiesce failure"):
            controller.run()
        failed = [payload for phase, status, payload in journal.records
                  if (phase, status) == ("old_restore", "failed")]
        terminal = [payload for phase, status, payload in journal.records
                    if (phase, status) == ("terminal", "started")]
        self.assert_action_evidence(
            failed[0]["process_signal_actions"], [action], signal_failure)
        self.assert_action_evidence(
            terminal[0]["process_signal_actions"], [action], signal_failure)
        self.assertIsNone(terminal[0]["old_sampler_restored"])

    def test_quiesce_actions_survive_distinct_archive_and_launch_failures(self):
        plan = transition.TransitionPlan(controller_sha256="a" * 64)
        tools = transition.capture_tool_records()
        actions = [
            make_process_signal_receipt(
                plan, make_identity_request(
                    plan, allowed_absence="launcher"), tools,
                target="child", pidfd=101),
            make_process_signal_receipt(
                plan, make_identity_request(
                    plan, allowed_absence="child"), tools,
                target="launcher", pidfd=102),
        ]
        signal_failures = {
            "launcher": {
                "error": {
                    "message": "synthetic later launcher signal failure",
                    "type": "TransitionError",
                },
                "no_retry": True,
                "signal": int(signal.SIGTERM),
                "target": "launcher",
            },
        }
        for boundary in ("archive", "launch"):
            with self.subTest(boundary=boundary):
                backend, _plan = self.actual_signal_evidence_backend([])
                backend.candidate_cleanup_complete = True

                def quiesce():
                    backend.calls.append("_quiesce_old_for_recovery")
                    backend.original_signal_receipts.extend(actions)
                    backend.original_signal_failures.update(signal_failures)
                    return {
                        "absence_receipts": {
                            "child": {"state": "synthetic-absent"},
                            "launcher": {"state": "synthetic-absent"},
                        },
                        "i2c_readers_after": [],
                        "signal_attempt_failures": copy.deepcopy(
                            signal_failures),
                        "signal_receipts": copy.deepcopy(actions),
                        "signal_receipts_before_recovery": 0,
                    }

                backend._quiesce_old_for_recovery = mock.Mock(
                    side_effect=quiesce)
                sealed_archive = {
                    "binding": {"sha256": "a" * 64},
                    "forced_stop": False,
                    "path": "/archive",
                    "stale_pid": None,
                }
                if boundary == "archive":
                    backend._ensure_recovery_archive = mock.Mock(
                        side_effect=transition.TransitionError(
                            "injected _ensure_recovery_archive failure"))
                    backend._launch_old = mock.Mock()
                else:
                    backend._ensure_recovery_archive = mock.Mock(
                        return_value=sealed_archive)
                    backend._launch_old = mock.Mock(
                        side_effect=transition.TransitionError(
                            "injected _launch_old failure"))

                def restore_old(archive):
                    backend.calls.append("restore_old")
                    return transition.LiveBackend.restore_old(backend, archive)

                backend.restore_old = restore_old
                journal = FakeJournal()
                controller, _clock = self.make_controller(backend, journal)
                with self.assertRaisesRegex(
                        transition.TransitionError,
                        "injected _%s" % (
                            "ensure_recovery_archive failure"
                            if boundary == "archive" else
                            "launch_old failure")):
                    controller.run()
                backend._quiesce_old_for_recovery.assert_called_once_with()
                backend._ensure_recovery_archive.assert_called_once()
                if boundary == "archive":
                    backend._launch_old.assert_not_called()
                else:
                    backend._launch_old.assert_called_once_with()
                failed = [
                    payload for phase, status, payload in journal.records
                    if (phase, status) == ("old_restore", "failed")][0]
                terminal = [
                    payload for phase, status, payload in journal.records
                    if (phase, status) == ("terminal", "started")][0]
                self.assert_action_evidence(
                    failed["process_signal_actions"], actions,
                    signal_failures)
                self.assert_action_evidence(
                    terminal["process_signal_actions"], actions,
                    signal_failures)
                self.assertIsNone(terminal["old_sampler_restored"])
                finished = [
                    payload for phase, status, payload in journal.records
                    if (phase, status) ==
                    ("process_signal_actions", "completed") and
                    payload["checkpoint"] == "recovery-finished"]
                self.assert_action_evidence(
                    finished[-1]["evidence"], actions, signal_failures)

    def test_failed_terminal_receipt_failure_follows_completed_recovery(self):
        backend = FakeBackend(fail_at="accept_candidate")
        controller, _clock = self.make_controller(
            backend, FakeJournal(fail_at=("terminal", "started")))
        with self.assertRaises(transition.TransitionError):
            controller.run()
        self.assertIn("restore_old", backend.calls)
        self.assertIn("publish_audit_binding", backend.calls)

    def test_audit_binding_failure_is_terminal_after_restore(self):
        backend = FakeBackend(fail_at="publish_audit_binding")
        controller, _clock = self.make_controller(backend)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "publish_audit_binding"):
            controller.run()
        self.assertIn("restore_old", backend.calls)

    def test_signal_received_during_recovery_is_not_lost(self):
        guard = FakeSignalGuard()
        backend = FakeBackend()
        original = backend.publish_audit_binding

        def signal_during_audit(archive, restored, dry_accepted,
                                candidate_accept, receipt_chain_prefix):
            result = original(
                archive, restored, dry_accepted, candidate_accept,
                receipt_chain_prefix)
            guard.requested = True
            return result

        backend.publish_audit_binding = signal_during_audit
        journal = FakeJournal()
        controller, _clock = self.make_controller(
            backend, journal, signal_guard=guard)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "deferred controller signal"):
            controller.run()
        self.assertIn("restore_old", backend.calls)
        terminal = [record for record in journal.records
                    if record[0] == "terminal"]
        self.assertEqual(terminal[-1][1], "started")
        self.assertIsNone(terminal[-1][2]["transition_success"])

    def test_normal_deadline_failure_enters_recovery_reserve(self):
        clock = FakeClock()
        backend = FakeBackend(clock=clock, advance_at="stop_old")
        controller, _clock = self.make_controller(backend, clock=clock)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "normal deadline"):
            controller.run()
        self.assertIn("cleanup_candidate", backend.calls)
        self.assertIn("restore_old", backend.calls)

    def test_expired_scoring_deadline_fails_but_never_vetoes_recovery(self):
        clock = FakeClock()
        backend = FakeBackend(clock=clock, advance_at="stop_old")
        journal = FakeJournal()
        controller, _clock = self.make_controller(
            backend, journal=journal, clock=clock)
        clock.value = 200.0
        backend.advance_at = "stop_old"
        # Advance beyond the full 540 s at the first post-arm operation.
        backend._call = self._absolute_deadline_call(backend, clock)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "scoring absolute deadline exhausted"):
            controller.run()
        self.assertIn("cleanup_candidate", backend.calls)
        self.assertIn("restore_old", backend.calls)
        self.assertIn("publish_audit_binding", backend.calls)
        self.assertIsNotNone(backend.recovery_budget)
        self.assertEqual(backend.recovery_budget.started, clock.value)
        self.assertEqual(
            backend.recovery_budget.absolute,
            clock.value + controller.plan.emergency_recovery_s)
        deadline_receipts = [record for record in journal.records
                             if record[0] == "scoring_deadline"]
        self.assertEqual(len(deadline_receipts), 1)
        self.assertEqual(deadline_receipts[0][1], "failed")
        self.assertTrue(deadline_receipts[0][2]["recovery_actions_continue"])

    def test_terminal_checkpoint_write_crossing_absolute_forces_failure(self):
        clock = FakeClock()
        backend = FakeBackend()
        journal = FakeJournal(
            clock=clock, advance_at=("terminal", "started"))
        controller, _clock = self.make_controller(
            backend, journal=journal, clock=clock)
        with self.assertRaisesRegex(
                transition.TransitionError,
                "scoring absolute deadline exhausted"):
            controller.run()
        self.assertGreaterEqual(clock.value, controller.deadline.absolute)
        self.assertIn("restore_old", backend.calls)
        self.assertIn("final_replay", backend.calls)
        checkpoints = [payload for phase, _status, payload in journal.records
                       if phase == "terminal"]
        self.assertTrue(checkpoints)
        self.assertTrue(all("success" not in payload for payload in checkpoints))
        self.assertTrue(all(payload["transition_success"] is None
                            for payload in checkpoints))
        self.assertFalse(any(phase == "terminal_verdict"
                             for phase, _status, _payload in journal.records))

    def test_terminal_receipt_replay_crossing_absolute_forces_failure(self):
        clock = FakeClock()
        backend = FakeBackend()
        journal = FakeJournal(
            clock=clock, advance_at=("receipt_chain", "replay-2"))
        controller, _clock = self.make_controller(
            backend, journal=journal, clock=clock)
        with self.assertRaisesRegex(
                transition.TransitionError,
                "scoring absolute deadline exhausted"):
            controller.run()
        self.assertIn("restore_old", backend.calls)
        self.assertIn("final_replay", backend.calls)
        self.assertGreaterEqual(
            controller.scoring_evidence_completed_monotonic_s,
            controller.deadline.absolute)

    def test_terminal_uses_one_post_replay_scoring_observation(self):
        class ContradictionClock:
            def __init__(self):
                self.armed = False
                self.calls_after_arm = 0

            def __call__(self):
                if not self.armed:
                    return 100.0
                self.calls_after_arm += 1
                return 100.0 if self.calls_after_arm == 1 else 641.0

        class ArmAfterFinalReplayJournal(FakeJournal):
            def replay(inner_self):
                result = super().replay()
                if inner_self.replay_calls == 2:
                    clock.armed = True
                return result

        clock = ContradictionClock()
        journal = ArmAfterFinalReplayJournal()
        controller, _clock = self.make_controller(
            FakeBackend(), journal=journal, clock=clock)
        terminal = controller.run()
        self.assertTrue(terminal["success"])
        self.assertFalse(terminal["scoring_deadline"]["exhausted"])
        self.assertEqual(clock.calls_after_arm, 1)

    def test_emergency_budget_exhaustion_is_terminal_after_safe_attempts(self):
        clock = FakeClock()
        backend = FakeBackend(clock=clock)
        original = backend.restore_old

        def exhaust_during_restore(archive):
            result = original(archive)
            clock.value += 61.0
            return result

        backend.restore_old = exhaust_during_restore
        journal = FakeJournal()
        controller, _clock = self.make_controller(
            backend, journal=journal, clock=clock)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "emergency recovery budget exhausted"):
            controller.run()
        self.assertIn("restore_old", backend.calls)
        self.assertIn("publish_audit_binding", backend.calls)
        emergency_receipts = [record for record in journal.records
                              if record[0] == "emergency_recovery"]
        self.assertEqual([record[1] for record in emergency_receipts],
                         ["started", "failed"])
        self.assertTrue(emergency_receipts[-1][2]["exhausted"])

    def test_emergency_setup_failure_does_not_veto_cleanup_or_restore(self):
        backend = FakeBackend()

        def fail_after_install(budget):
            backend.recovery_budget = budget
            raise transition.TransitionError("injected recovery tool replay")

        backend.begin_emergency_recovery = fail_after_install
        controller, _clock = self.make_controller(backend)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "recovery tool replay"):
            controller.run()
        self.assertIn("cleanup_candidate", backend.calls)
        self.assertIn("restore_old", backend.calls)
        self.assertIn("publish_audit_binding", backend.calls)

    def test_emergency_and_scoring_failed_receipt_errors_do_not_veto_restore(self):
        cases = ("emergency_recovery", "scoring_deadline")
        for phase in cases:
            with self.subTest(phase=phase):
                clock = FakeClock()
                backend = FakeBackend()
                if phase == "emergency_recovery":
                    def fail_after_install(budget):
                        backend.recovery_budget = budget
                        raise transition.TransitionError(
                            "injected emergency setup")
                    backend.begin_emergency_recovery = fail_after_install
                else:
                    def expire_after_stop(name, value=None):
                        result = FakeBackend._call(backend, name, value)
                        if name == "stop_old":
                            clock.value += 600.0
                        return result
                    backend._call = expire_after_stop
                controller, _clock = self.make_controller(
                    backend, FakeJournal(fail_at=(phase, "failed")), clock)
                with self.assertRaises(transition.TransitionError):
                    controller.run()
                self.assertIn("cleanup_candidate", backend.calls)
                self.assertIn("restore_old", backend.calls)

    def test_live_backend_installs_emergency_budget_before_tool_replay(self):
        backend = object.__new__(transition.LiveBackend)
        backend.recovery_budget = None
        backend._verify_tools = mock.Mock(side_effect=
            transition.TransitionError("injected tool mismatch"))
        budget = transition.EmergencyRecoveryBudget(60.0)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "tool mismatch"):
            transition.LiveBackend.begin_emergency_recovery(backend, budget)
        self.assertIs(backend.recovery_budget, budget)

    def test_recovery_waits_never_use_the_expired_scoring_deadline(self):
        for method in (
                transition.LiveBackend._quiesce_old_for_recovery,
                transition.LiveBackend._launch_old):
            source = inspect.getsource(method)
            self.assertNotIn("deadline.absolute", source)
            self.assertIn("_recovery_wait_deadline", source)
        clock = FakeClock()
        scoring = transition.Deadline(10.0, 2.0, clock=clock)
        clock.value = scoring.absolute + 25.0
        emergency = scoring.start_emergency_recovery(60.0)
        backend = object.__new__(transition.LiveBackend)
        backend.recovery_budget = emergency
        self.assertGreater(
            transition.LiveBackend._recovery_wait_deadline(backend, 5.0),
            scoring.absolute)
        clock.value = emergency.absolute + 1.0
        extended = transition.LiveBackend._recovery_wait_deadline(
            backend, 20.0, minimum_safety_wait_s=5.0)
        self.assertEqual(extended, clock.value + 5.0)
        self.assertEqual(emergency.receipt()["safety_extension_count"], 1)
        self.assertTrue(emergency.receipt()["exhausted"])

    def test_emergency_receipt_uses_one_consistent_clock_observation(self):
        values = iter((100.0, 159.0, 161.0))
        calls = []

        def clock():
            value = next(values)
            calls.append(value)
            return value

        budget = transition.EmergencyRecoveryBudget(60.0, clock=clock)
        receipt = budget.receipt()
        self.assertFalse(receipt["exhausted"])
        self.assertEqual(receipt["observed_monotonic_s"], 159.0)
        self.assertEqual(receipt["remaining_s"], 1.0)
        self.assertEqual(calls, [100.0, 159.0])

    @staticmethod
    def _absolute_deadline_call(backend, clock):
        original = FakeBackend._call.__get__(backend, FakeBackend)

        def call(name, value=None):
            result = original(name, value)
            if name == "stop_old":
                clock.value += 600.0
            return result

        return call


class PrimitiveTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        os.chmod(self.root, 0o700)
        self.uid = os.geteuid()

    def tearDown(self):
        self.temp.cleanup()

    def test_write_new_is_exclusive_sealed_and_single_link(self):
        path = self.root / "receipt.json"
        binding = transition.write_new(
            path, b"evidence\n", owner_uid=self.uid)
        self.assertEqual(binding["sha256"], hashlib.sha256(b"evidence\n").hexdigest())
        self.assertEqual(binding["mode"], 0o444)
        self.assertEqual(binding["nlink"], 1)
        with self.assertRaises(FileExistsError):
            transition.write_new(path, b"replacement\n", owner_uid=self.uid)

    def test_receipt_journal_is_exclusive_canonical_and_hash_chained(self):
        clock = FakeClock()
        deadline = transition.Deadline(540.0, 60.0, clock=clock)
        journal = transition.ReceiptJournal(
            self.root, "test-transition", self.uid, deadline)
        first = journal.record("alpha", "started", {"value": 1})
        second = journal.record("alpha", "completed", {"value": 2})
        first_path = self.root / "0000-alpha-started.json"
        second_path = self.root / "0001-alpha-completed.json"
        self.assertEqual(
            transition.load_canonical(
                first_path, transition.PHASE_SCHEMA, "first phase"),
            first)
        self.assertEqual(
            transition.load_canonical(
                second_path, transition.PHASE_SCHEMA, "second phase"),
            second)
        self.assertEqual(second["previous_receipt_sha256"],
                         transition.file_binding(
                             first_path, with_hash=True)["sha256"])
        chain = journal.replay()
        self.assertEqual(chain["count"], 2)
        self.assertEqual(chain["head_sha256"], transition.file_binding(
            second_path, with_hash=True)["sha256"])
        os.chmod(first_path, 0o644)
        first_path.write_bytes(first_path.read_bytes().replace(b'"value":1',
                                                               b'"value":9'))
        os.chmod(first_path, 0o444)
        with self.assertRaises(transition.TransitionError):
            journal.replay()
        with self.assertRaises(FileExistsError):
            transition.ReceiptJournal(
                self.root, "test-transition", self.uid, deadline).record(
                    "alpha", "started", {})

    def test_receipt_replay_rejects_fifo_before_path_read(self):
        clock = FakeClock()
        journal = transition.ReceiptJournal(
            self.root, "test-transition", self.uid,
            transition.Deadline(540.0, 60.0, clock=clock))
        fifo = self.root / "0000-alpha-started.json"
        os.mkfifo(fifo)
        journal.sequence = 1
        with self.assertRaisesRegex(
                transition.TransitionError, "not a regular file"):
            journal.replay()

    def test_mpstat_receipt_requires_three_finite_exact_pair_intervals(self):
        def raw(intervals):
            return json.dumps({
                "sysstat": {"hosts": [{"statistics": [
                    {"cpu-load": loads} for loads in intervals
                ]}]}
            }).encode("ascii")

        good = [
            [{"cpu": "60", "idle": 99.0},
             {"cpu": "124", "idle": 100.0}]
            for _ in range(3)
        ]
        self.assertEqual(
            transition.parse_mpstat_idle_receipt(raw(good), (60, 124)),
            {60: 99.0, 124: 100.0})
        malformed = {
            "short": good[:2],
            "duplicate": [good[0] + [{"cpu": "60", "idle": 100.0}],
                          good[1], good[2]],
            "missing": [[{"cpu": "60", "idle": 100.0}], good[1], good[2]],
            "nonfinite-string": [
                [{"cpu": "60", "idle": "NaN"},
                 {"cpu": "124", "idle": 100.0}], good[1], good[2]],
            "too-busy": [
                [{"cpu": "60", "idle": 96.0},
                 {"cpu": "124", "idle": 100.0}], good[1], good[2]],
        }
        for name, intervals in malformed.items():
            with self.subTest(name=name):
                with self.assertRaises(transition.TransitionError):
                    transition.parse_mpstat_idle_receipt(
                        raw(intervals), (60, 124))
        nan_token = raw(good).replace(b"99.0", b"NaN", 1)
        with self.assertRaises(transition.TransitionError):
            transition.parse_mpstat_idle_receipt(nan_token, (60, 124))

    def test_rename_noreplace_preserves_exact_binding(self):
        source = self.root / "live.csv"
        source.write_bytes(b"raw\n")
        os.chmod(source, 0o444)
        binding = transition.file_binding(source, with_hash=True)
        destination = self.root / "archive.csv"
        observed = transition.rename_noreplace(
            source, destination, binding, parent_uid=self.uid)
        self.assertEqual(observed, binding)
        self.assertFalse(source.exists())

    def test_rename_noreplace_refuses_collision_without_mutation(self):
        source = self.root / "live.csv"
        destination = self.root / "archive.csv"
        source.write_bytes(b"raw\n")
        destination.write_bytes(b"old\n")
        os.chmod(source, 0o444)
        os.chmod(destination, 0o444)
        binding = transition.file_binding(source, with_hash=True)
        with self.assertRaises(FileExistsError):
            transition.rename_noreplace(
                source, destination, binding, parent_uid=self.uid)
        self.assertEqual(source.read_bytes(), b"raw\n")
        self.assertEqual(destination.read_bytes(), b"old\n")

    def test_rename_noreplace_refuses_changed_source(self):
        source = self.root / "live.csv"
        source.write_bytes(b"raw\n")
        os.chmod(source, 0o444)
        binding = transition.file_binding(source, with_hash=True)
        os.chmod(source, 0o644)
        source.write_bytes(b"bad\n")
        os.chmod(source, 0o444)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "source binding changed"):
            transition.rename_noreplace(
                source, self.root / "archive.csv", binding,
                parent_uid=self.uid)

    def test_file_binding_fifo_is_nonblocking_and_rejected(self):
        fifo = self.root / "fifo"
        os.mkfifo(fifo)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "not a regular file"):
            transition.file_binding(fifo, with_hash=True)

    def test_root_sealed_signal_capture_holds_validated_pidfd_through_send(self):
        source = inspect.getsource(
            transition.capture_process_identity_targets)
        self.assertLess(source.index("opened[name] = opener(pid, 0)"),
                        source.index("captured[name] = capture_one(name)"))
        self.assertLess(source.index("for name in (\"child\", \"launcher\"):\n"
                                     "            if name not in opened:"),
                        source.index("pidfd_send_signal(target_fd"))
        self.assertIn("pidfd_held_from_before_target_proc_through_signal",
                      source)
        self.assertNotIn("os.kill", source)

    def test_background_launcher_allows_only_absent_or_exact_session_leader(self):
        proc_root = self.root / "proc"
        proc_root.mkdir()
        self.assertEqual(
            transition.capture_owned_session_leader(
                7001, proc_root=proc_root),
            0)
        leader = proc_root / "7001"
        leader.mkdir()
        (leader / "stat").write_bytes(b"synthetic")
        with mock.patch.object(transition, "_parse_proc_stat", return_value={
                "process_group": 7001, "session": 7001,
                "start_tick": 123}):
            self.assertEqual(
                transition.capture_owned_session_leader(
                    7001, proc_root=proc_root),
                123)
        with mock.patch.object(transition, "_parse_proc_stat", return_value={
                "process_group": 7002, "session": 7001,
                "start_tick": 123}):
            with self.assertRaisesRegex(
                    transition.TransitionError, "exact session"):
                transition.capture_owned_session_leader(
                    7001, proc_root=proc_root)

    def test_prepare_freezes_exact_bytes_and_plan_replays(self):
        plan = make_test_code_seal(self.root, output_name="dry")
        with exact_prepare_runtime(plan):
            prepared = transition._prepare_transition_source_test_model(
                plan, controller_source=plan.controller,
                candidate_source=plan.sampler, p32_source=plan.p32)
        self.assertEqual(prepared["value"]["transition_id"], plan.transition_id)
        self.assertEqual(prepared["value"]["emergency_recovery_s"], 60.0)
        self.assertEqual(
            prepared["value"]["prepare_runtime"],
            transition.prepare_runtime_contract(
                plan, prepared["value"]["tools"]))
        self.assertEqual(
            prepared["value"]["execution"], {
                "enabled": False,
                "retirement": transition.SOURCE_ONLY_RETIREMENT,
            })
        self.assertNotIn("execute_contract", prepared["value"])
        self.assertNotIn("identity_inspection", prepared["value"])
        self.assertNotIn("python_isolation", prepared["value"])
        self.assertEqual(
            set(prepared["value"]["tools"]),
            {"bash", "chmod", "env", "fuser", "install", "mkdir",
             "mpstat", "perl", "python3", "sha256sum", "stat", "sudo",
             "taskset", "timeout"})
        self.assertNotEqual(plan.root, plan.code_seal_root)
        self.assertEqual(
            prepared["value"]["candidate"]["uid"], plan.code_owner_uid)
        self.assertEqual(
            transition.directory_binding(plan.root)["uid"],
            plan.controller_uid)
        replay = transition.verify_transition_plan(plan)
        self.assertEqual(replay, prepared["value"])
        with self.assertRaisesRegex(transition.TransitionError,
                                    "plan contract changed"):
            transition.verify_transition_plan(
                replace(plan, candidate_cpu=plan.candidate_cpu + 1))
        with self.assertRaisesRegex(transition.TransitionError,
                                    "old-sampler plan changed"):
            transition.verify_transition_plan(
                replace(plan, old_csv_inode=plan.old_csv_inode + 1))
        changed_tools = json.loads(json.dumps(prepared["value"]["tools"]))
        changed_tools["env"]["binding"]["inode"] += 1
        with mock.patch.object(
                transition, "capture_tool_records",
                return_value=changed_tools):
            with self.assertRaisesRegex(
                    transition.TransitionError, "tool identity changed"):
                transition.verify_transition_plan(plan)

        payload = dict(prepared["value"])
        del payload["schema"]
        del payload["self_sha256_excluding_field"]
        payload["execution"] = {
            "enabled": True,
            "retirement": transition.SOURCE_ONLY_RETIREMENT,
        }
        tampered = transition.sealed_record(transition.PLAN_SCHEMA, payload)
        os.chmod(plan.plan_receipt, 0o644)
        plan.plan_receipt.write_bytes(transition.canonical_json(tampered))
        os.chmod(plan.plan_receipt, 0o444)
        with self.assertRaisesRegex(
                transition.TransitionError, "execution retirement changed"):
            transition.verify_transition_plan(plan)
        for forbidden in (
                "execute_contract", "identity_inspection",
                "python_isolation"):
            with self.subTest(forbidden=forbidden):
                payload = dict(prepared["value"])
                del payload["schema"]
                del payload["self_sha256_excluding_field"]
                payload[forbidden] = {"command": ["/forbidden/live"]}
                tampered = transition.sealed_record(
                    transition.PLAN_SCHEMA, payload)
                os.chmod(plan.plan_receipt, 0o644)
                plan.plan_receipt.write_bytes(
                    transition.canonical_json(tampered))
                os.chmod(plan.plan_receipt, 0o444)
                with self.assertRaisesRegex(
                        transition.TransitionError,
                        "execution retirement changed"):
                    transition.verify_transition_plan(plan)
        os.chmod(plan.plan_receipt, 0o644)
        plan.plan_receipt.write_bytes(
            transition.canonical_json(prepared["value"]))
        os.chmod(plan.plan_receipt, 0o444)

        payload = dict(prepared["value"])
        del payload["schema"]
        del payload["self_sha256_excluding_field"]
        payload["prepare_runtime"] = dict(payload["prepare_runtime"])
        payload["prepare_runtime"]["command_nul_sha256"] = "f" * 64
        tampered = transition.sealed_record(transition.PLAN_SCHEMA, payload)
        os.chmod(plan.plan_receipt, 0o644)
        plan.plan_receipt.write_bytes(transition.canonical_json(tampered))
        os.chmod(plan.plan_receipt, 0o444)
        with self.assertRaisesRegex(
                transition.TransitionError, "prepare runtime contract"):
            transition.verify_transition_plan(plan)

    def test_production_plan_splits_root_code_from_uid1000_output_custody(self):
        plan = transition.TransitionPlan(controller_sha256="a" * 64)
        self.assertEqual(
            (plan.code_owner_uid, plan.code_owner_gid), (0, 0))
        self.assertEqual(
            (plan.controller_uid, plan.controller_gid), (1000, 1000))
        self.assertNotEqual(plan.code_owner_uid, plan.controller_uid)
        self.assertEqual(plan.controller.parent, plan.code_seal_root)
        self.assertEqual(plan.sampler.parent, plan.code_seal_root)
        self.assertEqual(plan.p32.parent, plan.code_seal_root)
        self.assertEqual(plan.legacy_sampler.parent, plan.code_seal_root)
        self.assertEqual(plan.receipts.parent, plan.root)
        self.assertNotEqual(plan.root, plan.code_seal_root)

    def test_v6_namespace_retires_every_v2_through_v5_transition_destination(self):
        plan = transition.TransitionPlan(controller_sha256="a" * 64)
        self.assertEqual(
            plan.transition_id,
            "wirehair-g8iv.6.22.1-live-dry-3c291790-v6")
        self.assertEqual(plan.root, Path(
            "/dev/shm/wirehair-g8iv.6.22.1-live-dry-3c291790-v6"))
        self.assertEqual(plan.code_seal_parent, Path(
            "/dev/shm/wirehair-g8iv.6.22.1-live-dry-3c291790-v6-root-code"))
        self.assertEqual(plan.code_seal_root,
                         plan.code_seal_parent / "seal-v6")
        self.assertEqual(
            plan.root_code_stage_receipt,
            plan.code_seal_parent /
                "root-code-seal-stage-receipt.v6.json")
        self.assertEqual(
            plan.code_stage_receipt,
            plan.root / "root-code-seal-stage-receipt.v6.json")
        self.assertEqual(plan.old_archive.name,
                         "thermal.pre-g8iv.6.22.1.3c291790.v6.csv")
        self.assertEqual(
            plan.old_unclean_archive.name,
            "thermal.pre-g8iv.6.22.1.3c291790.v6.unclean.csv")
        self.assertEqual(
            plan.old_stale_pid_archive.name,
            "sampler.pre-g8iv.6.22.1.3c291790.v6.unclean.pid")
        self.assertEqual(
            plan.audit_binding.name,
            "future-audit-binding.g8iv.6.22.1.3c291790.v6.json")
        self.assertEqual(transition.TRANSITION_SCHEMA,
                         "wirehair.wh2.thermal_transition.v8")
        self.assertEqual(transition.PHASE_SCHEMA,
                         "wirehair.wh2.thermal_transition.phase.v8")
        self.assertEqual(
            transition.AUDIT_BINDING_SCHEMA,
            "wirehair.wh2.thermal_transition.audit_binding.v8")
        self.assertEqual(transition.PLAN_SCHEMA,
                         "wirehair.wh2.thermal_transition.plan.v9")
        self.assertEqual(
            transition.CODE_SEAL_STAGE_SCHEMA,
            "wirehair.wh2.thermal_transition.code_stage.v6")
        self.assertEqual(
            transition.ROOT_CODE_STAGE_RECEIPT_SCHEMA,
            "wirehair.wh2.thermal_transition.root_code_stage_receipt.v6")
        self.assertEqual(
            transition.CODE_SEAL_RECEIPT_SCHEMA,
            "wirehair.wh2.thermal_transition.code_seal.v6")
        self.assertEqual(
            transition.FD_INSPECTION_SCHEMA,
            "wirehair.wh2.thermal_transition.fd_inspection.v6")
        self.assertEqual(
            transition.IDENTITY_INSPECTION_SCHEMA,
            "wirehair.wh2.thermal_transition.identity_inspection.v6")
        self.assertEqual(
            transition.PROCESS_SIGNAL_SCHEMA,
            "wirehair.wh2.thermal_transition.process_signal.v6")
        self.assertEqual(
            transition.RECOVERY_ACTION_EVIDENCE_SCHEMA,
            "wirehair.wh2.thermal_transition.recovery_action_evidence.v6")
        replacement_hash = hashlib.sha256(
            b"\0".join(value.encode("ascii")
                        for value in plan.replacement_old_argv) +
            b"\0").hexdigest()
        self.assertEqual(
            replacement_hash,
            "c04bcd2331f6d580d4ad7f3d0d27c0ba2692046395a8c0b67b1084b57ff19c8f")
        self.assertEqual(
            plan.replacement_old_cmdline_sha256, replacement_hash)
        retired_destinations = set()
        for version in (2, 3, 4, 5):
            retired_id = (
                "wirehair-g8iv.6.22.1-live-dry-3c291790-v%d" % version)
            retired_output = Path("/dev/shm") / retired_id
            retired_parent = Path("/dev/shm") / (retired_id + "-root-code")
            retired_destinations.update({
                retired_output, retired_parent,
                retired_parent / ("seal-v%d" % version),
                retired_parent /
                    ("root-code-seal-stage-receipt.v%d.json" % version),
                retired_output /
                    ("root-code-seal-stage-receipt.v%d.json" % version),
                transition.OLD_THERMAL_DIR /
                    ("thermal.pre-g8iv.6.22.1.3c291790.v%d.csv" % version),
                transition.OLD_THERMAL_DIR /
                    ("thermal.pre-g8iv.6.22.1.3c291790.v%d.unclean.csv" %
                     version),
                transition.OLD_THERMAL_DIR /
                    ("sampler.pre-g8iv.6.22.1.3c291790.v%d.unclean.pid" %
                     version),
                transition.OLD_THERMAL_DIR /
                    ("future-audit-binding.g8iv.6.22.1.3c291790.v%d.json" %
                     version),
            })
        self.assertTrue(retired_destinations.isdisjoint({
            plan.root, plan.code_seal_parent, plan.code_seal_root,
            plan.root_code_stage_receipt, plan.code_stage_receipt,
            plan.old_archive, plan.old_unclean_archive,
            plan.old_stale_pid_archive, plan.audit_binding,
        }))

    def test_prepare_rejects_unreviewed_controller_digest_before_creating_root(self):
        plan = make_test_code_seal(
            self.root, output_name="unreviewed-dry")
        dry = plan.root
        plan = replace(plan, controller_sha256="f" * 64)
        with self.assertRaisesRegex(
                transition.TransitionError, "reviewed root code seal"):
            transition._prepare_transition_source_test_model(
                plan, controller_source=plan.controller,
                candidate_source=plan.sampler, p32_source=plan.p32)
        self.assertFalse(dry.exists())

    def test_public_prepare_and_synthetic_model_reject_production_before_mutation(self):
        plan = transition.TransitionPlan(controller_sha256="a" * 64)
        poison_names = (
            "verify_code_seal", "capture_tool_records", "_mkdir_exact",
            "write_new")
        with ExitStack() as stack:
            poisons = {
                name: stack.enter_context(mock.patch.object(
                    transition, name, side_effect=AssertionError(
                        "retired prepare reached " + name)))
                for name in poison_names
            }
            with self.assertRaisesRegex(
                    transition.TransitionError, "source/test-only"):
                transition.prepare_transition(
                    plan, controller_source=plan.controller,
                    candidate_source=plan.sampler, p32_source=plan.p32)
            with self.assertRaisesRegex(
                    transition.TransitionError, "production V6 namespace"):
                transition._prepare_transition_source_test_model(
                    plan, controller_source=plan.controller,
                    candidate_source=plan.sampler, p32_source=plan.p32)
            source_test_id = transition.SOURCE_TEST_TRANSITION_PREFIX + \
                "production-descendant-negative"
            synthetic_parent = self.root / "synthetic-code-parent"
            unsafe_plans = (
                replace(
                    plan, transition_id=source_test_id,
                    controller_uid=os.geteuid(),
                    controller_gid=os.getegid(),
                    code_owner_uid=os.geteuid(),
                    code_owner_gid=os.getegid(),
                    root=transition.DRY_ROOT / "nested-source-test-output",
                    code_seal_parent=synthetic_parent,
                    code_seal_root=synthetic_parent / "seal-v6"),
                replace(
                    plan, transition_id=source_test_id,
                    controller_uid=os.geteuid(),
                    controller_gid=os.getegid(),
                    code_owner_uid=os.geteuid(),
                    code_owner_gid=os.getegid(),
                    root=self.root / "synthetic-output",
                    code_seal_parent=
                        transition.CODE_SEAL_PARENT / "nested-source-test",
                    code_seal_root=transition.CODE_SEAL_PARENT /
                        "nested-source-test/seal-v6"),
            )
            for unsafe in unsafe_plans:
                with self.subTest(unsafe=unsafe), self.assertRaisesRegex(
                        transition.TransitionError,
                        "production V6 namespace"):
                    transition._prepare_transition_source_test_model(
                        unsafe, controller_source=unsafe.controller,
                        candidate_source=unsafe.sampler,
                        p32_source=unsafe.p32)
                with self.subTest(stage_unsafe=unsafe), self.assertRaisesRegex(
                        transition.TransitionError,
                        "production V6 namespace"):
                    transition._code_seal_stage_plan_value(
                        unsafe, controller_source=unsafe.controller,
                        candidate_source=unsafe.sampler,
                        p32_source=unsafe.p32,
                        legacy_source=unsafe.legacy_sampler)
        for poison in poisons.values():
            poison.assert_not_called()

    def test_execute_runtime_requires_exact_env_flags_and_orig_argv(self):
        plan = replace(
            transition.TransitionPlan(), root=self.root / "runtime",
            controller_sha256="a" * 64)
        environment = transition.execute_environment(plan)
        orig_argv = transition.expected_execute_orig_argv(plan)
        flags = dict(transition.EXECUTE_FLAG_CONTRACT)
        receipt = transition.verify_execute_runtime(
            plan, observed_environment=environment,
            observed_orig_argv=orig_argv, observed_flags=flags)
        self.assertEqual(receipt["command"][:2], ["/usr/bin/env", "-i"])
        self.assertIn("/usr/bin/python3.12", receipt["command"])
        self.assertEqual(receipt["sys_orig_argv"][1:4], ["-I", "-S", "-B"])
        bad_environment = dict(environment, PYTHONPATH="/attacker")
        with self.assertRaisesRegex(transition.TransitionError, "environment"):
            transition.verify_execute_runtime(
                plan, observed_environment=bad_environment,
                observed_orig_argv=orig_argv, observed_flags=flags)
        bad_flags = dict(flags, no_site=0)
        with self.assertRaisesRegex(transition.TransitionError, "sys.flags"):
            transition.verify_execute_runtime(
                plan, observed_environment=environment,
                observed_orig_argv=orig_argv, observed_flags=bad_flags)
        with self.assertRaisesRegex(transition.TransitionError, "orig_argv"):
            transition.verify_execute_runtime(
                plan, observed_environment=environment,
                observed_orig_argv=orig_argv[:2] + orig_argv[3:],
                observed_flags=flags)

    def test_prepare_command_and_hash_are_retired_nonexecutable_fixtures(self):
        plan = transition.TransitionPlan(controller_sha256="a" * 64)
        command = transition.expected_prepare_command(plan)
        self.assertEqual(command, (
            "/usr/bin/env", "-i",
            "HOME=/dev/shm/wirehair-g8iv.6.22.1-live-dry-3c291790-v6/"
            "runtime-home",
            "PATH=/usr/bin:/bin", "LANG=C", "LC_ALL=C", "TZ=UTC",
            "/usr/bin/python3.12", "-I", "-S", "-B",
            "/dev/shm/wirehair-g8iv.6.22.1-live-dry-3c291790-v6-"
            "root-code/seal-v6/wh2_thermal_sampler_transition.py",
            "--prepare-sealed-transition",
            "--expected-controller-sha256", "a" * 64,
        ))
        command_hash = hashlib.sha256(
            b"\0".join(value.encode("ascii") for value in command) +
            b"\0").hexdigest()
        self.assertEqual(
            command_hash,
            "e65d53b5697de4d379abd7b77cfbebf51e3c02743a58a807bf80b147a8734696")
        self.assertEqual(
            dict(transition.RETIRED_NONEXECUTABLE_AUTHORITY_HASHES)
                ["prepare_argv_nul_sha256"],
            command_hash)
        tools = transition.capture_tool_records()
        contract = transition.prepare_runtime_contract(plan, tools)
        self.assertEqual(contract["command"], list(command))
        self.assertEqual(contract["command_nul_sha256"], command_hash)
        self.assertEqual(
            contract["controller_interpreter"], {
                "argv_path": "/usr/bin/python3.12",
                "device": tools["python3"]["binding"]["device"],
                "inode": tools["python3"]["binding"]["inode"],
                "resolved_path": "/usr/bin/python3.12",
            })

    def test_prepare_runtime_rejections_never_create_output_root(self):
        cases = []
        for name, value in (
                ("extra_environment", "/unexpected"),
                ("pythonpath_environment", "/attacker")):
            cases.append((
                name, "environment",
                lambda plan, name=name, value=value: {
                    "environment": {
                        **transition.prepare_environment(plan),
                        ("EXTRA" if name == "extra_environment"
                         else "PYTHONPATH"): value,
                    }}))
        for flag_name in ("no_site", "isolated"):
            cases.append((
                "missing_" + flag_name, "sys.flags",
                lambda _plan, flag_name=flag_name: {
                    "flags": {
                        **transition.PREPARE_FLAG_CONTRACT,
                        flag_name: 0,
                    }}))
        cases.extend((
            ("wrong_orig_argv", "orig_argv",
             lambda plan: {
                 "orig_argv": transition.expected_prepare_orig_argv(plan)[:-1]}),
            ("wrong_interpreter", "interpreter",
             lambda _plan: {"executable": "/usr/bin/python3"}),
        ))
        for label, error, observed in cases:
            with self.subTest(case=label):
                base = self.root / label
                base.mkdir()
                plan = make_test_code_seal(base, output_name="prepare-output")
                with exact_prepare_runtime(plan, **observed(plan)), \
                        mock.patch.object(
                            transition, "_mkdir_exact",
                            side_effect=AssertionError(
                                "rejected prepare attempted output mutation")) \
                            as mkdir_exact:
                    with self.assertRaisesRegex(
                            transition.TransitionError, error):
                        transition._prepare_transition_source_test_model(
                            plan, controller_source=plan.controller,
                            candidate_source=plan.sampler,
                            p32_source=plan.p32)
                mkdir_exact.assert_not_called()
                self.assertFalse(plan.root.exists())

    def test_hardcoded_tool_contract_matches_current_exact_binaries(self):
        records = transition.capture_tool_records()
        self.assertEqual(set(records), {
            "bash", "chmod", "env", "fuser", "install", "mkdir",
            "mpstat", "perl", "python3", "sha256sum", "stat", "sudo",
            "taskset", "timeout"})
        expected = {name: (path, digest)
                    for name, path, digest in transition.TOOL_CONTRACT}
        for name, record in records.items():
            self.assertEqual(record["path"], expected[name][0])
            self.assertEqual(record["sha256"], expected[name][1])
            self.assertEqual(record["binding"]["sha256"], expected[name][1])
            self.assertEqual(record["binding"]["uid"], 0)
            self.assertFalse(record["binding"]["mode"] & 0o022)

    def test_tool_replay_rejects_same_path_binding_change(self):
        expected = transition.capture_tool_records()
        changed = json.loads(json.dumps(expected))
        changed["taskset"]["binding"]["inode"] += 1
        with mock.patch.object(
                transition, "capture_tool_records", return_value=changed):
            with self.assertRaisesRegex(
                    transition.TransitionError,
                    "sealed tool identity changed: taskset"):
                transition.verify_tool_records(expected)
        wrong_python = json.loads(json.dumps(expected))
        wrong_python["python3"]["binding"]["inode"] += 1
        with self.assertRaisesRegex(
                transition.TransitionError,
                "controller interpreter differs"):
            transition.verify_running_interpreter(wrong_python)
        with mock.patch.object(sys, "executable", "/usr/bin/python3"):
            with self.assertRaisesRegex(
                    transition.TransitionError,
                    "controller interpreter differs"):
                transition.verify_running_interpreter(
                    expected, require_exact_path=True)

    def test_v6_every_cli_mode_is_retired_before_any_dispatch(self):
        digest = "a" * 64
        exact_modes = (
            ["--print-root-code-seal-stage-plan"],
            ["--prepare-sealed-transition"],
            ["--execute-sealed-transition", transition.TRANSITION_ID,
             "--confirmation", transition.EXECUTE_CONFIRMATION],
            ["--inspect-sealed-sampler-fds", "candidate"],
            ["--inspect-sealed-process-identities", "original"],
            ["--signal-sealed-process-identities", "original"],
        )
        with redirect_stderr(io.StringIO()):
            with self.assertRaises(SystemExit):
                transition.parse_arguments(["--execute-sealed-transition",
                                            transition.TRANSITION_ID])
            with self.assertRaises(SystemExit):
                transition.parse_arguments([
                    "--stage-root-code-seal",
                    "--expected-controller-sha256", "a" * 64])
            for mode in exact_modes:
                with self.subTest(mode=mode), self.assertRaises(SystemExit):
                    transition.parse_arguments([
                        *mode, "--expected-controller-sha256", digest])
        poison_names = (
            "parse_arguments", "code_seal_stage_plan", "prepare_transition",
            "execute_transition", "execute_fd_inspection_mode",
            "execute_identity_inspection_mode",
            "execute_process_signal_mode")
        with ExitStack() as stack:
            poisons = {
                name: stack.enter_context(mock.patch.object(
                    transition, name, side_effect=AssertionError(
                        "retired V6 dispatched " + name)))
                for name in poison_names
            }
            with self.assertRaisesRegex(
                    transition.TransitionError, "source/test-only"):
                transition.main([
                    "--prepare-sealed-transition",
                    "--expected-controller-sha256", digest])
        for poison in poisons.values():
            poison.assert_not_called()

    def test_v6_direct_execute_retires_before_every_runtime_or_mutation(self):
        poison_names = (
            "verify_execute_runtime", "verify_transition_plan",
            "load_verified_p32", "ReceiptJournal", "LiveBackend", "Deadline")
        with ExitStack() as stack:
            poisons = {
                name: stack.enter_context(mock.patch.object(
                    transition, name, side_effect=AssertionError(
                        "retired transition reached " + name)))
                for name in poison_names
            }
            opened = stack.enter_context(mock.patch.object(
                transition.os, "open", side_effect=AssertionError(
                    "retired transition created a lock")))
            with self.assertRaisesRegex(
                    transition.TransitionError, "source/test-only"):
                transition.execute_transition(transition.TransitionPlan())
        for poison in poisons.values():
            poison.assert_not_called()
        opened.assert_not_called()

    def test_production_source_has_no_stale_per_attempt_v5_literal(self):
        source = Path(transition.__file__).read_text(encoding="utf-8")
        for stale in (".v" + "5", "-v" + "5", "seal-v" + "5"):
            with self.subTest(stale=stale):
                self.assertNotIn(stale, source)

    def test_every_live_entrypoint_first_statement_is_retirement_raise(self):
        tree = ast.parse(Path(transition.__file__).read_text(encoding="utf-8"))
        functions = {
            node.name: node for node in tree.body
            if isinstance(node, ast.FunctionDef)}
        for name in (
                "code_seal_stage_plan", "prepare_transition", "main",
                "execute_transition", "execute_fd_inspection_mode",
                "execute_identity_inspection_mode",
                "execute_process_signal_mode",
                "inspect_sampler_fd_provenance"):
            with self.subTest(name=name):
                first = functions[name].body[0]
                self.assertIsInstance(first, ast.Raise)
                self.assertIsInstance(first.exc, ast.Call)
                self.assertEqual(first.exc.func.id, "TransitionError")
                self.assertEqual(first.exc.args[0].id, "SOURCE_ONLY_RETIREMENT")

    def test_thermal_parent_requires_the_expected_controller_owner(self):
        paths = {
            "old_csv": self.root / "thermal.csv",
            "old_pid_file": self.root / "sampler.pid",
            "old_archive": self.root / "archive.csv",
            "old_unclean_archive": self.root / "unclean.csv",
            "old_stale_pid_archive": self.root / "stale.pid",
            "audit_binding": self.root / "audit.json",
        }
        backend = object.__new__(transition.LiveBackend)
        backend.plan = replace(
            transition.TransitionPlan(), controller_uid=self.uid + 1,
            **paths)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "trust boundary"):
            transition.LiveBackend._validate_thermal_parent(backend)
        backend.plan = replace(backend.plan, controller_uid=self.uid)
        os.chmod(self.root, 0o750)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "trust boundary"):
            transition.LiveBackend._validate_thermal_parent(backend)


class RootSealAndFdProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        os.chmod(self.root, 0o700)

    def tearDown(self):
        self.temp.cleanup()

    def test_root_code_seal_exact_roster_and_nonwritable_ancestors(self):
        plan = make_test_code_seal(self.root)
        receipt = transition.verify_code_seal(plan)
        self.assertEqual(set(receipt["files"]), {
            "candidate", "controller", "legacy", "p32"})
        self.assertTrue(all(
            value["binding"]["mode"] == 0o444 and
            value["binding"]["nlink"] == 1
            for value in receipt["files"].values()))
        self.assertEqual(receipt["parent"]["binding"]["mode"], 0o555)
        self.assertEqual(receipt["root"]["binding"]["mode"], 0o555)
        os.chmod(plan.code_seal_parent, 0o755)
        with self.assertRaisesRegex(
                transition.TransitionError, "parent binding"):
            transition.verify_code_seal(plan)
        os.chmod(plan.code_seal_parent, 0o555)
        os.chmod(plan.code_seal_root, 0o755)
        extra = plan.code_seal_root / "unreviewed.py"
        extra.write_bytes(b"unreviewed\n")
        os.chmod(extra, 0o444)
        os.chmod(plan.code_seal_root, 0o555)
        with self.assertRaisesRegex(transition.TransitionError, "roster"):
            transition.verify_code_seal(plan)

    def test_shared_anchor_direct_child_nlink_change_replays(self):
        anchor = self.root / "shared-anchor"
        anchor.mkdir()
        plan = make_test_code_seal(anchor)
        before = transition.verify_code_seal(plan)
        recorded_anchor = before["root_stage_receipt"]["value"] \
            ["directories"]["anchor"]
        observed_before = before["anchor"]["binding"]
        self.assertNotIn("nlink", observed_before)
        output_root = plan.code_anchor / "separate-uid-output-root"
        output_root.mkdir()
        after = transition.verify_code_seal(plan)
        observed_after = after["anchor"]["binding"]
        self.assertEqual(
            transition.directory_binding(plan.code_anchor)["nlink"],
            recorded_anchor["nlink"] + 1)
        self.assertEqual(
            transition._stable_anchor_binding(observed_after),
            transition._stable_anchor_binding(recorded_anchor))
        self.assertEqual(
            after["root_stage_receipt"]["value"],
            before["root_stage_receipt"]["value"])
        self.assertEqual(after, before)

        ordered_base = self.root / "real-order"
        ordered_base.mkdir()
        ordered = make_test_code_seal(ordered_base)
        ordered = replace(
            ordered, root=ordered.code_anchor / "uid-output-root")
        with exact_prepare_runtime(ordered):
            prepared = transition._prepare_transition_source_test_model(
                ordered, controller_source=ordered.controller,
                candidate_source=ordered.sampler, p32_source=ordered.p32)
        self.assertEqual(
            transition.verify_transition_plan(ordered), prepared["value"])

    def test_root_stage_receipt_tamper_and_partial_residue_fail_closed(self):
        tamper_base = self.root / "tamper"
        tamper_base.mkdir()
        tampered = make_test_code_seal(tamper_base)
        value = transition.load_canonical(
            tampered.root_code_stage_receipt,
            transition.ROOT_CODE_STAGE_RECEIPT_SCHEMA,
            "synthetic root stage receipt")
        payload = dict(value)
        del payload["schema"]
        del payload["self_sha256_excluding_field"]
        payload["authority"] = dict(payload["authority"])
        payload["authority"]["stage_program_sha256"] = "f" * 64
        replacement = transition.sealed_record(
            transition.ROOT_CODE_STAGE_RECEIPT_SCHEMA, payload)
        os.chmod(tampered.root_code_stage_receipt, 0o644)
        tampered.root_code_stage_receipt.write_bytes(
            transition.canonical_json(replacement))
        os.chmod(tampered.root_code_stage_receipt, 0o444)
        with self.assertRaisesRegex(
                transition.TransitionError, "receipt contract changed"):
            transition.verify_code_seal(tampered)

        partial_base = self.root / "partial"
        partial_base.mkdir()
        partial = make_test_code_seal(partial_base)
        os.chmod(partial.code_seal_parent, 0o755)
        partial.root_code_stage_receipt.unlink()
        os.chmod(partial.code_seal_parent, 0o555)
        with self.assertRaisesRegex(
                transition.TransitionError, "roster"):
            transition.verify_code_seal(partial)

    def test_transition_replay_rejects_equivalent_ancestor_swap(self):
        plan = make_test_code_seal(self.root)
        with exact_prepare_runtime(plan):
            transition._prepare_transition_source_test_model(
                plan, controller_source=plan.controller,
                candidate_source=plan.sampler, p32_source=plan.p32)
        original_parent = plan.code_seal_parent.with_name("displaced-code")
        plan.code_seal_parent.rename(original_parent)
        plan.code_seal_parent.mkdir()
        plan.code_seal_root.mkdir()
        source_root = original_parent / "seal-v6"
        for destination in (
                plan.sampler, plan.controller, plan.legacy_sampler, plan.p32):
            source = source_root / destination.name
            destination.write_bytes(source.read_bytes())
            os.chmod(destination, 0o444)
        os.chmod(plan.code_seal_root, 0o555)
        os.chmod(plan.code_seal_parent, 0o555)
        with self.assertRaisesRegex(
                transition.TransitionError, "roster|code-seal receipt changed"):
            transition.verify_transition_plan(plan)

    def test_synthetic_stage_model_is_pure_and_production_emitter_retired(self):
        sources = self.root / "sources"
        sources.mkdir()
        raws = {
            "candidate": b"candidate-stage\n",
            "controller": b"controller-stage\n",
            "legacy": b"legacy-stage\n",
            "p32": b"p32-stage\n",
        }
        paths = {name: sources / (name + ".py") for name in raws}
        for name, path in paths.items():
            path.write_bytes(raws[name])
            os.chmod(path, 0o444)
        test_id = transition.SOURCE_TEST_TRANSITION_PREFIX + \
            hashlib.sha256(str(self.root).encode("ascii")).hexdigest()[:16]
        stage_parent = Path("/dev/shm") / (test_id + "-root-code")
        plan = replace(
            transition.TransitionPlan(), transition_id=test_id,
            root=Path("/dev/shm") / (test_id + "-output"),
            controller_uid=os.geteuid(), controller_gid=os.getegid(),
            code_owner_uid=os.geteuid(), code_owner_gid=os.getegid(),
            code_seal_parent=stage_parent,
            code_seal_root=stage_parent / "seal-v6",
            candidate_sha256=hashlib.sha256(raws["candidate"]).hexdigest(),
            controller_sha256=hashlib.sha256(raws["controller"]).hexdigest(),
            old_source_sha256=hashlib.sha256(raws["legacy"]).hexdigest(),
            p32_sha256=hashlib.sha256(raws["p32"]).hexdigest(),
            stage_candidate_source=paths["candidate"],
            stage_controller_source=paths["controller"],
            stage_legacy_source=paths["legacy"],
            stage_p32_source=paths["p32"],
            old_source=paths["legacy"])

        # Model a swap before any hypothetical Python open.  If the mutable
        # planner attempts a source/tool read, the trap restores the good bytes
        # before its hypothetical path hash and fails the test.  The pure
        # generator never calls the trap, and its literal plan is byte-identical
        # before and after restoration.
        os.chmod(paths["candidate"], 0o644)
        paths["candidate"].write_bytes(b"attacker swap\n")

        def forbidden_python_open(*_args, **_kwargs):
            paths["candidate"].write_bytes(raws["candidate"])
            raise AssertionError("mutable Python attempted bootstrap authority")

        with mock.patch.object(
                transition, "file_binding",
                side_effect=forbidden_python_open) as opened, \
                mock.patch.object(
                    transition, "capture_tool_records",
                    side_effect=AssertionError("mutable Python read tools")), \
                mock.patch.object(
                    transition.subprocess, "run",
                    side_effect=AssertionError("mutable Python executed stage")) \
                    as run:
            with self.assertRaisesRegex(
                    transition.TransitionError, "source/test-only"):
                transition.code_seal_stage_plan(
                    plan, controller_source=paths["controller"],
                    candidate_source=paths["candidate"],
                    p32_source=paths["p32"],
                    legacy_source=paths["legacy"])
            swapped = transition._code_seal_stage_plan_value(
                plan, controller_source=paths["controller"],
                candidate_source=paths["candidate"], p32_source=paths["p32"],
                legacy_source=paths["legacy"])
        opened.assert_not_called()
        run.assert_not_called()
        paths["candidate"].write_bytes(raws["candidate"])
        os.chmod(paths["candidate"], 0o444)
        restored = transition._code_seal_stage_plan_value(
            plan, controller_source=paths["controller"],
            candidate_source=paths["candidate"], p32_source=paths["p32"],
            legacy_source=paths["legacy"])
        self.assertEqual(swapped, restored)
        stage = restored
        self.assertTrue(stage["no_live_state_or_workload"])
        self.assertEqual(stage["destinations"], {
            "parent": str(plan.code_seal_parent),
            "receipt": str(plan.root_code_stage_receipt),
            "root": str(plan.code_seal_root)})
        argv = stage["authority_argv"]
        self.assertEqual(argv[:3], ["/usr/bin/sudo", "-n", "--"])
        self.assertIn("/usr/bin/perl", argv)
        self.assertIn("/usr/bin/bash", argv)
        self.assertNotIn("/usr/bin/python3.12", argv)
        self.assertFalse(hasattr(transition, "execute_code_seal_stage"))
        self.assertEqual(
            stage["authority_argv_canonical_sha256"],
            hashlib.sha256(transition.canonical_json(argv)).hexdigest())
        self.assertEqual(
            stage["authority_argv_nul_sha256"],
            hashlib.sha256(
                b"\0".join(value.encode("ascii") for value in argv) +
                b"\0").hexdigest())
        self.assertEqual(stage["programs"], {
            "perl_opener_sha256": hashlib.sha256(
                transition.ROOT_STAGE_PERL_OPENER.encode("ascii")).hexdigest(),
            "stage_program_sha256": hashlib.sha256(
                transition.ROOT_STAGE_BASH_PROGRAM.encode("ascii")).hexdigest(),
        })
        self.assertEqual(set(stage["tools"]), set(transition.ROOT_STAGE_TOOL_NAMES))
        self.assertIn("$open_flags = 131072",
                      transition.ROOT_STAGE_PERL_OPENER)
        self.assertLess(
            transition.ROOT_STAGE_PERL_OPENER.index("sysopen"),
            transition.ROOT_STAGE_PERL_OPENER.index("stat("))
        self.assertNotIn("rm ", transition.ROOT_STAGE_BASH_PROGRAM)
        residue = transition.ROOT_STAGE_BASH_PROGRAM.index(
            "stage residue requires forensic review")
        first_mkdir = transition.ROOT_STAGE_BASH_PROGRAM.index(
            '"$mkdir_tool" --mode=0700')
        self.assertLess(residue, first_mkdir)

    def test_exact_perl_to_bash_authority_builds_complete_v6_seal(self):
        with unprivileged_root_stage_fixture("success") as (
                plan, paths, _raws):
            self.assertNotEqual(plan.code_owner_uid, 0)
            with self.assertRaisesRegex(
                    transition.TransitionError, "source/test-only"):
                transition.code_seal_stage_plan(
                    plan, controller_source=paths["controller"],
                    candidate_source=paths["candidate"],
                    p32_source=paths["p32"],
                    legacy_source=paths["legacy"])
            stage = transition._code_seal_stage_plan_value(
                plan, controller_source=paths["controller"],
                candidate_source=paths["candidate"],
                p32_source=paths["p32"],
                legacy_source=paths["legacy"])
            argv = stage["authority_argv"]
            bash_program_index = argv.index(
                transition.ROOT_STAGE_BASH_PROGRAM)
            receipt_index = argv.index(
                str(plan.root_code_stage_receipt), bash_program_index)
            self.assertEqual(
                argv[receipt_index + 1:receipt_index + 3],
                [str(os.geteuid()), str(os.getegid())])

            # The sudo prefix is intentionally not executed.  Starting at the
            # exact emitted timeout suffix still runs the literal env, Perl
            # O_NOFOLLOW opener, and Bash authority without mutable Python in
            # the child authority chain.
            timeout_index = argv.index("/usr/bin/timeout")
            result = subprocess.run(
                argv[timeout_index:], stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                check=False, timeout=10.0, close_fds=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            receipt_raw = plan.root_code_stage_receipt.read_bytes()
            self.assertEqual(result.stdout, receipt_raw)
            receipt = transition.load_canonical(
                plan.root_code_stage_receipt,
                transition.ROOT_CODE_STAGE_RECEIPT_SCHEMA,
                "unprivileged v6 stage receipt")
            self.assertEqual(receipt_raw, transition.canonical_json(receipt))
            self.assertEqual(receipt["authority"], stage["programs"])
            self.assertEqual(
                set(receipt["sources"]),
                {"candidate", "controller", "legacy", "p32"})
            self.assertEqual(
                {record["fd"] for record in receipt["sources"].values()},
                {3, 4, 5, 6})
            self.assertTrue(all(
                record["stability_observations"] == 2 and
                record["binding"]["nlink"] == 1
                for record in receipt["sources"].values()))

            replay = transition.verify_code_seal(plan)
            self.assertEqual(
                set(replay["files"]),
                {"candidate", "controller", "legacy", "p32"})
            self.assertTrue(all(
                value["binding"]["uid"] == os.geteuid() and
                value["binding"]["gid"] == os.getegid() and
                value["binding"]["mode"] == 0o444 and
                value["binding"]["nlink"] == 1
                for value in replay["files"].values()))
            self.assertEqual(
                replay["parent"]["binding"]["mode"], 0o555)
            self.assertEqual(replay["root"]["binding"]["mode"], 0o555)
            self.assertEqual(
                sorted(entry.name for entry in plan.code_seal_parent.iterdir()),
                ["root-code-seal-stage-receipt.v6.json", "seal-v6"])

    def test_exact_authority_failure_retains_terminal_partial_residue(self):
        with unprivileged_root_stage_fixture("partial") as (
                plan, paths, raws):
            stage = transition._code_seal_stage_plan_value(
                plan, controller_source=paths["controller"],
                candidate_source=paths["candidate"],
                p32_source=paths["p32"],
                legacy_source=paths["legacy"])
            argv = list(stage["authority_argv"])
            bash_program_index = argv.index(
                transition.ROOT_STAGE_BASH_PROGRAM)
            controller_destination_index = argv.index(
                str(plan.controller), bash_program_index)
            argv[controller_destination_index] = str(
                plan.code_seal_root / "wrong-controller.py")
            timeout_index = argv.index("/usr/bin/timeout")
            command = argv[timeout_index:]
            failed = subprocess.run(
                command, stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                check=False, timeout=10.0, close_fds=True)
            self.assertNotEqual(failed.returncode, 0)
            self.assertEqual(failed.stdout, b"")
            self.assertIn(b"destination roster changed", failed.stderr)
            self.assertFalse(plan.root_code_stage_receipt.exists())
            self.assertEqual(
                stat.S_IMODE(plan.code_seal_parent.stat().st_mode), 0o700)
            self.assertEqual(
                stat.S_IMODE(plan.code_seal_root.stat().st_mode), 0o700)
            self.assertEqual(
                [entry.name for entry in plan.code_seal_root.iterdir()],
                [plan.sampler.name])
            candidate = transition.file_binding(
                plan.sampler, with_hash=True)
            self.assertEqual(candidate["sha256"], hashlib.sha256(
                raws["candidate"]).hexdigest())
            self.assertEqual(candidate["mode"], 0o444)

            # A repeated synthetic invocation must hard-stop on the retained
            # parent before any attempt to reinterpret or clean the residue.
            rejected_retry = subprocess.run(
                command, stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                check=False, timeout=10.0, close_fds=True)
            self.assertNotEqual(rejected_retry.returncode, 0)
            self.assertEqual(rejected_retry.stdout, b"")
            self.assertIn(
                b"stage residue requires forensic review",
                rejected_retry.stderr)
            self.assertEqual(
                transition.file_binding(plan.sampler, with_hash=True),
                candidate)
            self.assertFalse(plan.root_code_stage_receipt.exists())

    def test_v3_unquoted_array_write_failure_is_reproduced_and_retired(self):
        retired = r'''set -Eeuo pipefail
cd "$1"
shopt -s nullglob dotglob
file_records=()
i=0
printf -v file_records[$i] '"%s":"%s"' candidate value
'''
        reproduced = subprocess.run(
            ["/usr/bin/bash", "--noprofile", "--norc", "-c", retired,
             "retired-v3-array-write", str(self.root)],
            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=False, timeout=5.0)
        self.assertNotEqual(reproduced.returncode, 0)
        self.assertIn(b"not a valid identifier", reproduced.stderr)

        corrected = retired.replace(
            "printf -v file_records[$i]",
            'printf -v "file_records[$i]"') + \
            'printf \'%s\\n\' "${file_records[0]}"\n'
        accepted = subprocess.run(
            ["/usr/bin/bash", "--noprofile", "--norc", "-c", corrected,
             "quoted-v6-array-write", str(self.root)],
            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=False, timeout=5.0)
        self.assertEqual(accepted.returncode, 0, accepted.stderr)
        self.assertEqual(accepted.stdout, b'"candidate":"value"\n')

        post_glob = transition.ROOT_STAGE_BASH_PROGRAM.split(
            "shopt -s nullglob dotglob", 1)[1]
        self.assertIn('printf -v "file_records[$i]"', post_glob)
        self.assertIn('printf -v "source_records[$i]"', post_glob)
        self.assertNotIn("printf -v file_records[$i]", post_glob)
        self.assertNotIn("printf -v source_records[$i]", post_glob)

    def test_retired_authority_hashes_are_nonexecutable_evidence_only(self):
        plan = transition.TransitionPlan(controller_sha256="a" * 64)
        expected = {
            "root_stage_perl_opener_sha256":
                "cbccd80040d8d3ea50f2404a6185b6e04836c8a652c1b03da825157a0c7a01ef",
            "root_stage_bash_program_sha256":
                "3f17a3fdb3e4dbdabcb2a8e03d27b02e33f1cc96bb330b188bef49e4b75be686",
            "root_stage_48_argv_canonical_sha256":
                "792ada77045bc7e971dfe42aee1520b0ff233bea8e6971cf5a8a2d4a64a88fcb",
            "root_stage_48_argv_nul_sha256":
                "bba1a91a967a9ef76ac32034184010a3503a62ccc1288284cbd7ace68a506113",
            "prepare_argv_nul_sha256":
                "e65d53b5697de4d379abd7b77cfbebf51e3c02743a58a807bf80b147a8734696",
            "identity_helper_argv_nul_sha256":
                "6261bae1a8f92fde22d518f57ea754f26831ea18fd98f960a951a88719a035d0",
            "signal_helper_argv_nul_sha256":
                "50ab01a730795d53a877ac9da1dac996d59bda27d90329e5d3c3f8c1ad936d7d",
        }
        self.assertEqual(
            dict(transition.RETIRED_NONEXECUTABLE_AUTHORITY_HASHES),
            expected)
        self.assertTrue(all(
            name.endswith("sha256") and
            transition.SHA256_RE.fullmatch(digest) is not None
            for name, digest in
            transition.RETIRED_NONEXECUTABLE_AUTHORITY_HASHES))
        with self.assertRaisesRegex(
                transition.TransitionError, "source/test-only"):
            transition.code_seal_stage_plan(
                plan, controller_source=transition.REPO_CONTROLLER)
        with self.assertRaisesRegex(
                transition.TransitionError, "production V6 namespace"):
            transition._code_seal_stage_plan_value(
                plan, controller_source=transition.REPO_CONTROLLER)

    def test_root_stage_timeout_grammar_executes_env_and_rejects_v2_form(self):
        # Historical grammar is tested as data only; no production authority
        # builder or stage destination is invoked.
        timeout_prefix = [
            "/usr/bin/timeout", "--signal=KILL", "--kill-after=1.000s",
            "30.000s",
        ]
        env_prefix = [
            "/usr/bin/env", "-i", "HOME=/root", "LANG=C", "LC_ALL=C",
            "PATH=/usr/bin:/bin", "TZ=UTC",
        ]
        self.assertEqual(timeout_prefix, [
            "/usr/bin/timeout", "--signal=KILL", "--kill-after=1.000s",
            "30.000s",
        ])
        self.assertEqual(env_prefix, [
            "/usr/bin/env", "-i", "HOME=/root", "LANG=C", "LC_ALL=C",
            "PATH=/usr/bin:/bin", "TZ=UTC",
        ])
        harmless = [
            *timeout_prefix, *env_prefix,
            "/usr/bin/perl", "-e", "print qq(timeout-env-ok\\n)",
        ]
        accepted = subprocess.run(
            harmless, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=False, timeout=5.0)
        self.assertEqual(accepted.returncode, 0, accepted.stderr)
        self.assertEqual(
            (accepted.stdout, accepted.stderr), (b"timeout-env-ok\n", b""))

        retired_v2_grammar = [
            *timeout_prefix, "--", *env_prefix,
            "/usr/bin/perl", "-e", "print qq(timeout-env-ok\\n)",
        ]
        rejected = subprocess.run(
            retired_v2_grammar, stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            check=False, timeout=5.0)
        self.assertEqual(rejected.returncode, 127, rejected.stderr)
        self.assertEqual(rejected.stdout, b"")
        self.assertIn(b"failed to run command", rejected.stderr)
        self.assertIn(b"--", rejected.stderr)

    def test_perl_opener_nofollow_and_exec_fd_survival(self):
        syntax = subprocess.run(
            ["/usr/bin/bash", "--noprofile", "--norc", "-n", "-c",
             transition.ROOT_STAGE_BASH_PROGRAM],
            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=False, timeout=5.0)
        self.assertEqual(syntax.returncode, 0, syntax.stderr)
        raws = [b"source-%d\n" % index for index in range(4)]
        paths = []
        for index, raw in enumerate(raws):
            path = self.root / ("source-%d.py" % index)
            path.write_bytes(raw)
            paths.append(path)
        probe = r'''set -Eeuo pipefail
[[ $# -eq 4 ]]
for fd in "$@"; do
    /usr/bin/stat -L -c '%d:%i:%s' -- "/proc/self/fd/$fd"
    /usr/bin/sha256sum -- "/proc/self/fd/$fd"
done
'''
        command = [
            "/usr/bin/perl", "-e", transition.ROOT_STAGE_PERL_OPENER,
            *(str(path) for path in paths),
            "/usr/bin/bash", "--noprofile", "--norc", "-c", probe,
            "fd-survival-probe",
        ]
        result = subprocess.run(
            command, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=False, timeout=5.0,
            close_fds=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = result.stdout.decode("ascii").splitlines()
        self.assertEqual(len(lines), 8)
        for index, (path, raw) in enumerate(zip(paths, raws)):
            info = path.stat()
            self.assertEqual(
                lines[index * 2],
                "%d:%d:%d" % (info.st_dev, info.st_ino, info.st_size))
            self.assertEqual(
                lines[index * 2 + 1].split()[0],
                hashlib.sha256(raw).hexdigest())

        symlink = self.root / "source-link.py"
        symlink.symlink_to(paths[0])
        rejected = subprocess.run(
            [*command[:3], str(symlink), *command[4:]],
            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=False, timeout=5.0,
            close_fds=True)
        self.assertNotEqual(rejected.returncode, 0)
        self.assertIn(b"source open failed errno=40", rejected.stderr)

    def test_split_root_routes_only_frozen_namespace(self):
        split = transition.SplitCustodyRoot(
            Path("/output"), Path("/root-code/seal-v6"))
        self.assertEqual(
            split / "frozen/wirehair_expo_thermal_sampler.py",
            Path("/root-code/seal-v6/wirehair_expo_thermal_sampler.py"))
        self.assertEqual(split / "frozen", Path("/root-code/seal-v6"))
        self.assertEqual(split / "segments", Path("/output/segments"))
        self.assertEqual(split / "receipts/x.json",
                         Path("/output/receipts/x.json"))
        for child in ("", "/absolute", "../escape", "frozen/../escape"):
            with self.subTest(child=child):
                with self.assertRaises(transition.TransitionError):
                    split / child

    def test_candidate_and_legacy_fd_receipts_are_nonvacuous(self):
        candidate, candidate_contract = make_fd_snapshot(kind="candidate")
        candidate_receipt = transition.validate_sampler_fd_snapshot(
            candidate, **candidate_contract)
        transition.verify_sealed(
            candidate_receipt, transition.FD_INSPECTION_SCHEMA,
            "candidate receipt")
        self.assertTrue(
            candidate_receipt["loaded_source"]["retained_source_fd"])
        self.assertEqual(
            {value["fd"]["fd"] for value in
             candidate_receipt["i2c"].values()}, {5, 6})
        self.assertEqual(
            candidate_receipt["csv"]["fd"]["access_mode"], os.O_RDWR)
        self.assertEqual(
            candidate_receipt["csv"]["fd"]["flags"],
            next(record["flags"] for record in candidate["fds"]
                 if record["fd"] == 4))
        self.assertIn("no-continuous-host-wide",
                      candidate_receipt["ownership_scope"])

        legacy, legacy_contract = make_fd_snapshot(kind="legacy")
        legacy_receipt = transition.validate_sampler_fd_snapshot(
            legacy, **legacy_contract)
        self.assertFalse(legacy_receipt["loaded_source"]["retained_source_fd"])
        self.assertEqual(
            legacy_receipt["loaded_source"]["basis"],
            "root-custody-normal-path-and-exact-argv")
        self.assertNotIn("fd", legacy_receipt["loaded_source"])
        self.assertEqual(
            legacy_receipt["csv"]["fd"]["access_mode"], os.O_WRONLY)
        wrong_legacy_csv = copy.deepcopy(legacy)
        legacy_csv_fd = next(
            record for record in wrong_legacy_csv["fds"]
            if record["fd"] == 4)
        legacy_csv_fd["access_mode"] = os.O_RDWR
        legacy_csv_fd["flags"] = os.O_RDWR
        with self.assertRaisesRegex(
                transition.TransitionError, "growing CSV"):
            transition.validate_sampler_fd_snapshot(
                wrong_legacy_csv, **legacy_contract)
        candidate_source_fd = make_fd_snapshot(kind="candidate")[0]["fds"][0]
        legacy_with_unclaimed_source = copy.deepcopy(legacy)
        retained = copy.deepcopy(candidate_source_fd)
        retained.update({
            "device": legacy["source"]["binding"]["device"],
            "inode": legacy["source"]["binding"]["inode"],
            "link_path": legacy["source"]["path"],
            "size": legacy["source"]["binding"]["size"],
        })
        legacy_with_unclaimed_source["fds"].append(retained)
        with self.assertRaisesRegex(
                transition.TransitionError, "unclaimed"):
            transition.validate_sampler_fd_snapshot(
                legacy_with_unclaimed_source, **legacy_contract)

    def test_fd_receipt_rejects_source_csv_and_i2c_ambiguity_matrix(self):
        def remove_fd(snapshot, fd_number):
            snapshot["fds"] = [
                record for record in snapshot["fds"]
                if record["fd"] != fd_number]

        def duplicate_fd(snapshot, fd_number, new_number):
            record = copy.deepcopy(next(
                value for value in snapshot["fds"]
                if value["fd"] == fd_number))
            record["fd"] = new_number
            snapshot["fds"].append(record)

        def mutate_fd(snapshot, fd_number, key, value):
            next(record for record in snapshot["fds"]
                 if record["fd"] == fd_number)[key] = value

        def mutate_access(snapshot, fd_number, value):
            record = next(record for record in snapshot["fds"]
                          if record["fd"] == fd_number)
            record["access_mode"] = value
            record["flags"] = (record["flags"] & ~os.O_ACCMODE) | value

        mutations = {
            "missing-source": lambda value: remove_fd(value, 3),
            "duplicate-source": lambda value: duplicate_fd(value, 3, 7),
            "wrong-source-flags": lambda value: mutate_access(
                value, 3, os.O_WRONLY),
            "missing-source-nonblock": lambda value: mutate_fd(
                value, 3, "flags", value["fds"][0]["flags"] & ~os.O_NONBLOCK),
            "missing-csv": lambda value: remove_fd(value, 4),
            "duplicate-csv": lambda value: duplicate_fd(value, 4, 7),
            "wrong-csv-flags": lambda value: mutate_access(
                value, 4, os.O_RDONLY),
            "append-csv-flags": lambda value: mutate_fd(
                value, 4, "flags", value["fds"][1]["flags"] | os.O_APPEND),
            "csv-path-inode": lambda value: value["csv"]["binding"].__setitem__(
                "inode", 999),
            "missing-i2c": lambda value: remove_fd(value, 5),
            "duplicate-i2c": lambda value: duplicate_fd(value, 5, 7),
            "wrong-i2c-flags": lambda value: mutate_access(
                value, 5, os.O_RDONLY),
            "nonblocking-i2c-flags": lambda value: mutate_fd(
                value, 5, "flags", value["fds"][2]["flags"] | os.O_NONBLOCK),
            "wrong-i2c-path-inode": lambda value: mutate_fd(
                value, 5, "inode", 999),
            "wrong-i2c-rdev": lambda value: mutate_fd(
                value, 5, "rdev", os.makedev(89, 2)),
            "wrong-i2c-major": lambda value: mutate_fd(
                value, 5, "major", 90),
            "wrong-i2c-minor": lambda value: mutate_fd(
                value, 5, "minor", 2),
            "additional-i2c": lambda value: duplicate_fd(value, 5, 7),
            "unsealed-source": lambda value: value["source"]["binding"].__setitem__(
                "mode", 0o644),
        }
        for name, mutate in mutations.items():
            snapshot, contract = make_fd_snapshot(kind="candidate")
            if name == "additional-i2c":
                mutate(snapshot)
                extra = snapshot["fds"][-1]
                extra.update({
                    "link_path": "/dev/i2c-3", "minor": 3,
                    "rdev": os.makedev(89, 3), "inode": 223})
            else:
                mutate(snapshot)
            with self.subTest(name=name):
                with self.assertRaises(transition.TransitionError):
                    transition.validate_sampler_fd_snapshot(
                        snapshot, **contract)

    def test_fd_receipt_rejects_pid_start_and_boot_churn(self):
        for field in ("pid", "start_tick"):
            snapshot, contract = make_fd_snapshot(kind="candidate")
            snapshot["identity"][field] += 1
            with self.subTest(field=field):
                with self.assertRaisesRegex(
                        transition.TransitionError, "identity changed"):
                    transition.validate_sampler_fd_snapshot(
                        snapshot, **contract)
        snapshot, contract = make_fd_snapshot(kind="candidate")
        snapshot["boot_id"] = "00000000-0000-0000-0000-000000000002"
        with self.assertRaisesRegex(
                transition.TransitionError, "identity changed"):
            transition.validate_sampler_fd_snapshot(snapshot, **contract)

    def test_csv_fd_replay_requires_same_inode_and_strict_growth(self):
        snapshot, contract = make_fd_snapshot(kind="candidate")
        first_receipt = transition.validate_sampler_fd_snapshot(
            snapshot, **contract)
        first = {"fd_receipt": first_receipt}
        second = copy.deepcopy(first)
        second["fd_receipt"]["csv"]["binding"]["size"] += 1
        second["fd_receipt"]["csv"]["fd"]["size"] += 1
        transition.LiveBackend._replay_fd_growth(
            first, second, description="synthetic")
        with self.assertRaisesRegex(
                transition.TransitionError, "did not remain exact and grow"):
            transition.LiveBackend._replay_fd_growth(
                first, first, description="synthetic")
        changed = copy.deepcopy(second)
        changed["fd_receipt"]["csv"]["binding"]["inode"] += 1
        with self.assertRaisesRegex(
                transition.TransitionError, "did not remain exact and grow"):
            transition.LiveBackend._replay_fd_growth(
                first, changed, description="synthetic")


class PrivilegedIdentityInspectionTests(unittest.TestCase):
    def setUp(self):
        self.plan = transition.TransitionPlan(controller_sha256="a" * 64)
        self.tools = transition.capture_tool_records()
        request = make_identity_request(self.plan)
        contracts = transition.identity_target_contracts(
            self.plan, request, self.tools)
        stats = {
            contract["pid"]: transition._identity_proc_stat(
                identity_from_contract(contract))
            for contract in contracts.values()
        }
        self.proc_stat_patch = mock.patch.object(
            transition, "_read_proc_stat_at",
            side_effect=lambda _root, pid: copy.deepcopy(stats[pid]))
        self.proc_stat_patch.start()
        self.addCleanup(self.proc_stat_patch.stop)

    def direct_entry_plan(self):
        controller = Path(transition.__file__).resolve()
        return replace(
            self.plan, code_seal_root=controller.parent,
            controller_sha256=hashlib.sha256(
                controller.read_bytes()).hexdigest())

    def test_exact_helper_command_hash_and_plan_authority_contract(self):
        request = make_identity_request(self.plan)
        command = transition.identity_inspection_command(
            self.plan, request, self.tools)
        self.assertEqual(command, (
            "/usr/bin/env", "-i",
            "HOME=/dev/shm/wirehair-g8iv.6.22.1-live-dry-3c291790-v6/"
            "runtime-home",
            "PATH=/usr/bin:/bin", "LC_ALL=C", "LANG=C", "TZ=UTC",
            "PYTHONDONTWRITEBYTECODE=1",
            "/usr/bin/python3.12", "-I", "-S", "-B",
            "/dev/shm/wirehair-g8iv.6.22.1-live-dry-3c291790-v6-"
            "root-code/seal-v6/wh2_thermal_sampler_transition.py",
            "--inspect-sealed-process-identities", "original",
            "--identity-child-pid", "3320493",
            "--identity-child-start-tick", "160912119",
            "--identity-launcher-pid", "3320490",
            "--identity-launcher-start-tick", "160912119",
            "--identity-controller-pid", "8123",
            "--identity-allowed-absence", "none",
            "--expected-controller-sha256", "a" * 64,
        ))
        command_hash = hashlib.sha256(
            b"\0".join(item.encode("ascii") for item in command) +
            b"\0").hexdigest()
        self.assertEqual(
            command_hash,
            "6261bae1a8f92fde22d518f57ea754f26831ea18fd98f960a951a88719a035d0")
        runtime = transition.identity_inspection_runtime_contract(
            self.plan, request, self.tools)
        self.assertEqual(runtime["command"], list(command))
        self.assertEqual(runtime["command_nul_sha256"], command_hash)
        authority = transition.identity_inspection_plan_contract(
            self.plan, self.tools)
        self.assertEqual(authority["receipt_schema"],
                         transition.IDENTITY_INSPECTION_SCHEMA)
        self.assertEqual(authority["pidfd_policy"],
                         transition.IDENTITY_PIDFD_POLICY)
        self.assertEqual(authority["code_seal"], {
            "full_receipt_in_helper_runtime_and_receipt": True,
            "replayed_before_target_proc": True,
            "schema": transition.CODE_SEAL_RECEIPT_SCHEMA,
        })
        self.assertEqual(authority["absence_policies"], {
            "original": ["both", "child", "launcher", "none"],
            "replacement": ["launcher", "none"],
        })
        self.assertEqual(
            authority["controller_interpreter"],
            runtime["controller_interpreter"])
        self.assertEqual(
            authority["argv_contract"]["python_prefix"],
            list(command[8:13]))
        self.assertEqual(
            authority["process_signal"]["argv_contract"]["python_prefix"],
            list(command[8:13]))
        self.assertEqual(
            authority["process_signal"]["argv_contract"]
                     ["ordered_request_options"][-3:],
            ["--signal-target", "--signal-number",
             "--expected-controller-sha256"])
        self.assertEqual(
            authority["process_signal"]["recovery_action_evidence_schema"],
            transition.RECOVERY_ACTION_EVIDENCE_SCHEMA)

    def test_direct_live_entrypoints_retire_before_runtime_or_target_access(self):
        plan = self.direct_entry_plan()
        identity = make_identity_request(plan)
        signal_request = transition.process_signal_request(
            plan, make_identity_request(plan, allowed_absence="launcher"),
            target="child", signum=signal.SIGTERM)
        poison_names = (
            "capture_tool_records", "verify_running_interpreter",
            "verify_code_seal", "capture_process_identity_targets_under_seal",
            "capture_sampler_fd_snapshot", "sha256_file")
        with ExitStack() as stack:
            poisons = {
                name: stack.enter_context(mock.patch.object(
                    transition, name, side_effect=AssertionError(
                        "retired helper reached " + name)))
                for name in poison_names
            }
            calls = (
                lambda: transition.inspect_sampler_fd_provenance(
                    plan, kind="candidate", pid=7001, start_tick=123,
                    boot_id="boot", argv=("/synthetic",),
                    csv_path=Path("/synthetic.csv")),
                lambda: transition.execute_fd_inspection_mode(
                    plan, kind="candidate", pid=7001, start_tick=123,
                    boot_id="boot", csv_path=Path("/synthetic.csv"),
                    target_argv_hex="00"),
                lambda: transition.execute_identity_inspection_mode(
                    plan, identity),
                lambda: transition.execute_process_signal_mode(
                    plan, signal_request),
            )
            for call in calls:
                with self.subTest(call=call), self.assertRaisesRegex(
                        transition.TransitionError, "source/test-only"):
                    call()
        for poison in poisons.values():
            poison.assert_not_called()

    def test_main_rejects_exact_privileged_requests_before_dispatch(self):
        digest = "a" * 64
        fd_argv = [
            "--inspect-sealed-sampler-fds", "candidate",
            "--target-pid", "7001", "--target-start-tick", "123",
            "--target-boot-id", "boot", "--target-csv", "/synthetic.csv",
            "--target-argv-json-hex", "5b222f62696e2f74727565225d0a",
            "--expected-controller-sha256", digest,
        ]
        identity_argv = [
            "--inspect-sealed-process-identities", "original",
            "--identity-child-pid", "3320493",
            "--identity-child-start-tick", "160912119",
            "--identity-launcher-pid", "3320490",
            "--identity-launcher-start-tick", "160912119",
            "--identity-controller-pid", "8123",
            "--identity-allowed-absence", "none",
            "--expected-controller-sha256", digest,
        ]
        signal_argv = [
            "--signal-sealed-process-identities", "original",
            "--identity-child-pid", "3320493",
            "--identity-child-start-tick", "160912119",
            "--identity-launcher-pid", "3320490",
            "--identity-launcher-start-tick", "160912119",
            "--identity-controller-pid", "8123",
            "--identity-allowed-absence", "launcher",
            "--signal-target", "child", "--signal-number", "15",
            "--expected-controller-sha256", digest,
        ]
        names = (
            "execute_fd_inspection_mode", "execute_identity_inspection_mode",
            "execute_process_signal_mode", "execute_transition")
        with ExitStack() as stack:
            dispatches = {
                name: stack.enter_context(mock.patch.object(
                    transition, name, side_effect=AssertionError(
                        "retired CLI dispatched " + name)))
                for name in names
            }
            for argv in (fd_argv, identity_argv, signal_argv):
                with self.subTest(argv=argv), self.assertRaisesRegex(
                        transition.TransitionError, "source/test-only"):
                    transition.main(argv)
        for dispatch in dispatches.values():
            dispatch.assert_not_called()

    def test_code_seal_failure_precedes_every_target_process_capture(self):
        request = make_identity_request(self.plan)
        capture = mock.Mock(side_effect=AssertionError(
            "target proc capture ran after rejected code seal"))
        for defect in ("altered seal ancestor", "altered seal file metadata"):
            with self.subTest(defect=defect), \
                    mock.patch.object(
                        transition, "verify_code_seal",
                        side_effect=transition.TransitionError(defect)), \
                    mock.patch.object(
                        transition, "capture_process_identity_targets",
                        capture), \
                    self.assertRaisesRegex(
                        transition.TransitionError, defect):
                transition.capture_process_identity_targets_under_seal(
                    self.plan, request, self.tools)
        capture.assert_not_called()

    def test_zombie_is_exact_exiting_evidence_and_never_absence(self):
        request = make_identity_request(
            self.plan, allowed_absence="both")
        contracts = transition.identity_target_contracts(
            self.plan, request, self.tools)

        def identity(_root, pid):
            name = "child" if pid == self.plan.old_pid else "launcher"
            value = identity_from_contract(contracts[name])
            if name == "child":
                value["state"] = "Z"
            return value

        with tempfile.TemporaryDirectory() as name:
            boot = Path(name) / "boot_id"
            boot.write_text(self.plan.old_boot_id + "\n", encoding="ascii")
            with mock.patch.object(
                    transition, "_proc_identity_at", side_effect=identity):
                captured = transition.capture_process_identity_targets(
                    self.plan, request, self.tools,
                    proc_root=Path(name) / "proc", boot_id_path=boot,
                    pidfd_open=lambda pid, _flags: 101
                    if pid == self.plan.old_pid else 102,
                    pidfd_ready=lambda fd: fd == 101,
                    pidfd_close=lambda _fd: None,
                    proc_stat_reader=lambda _root, pid:
                        transition._identity_proc_stat(identity(_root, pid)))
        child = captured["targets"]["child"]
        self.assertEqual(child["state"], "exiting")
        self.assertEqual(child["snapshot_kind"], "pidfd-ready-with-stat")
        self.assertTrue(child["start_tick_proven"])
        self.assertNotEqual(child["state"], "absent")

        with tempfile.TemporaryDirectory() as name:
            boot = Path(name) / "boot_id"
            boot.write_text(self.plan.old_boot_id + "\n", encoding="ascii")
            with mock.patch.object(
                    transition, "_proc_identity_at", side_effect=identity):
                pre_ready = transition.capture_process_identity_targets(
                    self.plan, request, self.tools,
                    proc_root=Path(name) / "proc", boot_id_path=boot,
                    pidfd_open=lambda pid, _flags: 101
                    if pid == self.plan.old_pid else 102,
                    pidfd_ready=lambda _fd: False,
                    pidfd_close=lambda _fd: None,
                    proc_stat_reader=lambda _root, pid:
                        transition._identity_proc_stat(identity(_root, pid)))
        child = pre_ready["targets"]["child"]
        self.assertEqual(child["state"], "exiting")
        self.assertEqual(
            child["snapshot_kind"], "exit-stat-before-pidfd-ready")
        self.assertFalse(child["pidfd_ready_before_receipt"])
        self.assertTrue(child["start_tick_proven"])

        # Only a separate invocation whose fresh pidfd_open gets ESRCH may
        # advance the exact same requested PID from exiting to absent.
        attempts = []

        def exited_child(pid, _flags):
            attempts.append(pid)
            if pid == self.plan.old_pid:
                raise ProcessLookupError(errno.ESRCH, "synthetic reaped child")
            return 102

        launcher = identity_from_contract(contracts["launcher"])
        with tempfile.TemporaryDirectory() as name:
            boot = Path(name) / "boot_id"
            boot.write_text(self.plan.old_boot_id + "\n", encoding="ascii")
            with mock.patch.object(
                    transition, "_proc_identity_at",
                    return_value=launcher):
                reaped = transition.capture_process_identity_targets(
                    self.plan, request, self.tools,
                    proc_root=Path(name) / "proc", boot_id_path=boot,
                    pidfd_open=exited_child,
                    pidfd_ready=lambda _fd: False,
                    pidfd_close=lambda _fd: None)
        self.assertEqual(reaped["targets"]["child"]["state"], "absent")
        self.assertEqual(attempts.count(self.plan.old_pid), 2)

    def test_signal_mode_holds_same_pidfd_and_receipt_binds_action(self):
        identity_request = make_identity_request(
            self.plan, allowed_absence="launcher")
        request = transition.process_signal_request(
            self.plan, identity_request, target="child",
            signum=signal.SIGTERM)
        contracts = transition.identity_target_contracts(
            self.plan, identity_request, self.tools)
        identities = {
            name: identity_from_contract(contract)
            for name, contract in contracts.items()
        }
        opened = []
        closed = []
        sent = []

        def pidfd_open(pid, flags):
            self.assertEqual(flags, 0)
            opened.append(pid)
            return 101 if pid == self.plan.old_pid else 102

        def proc_identity(_root, pid):
            self.assertEqual(opened, [
                self.plan.old_pid, self.plan.old_launcher_pid])
            name = "child" if pid == self.plan.old_pid else "launcher"
            return copy.deepcopy(identities[name])

        def send(fd, signum):
            self.assertEqual(opened, [
                self.plan.old_pid, self.plan.old_launcher_pid])
            self.assertEqual(closed, [])
            sent.append((fd, signum))

        with tempfile.TemporaryDirectory() as name:
            boot = Path(name) / "boot_id"
            boot.write_text(self.plan.old_boot_id + "\n", encoding="ascii")
            with mock.patch.object(
                    transition, "_proc_identity_at",
                    side_effect=proc_identity):
                captured = transition.capture_process_identity_targets(
                    self.plan, identity_request, self.tools,
                    proc_root=Path(name) / "proc", boot_id_path=boot,
                    pidfd_open=pidfd_open,
                    pidfd_ready=lambda _fd: False,
                    pidfd_close=closed.append,
                    signal_target="child", signum=signal.SIGTERM,
                    pidfd_send_signal=send)
        self.assertEqual(sent, [(101, signal.SIGTERM)])
        self.assertEqual(closed, [101, 102])
        self.assertEqual(captured["signal_action"]["pidfd"], 101)
        code_seal = synthetic_code_seal_receipt(self.plan)
        receipt = transition.sealed_record(
            transition.PROCESS_SIGNAL_SCHEMA, {
                **captured,
                "code_seal": code_seal,
                "helper_runtime": transition._runtime_with_code_seal(
                    transition.process_signal_runtime_contract(
                        self.plan, request, self.tools), code_seal),
                "pidfd_policy": transition.PROCESS_SIGNAL_PIDFD_POLICY,
                "request": request,
                "tools": self.tools,
                "transition_id": self.plan.transition_id,
            })
        with mock.patch.object(
                transition, "verify_code_seal", return_value=code_seal):
            replay = transition.verify_process_signal_receipt(
                self.plan, receipt, request, self.tools)
        self.assertEqual(replay, receipt)
        command = transition.process_signal_command(
            self.plan, request, self.tools)
        self.assertEqual(len(command), 33)
        self.assertEqual(hashlib.sha256(
            b"\0".join(item.encode("ascii") for item in command) +
            b"\0").hexdigest(),
            "50ab01a730795d53a877ac9da1dac996d59bda27d90329e5d3c3f8c1ad936d7d")
        self.assertNotIn("-c", command)
        self.assertIn(str(self.plan.controller), command)
        self.assertEqual(command[-6:-2], (
            "--signal-target", "child", "--signal-number", "15"))

        changed = copy.deepcopy(receipt)
        del changed["schema"]
        del changed["self_sha256_excluding_field"]
        changed["signal_action"]["pidfd_send_signal_result"] = 1
        changed = transition.sealed_record(
            transition.PROCESS_SIGNAL_SCHEMA, changed)
        with mock.patch.object(
                transition, "verify_code_seal", return_value=code_seal), \
                self.assertRaisesRegex(
                    transition.TransitionError, "action receipt changed"):
            transition.verify_process_signal_receipt(
                self.plan, changed, request, self.tools)

    def test_require_absence_retries_only_exiting_until_fresh_esrch_receipt(self):
        request = make_identity_request(
            self.plan, allowed_absence="both")
        exiting = make_identity_receipt(
            self.plan, request, self.tools, exiting=("child",))
        absent = make_identity_receipt(
            self.plan, request, self.tools, absent=("child",))
        backend = object.__new__(transition.LiveBackend)
        backend.deadline = transition.Deadline(60.0, 10.0)
        backend.recovery_budget = None
        backend._inspect_original_identities = mock.Mock(
            side_effect=(exiting, absent))
        with mock.patch.object(transition.time, "sleep") as sleep:
            proof = transition.LiveBackend._require_old_pid_absent(backend)
        self.assertEqual(proof["absence_receipt"], absent)
        self.assertEqual(proof["exiting_receipts"], [exiting])
        self.assertEqual(backend._inspect_original_identities.call_count, 2)
        sleep.assert_called_once_with(0.05)

        present = make_identity_receipt(
            self.plan, request, self.tools)
        backend._inspect_original_identities = mock.Mock(
            side_effect=(present, absent))
        with mock.patch.object(transition.time, "sleep") as sleep, \
                self.assertRaisesRegex(
                    transition.TransitionError, "non-exit identity"):
            transition.LiveBackend._require_old_pid_absent(backend)
        backend._inspect_original_identities.assert_called_once()
        sleep.assert_not_called()

    def test_stop_accepts_zombie_then_esrch_and_retains_signal_receipt(self):
        request = make_identity_request(
            self.plan, allowed_absence="both")
        exiting = make_identity_receipt(
            self.plan, request, self.tools, exiting=("child",))
        child_absent = make_identity_receipt(
            self.plan, request, self.tools, absent=("child",))
        launcher_absent = make_identity_receipt(
            self.plan, request, self.tools, absent=("launcher",))
        with tempfile.TemporaryDirectory() as name:
            backend = object.__new__(transition.LiveBackend)
            backend.plan = replace(
                self.plan, old_pid_file=Path(name) / "old.pid")
            backend.tools = self.tools
            backend.deadline = transition.Deadline(60.0, 10.0)
            backend.recovery_budget = None
            backend.old_preflight = {"sealed": True}
            backend.original_signal_receipts = []
            backend._verify_tools = lambda: self.tools
            backend._require_exclusive_old_readers = lambda: {}
            backend._i2c_readers = lambda: ()
            backend._old_identity_lives = lambda: False
            backend._old_launcher_identity_lives = lambda: False
            backend._inspect_original_identities = mock.Mock(
                side_effect=(exiting, child_absent, launcher_absent))
            action = {"exact_signal_action": "SIGTERM"}

            def signal_child(signum):
                self.assertEqual(signum, signal.SIGTERM)
                backend.original_signal_receipts.append(action)
                return action

            backend._old_pidfd_signal = signal_child
            with mock.patch.object(
                    transition, "verify_running_interpreter",
                    return_value={"exact": True}), \
                    mock.patch.object(transition.time, "sleep"):
                stopped = transition.LiveBackend.stop_old(backend)
        self.assertFalse(stopped["forced"])
        self.assertEqual(stopped["signal_receipts"], [action])
        self.assertEqual(
            stopped["absence_receipts"]["child"]["exiting_receipts"],
            [exiting])

    def test_recovery_retains_preexisting_signal_after_stop_or_journal_failure(self):
        request = make_identity_request(
            self.plan, allowed_absence="both")
        exiting = make_identity_receipt(
            self.plan, request, self.tools, exiting=("child",))
        child_absent = make_identity_receipt(
            self.plan, request, self.tools, absent=("child",))
        launcher_absent = make_identity_receipt(
            self.plan, request, self.tools, absent=("launcher",))
        backend = object.__new__(transition.LiveBackend)
        backend.plan = self.plan
        backend.recovery_budget = transition.EmergencyRecoveryBudget(60.0)
        preexisting_action = {"exact_signal_action": "pre-recovery-SIGTERM"}
        backend.original_signal_receipts = [preexisting_action]
        backend._old_identity_lives = lambda: False
        backend._old_launcher_identity_lives = lambda: False
        backend._inspect_original_identities = mock.Mock(
            side_effect=(exiting, child_absent, launcher_absent))
        backend._i2c_readers = lambda: ()
        with mock.patch.object(transition.time, "sleep"):
            proof = transition.LiveBackend._quiesce_old_for_recovery(backend)
        self.assertEqual(proof["signal_receipts"], [preexisting_action])
        self.assertEqual(proof["signal_receipts_before_recovery"], 1)
        self.assertEqual(
            proof["absence_receipts"]["child"]["exiting_receipts"],
            [exiting])
        self.assertEqual(
            proof["absence_receipts"]["launcher"]["absence_receipt"],
            launcher_absent)

    def test_capture_opens_complete_pidfd_roster_before_proc_and_receipt_replays(self):
        request = make_identity_request(self.plan)
        contracts = transition.identity_target_contracts(
            self.plan, request, self.tools)
        identities = {
            name: identity_from_contract(contract)
            for name, contract in contracts.items()
        }
        opened = []
        closed = []
        ready = []
        identity_calls = []

        def pidfd_open(pid, flags):
            self.assertEqual(flags, 0)
            opened.append(pid)
            return 100 + len(opened)

        def proc_identity(_root, pid):
            self.assertEqual(opened, [
                self.plan.old_pid, self.plan.old_launcher_pid])
            identity_calls.append(pid)
            name = "child" if pid == self.plan.old_pid else "launcher"
            return copy.deepcopy(identities[name])

        with tempfile.TemporaryDirectory() as name:
            boot = Path(name) / "boot_id"
            boot.write_text(self.plan.old_boot_id + "\n", encoding="ascii")
            with mock.patch.object(
                    transition, "_proc_identity_at",
                    side_effect=proc_identity):
                captured = transition.capture_process_identity_targets(
                    self.plan, request, self.tools,
                    proc_root=Path(name) / "proc", boot_id_path=boot,
                    pidfd_open=pidfd_open,
                    pidfd_ready=lambda fd: ready.append(fd) or False,
                    pidfd_close=closed.append)
        self.assertEqual(opened, [
            self.plan.old_pid, self.plan.old_launcher_pid])
        self.assertEqual(identity_calls, [
            self.plan.old_pid, self.plan.old_launcher_pid,
            self.plan.old_pid, self.plan.old_launcher_pid])
        self.assertEqual(
            ready, [101, 101, 102, 102, 101, 101, 102, 102, 101, 102])
        self.assertEqual(closed, [101, 102])
        code_seal = synthetic_code_seal_receipt(self.plan)
        receipt = transition.sealed_record(
            transition.IDENTITY_INSPECTION_SCHEMA, {
                **captured,
                "code_seal": code_seal,
                "helper_runtime":
                    transition._runtime_with_code_seal(
                        transition.identity_inspection_runtime_contract(
                            self.plan, request, self.tools), code_seal),
                "pidfd_policy": transition.IDENTITY_PIDFD_POLICY,
                "request": request, "tools": self.tools,
                "transition_id": self.plan.transition_id,
            })
        replay = verify_synthetic_identity_receipt(
            self.plan, receipt, request, self.tools, code_seal)
        self.assertEqual(
            transition.canonical_json(replay),
            transition.canonical_json(receipt))
        self.assertTrue(all(
            replay["targets"][target]["pidfd_opened_before_target_proc"]
            for target in ("child", "launcher")))

    def test_launcher_absence_is_explicit_and_eacces_is_never_absence(self):
        request = make_identity_request(
            self.plan, allowed_absence="launcher")
        child_contract = transition.identity_target_contracts(
            self.plan, request, self.tools)["child"]
        attempts = []

        def absent_launcher(pid, _flags):
            attempts.append(pid)
            if pid == self.plan.old_launcher_pid:
                raise ProcessLookupError(errno.ESRCH, "synthetic vanished launcher")
            return 101

        with tempfile.TemporaryDirectory() as name:
            boot = Path(name) / "boot_id"
            boot.write_text(self.plan.old_boot_id + "\n", encoding="ascii")
            with mock.patch.object(
                    transition, "_proc_identity_at",
                    return_value=identity_from_contract(child_contract)):
                captured = transition.capture_process_identity_targets(
                    self.plan, request, self.tools,
                    proc_root=Path(name) / "proc", boot_id_path=boot,
                    pidfd_open=absent_launcher,
                    pidfd_ready=lambda _fd: False,
                    pidfd_close=lambda _fd: None)
            self.assertEqual(captured["targets"]["launcher"], {
                "expected_pid": self.plan.old_launcher_pid,
                "expected_start_tick": self.plan.old_launcher_start_tick,
                "pidfd_absence_observations": 2,
                "state": "absent",
            })
            self.assertEqual(attempts.count(self.plan.old_launcher_pid), 2)

            denied_attempts = []

            def denied_launcher(pid, _flags):
                denied_attempts.append(pid)
                if pid == self.plan.old_launcher_pid:
                    raise PermissionError(errno.EACCES, "synthetic EACCES")
                return 102

            with self.assertRaisesRegex(
                    transition.TransitionError,
                    "pidfd open failed without allowed absence"):
                transition.capture_process_identity_targets(
                    self.plan, request, self.tools,
                    proc_root=Path(name) / "proc", boot_id_path=boot,
                    pidfd_open=denied_launcher,
                    pidfd_ready=lambda _fd: False,
                    pidfd_close=lambda _fd: None)
            self.assertEqual(
                denied_attempts,
                [self.plan.old_pid, self.plan.old_launcher_pid])

            required = make_identity_request(
                self.plan, allowed_absence="none")
            with self.assertRaisesRegex(
                    transition.TransitionError,
                    "pidfd open failed without allowed absence"):
                transition.capture_process_identity_targets(
                    self.plan, required, self.tools,
                    proc_root=Path(name) / "proc", boot_id_path=boot,
                    pidfd_open=absent_launcher,
                    pidfd_ready=lambda _fd: False,
                    pidfd_close=lambda _fd: None)

    def test_pidfd_readiness_is_tagged_exiting_and_identity_churn_is_terminal(self):
        request = make_identity_request(self.plan)
        contracts = transition.identity_target_contracts(
            self.plan, request, self.tools)
        identities = {
            name: identity_from_contract(contract)
            for name, contract in contracts.items()
        }
        pid_to_name = {
            contracts[name]["pid"]: name for name in contracts}

        def identity(_root, pid):
            return copy.deepcopy(identities[pid_to_name[pid]])

        with tempfile.TemporaryDirectory() as name:
            boot = Path(name) / "boot_id"
            boot.write_text(self.plan.old_boot_id + "\n", encoding="ascii")
            closed = []
            with mock.patch.object(
                    transition, "_proc_identity_at", side_effect=identity):
                captured = transition.capture_process_identity_targets(
                    self.plan, request, self.tools,
                    proc_root=Path(name) / "proc", boot_id_path=boot,
                    pidfd_open=lambda pid, _flags: 101
                    if pid == self.plan.old_pid else 102,
                    pidfd_ready=lambda fd: fd == 101,
                    pidfd_close=closed.append)
            self.assertEqual(closed, [101, 102])
            child = captured["targets"]["child"]
            self.assertEqual(child["state"], "exiting")
            self.assertEqual(
                child["snapshot_kind"], "pidfd-ready-with-stat")
            self.assertTrue(child["start_tick_proven"])

            observations = {self.plan.old_pid: 0,
                            self.plan.old_launcher_pid: 0}

            def churn(_root, pid):
                observations[pid] += 1
                value = copy.deepcopy(identities[pid_to_name[pid]])
                if pid == self.plan.old_pid and observations[pid] == 2:
                    # Both PPIDs are contract-admitted, but changing between
                    # observations is still an identity race.
                    value["ppid"] = 1
                return value

            with mock.patch.object(
                    transition, "_proc_identity_at", side_effect=churn), \
                    self.assertRaisesRegex(
                        transition.TransitionError,
                        "changed during (?:identity capture|inspection)"):
                transition.capture_process_identity_targets(
                    self.plan, request, self.tools,
                    proc_root=Path(name) / "proc", boot_id_path=boot,
                    pidfd_open=lambda pid, _flags: 101
                    if pid == self.plan.old_pid else 102,
                    pidfd_ready=lambda _fd: False,
                    pidfd_close=lambda _fd: None)

    def test_scheduler_state_and_processor_are_evidence_not_identity(self):
        request = make_identity_request(self.plan)
        contracts = transition.identity_target_contracts(
            self.plan, request, self.tools)
        observations = {self.plan.old_pid: 0,
                        self.plan.old_launcher_pid: 0}

        def scheduled(_root, pid):
            name = "child" if pid == self.plan.old_pid else "launcher"
            observations[pid] += 1
            value = identity_from_contract(contracts[name])
            if observations[pid] == 2:
                value["state"] = "R"
                if name == "launcher":
                    value["processor"] += 1
            return value

        with tempfile.TemporaryDirectory() as name:
            boot = Path(name) / "boot_id"
            boot.write_text(self.plan.old_boot_id + "\n", encoding="ascii")
            with mock.patch.object(
                    transition, "_proc_identity_at", side_effect=scheduled):
                captured = transition.capture_process_identity_targets(
                    self.plan, request, self.tools,
                    proc_root=Path(name) / "proc", boot_id_path=boot,
                    pidfd_open=lambda pid, _flags: 101
                    if pid == self.plan.old_pid else 102,
                    pidfd_ready=lambda _fd: False,
                    pidfd_close=lambda _fd: None)
        self.assertEqual(captured["targets"]["child"]["identity"]["state"], "R")
        self.assertEqual(
            captured["targets"]["launcher"]["identity"]["processor"], 17)

    def test_replacement_profile_captures_dynamic_child_and_launcher_absence(self):
        for allowed_absence in ("none", "launcher"):
            with self.subTest(allowed_absence=allowed_absence):
                request = make_identity_request(
                    self.plan, profile="replacement",
                    child_start_tick=0,
                    allowed_absence=allowed_absence)
                contracts = transition.identity_target_contracts(
                    self.plan, request, self.tools)
                identities = {
                    name: identity_from_contract(
                        contract, captured_start_tick=456)
                    for name, contract in contracts.items()
                }
                if allowed_absence == "launcher":
                    identities["child"]["ppid"] = 1
                pid_to_name = {
                    contract["pid"]: name
                    for name, contract in contracts.items()}
                opened = []
                closed = []

                def pidfd_open(pid, flags):
                    self.assertEqual(flags, 0)
                    opened.append(pid)
                    name = pid_to_name[pid]
                    if name == "launcher" and \
                            allowed_absence == "launcher":
                        raise ProcessLookupError(
                            errno.ESRCH, "synthetic exited sudo launcher")
                    return 101 if name == "child" else 102

                def identity(_root, pid):
                    return copy.deepcopy(identities[pid_to_name[pid]])

                def proc_stat(_root, pid):
                    return transition._identity_proc_stat(
                        identities[pid_to_name[pid]])

                with tempfile.TemporaryDirectory() as name:
                    boot = Path(name) / "boot_id"
                    boot.write_text(
                        self.plan.old_boot_id + "\n", encoding="ascii")
                    with mock.patch.object(
                            transition, "_proc_identity_at",
                            side_effect=identity):
                        captured = transition.capture_process_identity_targets(
                            self.plan, request, self.tools,
                            proc_root=Path(name) / "proc",
                            boot_id_path=boot, pidfd_open=pidfd_open,
                            pidfd_ready=lambda _fd: False,
                            pidfd_close=closed.append,
                            proc_stat_reader=proc_stat)
                receipt = make_identity_receipt_from_capture(
                    self.plan, request, self.tools, captured)
                self.assertEqual(
                    receipt["targets"]["child"]["identity"]["start_tick"],
                    456)
                self.assertEqual(
                    receipt["targets"]["child"]
                           ["proc_identity_observations"], 2)
                if allowed_absence == "none":
                    self.assertEqual(
                        receipt["targets"]["launcher"]["state"], "present")
                    self.assertEqual(closed, [101, 102])
                    self.assertEqual(opened, [7002, 7001])
                else:
                    self.assertEqual(
                        receipt["targets"]["launcher"], {
                            "expected_pid": 7001,
                            "expected_start_tick": 123,
                            "pidfd_absence_observations": 2,
                            "state": "absent",
                        })
                    self.assertEqual(closed, [101])
                    self.assertEqual(opened, [7002, 7001, 7001])

    def test_replacement_roster_uses_exact_privileged_boundary_receipt(self):
        backend = object.__new__(transition.LiveBackend)
        backend.plan = self.plan
        backend.tools = self.tools
        backend.p32 = mock.Mock()
        backend.p32.sanitized_environment.return_value = {
            "SANITIZED": "synthetic"}
        backend._owned_session_members = mock.Mock(
            return_value=(7001, 7002))
        request = make_identity_request(
            self.plan, profile="replacement", child_pid=7002,
            child_start_tick=456, launcher_pid=7001,
            launcher_start_tick=123, controller_pid=os.getpid())
        receipt = make_identity_receipt(
            self.plan, request, self.tools)
        backend.p32.run_privileged_bounded.return_value = (
            0, transition.canonical_json(receipt), b"")
        command = transition.replacement_launch_command(
            self.plan, self.tools)
        restored = {
            "launcher_command": list(command),
            "launcher_command_sha256": hashlib.sha256(
                b"\0".join(value.encode("ascii") for value in command) +
                b"\0").hexdigest(),
            "launcher_identity": copy.deepcopy(
                receipt["targets"]["launcher"]["identity"]),
            "launcher_session": 7001,
            "launcher_start_tick": 123,
            "pid": 7002,
            "start_tick": 456,
        }
        with mock.patch.object(
                transition, "verify_code_seal",
                return_value=receipt["code_seal"]):
            roster = transition.LiveBackend._replacement_launcher_roster(
                backend, restored)
        self.assertEqual(roster, [
            receipt["targets"]["launcher"]["identity"],
            receipt["targets"]["child"]["identity"],
        ])
        expected_helper = transition.identity_inspection_command(
            self.plan, request, self.tools)
        backend.p32.run_privileged_bounded.assert_called_once_with(
            Path(str(self.tools["sudo"]["path"])),
            Path(str(self.tools["timeout"]["path"])), expected_helper,
            {"SANITIZED": "synthetic"}, stdout_limit=65536,
            stderr_limit=4096)
        backend.p32.sanitized_environment.assert_called_once_with(
            self.plan.root / "runtime-home", allocator=False)
        backend._owned_session_members.assert_called_once_with(7001)

        backend.p32.run_privileged_bounded.reset_mock()
        backend.p32.sanitized_environment.reset_mock()
        backend._owned_session_members.reset_mock()
        backend._owned_session_members.return_value = (7002,)
        absent_request = make_identity_request(
            self.plan, profile="replacement", child_pid=7002,
            child_start_tick=456, launcher_pid=7001,
            launcher_start_tick=123, controller_pid=os.getpid(),
            allowed_absence="launcher")
        absent_contracts = transition.identity_target_contracts(
            self.plan, absent_request, self.tools)
        absent_child = identity_from_contract(absent_contracts["child"])
        absent_child["ppid"] = 1
        absent_receipt = make_identity_receipt_from_capture(
            self.plan, absent_request, self.tools, {
                "boot_id": self.plan.old_boot_id,
                "targets": {
                    "child": {
                        "identity": absent_child,
                        "pidfd_opened_before_target_proc": True,
                        "pidfd_unready_after_replay": True,
                        "proc_identity_observations": 2,
                        "state": "present",
                    },
                    "launcher": {
                        "expected_pid": 7001,
                        "expected_start_tick": 123,
                        "pidfd_absence_observations": 2,
                        "state": "absent",
                    },
                },
            })
        backend.p32.run_privileged_bounded.return_value = (
            0, transition.canonical_json(absent_receipt), b"")
        absent_restored = {
            **restored,
            "launcher_identity": None,
        }
        with mock.patch.object(
                transition, "verify_code_seal",
                return_value=absent_receipt["code_seal"]):
            absent_roster = \
                transition.LiveBackend._replacement_launcher_roster(
                    backend, absent_restored)
        self.assertEqual(
            absent_roster,
            [absent_receipt["targets"]["child"]["identity"]])
        absent_helper = transition.identity_inspection_command(
            self.plan, absent_request, self.tools)
        backend.p32.run_privileged_bounded.assert_called_once_with(
            Path(str(self.tools["sudo"]["path"])),
            Path(str(self.tools["timeout"]["path"])), absent_helper,
            {"SANITIZED": "synthetic"}, stdout_limit=65536,
            stderr_limit=4096)
        backend.p32.sanitized_environment.assert_called_once_with(
            self.plan.root / "runtime-home", allocator=False)
        backend._owned_session_members.assert_called_once_with(7001)

    def test_exact_receipt_rejects_runtime_tools_identity_and_field_tamper(self):
        request = make_identity_request(self.plan)
        receipt = make_identity_receipt(
            self.plan, request, self.tools)
        cases = {
            "runtime": lambda value: value["helper_runtime"]["flags"].update(
                {"isolated": 0}),
            "tools": lambda value: value["tools"]["python3"]["binding"].update(
                {"inode": 1}),
            "identity": lambda value: value["targets"]["child"][
                "identity"]["executable"].update({"inode": 1}),
            "field-roster": lambda value: value["targets"]["child"][
                "identity"].update({"unbound": True}),
        }
        for label, mutate in cases.items():
            with self.subTest(case=label):
                changed = copy.deepcopy(receipt)
                del changed["self_sha256_excluding_field"]
                del changed["schema"]
                mutate(changed)
                changed = transition.sealed_record(
                    transition.IDENTITY_INSPECTION_SCHEMA, changed)
                with self.assertRaises(transition.TransitionError):
                    verify_synthetic_identity_receipt(
                        self.plan, changed, request, self.tools,
                        receipt["code_seal"])

    def test_backend_semantic_mismatch_is_single_shot_and_signal_precheck_vetoes(self):
        backend = object.__new__(transition.LiveBackend)
        backend.plan = self.plan
        backend.tools = self.tools
        backend.p32 = mock.Mock()
        backend.p32.sanitized_environment.return_value = {}
        request = make_identity_request(
            self.plan, controller_pid=os.getpid())
        receipt = make_identity_receipt(
            self.plan, request, self.tools)
        changed = copy.deepcopy(receipt)
        del changed["self_sha256_excluding_field"]
        del changed["schema"]
        changed["targets"]["child"]["identity"]["executable"]["inode"] = 1
        changed = transition.sealed_record(
            transition.IDENTITY_INSPECTION_SCHEMA, changed)
        backend.p32.run_privileged_bounded.return_value = (
            0, transition.canonical_json(changed), b"")
        with mock.patch.object(
                transition, "verify_code_seal",
                return_value=receipt["code_seal"]), \
                self.assertRaisesRegex(
                    transition.TransitionError,
                    "identity changed: executable"):
            transition.LiveBackend._inspect_original_identities(backend)
        backend.p32.run_privileged_bounded.assert_called_once()

        backend._signal_original_target = mock.Mock(side_effect=
            transition.TransitionError("injected root helper failure"))
        backend.p32.run_privileged_bounded.reset_mock()
        with self.assertRaisesRegex(
                transition.TransitionError, "root helper failure"):
            transition.LiveBackend._old_pidfd_signal(
                backend, signal.SIGTERM)
        backend._signal_original_target.assert_called_once_with(
            "child", signal.SIGTERM)
        backend.p32.run_privileged_bounded.assert_not_called()

    def test_backend_uses_one_combined_signal_helper_and_retains_receipt(self):
        backend = object.__new__(transition.LiveBackend)
        backend.plan = self.plan
        backend.tools = self.tools
        backend.p32 = mock.Mock()
        backend.p32.sanitized_environment.return_value = {}
        backend.original_signal_receipts = []
        identity_request = make_identity_request(
            self.plan, controller_pid=os.getpid(),
            allowed_absence="launcher")
        receipt = make_process_signal_receipt(
            self.plan, identity_request, self.tools, target="child")
        backend.p32.run_privileged_bounded.return_value = (
            0, transition.canonical_json(receipt), b"")
        backend._inspect_original_identities = mock.Mock(
            side_effect=AssertionError("separate pre-signal helper reopened pidfd"))
        with mock.patch.object(
                transition, "verify_code_seal",
                return_value=receipt["code_seal"]):
            observed = transition.LiveBackend._old_pidfd_signal(
                backend, signal.SIGTERM)
        self.assertEqual(observed, receipt)
        self.assertEqual(backend.original_signal_receipts, [receipt])
        backend._inspect_original_identities.assert_not_called()
        command = transition.process_signal_command(
            backend.plan,
            transition.process_signal_request(
                backend.plan, identity_request, target="child",
                signum=signal.SIGTERM),
            backend.tools)
        backend.p32.run_privileged_bounded.assert_called_once_with(
            Path(str(backend.tools["sudo"]["path"])),
            Path(str(backend.tools["timeout"]["path"])), command, {},
            stdout_limit=65536, stderr_limit=4096)
        backend.p32.sanitized_environment.assert_called_once_with(
            backend.plan.root / "runtime-home", allocator=False)
        self.assertNotIn("-c", command)
        self.assertIn("--signal-sealed-process-identities", command)

        failed = object.__new__(transition.LiveBackend)
        failed.plan = self.plan
        failed.tools = self.tools
        failed.p32 = mock.Mock()
        failed.p32.sanitized_environment.return_value = {}
        failed.original_signal_receipts = []
        failed.original_signal_failures = {}
        tampered = copy.deepcopy(receipt)
        del tampered["schema"]
        del tampered["self_sha256_excluding_field"]
        tampered["signal_action"]["pidfd_send_signal_result"] = 1
        tampered = transition.sealed_record(
            transition.PROCESS_SIGNAL_SCHEMA, tampered)
        failed.p32.run_privileged_bounded.return_value = (
            0, transition.canonical_json(tampered), b"")
        with mock.patch.object(
                transition, "verify_code_seal",
                return_value=receipt["code_seal"]), \
                self.assertRaisesRegex(
                    transition.TransitionError, "action receipt changed"):
            transition.LiveBackend._old_pidfd_signal(
                failed, signal.SIGTERM)
        self.assertEqual(failed.p32.run_privileged_bounded.call_count, 1)
        with self.assertRaisesRegex(
                transition.TransitionError,
                "prior child signal helper failure forbids identity retry"):
            transition.LiveBackend._old_pidfd_signal(
                failed, signal.SIGKILL)
        self.assertEqual(failed.p32.run_privileged_bounded.call_count, 1)
        self.assertTrue(failed.original_signal_failures["child"]["no_retry"])

        with mock.patch.object(
                transition, "verify_code_seal",
                return_value=receipt["code_seal"]):
            evidence = transition.LiveBackend.recovery_action_evidence(
                backend)
        self.assertEqual(evidence["signal_receipts"], [receipt])
        self.assertEqual(evidence["receipt_count"], 1)
        altered_log = copy.deepcopy(backend.original_signal_receipts)
        altered_log[0]["signal_action"]["pidfd_send_signal_result"] = 1
        backend.original_signal_receipts = altered_log
        with mock.patch.object(
                transition, "verify_code_seal",
                return_value=receipt["code_seal"]), \
                self.assertRaises(transition.TransitionError):
            transition.LiveBackend.recovery_action_evidence(backend)

    def test_root_helper_failure_during_arm_cannot_reach_signal(self):
        backend = object.__new__(transition.LiveBackend)
        backend.plan = self.plan
        backend.tools = self.tools
        backend.execute_runtime = transition.verify_execute_runtime(
            self.plan,
            observed_environment=transition.execute_environment(self.plan),
            observed_orig_argv=transition.expected_execute_orig_argv(
                self.plan),
            observed_flags=transition.EXECUTE_FLAG_CONTRACT)
        backend.old_preflight = {"thermal_parent": {"sealed": True}}
        backend._verify_tools = mock.Mock(return_value=self.tools)
        backend._validate_old_identity = mock.Mock(side_effect=
            transition.TransitionError("injected root identity helper failure"))
        backend._old_pidfd_signal = mock.Mock()
        backend.hard_preflight = mock.Mock(return_value={"sealed": True})
        backend.arm_recovery = lambda: \
            transition.LiveBackend.arm_recovery(backend)
        backend.stop_old = mock.Mock(side_effect=AssertionError(
            "stop_old ran after rejected recovery arm"))
        controller = transition.TransitionController(
            self.plan, backend, FakeJournal(),
            transition.Deadline(540.0, 60.0, clock=FakeClock()))
        with mock.patch.object(
                transition, "verify_running_interpreter",
                return_value={"exact": True}), \
                mock.patch.object(
                    transition, "verify_execute_runtime",
                    return_value=backend.execute_runtime):
            with self.assertRaisesRegex(
                    transition.TransitionError, "identity helper failure"):
                controller.run()
        backend.hard_preflight.assert_called_once_with()
        backend._validate_old_identity.assert_called_once_with()
        backend._old_pidfd_signal.assert_not_called()
        backend.stop_old.assert_not_called()
        self.assertFalse(controller.recovery_armed)


class ExclusiveReaderAndAuditTests(unittest.TestCase):
    def test_candidate_cleanup_before_start_allows_only_exact_old_reader(self):
        backend = object.__new__(transition.LiveBackend)
        backend.plan = transition.TransitionPlan()
        backend.candidate_owner = None
        backend._i2c_readers = lambda: (backend.plan.old_pid,)
        result = transition.LiveBackend.cleanup_candidate(backend)
        self.assertTrue(backend.candidate_cleanup_complete)
        self.assertEqual(result["i2c_readers_after"], [backend.plan.old_pid])

    def test_candidate_cleanup_allows_exact_old_reader_to_exit_between_probes(self):
        backend = object.__new__(transition.LiveBackend)
        backend.plan = transition.TransitionPlan()
        backend.candidate_owner = None
        readers = iter(((backend.plan.old_pid,), ()))
        backend._i2c_readers = lambda: next(readers)
        backend._cleanup_proof_pause = lambda: None
        result = transition.LiveBackend.cleanup_candidate(backend)
        self.assertTrue(backend.candidate_cleanup_complete)
        self.assertEqual(result["i2c_readers_after"], [])

    def test_candidate_cleanup_before_start_rejects_unknown_reader(self):
        backend = object.__new__(transition.LiveBackend)
        backend.plan = transition.TransitionPlan()
        backend.candidate_owner = None
        backend._i2c_readers = lambda: (99123,)
        with self.assertRaisesRegex(transition.TransitionError, "unowned"):
            transition.LiveBackend.cleanup_candidate(backend)

    @staticmethod
    def _owned_candidate_backend(*, kill_error=None, observations_empty=True):
        backend = object.__new__(transition.LiveBackend)
        backend.plan = transition.TransitionPlan()
        backend.design = {"test": True}
        backend.recovery_budget = None
        backend.candidate_cleanup_complete = False
        backend.p32 = mock.Mock()
        owner = mock.Mock()
        owner.process.pid = 7701
        owner.process.poll.side_effect = [0, 0]
        owner.pid = 7702
        owner.identity = {"pid": owner.pid}
        owner.csv_part = Path("/synthetic/candidate.csv")
        backend.candidate_owner = owner
        if kill_error is not None:
            backend.p32._kill_owned_sampler_session.side_effect = kill_error
        backend.p32._sampler_evidence_paths.return_value = {
            "csv": Path("/synthetic/candidate-final.csv")}
        backend.p32.process_identity_matches.return_value = \
            not observations_empty
        backend._owned_session_members = lambda _session: \
            () if observations_empty else (owner.pid,)
        backend._candidate_pid_present = lambda _pid: \
            not observations_empty
        backend._i2c_readers = lambda: \
            () if observations_empty else (owner.pid,)
        backend._fuser = lambda _path: \
            () if observations_empty else (owner.pid,)
        backend._cleanup_proof_pause = lambda: None
        return backend, owner

    def test_kill_helper_error_with_independent_empty_proof_allows_restore(self):
        helper_error = transition.TransitionError("injected kill transport error")
        backend, _owner = self._owned_candidate_backend(
            kill_error=helper_error, observations_empty=True)
        with self.assertRaisesRegex(
                transition.TransitionError,
                "independently proven empty ownership"):
            transition.LiveBackend.cleanup_candidate(backend)
        self.assertTrue(backend.candidate_cleanup_complete)

    def test_kill_helper_error_without_empty_proof_forbids_second_reader(self):
        helper_error = transition.TransitionError("injected kill transport error")
        backend, _owner = self._owned_candidate_backend(
            kill_error=helper_error, observations_empty=False)
        with self.assertRaisesRegex(
                transition.TransitionError, "independent empty-ownership proof failed"):
            transition.LiveBackend.cleanup_candidate(backend)
        self.assertFalse(backend.candidate_cleanup_complete)
        backend._quiesce_old_for_recovery = mock.Mock()
        with self.assertRaisesRegex(
                transition.TransitionError, "before candidate cleanup"):
            transition.LiveBackend.restore_old(backend, None)
        backend._quiesce_old_for_recovery.assert_not_called()

    def test_exited_candidate_launcher_still_invokes_exact_session_kill(self):
        backend, owner = self._owned_candidate_backend(
            observations_empty=True)
        result = transition.LiveBackend.cleanup_candidate(backend)
        backend.p32._kill_owned_sampler_session.assert_called_once_with(
            owner, backend.plan.root, backend.design)
        self.assertTrue(backend.candidate_cleanup_complete)
        self.assertEqual(len(result["empty_ownership_observations"]), 2)

    def test_cleanup_probe_separation_survives_expired_emergency_budget(self):
        backend, owner = self._owned_candidate_backend(
            observations_empty=True)
        clock = FakeClock()
        backend.recovery_budget = transition.EmergencyRecoveryBudget(
            1.0, clock=clock)
        clock.value += 2.0
        backend._cleanup_proof_pause = \
            transition.LiveBackend._cleanup_proof_pause.__get__(
                backend, transition.LiveBackend)
        with mock.patch.object(transition.time, "sleep") as sleep:
            result = transition.LiveBackend.cleanup_candidate(backend)
        sleep.assert_called_once_with(0.05)
        backend.p32._kill_owned_sampler_session.assert_called_once_with(
            owner, backend.plan.root, backend.design)
        self.assertEqual(len(result["empty_ownership_observations"]), 2)
        self.assertTrue(backend.candidate_cleanup_complete)

    def test_candidate_accept_replay_rejects_changed_final_artifacts(self):
        backend = object.__new__(transition.LiveBackend)
        backend.candidate_owner = object()
        backend.candidate_cleanup_complete = True
        backend.accept_candidate = lambda: {
            "raw_sha256": "b" * 64, "sample_count": 6}
        with self.assertRaisesRegex(
                transition.TransitionError, "changed after acceptance"):
            transition.LiveBackend._replay_candidate_accept(
                backend, {"raw_sha256": "a" * 64, "sample_count": 6},
                7001)

    def test_candidate_replay_allows_only_restored_old_i2c_reader(self):
        backend = object.__new__(transition.LiveBackend)
        backend.candidate_owner = object()
        backend.candidate_cleanup_complete = True
        accepted = {"raw_sha256": "a" * 64, "sample_count": 6}
        backend.accept_candidate = lambda: dict(accepted)
        base = {
            "candidate_pid_present": False,
            "csv_readers": [],
            "exact_identity_live": False,
            "i2c_readers": [7001],
            "launcher_returncode": 0,
            "session_members": [],
        }
        backend._candidate_cleanup_observation = lambda _owner: dict(base)
        replay = transition.LiveBackend._replay_candidate_accept(
            backend, accepted, 7001)
        self.assertEqual(replay["ownership"]["i2c_readers"], [7001])
        rejected_replay = transition.LiveBackend._replay_candidate_accept(
            backend, None, 7001)
        self.assertIsNone(rejected_replay["acceptance"])
        self.assertEqual(
            rejected_replay["ownership"]["i2c_readers"], [7001])
        for readers in ([], [7999], [7001, 7999]):
            with self.subTest(readers=readers):
                backend._candidate_cleanup_observation = lambda _owner, r=readers: {
                    **base, "i2c_readers": r}
                with self.assertRaisesRegex(
                        transition.TransitionError,
                        "ownership reappeared"):
                    transition.LiveBackend._replay_candidate_accept(
                        backend, accepted, 7001)

    def test_candidate_final_roster_rejects_unbound_extra_artifact(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            segments = root / "segments"
            segments.mkdir()
            paths = {
                "csv": segments / "segment0000.thermal.csv",
                "receipt": segments /
                    "segment0000.thermal.sampler-receipt.json",
                "validation": segments /
                    "segment0000.thermal.validation.jsonl",
            }
            for path in paths.values():
                path.write_bytes(b"sealed\n")
            backend = object.__new__(transition.LiveBackend)
            backend.plan = replace(transition.TransitionPlan(), root=root)
            self.assertEqual(
                set(transition.LiveBackend._candidate_artifact_roster(
                    backend, paths)), set(paths.values()))
            (segments / "unbound-extra").write_bytes(b"extra\n")
            with self.assertRaisesRegex(
                    transition.TransitionError, "roster changed"):
                transition.LiveBackend._candidate_artifact_roster(
                    backend, paths)

    def test_candidate_acceptance_replays_retained_thermal_summary(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            segments = root / "segments"
            segments.mkdir()
            paths = {
                "csv": segments / "segment0000.thermal.csv",
                "receipt": segments /
                    "segment0000.thermal.sampler-receipt.json",
                "validation": segments /
                    "segment0000.thermal.validation.jsonl",
            }
            for path in paths.values():
                path.write_bytes(b"sealed\n")
            backend = object.__new__(transition.LiveBackend)
            backend.plan = replace(transition.TransitionPlan(), root=root)
            backend.p32 = mock.Mock()
            backend.p32_root = root
            backend.design = {}
            backend.candidate_timing_start = 100.0
            backend.candidate_benchmark_end = 101.0
            backend.candidate_fd_bootstrap = {"nonvacuous": "bootstrap"}
            backend.candidate_fd_replay = {"nonvacuous": "replay"}
            backend.candidate_terminal = {
                "thermal_summary": {"cpu_max_c": 1.0}}
            backend.p32._sampler_evidence_paths.return_value = paths
            backend.p32.validate_sampler_terminal_evidence.return_value = (
                {}, b"validation\n")
            backend.p32._parse_thermal_csv.return_value = ({},)
            backend.p32._parse_thermal_validation.return_value = ({},)
            backend.p32.validate_thermal_interval.return_value = {
                "cpu_max_c": 2.0}
            with self.assertRaisesRegex(
                    transition.TransitionError, "summary did not replay"):
                transition.LiveBackend.accept_candidate(backend)

    def test_cleanup_inspection_error_never_authorizes_old_restart(self):
        backend, _owner = self._owned_candidate_backend(
            kill_error=transition.TransitionError("kill helper error"),
            observations_empty=True)
        backend._i2c_readers = mock.Mock(
            side_effect=transition.TransitionError("fuser inspection error"))
        with self.assertRaisesRegex(
                transition.TransitionError, "empty-ownership proof failed"):
            transition.LiveBackend.cleanup_candidate(backend)
        self.assertFalse(backend.candidate_cleanup_complete)

    def test_old_restart_real_p32_retries_header_and_partial_then_accepts_row(self):
        valid = frozen_legacy_thermal_bytes()
        header = valid[:valid.index(b"\n") + 1]
        incomplete = header[:-7]
        partial = valid[:-7]
        self.assertEqual(
            frozen_p32._bootstrap_thermal_csv_state(incomplete), "incomplete")
        self.assertEqual(
            frozen_p32._bootstrap_thermal_csv_state(header), "header")
        self.assertEqual(
            frozen_p32._bootstrap_thermal_csv_state(partial), "header")
        self.assertEqual(
            frozen_p32._bootstrap_thermal_csv_state(valid), "row")
        with self.assertRaisesRegex(
                frozen_p32.TimingError, "thermal CSV has no samples"):
            frozen_p32._parse_thermal_csv(header)
        with self.assertRaisesRegex(
                frozen_p32.TimingError,
                "thermal CSV is not canonical complete text"):
            frozen_p32._parse_thermal_csv(partial)
        with tempfile.TemporaryDirectory() as name, \
                real_p32_old_restart_fixture(
                    Path(name), (incomplete, header, partial, valid),
                    now_observations=(0.0, 0.1, 0.2, 0.3), success=True) as case:
            restored = transition.LiveBackend._launch_old(case.backend)
        self.assertEqual(restored["first_sample_monotonic_s"], "123.000000")
        self.assertEqual(case.classifier.call_count, 4)
        case.parser.assert_called_once_with(valid)
        self.assertEqual(case.sleep.call_args_list,
                         [mock.call(0.05), mock.call(0.05), mock.call(0.05)])
        case.cleanup.assert_not_called()
        case.backend._inspect_replacement_identities.assert_called_once()
        case.backend._inspect_sampler_fds.assert_called_once()

    def test_old_restart_real_p32_perpetual_header_times_out_and_cleans_once(self):
        valid = frozen_legacy_thermal_bytes()
        header = valid[:valid.index(b"\n") + 1]
        with tempfile.TemporaryDirectory() as name, \
                real_p32_old_restart_fixture(
                    Path(name), (header,), now_observations=(0.0, 1.0)) as case:
            with self.assertRaisesRegex(
                    transition.TransitionError,
                    "bootstrap failed: replacement old CSV has no complete "
                    "sample: header"):
                transition.LiveBackend._launch_old(case.backend)
        case.classifier.assert_called_once_with(header)
        case.parser.assert_not_called()
        case.sleep.assert_called_once_with(0.05)
        case.cleanup.assert_called_once_with(
            case.process, 123, "boot", case.plan.root, case.backend.design)
        case.backend._inspect_replacement_identities.assert_not_called()

    def test_old_restart_invalid_or_unknown_bootstrap_state_is_terminal(self):
        cases = (
            ("invalid-schema", b"not,the,frozen,header\n", None,
             "bootstrap state is terminal: invalid-row"),
            ("unknown-state", frozen_legacy_thermal_bytes(),
             lambda _raw: "future-state",
             "bootstrap state is terminal: future-state"),
            ("no-newline-nul", b"partial\0", None, "contains NUL"),
            ("no-newline-nonascii", b"partial\xff", None, "not ASCII"),
            ("no-newline-wrong-prefix", b"not-a-header", None,
             "incomplete header prefix is invalid"),
        )
        for label, raw, classifier, message in cases:
            with self.subTest(label=label), tempfile.TemporaryDirectory() as name, \
                    real_p32_old_restart_fixture(
                        Path(name), (raw,), now_observations=(0.0,),
                        classifier=classifier) as case:
                with self.assertRaisesRegex(
                        transition.TransitionError, message):
                    transition.LiveBackend._launch_old(case.backend)
                case.parser.assert_not_called()
                case.sleep.assert_not_called()
                case.cleanup.assert_called_once_with(
                    case.process, 123, "boot", case.plan.root,
                    case.backend.design)
                case.backend._inspect_replacement_identities.assert_not_called()
                if label in ("no-newline-nul", "no-newline-nonascii"):
                    case.classifier.assert_not_called()
                else:
                    case.classifier.assert_called_once_with(raw)

    def test_old_restart_complete_malformed_or_implausible_row_is_p32_terminal(self):
        cases = (
            ("malformed", "bad", "thermal utilization field missing"),
            ("implausible", "101.0",
             "thermal CPU utilization exceeds 100 percent"),
        )
        for label, busy, message in cases:
            malformed = frozen_legacy_thermal_bytes(cpu_busy_pct=busy)
            self.assertEqual(
                frozen_p32._bootstrap_thermal_csv_state(malformed), "row")
            with self.assertRaisesRegex(frozen_p32.TimingError, message):
                frozen_p32._parse_thermal_csv(malformed)
            with self.subTest(label=label), tempfile.TemporaryDirectory() as name, \
                    real_p32_old_restart_fixture(
                        Path(name), (malformed,),
                        now_observations=(0.0,)) as case:
                with self.assertRaisesRegex(frozen_p32.TimingError, message):
                    transition.LiveBackend._launch_old(case.backend)
                case.classifier.assert_called_once_with(malformed)
                case.parser.assert_called_once_with(malformed)
                case.sleep.assert_not_called()
                case.cleanup.assert_called_once_with(
                    case.process, 123, "boot", case.plan.root,
                    case.backend.design)
                case.backend._inspect_replacement_identities.assert_not_called()

    def test_old_restart_cleans_exact_session_after_unexpected_capture_error(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            source = root / "old.py"
            source.write_bytes(b"old\n")
            os.chmod(source, 0o444)
            plan = make_test_code_seal(root, legacy=b"old\n")
            plan = replace(
                plan, root=root,
                old_source=source, old_csv=root / "thermal.csv",
                old_pid_file=root / "sampler.pid",
                old_source_sha256=hashlib.sha256(b"old\n").hexdigest())
            backend = object.__new__(transition.LiveBackend)
            backend.plan = plan
            backend.tools = transition.capture_tool_records()
            backend.design = {"tools": backend.tools}
            backend.p32 = mock.Mock()
            backend.deadline = transition.Deadline(540.0, 60.0)
            backend.recovery_budget = \
                transition.EmergencyRecoveryBudget(60.0)
            backend.old_preflight = {
                "source": transition.file_binding(source, with_hash=True)}
            backend._i2c_readers = lambda: ()
            backend._environment = lambda: {"PATH": "/usr/bin:/bin"}
            backend._require_old_launcher_absent = lambda: None
            process = mock.Mock(pid=7701)
            process.poll.return_value = None
            with mock.patch.object(
                    transition.subprocess, "Popen", return_value=process), \
                    mock.patch.object(
                        transition, "capture_owned_session_leader",
                        side_effect=RuntimeError("injected capture error")):
                with self.assertRaisesRegex(RuntimeError,
                                            "injected capture error"):
                    transition.LiveBackend._launch_old(backend)
            backend.p32._kill_owned_process_session.assert_called_once_with(
                process, 0, mock.ANY, plan.root, backend.design)

    def test_old_restart_does_not_retry_failed_fd_provenance(self):
        with tempfile.TemporaryDirectory() as name:
            base = Path(name)
            plan = make_test_code_seal(base, legacy=b"old source\n")
            plan = replace(
                plan, old_csv=base / "legacy.csv",
                old_pid_file=base / "legacy.pid")
            source_binding = transition.file_binding(
                plan.legacy_sampler, with_hash=True)
            csv_binding = {
                "device": 31, "gid": 0, "inode": 310, "mode": 0o444,
                "nlink": 1, "size": 7, "uid": 0,
            }
            raw_pid = b"7702\n"
            pid_binding = {
                "device": 32, "gid": 0, "inode": 320, "mode": 0o444,
                "nlink": 1, "sha256": hashlib.sha256(raw_pid).hexdigest(),
                "size": len(raw_pid), "uid": 0,
            }
            backend = object.__new__(transition.LiveBackend)
            backend.plan = plan
            backend.tools = transition.capture_tool_records()
            backend.design = {}
            backend.recovery_budget = \
                transition.EmergencyRecoveryBudget(60.0)
            backend.p32 = mock.Mock()
            backend.p32.capture_process_identity.side_effect = AssertionError(
                "retryable readiness discovery read privileged process identity")
            helper_identity = {
                "affinity": str(plan.old_cpu),
                "argv": list(plan.replacement_old_argv),
                "cmdline_sha256": plan.replacement_old_cmdline_sha256,
                "executable": {
                    "device": backend.tools["python3"]["binding"]["device"],
                    "inode": backend.tools["python3"]["binding"]["inode"],
                },
                "pid": 7702, "ppid": 7701, "process_group": 7701,
                "processor": plan.old_cpu, "session": 7701,
                "start_tick": 456, "state": "S", "uids": [0, 0, 0, 0],
            }
            backend._inspect_replacement_identities = mock.Mock(return_value={
                "boot_id": "boot",
                "targets": {
                    "child": {"identity": helper_identity,
                              "state": "present"},
                    "launcher": {"state": "absent"},
                },
            })
            backend.p32._parse_thermal_csv.return_value = (
                {"monotonic_s": "1.000000"},)
            backend.p32._bootstrap_thermal_csv_state.return_value = "row"
            backend._require_old_launcher_absent = mock.Mock()
            backend._i2c_readers = mock.Mock(side_effect=((), (7702,)))
            backend._fuser = mock.Mock(return_value=(7702,))
            backend._replacement_launch_command = mock.Mock(
                return_value=("/usr/bin/false",))
            backend._inspect_sampler_fds = mock.Mock(side_effect=
                transition.TransitionError("injected provenance failure"))
            process = mock.Mock(pid=7701)
            process.poll.return_value = None

            def launch(*_args, **_kwargs):
                plan.old_csv.write_bytes(b"sample\n")
                plan.old_pid_file.write_bytes(raw_pid)
                return process

            def binding(path, *, with_hash):
                del with_hash
                if path == plan.legacy_sampler:
                    return source_binding
                if path == plan.old_csv:
                    return csv_binding
                if path == plan.old_pid_file:
                    return pid_binding
                raise AssertionError("unexpected binding path: %s" % path)

            code_seal = {"files": {"legacy": {"binding": source_binding}}}
            with mock.patch.object(
                    transition, "verify_code_seal", return_value=code_seal), \
                    mock.patch.object(
                        transition, "execute_environment", return_value={}), \
                    mock.patch.object(
                        transition, "file_binding", side_effect=binding), \
                    mock.patch.object(
                        transition.subprocess, "Popen", side_effect=launch), \
                    mock.patch.object(
                        transition, "capture_owned_session_leader",
                        return_value=123), \
                    mock.patch.object(
                        transition.Path, "read_text", return_value="boot\n"), \
                    mock.patch.object(transition.time, "sleep") as sleep:
                with self.assertRaisesRegex(
                        transition.TransitionError,
                        "injected provenance failure"):
                    transition.LiveBackend._launch_old(backend)
            backend._inspect_sampler_fds.assert_called_once_with(
                "legacy", pid=7702, start_tick=456, boot_id="boot",
                argv=plan.replacement_old_argv, csv_path=plan.old_csv)
            backend._inspect_replacement_identities.assert_called_once_with(
                child_pid=7702, child_start_tick=0, launcher_pid=7701,
                launcher_start_tick=123, allowed_absence="launcher")
            backend.p32.capture_process_identity.assert_not_called()
            backend.p32._kill_owned_process_session.assert_called_once_with(
                process, 123, "boot", plan.root, backend.design)
            sleep.assert_not_called()

    def test_old_restart_does_not_retry_privileged_identity_mismatch(self):
        with tempfile.TemporaryDirectory() as name:
            base = Path(name)
            plan = make_test_code_seal(base, legacy=b"old source\n")
            plan = replace(
                plan, old_csv=base / "legacy.csv",
                old_pid_file=base / "legacy.pid")
            source_binding = transition.file_binding(
                plan.legacy_sampler, with_hash=True)
            raw_pid = b"7702\n"
            csv_binding = {
                "device": 31, "gid": 0, "inode": 310, "mode": 0o444,
                "nlink": 1, "size": 7, "uid": 0,
            }
            pid_binding = {
                "device": 32, "gid": 0, "inode": 320, "mode": 0o444,
                "nlink": 1, "sha256": hashlib.sha256(raw_pid).hexdigest(),
                "size": len(raw_pid), "uid": 0,
            }
            backend = object.__new__(transition.LiveBackend)
            backend.plan = plan
            backend.tools = transition.capture_tool_records()
            backend.design = {}
            backend.recovery_budget = \
                transition.EmergencyRecoveryBudget(60.0)
            backend.p32 = mock.Mock()
            backend.p32._parse_thermal_csv.return_value = (
                {"monotonic_s": "1.000000"},)
            backend.p32._bootstrap_thermal_csv_state.return_value = "row"
            backend.p32.capture_process_identity.side_effect = AssertionError(
                "retryable discovery read privileged process identity")
            backend._require_old_launcher_absent = mock.Mock()
            backend._i2c_readers = mock.Mock(side_effect=((), (7702,)))
            backend._fuser = mock.Mock(return_value=(7702,))
            backend._replacement_launch_command = mock.Mock(
                return_value=("/usr/bin/false",))
            backend._inspect_replacement_identities = mock.Mock(side_effect=
                transition.TransitionError("injected identity mismatch"))
            backend._inspect_sampler_fds = mock.Mock()
            process = mock.Mock(pid=7701)
            process.poll.return_value = None

            def launch(*_args, **_kwargs):
                plan.old_csv.write_bytes(b"sample\n")
                plan.old_pid_file.write_bytes(raw_pid)
                return process

            def binding(path, *, with_hash):
                del with_hash
                if path == plan.legacy_sampler:
                    return source_binding
                if path == plan.old_csv:
                    return csv_binding
                if path == plan.old_pid_file:
                    return pid_binding
                raise AssertionError("unexpected binding path: %s" % path)

            with mock.patch.object(
                    transition, "verify_code_seal",
                    return_value={"files": {
                        "legacy": {"binding": source_binding}}}), \
                    mock.patch.object(
                        transition, "execute_environment", return_value={}), \
                    mock.patch.object(
                        transition, "file_binding", side_effect=binding), \
                    mock.patch.object(
                        transition.subprocess, "Popen", side_effect=launch), \
                    mock.patch.object(
                        transition, "capture_owned_session_leader",
                        return_value=123), \
                    mock.patch.object(
                        transition.Path, "read_text", return_value="boot\n"), \
                    mock.patch.object(transition.time, "sleep") as sleep:
                with self.assertRaisesRegex(
                        transition.TransitionError, "identity mismatch"):
                    transition.LiveBackend._launch_old(backend)
            backend._inspect_replacement_identities.assert_called_once_with(
                child_pid=7702, child_start_tick=0, launcher_pid=7701,
                launcher_start_tick=123, allowed_absence="launcher")
            backend.p32.capture_process_identity.assert_not_called()
            backend._inspect_sampler_fds.assert_not_called()
            backend.p32._kill_owned_process_session.assert_called_once_with(
                process, 123, "boot", plan.root, backend.design)
            sleep.assert_not_called()

    def test_recovery_quiesces_exact_original_launcher_before_restart(self):
        backend = object.__new__(transition.LiveBackend)
        backend.plan = transition.TransitionPlan()
        backend.recovery_budget = transition.EmergencyRecoveryBudget(60.0)
        backend._old_identity_lives = lambda: False
        backend._require_old_pid_absent = mock.Mock()
        launcher_liveness = iter((True, False, False))
        backend._old_launcher_identity_lives = lambda: next(
            launcher_liveness, False)
        backend._old_launcher_pidfd_signal = mock.Mock()
        backend._require_old_launcher_absent = mock.Mock()
        backend._raw_proc_pid_present = lambda _pid: False
        backend._i2c_readers = lambda: ()
        transition.LiveBackend._quiesce_old_for_recovery(backend)
        backend._old_launcher_pidfd_signal.assert_called_once_with(signal.SIGTERM)
        backend._require_old_launcher_absent.assert_called_once_with()

    def test_recovery_uses_privileged_tagged_absence_not_raw_proc_path(self):
        backend = object.__new__(transition.LiveBackend)
        backend.plan = transition.TransitionPlan()
        backend.recovery_budget = transition.EmergencyRecoveryBudget(60.0)
        backend._old_identity_lives = lambda: False
        backend._old_launcher_identity_lives = lambda: False
        backend._require_old_pid_absent = mock.Mock()
        backend._require_old_launcher_absent = mock.Mock()
        backend._raw_proc_pid_present = mock.Mock(side_effect=AssertionError(
            "recovery used unprivileged raw /proc presence"))
        backend._i2c_readers = lambda: ()
        with mock.patch.object(transition.time, "sleep") as sleep:
            transition.LiveBackend._quiesce_old_for_recovery(backend)
        sleep.assert_not_called()
        backend._raw_proc_pid_present.assert_not_called()
        backend._require_old_launcher_absent.assert_called_once_with()

    def test_old_sampler_identity_uses_privileged_helper_on_proc_exe_eacces(self):
        backend = object.__new__(transition.LiveBackend)
        backend.plan = transition.TransitionPlan(controller_sha256="a" * 64)
        backend.tools = transition.capture_tool_records()
        backend.p32 = mock.Mock()
        backend.p32.sanitized_environment.return_value = {}
        request = make_identity_request(
            backend.plan, controller_pid=os.getpid())
        receipt = make_identity_receipt(
            backend.plan, request, backend.tools)
        backend.p32.run_privileged_bounded.return_value = (
            0, transition.canonical_json(receipt), b"")
        with mock.patch.object(
                transition, "verify_code_seal",
                return_value=receipt["code_seal"]), \
                mock.patch.object(
                    transition, "_proc_identity",
                    side_effect=PermissionError(
                        13, "simulated UID1000 /proc/PID/exe EACCES")) as direct:
            observed = transition.LiveBackend._validate_old_identity(
                backend)
        self.assertEqual(
            observed["child"], receipt["targets"]["child"]["identity"])
        self.assertEqual(observed["identity_receipt"], receipt)
        command = transition.identity_inspection_command(
            backend.plan, request, backend.tools)
        backend.p32.run_privileged_bounded.assert_called_once_with(
            Path(str(backend.tools["sudo"]["path"])),
            Path(str(backend.tools["timeout"]["path"])), command, {},
            stdout_limit=65536, stderr_limit=4096)
        backend.p32.sanitized_environment.assert_called_once_with(
            backend.plan.root / "runtime-home", allocator=False)
        direct.assert_not_called()

    def test_replacement_launcher_roster_rejects_unknown_session_member(self):
        backend = object.__new__(transition.LiveBackend)
        backend.plan = transition.TransitionPlan()
        backend.tools = transition.capture_tool_records()
        backend._owned_session_members = lambda _session: (7001, 7002, 7003)
        command = transition.LiveBackend._replacement_launch_command(backend)
        with self.assertRaisesRegex(
                transition.TransitionError, "unknown member"):
            transition.LiveBackend._replacement_launcher_roster(backend, {
                "launcher_command": list(command),
                "launcher_command_sha256": hashlib.sha256(
                    b"\0".join(value.encode("ascii")
                                for value in command) + b"\0").hexdigest(),
                "launcher_session": 7001, "launcher_start_tick": 123,
                "pid": 7002, "start_tick": 456})

    def test_replacement_launcher_roster_rejects_unsealed_command(self):
        backend = object.__new__(transition.LiveBackend)
        backend.plan = transition.TransitionPlan()
        backend.tools = transition.capture_tool_records()
        backend._owned_session_members = mock.Mock(return_value=(7002,))
        command = transition.LiveBackend._replacement_launch_command(backend)
        with self.assertRaisesRegex(
                transition.TransitionError, "identity is malformed"):
            transition.LiveBackend._replacement_launcher_roster(backend, {
                "launcher_command": [*command[:-1], "--changed"],
                "launcher_command_sha256": hashlib.sha256(
                    b"\0".join(value.encode("ascii")
                                for value in command) + b"\0").hexdigest(),
                "launcher_session": 7001, "launcher_start_tick": 0,
                "pid": 7002, "start_tick": 456})
        backend._owned_session_members.assert_not_called()

    def test_live_candidate_start_defensively_rejects_forced_archive(self):
        backend = object.__new__(transition.LiveBackend)
        backend.plan = transition.TransitionPlan()
        backend.stop_record = {"forced": False}
        backend.archive_record = {"forced_stop": True}
        backend._i2c_readers = lambda: ()
        with self.assertRaisesRegex(
                transition.TransitionError, "lacks an empty-I2C archive seal"):
            transition.LiveBackend.start_candidate(backend)

    def test_csv_reader_preflight_rejects_any_extra_pid(self):
        backend = object.__new__(transition.LiveBackend)
        backend.plan = transition.TransitionPlan()
        backend._i2c_readers = lambda: (backend.plan.old_pid,)
        backend._fuser = lambda _path: (backend.plan.old_pid, 919688)
        with self.assertRaisesRegex(transition.TransitionError,
                                    "auditor or unknown"):
            transition.LiveBackend._require_exclusive_old_readers(backend)

    def test_csv_reader_preflight_accepts_only_old_writer(self):
        backend = object.__new__(transition.LiveBackend)
        backend.plan = transition.TransitionPlan()
        backend._i2c_readers = lambda: (backend.plan.old_pid,)
        backend._fuser = lambda _path: (backend.plan.old_pid,)
        self.assertEqual(
            transition.LiveBackend._require_exclusive_old_readers(backend),
            {"csv_readers": [backend.plan.old_pid],
             "i2c_readers": [backend.plan.old_pid]})

    def test_final_old_identity_replays_owned_session_and_process_group(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            os.chmod(root, 0o700)
            source = root / "old.py"
            csv_path = root / "thermal.csv"
            pid_path = root / "sampler.pid"
            source.write_bytes(b"old source\n")
            csv_path.write_bytes(b"sample\n")
            pid_path.write_bytes(b"7001\n")
            for path in (source, csv_path, pid_path):
                os.chmod(path, 0o444)
            plan = make_test_code_seal(root, legacy=b"old source\n")
            plan = replace(
                plan, controller_uid=os.geteuid(), controller_gid=os.getegid(),
                old_source=source, old_csv=csv_path, old_pid_file=pid_path,
                old_source_sha256=hashlib.sha256(b"old source\n").hexdigest())
            identity = {
                "argv": list(plan.replacement_old_argv),
                "cmdline_sha256": plan.replacement_old_cmdline_sha256,
                "pid": 7001,
                "process_group": 9000, "session": 9000,
                "start_tick": 123,
            }
            backend = object.__new__(transition.LiveBackend)
            backend.plan = plan
            backend.p32 = mock.Mock()
            backend.old_preflight = {
                "source": transition.file_binding(source, with_hash=True)}
            identity_receipt = {
                "boot_id": "boot",
                "targets": {
                    "child": {"identity": identity, "state": "present"},
                    "launcher": {"state": "absent"},
                },
            }
            backend._inspect_replacement_identities = mock.Mock(
                return_value=identity_receipt)
            backend.p32._parse_thermal_csv.return_value = (
                {"monotonic_s": "%.6f" % transition.time.monotonic()},)
            backend._i2c_readers = lambda: (7001,)
            backend._fuser = lambda _path: (7001,)
            csv_binding = transition.file_binding(csv_path, with_hash=False)
            source_binding = transition.file_binding(
                plan.legacy_sampler, with_hash=True)
            bootstrap = {
                "fd_receipt": {"csv": {
                    "binding": {**csv_binding, "size": 1},
                    "fd": {"fd": 9, "device": csv_binding["device"],
                           "inode": csv_binding["inode"],
                           "access_mode": os.O_WRONLY,
                           "link_path": str(csv_path), "size": 1}}}}
            replay = {
                "fd_receipt": {"csv": {
                    "binding": dict(csv_binding),
                    "fd": {"fd": 9, "device": csv_binding["device"],
                           "inode": csv_binding["inode"],
                           "access_mode": os.O_WRONLY,
                           "link_path": str(csv_path),
                           "size": csv_binding["size"]}}}}
            backend.legacy_fd_bootstrap = bootstrap
            backend.legacy_fd_replay = None
            backend._inspect_sampler_fds = lambda *_args, **_kwargs: replay
            restored = {
                "boot_id": "boot",
                "cmdline_sha256": plan.replacement_old_cmdline_sha256,
                "csv_initial_size": csv_binding["size"],
                "csv_live_identity": {
                    key: csv_binding[key] for key in (
                        "device", "gid", "inode", "mode", "nlink", "uid")},
                "launcher_session": 9000,
                "pid": 7001,
                "pid_binding": transition.file_binding(pid_path, with_hash=True),
                "source_binding": source_binding,
                "source_path": str(plan.legacy_sampler),
                "fd_provenance_bootstrap": bootstrap,
                "start_tick": 123,
            }
            backend._replacement_launcher_roster = lambda _restored: [
                {"pid": 7001}]
            result = transition.LiveBackend._revalidate_restored_old(
                backend, restored)
            self.assertEqual(result["identity"], identity)
            backend.p32._parse_thermal_csv.return_value = ()
            with self.assertRaisesRegex(
                    transition.TransitionError, "CSV content changed"):
                transition.LiveBackend._revalidate_restored_old(
                    backend, restored)
            backend.p32._parse_thermal_csv.return_value = (
                {"monotonic_s": "%.6f" % transition.time.monotonic()},)
            bad_identity = dict(identity)
            bad_identity["process_group"] = 9001
            backend._inspect_replacement_identities.return_value = {
                "boot_id": "boot",
                "targets": {
                    "child": {"identity": bad_identity, "state": "present"},
                    "launcher": {"state": "absent"},
                },
            }
            with self.assertRaisesRegex(
                    transition.TransitionError, "changed before audit"):
                transition.LiveBackend._revalidate_restored_old(
                    backend, restored)

    def test_future_audit_binding_preserves_archive_sha_and_new_inode(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            os.chmod(root, 0o700)
            archive = root / "archive.csv"
            archive.write_bytes(b"pre-dry\n")
            os.chmod(archive, 0o444)
            archive_binding = transition.file_binding(archive, with_hash=True)
            archive_binding["uid"] = 0
            plan = replace(
                transition.TransitionPlan(), controller_uid=os.geteuid(),
                old_archive=archive, audit_binding=root / "future.json")
            backend = object.__new__(transition.LiveBackend)
            backend.plan = plan
            backend.tools = transition.capture_tool_records()
            backend._replay_external_state = lambda *_args: {
                "prepublication": "replayed"}
            backend._revalidate_restored_old = lambda _restored: {
                "csv_current_size": 9, "identity": {"pid": 8123}}
            restored = {
                "pid": 8123,
                "csv_live_identity": {"device": 7, "inode": 99},
                "archived_pre_dry_sha256": archive_binding["sha256"],
                "archived_pre_dry_path": str(archive),
            }
            real_file_binding = transition.file_binding

            def root_archive_binding(path, *, with_hash):
                binding = real_file_binding(path, with_hash=with_hash)
                if path == archive:
                    binding["uid"] = 0
                return binding

            with mock.patch.object(
                    transition, "file_binding", side_effect=root_archive_binding), \
                    mock.patch.object(
                        transition, "verify_running_interpreter",
                        return_value={"device": 1, "inode": 2,
                                      "argv_path": "/sealed/python",
                                      "resolved_path": "/sealed/python"}):
                result = transition.LiveBackend.publish_audit_binding(
                    backend,
                    {"path": str(archive), "binding": archive_binding},
                    restored, False, None,
                    {"count": 1, "head_sha256": "c" * 64,
                     "roster": [{"sequence": 0}]})
            value = transition.load_canonical(
                plan.audit_binding, transition.AUDIT_BINDING_SCHEMA,
                "test audit binding")
            self.assertEqual(value, result["value"])
            self.assertEqual(value["tools"], backend.tools)
            self.assertEqual(value["receipt_chain_prefix"]["head_sha256"],
                             "c" * 64)
            self.assertEqual(result["tools_after_audit"], backend.tools)
            self.assertEqual(
                value["archived_pre_dry"]["binding"]["sha256"],
                archive_binding["sha256"])
            self.assertEqual(value["live_old_sampler"]
                             ["csv_live_identity"]["inode"],
                             99)

    def test_audit_rejects_success_without_retained_candidate_acceptance(self):
        backend = object.__new__(transition.LiveBackend)
        backend._replay_external_state = mock.Mock()
        with self.assertRaisesRegex(
                transition.TransitionError, "lacks retained candidate"):
            transition.LiveBackend.publish_audit_binding(
                backend, None, {}, True, None,
                {"count": 1, "head_sha256": "c" * 64,
                 "roster": [{"sequence": 0}]})
        backend._replay_external_state.assert_not_called()

    def test_future_audit_binding_rehashes_the_archive(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            os.chmod(root, 0o700)
            archive = root / "archive.csv"
            archive.write_bytes(b"pre-dry\n")
            os.chmod(archive, 0o444)
            stale = transition.file_binding(archive, with_hash=True)
            os.chmod(archive, 0o644)
            archive.write_bytes(b"changed!\n")
            os.chmod(archive, 0o444)
            plan = replace(
                transition.TransitionPlan(), controller_uid=os.geteuid(),
                old_archive=archive, audit_binding=root / "future.json")
            backend = object.__new__(transition.LiveBackend)
            backend.plan = plan
            backend._replay_external_state = lambda *_args: {
                "prepublication": "replayed"}
            backend._revalidate_restored_old = lambda _restored: {
                "csv_current_size": 9, "identity": {"pid": 8123}}
            restored = {
                "pid": 8123,
                "archived_pre_dry_sha256": stale["sha256"],
                "archived_pre_dry_path": str(archive),
            }
            with self.assertRaisesRegex(
                    transition.TransitionError, "archive changed"):
                transition.LiveBackend.publish_audit_binding(
                    backend, {"path": str(archive), "binding": stale},
                    restored, False, None,
                    {"count": 1, "head_sha256": "c" * 64,
                     "roster": [{"sequence": 0}]})
            self.assertFalse(plan.audit_binding.exists())

    def test_future_audit_binding_fails_if_post_write_tool_replay_changes(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            os.chmod(root, 0o700)
            archive = root / "archive.csv"
            archive.write_bytes(b"pre-dry\n")
            os.chmod(archive, 0o444)
            archive_binding = transition.file_binding(archive, with_hash=True)
            archive_binding["uid"] = 0
            plan = replace(
                transition.TransitionPlan(), controller_uid=os.geteuid(),
                old_archive=archive, audit_binding=root / "future.json")
            tools = transition.capture_tool_records()
            changed_tools = json.loads(json.dumps(tools))
            changed_tools["sudo"]["binding"]["inode"] += 1
            backend = object.__new__(transition.LiveBackend)
            backend.plan = plan
            backend.tools = tools
            backend._verify_tools = lambda: tools
            backend._replay_external_state = lambda *_args: {
                "prepublication": "replayed"}
            backend._revalidate_restored_old = lambda _restored: {
                "csv_current_size": 9, "identity": {"pid": 8123}}
            restored = {
                "pid": 8123,
                "archived_pre_dry_sha256": archive_binding["sha256"],
                "archived_pre_dry_path": str(archive),
            }
            real_file_binding = transition.file_binding

            def root_archive_binding(path, *, with_hash):
                binding = real_file_binding(path, with_hash=with_hash)
                if path == archive:
                    binding["uid"] = 0
                return binding

            with mock.patch.object(
                    transition, "file_binding",
                    side_effect=root_archive_binding), \
                    mock.patch.object(
                        transition, "verify_running_interpreter",
                        return_value={"device": 1, "inode": 2,
                                      "argv_path": "/sealed/python",
                                      "resolved_path": "/sealed/python"}), \
                    mock.patch.object(
                        transition, "verify_tool_records",
                        return_value=changed_tools):
                with self.assertRaisesRegex(
                        transition.TransitionError,
                        "tool identities changed after audit"):
                    transition.LiveBackend.publish_audit_binding(
                        backend,
                        {"path": str(archive), "binding": archive_binding},
                        restored, False, None,
                        {"count": 1, "head_sha256": "c" * 64,
                         "roster": [{"sequence": 0}]})
            self.assertTrue(plan.audit_binding.exists())

    def test_final_replay_rejects_postpublication_audit_replacement(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            os.chmod(root, 0o700)
            audit_path = root / "future.json"
            value = transition.sealed_record(
                transition.AUDIT_BINDING_SCHEMA, {
                    "candidate_accept": None,
                    "prepublication_replay": {"state": "clean"},
                })
            binding = transition.write_new(
                audit_path, transition.canonical_json(value),
                owner_uid=os.geteuid())
            backend = object.__new__(transition.LiveBackend)
            backend.plan = replace(
                transition.TransitionPlan(), controller_uid=os.geteuid(),
                audit_binding=audit_path)
            backend._replay_external_state = lambda *_args: {"state": "clean"}
            audit = {"binding": binding, "path": str(audit_path),
                     "value": value}
            replay = transition.LiveBackend.final_replay(
                backend, None, None, {}, audit)
            self.assertEqual(replay["audit"]["binding"], binding)
            replacement = transition.sealed_record(
                transition.AUDIT_BINDING_SCHEMA, {
                    "candidate_accept": None,
                    "prepublication_replay": {"state": "changed"},
                })
            os.chmod(audit_path, 0o644)
            audit_path.write_bytes(transition.canonical_json(replacement))
            os.chmod(audit_path, 0o444)
            with self.assertRaisesRegex(
                    transition.TransitionError, "changed after publication"):
                transition.LiveBackend.final_replay(
                    backend, None, None, {}, audit)

    def test_recovery_discovers_archive_after_post_rename_exception(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            os.chmod(root, 0o700)
            archive = root / "archive.csv"
            archive.write_bytes(b"pre-dry\n")
            os.chmod(archive, 0o444)
            plan = replace(
                transition.TransitionPlan(), controller_uid=os.geteuid(),
                old_csv=root / "thermal.csv", old_pid_file=root / "sampler.pid",
                old_archive=archive,
                old_unclean_archive=root / "unclean.csv",
                old_stale_pid_archive=root / "stale.pid")
            backend = object.__new__(transition.LiveBackend)
            backend.plan = plan
            backend.archive_record = None
            backend._fuser = lambda _path: ()
            record = transition.LiveBackend._ensure_recovery_archive(
                backend, None)
            self.assertEqual(record["path"], str(archive))
            self.assertEqual(record["binding"]["sha256"],
                             hashlib.sha256(b"pre-dry\n").hexdigest())

    def test_recovery_rejects_archive_with_wrong_preflight_inode(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            os.chmod(root, 0o700)
            archive = root / "archive.csv"
            archive.write_bytes(b"pre-dry\n")
            os.chmod(archive, 0o444)
            binding = transition.file_binding(archive, with_hash=True)
            expected = dict(binding)
            expected["inode"] += 1
            plan = replace(
                transition.TransitionPlan(), controller_uid=os.geteuid(),
                old_csv=root / "thermal.csv", old_pid_file=root / "sampler.pid",
                old_archive=archive,
                old_unclean_archive=root / "unclean.csv",
                old_stale_pid_archive=root / "stale.pid")
            backend = object.__new__(transition.LiveBackend)
            backend.plan = plan
            backend.old_preflight = {"csv": expected}
            backend.archive_record = None
            with self.assertRaisesRegex(
                    transition.TransitionError, "preflight CSV inode"):
                transition.LiveBackend._ensure_recovery_archive(backend, None)

    def test_recovery_rejects_alternate_or_unbound_stale_archive_names(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            os.chmod(root, 0o700)
            clean = root / "clean.csv"
            unclean = root / "unclean.csv"
            stale = root / "stale.pid"
            for path, raw in ((clean, b"pre-dry\n"),
                              (unclean, b"other\n")):
                path.write_bytes(raw)
                os.chmod(path, 0o444)
            plan = replace(
                transition.TransitionPlan(), controller_uid=os.geteuid(),
                old_csv=root / "thermal.csv",
                old_pid_file=root / "sampler.pid", old_archive=clean,
                old_unclean_archive=unclean,
                old_stale_pid_archive=stale)
            backend = object.__new__(transition.LiveBackend)
            backend.plan = plan
            backend.archive_record = None
            with self.assertRaisesRegex(
                    transition.TransitionError, "multiple pre-dry archives"):
                transition.LiveBackend._ensure_recovery_archive(backend, None)
            unclean.unlink()
            stale.write_bytes(b"123\n")
            os.chmod(stale, 0o444)
            with self.assertRaisesRegex(
                    transition.TransitionError, "stale PID evidence"):
                transition.LiveBackend._ensure_recovery_archive(backend, None)

    def test_recovery_archives_canonical_csv_and_stale_pid(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            os.chmod(root, 0o700)
            csv_path = root / "thermal.csv"
            pid_path = root / "sampler.pid"
            csv_path.write_bytes(b"pre-dry\n")
            pid_path.write_bytes(b"123\n")
            os.chmod(csv_path, 0o444)
            os.chmod(pid_path, 0o444)
            plan = replace(
                transition.TransitionPlan(), controller_uid=os.geteuid(),
                old_csv=csv_path, old_pid_file=pid_path,
                old_archive=root / "archive.csv",
                old_unclean_archive=root / "unclean.csv",
                old_stale_pid_archive=root / "stale.pid")
            backend = object.__new__(transition.LiveBackend)
            backend.plan = plan
            backend.archive_record = None
            backend._fuser = lambda _path: ()
            record = transition.LiveBackend._ensure_recovery_archive(
                backend, None)
            self.assertEqual(record["path"], str(plan.old_unclean_archive))
            self.assertEqual(record["stale_pid"]["path"],
                             str(plan.old_stale_pid_archive))
            self.assertFalse(csv_path.exists())
            self.assertFalse(pid_path.exists())


if __name__ == "__main__":
    unittest.main()
