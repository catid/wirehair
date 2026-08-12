#!/usr/bin/env python3
"""Run and seal the bounded WH2 n16 decoder phase-attribution screen.

The phase screen is diagnostic evidence, not canonical timing evidence.  It
consumes one already-complete Two07 recovery campaign and binds that campaign
to an exact same-source/same-binary 24-coordinate n16 execution.  Completion
is a publish-wins transaction: every named dependency remains pinned while an
unnamed ``O_TMPFILE`` summary is validated, the hard wall is torn down, and
the summary is linked exactly once.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import os
from pathlib import Path
import selectors
import signal
import stat
import subprocess
import sys
import time
from typing import (Any, Callable, Dict, List, Mapping, Optional, Sequence,
                    Set, Tuple)

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_benchmark_contract as contract_api
import wh2_native_short_screen as native_api
import wh2_phase_attribution as phase_api
import wh2_run_native_recovery_screen as recovery_api
import wh2_run_native_short_screen as runner_api


PHASE_FREEZE_SCHEMA = "wirehair.wh2.native-phase-attribution-freeze.v1"
PHASE_EXECUTION_SCHEMA = \
    "wirehair.wh2.native-phase-attribution-execution.v1"
PHASE_SUMMARY_SCHEMA = "wirehair.wh2.native-phase-attribution-run.v1"
PHASE_RECOVERY_BINDING_SCHEMA = \
    "wirehair.wh2.phase-attribution-recovery-binding.v1"

PHASE_WORKER_COUNT = 8
PHASE_RECORD_COUNT = 24
PHASE_WAVE_COUNT = 3
PHASE_PROFILE_ORDINAL = 0
MAX_PHASE_RECORD_BYTES = phase_api.MAX_NATIVE_LINE_BYTES
MAX_PHASE_METADATA_BYTES = 1024 * 1024
MAX_PHASE_TRACE_STDERR_BYTES = 8192
PHASE_NATIVE_STREAM_BYTE_CAP = PHASE_RECORD_COUNT * MAX_PHASE_RECORD_BYTES

PHASE_SOURCE_NAMES = {
    "summary": "run-summary.json",
    "freeze": "phase-freeze.json",
    "trace": "phase-traces.jsonl",
    "native": "phase-native-results.jsonl",
    "analysis": "phase-analysis.json",
    "execution": "phase-execution.json",
    "sampler": "sampler-attestation.json",
}
PHASE_DEPENDENCY_NAMES = {
    key: name for key, name in PHASE_SOURCE_NAMES.items()
    if key != "summary"
}
PHASE_ARTIFACT_BYTE_LIMITS = {
    "run-summary.json": MAX_PHASE_METADATA_BYTES,
    "phase-freeze.json": MAX_PHASE_METADATA_BYTES,
    "phase-traces.jsonl": phase_api.MAX_TRACE_BYTES,
    "phase-native-results.jsonl": PHASE_NATIVE_STREAM_BYTE_CAP,
    "phase-analysis.json": MAX_PHASE_METADATA_BYTES,
    "phase-execution.json": MAX_PHASE_METADATA_BYTES,
    "sampler-attestation.json": MAX_PHASE_METADATA_BYTES,
}
MAX_PHASE_BUNDLE_BYTES = sum(PHASE_ARTIFACT_BYTE_LIMITS.values())

RECOVERY_BINDING_FIELDS = frozenset((
    "schema", "campaign_directory", "campaign_directory_device",
    "campaign_directory_inode", "candidate_id", "candidate_arm",
    "candidate_arm_descriptor_sha256", "source_git_commit",
    "worker_binary_sha256", "timing_proxy_witness_sha256",
    "recovery_run_summary_sha256", "recovery_freeze_sha256",
    "recovery_result_sha256", "recovery_execution_receipt_sha256",
    "binding_sha256",
))
PHASE_FREEZE_FIELDS = frozenset((
    "schema", "evidence_kind", "phase", "profile", "profile_ordinal",
    "contract_sha256", "source_git_commit", "worker_binary_sha256",
    "production_head_descriptor_sha256", "candidate_arm",
    "candidate_arm_descriptor_sha256", "timing_proxy_witness_sha256",
    "recovery_binding", "recovery_binding_sha256",
    "trace_manifest_sha256", "trace_file_sha256", "worker_cpus",
    "coordinate_cpus", "controller_cpu", "host_identity", "commands",
))
PHASE_COMMAND_FIELDS = frozenset((
    "description", "trace", "workers", "jobs", "job_ordinals",
))
PHASE_EXECUTION_FIELDS = frozenset((
    "schema", "status", "source_git_commit", "worker_binary_sha256",
    "contract_sha256", "recovery_binding_sha256", "freeze_sha256",
    "trace_manifest_sha256", "trace_file_sha256", "native_stream_sha256",
    "analysis_sha256", "sampler_attestation_sha256", "record_count",
    "worker_cpus", "coordinate_cpus",
    "workers", "controller_cpu", "worker_start_monotonic_ns",
    "worker_end_monotonic_ns", "thermal", "receipt_sha256",
))
PHASE_SUMMARY_FIELDS = frozenset((
    "schema", "status", "output_dir", "source_git_commit",
    "worker_binary_sha256", "contract_sha256", "candidate_arm",
    "candidate_arm_descriptor_sha256", "timing_proxy_witness_sha256",
    "recovery_binding_sha256", "recovery_run_summary_sha256",
    "recovery_freeze_sha256", "recovery_execution_receipt_sha256",
    "phase_freeze_sha256", "trace_manifest_sha256", "trace_file_sha256",
    "native_stream_sha256", "analysis_sha256",
    "sampler_attestation_sha256", "execution_receipt_sha256",
    "phase_records", "controller_cpu", "worker_cpus", "thermal_samples",
    "cpu_tctl_max_millic", "dimm_max_millic", "decision_outcome",
    "summary_sha256",
))
PHASE_WORKER_FIELDS = frozenset((
    "cpu", "pid", "process_start_ticks",
))
HOST_IDENTITY_FIELDS = frozenset((
    "controller_cpu", "hostname", "kernel_release", "machine", "system",
))


class PhaseNativeRunnerError(RuntimeError):
    """The sealed native phase-attribution transaction is invalid."""


def fail(message: str) -> None:
    raise PhaseNativeRunnerError(message)


def _is_sha256(value: Any) -> bool:
    return isinstance(value, str) and \
        contract_api.SHA256.fullmatch(value) is not None


def _hash_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _self_hash(value: Mapping[str, Any], field: str) -> str:
    return contract_api.sha256_json({
        key: item for key, item in value.items() if key != field
    })


def _exact_mapping(value: Any, fields: Sequence[str], context: str) \
        -> Mapping[str, Any]:
    if not isinstance(value, dict) or set(value) != set(fields):
        fail("{} has an unexpected exact schema".format(context))
    return value


def _parse_object_bytes(data: bytes, context: str) -> Mapping[str, Any]:
    try:
        return runner_api._parse_canonical_line(data, context)
    except runner_api.RunnerError as exc:
        fail(str(exc))
    raise AssertionError("unreachable")


def _canonical_object_bytes(value: Mapping[str, Any]) -> bytes:
    return (contract_api.canonical_json(value) + "\n").encode("utf-8")


def _public_phase_description(
        description: Mapping[str, Any]) -> Mapping[str, Any]:
    if (not isinstance(description, Mapping) or
            not phase_api.DESCRIPTION_FIELDS.issubset(description)):
        fail("phase worker description lacks its public identity")
    return {
        field: description[field] for field in phase_api.DESCRIPTION_FIELDS
    }


def _campaign_binding(campaign: Mapping[str, Any]) -> Mapping[str, Any]:
    """Extract the exact immutable recovery identity needed by this screen."""
    summary = campaign.get("summary")
    freeze = campaign.get("freeze")
    receipt = campaign.get("receipt")
    witness = campaign.get("timing_proxy_witness")
    identity = campaign.get("directory_identity")
    if (not isinstance(summary, Mapping) or
            not isinstance(freeze, Mapping) or
            not isinstance(receipt, Mapping) or
            not isinstance(witness, Mapping) or
            not isinstance(identity, tuple) or len(identity) != 2 or
            any(type(item) is not int or item < 0 for item in identity)):
        fail("completed recovery campaign lacks its validated exact identity")
    arms = freeze.get("arms")
    if not isinstance(arms, list) or len(arms) != 4:
        fail("completed recovery campaign lacks its exact four-arm freeze")
    candidate = arms[3]
    if (not isinstance(candidate, Mapping) or
            candidate.get("arm") != phase_api.LEFT_ARM or
            candidate.get("arm_descriptor_sha256") !=
                phase_api.EXPECTED_ARMS[2][2]):
        fail("completed recovery campaign does not bind the Two07 candidate")
    values = {
        "schema": PHASE_RECOVERY_BINDING_SCHEMA,
        "campaign_directory": campaign.get("directory"),
        "campaign_directory_device": identity[0],
        "campaign_directory_inode": identity[1],
        "candidate_id": campaign.get("candidate_id"),
        "candidate_arm": campaign.get("candidate_arm"),
        "candidate_arm_descriptor_sha256":
            candidate.get("arm_descriptor_sha256"),
        "source_git_commit": summary.get("source_git_commit"),
        "worker_binary_sha256": summary.get("worker_binary_sha256"),
        "timing_proxy_witness_sha256":
            summary.get("timing_proxy_witness_sha256"),
        "recovery_run_summary_sha256": summary.get("summary_sha256"),
        "recovery_freeze_sha256":
            summary.get("recovery_freeze_sha256"),
        "recovery_result_sha256":
            summary.get("recovery_result_sha256"),
        "recovery_execution_receipt_sha256":
            summary.get("recovery_execution_receipt_sha256"),
    }
    if (not isinstance(values["campaign_directory"], str) or
            not values["campaign_directory"] or
            values["candidate_id"] != "two07" or
            values["candidate_arm"] != phase_api.LEFT_ARM or
            values["source_git_commit"] != freeze.get("source_git_commit") or
            values["worker_binary_sha256"] !=
                freeze.get("arms", [{}])[0].get("binary_sha256") or
            values["timing_proxy_witness_sha256"] !=
                freeze.get("timing_proxy_witness_sha256") or
            values["timing_proxy_witness_sha256"] !=
                contract_api.sha256_json(witness) or
            values["recovery_freeze_sha256"] !=
                receipt.get("freeze_manifest_sha256") or
            values["recovery_result_sha256"] !=
                receipt.get("result_stream_sha256") or
            values["recovery_execution_receipt_sha256"] !=
                receipt.get("receipt_sha256") or
            any(not _is_sha256(values[field]) for field in (
                "candidate_arm_descriptor_sha256", "worker_binary_sha256",
                "timing_proxy_witness_sha256",
                "recovery_run_summary_sha256", "recovery_freeze_sha256",
                "recovery_result_sha256",
                "recovery_execution_receipt_sha256"))):
        fail("completed recovery campaign binding is internally inconsistent")
    values["binding_sha256"] = contract_api.sha256_json(values)
    return values


def _load_recovery_binding(
        contract: Mapping[str, Any], recovery_dir: Path,
        ) -> Tuple[Mapping[str, Any], Mapping[str, Any]]:
    try:
        campaign = recovery_api.load_completed_campaign(
            contract, recovery_dir)
    except (recovery_api.RecoveryRunnerError,
            runner_api.RunnerError) as exc:
        fail(str(exc))
    return campaign, _campaign_binding(campaign)


def _validate_recovery_binding(
        actual: Any, expected: Mapping[str, Any]) -> Mapping[str, Any]:
    binding = _exact_mapping(
        actual, RECOVERY_BINDING_FIELDS, "phase recovery binding")
    if (binding.get("schema") != PHASE_RECOVERY_BINDING_SCHEMA or
            binding.get("binding_sha256") !=
                _self_hash(binding, "binding_sha256") or
            contract_api.canonical_json(binding) !=
                contract_api.canonical_json(expected)):
        fail("phase evidence binds another completed recovery campaign")
    return binding


class PhaseJob:
    """One exact diagnostic P command; the worker-domain ordinal is even."""

    def __init__(self, coordinate: int) -> None:
        if type(coordinate) is not int or not 0 <= coordinate < \
                PHASE_RECORD_COUNT:
            fail("phase job coordinate is outside the exact n16 domain")
        self.kind = "phase_attribution"
        self.cell = coordinate
        self.item = PHASE_PROFILE_ORDINAL
        self.ordinal = coordinate * 2 + PHASE_PROFILE_ORDINAL

    def command(self) -> bytes:
        return "P {} {}\n".format(self.cell, self.item).encode("ascii")


def _phase_job_waves() -> List[List[PhaseJob]]:
    jobs = [PhaseJob(coordinate) for coordinate in range(PHASE_RECORD_COUNT)]
    waves = [
        jobs[index:index + PHASE_WORKER_COUNT]
        for index in range(0, PHASE_RECORD_COUNT, PHASE_WORKER_COUNT)
    ]
    if (len(waves) != PHASE_WAVE_COUNT or
            any(len(wave) != PHASE_WORKER_COUNT for wave in waves) or
            [job.ordinal for wave in waves for job in wave] !=
                list(range(0, 48, 2)) or
            [job.command() for wave in waves for job in wave] != [
                "P {} 0\n".format(value).encode("ascii")
                for value in range(PHASE_RECORD_COUNT)]):
        fail("phase job waves differ from the exact 24-coordinate n16 domain")
    return waves


def _phase_commands(
        description: Mapping[str, Any], cpus: Sequence[int],
        ) -> Mapping[str, Any]:
    worker = description.get("resolved_path")
    if not isinstance(worker, str) or not worker:
        fail("phase worker description lacks its resolved executable path")
    waves = _phase_job_waves()
    return {
        "description": [worker, "--describe"],
        "trace": [worker, "--emit-traces", "phase-attribution"],
        "workers": [[worker, "--worker", str(cpu)] for cpu in cpus],
        "jobs": [
            ["P", str(job.cell), str(job.item)]
            for wave in waves for job in wave
        ],
        "job_ordinals": [job.ordinal for wave in waves for job in wave],
    }


def _phase_freeze(
        contract: Mapping[str, Any], description: Mapping[str, Any],
        recovery_binding: Mapping[str, Any], cpus: Sequence[int],
        controller_cpu: int, trace_manifest_sha256: str,
        trace_file_sha256: str) -> Mapping[str, Any]:
    worker_cpus = list(cpus)
    coordinate_cpus = worker_cpus * PHASE_WAVE_COUNT
    if (len(worker_cpus) != PHASE_WORKER_COUNT or
            worker_cpus != sorted(set(worker_cpus)) or
            any(type(cpu) is not int or cpu < 0 for cpu in worker_cpus) or
            type(controller_cpu) is not int or controller_cpu < 0 or
            controller_cpu in worker_cpus):
        fail("phase freeze requires eight sorted unique workers and a controller")
    try:
        phase_api._validate_description(
            _public_phase_description(description),
            recovery_binding["worker_binary_sha256"])
    except phase_api.PhaseRunnerError as exc:
        fail(str(exc))
    if (description.get("source_git_commit") !=
            recovery_binding["source_git_commit"] or
            description["arms"][2]["arm_descriptor_sha256"] !=
                recovery_binding["candidate_arm_descriptor_sha256"] or
            not _is_sha256(trace_manifest_sha256) or
            not _is_sha256(trace_file_sha256)):
        fail("phase freeze differs from its recovery/source/trace identity")
    return {
        "schema": PHASE_FREEZE_SCHEMA,
        "evidence_kind": "phase_attribution",
        "phase": "development",
        "profile": phase_api.PROFILE_NAME,
        "profile_ordinal": PHASE_PROFILE_ORDINAL,
        "contract_sha256": contract_api.contract_sha256(contract),
        "source_git_commit": description["source_git_commit"],
        "worker_binary_sha256": description["binary_sha256"],
        "production_head_descriptor_sha256":
            description["arms"][0]["arm_descriptor_sha256"],
        "candidate_arm": phase_api.LEFT_ARM,
        "candidate_arm_descriptor_sha256":
            description["arms"][2]["arm_descriptor_sha256"],
        "timing_proxy_witness_sha256":
            recovery_binding["timing_proxy_witness_sha256"],
        "recovery_binding": dict(recovery_binding),
        "recovery_binding_sha256": recovery_binding["binding_sha256"],
        "trace_manifest_sha256": trace_manifest_sha256,
        "trace_file_sha256": trace_file_sha256,
        "worker_cpus": worker_cpus,
        "coordinate_cpus": coordinate_cpus,
        "controller_cpu": controller_cpu,
        "host_identity": runner_api._host_identity(controller_cpu),
        "commands": _phase_commands(description, worker_cpus),
    }


def _validate_physical_layout(
        worker_cpus: Sequence[int], controller_cpu: int) -> None:
    try:
        worker_cores = [
            runner_api._cpu_topology(cpu, runner_api.CPU_SYSFS_ROOT)
            for cpu in worker_cpus
        ]
        controller_core = runner_api._cpu_topology(
            controller_cpu, runner_api.CPU_SYSFS_ROOT)
    except runner_api.RunnerError as exc:
        fail(str(exc))
    if (len(worker_cores) != PHASE_WORKER_COUNT or
            len(set(worker_cores)) != PHASE_WORKER_COUNT or
            controller_core in set(worker_cores)):
        fail("phase worker/controller roster is not nine isolated physical cores")


def _validate_phase_freeze(
        contract: Mapping[str, Any], freeze: Any,
        recovery_binding: Mapping[str, Any],
        worker_path: Optional[str] = None) -> Mapping[str, Any]:
    value = _exact_mapping(freeze, PHASE_FREEZE_FIELDS, "phase freeze")
    worker_cpus = value.get("worker_cpus")
    coordinate_cpus = value.get("coordinate_cpus")
    if (value.get("schema") != PHASE_FREEZE_SCHEMA or
            value.get("evidence_kind") != "phase_attribution" or
            value.get("phase") != "development" or
            value.get("profile") != phase_api.PROFILE_NAME or
            type(value.get("profile_ordinal")) is not int or
            value.get("profile_ordinal") != PHASE_PROFILE_ORDINAL or
            value.get("contract_sha256") !=
                contract_api.contract_sha256(contract) or
            value.get("source_git_commit") !=
                recovery_binding["source_git_commit"] or
            value.get("worker_binary_sha256") !=
                recovery_binding["worker_binary_sha256"] or
            value.get("production_head_descriptor_sha256") !=
                phase_api.EXPECTED_ARMS[0][2] or
            value.get("candidate_arm") != phase_api.LEFT_ARM or
            value.get("candidate_arm_descriptor_sha256") !=
                recovery_binding["candidate_arm_descriptor_sha256"] or
            value.get("timing_proxy_witness_sha256") !=
                recovery_binding["timing_proxy_witness_sha256"] or
            value.get("recovery_binding_sha256") !=
                recovery_binding["binding_sha256"] or
            not _is_sha256(value.get("trace_manifest_sha256")) or
            not _is_sha256(value.get("trace_file_sha256")) or
            not isinstance(worker_cpus, list) or
            len(worker_cpus) != PHASE_WORKER_COUNT or
            worker_cpus != sorted(set(worker_cpus)) or
            any(type(cpu) is not int or cpu < 0 for cpu in worker_cpus) or
            coordinate_cpus != worker_cpus * PHASE_WAVE_COUNT or
            type(value.get("controller_cpu")) is not int or
            value["controller_cpu"] < 0 or
            value["controller_cpu"] in worker_cpus or
            not isinstance(value.get("host_identity"), dict) or
            set(value["host_identity"]) != HOST_IDENTITY_FIELDS or
            value["host_identity"].get("controller_cpu") !=
                value["controller_cpu"]):
        fail("phase freeze differs from its exact sealed geometry")
    _validate_physical_layout(worker_cpus, value["controller_cpu"])
    _validate_recovery_binding(value.get("recovery_binding"), recovery_binding)
    commands = _exact_mapping(
        value.get("commands"), PHASE_COMMAND_FIELDS, "phase command ledger")
    if worker_path is None:
        try:
            worker_path = commands["description"][0]
        except (IndexError, KeyError, TypeError):
            fail("phase command ledger lacks its worker path")
    expected_commands = _phase_commands(
        {"resolved_path": worker_path}, worker_cpus)
    if commands != expected_commands:
        fail("phase command ledger differs from its exact P schedule")
    return value


def _finish_capped_phase_process(
        process: Optional[subprocess.Popen],
        selector: Optional[selectors.BaseSelector], successful: bool,
        reaped: bool, primary: Optional[BaseException]) -> None:
    """Exhaust group shutdown/reap/pipe cleanup and preserve control flow."""
    failures: List[Tuple[str, BaseException]] = []
    terminal = []
    leader_reaped = reaped
    if process is not None:
        try:
            # Popen.wait() stores returncode before returning.  Reading it
            # closes the async-interrupt sliver between that store and the
            # caller's subsequent `reaped = True` assignment without invoking
            # poll()/wait() and therefore without reaping here.
            leader_reaped = leader_reaped or process.returncode is not None
        except BaseException as cleanup:
            failures.append(("phase trace reap-state inspection failed",
                             cleanup))
            # Destructive signaling is never safe when ownership is uncertain.
            leader_reaped = True
    if process is not None and not successful and not leader_reaped:
        # start_new_session=True makes the child PID its process-group ID.
        # Keep the leader unreaped until both signals have been attempted.  A
        # dead leader remains a zombie, preventing its numeric PID/PGID from
        # being recycled between TERM and KILL.  Never signal after wait()
        # returned: at that point the same number could name an unrelated
        # same-UID process group.
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
        except BaseException as cleanup:
            failures.append(("phase trace process-group terminate failed",
                             cleanup))
        # Bounded grace without wait() or poll(): either could reap the leader
        # and open a destructive PGID-reuse window before SIGKILL.
        try:
            time.sleep(0.05)
        except BaseException as cleanup:
            failures.append(("phase trace termination grace failed", cleanup))
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        except BaseException as cleanup:
            failures.append(("phase trace process-group kill failed", cleanup))
        try:
            process.wait(timeout=2.0)
        except subprocess.TimeoutExpired as cleanup:
            failures.append(("phase trace terminal reap timed out", cleanup))
        except BaseException as cleanup:
            failures.append(("phase trace terminal reap failed", cleanup))
    if selector is not None:
        try:
            selector.close()
        except BaseException as cleanup:
            failures.append(("phase trace selector close failed", cleanup))
    if process is not None:
        for label, stream in (
                ("stdin", process.stdin), ("stdout", process.stdout),
                ("stderr", process.stderr)):
            if stream is not None:
                try:
                    stream.close()
                except BaseException as cleanup:
                    failures.append((
                        "phase trace {} close failed".format(label), cleanup))
        try:
            alive = process.poll() is None
        except BaseException as cleanup:
            failures.append(("phase trace terminal poll failed", cleanup))
            alive = True
        if alive:
            terminal.append("phase trace process was not reaped")
    runner_api._raise_after_cleanup(primary, failures, terminal)


def _phase_remaining(deadline: float, context: str) -> float:
    try:
        return runner_api._remaining(deadline, context)
    except runner_api.RunnerError as exc:
        fail(str(exc))
    raise AssertionError("unreachable")


def _run_capped_phase_command(
        argv: Sequence[str], deadline: float, context: str) -> bytes:
    """Concurrently drain one trace command under exact byte/deadline caps."""
    if (not isinstance(argv, (list, tuple)) or not argv or
            any(not isinstance(value, str) or not value for value in argv)):
        fail("{} has an invalid argv".format(context))
    process: Optional[subprocess.Popen] = None
    selector: Optional[selectors.BaseSelector] = None
    successful = False
    reaped = False
    output = bytearray()
    diagnostic = bytearray()
    try:
        try:
            process = subprocess.Popen(
                list(argv), stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=0,
                close_fds=True, start_new_session=True)
        except OSError as exc:
            fail("{} failed: {}".format(context, exc))
        if process.stdout is None or process.stderr is None:
            fail("{} did not establish both bounded output pipes".format(
                context))
        selector = selectors.DefaultSelector()
        streams = (
            (process.stdout, "stdout", output, phase_api.MAX_TRACE_BYTES),
            (process.stderr, "stderr", diagnostic,
             MAX_PHASE_TRACE_STDERR_BYTES),
        )
        for stream, label, buffer, cap in streams:
            os.set_blocking(stream.fileno(), False)
            selector.register(
                stream, selectors.EVENT_READ, (label, buffer, cap))
        open_streams = len(streams)
        while open_streams:
            timeout = min(
                1.0, _phase_remaining(deadline, context))
            events = selector.select(timeout)
            if not events:
                continue
            for key, _mask in events:
                label, buffer, cap = key.data
                remaining = cap + 1 - len(buffer)
                if remaining <= 0:
                    fail("{} {} exceeds its exact byte cap".format(
                        context, label))
                try:
                    block = os.read(key.fileobj.fileno(), min(65536, remaining))
                except BlockingIOError:
                    continue
                except OSError as exc:
                    fail("cannot read {} {}: {}".format(context, label, exc))
                if not block:
                    selector.unregister(key.fileobj)
                    open_streams -= 1
                    continue
                buffer.extend(block)
                if len(buffer) > cap:
                    fail("{} {} exceeds its exact byte cap".format(
                        context, label))
        try:
            returncode = process.wait(
                timeout=_phase_remaining(deadline, context))
            reaped = True
        except subprocess.TimeoutExpired as exc:
            fail("{} timed out while reaping: {}".format(context, exc))
        if returncode != 0:
            text = bytes(diagnostic).decode("utf-8", "replace").strip()
            fail("{} exited {}{}".format(
                context, returncode, ": " + text if text else ""))
        if diagnostic:
            fail("{} produced unexpected stderr: {}".format(
                context,
                bytes(diagnostic).decode("utf-8", "replace").strip()))
        successful = True
        result = bytes(output)
    except BaseException as primary:
        _finish_capped_phase_process(
            process, selector, successful, reaped, primary)
        raise
    _finish_capped_phase_process(
        process, selector, successful, reaped, None)
    return result


def _emit_phase_traces(
        contract: Mapping[str, Any], description: Mapping[str, Any],
        output_dir: Path, deadline: float) -> Tuple[Path, bytes, str]:
    try:
        data = _run_capped_phase_command(
            [description["resolved_path"], "--emit-traces",
             "phase-attribution"], deadline, "phase-attribution traces")
        _rows, manifest_sha256 = phase_api.validate_trace_manifest(
            contract, _public_phase_description(description), data)
    except (runner_api.RunnerError, phase_api.PhaseRunnerError) as exc:
        fail(str(exc))
    path = output_dir / PHASE_SOURCE_NAMES["trace"]
    try:
        runner_api._atomic_write_bytes(path, data)
    except runner_api.RunnerError as exc:
        fail(str(exc))
    return path, data, manifest_sha256


def spawn_phase_workers(
        description: Mapping[str, Any], cpus: Sequence[int], deadline: float,
        ) -> List[runner_api.PersistentWorker]:
    """Use the hardened generic worker handoff for the P-capable process."""
    try:
        return runner_api.spawn_workers(description, cpus, deadline)
    except runner_api.RunnerError as exc:
        fail(str(exc))
    raise AssertionError("unreachable")


def _strict_phase_response_validator(
        contract: Mapping[str, Any], description: Mapping[str, Any],
        trace_data: bytes, window_start_ns: int,
        ) -> runner_api.ResponseValidator:
    """Build the exact post-parse validator used after the raw byte cap."""
    try:
        public_description = _public_phase_description(description)
        traces, _digest = phase_api.validate_trace_manifest(
            contract, public_description, trace_data)
        coordinates = phase_api._phase_coordinates(contract)
    except phase_api.PhaseRunnerError as exc:
        fail(str(exc))
    if type(window_start_ns) is not int or window_start_ns <= 0:
        fail("phase worker window start is invalid")

    def validate(
            value: Mapping[str, Any], line: bytes,
            worker: runner_api.PersistentWorker, job: Any) -> int:
        # run_job_batch enforces this bound before calling its JSON parser.
        # Keep the assertion here as a defense against another dispatcher.
        if len(line) > MAX_PHASE_RECORD_BYTES:
            fail("phase native response exceeds its exact raw record cap")
        if not isinstance(job, PhaseJob):
            fail("phase worker response is attached to a non-P command")
        if (not isinstance(value, dict) or
                set(value) != phase_api.NATIVE_RECORD_FIELDS or
                value.get("schema") != phase_api.PHASE_RECORD_SCHEMA):
            fail("phase native response has an unexpected envelope schema")
        integer_fields = (
            "cpu", "ordinal", "worker_pid", "worker_process_start_ticks",
            "started_monotonic_ns", "finished_monotonic_ns",
        )
        if any(type(value.get(field)) is not int for field in integer_fields):
            fail("phase native response has a noncanonical integer field")
        start = value["started_monotonic_ns"]
        end = value["finished_monotonic_ns"]
        if (value["cpu"] != worker.cpu or
                value["ordinal"] != job.ordinal or
                value["worker_pid"] != worker.process.pid or
                value["worker_process_start_ticks"] != worker.start_ticks or
                value.get("worker_binary_sha256") !=
                    description["binary_sha256"] or
                not window_start_ns <= start < end):
            fail("phase native response differs from its live worker/P command")
        try:
            native_api._require_process_predates_window(
                worker.start_ticks, window_start_ns, "phase native worker")
        except native_api.NativeEvidenceError as exc:
            fail(str(exc))
        for field in ("message_sha256", "work_sha256",
                      "worker_binary_sha256"):
            if not _is_sha256(value.get(field)):
                fail("phase native response has an invalid {}".format(field))
        coordinate = coordinates[job.cell]
        if value["message_sha256"] != \
                phase_api._deterministic_message_sha256(coordinate["base"]):
            fail("phase native response substitutes its deterministic message")
        try:
            payload = phase_api._validate_payload(
                value.get("payload"), coordinate, traces[job.cell],
                public_description)
        except phase_api.PhaseRunnerError as exc:
            fail(str(exc))
        expected_work = phase_api._expected_work_sha256(
            payload["cell_sha256"], job.ordinal)
        if value["work_sha256"] != expected_work:
            fail("phase native response work hash differs from its P command")
        return end

    return validate


def _validate_phase_wave_barriers(data: bytes) -> None:
    """Prove all lanes complete wave n before any lane starts wave n+1."""
    try:
        rows = phase_api._parse_jsonl_bytes(
            data, "phase native barrier stream", PHASE_RECORD_COUNT,
            PHASE_NATIVE_STREAM_BYTE_CAP, MAX_PHASE_RECORD_BYTES)
    except phase_api.PhaseRunnerError as exc:
        fail(str(exc))
    by_coordinate: Dict[int, Mapping[str, Any]] = {}
    for row in rows:
        ordinal = row.get("ordinal") if isinstance(row, Mapping) else None
        if (type(ordinal) is not int or ordinal % 2 != 0 or
                not 0 <= ordinal < PHASE_RECORD_COUNT * 2):
            fail("phase barrier stream has an invalid n16 ordinal")
        coordinate = ordinal // 2
        if coordinate in by_coordinate:
            fail("phase barrier stream repeats a coordinate")
        by_coordinate[coordinate] = row
    if set(by_coordinate) != set(range(PHASE_RECORD_COUNT)):
        fail("phase barrier stream omits a coordinate")
    for wave in range(1, PHASE_WAVE_COUNT):
        prior = range((wave - 1) * PHASE_WORKER_COUNT,
                      wave * PHASE_WORKER_COUNT)
        current = range(wave * PHASE_WORKER_COUNT,
                        (wave + 1) * PHASE_WORKER_COUNT)
        prior_end = max(by_coordinate[index].get("finished_monotonic_ns", -1)
                        for index in prior)
        current_start = min(
            by_coordinate[index].get("started_monotonic_ns", -1)
            for index in current)
        if (type(prior_end) is not int or type(current_start) is not int or
                current_start < prior_end):
            fail("phase native intervals violate an exact wave barrier")


def _run_phase_jobs(
        contract: Mapping[str, Any], freeze: Mapping[str, Any],
        description: Mapping[str, Any], trace_data: bytes,
        workers: Sequence[runner_api.PersistentWorker], output_dir: Path,
        window_start_ns: int, deadline: float) -> Tuple[Path, int]:
    if (len(workers) != PHASE_WORKER_COUNT or
            [worker.cpu for worker in workers] != freeze.get("worker_cpus")):
        fail("phase execution does not use its exact frozen worker roster")
    # Construct every possibly-raising validator before creating the staged
    # sink so a validation-setup failure cannot leave a .tmp artifact.
    validator = _strict_phase_response_validator(
        contract, description, trace_data, window_start_ns)
    waves = _phase_job_waves()
    path = output_dir / PHASE_SOURCE_NAMES["native"]
    sink: Optional[runner_api.AtomicLineSink] = None
    maximum_end = 0
    try:
        sink = runner_api.AtomicLineSink(
            path, maximum_bytes=PHASE_NATIVE_STREAM_BYTE_CAP)
        expected_cpus = set(freeze["worker_cpus"])
        for wave_ordinal, jobs in enumerate(waves):
            try:
                end, used = runner_api.run_job_batch(
                    workers, jobs, 0, sink, deadline, validator,
                    maximum_response_bytes=MAX_PHASE_RECORD_BYTES)
            except runner_api.RunnerError as exc:
                fail(str(exc))
            if used != expected_cpus:
                fail("phase wave {} did not use every frozen CPU".format(
                    wave_ordinal))
            maximum_end = max(maximum_end, end)
        sink.publish()
        return path, maximum_end
    finally:
        if sink is not None:
            sink.abort()


def _finish_phase_cleanup(
        workers: Sequence[runner_api.PersistentWorker], clean_shutdown: bool,
        controller_pinned: bool, original_affinity: Set[int],
        primary: Optional[BaseException]) -> None:
    """Exhaust worker and affinity cleanup while preserving control flow."""
    failures: List[Tuple[str, BaseException]] = []
    if workers and not clean_shutdown:
        try:
            runner_api.terminate_workers(workers)
        except BaseException as cleanup:
            failures.append(("phase worker cleanup failed", cleanup))
    if controller_pinned:
        try:
            runner_api._restore_controller_affinity(original_affinity)
        except BaseException as cleanup:
            failures.append(("controller affinity cleanup failed", cleanup))
    if not failures:
        return
    details = []
    if primary is not None:
        details.append("primary failure: {}".format(primary))
    details.extend("{}: {}: {}".format(label, type(exc).__name__, exc)
                   for label, exc in failures)
    diagnostic = PhaseNativeRunnerError("; ".join(details))
    control: Optional[BaseException]
    if primary is not None and not isinstance(primary, Exception):
        control = primary
    else:
        control = next((exc for _label, exc in failures
                        if not isinstance(exc, Exception)), None)
    if control is not None:
        raise control from diagnostic
    if primary is not None:
        raise diagnostic from primary
    raise diagnostic from failures[0][1]


def _read_bounded_regular(path: Path, context: str, byte_limit: int) -> bytes:
    owner: Dict[str, int] = {}
    try:
        data, _fingerprint = runner_api._open_completed_regular_snapshot(
            path, context, getattr(os, "AT_FDCWD", -100), byte_limit,
            owner, "artifact")
        return data
    except runner_api.RunnerError as exc:
        fail(str(exc))
    finally:
        descriptor = owner.pop("artifact", -1)
        if descriptor >= 0:
            runner_api._drain_descriptor_closes((
                ("phase artifact descriptor cleanup failed", descriptor),
            ), sys.exc_info()[1])
    raise AssertionError("unreachable")


def _runtime_worker_map(
        workers: Sequence[runner_api.PersistentWorker],
        description: Mapping[str, Any], verify_live: bool,
        ) -> Tuple[Mapping[int, Tuple[int, int]], List[Mapping[str, Any]]]:
    runtime: Dict[int, Tuple[int, int]] = {}
    rows = []
    for worker in workers:
        identity = (worker.process.pid, worker.start_ticks)
        if worker.cpu in runtime:
            fail("phase runtime worker roster reuses a CPU")
        if verify_live:
            try:
                native_api._verify_live_worker_process(
                    worker.process.pid, worker.cpu, worker.start_ticks,
                    description["binary_sha256"])
            except native_api.NativeEvidenceError as exc:
                fail(str(exc))
        runtime[worker.cpu] = identity
        rows.append({
            "cpu": worker.cpu,
            "pid": worker.process.pid,
            "process_start_ticks": worker.start_ticks,
        })
    if (len(runtime) != PHASE_WORKER_COUNT or
            len({identity[0] for identity in runtime.values()}) !=
                PHASE_WORKER_COUNT):
        fail("phase runtime worker roster is not eight unique live PIDs")
    return runtime, rows


def _validate_thermal(
        sampler: Mapping[str, Any], worker_start_ns: int, worker_end_ns: int,
        worker_cpus: Sequence[int], controller_cpu: int,
        verify_live_sampler: bool) -> Mapping[str, Any]:
    try:
        return native_api._thermal_window(
            sampler, worker_start_ns, worker_end_ns, worker_cpus,
            verify_live_sampler, controller_cpu=controller_cpu)
    except native_api.NativeEvidenceError as exc:
        fail(str(exc))
    raise AssertionError("unreachable")


def _description_from_freeze(freeze: Mapping[str, Any]) -> Mapping[str, Any]:
    return {
        "schema": runner_api.DESCRIPTION_SCHEMA,
        "source_git_commit": freeze["source_git_commit"],
        "binary_sha256": freeze["worker_binary_sha256"],
        "arms": [
            {
                "arm": arm,
                "codec": codec,
                "arm_descriptor_sha256": descriptor,
            }
            for arm, codec, descriptor in phase_api.EXPECTED_ARMS
        ],
    }


def _execution_receipt(
        contract: Mapping[str, Any], freeze: Mapping[str, Any],
        freeze_bytes: bytes, trace_data: bytes, native_data: bytes,
        analysis_bytes: bytes, sampler_bytes: bytes,
        metadata: Mapping[str, Any], thermal: Mapping[str, Any],
        controller_cpu: int) -> Mapping[str, Any]:
    value: Dict[str, Any] = {
        "schema": PHASE_EXECUTION_SCHEMA,
        "status": "complete",
        "source_git_commit": freeze["source_git_commit"],
        "worker_binary_sha256": freeze["worker_binary_sha256"],
        "contract_sha256": contract_api.contract_sha256(contract),
        "recovery_binding_sha256": freeze["recovery_binding_sha256"],
        "freeze_sha256": _hash_bytes(freeze_bytes),
        "trace_manifest_sha256": freeze["trace_manifest_sha256"],
        "trace_file_sha256": _hash_bytes(trace_data),
        "native_stream_sha256": _hash_bytes(native_data),
        "analysis_sha256": _hash_bytes(analysis_bytes),
        "sampler_attestation_sha256": _hash_bytes(sampler_bytes),
        "record_count": metadata["record_count"],
        "worker_cpus": list(metadata["worker_cpus"]),
        "coordinate_cpus": list(metadata["coordinate_cpus"]),
        "workers": list(metadata["workers"]),
        "controller_cpu": controller_cpu,
        "worker_start_monotonic_ns":
            metadata["worker_start_monotonic_ns"],
        "worker_end_monotonic_ns": metadata["worker_end_monotonic_ns"],
        "thermal": dict(thermal),
    }
    value["receipt_sha256"] = contract_api.sha256_json(value)
    return value


def _validate_execution(
        contract: Mapping[str, Any], execution: Any,
        freeze: Mapping[str, Any], snapshots: Mapping[str, bytes],
        metadata: Mapping[str, Any], thermal: Mapping[str, Any],
        ) -> Mapping[str, Any]:
    value = _exact_mapping(
        execution, PHASE_EXECUTION_FIELDS, "phase execution receipt")
    workers = value.get("workers")
    if (value.get("schema") != PHASE_EXECUTION_SCHEMA or
            value.get("status") != "complete" or
            value.get("source_git_commit") != freeze["source_git_commit"] or
            value.get("worker_binary_sha256") !=
                freeze["worker_binary_sha256"] or
            value.get("contract_sha256") !=
                contract_api.contract_sha256(contract) or
            value.get("recovery_binding_sha256") !=
                freeze["recovery_binding_sha256"] or
            value.get("freeze_sha256") != _hash_bytes(snapshots["freeze"]) or
            value.get("trace_manifest_sha256") !=
                freeze["trace_manifest_sha256"] or
            value.get("trace_file_sha256") != _hash_bytes(snapshots["trace"]) or
            value.get("native_stream_sha256") !=
                _hash_bytes(snapshots["native"]) or
            value.get("analysis_sha256") !=
                _hash_bytes(snapshots["analysis"]) or
            value.get("sampler_attestation_sha256") !=
                _hash_bytes(snapshots["sampler"]) or
            value.get("record_count") != PHASE_RECORD_COUNT or
            value.get("record_count") != metadata.get("record_count") or
            value.get("worker_cpus") != freeze["worker_cpus"] or
            value.get("worker_cpus") != metadata.get("worker_cpus") or
            value.get("coordinate_cpus") != freeze["coordinate_cpus"] or
            value.get("coordinate_cpus") != metadata.get("coordinate_cpus") or
            value.get("controller_cpu") != freeze["controller_cpu"] or
            value.get("worker_start_monotonic_ns") !=
                metadata.get("worker_start_monotonic_ns") or
            value.get("worker_end_monotonic_ns") !=
                metadata.get("worker_end_monotonic_ns") or
            contract_api.canonical_json(value.get("thermal")) !=
                contract_api.canonical_json(thermal) or
            value.get("receipt_sha256") !=
                _self_hash(value, "receipt_sha256") or
            not isinstance(workers, list) or
            len(workers) != PHASE_WORKER_COUNT or
            any(not isinstance(worker, dict) or
                set(worker) != PHASE_WORKER_FIELDS for worker in workers) or
            workers != metadata.get("workers")):
        fail("phase execution receipt differs from validated native evidence")
    return value


def _validate_summary(
        contract: Mapping[str, Any], resolved: Path,
        summary: Any, freeze: Mapping[str, Any],
        execution: Mapping[str, Any], analysis: Mapping[str, Any],
        snapshots: Mapping[str, bytes],
        recovery_binding: Mapping[str, Any]) -> Mapping[str, Any]:
    value = _exact_mapping(summary, PHASE_SUMMARY_FIELDS, "phase run summary")
    thermal = execution["thermal"]
    decision = analysis.get("decision")
    expected = {
        "schema": PHASE_SUMMARY_SCHEMA,
        "status": "complete",
        "output_dir": str(resolved),
        "source_git_commit": freeze["source_git_commit"],
        "worker_binary_sha256": freeze["worker_binary_sha256"],
        "contract_sha256": contract_api.contract_sha256(contract),
        "candidate_arm": phase_api.LEFT_ARM,
        "candidate_arm_descriptor_sha256":
            freeze["candidate_arm_descriptor_sha256"],
        "timing_proxy_witness_sha256":
            freeze["timing_proxy_witness_sha256"],
        "recovery_binding_sha256": recovery_binding["binding_sha256"],
        "recovery_run_summary_sha256":
            recovery_binding["recovery_run_summary_sha256"],
        "recovery_freeze_sha256":
            recovery_binding["recovery_freeze_sha256"],
        "recovery_execution_receipt_sha256":
            recovery_binding["recovery_execution_receipt_sha256"],
        "phase_freeze_sha256": _hash_bytes(snapshots["freeze"]),
        "trace_manifest_sha256": freeze["trace_manifest_sha256"],
        "trace_file_sha256": _hash_bytes(snapshots["trace"]),
        "native_stream_sha256": _hash_bytes(snapshots["native"]),
        "analysis_sha256": _hash_bytes(snapshots["analysis"]),
        "sampler_attestation_sha256": _hash_bytes(snapshots["sampler"]),
        "execution_receipt_sha256": _hash_bytes(snapshots["execution"]),
        "phase_records": PHASE_RECORD_COUNT,
        "controller_cpu": freeze["controller_cpu"],
        "worker_cpus": freeze["worker_cpus"],
        "thermal_samples": thermal.get("sample_count")
            if isinstance(thermal, Mapping) else None,
        "cpu_tctl_max_millic": thermal.get("cpu_tctl_max_millic")
            if isinstance(thermal, Mapping) else None,
        "dimm_max_millic": thermal.get("dimm_max_millic")
            if isinstance(thermal, Mapping) else None,
        "decision_outcome": decision.get("outcome")
            if isinstance(decision, Mapping) else None,
    }
    if (value.get("summary_sha256") !=
            _self_hash(value, "summary_sha256") or
            any(contract_api.canonical_json(value.get(field)) !=
                contract_api.canonical_json(expected_value)
                for field, expected_value in expected.items())):
        fail("phase run summary differs from validated phase evidence")
    return value


def _validate_phase_snapshots(
        contract: Mapping[str, Any], resolved: Path,
        directory_identity: Tuple[int, int],
        snapshots: Mapping[str, bytes],
        recovery_binding: Mapping[str, Any],
        verify_live_sampler: bool = False) -> Mapping[str, Any]:
    """Semantically rederive one complete exact seven-artifact snapshot."""
    if (set(snapshots) != set(PHASE_SOURCE_NAMES) or
            any(type(data) is not bytes for data in snapshots.values())):
        fail("completed phase snapshot is not the exact seven-artifact bundle")
    total = 0
    for key, name in PHASE_SOURCE_NAMES.items():
        data = snapshots[key]
        if len(data) > PHASE_ARTIFACT_BYTE_LIMITS[name]:
            fail("completed phase {} exceeds its exact byte cap".format(key))
        total += len(data)
    if total > MAX_PHASE_BUNDLE_BYTES:
        fail("completed phase evidence exceeds its aggregate byte cap")

    summary = _parse_object_bytes(snapshots["summary"], "phase run summary")
    freeze = _parse_object_bytes(snapshots["freeze"], "phase freeze")
    analysis = _parse_object_bytes(
        snapshots["analysis"], "phase analysis")
    execution = _parse_object_bytes(
        snapshots["execution"], "phase execution receipt")
    sampler = _parse_object_bytes(
        snapshots["sampler"], "phase sampler attestation")
    freeze = _validate_phase_freeze(
        contract, freeze, recovery_binding)
    description = _description_from_freeze(freeze)
    try:
        _traces, trace_manifest_sha256 = phase_api.validate_trace_manifest(
            contract, description, snapshots["trace"])
    except phase_api.PhaseRunnerError as exc:
        fail(str(exc))
    if (trace_manifest_sha256 != freeze["trace_manifest_sha256"] or
            _hash_bytes(snapshots["trace"]) != freeze["trace_file_sha256"]):
        fail("phase trace bytes differ from their frozen identity")

    if not isinstance(execution, Mapping):
        fail("phase execution receipt is malformed")
    raw_workers = execution.get("workers")
    if (not isinstance(raw_workers, list) or
            len(raw_workers) != PHASE_WORKER_COUNT or
            any(not isinstance(worker, dict) or
                set(worker) != PHASE_WORKER_FIELDS for worker in raw_workers)):
        fail("phase execution receipt has an invalid worker roster")
    runtime_workers: Dict[int, Tuple[int, int]] = {}
    for worker in raw_workers:
        cpu = worker.get("cpu")
        pid = worker.get("pid")
        ticks = worker.get("process_start_ticks")
        if (type(cpu) is not int or cpu < 0 or type(pid) is not int or
                pid <= 0 or type(ticks) is not int or ticks <= 0 or
                cpu in runtime_workers):
            fail("phase execution worker identity is noncanonical")
        runtime_workers[cpu] = (pid, ticks)
    start_ns = execution.get("worker_start_monotonic_ns")
    end_ns = execution.get("worker_end_monotonic_ns")
    if (type(start_ns) is not int or start_ns <= 0 or
            type(end_ns) is not int or end_ns <= start_ns):
        fail("phase execution receipt has an invalid worker window")
    try:
        payloads, metadata = phase_api.validate_native_records(
            contract, description, snapshots["trace"], snapshots["native"],
            freeze["coordinate_cpus"], runtime_workers, start_ns, end_ns)
        expected_analysis = phase_api.build_phase_analysis(payloads)
    except phase_api.PhaseRunnerError as exc:
        fail(str(exc))
    _validate_phase_wave_barriers(snapshots["native"])
    if contract_api.canonical_json(analysis) != \
            contract_api.canonical_json(expected_analysis):
        fail("phase analysis differs from recomputed native phase attribution")
    thermal = _validate_thermal(
        sampler, metadata["worker_start_monotonic_ns"],
        metadata["worker_end_monotonic_ns"], metadata["worker_cpus"],
        freeze["controller_cpu"], verify_live_sampler)
    execution = _validate_execution(
        contract, execution, freeze, snapshots, metadata, thermal)
    summary = _validate_summary(
        contract, resolved, summary, freeze, execution, analysis, snapshots,
        recovery_binding)
    return {
        "directory": str(resolved),
        "directory_identity": directory_identity,
        "summary": summary,
        "freeze": freeze,
        "execution": execution,
        "analysis": analysis,
        "sampler": sampler,
        "payloads": payloads,
        "recovery_binding": recovery_binding,
    }


def _open_pinned_phase_bundle(
        source_names: Mapping[str, str], directory_fd: int,
        descriptors: Dict[str, int],
        ) -> Tuple[Dict[str, bytes],
                   Dict[str, runner_api.CompletedFingerprint]]:
    if (not isinstance(descriptors, dict) or descriptors or
            not source_names or
            not set(source_names).issubset(PHASE_SOURCE_NAMES) or
            any(PHASE_SOURCE_NAMES[key] != name
                for key, name in source_names.items())):
        fail("pinned completed phase bundle is not an exact known subset")
    snapshots: Dict[str, bytes] = {}
    fingerprints: Dict[str, runner_api.CompletedFingerprint] = {}
    total = 0
    try:
        for key, name in source_names.items():
            try:
                data, fingerprint = runner_api._open_completed_regular_snapshot(
                    Path(name), "completed phase {}".format(key), directory_fd,
                    PHASE_ARTIFACT_BYTE_LIMITS[name], descriptors, key)
            except runner_api.RunnerError as exc:
                fail(str(exc))
            snapshots[key] = data
            fingerprints[key] = fingerprint
            total += len(data)
            if total > MAX_PHASE_BUNDLE_BYTES:
                fail("completed phase evidence exceeds its aggregate byte cap")
        return snapshots, fingerprints
    except BaseException as primary:
        pending = [
            ("phase dependency {} descriptor cleanup failed".format(key), fd)
            for key, fd in descriptors.items()
        ]
        descriptors.clear()
        runner_api._drain_descriptor_closes(pending, primary)
        raise


def _reread_pinned_phase_bundle(
        source_names: Mapping[str, str], directory_fd: int,
        descriptors: Mapping[str, int],
        expected: Mapping[str, runner_api.CompletedFingerprint],
        ) -> Tuple[Dict[str, bytes],
                   Dict[str, runner_api.CompletedFingerprint]]:
    if set(descriptors) != set(source_names) or set(expected) != set(source_names):
        fail("pinned completed phase bundle is incomplete")
    snapshots: Dict[str, bytes] = {}
    fingerprints: Dict[str, runner_api.CompletedFingerprint] = {}
    total = 0
    for key, name in source_names.items():
        try:
            data, fingerprint = runner_api._read_completed_descriptor_bytes(
                descriptors[key], "completed phase {}".format(key),
                PHASE_ARTIFACT_BYTE_LIMITS[name])
            named = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
        except runner_api.RunnerError as exc:
            fail(str(exc))
        except OSError as exc:
            fail("cannot terminally inspect completed phase {}: {}".format(
                key, exc))
        if (not stat.S_ISREG(named.st_mode) or
                fingerprint != expected[key] or
                runner_api._completed_fingerprint(named) != fingerprint):
            fail("completed phase {} changed during semantic validation".format(
                key))
        snapshots[key] = data
        fingerprints[key] = fingerprint
        total += len(data)
        if total > MAX_PHASE_BUNDLE_BYTES:
            fail("completed phase evidence exceeds its aggregate byte cap")
    return snapshots, fingerprints


def _require_pinned_phase_unchanged(
        source_names: Mapping[str, str], directory_fd: int,
        descriptors: Mapping[str, int],
        expected: Mapping[str, runner_api.CompletedFingerprint]) -> None:
    if set(descriptors) != set(source_names) or set(expected) != set(source_names):
        fail("pinned completed phase bundle is incomplete")
    for key, name in source_names.items():
        try:
            retained = os.fstat(descriptors[key])
            named = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
        except OSError as exc:
            fail("cannot recheck completed phase {}: {}".format(key, exc))
        if (not stat.S_ISREG(retained.st_mode) or
                not stat.S_ISREG(named.st_mode) or
                runner_api._completed_fingerprint(retained) != expected[key] or
                runner_api._completed_fingerprint(named) != expected[key]):
            fail("completed phase {} changed before publication".format(key))


def _require_exact_phase_names(directory_fd: int, complete: bool) -> None:
    expected = set(PHASE_SOURCE_NAMES.values()) if complete else \
        set(PHASE_DEPENDENCY_NAMES.values())
    names: Set[str] = set()
    try:
        # This helper also runs after hard-wall teardown.  Stream at most one
        # entry beyond the exact six/seven-name roster; never materialize an
        # attacker-controlled directory with listdir().
        with os.scandir(directory_fd) as entries:
            for entry in entries:
                if len(names) >= len(expected):
                    fail("phase directory contains too many artifacts")
                name = entry.name
                if (not isinstance(name, str) or not name or
                        name in names or
                        not entry.is_file(follow_symlinks=False)):
                    fail("phase directory has an invalid artifact entry")
                names.add(name)
    except OSError as exc:
        fail("cannot enumerate completed phase artifact roster: {}".format(exc))
    if names != expected:
        fail("phase directory does not contain its exact artifact roster")


def _open_phase_directory(
        phase_dir: Path, descriptor_holder: List[int],
        expected_identity: Optional[Tuple[int, int]] = None,
        ) -> Tuple[int, Path, Tuple[int, int]]:
    if not isinstance(descriptor_holder, list) or descriptor_holder != [-1]:
        fail("phase directory descriptor holder is invalid")
    marker = phase_dir / PHASE_SOURCE_NAMES["summary"]
    try:
        opened_fd, identity = runner_api._open_completion_parent(
            marker, descriptor_holder)
        resolved = phase_dir.resolve(strict=True)
        runner_api._verify_completed_directory_path(resolved, identity)
    except runner_api.RunnerError as exc:
        fail(str(exc))
    if expected_identity is not None and identity != expected_identity:
        owned = descriptor_holder[0]
        descriptor_holder[0] = -1
        runner_api._drain_descriptor_closes((
            ("substituted phase directory cleanup failed", owned),
        ), None)
        fail("completed phase directory changed before terminal read")
    return opened_fd, resolved, identity


def load_completed_phase_screen(
        contract: Mapping[str, Any], phase_dir: Path, recovery_dir: Path,
        verify_live_sampler: bool = False) -> Mapping[str, Any]:
    """Pin, semantically derive twice, and recovery-rebind a phase screen."""
    _campaign, recovery_binding = _load_recovery_binding(
        contract, recovery_dir)
    directory_holder = [-1]
    descriptors: Dict[str, int] = {}
    try:
        opened_fd, resolved, identity = _open_phase_directory(
            phase_dir, directory_holder)
        directory_fd = directory_holder[0]
        if opened_fd != directory_fd:
            fail("phase directory ownership handoff failed")
        _require_exact_phase_names(directory_fd, complete=True)
        first_snapshots, first_fingerprints = _open_pinned_phase_bundle(
            PHASE_SOURCE_NAMES, directory_fd, descriptors)
        runner_api._verify_completed_directory_path(resolved, identity)
        first = _validate_phase_snapshots(
            contract, resolved, identity, first_snapshots, recovery_binding,
            verify_live_sampler)
        second_snapshots, second_fingerprints = _reread_pinned_phase_bundle(
            PHASE_SOURCE_NAMES, directory_fd, descriptors,
            first_fingerprints)
        runner_api._verify_completed_directory_path(resolved, identity)
        second = _validate_phase_snapshots(
            contract, resolved, identity, second_snapshots, recovery_binding,
            verify_live_sampler)
        if (first_snapshots != second_snapshots or
                first_fingerprints != second_fingerprints or
                contract_api.canonical_json(first) !=
                    contract_api.canonical_json(second)):
            fail("completed phase evidence changed between semantic reads")
        _terminal_campaign, terminal_binding = _load_recovery_binding(
            contract, recovery_dir)
        if contract_api.canonical_json(terminal_binding) != \
                contract_api.canonical_json(recovery_binding):
            fail("completed recovery input changed during phase validation")
        final_snapshots, final_fingerprints = _reread_pinned_phase_bundle(
            PHASE_SOURCE_NAMES, directory_fd, descriptors,
            first_fingerprints)
        if (final_snapshots != first_snapshots or
                final_fingerprints != first_fingerprints):
            fail("completed phase evidence changed during recovery revalidation")
        _require_exact_phase_names(directory_fd, complete=True)
        # This retained-inode/public-name check is intentionally the last
        # operation before return, closing the recovery-reload substitution
        # window without another fallible dependency read afterward.
        _require_pinned_phase_unchanged(
            PHASE_SOURCE_NAMES, directory_fd, descriptors, first_fingerprints)
        return first
    finally:
        pending = [
            ("phase {} descriptor cleanup failed".format(key), descriptor)
            for key, descriptor in descriptors.items()
        ]
        descriptors.clear()
        owned_directory = directory_holder[0]
        directory_holder[0] = -1
        pending.append(("phase directory descriptor cleanup failed",
                        owned_directory))
        runner_api._drain_descriptor_closes(pending, sys.exc_info()[1])


def _phase_summary(
        contract: Mapping[str, Any], output_dir: Path,
        freeze: Mapping[str, Any], execution: Mapping[str, Any],
        analysis: Mapping[str, Any], recovery_binding: Mapping[str, Any],
        freeze_bytes: bytes, trace_data: bytes, native_data: bytes,
        analysis_bytes: bytes, execution_bytes: bytes,
        sampler_bytes: bytes) -> Mapping[str, Any]:
    thermal = execution["thermal"]
    value: Dict[str, Any] = {
        "schema": PHASE_SUMMARY_SCHEMA,
        "status": "complete",
        "output_dir": str(output_dir),
        "source_git_commit": freeze["source_git_commit"],
        "worker_binary_sha256": freeze["worker_binary_sha256"],
        "contract_sha256": contract_api.contract_sha256(contract),
        "candidate_arm": freeze["candidate_arm"],
        "candidate_arm_descriptor_sha256":
            freeze["candidate_arm_descriptor_sha256"],
        "timing_proxy_witness_sha256":
            freeze["timing_proxy_witness_sha256"],
        "recovery_binding_sha256": recovery_binding["binding_sha256"],
        "recovery_run_summary_sha256":
            recovery_binding["recovery_run_summary_sha256"],
        "recovery_freeze_sha256":
            recovery_binding["recovery_freeze_sha256"],
        "recovery_execution_receipt_sha256":
            recovery_binding["recovery_execution_receipt_sha256"],
        "phase_freeze_sha256": _hash_bytes(freeze_bytes),
        "trace_manifest_sha256": freeze["trace_manifest_sha256"],
        "trace_file_sha256": _hash_bytes(trace_data),
        "native_stream_sha256": _hash_bytes(native_data),
        "analysis_sha256": _hash_bytes(analysis_bytes),
        "sampler_attestation_sha256": _hash_bytes(sampler_bytes),
        "execution_receipt_sha256": _hash_bytes(execution_bytes),
        "phase_records": PHASE_RECORD_COUNT,
        "controller_cpu": freeze["controller_cpu"],
        "worker_cpus": list(freeze["worker_cpus"]),
        "thermal_samples": thermal["sample_count"],
        "cpu_tctl_max_millic": thermal["cpu_tctl_max_millic"],
        "dimm_max_millic": thermal["dimm_max_millic"],
        "decision_outcome": analysis["decision"]["outcome"],
    }
    value["summary_sha256"] = contract_api.sha256_json(value)
    return value


def _commit_completed_phase_screen(
        contract: Mapping[str, Any], output_dir: Path,
        recovery_dir: Path, recovery_binding: Mapping[str, Any],
        summary: Mapping[str, Any], source_commit: str, deadline: float,
        finish_hard_wall: Optional[Callable[[], None]] = None) -> None:
    """Double-derive pinned dependencies, then link one unnamed summary."""
    summary_path = output_dir / PHASE_SOURCE_NAMES["summary"]
    summary_bytes = _canonical_object_bytes(summary)
    if len(summary_bytes) > PHASE_ARTIFACT_BYTE_LIMITS["run-summary.json"]:
        fail("phase completion marker exceeds its exact byte cap")
    parent_holder = [-1]
    summary_holder = [-1]
    dependency_fds: Dict[str, int] = {}
    try:
        try:
            opened_parent, parent_identity = runner_api._open_completion_parent(
                summary_path, parent_holder)
            parent_fd = parent_holder[0]
            if opened_parent != parent_fd:
                fail("phase completion parent ownership handoff failed")
            runner_api._require_completion_parent_descriptor(
                parent_fd, parent_identity)
            runner_api._verify_completed_directory_path(
                output_dir, parent_identity)
            runner_api._require_completion_marker_absent(summary_path, parent_fd)
            summary_fingerprint = runner_api._open_unnamed_completion_marker(
                parent_fd, summary_bytes, summary_holder)
        except runner_api.RunnerError as exc:
            fail(str(exc))
        summary_fd = summary_holder[0]
        _require_exact_phase_names(parent_fd, complete=False)
        dependency_snapshots, dependency_fingerprints = \
            _open_pinned_phase_bundle(
                PHASE_DEPENDENCY_NAMES, parent_fd, dependency_fds)
        if len(summary_bytes) + sum(
                len(data) for data in dependency_snapshots.values()) > \
                MAX_PHASE_BUNDLE_BYTES:
            fail("completed phase evidence exceeds its aggregate byte cap")
        prospective = dict(dependency_snapshots)
        prospective["summary"] = summary_bytes
        resolved = output_dir.resolve(strict=True)
        first = _validate_phase_snapshots(
            contract, resolved, parent_identity, prospective,
            recovery_binding, verify_live_sampler=False)

        try:
            final_summary_bytes, final_summary_fingerprint = \
                runner_api._read_unnamed_completion_marker(
                    summary_fd,
                    PHASE_ARTIFACT_BYTE_LIMITS["run-summary.json"])
        except runner_api.RunnerError as exc:
            fail(str(exc))
        final_dependencies, final_fingerprints = \
            _reread_pinned_phase_bundle(
                PHASE_DEPENDENCY_NAMES, parent_fd, dependency_fds,
                dependency_fingerprints)
        final_prospective = dict(final_dependencies)
        final_prospective["summary"] = final_summary_bytes
        second = _validate_phase_snapshots(
            contract, resolved, parent_identity, final_prospective,
            recovery_binding, verify_live_sampler=False)
        if (final_summary_bytes != summary_bytes or
                final_summary_fingerprint != summary_fingerprint or
                final_dependencies != dependency_snapshots or
                final_fingerprints != dependency_fingerprints or
                contract_api.canonical_json(first) !=
                    contract_api.canonical_json(second) or
                contract_api.canonical_json(first.get("summary")) !=
                    contract_api.canonical_json(summary)):
            fail("prospective phase evidence changed between semantic reads")

        try:
            if runner_api._git_head(deadline) != source_commit:
                fail("codec source commit changed before phase publication")
        except runner_api.RunnerError as exc:
            fail(str(exc))
        # This is the terminal input-transaction reopen.  Source validation is
        # deliberately before it, leaving recovery as the final external
        # semantic input read before local pinned checks and publication.
        _terminal_campaign, terminal_binding = _load_recovery_binding(
            contract, recovery_dir)
        if contract_api.canonical_json(terminal_binding) != \
                contract_api.canonical_json(recovery_binding):
            fail("completed recovery input changed before phase publication")
        try:
            runner_api._require_completion_parent_descriptor(
                parent_fd, parent_identity)
            runner_api._verify_completed_directory_path(
                output_dir, parent_identity)
            runner_api._require_completion_marker_absent(summary_path, parent_fd)
        except runner_api.RunnerError as exc:
            fail(str(exc))

        if finish_hard_wall is not None:
            finish_hard_wall()
        try:
            marker_info = os.fstat(summary_fd)
        except OSError as exc:
            fail("cannot recheck unnamed phase completion marker: {}".format(
                exc))
        if runner_api._completed_fingerprint(marker_info) != summary_fingerprint:
            fail("prospective phase completion marker changed before publication")
        _require_pinned_phase_unchanged(
            PHASE_DEPENDENCY_NAMES, parent_fd, dependency_fds,
            dependency_fingerprints)
        _require_exact_phase_names(parent_fd, complete=False)
        try:
            runner_api._require_completion_parent_descriptor(
                parent_fd, parent_identity)
            runner_api._verify_completed_directory_path(
                output_dir, parent_identity)
            runner_api._require_completion_marker_absent(summary_path, parent_fd)
            runner_api._link_unnamed_completion_marker(
                summary_path, summary_fd, parent_fd, parent_identity)
        except runner_api.RunnerError as exc:
            fail(str(exc))
    finally:
        pending = [
            ("phase dependency {} descriptor cleanup failed".format(key), fd)
            for key, fd in dependency_fds.items()
        ]
        dependency_fds.clear()
        pending.extend((
            ("unnamed phase summary descriptor cleanup failed",
             summary_holder[0]),
            ("phase completion parent descriptor cleanup failed",
             parent_holder[0]),
        ))
        summary_holder[0] = -1
        parent_holder[0] = -1
        runner_api._drain_descriptor_closes(pending, sys.exc_info()[1])


def run_phase_screen(
        args: argparse.Namespace,
        finish_hard_wall: Optional[Callable[[], None]] = None,
        ) -> Mapping[str, Any]:
    if (not math.isfinite(args.deadline_seconds) or
            not 0.0 < args.deadline_seconds <= runner_api.MAX_WALL_SECONDS):
        fail("--deadline-seconds must be in (0,7200]")
    deadline = time.monotonic() + args.deadline_seconds
    contract = contract_api.load_contract(args.contract)
    _campaign, recovery_binding = _load_recovery_binding(
        contract, args.recovery_dir)
    try:
        description = runner_api.describe_worker(args.worker, deadline)
        original_affinity = set(os.sched_getaffinity(0))
    except (runner_api.RunnerError, AttributeError, OSError) as exc:
        fail(str(exc))
    explicit = runner_api.parse_cpu_list(args.cpus) \
        if args.cpus is not None else None
    try:
        cpus, controller_cpu = runner_api.select_cpu_layout(
            args.sampler_cpu, explicit, args.controller_cpu,
            affinity=original_affinity)
        runner_api._preflight_sampler(
            args.sampler_pid, args.sampler_cpu, args.sampler_script,
            args.sampler_csv)
        source_commit = runner_api._git_head(deadline)
        runner_api._require_worker_source_commit(description, source_commit)
    except runner_api.RunnerError as exc:
        fail(str(exc))
    if (source_commit != recovery_binding["source_git_commit"] or
            description["binary_sha256"] !=
                recovery_binding["worker_binary_sha256"] or
            description["arms"][2]["arm_descriptor_sha256"] !=
                recovery_binding["candidate_arm_descriptor_sha256"]):
        fail("phase worker is not the recovery campaign's exact source/binary/Two07")
    try:
        output_dir = runner_api._create_output_dir(args.output_dir)
    except runner_api.RunnerError as exc:
        fail(str(exc))
    trace_path, trace_data, trace_manifest_sha256 = _emit_phase_traces(
        contract, description, output_dir, deadline)
    try:
        if runner_api._git_head(deadline) != source_commit:
            fail("codec source commit changed before the phase freeze")
    except runner_api.RunnerError as exc:
        fail(str(exc))
    freeze = _phase_freeze(
        contract, description, recovery_binding, cpus, controller_cpu,
        trace_manifest_sha256, _hash_bytes(trace_data))
    freeze_path = output_dir / PHASE_SOURCE_NAMES["freeze"]
    try:
        runner_api._atomic_write_object(freeze_path, freeze)
    except runner_api.RunnerError as exc:
        fail(str(exc))
    freeze_bytes = _read_bounded_regular(
        freeze_path, "phase freeze", PHASE_ARTIFACT_BYTE_LIMITS[
            PHASE_SOURCE_NAMES["freeze"]])

    workers: List[runner_api.PersistentWorker] = []
    clean_shutdown = False
    controller_pinned = False
    completed: Optional[Tuple[Mapping[str, Any], Mapping[str, Any]]] = None
    try:
        workers = spawn_phase_workers(description, cpus, deadline)
        controller_pinned = True
        runner_api._pin_controller(controller_cpu)
        window_start_ns = runner_api.choose_new_sampler_start(
            args.sampler_csv, deadline)
        native_path, maximum_worker_end = _run_phase_jobs(
            contract, freeze, description, trace_data, workers, output_dir,
            window_start_ns, deadline)
        window_end_ns = runner_api._wait_for_sampler_sample(
            args.sampler_csv, deadline, at_or_after_ns=maximum_worker_end)
        sampler_path = output_dir / PHASE_SOURCE_NAMES["sampler"]
        try:
            native_api.write_sampler_attestation(
                args.sampler_pid, args.sampler_cpu, args.sampler_script,
                args.sampler_csv, window_start_ns, window_end_ns,
                sampler_path)
        except native_api.NativeEvidenceError as exc:
            fail(str(exc))
        runtime_workers, runtime_rows = _runtime_worker_map(
            workers, description, verify_live=True)
        native_data = _read_bounded_regular(
            native_path, "phase native stream", PHASE_NATIVE_STREAM_BYTE_CAP)
        try:
            payloads, metadata = phase_api.validate_native_records(
                contract, _public_phase_description(description),
                trace_data, native_data,
                freeze["coordinate_cpus"], runtime_workers,
                window_start_ns, window_end_ns)
            analysis = phase_api.build_phase_analysis(payloads)
        except phase_api.PhaseRunnerError as exc:
            fail(str(exc))
        if metadata["workers"] != runtime_rows:
            fail("phase native metadata differs from the live worker roster")
        _validate_phase_wave_barriers(native_data)
        analysis_path = output_dir / PHASE_SOURCE_NAMES["analysis"]
        try:
            runner_api._atomic_write_object(analysis_path, analysis)
        except runner_api.RunnerError as exc:
            fail(str(exc))
        analysis_bytes = _read_bounded_regular(
            analysis_path, "phase analysis", PHASE_ARTIFACT_BYTE_LIMITS[
                PHASE_SOURCE_NAMES["analysis"]])
        sampler_bytes = _read_bounded_regular(
            sampler_path, "phase sampler attestation",
            PHASE_ARTIFACT_BYTE_LIMITS[PHASE_SOURCE_NAMES["sampler"]])
        sampler = _parse_object_bytes(
            sampler_bytes, "phase sampler attestation")
        thermal = _validate_thermal(
            sampler, metadata["worker_start_monotonic_ns"],
            metadata["worker_end_monotonic_ns"], metadata["worker_cpus"],
            controller_cpu, verify_live_sampler=True)
        execution = _execution_receipt(
            contract, freeze, freeze_bytes, trace_data, native_data,
            analysis_bytes, sampler_bytes, metadata, thermal,
            controller_cpu)
        execution_path = output_dir / PHASE_SOURCE_NAMES["execution"]
        try:
            runner_api._atomic_write_object(execution_path, execution)
        except runner_api.RunnerError as exc:
            fail(str(exc))
        execution_bytes = _read_bounded_regular(
            execution_path, "phase execution receipt",
            PHASE_ARTIFACT_BYTE_LIMITS[PHASE_SOURCE_NAMES["execution"]])
        runner_api.quit_workers(workers, deadline)
        clean_shutdown = True
        summary = _phase_summary(
            contract, output_dir, freeze, execution, analysis,
            recovery_binding, freeze_bytes, trace_data, native_data,
            analysis_bytes, execution_bytes, sampler_bytes)
        completed = (summary, execution)
    finally:
        _finish_phase_cleanup(
            workers, clean_shutdown, controller_pinned, original_affinity,
            sys.exc_info()[1])
    if completed is None:
        fail("phase screen reached cleanup without complete evidence")
    summary, _execution = completed

    # Q/reap and exact controller-affinity restoration deliberately precede
    # both terminal source validation and terminal recovery-input reopening.
    try:
        if runner_api._git_head(deadline) != source_commit:
            fail("codec source commit changed during the phase screen")
    except runner_api.RunnerError as exc:
        fail(str(exc))
    _terminal_campaign, terminal_binding = _load_recovery_binding(
        contract, args.recovery_dir)
    if contract_api.canonical_json(terminal_binding) != \
            contract_api.canonical_json(recovery_binding):
        fail("completed recovery input changed during the phase screen")
    _commit_completed_phase_screen(
        contract, output_dir, args.recovery_dir, recovery_binding, summary,
        source_commit, deadline, finish_hard_wall)
    return summary


def main(argv: Sequence[str] = ()) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--contract", type=Path, default=contract_api.DEFAULT_CONTRACT)
    parser.add_argument(
        "--worker", type=Path,
        default=Path("build/wirehair_wh2_contract_worker"))
    parser.add_argument("--recovery-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--sampler-pid", type=int, required=True)
    parser.add_argument("--sampler-cpu", type=int, required=True)
    parser.add_argument("--sampler-script", type=Path, required=True)
    parser.add_argument("--sampler-csv", type=Path, required=True)
    parser.add_argument(
        "--cpus", help="strictly increasing eight logical worker CPUs")
    parser.add_argument("--controller-cpu", type=int)
    parser.add_argument(
        "--deadline-seconds", type=float,
        default=runner_api.MAX_WALL_SECONDS)
    args = parser.parse_args(argv if argv else None)
    try:
        with runner_api._hard_wall(
                args.deadline_seconds) as finish_hard_wall:
            summary = run_phase_screen(args, finish_hard_wall)
    except (PhaseNativeRunnerError, runner_api.RunnerError,
            recovery_api.RecoveryRunnerError, phase_api.PhaseRunnerError,
            native_api.NativeEvidenceError, contract_api.ContractError,
            OSError, UnicodeError) as exc:
        print("wh2 native phase attribution: {}".format(exc), file=sys.stderr)
        return 1
    print(contract_api.canonical_json(summary))
    return 0


if __name__ == "__main__":
    sys.exit(main())
