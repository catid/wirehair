#!/usr/bin/env python3
"""Run and combine the three bounded WH2 recovery-only candidate screens.

Each ``run`` invocation freezes exactly one closed candidate with the WH2-head,
Wirehair1, and uniform-raw D12/H12 controls, emits no timing work, and publishes
an independently validated native recovery execution receipt.  ``combine``
revalidates three completed runs and produces a six-arm *logical* recovery
ledger and summary.
The combination deliberately is not, and must not be presented as, one native
execution receipt: its provenance is the three input receipt hashes.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import os
from pathlib import Path
import re
import stat
import subprocess
import sys
import tempfile
import time
from typing import Any, Dict, List, Mapping, Optional, Sequence, Set, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_benchmark_contract as contract_api
import wh2_native_short_screen as native_api
import wh2_run_native_short_screen as runner_api


CANDIDATE_SPECS = (
    ("d12-h11-periodic", "wirehair2_raw_d12_h11_periodic",
     "80f57b83eac9b8e3a19d8235cc067e39990980510e46588ddefa16f9561e1c38"),
    ("d12-h13-periodic", "wirehair2_raw_d12_h13_periodic",
     "0eb3aef0602b5e7de15c822de84a5dbfc5dfdd99b76fbfd41538f7a13248c3a5"),
    ("d13-h12-periodic", "wirehair2_raw_d13_h12_periodic",
     "2dc244661b3b073569319377ee3e55333a82ddad7bd328e1b0fef67395174614"),
)
CANDIDATE_BY_ID = {
    candidate_id: (arm, descriptor)
    for candidate_id, arm, descriptor in CANDIDATE_SPECS
}
CONTROL_ARMS = (
    ("wirehair2_head", "wirehair2_certified"),
    ("wirehair1", "wirehair1"),
)
RAW_CONTROL_ARM = (
    "wirehair2_raw_d12_h12_periodic", "wirehair2_experiment")
RAW_CONTROL_DESCRIPTOR_SHA256 = \
    "0550e0ed0c62d5491ff6915652fd96ed25f3c7782462da8c551636ec2e0294dd"
RAW_SEED_BASIS = "uniform-raw-v1"
RAW_SEED_SCHEDULE_SHA256 = \
    "90a98a3db207852dabdf5fb27573ef48bce52e0228cee4e291d96fa44ed509a7"
CONTROL_DESCRIPTOR_SHA256S = (
    "4cafe27a8fb388ca9a4249b2c279b1406e7a0a86bcf14e98246988c7c503fa7a",
    "d5a24d404e69efeb439907cd8271eba98d6af86b58efe159a820fb7aea08883d",
)

CAMPAIGN_SUMMARY_SCHEMA = "wirehair.wh2.native-recovery-screen-run.v1"
COMBINATION_SUMMARY_SCHEMA = \
    "wirehair.wh2.logical-recovery-combination.v1"
CAMPAIGN_SUMMARY_FIELDS = frozenset((
    "schema", "status", "output_dir", "candidate_id", "candidate_arm",
    "source_git_commit", "contract_sha256", "domain_sha256",
    "trace_manifest_sha256", "worker_binary_sha256", "controller_cpu",
    "worker_cpus", "recovery_records", "recovery_freeze_sha256",
    "recovery_result_sha256", "recovery_execution_receipt_sha256",
    "thermal_samples", "cpu_tctl_max_millic", "dimm_max_millic",
    "seed_basis", "seed_schedule_sha256", "summary_sha256",
))
COMBINATION_SUMMARY_FIELDS = frozenset((
    "schema", "status", "artifact_kind", "is_execution_receipt", "phase",
    "contract_sha256", "domain_sha256", "trace_manifest_sha256",
    "trace_file_sha256", "source_git_commit", "worker_binary_sha256",
    "candidate_roster", "arm_roster", "campaigns",
    "logical_freeze_manifest_sha256", "logical_result_stream_sha256",
    "seed_basis", "seed_schedule_sha256", "validator_summary",
    "validator_summary_sha256", "combination_sha256",
))
CAMPAIGN_BINDING_FIELDS = frozenset((
    "candidate_id", "candidate_arm", "run_summary_sha256",
    "freeze_manifest_sha256", "result_stream_sha256",
    "execution_receipt_sha256",
))
ZERO_SHA256 = "0" * 64
RECOVERY_RECORDS = 1440
LOGICAL_RECORDS = 2160
RECOVERY_WORKER_COUNT = 8
MAX_COMPLETED_ARTIFACT_BYTES = 64 * 1024 * 1024


class RecoveryRunnerError(RuntimeError):
    """The recovery-only campaign or logical combination is invalid."""


def fail(message: str) -> None:
    raise RecoveryRunnerError(message)


def _candidate(candidate_id: Any) -> Tuple[str, str]:
    if not isinstance(candidate_id, str) or candidate_id not in CANDIDATE_BY_ID:
        fail("candidate ID is outside the closed recovery-only roster")
    return CANDIDATE_BY_ID[candidate_id]


def _expected_arms(candidate_id: str) \
        -> Tuple[Tuple[str, str], ...]:
    candidate_arm, _ = _candidate(candidate_id)
    return CONTROL_ARMS + (RAW_CONTROL_ARM,) + (
        (candidate_arm, "wirehair2_experiment"),)


def _is_sha256(value: Any) -> bool:
    return isinstance(value, str) and \
        contract_api.SHA256.fullmatch(value) is not None


def _hash_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _hash_jsonl(values: Sequence[Mapping[str, Any]]) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(contract_api.canonical_json(value).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def _self_hash(value: Mapping[str, Any], field: str) -> str:
    return contract_api.sha256_json({
        key: item for key, item in value.items() if key != field
    })


def describe_candidate_worker(
        worker_path: Path, candidate_id: str, deadline: float,
        ) -> Mapping[str, Any]:
    """Describe one exact four-arm worker configuration and bind its file."""
    candidate_arm, expected_descriptor = _candidate(candidate_id)
    try:
        resolved = worker_path.resolve(strict=True)
        info = resolved.stat()
    except OSError as exc:
        fail("cannot resolve native worker {}: {}".format(worker_path, exc))
    if not stat.S_ISREG(info.st_mode) or not os.access(str(resolved), os.X_OK):
        fail("native worker must be an executable regular file")
    binary_sha256 = runner_api._hash_file(resolved)
    description = runner_api._parse_canonical_line(
        runner_api._run_command(
            [str(resolved), "--describe-recovery-candidate", candidate_id],
            deadline, "recovery candidate worker description"),
        "recovery candidate worker description")
    if (set(description) != runner_api.DESCRIPTION_FIELDS or
            description.get("schema") != runner_api.DESCRIPTION_SCHEMA or
            not isinstance(description.get("source_git_commit"), str) or
            re.fullmatch(r"[0-9a-f]{40}",
                         description["source_git_commit"]) is None or
            description.get("binary_sha256") != binary_sha256):
        fail("recovery candidate description does not bind its executable")
    arms = description.get("arms")
    expected = _expected_arms(candidate_id)
    if not isinstance(arms, list) or len(arms) != len(expected):
        fail("recovery candidate description is not an exact four-arm roster")
    for index, expected_arm in enumerate(expected):
        arm = arms[index]
        if (not isinstance(arm, dict) or
                set(arm) != runner_api.DESCRIPTION_ARM_FIELDS or
                arm.get("arm") != expected_arm[0] or
                arm.get("codec") != expected_arm[1] or
                not _is_sha256(arm.get("arm_descriptor_sha256"))):
            fail("recovery candidate description has an invalid arm at index {}"
                 .format(index))
    if tuple(arms[index]["arm_descriptor_sha256"] for index in (0, 1)) != \
            CONTROL_DESCRIPTOR_SHA256S:
        fail("recovery candidate description substitutes a control descriptor")
    if (arms[2]["arm"] != RAW_CONTROL_ARM[0] or
            arms[2]["arm_descriptor_sha256"] !=
            RAW_CONTROL_DESCRIPTOR_SHA256):
        fail("recovery candidate description substitutes its raw control")
    if arms[3]["arm"] != candidate_arm or \
            arms[3]["arm_descriptor_sha256"] != expected_descriptor:
        fail("candidate descriptor does not match its closed candidate ID")
    if len({arm["arm_descriptor_sha256"] for arm in arms}) != len(arms):
        fail("recovery candidate description reuses an arm descriptor")
    result = dict(description)
    result["resolved_path"] = str(resolved)
    result["candidate_id"] = candidate_id
    return result


def _candidate_commands(
        description: Mapping[str, Any], candidate_id: str,
        cpus: Sequence[int]) -> List[List[str]]:
    _candidate(candidate_id)
    worker = description.get("resolved_path")
    if not isinstance(worker, str) or not worker:
        fail("worker description lacks its resolved executable path")
    return [
        [worker, "--describe-recovery-candidate", candidate_id],
        [worker, "--emit-traces", "recovery"],
    ] + [
        [worker, "--recovery-candidate-worker", candidate_id, str(cpu)]
        for cpu in cpus
    ]


def _candidate_freeze(
        contract: Mapping[str, Any], description: Mapping[str, Any],
        candidate_id: str, cpus: Sequence[int], controller_cpu: int,
        source_commit: str, trace_sha256: str) -> Mapping[str, Any]:
    expected = _expected_arms(candidate_id)
    if len(cpus) != RECOVERY_WORKER_COUNT or \
            list(cpus) != sorted(set(cpus)):
        fail("recovery-only screen requires eight sorted unique worker CPUs")
    roster = [arm for arm, _ in expected]
    arms = []
    for index, value in enumerate(description["arms"]):
        if value["arm"] != roster[index]:
            fail("worker description changed before recovery freeze")
        arms.append({
            "arm": value["arm"],
            "codec": value["codec"],
            "binary_sha256": description["binary_sha256"],
            "arm_descriptor_sha256": value["arm_descriptor_sha256"],
            "construction_policy": "not_applicable"
                if value["codec"] == "wirehair1" else "raw_base",
            "repair_map_sha256": ZERO_SHA256,
        })
    return {
        "schema": contract_api.FREEZE_SCHEMA,
        "contract_sha256": contract_api.contract_sha256(contract),
        "evidence_kind": "recovery",
        "phase": "development",
        "domain_sha256": contract["recovery"]["domains"]["development"][
            "domain_sha256"],
        "source_git_commit": source_commit,
        "arm_roster": roster,
        "arm_roster_sha256": contract_api.arm_roster_sha256(roster),
        "trace_manifest_sha256": trace_sha256,
        "repair_training_trace_manifest_sha256": ZERO_SHA256,
        "commands": _candidate_commands(description, candidate_id, cpus),
        "cpu_affinity": list(cpus),
        "host_identity": runner_api._host_identity(controller_cpu),
        "arms": arms,
    }


def _validate_candidate_freeze(
        freeze: Mapping[str, Any], candidate_id: str,
        worker_path: Optional[str] = None) -> None:
    candidate_arm, expected_descriptor = _candidate(candidate_id)
    expected_roster = [arm for arm, _ in _expected_arms(candidate_id)]
    if freeze.get("arm_roster") != expected_roster or \
            freeze.get("evidence_kind") != "recovery" or \
            freeze.get("phase") != "development":
        fail("candidate freeze does not bind the exact four-arm recovery roster")
    arms = freeze.get("arms")
    if (not isinstance(arms, list) or len(arms) != 4 or
            any(not isinstance(arm, dict) for arm in arms) or
            tuple((arm.get("arm"), arm.get("codec")) for arm in arms) !=
            _expected_arms(candidate_id) or
            arms[2].get("arm") != RAW_CONTROL_ARM[0] or
            arms[2].get("arm_descriptor_sha256") !=
            RAW_CONTROL_DESCRIPTOR_SHA256 or
            arms[3].get("arm") != candidate_arm or
            arms[3].get("arm_descriptor_sha256") != expected_descriptor):
        fail("candidate freeze substitutes its candidate descriptor")
    if tuple(arms[index].get("arm_descriptor_sha256")
             for index in (0, 1)) != CONTROL_DESCRIPTOR_SHA256S:
        fail("candidate freeze substitutes a control descriptor")
    binaries = {arm.get("binary_sha256") for arm in arms}
    if len(binaries) != 1 or not _is_sha256(next(iter(binaries))):
        fail("candidate freeze does not use one bound worker binary")
    commands = freeze.get("commands")
    cpus = freeze.get("cpu_affinity")
    if (not isinstance(commands, list) or not isinstance(cpus, list) or
            len(cpus) != RECOVERY_WORKER_COUNT):
        fail("candidate freeze lacks its exact eight-worker argv roster")
    if worker_path is None:
        try:
            worker_path = commands[0][0]
        except (IndexError, TypeError):
            fail("candidate freeze lacks its description argv")
    description_stub = {"resolved_path": worker_path}
    if commands != _candidate_commands(description_stub, candidate_id, cpus):
        fail("candidate freeze commands differ from the exact candidate argv")


def write_recovery_freeze(
        contract: Mapping[str, Any], description: Mapping[str, Any],
        candidate_id: str, cpus: Sequence[int], controller_cpu: int,
        source_commit: str, trace_sha256: str, output_dir: Path,
        ) -> Mapping[str, Any]:
    path = output_dir / "recovery-freeze.json"
    runner_api._atomic_write_object(path, _candidate_freeze(
        contract, description, candidate_id, cpus, controller_cpu,
        source_commit, trace_sha256))
    try:
        loaded = contract_api.load_freeze_manifest(
            contract, "development", path, "recovery")
    except contract_api.ContractError as exc:
        fail(str(exc))
    _validate_candidate_freeze(
        loaded, candidate_id, description["resolved_path"])
    return loaded


def spawn_candidate_workers(
        description: Mapping[str, Any], candidate_id: str,
        cpus: Sequence[int], deadline: float,
        ) -> List[runner_api.PersistentWorker]:
    """Spawn workers in the recovery-only mode; that mode rejects every T job."""
    _candidate(candidate_id)
    workers: List[runner_api.PersistentWorker] = []
    try:
        for cpu in cpus:
            argv = [description["resolved_path"],
                    "--recovery-candidate-worker", candidate_id, str(cpu)]
            try:
                process = subprocess.Popen(
                    argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE, bufsize=0)
            except OSError as exc:
                fail("cannot spawn recovery worker on CPU {}: {}".format(
                    cpu, exc))
            workers.append(runner_api.PersistentWorker(cpu, process, 0))
        pending = set(range(len(workers)))
        ready_deadline = min(deadline, time.monotonic() + 10.0)
        while pending:
            if time.monotonic() >= ready_deadline:
                fail("recovery workers did not establish singleton affinity")
            for index in list(pending):
                worker = workers[index]
                if worker.process.poll() is not None:
                    fail("recovery worker on CPU {} exited during startup: {}"
                         .format(worker.cpu,
                                 runner_api._worker_stderr(worker)))
                try:
                    ticks = native_api._parse_proc_start_ticks(
                        worker.process.pid)
                    affinity = os.sched_getaffinity(worker.process.pid)
                except (native_api.NativeEvidenceError,
                        AttributeError, OSError):
                    continue
                if affinity == {worker.cpu}:
                    worker.start_ticks = ticks
                    os.set_blocking(worker.process.stdout.fileno(), False)
                    os.set_blocking(worker.process.stderr.fileno(), False)
                    pending.remove(index)
            if pending:
                time.sleep(0.01)
        return workers
    except BaseException:
        runner_api.terminate_workers(workers)
        raise


def _candidate_recovery_jobs(
        contract: Mapping[str, Any]) -> List[runner_api.Job]:
    """Build the exact four-arm job ledger without weakening legacy gates."""
    cells = contract["recovery"]["domains"]["development"][
        "expected_cells_per_arm"]
    if type(cells) is not int or cells != 360:
        fail("development recovery domain is not exactly 360 frozen cells")
    jobs = [
        runner_api.Job("recovery", cell, arm, cell * 4 + arm)
        for cell in range(cells) for arm in range(4)
    ]
    if len(jobs) != RECOVERY_RECORDS or any(
            job.ordinal != index for index, job in enumerate(jobs)):
        fail("candidate recovery jobs differ from the four-arm frozen roster")
    return jobs


def _run_recovery_jobs(
        contract: Mapping[str, Any], freeze: Mapping[str, Any],
        description: Mapping[str, Any],
        workers: Sequence[runner_api.PersistentWorker], output_dir: Path,
        window_start_ns: int, deadline: float) -> Tuple[Path, int]:
    if len(workers) != RECOVERY_WORKER_COUNT:
        fail("recovery-only screen requires exactly eight workers")
    jobs = _candidate_recovery_jobs(contract)
    if len(jobs) != RECOVERY_RECORDS or \
            any(not job.command().startswith(b"R ") for job in jobs):
        fail("recovery-only job roster is not exactly 1440 R commands")
    path = output_dir / "recovery-native-results.jsonl"
    sink = runner_api.AtomicLineSink(path)
    try:
        maximum_end, used = runner_api.run_job_batch(
            workers, jobs, 0, sink, deadline,
            runner_api._strict_response_validator(
                contract, freeze, "recovery", description, window_start_ns))
        if used != {worker.cpu for worker in workers}:
            fail("recovery campaign did not exercise every frozen CPU")
        sink.publish()
        return path, maximum_end
    finally:
        sink.abort()


def run_recovery_screen(args: argparse.Namespace) -> Mapping[str, Any]:
    candidate_id = args.candidate
    candidate_arm, _ = _candidate(candidate_id)
    if (not math.isfinite(args.deadline_seconds) or
            not 0.0 < args.deadline_seconds <= runner_api.MAX_WALL_SECONDS):
        fail("--deadline-seconds must be in (0,7200]")
    deadline = time.monotonic() + args.deadline_seconds
    contract = contract_api.load_contract(args.contract)
    description = describe_candidate_worker(args.worker, candidate_id, deadline)
    try:
        original_affinity = set(os.sched_getaffinity(0))
    except (AttributeError, OSError) as exc:
        fail("cannot inspect initial controller affinity: {}".format(exc))
    explicit = runner_api.parse_cpu_list(args.cpus) \
        if args.cpus is not None else None
    cpus, controller_cpu = runner_api.select_cpu_layout(
        args.sampler_cpu, explicit, args.controller_cpu,
        affinity=original_affinity)
    runner_api._preflight_sampler(
        args.sampler_pid, args.sampler_cpu,
        args.sampler_script, args.sampler_csv)
    source_commit = runner_api._git_head(deadline)
    runner_api._require_worker_source_commit(description, source_commit)
    output_dir = runner_api._create_output_dir(args.output_dir)

    trace_path, trace_sha256 = runner_api._emit_and_assemble_trace(
        contract, "recovery", description, output_dir, deadline)
    if runner_api._git_head(deadline) != source_commit:
        fail("codec source commit changed before the recovery freeze")
    freeze = write_recovery_freeze(
        contract, description, candidate_id, cpus, controller_cpu,
        source_commit, trace_sha256, output_dir)
    freeze_path = output_dir / "recovery-freeze.json"

    workers: List[runner_api.PersistentWorker] = []
    clean_shutdown = False
    controller_pinned = False
    completed_summary: Optional[Mapping[str, Any]] = None
    try:
        workers = spawn_candidate_workers(
            description, candidate_id, cpus, deadline)
        controller_pinned = True
        runner_api._pin_controller(controller_cpu)
        window_start_ns = runner_api.choose_new_sampler_start(
            args.sampler_csv, deadline)
        native_path, maximum_worker_end = _run_recovery_jobs(
            contract, freeze, description, workers, output_dir,
            window_start_ns, deadline)
        window_end_ns = runner_api._wait_for_sampler_sample(
            args.sampler_csv, deadline,
            at_or_after_ns=maximum_worker_end)
        if runner_api._git_head(deadline) != source_commit:
            fail("codec source commit changed during the recovery screen")
        sampler_path = output_dir / "sampler-attestation.json"
        try:
            native_api.write_sampler_attestation(
                args.sampler_pid, args.sampler_cpu,
                args.sampler_script, args.sampler_csv,
                window_start_ns, window_end_ns, sampler_path)
        except native_api.NativeEvidenceError as exc:
            fail(str(exc))

        result_path = output_dir / "recovery-results.jsonl"
        receipt_path = output_dir / "recovery-execution.json"
        try:
            assembled = native_api.assemble_results(
                contract, "recovery", "development", freeze_path,
                trace_path, native_path, sampler_path,
                result_path, receipt_path, verify_live_sampler=True)
            runner_api._remaining(deadline, "assembling recovery results")
            validated = native_api.validate_execution_receipt(
                contract, "recovery", "development", freeze_path,
                trace_path, native_path, result_path, receipt_path,
                verify_live_sampler=True)
            runner_api._remaining(deadline, "validating recovery execution")
        except (native_api.NativeEvidenceError,
                contract_api.ContractError) as exc:
            fail(str(exc))
        if assembled != validated:
            fail("recovery assembly and terminal validation disagree")
        runner_api.quit_workers(workers, deadline)
        clean_shutdown = True

        receipt = assembled["execution_receipt"]
        if receipt.get("record_count") != RECOVERY_RECORDS:
            fail("recovery execution receipt does not contain 1440 records")
        thermal = receipt["thermal"]
        summary: Dict[str, Any] = {
            "schema": CAMPAIGN_SUMMARY_SCHEMA,
            "status": "complete",
            "output_dir": str(output_dir),
            "candidate_id": candidate_id,
            "candidate_arm": candidate_arm,
            "source_git_commit": source_commit,
            "contract_sha256": contract_api.contract_sha256(contract),
            "domain_sha256": contract["recovery"]["domains"]
                ["development"]["domain_sha256"],
            "trace_manifest_sha256": trace_sha256,
            "worker_binary_sha256": description["binary_sha256"],
            "controller_cpu": controller_cpu,
            "worker_cpus": list(cpus),
            "recovery_records": receipt["record_count"],
            "recovery_freeze_sha256": receipt["freeze_manifest_sha256"],
            "recovery_result_sha256": receipt["result_stream_sha256"],
            "recovery_execution_receipt_sha256": receipt["receipt_sha256"],
            "thermal_samples": thermal["sample_count"],
            "cpu_tctl_max_millic": thermal["cpu_tctl_max_millic"],
            "dimm_max_millic": thermal["dimm_max_millic"],
            "seed_basis": RAW_SEED_BASIS,
            "seed_schedule_sha256": RAW_SEED_SCHEDULE_SHA256,
        }
        summary["summary_sha256"] = contract_api.sha256_json(summary)
        completed_summary = summary
    finally:
        if workers and not clean_shutdown:
            runner_api.terminate_workers(workers)
        if controller_pinned:
            runner_api._restore_controller_affinity(original_affinity)
    if completed_summary is None:
        fail("recovery screen reached cleanup without a completed summary")
    runner_api._atomic_write_object(
        output_dir / "run-summary.json", completed_summary)
    return completed_summary


def _load_canonical_object(path: Path, context: str) -> Mapping[str, Any]:
    return _parse_canonical_object(_read_regular_bytes(path, context), context)


def _parse_canonical_object(data: bytes, context: str) -> Mapping[str, Any]:
    try:
        value = contract_api._load_json_bytes(data, context)
    except contract_api.ContractError as exc:
        fail(str(exc))
    if not isinstance(value, dict):
        fail("{} must be a JSON object".format(context))
    logical = contract_api.canonical_json(value).encode("utf-8")
    if data != logical + b"\n":
        fail("{} must be exact canonical JSON followed by LF".format(context))
    return value


def _read_regular_bytes(path: Path, context: str) -> bytes:
    """Read one bounded, already-opened regular artifact without symlink races."""
    descriptor = -1
    try:
        nofollow = getattr(os, "O_NOFOLLOW", 0)
        if nofollow == 0:
            fail("{} cannot be opened fail-closed without O_NOFOLLOW".format(
                context))
        flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | nofollow | \
            getattr(os, "O_NONBLOCK", 0)
        descriptor = os.open(str(path), flags)
        info = os.fstat(descriptor)
        if not stat.S_ISREG(info.st_mode):
            fail("{} must be a regular non-symlink file".format(context))
        if info.st_size > MAX_COMPLETED_ARTIFACT_BYTES:
            fail("{} exceeds the bounded artifact size".format(context))
        chunks = []
        size = 0
        while True:
            block = os.read(descriptor, 1024 * 1024)
            if not block:
                break
            size += len(block)
            if size > MAX_COMPLETED_ARTIFACT_BYTES:
                fail("{} exceeds the bounded artifact size".format(context))
            chunks.append(block)
        return b"".join(chunks)
    except OSError as exc:
        fail("cannot read {} {}: {}".format(context, path, exc))
    finally:
        if descriptor >= 0:
            os.close(descriptor)
    return b""


def _parse_exact_jsonl(
        data: bytes, context: str) -> List[Mapping[str, Any]]:
    if not data or not data.endswith(b"\n"):
        fail("{} must be a nonempty newline-terminated JSONL stream".format(
            context))
    rows = []
    for index, line in enumerate(data.splitlines(keepends=True), 1):
        try:
            rows.append(runner_api._parse_canonical_line(
                line, "{} line {}".format(context, index)))
        except runner_api.RunnerError as exc:
            fail(str(exc))
    return rows


def _validate_bound_recovery_trace_bytes(
        contract: Mapping[str, Any], freeze: Mapping[str, Any], data: bytes,
        ) -> None:
    rows = _parse_exact_jsonl(data, "campaign trace manifest")
    cells = list(contract_api.iter_recovery_cells(contract, "development"))
    expected_count = contract["recovery"]["domains"]["development"][
        "expected_cells_per_arm"]
    if len(cells) != expected_count or len(rows) != expected_count:
        fail("campaign trace manifest has the wrong frozen cardinality")
    digest = contract_api._trace_manifest_hasher(
        contract, "recovery", "development")
    for ordinal, (row, cell) in enumerate(zip(rows, cells)):
        if (set(row) != contract_api.TRACE_FIELDS or
                type(row.get("ordinal")) is not int or
                row.get("ordinal") != ordinal or
                not _is_sha256(row.get("cell_sha256")) or
                not _is_sha256(row.get("trace_sha256")) or
                row.get("cell_sha256") != contract_api.sha256_json(cell)):
            fail("campaign trace manifest row {} differs from its frozen cell"
                 .format(ordinal))
        contract_api._hash_trace_manifest_row(digest, row)
    if digest.hexdigest() != freeze.get("trace_manifest_sha256"):
        fail("campaign trace bytes differ from the frozen trace hash")


def _write_snapshot(path: Path, data: bytes) -> None:
    try:
        with path.open("xb") as output:
            output.write(data)
            output.flush()
            os.fsync(output.fileno())
    except OSError as exc:
        fail("cannot write private campaign snapshot {}: {}".format(path, exc))


def load_completed_campaign(
        contract: Mapping[str, Any], campaign_dir: Path,
        ) -> Mapping[str, Any]:
    """Revalidate one terminal four-arm campaign without claiming liveness."""
    try:
        info = os.lstat(str(campaign_dir))
        if not stat.S_ISDIR(info.st_mode):
            fail("campaign path must be a real directory, not a symlink")
        resolved = campaign_dir.resolve(strict=True)
    except OSError as exc:
        fail("cannot open campaign directory {}: {}".format(campaign_dir, exc))
    source_paths = {
        "summary": resolved / "run-summary.json",
        "freeze": resolved / "recovery-freeze.json",
        "trace": resolved / "recovery-traces.jsonl",
        "native": resolved / "recovery-native-results.jsonl",
        "result": resolved / "recovery-results.jsonl",
        "receipt": resolved / "recovery-execution.json",
    }
    contexts = {
        "summary": "campaign run summary",
        "freeze": "campaign recovery freeze",
        "trace": "campaign trace manifest",
        "native": "campaign native recovery stream",
        "result": "campaign recovery ledger",
        "receipt": "campaign recovery execution receipt",
    }
    snapshots = {
        key: _read_regular_bytes(path, contexts[key])
        for key, path in source_paths.items()
    }
    with tempfile.TemporaryDirectory(prefix="wh2-recovery-snapshot-") as raw:
        snapshot_root = Path(raw)
        paths = {
            key: snapshot_root / path.name
            for key, path in source_paths.items()
        }
        for key, path in paths.items():
            _write_snapshot(path, snapshots[key])

        summary = _parse_canonical_object(
            snapshots["summary"], "campaign run summary")
        if set(summary) != CAMPAIGN_SUMMARY_FIELDS or \
                summary.get("schema") != CAMPAIGN_SUMMARY_SCHEMA or \
                summary.get("status") != "complete" or \
                summary.get("summary_sha256") != _self_hash(
                    summary, "summary_sha256"):
            fail("campaign run summary is incomplete or has an invalid self-hash")
        candidate_id = summary.get("candidate_id")
        candidate_arm, _ = _candidate(candidate_id)
        if summary.get("candidate_arm") != candidate_arm or \
                summary.get("output_dir") != str(resolved):
            fail("campaign run summary is bound to another candidate or directory")
        try:
            validated = native_api.validate_execution_receipt(
                contract, "recovery", "development", paths["freeze"],
                paths["trace"], paths["native"], paths["result"],
                paths["receipt"], verify_live_sampler=False)
            freeze = contract_api.load_freeze_manifest(
                contract, "development", paths["freeze"], "recovery")
        except (native_api.NativeEvidenceError,
                contract_api.ContractError) as exc:
            fail(str(exc))
        _validate_candidate_freeze(freeze, candidate_id)
        receipt = validated["execution_receipt"]
        if contract_api.freeze_manifest_sha256(freeze) != \
                receipt.get("freeze_manifest_sha256"):
            fail("reopened candidate freeze differs from its execution receipt")
        expected_summary = {
            "source_git_commit": freeze["source_git_commit"],
            "contract_sha256": contract_api.contract_sha256(contract),
            "domain_sha256": freeze["domain_sha256"],
            "trace_manifest_sha256": freeze["trace_manifest_sha256"],
            "worker_binary_sha256": freeze["arms"][0]["binary_sha256"],
            "recovery_records": receipt["record_count"],
            "recovery_freeze_sha256": receipt["freeze_manifest_sha256"],
            "recovery_result_sha256": receipt["result_stream_sha256"],
            "recovery_execution_receipt_sha256": receipt["receipt_sha256"],
            "worker_cpus": receipt["worker_cpus"],
            "controller_cpu": freeze.get("host_identity", {}).get(
                "controller_cpu"),
            "thermal_samples": receipt["thermal"]["sample_count"],
            "cpu_tctl_max_millic": receipt["thermal"]["cpu_tctl_max_millic"],
            "dimm_max_millic": receipt["thermal"]["dimm_max_millic"],
            "seed_basis": RAW_SEED_BASIS,
            "seed_schedule_sha256": RAW_SEED_SCHEDULE_SHA256,
        }
        if (type(expected_summary["controller_cpu"]) is not int or
                expected_summary["controller_cpu"] < 0 or
                any(summary.get(field) != value
                    for field, value in expected_summary.items()) or
                summary.get("recovery_records") != RECOVERY_RECORDS):
            fail("campaign run summary differs from its validated execution")
        rows = _parse_exact_jsonl(
            snapshots["result"], "campaign recovery ledger")
        if len(rows) != RECOVERY_RECORDS or \
                _hash_bytes(snapshots["result"]) != \
                receipt["result_stream_sha256"]:
            fail("campaign recovery ledger differs from its receipted result hash")
        _validate_bound_recovery_trace_bytes(
            contract, freeze, snapshots["trace"])
        return {
            "directory": str(resolved),
            "directory_identity": (info.st_dev, info.st_ino),
            "candidate_id": candidate_id,
            "candidate_arm": candidate_arm,
            "summary": summary,
            "freeze": freeze,
            "receipt": receipt,
            "rows": rows,
            "trace_bytes": snapshots["trace"],
        }


def _combine_loaded_campaigns(
        contract: Mapping[str, Any], campaigns: Sequence[Mapping[str, Any]],
        ) -> Tuple[Mapping[str, Any], List[Mapping[str, Any]],
                   List[Mapping[str, Any]], bytes]:
    """Construct a logical six-arm freeze/ledger from validated payloads."""
    if len(campaigns) != len(CANDIDATE_SPECS):
        fail("combination requires exactly three completed campaign directories")
    by_id: Dict[str, Mapping[str, Any]] = {}
    identities = set()
    for campaign in campaigns:
        candidate_id = campaign.get("candidate_id")
        _candidate(candidate_id)
        if candidate_id in by_id:
            fail("combination contains a duplicate candidate campaign")
        identity = campaign.get("directory_identity")
        if identity is not None and identity in identities:
            fail("combination reuses the same campaign directory")
        identities.add(identity)
        by_id[candidate_id] = campaign
    expected_ids = [value[0] for value in CANDIDATE_SPECS]
    if set(by_id) != set(expected_ids):
        fail("combination does not cover the closed candidate roster")
    ordered = [by_id[candidate_id] for candidate_id in expected_ids]
    first = ordered[0]

    common_fields = (
        "contract_sha256", "domain_sha256", "trace_manifest_sha256",
        "source_git_commit",
    )
    for campaign in ordered[1:]:
        for field in common_fields:
            if campaign["freeze"].get(field) != first["freeze"].get(field):
                fail("campaigns differ in their {}".format(field))
        if campaign.get("trace_bytes") != first.get("trace_bytes"):
            fail("campaign trace manifests are not byte-identical")
    expected_contract = contract_api.contract_sha256(contract)
    expected_domain = contract["recovery"]["domains"]["development"][
        "domain_sha256"]
    if first["freeze"].get("contract_sha256") != expected_contract or \
            first["freeze"].get("domain_sha256") != expected_domain:
        fail("campaigns do not bind the loaded recovery contract/domain")

    binaries = set()
    for campaign in ordered:
        if contract_api.freeze_manifest_sha256(campaign["freeze"]) != \
                campaign["receipt"].get("freeze_manifest_sha256"):
            fail("candidate campaign freeze differs from its receipt hash")
        freeze_arms = campaign["freeze"].get("arms")
        if not isinstance(freeze_arms, list):
            fail("campaign freeze arm records are malformed")
        binaries.update(arm.get("binary_sha256") for arm in freeze_arms)
    if len(binaries) != 1 or not _is_sha256(next(iter(binaries))):
        fail("campaigns were not produced by one identical worker binary")
    first_common = first["freeze"]["arms"][:3]
    for campaign in ordered[1:]:
        if campaign["freeze"]["arms"][:3] != first_common:
            fail("campaigns substitute a common control arm descriptor")
    cell_count = contract["recovery"]["domains"]["development"][
        "expected_cells_per_arm"]
    if cell_count * 4 != RECOVERY_RECORDS:
        fail("development recovery cardinality is not 360 cells per arm")
    for campaign in ordered:
        rows = campaign.get("rows")
        if not isinstance(rows, list) or len(rows) != RECOVERY_RECORDS:
            fail("candidate campaign does not contain exactly 1440 payload rows")
        if _hash_jsonl(rows) != campaign["receipt"].get(
                "result_stream_sha256"):
            fail("candidate campaign rows differ from their result-stream hash")
        expected_roster = campaign["freeze"]["arm_roster"]
        for cell in range(cell_count):
            for arm_index, arm in enumerate(expected_roster):
                if rows[cell * 4 + arm_index].get("arm") != arm:
                    fail("candidate ledger rows are not in frozen ordinal order")
    for campaign in ordered[1:]:
        for cell in range(cell_count):
            for common_index in (0, 1, 2):
                if campaign["rows"][cell * 4 + common_index] != \
                        first["rows"][cell * 4 + common_index]:
                    fail("campaign common payloads are not cell-identical")

    roster = [value[0] for value in CONTROL_ARMS] + [RAW_CONTROL_ARM[0]] + [
        value[1] for value in CANDIDATE_SPECS
    ]
    logical_arms = [dict(value) for value in first_common] + [
        dict(campaign["freeze"]["arms"][3]) for campaign in ordered
    ]
    commands: List[List[str]] = []
    cpu_affinity: Set[int] = set()
    input_freezes = []
    logical_rows: List[Mapping[str, Any]] = []
    for campaign in ordered:
        commands.extend(campaign["freeze"]["commands"])
        cpu_affinity.update(campaign["freeze"]["cpu_affinity"])
        input_freezes.append(contract_api.freeze_manifest_sha256(
            campaign["freeze"]))
    for cell in range(cell_count):
        logical_rows.extend((first["rows"][cell * 4],
                             first["rows"][cell * 4 + 1],
                             first["rows"][cell * 4 + 2]))
        logical_rows.extend(
            campaign["rows"][cell * 4 + 3] for campaign in ordered)
    if len(logical_rows) != LOGICAL_RECORDS:
        fail("logical recovery ledger is not exactly six arms by 360 cells")

    logical_freeze = {
        "schema": contract_api.FREEZE_SCHEMA,
        "contract_sha256": expected_contract,
        "evidence_kind": "recovery",
        "phase": "development",
        "domain_sha256": expected_domain,
        "source_git_commit": first["freeze"]["source_git_commit"],
        "arm_roster": roster,
        "arm_roster_sha256": contract_api.arm_roster_sha256(roster),
        "trace_manifest_sha256": first["freeze"]["trace_manifest_sha256"],
        "repair_training_trace_manifest_sha256": ZERO_SHA256,
        "commands": commands,
        "cpu_affinity": sorted(cpu_affinity),
        "host_identity": {
            "artifact_kind": "logical_recovery_combination",
            "input_freeze_manifest_sha256s": input_freezes,
        },
        "arms": logical_arms,
    }
    try:
        # Validate the synthesized structure before allowing publication.
        contract_api._exact_keys(
            logical_freeze, contract_api.FREEZE_FIELDS,
            "logical recovery freeze")
    except contract_api.ContractError as exc:
        fail(str(exc))
    bindings = [{
        "candidate_id": campaign["candidate_id"],
        "candidate_arm": campaign["candidate_arm"],
        "run_summary_sha256": campaign["summary"]["summary_sha256"],
        "freeze_manifest_sha256": campaign["receipt"]
            ["freeze_manifest_sha256"],
        "result_stream_sha256": campaign["receipt"]["result_stream_sha256"],
        "execution_receipt_sha256": campaign["receipt"]["receipt_sha256"],
    } for campaign in ordered]
    return logical_freeze, logical_rows, bindings, first["trace_bytes"]


def combine_recovery_screens(args: argparse.Namespace) -> Mapping[str, Any]:
    contract = contract_api.load_contract(args.contract)
    campaigns = [load_completed_campaign(contract, path)
                 for path in args.campaign_dir]
    logical_freeze, logical_rows, bindings, trace_bytes = \
        _combine_loaded_campaigns(contract, campaigns)
    output_dir = runner_api._create_output_dir(args.output_dir)
    trace_path = output_dir / "logical-recovery-traces.jsonl"
    freeze_path = output_dir / "logical-recovery-freeze.json"
    result_path = output_dir / "logical-recovery-results.jsonl"
    runner_api._atomic_write_bytes(trace_path, trace_bytes)
    runner_api._atomic_write_object(freeze_path, logical_freeze)
    sink = runner_api.AtomicLineSink(result_path)
    try:
        for row in logical_rows:
            sink.write((contract_api.canonical_json(row) + "\n").encode(
                "utf-8"))
        sink.publish()
    finally:
        sink.abort()
    try:
        loaded_freeze = contract_api.load_freeze_manifest(
            contract, "development", freeze_path, "recovery")
        validator_summary = contract_api.validate_ledger(
            contract, "development", result_path, freeze_path, trace_path)
    except contract_api.ContractError as exc:
        fail(str(exc))
    if contract_api.freeze_manifest_sha256(loaded_freeze) != \
            validator_summary.get("freeze_manifest_sha256"):
        fail("logical validator summary does not bind its synthesized freeze")
    result_hash = _hash_jsonl(logical_rows)
    try:
        if runner_api._hash_file(result_path) != result_hash:
            fail("published logical ledger differs from its canonical hash")
    except runner_api.RunnerError as exc:
        fail(str(exc))
    summary: Dict[str, Any] = {
        "schema": COMBINATION_SUMMARY_SCHEMA,
        "status": "complete",
        "artifact_kind": "logical_recovery_combination",
        "is_execution_receipt": False,
        "phase": "development",
        "contract_sha256": contract_api.contract_sha256(contract),
        "domain_sha256": logical_freeze["domain_sha256"],
        "trace_manifest_sha256": logical_freeze["trace_manifest_sha256"],
        "trace_file_sha256": _hash_bytes(trace_bytes),
        "source_git_commit": logical_freeze["source_git_commit"],
        "worker_binary_sha256": logical_freeze["arms"][0]["binary_sha256"],
        "candidate_roster": [value[0] for value in CANDIDATE_SPECS],
        "arm_roster": logical_freeze["arm_roster"],
        "campaigns": bindings,
        "seed_basis": RAW_SEED_BASIS,
        "seed_schedule_sha256": RAW_SEED_SCHEDULE_SHA256,
        "logical_freeze_manifest_sha256":
            contract_api.freeze_manifest_sha256(loaded_freeze),
        "logical_result_stream_sha256": result_hash,
        "validator_summary": validator_summary,
        "validator_summary_sha256":
            contract_api.sha256_json(validator_summary),
    }
    summary["combination_sha256"] = contract_api.sha256_json(summary)
    if set(summary) != COMBINATION_SUMMARY_FIELDS or \
            any(set(value) != CAMPAIGN_BINDING_FIELDS
                for value in summary["campaigns"]):
        fail("internal logical-combination summary schema mismatch")
    runner_api._atomic_write_object(
        output_dir / "logical-recovery-summary.json", summary)
    return summary


def _add_common_run_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--contract", type=Path, default=contract_api.DEFAULT_CONTRACT)
    parser.add_argument(
        "--worker", type=Path,
        default=Path("build/wirehair_wh2_contract_worker"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--candidate", choices=tuple(CANDIDATE_BY_ID),
                        required=True)
    parser.add_argument("--sampler-pid", type=int, required=True)
    parser.add_argument("--sampler-cpu", type=int, required=True)
    parser.add_argument("--sampler-script", type=Path, required=True)
    parser.add_argument("--sampler-csv", type=Path, required=True)
    parser.add_argument(
        "--cpus", help="strictly increasing comma-separated logical CPUs")
    parser.add_argument(
        "--controller-cpu", type=int,
        help="dedicated non-worker, non-sampler physical core")
    parser.add_argument(
        "--deadline-seconds", type=float,
        default=runner_api.MAX_WALL_SECONDS)


def main(argv: Sequence[str] = ()) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser(
        "run", help="run one exact four-arm native recovery campaign")
    _add_common_run_arguments(run)
    combine = commands.add_parser(
        "combine", help="combine three completed runs into a logical ledger")
    combine.add_argument(
        "--contract", type=Path, default=contract_api.DEFAULT_CONTRACT)
    combine.add_argument(
        "--campaign-dir", type=Path, action="append", required=True,
        help="repeat exactly three times, once per closed candidate")
    combine.add_argument("--output-dir", type=Path, required=True)
    combine.add_argument(
        "--deadline-seconds", type=float,
        default=runner_api.MAX_WALL_SECONDS)
    args = parser.parse_args(argv if argv else None)
    try:
        if args.command == "run":
            if (not math.isfinite(args.deadline_seconds) or
                    not 0.0 < args.deadline_seconds <=
                    runner_api.MAX_WALL_SECONDS):
                fail("--deadline-seconds must be in (0,7200]")
            with runner_api._hard_wall(args.deadline_seconds):
                summary = run_recovery_screen(args)
        else:
            if len(args.campaign_dir) != len(CANDIDATE_SPECS):
                fail("--campaign-dir must be supplied exactly three times")
            if (not math.isfinite(args.deadline_seconds) or
                    not 0.0 < args.deadline_seconds <=
                    runner_api.MAX_WALL_SECONDS):
                fail("--deadline-seconds must be in (0,7200]")
            with runner_api._hard_wall(args.deadline_seconds):
                summary = combine_recovery_screens(args)
    except (RecoveryRunnerError, runner_api.RunnerError,
            native_api.NativeEvidenceError, contract_api.ContractError,
            OSError, UnicodeError) as exc:
        print("wh2 native recovery screen: {}".format(exc), file=sys.stderr)
        return 1
    print(contract_api.canonical_json(summary))
    return 0


if __name__ == "__main__":
    sys.exit(main())
