#!/usr/bin/env python3
"""Parallel strict verifier for an immutable WH2 za5v campaign runner.

This wrapper source-loads the manifest-bound runner and dependencies from exact
snapshotted bytes, verifies those runtime bindings before and after replay, and
monkey-patches only ``build_summaries`` in memory.  Every result passes a
source-pinned reproduction of the frozen runner's complete result validation.
Only independent job replay is parallel; aggregation and artifact generation
retain manifest order.
"""

import argparse
from contextlib import contextmanager, nullcontext
import gzip
import hashlib
import multiprocessing
import os
from pathlib import Path
import queue
import signal
import stat
import sys
import time
import types
import zlib


MAX_REPLAY_WORKERS = 32
# These three source hashes are checked before executing any campaign-provided
# Python.  The manifest binding is then checked independently after the exact
# source snapshots are running.
SUPPORTED_RUNNER_SHA256 = frozenset({
    "bae4336508d9d5fd9dd514950da709cb01cc129dcf9e0b36f3e8fe96e17346ef",
})
SUPPORTED_PARSER_SHA256 = frozenset({
    "ceeebe3cab0507fdfe86db9067d2e1aea1741aaeb1629cdcaafb589bb568e519",
})
SUPPORTED_CONTEXT_TOOL_SHA256 = frozenset({
    "07cd2dd509e0847caaaab813ba2d7fe6f4cf57a4fe2d78e2c9400c2545f87c82",
})
_PROJECTION_FIELDS = (
    "arms",
    "contrasts",
    "replicates",
    "stream_sha256",
    "context",
    "evidence",
)
_WORKER_CONTEXT = None
_LOADED_CAMPAIGNS = {}
RESULT_COMPRESSED_FIXED_BYTES = 64 * 1024
RESULT_COMPRESSED_BYTES_PER_REPLICATE = 8 * 1024
RESULT_UNCOMPRESSED_FIXED_BYTES = 512 * 1024
RESULT_UNCOMPRESSED_BYTES_PER_REPLICATE = 64 * 1024
REPLAY_QUEUE_POLL_SECONDS = 0.05
REPLAY_TERMINATE_GRACE_SECONDS = 0.5
REPLAY_KILL_GRACE_SECONDS = 1.0
SOURCE_BYTE_LIMITS = {
    "runner": 512 * 1024,
    "parser": 256 * 1024,
    "context_tool": 512 * 1024,
}


def _source_snapshot(path, kind):
    """Return the exact stable bytes that will be compiled and executed."""
    if kind not in SOURCE_BYTE_LIMITS:
        raise RuntimeError("frozen source kind is invalid")
    label = kind.replace("_", " ")
    limit = SOURCE_BYTE_LIMITS[kind]
    path = Path(path).resolve()
    try:
        with open(path, "rb") as source:
            before = os.fstat(source.fileno())
            if not stat.S_ISREG(before.st_mode):
                raise RuntimeError(
                    f"frozen {label} is not regular: {path}")
            if before.st_size > limit:
                raise RuntimeError(
                    f"frozen {label} exceeds source byte limit: {path}")
            chunks = []
            digest = hashlib.sha256()
            payload_bytes = 0
            while payload_bytes <= limit:
                chunk = source.read(
                    min(64 * 1024, limit - payload_bytes + 1))
                if not chunk:
                    break
                chunks.append(chunk)
                digest.update(chunk)
                payload_bytes += len(chunk)
            if (
                payload_bytes > limit
                or payload_bytes != before.st_size
            ):
                raise RuntimeError(
                    f"frozen {label} changed while being snapshotted")
            payload = b"".join(chunks)
            after = os.fstat(source.fileno())
        current = os.stat(path, follow_symlinks=False)
    except OSError as error:
        raise RuntimeError(
            f"could not snapshot frozen {label} {path}: {error}")
    stable = (
        "st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns",
        "st_ctime_ns",
    )
    if (
        any(
            getattr(before, name) != getattr(after, name)
            for name in stable
        )
        or any(
            getattr(after, name) != getattr(current, name)
            for name in stable
        )
    ):
        raise RuntimeError(
            f"frozen {label} changed while being snapshotted")
    return {
        "kind": kind,
        "path": path,
        "identity": tuple(getattr(after, name) for name in stable),
        "sha256": digest.hexdigest(),
        "payload": payload,
    }


def _snapshot_attestation(snapshot):
    return {
        "kind": snapshot["kind"],
        "path": snapshot["path"],
        "identity": snapshot["identity"],
        "sha256": snapshot["sha256"],
    }


def _execute_source_snapshot(module_name, snapshot):
    """Execute source bytes directly, never a timestamp-selected .pyc."""
    module = types.ModuleType(module_name)
    module.__file__ = str(snapshot["path"])
    module.__package__ = ""
    module.__loader__ = None
    module.__spec__ = None
    sys.modules[module_name] = module
    try:
        code = compile(
            snapshot["payload"],
            str(snapshot["path"]),
            "exec",
            dont_inherit=True,
        )
        exec(code, module.__dict__)
        if _snapshot_attestation(
                _source_snapshot(
                    snapshot["path"], snapshot["kind"])) != \
                _snapshot_attestation(snapshot):
            raise RuntimeError(
                f"frozen {module_name} changed while being executed")
    except BaseException:
        if sys.modules.get(module_name) is module:
            del sys.modules[module_name]
        raise
    return module


def _campaign_attestation(campaign):
    for attestation in _LOADED_CAMPAIGNS.values():
        if attestation["modules"]["runner"] is campaign:
            return attestation
    return None


def _assert_source_attestation(campaign):
    attestation = _campaign_attestation(campaign)
    if attestation is None:
        raise RuntimeError(
            "campaign module was not source-loaded by this verifier")
    expected_names = {
        "runner": attestation["module_name"],
        "parser": "band_timing_codec",
        "context_tool": "peel_codec",
    }
    for name, module_name in expected_names.items():
        module = attestation["modules"][name]
        if sys.modules.get(module_name) is not module:
            raise RuntimeError(
                f"source-attested {name} module was substituted")
        current = _source_snapshot(
            attestation["sources"][name]["path"], name)
        if _snapshot_attestation(current) != \
                attestation["sources"][name]:
            raise RuntimeError(
                f"source-attested {name} changed after import")
    if (
        campaign.band is not attestation["modules"]["parser"]
        or campaign.peel_codec
            is not attestation["modules"]["context_tool"]
    ):
        raise RuntimeError(
            "campaign dependency objects were substituted")
    return attestation


def _load_frozen_campaign(path):
    path = Path(path).resolve()
    repository = path.parents[1]
    snapshots = {
        "runner": _source_snapshot(path, "runner"),
        "parser": _source_snapshot(
            repository / "tools" / "band_timing_codec.py",
            "parser",
        ),
        "context_tool": _source_snapshot(
            repository / "tools" / "peel_codec.py",
            "context_tool",
        ),
    }
    supported = {
        "runner": SUPPORTED_RUNNER_SHA256,
        "parser": SUPPORTED_PARSER_SHA256,
        "context_tool": SUPPORTED_CONTEXT_TOOL_SHA256,
    }
    unsupported = [
        name for name, snapshot in snapshots.items()
        if snapshot["sha256"] not in supported[name]
    ]
    if unsupported:
        raise RuntimeError(
            "parallel wrapper does not support source snapshot(s): "
            + ",".join(unsupported))
    module_name = (
        "_wirehair_frozen_wh2_za5v_"
        + hashlib.sha256(str(path).encode("utf-8")).hexdigest()[:12]
        + "_"
        + snapshots["runner"]["sha256"]
    )
    cached = _LOADED_CAMPAIGNS.get(module_name)
    if cached is not None:
        campaign = cached["modules"]["runner"]
        _assert_source_attestation(campaign)
        return campaign
    occupied = [
        name for name in ("peel_codec", "band_timing_codec", module_name)
        if name in sys.modules
    ]
    if occupied:
        raise RuntimeError(
            "refusing preloaded executable campaign modules: "
            + ",".join(occupied))
    loaded = {}
    try:
        loaded["context_tool"] = _execute_source_snapshot(
            "peel_codec", snapshots["context_tool"])
        loaded["parser"] = _execute_source_snapshot(
            "band_timing_codec", snapshots["parser"])
        loaded["runner"] = _execute_source_snapshot(
            module_name, snapshots["runner"])
    except BaseException:
        for name in (
                module_name, "band_timing_codec", "peel_codec"):
            sys.modules.pop(name, None)
        raise
    attestation = {
        "module_name": module_name,
        "modules": loaded,
        "sources": {
            name: _snapshot_attestation(snapshot)
            for name, snapshot in snapshots.items()
        },
    }
    _LOADED_CAMPAIGNS[module_name] = attestation
    try:
        _assert_source_attestation(loaded["runner"])
    except BaseException:
        _LOADED_CAMPAIGNS.pop(module_name, None)
        for name in (
                module_name, "band_timing_codec", "peel_codec"):
            sys.modules.pop(name, None)
        raise
    return loaded["runner"]


def _assert_imported_runner_binding(campaign, runtime_bindings):
    campaign._validate_runtime_bindings_schema(runtime_bindings)
    attestation = _campaign_attestation(campaign)
    if attestation is not None:
        _assert_source_attestation(campaign)
    imported = {
        "runner": campaign,
        "parser": campaign.band,
        "context_tool": campaign.peel_codec,
    }
    for name, module in imported.items():
        expected = Path(runtime_bindings[name]["path"]).resolve()
        actual = Path(module.__file__).resolve()
        if actual != expected:
            raise campaign.CampaignError(
                f"imported {name} is not manifest-bound")
        if (
            attestation is not None
            and (
                attestation["sources"][name]["path"] != expected
                or attestation["sources"][name]["sha256"]
                    != runtime_bindings[name]["sha256"]
            )
        ):
            raise campaign.CampaignError(
                f"executed {name} source is not manifest-bound")
    if runtime_bindings["runner"]["sha256"] not in \
            SUPPORTED_RUNNER_SHA256:
        raise campaign.CampaignError(
            "parallel wrapper does not support this campaign runner hash")


class _EvidenceByteLimitExceeded(Exception):
    pass


class _CappedHashReader:
    def __init__(self, source, limit):
        self.source = source
        self.limit = limit
        self.count = 0
        self.digest = hashlib.sha256()

    def read(self, size=-1):
        remaining_with_probe = self.limit - self.count + 1
        if size is None or size < 0 or size > remaining_with_probe:
            size = remaining_with_probe
        data = self.source.read(size)
        if self.count + len(data) > self.limit:
            raise _EvidenceByteLimitExceeded
        self.count += len(data)
        self.digest.update(data)
        return data

    def seek(self, offset, whence=os.SEEK_SET):
        if offset != self.count or whence != os.SEEK_SET:
            raise OSError("capped gzip evidence reader is forward-only")
        return self.count

    def tell(self):
        return self.count


def _result_evidence_byte_limits(campaign, job):
    if not isinstance(job, dict):
        raise campaign.CampaignError("result evidence job is malformed")
    replicates = job.get("replicates")
    warmups = job.get("warmup_replicates")
    if (
        type(replicates) is not int
        or replicates < 1
        or type(warmups) is not int
        or warmups < 0
    ):
        raise campaign.CampaignError(
            "result evidence replicate count is invalid")
    records = replicates + warmups
    return {
        "compressed":
            RESULT_COMPRESSED_FIXED_BYTES
            + records * RESULT_COMPRESSED_BYTES_PER_REPLICATE,
        "uncompressed":
            RESULT_UNCOMPRESSED_FIXED_BYTES
            + records * RESULT_UNCOMPRESSED_BYTES_PER_REPLICATE,
    }


def _read_canonical_gzip_json_capped(
        campaign, path, *, compressed_limit, uncompressed_limit):
    path = Path(path)
    stable = (
        "st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns",
        "st_ctime_ns",
    )
    try:
        with open(path, "rb") as source:
            before = os.fstat(source.fileno())
            if not stat.S_ISREG(before.st_mode):
                raise campaign.CampaignError(
                    f"gzip evidence is not a regular file: {path}")
            if before.st_size > compressed_limit:
                raise campaign.CampaignError(
                    f"gzip evidence exceeds compressed byte limit: {path}")
            reader = _CappedHashReader(source, compressed_limit)
            chunks = []
            uncompressed_bytes = 0
            with gzip.GzipFile(fileobj=reader, mode="rb") as archive:
                while True:
                    remaining_with_probe = (
                        uncompressed_limit - uncompressed_bytes + 1)
                    chunk = archive.read(
                        min(64 * 1024, remaining_with_probe))
                    if not chunk:
                        break
                    uncompressed_bytes += len(chunk)
                    if uncompressed_bytes > uncompressed_limit:
                        raise campaign.CampaignError(
                            "gzip evidence exceeds uncompressed byte "
                            f"limit: {path}")
                    chunks.append(chunk)
            while reader.read(64 * 1024):
                pass
            after = os.fstat(source.fileno())
        current = os.stat(path, follow_symlinks=False)
        if (
            any(
                getattr(before, name) != getattr(after, name)
                for name in stable
            )
            or any(
                getattr(after, name) != getattr(current, name)
                for name in stable
            )
            or reader.count != after.st_size
        ):
            raise campaign.CampaignError(
                f"gzip evidence changed while being read: {path}")
        payload = b"".join(chunks)
    except _EvidenceByteLimitExceeded:
        raise campaign.CampaignError(
            f"gzip evidence exceeds compressed byte limit: {path}")
    except (OSError, EOFError, gzip.BadGzipFile, zlib.error) as error:
        raise campaign.CampaignError(
            f"could not read gzip evidence {path}: {error}")
    return (
        campaign._strict_json_payload(payload, str(path)),
        {
            "compressed_sha256": reader.digest.hexdigest(),
            "uncompressed_sha256": campaign.sha256_bytes(payload),
            "compressed_bytes": reader.count,
            "uncompressed_bytes": uncompressed_bytes,
        },
    )


def _replay_result_capped(
        campaign, directory, manifest, manifest_sha256, assignment, job):
    path = Path(directory) / assignment["output"]
    limits = _result_evidence_byte_limits(campaign, job)
    envelope, evidence = _read_canonical_gzip_json_capped(
        campaign,
        path,
        compressed_limit=limits["compressed"],
        uncompressed_limit=limits["uncompressed"],
    )
    required = {
        "schema", "manifest_sha256", "job", "assignment",
        "wall_started_unix_ns", "wall_finished_unix_ns",
        "runtime_bindings_before", "runtime_bindings_after", "receipt",
    }
    if (
        not isinstance(envelope, dict)
        or set(envelope) != required
        or envelope["schema"] != campaign.RESULT_SCHEMA
        or envelope["manifest_sha256"] != manifest_sha256
        or not campaign.peel_codec._same_typed_json(
            envelope["job"], job)
        or not campaign.peel_codec._same_typed_json(
            envelope["assignment"], assignment)
        or not campaign.peel_codec._same_typed_json(
            envelope["runtime_bindings_before"],
            manifest["runtime_bindings"])
        or not campaign.peel_codec._same_typed_json(
            envelope["runtime_bindings_after"],
            manifest["runtime_bindings"])
        or type(envelope["wall_started_unix_ns"]) is not int
        or type(envelope["wall_finished_unix_ns"]) is not int
        or envelope["wall_finished_unix_ns"]
            <= envelope["wall_started_unix_ns"]
    ):
        raise campaign.CampaignError(
            f"result envelope is inconsistent: {path}")
    receipt = envelope["receipt"]
    request = campaign.expected_request(job)
    replayed = campaign.band.replay_bandtiming_receipt(
        receipt, expected_request=request)
    if (
        replayed.protocol != manifest["bandtiming_protocol"]
        or replayed.protocol != campaign.REQUIRED_BANDTIMING_PROTOCOL
        or not replayed.valid_for_promotion
    ):
        raise campaign.CampaignError(
            f"result receipt is non-promotable: {path}")
    if (
        receipt["candidate_descriptor"] != job["candidate_descriptor"]
        or receipt["dispatch_descriptor"] != job["dispatch_descriptor"]
        or len(receipt["replicates"]) != job["replicates"]
    ):
        raise campaign.CampaignError(
            f"native receipt changed its job: {path}")
    bound = receipt["context"]["bound"]
    if (
        bound["cpu_affinity"] != [assignment["cpu"]]
        or bound["thermal_device"]
            != manifest["runtime_bindings"]["thermal"]["device"]
        or bound["thermal_inode"]
            != manifest["runtime_bindings"]["thermal"]["inode"]
    ):
        raise campaign.CampaignError(
            f"native runtime binding is inconsistent: {path}")
    return receipt, evidence


def _projection(campaign, receipt):
    if (
        not isinstance(receipt, dict)
        or any(name not in receipt for name in _PROJECTION_FIELDS)
    ):
        raise campaign.CampaignError(
            "authenticated receipt projection is incomplete")
    return {name: receipt[name] for name in _PROJECTION_FIELDS}


def _initialize_worker(
        runner_path, directory, replay_manifest, manifest_sha256):
    global _WORKER_CONTEXT
    if _WORKER_CONTEXT is not None:
        raise RuntimeError("parallel verifier worker initialized twice")
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGTERM, signal.SIG_DFL)
    signal.signal(signal.SIGHUP, signal.SIG_DFL)
    campaign = _load_frozen_campaign(runner_path)
    campaign._PARENT_SIGNAL_STATE = None
    _assert_source_attestation(campaign)
    _assert_imported_runner_binding(
        campaign, replay_manifest["runtime_bindings"])
    campaign.check_runtime_bindings(
        replay_manifest["runtime_bindings"], full_hash=True)
    _WORKER_CONTEXT = (
        campaign,
        Path(directory),
        replay_manifest,
        manifest_sha256,
    )


def _replay_task(task):
    if _WORKER_CONTEXT is None:
        raise RuntimeError("parallel verifier worker is not initialized")
    campaign, directory, replay_manifest, manifest_sha256 = _WORKER_CONTEXT
    assignment, job = task
    try:
        receipt, evidence = _replay_result_capped(
            campaign,
            directory,
            replay_manifest,
            manifest_sha256,
            assignment,
            job,
        )
        return {
            "status": "ok",
            "order": assignment["order"],
            "job_id": job["job_id"],
            "receipt": _projection(campaign, receipt),
            "evidence": evidence,
        }
    except (
            campaign.CampaignError,
            campaign.peel_codec.MeasurementError,
            ValueError,
    ) as error:
        if isinstance(error, campaign.CampaignError):
            kind = "campaign"
        elif isinstance(error, campaign.peel_codec.MeasurementError):
            kind = "measurement"
        else:
            kind = "value"
        return {
            "status": "error",
            "order": assignment["order"],
            "job_id": job["job_id"],
            "exception_kind": kind,
            "message": str(error),
        }


def _worker_main(
        task_queue, result_queue, runner_path, directory,
        replay_manifest, manifest_sha256):
    _initialize_worker(
        runner_path, directory, replay_manifest, manifest_sha256)
    while True:
        task = task_queue.get()
        if task is None:
            return
        index, assignment, job = task
        result_queue.put((index, _replay_task((assignment, job))))


def _worker_count(campaign, manifest, replay_workers, job_count):
    value = manifest.get("workers") if replay_workers is None \
        else replay_workers
    if (
        type(value) is not int
        or not 1 <= value <= MAX_REPLAY_WORKERS
        or type(job_count) is not int
        or job_count < 1
    ):
        raise campaign.CampaignError("replay worker count is invalid")
    return min(value, job_count)


def _alive_processes(processes):
    return tuple(process for process in processes if process.is_alive())


def _join_processes_until(processes, deadline):
    for process in processes:
        if process.pid is None:
            continue
        remaining = deadline - time.monotonic()
        if remaining <= 0.0:
            break
        process.join(min(REPLAY_QUEUE_POLL_SECONDS, remaining))


def _terminate_processes(campaign, processes):
    for process in _alive_processes(processes):
        process.terminate()
    _join_processes_until(
        processes, time.monotonic() + REPLAY_TERMINATE_GRACE_SECONDS)
    for process in _alive_processes(processes):
        process.kill()
    _join_processes_until(
        processes, time.monotonic() + REPLAY_KILL_GRACE_SECONDS)
    survivors = _alive_processes(processes)
    if survivors:
        raise campaign.CampaignError(
            "parallel replay workers survived SIGKILL: "
            + ",".join(str(process.pid) for process in survivors))


def _close_queues(*queues):
    for work_queue in queues:
        try:
            work_queue.cancel_join_thread()
        except (AttributeError, ValueError):
            pass
        try:
            work_queue.close()
        except (AttributeError, ValueError):
            pass


def _worker_failure(
        campaign, processes, *, allow_clean_exit=False):
    failures = [
        (process.pid, process.exitcode)
        for process in processes
        if (
            process.exitcode is not None
            and (process.exitcode != 0 or not allow_clean_exit)
        )
    ]
    if failures:
        raise campaign.CampaignError(
            f"parallel replay worker failed: {failures}")


@contextmanager
def _parent_termination_handlers(campaign):
    """Make INT/TERM/HUP record-only until replay children are reaped."""
    handled = (signal.SIGINT, signal.SIGTERM, signal.SIGHUP)
    if campaign._PARENT_SIGNAL_STATE is not None:
        yield
        return
    previous = {
        item: signal.getsignal(item)
        for item in handled
    }
    state = {"pending": None}

    def interrupt(signum, unused_frame):
        if state["pending"] is None:
            state["pending"] = signum

    completed = False
    try:
        campaign._PARENT_SIGNAL_STATE = state
        for item in handled:
            signal.signal(item, interrupt)
        yield
        completed = True
    finally:
        old_mask = signal.pthread_sigmask(signal.SIG_BLOCK, handled)
        try:
            pending = state["pending"]
            for item in handled:
                signal.signal(item, previous[item])
            campaign._PARENT_SIGNAL_STATE = None
        finally:
            signal.pthread_sigmask(signal.SIG_SETMASK, old_mask)
        if completed and pending is not None:
            raise campaign._parent_signal_error(pending)


def _put_task(
        campaign, task_queue, task, processes, *,
        allow_clean_exit=False):
    while True:
        campaign._parent_signal_safe_point()
        _worker_failure(
            campaign,
            processes,
            allow_clean_exit=allow_clean_exit,
        )
        try:
            task_queue.put(
                task, timeout=REPLAY_QUEUE_POLL_SECONDS)
            return
        except queue.Full:
            continue


def _finish_processes(campaign, task_queue, processes):
    try:
        for unused in processes:
            _put_task(
                campaign,
                task_queue,
                None,
                processes,
                allow_clean_exit=True,
            )
        deadline = time.monotonic() + REPLAY_TERMINATE_GRACE_SECONDS
        while _alive_processes(processes):
            campaign._parent_signal_safe_point()
            _worker_failure(
                campaign,
                processes,
                allow_clean_exit=True,
            )
            if time.monotonic() >= deadline:
                break
            _join_processes_until(
                processes,
                min(
                    deadline,
                    time.monotonic() + REPLAY_QUEUE_POLL_SECONDS,
                ),
            )
    finally:
        _terminate_processes(campaign, processes)


def _validate_result(campaign, result, assignment, job):
    if (
        not isinstance(result, dict)
        or type(result.get("status")) is not str
        or result["status"] not in ("ok", "error")
        or type(result.get("order")) is not int
        or result.get("order") != assignment["order"]
        or type(result.get("job_id")) is not str
        or result.get("job_id") != job["job_id"]
    ):
        raise campaign.CampaignError(
            "parallel replay result is inconsistent")
    if result["status"] == "error":
        if (
            set(result) != {
                "status", "order", "job_id", "exception_kind", "message",
            }
            or result["exception_kind"]
                not in ("campaign", "measurement", "value")
            or type(result["message"]) is not str
        ):
            raise campaign.CampaignError(
                "parallel replay result is inconsistent")
        exception_type = {
            "campaign": campaign.CampaignError,
            "measurement": campaign.peel_codec.MeasurementError,
            "value": ValueError,
        }[result["exception_kind"]]
        raise exception_type(result["message"])
    if (
        set(result) != {
            "status", "order", "job_id", "receipt", "evidence",
        }
        or not isinstance(result["receipt"], dict)
        or set(result["receipt"]) != set(_PROJECTION_FIELDS)
        or not isinstance(result["evidence"], dict)
        or set(result["evidence"]) != {
            "compressed_sha256", "uncompressed_sha256",
            "compressed_bytes", "uncompressed_bytes",
        }
    ):
        raise campaign.CampaignError(
            "parallel replay result is inconsistent")
    return result["receipt"], result["evidence"]


@contextmanager
def _ordered_results(
        campaign, runner_path, directory, manifest, manifest_sha256, *,
        replay_workers=None):
    assignments = tuple(manifest["assignments"])
    jobs = tuple(manifest["pre_cpu_jobs"])
    if len(assignments) != len(jobs) or not jobs:
        raise campaign.CampaignError(
            "campaign replay population is incomplete")
    workers = _worker_count(
        campaign, manifest, replay_workers, len(jobs))
    if workers == 1:
        def serial_results():
            for assignment, job in zip(assignments, jobs):
                campaign._parent_signal_safe_point()
                receipt, evidence = _replay_result_capped(
                    campaign,
                    directory,
                    manifest,
                    manifest_sha256,
                    assignment,
                    job,
                )
                yield _projection(campaign, receipt), evidence

        yield serial_results()
        return

    replay_manifest = {
        "bandtiming_protocol": manifest["bandtiming_protocol"],
        "runtime_bindings": manifest["runtime_bindings"],
    }
    context = multiprocessing.get_context("fork")
    capacity = min(len(jobs), workers * 2)
    task_queue = context.Queue(maxsize=capacity)
    result_queue = context.Queue(maxsize=capacity)
    processes = []
    completed = [False]
    handler = (
        _parent_termination_handlers(campaign)
        if campaign._PARENT_SIGNAL_STATE is None
        else nullcontext()
    )
    try:
        with handler:
            try:
                for index in range(workers):
                    campaign._parent_signal_safe_point()
                    process = context.Process(
                        target=_worker_main,
                        args=(
                            task_queue,
                            result_queue,
                            str(Path(runner_path).resolve()),
                            str(Path(directory).resolve()),
                            replay_manifest,
                            manifest_sha256,
                        ),
                        name=f"za5v-wrapper-replay-{index:02d}",
                    )
                    processes.append(process)
                    process.start()
                    campaign._parent_signal_safe_point()

                next_submit = 0
                for unused in range(capacity):
                    _put_task(
                        campaign,
                        task_queue,
                        (
                            next_submit,
                            assignments[next_submit],
                            jobs[next_submit],
                        ),
                        processes,
                    )
                    next_submit += 1

                def parallel_results():
                    nonlocal next_submit
                    buffered = {}
                    next_expected = 0
                    while next_expected < len(jobs):
                        while next_expected not in buffered:
                            campaign._parent_signal_safe_point()
                            _worker_failure(campaign, processes)
                            try:
                                index, result = result_queue.get(
                                    timeout=REPLAY_QUEUE_POLL_SECONDS)
                            except queue.Empty:
                                continue
                            if (
                                type(index) is not int
                                or not 0 <= index < next_submit
                                or index in buffered
                                or index < next_expected
                            ):
                                raise campaign.CampaignError(
                                    "parallel replay result order is "
                                    "inconsistent")
                            buffered[index] = _validate_result(
                                campaign,
                                result,
                                assignments[index],
                                jobs[index],
                            )
                        value = buffered.pop(next_expected)
                        next_expected += 1
                        yield value
                        if next_submit < len(jobs):
                            _put_task(
                                campaign,
                                task_queue,
                                (
                                    next_submit,
                                    assignments[next_submit],
                                    jobs[next_submit],
                                ),
                                processes,
                            )
                            next_submit += 1
                    completed[0] = True

                yield parallel_results()
                campaign._parent_signal_safe_point()
            finally:
                if completed[0]:
                    _finish_processes(
                        campaign, task_queue, processes)
                else:
                    _terminate_processes(campaign, processes)
    finally:
        _close_queues(task_queue, result_queue)


def parallel_build_summaries(
        campaign, runner_path, directory, manifest, output_directory, *,
        replay_workers=None):
    """Reproduce the frozen builder with only full job replay parallelized."""
    _assert_imported_runner_binding(
        campaign, manifest["runtime_bindings"])
    campaign.check_runtime_bindings(
        manifest["runtime_bindings"], full_hash=True)
    output_directory = Path(output_directory)
    output_directory.mkdir(parents=True, exist_ok=True)
    ledger_path = output_directory / "cell-ledger.jsonl.gz"
    summary_path = output_directory / "summary.json.gz"
    jobs = tuple(manifest["pre_cpu_jobs"])
    manifest_sha256 = campaign.canonical_sha256(manifest)
    job_files = campaign._verify_job_files(
        directory, manifest, manifest_sha256)
    aggregator = campaign.RecoveryAggregator(manifest["phase"], jobs)
    job_evidence = campaign.JobEvidenceAggregator(
        manifest["phase"], jobs)
    result_files = []
    stream_hashes = set()
    with _ordered_results(
            campaign,
            runner_path,
            directory,
            manifest,
            manifest_sha256,
            replay_workers=replay_workers,
    ) as replayed:
        with campaign.AtomicGzipJsonLines(ledger_path) as ledger:
            for assignment, job in zip(manifest["assignments"], jobs):
                campaign._parent_signal_safe_point()
                try:
                    receipt, evidence = next(replayed)
                except StopIteration:
                    raise campaign.CampaignError(
                        "parallel replay omitted a result")
                stream = receipt["stream_sha256"]
                if stream in stream_hashes:
                    raise campaign.CampaignError(
                        "duplicate native stream hash")
                stream_hashes.add(stream)
                job_evidence.add(job, receipt)
                result_files.append({
                    "path": assignment["output"],
                    **evidence,
                    "stream_sha256": stream,
                })
                for replicate in receipt["replicates"]:
                    ledger.write(
                        aggregator.add_replicate(job, replicate))
            try:
                next(replayed)
            except StopIteration:
                pass
            else:
                raise campaign.CampaignError(
                    "parallel replay produced extra results")
    ledger_evidence = ledger.evidence()
    aggregates = aggregator.finish()
    timing_and_thermal = job_evidence.finish()
    decision = campaign.derive_campaign_decision(
        manifest["phase"],
        aggregates,
        timing_and_thermal,
        survivor=manifest["survivor"],
        selection_production_blockers=(
            manifest["selection_contract"][
                "selected_candidate_production_blockers"]
            if manifest["phase"] == "holdout" else None
        ),
    )
    if timing_and_thermal["thermal_context"]["cpu_assignments"] != \
            manifest["worker_cpus"]:
        raise campaign.CampaignError(
            "thermal evidence did not cover every assigned worker CPU")
    summary = {
        "schema": campaign.SUMMARY_SCHEMA,
        "phase": manifest["phase"],
        "survivor": manifest["survivor"],
        "frozen_roster_sha256": manifest["frozen_roster_sha256"],
        "pre_cpu_job_list_sha256":
            manifest["pre_cpu_job_list_sha256"],
        "job_count": len(jobs),
        "cell_ledger": {
            "path": "cell-ledger.jsonl.gz",
            **ledger_evidence,
        },
        "decision": decision,
        **aggregates,
        **timing_and_thermal,
    }
    summary_evidence = campaign.atomic_gzip_json(
        summary_path, summary)
    return {
        "job_files": job_files,
        "result_files": result_files,
        "cell_ledger": {
            "path": "cell-ledger.jsonl.gz",
            **ledger_evidence,
        },
        "summary": {
            "path": "summary.json.gz",
            **summary_evidence,
        },
    }


def verify_with_parallel_replay(
        runner_path, directory, *, replay_workers=None):
    """Call the frozen verifier after one fail-closed in-memory patch."""
    runner_path = Path(runner_path).resolve()
    directory = Path(directory).resolve()
    campaign = _load_frozen_campaign(runner_path)
    _assert_source_attestation(campaign)
    manifest, unused_manifest_sha256 = campaign.read_canonical_json(
        directory / "manifest.json")
    campaign._validate_manifest(manifest)
    _assert_imported_runner_binding(
        campaign, manifest["runtime_bindings"])
    campaign.check_runtime_bindings(
        manifest["runtime_bindings"], full_hash=True)
    original_build_summaries = campaign.build_summaries

    def patched_build_summaries(
            replay_directory, replay_manifest, output_directory):
        if (
            Path(replay_directory).resolve() != directory
            or not campaign.peel_codec._same_typed_json(
                replay_manifest, manifest)
        ):
            raise campaign.CampaignError(
                "verified manifest changed before parallel replay")
        return parallel_build_summaries(
            campaign,
            runner_path,
            replay_directory,
            replay_manifest,
            output_directory,
            replay_workers=replay_workers,
        )

    campaign.build_summaries = patched_build_summaries
    try:
        with _parent_termination_handlers(campaign):
            verified = campaign.verify_campaign(directory)
            campaign._parent_signal_safe_point()
    finally:
        campaign.build_summaries = original_build_summaries
    _assert_source_attestation(campaign)
    _assert_imported_runner_binding(
        campaign, manifest["runtime_bindings"])
    campaign.check_runtime_bindings(
        manifest["runtime_bindings"], full_hash=True)
    return campaign, verified


def _parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--runner", required=True, type=Path,
        help="exact manifest-bound wh2_za5v_campaign.py")
    parser.add_argument(
        "--directory", required=True, type=Path,
        help="completed campaign directory")
    parser.add_argument(
        "--replay-workers", type=int, default=None,
        help="bounded replay processes (default: manifest worker count)")
    return parser


def main(argv=None):
    args = _parser().parse_args(argv)
    campaign, verified = verify_with_parallel_replay(
        args.runner,
        args.directory,
        replay_workers=args.replay_workers,
    )
    print(
        campaign.canonical_json_bytes(verified).decode("ascii"),
        end="",
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(2)
