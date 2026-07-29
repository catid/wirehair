#!/usr/bin/env python3
"""Parallel strict verifier for an immutable WH2 za5v campaign runner.

This wrapper imports the manifest-bound ``wh2_za5v_campaign.py`` by its exact
path, verifies that runtime binding before and after replay, and monkey-patches
only ``build_summaries`` in memory.  Every result still passes the frozen
runner's complete ``_replay_result`` validation.  Only independent job replay
is parallel; aggregation and artifact generation retain manifest order.
"""

import argparse
from collections import deque
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
import hashlib
import importlib.util
import multiprocessing
import os
from pathlib import Path
import stat
import sys


MAX_REPLAY_WORKERS = 32
SUPPORTED_RUNNER_SHA256 = frozenset({
    "bae4336508d9d5fd9dd514950da709cb01cc129dcf9e0b36f3e8fe96e17346ef",
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


def _runner_snapshot_identity(path):
    """Hash one stable regular-file snapshot before importing executable code."""
    path = Path(path).resolve()
    try:
        with open(path, "rb") as source:
            before = os.fstat(source.fileno())
            if not stat.S_ISREG(before.st_mode):
                raise RuntimeError(
                    f"frozen campaign runner is not regular: {path}")
            digest = hashlib.sha256()
            while True:
                chunk = source.read(1024 * 1024)
                if not chunk:
                    break
                digest.update(chunk)
            after = os.fstat(source.fileno())
        current = os.stat(path, follow_symlinks=False)
    except OSError as error:
        raise RuntimeError(
            f"could not snapshot frozen campaign runner {path}: {error}")
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
            "frozen campaign runner changed while being snapshotted")
    return tuple(getattr(after, name) for name in stable), digest.hexdigest()


def _load_frozen_campaign(path):
    path = Path(path).resolve()
    before = _runner_snapshot_identity(path)
    module_name = (
        "_wirehair_frozen_wh2_za5v_"
        + hashlib.sha256(str(path).encode("utf-8")).hexdigest()[:12]
        + "_"
        + before[1]
    )
    existing = sys.modules.get(module_name)
    if existing is not None:
        if Path(existing.__file__).resolve() != path:
            raise RuntimeError("frozen campaign module name collision")
        if _runner_snapshot_identity(path) != before:
            raise RuntimeError(
                "frozen campaign runner changed before cached import")
        return existing
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not import frozen campaign runner {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    previous_dont_write_bytecode = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        spec.loader.exec_module(module)
    except BaseException:
        del sys.modules[module_name]
        raise
    finally:
        sys.dont_write_bytecode = previous_dont_write_bytecode
    if (
        Path(module.__file__).resolve() != path
        or _runner_snapshot_identity(path) != before
    ):
        raise RuntimeError("imported campaign runner path changed")
    return module


def _assert_imported_runner_binding(campaign, runtime_bindings):
    campaign._validate_runtime_bindings_schema(runtime_bindings)
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
    if runtime_bindings["runner"]["sha256"] not in \
            SUPPORTED_RUNNER_SHA256:
        raise campaign.CampaignError(
            "parallel wrapper does not support this campaign runner hash")


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
    campaign = _load_frozen_campaign(runner_path)
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
        receipt, evidence = campaign._replay_result(
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
                receipt, evidence = campaign._replay_result(
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
    executor = ProcessPoolExecutor(
        max_workers=workers,
        mp_context=multiprocessing.get_context("fork"),
        initializer=_initialize_worker,
        initargs=(
            str(Path(runner_path).resolve()),
            str(Path(directory).resolve()),
            replay_manifest,
            manifest_sha256,
        ),
    )
    pending = deque()
    population = iter(zip(assignments, jobs))

    def submit_one():
        try:
            assignment, job = next(population)
        except StopIteration:
            return False
        pending.append((
            assignment,
            job,
            executor.submit(_replay_task, (assignment, job)),
        ))
        return True

    try:
        for unused in range(min(len(jobs), workers * 2)):
            if not submit_one():
                break

        def parallel_results():
            while pending:
                campaign._parent_signal_safe_point()
                assignment, job, future = pending.popleft()
                result = future.result()
                campaign._parent_signal_safe_point()
                yield _validate_result(
                    campaign, result, assignment, job)
                submit_one()

        yield parallel_results()
    finally:
        for unused_assignment, unused_job, future in pending:
            future.cancel()
        executor.shutdown(wait=True, cancel_futures=True)


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
        verified = campaign.verify_campaign(directory)
    finally:
        campaign.build_summaries = original_build_summaries
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
