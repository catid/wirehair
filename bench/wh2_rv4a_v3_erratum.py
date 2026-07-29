#!/usr/bin/env python3
"""Derive an authenticated repair-v1 classification erratum from v3 evidence."""

import argparse
import collections
import concurrent.futures
import hashlib
import multiprocessing
import os
from pathlib import Path
import stat
import subprocess
import sys
import types


SOURCE_BYTE_LIMIT = 16 * 1024 * 1024
REPOSITORY = Path(__file__).resolve().parents[1]
_BOOTSTRAP_ERRATUM_SOURCE_SHA256 = globals().get(
    "_BOOTSTRAP_ERRATUM_SOURCE_SHA256")


def _source_bytes(path):
    path = Path(path).resolve()
    flags = os.O_RDONLY | os.O_CLOEXEC | getattr(os, "O_NONBLOCK", 0)
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        descriptor = os.open(str(path), flags)
        try:
            before = os.fstat(descriptor)
            if (
                not stat.S_ISREG(before.st_mode) or
                before.st_size > SOURCE_BYTE_LIMIT
            ):
                raise RuntimeError(
                    f"Python source is invalid or oversized: {path}")
            chunks = []
            size = 0
            while True:
                chunk = os.read(
                    descriptor,
                    min(1024 * 1024, SOURCE_BYTE_LIMIT - size + 1),
                )
                if not chunk:
                    break
                size += len(chunk)
                if size > SOURCE_BYTE_LIMIT:
                    raise RuntimeError(
                        f"Python source is oversized: {path}")
                chunks.append(chunk)
            after = os.fstat(descriptor)
        finally:
            os.close(descriptor)
        current = os.stat(path, follow_symlinks=False)
    except OSError as error:
        raise RuntimeError(f"could not read Python source {path}: {error}")
    stable = (
        "st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns",
        "st_ctime_ns",
    )
    payload = b"".join(chunks)
    if (
        any(getattr(before, name) != getattr(after, name)
            for name in stable) or
        any(getattr(after, name) != getattr(current, name)
            for name in stable) or
        size != after.st_size or len(payload) != size
    ):
        raise RuntimeError(f"Python source changed while read: {path}")
    return payload


def _source_module(module_name, path):
    path = Path(path).resolve()
    payload = _source_bytes(path)
    digest = hashlib.sha256(payload).hexdigest()
    module = types.ModuleType(module_name)
    module.__file__ = str(path)
    module.__package__ = ""
    module.__loader__ = None
    module.__spec__ = None
    sys.modules[module_name] = module
    try:
        code = compile(payload, str(path), "exec", dont_inherit=True)
        exec(code, module.__dict__)
        module.__rv4a_source_sha256__ = digest
        if _source_bytes(path) != payload:
            raise RuntimeError(
                f"Python source changed during import: {path}")
    except BaseException:
        sys.modules.pop(module_name, None)
        raise
    return module


campaign = _source_module(
    "wh2_rv4a_campaign",
    REPOSITORY / "bench/wh2_rv4a_campaign.py",
)


ERRATUM_SCHEMA = "wirehair.wh2.rv4a.repair-classification-erratum.v1"
POLICY_SCHEMA = "wirehair.wh2.rv4a.repair-classification-policy.v1"
_WORKER_STATE = None
classification = campaign._bootstrap_source_module(
    "repair_v1_classification",
    campaign.REPOSITORY / "tools/repair_v1_classification.py",
)


class ErratumError(RuntimeError):
    """Raised when immutable evidence cannot support the erratum."""


def classification_policy():
    """Return the versioned additive classification contract."""
    return {
        "schema": POLICY_SCHEMA,
        "source": {
            "repairtiming_protocol":
                campaign.REQUIRED_REPAIRTIMING_PROTOCOL,
            "repairtiming_schema":
                campaign.REQUIRED_REPAIRTIMING_SCHEMA,
            "rule":
                "strictly-replay-immutable-v3-receipts-without-rewriting-"
                "the-v1-summary-decision-ledger-or-completion",
        },
        "selector_weighting": {
            "key": ["arm_provisional_id", "K", "construction_root"],
            "witness":
                "lowest-replicate-for-the-key-in-the-first-frozen-"
                "schedule",
            "weak_count": "exactly-once-per-selector-key",
        },
        "raw_attempt_zero_structural_weakness": {
            "explicit_weak":
                "raw-and-attempt0-real-result-are-identical-code3-or-code4",
            "need_more":
                "raw-and-attempt0-real-result-are-identical-code1",
            "corroborated_error":
                "raw-and-attempt0-real-result-are-identical-code8-and-"
                "attempt0-structural-probe-executed-with-code1-and-no-"
                "fatal-attempt0-mismatch-and-selected-attempt1-through7-"
                "has-matching-probe-and-real-code0-and-the-selector-is-"
                "committed-with-bound-descriptor-and-cap-fatal-oom-zero-"
                "and-repaired-recovery-is-successful",
            "other_error":
                "code8-without-that-exact-corroboration-is-not-a-"
                "structural-weakness",
        },
        "runtime_error": {
            "weight": "once-per-physical-recovery-cell",
            "expected_code8_exemptions": sorted(
                classification.EXPECTED_ERROR_LOCATIONS),
            "exemption_rule":
                "exempt-only-both-named-duplicate-observations-of-one-"
                "corroborated-attempt0-error-with-a-bound-later-success",
            "fail_closed":
                "every-other-code8-and-every-result-code-outside-"
                "minus1-0-1-3-4-7-9-is-a-runtime-error",
        },
        "repaired_final_weak": {
            "unique_selectors":
                "count-once-per-verified-selector-key",
            "physical_cells":
                "count-once-per-fully-replayed-recovery-cell",
            "historical_check":
                "both-counts-must-equal-the-corresponding-v1-summary-"
                "values-because-this-erratum-does-not-reclassify-final-"
                "outcomes",
        },
        "outcome_use":
            "reporting-erratum-only-this-historical-v1-no-survivor-"
            "decision-and-this-erratum-grant-no-sealed-or-public-"
            "promotion-authorization",
    }


def _strict_replay_job(directory, manifest, manifest_sha256, assignment, job):
    parser_module = campaign._require_v3_parser()
    limits = campaign._result_evidence_byte_limits(job)
    path = directory / assignment["output"]
    envelope, evidence = campaign.read_canonical_gzip_json(
        path,
        compressed_limit=limits["compressed"],
        uncompressed_limit=limits["uncompressed"],
    )
    expected_fields = {
        "schema", "manifest_sha256", "job", "assignment",
        "wall_started_unix_ns", "wall_finished_unix_ns",
        "runtime_bindings_before", "runtime_bindings_after", "receipt",
    }
    if (
        not isinstance(envelope, dict) or
        set(envelope) != expected_fields or
        envelope.get("schema") != campaign.RESULT_SCHEMA or
        envelope.get("manifest_sha256") != manifest_sha256 or
        not campaign.peel_codec._same_typed_json(
            envelope.get("job"), job) or
        not campaign.peel_codec._same_typed_json(
            envelope.get("assignment"), assignment) or
        type(envelope.get("wall_started_unix_ns")) is not int or
        type(envelope.get("wall_finished_unix_ns")) is not int or
        not 0 < envelope["wall_started_unix_ns"] <= \
            envelope["wall_finished_unix_ns"] or
        not campaign.peel_codec._same_typed_json(
            envelope.get("runtime_bindings_before"),
            manifest["runtime_bindings"]) or
        not campaign.peel_codec._same_typed_json(
            envelope.get("runtime_bindings_after"),
            manifest["runtime_bindings"])
    ):
        raise ErratumError(f"result envelope is invalid: {path}")
    receipt = envelope["receipt"]
    replayed = parser_module.replay_repairtiming_receipt(
        receipt, expected_request=campaign.expected_request(job))
    if (
        getattr(replayed, "protocol", None) !=
            campaign.REQUIRED_REPAIRTIMING_PROTOCOL or
        not campaign.peel_codec._same_typed_json(
            campaign._measurement_as_dict(replayed), receipt)
    ):
        raise ErratumError("strict post-v3 replay changed the receipt")
    thermal = campaign._receipt_context_projection(receipt, assignment)
    if (
        thermal["thermal_device"] !=
            manifest["runtime_bindings"]["thermal"]["device"] or
        thermal["thermal_inode"] !=
            manifest["runtime_bindings"]["thermal"]["inode"]
    ):
        raise ErratumError("result receipt changed thermal identity")

    cells = getattr(replayed, "cells", None)
    if not isinstance(cells, tuple) or len(cells) != job["replicates"]:
        raise ErratumError("post-v3 replay cell count is invalid")
    selectors = {}
    runtime_cells = 0
    repaired_final_weak_cells = 0
    runtime_observations = collections.Counter()
    for replicate, cell in enumerate(cells):
        if (
            type(getattr(cell, "replicate", None)) is not int or
            cell.replicate != replicate or
            type(getattr(cell, "construction_root", None)) is not int or
            cell.construction_root !=
                campaign._expected_root(job, replicate) or
            type(getattr(cell, "loss_seed", None)) is not int or
            cell.loss_seed != campaign._expected_loss(job, replicate)
        ):
            raise ErratumError("post-v3 cell seed/order changed")
        key = (
            job["arm_provisional_id"], job["K"], cell.construction_root)
        if getattr(cell, "selector_key", None) != key:
            raise ErratumError("post-v3 selector key changed")
        selector = parser_module.selector_projection(cell)
        real = parser_module.real_projection(cell)
        cell_classification = classification.classify_repair_v1_cell(
            selector, real, parser_module=parser_module)
        prior = selectors.get(key)
        if prior is None:
            selectors[key] = {
                "projection": selector,
                "classification": cell_classification,
                "replicate": replicate,
            }
        elif not campaign.peel_codec._same_typed_json(
                prior["projection"], selector):
            raise ErratumError(
                "duplicate selector key changed structural projection")
        runtime_cells += \
            cell_classification["candidate_runtime_error"]
        repaired_final_weak_cells += \
            cell_classification["repaired_final_weak"]
        runtime_observations.update(
            cell_classification["runtime_error_observations"])

    selector_rows = [
        {
            "selector_key": list(key),
            "selector_projection": selectors[key]["projection"],
        }
        for key in sorted(selectors)
    ]
    selector_classifications = [
        {
            "selector_key": list(key),
            "replicate": selectors[key]["replicate"],
            **selectors[key]["classification"],
        }
        for key in sorted(selectors)
    ]
    result_evidence = {
        "order": assignment["order"],
        "job_id": job["job_id"],
        "path": assignment["output"],
        **evidence,
        "wall_started_unix_ns": envelope["wall_started_unix_ns"],
        "wall_finished_unix_ns": envelope["wall_finished_unix_ns"],
        "thermal": thermal,
    }
    return {
        "order": assignment["order"],
        "job_id": job["job_id"],
        "arm": job["arm"],
        "K": job["K"],
        "schedule": job["schedule"],
        "lane": job["lane"],
        "cells": len(cells),
        "result_evidence": result_evidence,
        "selector_set_sha256": campaign.canonical_sha256(selector_rows),
        "unique_selectors": len(selector_rows),
        "selector_classifications": selector_classifications,
        "candidate_runtime_error": runtime_cells,
        "repaired_final_weak_cells": repaired_final_weak_cells,
        "runtime_error_observations": dict(
            sorted(runtime_observations.items())),
    }


def _initialize_worker(directory, manifest, manifest_sha256):
    global _WORKER_STATE
    _WORKER_STATE = (Path(directory), manifest, manifest_sha256)


def _worker(task):
    if _WORKER_STATE is None:
        raise ErratumError("post-v3 replay worker was not initialized")
    assignment, job = task
    return _strict_replay_job(
        *_WORKER_STATE, assignment, job)


def _read_completion_source(directory):
    directory = Path(directory).resolve()
    manifest, manifest_sha256 = campaign.read_canonical_json(
        directory / "manifest.json")
    jobs = campaign._validate_manifest(manifest)
    if (
        manifest.get("phase") != "training" or
        manifest.get("survivor") is not None
    ):
        raise ErratumError(
            "classification erratum requires a training completion")
    completion, completion_sha256 = campaign.read_canonical_json(
        directory / "completion.json")
    checksum = campaign._read_stable_bytes(
        directory / "completion.sha256", byte_limit=256)
    if checksum != (
            f"{completion_sha256}  completion.json\n".encode("ascii")):
        raise ErratumError("completion checksum does not match")
    completion_fields = {
        "schema", "phase", "survivor", "manifest_sha256",
        "frozen_roster_sha256", "pre_cpu_job_list_sha256",
        "cell_set_sha256", "policy_sha256", "preflight_pin_sha256", "jobs",
        "runtime_binding_monitor", "runtime_bindings_final",
        "job_files", "log_files", "result_files", "cell_ledger", "summary",
        "training_decision_sha256", "sealed_confirmation_sha256",
        "completed_unix_ns",
    }
    if (
        not isinstance(completion, dict) or
        set(completion) != completion_fields or
        completion.get("schema") != campaign.COMPLETION_SCHEMA or
        completion.get("phase") != manifest["phase"] or
        completion.get("survivor") != manifest["survivor"] or
        completion.get("manifest_sha256") != manifest_sha256 or
        completion.get("frozen_roster_sha256") !=
            manifest["frozen_roster_sha256"] or
        completion.get("pre_cpu_job_list_sha256") !=
            manifest["pre_cpu_job_list_sha256"] or
        completion.get("cell_set_sha256") !=
            manifest["cell_set_sha256"] or
        completion.get("policy_sha256") != manifest["policy_sha256"] or
        completion.get("preflight_pin_sha256") !=
            manifest["preflight_pin_sha256"] or
        completion.get("jobs") != len(jobs) or
        not campaign.peel_codec._same_typed_json(
            completion.get("runtime_bindings_final"),
            manifest["runtime_bindings"]) or
        type(completion.get("completed_unix_ns")) is not int or
        completion["completed_unix_ns"] < manifest["created_unix_ns"]
    ):
        raise ErratumError("completion does not match manifest")
    campaign._validate_runtime_monitor(
        completion["runtime_binding_monitor"], manifest)
    summary = campaign._verify_stored_artifacts(directory, completion)
    decision = summary.get("training_decision")
    confirmation = summary.get("sealed_confirmation")
    if (
        not isinstance(decision, dict) or
        decision.get("schema") != campaign.TRAINING_DECISION_SCHEMA or
        decision.get("policy_sha256") != manifest["policy_sha256"] or
        decision.get("status") != "no-survivor" or
        decision.get("selected_survivor") is not None or
        decision.get("sealed_authorization") != "forbidden" or
        confirmation is not None or
        completion["training_decision_sha256"] != (
            campaign.canonical_sha256(decision)) or
        completion["sealed_confirmation_sha256"] is not None
    ):
        raise ErratumError(
            "stored v1 no-survivor decision is invalid")
    if not campaign.peel_codec._same_typed_json(
            campaign._verify_job_files(
                directory, manifest, manifest_sha256),
            completion["job_files"]):
        raise ErratumError("stored job files changed")
    if not campaign.peel_codec._same_typed_json(
            campaign._verify_empty_logs(directory, manifest),
            completion["log_files"]):
        raise ErratumError("stored worker logs changed")
    ledger_rows, ledger_evidence = campaign.read_canonical_gzip_json_lines(
        directory / completion["cell_ledger"]["path"],
        compressed_limit=campaign.SUMMARY_COMPRESSED_BYTE_LIMIT,
        uncompressed_limit=campaign.SUMMARY_UNCOMPRESSED_BYTE_LIMIT,
        row_limit=len(jobs),
    )
    expected_ledger = {
        key: completion["cell_ledger"][key]
        for key in (
            "compressed_sha256", "uncompressed_sha256",
            "compressed_bytes", "uncompressed_bytes", "rows")
    }
    if (
        not campaign.peel_codec._same_typed_json(
            ledger_evidence, expected_ledger) or
        len(ledger_rows) != len(jobs)
    ):
        raise ErratumError("stored ledger changed")
    return {
        "directory": directory,
        "manifest": manifest,
        "manifest_sha256": manifest_sha256,
        "jobs": jobs,
        "completion": completion,
        "completion_sha256": completion_sha256,
        "summary": summary,
        "ledger_rows": ledger_rows,
    }


def _recheck_result_file(task):
    directory, expected = task
    if (
        not isinstance(expected, dict) or
        not isinstance(expected.get("path"), str) or
        type(expected.get("compressed_bytes")) is not int or
        expected["compressed_bytes"] < 1 or
        not campaign._is_sha256(expected.get("compressed_sha256"))
    ):
        raise ErratumError("result-file evidence is malformed")
    relative = Path(expected["path"])
    if relative.is_absolute() or ".." in relative.parts:
        raise ErratumError("result-file evidence path is invalid")
    digest, size = campaign._stream_file_sha256(
        Path(directory) / relative,
        byte_limit=expected["compressed_bytes"],
    )
    if (
        digest != expected["compressed_sha256"] or
        size != expected["compressed_bytes"]
    ):
        raise ErratumError(
            f"result file changed after replay: {relative}")
    return expected["path"]


def _recheck_completed_source(source, *, workers):
    current = _read_completion_source(source["directory"])
    for field in (
            "manifest", "manifest_sha256", "jobs", "completion",
            "completion_sha256", "summary", "ledger_rows"):
        if not campaign.peel_codec._same_typed_json(
                current[field], source[field]):
            raise ErratumError(
                f"completed evidence changed after replay: {field}")
    expected_results = source["completion"]["result_files"]
    tasks = [
        (str(source["directory"]), expected)
        for expected in expected_results
    ]
    with concurrent.futures.ThreadPoolExecutor(
            max_workers=min(workers, len(tasks))) as executor:
        paths = list(executor.map(_recheck_result_file, tasks))
    if paths != [
            expected["path"] for expected in expected_results]:
        raise ErratumError("result-file recheck changed order")
    return current


def _git_source_receipt(commit, relative_paths):
    repository = campaign.REPOSITORY.resolve()
    try:
        resolved = subprocess.run(
            [
                "git", "-C", str(repository), "rev-parse", "--verify",
                f"{commit}^{{commit}}",
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=30,
        ).stdout.strip()
        tree = subprocess.run(
            [
                "git", "-C", str(repository), "rev-parse", "--verify",
                f"{resolved}^{{tree}}",
            ],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=30,
        ).stdout.strip()
    except (OSError, subprocess.SubprocessError) as error:
        raise ErratumError(f"could not resolve source commit: {error}")
    if (
        len(resolved) != 40 or
        any(character not in "0123456789abcdef" for character in resolved) or
        len(tree) != 40 or
        any(character not in "0123456789abcdef" for character in tree)
    ):
        raise ErratumError("git returned an invalid source commit")
    files = {}
    for relative in relative_paths:
        if (
            not isinstance(relative, str) or
            Path(relative).is_absolute() or ".." in Path(relative).parts
        ):
            raise ErratumError("source receipt path is invalid")
        object_name = f"{resolved}:{relative}"
        try:
            size_text = subprocess.run(
                [
                    "git", "-C", str(repository), "cat-file", "-s",
                    object_name,
                ],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                timeout=30,
            ).stdout.strip()
            size = int(size_text, 10)
            if size < 0 or size > SOURCE_BYTE_LIMIT:
                raise ErratumError(
                    "committed source exceeds byte limit")
            payload = subprocess.run(
                [
                    "git", "-C", str(repository), "show",
                    object_name,
                ],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                timeout=30,
            ).stdout
        except (
                OSError, subprocess.SubprocessError, ValueError
        ) as error:
            raise ErratumError(
                f"could not read {relative} at source commit: {error}")
        if len(payload) != size:
            raise ErratumError("committed source size changed")
        files[relative] = {
            "bytes": len(payload),
            "sha256": hashlib.sha256(payload).hexdigest(),
        }
    return {
        "commit": resolved,
        "tree": tree,
        "files": files,
    }


def _validate_source_receipts(
        source, historical_source_commit, classifier_source_commit):
    manifest = source["manifest"]
    pinned_runner = manifest["runtime_bindings"]["sources"]["runner"]
    historical = _git_source_receipt(
        historical_source_commit, ["bench/wh2_rv4a_campaign.py"])
    historical_runner = historical["files"]["bench/wh2_rv4a_campaign.py"]
    if (
        historical_runner["sha256"] != pinned_runner["sha256"] or
        historical_runner["bytes"] != pinned_runner["size"]
    ):
        raise ErratumError(
            "historical commit does not contain the pinned v1 runner")

    classifier_paths = [
        "bench/wh2_rv4a_campaign.py",
        "bench/wh2_rv4a_v3_erratum.py",
        "tools/repair_v1_classification.py",
    ]
    classifier = _git_source_receipt(
        classifier_source_commit, classifier_paths)
    loaded_hashes = {
        "bench/wh2_rv4a_campaign.py":
            getattr(campaign, "__rv4a_source_sha256__", None),
        "bench/wh2_rv4a_v3_erratum.py":
            _BOOTSTRAP_ERRATUM_SOURCE_SHA256,
        "tools/repair_v1_classification.py":
            getattr(classification, "__rv4a_source_sha256__", None),
    }
    for relative in classifier_paths:
        disk = campaign._stable_file_binding(
            campaign.REPOSITORY / relative,
            byte_limit=SOURCE_BYTE_LIMIT,
        )
        committed = classifier["files"][relative]
        if (
            disk["sha256"] != committed["sha256"] or
            disk["size"] != committed["bytes"] or
            loaded_hashes[relative] != committed["sha256"]
        ):
            raise ErratumError(
                "loaded or on-disk classification source differs from "
                f"commit: {relative}")
    classifier_campaign = classifier["files"][
        "bench/wh2_rv4a_campaign.py"]
    if (
        classifier_campaign != historical_runner or
        classifier_campaign["sha256"] != pinned_runner["sha256"] or
        classifier_campaign["bytes"] != pinned_runner["size"]
    ):
        raise ErratumError(
            "classification runtime campaign helper is not the exact "
            "historical pinned runner")

    parser_module = campaign._require_v3_parser()
    for name, module in (
            ("parser", parser_module),
            ("context_tool", campaign.peel_codec)):
        expected = manifest["runtime_bindings"]["sources"][name]
        actual = campaign._stable_file_binding(
            Path(module.__file__), byte_limit=expected["byte_limit"])
        loaded_sha256 = getattr(
            module, "__rv4a_source_sha256__", None)
        if (
            actual["sha256"] != expected["sha256"] or
            actual["size"] != expected["size"] or
            loaded_sha256 != expected["sha256"]
        ):
            raise ErratumError(
                f"loaded v3 {name} differs from the runtime pin")
    return historical, classifier


def _new_arm_bucket(name):
    return {
        "name": name,
        "jobs": 0,
        "cells": 0,
        "selector_records": [],
        "candidate_runtime_error": 0,
        "repaired_final_weak_cells": 0,
        "runtime_error_observations": collections.Counter(),
    }


def _accumulate_job(bucket, replayed):
    bucket["jobs"] += 1
    bucket["cells"] += replayed["cells"]
    bucket["candidate_runtime_error"] += \
        replayed["candidate_runtime_error"]
    bucket["repaired_final_weak_cells"] += \
        replayed["repaired_final_weak_cells"]
    bucket["runtime_error_observations"].update(
        replayed["runtime_error_observations"])
    for selector in replayed["selector_classifications"]:
        bucket["selector_records"].append({
            "selector_key": selector["selector_key"],
            "schedule": replayed["schedule"],
            "replicate": selector["replicate"],
            "classification": {
                key: value for key, value in selector.items()
                if key not in ("selector_key", "replicate")
            },
        })


def _finish_bucket(bucket, legacy):
    selector_aggregate = \
        classification.aggregate_selector_classifications(
            bucket["selector_records"],
            schedules=campaign.SCHEDULES,
            first_schedule=campaign.SCHEDULES[0],
        )
    selector_aggregate["repaired_final_weak"]["physical_cells"] = \
        bucket["repaired_final_weak_cells"]
    return {
        **{
            key: value for key, value in bucket.items()
            if key not in (
                "runtime_error_observations", "selector_records",
                "repaired_final_weak_cells",
            )
        },
        **selector_aggregate,
        "runtime_error_observations": dict(
            sorted(bucket["runtime_error_observations"].items())),
        "legacy_v1_reported": legacy,
    }


def derive_erratum(
        directory, *, historical_source_commit,
        classifier_source_commit, workers):
    if type(workers) is not int or not 1 <= workers <= campaign.MAX_WORKERS:
        raise ErratumError("worker count is invalid")
    source = _read_completion_source(directory)
    historical_source, classifier_source = _validate_source_receipts(
        source, historical_source_commit, classifier_source_commit)
    manifest = source["manifest"]
    completion = source["completion"]
    summary = source["summary"]
    jobs = source["jobs"]
    assignments = manifest["assignments"]
    if len(completion["result_files"]) != len(jobs):
        raise ErratumError("completion result roster is incomplete")

    arm_buckets = {
        name: _new_arm_bucket(name)
        for name in summary["arms"]
    }
    overall = _new_arm_bucket("overall")
    selector_groups = {}
    context = multiprocessing.get_context("fork")
    tasks = list(zip(assignments, jobs))
    with context.Pool(
            min(workers, len(tasks)),
            initializer=_initialize_worker,
            initargs=(
                str(source["directory"]), manifest,
                source["manifest_sha256"],
            ),
            maxtasksperchild=32,
    ) as pool:
        replayed_jobs = pool.imap(_worker, tasks, chunksize=1)
        for index, replayed in enumerate(replayed_jobs):
            expected_result = completion["result_files"][index]
            ledger = source["ledger_rows"][index]
            if (
                replayed["order"] != index or
                replayed["job_id"] != jobs[index]["job_id"] or
                not campaign.peel_codec._same_typed_json(
                    replayed["result_evidence"], expected_result) or
                ledger.get("order") != index or
                ledger.get("job_id") != replayed["job_id"] or
                ledger.get("selector_set_sha256") !=
                    replayed["selector_set_sha256"] or
                ledger.get("unique_selectors") !=
                    replayed["unique_selectors"]
            ):
                raise ErratumError(
                    f"strict replay disagrees with job {index}")
            selector_key = (
                replayed["arm"], replayed["lane"], replayed["K"])
            value = (
                replayed["selector_set_sha256"],
                replayed["unique_selectors"],
            )
            prior = selector_groups.setdefault(
                selector_key, {
                    "value": value,
                    "schedules": set(),
                })
            if (
                prior["value"] != value or
                replayed["schedule"] in prior["schedules"]
            ):
                raise ErratumError(
                    "selector projection changed across schedules")
            prior["schedules"].add(replayed["schedule"])
            try:
                arm_bucket = arm_buckets[replayed["arm"]]
            except KeyError:
                raise ErratumError("replayed result has an unknown arm")
            _accumulate_job(arm_bucket, replayed)
            _accumulate_job(overall, replayed)

    if any(
            receipt["schedules"] != set(campaign.SCHEDULES)
            for receipt in selector_groups.values()):
        raise ErratumError("selector group lacks a frozen schedule")
    expected_unique = sum(
        receipt["value"][1] for receipt in selector_groups.values())
    if (
        expected_unique != summary["unique_selectors"] or
        overall["cells"] != summary["cells"]
    ):
        raise ErratumError("post-v3 aggregate cardinality changed")

    arms = {}
    for name, bucket in arm_buckets.items():
        legacy_summary = summary["arms"][name]
        finished = _finish_bucket(bucket, {
            "raw_attempt0_structural_weak":
                legacy_summary["selector"]["weak_constructions"][
                    "raw_attempt0"],
            "candidate_runtime_error":
                legacy_summary["candidate_runtime_error"],
            "repaired_final_weak_unique_selectors":
                legacy_summary["selector"]["weak_constructions"][
                    "repaired_final"],
            "repaired_final_weak_physical_cells":
                legacy_summary["candidate_final_weak"],
        })
        if (
            finished["cells"] != legacy_summary["cells"] or
            finished["unique_selectors"] !=
                legacy_summary["selector"]["observations"] or
            finished["repaired_final_weak"] != {
                "unique_selectors":
                    legacy_summary["selector"]["weak_constructions"][
                        "repaired_final"],
                "physical_cells":
                    legacy_summary["candidate_final_weak"],
            }
        ):
            raise ErratumError(
                f"post-v3 arm cardinality changed for {name}")
        arms[name] = finished
    finished_overall = _finish_bucket(overall, {
        "raw_attempt0_structural_weak":
            summary["overall"]["selector"]["weak_constructions"][
                "raw_attempt0"],
        "candidate_runtime_error":
            summary["overall"]["candidate_runtime_error"],
        "repaired_final_weak_unique_selectors":
            summary["overall"]["selector"]["weak_constructions"][
                "repaired_final"],
        "repaired_final_weak_physical_cells":
            summary["overall"]["candidate_final_weak"],
    })
    if (
        finished_overall["unique_selectors"] != summary["unique_selectors"] or
        finished_overall["repaired_final_weak"] != {
            "unique_selectors":
                summary["overall"]["selector"]["weak_constructions"][
                    "repaired_final"],
            "physical_cells": summary["overall"]["candidate_final_weak"],
        }
    ):
        raise ErratumError("post-v3 selector cardinality changed")
    _recheck_completed_source(source, workers=workers)
    historical_recheck, classifier_recheck = _validate_source_receipts(
        source, historical_source_commit, classifier_source_commit)
    if (
        not campaign.peel_codec._same_typed_json(
            historical_recheck, historical_source) or
        not campaign.peel_codec._same_typed_json(
            classifier_recheck, classifier_source)
    ):
        raise ErratumError("classification source receipt changed")
    policy = classification_policy()
    return {
        "schema": ERRATUM_SCHEMA,
        "policy": policy,
        "policy_sha256": campaign.canonical_sha256(policy),
        "source_evidence": {
            "manifest_schema": manifest["schema"],
            "manifest_sha256": source["manifest_sha256"],
            "completion_schema": completion["schema"],
            "completion_sha256": source["completion_sha256"],
            "summary_schema": summary["schema"],
            "summary": dict(completion["summary"]),
            "summary_evidence_sha256":
                campaign.canonical_sha256(completion["summary"]),
            "historical_training_decision": {
                "schema":
                    summary["training_decision"]["schema"],
                "sha256":
                    completion["training_decision_sha256"],
                "status":
                    summary["training_decision"]["status"],
                "selected_survivor":
                    summary["training_decision"][
                        "selected_survivor"],
                "sealed_authorization":
                    summary["training_decision"][
                        "sealed_authorization"],
            },
            "result_schema": campaign.RESULT_SCHEMA,
            "result_files": len(completion["result_files"]),
            "result_file_evidence_sha256":
                campaign.canonical_sha256(completion["result_files"]),
            "cell_ledger": dict(completion["cell_ledger"]),
            "cell_ledger_evidence_sha256":
                campaign.canonical_sha256(completion["cell_ledger"]),
            "historical_source": historical_source,
            "classification_source": classifier_source,
            "pinned_runner_sha256":
                manifest["runtime_bindings"]["sources"]["runner"]["sha256"],
            "pinned_parser_sha256":
                manifest["runtime_bindings"]["sources"]["parser"]["sha256"],
            "pinned_context_tool_sha256":
                manifest["runtime_bindings"]["sources"][
                    "context_tool"]["sha256"],
        },
        "phase": manifest["phase"],
        "survivor": manifest["survivor"],
        "jobs": len(jobs),
        "cells": summary["cells"],
        "unique_selectors": summary["unique_selectors"],
        "arms": arms,
        "overall": finished_overall,
        "outcome_authorization":
            "forbidden-v1-decision-retained-only-as-historical-evidence",
        "promotion_authorization":
            "forbidden-reporting-erratum-is-not-a-promotion-decision",
    }


def _parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--historical-source-commit", required=True)
    parser.add_argument("--classifier-source-commit", required=True)
    parser.add_argument(
        "--workers", type=int, default=campaign.MAX_WORKERS)
    return parser


def main(argv=None):
    args = _parser().parse_args(argv)
    receipt = derive_erratum(
        args.directory,
        historical_source_commit=args.historical_source_commit,
        classifier_source_commit=args.classifier_source_commit,
        workers=args.workers,
    )
    output = args.output.resolve()
    digest = campaign.atomic_json(output, receipt)
    print(campaign.canonical_json_bytes({
        "output": str(output),
        "schema": ERRATUM_SCHEMA,
        "sha256": digest,
    }).decode("ascii"), end="")


def _main_entry(argv):
    try:
        main(argv)
    except (
            ErratumError,
            campaign.CampaignError,
            classification.ClassificationError,
            ValueError,
    ) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    return 0


def _source_forced_cli(argv):
    path = Path(__file__).resolve()
    payload = _source_bytes(path)
    digest = hashlib.sha256(payload).hexdigest()
    module_name = f"_wirehair_rv4a_v3_erratum_{digest[:16]}"
    module = types.ModuleType(module_name)
    module.__file__ = str(path)
    module.__package__ = ""
    module.__loader__ = None
    module.__spec__ = None
    module._BOOTSTRAP_ERRATUM_SOURCE_SHA256 = digest
    sys.modules[module_name] = module
    try:
        code = compile(payload, str(path), "exec", dont_inherit=True)
        exec(code, module.__dict__)
        if module._source_bytes(path) != payload:
            raise RuntimeError(
                "erratum source changed during source-forced launch")
        result = module._main_entry(argv)
        if module._source_bytes(path) != payload:
            raise RuntimeError(
                "erratum source changed during source-forced execution")
        return result
    finally:
        sys.modules.pop(module_name, None)


if __name__ == "__main__":
    raise SystemExit(_source_forced_cli(sys.argv[1:]))
