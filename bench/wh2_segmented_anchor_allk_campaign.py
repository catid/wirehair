#!/usr/bin/env python3
"""Sealed raw all-K recovery comparison for segmented D12 anchor layouts.

This controller compares one raw D12/two-anchor baseline against the
experiment-only three-048, three-059, and four-0369 layouts.  It fixes the
mixed 10+2, P244, q0 geometry, uses no K-specific seed fixes, and covers every
K in 2..64000 under three independent seeds and the burst, adversarial, and
repair-only 50% loss schedules.

The campaign is recovery/work evidence, not quiet timing evidence.  ``run``
requires the normal nice-19 filler load and receipts the already-running
external thermal sampler without starting, stopping, or signalling it.
"""

from __future__ import annotations

import argparse
from collections import Counter
import concurrent.futures
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation, getcontext
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import threading
import time
from typing import Dict, List, Mapping, Optional, Sequence, Tuple


sys.dont_write_bytecode = True
getcontext().prec = 50

SCHEMA = "wirehair.wh2.segmented_anchor_allk.v1"
PROTOTYPE_COMMIT = "57d79903716bfb7bb815106040e478ed6996a029"
PROTOTYPE_PARENT = "a3f5c660b1d0de7be0f116c65adfd7eb23a9f617"
CONTROLLER_NAME = "wh2_segmented_anchor_allk_campaign.py"
BINARY_NAME = "wirehair_v2_bench.segmented_anchor"
K_LO = 2
K_HI = 64_000
CHUNK_MAX = 256
WORKERS = 56
LOSS = "0.5"
SEEDS = (
    "0x475b1dd3864534a3",
    "0x65eda454ed6be8ec",
    "0x6fd9d0c558cad523",
)
SCHEDULES = ("burst", "adversarial", "repair-only")
ARMS = ("baseline", "three-048", "three-059", "four-0369")
ARM_OPTIONS = {
    "baseline": ("--binary-dense-two-anchor",),
    "three-048": ("--binary-dense-segmented-anchors", "three-048"),
    "three-059": ("--binary-dense-segmented-anchors", "three-059"),
    "four-0369": ("--binary-dense-segmented-anchors", "four-0369"),
}
EXPECTED_PROTOTYPE_PATHS = (
    "codec/PrecodeEncodeTest.cpp",
    "codec/PrecodeSolveTest.cpp",
    "codec/PrecodeTest.cpp",
    "codec/README.md",
    "codec/V2BenchCliTest.cmake",
    "codec/V2NoInitTest.cpp",
    "codec/WirehairV2Bench.cpp",
    "codec/WirehairV2Precode.cpp",
    "codec/WirehairV2Precode.h",
    "codec/WirehairV2PrecodeEncode.cpp",
    "codec/fuzz/V2ProfileFuzz.cpp",
)
MAX_STDOUT_BYTES = 256 * 1024
MAX_STDERR_BYTES = 16 * 1024
TASK_TIMEOUT_SECONDS = 600
MIN_FILLER_COUNT = 56
FILLER_PID_FILE = Path("/tmp/wirehair-load-fillers.pid")
SAMPLER_BASENAME = "wirehair_expo_thermal_sampler.py"

CSV_FIELDS = (
    "N", "bb", "heavy_family", "mix_count", "overhead", "trials",
    "success", "rank_fail", "error", "fail_rate", "inact_mu",
    "inact_max", "binary_def_mu", "binary_def_max", "heavy_gain_mu",
    "heavy_gain_min", "heavy_shortfall", "solve_ms_mu", "build_ms_mu",
    "peel_ms_mu", "project_ms_mu", "residual_ms_mu", "backsub_ms_mu",
    "seed_attempt", "block_xors_mu", "block_muladds_mu",
    "first_rank_fail", "binary_def_hist", "heavy_gain_hist",
    "failure_trials", "active_packet_peel_seed_xor",
    "mixed_joint_source_xors_mu", "mixed_joint_marginal_xors_mu",
    "mixed_joint_marginal_copies_mu", "mixed_joint_active_deltas_mu",
    "mixed_joint_scratch_bytes_mu", "mixed_dual_source_columns_mu",
)
TIMING_FIELDS = (
    "solve_ms_mu", "build_ms_mu", "peel_ms_mu", "project_ms_mu",
    "residual_ms_mu", "backsub_ms_mu",
)
ZERO_JOINT_FIELDS = (
    "mixed_joint_source_xors_mu", "mixed_joint_marginal_xors_mu",
    "mixed_joint_marginal_copies_mu", "mixed_joint_active_deltas_mu",
    "mixed_joint_scratch_bytes_mu", "mixed_dual_source_columns_mu",
)
UINT_RE = re.compile(r"0|[1-9][0-9]*\Z")
MILLI_RE = re.compile(r"(0|[1-9][0-9]*)\.000\Z")
SHA_RE = re.compile(r"[0-9a-f]{64}\Z")
JOB_RECEIPT_FIELDS = frozenset((
    "schema", "self_sha256_excluding_field", "job", "task_sha256", "argv",
    "output_name", "started_utc", "start_monotonic_ns", "end_monotonic_ns",
    "duration_ns", "returncode", "row_count", "stdout_sha256",
    "stdout_bytes", "stderr_sha256", "stderr_bytes", "exit_sha256",
))
LAUNCH_RECEIPT_FIELDS = frozenset((
    "schema", "self_sha256_excluding_field", "started_utc", "ended_utc",
    "start_monotonic_ns", "end_monotonic_ns", "duration_ns", "design_sha256",
    "prepare_receipt_sha256", "tasks_manifest_sha256", "workers", "task_count",
    "job_receipts", "filler_start", "filler_end", "external_sampler_start",
    "external_sampler_end", "controller_started_sampler",
    "controller_signalled_sampler", "quiet_timing",
))


class CampaignError(RuntimeError):
    pass


def die(message: str) -> None:
    raise CampaignError(message)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds").replace(
        "+00:00", "Z")


def canonical_json(value: object) -> bytes:
    return (json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False,
    ) + "\n").encode("ascii")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_new(path: Path, data: bytes, mode: int = 0o444) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() or path.is_symlink():
        die("refusing to replace {}".format(path))
    part = path.with_name(
        path.name + ".part.{}.{}".format(os.getpid(), threading.get_ident()))
    try:
        with part.open("xb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.chmod(part, mode)
        os.link(part, path)
        directory = os.open(str(path.parent), os.O_RDONLY | os.O_DIRECTORY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    finally:
        try:
            part.unlink()
        except FileNotFoundError:
            pass


def sealed(schema: str, payload: Mapping[str, object]) -> Dict[str, object]:
    if "schema" in payload or "self_sha256_excluding_field" in payload:
        die("invalid sealed record construction")
    value: Dict[str, object] = {"schema": schema, **payload}
    value["self_sha256_excluding_field"] = sha256_bytes(canonical_json(value))
    return value


def verify_sealed(value: object, schema: str, name: str) -> Dict[str, object]:
    if not isinstance(value, dict) or value.get("schema") != schema:
        die("{} schema mismatch".format(name))
    claimed = value.get("self_sha256_excluding_field")
    if not isinstance(claimed, str) or SHA_RE.fullmatch(claimed) is None:
        die("{} self hash is malformed".format(name))
    unhashed = dict(value)
    del unhashed["self_sha256_excluding_field"]
    if sha256_bytes(canonical_json(unhashed)) != claimed:
        die("{} self hash mismatch".format(name))
    return value


def load_canonical(path: Path, name: str) -> Dict[str, object]:
    raw = path.read_bytes()
    try:
        value = json.loads(raw.decode("ascii"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        die("{} is not ASCII JSON: {}".format(name, exc))
    if not isinstance(value, dict) or canonical_json(value) != raw:
        die("{} is not a canonical JSON object".format(name))
    return value


def checked(
    argv: Sequence[str], *, cwd: Optional[Path] = None,
    environment: Optional[Mapping[str, str]] = None, timeout: float = 600,
) -> Tuple[bytes, bytes]:
    try:
        result = subprocess.run(
            list(argv), cwd=str(cwd) if cwd is not None else None,
            env=dict(environment) if environment is not None else None,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            timeout=timeout,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        die("command failed to run {!r}: {}".format(list(argv), exc))
    if result.returncode != 0:
        die("command failed exit={} argv={!r} stderr={!r}".format(
            result.returncode, list(argv), result.stderr[-2000:]))
    return result.stdout, result.stderr


def tool(name: str) -> Path:
    value = shutil.which(name)
    if value is None:
        die("required tool is unavailable: {}".format(name))
    path = Path(value).resolve()
    if not path.is_file() or not os.access(path, os.X_OK):
        die("required tool is not executable: {}".format(path))
    return path


def parse_uint(text: object, name: str, maximum: int = (1 << 64) - 1) -> int:
    if not isinstance(text, str) or UINT_RE.fullmatch(text) is None:
        die("{} is not a canonical uint".format(name))
    value = int(text)
    if value > maximum:
        die("{} exceeds its integer domain".format(name))
    return value


def parse_milli(text: object, name: str) -> int:
    if not isinstance(text, str):
        die("{} is not a fixed integer mean".format(name))
    match = MILLI_RE.fullmatch(text)
    if match is None:
        die("{} is not an exact one-trial integer mean".format(name))
    return int(match.group(1))


def process_start_ticks(pid: int) -> int:
    try:
        raw = Path("/proc/{}/stat".format(pid)).read_bytes()
        return int(raw.rsplit(b") ", 1)[1].split()[19])
    except (OSError, IndexError, ValueError) as exc:
        die("cannot receipt process {} start time: {}".format(pid, exc))


def process_cmdline(pid: int) -> Tuple[str, ...]:
    try:
        raw = Path("/proc/{}/cmdline".format(pid)).read_bytes()
        return tuple(token.decode("utf-8") for token in raw.split(b"\0") if token)
    except (OSError, UnicodeDecodeError) as exc:
        die("cannot receipt process {} command line: {}".format(pid, exc))


def filler_snapshot() -> Dict[str, object]:
    if FILLER_PID_FILE.is_symlink() or not FILLER_PID_FILE.is_file():
        die("normal filler PID file is unavailable")
    text = FILLER_PID_FILE.read_text(encoding="ascii")
    if re.fullmatch(r"[1-9][0-9]*\n", text) is None:
        die("normal filler PID file is malformed")
    parent = int(text)
    parent_tokens = process_cmdline(parent)
    if not any(Path(token).name == "wirehair_load_fillers.sh"
               for token in parent_tokens):
        die("filler PID does not identify the normal filler controller")
    workers = []
    for status_path in Path("/proc").glob("[0-9]*/status"):
        pid = int(status_path.parent.name)
        try:
            status = status_path.read_text(encoding="ascii")
            ppid_match = re.search(r"^PPid:\s*([0-9]+)$", status, re.MULTILINE)
            nice_match = re.search(r"^NSpid:\s*(.+)$", status, re.MULTILINE)
            del nice_match  # Namespace depth is irrelevant; keep parsing bounded.
            if ppid_match is None or int(ppid_match.group(1)) != parent:
                continue
            tokens = process_cmdline(pid)
            if tokens[-3:] != ("-c", "while :; do :; done") and not any(
                    "while :; do :; done" in token for token in tokens):
                continue
            affinity = sorted(os.sched_getaffinity(pid))
            priority = os.getpriority(os.PRIO_PROCESS, pid)
        except (OSError, IndexError, CampaignError):
            continue
        if len(affinity) != 1 or priority < 19:
            die("normal filler worker is not singleton-pinned at nice 19")
        workers.append({
            "pid": pid, "start_ticks": process_start_ticks(pid),
            "cpu": affinity[0], "nice": priority,
            "cmdline_sha256": sha256_bytes(b"\0".join(
                token.encode("utf-8") for token in tokens) + b"\0"),
        })
    workers.sort(key=lambda value: int(value["cpu"]))
    cpus = [int(value["cpu"]) for value in workers]
    if len(workers) < MIN_FILLER_COUNT or len(cpus) != len(set(cpus)):
        die("normal filler load is incomplete or duplicates CPUs")
    return {
        "pid_file": str(FILLER_PID_FILE), "parent_pid": parent,
        "parent_start_ticks": process_start_ticks(parent),
        "parent_cmdline_sha256": sha256_bytes(b"\0".join(
            token.encode("utf-8") for token in parent_tokens) + b"\0"),
        "worker_count": len(workers), "cpus": cpus, "workers": workers,
    }


def sampler_snapshot() -> List[Dict[str, object]]:
    records = []
    for path in Path("/proc").glob("[0-9]*/cmdline"):
        pid = int(path.parent.name)
        try:
            tokens = process_cmdline(pid)
            matching = [token for token in tokens
                        if Path(token).name == SAMPLER_BASENAME]
            if len(matching) != 1:
                continue
            affinity = sorted(os.sched_getaffinity(pid))
        except (OSError, CampaignError):
            continue
        records.append({
            "pid": pid, "start_ticks": process_start_ticks(pid),
            "affinity": affinity,
            "cmdline_sha256": sha256_bytes(b"\0".join(
                token.encode("utf-8") for token in tokens) + b"\0"),
        })
    records.sort(key=lambda value: int(value["pid"]))
    if not records:
        die("the externally managed thermal sampler is not running")
    return records


def validate_load_receipts(filler: object, sampler: object) -> None:
    filler_fields = {
        "pid_file", "parent_pid", "parent_start_ticks",
        "parent_cmdline_sha256", "worker_count", "cpus", "workers",
    }
    worker_fields = {"pid", "start_ticks", "cpu", "nice", "cmdline_sha256"}
    if not isinstance(filler, dict) or set(filler) != filler_fields:
        die("filler receipt schema is malformed")
    workers = filler.get("workers")
    cpus = filler.get("cpus")
    count = filler.get("worker_count")
    if (filler.get("pid_file") != str(FILLER_PID_FILE) or
            not isinstance(workers, list) or not isinstance(cpus, list) or
            not isinstance(count, int) or isinstance(count, bool) or
            count < MIN_FILLER_COUNT or len(workers) != count or
            len(cpus) != count or cpus != sorted(set(cpus)) or
            not isinstance(filler.get("parent_pid"), int) or
            filler["parent_pid"] <= 0 or
            not isinstance(filler.get("parent_start_ticks"), int) or
            filler["parent_start_ticks"] < 0 or
            not isinstance(filler.get("parent_cmdline_sha256"), str) or
            SHA_RE.fullmatch(filler["parent_cmdline_sha256"]) is None):
        die("filler receipt values are malformed")
    for index, worker in enumerate(workers):
        if (not isinstance(worker, dict) or set(worker) != worker_fields or
                worker.get("cpu") != cpus[index] or worker.get("nice") != 19 or
                not isinstance(worker.get("pid"), int) or worker["pid"] <= 0 or
                not isinstance(worker.get("start_ticks"), int) or
                worker["start_ticks"] < 0 or
                not isinstance(worker.get("cmdline_sha256"), str) or
                SHA_RE.fullmatch(worker["cmdline_sha256"]) is None):
            die("filler worker receipt is malformed")
    sampler_fields = {"pid", "start_ticks", "affinity", "cmdline_sha256"}
    if not isinstance(sampler, list) or not sampler:
        die("external sampler receipt is malformed")
    previous_pid = 0
    for record in sampler:
        if (not isinstance(record, dict) or set(record) != sampler_fields or
                not isinstance(record.get("pid"), int) or
                record["pid"] <= previous_pid or
                not isinstance(record.get("start_ticks"), int) or
                record["start_ticks"] < 0 or
                not isinstance(record.get("affinity"), list) or
                not record["affinity"] or record["affinity"] !=
                    sorted(set(record["affinity"])) or
                not isinstance(record.get("cmdline_sha256"), str) or
                SHA_RE.fullmatch(record["cmdline_sha256"]) is None):
            die("external sampler process receipt is malformed")
        previous_pid = record["pid"]


def arm_preamble(arm: str, seed_index: int, schedule: str) -> Dict[str, str]:
    if arm not in ARMS or not 0 <= seed_index < len(SEEDS) or \
            schedule not in SCHEDULES:
        die("invalid preamble coordinate")
    return {
        "trials": "1", "threads": "1", "loss": LOSS,
        "seed": SEEDS[seed_index], "completion": "mixed",
        "mixed_period": "244", "mixed_gf256_rows": "10",
        "mixed_gf16_rows": "2", "mixed_geometry": "frozen",
        "mixed_residue_skew": "0", "mixed_residue_schedule": "constant",
        "mixed_residue_hash_seed": "0x0", "mixed_residue_hash_keyed": "0",
        "mixed_independent_extension_residues": "0",
        "mixed_grouped_gf256_rows": "0",
        "mixed_grouped_gf256_hash_seed": "0x0",
        "mixed_grouped_final_h_a_columns": "0",
        "mixed_residue_buckets_requested": "auto",
        "mixed_extension_residue_seed_xor": "0x4e",
        "source_hits_override": "0", "packet_peel_seed_xor": "0x0",
        "packet_peel_seed_table": "none", "binary_dense_rows_override": "0",
        "binary_dense_two_anchor": "1" if arm == "baseline" else "0",
        "binary_dense_two_anchor_phase": "0",
        "binary_dense_segmented_anchors": "none" if arm == "baseline" else arm,
        "gf256_heavy_rows_override": "0", "odd_packet_peel_seed_xor": "0x0",
        "packet_row_seed_multiplier": "0x1",
        "packet_row_seed_avalanche": "0", "seed_block_bytes_override": "0",
        "overhead_stream": "salted", "full_payload_solve": "0",
        "schedule": schedule,
    }


def parse_preamble(line: str) -> Dict[str, str]:
    prefix = "# precodefail: "
    if not line.startswith(prefix):
        die("missing precodefail preamble")
    result: Dict[str, str] = {}
    for token in line[len(prefix):].split(" "):
        if token.count("=") != 1:
            die("malformed precodefail preamble token")
        key, value = token.split("=", 1)
        if key in result or not key or not value:
            die("duplicate or empty precodefail preamble token")
        result[key] = value
    return result


def task_argv(
    binary: str, arm: str, seed_index: int, schedule: str,
    Ks: Sequence[int],
) -> List[str]:
    if arm not in ARMS or schedule not in SCHEDULES or not Ks:
        die("invalid task command coordinate")
    return [
        binary, "precodefail", "--bb-list", "64", "--overhead", "0",
        "--trials", "1", "--threads", "1", "--loss", "0.50",
        "--completion", "mixed", "--mix-count", "2",
        "--packet-peel-seed-xor", "0", "--N", ",".join(map(str, Ks)),
        "--seed", SEEDS[seed_index], "--schedule", schedule,
        "--mixed-period", "244", "--mixed-geometry", "frozen",
        "--mixed-gf256-rows", "10", "--mixed-gf16-rows", "2",
        "--mixed-grouped-gf256-rows", "0", "--mixed-residue-buckets", "auto",
        *ARM_OPTIONS[arm],
    ]


def build_tasks(
    binary: str, *, k_lo: int = K_LO, k_hi: int = K_HI,
    chunk_max: int = CHUNK_MAX,
) -> Tuple[Dict[str, object], ...]:
    if not K_LO <= k_lo <= k_hi <= K_HI or not 1 <= chunk_max <= CHUNK_MAX:
        die("invalid task-grid bounds")
    tasks: List[Dict[str, object]] = []
    for seed_index, seed in enumerate(SEEDS):
        for schedule in SCHEDULES:
            for chunk_index, lo in enumerate(range(k_lo, k_hi + 1, chunk_max)):
                Ks = list(range(lo, min(lo + chunk_max - 1, k_hi) + 1))
                for arm in ARMS:
                    job = len(tasks)
                    name = "{:05d}.{}.seed{}.{}.chunk{:03d}.csv".format(
                        job, arm, seed_index, schedule, chunk_index)
                    task: Dict[str, object] = {
                        "job": job, "arm": arm, "seed_index": seed_index,
                        "seed": seed, "schedule": schedule,
                        "chunk_index": chunk_index, "K_lo": Ks[0],
                        "K_hi": Ks[-1], "Ks": Ks, "output_name": name,
                        "stdout_max_bytes": MAX_STDOUT_BYTES,
                        "stderr_max_bytes": MAX_STDERR_BYTES,
                        "timeout_seconds": TASK_TIMEOUT_SECONDS,
                    }
                    task["argv"] = task_argv(
                        binary, arm, seed_index, schedule, Ks)
                    tasks.append(task)
    return tuple(tasks)


def audit_tasks(
    tasks: Sequence[Mapping[str, object]], *, k_lo: int = K_LO,
    k_hi: int = K_HI,
) -> Dict[str, object]:
    expected_Ks = list(range(k_lo, k_hi + 1))
    strata = {(arm, seed_index, schedule): []
              for arm in ARMS for seed_index in range(len(SEEDS))
              for schedule in SCHEDULES}
    names = set()
    for expected_job, task in enumerate(tasks):
        arm = str(task.get("arm"))
        seed_index = task.get("seed_index")
        schedule = str(task.get("schedule"))
        Ks = task.get("Ks")
        if (task.get("job") != expected_job or arm not in ARMS or
                not isinstance(seed_index, int) or isinstance(seed_index, bool) or
                not 0 <= seed_index < len(SEEDS) or schedule not in SCHEDULES or
                task.get("seed") != SEEDS[seed_index] or
                not isinstance(Ks, list) or not Ks or
                any(not isinstance(K, int) or isinstance(K, bool) for K in Ks) or
                Ks != list(range(Ks[0], Ks[-1] + 1)) or len(Ks) > CHUNK_MAX or
                task.get("K_lo") != Ks[0] or task.get("K_hi") != Ks[-1] or
                task.get("stdout_max_bytes") != MAX_STDOUT_BYTES or
                task.get("stderr_max_bytes") != MAX_STDERR_BYTES or
                task.get("timeout_seconds") != TASK_TIMEOUT_SECONDS):
            die("task semantic ledger mismatch at job {}".format(expected_job))
        argv = task.get("argv")
        if not isinstance(argv, list) or not argv or argv != task_argv(
                str(argv[0]), arm, seed_index, schedule, Ks):
            die("task command ledger mismatch at job {}".format(expected_job))
        name = task.get("output_name")
        if not isinstance(name, str) or "/" in name or name in names:
            die("task output ledger is unsafe or duplicate")
        names.add(name)
        strata[(arm, seed_index, schedule)].extend(Ks)
    for offset in range(0, len(tasks), len(ARMS)):
        group = tasks[offset:offset + len(ARMS)]
        if tuple(str(task["arm"]) for task in group) != ARMS or any(
                group[0][key] != task[key] for task in group[1:]
                for key in ("seed_index", "seed", "schedule", "chunk_index",
                            "K_lo", "K_hi", "Ks")):
            die("task architecture group is not exactly paired")
    if any(values != expected_Ks for values in strata.values()):
        die("task ledger does not cover every all-K stratum exactly")
    cells = len(expected_Ks) * len(SEEDS) * len(SCHEDULES)
    return {
        "audit": "OK", "K_range": [k_lo, k_hi],
        "K_count": len(expected_Ks), "arms": list(ARMS),
        "seeds": list(SEEDS), "schedules": list(SCHEDULES),
        "chunk_max": CHUNK_MAX, "task_count": len(tasks),
        "task_groups": len(tasks) // len(ARMS),
        "cells_per_arm": cells, "total_cells": cells * len(ARMS),
    }


def validate_row(row: Mapping[str, str]) -> Dict[str, object]:
    for field in CSV_FIELDS:
        if field not in row or not isinstance(row[field], str):
            die("benchmark row is missing {}".format(field))
    K = parse_uint(row["N"], "K", K_HI)
    if K < K_LO:
        die("benchmark K is below the codec domain")
    fixed = {
        "bb": "64", "heavy_family": "periodic", "mix_count": "2",
        "overhead": "0", "trials": "1",
        "active_packet_peel_seed_xor": "0x0",
    }
    if any(row[key] != value for key, value in fixed.items()):
        die("benchmark row changed a fixed architecture field")
    success = parse_uint(row["success"], "success", 1)
    rank_fail = parse_uint(row["rank_fail"], "rank_fail", 1)
    error = parse_uint(row["error"], "error", 1)
    if success + rank_fail + error != 1 or error:
        die("benchmark one-trial outcome is inconsistent")
    expected_rate = "0.00000000" if success else "1.00000000"
    if row["fail_rate"] != expected_rate:
        die("benchmark failure rate is inconsistent")
    inactivated = parse_milli(row["inact_mu"], "inactivation mean")
    q = parse_milli(row["binary_def_mu"], "binary deficit mean")
    heavy_gain = parse_milli(row["heavy_gain_mu"], "heavy gain mean")
    xors = parse_milli(row["block_xors_mu"], "block XOR mean")
    muladds = parse_milli(row["block_muladds_mu"], "block muladd mean")
    heavy_shortfall = parse_uint(row["heavy_shortfall"], "heavy shortfall", 1)
    expected_shortfall = int(
        rank_fail == 1 and q <= 12 and heavy_gain < q)
    if (parse_uint(row["inact_max"], "inactivation max") != inactivated or
            parse_uint(row["binary_def_max"], "binary deficit max") != q or
            parse_uint(row["heavy_gain_min"], "heavy gain min") != heavy_gain or
            heavy_shortfall != expected_shortfall or
            row["binary_def_hist"] != "{}:1".format(q) or
            row["heavy_gain_hist"] != "{}:1".format(heavy_gain)):
        die("benchmark deterministic work receipt is inconsistent")
    seed_attempt = parse_uint(row["seed_attempt"], "seed attempt", 255)
    if success:
        if (row["first_rank_fail"] != "-1" or row["failure_trials"] != ""):
            die("successful benchmark row has a failure receipt")
    else:
        if row["first_rank_fail"] != "0" or row["failure_trials"] != "0":
            die("failed benchmark row lacks the one-trial failure receipt")
    for field in TIMING_FIELDS:
        try:
            value = Decimal(row[field])
        except InvalidOperation as exc:
            die("benchmark timing field is malformed: {}".format(field))
        if not value.is_finite() or value < 0:
            die("benchmark timing field is negative or nonfinite")
    for field in ZERO_JOINT_FIELDS:
        if row[field] != "0.000":
            die("grouped/dual work unexpectedly entered the q0 campaign")
    outcome = (
        "success" if success else "q>H" if q > 12 else
        "field_shortfall" if heavy_shortfall else "other_rank_fail")
    return {
        "K": K, "outcome": outcome, "inactivated": inactivated,
        "q": q, "heavy_gain": heavy_gain, "seed_attempt": seed_attempt,
        "xors": xors, "muladds": muladds,
    }


def parse_output(
    raw: bytes, task: Mapping[str, object],
) -> Tuple[Dict[str, str], ...]:
    maximum = task.get("stdout_max_bytes")
    if not isinstance(maximum, int) or len(raw) > maximum:
        die("benchmark stdout exceeded its sealed bound")
    if (not raw or not raw.endswith(b"\n") or b"\r" in raw or b"\0" in raw or
            b'"' in raw):
        die("benchmark stdout is not canonical LF text")
    try:
        lines = raw.decode("ascii").splitlines()
    except UnicodeDecodeError as exc:
        die("benchmark stdout is not ASCII: {}".format(exc))
    if len(lines) != len(task["Ks"]) + 2:
        die("benchmark stdout row count changed")
    preamble = parse_preamble(lines[0])
    expected = arm_preamble(
        str(task["arm"]), int(task["seed_index"]), str(task["schedule"]))
    if preamble != expected:
        die("benchmark preamble differs from the raw architecture contract")
    reader = csv.DictReader(io.StringIO("\n".join(lines[1:]) + "\n"))
    if tuple(reader.fieldnames or ()) != CSV_FIELDS:
        die("benchmark CSV header changed")
    if any(line.count(",") != len(CSV_FIELDS) - 1 for line in lines[1:]):
        die("benchmark CSV field count changed")
    rows = tuple(dict(row) for row in reader)
    if any(set(row) != set(CSV_FIELDS) for row in rows):
        die("benchmark CSV has missing or extra fields")
    if [parse_uint(row["N"], "row K", K_HI) for row in rows] != task["Ks"]:
        die("benchmark K ledger differs from its task")
    for row in rows:
        validate_row(row)
    return rows


def git_text(git: Path, repo: Path, *args: str) -> str:
    stdout, stderr = checked((str(git), "-C", str(repo), *args), timeout=120)
    if stderr:
        die("git command produced stderr")
    return stdout.decode("utf-8").strip()


def prepare(args: argparse.Namespace) -> None:
    root = Path(args.root).resolve()
    repo = Path(args.repo).resolve()
    if root.exists():
        die("campaign root already exists")
    if not repo.is_dir() or args.workers != WORKERS:
        die("prepare requires the frozen repository and 56 workers")
    tools = {name: tool(name) for name in (
        "git", "cmake", "ninja", "env", "python3", "c++")}
    if Path(sys.executable).resolve() != tools["python3"]:
        die("prepare is not using the receipted Python interpreter")
    resolved = git_text(tools["git"], repo, "rev-parse", PROTOTYPE_COMMIT + "^{commit}")
    parent = git_text(tools["git"], repo, "rev-parse", PROTOTYPE_COMMIT + "^")
    if resolved != PROTOTYPE_COMMIT or parent != PROTOTYPE_PARENT:
        die("segmented-anchor prototype identity or parent changed")
    paths = tuple(git_text(
        tools["git"], repo, "diff", "--name-only", PROTOTYPE_PARENT,
        PROTOTYPE_COMMIT).splitlines())
    if paths != EXPECTED_PROTOTYPE_PATHS:
        die("segmented-anchor prototype changed-path audit failed")
    diff_stdout, diff_stderr = checked((
        str(tools["git"]), "-C", str(repo), "diff", "--binary",
        PROTOTYPE_PARENT, PROTOTYPE_COMMIT), timeout=120)
    if diff_stderr:
        die("prototype diff produced stderr")

    staging = root.with_name(root.name + ".prepare.{}".format(os.getpid()))
    workspace = root.with_name(root.name + ".build.{}".format(os.getpid()))
    if staging.exists() or workspace.exists():
        die("stale prepare workspace exists")
    staging.mkdir(parents=True)
    workspace.mkdir(parents=True)
    try:
        frozen = staging / "frozen"
        provenance = staging / "provenance"
        for directory in (
                frozen, provenance, staging / "raw", staging / "stderr",
                staging / "exit", staging / "receipts"):
            directory.mkdir()
        source = workspace / "source"
        build = workspace / "build"
        home = workspace / "home"
        add_out, add_err = checked((
            str(tools["git"]), "-C", str(repo), "worktree", "add", "--detach",
            str(source), PROTOTYPE_COMMIT), timeout=120)
        atomic_new(provenance / "worktree.stdout", add_out)
        atomic_new(provenance / "worktree.stderr", add_err)
        try:
            home.mkdir()
            environment = {
                "HOME": str(home), "PATH": "/usr/bin:/bin", "LC_ALL": "C",
                "LANG": "C", "TZ": "UTC",
            }
            configure_argv = (
                str(tools["cmake"]), "-S", str(source), "-B", str(build),
                "-G", "Ninja", "-DCMAKE_BUILD_TYPE=Release",
                "-DBUILD_TESTS=ON", "-DBUILD_CODEC_V2=ON",
                "-DWIREHAIR_BUILD_BENCHMARKS=ON", "-DMARCH_NATIVE=ON",
                "-DWIREHAIR_STRICT_WARNINGS=ON", "-DWH_LTO=OFF",
                "-DWH_PGO_MODE=OFF", "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
            )
            configure_out, configure_err = checked(
                configure_argv, environment=environment, timeout=600)
            atomic_new(provenance / "configure.stdout", configure_out)
            atomic_new(provenance / "configure.stderr", configure_err)
            build_argv = (
                str(tools["cmake"]), "--build", str(build), "--target",
                "wirehair_v2_bench", "--parallel", str(WORKERS), "--verbose",
            )
            build_out, build_err = checked(
                build_argv, environment=environment, timeout=1800)
            atomic_new(provenance / "build.stdout", build_out)
            atomic_new(provenance / "build.stderr", build_err)
            binary_source = build / "codec/wirehair_v2_bench"
            if not binary_source.is_file() or not os.access(binary_source, os.X_OK):
                die("fresh prototype build did not produce the benchmark binary")
            binary = frozen / BINARY_NAME
            shutil.copyfile(binary_source, binary)
            os.chmod(binary, 0o555)
            for path, name in (
                    (build / "CMakeCache.txt", "CMakeCache.txt"),
                    (build / "compile_commands.json", "compile_commands.json")):
                if not path.is_file():
                    die("fresh build provenance is incomplete")
                shutil.copyfile(path, provenance / name)
                os.chmod(provenance / name, 0o444)
        finally:
            removal = subprocess.run((
                str(tools["git"]), "-C", str(repo), "worktree", "remove",
                "--force", str(source)), stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, check=False)
            if removal.returncode != 0 and source.exists():
                die("temporary build worktree could not be removed")

        controller_source = Path(__file__).resolve()
        controller = frozen / CONTROLLER_NAME
        shutil.copyfile(controller_source, controller)
        os.chmod(controller, 0o555)
        binary = frozen / BINARY_NAME
        smoke_records = []
        smoke_Ks = [2, 257, 4096, 64_000]
        for schedule in SCHEDULES:
            for arm in ARMS:
                task = {
                    "arm": arm, "seed_index": 0, "schedule": schedule,
                    "Ks": smoke_Ks, "stdout_max_bytes": MAX_STDOUT_BYTES,
                }
                argv = task_argv(str(binary), arm, 0, schedule, smoke_Ks)
                stdout, stderr = checked(argv, cwd=staging, timeout=300)
                if stderr:
                    die("bounded segmented-anchor smoke produced stderr")
                rows = parse_output(stdout, task)
                name = "smoke.{}.{}.csv".format(schedule, arm)
                atomic_new(provenance / name, stdout)
                smoke_records.append({
                    "arm": arm, "schedule": schedule, "seed_index": 0,
                    "Ks": smoke_Ks, "argv": ["frozen/" + BINARY_NAME, *argv[1:]],
                    "path": "provenance/" + name,
                    "sha256": sha256_bytes(stdout), "row_count": len(rows),
                })

        final_binary = root / "frozen" / BINARY_NAME
        tasks = build_tasks(str(final_binary))
        audit = audit_tasks(tasks)
        manifest = b"".join(canonical_json(task) for task in tasks)
        atomic_new(staging / "tasks_manifest.jsonl", manifest)
        immutable = {}
        for directory in (frozen, provenance):
            for path in sorted(directory.iterdir()):
                if not path.is_file() or path.is_symlink():
                    die("immutable prepare inventory contains a non-regular file")
                immutable[str(path.relative_to(staging))] = sha256_file(path)
        tool_records = {
            name: {"path": str(path), "sha256": sha256_file(path)}
            for name, path in sorted(tools.items())
        }
        design = sealed(SCHEMA + ".design", {
            "root": str(root), "prototype_commit": PROTOTYPE_COMMIT,
            "prototype_parent": PROTOTYPE_PARENT,
            "prototype_diff_sha256": sha256_bytes(diff_stdout),
            "prototype_changed_paths": list(paths),
            "binary_sha256": sha256_file(binary),
            "controller_sha256": sha256_file(controller),
            "architecture": {
                "baseline": "D12 two-anchor phase0 q0",
                "candidates": ["three-048", "three-059", "four-0369"],
                "dense_rows": 12, "gf256_rows": 10, "gf16_rows": 2,
                "period": 244, "mix_count": 2, "overhead": 0,
                "loss": LOSS, "block_bytes": 64,
                "seed_fixes_applied": False,
                "packet_peel_seed_xor": 0,
                "packet_peel_seed_table": "none",
                "segmented_layout_is_only_architecture_selector": True,
            },
            "K_range": [K_LO, K_HI], "seeds": list(SEEDS),
            "schedules": list(SCHEDULES), "arms": list(ARMS),
            "workers": WORKERS, "task_audit": audit,
            "tasks_manifest_sha256": sha256_bytes(manifest),
            "smoke": smoke_records, "immutable_files": immutable,
            "tools": tool_records,
            "execution_environment": {
                "evidentiary_scope": "deterministic raw recovery and work; not timing",
                "normal_nice19_fillers_required": True,
                "minimum_filler_workers": MIN_FILLER_COUNT,
                "external_thermal_sampler_required": True,
                "controller_starts_sampler": False,
                "controller_signals_sampler": False,
                "quiet_timing": False,
            },
            "python": {
                "executable": str(Path(sys.executable).resolve()),
                "version": sys.version,
            },
        })
        atomic_new(staging / "design.json", canonical_json(design))
        receipt = sealed(SCHEMA + ".prepare", {
            "prepared_utc": utc_now(),
            "design_sha256": sha256_file(staging / "design.json"),
            "tasks_manifest_sha256": sha256_bytes(manifest),
            "binary_sha256": sha256_file(binary),
            "controller_sha256": sha256_file(controller),
            "immutable_files": immutable,
        })
        atomic_new(staging / "prepare_receipt.json", canonical_json(receipt))
        for directory in (frozen, provenance):
            os.chmod(directory, 0o555)
        os.replace(staging, root)
        print(json.dumps({
            "prepared": True, "root": str(root),
            "design_sha256": sha256_file(root / "design.json"),
            "prepare_receipt_sha256": sha256_file(root / "prepare_receipt.json"),
            "binary_sha256": sha256_file(root / "frozen" / BINARY_NAME),
            **audit,
        }, sort_keys=True))
    finally:
        if staging.exists():
            shutil.rmtree(staging)
        if workspace.exists():
            shutil.rmtree(workspace)


def load_prepared(root: Path, *, require_fresh: bool = False) -> Tuple[
        Dict[str, object], Tuple[Dict[str, object], ...]]:
    root = root.resolve()
    design = load_canonical(root / "design.json", "campaign design")
    verify_sealed(design, SCHEMA + ".design", "campaign design")
    if (design.get("root") != str(root) or
            design.get("prototype_commit") != PROTOTYPE_COMMIT or
            design.get("prototype_parent") != PROTOTYPE_PARENT or
            design.get("prototype_changed_paths") != list(EXPECTED_PROTOTYPE_PATHS) or
            design.get("K_range") != [K_LO, K_HI] or
            design.get("seeds") != list(SEEDS) or
            design.get("schedules") != list(SCHEDULES) or
            design.get("arms") != list(ARMS) or design.get("workers") != WORKERS):
        die("campaign design differs from the frozen experiment")
    expected_architecture = {
        "baseline": "D12 two-anchor phase0 q0",
        "candidates": ["three-048", "three-059", "four-0369"],
        "dense_rows": 12, "gf256_rows": 10, "gf16_rows": 2,
        "period": 244, "mix_count": 2, "overhead": 0,
        "loss": LOSS, "block_bytes": 64, "seed_fixes_applied": False,
        "packet_peel_seed_xor": 0, "packet_peel_seed_table": "none",
        "segmented_layout_is_only_architecture_selector": True,
    }
    expected_environment = {
        "evidentiary_scope": "deterministic raw recovery and work; not timing",
        "normal_nice19_fillers_required": True,
        "minimum_filler_workers": MIN_FILLER_COUNT,
        "external_thermal_sampler_required": True,
        "controller_starts_sampler": False,
        "controller_signals_sampler": False,
        "quiet_timing": False,
    }
    if (design.get("architecture") != expected_architecture or
            design.get("execution_environment") != expected_environment):
        die("campaign architecture or environment policy changed")
    python = design.get("python")
    if python != {
            "executable": str(Path(sys.executable).resolve()),
            "version": sys.version}:
        die("campaign Python runtime changed")
    immutable = design.get("immutable_files")
    if not isinstance(immutable, dict) or not immutable:
        die("immutable file ledger is malformed")
    actual_immutable = set()
    for dirname in ("frozen", "provenance"):
        directory = root / dirname
        if directory.is_symlink() or not directory.is_dir():
            die("immutable campaign directory is unavailable")
        for path in directory.iterdir():
            if path.is_symlink() or not path.is_file():
                die("immutable campaign inventory contains a non-regular file")
            actual_immutable.add(str(path.relative_to(root)))
    if actual_immutable != set(immutable):
        die("immutable campaign inventory changed")
    for relative, digest in immutable.items():
        if (not isinstance(relative, str) or not isinstance(digest, str) or
                SHA_RE.fullmatch(digest) is None or
                sha256_file(root / relative) != digest):
            die("immutable campaign file changed: {}".format(relative))
    controller = root / "frozen" / CONTROLLER_NAME
    binary = root / "frozen" / BINARY_NAME
    if (design.get("controller_sha256") != sha256_file(controller) or
            design.get("binary_sha256") != sha256_file(binary) or
            Path(__file__).resolve().read_bytes() != controller.read_bytes()):
        die("active controller or frozen binary identity changed")
    tools = design.get("tools")
    if not isinstance(tools, dict):
        die("tool identity ledger is malformed")
    for name, record in tools.items():
        if (not isinstance(record, dict) or set(record) != {"path", "sha256"} or
                not isinstance(record.get("path"), str) or
                not isinstance(record.get("sha256"), str) or
                sha256_file(Path(record["path"])) != record["sha256"]):
            die("campaign tool changed: {}".format(name))
    raw_manifest = (root / "tasks_manifest.jsonl").read_bytes()
    if sha256_bytes(raw_manifest) != design.get("tasks_manifest_sha256"):
        die("task manifest hash changed")
    tasks: List[Dict[str, object]] = []
    for line in raw_manifest.splitlines(keepends=True):
        try:
            value = json.loads(line.decode("ascii"))
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            die("task manifest line is malformed: {}".format(exc))
        if not isinstance(value, dict) or canonical_json(value) != line:
            die("task manifest line is not canonical")
        tasks.append(value)
    expected_tasks = build_tasks(str(binary))
    if tuple(tasks) != expected_tasks:
        die("task manifest differs from the generated all-K cube")
    audit = audit_tasks(tasks)
    if design.get("task_audit") != audit:
        die("task audit receipt does not replay")
    smoke = design.get("smoke")
    if not isinstance(smoke, list) or len(smoke) != len(ARMS) * len(SCHEDULES):
        die("bounded smoke ledger is malformed")
    smoke_Ks = [2, 257, 4096, 64_000]
    expected_smoke = []
    for schedule in SCHEDULES:
        for arm in ARMS:
            name = "smoke.{}.{}.csv".format(schedule, arm)
            task = {
                "arm": arm, "seed_index": 0, "schedule": schedule,
                "Ks": smoke_Ks, "stdout_max_bytes": MAX_STDOUT_BYTES,
            }
            raw = (root / "provenance" / name).read_bytes()
            rows = parse_output(raw, task)
            argv = task_argv(str(binary), arm, 0, schedule, smoke_Ks)
            expected_smoke.append({
                "arm": arm, "schedule": schedule, "seed_index": 0,
                "Ks": smoke_Ks, "argv": ["frozen/" + BINARY_NAME, *argv[1:]],
                "path": "provenance/" + name,
                "sha256": sha256_bytes(raw), "row_count": len(rows),
            })
    if smoke != expected_smoke:
        die("bounded smoke ledger does not replay")
    prepare_receipt = load_canonical(
        root / "prepare_receipt.json", "prepare receipt")
    verify_sealed(prepare_receipt, SCHEMA + ".prepare", "prepare receipt")
    expected_prepare = {
        "design_sha256": sha256_file(root / "design.json"),
        "tasks_manifest_sha256": sha256_bytes(raw_manifest),
        "binary_sha256": sha256_file(binary),
        "controller_sha256": sha256_file(controller),
        "immutable_files": immutable,
    }
    if any(prepare_receipt.get(key) != value
           for key, value in expected_prepare.items()):
        die("prepare receipt does not bind the campaign")
    if not isinstance(prepare_receipt.get("prepared_utc"), str):
        die("prepare receipt timestamp is malformed")
    expected_directories = {"frozen", "provenance", "raw", "stderr", "exit", "receipts"}
    actual_directories = {
        str(path.relative_to(root)) for path in root.rglob("*")
        if path.is_dir() and not path.is_symlink()}
    if actual_directories != expected_directories or any(
            path.is_symlink() for path in root.rglob("*")):
        die("campaign directory inventory changed")
    if require_fresh:
        for dirname in ("raw", "stderr", "exit", "receipts"):
            if any((root / dirname).iterdir()):
                die("fresh launch found prior output in {}".format(dirname))
        for name in ("launch_receipt.json", "data_manifest.sha256",
                     "validated_summary.json"):
            if (root / name).exists():
                die("fresh launch found prior {}".format(name))
    return design, tuple(tasks)


def verify_prepared(args: argparse.Namespace) -> None:
    root = Path(args.root).resolve()
    design, tasks = load_prepared(root, require_fresh=True)
    print(json.dumps({
        "verified_prepared": True, "root": str(root),
        "design_sha256": sha256_file(root / "design.json"),
        "prepare_receipt_sha256": sha256_file(root / "prepare_receipt.json"),
        "binary_sha256": design["binary_sha256"], "task_count": len(tasks),
        "total_cells": design["task_audit"]["total_cells"],
    }, sort_keys=True))


def run_one(root: Path, task: Mapping[str, object]) -> Dict[str, object]:
    started_utc = utc_now()
    start_ns = time.monotonic_ns()
    try:
        result = subprocess.run(
            list(task["argv"]), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            check=False, timeout=int(task["timeout_seconds"]),
            start_new_session=True,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        die("task {} failed to run: {}".format(task["job"], exc))
    end_ns = time.monotonic_ns()
    if result.returncode != 0:
        die("task {} exited {} stderr={!r}".format(
            task["job"], result.returncode, result.stderr[-2000:]))
    if result.stderr or len(result.stderr) > int(task["stderr_max_bytes"]):
        die("task {} produced stderr".format(task["job"]))
    rows = parse_output(result.stdout, task)
    name = str(task["output_name"])
    raw_path = root / "raw" / name
    stderr_path = root / "stderr" / (name + ".stderr")
    exit_path = root / "exit" / (name + ".exit")
    atomic_new(raw_path, result.stdout)
    atomic_new(stderr_path, result.stderr)
    atomic_new(exit_path, b"0\n")
    receipt = sealed(SCHEMA + ".job", {
        "job": task["job"], "task_sha256": sha256_bytes(canonical_json(task)),
        "argv": task["argv"], "output_name": name,
        "started_utc": started_utc, "start_monotonic_ns": start_ns,
        "end_monotonic_ns": end_ns, "duration_ns": end_ns - start_ns,
        "returncode": 0, "row_count": len(rows),
        "stdout_sha256": sha256_bytes(result.stdout),
        "stdout_bytes": len(result.stdout),
        "stderr_sha256": sha256_bytes(result.stderr),
        "stderr_bytes": len(result.stderr),
        "exit_sha256": sha256_bytes(b"0\n"),
    })
    receipt_path = root / "receipts" / (name + ".json")
    atomic_new(receipt_path, canonical_json(receipt))
    return {
        "job": task["job"], "name": name,
        "receipt_sha256": sha256_file(receipt_path),
    }


def run_campaign(args: argparse.Namespace) -> None:
    root = Path(args.root).resolve()
    design, tasks = load_prepared(root, require_fresh=True)
    if args.workers != WORKERS:
        die("launch worker count differs from the frozen campaign")
    filler_start = filler_snapshot()
    sampler_start = sampler_snapshot()
    validate_load_receipts(filler_start, sampler_start)
    started_utc = utc_now()
    start_ns = time.monotonic_ns()
    records: List[Dict[str, object]] = []
    iterator = iter(tasks)
    with concurrent.futures.ThreadPoolExecutor(max_workers=WORKERS) as pool:
        pending: Dict[concurrent.futures.Future[Dict[str, object]], int] = {}
        for _ in range(WORKERS):
            try:
                task = next(iterator)
            except StopIteration:
                break
            pending[pool.submit(run_one, root, task)] = int(task["job"])
        completed = 0
        while pending:
            done, _not_done = concurrent.futures.wait(
                pending, return_when=concurrent.futures.FIRST_COMPLETED)
            for future in done:
                job = pending.pop(future)
                try:
                    records.append(future.result())
                except BaseException:
                    for active in pending:
                        active.cancel()
                    raise
                completed += 1
                if completed % 100 == 0 or completed == len(tasks):
                    print("progress={}/{}".format(completed, len(tasks)), flush=True)
                try:
                    task = next(iterator)
                except StopIteration:
                    continue
                pending[pool.submit(run_one, root, task)] = int(task["job"])
    end_ns = time.monotonic_ns()
    ended_utc = utc_now()
    records.sort(key=lambda value: int(value["job"]))
    if [record["job"] for record in records] != list(range(len(tasks))):
        die("launch job receipt ledger is incomplete")
    filler_end = filler_snapshot()
    sampler_end = sampler_snapshot()
    validate_load_receipts(filler_end, sampler_end)
    if filler_end != filler_start:
        die("normal filler identity changed during the campaign")
    if sampler_end != sampler_start:
        die("external thermal sampler identity changed during the campaign")
    launch = sealed(SCHEMA + ".launch", {
        "started_utc": started_utc, "ended_utc": ended_utc,
        "start_monotonic_ns": start_ns, "end_monotonic_ns": end_ns,
        "duration_ns": end_ns - start_ns,
        "design_sha256": sha256_file(root / "design.json"),
        "prepare_receipt_sha256": sha256_file(root / "prepare_receipt.json"),
        "tasks_manifest_sha256": sha256_file(root / "tasks_manifest.jsonl"),
        "workers": WORKERS, "task_count": len(tasks),
        "job_receipts": records, "filler_start": filler_start,
        "filler_end": filler_end, "external_sampler_start": sampler_start,
        "external_sampler_end": sampler_end,
        "controller_started_sampler": False,
        "controller_signalled_sampler": False, "quiet_timing": False,
    })
    atomic_new(root / "launch_receipt.json", canonical_json(launch))
    print(json.dumps({
        "launched": True, "tasks": len(tasks),
        "duration_s": (end_ns - start_ns) / 1e9,
        "launch_receipt_sha256": sha256_file(root / "launch_receipt.json"),
        "filler_workers": filler_start["worker_count"],
        "external_sampler_processes": len(sampler_start),
    }, sort_keys=True))


def verify_runtime(
    root: Path, tasks: Sequence[Mapping[str, object]],
) -> Tuple[Dict[str, object], set[str]]:
    launch = load_canonical(root / "launch_receipt.json", "launch receipt")
    verify_sealed(launch, SCHEMA + ".launch", "launch receipt")
    if (set(launch) != LAUNCH_RECEIPT_FIELDS or
            launch.get("design_sha256") != sha256_file(root / "design.json") or
            launch.get("prepare_receipt_sha256") !=
                sha256_file(root / "prepare_receipt.json") or
            launch.get("tasks_manifest_sha256") !=
                sha256_file(root / "tasks_manifest.jsonl") or
            launch.get("workers") != WORKERS or
            launch.get("task_count") != len(tasks) or
            launch.get("controller_started_sampler") is not False or
            launch.get("controller_signalled_sampler") is not False or
            launch.get("quiet_timing") is not False or
            launch.get("filler_start") != launch.get("filler_end") or
            launch.get("external_sampler_start") !=
                launch.get("external_sampler_end")):
        die("launch receipt does not bind the recovery environment")
    start_ns = launch.get("start_monotonic_ns")
    end_ns = launch.get("end_monotonic_ns")
    if (not isinstance(start_ns, int) or isinstance(start_ns, bool) or
            not isinstance(end_ns, int) or isinstance(end_ns, bool) or
            end_ns <= start_ns or launch.get("duration_ns") != end_ns - start_ns):
        die("launch monotonic interval is malformed")
    filler = launch.get("filler_start")
    sampler = launch.get("external_sampler_start")
    validate_load_receipts(filler, sampler)
    ledger = launch.get("job_receipts")
    if not isinstance(ledger, list) or len(ledger) != len(tasks):
        die("launch job ledger is malformed")
    expected_files = {
        "design.json", "prepare_receipt.json", "tasks_manifest.jsonl",
        "launch_receipt.json",
    }
    design = load_canonical(root / "design.json", "campaign design")
    expected_files.update(str(value) for value in design["immutable_files"])
    expected_ledger = []
    for task in tasks:
        name = str(task["output_name"])
        raw_path = root / "raw" / name
        stderr_path = root / "stderr" / (name + ".stderr")
        exit_path = root / "exit" / (name + ".exit")
        receipt_path = root / "receipts" / (name + ".json")
        for path in (raw_path, stderr_path, exit_path, receipt_path):
            if path.is_symlink() or not path.is_file():
                die("runtime artifact is missing or indirect: {}".format(path))
        raw = raw_path.read_bytes()
        stderr = stderr_path.read_bytes()
        exit_raw = exit_path.read_bytes()
        if stderr or exit_raw != b"0\n":
            die("runtime stderr or exit evidence changed")
        rows = parse_output(raw, task)
        receipt = load_canonical(receipt_path, "job receipt")
        verify_sealed(receipt, SCHEMA + ".job", "job receipt")
        expected = {
            "job": task["job"], "task_sha256": sha256_bytes(canonical_json(task)),
            "argv": task["argv"], "output_name": name, "returncode": 0,
            "row_count": len(rows), "stdout_sha256": sha256_bytes(raw),
            "stdout_bytes": len(raw), "stderr_sha256": sha256_bytes(stderr),
            "stderr_bytes": len(stderr), "exit_sha256": sha256_bytes(exit_raw),
        }
        if (set(receipt) != JOB_RECEIPT_FIELDS or
                any(receipt.get(key) != value for key, value in expected.items())):
            die("job receipt does not replay: {}".format(task["job"]))
        job_start = receipt.get("start_monotonic_ns")
        job_end = receipt.get("end_monotonic_ns")
        if (not isinstance(job_start, int) or isinstance(job_start, bool) or
                not isinstance(job_end, int) or isinstance(job_end, bool) or
                not start_ns <= job_start < job_end <= end_ns or
                receipt.get("duration_ns") != job_end - job_start or
                not isinstance(receipt.get("started_utc"), str)):
            die("job timing receipt is malformed")
        ledger_record = {
            "job": task["job"], "name": name,
            "receipt_sha256": sha256_file(receipt_path),
        }
        expected_ledger.append(ledger_record)
        expected_files.update((
            "raw/" + name, "stderr/" + name + ".stderr",
            "exit/" + name + ".exit", "receipts/" + name + ".json"))
    if ledger != expected_ledger:
        die("launch job ledger does not replay")
    actual_files = {
        str(path.relative_to(root)) for path in root.rglob("*") if path.is_file()}
    optional = {"data_manifest.sha256", "validated_summary.json"}
    if not expected_files <= actual_files or actual_files - expected_files - optional:
        die("runtime campaign inventory contains missing or extra files")
    return launch, expected_files


def new_work() -> Dict[str, int]:
    return {"cells": 0, "xors": 0, "muladds": 0, "inactivated": 0}


def add_work(work: Dict[str, int], record: Mapping[str, object]) -> None:
    work["cells"] += 1
    for metric in ("xors", "muladds", "inactivated"):
        work[metric] += int(record[metric])


def finish_work(work: Mapping[str, int]) -> Dict[str, object]:
    cells = int(work["cells"])
    result: Dict[str, object] = {"cells": cells}
    for metric in ("xors", "muladds", "inactivated"):
        total = int(work[metric])
        result[metric + "_sum"] = total
        result[metric + "_mean"] = (
            str(Decimal(total) / Decimal(cells)) if cells else None)
    return result


def ratio(numerator: int, denominator: int) -> Optional[str]:
    return str(Decimal(numerator) / Decimal(denominator)) if denominator else None


def new_arm_accumulator() -> Dict[str, object]:
    return {
        "cells": 0, "outcomes": Counter(), "weak": Counter(),
        "by_seed": Counter(), "by_schedule": Counter(),
        "q_all": Counter(), "q_failure": Counter(), "attempts": Counter(),
        "failures": [], "work_all": new_work(), "work_success": new_work(),
        "work_failure": new_work(),
    }


def add_arm(
    accumulator: Dict[str, object], record: Mapping[str, object],
    seed_index: int, schedule: str,
) -> None:
    accumulator["cells"] = int(accumulator["cells"]) + 1
    outcome = str(record["outcome"])
    K = int(record["K"])
    q = int(record["q"])
    accumulator["outcomes"][outcome] += 1
    accumulator["q_all"][q] += 1
    accumulator["attempts"][int(record["seed_attempt"])] += 1
    add_work(accumulator["work_all"], record)
    if outcome == "success":
        add_work(accumulator["work_success"], record)
    else:
        add_work(accumulator["work_failure"], record)
        accumulator["weak"][K] += 1
        accumulator["by_seed"][seed_index] += 1
        accumulator["by_schedule"][schedule] += 1
        accumulator["q_failure"][q] += 1
        accumulator["failures"].append({
            "K": K, "seed_index": seed_index, "seed": SEEDS[seed_index],
            "schedule": schedule, "cause": outcome, "q": q,
            "seed_attempt": int(record["seed_attempt"]),
            "heavy_gain": int(record["heavy_gain"]),
        })


def counter_dict(counter: Counter) -> Dict[str, int]:
    return {str(key): int(counter[key]) for key in sorted(counter)}


def finish_arm(accumulator: Mapping[str, object]) -> Dict[str, object]:
    weak = accumulator["weak"]
    multiplicity_hist = Counter(weak.values())
    q_failure = accumulator["q_failure"]
    return {
        "cells": accumulator["cells"],
        "outcomes": counter_dict(accumulator["outcomes"]),
        "failures": len(accumulator["failures"]), "weak_K": len(weak),
        "maximum_weak_K_multiplicity": max(weak.values(), default=0),
        "weak_K_multiplicity_histogram": counter_dict(multiplicity_hist),
        "failures_by_seed_index": counter_dict(accumulator["by_seed"]),
        "failures_by_schedule": counter_dict(accumulator["by_schedule"]),
        "q_histogram_all": counter_dict(accumulator["q_all"]),
        "q_histogram_failures": counter_dict(q_failure),
        "q_gt_12_failures": sum(
            count for q, count in q_failure.items() if int(q) > 12),
        "seed_attempt_histogram": counter_dict(accumulator["attempts"]),
        "work_all": finish_work(accumulator["work_all"]),
        "work_success": finish_work(accumulator["work_success"]),
        "work_failure": finish_work(accumulator["work_failure"]),
        "failure_records": accumulator["failures"],
    }


def new_comparison() -> Dict[str, object]:
    return {
        "quadrants": Counter(), "repairs": [], "introductions": [],
        "both_fail": [], "q_transitions": Counter(),
        "common_base": new_work(), "common_candidate": new_work(),
        "all_base": new_work(), "all_candidate": new_work(),
        "seed_attempt_delta_cells": 0, "seed_attempt_delta_sum": 0,
        "q_delta_cells": 0, "q_delta_sum": 0,
    }


def comparison_key(
    K: int, seed_index: int, schedule: str,
    base: Mapping[str, object], candidate: Mapping[str, object],
) -> Dict[str, object]:
    return {
        "K": K, "seed_index": seed_index, "seed": SEEDS[seed_index],
        "schedule": schedule, "baseline_cause": base["outcome"],
        "candidate_cause": candidate["outcome"],
        "baseline_q": base["q"], "candidate_q": candidate["q"],
        "baseline_seed_attempt": base["seed_attempt"],
        "candidate_seed_attempt": candidate["seed_attempt"],
    }


def add_comparison(
    accumulator: Dict[str, object], base: Mapping[str, object],
    candidate: Mapping[str, object], seed_index: int, schedule: str,
) -> None:
    if base["K"] != candidate["K"]:
        die("paired architecture K coordinate changed")
    base_ok = base["outcome"] == "success"
    candidate_ok = candidate["outcome"] == "success"
    quadrant = (
        "both_success" if base_ok and candidate_ok else
        "repair" if not base_ok and candidate_ok else
        "introduction" if base_ok and not candidate_ok else "both_fail")
    accumulator["quadrants"][quadrant] += 1
    key = comparison_key(
        int(base["K"]), seed_index, schedule, base, candidate)
    if quadrant == "repair":
        accumulator["repairs"].append(key)
    elif quadrant == "introduction":
        accumulator["introductions"].append(key)
    elif quadrant == "both_fail":
        accumulator["both_fail"].append(key)
    accumulator["q_transitions"][(int(base["q"]), int(candidate["q"]))] += 1
    add_work(accumulator["all_base"], base)
    add_work(accumulator["all_candidate"], candidate)
    if quadrant == "both_success":
        add_work(accumulator["common_base"], base)
        add_work(accumulator["common_candidate"], candidate)
    attempt_delta = int(candidate["seed_attempt"]) - int(base["seed_attempt"])
    if attempt_delta:
        accumulator["seed_attempt_delta_cells"] += 1
        accumulator["seed_attempt_delta_sum"] += attempt_delta
    q_delta = int(candidate["q"]) - int(base["q"])
    if q_delta:
        accumulator["q_delta_cells"] += 1
        accumulator["q_delta_sum"] += q_delta


def work_pair(
    base: Mapping[str, int], candidate: Mapping[str, int],
) -> Dict[str, object]:
    base_finished = finish_work(base)
    candidate_finished = finish_work(candidate)
    return {
        "baseline": base_finished, "candidate": candidate_finished,
        "ratios": {
            metric: ratio(
                int(candidate_finished[metric + "_sum"]),
                int(base_finished[metric + "_sum"]))
            for metric in ("xors", "muladds", "inactivated")
        },
        "candidate_minus_baseline_sums": {
            metric: int(candidate_finished[metric + "_sum"]) -
                int(base_finished[metric + "_sum"])
            for metric in ("xors", "muladds", "inactivated")
        },
    }


def finish_comparison(accumulator: Mapping[str, object]) -> Dict[str, object]:
    quadrants = accumulator["quadrants"]
    repairs = accumulator["repairs"]
    introductions = accumulator["introductions"]
    return {
        "repairs": len(repairs), "introductions": len(introductions),
        "net_failure_change": len(introductions) - len(repairs),
        "both_success": quadrants["both_success"],
        "both_fail": quadrants["both_fail"],
        "repair_q_gt_12": sum(int(item["baseline_q"]) > 12 for item in repairs),
        "introduction_q_gt_12": sum(
            int(item["candidate_q"]) > 12 for item in introductions),
        "repair_records": repairs, "introduction_records": introductions,
        "both_fail_records": accumulator["both_fail"],
        "q_transition_histogram": {
            "{}->{}".format(before, after): count
            for (before, after), count in sorted(
                accumulator["q_transitions"].items())},
        "common_success_work": work_pair(
            accumulator["common_base"], accumulator["common_candidate"]),
        "all_cell_work": work_pair(
            accumulator["all_base"], accumulator["all_candidate"]),
        "seed_attempt_delta": {
            "changed_cells": accumulator["seed_attempt_delta_cells"],
            "candidate_minus_baseline_sum": accumulator["seed_attempt_delta_sum"],
        },
        "q_delta": {
            "changed_cells": accumulator["q_delta_cells"],
            "candidate_minus_baseline_sum": accumulator["q_delta_sum"],
        },
    }


def replay_summary(root: Path) -> Tuple[Dict[str, object], set[str], bytes]:
    design, tasks = load_prepared(root)
    launch, expected_files = verify_runtime(root, tasks)
    arms = {arm: new_arm_accumulator() for arm in ARMS}
    comparisons = {arm: new_comparison() for arm in ARMS[1:]}
    if len(tasks) % len(ARMS):
        die("parsed task ledger is not architecture-grouped")
    for offset in range(0, len(tasks), len(ARMS)):
        task_group = tasks[offset:offset + len(ARMS)]
        if tuple(str(task["arm"]) for task in task_group) != ARMS:
            die("parsed task architecture order changed")
        base_task = task_group[0]
        coordinates = (
            base_task["seed_index"], base_task["schedule"], base_task["Ks"])
        if any((task["seed_index"], task["schedule"], task["Ks"]) != coordinates
               for task in task_group[1:]):
            die("parsed task architecture coordinates differ")
        simplified = []
        for task in task_group:
            raw = (root / "raw" / str(task["output_name"])).read_bytes()
            rows = parse_output(raw, task)
            values = [validate_row(row) for row in rows]
            simplified.append(values)
            for record in values:
                add_arm(
                    arms[str(task["arm"])], record,
                    int(task["seed_index"]), str(task["schedule"]))
        for candidate_index, candidate_arm in enumerate(ARMS[1:], start=1):
            if len(simplified[0]) != len(simplified[candidate_index]):
                die("paired architecture row count changed")
            for base, candidate in zip(
                    simplified[0], simplified[candidate_index]):
                add_comparison(
                    comparisons[candidate_arm], base, candidate,
                    int(base_task["seed_index"]), str(base_task["schedule"]))
    finished_arms = {arm: finish_arm(arms[arm]) for arm in ARMS}
    finished_comparisons = {
        arm: finish_comparison(comparisons[arm]) for arm in ARMS[1:]}
    cells_per_arm = int(design["task_audit"]["cells_per_arm"])
    if any(summary["cells"] != cells_per_arm for summary in finished_arms.values()):
        die("reduced arm cell count differs from the frozen Cartesian cube")
    manifest = bytearray()
    for relative in sorted(expected_files):
        path = root / relative
        if path.is_symlink() or not path.is_file():
            die("data-manifest input is missing or indirect")
        manifest.extend("{}  {}\n".format(
            sha256_file(path), relative).encode("ascii"))
    manifest_bytes = bytes(manifest)
    payload = {
        "validation_issue_count": 0, "raw_architectures": True,
        "seed_fixes_applied": False,
        "prototype_commit": PROTOTYPE_COMMIT,
        "prototype_parent": PROTOTYPE_PARENT,
        "design_sha256": sha256_file(root / "design.json"),
        "prepare_receipt_sha256": sha256_file(root / "prepare_receipt.json"),
        "launch_receipt_sha256": sha256_file(root / "launch_receipt.json"),
        "tasks_manifest_sha256": sha256_file(root / "tasks_manifest.jsonl"),
        "data_manifest_sha256": sha256_bytes(manifest_bytes),
        "architecture": design["architecture"],
        "K_range": [K_LO, K_HI], "seeds": list(SEEDS),
        "schedules": list(SCHEDULES), "arms_order": list(ARMS),
        "cells_per_arm": cells_per_arm,
        "total_cells": int(design["task_audit"]["total_cells"]),
        "task_count": len(tasks),
        "execution_environment": design["execution_environment"],
        "load_receipt": {
            "filler_worker_count": launch["filler_start"]["worker_count"],
            "filler_cpus": launch["filler_start"]["cpus"],
            "external_sampler_processes": len(launch["external_sampler_start"]),
            "controller_started_sampler": False,
            "controller_signalled_sampler": False,
            "quiet_timing": False,
        },
        "arms": finished_arms, "comparisons_to_baseline": finished_comparisons,
        "selection_methodology": {
            "weak_K_and_multiplicity_are_descriptive_not_rejection_gates": True,
            "metrics": [
                "failures", "weak_K", "maximum_weak_K_multiplicity",
                "repairs", "introductions", "q>12 failure movement",
                "common-success XOR ratio", "common-success muladd ratio",
                "common-success inactivation ratio",
            ],
            "timing_claim": False,
            "seed_fix_evaluation_deferred_until_after_architecture_selection": True,
        },
    }
    return payload, expected_files, manifest_bytes


def reduce_campaign(args: argparse.Namespace) -> None:
    root = Path(args.root).resolve()
    for name in ("data_manifest.sha256", "validated_summary.json"):
        if (root / name).exists():
            die("reducer refuses to replace {}".format(name))
    payload, _expected_files, manifest = replay_summary(root)
    atomic_new(root / "data_manifest.sha256", manifest)
    summary = sealed(SCHEMA + ".validated", payload)
    encoded = canonical_json(summary)
    if len(encoded) > 16 * 1024 * 1024:
        die("validated summary exceeds its sealed size bound")
    atomic_new(root / "validated_summary.json", encoded)
    comparison = summary["comparisons_to_baseline"]
    print(json.dumps({
        "validated": True, "cells": summary["total_cells"],
        "failures": {
            arm: summary["arms"][arm]["failures"] for arm in ARMS},
        "repairs_introductions": {
            arm: [comparison[arm]["repairs"], comparison[arm]["introductions"]]
            for arm in ARMS[1:]},
        "validated_summary_sha256": sha256_file(
            root / "validated_summary.json"),
        "data_manifest_sha256": sha256_bytes(manifest),
    }, sort_keys=True))


def verify_campaign(args: argparse.Namespace) -> None:
    root = Path(args.root).resolve()
    payload, expected_files, manifest = replay_summary(root)
    if (root / "data_manifest.sha256").read_bytes() != manifest:
        die("data manifest does not replay")
    summary = load_canonical(root / "validated_summary.json", "validated summary")
    verify_sealed(summary, SCHEMA + ".validated", "validated summary")
    if summary != sealed(SCHEMA + ".validated", payload):
        die("validated summary does not replay")
    actual_files = {
        str(path.relative_to(root)) for path in root.rglob("*") if path.is_file()}
    allowed = set(expected_files) | {"data_manifest.sha256", "validated_summary.json"}
    if actual_files != allowed:
        die("reduced campaign inventory changed")
    print(json.dumps({
        "verified": True,
        "validated_summary_sha256": sha256_file(
            root / "validated_summary.json"),
        "data_manifest_sha256": sha256_bytes(manifest),
        "file_count": len(actual_files),
    }, sort_keys=True))


def selftest(_args: argparse.Namespace) -> None:
    expected_seeds = tuple(
        "0x" + hashlib.sha256(
            "wirehair.wh2.degree-balanced.allk.holdout.v1|seed|{}".format(index)
            .encode("ascii")).hexdigest()[:16]
        for index in range(3))
    if SEEDS != expected_seeds or len(set(SEEDS)) != len(SEEDS):
        die("sealed seed derivation selftest failed")
    tasks = build_tasks("/frozen/binary", k_lo=2, k_hi=19, chunk_max=5)
    audit = audit_tasks(tasks, k_lo=2, k_hi=19)
    if (len(tasks) != 144 or audit["task_groups"] != 36 or
            audit["cells_per_arm"] != 162 or audit["total_cells"] != 648):
        die("small task-cube cardinality selftest failed")
    for offset in range(0, len(tasks), len(ARMS)):
        group = tasks[offset:offset + len(ARMS)]
        if ("--binary-dense-two-anchor" not in group[0]["argv"] or
                any("--binary-dense-two-anchor" in task["argv"]
                    for task in group[1:]) or
                [task["argv"][-1] for task in group[1:]] != list(ARMS[1:])):
            die("raw architecture isolation selftest failed")
    if parse_preamble("# precodefail: a=1 b=two") != {"a": "1", "b": "two"}:
        die("preamble parser selftest failed")
    try:
        parse_preamble("# precodefail: a=1 a=2")
    except CampaignError:
        pass
    else:
        die("duplicate preamble selftest failed")
    base = {
        "K": 10, "outcome": "rank_fail", "inactivated": 20, "q": 13,
        "heavy_gain": 12, "seed_attempt": 0, "xors": 100, "muladds": 10,
    }
    candidate = dict(
        base, outcome="success", inactivated=19, q=12, xors=99, muladds=9)
    comparison = new_comparison()
    add_comparison(comparison, base, candidate, 0, "burst")
    finished = finish_comparison(comparison)
    if (finished["repairs"] != 1 or finished["repair_q_gt_12"] != 1 or
            finished["introductions"] != 0):
        die("paired reducer selftest failed")
    print(json.dumps({"selftest": "PASS", "small_task_audit": audit}, sort_keys=True))


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("--repo", required=True)
    prepare_parser.add_argument("--root", required=True)
    prepare_parser.add_argument("--workers", type=int, default=WORKERS)
    prepare_parser.set_defaults(function=prepare)
    verify_prepared_parser = subparsers.add_parser("verify-prepared")
    verify_prepared_parser.add_argument("--root", required=True)
    verify_prepared_parser.set_defaults(function=verify_prepared)
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--root", required=True)
    run_parser.add_argument("--workers", type=int, default=WORKERS)
    run_parser.set_defaults(function=run_campaign)
    reduce_parser = subparsers.add_parser("reduce")
    reduce_parser.add_argument("--root", required=True)
    reduce_parser.set_defaults(function=reduce_campaign)
    verify_parser = subparsers.add_parser("verify")
    verify_parser.add_argument("--root", required=True)
    verify_parser.set_defaults(function=verify_campaign)
    selftest_parser = subparsers.add_parser("selftest")
    selftest_parser.set_defaults(function=selftest)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = make_parser().parse_args(argv)
    try:
        args.function(args)
    except (CampaignError, OSError, ValueError) as exc:
        print("segmented-anchor campaign failed: {}".format(exc),
              file=sys.stderr, flush=True)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
