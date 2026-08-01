#!/usr/bin/env python3
"""Uniform-construction-root H12 versus H13 Stage-A campaign controller.

The controller deliberately has a small public surface:

    prepare  freeze this script, benchmark, roots, and every OH0 plan
    run      screen C0, then confirm C1+C2 and escalate paired failures to OH1024
    reduce   verify every seal again and emit the campaign summary

This is a deterministic census over named construction roots, loss roots, and
schedules.  McNemar values and K-cluster summaries are descriptive diagnostics,
not IID population claims.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from decimal import Decimal
import fcntl
from functools import wraps
import hashlib
import io
import json
import math
import os
from pathlib import Path
import queue
import re
import shutil
import signal
import stat
import subprocess
import sys
import tempfile
import threading
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from typing import Any, Iterable, Mapping, Sequence


SCHEMA = "wirehair.wh2.h13_stage_a.v2.uniform_construction_seed_v1"
K_MIN = 2
K_MAX = 64000
K_COUNT = K_MAX - K_MIN + 1
CHUNK_SIZE = 250
MAX_OVERHEAD = 1024
ARMS = ("h12", "h13")
ARM_GF256_ROWS = {"h12": 10, "h13": 11}
LOSS_ROOTS = (
    "0xd1b54a32d192ed03",
    "0x94d049bb133111eb",
    "0x8538ecb5bd456ea3",
)
CONSTRUCTION_SEEDS = (
    0x243F6A8885A308D3,
    0x13198A2E03707344,
    0xA4093822299F31D0,
)
CONSTRUCTION_SEED_POLICY = "matrix-c-peel-lo32-xor-hi32-v1"
SCHEDULES = ("burst", "adversarial", "repair-only")
PAIRED_CELLS_PER_ROOT = K_COUNT * len(LOSS_ROOTS) * len(SCHEDULES)
R0_PAIRED_CELLS = PAIRED_CELLS_PER_ROOT
R1_PAIRED_CELLS = PAIRED_CELLS_PER_ROOT * 2
FULL_PAIRED_CELLS = PAIRED_CELLS_PER_ROOT * len(CONSTRUCTION_SEEDS)
TOTAL_ARM_OUTCOMES_PER_ROOT = PAIRED_CELLS_PER_ROOT * len(ARMS)
FULL_ARM_OUTCOMES = FULL_PAIRED_CELLS * len(ARMS)
OH0_PAIR_JOBS_PER_ROOT = (
    math.ceil(K_COUNT / CHUNK_SIZE) * len(LOSS_ROOTS) * len(SCHEDULES)
)
OH0_JOBS_PER_ROOT = OH0_PAIR_JOBS_PER_ROOT * len(ARMS)
FULL_OH0_JOBS = OH0_JOBS_PER_ROOT * len(CONSTRUCTION_SEEDS)
# Compatibility aliases intentionally retain the per-root meaning used by the
# stage machinery.  Serialized contracts use the unambiguous names above.
SEEDS = LOSS_ROOTS
PAIRED_CELLS = PAIRED_CELLS_PER_ROOT
TOTAL_ARM_OUTCOMES = TOTAL_ARM_OUTCOMES_PER_ROOT
OH0_PAIR_JOBS = OH0_PAIR_JOBS_PER_ROOT
OH0_JOBS = OH0_JOBS_PER_ROOT
K_BANDS = (
    ("2-63", 2, 63),
    ("64-255", 64, 255),
    ("256-1023", 256, 1023),
    ("1024-4095", 1024, 4095),
    ("4096-16383", 4096, 16383),
    ("16384-64000", 16384, 64000),
)
MASK64 = (1 << 64) - 1
MAX_I63 = (1 << 63) - 1
OUTPUT_LIMIT = 8 << 20
EXECUTION_MODE = "authenticated-fd-taskset-v1"
LOSS_SEED_FORMULA = (
    "(external_seed ^ (K * 0x9e3779b97f4a7c15) ^ "
    "(64 * 0xbf58476d1ce4e5b9)) mod 2^64; paired overhead salt = 0"
)
THERMAL_FIELDS = (
    "utc", "monotonic_s", "cpu_busy_pct", "cpu_avg_mhz", "cpu_tctl_c",
    "dimm_i2c1_50_c", "dimm_i2c1_51_c", "dimm_i2c1_52_c",
    "dimm_i2c1_53_c", "dimm_i2c2_50_c", "dimm_i2c2_51_c",
    "dimm_i2c2_52_c", "dimm_i2c2_53_c", "dimm_read_errors",
    "load1", "load5", "load15", "edac_ce", "edac_ue",
)
THERMAL_DIMM_FIELDS = THERMAL_FIELDS[5:13]
TELEMETRY_POLICY = {
    "max_gap_seconds": 5.0, "min_cpu_busy_pct": 95.0,
    "max_temperature_c": 90.0, "require_zero_dimm_read_errors": True,
    "require_constant_edac_counters": True,
}
TELEMETRY_CONTINUITY = (
    "single append-only interval; sampler and independent >=95% filler remain "
    "continuously active from first run through terminal, including controller "
    "downtime and serial verification; any low row, >5s gap, rotation, or "
    "stopped filler irreversibly invalidates this result directory and requires "
    "a fresh campaign; shard restartability does not restart thermal provenance"
)

BENCH_HEADER = (
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
    "construction_attempt", "construction_seed_policy", "construction_seed",
    "production_seed_fixups_applied", "base_matrix_seed", "base_peel_seed",
    "matrix_seed", "peel_seed", "actual_staircase_rows",
    "actual_dense_rows", "actual_heavy_rows", "actual_source_hits",
    "actual_dense_two_anchor", "actual_dense_two_anchor_phase",
    "systematic_probe_result",
    "precode_count", "packet_trace_seed", "packet_trace_sha256",
)

THREE_DECIMAL_FIELDS = (
    "inact_mu", "binary_def_mu", "heavy_gain_mu", "solve_ms_mu",
    "build_ms_mu", "peel_ms_mu", "project_ms_mu", "residual_ms_mu",
    "backsub_ms_mu", "block_xors_mu", "block_muladds_mu",
    "mixed_joint_source_xors_mu", "mixed_joint_marginal_xors_mu",
    "mixed_joint_marginal_copies_mu", "mixed_joint_active_deltas_mu",
    "mixed_joint_scratch_bytes_mu", "mixed_dual_source_columns_mu",
)
WORK_FIELDS = (
    "block_xors_mu", "block_muladds_mu", "mixed_joint_source_xors_mu",
    "mixed_joint_marginal_xors_mu", "mixed_joint_marginal_copies_mu",
    "mixed_joint_active_deltas_mu", "mixed_joint_scratch_bytes_mu",
    "mixed_dual_source_columns_mu",
)


class CampaignError(RuntimeError):
    """A fail-closed campaign contract violation."""


class CampaignInterrupted(CampaignError):
    """A controlled SIGINT/SIGTERM shutdown."""

    def __init__(self, signum: int) -> None:
        self.signum = signum
        super().__init__(f"campaign interrupted by signal {signum}")


class ActiveChildren:
    """Own all benchmark process groups launched by one run invocation."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._stopped = threading.Event()
        self._processes: dict[int, subprocess.Popen[bytes]] = {}
        self._signal: int | None = None

    @staticmethod
    def _kill(process: subprocess.Popen[bytes]) -> None:
        if process.poll() is not None:
            return
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass

    def register(self, process: subprocess.Popen[bytes]) -> None:
        with self._lock:
            self._processes[process.pid] = process
            stopped = self._stopped.is_set()
        # Close the launch-versus-signal race: a process registered after the
        # stop snapshot must kill itself before it can produce a seal.
        if stopped:
            self._kill(process)

    def unregister(self, process: subprocess.Popen[bytes]) -> None:
        with self._lock:
            self._processes.pop(process.pid, None)

    def stop(self, signum: int | None = None) -> None:
        with self._lock:
            if signum is not None and self._signal is None:
                self._signal = signum
            self._stopped.set()
            processes = list(self._processes.values())
        for process in processes:
            self._kill(process)

    def check(self) -> None:
        if self._stopped.is_set():
            if self._signal is not None:
                raise CampaignInterrupted(self._signal)
            die("campaign dispatch cancelled after a job failure")

    @property
    def interrupted_signal(self) -> int | None:
        with self._lock:
            return self._signal

    @property
    def active_count(self) -> int:
        with self._lock:
            return len(self._processes)


def die(message: str) -> None:
    raise CampaignError(message)


def pdeath_exec(expected_parent: int, command: Sequence[str]) -> int:
    """Exec a command that the Linux kernel kills if this controller dies."""
    if (sys.platform != "linux" or type(expected_parent) is not int or
            expected_parent <= 1 or not command):
        die("parent-death launcher requires Linux and a nonempty command")
    # This runs in a fresh, single-threaded Python launcher rather than in a
    # preexec_fn after fork from the multithreaded controller.
    import ctypes

    parent = os.getppid()
    if parent != expected_parent:
        os.kill(os.getpid(), signal.SIGKILL)
    libc = ctypes.CDLL(None, use_errno=True)
    if libc.prctl(1, signal.SIGKILL, 0, 0, 0) != 0:  # PR_SET_PDEATHSIG
        error = ctypes.get_errno()
        die(f"cannot arm benchmark parent-death signal: errno {error}")
    if os.getppid() != expected_parent:
        os.kill(os.getpid(), signal.SIGKILL)
    try:
        os.execv(command[0], list(command))
    except OSError as exc:
        die(f"cannot exec benchmark command: {exc}")
    return 125


def canonical_json(value: Any) -> str:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False,
    )


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def stable_source_bytes(path: Path, *, require_unique: bool) -> bytes:
    flags = (os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
             getattr(os, "O_NOFOLLOW", 0))
    try:
        descriptor = os.open(path, flags)
    except OSError as exc:
        die(f"cannot open sealed regular file {path}: {exc}")
    try:
        before = os.fstat(descriptor)
        if (not stat.S_ISREG(before.st_mode) or
                (require_unique and before.st_nlink != 1)):
            die(f"sealed input is not a unique regular file: {path}")
        chunks: list[bytes] = []
        while True:
            chunk = os.read(descriptor, 1 << 20)
            if not chunk:
                break
            chunks.append(chunk)
        midpoint = os.fstat(descriptor)
        os.lseek(descriptor, 0, os.SEEK_SET)
        confirmation_chunks: list[bytes] = []
        while True:
            chunk = os.read(descriptor, 1 << 20)
            if not chunk:
                break
            confirmation_chunks.append(chunk)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    identity = lambda value: (
        value.st_dev, value.st_ino, value.st_mode, value.st_nlink,
        value.st_size, value.st_mtime_ns, value.st_ctime_ns,
    )
    try:
        named = os.stat(path, follow_symlinks=False)
    except OSError as exc:
        die(f"sealed input pathname disappeared: {path}: {exc}")
    if (identity(before) != identity(midpoint) or
            identity(midpoint) != identity(after) or
            identity(after) != identity(named) or
            not stat.S_ISREG(named.st_mode) or
            (require_unique and named.st_nlink != 1)):
        die(f"sealed input changed while reading: {path}")
    data = b"".join(chunks)
    confirmation = b"".join(confirmation_chunks)
    if len(data) != before.st_size or data != confirmation:
        die(f"short read from sealed input: {path}")
    return data


def stable_bytes(path: Path) -> bytes:
    return stable_source_bytes(path, require_unique=True)


def sha256_file(path: Path) -> str:
    return sha256_bytes(stable_bytes(path))


def external_file_mark(path: Path) -> dict[str, Any]:
    try:
        metadata = os.stat(path, follow_symlinks=False)
    except OSError as exc:
        die(f"cannot stat external telemetry log {path}: {exc}")
    if not stat.S_ISREG(metadata.st_mode) or metadata.st_nlink != 1:
        die("external telemetry log is not a unique regular file")
    return {
        "path": str(path), "dev": metadata.st_dev, "ino": metadata.st_ino,
        "offset": metadata.st_size, "mtime_ns": metadata.st_mtime_ns,
    }


def validate_external_mark(
    mark: Any, expected_path: str, context: str,
) -> None:
    if (not isinstance(mark, dict) or set(mark) != {
            "path", "dev", "ino", "offset", "mtime_ns"} or
            type(mark["path"]) is not str or mark["path"] != expected_path or
            any(type(mark[key]) is not int for key in (
                "dev", "ino", "offset", "mtime_ns")) or
            mark["dev"] < 0 or mark["ino"] <= 0 or mark["offset"] < 0):
        die(f"{context}: external telemetry mark is malformed")


def external_range_bytes(
    path: Path, start: int, end: int,
    expected_identity: tuple[int, int] | None = None,
) -> bytes:
    if start < 0 or end < start:
        die("external telemetry interval is malformed")
    if end - start > 256 << 20:
        die("external telemetry interval exceeds 256 MiB")
    flags = (os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
             getattr(os, "O_NOFOLLOW", 0))
    try:
        descriptor = os.open(path, flags)
    except OSError as exc:
        die(f"cannot open external telemetry interval: {exc}")
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode) or before.st_size < end:
            die("external telemetry interval is truncated")
        chunks: list[bytes] = []
        offset = start
        while offset < end:
            block = os.pread(descriptor, min(1 << 20, end - offset), offset)
            if not block:
                die("short read from external telemetry interval")
            chunks.append(block)
            offset += len(block)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    try:
        named = os.stat(path, follow_symlinks=False)
    except OSError as exc:
        die(f"external telemetry pathname disappeared: {exc}")
    identity = lambda value: (
        value.st_dev, value.st_ino, value.st_mode, value.st_nlink,
    )
    if (identity(before) != identity(after) or
            identity(after) != identity(named) or
            (expected_identity is not None and
             (before.st_dev, before.st_ino) != expected_identity) or
            not stat.S_ISREG(named.st_mode) or named.st_nlink != 1 or
            after.st_size < end or named.st_size < end):
        die("external telemetry changed identity while hashing")
    return b"".join(chunks)


def _telemetry_scalar(text: str, context: str) -> float:
    if not re.fullmatch(r"[0-9]+(?:\.[0-9]+)?", text):
        die(f"{context}: noncanonical telemetry scalar")
    value = float(text)
    if not math.isfinite(value):
        die(f"{context}: nonfinite telemetry scalar")
    return value


def _telemetry_counter(text: str, context: str) -> int:
    if not re.fullmatch(r"0|[1-9][0-9]*", text):
        die(f"{context}: noncanonical telemetry counter")
    return int(text)


def system_uptime_s() -> float:
    try:
        fields = Path("/proc/uptime").read_text(encoding="ascii").split()
    except OSError as exc:
        die(f"cannot read canonical system uptime: {exc}")
    if len(fields) != 2:
        die("canonical system uptime has the wrong schema")
    return _telemetry_scalar(fields[0], "system uptime")


def parse_telemetry_row(line: bytes, context: str) -> dict[str, Any]:
    if len(line) > 16384 or not line.endswith(b"\n"):
        die(f"{context}: incomplete or oversized telemetry row")
    try:
        text = line.decode("utf-8")[:-1]
    except UnicodeDecodeError as exc:
        die(f"{context}: telemetry row is not UTF-8: {exc}")
    if text.endswith("\r"):
        text = text[:-1]
    values = text.split(",")
    if len(values) != len(THERMAL_FIELDS):
        die(f"{context}: telemetry field count mismatch")
    row = dict(zip(THERMAL_FIELDS, values))
    if not row["utc"]:
        die(f"{context}: empty UTC field")
    monotonic = _telemetry_scalar(row["monotonic_s"], f"{context}:monotonic")
    cpu_busy = _telemetry_scalar(row["cpu_busy_pct"], f"{context}:busy")
    cpu_mhz = _telemetry_scalar(row["cpu_avg_mhz"], f"{context}:mhz")
    cpu_temp = _telemetry_scalar(row["cpu_tctl_c"], f"{context}:cpu_temp")
    if any(row[field] == "" for field in THERMAL_DIMM_FIELDS):
        die(f"{context}: missing DIMM temperature")
    dimm_temperatures = [
        _telemetry_scalar(row[field], f"{context}:{field}")
        for field in THERMAL_DIMM_FIELDS
    ]
    dimm_errors = _telemetry_counter(
        row["dimm_read_errors"], f"{context}:dimm errors",
    )
    loads = [
        _telemetry_scalar(row[field], f"{context}:{field}")
        for field in ("load1", "load5", "load15")
    ]
    edac_ce = _telemetry_counter(row["edac_ce"], f"{context}:edac_ce")
    edac_ue = _telemetry_counter(row["edac_ue"], f"{context}:edac_ue")
    if (not 0 <= cpu_busy <= 100 or cpu_mhz <= 0 or
            not 0 <= cpu_temp <= 120 or
            any(not 0 <= value <= 100 for value in dimm_temperatures)):
        die(f"{context}: telemetry value outside physical bounds")
    return {
        "utc": row["utc"], "monotonic_s": monotonic,
        "cpu_busy_pct": cpu_busy, "cpu_avg_mhz": cpu_mhz,
        "cpu_tctl_c": cpu_temp, "dimm_temperatures": dimm_temperatures,
        "dimm_read_errors": dimm_errors, "loads": loads,
        "edac_ce": edac_ce, "edac_ue": edac_ue,
    }


def telemetry_baseline(
    path: Path, end: int, expected_identity: tuple[int, int],
) -> dict[str, Any]:
    expected_header = (",".join(THERMAL_FIELDS) + "\n").encode("ascii")
    if external_range_bytes(
            path, 0, len(expected_header), expected_identity) != expected_header:
        die("external telemetry header mismatch")
    start = max(len(expected_header), end - 16384)
    suffix = external_range_bytes(path, start, end, expected_identity)
    if not suffix.endswith(b"\n"):
        die("external telemetry start is not at a complete row boundary")
    previous_newline = suffix[:-1].rfind(b"\n")
    line = suffix[previous_newline + 1:]
    if not line or line == expected_header:
        die("external telemetry has no baseline sample")
    return parse_telemetry_row(line, "telemetry baseline")


def audit_telemetry_interval(
    path: Path, start: Mapping[str, Any], end: Mapping[str, Any],
    audit_uptime_s: float | None = None,
) -> dict[str, Any]:
    expected_identity = (start["dev"], start["ino"])
    if (end["dev"], end["ino"]) != expected_identity:
        die("external telemetry interval endpoints differ in identity")
    suffix = external_range_bytes(
        path, start["offset"], end["offset"], expected_identity,
    )
    if not suffix or not suffix.endswith(b"\n"):
        die("external telemetry interval is empty or ends in a partial row")
    rows = suffix.splitlines(keepends=True)
    samples = [start["baseline"]] + [
        parse_telemetry_row(line, f"telemetry interval row {index}")
        for index, line in enumerate(rows, 1)
    ]
    gaps = [
        current["monotonic_s"] - previous["monotonic_s"]
        for previous, current in zip(samples, samples[1:])
    ]
    if any(delta <= 0 or delta > TELEMETRY_POLICY["max_gap_seconds"] for delta in gaps):
        die("external telemetry interval is nonmonotonic or has a sample gap")
    interval = samples[1:]
    if audit_uptime_s is None:
        audit_uptime_s = system_uptime_s()
    if type(audit_uptime_s) is not float or not math.isfinite(audit_uptime_s):
        die("external telemetry audit clock is malformed")
    tail_age_seconds = audit_uptime_s - interval[-1]["monotonic_s"]
    if (tail_age_seconds < 0 or
            tail_age_seconds > TELEMETRY_POLICY["max_gap_seconds"]):
        die("external telemetry sampler was not live at interval sealing")
    if min(sample["cpu_busy_pct"] for sample in interval) < \
            TELEMETRY_POLICY["min_cpu_busy_pct"]:
        die("external telemetry CPU busy floor was not maintained")
    if max(
            max(sample["cpu_tctl_c"], *sample["dimm_temperatures"])
            for sample in samples) >= TELEMETRY_POLICY["max_temperature_c"]:
        die("external telemetry temperature ceiling was reached")
    if any(sample["dimm_read_errors"] != 0 for sample in samples):
        die("external telemetry reports a DIMM read error")
    if (len({sample["edac_ce"] for sample in samples}) != 1 or
            len({sample["edac_ue"] for sample in samples}) != 1):
        die("external telemetry EDAC counters changed or rolled back")
    return {
        "samples": len(interval), "max_gap_seconds": max(gaps),
        "cpu_busy_min_pct": min(sample["cpu_busy_pct"] for sample in interval),
        "cpu_tctl_max_c": max(sample["cpu_tctl_c"] for sample in samples),
        "dimm_max_c": max(
            value for sample in samples for value in sample["dimm_temperatures"]
        ),
        "dimm_read_errors_max": 0,
        "edac_ce": samples[0]["edac_ce"],
        "edac_ue": samples[0]["edac_ue"],
        "edac_ce_delta": 0, "edac_ue_delta": 0,
        "audit_uptime_s": audit_uptime_s,
        "tail_age_seconds": tail_age_seconds,
        "interval_sha256": sha256_bytes(suffix),
    }


def sealed_record(schema: str, payload: Any) -> dict[str, Any]:
    body = {"schema": schema, "payload": payload}
    return {**body, "seal_sha256": sha256_bytes(canonical_json(body).encode())}


def verify_sealed(record: Any, schema: str) -> Any:
    if not isinstance(record, dict) or set(record) != {
            "schema", "payload", "seal_sha256"}:
        die(f"{schema}: malformed sealed record")
    if record["schema"] != schema:
        die(f"sealed schema {record['schema']!r}, want {schema!r}")
    seal = record["seal_sha256"]
    if not isinstance(seal, str) or not re.fullmatch(r"[0-9a-f]{64}", seal):
        die(f"{schema}: malformed seal")
    body = {"schema": schema, "payload": record["payload"]}
    if sha256_bytes(canonical_json(body).encode()) != seal:
        die(f"{schema}: seal mismatch")
    return record["payload"]


def load_sealed(path: Path, schema: str) -> Any:
    try:
        record = json.loads(stable_bytes(path).decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        die(f"{path}: invalid sealed JSON: {exc}")
    if stable_bytes(path) != (canonical_json(record) + "\n").encode():
        die(f"{path}: JSON is not canonical")
    return verify_sealed(record, schema)


def atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", dir=path.parent,
    )
    temporary_path = Path(temporary)
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(data)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary_path, path)
    finally:
        try:
            temporary_path.unlink()
        except FileNotFoundError:
            pass


def write_sealed_once(path: Path, schema: str, payload: Any) -> None:
    data = (canonical_json(sealed_record(schema, payload)) + "\n").encode()
    if path.exists():
        if stable_bytes(path) != data:
            die(f"refusing to replace a different sealed record: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.create-", dir=path.parent,
    )
    temporary_path = Path(temporary)
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(data)
            output.flush()
            os.fsync(output.fileno())
        try:
            # A same-directory hard link is atomic and, unlike os.replace(),
            # cannot overwrite a different seal created after the exists()
            # check above.  The temporary name is removed immediately below.
            os.link(temporary_path, path, follow_symlinks=False)
        except FileExistsError:
            if stable_bytes(path) != data:
                die(f"refusing to replace a different sealed record: {path}")
    finally:
        try:
            temporary_path.unlink()
        except FileNotFoundError:
            pass


def checked_product(*values: int) -> int:
    result = 1
    for value in values:
        if (type(value) is not int or value < 0 or
                result > MAX_I63 // max(value, 1)):
            die("campaign cardinality exceeds signed 63-bit accounting")
        result *= value
    return result


def strict_uint(text: str, context: str, maximum: int | None = None) -> int:
    if (type(text) is not str or len(text) > 20 or
            not re.fullmatch(r"0|[1-9][0-9]*", text)):
        die(f"{context}: noncanonical unsigned integer {text!r}")
    value = int(text)
    if maximum is not None and value > maximum:
        die(f"{context}: integer exceeds {maximum}")
    return value


def strict_hex(text: str, context: str, bits: int) -> int:
    if (type(text) is not str or type(bits) is not int or bits < 1 or
            len(text) > 2 + (bits + 3) // 4 or
            not re.fullmatch(r"0x(?:0|[1-9a-f][0-9a-f]*)", text)):
        die(f"{context}: noncanonical hexadecimal integer {text!r}")
    value = int(text, 16)
    if value >= 1 << bits:
        die(f"{context}: integer does not fit u{bits}")
    return value


def strict_fixed(text: str, context: str, places: int) -> Decimal:
    if (type(text) is not str or type(places) is not int or
            not 0 <= places <= 9 or len(text) > 32 or
            not re.fullmatch(rf"(?:0|[1-9][0-9]*)\.[0-9]{{{places}}}", text)):
        die(f"{context}: noncanonical fixed-point value {text!r}")
    return Decimal(text)


def validate_histogram(text: str, context: str) -> int:
    if not text:
        die(f"{context}: empty histogram")
    previous = -1
    total = 0
    for token in text.split(","):
        value_text, separator, count_text = token.partition(":")
        if not separator:
            die(f"{context}: malformed histogram token")
        value = strict_uint(value_text, context)
        count = strict_uint(count_text, context)
        if count == 0 or value <= previous:
            die(f"{context}: histogram is empty, duplicate, or unsorted")
        previous = value
        total += count
    if total != 1:
        die(f"{context}: one-trial histogram count is not one")
    return previous


def parse_preamble(line: str) -> dict[str, str]:
    prefix = "# precodefail: "
    if not line.startswith(prefix):
        die("missing raw precodefail preamble")
    result: dict[str, str] = {}
    for token in line[len(prefix):].split(" "):
        key, separator, value = token.partition("=")
        if not separator or not key or not value or key in result:
            die(f"malformed or duplicate preamble token {token!r}")
        result[key] = value
    return result


def expected_preamble(task: "Task") -> dict[str, str]:
    return {
        "trials": "1", "threads": "1", "loss": "0.5",
        "seed": task.seed, "seed_role": "loss-trace-root",
        "construction_seed_policy": CONSTRUCTION_SEED_POLICY,
        "construction_seed": str(task.construction_seed),
        "production_seed_fixups_applied": "0",
        "raw_attempt0": "1",
        "completion": "mixed", "mixed_period": "48",
        "mixed_gf256_rows": str(ARM_GF256_ROWS[task.arm]),
        "mixed_gf16_rows": "2", "mixed_geometry": "shared-x",
        "mixed_residue_skew": "0", "mixed_residue_schedule": "constant",
        "mixed_residue_hash_seed": "0x0", "mixed_residue_hash_keyed": "0",
        "mixed_independent_extension_residues": "0",
        "mixed_grouped_gf256_rows": "0",
        "mixed_grouped_gf256_row_mask": "0x0",
        "mixed_grouped_gf256_hash_seed": "0x0",
        "mixed_grouped_final_h_a_columns": "0",
        "mixed_residue_buckets_requested": "auto",
        "mixed_extension_residue_seed_xor": "0x4e",
        "source_hits_override": "0", "packet_peel_seed_xor": "0x0",
        "packet_peel_seed_table": "none", "binary_dense_rows_override": "0",
        "binary_dense_two_anchor": "1", "binary_dense_two_anchor_phase": "0",
        "gf256_heavy_rows_override": "0", "odd_packet_peel_seed_xor": "0x0",
        "packet_row_seed_multiplier": "0x1",
        "packet_row_seed_avalanche": "0",
        "seed_block_bytes_override": "0", "overhead_stream": "paired",
        "full_payload_solve": "0", "schedule": task.schedule,
    }


def expected_preamble_line(task: "Task") -> str:
    return "# precodefail: " + " ".join(
        f"{key}={value}" for key, value in expected_preamble(task).items()
    )


def loss_seed(seed: str, K: int) -> int:
    return (
        int(seed, 16) ^ (K * 0x9E3779B97F4A7C15) ^
        (64 * 0xBF58476D1CE4E5B9)
    ) & MASK64


def construction_peel_seed(construction_seed: int) -> int:
    if (type(construction_seed) is not int or construction_seed <= 0 or
            construction_seed > MASK64):
        die("construction seed is not a nonzero u64")
    return ((construction_seed & 0xffffffff) ^
            (construction_seed >> 32)) & 0xffffffff


@dataclass(frozen=True)
class Task:
    overhead: int
    job: int
    pair: int
    arm: str
    construction_index: int
    construction_seed: int
    seed_index: int
    seed: str
    schedule: str
    ks: tuple[int, ...]
    task_id: str

    @property
    def stem(self) -> str:
        return f"job{self.job:07d}.{self.arm}.{self.task_id}"

    def core(self) -> dict[str, Any]:
        return {
            "overhead": self.overhead, "job": self.job, "pair": self.pair,
            "arm": self.arm,
            "construction_index": self.construction_index,
            "construction_seed": self.construction_seed,
            "seed_index": self.seed_index,
            "seed": self.seed, "schedule": self.schedule,
            "ks": list(self.ks),
        }

    def payload(self) -> dict[str, Any]:
        return {**self.core(), "task_id": self.task_id}


def task_from_payload(payload: Any) -> Task:
    required = {
        "overhead", "job", "pair", "arm", "construction_index",
        "construction_seed", "seed_index", "seed", "schedule", "ks",
        "task_id",
    }
    if not isinstance(payload, dict) or set(payload) != required:
        die("task has an unexpected schema")
    if (any(type(payload[key]) is not int for key in (
            "overhead", "job", "pair", "construction_index",
            "construction_seed", "seed_index")) or
            any(type(payload[key]) is not str for key in (
                "arm", "seed", "schedule", "task_id")) or
            type(payload["ks"]) is not list or
            any(type(value) is not int for value in payload["ks"])):
        die("task scalar types are not canonical")
    task = Task(
        overhead=payload["overhead"], job=payload["job"],
        pair=payload["pair"], arm=payload["arm"],
        construction_index=payload["construction_index"],
        construction_seed=payload["construction_seed"],
        seed_index=payload["seed_index"], seed=payload["seed"],
        schedule=payload["schedule"], ks=tuple(payload["ks"]),
        task_id=payload["task_id"],
    )
    validate_task(task)
    return task


def task_digest(core: Mapping[str, Any]) -> str:
    return sha256_bytes(canonical_json(core).encode())[:24]


def validate_task(task: Task) -> None:
    if (any(type(value) is not int for value in (
            task.overhead, task.job, task.pair, task.construction_index,
            task.construction_seed, task.seed_index, *task.ks)) or
            any(type(value) is not str for value in (
                task.arm, task.seed, task.schedule, task.task_id)) or
            task.overhead < 0 or task.overhead > MAX_OVERHEAD or task.job < 0 or
            task.pair < 0 or task.arm not in ARMS or
            task.construction_index not in range(len(CONSTRUCTION_SEEDS)) or
            task.construction_seed != CONSTRUCTION_SEEDS[task.construction_index] or
            task.seed_index not in range(len(SEEDS)) or
            task.seed != SEEDS[task.seed_index] or
            task.schedule not in SCHEDULES or not task.ks or
            len(task.ks) > CHUNK_SIZE or tuple(sorted(task.ks)) != task.ks or
            len(set(task.ks)) != len(task.ks) or
            any(K < K_MIN or K > K_MAX for K in task.ks)):
        die("task violates the Stage-A geometry")
    if task.task_id != task_digest(task.core()):
        die("task ID does not bind its exact payload")


def _new_task(
    overhead: int, job: int, pair: int, arm: str, construction_index: int,
    seed_index: int, schedule: str, ks: Sequence[int],
) -> Task:
    provisional = Task(
        overhead, job, pair, arm, construction_index,
        CONSTRUCTION_SEEDS[construction_index], seed_index, SEEDS[seed_index],
        schedule, tuple(ks), "",
    )
    return Task(**{**provisional.__dict__, "task_id": task_digest(provisional.core())})


def build_tasks(
    overhead: int, construction_index: int = 0,
    cohort: Iterable[tuple[int, str, int]] | None = None,
) -> list[Task]:
    if type(overhead) is not int or overhead < 0 or overhead > MAX_OVERHEAD:
        die("overhead outside 0..1024")
    if (type(construction_index) is not int or
            construction_index not in range(len(CONSTRUCTION_SEEDS))):
        die("construction root index outside the sealed domain")
    grouped: dict[tuple[int, str], list[int]] = {}
    if cohort is None:
        if overhead != 0:
            die("only OH0 may use the full implicit cohort")
        for seed_index in range(len(SEEDS)):
            for schedule in SCHEDULES:
                grouped[(seed_index, schedule)] = list(range(K_MIN, K_MAX + 1))
    else:
        seen: set[tuple[int, str, int]] = set()
        for seed_index, schedule, K in cohort:
            key = (seed_index, schedule, K)
            if key in seen or seed_index not in range(len(SEEDS)) or \
                    schedule not in SCHEDULES or not K_MIN <= K <= K_MAX:
                die("follow-up cohort is duplicated or malformed")
            seen.add(key)
            grouped.setdefault((seed_index, schedule), []).append(K)
        for values in grouped.values():
            values.sort()

    tasks: list[Task] = []
    pair = 0
    for seed_index in range(len(SEEDS)):
        for schedule in SCHEDULES:
            values = grouped.get((seed_index, schedule), [])
            for offset in range(0, len(values), CHUNK_SIZE):
                ks = values[offset:offset + CHUNK_SIZE]
                for arm in ARMS:
                    tasks.append(_new_task(
                        overhead, len(tasks), pair, arm, construction_index,
                        seed_index,
                        schedule, ks,
                    ))
                pair += 1
    def paired_core(value: Task) -> dict[str, Any]:
        core = value.core()
        core["arm"] = "paired"
        core["job"] = -1
        return core

    if len(tasks) % 2 or any(
            paired_core(tasks[index]) != paired_core(tasks[index + 1])
            for index in range(0, len(tasks), 2)):
        die("task planner did not emit adjacent exact arm pairs")
    if cohort is None:
        if (len(tasks) != OH0_JOBS or pair != OH0_PAIR_JOBS or
                sum(len(task.ks) for task in tasks) != TOTAL_ARM_OUTCOMES):
            die("OH0 task cardinality is not 4,608 jobs / 575,991 cells per arm")
    return tasks


def make_benchmark_argv(binary: Path, task: Task) -> list[str]:
    validate_task(task)
    return [
        str(binary), "precodefail", "--N", ",".join(map(str, task.ks)),
        "--bb-list", "64", "--overhead", str(task.overhead),
        "--trials", "1", "--threads", "1", "--loss", "0.5",
        "--completion", "mixed", "--heavy-family", "periodic",
        "--mix-count", "2", "--mixed-period", "48",
        "--mixed-gf256-rows", str(ARM_GF256_ROWS[task.arm]),
        "--mixed-gf16-rows", "2", "--mixed-geometry", "shared-x",
        "--mixed-residue-skew", "0", "--mixed-residue-schedule", "constant",
        "--binary-dense-two-anchor", "--binary-dense-two-anchor-phase", "0",
        "--packet-peel-seed-xor", "0",
        "--paired-overhead-stream", "--seed", task.seed,
        "--schedule", task.schedule, "--raw-attempt0",
        "--construction-seed", str(task.construction_seed),
    ]


def parse_output(output: str, task: Task) -> list[dict[str, str]]:
    if "\r" in output or not output.endswith("\n"):
        die("benchmark output is not canonical LF-terminated text")
    lines = output.splitlines()
    if len(lines) != len(task.ks) + 2:
        die("benchmark output line count does not match the sealed K list")
    if (parse_preamble(lines[0]) != expected_preamble(task) or
            lines[0] != expected_preamble_line(task)):
        die("benchmark preamble differs from the sealed Stage-A task")
    try:
        header = tuple(next(csv.reader([lines[1]], strict=True)))
    except csv.Error as exc:
        die(f"malformed benchmark header: {exc}")
    if (lines[1] != ",".join(BENCH_HEADER) or header != BENCH_HEADER):
        die("benchmark CSV header differs from the raw receipt schema")
    rows: list[dict[str, str]] = []
    for index, line in enumerate(lines[2:]):
        try:
            values = next(csv.reader([line], strict=True))
        except csv.Error as exc:
            die(f"malformed benchmark row: {exc}")
        if (line != ",".join(values) or len(values) != len(BENCH_HEADER)):
            die("benchmark row has the wrong field count")
        row = dict(zip(BENCH_HEADER, values))
        K = strict_uint(row["N"], "N", K_MAX)
        if K != task.ks[index]:
            die("benchmark K order differs from the sealed task")
        fixed = {
            "bb": "64", "heavy_family": "periodic", "mix_count": "2",
            "overhead": str(task.overhead), "trials": "1",
            "seed_attempt": "0", "active_packet_peel_seed_xor": "0x0",
            "construction_attempt": "0",
            "construction_seed_policy": CONSTRUCTION_SEED_POLICY,
            "construction_seed": str(task.construction_seed),
            "production_seed_fixups_applied": "0",
        }
        if any(row[key] != value for key, value in fixed.items()):
            die(f"benchmark row geometry mismatch at K={K}")
        outcomes = [strict_uint(row[key], f"K={K}:{key}", 1)
                    for key in ("success", "rank_fail", "error")]
        if sum(outcomes) != 1:
            die(f"benchmark outcome is not exactly one trial at K={K}")
        if outcomes[2]:
            die(f"benchmark reported a codec/internal error at K={K}")
        fail_rate = strict_fixed(row["fail_rate"], f"K={K}:fail_rate", 8)
        if fail_rate != Decimal(outcomes[1] + outcomes[2]):
            die(f"benchmark failure rate disagrees with outcome at K={K}")
        fixed_values = {
            key: strict_fixed(row[key], f"K={K}:{key}", 3)
            for key in THREE_DECIMAL_FIELDS
        }
        for key in (
            "inact_max", "binary_def_max", "heavy_gain_min",
            "heavy_shortfall", "actual_staircase_rows", "actual_dense_rows",
            "actual_heavy_rows", "actual_source_hits",
            "actual_dense_two_anchor", "actual_dense_two_anchor_phase",
            "precode_count",
        ):
            maximum = 1 if key == "heavy_shortfall" else None
            strict_uint(row[key], f"K={K}:{key}", maximum)
        staircase = int(row["actual_staircase_rows"])
        dense_rows = int(row["actual_dense_rows"])
        heavy_rows = int(row["actual_heavy_rows"])
        if (dense_rows != 12 or
                heavy_rows != ARM_GF256_ROWS[task.arm] + 2 or
                row["actual_dense_two_anchor"] != "1" or
                row["actual_dense_two_anchor_phase"] != "0" or
                int(row["actual_source_hits"]) != (2 if K < 10000 else 3) or
                int(row["precode_count"]) != staircase + dense_rows + heavy_rows):
            die(f"actual attempt-0 construction geometry mismatch at K={K}")
        first_rank = row["first_rank_fail"]
        if first_rank not in ("-1", "0") or (outcomes[1] == 1) != (first_rank == "0"):
            die(f"first-rank-failure field disagrees with outcome at K={K}")
        binary_hist_value = validate_histogram(
            row["binary_def_hist"], f"K={K}:binary histogram",
        )
        heavy_hist_value = validate_histogram(
            row["heavy_gain_hist"], f"K={K}:heavy histogram",
        )
        inact = int(row["inact_max"])
        binary_def = int(row["binary_def_max"])
        heavy_gain = int(row["heavy_gain_min"])
        if (fixed_values["inact_mu"] != inact or
                fixed_values["binary_def_mu"] != binary_def or
                fixed_values["heavy_gain_mu"] != heavy_gain or
                binary_hist_value != binary_def or
                heavy_hist_value != heavy_gain):
            die(f"one-trial scalar/histogram metrics disagree at K={K}")
        expected_shortfall = int(
            outcomes[1] == 1 and binary_def <= heavy_rows and
            heavy_gain < binary_def
        )
        if (outcomes[1] and heavy_gain >= binary_def) or \
                int(row["heavy_shortfall"]) != expected_shortfall:
            die(f"rank-failure cause metrics disagree at K={K}")
        expected_failure_trials = "" if outcomes[0] else "0"
        if row["failure_trials"] != expected_failure_trials:
            die(f"failure-trial list disagrees with outcome at K={K}")
        base_matrix = strict_hex(row["base_matrix_seed"], "base matrix", 64)
        base_peel = strict_hex(row["base_peel_seed"], "base peel", 32)
        if (base_matrix != task.construction_seed or
                base_peel != construction_peel_seed(task.construction_seed) or
                strict_hex(row["matrix_seed"], "matrix", 64) != base_matrix or
                strict_hex(row["peel_seed"], "peel", 32) != base_peel):
            die(f"uniform construction-root seed mapping mismatch at K={K}")
        if strict_uint(row["systematic_probe_result"], "systematic probe", 1) not in (0, 1):
            die(f"unexpected systematic probe result at K={K}")
        expected_loss_seed = loss_seed(task.seed, K)
        if strict_hex(row["packet_trace_seed"], "packet trace seed", 64) != expected_loss_seed:
            die(f"packet trace seed formula mismatch at K={K}")
        if not re.fullmatch(r"[0-9a-f]{64}", row["packet_trace_sha256"]):
            die(f"malformed packet trace SHA256 at K={K}")
        rows.append(row)
    return rows


def row_failed(row: Mapping[str, str]) -> bool:
    return row["rank_fail"] == "1" or row["error"] == "1"


def pair_rows(
    tasks_and_rows: Iterable[tuple[Task, Sequence[dict[str, str]]]],
) -> dict[tuple[int, int, str, int], dict[str, dict[str, str]]]:
    cells: dict[tuple[int, int, str, int], dict[str, dict[str, str]]] = {}
    for task, rows in tasks_and_rows:
        if len(rows) != len(task.ks):
            die("parsed task row count mismatch")
        for K, row in zip(task.ks, rows):
            key = (task.construction_index, task.seed_index, task.schedule, K)
            arms = cells.setdefault(key, {})
            if task.arm in arms:
                die(f"duplicate arm result for cell {key}")
            arms[task.arm] = row
    for key, arms in cells.items():
        validate_pair_receipt(key, arms)
    return cells


def validate_pair_receipt(
    key: tuple[int, int, str, int], arms: Mapping[str, Mapping[str, str]],
) -> None:
    if set(arms) != set(ARMS):
        die(f"unpaired Stage-A cell {key}")
    h12, h13 = arms["h12"], arms["h13"]
    for field in (
        "base_matrix_seed", "base_peel_seed", "matrix_seed", "peel_seed",
        "construction_seed_policy", "construction_seed",
        "production_seed_fixups_applied",
        "actual_staircase_rows", "actual_dense_rows", "actual_source_hits",
        "actual_dense_two_anchor", "actual_dense_two_anchor_phase",
        "packet_trace_seed", "packet_trace_sha256",
    ):
        if h12[field] != h13[field]:
            die(f"paired raw receipt {field} differs at {key}")
    if strict_uint(h13["precode_count"], "H13 precode count") != \
            strict_uint(h12["precode_count"], "H12 precode count") + 1:
        die(f"H13 precode count is not H12 + 1 at {key}")


def exact_mcnemar(repairs: int, introductions: int) -> str:
    discordant = repairs + introductions
    if discordant == 0:
        return "1"
    lower = min(repairs, introductions)
    numerator = sum(math.comb(discordant, index) for index in range(lower + 1))
    value = min(Decimal(1), Decimal(2 * numerator) / Decimal(1 << discordant))
    return str(value)


def _milli(row: Mapping[str, str], field: str) -> int:
    return int(Decimal(row[field]) * 1000)


def _ratio(numerator: int, denominator: int) -> str | None:
    return None if denominator == 0 else str(Decimal(numerator) / Decimal(denominator))


def _empty_scope() -> dict[str, Any]:
    return {
        "paired_cells": 0,
        "arms": {
            arm: {
                "success": 0, "rank_fail": 0, "error": 0,
                "heavy_shortfall_sum": 0, "inact_milli_sum": 0,
                "binary_def_milli_sum": 0, "heavy_gain_milli_sum": 0,
                **{f"{field}_milli_sum": 0 for field in WORK_FIELDS},
            } for arm in ARMS
        },
        "comparison": {
            "repairs": 0, "introductions": 0, "both_fail": 0,
            "common_success": 0,
            "common_success_work": {
                arm: {f"{field}_milli_sum": 0 for field in WORK_FIELDS}
                for arm in ARMS
            },
        },
    }


def _update_scope(
    scope: dict[str, Any], pair: Mapping[str, Mapping[str, str]],
) -> None:
    scope["paired_cells"] += 1
    base_fail = row_failed(pair["h12"])
    candidate_fail = row_failed(pair["h13"])
    comparison = scope["comparison"]
    comparison["repairs"] += base_fail and not candidate_fail
    comparison["introductions"] += not base_fail and candidate_fail
    comparison["both_fail"] += base_fail and candidate_fail
    comparison["common_success"] += not base_fail and not candidate_fail
    for arm in ARMS:
        row = pair[arm]
        arm_scope = scope["arms"][arm]
        for field in ("success", "rank_fail", "error"):
            arm_scope[field] += int(row[field])
        arm_scope["heavy_shortfall_sum"] += int(row["heavy_shortfall"])
        for field in ("inact_mu", "binary_def_mu", "heavy_gain_mu"):
            arm_scope[f"{field[:-3]}_milli_sum"] += _milli(row, field)
        for field in WORK_FIELDS:
            arm_scope[f"{field}_milli_sum"] += _milli(row, field)
            if not base_fail and not candidate_fail:
                comparison["common_success_work"][arm][
                    f"{field}_milli_sum"
                ] += _milli(row, field)


def _finish_scope(scope: dict[str, Any]) -> dict[str, Any]:
    comparison = scope["comparison"]
    comparison["exact_two_sided_mcnemar"] = exact_mcnemar(
        comparison["repairs"], comparison["introductions"],
    )
    comparison["h13_over_h12_common_success_work_ratios"] = {
        field: _ratio(
            comparison["common_success_work"]["h13"][f"{field}_milli_sum"],
            comparison["common_success_work"]["h12"][f"{field}_milli_sum"],
        ) for field in WORK_FIELDS
    }
    return scope


def aggregate_pair_stream(
    pairs: Iterable[
        tuple[tuple[int, int, str, int], Mapping[str, Mapping[str, str]]]
    ],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    overall = _empty_scope()
    by_construction = {
        str(i): _empty_scope() for i in range(len(CONSTRUCTION_SEEDS))
    }
    by_seed = {str(i): _empty_scope() for i in range(len(LOSS_ROOTS))}
    by_schedule = {name: _empty_scope() for name in SCHEDULES}
    by_band = {name: _empty_scope() for name, _, _ in K_BANDS}
    by_joint = {
        f"construction{construction_index}:loss{seed_index}:{schedule}":
            _empty_scope()
        for construction_index in range(len(CONSTRUCTION_SEEDS))
        for seed_index in range(len(LOSS_ROOTS)) for schedule in SCHEDULES
    }
    failure_by_k: dict[str, dict[int, int]] = {arm: {} for arm in ARMS}
    failure_by_kc: dict[str, dict[tuple[int, int], int]] = {
        arm: {} for arm in ARMS
    }
    failure_records_kc: dict[
        str, dict[tuple[int, int], list[dict[str, Any]]]
    ] = {
        arm: {} for arm in ARMS
    }
    sampled_ks: set[int] = set()
    sampled_kcs: set[tuple[int, int]] = set()
    union_failures: list[dict[str, Any]] = []
    previous_ordinal = -1
    for (construction_index, seed_index, schedule, K), pair in pairs:
        ordinal = ((((construction_index * len(LOSS_ROOTS)) + seed_index) *
                    len(SCHEDULES) + SCHEDULES.index(schedule)) * K_COUNT +
                   K - K_MIN)
        if ordinal <= previous_ordinal:
            die("stage cells are duplicated or not in canonical order")
        previous_ordinal = ordinal
        sampled_ks.add(K)
        sampled_kcs.add((K, construction_index))
        _update_scope(overall, pair)
        _update_scope(by_construction[str(construction_index)], pair)
        _update_scope(by_seed[str(seed_index)], pair)
        _update_scope(by_schedule[schedule], pair)
        band_names = [
            name for name, first, last in K_BANDS if first <= K <= last
        ]
        if len(band_names) != 1:
            die("K band partition does not cover each sampled K exactly once")
        _update_scope(by_band[band_names[0]], pair)
        _update_scope(by_joint[
            f"construction{construction_index}:loss{seed_index}:{schedule}"
        ], pair)
        any_failure = False
        for arm in ARMS:
            row = pair[arm]
            failed = row_failed(row)
            any_failure |= failed
            failure_by_k[arm][K] = failure_by_k[arm].get(K, 0) + failed
            kc = (K, construction_index)
            failure_by_kc[arm][kc] = failure_by_kc[arm].get(kc, 0) + failed
            if failed:
                failure_records_kc[arm].setdefault(kc, []).append({
                    "construction_seed_index": construction_index,
                    "construction_seed": str(
                        CONSTRUCTION_SEEDS[construction_index]),
                    "packet_trace_root_index": seed_index,
                    "packet_trace_root": LOSS_ROOTS[seed_index],
                    "schedule": schedule,
                    "cause": "error" if row["error"] == "1" else "rank_fail",
                    "systematic_probe_result": int(row["systematic_probe_result"]),
                    "inact_mu": row["inact_mu"],
                    "binary_def_mu": row["binary_def_mu"],
                    "heavy_gain_mu": row["heavy_gain_mu"],
                    "heavy_shortfall": int(row["heavy_shortfall"]),
                    "work": {field: row[field] for field in WORK_FIELDS},
                    "loss_seed": row["packet_trace_seed"],
                    "packet_trace_sha256": row["packet_trace_sha256"],
                    "base_matrix_seed": row["base_matrix_seed"],
                    "base_peel_seed": row["base_peel_seed"],
                })
        if any_failure:
            union_failures.append({
                "construction_seed_index": construction_index,
                "construction_seed": str(CONSTRUCTION_SEEDS[construction_index]),
                "packet_trace_root_index": seed_index,
                "packet_trace_root": LOSS_ROOTS[seed_index],
                "schedule": schedule, "K": K,
            })

    arms = overall["arms"]
    for arm in ARMS:
        kc_histogram = {str(count): 0 for count in range(10)}
        for kc in sampled_kcs:
            count = failure_by_kc[arm].get(kc, 0)
            if count > 9:
                die("one (K,construction) has more than nine failures")
            kc_histogram[str(count)] += 1
        k_histogram = {str(count): 0 for count in range(28)}
        for K in sampled_ks:
            count = failure_by_k[arm].get(K, 0)
            if count > 27:
                die("one K has more than 27 named-stratum failures")
            k_histogram[str(count)] += 1
        arms[arm]["sampled_unique_K"] = len(sampled_ks)
        arms[arm]["sampled_K_construction_pairs"] = len(sampled_kcs)
        arms[arm]["weak_K_construction_pairs"] = (
            len(sampled_kcs) - kc_histogram["0"]
        )
        arms[arm]["weak_K"] = len(sampled_ks) - k_histogram["0"]
        arms[arm]["multi_failure_K_construction_pairs"] = sum(
            kc_histogram[str(index)] for index in range(2, 10)
        )
        arms[arm]["maximum_failure_strata_per_K_construction"] = max(
            failure_by_kc[arm].values(), default=0,
        )
        arms[arm]["maximum_failure_strata_per_K"] = max(
            failure_by_k[arm].values(), default=0,
        )
        arms[arm]["failure_strata_histogram_per_K_construction_0_to_9"] = (
            kc_histogram
        )
        arms[arm]["failure_strata_histogram_per_K_0_to_27"] = k_histogram
        arms[arm]["weak_K_construction_records"] = [
            {"K": K, "construction_seed_index": construction_index,
             "construction_seed": str(CONSTRUCTION_SEEDS[construction_index]),
             "failure_strata": len(records), "failures": records}
            for (K, construction_index), records in sorted(
                failure_records_kc[arm].items())
        ]
        arms[arm]["weak_K_records"] = [{
            "K": K,
            "failure_strata": failure_by_k[arm][K],
            "construction_seed_multiplicities": [{
                "construction_seed_index": construction_index,
                "construction_seed": str(CONSTRUCTION_SEEDS[construction_index]),
                "failure_strata": failure_by_kc[arm].get(
                    (K, construction_index), 0),
            } for construction_index in range(len(CONSTRUCTION_SEEDS))
              if (K, construction_index) in sampled_kcs],
        } for K in sorted(failure_by_k[arm]) if failure_by_k[arm][K] > 0]
    report = _finish_scope(overall)
    report["descriptive_scopes"] = {
        "by_construction_seed": {
            key: _finish_scope(value) for key, value in by_construction.items()
            if value["paired_cells"]
        },
        "by_packet_trace_root": {
            key: _finish_scope(value) for key, value in by_seed.items()
            if value["paired_cells"]
        },
        "by_schedule": {
            key: _finish_scope(value) for key, value in by_schedule.items()
            if value["paired_cells"]
        },
        "by_K_band": {
            key: _finish_scope(value) for key, value in by_band.items()
            if value["paired_cells"]
        },
        "by_construction_loss_schedule": {
            key: _finish_scope(value) for key, value in by_joint.items()
            if value["paired_cells"]
        },
    }
    cluster_differences = [
        failure_by_k["h12"].get(K, 0) - failure_by_k["h13"].get(K, 0)
        for K in sorted(sampled_ks)
    ]
    cluster_count = len(cluster_differences)
    if cluster_count > 1:
        mean = Decimal(sum(cluster_differences)) / Decimal(cluster_count)
        squared = sum((Decimal(value) - mean) ** 2
                      for value in cluster_differences)
        standard_error = (
            squared / Decimal(cluster_count * (cluster_count - 1))
        ).sqrt()
        lower = mean - Decimal("1.96") * standard_error
        upper = mean + Decimal("1.96") * standard_error
    else:
        mean = Decimal(cluster_differences[0]) if cluster_differences else Decimal(0)
        standard_error = None
        lower = upper = None
    report["comparison"]["K_cluster_descriptive"] = {
        "clusters": cluster_count,
        "estimand": "mean per-K (H12 failures minus H13 failures)",
        "mean": str(mean),
        "delete_one_jackknife_standard_error": (
            None if standard_error is None else str(standard_error)
        ),
        "normal_95_interval": (
            None if lower is None else [str(lower), str(upper)]
        ),
        "interpretation": (
            "descriptive sensitivity across the complete ordered K domain; "
            "K clusters are not asserted IID and this interval is not a "
            "population-level acceptance gate"
        ),
    }
    report["comparison"]["union_failures"] = len(union_failures)
    return report, union_failures


def aggregate_cells(
    cells: Mapping[tuple[int, int, str, int], Mapping[str, Mapping[str, str]]],
) -> dict[str, Any]:
    return aggregate_pair_stream(cells.items())[0]


def campaign_round(construction_index: int) -> str:
    if construction_index == 0:
        return "R0"
    if construction_index in (1, 2):
        return "R1"
    die("construction index has no campaign round")


def manifest_payload(
    construction_index: int, overhead: int, tasks: Sequence[Task],
) -> dict[str, Any]:
    if (type(overhead) is not int or not 0 <= overhead <= MAX_OVERHEAD or
            any(task.construction_index != construction_index or
                task.overhead != overhead for task in tasks)):
        die("cannot seal a manifest with mismatched task geometry")
    return {
        "campaign_round": campaign_round(construction_index),
        "construction_seed_index": construction_index,
        "construction_seed": str(CONSTRUCTION_SEEDS[construction_index]),
        "overhead": overhead, "chunk_size": CHUNK_SIZE,
        "paired_overhead_stream": True, "tasks": [task.payload() for task in tasks],
    }


def load_manifest(
    stage: Path, construction_index: int, overhead: int,
) -> list[Task]:
    payload = load_sealed(
        stage / "manifest.json", f"{SCHEMA}.stage_manifest",
    )
    if (not isinstance(payload, dict) or set(payload) != {
            "campaign_round", "construction_seed_index", "construction_seed",
            "overhead", "chunk_size", "paired_overhead_stream", "tasks"} or
            payload["campaign_round"] != campaign_round(construction_index) or
            type(payload["construction_seed_index"]) is not int or
            payload["construction_seed_index"] != construction_index or
            type(payload["construction_seed"]) is not str or
            payload["construction_seed"] != str(
                CONSTRUCTION_SEEDS[construction_index]) or
            type(payload["overhead"]) is not int or
            type(payload["chunk_size"]) is not int or
            payload["overhead"] != overhead or payload["chunk_size"] != CHUNK_SIZE or
            payload["paired_overhead_stream"] is not True or
            not isinstance(payload["tasks"], list)):
        die("stage manifest contract mismatch")
    tasks = [task_from_payload(value) for value in payload["tasks"]]
    if any(task.construction_index != construction_index for task in tasks):
        die("stage manifest task construction root mismatch")
    if [task.job for task in tasks] != list(range(len(tasks))):
        die("stage task IDs are not contiguous unbounded integers")
    if (overhead == 0 and
            (len(tasks) != OH0_JOBS or
             sum(len(task.ks) for task in tasks) != TOTAL_ARM_OUTCOMES)):
        die("OH0 manifest does not contain the exact full campaign census")
    if overhead == 0:
        expected_tasks = build_tasks(0, construction_index)
    else:
        if not tasks:
            die("nested stage manifest may not encode an empty cohort")
        cohort = {
            (task.seed_index, task.schedule, K)
            for task in tasks[::2] for K in task.ks
        }
        expected_tasks = build_tasks(overhead, construction_index, cohort)
    if tasks != expected_tasks:
        die("stage manifest is not the canonical planner output")
    return tasks


def construction_root_path(result_dir: Path, construction_index: int) -> Path:
    if (type(construction_index) is not int or
            construction_index not in range(len(CONSTRUCTION_SEEDS))):
        die("construction root index outside sealed domain")
    return result_dir / "roots" / f"c{construction_index}"


def stage_path(
    result_dir: Path, construction_index: int, overhead: int,
) -> Path:
    return (construction_root_path(result_dir, construction_index) /
            "stages" / f"oh{overhead:04d}")


def shard_path(stage: Path, task: Task) -> Path:
    return stage / "shards" / task.stem


def receipt_payload(
    task: Task, benchmark_argv: Sequence[str], command: Sequence[str], cpu: int,
    binary_sha256: str, stdout: bytes, stderr: bytes, row_count: int,
) -> dict[str, Any]:
    return {
        "task": task.payload(), "benchmark_argv": list(benchmark_argv),
        "command": list(command), "cpu": cpu, "binary_sha256": binary_sha256,
        "execution_mode": EXECUTION_MODE,
        "returncode": 0, "stdout_sha256": sha256_bytes(stdout),
        "stderr_sha256": sha256_bytes(stderr), "row_count": row_count,
        "loss_seed_formula": LOSS_SEED_FORMULA,
        "loss_seeds": [{"K": K, "hex": hex(loss_seed(task.seed, K))} for K in task.ks],
    }


def validate_shard(
    stage: Path, task: Task, binary: Path, binary_sha256: str, taskset: Path,
    allowed_cpus: set[int],
) -> list[dict[str, str]]:
    root = shard_path(stage, task)
    if not root.is_dir() or root.is_symlink():
        die(f"missing sealed shard for {task.task_id}")
    inventory = {path.name: path for path in root.iterdir()}
    if (set(inventory) != {"stdout.csv", "stderr.txt", "receipt.json"} or
            any(path.is_symlink() or not path.is_file()
                for path in inventory.values())):
        die(f"sealed shard inventory mismatch for {task.task_id}")
    stdout = stable_bytes(root / "stdout.csv")
    stderr = stable_bytes(root / "stderr.txt")
    receipt = load_sealed(root / "receipt.json", f"{SCHEMA}.job_receipt")
    if not isinstance(receipt, dict):
        die("job receipt payload is not an object")
    expected_argv = make_benchmark_argv(binary, task)
    cpu = receipt.get("cpu")
    if type(cpu) is not int or cpu == 127 or cpu not in allowed_cpus:
        die(f"job receipt CPU is outside the sealed controller affinity for {task.task_id}")
    expected_command = [str(taskset), "-c", str(cpu), *expected_argv]
    expected_receipt = receipt_payload(
        task, expected_argv, expected_command, cpu, binary_sha256,
        stdout, stderr, len(task.ks),
    )
    if stderr or canonical_json(receipt) != canonical_json(expected_receipt):
        die(f"job receipt substitution or content mismatch for {task.task_id}")
    try:
        text = stdout.decode("ascii")
    except UnicodeDecodeError as exc:
        die(f"benchmark output is not ASCII: {exc}")
    return parse_output(text, task)


def _stat_identity(value: os.stat_result) -> tuple[int, ...]:
    return (
        value.st_dev, value.st_ino, value.st_mode, value.st_nlink,
        value.st_size, value.st_mtime_ns, value.st_ctime_ns,
    )


def _stable_fd_sha256(descriptor: int, size: int) -> str:
    digests: list[str] = []
    for _ in range(2):
        digest = hashlib.sha256()
        offset = 0
        while offset < size:
            block = os.pread(descriptor, min(1 << 20, size - offset), offset)
            if not block:
                die("authenticated executable ended before its sealed size")
            digest.update(block)
            offset += len(block)
        if os.pread(descriptor, 1, size):
            die("authenticated executable grew while hashing")
        digests.append(digest.hexdigest())
    if digests[0] != digests[1]:
        die("authenticated executable changed between hash passes")
    return digests[0]


def open_authenticated_executable(path: Path, expected_sha256: str) -> int:
    """Open one immutable frozen executable and bind its bytes to an fd."""
    if not re.fullmatch(r"[0-9a-f]{64}", expected_sha256):
        die("authenticated executable has a malformed expected SHA256")
    flags = (os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
             getattr(os, "O_NOFOLLOW", 0))
    try:
        descriptor = os.open(path, flags)
    except OSError as exc:
        die(f"cannot open frozen executable {path}: {exc}")
    try:
        before = os.fstat(descriptor)
        if (not stat.S_ISREG(before.st_mode) or before.st_nlink != 1 or
                stat.S_IMODE(before.st_mode) != 0o555):
            die("frozen executable is not a unique mode-0555 regular file")
        digest = _stable_fd_sha256(descriptor, before.st_size)
        after = os.fstat(descriptor)
        try:
            named = os.stat(path, follow_symlinks=False)
        except OSError as exc:
            die(f"frozen executable pathname changed while opening: {exc}")
        if (_stat_identity(before) != _stat_identity(after) or
                _stat_identity(after) != _stat_identity(named) or
                digest != expected_sha256):
            die("frozen executable identity/hash changed while opening")
        return descriptor
    except BaseException:
        os.close(descriptor)
        raise


# The child is launched by the already-running interpreter and authenticates
# inherited taskset/benchmark fds immediately before exec.  No campaign tool is
# reopened through its mutable visible pathname after parent authentication.
AUTHENTICATED_FD_LAUNCHER = r'''import ctypes, hashlib, os, signal, stat, sys
def fail(message):
    os.write(2, ("authenticated launcher: " + message + "\n").encode("ascii", "replace"))
    os._exit(125)
def identity(value):
    return (value.st_dev, value.st_ino, value.st_mode, value.st_nlink,
            value.st_size, value.st_mtime_ns, value.st_ctime_ns)
def authenticate(descriptor, expected):
    before = os.fstat(descriptor)
    if (not stat.S_ISREG(before.st_mode) or before.st_nlink != 1 or
            stat.S_IMODE(before.st_mode) != 0o555):
        fail("executable fd is not a unique mode-0555 regular file")
    found = []
    for unused in range(2):
        digest = hashlib.sha256()
        offset = 0
        while offset < before.st_size:
            block = os.pread(descriptor, min(1048576, before.st_size - offset), offset)
            if not block:
                fail("executable fd ended before its sealed size")
            digest.update(block)
            offset += len(block)
        if os.pread(descriptor, 1, before.st_size):
            fail("executable fd grew while hashing")
        found.append(digest.hexdigest())
    after = os.fstat(descriptor)
    if identity(before) != identity(after) or found != [expected, expected]:
        fail("executable fd identity/hash mismatch")
try:
    if len(sys.argv) < 9 or sys.argv[0] != "-c":
        fail("malformed launcher arguments")
    parent = int(sys.argv[1]); cpu = int(sys.argv[2])
    taskset_fd = int(sys.argv[3]); taskset_sha256 = sys.argv[4]
    binary_fd = int(sys.argv[5]); binary_sha256 = sys.argv[6]
    taskset_argv0 = sys.argv[7]; benchmark_tail = sys.argv[8:]
    if parent <= 1 or cpu < 0 or taskset_fd <= 2 or binary_fd <= 2 or not benchmark_tail:
        fail("launcher scalar argument is outside its domain")
    if os.getppid() != parent:
        os.kill(os.getpid(), signal.SIGKILL)
    libc = ctypes.CDLL(None, use_errno=True)
    if libc.prctl(1, signal.SIGKILL, 0, 0, 0) != 0:
        fail("cannot arm PR_SET_PDEATHSIG")
    if os.getppid() != parent:
        os.kill(os.getpid(), signal.SIGKILL)
    os.sched_setaffinity(0, {cpu})
    if os.sched_getaffinity(0) != {cpu}:
        fail("child affinity did not bind to its exact assigned CPU")
    authenticate(taskset_fd, taskset_sha256)
    authenticate(binary_fd, binary_sha256)
    os.set_inheritable(binary_fd, True)
    command = [taskset_argv0, "-c", str(cpu),
               "/proc/self/fd/" + str(binary_fd)] + benchmark_tail
    os.execve(taskset_fd, command, os.environ)
except BaseException as exc:
    fail("launch failed: " + str(exc))
'''


def _run_process(
    command: Sequence[str], timeout: float, children: ActiveChildren, *,
    cpu: int, binary_fd: int, binary_sha256: str,
    taskset_fd: int, taskset_sha256: str,
) -> tuple[int, bytes, bytes]:
    children.check()
    if (sys.platform != "linux" or type(cpu) is not int or cpu < 0 or
            len(command) < 5 or command[1:3] != ["-c", str(cpu)] or
            any(type(value) is not str for value in command) or
            any(type(value) is not int or value <= 2
                for value in (binary_fd, taskset_fd))):
        die("authenticated benchmark launcher arguments are malformed")
    launcher = [
        "/proc/self/exe", "-I", "-S", "-c", AUTHENTICATED_FD_LAUNCHER,
        str(os.getpid()), str(cpu), str(taskset_fd), taskset_sha256,
        str(binary_fd), binary_sha256, command[0], *command[4:],
    ]
    process = subprocess.Popen(
        launcher, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        start_new_session=True, pass_fds=(taskset_fd, binary_fd),
    )
    children.register(process)
    try:
        try:
            stdout, stderr = process.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            ActiveChildren._kill(process)
            stdout, stderr = process.communicate()
            die(f"benchmark timed out after {timeout:g} seconds")
        except BaseException:
            ActiveChildren._kill(process)
            process.communicate()
            raise
    finally:
        children.unregister(process)
    children.check()
    if len(stdout) > OUTPUT_LIMIT or len(stderr) > OUTPUT_LIMIT:
        die("benchmark output exceeded the sealed capture limit")
    return process.returncode, stdout, stderr


def execute_task(
    result_dir: Path, stage: Path, task: Task, binary: Path,
    binary_sha256: str, taskset: Path, taskset_sha256: str,
    cpus: "queue.Queue[int]", timeout: float, children: ActiveChildren,
    allowed_cpus: set[int], binary_fd: int, taskset_fd: int,
) -> None:
    children.check()
    final = shard_path(stage, task)
    if final.exists():
        validate_shard(
            stage, task, binary, binary_sha256, taskset, allowed_cpus,
        )
        return
    cpu = cpus.get()
    try:
        if cpu == 127:
            die("CPU 127 is reserved for the thermal sampler")
        children.check()
        benchmark_argv = make_benchmark_argv(binary, task)
        command = [str(taskset), "-c", str(cpu), *benchmark_argv]
        returncode, stdout, stderr = _run_process(
            command, timeout, children, cpu=cpu, binary_fd=binary_fd,
            binary_sha256=binary_sha256, taskset_fd=taskset_fd,
            taskset_sha256=taskset_sha256,
        )
        if returncode != 0 or stderr:
            tail = stderr.decode("utf-8", errors="replace")[-1000:]
            die(f"job {task.job} failed rc={returncode}: {tail}")
        try:
            rows = parse_output(stdout.decode("ascii"), task)
        except UnicodeDecodeError as exc:
            die(f"job {task.job} emitted non-ASCII output: {exc}")
        children.check()
        staging = final.with_name(f".staging-{task.task_id}-{os.getpid()}")
        if staging.exists():
            if staging.is_symlink() or not staging.is_dir():
                die("orphan shard staging path has an unsafe type")
            shutil.rmtree(staging)
        staging.mkdir(parents=True)
        atomic_write(staging / "stdout.csv", stdout)
        atomic_write(staging / "stderr.txt", stderr)
        write_sealed_once(
            staging / "receipt.json", f"{SCHEMA}.job_receipt",
            receipt_payload(
                task, benchmark_argv, command, cpu, binary_sha256,
                stdout, stderr, len(rows),
            ),
        )
        children.check()
        final.parent.mkdir(parents=True, exist_ok=True)
        try:
            staging.rename(final)
        except FileExistsError:
            shutil.rmtree(staging)
            validate_shard(
                stage, task, binary, binary_sha256, taskset, allowed_cpus,
            )
    finally:
        cpus.put(cpu)


def iter_stage_pairs(
    stage: Path, tasks: Sequence[Task], binary: Path, binary_sha256: str,
    taskset: Path, allowed_cpus: set[int], *, progress: bool = False,
) -> Iterable[
    tuple[tuple[int, int, str, int], dict[str, dict[str, str]]]
]:
    if len(tasks) % 2:
        die("stage manifest has an odd task count")
    pair_count = len(tasks) // 2
    for pair_index in range(pair_count):
        left, right = tasks[2 * pair_index:2 * pair_index + 2]
        if (left.arm != "h12" or right.arm != "h13" or
                left.pair != pair_index or right.pair != pair_index or
                left.overhead != right.overhead or
                left.construction_index != right.construction_index or
                left.construction_seed != right.construction_seed or
                left.seed_index != right.seed_index or left.seed != right.seed or
                left.schedule != right.schedule or left.ks != right.ks):
            die("stage manifest arm pairing/order mismatch")
        left_rows = validate_shard(
            stage, left, binary, binary_sha256, taskset, allowed_cpus,
        )
        right_rows = validate_shard(
            stage, right, binary, binary_sha256, taskset, allowed_cpus,
        )
        for K, left_row, right_row in zip(left.ks, left_rows, right_rows):
            key = (left.construction_index, left.seed_index, left.schedule, K)
            arms = {"h12": left_row, "h13": right_row}
            validate_pair_receipt(key, arms)
            yield key, arms
        if progress and ((pair_index + 1) % 256 == 0 or pair_index + 1 == pair_count):
            print(
                f"# Stage-A OH{left.overhead}: verified paired shard "
                f"{pair_index + 1}/{pair_count}",
                file=sys.stderr, flush=True,
            )


def validate_stage_shard_inventory(stage: Path, tasks: Sequence[Task]) -> None:
    root = stage / "shards"
    if not root.is_dir() or root.is_symlink():
        die("stage shard root is missing, symlinked, or not a directory")
    expected = {task.stem for task in tasks}
    if len(expected) != len(tasks):
        die("stage manifest contains duplicate shard stems")
    entries = list(root.iterdir())
    if ({entry.name for entry in entries} != expected or
            any(not entry.is_dir() or entry.is_symlink() for entry in entries)):
        die("stage shard inventory differs from the sealed manifest")


def validate_stage_file_inventory(
    stage: Path, *, complete_required: bool,
) -> None:
    """Reject out-of-ledger files beside a stage manifest and shard set."""
    if not stage.is_dir() or stage.is_symlink():
        die("campaign stage is missing, symlinked, or not a directory")
    expected = {"manifest.json", "shards"}
    complete = stage / "complete.json"
    if complete_required or complete.exists():
        expected.add("complete.json")
    entries = {entry.name: entry for entry in stage.iterdir()}
    if set(entries) != expected:
        die("stage file inventory differs from the sealed campaign layout")
    if (not entries["manifest.json"].is_file() or
            entries["manifest.json"].is_symlink() or
            not entries["shards"].is_dir() or
            entries["shards"].is_symlink() or
            ("complete.json" in entries and
             (not entries["complete.json"].is_file() or
              entries["complete.json"].is_symlink()))):
        die("stage file inventory contains a wrong-type or symlinked entry")


def compute_stage_payload(
    stage: Path, overhead: int, tasks: Sequence[Task], binary: Path,
    binary_sha256: str, taskset: Path, allowed_cpus: set[int],
    *, progress: bool = False,
) -> dict[str, Any]:
    validate_stage_shard_inventory(stage, tasks)
    summary, failures = aggregate_pair_stream(iter_stage_pairs(
        stage, tasks, binary, binary_sha256, taskset, allowed_cpus,
        progress=progress,
    ))
    return {
        "overhead": overhead, "job_count": len(tasks),
        "paired_cell_count": summary["paired_cells"], "metrics": summary,
        "union_failures": failures,
    }


def complete_stage(
    stage: Path, overhead: int, tasks: Sequence[Task], binary: Path,
    binary_sha256: str, taskset: Path, allowed_cpus: set[int],
) -> dict[str, Any]:
    validate_stage_file_inventory(stage, complete_required=False)
    payload = compute_stage_payload(
        stage, overhead, tasks, binary, binary_sha256, taskset, allowed_cpus,
        progress=True,
    )
    write_sealed_once(stage / "complete.json", f"{SCHEMA}.stage_complete", payload)
    validate_stage_file_inventory(stage, complete_required=True)
    stored = load_sealed(stage / "complete.json", f"{SCHEMA}.stage_complete")
    if canonical_json(stored) != canonical_json(payload):
        die("stored stage completion differs from recomputed results")
    return payload


def verify_complete_stage(
    stage: Path, overhead: int, tasks: Sequence[Task], binary: Path,
    binary_sha256: str, taskset: Path, allowed_cpus: set[int],
) -> dict[str, Any]:
    validate_stage_file_inventory(stage, complete_required=True)
    stored = load_sealed(stage / "complete.json", f"{SCHEMA}.stage_complete")
    recomputed = compute_stage_payload(
        stage, overhead, tasks, binary, binary_sha256, taskset, allowed_cpus,
        progress=True,
    )
    if canonical_json(stored) != canonical_json(recomputed):
        die("sealed stage completion differs from recomputed shards")
    return stored


def parse_cpu_list(text: str) -> list[int]:
    cpus: set[int] = set()
    for token in text.split(","):
        if "-" in token:
            first_text, last_text = token.split("-", 1)
            first = strict_uint(first_text, "CPU")
            last = strict_uint(last_text, "CPU")
            if first > last:
                die("descending CPU range")
            cpus.update(range(first, last + 1))
        else:
            cpus.add(strict_uint(token, "CPU"))
    if not cpus or 127 in cpus:
        die("CPU list is empty or includes reserved CPU 127")
    available = os.sched_getaffinity(0)
    if not cpus.issubset(available):
        die("requested CPU is outside the controller affinity")
    return sorted(cpus)


def resolve_worker_count(requested: int | None, cpus: Sequence[int]) -> int:
    workers = len(cpus) if requested is None else requested
    if type(workers) is not int or workers < 1 or workers > len(cpus):
        die("workers must be in 1..len(cpus)")
    return workers


def exclusive_campaign(handler):
    @wraps(handler)
    def locked(args: argparse.Namespace):
        result_dir = args.result_dir.resolve(strict=True)
        lock_path = result_dir / "controller.lock"
        flags = (os.O_RDWR | os.O_CREAT | getattr(os, "O_CLOEXEC", 0) |
                 getattr(os, "O_NOFOLLOW", 0))
        try:
            descriptor = os.open(lock_path, flags, 0o600)
        except OSError as exc:
            die(f"cannot open campaign ownership lock: {exc}")
        try:
            metadata = os.fstat(descriptor)
            if not stat.S_ISREG(metadata.st_mode) or metadata.st_nlink != 1:
                die("campaign ownership lock is not a unique regular file")
            try:
                fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                die("another run/reduce controller owns this campaign")
            os.ftruncate(descriptor, 0)
            os.write(descriptor, f"pid={os.getpid()}\n".encode("ascii"))
            os.fsync(descriptor)
            return handler(args, result_dir)
        finally:
            os.close(descriptor)
    return locked


def _copy_unique(
    source: Path, destination: Path, executable: bool,
    *, allow_source_hardlinks: bool = False,
) -> str:
    if destination.exists():
        die(f"freeze target already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    source_data = stable_source_bytes(
        source, require_unique=not allow_source_hardlinks,
    )
    atomic_write(destination, source_data)
    os.chmod(destination, 0o555 if executable else 0o444)
    digest = sha256_bytes(source_data)
    if sha256_file(destination) != digest:
        die(f"frozen copy hash mismatch: {destination}")
    return digest


def prepare(args: argparse.Namespace) -> None:
    result_dir = args.result_dir.resolve()
    binary_source = args.binary.resolve(strict=True)
    controller_source = Path(__file__).resolve(strict=True)
    taskset_source_text = shutil.which("taskset")
    if not taskset_source_text:
        die("taskset is required")
    taskset_source = Path(taskset_source_text).resolve(strict=True)
    telemetry_path = (
        args.telemetry_log.resolve(strict=True) if args.telemetry_log else None
    )
    if (checked_product(K_COUNT, len(LOSS_ROOTS), len(SCHEDULES)) !=
            PAIRED_CELLS_PER_ROOT or
            checked_product(PAIRED_CELLS_PER_ROOT,
                            len(CONSTRUCTION_SEEDS)) != FULL_PAIRED_CELLS or
            checked_product(FULL_PAIRED_CELLS, len(ARMS)) !=
            FULL_ARM_OUTCOMES or any(seed == 0 for seed in CONSTRUCTION_SEEDS) or
            len(set(CONSTRUCTION_SEEDS)) != len(CONSTRUCTION_SEEDS)):
        die("campaign cardinality constants are inconsistent")
    telemetry_contract = {
        "managed_by_controller": False,
        "path": str(telemetry_path) if telemetry_path else None,
        "prepare_mark": external_file_mark(telemetry_path)
            if telemetry_path else None,
        "policy": TELEMETRY_POLICY,
        "sampler_identity": (
            "external operator responsibility; log inode and semantic interval "
            "are bound, sampler process identity is not authenticated"
        ),
        "continuity": TELEMETRY_CONTINUITY,
    }
    if result_dir.exists():
        die("result directory already exists")
    staging = result_dir.with_name(f".{result_dir.name}.prepare-{os.getpid()}")
    if staging.exists():
        die("prepare staging directory already exists")
    staging.mkdir(parents=True)
    try:
        frozen = staging / "frozen"
        frozen_binary = frozen / "wirehair_v2_bench"
        frozen_controller = frozen / Path(__file__).name
        frozen_taskset = frozen / "taskset"
        binary_sha256 = _copy_unique(binary_source, frozen_binary, True)
        controller_sha256 = _copy_unique(controller_source, frozen_controller, True)
        taskset_sha256 = _copy_unique(
            taskset_source, frozen_taskset, True, allow_source_hardlinks=True,
        )
        final_binary = result_dir / "frozen" / frozen_binary.name
        final_controller = result_dir / "frozen" / frozen_controller.name
        final_taskset = result_dir / "frozen" / frozen_taskset.name
        contract = {
            "K_domain": [K_MIN, K_MAX], "chunk_size": CHUNK_SIZE,
            "max_overhead": MAX_OVERHEAD, "arms": list(ARMS),
            "arm_gf256_rows": ARM_GF256_ROWS, "gf16_rows": 2,
            "period": 48, "geometry": "shared-x", "residue_schedule": "constant",
            "source_hits_rule": {"K_lt_10000": 2, "K_ge_10000": 3},
            "binary_dense_two_anchor": True, "binary_dense_two_anchor_phase": 0,
            "block_bytes": 64, "loss": "0.5", "trials": 1,
            "paired_overhead_stream": True, "raw_attempt": 0,
            "packet_trace_roots": list(LOSS_ROOTS),
            "packet_trace_seed_role": "loss-trace-root",
            "construction_seed_policy": CONSTRUCTION_SEED_POLICY,
            "construction_seeds": [str(seed) for seed in CONSTRUCTION_SEEDS],
            "construction_seed_assignment": (
                "same predeclared root for every K, loss root, schedule, "
                "overhead, and arm within one construction-root campaign"
            ),
            "base_matrix_seed_rule": "base_matrix_seed = construction_seed",
            "base_peel_seed_rule": (
                "base_peel_seed = low32(construction_seed) XOR "
                "high32(construction_seed)"
            ),
            "production_seed_fixups_applied": 0,
            "construction_retry_policy": "none; construction_attempt=0",
            "schedules": list(SCHEDULES),
            "K_bands": [
                {"name": name, "K_min": first, "K_max": last}
                for name, first, last in K_BANDS
            ],
            "rounds": [
                {"name": "R0", "construction_seed_indices": [0],
                 "paired_cells": R0_PAIRED_CELLS,
                 "gate": "H13 OH0 failures strictly fewer than H12 OH0 failures"},
                {"name": "R1", "construction_seed_indices": [1, 2],
                 "paired_cells": R1_PAIRED_CELLS,
                 "gate": "held-out pooled H13 OH0 failures strictly fewer than H12"},
            ],
            "round_execution_policy": (
                "R0 OH0 then fail-closed stop when gate fails; when it passes, "
                "finish the exact R0 union-only nested stream, run both R1 OH0 "
                "plans before either R1 nested stream, then finish each exact "
                "R1 union-only nested stream"
            ),
            "predeclaration_policy": (
                "all three construction/loss/schedule/K roots and all three "
                "OH0 manifests are sealed by prepare before any outcome"
            ),
            "paired_cells_per_construction_root": PAIRED_CELLS_PER_ROOT,
            "r0_paired_cells": R0_PAIRED_CELLS,
            "r1_paired_cells": R1_PAIRED_CELLS,
            "full_paired_cells": FULL_PAIRED_CELLS,
            "outcomes_per_arm_full": FULL_PAIRED_CELLS,
            "full_arm_outcomes_at_OH0": FULL_ARM_OUTCOMES,
            "oh0_jobs_per_construction_root": OH0_JOBS_PER_ROOT,
            "full_oh0_jobs": FULL_OH0_JOBS,
            "execution_mode": EXECUTION_MODE,
            "child_lifetime_policy": (
                "dedicated process group; registered kill/reap on SIGINT or "
                "SIGTERM; Linux PR_SET_PDEATHSIG=SIGKILL; inherited "
                "authenticated taskset/benchmark fds rehashed in child"
            ),
            "binary": str(final_binary), "binary_sha256": binary_sha256,
            "controller": str(final_controller),
            "controller_sha256": controller_sha256,
            "taskset": str(final_taskset), "taskset_sha256": taskset_sha256,
            "result_dir": str(result_dir), "loss_seed_formula": LOSS_SEED_FORMULA,
            "external_telemetry": telemetry_contract,
        }
        write_sealed_once(
            staging / "contract.json", f"{SCHEMA}.contract", contract,
        )
        for construction_index in range(len(CONSTRUCTION_SEEDS)):
            tasks = build_tasks(0, construction_index)
            initial_stage = (
                staging / "roots" / f"c{construction_index}" /
                "stages" / "oh0000"
            )
            write_sealed_once(
                initial_stage / "manifest.json", f"{SCHEMA}.stage_manifest",
                manifest_payload(construction_index, 0, tasks),
            )
        staging.rename(result_dir)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    record = {
        "result_dir": str(result_dir),
        "oh0_jobs_per_construction_root": OH0_JOBS_PER_ROOT,
        "full_oh0_jobs": FULL_OH0_JOBS,
        "r0_paired_cells": R0_PAIRED_CELLS,
        "r1_paired_cells": R1_PAIRED_CELLS,
        "full_paired_cells": FULL_PAIRED_CELLS,
        "telemetry_continuity": TELEMETRY_CONTINUITY,
        "run_command": [
            str(result_dir / "frozen" / Path(__file__).name), "run",
            "--result-dir", str(result_dir), "--cpus", "<CPU-LIST>",
        ],
        "reduce_command": [
            str(result_dir / "frozen" / Path(__file__).name), "reduce",
            "--result-dir", str(result_dir),
        ],
    }
    print(canonical_json(record))


def load_contract(result_dir: Path, *, require_frozen_controller: bool) -> dict[str, Any]:
    contract = load_sealed(result_dir / "contract.json", f"{SCHEMA}.contract")
    fixed = {
        "K_domain": [K_MIN, K_MAX], "chunk_size": CHUNK_SIZE,
        "max_overhead": MAX_OVERHEAD, "arms": list(ARMS),
        "arm_gf256_rows": ARM_GF256_ROWS, "gf16_rows": 2, "period": 48,
        "geometry": "shared-x", "residue_schedule": "constant",
        "source_hits_rule": {"K_lt_10000": 2, "K_ge_10000": 3},
        "binary_dense_two_anchor": True, "binary_dense_two_anchor_phase": 0,
        "block_bytes": 64, "loss": "0.5", "trials": 1,
        "paired_overhead_stream": True, "raw_attempt": 0,
        "packet_trace_roots": list(LOSS_ROOTS),
        "packet_trace_seed_role": "loss-trace-root",
        "construction_seed_policy": CONSTRUCTION_SEED_POLICY,
        "construction_seeds": [str(seed) for seed in CONSTRUCTION_SEEDS],
        "construction_seed_assignment": (
            "same predeclared root for every K, loss root, schedule, "
            "overhead, and arm within one construction-root campaign"
        ),
        "base_matrix_seed_rule": "base_matrix_seed = construction_seed",
        "base_peel_seed_rule": (
            "base_peel_seed = low32(construction_seed) XOR "
            "high32(construction_seed)"
        ),
        "production_seed_fixups_applied": 0,
        "construction_retry_policy": "none; construction_attempt=0",
        "schedules": list(SCHEDULES),
        "K_bands": [
            {"name": name, "K_min": first, "K_max": last}
            for name, first, last in K_BANDS
        ],
        "rounds": [
            {"name": "R0", "construction_seed_indices": [0],
             "paired_cells": R0_PAIRED_CELLS,
             "gate": "H13 OH0 failures strictly fewer than H12 OH0 failures"},
            {"name": "R1", "construction_seed_indices": [1, 2],
             "paired_cells": R1_PAIRED_CELLS,
             "gate": "held-out pooled H13 OH0 failures strictly fewer than H12"},
        ],
        "round_execution_policy": (
            "R0 OH0 then fail-closed stop when gate fails; when it passes, "
            "finish the exact R0 union-only nested stream, run both R1 OH0 "
            "plans before either R1 nested stream, then finish each exact "
            "R1 union-only nested stream"
        ),
        "predeclaration_policy": (
            "all three construction/loss/schedule/K roots and all three "
            "OH0 manifests are sealed by prepare before any outcome"
        ),
        "paired_cells_per_construction_root": PAIRED_CELLS_PER_ROOT,
        "r0_paired_cells": R0_PAIRED_CELLS,
        "r1_paired_cells": R1_PAIRED_CELLS,
        "full_paired_cells": FULL_PAIRED_CELLS,
        "outcomes_per_arm_full": FULL_PAIRED_CELLS,
        "full_arm_outcomes_at_OH0": FULL_ARM_OUTCOMES,
        "oh0_jobs_per_construction_root": OH0_JOBS_PER_ROOT,
        "full_oh0_jobs": FULL_OH0_JOBS,
        "loss_seed_formula": LOSS_SEED_FORMULA,
        "execution_mode": EXECUTION_MODE,
        "child_lifetime_policy": (
            "dedicated process group; registered kill/reap on SIGINT or "
            "SIGTERM; Linux PR_SET_PDEATHSIG=SIGKILL; inherited authenticated "
            "taskset/benchmark fds rehashed in child"
        ),
    }
    dynamic_keys = {
        "binary", "binary_sha256", "controller", "controller_sha256",
        "taskset", "taskset_sha256", "result_dir", "external_telemetry",
    }
    if (not isinstance(contract, dict) or
            set(contract) != set(fixed) | dynamic_keys or
            contract.get("result_dir") != str(result_dir) or
            any(type(contract.get(key)) is not str for key in (
                "binary", "binary_sha256", "controller", "controller_sha256",
                "taskset", "taskset_sha256", "result_dir",
            )) or
            any(not re.fullmatch(r"[0-9a-f]{64}", contract[key]) for key in (
                "binary_sha256", "controller_sha256", "taskset_sha256",
            ))):
        die("campaign contract schema/result directory mismatch")
    if any(
            canonical_json(contract.get(key)) != canonical_json(value)
            for key, value in fixed.items()):
        die("campaign contract geometry mismatch")
    telemetry = contract.get("external_telemetry")
    if (not isinstance(telemetry, dict) or set(telemetry) != {
            "managed_by_controller", "path", "prepare_mark", "policy",
            "sampler_identity", "continuity"} or
            telemetry["managed_by_controller"] is not False):
        die("external telemetry contract is malformed")
    if telemetry["path"] is None:
        if telemetry["prepare_mark"] is not None:
            die("external telemetry null path has a prepare mark")
    elif (not isinstance(telemetry["path"], str) or
            not isinstance(telemetry["prepare_mark"], dict) or
            telemetry["prepare_mark"].get("path") != telemetry["path"]):
        die("external telemetry path/prepare mark mismatch")
    if telemetry["path"] is not None:
        validate_external_mark(
            telemetry["prepare_mark"], telemetry["path"],
            "campaign contract",
        )
    if (canonical_json(telemetry["policy"]) != canonical_json(TELEMETRY_POLICY) or
            telemetry["sampler_identity"] !=
            "external operator responsibility; log inode and semantic interval "
            "are bound, sampler process identity is not authenticated" or
            telemetry["continuity"] != TELEMETRY_CONTINUITY):
        die("external telemetry policy/identity statement mismatch")
    for name in ("binary", "controller", "taskset"):
        path = Path(contract[name])
        if path.parent != result_dir / "frozen":
            die(f"frozen {name} substitution detected")
        descriptor = open_authenticated_executable(
            path, contract[f"{name}_sha256"],
        )
        os.close(descriptor)
    if require_frozen_controller and Path(__file__).resolve() != Path(contract["controller"]):
        die("run/reduce must use the frozen controller path")
    return contract


def telemetry_start(result_dir: Path, contract: Mapping[str, Any]) -> dict[str, Any] | None:
    telemetry = contract["external_telemetry"]
    if telemetry["path"] is None:
        return None
    path = Path(telemetry["path"])
    current = external_file_mark(path)
    prepare_mark = telemetry["prepare_mark"]
    if ((current["dev"], current["ino"]) !=
            (prepare_mark["dev"], prepare_mark["ino"]) or
            current["offset"] < prepare_mark["offset"]):
        die("external telemetry rotated or truncated after prepare")
    receipt_path = result_dir / "telemetry_start.json"
    if receipt_path.exists():
        start = load_sealed(receipt_path, f"{SCHEMA}.telemetry_start")
    else:
        start = {**current, "baseline": telemetry_baseline(
            path, current["offset"], (current["dev"], current["ino"]),
        )}
        write_sealed_once(receipt_path, f"{SCHEMA}.telemetry_start", start)
    confirmation = external_file_mark(path)
    if (not isinstance(start, dict) or set(start) != {
            "path", "dev", "ino", "offset", "mtime_ns", "baseline"} or
            not isinstance(start["baseline"], dict) or
            any(type(start[key]) is not int for key in (
                "dev", "ino", "offset", "mtime_ns")) or
            start["path"] != str(path) or
            (confirmation["dev"], confirmation["ino"]) !=
            (start["dev"], start["ino"]) or
            confirmation["offset"] < start["offset"] or
            canonical_json(start["baseline"]) != canonical_json(
                telemetry_baseline(
                    path, start["offset"], (start["dev"], start["ino"])))):
        die("external telemetry rotated or truncated after campaign start")
    return start


def telemetry_finish(
    result_dir: Path, contract: Mapping[str, Any], start: Mapping[str, Any] | None,
) -> dict[str, Any] | None:
    if start is None:
        return None
    path = Path(contract["external_telemetry"]["path"])
    receipt_path = result_dir / "telemetry_end.json"
    if receipt_path.exists():
        end = load_sealed(receipt_path, f"{SCHEMA}.telemetry_end")
    else:
        mark = external_file_mark(path)
        if ((mark["dev"], mark["ino"]) != (start["dev"], start["ino"]) or
                mark["offset"] < start["offset"]):
            die("external telemetry rotated or truncated during campaign")
        audit = audit_telemetry_interval(path, start, mark)
        confirmation = external_file_mark(path)
        if ((confirmation["dev"], confirmation["ino"]) !=
                (mark["dev"], mark["ino"]) or
                confirmation["offset"] < mark["offset"]):
            die("external telemetry rotated or truncated while sealing")
        end = {**mark, "start_offset": start["offset"],
               "interval_sha256": audit["interval_sha256"],
               "semantic_audit": audit}
        write_sealed_once(receipt_path, f"{SCHEMA}.telemetry_end", end)
    if (not isinstance(end, dict) or set(end) != {
            "path", "dev", "ino", "offset", "mtime_ns", "start_offset",
            "interval_sha256", "semantic_audit"} or
            any(type(end[key]) is not int for key in (
                "dev", "ino", "offset", "mtime_ns", "start_offset")) or
            end["path"] != str(path) or
            not isinstance(end["semantic_audit"], dict) or
            not isinstance(end["interval_sha256"], str) or
            not re.fullmatch(r"[0-9a-f]{64}", end["interval_sha256"]) or
            (end["dev"], end["ino"]) != (start["dev"], start["ino"]) or
            end["start_offset"] != start["offset"] or
            canonical_json(end["semantic_audit"]) != canonical_json(
                audit_telemetry_interval(
                    path, start, end,
                    end["semantic_audit"].get("audit_uptime_s"),
                )) or
            end["interval_sha256"] != end["semantic_audit"]["interval_sha256"]):
        die("external telemetry interval receipt mismatch")
    return end


def seal_campaign_completion(
    result_dir: Path, telemetry_end: Mapping[str, Any] | None,
) -> dict[str, Any]:
    """Cross-seal terminal results with successful telemetry disposition."""
    terminal_path = result_dir / "terminal.json"
    if not terminal_path.exists():
        die("cannot complete a campaign without terminal.json")
    null_path = result_dir / "telemetry_null.json"
    start_path = result_dir / "telemetry_start.json"
    end_path = result_dir / "telemetry_end.json"
    if telemetry_end is None:
        if start_path.exists() or end_path.exists():
            die("unexpected telemetry receipt for a null-telemetry campaign")
        write_sealed_once(
            null_path, f"{SCHEMA}.telemetry_null",
            {"external_telemetry_configured": False},
        )
        telemetry_name = null_path.name
        telemetry_path = null_path
    else:
        if (null_path.exists() or not start_path.exists() or
                not end_path.exists()):
            die("external telemetry completion receipt inventory mismatch")
        telemetry_name = end_path.name
        telemetry_path = end_path
    terminal = load_sealed(terminal_path, f"{SCHEMA}.terminal")
    stage_chain = campaign_stage_chain(result_dir, terminal)
    payload = {
        "contract_sha256": sha256_file(result_dir / "contract.json"),
        "controller_affinity_sha256": sha256_file(
            result_dir / "controller_affinity.json"),
        "terminal_sha256": sha256_file(terminal_path),
        "stage_entry_count": len(stage_chain),
        "stage_chain_sha256": sha256_bytes(
            canonical_json(stage_chain).encode()),
        "telemetry_receipt": telemetry_name,
        "telemetry_receipt_sha256": sha256_file(telemetry_path),
        "external_telemetry_bound": telemetry_end is not None,
    }
    write_sealed_once(
        result_dir / "campaign_complete.json",
        f"{SCHEMA}.campaign_complete", payload,
    )
    stored = load_sealed(
        result_dir / "campaign_complete.json",
        f"{SCHEMA}.campaign_complete",
    )
    if canonical_json(stored) != canonical_json(payload):
        die("campaign completion cross-seal mismatch")
    return payload


def verify_campaign_completion(
    result_dir: Path, contract: Mapping[str, Any],
    telemetry_start_receipt: Mapping[str, Any] | None,
    telemetry_end_receipt: Mapping[str, Any] | None,
) -> dict[str, Any]:
    complete_path = result_dir / "campaign_complete.json"
    if not complete_path.exists():
        die("campaign_complete.json is required before reduction")
    if ((contract["external_telemetry"]["path"] is None) !=
            (telemetry_start_receipt is None) or
            (telemetry_start_receipt is None) != (telemetry_end_receipt is None)):
        die("campaign completion telemetry disposition mismatch")
    null_path = result_dir / "telemetry_null.json"
    start_path = result_dir / "telemetry_start.json"
    end_path = result_dir / "telemetry_end.json"
    if telemetry_end_receipt is None:
        if (not null_path.exists() or start_path.exists() or
                end_path.exists()):
            die("null-telemetry completion receipt inventory mismatch")
        null = load_sealed(null_path, f"{SCHEMA}.telemetry_null")
        if canonical_json(null) != canonical_json({
                "external_telemetry_configured": False}):
            die("null-telemetry completion receipt mismatch")
        telemetry_path = null_path
    else:
        if (null_path.exists() or not start_path.exists() or
                not end_path.exists()):
            die("external telemetry completion receipt inventory mismatch")
        telemetry_path = end_path
    terminal = load_sealed(
        result_dir / "terminal.json", f"{SCHEMA}.terminal",
    )
    stage_chain = campaign_stage_chain(result_dir, terminal)
    expected = {
        "contract_sha256": sha256_file(result_dir / "contract.json"),
        "controller_affinity_sha256": sha256_file(
            result_dir / "controller_affinity.json"),
        "terminal_sha256": sha256_file(result_dir / "terminal.json"),
        "stage_entry_count": len(stage_chain),
        "stage_chain_sha256": sha256_bytes(
            canonical_json(stage_chain).encode()),
        "telemetry_receipt": telemetry_path.name,
        "telemetry_receipt_sha256": sha256_file(telemetry_path),
        "external_telemetry_bound": telemetry_end_receipt is not None,
    }
    stored = load_sealed(complete_path, f"{SCHEMA}.campaign_complete")
    if canonical_json(stored) != canonical_json(expected):
        die("campaign completion cross-seal mismatch")
    return stored


def bind_controller_affinity(result_dir: Path, cpus: Sequence[int]) -> None:
    selected = sorted(cpus)
    if not selected or 127 in selected:
        die("controller affinity is empty or includes reserved CPU 127")
    try:
        os.sched_setaffinity(0, set(selected))
    except OSError as exc:
        die(f"cannot bind controller affinity: {exc}")
    actual = sorted(os.sched_getaffinity(0))
    if actual != selected:
        die("controller affinity differs from the selected non-127 CPU set")
    write_sealed_once(
        result_dir / "controller_affinity.json",
        f"{SCHEMA}.controller_affinity",
        {"selected_cpus": selected, "actual_cpus": actual,
         "reserved_cpu_127_excluded": True},
    )


def apply_sealed_controller_affinity(result_dir: Path) -> list[int]:
    receipt = load_sealed(
        result_dir / "controller_affinity.json",
        f"{SCHEMA}.controller_affinity",
    )
    if (not isinstance(receipt, dict) or set(receipt) != {
            "selected_cpus", "actual_cpus", "reserved_cpu_127_excluded"} or
            type(receipt["selected_cpus"]) is not list or
            any(type(cpu) is not int for cpu in receipt["selected_cpus"]) or
            type(receipt["actual_cpus"]) is not list or
            any(type(cpu) is not int for cpu in receipt["actual_cpus"]) or
            receipt["actual_cpus"] != receipt["selected_cpus"] or
            receipt["reserved_cpu_127_excluded"] is not True or
            not receipt["selected_cpus"] or 127 in receipt["selected_cpus"] or
            receipt["selected_cpus"] != sorted(set(receipt["selected_cpus"]))):
        die("sealed controller affinity receipt is malformed")
    bind_controller_affinity(result_dir, receipt["selected_cpus"])
    return receipt["selected_cpus"]


def root_terminal_payload(
    construction_index: int, stage: Path, overhead: int, unresolved: int,
) -> dict[str, Any]:
    return {
        "campaign_round": campaign_round(construction_index),
        "construction_seed_index": construction_index,
        "construction_seed": str(CONSTRUCTION_SEEDS[construction_index]),
        "terminal_overhead": overhead,
        "unresolved_paired_cells": unresolved,
        "last_stage_complete_sha256": sha256_file(stage / "complete.json"),
    }


def failure_gate(
    name: str, construction_indices: Sequence[int],
    oh0_metrics_by_root: Mapping[int, Mapping[str, Any]],
) -> dict[str, Any]:
    expected = list(construction_indices)
    if (not expected or any(type(index) is not int for index in expected) or
            expected != sorted(set(expected)) or
            any(index not in range(len(CONSTRUCTION_SEEDS)) for index in expected) or
            set(expected) - set(oh0_metrics_by_root)):
        die("failure gate construction-root set is malformed")
    failures = {
        arm: sum(
            int(oh0_metrics_by_root[index]["arms"][arm]["rank_fail"]) +
            int(oh0_metrics_by_root[index]["arms"][arm]["error"])
            for index in expected
        ) for arm in ARMS
    }
    return {
        "name": name,
        "construction_seed_indices": expected,
        "paired_cells": PAIRED_CELLS_PER_ROOT * len(expected),
        "h12_failures": failures["h12"],
        "h13_failures": failures["h13"],
        "h12_minus_h13_failures": failures["h12"] - failures["h13"],
        "h13_strictly_fewer_failures": failures["h13"] < failures["h12"],
        "gate_uses_seed_fixes": False,
        "weak_seed_counts_are_descriptive": True,
    }


def terminal_payload(
    disposition: str, root_terminals: Sequence[Mapping[str, Any]],
    oh0_metrics_by_root: Mapping[int, Mapping[str, Any]],
) -> dict[str, Any]:
    if disposition not in ("r0_rejected", "confirmation_complete"):
        die("unknown terminal campaign disposition")
    expected_indices = [0] if disposition == "r0_rejected" else [0, 1, 2]
    if (len(root_terminals) != len(expected_indices) or
            any(not isinstance(item, Mapping) for item in root_terminals) or
            [item.get("construction_seed_index") for item in root_terminals] !=
                expected_indices):
        die("terminal root records are not in exact campaign order")
    r0 = failure_gate("R0", [0], oh0_metrics_by_root)
    if disposition == "r0_rejected":
        if r0["h13_strictly_fewer_failures"] or len(root_terminals) != 1:
            die("R0 rejection terminal contradicts the screening gate")
        r1 = None
        full = None
        accepted = False
    else:
        if not r0["h13_strictly_fewer_failures"] or len(root_terminals) != 3:
            die("confirmation terminal contradicts the R0 screening gate")
        r1 = failure_gate("R1", [1, 2], oh0_metrics_by_root)
        full = failure_gate("FULL", [0, 1, 2], oh0_metrics_by_root)
        accepted = bool(
            r1["h13_strictly_fewer_failures"] and
            full["h13_strictly_fewer_failures"]
        )
    return {
        "disposition": disposition,
        "architecture_accepted": accepted,
        "executed_construction_seed_indices": [
            item["construction_seed_index"] for item in root_terminals
        ],
        "root_terminals": list(root_terminals),
        "gates": {"R0": r0, "R1": r1, "FULL": full},
        "production_seed_fixups_applied": 0,
    }


def campaign_stage_chain(
    result_dir: Path, terminal: Any,
) -> list[dict[str, Any]]:
    if not isinstance(terminal, dict) or not isinstance(
            terminal.get("root_terminals"), list):
        die("campaign terminal lacks a root-terminal list")
    executed: dict[int, int] = {}
    for item in terminal["root_terminals"]:
        if (not isinstance(item, dict) or
                type(item.get("construction_seed_index")) is not int or
                type(item.get("terminal_overhead")) is not int):
            die("campaign root-terminal entry is malformed")
        construction_index = item["construction_seed_index"]
        terminal_overhead = item["terminal_overhead"]
        if (construction_index in executed or
                construction_index not in range(len(CONSTRUCTION_SEEDS)) or
                not 0 <= terminal_overhead <= MAX_OVERHEAD):
            die("campaign root-terminal domain is malformed")
        executed[construction_index] = terminal_overhead
    entries: list[dict[str, Any]] = []
    for construction_index in range(len(CONSTRUCTION_SEEDS)):
        last = executed.get(construction_index, -1)
        # Every OH0 plan is predeclared.  An unexecuted held-out root therefore
        # contributes its manifest, but deliberately has no completion seal.
        for overhead in range(max(last, 0) + 1):
            stage = stage_path(result_dir, construction_index, overhead)
            complete = stage / "complete.json"
            entries.append({
                "construction_seed_index": construction_index,
                "overhead": overhead,
                "manifest_sha256": sha256_file(stage / "manifest.json"),
                "complete_sha256": (
                    sha256_file(complete) if overhead <= last else None
                ),
            })
    return entries


def validate_campaign_root_inventory(
    result_dir: Path, contract: Mapping[str, Any],
) -> None:
    """Reject files and directories outside the finite campaign ledger."""
    allowed = {
        "contract.json", "controller.lock", "frozen", "roots",
        "controller_affinity.json", "telemetry_start.json",
        "telemetry_end.json", "telemetry_null.json", "terminal.json",
        "campaign_complete.json", "analysis.json",
    }
    # Tiny unit-test contracts keep their stand-in executables at the root.
    # A real sealed contract places all three tools under frozen/ instead.
    for name in ("binary", "controller", "taskset"):
        value = contract.get(name)
        if isinstance(value, str) and Path(value).parent == result_dir:
            allowed.add(Path(value).name)
    entries = {entry.name: entry for entry in result_dir.iterdir()}
    unexpected = set(entries) - allowed
    if unexpected:
        die("campaign root contains out-of-ledger entries")
    for name, entry in entries.items():
        if entry.is_symlink():
            die("campaign root contains a symlinked entry")
        if name in ("frozen", "roots"):
            if not entry.is_dir():
                die("campaign root directory entry has the wrong type")
        elif not entry.is_file():
            die("campaign root file entry has the wrong type")

    frozen = result_dir / "frozen"
    expected_frozen = {
        Path(contract[name]).name
        for name in ("binary", "controller", "taskset")
        if isinstance(contract.get(name), str) and
        Path(contract[name]).parent == frozen
    }
    if expected_frozen:
        if (not frozen.is_dir() or frozen.is_symlink() or
                {entry.name for entry in frozen.iterdir()} != expected_frozen or
                any(not entry.is_file() or entry.is_symlink()
                    for entry in frozen.iterdir())):
            die("frozen tool inventory differs from the sealed contract")


def validate_stage_root_inventory(
    result_dir: Path, construction_index: int,
    terminal_overhead: int | None,
) -> None:
    roots = result_dir / "roots"
    expected_roots = {f"c{index}" for index in range(len(CONSTRUCTION_SEEDS))}
    if (not roots.is_dir() or roots.is_symlink() or
            {entry.name for entry in roots.iterdir()} != expected_roots or
            any(not entry.is_dir() or entry.is_symlink()
                for entry in roots.iterdir())):
        die("campaign construction-root inventory mismatch")
    construction_root = construction_root_path(result_dir, construction_index)
    if ({entry.name for entry in construction_root.iterdir()} != {"stages"} or
            any(not entry.is_dir() or entry.is_symlink()
                for entry in construction_root.iterdir())):
        die("construction-root layout differs from its sealed plan")
    root = construction_root / "stages"
    entries = list(root.iterdir())
    last = 0 if terminal_overhead is None else terminal_overhead
    expected = {f"oh{overhead:04d}" for overhead in range(last + 1)}
    if ({entry.name for entry in entries} != expected or
            any(not entry.is_dir() or entry.is_symlink() for entry in entries)):
        die("campaign stage inventory differs from the terminal chain")
    if terminal_overhead is None:
        stage = stage_path(result_dir, construction_index, 0)
        if ({entry.name for entry in stage.iterdir()} != {"manifest.json"} or
                any(not entry.is_file() or entry.is_symlink()
                    for entry in stage.iterdir())):
            die("unexecuted predeclared OH0 plan contains outcome artifacts")


def validate_nested_arm_monotonicity(
    stage: Path, tasks: Sequence[Task], binary: Path, binary_sha256: str,
    taskset: Path, allowed_cpus: set[int],
    previous: Mapping[tuple[int, int, str, int], Mapping[str, bool]] | None,
) -> dict[tuple[int, int, str, int], dict[str, bool]]:
    """Return union failures while rejecting per-arm success reversals."""
    current: dict[tuple[int, int, str, int], dict[str, bool]] = {}
    seen: set[tuple[int, int, str, int]] = set()
    for key, pair in iter_stage_pairs(
            stage, tasks, binary, binary_sha256, taskset, allowed_cpus):
        failures = {arm: row_failed(pair[arm]) for arm in ARMS}
        if previous is not None:
            if key not in previous or key in seen:
                die("nested stage has an unexpected or duplicate arm state")
            for arm in ARMS:
                if not previous[key][arm] and failures[arm]:
                    die(f"nested packet stream reverted {arm} success at {key}")
            seen.add(key)
        if any(failures.values()):
            current[key] = failures
    if previous is not None and seen != set(previous):
        die("nested stage omitted a prior union-failure arm state")
    return current


def verify_terminal_campaign(
    result_dir: Path, contract: Mapping[str, Any], allowed_cpus: set[int],
) -> dict[str, Any]:
    validate_campaign_root_inventory(result_dir, contract)
    terminal = load_sealed(
        result_dir / "terminal.json", f"{SCHEMA}.terminal",
    )
    if (not isinstance(terminal, dict) or set(terminal) != {
            "disposition", "architecture_accepted",
            "executed_construction_seed_indices", "root_terminals", "gates",
            "production_seed_fixups_applied"} or
            terminal["disposition"] not in
                ("r0_rejected", "confirmation_complete") or
            type(terminal["architecture_accepted"]) is not bool or
            type(terminal["executed_construction_seed_indices"]) is not list or
            type(terminal["root_terminals"]) is not list or
            not isinstance(terminal["gates"], dict) or
            type(terminal["production_seed_fixups_applied"]) is not int or
            terminal["production_seed_fixups_applied"] != 0):
        die("terminal campaign receipt is malformed")
    binary = Path(contract["binary"])
    taskset = Path(contract["taskset"])
    root_records = terminal["root_terminals"]
    executed_indices = terminal["executed_construction_seed_indices"]
    if (any(not isinstance(item, dict) for item in root_records) or
            executed_indices != [
                item.get("construction_seed_index") for item in root_records
            ] or executed_indices not in ([0], [0, 1, 2])):
        die("terminal executed construction-root ordering mismatch")
    oh0_metrics: dict[int, Mapping[str, Any]] = {}
    recomputed_root_records: list[dict[str, Any]] = []
    for construction_index in range(len(CONSTRUCTION_SEEDS)):
        if construction_index not in executed_indices:
            validate_stage_root_inventory(result_dir, construction_index, None)
            tasks = load_manifest(
                stage_path(result_dir, construction_index, 0),
                construction_index, 0,
            )
            if len(tasks) != OH0_JOBS_PER_ROOT:
                die("unexecuted OH0 plan cardinality mismatch")
            continue
        record = root_records[executed_indices.index(construction_index)]
        terminal_overhead = record.get("terminal_overhead")
        if (type(terminal_overhead) is not int or
                not 0 <= terminal_overhead <= MAX_OVERHEAD):
            die("root terminal overhead is malformed")
        validate_stage_root_inventory(
            result_dir, construction_index, terminal_overhead,
        )
        previous_union: list[dict[str, Any]] | None = None
        previous_arm_failures: dict[
            tuple[int, int, str, int], dict[str, bool]
        ] | None = None
        for overhead in range(terminal_overhead + 1):
            if previous_union is not None and not previous_union:
                die("root chain continues after union failures resolved")
            stage = stage_path(result_dir, construction_index, overhead)
            tasks = load_manifest(stage, construction_index, overhead)
            if previous_union is not None:
                expected = {
                    (item["construction_seed_index"],
                     item["packet_trace_root_index"], item["schedule"], item["K"])
                    for item in previous_union
                }
                actual = {
                    (task.construction_index, task.seed_index, task.schedule, K)
                    for task in tasks[::2] for K in task.ks
                }
                if actual != expected:
                    die("terminal nested manifest is not the prior union cohort")
            completed = verify_complete_stage(
                stage, overhead, tasks, binary, contract["binary_sha256"],
                taskset, allowed_cpus,
            )
            if overhead == 0:
                oh0_metrics[construction_index] = completed["metrics"]
            current_arm_failures = validate_nested_arm_monotonicity(
                stage, tasks, binary, contract["binary_sha256"], taskset,
                allowed_cpus, previous_arm_failures,
            )
            sealed_union_keys = {
                (item["construction_seed_index"],
                 item["packet_trace_root_index"], item["schedule"], item["K"])
                for item in completed["union_failures"]
            }
            if set(current_arm_failures) != sealed_union_keys:
                die("terminal arm states differ from the sealed union cohort")
            previous_arm_failures = current_arm_failures
            previous_union = completed["union_failures"]
        final_stage = stage_path(
            result_dir, construction_index, terminal_overhead,
        )
        recomputed_root_records.append(root_terminal_payload(
            construction_index, final_stage, terminal_overhead,
            len(previous_union or []),
        ))
        screened_stop = (
            terminal["disposition"] == "r0_rejected" and
            construction_index == 0 and terminal_overhead == 0
        )
        if ((previous_union and terminal_overhead != MAX_OVERHEAD and
             not screened_stop) or
                stage_path(result_dir, construction_index,
                           terminal_overhead + 1).exists()):
            die("terminal construction-root chain is incomplete or overlong")
    recomputed = terminal_payload(
        terminal["disposition"], recomputed_root_records, oh0_metrics,
    )
    if canonical_json(terminal) != canonical_json(recomputed):
        die("terminal campaign receipt/chain mismatch")
    return terminal


def dispatch_tasks(
    result_dir: Path, stage: Path, tasks: Sequence[Task], binary: Path,
    binary_sha256: str, taskset: Path, taskset_sha256: str,
    cpu_queue: "queue.Queue[int]", workers: int, timeout: float,
    allowed_cpus: set[int],
) -> None:
    """Run one manifest with fail-fast, signal-safe child ownership."""
    children = ActiveChildren()
    previous_handlers: dict[int, Any] = {}
    received_signal: list[int | None] = [None]
    binary_fd = open_authenticated_executable(binary, binary_sha256)
    try:
        taskset_fd = open_authenticated_executable(taskset, taskset_sha256)
    except BaseException:
        os.close(binary_fd)
        raise
    executor = ThreadPoolExecutor(max_workers=workers)

    def handle_signal(signum: int, _frame: Any) -> None:
        # Python delivers this on the main thread.  Keep the handler
        # lock-free/reentrant: the dispatch loop notices it within 100 ms and
        # performs all cancellation, kill, and reap work in ordinary control.
        if received_signal[0] is None:
            received_signal[0] = signum

    for signum in (signal.SIGINT, signal.SIGTERM):
        previous_handlers[signum] = signal.signal(signum, handle_signal)

    futures = []
    failure: BaseException | None = None
    try:
        for task in tasks:
            if received_signal[0] is not None:
                failure = CampaignInterrupted(received_signal[0])
                break
            futures.append(executor.submit(
                execute_task, result_dir, stage, task, binary,
                binary_sha256, taskset, taskset_sha256, cpu_queue,
                timeout, children, allowed_cpus, binary_fd, taskset_fd,
            ))
        pending = set(futures)
        while failure is None and pending:
            if received_signal[0] is not None:
                failure = CampaignInterrupted(received_signal[0])
                break
            done, pending = wait(
                pending, timeout=0.1, return_when=FIRST_COMPLETED,
            )
            for future in done:
                try:
                    future.result()
                except BaseException as exc:
                    failure = exc
                    break
    except BaseException as exc:
        failure = exc
    finally:
        # Ignore reentrant signals before touching registry locks or waiting
        # for workers.  The first signal has already been recorded.
        for signum in previous_handlers:
            signal.signal(signum, signal.SIG_IGN)
        if received_signal[0] is not None:
            failure = CampaignInterrupted(received_signal[0])
        if failure is not None:
            children.stop(
                failure.signum
                if isinstance(failure, CampaignInterrupted) else None,
            )
            for future in futures:
                future.cancel()
        # Python 3.8 lacks Executor.shutdown(cancel_futures=True).  Pending
        # futures were explicitly cancelled above; running workers return only
        # after their registered process groups have been killed and reaped.
        try:
            executor.shutdown(wait=True)
        finally:
            try:
                for signum, previous in previous_handlers.items():
                    signal.signal(signum, previous)
            finally:
                os.close(taskset_fd)
                os.close(binary_fd)
    if children.active_count != 0:
        die("benchmark child registry was not empty after executor shutdown")
    if received_signal[0] is not None:
        raise CampaignInterrupted(received_signal[0])
    if failure is not None:
        raise failure


def cleanup_orphan_seal_temporaries(result_dir: Path) -> None:
    """Remove only recognizable create-only seal artifacts after a crash."""
    root_targets = {
        "contract.json", "controller_affinity.json", "telemetry_start.json",
        "telemetry_end.json", "telemetry_null.json", "terminal.json",
        "campaign_complete.json", "analysis.json",
    }

    def cleanup_directory(directory: Path, targets: set[str]) -> None:
        if not directory.is_dir() or directory.is_symlink():
            die("seal-temporary parent is missing, symlinked, or not a directory")
        for entry in list(directory.iterdir()):
            target_name = next((
                target for target in targets
                if entry.name.startswith(f".{target}.create-")
            ), None)
            if target_name is None:
                continue
            metadata = entry.lstat()
            if (entry.is_symlink() or not stat.S_ISREG(metadata.st_mode) or
                    metadata.st_nlink not in (1, 2)):
                die("orphan seal temporary has an unsafe type or link count")
            if metadata.st_nlink == 2:
                target = directory / target_name
                try:
                    target_metadata = target.lstat()
                except OSError as exc:
                    die(f"linked seal temporary lacks its destination: {exc}")
                if (target.is_symlink() or
                        not stat.S_ISREG(target_metadata.st_mode) or
                        (metadata.st_dev, metadata.st_ino) !=
                        (target_metadata.st_dev, target_metadata.st_ino)):
                    die("orphan seal temporary is linked outside its destination")
            entry.unlink()

    cleanup_directory(result_dir, root_targets)
    roots = result_dir / "roots"
    if roots.exists():
        if not roots.is_dir() or roots.is_symlink():
            die("campaign construction-root directory is unsafe")
        for construction_root in roots.iterdir():
            if (not re.fullmatch(r"c[0-9]+", construction_root.name) or
                    not construction_root.is_dir() or
                    construction_root.is_symlink()):
                continue
            stages = construction_root / "stages"
            if not stages.is_dir() or stages.is_symlink():
                die("campaign stage root is symlinked or not a directory")
            for stage in stages.iterdir():
                if (re.fullmatch(r"oh[0-9]{4}", stage.name) and
                        stage.is_dir() and not stage.is_symlink()):
                    cleanup_directory(stage, {"manifest.json", "complete.json"})


@exclusive_campaign
def run(args: argparse.Namespace, result_dir: Path) -> None:
    cleanup_orphan_seal_temporaries(result_dir)
    contract = load_contract(result_dir, require_frozen_controller=True)
    cpus_values = parse_cpu_list(args.cpus)
    workers = resolve_worker_count(args.workers, cpus_values)
    bind_controller_affinity(result_dir, cpus_values)
    terminal_path = result_dir / "terminal.json"
    if ((result_dir / "campaign_complete.json").exists() and
            not terminal_path.exists()):
        die("campaign completion exists without terminal.json")
    if (result_dir / "telemetry_end.json").exists() and not terminal_path.exists():
        die("telemetry end exists without an immutable terminal campaign seal")
    if (terminal_path.exists() and
            contract["external_telemetry"]["path"] is not None and
            not (result_dir / "telemetry_start.json").exists()):
        die("terminal campaign lacks its external telemetry start seal")
    telemetry_receipt_start = telemetry_start(result_dir, contract)
    if terminal_path.exists():
        terminal = verify_terminal_campaign(
            result_dir, contract, set(cpus_values),
        )
        telemetry_receipt_end = telemetry_finish(
            result_dir, contract, telemetry_receipt_start,
        )
        seal_campaign_completion(result_dir, telemetry_receipt_end)
        print(canonical_json({
            **terminal, "already_terminal": True,
            "external_telemetry_bound": telemetry_receipt_end is not None,
            "telemetry_continuity": contract["external_telemetry"]["continuity"],
            "campaign_complete_sha256": sha256_file(
                result_dir / "campaign_complete.json"),
        }))
        return
    cpu_queue: queue.Queue[int] = queue.Queue()
    for cpu in cpus_values[:workers]:
        cpu_queue.put(cpu)
    binary = Path(contract["binary"])
    taskset = Path(contract["taskset"])

    def execute_stage(construction_index: int, overhead: int) -> dict[str, Any]:
        stage = stage_path(result_dir, construction_index, overhead)
        tasks = load_manifest(stage, construction_index, overhead)
        shards = stage / "shards"
        if shards.exists():
            if shards.is_symlink() or not shards.is_dir():
                die("stage shard root has an unsafe type")
            task_ids = {task.task_id for task in tasks}
            task_stems = {task.stem for task in tasks}
            for path in shards.iterdir():
                if path.name.startswith(".staging-"):
                    match = re.fullmatch(
                        r"\.staging-([0-9a-f]{24})-([1-9][0-9]*)",
                        path.name,
                    )
                    if (match is None or match.group(1) not in task_ids or
                            path.is_symlink() or not path.is_dir()):
                        die("stage contains an unsafe shard staging artifact")
                    shutil.rmtree(path)
                elif (path.name not in task_stems or path.is_symlink() or
                      not path.is_dir()):
                    die("stage shard root contains an out-of-ledger entry")
        dispatch_tasks(
            result_dir, stage, tasks, binary, contract["binary_sha256"],
            taskset, contract["taskset_sha256"], cpu_queue, workers,
            args.timeout, set(cpus_values),
        )
        return complete_stage(
            stage, overhead, tasks, binary, contract["binary_sha256"], taskset,
            set(cpus_values),
        )

    def finish_root(
        construction_index: int, oh0_completed: Mapping[str, Any],
    ) -> dict[str, Any]:
        overhead = 0
        completed = oh0_completed
        while True:
            union = completed["union_failures"]
            stage = stage_path(result_dir, construction_index, overhead)
            if not union or overhead == MAX_OVERHEAD:
                return root_terminal_payload(
                    construction_index, stage, overhead, len(union),
                )
            next_overhead = overhead + 1
            cohort = [
                (item["packet_trace_root_index"], item["schedule"], item["K"])
                for item in union
            ]
            next_tasks = build_tasks(
                next_overhead, construction_index, cohort,
            )
            write_sealed_once(
                stage_path(result_dir, construction_index, next_overhead) /
                    "manifest.json",
                f"{SCHEMA}.stage_manifest",
                manifest_payload(construction_index, next_overhead, next_tasks),
            )
            overhead = next_overhead
            completed = execute_stage(construction_index, overhead)

    oh0_completed: dict[int, dict[str, Any]] = {}
    oh0_completed[0] = execute_stage(0, 0)
    r0_gate = failure_gate("R0", [0], {0: oh0_completed[0]["metrics"]})
    if not r0_gate["h13_strictly_fewer_failures"]:
        root_terminals = [root_terminal_payload(
            0, stage_path(result_dir, 0, 0), 0,
            len(oh0_completed[0]["union_failures"]),
        )]
        terminal = terminal_payload(
            "r0_rejected", root_terminals, {0: oh0_completed[0]["metrics"]},
        )
    else:
        root0_terminal = finish_root(0, oh0_completed[0])
        # Both held-out OH0 censuses are completed before either held-out root
        # can consume outcome-derived nested manifests.
        oh0_completed[1] = execute_stage(1, 0)
        oh0_completed[2] = execute_stage(2, 0)
        root1_terminal = finish_root(1, oh0_completed[1])
        root2_terminal = finish_root(2, oh0_completed[2])
        terminal = terminal_payload(
            "confirmation_complete",
            [root0_terminal, root1_terminal, root2_terminal],
            {index: value["metrics"] for index, value in oh0_completed.items()},
        )
    write_sealed_once(terminal_path, f"{SCHEMA}.terminal", terminal)
    terminal = verify_terminal_campaign(result_dir, contract, set(cpus_values))
    telemetry_receipt_end = telemetry_finish(
        result_dir, contract, telemetry_receipt_start,
    )
    seal_campaign_completion(result_dir, telemetry_receipt_end)
    print(canonical_json({
        **terminal,
        "external_telemetry_bound": telemetry_receipt_end is not None,
        "telemetry_continuity": contract["external_telemetry"]["continuity"],
        "campaign_complete_sha256": sha256_file(
            result_dir / "campaign_complete.json"),
        "reduce_command": [
            str(contract["controller"]), "reduce", "--result-dir",
            str(result_dir),
        ],
    }))
    return


def nearest_rank_p99(values: Sequence[int | None]) -> int | None:
    if not values:
        die("cannot compute p99 of an empty cohort")
    ordered = sorted(MAX_OVERHEAD + 1 if value is None else value for value in values)
    value = ordered[math.ceil(Decimal("0.99") * len(ordered)) - 1]
    return None if value > MAX_OVERHEAD else value


def overhead_report(values: Sequence[int | None]) -> dict[str, Any]:
    histogram: dict[str, int] = {}
    for value in values:
        key = "censored" if value is None else str(value)
        histogram[key] = histogram.get(key, 0) + 1
    observed_thresholds = sorted({
        0, *(value for value in values if value is not None)
    })
    return {
        "cells": len(values),
        "p99": nearest_rank_p99(values) if values else None,
        "censored": sum(value is None for value in values),
        "maximum_observed": max(
            (value for value in values if value is not None), default=None,
        ),
        "histogram": dict(sorted(
            histogram.items(),
            key=lambda item: (
                MAX_OVERHEAD + 1 if item[0] == "censored" else int(item[0])
            ),
        )),
        "survivors_after_observed_overhead": {
            str(overhead): sum(
                value is None or value > overhead for value in values
            ) for overhead in observed_thresholds
        },
        "right_censored_values_are_null": True,
    }


@exclusive_campaign
def reduce_campaign(args: argparse.Namespace, result_dir: Path) -> None:
    cleanup_orphan_seal_temporaries(result_dir)
    contract = load_contract(result_dir, require_frozen_controller=True)
    allowed_cpus = set(apply_sealed_controller_affinity(result_dir))
    if not (result_dir / "campaign_complete.json").exists():
        die("campaign_complete.json is required before reduction")
    if (contract["external_telemetry"]["path"] is not None and
            (not (result_dir / "telemetry_start.json").exists() or
             not (result_dir / "telemetry_end.json").exists())):
        die("campaign external telemetry interval is not terminally sealed")
    telemetry_receipt_start = telemetry_start(result_dir, contract)
    terminal = verify_terminal_campaign(result_dir, contract, allowed_cpus)
    telemetry_receipt_end = telemetry_finish(
        result_dir, contract, telemetry_receipt_start,
    )
    campaign_completion = verify_campaign_completion(
        result_dir, contract, telemetry_receipt_start, telemetry_receipt_end,
    )
    binary = Path(contract["binary"])
    taskset = Path(contract["taskset"])
    root_terminal_by_index = {
        item["construction_seed_index"]: item
        for item in terminal["root_terminals"]
    }

    def reduce_root(construction_index: int) -> dict[str, Any]:
        terminal_overhead = root_terminal_by_index[construction_index][
            "terminal_overhead"
        ]
        states: dict[tuple[int, int, str, int], dict[str, Any]] = {}
        arm_oh0_successes = {arm: 0 for arm in ARMS}
        arm_oh0_failure_minimums: dict[str, list[int | None]] = {
            arm: [] for arm in ARMS
        }
        union_oh0_minimums: dict[str, list[int | None]] = {
            arm: [] for arm in ARMS
        }
        common_failure_oh0_minimums: dict[str, list[int | None]] = {
            arm: [] for arm in ARMS
        }

        def finalize_state(state: Mapping[str, Any]) -> None:
            for arm in ARMS:
                minimum = state["minimum"][arm]
                union_oh0_minimums[arm].append(minimum)
                if state["oh0_failed"][arm]:
                    arm_oh0_failure_minimums[arm].append(minimum)
                if all(state["oh0_failed"].values()):
                    common_failure_oh0_minimums[arm].append(minimum)

        previous_union: list[dict[str, Any]] | None = None
        oh0_metrics: dict[str, Any] | None = None
        for overhead in range(terminal_overhead + 1):
            stage = stage_path(result_dir, construction_index, overhead)
            tasks = load_manifest(stage, construction_index, overhead)
            if previous_union is not None:
                expected = {
                    (item["construction_seed_index"],
                     item["packet_trace_root_index"], item["schedule"], item["K"])
                    for item in previous_union
                }
                actual = {
                    (task.construction_index, task.seed_index, task.schedule, K)
                    for task in tasks[::2] for K in task.ks
                }
                if actual != expected:
                    die("nested manifest is not the prior union-failure cohort")
            completed = load_sealed(
                stage / "complete.json", f"{SCHEMA}.stage_complete",
            )
            if overhead == 0:
                oh0_metrics = completed["metrics"]
                if completed["paired_cell_count"] != PAIRED_CELLS_PER_ROOT:
                    die("construction-root OH0 cardinality is not 575,991")
            next_states: dict[tuple[int, int, str, int], dict[str, Any]] = {}
            seen: set[tuple[int, int, str, int]] = set()
            paired_count = 0
            for key, pair in iter_stage_pairs(
                    stage, tasks, binary, contract["binary_sha256"], taskset,
                    allowed_cpus, progress=True):
                paired_count += 1
                if overhead == 0:
                    failed = {arm: row_failed(pair[arm]) for arm in ARMS}
                    for arm in ARMS:
                        arm_oh0_successes[arm] += not failed[arm]
                    if not any(failed.values()):
                        continue
                    state: dict[str, Any] = {
                        "minimum": {
                            arm: None if failed[arm] else 0 for arm in ARMS
                        },
                        "oh0_failed": failed,
                    }
                else:
                    if key not in states or key in seen:
                        die("nested stage has an unexpected or duplicate cell")
                    seen.add(key)
                    state = states[key]
                    for arm in ARMS:
                        failed = row_failed(pair[arm])
                        if state["minimum"][arm] is not None and failed:
                            die(f"nested packet stream reverted success at {key}")
                        if state["minimum"][arm] is None and not failed:
                            state["minimum"][arm] = overhead
                if any(row_failed(pair[arm]) for arm in ARMS):
                    next_states[key] = state
                else:
                    finalize_state(state)
            if paired_count != completed["paired_cell_count"]:
                die("stage streaming reduction cardinality mismatch")
            if overhead > 0 and seen != set(states):
                die("nested stage omitted a prior union-failure cell")
            states = next_states
            previous_union = completed["union_failures"]
            expected_union = {
                (item["construction_seed_index"],
                 item["packet_trace_root_index"], item["schedule"], item["K"])
                for item in previous_union
            }
            if set(states) != expected_union:
                die("computed union failures differ from sealed stage receipt")
        for state in states.values():
            finalize_state(state)
        minimums = {
            arm: [0] * arm_oh0_successes[arm] + arm_oh0_failure_minimums[arm]
            for arm in ARMS
        }
        if oh0_metrics is None:
            die("construction-root reduction lacks its OH0 metrics")
        if (any(len(minimums[arm]) != PAIRED_CELLS_PER_ROOT for arm in ARMS) or
                any(len(union_oh0_minimums[arm]) !=
                    oh0_metrics["comparison"]["union_failures"] for arm in ARMS) or
                any(len(common_failure_oh0_minimums[arm]) !=
                    oh0_metrics["comparison"]["both_fail"] for arm in ARMS)):
            die("minimum-overhead cohort accounting mismatch")
        return {
            "oh0": oh0_metrics,
            "terminal_overhead": terminal_overhead,
            "unresolved_union_failures": len(previous_union or []),
            "minimum_values": {
                "all_cells": minimums,
                "among_arm_OH0_failures": arm_oh0_failure_minimums,
                "among_OH0_union_failures": union_oh0_minimums,
                "among_common_OH0_failures": common_failure_oh0_minimums,
            },
        }

    root_results = {
        construction_index: reduce_root(construction_index)
        for construction_index in terminal["executed_construction_seed_indices"]
    }

    def combined_oh0(indices: Sequence[int]) -> dict[str, Any]:
        def pairs() -> Iterable[
            tuple[tuple[int, int, str, int], dict[str, dict[str, str]]]
        ]:
            for construction_index in indices:
                stage = stage_path(result_dir, construction_index, 0)
                tasks = load_manifest(stage, construction_index, 0)
                yield from iter_stage_pairs(
                    stage, tasks, binary, contract["binary_sha256"], taskset,
                    allowed_cpus, progress=True,
                )
        return aggregate_pair_stream(pairs())[0]

    def tail_summary(indices: Sequence[int]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for arm in ARMS:
            result[arm] = {}
            for cohort in (
                    "all_cells", "among_arm_OH0_failures",
                    "among_OH0_union_failures", "among_common_OH0_failures"):
                values: list[int | None] = []
                for construction_index in indices:
                    values.extend(root_results[construction_index][
                        "minimum_values"][cohort][arm])
                result[arm][cohort] = overhead_report(values)
        return result

    def censoring_summary(indices: Sequence[int]) -> list[dict[str, Any]]:
        records = []
        for construction_index in indices:
            unresolved = root_results[construction_index][
                "unresolved_union_failures"
            ]
            if unresolved:
                rejected = terminal["disposition"] == "r0_rejected"
                records.append({
                    "construction_seed_index": construction_index,
                    "right_censoring_overhead": root_results[
                        construction_index]["terminal_overhead"],
                    "right_censored_union_cells": unresolved,
                    "reason": (
                        "R0 strict-better screen rejected; no OH1+ outcomes "
                        "were collected" if rejected else
                        "success was not observed through the OH1024 cap"
                    ),
                })
        return records

    r0_oh0 = combined_oh0([0])
    if canonical_json(r0_oh0) != canonical_json(root_results[0]["oh0"]):
        die("R0 pooled OH0 metrics differ from the construction-root receipt")
    rounds: dict[str, Any] = {
        "R0": {
            "construction_seed_indices": [0],
            "paired_cells": R0_PAIRED_CELLS,
            "oh0": r0_oh0,
            "minimum_success_overhead": tail_summary([0]),
            "right_censoring": censoring_summary([0]),
        },
        "R1": None,
        "FULL": None,
    }
    if terminal["disposition"] == "confirmation_complete":
        rounds["R1"] = {
            "construction_seed_indices": [1, 2],
            "paired_cells": R1_PAIRED_CELLS,
            "oh0": combined_oh0([1, 2]),
            "minimum_success_overhead": tail_summary([1, 2]),
            "right_censoring": censoring_summary([1, 2]),
        }
        rounds["FULL"] = {
            "construction_seed_indices": [0, 1, 2],
            "paired_cells": FULL_PAIRED_CELLS,
            "oh0": combined_oh0([0, 1, 2]),
            "minimum_success_overhead": tail_summary([0, 1, 2]),
            "right_censoring": censoring_summary([0, 1, 2]),
        }
    summary = {
        "schema": f"{SCHEMA}.analysis",
        "coverage": {
            "K_min": K_MIN, "K_max": K_MAX, "unique_K": K_COUNT,
            "construction_seed_policy": CONSTRUCTION_SEED_POLICY,
            "construction_seeds_predeclared": [
                str(seed) for seed in CONSTRUCTION_SEEDS
            ],
            "construction_seeds_executed": terminal[
                "executed_construction_seed_indices"],
            "packet_trace_roots": len(LOSS_ROOTS),
            "schedules": len(SCHEDULES),
            "paired_cells_per_construction_root": PAIRED_CELLS_PER_ROOT,
            "r0_paired_cells": R0_PAIRED_CELLS,
            "r1_paired_cells": (
                R1_PAIRED_CELLS
                if terminal["disposition"] == "confirmation_complete" else 0
            ),
            "full_paired_cells": (
                FULL_PAIRED_CELLS
                if terminal["disposition"] == "confirmation_complete"
                else R0_PAIRED_CELLS
            ),
            "production_seed_fixups_applied": 0,
        },
        "terminal": terminal,
        "rounds": rounds,
        "by_construction_seed": {
            str(index): {
                "construction_seed": str(CONSTRUCTION_SEEDS[index]),
                "campaign_round": campaign_round(index),
                "terminal_overhead": result["terminal_overhead"],
                "unresolved_union_failures": result[
                    "unresolved_union_failures"],
                "oh0": result["oh0"],
                "minimum_success_overhead": tail_summary([index]),
                "right_censoring": censoring_summary([index]),
            } for index, result in root_results.items()
        },
        "interpretation": (
            "Deterministic uniform-construction-root/loss-root/schedule census; "
            "raw weak (K,C), weak-K, McNemar, work, and K-cluster values are "
            "descriptive. Architecture gates use only strictly fewer H13 OH0 "
            "failures in R0, held-out pooled R1, and full totals; no seed fixes "
            "or introductions gate was applied. A censored minimum in a full "
            "confirmation means no success was observed through OH1024. In an "
            "R0-rejected campaign it means the cell was deliberately not "
            "followed beyond OH0, never an observed OH1025."
        ),
        "loss_seed_formula": LOSS_SEED_FORMULA,
        "external_telemetry": {
            "managed_by_controller": False,
            "policy": contract["external_telemetry"]["policy"],
            "sampler_identity": contract["external_telemetry"]
                ["sampler_identity"],
            "continuity": contract["external_telemetry"]["continuity"],
            "start": telemetry_receipt_start,
            "end": telemetry_receipt_end,
        },
        "campaign_completion": campaign_completion,
    }
    write_sealed_once(
        result_dir / "analysis.json", f"{SCHEMA}.analysis_record", summary,
    )
    print(canonical_json(summary))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument("--binary", type=Path, required=True)
    prepare_parser.add_argument("--result-dir", type=Path, required=True)
    prepare_parser.add_argument(
        "--telemetry-log", type=Path,
        help="existing externally managed append-only thermal/host telemetry log",
    )
    prepare_parser.set_defaults(handler=prepare)
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--result-dir", type=Path, required=True)
    run_parser.add_argument("--cpus", required=True)
    run_parser.add_argument("--workers", type=int)
    run_parser.add_argument("--timeout", type=float, default=300.0)
    run_parser.set_defaults(handler=run)
    reduce_parser = subparsers.add_parser("reduce")
    reduce_parser.add_argument("--result-dir", type=Path, required=True)
    reduce_parser.set_defaults(handler=reduce_campaign)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    try:
        args = build_parser().parse_args(argv)
        timeout = getattr(args, "timeout", 1.0)
        if not math.isfinite(timeout) or timeout <= 0:
            die("timeout must be finite and positive")
        args.handler(args)
        return 0
    except CampaignInterrupted as exc:
        print(f"wh2_h13_stage_a: {exc}", file=sys.stderr)
        return 128 + exc.signum
    except CampaignError as exc:
        print(f"wh2_h13_stage_a: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    if len(sys.argv) >= 2 and sys.argv[1] == "__pdeath_exec":
        try:
            if len(sys.argv) < 4 or not re.fullmatch(r"[1-9][0-9]*", sys.argv[2]):
                die("parent-death launcher has a malformed parent PID")
            raise SystemExit(pdeath_exec(int(sys.argv[2]), sys.argv[3:]))
        except CampaignError as exc:
            print(f"wh2_h13_stage_a: {exc}", file=sys.stderr)
            raise SystemExit(125)
    raise SystemExit(main())
