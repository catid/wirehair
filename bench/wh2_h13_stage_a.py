#!/usr/bin/env python3
"""Raw attempt-0 H12 versus H13 Stage-A campaign controller.

The controller deliberately has a small public surface:

    prepare  freeze this script, the benchmark, and the 4,608-job OH0 plan
    run      execute/restart sealed jobs and escalate paired failures to OH1024
    reduce   verify every seal again and emit the campaign summary

This is a deterministic census over named seeds and schedules.  McNemar values
and weak-K summaries are descriptive diagnostics, not IID population claims.
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


SCHEMA = "wirehair.wh2.h13_stage_a.v1"
K_MIN = 2
K_MAX = 64000
K_COUNT = K_MAX - K_MIN + 1
CHUNK_SIZE = 250
MAX_OVERHEAD = 1024
ARMS = ("h12", "h13")
ARM_GF256_ROWS = {"h12": 10, "h13": 11}
SEEDS = (
    "0xd1b54a32d192ed03",
    "0x94d049bb133111eb",
    "0x8538ecb5bd456ea3",
)
SCHEDULES = ("burst", "adversarial", "repair-only")
PAIRED_CELLS = K_COUNT * len(SEEDS) * len(SCHEDULES)
TOTAL_ARM_OUTCOMES = PAIRED_CELLS * len(ARMS)
OH0_PAIR_JOBS = math.ceil(K_COUNT / CHUNK_SIZE) * len(SEEDS) * len(SCHEDULES)
OH0_JOBS = OH0_PAIR_JOBS * len(ARMS)
MASK64 = (1 << 64) - 1
MAX_I63 = (1 << 63) - 1
OUTPUT_LIMIT = 8 << 20
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
    "construction_attempt", "base_matrix_seed", "base_peel_seed",
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
        if value < 0 or result > MAX_I63 // max(value, 1):
            die("campaign cardinality exceeds signed 63-bit accounting")
        result *= value
    return result


def strict_uint(text: str, context: str, maximum: int | None = None) -> int:
    if not re.fullmatch(r"0|[1-9][0-9]*", text):
        die(f"{context}: noncanonical unsigned integer {text!r}")
    value = int(text)
    if maximum is not None and value > maximum:
        die(f"{context}: integer exceeds {maximum}")
    return value


def strict_hex(text: str, context: str, bits: int) -> int:
    if not re.fullmatch(r"0x(?:0|[1-9a-f][0-9a-f]*)", text):
        die(f"{context}: noncanonical hexadecimal integer {text!r}")
    value = int(text, 16)
    if value >= 1 << bits:
        die(f"{context}: integer does not fit u{bits}")
    return value


def strict_fixed(text: str, context: str, places: int) -> Decimal:
    if not re.fullmatch(rf"(?:0|[1-9][0-9]*)\.[0-9]{{{places}}}", text):
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
        "seed": task.seed, "completion": "mixed", "mixed_period": "48",
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
        "seed_block_bytes_override": "64", "overhead_stream": "paired",
        "full_payload_solve": "0", "schedule": task.schedule,
        "raw_attempt0": "1",
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


@dataclass(frozen=True)
class Task:
    overhead: int
    job: int
    pair: int
    arm: str
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
            "arm": self.arm, "seed_index": self.seed_index,
            "seed": self.seed, "schedule": self.schedule,
            "ks": list(self.ks),
        }

    def payload(self) -> dict[str, Any]:
        return {**self.core(), "task_id": self.task_id}


def task_from_payload(payload: Any) -> Task:
    required = {
        "overhead", "job", "pair", "arm", "seed_index", "seed",
        "schedule", "ks", "task_id",
    }
    if not isinstance(payload, dict) or set(payload) != required:
        die("task has an unexpected schema")
    if (any(type(payload[key]) is not int for key in (
            "overhead", "job", "pair", "seed_index")) or
            any(type(payload[key]) is not str for key in (
                "arm", "seed", "schedule", "task_id")) or
            type(payload["ks"]) is not list or
            any(type(value) is not int for value in payload["ks"])):
        die("task scalar types are not canonical")
    task = Task(
        overhead=payload["overhead"], job=payload["job"],
        pair=payload["pair"], arm=payload["arm"],
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
            task.overhead, task.job, task.pair, task.seed_index, *task.ks)) or
            any(type(value) is not str for value in (
                task.arm, task.seed, task.schedule, task.task_id)) or
            task.overhead < 0 or task.overhead > MAX_OVERHEAD or task.job < 0 or
            task.pair < 0 or task.arm not in ARMS or
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
    overhead: int, job: int, pair: int, arm: str, seed_index: int,
    schedule: str, ks: Sequence[int],
) -> Task:
    provisional = Task(
        overhead, job, pair, arm, seed_index, SEEDS[seed_index], schedule,
        tuple(ks), "",
    )
    return Task(**{**provisional.__dict__, "task_id": task_digest(provisional.core())})


def build_tasks(
    overhead: int,
    cohort: Iterable[tuple[int, str, int]] | None = None,
) -> list[Task]:
    if overhead < 0 or overhead > MAX_OVERHEAD:
        die("overhead outside 0..1024")
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
                        overhead, len(tasks), pair, arm, seed_index,
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
        "--seed-block-bytes", "64", "--packet-peel-seed-xor", "0",
        "--paired-overhead-stream", "--seed", task.seed,
        "--schedule", task.schedule, "--raw-attempt0",
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
        if (strict_hex(row["matrix_seed"], "matrix", 64) != base_matrix or
                strict_hex(row["peel_seed"], "peel", 32) != base_peel):
            die(f"attempt-0 active seeds differ from base seeds at K={K}")
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
) -> dict[tuple[int, str, int], dict[str, dict[str, str]]]:
    cells: dict[tuple[int, str, int], dict[str, dict[str, str]]] = {}
    for task, rows in tasks_and_rows:
        if len(rows) != len(task.ks):
            die("parsed task row count mismatch")
        for K, row in zip(task.ks, rows):
            key = (task.seed_index, task.schedule, K)
            arms = cells.setdefault(key, {})
            if task.arm in arms:
                die(f"duplicate arm result for cell {key}")
            arms[task.arm] = row
    for key, arms in cells.items():
        validate_pair_receipt(key, arms)
    return cells


def validate_pair_receipt(
    key: tuple[int, str, int], arms: Mapping[str, Mapping[str, str]],
) -> None:
    if set(arms) != set(ARMS):
        die(f"unpaired Stage-A cell {key}")
    h12, h13 = arms["h12"], arms["h13"]
    for field in (
        "base_matrix_seed", "base_peel_seed", "matrix_seed", "peel_seed",
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
        tuple[tuple[int, str, int], Mapping[str, Mapping[str, str]]]
    ],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    overall = _empty_scope()
    by_seed = {str(i): _empty_scope() for i in range(len(SEEDS))}
    by_schedule = {name: _empty_scope() for name in SCHEDULES}
    by_seed_schedule = {
        f"seed{seed_index}:{schedule}": _empty_scope()
        for seed_index in range(len(SEEDS)) for schedule in SCHEDULES
    }
    failure_by_k: dict[str, dict[int, int]] = {arm: {} for arm in ARMS}
    failure_records: dict[str, dict[int, list[dict[str, Any]]]] = {
        arm: {} for arm in ARMS
    }
    sampled_ks: set[int] = set()
    union_failures: list[dict[str, Any]] = []
    previous_ordinal = -1
    for (seed_index, schedule, K), pair in pairs:
        ordinal = ((seed_index * len(SCHEDULES) + SCHEDULES.index(schedule)) *
                   K_COUNT + K - K_MIN)
        if ordinal <= previous_ordinal:
            die("stage cells are duplicated or not in canonical order")
        previous_ordinal = ordinal
        sampled_ks.add(K)
        _update_scope(overall, pair)
        _update_scope(by_seed[str(seed_index)], pair)
        _update_scope(by_schedule[schedule], pair)
        _update_scope(by_seed_schedule[f"seed{seed_index}:{schedule}"], pair)
        any_failure = False
        for arm in ARMS:
            row = pair[arm]
            failed = row_failed(row)
            any_failure |= failed
            failure_by_k[arm][K] = failure_by_k[arm].get(K, 0) + failed
            if failed:
                failure_records[arm].setdefault(K, []).append({
                    "packet_trace_root_index": seed_index,
                    "packet_trace_root": SEEDS[seed_index],
                    "schedule": schedule,
                    "cause": "rank_fail",
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
                "seed_index": seed_index, "seed": SEEDS[seed_index],
                "schedule": schedule, "K": K,
            })

    arms = overall["arms"]
    for arm in ARMS:
        histogram = {str(count): 0 for count in range(10)}
        for K in sampled_ks:
            count = failure_by_k[arm].get(K, 0)
            if count > 9:
                die("one K has more than nine named-stratum failures")
            histogram[str(count)] += 1
        arms[arm]["sampled_unique_K"] = len(sampled_ks)
        arms[arm]["weak_K"] = len(sampled_ks) - histogram["0"]
        arms[arm]["multi_failure_K"] = sum(
            histogram[str(index)] for index in range(2, 10)
        )
        arms[arm]["maximum_failure_strata"] = max(
            failure_by_k[arm].values(), default=0,
        )
        arms[arm]["failure_strata_histogram_0_to_9"] = histogram
        arms[arm]["weak_K_records"] = [
            {"K": K, "failure_strata": len(records), "failures": records}
            for K, records in sorted(failure_records[arm].items())
        ]
    report = _finish_scope(overall)
    report["descriptive_scopes"] = {
        "by_packet_trace_root": {
            key: _finish_scope(value) for key, value in by_seed.items()
        },
        "by_schedule": {
            key: _finish_scope(value) for key, value in by_schedule.items()
        },
        "by_packet_trace_root_schedule": {
            key: _finish_scope(value) for key, value in by_seed_schedule.items()
        },
    }
    report["comparison"]["union_failures"] = len(union_failures)
    return report, union_failures


def aggregate_cells(
    cells: Mapping[tuple[int, str, int], Mapping[str, Mapping[str, str]]],
) -> dict[str, Any]:
    return aggregate_pair_stream(cells.items())[0]


def manifest_payload(overhead: int, tasks: Sequence[Task]) -> dict[str, Any]:
    return {
        "overhead": overhead, "chunk_size": CHUNK_SIZE,
        "paired_overhead_stream": True, "tasks": [task.payload() for task in tasks],
    }


def load_manifest(stage: Path, overhead: int) -> list[Task]:
    payload = load_sealed(
        stage / "manifest.json", f"{SCHEMA}.stage_manifest",
    )
    if (not isinstance(payload, dict) or set(payload) != {
            "overhead", "chunk_size", "paired_overhead_stream", "tasks"} or
            type(payload["overhead"]) is not int or
            type(payload["chunk_size"]) is not int or
            payload["overhead"] != overhead or payload["chunk_size"] != CHUNK_SIZE or
            payload["paired_overhead_stream"] is not True or
            not isinstance(payload["tasks"], list)):
        die("stage manifest contract mismatch")
    tasks = [task_from_payload(value) for value in payload["tasks"]]
    if [task.job for task in tasks] != list(range(len(tasks))):
        die("stage task IDs are not contiguous unbounded integers")
    if (overhead == 0 and
            (len(tasks) != OH0_JOBS or
             sum(len(task.ks) for task in tasks) != TOTAL_ARM_OUTCOMES)):
        die("OH0 manifest does not contain the exact full campaign census")
    return tasks


def stage_path(result_dir: Path, overhead: int) -> Path:
    return result_dir / "stages" / f"oh{overhead:04d}"


def shard_path(stage: Path, task: Task) -> Path:
    return stage / "shards" / task.stem


def receipt_payload(
    task: Task, benchmark_argv: Sequence[str], command: Sequence[str], cpu: int,
    binary_sha256: str, stdout: bytes, stderr: bytes, row_count: int,
) -> dict[str, Any]:
    return {
        "task": task.payload(), "benchmark_argv": list(benchmark_argv),
        "command": list(command), "cpu": cpu, "binary_sha256": binary_sha256,
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
    inventory = {path.name for path in root.iterdir()}
    if inventory != {"stdout.csv", "stderr.txt", "receipt.json"}:
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


def _run_process(
    command: Sequence[str], timeout: float, children: ActiveChildren,
) -> tuple[int, bytes, bytes]:
    children.check()
    launcher = [
        str(Path(__file__).resolve()), "__pdeath_exec", str(os.getpid()),
        *command,
    ]
    process = subprocess.Popen(
        launcher, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        start_new_session=True,
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
    allowed_cpus: set[int],
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
        if sha256_file(binary) != binary_sha256 or sha256_file(taskset) != taskset_sha256:
            die("frozen executable changed before job launch")
        children.check()
        benchmark_argv = make_benchmark_argv(binary, task)
        command = [str(taskset), "-c", str(cpu), *benchmark_argv]
        returncode, stdout, stderr = _run_process(command, timeout, children)
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
    tuple[tuple[int, str, int], dict[str, dict[str, str]]]
]:
    if len(tasks) % 2:
        die("stage manifest has an odd task count")
    pair_count = len(tasks) // 2
    for pair_index in range(pair_count):
        left, right = tasks[2 * pair_index:2 * pair_index + 2]
        if (left.arm != "h12" or right.arm != "h13" or
                left.pair != pair_index or right.pair != pair_index or
                left.overhead != right.overhead or
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
            key = (left.seed_index, left.schedule, K)
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
    if (checked_product(K_COUNT, len(SEEDS), len(SCHEDULES)) != PAIRED_CELLS or
            checked_product(PAIRED_CELLS, len(ARMS)) != TOTAL_ARM_OUTCOMES):
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
            "packet_trace_roots": list(SEEDS),
            "seed_role": (
                "external packet/loss-trace root; not a construction seed"
            ),
            "schedules": list(SCHEDULES),
            "paired_cells": PAIRED_CELLS, "outcomes_per_arm": PAIRED_CELLS,
            "arm_outcomes_at_OH0": TOTAL_ARM_OUTCOMES,
            "oh0_jobs": OH0_JOBS,
            "child_lifetime_policy": (
                "dedicated process group; registered kill/reap on SIGINT or "
                "SIGTERM; Linux PR_SET_PDEATHSIG=SIGKILL"
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
        tasks = build_tasks(0)
        initial_stage = staging / "stages" / "oh0000"
        write_sealed_once(
            initial_stage / "manifest.json", f"{SCHEMA}.stage_manifest",
            manifest_payload(0, tasks),
        )
        staging.rename(result_dir)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    record = {
        "result_dir": str(result_dir), "oh0_jobs": OH0_JOBS,
        "paired_cells": PAIRED_CELLS, "outcomes_per_arm": PAIRED_CELLS,
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
        "packet_trace_roots": list(SEEDS),
        "seed_role": "external packet/loss-trace root; not a construction seed",
        "schedules": list(SCHEDULES),
        "paired_cells": PAIRED_CELLS, "outcomes_per_arm": PAIRED_CELLS,
        "arm_outcomes_at_OH0": TOTAL_ARM_OUTCOMES,
        "oh0_jobs": OH0_JOBS, "loss_seed_formula": LOSS_SEED_FORMULA,
        "child_lifetime_policy": (
            "dedicated process group; registered kill/reap on SIGINT or "
            "SIGTERM; Linux PR_SET_PDEATHSIG=SIGKILL"
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
        if path.parent != result_dir / "frozen" or \
                sha256_file(path) != contract[f"{name}_sha256"]:
            die(f"frozen {name} substitution detected")
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
    terminal_overhead = terminal.get("terminal_overhead") \
        if isinstance(terminal, dict) else None
    if type(terminal_overhead) is not int or not 0 <= terminal_overhead <= MAX_OVERHEAD:
        die("cannot build campaign completion from malformed terminal receipt")
    stage_chain = [{
        "overhead": overhead,
        "manifest_sha256": sha256_file(
            stage_path(result_dir, overhead) / "manifest.json"),
        "complete_sha256": sha256_file(
            stage_path(result_dir, overhead) / "complete.json"),
    } for overhead in range(terminal_overhead + 1)]
    payload = {
        "contract_sha256": sha256_file(result_dir / "contract.json"),
        "controller_affinity_sha256": sha256_file(
            result_dir / "controller_affinity.json"),
        "terminal_sha256": sha256_file(terminal_path),
        "stage_count": len(stage_chain),
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
    terminal_overhead = terminal.get("terminal_overhead") \
        if isinstance(terminal, dict) else None
    if type(terminal_overhead) is not int or not 0 <= terminal_overhead <= MAX_OVERHEAD:
        die("campaign completion terminal receipt is malformed")
    stage_chain = [{
        "overhead": overhead,
        "manifest_sha256": sha256_file(
            stage_path(result_dir, overhead) / "manifest.json"),
        "complete_sha256": sha256_file(
            stage_path(result_dir, overhead) / "complete.json"),
    } for overhead in range(terminal_overhead + 1)]
    expected = {
        "contract_sha256": sha256_file(result_dir / "contract.json"),
        "controller_affinity_sha256": sha256_file(
            result_dir / "controller_affinity.json"),
        "terminal_sha256": sha256_file(result_dir / "terminal.json"),
        "stage_count": len(stage_chain),
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


def terminal_payload(
    stage: Path, overhead: int, unresolved: int,
) -> dict[str, Any]:
    return {
        "terminal_overhead": overhead,
        "unresolved_paired_cells": unresolved,
        "last_stage_complete_sha256": sha256_file(stage / "complete.json"),
    }


def validate_campaign_root_inventory(
    result_dir: Path, contract: Mapping[str, Any],
) -> None:
    """Reject files and directories outside the finite campaign ledger."""
    allowed = {
        "contract.json", "controller.lock", "frozen", "stages",
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
        if name in ("frozen", "stages"):
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


def validate_stage_root_inventory(result_dir: Path, terminal_overhead: int) -> None:
    root = result_dir / "stages"
    if not root.is_dir() or root.is_symlink():
        die("campaign stage root is missing, symlinked, or not a directory")
    entries = list(root.iterdir())
    expected = {f"oh{overhead:04d}" for overhead in range(terminal_overhead + 1)}
    if ({entry.name for entry in entries} != expected or
            any(not entry.is_dir() or entry.is_symlink() for entry in entries)):
        die("campaign stage inventory differs from the terminal chain")


def validate_nested_arm_monotonicity(
    stage: Path, tasks: Sequence[Task], binary: Path, binary_sha256: str,
    taskset: Path, allowed_cpus: set[int],
    previous: Mapping[tuple[int, str, int], Mapping[str, bool]] | None,
) -> dict[tuple[int, str, int], dict[str, bool]]:
    """Return union failures while rejecting per-arm success reversals."""
    current: dict[tuple[int, str, int], dict[str, bool]] = {}
    seen: set[tuple[int, str, int]] = set()
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
            "terminal_overhead", "unresolved_paired_cells",
            "last_stage_complete_sha256"} or
            type(terminal["terminal_overhead"]) is not int or
            type(terminal["unresolved_paired_cells"]) is not int or
            not 0 <= terminal["terminal_overhead"] <= MAX_OVERHEAD or
            terminal["unresolved_paired_cells"] < 0 or
            not isinstance(terminal["last_stage_complete_sha256"], str) or
            not re.fullmatch(r"[0-9a-f]{64}",
                             terminal["last_stage_complete_sha256"])):
        die("terminal campaign receipt is malformed")
    binary = Path(contract["binary"])
    taskset = Path(contract["taskset"])
    validate_stage_root_inventory(result_dir, terminal["terminal_overhead"])
    previous_union: list[dict[str, Any]] | None = None
    previous_arm_failures: dict[
        tuple[int, str, int], dict[str, bool]
    ] | None = None
    for overhead in range(terminal["terminal_overhead"] + 1):
        if previous_union is not None and not previous_union:
            die("terminal chain continues after the union failure cohort resolved")
        stage = stage_path(result_dir, overhead)
        tasks = load_manifest(stage, overhead)
        if previous_union is not None:
            expected = {
                (item["seed_index"], item["schedule"], item["K"])
                for item in previous_union
            }
            actual = {
                (task.seed_index, task.schedule, K)
                for task in tasks[::2] for K in task.ks
            }
            if actual != expected:
                die("terminal nested manifest is not the prior union cohort")
        completed = verify_complete_stage(
            stage, overhead, tasks, binary, contract["binary_sha256"], taskset,
            allowed_cpus,
        )
        current_arm_failures = validate_nested_arm_monotonicity(
            stage, tasks, binary, contract["binary_sha256"], taskset,
            allowed_cpus, previous_arm_failures,
        )
        sealed_union_keys = {
            (item["seed_index"], item["schedule"], item["K"])
            for item in completed["union_failures"]
        }
        if set(current_arm_failures) != sealed_union_keys:
            die("terminal arm states differ from the sealed union cohort")
        previous_arm_failures = current_arm_failures
        previous_union = completed["union_failures"]
    final_stage = stage_path(result_dir, terminal["terminal_overhead"])
    if (canonical_json(terminal) != canonical_json(terminal_payload(
            final_stage, terminal["terminal_overhead"],
            len(previous_union or []))) or
            (previous_union and terminal["terminal_overhead"] != MAX_OVERHEAD) or
            stage_path(result_dir, terminal["terminal_overhead"] + 1).exists()):
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
                timeout, children, allowed_cpus,
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
            for signum, previous in previous_handlers.items():
                signal.signal(signum, previous)
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
    stages = result_dir / "stages"
    if stages.exists():
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
    overhead = 0
    while True:
        stage = stage_path(result_dir, overhead)
        tasks = load_manifest(stage, overhead)
        shards = stage / "shards"
        if shards.exists():
            for path in shards.iterdir():
                if path.name.startswith(".staging-"):
                    shutil.rmtree(path)
        dispatch_tasks(
            result_dir, stage, tasks, binary, contract["binary_sha256"],
            taskset, contract["taskset_sha256"], cpu_queue, workers,
            args.timeout, set(cpus_values),
        )
        completed = complete_stage(
            stage, overhead, tasks, binary, contract["binary_sha256"], taskset,
            set(cpus_values),
        )
        union = completed["union_failures"]
        if not union or overhead == MAX_OVERHEAD:
            terminal = terminal_payload(stage, overhead, len(union))
            write_sealed_once(
                terminal_path, f"{SCHEMA}.terminal", terminal,
            )
            terminal = verify_terminal_campaign(
                result_dir, contract, set(cpus_values),
            )
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
        next_overhead = overhead + 1
        cohort = [
            (item["seed_index"], item["schedule"], item["K"]) for item in union
        ]
        next_tasks = build_tasks(next_overhead, cohort)
        write_sealed_once(
            stage_path(result_dir, next_overhead) / "manifest.json",
            f"{SCHEMA}.stage_manifest",
            manifest_payload(next_overhead, next_tasks),
        )
        overhead = next_overhead


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
    states: dict[tuple[int, str, int], dict[str, Any]] = {}
    arm_oh0_successes = {arm: 0 for arm in ARMS}
    arm_oh0_failure_minimums: dict[str, list[int | None]] = {
        arm: [] for arm in ARMS
    }
    union_oh0_minimums: dict[str, list[int | None]] = {arm: [] for arm in ARMS}
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
    overhead = 0
    oh0_metrics: dict[str, Any] | None = None
    while True:
        stage = stage_path(result_dir, overhead)
        tasks = load_manifest(stage, overhead)
        if previous_union is not None:
            expected = {
                (item["seed_index"], item["schedule"], item["K"])
                for item in previous_union
            }
            actual = {(task.seed_index, task.schedule, K) for task in tasks for K in task.ks}
            if actual != expected:
                die("nested overhead manifest is not the prior union-failure cohort")
        completed = load_sealed(
            stage / "complete.json", f"{SCHEMA}.stage_complete",
        )
        if overhead == 0:
            oh0_metrics = completed["metrics"]
            if completed["paired_cell_count"] != PAIRED_CELLS:
                die("OH0 reduction does not contain exactly 575,991 paired cells")
        next_states: dict[tuple[int, str, int], dict[str, Any]] = {}
        seen: set[tuple[int, str, int]] = set()
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
            (item["seed_index"], item["schedule"], item["K"])
            for item in previous_union
        }
        if set(states) != expected_union:
            die("computed union failures differ from the sealed stage receipt")
        if not previous_union or overhead == MAX_OVERHEAD:
            break
        overhead += 1
        if not stage_path(result_dir, overhead).exists():
            die("campaign is incomplete: next overhead manifest is missing")
    for state in states.values():
        finalize_state(state)
    minimums = {
        arm: [0] * arm_oh0_successes[arm] + arm_oh0_failure_minimums[arm]
        for arm in ARMS
    }
    if (any(len(minimums[arm]) != PAIRED_CELLS for arm in ARMS) or
            any(len(union_oh0_minimums[arm]) !=
                oh0_metrics["comparison"]["union_failures"] for arm in ARMS) or
            any(len(common_failure_oh0_minimums[arm]) !=
                oh0_metrics["comparison"]["both_fail"] for arm in ARMS)):
        die("minimum-overhead cohort accounting mismatch")
    summary = {
        "schema": f"{SCHEMA}.analysis",
        "coverage": {
            "K_min": K_MIN, "K_max": K_MAX, "unique_K": K_COUNT,
            "packet_trace_roots": len(SEEDS), "schedules": len(SCHEDULES),
            "paired_cells": PAIRED_CELLS, "outcomes_per_arm": PAIRED_CELLS,
            "arm_outcomes_at_OH0": TOTAL_ARM_OUTCOMES,
            "oh0_jobs": OH0_JOBS,
            "terminal_overhead": terminal["terminal_overhead"],
        },
        "oh0": oh0_metrics,
        "minimum_success_overhead": {
            arm: {
                "all_cells": overhead_report(minimums[arm]),
                "among_arm_OH0_failures": overhead_report(
                    arm_oh0_failure_minimums[arm]),
                "among_OH0_union_failures": overhead_report(
                    union_oh0_minimums[arm]),
                "among_common_OH0_failures": overhead_report(
                    common_failure_oh0_minimums[arm]),
            } for arm in ARMS
        },
        "interpretation": (
            "Deterministic named-seed/schedule census; weak-K counts and exact "
            "McNemar values are descriptive and are not rejection gates. "
            "A censored minimum means success was not observed through OH1024; "
            "it is unknown/>1024 and is never encoded as observed OH1025."
        ),
        "unresolved_union_failures_at_cap": len(previous_union or []),
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
