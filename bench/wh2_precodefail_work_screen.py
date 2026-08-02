#!/usr/bin/env python3
"""Run the fixed pure-GF256 WH2 exact-attempt work/rank screen.

This controller deliberately does not measure or select on time.  It invokes
the test-hook ``precodefail`` worker on the 360-cell development recovery
domain using the K-independent ``uniform-raw-v1`` construction seed basis,
validates every CSV field that identifies the requested experiment, and
publishes only rank and algebraic-work counters.

Source inspection establishes equivalence to the frozen stream when each
invocation has ``trials=1``, ``bb=2``, the frozen loss root as ``--seed``, and
``--paired-overhead-stream``.  Consequently precodefail's local trial and
overhead salts are zero.  Its IID path consumes that seed directly; the other
three paths pass it through ``BuildPacketSchedule``, whose sole seed transform
is the same xor-0x10fade used by ``Wh2FrozenTrace``.  The generator below also
checks the aggregate oracle in ``Wh2FrozenTraceTest.cpp`` before launch.

That is audited source equivalence, not runtime packet-ID attestation: the
precodefail CSV does not emit IDs or their hash.  Every artifact binds the
worker and all audited source files (including the raw-seed schedule), labels
trace hashes as frozen expectations, and states that the native trace ledger
remains authoritative.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
from decimal import Decimal, InvalidOperation
import hashlib
import math
import os
from pathlib import Path
import re
import selectors
import signal
import stat
import subprocess
import sys
import tempfile
import threading
import time
from dataclasses import dataclass
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_benchmark_contract as contract_api


RECORD_SCHEMA = "wirehair.wh2.precodefail-work-record.v1"
SUMMARY_SCHEMA = "wirehair.wh2.precodefail-work-summary.v1"
TRACE_IDENTITY_SCHEMA = "wirehair.wh2.precodefail-trace-identity.v1"
RESULT_NAME = "work-rank-results.jsonl"
SUMMARY_NAME = "work-rank-summary.json"
MAX_WALL_SECONDS = 7200.0
MAX_JOBS = 48
MAX_STDOUT_BYTES = 2 * 1024 * 1024
MAX_STDERR_BYTES = 64 * 1024
EXPECTED_RECORDS = 5760
EXPECTED_INVOCATIONS = 48
EXPECTED_FROZEN_TRACE_AGGREGATE_SHA256 = (
    "ec6f5ed3976e5e9664bb5db927c525a7"
    "3e8ab952d080492aaecd99c366db41a9"
)
EXPECTED_FROZEN_CANDIDATES = 5345254
EXPECTED_FROZEN_MAX_CANDIDATES = 128334
REPO_ROOT = Path(__file__).resolve().parent.parent
BENCH_SOURCE = REPO_ROOT / "codec" / "WirehairV2Bench.cpp"
RAW_SEED_SOURCE = REPO_ROOT / "codec" / "WirehairV2RawSeed.h"
FROZEN_TRACE_SOURCE = REPO_ROOT / "bench" / "Wh2FrozenTrace.cpp"
CONTROLLER_SOURCE = REPO_ROOT / "bench" / "wh2_precodefail_work_screen.py"
MASK32 = (1 << 32) - 1
MASK64 = (1 << 64) - 1
UINT_TOKEN = re.compile(r"(?:0|[1-9][0-9]*)\Z")
HEX_TOKEN = re.compile(r"0x(?:0|[1-9a-f][0-9a-f]*)\Z")

RAW_SEED_BASIS = "uniform-raw-v1"
RAW_SEED_SCHEDULE_SHA256 = (
    "90a98a3db207852dabdf5fb27573ef48b"
    "ce52e0228cee4e291d96fa44ed509a7"
)
RAW_PRECODE_BASE_SEED = 0x487468302aad7105
RAW_PACKET_BASE_SEED = 0x4ec72102
RAW_PRECODE_ATTEMPT_STRIDE = 0x9e3779b97f4a7c15
RAW_PACKET_ATTEMPT_STRIDE = 0x9e3779b9

K_VALUES = tuple(contract_api.EXPECTED_SHORT_K)
OVERHEADS = (0, 1, 2, 4)
STRATA = (
    ("iid", 100000),
    ("burst", 500000),
    ("adversarial", 500000),
    ("repair-only", 500000),
)
ROOTS = tuple(contract_api.EXPECTED_TRAINING_LOSS_ROOTS)
ARMS = (
    ("wirehair2_raw_d12_h11_periodic", 12, 11),
    ("wirehair2_raw_d12_h12_periodic", 12, 12),
    ("wirehair2_raw_d12_h13_periodic", 12, 13),
    ("wirehair2_raw_d13_h12_periodic", 13, 12),
)

CSV_HEADER = (
    "N", "bb", "heavy_family", "mix_count", "overhead", "trials",
    "success", "rank_fail", "error", "fail_rate", "inact_mu",
    "inact_max", "binary_def_mu", "binary_def_max", "heavy_gain_mu",
    "heavy_gain_min", "heavy_shortfall", "solve_ms_mu", "build_ms_mu",
    "peel_ms_mu", "project_ms_mu", "residual_ms_mu", "backsub_ms_mu",
    "seed_attempt", "block_xors_mu", "block_muladds_mu",
    "first_rank_fail", "binary_def_hist", "heavy_gain_hist",
    "failure_trials", "active_packet_peel_seed_xor", "precode_attempt",
    "packet_attempt", "attempt_mode", "construction_seed_basis",
    "seed_schedule_sha256", "effective_precode_seed",
    "effective_packet_seed",
)
METADATA_KEYS = frozenset((
    "trials", "threads", "loss", "seed", "source_hits_override",
    "packet_peel_seed_xor", "binary_dense_rows_override",
    "gf256_heavy_rows_override", "odd_packet_peel_seed_xor",
    "packet_row_seed_multiplier", "packet_row_seed_avalanche",
    "seed_block_bytes_override", "overhead_stream", "full_payload_solve",
    "schedule", "exact_attempt_mode", "exact_precode_attempt",
    "exact_packet_attempt", "construction_seed_basis",
    "seed_schedule_sha256",
))


class WorkScreenError(RuntimeError):
    """The fixed work/rank screen cannot be published safely."""


def fail(message: str) -> None:
    raise WorkScreenError(message)


@dataclass(frozen=True)
class Invocation:
    ordinal: int
    arm: str
    dense_rows: int
    heavy_rows: int
    trial: int
    root: str
    schedule: str
    loss_ppm: int

    def argv(self, worker: Path) -> List[str]:
        return [
            str(worker), "precodefail",
            "--N", ",".join(str(value) for value in K_VALUES),
            "--bb-list", "2",
            "--overhead", ",".join(str(value) for value in OVERHEADS),
            "--trials", "1", "--threads", "1",
            "--loss", _loss_arg(self.loss_ppm),
            "--seed", self.root,
            "--schedule", self.schedule,
            "--heavy-family", "periodic", "--mix-count", "3",
            "--binary-dense-rows", str(self.dense_rows),
            "--gf256-heavy-rows", str(self.heavy_rows),
            "--paired-overhead-stream", "--full-payload-solve",
            "--exact-precode-attempt", str(self.trial),
            "--exact-packet-attempt", str(self.trial),
            "--construction-seed-basis", RAW_SEED_BASIS,
        ]

    def identity(self) -> Mapping[str, Any]:
        return {
            "arm": self.arm,
            "binary_dense_rows": self.dense_rows,
            "gf256_heavy_rows": self.heavy_rows,
            "heavy_family": "periodic",
            "loss_ppm": self.loss_ppm,
            "loss_seed": self.root,
            "mix_count": 3,
            "attempt_mode": "exact",
            "packet_attempt": self.trial,
            "precode_attempt": self.trial,
            "schedule": self.schedule,
            "seed_basis": RAW_SEED_BASIS,
            "seed_schedule_sha256": RAW_SEED_SCHEDULE_SHA256,
            "effective_precode_seed": _effective_precode_seed(self.trial),
            "effective_packet_seed": _effective_packet_seed(self.trial),
            "trial": self.trial,
        }


@dataclass(frozen=True)
class FrozenTrace:
    packet_ids: Tuple[int, ...]
    attempted_candidates: int
    trace_sha256: str
    prefix_attempted_candidates: Tuple[int, int, int, int]


@dataclass(frozen=True)
class InvocationResult:
    invocation: Invocation
    command_sha256: str
    stdout_sha256: str
    rows: Tuple[Mapping[str, Any], ...]


class _SplitMix64:
    def __init__(self, seed: int) -> None:
        self.state = seed & MASK64

    def next(self) -> int:
        self.state = (self.state + 0x9e3779b97f4a7c15) & MASK64
        value = self.state
        value = ((value ^ (value >> 30)) * 0xbf58476d1ce4e5b9) & MASK64
        value = ((value ^ (value >> 27)) * 0x94d049bb133111eb) & MASK64
        return (value ^ (value >> 31)) & MASK64

    def unit(self) -> float:
        return float(self.next() >> 11) * (1.0 / 9007199254740992.0)


class _ProcessRegistry:
    def __init__(self) -> None:
        self.lock = threading.Lock()
        self.processes: Dict[int, subprocess.Popen] = {}
        self.cancelled = threading.Event()

    def add(self, process: subprocess.Popen) -> None:
        with self.lock:
            if self.cancelled.is_set():
                _signal_process_group(process, signal.SIGKILL)
                try:
                    process.wait(timeout=1.0)
                except (OSError, subprocess.TimeoutExpired):
                    pass
                fail("work/rank campaign was cancelled")
            self.processes[process.pid] = process

    def remove(self, process: subprocess.Popen) -> None:
        with self.lock:
            self.processes.pop(process.pid, None)

    def cancel(self) -> None:
        self.cancelled.set()
        with self.lock:
            processes = list(self.processes.values())
        for process in processes:
            _signal_process_group(process, signal.SIGKILL)


def _signal_process_group(process: subprocess.Popen, sig: int) -> None:
    if process.poll() is not None:
        return
    try:
        os.killpg(process.pid, sig)
    except ProcessLookupError:
        pass
    except OSError:
        try:
            process.send_signal(sig)
        except OSError:
            pass


def _kill_and_reap(process: subprocess.Popen) -> None:
    _signal_process_group(process, signal.SIGKILL)
    try:
        process.wait(timeout=2.0)
    except (OSError, subprocess.TimeoutExpired):
        _signal_process_group(process, signal.SIGKILL)
        try:
            process.wait(timeout=2.0)
        except (OSError, subprocess.TimeoutExpired):
            pass


def _capture_bounded(
        process: subprocess.Popen, deadline: float,
        registry: _ProcessRegistry) -> Tuple[bytes, bytes]:
    if process.stdout is None or process.stderr is None:
        registry.cancel()
        fail("precodefail worker pipes were not created")
    selector = selectors.DefaultSelector()
    streams = {
        process.stdout.fileno(): ("stdout", process.stdout, MAX_STDOUT_BYTES),
        process.stderr.fileno(): ("stderr", process.stderr, MAX_STDERR_BYTES),
    }
    buffers = {"stdout": bytearray(), "stderr": bytearray()}
    try:
        for descriptor in streams:
            selector.register(descriptor, selectors.EVENT_READ)
        while selector.get_map():
            remaining = deadline - time.monotonic()
            if remaining <= 0.0:
                registry.cancel()
                fail("work/rank campaign hard wall expired")
            events = selector.select(timeout=remaining)
            if not events:
                registry.cancel()
                fail("work/rank campaign hard wall expired")
            for key, _mask in events:
                descriptor = key.fd
                name, _stream, cap = streams[descriptor]
                try:
                    data = os.read(descriptor, 65536)
                except OSError as exc:
                    registry.cancel()
                    fail("cannot read precodefail {}: {}".format(name, exc))
                if not data:
                    selector.unregister(descriptor)
                    continue
                if len(buffers[name]) + len(data) > cap:
                    registry.cancel()
                    fail("precodefail {} exceeded its fixed byte cap".format(
                        name))
                buffers[name].extend(data)
        remaining = deadline - time.monotonic()
        if remaining <= 0.0:
            registry.cancel()
            fail("work/rank campaign hard wall expired")
        try:
            process.wait(timeout=remaining)
        except subprocess.TimeoutExpired:
            registry.cancel()
            fail("work/rank campaign hard wall expired")
        return bytes(buffers["stdout"]), bytes(buffers["stderr"])
    finally:
        selector.close()
        process.stdout.close()
        process.stderr.close()


def _loss_arg(loss_ppm: int) -> str:
    if loss_ppm == 100000:
        return "0.1"
    if loss_ppm == 500000:
        return "0.5"
    fail("unsupported frozen loss stratum {}".format(loss_ppm))
    return ""  # Unreachable; keeps Python 3.8 type checkers satisfied.


def _canonical(value: Any) -> str:
    return contract_api.canonical_json(value)


def _sha256_json(value: Any) -> str:
    return hashlib.sha256(_canonical(value).encode("utf-8")).hexdigest()


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        while True:
            block = source.read(1024 * 1024)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def _effective_precode_seed(attempt: int) -> str:
    value = (RAW_PRECODE_BASE_SEED +
             attempt * RAW_PRECODE_ATTEMPT_STRIDE) & MASK64
    return "0x{:016x}".format(value)


def _effective_packet_seed(attempt: int) -> str:
    value = (RAW_PACKET_BASE_SEED +
             attempt * RAW_PACKET_ATTEMPT_STRIDE) & MASK32
    return "0x{:08x}".format(value)


def _source_provenance(worker: Path) -> Mapping[str, Any]:
    for path in (
            BENCH_SOURCE, RAW_SEED_SOURCE, FROZEN_TRACE_SOURCE,
            CONTROLLER_SOURCE):
        try:
            info = path.stat()
        except OSError as exc:
            fail("cannot inspect audited source {}: {}".format(path, exc))
        if not stat.S_ISREG(info.st_mode):
            fail("audited source must be a regular file: {}".format(path))
    try:
        head = subprocess.run(
            ["git", "rev-parse", "--verify", "HEAD"], cwd=str(REPO_ROOT),
            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, timeout=10.0, check=False)
        status = subprocess.run(
            ["git", "status", "--porcelain=v1", "--untracked-files=no"],
            cwd=str(REPO_ROOT), stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=10.0,
            check=False)
    except (OSError, subprocess.TimeoutExpired) as exc:
        fail("cannot bind source git provenance: {}".format(exc))
    if head.returncode != 0 or status.returncode != 0:
        fail("cannot bind source git provenance")
    try:
        commit = head.stdout.decode("ascii").strip()
        status.stdout.decode("utf-8")
    except UnicodeDecodeError:
        fail("source git provenance is not valid text")
    if contract_api.GIT_COMMIT.fullmatch(commit) is None:
        fail("source git commit is malformed")
    return {
        "source_git_commit": commit,
        "tracked_worktree_clean": not bool(status.stdout),
        "tracked_status_sha256": hashlib.sha256(status.stdout).hexdigest(),
        "worker_binary_sha256": _sha256_file(worker),
        "controller_source_sha256": _sha256_file(CONTROLLER_SOURCE),
        "wirehair_v2_bench_source_sha256": _sha256_file(BENCH_SOURCE),
        "wirehair_v2_raw_seed_source_sha256": _sha256_file(RAW_SEED_SOURCE),
        "wh2_frozen_trace_source_sha256": _sha256_file(FROZEN_TRACE_SOURCE),
    }


def _packet_hash(packet_ids: Sequence[int]) -> str:
    digest = hashlib.sha256()
    for packet_id in packet_ids:
        digest.update(int(packet_id).to_bytes(4, "little"))
    return digest.hexdigest()


def _frozen_trace(cell: Mapping[str, Any]) -> FrozenTrace:
    K = cell["K"]
    block_bytes = cell["block_bytes"]
    loss_ppm = cell["loss_ppm"]
    schedule = cell["schedule"]
    loss_seed = int(cell["loss_seed"], 16)
    delivered_count = K + 4
    candidate_limit = delivered_count * 256 + 65536
    trace_seed = loss_seed ^ ((K * 0x9e3779b97f4a7c15) & MASK64) ^ \
        ((block_bytes * 0xbf58476d1ce4e5b9) & MASK64)
    rng = _SplitMix64(trace_seed if schedule == "iid" else
                      trace_seed ^ 0x10fade)
    loss_rate = float(loss_ppm) / 1000000.0
    burst_remaining = 0
    packet_ids: List[int] = []
    prefix_attempts: Dict[int, int] = {}
    attempted = 0
    while len(packet_ids) < delivered_count and attempted < candidate_limit:
        candidate = attempted
        attempted += 1
        if schedule == "adversarial":
            packet_id = (MASK32 - ((candidate * 2) & MASK32)) & MASK32
        elif schedule == "repair-only":
            packet_id = (K + candidate) & MASK32
        else:
            packet_id = candidate & MASK32

        if schedule == "burst":
            if burst_remaining:
                burst_remaining -= 1
                drop = True
            else:
                start_probability = loss_rate / (8.0 - 7.0 * loss_rate)
                drop = rng.unit() < start_probability
                if drop:
                    burst_remaining = 7
        else:
            drop = rng.unit() < loss_rate
        if not drop:
            packet_ids.append(packet_id)
            delivered_overhead = len(packet_ids) - K
            if delivered_overhead in OVERHEADS:
                prefix_attempts[delivered_overhead] = attempted
    if len(packet_ids) != delivered_count:
        fail("frozen packet generator hit its candidate cap")
    if set(prefix_attempts) != set(OVERHEADS):
        fail("frozen packet generator omitted a nested threshold")
    for overhead in OVERHEADS:
        if prefix_attempts[overhead] > (K + overhead) * 256 + 65536:
            fail("precodefail's per-prefix candidate cap is not frozen-equivalent")
    return FrozenTrace(
        tuple(packet_ids), attempted, _packet_hash(packet_ids),
        tuple(prefix_attempts[value] for value in OVERHEADS))


def build_frozen_domain(contract: Mapping[str, Any]) \
        -> Tuple[Mapping[Tuple[int, str, int], Mapping[str, Any]],
                 Mapping[Tuple[int, str, int], FrozenTrace], Mapping[str, Any]]:
    cells = list(contract_api.iter_recovery_cells(contract, "development"))
    if len(cells) != 360:
        fail("development recovery domain is not the frozen 360 cells")
    cell_map: Dict[Tuple[int, str, int], Mapping[str, Any]] = {}
    trace_map: Dict[Tuple[int, str, int], FrozenTrace] = {}
    aggregate = hashlib.sha256()
    total_candidates = 0
    maximum_candidates = 0
    for cell in cells:
        key = (cell["trial"], cell["schedule"], cell["K"])
        if key in cell_map:
            fail("development recovery domain contains a duplicate cell")
        trace = _frozen_trace(cell)
        aggregate.update((trace.trace_sha256 + ":" +
                          str(trace.attempted_candidates) + "\n").encode(
                              "ascii"))
        total_candidates += trace.attempted_candidates
        maximum_candidates = max(maximum_candidates,
                                 trace.attempted_candidates)
        cell_map[key] = cell
        trace_map[key] = trace
    aggregate_sha256 = aggregate.hexdigest()
    if (aggregate_sha256 != EXPECTED_FROZEN_TRACE_AGGREGATE_SHA256 or
            total_candidates != EXPECTED_FROZEN_CANDIDATES or
            maximum_candidates != EXPECTED_FROZEN_MAX_CANDIDATES):
        fail("independent packet generator disagrees with Wh2FrozenTrace")
    proof = {
        "schema": TRACE_IDENTITY_SCHEMA,
        "status": "audited_source_equivalence",
        "runtime_packet_ids_observed": False,
        "authoritative_trace_evidence": "native_frozen_trace_ledger",
        "runtime_fingerprint_seam": (
            "precodefail must emit the SHA-256 of each delivered-ID prefix"
        ),
        "cell_count": len(cells),
        "root_stratum_case_count": len(ROOTS) * len(STRATA),
        "nested_prefix_count": len(cells) * len(OVERHEADS),
        "trace_aggregate_sha256": aggregate_sha256,
        "total_attempted_candidates": total_candidates,
        "maximum_attempted_candidates": maximum_candidates,
        "precodefail_local_trial": 0,
        "precodefail_overhead_salt": 0,
        "precodefail_block_bytes": 2,
        "non_iid_seed_xor": "0x10fade",
    }
    return cell_map, trace_map, proof


def make_invocations() -> Tuple[Invocation, ...]:
    values: List[Invocation] = []
    for arm, dense_rows, heavy_rows in ARMS:
        for trial, root in enumerate(ROOTS):
            for schedule, loss_ppm in STRATA:
                values.append(Invocation(
                    len(values), arm, dense_rows, heavy_rows, trial, root,
                    schedule, loss_ppm))
    if len(values) != EXPECTED_INVOCATIONS:
        fail("internal invocation cardinality is inconsistent")
    return tuple(values)


def _parse_uint(text: str, field: str, maximum: int = MASK64) -> int:
    if UINT_TOKEN.fullmatch(text) is None:
        fail("{} is not a canonical unsigned integer".format(field))
    value = int(text)
    if value > maximum:
        fail("{} exceeds its supported range".format(field))
    return value


def _parse_integral_decimal(text: str, field: str) -> int:
    try:
        value = Decimal(text)
    except InvalidOperation:
        fail("{} is not a decimal number".format(field))
    if not value.is_finite() or value < 0 or value != value.to_integral_value():
        fail("{} is not a nonnegative integral counter".format(field))
    integer = int(value)
    if integer > MASK64:
        fail("{} exceeds its supported range".format(field))
    return integer


def _parse_nonnegative_decimal(text: str, field: str) -> Decimal:
    try:
        value = Decimal(text)
    except InvalidOperation:
        fail("{} is not a decimal number".format(field))
    if not value.is_finite() or value < 0:
        fail("{} is not a finite nonnegative value".format(field))
    return value


def _parse_metadata(line: str, invocation: Invocation) -> Mapping[str, str]:
    prefix = "# precodefail: "
    if not line.startswith(prefix):
        fail("precodefail metadata line is missing")
    metadata: Dict[str, str] = {}
    for token in line[len(prefix):].split(" "):
        if token.count("=") != 1:
            fail("precodefail metadata token is malformed")
        key, value = token.split("=", 1)
        if not key or not value or key in metadata:
            fail("precodefail metadata has a duplicate or empty field")
        metadata[key] = value
    if set(metadata) != set(METADATA_KEYS):
        fail("precodefail metadata schema changed")
    expected = {
        "trials": "1", "threads": "1",
        "source_hits_override": "0", "packet_peel_seed_xor": "0x0",
        "binary_dense_rows_override": str(invocation.dense_rows),
        "gf256_heavy_rows_override": str(invocation.heavy_rows),
        "odd_packet_peel_seed_xor": "0x0",
        "packet_row_seed_multiplier": "0x1",
        "packet_row_seed_avalanche": "0",
        "seed_block_bytes_override": "0", "overhead_stream": "paired",
        "full_payload_solve": "1", "schedule": invocation.schedule,
        "exact_attempt_mode": "1",
        "exact_precode_attempt": str(invocation.trial),
        "exact_packet_attempt": str(invocation.trial),
        "construction_seed_basis": RAW_SEED_BASIS,
        "seed_schedule_sha256": RAW_SEED_SCHEDULE_SHA256,
    }
    for key, expected_value in expected.items():
        if metadata[key] != expected_value:
            fail("precodefail metadata {} disagrees with command".format(key))
    if HEX_TOKEN.fullmatch(metadata["seed"]) is None or \
            int(metadata["seed"], 16) != int(invocation.root, 16):
        fail("precodefail metadata seed disagrees with frozen root")
    expected_loss = format(float(invocation.loss_ppm) / 1000000.0, ".17g")
    if metadata["loss"] != expected_loss:
        fail("precodefail metadata loss disagrees with frozen stratum")
    return metadata


def _parse_csv_line(line: str) -> List[str]:
    try:
        rows = list(csv.reader([line], strict=True))
    except csv.Error as exc:
        fail("precodefail CSV is malformed: {}".format(exc))
    if len(rows) != 1 or len(rows[0]) != len(CSV_HEADER):
        fail("precodefail CSV row is not {} aligned fields".format(
            len(CSV_HEADER)))
    if line != ",".join(rows[0]):
        fail("precodefail CSV row is not in its unquoted canonical form")
    return rows[0]


def parse_invocation_output(
        invocation: Invocation, command: Sequence[str], stdout: bytes,
        stderr: bytes,
        cell_map: Mapping[Tuple[int, str, int], Mapping[str, Any]],
        trace_map: Mapping[Tuple[int, str, int], FrozenTrace],
        worker_binary_sha256: str) -> InvocationResult:
    if len(stdout) > MAX_STDOUT_BYTES or len(stderr) > MAX_STDERR_BYTES:
        fail("precodefail output exceeded its fixed byte cap")
    if stderr:
        fail("precodefail emitted unexpected stderr")
    if not stdout.endswith(b"\n") or b"\r" in stdout:
        fail("precodefail stdout is not newline terminated")
    try:
        text = stdout.decode("ascii")
    except UnicodeDecodeError:
        fail("precodefail stdout is not ASCII")
    lines = text.splitlines()
    expected_rows = len(K_VALUES) * len(OVERHEADS)
    if len(lines) != expected_rows + 2:
        fail("precodefail invocation did not emit exactly 120 rows")
    _parse_metadata(lines[0], invocation)
    if tuple(_parse_csv_line(lines[1])) != CSV_HEADER:
        fail("precodefail CSV header changed")

    command_identity = {
        "argv": list(command[1:]),
        "invocation": invocation.identity(),
    }
    command_sha256 = _sha256_json(command_identity)
    stdout_sha256 = hashlib.sha256(stdout).hexdigest()
    expected_precode_seed = _effective_precode_seed(invocation.trial)
    expected_packet_seed = _effective_packet_seed(invocation.trial)
    seen = set()
    parsed: List[Mapping[str, Any]] = []
    for line_number, line in enumerate(lines[2:], 3):
        values = _parse_csv_line(line)
        row = dict(zip(CSV_HEADER, values))
        K = _parse_uint(row["N"], "N", 64000)
        overhead = _parse_uint(row["overhead"], "overhead", 4)
        key = (K, overhead)
        if K not in K_VALUES or overhead not in OVERHEADS or key in seen:
            fail("precodefail row {} has an unexpected or duplicate key".format(
                line_number))
        seen.add(key)
        exact_text = {
            "bb": "2", "heavy_family": "periodic", "mix_count": "3",
            "trials": "1", "seed_attempt": "",
            "active_packet_peel_seed_xor": "0x0",
            "precode_attempt": str(invocation.trial),
            "packet_attempt": str(invocation.trial),
            "attempt_mode": "exact",
            "construction_seed_basis": RAW_SEED_BASIS,
            "seed_schedule_sha256": RAW_SEED_SCHEDULE_SHA256,
            "effective_precode_seed": expected_precode_seed,
            "effective_packet_seed": expected_packet_seed,
        }
        for field, expected in exact_text.items():
            if row[field] != expected:
                fail("precodefail row {} field {} disagrees with the fixed "
                     "command".format(line_number, field))

        success = _parse_uint(row["success"], "success", 1)
        rank_fail = _parse_uint(row["rank_fail"], "rank_fail", 1)
        errors = _parse_uint(row["error"], "error", 1)
        if success + rank_fail + errors != 1 or errors != 0:
            fail("precodefail row has an error or invalid terminal counts")
        fail_rate = _parse_nonnegative_decimal(row["fail_rate"], "fail_rate")
        if fail_rate != Decimal(rank_fail):
            fail("precodefail fail_rate disagrees with terminal counts")
        inactivated = _parse_integral_decimal(row["inact_mu"], "inact_mu")
        if _parse_uint(row["inact_max"], "inact_max") != inactivated:
            fail("precodefail inactivation mean/max disagree for one trial")
        binary_deficiency = _parse_integral_decimal(
            row["binary_def_mu"], "binary_def_mu")
        if _parse_uint(row["binary_def_max"], "binary_def_max") != \
                binary_deficiency:
            fail("precodefail binary deficiency mean/max disagree")
        gf256_rank_gain = _parse_integral_decimal(
            row["heavy_gain_mu"], "heavy_gain_mu")
        if _parse_uint(row["heavy_gain_min"], "heavy_gain_min") != \
                gf256_rank_gain:
            fail("precodefail GF256 rank gain mean/min disagree")
        heavy_shortfall = _parse_uint(
            row["heavy_shortfall"], "heavy_shortfall", 1)
        if (binary_deficiency > inactivated or
                gf256_rank_gain > binary_deficiency or
                bool(success) != (gf256_rank_gain == binary_deficiency) or
                bool(rank_fail) != (gf256_rank_gain < binary_deficiency)):
            fail("precodefail rank counters disagree with terminal result")
        expected_shortfall = int(
            bool(rank_fail) and binary_deficiency <= invocation.heavy_rows and
            gf256_rank_gain < binary_deficiency)
        if heavy_shortfall != expected_shortfall:
            fail("precodefail heavy shortfall disagrees with rank counters")
        block_xors = _parse_integral_decimal(
            row["block_xors_mu"], "block_xors_mu")
        gf256_muladds = _parse_integral_decimal(
            row["block_muladds_mu"], "block_muladds_mu")
        for field in (
                "solve_ms_mu", "build_ms_mu", "peel_ms_mu",
                "project_ms_mu", "residual_ms_mu", "backsub_ms_mu"):
            _parse_nonnegative_decimal(row[field], field)
        expected_first = "0" if rank_fail else "-1"
        expected_failures = "0" if rank_fail else ""
        if (row["first_rank_fail"] != expected_first or
                row["failure_trials"] != expected_failures or
                row["binary_def_hist"] != "{}:1".format(binary_deficiency) or
                row["heavy_gain_hist"] != "{}:1".format(gf256_rank_gain)):
            fail("precodefail one-trial diagnostics are internally inconsistent")

        frozen_key = (invocation.trial, invocation.schedule, K)
        cell = cell_map.get(frozen_key)
        trace = trace_map.get(frozen_key)
        if cell is None or trace is None:
            fail("precodefail row does not map to a frozen recovery cell")
        if (cell["loss_seed"] != invocation.root or
                cell["loss_ppm"] != invocation.loss_ppm or
                cell["base_seed_attempt"] != invocation.trial or
                cell["block_bytes"] != 2):
            fail("precodefail invocation disagrees with its frozen cell")
        prefix = trace.packet_ids[:K + overhead]
        parsed.append({
            "schema": RECORD_SCHEMA,
            "arm": invocation.arm,
            "binary_dense_rows": invocation.dense_rows,
            "gf256_heavy_rows": invocation.heavy_rows,
            "heavy_family": "periodic",
            "mix_count": 3,
            "phase": "development",
            "band": cell["band"],
            "K": K,
            "block_bytes": 2,
            "loss_ppm": invocation.loss_ppm,
            "schedule": invocation.schedule,
            "trial": invocation.trial,
            "base_seed_attempt": invocation.trial,
            "loss_seed": invocation.root,
            "overhead": overhead,
            "precode_attempt": invocation.trial,
            "packet_attempt": invocation.trial,
            "attempt_mode": "exact",
            "seed_basis": RAW_SEED_BASIS,
            "seed_schedule_sha256": RAW_SEED_SCHEDULE_SHA256,
            "effective_precode_seed": expected_precode_seed,
            "effective_packet_seed": expected_packet_seed,
            "cell_sha256": _sha256_json(cell),
            "frozen_trace_sha256": trace.trace_sha256,
            "frozen_packet_prefix_sha256": _packet_hash(prefix),
            "packet_count": len(prefix),
            "frozen_trace_attempted_candidates": trace.attempted_candidates,
            "frozen_prefix_attempted_candidates":
                trace.prefix_attempted_candidates[OVERHEADS.index(overhead)],
            "success": success,
            "rank_fail": rank_fail,
            "error": errors,
            "inactivated_columns": inactivated,
            "binary_deficiency": binary_deficiency,
            "gf256_rank_gain": gf256_rank_gain,
            "heavy_shortfall": heavy_shortfall,
            "block_xors": block_xors,
            "gf256_muladds": gf256_muladds,
            "command_sha256": command_sha256,
            "worker_stdout_sha256": stdout_sha256,
            "worker_binary_sha256": worker_binary_sha256,
        })
    if seen != set((K, overhead) for K in K_VALUES for overhead in OVERHEADS):
        fail("precodefail invocation omitted a frozen K/overhead row")
    return InvocationResult(
        invocation, command_sha256, stdout_sha256, tuple(parsed))


def _run_invocation(
        invocation: Invocation, worker: Path, deadline: float,
        registry: _ProcessRegistry,
        cell_map: Mapping[Tuple[int, str, int], Mapping[str, Any]],
        trace_map: Mapping[Tuple[int, str, int], FrozenTrace],
        worker_binary_sha256: str) -> InvocationResult:
    if registry.cancelled.is_set():
        fail("work/rank campaign was cancelled")
    command = invocation.argv(worker)
    try:
        process = subprocess.Popen(
            command, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, start_new_session=True)
    except OSError as exc:
        registry.cancel()
        fail("cannot start precodefail worker: {}".format(exc))
    registry.add(process)
    try:
        stdout, stderr = _capture_bounded(process, deadline, registry)
        if process.returncode != 0:
            registry.cancel()
            detail = stderr[:MAX_STDERR_BYTES].decode("utf-8", "replace").strip()
            fail("precodefail worker exited {}{}".format(
                process.returncode, ": " + detail if detail else ""))
        return parse_invocation_output(
            invocation, command, stdout, stderr, cell_map, trace_map,
            worker_binary_sha256)
    except BaseException:
        registry.cancel()
        _kill_and_reap(process)
        raise
    finally:
        registry.remove(process)


def _domain_sha256() -> str:
    digest = hashlib.sha256()
    for arm, dense_rows, heavy_rows in ARMS:
        for trial, root in enumerate(ROOTS):
            for schedule, loss_ppm in STRATA:
                for K in K_VALUES:
                    for overhead in OVERHEADS:
                        identity = {
                            "K": K, "arm": arm,
                            "binary_dense_rows": dense_rows,
                            "block_bytes": 2,
                            "gf256_heavy_rows": heavy_rows,
                            "heavy_family": "periodic",
                            "loss_ppm": loss_ppm, "loss_seed": root,
                            "mix_count": 3, "overhead": overhead,
                            "attempt_mode": "exact",
                            "packet_attempt": trial,
                            "precode_attempt": trial,
                            "schedule": schedule,
                            "seed_basis": RAW_SEED_BASIS,
                            "seed_schedule_sha256":
                                RAW_SEED_SCHEDULE_SHA256,
                            "effective_precode_seed":
                                _effective_precode_seed(trial),
                            "effective_packet_seed":
                                _effective_packet_seed(trial),
                            "trial": trial,
                        }
                        digest.update(_canonical(identity).encode("utf-8"))
                        digest.update(b"\n")
    return digest.hexdigest()


def _aggregate(records: Sequence[Mapping[str, Any]]) -> List[Mapping[str, Any]]:
    control_arm = "wirehair2_raw_d12_h12_periodic"
    control_rows = {
        (row["trial"], row["schedule"], row["K"], row["overhead"]): row
        for row in records if row["arm"] == control_arm
    }
    if len(control_rows) != 1440:
        fail("control arm aggregation has the wrong cardinality")
    summaries: List[Mapping[str, Any]] = []
    for arm, dense_rows, heavy_rows in ARMS:
        arm_rows = [row for row in records if row["arm"] == arm]
        arm_index = {
            (row["trial"], row["schedule"], row["K"], row["overhead"]): row
            for row in arm_rows
        }
        if len(arm_rows) != 1440 or len(arm_index) != 1440:
            fail("arm aggregation has the wrong cardinality")
        by_overhead = []
        for overhead in OVERHEADS:
            rows = [row for row in arm_rows if row["overhead"] == overhead]
            if len(rows) != 360:
                fail("arm/overhead aggregation has the wrong cardinality")
            by_overhead.append({
                "overhead": overhead,
                "rows": len(rows),
                "successes": sum(row["success"] for row in rows),
                "rank_failures": sum(row["rank_fail"] for row in rows),
                "errors": sum(row["error"] for row in rows),
                "inactivated_columns_sum": sum(
                    row["inactivated_columns"] for row in rows),
                "binary_deficiency_sum": sum(
                    row["binary_deficiency"] for row in rows),
                "gf256_rank_gain_sum": sum(
                    row["gf256_rank_gain"] for row in rows),
                "block_xors_sum": sum(row["block_xors"] for row in rows),
                "gf256_muladds_sum": sum(
                    row["gf256_muladds"] for row in rows),
            })

        weak_units = []
        multiplicities: Dict[int, int] = {}
        for trial in range(3):
            for K in K_VALUES:
                failures = [
                    arm_index[(trial, schedule, K, 0)]
                    for schedule, _ in STRATA
                    if arm_index[(trial, schedule, K, 0)]["rank_fail"] == 1
                ]
                if not failures:
                    continue
                multiplicity = len(failures)
                multiplicities[multiplicity] = \
                    multiplicities.get(multiplicity, 0) + 1
                weak_units.append({
                    "K": K,
                    "block_bytes": 2,
                    "precode_attempt": trial,
                    "packet_attempt": trial,
                    "loss_seed": ROOTS[trial],
                    "failed_strata": multiplicity,
                    "failed_cells": [{
                        "cell_sha256": row["cell_sha256"],
                        "frozen_trace_sha256": row["frozen_trace_sha256"],
                        "loss_ppm": row["loss_ppm"],
                        "schedule": row["schedule"],
                    } for row in failures],
                })

        repairs = []
        introductions = []
        for key in sorted(arm_index):
            row = arm_index[key]
            control = control_rows[key]
            if row["rank_fail"] == control["rank_fail"]:
                continue
            transition = {
                "K": row["K"],
                "cell_sha256": row["cell_sha256"],
                "frozen_trace_sha256": row["frozen_trace_sha256"],
                "overhead": row["overhead"],
                "schedule": row["schedule"],
                "trial": row["trial"],
                "record_ordinal": row["ordinal"],
                "control_record_ordinal": control["ordinal"],
            }
            if control["rank_fail"] and not row["rank_fail"]:
                repairs.append(transition)
            elif row["rank_fail"] and not control["rank_fail"]:
                introductions.append(transition)
            else:
                fail("rank transition classification is inconsistent")
        summaries.append({
            "arm": arm,
            "binary_dense_rows": dense_rows,
            "gf256_heavy_rows": heavy_rows,
            "heavy_family": "periodic",
            "mix_count": 3,
            "rows": len(arm_rows),
            "by_overhead": by_overhead,
            "overhead_zero_weak_cell_count": sum(
                item["failed_strata"] for item in weak_units),
            "overhead_zero_weak_unit_count": len(weak_units),
            "weak_unit_multiplicity": [
                {"failed_strata": value, "unit_count": multiplicities[value]}
                for value in sorted(multiplicities)
            ],
            "weak_units": weak_units,
            "transitions_vs_control": {
                "control_arm": control_arm,
                "repair_count": len(repairs),
                "introduction_count": len(introductions),
                "repairs": repairs,
                "introductions": introductions,
            },
        })
    return summaries


def _atomic_write(path: Path, data: bytes) -> None:
    descriptor, raw_path = tempfile.mkstemp(
        prefix="." + path.name + ".", suffix=".tmp", dir=str(path.parent))
    staged = Path(raw_path)
    linked = False
    staged_info: Optional[os.stat_result] = None
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(data)
            output.flush()
            os.fsync(output.fileno())
        staged_info = staged.stat()
        try:
            os.link(str(staged), str(path))
            linked = True
        except FileExistsError:
            fail("refusing to replace existing artifact {}".format(path))
        directory = os.open(
            str(path.parent), os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    except BaseException:
        if linked and staged_info is not None:
            try:
                destination_info = os.lstat(str(path))
                if (destination_info.st_dev == staged_info.st_dev and
                        destination_info.st_ino == staged_info.st_ino):
                    os.unlink(str(path))
            except OSError:
                pass
        raise
    finally:
        try:
            staged.unlink()
        except FileNotFoundError:
            pass


def run_campaign(
        worker: Path, output_dir: Path, jobs: int,
        deadline_seconds: float) -> Mapping[str, Any]:
    if type(jobs) is not int or not 1 <= jobs <= MAX_JOBS:
        fail("jobs must be in [1,{}]".format(MAX_JOBS))
    if (not math.isfinite(deadline_seconds) or deadline_seconds <= 0.0 or
            deadline_seconds > MAX_WALL_SECONDS):
        fail("deadline must be finite and in (0,{}]".format(
            int(MAX_WALL_SECONDS)))
    worker = worker.resolve(strict=True)
    info = worker.stat()
    if not stat.S_ISREG(info.st_mode) or not os.access(str(worker), os.X_OK):
        fail("worker must be an executable regular file")
    output_dir.mkdir(parents=True, exist_ok=True)
    if output_dir.is_symlink() or not output_dir.is_dir():
        fail("output directory must be a real directory")
    result_path = output_dir / RESULT_NAME
    summary_path = output_dir / SUMMARY_NAME
    if (result_path.exists() or result_path.is_symlink() or
            summary_path.exists() or summary_path.is_symlink()):
        fail("refusing to replace an existing work/rank artifact")

    contract = contract_api.load_contract()
    cell_map, trace_map, trace_identity = build_frozen_domain(contract)
    source_provenance = _source_provenance(worker)
    source_provenance_sha256 = _sha256_json(source_provenance)
    worker_binary_sha256 = source_provenance["worker_binary_sha256"]
    invocations = make_invocations()
    registry = _ProcessRegistry()
    deadline = time.monotonic() + deadline_seconds
    results: Dict[int, InvocationResult] = {}
    executor: Optional[concurrent.futures.ThreadPoolExecutor] = None
    futures: Dict[concurrent.futures.Future, int] = {}
    try:
        executor = concurrent.futures.ThreadPoolExecutor(max_workers=jobs)
        for invocation in invocations:
            future = executor.submit(
                _run_invocation, invocation, worker, deadline, registry,
                cell_map, trace_map, worker_binary_sha256)
            futures[future] = invocation.ordinal
        remaining = deadline - time.monotonic()
        if remaining <= 0.0:
            registry.cancel()
            fail("work/rank campaign hard wall expired")
        for future in concurrent.futures.as_completed(
                futures, timeout=remaining):
            result = future.result()
            if result.invocation.ordinal in results:
                fail("duplicate invocation result")
            results[result.invocation.ordinal] = result
    except concurrent.futures.TimeoutError:
        registry.cancel()
        fail("work/rank campaign hard wall expired")
    except BaseException:
        registry.cancel()
        raise
    finally:
        for future in futures:
            future.cancel()
        if executor is not None:
            executor.shutdown(wait=True)
    if len(results) != EXPECTED_INVOCATIONS:
        fail("work/rank campaign did not complete all invocations")
    if _source_provenance(worker) != source_provenance:
        fail("worker, audited source, or tracked git state changed during run")

    rows: List[Mapping[str, Any]] = []
    receipts: List[Mapping[str, Any]] = []
    for ordinal in range(EXPECTED_INVOCATIONS):
        result = results[ordinal]
        rows.extend(result.rows)
        receipts.append({
            "ordinal": ordinal,
            "arm": result.invocation.arm,
            "trial": result.invocation.trial,
            "schedule": result.invocation.schedule,
            "row_count": len(result.rows),
            "command_sha256": result.command_sha256,
            "worker_stdout_sha256": result.stdout_sha256,
        })
    expected_keys = set()
    for arm, _, _ in ARMS:
        for trial in range(3):
            for schedule, _ in STRATA:
                for K in K_VALUES:
                    for overhead in OVERHEADS:
                        expected_keys.add((arm, trial, schedule, K, overhead))
    actual_keys = {
        (row["arm"], row["trial"], row["schedule"], row["K"],
         row["overhead"]) for row in rows
    }
    if (len(rows) != EXPECTED_RECORDS or len(actual_keys) != EXPECTED_RECORDS or
            actual_keys != expected_keys):
        fail("work/rank campaign is not exactly 5760 unique frozen rows")
    records: List[Mapping[str, Any]] = []
    for ordinal, row in enumerate(rows):
        record = dict(row)
        record["ordinal"] = ordinal
        record["source_provenance_sha256"] = source_provenance_sha256
        records.append(record)
    result_bytes = b"".join(
        (_canonical(record) + "\n").encode("utf-8") for record in records)
    result_sha256 = hashlib.sha256(result_bytes).hexdigest()
    unsigned_summary = {
        "schema": SUMMARY_SCHEMA,
        "contract_sha256": contract_api.contract_sha256(contract),
        "recovery_domain_sha256": contract_api.recovery_domain_sha256(
            contract, "development"),
        "work_domain_sha256": _domain_sha256(),
        "worker_binary_sha256": worker_binary_sha256,
        "source_provenance": source_provenance,
        "source_provenance_sha256": source_provenance_sha256,
        "result_stream_sha256": result_sha256,
        "record_count": len(records),
        "invocation_count": len(receipts),
        "invocations": receipts,
        "trace_identity": trace_identity,
        "arms": _aggregate(records),
        "seed_basis": RAW_SEED_BASIS,
        "seed_schedule_sha256": RAW_SEED_SCHEDULE_SHA256,
        "selection_performed": False,
        "timing_metrics_included": False,
    }
    summary = dict(unsigned_summary)
    summary["summary_sha256"] = _sha256_json(unsigned_summary)
    summary_bytes = (_canonical(summary) + "\n").encode("utf-8")
    _atomic_write(result_path, result_bytes)
    _atomic_write(summary_path, summary_bytes)
    return summary


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=min(12, os.cpu_count() or 1))
    parser.add_argument("--deadline-seconds", type=float,
                        default=MAX_WALL_SECONDS)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    try:
        summary = run_campaign(
            args.worker, args.output_dir, args.jobs, args.deadline_seconds)
    except (WorkScreenError, contract_api.ContractError, OSError) as exc:
        print("wh2 precodefail work screen: {}".format(exc), file=sys.stderr)
        return 1
    print(_canonical(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
