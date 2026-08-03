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
import ctypes
import csv
from decimal import Decimal, InvalidOperation
import hashlib
import math
import os
from pathlib import Path
import re
import selectors
import secrets
import signal
import stat
import subprocess
import sys
import threading
import time
from dataclasses import dataclass
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_benchmark_contract as contract_api


RECORD_SCHEMA = "wirehair.wh2.precodefail-work-record.v3"
SUMMARY_SCHEMA = "wirehair.wh2.precodefail-work-summary.v3"
TRACE_IDENTITY_SCHEMA = "wirehair.wh2.precodefail-trace-identity.v1"
NATIVE_DESCRIPTOR_SCHEMA = "wirehair.wh2.native-arm-descriptor.v1"
RESULT_NAME = "work-rank-results.jsonl"
SUMMARY_NAME = "work-rank-summary.json"
MAX_WALL_SECONDS = 7200.0
MAX_JOBS = 48
MAX_STDOUT_BYTES = 2 * 1024 * 1024
MAX_STDERR_BYTES = 64 * 1024
EXPECTED_RECORDS = 2880
EXPECTED_INVOCATIONS = 24
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
_LIBC = ctypes.CDLL(None, use_errno=True)
_RENAMEAT2 = getattr(_LIBC, "renameat2", None)
if _RENAMEAT2 is not None:
    _RENAMEAT2.argtypes = (
        ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p,
        ctypes.c_uint)
    _RENAMEAT2.restype = ctypes.c_int
_RENAME_NOREPLACE = 1

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
    ("wirehair2_raw_d12_h12_periodic", 12, 12, "disabled"),
    ("wirehair2_dense_two07_basis_v1", 12, 12, "two07"),
)
ARM_TRANSFORMS = {
    "wirehair2_raw_d12_h12_periodic": "d12-h12-periodic",
    "wirehair2_dense_two07_basis_v1":
        "dense-anchor-two07+basis-segment-adjacent-symdiff-v1",
}
# These are SHA-256 hashes of the closed, structure-only native arm
# descriptors.  The uniform raw seed schedule is bound separately by every
# realized-construction receipt, so changing the schedule cannot masquerade
# as changing the equation structure (or vice versa).
ARM_DESCRIPTOR_SHA256 = {
    "wirehair2_raw_d12_h12_periodic":
        "739092a7824449e6168f08b46661dfbe8ad5495ea4166b36073c79cd3bacdd11",
    "wirehair2_dense_two07_basis_v1":
        "9527f200ad38c7eec6502b2f768fdd67b92787fb227eed3d7616274ffc2df388",
}

CSV_HEADER = (
    "N", "bb", "heavy_family", "mix_count", "staircase",
    "binary_dense_rows", "gf256_heavy_rows", "source_hits",
    "dense_identity_corner", "overhead", "trials",
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
    "gf256_heavy_rows_override", "dense_anchor_layout",
    "odd_packet_peel_seed_xor",
    "packet_row_seed_multiplier", "packet_row_seed_avalanche",
    "seed_block_bytes_override", "overhead_stream", "full_payload_solve",
    "schedule", "exact_attempt_mode", "exact_precode_attempt",
    "exact_packet_attempt", "construction_seed_basis",
    "seed_schedule_sha256", "source_git_commit",
))
SUMMARY_FIELDS = frozenset((
    "schema", "contract_sha256", "recovery_domain_sha256",
    "work_domain_sha256", "worker_binary_sha256", "source_provenance",
    "source_provenance_sha256", "result_stream_sha256", "record_count",
    "invocation_count", "invocations", "trace_identity", "arms",
    "construction_seed_basis", "seed_schedule_sha256",
    "selection_performed", "timing_metrics_included", "summary_sha256",
))
SOURCE_PROVENANCE_FIELDS = frozenset((
    "source_git_commit", "tracked_worktree_clean", "tracked_status_sha256",
    "untracked_source_clean", "untracked_source_sha256",
    "worker_source_git_commit",
    "worker_binary_sha256", "controller_source_sha256",
    "wirehair_v2_bench_source_sha256",
    "wirehair_v2_raw_seed_source_sha256",
    "wh2_frozen_trace_source_sha256",
))
INVOCATION_RECEIPT_FIELDS = frozenset((
    "ordinal", "arm", "trial", "schedule", "row_count", "command_sha256",
    "worker_stdout_sha256",
))
RECORD_FIELDS = frozenset((
    "schema", "arm", "arm_descriptor_sha256",
    "realized_construction_sha256", "construction_seed_basis",
    "seed_schedule_sha256", "precode_attempt", "packet_attempt",
    "effective_precode_seed", "effective_packet_seed", "staircase",
    "binary_dense_rows", "gf256_heavy_rows", "source_hits",
    "dense_anchor_layout", "dense_identity_corner", "heavy_family",
    "mix_count", "phase", "band",
    "K", "block_bytes", "loss_ppm", "schedule", "trial",
    "base_seed_attempt", "loss_seed", "overhead", "attempt_mode",
    "cell_sha256", "frozen_trace_sha256", "frozen_packet_prefix_sha256",
    "packet_count", "frozen_trace_attempted_candidates",
    "frozen_prefix_attempted_candidates", "success", "rank_fail", "error",
    "inactivated_columns", "binary_deficiency", "gf256_rank_gain",
    "heavy_shortfall", "block_xors", "gf256_muladds", "command_sha256",
    "worker_stdout_sha256", "worker_binary_sha256", "ordinal",
    "source_provenance_sha256",
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
    dense_anchor_layout: str
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
            "--dense-anchors", self.dense_anchor_layout,
            "--paired-overhead-stream", "--full-payload-solve",
            "--exact-precode-attempt", str(self.trial),
            "--exact-packet-attempt", str(self.trial),
            "--construction-seed-basis", RAW_SEED_BASIS,
        ]

    def identity(self) -> Mapping[str, Any]:
        return {
            "arm": self.arm,
            "binary_dense_rows": self.dense_rows,
            "dense_anchor_layout": self.dense_anchor_layout,
            "gf256_heavy_rows": self.heavy_rows,
            "heavy_family": "periodic",
            "loss_ppm": self.loss_ppm,
            "loss_seed": self.root,
            "mix_count": 3,
            "attempt_mode": "exact",
            "packet_attempt": self.trial,
            "precode_attempt": self.trial,
            "schedule": self.schedule,
            "construction_seed_basis": RAW_SEED_BASIS,
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


def _sha256_fd(descriptor: int, context: str) -> str:
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            fail("{} must be a regular file".format(context))
        digest = hashlib.sha256()
        offset = 0
        while offset < before.st_size:
            block = os.pread(
                descriptor, min(1024 * 1024, before.st_size - offset),
                offset)
            if not block:
                fail("{} was truncated while hashing".format(context))
            digest.update(block)
            offset += len(block)
        if os.pread(descriptor, 1, offset):
            fail("{} grew while hashing".format(context))
        after = os.fstat(descriptor)
    except OSError as exc:
        fail("cannot hash {}: {}".format(context, exc))
    identity_before = (
        before.st_dev, before.st_ino, before.st_size,
        before.st_mtime_ns, before.st_ctime_ns)
    identity_after = (
        after.st_dev, after.st_ino, after.st_size,
        after.st_mtime_ns, after.st_ctime_ns)
    if identity_before != identity_after:
        fail("{} changed while hashing".format(context))
    return digest.hexdigest()


def _open_pinned_worker(path: Path) -> Tuple[Path, int]:
    try:
        resolved = path.resolve(strict=True)
        flags = (os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0) |
                 getattr(os, "O_CLOEXEC", 0) |
                 getattr(os, "O_NONBLOCK", 0))
        descriptor = os.open(str(resolved), flags)
        opened = os.fstat(descriptor)
        current = os.stat(str(resolved), follow_symlinks=False)
        if (not stat.S_ISREG(opened.st_mode) or
                not opened.st_mode & 0o111 or
                (opened.st_dev, opened.st_ino) !=
                (current.st_dev, current.st_ino)):
            fail("worker must be one stable executable regular file")
        _sha256_fd(descriptor, "precodefail worker")
        return resolved, descriptor
    except BaseException:
        if "descriptor" in locals():
            try:
                os.close(descriptor)
            except BaseException:
                pass
        raise


def _effective_precode_seed(attempt: int) -> str:
    value = (RAW_PRECODE_BASE_SEED +
             attempt * RAW_PRECODE_ATTEMPT_STRIDE) & MASK64
    return "0x{:016x}".format(value)


def _effective_packet_seed(attempt: int) -> str:
    value = (RAW_PACKET_BASE_SEED +
             attempt * RAW_PACKET_ATTEMPT_STRIDE) & MASK32
    return "0x{:08x}".format(value)


def _arm_descriptor_sha256(arm: str) -> str:
    transform = ARM_TRANSFORMS.get(arm)
    expected = ARM_DESCRIPTOR_SHA256.get(arm)
    if transform is None or expected is None or \
            set(ARM_TRANSFORMS) != set(ARM_DESCRIPTOR_SHA256):
        fail("raw arm is not in the closed structure descriptor roster")
    descriptor = {
        "arm": arm,
        "codec": "wirehair2_experiment",
        "equation_transform": transform,
        "schema": NATIVE_DESCRIPTOR_SCHEMA,
    }
    actual = _sha256_json(descriptor)
    if actual != expected:
        fail("raw structure descriptor golden is inconsistent")
    return expected


def _raw_realized_construction_sha256(
        arm: str, K: int, block_bytes: int,
        raw_fields: Mapping[str, Any]) -> str:
    try:
        return contract_api.raw_realized_construction_sha256(
            "wirehair2_experiment", arm, _arm_descriptor_sha256(arm), K,
            block_bytes, raw_fields)
    except contract_api.ContractError as exc:
        fail(str(exc))
    return ""


def _git_source_state(repo_root: Path = REPO_ROOT) -> Mapping[str, Any]:
    """Bind one stable HEAD and reject relevant untracked build inputs."""
    commands = (
        ["git", "rev-parse", "--verify", "HEAD^{commit}"],
        ["git", "status", "--porcelain=v1", "--untracked-files=no"],
        ["git", "ls-files", "--others", "--exclude-standard", "--",
         "CMakeLists.txt", "codec", "bench", "include", "cmake",
         ":(top,glob)*.c", ":(top,glob)*.cc", ":(top,glob)*.cpp",
         ":(top,glob)*.h", ":(top,glob)*.hpp", ":(top,glob)*.inc"],
        ["git", "rev-parse", "--verify", "HEAD^{commit}"],
    )
    outputs = []
    try:
        for command in commands:
            result = subprocess.run(
                command, cwd=str(repo_root), stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=10.0,
                check=False)
            if result.returncode != 0:
                fail("cannot bind source git provenance")
            outputs.append(result.stdout)
    except (OSError, subprocess.TimeoutExpired) as exc:
        fail("cannot bind source git provenance: {}".format(exc))
    try:
        first = outputs[0].decode("ascii").strip()
        status = outputs[1].decode("utf-8")
        untracked = outputs[2].decode("utf-8")
        second = outputs[3].decode("ascii").strip()
    except UnicodeDecodeError:
        fail("source git provenance is not valid text")
    if re.fullmatch(r"[0-9a-f]{40}", first) is None or first != second:
        fail("source git commit is malformed or changed during inspection")
    return {
        "source_git_commit": first,
        "tracked_worktree_clean": not bool(status),
        "tracked_status_sha256":
            hashlib.sha256(status.encode("utf-8")).hexdigest(),
        "untracked_source_clean": not bool(untracked),
        "untracked_source_sha256":
            hashlib.sha256(untracked.encode("utf-8")).hexdigest(),
    }


def _source_provenance(
        worker: Path, worker_fd: Optional[int] = None) -> Mapping[str, Any]:
    for path in (
            BENCH_SOURCE, RAW_SEED_SOURCE, FROZEN_TRACE_SOURCE,
            CONTROLLER_SOURCE):
        try:
            info = path.stat()
        except OSError as exc:
            fail("cannot inspect audited source {}: {}".format(path, exc))
        if not stat.S_ISREG(info.st_mode):
            fail("audited source must be a regular file: {}".format(path))
    git_state = _git_source_state()
    if worker_fd is None:
        worker_sha256 = _sha256_file(worker)
    else:
        try:
            opened = os.fstat(worker_fd)
            current = os.stat(str(worker), follow_symlinks=False)
        except OSError as exc:
            fail("cannot revalidate pinned worker: {}".format(exc))
        if ((opened.st_dev, opened.st_ino) !=
                (current.st_dev, current.st_ino) or
                not stat.S_ISREG(opened.st_mode)):
            fail("precodefail worker path no longer names the pinned binary")
        worker_sha256 = _sha256_fd(worker_fd, "precodefail worker")
    return {
        **git_state,
        # Every invocation must independently emit this same embedded commit
        # before the controller can publish the campaign.
        "worker_source_git_commit": git_state["source_git_commit"],
        "worker_binary_sha256": worker_sha256,
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
    for arm, dense_rows, heavy_rows, dense_anchor_layout in ARMS:
        for trial, root in enumerate(ROOTS):
            for schedule, loss_ppm in STRATA:
                values.append(Invocation(
                    len(values), arm, dense_rows, heavy_rows,
                    dense_anchor_layout, trial, root, schedule, loss_ppm))
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


def _parse_metadata(
        line: str, invocation: Invocation,
        expected_source_git_commit: str) -> Mapping[str, str]:
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
        "dense_anchor_layout": invocation.dense_anchor_layout,
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
        "source_git_commit": expected_source_git_commit,
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
        worker_binary_sha256: str,
        expected_source_git_commit: str) -> InvocationResult:
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
    _parse_metadata(lines[0], invocation, expected_source_git_commit)
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
            "binary_dense_rows": str(invocation.dense_rows),
            "gf256_heavy_rows": str(invocation.heavy_rows),
            "dense_identity_corner": "0",
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

        staircase = _parse_uint(row["staircase"], "staircase", 400)
        binary_dense_rows = _parse_uint(
            row["binary_dense_rows"], "binary_dense_rows", 64)
        gf256_heavy_rows = _parse_uint(
            row["gf256_heavy_rows"], "gf256_heavy_rows", 128)
        source_hits = _parse_uint(row["source_hits"], "source_hits", 8)
        dense_identity_corner = bool(_parse_uint(
            row["dense_identity_corner"], "dense_identity_corner", 1))
        expected_source_hits = 3 if K >= 10000 else 2
        if (staircase == 0 or binary_dense_rows != invocation.dense_rows or
                gf256_heavy_rows != invocation.heavy_rows or
                source_hits != expected_source_hits or dense_identity_corner or
                K + staircase + binary_dense_rows + gf256_heavy_rows > 65535):
            fail("precodefail row has an invalid realized raw structure")

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
        raw_fields = {
            "construction_seed_basis": RAW_SEED_BASIS,
            "seed_schedule_sha256": RAW_SEED_SCHEDULE_SHA256,
            "precode_attempt": invocation.trial,
            "packet_attempt": invocation.trial,
            "effective_precode_seed": expected_precode_seed,
            "effective_packet_seed": expected_packet_seed,
            "staircase": staircase,
            "binary_dense_rows": binary_dense_rows,
            "gf256_heavy_rows": gf256_heavy_rows,
            "source_hits": source_hits,
            "dense_identity_corner": dense_identity_corner,
            "dense_anchor_layout": invocation.dense_anchor_layout,
            "heavy_family": "periodic-cauchy",
            "mix_count": 3,
        }
        arm_descriptor_sha256 = _arm_descriptor_sha256(invocation.arm)
        parsed.append({
            "schema": RECORD_SCHEMA,
            "arm": invocation.arm,
            "arm_descriptor_sha256": arm_descriptor_sha256,
            "realized_construction_sha256":
                _raw_realized_construction_sha256(
                    invocation.arm, K, 2, raw_fields),
            **raw_fields,
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
            "attempt_mode": "exact",
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
        worker_binary_sha256: str,
        expected_source_git_commit: str,
        worker_fd: Optional[int] = None) -> InvocationResult:
    if registry.cancelled.is_set():
        fail("work/rank campaign was cancelled")
    command = invocation.argv(worker)
    execution_command = list(command)
    pass_fds: Tuple[int, ...] = ()
    if worker_fd is not None:
        execution_command[0] = "/proc/self/fd/{}".format(worker_fd)
        pass_fds = (worker_fd,)
    try:
        process = subprocess.Popen(
            execution_command, stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            start_new_session=True, pass_fds=pass_fds)
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
            worker_binary_sha256, expected_source_git_commit)
    except BaseException:
        registry.cancel()
        _kill_and_reap(process)
        raise
    finally:
        registry.remove(process)


def _domain_sha256(records: Sequence[Mapping[str, Any]]) -> str:
    if len(records) != EXPECTED_RECORDS:
        fail("work domain requires exactly {} records".format(
            EXPECTED_RECORDS))
    digest = hashlib.sha256()
    identity_fields = (
        "K", "arm", "arm_descriptor_sha256", "binary_dense_rows",
        "block_bytes", "construction_seed_basis", "dense_anchor_layout",
        "dense_identity_corner",
        "effective_packet_seed", "effective_precode_seed",
        "gf256_heavy_rows", "heavy_family", "loss_ppm", "loss_seed",
        "mix_count", "overhead", "packet_attempt", "precode_attempt",
        "realized_construction_sha256", "schedule",
        "seed_schedule_sha256", "source_hits", "staircase", "trial",
    )
    ordered = sorted(records, key=lambda row: (
        row["arm"], row["trial"], row["schedule"], row["K"],
        row["overhead"]))
    for row in ordered:
        if not all(field in row for field in identity_fields):
            fail("work domain record is missing raw construction identity")
        identity = {field: row[field] for field in identity_fields}
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
    for arm, dense_rows, heavy_rows, dense_anchor_layout in ARMS:
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
            "arm_descriptor_sha256": _arm_descriptor_sha256(arm),
            "binary_dense_rows": dense_rows,
            "dense_anchor_layout": dense_anchor_layout,
            "gf256_heavy_rows": heavy_rows,
            "heavy_family": "periodic-cauchy",
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


def _close_preserving_exception(descriptor: int) -> None:
    initiating_exception = sys.exc_info()[0] is not None
    try:
        os.close(descriptor)
    except BaseException:
        if not initiating_exception:
            raise


def _open_pinned_directory(path: Path) -> Tuple[int, Tuple[int, int]]:
    flags = (os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) |
             getattr(os, "O_NOFOLLOW", 0))
    try:
        directory_fd = os.open(str(path), flags)
        opened = os.fstat(directory_fd)
        current = os.stat(str(path), follow_symlinks=False)
        identity = (opened.st_dev, opened.st_ino)
        if not stat.S_ISDIR(opened.st_mode) or identity != (
                current.st_dev, current.st_ino):
            fail("artifact output directory identity changed while opening")
        return directory_fd, identity
    except BaseException:
        if "directory_fd" in locals():
            _close_preserving_exception(directory_fd)
        raise


def _verify_pinned_directory(
        path: Path, identity: Tuple[int, int]) -> None:
    try:
        current = os.stat(str(path), follow_symlinks=False)
    except OSError as exc:
        fail("artifact output directory became unavailable: {}".format(exc))
    if not stat.S_ISDIR(current.st_mode) or (
            current.st_dev, current.st_ino) != identity:
        fail("artifact output directory identity changed during publication")


def _remove_if_identity_at(
        directory_fd: int, name: str, identity: Tuple[int, int],
        display_path: Path) -> None:
    """Rollback only the artifact inode published by this process."""
    try:
        if _quarantine_unlink_if_identity_at(
                directory_fd, name, identity, display_path):
            os.fsync(directory_fd)
    except OSError as exc:
        fail("cannot roll back partial artifact publication {}: {}".format(
            display_path, exc))


def _rename_noreplace_at(
        directory_fd: int, source: str, destination: str) -> None:
    if _RENAMEAT2 is None:
        fail("artifact publication requires Linux renameat2")
    result = _RENAMEAT2(
        directory_fd, source.encode("utf-8"),
        directory_fd, destination.encode("utf-8"), _RENAME_NOREPLACE)
    if result != 0:
        error_number = ctypes.get_errno()
        raise OSError(error_number, os.strerror(error_number))


def _probe_rename_noreplace(directory_fd: int) -> None:
    source = ".renameat2-probe.{}.src".format(secrets.token_hex(16))
    destination = ".renameat2-probe.{}.dst".format(secrets.token_hex(16))
    descriptor = -1
    try:
        flags = (os.O_WRONLY | os.O_CREAT | os.O_EXCL |
                 getattr(os, "O_NOFOLLOW", 0) |
                 getattr(os, "O_CLOEXEC", 0))
        descriptor = os.open(source, flags, 0o600, dir_fd=directory_fd)
        closing_descriptor = descriptor
        descriptor = -1
        os.close(closing_descriptor)
        _rename_noreplace_at(directory_fd, source, destination)
        moved = os.stat(
            destination, dir_fd=directory_fd, follow_symlinks=False)
        if not stat.S_ISREG(moved.st_mode):
            fail("renameat2 capability probe produced a non-regular file")
    except OSError as exc:
        fail("output directory does not support renameat2 NOREPLACE: {}".format(
            exc))
    finally:
        if descriptor >= 0:
            try:
                os.close(descriptor)
            except BaseException:
                pass
        for name in (source, destination):
            try:
                os.unlink(name, dir_fd=directory_fd)
            except FileNotFoundError:
                pass
        os.fsync(directory_fd)


def _quarantine_unlink_if_identity_at(
        directory_fd: int, name: str, identity: Tuple[int, int],
        display_path: Path) -> bool:
    """Remove one matching name without ever unlinking a path replacement."""
    quarantine = ".{}.rollback.{}.tmp".format(name, secrets.token_hex(16))
    move_attempted = False
    try:
        move_attempted = True
        _rename_noreplace_at(directory_fd, name, quarantine)
    except FileNotFoundError:
        return False
    except BaseException:
        if move_attempted:
            try:
                _rename_noreplace_at(directory_fd, quarantine, name)
            except BaseException:
                pass
        raise
    try:
        moved = os.stat(
            quarantine, dir_fd=directory_fd, follow_symlinks=False)
        if (moved.st_dev, moved.st_ino) == identity:
            os.unlink(quarantine, dir_fd=directory_fd)
            return True

        # A foreign replacement won the race.  Restore it without overwriting
        # a second replacement.  If the original name is occupied, retain both
        # foreign inodes (one under the unpredictable quarantine name).
        try:
            _rename_noreplace_at(directory_fd, quarantine, name)
        except FileExistsError:
            pass
        return False
    except BaseException:
        # Put the atomically quarantined name back whenever possible so the
        # outer transaction can retry identity-checked cleanup.  This also
        # covers an asynchronous exception immediately after renameat2.
        if move_attempted:
            try:
                _rename_noreplace_at(directory_fd, quarantine, name)
            except (FileNotFoundError, FileExistsError):
                pass
            except BaseException:
                pass
        raise


@dataclass(frozen=True)
class _PublishedArtifact:
    device: int
    inode: int
    descriptor: int
    staged_name: str
    byte_count: int
    sha256: str

    @property
    def identity(self) -> Tuple[int, int]:
        return self.device, self.inode


def _verify_identity_at(
        directory_fd: int, name: str, expected: _PublishedArtifact,
        display_path: Path, require_staged: bool = True) -> None:
    try:
        retained = os.fstat(expected.descriptor)
        current = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
        if require_staged:
            staged = os.stat(
                expected.staged_name, dir_fd=directory_fd,
                follow_symlinks=False)
    except OSError as exc:
        fail("published artifact became unavailable {}: {}".format(
            display_path, exc))
    if (not stat.S_ISREG(retained.st_mode) or
            (retained.st_dev, retained.st_ino) != expected.identity or
            (current.st_dev, current.st_ino) != expected.identity or
            (require_staged and
             (staged.st_dev, staged.st_ino) != expected.identity)):
        fail("published artifact identity changed for {}".format(
            display_path))
    flags = (os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0) |
             getattr(os, "O_NONBLOCK", 0) |
             getattr(os, "O_CLOEXEC", 0))
    try:
        descriptor = os.open(name, flags, dir_fd=directory_fd)
    except OSError as exc:
        fail("cannot open published artifact {}: {}".format(
            display_path, exc))
    try:
        before = os.fstat(descriptor)
        if (not stat.S_ISREG(before.st_mode) or
                (before.st_dev, before.st_ino) != expected.identity):
            fail("published artifact identity changed while opening {}".format(
                display_path))
        digest = hashlib.sha256()
        byte_count = 0
        while True:
            remaining = expected.byte_count - byte_count
            if remaining < 0:
                fail("published artifact exceeded its fixed size {}".format(
                    display_path))
            chunk = os.read(descriptor, min(1024 * 1024, remaining + 1))
            if not chunk:
                break
            digest.update(chunk)
            byte_count += len(chunk)
        after = os.fstat(descriptor)
    except OSError as exc:
        fail("cannot verify published artifact {}: {}".format(
            display_path, exc))
    finally:
        _close_preserving_exception(descriptor)
    stable_fields_before = (
        before.st_dev, before.st_ino, before.st_size,
        before.st_mtime_ns, before.st_ctime_ns)
    stable_fields_after = (
        after.st_dev, after.st_ino, after.st_size,
        after.st_mtime_ns, after.st_ctime_ns)
    if stable_fields_before != stable_fields_after:
        fail("published artifact changed while verifying {}".format(
            display_path))
    if (byte_count != expected.byte_count or
            digest.hexdigest() != expected.sha256):
        fail("published artifact content changed for {}".format(display_path))


def _discard_staged_at(
        directory_fd: int, expected: _PublishedArtifact,
        display_path: Path) -> None:
    if not _quarantine_unlink_if_identity_at(
            directory_fd, expected.staged_name, expected.identity,
            display_path):
        fail("staged artifact disappeared for {}".format(display_path))


def _rollback_publications_at(
        directory_fd: int, paths: Sequence[Path],
        expected: Mapping[str, _PublishedArtifact]) -> None:
    """Best-effort rollback that never masks the initiating exception."""
    for path in paths:
        publication = expected.get(path.name)
        if publication is None:
            continue
        try:
            _remove_if_identity_at(
                directory_fd, path.name, publication.identity, path)
        except BaseException:
            pass
    for path in paths:
        publication = expected.get(path.name)
        if publication is None:
            continue
        try:
            _quarantine_unlink_if_identity_at(
                directory_fd, publication.staged_name,
                publication.identity, path)
        except BaseException:
            pass
    try:
        os.fsync(directory_fd)
    except BaseException:
        pass


def _close_publication_descriptors(
        expected: Mapping[str, _PublishedArtifact]) -> None:
    for publication in expected.values():
        try:
            os.close(publication.descriptor)
        except BaseException:
            # Descriptor closure occurs only after either successful commit or
            # completed best-effort rollback.  It cannot change that outcome.
            pass


def _atomic_write_at(
        directory_fd: int, name: str, display_path: Path, data: bytes,
        expected_identities: Dict[str, _PublishedArtifact]
        ) -> Tuple[int, int]:
    if not name or name in (".", "..") or "/" in name:
        fail("artifact destination must be a plain filename")
    staged_name = ".{}.{}.{}.{}.tmp".format(
        name, os.getpid(), threading.get_ident(), secrets.token_hex(16))
    flags = (os.O_WRONLY | os.O_CREAT | os.O_EXCL |
             getattr(os, "O_NOFOLLOW", 0) |
             getattr(os, "O_CLOEXEC", 0))
    descriptor = -1
    retained_descriptor = -1
    staged_identity: Optional[Tuple[int, int]] = None
    caller_owns_staged = False
    try:
        descriptor = os.open(
            staged_name, flags, 0o600, dir_fd=directory_fd)
        opened_info = os.fstat(descriptor)
        staged_identity = (opened_info.st_dev, opened_info.st_ino)
        remaining = memoryview(data)
        while remaining:
            written = os.write(descriptor, remaining)
            if written <= 0:
                fail("short write while staging artifact {}".format(
                    display_path))
            remaining = remaining[written:]
        os.fsync(descriptor)
        staged_info = os.fstat(descriptor)
        if (staged_info.st_dev, staged_info.st_ino) != staged_identity:
            fail("staged artifact identity changed while writing {}".format(
                display_path))

        # Publish the expected inode before linking.  The pair-level rollback
        # can therefore identify our destination even if an asynchronous
        # BaseException lands after link(2) but before this function returns.
        retained_descriptor = os.dup(descriptor)
        expected_identities[name] = _PublishedArtifact(
            device=staged_info.st_dev,
            inode=staged_info.st_ino,
            descriptor=retained_descriptor,
            staged_name=staged_name,
            byte_count=len(data),
            sha256=hashlib.sha256(data).hexdigest())
        caller_owns_staged = True
        try:
            os.link(
                staged_name, name,
                src_dir_fd=directory_fd, dst_dir_fd=directory_fd,
                follow_symlinks=False)
        except FileExistsError:
            fail("refusing to replace existing artifact {}".format(
                display_path))
        destination_info = os.stat(
            name, dir_fd=directory_fd, follow_symlinks=False)
        if (destination_info.st_dev, destination_info.st_ino) != (
                staged_identity):
            fail("published artifact identity mismatch for {}".format(
                display_path))
        os.fsync(directory_fd)
        return staged_identity
    except BaseException:
        if staged_identity is not None:
            try:
                _remove_if_identity_at(
                    directory_fd, name, staged_identity, display_path)
            except BaseException:
                pass
        raise
    finally:
        initiating_exception = sys.exc_info()[0] is not None
        cleanup_error: Optional[BaseException] = None
        if descriptor >= 0:
            if staged_identity is None:
                try:
                    fallback_info = os.stat(descriptor)
                    staged_identity = (
                        fallback_info.st_dev, fallback_info.st_ino)
                except BaseException:
                    pass
            try:
                os.close(descriptor)
            except BaseException as exc:
                cleanup_error = exc
        if not caller_owns_staged and staged_identity is not None:
            try:
                _quarantine_unlink_if_identity_at(
                    directory_fd, staged_name, staged_identity, display_path)
            except BaseException as exc:
                if cleanup_error is None:
                    cleanup_error = exc
        if cleanup_error is not None and not initiating_exception:
            raise cleanup_error
        if name not in expected_identities and retained_descriptor >= 0:
            try:
                os.close(retained_descriptor)
            except BaseException:
                pass


def _atomic_write(path: Path, data: bytes) -> Tuple[int, int]:
    directory_fd, directory_identity = _open_pinned_directory(path.parent)
    expected: Dict[str, _PublishedArtifact] = {}
    try:
        _probe_rename_noreplace(directory_fd)
        identity = _atomic_write_at(
            directory_fd, path.name, path, data, expected)
        _verify_identity_at(
            directory_fd, path.name, expected[path.name], path)
        _verify_pinned_directory(path.parent, directory_identity)
        _discard_staged_at(directory_fd, expected[path.name], path)
        os.fsync(directory_fd)
        _verify_identity_at(
            directory_fd, path.name, expected[path.name], path,
            require_staged=False)
        _verify_pinned_directory(path.parent, directory_identity)
        return identity
    except BaseException:
        _rollback_publications_at(directory_fd, (path,), expected)
        raise
    finally:
        try:
            _close_publication_descriptors(expected)
        finally:
            try:
                os.close(directory_fd)
            except BaseException:
                pass


def _publish_artifact_pair(
        result_path: Path, result_bytes: bytes,
        summary_path: Path, summary_bytes: bytes) -> None:
    if result_path.parent != summary_path.parent:
        fail("artifact pair must share one output directory")
    if result_path.name == summary_path.name:
        fail("artifact pair destinations must be distinct")
    directory_fd, directory_identity = _open_pinned_directory(
        result_path.parent)
    expected: Dict[str, _PublishedArtifact] = {}
    try:
        _probe_rename_noreplace(directory_fd)
        _atomic_write_at(
            directory_fd, result_path.name, result_path, result_bytes,
            expected)
        _atomic_write_at(
            directory_fd, summary_path.name, summary_path, summary_bytes,
            expected)
        for path in (result_path, summary_path):
            _verify_identity_at(
                directory_fd, path.name, expected[path.name], path)
        _verify_pinned_directory(result_path.parent, directory_identity)
        for path in (summary_path, result_path):
            _discard_staged_at(directory_fd, expected[path.name], path)
        os.fsync(directory_fd)
        for path in (result_path, summary_path):
            _verify_identity_at(
                directory_fd, path.name, expected[path.name], path,
                require_staged=False)
        _verify_pinned_directory(result_path.parent, directory_identity)
        return
    except BaseException:
        _rollback_publications_at(
            directory_fd, (summary_path, result_path), expected)
        raise
    finally:
        try:
            _close_publication_descriptors(expected)
        finally:
            try:
                os.close(directory_fd)
            except BaseException:
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
    worker, worker_fd = _open_pinned_worker(worker)
    try:
        return _run_campaign_pinned(
            worker, worker_fd, output_dir, jobs, deadline_seconds)
    finally:
        try:
            os.close(worker_fd)
        except BaseException:
            pass


def _run_campaign_pinned(
        worker: Path, worker_fd: int, output_dir: Path, jobs: int,
        deadline_seconds: float) -> Mapping[str, Any]:
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
    source_provenance = _source_provenance(worker, worker_fd)
    if (source_provenance["tracked_worktree_clean"] is not True or
            source_provenance["untracked_source_clean"] is not True):
        fail("work/rank campaign requires a clean codec source worktree")
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
                cell_map, trace_map, worker_binary_sha256,
                source_provenance["source_git_commit"], worker_fd)
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
    if _source_provenance(worker, worker_fd) != source_provenance:
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
    for arm, _, _, _ in ARMS:
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
        fail("work/rank campaign is not exactly 2880 unique frozen rows")
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
        "work_domain_sha256": _domain_sha256(records),
        "worker_binary_sha256": worker_binary_sha256,
        "source_provenance": source_provenance,
        "source_provenance_sha256": source_provenance_sha256,
        "result_stream_sha256": result_sha256,
        "record_count": len(records),
        "invocation_count": len(receipts),
        "invocations": receipts,
        "trace_identity": trace_identity,
        "arms": _aggregate(records),
        "construction_seed_basis": RAW_SEED_BASIS,
        "seed_schedule_sha256": RAW_SEED_SCHEDULE_SHA256,
        "selection_performed": False,
        "timing_metrics_included": False,
    }
    summary = dict(unsigned_summary)
    summary["summary_sha256"] = _sha256_json(unsigned_summary)
    summary_bytes = (_canonical(summary) + "\n").encode("utf-8")
    _publish_artifact_pair(
        result_path, result_bytes, summary_path, summary_bytes)
    return summary


def _is_sha256(value: Any) -> bool:
    return isinstance(value, str) and re.fullmatch(
        r"[0-9a-f]{64}", value) is not None


def _read_regular_bytes(
        path: Path, byte_cap: int, context: str,
        directory_fd: Optional[int] = None) -> bytes:
    descriptor = -1
    try:
        nofollow = getattr(os, "O_NOFOLLOW", 0)
        if nofollow == 0:
            fail("{} cannot be opened fail-closed without O_NOFOLLOW".format(
                context))
        descriptor = os.open(
            str(path), os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
            getattr(os, "O_NONBLOCK", 0) | nofollow, dir_fd=directory_fd)
        before = os.fstat(descriptor)
        if (not stat.S_ISREG(before.st_mode) or before.st_size <= 0 or
                before.st_size > byte_cap):
            fail("{} has an invalid type or size".format(context))
        blocks = []
        size = 0
        while True:
            block = os.read(descriptor, min(1024 * 1024, byte_cap + 1 - size))
            if not block:
                break
            blocks.append(block)
            size += len(block)
            if size > byte_cap:
                fail("{} exceeds its bounded size".format(context))
        after = os.fstat(descriptor)
        stable_before = (
            before.st_dev, before.st_ino, before.st_size,
            getattr(before, "st_mtime_ns", None),
            getattr(before, "st_ctime_ns", None))
        stable_after = (
            after.st_dev, after.st_ino, after.st_size,
            getattr(after, "st_mtime_ns", None),
            getattr(after, "st_ctime_ns", None))
        data = b"".join(blocks)
        if stable_before != stable_after or len(data) != before.st_size:
            fail("{} changed while it was being read".format(context))
        return data
    except OSError as exc:
        fail("cannot read {} {}: {}".format(context, path, exc))
    finally:
        if descriptor >= 0:
            os.close(descriptor)
    return b""


def _read_canonical_json(
        path: Path, byte_cap: int,
        directory_fd: Optional[int] = None) -> Mapping[str, Any]:
    data = _read_regular_bytes(
        path, byte_cap, "work/rank summary", directory_fd)
    try:
        value = contract_api._load_json_bytes(data, "work/rank summary")
    except contract_api.ContractError as exc:
        fail("cannot parse artifact {}: {}".format(path, exc))
    if not isinstance(value, dict) or \
            data != (_canonical(value) + "\n").encode("utf-8"):
        fail("artifact is not canonical JSON: {}".format(path))
    return value


def load_completed_work_screen(
        contract: Mapping[str, Any], artifact_dir: Path) -> Mapping[str, Any]:
    """Strictly load one complete v3 work/rank campaign.

    This is the only supported ingestion seam for joining sidecar algebraic
    work with native recovery records.  It deliberately replays every raw
    construction identity instead of trusting the summary's self-hash.
    """
    artifact_dir = Path(artifact_dir)
    directory_fd = -1
    try:
        nofollow = getattr(os, "O_NOFOLLOW", 0)
        directory_flag = getattr(os, "O_DIRECTORY", 0)
        if nofollow == 0 or directory_flag == 0:
            fail("work/rank artifact directory cannot be opened fail-closed")
        directory_fd = os.open(
            str(artifact_dir), os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
            nofollow | directory_flag)
        directory_info = os.fstat(directory_fd)
        if not stat.S_ISDIR(directory_info.st_mode):
            fail("work/rank artifact directory must be a real directory")
        resolved_dir = artifact_dir.resolve(strict=True)
        resolved_info = resolved_dir.stat()
        if (directory_info.st_dev, directory_info.st_ino) != \
                (resolved_info.st_dev, resolved_info.st_ino):
            fail("work/rank artifact directory identity changed")
        summary = _read_canonical_json(
            Path(SUMMARY_NAME), 16 * 1024 * 1024, directory_fd)
        result_bytes = _read_regular_bytes(
            Path(RESULT_NAME), 64 * 1024 * 1024,
            "work/rank result stream", directory_fd)
    except OSError as exc:
        fail("work/rank artifact directory must be a real directory: {}"
             .format(exc))
    finally:
        if directory_fd >= 0:
            os.close(directory_fd)
    if set(summary) != set(SUMMARY_FIELDS) or \
            summary.get("schema") != SUMMARY_SCHEMA:
        fail("work/rank summary schema or fields changed")
    unsigned_summary = dict(summary)
    summary_sha256 = unsigned_summary.pop("summary_sha256", None)
    if not _is_sha256(summary_sha256) or \
            summary_sha256 != _sha256_json(unsigned_summary):
        fail("work/rank summary self-hash is invalid")
    expected_contract_sha256 = contract_api.contract_sha256(contract)
    expected_domain_sha256 = contract_api.recovery_domain_sha256(
        contract, "development")
    if (summary.get("contract_sha256") != expected_contract_sha256 or
            summary.get("recovery_domain_sha256") != expected_domain_sha256):
        fail("work/rank summary is bound to a different contract or domain")
    if (type(summary.get("record_count")) is not int or
            summary["record_count"] != EXPECTED_RECORDS or
            type(summary.get("invocation_count")) is not int or
            summary["invocation_count"] != EXPECTED_INVOCATIONS or
            summary.get("construction_seed_basis") != RAW_SEED_BASIS or
            summary.get("seed_schedule_sha256") !=
                RAW_SEED_SCHEDULE_SHA256 or
            summary.get("selection_performed") is not False or
            summary.get("timing_metrics_included") is not False):
        fail("work/rank summary policy or cardinality is invalid")

    provenance = summary.get("source_provenance")
    if not isinstance(provenance, dict) or \
            set(provenance) != set(SOURCE_PROVENANCE_FIELDS) or \
            summary.get("source_provenance_sha256") != \
                _sha256_json(provenance):
        fail("work/rank source provenance is invalid")
    source_git_commit = provenance.get("source_git_commit")
    if (not isinstance(source_git_commit, str) or
            re.fullmatch(r"[0-9a-f]{40}", source_git_commit) is None or
            provenance.get("worker_source_git_commit") !=
                source_git_commit or
            provenance.get("tracked_worktree_clean") is not True or
            provenance.get("untracked_source_clean") is not True or
            any(not _is_sha256(provenance.get(field)) for field in
                SOURCE_PROVENANCE_FIELDS - {
                    "source_git_commit", "tracked_worktree_clean",
                    "untracked_source_clean", "worker_source_git_commit"}) or
            provenance["tracked_status_sha256"] !=
                hashlib.sha256(b"").hexdigest() or
            provenance["untracked_source_sha256"] !=
                hashlib.sha256(b"").hexdigest() or
            provenance["worker_binary_sha256"] !=
                summary.get("worker_binary_sha256")):
        fail("work/rank source provenance fields are malformed")

    if (not result_bytes.endswith(b"\n") or b"\r" in result_bytes or
            hashlib.sha256(result_bytes).hexdigest() !=
                summary.get("result_stream_sha256")):
        fail("work/rank result stream hash, type, or framing is invalid")
    raw_lines = result_bytes.splitlines()
    if len(raw_lines) != EXPECTED_RECORDS:
        fail("work/rank result stream has the wrong cardinality")

    cell_map, trace_map, trace_identity = build_frozen_domain(contract)
    if _canonical(summary.get("trace_identity")) != \
            _canonical(trace_identity):
        fail("work/rank trace identity is invalid")
    arm_shapes = {
        arm: (dense, heavy, layout)
        for arm, dense, heavy, layout in ARMS
    }
    rows: List[Mapping[str, Any]] = []
    seen = set()
    hash_fields = (
        "arm_descriptor_sha256", "realized_construction_sha256",
        "seed_schedule_sha256", "cell_sha256", "frozen_trace_sha256",
        "frozen_packet_prefix_sha256", "command_sha256",
        "worker_stdout_sha256", "worker_binary_sha256",
        "source_provenance_sha256",
    )
    counter_fields = (
        "inactivated_columns", "binary_deficiency", "gf256_rank_gain",
        "heavy_shortfall", "block_xors", "gf256_muladds",
        "frozen_trace_attempted_candidates",
        "frozen_prefix_attempted_candidates",
    )
    for ordinal, raw_line in enumerate(raw_lines):
        try:
            row = contract_api._load_json_bytes(
                raw_line, "work/rank record {}".format(ordinal))
        except contract_api.ContractError as exc:
            fail("work/rank record is not JSON: {}".format(exc))
        if (not isinstance(row, dict) or set(row) != set(RECORD_FIELDS) or
                raw_line.decode("utf-8") != _canonical(row) or
                row.get("schema") != RECORD_SCHEMA or
                type(row.get("ordinal")) is not int or
                row.get("ordinal") != ordinal):
            fail("work/rank record schema, canonical form, or ordinal changed")
        if any(not _is_sha256(row.get(field)) for field in hash_fields):
            fail("work/rank record contains a malformed hash")
        arm = row.get("arm")
        if not isinstance(arm, str) or arm not in arm_shapes or \
                row.get("arm_descriptor_sha256") != \
                _arm_descriptor_sha256(arm):
            fail("work/rank record has an unknown raw arm descriptor")
        dense_rows, heavy_rows, dense_anchor_layout = arm_shapes[arm]
        K = row.get("K")
        overhead = row.get("overhead")
        trial = row.get("trial")
        schedule = row.get("schedule")
        if (type(K) is not int or K not in K_VALUES or
                type(overhead) is not int or overhead not in OVERHEADS or
                type(trial) is not int or not 0 <= trial < len(ROOTS) or
                not isinstance(schedule, str) or
                schedule not in {value[0] for value in STRATA}):
            fail("work/rank record has an invalid frozen cell coordinate")
        expected_loss_ppm = dict(STRATA)[schedule]
        staircase = row.get("staircase")
        source_hits = row.get("source_hits")
        if (row.get("phase") != "development" or
                type(row.get("block_bytes")) is not int or
                row["block_bytes"] != 2 or
                type(row.get("loss_ppm")) is not int or
                row["loss_ppm"] != expected_loss_ppm or
                row.get("loss_seed") != ROOTS[trial] or
                type(row.get("base_seed_attempt")) is not int or
                row["base_seed_attempt"] != trial or
                row.get("precode_attempt") != trial or
                row.get("packet_attempt") != trial or
                row.get("attempt_mode") != "exact" or
                row.get("construction_seed_basis") != RAW_SEED_BASIS or
                row.get("seed_schedule_sha256") !=
                    RAW_SEED_SCHEDULE_SHA256 or
                row.get("effective_precode_seed") !=
                    _effective_precode_seed(trial) or
                row.get("effective_packet_seed") !=
                    _effective_packet_seed(trial) or
                type(staircase) is not int or not 1 <= staircase <= 400 or
                row.get("binary_dense_rows") != dense_rows or
                row.get("gf256_heavy_rows") != heavy_rows or
                row.get("dense_anchor_layout") != dense_anchor_layout or
                type(source_hits) is not int or
                source_hits != (3 if K >= 10000 else 2) or
                row.get("dense_identity_corner") is not False or
                row.get("heavy_family") != "periodic-cauchy" or
                row.get("mix_count") != 3 or
                K + staircase + dense_rows + heavy_rows > 65535):
            fail("work/rank record has an invalid raw construction payload")
        raw_fields = {
            field: row[field] for field in (
                "construction_seed_basis", "seed_schedule_sha256",
                "precode_attempt", "packet_attempt",
                "effective_precode_seed", "effective_packet_seed",
                "staircase", "binary_dense_rows", "gf256_heavy_rows",
                "source_hits", "dense_anchor_layout",
                "dense_identity_corner", "heavy_family",
                "mix_count")
        }
        if row["realized_construction_sha256"] != \
                _raw_realized_construction_sha256(arm, K, 2, raw_fields):
            fail("work/rank realized construction hash is invalid")

        cell_key = (trial, schedule, K)
        cell = cell_map.get(cell_key)
        trace = trace_map.get(cell_key)
        if cell is None or trace is None:
            fail("work/rank record is outside the frozen domain")
        prefix = trace.packet_ids[:K + overhead]
        if (row.get("band") != cell["band"] or
                row.get("cell_sha256") != _sha256_json(cell) or
                row.get("frozen_trace_sha256") != trace.trace_sha256 or
                row.get("frozen_packet_prefix_sha256") !=
                    _packet_hash(prefix) or
                type(row.get("packet_count")) is not int or
                row["packet_count"] != K + overhead or
                row.get("frozen_trace_attempted_candidates") !=
                    trace.attempted_candidates or
                row.get("frozen_prefix_attempted_candidates") !=
                    trace.prefix_attempted_candidates[
                        OVERHEADS.index(overhead)]):
            fail("work/rank record trace identity is invalid")
        if (row.get("worker_binary_sha256") !=
                summary["worker_binary_sha256"] or
                row.get("source_provenance_sha256") !=
                summary["source_provenance_sha256"] or
                any(type(row.get(field)) is not int or
                    not 0 <= row[field] <= MASK64
                    for field in counter_fields) or
                any(type(row.get(field)) is not int or row[field] not in (0, 1)
                    for field in ("success", "rank_fail", "error")) or
                row["success"] + row["rank_fail"] + row["error"] != 1 or
                row["error"] != 0 or
                row["binary_deficiency"] > row["inactivated_columns"] or
                row["gf256_rank_gain"] > row["binary_deficiency"] or
                bool(row["success"]) !=
                    (row["gf256_rank_gain"] == row["binary_deficiency"]) or
                bool(row["rank_fail"]) !=
                    (row["gf256_rank_gain"] < row["binary_deficiency"])):
            fail("work/rank record counters or provenance are invalid")
        expected_shortfall = int(
            bool(row["rank_fail"]) and
            row["binary_deficiency"] <= heavy_rows and
            row["gf256_rank_gain"] < row["binary_deficiency"])
        if row["heavy_shortfall"] != expected_shortfall:
            fail("work/rank record heavy shortfall is invalid")
        coordinate = (arm, trial, schedule, K, overhead)
        if coordinate in seen:
            fail("work/rank result stream contains a duplicate coordinate")
        seen.add(coordinate)
        rows.append(row)

    expected_coordinates = {
        (arm, trial, schedule, K, overhead)
        for arm, _dense, _heavy, _layout in ARMS
        for trial in range(len(ROOTS))
        for schedule, _loss_ppm in STRATA
        for K in K_VALUES for overhead in OVERHEADS
    }
    if seen != expected_coordinates:
        fail("work/rank result stream omitted a frozen coordinate")
    if summary.get("work_domain_sha256") != _domain_sha256(rows):
        fail("work/rank domain hash is invalid")
    if _canonical(summary.get("arms")) != _canonical(_aggregate(rows)):
        fail("work/rank arm summary is invalid")

    receipts = []
    for invocation in make_invocations():
        invocation_rows = [
            row for row in rows
            if row["arm"] == invocation.arm and
               row["trial"] == invocation.trial and
               row["schedule"] == invocation.schedule
        ]
        commands = {row["command_sha256"] for row in invocation_rows}
        outputs = {row["worker_stdout_sha256"] for row in invocation_rows}
        if len(invocation_rows) != len(K_VALUES) * len(OVERHEADS) or \
                len(commands) != 1 or len(outputs) != 1:
            fail("work/rank invocation receipt cannot be reconstructed")
        receipts.append({
            "ordinal": invocation.ordinal,
            "arm": invocation.arm,
            "trial": invocation.trial,
            "schedule": invocation.schedule,
            "row_count": len(invocation_rows),
            "command_sha256": next(iter(commands)),
            "worker_stdout_sha256": next(iter(outputs)),
        })
    summary_receipts = summary.get("invocations")
    if (not isinstance(summary_receipts, list) or
            any(not isinstance(value, dict) or
                set(value) != set(INVOCATION_RECEIPT_FIELDS)
                for value in summary_receipts) or
            _canonical(summary_receipts) != _canonical(receipts)):
        fail("work/rank invocation receipts are invalid")
    return {
        "summary": summary,
        "rows": tuple(rows),
        "summary_sha256": summary_sha256,
        "result_stream_sha256": summary["result_stream_sha256"],
        "work_domain_sha256": summary["work_domain_sha256"],
        "source_git_commit": provenance["source_git_commit"],
    }


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
