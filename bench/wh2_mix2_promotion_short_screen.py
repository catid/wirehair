#!/usr/bin/env python3
"""Run the frozen v7 production-basis MIX2 bounded promotion screen.

This is deliberately a short screen, not the all-K recovery census and not
the final Wirehair1 timing gate.  The v7 repair worker first selects one uint8
construction attempt per preregistered K on the 18 frozen selection cells.
The ``precodefail`` test hook then byte-verifies overhead-zero recovery and
measures deterministic work on three fresh roots.  Each benchmark coordinate
uses two observations per arm in a deterministic ABBA or BAAB order.

The controller owns no equation implementation.  It validates the native
binary protocols, freezes the exact selection stream and input roster before
the first benchmark solve, and writes canonical JSON artifacts
transactionally.  Worker validation commands are intentionally reserved for
the untouched final all-K holdout and are never issued by this bounded screen.
"""

from __future__ import annotations

import argparse
import ctypes
import csv
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation, getcontext
import errno
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import time
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


getcontext().prec = 50

SCREEN_SCHEMA = "wirehair.wh2.mix2-promotion-short-screen.v2"
INPUT_SCHEMA = "wirehair.wh2.mix2-promotion-short-screen-input.v2"
RESULT_SCHEMA = "wirehair.wh2.mix2-promotion-short-screen-record.v2"
SUMMARY_SCHEMA = "wirehair.wh2.mix2-promotion-short-screen-summary.v2"
ATTEMPT_STREAM_SCHEMA = \
    "wirehair.wh2.mix2-promotion-short-screen-attempt-stream.v2"

WORKER_DESCRIPTION_SCHEMA = \
    "wirehair.wh2.mix2-seed-repair-worker-description.v3"
WORKER_DERIVATION_SCHEMA = \
    "wirehair.wh2.mix2-seed-repair-derivation-record.v3"
WORKER_VALIDATION_SCHEMA = \
    "wirehair.wh2.mix2-seed-repair-validation-record.v3"
WORKER_CONTRACT_SCHEMA = "wirehair.wh2.mix2-seed-repair-contract.v7"
WORKER_SCHEMA = "wirehair.wh2.mix2-seed-repair-worker.v3"
CANDIDATE_PROFILE_SCHEMA = "wirehair.wh2.mix2-production-profile.v1"

RESULT_NAME = "promotion-short-screen-results.jsonl"
SUMMARY_NAME = "promotion-short-screen-summary.json"
INPUT_NAME = "promotion-short-screen-input.json"
ATTEMPT_NAME = "promotion-short-screen-attempts.jsonl"

K_VALUES = (
    2, 3, 4, 5, 6, 8, 16, 32, 64, 100, 101, 128, 256, 512, 513,
    1000, 1001, 2048, 4096, 5000, 5001, 8192, 10000, 10001,
    16384, 20000, 20001, 32768, 49152, 64000,
)
SELECTION_ROOTS = (
    "0xd1b54a32d192ed03",
    "0x94d049bb133111eb",
    "0x8538ecb5bd456ea3",
    "0xc0ac29b7c97c50dd",
    "0x3f84d5b5b5470917",
    "0x9216d5d98979fb1b",
)
ROOTS = (
    "0x22ef7f82b3d08e8d",
    "0x9e5241defc95255c",
    "0xdd0a4e8205da8ed0",
)
BENCHMARK_ROOT_NAMESPACE_TEMPLATE = \
    "wirehair2-two07-mix2-graph-b2-short-v7:holdout-root:{i}"
BENCHMARK_ROOT_FULL_SHA256 = (
    "22ef7f82b3d08e8da867549574a58821f6992e1ec57f6aa8b1015a65b0b1eb6f",
    "9e5241defc95255cc27fcabfa4485097fd25e4f7f44e47daaa1d2b80a21d51f7",
    "dd0a4e8205da8ed04035547611843b364b52d969c7b00bf87b87230cb8c98695",
)
DISCARDED_V6_ROOTS = (
    "0x343c26b8a0a06468",
    "0x5b6e50735153e7db",
    "0x738c94f4e3310746",
)
FINAL_VALIDATION_ROOTS = (
    "0xb501025fdce63900",
    "0x7fb960494dece7de",
    "0x6ad0017d0069e483",
)
SCHEDULES = ("burst", "adversarial", "repair-only")
VALIDATION_ROSTER_SCHEMA = \
    "wirehair.wh2.mix2-seed-repair-validation-roster.v1"
VALIDATION_ROSTER_SHA256 = \
    "030bb1c51e21777266edd4c2349d4a81ccf6e79e2fe4ed9eb75856e16f3387c7"
ARMS = ("current_disabled_mix3", "candidate_two07_mix2")
TIMING_ORDERS = ("ABBA", "BAAB")
OBSERVATIONS_PER_ARM = 2
BLOCK_BYTES = 2
LOSS_PPM = 500000
OVERHEAD = 0
BINARY_DENSE_ROWS = 12
GF256_HEAVY_ROWS = 12
EXPECTED_CELL_COUNT = len(K_VALUES) * len(ROOTS) * len(SCHEDULES)
EXPECTED_RECORD_COUNT = (
    EXPECTED_CELL_COUNT * len(ARMS) * OBSERVATIONS_PER_ARM)

XOR_RATIO_MAX = Decimal("0.9829453304")
GF256_MULADD_RATIO_MAX = Decimal("1")
INACTIVATION_RATIO_MAX = Decimal("1")
SOLVE_TIME_RATIO_MAX = Decimal("1")
BUILD_TIME_RATIO_MAX = Decimal("1.05")

PRECODE_ATTEMPT_STRIDE = 0x9E3779B97F4A7C15
PACKET_ATTEMPT_STRIDE = 0x9E3779B9
MASK64 = (1 << 64) - 1
MASK32 = (1 << 32) - 1
SEED_SCHEDULE_SHA256 = (
    "fe8101781024b4a30797d66cb1512640"
    "e6a939586ad5e32f73c5b7ff6f411294"
)
ZERO_SHA256 = "0" * 64
EXPECTED_CANDIDATE_PROFILE_SHA256 = (
    "90233b44a0893f96c1a18c19aa61ada052c935a48c6bf7d6a2813065856651f0"
)
# Updated only when the preregistered contract object below intentionally
# changes.  Keeping the digest literal prevents a silent roster/gate edit from
# merely blessing itself with a newly computed hash.
EXPECTED_CONTRACT_SHA256 = (
    "574bef7638a51e34c48b118613833df554e3276ae0cea952db787ed0d301af8e"
)

REPO_ROOT = Path(__file__).resolve().parent.parent
SOURCE_PATHS = (
    "CMakeLists.txt",
    "bench/wh2_mix2_promotion_short_screen.py",
    "bench/Wh2Mix2SeedRepairWorker.cpp",
    "bench/Wh2FrozenTrace.cpp",
    "bench/Wh2FrozenTrace.h",
    "bench/Wh2NativeCodec.cpp",
    "bench/Wh2NativeCodec.h",
    "codec/WirehairV2Bench.cpp",
    "codec/WirehairV2Precode.cpp",
    "codec/WirehairV2Precode.h",
    "codec/WirehairV2PrecodeEncode.cpp",
    "codec/WirehairV2PrecodeEncode.h",
    "codec/WirehairV2Profile.cpp",
    "codec/WirehairV2Seeds.h",
    "include/wirehair/wirehair.h",
    "WirehairDenseFixups.inc",
    "WirehairPeelFixups.inc",
)
SOURCE_STATUS_PATHS = (
    "CMakeLists.txt", "cmake", "codec", "include", "tables",
    "wirehair.cpp", "gf256.cpp", "gf256.h", "WirehairCodec.cpp",
    "WirehairCodec.h", "WirehairCompiler.h", "WirehairEnvironment.h",
    "WirehairHeavy.h", "WirehairTools.cpp", "WirehairTools.h",
    "WirehairDenseFixups.inc", "WirehairPeelFixups.inc",
    "bench/Wh2Mix2SeedRepairWorker.cpp", "bench/Wh2FrozenTrace.cpp",
    "bench/Wh2FrozenTrace.h", "bench/Wh2NativeCodec.cpp",
    "bench/Wh2NativeCodec.h",
    "bench/wh2_mix2_promotion_short_screen.py",
)

MAX_DESCRIPTION_BYTES = 64 * 1024
MAX_WORKER_STDOUT_BYTES = 32 * 1024 * 1024
MAX_BENCH_STDOUT_BYTES = 256 * 1024
MAX_STDERR_BYTES = 64 * 1024
WORKER_TIMEOUT_SECONDS = 3600
BENCH_TIMEOUT_SECONDS = 600
TOTAL_WALL_SECONDS = 4 * 60 * 60

SHA256_TOKEN = re.compile(r"[0-9a-f]{64}\Z")
COMMIT_TOKEN = re.compile(r"[0-9a-f]{40}\Z")
UINT_TOKEN = re.compile(r"(?:0|[1-9][0-9]*)\Z")
HEX32_TOKEN = re.compile(r"0x[0-9a-f]{8}\Z")
HEX64_TOKEN = re.compile(r"0x[0-9a-f]{16}\Z")
DECIMAL_TOKEN = re.compile(r"(?:0|[1-9][0-9]*)\.[0-9]{3}\Z")
FAIL_RATE_TOKEN = re.compile(r"(?:0|1)\.[0-9]{8}\Z")

_LIBC = ctypes.CDLL(None, use_errno=True)
_RENAMEAT2 = getattr(_LIBC, "renameat2", None)
if _RENAMEAT2 is not None:
    _RENAMEAT2.argtypes = (
        ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p,
        ctypes.c_uint)
    _RENAMEAT2.restype = ctypes.c_int
_AT_FDCWD = -100
_RENAME_NOREPLACE = 1

CANDIDATE_PROFILE = {
    "arm": "wirehair2_two07_mix2_graph_b2_v1",
    "binary_dense_rows": 12,
    "codec": "wirehair2_candidate",
    "construction_seed_basis": "production-normalized-b2-v1",
    "dense_anchor_layout": "two07",
    "dense_identity_corner": False,
    "field": "GF(256)",
    "gf256_heavy_rows": 12,
    "graph_seed_block_bytes": 2,
    "heavy_family": "periodic-cauchy",
    "mix_count": 2,
    "packet_seed_salt": "0x76327265636f7665",
    "precode_seed_salt": "0x763263707265636f",
    "schema": CANDIDATE_PROFILE_SCHEMA,
    "seed_schedule_sha256": SEED_SCHEDULE_SHA256,
    "source_hits": "certified-by-K",
}

CELL_FIELDS = frozenset((
    "attempted_candidates", "cell_ordinal", "construction_attempt",
    "decoded_extra", "loss_ppm", "loss_root", "outcome", "result_code",
    "root_index", "schedule", "trace_sha256",
))
DERIVATION_FIELDS = frozenset((
    "K", "base_packet_seed", "base_precode_seed",
    "candidate_profile_sha256", "effective_packet_seed",
    "effective_precode_seed", "lower_attempt_failure_witnesses", "mode",
    "ordinal", "schema", "selected_attempt", "selected_successes",
    "source_sha256", "worker_binary_sha256",
))
DESCRIPTION_FIELDS = frozenset((
    "binary_sha256", "candidate_profile", "candidate_profile_sha256",
    "contract_schema", "derivation_schema", "protocol", "schema",
    "source_git_commit", "validation_roster_schema",
    "validation_roster_sha256", "validation_schema", "worker_schema",
))

CSV_HEADER = (
    "N", "bb", "heavy_family", "mix_count", "staircase",
    "binary_dense_rows", "gf256_heavy_rows", "source_hits",
    "dense_identity_corner", "overhead", "trials", "success",
    "rank_fail", "error", "fail_rate", "inact_mu", "inact_max",
    "binary_def_mu", "binary_def_max", "heavy_gain_mu", "heavy_gain_min",
    "heavy_shortfall", "solve_ms_mu", "build_ms_mu", "peel_ms_mu",
    "project_ms_mu", "residual_ms_mu", "backsub_ms_mu", "seed_attempt",
    "block_xors_mu", "block_muladds_mu", "first_rank_fail",
    "binary_def_hist", "heavy_gain_hist", "failure_trials",
    "active_packet_peel_seed_xor", "precode_attempt", "packet_attempt",
    "attempt_mode", "construction_seed_basis", "seed_schedule_sha256",
    "effective_precode_seed", "effective_packet_seed",
)
METADATA_KEYS = frozenset((
    "trials", "threads", "loss", "seed", "source_hits_override",
    "packet_peel_seed_xor", "binary_dense_rows_override",
    "gf256_heavy_rows_override", "dense_anchor_layout",
    "odd_packet_peel_seed_xor", "packet_row_seed_multiplier",
    "packet_row_seed_avalanche", "seed_block_bytes_override",
    "overhead_stream", "full_payload_solve", "schedule",
    "cold_solve_wide_xor", "exact_attempt_mode",
    "exact_precode_attempt", "exact_packet_attempt",
    "construction_seed_basis", "seed_schedule_sha256",
    "source_git_commit",
))


class ShortScreenError(RuntimeError):
    """The frozen confirmation cannot be executed or admitted safely."""


def fail(message: str) -> None:
    raise ShortScreenError(message)


def canonical_json(value: Any) -> str:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False)


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_text(value: str) -> str:
    return sha256_bytes(value.encode("utf-8"))


def validation_roster_sha256() -> str:
    return sha256_text(canonical_json({
        "schema": VALIDATION_ROSTER_SCHEMA,
        "cell_order": "root-major-then-schedule",
        "root_count": len(FINAL_VALIDATION_ROOTS),
        "roots": list(FINAL_VALIDATION_ROOTS),
        "schedule_count": len(SCHEDULES),
        "schedules": list(SCHEDULES),
        "cell_count": len(FINAL_VALIDATION_ROOTS) * len(SCHEDULES),
    }))


def self_hash(value: Mapping[str, Any], field: str) -> str:
    body = dict(value)
    body.pop(field, None)
    return sha256_text(canonical_json(body))


def _exact_fields(value: Mapping[str, Any], expected: Iterable[str],
                  context: str) -> None:
    if set(value) != set(expected):
        fail("{} fields differ from the frozen schema".format(context))


def _uint(value: Any, maximum: int, context: str) -> int:
    if type(value) is not int or value < 0 or value > maximum:
        fail("{} is not a bounded unsigned integer".format(context))
    return value


def _uint_text(value: str, maximum: int, context: str) -> int:
    if not UINT_TOKEN.fullmatch(value):
        fail("{} is not a canonical unsigned integer".format(context))
    parsed = int(value, 10)
    if parsed > maximum:
        fail("{} exceeds its frozen bound".format(context))
    return parsed


def _hex_value(value: Any, regex: re.Pattern[str], context: str) -> int:
    if type(value) is not str or not regex.fullmatch(value):
        fail("{} is not a canonical fixed-width hexadecimal value".format(
            context))
    return int(value, 16)


def _sha(value: Any, context: str) -> str:
    if type(value) is not str or not SHA256_TOKEN.fullmatch(value):
        fail("{} is not a lowercase SHA-256".format(context))
    return value


def _decimal(value: str, context: str,
             pattern: re.Pattern[str] = DECIMAL_TOKEN) -> Decimal:
    if not pattern.fullmatch(value):
        fail("{} is not a canonical decimal".format(context))
    try:
        parsed = Decimal(value)
    except InvalidOperation:
        fail("{} is not a finite decimal".format(context))
    if not parsed.is_finite() or parsed < 0:
        fail("{} is not a nonnegative finite decimal".format(context))
    return parsed


def _integral_decimal(value: str, context: str) -> int:
    parsed = _decimal(value, context)
    integral = parsed.to_integral_value()
    if parsed != integral:
        fail("{} is not integral for a one-trial cell".format(context))
    return int(integral)


def _stable_file_bytes(path: Path, maximum: Optional[int] = None) -> bytes:
    flags = os.O_RDONLY
    if hasattr(os, "O_CLOEXEC"):
        flags |= os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        descriptor = os.open(str(path), flags)
    except OSError as exc:
        fail("cannot open {} safely: {}".format(path, exc))
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            fail("{} is not a regular file".format(path))
        if maximum is not None and before.st_size > maximum:
            fail("{} exceeds its size bound".format(path))
        chunks: List[bytes] = []
        total = 0
        while True:
            chunk = os.read(descriptor, 1024 * 1024)
            if not chunk:
                break
            total += len(chunk)
            if maximum is not None and total > maximum:
                fail("{} grew beyond its size bound".format(path))
            chunks.append(chunk)
        after = os.fstat(descriptor)
        if (before.st_dev, before.st_ino, before.st_size,
                before.st_mtime_ns) != (after.st_dev, after.st_ino,
                                        after.st_size, after.st_mtime_ns):
            fail("{} changed while it was read".format(path))
        return b"".join(chunks)
    finally:
        os.close(descriptor)


def _hash_open_file(descriptor: int, context: str) -> Tuple[str, int]:
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            fail("{} is no longer a regular file".format(context))
        digest = hashlib.sha256()
        offset = 0
        while offset < before.st_size:
            chunk = os.pread(
                descriptor, min(1024 * 1024, before.st_size - offset), offset)
            if not chunk:
                fail("{} was truncated while hashing".format(context))
            digest.update(chunk)
            offset += len(chunk)
        if os.pread(descriptor, 1, offset):
            fail("{} grew while hashing".format(context))
        after = os.fstat(descriptor)
        identity = lambda value: (
            value.st_dev, value.st_ino, value.st_size,
            value.st_mtime_ns, value.st_ctime_ns)
        if identity(before) != identity(after):
            fail("{} changed while hashing".format(context))
        return digest.hexdigest(), before.st_size
    except OSError as exc:
        fail("cannot hash {}: {}".format(context, exc))
    return "", 0


@dataclass
class PinnedExecutable:
    path: Path
    descriptor: int
    sha256: str
    size: int
    label: str

    def receipt(self) -> Mapping[str, Any]:
        if self.descriptor < 0:
            fail("{} executable is already closed".format(self.label))
        digest, size = _hash_open_file(
            self.descriptor, "{} executable".format(self.label))
        if digest != self.sha256 or size != self.size:
            fail("{} executable changed after it was pinned".format(self.label))
        return {
            "path": str(self.path),
            "sha256": digest,
            "size": size,
        }

    def close(self) -> None:
        if self.descriptor >= 0:
            os.close(self.descriptor)
            self.descriptor = -1

    def __enter__(self) -> "PinnedExecutable":
        return self

    def __exit__(self, exc_type: Any, exc: Any, traceback: Any) -> None:
        del exc_type, exc, traceback
        self.close()


def _open_binary(path: Path, label: str) -> PinnedExecutable:
    descriptor = -1
    try:
        original = path.lstat()
        if stat.S_ISLNK(original.st_mode) or not stat.S_ISREG(original.st_mode):
            fail("{} binary must be a non-symlink regular file".format(label))
        resolved = path.resolve(strict=True)
        nofollow = getattr(os, "O_NOFOLLOW", 0)
        if nofollow == 0:
            fail("{} binary cannot be opened fail-closed".format(label))
        descriptor = os.open(
            str(resolved), os.O_RDONLY | os.O_CLOEXEC | nofollow)
        opened = os.fstat(descriptor)
        current = os.stat(str(resolved), follow_symlinks=False)
        if (not stat.S_ISREG(opened.st_mode) or
                not opened.st_mode & 0o111 or
                (original.st_dev, original.st_ino) !=
                    (opened.st_dev, opened.st_ino) or
                (current.st_dev, current.st_ino) !=
                    (opened.st_dev, opened.st_ino)):
            fail("{} binary is not one stable executable regular file".format(
                label))
        digest, size = _hash_open_file(
            descriptor, "{} executable".format(label))
        if size <= 0:
            fail("{} binary is empty".format(label))
        pinned = PinnedExecutable(resolved, descriptor, digest, size, label)
        descriptor = -1
        return pinned
    except OSError as exc:
        fail("cannot inspect {} binary: {}".format(label, exc))
    finally:
        if descriptor >= 0:
            os.close(descriptor)
    raise AssertionError("unreachable")


def _run_process(argv: Sequence[str], input_bytes: Optional[bytes],
                 timeout: int, stdout_limit: int,
                 context: str,
                 executable_fd: Optional[int] = None) -> Tuple[bytes, bytes]:
    environment = dict(os.environ)
    environment["LC_ALL"] = "C"
    environment["LANG"] = "C"
    execution_options: Dict[str, Any] = {}
    if executable_fd is not None:
        if type(executable_fd) is not int or executable_fd < 0:
            fail("{} executable descriptor is invalid".format(context))
        execution_options = {
            "executable": "/proc/self/fd/{}".format(executable_fd),
            "pass_fds": (executable_fd,),
        }
    try:
        completed = subprocess.run(
            list(argv), input=input_bytes, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, timeout=timeout, check=False,
            env=environment, **execution_options)
    except (OSError, subprocess.SubprocessError) as exc:
        fail("{} failed to execute: {}".format(context, exc))
    if len(completed.stdout) > stdout_limit:
        fail("{} stdout exceeds its frozen bound".format(context))
    if len(completed.stderr) > MAX_STDERR_BYTES:
        fail("{} stderr exceeds its frozen bound".format(context))
    if completed.returncode != 0:
        diagnostic = completed.stderr.decode("utf-8", "replace").strip()
        fail("{} exited {}: {}".format(
            context, completed.returncode, diagnostic[:1000]))
    return completed.stdout, completed.stderr


def candidate_profile_sha256() -> str:
    digest = sha256_text(canonical_json(CANDIDATE_PROFILE))
    if digest != EXPECTED_CANDIDATE_PROFILE_SHA256:
        fail("candidate profile object differs from its frozen digest")
    return digest


def _validate_benchmark_root_derivation() -> None:
    if len(ROOTS) != 3 or len(BENCHMARK_ROOT_FULL_SHA256) != len(ROOTS):
        fail("benchmark-root derivation roster has the wrong size")
    for index, (root, expected_digest) in enumerate(zip(
            ROOTS, BENCHMARK_ROOT_FULL_SHA256)):
        preimage = BENCHMARK_ROOT_NAMESPACE_TEMPLATE.format(i=index)
        try:
            encoded = preimage.encode("ascii")
        except UnicodeEncodeError:
            fail("benchmark-root namespace is not raw ASCII")
        digest = sha256_bytes(encoded)
        if digest != expected_digest or root != "0x" + digest[:16]:
            fail("benchmark root {} differs from its SHA-256 preimage".format(
                index))
    root_sets = (
        set(ROOTS), set(SELECTION_ROOTS), set(DISCARDED_V6_ROOTS),
        set(FINAL_VALIDATION_ROOTS),
    )
    if any(len(values) != len(tuple_values) for values, tuple_values in zip(
            root_sets,
            (ROOTS, SELECTION_ROOTS, DISCARDED_V6_ROOTS,
             FINAL_VALIDATION_ROOTS))):
        fail("a frozen root roster contains duplicates")
    for left in range(len(root_sets)):
        for right in range(left + 1, len(root_sets)):
            if root_sets[left] & root_sets[right]:
                fail("frozen selection, short-screen, discarded-v6, and "
                     "final-validation roots are not pairwise disjoint")


def _contract_body() -> Mapping[str, Any]:
    _validate_benchmark_root_derivation()
    return {
        "schema": SCREEN_SCHEMA,
        "preregistration": "wirehair-sxvz.16.1.20.38@2026-08-30T22:23Z",
        "K": list(K_VALUES),
        "block_bytes": BLOCK_BYTES,
        "loss_ppm": LOSS_PPM,
        "overhead": OVERHEAD,
        "selection_roots": list(SELECTION_ROOTS),
        "benchmark_roots": list(ROOTS),
        "benchmark_root_derivation": {
            "namespace_template": BENCHMARK_ROOT_NAMESPACE_TEMPLATE,
            "indices": [0, 1, 2],
            "preimage_encoding": "raw ASCII",
            "full_sha256": list(BENCHMARK_ROOT_FULL_SHA256),
            "root_derivation":
                "first 8 SHA-256 digest bytes interpreted big-endian and "
                "formatted as lowercase fixed-width 0x hexadecimal",
            "pairwise_distinct": True,
            "disjoint_from_selection_roots": True,
            "disjoint_from_final_validation_roots": True,
            "discarded_v6_roots": list(DISCARDED_V6_ROOTS),
            "discarded_v6_roots_excluded": True,
        },
        "schedules": list(SCHEDULES),
        "arms": list(ARMS),
        "arm_definitions": {
            "current_disabled_mix3": {
                "configuration": "exact current production selection",
                "dense_anchor_layout": "disabled",
                "mix_count": 3,
                "construction_attempt": "SelectSystematicConfiguration",
            },
            "candidate_two07_mix2": {
                "configuration": "normalized production candidate",
                "graph_seed_basis": "SelectSeedProfile(K,2)",
                "dense_anchor_layout": "two07",
                "mix_count": 2,
                "construction_attempt":
                    "exact v7 18-cell-selected uint8 attempt indexed by K",
            },
        },
        "candidate_profile": dict(CANDIDATE_PROFILE),
        "candidate_profile_sha256": candidate_profile_sha256(),
        "attempt_selection": {
            "commands": 30,
            "cells_per_K": 18,
            "cell_order": "selection-root-major-then-schedule",
            "rule":
                "first uint8 attempt passing all 18 frozen selection cells; "
                "exact selected attempt only in the bounded screen",
            "stream": "canonical worker derivation JSONL frozen before the "
                      "first benchmark invocation",
        },
        "aggregation": {
            "deterministic_work":
                "one observation per arm/coordinate after exact equality "
                "with the arm's repeated observation",
            "timing":
                "sums both observations per arm over paired common-success "
                "coordinates",
        },
        "payload_work":
            "two two-byte full-payload solves per arm/coordinate; each "
            "benchmark solve byte-verifies overhead-zero recovery",
        "evidence_semantics": {
            "attempt_selection_provenance":
                "the v7 worker's 18 selection-cell receipts are selection "
                "provenance only and contribute no promotion recovery or "
                "work evidence",
            "benchmark_loss_trace":
                "precodefail emits no packet IDs or trace hash; benchmark "
                "arms are paired by identical fresh requested coordinates, "
                "not by an observed or hashed benchmark trace",
            "final_validation_holdout":
                "the worker V command and its final all-K roots are reserved "
                "for the later final campaign and are not invoked here",
        },
        "timing_protocol": {
            "orders": list(TIMING_ORDERS),
            "order_rule":
                "ABBA for even canonical benchmark-coordinate ordinals and "
                "BAAB for odd ordinals; A=current, B=candidate",
            "observations_per_arm_per_coordinate": OBSERVATIONS_PER_ARM,
            "cpu_affinity_required": False,
            "scope":
                "bounded same-process-order screen only; not the final "
                "Wirehair1 timing gate",
        },
        "gates": {
            "candidate_benchmark_weak_observations_max": 0,
            "candidate_benchmark_weak_observations_not_above_control": True,
            "repeated_deterministic_work_equal": True,
            "block_xor_ratio_max": format(XOR_RATIO_MAX, "f"),
            "gf256_muladd_ratio_max": format(
                GF256_MULADD_RATIO_MAX, "f"),
            "inactivation_ratio_max": format(
                INACTIVATION_RATIO_MAX, "f"),
            "solve_time_ratio_max": format(SOLVE_TIME_RATIO_MAX, "f"),
            "build_time_ratio_max": format(BUILD_TIME_RATIO_MAX, "f"),
        },
        "complexity": {
            "allocation_count_delta": 0,
            "asymptotic_construction": "O(K)",
            "asymptotic_storage": "O(K)",
            "production_map_lookups_per_construction": 1,
            "online_retries": 0,
            "extra_dense_basis_capacity_entries":
                "ceil((K+S+12)/2)-2 uint32 entries",
        },
        "scope_exclusions": [
            "all-K recovery gate",
            "final Wirehair1 timing gate",
            "worker final-holdout validation",
        ],
    }


def contract_description() -> Mapping[str, Any]:
    value = dict(_contract_body())
    digest = sha256_text(canonical_json(value))
    if digest != EXPECTED_CONTRACT_SHA256:
        fail("short-screen contract differs from its frozen digest")
    value["contract_sha256"] = digest
    return value


def _validate_description(value: Any, binary_sha256: str,
                          source_commit: str) -> Mapping[str, Any]:
    if type(value) is not dict:
        fail("repair-worker description is not an object")
    _exact_fields(value, DESCRIPTION_FIELDS, "repair-worker description")
    if value["schema"] != WORKER_DESCRIPTION_SCHEMA or \
            value["contract_schema"] != WORKER_CONTRACT_SCHEMA or \
            value["derivation_schema"] != WORKER_DERIVATION_SCHEMA or \
            value["validation_schema"] != WORKER_VALIDATION_SCHEMA or \
            value["validation_roster_schema"] != VALIDATION_ROSTER_SCHEMA or \
            value["validation_roster_sha256"] != \
                VALIDATION_ROSTER_SHA256 or \
            validation_roster_sha256() != VALIDATION_ROSTER_SHA256 or \
            value["worker_schema"] != WORKER_SCHEMA or \
            value["protocol"] != "D ordinal K | V ordinal K attempt | Q":
        fail("repair-worker description protocol differs from v7")
    if value["binary_sha256"] != binary_sha256:
        fail("repair-worker self hash differs from the executable hash")
    if value["source_git_commit"] != source_commit:
        fail("repair-worker source commit differs from the frozen source")
    if (type(value["candidate_profile"]) is not dict or
            canonical_json(value["candidate_profile"]) !=
            canonical_json(CANDIDATE_PROFILE)):
        fail("repair-worker candidate profile differs from preregistration")
    expected_profile_sha = candidate_profile_sha256()
    if value["candidate_profile_sha256"] != expected_profile_sha:
        fail("repair-worker candidate profile hash is inconsistent")
    return value


def read_worker_description(worker: Path, binary_sha256: str,
                            source_commit: str,
                            executable_fd: Optional[int] = None) \
        -> Mapping[str, Any]:
    stdout, stderr = _run_process(
        [str(worker), "--describe"], None, BENCH_TIMEOUT_SECONDS,
        MAX_DESCRIPTION_BYTES, "repair-worker description", executable_fd)
    if stderr:
        fail("repair-worker description wrote stderr")
    try:
        text = stdout.decode("utf-8")
    except UnicodeDecodeError:
        fail("repair-worker description is not UTF-8")
    lines = text.splitlines()
    if len(lines) != 1 or not lines[0]:
        fail("repair-worker description must contain exactly one JSON line")
    if not text.endswith("\n") or "\r" in text:
        fail("repair-worker description is not newline-canonical")
    try:
        value = json.loads(lines[0])
    except (TypeError, ValueError) as exc:
        fail("repair-worker description is invalid JSON: {}".format(exc))
    if lines[0] != canonical_json(value):
        fail("repair-worker description is not canonical JSON")
    return _validate_description(value, binary_sha256, source_commit)


def _effective_precode_seed(base: str, attempt: int) -> str:
    value = (_hex_value(base, HEX64_TOKEN, "base precode seed") +
             attempt * PRECODE_ATTEMPT_STRIDE) & MASK64
    return "0x{:016x}".format(value)


def _effective_packet_seed(base: str, attempt: int) -> str:
    value = (_hex_value(base, HEX32_TOKEN, "base packet seed") +
             attempt * PACKET_ATTEMPT_STRIDE) & MASK32
    return "0x{:08x}".format(value)


def _validate_worker_cell(cell: Any, attempt: int, context: str,
                          require_success: bool,
                          roots: Sequence[str] = SELECTION_ROOTS) \
        -> Mapping[str, Any]:
    if type(cell) is not dict:
        fail("{} is not an object".format(context))
    _exact_fields(cell, CELL_FIELDS, context)
    cell_count = len(roots) * len(SCHEDULES)
    ordinal = _uint(
        cell["cell_ordinal"], cell_count - 1, context + " cell ordinal")
    root_index = _uint(
        cell["root_index"], len(roots) - 1, context + " root index")
    schedule_index = ordinal % len(SCHEDULES)
    if root_index != ordinal // len(SCHEDULES) or \
            cell["loss_root"] != roots[root_index] or \
            cell["schedule"] != SCHEDULES[schedule_index]:
        fail("{} is outside the frozen root/schedule roster".format(context))
    loss_ppm = _uint(cell["loss_ppm"], 1000000,
                     context + " loss ppm")
    construction_attempt = _uint(
        cell["construction_attempt"], 255,
        context + " construction attempt")
    if loss_ppm != LOSS_PPM or construction_attempt != attempt:
        fail("{} has the wrong loss or construction attempt".format(context))
    attempted_candidates = _uint(
        cell["attempted_candidates"], (1 << 63) - 1,
        context + " attempted candidates")
    if attempted_candidates == 0:
        fail("{} has an empty packet trace".format(context))
    _sha(cell["trace_sha256"],
         context + " attempt-selection trace hash")
    decoded = cell["decoded_extra"]
    if decoded is not None:
        _uint(decoded, 4, context + " decoded extra")
    result_code = _uint(
        cell["result_code"], (1 << 31) - 1, context + " result code")
    if require_success:
        if cell["outcome"] != "success" or result_code != 0 or \
                decoded != 0:
            fail("{} is not an overhead-zero success".format(context))
    elif cell["outcome"] == "need_more_at_oh0":
        if decoded is None or decoded == 0 or result_code != 0:
            fail("{} has an inconsistent overhead-zero failure".format(
                context))
    elif cell["outcome"] == "need_more_at_cap":
        if decoded is not None or result_code != 1:
            fail("{} has an inconsistent cap failure".format(context))
    elif cell["outcome"] == "construct_failed":
        if decoded is not None or result_code != 4:
            fail("{} has an inconsistent construction/cap failure".format(
                context))
    else:
        fail("{} is not a worker-emittable bounded failure".format(context))
    return cell


def validate_derivation_record(value: Any, expected_K: int,
                               worker_binary_sha256: str) -> Mapping[str, Any]:
    context = "repair derivation K={}".format(expected_K)
    if type(value) is not dict:
        fail("{} is not an object".format(context))
    _exact_fields(value, DERIVATION_FIELDS, context)
    K = _uint(value["K"], 64000, context + " K")
    ordinal = _uint(value["ordinal"], 63998, context + " ordinal")
    if value["schema"] != WORKER_DERIVATION_SCHEMA or \
            value["mode"] != "derive" or K != expected_K or \
            ordinal != expected_K - 2:
        fail("{} has the wrong identity".format(context))
    if value["worker_binary_sha256"] != worker_binary_sha256:
        fail("{} has the wrong worker binary hash".format(context))
    if value["candidate_profile_sha256"] != candidate_profile_sha256():
        fail("{} has the wrong candidate profile hash".format(context))
    _sha(value["source_sha256"], context + " source hash")
    _hex_value(value["base_precode_seed"], HEX64_TOKEN,
               context + " base precode seed")
    _hex_value(value["base_packet_seed"], HEX32_TOKEN,
               context + " base packet seed")
    attempt = _uint(value["selected_attempt"], 255,
                    context + " selected attempt")
    if value["effective_precode_seed"] != _effective_precode_seed(
            value["base_precode_seed"], attempt) or \
            value["effective_packet_seed"] != _effective_packet_seed(
                value["base_packet_seed"], attempt):
        fail("{} effective seeds do not implement the frozen schedule".format(
            context))
    lower = value["lower_attempt_failure_witnesses"]
    successes = value["selected_successes"]
    if type(lower) is not list or len(lower) != attempt:
        fail("{} omits a lower-attempt witness".format(context))
    if type(successes) is not list or len(successes) != 18:
        fail("{} does not contain all 18 selected successes".format(
            context))
    for ordinal, cell in enumerate(successes):
        _validate_worker_cell(
            cell, attempt, "{} success {}".format(context, ordinal), True,
            SELECTION_ROOTS)
        if cell["cell_ordinal"] != ordinal:
            fail("{} selected successes are not in frozen order".format(
                context))
    for lower_attempt, cell in enumerate(lower):
        _validate_worker_cell(
            cell, lower_attempt,
            "{} lower attempt {}".format(context, lower_attempt), False,
            SELECTION_ROOTS)
        selected_cell = successes[cell["cell_ordinal"]]
        if (cell["trace_sha256"] != selected_cell["trace_sha256"] or
                cell["attempted_candidates"] !=
                selected_cell["attempted_candidates"]):
            fail("{} lower-attempt witness does not match the frozen "
                 "selection trace".format(context))
    return value


def parse_derivation_stream(stdout: bytes, worker_binary_sha256: str) \
        -> Tuple[Mapping[str, Any], ...]:
    try:
        text = stdout.decode("utf-8")
    except UnicodeDecodeError:
        fail("repair derivation stream is not UTF-8")
    lines = text.splitlines()
    if len(lines) != len(K_VALUES) or any(not line for line in lines):
        fail("repair derivation stream has the wrong record count")
    records: List[Mapping[str, Any]] = []
    for line, K in zip(lines, K_VALUES):
        try:
            value = json.loads(line)
        except (TypeError, ValueError) as exc:
            fail("repair derivation K={} is invalid JSON: {}".format(K, exc))
        if line != canonical_json(value):
            fail("repair derivation K={} is not canonical JSON".format(K))
        records.append(validate_derivation_record(
            value, K, worker_binary_sha256))
    return tuple(records)


def derive_attempts(worker: Path, worker_binary_sha256: str,
                    executable_fd: Optional[int] = None) \
        -> Tuple[Tuple[Mapping[str, Any], ...], bytes, str]:
    commands = "".join(
        "D {} {}\n".format(K - 2, K) for K in K_VALUES) + "Q\n"
    command_bytes = commands.encode("ascii")
    stdout, stderr = _run_process(
        [str(worker), "--worker"], command_bytes, WORKER_TIMEOUT_SECONDS,
        MAX_WORKER_STDOUT_BYTES, "repair-worker derivation", executable_fd)
    if stderr:
        fail("repair-worker derivation wrote stderr")
    records = parse_derivation_stream(stdout, worker_binary_sha256)
    canonical = ("\n".join(canonical_json(record) for record in records) +
                 "\n").encode("utf-8")
    if stdout != canonical:
        fail("repair-worker derivation stream is not byte-canonical")
    return records, canonical, sha256_bytes(command_bytes)


@dataclass(frozen=True)
class Invocation:
    ordinal: int
    coordinate_ordinal: int
    K: int
    root_index: int
    schedule_index: int
    arm: str
    attempt: int
    timing_order: str
    timing_slot: int
    observation_index: int

    @property
    def root(self) -> str:
        return ROOTS[self.root_index]

    @property
    def schedule(self) -> str:
        return SCHEDULES[self.schedule_index]

    @property
    def cell_ordinal(self) -> int:
        return self.root_index * len(SCHEDULES) + self.schedule_index

    def argv(self, bench: Path) -> List[str]:
        candidate = self.arm == "candidate_two07_mix2"
        command = [
            str(bench), "precodefail", "--N", str(self.K),
            "--bb-list", str(BLOCK_BYTES), "--overhead", str(OVERHEAD),
            "--trials", "1", "--threads", "1", "--loss", "0.5",
            "--seed", self.root, "--schedule", self.schedule,
            "--heavy-family", "periodic", "--mix-count",
            "2" if candidate else "3", "--binary-dense-rows", "12",
            "--gf256-heavy-rows", "12", "--dense-anchors",
            "two07" if candidate else "disabled",
            "--paired-overhead-stream", "--full-payload-solve",
            "--cold-solve-wide-xor", "policy",
        ]
        if candidate:
            command.extend([
                "--exact-precode-attempt", str(self.attempt),
                "--exact-packet-attempt", str(self.attempt),
                "--construction-seed-basis", "production-profile",
            ])
        return command

    def identity(self, bench: Path) -> Mapping[str, Any]:
        command = self.argv(bench)
        return {
            "ordinal": self.ordinal,
            "coordinate_ordinal": self.coordinate_ordinal,
            "K": self.K,
            "root_index": self.root_index,
            "loss_root": self.root,
            "schedule": self.schedule,
            "cell_ordinal": self.cell_ordinal,
            "arm": self.arm,
            "construction_attempt": self.attempt,
            "timing_order": self.timing_order,
            "timing_slot": self.timing_slot,
            "observation_index": self.observation_index,
            "argv": command,
            "command_sha256": sha256_text(canonical_json(command)),
        }


def make_invocations(attempt_by_K: Mapping[int, int]) \
        -> Tuple[Invocation, ...]:
    if set(attempt_by_K) != set(K_VALUES):
        fail("attempt map does not cover the frozen K roster exactly")
    invocations: List[Invocation] = []
    coordinate_ordinal = 0
    for K in K_VALUES:
        attempt = _uint(attempt_by_K[K], 255, "attempt map K={}".format(K))
        for root_index in range(len(ROOTS)):
            for schedule_index in range(len(SCHEDULES)):
                timing_order = TIMING_ORDERS[coordinate_ordinal % 2]
                if timing_order == "ABBA":
                    arm_order = (ARMS[0], ARMS[1], ARMS[1], ARMS[0])
                else:
                    arm_order = (ARMS[1], ARMS[0], ARMS[0], ARMS[1])
                observation_count = {arm: 0 for arm in ARMS}
                for timing_slot, arm in enumerate(arm_order):
                    observation_index = observation_count[arm]
                    observation_count[arm] += 1
                    invocations.append(Invocation(
                        len(invocations), coordinate_ordinal, K, root_index,
                        schedule_index, arm,
                        attempt if arm == "candidate_two07_mix2" else 0,
                        timing_order, timing_slot, observation_index))
                if any(count != OBSERVATIONS_PER_ARM
                       for count in observation_count.values()):
                    fail("internal timing order has the wrong arm counts")
                coordinate_ordinal += 1
    if coordinate_ordinal != EXPECTED_CELL_COUNT:
        fail("internal benchmark coordinate roster has the wrong size")
    if len(invocations) != EXPECTED_RECORD_COUNT:
        fail("internal invocation roster has the wrong size")
    return tuple(invocations)


def _parse_metadata(line: str, invocation: Invocation,
                    source_commit: str) -> Mapping[str, str]:
    prefix = "# precodefail: "
    if not line.startswith(prefix):
        fail("bench cell {} omits precodefail metadata".format(
            invocation.ordinal))
    values: Dict[str, str] = {}
    for token in line[len(prefix):].split(" "):
        if token.count("=") != 1:
            fail("bench cell {} metadata is malformed".format(
                invocation.ordinal))
        key, value = token.split("=", 1)
        if not key or not value or key in values:
            fail("bench cell {} metadata is ambiguous".format(
                invocation.ordinal))
        values[key] = value
    if set(values) != METADATA_KEYS:
        fail("bench cell {} metadata fields differ from the protocol".format(
            invocation.ordinal))
    candidate = invocation.arm == "candidate_two07_mix2"
    expected = {
        "trials": "1", "threads": "1", "loss": "0.5",
        "source_hits_override": "0", "packet_peel_seed_xor": "0x0",
        "binary_dense_rows_override": "12",
        "gf256_heavy_rows_override": "12",
        "dense_anchor_layout": "two07" if candidate else "disabled",
        "odd_packet_peel_seed_xor": "0x0",
        "packet_row_seed_multiplier": "0x1",
        "packet_row_seed_avalanche": "0",
        "seed_block_bytes_override": "0", "overhead_stream": "paired",
        "full_payload_solve": "1", "schedule": invocation.schedule,
        "cold_solve_wide_xor": "policy",
        "exact_attempt_mode": "1" if candidate else "0",
        "exact_precode_attempt": str(invocation.attempt if candidate else 0),
        "exact_packet_attempt": str(invocation.attempt if candidate else 0),
        "construction_seed_basis": "production-profile",
        "seed_schedule_sha256": ZERO_SHA256,
        "source_git_commit": source_commit,
    }
    for key, expected_value in expected.items():
        if values[key] != expected_value:
            fail("bench cell {} metadata {} differs from the frozen input".
                 format(invocation.ordinal, key))
    try:
        seed = int(values["seed"], 0)
    except ValueError:
        fail("bench cell {} metadata seed is invalid".format(
            invocation.ordinal))
    if seed != int(invocation.root, 16):
        fail("bench cell {} metadata seed differs from its loss root".format(
            invocation.ordinal))
    return values


def _parse_histogram(value: str, trials: int, context: str) -> None:
    total = 0
    previous = -1
    if not value:
        fail("{} histogram is empty".format(context))
    for item in value.split("|"):
        parts = item.split(":")
        if len(parts) != 2:
            fail("{} histogram is malformed".format(context))
        key = _uint_text(parts[0], (1 << 32) - 1, context + " key")
        count = _uint_text(parts[1], trials, context + " count")
        if key <= previous or count == 0:
            fail("{} histogram is noncanonical".format(context))
        previous = key
        total += count
    if total != trials:
        fail("{} histogram count differs from trials".format(context))


def parse_bench_stdout(stdout: bytes, invocation: Invocation,
                       source_commit: str, bench_binary_sha256: str,
                       derivation: Mapping[str, Any], command_sha256: str,
                       source_receipt_sha256: str,
                       attempt_stream_sha256: str) -> Mapping[str, Any]:
    try:
        text = stdout.decode("utf-8")
    except UnicodeDecodeError:
        fail("bench cell {} stdout is not UTF-8".format(invocation.ordinal))
    lines = text.splitlines()
    if len(lines) != 3 or any(not line for line in lines):
        fail("bench cell {} must emit metadata, header, and one row".format(
            invocation.ordinal))
    if not text.endswith("\n") or "\r" in text:
        fail("bench cell {} stdout is not newline-canonical".format(
            invocation.ordinal))
    _parse_metadata(lines[0], invocation, source_commit)
    if tuple(lines[1].split(",")) != CSV_HEADER:
        fail("bench cell {} CSV header differs from the protocol".format(
            invocation.ordinal))
    try:
        parsed_rows = list(csv.reader(io.StringIO(lines[2])))
    except csv.Error as exc:
        fail("bench cell {} CSV is invalid: {}".format(
            invocation.ordinal, exc))
    if len(parsed_rows) != 1 or len(parsed_rows[0]) != len(CSV_HEADER):
        fail("bench cell {} CSV row has the wrong width".format(
            invocation.ordinal))
    if ",".join(parsed_rows[0]) != lines[2]:
        fail("bench cell {} CSV row is not canonical".format(
            invocation.ordinal))
    row = dict(zip(CSV_HEADER, parsed_rows[0]))
    context = "bench cell {}".format(invocation.ordinal)
    candidate = invocation.arm == "candidate_two07_mix2"

    K = _uint_text(row["N"], 64000, context + " K")
    bb = _uint_text(row["bb"], (1 << 31) - 1, context + " block bytes")
    overhead = _uint_text(row["overhead"], 1024, context + " overhead")
    trials = _uint_text(row["trials"], 1, context + " trials")
    if (K, bb, overhead, trials) != \
            (invocation.K, BLOCK_BYTES, OVERHEAD, 1):
        fail("{} coordinates differ from the frozen invocation".format(
            context))
    if row["heavy_family"] != "periodic" or \
            _uint_text(row["mix_count"], 3, context + " mix count") != \
            (2 if candidate else 3) or \
            _uint_text(row["binary_dense_rows"], 64,
                       context + " binary dense rows") != 12 or \
            _uint_text(row["gf256_heavy_rows"], 128,
                       context + " GF256 heavy rows") != 12 or \
            row["dense_identity_corner"] != "0":
        fail("{} realized structure differs from the frozen arm".format(
            context))
    staircase = _uint_text(row["staircase"], 4096,
                           context + " staircase")
    source_hits = _uint_text(row["source_hits"], 8,
                            context + " source hits")
    if staircase == 0 or source_hits == 0:
        fail("{} realized profile geometry is empty".format(context))

    success = _uint_text(row["success"], 1, context + " successes")
    rank_fail = _uint_text(row["rank_fail"], 1, context + " rank failures")
    error = _uint_text(row["error"], 1, context + " errors")
    if success + rank_fail + error != 1:
        fail("{} outcome counts do not sum to one".format(context))
    weak = success != 1
    expected_fail_rate = "1.00000000" if weak else "0.00000000"
    _decimal(row["fail_rate"], context + " failure rate", FAIL_RATE_TOKEN)
    if row["fail_rate"] != expected_fail_rate:
        fail("{} failure rate disagrees with its outcome".format(context))
    if row["first_rank_fail"] not in ("-1", "0"):
        fail("{} first-rank-failure receipt is noncanonical".format(context))
    first_rank_fail = int(row["first_rank_fail"])
    if first_rank_fail != (0 if rank_fail else -1):
        fail("{} first-rank-failure receipt is inconsistent".format(context))
    if row["failure_trials"] != ("0" if weak else ""):
        fail("{} failure trial receipt is inconsistent".format(context))
    if row["active_packet_peel_seed_xor"] != "0x0":
        fail("{} has an unexpected packet seed perturbation".format(context))

    inactivated = _integral_decimal(row["inact_mu"],
                                    context + " inactivation")
    if _uint_text(row["inact_max"], (1 << 32) - 1,
                  context + " maximum inactivation") != inactivated:
        fail("{} inactivation mean/max disagree for one trial".format(context))
    for key in ("binary_def_mu", "heavy_gain_mu"):
        _integral_decimal(row[key], context + " " + key)
    for key in ("binary_def_max", "heavy_gain_min", "heavy_shortfall"):
        _uint_text(row[key], (1 << 32) - 1, context + " " + key)
    _parse_histogram(row["binary_def_hist"], 1,
                     context + " binary deficiency")
    _parse_histogram(row["heavy_gain_hist"], 1,
                     context + " heavy gain")
    timing = {}
    for key in ("solve_ms_mu", "build_ms_mu", "peel_ms_mu",
                "project_ms_mu", "residual_ms_mu", "backsub_ms_mu"):
        _decimal(row[key], context + " " + key)
        timing[key] = row[key]
    block_xors = _integral_decimal(row["block_xors_mu"],
                                   context + " block XORs")
    block_muladds = _integral_decimal(row["block_muladds_mu"],
                                      context + " GF256 muladds")

    precode_attempt = _uint_text(row["precode_attempt"], 255,
                                 context + " precode attempt")
    packet_attempt = _uint_text(row["packet_attempt"], 255,
                                context + " packet attempt")
    if precode_attempt != packet_attempt:
        fail("{} uses different precode and packet attempts".format(context))
    if candidate:
        if row["seed_attempt"] != "" or row["attempt_mode"] != "exact" or \
                precode_attempt != invocation.attempt:
            fail("{} did not execute the exact repaired attempt".format(
                context))
    else:
        selected = _uint_text(row["seed_attempt"], 255,
                              context + " selected attempt")
        if row["attempt_mode"] != "selected" or \
                precode_attempt != selected:
            fail("{} did not execute current production selection".format(
                context))
    if row["construction_seed_basis"] != "production-profile" or \
            row["seed_schedule_sha256"] != ZERO_SHA256:
        fail("{} did not execute the production-profile seed basis".format(
            context))
    _hex_value(row["effective_precode_seed"], HEX64_TOKEN,
               context + " effective precode seed")
    _hex_value(row["effective_packet_seed"], HEX32_TOKEN,
               context + " effective packet seed")
    if candidate and (row["effective_precode_seed"] !=
                      derivation["effective_precode_seed"] or
                      row["effective_packet_seed"] !=
                      derivation["effective_packet_seed"]):
        fail("{} effective seeds differ from the v7-derived attempt".format(
            context))

    _sha(attempt_stream_sha256, context + " attempt stream hash")
    extra_entries = ((K + staircase + 12 + 1) // 2) - 2
    if extra_entries < 0:
        fail("{} dense-basis capacity formula underflowed".format(context))
    record: Dict[str, Any] = {
        "schema": RESULT_SCHEMA,
        "ordinal": invocation.ordinal,
        "coordinate_ordinal": invocation.coordinate_ordinal,
        "arm": invocation.arm,
        "K": K,
        "block_bytes": bb,
        "loss_ppm": LOSS_PPM,
        "overhead": overhead,
        "root_index": invocation.root_index,
        "loss_root": invocation.root,
        "schedule": invocation.schedule,
        "cell_ordinal": invocation.cell_ordinal,
        "timing_order": invocation.timing_order,
        "timing_slot": invocation.timing_slot,
        "observation_index": invocation.observation_index,
        "attempt_selection_stream_sha256": attempt_stream_sha256,
        "attempt_selection_cell_receipts_used_as_promotion_evidence": False,
        "benchmark_loss_trace_hash_recorded": False,
        "full_payload_byte_recovery_verified": not weak,
        "candidate_profile_sha256": candidate_profile_sha256(),
        "construction_attempt": precode_attempt,
        "attempt_mode": row["attempt_mode"],
        "effective_precode_seed": row["effective_precode_seed"],
        "effective_packet_seed": row["effective_packet_seed"],
        "staircase": staircase,
        "binary_dense_rows": 12,
        "gf256_heavy_rows": 12,
        "source_hits": source_hits,
        "dense_anchor_layout": "two07" if candidate else "disabled",
        "mix_count": 2 if candidate else 3,
        "success": success,
        "rank_fail": rank_fail,
        "error": error,
        "weak": weak,
        "inactivated_columns": inactivated,
        "block_xors": block_xors,
        "gf256_muladds": block_muladds,
        "solve_ms": timing["solve_ms_mu"],
        "build_ms": timing["build_ms_mu"],
        "peel_ms": timing["peel_ms_mu"],
        "project_ms": timing["project_ms_mu"],
        "residual_ms": timing["residual_ms_mu"],
        "backsub_ms": timing["backsub_ms_mu"],
        "extra_dense_basis_capacity_entries":
            extra_entries if candidate else 0,
        "extra_dense_basis_capacity_bytes":
            extra_entries * 4 if candidate else 0,
        "command_sha256": command_sha256,
        "bench_stdout_sha256": sha256_bytes(stdout),
        "bench_binary_sha256": bench_binary_sha256,
        "bench_source_git_commit": source_commit,
        "source_receipt_sha256": source_receipt_sha256,
    }
    return record


def run_bench_cell(bench: Path, invocation: Invocation,
                   source_commit: str, bench_binary_sha256: str,
                   derivation: Mapping[str, Any],
                   source_receipt_sha256: str,
                   attempt_stream_sha256: str,
                   executable_fd: Optional[int] = None) \
        -> Mapping[str, Any]:
    command = invocation.argv(bench)
    command_sha = sha256_text(canonical_json(command))
    stdout, stderr = _run_process(
        command, None, BENCH_TIMEOUT_SECONDS, MAX_BENCH_STDOUT_BYTES,
        "bench cell {}".format(invocation.ordinal), executable_fd)
    if stderr:
        fail("bench cell {} wrote stderr".format(invocation.ordinal))
    return parse_bench_stdout(
        stdout, invocation, source_commit, bench_binary_sha256, derivation,
        command_sha, source_receipt_sha256, attempt_stream_sha256)


def extra_dense_basis_capacity_entries(K: int, staircase: int) -> int:
    _uint(K, 64000, "complexity K")
    _uint(staircase, 4096, "complexity staircase")
    if K < 2 or staircase == 0:
        fail("complexity coordinates are outside the production domain")
    return ((K + staircase + 12 + 1) // 2) - 2


def _ratio(numerator: Decimal, denominator: Decimal,
           context: str) -> Decimal:
    if denominator <= 0:
        fail("{} control aggregate is not positive".format(context))
    return numerator / denominator


def _ratio_text(value: Decimal) -> str:
    return format(value, ".12f")


DETERMINISTIC_OBSERVATION_FIELDS = (
    "arm", "K", "block_bytes", "loss_ppm", "overhead", "root_index",
    "loss_root", "schedule", "cell_ordinal",
    "attempt_selection_stream_sha256",
    "attempt_selection_cell_receipts_used_as_promotion_evidence",
    "benchmark_loss_trace_hash_recorded", "full_payload_byte_recovery_verified",
    "candidate_profile_sha256", "construction_attempt", "attempt_mode",
    "effective_precode_seed", "effective_packet_seed", "staircase",
    "binary_dense_rows", "gf256_heavy_rows", "source_hits",
    "dense_anchor_layout", "mix_count", "success", "rank_fail", "error",
    "weak", "inactivated_columns", "block_xors", "gf256_muladds",
    "extra_dense_basis_capacity_entries",
    "extra_dense_basis_capacity_bytes", "command_sha256",
    "bench_binary_sha256", "bench_source_git_commit",
    "source_receipt_sha256",
)


def _deterministic_projection(record: Mapping[str, Any]) \
        -> Mapping[str, Any]:
    try:
        return {field: record[field]
                for field in DETERMINISTIC_OBSERVATION_FIELDS}
    except KeyError as exc:
        fail("result record omits deterministic field {}".format(exc))


def adjudicate(records: Sequence[Mapping[str, Any]],
               expected_attempt_stream_sha256: Optional[str] = None) \
        -> Mapping[str, Any]:
    if len(records) != EXPECTED_RECORD_COUNT:
        fail("result stream has the wrong record count")
    if expected_attempt_stream_sha256 is not None:
        _sha(expected_attempt_stream_sha256, "expected attempt stream hash")
    keyed: Dict[Tuple[int, int, str, str, int], Mapping[str, Any]] = {}
    coordinate_ordinal = 0
    ordinal = 0
    observed_attempt_stream_sha: Optional[str] = None
    for K in K_VALUES:
        for root_index, root in enumerate(ROOTS):
            for schedule_index, schedule in enumerate(SCHEDULES):
                timing_order = TIMING_ORDERS[coordinate_ordinal % 2]
                arm_order = ((ARMS[0], ARMS[1], ARMS[1], ARMS[0])
                             if timing_order == "ABBA" else
                             (ARMS[1], ARMS[0], ARMS[0], ARMS[1]))
                observation_count = {arm: 0 for arm in ARMS}
                for timing_slot, arm in enumerate(arm_order):
                    record = records[ordinal]
                    observation_index = observation_count[arm]
                    observation_count[arm] += 1
                    if (type(record) is not dict or
                            record.get("schema") != RESULT_SCHEMA or
                            record.get("ordinal") != ordinal or
                            record.get("coordinate_ordinal") !=
                            coordinate_ordinal or record.get("K") != K or
                            record.get("root_index") != root_index or
                            record.get("loss_root") != root or
                            record.get("schedule") != schedule or
                            record.get("cell_ordinal") !=
                            root_index * len(SCHEDULES) + schedule_index or
                            record.get("arm") != arm or
                            record.get("timing_order") != timing_order or
                            record.get("timing_slot") != timing_slot or
                            record.get("observation_index") !=
                            observation_index):
                        fail("result stream is not in canonical v7 ABBA/BAAB "
                             "roster order")
                    stream_sha = _sha(
                        record.get("attempt_selection_stream_sha256"),
                        "result attempt stream hash")
                    if observed_attempt_stream_sha is None:
                        observed_attempt_stream_sha = stream_sha
                    if stream_sha != observed_attempt_stream_sha or (
                            expected_attempt_stream_sha256 is not None and
                            stream_sha != expected_attempt_stream_sha256):
                        fail("result records do not bind one frozen attempt "
                             "stream")
                    weak = record.get("weak")
                    recovery_verified = record.get(
                        "full_payload_byte_recovery_verified")
                    if (type(weak) is not bool or record.get(
                            "attempt_selection_cell_receipts_used_as_"
                            "promotion_evidence") is not False or
                            record.get("benchmark_loss_trace_hash_recorded")
                            is not False or type(recovery_verified) is not
                            bool or recovery_verified != (not weak)):
                        fail("result recovery-evidence semantics are "
                             "inconsistent")
                    key = (K, root_index, schedule, arm, observation_index)
                    if key in keyed:
                        fail("result stream repeats an arm observation")
                    keyed[key] = record
                    ordinal += 1
                coordinate_ordinal += 1
    if ordinal != EXPECTED_RECORD_COUNT or len(keyed) != EXPECTED_RECORD_COUNT:
        fail("result stream contains observations outside the frozen roster")

    candidate_weak = 0
    control_weak = 0
    candidate_weak_coordinates = 0
    control_weak_coordinates = 0
    common = 0
    common_timing_pairs = 0
    old_xor = new_xor = 0
    old_muladd = new_muladd = 0
    old_inact = new_inact = 0
    old_solve = new_solve = Decimal(0)
    old_build = new_build = Decimal(0)
    maximum_entries = -1
    maximum_bytes = -1
    for K in K_VALUES:
        for root_index in range(len(ROOTS)):
            for schedule in SCHEDULES:
                coordinate = (K, root_index, schedule)
                control = [keyed[coordinate + (ARMS[0], index)]
                           for index in range(OBSERVATIONS_PER_ARM)]
                candidate = [keyed[coordinate + (ARMS[1], index)]
                             for index in range(OBSERVATIONS_PER_ARM)]
                if (_deterministic_projection(control[0]) !=
                        _deterministic_projection(control[1]) or
                        _deterministic_projection(candidate[0]) !=
                        _deterministic_projection(candidate[1])):
                    fail("repeated observations disagree on deterministic "
                         "work or recovery")
                if (control[0]["staircase"] != candidate[0]["staircase"] or
                        control[0]["source_hits"] !=
                        candidate[0]["source_hits"]):
                    fail("paired arms do not share production geometry")
                expected_entries = extra_dense_basis_capacity_entries(
                    K, candidate[0]["staircase"])
                candidate_entries = candidate[0][
                    "extra_dense_basis_capacity_entries"]
                candidate_bytes = candidate[0][
                    "extra_dense_basis_capacity_bytes"]
                if (control[0]["extra_dense_basis_capacity_entries"] != 0 or
                        control[0]["extra_dense_basis_capacity_bytes"] != 0 or
                        candidate_entries != expected_entries or
                        candidate_bytes != expected_entries * 4):
                    fail("dense-basis extra capacity differs from the formula")
                maximum_entries = max(maximum_entries, expected_entries)
                maximum_bytes = max(maximum_bytes, expected_entries * 4)
                candidate_cell_weak = bool(candidate[0]["weak"])
                control_cell_weak = bool(control[0]["weak"])
                candidate_weak += sum(1 for item in candidate if item["weak"])
                control_weak += sum(1 for item in control if item["weak"])
                candidate_weak_coordinates += 1 if candidate_cell_weak else 0
                control_weak_coordinates += 1 if control_cell_weak else 0
                if candidate_cell_weak or control_cell_weak:
                    continue
                common += 1
                old_xor += control[0]["block_xors"]
                new_xor += candidate[0]["block_xors"]
                old_muladd += control[0]["gf256_muladds"]
                new_muladd += candidate[0]["gf256_muladds"]
                old_inact += control[0]["inactivated_columns"]
                new_inact += candidate[0]["inactivated_columns"]
                for observation_index in range(OBSERVATIONS_PER_ARM):
                    old_solve += Decimal(control[observation_index]["solve_ms"])
                    new_solve += Decimal(
                        candidate[observation_index]["solve_ms"])
                    old_build += Decimal(control[observation_index]["build_ms"])
                    new_build += Decimal(
                        candidate[observation_index]["build_ms"])
                    common_timing_pairs += 1

    K64000 = [record for record in records
              if record["K"] == 64000 and
              record["arm"] == "candidate_two07_mix2"]
    if not K64000 or any(record["staircase"] != 346 for record in K64000):
        fail("realized K64000 staircase is not the preregistered S346")
    if maximum_entries != 32177 or maximum_bytes != 128708:
        fail("maximum dense-basis capacity is not 32177 entries/128708 bytes")

    xor_ratio = _ratio(Decimal(new_xor), Decimal(old_xor), "block XOR")
    muladd_ratio = _ratio(
        Decimal(new_muladd), Decimal(old_muladd), "GF256 muladd")
    inact_ratio = _ratio(
        Decimal(new_inact), Decimal(old_inact), "inactivation")
    solve_ratio = _ratio(new_solve, old_solve, "solve time")
    build_ratio = _ratio(new_build, old_build, "build time")
    gates = {
        "candidate_benchmark_weak_observations_zero": candidate_weak == 0,
        "candidate_benchmark_weak_observations_not_above_control":
            candidate_weak <= control_weak,
        "repeated_deterministic_work_equal": True,
        "block_xor_ratio_at_most_0.9829453304":
            xor_ratio <= XOR_RATIO_MAX,
        "gf256_muladd_ratio_at_most_1":
            muladd_ratio <= GF256_MULADD_RATIO_MAX,
        "inactivation_ratio_at_most_1":
            inact_ratio <= INACTIVATION_RATIO_MAX,
        "aggregate_b2_solve_time_ratio_at_most_1":
            solve_ratio <= SOLVE_TIME_RATIO_MAX,
        "aggregate_b2_build_time_ratio_at_most_1.05":
            build_ratio <= BUILD_TIME_RATIO_MAX,
    }
    return {
        "disposition": "PASS" if all(gates.values()) else "REJECT",
        "gates": gates,
        "candidate_weak_observations": candidate_weak,
        "control_weak_observations": control_weak,
        "candidate_weak_coordinates": candidate_weak_coordinates,
        "control_weak_coordinates": control_weak_coordinates,
        "common_success_cells": common,
        "common_success_timing_pairs": common_timing_pairs,
        "aggregates": {
            "control": {
                "block_xors": old_xor,
                "gf256_muladds": old_muladd,
                "inactivated_columns": old_inact,
                "solve_ms": format(old_solve, "f"),
                "build_ms": format(old_build, "f"),
            },
            "candidate": {
                "block_xors": new_xor,
                "gf256_muladds": new_muladd,
                "inactivated_columns": new_inact,
                "solve_ms": format(new_solve, "f"),
                "build_ms": format(new_build, "f"),
            },
        },
        "ratios": {
            "block_xors": _ratio_text(xor_ratio),
            "gf256_muladds": _ratio_text(muladd_ratio),
            "inactivated_columns": _ratio_text(inact_ratio),
            "solve_time": _ratio_text(solve_ratio),
            "build_time": _ratio_text(build_ratio),
        },
        "complexity_receipt": {
            "asymptotic_construction_control": "O(K)",
            "asymptotic_construction_candidate": "O(K)",
            "asymptotic_storage_control": "O(K)",
            "asymptotic_storage_candidate": "O(K)",
            "allocation_count_delta": 0,
            "runtime_online_retries": 0,
            "production_map_lookups_per_construction": 1,
            "new_equations": 0,
            "extra_dense_basis_capacity_formula":
                "ceil((K+S+12)/2)-2 uint32 entries",
            "maximum_realized_K": 64000,
            "maximum_realized_staircase": 346,
            "maximum_extra_entries": maximum_entries,
            "maximum_extra_bytes": maximum_bytes,
            "scope":
                "static promotion complexity claim plus exact realized "
                "capacity arithmetic; allocation claim requires source review",
        },
    }


def _source_receipt() -> Mapping[str, Any]:
    hash_paths = list(SOURCE_PATHS)
    for item in hash_paths:
        path = REPO_ROOT / item
        if not path.is_file():
            fail("required source is missing: {}".format(item))
    try:
        head = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=str(REPO_ROOT),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True,
            timeout=30).stdout.decode("ascii").strip()
        status_output = subprocess.run(
            ["git", "status", "--porcelain=v1", "--untracked-files=all",
             "--"] + list(SOURCE_STATUS_PATHS), cwd=str(REPO_ROOT),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=True, timeout=30).stdout
    except (OSError, subprocess.SubprocessError, UnicodeDecodeError) as exc:
        fail("cannot attest source checkout: {}".format(exc))
    if not COMMIT_TOKEN.fullmatch(head):
        fail("source HEAD is not a full lowercase commit")
    if status_output:
        fail("linked short-screen sources are not clean at source HEAD")
    hashes = {
        item: sha256_bytes(_stable_file_bytes(REPO_ROOT / item))
        for item in hash_paths
    }
    receipt: Dict[str, Any] = {
        "source_git_commit": head,
        "tracked_and_untracked_linked_sources_clean": True,
        "clean_status_scope": list(SOURCE_STATUS_PATHS),
        "source_files": hashes,
    }
    receipt["source_receipt_sha256"] = self_hash(
        receipt, "source_receipt_sha256")
    return receipt


def _write_exclusive(path: Path, data: bytes) -> None:
    try:
        descriptor = os.open(
            str(path), os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    except OSError as exc:
        fail("cannot create {}: {}".format(path, exc))
    try:
        offset = 0
        while offset < len(data):
            written = os.write(descriptor, data[offset:])
            if written <= 0:
                fail("short write while creating {}".format(path))
            offset += written
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _canonical_json_file(value: Mapping[str, Any]) -> bytes:
    return (canonical_json(value) + "\n").encode("utf-8")


def _publish_directory_noreplace(source: Path, destination: Path) -> None:
    if _RENAMEAT2 is None:
        fail("atomic no-replace publication requires renameat2")
    source_bytes = os.fsencode(str(source))
    destination_bytes = os.fsencode(str(destination))
    ctypes.set_errno(0)
    result = _RENAMEAT2(
        _AT_FDCWD, source_bytes, _AT_FDCWD, destination_bytes,
        _RENAME_NOREPLACE)
    if result == 0:
        return
    error_number = ctypes.get_errno()
    if error_number == errno.EEXIST:
        fail("output directory appeared before publication")
    fail("cannot publish output directory atomically: {}".format(
        os.strerror(error_number)))


def run_screen(bench: Path, worker: Path, output_dir: Path) \
        -> Mapping[str, Any]:
    started = time.monotonic()
    if output_dir.exists():
        fail("output directory already exists")
    parent = output_dir.parent
    if not parent.is_dir():
        fail("output parent directory does not exist")
    with _open_binary(bench, "wirehair_v2_bench") as pinned_bench, \
            _open_binary(worker, "repair worker") as pinned_worker:
        return _run_screen_pinned(
            pinned_bench, pinned_worker, output_dir, parent, started)


def _run_screen_pinned(
        pinned_bench: PinnedExecutable, pinned_worker: PinnedExecutable,
        output_dir: Path, parent: Path, started: float) \
        -> Mapping[str, Any]:
    bench_receipt = pinned_bench.receipt()
    worker_receipt = pinned_worker.receipt()
    bench = pinned_bench.path
    worker = pinned_worker.path
    source = _source_receipt()
    source_commit = source["source_git_commit"]
    source_receipt_sha = source["source_receipt_sha256"]
    description = read_worker_description(
        worker, worker_receipt["sha256"], source_commit,
        pinned_worker.descriptor)
    derivations, attempt_bytes, worker_command_sha = derive_attempts(
        worker, worker_receipt["sha256"], pinned_worker.descriptor)
    attempt_stream_sha = sha256_bytes(attempt_bytes)
    attempt_by_K = {
        record["K"]: record["selected_attempt"] for record in derivations
    }
    derivation_by_K = {record["K"]: record for record in derivations}
    invocations = make_invocations(attempt_by_K)
    contract = contract_description()
    input_receipt: Dict[str, Any] = {
        "schema": INPUT_SCHEMA,
        "contract": contract,
        "contract_sha256": contract["contract_sha256"],
        "source_receipt": source,
        "source_receipt_sha256": source_receipt_sha,
        "bench_binary": bench_receipt,
        "repair_worker_binary": worker_receipt,
        "repair_worker_description": description,
        "repair_worker_description_sha256":
            sha256_text(canonical_json(description)),
        "attempt_selection_worker_command_sha256": worker_command_sha,
        "attempt_selection_stream_schema": ATTEMPT_STREAM_SCHEMA,
        "attempt_selection_stream_sha256": attempt_stream_sha,
        "attempt_selection_record_count": len(derivations),
        "selected_attempts": [attempt_by_K[K] for K in K_VALUES],
        "worker_validation_commands_issued": 0,
        "worker_final_validation_stream_present": False,
        "invocations": [item.identity(bench) for item in invocations],
    }
    input_receipt["input_sha256"] = self_hash(input_receipt, "input_sha256")
    input_bytes = _canonical_json_file(input_receipt)

    temporary = Path(tempfile.mkdtemp(
        prefix=".{}-".format(output_dir.name), dir=str(parent)))
    published = False
    try:
        _write_exclusive(temporary / ATTEMPT_NAME, attempt_bytes)
        _write_exclusive(temporary / INPUT_NAME, input_bytes)
        records: List[Mapping[str, Any]] = []
        for invocation in invocations:
            if time.monotonic() - started > TOTAL_WALL_SECONDS:
                fail("promotion short-screen hard wall expired")
            records.append(run_bench_cell(
                bench, invocation, source_commit, bench_receipt["sha256"],
                derivation_by_K[invocation.K], source_receipt_sha,
                attempt_stream_sha, pinned_bench.descriptor))
        result_bytes = ("\n".join(canonical_json(record)
                                  for record in records) + "\n").encode(
                                      "utf-8")
        verdict = adjudicate(records, attempt_stream_sha)
        if pinned_bench.receipt() != bench_receipt or \
                pinned_worker.receipt() != worker_receipt:
            fail("a measured binary changed during the short screen")
        if _source_receipt() != source:
            fail("linked source provenance changed during the short screen")
        summary: Dict[str, Any] = {
            "schema": SUMMARY_SCHEMA,
            "contract_sha256": contract["contract_sha256"],
            "input_sha256": input_receipt["input_sha256"],
            "input_file_sha256": sha256_bytes(input_bytes),
            "attempt_selection_stream_schema": ATTEMPT_STREAM_SCHEMA,
            "attempt_selection_stream_sha256": attempt_stream_sha,
            "attempt_selection_record_count": len(derivations),
            "worker_validation_commands_issued": 0,
            "worker_final_validation_stream_present": False,
            "result_stream_sha256": sha256_bytes(result_bytes),
            "record_count": len(records),
            "invocation_count": len(invocations),
            "bench_binary_sha256": bench_receipt["sha256"],
            "repair_worker_binary_sha256": worker_receipt["sha256"],
            "source_git_commit": source_commit,
            "source_receipt": source,
            "source_receipt_sha256": source_receipt_sha,
            "candidate_profile_sha256": candidate_profile_sha256(),
            "architecture_selection_performed": False,
            "offline_attempt_derivation_performed": True,
            "official_scope": "v7 production-basis bounded promotion screen",
            "all_K_recovery_claimed": False,
            "wirehair1_timing_claimed": False,
            **verdict,
        }
        summary["summary_sha256"] = self_hash(summary, "summary_sha256")
        summary_bytes = _canonical_json_file(summary)
        _write_exclusive(temporary / RESULT_NAME, result_bytes)
        _write_exclusive(temporary / SUMMARY_NAME, summary_bytes)
        _publish_directory_noreplace(temporary, output_dir)
        published = True
        return summary
    finally:
        if not published:
            shutil.rmtree(temporary, ignore_errors=True)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="action", required=True)
    subparsers.add_parser("describe", help="print the frozen contract")
    run = subparsers.add_parser("run", help="run and publish one confirmation")
    run.add_argument("--bench", required=True, type=Path)
    run.add_argument("--repair-worker", required=True, type=Path)
    run.add_argument("--output-dir", required=True, type=Path)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _build_parser().parse_args(argv)
    try:
        if args.action == "describe":
            print(canonical_json(contract_description()))
            return 0
        summary = run_screen(args.bench, args.repair_worker, args.output_dir)
        print(canonical_json(summary))
        return 0 if summary["disposition"] == "PASS" else 2
    except ShortScreenError as exc:
        print("mix2 promotion short screen: {}".format(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
