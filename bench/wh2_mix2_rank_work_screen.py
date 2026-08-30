#!/usr/bin/env python3
"""Launch and reduce the fixed pure-GF256 mix2 rank/work screen.

The ``manifest`` command is a read-only dry run.  The ``run`` command is the
only execution seam: it verifies one explicit benchmark artifact, runs the 96
closed precodefail commands sequentially without a shell, validates every
one-trial row, reduces the paired A/B/C/D experiment, and publishes two
no-clobber canonical artifacts.  Wall-clock fields are retained as diagnostic
raw evidence but are absent from the scientific decision projection.  A run
must present the exact manifest SHA-256 emitted by its preregistered dry run.
"""

import sys
sys.dont_write_bytecode = True

import argparse
import csv
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
import hashlib
import json
import math
import os
from pathlib import Path
import re
import selectors
import signal
import stat
import struct
import subprocess
import tempfile
import time
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple


MANIFEST_SCHEMA = "wirehair.wh2.mix2-rank-work-manifest.v1"
RECORD_SCHEMA = "wirehair.wh2.mix2-rank-work-record.v1"
DECISION_SCHEMA = "wirehair.wh2.mix2-rank-work-decision.v1"
DECISION_INPUT_SCHEMA = "wirehair.wh2.mix2-rank-work-decision-input.v1"
SUMMARY_SCHEMA = "wirehair.wh2.mix2-rank-work-summary.v1"
RESULT_NAME = "mix2-rank-work-results.jsonl"
SUMMARY_NAME = "mix2-rank-work-summary.json"

EXPECTED_INVOCATIONS = 96
ROWS_PER_INVOCATION = 64
EXPECTED_RECORDS = EXPECTED_INVOCATIONS * ROWS_PER_INVOCATION
MAX_STDOUT_BYTES = 512 * 1024
MAX_STDERR_BYTES = 64 * 1024
MAX_BINARY_BYTES = 512 * 1024 * 1024
MAX_TIMEOUT_SECONDS = 3600.0
MASK32 = (1 << 32) - 1
MASK64 = (1 << 64) - 1
TERMINATION_SIGNALS = (signal.SIGHUP, signal.SIGTERM)
CHILD_SPAWN_SIGNALS = (signal.SIGHUP, signal.SIGINT, signal.SIGTERM)
CHILD_ENVIRONMENT = {
    "LANG": "C",
    "LC_ALL": "C",
    "TZ": "UTC",
}

K_VALUES = (10000, 10001, 16384, 20000, 20001, 32768, 49152, 64000)
STAIRCASE_BY_K = {
    10000: 86, 10001: 86, 16384: 114, 20000: 134,
    20001: 134, 32768: 194, 49152: 378, 64000: 346,
}
OVERHEADS = (0, 1, 2, 4)
MIX_COUNTS = (3, 2)
LAYOUTS = ("two07", "four0369")
ATTEMPTS = (0, 1, 2)
ROOTS = (
    "0x000000005eedf411",
    "0xd1b54a32d192ed03",
    "0x94d049bb133111eb",
    "0x8538ecb5bd456ea3",
)
SCHEDULES = (
    ("iid", "0.1", 100000),
    ("burst", "0.5", 500000),
    ("adversarial", "0.5", 500000),
    ("repair-only", "0.5", 500000),
)

RAW_SEED_BASIS = "uniform-raw-v1"
RAW_SEED_SCHEDULE_SHA256 = (
    "90a98a3db207852dabdf5fb27573ef48b"
    "ce52e0228cee4e291d96fa44ed509a7"
)
RAW_PRECODE_BASE_SEED = 0x487468302aad7105
RAW_PACKET_BASE_SEED = 0x4ec72102
RAW_PRECODE_ATTEMPT_STRIDE = 0x9e3779b97f4a7c15
RAW_PACKET_ATTEMPT_STRIDE = 0x9e3779b9
REQUIRED_TEST_HOOK_MARKERS = (
    "--binary-dense-rows",
    "--gf256-heavy-rows",
    "--dense-anchors",
    "--exact-precode-attempt",
    "--exact-packet-attempt",
    "--construction-seed-basis",
    "uniform-raw-v1",
)

ARM_BY_LAYOUT_MIX = {
    ("two07", 3): "A",
    ("two07", 2): "B",
    ("four0369", 3): "C",
    ("four0369", 2): "D",
}
CSV_HEADER = (
    "N", "bb", "heavy_family", "mix_count", "staircase",
    "binary_dense_rows", "gf256_heavy_rows", "source_hits",
    "dense_identity_corner", "overhead", "trials", "success",
    "rank_fail", "error", "fail_rate", "inact_mu", "inact_max",
    "binary_def_mu", "binary_def_max", "heavy_gain_mu",
    "heavy_gain_min", "heavy_shortfall", "solve_ms_mu", "build_ms_mu",
    "peel_ms_mu", "project_ms_mu", "residual_ms_mu", "backsub_ms_mu",
    "seed_attempt", "block_xors_mu", "block_muladds_mu",
    "first_rank_fail", "binary_def_hist", "heavy_gain_hist",
    "failure_trials", "active_packet_peel_seed_xor", "precode_attempt",
    "packet_attempt", "attempt_mode", "construction_seed_basis",
    "seed_schedule_sha256", "effective_precode_seed",
    "effective_packet_seed",
)
TIMING_FIELDS = (
    "solve_ms_mu", "build_ms_mu", "peel_ms_mu", "project_ms_mu",
    "residual_ms_mu", "backsub_ms_mu",
)
RECORD_FIELDS = frozenset((
    "active_packet_peel_seed_xor", "arm", "argv_sha256", "attempt_mode",
    "binary_build_id", "binary_deficiency", "binary_dense_rows",
    "binary_sha256", "block_bytes", "block_xors", "construction_attempt",
    "construction_seed_basis", "controller_sha256", "csv_row_sha256",
    "dense_anchor_layout", "dense_identity_corner", "diagnostic_timings",
    "effective_packet_seed", "effective_precode_seed", "error",
    "fail_rate", "failure_trials", "first_rank_fail", "gf256_heavy_rows",
    "gf256_muladds", "gf256_rank_gain", "heavy_family",
    "heavy_shortfall", "inactivated_columns", "invocation_ordinal", "K",
    "legacy_seed_attempt", "loss_ppm", "loss_root", "metadata_sha256",
    "mix_count", "overhead", "packet_attempt", "precode_attempt",
    "rank_fail", "record_ordinal", "root_index", "row_ordinal", "schedule", "schema",
    "seed_schedule_sha256", "source_commit", "source_hits", "staircase",
    "stdout_sha256", "success", "trials",
))
METADATA_KEYS = (
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
)

HEX40 = re.compile(r"[0-9a-f]{40}\Z")
HEX64 = re.compile(r"[0-9a-f]{64}\Z")
BUILD_ID = re.compile(r"[0-9a-f]{8,128}\Z")
UINT_TEXT = re.compile(r"(?:0|[1-9][0-9]*)\Z")
THREE_DECIMAL = re.compile(r"(?:0|[1-9][0-9]*)\.[0-9]{3}\Z")


class ScreenError(RuntimeError):
    """The closed screen cannot be executed or reduced safely."""


def fail(message: str) -> None:
    raise ScreenError(message)


def _termination_handler(signum: int, _frame: Any) -> None:
    raise ScreenError("campaign received {}".format(signal.Signals(signum).name))


def _install_termination_handlers() -> List[Tuple[int, Any]]:
    previous: List[Tuple[int, Any]] = []
    try:
        for signum in TERMINATION_SIGNALS:
            previous.append((signum, signal.getsignal(signum)))
            signal.signal(signum, _termination_handler)
    except BaseException:
        _restore_termination_handlers(previous)
        raise
    return previous


def _restore_termination_handlers(previous: Sequence[Tuple[int, Any]]) -> None:
    for signum, handler in reversed(previous):
        signal.signal(signum, handler)


def canonical_json(value: Any) -> str:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_json(value: Any) -> str:
    return sha256_bytes(canonical_json(value).encode("ascii"))


def _require_sha256(value: str, name: str) -> str:
    if not isinstance(value, str) or HEX64.fullmatch(value) is None:
        fail("{} must be exactly 64 lowercase hexadecimal digits".format(name))
    return value


def _require_commit(value: str) -> str:
    if (not isinstance(value, str) or HEX40.fullmatch(value) is None or
            value == "0" * 40):
        fail("source commit must be one nonzero lowercase 40-hex identity")
    return value


def _require_build_id(value: str) -> str:
    if (not isinstance(value, str) or BUILD_ID.fullmatch(value) is None or
            len(value) % 2 or set(value) == {"0"}):
        fail("ELF build ID must be nonzero lowercase even-length hexadecimal")
    return value


def _effective_precode_seed(attempt: int) -> str:
    value = (RAW_PRECODE_BASE_SEED +
             attempt * RAW_PRECODE_ATTEMPT_STRIDE) & MASK64
    return "0x{:016x}".format(value)


def _effective_packet_seed(attempt: int) -> str:
    value = (RAW_PACKET_BASE_SEED +
             attempt * RAW_PACKET_ATTEMPT_STRIDE) & MASK32
    return "0x{:08x}".format(value)


@dataclass(frozen=True)
class Invocation:
    ordinal: int
    layout: str
    attempt: int
    root_index: int
    root: str
    schedule: str
    loss_arg: str
    loss_ppm: int

    def identity(self) -> Mapping[str, Any]:
        return {
            "construction_attempt": self.attempt,
            "dense_anchor_layout": self.layout,
            "loss_ppm": self.loss_ppm,
            "loss_root": self.root,
            "root_index": self.root_index,
            "schedule": self.schedule,
        }

    def argv(self, benchmark: Path) -> List[str]:
        return [
            str(benchmark), "precodefail",
            "--N", ",".join(str(value) for value in K_VALUES),
            "--bb-list", "1280",
            "--overhead", ",".join(str(value) for value in OVERHEADS),
            "--trials", "1", "--threads", "1",
            "--loss", self.loss_arg,
            "--seed", self.root,
            "--schedule", self.schedule,
            "--heavy-family", "periodic",
            "--mix-count", "3,2",
            "--binary-dense-rows", "12",
            "--gf256-heavy-rows", "12",
            "--dense-anchors", self.layout,
            "--paired-overhead-stream",
            "--cold-solve-wide-xor", "policy",
            "--exact-precode-attempt", str(self.attempt),
            "--exact-packet-attempt", str(self.attempt),
            "--construction-seed-basis", RAW_SEED_BASIS,
        ]


@dataclass(frozen=True)
class VerifiedBinary:
    path: Path
    descriptor: int
    device: int
    inode: int
    byte_count: int
    sha256: str
    build_id: str
    source_commit: str
    test_hooks_verified: bool


@dataclass(frozen=True)
class InvocationResult:
    invocation: Invocation
    argv_sha256: str
    stdout_sha256: str
    metadata_line: str
    metadata_sha256: str
    rows: Tuple[Mapping[str, Any], ...]


def make_invocations() -> Tuple[Invocation, ...]:
    result: List[Invocation] = []
    for layout in LAYOUTS:
        for attempt in ATTEMPTS:
            for root_index, root in enumerate(ROOTS):
                for schedule, loss_arg, loss_ppm in SCHEDULES:
                    result.append(Invocation(
                        len(result), layout, attempt, root_index, root,
                        schedule, loss_arg, loss_ppm))
    if len(result) != EXPECTED_INVOCATIONS:
        fail("internal invocation roster is not exactly 96 commands")
    return tuple(result)


def _read_stable_fd(descriptor: int, cap: int, context: str) -> bytes:
    try:
        before = os.fstat(descriptor)
        if (not stat.S_ISREG(before.st_mode) or before.st_size <= 0 or
                before.st_size > cap):
            fail("{} is not a bounded nonempty regular file".format(context))
        blocks: List[bytes] = []
        offset = 0
        while offset < before.st_size:
            block = os.pread(
                descriptor, min(1024 * 1024, before.st_size - offset), offset)
            if not block:
                fail("{} was truncated while reading".format(context))
            blocks.append(block)
            offset += len(block)
        if os.pread(descriptor, 1, offset):
            fail("{} grew while reading".format(context))
        after = os.fstat(descriptor)
    except OSError as exc:
        fail("cannot read {}: {}".format(context, exc))
    stable_before = (
        before.st_dev, before.st_ino, before.st_size,
        before.st_mtime_ns, before.st_ctime_ns)
    stable_after = (
        after.st_dev, after.st_ino, after.st_size,
        after.st_mtime_ns, after.st_ctime_ns)
    if stable_before != stable_after:
        fail("{} changed while reading".format(context))
    return b"".join(blocks)


def _gnu_build_id(data: bytes) -> str:
    if len(data) < 64 or data[:4] != b"\x7fELF":
        fail("benchmark is not an ELF file")
    elf_class = data[4]
    byte_order = data[5]
    if elf_class not in (1, 2) or byte_order not in (1, 2):
        fail("benchmark has an unsupported ELF class or byte order")
    endian = "<" if byte_order == 1 else ">"
    if elf_class == 2:
        if len(data) < 64:
            fail("benchmark has a truncated ELF64 header")
        phoff = struct.unpack_from(endian + "Q", data, 32)[0]
        phentsize = struct.unpack_from(endian + "H", data, 54)[0]
        phnum = struct.unpack_from(endian + "H", data, 56)[0]
        minimum = 56
    else:
        if len(data) < 52:
            fail("benchmark has a truncated ELF32 header")
        phoff = struct.unpack_from(endian + "I", data, 28)[0]
        phentsize = struct.unpack_from(endian + "H", data, 42)[0]
        phnum = struct.unpack_from(endian + "H", data, 44)[0]
        minimum = 32
    if phnum == 0 or phentsize < minimum or phnum > 65535:
        fail("benchmark has no bounded ELF program-header table")
    if phoff > len(data) or phnum * phentsize > len(data) - phoff:
        fail("benchmark ELF program-header table is out of bounds")
    found = set()
    for index in range(phnum):
        offset = phoff + index * phentsize
        p_type = struct.unpack_from(endian + "I", data, offset)[0]
        if p_type != 4:  # PT_NOTE
            continue
        if elf_class == 2:
            note_offset = struct.unpack_from(endian + "Q", data, offset + 8)[0]
            note_size = struct.unpack_from(endian + "Q", data, offset + 32)[0]
        else:
            note_offset = struct.unpack_from(endian + "I", data, offset + 4)[0]
            note_size = struct.unpack_from(endian + "I", data, offset + 16)[0]
        if note_offset > len(data) or note_size > len(data) - note_offset:
            fail("benchmark ELF note segment is out of bounds")
        cursor = note_offset
        end = note_offset + note_size
        while cursor < end:
            if end - cursor < 12:
                fail("benchmark ELF note segment is truncated")
            namesz, descsz, note_type = struct.unpack_from(
                endian + "III", data, cursor)
            cursor += 12
            name_end = cursor + namesz
            name_padded = cursor + ((namesz + 3) & ~3)
            desc_end = name_padded + descsz
            desc_padded = name_padded + ((descsz + 3) & ~3)
            if (name_end > end or name_padded > end or desc_end > end or
                    desc_padded > end):
                fail("benchmark ELF note entry is out of bounds")
            name = data[cursor:name_end]
            desc = data[name_padded:desc_end]
            if note_type == 3 and name.rstrip(b"\0") == b"GNU" and desc:
                found.add(desc.hex())
            cursor = desc_padded
    if len(found) != 1:
        fail("benchmark must contain one unique nonempty GNU build ID")
    return _require_build_id(next(iter(found)))


def _missing_test_hook_markers(data: bytes) -> List[str]:
    return [marker for marker in REQUIRED_TEST_HOOK_MARKERS
            if marker.encode("ascii") not in data]


def verify_binary(
        path: Path, expected_sha256: str, expected_build_id: str,
        expected_source_commit: str) -> VerifiedBinary:
    expected_sha256 = _require_sha256(expected_sha256, "binary SHA-256")
    expected_build_id = _require_build_id(expected_build_id)
    expected_source_commit = _require_commit(expected_source_commit)
    path = Path(path).absolute()
    try:
        original = os.lstat(str(path))
        if stat.S_ISLNK(original.st_mode) or not stat.S_ISREG(original.st_mode):
            fail("benchmark path must name a nonsymlink regular file")
        flags = (os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0) |
                 getattr(os, "O_CLOEXEC", 0))
        descriptor = os.open(str(path), flags)
        opened = os.fstat(descriptor)
        current = os.stat(str(path), follow_symlinks=False)
        if (not opened.st_mode & 0o111 or
                (opened.st_dev, opened.st_ino) !=
                (current.st_dev, current.st_ino)):
            fail("benchmark must be one stable executable regular file")
        data = _read_stable_fd(descriptor, MAX_BINARY_BYTES, "benchmark")
        actual_sha256 = sha256_bytes(data)
        actual_build_id = _gnu_build_id(data)
        if actual_sha256 != expected_sha256:
            fail("benchmark SHA-256 disagrees with the explicit identity")
        if actual_build_id != expected_build_id:
            fail("benchmark build ID disagrees with the explicit identity")
        if expected_source_commit.encode("ascii") not in data:
            fail("benchmark does not embed the explicit source commit")
        missing_markers = _missing_test_hook_markers(data)
        if missing_markers:
            fail("benchmark lacks required BUILD_TESTS/test-hook markers: {}"
                 .format(",".join(missing_markers)))
        return VerifiedBinary(
            path, descriptor, opened.st_dev, opened.st_ino, opened.st_size,
            actual_sha256, actual_build_id, expected_source_commit, True)
    except BaseException:
        if "descriptor" in locals():
            try:
                os.close(descriptor)
            except OSError:
                pass
        raise


def _controller_identity() -> Mapping[str, str]:
    path = Path(__file__).resolve(strict=True)
    try:
        before = path.stat()
        data = path.read_bytes()
        after = path.stat()
    except OSError as exc:
        fail("cannot bind controller source identity: {}".format(exc))
    if (not stat.S_ISREG(before.st_mode) or
            (before.st_dev, before.st_ino, before.st_size,
             before.st_mtime_ns, before.st_ctime_ns) !=
            (after.st_dev, after.st_ino, after.st_size,
             after.st_mtime_ns, after.st_ctime_ns) or
            len(data) != before.st_size):
        fail("controller source changed while hashing")
    return {"path": str(path), "sha256": sha256_bytes(data)}


def build_manifest(
        benchmark_identity: Mapping[str, Any],
        controller_identity: Mapping[str, Any]) -> Mapping[str, Any]:
    benchmark = Path(str(benchmark_identity["path"])).absolute()
    invocations = []
    for invocation in make_invocations():
        argv = invocation.argv(benchmark)
        invocations.append({
            "argv": argv,
            "argv_sha256": sha256_json(argv),
            "identity": invocation.identity(),
            "ordinal": invocation.ordinal,
        })
    unsigned = {
        "benchmark": {
            "build_id": _require_build_id(
                str(benchmark_identity["build_id"])),
            "embedded_source_commit": _require_commit(
                str(benchmark_identity["embedded_source_commit"])),
            "path": str(benchmark),
            "sha256": _require_sha256(
                str(benchmark_identity["sha256"]), "binary SHA-256"),
            "test_hooks_verified":
                benchmark_identity.get("test_hooks_verified") is True,
        },
        "controller": {
            "path": str(Path(str(controller_identity["path"])).absolute()),
            "sha256": _require_sha256(
                str(controller_identity["sha256"]), "controller SHA-256"),
        },
        "expected_invocations": EXPECTED_INVOCATIONS,
        "expected_records": EXPECTED_RECORDS,
        "experiment": {
            "attempts": list(ATTEMPTS),
            "block_bytes": 1280,
            "binary_dense_rows": 12,
            "construction_seed_basis": RAW_SEED_BASIS,
            "gf256_heavy_rows": 12,
            "heavy_family": "periodic",
            "K": list(K_VALUES),
            "layouts": list(LAYOUTS),
            "mix_count_cli_order": list(MIX_COUNTS),
            "overheads": list(OVERHEADS),
            "paired_overhead_stream": True,
            "roots": list(ROOTS),
            "schedules": [
                {"loss_arg": loss, "loss_ppm": ppm, "name": name}
                for name, loss, ppm in SCHEDULES
            ],
            "seed_schedule_sha256": RAW_SEED_SCHEDULE_SHA256,
            "threads": 1,
            "trials": 1,
        },
        "invocations": invocations,
        "policy": {
            "cold_solve_wide_xor": "policy",
            "child_environment": dict(CHILD_ENVIRONMENT),
            "diagnostic_timings_in_decision": False,
            "full_payload_solve": False,
            "payload_e2e": False,
            "runtime_delivered_id_hash_claimed": False,
            "runtime_trace_hash_claimed": False,
            "scientific_pairing_basis":
                "frozen-source+binary+argv-derived",
            "seed_perturbations": False,
            "sequential_execution": True,
            "required_test_hook_markers": list(REQUIRED_TEST_HOOK_MARKERS),
            "required_test_hook_markers_sha256":
                sha256_json(list(REQUIRED_TEST_HOOK_MARKERS)),
            "required_compile_capability":
                "WIREHAIR_V2_ENABLE_TEST_HOOKS (BUILD_TESTS)",
        },
        "schema": MANIFEST_SCHEMA,
    }
    if not unsigned["benchmark"]["test_hooks_verified"]:
        fail("manifest requires a verified BUILD_TESTS/test-hook benchmark")
    result = dict(unsigned)
    result["manifest_sha256"] = sha256_json(unsigned)
    return result


def _manifest_for_verified(binary: VerifiedBinary) -> Mapping[str, Any]:
    return build_manifest({
        "path": str(binary.path),
        "sha256": binary.sha256,
        "build_id": binary.build_id,
        "embedded_source_commit": binary.source_commit,
        "test_hooks_verified": binary.test_hooks_verified,
    }, _controller_identity())


def _parse_uint(text: str, field: str, maximum: int = MASK64) -> int:
    if UINT_TEXT.fullmatch(text) is None:
        fail("{} is not a canonical unsigned integer".format(field))
    value = int(text)
    if value > maximum:
        fail("{} exceeds its supported range".format(field))
    return value


def _parse_counter_mean(text: str, field: str) -> int:
    if THREE_DECIMAL.fullmatch(text) is None:
        fail("{} is not a canonical three-decimal mean".format(field))
    try:
        value = Decimal(text)
    except InvalidOperation:
        fail("{} is not a decimal".format(field))
    if value != value.to_integral_value() or value > MASK64:
        fail("{} is not an integral one-trial counter".format(field))
    return int(value)


def _parse_timing(text: str, field: str) -> str:
    if THREE_DECIMAL.fullmatch(text) is None:
        fail("{} is not a canonical three-decimal timing".format(field))
    try:
        value = Decimal(text)
    except InvalidOperation:
        fail("{} is not a decimal timing".format(field))
    if not value.is_finite() or value < 0:
        fail("{} is not a finite nonnegative timing".format(field))
    return text


def _parse_csv_line(line: str) -> List[str]:
    try:
        rows = list(csv.reader([line], strict=True))
    except csv.Error as exc:
        fail("precodefail CSV is malformed: {}".format(exc))
    if len(rows) != 1 or len(rows[0]) != len(CSV_HEADER):
        fail("precodefail CSV row is not 43 aligned fields")
    if line != ",".join(rows[0]):
        fail("precodefail CSV is not in canonical unquoted form")
    return rows[0]


def _parse_metadata(
        line: str, invocation: Invocation,
        expected_source_commit: str) -> Mapping[str, str]:
    prefix = "# precodefail: "
    if not line.startswith(prefix):
        fail("precodefail metadata line is missing")
    values: Dict[str, str] = {}
    order: List[str] = []
    for token in line[len(prefix):].split(" "):
        if token.count("=") != 1:
            fail("precodefail metadata token is malformed")
        key, value = token.split("=", 1)
        if not key or not value or key in values:
            fail("precodefail metadata has a duplicate or empty field")
        order.append(key)
        values[key] = value
    if tuple(order) != METADATA_KEYS:
        fail("precodefail metadata schema or field order changed")
    expected = {
        "trials": "1", "threads": "1",
        "loss": "0.10000000000000001" if invocation.loss_ppm == 100000
                else "0.5",
        "seed": "0x{:x}".format(int(invocation.root, 16)),
        "source_hits_override": "0", "packet_peel_seed_xor": "0x0",
        "binary_dense_rows_override": "12",
        "gf256_heavy_rows_override": "12",
        "dense_anchor_layout": invocation.layout,
        "odd_packet_peel_seed_xor": "0x0",
        "packet_row_seed_multiplier": "0x1",
        "packet_row_seed_avalanche": "0",
        "seed_block_bytes_override": "0", "overhead_stream": "paired",
        "full_payload_solve": "0", "schedule": invocation.schedule,
        "cold_solve_wide_xor": "policy", "exact_attempt_mode": "1",
        "exact_precode_attempt": str(invocation.attempt),
        "exact_packet_attempt": str(invocation.attempt),
        "construction_seed_basis": RAW_SEED_BASIS,
        "seed_schedule_sha256": RAW_SEED_SCHEDULE_SHA256,
        "source_git_commit": expected_source_commit,
    }
    if values != expected:
        differing = [key for key in METADATA_KEYS
                     if values.get(key) != expected.get(key)]
        fail("precodefail metadata disagrees with invocation: {}".format(
            ",".join(differing)))
    return values


def parse_invocation_output(
        invocation: Invocation, command: Sequence[str], stdout: bytes,
        stderr: bytes, binary_identity: Mapping[str, str],
        controller_sha256: str) -> InvocationResult:
    if len(stdout) > MAX_STDOUT_BYTES or len(stderr) > MAX_STDERR_BYTES:
        fail("precodefail output exceeded its fixed byte cap")
    if stderr:
        fail("precodefail emitted unexpected stderr")
    if not stdout.endswith(b"\n") or b"\r" in stdout or b"\0" in stdout:
        fail("precodefail stdout has invalid framing")
    try:
        lines = stdout.decode("ascii").splitlines()
    except UnicodeDecodeError:
        fail("precodefail stdout is not ASCII")
    if len(lines) != ROWS_PER_INVOCATION + 2:
        fail("precodefail invocation must emit metadata, header, and 64 rows")
    _parse_metadata(
        lines[0], invocation, binary_identity["embedded_source_commit"])
    if tuple(_parse_csv_line(lines[1])) != CSV_HEADER:
        fail("precodefail CSV header changed")
    expected_command = invocation.argv(Path(command[0]))
    if list(command) != expected_command:
        fail("precodefail command escaped the closed argv allowlist")
    argv_sha256 = sha256_json(list(command))
    stdout_sha256 = sha256_bytes(stdout)
    metadata_sha256 = sha256_bytes((lines[0] + "\n").encode("ascii"))
    expected_precode_seed = _effective_precode_seed(invocation.attempt)
    expected_packet_seed = _effective_packet_seed(invocation.attempt)
    expected_order = [
        (K, mix_count, overhead)
        for K in K_VALUES for mix_count in MIX_COUNTS for overhead in OVERHEADS
    ]
    seen = set()
    parsed: List[Mapping[str, Any]] = []
    for row_ordinal, (line, expected_coordinate) in enumerate(
            zip(lines[2:], expected_order)):
        row = dict(zip(CSV_HEADER, _parse_csv_line(line)))
        K = _parse_uint(row["N"], "N", 64000)
        mix_count = _parse_uint(row["mix_count"], "mix_count", 3)
        overhead = _parse_uint(row["overhead"], "overhead", 4)
        coordinate = (K, mix_count, overhead)
        if coordinate != expected_coordinate or coordinate in seen:
            fail("precodefail row order, coordinate, or uniqueness changed")
        seen.add(coordinate)
        arm = ARM_BY_LAYOUT_MIX.get((invocation.layout, mix_count))
        if arm is None:
            fail("precodefail row does not derive one closed A/B/C/D arm")
        exact = {
            "bb": "1280", "heavy_family": "periodic",
            "staircase": str(STAIRCASE_BY_K[K]),
            "binary_dense_rows": "12", "gf256_heavy_rows": "12",
            "source_hits": "3", "dense_identity_corner": "0",
            "trials": "1", "seed_attempt": "",
            "active_packet_peel_seed_xor": "0x0",
            "precode_attempt": str(invocation.attempt),
            "packet_attempt": str(invocation.attempt),
            "attempt_mode": "exact",
            "construction_seed_basis": RAW_SEED_BASIS,
            "seed_schedule_sha256": RAW_SEED_SCHEDULE_SHA256,
            "effective_precode_seed": expected_precode_seed,
            "effective_packet_seed": expected_packet_seed,
        }
        for field, expected in exact.items():
            if row[field] != expected:
                fail("precodefail row field {} disagrees with H12/current "
                     "configuration".format(field))
        success = _parse_uint(row["success"], "success", 1)
        rank_fail = _parse_uint(row["rank_fail"], "rank_fail", 1)
        errors = _parse_uint(row["error"], "error", 1)
        if success + rank_fail + errors != 1 or errors != 0:
            fail("precodefail row has an error or invalid terminal counts")
        if row["fail_rate"] != ("1.00000000" if rank_fail else "0.00000000"):
            fail("precodefail fail_rate disagrees with one-trial result")
        inactivated = _parse_counter_mean(row["inact_mu"], "inact_mu")
        if (_parse_uint(row["inact_max"], "inact_max") != inactivated or
                inactivated > 1024):
            fail("precodefail inactivation mean/max/cap is invalid")
        binary_deficiency = _parse_counter_mean(
            row["binary_def_mu"], "binary_def_mu")
        gf256_rank_gain = _parse_counter_mean(
            row["heavy_gain_mu"], "heavy_gain_mu")
        if (_parse_uint(row["binary_def_max"], "binary_def_max") !=
                binary_deficiency or
                _parse_uint(row["heavy_gain_min"], "heavy_gain_min") !=
                gf256_rank_gain or binary_deficiency > inactivated or
                gf256_rank_gain > binary_deficiency or
                bool(success) != (gf256_rank_gain == binary_deficiency) or
                bool(rank_fail) != (gf256_rank_gain < binary_deficiency)):
            fail("precodefail rank counters disagree with terminal result")
        heavy_shortfall = _parse_uint(
            row["heavy_shortfall"], "heavy_shortfall", 1)
        expected_shortfall = int(
            bool(rank_fail) and binary_deficiency <= 12 and
            gf256_rank_gain < binary_deficiency)
        if heavy_shortfall != expected_shortfall:
            fail("precodefail heavy_shortfall is inconsistent")
        block_xors = _parse_counter_mean(row["block_xors_mu"], "block_xors_mu")
        gf256_muladds = _parse_counter_mean(
            row["block_muladds_mu"], "block_muladds_mu")
        timings = {field: _parse_timing(row[field], field)
                   for field in TIMING_FIELDS}
        expected_first = "0" if rank_fail else "-1"
        expected_failure_trials = "0" if rank_fail else ""
        if (row["first_rank_fail"] != expected_first or
                row["failure_trials"] != expected_failure_trials or
                row["binary_def_hist"] != "{}:1".format(binary_deficiency) or
                row["heavy_gain_hist"] != "{}:1".format(gf256_rank_gain)):
            fail("precodefail one-trial diagnostics are inconsistent")
        parsed.append({
            "active_packet_peel_seed_xor": "0x0",
            "arm": arm,
            "attempt_mode": "exact",
            "binary_build_id": binary_identity["build_id"],
            "binary_deficiency": binary_deficiency,
            "binary_dense_rows": 12,
            "binary_sha256": binary_identity["sha256"],
            "block_bytes": 1280,
            "block_xors": block_xors,
            "construction_attempt": invocation.attempt,
            "construction_seed_basis": RAW_SEED_BASIS,
            "controller_sha256": controller_sha256,
            "csv_row_sha256": sha256_bytes((line + "\n").encode("ascii")),
            "dense_anchor_layout": invocation.layout,
            "dense_identity_corner": False,
            "diagnostic_timings": timings,
            "effective_packet_seed": expected_packet_seed,
            "effective_precode_seed": expected_precode_seed,
            "error": errors,
            "fail_rate": row["fail_rate"],
            "failure_trials": row["failure_trials"],
            "first_rank_fail": int(row["first_rank_fail"]),
            "gf256_heavy_rows": 12,
            "gf256_muladds": gf256_muladds,
            "gf256_rank_gain": gf256_rank_gain,
            "heavy_family": "periodic",
            "heavy_shortfall": heavy_shortfall,
            "inactivated_columns": inactivated,
            "invocation_ordinal": invocation.ordinal,
            "K": K,
            "legacy_seed_attempt": "",
            "loss_ppm": invocation.loss_ppm,
            "loss_root": invocation.root,
            "metadata_sha256": metadata_sha256,
            "mix_count": mix_count,
            "overhead": overhead,
            "packet_attempt": invocation.attempt,
            "precode_attempt": invocation.attempt,
            "rank_fail": rank_fail,
            "record_ordinal":
                invocation.ordinal * ROWS_PER_INVOCATION + row_ordinal,
            "root_index": invocation.root_index,
            "row_ordinal": row_ordinal,
            "schedule": invocation.schedule,
            "schema": RECORD_SCHEMA,
            "seed_schedule_sha256": RAW_SEED_SCHEDULE_SHA256,
            "source_commit": binary_identity["embedded_source_commit"],
            "source_hits": 3,
            "staircase": STAIRCASE_BY_K[K],
            "stdout_sha256": stdout_sha256,
            "success": success,
            "trials": 1,
            "argv_sha256": argv_sha256,
        })
    if seen != set(expected_order):
        fail("precodefail invocation omitted one or more closed rows")
    return InvocationResult(
        invocation, argv_sha256, stdout_sha256, lines[0], metadata_sha256,
        tuple(parsed))


def _ratio(
        candidate: int, control: int, limit_numerator: int,
        limit_denominator: int, common_count: int) -> Mapping[str, Any]:
    if min(candidate, control, common_count) < 0:
        fail("work ratio received a negative counter")
    if common_count == 0:
        passed = False
        zero_rule = "no-common-success"
    elif control == 0:
        passed = candidate == 0
        zero_rule = "both-zero-pass" if passed else "positive-over-zero-fail"
    else:
        passed = (limit_denominator * candidate <=
                  limit_numerator * control)
        zero_rule = "ordinary"
    return {
        "candidate_sum": candidate,
        "common_success_count": common_count,
        "control_sum": control,
        "limit_denominator": limit_denominator,
        "limit_numerator": limit_numerator,
        "pass": passed,
        "zero_denominator_rule": zero_rule,
    }


def _resign_manifest(manifest: Mapping[str, Any]) -> Mapping[str, Any]:
    unsigned = dict(manifest)
    unsigned.pop("manifest_sha256", None)
    result = dict(unsigned)
    result["manifest_sha256"] = sha256_json(unsigned)
    return result


def _validate_manifest(manifest: Mapping[str, Any]) -> Tuple[Invocation, ...]:
    if not isinstance(manifest, dict) or manifest.get("schema") != MANIFEST_SCHEMA:
        fail("manifest schema changed")
    supplied_hash = manifest.get("manifest_sha256")
    unsigned = dict(manifest)
    unsigned.pop("manifest_sha256", None)
    if supplied_hash != sha256_json(unsigned):
        fail("manifest self-hash is invalid")
    benchmark = manifest.get("benchmark")
    controller = manifest.get("controller")
    if not isinstance(benchmark, dict) or not isinstance(controller, dict):
        fail("manifest identities are missing")
    expected = build_manifest(benchmark, controller)
    if canonical_json(expected) != canonical_json(manifest):
        fail("manifest is not the exact closed 96-command roster")
    return make_invocations()


def _validate_records(
        records: Sequence[Mapping[str, Any]], manifest: Mapping[str, Any]
        ) -> List[Mapping[str, Any]]:
    invocations = _validate_manifest(manifest)
    if len(records) != EXPECTED_RECORDS:
        fail("result set must contain exactly 6144 records")
    manifest_entries = {entry["ordinal"]: entry
                        for entry in manifest["invocations"]}
    expected_coordinates = {
        (ARM_BY_LAYOUT_MIX[(inv.layout, mix)], inv.attempt, inv.root_index,
         inv.schedule, K, overhead)
        for inv in invocations for mix in MIX_COUNTS
        for K in K_VALUES for overhead in OVERHEADS
    }
    seen = set()
    normalized: List[Mapping[str, Any]] = []
    sha_fields = (
        "argv_sha256", "stdout_sha256", "metadata_sha256",
        "csv_row_sha256", "binary_sha256", "controller_sha256",
        "seed_schedule_sha256",
    )
    invocation_receipts: Dict[int, Tuple[str, str]] = {}
    for record_ordinal, source in enumerate(records):
        if (not isinstance(source, dict) or set(source) != RECORD_FIELDS or
                source.get("schema") != RECORD_SCHEMA):
            fail("record schema changed")
        record = dict(source)
        invocation_ordinal = record.get("invocation_ordinal")
        row_ordinal = record.get("row_ordinal")
        if (type(invocation_ordinal) is not int or
                not 0 <= invocation_ordinal < EXPECTED_INVOCATIONS or
                type(row_ordinal) is not int or
                row_ordinal != record_ordinal % ROWS_PER_INVOCATION or
                invocation_ordinal != record_ordinal // ROWS_PER_INVOCATION or
                type(record.get("record_ordinal")) is not int or
                record.get("record_ordinal") != record_ordinal):
            fail("record/invocation ordinals are not canonical")
        invocation = invocations[invocation_ordinal]
        entry = manifest_entries[invocation_ordinal]
        integer_identity_fields = (
            "block_bytes", "construction_attempt", "gf256_heavy_rows",
            "binary_dense_rows", "first_rank_fail", "loss_ppm",
            "mix_count", "overhead", "packet_attempt", "precode_attempt",
            "record_ordinal", "root_index", "source_hits", "staircase",
            "trials",
        )
        if any(type(record.get(field)) is not int
               for field in integer_identity_fields):
            fail("record identity/configuration integers are not canonical")
        expected_arm = ARM_BY_LAYOUT_MIX.get(
            (invocation.layout, record.get("mix_count")))
        coordinate = (
            record.get("arm"), record.get("construction_attempt"),
            record.get("root_index"), record.get("schedule"),
            record.get("K"), record.get("overhead"))
        if expected_arm is None or record.get("arm") != expected_arm:
            fail("record arm is not derived from layout and mix count")
        if coordinate in seen or coordinate not in expected_coordinates:
            fail("record coordinate is missing, duplicate, or outside roster")
        seen.add(coordinate)
        if (record.get("dense_anchor_layout") != invocation.layout or
                record.get("construction_attempt") != invocation.attempt or
                record.get("precode_attempt") != invocation.attempt or
                record.get("packet_attempt") != invocation.attempt or
                record.get("root_index") != invocation.root_index or
                record.get("loss_root") != invocation.root or
                record.get("schedule") != invocation.schedule or
                record.get("loss_ppm") != invocation.loss_ppm or
                record.get("argv_sha256") != entry["argv_sha256"]):
            fail("record trusted invocation identity was mutated")
        K = record.get("K")
        if (type(K) is not int or K not in K_VALUES or
                record.get("staircase") != STAIRCASE_BY_K[K] or
                record.get("block_bytes") != 1280 or
                record.get("binary_dense_rows") != 12 or
                record.get("gf256_heavy_rows") != 12 or
                record.get("source_hits") != 3 or
                record.get("dense_identity_corner") is not False or
                record.get("heavy_family") != "periodic" or
                record.get("trials") != 1 or
                record.get("legacy_seed_attempt") != "" or
                record.get("active_packet_peel_seed_xor") != "0x0" or
                record.get("attempt_mode") != "exact" or
                record.get("construction_seed_basis") != RAW_SEED_BASIS or
                record.get("seed_schedule_sha256") !=
                    RAW_SEED_SCHEDULE_SHA256 or
                record.get("effective_precode_seed") !=
                    _effective_precode_seed(invocation.attempt) or
                record.get("effective_packet_seed") !=
                    _effective_packet_seed(invocation.attempt)):
            fail("record is not the current pure-GF256 H12 configuration")
        if (record.get("binary_sha256") != manifest["benchmark"]["sha256"] or
                record.get("binary_build_id") !=
                    manifest["benchmark"]["build_id"] or
                record.get("source_commit") !=
                    manifest["benchmark"]["embedded_source_commit"] or
                record.get("controller_sha256") !=
                    manifest["controller"]["sha256"] or
                any(not isinstance(record.get(field), str) or
                    HEX64.fullmatch(record[field]) is None
                    for field in sha_fields)):
            fail("record binary/controller/hash provenance is invalid")
        receipt = (record["stdout_sha256"], record["metadata_sha256"])
        previous = invocation_receipts.setdefault(invocation_ordinal, receipt)
        if previous != receipt:
            fail("rows from one invocation disagree on stdout/metadata binding")
        terminal = tuple(record.get(field)
                         for field in ("success", "rank_fail", "error"))
        if (any(type(value) is not int or value not in (0, 1)
                for value in terminal) or sum(terminal) != 1 or
                record.get("error") != 0):
            fail("record terminal counts are invalid")
        counters = tuple(record.get(field) for field in (
            "inactivated_columns", "binary_deficiency", "gf256_rank_gain",
            "heavy_shortfall", "block_xors", "gf256_muladds"))
        if any(type(value) is not int or not 0 <= value <= MASK64
               for value in counters):
            fail("record work/rank counters are invalid")
        if (record["inactivated_columns"] > 1024 or
                record["binary_deficiency"] > record["inactivated_columns"] or
                record["gf256_rank_gain"] > record["binary_deficiency"] or
                bool(record["success"]) !=
                    (record["gf256_rank_gain"] ==
                     record["binary_deficiency"]) or
                bool(record["rank_fail"]) !=
                    (record["gf256_rank_gain"] <
                     record["binary_deficiency"])):
            fail("record rank result and counters disagree")
        expected_shortfall = int(
            bool(record["rank_fail"]) and
            record["binary_deficiency"] <= 12 and
            record["gf256_rank_gain"] < record["binary_deficiency"])
        if record["heavy_shortfall"] != expected_shortfall:
            fail("record heavy_shortfall disagrees with rank counters")
        expected_failure = "0" if record["rank_fail"] else ""
        expected_first = 0 if record["rank_fail"] else -1
        if (record.get("failure_trials") != expected_failure or
                record.get("first_rank_fail") != expected_first or
                record.get("fail_rate") !=
                    ("1.00000000" if record["rank_fail"] else "0.00000000")):
            fail("record one-trial diagnostics disagree")
        timings = record.get("diagnostic_timings")
        if (not isinstance(timings, dict) or set(timings) != set(TIMING_FIELDS) or
                any(_parse_timing(value, key) != value
                    for key, value in timings.items())):
            fail("record diagnostic timing schema is invalid")
        normalized.append(record)
    if seen != expected_coordinates or len(invocation_receipts) != 96:
        fail("record set does not exactly cover the closed paired domain")
    return normalized


def _decision_projection(records: Sequence[Mapping[str, Any]]) -> Mapping[str, Any]:
    fields = (
        "arm", "construction_attempt", "root_index", "loss_root",
        "schedule", "loss_ppm", "K", "overhead", "mix_count",
        "staircase", "binary_dense_rows", "gf256_heavy_rows",
        "source_hits", "success", "rank_fail", "error",
        "inactivated_columns", "binary_deficiency", "gf256_rank_gain",
        "heavy_shortfall", "block_xors", "gf256_muladds",
    )
    return {
        "gates": {
            "global_inact": [100, 100],
            "global_muladd": [99, 100],
            "global_xor": [95, 100],
            "per_K_and_schedule": [101, 100],
            "recovery_reduction": [95, 100],
        },
        "records": [{field: row[field] for field in fields}
                    for row in records],
        "schema": DECISION_INPUT_SCHEMA,
    }


def _failure_count_violations(
        rows: Sequence[Mapping[str, Any]]) -> List[Mapping[str, Any]]:
    dimensions: Tuple[Tuple[str, ...], ...] = (
        (), ("schedule",), ("construction_attempt",), ("loss_root",),
        ("K",), ("overhead",),
    )
    violations = []
    for dimensions_used in dimensions:
        keys = {tuple(row[field] for field in dimensions_used) for row in rows}
        for values in sorted(keys, key=lambda item: tuple(str(x) for x in item)):
            subset = [row for row in rows if all(
                row[field] == value
                for field, value in zip(dimensions_used, values))]
            counts = {arm: sum(row["rank_fail"] for row in subset
                               if row["arm"] == arm)
                      for arm in ("A", "C", "D")}
            for reference in ("A", "C"):
                if counts["D"] > counts[reference]:
                    violations.append({
                        "candidate_failures": counts["D"],
                        "dimension": "global" if not dimensions_used else
                            "+".join(dimensions_used),
                        "reference": reference,
                        "reference_failures": counts[reference],
                        "values": list(values),
                    })
    return violations


def reduce_records(
        source_records: Sequence[Mapping[str, Any]],
        manifest: Mapping[str, Any]) -> Mapping[str, Any]:
    records = _validate_records(source_records, manifest)
    by_key = {}
    for row in records:
        key = (row["construction_attempt"], row["root_index"],
               row["schedule"], row["K"], row["overhead"])
        by_key.setdefault(key, {})[row["arm"]] = row
    if (len(by_key) * 4 != EXPECTED_RECORDS or
            any(set(arms) != {"A", "B", "C", "D"}
                for arms in by_key.values())):
        fail("A/B/C/D rows are not exactly paired")

    introductions = {"A": [], "C": []}
    repairs_a = []
    for key, arms in by_key.items():
        if arms["D"]["rank_fail"] and not arms["A"]["rank_fail"]:
            introductions["A"].append(key)
        if arms["D"]["rank_fail"] and not arms["C"]["rank_fail"]:
            introductions["C"].append(key)
        if arms["A"]["rank_fail"] and not arms["D"]["rank_fail"]:
            repairs_a.append(key)
    failure_totals = {
        arm: sum(row["rank_fail"] for row in records if row["arm"] == arm)
        for arm in ("A", "B", "C", "D")
    }
    no_introductions = not introductions["A"] and not introductions["C"]
    recovery_reduction = (
        failure_totals["A"] > 0 and
        100 * failure_totals["D"] <= 95 * failure_totals["A"])
    no_regression_violations = _failure_count_violations(records)

    nesting_violations = []
    for arm in ("A", "B", "C", "D"):
        grouped: Dict[Tuple[Any, ...], Dict[int, int]] = {}
        for row in records:
            if row["arm"] != arm:
                continue
            key = (row["construction_attempt"], row["root_index"],
                   row["schedule"], row["K"])
            grouped.setdefault(key, {})[row["overhead"]] = row["rank_fail"]
        for key, failures in grouped.items():
            for lower, higher in zip(OVERHEADS, OVERHEADS[1:]):
                if failures[higher] > failures[lower]:
                    nesting_violations.append({
                        "arm": arm, "base_key": list(key),
                        "higher_overhead": higher, "lower_overhead": lower,
                    })

    common = [(arms["A"], arms["D"]) for arms in by_key.values()
              if arms["A"]["success"] and arms["D"]["success"]]

    def work_summary(
            pairs: Sequence[Tuple[Mapping[str, Any], Mapping[str, Any]]],
            limits: Mapping[str, int]) -> Mapping[str, Any]:
        return {
            field: _ratio(
                sum(candidate[field] for _, candidate in pairs),
                sum(control[field] for control, _ in pairs),
                numerator, 100, len(pairs))
            for field, numerator in limits.items()
        }

    global_work = work_summary(common, {
        "block_xors": 95,
        "gf256_muladds": 99,
        "inactivated_columns": 100,
    })
    global_work_pass = all(value["pass"] for value in global_work.values())
    stratum_work = {"K": [], "schedule": []}
    stratum_work_pass = True
    for dimension, values in (("K", K_VALUES),
                              ("schedule", tuple(x[0] for x in SCHEDULES))):
        for value in values:
            pairs = [(control, candidate) for control, candidate in common
                     if control[dimension] == value]
            metrics = work_summary(pairs, {
                "block_xors": 101,
                "gf256_muladds": 101,
                "inactivated_columns": 101,
            })
            passed = all(metric["pass"] for metric in metrics.values())
            stratum_work_pass = stratum_work_pass and passed
            stratum_work[dimension].append({
                "metrics": metrics, "pass": passed, "value": value,
            })

    hard_fail = (
        not no_introductions or bool(no_regression_violations) or
        bool(nesting_violations) or not global_work_pass or
        not stratum_work_pass)
    if failure_totals["A"] == 0:
        disposition = "REJECT" if hard_fail else "INCONCLUSIVE"
    else:
        passed = (not hard_fail and bool(repairs_a) and recovery_reduction)
        disposition = "PASS" if passed else "REJECT"

    decision_projection = _decision_projection(records)
    decision = {
        "common_success_work": {
            "global": global_work,
            "per_stratum": stratum_work,
        },
        "decision_sha256": sha256_json(decision_projection),
        "diagnostic_timing_fields_included": False,
        "disposition": disposition,
        "failure_count_regression_violations": no_regression_violations,
        "failure_totals": failure_totals,
        "gates": {
            "common_success_global": global_work_pass,
            "common_success_per_K_schedule": stratum_work_pass,
            "composite_repairs_A": bool(repairs_a),
            "five_percent_reduction_vs_A": recovery_reduction,
            "no_failure_count_regression_vs_A_C":
                not no_regression_violations,
            "no_introductions_vs_A_C": no_introductions,
            "overhead_nesting": not nesting_violations,
        },
        "introduction_counts": {
            "vs_A": len(introductions["A"]),
            "vs_C": len(introductions["C"]),
        },
        "overhead_nesting_violations": nesting_violations,
        "record_count": len(records),
        "repair_count_vs_A": len(repairs_a),
        "schema": DECISION_SCHEMA,
    }
    return decision


def _kill_process(process: subprocess.Popen) -> None:
    if process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except OSError:
        if process.poll() is None:
            try:
                process.kill()
            except OSError:
                pass


def _capture_bounded(
        process: subprocess.Popen, timeout_seconds: float) -> Tuple[bytes, bytes]:
    if process.stdout is None or process.stderr is None:
        fail("benchmark pipes were not created")
    selector = selectors.DefaultSelector()
    streams = {
        process.stdout.fileno(): ("stdout", MAX_STDOUT_BYTES),
        process.stderr.fileno(): ("stderr", MAX_STDERR_BYTES),
    }
    buffers = {"stdout": bytearray(), "stderr": bytearray()}
    deadline = time.monotonic() + timeout_seconds
    try:
        for descriptor in streams:
            selector.register(descriptor, selectors.EVENT_READ)
        while selector.get_map():
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                _kill_process(process)
                fail("benchmark invocation timed out")
            events = selector.select(remaining)
            if not events:
                _kill_process(process)
                fail("benchmark invocation timed out")
            for key, _mask in events:
                descriptor = key.fd
                name, cap = streams[descriptor]
                try:
                    data = os.read(descriptor, 65536)
                except OSError as exc:
                    _kill_process(process)
                    fail("cannot read benchmark {}: {}".format(name, exc))
                if not data:
                    selector.unregister(descriptor)
                elif len(buffers[name]) + len(data) > cap:
                    _kill_process(process)
                    fail("benchmark {} exceeded its byte cap".format(name))
                else:
                    buffers[name].extend(data)
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            _kill_process(process)
            fail("benchmark invocation timed out")
        try:
            process.wait(timeout=remaining)
        except subprocess.TimeoutExpired:
            _kill_process(process)
            fail("benchmark invocation timed out")
        return bytes(buffers["stdout"]), bytes(buffers["stderr"])
    finally:
        selector.close()
        process.stdout.close()
        process.stderr.close()


def _execute_invocation(
        binary: VerifiedBinary, invocation: Invocation,
        timeout_seconds: float) -> Tuple[Sequence[str], bytes, bytes]:
    command = invocation.argv(binary.path)
    process: Optional[subprocess.Popen] = None
    try:
        pthread_sigmask = getattr(signal, "pthread_sigmask", None)
        if pthread_sigmask is None:
            fail("platform cannot mask termination signals during child spawn")
        previous_mask = pthread_sigmask(
            signal.SIG_BLOCK, CHILD_SPAWN_SIGNALS)
        try:
            try:
                process = subprocess.Popen(
                    command,
                    executable="/proc/self/fd/{}".format(binary.descriptor),
                    stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE, shell=False, close_fds=True,
                    env=CHILD_ENVIRONMENT,
                    pass_fds=(binary.descriptor,), start_new_session=True)
            except OSError as exc:
                fail("cannot start benchmark invocation: {}".format(exc))
        finally:
            pthread_sigmask(signal.SIG_SETMASK, previous_mask)
        stdout, stderr = _capture_bounded(process, timeout_seconds)
        if process.returncode != 0:
            fail("benchmark invocation {} exited {}".format(
                invocation.ordinal, process.returncode))
        return command, stdout, stderr
    except BaseException:
        if process is not None:
            _kill_process(process)
            try:
                process.wait(timeout=2.0)
            except (OSError, subprocess.TimeoutExpired):
                pass
        raise


def _verify_binary_still_pinned(binary: VerifiedBinary) -> None:
    try:
        opened = os.fstat(binary.descriptor)
        current = os.stat(str(binary.path), follow_symlinks=False)
    except OSError as exc:
        fail("benchmark became unavailable: {}".format(exc))
    if ((opened.st_dev, opened.st_ino, opened.st_size) !=
            (binary.device, binary.inode, binary.byte_count) or
            (current.st_dev, current.st_ino) != (binary.device, binary.inode)):
        fail("benchmark path no longer names the verified inode")
    data = _read_stable_fd(binary.descriptor, MAX_BINARY_BYTES, "benchmark")
    if (sha256_bytes(data) != binary.sha256 or
            _gnu_build_id(data) != binary.build_id or
            binary.source_commit.encode("ascii") not in data):
        fail("benchmark identity changed during campaign")


def _open_output_directory(path: Path) -> Tuple[int, Path, Tuple[int, int]]:
    path = Path(path).absolute()
    try:
        info = os.lstat(str(path))
        if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
            fail("output directory must be an existing nonsymlink directory")
        directory_flag = getattr(os, "O_DIRECTORY", 0)
        nofollow = getattr(os, "O_NOFOLLOW", 0)
        if not directory_flag or not nofollow:
            fail("platform cannot pin output directories without symlinks")
        flags = (os.O_RDONLY | directory_flag | nofollow |
                 getattr(os, "O_CLOEXEC", 0))
        descriptor = os.open(str(path), flags)
        opened = os.fstat(descriptor)
        if (opened.st_dev, opened.st_ino) != (info.st_dev, info.st_ino):
            fail("output directory identity changed while opening")
        for name in (RESULT_NAME, SUMMARY_NAME):
            try:
                os.stat(name, dir_fd=descriptor, follow_symlinks=False)
            except FileNotFoundError:
                continue
            fail("refusing to replace existing artifact {}".format(path / name))
        return descriptor, path, (opened.st_dev, opened.st_ino)
    except BaseException:
        if "descriptor" in locals():
            try:
                os.close(descriptor)
            except OSError:
                pass
        raise


def _stage_file(directory_fd: int, final_name: str, data: bytes) -> Tuple[str, Tuple[int, int]]:
    staged = ".{}.{}.tmp".format(final_name, next(tempfile._get_candidate_names()))
    flags = (os.O_WRONLY | os.O_CREAT | os.O_EXCL |
             getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_CLOEXEC", 0))
    descriptor = os.open(staged, flags, 0o600, dir_fd=directory_fd)
    try:
        info = os.fstat(descriptor)
        identity = (info.st_dev, info.st_ino)
        view = memoryview(data)
        while view:
            count = os.write(descriptor, view)
            if count <= 0:
                fail("short write while staging {}".format(final_name))
            view = view[count:]
        os.fsync(descriptor)
    except BaseException:
        if "identity" in locals():
            _unlink_if_identity(directory_fd, staged, identity)
        try:
            os.close(descriptor)
        except OSError:
            pass
        raise
    try:
        os.close(descriptor)
    except BaseException:
        _unlink_if_identity(directory_fd, staged, identity)
        raise
    return staged, identity


def _unlink_if_identity(
        directory_fd: int, name: str, identity: Tuple[int, int]) -> None:
    try:
        info = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
        if (info.st_dev, info.st_ino) == identity:
            os.unlink(name, dir_fd=directory_fd)
    except FileNotFoundError:
        pass


def _verify_output_directory(
        directory_fd: int, output_dir: Path,
        identity: Tuple[int, int]) -> None:
    try:
        retained = os.fstat(directory_fd)
        current = os.stat(str(output_dir), follow_symlinks=False)
    except OSError as exc:
        fail("output directory became unavailable: {}".format(exc))
    if (not stat.S_ISDIR(retained.st_mode) or
            (retained.st_dev, retained.st_ino) != identity or
            (current.st_dev, current.st_ino) != identity):
        fail("output directory identity changed during campaign")


def _publish_pair_at(
        directory_fd: int, output_dir: Path,
        directory_identity: Tuple[int, int], result: bytes,
        summary: bytes) -> None:
    # Results are linked first and the summary last as the completion marker.
    # Catchable exceptions roll both back.  Two filesystem names cannot be
    # crash-atomic: SIGKILL between links may leave results without a summary,
    # which is intentionally incomplete and makes this directory non-reusable.
    _verify_output_directory(directory_fd, output_dir, directory_identity)
    staged: List[Tuple[str, str, Tuple[int, int]]] = []
    try:
        for name, data in ((RESULT_NAME, result), (SUMMARY_NAME, summary)):
            staged_name, identity = _stage_file(directory_fd, name, data)
            staged.append((staged_name, name, identity))
        for staged_name, name, identity in staged:
            try:
                os.link(staged_name, name, src_dir_fd=directory_fd,
                        dst_dir_fd=directory_fd, follow_symlinks=False)
            except FileExistsError:
                fail("refusing to replace existing artifact {}".format(
                    Path(output_dir) / name))
        os.fsync(directory_fd)
        for staged_name, _name, identity in staged:
            _unlink_if_identity(directory_fd, staged_name, identity)
        os.fsync(directory_fd)
        _verify_output_directory(
            directory_fd, output_dir, directory_identity)
    except BaseException:
        # Sweep every expected final name, not merely links whose Python-side
        # bookkeeping completed.  An asynchronous exception can arrive after
        # link(2) succeeds but before the next bytecode records that success.
        for _staged_name, name, identity in reversed(staged):
            _unlink_if_identity(directory_fd, name, identity)
        raise
    finally:
        for staged_name, _name, identity in staged:
            _unlink_if_identity(directory_fd, staged_name, identity)


def _publish_pair(output_dir: Path, result: bytes, summary: bytes) -> None:
    directory_fd, output_dir, identity = _open_output_directory(output_dir)
    try:
        _publish_pair_at(
            directory_fd, output_dir, identity, result, summary)
    finally:
        os.close(directory_fd)


def _run_campaign(
        benchmark: Path, output_dir: Path, expected_sha256: str,
        expected_build_id: str, expected_source_commit: str,
        expected_manifest_sha256: str,
        timeout_seconds: float = 300.0) -> Mapping[str, Any]:
    if (not math.isfinite(timeout_seconds) or timeout_seconds <= 0 or
            timeout_seconds > MAX_TIMEOUT_SECONDS):
        fail("invocation timeout must be finite and in (0,3600]")
    # Destination absence is checked before any child can execute and again
    # during publication through no-overwrite hard links.
    directory_fd, pinned_output_dir, directory_identity = \
        _open_output_directory(output_dir)
    binary: Optional[VerifiedBinary] = None
    try:
        binary = verify_binary(
            benchmark, expected_sha256, expected_build_id,
            expected_source_commit)
        manifest = _manifest_for_verified(binary)
        if (manifest["manifest_sha256"] != _require_sha256(
                expected_manifest_sha256, "manifest SHA-256")):
            fail("fresh manifest disagrees with the preregistered identity")
        controller_sha256 = manifest["controller"]["sha256"]
        binary_identity = {
            "build_id": binary.build_id,
            "embedded_source_commit": binary.source_commit,
            "sha256": binary.sha256,
        }
        rows: List[Mapping[str, Any]] = []
        receipts = []
        for invocation in make_invocations():
            command, stdout, stderr = _execute_invocation(
                binary, invocation, timeout_seconds)
            result = parse_invocation_output(
                invocation, command, stdout, stderr, binary_identity,
                controller_sha256)
            rows.extend(result.rows)
            receipts.append({
                "argv_sha256": result.argv_sha256,
                "metadata_line": result.metadata_line,
                "metadata_sha256": result.metadata_sha256,
                "ordinal": invocation.ordinal,
                "row_count": len(result.rows),
                "stdout_bytes": len(stdout),
                "stdout_sha256": result.stdout_sha256,
            })
        _verify_binary_still_pinned(binary)
        if _controller_identity()["sha256"] != controller_sha256:
            fail("controller source changed during campaign")
        decision = reduce_records(rows, manifest)
        result_bytes = b"".join(
            (canonical_json(row) + "\n").encode("ascii") for row in rows)
        result_sha256 = sha256_bytes(result_bytes)
        unsigned_summary = {
            "benchmark": manifest["benchmark"],
            "controller": manifest["controller"],
            "decision": decision,
            "invocation_count": len(receipts),
            "invocations": receipts,
            "manifest_sha256": manifest["manifest_sha256"],
            "record_count": len(rows),
            "result_stream_bytes": len(result_bytes),
            "result_stream_sha256": result_sha256,
            "schema": SUMMARY_SCHEMA,
            "timing_metrics_retained_raw_only": True,
        }
        summary = dict(unsigned_summary)
        summary["summary_sha256"] = sha256_json(unsigned_summary)
        summary_bytes = (canonical_json(summary) + "\n").encode("ascii")
        _publish_pair_at(
            directory_fd, pinned_output_dir, directory_identity,
            result_bytes, summary_bytes)
        return summary
    finally:
        if binary is not None:
            os.close(binary.descriptor)
        os.close(directory_fd)


def run_campaign(
        benchmark: Path, output_dir: Path, expected_sha256: str,
        expected_build_id: str, expected_source_commit: str,
        expected_manifest_sha256: str,
        timeout_seconds: float = 300.0) -> Mapping[str, Any]:
    previous_handlers = _install_termination_handlers()
    try:
        return _run_campaign(
            benchmark, output_dir, expected_sha256, expected_build_id,
            expected_source_commit, expected_manifest_sha256,
            timeout_seconds)
    finally:
        _restore_termination_handlers(previous_handlers)


def create_manifest(
        benchmark: Path, expected_sha256: str, expected_build_id: str,
        expected_source_commit: str) -> Mapping[str, Any]:
    binary = verify_binary(
        benchmark, expected_sha256, expected_build_id, expected_source_commit)
    try:
        return _manifest_for_verified(binary)
    finally:
        os.close(binary.descriptor)


def _add_identity_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--benchmark", type=Path, required=True)
    parser.add_argument("--binary-sha256", required=True)
    parser.add_argument("--binary-build-id", required=True)
    parser.add_argument("--source-commit", required=True)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    manifest = commands.add_parser(
        "manifest", aliases=("dry-run",),
        help="verify identities and print the exact 96-command receipt")
    _add_identity_arguments(manifest)
    run = commands.add_parser("run", help="execute and reduce the closed screen")
    _add_identity_arguments(run)
    run.add_argument(
        "--manifest-sha256", required=True,
        help="exact manifest_sha256 from the preregistered dry run")
    run.add_argument("--output-dir", type=Path, required=True)
    run.add_argument("--invocation-timeout-seconds", type=float, default=300.0)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.command in ("manifest", "dry-run"):
            result = create_manifest(
                args.benchmark, args.binary_sha256, args.binary_build_id,
                args.source_commit)
        else:
            result = run_campaign(
                args.benchmark, args.output_dir, args.binary_sha256,
                args.binary_build_id, args.source_commit,
                args.manifest_sha256,
                args.invocation_timeout_seconds)
    except (ScreenError, OSError, ValueError, KeyError) as exc:
        invalid = {
            "diagnostic": str(exc), "disposition": "INVALID",
            "schema": "wirehair.wh2.mix2-rank-work-invalid.v1",
        }
        print(canonical_json(invalid), file=sys.stderr)
        return 1
    print(canonical_json(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
