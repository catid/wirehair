#!/usr/bin/python3
"""Fail-closed one-shot controller for the WH2 fixed-five packet screen.

The controller never calibrates or retries.  It authenticates one clean source
commit and one binary, invokes the six frozen cells once in their registered
order through ``taskset``, independently adjudicates the emitted paired data,
and atomically publishes a terminal evidence directory.  Exit 10 from a cell
is a scientific rejection and does not stop the remaining frozen cells.
"""

from __future__ import annotations

import argparse
import ctypes
from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import tempfile
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple


SCHEMA_VERSION = 1
T31_95 = 2.039513446
SAMPLES = 32
ROWS = 256
TARGET_CPU = 50
TARGET_FULL_APIC_ID = "00000064"
TARGET_IDENTITY_CANONICAL_BYTES = 617
TARGET_IDENTITY_SHA256 = (
    "3288e0ef61cf3e628dcd827f9cf003c9d6ec6b5a12169e7a8bfc796baacddba7"
)
CELL_TIMEOUT_SECONDS = 1800.0
MAX_OUTPUT_BYTES = 2 * 1024 * 1024
LOWER40 = re.compile(r"[0-9a-f]{40}\Z")
LOWER64 = re.compile(r"[0-9a-f]{64}\Z")
LOWER16 = re.compile(r"[0-9a-f]{16}\Z")
UINT = re.compile(r"(?:0|[1-9][0-9]*)\Z")
DECIMAL = re.compile(r"(?:0|[1-9][0-9]*)\.[0-9]+\Z")

FROZEN_CELLS: Tuple[Tuple[int, int, int, int], ...] = (
    (0, 10000, 1280, 256),
    (1, 64000, 32768, 10),
    (2, 64000, 1280, 256),
    (3, 10000, 32768, 10),
    (4, 10000, 4096, 80),
    (5, 64000, 4096, 80),
)

PANELS: Tuple[Tuple[str, str, int], ...] = (
    ("active", "direct_fixed5", 2),
    ("fallback", "direct_fixed5_dispatch", 3),
    ("aa", "legacy_copy", 2),
    ("split", "wide_split2_plus3", 2),
)
EXPECTED_P = {10000: 110, 64000: 370}

DESIGN_FIELDS = (
    "version", "source_commit", "rows", "samples", "mix_count",
    "peel_seed", "precode_seed", "target_cpu", "target_full_apic_id",
    "target_identity_sha256", "candidate", "primary_gate", "control_band",
    "student_t_df31", "timed_scope", "timed_counters", "binary_sha256",
)
DESCRIBE_CELL_FIELDS = (
    "version", "ordinal", "K", "block_bytes", "repetitions",
)
META_FIELDS = (
    "version", "source_commit", "ordinal", "K", "P", "block_bytes",
    "repetitions", "rows", "samples", "mix_count", "peel_seed",
    "precode_seed",
    "active_degree", "fallback_degree", "active_terms", "fallback_terms",
    "active_rows_hash", "fallback_rows_hash", "input_hash",
    "active_first_id", "active_last_id", "fallback_first_id",
    "fallback_last_id", "wide_build", "avx2", "avx512", "candidate",
    "direct_wide_available", "timed_row_generation", "timed_pointer_gather",
    "timed_counters", "panel_order",
)
RAW_FIELDS = (
    "version", "cell", "candidate", "ordinal", "K", "block_bytes",
    "degree", "rows", "repetitions", "sample", "order", "legacy_ns",
    "candidate_ns", "ratio",
)
SUMMARY_FIELDS = (
    "version", "cell", "candidate", "ordinal", "K", "block_bytes",
    "degree", "samples", "geomean", "lower95", "upper95",
)
TARGET_FIELDS = (
    "version", "source_commit", "ordinal", "K", "block_bytes",
    "target_cpu", "full_apic_id", "pre_identity_sha256",
    "post_identity_sha256", "canonical_bytes",
    "pre_before_cpu", "pre_after_cpu", "post_before_cpu", "post_after_cpu",
    "pre_before_affinity_count", "pre_after_affinity_count",
    "post_before_affinity_count", "post_after_affinity_count",
    "pre_voluntary_delta", "pre_involuntary_delta", "post_voluntary_delta",
    "post_involuntary_delta", "gate",
)
RESULT_FIELDS = (
    "version", "ordinal", "K", "block_bytes", "primary",
    "fallback_control", "aa_control", "controls", "split_promotional",
    "mismatch_sink", "status",
)

RUN_ENV = {
    "LANG": "C",
    "LC_ALL": "C",
    "PATH": "/usr/bin:/bin",
    "PYTHONDONTWRITEBYTECODE": "1",
    "TZ": "UTC",
}
GIT_ENV = dict(RUN_ENV, GIT_CONFIG_GLOBAL="/dev/null",
               GIT_CONFIG_NOSYSTEM="1", GIT_NO_REPLACE_OBJECTS="1",
               GIT_OPTIONAL_LOCKS="0", GIT_TERMINAL_PROMPT="0")


class ScreenError(RuntimeError):
    """A fail-closed provenance, protocol, or publication error."""


@dataclass(frozen=True)
class Config:
    benchmark: Path
    output_dir: Path
    source_root: Path
    expected_commit: str
    expected_binary_sha256: str
    cpu: int
    taskset: Path = Path("/usr/bin/taskset")
    nm: Path = Path("/usr/bin/nm")
    git: Path = Path("/usr/bin/git")


@dataclass(frozen=True)
class Preflight:
    binary_sha256: str
    describe_sha256: str
    describe_text: str
    nm_sha256: str
    taskset_sha256: str
    controller_sha256: str


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(
        value, allow_nan=False, ensure_ascii=True,
        separators=(",", ":"), sort_keys=True,
    ) + "\n").encode("ascii")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        while True:
            chunk = source.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def require_hash(value: str, pattern: re.Pattern[str], name: str) -> str:
    if pattern.fullmatch(value) is None:
        raise ScreenError(name + " must be canonical lowercase hex")
    return value


def parse_uint(value: str, name: str, positive: bool = False) -> int:
    if UINT.fullmatch(value) is None:
        raise ScreenError(name + " is not a canonical unsigned integer")
    result = int(value)
    if positive and result == 0:
        raise ScreenError(name + " must be positive")
    return result


def parse_float(value: str, name: str) -> float:
    if DECIMAL.fullmatch(value) is None:
        raise ScreenError(name + " is not a canonical fixed decimal")
    result = float(value)
    if not math.isfinite(result):
        raise ScreenError(name + " is not finite")
    return result


def parse_csv_line(
        line: str, prefix: str,
        expected_fields: Sequence[str]) -> Dict[str, str]:
    pieces = line.split(",")
    if not pieces or pieces[0] != prefix:
        raise ScreenError("expected " + prefix + " line")
    fields: Dict[str, str] = {}
    for piece in pieces[1:]:
        if piece.count("=") != 1:
            raise ScreenError(prefix + " contains a malformed field")
        key, value = piece.split("=", 1)
        if not key or not value or key in fields:
            raise ScreenError(prefix + " contains an empty/duplicate field")
        fields[key] = value
    if tuple(fields) != tuple(expected_fields):
        expected_set = set(expected_fields)
        missing = sorted(expected_set - set(fields))
        extra = sorted(set(fields) - expected_set)
        raise ScreenError(
            "%s schema/order mismatch missing=%r extra=%r" %
            (prefix, missing, extra)
        )
    return fields


def decode_output(data: bytes, name: str) -> List[str]:
    if len(data) > MAX_OUTPUT_BYTES:
        raise ScreenError(name + " exceeds the bounded output size")
    if not data or not data.endswith(b"\n") or b"\x00" in data:
        raise ScreenError(name + " is empty, unterminated, or contains NUL")
    try:
        text = data.decode("ascii")
    except UnicodeDecodeError as error:
        raise ScreenError(name + " is not ASCII") from error
    lines = text.splitlines()
    if any(not line or line != line.strip() for line in lines):
        raise ScreenError(name + " contains empty/padded lines")
    return lines


def run_checked(
        argv: Sequence[str], timeout: float, name: str,
        cwd: Optional[Path] = None, env: Mapping[str, str] = RUN_ENV,
        allowed_codes: Tuple[int, ...] = (0,)) -> subprocess.CompletedProcess[bytes]:
    try:
        result = subprocess.run(
            list(argv), cwd=None if cwd is None else str(cwd), env=dict(env),
            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, timeout=timeout, check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise ScreenError(name + " execution failed: " + str(error)) from error
    if len(result.stdout) > MAX_OUTPUT_BYTES or len(result.stderr) > MAX_OUTPUT_BYTES:
        raise ScreenError(name + " output exceeded limit")
    if result.returncode not in allowed_codes:
        raise ScreenError("%s exit %d" % (name, result.returncode))
    return result


def require_regular_executable(path: Path, name: str) -> Path:
    try:
        info = path.lstat()
    except OSError as error:
        raise ScreenError(name + " is unavailable: " + str(error)) from error
    if not stat.S_ISREG(info.st_mode) or stat.S_ISLNK(info.st_mode) or \
            info.st_mode & 0o111 == 0:
        raise ScreenError(name + " must be a non-symlink executable file")
    return path.resolve(strict=True)


def git_output(config: Config, *arguments: str) -> bytes:
    result = run_checked(
        [str(config.git), "-C", str(config.source_root), *arguments],
        30.0, "git " + arguments[0], env=GIT_ENV,
    )
    if result.stderr:
        raise ScreenError("git " + arguments[0] + " wrote stderr")
    return result.stdout


def require_clean_commit(config: Config) -> None:
    head = git_output(config, "rev-parse", "--verify", "HEAD^{commit}")
    try:
        head_text = head.decode("ascii").strip()
    except UnicodeDecodeError as error:
        raise ScreenError("git HEAD is not ASCII") from error
    if head != (config.expected_commit + "\n").encode("ascii") or \
            head_text != config.expected_commit:
        raise ScreenError("source HEAD differs from expected commit")
    status_output = git_output(
        config, "status", "--porcelain=v1", "--untracked-files=all")
    if status_output:
        raise ScreenError("source worktree is not clean")


def parse_describe(raw: bytes, expected_commit: str) -> None:
    lines = decode_output(raw, "--describe stdout")
    if len(lines) != 1 + len(FROZEN_CELLS):
        raise ScreenError("--describe emitted the wrong line count")
    design = parse_csv_line(lines[0], "wh2_fixed5_design", DESIGN_FIELDS)
    expected = {
        "version": "1", "source_commit": expected_commit,
        "rows": "256", "samples": "32", "mix_count": "3",
        "peel_seed": "4d241359",
        "precode_seed": "786f72667573696f",
        "target_cpu": str(TARGET_CPU),
        "target_full_apic_id": TARGET_FULL_APIC_ID,
        "target_identity_sha256": TARGET_IDENTITY_SHA256,
        "candidate": "direct_no_write_fixed5_try",
        "primary_gate": "upper95_lt_0.99",
        "control_band": "0.99_to_1.01",
        "student_t_df31": "2.039513446",
        "timed_scope": "xor_schedules_only", "timed_counters": "0",
        "binary_sha256": "controller_required",
    }
    if design != expected:
        raise ScreenError("--describe design identity mismatch")
    for line, cell in zip(lines[1:], FROZEN_CELLS):
        fields = parse_csv_line(
            line, "wh2_fixed5_cell", DESCRIBE_CELL_FIELDS)
        ordinal, K, block_bytes, repetitions = cell
        if fields != {
            "version": "1", "ordinal": str(ordinal), "K": str(K),
            "block_bytes": str(block_bytes),
            "repetitions": str(repetitions),
        }:
            raise ScreenError("--describe frozen cell mismatch")


def audit_symbols(config: Config) -> str:
    result = run_checked(
        [str(config.nm), "-g", "--defined-only", str(config.benchmark)],
        30.0, "nm symbol audit",
    )
    if result.stderr:
        raise ScreenError("nm symbol audit wrote stderr")
    lines = decode_output(result.stdout, "nm stdout")
    symbols: List[str] = []
    for line in lines:
        pieces = line.split()
        if len(pieces) < 2:
            raise ScreenError("unrecognized nm output")
        symbols.append(pieces[-1])
    if any(symbol.startswith("gf256_count_") for symbol in symbols):
        raise ScreenError("timed GF256 counter symbols are present")
    if "gf256_try_addset5_wide_mem" not in symbols:
        raise ScreenError("direct fixed-five symbol is absent")
    return sha256_bytes(result.stdout)


def preflight(config: Config) -> Preflight:
    require_hash(config.expected_commit, LOWER40, "source commit")
    require_hash(config.expected_binary_sha256, LOWER64, "binary SHA-256")
    if config.cpu != TARGET_CPU:
        raise ScreenError("requested CPU differs from frozen CPU 50")
    if config.cpu not in os.sched_getaffinity(0):
        raise ScreenError("requested CPU is outside controller affinity")
    source_root = config.source_root.resolve(strict=True)
    output_absolute = config.output_dir.absolute()
    if not config.output_dir.is_absolute() or \
            config.output_dir.parent.resolve(strict=True) != \
            config.output_dir.parent or \
            config.output_dir.name in ("", ".", ".."):
        raise ScreenError("output path must have a canonical real parent")
    try:
        output_absolute.relative_to(source_root)
    except ValueError:
        pass
    else:
        raise ScreenError("output directory must be outside the source tree")
    if config.output_dir.exists() or config.output_dir.is_symlink():
        raise ScreenError("output directory already exists")
    if not config.output_dir.parent.is_dir():
        raise ScreenError("output parent does not exist")
    claim = attempt_path(config.output_dir)
    if claim.exists() or claim.is_symlink():
        raise ScreenError("output namespace has already been attempted")
    require_regular_executable(config.benchmark, "benchmark")
    require_regular_executable(config.taskset, "taskset")
    require_regular_executable(config.nm, "nm")
    require_regular_executable(config.git, "git")
    require_clean_commit(config)
    actual_binary_hash = sha256_file(config.benchmark)
    if actual_binary_hash != config.expected_binary_sha256:
        raise ScreenError("binary SHA-256 mismatch")
    nm_sha256 = audit_symbols(config)
    if sha256_file(config.benchmark) != actual_binary_hash:
        raise ScreenError("binary changed during symbol audit")
    describe = run_checked(
        [str(config.benchmark), "--describe"], 30.0,
        "benchmark --describe", cwd=config.source_root,
    )
    if describe.stderr:
        raise ScreenError("benchmark --describe wrote stderr")
    parse_describe(describe.stdout, config.expected_commit)
    require_clean_commit(config)
    if sha256_file(config.benchmark) != actual_binary_hash:
        raise ScreenError("binary changed during preflight")
    return Preflight(
        binary_sha256=actual_binary_hash,
        describe_sha256=sha256_bytes(describe.stdout),
        describe_text=describe.stdout.decode("ascii"),
        nm_sha256=nm_sha256,
        taskset_sha256=sha256_file(config.taskset),
        controller_sha256=sha256_file(Path(__file__).resolve()),
    )


def close_enough(emitted: float, expected: float) -> bool:
    # The producer emits exactly 12 digits after the decimal point.
    return math.isclose(emitted, expected, rel_tol=5e-15, abs_tol=5.1e-13)


def compute_stats(pairs: Sequence[Tuple[int, int]]) -> Dict[str, float]:
    if len(pairs) != SAMPLES:
        raise ScreenError("panel does not contain 32 pairs")
    logs = [math.log(candidate / legacy) for legacy, candidate in pairs]
    mean = sum(logs) / len(logs)
    variance = sum((value - mean) ** 2 for value in logs) / (len(logs) - 1)
    half_width = T31_95 * math.sqrt(variance / len(logs))
    return {
        "geomean": math.exp(mean),
        "lower95": math.exp(mean - half_width),
        "upper95": math.exp(mean + half_width),
    }


def require_cell_identity(
        fields: Mapping[str, str], cell: Tuple[int, int, int, int]) -> None:
    ordinal, K, block_bytes, repetitions = cell
    expected = {
        "version": "1", "ordinal": str(ordinal), "K": str(K),
        "block_bytes": str(block_bytes),
    }
    for key, value in expected.items():
        if fields.get(key) != value:
            raise ScreenError("cell identity mismatch in field " + key)
    if "repetitions" in fields and fields["repetitions"] != str(repetitions):
        raise ScreenError("cell repetition identity mismatch")


def parse_cell_output(
        stdout: bytes, returncode: int, cell: Tuple[int, int, int, int],
        expected_commit: str) -> Dict[str, Any]:
    lines = decode_output(stdout, "cell stdout")
    expected_line_count = 1 + len(PANELS) * (SAMPLES + 1) + 2
    if len(lines) != expected_line_count:
        raise ScreenError("cell emitted the wrong line count")
    cursor = 0
    meta = parse_csv_line(lines[cursor], "wh2_fixed5_meta", META_FIELDS)
    cursor += 1
    require_cell_identity(meta, cell)
    fixed_meta = {
        "version": "1", "source_commit": expected_commit,
        "rows": "256", "samples": "32", "mix_count": "3",
        "peel_seed": "4d241359", "active_degree": "2",
        "precode_seed": "786f72667573696f",
        "fallback_degree": "3", "active_terms": "5",
        "fallback_terms": "6", "wide_build": "1", "avx2": "1",
        "avx512": "1",
        "candidate": "direct_no_write_fixed5_try",
        "direct_wide_available": "1", "timed_row_generation": "0",
        "timed_pointer_gather": "0", "timed_counters": "0",
        "panel_order": "active_fallback_aa_split",
    }
    for key, value in fixed_meta.items():
        if meta[key] != value:
            raise ScreenError("metadata mismatch in field " + key)
    if parse_uint(meta["P"], "P", positive=True) != EXPECTED_P[cell[1]]:
        raise ScreenError("certified precode width mismatch")
    for key in ("active_first_id", "active_last_id",
                "fallback_first_id", "fallback_last_id"):
        if parse_uint(meta[key], key) >= 1000000:
            raise ScreenError(key + " is outside the frozen scan")
    for key in ("active_rows_hash", "fallback_rows_hash", "input_hash"):
        require_hash(meta[key], LOWER16, key)
    parsed_panels: Dict[str, Dict[str, Any]] = {}
    for panel_name, candidate_name, degree in PANELS:
        pairs: List[Tuple[int, int]] = []
        for sample in range(SAMPLES):
            fields = parse_csv_line(
                lines[cursor], "wh2_fixed5_raw", RAW_FIELDS)
            cursor += 1
            require_cell_identity(fields, cell)
            expected_order = "AB" if sample % 2 == 0 else "BA"
            expected_values = {
                "cell": panel_name, "candidate": candidate_name,
                "degree": str(degree), "rows": "256",
                "sample": str(sample), "order": expected_order,
            }
            for key, value in expected_values.items():
                if fields[key] != value:
                    raise ScreenError(
                        "%s raw mismatch in field %s" % (panel_name, key))
            legacy = parse_uint(fields["legacy_ns"], "legacy_ns", positive=True)
            candidate = parse_uint(
                fields["candidate_ns"], "candidate_ns", positive=True)
            ratio = parse_float(fields["ratio"], "ratio")
            if not close_enough(ratio, candidate / legacy):
                raise ScreenError(panel_name + " raw ratio mismatch")
            pairs.append((legacy, candidate))
        emitted = parse_csv_line(
            lines[cursor], "wh2_fixed5_summary", SUMMARY_FIELDS)
        cursor += 1
        require_cell_identity(emitted, cell)
        expected_summary = {
            "cell": panel_name, "candidate": candidate_name,
            "degree": str(degree), "samples": "32",
        }
        for key, value in expected_summary.items():
            if emitted[key] != value:
                raise ScreenError(
                    "%s summary mismatch in field %s" % (panel_name, key))
        stats = compute_stats(pairs)
        for key in ("geomean", "lower95", "upper95"):
            if not close_enough(parse_float(emitted[key], key), stats[key]):
                raise ScreenError(panel_name + " independently computed CI mismatch")
        parsed_panels[panel_name] = {
            "pairs": [[legacy, candidate] for legacy, candidate in pairs],
            **stats,
        }

    target = parse_csv_line(
        lines[cursor], "wh2_fixed5_target", TARGET_FIELDS)
    cursor += 1
    require_cell_identity(target, cell)
    expected_target = {
        "source_commit": expected_commit,
        "target_cpu": str(TARGET_CPU),
        "full_apic_id": TARGET_FULL_APIC_ID,
        "pre_identity_sha256": TARGET_IDENTITY_SHA256,
        "post_identity_sha256": TARGET_IDENTITY_SHA256,
        "canonical_bytes": str(TARGET_IDENTITY_CANONICAL_BYTES),
        "pre_before_cpu": str(TARGET_CPU),
        "pre_after_cpu": str(TARGET_CPU),
        "post_before_cpu": str(TARGET_CPU),
        "post_after_cpu": str(TARGET_CPU),
        "pre_before_affinity_count": "1",
        "pre_after_affinity_count": "1",
        "post_before_affinity_count": "1",
        "post_after_affinity_count": "1",
        "pre_voluntary_delta": "0", "pre_involuntary_delta": "0",
        "post_voluntary_delta": "0", "post_involuntary_delta": "0",
        "gate": "pass",
    }
    for key, value in expected_target.items():
        if target[key] != value:
            raise ScreenError("target receipt mismatch in field " + key)

    result = parse_csv_line(lines[cursor], "wh2_fixed5_result", RESULT_FIELDS)
    cursor += 1
    if cursor != len(lines):
        raise ScreenError("trailing cell output")
    require_cell_identity(result, cell)
    primary_pass = parsed_panels["active"]["upper95"] < 0.99
    fallback_pass = (
        parsed_panels["fallback"]["lower95"] >= 0.99 and
        parsed_panels["fallback"]["upper95"] <= 1.01
    )
    aa_pass = (
        parsed_panels["aa"]["lower95"] >= 0.99 and
        parsed_panels["aa"]["upper95"] <= 1.01 and
        parsed_panels["aa"]["lower95"] <= 1.0 <=
        parsed_panels["aa"]["upper95"]
    )
    controls_pass = fallback_pass and aa_pass
    expected_result = {
        "primary": "pass" if primary_pass else "reject",
        "fallback_control": "pass" if fallback_pass else "invalid",
        "aa_control": "pass" if aa_pass else "invalid",
        "controls": "pass" if controls_pass else "invalid",
        "split_promotional": "0", "mismatch_sink": "0",
        "status": ("invalid" if not controls_pass else
                   ("pass" if primary_pass else "reject")),
    }
    for key, value in expected_result.items():
        if result[key] != value:
            raise ScreenError("result mismatch in field " + key)
    expected_exit = 0 if controls_pass and primary_pass else \
        (10 if controls_pass else 11)
    if returncode != expected_exit or returncode not in (0, 10):
        raise ScreenError("cell exit is not a scientifically accepted 0/10")
    return {
        "meta": dict(meta), "panels": parsed_panels, "target": dict(target),
        "outcome": "pass" if primary_pass else "reject",
    }


def attempt_path(output_dir: Path) -> Path:
    return output_dir.with_name(output_dir.name + ".ATTEMPT")


def write_exclusive(path: Path, data: bytes, mode: int = 0o444) -> None:
    descriptor = os.open(
        str(path), os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
        mode,
    )
    try:
        offset = 0
        while offset < len(data):
            written = os.write(descriptor, data[offset:])
            if written <= 0:
                raise ScreenError("short evidence write")
            offset += written
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def fsync_directory(path: Path) -> None:
    descriptor = os.open(str(path), os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def claim_attempt(
        config: Config, preflight_value: Preflight,
        inputs: Dict[str, Any]) -> Path:
    path = attempt_path(config.output_dir)
    inputs_bytes = canonical_bytes(inputs)
    value = {
        "schema": "wirehair.wh2.fixed5-packet-screen.attempt.v1",
        "binary_sha256": preflight_value.binary_sha256,
        "controller_sha256": preflight_value.controller_sha256,
        "source_commit": config.expected_commit,
        "output_dir": str(config.output_dir.absolute()),
        "inputs": inputs,
        "inputs_sha256": sha256_bytes(inputs_bytes),
    }
    write_exclusive(path, canonical_bytes(value))
    fsync_directory(path.parent)
    return path


def rename_noreplace(source: Path, destination: Path) -> None:
    libc = ctypes.CDLL(None, use_errno=True)
    function = getattr(libc, "renameat2", None)
    if function is None:
        raise ScreenError("renameat2(RENAME_NOREPLACE) is unavailable")
    function.argtypes = [
        ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p,
        ctypes.c_uint,
    ]
    function.restype = ctypes.c_int
    result = function(
        -100, os.fsencode(source), -100, os.fsencode(destination), 1,
    )
    if result != 0:
        error = ctypes.get_errno()
        raise ScreenError(
            "atomic no-replace publication failed: " + os.strerror(error))


def self_hash(schema: str, value: Dict[str, Any], field: str) -> Dict[str, Any]:
    result = dict(value)
    result["schema"] = schema
    result[field] = sha256_bytes(canonical_bytes(result))
    return result


def publish(
        config: Config, inputs: Dict[str, Any], records: List[Dict[str, Any]],
        raw_outputs: Sequence[Tuple[str, bytes]], summary: Dict[str, Any]) -> None:
    parent = config.output_dir.parent
    staging = Path(tempfile.mkdtemp(
        prefix="." + config.output_dir.name + ".staging.", dir=str(parent)))
    try:
        inputs_bytes = canonical_bytes(inputs)
        results_bytes = b"".join(canonical_bytes(record) for record in records)
        summary_value = self_hash(
            "wirehair.wh2.fixed5-packet-screen.summary.v1",
            summary, "summary_sha256")
        summary_bytes = canonical_bytes(summary_value)
        write_exclusive(staging / "inputs.json", inputs_bytes)
        write_exclusive(staging / "results.jsonl", results_bytes)
        write_exclusive(staging / "summary.json", summary_bytes)
        raw_hashes: Dict[str, str] = {}
        for name, data in raw_outputs:
            if name in raw_hashes:
                raise ScreenError("duplicate raw artifact name")
            write_exclusive(staging / name, data)
            raw_hashes[name] = sha256_bytes(data)
        complete = self_hash(
            "wirehair.wh2.fixed5-packet-screen.complete.v1", {
                "inputs_sha256": sha256_bytes(inputs_bytes),
                "results_sha256": sha256_bytes(results_bytes),
                "summary_sha256": sha256_bytes(summary_bytes),
                "raw_artifacts": raw_hashes,
                "status": summary_value["status"],
            }, "receipt_sha256")
        write_exclusive(staging / "COMPLETE", canonical_bytes(complete))
        fsync_directory(staging)
        rename_noreplace(staging, config.output_dir)
        fsync_directory(parent)
    except BaseException:
        if staging.exists():
            shutil.rmtree(staging)
        raise


def make_inputs(config: Config, preflight_value: Preflight) -> Dict[str, Any]:
    commands = []
    for ordinal, K, block_bytes, _repetitions in FROZEN_CELLS:
        commands.append([
            str(config.taskset), "-c", str(config.cpu),
            str(config.benchmark), "--run-cell", str(ordinal), str(K),
            str(block_bytes), config.expected_commit, str(config.cpu),
        ])
    return {
        "schema": "wirehair.wh2.fixed5-packet-screen.inputs.v1",
        "source_commit": config.expected_commit,
        "binary": {
            "path": str(config.benchmark),
            "sha256": preflight_value.binary_sha256,
            "nm_stdout_sha256": preflight_value.nm_sha256,
        },
        "controller_sha256": preflight_value.controller_sha256,
        "taskset": {
            "path": str(config.taskset),
            "sha256": preflight_value.taskset_sha256,
            "cpu": config.cpu,
        },
        "describe_sha256": preflight_value.describe_sha256,
        "describe": preflight_value.describe_text,
        "cells": [list(cell) for cell in FROZEN_CELLS],
        "commands": commands,
        "samples": SAMPLES, "rows": ROWS, "student_t_df31": T31_95,
        "target_identity": {
            "cpu": TARGET_CPU, "full_apic_id": TARGET_FULL_APIC_ID,
            "canonical_bytes": TARGET_IDENTITY_CANONICAL_BYTES,
            "sha256": TARGET_IDENTITY_SHA256,
        },
        "cell_timeout_seconds": CELL_TIMEOUT_SECONDS,
    }


def run_campaign(config: Config) -> int:
    preflight_value = preflight(config)
    inputs = make_inputs(config, preflight_value)
    claim = claim_attempt(config, preflight_value, inputs)
    records: List[Dict[str, Any]] = []
    raw_outputs: List[Tuple[str, bytes]] = []
    invalid_reason: Optional[str] = None
    parsed_cells: List[Dict[str, Any]] = []
    for cell in FROZEN_CELLS:
        ordinal, K, block_bytes, _repetitions = cell
        command = inputs["commands"][ordinal]
        stdout = b""
        stderr = b""
        returncode: Optional[int] = None
        error_text: Optional[str] = None
        try:
            require_clean_commit(config)
            if sha256_file(config.benchmark) != preflight_value.binary_sha256:
                raise ScreenError("binary changed before cell")
            if sha256_file(config.taskset) != preflight_value.taskset_sha256:
                raise ScreenError("taskset changed before cell")
            result = subprocess.run(
                command, cwd=str(config.source_root), env=dict(RUN_ENV),
                stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, timeout=CELL_TIMEOUT_SECONDS,
                check=False,
            )
            stdout, stderr, returncode = (
                result.stdout, result.stderr, result.returncode)
            require_clean_commit(config)
            if sha256_file(config.benchmark) != preflight_value.binary_sha256:
                raise ScreenError("binary changed after cell")
            if sha256_file(config.taskset) != preflight_value.taskset_sha256:
                raise ScreenError("taskset changed after cell")
            if len(stdout) > MAX_OUTPUT_BYTES or len(stderr) > MAX_OUTPUT_BYTES:
                raise ScreenError("cell output exceeded limit")
            if stderr:
                raise ScreenError("cell wrote stderr")
            if returncode not in (0, 10):
                raise ScreenError("cell exit %d is not accepted" % returncode)
            parsed = parse_cell_output(
                stdout, returncode, cell, config.expected_commit)
            parsed_cells.append(parsed)
        except subprocess.TimeoutExpired as error:
            stdout = error.stdout or b""
            stderr = error.stderr or b""
            oversized: List[str] = []
            if len(stdout) > MAX_OUTPUT_BYTES:
                stdout = stdout[:MAX_OUTPUT_BYTES]
                oversized.append("stdout")
            if len(stderr) > MAX_OUTPUT_BYTES:
                stderr = stderr[:MAX_OUTPUT_BYTES]
                oversized.append("stderr")
            error_text = "cell timed out"
            if oversized:
                error_text += "; truncated oversized " + "/".join(oversized)
        except (OSError, ScreenError) as error:
            error_text = str(error)
        stdout_name = "cell-%02d.stdout" % ordinal
        stderr_name = "cell-%02d.stderr" % ordinal
        raw_outputs.extend(((stdout_name, stdout), (stderr_name, stderr)))
        record: Dict[str, Any] = {
            "schema": "wirehair.wh2.fixed5-packet-screen.result.v1",
            "ordinal": ordinal, "K": K, "block_bytes": block_bytes,
            "argv": command, "returncode": returncode,
            "stdout_file": stdout_name, "stdout_sha256": sha256_bytes(stdout),
            "stderr_file": stderr_name, "stderr_sha256": sha256_bytes(stderr),
            "error": error_text,
        }
        if error_text is None:
            record["parsed"] = parsed_cells[-1]
        records.append(record)
        if error_text is not None:
            invalid_reason = "cell-%02d: %s" % (ordinal, error_text)
            break

    if invalid_reason is None:
        by_K: Dict[int, Tuple[str, str, int, str]] = {}
        for parsed in parsed_cells:
            meta = parsed["meta"]
            K = int(meta["K"])
            identity = (
                meta["active_rows_hash"], meta["fallback_rows_hash"],
                int(meta["P"]), meta["avx512"],
            )
            if K in by_K and by_K[K] != identity:
                invalid_reason = "same-K row/P/feature identity changed across B"
                break
            by_K[K] = identity
        feature_rows = {
            (cell["meta"]["wide_build"], cell["meta"]["avx2"],
             cell["meta"]["avx512"])
            for cell in parsed_cells
        }
        if invalid_reason is None and len(feature_rows) != 1:
            invalid_reason = "CPU feature stratum changed across cells"
    if invalid_reason is not None:
        status = "invalid"
    elif all(cell["outcome"] == "pass" for cell in parsed_cells):
        status = "pass"
    else:
        status = "reject"
    summary = {
        "status": status,
        "invalid_reason": invalid_reason,
        "completed_cells": len(records),
        "required_cells": len(FROZEN_CELLS),
        "outcomes": [cell["outcome"] for cell in parsed_cells],
        "features": ([
            parsed_cells[0]["meta"]["wide_build"],
            parsed_cells[0]["meta"]["avx2"],
            parsed_cells[0]["meta"]["avx512"],
        ] if parsed_cells else None),
        "cpu": config.cpu,
        "source_commit": config.expected_commit,
        "binary_sha256": preflight_value.binary_sha256,
    }
    publish(config, inputs, records, raw_outputs, summary)
    try:
        claim.unlink()
        fsync_directory(claim.parent)
    except OSError:
        # The published no-replace directory is already the durable consumed
        # namespace.  A retained claim is conservative and still forbids rerun.
        pass
    return 0 if status == "pass" else (10 if status == "reject" else 1)


def parse_args(argv: Optional[Sequence[str]] = None) -> Config:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--expected-source-commit", required=True)
    parser.add_argument("--binary-sha256", required=True)
    parser.add_argument("--cpu", required=True, type=int)
    parser.add_argument("--taskset", type=Path, default=Path("/usr/bin/taskset"))
    parser.add_argument("--nm", type=Path, default=Path("/usr/bin/nm"))
    parser.add_argument("--git", type=Path, default=Path("/usr/bin/git"))
    arguments = parser.parse_args(argv)
    return Config(
        benchmark=arguments.benchmark.resolve(),
        output_dir=arguments.output_dir.absolute(),
        source_root=arguments.source_root.resolve(),
        expected_commit=arguments.expected_source_commit,
        expected_binary_sha256=arguments.binary_sha256,
        cpu=arguments.cpu,
        taskset=arguments.taskset.resolve(), nm=arguments.nm.resolve(),
        git=arguments.git.resolve(),
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    try:
        return run_campaign(parse_args(argv))
    except ScreenError as error:
        print("fixed5 packet screen: " + str(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
