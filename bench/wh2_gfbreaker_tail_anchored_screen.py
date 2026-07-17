#!/usr/bin/env python3
"""Paired screen for the experiment-only tail-anchored GF(256) R2 schedule.

The screen is intentionally limited to the already revealed, sealed R2
confirmation artifacts.  It runs the current and nested binaries at R1 and
R2 on every one of the nine hard strata for all 30 revealed discordance K,
plus an equally sized deterministic set of clean controls.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import shutil
import stat
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation, getcontext
from pathlib import Path
from typing import Iterable, Sequence


getcontext().prec = 50

SCHEMA = "wirehair.wh2.gfbreaker_tail_anchored_screen.v1"
DEFAULT_REVEALED_ROOT = Path(
    "/tmp/wh2-gfbreaker-r2-confirm-20260717T142410Z-09d44f1-retry2"
)
DEFAULT_CURRENT_BINARY = Path(
    "/tmp/wirehair-gf256-r2-confirmation/"
    "build-r2-native/codec/wirehair_v2_bench"
)
EXPECTED_CURRENT_BINARY_SHA256 = (
    "7e5ff7ca4fc276b6b10d4ed3c1a9b1256646543ec4b7bded2b5c4a1597b2da9c"
)
EXPECTED_INPUT_SHA256 = {
    "analysis_complete.sha256":
        "0970e2dc1fd1adda652a94ab7c059d3a53947e6af96307f04a0517582c3ad7b7",
    "cells.csv":
        "50cd14e55eb1e7696639f3795e3bf34322983dabea5e06f8215db393c8e9e267",
    "exact_k_changes.csv":
        "d4d9efc78cf3baafa620db6ed15398197f34c5800206d9b108f23b372dd40699",
}

SCHEDULES = ("burst", "adversarial", "repair-only")
SEEDS = (
    "0xa6b527b9ae8de8d7",
    "0x75446e1619e81d8a",
    "0x007e9dd892a319f5",
)
CSV_HEADER = (
    "N", "bb", "heavy_family", "mix_count", "overhead", "trials",
    "success", "rank_fail", "error", "fail_rate", "inact_mu",
    "inact_max", "binary_def_mu", "binary_def_max", "heavy_gain_mu",
    "heavy_gain_min", "heavy_shortfall", "solve_ms_mu", "build_ms_mu",
    "peel_ms_mu", "project_ms_mu", "residual_ms_mu", "backsub_ms_mu",
    "seed_attempt", "block_xors_mu", "block_muladds_mu",
    "first_rank_fail", "binary_def_hist", "heavy_gain_hist",
    "failure_trials", "active_packet_peel_seed_xor",
)
NON_TIMING_INDICES = tuple(i for i in range(len(CSV_HEADER)) if i not in range(17, 23))
SEALED_CELL_HEADER = (
    "K", "salt", "schedule", "seed_index", "seed", "r0_failed",
    "r1_failed", "r2_failed", "r0_error", "r1_error", "r2_error",
    "r0_heavy_shortfall", "r1_heavy_shortfall", "r2_heavy_shortfall",
    "r0_block_xors_milli", "r1_block_xors_milli", "r2_block_xors_milli",
    "r0_block_muladds_milli", "r1_block_muladds_milli",
    "r2_block_muladds_milli", "r2_source",
)
EXACT_K_HEADER = (
    "K", "cells", "r0_fail", "r1_fail", "r2_fail", "r2_minus_r0",
    "r2_minus_r1", "r2_intro_r0", "r2_intro_r1", "r2_repair_r0",
    "r2_repair_r1",
)
ARMS = (
    ("current_r1", "current", 1),
    ("nested_r1", "nested", 1),
    ("current_r2", "current", 2),
    ("nested_r2", "nested", 2),
)
UINT_RE = re.compile(r"^(0|[1-9][0-9]*)$")
INT_RE = re.compile(r"^(0|-?[1-9][0-9]*)$")
DIAGNOSTIC_RE = re.compile(
    r"^# precodefail_gf256_breaker: N=(0|[1-9][0-9]*) "
    r"base_hash_seed=(0x[0-9a-f]+) breaker_rows=([12]) "
    r"breaker_hash_seed0=(0x[0-9a-f]+) seed_xor=0xb7e15162 "
    r"heavy_tail=A$"
)


class ScreenError(RuntimeError):
    pass


def die(message: str) -> None:
    raise ScreenError(message)


def strict_uint(text: str, context: str) -> int:
    if not UINT_RE.fullmatch(text):
        die(f"{context}: invalid canonical unsigned integer {text!r}")
    return int(text)


def strict_int(text: str, context: str) -> int:
    if not INT_RE.fullmatch(text):
        die(f"{context}: invalid canonical integer {text!r}")
    return int(text)


def strict_decimal(text: str, context: str) -> Decimal:
    try:
        value = Decimal(text)
    except InvalidOperation:
        die(f"{context}: invalid decimal {text!r}")
    if not value.is_finite() or value < 0:
        die(f"{context}: non-finite/negative decimal {text!r}")
    return value


def canonical_u64(text: str, context: str) -> int:
    try:
        value = int(text, 0)
    except ValueError:
        die(f"{context}: invalid integer {text!r}")
    if value < 0 or value > (1 << 64) - 1:
        die(f"{context}: integer outside uint64")
    return value


def regular_file(path: Path, context: str) -> Path:
    if path.is_symlink():
        die(f"{context}: symlink is not accepted: {path}")
    try:
        status = path.stat()
    except OSError as exc:
        die(f"{context}: cannot stat {path}: {exc}")
    if not stat.S_ISREG(status.st_mode):
        die(f"{context}: not a regular file: {path}")
    return path.resolve()


def sha256(path: Path) -> str:
    resolved = regular_file(path, "SHA256 input")
    before = resolved.stat()
    digest = hashlib.sha256()
    with resolved.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    after = resolved.stat()
    identity = lambda value: (
        value.st_dev, value.st_ino, value.st_size,
        value.st_mtime_ns, value.st_ctime_ns,
    )
    if identity(before) != identity(after):
        die(f"SHA256 input changed while being read: {resolved}")
    return digest.hexdigest()


def file_identity(path: Path) -> tuple[int, int, int, int, int]:
    status = regular_file(path, "file identity").stat()
    return (
        status.st_dev, status.st_ino, status.st_size,
        status.st_mtime_ns, status.st_ctime_ns,
    )


def atomic_bytes(path: Path, value: bytes) -> None:
    temporary = path.with_name(path.name + ".partial")
    if temporary.exists() or temporary.is_symlink():
        die(f"refusing stale temporary output {temporary}")
    with temporary.open("xb") as stream:
        stream.write(value)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def atomic_text(path: Path, value: str) -> None:
    atomic_bytes(path, value.encode("utf-8"))


def atomic_json(path: Path, value: object) -> None:
    atomic_text(path, json.dumps(
        value, sort_keys=True, indent=2, allow_nan=False,
        ensure_ascii=True,
    ) + "\n")


def parse_preamble(line: str, context: str) -> dict[str, str]:
    prefix = "# precodefail: "
    if not line.startswith(prefix):
        die(f"{context}: missing preamble")
    result: dict[str, str] = {}
    for token in line[len(prefix):].split():
        if token.count("=") != 1:
            die(f"{context}: malformed preamble token {token!r}")
        key, value = token.split("=", 1)
        if not key or key in result:
            die(f"{context}: duplicate/empty preamble key {key!r}")
        result[key] = value
    return result


def parse_histogram(text: str, context: str) -> tuple[int, int]:
    if not text:
        die(f"{context}: empty histogram")
    count_sum = 0
    weighted_sum = 0
    previous = -1
    for item in text.split("|"):
        if item.count(":") != 1:
            die(f"{context}: malformed histogram")
        raw_value, raw_count = item.split(":", 1)
        value = strict_uint(raw_value, context)
        count = strict_uint(raw_count, context)
        if count == 0 or value <= previous:
            die(f"{context}: zero count/nonascending histogram")
        previous = value
        count_sum += count
        weighted_sum += value * count
    return count_sum, weighted_sum


@dataclass(frozen=True)
class SealedCell:
    values: dict[str, str]

    @property
    def key(self) -> tuple[int, str, int]:
        return (
            int(self.values["K"]), self.values["schedule"],
            int(self.values["seed_index"]),
        )


@dataclass(frozen=True)
class Task:
    task_id: int
    kind: str
    k: int
    schedule: str
    seed_index: int
    seed: str
    salt: str
    arm: str
    binary_name: str
    binary: Path
    rows: int
    cpu: int = -1

    @property
    def stem(self) -> str:
        return f"task{self.task_id:04d}"


@dataclass(frozen=True)
class ParsedOutput:
    preamble: dict[str, str]
    diagnostic: tuple[str, str]
    row: tuple[str, ...]

    @property
    def failed(self) -> bool:
        return self.row[7] != "0" or self.row[8] != "0"

    @property
    def error(self) -> int:
        return int(self.row[8])

    @property
    def heavy_shortfall(self) -> int:
        return int(self.row[16])

    @property
    def block_xors(self) -> Decimal:
        return Decimal(self.row[24])

    @property
    def block_muladds(self) -> Decimal:
        return Decimal(self.row[25])

    @property
    def semantic(self) -> tuple[str, ...]:
        return tuple(self.row[index] for index in NON_TIMING_INDICES)


def verify_revealed_inputs(root: Path) -> tuple[Path, Path, dict[str, str]]:
    analysis = root.resolve() / "analysis/r2"
    observed: dict[str, str] = {}
    for name, expected in EXPECTED_INPUT_SHA256.items():
        path = regular_file(analysis / name, f"revealed {name}")
        digest = sha256(path)
        if digest != expected:
            die(f"revealed {name} SHA256 {digest}, want {expected}")
        observed[name] = digest
    seal_path = analysis / "analysis_complete.sha256"
    sealed: dict[str, str] = {}
    for number, line in enumerate(seal_path.read_text(encoding="utf-8").splitlines(), 1):
        match = re.fullmatch(r"([0-9a-f]{64})  (/.+)", line)
        if not match:
            die(f"revealed analysis seal line {number}: malformed")
        name = Path(match.group(2)).name
        if name in sealed:
            die(f"revealed analysis seal duplicates {name}")
        sealed[name] = match.group(1)
    if set(sealed) != {
        "cells.csv", "exact_k_changes.csv", "scopes.csv",
        "summary.json", "summary.txt",
    }:
        die("revealed analysis seal does not name the exact five artifacts")
    for name in ("cells.csv", "exact_k_changes.csv"):
        if sealed[name] != observed[name]:
            die(f"revealed analysis seal/hash mismatch for {name}")
    return analysis / "cells.csv", analysis / "exact_k_changes.csv", observed


def load_changed_ks(path: Path) -> tuple[int, ...]:
    before = file_identity(path)
    result: list[int] = []
    totals = {name: 0 for name in EXACT_K_HEADER[2:]}
    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream)
        try:
            header = tuple(next(reader))
        except StopIteration:
            die("exact K changes is empty")
        if header != EXACT_K_HEADER:
            die("exact K changes header mismatch")
        for number, values in enumerate(reader, 2):
            if len(values) != len(EXACT_K_HEADER):
                die(f"exact K changes:{number}: field count")
            parsed = [
                (strict_int if index in (5, 6) else strict_uint)(
                    value, f"exact K changes:{number}:{EXACT_K_HEADER[index]}",
                )
                for index, value in enumerate(values)
            ]
            k = parsed[0]
            if not 2 <= k <= 64000 or parsed[1] != 9:
                die(f"exact K changes:{number}: K/cell invariant")
            if result and k <= result[-1]:
                die("exact K changes is not strictly sorted")
            result.append(k)
            for name, value in zip(EXACT_K_HEADER[2:], parsed[2:]):
                totals[name] += value
    expected = {
        "r0_fail": 11, "r1_fail": 8, "r2_fail": 11,
        "r2_minus_r0": 0, "r2_minus_r1": 3,
        "r2_intro_r0": 11, "r2_intro_r1": 11,
        "r2_repair_r0": 11, "r2_repair_r1": 8,
    }
    if file_identity(path) != before:
        die("exact K changes mutated while being parsed")
    if len(result) != 30 or totals != expected:
        die(f"revealed exact K ledger geometry/totals mismatch: {len(result)}, {totals}")
    return tuple(result)


def parse_sealed_cell(values: Sequence[str], number: int) -> SealedCell:
    if len(values) != len(SEALED_CELL_HEADER):
        die(f"sealed cells:{number}: field count")
    row = dict(zip(SEALED_CELL_HEADER, values))
    k = strict_uint(row["K"], f"sealed cells:{number}:K")
    seed_index = strict_uint(row["seed_index"], f"sealed cells:{number}:seed")
    if not 2 <= k <= 64000 or row["schedule"] not in SCHEDULES or seed_index >= 3:
        die(f"sealed cells:{number}: key invariant")
    if row["seed"] != SEEDS[seed_index]:
        die(f"sealed cells:{number}: seed value mismatch")
    canonical_u64(row["salt"], f"sealed cells:{number}:salt")
    for field in (
        "r0_failed", "r1_failed", "r2_failed", "r0_error", "r1_error",
        "r2_error", "r0_heavy_shortfall", "r1_heavy_shortfall",
        "r2_heavy_shortfall", "r0_block_xors_milli",
        "r1_block_xors_milli", "r2_block_xors_milli",
        "r0_block_muladds_milli", "r1_block_muladds_milli",
        "r2_block_muladds_milli",
    ):
        strict_uint(row[field], f"sealed cells:{number}:{field}")
    return SealedCell(row)


def read_sealed_selection(
    cells_path: Path,
    changed_ks: Sequence[int],
) -> tuple[tuple[int, ...], dict[tuple[int, str, int], SealedCell]]:
    before = file_identity(cells_path)
    changed = set(changed_ks)
    salt0_count = [0] * 64001
    clean = [True] * 64001
    changed_rows: dict[tuple[int, str, int], SealedCell] = {}
    total_rows = 0
    salt0_rows = 0
    with cells_path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream)
        if tuple(next(reader, ())) != SEALED_CELL_HEADER:
            die("sealed cells header mismatch")
        for number, values in enumerate(reader, 2):
            cell = parse_sealed_cell(values, number)
            total_rows += 1
            if canonical_u64(cell.values["salt"], "sealed salt") != 0:
                continue
            salt0_rows += 1
            k, schedule, seed_index = cell.key
            salt0_count[k] += 1
            if any(cell.values[field] != "0" for field in (
                "r0_failed", "r1_failed", "r2_failed", "r0_error",
                "r1_error", "r2_error", "r0_heavy_shortfall",
                "r1_heavy_shortfall", "r2_heavy_shortfall",
            )):
                clean[k] = False
            if k in changed:
                key = (k, schedule, seed_index)
                if key in changed_rows:
                    die(f"sealed cells duplicate changed key {key}")
                changed_rows[key] = cell
    if total_rows != 575991 or salt0_rows != 575361:
        die(f"sealed cells cardinality {total_rows}/{salt0_rows} mismatch")
    salt0_geometry = {
        count: sum(salt0_count[k] == count for k in range(2, 64001))
        for count in set(salt0_count[2:])
    }
    if salt0_geometry != {0: 70, 9: 63929}:
        die(f"sealed salt-0 K geometry mismatch: {salt0_geometry}")
    expected_changed = {
        (k, schedule, seed_index)
        for k in changed for schedule in SCHEDULES for seed_index in range(3)
    }
    if set(changed_rows) != expected_changed:
        die("sealed cells do not cover exact changed-K strata")

    eligible = {
        k for k in range(2, 64001)
        if clean[k] and salt0_count[k] == 9 and k not in changed
    }
    controls: list[int] = []
    used: set[int] = set()
    for target in changed_ks:
        candidates = sorted(eligible - used, key=lambda k: (abs(k - target), k))
        if not candidates:
            die(f"cannot select a clean control for K={target}")
        controls.append(candidates[0])
        used.add(candidates[0])
    if len(controls) != 30 or len(set(controls)) != 30:
        die("clean control selection cardinality mismatch")

    selected_ks = changed | set(controls)
    selected = dict(changed_rows)
    with cells_path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream)
        if tuple(next(reader, ())) != SEALED_CELL_HEADER:
            die("sealed cells changed before second pass")
        for number, values in enumerate(reader, 2):
            # The canonical K field is cheap to inspect before full validation.
            if not values or not UINT_RE.fullmatch(values[0]):
                die(f"sealed cells:{number}: malformed K on second pass")
            if int(values[0]) not in selected_ks:
                continue
            cell = parse_sealed_cell(values, number)
            if canonical_u64(cell.values["salt"], "selected salt") != 0:
                continue
            selected[cell.key] = cell
    expected_selected = {
        (k, schedule, seed_index)
        for k in selected_ks for schedule in SCHEDULES for seed_index in range(3)
    }
    if set(selected) != expected_selected or len(selected) != 540:
        die("sealed selected-cell set is not exact 60 K x nine strata")
    if file_identity(cells_path) != before:
        die("sealed cells mutated while being parsed")
    return tuple(controls), selected


def make_command(task: Task) -> list[str]:
    return [
        str(task.binary), "precodefail", "--N", str(task.k),
        "--bb-list", "64", "--seed-block-bytes", "1280",
        "--overhead", "0", "--trials", "1", "--threads", "1",
        "--loss", "0.50", "--seed", task.seed,
        "--schedule", task.schedule, "--completion", "mixed",
        "--mix-count", "2", "--packet-peel-seed-xor", "0",
        "--mixed-gf256-rows", "11", "--mixed-gf16-rows", "4",
        "--mixed-period", "32", "--mixed-geometry", "shared-x",
        "--mixed-residue-schedule", "hashed",
        "--mixed-residue-hash-seed", "68", "--mixed-residue-hash-keyed",
        "--mixed-independent-extension-residues",
        "--mixed-extension-residue-seed-xor", "78",
        "--mixed-independent-gf256-breaker-rows", str(task.rows),
        "--mixed-fused-schedule-buckets", "auto",
    ]


def parse_output(data: bytes, task: Task) -> ParsedOutput:
    context = task.stem
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        die(f"{context}: stdout is not UTF-8: {exc}")
    lines = text.splitlines()
    if len(lines) != 4:
        die(f"{context}: expected exactly four output lines, got {len(lines)}")
    preamble = parse_preamble(lines[0], context)
    expected_preamble = {
        "trials": "1", "threads": "1", "loss": "0.5",
        "seed": hex(canonical_u64(task.seed, context)),
        "completion": "mixed", "mixed_period": "32",
        "mixed_gf256_rows": "11", "mixed_gf16_rows": "4",
        "mixed_geometry": "shared-x", "mixed_residue_skew": "0",
        "mixed_residue_schedule": "hashed",
        "mixed_residue_hash_seed": "0x44",
        "mixed_residue_hash_keyed": "1",
        "mixed_independent_extension_residues": "1",
        "mixed_extension_residue_seed_xor": "0x4e",
        "mixed_fused_schedule_buckets_mode": "-1",
        "source_hits_override": "0", "packet_peel_seed_xor": "0x0",
        "packet_peel_seed_table": "none", "binary_dense_rows_override": "0",
        "gf256_heavy_rows_override": "0", "odd_packet_peel_seed_xor": "0x0",
        "packet_row_seed_multiplier": "0x1",
        "packet_row_seed_avalanche": "0",
        "seed_block_bytes_override": "1280", "overhead_stream": "salted",
        "full_payload_solve": "0", "schedule": task.schedule,
        "mixed_independent_gf256_breaker_residues": "1",
        "mixed_independent_gf256_breaker_rows": str(task.rows),
        "mixed_gf256_breaker_residue_seed_xor": "0xb7e15162",
    }
    if preamble != expected_preamble:
        die(f"{context}: exact preamble mismatch")
    if tuple(next(csv.reader([lines[1]]))) != CSV_HEADER:
        die(f"{context}: CSV header mismatch")
    diagnostic = DIAGNOSTIC_RE.fullmatch(lines[2])
    if not diagnostic or int(diagnostic.group(1)) != task.k or int(diagnostic.group(3)) != task.rows:
        die(f"{context}: diagnostic mismatch")
    row = tuple(next(csv.reader([lines[3]])))
    if len(row) != len(CSV_HEADER):
        die(f"{context}: CSV row field count")
    if (row[0], row[1], row[2], row[3], row[4], row[5]) != (
        str(task.k), "64", "periodic", "2", "0", "1",
    ):
        die(f"{context}: fixed row geometry mismatch")
    for index in (1, 3, 4, 5, 6, 7, 8, 11, 13, 15, 16, 23):
        strict_uint(row[index], f"{context}:{CSV_HEADER[index]}")
    for index in (9, 10, 12, 14, 17, 18, 19, 20, 21, 22, 24, 25):
        strict_decimal(row[index], f"{context}:{CSV_HEADER[index]}")
    success, rank_fail, error = map(int, row[6:9])
    if success + rank_fail + error != 1 or int(row[16]) > rank_fail:
        die(f"{context}: outcome count mismatch")
    if row[9] != ("1.00000000" if rank_fail + error else "0.00000000"):
        die(f"{context}: fail-rate mismatch")
    binary_count, binary_sum = parse_histogram(row[27], f"{context}:binary histogram")
    heavy_count, heavy_sum = parse_histogram(row[28], f"{context}:heavy histogram")
    if binary_count != 1 or heavy_count != 1:
        die(f"{context}: histogram trial count")
    if Decimal(row[12]) != binary_sum or int(row[13]) != binary_sum:
        die(f"{context}: binary histogram summary")
    if Decimal(row[14]) != heavy_sum:
        die(f"{context}: heavy histogram summary")
    if (rank_fail == 0 and row[26] != "-1") or (rank_fail == 1 and row[26] != "0"):
        die(f"{context}: first-rank-fail mismatch")
    failures = () if row[29] == "" else tuple(row[29].split("|"))
    if len(failures) != rank_fail + error or any(item != "0" for item in failures):
        die(f"{context}: failure-trial ledger mismatch")
    if canonical_u64(row[30], context) != 0:
        die(f"{context}: active salt is not zero")
    return ParsedOutput(
        preamble, (diagnostic.group(2), diagnostic.group(4)), row,
    )


def milli(value: Decimal, context: str) -> int:
    scaled = value * 1000
    if scaled != scaled.to_integral_value():
        die(f"{context}: work value is not exact milli-unit")
    return int(scaled)


def check_sealed_reproduction(
    cell: SealedCell,
    output: ParsedOutput,
    prefix: str,
    context: str,
) -> None:
    expected_failed = cell.values[f"{prefix}_failed"] == "1"
    expected_error = int(cell.values[f"{prefix}_error"])
    expected_shortfall = int(cell.values[f"{prefix}_heavy_shortfall"])
    expected_xors = int(cell.values[f"{prefix}_block_xors_milli"])
    expected_muladds = int(cell.values[f"{prefix}_block_muladds_milli"])
    observed = (
        output.failed, output.error, output.heavy_shortfall,
        milli(output.block_xors, context),
        milli(output.block_muladds, context),
    )
    expected = (
        expected_failed, expected_error, expected_shortfall,
        expected_xors, expected_muladds,
    )
    if observed != expected:
        die(f"{context}: rerun does not reproduce sealed {prefix}: {observed} != {expected}")


def run_task(
    task: Task,
    raw: Path,
    taskset: str,
    timeout: float,
    environment: dict[str, str],
) -> ParsedOutput:
    command = [taskset, "-c", str(task.cpu), *make_command(task)]
    try:
        completed = subprocess.run(
            command, cwd=Path(__file__).resolve().parents[1],
            env=environment, capture_output=True, check=False,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired as exc:
        stdout = exc.stdout or b""
        stderr = exc.stderr or b""
        atomic_bytes(raw / f"{task.stem}.stdout", stdout)
        atomic_bytes(raw / f"{task.stem}.stderr", stderr)
        die(f"{task.stem}: timed out after {timeout} seconds")
    atomic_bytes(raw / f"{task.stem}.stdout", completed.stdout)
    atomic_bytes(raw / f"{task.stem}.stderr", completed.stderr)
    if completed.returncode != 0:
        die(f"{task.stem}: exit status {completed.returncode}")
    if completed.stderr:
        die(f"{task.stem}: non-empty stderr")
    return parse_output(completed.stdout, task)


def assign_cpus(tasks: Sequence[Task], cpus: Sequence[int]) -> tuple[list[Task], list[list[Task]]]:
    bins: list[list[Task]] = [[] for _ in cpus]
    loads = [0] * len(cpus)
    assigned: dict[int, Task] = {}
    for task in sorted(tasks, key=lambda value: (-value.k, value.task_id)):
        slot = min(range(len(cpus)), key=lambda index: (loads[index], index))
        bound = replace(task, cpu=cpus[slot])
        bins[slot].append(bound)
        assigned[task.task_id] = bound
        loads[slot] += task.k + 1000
    return [assigned[index] for index in range(len(tasks))], bins


def run_screen_tasks(
    tasks: Sequence[Task],
    bins: Sequence[Sequence[Task]],
    raw: Path,
    taskset: str,
    timeout: float,
) -> dict[int, ParsedOutput]:
    environment = dict(os.environ)
    for key in list(environment):
        if key.startswith("WIREHAIR_") or key in {
            "ASAN_OPTIONS", "UBSAN_OPTIONS", "LSAN_OPTIONS", "TSAN_OPTIONS",
        }:
            environment.pop(key, None)
    environment["LC_ALL"] = "C"
    environment["TZ"] = "UTC"

    def run_bin(batch: Sequence[Task]) -> list[tuple[int, ParsedOutput]]:
        return [
            (task.task_id, run_task(task, raw, taskset, timeout, environment))
            for task in batch
        ]

    outputs: dict[int, ParsedOutput] = {}
    with ThreadPoolExecutor(max_workers=len(bins)) as executor:
        futures = [executor.submit(run_bin, batch) for batch in bins if batch]
        for future in as_completed(futures):
            for task_id, output in future.result():
                if task_id in outputs:
                    die(f"duplicate task result {task_id}")
                outputs[task_id] = output
    if set(outputs) != {task.task_id for task in tasks}:
        die("task result set is incomplete")
    return outputs


def decimal_ratio(numerator: Decimal, denominator: Decimal) -> str | None:
    if denominator == 0:
        return None
    return str(numerator / denominator)


def write_csv(path: Path, rows: Iterable[Sequence[object]]) -> None:
    temporary = path.with_name(path.name + ".partial")
    if temporary.exists() or temporary.is_symlink():
        die(f"refusing stale temporary output {temporary}")
    with temporary.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerows(rows)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def summarize_scope(records: Sequence[dict[str, object]]) -> dict[str, object]:
    totals: dict[str, object] = {
        "cells": len(records),
        "r0_fail": sum(bool(row["r0_failed"]) for row in records),
        "r1_fail": sum(bool(row["r1_failed"]) for row in records),
    }
    for arm in ("current_r2", "nested_r2"):
        arm_fail = [bool(row[f"{arm}_failed"]) for row in records]
        r0_fail = [bool(row["r0_failed"]) for row in records]
        r1_fail = [bool(row["r1_failed"]) for row in records]
        totals[f"{arm}_fail"] = sum(arm_fail)
        totals[f"{arm}_repair_r0"] = sum(a and not b for a, b in zip(r0_fail, arm_fail))
        totals[f"{arm}_intro_r0"] = sum(not a and b for a, b in zip(r0_fail, arm_fail))
        totals[f"{arm}_repair_r1"] = sum(a and not b for a, b in zip(r1_fail, arm_fail))
        totals[f"{arm}_intro_r1"] = sum(not a and b for a, b in zip(r1_fail, arm_fail))
    totals["nested_repairs_current"] = sum(
        bool(row["current_r2_failed"]) and not bool(row["nested_r2_failed"])
        for row in records
    )
    totals["nested_intros_current"] = sum(
        not bool(row["current_r2_failed"]) and bool(row["nested_r2_failed"])
        for row in records
    )
    common = [
        row for row in records
        if not bool(row["current_r2_failed"]) and not bool(row["nested_r2_failed"])
    ]
    totals["common_success"] = len(common)
    for metric in ("block_xors", "block_muladds"):
        current_all = sum((Decimal(str(row[f"current_r2_{metric}"])) for row in records), Decimal(0))
        nested_all = sum((Decimal(str(row[f"nested_r2_{metric}"])) for row in records), Decimal(0))
        current_common = sum((Decimal(str(row[f"current_r2_{metric}"])) for row in common), Decimal(0))
        nested_common = sum((Decimal(str(row[f"nested_r2_{metric}"])) for row in common), Decimal(0))
        totals[f"nested_over_current_{metric}_all"] = decimal_ratio(nested_all, current_all)
        totals[f"nested_over_current_{metric}_common_success"] = decimal_ratio(
            nested_common, current_common,
        )
    return totals


def derive_results(
    out: Path,
    tasks: Sequence[Task],
    outputs: dict[int, ParsedOutput],
    selected: dict[tuple[int, str, int], SealedCell],
    changed_ks: Sequence[int],
    control_ks: Sequence[int],
) -> dict[str, object]:
    by_cell: dict[tuple[int, str, int], dict[str, ParsedOutput]] = {}
    for task in tasks:
        key = (task.k, task.schedule, task.seed_index)
        by_cell.setdefault(key, {})[task.arm] = outputs[task.task_id]
    if any(set(arms) != {arm[0] for arm in ARMS} for arms in by_cell.values()):
        die("per-cell arm set is incomplete")

    records: list[dict[str, object]] = []
    identity_mismatches = 0
    changed_set = set(changed_ks)
    for key in sorted(by_cell, key=lambda value: (value[0], SCHEDULES.index(value[1]), value[2])):
        k, schedule, seed_index = key
        sealed = selected[key]
        arms = by_cell[key]
        current_r1 = arms["current_r1"]
        nested_r1 = arms["nested_r1"]
        current_r2 = arms["current_r2"]
        nested_r2 = arms["nested_r2"]
        context = f"K={k}/{schedule}/seed{seed_index}"
        check_sealed_reproduction(sealed, current_r1, "r1", context + "/R1")
        check_sealed_reproduction(sealed, current_r2, "r2", context + "/R2")
        r1_identity = (
            current_r1.preamble == nested_r1.preamble and
            current_r1.diagnostic == nested_r1.diagnostic and
            current_r1.semantic == nested_r1.semantic
        )
        if not r1_identity:
            identity_mismatches += 1
        record: dict[str, object] = {
            "kind": "hard" if k in changed_set else "control",
            "K": k, "schedule": schedule, "seed_index": seed_index,
            "seed": SEEDS[seed_index], "salt": "0x0",
            "r0_failed": sealed.values["r0_failed"] == "1",
            "r1_failed": sealed.values["r1_failed"] == "1",
            "current_r2_failed": current_r2.failed,
            "nested_r2_failed": nested_r2.failed,
            "current_r2_error": current_r2.error,
            "nested_r2_error": nested_r2.error,
            "current_r2_heavy_shortfall": current_r2.heavy_shortfall,
            "nested_r2_heavy_shortfall": nested_r2.heavy_shortfall,
            "current_r2_block_xors": str(current_r2.block_xors),
            "nested_r2_block_xors": str(nested_r2.block_xors),
            "current_r2_block_muladds": str(current_r2.block_muladds),
            "nested_r2_block_muladds": str(nested_r2.block_muladds),
            "current_r2_breaker_seed0": current_r2.diagnostic[1],
            "nested_r2_breaker_seed0": nested_r2.diagnostic[1],
            "r1_exact_identity": r1_identity,
        }
        records.append(record)

    cell_header = list(records[0])
    write_csv(out / "cells.csv", [cell_header] + [
        [int(value) if isinstance(value, bool) else value for value in (row[name] for name in cell_header)]
        for row in records
    ])
    changes = [
        row for row in records
        if bool(row["current_r2_failed"]) != bool(row["nested_r2_failed"])
    ]
    change_header = [
        "kind", "K", "schedule", "seed_index", "r0_failed", "r1_failed",
        "current_r2_failed", "nested_r2_failed", "current_r2_heavy_shortfall",
        "nested_r2_heavy_shortfall",
    ]
    write_csv(out / "outcome_changes.csv", [change_header] + [
        [int(value) if isinstance(value, bool) else value for value in (row[name] for name in change_header)]
        for row in changes
    ])

    hard = [row for row in records if row["kind"] == "hard"]
    controls = [row for row in records if row["kind"] == "control"]
    hard_summary = summarize_scope(hard)
    control_summary = summarize_scope(controls)
    gates = {
        "exact_540_cell_geometry": len(records) == 540 and len(hard) == 270 and len(controls) == 270,
        "r1_exact_identity": identity_mismatches == 0,
        "zero_runtime_errors": all(
            int(row["current_r2_error"]) == 0 and int(row["nested_r2_error"]) == 0
            for row in records
        ),
        "current_hard_reproduces_11_repairs_11_introductions":
            hard_summary["current_r2_repair_r0"] == 11 and
            hard_summary["current_r2_intro_r0"] == 11,
        "nested_retains_all_11_r0_repairs": hard_summary["nested_r2_repair_r0"] == 11,
        "nested_reduces_r0_introductions": hard_summary["nested_r2_intro_r0"] < 11,
        "nested_controls_clean": control_summary["nested_r2_fail"] == 0,
    }
    summary: dict[str, object] = {
        "schema": SCHEMA,
        "hard": hard_summary,
        "controls": control_summary,
        "r1_identity_mismatches": identity_mismatches,
        "current_to_nested_outcome_changes": len(changes),
        "gates": gates,
        "promotion_screen_passed": all(gates.values()),
        "timing_claim_valid": False,
        "timing_note": "Parallel saturated screen; structural work only, no timing claim.",
        "changed_ks": list(changed_ks),
        "control_ks": list(control_ks),
    }
    atomic_json(out / "summary.json", summary)
    lines = [
        f"promotion_screen_passed={int(bool(summary['promotion_screen_passed']))}",
        f"r1_identity_mismatches={identity_mismatches}",
        f"hard_current_fail={hard_summary['current_r2_fail']}",
        f"hard_nested_fail={hard_summary['nested_r2_fail']}",
        f"hard_current_repairs_r0={hard_summary['current_r2_repair_r0']}",
        f"hard_current_intros_r0={hard_summary['current_r2_intro_r0']}",
        f"hard_nested_repairs_r0={hard_summary['nested_r2_repair_r0']}",
        f"hard_nested_intros_r0={hard_summary['nested_r2_intro_r0']}",
        f"hard_nested_repairs_current={hard_summary['nested_repairs_current']}",
        f"hard_nested_intros_current={hard_summary['nested_intros_current']}",
        f"control_nested_fail={control_summary['nested_r2_fail']}",
        f"hard_xor_nested_over_current_common={hard_summary['nested_over_current_block_xors_common_success']}",
        f"hard_muladd_nested_over_current_common={hard_summary['nested_over_current_block_muladds_common_success']}",
        "timing_claim_valid=0",
    ]
    atomic_text(out / "summary.txt", "\n".join(lines) + "\n")
    return summary


def seal_output(out: Path) -> str:
    seal = out / "artifact.sha256"
    paths = sorted(
        path for path in out.rglob("*")
        if path.is_file() and path != seal and path.name != "INCOMPLETE"
    )
    rows: list[str] = []
    for path in paths:
        if path.is_symlink():
            die(f"output seal refuses symlink {path}")
        rows.append(f"{sha256(path)}  {path.relative_to(out)}\n")
    atomic_text(seal, "".join(rows))
    for number, line in enumerate(seal.read_text(encoding="utf-8").splitlines(), 1):
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
        if not match:
            die(f"output seal:{number}: malformed")
        path = out / match.group(2)
        if sha256(path) != match.group(1):
            die(f"output seal:{number}: verification failed")
    return sha256(seal)


def git_provenance(worktree: Path) -> dict[str, str]:
    def capture(*args: str) -> bytes:
        return subprocess.check_output(
            ["git", "-C", str(worktree), *args], stderr=subprocess.DEVNULL,
        )
    try:
        head = capture("rev-parse", "HEAD").decode("ascii").strip()
        diff = capture("diff", "--binary", "--", "codec", "bench")
        status = capture("status", "--short", "--", "codec", "bench").decode("utf-8")
    except (OSError, subprocess.CalledProcessError, UnicodeError) as exc:
        die(f"cannot capture candidate Git provenance: {exc}")
    return {
        "head": head,
        "codec_bench_diff_sha256": hashlib.sha256(diff).hexdigest(),
        "codec_bench_status": status,
    }


def self_test() -> None:
    parsed = parse_preamble("# precodefail: a=1 b=2", "self-test")
    if parsed != {"a": "1", "b": "2"}:
        die("self-test preamble parse")
    try:
        parse_preamble("# precodefail: a=1 a=2", "self-test")
    except ScreenError:
        pass
    else:
        die("self-test duplicate preamble acceptance")
    if decimal_ratio(Decimal(3), Decimal(2)) != "1.5":
        die("self-test decimal ratio")
    if parse_histogram("2:1", "self-test") != (1, 2):
        die("self-test histogram")
    print("self-test: PASS")


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--revealed-root", type=Path, default=DEFAULT_REVEALED_ROOT)
    parser.add_argument("--current-binary", type=Path, default=DEFAULT_CURRENT_BINARY)
    parser.add_argument(
        "--nested-binary", type=Path,
        default=Path(__file__).resolve().parents[1] /
            "build-tail/codec/wirehair_v2_bench",
    )
    parser.add_argument("--out", type=Path)
    parser.add_argument("--workers", type=int, default=120)
    parser.add_argument("--timeout-seconds", type=float, default=120.0)
    parser.add_argument("--self-test", action="store_true")
    return parser.parse_args(argv)


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)
    if args.self_test:
        self_test()
        return 0
    if args.out is None:
        die("--out is required unless --self-test is used")
    if args.workers <= 0 or args.timeout_seconds <= 0:
        die("workers and timeout must be positive")
    taskset = shutil.which("taskset")
    if not taskset:
        die("taskset is required for one-worker-per-CPU scheduling")
    allowed = sorted(os.sched_getaffinity(0))
    cpus = allowed[:min(args.workers, len(allowed))]
    if not cpus:
        die("empty CPU affinity set")

    cells_path, exact_path, input_hashes = verify_revealed_inputs(args.revealed_root)
    changed_ks = load_changed_ks(exact_path)
    controls, selected = read_sealed_selection(cells_path, changed_ks)
    current_binary = regular_file(args.current_binary, "current binary")
    nested_binary = regular_file(args.nested_binary, "nested binary")
    current_sha = sha256(current_binary)
    nested_sha = sha256(nested_binary)
    binary_identities = {
        "current": file_identity(current_binary),
        "nested": file_identity(nested_binary),
    }
    script_path = Path(__file__).resolve()
    script_sha = sha256(script_path)
    script_identity = file_identity(script_path)
    if current_sha != EXPECTED_CURRENT_BINARY_SHA256:
        die(f"current binary SHA256 {current_sha}, want {EXPECTED_CURRENT_BINARY_SHA256}")
    if current_sha == nested_sha:
        die("current and nested binaries are byte-identical")

    out = args.out.resolve()
    if out.exists() or out.is_symlink():
        die(f"output already exists: {out}")
    out.mkdir(parents=True)
    raw = out / "raw"
    raw.mkdir()
    atomic_text(out / "INCOMPLETE", "tail-anchored screen in progress\n")

    tasks: list[Task] = []
    binaries = {"current": current_binary, "nested": nested_binary}
    kinds = (("hard", changed_ks), ("control", controls))
    for kind, ks in kinds:
        for k in sorted(ks):
            for schedule in SCHEDULES:
                for seed_index, seed in enumerate(SEEDS):
                    for arm, binary_name, rows in ARMS:
                        tasks.append(Task(
                            len(tasks), kind, k, schedule, seed_index, seed,
                            "0x0", arm, binary_name, binaries[binary_name], rows,
                        ))
    if len(tasks) != 2160:
        die(f"task ledger has {len(tasks)} entries, want 2160")
    tasks, bins = assign_cpus(tasks, cpus)

    selection_rows: list[list[object]] = [[
        "kind", "K", "schedule", "seed_index", "seed", "salt",
        "sealed_r0_failed", "sealed_r1_failed", "sealed_current_r2_failed",
    ]]
    for kind, ks in kinds:
        for k in sorted(ks):
            for schedule in SCHEDULES:
                for seed_index, seed in enumerate(SEEDS):
                    cell = selected[(k, schedule, seed_index)]
                    selection_rows.append([
                        kind, k, schedule, seed_index, seed, "0x0",
                        cell.values["r0_failed"], cell.values["r1_failed"],
                        cell.values["r2_failed"],
                    ])
    write_csv(out / "selection.csv", selection_rows)
    write_csv(out / "jobs.csv", [[
        "task", "kind", "K", "schedule", "seed_index", "seed", "salt",
        "arm", "binary", "rows", "cpu", "stdout", "stderr", "command",
    ]] + [[
        task.task_id, task.kind, task.k, task.schedule, task.seed_index,
        task.seed, task.salt, task.arm, task.binary_name, task.rows, task.cpu,
        f"raw/{task.stem}.stdout", f"raw/{task.stem}.stderr",
        json.dumps(make_command(task), separators=(",", ":")),
    ] for task in tasks])
    run_metadata = {
        "schema": SCHEMA,
        "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "script": str(script_path),
        "script_sha256": script_sha,
        "revealed_root": str(args.revealed_root.resolve()),
        "revealed_input_sha256": input_hashes,
        "current_binary": str(current_binary),
        "current_binary_sha256": current_sha,
        "nested_binary": str(nested_binary),
        "nested_binary_sha256": nested_sha,
        "candidate_git": git_provenance(Path(__file__).resolve().parents[1]),
        "workers": len(cpus), "cpus": cpus, "tasks": len(tasks),
        "hard_cells": 270, "control_cells": 270,
        "control_selection": "nearest unique all-R0/R1/R2-clean salt0 K; tie to lower K",
        "timing_claim_valid": False,
        "environment_policy": "inherit, strip WIREHAIR_* and sanitizer options, LC_ALL=C,TZ=UTC",
    }
    atomic_json(out / "run.json", run_metadata)

    outputs = run_screen_tasks(
        tasks, bins, raw, taskset, args.timeout_seconds,
    )
    if (file_identity(current_binary) != binary_identities["current"] or
        file_identity(nested_binary) != binary_identities["nested"] or
        sha256(current_binary) != current_sha or
        sha256(nested_binary) != nested_sha):
        die("a benchmark binary mutated during the screen")
    if file_identity(script_path) != script_identity or sha256(script_path) != script_sha:
        die("the screen harness mutated during the screen")
    summary = derive_results(
        out, tasks, outputs, selected, changed_ks, controls,
    )
    (out / "INCOMPLETE").unlink()
    seal_digest = seal_output(out)
    print(json.dumps({
        "output": str(out),
        "artifact_seal_sha256": seal_digest,
        "promotion_screen_passed": summary["promotion_screen_passed"],
        "hard": summary["hard"], "controls": summary["controls"],
        "r1_identity_mismatches": summary["r1_identity_mismatches"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main(sys.argv[1:]))
    except ScreenError as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2)
