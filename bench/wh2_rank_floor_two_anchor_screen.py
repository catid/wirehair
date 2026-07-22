#!/usr/bin/env python3
"""Shared parsing, statistics, and thermal helpers for WH2 campaigns.

The one-off standalone rank-floor screen that originally lived in this module
is retired: it predates the frozen-executable and process-group guarantees of
the supported controllers and cannot produce supported evidence.  The helper
APIs remain because the hardened all-K, phase, and preferred-attempt
controllers freeze and import this module.
"""

from __future__ import annotations

import sys

# This module is imported from exact-inventory frozen campaign directories.
sys.dont_write_bytecode = True

import csv
import errno
import hashlib
import io
import json
import math
import os
from pathlib import Path
import queue
import re
import stat
import time
from collections import defaultdict
from decimal import Decimal
from typing import Any, BinaryIO, Iterable, Sequence


SOURCE_CELLS_SHA256 = (
    "afbed8f0c6fdd49f1644fc0e12b827097cbc6752fddf6e6ce4484b118b93659f"
)
BASELINE_BINARY_SHA256 = (
    "7e5ff7ca4fc276b6b10d4ed3c1a9b1256646543ec4b7bded2b5c4a1597b2da9c"
)
SEEDS = (
    "0xa6b527b9ae8de8d7",
    "0x75446e1619e81d8a",
    "0x007e9dd892a319f5",
)
SCHEDULES = ("burst", "adversarial", "repair-only")
ARMS = ("d12", "two_anchor_adaptive", "d13_adaptive", "d14")
CUTOFF = 4096
CUTOFF_NEIGHBORS = (4094, 4095, 4096, 4097, 4098)
IDENTITY_KS = (2, 64, 945, 4095, 4096, 4097, 8192, 32000, 64000)
TIMING_FIELDS = {
    "solve_ms_mu", "build_ms_mu", "peel_ms_mu", "project_ms_mu",
    "residual_ms_mu", "backsub_ms_mu",
}
THERMAL_FIELDS = (
    "utc", "monotonic_s", "cpu_busy_pct", "cpu_avg_mhz", "cpu_tctl_c",
    "dimm_i2c1_50_c", "dimm_i2c1_51_c", "dimm_i2c1_52_c",
    "dimm_i2c1_53_c", "dimm_i2c2_50_c", "dimm_i2c2_51_c",
    "dimm_i2c2_52_c", "dimm_i2c2_53_c", "dimm_read_errors",
    "load1", "load5", "load15", "edac_ce", "edac_ue",
)
THERMAL_DIMM_FIELDS = THERMAL_FIELDS[5:13]
THERMAL_ROW_MAX_BYTES = 16384
THERMAL_SCALAR_PATTERN = re.compile(r"[0-9]+(?:\.[0-9]+)?\Z")
THERMAL_COUNTER_PATTERN = re.compile(r"[0-9]+\Z")
BENCH_HEADER_V1 = (
    "N", "bb", "heavy_family", "mix_count", "overhead", "trials",
    "success", "rank_fail", "error", "fail_rate", "inact_mu",
    "inact_max", "binary_def_mu", "binary_def_max", "heavy_gain_mu",
    "heavy_gain_min", "heavy_shortfall", "solve_ms_mu", "build_ms_mu",
    "peel_ms_mu", "project_ms_mu", "residual_ms_mu", "backsub_ms_mu",
    "seed_attempt", "block_xors_mu", "block_muladds_mu",
    "first_rank_fail", "binary_def_hist", "heavy_gain_hist",
    "failure_trials", "active_packet_peel_seed_xor",
)
BENCH_HEADER_JOINT_DELTA = BENCH_HEADER_V1 + (
    "mixed_joint_source_xors_mu", "mixed_joint_marginal_xors_mu",
    "mixed_joint_marginal_copies_mu", "mixed_joint_active_deltas_mu",
    "mixed_joint_scratch_bytes_mu", "mixed_dual_source_columns_mu",
)
BASE_OPTIONS = (
    "precodefail",
    "--bb-list", "64",
    "--overhead", "0",
    "--trials", "1",
    "--threads", "1",
    "--loss", "0.50",
    "--completion", "mixed",
    "--mix-count", "2",
)

STANDALONE_RETIREMENT_MESSAGE = (
    "the standalone WH2 rank-floor screen is retired and cannot produce "
    "supported evidence; use wh2_rank_floor_two_anchor_allk.py prepare, then "
    "execute the exact frozen run_command recorded in prepare.json"
)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path, *, require_unique: bool = True) -> str:
    flags = (os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
             getattr(os, "O_NOFOLLOW", 0) |
             getattr(os, "O_NONBLOCK", 0))
    try:
        descriptor = os.open(path, flags)
    except OSError as exc:
        if exc.errno == errno.ELOOP:
            raise ValueError("refusing symlink input: %s" % path) from exc
        raise ValueError("unable to open unique regular input: %s" % path) from exc
    try:
        before = os.fstat(descriptor)
        if (not stat.S_ISREG(before.st_mode) or
                (require_unique and before.st_nlink != 1)):
            raise ValueError("refusing nonunique regular input: %s" % path)
        with os.fdopen(descriptor, "rb", closefd=False) as source:
            digest = hashlib.sha256()
            for block in iter(lambda: source.read(1 << 20), b""):
                digest.update(block)
            midpoint = os.fstat(descriptor)
            pathname_midpoint = os.stat(path, follow_symlinks=False)
            source.seek(0)
            confirmation = hashlib.sha256()
            for block in iter(lambda: source.read(1 << 20), b""):
                confirmation.update(block)
        after = os.fstat(descriptor)
        pathname_after = os.stat(path, follow_symlinks=False)
    finally:
        os.close(descriptor)
    identity = lambda value: (
        value.st_dev, value.st_ino, value.st_mode, value.st_nlink,
        value.st_size, value.st_mtime_ns, value.st_ctime_ns,
    )
    if (identity(before) != identity(midpoint) or
            not stat.S_ISREG(pathname_midpoint.st_mode) or
            (require_unique and pathname_midpoint.st_nlink != 1) or
            identity(pathname_midpoint) != identity(midpoint) or
            digest.digest() != confirmation.digest() or
            identity(midpoint) != identity(after) or
            not stat.S_ISREG(pathname_after.st_mode) or
            (require_unique and pathname_after.st_nlink != 1) or
            identity(pathname_after) != identity(after)):
        raise ValueError("input changed while hashing: %s" % path)
    return digest.hexdigest()


def stable_file_bytes(path: Path) -> bytes:
    flags = (os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
             getattr(os, "O_NOFOLLOW", 0) |
             getattr(os, "O_NONBLOCK", 0))
    try:
        descriptor = os.open(path, flags)
    except OSError as exc:
        if exc.errno == errno.ELOOP:
            raise ValueError("refusing symlink input: %s" % path) from exc
        raise ValueError("unable to open unique regular input: %s" % path) from exc
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode) or before.st_nlink != 1:
            raise ValueError("refusing nonunique regular input: %s" % path)
        with os.fdopen(descriptor, "rb", closefd=False) as source:
            data = source.read()
            midpoint = os.fstat(descriptor)
            pathname_midpoint = os.stat(path, follow_symlinks=False)
            source.seek(0)
            confirmation = source.read()
        after = os.fstat(descriptor)
        pathname_after = os.stat(path, follow_symlinks=False)
    finally:
        os.close(descriptor)
    identity = lambda value: (
        value.st_dev, value.st_ino, value.st_mode, value.st_nlink,
        value.st_size, value.st_mtime_ns, value.st_ctime_ns,
    )
    if (identity(before) != identity(midpoint) or
            len(data) != before.st_size or
            not stat.S_ISREG(pathname_midpoint.st_mode) or
            pathname_midpoint.st_nlink != 1 or
            identity(pathname_midpoint) != identity(midpoint) or
            data != confirmation or
            identity(midpoint) != identity(after) or
            not stat.S_ISREG(pathname_after.st_mode) or
            pathname_after.st_nlink != 1 or
            identity(pathname_after) != identity(after)):
        raise ValueError("input changed while reading: %s" % path)
    return data


def canonical_json(value: Any) -> str:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    )


def stratum_key(row: dict[str, str]) -> tuple[int, str, str]:
    return int(row["seed_index"]), row["seed"], row["schedule"]


def expected_strata() -> tuple[tuple[int, str, str], ...]:
    return tuple(
        (seed_index, seed, schedule)
        for seed_index, seed in enumerate(SEEDS)
        for schedule in SCHEDULES
    )


def load_source_cells(path: Path) -> dict[str, Any]:
    raw = stable_file_bytes(path)
    if sha256_bytes(raw) != SOURCE_CELLS_SHA256:
        raise ValueError("source cells SHA256 mismatch")
    required = {
        "K", "bb", "salt", "schedule", "loss", "seed_index", "seed",
        "base_rank_fail", "candidate_rank_fail", "base_error",
        "candidate_error", "base_failed", "candidate_failed",
        "base_heavy_shortfall", "base_inact_mu", "base_inact_max",
        "base_seed_attempt", "base_binary_def_mu", "base_binary_def_max",
        "base_binary_def_hist", "base_block_xors_mu",
        "base_block_muladds_mu",
    }
    states: dict[int, dict[str, Any]] = {}
    hard_cases: list[dict[str, str]] = []
    seen_cells: set[tuple[int, int, str]] = set()
    strata: set[tuple[int, str, str]] = set()
    row_count = 0
    try:
        source_text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError("source cells CSV is not UTF-8") from exc
    with io.StringIO(source_text, newline="") as source:
        reader = csv.DictReader(source)
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError("source cells CSV has an unexpected schema")
        for row in reader:
            row_count += 1
            K = int(row["K"])
            seed_index = int(row["seed_index"])
            schedule = row["schedule"]
            key = (K, seed_index, schedule)
            if key in seen_cells:
                raise ValueError(f"duplicate source cell {key}")
            seen_cells.add(key)
            if not 2 <= K <= 64000 or row["bb"] != "64" or row["loss"] != "0.50":
                raise ValueError(f"source cell geometry mismatch at {key}")
            if not 0 <= seed_index < len(SEEDS) or row["seed"] != SEEDS[seed_index]:
                raise ValueError(f"source seed mismatch at {key}")
            if schedule not in SCHEDULES:
                raise ValueError(f"source schedule mismatch at {key}")
            salt_value = int(row["salt"], 0)
            if not 0 <= salt_value <= 0xFFFFFFFF:
                raise ValueError(f"source salt does not fit u32 at {key}")
            strata.add(stratum_key(row))
            state = states.setdefault(K, {
                "salt": row["salt"], "strata": set(), "shared_success": True,
            })
            if state["salt"] != row["salt"]:
                raise ValueError(f"inconsistent source salt at K={K}")
            state["strata"].add((seed_index, schedule))
            base_failed = int(row["base_failed"])
            candidate_failed = int(row["candidate_failed"])
            base_error = int(row["base_error"])
            candidate_error = int(row["candidate_error"])
            if any(value not in (0, 1) for value in (
                base_failed, candidate_failed, base_error, candidate_error
            )):
                raise ValueError(f"non-binary source outcome at {key}")
            if base_error or candidate_error:
                raise ValueError(f"source error at {key}")
            if (base_failed != int(row["base_rank_fail"]) + base_error or
                    candidate_failed !=
                    int(row["candidate_rank_fail"]) + candidate_error):
                raise ValueError(f"source failed/rank/error mismatch at {key}")
            state["shared_success"] &= not (base_failed or candidate_failed)
            if base_failed or candidate_failed:
                hard_cases.append(dict(row))

    if row_count != 575991 or len(seen_cells) != 575991:
        raise ValueError(
            f"source cell cardinality {row_count}/{len(seen_cells)}, want 575991"
        )
    if set(states) != set(range(2, 64001)):
        raise ValueError("source K domain is not exactly 2..64000")
    if strata != set(expected_strata()):
        raise ValueError("source strata are not the expected three-by-three grid")
    if any(len(state["strata"]) != 9 for state in states.values()):
        raise ValueError("one or more source K values lack nine strata")
    if len(hard_cases) != 485:
        raise ValueError(f"expected 485 hard source cases, found {len(hard_cases)}")
    hard_ks = {int(row["K"]) for row in hard_cases}
    if len(hard_ks) != 474:
        raise ValueError(f"expected 474 unique hard K values, found {len(hard_ks)}")
    return {
        "states": states,
        "hard_cases": hard_cases,
        "hard_ks": hard_ks,
    }


def select_scopes(source: dict[str, Any], control_count: int) -> dict[str, set[int]]:
    hard_ks = set(source["hard_ks"])
    eligible = [
        K for K, state in source["states"].items()
        if state["shared_success"] and state["salt"] == "0x0" and
        K not in hard_ks and K not in CUTOFF_NEIGHBORS
    ]
    eligible.sort(key=lambda K: hashlib.sha256(
        f"wh2-rank-floor-control-v1|{K}".encode("ascii")
    ).digest())
    if len(eligible) < control_count:
        raise ValueError("not enough shared-success controls")
    controls = set(eligible[:control_count])
    scopes = {
        "revealed_hard_k": hard_ks,
        "shared_success_control_k": controls,
        "cutoff_neighbors": set(CUTOFF_NEIGHBORS),
    }
    if len(controls) != control_count or controls & hard_ks:
        raise ValueError("control selection collision/cardinality mismatch")
    return scopes


def arm_options(arm: str, band: str) -> tuple[str, ...]:
    if arm == "d12":
        return ()
    if arm == "two_anchor_adaptive":
        return ("--binary-dense-two-anchor",) if band == "large" else ()
    if arm == "d13_adaptive":
        return ("--binary-dense-rows", "13") if band == "large" else ()
    if arm == "d14":
        return ("--binary-dense-rows", "14")
    raise ValueError(f"unknown arm {arm}")


def make_command(
    binary: Path,
    Ks: Sequence[int],
    seed: str,
    schedule: str,
    salt: str,
    options: Sequence[str] = (),
) -> list[str]:
    return [
        str(binary), *BASE_OPTIONS,
        "--N", ",".join(map(str, Ks)),
        "--seed", seed,
        "--schedule", schedule,
        "--packet-peel-seed-xor", str(int(salt, 0)),
        *options,
    ]


def build_tasks(
    selected: set[int], states: dict[int, dict[str, Any]], chunk_size: int
) -> list[dict[str, Any]]:
    tasks: list[dict[str, Any]] = []
    for arm in ARMS:
        for seed_index, seed, schedule in expected_strata():
            for band, predicate in (
                ("small", lambda K: K < CUTOFF),
                ("large", lambda K: K >= CUTOFF),
            ):
                by_salt: dict[str, list[int]] = defaultdict(list)
                for K in sorted(selected):
                    if predicate(K):
                        by_salt[states[K]["salt"]].append(K)
                for salt in sorted(by_salt, key=lambda value: (int(value, 0), value)):
                    values = by_salt[salt]
                    for offset in range(0, len(values), chunk_size):
                        tasks.append({
                            "arm": arm,
                            "band": band,
                            "seed_index": seed_index,
                            "seed": seed,
                            "schedule": schedule,
                            "salt": salt,
                            "Ks": values[offset:offset + chunk_size],
                        })
    for job, task in enumerate(tasks):
        task["job"] = job
    keys = [
        (task["arm"], task["seed_index"], task["schedule"], K)
        for task in tasks for K in task["Ks"]
    ]
    expected = len(ARMS) * len(expected_strata()) * len(selected)
    if len(keys) != expected or len(keys) != len(set(keys)):
        raise ValueError("task ledger does not exactly cover selected arms/strata/K")
    return tasks


def parse_preamble(line: str) -> dict[str, str]:
    if not line.startswith("# precodefail:"):
        raise ValueError("missing precodefail preamble")
    result: dict[str, str] = {}
    for token in line.split()[2:]:
        if "=" not in token:
            raise ValueError(f"malformed preamble token {token!r}")
        key, value = token.split("=", 1)
        if key in result:
            raise ValueError(f"duplicate preamble token {key}")
        result[key] = value
    return result


def parse_bench_output(
    output: str,
    expected_Ks: Sequence[int],
    expected_salt: str,
    expected_seed: str,
    expected_schedule: str,
    expected_options: Sequence[str],
    candidate: bool,
    expected_trials: int = 1,
    expected_bb: int = 64,
) -> list[dict[str, str]]:
    def option_value(flag: str, default: str) -> str:
        try:
            index = expected_options.index(flag)
        except ValueError:
            return default
        if index + 1 >= len(expected_options):
            raise ValueError(f"benchmark option {flag} lacks a value")
        return expected_options[index + 1]

    lines = output.splitlines()
    if len(lines) != len(expected_Ks) + 2:
        raise ValueError(
            f"benchmark line count {len(lines)}, want {len(expected_Ks) + 2}"
        )
    preamble = parse_preamble(lines[0])
    expected_anchor = "1" if "--binary-dense-two-anchor" in expected_options else "0"
    if candidate and preamble.get("binary_dense_two_anchor") != expected_anchor:
        raise ValueError("binary_dense_two_anchor preamble mismatch")
    if not candidate and "binary_dense_two_anchor" in preamble:
        raise ValueError("baseline binary unexpectedly reports two-anchor hook")
    dense_rows = option_value("--binary-dense-rows", "0")
    hash_seed = int(option_value("--mixed-residue-hash-seed", "0"), 0)
    extension_xor = int(
        option_value("--mixed-extension-residue-seed-xor", "78"), 0
    )
    fixed = {
        "trials": str(expected_trials), "threads": "1", "loss": "0.5",
        "completion": option_value("--completion", "mixed"),
        "mixed_period": option_value("--mixed-period", "244"),
        "mixed_gf256_rows": option_value("--mixed-gf256-rows", "10"),
        "mixed_gf16_rows": option_value("--mixed-gf16-rows", "2"),
        "mixed_geometry": option_value("--mixed-geometry", "frozen"),
        "mixed_residue_skew": option_value("--mixed-residue-skew", "0"),
        "mixed_residue_schedule": option_value(
            "--mixed-residue-schedule", "constant"
        ),
        "mixed_residue_hash_seed": hex(hash_seed),
        "mixed_residue_hash_keyed": (
            "1" if "--mixed-residue-hash-keyed" in expected_options else "0"
        ),
        "mixed_independent_extension_residues": (
            "1" if "--mixed-independent-extension-residues" in expected_options
            else "0"
        ),
        "mixed_residue_buckets_requested": option_value(
            "--mixed-fused-schedule-buckets", "auto"
        ),
        "mixed_extension_residue_seed_xor": hex(extension_xor),
        "source_hits_override": option_value("--source-hits", "0"),
        "packet_peel_seed_table": option_value(
            "--packet-peel-seed-table", "none"
        ),
        "binary_dense_rows_override": dense_rows,
        "gf256_heavy_rows_override": option_value(
            "--gf256-heavy-rows", "0"
        ),
        "odd_packet_peel_seed_xor": hex(int(option_value(
            "--odd-packet-peel-seed-xor", "0"
        ), 0)),
        "packet_row_seed_multiplier": hex(int(option_value(
            "--packet-row-seed-multiplier", "1"
        ), 0)),
        "packet_row_seed_avalanche": (
            "1" if "--packet-row-seed-avalanche" in expected_options else "0"
        ),
        "seed_block_bytes_override": option_value("--seed-block-bytes", "0"),
        "overhead_stream": (
            "paired" if "--paired-overhead-stream" in expected_options
            else "salted"
        ),
        "full_payload_solve": (
            "1" if "--full-payload-solve" in expected_options else "0"
        ),
        "schedule": expected_schedule,
    }
    for key, expected in fixed.items():
        if preamble.get(key) != expected:
            raise ValueError(
                f"preamble {key}={preamble.get(key)!r}, want {expected!r}"
            )
    if int(preamble.get("packet_peel_seed_xor", "-1"), 0) != int(expected_salt, 0):
        raise ValueError("preamble salt mismatch")
    if int(preamble.get("seed", "-1"), 0) != int(expected_seed, 0):
        raise ValueError("preamble seed mismatch")
    header = next(csv.reader([lines[1]]))
    if tuple(header) not in (BENCH_HEADER_V1, BENCH_HEADER_JOINT_DELTA):
        raise ValueError("benchmark CSV header mismatch")
    records: list[dict[str, str]] = []
    order: list[int] = []
    for line in lines[2:]:
        values = next(csv.reader([line]))
        if len(values) != len(header):
            raise ValueError("benchmark CSV field count mismatch")
        row = dict(zip(header, values))
        K = int(row["N"])
        order.append(K)
        if (row["bb"], row["heavy_family"], row["mix_count"],
                row["overhead"], row["trials"]) != (
                    str(expected_bb), "periodic", "2", "0",
                    str(expected_trials)):
            raise ValueError(f"result geometry mismatch K={K}")
        success = int(row["success"])
        rank_fail = int(row["rank_fail"])
        error = int(row["error"])
        if success + rank_fail + error != expected_trials:
            raise ValueError(f"result outcome total mismatch K={K}")
        if int(row["active_packet_peel_seed_xor"], 0) != int(expected_salt, 0):
            raise ValueError(f"result salt mismatch K={K}")
        records.append(row)
    if order != list(expected_Ks):
        raise ValueError("benchmark K order/set mismatch")
    return records


class CpuPool:
    def __init__(self, cpus: Sequence[int], workers: int) -> None:
        self._queue: queue.Queue[int] = queue.Queue()
        for cpu in list(cpus)[:workers]:
            self._queue.put(cpu)

    def acquire(self) -> int:
        return self._queue.get()

    def release(self, cpu: int) -> None:
        self._queue.put(cpu)


def execute(
    command: Sequence[str], cpu_pool: CpuPool, timeout: float
) -> tuple[list[str], int, str, str, int, int]:
    del command, cpu_pool, timeout
    raise RuntimeError(STANDALONE_RETIREMENT_MESSAGE)


def run_identity(
    baseline: Path,
    candidate: Path,
    result_dir: Path,
    cpu_pool: CpuPool,
    timeout: float,
) -> dict[str, Any]:
    del baseline, candidate, result_dir, cpu_pool, timeout
    raise RuntimeError(STANDALONE_RETIREMENT_MESSAGE)


def run_task(
    task: dict[str, Any],
    binary: Path,
    result_dir: Path,
    cpu_pool: CpuPool,
    timeout: float,
) -> dict[str, Any]:
    del task, binary, result_dir, cpu_pool, timeout
    raise RuntimeError(STANDALONE_RETIREMENT_MESSAGE)


def binomial_one_sided(repairs: int, introductions: int) -> dict[str, Any]:
    if repairs < 0 or introductions < 0:
        raise ValueError("paired outcome counts must be nonnegative")
    discordant = repairs + introductions
    if discordant == 0:
        return {"numerator": 1, "denominator": 1, "decimal": "1"}

    def coefficient_sum(begin: int, end: int) -> int:
        if begin >= end:
            return 0
        term = math.comb(discordant, begin)
        total = term
        for k in range(begin + 1, end):
            term = term * (discordant - k + 1) // k
            total += term
        return total

    denominator = 1 << discordant
    center = discordant // 2 + (discordant & 1)
    choices = (
        (discordant - repairs + 1, "upper"),
        (repairs, "lower"),
        (abs(repairs - center), "center"),
    )
    representation = min(choices)[1]
    if representation == "upper":
        numerator = coefficient_sum(repairs, discordant + 1)
    elif representation == "lower":
        numerator = denominator - coefficient_sum(0, repairs)
    else:
        numerator = denominator >> 1
        if (discordant & 1) == 0:
            numerator += math.comb(discordant, center) >> 1
        if repairs < center:
            numerator += coefficient_sum(repairs, center)
        elif repairs > center:
            numerator -= coefficient_sum(center, repairs)
    divisor = math.gcd(numerator, denominator)
    numerator //= divisor
    denominator //= divisor
    return {
        "numerator": numerator, "denominator": denominator,
        "decimal": str(Decimal(numerator) / Decimal(denominator)),
    }


def failed(row: dict[str, str]) -> bool:
    return int(row["rank_fail"]) != 0 or int(row["error"]) != 0


def deterministic_equal(a: dict[str, str], b: dict[str, str]) -> bool:
    if a.keys() != b.keys():
        return False
    fields = set(a)
    ignored = TIMING_FIELDS | {"arm", "band"}
    return all(a[field] == b[field] for field in fields - ignored)


def summarize_scope(
    Ks: set[int],
    indexed: dict[tuple[str, int, str, int], dict[str, str]],
) -> dict[str, Any]:
    keys = [
        (seed_index, schedule, K)
        for seed_index, _, schedule in expected_strata() for K in sorted(Ks)
    ]
    arms: dict[str, Any] = {}
    for arm in ARMS:
        rows = [indexed[(arm, seed_index, schedule, K)] for seed_index, schedule, K in keys]
        success_rows = [row for row in rows if not failed(row)]
        failure_rows = [row for row in rows if failed(row)]
        arms[arm] = {
            "cells": len(rows),
            "successful_cells": len(success_rows),
            "failure_path_cells": len(failure_rows),
            "failures": len(failure_rows),
            "errors": sum(int(row["error"]) for row in rows),
            "heavy_shortfalls": sum(int(row["heavy_shortfall"]) for row in rows),
            **{
                f"{metric}_{outcome}_sum": str(sum(
                    Decimal(row[field]) for row in selected_rows
                ))
                for metric, field in (
                    ("inact", "inact_mu"),
                    ("block_xors", "block_xors_mu"),
                    ("block_muladds", "block_muladds_mu"),
                )
                for outcome, selected_rows in (
                    ("all_cells", rows),
                    ("success", success_rows),
                    ("failure_path", failure_rows),
                )
            },
        }
    comparisons: dict[str, Any] = {}
    for arm in ARMS[1:]:
        repairs = introductions = both_fail = 0
        new_shortfall = cleared_shortfall = 0
        base_inact = cand_inact = Decimal(0)
        base_xors = cand_xors = Decimal(0)
        base_muladds = cand_muladds = Decimal(0)
        paired_success = 0
        for seed_index, schedule, K in keys:
            base = indexed[("d12", seed_index, schedule, K)]
            cand = indexed[(arm, seed_index, schedule, K)]
            bf, cf = failed(base), failed(cand)
            repairs += bf and not cf
            introductions += not bf and cf
            both_fail += bf and cf
            bs = int(base["heavy_shortfall"]) != 0
            cs = int(cand["heavy_shortfall"]) != 0
            new_shortfall += not bs and cs
            cleared_shortfall += bs and not cs
            if not bf and not cf:
                paired_success += 1
                base_inact += Decimal(base["inact_mu"])
                cand_inact += Decimal(cand["inact_mu"])
                base_xors += Decimal(base["block_xors_mu"])
                cand_xors += Decimal(cand["block_xors_mu"])
                base_muladds += Decimal(base["block_muladds_mu"])
                cand_muladds += Decimal(cand["block_muladds_mu"])
        if paired_success + repairs + introductions + both_fail != len(keys):
            raise ValueError(f"{arm} paired outcome partition is incomplete")
        comparisons[f"{arm}_vs_d12"] = {
            "repairs": repairs, "introductions": introductions,
            "both_fail": both_fail, "new_heavy_shortfall": new_shortfall,
            "cleared_heavy_shortfall": cleared_shortfall,
            "one_sided_p": binomial_one_sided(repairs, introductions),
            "paired_success_cells": paired_success,
            "paired_success_base_inact_sum": str(base_inact),
            "paired_success_candidate_inact_sum": str(cand_inact),
            "inact_ratio_paired_success": (
                str(cand_inact / base_inact) if base_inact else None
            ),
            "paired_success_base_block_xors_sum": str(base_xors),
            "paired_success_candidate_block_xors_sum": str(cand_xors),
            "block_xor_ratio_paired_success": (
                str(cand_xors / base_xors) if base_xors else None
            ),
            "paired_success_base_block_muladds_sum": str(base_muladds),
            "paired_success_candidate_block_muladds_sum": str(cand_muladds),
            "block_muladd_ratio_paired_success": (
                str(cand_muladds / base_muladds) if base_muladds else None
            ),
            "failure_path_costs_excluded_from_ratios": True,
        }
    consistency: dict[str, Any] = {}
    base_counts = {
        K: sum(failed(indexed[("d12", seed_index, schedule, K)])
               for seed_index, _, schedule in expected_strata())
        for K in Ks
    }
    for arm in ARMS[1:]:
        candidate_counts = {
            K: sum(failed(indexed[(arm, seed_index, schedule, K)])
                   for seed_index, _, schedule in expected_strata())
            for K in Ks
        }
        consistency[f"{arm}_vs_d12"] = {
            "K_better": sum(candidate_counts[K] < base_counts[K] for K in Ks),
            "K_worse": sum(candidate_counts[K] > base_counts[K] for K in Ks),
            "K_equal": sum(candidate_counts[K] == base_counts[K] for K in Ks),
            "max_fail_strata": max(candidate_counts.values(), default=0),
            "worse_K": sorted(K for K in Ks if candidate_counts[K] > base_counts[K]),
        }
    return {
        "K_count": len(Ks), "cells_per_arm": len(Ks) * 9,
        "arms": arms, "comparisons": comparisons, "consistency": consistency,
    }


def analyze(
    scopes: dict[str, set[int]],
    source: dict[str, Any],
    records: list[dict[str, str]],
) -> tuple[dict[str, Any], list[dict[str, str]]]:
    indexed: dict[tuple[str, int, str, int], dict[str, str]] = {}
    for row in records:
        key = (row["arm"], int(row["seed_index"]), row["schedule"], int(row["N"]))
        if key in indexed:
            raise ValueError(f"duplicate result {key}")
        indexed[key] = row
    selected = set().union(*scopes.values())
    expected = len(ARMS) * 9 * len(selected)
    if len(indexed) != expected:
        raise ValueError(f"result count {len(indexed)}, want {expected}")

    low_identity_mismatches = 0
    for seed_index, _, schedule in expected_strata():
        for K in selected:
            if K >= CUTOFF:
                continue
            base = indexed[("d12", seed_index, schedule, K)]
            for arm in ("two_anchor_adaptive", "d13_adaptive"):
                if not deterministic_equal(base, indexed[(arm, seed_index, schedule, K)]):
                    low_identity_mismatches += 1
    if low_identity_mismatches:
        raise ValueError(f"adaptive low-K identity mismatches={low_identity_mismatches}")

    archive_mismatches = 0
    for source_row in source["hard_cases"]:
        key = (
            "d12", int(source_row["seed_index"]), source_row["schedule"],
            int(source_row["K"]),
        )
        row = indexed[key]
        checks = {
            "rank_fail": source_row["base_rank_fail"],
            "error": source_row["base_error"],
            "heavy_shortfall": source_row["base_heavy_shortfall"],
            "inact_mu": source_row["base_inact_mu"],
            "inact_max": source_row["base_inact_max"],
            "seed_attempt": source_row["base_seed_attempt"],
            "binary_def_mu": source_row["base_binary_def_mu"],
            "binary_def_max": source_row["base_binary_def_max"],
            "binary_def_hist": source_row["base_binary_def_hist"],
            "block_xors_mu": source_row["base_block_xors_mu"],
            "block_muladds_mu": source_row["base_block_muladds_mu"],
        }
        archive_mismatches += any(row[field] != value for field, value in checks.items())
    if archive_mismatches:
        raise ValueError(f"D12 archive reproduction mismatches={archive_mismatches}")

    summary = {
        name: summarize_scope(Ks, indexed) for name, Ks in scopes.items()
    }
    summary["all_selected_unique"] = summarize_scope(selected, indexed)
    summary["provenance"] = {
        "revealed_hard_source_cases": len(source["hard_cases"]),
        "revealed_hard_unique_K": len(source["hard_ks"]),
        "hard_source_case_x_stratum_mappings": len(source["hard_cases"]) * 9,
        "hard_unique_execution_cells_per_arm": len(source["hard_ks"]) * 9,
        "adaptive_low_K_identity_mismatches": low_identity_mismatches,
        "d12_archive_reproduction_mismatches": archive_mismatches,
    }

    per_k: list[dict[str, str]] = []
    for K in sorted(selected):
        row = {
            "K": str(K), "salt": source["states"][K]["salt"],
            "revealed_hard_k": str(int(K in scopes["revealed_hard_k"])),
            "shared_success_control_k": str(int(K in scopes["shared_success_control_k"])),
            "cutoff_neighbor": str(int(K in scopes["cutoff_neighbors"])),
        }
        for arm in ARMS:
            arm_rows = [
                indexed[(arm, seed_index, schedule, K)]
                for seed_index, _, schedule in expected_strata()
            ]
            row[f"{arm}_fail_strata"] = str(sum(failed(value) for value in arm_rows))
            row[f"{arm}_heavy_shortfall_strata"] = str(
                sum(int(value["heavy_shortfall"]) for value in arm_rows)
            )
        per_k.append(row)
    return summary, per_k


def write_csv(path: Path, rows: Iterable[dict[str, Any]], fields: Sequence[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def thermal_scalar(text: str, context: str) -> float:
    if THERMAL_SCALAR_PATTERN.fullmatch(text) is None:
        raise ValueError(f"{context} is not a canonical nonnegative scalar")
    value = float(text)
    if not math.isfinite(value):
        raise ValueError(f"{context} is not finite")
    return value


def thermal_counter(text: str, context: str) -> int:
    if THERMAL_COUNTER_PATTERN.fullmatch(text) is None:
        raise ValueError(f"{context} is not a canonical unsigned counter")
    return int(text, 10)


def parse_thermal_sample(
    line: bytes,
    context: str,
    allow_dimm_read_errors: bool = False,
) -> dict[str, Any]:
    if len(line) > THERMAL_ROW_MAX_BYTES:
        raise ValueError(f"{context} exceeds the thermal row size limit")
    if not line.endswith(b"\n"):
        raise ValueError(f"{context} is not a complete CSV row")
    try:
        text = line.decode("utf-8")[:-1]
    except UnicodeDecodeError as exc:
        raise ValueError(f"{context} is not UTF-8") from exc
    if text.endswith("\r"):
        text = text[:-1]
    values = text.split(",")
    if len(values) != len(THERMAL_FIELDS):
        raise ValueError(f"{context} has the wrong field count")
    row = dict(zip(THERMAL_FIELDS, values))
    if not row["utc"]:
        raise ValueError(f"{context} has an empty UTC timestamp")
    monotonic_s = thermal_scalar(row["monotonic_s"], f"{context}:monotonic_s")
    cpu_busy_pct = thermal_scalar(row["cpu_busy_pct"], f"{context}:cpu_busy_pct")
    cpu_avg_mhz = thermal_scalar(row["cpu_avg_mhz"], f"{context}:cpu_avg_mhz")
    cpu_tctl_c = thermal_scalar(row["cpu_tctl_c"], f"{context}:cpu_tctl_c")
    dimm_read_errors = thermal_counter(
        row["dimm_read_errors"], f"{context}:dimm_read_errors"
    )
    missing_dimm_fields = tuple(
        field for field in THERMAL_DIMM_FIELDS if row[field] == ""
    )
    if allow_dimm_read_errors:
        if len(missing_dimm_fields) != dimm_read_errors:
            raise ValueError(
                f"{context} DIMM read-error count does not match missing "
                "temperature fields"
            )
        if len(missing_dimm_fields) == len(THERMAL_DIMM_FIELDS):
            raise ValueError(f"{context} has no readable DIMM temperatures")
    elif missing_dimm_fields:
        raise ValueError(f"{context} has a missing DIMM temperature")
    dimm_temperatures = tuple(
        thermal_scalar(row[field], f"{context}:{field}")
        for field in THERMAL_DIMM_FIELDS if row[field] != ""
    )
    loads = tuple(
        thermal_scalar(row[field], f"{context}:{field}")
        for field in ("load1", "load5", "load15")
    )
    edac_ce = thermal_counter(row["edac_ce"], f"{context}:edac_ce")
    edac_ue = thermal_counter(row["edac_ue"], f"{context}:edac_ue")
    if not 0.0 <= cpu_busy_pct <= 100.0:
        raise ValueError(f"{context} CPU busy percentage is out of range")
    if cpu_avg_mhz <= 0.0:
        raise ValueError(f"{context} CPU frequency is not positive")
    if not 0.0 <= cpu_tctl_c <= 120.0:
        raise ValueError(f"{context} CPU temperature is out of range")
    if any(not 0.0 <= value <= 100.0 for value in dimm_temperatures):
        raise ValueError(f"{context} DIMM temperature is out of range")
    return {
        "raw": line,
        "monotonic_s": monotonic_s,
        "cpu_busy_pct": cpu_busy_pct,
        "cpu_avg_mhz": cpu_avg_mhz,
        "cpu_tctl_c": cpu_tctl_c,
        "dimm_temperatures": dimm_temperatures,
        "dimm_read_errors": dimm_read_errors,
        "loads": loads,
        "edac_ce": edac_ce,
        "edac_ue": edac_ue,
        "max_temperature_c": max((cpu_tctl_c, *dimm_temperatures)),
    }


def validate_thermal_current(
    sample: dict[str, Any], stale_seconds: float, context: str,
) -> None:
    if not math.isfinite(stale_seconds) or stale_seconds <= 0.0:
        raise ValueError("thermal stale interval must be positive and finite")
    try:
        with Path("/proc/uptime").open("r", encoding="ascii") as source:
            uptime_fields = source.read().split()
    except OSError as exc:
        raise ValueError("cannot read canonical system uptime") from exc
    if len(uptime_fields) != 2:
        raise ValueError("canonical system uptime has the wrong schema")
    now = thermal_scalar(uptime_fields[0], "system uptime")
    if (now - sample["monotonic_s"] > stale_seconds or
            sample["monotonic_s"] - now > 1.0):
        raise ValueError(f"{context} is stale or future-dated")


def validate_thermal_health(
    sample: dict[str, Any], context: str, require_zero_edac: bool = True,
    require_zero_dimm_errors: bool = True,
) -> None:
    if require_zero_dimm_errors and sample["dimm_read_errors"] != 0:
        raise ValueError(f"{context} reports a DIMM read error")
    if (require_zero_edac and
            (sample["edac_ce"] != 0 or sample["edac_ue"] != 0)):
        raise ValueError(f"{context} reports a nonzero EDAC counter")


def validate_thermal_log_descriptor(
    source: BinaryIO,
    path: Path,
) -> os.stat_result:
    descriptor = os.fstat(source.fileno())
    try:
        named = os.stat(path, follow_symlinks=False)
    except OSError as exc:
        raise ValueError(f"thermal log pathname is unavailable: {path}") from exc
    identity = lambda value: (
        value.st_dev, value.st_ino, value.st_mode, value.st_nlink,
    )
    if (not stat.S_ISREG(descriptor.st_mode) or descriptor.st_nlink != 1 or
            not stat.S_ISREG(named.st_mode) or named.st_nlink != 1 or
            identity(descriptor) != identity(named)):
        raise ValueError(f"thermal log is not a unique regular file: {path}")
    return descriptor


def open_thermal_log(path: Path) -> BinaryIO:
    flags = (os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
             getattr(os, "O_NOFOLLOW", 0) |
             getattr(os, "O_NONBLOCK", 0))
    try:
        descriptor = os.open(path, flags)
    except OSError as exc:
        if exc.errno == errno.ELOOP:
            raise ValueError(f"thermal log is a symlink: {path}") from exc
        raise ValueError(f"thermal log is unavailable: {path}") from exc
    source: BinaryIO | None = None
    try:
        # Critical header/baseline rereads must observe same-inode writes;
        # BufferedReader may otherwise serve stale bytes after seek().
        source = os.fdopen(descriptor, "rb", buffering=0)
        validate_thermal_log_descriptor(source, path)
    except BaseException:
        try:
            if source is None:
                os.close(descriptor)
            else:
                source.close()
        except OSError:
            pass
        raise
    if source is None:
        raise AssertionError("thermal log descriptor wrapper is unavailable")
    return source


def thermal_start(
    path: Path,
    stale_seconds: float = 5.0,
    require_zero_edac: bool = True,
    require_zero_dimm_errors: bool = True,
) -> dict[str, Any]:
    deadline = time.monotonic() + 5.0
    while True:
        with open_thermal_log(path) as source:
            header = source.readline()
            if not header.endswith(b"\n"):
                raise ValueError("thermal log has an invalid header")
            try:
                header_text = header.decode("utf-8")[:-1]
            except UnicodeDecodeError as exc:
                raise ValueError("thermal log has an invalid header") from exc
            if header_text.endswith("\r"):
                header_text = header_text[:-1]
            if header_text != ",".join(THERMAL_FIELDS):
                raise ValueError("thermal log does not use the canonical schema")
            source.seek(0, os.SEEK_END)
            offset = source.tell()
            if offset <= len(header):
                raise ValueError("thermal log has no data rows")
            source.seek(offset - 1)
            complete = source.read(1) == b"\n"
            stat_before = validate_thermal_log_descriptor(source, path)
            if complete:
                source.seek(max(len(header), offset - THERMAL_ROW_MAX_BYTES))
                tail = source.read(offset - source.tell())
                lines = tail.splitlines(keepends=True)
                if not lines:
                    raise ValueError("thermal log has no complete data row")
                baseline_row = lines[-1]
                row_start = offset - len(baseline_row)
                if row_start < len(header):
                    raise ValueError("thermal baseline overlaps its header")
                if row_start > len(header):
                    source.seek(row_start - 1)
                    if source.read(1) != b"\n":
                        raise ValueError(
                            "thermal baseline exceeds the row size limit")
                stat_before = validate_thermal_log_descriptor(source, path)
        if complete:
            break
        if time.monotonic() >= deadline:
            raise ValueError("thermal log ends in a persistent partial row")
        time.sleep(0.05)
    baseline = parse_thermal_sample(
        baseline_row, "thermal baseline",
        allow_dimm_read_errors=not require_zero_dimm_errors)
    wall_age = time.time() - stat_before.st_mtime
    if wall_age > stale_seconds or wall_age < -1.0:
        raise ValueError("thermal baseline file timestamp is stale or future-dated")
    validate_thermal_current(baseline, stale_seconds, "thermal baseline")
    validate_thermal_health(
        baseline, "thermal baseline", require_zero_edac=require_zero_edac,
        require_zero_dimm_errors=require_zero_dimm_errors)
    return {
        "dev": stat_before.st_dev, "ino": stat_before.st_ino,
        "offset": offset, "header": header,
        "baseline_row": baseline_row, "baseline": baseline,
        "edac_ce": baseline["edac_ce"], "edac_ue": baseline["edac_ue"],
        "monotonic_s": baseline["monotonic_s"],
        "max_temperature_c": baseline["max_temperature_c"],
    }


def thermal_finish(
    path: Path,
    start: dict[str, Any],
    output: Path,
    limit_c: float = 90.0,
    stale_seconds: float = 5.0,
    require_zero_edac: bool = True,
    require_zero_dimm_errors: bool = True,
    dimm_limit_c: float | None = None,
) -> dict[str, Any]:
    if dimm_limit_c is None:
        dimm_limit_c = limit_c
    if not math.isfinite(limit_c) or not 0.0 <= limit_c <= 120.0:
        raise ValueError("thermal limit is outside the canonical range")
    if (not math.isfinite(dimm_limit_c) or
            not 0.0 <= dimm_limit_c <= 120.0):
        raise ValueError("DIMM thermal limit is outside the canonical range")
    suffix = b""
    deadline = time.monotonic() + 5.0
    while True:
        with open_thermal_log(path) as source:
            stat_before_read = validate_thermal_log_descriptor(source, path)
            if ((stat_before_read.st_dev, stat_before_read.st_ino) !=
                    (start["dev"], start["ino"])):
                raise ValueError("thermal log inode changed")
            if stat_before_read.st_size < start["offset"]:
                raise ValueError("thermal log shrank")
            header = start.get("header")
            baseline_row = start.get("baseline_row")
            if (not isinstance(header, bytes) or
                    not isinstance(baseline_row, bytes) or
                    start["offset"] < len(header) + len(baseline_row)):
                raise ValueError("thermal start record is malformed")
            baseline_start = start["offset"] - len(baseline_row)
            source.seek(0)
            if source.read(len(header)) != header:
                raise ValueError("thermal log header changed after the mark")
            source.seek(baseline_start)
            if source.read(len(baseline_row)) != baseline_row:
                raise ValueError("thermal baseline row changed after the mark")
            source.seek(start["offset"])
            suffix = source.read()
            source.seek(0)
            if source.read(len(header)) != header:
                raise ValueError("thermal log header changed during finish")
            source.seek(baseline_start)
            if source.read(len(baseline_row)) != baseline_row:
                raise ValueError("thermal baseline row changed during finish")
            stat_after = validate_thermal_log_descriptor(source, path)
            if stat_after.st_size < start["offset"]:
                raise ValueError("thermal log shrank")
            # Reread after the post-read stat as well.  Metadata timestamps
            # can have coarser granularity than a same-length in-place write.
            source.seek(0)
            if source.read(len(header)) != header:
                raise ValueError("thermal log header changed at finish")
            source.seek(baseline_start)
            if source.read(len(baseline_row)) != baseline_row:
                raise ValueError("thermal baseline row changed at finish")
            stat_final = validate_thermal_log_descriptor(source, path)
            stable_identity = lambda value: (
                value.st_dev, value.st_ino, value.st_mode, value.st_nlink,
                value.st_size, value.st_mtime_ns, value.st_ctime_ns,
            )
            captured_to_stable_end = (
                stable_identity(stat_before_read) == stable_identity(stat_after) ==
                stable_identity(stat_final) and
                start["offset"] + len(suffix) == stat_final.st_size)
            stat_after = stat_final
        if captured_to_stable_end and suffix.endswith(b"\n"):
            break
        if time.monotonic() >= deadline:
            raise ValueError("thermal log ends in a persistent partial interval row")
        time.sleep(0.05)
    wall_age = time.time() - stat_after.st_mtime
    if wall_age > stale_seconds or wall_age < -1.0:
        raise ValueError("thermal interval file timestamp is stale or future-dated")
    interval_lines = suffix.splitlines(keepends=True)
    if not interval_lines:
        raise ValueError("thermal interval contains no samples")
    samples = [start["baseline"]]
    for number, line in enumerate(interval_lines, 1):
        samples.append(parse_thermal_sample(
            line, f"thermal interval row {number}",
            allow_dimm_read_errors=not require_zero_dimm_errors))
    for number, (previous, sample) in enumerate(zip(samples, samples[1:]), 1):
        delta = sample["monotonic_s"] - previous["monotonic_s"]
        if delta <= 0.0 or delta > stale_seconds:
            raise ValueError(
                f"thermal interval row {number} is nonmonotonic or follows a gap"
            )
    validate_thermal_current(samples[-1], stale_seconds, "thermal interval tail")
    for number, sample in enumerate(samples):
        validate_thermal_health(
            sample, f"thermal sealed row {number}",
            require_zero_edac=require_zero_edac,
            require_zero_dimm_errors=require_zero_dimm_errors)
    output.write_bytes(start["header"] + start["baseline_row"] + suffix)
    interval_samples = samples[1:]
    cpu_busy_values = [sample["cpu_busy_pct"] for sample in interval_samples]
    cpu_temperature_values = [sample["cpu_tctl_c"] for sample in samples]
    dimm_temperature_values = [
        value for sample in samples for value in sample["dimm_temperatures"]
    ]
    consecutive_high = 0
    max_consecutive_high = consecutive_high
    high_samples = 0
    for sample in samples:
        if (sample["cpu_tctl_c"] >= limit_c or
                max(sample["dimm_temperatures"]) >= dimm_limit_c):
            consecutive_high += 1
            high_samples += 1
            max_consecutive_high = max(max_consecutive_high, consecutive_high)
        else:
            consecutive_high = 0
    edac_ce_values = [sample["edac_ce"] for sample in samples]
    edac_ue_values = [sample["edac_ue"] for sample in samples]
    def counter_delta(values: Sequence[int]) -> int:
        """Return the net delta, or a negative rollback step as a gap marker."""
        steps = [current - previous
                 for previous, current in zip(values, values[1:])]
        rollback = min(steps, default=0)
        if rollback < 0:
            # The raw sealed rows remain the authority.  A negative result is
            # deliberately propagated so permissive timing callers classify
            # the counter discontinuity as a telemetry gap instead of hiding
            # it behind a nonnegative net change.
            return rollback
        return values[-1] - values[0]

    edac_ce_delta = counter_delta(edac_ce_values)
    edac_ue_delta = counter_delta(edac_ue_values)
    return {
        "samples": len(interval_samples),
        "sealed_samples_including_baseline": len(samples),
        "cpu_busy_min_pct": min(cpu_busy_values),
        "cpu_tctl_max_c": max(cpu_temperature_values),
        "dimm_max_c": max(dimm_temperature_values),
        "dimm_read_errors_max": max(
            sample["dimm_read_errors"] for sample in samples
        ),
        "edac_ce_delta": edac_ce_delta,
        "edac_ue_delta": edac_ue_delta,
        "thermal_limit_c": limit_c,
        "thermal_high_samples": high_samples,
        "thermal_high_max_consecutive_samples": max_consecutive_high,
    }


def main() -> int:
    """Reject the retired standalone runner before touching any input."""
    print(STANDALONE_RETIREMENT_MESSAGE, file=sys.stderr)
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
