#!/usr/bin/env python3
"""Paired production-H12 two-anchor screen over revealed development inputs.

This one-off development screen consumes only the published gen2b all-K
analysis, but evaluates the exact mixed10-GF(256)/2-GF(65536), period-244,
mix2 production geometry.  Its selected inputs are not independent evidence;
the companion all-K campaign provides the independent confirmation.  This
script does not locate, inspect, or depend on any H15-v5 holdout material.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import queue
import re
import subprocess
import sys
import time
from collections import defaultdict
from decimal import Decimal, getcontext
from typing import Any, Iterable, Sequence


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


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


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
    if sha256_file(path) != SOURCE_CELLS_SHA256:
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
    with path.open("r", encoding="utf-8", newline="") as source:
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
    cpu = cpu_pool.acquire()
    pinned = ["taskset", "-c", str(cpu), *command]
    start_ns = time.time_ns()
    try:
        process = subprocess.run(
            pinned, check=False, text=True, encoding="utf-8",
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            timeout=timeout,
        )
    finally:
        cpu_pool.release(cpu)
    return (
        pinned, process.returncode, process.stdout, process.stderr,
        start_ns, time.time_ns(),
    )


def run_identity(
    baseline: Path,
    candidate: Path,
    result_dir: Path,
    cpu_pool: CpuPool,
    timeout: float,
) -> dict[str, Any]:
    identity_dir = result_dir / "identity"
    identity_dir.mkdir()
    jobs: list[tuple[str, Path, int, str, str]] = []
    for role, binary in (("baseline", baseline), ("candidate", candidate)):
        for seed_index, seed, schedule in expected_strata():
            jobs.append((role, binary, seed_index, seed, schedule))

    def one(job: tuple[str, Path, int, str, str]) -> dict[str, Any]:
        role, binary, seed_index, seed, schedule = job
        command = make_command(binary, IDENTITY_KS, seed, schedule, "0x0")
        pinned, returncode, stdout, stderr, start_ns, end_ns = execute(
            command, cpu_pool, timeout
        )
        stem = f"{role}.seed{seed_index}.{schedule}"
        out_path = identity_dir / f"{stem}.stdout"
        err_path = identity_dir / f"{stem}.stderr"
        out_path.write_text(stdout, encoding="utf-8", newline="\n")
        err_path.write_text(stderr, encoding="utf-8", newline="\n")
        if returncode != 0 or stderr:
            raise RuntimeError(f"identity {stem} failed rc={returncode}: {stderr[-1000:]}")
        rows = parse_bench_output(
            stdout, IDENTITY_KS, "0x0", seed, schedule, BASE_OPTIONS,
            role == "candidate"
        )
        return {
            "role": role, "seed_index": seed_index, "schedule": schedule,
            "command": pinned, "start_ns": start_ns, "end_ns": end_ns,
            "stdout": str(out_path.relative_to(result_dir)),
            "stdout_sha256": sha256_file(out_path),
            "stderr": str(err_path.relative_to(result_dir)),
            "stderr_sha256": sha256_file(err_path), "rows": rows,
        }

    results: list[dict[str, Any]] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=len(jobs)) as executor:
        futures = [executor.submit(one, job) for job in jobs]
        for future in concurrent.futures.as_completed(futures):
            results.append(future.result())
    indexed: dict[tuple[str, int, str, int], dict[str, str]] = {}
    for result in results:
        for row in result.pop("rows"):
            indexed[(
                result["role"], result["seed_index"], result["schedule"],
                int(row["N"]),
            )] = row
    comparisons = 0
    for seed_index, _, schedule in expected_strata():
        for K in IDENTITY_KS:
            base = indexed[("baseline", seed_index, schedule, K)]
            cand = indexed[("candidate", seed_index, schedule, K)]
            if base.keys() != cand.keys():
                raise ValueError("identity CSV schemas differ")
            for field in base:
                if field in TIMING_FIELDS:
                    continue
                comparisons += 1
                if base[field] != cand[field]:
                    raise ValueError(
                        f"disabled-hook identity mismatch seed={seed_index} "
                        f"schedule={schedule} K={K} field={field}: "
                        f"{base[field]} != {cand[field]}"
                    )
    results.sort(key=lambda value: (
        value["role"], value["seed_index"], value["schedule"]
    ))
    return {
        "passed": True,
        "cells_per_binary": len(expected_strata()) * len(IDENTITY_KS),
        "non_timing_field_comparisons": comparisons,
        "baseline_binary_sha256": sha256_file(baseline),
        "candidate_binary_sha256": sha256_file(candidate),
        "jobs": results,
    }


def run_task(
    task: dict[str, Any],
    binary: Path,
    result_dir: Path,
    cpu_pool: CpuPool,
    timeout: float,
) -> dict[str, Any]:
    options = arm_options(task["arm"], task["band"])
    command = make_command(
        binary, task["Ks"], task["seed"], task["schedule"],
        task["salt"], options,
    )
    pinned, returncode, stdout, stderr, start_ns, end_ns = execute(
        command, cpu_pool, timeout
    )
    stem = f"job{task['job']:04d}.{task['arm']}.seed{task['seed_index']}.{task['schedule']}.{task['band']}"
    out_path = result_dir / "stdout" / f"{stem}.csv"
    err_path = result_dir / "stderr" / f"{stem}.txt"
    out_path.write_text(stdout, encoding="utf-8", newline="\n")
    err_path.write_text(stderr, encoding="utf-8", newline="\n")
    if returncode != 0 or stderr:
        raise RuntimeError(f"{stem} failed rc={returncode}: {stderr[-1000:]}")
    rows = parse_bench_output(
        stdout, task["Ks"], task["salt"], task["seed"],
        task["schedule"], options, True
    )
    for row in rows:
        row.update({
            "arm": task["arm"], "band": task["band"],
            "seed_index": str(task["seed_index"]), "seed": task["seed"],
            "schedule": task["schedule"], "salt": task["salt"],
        })
    return {
        "job": task["job"], "arm": task["arm"], "band": task["band"],
        "seed_index": task["seed_index"], "seed": task["seed"],
        "schedule": task["schedule"], "salt": task["salt"],
        "Ks": task["Ks"], "command": pinned,
        "start_ns": start_ns, "end_ns": end_ns, "returncode": returncode,
        "stdout": str(out_path.relative_to(result_dir)),
        "stdout_sha256": sha256_file(out_path),
        "stderr": str(err_path.relative_to(result_dir)),
        "stderr_sha256": sha256_file(err_path), "records": rows,
    }


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
        arms[arm] = {
            "cells": len(rows),
            "failures": sum(failed(row) for row in rows),
            "errors": sum(int(row["error"]) for row in rows),
            "heavy_shortfalls": sum(int(row["heavy_shortfall"]) for row in rows),
            "block_xors_all": str(sum(Decimal(row["block_xors_mu"]) for row in rows)),
            "block_muladds_all": str(sum(Decimal(row["block_muladds_mu"]) for row in rows)),
        }
    comparisons: dict[str, Any] = {}
    for arm in ARMS[1:]:
        repairs = introductions = both_fail = 0
        new_shortfall = cleared_shortfall = 0
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
                base_xors += Decimal(base["block_xors_mu"])
                cand_xors += Decimal(cand["block_xors_mu"])
                base_muladds += Decimal(base["block_muladds_mu"])
                cand_muladds += Decimal(cand["block_muladds_mu"])
        comparisons[f"{arm}_vs_d12"] = {
            "repairs": repairs, "introductions": introductions,
            "both_fail": both_fail, "new_heavy_shortfall": new_shortfall,
            "cleared_heavy_shortfall": cleared_shortfall,
            "one_sided_p": binomial_one_sided(repairs, introductions),
            "paired_success_cells": paired_success,
            "block_xor_ratio_paired_success": str(cand_xors / base_xors),
            "block_muladd_ratio_paired_success": str(cand_muladds / base_muladds),
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


def thermal_start(
    path: Path,
    stale_seconds: float = 5.0,
    require_zero_edac: bool = True,
    require_zero_dimm_errors: bool = True,
) -> dict[str, Any]:
    deadline = time.monotonic() + 5.0
    while True:
        with path.open("rb") as source:
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
            stat_before = os.fstat(source.fileno())
            if complete:
                source.seek(max(len(header), offset - 16384))
                tail = source.read(offset - source.tell())
        if complete:
            break
        if time.monotonic() >= deadline:
            raise ValueError("thermal log ends in a persistent partial row")
        time.sleep(0.05)
    lines = tail.splitlines(keepends=True)
    if not lines:
        raise ValueError("thermal log has no complete data row")
    baseline_row = lines[-1]
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
        with path.open("rb") as source:
            stat_after = os.fstat(source.fileno())
            if ((stat_after.st_dev, stat_after.st_ino) !=
                    (start["dev"], start["ino"])):
                raise ValueError("thermal log inode changed")
            if stat_after.st_size < start["offset"]:
                raise ValueError("thermal log shrank")
            source.seek(start["offset"])
            suffix = source.read()
        if suffix.endswith(b"\n"):
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
    parser = argparse.ArgumentParser()
    parser.add_argument("--cells", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--baseline-binary", type=Path, required=True)
    parser.add_argument("--thermal", type=Path, required=True)
    parser.add_argument("--result-dir", type=Path, required=True)
    parser.add_argument("--controls", type=int, default=512)
    parser.add_argument("--workers", type=int, default=120)
    parser.add_argument("--chunk-size", type=int, default=16)
    parser.add_argument("--timeout", type=float, default=600.0)
    args = parser.parse_args()
    if args.controls < 512 or args.workers <= 0 or args.chunk_size <= 0 or args.timeout <= 0:
        raise ValueError("controls>=512 and positive workers/chunk-size/timeout required")

    cells = args.cells.resolve(strict=True)
    binary = args.binary.resolve(strict=True)
    baseline = args.baseline_binary.resolve(strict=True)
    thermal = args.thermal.resolve(strict=True)
    for name, path in (("binary", binary), ("baseline binary", baseline)):
        if not path.is_file() or not os.access(path, os.X_OK):
            raise ValueError(f"{name} is not an executable regular file")
    if sha256_file(baseline) != BASELINE_BINARY_SHA256:
        raise ValueError("baseline binary is not the exact preserved 09d44f1 build")

    result_dir = args.result_dir.resolve()
    result_dir.mkdir(parents=True, exist_ok=False)
    (result_dir / "stdout").mkdir()
    (result_dir / "stderr").mkdir()
    thermal_mark = thermal_start(thermal)

    source = load_source_cells(cells)
    scopes = select_scopes(source, args.controls)
    selected = set().union(*scopes.values())
    tasks = build_tasks(selected, source["states"], args.chunk_size)
    cpus = sorted(os.sched_getaffinity(0))
    workers = min(args.workers, len(cpus), len(tasks))
    cpu_pool = CpuPool(cpus, workers)
    identity = run_identity(
        baseline, binary, result_dir, cpu_pool, args.timeout
    )

    task_results: list[dict[str, Any]] = []
    completed = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [
            executor.submit(
                run_task, task, binary, result_dir, cpu_pool, args.timeout
            ) for task in tasks
        ]
        for future in concurrent.futures.as_completed(futures):
            task_results.append(future.result())
            completed += 1
            if completed % 100 == 0 or completed == len(tasks):
                print(f"screen progress {completed}/{len(tasks)}", flush=True)
    task_results.sort(key=lambda result: result["job"])

    records: list[dict[str, str]] = []
    command_records: list[dict[str, Any]] = []
    for task_result in task_results:
        records.extend(task_result.pop("records"))
        command_records.append(task_result)
    records.sort(key=lambda row: (
        row["arm"], int(row["seed_index"]), row["schedule"], int(row["N"])
    ))
    record_fields = [
        "arm", "band", "seed_index", "seed", "schedule", "salt",
        *[field for field in records[0] if field not in {
            "arm", "band", "seed_index", "seed", "schedule", "salt"
        }],
    ]
    write_csv(result_dir / "results.csv", records, record_fields)
    with (result_dir / "commands.jsonl").open(
        "w", encoding="utf-8", newline="\n"
    ) as output:
        for record in command_records:
            output.write(canonical_json(record) + "\n")
    (result_dir / "identity.json").write_text(
        json.dumps(identity, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8", newline="\n",
    )

    hard_fields = [
        "case", "K", "salt", "schedule", "seed_index", "seed",
        "base_failed", "candidate_failed", "base_rank_fail",
        "candidate_rank_fail", "base_error", "candidate_error",
    ]
    hard_rows = []
    for case, row in enumerate(source["hard_cases"]):
        value = dict(row)
        value["case"] = case
        hard_rows.append(value)
    write_csv(result_dir / "revealed_hard_cases.csv", hard_rows, hard_fields)
    selection_rows = []
    for K in sorted(selected):
        selection_rows.append({
            "K": K, "salt": source["states"][K]["salt"],
            "revealed_hard_k": int(K in scopes["revealed_hard_k"]),
            "shared_success_control_k": int(K in scopes["shared_success_control_k"]),
            "cutoff_neighbor": int(K in scopes["cutoff_neighbors"]),
        })
    write_csv(
        result_dir / "selected_k.csv", selection_rows,
        ("K", "salt", "revealed_hard_k", "shared_success_control_k", "cutoff_neighbor"),
    )

    getcontext().prec = 80
    summary, per_k = analyze(scopes, source, records)
    per_k_fields = list(per_k[0])
    write_csv(result_dir / "per_k.csv", per_k, per_k_fields)
    thermal_summary = thermal_finish(
        thermal, thermal_mark, result_dir / "thermal_interval.csv",
        dimm_limit_c=90.0,
    )
    summary["thermal"] = thermal_summary
    summary_path = result_dir / "summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8", newline="\n",
    )

    repo = Path(__file__).resolve().parents[1]
    head = subprocess.check_output(
        ("git", "rev-parse", "HEAD"), cwd=repo
    ).strip().decode("ascii")
    diff = subprocess.check_output(
        ("git", "diff", "--binary", "--no-ext-diff"), cwd=repo
    )
    status = subprocess.check_output(
        ("git", "status", "--short"), cwd=repo
    ).decode("utf-8")
    artifact_hashes = {
        str(path.relative_to(result_dir)): sha256_file(path)
        for path in sorted(result_dir.rglob("*")) if path.is_file()
    }
    manifest = {
        "schema": "wirehair.wh2.rank_floor_two_anchor_screen.v1",
        "source_cells": str(cells),
        "source_cells_sha256": sha256_file(cells),
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256_file(Path(__file__).resolve()),
        "binary": str(binary), "binary_sha256": sha256_file(binary),
        "baseline_binary": str(baseline),
        "baseline_binary_sha256": sha256_file(baseline),
        "git_head": head, "git_diff_sha256": sha256_bytes(diff),
        "git_status": status, "workers": workers,
        "chunk_size": args.chunk_size, "controls": args.controls,
        "hard_source_cases": len(source["hard_cases"]),
        "hard_unique_K": len(source["hard_ks"]),
        "selected_unique_K": len(selected), "jobs": len(tasks),
        "arms": ARMS, "cutoff": CUTOFF,
        "artifacts": artifact_hashes,
    }
    manifest_path = result_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8", newline="\n",
    )
    print(canonical_json({
        "result_dir": str(result_dir),
        "manifest_sha256": sha256_file(manifest_path),
        "summary": summary,
    }))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"rank-floor two-anchor screen: {exc}", file=sys.stderr)
        raise
