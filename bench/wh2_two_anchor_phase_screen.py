#!/usr/bin/env python3
"""Freeze and run the focused WH2 second-anchor q-phase screen.

The cohort is derived only from the sealed production-H12 all-K v3 result:
all 197 K values where D12 and q0 differed, all 26 K values where either arm
failed in at least two of the nine sealed strata, and seven explicit boundary
controls.  Four adaptive arms (D12, q0, q1, q2) run over that exact 226-K
cohort under the three revealed development seeds plus the three sealed fresh
seeds and three hard loss schedules.  Every packet-peel xor is zero.

Prepare must run from a clean immutable commit.  Run must execute the frozen
script, support modules, binary, cohort, and contract produced by prepare.
Saturated-run timings are provenance only and are never a speed claim.

This is deliberately post-selection tooling.  Its outcome-derived K cohort
and q-phase variants are forbidden inputs to raw architecture selection.  The
prepare command requires an explicit acknowledgement so an architecture
campaign cannot accidentally use this development screen.
"""

from __future__ import annotations

from collections import defaultdict
import argparse
import concurrent.futures
import csv
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
import heapq
import io
import json
import math
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import threading
import time
from typing import Any, Sequence

import wh2_rank_floor_two_anchor_allk as common
from wh2_rank_floor_two_anchor_screen import (
    THERMAL_FIELDS,
    TIMING_FIELDS,
    canonical_json,
    parse_bench_output,
    parse_preamble,
    thermal_start,
)


SOURCE_CELLS_SHA256 = (
    "ab4a1f1e31b59a302f9803c99bed29bb1bfbb1715a2a38d470196c5356d41587"
)
COHORT_KS_SHA256 = (
    "705098e645f3bdbe96fcc8d6589549d52aad2ffaf81cdc4997cddf7f84664db3"
)
CUTOFF = 4096
SOURCE_SEEDS = (
    "0xd1b54a32d192ed03",
    "0x94d049bb133111eb",
    "0x8538ecb5bd456ea3",
)
SEEDS = (
    "0xa6b527b9ae8de8d7",
    "0x75446e1619e81d8a",
    "0x007e9dd892a319f5",
    *SOURCE_SEEDS,
)
SCHEDULES = ("burst", "adversarial", "repair-only")
ARMS = ("d12", "q0", "q1", "q2")
EXPLICIT_KS = frozenset((4095, 4096, 46008, 46251, 46253, 46496, 64000))
RETAINED_CLEAR_KS = (6841, 44189, 53732)
TARGET_K = 46252
LARGE_PARTITIONS = 8
STATS = (
    "rank_fail", "error", "heavy_shortfall", "seed_attempt",
    "inact_milli", "binary_def_milli", "heavy_gain_milli",
    "block_xors_milli", "block_muladds_milli",
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


def die(message: str) -> None:
    raise common.CampaignError(message)


def sha256_bytes(data: bytes) -> str:
    return common.sha256_bytes(data)


def json_bytes(value: Any) -> bytes:
    return common.json_bytes(value)


def source_stat(row: dict[str, str], prefix: str, field: str) -> int:
    source_field = f"{prefix}_{field}"
    if source_field not in row:
        die(f"source cells missing {source_field}")
    text = row[source_field]
    try:
        value = int(text, 10)
    except ValueError:
        die(f"source cells {source_field} is not an integer: {text!r}")
    if value < 0:
        die(f"source cells {source_field} is negative")
    return value


def load_source_cohort(path: Path) -> dict[str, Any]:
    data = common.stable_bytes(path)
    if sha256_bytes(data) != SOURCE_CELLS_SHA256:
        die("sealed source paired_cells.csv SHA256 mismatch")
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as exc:
        die(f"sealed source cells are not UTF-8: {exc}")
    required = {
        "K", "active_packet_peel_seed_xor", "schedule", "seed_index",
        "seed",
        *(f"{prefix}_{field}" for prefix in ("d12", "two_anchor_adaptive")
          for field in STATS),
    }
    counts: dict[int, list[int]] = defaultdict(lambda: [0, 0, 0])
    seen: set[tuple[int, int, str]] = set()
    source_rows = 0
    reader = csv.DictReader(io.StringIO(text, newline=""))
    if reader.fieldnames is None or not required.issubset(reader.fieldnames):
        die("sealed source cells have an unexpected schema")
    for row in reader:
        source_rows += 1
        try:
            K = int(row["K"], 10)
            seed_index = int(row["seed_index"], 10)
        except ValueError:
            die("sealed source cells contain a malformed K or seed index")
        schedule = row["schedule"]
        key = (K, seed_index, schedule)
        if key in seen:
            die(f"duplicate sealed source cell {key}")
        seen.add(key)
        if not 2 <= K <= 64000:
            die(f"sealed source K outside 2..64000: {K}")
        if (not 0 <= seed_index < len(SOURCE_SEEDS) or
                row["seed"] != SOURCE_SEEDS[seed_index] or
                schedule not in SCHEDULES or
                int(row["active_packet_peel_seed_xor"], 0) != 0):
            die(f"sealed source stratum/packet-xor mismatch at {key}")
        d12_error = source_stat(row, "d12", "error")
        q0_error = source_stat(row, "two_anchor_adaptive", "error")
        d12_fail = source_stat(row, "d12", "rank_fail") + d12_error
        q0_fail = (
            source_stat(row, "two_anchor_adaptive", "rank_fail") + q0_error
        )
        if d12_error or q0_error or d12_fail not in (0, 1) or q0_fail not in (0, 1):
            die(f"sealed source has error or nonbinary outcome at {key}")
        counts[K][0] += d12_fail
        counts[K][1] += q0_fail
        counts[K][2] += d12_fail != q0_fail
    expected_rows = 63999 * len(SOURCE_SEEDS) * len(SCHEDULES)
    if source_rows != expected_rows or len(seen) != expected_rows:
        die(f"sealed source cardinality {source_rows}/{len(seen)}, want {expected_rows}")
    if set(counts) != set(range(2, 64001)):
        die("sealed source K domain is not exactly 2..64000")

    sensitive = {K for K, values in counts.items() if values[2] != 0}
    multi = {
        K for K, values in counts.items()
        if values[0] >= 2 or values[1] >= 2
    }
    cohort = sorted(sensitive | multi | EXPLICIT_KS)
    cohort_data = canonical_json(cohort).encode("utf-8")
    if (len(sensitive), len(multi), len(sensitive | multi), len(cohort)) != (
            197, 26, 219, 226):
        die("sealed source cohort component cardinalities changed")
    if sha256_bytes(cohort_data) != COHORT_KS_SHA256:
        die("derived 226-K cohort fingerprint mismatch")

    references: list[dict[str, Any]] = []
    cohort_set = set(cohort)
    reader = csv.DictReader(io.StringIO(text, newline=""))
    for row in reader:
        K = int(row["K"], 10)
        if K not in cohort_set:
            continue
        reference: dict[str, Any] = {
            "K": K,
            "source_seed_index": int(row["seed_index"], 10),
            "seed": row["seed"],
            "schedule": row["schedule"],
        }
        for output_arm, prefix in (
                ("d12", "d12"), ("q0", "two_anchor_adaptive")):
            for field in STATS:
                reference[f"{output_arm}_{field}"] = source_stat(
                    row, prefix, field
                )
        references.append(reference)
    references.sort(key=lambda row: (
        row["source_seed_index"], SCHEDULES.index(row["schedule"]), row["K"]
    ))
    if len(references) != len(cohort) * len(SOURCE_SEEDS) * len(SCHEDULES):
        die("sealed source reference extraction cardinality mismatch")
    return {
        "schema": "wirehair.wh2.two_anchor_phase.cohort.v1",
        "source_cells_sha256": SOURCE_CELLS_SHA256,
        "ks_sha256": COHORT_KS_SHA256,
        "ks": cohort,
        "selection": {
            "sensitive_k": len(sensitive),
            "multi_stratum_k": len(multi),
            "sensitive_or_multi_k": len(sensitive | multi),
            "explicit_k": sorted(EXPLICIT_KS),
            "cohort_k": len(cohort),
            "failure_resonance_threshold_strata": 2,
        },
        "source_references": references,
    }


def arm_options(arm: str, band: str) -> tuple[str, ...]:
    if arm not in ARMS or band not in ("small", "large"):
        die(f"unknown phase-screen arm/band {arm}/{band}")
    if arm == "d12" or band == "small":
        return ()
    if arm == "q0":
        phase = "0"
    elif arm == "q1":
        phase = "1"
    elif arm == "q2":
        phase = "2"
    return (
        "--binary-dense-two-anchor",
        "--binary-dense-two-anchor-phase", phase,
    )


@dataclass(frozen=True)
class Job:
    job: int
    arm: str
    band: str
    seed_index: int
    seed: str
    schedule: str
    chunk: int
    ks: tuple[int, ...]

    @property
    def stem(self) -> str:
        return (
            f"job{self.job:05d}.{self.arm}.seed{self.seed_index}."
            f"{self.schedule}.chunk{self.chunk:02d}.{self.band}"
        )


def build_jobs(ks: Sequence[int]) -> list[Job]:
    if (len(ks) != 226 or list(ks) != sorted(set(ks)) or
            sha256_bytes(canonical_json(list(ks)).encode("utf-8")) !=
            COHORT_KS_SHA256):
        die("frozen cohort is not the exact ordered 226-K set")
    small = tuple(K for K in ks if K < CUTOFF)
    large = tuple(K for K in ks if K >= CUTOFF)
    # Process-launch microbenchmarks on the complete 16,272-cell workload
    # found two coarse large-K partitions (216 jobs) 23.7% slower than the
    # 4-K reference, while eight K-sum-balanced partitions (648 jobs) were
    # within 0.2% and removed 84% of the process/manifest overhead.  K is a
    # stable proxy for construction work here.  The heap's partition id is a
    # deterministic tie-breaker, and each final command retains exact per-K
    # rows and attribution.
    partition_heap: list[tuple[int, int, list[int]]] = [
        (0, partition, []) for partition in range(LARGE_PARTITIONS)
    ]
    heapq.heapify(partition_heap)
    for K in reversed(large):
        weight, partition, values = heapq.heappop(partition_heap)
        values.append(K)
        heapq.heappush(partition_heap, (weight + K, partition, values))
    balanced = sorted(partition_heap, key=lambda item: item[1])
    chunks: list[tuple[str, int, tuple[int, ...]]] = [
        ("small", 0, small)
    ]
    chunks.extend(
        ("large", partition + 1, tuple(sorted(values)))
        for _weight, partition, values in balanced
    )
    if ([len(values) for _band, _chunk, values in chunks] !=
            [10, 27, 27, 27, 27, 27, 27, 27, 27]):
        die("focused balanced-partition cardinalities changed")
    jobs: list[Job] = []
    for arm in ARMS:
        for seed_index, seed in enumerate(SEEDS):
            for schedule in SCHEDULES:
                for band, chunk, values in chunks:
                    jobs.append(Job(
                        len(jobs), arm, band, seed_index, seed, schedule,
                        chunk, values,
                    ))
    expected_cells = len(ARMS) * len(SEEDS) * len(SCHEDULES) * len(ks)
    keys = {
        (job.arm, job.seed_index, job.schedule, K)
        for job in jobs for K in job.ks
    }
    if len(jobs) != 648 or len(keys) != expected_cells or expected_cells != 16272:
        die("focused job ledger cardinality mismatch")
    return jobs


def make_command(binary: Path, job: Job) -> list[str]:
    return [
        str(binary), *BASE_OPTIONS,
        "--N", ",".join(map(str, job.ks)),
        "--seed", job.seed,
        "--schedule", job.schedule,
        "--packet-peel-seed-xor", "0",
        *arm_options(job.arm, job.band),
    ]


def expected_phase(options: Sequence[str]) -> int:
    if "--binary-dense-two-anchor-phase" not in options:
        return 0
    index = options.index("--binary-dense-two-anchor-phase")
    if index + 1 >= len(options):
        die("phase option lacks a value")
    return int(options[index + 1], 10)


def run_job(
    job: Job,
    binary: Path,
    binary_sha256: str,
    taskset: Path,
    taskset_record: dict[str, str],
    result_root: Path,
    pool: common.CpuPool,
    abort: threading.Event,
    registry: common.ProcessRegistry,
    timeout: float,
) -> dict[str, Any]:
    if abort.is_set():
        die(f"campaign abort set before job {job.job}")
    cpu = pool.acquire()
    try:
        if abort.is_set():
            die(f"campaign abort set before job {job.job} launch")
        common.verify_frozen_binary(binary, binary_sha256)
        if common.frozen_executable_path(taskset_record, "taskset") != taskset:
            die("frozen taskset path changed before job launch")
        command = [str(taskset), "-c", str(cpu), *make_command(binary, job)]
        start_ns = time.time_ns()
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            start_new_session=True,
        )
        registered = False
        stdout_bytes = stderr_bytes = b""
        try:
            registry.add(process)
            registered = True
            deadline = time.monotonic() + timeout
            while True:
                remaining = deadline - time.monotonic()
                try:
                    stdout_bytes, stderr_bytes = process.communicate(
                        timeout=max(0.001, min(0.25, remaining))
                    )
                    break
                except subprocess.TimeoutExpired:
                    reason = None
                    if abort.is_set():
                        reason = "campaign abort"
                    elif remaining <= 0:
                        reason = f"timeout after {timeout:g}s"
                    if reason is None:
                        continue
                    common.stop_and_reap_process_group(process)
                    die(f"job {job.job} terminated on {reason}")
        finally:
            leader_live = process.poll() is None
            descendants_live = common.process_group_exists(process)
            if leader_live or descendants_live:
                common.stop_and_reap_process_group(process)
            else:
                common.close_process_streams(process)
            if process.poll() is None or common.process_group_exists(process):
                die(f"job {job.job} child cleanup was not proven")
            if registered:
                registry.remove(process)
            if not leader_live and descendants_live:
                die(f"job {job.job} left an unexpected child process")
        end_ns = time.time_ns()
    finally:
        pool.release(cpu)

    stdout_path = result_root / "stdout" / f"{job.stem}.csv"
    stderr_path = result_root / "stderr" / f"{job.stem}.txt"
    command_path = result_root / "commands" / f"{job.stem}.json"
    common.atomic_write(stdout_path, stdout_bytes)
    common.atomic_write(stderr_path, stderr_bytes)
    try:
        stdout = stdout_bytes.decode("utf-8")
        stderr = stderr_bytes.decode("utf-8")
    except UnicodeDecodeError as exc:
        die(f"job {job.job} emitted non-UTF-8 output: {exc}")
    if process.returncode != 0 or stderr_bytes:
        die(f"job {job.job} failed rc={process.returncode}: {stderr[-1000:]}")
    options = arm_options(job.arm, job.band)
    rows = parse_bench_output(
        stdout, job.ks, "0x0", job.seed, job.schedule, options, True,
    )
    preamble = parse_preamble(stdout.splitlines()[0])
    if int(preamble.get("binary_dense_two_anchor_phase", "-1"), 10) != \
            expected_phase(options):
        die(f"job {job.job} two-anchor phase preamble mismatch")
    record = {
        "job": job.job, "arm": job.arm, "band": job.band,
        "seed_index": job.seed_index, "seed": job.seed,
        "schedule": job.schedule, "chunk": job.chunk,
        "active_packet_peel_seed_xor": "0x0",
        "K_count": len(job.ks), "cpu": cpu, "command": command,
        "start_ns": start_ns, "end_ns": end_ns,
        "elapsed_ns": end_ns - start_ns, "returncode": process.returncode,
        "stdout": str(stdout_path.relative_to(result_root)),
        "stdout_sha256": common.sha256_file(stdout_path),
        "stderr": str(stderr_path.relative_to(result_root)),
        "stderr_sha256": common.sha256_file(stderr_path),
        "saturated_timing_speed_claim_valid": False,
    }
    common.atomic_json(command_path, record)
    return record


def tracked_clean_sources(
    paths: Sequence[Path],
    git: Path,
) -> tuple[Path, str]:
    if not paths:
        die("no tracked sources supplied")
    repo = Path(subprocess.check_output(
        (
            str(git), "--no-replace-objects", "-C", str(paths[0].parent),
            "rev-parse", "--show-toplevel",
        ),
        text=True,
    ).strip()).resolve(strict=True)
    for args, name in (
        (("diff", "--quiet", "HEAD", "--"), "worktree"),
        (("diff", "--cached", "--quiet", "HEAD", "--"), "index"),
    ):
        if subprocess.run(
                (str(git), "--no-replace-objects", "-C", str(repo), *args),
                check=False).returncode:
            die(f"tracked source {name} is dirty")
    head = subprocess.check_output(
        (
            str(git), "--no-replace-objects", "-C", str(repo),
            "rev-parse", "HEAD",
        ), text=True,
    ).strip()
    for path in paths:
        relative = path.relative_to(repo)
        tracked = subprocess.check_output(
            (
                str(git), "--no-replace-objects", "-C", str(repo),
                "cat-file", "blob", f"{head}:{relative}",
            ),
        )
        if tracked != common.stable_bytes(path):
            die(f"{relative} does not exactly match immutable HEAD")
    return repo, head


def q0_identity(
    binary: Path,
    binary_sha256: str,
    taskset: Path,
    taskset_record: dict[str, str],
) -> dict[str, Any]:
    ks = (4096, TARGET_K, 64000)
    fixed = [
        str(binary), *BASE_OPTIONS, "--N", ",".join(map(str, ks)),
        "--seed", SOURCE_SEEDS[1], "--schedule", "burst",
        "--packet-peel-seed-xor", "0", "--binary-dense-two-anchor",
    ]
    outputs: list[str] = []
    commands = (fixed, [*fixed, "--binary-dense-two-anchor-phase", "0"])
    pinned_commands = tuple(
        [str(taskset), "-c", "0", *command] for command in commands)
    options = (
        ("--binary-dense-two-anchor",),
        ("--binary-dense-two-anchor", "--binary-dense-two-anchor-phase", "0"),
    )
    rows: list[list[dict[str, str]]] = []
    for command, expected_options in zip(pinned_commands, options):
        common.verify_frozen_binary(binary, binary_sha256)
        if common.frozen_executable_path(
                taskset_record, "taskset") != taskset:
            die("frozen taskset changed before q0 identity launch")
        completed = subprocess.run(
            command, capture_output=True,
            check=False, timeout=300,
        )
        if completed.returncode != 0 or completed.stderr:
            die("q0 implicit/explicit identity command failed")
        output = completed.stdout.decode("utf-8")
        outputs.append(output)
        preamble = parse_preamble(output.splitlines()[0])
        if preamble.get("binary_dense_two_anchor_phase") != "0":
            die("q0 identity command did not report phase zero")
        rows.append(parse_bench_output(
            output, ks, "0x0", SOURCE_SEEDS[1], "burst",
            expected_options, True,
        ))
    comparisons = 0
    canonical_rows: list[dict[str, str]] = []
    for implicit, explicit in zip(*rows):
        if implicit.keys() != explicit.keys():
            die("q0 identity schemas differ")
        for field in implicit:
            if field in TIMING_FIELDS:
                continue
            comparisons += 1
            if implicit[field] != explicit[field]:
                die(f"q0 implicit/explicit identity differs at {field}")
        canonical_rows.append({
            field: implicit[field]
            for field in implicit if field not in TIMING_FIELDS
        })
    common.verify_frozen_binary(binary, binary_sha256)
    if common.frozen_executable_path(
            taskset_record, "taskset") != taskset:
        die("frozen taskset changed during q0 identity validation")
    return {
        "passed": True,
        "K": list(ks),
        "non_timing_field_comparisons": comparisons,
        "canonical_rows_sha256": sha256_bytes(json_bytes(canonical_rows)),
        "implicit_stdout_sha256": sha256_bytes(outputs[0].encode("utf-8")),
        "explicit_stdout_sha256": sha256_bytes(outputs[1].encode("utf-8")),
        "implicit_command": pinned_commands[0],
        "explicit_command": pinned_commands[1],
    }


def _prepare_with_fresh_pinned_benchmark_build(
    args: argparse.Namespace,
    final_result_dir: Path,
) -> dict[str, Any]:
    if not args.acknowledge_post_selection_only:
        die(
            "prepare requires --acknowledge-post-selection-only; this "
            "outcome-derived screen must not select a raw architecture"
        )
    script = Path(__file__).resolve(strict=True)
    allk = script.with_name("wh2_rank_floor_two_anchor_allk.py")
    helper = script.with_name("wh2_rank_floor_two_anchor_screen.py")
    git_record = common._trusted_executable("git")
    git = common.frozen_executable_path(git_record, "git")
    repo, head = tracked_clean_sources((script, allk, helper), git)
    taskset_record = common._trusted_executable("taskset")
    cmake_record = common._trusted_executable("cmake")
    python_runtime = common.python_runtime_identity()
    taskset = common.frozen_executable_path(taskset_record, "taskset")
    cmake = common.frozen_executable_path(cmake_record, "cmake")
    if args.source_cells.absolute().is_symlink():
        die("--source-cells must not be a symlink")
    cohort = load_source_cohort(args.source_cells.resolve(strict=True))
    if args.workers != 64 or args.build_workers < 1 or args.build_workers > 64:
        die("focused screen requires 64 workers and 1..64 build workers")
    thermal = args.thermal.resolve(strict=True)
    result_dir = args.result_dir.resolve()
    frozen = result_dir / "frozen"
    destinations = {
        "script": frozen / script.name,
        "allk": frozen / allk.name,
        "helper": frozen / helper.name,
        "binary": frozen / "wirehair_v2_bench",
        "cmake_cache": frozen / "CMakeCache.txt",
        "cohort": frozen / "cohort.json",
        "q0_identity": frozen / "q0_identity.json",
        "pinned_build": frozen / "pinned_build.json",
    }
    with common.fresh_pinned_benchmark_build(
            repo, head, args.binary, git_record, cmake_record,
            taskset_record, build_workers=args.build_workers,
            cpu_set="0-63") as fresh_build:
        build_policy = fresh_build.record["build_policy"]
        build_command = fresh_build.record["build_command"]
        live_binary_sha256 = fresh_build.record["binary_sha256"]
        identity = q0_identity(
            fresh_build.binary, live_binary_sha256, taskset,
            taskset_record)
        thermal_mark = thermal_start(thermal)
        post_repo, post_head = tracked_clean_sources(
            (script, allk, helper), git)
        if post_repo != repo or post_head != head:
            die("source repository changed during fresh benchmark build")
        if (common.frozen_executable_path(git_record, "git") != git or
                common.frozen_executable_path(
                    taskset_record, "taskset") != taskset or
                common.frozen_executable_path(
                    cmake_record, "cmake") != cmake):
            die("frozen build-tool path changed during fresh benchmark build")
        binary_bytes = common.stable_bytes(fresh_build.binary)
        cache_bytes = common.stable_bytes(fresh_build.cache)
        pinned_build_bytes = common.json_bytes(fresh_build.record)
        if common.sha256_bytes(binary_bytes) != live_binary_sha256:
            die("q0-tested fresh benchmark changed while it was captured")
        common.validate_pinned_build_record(
            fresh_build.record, fresh_build.cache)
    # Publish only after the fresh-build context's post-yield graph, source,
    # compiler, and tool checks have all succeeded.
    frozen.mkdir()
    common.atomic_write(destinations["binary"], binary_bytes)
    common.atomic_write(destinations["cmake_cache"], cache_bytes)
    common.atomic_write(destinations["pinned_build"], pinned_build_bytes)
    if common.sha256_file(destinations["binary"]) != live_binary_sha256:
        die("captured q0-tested benchmark changed while it was frozen")
    common.validate_pinned_build_record(
        json.loads(pinned_build_bytes), destinations["cmake_cache"])
    for source, destination in (
        (script, destinations["script"]),
        (allk, destinations["allk"]),
        (helper, destinations["helper"]),
    ):
        shutil.copyfile(source, destination)
    common.atomic_json(destinations["cohort"], cohort)
    common.atomic_json(destinations["q0_identity"], identity)
    for name in ("script", "binary"):
        destinations[name].chmod(0o755)
    common.verify_frozen_sources_at_commit(
        repo, head,
        (
            (script, destinations["script"]),
            (allk, destinations["allk"]),
            (helper, destinations["helper"]),
        ),
        git,
    )
    if common.sha256_file(destinations["binary"]) != live_binary_sha256:
        die("q0 identity binary changed while it was being frozen")
    pinned_build_sha256 = common.sha256_file(destinations["pinned_build"])
    common.validate_pinned_build_record(
        json.loads(common.stable_bytes(destinations["pinned_build"])),
        destinations["cmake_cache"],
    )
    contract = {
        "schema": "wirehair.wh2.two_anchor_phase.contract.v3",
        "purpose": "post-selection tuning only; forbidden for architecture selection",
        "source_commit": head,
        "source_repo": str(repo),
        "build_command": build_command,
        "binary_sha256": common.sha256_file(destinations["binary"]),
        "cmake_cache_sha256": common.sha256_file(destinations["cmake_cache"]),
        "pinned_build_sha256": pinned_build_sha256,
        "build_policy": build_policy,
        "taskset": taskset_record,
        "cmake": cmake_record,
        "git": git_record,
        "python": python_runtime,
        "script_sha256": common.sha256_file(destinations["script"]),
        "allk_sha256": common.sha256_file(destinations["allk"]),
        "helper_sha256": common.sha256_file(destinations["helper"]),
        "cohort_sha256": common.sha256_file(destinations["cohort"]),
        "q0_identity_sha256": common.sha256_file(destinations["q0_identity"]),
        "source_cells_sha256": SOURCE_CELLS_SHA256,
        "cohort_ks_sha256": COHORT_KS_SHA256,
        "thermal": str(thermal),
        "thermal_baseline": {
            "dev": thermal_mark["dev"], "ino": thermal_mark["ino"],
            "offset": thermal_mark["offset"],
            "edac_ce": thermal_mark["edac_ce"],
            "edac_ue": thermal_mark["edac_ue"],
            "monotonic_s": thermal_mark["monotonic_s"],
            "max_temperature_c": thermal_mark["max_temperature_c"],
            "row_sha256": sha256_bytes(thermal_mark["baseline_row"]),
        },
        "workers": args.workers,
        "timeout_seconds": args.timeout,
        "cpu_set": list(range(64)),
        "K_count": 226,
        "cutoff": CUTOFF,
        "arms": list(ARMS),
        "seeds": list(SEEDS),
        "schedules": list(SCHEDULES),
        "loss": "0.50",
        "trials": 1,
        "block_bytes": 64,
        "packet_peel_seed_xor": "0x0",
        "expected_jobs": 648,
        "expected_cells": 16272,
        "batching": {
            "small_partitions": 1,
            "large_partitions": LARGE_PARTITIONS,
            "balance_proxy": "sum(K)",
            "max_K_per_job": 27,
        },
        "saturated_timing_speed_claim_valid": False,
        "thermal_fields": list(THERMAL_FIELDS),
        "thermal_policy": {
            "limit_c": 90.0,
            "consecutive_samples": 3,
            "stale_seconds": 5.0,
            "min_cpu_busy_pct": 95.0,
        },
    }
    contract_path = frozen / "contract.json"
    common.atomic_json(contract_path, contract)
    staged_files = [*destinations.values(), contract_path]
    staged = frozen / "staged.sha256"
    common.atomic_write(staged, common.sha_manifest(result_dir, staged_files))
    verified = common.verify_sha_manifest(result_dir, staged)
    if set(verified) != {path.resolve(strict=True) for path in staged_files}:
        die("staged seal does not cover the exact frozen phase inputs")
    prepare_record = {
        "schema": "wirehair.wh2.two_anchor_phase.prepare.v2",
        "source_commit": head,
        "binary_sha256": contract["binary_sha256"],
        "pinned_build_sha256": pinned_build_sha256,
        "cohort_sha256": contract["cohort_sha256"],
        "q0_identity": identity,
        "staged_seal_sha256": common.sha256_file(staged),
        "result_dir": str(final_result_dir),
        "run_command": [
            python_runtime["path"],
            str(final_result_dir / "frozen" / script.name), "run",
            "--result-dir", str(final_result_dir),
        ],
    }
    common.atomic_json(result_dir / "prepare.json", prepare_record)
    return prepare_record


def prepare(args: argparse.Namespace) -> int:
    with common.prepare_result_staging(args.result_dir) as (staging, final):
        staged_args = argparse.Namespace(**vars(args))
        staged_args.result_dir = staging
        prepare_record = _prepare_with_fresh_pinned_benchmark_build(
            staged_args, final)
    print(canonical_json(prepare_record))
    return 0


def load_frozen(
    result_dir: Path,
) -> tuple[
    dict[str, Any], dict[str, Any], Path, Path, Path,
    tuple[int, int, int, int],
]:
    frozen = result_dir / "frozen"
    script = Path(__file__).resolve(strict=True)
    if script.parent != frozen.resolve(strict=True):
        die("run must use the frozen phase-screen script")
    expected_names = {
        script.name, "wh2_rank_floor_two_anchor_allk.py",
        "wh2_rank_floor_two_anchor_screen.py", "wirehair_v2_bench",
        "CMakeCache.txt", "cohort.json", "q0_identity.json",
        "pinned_build.json", "contract.json", "staged.sha256",
    }
    actual_names = {path.name for path in frozen.iterdir()}
    if (actual_names != expected_names or
            any(path.is_symlink() or not path.is_file()
                for path in frozen.iterdir())):
        die("frozen tree names or file types changed")
    prepare_anchor = common.validate_prepare_anchor(
        result_dir, "wirehair.wh2.two_anchor_phase.prepare.v2",
        ("source_commit", "binary_sha256", "pinned_build_sha256",
         "cohort_sha256", "q0_identity", "staged_seal_sha256",
         "run_command"),
        "staged_seal_sha256")
    staged = frozen / "staged.sha256"
    verified = common.verify_sha_manifest(result_dir, staged)
    expected_targets = {
        (frozen / name).resolve(strict=True)
        for name in expected_names if name != "staged.sha256"
    }
    if set(verified) != expected_targets:
        die("frozen staged seal target set mismatch")
    contract = json.loads(common.stable_bytes(frozen / "contract.json"))
    cohort = json.loads(common.stable_bytes(frozen / "cohort.json"))
    identity = json.loads(common.stable_bytes(frozen / "q0_identity.json"))
    pinned_build = json.loads(
        common.stable_bytes(frozen / "pinned_build.json"))
    common.validate_pinned_build_record(
        pinned_build, frozen / "CMakeCache.txt")
    if (not isinstance(contract, dict) or not isinstance(cohort, dict) or
            not isinstance(identity, dict)):
        die("frozen contract, cohort, and identity must be objects")
    python_runtime = contract.get("python")
    if python_runtime != common.python_runtime_identity():
        die("frozen phase-screen Python runtime changed")
    if prepare_anchor.get("run_command") != [
            python_runtime["path"], str(script), "run",
            "--result-dir", str(result_dir)]:
        die("phase-screen prepare run command changed")
    if (prepare_anchor.get("source_commit") != contract.get("source_commit") or
            prepare_anchor.get("binary_sha256") !=
            contract.get("binary_sha256") or
            prepare_anchor.get("pinned_build_sha256") !=
            contract.get("pinned_build_sha256") or
            prepare_anchor.get("cohort_sha256") !=
            contract.get("cohort_sha256") or
            prepare_anchor.get("q0_identity") != identity):
        die("prepare anchor differs from the frozen phase-screen inputs")
    if (contract.get("schema") != "wirehair.wh2.two_anchor_phase.contract.v3" or
            contract.get("purpose") !=
            "post-selection tuning only; forbidden for architecture selection" or
            contract.get("workers") != 64 or
            contract.get("cpu_set") != list(range(64)) or
            contract.get("K_count") != 226 or
            contract.get("expected_jobs") != 648 or
            contract.get("expected_cells") != 16272 or
            contract.get("source_cells_sha256") != SOURCE_CELLS_SHA256 or
            contract.get("cutoff") != CUTOFF or
            contract.get("arms") != list(ARMS) or
            contract.get("seeds") != list(SEEDS) or
            contract.get("schedules") != list(SCHEDULES) or
            contract.get("loss") != "0.50" or
            contract.get("trials") != 1 or
            contract.get("block_bytes") != 64 or
            contract.get("packet_peel_seed_xor") != "0x0" or
            contract.get("cohort_ks_sha256") != COHORT_KS_SHA256 or
            contract.get("batching") != {
                "small_partitions": 1,
                "large_partitions": LARGE_PARTITIONS,
                "balance_proxy": "sum(K)",
                "max_K_per_job": 27,
            } or
            contract.get("saturated_timing_speed_claim_valid") is not False or
            contract.get("thermal_fields") != list(THERMAL_FIELDS) or
            contract.get("thermal_policy") != {
                "limit_c": 90.0,
                "consecutive_samples": 3,
                "stale_seconds": 5.0,
                "min_cpu_busy_pct": 95.0,
            } or
            type(contract.get("timeout_seconds")) not in (int, float) or
            not math.isfinite(float(contract["timeout_seconds"])) or
            float(contract["timeout_seconds"]) <= 0 or
            not identity.get("passed")):
        die("frozen contract or q0 identity record mismatch")
    hash_fields = {
        "binary_sha256": "wirehair_v2_bench",
        "cmake_cache_sha256": "CMakeCache.txt",
        "pinned_build_sha256": "pinned_build.json",
        "script_sha256": script.name,
        "allk_sha256": "wh2_rank_floor_two_anchor_allk.py",
        "helper_sha256": "wh2_rank_floor_two_anchor_screen.py",
        "cohort_sha256": "cohort.json",
        "q0_identity_sha256": "q0_identity.json",
    }
    for field, name in hash_fields.items():
        if contract.get(field) != common.sha256_file(frozen / name):
            die(f"frozen contract {field} mismatch")
    common.validate_pinned_build_contract_binding(
        contract, pinned_build, "frozen phase-screen contract")
    if (cohort.get("schema") != "wirehair.wh2.two_anchor_phase.cohort.v1" or
            cohort.get("source_cells_sha256") != SOURCE_CELLS_SHA256 or
            cohort.get("ks_sha256") != COHORT_KS_SHA256):
        die("frozen cohort contract mismatch")
    binary = frozen / "wirehair_v2_bench"
    common.verify_frozen_binary(binary, contract.get("binary_sha256"))
    taskset = common.frozen_executable_path(
        contract.get("taskset"), "taskset")
    common.frozen_executable_path(contract.get("cmake"), "cmake")
    common.frozen_executable_path(contract.get("git"), "git")
    thermal = Path(contract["thermal"]).resolve(strict=True)
    thermal_identity = common.validate_frozen_thermal_source(
        thermal, contract.get("thermal_baseline"),
        float(contract["thermal_policy"]["stale_seconds"]),
        cpu_limit_c=float(contract["thermal_policy"]["limit_c"]),
        dimm_limit_c=float(contract["thermal_policy"]["limit_c"]),
        consecutive_limit=int(
            contract["thermal_policy"]["consecutive_samples"]))
    return contract, cohort, binary, taskset, thermal, thermal_identity


def milli(text: str, context: str) -> int:
    try:
        value = Decimal(text) * 1000
    except (InvalidOperation, ValueError):
        die(f"{context}: invalid decimal {text!r}")
    if not value.is_finite() or value != value.to_integral_value() or value < 0:
        die(f"{context}: decimal is not a nonnegative exact milli-value")
    return int(value)


def row_stats(row: dict[str, str], context: str) -> dict[str, int]:
    return {
        "rank_fail": int(row["rank_fail"], 10),
        "error": int(row["error"], 10),
        "heavy_shortfall": int(row["heavy_shortfall"], 10),
        "seed_attempt": int(row["seed_attempt"], 10),
        "inact_milli": milli(row["inact_mu"], f"{context}:inact"),
        "binary_def_milli": milli(
            row["binary_def_mu"], f"{context}:binary_def"
        ),
        "heavy_gain_milli": milli(
            row["heavy_gain_mu"], f"{context}:heavy_gain"
        ),
        "block_xors_milli": milli(
            row["block_xors_mu"], f"{context}:block_xors"
        ),
        "block_muladds_milli": milli(
            row["block_muladds_mu"], f"{context}:block_muladds"
        ),
    }


def load_results(
    result_dir: Path,
    contract: dict[str, Any],
    jobs: Sequence[Job],
    phase_hashes: dict[Path, str],
) -> tuple[dict[tuple[str, int, str, int], dict[str, str]], list[dict[str, Any]]]:
    result_root = result_dir / "results"
    taskset = common.frozen_executable_path(
        contract.get("taskset"), "taskset")
    expected_phase_targets = {
        (result_root / directory / filename).resolve(strict=True)
        for job in jobs
        for directory, filename in (
            ("commands", f"{job.stem}.json"),
            ("stdout", f"{job.stem}.csv"),
            ("stderr", f"{job.stem}.txt"),
        )
    }
    expected_phase_targets.update(
        (result_root / name).resolve(strict=True)
        for name in ("tasks.jsonl", "thermal_interval.csv", "run.json")
    )
    if set(phase_hashes) != expected_phase_targets:
        die("phase seal does not cover the exact focused artifact set")
    indexed: dict[tuple[str, int, str, int], dict[str, str]] = {}
    records: list[dict[str, Any]] = []
    for job in jobs:
        command_path = result_root / "commands" / f"{job.stem}.json"
        stdout_path = result_root / "stdout" / f"{job.stem}.csv"
        stderr_path = result_root / "stderr" / f"{job.stem}.txt"
        for path in (command_path, stdout_path, stderr_path):
            if path.resolve(strict=True) not in phase_hashes:
                die(f"job {job.job} artifact is outside phase seal")
        record = json.loads(common.stable_bytes(command_path))
        if not isinstance(record, dict):
            die(f"job {job.job} command record is not an object")
        command = [str(taskset), "-c", str(record.get("cpu")), *make_command(
            result_root.parent / "frozen/wirehair_v2_bench", job
        )]
        expected_stdout = str(stdout_path.relative_to(result_root))
        expected_stderr = str(stderr_path.relative_to(result_root))
        start_ns = record.get("start_ns")
        end_ns = record.get("end_ns")
        if (record.get("job") != job.job or record.get("arm") != job.arm or
                record.get("band") != job.band or
                record.get("seed_index") != job.seed_index or
                record.get("seed") != job.seed or
                record.get("schedule") != job.schedule or
                record.get("chunk") != job.chunk or
                record.get("K_count") != len(job.ks) or
                record.get("command") != command or
                record.get("returncode") != 0 or
                record.get("active_packet_peel_seed_xor") != "0x0" or
                type(record.get("cpu")) is not int or
                record.get("cpu") not in range(64) or
                record.get("stdout") != expected_stdout or
                record.get("stderr") != expected_stderr or
                record.get("stdout_sha256") != common.sha256_file(stdout_path) or
                record.get("stderr_sha256") != common.sha256_file(stderr_path) or
                type(start_ns) is not int or type(end_ns) is not int or
                start_ns <= 0 or end_ns < start_ns or
                record.get("elapsed_ns") != end_ns - start_ns or
                record.get("saturated_timing_speed_claim_valid") is not False or
                common.stable_bytes(stderr_path) != b""):
            die(f"job {job.job} command/provenance mismatch")
        stdout = common.stable_bytes(stdout_path).decode("utf-8")
        options = arm_options(job.arm, job.band)
        preamble = parse_preamble(stdout.splitlines()[0])
        if int(preamble.get("binary_dense_two_anchor_phase", "-1"), 10) != \
                expected_phase(options):
            die(f"job {job.job} sealed phase preamble mismatch")
        rows = parse_bench_output(
            stdout, job.ks, "0x0", job.seed, job.schedule, options, True,
        )
        for row in rows:
            key = (job.arm, job.seed_index, job.schedule, int(row["N"], 10))
            if key in indexed:
                die(f"duplicate result cell {key}")
            indexed[key] = row
        records.append(record)
    expected = len(ARMS) * len(SEEDS) * len(SCHEDULES) * 226
    if len(indexed) != expected or len(records) != 648:
        die("sealed result cell/job cardinality mismatch")
    records.sort(key=lambda record: record["job"])
    expected_tasks = b"".join(
        (canonical_json(record) + "\n").encode("utf-8") for record in records
    )
    if common.stable_bytes(result_root / "tasks.jsonl") != expected_tasks:
        die("tasks ledger is not the exact canonical ordered command ledger")
    run_record = json.loads(common.stable_bytes(result_root / "run.json"))
    if not isinstance(run_record, dict):
        die("sealed run summary is not an object")
    if (
        run_record.get("schema") != "wirehair.wh2.two_anchor_phase.run.v1" or
        run_record.get("source_commit") != contract["source_commit"] or
        run_record.get("binary_sha256") != contract["binary_sha256"] or
        run_record.get("jobs") != len(jobs) or
        run_record.get("workers") != 64 or
        run_record.get("recovery_cells") != 16272 or
        not isinstance(run_record.get("thermal"), dict) or
        run_record.get("saturated_timing_speed_claim_valid") is not False
    ):
        die("sealed run summary contract mismatch")
    return indexed, records


def is_failed(row: dict[str, str]) -> bool:
    return int(row["rank_fail"], 10) != 0 or int(row["error"], 10) != 0


def comparison(
    indexed: dict[tuple[str, int, str, int], dict[str, str]],
    ks: Sequence[int], base_arm: str, candidate_arm: str,
) -> dict[str, Any]:
    repairs = introductions = both_fail = paired_success = 0
    base_xors = cand_xors = base_muladds = cand_muladds = 0
    for seed_index in range(len(SEEDS)):
        for schedule in SCHEDULES:
            for K in ks:
                base = indexed[(base_arm, seed_index, schedule, K)]
                cand = indexed[(candidate_arm, seed_index, schedule, K)]
                bf, cf = is_failed(base), is_failed(cand)
                repairs += bf and not cf
                introductions += not bf and cf
                both_fail += bf and cf
                if not bf and not cf:
                    paired_success += 1
                    base_stats = row_stats(base, "comparison base")
                    cand_stats = row_stats(cand, "comparison candidate")
                    base_xors += base_stats["block_xors_milli"]
                    cand_xors += cand_stats["block_xors_milli"]
                    base_muladds += base_stats["block_muladds_milli"]
                    cand_muladds += cand_stats["block_muladds_milli"]
    return {
        "repairs": repairs,
        "introductions": introductions,
        "net_failures_removed": repairs - introductions,
        "both_fail": both_fail,
        "paired_success_cells": paired_success,
        "block_xor_ratio_paired_success": (
            str(Decimal(cand_xors) / Decimal(base_xors)) if base_xors else None
        ),
        "block_muladd_ratio_paired_success": (
            str(Decimal(cand_muladds) / Decimal(base_muladds))
            if base_muladds else None
        ),
    }


def fail_counts(
    indexed: dict[tuple[str, int, str, int], dict[str, str]],
    ks: Sequence[int], arm: str, seed_indices: Sequence[int],
) -> dict[int, int]:
    return {
        K: sum(
            is_failed(indexed[(arm, seed_index, schedule, K)])
            for seed_index in seed_indices for schedule in SCHEDULES
        )
        for K in ks
    }


def analyze(
    result_dir: Path,
    contract: dict[str, Any],
    cohort: dict[str, Any],
    jobs: Sequence[Job],
    phase_hashes: dict[Path, str],
) -> dict[str, Any]:
    ks = [int(K) for K in cohort["ks"]]
    indexed, records = load_results(result_dir, contract, jobs, phase_hashes)
    if len(records) != len(jobs):
        die("analysis did not load the exact job record count")
    source_comparisons = 0
    seen_references: set[tuple[int, int, str]] = set()
    for reference in cohort["source_references"]:
        source_seed_index = int(reference["source_seed_index"])
        if (not 0 <= source_seed_index < len(SOURCE_SEEDS) or
                reference.get("seed") != SOURCE_SEEDS[source_seed_index] or
                reference.get("schedule") not in SCHEDULES):
            die("source replay reference stratum mismatch")
        seed_index = len(SEEDS) - len(SOURCE_SEEDS) + source_seed_index
        K = int(reference["K"])
        schedule = reference["schedule"]
        reference_key = (K, source_seed_index, schedule)
        if K not in ks or reference_key in seen_references:
            die("source replay reference K is absent or duplicated")
        seen_references.add(reference_key)
        for arm in ("d12", "q0"):
            stats = row_stats(indexed[(arm, seed_index, schedule, K)], "source replay")
            for field in STATS:
                source_comparisons += 1
                if stats[field] != int(reference[f"{arm}_{field}"]):
                    die(
                        f"source replay mismatch arm={arm} seed={seed_index} "
                        f"schedule={schedule} K={K} field={field}"
                    )
    expected_references = len(ks) * len(SOURCE_SEEDS) * len(SCHEDULES)
    if len(seen_references) != expected_references:
        die("source replay reference grid is incomplete")

    low_ks = [K for K in ks if K < CUTOFF]
    low_comparisons = 0
    for seed_index in range(len(SEEDS)):
        for schedule in SCHEDULES:
            for K in low_ks:
                base = indexed[("d12", seed_index, schedule, K)]
                for arm in ARMS[1:]:
                    candidate = indexed[(arm, seed_index, schedule, K)]
                    if base.keys() != candidate.keys():
                        die("low-K phase identity schemas differ")
                    for field in base:
                        if field in TIMING_FIELDS:
                            continue
                        low_comparisons += 1
                        if base[field] != candidate[field]:
                            die(
                                f"low-K identity mismatch arm={arm} "
                                f"seed={seed_index} schedule={schedule} "
                                f"K={K} field={field}"
                            )

    analysis = result_dir / "analysis"
    analysis.mkdir(exist_ok=False)
    paired_header = ["K", "schedule", "seed_index", "seed"] + [
        f"{arm}_{field}" for arm in ARMS for field in STATS
    ]
    paired_buffer = io.StringIO(newline="")
    paired_writer = csv.DictWriter(
        paired_buffer, fieldnames=paired_header, lineterminator="\n"
    )
    paired_writer.writeheader()
    for seed_index, seed in enumerate(SEEDS):
        for schedule in SCHEDULES:
            for K in ks:
                output: dict[str, Any] = {
                    "K": K, "schedule": schedule,
                    "seed_index": seed_index, "seed": seed,
                }
                for arm in ARMS:
                    stats = row_stats(indexed[(arm, seed_index, schedule, K)], "paired")
                    output.update({f"{arm}_{field}": stats[field] for field in STATS})
                paired_writer.writerow(output)
    paired_path = analysis / "paired_cells.csv"
    common.atomic_write(paired_path, paired_buffer.getvalue().encode("utf-8"))

    all_indices = tuple(range(len(SEEDS)))
    # SEEDS deliberately puts the three previously unseen seeds first and
    # the three source-campaign (already revealed) seeds last.  Keep these
    # labels tied to that frozen ordering: reversing them would silently use
    # the development outcomes as the nominal fresh-seed holdout.
    fresh_indices = tuple(range(3))
    revealed_indices = tuple(range(3, 6))
    cohorts = {
        "all6": all_indices,
        "fresh3": fresh_indices,
        "source_revealed3": revealed_indices,
    }
    counts = {
        cohort_name: {
            arm: fail_counts(indexed, ks, arm, indices)
            for arm in ARMS
        }
        for cohort_name, indices in cohorts.items()
    }
    resonance: dict[str, Any] = {}
    for cohort_name in cohorts:
        cohort_report: dict[str, Any] = {}
        base_counts = counts[cohort_name]["d12"]
        for arm in ARMS:
            arm_counts = counts[cohort_name][arm]
            cohort_report[arm] = {
                "multi_stratum_k": [K for K in ks if arm_counts[K] >= 2],
                "new_vs_d12_multi_stratum_k": [
                    K for K in ks if base_counts[K] < 2 and arm_counts[K] >= 2
                ],
            }
        resonance[cohort_name] = cohort_report

    arm_summary: dict[str, Any] = {}
    for arm in ARMS:
        rows = [
            indexed[(arm, seed_index, schedule, K)]
            for seed_index in range(len(SEEDS))
            for schedule in SCHEDULES for K in ks
        ]
        arm_summary[arm] = {
            "cells": len(rows),
            "failures": sum(is_failed(row) for row in rows),
            "errors": sum(int(row["error"], 10) for row in rows),
        }

    post_selection_candidate_gates: dict[str, Any] = {}
    for arm in ("q1", "q2"):
        all_counts = counts["all6"][arm]
        fresh_counts = counts["fresh3"][arm]
        retained = {
            str(K): fresh_counts[K] == 0 for K in RETAINED_CLEAR_KS
        }
        q0_comparison = comparison(indexed, ks, "q0", arm)
        post_selection_candidate_gates[arm] = {
            "K46252_all6_fail_strata": all_counts[TARGET_K],
            "K46252_fresh3_fail_strata": fresh_counts[TARGET_K],
            "K46252_lt2": all_counts[TARGET_K] < 2,
            "retained_q0_clears_on_fresh3": retained,
            "no_new_all6_multi_stratum_resonance_vs_d12": not resonance[
                "all6"
            ][arm]["new_vs_d12_multi_stratum_k"],
            "no_new_fresh3_multi_stratum_resonance_vs_d12": not resonance[
                "fresh3"
            ][arm]["new_vs_d12_multi_stratum_k"],
            "net_not_worse_vs_q0": q0_comparison["net_failures_removed"] >= 0,
            "all_recovery_gates": (
                all_counts[TARGET_K] < 2 and all(retained.values()) and
                not resonance["all6"][arm]["new_vs_d12_multi_stratum_k"] and
                not resonance["fresh3"][arm]["new_vs_d12_multi_stratum_k"] and
                q0_comparison["net_failures_removed"] >= 0
            ),
        }

    summary = {
        "schema": "wirehair.wh2.two_anchor_phase.analysis.v1",
        "source_commit": contract["source_commit"],
        "binary_sha256": contract["binary_sha256"],
        "cohort_ks_sha256": COHORT_KS_SHA256,
        "K_count": len(ks),
        "jobs": len(jobs),
        "cells": len(ARMS) * len(SEEDS) * len(SCHEDULES) * len(ks),
        "arms": arm_summary,
        "comparisons": {
            f"{candidate}_vs_{base}": comparison(
                indexed, ks, base, candidate
            )
            for base in ("d12", "q0") for candidate in ARMS[1:]
            if candidate != base
        },
        "source_replay_identity": {
            "passed": True,
            "cells_per_arm": len(cohort["source_references"]),
            "non_timing_stat_comparisons": source_comparisons,
        },
        "low_K_identity": {
            "passed": True,
            "K": low_ks,
            "non_timing_field_comparisons": low_comparisons,
        },
        "resonance_threshold_strata": 2,
        "resonance": resonance,
        "post_selection_candidate_gates": post_selection_candidate_gates,
        "saturated_timing_speed_claim_valid": False,
        "phase_seal_sha256": common.sha256_file(
            result_dir / "results/phase_complete.sha256"
        ),
    }
    summary_path = analysis / "summary.json"
    common.atomic_json(summary_path, summary)
    analysis_files = [paired_path, summary_path]
    seal = analysis / "analysis_complete.sha256"
    common.atomic_write(seal, common.sha_manifest(result_dir, analysis_files))
    verified = common.verify_sha_manifest(result_dir, seal)
    if set(verified) != {path.resolve(strict=True) for path in analysis_files}:
        die("analysis seal does not cover exact phase-screen artifacts")
    return summary


def run(args: argparse.Namespace) -> int:
    result_dir = args.result_dir.resolve(strict=True)
    if set(path.name for path in result_dir.iterdir()) != {"frozen", "prepare.json"}:
        die("result directory is not pristine after prepare")
    contract, cohort, binary, taskset, thermal, thermal_identity = \
        load_frozen(result_dir)
    jobs = build_jobs([int(K) for K in cohort["ks"]])
    cpus = sorted(os.sched_getaffinity(0))
    if cpus != list(range(64)):
        die(f"focused run requires exact process affinity 0..63, got {cpus}")
    result_root = result_dir / "results"
    result_root.mkdir(exist_ok=False)
    for name in ("stdout", "stderr", "commands"):
        (result_root / name).mkdir()
    pool = common.CpuPool(cpus, 64)
    abort = threading.Event()
    guard = common.ThermalGuard(
        thermal, abort,
        stale_seconds=float(contract["thermal_policy"]["stale_seconds"]),
        limit_c=float(contract["thermal_policy"]["limit_c"]),
        consecutive_limit=int(contract["thermal_policy"]["consecutive_samples"]),
        expected_dev=thermal_identity[0], expected_ino=thermal_identity[1],
        expected_edac_ce=thermal_identity[2],
        expected_edac_ue=thermal_identity[3],
    )
    if common.validate_frozen_thermal_source(
            thermal, contract.get("thermal_baseline"),
            float(contract["thermal_policy"]["stale_seconds"]),
            cpu_limit_c=float(contract["thermal_policy"]["limit_c"]),
            dimm_limit_c=float(contract["thermal_policy"]["limit_c"]),
            consecutive_limit=int(
                contract["thermal_policy"]["consecutive_samples"])) != \
            thermal_identity:
        die("frozen thermal baseline changed before phase guard start")
    registry = common.ProcessRegistry()
    records: list[dict[str, Any]] = []
    campaign_error: BaseException | None = None
    thermal_summary: dict[str, Any] | None = None
    with common.CampaignSignalGuard(abort, registry):
        try:
            guard.start()
            with concurrent.futures.ThreadPoolExecutor(max_workers=64) as executor:
                futures: list[concurrent.futures.Future[dict[str, Any]]] = []
                try:
                    for job in jobs:
                        if abort.is_set():
                            die("campaign aborted during job submission")
                        futures.append(executor.submit(
                            run_job, job, binary, contract["binary_sha256"],
                            taskset, contract["taskset"], result_root, pool,
                            abort, registry,
                            float(contract["timeout_seconds"]),
                        ))
                    for completed, future in enumerate(
                            concurrent.futures.as_completed(futures), 1):
                        records.append(future.result())
                        if completed % 100 == 0 or completed == len(jobs):
                            print(
                                f"phase-screen progress {completed}/{len(jobs)}",
                                flush=True,
                            )
                        if guard.error:
                            die(f"thermal guard failed: {guard.error}")
                except BaseException as exc:
                    campaign_error = exc
                    abort.set()
                    registry.signal_all(signal.SIGTERM)
                    for future in futures:
                        future.cancel()
        finally:
            if guard.started:
                thermal_summary = guard.finish(result_root / "thermal_interval.csv")
    if guard.error:
        die(f"thermal guard failed: {guard.error}")
    if thermal_summary is None:
        if campaign_error is not None:
            raise campaign_error
        die("thermal guard did not produce a summary")
    if (
        thermal_summary["thermal_high_max_consecutive_samples"] >=
            int(contract["thermal_policy"]["consecutive_samples"]) or
        thermal_summary["dimm_read_errors_max"] != 0 or
        thermal_summary["edac_ce_delta"] != 0 or
        thermal_summary["edac_ue_delta"] != 0
    ):
        die(f"memory telemetry changed during campaign: {thermal_summary}")
    if campaign_error is not None:
        raise campaign_error
    # Utilization is a validity gate for a complete saturated campaign, not
    # for an intentionally aborted run.  Checking it after campaign_error
    # preserves the causal job failure while the safety telemetry above still
    # takes priority.
    if thermal_summary["cpu_busy_min_pct"] < \
            float(contract["thermal_policy"]["min_cpu_busy_pct"]):
        die(f"CPU utilization fell below campaign gate: {thermal_summary}")
    common.verify_frozen_binary(binary, contract.get("binary_sha256"))
    if common.frozen_executable_path(
            contract.get("taskset"), "taskset") != taskset:
        die("frozen taskset path changed during the phase screen")
    if common.validate_frozen_thermal_source(
            thermal, contract.get("thermal_baseline"),
            float(contract["thermal_policy"]["stale_seconds"]),
            cpu_limit_c=float(contract["thermal_policy"]["limit_c"]),
            dimm_limit_c=float(contract["thermal_policy"]["limit_c"]),
            consecutive_limit=int(
                contract["thermal_policy"]["consecutive_samples"])) != \
            thermal_identity:
        die("frozen thermal baseline changed during the phase screen")
    records.sort(key=lambda record: record["job"])
    if len(records) != len(jobs) or registry.count() != 0:
        die("focused campaign did not reap exact job set")
    tasks_data = b"".join(
        (canonical_json(record) + "\n").encode("utf-8") for record in records
    )
    common.atomic_write(result_root / "tasks.jsonl", tasks_data)
    run_summary = {
        "schema": "wirehair.wh2.two_anchor_phase.run.v1",
        "source_commit": contract["source_commit"],
        "binary_sha256": contract["binary_sha256"],
        "jobs": len(jobs),
        "workers": 64,
        "recovery_cells": 16272,
        "thermal": thermal_summary,
        "saturated_timing_speed_claim_valid": False,
    }
    common.atomic_json(result_root / "run.json", run_summary)
    phase_files = [
        path for path in result_root.rglob("*")
        if path.is_file() and path.name != "phase_complete.sha256"
    ]
    phase_seal = result_root / "phase_complete.sha256"
    common.atomic_write(phase_seal, common.sha_manifest(result_dir, phase_files))
    phase_hashes = common.verify_sha_manifest(result_dir, phase_seal)
    summary = analyze(result_dir, contract, cohort, jobs, phase_hashes)
    complete = {
        "result_dir": str(result_dir),
        "source_commit": contract["source_commit"],
        "binary_sha256": contract["binary_sha256"],
        "staged_seal_sha256": common.sha256_file(
            result_dir / "frozen/staged.sha256"
        ),
        "phase_seal_sha256": common.sha256_file(phase_seal),
        "analysis_seal_sha256": common.sha256_file(
            result_dir / "analysis/analysis_complete.sha256"
        ),
        "summary": summary,
    }
    common.atomic_json(result_dir / "complete.json", complete)
    print(canonical_json(complete))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="mode", required=True)
    prepare_parser = subparsers.add_parser("prepare")
    prepare_parser.add_argument(
        "--binary", type=Path, required=True,
        help=("accepted local benchmark used only to confirm the system "
              "toolchain; its binary and generated graph are never reused"))
    prepare_parser.add_argument("--source-cells", type=Path, required=True)
    prepare_parser.add_argument("--thermal", type=Path, required=True)
    prepare_parser.add_argument("--result-dir", type=Path, required=True)
    prepare_parser.add_argument("--workers", type=int, default=64)
    prepare_parser.add_argument("--build-workers", type=int, default=64)
    prepare_parser.add_argument("--timeout", type=float, default=900.0)
    prepare_parser.add_argument(
        "--acknowledge-post-selection-only", action="store_true"
    )
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--result-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.mode == "prepare":
        if not math.isfinite(args.timeout) or args.timeout <= 0:
            die("--timeout must be finite and positive")
        return prepare(args)
    return run(args)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except common.CampaignError as exc:
        print(f"phase-screen error: {exc}", file=sys.stderr)
        raise SystemExit(1)
