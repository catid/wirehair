#!/usr/bin/env python3
"""Prepare, execute, and strictly reduce the bounded WH2 row-mask screen.

The preparation path is intentionally separate from execution.  It freezes an
already-built test-hook benchmark, the sole CPU/DIMM thermal sampler, the
validated all-K source receipt, and a complete task ledger.  ``launch
--preflight-only`` performs no benchmark work.  ``launch --execute`` is the
only path that starts jobs, and it refuses a second SPD/I2C sampler.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import hashlib
import heapq
import io
import itertools
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import threading
import time
from collections import defaultdict
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation, getcontext
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Set, Tuple


sys.dont_write_bytecode = True
getcontext().prec = 50

SCHEMA = "wirehair.wh2.grouped_row_mask.v1"
CODEC_IMPLEMENTATION_HEAD = "42fb578b3df159970708c2db51684a9ad1abf93c"
SOURCE_HEAD = "0978602bb535712f6136b94c298ef428e4883fbc"
SOURCE_SUMMARY_SHA256 = (
    "c09d1e9b74a0e6cc38fea67d7038a3a77281736d7fd1be4d2867c48d6e3a2e72"
)
SOURCE_SUMMARY_NAME = "validated_summary.json"
SAMPLER_NAME = "wirehair_expo_thermal_sampler.py"
BINARY_NAME = "wirehair_v2_bench"
CONTROLLER_NAME = "wh2_grouped_row_mask_campaign.py"

SEEDS = (
    "0xd1b54a32d192ed03",
    "0x94d049bb133111eb",
    "0x8538ecb5bd456ea3",
)
SCHEDULES = ("burst", "adversarial", "repair-only")
PANELS = (
    {"name": "p48_r3", "period": 48, "rows": 3, "suffix_mask": 0x380},
    {"name": "p32_r7", "period": 32, "rows": 7, "suffix_mask": 0x3F8},
)
SOURCE_ARMS = ("p48_r3", "p32_r7", "p48_r0", "prod244")
FINALIST_ARMS = ("p48_r3", "p32_r7")
K_LO = 2
K_HI = 64_000
CONTROL_COUNT = 2048
MASK_COUNT = 120
CHUNK_MAX = 256
DEFAULT_WORKERS = 126
# The saturation floor is a coarse tripwire for a collapsed or undersized
# worker pool, not a precision gate: at 126 workers on a 128-CPU host the
# sustained ceiling is 98.4% and the measured whole-segment mean lands
# near 95% once the ramp and drain tails are included, so a 95 floor is
# flaky by construction while 90 still trips when saturation is lost.
CPU_BUSY_FLOOR = Decimal("90")
# The floor is a control-stage guarantee: hard segments are sub-minute
# bursts of one- and two-cell jobs whose interval busy mean is dominated
# by worker-pool spawn latency.  Every segment still passes the coverage,
# cadence, DIMM, EDAC and sampler-exit gates.
CPU_BUSY_FLOOR_STAGES = ("control",)
THERMAL_INTERVAL_SECONDS = Decimal("1")
THERMAL_MIN_GAP_SECONDS = Decimal("0.5")
THERMAL_MAX_GAP_SECONDS = Decimal("2.5")
THERMAL_READY_SAMPLES = 2
THERMAL_READY_TIMEOUT_SECONDS = 12.0
THERMAL_END_TIMEOUT_SECONDS = 5.0
BUILD_TARGET = "wirehair_v2_bench"
BUILD_CONFIGURE_ARGS = (
    "-G", "Ninja",
    "-DCMAKE_BUILD_TYPE=Release",
    "-DBUILD_TESTS=ON",
    "-DBUILD_CODEC_V2=ON",
    "-DMARCH_NATIVE=ON",
    "-DWIREHAIR_BUILD_BENCHMARKS=ON",
    "-DWIREHAIR_STRICT_WARNINGS=ON",
    "-DWH_LTO=OFF",
    "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
)

CSV_FIELDS = (
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
)
THERMAL_FIELDS = (
    "utc", "monotonic_s", "cpu_busy_pct", "cpu_avg_mhz", "cpu_tctl_c",
    "dimm_i2c1_50_c", "dimm_i2c1_51_c", "dimm_i2c1_52_c",
    "dimm_i2c1_53_c", "dimm_i2c2_50_c", "dimm_i2c2_51_c",
    "dimm_i2c2_52_c", "dimm_i2c2_53_c", "dimm_read_errors",
    "load1", "load5", "load15", "edac_ce", "edac_ue",
)


class CampaignError(RuntimeError):
    pass


def die(message: str) -> None:
    raise CampaignError(message)


def canonical_json(value: object) -> bytes:
    return (json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False,
    ) + "\n").encode("ascii")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_once(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_symlink():
        die("immutable artifact is a symlink {}".format(path))
    if path.exists():
        if path.read_bytes() != data:
            die("refusing to replace immutable artifact {}".format(path))
        return
    temporary = Path(
        str(path) + ".part.{}.{}".format(os.getpid(), threading.get_ident())
    )
    try:
        with temporary.open("xb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        try:
            os.link(str(temporary), str(path))
        except FileExistsError:
            if path.read_bytes() != data:
                die("refusing to replace immutable artifact {}".format(path))
    finally:
        try:
            temporary.unlink()
        except FileNotFoundError:
            pass


def json_lines(values: Iterable[object]) -> bytes:
    return b"".join(canonical_json(value) for value in values)


def load_json(path: Path) -> Dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="ascii"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        die("cannot read {}: {}".format(path, exc))
    if not isinstance(value, dict):
        die("{} is not a JSON object".format(path))
    return value


def run_captured(
    argv: Sequence[str],
    *,
    cwd: Path,
    timeout: float,
    environment: Optional[Mapping[str, str]] = None,
    label: str,
) -> Tuple[bytes, bytes]:
    try:
        result = subprocess.run(
            list(argv), cwd=str(cwd), env=None if environment is None else dict(environment),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        die("{} failed to run: {}".format(label, exc))
    if result.returncode:
        detail = result.stderr[-4096:].decode("utf-8", "replace")
        die("{} exited {}: {}".format(label, result.returncode, detail))
    return result.stdout, result.stderr


def git_output(repo: Path, *args: str) -> bytes:
    stdout, stderr = run_captured(
        ("git", *args), cwd=repo, timeout=60.0, label="git " + " ".join(args),
    )
    if stderr:
        die("git {} produced stderr".format(" ".join(args)))
    return stdout


def clean_source_identity(repo: Path) -> Dict[str, object]:
    status = git_output(repo, "status", "--porcelain=v1", "--untracked-files=all")
    if status:
        die("campaign preparation requires a completely clean source worktree")
    head = git_output(repo, "rev-parse", "HEAD").decode("ascii").strip()
    tree = git_output(repo, "rev-parse", "HEAD^{tree}").decode("ascii").strip()
    if re.fullmatch(r"[0-9a-f]{40}", head) is None or \
            re.fullmatch(r"[0-9a-f]{40}", tree) is None:
        die("source commit/tree identity is malformed")
    ancestor = subprocess.run(
        ["git", "merge-base", "--is-ancestor", CODEC_IMPLEMENTATION_HEAD, head],
        cwd=str(repo), stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    if ancestor.returncode != 0:
        die("campaign source does not contain the reviewed row-mask implementation")
    manifest = git_output(repo, "ls-tree", "-r", "--full-tree", head)
    if not manifest or not manifest.endswith(b"\n"):
        die("source tree manifest is empty or noncanonical")
    return {
        "commit": head,
        "tree_oid": tree,
        "status_porcelain_sha256": sha256_bytes(status),
        "tree_manifest_sha256": sha256_bytes(manifest),
        "tree_manifest_bytes": len(manifest),
    }


def executable_identity(name: str, version_args: Sequence[str]) -> Dict[str, object]:
    located = shutil.which(name)
    if not located:
        die("required build tool is unavailable: {}".format(name))
    path = Path(located).resolve(strict=True)
    stdout, stderr = run_captured(
        (str(path), *version_args), cwd=Path.cwd(), timeout=30.0,
        label="{} version".format(name),
    )
    return {
        "path": str(path), "sha256": sha256_file(path),
        "version_stdout_sha256": sha256_bytes(stdout),
        "version_stderr_sha256": sha256_bytes(stderr),
    }


def parse_cmake_cache_value(cache: Path, key: str) -> str:
    prefix = key + ":"
    matches: List[str] = []
    for line in cache.read_text(encoding="utf-8").splitlines():
        if line.startswith(prefix) and "=" in line:
            matches.append(line.split("=", 1)[1])
    if len(matches) != 1:
        die("CMake cache lacks one {} entry".format(key))
    return matches[0]


def fresh_build_binary(repo: Path, frozen: Path, workers: int) -> Dict[str, object]:
    """Build the benchmark from a clean committed tree and freeze provenance."""
    before = clean_source_identity(repo)
    tools = {
        "cmake": executable_identity("cmake", ("--version",)),
        "ninja": executable_identity("ninja", ("--version",)),
        "readelf": executable_identity("readelf", ("--version",)),
    }
    manifest = git_output(repo, "ls-tree", "-r", "--full-tree", str(before["commit"]))
    archive = git_output(repo, "archive", "--format=tar", str(before["commit"]))
    source_record = dict(before)
    source_record["source_archive_sha256"] = sha256_bytes(archive)
    build_evidence = frozen / "build"
    build_evidence.mkdir()
    (build_evidence / "source_tree.txt").write_bytes(manifest)
    with tempfile.TemporaryDirectory(prefix="wirehair-row-mask-build-") as temporary:
        source = Path(temporary) / "source"
        source.mkdir()
        with tarfile.open(fileobj=io.BytesIO(archive), mode="r:") as stream:
            stream.extractall(path=source, filter="data")
        build = Path(temporary) / "build"
        configure = [
            str(tools["cmake"]["path"]), "-S", str(source), "-B", str(build),
            *BUILD_CONFIGURE_ARGS,
        ]
        environment = dict(os.environ)
        environment.update({"LANG": "C", "LC_ALL": "C", "TZ": "UTC"})
        configure_stdout, configure_stderr = run_captured(
            configure, cwd=Path(temporary), timeout=600.0,
            environment=environment, label="fresh campaign CMake configure",
        )
        build_command = [
            str(tools["cmake"]["path"]), "--build", str(build),
            "--target", BUILD_TARGET, "--parallel", str(workers), "--verbose",
        ]
        build_stdout, build_stderr = run_captured(
            build_command, cwd=Path(temporary), timeout=1800.0,
            environment=environment, label="fresh campaign benchmark build",
        )
        candidate = build / "codec" / BINARY_NAME
        cache = build / "CMakeCache.txt"
        graph = build / "build.ninja"
        compile_commands = build / "compile_commands.json"
        for path, label in ((candidate, "built binary"), (cache, "CMake cache"),
                            (graph, "Ninja graph"),
                            (compile_commands, "compile commands")):
            if not path.is_file() or path.is_symlink():
                die("fresh build lacks a regular {}".format(label))
        if not os.access(str(candidate), os.X_OK):
            die("fresh benchmark is not executable")
        if Path(parse_cmake_cache_value(cache, "CMAKE_HOME_DIRECTORY")).resolve() != source:
            die("fresh CMake cache is not bound to the campaign source")
        compiler = Path(parse_cmake_cache_value(cache, "CMAKE_CXX_COMPILER")).resolve(strict=True)
        compiler_stdout, compiler_stderr = run_captured(
            (str(compiler), "--version"), cwd=Path(temporary), timeout=30.0,
            label="C++ compiler version",
        )
        notes_stdout, notes_stderr = run_captured(
            (str(tools["readelf"]["path"]), "--notes", str(candidate)),
            cwd=Path(temporary), timeout=30.0, label="ELF build-ID inspection",
        )
        build_ids = re.findall(
            rb"^[ \t]*Build ID: ([0-9a-f]+)[ \t]*$", notes_stdout, re.MULTILINE,
        )
        if notes_stderr or len(build_ids) != 1:
            die("fresh benchmark lacks exactly one lowercase ELF build ID")
        evidence_sources = {
            "CMakeCache.txt": cache,
            "build.ninja": graph,
            "compile_commands.json": compile_commands,
        }
        for name, source in evidence_sources.items():
            shutil.copyfile(str(source), str(build_evidence / name))
        log_bytes = {
            "configure.stdout": configure_stdout,
            "configure.stderr": configure_stderr,
            "build.stdout": build_stdout,
            "build.stderr": build_stderr,
            "compiler.stdout": compiler_stdout,
            "compiler.stderr": compiler_stderr,
            "readelf-notes.stdout": notes_stdout,
            "readelf-notes.stderr": notes_stderr,
        }
        for name, data in log_bytes.items():
            (build_evidence / name).write_bytes(data)
        shutil.copyfile(str(candidate), str(frozen / BINARY_NAME))
        os.chmod(str(frozen / BINARY_NAME), 0o555)

    after = clean_source_identity(repo)
    if after != before:
        die("source commit/tree/worktree changed during the fresh build")
    evidence_files = sorted(path for path in build_evidence.iterdir() if path.is_file())
    record = {
        "schema": SCHEMA + ".fresh_build",
        "source": source_record,
        "configure_command": configure,
        "build_command": build_command,
        "workers": workers,
        "tools": tools,
        "compiler": {
            "path": str(compiler), "sha256": sha256_file(compiler),
            "version_stdout_sha256": sha256_bytes(compiler_stdout),
            "version_stderr_sha256": sha256_bytes(compiler_stderr),
        },
        "binary_sha256": sha256_file(frozen / BINARY_NAME),
        "binary_elf_build_id": build_ids[0].decode("ascii"),
        "evidence": {
            path.name: {"sha256": sha256_file(path), "bytes": path.stat().st_size}
            for path in evidence_files
        },
        "fresh_build": True,
    }
    write_once(frozen / "build_receipt.json", canonical_json(record))
    return record


def cell_tuple(record: Mapping[str, object]) -> Tuple[int, int, str]:
    try:
        key = (int(record["K"]), int(record["seed_index"]), str(record["schedule"]))
    except (KeyError, TypeError, ValueError) as exc:
        die("malformed all-K failure record: {}".format(exc))
    if not (K_LO <= key[0] <= K_HI and 0 <= key[1] < len(SEEDS) and
            key[2] in SCHEDULES):
        die("out-of-domain all-K failure key {}".format(key))
    return key


def cell_record(key: Tuple[int, int, str]) -> Dict[str, object]:
    return {"K": key[0], "seed_index": key[1], "schedule": key[2]}


def validate_source_summary(summary: Mapping[str, object]) -> Dict[str, List[dict]]:
    exact = {
        "schema": 3,
        "head": SOURCE_HEAD,
        "K_range": [K_LO, K_HI],
        "K_count": K_HI - K_LO + 1,
        "cells_per_arm": (K_HI - K_LO + 1) * len(SEEDS) * len(SCHEDULES),
        "validation_issue_count": 0,
        "timing_promotional": False,
    }
    for key, expected in exact.items():
        if summary.get(key) != expected:
            die("all-K source receipt {} changed".format(key))
    arms = summary.get("arms")
    if not isinstance(arms, dict) or set(arms) != set(SOURCE_ARMS):
        die("all-K source arm set changed")
    records: Dict[str, List[dict]] = {}
    for arm in SOURCE_ARMS:
        value = arms.get(arm)
        if not isinstance(value, dict) or not isinstance(value.get("failure_records"), list):
            die("all-K source arm {} lacks failure records".format(arm))
        seen: Set[Tuple[int, int, str]] = set()
        parsed: List[dict] = []
        for raw in value["failure_records"]:
            if not isinstance(raw, dict):
                die("all-K source failure record is not an object")
            key = cell_tuple(raw)
            if key in seen:
                die("duplicate all-K source failure {} {}".format(arm, key))
            seen.add(key)
            cause = raw.get("cause")
            if cause not in ("field_shortfall", "q>H"):
                die("unknown all-K failure cause {}".format(cause))
            if cause == "field_shortfall" and (
                    raw.get("binary_def") != 12 or raw.get("heavy_gain") != 11):
                die("field-shortfall receipt changed for {} {}".format(arm, key))
            parsed.append(dict(raw))
        if int(value.get("failures", -1)) != len(parsed):
            die("all-K failure cardinality changed for {}".format(arm))
        records[arm] = parsed
    expected_field_counts = {"p48_r3": 7, "p32_r7": 9}
    for arm, count in expected_field_counts.items():
        if sum(row["cause"] == "field_shortfall" for row in records[arm]) != count:
            die("all-K finalist field-shortfall count changed for {}".format(arm))
    return records


def derive_hard_keys(records: Mapping[str, Sequence[dict]]) -> List[dict]:
    sources: Dict[Tuple[int, int, str], List[str]] = defaultdict(list)
    for arm in FINALIST_ARMS:
        for row in records[arm]:
            if row["cause"] == "field_shortfall":
                sources[cell_tuple(row)].append(arm)
    if len(sources) != 16:
        die("expected exact 16-cell finalist field-shortfall union")
    return [
        {**cell_record(key), "source_arms": sorted(sources[key])}
        for key in sorted(sources)
    ]


def common_failure_keys(records: Mapping[str, Sequence[dict]]) -> Set[Tuple[int, int, str]]:
    return {
        cell_tuple(row)
        for arm in SOURCE_ARMS
        for row in records[arm]
    }


def deterministic_controls(
    period: int,
    excluded: Set[Tuple[int, int, str]],
    count: int = CONTROL_COUNT,
    k_lo: int = K_LO,
    k_hi: int = K_HI,
) -> List[dict]:
    domain = "{}|common-success|P{}|".format(SCHEMA, period).encode("ascii")

    def scored() -> Iterator[Tuple[bytes, int, int, str]]:
        for K in range(k_lo, k_hi + 1):
            for seed_index in range(len(SEEDS)):
                for schedule in SCHEDULES:
                    key = (K, seed_index, schedule)
                    if key in excluded:
                        continue
                    suffix = "{}|{}|{}".format(*key).encode("ascii")
                    yield (hashlib.sha256(domain + suffix).digest(), *key)

    selected = heapq.nsmallest(count, scored())
    if len(selected) != count:
        die("insufficient common-success controls for P{}".format(period))
    keys = sorted((item[1], item[2], item[3]) for item in selected)
    if len(set(keys)) != count or any(key in excluded for key in keys):
        die("invalid deterministic control sample")
    return [cell_record(key) for key in keys]


def enumerate_masks(period: int, rows: int, suffix_mask: int) -> List[dict]:
    result: List[dict] = []
    for index, selected in enumerate(itertools.combinations(range(10), rows)):
        mask = sum(1 << row for row in selected)
        result.append({
            "period": period,
            "rows": rows,
            "mask_index": index,
            "mask": mask,
            "mask_hex": "0x{:03x}".format(mask),
            "selected_y": list(selected),
            "canonical_suffix": mask == suffix_mask,
        })
    if len(result) != MASK_COUNT or sum(row["canonical_suffix"] for row in result) != 1:
        die("mask enumeration/canonical suffix invariant failed")
    return result


def expected_preamble(task: Mapping[str, object]) -> Dict[str, str]:
    return {
        "trials": "1", "threads": "1", "loss": "0.5",
        "seed": SEEDS[int(task["seed_index"])], "completion": "mixed",
        "mixed_period": str(task["period"]), "mixed_gf256_rows": "10",
        "mixed_gf16_rows": "2", "mixed_geometry": "shared-x",
        "mixed_residue_skew": "0", "mixed_residue_schedule": "constant",
        "mixed_residue_hash_seed": "0x0", "mixed_residue_hash_keyed": "0",
        "mixed_independent_extension_residues": "0",
        "mixed_grouped_gf256_rows": str(task["rows"]),
        "mixed_grouped_gf256_row_mask": "0x{:x}".format(int(task["mask"])),
        "mixed_grouped_gf256_hash_seed": "0xb7e15162",
        "mixed_grouped_final_h_a_columns": "12",
        "mixed_residue_buckets_requested": "auto",
        "mixed_extension_residue_seed_xor": "0x4e",
        "source_hits_override": "0", "packet_peel_seed_xor": "0x0",
        "packet_peel_seed_table": "none", "binary_dense_rows_override": "0",
        "binary_dense_two_anchor": "1", "binary_dense_two_anchor_phase": "0",
        "gf256_heavy_rows_override": "0", "odd_packet_peel_seed_xor": "0x0",
        "packet_row_seed_multiplier": "0x1", "packet_row_seed_avalanche": "0",
        "seed_block_bytes_override": "0", "overhead_stream": "salted",
        "full_payload_solve": "0", "schedule": str(task["schedule"]),
    }


def parse_preamble(line: str) -> Dict[str, str]:
    prefix = "# precodefail: "
    if not line.startswith(prefix):
        die("missing precodefail preamble")
    result: Dict[str, str] = {}
    for token in line[len(prefix):].split():
        if "=" not in token:
            die("malformed preamble token")
        key, value = token.split("=", 1)
        if key in result:
            die("duplicate preamble token {}".format(key))
        result[key] = value
    return result


def parse_benchmark_csv(data: bytes, task: Mapping[str, object]) -> List[Dict[str, str]]:
    try:
        text = data.decode("ascii")
    except UnicodeDecodeError as exc:
        die("benchmark output is not ASCII: {}".format(exc))
    lines = text.splitlines(keepends=True)
    if not lines or not all(line.endswith("\n") for line in lines) or \
            any("\r" in line or line == "\n" for line in lines):
        die("benchmark output is empty or lacks final LF")
    preamble = parse_preamble(lines[0][:-1])
    expected = expected_preamble(task)
    if set(preamble) != set(expected) or preamble != expected:
        die("benchmark preamble differs from the sealed task")
    reader = csv.DictReader(lines[1:])
    if tuple(reader.fieldnames or ()) != CSV_FIELDS:
        die("benchmark CSV header mismatch")
    rows = list(reader)
    try:
        observed_Ks = [int(row["N"]) for row in rows]
    except (KeyError, TypeError, ValueError) as exc:
        die("benchmark CSV has an invalid K ledger: {}".format(exc))
    if observed_Ks != list(task["Ks"]):
        die("benchmark CSV K ledger mismatch")
    fixed = {
        "bb": "64", "heavy_family": "periodic", "mix_count": "2",
        "overhead": "0", "trials": "1", "active_packet_peel_seed_xor": "0x0",
    }
    for row in rows:
        if any(row.get(key) != value for key, value in fixed.items()):
            die("benchmark fixed CSV field mismatch")
        classify(row)
    return rows


def smoke_frozen_binary(binary: Path, frozen: Path) -> Dict[str, object]:
    smoke_dir = frozen / "smoke"
    smoke_dir.mkdir()
    cases = (
        (48, 3, 0x380), (48, 3, 0x049),
        (32, 7, 0x3F8), (32, 7, 0x1DD),
    )
    records: List[dict] = []
    for index, (period, rows, mask) in enumerate(cases):
        task = {
            "period": period, "rows": rows, "mask": mask,
            "seed_index": 0, "schedule": "adversarial", "Ks": [257],
        }
        argv = [
            str(binary), "precodefail", "--bb-list", "64", "--overhead", "0",
            "--trials", "1", "--threads", "1", "--loss", "0.50",
            "--completion", "mixed", "--mix-count", "2",
            "--packet-peel-seed-xor", "0", "--N", "257",
            "--seed", SEEDS[0], "--schedule", "adversarial",
            "--mixed-period", str(period), "--mixed-geometry", "shared-x",
            "--mixed-gf256-rows", "10", "--mixed-gf16-rows", "2",
            "--mixed-grouped-gf256-rows", str(rows),
            "--mixed-grouped-gf256-row-mask", "0x{:x}".format(mask),
            "--mixed-residue-buckets", "auto", "--binary-dense-two-anchor",
        ]
        stdout, stderr = run_captured(
            argv, cwd=frozen, timeout=60.0,
            label="bounded row-mask smoke P{} mask 0x{:x}".format(period, mask),
        )
        parsed = parse_benchmark_csv(stdout, task)
        if stderr:
            die("bounded row-mask smoke produced stderr")
        stem = "case{:02d}.P{}.r{}.mask{:03x}".format(index, period, rows, mask)
        stdout_path = smoke_dir / (stem + ".stdout")
        stderr_path = smoke_dir / (stem + ".stderr")
        stdout_path.write_bytes(stdout)
        stderr_path.write_bytes(stderr)
        records.append({
            "case": index, "period": period, "rows": rows,
            "mask": mask, "mask_hex": "0x{:03x}".format(mask), "K": 257,
            "argv": ["frozen/" + BINARY_NAME, *argv[1:]],
            "stdout": str(stdout_path.relative_to(frozen.parent)),
            "stdout_sha256": sha256_bytes(stdout),
            "stderr": str(stderr_path.relative_to(frozen.parent)),
            "stderr_sha256": sha256_bytes(stderr),
            "outcome": classify(parsed[0]),
        })
    record = {
        "schema": SCHEMA + ".binary_smoke", "bounded": True,
        "binary_sha256": sha256_file(binary), "cases": records,
    }
    write_once(frozen / "binary_smoke.json", canonical_json(record))
    return record


def build_tasks(
    binary: str,
    masks: Mapping[int, Sequence[dict]],
    hard: Sequence[dict],
    controls: Mapping[int, Sequence[dict]],
    control_count: int = CONTROL_COUNT,
) -> List[dict]:
    tasks: List[dict] = []
    common = [
        binary, "precodefail", "--bb-list", "64", "--overhead", "0",
        "--trials", "1", "--threads", "1", "--loss", "0.50",
        "--completion", "mixed", "--mix-count", "2",
        "--packet-peel-seed-xor", "0",
    ]
    # Seal the hard phase before the control phase.  The launcher also places
    # a process-pool barrier between these phases, so no control task can start
    # while a hard task is still running.
    for stage in ("hard", "control"):
        for panel in PANELS:
            period = int(panel["period"])
            source = hard if stage == "hard" else controls[period]
            for mask_entry in masks[period]:
                mask = int(mask_entry["mask"])
                strata: Dict[Tuple[int, str], List[int]] = defaultdict(list)
                for record in source:
                    key = cell_tuple(record)
                    strata[(key[1], key[2])].append(key[0])
                for seed_index in range(len(SEEDS)):
                    for schedule in SCHEDULES:
                        Ks = sorted(strata.get((seed_index, schedule), []))
                        for chunk_index, offset in enumerate(range(0, len(Ks), CHUNK_MAX)):
                            chunk = Ks[offset:offset + CHUNK_MAX]
                            job = len(tasks)
                            name = (
                                "{:05d}.P{}.r{}.mask{:03x}.{}.seed{}.{}.chunk{:02d}.csv"
                                .format(job, period, panel["rows"], mask, stage,
                                        seed_index, schedule, chunk_index)
                            )
                            task = {
                                "job": job, "panel": panel["name"],
                                "period": period, "rows": panel["rows"],
                                "mask_index": mask_entry["mask_index"], "mask": mask,
                                "mask_hex": mask_entry["mask_hex"],
                                "canonical_suffix": mask_entry["canonical_suffix"],
                                "stage": stage, "seed_index": seed_index,
                                "seed": SEEDS[seed_index], "schedule": schedule,
                                "chunk_index": chunk_index, "Ks": chunk,
                                "output_name": name,
                            }
                            task["argv"] = [
                                *common, "--N", ",".join(map(str, chunk)),
                                "--seed", SEEDS[seed_index], "--schedule", schedule,
                                "--mixed-period", str(period),
                                "--mixed-geometry", "shared-x",
                                "--mixed-gf256-rows", "10", "--mixed-gf16-rows", "2",
                                "--mixed-grouped-gf256-rows", str(panel["rows"]),
                                "--mixed-grouped-gf256-row-mask", "0x{:x}".format(mask),
                                "--mixed-residue-buckets", "auto",
                                "--binary-dense-two-anchor",
                            ]
                            tasks.append(task)
    expected_cells = len(PANELS) * MASK_COUNT * (len(hard) + control_count)
    if sum(len(task["Ks"]) for task in tasks) != expected_cells:
        die("task ledger does not cover the exact Cartesian product")
    if any(not task["Ks"] or len(task["Ks"]) > CHUNK_MAX for task in tasks):
        die("empty or oversized task generated")
    return tasks


def plan_from_summary(summary: Mapping[str, object], binary: str) -> Dict[str, object]:
    records = validate_source_summary(summary)
    hard = derive_hard_keys(records)
    failures = common_failure_keys(records)
    controls = {
        int(panel["period"]): deterministic_controls(int(panel["period"]), failures)
        for panel in PANELS
    }
    hard_set = {cell_tuple(record) for record in hard}
    if any(hard_set & {cell_tuple(record) for record in values}
           for values in controls.values()):
        die("hard/control ledgers overlap")
    masks = {
        int(panel["period"]): enumerate_masks(
            int(panel["period"]), int(panel["rows"]), int(panel["suffix_mask"]))
        for panel in PANELS
    }
    tasks = build_tasks(binary, masks, hard, controls)
    return {"hard": hard, "controls": controls, "masks": masks, "tasks": tasks}


def artifact_map(root: Path) -> Dict[str, Path]:
    result = {
        "controller": root / "frozen" / CONTROLLER_NAME,
        "binary": root / "frozen" / BINARY_NAME,
        "sampler": root / "frozen" / SAMPLER_NAME,
        "source_summary": root / "frozen" / SOURCE_SUMMARY_NAME,
        "build_receipt": root / "frozen" / "build_receipt.json",
        "binary_smoke": root / "frozen" / "binary_smoke.json",
        "design": root / "design.json",
        "hard_keys": root / "hard_keys.jsonl",
        "p48_controls": root / "controls_p48.jsonl",
        "p32_controls": root / "controls_p32.jsonl",
        "p48_masks": root / "masks_p48_r3.jsonl",
        "p32_masks": root / "masks_p32_r7.jsonl",
        "tasks": root / "tasks_manifest.jsonl",
    }
    for name in (
        "source_tree.txt", "CMakeCache.txt", "build.ninja",
        "compile_commands.json", "configure.stdout", "configure.stderr",
        "build.stdout", "build.stderr", "compiler.stdout", "compiler.stderr",
        "readelf-notes.stdout", "readelf-notes.stderr",
    ):
        result["build_" + name] = root / "frozen" / "build" / name
    for index, (period, rows, mask) in enumerate((
            (48, 3, 0x380), (48, 3, 0x049),
            (32, 7, 0x3F8), (32, 7, 0x1DD))):
        stem = "case{:02d}.P{}.r{}.mask{:03x}".format(index, period, rows, mask)
        result["smoke_{}_stdout".format(index)] = \
            root / "frozen" / "smoke" / (stem + ".stdout")
        result["smoke_{}_stderr".format(index)] = \
            root / "frozen" / "smoke" / (stem + ".stderr")
    return result


def command_prepare(args: argparse.Namespace) -> None:
    root = Path(args.output_root).resolve()
    source_path = Path(args.allk_root).resolve() / SOURCE_SUMMARY_NAME
    sampler = Path(args.sampler).resolve()
    controller = Path(__file__).resolve()
    repo = controller.parent.parent.resolve()
    if root.exists():
        die("output root already exists: {}".format(root))
    staging = root.parent / (root.name + ".prepare.{}".format(os.getpid()))
    if staging.exists():
        die("preparation staging path already exists: {}".format(staging))
    if not source_path.is_file() or source_path.is_symlink() or \
            sha256_file(source_path) != SOURCE_SUMMARY_SHA256:
        die("all-K validated-summary hash mismatch")
    for path, label in ((sampler, "sampler"), (controller, "controller")):
        if not path.is_file() or path.is_symlink():
            die("{} is not a regular file: {}".format(label, path))
    try:
        staging.mkdir(parents=True)
        for directory in (
                "frozen", "raw", "stderr", "exit", "thermal", "segments",
                "attempts", "job_receipts"):
            (staging / directory).mkdir()
        frozen = staging / "frozen"
        shutil.copyfile(str(controller), str(frozen / CONTROLLER_NAME))
        shutil.copyfile(str(sampler), str(frozen / SAMPLER_NAME))
        shutil.copyfile(str(source_path), str(frozen / SOURCE_SUMMARY_NAME))
        for name in (CONTROLLER_NAME, SAMPLER_NAME):
            os.chmod(str(frozen / name), 0o555)

        build = fresh_build_binary(repo, frozen, int(args.workers))
        smoke = smoke_frozen_binary(frozen / BINARY_NAME, frozen)
        if smoke.get("binary_sha256") != build.get("binary_sha256"):
            die("binary smoke is not bound to the fresh build")
        summary = load_json(frozen / SOURCE_SUMMARY_NAME)
        # The task argv is sealed to the final root, not the transient staging
        # directory that is atomically renamed after every receipt is written.
        frozen_binary = str(root / "frozen" / BINARY_NAME)
        plan = plan_from_summary(summary, frozen_binary)
        write_once(staging / "hard_keys.jsonl", json_lines(plan["hard"]))
        for period in (48, 32):
            write_once(staging / "controls_p{}.jsonl".format(period),
                       json_lines(plan["controls"][period]))
        write_once(staging / "masks_p48_r3.jsonl", json_lines(plan["masks"][48]))
        write_once(staging / "masks_p32_r7.jsonl", json_lines(plan["masks"][32]))
        write_once(staging / "tasks_manifest.jsonl", json_lines(plan["tasks"]))

        design = {
            "schema": SCHEMA + ".design", "root": str(root),
            "git_head": build["source"]["commit"],
            "git_tree_oid": build["source"]["tree_oid"],
            "build_receipt_sha256": sha256_file(frozen / "build_receipt.json"),
            "binary_sha256": build["binary_sha256"],
            "binary_smoke_sha256": sha256_file(frozen / "binary_smoke.json"),
            "codec_implementation_head": CODEC_IMPLEMENTATION_HEAD,
            "source_head": SOURCE_HEAD, "source_summary_sha256": SOURCE_SUMMARY_SHA256,
            "source_policy": "exact finalist field-shortfall union; controls common-success in all four source arms",
            "hard_cells": len(plan["hard"]), "controls_per_period": CONTROL_COUNT,
            "masks_per_period": MASK_COUNT, "periods": [48, 32], "rows": [3, 7],
            "canonical_suffix_masks": {"48": "0x380", "32": "0x3f8"},
            "seeds": list(SEEDS), "schedules": list(SCHEDULES),
            "K_range": [K_LO, K_HI], "bb": 64, "loss": 0.5,
            "overhead": 0, "trials": 1, "mix_count": 2,
            "seed_fixes": "none", "binary_dense_two_anchor": True,
            "cells_per_period": MASK_COUNT * (len(plan["hard"]) + CONTROL_COUNT),
            "total_cells": len(PANELS) * MASK_COUNT * (len(plan["hard"]) + CONTROL_COUNT),
            "task_count": len(plan["tasks"]), "chunk_max": CHUNK_MAX,
            "stage_order": ["hard", "control"],
            "worker_count": args.workers, "thermal_single_sampler_required": True,
            "thermal_cpu_busy_floor_pct": str(CPU_BUSY_FLOOR),
            "thermal_busy_floor_stages": list(CPU_BUSY_FLOOR_STAGES),
            "thermal_interval_seconds": str(THERMAL_INTERVAL_SECONDS),
            "thermal_min_gap_seconds": str(THERMAL_MIN_GAP_SECONDS),
            "thermal_max_gap_seconds": str(THERMAL_MAX_GAP_SECONDS),
            "thermal_ready_samples": THERMAL_READY_SAMPLES,
            "retry_policy": "stage-atomic non-selective retry after failed/interrupted segment",
            "timing_promotional": False,
            "timing_policy": "NONPROMOTIONAL saturated rank-proxy screen; recovery and work receipts only",
        }
        write_once(staging / "design.json", canonical_json(design))
        paths = artifact_map(staging)
        receipts = {
            "schema": SCHEMA + ".prelaunch", "artifacts": {
                key: {"path": str(path.relative_to(staging)), "sha256": sha256_file(path)}
                for key, path in sorted(paths.items())
            },
        }
        write_once(staging / "prelaunch_receipts.json", canonical_json(receipts))
        os.replace(str(staging), str(root))
    except BaseException:
        if staging.exists():
            shutil.rmtree(staging)
        raise
    print(json.dumps({
        "status": "PREPARED_NOT_LAUNCHED", "root": str(root),
        "prelaunch_receipts_sha256": sha256_file(root / "prelaunch_receipts.json"),
        "tasks": len(plan["tasks"]), "cells": design["total_cells"],
        "source_commit": design["git_head"], "source_tree": design["git_tree_oid"],
        "binary_sha256": design["binary_sha256"], "smoke_cases": len(smoke["cases"]),
    }, indent=2, sort_keys=True))


def read_jsonl(path: Path) -> List[dict]:
    result: List[dict] = []
    try:
        raw_lines = path.read_bytes().splitlines(keepends=True)
    except OSError as exc:
        die("cannot read {}: {}".format(path, exc))
    for index, raw in enumerate(raw_lines):
        if not raw.endswith(b"\n"):
            die("{} line {} lacks LF".format(path, index + 1))
        try:
            value = json.loads(raw.decode("ascii"))
        except (UnicodeError, json.JSONDecodeError) as exc:
            die("bad JSONL {} line {}: {}".format(path, index + 1, exc))
        if not isinstance(value, dict) or canonical_json(value) != raw:
            die("noncanonical JSONL {} line {}".format(path, index + 1))
        result.append(value)
    return result


def verify_root(
    root: Path,
    expected_receipts: str,
    *,
    runtime_tolerant: bool = False,
) -> Tuple[dict, List[dict]]:
    receipt_path = root / "prelaunch_receipts.json"
    if not re.fullmatch(r"[0-9a-f]{64}", expected_receipts or ""):
        die("expected receipt hash must be 64 lowercase hex digits")
    if sha256_file(receipt_path) != expected_receipts:
        die("external prelaunch receipt hash mismatch")
    receipts = load_json(receipt_path)
    if receipts.get("schema") != SCHEMA + ".prelaunch":
        die("prelaunch schema mismatch")
    artifacts = receipts.get("artifacts")
    if not isinstance(artifacts, dict) or set(artifacts) != set(artifact_map(root)):
        die("prelaunch artifact set mismatch")
    for key, expected_path in artifact_map(root).items():
        record = artifacts[key]
        if not isinstance(record, dict) or record.get("path") != str(expected_path.relative_to(root)):
            die("artifact path mismatch: {}".format(key))
        if expected_path.is_symlink() or not expected_path.is_file():
            die("artifact is not a regular file: {}".format(key))
        if sha256_file(expected_path) != record.get("sha256"):
            die("artifact hash mismatch: {}".format(key))
    design = load_json(root / "design.json")
    if design.get("schema") != SCHEMA + ".design" or design.get("root") != str(root):
        die("design identity mismatch")
    if design.get("codec_implementation_head") != CODEC_IMPLEMENTATION_HEAD or \
            design.get("source_summary_sha256") != SOURCE_SUMMARY_SHA256 or \
            design.get("thermal_cpu_busy_floor_pct") != str(CPU_BUSY_FLOOR) or \
            design.get("thermal_busy_floor_stages") != list(CPU_BUSY_FLOOR_STAGES) or \
            design.get("thermal_interval_seconds") != str(THERMAL_INTERVAL_SECONDS) or \
            design.get("thermal_min_gap_seconds") != str(THERMAL_MIN_GAP_SECONDS) or \
            design.get("thermal_max_gap_seconds") != str(THERMAL_MAX_GAP_SECONDS) or \
            design.get("thermal_ready_samples") != THERMAL_READY_SAMPLES or \
            design.get("retry_policy") != \
                "stage-atomic non-selective retry after failed/interrupted segment":
        die("design trust anchor mismatch")
    build = load_json(root / "frozen" / "build_receipt.json")
    smoke = load_json(root / "frozen" / "binary_smoke.json")
    source = build.get("source")
    if not isinstance(source, dict) or \
            build.get("schema") != SCHEMA + ".fresh_build" or \
            build.get("fresh_build") is not True or \
            build.get("binary_sha256") != sha256_file(root / "frozen" / BINARY_NAME) or \
            source.get("commit") != design.get("git_head") or \
            source.get("tree_oid") != design.get("git_tree_oid") or \
            design.get("binary_sha256") != build.get("binary_sha256") or \
            design.get("build_receipt_sha256") != sha256_file(
                root / "frozen" / "build_receipt.json"):
        die("fresh-build provenance binding mismatch")
    if re.fullmatch(r"[0-9a-f]{40}", str(source.get("commit", ""))) is None or \
            re.fullmatch(r"[0-9a-f]{40}", str(source.get("tree_oid", ""))) is None or \
            re.fullmatch(r"[0-9a-f]{64}", str(source.get("source_archive_sha256", ""))) is None or \
            source.get("status_porcelain_sha256") != sha256_bytes(b"") or \
            source.get("tree_manifest_sha256") != sha256_file(
                root / "frozen" / "build" / "source_tree.txt") or \
            source.get("tree_manifest_bytes") != (
                root / "frozen" / "build" / "source_tree.txt").stat().st_size:
        die("fresh-build clean source/tree receipt mismatch")
    if smoke.get("schema") != SCHEMA + ".binary_smoke" or \
            smoke.get("bounded") is not True or \
            smoke.get("binary_sha256") != build.get("binary_sha256") or \
            not isinstance(smoke.get("cases"), list) or len(smoke["cases"]) != 4 or \
            design.get("binary_smoke_sha256") != sha256_file(
                root / "frozen" / "binary_smoke.json"):
        die("bounded binary-smoke provenance binding mismatch")
    expected_build_names = {
        path.name for key, path in artifact_map(root).items()
        if key.startswith("build_") and key != "build_receipt"
    }
    expected_smoke_names = {
        path.name for key, path in artifact_map(root).items()
        if key.startswith("smoke_")
    }
    for directory, expected_names in (
            (root / "frozen" / "build", expected_build_names),
            (root / "frozen" / "smoke", expected_smoke_names)):
        entries = list(directory.iterdir())
        actual_names = {path.name for path in entries}
        if actual_names != expected_names or any(
                path.is_symlink() or not path.is_file() for path in entries):
            die("frozen evidence inventory mismatch: {}".format(directory.name))
    evidence = build.get("evidence")
    if not isinstance(evidence, dict) or set(evidence) != expected_build_names:
        die("fresh-build evidence receipt set mismatch")
    for name in expected_build_names:
        path = root / "frozen" / "build" / name
        record = evidence[name]
        if not isinstance(record, dict) or record.get("sha256") != sha256_file(path) or \
                record.get("bytes") != path.stat().st_size:
            die("fresh-build evidence receipt mismatch: {}".format(name))
    expected_cases = ((48, 3, 0x380), (48, 3, 0x049),
                      (32, 7, 0x3F8), (32, 7, 0x1DD))
    for index, ((period, rows, mask), record) in enumerate(zip(expected_cases, smoke["cases"])):
        stem = "case{:02d}.P{}.r{}.mask{:03x}".format(index, period, rows, mask)
        stdout_path = root / "frozen" / "smoke" / (stem + ".stdout")
        stderr_path = root / "frozen" / "smoke" / (stem + ".stderr")
        task = {"period": period, "rows": rows, "mask": mask,
                "seed_index": 0, "schedule": "adversarial", "Ks": [257]}
        if not isinstance(record, dict) or record.get("case") != index or \
                record.get("period") != period or record.get("rows") != rows or \
                record.get("mask") != mask or record.get("K") != 257 or \
                record.get("stdout") != str(stdout_path.relative_to(root)) or \
                record.get("stderr") != str(stderr_path.relative_to(root)) or \
                record.get("stdout_sha256") != sha256_file(stdout_path) or \
                record.get("stderr_sha256") != sha256_file(stderr_path) or \
                stderr_path.stat().st_size:
            die("bounded row-mask smoke receipt mismatch")
        parsed = parse_benchmark_csv(stdout_path.read_bytes(), task)
        if record.get("outcome") != classify(parsed[0]):
            die("bounded row-mask smoke outcome receipt mismatch")
    summary = load_json(root / "frozen" / SOURCE_SUMMARY_NAME)
    plan = plan_from_summary(summary, str(root / "frozen" / BINARY_NAME))
    expected_files = {
        root / "hard_keys.jsonl": json_lines(plan["hard"]),
        root / "controls_p48.jsonl": json_lines(plan["controls"][48]),
        root / "controls_p32.jsonl": json_lines(plan["controls"][32]),
        root / "masks_p48_r3.jsonl": json_lines(plan["masks"][48]),
        root / "masks_p32_r7.jsonl": json_lines(plan["masks"][32]),
        root / "tasks_manifest.jsonl": json_lines(plan["tasks"]),
    }
    for path, expected in expected_files.items():
        if path.read_bytes() != expected:
            die("regenerated ledger mismatch: {}".format(path.name))
    tasks = plan["tasks"]
    if design.get("task_count") != len(tasks) or \
            design.get("total_cells") != sum(len(task["Ks"]) for task in tasks):
        die("design Cartesian receipt mismatch")
    expected_outputs = {task["output_name"] for task in tasks}
    for directory, suffix in (("raw", ""), ("stderr", ".stderr"), ("exit", ".exit")):
        entries = list((root / directory).iterdir())
        if any(path.is_symlink() or not path.is_file() for path in entries):
            die("invalid {} runtime artifact".format(directory))
        actual = {path.name for path in entries}
        expected = {name + suffix for name in expected_outputs}
        partial = {name for name in actual if ".part." in name}
        if not (actual - partial) <= expected:
            die("unexpected {} output(s)".format(directory))
        if partial and not runtime_tolerant:
            die("partial {} output remains".format(directory))
        if any(not any(name.startswith(prefix + ".part.") for prefix in expected)
               for name in partial):
            die("unexpected partial {} output(s)".format(directory))
    return design, tasks


def job_receipt_path(root: Path, job: int) -> Path:
    return root / "job_receipts" / "job{:05d}.json".format(job)


def completed_jobs(root: Path, tasks: Sequence[dict]) -> Set[int]:
    completed: Set[int] = set()
    for task in tasks:
        name = task["output_name"]
        paths = (
            root / "raw" / name,
            root / "stderr" / (name + ".stderr"),
            root / "exit" / (name + ".exit"),
            job_receipt_path(root, int(task["job"])),
        )
        present = tuple(path.exists() for path in paths)
        if any(present) and not all(present):
            die("partial output triplet for job {}".format(task["job"]))
        if all(present):
            if any(path.is_symlink() or not path.is_file() for path in paths):
                die("completed job contains an indirect artifact")
            receipt = load_json(paths[3])
            if not paths[0].stat().st_size or paths[1].stat().st_size or \
                    paths[2].read_text(encoding="ascii") != "0\n" or \
                    receipt.get("schema") != SCHEMA + ".job" or \
                    receipt.get("job") != task["job"] or \
                    receipt.get("stage") != task["stage"] or \
                    receipt.get("output_name") != name or \
                    receipt.get("returncode") != 0 or \
                    receipt.get("stdout_sha256") != sha256_file(paths[0]) or \
                    receipt.get("stderr_sha256") != sha256_file(paths[1]) or \
                    receipt.get("exit_sha256") != sha256_file(paths[2]) or \
                    not isinstance(receipt.get("segment"), int) or \
                    isinstance(receipt.get("segment"), bool) or \
                    not isinstance(receipt.get("started_monotonic_s"), (int, float)) or \
                    isinstance(receipt.get("started_monotonic_s"), bool) or \
                    not isinstance(receipt.get("ended_monotonic_s"), (int, float)) or \
                    isinstance(receipt.get("ended_monotonic_s"), bool) or \
                    receipt["ended_monotonic_s"] < receipt["started_monotonic_s"]:
                die("unclean completed job {}".format(task["job"]))
            completed.add(int(task["job"]))
    return completed


def process_state(pid: int) -> Optional[str]:
    path = Path("/proc") / str(pid) / "stat"
    try:
        raw = path.read_bytes()
    except PermissionError:
        result = subprocess.run(
            ["sudo", "-n", "cat", str(path)], stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
        )
        if result.returncode:
            return None
        raw = result.stdout
    except OSError:
        return None
    try:
        return raw.rsplit(b") ", 1)[1].split(None, 1)[0].decode("ascii")
    except (IndexError, UnicodeError):
        return None


def process_start_ticks(pid: int) -> Optional[int]:
    path = Path("/proc") / str(pid) / "stat"
    try:
        raw = path.read_bytes()
    except PermissionError:
        result = subprocess.run(
            ["sudo", "-n", "cat", str(path)], stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
        )
        if result.returncode:
            return None
        raw = result.stdout
    except OSError:
        return None
    try:
        # Field 22 is the process start time; after removing PID and comm it is
        # zero-based field 19 in the remainder beginning with the state.
        return int(raw.rsplit(b") ", 1)[1].split()[19])
    except (IndexError, ValueError):
        return None


def process_alive(pid: int) -> bool:
    if pid <= 1:
        return False
    try:
        os.kill(pid, 0)
        return process_state(pid) != "Z"
    except ProcessLookupError:
        return False
    except PermissionError:
        # Root-owned samplers return EPERM to a non-root kill(2).  A plain
        # exception must never be interpreted as process death.
        result = subprocess.run(
            ["sudo", "-n", "kill", "-0", str(pid)],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
        if result.returncode == 0:
            return process_state(pid) != "Z"
        if (Path("/proc") / str(pid)).exists():
            die("cannot determine liveness of permission-protected PID {}".format(pid))
        return False


def process_tokens(pid: int) -> Optional[List[str]]:
    path = Path("/proc") / str(pid) / "cmdline"
    try:
        raw = path.read_bytes()
    except PermissionError:
        result = subprocess.run(
            ["sudo", "-n", "cat", str(path)], stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
        )
        if result.returncode:
            return None
        raw = result.stdout
    except OSError:
        return None
    return [token.decode("utf-8", "replace") for token in raw.split(b"\0") if token]


def other_samplers(excluded_pids: Iterable[int] = ()) -> List[dict]:
    excluded = {os.getpid(), *excluded_pids}
    found: List[dict] = []
    for proc in Path("/proc").iterdir():
        if not proc.name.isdigit() or int(proc.name) in excluded:
            continue
        pid = int(proc.name)
        tokens = process_tokens(pid)
        if tokens is not None and any(
                Path(token).name == SAMPLER_NAME for token in tokens) and \
                process_alive(pid):
            found.append({"pid": pid, "command": " ".join(tokens)})
    return sorted(found, key=lambda item: item["pid"])


def atomic_result(path: Path, data: bytes) -> None:
    if path.exists() or path.is_symlink():
        die("refusing to replace runtime artifact {}".format(path))
    temporary = Path(str(path) + ".part.{}.{}".format(
        os.getpid(), threading.get_ident()))
    with temporary.open("xb") as stream:
        stream.write(data)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(str(temporary), str(path))


def segment_path(root: Path, segment: int, kind: str) -> Path:
    return root / "segments" / "segment{:03d}.{}.json".format(segment, kind)


def segment_indices(root: Path) -> List[int]:
    indices: Set[int] = set()
    patterns = (
        (root / "segments", re.compile(r"segment([0-9]{3})\.(?:intent|ready|final)\.json$")),
        (root / "thermal", re.compile(r"segment([0-9]{3})\.(?:csv|pid|stderr)$")),
        (root / "attempts", re.compile(r"segment([0-9]{3})$")),
    )
    for directory, pattern in patterns:
        for path in directory.iterdir():
            match = pattern.fullmatch(path.name)
            if match:
                indices.add(int(match.group(1)))
    ordered = sorted(indices)
    if ordered and ordered != list(range(ordered[-1] + 1)):
        die("campaign segment numbering has a gap")
    return ordered


def segment_jobs(intent: Mapping[str, object], tasks: Sequence[dict]) -> List[dict]:
    raw = intent.get("jobs")
    if not isinstance(raw, list) or not raw or any(
            not isinstance(job, int) or isinstance(job, bool) for job in raw):
        die("segment intent has an invalid job ledger")
    if raw != sorted(set(raw)):
        die("segment intent job ledger is not canonical")
    task_by_job = {int(task["job"]): task for task in tasks}
    if any(job not in task_by_job for job in raw):
        die("segment intent references an unknown job")
    selected = [task_by_job[job] for job in raw]
    stage = intent.get("stage")
    if stage not in ("hard", "control") or any(task["stage"] != stage for task in selected):
        die("segment intent crosses campaign stages")
    if intent.get("jobs_sha256") != sha256_bytes(json_lines(raw)):
        die("segment intent job hash mismatch")
    return selected


def terminate_verified_process(pid: int, expected: Sequence[str]) -> str:
    if not process_alive(pid):
        return "already_exited"
    tokens = process_tokens(pid)
    if tokens != list(expected):
        die("refusing to terminate PID {} with changed command identity".format(pid))
    subprocess.run(
        ["sudo", "-n", "kill", "-TERM", str(pid)],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False,
    )
    deadline = time.monotonic() + 10.0
    while time.monotonic() < deadline and process_alive(pid):
        time.sleep(0.05)
    if process_alive(pid):
        subprocess.run(
            ["sudo", "-n", "kill", "-KILL", str(pid)],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False,
        )
        deadline = time.monotonic() + 5.0
        while time.monotonic() < deadline and process_alive(pid):
            time.sleep(0.05)
        if process_alive(pid):
            die("verified owned process {} did not stop".format(pid))
        return "killed"
    return "terminated"


def terminate_owned_segment_processes(root: Path, segment: int, tasks: Sequence[dict]) -> List[dict]:
    actions: List[dict] = []
    thermal_pid = root / "thermal" / "segment{:03d}.pid".format(segment)
    thermal_csv = root / "thermal" / "segment{:03d}.csv".format(segment)
    if thermal_pid.is_symlink():
        die("owned thermal PID file is a symlink")
    if thermal_pid.exists():
        try:
            pid = int(thermal_pid.read_text(encoding="ascii").strip())
        except (OSError, ValueError):
            die("malformed owned thermal PID file for segment {}".format(segment))
        ready_path = segment_path(root, segment, "ready")
        expected_start: Optional[int] = None
        if ready_path.is_file() and not ready_path.is_symlink():
            ready = load_json(ready_path)
            if ready.get("sampler_pid") == pid and \
                    isinstance(ready.get("sampler_start_ticks"), int) and \
                    not isinstance(ready.get("sampler_start_ticks"), bool):
                expected_start = int(ready["sampler_start_ticks"])
        if process_alive(pid) and expected_start is not None and \
                process_start_ticks(pid) != expected_start:
            actions.append({"kind": "thermal", "pid": pid,
                            "action": "stale_reused_pid"})
            try:
                thermal_pid.unlink()
            except FileNotFoundError:
                pass
            pid = -1
        tokens = process_tokens(pid) if process_alive(pid) else None
        if process_alive(pid) and tokens is None:
            die("cannot inspect the live owned thermal sampler")
        if tokens is not None:
            frozen_sampler = str(root / "frozen" / SAMPLER_NAME)
            if frozen_sampler not in tokens or str(thermal_csv) not in tokens:
                actions.append({"kind": "thermal", "pid": pid,
                                "action": "stale_changed_command"})
                try:
                    thermal_pid.unlink()
                except FileNotFoundError:
                    pass
                tokens = None
            actions.append({"kind": "thermal", "pid": pid,
                            "action": terminate_verified_process(pid, tokens)})
        if not process_alive(pid):
            try:
                thermal_pid.unlink()
            except FileNotFoundError:
                pass
    task_by_job = {int(task["job"]): task for task in tasks}
    attempt_dir = root / "attempts" / "segment{:03d}".format(segment)
    if attempt_dir.is_dir():
        for pid_path in sorted(attempt_dir.glob("job*.pid")):
            if pid_path.is_symlink() or not pid_path.is_file():
                die("benchmark PID receipt is indirect")
            match = re.fullmatch(r"job([0-9]{5})\.pid", pid_path.name)
            if not match:
                die("malformed benchmark PID receipt")
            job = int(match.group(1))
            if job not in task_by_job:
                die("benchmark PID receipt names an unknown job")
            identity = load_json(pid_path)
            pid = identity.get("pid")
            start_ticks = identity.get("start_ticks")
            if not isinstance(pid, int) or isinstance(pid, bool) or pid <= 1 or \
                    not isinstance(start_ticks, int) or isinstance(start_ticks, bool) or \
                    start_ticks < 0:
                die("malformed benchmark PID identity")
            if not process_alive(pid):
                continue
            if process_start_ticks(pid) != start_ticks:
                actions.append({"kind": "benchmark", "job": job, "pid": pid,
                                "action": "stale_reused_pid"})
                continue
            tokens = process_tokens(pid)
            if tokens is None:
                die("cannot inspect a live owned benchmark PID")
            if tokens != task_by_job[job]["argv"]:
                actions.append({"kind": "benchmark", "job": job, "pid": pid,
                                "action": "stale_changed_command"})
                continue
            actions.append({"kind": "benchmark", "job": job, "pid": pid,
                            "action": terminate_verified_process(pid, tokens)})
    return actions


def rollback_segment_outputs(root: Path, selected: Sequence[dict]) -> List[int]:
    rolled_back: List[int] = []
    for task in selected:
        job = int(task["job"])
        name = task["output_name"]
        paths = (
            root / "raw" / name,
            root / "stderr" / (name + ".stderr"),
            root / "exit" / (name + ".exit"),
            job_receipt_path(root, job),
        )
        if any(path.exists() or path.is_symlink() for path in paths):
            rolled_back.append(job)
        for path in paths:
            try:
                path.unlink()
            except FileNotFoundError:
                pass
        for directory, prefix in (
                (root / "raw", name), (root / "stderr", name + ".stderr"),
                (root / "exit", name + ".exit")):
            for partial in directory.glob(prefix + ".part.*"):
                partial.unlink()
    return rolled_back


def seal_attempt_manifest(root: Path, segment: int) -> Tuple[str, int]:
    directory = root / "attempts" / "segment{:03d}".format(segment)
    if not directory.exists():
        directory.mkdir()
    if directory.is_symlink() or not directory.is_dir():
        die("attempt evidence directory is indirect")
    for partial in directory.glob("*.part.*"):
        partial.unlink()
    files = sorted(path for path in directory.iterdir() if path.is_file())
    if any(path.is_symlink() for path in files):
        die("attempt evidence contains a symlink")
    data = b"".join(
        "{}  {}\n".format(sha256_file(path), path.relative_to(root)).encode("ascii")
        for path in files
    )
    manifest = root / "segments" / "segment{:03d}.attempts.sha256".format(segment)
    write_once(manifest, data)
    return sha256_bytes(data), len(files)


def write_segment_final(
    root: Path,
    segment: int,
    intent: Mapping[str, object],
    state: str,
    *,
    failures: Sequence[dict],
    rolled_back: Sequence[int],
    process_actions: Sequence[dict],
    jobs_ended_monotonic_s: Optional[float],
    sampler_returncode: Optional[int],
) -> Dict[str, object]:
    thermal_csv = root / "thermal" / "segment{:03d}.csv".format(segment)
    thermal_stderr = root / "thermal" / "segment{:03d}.stderr".format(segment)
    attempts_sha256, attempt_file_count = seal_attempt_manifest(root, segment)
    record = {
        "schema": SCHEMA + ".segment_final", "segment": segment,
        "state": state, "stage": intent.get("stage"),
        "intent_sha256": sha256_file(segment_path(root, segment, "intent")),
        "ready_sha256": (
            sha256_file(segment_path(root, segment, "ready"))
            if segment_path(root, segment, "ready").is_file() else None),
        "ended_utc": datetime.now(timezone.utc).isoformat(),
        "jobs_ended_monotonic_s": jobs_ended_monotonic_s,
        "jobs": intent.get("jobs"), "jobs_sha256": intent.get("jobs_sha256"),
        "failures": list(failures), "rolled_back_jobs": sorted(rolled_back),
        "published_jobs": list(intent.get("jobs", [])) if state == "success" else [],
        "process_actions": list(process_actions),
        "sampler_returncode": sampler_returncode,
        "thermal_csv": thermal_csv.name,
        "thermal_csv_sha256": sha256_file(thermal_csv) if thermal_csv.is_file() else None,
        "thermal_stderr": thermal_stderr.name,
        "thermal_stderr_sha256": (
            sha256_file(thermal_stderr) if thermal_stderr.is_file() else None),
        "attempt_manifest": "segment{:03d}.attempts.sha256".format(segment),
        "attempt_manifest_sha256": attempts_sha256,
        "attempt_file_count": attempt_file_count,
        "retry_policy": "entire stage segment; no successful survivor is retained on failure",
    }
    write_once(segment_path(root, segment, "final"), canonical_json(record))
    return record


def reconcile_incomplete_segments(root: Path, tasks: Sequence[dict]) -> List[dict]:
    reconciled: List[dict] = []
    for segment in segment_indices(root):
        intent_path = segment_path(root, segment, "intent")
        final_path = segment_path(root, segment, "final")
        if final_path.exists():
            if not intent_path.exists():
                die("finalized segment lacks its intent")
            continue
        if not intent_path.exists():
            die("campaign segment artifacts lack an immutable intent")
        intent = load_json(intent_path)
        if intent.get("schema") != SCHEMA + ".segment_intent" or \
                intent.get("segment") != segment:
            die("incomplete segment intent is malformed")
        selected = segment_jobs(intent, tasks)
        actions = terminate_owned_segment_processes(root, segment, selected)
        rolled_back = rollback_segment_outputs(root, selected)
        final = write_segment_final(
            root, segment, intent, "interrupted", failures=[],
            rolled_back=rolled_back, process_actions=actions,
            jobs_ended_monotonic_s=None, sampler_returncode=None,
        )
        reconciled.append(final)
    return reconciled


def read_thermal_rows(path: Path) -> List[Dict[str, str]]:
    try:
        with path.open(newline="", encoding="ascii") as stream:
            reader = csv.DictReader(stream)
            if tuple(reader.fieldnames or ()) != THERMAL_FIELDS:
                die("thermal header mismatch: {}".format(path))
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as exc:
        die("cannot parse thermal stream {}: {}".format(path, exc))
    if any(set(row) != set(THERMAL_FIELDS) or any(value is None for value in row.values())
           for row in rows):
        die("thermal stream contains a malformed row")
    return rows


def sampler_identity(pid_file: Path, sampler_path: Path, csv_path: Path) -> Optional[int]:
    if pid_file.is_symlink():
        die("thermal sampler PID file is a symlink")
    if not pid_file.is_file():
        return None
    try:
        pid = int(pid_file.read_text(encoding="ascii").strip())
    except (OSError, ValueError):
        die("thermal sampler PID file is malformed")
    if not process_alive(pid):
        return None
    tokens = process_tokens(pid)
    if tokens is None or str(sampler_path) not in tokens or str(csv_path) not in tokens:
        die("thermal sampler PID is not bound to this segment")
    return pid


def wait_sampler_ready(
    sampler: subprocess.Popen,
    pid_file: Path,
    sampler_path: Path,
    csv_path: Path,
) -> Tuple[int, List[Dict[str, str]]]:
    deadline = time.monotonic() + THERMAL_READY_TIMEOUT_SECONDS
    while time.monotonic() < deadline:
        if sampler.poll() is not None:
            die("thermal sampler exited before readiness")
        pid = sampler_identity(pid_file, sampler_path, csv_path)
        if pid is not None and csv_path.is_file():
            rows = read_thermal_rows(csv_path)
            if len(rows) >= THERMAL_READY_SAMPLES:
                last = decimal_field(rows[-1], "monotonic_s")
                age = Decimal(str(time.monotonic())) - last
                if 0 <= age <= THERMAL_MAX_GAP_SECONDS:
                    return pid, rows
        time.sleep(0.05)
    die("thermal sampler did not provide two live data rows")


def wait_thermal_coverage(
    sampler: subprocess.Popen,
    csv_path: Path,
    end_monotonic_s: float,
) -> bool:
    deadline = time.monotonic() + THERMAL_END_TIMEOUT_SECONDS
    target = Decimal(str(end_monotonic_s))
    while time.monotonic() < deadline:
        if sampler.poll() is not None:
            return False
        rows = read_thermal_rows(csv_path)
        if rows and decimal_field(rows[-1], "monotonic_s") >= target:
            return True
        time.sleep(0.05)
    return False


def stop_launched_sampler(
    sampler: subprocess.Popen,
    pid_file: Path,
    sampler_path: Path,
    csv_path: Path,
) -> Tuple[Optional[int], List[dict]]:
    actions: List[dict] = []
    pid = sampler_identity(pid_file, sampler_path, csv_path)
    if pid is None and sampler.poll() is None:
        tokens = process_tokens(sampler.pid)
        if tokens is None or str(sampler_path) not in tokens or str(csv_path) not in tokens:
            die("cannot prove ownership of launched thermal sampler")
        pid = sampler.pid
    if pid is not None:
        tokens = process_tokens(pid) if process_alive(pid) else None
        if tokens is not None:
            actions.append({"kind": "thermal", "pid": pid,
                            "action": terminate_verified_process(pid, tokens)})
    try:
        returncode = sampler.wait(timeout=15.0)
    except subprocess.TimeoutExpired:
        tokens = process_tokens(sampler.pid)
        if tokens is None:
            die("launched sampler wrapper did not exit and cannot be inspected")
        actions.append({"kind": "thermal_wrapper", "pid": sampler.pid,
                        "action": terminate_verified_process(sampler.pid, tokens)})
        returncode = sampler.wait(timeout=5.0)
    try:
        pid_file.unlink()
    except FileNotFoundError:
        pass
    return returncode, actions


def run_segment(
    root: Path,
    design: Mapping[str, object],
    tasks: Sequence[dict],
    stage: str,
    previously_complete: int,
    resume: bool,
) -> Dict[str, object]:
    indices = segment_indices(root)
    segment = indices[-1] + 1 if indices else 0
    jobs = sorted(int(task["job"]) for task in tasks)
    intent = {
        "schema": SCHEMA + ".segment_intent", "segment": segment,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "created_monotonic_s": time.monotonic(),
        "stage": stage, "jobs": jobs,
        "jobs_sha256": sha256_bytes(json_lines(jobs)),
        "workers": design["worker_count"], "previously_complete": previously_complete,
        "resume": bool(resume),
        "retry_policy": "stage-atomic non-selective retry",
    }
    write_once(segment_path(root, segment, "intent"), canonical_json(intent))
    attempt_dir = root / "attempts" / "segment{:03d}".format(segment)
    attempt_dir.mkdir()
    thermal_csv = root / "thermal" / "segment{:03d}.csv".format(segment)
    thermal_pid = root / "thermal" / "segment{:03d}.pid".format(segment)
    thermal_stderr = root / "thermal" / "segment{:03d}.stderr".format(segment)
    sampler_path = root / "frozen" / SAMPLER_NAME
    sampler_error = thermal_stderr.open("xb")
    sampler_command = [
        "sudo", "-n", sys.executable, str(sampler_path),
        "--csv", str(thermal_csv), "--pid-file", str(thermal_pid),
        "--interval", str(THERMAL_INTERVAL_SECONDS),
        "--dimm-attempts", "5", "--dimm-retry-delay", "0.01",
    ]
    try:
        sampler = subprocess.Popen(
            sampler_command, stdout=subprocess.DEVNULL, stderr=sampler_error,
        )
    except OSError as exc:
        sampler_error.close()
        return write_segment_final(
            root, segment, intent, "failed",
            failures=[{"status": "sampler_spawn_error", "detail": str(exc)}],
            rolled_back=[], process_actions=[], jobs_ended_monotonic_s=None,
            sampler_returncode=None,
        )
    abort = threading.Event()
    active: Dict[int, subprocess.Popen] = {}
    active_lock = threading.Lock()
    failures: List[dict] = []
    process_actions: List[dict] = []
    jobs_ended: Optional[float] = None
    sampler_returncode: Optional[int] = None
    ready_written = False

    def terminate_active() -> None:
        with active_lock:
            snapshot = list(active.items())
        for job, process in snapshot:
            if process.poll() is not None:
                continue
            tokens = process_tokens(process.pid)
            if tokens is None:
                try:
                    process.terminate()
                    action = "terminated_via_owned_popen"
                except ProcessLookupError:
                    action = "already_exited"
                process_actions.append({
                    "kind": "benchmark", "job": job, "pid": process.pid,
                    "action": action,
                })
                continue
            process_actions.append({
                "kind": "benchmark", "job": job, "pid": process.pid,
                "action": terminate_verified_process(process.pid, tokens),
            })

    def run_one(task: dict) -> dict:
        job = int(task["job"])
        if abort.is_set():
            return {"job": job, "status": "cancelled_before_start"}
        started_monotonic = time.monotonic()
        try:
            process = subprocess.Popen(
                task["argv"], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                start_new_session=True,
            )
        except OSError as exc:
            return {"job": job, "status": "spawn_error", "detail": str(exc)}
        start_ticks = process_start_ticks(process.pid)
        if start_ticks is None:
            process.terminate()
            process.wait()
            return {"job": job, "status": "pid_identity_error"}
        atomic_result(
            attempt_dir / "job{:05d}.pid".format(job),
            canonical_json({"pid": process.pid, "start_ticks": start_ticks}),
        )
        with active_lock:
            active[job] = process
        if abort.is_set() and process.poll() is None:
            tokens = process_tokens(process.pid)
            if tokens is not None:
                terminate_verified_process(process.pid, tokens)
            else:
                process.terminate()
        stdout, stderr = process.communicate()
        with active_lock:
            active.pop(job, None)
        ended_monotonic = time.monotonic()
        exit_bytes = (str(process.returncode) + "\n").encode("ascii")
        for suffix, data in (("stdout", stdout), ("stderr", stderr), ("exit", exit_bytes)):
            atomic_result(attempt_dir / "job{:05d}.{}".format(job, suffix), data)
        if process.returncode != 0:
            return {"job": job, "status": "exit", "returncode": process.returncode}
        try:
            parse_benchmark_csv(stdout, task)
            if stderr:
                die("successful benchmark wrote stderr")
        except CampaignError as exc:
            return {"job": job, "status": "validation_error", "detail": str(exc)}
        if abort.is_set():
            return {"job": job, "status": "cancelled_after_run"}
        name = task["output_name"]
        raw_path = root / "raw" / name
        stderr_path = root / "stderr" / (name + ".stderr")
        exit_path = root / "exit" / (name + ".exit")
        atomic_result(raw_path, stdout)
        atomic_result(stderr_path, stderr)
        atomic_result(exit_path, exit_bytes)
        receipt = {
            "schema": SCHEMA + ".job", "job": job, "segment": segment,
            "stage": stage, "output_name": name, "returncode": 0,
            "started_monotonic_s": started_monotonic,
            "ended_monotonic_s": ended_monotonic,
            "stdout_sha256": sha256_bytes(stdout),
            "stderr_sha256": sha256_bytes(stderr),
            "exit_sha256": sha256_bytes(exit_bytes),
        }
        write_once(job_receipt_path(root, job), canonical_json(receipt))
        return {"job": job, "status": "success"}

    state = "failed"
    caught: Optional[BaseException] = None
    try:
        sampled_pid, ready_rows = wait_sampler_ready(
            sampler, thermal_pid, sampler_path, thermal_csv,
        )
        sampled_start_ticks = process_start_ticks(sampled_pid)
        if sampled_start_ticks is None:
            die("cannot bind thermal sampler to its process start time")
        jobs_started = time.monotonic()
        ready = {
            "schema": SCHEMA + ".segment_ready", "segment": segment,
            "intent_sha256": sha256_file(segment_path(root, segment, "intent")),
            "ready_utc": datetime.now(timezone.utc).isoformat(),
            "jobs_started_monotonic_s": jobs_started,
            "sampler_pid": sampled_pid, "samples_at_ready": len(ready_rows),
            "sampler_start_ticks": sampled_start_ticks,
            "last_sample_monotonic_s": ready_rows[-1]["monotonic_s"],
        }
        write_once(segment_path(root, segment, "ready"), canonical_json(ready))
        ready_written = True
        future_results: List[dict] = []
        with concurrent.futures.ThreadPoolExecutor(
                max_workers=int(design["worker_count"])) as pool:
            pending = {pool.submit(run_one, task) for task in tasks}
            while pending:
                if sampler.poll() is not None and not abort.is_set():
                    failures.append({"status": "sampler_exited", "returncode": sampler.returncode})
                    abort.set()
                    terminate_active()
                done, pending = concurrent.futures.wait(
                    pending, timeout=0.2,
                    return_when=concurrent.futures.FIRST_COMPLETED,
                )
                for future in done:
                    if future.cancelled():
                        result = {"status": "future_cancelled"}
                    else:
                        try:
                            result = future.result()
                        except BaseException as exc:
                            result = {"status": "worker_exception", "detail": repr(exc)}
                    future_results.append(result)
                    if result.get("status") != "success" and not abort.is_set():
                        failures.append(result)
                        abort.set()
                        terminate_active()
                if abort.is_set():
                    terminate_active()
                    for future in pending:
                        future.cancel()
            for result in future_results:
                if result.get("status") != "success" and result not in failures:
                    failures.append(result)
        jobs_ended = time.monotonic()
        if not failures and not wait_thermal_coverage(sampler, thermal_csv, jobs_ended):
            failures.append({"status": "thermal_end_coverage_failed"})
        if not failures:
            state = "success"
    except BaseException as exc:
        caught = exc
        abort.set()
        terminate_active()
        failures.append({"status": "launcher_exception", "detail": repr(exc)})
        jobs_ended = time.monotonic() if ready_written else None
    finally:
        try:
            sampler_returncode, sampler_actions = stop_launched_sampler(
                sampler, thermal_pid, sampler_path, thermal_csv,
            )
            process_actions.extend(sampler_actions)
        except BaseException as exc:
            if caught is None:
                caught = exc
            failures.append({"status": "sampler_stop_error", "detail": repr(exc)})
            state = "failed"
        sampler_error.close()
    if sampler_returncode != 0:
        failures.append({"status": "sampler_returncode", "returncode": sampler_returncode})
        state = "failed"
    if not thermal_stderr.is_file() or thermal_stderr.stat().st_size:
        failures.append({"status": "sampler_stderr_nonempty_or_missing"})
        state = "failed"
    rolled_back: List[int] = []
    if state != "success":
        rolled_back = rollback_segment_outputs(root, tasks)
    final = write_segment_final(
        root, segment, intent, state, failures=failures,
        rolled_back=rolled_back, process_actions=process_actions,
        jobs_ended_monotonic_s=jobs_ended,
        sampler_returncode=sampler_returncode,
    )
    if caught is not None:
        raise caught
    return final


def command_launch(args: argparse.Namespace) -> None:
    root = Path(args.root).resolve()
    design, tasks = verify_root(
        root, args.expected_receipts_sha256, runtime_tolerant=True,
    )
    receipt_partials = sorted(
        path for directory in (
            root / "segments", root / "job_receipts", root / "raw",
            root / "stderr", root / "exit", root / "attempts",
        )
        for path in directory.rglob("*.part.*")
    )
    indices = segment_indices(root)
    incomplete = [
        segment for segment in indices
        if not segment_path(root, segment, "final").is_file()
    ]
    needs_reconciliation = bool(incomplete or receipt_partials)
    if args.preflight_only:
        samplers = other_samplers()
        complete_count: Optional[int] = None
        if not needs_reconciliation:
            complete_count = len(completed_jobs(root, tasks))
        print(json.dumps({
            "status": (
                "RESUME_RECONCILIATION_REQUIRED" if needs_reconciliation
                else "BLOCKED_CONCURRENT_SAMPLER" if samplers
                else "READY_NOT_LAUNCHED"),
            "root": str(root), "tasks": len(tasks), "complete": complete_count,
            "remaining": None if complete_count is None else len(tasks) - complete_count,
            "workers": design["worker_count"], "requested_stage": args.stage,
            "single_sampler": not bool(samplers), "other_samplers": samplers,
            "incomplete_segments": incomplete,
            "orphan_receipt_partials": [str(path.relative_to(root)) for path in receipt_partials],
        }, indent=2, sort_keys=True))
        return
    if needs_reconciliation and not args.resume:
        die("interrupted segment(s) require --resume reconciliation")
    if indices and not args.resume:
        die("existing campaign segment history requires --resume")
    reconciled: List[dict] = []
    if args.resume:
        for path in receipt_partials:
            path.unlink()
        reconciled = reconcile_incomplete_segments(root, tasks)
    design, tasks = verify_root(root, args.expected_receipts_sha256)
    complete = completed_jobs(root, tasks)
    remaining = [task for task in tasks if int(task["job"]) not in complete]
    by_stage = {
        stage: [task for task in remaining if task["stage"] == stage]
        for stage in ("hard", "control")
    }
    samplers = other_samplers()
    if samplers:
        die("refusing concurrent SPD/I2C reader(s): {}".format(
            json.dumps(samplers, sort_keys=True)))
    if complete and not args.resume:
        die("outputs already exist; inspect then pass --resume")
    if not remaining:
        die("campaign ledger is already complete")
    hard_incomplete = bool(by_stage["hard"])
    if args.stage == "control" and hard_incomplete:
        die("control stage requires the complete hard-stage Cartesian product")
    requested_stages = ("hard", "control") if args.stage == "all" else (args.stage,)
    if not any(by_stage[stage] for stage in requested_stages):
        die("requested campaign stage is already complete")
    launched: List[dict] = []
    for stage in requested_stages:
        stage_tasks = by_stage[stage]
        if not stage_tasks:
            continue
        # Recheck immediately before each privileged I2C sampler launch.
        samplers = other_samplers()
        if samplers:
            die("refusing concurrent SPD/I2C reader(s): {}".format(
                json.dumps(samplers, sort_keys=True)))
        final = run_segment(
            root, design, stage_tasks, stage, len(complete), bool(args.resume),
        )
        launched.append(final)
        print(json.dumps(final, indent=2, sort_keys=True))
        if final.get("state") != "success":
            raise SystemExit(1)
        complete.update(int(task["job"]) for task in stage_tasks)
    print(json.dumps({
        "status": "LAUNCH_OK", "reconciled_segments": len(reconciled),
        "launched_segments": [record["segment"] for record in launched],
        "complete": len(complete), "remaining": len(tasks) - len(complete),
    }, indent=2, sort_keys=True))


def decimal_field(row: Mapping[str, str], key: str) -> Decimal:
    try:
        value = Decimal(row[key])
    except (KeyError, InvalidOperation, TypeError) as exc:
        die("invalid decimal field {}: {}".format(key, exc))
    if not value.is_finite():
        die("nonfinite decimal field {}".format(key))
    return value


def classify(row: Mapping[str, str]) -> str:
    try:
        success, rank_fail, error = (
            int(row[key]) for key in ("success", "rank_fail", "error")
        )
    except (KeyError, TypeError, ValueError) as exc:
        die("invalid benchmark outcome field: {}".format(exc))
    if any(value not in (0, 1) for value in (success, rank_fail, error)) or \
            success + rank_fail + error != 1 or error:
        die("invalid benchmark outcome")
    if success:
        return "success"
    try:
        return "q>H" if int(row["binary_def_max"]) > 12 else "field_shortfall"
    except (KeyError, TypeError, ValueError) as exc:
        die("invalid benchmark binary deficit: {}".format(exc))


def validate_attempt_manifest(
    root: Path,
    segment: int,
    final: Mapping[str, object],
    allowed_jobs: Set[int],
    *,
    require_complete: bool,
) -> None:
    attempt_dir = root / "attempts" / "segment{:03d}".format(segment)
    if not attempt_dir.is_dir() or attempt_dir.is_symlink():
        die("segment attempt directory is missing or indirect")
    files = sorted(path for path in attempt_dir.iterdir() if path.is_file())
    entries = list(attempt_dir.iterdir())
    if any(path.is_symlink() or not path.is_file() for path in entries):
        die("segment attempt inventory contains a symlink")
    allowed = re.compile(r"job([0-9]{5})\.(?:pid|stdout|stderr|exit)$")
    matches = [allowed.fullmatch(path.name) for path in files]
    if any(match is None for match in matches):
        die("segment attempt inventory contains an unexpected artifact")
    groups: Dict[int, Set[str]] = defaultdict(set)
    for path, match in zip(files, matches):
        assert match is not None
        job = int(match.group(1))
        if job not in allowed_jobs:
            die("segment attempt inventory names a job outside its intent")
        suffix = path.suffix[1:]
        groups[job].add(suffix)
        if suffix == "pid":
            identity = load_json(path)
            value = identity.get("pid")
            start_ticks = identity.get("start_ticks")
            if not isinstance(value, int) or isinstance(value, bool) or value <= 1 or \
                    not isinstance(start_ticks, int) or isinstance(start_ticks, bool) or \
                    start_ticks < 0:
                die("segment attempt PID receipt is out of range")
    complete_suffixes = {"pid", "stdout", "stderr", "exit"}
    for job, suffixes in groups.items():
        if suffixes != {"pid"} and suffixes != complete_suffixes:
            die("segment attempt contains a partial completed artifact set")
        if suffixes == complete_suffixes:
            exit_path = attempt_dir / "job{:05d}.exit".format(job)
            try:
                int(exit_path.read_text(encoding="ascii").strip())
            except (OSError, ValueError):
                die("segment attempt exit receipt is malformed")
    if require_complete and (set(groups) != allowed_jobs or
                             any(value != complete_suffixes for value in groups.values())):
        die("successful segment lacks complete attempt evidence for every job")
    data = b"".join(
        "{}  {}\n".format(sha256_file(path), path.relative_to(root)).encode("ascii")
        for path in files
    )
    manifest = root / "segments" / "segment{:03d}.attempts.sha256".format(segment)
    if not manifest.is_file() or manifest.is_symlink() or manifest.read_bytes() != data or \
            final.get("attempt_manifest") != manifest.name or \
            final.get("attempt_manifest_sha256") != sha256_bytes(data) or \
            final.get("attempt_file_count") != len(files):
        die("segment attempt manifest mismatch")


def validate_thermal_rows(
    rows: Sequence[Mapping[str, str]],
    *,
    segment: int,
) -> Dict[str, object]:
    if len(rows) < THERMAL_READY_SAMPLES:
        die("ready segment {} has fewer than two thermal rows".format(segment))
    monotonic: List[Decimal] = []
    busy: List[Decimal] = []
    cpu: List[Decimal] = []
    dimm_by_field: Dict[str, List[Decimal]] = {
        key: [] for key in THERMAL_FIELDS[5:13]
    }
    utc_previous: Optional[datetime] = None
    for row in rows:
        try:
            utc = datetime.fromisoformat(row["utc"].replace("Z", "+00:00"))
        except (KeyError, ValueError) as exc:
            die("invalid thermal UTC timestamp: {}".format(exc))
        if utc.tzinfo is None or (utc_previous is not None and utc <= utc_previous):
            die("thermal UTC timestamps are not strictly increasing")
        utc_previous = utc
        current = decimal_field(row, "monotonic_s")
        if monotonic and current <= monotonic[-1]:
            die("thermal monotonic timestamps are not strictly increasing")
        if monotonic and not (
                THERMAL_MIN_GAP_SECONDS <= current - monotonic[-1] <=
                THERMAL_MAX_GAP_SECONDS):
            die("thermal cadence gap is outside policy in segment {}".format(segment))
        monotonic.append(current)
        if row["cpu_busy_pct"]:
            value = decimal_field(row, "cpu_busy_pct")
            if value < 0 or value > 100:
                die("thermal CPU busy percentage is out of range")
            busy.append(value)
        if not row["cpu_tctl_c"]:
            die("thermal row lacks CPU Tctl")
        cpu.append(decimal_field(row, "cpu_tctl_c"))
        for key in THERMAL_FIELDS[5:13]:
            if not row[key]:
                die("thermal row lacks DIMM field {}".format(key))
            dimm_by_field[key].append(decimal_field(row, key))
        try:
            dimm_errors = int(row["dimm_read_errors"])
            edac_ce = int(row["edac_ce"])
            edac_ue = int(row["edac_ue"])
        except (KeyError, ValueError) as exc:
            die("invalid thermal hardware counter: {}".format(exc))
        if dimm_errors != 0 or edac_ce != 0 or edac_ue != 0:
            die("thermal hardware error receipt is nonzero")
    return {
        "samples": len(rows), "monotonic": monotonic, "busy": busy,
        "cpu": cpu, "dimm_by_field": dimm_by_field,
    }


def validate_thermal(root: Path, design: Mapping[str, object], tasks: Sequence[dict]) -> Dict[str, object]:
    indices = segment_indices(root)
    if not indices:
        die("campaign has no launch segments")
    expected_segment_files: Set[str] = set()
    successful_jobs: Set[int] = set()
    successful_segments = failed_segments = interrupted_segments = 0
    samples = 0
    all_busy: List[Decimal] = []
    all_cpu: List[Decimal] = []
    dimm_maxima: Dict[str, Decimal] = {}
    seen_control_success = False
    expected_thermal_files: Set[str] = set()
    task_by_job = {int(task["job"]): task for task in tasks}
    for segment in indices:
        intent_path = segment_path(root, segment, "intent")
        ready_path = segment_path(root, segment, "ready")
        final_path = segment_path(root, segment, "final")
        manifest_path = root / "segments" / "segment{:03d}.attempts.sha256".format(segment)
        expected_segment_files.update((intent_path.name, final_path.name, manifest_path.name))
        if not intent_path.is_file() or intent_path.is_symlink() or \
                not final_path.is_file() or final_path.is_symlink():
            die("campaign contains an unreconciled segment")
        intent = load_json(intent_path)
        final = load_json(final_path)
        selected = segment_jobs(intent, tasks)
        jobs = [int(task["job"]) for task in selected]
        stage = str(intent.get("stage"))
        expected_jobs = [
            int(task["job"]) for task in tasks
            if task["stage"] == stage and int(task["job"]) not in successful_jobs
        ]
        hard_remaining = any(
            task["stage"] == "hard" and int(task["job"]) not in successful_jobs
            for task in tasks
        )
        if intent.get("schema") != SCHEMA + ".segment_intent" or \
                intent.get("segment") != segment or \
                intent.get("workers") != design.get("worker_count") or \
                intent.get("previously_complete") != len(successful_jobs) or \
                jobs != expected_jobs or (stage == "control" and hard_remaining) or \
                intent.get("retry_policy") != "stage-atomic non-selective retry" or \
                final.get("schema") != SCHEMA + ".segment_final" or \
                final.get("segment") != segment or \
                final.get("intent_sha256") != sha256_file(intent_path) or \
                final.get("jobs") != jobs or final.get("jobs_sha256") != intent.get("jobs_sha256"):
            die("segment intent/final binding mismatch")
        state = final.get("state")
        if state not in ("success", "failed", "interrupted"):
            die("unknown segment terminal state")
        rolled_back = final.get("rolled_back_jobs")
        published = final.get("published_jobs")
        failures = final.get("failures")
        if not isinstance(rolled_back, list) or \
                any(not isinstance(job, int) or isinstance(job, bool) or job not in jobs
                    for job in rolled_back) or rolled_back != sorted(set(rolled_back)) or \
                not isinstance(published, list) or \
                not isinstance(failures, list):
            die("segment terminal retry ledger is malformed")
        validate_attempt_manifest(
            root, segment, final, set(jobs), require_complete=state == "success",
        )
        ready: Optional[Dict[str, object]] = None
        if ready_path.exists():
            expected_segment_files.add(ready_path.name)
            if ready_path.is_symlink() or not ready_path.is_file():
                die("segment ready receipt is indirect")
            ready = load_json(ready_path)
            samples_at_ready = ready.get("samples_at_ready")
            if ready.get("schema") != SCHEMA + ".segment_ready" or \
                    ready.get("segment") != segment or \
                    ready.get("intent_sha256") != sha256_file(intent_path) or \
                    final.get("ready_sha256") != sha256_file(ready_path) or \
                    not isinstance(samples_at_ready, int) or isinstance(samples_at_ready, bool) or \
                    samples_at_ready < THERMAL_READY_SAMPLES or \
                    not isinstance(ready.get("sampler_pid"), int) or \
                    isinstance(ready.get("sampler_pid"), bool) or ready["sampler_pid"] <= 1 or \
                    not isinstance(ready.get("sampler_start_ticks"), int) or \
                    isinstance(ready.get("sampler_start_ticks"), bool) or \
                    ready["sampler_start_ticks"] < 0 or \
                    not isinstance(ready.get("jobs_started_monotonic_s"), (int, float)) or \
                    isinstance(ready.get("jobs_started_monotonic_s"), bool):
                die("segment readiness receipt mismatch")
        elif final.get("ready_sha256") is not None:
            die("segment final references a missing readiness receipt")
        csv_path = root / "thermal" / "segment{:03d}.csv".format(segment)
        stderr_path = root / "thermal" / "segment{:03d}.stderr".format(segment)
        pid_path = root / "thermal" / "segment{:03d}.pid".format(segment)
        if pid_path.exists() or pid_path.is_symlink():
            die("finalized segment retains a thermal PID file")
        csv_exists = csv_path.is_file() and not csv_path.is_symlink()
        stderr_exists = stderr_path.is_file() and not stderr_path.is_symlink()
        if csv_exists != (final.get("thermal_csv_sha256") is not None) or \
                stderr_exists != (final.get("thermal_stderr_sha256") is not None):
            die("segment thermal artifact cardinality mismatch")
        if ready is not None and (not csv_exists or not stderr_exists):
            die("ready segment lacks complete thermal artifacts")
        if csv_exists:
            expected_thermal_files.add(csv_path.name)
        if stderr_exists:
            expected_thermal_files.add(stderr_path.name)
        thermal: Optional[Dict[str, object]] = None
        rows: List[Dict[str, str]] = []
        if csv_exists:
            if final.get("thermal_csv") != csv_path.name or \
                    final.get("thermal_csv_sha256") != sha256_file(csv_path) or \
                    (stderr_exists and final.get("thermal_stderr") != stderr_path.name):
                die("segment thermal hash binding mismatch")
            if ready is not None and stderr_path.stat().st_size:
                die("thermal sampler stderr is nonempty")
            rows = read_thermal_rows(csv_path)
            if ready is not None:
                thermal = validate_thermal_rows(rows, segment=segment)
                count = int(ready["samples_at_ready"])
                if count > len(rows) or \
                        ready.get("last_sample_monotonic_s") != rows[count - 1]["monotonic_s"] or \
                        Decimal(str(ready["jobs_started_monotonic_s"])) < \
                            decimal_field(rows[count - 1], "monotonic_s") or \
                        Decimal(str(ready["jobs_started_monotonic_s"])) - \
                            decimal_field(rows[count - 1], "monotonic_s") > THERMAL_MAX_GAP_SECONDS:
                    die("segment readiness is not bound to its live thermal prefix")
                samples += int(thermal["samples"])
                all_cpu.extend(thermal["cpu"])
                for key, values in thermal["dimm_by_field"].items():
                    dimm_maxima[key] = max(dimm_maxima.get(key, values[0]), max(values))
        if state == "success":
            successful_segments += 1
            end_value = final.get("jobs_ended_monotonic_s")
            if ready is None or thermal is None or \
                    not isinstance(final.get("sampler_returncode"), int) or \
                    isinstance(final.get("sampler_returncode"), bool) or \
                    final.get("sampler_returncode") != 0 or \
                    final.get("failures") != [] or final.get("rolled_back_jobs") != [] or \
                    final.get("published_jobs") != jobs or \
                    not isinstance(end_value, (int, float)) or isinstance(end_value, bool):
                die("successful segment terminal receipt is incomplete")
            start = Decimal(str(ready["jobs_started_monotonic_s"]))
            end = Decimal(str(final["jobs_ended_monotonic_s"]))
            monotonic = thermal["monotonic"]
            if end < start or monotonic[0] > start or monotonic[-1] < end:
                die("successful thermal segment lacks full launch start/end coverage")
            # Include the first sample at/after the end boundary and the sample
            # immediately before start, since CPU busy is interval-derived.
            first_after_end = next(index for index, value in enumerate(monotonic) if value >= end)
            first_after_start = next(index for index, value in enumerate(monotonic) if value >= start)
            coverage_set = set(range(max(0, first_after_start - 1), first_after_end + 1))
            coverage_busy = [
                decimal_field(rows[index], "cpu_busy_pct")
                for index in sorted(coverage_set)
                if rows[index]["cpu_busy_pct"]
            ]
            if not coverage_busy:
                die("successful segment has no CPU busy coverage")
            if stage in CPU_BUSY_FLOOR_STAGES and \
                    sum(coverage_busy) / len(coverage_busy) < CPU_BUSY_FLOOR:
                die("successful segment CPU busy mean is below the sealed floor")
            all_busy.extend(coverage_busy)
            if intent.get("stage") == "control":
                seen_control_success = True
            elif seen_control_success:
                die("successful hard segment appears after a successful control segment")
            overlap = successful_jobs.intersection(jobs)
            if overlap:
                die("successful job appears in more than one segment")
            for job in jobs:
                receipt = load_json(job_receipt_path(root, job))
                task = task_by_job[job]
                name = task["output_name"]
                if receipt.get("segment") != segment or \
                        Decimal(str(receipt["started_monotonic_s"])) < start or \
                        Decimal(str(receipt["ended_monotonic_s"])) > end or \
                        sha256_file(root / "attempts" / "segment{:03d}".format(segment) /
                                    "job{:05d}.stdout".format(job)) != \
                            sha256_file(root / "raw" / name) or \
                        sha256_file(root / "attempts" / "segment{:03d}".format(segment) /
                                    "job{:05d}.stderr".format(job)) != \
                            sha256_file(root / "stderr" / (name + ".stderr")) or \
                        sha256_file(root / "attempts" / "segment{:03d}".format(segment) /
                                    "job{:05d}.exit".format(job)) != \
                            sha256_file(root / "exit" / (name + ".exit")):
                    die("job receipt is bound to the wrong successful segment")
            successful_jobs.update(jobs)
        else:
            if state == "failed":
                failed_segments += 1
            else:
                interrupted_segments += 1
            if final.get("published_jobs") != []:
                die("non-success segment retained published jobs")
            if state == "failed" and not failures:
                die("failed segment lacks a failure receipt")
            for job in jobs:
                receipt_path = job_receipt_path(root, job)
                if receipt_path.is_file() and load_json(receipt_path).get("segment") == segment:
                    die("non-success segment retained a job receipt")
        if stderr_exists and final.get("thermal_stderr") != stderr_path.name or \
                (stderr_exists and final.get("thermal_stderr_sha256") != sha256_file(stderr_path)):
            die("segment thermal stderr hash binding mismatch")
    actual_segment_files = {path.name for path in (root / "segments").iterdir()}
    if actual_segment_files != expected_segment_files:
        die("segment receipt inventory contains missing or unexpected files")
    actual_thermal_files = {path.name for path in (root / "thermal").iterdir()}
    if actual_thermal_files != expected_thermal_files:
        die("thermal segment inventory contains missing or unexpected files")
    expected_attempt_dirs = {"segment{:03d}".format(index) for index in indices}
    attempt_entries = list((root / "attempts").iterdir())
    if {path.name for path in attempt_entries} != expected_attempt_dirs or \
            any(not path.is_dir() or path.is_symlink() for path in attempt_entries):
        die("attempt segment directory inventory mismatch")
    completed = completed_jobs(root, tasks)
    if successful_jobs != completed:
        die("successful segment/job receipt coverage mismatch")
    if completed != set(task_by_job):
        die("successful segment receipts do not cover the complete campaign")
    expected_job_receipts = {
        job_receipt_path(root, job).name for job in completed
    }
    job_receipt_entries = list((root / "job_receipts").iterdir())
    if {path.name for path in job_receipt_entries} != expected_job_receipts or \
            any(not path.is_file() or path.is_symlink() for path in job_receipt_entries):
        die("job receipt inventory mismatch")
    if not all_busy:
        die("campaign has no successful thermal load samples")
    return {
        "segments": len(indices), "successful_segments": successful_segments,
        "failed_segments": failed_segments,
        "interrupted_segments": interrupted_segments, "samples": samples,
        "successful_busy_mean": str(sum(all_busy) / len(all_busy)),
        "cpu_max_c": str(max(all_cpu)) if all_cpu else None,
        "dimm_max_c_by_field": {
            key: str(dimm_maxima[key]) for key in sorted(dimm_maxima)
        },
    }


def command_reduce(args: argparse.Namespace) -> None:
    root = Path(args.root).resolve()
    design, tasks = verify_root(root, args.expected_receipts_sha256)
    if len(completed_jobs(root, tasks)) != len(tasks):
        die("campaign is incomplete")
    thermal = validate_thermal(root, design, tasks)
    outcomes: Dict[Tuple[int, int, str, int, int, str], dict] = {}
    manifest_lines = bytearray()
    runtime_paths: List[Path] = []
    for directory_name in ("segments", "thermal", "attempts", "job_receipts"):
        directory = root / directory_name
        for path in directory.rglob("*"):
            if path.is_symlink():
                die("runtime evidence contains a symlink")
            if path.is_file():
                runtime_paths.append(path)
    for path in sorted(runtime_paths, key=lambda value: str(value.relative_to(root))):
        manifest_lines.extend("{}  {}\n".format(
            sha256_file(path), path.relative_to(root)).encode("ascii"))
    for task in tasks:
        name = task["output_name"]
        raw_path = root / "raw" / name
        stderr_path = root / "stderr" / (name + ".stderr")
        exit_path = root / "exit" / (name + ".exit")
        for path in (raw_path, stderr_path, exit_path):
            manifest_lines.extend("{}  {}\n".format(
                sha256_file(path), path.relative_to(root)).encode("ascii"))
        rows = parse_benchmark_csv(raw_path.read_bytes(), task)
        for row in rows:
            key = (int(task["period"]), int(task["mask"]), str(task["stage"]),
                   int(row["N"]), int(task["seed_index"]), str(task["schedule"]))
            if key in outcomes:
                die("duplicate Cartesian cell {}".format(key))
            outcomes[key] = {
                "outcome": classify(row),
                "xors": decimal_field(row, "block_xors_mu"),
                "muladds": decimal_field(row, "block_muladds_mu"),
                "seed_attempt": int(row["seed_attempt"]),
                "inactivated": int(row["inact_max"]),
                "binary_deficit": int(row["binary_def_max"]),
            }
    if len(outcomes) != int(design["total_cells"]):
        die("Cartesian result count mismatch")

    hard = read_jsonl(root / "hard_keys.jsonl")
    controls = {48: read_jsonl(root / "controls_p48.jsonl"),
                32: read_jsonl(root / "controls_p32.jsonl")}
    # The named suffix controls must reproduce the source campaign exactly.
    # This catches an accidental graph, seed-attempt, schedule, or equation
    # change before any alternate mask is scored.
    for panel in PANELS:
        period = int(panel["period"])
        suffix = int(panel["suffix_mask"])
        for cell in hard:
            K, seed_index, schedule = cell_tuple(cell)
            observed = outcomes[(period, suffix, "hard", K, seed_index, schedule)]["outcome"]
            expected = "field_shortfall" if panel["name"] in cell["source_arms"] else "success"
            if observed != expected:
                die("canonical suffix failed exact all-K hard-cell replay")
        for cell in controls[period]:
            K, seed_index, schedule = cell_tuple(cell)
            if outcomes[(period, suffix, "control", K, seed_index, schedule)]["outcome"] != "success":
                die("canonical suffix failed a sealed common-success control")
    results: Dict[str, List[dict]] = {}
    for panel in PANELS:
        period = int(panel["period"])
        suffix = int(panel["suffix_mask"])
        panel_rows: List[dict] = []
        for mask in read_jsonl(root / "masks_p{}_r{}.jsonl".format(period, panel["rows"])):
            mask_value = int(mask["mask"])
            counts = defaultdict(int)
            repairs = introductions = 0
            candidate_work = [Decimal(0), Decimal(0)]
            suffix_work = [Decimal(0), Decimal(0)]
            common_success = 0
            weak_K: Set[int] = set()
            weak_multiplicity: Dict[int, int] = defaultdict(int)
            by_seed: Dict[int, int] = defaultdict(int)
            by_schedule: Dict[str, int] = defaultdict(int)
            failures: List[dict] = []
            repair_keys: List[dict] = []
            introduction_keys: List[dict] = []
            receipt_deltas = {
                field: {"cells": 0, "candidate_minus_suffix_sum": 0, "max_abs": 0}
                for field in ("seed_attempt", "inactivated", "binary_deficit")
            }
            for stage, source in (("hard", hard), ("control", controls[period])):
                for cell in source:
                    K, seed_index, schedule = cell_tuple(cell)
                    candidate = outcomes[(period, mask_value, stage, K, seed_index, schedule)]
                    canonical = outcomes[(period, suffix, stage, K, seed_index, schedule)]
                    # Configuration selection is inherited from the codec, with
                    # no experiment-specific seed patches.  A mask can therefore
                    # alter the selected attempt or solve dimensions; preserve
                    # those paired deltas as results instead of rejecting them.
                    for field, stats in receipt_deltas.items():
                        delta = candidate[field] - canonical[field]
                        if delta:
                            stats["cells"] += 1
                            stats["candidate_minus_suffix_sum"] += delta
                            stats["max_abs"] = max(stats["max_abs"], abs(delta))
                    counts[stage + "_" + candidate["outcome"]] += 1
                    if candidate["outcome"] != "success":
                        weak_K.add(K)
                        weak_multiplicity[K] += 1
                        by_seed[seed_index] += 1
                        by_schedule[schedule] += 1
                        failures.append({
                            "stage": stage, **cell_record((K, seed_index, schedule)),
                            "cause": candidate["outcome"],
                        })
                    if canonical["outcome"] != "success" and candidate["outcome"] == "success":
                        repairs += 1
                        repair_keys.append({"stage": stage, **cell_record((K, seed_index, schedule))})
                    if canonical["outcome"] == "success" and candidate["outcome"] != "success":
                        introductions += 1
                        introduction_keys.append({
                            "stage": stage, **cell_record((K, seed_index, schedule))})
                    if candidate["outcome"] == canonical["outcome"] == "success":
                        common_success += 1
                        candidate_work[0] += candidate["xors"]
                        candidate_work[1] += candidate["muladds"]
                        suffix_work[0] += canonical["xors"]
                        suffix_work[1] += canonical["muladds"]
            panel_rows.append({
                "mask": mask_value, "mask_hex": mask["mask_hex"],
                "canonical_suffix": bool(mask["canonical_suffix"]),
                "counts": dict(sorted(counts.items())), "weak_K": len(weak_K),
                "weak_K_multiplicity": {
                    str(K): weak_multiplicity[K] for K in sorted(weak_multiplicity)
                },
                "failures_by_seed": {
                    str(seed): by_seed[seed] for seed in sorted(by_seed)
                },
                "failures_by_schedule": {
                    schedule: by_schedule[schedule] for schedule in sorted(by_schedule)
                },
                "failure_records": failures,
                "repairs_vs_suffix": repairs, "introductions_vs_suffix": introductions,
                "repair_keys_vs_suffix": repair_keys,
                "introduction_keys_vs_suffix": introduction_keys,
                "paired_receipt_deltas_vs_suffix": receipt_deltas,
                "common_success": common_success,
                "common_success_xor_ratio": (
                    str(candidate_work[0] / suffix_work[0]) if suffix_work[0] else None),
                "common_success_muladd_ratio": (
                    str(candidate_work[1] / suffix_work[1]) if suffix_work[1] else None),
            })
        if len(panel_rows) != MASK_COUNT:
            die("reducer mask cardinality mismatch")
        results[str(period)] = panel_rows
    data_manifest = bytes(manifest_lines)
    data_manifest_sha256 = sha256_bytes(data_manifest)
    write_once(root / "data_manifest.sha256", data_manifest)
    summary = {
        "schema": SCHEMA + ".validated", "validation_issue_count": 0,
        "root": str(root), "prelaunch_receipts_sha256": args.expected_receipts_sha256,
        "data_manifest_sha256": data_manifest_sha256, "design": design,
        "hard_keys": hard, "thermal": thermal, "results": results,
    }
    write_once(root / "validated_summary.json", canonical_json(summary))
    print(json.dumps({
        "status": "VALIDATION_OK", "cells": len(outcomes),
        "validated_summary_sha256": sha256_file(root / "validated_summary.json"),
        "data_manifest_sha256": data_manifest_sha256,
    }, indent=2, sort_keys=True))


def command_selftest(_: argparse.Namespace) -> None:
    masks3 = enumerate_masks(48, 3, 0x380)
    masks7 = enumerate_masks(32, 7, 0x3F8)
    if {int(row["mask"]) for row in masks3} != {
            sum(1 << bit for bit in choice) for choice in itertools.combinations(range(10), 3)}:
        die("r3 mask selftest failed")
    if {int(row["mask"]) for row in masks7} != {
            sum(1 << bit for bit in choice) for choice in itertools.combinations(range(10), 7)}:
        die("r7 mask selftest failed")
    excluded = {(2, 0, "burst"), (3, 1, "adversarial")}
    a = deterministic_controls(48, excluded, count=32, k_lo=2, k_hi=100)
    b = deterministic_controls(48, excluded, count=32, k_lo=2, k_hi=100)
    c = deterministic_controls(32, excluded, count=32, k_lo=2, k_hi=100)
    if a != b or a == c or any(cell_tuple(row) in excluded for row in a + c):
        die("deterministic control selftest failed")
    hard = [cell_record((K, K % 3, SCHEDULES[K % 3])) for K in range(2, 18)]
    hard_set = {cell_tuple(row) for row in hard}
    controls = {
        48: deterministic_controls(48, hard_set, count=32, k_lo=20, k_hi=100),
        32: deterministic_controls(32, hard_set, count=32, k_lo=20, k_hi=100),
    }
    tasks = build_tasks(
        "/frozen/wirehair_v2_bench", {48: masks3, 32: masks7},
        hard, controls, control_count=32,
    )
    if [task["job"] for task in tasks] != list(range(len(tasks))) or \
            len({task["output_name"] for task in tasks}) != len(tasks) or \
            sum(len(task["Ks"]) for task in tasks) != 2 * 120 * 48:
        die("task-ledger identity/cardinality selftest failed")
    stages = [task["stage"] for task in tasks]
    if stages != sorted(stages, key=("hard", "control").index):
        die("task-ledger stage ordering selftest failed")
    cartesian = {
        (task["period"], task["mask"], task["stage"], K,
         task["seed_index"], task["schedule"])
        for task in tasks for K in task["Ks"]
    }
    if len(cartesian) != 2 * 120 * 48:
        die("task-ledger Cartesian uniqueness selftest failed")
    if parse_preamble("# precodefail: a=1 b=two") != {"a": "1", "b": "two"}:
        die("preamble selftest failed")
    try:
        parse_preamble("# precodefail: duplicate=1 duplicate=2")
    except CampaignError:
        pass
    else:
        die("duplicate preamble rejection selftest failed")

    def expect_reject(function: Any, label: str) -> None:
        try:
            function()
        except CampaignError:
            return
        die("{} tamper selftest was accepted".format(label))

    def thermal_bytes(monotonic_values: Sequence[str], *, edac_ce: str = "0",
                      blank_dimm: bool = False, busy: str = "99.5") -> bytes:
        output = io.StringIO(newline="")
        writer = csv.DictWriter(output, fieldnames=THERMAL_FIELDS, lineterminator="\n")
        writer.writeheader()
        for index, monotonic_value in enumerate(monotonic_values):
            row = {key: "0" for key in THERMAL_FIELDS}
            row.update({
                "utc": "2026-07-18T00:00:{:02d}.000Z".format(index),
                "monotonic_s": monotonic_value, "cpu_busy_pct": busy,
                "cpu_avg_mhz": "5000", "cpu_tctl_c": "60",
                "dimm_read_errors": "0", "load1": "128", "load5": "128",
                "load15": "128", "edac_ce": edac_ce, "edac_ue": "0",
            })
            for key in THERMAL_FIELDS[5:13]:
                row[key] = "45"
            if blank_dimm:
                row[THERMAL_FIELDS[5]] = ""
            writer.writerow(row)
        return output.getvalue().encode("ascii")

    def make_success_fixture(
            root: Path, *, stage: str = "hard",
            busy: str = "99.5") -> Tuple[List[dict], dict]:
        for directory in (
                "raw", "stderr", "exit", "thermal", "segments",
                "attempts", "job_receipts"):
            (root / directory).mkdir(parents=True, exist_ok=True)
        task = {"job": 0, "stage": stage, "output_name": "job.csv", "argv": ["/bin/true"]}
        intent = {
            "schema": SCHEMA + ".segment_intent", "segment": 0,
            "created_utc": "2026-07-18T00:00:00+00:00",
            "created_monotonic_s": 99.5, "stage": stage, "jobs": [0],
            "jobs_sha256": sha256_bytes(json_lines([0])), "workers": 1,
            "previously_complete": 0, "resume": False,
            "retry_policy": "stage-atomic non-selective retry",
        }
        write_once(segment_path(root, 0, "intent"), canonical_json(intent))
        (root / "attempts" / "segment000").mkdir()
        (root / "attempts" / "segment000" / "job00000.pid").write_bytes(
            canonical_json({"pid": 99999999, "start_ticks": 1}))
        ready = {
            "schema": SCHEMA + ".segment_ready", "segment": 0,
            "intent_sha256": sha256_file(segment_path(root, 0, "intent")),
            "ready_utc": "2026-07-18T00:00:00+00:00",
            "jobs_started_monotonic_s": 101.1, "sampler_pid": 99999999,
            "sampler_start_ticks": 1,
            "samples_at_ready": 2, "last_sample_monotonic_s": "101.0",
        }
        write_once(segment_path(root, 0, "ready"), canonical_json(ready))
        thermal_csv = root / "thermal" / "segment000.csv"
        thermal_stderr = root / "thermal" / "segment000.stderr"
        thermal_csv.write_bytes(thermal_bytes(("100.0", "101.0", "102.0"), busy=busy))
        thermal_stderr.write_bytes(b"")
        raw = root / "raw" / "job.csv"
        stderr = root / "stderr" / "job.csv.stderr"
        exit_path = root / "exit" / "job.csv.exit"
        raw.write_bytes(b"synthetic output\n")
        stderr.write_bytes(b"")
        exit_path.write_bytes(b"0\n")
        (root / "attempts" / "segment000" / "job00000.stdout").write_bytes(raw.read_bytes())
        (root / "attempts" / "segment000" / "job00000.stderr").write_bytes(b"")
        (root / "attempts" / "segment000" / "job00000.exit").write_bytes(b"0\n")
        write_once(job_receipt_path(root, 0), canonical_json({
            "schema": SCHEMA + ".job", "job": 0, "segment": 0,
            "stage": stage, "output_name": "job.csv", "returncode": 0,
            "started_monotonic_s": 101.2, "ended_monotonic_s": 101.4,
            "stdout_sha256": sha256_file(raw), "stderr_sha256": sha256_file(stderr),
            "exit_sha256": sha256_file(exit_path),
        }))
        write_segment_final(
            root, 0, intent, "success", failures=[], rolled_back=[],
            process_actions=[], jobs_ended_monotonic_s=101.5,
            sampler_returncode=0,
        )
        return [task], {"worker_count": 1}

    with tempfile.TemporaryDirectory(prefix="wirehair-row-mask-selftest-") as temporary:
        fixture = Path(temporary) / "success"
        fixture.mkdir()
        fixture_tasks, fixture_design = make_success_fixture(fixture)
        receipt = validate_thermal(fixture, fixture_design, fixture_tasks)
        if receipt["successful_segments"] != 1 or receipt["samples"] != 3:
            die("synthetic segment validation selftest failed")

        def tampered_copy(name: str) -> Path:
            target = Path(temporary) / name
            shutil.copytree(fixture, target)
            return target

        cadence = tampered_copy("cadence")
        cadence_csv = cadence / "thermal" / "segment000.csv"
        cadence_csv.write_bytes(thermal_bytes(("100.0", "101.0", "104.0")))
        cadence_final = load_json(segment_path(cadence, 0, "final"))
        cadence_final["thermal_csv_sha256"] = sha256_file(cadence_csv)
        segment_path(cadence, 0, "final").write_bytes(canonical_json(cadence_final))
        expect_reject(
            lambda: validate_thermal(cadence, fixture_design, fixture_tasks),
            "thermal cadence",
        )

        dimm = tampered_copy("dimm")
        dimm_csv = dimm / "thermal" / "segment000.csv"
        dimm_csv.write_bytes(thermal_bytes(("100.0", "101.0", "102.0"), blank_dimm=True))
        dimm_final = load_json(segment_path(dimm, 0, "final"))
        dimm_final["thermal_csv_sha256"] = sha256_file(dimm_csv)
        segment_path(dimm, 0, "final").write_bytes(canonical_json(dimm_final))
        expect_reject(
            lambda: validate_thermal(dimm, fixture_design, fixture_tasks),
            "missing DIMM",
        )

        edac = tampered_copy("edac")
        edac_csv = edac / "thermal" / "segment000.csv"
        edac_csv.write_bytes(thermal_bytes(("100.0", "101.0", "102.0"), edac_ce="1"))
        edac_final = load_json(segment_path(edac, 0, "final"))
        edac_final["thermal_csv_sha256"] = sha256_file(edac_csv)
        segment_path(edac, 0, "final").write_bytes(canonical_json(edac_final))
        expect_reject(
            lambda: validate_thermal(edac, fixture_design, fixture_tasks),
            "EDAC",
        )

        coverage = tampered_copy("coverage")
        coverage_ready = load_json(segment_path(coverage, 0, "ready"))
        coverage_ready["jobs_started_monotonic_s"] = 99.0
        segment_path(coverage, 0, "ready").write_bytes(canonical_json(coverage_ready))
        coverage_final = load_json(segment_path(coverage, 0, "final"))
        coverage_final["ready_sha256"] = sha256_file(segment_path(coverage, 0, "ready"))
        segment_path(coverage, 0, "final").write_bytes(canonical_json(coverage_final))
        expect_reject(
            lambda: validate_thermal(coverage, fixture_design, fixture_tasks),
            "thermal start coverage",
        )

        sampler_exit = tampered_copy("sampler-exit")
        sampler_final = load_json(segment_path(sampler_exit, 0, "final"))
        sampler_final["sampler_returncode"] = 9
        segment_path(sampler_exit, 0, "final").write_bytes(canonical_json(sampler_final))
        expect_reject(
            lambda: validate_thermal(sampler_exit, fixture_design, fixture_tasks),
            "sampler exit",
        )

        hard_burst = Path(temporary) / "hard-burst"
        hard_burst.mkdir()
        hard_tasks, hard_design = make_success_fixture(hard_burst, busy="20.0")
        hard_receipt = validate_thermal(hard_burst, hard_design, hard_tasks)
        if hard_receipt["successful_segments"] != 1:
            die("hard burst segment below the control floor was not accepted")

        control_floor = Path(temporary) / "control-floor"
        control_floor.mkdir()
        control_tasks, control_design = make_success_fixture(
            control_floor, stage="control", busy=str(CPU_BUSY_FLOOR))
        control_receipt = validate_thermal(
            control_floor, control_design, control_tasks)
        if control_receipt["successful_segments"] != 1:
            die("control segment at the sealed floor was not accepted")

        control_low = Path(temporary) / "control-low"
        control_low.mkdir()
        low_tasks, low_design = make_success_fixture(
            control_low, stage="control", busy="89.9")
        expect_reject(
            lambda: validate_thermal(control_low, low_design, low_tasks),
            "control busy floor",
        )

        interrupted = Path(temporary) / "interrupted"
        for directory in (
                "raw", "stderr", "exit", "thermal", "segments",
                "attempts", "job_receipts"):
            (interrupted / directory).mkdir(parents=True, exist_ok=True)
        retry_task = {"job": 0, "stage": "hard", "output_name": "retry.csv", "argv": ["/bin/true"]}
        retry_intent = {
            "schema": SCHEMA + ".segment_intent", "segment": 0,
            "stage": "hard", "jobs": [0],
            "jobs_sha256": sha256_bytes(json_lines([0])), "workers": 1,
            "retry_policy": "stage-atomic non-selective retry",
        }
        write_once(segment_path(interrupted, 0, "intent"), canonical_json(retry_intent))
        (interrupted / "attempts" / "segment000").mkdir()
        (interrupted / "raw" / "retry.csv").write_bytes(b"failed output\n")
        (interrupted / "stderr" / "retry.csv.stderr").write_bytes(b"failure\n")
        (interrupted / "exit" / "retry.csv.exit").write_bytes(b"1\n")
        reconciled = reconcile_incomplete_segments(interrupted, [retry_task])
        if len(reconciled) != 1 or reconciled[0]["state"] != "interrupted" or \
                completed_jobs(interrupted, [retry_task]) or \
                any((interrupted / directory / name).exists() for directory, name in (
                    ("raw", "retry.csv"), ("stderr", "retry.csv.stderr"),
                    ("exit", "retry.csv.exit"))) or \
                reconcile_incomplete_segments(interrupted, [retry_task]):
            die("interrupted failed-triplet resume selftest failed")
        expect_reject(lambda: validate_thermal_rows([], segment=0), "header-only thermal")
    print("wh2_grouped_row_mask_campaign selftest: PASS")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prepare = sub.add_parser("prepare", help="freeze a new campaign; never launches jobs")
    prepare.add_argument("--output-root", required=True)
    prepare.add_argument("--allk-root", required=True)
    prepare.add_argument("--sampler", required=True)
    prepare.add_argument("--workers", type=int, default=DEFAULT_WORKERS, choices=range(1, 257))
    prepare.set_defaults(function=command_prepare)
    launch = sub.add_parser("launch", help="preflight or execute a frozen campaign")
    launch.add_argument("--root", required=True)
    launch.add_argument("--expected-receipts-sha256", required=True)
    mode = launch.add_mutually_exclusive_group(required=True)
    mode.add_argument("--preflight-only", action="store_true")
    mode.add_argument("--execute", action="store_true")
    launch.add_argument("--resume", action="store_true")
    launch.add_argument("--stage", choices=("hard", "control", "all"), default="all")
    launch.set_defaults(function=command_launch)
    reduce_parser = sub.add_parser("reduce", help="strictly validate and reduce all outputs")
    reduce_parser.add_argument("--root", required=True)
    reduce_parser.add_argument("--expected-receipts-sha256", required=True)
    reduce_parser.set_defaults(function=command_reduce)
    selftest = sub.add_parser("selftest", help="run bounded pure-Python invariants")
    selftest.set_defaults(function=command_selftest)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if getattr(args, "resume", False) and not getattr(args, "execute", False):
        parser.error("--resume requires --execute")
    try:
        args.function(args)
    except CampaignError as exc:
        print("campaign error: {}".format(exc), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
