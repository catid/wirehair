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
import itertools
import json
import os
import re
import shutil
import subprocess
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation, getcontext
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Set, Tuple


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
DEFAULT_WORKERS = 120
CPU_BUSY_FLOOR = Decimal("95")

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
    try:
        with path.open("xb") as stream:
            stream.write(data)
    except FileExistsError:
        if path.read_bytes() != data:
            die("refusing to replace immutable artifact {}".format(path))


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
    return {
        "controller": root / "frozen" / CONTROLLER_NAME,
        "binary": root / "frozen" / BINARY_NAME,
        "sampler": root / "frozen" / SAMPLER_NAME,
        "source_summary": root / "frozen" / SOURCE_SUMMARY_NAME,
        "design": root / "design.json",
        "hard_keys": root / "hard_keys.jsonl",
        "p48_controls": root / "controls_p48.jsonl",
        "p32_controls": root / "controls_p32.jsonl",
        "p48_masks": root / "masks_p48_r3.jsonl",
        "p32_masks": root / "masks_p32_r7.jsonl",
        "tasks": root / "tasks_manifest.jsonl",
    }


def command_prepare(args: argparse.Namespace) -> None:
    root = Path(args.output_root).resolve()
    source_path = Path(args.allk_root).resolve() / SOURCE_SUMMARY_NAME
    binary = Path(args.binary).resolve()
    sampler = Path(args.sampler).resolve()
    controller = Path(__file__).resolve()
    if root.exists():
        die("output root already exists: {}".format(root))
    if sha256_file(source_path) != SOURCE_SUMMARY_SHA256:
        die("all-K validated-summary hash mismatch")
    for path, label in ((binary, "binary"), (sampler, "sampler"), (controller, "controller")):
        if not path.is_file() or path.is_symlink():
            die("{} is not a regular file: {}".format(label, path))
    if not os.access(str(binary), os.X_OK):
        die("benchmark binary is not executable")
    git_head = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=str(controller.parent.parent), text=True,
    ).strip()
    ancestor = subprocess.run(
        ["git", "merge-base", "--is-ancestor", CODEC_IMPLEMENTATION_HEAD, git_head],
        cwd=str(controller.parent.parent), stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    if ancestor.returncode != 0:
        die("campaign source does not contain the reviewed row-mask implementation")

    root.mkdir(parents=True)
    for directory in ("frozen", "raw", "stderr", "exit", "thermal"):
        (root / directory).mkdir()
    shutil.copyfile(str(controller), str(root / "frozen" / CONTROLLER_NAME))
    shutil.copyfile(str(binary), str(root / "frozen" / BINARY_NAME))
    shutil.copyfile(str(sampler), str(root / "frozen" / SAMPLER_NAME))
    shutil.copyfile(str(source_path), str(root / "frozen" / SOURCE_SUMMARY_NAME))
    os.chmod(str(root / "frozen" / CONTROLLER_NAME), 0o555)
    os.chmod(str(root / "frozen" / BINARY_NAME), 0o555)
    os.chmod(str(root / "frozen" / SAMPLER_NAME), 0o555)

    summary = load_json(root / "frozen" / SOURCE_SUMMARY_NAME)
    frozen_binary = str(root / "frozen" / BINARY_NAME)
    plan = plan_from_summary(summary, frozen_binary)
    write_once(root / "hard_keys.jsonl", json_lines(plan["hard"]))
    for period in (48, 32):
        write_once(root / "controls_p{}.jsonl".format(period),
                   json_lines(plan["controls"][period]))
    write_once(root / "masks_p48_r3.jsonl", json_lines(plan["masks"][48]))
    write_once(root / "masks_p32_r7.jsonl", json_lines(plan["masks"][32]))
    write_once(root / "tasks_manifest.jsonl", json_lines(plan["tasks"]))

    design = {
        "schema": SCHEMA + ".design", "root": str(root),
        "git_head": git_head, "codec_implementation_head": CODEC_IMPLEMENTATION_HEAD,
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
        "thermal_cpu_busy_floor_pct": float(CPU_BUSY_FLOOR),
        "timing_promotional": False,
        "timing_policy": "NONPROMOTIONAL saturated rank-proxy screen; recovery and work receipts only",
    }
    write_once(root / "design.json", canonical_json(design))
    paths = artifact_map(root)
    receipts = {
        "schema": SCHEMA + ".prelaunch", "artifacts": {
            key: {"path": str(path.relative_to(root)), "sha256": sha256_file(path)}
            for key, path in sorted(paths.items())
        },
    }
    write_once(root / "prelaunch_receipts.json", canonical_json(receipts))
    print(json.dumps({
        "status": "PREPARED_NOT_LAUNCHED", "root": str(root),
        "prelaunch_receipts_sha256": sha256_file(root / "prelaunch_receipts.json"),
        "tasks": len(plan["tasks"]), "cells": design["total_cells"],
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


def verify_root(root: Path, expected_receipts: str) -> Tuple[dict, List[dict]]:
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
            design.get("source_summary_sha256") != SOURCE_SUMMARY_SHA256:
        die("design trust anchor mismatch")
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
        actual = {path.name for path in (root / directory).iterdir() if path.is_file()}
        expected = {name + suffix for name in expected_outputs}
        if not actual <= expected:
            die("unexpected {} output(s)".format(directory))
        if any(".part." in name for name in actual):
            die("partial {} output remains".format(directory))
    return design, tasks


def completed_jobs(root: Path, tasks: Sequence[dict]) -> Set[int]:
    completed: Set[int] = set()
    for task in tasks:
        name = task["output_name"]
        paths = (
            root / "raw" / name,
            root / "stderr" / (name + ".stderr"),
            root / "exit" / (name + ".exit"),
        )
        present = tuple(path.exists() for path in paths)
        if any(present) and not all(present):
            die("partial output triplet for job {}".format(task["job"]))
        if all(present):
            if not paths[0].stat().st_size or paths[1].stat().st_size or \
                    paths[2].read_text(encoding="ascii").strip() != "0":
                die("unclean completed job {}".format(task["job"]))
            completed.add(int(task["job"]))
    return completed


def other_samplers() -> List[dict]:
    found: List[dict] = []
    for proc in Path("/proc").iterdir():
        if not proc.name.isdigit() or int(proc.name) == os.getpid():
            continue
        try:
            tokens = [token.decode("utf-8", "replace") for token in
                      (proc / "cmdline").read_bytes().split(b"\0") if token]
        except OSError:
            continue
        if any(Path(token).name == SAMPLER_NAME for token in tokens):
            found.append({"pid": int(proc.name), "command": " ".join(tokens)})
    return sorted(found, key=lambda item: item["pid"])


def atomic_result(path: Path, data: bytes) -> None:
    temporary = Path(str(path) + ".part.{}".format(os.getpid()))
    with temporary.open("xb") as stream:
        stream.write(data)
    os.replace(str(temporary), str(path))


def command_launch(args: argparse.Namespace) -> None:
    root = Path(args.root).resolve()
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
    if args.preflight_only:
        print(json.dumps({
            "status": "READY_NOT_LAUNCHED", "root": str(root),
            "tasks": len(tasks), "complete": len(complete),
            "remaining": len(remaining), "workers": design["worker_count"],
            "remaining_by_stage": {
                stage: len(values) for stage, values in by_stage.items()
            },
            "requested_stage": args.stage, "single_sampler": True,
        }, indent=2, sort_keys=True))
        return
    if complete and not args.resume:
        die("outputs already exist; inspect then pass --resume")
    if not remaining:
        die("campaign ledger is already complete")
    hard_incomplete = bool(by_stage["hard"])
    if args.stage == "control" and hard_incomplete:
        die("control stage requires the complete hard-stage Cartesian product")
    requested_stages = ("hard", "control") if args.stage == "all" else (args.stage,)
    requested = [task for stage in requested_stages for task in by_stage[stage]]
    if not requested:
        die("requested campaign stage is already complete")

    segment = len(list((root / "thermal").glob("segment*.csv")))
    thermal_csv = root / "thermal" / "segment{:03d}.csv".format(segment)
    thermal_pid = root / "thermal" / "segment{:03d}.pid".format(segment)
    thermal_stderr = root / "thermal" / "segment{:03d}.stderr".format(segment)
    sampler_error = thermal_stderr.open("xb")
    sampler = subprocess.Popen([
        "sudo", "-n", sys.executable, str(root / "frozen" / SAMPLER_NAME),
        "--csv", str(thermal_csv), "--pid-file", str(thermal_pid),
        "--interval", "1", "--dimm-attempts", "5", "--dimm-retry-delay", "0.01",
    ], stdout=subprocess.DEVNULL, stderr=sampler_error)
    deadline = time.monotonic() + 8.0
    while time.monotonic() < deadline:
        if sampler.poll() is not None:
            sampler_error.close()
            die("thermal sampler exited early")
        if thermal_csv.exists() and thermal_csv.stat().st_size and thermal_pid.exists():
            break
        time.sleep(0.05)
    else:
        subprocess.run(["sudo", "-n", "kill", "-TERM", str(sampler.pid)], check=False)
        sampler.wait(timeout=10)
        sampler_error.close()
        die("thermal sampler did not become ready")

    started = datetime.now(timezone.utc).isoformat()
    failures: List[dict] = []

    def run(task: dict) -> dict:
        result = subprocess.run(task["argv"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        name = task["output_name"]
        atomic_result(root / "raw" / name, result.stdout)
        atomic_result(root / "stderr" / (name + ".stderr"), result.stderr)
        atomic_result(root / "exit" / (name + ".exit"),
                      (str(result.returncode) + "\n").encode("ascii"))
        return {"job": task["job"], "returncode": result.returncode}

    try:
        for stage in requested_stages:
            stage_tasks = by_stage[stage]
            if not stage_tasks:
                continue
            with concurrent.futures.ThreadPoolExecutor(
                    max_workers=int(design["worker_count"])) as pool:
                for result in pool.map(run, stage_tasks):
                    if result["returncode"]:
                        failures.append(result)
            if failures:
                break
    finally:
        try:
            sampled_pid = int(thermal_pid.read_text(encoding="ascii").strip())
            subprocess.run(["sudo", "-n", "kill", "-TERM", str(sampled_pid)], check=False)
        except (OSError, ValueError):
            sampled_pid = None
        try:
            sampler.wait(timeout=15)
        except subprocess.TimeoutExpired:
            subprocess.run(["sudo", "-n", "kill", "-KILL", str(sampler.pid)], check=False)
            sampler.wait(timeout=5)
        sampler_error.close()
    receipt = {
        "schema": SCHEMA + ".launch", "segment": segment,
        "started_utc": started, "ended_utc": datetime.now(timezone.utc).isoformat(),
        "resume": bool(args.resume), "workers": design["worker_count"],
        "requested_stage": args.stage, "stage_order": list(requested_stages),
        "jobs_requested": len(requested),
        "jobs_requested_by_stage": {
            stage: len(by_stage[stage]) for stage in requested_stages
        },
        "jobs_sha256": sha256_bytes(json_lines([task["job"] for task in requested])),
        "previously_complete": len(complete), "failures": failures,
        "thermal_csv": thermal_csv.name, "thermal_stderr": thermal_stderr.name,
        "single_sampler": True,
    }
    write_once(root / "launch_receipt.segment{:03d}.json".format(segment),
               canonical_json(receipt))
    print(json.dumps(receipt, indent=2, sort_keys=True))
    if failures:
        raise SystemExit(1)


def decimal_field(row: Mapping[str, str], key: str) -> Decimal:
    try:
        value = Decimal(row[key])
    except (KeyError, InvalidOperation) as exc:
        die("invalid decimal field {}: {}".format(key, exc))
    if not value.is_finite():
        die("nonfinite decimal field {}".format(key))
    return value


def classify(row: Mapping[str, str]) -> str:
    success, rank_fail, error = (int(row[key]) for key in ("success", "rank_fail", "error"))
    if any(value not in (0, 1) for value in (success, rank_fail, error)) or \
            success + rank_fail + error != 1 or error:
        die("invalid benchmark outcome")
    if success:
        return "success"
    return "q>H" if int(row["binary_def_max"]) > 12 else "field_shortfall"


def validate_thermal(root: Path, design: Mapping[str, object]) -> Dict[str, object]:
    segments = sorted((root / "thermal").glob("segment*.csv"))
    receipts = sorted(root.glob("launch_receipt.segment*.json"))
    if not segments or len(segments) != len(receipts):
        die("thermal/launch segment cardinality mismatch")
    busy: List[Decimal] = []
    cpu: List[Decimal] = []
    dimm: List[Decimal] = []
    samples = 0
    for index, (path, receipt_path) in enumerate(zip(segments, receipts)):
        receipt = load_json(receipt_path)
        if receipt.get("schema") != SCHEMA + ".launch" or \
                receipt.get("segment") != index or receipt.get("single_sampler") is not True or \
                receipt.get("workers") != design.get("worker_count") or receipt.get("failures") != [] or \
                receipt.get("thermal_csv") != path.name or \
                receipt.get("thermal_stderr") != path.with_suffix(".stderr").name:
            die("launch/thermal receipt mismatch for segment {}".format(index))
        stderr = path.with_suffix(".stderr")
        if not stderr.exists() or stderr.stat().st_size:
            die("missing/nonempty thermal sampler stderr")
        with path.open(newline="", encoding="ascii") as stream:
            reader = csv.DictReader(stream)
            if tuple(reader.fieldnames or ()) != THERMAL_FIELDS:
                die("thermal header mismatch")
            for row in reader:
                samples += 1
                if row["cpu_busy_pct"]:
                    busy.append(decimal_field(row, "cpu_busy_pct"))
                if row["cpu_tctl_c"]:
                    cpu.append(decimal_field(row, "cpu_tctl_c"))
                for key in THERMAL_FIELDS[5:13]:
                    if row[key]:
                        dimm.append(decimal_field(row, key))
                if int(row["dimm_read_errors"]) != 0 or \
                        int(row["edac_ce"]) != 0 or int(row["edac_ue"]) != 0:
                    die("thermal hardware error receipt is nonzero")
    if samples < 2 or not busy or sum(busy) / len(busy) < CPU_BUSY_FLOOR:
        die("thermal load receipt is incomplete or below CPU busy floor")
    return {
        "segments": len(segments), "samples": samples,
        "busy_mean": str(sum(busy) / len(busy)),
        "cpu_max_c": str(max(cpu)) if cpu else None,
        "dimm_max_c": str(max(dimm)) if dimm else None,
    }


def command_reduce(args: argparse.Namespace) -> None:
    root = Path(args.root).resolve()
    design, tasks = verify_root(root, args.expected_receipts_sha256)
    if len(completed_jobs(root, tasks)) != len(tasks):
        die("campaign is incomplete")
    thermal = validate_thermal(root, design)
    outcomes: Dict[Tuple[int, int, str, int, int, str], dict] = {}
    manifest_lines = bytearray()
    for task in tasks:
        name = task["output_name"]
        raw_path = root / "raw" / name
        stderr_path = root / "stderr" / (name + ".stderr")
        exit_path = root / "exit" / (name + ".exit")
        for path in (raw_path, stderr_path, exit_path):
            manifest_lines.extend("{}  {}\n".format(
                sha256_file(path), path.relative_to(root)).encode("ascii"))
        with raw_path.open(newline="", encoding="ascii") as stream:
            preamble = parse_preamble(stream.readline().rstrip("\n"))
            expected = expected_preamble(task)
            if any(preamble.get(key) != value for key, value in expected.items()):
                die("preamble mismatch job {}".format(task["job"]))
            reader = csv.DictReader(stream)
            if tuple(reader.fieldnames or ()) != CSV_FIELDS:
                die("CSV header mismatch job {}".format(task["job"]))
            rows = list(reader)
        if [int(row["N"]) for row in rows] != task["Ks"]:
            die("CSV K ledger mismatch job {}".format(task["job"]))
        for row in rows:
            fixed = {
                "bb": "64", "heavy_family": "periodic", "mix_count": "2",
                "overhead": "0", "trials": "1", "active_packet_peel_seed_xor": "0x0",
            }
            if any(row.get(key) != value for key, value in fixed.items()):
                die("fixed CSV field mismatch job {}".format(task["job"]))
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
    print("wh2_grouped_row_mask_campaign selftest: PASS")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prepare = sub.add_parser("prepare", help="freeze a new campaign; never launches jobs")
    prepare.add_argument("--output-root", required=True)
    prepare.add_argument("--allk-root", required=True)
    prepare.add_argument("--binary", required=True)
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
