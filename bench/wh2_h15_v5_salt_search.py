#!/usr/bin/env python3
"""Prepare and audit bounded WH2 normalized-H15-v5 peel-salt experiments.

The discovery set is intentionally limited to independently audited baseline
failures from the frozen gen2b campaign.  Candidate-breaker failures are never
used as peel-salt training input.  This tool writes only to a caller-provided
fresh experiment directory; the frozen campaign remains read-only.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import tempfile
from collections import defaultdict
from decimal import Decimal
from pathlib import Path


BASE_COMMIT = "709a3c5060a7434ff2652e64528cd869bcfd874d"
BASE_BINARY_SHA256 = (
    "4340fd9444b6dc2adbd94304b55e987d9887efd8f6d86f205511e7fda75cb41c"
)
SOURCE_ARTIFACTS_SHA256 = (
    "e40272946f821e740cd6666bb611427e1c56f3d0ce0996123ad540c46bf7c171"
)
ANALYSIS_COMPLETE_SHA256 = (
    "49d50aacb2722488abc0122bae0d15a0aa4e38daad7bd6f033871ff809db02aa"
)
CAMPAIGN_SCHEMA = "wirehair.wh2.gfbreaker_exhaustive_k2_64000.v2"
DISCOVERY_SEEDS = (
    "0xa6b527b9ae8de8d7",
    "0x75446e1619e81d8a",
    "0x007e9dd892a319f5",
)
SCHEDULES = ("burst", "adversarial", "repair-only")
LOSS = "0.50"
OUTPUT_LOSS = "0.5"
EXPECTED_FAILURE_CELLS = 477
EXPECTED_FAILURE_K = 466

RAW_HEADER = (
    "N", "bb", "heavy_family", "mix_count", "overhead", "trials",
    "success", "rank_fail", "error", "fail_rate", "inact_mu",
    "inact_max", "binary_def_mu", "binary_def_max", "heavy_gain_mu",
    "heavy_gain_min", "heavy_shortfall", "solve_ms_mu", "build_ms_mu",
    "peel_ms_mu", "project_ms_mu", "residual_ms_mu", "backsub_ms_mu",
    "seed_attempt", "block_xors_mu", "block_muladds_mu",
    "first_rank_fail", "binary_def_hist", "heavy_gain_hist",
    "failure_trials", "active_packet_peel_seed_xor",
)

FAILURE_HEADER = (
    "witness_id", "K", "bb", "certified_v4_salt", "schedule", "loss",
    "seed_index", "seed", "trial", "rank_fail", "error",
    "failure_trials", "heavy_shortfall", "inact_mu", "inact_max",
    "binary_def_mu", "binary_def_max", "seed_attempt", "block_xors_mu",
    "block_muladds_mu", "source_base", "source_sha256",
)

EXACT_JOB_HEADER = (
    "job_id", "stem", "schedule", "seed_index", "seed", "loss", "salt",
    "k_count", "k_csv",
)


class AuditError(RuntimeError):
    pass


def fail(message: str) -> None:
    raise AuditError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_kv(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for number, line in enumerate(path.read_text().splitlines(), 1):
        if not line or "=" not in line:
            fail(f"{path}:{number}: malformed key/value line")
        key, value = line.split("=", 1)
        if not key or key in result:
            fail(f"{path}:{number}: duplicate or empty key")
        result[key] = value
    return result


def read_source_artifact_seal(
    source_seal: Path,
    analysis_complete: Path,
) -> dict[str, str]:
    if sha256(source_seal) != SOURCE_ARTIFACTS_SHA256:
        fail(f"{source_seal}: frozen source-artifact seal hash mismatch")
    if sha256(analysis_complete) != ANALYSIS_COMPLETE_SHA256:
        fail(f"{analysis_complete}: frozen analysis-complete hash mismatch")
    complete_lines = analysis_complete.read_text().splitlines()
    expected_line = f"{SOURCE_ARTIFACTS_SHA256}  source_artifacts.sha256"
    if expected_line not in complete_lines:
        fail(f"{analysis_complete}: source-artifact seal is not linked")

    result: dict[str, str] = {}
    for number, line in enumerate(source_seal.read_text().splitlines(), 1):
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
        if not match:
            fail(f"{source_seal}:{number}: malformed SHA-256 record")
        digest, relative_text = match.groups()
        relative = Path(relative_text)
        if relative.is_absolute() or ".." in relative.parts or relative_text in result:
            fail(f"{source_seal}:{number}: unsafe or duplicate relative path")
        result[relative_text] = digest
    return result


def parse_preamble(line: str, path: Path) -> dict[str, str]:
    prefix = "# precodefail: "
    if not line.startswith(prefix):
        fail(f"{path}: missing precodefail preamble")
    result: dict[str, str] = {}
    for token in line[len(prefix):].split():
        if "=" not in token:
            fail(f"{path}: malformed preamble token {token!r}")
        key, value = token.split("=", 1)
        if not key or key in result:
            fail(f"{path}: duplicate or empty preamble key")
        result[key] = value
    return result


def strict_uint(text: str, context: str) -> int:
    if not re.fullmatch(r"(?:0|[1-9][0-9]*)", text):
        fail(f"{context}: noncanonical decimal unsigned integer {text!r}")
    return int(text)


def parse_u64(text: str, context: str) -> int:
    if not re.fullmatch(r"(?:0[xX][0-9a-fA-F]+|[0-9]+)", text):
        fail(f"{context}: invalid unsigned integer {text!r}")
    value = int(text, 0 if text.lower().startswith("0x") else 10)
    if value < 0 or value > (1 << 64) - 1:
        fail(f"{context}: out of u64 range")
    return value


def canonical_salt(text: str, context: str) -> str:
    value = parse_u64(text, context)
    if value > 255:
        fail(f"{context}: salt exceeds uint8 range")
    return f"0x{value:x}"


def nonnegative_decimal(text: str, context: str) -> Decimal:
    try:
        value = Decimal(text)
    except Exception as exc:
        fail(f"{context}: invalid decimal {text!r}: {exc}")
    if not value.is_finite() or value < 0:
        fail(f"{context}: nonfinite or negative decimal {text!r}")
    return value


def require_fresh_dir(path: Path) -> None:
    if not path.is_dir():
        fail(f"output directory does not exist: {path}")
    if any(path.iterdir()):
        fail(f"output directory must be empty: {path}")


def write_csv(
    path: Path,
    header: tuple[str, ...],
    rows: list[dict[str, object]],
    *,
    delimiter: str = ",",
) -> None:
    with path.open("x", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=header,
            delimiter=delimiter,
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_requested_map(path: Path) -> dict[int, str]:
    result: dict[int, str] = {}
    with path.open(newline="") as stream:
        for number, fields in enumerate(csv.reader(stream, delimiter="\t"), 1):
            if len(fields) != 2:
                fail(f"{path}:{number}: expected two tab-separated fields")
            k = strict_uint(fields[0], f"{path}:{number}:K")
            salt = canonical_salt(fields[1], f"{path}:{number}")
            if k in result:
                fail(f"{path}:{number}: duplicate K={k}")
            result[k] = salt
    if set(result) != set(range(2, 64001)):
        fail(f"{path}: expected exact K range 2..64000")
    return result


def read_groups(path: Path) -> dict[int, tuple[str, list[int]]]:
    groups: dict[int, tuple[str, list[int]]] = {}
    with path.open(newline="") as stream:
        for number, fields in enumerate(csv.reader(stream, delimiter="\t"), 1):
            if len(fields) != 3:
                fail(f"{path}:{number}: expected three tab-separated fields")
            group_id = strict_uint(fields[0], f"{path}:{number}:group_id")
            salt = canonical_salt(fields[1], f"{path}:{number}:salt")
            ks = [strict_uint(value, f"{path}:{number}:K") for value in fields[2].split(",")]
            if not ks or ks != sorted(set(ks)):
                fail(f"{path}:{number}: K list is empty, duplicate, or unsorted")
            if group_id in groups:
                fail(f"{path}:{number}: duplicate group {group_id}")
            groups[group_id] = (salt, ks)
    if set(groups) != set(range(186)):
        fail(f"{path}: expected exact group IDs 0..185")
    flattened = [k for group_id in range(186) for k in groups[group_id][1]]
    if len(flattened) != 63999 or sorted(flattened) != list(range(2, 64001)):
        fail(f"{path}: groups do not form a disjoint exact K partition 2..64000")
    return groups


def read_campaign_jobs(path: Path) -> list[dict[str, object]]:
    jobs: list[dict[str, object]] = []
    seen: set[tuple[int, str, int]] = set()
    with path.open(newline="") as stream:
        for number, fields in enumerate(csv.reader(stream, delimiter="\t"), 1):
            if len(fields) != 6:
                fail(f"{path}:{number}: expected six tab-separated fields")
            job_id = strict_uint(fields[0], f"{path}:{number}:job_id")
            group_id = strict_uint(fields[1], f"{path}:{number}:group_id")
            schedule = fields[2]
            loss = fields[3]
            seed_index = strict_uint(fields[4], f"{path}:{number}:seed_index")
            seed = fields[5]
            if job_id != len(jobs):
                fail(f"{path}:{number}: noncontiguous job ID")
            if group_id not in range(186):
                fail(f"{path}:{number}: group ID outside 0..185")
            if schedule not in SCHEDULES or loss != LOSS:
                fail(f"{path}:{number}: schedule/loss mismatch")
            if seed_index not in range(len(DISCOVERY_SEEDS)) or seed != DISCOVERY_SEEDS[seed_index]:
                fail(f"{path}:{number}: seed mismatch")
            key = (group_id, schedule, seed_index)
            if key in seen:
                fail(f"{path}:{number}: duplicate group/schedule/seed cell")
            seen.add(key)
            jobs.append({
                "job_id": job_id,
                "group_id": group_id,
                "schedule": schedule,
                "loss": loss,
                "seed_index": seed_index,
                "seed": seed,
            })
    expected = {
        (group_id, schedule, seed_index)
        for schedule in SCHEDULES
        for seed_index in range(len(DISCOVERY_SEEDS))
        for group_id in range(186)
    }
    if len(jobs) != 1674 or seen != expected:
        fail(f"{path}: expected complete 186x3x3 baseline job grid")
    return jobs


def load_raw_base(path: Path) -> tuple[dict[str, str], dict[int, dict[str, str]]]:
    with path.open(newline="") as stream:
        preamble_line = stream.readline().rstrip("\r\n")
        preamble = parse_preamble(preamble_line, path)
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != RAW_HEADER:
            fail(f"{path}: unexpected CSV header")
        rows: dict[int, dict[str, str]] = {}
        for number, row in enumerate(reader, 3):
            if None in row:
                fail(f"{path}:{number}: malformed CSV row")
            k = strict_uint(row["N"], f"{path}:{number}:N")
            if k in rows:
                fail(f"{path}:{number}: duplicate K={k}")
            rows[k] = row
    return preamble, rows


def validate_raw_geometry(
    path: Path,
    preamble: dict[str, str],
    schedule: str,
    seed: str,
    salt: str,
) -> None:
    expected = {
        "trials": "1",
        "threads": "1",
        "loss": OUTPUT_LOSS,
        "completion": "mixed",
        "mixed_period": "32",
        "mixed_gf256_rows": "11",
        "mixed_gf16_rows": "4",
        "mixed_geometry": "shared-x",
        "mixed_residue_skew": "0",
        "mixed_residue_schedule": "hashed",
        "mixed_residue_hash_seed": "0x44",
        "mixed_residue_hash_keyed": "1",
        "mixed_independent_extension_residues": "1",
        "mixed_extension_residue_seed_xor": "0x4e",
        "source_hits_override": "0",
        "packet_peel_seed_table": "none",
        "binary_dense_rows_override": "0",
        "gf256_heavy_rows_override": "0",
        "odd_packet_peel_seed_xor": "0x0",
        "packet_row_seed_multiplier": "0x1",
        "packet_row_seed_avalanche": "0",
        "seed_block_bytes_override": "1280",
        "overhead_stream": "salted",
        "full_payload_solve": "0",
        "schedule": schedule,
    }
    for key, value in expected.items():
        if preamble.get(key) != value:
            fail(f"{path}: {key}={preamble.get(key)!r}, expected {value!r}")
    if parse_u64(preamble.get("seed", ""), f"{path}:seed") != parse_u64(seed, "cell seed"):
        fail(f"{path}: seed mismatch")
    if canonical_salt(preamble.get("packet_peel_seed_xor", ""), f"{path}:salt") != salt:
        fail(f"{path}: packet peel salt mismatch")
    if "mixed_independent_gf256_breaker_residues" in preamble:
        fail(f"{path}: candidate-only breaker flag present in baseline output")


def extract(args: argparse.Namespace) -> None:
    campaign = args.campaign.resolve()
    source_seal_path = args.source_seal.resolve()
    analysis_complete_path = args.analysis_complete.resolve()
    out = args.output.resolve()
    require_fresh_dir(out)

    sealed_hashes = read_source_artifact_seal(source_seal_path, analysis_complete_path)

    def verify_frozen_source(relative_text: str) -> tuple[Path, str]:
        relative = Path(relative_text)
        if relative.is_absolute() or ".." in relative.parts:
            fail(f"unsafe frozen source path: {relative_text}")
        path = (campaign / relative).resolve()
        try:
            path.relative_to(campaign)
        except ValueError:
            fail(f"frozen source escaped campaign: {relative_text}")
        expected = sealed_hashes.get(relative_text)
        if expected is None:
            fail(f"frozen source is absent from pinned seal: {relative_text}")
        actual = sha256(path)
        if actual != expected:
            fail(f"frozen source hash mismatch: {relative_text}")
        return path, actual

    for relative in (
        "manifest.txt",
        "meta/requested_k_salt.tsv",
        "meta/groups.tsv",
        "meta/jobs.tsv",
    ):
        verify_frozen_source(relative)

    manifest = read_kv(campaign / "manifest.txt")
    if manifest.get("schema") != CAMPAIGN_SCHEMA:
        fail("campaign schema mismatch")
    if manifest.get("base_commit") != BASE_COMMIT:
        fail("campaign baseline commit mismatch")
    if manifest.get("baseline_sha256") != BASE_BINARY_SHA256:
        fail("campaign baseline binary hash mismatch")
    if manifest.get("seeds", "").split(",") != list(DISCOVERY_SEEDS):
        fail("campaign discovery seeds mismatch")
    if manifest.get("schedules", "").split(",") != list(SCHEDULES):
        fail("campaign schedules mismatch")
    if manifest.get("losses") != LOSS or manifest.get("trials_per_cell") != "1":
        fail("campaign loss/trial contract mismatch")

    requested_map = read_requested_map(campaign / "meta/requested_k_salt.tsv")
    groups = read_groups(campaign / "meta/groups.tsv")
    campaign_jobs = read_campaign_jobs(campaign / "meta/jobs.tsv")
    raw_hashes: dict[Path, str] = {}
    failures: list[dict[str, object]] = []

    raw_root = (campaign / "raw/stdout").resolve()
    for campaign_job in campaign_jobs:
        job_id = int(campaign_job["job_id"])
        group_id = int(campaign_job["group_id"])
        schedule = str(campaign_job["schedule"])
        seed_index = int(campaign_job["seed_index"])
        seed = str(campaign_job["seed"])
        salt, expected_k = groups[group_id]
        stem = f"job{job_id:04d}_seed{seed_index}_{schedule}_l050_group{group_id:03d}.base.csv"
        source_rel = Path("raw/stdout") / stem
        if source_rel.is_absolute() or ".." in source_rel.parts:
            fail(f"unsafe baseline source path: {source_rel}")
        source = (campaign / source_rel).resolve()
        if source.parent != raw_root:
            fail(f"baseline source escaped raw/stdout: {source}")
        verified_source, verified_hash = verify_frozen_source(str(source_rel))
        if verified_source != source:
            fail(f"baseline source resolution mismatch: {source_rel}")
        preamble, raw_rows = load_raw_base(source)
        raw_hashes[source] = verified_hash
        validate_raw_geometry(source, preamble, schedule, seed, salt)
        if list(raw_rows) != expected_k:
            fail(f"{source}: K row order/content mismatch")

        for k, raw in raw_rows.items():
            if requested_map.get(k) != salt:
                fail(f"{source}: K={k} certified v4 salt mismatch")
            integer_fields = (
                "bb", "mix_count", "overhead", "trials", "success",
                "rank_fail", "error", "inact_max", "binary_def_max",
                "heavy_shortfall", "seed_attempt",
            )
            parsed = {key: strict_uint(raw[key], f"{source}:K={k}:{key}") for key in integer_fields}
            if parsed["bb"] != 64 or parsed["mix_count"] != 2 or parsed["overhead"] != 0:
                fail(f"{source}: K={k} row geometry mismatch")
            if parsed["trials"] != 1:
                fail(f"{source}: K={k} trial count mismatch")
            if any(parsed[key] not in (0, 1) for key in ("success", "rank_fail", "error")):
                fail(f"{source}: K={k} non-Boolean result count")
            if parsed["success"] + parsed["rank_fail"] + parsed["error"] != 1:
                fail(f"{source}: K={k} trial ledger mismatch")
            failure_ids = parse_failure_ids(raw["failure_trials"], 1, f"{source}:K={k}")
            if len(failure_ids) != parsed["rank_fail"] + parsed["error"]:
                fail(f"{source}: K={k} failure list/count mismatch")
            if canonical_salt(raw["active_packet_peel_seed_xor"], f"{source}:K={k}") != salt:
                fail(f"{source}: K={k} active salt mismatch")
            if parsed["rank_fail"] + parsed["error"] == 0:
                continue

            failures.append({
                "witness_id": "",
                "K": k,
                "bb": raw["bb"],
                "certified_v4_salt": salt,
                "schedule": schedule,
                "loss": LOSS,
                "seed_index": seed_index,
                "seed": seed,
                "trial": 0,
                "rank_fail": raw["rank_fail"],
                "error": raw["error"],
                "failure_trials": raw["failure_trials"],
                "heavy_shortfall": raw["heavy_shortfall"],
                "inact_mu": raw["inact_mu"],
                "inact_max": raw["inact_max"],
                "binary_def_mu": raw["binary_def_mu"],
                "binary_def_max": raw["binary_def_max"],
                "seed_attempt": raw["seed_attempt"],
                "block_xors_mu": raw["block_xors_mu"],
                "block_muladds_mu": raw["block_muladds_mu"],
                "source_base": str(source_rel),
                "source_sha256": raw_hashes[source],
            })

    schedule_index = {value: index for index, value in enumerate(SCHEDULES)}
    failures.sort(key=lambda row: (
        schedule_index[str(row["schedule"])], int(row["seed_index"]), int(row["K"])
    ))
    for index, row in enumerate(failures):
        row["witness_id"] = f"w{index:04d}"
    unique_k = sorted({int(row["K"]) for row in failures})
    if len(failures) != EXPECTED_FAILURE_CELLS or len(unique_k) != EXPECTED_FAILURE_K:
        fail(
            f"baseline failure support mismatch: cells={len(failures)} K={len(unique_k)}"
        )

    write_csv(out / "discovery_failures.csv", FAILURE_HEADER, failures)

    by_k: dict[int, list[dict[str, object]]] = defaultdict(list)
    for row in failures:
        by_k[int(row["K"])].append(row)
    k_rows: list[dict[str, object]] = []
    for k in unique_k:
        rows = by_k[k]
        salts = {str(row["certified_v4_salt"]) for row in rows}
        if len(salts) != 1:
            fail(f"K={k}: inconsistent certified salt across witnesses")
        k_rows.append({
            "K": k,
            "certified_v4_salt": next(iter(salts)),
            "witness_count": len(rows),
            "witness_ids": "|".join(str(row["witness_id"]) for row in rows),
        })
    write_csv(
        out / "discovery_k.csv",
        ("K", "certified_v4_salt", "witness_count", "witness_ids"),
        k_rows,
    )

    strata: dict[tuple[str, int, str], list[int]] = defaultdict(list)
    for row in failures:
        strata[(str(row["schedule"]), int(row["seed_index"]), str(row["seed"]))].append(int(row["K"]))
    jobs: list[dict[str, object]] = []
    for schedule in SCHEDULES:
        for seed_index, seed in enumerate(DISCOVERY_SEEDS):
            ks = sorted(strata[(schedule, seed_index, seed)])
            if not ks:
                fail(f"empty discovery stratum: {schedule}/seed{seed_index}")
            for salt in range(256):
                job_id = len(jobs)
                stem = f"job{job_id:04d}_{schedule}_seed{seed_index}_salt{salt:03d}"
                jobs.append({
                    "job_id": job_id,
                    "stem": stem,
                    "schedule": schedule,
                    "seed_index": seed_index,
                    "seed": seed,
                    "loss": LOSS,
                    "salt": salt,
                    "k_count": len(ks),
                    "k_csv": ",".join(str(k) for k in ks),
                })
    if len(jobs) != len(SCHEDULES) * len(DISCOVERY_SEEDS) * 256:
        fail("exact job count mismatch")
    write_csv(out / "exact_jobs.tsv", EXACT_JOB_HEADER, jobs, delimiter="\t")

    with (out / "source_files.sha256").open("x") as stream:
        source_paths = sorted(raw_hashes)
        extra = [
            campaign / "manifest.txt",
            campaign / "meta/requested_k_salt.tsv",
            campaign / "meta/groups.tsv",
            campaign / "meta/jobs.tsv",
            source_seal_path,
            analysis_complete_path,
        ]
        for path in extra + source_paths:
            stream.write(f"{sha256(path)}  {path}\n")

    nonzero_cells = sum(int(str(row["certified_v4_salt"]), 16) != 0 for row in failures)
    nonzero_k = sum(int(str(row["certified_v4_salt"]), 16) != 0 for row in k_rows)
    summary = {
        "schema": "wirehair.wh2.normalized_h15_v5.discovery.v1",
        "base_commit": BASE_COMMIT,
        "base_binary_sha256": BASE_BINARY_SHA256,
        "campaign": str(campaign),
        "source_artifacts_sha256": SOURCE_ARTIFACTS_SHA256,
        "analysis_complete_sha256": ANALYSIS_COMPLETE_SHA256,
        "discovery_policy": "frozen gen2b raw baseline outputs only; candidate files never opened",
        "failure_cells": len(failures),
        "failure_k": len(unique_k),
        "nonzero_v4_salt_failure_cells": nonzero_cells,
        "nonzero_v4_salt_failure_k": nonzero_k,
        "exact_jobs": len(jobs),
        "exact_logical_cells": sum(int(row["k_count"]) for row in jobs),
        "salt_domain": "all uint8 values 0..255",
        "discovery_seeds": list(DISCOVERY_SEEDS),
        "schedules": list(SCHEDULES),
        "loss": LOSS,
        "trials_per_exact_cell": 1,
    }
    with (out / "discovery_summary.json").open("x") as stream:
        json.dump(summary, stream, indent=2, sort_keys=True)
        stream.write("\n")


def read_tsv_jobs(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if tuple(reader.fieldnames or ()) != EXACT_JOB_HEADER:
            fail(f"{path}: exact job header mismatch")
        jobs = list(reader)
    seen: set[tuple[str, int, int]] = set()
    stratum_k: dict[tuple[str, int], tuple[int, ...]] = {}
    for row_number, row in enumerate(jobs, 2):
        if None in row or any(value is None for value in row.values()):
            fail(f"{path}:{row_number}: malformed or extra TSV fields")
        job_id = strict_uint(row["job_id"], f"{path}:{row_number}:job_id")
        seed_index = strict_uint(row["seed_index"], f"{path}:{row_number}:seed_index")
        salt = strict_uint(row["salt"], f"{path}:{row_number}:salt")
        k_count = strict_uint(row["k_count"], f"{path}:{row_number}:k_count")
        schedule = row["schedule"]
        if job_id != row_number - 2:
            fail(f"{path}:{row_number}: noncontiguous job ID")
        if schedule not in SCHEDULES or seed_index not in range(len(DISCOVERY_SEEDS)):
            fail(f"{path}:{row_number}: schedule/seed index mismatch")
        if row["seed"] != DISCOVERY_SEEDS[seed_index] or row["loss"] != LOSS:
            fail(f"{path}:{row_number}: seed/loss mismatch")
        if salt not in range(256):
            fail(f"{path}:{row_number}: salt outside uint8 range")
        expected_stem = f"job{job_id:04d}_{schedule}_seed{seed_index}_salt{salt:03d}"
        if row["stem"] != expected_stem:
            fail(f"{path}:{row_number}: stem mismatch")
        ks = tuple(strict_uint(value, f"{path}:{row_number}:K") for value in row["k_csv"].split(","))
        if not ks or ks != tuple(sorted(set(ks))) or len(ks) != k_count:
            fail(f"{path}:{row_number}: K list/count mismatch")
        stratum = (schedule, seed_index)
        if stratum in stratum_k and stratum_k[stratum] != ks:
            fail(f"{path}:{row_number}: inconsistent K coverage within stratum")
        stratum_k[stratum] = ks
        key = (schedule, seed_index, salt)
        if key in seen:
            fail(f"{path}:{row_number}: duplicate schedule/seed/salt job")
        seen.add(key)
    expected = {
        (schedule, seed_index, salt)
        for schedule in SCHEDULES
        for seed_index in range(len(DISCOVERY_SEEDS))
        for salt in range(256)
    }
    if len(jobs) != 2304 or seen != expected or len(stratum_k) != 9:
        fail(f"{path}: incomplete exact discovery grid")
    return jobs


def parse_failure_ids(text: str, trials: int, context: str) -> tuple[int, ...]:
    if not text:
        return ()
    values = tuple(strict_uint(part, context) for part in text.split("|"))
    if any(value < 0 or value >= trials for value in values):
        fail(f"{context}: failure trial outside range")
    if any(left >= right for left, right in zip(values, values[1:])):
        fail(f"{context}: failure trials are duplicate or unsorted")
    return values


def validate_exact_output(path: Path, job: dict[str, str]) -> dict[int, dict[str, str]]:
    preamble, rows = load_raw_base(path)
    salt = canonical_salt(job["salt"], f"{path}:job salt")
    validate_raw_geometry(path, preamble, job["schedule"], job["seed"], salt)
    expected_k = [int(value) for value in job["k_csv"].split(",")]
    if list(rows) != expected_k:
        fail(f"{path}: K row order/content mismatch")
    if len(rows) != int(job["k_count"]):
        fail(f"{path}: row count mismatch")
    for k, row in rows.items():
        checks = {
            "bb": "64",
            "heavy_family": "periodic",
            "mix_count": "2",
            "overhead": "0",
            "trials": "1",
        }
        for key, expected in checks.items():
            if row[key] != expected:
                fail(f"{path}: K={k} {key} mismatch")
        success = strict_uint(row["success"], f"{path}:K={k}:success")
        rank_fail = strict_uint(row["rank_fail"], f"{path}:K={k}:rank_fail")
        errors = strict_uint(row["error"], f"{path}:K={k}:error")
        if any(value not in (0, 1) for value in (success, rank_fail, errors)):
            fail(f"{path}: K={k} non-Boolean result count")
        if success + rank_fail + errors != 1:
            fail(f"{path}: K={k} trial ledger mismatch")
        failures = parse_failure_ids(row["failure_trials"], 1, f"{path}:K={k}")
        if len(failures) != rank_fail + errors:
            fail(f"{path}: K={k} failure list/count mismatch")
        expected_first_rank_fail = "0" if rank_fail else "-1"
        if row["first_rank_fail"] != expected_first_rank_fail:
            fail(f"{path}: K={k} first_rank_fail mismatch")
        for key in (
            "inact_mu", "binary_def_mu", "heavy_gain_mu", "solve_ms_mu",
            "build_ms_mu", "peel_ms_mu", "project_ms_mu", "residual_ms_mu",
            "backsub_ms_mu", "block_xors_mu", "block_muladds_mu",
        ):
            nonnegative_decimal(row[key], f"{path}:K={k}:{key}")
        if canonical_salt(row["active_packet_peel_seed_xor"], f"{path}:K={k}") != salt:
            fail(f"{path}: K={k} active salt mismatch")
    return rows


def validate_exact_job(args: argparse.Namespace) -> None:
    experiment = args.experiment.resolve()
    jobs = read_tsv_jobs(experiment / "meta/exact_jobs.tsv")
    if args.job_id < 0 or args.job_id >= len(jobs):
        fail(f"job ID {args.job_id} is outside [0,{len(jobs)})")
    job = jobs[args.job_id]
    if int(job["job_id"]) != args.job_id:
        fail(f"job row {args.job_id} has mismatched ID {job['job_id']}")
    stem = job["stem"]
    stderr = experiment / f"raw/exact/{stem}.stderr"
    if not stderr.is_file() or stderr.stat().st_size != 0:
        fail(f"nonempty or missing stderr for {stem}")
    validate_exact_output(experiment / f"raw/exact/{stem}.csv", job)


def expect_audit_failure(action: object, label: str) -> None:
    try:
        action()  # type: ignore[operator]
    except AuditError:
        return
    fail(f"selftest mutation unexpectedly passed: {label}")


def selftest(_args: argparse.Namespace) -> None:
    canonical_job = {
        "job_id": "0",
        "stem": "selftest",
        "schedule": "burst",
        "seed_index": "0",
        "seed": DISCOVERY_SEEDS[0],
        "loss": LOSS,
        "salt": "57",
        "k_count": "1",
        "k_csv": "7",
    }
    preamble = (
        "# precodefail: trials=1 threads=1 loss=0.5 "
        f"seed={int(DISCOVERY_SEEDS[0], 16)} "
        "completion=mixed mixed_period=32 mixed_gf256_rows=11 "
        "mixed_gf16_rows=4 mixed_geometry=shared-x mixed_residue_skew=0 "
        "mixed_residue_schedule=hashed mixed_residue_hash_seed=0x44 "
        "mixed_residue_hash_keyed=1 mixed_independent_extension_residues=1 "
        "mixed_extension_residue_seed_xor=0x4e source_hits_override=0 "
        "packet_peel_seed_xor=0x39 packet_peel_seed_table=none "
        "binary_dense_rows_override=0 gf256_heavy_rows_override=0 "
        "odd_packet_peel_seed_xor=0x0 packet_row_seed_multiplier=0x1 "
        "packet_row_seed_avalanche=0 seed_block_bytes_override=1280 "
        "overhead_stream=salted full_payload_solve=0 schedule=burst\n"
    )
    values = {
        "N": "7",
        "bb": "64",
        "heavy_family": "periodic",
        "mix_count": "2",
        "overhead": "0",
        "trials": "1",
        "success": "1",
        "rank_fail": "0",
        "error": "0",
        "fail_rate": "0.00000000",
        "inact_mu": "31.000",
        "inact_max": "31",
        "binary_def_mu": "15.000",
        "binary_def_max": "15",
        "heavy_gain_mu": "15.000",
        "heavy_gain_min": "15",
        "heavy_shortfall": "0",
        "solve_ms_mu": "0.010",
        "build_ms_mu": "0.001",
        "peel_ms_mu": "0.001",
        "project_ms_mu": "0.001",
        "residual_ms_mu": "0.001",
        "backsub_ms_mu": "0.001",
        "seed_attempt": "0",
        "block_xors_mu": "123.000",
        "block_muladds_mu": "456.000",
        "first_rank_fail": "-1",
        "binary_def_hist": "15:1",
        "heavy_gain_hist": "15:1",
        "failure_trials": "",
        "active_packet_peel_seed_xor": "57",
    }

    with tempfile.TemporaryDirectory(prefix="wh2-h15-v5-selftest-") as temp:
        valid = Path(temp) / "valid.csv"
        with valid.open("w", newline="") as stream:
            stream.write(preamble)
            writer = csv.DictWriter(stream, fieldnames=RAW_HEADER, lineterminator="\n")
            writer.writeheader()
            writer.writerow(values)
        validate_exact_output(valid, canonical_job)

        def mutated(name: str, old: str, new: str) -> Path:
            path = Path(temp) / name
            path.write_text(valid.read_text().replace(old, new, 1))
            return path

        bad_seed = mutated("bad-seed.csv", f"seed={int(DISCOVERY_SEEDS[0], 16)}", "seed=1")
        expect_audit_failure(lambda: validate_exact_output(bad_seed, canonical_job), "seed")
        bad_loss = mutated("bad-loss.csv", "loss=0.5", "loss=0.50")
        expect_audit_failure(lambda: validate_exact_output(bad_loss, canonical_job), "loss token")
        bad_salt = mutated("bad-salt.csv", "packet_peel_seed_xor=0x39", "packet_peel_seed_xor=0x38")
        expect_audit_failure(lambda: validate_exact_output(bad_salt, canonical_job), "salt")
        bad_k = mutated("bad-k.csv", "\n7,64,", "\n8,64,")
        expect_audit_failure(lambda: validate_exact_output(bad_k, canonical_job), "K coverage")
        bad_arm = mutated(
            "bad-arm.csv",
            "source_hits_override=0",
            "mixed_independent_gf256_breaker_residues=1 source_hits_override=0",
        )
        expect_audit_failure(lambda: validate_exact_output(bad_arm, canonical_job), "candidate arm")
        bad_count = mutated(
            "bad-count.csv",
            "\n7,64,periodic,2,0,1,1,0,0,",
            "\n7,64,periodic,2,0,1,1,-1,1,",
        )
        expect_audit_failure(lambda: validate_exact_output(bad_count, canonical_job), "negative result count")
        bad_first = mutated("bad-first.csv", ",-1,15:1,", ",0,15:1,")
        expect_audit_failure(lambda: validate_exact_output(bad_first, canonical_job), "first rank failure")

        valid_jobs = Path(temp) / "jobs.tsv"
        job_rows: list[dict[str, object]] = []
        for schedule in SCHEDULES:
            for seed_index, seed in enumerate(DISCOVERY_SEEDS):
                for salt in range(256):
                    job_id = len(job_rows)
                    job_rows.append({
                        "job_id": job_id,
                        "stem": f"job{job_id:04d}_{schedule}_seed{seed_index}_salt{salt:03d}",
                        "schedule": schedule,
                        "seed_index": seed_index,
                        "seed": seed,
                        "loss": LOSS,
                        "salt": salt,
                        "k_count": 1,
                        "k_csv": "7",
                    })
        write_csv(valid_jobs, EXACT_JOB_HEADER, job_rows, delimiter="\t")
        read_tsv_jobs(valid_jobs)
        bad_jobs = Path(temp) / "bad-jobs.tsv"
        bad_jobs.write_text(valid_jobs.read_text().replace("job0000_burst_seed0_salt000", "wrong", 1))
        expect_audit_failure(lambda: read_tsv_jobs(bad_jobs), "job stem")
        extra_jobs = Path(temp) / "extra-jobs.tsv"
        lines = valid_jobs.read_text().splitlines()
        lines[1] += "\textra"
        extra_jobs.write_text("\n".join(lines) + "\n")
        expect_audit_failure(lambda: read_tsv_jobs(extra_jobs), "extra TSV field")

    if parse_u64("57", "selftest") != parse_u64("0x39", "selftest"):
        fail("selftest numeric salt equivalence failed")
    expect_audit_failure(lambda: canonical_salt("256", "selftest"), "uint8 salt bound")
    expect_audit_failure(lambda: strict_uint("+1", "selftest"), "signed decimal")
    expect_audit_failure(lambda: strict_uint(" 1", "selftest"), "whitespace decimal")
    expect_audit_failure(lambda: parse_failure_ids("0|0", 1, "selftest"), "failure list duplicate")
    expect_audit_failure(lambda: parse_failure_ids("1", 1, "selftest"), "failure list range")
    expect_audit_failure(lambda: parse_failure_ids("+0", 1, "selftest"), "signed failure ID")
    print("selftest: PASS (valid fixtures + 15 rejected mutations)")


def analyze_exact(args: argparse.Namespace) -> None:
    experiment = args.experiment.resolve()
    jobs = read_tsv_jobs(experiment / "meta/exact_jobs.tsv")
    with (experiment / "meta/discovery_failures.csv").open(newline="") as stream:
        failures = list(csv.DictReader(stream))
    by_witness: dict[tuple[int, str, int], dict[str, str]] = {}
    by_k_witnesses: dict[int, list[dict[str, str]]] = defaultdict(list)
    for row in failures:
        key = (int(row["K"]), row["schedule"], int(row["seed_index"]))
        if key in by_witness:
            fail(f"duplicate discovery witness key {key}")
        by_witness[key] = row
        by_k_witnesses[int(row["K"])].append(row)

    exact_rows: dict[tuple[int, str, int, int], dict[str, str]] = {}
    for index, job in enumerate(jobs):
        if int(job["job_id"]) != index:
            fail("exact job IDs are not contiguous")
        stem = job["stem"]
        status = experiment / f"meta/exact_status/{stem}.ok"
        if not status.is_file() or status.read_text().strip() != "ok":
            fail(f"missing successful status for {stem}")
        stderr = experiment / f"raw/exact/{stem}.stderr"
        if not stderr.is_file() or stderr.stat().st_size != 0:
            fail(f"nonempty/missing stderr for {stem}")
        output = experiment / f"raw/exact/{stem}.csv"
        rows = validate_exact_output(output, job)
        salt = int(job["salt"])
        for k, row in rows.items():
            key = (k, job["schedule"], int(job["seed_index"]), salt)
            if key in exact_rows:
                fail(f"duplicate exact result {key}")
            exact_rows[key] = row

    for key, witness in by_witness.items():
        k, schedule, seed_index = key
        certified = int(witness["certified_v4_salt"], 16)
        row = exact_rows.get((k, schedule, seed_index, certified))
        if row is None:
            fail(f"missing certified-salt reproduction for witness {key}")
        if int(row["rank_fail"]) + int(row["error"]) <= 0:
            fail(f"certified-salt baseline failure did not reproduce for witness {key}")
        raw_fields = (
            ("rank_fail", "rank_fail"),
            ("error", "error"),
            ("failure_trials", "failure_trials"),
            ("heavy_shortfall", "heavy_shortfall"),
            ("inact_mu", "inact_mu"),
            ("inact_max", "inact_max"),
            ("binary_def_mu", "binary_def_mu"),
            ("binary_def_max", "binary_def_max"),
            ("seed_attempt", "seed_attempt"),
            ("block_xors_mu", "block_xors_mu"),
            ("block_muladds_mu", "block_muladds_mu"),
        )
        for result_key, witness_key in raw_fields:
            if row[result_key] != witness[witness_key]:
                fail(f"certified-salt reproduction field mismatch: {key} {result_key}")

    scores: list[dict[str, object]] = []
    shortlist: list[dict[str, object]] = []
    for k in sorted(by_k_witnesses):
        witnesses = by_k_witnesses[k]
        certified = int(witnesses[0]["certified_v4_salt"], 16)
        if any(int(row["certified_v4_salt"], 16) != certified for row in witnesses):
            fail(f"K={k}: inconsistent certified salt")
        eligible_results: list[dict[str, object]] = []
        for salt in range(256):
            rows = [
                exact_rows[(k, row["schedule"], int(row["seed_index"]), salt)]
                for row in witnesses
            ]
            rank_fail = sum(int(row["rank_fail"]) for row in rows)
            errors = sum(int(row["error"]) for row in rows)
            heavy_shortfall = sum(int(row["heavy_shortfall"]) for row in rows)
            block_xors = sum(Decimal(row["block_xors_mu"]) for row in rows)
            block_muladds = sum(Decimal(row["block_muladds_mu"]) for row in rows)
            inact_total = sum(Decimal(row["inact_mu"]) for row in rows)
            inact_max = max(int(row["inact_max"]) for row in rows)
            binary_def_total = sum(Decimal(row["binary_def_mu"]) for row in rows)
            binary_def_max = max(int(row["binary_def_max"]) for row in rows)
            eligible = rank_fail == 0 and errors == 0
            result: dict[str, object] = {
                "K": k,
                "certified_v4_salt": f"0x{certified:x}",
                "salt": f"0x{salt:x}",
                "witness_count": len(rows),
                "rank_fail": rank_fail,
                "error": errors,
                "heavy_shortfall": heavy_shortfall,
                "binary_def_mu_sum": str(binary_def_total),
                "binary_def_max": binary_def_max,
                "inact_mu_sum": str(inact_total),
                "inact_max": inact_max,
                "block_xors_mu_sum": str(block_xors),
                "block_muladds_mu_sum": str(block_muladds),
                "exact_eligible": int(eligible),
                "existing_v4_entry": int(certified != 0),
            }
            scores.append(result)
            if eligible:
                eligible_results.append(result)
        if len(eligible_results) < args.shortlist:
            fail(f"K={k}: only {len(eligible_results)} exact-eligible salts")

        def numeric(result: dict[str, object], key: str) -> Decimal:
            return Decimal(str(result[key]))

        def salt_number(result: dict[str, object]) -> int:
            return int(str(result["salt"]), 16)

        rankings = (
            ("structural", sorted(eligible_results, key=lambda result: (
                int(result["heavy_shortfall"]), int(result["binary_def_max"]),
                numeric(result, "binary_def_mu_sum"), int(result["inact_max"]),
                numeric(result, "inact_mu_sum"), numeric(result, "block_xors_mu_sum"),
                numeric(result, "block_muladds_mu_sum"), salt_number(result),
            ))),
            ("inactivation", sorted(eligible_results, key=lambda result: (
                int(result["heavy_shortfall"]), int(result["inact_max"]),
                numeric(result, "inact_mu_sum"), int(result["binary_def_max"]),
                numeric(result, "binary_def_mu_sum"), numeric(result, "block_xors_mu_sum"),
                numeric(result, "block_muladds_mu_sum"), salt_number(result),
            ))),
            ("work", sorted(eligible_results, key=lambda result: (
                int(result["heavy_shortfall"]), numeric(result, "block_xors_mu_sum"),
                numeric(result, "block_muladds_mu_sum"), int(result["binary_def_max"]),
                numeric(result, "binary_def_mu_sum"), int(result["inact_max"]),
                numeric(result, "inact_mu_sum"), salt_number(result),
            ))),
        )
        selected: list[tuple[str, dict[str, object]]] = []
        selected_salts: set[int] = set()
        depth = 0
        while len(selected) < args.shortlist:
            for basis, ranking in rankings:
                result = ranking[depth]
                salt = salt_number(result)
                if salt not in selected_salts:
                    selected_salts.add(salt)
                    selected.append((basis, result))
                    if len(selected) == args.shortlist:
                        break
            depth += 1
        for rank, (basis, result) in enumerate(selected, 1):
            shortlisted = dict(result)
            shortlisted["exact_rank"] = rank
            shortlisted["selection_basis"] = basis
            shortlist.append(shortlisted)

    score_header = (
        "K", "certified_v4_salt", "salt", "witness_count", "rank_fail",
        "error", "heavy_shortfall", "binary_def_mu_sum", "binary_def_max",
        "inact_mu_sum", "inact_max",
        "block_xors_mu_sum", "block_muladds_mu_sum", "exact_eligible",
        "existing_v4_entry",
    )
    shortlist_header = score_header + ("exact_rank", "selection_basis")
    write_csv(experiment / "analysis/exact_salt_scores.csv", score_header, scores)
    write_csv(experiment / "analysis/exact_shortlist.csv", shortlist_header, shortlist)

    with (experiment / "analysis/exact_summary.json").open("x") as stream:
        json.dump({
            "schema": "wirehair.wh2.normalized_h15_v5.exact_salt_search.v1",
            "jobs": len(jobs),
            "logical_cells": len(exact_rows),
            "discovery_witnesses": len(failures),
            "discovery_k": len(by_k_witnesses),
            "certified_v4_failures_reproduced": len(failures),
            "shortlist_per_k": args.shortlist,
            "shortlist_rows": len(shortlist),
            "ranking": [
                "zero exact-witness rank_fail/error required",
                "heavy_shortfall",
                "round-robin diversity across structural binary deficit, inactivation, and work rankings",
                "numeric salt is the final deterministic tie-break",
            ],
        }, stream, indent=2, sort_keys=True)
        stream.write("\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    extract_parser = subparsers.add_parser("extract")
    extract_parser.add_argument("--campaign", type=Path, required=True)
    extract_parser.add_argument("--source-seal", type=Path, required=True)
    extract_parser.add_argument("--analysis-complete", type=Path, required=True)
    extract_parser.add_argument("--output", type=Path, required=True)
    extract_parser.set_defaults(func=extract)

    exact_parser = subparsers.add_parser("analyze-exact")
    exact_parser.add_argument("--experiment", type=Path, required=True)
    exact_parser.add_argument("--shortlist", type=int, default=16)
    exact_parser.set_defaults(func=analyze_exact)

    validate_parser = subparsers.add_parser("validate-exact-job")
    validate_parser.add_argument("--experiment", type=Path, required=True)
    validate_parser.add_argument("--job-id", type=int, required=True)
    validate_parser.set_defaults(func=validate_exact_job)

    selftest_parser = subparsers.add_parser("selftest")
    selftest_parser.set_defaults(func=selftest)

    args = parser.parse_args()
    try:
        args.func(args)
    except (AuditError, OSError, ValueError, csv.Error) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
