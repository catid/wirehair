#!/usr/bin/env python3
"""Validate, monitor, seal, and reduce normalized-H15-v5 holdout results."""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import io
import json
import math
import os
import re
import shlex
import stat
import sys
import tempfile
import time
from collections import defaultdict
from collections.abc import Callable
from decimal import Decimal, InvalidOperation
from pathlib import Path

sys.dont_write_bytecode = True
import wh2_h15_v5_holdout_prepare as contract


RAW_HEADER = (
    "N", "bb", "heavy_family", "mix_count", "overhead", "trials", "success",
    "rank_fail", "error", "fail_rate", "inact_mu", "inact_max", "binary_def_mu",
    "binary_def_max", "heavy_gain_mu", "heavy_gain_min", "heavy_shortfall",
    "solve_ms_mu", "build_ms_mu", "peel_ms_mu", "project_ms_mu", "residual_ms_mu",
    "backsub_ms_mu", "seed_attempt", "block_xors_mu", "block_muladds_mu",
    "first_rank_fail", "binary_def_hist", "heavy_gain_hist", "failure_trials",
    "active_packet_peel_seed_xor",
    "mixed_joint_source_xors_mu", "mixed_joint_marginal_xors_mu",
    "mixed_joint_marginal_copies_mu", "mixed_joint_active_deltas_mu",
    "mixed_joint_scratch_bytes_mu", "mixed_dual_source_columns_mu",
)
JOINT_METRIC_FIELDS = (
    "mixed_joint_source_xors_mu", "mixed_joint_marginal_xors_mu",
    "mixed_joint_marginal_copies_mu", "mixed_joint_active_deltas_mu",
    "mixed_joint_scratch_bytes_mu", "mixed_dual_source_columns_mu",
)
THERMAL_HEADER = (
    "utc", "monotonic_s", "cpu_busy_pct", "cpu_avg_mhz", "cpu_tctl_c",
    "dimm_i2c1_50_c", "dimm_i2c1_51_c", "dimm_i2c1_52_c", "dimm_i2c1_53_c",
    "dimm_i2c2_50_c", "dimm_i2c2_51_c", "dimm_i2c2_52_c", "dimm_i2c2_53_c",
    "dimm_read_errors", "load1", "load5", "load15", "edac_ce", "edac_ue",
)
TIMING_RE = re.compile(
    r"elapsed_s=([0-9]+(?:\.[0-9]+)?) user_s=([0-9]+(?:\.[0-9]+)?) "
    r"sys_s=([0-9]+(?:\.[0-9]+)?) cpu_pct=([0-9]+)% "
    r"max_rss_kb=(0|[1-9][0-9]*) exit=0\n"
)
SEMANTIC_FIELDS = tuple(field for field in RAW_HEADER if field not in (
    "solve_ms_mu", "build_ms_mu", "peel_ms_mu", "project_ms_mu",
    "residual_ms_mu", "backsub_ms_mu",
))


AuditError = contract.AuditError
fail = contract.fail


def decimal(text: str, label: str) -> Decimal:
    try:
        value = Decimal(text)
    except InvalidOperation:
        fail(f"{label}: invalid decimal {text!r}")
    if not value.is_finite() or value < 0:
        fail(f"{label}: decimal must be finite and nonnegative")
    return value


def parse_failures(text: str, trials: int, label: str) -> tuple[int, ...]:
    if not text:
        return ()
    result = tuple(contract.uint(value, label, trials - 1) for value in text.split("|"))
    if any(left >= right for left, right in zip(result, result[1:])):
        fail(f"{label}: failure trial IDs are duplicate or unsorted")
    return result


def parse_histogram(text: str, trials: int, label: str) -> dict[int, int]:
    result: dict[int, int] = {}
    for item in text.split("|"):
        fields = item.split(":")
        if len(fields) != 2:
            fail(f"{label}: malformed histogram")
        key = contract.uint(fields[0], label)
        count = contract.uint(fields[1], label)
        if key in result or count == 0 or (result and key <= max(result)):
            fail(f"{label}: histogram is not canonical")
        result[key] = count
    if sum(result.values()) != trials:
        fail(f"{label}: histogram count does not equal trials")
    return result


def recover_integer_total(text: str, trials: int, label: str) -> int:
    if not re.fullmatch(r"[0-9]+\.[0-9]{3}", text):
        fail(f"{label}: work mean must use the CLI's exact three-decimal format")
    mean = decimal(text, label)
    scaled = mean * trials
    nearest = scaled.to_integral_value()
    tolerance = Decimal("0.000500001") * trials
    if abs(scaled - nearest) > tolerance:
        fail(f"{label}: rounded mean cannot represent a unique integer trial sum")
    total = int(nearest)
    if f"{Decimal(total) / trials:.3f}" != text:
        fail(f"{label}: mean does not round back from the recovered integer sum")
    return total


def parse_preamble(line: str, label: str) -> dict[str, str]:
    prefix = "# precodefail: "
    if not line.startswith(prefix):
        fail(f"{label}: missing precodefail preamble")
    result: dict[str, str] = {}
    for item in line[len(prefix):].split():
        if "=" not in item:
            fail(f"{label}: malformed preamble token")
        key, value = item.split("=", 1)
        if not key or key in result:
            fail(f"{label}: duplicate/empty preamble key")
        result[key] = value
    return result


def load_jobs(root: Path) -> tuple[list[dict[str, object]], dict[int, dict[str, str]]]:
    contract.ensure_private_root(root, fresh=False)
    contract.verify_seal(root, "meta/pre_reveal_seal.sha256", contract.PRE_REVEAL_ROLES)
    contract.verify_seal(root, "recovery/input_seal.sha256", contract.RECOVERY_INPUT_ROLES)
    manifest_data = contract.read_relative(root, "recovery/manifest.json")
    manifest = contract.parse_json(manifest_data, "recovery manifest")
    if not isinstance(manifest, dict):
        fail("recovery manifest root must be an object")
    reveal_data = contract.read_relative(root, "recovery/holdout_reveal.tsv")
    reveals = contract.parse_reveal(reveal_data)
    groups = contract.read_groups(contract.read_relative(root, "freeze/groups.tsv"))
    table = contract.parse_table(contract.read_relative(root, "freeze/selection_table.tsv"))
    changed = tuple(k for k, row in table.items() if row["changed"] == "1")
    expected = contract.make_recovery_jobs(groups, reveals, changed)
    jobs_data = contract.read_relative(root, "recovery/jobs.tsv")
    if jobs_data != expected or contract.sha256_bytes(jobs_data) != manifest.get("jobs_sha256"):
        fail("recovery jobs differ from exact seed/table-derived reconstruction")
    if manifest != {
        "artifact_records": contract.EXPECTED_RECOVERY_ARTIFACTS,
        "binary_sha256": contract.sha256_bytes(contract.read_relative(root, "freeze/wirehair_v2_bench")),
        "holdout_reveal_sha256": contract.sha256_bytes(reveal_data),
        "jobs": contract.EXPECTED_RECOVERY_JOBS,
        "jobs_sha256": contract.sha256_bytes(jobs_data),
        "paired_rows": contract.EXPECTED_PAIRED_ROWS,
        "pre_reveal_seal_sha256": contract.sha256_bytes(contract.read_relative(root, "meta/pre_reveal_seal.sha256")),
        "schema": "wirehair.wh2.normalized_h15_v5.holdout.recovery.v1",
        "seed_roles": {"breadth": "raw holdout derived_u64", "depth": "domain-separated SHA-256 derivative"},
    }:
        fail("recovery manifest contract mismatch")
    raw_rows = contract.parse_tsv(jobs_data, contract.RECOVERY_JOB_HEADER, "recovery jobs")
    jobs: list[dict[str, object]] = []
    for line, row in enumerate(raw_rows, 2):
        job: dict[str, object] = dict(row)
        job["_id"] = contract.uint(row["job_id"], f"jobs:{line}:id", 17423)
        job["_trials"] = contract.uint(row["trials"], f"jobs:{line}:trials", 16)
        job["_ks"] = tuple(contract.uint(k, f"jobs:{line}:K", 64000) for k in row["k_csv"].split(","))
        if job["_id"] != len(jobs) or len(job["_ks"]) != int(row["k_count"]):
            fail(f"jobs:{line}: ID/K count mismatch")
        jobs.append(job)
    return jobs, table


def expected_preamble(job: dict[str, object]) -> dict[str, str]:
    return {
        "trials": str(job["trials"]),
        "threads": "1",
        "loss": format(float(str(job["loss"])), ".17g"),
        "seed": f"0x{int(str(job['seed']), 16):x}",
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
        "mixed_residue_buckets_requested": "auto",
        "mixed_extension_residue_seed_xor": "0x4e",
        "source_hits_override": "0",
        "packet_peel_seed_table": str(job["profile"]),
        "binary_dense_rows_override": "0",
        "gf256_heavy_rows_override": "0",
        "odd_packet_peel_seed_xor": "0x0",
        "packet_row_seed_multiplier": "0x1",
        "packet_row_seed_avalanche": "0",
        "seed_block_bytes_override": "1280",
        "overhead_stream": "salted",
        "full_payload_solve": "0",
        "schedule": str(job["schedule"]),
    }


def validate_output_data(data: bytes, stderr: bytes, job: dict[str, object],
                         table: dict[int, dict[str, str]] | None,
                         label: str) -> dict[int, dict[str, object]]:
    if stderr:
        fail(f"{label}: stderr must be empty")
    if b"\r" in data or not data.endswith(b"\n") or b"\n\n" in data:
        fail(f"{label}: output framing is noncanonical")
    try:
        text = data.decode("ascii")
    except UnicodeError:
        fail(f"{label}: output must be ASCII")
    lines = text.splitlines()
    if len(lines) != len(job["_ks"]) + 2:
        fail(f"{label}: expected one row per K")
    preamble = parse_preamble(lines[0], label)
    if set(preamble) != set(expected_preamble(job)) | {"packet_peel_seed_xor"}:
        fail(f"{label}: unexpected or missing preamble keys")
    for key, expected in expected_preamble(job).items():
        if preamble.get(key) != expected:
            fail(f"{label}: preamble {key} expected {expected!r}, found {preamble.get(key)!r}")
    if preamble.get("packet_peel_seed_xor") not in (None, "0x0"):
        fail(f"{label}: explicit packet peel XOR must remain zero")
    reader = csv.DictReader(io.StringIO("\n".join(lines[1:]) + "\n", newline=""))
    if tuple(reader.fieldnames or ()) != RAW_HEADER:
        fail(f"{label}: raw header mismatch")
    parsed: dict[int, dict[str, object]] = {}
    for offset, (row, expected_k) in enumerate(zip(reader, job["_ks"]), 3):
        if None in row or any(value is None for value in row.values()):
            fail(f"{label}:{offset}: malformed CSV")
        k = contract.uint(row["N"], f"{label}:{offset}:N", 64000)
        trials = contract.uint(row["trials"], f"{label}:{offset}:trials", 16)
        success = contract.uint(row["success"], f"{label}:{offset}:success", trials)
        rank_fail = contract.uint(row["rank_fail"], f"{label}:{offset}:rank_fail", trials)
        errors = contract.uint(row["error"], f"{label}:{offset}:error", trials)
        failures = parse_failures(row["failure_trials"], trials, f"{label}:{offset}:failures")
        if (k != expected_k or trials != job["_trials"] or
                row["bb"] != "64" or row["heavy_family"] != "periodic" or
                row["mix_count"] != "2" or row["overhead"] != "0" or
                success + rank_fail + errors != trials or len(failures) != rank_fail + errors):
            fail(f"{label}:{offset}: geometry/outcome count mismatch")
        if table is not None:
            expected_salt = table[k]["v4_salt" if job["arm"] == "v4" else "v5_salt"]
            if row["active_packet_peel_seed_xor"] != expected_salt:
                fail(f"{label}:{offset}: active peel salt differs from frozen table")
        if decimal(row["fail_rate"], f"{label}:{offset}:fail_rate") != \
                Decimal(rank_fail + errors) / Decimal(trials):
            # The CLI prints eight decimals, so compare its exact rounding.
            expected_rate = Decimal(f"{(rank_fail + errors) / trials:.8f}")
            if decimal(row["fail_rate"], f"{label}:{offset}:fail_rate") != expected_rate:
                fail(f"{label}:{offset}: fail rate mismatch")
        for field in ("inact_mu", "binary_def_mu", "heavy_gain_mu", "solve_ms_mu",
                      "build_ms_mu", "peel_ms_mu", "project_ms_mu", "residual_ms_mu",
                      "backsub_ms_mu", "block_xors_mu", "block_muladds_mu",
                      *JOINT_METRIC_FIELDS):
            decimal(row[field], f"{label}:{offset}:{field}")
        for field in ("inact_max", "binary_def_max", "heavy_gain_min", "seed_attempt"):
            contract.uint(row[field], f"{label}:{offset}:{field}")
        contract.uint(row["heavy_shortfall"], f"{label}:{offset}:heavy_shortfall", trials)
        first = row["first_rank_fail"]
        if rank_fail == 0:
            if first != "-1":
                fail(f"{label}:{offset}: first rank failure must be -1")
        elif contract.uint(first, f"{label}:{offset}:first failure", trials - 1) not in failures:
            fail(f"{label}:{offset}: first rank failure is absent from failure IDs")
        binary_hist = parse_histogram(
            row["binary_def_hist"], trials, f"{label}:{offset}:binary histogram")
        heavy_hist = parse_histogram(
            row["heavy_gain_hist"], trials, f"{label}:{offset}:gain histogram")
        binary_mean = sum(Decimal(key * count) for key, count in binary_hist.items()) / trials
        heavy_mean = sum(Decimal(key * count) for key, count in heavy_hist.items()) / trials
        if (max(binary_hist) != int(row["binary_def_max"]) or
                min(heavy_hist) != int(row["heavy_gain_min"]) or
                f"{binary_mean:.3f}" != row["binary_def_mu"] or
                f"{heavy_mean:.3f}" != row["heavy_gain_mu"]):
            fail(f"{label}:{offset}: histogram summary mismatch")
        xor_total = recover_integer_total(
            row["block_xors_mu"], trials, f"{label}:{offset}:block_xors_mu")
        mul_total = recover_integer_total(
            row["block_muladds_mu"], trials, f"{label}:{offset}:block_muladds_mu")
        for field in JOINT_METRIC_FIELDS:
            recover_integer_total(row[field], trials, f"{label}:{offset}:{field}")
        cooked: dict[str, object] = dict(row)
        cooked.update({"_trials": trials, "_rank_fail": rank_fail, "_errors": errors,
                       "_failures": failures,
                       "_xor": decimal(row["block_xors_mu"], f"{label}:{offset}:xor"),
                       "_mul": decimal(row["block_muladds_mu"], f"{label}:{offset}:mul"),
                       "_xor_total": xor_total, "_mul_total": mul_total})
        parsed[k] = cooked
    if tuple(parsed) != job["_ks"]:
        fail(f"{label}: K output order differs from job ledger")
    return parsed


def output_paths(root: Path, job: dict[str, object]) -> tuple[str, str, str, str, str]:
    stem = str(job["stem"])
    return (f"recovery/raw/{stem}.csv", f"recovery/raw/{stem}.stderr",
            f"recovery/raw/{stem}.time", f"recovery/status/{stem}.ok",
            f"recovery/commands/{stem}.txt")


def validate_output(args: argparse.Namespace) -> None:
    jobs, table = load_jobs(args.root)
    job_id = contract.uint(args.job_id, "job ID", 17423)
    data = contract.read_once(args.stdout)
    stderr = contract.read_once(args.stderr)
    validate_output_data(data, stderr, jobs[job_id], table, str(args.stdout))


def parse_job_argument(text: str) -> dict[str, object]:
    if "\n" in text or "\r" in text:
        fail("temporary job row contains a line break")
    data = contract.canonical_tsv(contract.RECOVERY_JOB_HEADER, (text.split("\t"),))
    rows = contract.parse_tsv(data, contract.RECOVERY_JOB_HEADER, "temporary job row")
    if len(rows) != 1:
        fail("temporary job row count mismatch")
    row = rows[0]
    job: dict[str, object] = dict(row)
    job["_id"] = contract.uint(row["job_id"], "temporary job ID", 17423)
    job["_trials"] = contract.uint(row["trials"], "temporary trials", 16)
    job["_ks"] = tuple(contract.uint(k, "temporary K", 64000) for k in row["k_csv"].split(","))
    if len(job["_ks"]) != contract.uint(row["k_count"], "temporary K count"):
        fail("temporary job K count mismatch")
    if (row["arm"], row["profile"]) not in contract.ARMS:
        fail("temporary arm/profile mismatch")
    return job


def validate_temp(args: argparse.Namespace) -> None:
    job = parse_job_argument(args.job_row)
    data = contract.read_once(args.stdout)
    stderr = contract.read_once(args.stderr)
    validate_output_data(data, stderr, job, None, str(args.stdout))


def expected_command(job: dict[str, object], binary: str) -> list[str]:
    return [binary, "precodefail", "--N", str(job["k_csv"]), "--bb-list", "64",
            "--seed-block-bytes", "1280", "--overhead", "0", "--trials", str(job["trials"]),
            "--threads", "1", "--loss", str(job["loss"]), "--seed", str(job["seed"]),
            "--schedule", str(job["schedule"]), "--completion", "mixed", "--mix-count", "2",
            "--packet-peel-seed-table", str(job["profile"]), "--mixed-gf256-rows", "11",
            "--mixed-gf16-rows", "4", "--mixed-period", "32", "--mixed-geometry", "shared-x",
            "--mixed-residue-schedule", "hashed", "--mixed-residue-hash-seed", "68",
            "--mixed-residue-hash-keyed", "--mixed-independent-extension-residues",
            "--mixed-extension-residue-seed-xor", "78"]


def validate_job_files(root: Path, job: dict[str, object], table: dict[int, dict[str, str]]) -> None:
    stdout_path, stderr_path, timing_path, status_path, command_path = output_paths(root, job)
    stdout = contract.read_relative(root, stdout_path)
    stderr = contract.read_relative(root, stderr_path)
    validate_output_data(stdout, stderr, job, table, stdout_path)
    timing = contract.read_relative(root, timing_path).decode("ascii")
    if not TIMING_RE.fullmatch(timing):
        fail(f"{timing_path}: malformed /usr/bin/time record")
    if contract.read_relative(root, status_path) != b"ok\n":
        fail(f"{status_path}: status is not exact ok")
    command_data = contract.read_relative(root, command_path)
    if b"\r" in command_data or not command_data.endswith(b"\n") or command_data.count(b"\n") != 1:
        fail(f"{command_path}: noncanonical command record")
    text = command_data.decode("ascii").rstrip("\n")
    prefix = (f"job_id={job['job_id']}\tphase={job['phase']}\tarm={job['arm']}\t"
              f"seed_index={job['seed_index']}\tk_count={job['k_count']}\tcommand=")
    if not text.startswith(prefix):
        fail(f"{command_path}: command provenance prefix mismatch")
    actual = shlex.split(text[len(prefix):])
    binary = str((root / "freeze/wirehair_v2_bench").resolve())
    if actual != expected_command(job, binary):
        fail(f"{command_path}: command differs from exact reconstruction")


def validate_job(args: argparse.Namespace) -> None:
    jobs, table = load_jobs(args.root)
    job_id = contract.uint(args.job_id, "job ID", 17423)
    validate_job_files(args.root, jobs[job_id], table)


def audit_inputs(args: argparse.Namespace) -> None:
    jobs, _ = load_jobs(args.root)
    if len(jobs) != contract.EXPECTED_RECOVERY_JOBS:
        fail("input audit job count mismatch")
    print("input audit: PASS (19-role pre-reveal + 6-role recovery-input seals)")


def expected_artifact_paths(jobs: list[dict[str, object]]) -> list[str]:
    paths: list[str] = []
    for job in jobs:
        paths.extend(output_paths(Path("."), job))
    paths.extend(("recovery/run_start.json", "recovery/run_finish.json",
                  "recovery/thermal_interval.csv"))
    if len(paths) != contract.EXPECTED_RECOVERY_ARTIFACTS or len(set(paths)) != len(paths):
        fail("internal recovery artifact contract mismatch")
    return paths


def verify_recovery_filename_paths(root: Path, paths: list[str], sealed: bool) -> None:
    for relative in ("recovery/raw", "recovery/status", "recovery/commands"):
        expected = [Path(path).name for path in paths
                    if Path(path).parent.as_posix() == relative]
        contract.verify_exact_directory_entries(root, relative, expected)
    root_files = [Path(path).name for path in paths
                  if Path(path).parent.as_posix() == "recovery"]
    root_files.extend(("holdout_reveal.tsv", "reveal_provenance.json", "manifest.json",
                       "jobs.tsv", "input_seal.sha256"))
    if sealed:
        root_files.append("artifact_manifest.sha256")
    contract.verify_exact_directory_entries(
        root, "recovery", root_files, ("raw", "status", "commands"))


def verify_artifact_filename_sets(
        root: Path, jobs: list[dict[str, object]], sealed: bool) -> None:
    verify_recovery_filename_paths(root, expected_artifact_paths(jobs), sealed)


def audit_artifacts(args: argparse.Namespace) -> None:
    jobs, table = load_jobs(args.root)
    for job in jobs:
        validate_job_files(args.root, job, table)
    verify_artifact_filename_sets(args.root, jobs, sealed=False)
    print(f"artifact audit: PASS ({len(expected_artifact_paths(jobs))} exact records)")


def make_artifact_manifest(root: Path, jobs: list[dict[str, object]]) -> bytes:
    return contract.canonical_tsv(("sha256", "path"), (
        (contract.sha256_bytes(contract.read_relative(root, path)), path)
        for path in expected_artifact_paths(jobs)
    ))


def seal_artifacts(args: argparse.Namespace) -> None:
    jobs, _ = load_jobs(args.root)
    verify_artifact_filename_sets(args.root, jobs, sealed=False)
    data = make_artifact_manifest(args.root, jobs)
    contract.write_exclusive(args.root, "recovery/artifact_manifest.sha256", data)
    contract.verify_hash_manifest(args.root, "recovery/artifact_manifest.sha256",
                                  contract.EXPECTED_RECOVERY_ARTIFACTS)
    verify_artifact_filename_sets(args.root, jobs, sealed=True)
    print(f"sealed {contract.EXPECTED_RECOVERY_ARTIFACTS} recovery artifacts")


def verify_artifacts(args: argparse.Namespace) -> None:
    jobs, _ = load_jobs(args.root)
    verify_artifact_filename_sets(args.root, jobs, sealed=True)
    expected = make_artifact_manifest(args.root, jobs)
    actual = contract.read_relative(args.root, "recovery/artifact_manifest.sha256")
    if actual != expected:
        fail("recovery artifacts changed after sealing")
    contract.verify_hash_manifest(args.root, "recovery/artifact_manifest.sha256",
                                  contract.EXPECTED_RECOVERY_ARTIFACTS)
    verify_artifact_filename_sets(args.root, jobs, sealed=True)
    print("recovery artifact seal: PASS")


def read_edac() -> tuple[int, int]:
    def total(kind: str) -> int:
        base = Path("/sys/devices/system/edac/mc")
        paths = sorted(base.glob(f"mc*/{kind}_count"))
        if not paths:
            fail(f"EDAC {kind.upper()} counters are unavailable")
        result = 0
        for path in paths:
            flags = os.O_RDONLY | os.O_CLOEXEC
            if hasattr(os, "O_NOFOLLOW"):
                flags |= os.O_NOFOLLOW
            fd = os.open(path, flags)
            try:
                data = os.read(fd, 128).decode("ascii").strip()
            finally:
                os.close(fd)
            result += contract.uint(data, str(path))
        return result
    return total("ce"), total("ue")


def read_uptime() -> float:
    flags = os.O_RDONLY | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    fd = os.open("/proc/uptime", flags)
    try:
        data = os.read(fd, 256)
    finally:
        os.close(fd)
    fields = data.split()
    if len(fields) != 2 or not re.fullmatch(rb"[0-9]+(?:\.[0-9]+)?", fields[0]):
        fail("/proc/uptime returned a malformed record")
    return float(fields[0])


def parse_thermal_line(data: bytes, line_number: int) -> tuple[list[str], list[float | None]]:
    if data.endswith(b"\r\n"):
        data = data[:-2] + b"\n"
    if not data.endswith(b"\n") or b"\r" in data:
        fail(f"thermal:{line_number}: incomplete/noncanonical line")
    try:
        fields = next(csv.reader(io.StringIO(data.decode("ascii"), newline="")))
    except (UnicodeError, csv.Error, StopIteration):
        fail(f"thermal:{line_number}: malformed CSV")
    if len(fields) != len(THERMAL_HEADER):
        fail(f"thermal:{line_number}: expected {len(THERMAL_HEADER)} fields")
    if not re.fullmatch(r"[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}(?:\.[0-9]+)?Z", fields[0]):
        fail(f"thermal:{line_number}: UTC timestamp is noncanonical")
    numeric: list[float | None] = []
    for index, value in enumerate(fields):
        if index == 0:
            numeric.append(None)
            continue
        if 5 <= index <= 12 and value == "":
            numeric.append(None)
            continue
        if not re.fullmatch(r"-?[0-9]+(?:\.[0-9]+)?", value):
            fail(f"thermal:{line_number}: nonnumeric field {THERMAL_HEADER[index]}")
        number = float(value)
        if not math.isfinite(number):
            fail(f"thermal:{line_number}: nonfinite field")
        numeric.append(number)
    if not (numbers := numeric):
        fail(f"thermal:{line_number}: internal numeric parse failure")
    if (numbers[1] is None or numbers[1] <= 0 or numbers[2] is None or
            not 0 <= numbers[2] <= 100 or numbers[3] is None or numbers[3] < 0 or
            numbers[4] is None or not 0 <= numbers[4] < 150):
        fail(f"thermal:{line_number}: CPU monotonic/load/temperature range violation")
    for index in range(5, 13):
        if numbers[index] is not None and not 0 <= numbers[index] < 150:
            fail(f"thermal:{line_number}: DIMM temperature range violation")
    if any(numbers[index] is None or numbers[index] < 0 for index in range(14, 17)):
        fail(f"thermal:{line_number}: load averages must be nonnegative")
    contract.uint(fields[13], f"thermal:{line_number}:dimm_read_errors", 8)
    contract.uint(fields[17], f"thermal:{line_number}:edac_ce")
    contract.uint(fields[18], f"thermal:{line_number}:edac_ue")
    return fields, numeric


def watchdog(args: argparse.Namespace) -> None:
    source = args.thermal_source
    flags = os.O_RDONLY | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    fd = os.open(source, flags)
    try:
        st = os.fstat(fd)
        if not stat.S_ISREG(st.st_mode):
            fail("thermal source must be a regular file")
        stream = os.fdopen(fd, "rb", buffering=1024 * 1024, closefd=False)
        header = stream.readline()
        normalized_header = header.replace(b"\r\n", b"\n")
        if normalized_header != (",".join(THERMAL_HEADER) + "\n").encode("ascii"):
            fail("thermal source header mismatch")
        for line_number in range(2, args.start_line + 1):
            line = stream.readline()
            if not line or not line.endswith(b"\n"):
                fail("thermal source ended before the frozen start line")
        ce_start, ue_start = read_edac()
        wall_start = time.monotonic()
        uptime_start = read_uptime()
        samples: list[tuple[bytes, list[str], list[float | None]]] = []
        last_sample_wall = wall_start
        hot_streak = 0
        missing_streak = [0] * 8
        max_missing_streak = [0] * 8
        line_number = args.start_line
        while True:
            line = stream.readline()
            if line:
                line_number += 1
                fields, numbers = parse_thermal_line(line, line_number)
                normalized = line.replace(b"\r\n", b"\n")
                last_sample_wall = time.monotonic()
                dims = numbers[5:13]
                missing = sum(value is None for value in dims)
                if missing == 8:
                    fail("live watchdog observed an all-DIMM blackout")
                if int(numbers[13] or 0) != missing:
                    fail("live watchdog DIMM read-error count differs from blank sensors")
                for index, value in enumerate(dims):
                    missing_streak[index] = missing_streak[index] + 1 if value is None else 0
                    max_missing_streak[index] = max(max_missing_streak[index], missing_streak[index])
                    if missing_streak[index] > 5:
                        fail(f"live watchdog sensor {index} exceeded five consecutive misses")
                temperatures = [numbers[4]] + [value for value in dims if value is not None]
                hot_streak = hot_streak + 1 if any(value is not None and value >= 90.0 for value in temperatures) else 0
                if hot_streak >= 3:
                    fail("live watchdog observed three consecutive >=90C samples")
                if int(numbers[17] or 0) != ce_start or int(numbers[18] or 0) != ue_start:
                    fail("live watchdog observed an EDAC counter change")
                samples.append((normalized, fields, numbers))
                if args.done_file.exists():
                    break
            else:
                if time.monotonic() - last_sample_wall > 5.0:
                    fail("live watchdog received no complete thermal sample for five seconds")
                if args.done_file.exists():
                    break
                time.sleep(0.1)
        wall_end = time.monotonic()
        uptime_end = read_uptime()
        ce_end, ue_end = read_edac()
        if (ce_end, ue_end) != (ce_start, ue_start):
            fail("EDAC counters changed across watchdog interval")
        if len(samples) < 2:
            fail("watchdog captured fewer than two campaign samples")
        monotonic = [float(sample[1][1]) for sample in samples]
        busy = [float(sample[1][2]) for sample in samples]
        gaps = [right - left for left, right in zip(monotonic, monotonic[1:])]
        if (any(gap <= 0 or gap > 2.5 for gap in gaps) or
                abs(monotonic[0] - uptime_start) > 2.5 or abs(uptime_end - monotonic[-1]) > 2.5):
            fail("thermal samples violate monotonic gap/endpoint gates")
        if min(busy) < 90.0 or sum(busy) / len(busy) < 98.0:
            fail("CPU utilization violates the 90% sample / 98% mean gates")
        coverage = []
        dimm_max = []
        for index in range(8):
            valid = [sample[2][5 + index] for sample in samples if sample[2][5 + index] is not None]
            coverage.append(len(valid) / len(samples))
            dimm_max.append(max(valid) if valid else None)
        if any(value < 0.99 for value in coverage):
            fail("one or more DIMM sensors fell below 99% coverage")
        output = (",".join(THERMAL_HEADER) + "\n").encode("ascii") + b"".join(sample[0] for sample in samples)
        report = {
            "cpu_busy_mean_pct": sum(busy) / len(busy),
            "cpu_busy_min_pct": min(busy),
            "cpu_tctl_max_c": max(float(sample[1][4]) for sample in samples),
            "dimm_coverage": coverage,
            "dimm_max_c": dimm_max,
            "edac_ce_end": ce_end,
            "edac_ce_start": ce_start,
            "edac_ue_end": ue_end,
            "edac_ue_start": ue_start,
            "gate": "PASS",
            "max_consecutive_dimm_misses": max_missing_streak,
            "max_gap_s": max(gaps),
            "samples": len(samples),
            "schema": "wirehair.wh2.holdout.environment.v1",
            "uptime_end": uptime_end,
            "uptime_start": uptime_start,
            "wall_elapsed_s": wall_end - wall_start,
        }
        args.output.write_bytes(output)
        args.report.write_bytes(contract.canonical_json(report))
    finally:
        os.close(fd)


def run_record(args: argparse.Namespace) -> None:
    now = dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    ce, ue = read_edac()
    value: dict[str, object] = {
        "edac_ce": ce, "edac_ue": ue, "kind": args.kind,
        "schema": "wirehair.wh2.holdout.run_record.v1", "utc": now,
        "uptime_s": read_uptime(),
    }
    if args.kind == "start":
        value.update({"job_count": contract.EXPECTED_RECOVERY_JOBS,
                      "launcher_nice": args.launcher_nice, "thermal_line": args.thermal_line,
                      "thermal_source": str(args.thermal_source), "workers": args.workers})
    else:
        report = contract.parse_json(contract.read_once(args.watchdog_report), "watchdog report")
        value.update({"watchdog": report, "xargs_rc": args.xargs_rc})
    args.output.write_bytes(contract.canonical_json(value))


def audit_environment_data(root: Path, panel: str) -> dict[str, object]:
    if panel == "recovery":
        prefix = "recovery"
    elif panel in ("payload", "lookup"):
        prefix = f"performance/{panel}"
    else:
        fail(f"unknown environment panel {panel!r}")
    interval_data = contract.read_relative(root, f"{prefix}/thermal_interval.csv")
    expected_header = (",".join(THERMAL_HEADER) + "\n").encode("ascii")
    if not interval_data.startswith(expected_header) or b"\r" in interval_data:
        fail(f"{panel}: thermal interval header/framing mismatch")
    lines = interval_data[len(expected_header):].splitlines(keepends=True)
    if len(lines) < 2 or b"".join(lines) != interval_data[len(expected_header):]:
        fail(f"{panel}: thermal interval has incomplete or insufficient samples")
    samples = [parse_thermal_line(line, index + 2) for index, line in enumerate(lines)]
    finish = contract.parse_json(
        contract.read_relative(root, f"{prefix}/run_finish.json"), f"{panel} run finish")
    start = contract.parse_json(
        contract.read_relative(root, f"{prefix}/run_start.json"), f"{panel} run start")
    if not isinstance(finish, dict) or not isinstance(start, dict):
        fail(f"{panel}: run records must be JSON objects")
    if panel == "recovery":
        start_keys = {"edac_ce", "edac_ue", "job_count", "kind", "launcher_nice", "schema",
                      "thermal_line", "thermal_source", "uptime_s", "utc", "workers"}
        finish_keys = {"edac_ce", "edac_ue", "kind", "schema", "uptime_s", "utc",
                       "watchdog", "xargs_rc"}
        report = finish.get("watchdog")
        if (set(start) != start_keys or set(finish) != finish_keys or
                start.get("schema") != "wirehair.wh2.holdout.run_record.v1" or
                finish.get("schema") != "wirehair.wh2.holdout.run_record.v1" or
                start.get("kind") != "start" or finish.get("kind") != "finish" or
                start.get("job_count") != contract.EXPECTED_RECOVERY_JOBS or
                start.get("workers") != 120 or start.get("launcher_nice") != 0 or
                finish.get("xargs_rc") != 0):
            fail("recovery run-record contract mismatch")
    else:
        start_keys = {"jobs", "panel", "schema", "thermal_line_start", "uptime_s", "utc"}
        finish_keys = {"environment", "jobs", "panel", "schema", "uptime_s", "utc"}
        report = finish.get("environment")
        expected_jobs = 32 if panel == "payload" else 24
        if (set(start) != start_keys or set(finish) != finish_keys or
                start.get("schema") != "wirehair.wh2.holdout.performance_run.v1" or
                finish.get("schema") != "wirehair.wh2.holdout.performance_run.v1" or
                start.get("panel") != panel or finish.get("panel") != panel or
                start.get("jobs") != expected_jobs or finish.get("jobs") != expected_jobs):
            fail(f"{panel}: performance run-record contract mismatch")
    report_keys = {
        "cpu_busy_mean_pct", "cpu_busy_min_pct", "cpu_tctl_max_c", "dimm_coverage",
        "dimm_max_c", "edac_ce_end", "edac_ce_start", "edac_ue_end", "edac_ue_start",
        "gate", "max_consecutive_dimm_misses", "max_gap_s", "samples", "schema",
        "uptime_end", "uptime_start", "wall_elapsed_s",
    }
    if not isinstance(report, dict) or set(report) != report_keys or \
            report.get("schema") != "wirehair.wh2.holdout.environment.v1" or \
            report.get("gate") != "PASS":
        fail(f"{panel}: watchdog report schema/gate mismatch")
    numeric_report = ("cpu_busy_mean_pct", "cpu_busy_min_pct", "cpu_tctl_max_c",
                      "max_gap_s", "uptime_end", "uptime_start", "wall_elapsed_s")
    if any(type(report.get(key)) not in (int, float) or
           not math.isfinite(float(report[key])) for key in numeric_report):
        fail(f"{panel}: watchdog numeric report field is malformed")
    if (type(start.get("uptime_s")) not in (int, float) or
            type(finish.get("uptime_s")) not in (int, float) or
            not math.isfinite(float(start["uptime_s"])) or
            not math.isfinite(float(finish["uptime_s"]))):
        fail(f"{panel}: run-record uptime is malformed")
    monotonic = [float(fields[1]) for fields, _ in samples]
    busy = [float(fields[2]) for fields, _ in samples]
    gaps = [right - left for left, right in zip(monotonic, monotonic[1:])]
    if (any(gap <= 0 or gap > 2.5 for gap in gaps) or min(busy) < 90.0 or
            sum(busy) / len(busy) < 98.0):
        fail(f"{panel}: independently reduced gap/load gates failed")
    if (abs(monotonic[0] - float(report["uptime_start"])) > 2.5 or
            abs(float(report["uptime_end"]) - monotonic[-1]) > 2.5 or
            monotonic[0] - float(start["uptime_s"]) < 0 or
            monotonic[0] - float(start["uptime_s"]) > 2.5 or
            float(finish["uptime_s"]) - monotonic[-1] < 0 or
            float(finish["uptime_s"]) - monotonic[-1] > 2.5):
        fail(f"{panel}: independently reduced endpoint gates failed")
    coverage: list[float] = []
    dimm_max: list[float | None] = []
    missing_streak = [0] * 8
    max_missing = [0] * 8
    hot_streak = 0
    max_hot_streak = 0
    ce_values: set[int] = set()
    ue_values: set[int] = set()
    for fields, numbers in samples:
        dims = numbers[5:13]
        missing = sum(value is None for value in dims)
        if missing == 8 or contract.uint(fields[13], "thermal read errors", 8) != missing:
            fail(f"{panel}: independently reduced DIMM read accounting failed")
        for index, value in enumerate(dims):
            missing_streak[index] = missing_streak[index] + 1 if value is None else 0
            max_missing[index] = max(max_missing[index], missing_streak[index])
            if missing_streak[index] > 5:
                fail(f"{panel}: independently reduced consecutive DIMM-miss gate failed")
        temperatures = [numbers[4]] + [value for value in dims if value is not None]
        hot_streak = hot_streak + 1 if any(value is not None and value >= 90.0 for value in temperatures) else 0
        max_hot_streak = max(max_hot_streak, hot_streak)
        ce_values.add(contract.uint(fields[17], "thermal EDAC CE"))
        ue_values.add(contract.uint(fields[18], "thermal EDAC UE"))
    if max_hot_streak >= 3 or len(ce_values) != 1 or len(ue_values) != 1:
        fail(f"{panel}: independently reduced temperature/EDAC gate failed")
    for index in range(8):
        values = [numbers[5 + index] for _, numbers in samples if numbers[5 + index] is not None]
        coverage.append(len(values) / len(samples))
        dimm_max.append(max(values) if values else None)
    if any(value < 0.99 for value in coverage):
        fail(f"{panel}: independently reduced DIMM coverage gate failed")
    ce = next(iter(ce_values))
    ue = next(iter(ue_values))
    recomputed = dict(report)
    recomputed.update({
        "cpu_busy_mean_pct": sum(busy) / len(busy),
        "cpu_busy_min_pct": min(busy),
        "cpu_tctl_max_c": max(float(fields[4]) for fields, _ in samples),
        "dimm_coverage": coverage,
        "dimm_max_c": dimm_max,
        "edac_ce_end": ce, "edac_ce_start": ce, "edac_ue_end": ue, "edac_ue_start": ue,
        "max_consecutive_dimm_misses": max_missing,
        "max_gap_s": max(gaps), "samples": len(samples),
    })
    if recomputed != report:
        fail(f"{panel}: watchdog report differs from independent interval reduction")
    if (start.get("edac_ce", ce) != ce or start.get("edac_ue", ue) != ue or
            finish.get("edac_ce", ce) != ce or finish.get("edac_ue", ue) != ue):
        fail(f"{panel}: run-record EDAC values differ from thermal interval")
    if (float(report["wall_elapsed_s"]) < 0 or
            abs((float(report["uptime_end"]) - float(report["uptime_start"])) -
                float(report["wall_elapsed_s"])) > 1.0):
        fail(f"{panel}: watchdog wall/uptime duration mismatch")
    return report


def audit_environment(args: argparse.Namespace) -> None:
    report = audit_environment_data(args.root, args.panel)
    print(f"environment audit: PASS ({args.panel}, samples={report['samples']})")


def exact_tail(successes: int, trials: int) -> float:
    if trials == 0:
        return 1.0
    numerator = sum(math.comb(trials, value) for value in range(successes, trials + 1))
    return numerator / (2 ** trials)


def paired_p_values(repairs: int, introductions: int) -> tuple[float, float, float]:
    discordant = repairs + introductions
    improve = exact_tail(repairs, discordant)
    regress = exact_tail(introductions, discordant)
    two_sided = min(1.0, 2.0 * min(improve, regress))
    return improve, regress, two_sided


def holm_adjust(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda index: (values[index], index))
    adjusted = [1.0] * len(values)
    running = 0.0
    size = len(values)
    for rank, index in enumerate(order):
        running = max(running, min(1.0, (size - rank) * values[index]))
        adjusted[index] = running
    return adjusted


def scope_names(table: dict[int, dict[str, str]]) -> list[str]:
    result = ["changed"]
    result.extend(f"seed:{index}" for index in range(3))
    result.extend(f"schedule:{value}" for value in contract.SCHEDULES)
    result.extend(f"loss:{value}" for value in contract.LOSSES)
    result.extend(f"schedule_loss:{schedule}:{loss}"
                  for schedule in contract.SCHEDULES for loss in contract.LOSSES)
    result.extend(f"band:{index:02d}" for index in range(16))
    result.extend(("selected_arm:c27", "selected_arm:c79", "selected_arm:c6f"))
    if len(result) != 57:
        fail("internal Holm scope family is not exactly 57")
    return result


def scopes_for(k: int, job: dict[str, object], table: dict[int, dict[str, str]]) -> tuple[str, ...]:
    scopes = [f"seed:{job['seed_index']}", f"schedule:{job['schedule']}",
              f"loss:{job['loss']}", f"schedule_loss:{job['schedule']}:{job['loss']}",
              f"band:{(k - 2) // 4000:02d}"]
    if table[k]["changed"] == "1":
        scopes.append("changed")
        selected = table[k]["v5_salt"][2:]
        scopes.append(f"selected_arm:c{selected}")
    return tuple(scopes)


def empty_metric() -> dict[str, object]:
    return {"trials": 0, "v4_fail": 0, "v5_fail": 0, "repairs": 0,
            "introductions": 0, "overlap": 0, "v4_errors": 0, "v5_errors": 0,
            "paired_success_trials": 0, "v4_xor": Decimal(0), "v5_xor": Decimal(0),
            "v4_mul": Decimal(0), "v5_mul": Decimal(0)}


def add_metric(metric: dict[str, object], v4: dict[str, object], v5: dict[str, object]) -> None:
    trials = int(v4["_trials"])
    a = set(v4["_failures"])
    b = set(v5["_failures"])
    metric["trials"] = int(metric["trials"]) + trials
    metric["v4_fail"] = int(metric["v4_fail"]) + len(a)
    metric["v5_fail"] = int(metric["v5_fail"]) + len(b)
    metric["repairs"] = int(metric["repairs"]) + len(a - b)
    metric["introductions"] = int(metric["introductions"]) + len(b - a)
    metric["overlap"] = int(metric["overlap"]) + len(a & b)
    metric["v4_errors"] = int(metric["v4_errors"]) + int(v4["_errors"])
    metric["v5_errors"] = int(metric["v5_errors"]) + int(v5["_errors"])
    if not a and not b:
        metric["paired_success_trials"] = int(metric["paired_success_trials"]) + trials
        metric["v4_xor"] = Decimal(metric["v4_xor"]) + int(v4["_xor_total"])
        metric["v5_xor"] = Decimal(metric["v5_xor"]) + int(v5["_xor_total"])
        metric["v4_mul"] = Decimal(metric["v4_mul"]) + int(v4["_mul_total"])
        metric["v5_mul"] = Decimal(metric["v5_mul"]) + int(v5["_mul_total"])


def material_regression(metric: dict[str, object]) -> bool:
    delta = int(metric["v5_fail"]) - int(metric["v4_fail"])
    if delta < 3 or int(metric["trials"]) == 0:
        return False
    ratio_material = (int(metric["v4_fail"]) == 0 or
                      10 * int(metric["v5_fail"]) >= 11 * int(metric["v4_fail"]))
    risk_material = 2000 * delta >= int(metric["trials"])
    return ratio_material or risk_material


def ratio(numerator: Decimal, denominator: Decimal) -> str:
    if denominator == 0:
        return "1.000000000" if numerator == 0 else "inf"
    return f"{numerator / denominator:.9f}"


def ratio_value(numerator: Decimal, denominator: Decimal) -> Decimal:
    if denominator == 0:
        return Decimal(1) if numerator == 0 else Decimal("Infinity")
    return numerator / denominator


def p_text(value: float) -> str:
    return f"{value:.17g}"


def analyze(args: argparse.Namespace) -> None:
    jobs, table = load_jobs(args.root)
    verify_artifact_filename_sets(args.root, jobs, sealed=True)
    expected_manifest = make_artifact_manifest(args.root, jobs)
    if contract.read_relative(args.root, "recovery/artifact_manifest.sha256") != expected_manifest:
        fail("analysis inputs differ from sealed recovery artifact manifest")
    names = scope_names(table)
    scopes = {name: empty_metric() for name in names}
    per_k = {k: empty_metric() for k in range(2, 64001)}
    overall = empty_metric()
    semantic_comparisons = 0
    semantic_mismatch_examples: list[dict[str, object]] = []
    semantic_mismatch_count = 0
    for offset in range(0, len(jobs), 2):
        left, right = jobs[offset:offset + 2]
        if left["arm"] != "v4" or right["arm"] != "v5" or any(
                left[key] != right[key] for key in contract.RECOVERY_JOB_HEADER
                if key not in ("job_id", "stem", "arm", "profile")):
            fail(f"jobs {offset}/{offset + 1} are not an exact v4/v5 pair")
        v4 = validate_output_data(
            contract.read_relative(args.root, output_paths(args.root, left)[0]),
            contract.read_relative(args.root, output_paths(args.root, left)[1]), left, table,
            str(left["stem"]))
        v5 = validate_output_data(
            contract.read_relative(args.root, output_paths(args.root, right)[0]),
            contract.read_relative(args.root, output_paths(args.root, right)[1]), right, table,
            str(right["stem"]))
        for k in left["_ks"]:
            a, b = v4[k], v5[k]
            add_metric(overall, a, b)
            add_metric(per_k[k], a, b)
            for name in scopes_for(k, left, table):
                add_metric(scopes[name], a, b)
            if table[k]["changed"] == "0":
                semantic_comparisons += 1
                differing = [field for field in SEMANTIC_FIELDS if a[field] != b[field]]
                if differing:
                    semantic_mismatch_count += 1
                    if len(semantic_mismatch_examples) < 100:
                        semantic_mismatch_examples.append(
                            {"K": k, "job": offset, "fields": differing})
    if int(overall["trials"]) != contract.EXPECTED_PAIRED_ROWS:
        fail("analysis did not consume the exact paired-row geometry")
    scope_regression = [paired_p_values(int(scopes[name]["repairs"]),
                                        int(scopes[name]["introductions"]))[1] for name in names]
    scope_adjusted = holm_adjust(scope_regression)
    k_names = [k for k, row in table.items() if row["changed"] == "1"]
    if len(k_names) != 98:
        fail("changed-K Holm family is not exactly 98")
    k_regression = [paired_p_values(int(per_k[k]["repairs"]),
                                    int(per_k[k]["introductions"]))[1] for k in k_names]
    k_adjusted = dict(zip(k_names, holm_adjust(k_regression)))
    recovery_header = ("scope", "trials", "v4_fail", "v5_fail", "repairs",
                       "introductions", "overlap", "v4_errors", "v5_errors",
                       "p_improvement", "p_regression", "p_two_sided",
                       "holm_regression", "material_regression", "gate")
    recovery_rows = []
    regression_failures = 0
    for index, name in enumerate(names):
        metric = scopes[name]
        improve, regress, two = paired_p_values(int(metric["repairs"]), int(metric["introductions"]))
        material = material_regression(metric)
        gate = "FAIL" if material and scope_adjusted[index] < 0.05 else "PASS"
        regression_failures += gate == "FAIL"
        recovery_rows.append((name, metric["trials"], metric["v4_fail"], metric["v5_fail"],
                              metric["repairs"], metric["introductions"], metric["overlap"],
                              metric["v4_errors"], metric["v5_errors"], p_text(improve),
                              p_text(regress), p_text(two), p_text(scope_adjusted[index]),
                              int(material), gate))
    per_k_header = ("K", "changed", "v5_salt", "trials", "v4_fail", "v5_fail",
                    "repairs", "introductions", "overlap", "v4_errors", "v5_errors",
                    "p_improvement", "p_regression", "p_two_sided",
                    "holm_changed_k_regression", "material_regression", "gate")
    per_k_rows = []
    for k in range(2, 64001):
        metric = per_k[k]
        improve, regress, two = paired_p_values(int(metric["repairs"]), int(metric["introductions"]))
        material = material_regression(metric)
        adjusted = k_adjusted.get(k, 1.0)
        gate = "FAIL" if k in k_adjusted and material and adjusted < 0.05 else "PASS"
        regression_failures += gate == "FAIL"
        per_k_rows.append((k, table[k]["changed"], table[k]["v5_salt"], metric["trials"],
                           metric["v4_fail"], metric["v5_fail"], metric["repairs"],
                           metric["introductions"], metric["overlap"], metric["v4_errors"],
                           metric["v5_errors"], p_text(improve), p_text(regress), p_text(two),
                           p_text(adjusted), int(material), gate))
    work_header = ("scope", "paired_success_trials", "v4_xor", "v5_xor", "xor_ratio",
                   "v4_muladd", "v5_muladd", "muladd_ratio", "gate")
    work_rows = []
    for name, metric in (("overall", overall), *( (name, scopes[name]) for name in names )):
        xor_ratio = ratio(Decimal(metric["v5_xor"]), Decimal(metric["v4_xor"]))
        mul_ratio = ratio(Decimal(metric["v5_mul"]), Decimal(metric["v4_mul"]))
        gate = "PASS"
        if name == "overall" and (
                ratio_value(Decimal(metric["v5_xor"]), Decimal(metric["v4_xor"])) >
                Decimal("1.01") or
                ratio_value(Decimal(metric["v5_mul"]), Decimal(metric["v4_mul"])) >
                Decimal("1.0025")):
            gate = "FAIL"
            regression_failures += 1
        work_rows.append((name, metric["paired_success_trials"], metric["v4_xor"],
                          metric["v5_xor"], xor_ratio, metric["v4_mul"], metric["v5_mul"],
                          mul_ratio, gate))
    changed = scopes["changed"]
    improvement_p = paired_p_values(int(changed["repairs"]), int(changed["introductions"]))[0]
    inconclusive = int(changed["v4_fail"]) < 40
    efficacy = (not inconclusive and 4 * int(changed["v5_fail"]) <= 3 * int(changed["v4_fail"])
                and improvement_p < 0.01)
    hard_failure = (int(overall["v5_fail"]) > int(overall["v4_fail"]) or
                    int(overall["v4_errors"]) != 0 or int(overall["v5_errors"]) != 0 or
                    semantic_mismatch_count != 0 or regression_failures != 0 or
                    (not inconclusive and not efficacy))
    final_gate = "FAIL" if hard_failure else ("INCONCLUSIVE" if inconclusive else "PASS")
    semantic = {
        "comparisons": semantic_comparisons,
        "gate": "PASS" if semantic_mismatch_count == 0 else "FAIL",
        "mismatch_count": semantic_mismatch_count,
        "mismatches": semantic_mismatch_examples,
        "schema": "wirehair.wh2.holdout.semantic_identity.v1",
        "unchanged_k": 63901,
    }
    environment_report = audit_environment_data(args.root, "recovery")
    summary = {
        "changed": {key: changed[key] for key in
                    ("trials", "v4_fail", "v5_fail", "repairs", "introductions", "overlap")},
        "changed_efficacy_p": p_text(improvement_p),
        "changed_failure_reduction_fraction": (
            None if int(changed["v4_fail"]) == 0 else
            (int(changed["v4_fail"]) - int(changed["v5_fail"])) / int(changed["v4_fail"])),
        "gate": final_gate,
        "holm_changed_k": 98,
        "holm_scopes": 57,
        "inconclusive_low_changed_failure_support": inconclusive,
        "overall": {key: overall[key] for key in
                    ("trials", "v4_fail", "v5_fail", "repairs", "introductions", "overlap",
                     "v4_errors", "v5_errors")},
        "regression_gate_failures": regression_failures,
        "schema": "wirehair.wh2.normalized_h15_v5.holdout.summary.v1",
        "semantic_mismatches": semantic_mismatch_count,
    }
    outputs = {
        "analysis/recovery_scopes.csv": contract.canonical_tsv(recovery_header, recovery_rows),
        "analysis/per_k_recovery.csv": contract.canonical_tsv(per_k_header, per_k_rows),
        "analysis/work_scopes.csv": contract.canonical_tsv(work_header, work_rows),
        "analysis/semantic_identity.json": contract.canonical_json(semantic),
        "analysis/environment.json": contract.canonical_json(environment_report),
        "analysis/summary.json": contract.canonical_json(summary),
    }
    if getattr(args, "audit_existing", False):
        expected_names = {Path(relative).name for relative in outputs}
        contract.verify_exact_directory_entries(args.root, "analysis", expected_names)
        for relative, data in outputs.items():
            if contract.read_relative(args.root, relative) != data:
                fail(f"{relative}: differs from independent exact recomputation")
        verify_artifact_filename_sets(args.root, jobs, sealed=True)
        contract.verify_exact_directory_entries(args.root, "analysis", expected_names)
        print(f"analysis audit: PASS gate={final_gate}")
    else:
        for relative, data in outputs.items():
            contract.write_exclusive(args.root, relative, data)
        verify_artifact_filename_sets(args.root, jobs, sealed=True)
        contract.verify_exact_directory_entries(
            args.root, "analysis", (Path(relative).name for relative in outputs))
        print(f"analysis gate={final_gate} v4_fail={overall['v4_fail']} v5_fail={overall['v5_fail']}")


def audit_analysis(args: argparse.Namespace) -> None:
    args.audit_existing = True
    analyze(args)


def expect_failure(action: Callable[[], object], label: str) -> None:
    try:
        action()
    except (AuditError, OSError, UnicodeError, ValueError, csv.Error):
        return
    fail(f"selftest mutation unexpectedly passed: {label}")


def selftest(_args: argparse.Namespace) -> None:
    if paired_p_values(3, 0) != (0.125, 1.0, 0.25):
        fail("selftest exact binomial tails mismatch")
    if recover_integer_total("0.062", 16, "selftest") != 1:
        fail("selftest failed to recover an exact T16 integer work total")
    expect_failure(lambda: recover_integer_total("0.063", 16, "selftest"),
                   "wrong half-even T16 rounding")
    adjusted = holm_adjust([0.01, 0.04, 0.03])
    if adjusted != [0.03, 0.06, 0.06]:
        fail(f"selftest Holm adjustment mismatch: {adjusted}")
    metric = empty_metric()
    metric.update({"trials": 10000, "v4_fail": 0, "v5_fail": 5,
                   "repairs": 0, "introductions": 5})
    if not material_regression(metric):
        fail("selftest material regression was missed")
    metric["introductions"] = metric["v5_fail"] = 2
    if material_regression(metric):
        fail("selftest sparse nonmaterial delta was misclassified")
    header = (",".join(THERMAL_HEADER) + "\r\n").encode("ascii")
    if header.replace(b"\r\n", b"\n") != (",".join(THERMAL_HEADER) + "\n").encode("ascii"):
        fail("selftest thermal CRLF normalization mismatch")
    expect_failure(lambda: parse_thermal_line(b"bad\n", 2), "bad thermal field count")
    job: dict[str, object] = {
        "_ks": (2,), "_trials": 1, "arm": "v4", "loss": "0.10",
        "profile": "normalized-h15-v4", "schedule": "iid",
        "seed": "0x123456789abcdef0", "trials": "1",
    }
    meta = expected_preamble(job)
    meta["packet_peel_seed_xor"] = "0x0"
    raw_values = {
        "N": "2", "bb": "64", "heavy_family": "periodic", "mix_count": "2",
        "overhead": "0", "trials": "1", "success": "1", "rank_fail": "0",
        "error": "0", "fail_rate": "0.00000000", "inact_mu": "0.000",
        "inact_max": "0", "binary_def_mu": "0.000", "binary_def_max": "0",
        "heavy_gain_mu": "0.000", "heavy_gain_min": "0", "heavy_shortfall": "0",
        "solve_ms_mu": "0.000", "build_ms_mu": "0.000", "peel_ms_mu": "0.000",
        "project_ms_mu": "0.000", "residual_ms_mu": "0.000",
        "backsub_ms_mu": "0.000", "seed_attempt": "1", "block_xors_mu": "0.000",
        "block_muladds_mu": "0.000", "first_rank_fail": "-1",
        "binary_def_hist": "0:1", "heavy_gain_hist": "0:1", "failure_trials": "",
        "active_packet_peel_seed_xor": "0x0",
        "mixed_joint_source_xors_mu": "1.000",
        "mixed_joint_marginal_xors_mu": "2.000",
        "mixed_joint_marginal_copies_mu": "3.000",
        "mixed_joint_active_deltas_mu": "4.000",
        "mixed_joint_scratch_bytes_mu": "5.000",
        "mixed_dual_source_columns_mu": "6.000",
    }
    def synthetic_raw(values: dict[str, str], fields: tuple[str, ...] = RAW_HEADER) -> bytes:
        output = io.StringIO(newline="")
        output.write("# precodefail: " + " ".join(
            f"{key}={value}" for key, value in meta.items()) + "\n")
        writer = csv.writer(output, lineterminator="\n")
        writer.writerow(fields)
        writer.writerow([values[field] for field in fields])
        return output.getvalue().encode("ascii")
    parsed = validate_output_data(synthetic_raw(raw_values), b"", job, None, "selftest raw")
    if tuple(parsed) != (2,):
        fail("selftest integrated raw schema was not accepted")
    bad_joint = dict(raw_values)
    bad_joint["mixed_joint_source_xors_mu"] = "0.500"
    expect_failure(lambda: validate_output_data(
        synthetic_raw(bad_joint), b"", job, None, "selftest raw"),
        "nonintegral joint-delta T1 mean")
    expect_failure(lambda: validate_output_data(
        synthetic_raw(raw_values, RAW_HEADER[:-1]), b"", job, None, "selftest raw"),
        "missing integrated joint-delta column")
    with tempfile.TemporaryDirectory(prefix="wh2-recovery-set-selftest-") as text:
        root = Path(text) / "private"
        contract.ensure_private_root(root, fresh=True)
        paths = [
            "recovery/raw/job.csv", "recovery/status/job.ok",
            "recovery/commands/job.txt", "recovery/run_start.json",
            "recovery/run_finish.json", "recovery/thermal_interval.csv",
        ]
        for path in paths:
            contract.write_exclusive(root, path, b"fixture\n")
        for path in (
                "recovery/holdout_reveal.tsv", "recovery/reveal_provenance.json",
                "recovery/manifest.json", "recovery/jobs.tsv", "recovery/input_seal.sha256",
                "recovery/artifact_manifest.sha256"):
            contract.write_exclusive(root, path, b"fixture\n")
        verify_recovery_filename_paths(root, paths, sealed=True)
        contract.write_exclusive(root, "recovery/status/unmanifested.tmp", b"extra\n")
        expect_failure(lambda: verify_recovery_filename_paths(root, paths, sealed=True),
                       "post-recovery-seal unmanifested extra file")
        (root / "recovery/status/unmanifested.tmp").unlink()
        verify_recovery_filename_paths(root, paths, sealed=True)
    print("selftest: ok (tails, Holm, thermal, joint deltas, exact artifact sets)")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name, function in (("audit-inputs", audit_inputs), ("audit-artifacts", audit_artifacts),
                           ("seal-artifacts", seal_artifacts), ("verify-artifacts", verify_artifacts),
                           ("analyze", analyze), ("audit-analysis", audit_analysis)):
        command = commands.add_parser(name)
        command.add_argument("--root", type=Path, required=True)
        command.set_defaults(func=function)
    output = commands.add_parser("validate-output")
    output.add_argument("--root", type=Path, required=True)
    output.add_argument("--job-id", required=True)
    output.add_argument("--stdout", type=Path, required=True)
    output.add_argument("--stderr", type=Path, required=True)
    output.set_defaults(func=validate_output)
    temporary = commands.add_parser("validate-temp")
    temporary.add_argument("--job-row", required=True)
    temporary.add_argument("--stdout", type=Path, required=True)
    temporary.add_argument("--stderr", type=Path, required=True)
    temporary.set_defaults(func=validate_temp)
    job = commands.add_parser("validate-job")
    job.add_argument("--root", type=Path, required=True)
    job.add_argument("--job-id", required=True)
    job.set_defaults(func=validate_job)
    watch = commands.add_parser("watchdog")
    watch.add_argument("--thermal-source", type=Path, required=True)
    watch.add_argument("--start-line", type=int, required=True)
    watch.add_argument("--done-file", type=Path, required=True)
    watch.add_argument("--output", type=Path, required=True)
    watch.add_argument("--report", type=Path, required=True)
    watch.set_defaults(func=watchdog)
    record = commands.add_parser("run-record")
    record.add_argument("--root", type=Path, required=True)
    record.add_argument("--output", type=Path, required=True)
    record.add_argument("--kind", choices=("start", "finish"), required=True)
    record.add_argument("--workers", type=int)
    record.add_argument("--thermal-source", type=Path)
    record.add_argument("--thermal-line", type=int)
    record.add_argument("--launcher-nice", type=int)
    record.add_argument("--xargs-rc", type=int)
    record.add_argument("--watchdog-report", type=Path)
    record.set_defaults(func=run_record)
    environment = commands.add_parser("audit-environment")
    environment.add_argument("--root", type=Path, required=True)
    environment.add_argument("--panel", choices=("recovery", "payload", "lookup"), required=True)
    environment.set_defaults(func=audit_environment)
    test = commands.add_parser("selftest")
    test.set_defaults(func=selftest)
    args = parser.parse_args()
    try:
        args.func(args)
    except (AuditError, OSError, UnicodeError, ValueError, csv.Error) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
