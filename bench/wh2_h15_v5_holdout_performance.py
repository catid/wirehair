#!/usr/bin/env python3
"""Run and reduce the frozen H15-v5 payload, lookup, and static-size panels."""

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
import statistics
import subprocess
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
LOOKUP_RE = re.compile(
    r"profile=full-domain-permutation "
    r"arm=(normalized-h15-v[45]) "
    r"lookups=131069952 permutation_step=40501 domain_repeats=2048 "
    r"checksum=(0x[0-9a-f]+) elapsed_ns=([1-9][0-9]*) "
    r"ns_per_lookup=([0-9]+(?:\.[0-9]+)?) "
    r"selection_table_sha256=([0-9a-f]{64})\n"
)
AuditError = contract.AuditError
fail = contract.fail

# This is the shared integrated candidate/table-free-reference precodefail
# interface.  Holdout timing executes the frozen candidate for both table
# arms; the frozen reference is intentionally used only for static-size
# deltas, but its matching header is checked during the pre-freeze rehearsal.

def decimal(text: str, label: str) -> Decimal:
    try:
        value = Decimal(text)
    except InvalidOperation:
        fail(f"{label}: invalid decimal {text!r}")
    if not value.is_finite() or value < 0:
        fail(f"{label}: expected finite nonnegative decimal")
    return value


def integer_mean_t1(text: str, label: str) -> Decimal:
    if not re.fullmatch(r"[0-9]+\.[0-9]{3}", text):
        fail(f"{label}: T1 integer mean must use exact three-decimal formatting")
    value = decimal(text, label)
    if value != value.to_integral_value():
        fail(f"{label}: T1 integer mean is not integral")
    return value


def preamble(line: str, label: str) -> dict[str, str]:
    prefix = "# precodefail: "
    if not line.startswith(prefix):
        fail(f"{label}: missing precodefail preamble")
    result = {}
    for item in line[len(prefix):].split():
        if "=" not in item:
            fail(f"{label}: malformed preamble token")
        key, value = item.split("=", 1)
        if not key or key in result:
            fail(f"{label}: duplicate preamble key")
        result[key] = value
    return result


def parse_failures(text: str, trials: int, label: str) -> tuple[int, ...]:
    if not text:
        return ()
    result = tuple(contract.uint(value, label, trials - 1) for value in text.split("|"))
    if any(a >= b for a, b in zip(result, result[1:])):
        fail(f"{label}: failure trial IDs are duplicate/unsorted")
    return result


def expected_payload_preamble(job: dict[str, str]) -> dict[str, str]:
    return {
        "trials": "1", "threads": "1", "loss": "0.10000000000000001",
        "seed": f"0x{int(job['seed'], 16):x}", "completion": "mixed",
        "mixed_period": "32", "mixed_gf256_rows": "11", "mixed_gf16_rows": "4",
        "mixed_geometry": "shared-x", "mixed_residue_skew": "0",
        "mixed_residue_schedule": "hashed", "mixed_residue_hash_seed": "0x44",
        "mixed_residue_hash_keyed": "1", "mixed_independent_extension_residues": "1",
        "mixed_residue_buckets_requested": "auto",
        "mixed_extension_residue_seed_xor": "0x4e", "source_hits_override": "0",
        "packet_peel_seed_xor": "0x0", "packet_peel_seed_table": job["profile"],
        "binary_dense_rows_override": "0", "gf256_heavy_rows_override": "0",
        "odd_packet_peel_seed_xor": "0x0", "packet_row_seed_multiplier": "0x1",
        "packet_row_seed_avalanche": "0", "seed_block_bytes_override": "1280",
        "overhead_stream": "salted", "full_payload_solve": "1", "schedule": "repair-only",
    }


def parse_payload(data: bytes, stderr: bytes, job: dict[str, str],
                  table: dict[int, dict[str, str]], label: str) -> dict[int, dict[str, object]]:
    if stderr:
        fail(f"{label}: stderr must be empty")
    if b"\r" in data or not data.endswith(b"\n") or b"\n\n" in data:
        fail(f"{label}: noncanonical framing")
    text = data.decode("ascii")
    stream = io.StringIO(text, newline="")
    meta = preamble(stream.readline().rstrip("\n"), label)
    expected = expected_payload_preamble(job)
    if meta != expected:
        differing = sorted(key for key in set(meta) | set(expected) if meta.get(key) != expected.get(key))
        fail(f"{label}: payload preamble mismatch fields={differing}")
    reader = csv.DictReader(stream)
    if tuple(reader.fieldnames or ()) != RAW_HEADER:
        fail(f"{label}: raw header mismatch")
    expected_k = tuple(map(int, job["k_csv"].split(",")))
    rows = list(reader)
    if len(rows) != len(expected_k):
        fail(f"{label}: expected exactly one payload row per K")
    result: dict[int, dict[str, object]] = {}
    for line, (row, wanted) in enumerate(zip(rows, expected_k), 3):
        if None in row or any(value is None for value in row.values()):
            fail(f"{label}:{line}: malformed row")
        k = contract.uint(row["N"], f"{label}:{line}:K", 64000)
        trials = contract.uint(row["trials"], f"{label}:{line}:trials", 1)
        success = contract.uint(row["success"], f"{label}:{line}:success", 1)
        rank_fail = contract.uint(row["rank_fail"], f"{label}:{line}:rank_fail", 1)
        errors = contract.uint(row["error"], f"{label}:{line}:error", 1)
        failures = parse_failures(row["failure_trials"], trials, f"{label}:{line}:failures")
        expected_salt = table[k]["v4_salt" if job["arm"] == "v4" else "v5_salt"]
        if (k != wanted or row["bb"] != job["width"] or row["heavy_family"] != "periodic" or
                row["mix_count"] != "2" or row["overhead"] != "4" or trials != 1 or
                success + rank_fail + errors != 1 or len(failures) != rank_fail + errors or
                row["active_packet_peel_seed_xor"] != expected_salt):
            fail(f"{label}:{line}: payload geometry/outcome/salt mismatch")
        for field in JOINT_METRIC_FIELDS:
            integer_mean_t1(row[field], f"{label}:{line}:{field}")
        result[k] = {"failures": failures, "rank_fail": rank_fail, "errors": errors,
                     "solve": decimal(row["solve_ms_mu"], f"{label}:{line}:solve"),
                     "xor": integer_mean_t1(row["block_xors_mu"], f"{label}:{line}:xor"),
                     "mul": integer_mean_t1(row["block_muladds_mu"], f"{label}:{line}:mul")}
    if tuple(result) != expected_k:
        fail(f"{label}: payload K order/content mismatch")
    return result


def parse_lookup(data: bytes, stderr: bytes, job: dict[str, str], label: str) -> dict[str, object]:
    if stderr:
        fail(f"{label}: stderr must be empty")
    match = LOOKUP_RE.fullmatch(data.decode("ascii"))
    if not match:
        fail(f"{label}: lookup stdout contract mismatch")
    arm, checksum, elapsed, ns_per, table_sha = match.groups()
    if (arm != job["profile"] or checksum != job["expected_checksum"] or
            table_sha != contract.TABLE_SHA256):
        fail(f"{label}: lookup arm/checksum/table identity mismatch")
    elapsed_value = int(elapsed)
    ns_value = decimal(ns_per, f"{label}:ns_per_lookup")
    expected_ns = Decimal(elapsed_value) / Decimal(job["lookups"])
    if abs(ns_value - expected_ns) > Decimal("0.001"):
        fail(f"{label}: printed lookup time disagrees with elapsed/count")
    return {"elapsed_ns": elapsed_value, "ns": ns_value, "checksum": checksum}


def load_contract(root: Path) -> tuple[dict[str, object], dict[int, dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    manifest = contract.validate_freeze(root)
    if not isinstance(manifest, dict) or manifest.get("schema") != contract.CONTRACT_SCHEMA:
        fail("freeze contract schema mismatch")
    binary = contract.read_relative(root, "freeze/wirehair_v2_bench")
    reference = contract.read_relative(root, "freeze/reference_wirehair_v2_bench")
    if (contract.sha256_bytes(binary) != manifest.get("binary_sha256") or
            contract.sha256_bytes(reference) != manifest.get("reference_binary_sha256")):
        fail("candidate/reference binary differs from freeze contract")
    table = contract.parse_table(contract.read_relative(root, "freeze/selection_table.tsv"))
    payload_data = contract.read_relative(root, "freeze/payload_k.tsv")
    if contract.sha256_bytes(payload_data) != contract.PAYLOAD_K_SHA256:
        fail("payload panel SHA mismatch")
    expected_payload = contract.make_payload_jobs(payload_data)
    actual_payload = contract.read_relative(root, "freeze/payload_jobs.tsv")
    if actual_payload != expected_payload:
        fail("payload jobs differ from exact reconstruction")
    expected_lookup = contract.make_lookup_jobs()
    actual_lookup = contract.read_relative(root, "freeze/lookup_jobs.tsv")
    if actual_lookup != expected_lookup:
        fail("lookup jobs differ from exact reconstruction")
    payload_jobs = contract.parse_tsv(actual_payload, contract.PAYLOAD_JOB_HEADER, "payload jobs")
    lookup_jobs = contract.parse_tsv(actual_lookup, contract.LOOKUP_JOB_HEADER, "lookup jobs")
    return manifest, table, payload_jobs, lookup_jobs


def panel_paths(panel: str, job: dict[str, str]) -> tuple[str, ...]:
    stem = job["stem"]
    return (f"performance/{panel}/raw/{stem}.stdout",
            f"performance/{panel}/raw/{stem}.stderr",
            f"performance/{panel}/raw/{stem}.time",
            f"performance/{panel}/status/{stem}.ok",
            f"performance/{panel}/commands/{stem}.txt")


def terminate_and_drain(process: subprocess.Popen[bytes], timeout: float = 5.0) -> None:
    if process.poll() is None:
        try:
            process.terminate()
        except ProcessLookupError:
            pass
    try:
        process.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        if process.poll() is None:
            try:
                process.kill()
            except ProcessLookupError:
                pass
        process.communicate()


def finish_watchdog(watchdog: subprocess.Popen[bytes], done: Path,
                    graceful_timeout: float = 10.0,
                    terminate_timeout: float = 5.0) -> int:
    marker_error: BaseException | None = None
    wait_error: BaseException | None = None
    cleanup_error: BaseException | None = None
    try:
        done.write_bytes(b"done\n")
    except BaseException as exc:
        marker_error = exc
    if marker_error is None and watchdog.poll() is None:
        try:
            watchdog.wait(timeout=graceful_timeout)
        except subprocess.TimeoutExpired:
            pass
        except BaseException as exc:
            wait_error = exc
    if watchdog.poll() is None:
        try:
            terminate_and_drain(watchdog, terminate_timeout)
        except BaseException as exc:
            cleanup_error = exc
    # A local child must be reaped even when marker creation, graceful wait,
    # or termination itself failed.  SIGKILL plus an unbounded final wait is
    # the last-resort no-leak guarantee.
    if watchdog.poll() is None:
        try:
            watchdog.kill()
        except ProcessLookupError:
            pass
        watchdog.wait()
    if marker_error is not None:
        raise marker_error
    if wait_error is not None:
        raise wait_error
    if cleanup_error is not None:
        raise cleanup_error
    if watchdog.returncode is None:
        fail("environment watchdog was not reaped")
    return watchdog.returncode


def monitored_communicate(process: subprocess.Popen[bytes],
                          watchdog: subprocess.Popen[bytes]) -> tuple[bytes, bytes]:
    while True:
        if watchdog.poll() is not None:
            terminate_and_drain(process)
            fail("environment watchdog ended while a performance job was running")
        try:
            stdout, stderr = process.communicate(timeout=0.1)
        except subprocess.TimeoutExpired:
            # communicate() drains both pipes while it waits and is safe to
            # retry after TimeoutExpired; the eventual return contains every
            # byte accumulated across all retries.
            continue
        if watchdog.poll() is not None:
            fail("environment watchdog ended before performance output was accepted")
        return stdout, stderr


def run_child(command: list[str], watchdog: subprocess.Popen[bytes]) -> tuple[bytes, bytes, float]:
    started = time.perf_counter()
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    try:
        stdout, stderr = monitored_communicate(process, watchdog)
    except BaseException:
        terminate_and_drain(process)
        raise
    elapsed = time.perf_counter() - started
    if process.returncode != 0:
        fail(f"performance command exited {process.returncode}: {shlex.join(command)}")
    return stdout, stderr, elapsed


def thermal_line_count(path: Path) -> int:
    flags = os.O_RDONLY | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    fd = os.open(path, flags)
    try:
        count = 0
        while True:
            data = os.read(fd, 1024 * 1024)
            if not data:
                return count
            count += data.count(b"\n")
    finally:
        os.close(fd)


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


def write_panel_manifest(root: Path, panel: str, jobs: list[dict[str, str]], count: int) -> bytes:
    paths = [path for job in jobs for path in panel_paths(panel, job)]
    paths.extend((f"performance/{panel}/run_start.json", f"performance/{panel}/run_finish.json",
                  f"performance/{panel}/thermal_interval.csv"))
    if len(paths) != count:
        fail(f"{panel}: internal artifact count mismatch")
    return contract.canonical_tsv(("sha256", "path"), (
        (contract.sha256_bytes(contract.read_relative(root, path)), path) for path in paths
    ))


def panel_command(root: Path, panel: str, job: dict[str, str]) -> list[str]:
    binary = str(root / "freeze/wirehair_v2_bench")
    if panel == "payload":
        return [binary, "precodefail", "--N", job["k_csv"], "--bb-list", job["width"],
                "--seed-block-bytes", "1280", "--overhead", "4", "--trials", "1",
                "--threads", "1", "--loss", "0.10", "--seed", job["seed"],
                "--schedule", "repair-only", "--completion", "mixed", "--mix-count", "2",
                "--packet-peel-seed-table", job["profile"], "--full-payload-solve",
                "--mixed-gf256-rows", "11", "--mixed-gf16-rows", "4", "--mixed-period", "32",
                "--mixed-geometry", "shared-x", "--mixed-residue-schedule", "hashed",
                "--mixed-residue-hash-seed", "68", "--mixed-residue-hash-keyed",
                "--mixed-independent-extension-residues",
                "--mixed-extension-residue-seed-xor", "78"]
    if panel == "lookup":
        return ["taskset", "-c", job["cpu"], binary, "h15table",
                "--permutation-step", job["permutation_step"],
                "--domain-repeats", job["domain_repeats"], "--table", job["profile"]]
    fail(f"unknown performance panel {panel!r}")


def run_panel(root: Path, panel: str, jobs: list[dict[str, str]],
              table: dict[int, dict[str, str]], thermal_source: Path) -> dict[int, object]:
    panel_root = root / "performance" / panel
    if panel_root.exists() or panel_root.is_symlink():
        fail(f"refusing existing performance panel directory: {panel_root}")
    for relative in (f"performance/{panel}/raw", f"performance/{panel}/status",
                     f"performance/{panel}/commands"):
        contract.ensure_dir(root, relative)
    temp = Path(tempfile.mkdtemp(prefix=f".{panel}-", dir=root / "performance"))
    os.chmod(temp, 0o700)
    done = temp / "done"
    thermal = temp / "thermal.csv"
    report = temp / "thermal.json"
    analyzer = root / "freeze/tools/wh2_h15_v5_holdout_analyze.py"
    start_line = thermal_line_count(thermal_source)
    if start_line < 2:
        fail("live thermal source has no complete sample")
    start_uptime = read_uptime()
    now = dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")
    results: dict[int, object] = {}
    active_error: BaseException | None = None
    watchdog_rc: int | None = None
    watchdog: subprocess.Popen[bytes] | None = None
    try:
        watchdog = subprocess.Popen((sys.executable, str(analyzer), "watchdog",
                                     "--thermal-source", str(thermal_source),
                                     "--start-line", str(start_line),
                                     "--done-file", str(done), "--output", str(thermal),
                                     "--report", str(report)))
        start_record = {"jobs": len(jobs), "panel": panel,
                        "schema": "wirehair.wh2.holdout.performance_run.v1",
                        "thermal_line_start": start_line, "uptime_s": start_uptime, "utc": now}
        contract.write_exclusive(
            root, f"performance/{panel}/run_start.json", contract.canonical_json(start_record))
        for job in jobs:
            command = panel_command(root, panel, job)
            stdout, stderr, elapsed = run_child(command, watchdog)
            paths = panel_paths(panel, job)
            if panel == "payload":
                results[int(job["job_id"])] = parse_payload(stdout, stderr, job, table, job["stem"])
            else:
                results[int(job["job_id"])] = parse_lookup(stdout, stderr, job, job["stem"])
            artifacts = (stdout, stderr, f"elapsed_s={elapsed:.9f} exit=0\n".encode("ascii"),
                         b"ok\n", (shlex.join(command) + "\n").encode("ascii"))
            for path, data in zip(paths, artifacts):
                contract.write_exclusive(root, path, data)
    except BaseException as exc:
        active_error = exc
        raise
    finally:
        if watchdog is not None:
            try:
                watchdog_rc = finish_watchdog(watchdog, done)
            except BaseException:
                if active_error is None:
                    raise
    if watchdog_rc != 0:
        fail(f"{panel}: live environment watchdog failed")
    thermal_data = contract.read_once(thermal)
    watchdog_data = contract.parse_json(contract.read_once(report), f"{panel} watchdog")
    contract.write_exclusive(root, f"performance/{panel}/thermal_interval.csv", thermal_data)
    finish_record = {"environment": watchdog_data, "jobs": len(jobs), "panel": panel,
                     "schema": "wirehair.wh2.holdout.performance_run.v1",
                     "uptime_s": read_uptime(),
                     "utc": dt.datetime.now(dt.timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")}
    contract.write_exclusive(root, f"performance/{panel}/run_finish.json", contract.canonical_json(finish_record))
    for entry in temp.iterdir():
        entry.unlink()
    temp.rmdir()
    return results


def ratio(numerator: Decimal, denominator: Decimal) -> Decimal:
    if denominator == 0:
        return Decimal(1) if numerator == 0 else Decimal("Infinity")
    return numerator / denominator


def reduce_payload(jobs: list[dict[str, str]], results: dict[int, object],
                   table: dict[int, dict[str, str]],
                   enforce_panel_hash: bool = True) -> dict[str, object]:
    paired: dict[tuple[int, int], dict[str, tuple[dict[str, str], dict[int, dict[str, object]]]]] = defaultdict(dict)
    for job in jobs:
        paired[(int(job["phase"]), int(job["width"]))][job["arm"]] = (job, results[int(job["job_id"])])  # type: ignore[assignment]
    totals = {"all": [Decimal(0), Decimal(0)], "changed": [Decimal(0), Decimal(0)]}
    widths = {1280: [Decimal(0), Decimal(0)], 4096: [Decimal(0), Decimal(0)]}
    work = [Decimal(0), Decimal(0), Decimal(0), Decimal(0)]
    phase_solve = {phase: [Decimal(0), Decimal(0)] for phase in range(8)}
    errors = 0
    v5_only_rank = 0
    v4_only_repairs = 0
    included = 0
    excluded = 0
    cell_coverage: dict[str, dict[str, int]] = {}
    panel_rows = contract.parse_tsv(
        contract.make_payload_k(table, enforce_hash=enforce_panel_hash),
        contract.PAYLOAD_K_HEADER, "payload classes")
    classes = {int(row["K"]): row["class"] for row in panel_rows}
    class_coverage = {name: {"included": 0, "excluded": 0} for name in
                      ("changed-v5", "immutable-v4", "structural-control", "maximin-control")}
    for (phase, width), arms in paired.items():
        if set(arms) != {"v4", "v5"}:
            fail("payload phase/width lacks an exact arm pair")
        v4 = arms["v4"][1]
        v5 = arms["v5"][1]
        cell_included = 0
        for k in v4:
            a, b = v4[k], v5[k]
            errors += int(a["errors"]) + int(b["errors"])
            v5_only_rank += int(bool(b["rank_fail"]) and not bool(a["rank_fail"]))
            v4_only_repairs += int(bool(a["rank_fail"]) and not bool(b["rank_fail"]))
            both_success = not a["failures"] and not b["failures"]
            if not both_success:
                excluded += 1
                class_coverage[classes[k]]["excluded"] += 1
                continue
            included += 1
            cell_included += 1
            class_coverage[classes[k]]["included"] += 1
            totals["all"][0] += Decimal(a["solve"])
            totals["all"][1] += Decimal(b["solve"])
            widths[width][0] += Decimal(a["solve"])
            widths[width][1] += Decimal(b["solve"])
            phase_solve[phase][0] += Decimal(a["solve"])
            phase_solve[phase][1] += Decimal(b["solve"])
            if table[k]["changed"] == "1":
                totals["changed"][0] += Decimal(a["solve"])
                totals["changed"][1] += Decimal(b["solve"])
            work[0] += Decimal(a["xor"])
            work[1] += Decimal(b["xor"])
            work[2] += Decimal(a["mul"])
            work[3] += Decimal(b["mul"])
        cell_coverage[f"phase{phase}:width{width}"] = {
            "excluded": 256 - cell_included, "included": cell_included}
    all_ratio = ratio(totals["all"][1], totals["all"][0])
    changed_ratio = ratio(totals["changed"][1], totals["changed"][0])
    width_ratios = {str(width): ratio(value[1], value[0]) for width, value in widths.items()}
    phase_ratios = [ratio(value[1], value[0]) for value in phase_solve.values()]
    median_phase = statistics.median(phase_ratios)
    xor_ratio = ratio(work[1], work[0])
    mul_ratio = ratio(work[3], work[2])
    coverage_valid = all(value["included"] >= 231 for value in cell_coverage.values())
    hard_valid = errors == 0 and v5_only_rank == 0
    speed_pass = (all_ratio <= Decimal("1.01") and changed_ratio <= Decimal("1.02") and
                  all(value <= Decimal("1.02") for value in width_ratios.values()) and
                  median_phase <= Decimal("1.01") and xor_ratio <= Decimal("1.01") and
                  mul_ratio <= Decimal("1.0025"))
    gate = ("FAIL" if not hard_valid else
            "INCONCLUSIVE" if not coverage_valid else
            "PASS" if speed_pass else "FAIL")
    return {"all_256_solve_ratio": str(all_ratio), "changed_98_solve_ratio": str(changed_ratio),
            "class_coverage": class_coverage, "coverage_by_phase_width": cell_coverage,
            "coverage_valid": coverage_valid, "errors": errors, "excluded_observations": excluded,
            "gate": gate, "included_observations": included,
            "median_phase_ratio": str(median_phase), "muladd_ratio": str(mul_ratio),
            "phase_ratios": [str(value) for value in phase_ratios],
            "schema": "wirehair.wh2.holdout.payload_summary.v1",
            "v4_only_repairs": v4_only_repairs, "v5_only_rank_failures": v5_only_rank,
            "width_ratios": {k: str(v) for k, v in width_ratios.items()},
            "xor_ratio": str(xor_ratio)}


def reduce_lookup(jobs: list[dict[str, str]], results: dict[int, object]) -> dict[str, object]:
    phases: dict[int, dict[str, Decimal]] = defaultdict(dict)
    for job in jobs:
        value = results[int(job["job_id"])]  # type: ignore[assignment]
        # Reduce the integer timer reading, not the CLI's three-decimal display.
        # This keeps threshold decisions independent of presentation rounding.
        phases[int(job["phase"])][job["arm"]] = (
            Decimal(value["elapsed_ns"]) / Decimal(job["lookups"]))
    ratios = []
    deltas = []
    for phase in range(12):
        if set(phases[phase]) != {"v4", "v5"}:
            fail("lookup phase lacks an exact v4/v5 pair")
        ratios.append(ratio(phases[phase]["v5"], phases[phase]["v4"]))
        deltas.append(phases[phase]["v5"] - phases[phase]["v4"])
    median_ratio = statistics.median(ratios)
    median_delta = statistics.median(deltas)
    passed = (median_ratio <= Decimal("1.25") and median_delta <= Decimal("2") and
              max(ratios) <= Decimal("1.50"))
    return {"signed_median_delta_ns": str(median_delta), "checksums": contract.LOOKUP_CHECKSUMS,
            "gate": "PASS" if passed else "FAIL", "median_ratio": str(median_ratio),
            "per_phase_ratios": [str(value) for value in ratios],
            "schema": contract.LOOKUP_SUMMARY_SCHEMA}


def size_sections(path: Path) -> dict[str, int]:
    result = subprocess.run(("size", str(path)), check=False, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        fail(f"size failed for {path}: {result.stderr.strip()}")
    lines = result.stdout.splitlines()
    if len(lines) != 2:
        fail(f"size output shape mismatch for {path}")
    fields = lines[1].split()
    if (len(fields) < 6 or not fields[0].isdigit() or not fields[1].isdigit() or
            not fields[2].isdigit()):
        fail(f"size output fields mismatch for {path}")
    return {"text": int(fields[0]), "data": int(fields[1]), "bss": int(fields[2])}


def lookup_symbol_sizes(path: Path) -> tuple[int, int, dict[str, int]]:
    result = subprocess.run(("nm", "-S", "--size-sort", "-C", str(path)), check=False,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        fail(f"nm failed for {path}: {result.stderr.strip()}")
    # Code size is compiler/target dependent (the frozen native and portable
    # builds legitimately differ), so only the logical arrays have exact byte
    # sizes.  The two code symbols must remain separately emitted, non-empty
    # text symbols; this is the binary-level evidence that the benchmark entry
    # point did not disappear into an inline caller.
    wanted = {
        "helper": ("NormalizedH15V5PacketPeelSeedXor(unsigned int) [clone .part.0]", None, "text"),
        "benchmark_wrapper": ("BenchmarkNormalizedH15V5PacketPeelSeedXor(unsigned int)", None, "text"),
        "salts": ("kNormalizedH15V5Salts", 168, "rodata"),
        "offsets": ("kNormalizedH15V5BucketOffsets", 252, "rodata"),
        "keys": ("kNormalizedH15V5Keys", 336, "rodata"),
    }
    found: dict[str, int] = {}
    addresses: dict[str, int] = {}
    for line in result.stdout.splitlines():
        fields = line.split(maxsplit=3)
        if (len(fields) != 4 or not re.fullmatch(r"[0-9a-fA-F]+", fields[0]) or
                not re.fullmatch(r"[0-9a-fA-F]+", fields[1]) or len(fields[2]) != 1):
            continue
        name = fields[3]
        for role, (suffix, expected, section) in wanted.items():
            if name.endswith(suffix):
                if role in found:
                    fail(f"nm reported duplicate H15-v5 symbol role {role}")
                size = int(fields[1], 16)
                kind = fields[2]
                if size <= 0:
                    fail(f"H15-v5 symbol {role} is empty")
                if section == "text" and kind not in ("t", "T"):
                    fail(f"H15-v5 symbol {role} is not an emitted text symbol")
                if section == "rodata" and kind not in ("r", "R"):
                    fail(f"H15-v5 symbol {role} is not read-only data")
                if expected is not None and size != expected:
                    fail(f"H15-v5 symbol {role} size {size} differs from frozen {expected}")
                found[role] = size
                addresses[role] = int(fields[0], 16)
    if set(found) != set(wanted):
        fail(f"nm H15-v5 symbol set mismatch: found {sorted(found)}")
    if addresses["helper"] == addresses["benchmark_wrapper"]:
        fail("H15-v5 helper and benchmark wrapper were folded into one symbol")
    for role in ("helper", "benchmark_wrapper"):
        start = addresses[role]
        result = subprocess.run(("objdump", "-d", "-C", f"--start-address={start}",
                                 f"--stop-address={start + found[role]}", str(path)),
                                check=False, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True)
        if result.returncode != 0:
            fail(f"objdump failed for H15-v5 {role}: {result.stderr.strip()}")
        if not re.search(r"(?m)^\s*[0-9a-fA-F]+:\s+(?:[0-9a-fA-F]{2}\s+)+", result.stdout):
            fail(f"objdump found no instructions for H15-v5 {role}")
    logical = found["salts"] + found["offsets"] + found["keys"]
    return logical, sum(found.values()), found


def static_summary(root: Path) -> dict[str, object]:
    candidate = size_sections(root / "freeze/wirehair_v2_bench")
    reference = size_sections(root / "freeze/reference_wirehair_v2_bench")
    deltas = {key: candidate[key] - reference[key] for key in ("text", "data", "bss")}
    logical, lookup_symbols, symbols = lookup_symbol_sizes(root / "freeze/wirehair_v2_bench")
    # Whole-binary growth includes the untimed h15table dump/validation CLI and
    # is sealed as provenance.  The hard code/data cap applies to the exact
    # lookup-specific symbol set above, i.e. the production hot path.
    passed = logical <= 2048 and lookup_symbols <= 4096
    return {"candidate_sections": candidate, "gate": "PASS" if passed else "FAIL",
            "logical_table_bytes": logical,
            "lookup_symbol_bytes": lookup_symbols, "lookup_symbols": symbols,
            "reference_sections": reference,
            "schema": "wirehair.wh2.holdout.static_size.v1",
            "whole_binary_section_deltas": deltas}


def audit_static(args: argparse.Namespace) -> None:
    root = Path(os.path.abspath(args.root))
    contract.ensure_private_root(root, fresh=False)
    contract.validate_freeze(root)
    expected = contract.canonical_json(static_summary(root))
    actual = contract.read_relative(root, "performance/static_size_summary.json")
    if actual != expected:
        fail("static-size summary differs from an independent binary recomputation")
    print("static-size audit: PASS")


def audit_results(args: argparse.Namespace) -> None:
    root = Path(os.path.abspath(args.root))
    contract.ensure_private_root(root, fresh=False)
    _, table, payload_jobs, lookup_jobs = load_contract(root)
    contract.verify_performance_artifact_sets(root)
    for panel, jobs, count in (
            ("payload", payload_jobs, contract.EXPECTED_PAYLOAD_ARTIFACTS),
            ("lookup", lookup_jobs, contract.EXPECTED_LOOKUP_ARTIFACTS)):
        results: dict[int, object] = {}
        for job in jobs:
            paths = panel_paths(panel, job)
            stdout = contract.read_relative(root, paths[0])
            stderr = contract.read_relative(root, paths[1])
            if panel == "payload":
                results[int(job["job_id"])] = parse_payload(
                    stdout, stderr, job, table, paths[0])
            else:
                results[int(job["job_id"])] = parse_lookup(
                    stdout, stderr, job, paths[0])
            timing = contract.read_relative(root, paths[2])
            if not re.fullmatch(rb"elapsed_s=[0-9]+\.[0-9]{9} exit=0\n", timing):
                fail(f"{paths[2]}: malformed elapsed-time artifact")
            if contract.read_relative(root, paths[3]) != b"ok\n":
                fail(f"{paths[3]}: status is not exact ok")
            expected_command = (shlex.join(panel_command(root, panel, job)) + "\n").encode("ascii")
            if contract.read_relative(root, paths[4]) != expected_command:
                fail(f"{paths[4]}: command differs from exact reconstruction")
        expected_manifest = write_panel_manifest(root, panel, jobs, count)
        manifest_path = f"performance/{panel}_artifact_manifest.sha256"
        if contract.read_relative(root, manifest_path) != expected_manifest:
            fail(f"{panel}: artifact manifest differs from exact recomputation")
        expected_summary = (reduce_payload(jobs, results, table) if panel == "payload"
                            else reduce_lookup(jobs, results))
        summary_path = f"performance/{panel}_summary.json"
        if contract.read_relative(root, summary_path) != contract.canonical_json(expected_summary):
            fail(f"{panel}: summary differs from independent raw-result reduction")
    expected_static = contract.canonical_json(static_summary(root))
    if contract.read_relative(root, "performance/static_size_summary.json") != expected_static:
        fail("static-size summary differs from independent binary recomputation")
    contract.verify_performance_artifact_sets(root)
    print("performance result audit: PASS (32 payload + 24 lookup jobs)")


def run(args: argparse.Namespace) -> None:
    root = Path(os.path.abspath(args.root))
    contract.ensure_private_root(root, fresh=False)
    _, table, payload_jobs, lookup_jobs = load_contract(root)
    payload_results = run_panel(root, "payload", payload_jobs, table, args.thermal_source)
    payload_summary = reduce_payload(payload_jobs, payload_results, table)
    contract.write_exclusive(root, "performance/payload_summary.json", contract.canonical_json(payload_summary))
    payload_manifest = write_panel_manifest(root, "payload", payload_jobs, contract.EXPECTED_PAYLOAD_ARTIFACTS)
    contract.write_exclusive(root, "performance/payload_artifact_manifest.sha256", payload_manifest)
    contract.verify_hash_manifest(root, "performance/payload_artifact_manifest.sha256",
                                  contract.EXPECTED_PAYLOAD_ARTIFACTS)
    if payload_summary["gate"] != "PASS":
        fail("payload performance panel failed its frozen gates")
    lookup_results = run_panel(root, "lookup", lookup_jobs, table, args.thermal_source)
    lookup_summary = reduce_lookup(lookup_jobs, lookup_results)
    contract.write_exclusive(root, "performance/lookup_summary.json", contract.canonical_json(lookup_summary))
    lookup_manifest = write_panel_manifest(root, "lookup", lookup_jobs, contract.EXPECTED_LOOKUP_ARTIFACTS)
    contract.write_exclusive(root, "performance/lookup_artifact_manifest.sha256", lookup_manifest)
    contract.verify_hash_manifest(root, "performance/lookup_artifact_manifest.sha256",
                                  contract.EXPECTED_LOOKUP_ARTIFACTS)
    static = static_summary(root)
    contract.write_exclusive(root, "performance/static_size_summary.json", contract.canonical_json(static))
    if lookup_summary["gate"] != "PASS" or static["gate"] != "PASS":
        fail("lookup or static-size panel failed its frozen gates")
    audit_results(argparse.Namespace(root=root))
    print("performance gates: PASS (payload=32 jobs, lookup=24 jobs)")


def expect_failure(action: Callable[[], object], label: str) -> None:
    try:
        action()
    except (AuditError, OSError, UnicodeError, ValueError, csv.Error):
        return
    fail(f"selftest mutation unexpectedly passed: {label}")


def selftest(_args: argparse.Namespace) -> None:
    live_watchdog = subprocess.Popen(
        (sys.executable, "-c", "import time; time.sleep(30)"),
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    large_stdout = 2 * 1024 * 1024 + 17
    large_stderr = 1024 * 1024 + 31
    try:
        stdout, stderr, _ = run_child(
            [sys.executable, "-c",
             "import sys;"
             f"sys.stdout.buffer.write(b'x'*{large_stdout});"
             f"sys.stderr.buffer.write(b'y'*{large_stderr})"],
            live_watchdog)
    finally:
        if live_watchdog.poll() is None:
            live_watchdog.terminate()
        live_watchdog.wait(timeout=5)
    if stdout != b"x" * large_stdout or stderr != b"y" * large_stderr:
        fail("selftest monitored pipe drain lost or duplicated child output")
    dead_watchdog = subprocess.Popen(
        (sys.executable, "-c", "pass"),
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    dead_watchdog.wait(timeout=5)
    sleeping_child = subprocess.Popen(
        (sys.executable, "-c", "import time; time.sleep(30)"),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    expect_failure(lambda: monitored_communicate(sleeping_child, dead_watchdog),
                   "dead watchdog child cleanup")
    if sleeping_child.poll() is None:
        terminate_and_drain(sleeping_child)
        fail("selftest watchdog failure left its performance child running")
    with tempfile.TemporaryDirectory(prefix="wh2-watchdog-cleanup-") as text:
        cleanup_root = Path(text)
        preloop_watchdog = subprocess.Popen(
            (sys.executable, "-c", "import time; time.sleep(30)"),
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        active_error: BaseException | None = None
        observed_error = ""
        try:
            try:
                raise AuditError("forced failure before performance job loop")
            except BaseException as exc:
                active_error = exc
                raise
            finally:
                try:
                    finish_watchdog(
                        preloop_watchdog, cleanup_root / "missing/done", 0.1, 0.1)
                except BaseException:
                    if active_error is None:
                        raise
        except AuditError as exc:
            observed_error = str(exc)
        if (observed_error != "forced failure before performance job loop" or
                preloop_watchdog.poll() is None):
            if preloop_watchdog.poll() is None:
                terminate_and_drain(preloop_watchdog, 0.1)
            fail("selftest pre-loop failure did not preserve its error and reap watchdog")
        panel_root = cleanup_root / "actual-panel"
        panel_root.mkdir(mode=0o700)
        (panel_root / "performance").mkdir(mode=0o700)
        thermal_source = cleanup_root / "thermal-source.csv"
        thermal_source.write_bytes(b"header\nsample\n")
        original_popen = subprocess.Popen
        original_write = contract.write_exclusive
        spawned_watchdogs: list[subprocess.Popen[bytes]] = []
        def tracked_popen(command: object, *positional: object,
                          **keywords: object) -> subprocess.Popen[bytes]:
            command_fields = tuple(command)  # type: ignore[arg-type]
            done_path = command_fields[command_fields.index("--done-file") + 1]
            process = original_popen(
                (sys.executable, "-c",
                 "import pathlib,sys,time;"
                 "p=pathlib.Path(sys.argv[1]);"
                 "\nwhile not p.exists(): time.sleep(0.01)",
                 str(done_path)),
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            spawned_watchdogs.append(process)
            return process
        def fail_run_start(root: Path, relative: str, data: bytes,
                           executable: bool = False) -> None:
            if relative == "performance/payload/run_start.json":
                raise AuditError("forced actual run_panel pre-loop failure")
            original_write(root, relative, data, executable)
        subprocess.Popen = tracked_popen  # type: ignore[assignment]
        contract.write_exclusive = fail_run_start
        actual_panel_error = ""
        try:
            try:
                run_panel(panel_root, "payload", [], {}, thermal_source)
            except AuditError as exc:
                actual_panel_error = str(exc)
        finally:
            contract.write_exclusive = original_write
            subprocess.Popen = original_popen  # type: ignore[assignment]
        if (actual_panel_error != "forced actual run_panel pre-loop failure" or
                len(spawned_watchdogs) != 1 or spawned_watchdogs[0].poll() is None):
            for process in spawned_watchdogs:
                if process.poll() is None:
                    terminate_and_drain(process, 0.1)
            fail("selftest actual run_panel did not preserve failure and reap watchdog")
        stubborn_ready = cleanup_root / "stubborn.ready"
        stubborn_watchdog = subprocess.Popen(
            (sys.executable, "-c",
             "import pathlib,signal,sys,time;"
             "signal.signal(signal.SIGTERM, signal.SIG_IGN);"
             "pathlib.Path(sys.argv[1]).write_text('ready\\n');time.sleep(30)",
             str(stubborn_ready)),
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        ready_deadline = time.monotonic() + 5.0
        while time.monotonic() < ready_deadline:
            if stubborn_ready.exists():
                break
            if stubborn_watchdog.poll() is not None:
                break
            time.sleep(0.01)
        if not stubborn_ready.exists():
            terminate_and_drain(stubborn_watchdog, 0.1)
            fail("selftest stubborn watchdog did not reach its ready state")
        stubborn_rc = finish_watchdog(
            stubborn_watchdog, cleanup_root / "stubborn.done", 0.1, 0.1)
        if stubborn_watchdog.poll() is None or stubborn_rc == 0:
            if stubborn_watchdog.poll() is None:
                terminate_and_drain(stubborn_watchdog, 0.1)
            fail("selftest nonterminating watchdog escaped kill/reap fallback")
    def synthetic_lookup_reduction(v5_ns: int) -> dict[str, object]:
        fake_jobs = []
        fake_results: dict[int, object] = {}
        lookups = 131069952
        for phase, first in enumerate(contract.LOOKUP_FIRST_ARMS):
            arms = ("v4", "v5") if first == "v4" else ("v5", "v4")
            for arm in arms:
                job_id = len(fake_jobs)
                profile = f"normalized-h15-{arm}"
                ns = 100 if arm == "v4" else v5_ns
                fake_jobs.append({"job_id": str(job_id), "phase": str(phase), "arm": arm,
                                  "profile": profile, "lookups": str(lookups),
                                  "expected_checksum": contract.LOOKUP_CHECKSUMS[profile]})
                fake_results[job_id] = {
                    "elapsed_ns": lookups * ns, "ns": Decimal(ns)}
        return reduce_lookup(fake_jobs, fake_results)

    lookup_cases = (
        ("+1 ns slowdown", 101, "PASS", "1", "1.01"),
        ("+2 ns slowdown boundary", 102, "PASS", "2", "1.02"),
        ("+3 ns slowdown", 103, "FAIL", "3", "1.03"),
        ("-4 ns improvement", 96, "PASS", "-4", "0.96"),
    )
    for label, v5_ns, expected_gate, expected_delta, expected_ratio in lookup_cases:
        summary = synthetic_lookup_reduction(v5_ns)
        if (summary != {
                "checksums": contract.LOOKUP_CHECKSUMS,
                "gate": expected_gate,
                "median_ratio": expected_ratio,
                "per_phase_ratios": [expected_ratio] * 12,
                "schema": contract.LOOKUP_SUMMARY_SCHEMA,
                "signed_median_delta_ns": expected_delta}):
            fail(f"selftest lookup signed-delta reduction mismatch: {label}")
    if ratio(Decimal(0), Decimal(0)) != 1 or ratio(Decimal(1), Decimal(0)).is_finite():
        fail("selftest ratio edge cases mismatch")
    table_rows = contract.parse_tsv(
        contract.synthetic_table(), contract.TABLE_HEADER, "synthetic payload table")
    table = {int(row["K"]): row for row in table_rows}
    raw_job = {"arm": "v4", "k_csv": "2", "profile": "normalized-h15-v4",
               "seed": "0x123456789abcdef0", "width": "1280"}
    raw_values = {
        "N": "2", "bb": "1280", "heavy_family": "periodic", "mix_count": "2",
        "overhead": "4", "trials": "1", "success": "1", "rank_fail": "0",
        "error": "0", "fail_rate": "0.00000000", "inact_mu": "0.000",
        "inact_max": "0", "binary_def_mu": "0.000", "binary_def_max": "0",
        "heavy_gain_mu": "0.000", "heavy_gain_min": "0", "heavy_shortfall": "0",
        "solve_ms_mu": "1.000", "build_ms_mu": "0.000", "peel_ms_mu": "0.000",
        "project_ms_mu": "0.000", "residual_ms_mu": "0.000",
        "backsub_ms_mu": "0.000", "seed_attempt": "1", "block_xors_mu": "10.000",
        "block_muladds_mu": "20.000", "first_rank_fail": "-1",
        "binary_def_hist": "0:1", "heavy_gain_hist": "0:1", "failure_trials": "",
        "active_packet_peel_seed_xor": table[2]["v4_salt"],
        "mixed_joint_source_xors_mu": "1.000",
        "mixed_joint_marginal_xors_mu": "2.000",
        "mixed_joint_marginal_copies_mu": "3.000",
        "mixed_joint_active_deltas_mu": "4.000",
        "mixed_joint_scratch_bytes_mu": "5.000",
        "mixed_dual_source_columns_mu": "6.000",
    }
    def synthetic_payload(values: dict[str, str],
                          fields: tuple[str, ...] = RAW_HEADER) -> bytes:
        output = io.StringIO(newline="")
        output.write("# precodefail: " + " ".join(
            f"{key}={value}" for key, value in expected_payload_preamble(raw_job).items()) + "\n")
        writer = csv.writer(output, lineterminator="\n")
        writer.writerow(fields)
        writer.writerow([values[field] for field in fields])
        return output.getvalue().encode("ascii")
    if tuple(parse_payload(
            synthetic_payload(raw_values), b"", raw_job, table, "selftest payload")) != (2,):
        fail("selftest integrated payload schema was not accepted")
    bad_joint = dict(raw_values)
    bad_joint["mixed_joint_active_deltas_mu"] = "4.500"
    expect_failure(lambda: parse_payload(
        synthetic_payload(bad_joint), b"", raw_job, table, "selftest payload"),
        "nonintegral payload joint-delta T1 mean")
    expect_failure(lambda: parse_payload(
        synthetic_payload(raw_values, RAW_HEADER[:-1]), b"", raw_job, table,
        "selftest payload"), "missing payload joint-delta column")
    panel = contract.make_payload_k(table, enforce_hash=False)
    payload_jobs = contract.parse_tsv(
        contract.make_payload_jobs(panel), contract.PAYLOAD_JOB_HEADER, "synthetic payload jobs")
    payload_results: dict[int, object] = {}
    for job in payload_jobs:
        arm_solve = Decimal("100") if job["arm"] == "v4" else Decimal("99")
        payload_results[int(job["job_id"])] = {
            k: {"errors": 0, "failures": (), "mul": Decimal(100), "rank_fail": 0,
                "solve": arm_solve, "xor": Decimal(100)}
            for k in map(int, job["k_csv"].split(","))}
    payload_summary = reduce_payload(
        payload_jobs, payload_results, table, enforce_panel_hash=False)
    if (payload_summary["gate"] != "PASS" or
            payload_summary["included_observations"] != 4096 or
            payload_summary["excluded_observations"] != 0):
        fail("selftest paired-success payload reduction mismatch")
    expect_failure(lambda: parse_lookup(b"bad\n", b"", {
        "profile": "normalized-h15-v4", "expected_checksum": contract.LOOKUP_CHECKSUMS["normalized-h15-v4"],
        "lookups": "131069952"}, "bad"), "malformed lookup output")
    print("selftest: ok (watchdog lifecycle, drained pipes, signed lookup gate, reducers, and joint-delta schema)")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    run_parser = commands.add_parser("run")
    run_parser.add_argument("--root", type=Path, required=True)
    run_parser.add_argument("--thermal-source", type=Path,
                            default=Path("/tmp/wirehair-enoq-thermal.csv"))
    run_parser.set_defaults(func=run)
    audit_parser = commands.add_parser("audit-static")
    audit_parser.add_argument("--root", type=Path, required=True)
    audit_parser.set_defaults(func=audit_static)
    results_parser = commands.add_parser("audit-results")
    results_parser.add_argument("--root", type=Path, required=True)
    results_parser.set_defaults(func=audit_results)
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
