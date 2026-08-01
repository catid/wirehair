#!/usr/bin/env python3
"""Sealed H12-versus-H13 grouped-GF(256) Stage-B screen controller.

``prepare`` authenticates a completed Stage-A campaign, independently
replays its exact OH0 census, and derives a deterministic 1:1 matched-control
screen before freezing any Stage-B artifact.  ``run`` executes only that
screen and, on the two prespecified improvement gates, seals an exact all-K
escalation *plan*.  ``reduce`` revalidates the sealed shards and emits a
summary.  It never executes the all-K plan implicitly.

The controller is intentionally strict.  Benchmark text is a receipt, not a
friendly interchange format: field order, spelling, values, pair ordering,
construction attempt, packet trace, and architecture must all match exactly.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from decimal import Decimal
import fcntl
from functools import wraps
import hashlib
import heapq
import json
import math
import os
from pathlib import Path
import queue
import re
import shutil
import signal
import stat
import subprocess
import sys
import tempfile
import threading
import types
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from typing import Any, Iterable, Iterator, Mapping, Optional, Sequence, Tuple


SCHEMA = "wirehair.wh2.h13_stage_b.v2"
SOURCE_SCHEMA = "wirehair.wh2.h13_stage_a.v2.uniform_construction_seed_v1"
K_MIN = 2
K_MAX = 64000
K_COUNT = K_MAX - K_MIN + 1
CHUNK_SIZE = 250
SCREEN_MATCH_CHUNK = 8
ARMS = ("h12", "h13")
ARM_GF256_ROWS = {"h12": 10, "h13": 11}
ARM_HEAVY_ROWS = {"h12": 12, "h13": 13}
ARM_GROUP_HASH_SEED = {"h12": "0xb7e15162", "h13": "0xb7e15163"}
CONSTRUCTION_ROOTS = (
    "0x243f6a8885a308d3",
    "0x13198a2e03707344",
    "0xa4093822299f31d0",
)
LOSS_ROOTS = (
    "0xd1b54a32d192ed03",
    "0x94d049bb133111eb",
    "0x8538ecb5bd456ea3",
)
STAGE_A_ROOT_SETS = {"R0": (0,), "R1": (1, 2)}
SCHEDULES = ("burst", "adversarial", "repair-only")
BANDS = (
    (2, 4095), (4096, 8191), (8192, 9999), (10000, 16383),
    (16384, 32767), (32768, 49151), (49152, 64000),
)
SOURCE_PAIRED_CELLS = (
    K_COUNT * len(CONSTRUCTION_ROOTS) * len(LOSS_ROOTS) * len(SCHEDULES)
)
SOURCE_ARM_OUTCOMES = SOURCE_PAIRED_CELLS * len(ARMS)
FORMAL_STAGE_A_CONTROLLER_SHA256 = (
    "f989dd991de3159abf3cb0a8c42deb118c41b2a868c4759c599d38170e6fa823"
)
FORMAL_STAGE_A_BINARY_SHA256 = (
    "df8da0d908e172f3cfd2ecba367b39d95fa253370f681560d1dd927e39b66a06"
)
FORMAL_STAGE_A_CAMPAIGN_COMPLETE_FILE_SHA256 = (
    "38c3bcbfa392a96a95a08342a02c46d0404a557aff44401c7a0cba37b2237285"
)
FORMAL_STAGE_A_CAMPAIGN_COMPLETE_SEAL = (
    "3dee34a593203c6009d7f5ea3cbb913900eccbeec5ff40a0236d242f1f542cf0"
)
FORMAL_STAGE_A_ANALYSIS_FILE_SHA256 = (
    "61619d0539f72eea3b5f7045edbd11c2b78c3941bcb629bc5e9b092cc158d1a3"
)
FORMAL_STAGE_A_ANALYSIS_SEAL = (
    "ead9727fd0f1684fc9179fc481b5987d2eafe7e896b8eba508d262a4a24c5c34"
)
# The former pinned source above was later found to inherit production per-K
# construction-seed fixups.  It is useful for controller-mechanics fixtures,
# but it is not the required uniform-root experimental source.  Repin every
# formal identity and cohort invariant after that Stage-A rerun before setting
# this true.
FORMAL_STAGE_A_LAUNCHABLE = False
FORMAL_STAGE_A_IDENTITIES_REPINNED = False
FORMAL_COHORT_INVARIANTS_REPINNED = False
FORMAL_STAGE_A_BLOCKER = (
    "fresh uniform-root Stage A v2 has not completed; repin all formal "
    "identities, union/stratum/cohort invariants, then explicitly release "
    "the Stage-B launch interlock"
)
# Retained only as explicit provenance for the superseded v1 experiment.
FORMER_V1_SOURCE_UNION = 893
FORMER_V1_NONEMPTY_STRATA = 61
FORMER_V1_ZERO_STRATA = (
    (0, 0, "repair-only", 2),
    (0, 2, "burst", 2),
)
# Active formal cohort invariants are fail-closed sentinels until the fresh v2
# analysis is reviewed and pinned.  They must all be replaced together with
# FORMAL_COHORT_INVARIANTS_REPINNED=True.
EXPECTED_SOURCE_UNION = 0
EXPECTED_NONEMPTY_STRATA = 0
EXPECTED_ZERO_STRATA: tuple[tuple[int, int, str, int], ...] = ()
EXPECTED_SCREEN_MATCHES = 0
EXPECTED_SCREEN_CELLS = 0
EXPECTED_SCREEN_OUTCOMES = 0
EXPECTED_SCREEN_TASKS = 0
ALLK_TASKS_PER_CONSTRUCTION_ROOT = (
    math.ceil(K_COUNT / CHUNK_SIZE) * len(LOSS_ROOTS) * len(SCHEDULES)
)
ALLK_TASKS = ALLK_TASKS_PER_CONSTRUCTION_ROOT * len(CONSTRUCTION_ROOTS)
ALLK_CELLS_PER_CONSTRUCTION_ROOT = K_COUNT * len(LOSS_ROOTS) * len(SCHEDULES)
ALLK_CELLS = SOURCE_PAIRED_CELLS
ALLK_OUTCOMES = SOURCE_ARM_OUTCOMES
OUTPUT_LIMIT = 8 << 20
EXECUTION_MODE = "authenticated-fd-taskset-v1"
MAX_I63 = (1 << 63) - 1
MASK64 = (1 << 64) - 1
SELECTION_RULE = (
    "within each construction-root/packet-trace-root/schedule/K-band "
    "stratum, sort raw OH0 "
    "H12+H13 common-success non-union cells by SHA256 of "
    "wirehair.wh2.h13_stage_b.v2|matched-control|construction_seed_index|"
    "construction_seed|packet_trace_root_index|packet_trace_root|schedule|"
    "band_lo|band_hi|K, then pair the "
    "first m with the m union cells in ascending K order"
)
LOSS_SEED_FORMULA = (
    "(external_seed ^ (K * 0x9e3779b97f4a7c15) ^ "
    "(64 * 0xbf58476d1ce4e5b9)) mod 2^64; paired overhead salt = 0"
)
THERMAL_FIELDS = (
    "utc", "monotonic_s", "cpu_busy_pct", "cpu_avg_mhz", "cpu_tctl_c",
    "dimm_i2c1_50_c", "dimm_i2c1_51_c", "dimm_i2c1_52_c",
    "dimm_i2c1_53_c", "dimm_i2c2_50_c", "dimm_i2c2_51_c",
    "dimm_i2c2_52_c", "dimm_i2c2_53_c", "dimm_read_errors",
    "load1", "load5", "load15", "edac_ce", "edac_ue",
)
THERMAL_DIMM_FIELDS = THERMAL_FIELDS[5:13]
TELEMETRY_POLICY = {
    "max_gap_seconds": 5.0,
    "min_cpu_busy_pct": 95.0,
    "max_temperature_c": 90.0,
    "require_zero_dimm_read_errors": True,
    "require_constant_edac_counters": True,
}
TELEMETRY_CONTINUITY = (
    "single append-only interval; external sampler and independent >=95% "
    "load remain continuously active from run start through terminal sealing; "
    "a low row, >5s gap, rotation, temperature ceiling, DIMM read error, or "
    "EDAC counter change invalidates the result directory"
)


def architecture_contract() -> dict[str, Any]:
    """Exact fixed semantics emitted by the groupedrecovery preamble/rows."""
    return {
        "native_schema": "v2", "pair_order": list(ARMS), "arms_per_N": 2,
        "policy": "h12-h13-q0-grouped-v1", "overhead": 0,
        "period": 48, "geometry": "shared-x", "residue_skew": 0,
        "residue_schedule": "constant", "residue_hash_seed": "0x0",
        "extension_residue_seed_xor": "0x4e",
        "independent_extension_residues": False,
        "gf256_rows": ARM_GF256_ROWS, "gf16_rows": 2, "grouped_rows": 3,
        "grouped_gf256_row_mask": "0x380", "buckets": "separate",
        "grouped_hash_seed": ARM_GROUP_HASH_SEED,
        "grouped_final_h_a_columns": "arm-heavy-rows",
        "final_h_a_columns": ARM_HEAVY_ROWS,
        "dense_rows": 12, "dense_identity_corner": False,
        "dense_two_anchor": True, "dense_two_anchor_phase": 0,
        "source_hits": "canonical-K-rule",
        "staircase_rows": "GetDenseCount(K)",
        "field": "mixed-gf256-gf16", "heavy_family": "periodic-cauchy",
        "construction_attempt": 0, "systematic_probe": "direct-attempt0",
        "construction_seed_policy":
            "matrix-c-peel-lo32-xor-hi32-v1",
        "production_seed_fixups_applied": 0,
        "construction_roots": list(CONSTRUCTION_ROOTS),
        "construction_root_sets": {
            name: list(indices) for name, indices in STAGE_A_ROOT_SETS.items()
        },
        "external_seed_role": "loss-trace-root",
        "packet_trace_roots": list(LOSS_ROOTS),
        "mix_count": 2, "bb": 64,
        "solve_block_bytes": 2, "loss": "0.5",
        "overhead_stream": "paired", "packet_row_seed_multiplier": "0x1",
        "packet_row_seed_avalanche": 0, "odd_packet_peel_seed_xor": "0x0",
        "packet_vector": "one-shared-per-N", "payload": "shared-zero-v1",
        "packet_payload_bytes": 2,
        "packet_trace_schema":
            "wirehair-wh2-precodefail-raw-packet-trace-v1",
        "coefficient_layout": "materialized-before-solve",
        "grouped_schedule_prefix": "materialized-before-solve",
    }


class CampaignError(RuntimeError):
    """A fail-closed campaign contract violation."""


class CampaignInterrupted(CampaignError):
    """A controlled SIGINT/SIGTERM shutdown."""

    def __init__(self, signum: int) -> None:
        self.signum = signum
        super().__init__(f"campaign interrupted by signal {signum}")


def die(message: str) -> None:
    raise CampaignError(message)


def require_formal_stage_a_launchable() -> None:
    if (FORMAL_STAGE_A_LAUNCHABLE is not True or
            FORMAL_STAGE_A_IDENTITIES_REPINNED is not True or
            FORMAL_COHORT_INVARIANTS_REPINNED is not True):
        die(f"Stage-B launch is blocked: {FORMAL_STAGE_A_BLOCKER}")


def canonical_json(value: Any) -> str:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False,
    )


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def stable_source_bytes(path: Path, *, require_unique: bool = True) -> bytes:
    flags = (os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
             getattr(os, "O_NOFOLLOW", 0))
    try:
        descriptor = os.open(path, flags)
    except OSError as exc:
        die(f"cannot open sealed regular file {path}: {exc}")
    try:
        before = os.fstat(descriptor)
        if (not stat.S_ISREG(before.st_mode) or
                (require_unique and before.st_nlink != 1)):
            die(f"sealed input is not a unique regular file: {path}")
        first: list[bytes] = []
        while True:
            block = os.read(descriptor, 1 << 20)
            if not block:
                break
            first.append(block)
        middle = os.fstat(descriptor)
        os.lseek(descriptor, 0, os.SEEK_SET)
        second: list[bytes] = []
        while True:
            block = os.read(descriptor, 1 << 20)
            if not block:
                break
            second.append(block)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    try:
        named = os.stat(path, follow_symlinks=False)
    except OSError as exc:
        die(f"sealed input pathname disappeared: {path}: {exc}")
    identity = lambda value: (
        value.st_dev, value.st_ino, value.st_mode, value.st_nlink,
        value.st_size, value.st_mtime_ns, value.st_ctime_ns,
    )
    if (identity(before) != identity(middle) or
            identity(middle) != identity(after) or
            identity(after) != identity(named) or
            not stat.S_ISREG(named.st_mode) or
            (require_unique and named.st_nlink != 1)):
        die(f"sealed input changed while reading: {path}")
    data = b"".join(first)
    if data != b"".join(second) or len(data) != before.st_size:
        die(f"sealed input was not stable: {path}")
    return data


def stable_bytes(path: Path) -> bytes:
    return stable_source_bytes(path, require_unique=True)


def sha256_file(path: Path) -> str:
    return sha256_bytes(stable_bytes(path))


def sealed_record(schema: str, payload: Any) -> dict[str, Any]:
    body = {"schema": schema, "payload": payload}
    return {**body, "seal_sha256": sha256_bytes(canonical_json(body).encode())}


def verify_sealed(record: Any, schema: str) -> Any:
    if not isinstance(record, dict) or set(record) != {
            "schema", "payload", "seal_sha256"}:
        die(f"{schema}: malformed sealed record")
    if record["schema"] != schema:
        die(f"sealed schema {record['schema']!r}, want {schema!r}")
    seal = record["seal_sha256"]
    if not isinstance(seal, str) or not re.fullmatch(r"[0-9a-f]{64}", seal):
        die(f"{schema}: malformed seal")
    body = {"schema": schema, "payload": record["payload"]}
    if sha256_bytes(canonical_json(body).encode()) != seal:
        die(f"{schema}: seal mismatch")
    return record["payload"]


def load_sealed(path: Path, schema: str) -> Any:
    data = stable_bytes(path)
    try:
        record = json.loads(data.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        die(f"{path}: invalid sealed JSON: {exc}")
    if data != (canonical_json(record) + "\n").encode():
        die(f"{path}: JSON is not canonical")
    return verify_sealed(record, schema)


def require_exact_seal(
    path: Path, schema: str, expected_seal: str,
    expected_file_sha256: Optional[str] = None,
) -> Any:
    data = stable_bytes(path)
    try:
        record = json.loads(data.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        die(f"{path}: invalid pinned sealed JSON: {exc}")
    if (expected_file_sha256 is not None and
            sha256_bytes(data) != expected_file_sha256):
        die(f"{path}: file SHA256 is not the pinned formal Stage-A artifact")
    if (data != (canonical_json(record) + "\n").encode() or
            not isinstance(record, dict) or
            record.get("seal_sha256") != expected_seal):
        die(f"{path}: record is not the pinned formal Stage-A seal")
    return verify_sealed(record, schema)


def atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", dir=path.parent,
    )
    temporary_path = Path(temporary)
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(data)
            output.flush()
            os.fsync(output.fileno())
        os.replace(temporary_path, path)
    finally:
        try:
            temporary_path.unlink()
        except FileNotFoundError:
            pass


def write_sealed_once(path: Path, schema: str, payload: Any) -> None:
    data = (canonical_json(sealed_record(schema, payload)) + "\n").encode()
    if path.exists():
        if stable_bytes(path) != data:
            die(f"refusing to replace a different sealed record: {path}")
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.create-", dir=path.parent,
    )
    temporary_path = Path(temporary)
    try:
        with os.fdopen(descriptor, "wb") as output:
            output.write(data)
            output.flush()
            os.fsync(output.fileno())
        try:
            os.link(temporary_path, path, follow_symlinks=False)
        except FileExistsError:
            if stable_bytes(path) != data:
                die(f"refusing to replace a different sealed record: {path}")
    finally:
        try:
            temporary_path.unlink()
        except FileNotFoundError:
            pass


def strict_uint(text: str, context: str, maximum: Optional[int] = None) -> int:
    if (len(text) > 20 or
            not re.fullmatch(r"0|[1-9][0-9]*", text)):
        die(f"{context}: noncanonical unsigned integer {text!r}")
    value = int(text)
    limit = MASK64 if maximum is None else maximum
    if value > limit:
        die(f"{context}: integer exceeds {limit}")
    return value


def strict_hex(text: str, context: str, bits: int) -> int:
    if (type(bits) is not int or bits < 1 or
            len(text) > 2 + (bits + 3) // 4 or
            not re.fullmatch(r"0x(?:0|[1-9a-f][0-9a-f]*)", text)):
        die(f"{context}: noncanonical hexadecimal integer {text!r}")
    value = int(text, 16)
    if value >= 1 << bits:
        die(f"{context}: integer does not fit u{bits}")
    return value


def checked_product(*values: int) -> int:
    result = 1
    for value in values:
        if type(value) is not int or value < 0 or result > MAX_I63 // max(value, 1):
            die("campaign cardinality exceeds signed 63-bit accounting")
        result *= value
    return result


def loss_seed(seed: str, K: int) -> int:
    return (
        int(seed, 16) ^ (K * 0x9E3779B97F4A7C15) ^
        (64 * 0xBF58476D1CE4E5B9)
    ) & MASK64


def construction_seed_values(root: str) -> tuple[int, int]:
    """Return the exact g8iv matrix/peel seeds for one uniform root."""
    matrix = strict_hex(root, "construction root", 64)
    peel = ((matrix & 0xffffffff) ^ (matrix >> 32)) & 0xffffffff
    return matrix, peel


def construction_seed_text(index: int) -> str:
    """Canonical decimal spelling used by the Stage-A v2 JSON contract."""
    if type(index) is not int or index not in range(len(CONSTRUCTION_ROOTS)):
        die("construction root index is outside the sealed domain")
    return str(int(CONSTRUCTION_ROOTS[index], 16))


def band_index(K: int) -> int:
    for index, (low, high) in enumerate(BANDS):
        if low <= K <= high:
            return index
    die(f"K={K} is outside the sealed band partition")
    raise AssertionError


CellKey = Tuple[int, int, str, int]
StratumKey = Tuple[int, int, str, int]


def stratum_key(
    construction_seed_index: int, packet_trace_root_index: int,
    schedule: str, K: int,
) -> StratumKey:
    return (
        construction_seed_index, packet_trace_root_index, schedule,
        band_index(K),
    )


def cell_ordinal(key: CellKey) -> int:
    construction_seed_index, packet_trace_root_index, schedule, K = key
    return (((
        construction_seed_index * len(LOSS_ROOTS) + packet_trace_root_index
    ) * len(SCHEDULES) + SCHEDULES.index(schedule)) * K_COUNT + K - K_MIN)


def cell_key_payload(key: CellKey) -> dict[str, Any]:
    construction_seed_index, packet_trace_root_index, schedule, K = key
    return {
        "construction_seed_index": construction_seed_index,
        "construction_seed": construction_seed_text(construction_seed_index),
        "packet_trace_root_index": packet_trace_root_index,
        "packet_trace_root": LOSS_ROOTS[packet_trace_root_index],
        "schedule": schedule, "K": K,
    }


SOURCE_IDENTITY_FIELDS = (
    "construction_seed_policy", "production_seed_fixups_applied",
    "base_matrix_seed", "base_peel_seed", "matrix_seed", "peel_seed",
    "actual_staircase_rows", "actual_dense_rows", "actual_source_hits",
    "actual_dense_two_anchor", "actual_dense_two_anchor_phase",
    "packet_trace_seed", "packet_trace_sha256",
)


def source_cell_record(
    key: CellKey, arms: Mapping[str, Mapping[str, str]],
) -> dict[str, Any]:
    if set(arms) != set(ARMS):
        die(f"source cell {key} is not an exact arm pair")
    h12, h13 = arms["h12"], arms["h13"]
    for field in SOURCE_IDENTITY_FIELDS:
        if field not in h12 or field not in h13 or h12[field] != h13[field]:
            die(f"source pair identity {field} differs at {key}")
    failures: dict[str, bool] = {}
    for arm, row in arms.items():
        if row.get("success") not in ("0", "1") or \
                row.get("rank_fail") not in ("0", "1") or \
                row.get("error") not in ("0", "1"):
            die(f"source outcome is malformed at {key}/{arm}")
        if sum(int(row[name]) for name in ("success", "rank_fail", "error")) != 1:
            die(f"source outcome is not one-hot at {key}/{arm}")
        failures[arm] = row["rank_fail"] == "1" or row["error"] == "1"
    record: dict[str, Any] = {
        **cell_key_payload(key),
        "raw_failures": failures,
        "raw_systematic_probe_results": {},
    }
    for arm, row in arms.items():
        probe = strict_uint(
            row.get("systematic_probe_result", ""),
            f"source systematic probe {key}/{arm}", 1,
        )
        record["raw_systematic_probe_results"][arm] = probe
    integer_fields = {
        "production_seed_fixups_applied",
        "actual_staircase_rows", "actual_dense_rows", "actual_source_hits",
        "actual_dense_two_anchor", "actual_dense_two_anchor_phase",
    }
    for field in SOURCE_IDENTITY_FIELDS:
        value = h12[field]
        if type(value) is not str:
            die(f"source identity {field} has the wrong type at {key}")
        record[field] = (
            strict_uint(value, f"source {field}")
            if field in integer_fields else value
        )
    for name, bits in (
            ("base_matrix_seed", 64), ("matrix_seed", 64),
            ("base_peel_seed", 32), ("peel_seed", 32),
            ("packet_trace_seed", 64)):
        strict_hex(str(record[name]), f"source {name}", bits)
    if (type(record["packet_trace_sha256"]) is not str or
            not re.fullmatch(r"[0-9a-f]{64}", record["packet_trace_sha256"])):
        die(f"source packet trace digest is malformed at {key}")
    construction_seed_index, packet_trace_root_index, unused_schedule, K = key
    del unused_schedule
    matrix_seed, peel_seed = construction_seed_values(
        CONSTRUCTION_ROOTS[construction_seed_index],
    )
    if record["packet_trace_seed"] != hex(
            loss_seed(LOSS_ROOTS[packet_trace_root_index], K)):
        die(f"source packet trace seed differs from the frozen formula at {key}")
    if (record["construction_seed_policy"] !=
            "matrix-c-peel-lo32-xor-hi32-v1" or
            record["production_seed_fixups_applied"] != 0 or
            any(row.get("construction_seed") != str(matrix_seed)
                for row in arms.values()) or
            record["base_matrix_seed"] != hex(matrix_seed) or
            record["matrix_seed"] != hex(matrix_seed) or
            record["base_peel_seed"] != hex(peel_seed) or
            record["peel_seed"] != hex(peel_seed)):
        die(f"source construction seeds are not exact uniform-root v1 at {key}")
    if (record["actual_dense_rows"] != 12 or
            record["actual_source_hits"] != (2 if K < 10000 else 3) or
            record["actual_dense_two_anchor"] != 1 or
            record["actual_dense_two_anchor_phase"] != 0):
        die(f"source construction geometry differs at {key}")
    return record


def _control_digest(
    construction_seed_index: int, packet_trace_root_index: int,
    schedule: str, band: int, K: int,
) -> str:
    low, high = BANDS[band]
    material = (
        f"{SCHEMA}|matched-control|{construction_seed_index}|"
        f"{construction_seed_text(construction_seed_index)}|"
        f"{packet_trace_root_index}|{LOSS_ROOTS[packet_trace_root_index]}|"
        f"{schedule}|{low}|{high}|{K}"
    )
    return sha256_bytes(material.encode("ascii"))


class CohortDeriver:
    """Streaming exact-census cohort selection with bounded retained state."""

    def __init__(
        self, sealed_union: Sequence[Mapping[str, Any]], *,
        expected_cells: int = SOURCE_PAIRED_CELLS,
        expected_union: int = EXPECTED_SOURCE_UNION,
        require_formal_strata: bool = True,
    ) -> None:
        if type(expected_cells) is not int or expected_cells < 1:
            die("expected source cell count is malformed")
        if type(expected_union) is not int or expected_union < 1:
            die("expected source union count is malformed")
        if not isinstance(sealed_union, Sequence) or len(sealed_union) != expected_union:
            die("source OH0 union does not have the required exact cardinality")
        self.expected_cells = expected_cells
        self.expected_union = expected_union
        self.require_formal_strata = require_formal_strata
        self.union_keys: list[CellKey] = []
        self.union_set: set[CellKey] = set()
        self.targets: dict[StratumKey, int] = {}
        previous = -1
        for item in sealed_union:
            if not isinstance(item, Mapping) or set(item) != {
                    "construction_seed_index", "construction_seed",
                    "packet_trace_root_index", "packet_trace_root",
                    "schedule", "K"}:
                die("source OH0 union record schema mismatch")
            construction_seed_index, packet_trace_root_index, schedule, K = (
                item["construction_seed_index"],
                item["packet_trace_root_index"], item["schedule"], item["K"],
            )
            if (type(construction_seed_index) is not int or
                    type(packet_trace_root_index) is not int or
                    type(schedule) is not str or type(K) is not int or
                    construction_seed_index not in range(len(CONSTRUCTION_ROOTS)) or
                    packet_trace_root_index not in range(len(LOSS_ROOTS)) or
                    item["construction_seed"] !=
                        construction_seed_text(construction_seed_index) or
                    item["packet_trace_root"] !=
                        LOSS_ROOTS[packet_trace_root_index] or
                    schedule not in SCHEDULES or not K_MIN <= K <= K_MAX):
                die("source OH0 union record is malformed")
            key = (
                construction_seed_index, packet_trace_root_index, schedule, K,
            )
            ordinal = cell_ordinal(key)
            if ordinal <= previous:
                die("source OH0 union is duplicated or not canonically ordered")
            previous = ordinal
            self.union_keys.append(key)
            self.union_set.add(key)
            stratum = stratum_key(*key)
            self.targets[stratum] = self.targets.get(stratum, 0) + 1
        if require_formal_strata:
            nonempty = set(self.targets)
            zero = {
                (construction_seed_index, packet_trace_root_index, schedule, band)
                for construction_seed_index in range(len(CONSTRUCTION_ROOTS))
                for packet_trace_root_index in range(len(LOSS_ROOTS))
                for schedule in SCHEDULES
                for band in range(len(BANDS))
            } - nonempty
            if (len(nonempty) != EXPECTED_NONEMPTY_STRATA or
                    tuple(sorted(zero, key=_stratum_order)) !=
                    tuple(sorted(EXPECTED_ZERO_STRATA, key=_stratum_order))):
                die("source union stratum/zero-stratum invariant mismatch")
        self.identity = hashlib.sha256()
        self.count = 0
        self.previous_ordinal = -1
        self.seen_union: dict[CellKey, dict[str, Any]] = {}
        # Heap items use a negative digest integer so the largest retained
        # selection digest is at index zero and can be replaced in O(log m).
        self.controls: dict[
            StratumKey, list[tuple[int, int, str, dict[str, Any]]]
        ] = {key: [] for key in self.targets}

    def observe(
        self, key: CellKey, arms: Mapping[str, Mapping[str, str]],
    ) -> None:
        construction_seed_index, packet_trace_root_index, schedule, K = key
        if (construction_seed_index not in range(len(CONSTRUCTION_ROOTS)) or
                packet_trace_root_index not in range(len(LOSS_ROOTS)) or
                schedule not in SCHEDULES or
                not K_MIN <= K <= K_MAX):
            die("source stream contains an out-of-domain cell")
        ordinal = cell_ordinal(key)
        if ordinal <= self.previous_ordinal:
            die("source stream is duplicated or not canonically ordered")
        self.previous_ordinal = ordinal
        record = source_cell_record(key, arms)
        self.identity.update((canonical_json(record) + "\n").encode())
        self.count += 1
        if key in self.union_set:
            if key in self.seen_union or not any(record["raw_failures"].values()):
                die("source sealed union disagrees with replayed outcomes")
            self.seen_union[key] = record
            return
        stratum = stratum_key(*key)
        target = self.targets.get(stratum, 0)
        if target == 0 or any(record["raw_failures"].values()):
            return
        digest = _control_digest(
            construction_seed_index, packet_trace_root_index, schedule,
            stratum[3], K,
        )
        item = (-int(digest, 16), -K, digest, record)
        heap = self.controls[stratum]
        if len(heap) < target:
            heapq.heappush(heap, item)
        elif item > heap[0]:
            heapq.heapreplace(heap, item)

    def finish(self) -> dict[str, Any]:
        if self.count != self.expected_cells:
            die("source stream does not contain the exact OH0 census")
        if set(self.seen_union) != self.union_set:
            die("source stream omitted or changed a sealed union cell")
        matches: list[dict[str, Any]] = []
        counts: list[dict[str, Any]] = []
        match_id = 0
        for stratum in sorted(self.targets, key=_stratum_order):
            construction_seed_index, packet_trace_root_index, schedule, band = stratum
            target = self.targets[stratum]
            heap = self.controls[stratum]
            if len(heap) != target:
                die(f"insufficient matched controls in stratum {stratum}")
            unions = sorted(
                (self.seen_union[key] for key in self.union_keys
                 if stratum_key(*key) == stratum),
                key=lambda record: record["K"],
            )
            controls = sorted(
                ((item[2], item[3]) for item in heap),
                key=lambda value: (value[0], value[1]["K"]),
            )
            if len(unions) != target or len(controls) != target:
                die("matched-control stratum cardinality mismatch")
            low, high = BANDS[band]
            counts.append({
                "construction_seed_index": construction_seed_index,
                "construction_seed": construction_seed_text(
                    construction_seed_index),
                "packet_trace_root_index": packet_trace_root_index,
                "packet_trace_root": LOSS_ROOTS[packet_trace_root_index],
                "schedule": schedule, "band_index": band,
                "band": [low, high], "union": target, "controls": target,
            })
            for union, (selection_sha256, control) in zip(unions, controls):
                matches.append({
                    "match_id": match_id,
                    "stratum": {
                        "construction_seed_index": construction_seed_index,
                        "construction_seed": construction_seed_text(
                            construction_seed_index),
                        "packet_trace_root_index": packet_trace_root_index,
                        "packet_trace_root": LOSS_ROOTS[packet_trace_root_index],
                        "schedule": schedule, "band_index": band,
                        "band": [low, high],
                    },
                    "control_selection_sha256": selection_sha256,
                    "union": union, "control": control,
                })
                match_id += 1
        if len(matches) != self.expected_union:
            die("derived matched cohort cardinality mismatch")
        union_keys = {
            (match["union"]["construction_seed_index"],
             match["union"]["packet_trace_root_index"],
             match["union"]["schedule"], match["union"]["K"])
            for match in matches
        }
        control_keys = {
            (match["control"]["construction_seed_index"],
             match["control"]["packet_trace_root_index"],
             match["control"]["schedule"], match["control"]["K"])
            for match in matches
        }
        if (len(union_keys) != self.expected_union or
                len(control_keys) != self.expected_union or
                union_keys & control_keys):
            die("derived controls are duplicated or overlap the union")
        return {
            "selection_schema": f"{SCHEMA}.matched_control.v1",
            "selection_rule": SELECTION_RULE,
            "source_union_count": self.expected_union,
            "matched_control_count": self.expected_union,
            "paired_cells": checked_product(self.expected_union, 2),
            "arm_outcomes": checked_product(self.expected_union, 4),
            "bands": [list(value) for value in BANDS],
            "nonempty_strata": len(self.targets),
            "zero_strata": [
                {"construction_seed_index": construction_seed_index,
                 "construction_seed": construction_seed_text(
                    construction_seed_index),
                 "packet_trace_root_index": packet_trace_root_index,
                 "packet_trace_root": LOSS_ROOTS[packet_trace_root_index],
                 "schedule": schedule, "band_index": band,
                 "band": list(BANDS[band])}
                for construction_seed_index in range(len(CONSTRUCTION_ROOTS))
                for packet_trace_root_index in range(len(LOSS_ROOTS))
                for schedule in SCHEDULES
                for band in range(len(BANDS))
                if (construction_seed_index, packet_trace_root_index,
                    schedule, band) not in self.targets
            ],
            "source_identity_stream_sha256": self.identity.hexdigest(),
            "counts_by_stratum_band": counts,
            "matches": matches,
        }


def _stratum_order(value: StratumKey) -> tuple[int, int, int, int]:
    return value[0], value[1], SCHEDULES.index(value[2]), value[3]


@dataclass(frozen=True)
class Task:
    phase: str
    overhead: int
    job: int
    construction_seed_index: int
    construction_seed: str
    packet_trace_root_index: int
    packet_trace_root: str
    schedule: str
    match_ids: tuple[int, ...]
    ks: tuple[int, ...]
    task_id: str

    @property
    def stem(self) -> str:
        return f"job{self.job:07d}.{self.task_id}"

    def core(self) -> dict[str, Any]:
        return {
            "phase": self.phase, "overhead": self.overhead, "job": self.job,
            "construction_seed_index": self.construction_seed_index,
            "construction_seed": self.construction_seed,
            "packet_trace_root_index": self.packet_trace_root_index,
            "packet_trace_root": self.packet_trace_root,
            "schedule": self.schedule, "match_ids": list(self.match_ids),
            "ks": list(self.ks),
        }

    def payload(self) -> dict[str, Any]:
        return {**self.core(), "task_id": self.task_id}


def task_digest(core: Mapping[str, Any]) -> str:
    return sha256_bytes(canonical_json(core).encode())[:24]


def _new_task(
    phase: str, job: int, construction_seed_index: int,
    packet_trace_root_index: int, schedule: str,
    match_ids: Sequence[int], ks: Sequence[int],
) -> Task:
    provisional = Task(
        phase, 0, job,
        construction_seed_index, CONSTRUCTION_ROOTS[construction_seed_index],
        packet_trace_root_index, LOSS_ROOTS[packet_trace_root_index], schedule,
        tuple(match_ids), tuple(ks), "",
    )
    return Task(**{
        **provisional.__dict__, "task_id": task_digest(provisional.core()),
    })


def validate_task(task: Task) -> None:
    if (any(type(value) is not int for value in (
                task.overhead, task.job, task.construction_seed_index,
                task.packet_trace_root_index,
                *task.match_ids, *task.ks)) or
            any(type(value) is not str for value in (
                task.phase, task.construction_seed, task.packet_trace_root,
                task.schedule, task.task_id)) or
            task.phase not in ("screen", "allk") or task.overhead != 0 or
            task.job < 0 or
            task.construction_seed_index not in range(len(CONSTRUCTION_ROOTS)) or
            task.construction_seed !=
                CONSTRUCTION_ROOTS[task.construction_seed_index] or
            task.packet_trace_root_index not in range(len(LOSS_ROOTS)) or
            task.packet_trace_root != LOSS_ROOTS[task.packet_trace_root_index] or
            task.schedule not in SCHEDULES or not task.ks or
            tuple(sorted(task.ks)) != task.ks or len(set(task.ks)) != len(task.ks) or
            len(task.ks) > CHUNK_SIZE or
            any(not K_MIN <= K <= K_MAX for K in task.ks) or
            any(value < 0 for value in task.match_ids) or
            tuple(sorted(task.match_ids)) != task.match_ids or
            len(set(task.match_ids)) != len(task.match_ids) or
            (task.phase == "screen" and
             (not task.match_ids or len(task.match_ids) > SCREEN_MATCH_CHUNK or
              len(task.ks) != 2 * len(task.match_ids))) or
            (task.phase == "allk" and task.match_ids) or
            task.task_id != task_digest(task.core())):
        die("task violates the sealed Stage-B geometry")


def task_from_payload(payload: Any) -> Task:
    fields = {
        "phase", "overhead", "job", "construction_seed_index",
        "construction_seed", "packet_trace_root_index", "packet_trace_root",
        "schedule", "match_ids", "ks", "task_id",
    }
    if not isinstance(payload, dict) or set(payload) != fields:
        die("task has an unexpected schema")
    if (any(type(payload[name]) is not int for name in (
            "overhead", "job", "construction_seed_index",
            "packet_trace_root_index")) or
            any(type(payload[name]) is not str for name in (
                "phase", "construction_seed", "packet_trace_root",
                "schedule", "task_id")) or
            type(payload["match_ids"]) is not list or
            type(payload["ks"]) is not list or
            any(type(value) is not int for value in
                [*payload["match_ids"], *payload["ks"]])):
        die("task scalar types are not canonical")
    task = Task(
        payload["phase"], payload["overhead"], payload["job"],
        payload["construction_seed_index"], payload["construction_seed"],
        payload["packet_trace_root_index"], payload["packet_trace_root"],
        payload["schedule"],
        tuple(payload["match_ids"]), tuple(payload["ks"]), payload["task_id"],
    )
    validate_task(task)
    return task


def build_screen_tasks(cohort: Mapping[str, Any]) -> list[Task]:
    matches = cohort.get("matches") if isinstance(cohort, Mapping) else None
    if not isinstance(matches, list) or not matches:
        die("matched cohort has no matches")
    grouped: dict[StratumKey, list[Mapping[str, Any]]] = {}
    for expected_id, match in enumerate(matches):
        if (not isinstance(match, Mapping) or
                type(match.get("match_id")) is not int or
                match.get("match_id") != expected_id or
                not isinstance(match.get("stratum"), Mapping)):
            die("matched cohort has malformed or noncontiguous match IDs")
        stratum = match["stratum"]
        key = (
            stratum.get("construction_seed_index"),
            stratum.get("packet_trace_root_index"),
            stratum.get("schedule"), stratum.get("band_index"),
        )
        if (type(key[0]) is not int or type(key[1]) is not int or
                type(key[2]) is not str or type(key[3]) is not int or
                key[0] not in range(len(CONSTRUCTION_ROOTS)) or
                key[1] not in range(len(LOSS_ROOTS)) or
                key[2] not in SCHEDULES or
                key[3] not in range(len(BANDS)) or
                stratum.get("construction_seed") !=
                    construction_seed_text(key[0]) or
                stratum.get("packet_trace_root") != LOSS_ROOTS[key[1]] or
                stratum.get("band") != list(BANDS[key[3]])):
            die("matched cohort stratum is malformed")
        grouped.setdefault(key, []).append(match)
    tasks: list[Task] = []
    for stratum in sorted(grouped, key=_stratum_order):
        construction_seed_index, packet_trace_root_index, schedule, band = stratum
        if band not in range(len(BANDS)):
            die("matched cohort band is outside the sealed partition")
        values = grouped[stratum]
        for offset in range(0, len(values), SCREEN_MATCH_CHUNK):
            chunk = values[offset:offset + SCREEN_MATCH_CHUNK]
            match_ids = [int(match["match_id"]) for match in chunk]
            keys: list[int] = []
            for match in chunk:
                for role in ("union", "control"):
                    cell = match.get(role)
                    if (not isinstance(cell, Mapping) or
                            cell.get("construction_seed_index") !=
                                construction_seed_index or
                            cell.get("packet_trace_root_index") !=
                                packet_trace_root_index or
                            cell.get("construction_seed") !=
                                construction_seed_text(construction_seed_index) or
                            cell.get("packet_trace_root") !=
                                LOSS_ROOTS[packet_trace_root_index] or
                            cell.get("schedule") != schedule or
                            type(cell.get("K")) is not int or
                            band_index(cell["K"]) != band):
                        die("matched cohort cell differs from its stratum")
                    keys.append(cell["K"])
            if len(set(keys)) != len(keys):
                die("one screen task would contain a duplicate K")
            tasks.append(_new_task(
                "screen", len(tasks), construction_seed_index,
                packet_trace_root_index, schedule,
                match_ids, sorted(keys),
            ))
    for task in tasks:
        validate_task(task)
    return tasks


def build_allk_tasks() -> list[Task]:
    tasks: list[Task] = []
    for construction_seed_index in range(len(CONSTRUCTION_ROOTS)):
        for packet_trace_root_index in range(len(LOSS_ROOTS)):
            for schedule in SCHEDULES:
                for low in range(K_MIN, K_MAX + 1, CHUNK_SIZE):
                    high = min(low + CHUNK_SIZE - 1, K_MAX)
                    tasks.append(_new_task(
                        "allk", len(tasks), construction_seed_index,
                        packet_trace_root_index, schedule, (),
                        range(low, high + 1),
                    ))
    if (len(tasks) != ALLK_TASKS or
            sum(len(task.ks) for task in tasks) != ALLK_CELLS):
        die("all-K planner cardinality differs from the sealed census")
    return tasks


def manifest_payload(phase: str, tasks: Sequence[Task]) -> dict[str, Any]:
    if (phase not in ("screen", "allk") or
            any(task.phase != phase for task in tasks)):
        die("manifest payload phase differs from its tasks")
    for task in tasks:
        validate_task(task)
    return {
        "phase": phase, "overhead": 0,
        "pair_native": True, "arms_per_N": 2,
        "task_count": len(tasks),
        "paired_cell_count": sum(len(task.ks) for task in tasks),
        "arm_outcome_count": checked_product(
            sum(len(task.ks) for task in tasks), len(ARMS)),
        "tasks": [task.payload() for task in tasks],
    }


def load_manifest(stage: Path, phase: str, cohort: Mapping[str, Any]) -> list[Task]:
    if phase not in ("screen", "allk"):
        die("unknown Stage-B manifest phase")
    payload = load_sealed(stage / "manifest.json", f"{SCHEMA}.stage_manifest")
    if (not isinstance(payload, dict) or set(payload) != {
            "phase", "overhead", "pair_native", "arms_per_N", "task_count",
            "paired_cell_count", "arm_outcome_count", "tasks"} or
            payload["phase"] != phase or payload["overhead"] != 0 or
            payload["pair_native"] is not True or payload["arms_per_N"] != 2 or
            any(type(payload[name]) is not int for name in (
                "overhead", "arms_per_N", "task_count",
                "paired_cell_count", "arm_outcome_count")) or
            not isinstance(payload["tasks"], list)):
        die("Stage-B manifest schema mismatch")
    tasks = [task_from_payload(value) for value in payload["tasks"]]
    expected = build_screen_tasks(cohort) if phase == "screen" else build_allk_tasks()
    if (canonical_json([task.payload() for task in tasks]) !=
            canonical_json([task.payload() for task in expected]) or
            payload["task_count"] != len(tasks) or
            payload["paired_cell_count"] != sum(len(task.ks) for task in tasks) or
            payload["arm_outcome_count"] != 2 * payload["paired_cell_count"]):
        die("Stage-B manifest differs from its deterministic plan")
    return tasks


def _load_stage_a_module(controller: Path, digest: str) -> Any:
    source = stable_bytes(controller)
    if sha256_bytes(source) != digest:
        die("Stage-A frozen controller changed before import")
    module_name = f"_wirehair_stage_a_{digest}"
    existing = sys.modules.get(module_name)
    if existing is not None:
        if Path(existing.__file__).resolve() != controller.resolve():
            die("Stage-A module-name collision")
        return existing
    module = types.ModuleType(module_name)
    module.__file__ = str(controller)
    module.__package__ = ""
    module.__loader__ = None
    module.__spec__ = None
    sys.modules[module_name] = module
    try:
        code = compile(source, str(controller), "exec", dont_inherit=True)
        exec(code, module.__dict__)
    except BaseException:
        sys.modules.pop(module_name, None)
        raise
    if (getattr(module, "SCHEMA", None) != SOURCE_SCHEMA or
            getattr(module, "K_MIN", None) != K_MIN or
            getattr(module, "K_MAX", None) != K_MAX or
            tuple(getattr(module, "LOSS_ROOTS", ())) != LOSS_ROOTS or
            tuple(hex(value) for value in
                  getattr(module, "CONSTRUCTION_SEEDS", ())) !=
                CONSTRUCTION_ROOTS or
            getattr(module, "CONSTRUCTION_SEED_POLICY", None) !=
                "matrix-c-peel-lo32-xor-hi32-v1" or
            tuple(getattr(module, "SCHEDULES", ())) != SCHEDULES):
        die("authenticated Stage-A controller constants differ from Stage-B")
    return module


def _source_affinity(module: Any, result_dir: Path) -> list[int]:
    payload = module.load_sealed(
        result_dir / "controller_affinity.json",
        f"{SOURCE_SCHEMA}.controller_affinity",
    )
    if (not isinstance(payload, dict) or set(payload) != {
            "selected_cpus", "actual_cpus", "reserved_cpu_127_excluded"} or
            type(payload["selected_cpus"]) is not list or
            any(type(cpu) is not int for cpu in payload["selected_cpus"]) or
            payload["actual_cpus"] != payload["selected_cpus"] or
            payload["reserved_cpu_127_excluded"] is not True or
            not payload["selected_cpus"] or 127 in payload["selected_cpus"] or
            payload["selected_cpus"] != sorted(set(payload["selected_cpus"]))):
        die("Stage-A controller affinity receipt is malformed")
    return payload["selected_cpus"]


def _source_telemetry_receipts(
    module: Any, result_dir: Path, contract: Mapping[str, Any],
) -> tuple[Optional[Mapping[str, Any]], Optional[Mapping[str, Any]]]:
    configured = contract["external_telemetry"]["path"] is not None
    if not configured:
        if ((result_dir / "telemetry_start.json").exists() or
                (result_dir / "telemetry_end.json").exists()):
            die("Stage-A null-telemetry campaign has interval receipts")
        return None, None
    start = module.load_sealed(
        result_dir / "telemetry_start.json", f"{SOURCE_SCHEMA}.telemetry_start",
    )
    end = module.load_sealed(
        result_dir / "telemetry_end.json", f"{SOURCE_SCHEMA}.telemetry_end",
    )
    if not isinstance(start, dict) or not isinstance(end, dict):
        die("Stage-A telemetry receipt is not an object")
    return start, end


def source_shard_epoch_sha256(module: Any, stage: Path, task: Any) -> str:
    """Digest the exact three-file source shard bound to one Stage-A task."""
    root = module.shard_path(stage, task)
    if not root.is_dir() or root.is_symlink():
        die("Stage-A source shard epoch root is missing or unsafe")
    inventory = {entry.name: entry for entry in root.iterdir()}
    names = ("stdout.csv", "stderr.txt", "receipt.json")
    if (set(inventory) != set(names) or
            any(inventory[name].is_symlink() or
                not inventory[name].is_file() for name in names)):
        die("Stage-A source shard epoch inventory differs from its ledger")
    payload = {
        "task": task.payload(),
        "files": [
            {"name": name, "sha256": sha256_file(inventory[name])}
            for name in names
        ],
    }
    return sha256_bytes(canonical_json(payload).encode())



def authenticate_stage_a(result_dir: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    """Authenticate and independently replay all three uniform-root OH0 censuses."""
    result_dir = result_dir.resolve(strict=True)
    if not result_dir.is_dir() or result_dir.is_symlink():
        die("Stage-A result is not a real directory")
    independent = load_sealed(
        result_dir / "contract.json", f"{SOURCE_SCHEMA}.contract",
    )
    if (not isinstance(independent, dict) or
            independent.get("result_dir") != str(result_dir)):
        die("Stage-A contract does not bind its canonical result directory")
    for name in (
            "binary", "binary_sha256", "controller", "controller_sha256",
            "taskset", "taskset_sha256"):
        if type(independent.get(name)) is not str:
            die(f"Stage-A contract lacks canonical {name}")
    if (independent["controller_sha256"] !=
            FORMAL_STAGE_A_CONTROLLER_SHA256 or
            independent["binary_sha256"] != FORMAL_STAGE_A_BINARY_SHA256):
        die("Stage-A controller/binary is not the pinned formal build")
    for name in ("binary", "controller", "taskset"):
        path = Path(independent[name])
        digest = independent[f"{name}_sha256"]
        if (path.parent != result_dir / "frozen" or path.is_symlink() or
                not path.is_file() or not re.fullmatch(r"[0-9a-f]{64}", digest) or
                sha256_file(path) != digest):
            die(f"Stage-A frozen {name} path/hash mismatch")
    require_exact_seal(
        result_dir / "campaign_complete.json",
        f"{SOURCE_SCHEMA}.campaign_complete",
        FORMAL_STAGE_A_CAMPAIGN_COMPLETE_SEAL,
        FORMAL_STAGE_A_CAMPAIGN_COMPLETE_FILE_SHA256,
    )
    require_exact_seal(
        result_dir / "analysis.json", f"{SOURCE_SCHEMA}.analysis_record",
        FORMAL_STAGE_A_ANALYSIS_SEAL, FORMAL_STAGE_A_ANALYSIS_FILE_SHA256,
    )
    module = _load_stage_a_module(
        Path(independent["controller"]), independent["controller_sha256"],
    )
    try:
        contract = module.load_contract(
            result_dir, require_frozen_controller=False,
        )
        if canonical_json(contract) != canonical_json(independent):
            die("authenticated Stage-A controller parsed a different contract")
        allowed_cpus = set(_source_affinity(module, result_dir))
        common_paths = {
            "contract_sha256": result_dir / "contract.json",
            "controller_affinity_sha256":
                result_dir / "controller_affinity.json",
            "terminal_sha256": result_dir / "terminal.json",
        }
        root_stages = [
            module.stage_path(result_dir, index, 0)
            for index in range(len(CONSTRUCTION_ROOTS))
        ]
        root_paths = [{
            "construction_seed_index": index,
            "manifest": stage / "manifest.json",
            "complete": stage / "complete.json",
        } for index, stage in enumerate(root_stages)]
        initial_common = {
            name: sha256_file(path) for name, path in common_paths.items()
        }
        initial_root_hashes = [{
            "construction_seed_index": item["construction_seed_index"],
            "manifest_sha256": sha256_file(item["manifest"]),
            "complete_sha256": sha256_file(item["complete"]),
        } for item in root_paths]
        terminal = module.verify_terminal_campaign(
            result_dir, contract, allowed_cpus,
        )
        if (not isinstance(terminal, dict) or
                terminal.get("disposition") != "confirmation_complete" or
                terminal.get("architecture_accepted") is not True or
                terminal.get("executed_construction_seed_indices") != [0, 1, 2] or
                terminal.get("production_seed_fixups_applied") != 0):
            die("Stage-A did not complete the exact accepted R0+R1 campaign")
        start, end = _source_telemetry_receipts(module, result_dir, contract)
        completion = module.verify_campaign_completion(
            result_dir, contract, start, end,
        )
        if (not isinstance(completion, dict) or
                completion.get("contract_sha256") !=
                    initial_common["contract_sha256"] or
                completion.get("controller_affinity_sha256") !=
                    initial_common["controller_affinity_sha256"] or
                completion.get("terminal_sha256") !=
                    initial_common["terminal_sha256"]):
            die("Stage-A campaign completion/control-leaf mismatch")
        analysis = module.load_sealed(
            result_dir / "analysis.json", f"{SOURCE_SCHEMA}.analysis_record",
        )
        analysis_fields = {
            "schema", "coverage", "terminal", "rounds",
            "by_construction_seed", "interpretation", "loss_seed_formula",
            "external_telemetry", "campaign_completion",
        }
        if (not isinstance(analysis, dict) or set(analysis) != analysis_fields or
                analysis.get("schema") != f"{SOURCE_SCHEMA}.analysis" or
                type(analysis.get("interpretation")) is not str or
                analysis.get("loss_seed_formula") != LOSS_SEED_FORMULA or
                not isinstance(analysis.get("external_telemetry"), dict) or
                canonical_json(analysis.get("terminal")) !=
                    canonical_json(terminal) or
                canonical_json(analysis.get("campaign_completion")) !=
                    canonical_json(completion)):
            die("Stage-A analysis terminal/completion cross-seal mismatch")
        coverage = analysis.get("coverage")
        expected_construction_seeds = [
            str(int(root, 16)) for root in CONSTRUCTION_ROOTS
        ]
        coverage_fields = {
            "K_min", "K_max", "unique_K", "construction_seed_policy",
            "construction_seeds_predeclared", "construction_seeds_executed",
            "packet_trace_roots", "schedules",
            "paired_cells_per_construction_root", "r0_paired_cells",
            "r1_paired_cells", "full_paired_cells",
            "production_seed_fixups_applied",
        }
        if (not isinstance(coverage, dict) or set(coverage) != coverage_fields or
                any(type(coverage.get(name)) is not int for name in (
                    "K_min", "K_max", "unique_K", "packet_trace_roots",
                    "schedules", "paired_cells_per_construction_root",
                    "r0_paired_cells", "r1_paired_cells",
                    "full_paired_cells", "production_seed_fixups_applied",
                )) or
                coverage.get("K_min") != K_MIN or
                coverage.get("K_max") != K_MAX or
                coverage.get("unique_K") != K_COUNT or
                coverage.get("construction_seed_policy") !=
                    "matrix-c-peel-lo32-xor-hi32-v1" or
                coverage.get("construction_seeds_predeclared") !=
                    expected_construction_seeds or
                coverage.get("construction_seeds_executed") != [0, 1, 2] or
                coverage.get("packet_trace_roots") != len(LOSS_ROOTS) or
                coverage.get("schedules") != len(SCHEDULES) or
                coverage.get("paired_cells_per_construction_root") !=
                    ALLK_CELLS_PER_CONSTRUCTION_ROOT or
                coverage.get("r0_paired_cells") !=
                    ALLK_CELLS_PER_CONSTRUCTION_ROOT or
                coverage.get("r1_paired_cells") !=
                    2 * ALLK_CELLS_PER_CONSTRUCTION_ROOT or
                coverage.get("full_paired_cells") != SOURCE_PAIRED_CELLS or
                coverage.get("production_seed_fixups_applied") != 0):
            die("Stage-A analysis coverage is not the exact uniform-root census")
        rounds = analysis.get("rounds")
        by_construction = analysis.get("by_construction_seed")
        expected_rounds = {
            "R0": ([0], ALLK_CELLS_PER_CONSTRUCTION_ROOT),
            "R1": ([1, 2], 2 * ALLK_CELLS_PER_CONSTRUCTION_ROOT),
            "FULL": ([0, 1, 2], SOURCE_PAIRED_CELLS),
        }
        round_fields = {
            "construction_seed_indices", "paired_cells", "oh0",
            "minimum_success_overhead", "right_censoring",
        }
        root_fields = {
            "construction_seed", "campaign_round", "terminal_overhead",
            "unresolved_union_failures", "oh0", "minimum_success_overhead",
            "right_censoring",
        }
        if (not isinstance(rounds, dict) or set(rounds) != set(expected_rounds) or
                any(not isinstance(rounds[name], dict) or
                    set(rounds[name]) != round_fields or
                    rounds[name].get("construction_seed_indices") != indices or
                    type(rounds[name].get("paired_cells")) is not int or
                    rounds[name].get("paired_cells") != paired_cells
                    for name, (indices, paired_cells) in
                    expected_rounds.items()) or
                not isinstance(by_construction, dict) or
                set(by_construction) != {"0", "1", "2"} or
                any(not isinstance(by_construction[str(index)], dict) or
                    set(by_construction[str(index)]) != root_fields or
                    by_construction[str(index)].get("construction_seed") !=
                        construction_seed_text(index) or
                    by_construction[str(index)].get("campaign_round") !=
                        ("R0" if index == 0 else "R1")
                    for index in range(len(CONSTRUCTION_ROOTS)))):
            die("Stage-A analysis round/root schema mismatch")

        tasks_by_root: list[list[Any]] = []
        sealed_union: list[dict[str, Any]] = []
        for construction_index, stage in enumerate(root_stages):
            tasks = module.load_manifest(stage, construction_index, 0)
            module.validate_stage_file_inventory(stage, complete_required=True)
            module.validate_stage_shard_inventory(stage, tasks)
            complete = module.load_sealed(
                stage / "complete.json", f"{SOURCE_SCHEMA}.stage_complete",
            )
            if (not isinstance(complete, dict) or set(complete) != {
                    "overhead", "job_count", "paired_cell_count", "metrics",
                    "union_failures"} or complete["overhead"] != 0 or
                    complete["job_count"] != len(tasks) or
                    len(tasks) != 2 * ALLK_TASKS_PER_CONSTRUCTION_ROOT or
                    complete["paired_cell_count"] !=
                        ALLK_CELLS_PER_CONSTRUCTION_ROOT or
                    canonical_json(complete["metrics"]) != canonical_json(
                        by_construction[str(construction_index)].get("oh0")) or
                    not isinstance(complete["union_failures"], list) or
                    any(item.get("construction_seed_index") != construction_index
                        for item in complete["union_failures"]
                        if isinstance(item, Mapping)) or
                    any(not isinstance(item, Mapping)
                        for item in complete["union_failures"])):
                die("Stage-A OH0 root completion/analysis mismatch")
            tasks_by_root.append(tasks)
            sealed_union.extend(complete["union_failures"])
        if len(sealed_union) != EXPECTED_SOURCE_UNION:
            die("Stage-A OH0 union does not match the pinned v2 cardinality")
        deriver = CohortDeriver(sealed_union)
        shard_epoch_sha256s: list[str] = []

        def observed_pairs() -> Iterator[
            tuple[CellKey, Mapping[str, Mapping[str, str]]]
        ]:
            for construction_index, (stage, tasks) in enumerate(
                    zip(root_stages, tasks_by_root)):
                for pair_index in range(len(tasks) // 2):
                    left, right = tasks[2 * pair_index:2 * pair_index + 2]
                    if (left.arm != "h12" or right.arm != "h13" or
                            left.pair != pair_index or right.pair != pair_index or
                            left.overhead != 0 or right.overhead != 0 or
                            left.construction_index != construction_index or
                            right.construction_index != construction_index or
                            left.construction_seed !=
                                int(CONSTRUCTION_ROOTS[construction_index], 16) or
                            right.construction_seed != left.construction_seed or
                            left.seed_index != right.seed_index or
                            left.seed != right.seed or
                            left.schedule != right.schedule or
                            tuple(left.ks) != tuple(right.ks)):
                        die("Stage-A OH0 manifest arm/root pairing mismatch")
                    before = (
                        source_shard_epoch_sha256(module, stage, left),
                        source_shard_epoch_sha256(module, stage, right),
                    )
                    left_rows = module.validate_shard(
                        stage, left, Path(contract["binary"]),
                        contract["binary_sha256"], Path(contract["taskset"]),
                        allowed_cpus,
                    )
                    right_rows = module.validate_shard(
                        stage, right, Path(contract["binary"]),
                        contract["binary_sha256"], Path(contract["taskset"]),
                        allowed_cpus,
                    )
                    after = (
                        source_shard_epoch_sha256(module, stage, left),
                        source_shard_epoch_sha256(module, stage, right),
                    )
                    if before != after:
                        die("Stage-A source shard changed during authenticated replay")
                    shard_epoch_sha256s.extend(before)
                    if (len(left_rows) != len(left.ks) or
                            len(right_rows) != len(right.ks)):
                        die("Stage-A source shard row count mismatch")
                    for K, h12, h13 in zip(left.ks, left_rows, right_rows):
                        key: CellKey = (
                            construction_index, left.seed_index,
                            left.schedule, K,
                        )
                        arms = {"h12": h12, "h13": h13}
                        module.validate_pair_receipt(key, arms)
                        deriver.observe(key, arms)
                        yield key, arms
                print(
                    f"# Stage-A OH0: verified construction root "
                    f"{construction_index + 1}/{len(CONSTRUCTION_ROOTS)}",
                    file=sys.stderr, flush=True,
                )

        replayed_metrics, replayed_union = module.aggregate_pair_stream(
            observed_pairs(),
        )
        if (canonical_json(replayed_metrics) !=
                canonical_json(rounds["FULL"].get("oh0")) or
                canonical_json(replayed_union) != canonical_json(sealed_union)):
            die("Stage-A combined OH0 completion differs from streamed replay")
        cohort = deriver.finish()
        if (cohort["source_union_count"] != EXPECTED_SOURCE_UNION or
                cohort["matched_control_count"] != EXPECTED_SCREEN_MATCHES or
                cohort["paired_cells"] != EXPECTED_SCREEN_CELLS or
                cohort["arm_outcomes"] != EXPECTED_SCREEN_OUTCOMES or
                cohort["nonempty_strata"] != EXPECTED_NONEMPTY_STRATA or
                len(cohort["zero_strata"]) != len(EXPECTED_ZERO_STRATA) or
                len(build_screen_tasks(cohort)) != EXPECTED_SCREEN_TASKS):
            die("Stage-A-derived screen cohort misses pinned v2 invariants")
        if len(shard_epoch_sha256s) != 2 * ALLK_TASKS:
            die("Stage-A OH0 shard epoch cardinality mismatch")
        final_campaign_complete = require_exact_seal(
            result_dir / "campaign_complete.json",
            f"{SOURCE_SCHEMA}.campaign_complete",
            FORMAL_STAGE_A_CAMPAIGN_COMPLETE_SEAL,
            FORMAL_STAGE_A_CAMPAIGN_COMPLETE_FILE_SHA256,
        )
        final_analysis = require_exact_seal(
            result_dir / "analysis.json", f"{SOURCE_SCHEMA}.analysis_record",
            FORMAL_STAGE_A_ANALYSIS_SEAL,
            FORMAL_STAGE_A_ANALYSIS_FILE_SHA256,
        )
        final_common = {
            name: sha256_file(path) for name, path in common_paths.items()
        }
        final_root_hashes = [{
            "construction_seed_index": item["construction_seed_index"],
            "manifest_sha256": sha256_file(item["manifest"]),
            "complete_sha256": sha256_file(item["complete"]),
        } for item in root_paths]
        if (canonical_json(final_campaign_complete) != canonical_json(completion) or
                canonical_json(final_analysis) != canonical_json(analysis) or
                final_common != initial_common or
                final_root_hashes != initial_root_hashes):
            die("formal Stage-A seals changed during authenticated replay")
    except CampaignError:
        raise
    except BaseException as exc:
        if isinstance(exc, (KeyboardInterrupt, SystemExit)):
            raise
        die(f"authenticated Stage-A validation failed: {exc}")
    provenance = {
        "result_dir": str(result_dir),
        "contract_sha256": initial_common["contract_sha256"],
        "controller_sha256": contract["controller_sha256"],
        "binary_sha256": contract["binary_sha256"],
        "taskset_sha256": contract["taskset_sha256"],
        "controller_affinity_sha256":
            initial_common["controller_affinity_sha256"],
        "terminal_sha256": initial_common["terminal_sha256"],
        "oh0_root_control_sha256s": initial_root_hashes,
        "campaign_complete_sha256":
            FORMAL_STAGE_A_CAMPAIGN_COMPLETE_FILE_SHA256,
        "analysis_sha256": FORMAL_STAGE_A_ANALYSIS_FILE_SHA256,
        "source_identity_stream_sha256":
            cohort["source_identity_stream_sha256"],
        "oh0_shard_epoch_sha256s": shard_epoch_sha256s,
        "oh0_shard_epoch_root_sha256": sha256_bytes(
            canonical_json(shard_epoch_sha256s).encode()),
        "paired_cells_replayed": SOURCE_PAIRED_CELLS,
        "union_cells_replayed": EXPECTED_SOURCE_UNION,
        "construction_seed_policy": "matrix-c-peel-lo32-xor-hi32-v1",
        "production_seed_fixups_applied": 0,
        "formal_completion_authenticated": True,
    }
    return cohort, provenance


class StageAAllKOracle:
    """Bounded per-task view of the three pinned Stage-A OH0 root censuses."""

    def __init__(
        self, source_result: Path, expected_provenance: Mapping[str, Any],
    ) -> None:
        self.result_dir = source_result.resolve(strict=True)
        if (not isinstance(expected_provenance, Mapping) or
                expected_provenance.get("result_dir") != str(self.result_dir) or
                expected_provenance.get("controller_sha256") !=
                    FORMAL_STAGE_A_CONTROLLER_SHA256 or
                expected_provenance.get("binary_sha256") !=
                    FORMAL_STAGE_A_BINARY_SHA256 or
                expected_provenance.get("construction_seed_policy") !=
                    "matrix-c-peel-lo32-xor-hi32-v1" or
                expected_provenance.get("production_seed_fixups_applied") != 0 or
                type(expected_provenance.get("taskset_sha256")) is not str or
                not re.fullmatch(
                    r"[0-9a-f]{64}", expected_provenance["taskset_sha256"])):
            die("all-K oracle expected provenance is malformed")
        control_fields = (
            "contract_sha256", "controller_affinity_sha256", "terminal_sha256",
            "campaign_complete_sha256", "analysis_sha256",
        )
        if any(type(expected_provenance.get(name)) is not str or
               not re.fullmatch(r"[0-9a-f]{64}", expected_provenance[name])
               for name in control_fields):
            die("all-K oracle expected control-leaf hashes are malformed")
        contract_path = self.result_dir / "contract.json"
        if sha256_file(contract_path) != expected_provenance["contract_sha256"]:
            die("all-K oracle source contract differs from authenticated replay")
        independent = load_sealed(
            contract_path, f"{SOURCE_SCHEMA}.contract",
        )
        if (not isinstance(independent, dict) or
                independent.get("controller_sha256") !=
                    FORMAL_STAGE_A_CONTROLLER_SHA256 or
                independent.get("binary_sha256") !=
                    FORMAL_STAGE_A_BINARY_SHA256 or
                independent.get("taskset_sha256") !=
                    expected_provenance.get("taskset_sha256") or
                sha256_file(contract_path) !=
                    expected_provenance["contract_sha256"]):
            die("all-K oracle source is not the pinned formal Stage-A build")
        controller = Path(independent.get("controller", ""))
        if (controller.parent != self.result_dir / "frozen" or
                sha256_file(controller) != FORMAL_STAGE_A_CONTROLLER_SHA256):
            die("all-K oracle Stage-A controller path/hash mismatch")
        require_exact_seal(
            self.result_dir / "campaign_complete.json",
            f"{SOURCE_SCHEMA}.campaign_complete",
            FORMAL_STAGE_A_CAMPAIGN_COMPLETE_SEAL,
            FORMAL_STAGE_A_CAMPAIGN_COMPLETE_FILE_SHA256,
        )
        require_exact_seal(
            self.result_dir / "analysis.json",
            f"{SOURCE_SCHEMA}.analysis_record",
            FORMAL_STAGE_A_ANALYSIS_SEAL,
            FORMAL_STAGE_A_ANALYSIS_FILE_SHA256,
        )
        if (expected_provenance["campaign_complete_sha256"] !=
                FORMAL_STAGE_A_CAMPAIGN_COMPLETE_FILE_SHA256 or
                expected_provenance["analysis_sha256"] !=
                FORMAL_STAGE_A_ANALYSIS_FILE_SHA256):
            die("all-K oracle formal-anchor provenance differs")
        self.module = _load_stage_a_module(
            controller, FORMAL_STAGE_A_CONTROLLER_SHA256,
        )
        self.contract = self.module.load_contract(
            self.result_dir, require_frozen_controller=False,
        )
        if canonical_json(self.contract) != canonical_json(independent):
            die("all-K oracle authenticated controller parsed another contract")
        affinity_path = self.result_dir / "controller_affinity.json"
        if sha256_file(affinity_path) != \
                expected_provenance["controller_affinity_sha256"]:
            die("all-K oracle source affinity differs from authenticated replay")
        self.allowed_cpus = set(_source_affinity(self.module, self.result_dir))
        self.stages = [
            self.module.stage_path(self.result_dir, index, 0)
            for index in range(len(CONSTRUCTION_ROOTS))
        ]
        control_paths = {
            "controller_affinity_sha256": affinity_path,
            "terminal_sha256": self.result_dir / "terminal.json",
        }
        if any(sha256_file(path) != expected_provenance[name]
               for name, path in control_paths.items()):
            die("all-K oracle control leaf differs from authenticated replay")
        expected_roots = expected_provenance.get("oh0_root_control_sha256s")
        if (not isinstance(expected_roots, list) or
                len(expected_roots) != len(CONSTRUCTION_ROOTS) or
                any(not isinstance(item, dict) or set(item) != {
                    "construction_seed_index", "manifest_sha256",
                    "complete_sha256"} or
                    item["construction_seed_index"] != index or
                    any(type(item[name]) is not str or
                        not re.fullmatch(r"[0-9a-f]{64}", item[name])
                        for name in ("manifest_sha256", "complete_sha256"))
                    for index, item in enumerate(expected_roots))):
            die("all-K oracle root control provenance is malformed")
        self.tasks_by_root = []
        for index, stage in enumerate(self.stages):
            expected_root = expected_roots[index]
            if (sha256_file(stage / "manifest.json") !=
                    expected_root["manifest_sha256"] or
                    sha256_file(stage / "complete.json") !=
                    expected_root["complete_sha256"]):
                die("all-K oracle OH0 root leaf differs from authenticated replay")
            tasks = self.module.load_manifest(stage, index, 0)
            if len(tasks) != 2 * ALLK_TASKS_PER_CONSTRUCTION_ROOT:
                die("all-K oracle Stage-A OH0 root task count mismatch")
            self.tasks_by_root.append(tasks)
        if any(sha256_file(path) != expected_provenance[name]
               for name, path in control_paths.items()):
            die("all-K oracle control leaf changed while loading its manifest")
        for index, stage in enumerate(self.stages):
            if (sha256_file(stage / "manifest.json") !=
                    expected_roots[index]["manifest_sha256"] or
                    sha256_file(stage / "complete.json") !=
                    expected_roots[index]["complete_sha256"]):
                die("all-K oracle root control leaf changed while loading")
        expected_shard_epoch = expected_provenance.get(
            "oh0_shard_epoch_sha256s",
        )
        if (not isinstance(expected_shard_epoch, Sequence) or
                isinstance(expected_shard_epoch, (str, bytes)) or
                len(expected_shard_epoch) != 2 * ALLK_TASKS or
                any(type(value) is not str or
                    not re.fullmatch(r"[0-9a-f]{64}", value)
                    for value in expected_shard_epoch) or
                expected_provenance.get("oh0_shard_epoch_root_sha256") !=
                    sha256_bytes(canonical_json(expected_shard_epoch).encode()) or
                expected_provenance.get("paired_cells_replayed") !=
                    SOURCE_PAIRED_CELLS or
                expected_provenance.get("union_cells_replayed") !=
                    EXPECTED_SOURCE_UNION or
                expected_provenance.get("formal_completion_authenticated") is not
                    True):
            die("all-K oracle source shard epoch is malformed")
        self.expected_shard_epoch = tuple(expected_shard_epoch)

    def for_task(
        self, task: Task,
    ) -> dict[CellKey, dict[str, Any]]:
        validate_task(task)
        if task.phase != "allk" or task.job not in range(ALLK_TASKS):
            die("all-K oracle received a non-all-K task")
        construction_index = task.construction_seed_index
        local_job = task.job - (
            construction_index * ALLK_TASKS_PER_CONSTRUCTION_ROOT
        )
        if local_job not in range(ALLK_TASKS_PER_CONSTRUCTION_ROOT):
            die("all-K task job does not belong to its construction root")
        stage = self.stages[construction_index]
        tasks = self.tasks_by_root[construction_index]
        left, right = tasks[2 * local_job:2 * local_job + 2]
        if (left.arm != "h12" or right.arm != "h13" or
                left.overhead != 0 or right.overhead != 0 or
                left.pair != local_job or right.pair != local_job or
                left.construction_index != construction_index or
                right.construction_index != construction_index or
                left.construction_seed != int(task.construction_seed, 16) or
                right.construction_seed != int(task.construction_seed, 16) or
                left.seed_index != task.packet_trace_root_index or
                right.seed_index != task.packet_trace_root_index or
                left.seed != task.packet_trace_root or
                right.seed != task.packet_trace_root or
                left.schedule != task.schedule or right.schedule != task.schedule or
                tuple(left.ks) != task.ks or tuple(right.ks) != task.ks):
            die("all-K task does not align with its Stage-A OH0 oracle shards")
        epoch_indices = (2 * task.job, 2 * task.job + 1)
        before = (
            source_shard_epoch_sha256(self.module, stage, left),
            source_shard_epoch_sha256(self.module, stage, right),
        )
        expected_epoch = tuple(
            self.expected_shard_epoch[index] for index in epoch_indices
        )
        if before != expected_epoch:
            die("all-K oracle source shard differs from authenticated replay epoch")
        left_rows = self.module.validate_shard(
            stage, left, Path(self.contract["binary"]),
            self.contract["binary_sha256"], Path(self.contract["taskset"]),
            self.allowed_cpus,
        )
        right_rows = self.module.validate_shard(
            stage, right, Path(self.contract["binary"]),
            self.contract["binary_sha256"], Path(self.contract["taskset"]),
            self.allowed_cpus,
        )
        after = (
            source_shard_epoch_sha256(self.module, stage, left),
            source_shard_epoch_sha256(self.module, stage, right),
        )
        if after != expected_epoch:
            die("all-K oracle source shard changed while it was consumed")
        if len(left_rows) != len(task.ks) or len(right_rows) != len(task.ks):
            die("all-K oracle source shard row count mismatch")
        result: dict[CellKey, dict[str, Any]] = {}
        for K, h12, h13 in zip(task.ks, left_rows, right_rows):
            key: CellKey = (
                construction_index, task.packet_trace_root_index,
                task.schedule, K,
            )
            arms = {"h12": h12, "h13": h13}
            self.module.validate_pair_receipt(key, arms)
            result[key] = {"role": "allk", "match_id": task.job,
                           "source": source_cell_record(key, arms)}
        return result


def expected_cells_for_task(
    source: Any, task: Task,
) -> Mapping[CellKey, Mapping[str, Any]]:
    if isinstance(source, Mapping):
        return source
    loader = getattr(source, "for_task", None)
    if loader is None:
        die("task execution lacks an authenticated source-cell oracle")
    cells = loader(task)
    if not isinstance(cells, Mapping):
        die("source-cell oracle returned a non-mapping")
    return cells


def external_file_mark(path: Path) -> dict[str, Any]:
    try:
        metadata = os.stat(path, follow_symlinks=False)
    except OSError as exc:
        die(f"cannot stat external telemetry log {path}: {exc}")
    if not stat.S_ISREG(metadata.st_mode) or metadata.st_nlink != 1:
        die("external telemetry log is not a unique regular file")
    return {
        "path": str(path), "dev": metadata.st_dev, "ino": metadata.st_ino,
        "offset": metadata.st_size, "mtime_ns": metadata.st_mtime_ns,
    }


def validate_external_mark(mark: Any, expected_path: str, context: str) -> None:
    if (not isinstance(mark, dict) or set(mark) != {
            "path", "dev", "ino", "offset", "mtime_ns"} or
            mark.get("path") != expected_path or
            any(type(mark.get(name)) is not int for name in
                ("dev", "ino", "offset", "mtime_ns")) or
            mark["dev"] < 0 or mark["ino"] <= 0 or mark["offset"] < 0):
        die(f"{context}: external telemetry mark is malformed")


def external_range_bytes(
    path: Path, start: int, end: int,
    expected_identity: Optional[tuple[int, int]] = None,
) -> bytes:
    if (type(start) is not int or type(end) is not int or start < 0 or
            end < start or end - start > 256 << 20):
        die("external telemetry interval is malformed or too large")
    flags = (os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
             getattr(os, "O_NOFOLLOW", 0))
    try:
        descriptor = os.open(path, flags)
    except OSError as exc:
        die(f"cannot open external telemetry interval: {exc}")
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode) or before.st_nlink != 1 or \
                before.st_size < end:
            die("external telemetry interval is truncated")
        chunks: list[bytes] = []
        offset = start
        while offset < end:
            block = os.pread(descriptor, min(1 << 20, end - offset), offset)
            if not block:
                die("short read from external telemetry interval")
            chunks.append(block)
            offset += len(block)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    try:
        named = os.stat(path, follow_symlinks=False)
    except OSError as exc:
        die(f"external telemetry pathname disappeared: {exc}")
    identity = lambda value: (
        value.st_dev, value.st_ino, value.st_mode, value.st_nlink,
    )
    if (identity(before) != identity(after) or identity(after) != identity(named) or
            (expected_identity is not None and
             (before.st_dev, before.st_ino) != expected_identity) or
            not stat.S_ISREG(named.st_mode) or named.st_nlink != 1 or
            named.st_size < end):
        die("external telemetry changed identity while reading")
    return b"".join(chunks)


def _telemetry_scalar(text: str, context: str) -> float:
    if not re.fullmatch(r"[0-9]+(?:\.[0-9]+)?", text):
        die(f"{context}: noncanonical telemetry scalar")
    value = float(text)
    if not math.isfinite(value):
        die(f"{context}: nonfinite telemetry scalar")
    return value


def _telemetry_counter(text: str, context: str) -> int:
    if (len(text) > 20 or
            not re.fullmatch(r"0|[1-9][0-9]*", text)):
        die(f"{context}: noncanonical telemetry counter")
    value = int(text)
    if value > MASK64:
        die(f"{context}: telemetry counter exceeds u64")
    return value


def system_uptime_s() -> float:
    try:
        fields = Path("/proc/uptime").read_text(encoding="ascii").split()
    except OSError as exc:
        die(f"cannot read canonical system uptime: {exc}")
    if len(fields) != 2:
        die("canonical system uptime has the wrong schema")
    return _telemetry_scalar(fields[0], "system uptime")


def parse_telemetry_row(line: bytes, context: str) -> dict[str, Any]:
    if len(line) > 16384 or not line.endswith(b"\n"):
        die(f"{context}: incomplete or oversized telemetry row")
    try:
        text = line.decode("utf-8")[:-1]
    except UnicodeDecodeError as exc:
        die(f"{context}: telemetry row is not UTF-8: {exc}")
    if text.endswith("\r"):
        text = text[:-1]
    values = text.split(",")
    if len(values) != len(THERMAL_FIELDS):
        die(f"{context}: telemetry field count mismatch")
    row = dict(zip(THERMAL_FIELDS, values))
    if not row["utc"] or any(row[name] == "" for name in THERMAL_DIMM_FIELDS):
        die(f"{context}: missing UTC or DIMM temperature")
    monotonic = _telemetry_scalar(row["monotonic_s"], f"{context}:monotonic")
    busy = _telemetry_scalar(row["cpu_busy_pct"], f"{context}:busy")
    mhz = _telemetry_scalar(row["cpu_avg_mhz"], f"{context}:mhz")
    cpu_temp = _telemetry_scalar(row["cpu_tctl_c"], f"{context}:cpu_temp")
    dimm = [
        _telemetry_scalar(row[name], f"{context}:{name}")
        for name in THERMAL_DIMM_FIELDS
    ]
    dimm_errors = _telemetry_counter(
        row["dimm_read_errors"], f"{context}:dimm_errors",
    )
    loads = [
        _telemetry_scalar(row[name], f"{context}:{name}")
        for name in ("load1", "load5", "load15")
    ]
    edac_ce = _telemetry_counter(row["edac_ce"], f"{context}:edac_ce")
    edac_ue = _telemetry_counter(row["edac_ue"], f"{context}:edac_ue")
    if (not 0 <= busy <= 100 or mhz <= 0 or not 0 <= cpu_temp <= 120 or
            any(not 0 <= value <= 100 for value in dimm)):
        die(f"{context}: telemetry value outside physical bounds")
    return {
        "utc": row["utc"], "monotonic_s": monotonic,
        "cpu_busy_pct": busy, "cpu_avg_mhz": mhz, "cpu_tctl_c": cpu_temp,
        "dimm_temperatures": dimm, "dimm_read_errors": dimm_errors,
        "loads": loads, "edac_ce": edac_ce, "edac_ue": edac_ue,
    }


def telemetry_baseline(
    path: Path, end: int, expected_identity: tuple[int, int],
) -> dict[str, Any]:
    header = (",".join(THERMAL_FIELDS) + "\n").encode("ascii")
    if external_range_bytes(path, 0, len(header), expected_identity) != header:
        die("external telemetry header mismatch")
    start = max(len(header), end - 16384)
    suffix = external_range_bytes(path, start, end, expected_identity)
    if not suffix.endswith(b"\n"):
        die("external telemetry baseline ends in a partial row")
    previous = suffix[:-1].rfind(b"\n")
    row = suffix[previous + 1:]
    if not row or row == header:
        die("external telemetry has no baseline sample")
    return parse_telemetry_row(row, "telemetry baseline")


def audit_telemetry_interval(
    path: Path, start: Mapping[str, Any], end: Mapping[str, Any],
    audit_uptime_s: Optional[float] = None,
) -> dict[str, Any]:
    expected_identity = start["dev"], start["ino"]
    if (end["dev"], end["ino"]) != expected_identity:
        die("external telemetry interval endpoints differ in identity")
    suffix = external_range_bytes(
        path, start["offset"], end["offset"], expected_identity,
    )
    if not suffix or not suffix.endswith(b"\n"):
        die("external telemetry interval is empty or partial")
    samples = [start["baseline"]] + [
        parse_telemetry_row(line, f"telemetry interval row {index}")
        for index, line in enumerate(suffix.splitlines(keepends=True), 1)
    ]
    gaps = [
        current["monotonic_s"] - previous["monotonic_s"]
        for previous, current in zip(samples, samples[1:])
    ]
    if any(delta <= 0 or delta > TELEMETRY_POLICY["max_gap_seconds"]
           for delta in gaps):
        die("external telemetry is nonmonotonic or has a sample gap")
    interval = samples[1:]
    uptime = system_uptime_s() if audit_uptime_s is None else audit_uptime_s
    if type(uptime) is not float or not math.isfinite(uptime):
        die("external telemetry audit clock is malformed")
    tail_age = uptime - interval[-1]["monotonic_s"]
    if tail_age < 0 or tail_age > TELEMETRY_POLICY["max_gap_seconds"]:
        die("external telemetry sampler was not live at sealing")
    if min(row["cpu_busy_pct"] for row in interval) < \
            TELEMETRY_POLICY["min_cpu_busy_pct"]:
        die("external telemetry CPU busy floor was not maintained")
    if max(max(row["cpu_tctl_c"], *row["dimm_temperatures"])
           for row in samples) >= TELEMETRY_POLICY["max_temperature_c"]:
        die("external telemetry temperature ceiling was reached")
    if any(row["dimm_read_errors"] != 0 for row in samples):
        die("external telemetry reports a DIMM read error")
    if (len({row["edac_ce"] for row in samples}) != 1 or
            len({row["edac_ue"] for row in samples}) != 1):
        die("external telemetry EDAC counters changed")
    return {
        "samples": len(interval), "max_gap_seconds": max(gaps),
        "cpu_busy_min_pct": min(row["cpu_busy_pct"] for row in interval),
        "cpu_tctl_max_c": max(row["cpu_tctl_c"] for row in samples),
        "dimm_max_c": max(
            value for row in samples for value in row["dimm_temperatures"]),
        "dimm_read_errors_max": 0,
        "edac_ce": samples[0]["edac_ce"], "edac_ue": samples[0]["edac_ue"],
        "edac_ce_delta": 0, "edac_ue_delta": 0,
        "audit_uptime_s": uptime, "tail_age_seconds": tail_age,
        "interval_sha256": sha256_bytes(suffix),
    }


def telemetry_start(
    result_dir: Path, contract: Mapping[str, Any],
    phase: str = "screen",
) -> Optional[dict[str, Any]]:
    if phase not in ("screen", "allk"):
        die("unknown telemetry phase")
    stem = "telemetry" if phase == "screen" else "allk_telemetry"
    telemetry = contract["external_telemetry"]
    if telemetry["path"] is None:
        return None
    path = Path(telemetry["path"])
    current = external_file_mark(path)
    prepared = telemetry["prepare_mark"]
    if ((current["dev"], current["ino"]) !=
            (prepared["dev"], prepared["ino"]) or
            current["offset"] < prepared["offset"]):
        die("external telemetry rotated or truncated after prepare")
    receipt_path = result_dir / f"{stem}_start.json"
    if receipt_path.exists():
        start = load_sealed(receipt_path, f"{SCHEMA}.{stem}_start")
    else:
        start = {**current, "baseline": telemetry_baseline(
            path, current["offset"], (current["dev"], current["ino"]),
        )}
        write_sealed_once(receipt_path, f"{SCHEMA}.{stem}_start", start)
    confirmation = external_file_mark(path)
    if (not isinstance(start, dict) or set(start) != {
            "path", "dev", "ino", "offset", "mtime_ns", "baseline"} or
            start["path"] != str(path) or not isinstance(start["baseline"], dict) or
            any(type(start[name]) is not int for name in
                ("dev", "ino", "offset", "mtime_ns")) or
            (confirmation["dev"], confirmation["ino"]) !=
            (start["dev"], start["ino"]) or
            confirmation["offset"] < start["offset"] or
            canonical_json(start["baseline"]) != canonical_json(
                telemetry_baseline(
                    path, start["offset"], (start["dev"], start["ino"])))):
        die("external telemetry start receipt mismatch")
    return start


def telemetry_finish(
    result_dir: Path, contract: Mapping[str, Any],
    start: Optional[Mapping[str, Any]],
    phase: str = "screen",
) -> Optional[dict[str, Any]]:
    if phase not in ("screen", "allk"):
        die("unknown telemetry phase")
    stem = "telemetry" if phase == "screen" else "allk_telemetry"
    if start is None:
        return None
    path = Path(contract["external_telemetry"]["path"])
    receipt_path = result_dir / f"{stem}_end.json"
    if receipt_path.exists():
        end = load_sealed(receipt_path, f"{SCHEMA}.{stem}_end")
    else:
        mark = external_file_mark(path)
        if ((mark["dev"], mark["ino"]) != (start["dev"], start["ino"]) or
                mark["offset"] < start["offset"]):
            die("external telemetry rotated or truncated during screen")
        audit = audit_telemetry_interval(path, start, mark)
        end = {
            **mark, "start_offset": start["offset"],
            "interval_sha256": audit["interval_sha256"],
            "semantic_audit": audit,
        }
        write_sealed_once(receipt_path, f"{SCHEMA}.{stem}_end", end)
    if (not isinstance(end, dict) or set(end) != {
            "path", "dev", "ino", "offset", "mtime_ns", "start_offset",
            "interval_sha256", "semantic_audit"} or
            any(type(end[name]) is not int for name in
                ("dev", "ino", "offset", "mtime_ns", "start_offset")) or
            end["path"] != str(path) or end["start_offset"] != start["offset"] or
            (end["dev"], end["ino"]) != (start["dev"], start["ino"]) or
            not re.fullmatch(r"[0-9a-f]{64}", str(end["interval_sha256"])) or
            not isinstance(end["semantic_audit"], dict) or
            canonical_json(end["semantic_audit"]) != canonical_json(
                audit_telemetry_interval(
                    path, start, end,
                    end["semantic_audit"].get("audit_uptime_s"))) or
            end["interval_sha256"] !=
            end["semantic_audit"].get("interval_sha256")):
        die("external telemetry end receipt mismatch")
    return end


BENCH_HEADER = (
    "pair_id", "pair_index", "pair_order_index", "arm", "N", "bb",
    "solve_block_bytes", "overhead", "schedule", "external_seed", "loss",
    "period", "geometry", "residue_skew", "residue_schedule", "gf256_rows",
    "gf16_rows", "heavy_rows", "grouped_rows", "grouped_gf256_row_mask",
    "buckets_requested", "grouped_hash_seed", "final_h_a_columns",
    "construction_attempt", "construction_seed_policy", "construction_seed",
    "production_seed_fixups_applied", "systematic_probe_result",
    "base_matrix_seed",
    "base_peel_seed", "matrix_seed", "peel_seed", "staircase_rows",
    "dense_rows", "source_hits", "dense_identity_corner", "dense_two_anchor",
    "dense_two_anchor_phase", "field", "heavy_family", "mix_count",
    "precode_count", "packet_count", "packet_trace_seed",
    "packet_trace_sha256", "pair_class", "pair_results_equal", "result",
    "result_name", "packet_rows", "peeled_columns", "inactivated_columns",
    "residual_rows", "binary_residual_rank", "residual_rank",
    "binary_deficit", "heavy_gain", "binary_row_references",
    "binary_row_storage_bytes", "binary_adjacency_storage_bytes",
    "binary_row_storage_allocations", "binary_adjacency_storage_allocations",
    "block_xors", "block_muladds", "build_ns", "peel_ns", "project_ns",
    "residual_ns", "backsub_ns", "packet_seed_attempt", "joint_source_xors",
    "joint_marginal_xors", "joint_marginal_copies", "joint_active_deltas",
    "joint_scratch_bytes", "dual_source_columns", "intermediate_bytes",
)

UINT_FIELDS = (
    "pair_index", "pair_order_index", "N", "bb", "solve_block_bytes",
    "overhead", "external_seed", "period", "residue_skew", "gf256_rows",
    "gf16_rows", "heavy_rows", "grouped_rows", "final_h_a_columns",
    "construction_attempt", "construction_seed",
    "production_seed_fixups_applied", "systematic_probe_result",
    "staircase_rows",
    "dense_rows", "source_hits", "dense_identity_corner", "dense_two_anchor",
    "dense_two_anchor_phase", "mix_count", "precode_count", "packet_count",
    "pair_results_equal", "result", "packet_rows", "peeled_columns",
    "inactivated_columns", "residual_rows", "binary_residual_rank",
    "residual_rank", "binary_deficit", "heavy_gain", "binary_row_references",
    "binary_row_storage_bytes", "binary_adjacency_storage_bytes",
    "binary_row_storage_allocations", "binary_adjacency_storage_allocations",
    "block_xors", "block_muladds", "build_ns", "peel_ns", "project_ns",
    "residual_ns", "backsub_ns", "packet_seed_attempt", "joint_source_xors",
    "joint_marginal_xors", "joint_marginal_copies", "joint_active_deltas",
    "joint_scratch_bytes", "dual_source_columns", "intermediate_bytes",
)


def expected_preamble_items(task: Task) -> tuple[tuple[str, str], ...]:
    validate_task(task)
    return (
        ("schema", "v2"), ("pair_order", "h12,h13"),
        ("arms_per_N", "2"), ("N_count", str(len(task.ks))),
        ("N", ",".join(map(str, task.ks))), ("overhead", "0"),
        ("seed", str(int(task.packet_trace_root, 16))),
        ("seed_role", "loss-trace-root"),
        ("construction_seed_policy", "matrix-c-peel-lo32-xor-hi32-v1"),
        ("construction_seed", str(int(task.construction_seed, 16))),
        ("production_seed_fixups_applied", "0"),
        ("schedule", task.schedule),
        ("policy", "h12-h13-q0-grouped-v1"), ("period", "48"),
        ("geometry", "shared-x"), ("residue_skew", "0"),
        ("residue_schedule", "constant"), ("residue_hash_seed", "0x0"),
        ("extension_residue_seed_xor", "0x4e"),
        ("independent_extension_residues", "0"),
        ("gf256_rows", "h12:10|h13:11"), ("gf16_rows", "2"),
        ("grouped_rows", "3"), ("grouped_gf256_row_mask", "0x380"),
        ("buckets", "separate"),
        ("grouped_hash_seed", "h12:0xb7e15162|h13:0xb7e15163"),
        ("grouped_final_h_a_columns", "arm-heavy-rows"),
        ("dense_rows", "12"), ("dense_identity_corner", "0"),
        ("dense_two_anchor", "1"), ("dense_two_anchor_phase", "0"),
        ("source_hits", "canonical-K-rule"),
        ("staircase_rows", "GetDenseCount(K)"),
        ("field", "mixed-gf256-gf16"),
        ("heavy_family", "periodic-cauchy"),
        ("construction_attempt", "0"),
        ("systematic_probe", "direct-attempt0"),
        ("mix", "2"), ("bb", "64"),
        ("solve_block_bytes", "2"), ("loss", "0.5"),
        ("overhead_stream", "paired"),
        ("packet_row_seed_multiplier", "0x1"),
        ("packet_row_seed_avalanche", "0"),
        ("odd_packet_peel_seed_xor", "0x0"),
        ("packet_vector", "one-shared-per-N"),
        ("payload", "shared-zero-v1"), ("packet_payload_bytes", "2"),
        ("packet_trace_schema",
         "wirehair-wh2-precodefail-raw-packet-trace-v1"),
        ("coefficient_layout", "materialized-before-solve"),
        ("grouped_schedule_prefix", "materialized-before-solve"),
    )


def expected_preamble_line(task: Task) -> str:
    return "# groupedrecovery: " + " ".join(
        f"{name}={value}" for name, value in expected_preamble_items(task)
    )


def parse_preamble(line: str) -> dict[str, str]:
    prefix = "# groupedrecovery: "
    if not line.startswith(prefix):
        die("missing groupedrecovery schema-v2 preamble")
    values: dict[str, str] = {}
    for token in line[len(prefix):].split(" "):
        name, separator, value = token.partition("=")
        if not separator or not name or not value or name in values:
            die(f"malformed or duplicate groupedrecovery preamble token {token!r}")
        values[name] = value
    return values


def make_benchmark_argv(binary: Path, task: Task) -> list[str]:
    validate_task(task)
    return [
        str(binary), "groupedrecovery", "--N", ",".join(map(str, task.ks)),
        "--overhead", "0", "--seed", str(int(task.packet_trace_root, 16)),
        "--construction-seed", str(int(task.construction_seed, 16)),
        "--schedule", task.schedule,
    ]


def grouped_pair_id(task: Task, K: int, trace_seed: int, trace_sha256: str) -> str:
    matrix_seed, peel_seed = construction_seed_values(task.construction_seed)
    domain = (
        "wirehair-wh2-grouped-recovery-pair-v2\n"
        f"N={K}\nbb=64\nsolve_block_bytes=2\n"
        "overhead=0\nloss=0.5\n"
        f"external_seed={int(task.packet_trace_root, 16)}\n"
        "external_seed_role=loss-trace-root\n"
        "construction_seed_policy=matrix-c-peel-lo32-xor-hi32-v1\n"
        f"construction_seed={int(task.construction_seed, 16)}\n"
        "production_seed_fixups_applied=0\n"
        f"base_matrix_seed={matrix_seed}\nbase_peel_seed={peel_seed}\n"
        f"schedule={task.schedule}\n"
        "period=48\ngeometry=shared-x\nresidue_skew=0\n"
        "residue_schedule=constant\nresidue_hash_seed=0\n"
        "extension_residue_seed_xor=78\n"
        "independent_extension_residues=0\n"
        "gf256_rows=h12:10|h13:11\ngf16_rows=2\ngrouped_rows=3\n"
        "grouped_gf256_row_mask=0x380\n"
        "grouped_hash_seed=h12:0xb7e15162|h13:0xb7e15163\n"
        "buckets=separate\ndense_rows=12\ndense_identity_corner=0\n"
        "dense_two_anchor=1\ndense_two_anchor_phase=0\n"
        "source_hits=canonical-K-rule\nstaircase_rows=GetDenseCount(K)\n"
        "field=mixed-gf256-gf16\nheavy_family=periodic-cauchy\n"
        "construction_attempt=0\nmix=2\n"
        "overhead_stream=paired\npacket_row_seed_multiplier=1\n"
        "packet_row_seed_avalanche=0\nodd_packet_peel_seed_xor=0\n"
        "payload=shared-zero-v1\n"
        "packet_trace_schema=wirehair-wh2-precodefail-raw-packet-trace-v1\n"
        f"packet_trace_seed={trace_seed}\n"
        f"packet_trace_sha256={trace_sha256}\n"
    )
    return sha256_bytes(domain.encode("ascii"))


def cohort_cell_map(cohort: Mapping[str, Any]) -> dict[CellKey, dict[str, Any]]:
    matches = cohort.get("matches") if isinstance(cohort, Mapping) else None
    if not isinstance(matches, list):
        die("cohort matches are malformed")
    cells: dict[CellKey, dict[str, Any]] = {}
    stratum_fields = {
        "construction_seed_index", "construction_seed",
        "packet_trace_root_index", "packet_trace_root", "schedule",
        "band_index", "band",
    }
    source_fields = {
        "construction_seed_index", "construction_seed",
        "packet_trace_root_index", "packet_trace_root", "schedule", "K",
        "raw_failures", "raw_systematic_probe_results",
        *SOURCE_IDENTITY_FIELDS,
    }
    for expected_match_id, match in enumerate(matches):
        if (not isinstance(match, Mapping) or set(match) != {
                "match_id", "stratum", "control_selection_sha256", "union",
                "control"} or match.get("match_id") != expected_match_id):
            die("cohort match is malformed")
        stratum = match.get("stratum")
        if not isinstance(stratum, Mapping) or set(stratum) != stratum_fields:
            die("cohort match stratum schema mismatch")
        construction_seed_index = stratum["construction_seed_index"]
        packet_trace_root_index = stratum["packet_trace_root_index"]
        schedule = stratum["schedule"]
        band = stratum["band_index"]
        if (type(construction_seed_index) is not int or
                construction_seed_index not in range(len(CONSTRUCTION_ROOTS)) or
                stratum["construction_seed"] !=
                    construction_seed_text(construction_seed_index) or
                type(packet_trace_root_index) is not int or
                packet_trace_root_index not in range(len(LOSS_ROOTS)) or
                stratum["packet_trace_root"] !=
                    LOSS_ROOTS[packet_trace_root_index] or
                schedule not in SCHEDULES or type(band) is not int or
                band not in range(len(BANDS)) or
                stratum["band"] != list(BANDS[band])):
            die("cohort match stratum value mismatch")
        for role in ("union", "control"):
            source = match.get(role)
            if not isinstance(source, dict) or set(source) != source_fields:
                die("cohort source cell is malformed")
            key: CellKey = (
                source.get("construction_seed_index"),
                source.get("packet_trace_root_index"),
                source.get("schedule"), source.get("K"),
            )
            if (type(key[0]) is not int or type(key[1]) is not int or
                    type(key[2]) is not str or type(key[3]) is not int or
                    key in cells or key[0] != construction_seed_index or
                    key[1] != packet_trace_root_index or key[2] != schedule or
                    source.get("construction_seed") !=
                        construction_seed_text(construction_seed_index) or
                    source.get("packet_trace_root") !=
                        LOSS_ROOTS[packet_trace_root_index] or
                    not K_MIN <= key[3] <= K_MAX or
                    band_index(key[3]) != band):
                die("cohort source cell is duplicated or malformed")
            failures = source["raw_failures"]
            probes = source["raw_systematic_probe_results"]
            if (not isinstance(failures, dict) or set(failures) != set(ARMS) or
                    any(type(failures[arm]) is not bool for arm in ARMS) or
                    not isinstance(probes, dict) or set(probes) != set(ARMS) or
                    any(type(probes[arm]) is not int or probes[arm] not in (0, 1)
                        for arm in ARMS) or
                    (role == "union" and not any(failures.values())) or
                    (role == "control" and any(failures.values()))):
                die("cohort source outcome/probe receipt mismatch")
            for name, bits in (
                    ("base_matrix_seed", 64), ("matrix_seed", 64),
                    ("base_peel_seed", 32), ("peel_seed", 32),
                    ("packet_trace_seed", 64)):
                if type(source[name]) is not str:
                    die("cohort source seed receipt has the wrong type")
                strict_hex(source[name], f"cohort:{name}", bits)
            matrix_seed, peel_seed = construction_seed_values(
                CONSTRUCTION_ROOTS[construction_seed_index],
            )
            if (source["construction_seed_policy"] !=
                    "matrix-c-peel-lo32-xor-hi32-v1" or
                    source["production_seed_fixups_applied"] != 0 or
                    source["base_matrix_seed"] != hex(matrix_seed) or
                    source["matrix_seed"] != hex(matrix_seed) or
                    source["base_peel_seed"] != hex(peel_seed) or
                    source["peel_seed"] != hex(peel_seed) or
                    source["packet_trace_seed"] != hex(loss_seed(
                        LOSS_ROOTS[packet_trace_root_index], key[3])) or
                    type(source["packet_trace_sha256"]) is not str or
                    not re.fullmatch(
                        r"[0-9a-f]{64}", source["packet_trace_sha256"]) or
                    any(type(source[name]) is not int for name in (
                        "actual_staircase_rows", "actual_dense_rows",
                        "actual_source_hits", "actual_dense_two_anchor",
                        "actual_dense_two_anchor_phase")) or
                    source["actual_dense_rows"] != 12 or
                    source["actual_source_hits"] !=
                        (2 if key[3] < 10000 else 3) or
                    source["actual_dense_two_anchor"] != 1 or
                    source["actual_dense_two_anchor_phase"] != 0 or
                    type(source["actual_staircase_rows"]) is not int or
                    source["actual_staircase_rows"] < 1):
                die("cohort source construction/trace receipt mismatch")
            cells[key] = {
                "match_id": match["match_id"], "role": role, "source": source,
            }
        control = match["control"]
        expected_selection = _control_digest(
            construction_seed_index, packet_trace_root_index, schedule, band,
            control["K"],
        )
        if match["control_selection_sha256"] != expected_selection:
            die("cohort control selection digest mismatch")
    return cells


def _validate_row_geometry(
    row: Mapping[str, str], task: Task, K: int, arm: str,
    pair_index: int, expected: Optional[Mapping[str, Any]],
) -> None:
    context = f"K={K}/{arm}"
    parsed = {name: strict_uint(row[name], f"{context}:{name}") for name in UINT_FIELDS}
    fixed = {
        "pair_index": pair_index, "pair_order_index": ARMS.index(arm),
        "N": K, "bb": 64, "solve_block_bytes": 2, "overhead": 0,
        "external_seed": int(task.packet_trace_root, 16), "period": 48,
        "residue_skew": 0, "gf16_rows": 2,
        "gf256_rows": ARM_GF256_ROWS[arm],
        "heavy_rows": ARM_HEAVY_ROWS[arm], "grouped_rows": 3,
        "final_h_a_columns": ARM_HEAVY_ROWS[arm],
        "construction_attempt": 0,
        "construction_seed": int(task.construction_seed, 16),
        "production_seed_fixups_applied": 0, "dense_rows": 12,
        "source_hits": 2 if K < 10000 else 3,
        "dense_identity_corner": 0, "dense_two_anchor": 1,
        "dense_two_anchor_phase": 0, "mix_count": 2,
        "packet_count": K, "packet_rows": K, "packet_seed_attempt": 0,
    }
    for name, value in fixed.items():
        if parsed[name] != value:
            die(f"{context}: grouped recovery {name} geometry mismatch")
    text_fixed = {
        "schedule": task.schedule, "loss": "0.5", "geometry": "shared-x",
        "residue_schedule": "constant", "buckets_requested": "separate",
        "field": "mixed-gf256-gf16", "heavy_family": "periodic-cauchy",
        "construction_seed_policy": "matrix-c-peel-lo32-xor-hi32-v1",
    }
    for name, value in text_fixed.items():
        if row[name] != value:
            die(f"{context}: grouped recovery {name} text mismatch")
    if row["grouped_gf256_row_mask"] != "0x380":
        die(f"{context}: grouped row mask mismatch")
    if row["grouped_hash_seed"] != ARM_GROUP_HASH_SEED[arm]:
        die(f"{context}: grouped hash seed mismatch")
    if parsed["systematic_probe_result"] not in (0, 1):
        die(f"{context}: systematic probe result is not success/need-more")
    result = parsed["result"]
    if result not in (0, 1) or row["result_name"] != (
            "success" if result == 0 else "need-more"):
        die(f"{context}: codec result/result_name mismatch")
    for name, bits in (
            ("base_matrix_seed", 64), ("matrix_seed", 64),
            ("base_peel_seed", 32), ("peel_seed", 32),
            ("packet_trace_seed", 64)):
        strict_hex(row[name], f"{context}:{name}", bits)
    matrix_seed, peel_seed = construction_seed_values(task.construction_seed)
    if (row["base_matrix_seed"] != hex(matrix_seed) or
            row["matrix_seed"] != hex(matrix_seed) or
            row["base_peel_seed"] != hex(peel_seed) or
            row["peel_seed"] != hex(peel_seed)):
        die(f"{context}: seeds are not exact uniform-root v1")
    if expected is not None:
        source = expected["source"]
        source_names = {
            "base_matrix_seed": "base_matrix_seed",
            "base_peel_seed": "base_peel_seed",
            "matrix_seed": "matrix_seed", "peel_seed": "peel_seed",
            "staircase_rows": "actual_staircase_rows",
            "dense_rows": "actual_dense_rows",
            "source_hits": "actual_source_hits",
            "dense_two_anchor": "actual_dense_two_anchor",
            "dense_two_anchor_phase": "actual_dense_two_anchor_phase",
            "packet_trace_seed": "packet_trace_seed",
            "packet_trace_sha256": "packet_trace_sha256",
        }
        for output_name, source_name in source_names.items():
            if str(row[output_name]) != str(source[source_name]):
                die(
                    f"{context}: {output_name} differs from authenticated "
                    "Stage-A OH0"
                )
    if (not re.fullmatch(r"[0-9a-f]{64}", row["packet_trace_sha256"]) or
            int(row["packet_trace_seed"], 16) !=
                loss_seed(task.packet_trace_root, K)):
        die(f"{context}: packet trace receipt mismatch")
    if parsed["precode_count"] != (
            parsed["staircase_rows"] + parsed["dense_rows"] +
            parsed["heavy_rows"]):
        die(f"{context}: precode count arithmetic mismatch")
    total_columns = K + parsed["precode_count"]
    if (parsed["peeled_columns"] + parsed["inactivated_columns"] !=
            total_columns):
        die(f"{context}: peel column partition mismatch")
    if not (parsed["binary_residual_rank"] <= parsed["residual_rank"] <=
            parsed["inactivated_columns"]):
        die(f"{context}: residual rank ordering mismatch")
    if (parsed["residual_rank"] > parsed["residual_rows"] or
            parsed["residual_rows"] > total_columns):
        die(f"{context}: residual row/rank accounting mismatch")
    if parsed["binary_deficit"] != (
            parsed["inactivated_columns"] - parsed["binary_residual_rank"]):
        die(f"{context}: binary deficit arithmetic mismatch")
    if parsed["heavy_gain"] != (
            parsed["residual_rank"] - parsed["binary_residual_rank"]):
        die(f"{context}: heavy gain arithmetic mismatch")
    if parsed["heavy_gain"] > parsed["heavy_rows"]:
        die(f"{context}: heavy gain exceeds the heavy-row budget")
    if (result == 0) != (
            parsed["residual_rank"] == parsed["inactivated_columns"]):
        die(f"{context}: result differs from exact residual full rank")
    for name in (
            "joint_source_xors", "joint_marginal_xors",
            "joint_marginal_copies", "joint_active_deltas",
            "joint_scratch_bytes", "dual_source_columns"):
        if parsed[name] != 0:
            die(f"{context}: separate-bucket receipt has nonzero {name}")
    expected_intermediate = (
        2 * (K + parsed["precode_count"]) if result == 0 else 0
    )
    if parsed["intermediate_bytes"] != expected_intermediate:
        die(f"{context}: intermediate byte count mismatch")


def validate_output_pair(
    task: Task, K: int, pair_index: int,
    arms: Mapping[str, Mapping[str, str]],
    expected: Optional[Mapping[str, Any]],
) -> None:
    if set(arms) != set(ARMS):
        die(f"K={K}: output is not an exact H12/H13 pair")
    for arm in ARMS:
        _validate_row_geometry(arms[arm], task, K, arm, pair_index, expected)
    h12, h13 = arms["h12"], arms["h13"]
    h12_result, h13_result = int(h12["result"]), int(h13["result"])
    pair_class = (
        "both-success" if h12_result == 0 and h13_result == 0 else
        "h12-only" if h12_result == 0 else
        "h13-only" if h13_result == 0 else "both-need-more"
    )
    equal = "1" if h12_result == h13_result else "0"
    shared = (
        "pair_id", "pair_index", "N", "bb", "solve_block_bytes", "overhead",
        "schedule", "external_seed", "loss", "period", "geometry",
        "residue_skew", "residue_schedule", "gf16_rows", "grouped_rows",
        "grouped_gf256_row_mask", "buckets_requested", "construction_attempt",
        "construction_seed_policy", "construction_seed",
        "production_seed_fixups_applied",
        "base_matrix_seed", "base_peel_seed", "matrix_seed", "peel_seed",
        "staircase_rows", "dense_rows", "source_hits", "dense_identity_corner",
        "dense_two_anchor", "dense_two_anchor_phase", "field", "heavy_family",
        "mix_count", "packet_count", "packet_trace_seed", "packet_trace_sha256",
        "pair_class", "pair_results_equal", "packet_rows", "packet_seed_attempt",
    )
    for name in shared:
        if h12[name] != h13[name]:
            die(f"K={K}: paired receipt {name} differs between arms")
    if h12["pair_class"] != pair_class or h12["pair_results_equal"] != equal:
        die(f"K={K}: pair class/equality receipt mismatch")
    if int(h13["precode_count"]) != int(h12["precode_count"]) + 1:
        die(f"K={K}: H13 precode count is not H12 + 1")
    trace_seed = int(h12["packet_trace_seed"], 16)
    expected_pair_id = grouped_pair_id(
        task, K, trace_seed, h12["packet_trace_sha256"],
    )
    if (h12["pair_id"] != expected_pair_id or
            not re.fullmatch(r"[0-9a-f]{64}", h12["pair_id"])):
        die(f"K={K}: pair_id does not authenticate its declared pair domain")


def parse_output(
    output: str, task: Task,
    expected_cells: Mapping[CellKey, Mapping[str, Any]],
) -> list[dict[str, dict[str, str]]]:
    if "\r" in output or not output.endswith("\n"):
        die("benchmark output is not canonical LF-terminated text")
    lines = output.splitlines()
    if len(lines) != 2 + 2 * len(task.ks):
        die("benchmark output line count differs from pair-native task")
    expected_line = expected_preamble_line(task)
    if (lines[0] != expected_line or
            parse_preamble(lines[0]) != dict(expected_preamble_items(task))):
        die("benchmark preamble differs from the sealed Stage-B task")
    try:
        header = tuple(next(csv.reader([lines[1]], strict=True)))
    except csv.Error as exc:
        die(f"malformed grouped recovery CSV header: {exc}")
    if lines[1] != ",".join(BENCH_HEADER) or header != BENCH_HEADER:
        die("benchmark CSV header differs from groupedrecovery schema v2")
    pairs: list[dict[str, dict[str, str]]] = []
    seen_pair_ids: set[str] = set()
    for pair_index, K in enumerate(task.ks):
        expected = expected_cells.get((
            task.construction_seed_index, task.packet_trace_root_index,
            task.schedule, K,
        ))
        if expected is None:
            die(f"K={K}: task cell is absent from its authenticated source oracle")
        arms: dict[str, dict[str, str]] = {}
        for order_index, line in enumerate(lines[2 + 2 * pair_index:
                                                 4 + 2 * pair_index]):
            try:
                values = next(csv.reader([line], strict=True))
            except csv.Error as exc:
                die(f"malformed grouped recovery row: {exc}")
            if line != ",".join(values) or len(values) != len(BENCH_HEADER):
                die("grouped recovery row has noncanonical quoting/field count")
            row = dict(zip(BENCH_HEADER, values))
            arm = ARMS[order_index]
            if row["arm"] != arm or row["pair_order_index"] != str(order_index):
                die(f"K={K}: pair rows are not adjacent H12 then H13")
            if arm in arms:
                die(f"K={K}: duplicate output arm")
            arms[arm] = row
        validate_output_pair(task, K, pair_index, arms, expected)
        pair_id = arms["h12"]["pair_id"]
        if pair_id in seen_pair_ids:
            die("one task contains a duplicate pair_id")
        seen_pair_ids.add(pair_id)
        pairs.append(arms)
    return pairs


def _copy_unique(
    source: Path, destination: Path, executable: bool,
    *, allow_source_hardlinks: bool = False,
) -> str:
    if destination.exists():
        die(f"freeze target already exists: {destination}")
    data = stable_source_bytes(
        source, require_unique=not allow_source_hardlinks,
    )
    atomic_write(destination, data)
    os.chmod(destination, 0o555 if executable else 0o444)
    digest = sha256_bytes(data)
    if sha256_file(destination) != digest:
        die(f"frozen copy hash mismatch: {destination}")
    return digest


def prepare(args: argparse.Namespace) -> None:
    require_formal_stage_a_launchable()
    result_dir = args.result_dir.resolve()
    if result_dir.exists():
        die("result directory already exists")
    binary_source = args.binary.resolve(strict=True)
    controller_source = Path(__file__).resolve(strict=True)
    source_result = args.stage_a_result.resolve(strict=True)
    taskset_text = shutil.which("taskset")
    if not taskset_text:
        die("taskset is required")
    taskset_source = Path(taskset_text).resolve(strict=True)
    telemetry_path = (
        args.telemetry_log.resolve(strict=True) if args.telemetry_log else None
    )
    if (checked_product(
            K_COUNT, len(CONSTRUCTION_ROOTS), len(LOSS_ROOTS), len(SCHEDULES)) !=
            SOURCE_PAIRED_CELLS or
            checked_product(SOURCE_PAIRED_CELLS, len(ARMS)) !=
            SOURCE_ARM_OUTCOMES or ALLK_TASKS != 6912 or
            ALLK_CELLS != 1727973):
        die("Stage-B cardinality constants are inconsistent")

    # Expensive source validation deliberately precedes creation of any target
    # artifact.  A failed/partial Stage-A result can never leave something that
    # resembles a prepared Stage-B campaign.
    cohort, source_provenance = authenticate_stage_a(source_result)
    tasks = build_screen_tasks(cohort)
    if (len(tasks) != EXPECTED_SCREEN_TASKS or
            sum(len(task.ks) for task in tasks) != EXPECTED_SCREEN_CELLS):
        die("formal screen plan cardinality mismatch")
    telemetry_contract = {
        "managed_by_controller": False,
        "path": str(telemetry_path) if telemetry_path else None,
        "prepare_mark": external_file_mark(telemetry_path)
            if telemetry_path else None,
        "policy": TELEMETRY_POLICY,
        "sampler_identity": (
            "external operator responsibility; CPU and all eight DRAM TSOD "
            "temperature channels are sampled by the authenticated log worker"
        ),
        "continuity": TELEMETRY_CONTINUITY,
    }
    staging = result_dir.with_name(f".{result_dir.name}.prepare-{os.getpid()}")
    if staging.exists():
        die("prepare staging directory already exists")
    staging.mkdir(parents=True)
    try:
        frozen = staging / "frozen"
        frozen_binary = frozen / "wirehair_v2_bench"
        frozen_controller = frozen / Path(__file__).name
        frozen_taskset = frozen / "taskset"
        binary_sha256 = _copy_unique(binary_source, frozen_binary, True)
        controller_sha256 = _copy_unique(
            controller_source, frozen_controller, True,
        )
        taskset_sha256 = _copy_unique(
            taskset_source, frozen_taskset, True,
            allow_source_hardlinks=True,
        )
        write_sealed_once(
            staging / "cohort.json", f"{SCHEMA}.cohort", cohort,
        )
        cohort_sha256 = sha256_file(staging / "cohort.json")
        final_frozen = result_dir / "frozen"
        contract = {
            "K_domain": [K_MIN, K_MAX],
            "construction_roots": list(CONSTRUCTION_ROOTS),
            "construction_root_sets": {
                name: list(indices)
                for name, indices in STAGE_A_ROOT_SETS.items()
            },
            "packet_trace_roots": list(LOSS_ROOTS),
            "schedules": list(SCHEDULES),
            "bands": [list(value) for value in BANDS],
            "source_stage_a": source_provenance,
            "source_controller_sha256": FORMAL_STAGE_A_CONTROLLER_SHA256,
            "source_binary_sha256": FORMAL_STAGE_A_BINARY_SHA256,
            "source_campaign_complete_seal":
                FORMAL_STAGE_A_CAMPAIGN_COMPLETE_SEAL,
            "source_campaign_complete_file_sha256":
                FORMAL_STAGE_A_CAMPAIGN_COMPLETE_FILE_SHA256,
            "source_analysis_seal": FORMAL_STAGE_A_ANALYSIS_SEAL,
            "source_analysis_file_sha256":
                FORMAL_STAGE_A_ANALYSIS_FILE_SHA256,
            "cohort_sha256": cohort_sha256,
            "selection_rule": SELECTION_RULE,
            "screen": {
                "matches": EXPECTED_SCREEN_MATCHES,
                "union_cells": EXPECTED_SOURCE_UNION,
                "matched_controls": EXPECTED_SCREEN_MATCHES,
                "paired_cells": EXPECTED_SCREEN_CELLS,
                "arm_outcomes": EXPECTED_SCREEN_OUTCOMES,
                "nonempty_strata": EXPECTED_NONEMPTY_STRATA,
                "zero_strata": [
                    [construction_seed_index, packet_trace_root_index,
                     schedule, band]
                    for (construction_seed_index, packet_trace_root_index,
                         schedule, band) in EXPECTED_ZERO_STRATA
                ],
                "matches_per_task_max": SCREEN_MATCH_CHUNK,
                "pair_native_tasks": EXPECTED_SCREEN_TASKS,
            },
            "allk_escalation": {
                "conditional_plan_only": True, "chunk_size": CHUNK_SIZE,
                "pair_native_tasks": ALLK_TASKS,
                "paired_cells": ALLK_CELLS, "arm_outcomes": ALLK_OUTCOMES,
            },
            "architecture": architecture_contract(),
            "promotion_rule": {
                "mechanical_validation": "absolute",
                "union_scope": "H13 failures strictly less than H12",
                "full_matched_scope": "H13 failures strictly less than H12",
                "mcnemar_and_strata": "descriptive-only",
                "single_cell_regressions": "reported-not-vetoed",
            },
            "child_lifetime_policy": (
                "dedicated process group; registered kill/reap on SIGINT or "
                "SIGTERM; Linux PR_SET_PDEATHSIG=SIGKILL"
            ),
            "binary": str(final_frozen / frozen_binary.name),
            "binary_sha256": binary_sha256,
            "controller": str(final_frozen / frozen_controller.name),
            "controller_sha256": controller_sha256,
            "taskset": str(final_frozen / frozen_taskset.name),
            "taskset_sha256": taskset_sha256,
            "result_dir": str(result_dir),
            "loss_seed_formula": LOSS_SEED_FORMULA,
            "external_telemetry": telemetry_contract,
        }
        write_sealed_once(staging / "contract.json", f"{SCHEMA}.contract", contract)
        write_sealed_once(
            staging / "stages" / "screen" / "manifest.json",
            f"{SCHEMA}.stage_manifest", manifest_payload("screen", tasks),
        )
        staging.rename(result_dir)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    print(canonical_json({
        "result_dir": str(result_dir), "screen_tasks": len(tasks),
        "screen_paired_cells": EXPECTED_SCREEN_CELLS,
        "stage_a_identity_stream_sha256":
            cohort["source_identity_stream_sha256"],
        "run_command": [
            str(result_dir / "frozen" / Path(__file__).name), "run",
            "--result-dir", str(result_dir), "--cpus", "<CPU-LIST>",
        ],
    }))


def load_cohort(result_dir: Path, expected_sha256: str) -> dict[str, Any]:
    path = result_dir / "cohort.json"
    if sha256_file(path) != expected_sha256:
        die("sealed cohort file hash differs from the contract")
    cohort = load_sealed(path, f"{SCHEMA}.cohort")
    if (not isinstance(cohort, dict) or set(cohort) != {
            "selection_schema", "selection_rule", "source_union_count",
            "matched_control_count", "paired_cells", "arm_outcomes", "bands",
            "nonempty_strata", "zero_strata", "source_identity_stream_sha256",
            "counts_by_stratum_band", "matches"} or
            cohort["selection_schema"] != f"{SCHEMA}.matched_control.v1" or
            cohort["selection_rule"] != SELECTION_RULE or
            cohort["source_union_count"] != EXPECTED_SOURCE_UNION or
            cohort["matched_control_count"] != EXPECTED_SCREEN_MATCHES or
            cohort["paired_cells"] != EXPECTED_SCREEN_CELLS or
            cohort["arm_outcomes"] != EXPECTED_SCREEN_OUTCOMES or
            cohort["bands"] != [list(value) for value in BANDS] or
            cohort["nonempty_strata"] != EXPECTED_NONEMPTY_STRATA or
            not re.fullmatch(
                r"[0-9a-f]{64}",
                str(cohort["source_identity_stream_sha256"])) or
            not isinstance(cohort["matches"], list) or
            len(cohort["matches"]) != EXPECTED_SCREEN_MATCHES or
            not isinstance(cohort["counts_by_stratum_band"], list) or
            len(cohort["counts_by_stratum_band"]) != EXPECTED_NONEMPTY_STRATA or
            not isinstance(cohort["zero_strata"], list) or
            len(cohort["zero_strata"]) != len(EXPECTED_ZERO_STRATA)):
        die("sealed cohort formal invariants mismatch")
    cells = cohort_cell_map(cohort)
    if len(cells) != EXPECTED_SCREEN_CELLS:
        die("sealed cohort cell keys are not exact/disjoint")
    expected_zero = [
        {"construction_seed_index": construction_seed_index,
         "construction_seed": construction_seed_text(construction_seed_index),
         "packet_trace_root_index": packet_trace_root_index,
         "packet_trace_root": LOSS_ROOTS[packet_trace_root_index],
         "schedule": schedule, "band_index": band,
         "band": list(BANDS[band])}
        for construction_seed_index in range(len(CONSTRUCTION_ROOTS))
        for packet_trace_root_index in range(len(LOSS_ROOTS))
        for schedule in SCHEDULES
        for band in range(len(BANDS))
        if (construction_seed_index, packet_trace_root_index,
            schedule, band) in set(EXPECTED_ZERO_STRATA)
    ]
    if canonical_json(cohort["zero_strata"]) != canonical_json(expected_zero):
        die("sealed cohort zero-stratum identities mismatch")
    counts: dict[StratumKey, int] = {}
    for match in cohort["matches"]:
        stratum = match["stratum"]
        key = (
            stratum["construction_seed_index"],
            stratum["packet_trace_root_index"], stratum["schedule"],
            stratum["band_index"],
        )
        counts[key] = counts.get(key, 0) + 1
    expected_counts = [
        {"construction_seed_index": construction_seed_index,
         "construction_seed": construction_seed_text(construction_seed_index),
         "packet_trace_root_index": packet_trace_root_index,
         "packet_trace_root": LOSS_ROOTS[packet_trace_root_index],
         "schedule": schedule, "band_index": band,
         "band": list(BANDS[band]), "union": count, "controls": count}
        for (construction_seed_index, packet_trace_root_index,
             schedule, band), count in sorted(
            counts.items(), key=lambda item: _stratum_order(item[0]))
    ]
    if canonical_json(cohort["counts_by_stratum_band"]) != \
            canonical_json(expected_counts):
        die("sealed cohort count-by-stratum receipts mismatch")
    if len(build_screen_tasks(cohort)) != EXPECTED_SCREEN_TASKS:
        die("sealed cohort does not reproduce the pinned screen-task plan")
    return cohort


def load_contract(
    result_dir: Path, *, require_frozen_controller: bool,
) -> tuple[dict[str, Any], dict[str, Any]]:
    require_formal_stage_a_launchable()
    contract = load_sealed(result_dir / "contract.json", f"{SCHEMA}.contract")
    expected_keys = {
        "K_domain", "construction_roots", "construction_root_sets",
        "packet_trace_roots", "schedules", "bands",
        "source_stage_a", "source_controller_sha256", "source_binary_sha256",
        "source_campaign_complete_seal",
        "source_campaign_complete_file_sha256", "source_analysis_seal",
        "source_analysis_file_sha256",
        "cohort_sha256", "selection_rule", "screen", "allk_escalation",
        "architecture", "promotion_rule", "child_lifetime_policy", "binary",
        "binary_sha256", "controller", "controller_sha256", "taskset",
        "taskset_sha256", "result_dir", "loss_seed_formula",
        "external_telemetry",
    }
    if (not isinstance(contract, dict) or set(contract) != expected_keys or
            contract["K_domain"] != [K_MIN, K_MAX] or
            contract["construction_roots"] != list(CONSTRUCTION_ROOTS) or
            canonical_json(contract["construction_root_sets"]) !=
                canonical_json({
                    name: list(indices)
                    for name, indices in STAGE_A_ROOT_SETS.items()
                }) or
            contract["packet_trace_roots"] != list(LOSS_ROOTS) or
            contract["schedules"] != list(SCHEDULES) or
            contract["bands"] != [list(value) for value in BANDS] or
            contract["source_controller_sha256"] !=
                FORMAL_STAGE_A_CONTROLLER_SHA256 or
            contract["source_binary_sha256"] != FORMAL_STAGE_A_BINARY_SHA256 or
            contract["source_campaign_complete_seal"] !=
                FORMAL_STAGE_A_CAMPAIGN_COMPLETE_SEAL or
            contract["source_campaign_complete_file_sha256"] !=
                FORMAL_STAGE_A_CAMPAIGN_COMPLETE_FILE_SHA256 or
            contract["source_analysis_seal"] != FORMAL_STAGE_A_ANALYSIS_SEAL or
            contract["source_analysis_file_sha256"] !=
                FORMAL_STAGE_A_ANALYSIS_FILE_SHA256 or
            contract["selection_rule"] != SELECTION_RULE or
            contract["result_dir"] != str(result_dir) or
            contract["loss_seed_formula"] != LOSS_SEED_FORMULA or
            any(type(contract.get(name)) is not str for name in (
                "binary", "binary_sha256", "controller", "controller_sha256",
                "taskset", "taskset_sha256", "cohort_sha256")) or
            any(not re.fullmatch(r"[0-9a-f]{64}", contract[name]) for name in (
                "binary_sha256", "controller_sha256", "taskset_sha256",
                "cohort_sha256"))):
        die("Stage-B campaign contract schema/identity mismatch")
    fixed_screen = {
        "matches": EXPECTED_SCREEN_MATCHES,
        "union_cells": EXPECTED_SOURCE_UNION,
        "matched_controls": EXPECTED_SCREEN_MATCHES,
        "paired_cells": EXPECTED_SCREEN_CELLS,
        "arm_outcomes": EXPECTED_SCREEN_OUTCOMES,
        "nonempty_strata": EXPECTED_NONEMPTY_STRATA,
        "zero_strata": [list(value) for value in EXPECTED_ZERO_STRATA],
        "matches_per_task_max": SCREEN_MATCH_CHUNK,
        "pair_native_tasks": EXPECTED_SCREEN_TASKS,
    }
    fixed_allk = {
        "conditional_plan_only": True, "chunk_size": CHUNK_SIZE,
        "pair_native_tasks": ALLK_TASKS, "paired_cells": ALLK_CELLS,
        "arm_outcomes": ALLK_OUTCOMES,
    }
    fixed_architecture = architecture_contract()
    fixed_promotion = {
        "mechanical_validation": "absolute",
        "union_scope": "H13 failures strictly less than H12",
        "full_matched_scope": "H13 failures strictly less than H12",
        "mcnemar_and_strata": "descriptive-only",
        "single_cell_regressions": "reported-not-vetoed",
    }
    if (canonical_json(contract["screen"]) != canonical_json(fixed_screen) or
            canonical_json(contract["allk_escalation"]) !=
                canonical_json(fixed_allk) or
            canonical_json(contract["architecture"]) !=
                canonical_json(fixed_architecture) or
            canonical_json(contract["promotion_rule"]) !=
                canonical_json(fixed_promotion) or
            contract["child_lifetime_policy"] !=
                "dedicated process group; registered kill/reap on SIGINT or "
                "SIGTERM; Linux PR_SET_PDEATHSIG=SIGKILL"):
        die("Stage-B campaign geometry/promotion contract mismatch")
    provenance = contract["source_stage_a"]
    provenance_fields = {
        "result_dir", "contract_sha256", "controller_sha256",
        "binary_sha256", "taskset_sha256", "controller_affinity_sha256",
        "terminal_sha256", "oh0_root_control_sha256s",
        "campaign_complete_sha256", "analysis_sha256",
        "source_identity_stream_sha256", "oh0_shard_epoch_sha256s",
        "oh0_shard_epoch_root_sha256", "paired_cells_replayed",
        "union_cells_replayed", "construction_seed_policy",
        "production_seed_fixups_applied", "formal_completion_authenticated",
    }
    epoch = provenance.get("oh0_shard_epoch_sha256s") \
        if isinstance(provenance, dict) else None
    provenance_hashes = (
        "contract_sha256", "controller_sha256", "binary_sha256",
        "taskset_sha256", "controller_affinity_sha256",
        "terminal_sha256",
        "campaign_complete_sha256", "analysis_sha256",
        "source_identity_stream_sha256", "oh0_shard_epoch_root_sha256",
    )
    if (not isinstance(provenance, dict) or set(provenance) != provenance_fields or
            type(provenance.get("result_dir")) is not str or
            provenance.get("controller_sha256") !=
                FORMAL_STAGE_A_CONTROLLER_SHA256 or
            provenance.get("binary_sha256") != FORMAL_STAGE_A_BINARY_SHA256 or
            provenance.get("campaign_complete_sha256") !=
                FORMAL_STAGE_A_CAMPAIGN_COMPLETE_FILE_SHA256 or
            provenance.get("analysis_sha256") !=
                FORMAL_STAGE_A_ANALYSIS_FILE_SHA256 or
            provenance.get("paired_cells_replayed") != SOURCE_PAIRED_CELLS or
            provenance.get("union_cells_replayed") != EXPECTED_SOURCE_UNION or
            provenance.get("construction_seed_policy") !=
                "matrix-c-peel-lo32-xor-hi32-v1" or
            provenance.get("production_seed_fixups_applied") != 0 or
            provenance.get("formal_completion_authenticated") is not True or
            any(type(provenance.get(name)) is not str or
                not re.fullmatch(r"[0-9a-f]{64}", provenance[name])
                for name in provenance_hashes) or
            type(epoch) is not list or len(epoch) != 2 * ALLK_TASKS or
            any(type(value) is not str or
                not re.fullmatch(r"[0-9a-f]{64}", value) for value in epoch) or
            provenance.get("oh0_shard_epoch_root_sha256") != sha256_bytes(
                canonical_json(epoch).encode())):
        die("Stage-B source provenance is malformed")
    root_controls = provenance["oh0_root_control_sha256s"]
    if (type(root_controls) is not list or
            len(root_controls) != len(CONSTRUCTION_ROOTS) or
            any(not isinstance(item, dict) or set(item) != {
                "construction_seed_index", "manifest_sha256",
                "complete_sha256"} or
                item["construction_seed_index"] != index or
                any(type(item[name]) is not str or
                    not re.fullmatch(r"[0-9a-f]{64}", item[name])
                    for name in ("manifest_sha256", "complete_sha256"))
                for index, item in enumerate(root_controls))):
        die("Stage-B source root-control provenance is malformed")
    for name in ("binary", "controller", "taskset"):
        path = Path(contract[name])
        if (path.parent != result_dir / "frozen" or
                path.is_symlink() or not path.is_file() or
                stat.S_IMODE(path.stat().st_mode) != 0o555 or
                sha256_file(path) != contract[f"{name}_sha256"]):
            die(f"frozen Stage-B {name} substitution detected")
    if require_frozen_controller and Path(__file__).resolve() != Path(contract["controller"]):
        die("run/reduce must use the frozen Stage-B controller path")
    telemetry = contract["external_telemetry"]
    if (not isinstance(telemetry, dict) or set(telemetry) != {
            "managed_by_controller", "path", "prepare_mark", "policy",
            "sampler_identity", "continuity"} or
            telemetry["managed_by_controller"] is not False or
            canonical_json(telemetry["policy"]) != canonical_json(TELEMETRY_POLICY) or
            telemetry["continuity"] != TELEMETRY_CONTINUITY or
            telemetry["sampler_identity"] !=
                "external operator responsibility; CPU and all eight DRAM TSOD "
                "temperature channels are sampled by the authenticated log worker"):
        die("Stage-B external telemetry contract mismatch")
    if telemetry["path"] is None:
        if telemetry["prepare_mark"] is not None:
            die("null telemetry path has a prepare mark")
    else:
        if type(telemetry["path"]) is not str:
            die("external telemetry path is malformed")
        validate_external_mark(
            telemetry["prepare_mark"], telemetry["path"], "contract",
        )
    cohort = load_cohort(result_dir, contract["cohort_sha256"])
    if (provenance["source_identity_stream_sha256"] !=
            cohort["source_identity_stream_sha256"]):
        die("source identity stream digest differs between contract and cohort")
    return contract, cohort


class ActiveChildren:
    """Own all benchmark process groups launched by one run invocation."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._stopped = threading.Event()
        self._processes: dict[int, subprocess.Popen[bytes]] = {}
        self._signal: Optional[int] = None

    @staticmethod
    def _kill(process: subprocess.Popen[bytes]) -> None:
        if process.poll() is not None:
            return
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass

    def register(self, process: subprocess.Popen[bytes]) -> None:
        with self._lock:
            self._processes[process.pid] = process
            stopped = self._stopped.is_set()
        if stopped:
            self._kill(process)

    def unregister(self, process: subprocess.Popen[bytes]) -> None:
        with self._lock:
            self._processes.pop(process.pid, None)

    def stop(self, signum: Optional[int] = None) -> None:
        with self._lock:
            if signum is not None and self._signal is None:
                self._signal = signum
            self._stopped.set()
            processes = list(self._processes.values())
        for process in processes:
            self._kill(process)

    def check(self) -> None:
        if self._stopped.is_set():
            if self._signal is not None:
                raise CampaignInterrupted(self._signal)
            die("campaign dispatch cancelled after a job failure")


def _stat_identity(value: os.stat_result) -> tuple[int, ...]:
    return (
        value.st_dev, value.st_ino, value.st_mode, value.st_nlink,
        value.st_size, value.st_mtime_ns, value.st_ctime_ns,
    )


def _stable_fd_sha256(descriptor: int, size: int) -> str:
    digests: list[str] = []
    for _ in range(2):
        digest = hashlib.sha256()
        offset = 0
        while offset < size:
            block = os.pread(descriptor, min(1 << 20, size - offset), offset)
            if not block:
                die("authenticated executable ended before its sealed size")
            digest.update(block)
            offset += len(block)
        if os.pread(descriptor, 1, size):
            die("authenticated executable grew while hashing")
        digests.append(digest.hexdigest())
    if digests[0] != digests[1]:
        die("authenticated executable changed between hash passes")
    return digests[0]


def open_authenticated_executable(path: Path, expected_sha256: str) -> int:
    """Open one immutable frozen executable and bind the hash to its fd."""
    if not re.fullmatch(r"[0-9a-f]{64}", expected_sha256):
        die("authenticated executable has a malformed expected SHA256")
    flags = (os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
             getattr(os, "O_NOFOLLOW", 0))
    try:
        descriptor = os.open(path, flags)
    except OSError as exc:
        die(f"cannot open frozen executable {path}: {exc}")
    try:
        before = os.fstat(descriptor)
        if (not stat.S_ISREG(before.st_mode) or before.st_nlink != 1 or
                stat.S_IMODE(before.st_mode) != 0o555):
            die("frozen executable is not a unique mode-0555 regular file")
        digest = _stable_fd_sha256(descriptor, before.st_size)
        after = os.fstat(descriptor)
        try:
            named = os.stat(path, follow_symlinks=False)
        except OSError as exc:
            die(f"frozen executable pathname changed while opening: {exc}")
        if (_stat_identity(before) != _stat_identity(after) or
                _stat_identity(after) != _stat_identity(named) or
                digest != expected_sha256):
            die("frozen executable identity/hash changed while opening")
        return descriptor
    except BaseException:
        os.close(descriptor)
        raise


# This code is passed directly to the already-running interpreter through
# /proc/self/exe.  It never reopens the controller, taskset, or benchmark by
# pathname.  Hashing uses pread so concurrent children can safely share the
# two authenticated parent descriptors without sharing file offsets.
AUTHENTICATED_FD_LAUNCHER = r'''import ctypes, hashlib, os, signal, stat, sys
def fail(message):
    os.write(2, ("authenticated launcher: " + message + "\n").encode("ascii", "replace"))
    os._exit(125)
def identity(value):
    return (value.st_dev, value.st_ino, value.st_mode, value.st_nlink,
            value.st_size, value.st_mtime_ns, value.st_ctime_ns)
def authenticate(descriptor, expected):
    before = os.fstat(descriptor)
    if (not stat.S_ISREG(before.st_mode) or before.st_nlink != 1 or
            stat.S_IMODE(before.st_mode) != 0o555):
        fail("executable fd is not a unique mode-0555 regular file")
    found = []
    for unused in range(2):
        digest = hashlib.sha256()
        offset = 0
        while offset < before.st_size:
            block = os.pread(descriptor, min(1048576, before.st_size - offset), offset)
            if not block:
                fail("executable fd ended before its sealed size")
            digest.update(block)
            offset += len(block)
        if os.pread(descriptor, 1, before.st_size):
            fail("executable fd grew while hashing")
        found.append(digest.hexdigest())
    after = os.fstat(descriptor)
    if identity(before) != identity(after) or found != [expected, expected]:
        fail("executable fd identity/hash mismatch")
try:
    if len(sys.argv) < 9 or sys.argv[0] != "-c":
        fail("malformed launcher arguments")
    parent = int(sys.argv[1])
    cpu = int(sys.argv[2])
    taskset_fd = int(sys.argv[3])
    taskset_sha256 = sys.argv[4]
    binary_fd = int(sys.argv[5])
    binary_sha256 = sys.argv[6]
    taskset_argv0 = sys.argv[7]
    benchmark_tail = sys.argv[8:]
    if parent <= 1 or cpu < 0 or taskset_fd <= 2 or binary_fd <= 2 or not benchmark_tail:
        fail("launcher scalar argument is outside its domain")
    if os.getppid() != parent:
        os.kill(os.getpid(), signal.SIGKILL)
    libc = ctypes.CDLL(None, use_errno=True)
    if libc.prctl(1, signal.SIGKILL, 0, 0, 0) != 0:
        fail("cannot arm PR_SET_PDEATHSIG")
    if os.getppid() != parent:
        os.kill(os.getpid(), signal.SIGKILL)
    os.sched_setaffinity(0, {cpu})
    if os.sched_getaffinity(0) != {cpu}:
        fail("child affinity did not bind to its exact assigned CPU")
    authenticate(taskset_fd, taskset_sha256)
    authenticate(binary_fd, binary_sha256)
    os.set_inheritable(binary_fd, True)
    command = [taskset_argv0, "-c", str(cpu),
               "/proc/self/fd/" + str(binary_fd)] + benchmark_tail
    os.execve(taskset_fd, command, os.environ)
except BaseException as exc:
    fail("launch failed: " + str(exc))
'''


def parse_cpu_list(text: str) -> list[int]:
    cpus: set[int] = set()
    for token in text.split(","):
        if "-" in token:
            first_text, last_text = token.split("-", 1)
            first = strict_uint(first_text, "CPU")
            last = strict_uint(last_text, "CPU")
            if first > last:
                die("descending CPU range")
            cpus.update(range(first, last + 1))
        else:
            cpus.add(strict_uint(token, "CPU"))
    if not cpus or 127 in cpus:
        die("CPU list is empty or includes reserved CPU 127")
    if not cpus.issubset(os.sched_getaffinity(0)):
        die("requested CPU is outside the controller affinity")
    return sorted(cpus)


def bind_controller_affinity(result_dir: Path, cpus: Sequence[int]) -> None:
    selected = sorted(cpus)
    if not selected or 127 in selected or selected != sorted(set(selected)):
        die("controller affinity is malformed or includes CPU 127")
    try:
        os.sched_setaffinity(0, set(selected))
    except OSError as exc:
        die(f"cannot bind controller affinity: {exc}")
    actual = sorted(os.sched_getaffinity(0))
    if actual != selected:
        die("controller affinity differs from selected CPUs")
    write_sealed_once(
        result_dir / "controller_affinity.json",
        f"{SCHEMA}.controller_affinity",
        {"selected_cpus": selected, "actual_cpus": actual,
         "reserved_cpu_127_excluded": True},
    )


def apply_sealed_controller_affinity(result_dir: Path) -> list[int]:
    payload = load_sealed(
        result_dir / "controller_affinity.json",
        f"{SCHEMA}.controller_affinity",
    )
    if (not isinstance(payload, dict) or set(payload) != {
            "selected_cpus", "actual_cpus", "reserved_cpu_127_excluded"} or
            type(payload["selected_cpus"]) is not list or
            any(type(cpu) is not int for cpu in payload["selected_cpus"]) or
            payload["actual_cpus"] != payload["selected_cpus"] or
            payload["reserved_cpu_127_excluded"] is not True or
            not payload["selected_cpus"] or 127 in payload["selected_cpus"] or
            payload["selected_cpus"] != sorted(set(payload["selected_cpus"]))):
        die("sealed Stage-B controller affinity is malformed")
    bind_controller_affinity(result_dir, payload["selected_cpus"])
    return payload["selected_cpus"]


def resolve_worker_count(requested: Optional[int], cpus: Sequence[int]) -> int:
    workers = len(cpus) if requested is None else requested
    if type(workers) is not int or workers < 1 or workers > len(cpus):
        die("workers must be in 1..len(cpus)")
    return workers


def exclusive_campaign(handler):
    @wraps(handler)
    def locked(args: argparse.Namespace):
        result_dir = args.result_dir.resolve(strict=True)
        lock_path = result_dir / "controller.lock"
        flags = (os.O_RDWR | os.O_CREAT | getattr(os, "O_CLOEXEC", 0) |
                 getattr(os, "O_NOFOLLOW", 0))
        try:
            descriptor = os.open(lock_path, flags, 0o600)
        except OSError as exc:
            die(f"cannot open campaign ownership lock: {exc}")
        try:
            metadata = os.fstat(descriptor)
            if not stat.S_ISREG(metadata.st_mode) or metadata.st_nlink != 1:
                die("campaign ownership lock is not a unique regular file")
            try:
                fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                die("another run/reduce controller owns this campaign")
            os.ftruncate(descriptor, 0)
            os.write(descriptor, f"pid={os.getpid()}\n".encode("ascii"))
            os.fsync(descriptor)
            return handler(args, result_dir)
        finally:
            os.close(descriptor)
    return locked


def screen_stage(result_dir: Path) -> Path:
    return result_dir / "stages" / "screen"


def shard_path(stage: Path, task: Task) -> Path:
    return stage / "shards" / task.stem


def validate_root_inventory(result_dir: Path, contract: Mapping[str, Any]) -> None:
    allowed = {
        "contract.json", "cohort.json", "controller.lock", "frozen", "stages",
        "controller_affinity.json", "telemetry_start.json", "telemetry_end.json",
        "telemetry_null.json", "screen_decision.json", "terminal.json",
        "campaign_complete.json", "analysis.json",
        "allk_telemetry_start.json", "allk_telemetry_end.json",
        "allk_telemetry_null.json", "allk_terminal.json",
        "allk_campaign_complete.json", "allk_analysis.json",
    }
    entries = {entry.name: entry for entry in result_dir.iterdir()}
    if set(entries) - allowed:
        die("Stage-B root contains out-of-ledger entries")
    for name, entry in entries.items():
        if entry.is_symlink():
            die("Stage-B root contains a symlinked entry")
        if name in ("frozen", "stages"):
            if not entry.is_dir():
                die("Stage-B root directory entry has the wrong type")
        elif not entry.is_file():
            die("Stage-B root file entry has the wrong type")
    frozen = result_dir / "frozen"
    expected = {
        Path(contract[name]).name for name in ("binary", "controller", "taskset")
    }
    if (not frozen.is_dir() or frozen.is_symlink() or
            {entry.name for entry in frozen.iterdir()} != expected or
            any(not entry.is_file() or entry.is_symlink()
                for entry in frozen.iterdir())):
        die("frozen Stage-B tool inventory differs from contract")
    stages = result_dir / "stages"
    if not stages.is_dir() or stages.is_symlink():
        die("Stage-B stages root is missing or symlinked")
    stage_entries = list(stages.iterdir())
    stage_names = {entry.name for entry in stage_entries}
    if (not stage_names.issubset({"screen", "allk"}) or
            "screen" not in stage_names or
            any(not entry.is_dir() or entry.is_symlink()
                for entry in stage_entries)):
        die("Stage-B stage inventory is malformed")


def validate_screen_inventory(
    stage: Path, tasks: Sequence[Task], *, complete_required: bool,
) -> None:
    if not stage.is_dir() or stage.is_symlink():
        die("screen stage is missing or symlinked")
    expected = {"manifest.json", "shards"}
    if not (stage / "shards").exists() and not complete_required:
        expected.remove("shards")
    if complete_required or (stage / "complete.json").exists():
        expected.add("complete.json")
    entries = {entry.name: entry for entry in stage.iterdir()}
    if set(entries) != expected:
        die("screen stage file inventory differs from sealed layout")
    manifest = entries["manifest.json"]
    if manifest.is_symlink() or not manifest.is_file():
        die("screen manifest is not a real file")
    if "shards" in entries:
        root = entries["shards"]
        if root.is_symlink() or not root.is_dir():
            die("screen shard root is not a real directory")
        expected_shards = {task.stem for task in tasks}
        actual = list(root.iterdir())
        if complete_required:
            if ({entry.name for entry in actual} != expected_shards or
                    any(entry.is_symlink() or not entry.is_dir()
                        for entry in actual)):
                die("screen shard inventory differs from sealed manifest")
        else:
            task_ids = {task.task_id for task in tasks}
            for entry in actual:
                match = re.fullmatch(
                    r"\.staging-([0-9a-f]{24})-([1-9][0-9]*)", entry.name,
                )
                if (entry.name not in expected_shards and
                        (match is None or match.group(1) not in task_ids)):
                    die("screen shard root contains an out-of-ledger entry")
                if entry.is_symlink() or not entry.is_dir():
                    die("screen shard root contains an unsafe entry")
    if "complete.json" in entries:
        complete = entries["complete.json"]
        if complete.is_symlink() or not complete.is_file():
            die("screen completion is not a real file")


def validate_allk_inventory(
    result_dir: Path, cohort: Mapping[str, Any], *,
    complete_required: Optional[bool],
) -> list[Task]:
    stage = result_dir / "stages" / "allk"
    if not stage.is_dir() or stage.is_symlink():
        die("all-K escalation plan directory is missing or symlinked")
    entries = {entry.name: entry for entry in stage.iterdir()}
    allowed = {"manifest.json", "shards", "complete.json"}
    if set(entries) - allowed or "manifest.json" not in entries:
        die("all-K escalation contains an out-of-ledger entry")
    manifest = entries["manifest.json"]
    if manifest.is_symlink() or not manifest.is_file():
        die("all-K escalation manifest is not a real file")
    if complete_required is False and set(entries) != {"manifest.json"}:
        die("planned all-K escalation has execution artifacts")
    if complete_required is True and set(entries) != allowed:
        die("completed all-K escalation inventory is incomplete")
    tasks = load_manifest(stage, "allk", cohort)
    if (len(tasks) != ALLK_TASKS or
            sum(len(task.ks) for task in tasks) != ALLK_CELLS):
        die("all-K escalation manifest cardinality mismatch")
    if "shards" in entries:
        shards = entries["shards"]
        if shards.is_symlink() or not shards.is_dir():
            die("all-K shard root has unsafe type")
        expected_stems = {task.stem for task in tasks}
        task_ids = {task.task_id for task in tasks}
        shard_entries = list(shards.iterdir())
        for entry in shard_entries:
            staging = re.fullmatch(
                r"\.staging-([0-9a-f]{24})-([1-9][0-9]*)", entry.name,
            )
            if (entry.name not in expected_stems and
                    (staging is None or staging.group(1) not in task_ids)):
                die("all-K shard root contains an out-of-ledger entry")
            if entry.is_symlink() or not entry.is_dir():
                die("all-K shard root contains an unsafe entry")
        if complete_required is True and (
                {entry.name for entry in shard_entries} != expected_stems):
            die("completed all-K shard inventory differs from manifest")
    if "complete.json" in entries:
        complete = entries["complete.json"]
        if complete.is_symlink() or not complete.is_file():
            die("all-K completion is not a real file")
        if "shards" not in entries:
            die("all-K completion exists without its shard ledger")
        shard_entries = list(entries["shards"].iterdir())
        if ({entry.name for entry in shard_entries} !=
                {task.stem for task in tasks} or
                any(entry.is_symlink() or not entry.is_dir()
                    for entry in shard_entries)):
            die("all-K completion exists without every manifest shard")
    return tasks


def receipt_payload(
    task: Task, benchmark_argv: Sequence[str], command: Sequence[str], cpu: int,
    binary_sha256: str, stdout: bytes, stderr: bytes,
) -> dict[str, Any]:
    return {
        "task": task.payload(), "benchmark_argv": list(benchmark_argv),
        "command": list(command), "cpu": cpu,
        "execution_mode": EXECUTION_MODE,
        "binary_sha256": binary_sha256, "returncode": 0,
        "stdout_sha256": sha256_bytes(stdout),
        "stderr_sha256": sha256_bytes(stderr),
        "paired_cell_count": len(task.ks), "row_count": 2 * len(task.ks),
        "loss_seed_formula": LOSS_SEED_FORMULA,
        "loss_seeds": [
            {"K": K, "hex": hex(loss_seed(task.packet_trace_root, K))}
            for K in task.ks
        ],
    }


def validate_shard(
    stage: Path, task: Task, binary: Path, binary_sha256: str, taskset: Path,
    allowed_cpus: set[int], source_cells: Any,
) -> list[dict[str, dict[str, str]]]:
    root = shard_path(stage, task)
    if not root.is_dir() or root.is_symlink():
        die(f"missing sealed shard for {task.task_id}")
    inventory = {entry.name: entry for entry in root.iterdir()}
    if (set(inventory) != {"stdout.csv", "stderr.txt", "receipt.json"} or
            any(entry.is_symlink() or not entry.is_file()
                for entry in inventory.values())):
        die(f"sealed shard inventory mismatch for {task.task_id}")
    stdout = stable_bytes(root / "stdout.csv")
    stderr = stable_bytes(root / "stderr.txt")
    receipt = load_sealed(root / "receipt.json", f"{SCHEMA}.job_receipt")
    expected_argv = make_benchmark_argv(binary, task)
    cpu = receipt.get("cpu") if isinstance(receipt, dict) else None
    if type(cpu) is not int or cpu == 127 or cpu not in allowed_cpus:
        die(f"job receipt CPU is outside sealed affinity for {task.task_id}")
    command = [str(taskset), "-c", str(cpu), *expected_argv]
    expected_receipt = receipt_payload(
        task, expected_argv, command, cpu, binary_sha256, stdout, stderr,
    )
    if stderr or canonical_json(receipt) != canonical_json(expected_receipt):
        die(f"job receipt substitution/content mismatch for {task.task_id}")
    try:
        text = stdout.decode("ascii")
    except UnicodeDecodeError as exc:
        die(f"benchmark output is not ASCII: {exc}")
    return parse_output(text, task, expected_cells_for_task(source_cells, task))


def _run_process(
    command: Sequence[str], timeout: float, children: ActiveChildren, *,
    cpu: int, binary_fd: int, binary_sha256: str,
    taskset_fd: int, taskset_sha256: str,
) -> tuple[int, bytes, bytes]:
    children.check()
    if (sys.platform != "linux" or type(cpu) is not int or cpu < 0 or
            len(command) < 5 or command[1:3] != ["-c", str(cpu)] or
            any(type(value) is not str for value in command) or
            any(type(value) is not int or value <= 2 for value in
                (binary_fd, taskset_fd))):
        die("authenticated benchmark launcher arguments are malformed")
    launcher = [
        "/proc/self/exe", "-I", "-S", "-c", AUTHENTICATED_FD_LAUNCHER,
        str(os.getpid()), str(cpu), str(taskset_fd), taskset_sha256,
        str(binary_fd), binary_sha256, command[0], *command[4:],
    ]
    process = subprocess.Popen(
        launcher, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        start_new_session=True, pass_fds=(taskset_fd, binary_fd),
    )
    children.register(process)
    try:
        try:
            stdout, stderr = process.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            ActiveChildren._kill(process)
            stdout, stderr = process.communicate()
            die(f"benchmark timed out after {timeout:g} seconds")
        except BaseException:
            ActiveChildren._kill(process)
            process.communicate()
            raise
    finally:
        children.unregister(process)
    children.check()
    if len(stdout) > OUTPUT_LIMIT or len(stderr) > OUTPUT_LIMIT:
        die("benchmark output exceeded sealed capture limit")
    return process.returncode, stdout, stderr


def execute_task(
    stage: Path, task: Task, binary: Path, binary_sha256: str,
    taskset: Path, taskset_sha256: str, cpus: "queue.Queue[int]",
    timeout: float, children: ActiveChildren, allowed_cpus: set[int],
    source_cells: Any, binary_fd: int, taskset_fd: int,
) -> None:
    children.check()
    final = shard_path(stage, task)
    if final.exists():
        validate_shard(
            stage, task, binary, binary_sha256, taskset, allowed_cpus,
            source_cells,
        )
        return
    cpu = cpus.get()
    try:
        if cpu == 127:
            die("CPU 127 is reserved for the thermal sampler")
        benchmark_argv = make_benchmark_argv(binary, task)
        command = [str(taskset), "-c", str(cpu), *benchmark_argv]
        returncode, stdout, stderr = _run_process(
            command, timeout, children, cpu=cpu,
            binary_fd=binary_fd, binary_sha256=binary_sha256,
            taskset_fd=taskset_fd, taskset_sha256=taskset_sha256,
        )
        if returncode != 0 or stderr:
            tail = stderr.decode("utf-8", errors="replace")[-1000:]
            die(f"screen job {task.job} failed rc={returncode}: {tail}")
        try:
            parse_output(
                stdout.decode("ascii"), task,
                expected_cells_for_task(source_cells, task),
            )
        except UnicodeDecodeError as exc:
            die(f"screen job {task.job} emitted non-ASCII output: {exc}")
        children.check()
        staging = final.with_name(f".staging-{task.task_id}-{os.getpid()}")
        if staging.exists():
            if staging.is_symlink() or not staging.is_dir():
                die("orphan shard staging path has unsafe type")
            shutil.rmtree(staging)
        staging.mkdir(parents=True)
        atomic_write(staging / "stdout.csv", stdout)
        atomic_write(staging / "stderr.txt", stderr)
        write_sealed_once(
            staging / "receipt.json", f"{SCHEMA}.job_receipt",
            receipt_payload(
                task, benchmark_argv, command, cpu, binary_sha256,
                stdout, stderr,
            ),
        )
        children.check()
        final.parent.mkdir(parents=True, exist_ok=True)
        try:
            staging.rename(final)
        except FileExistsError:
            shutil.rmtree(staging)
            validate_shard(
                stage, task, binary, binary_sha256, taskset, allowed_cpus,
                source_cells,
            )
    finally:
        cpus.put(cpu)


def dispatch_tasks(
    stage: Path, tasks: Sequence[Task], binary: Path, binary_sha256: str,
    taskset: Path, taskset_sha256: str, cpus: "queue.Queue[int]", workers: int,
    timeout: float, allowed_cpus: set[int],
    source_cells: Any,
) -> None:
    children = ActiveChildren()
    previous_handlers: dict[int, Any] = {}
    binary_fd = open_authenticated_executable(binary, binary_sha256)
    try:
        taskset_fd = open_authenticated_executable(taskset, taskset_sha256)
    except BaseException:
        os.close(binary_fd)
        raise

    def stop_handler(signum: int, _frame: Any) -> None:
        children.stop(signum)

    for signum in (signal.SIGINT, signal.SIGTERM):
        previous_handlers[signum] = signal.signal(signum, stop_handler)
    try:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            pending: dict[Any, Task] = {}
            iterator = iter(tasks)
            failure: Optional[BaseException] = None
            completed_count = 0
            while True:
                while failure is None and len(pending) < workers:
                    try:
                        task = next(iterator)
                    except StopIteration:
                        break
                    future = executor.submit(
                        execute_task, stage, task, binary, binary_sha256,
                        taskset, taskset_sha256, cpus, timeout, children,
                        allowed_cpus, source_cells, binary_fd, taskset_fd,
                    )
                    pending[future] = task
                if not pending:
                    break
                done, _ = wait(pending, return_when=FIRST_COMPLETED)
                for future in done:
                    task = pending.pop(future)
                    try:
                        future.result()
                    except BaseException as exc:
                        if failure is None:
                            failure = exc
                            children.stop()
                    else:
                        completed_count += 1
                        if (completed_count % 16 == 0 or
                                completed_count == len(tasks)):
                            print(
                                f"# Stage-B {task.phase} completed "
                                f"{completed_count}/{len(tasks)}",
                                file=sys.stderr, flush=True,
                            )
                if failure is not None:
                    for future in pending:
                        future.cancel()
                    for future in list(pending):
                        try:
                            future.result()
                        except BaseException:
                            pass
                    raise failure
    finally:
        children.stop()
        for signum, handler in previous_handlers.items():
            signal.signal(signum, handler)
        os.close(taskset_fd)
        os.close(binary_fd)


WORK_RECEIPT_FIELDS = (
    "block_xors", "block_muladds", "binary_row_references",
    "binary_row_storage_bytes", "binary_adjacency_storage_bytes",
    "binary_row_storage_allocations", "binary_adjacency_storage_allocations",
    "joint_source_xors", "joint_marginal_xors", "joint_marginal_copies",
    "joint_active_deltas", "joint_scratch_bytes", "dual_source_columns",
    "build_ns", "peel_ns", "project_ns", "residual_ns", "backsub_ns",
    "intermediate_bytes",
)


def _empty_scope() -> dict[str, Any]:
    return {
        "paired_cells": 0,
        "arms": {
            arm: {
                "success": 0, "need_more": 0,
                **{f"{name}_sum": 0 for name in WORK_RECEIPT_FIELDS},
            } for arm in ARMS
        },
        "comparison": {
            "repairs": 0, "introductions": 0, "both_need_more": 0,
            "common_success": 0,
            "common_success_work": {
                arm: {f"{name}_sum": 0 for name in WORK_RECEIPT_FIELDS}
                for arm in ARMS
            },
        },
    }


def _update_scope(
    scope: dict[str, Any], arms: Mapping[str, Mapping[str, str]],
) -> None:
    scope["paired_cells"] += 1
    failures = {arm: arms[arm]["result"] == "1" for arm in ARMS}
    comparison = scope["comparison"]
    comparison["repairs"] += failures["h12"] and not failures["h13"]
    comparison["introductions"] += not failures["h12"] and failures["h13"]
    comparison["both_need_more"] += failures["h12"] and failures["h13"]
    comparison["common_success"] += not failures["h12"] and not failures["h13"]
    for arm in ARMS:
        target = scope["arms"][arm]
        target["need_more"] += failures[arm]
        target["success"] += not failures[arm]
        for name in WORK_RECEIPT_FIELDS:
            value = int(arms[arm][name])
            target[f"{name}_sum"] += value
            if not failures["h12"] and not failures["h13"]:
                comparison["common_success_work"][arm][f"{name}_sum"] += value


def exact_mcnemar(repairs: int, introductions: int) -> str:
    discordant = repairs + introductions
    if discordant == 0:
        return "1"
    lower = min(repairs, introductions)
    numerator = sum(math.comb(discordant, index) for index in range(lower + 1))
    return str(min(
        Decimal(1), Decimal(2 * numerator) / Decimal(1 << discordant),
    ))


def _finish_scope(scope: dict[str, Any]) -> dict[str, Any]:
    comparison = scope["comparison"]
    comparison["exact_two_sided_mcnemar_descriptive"] = exact_mcnemar(
        comparison["repairs"], comparison["introductions"],
    )
    comparison["h13_over_h12_common_success_work_ratios"] = {
        name: (
            None if comparison["common_success_work"]["h12"][f"{name}_sum"] == 0
            else str(
                Decimal(comparison["common_success_work"]["h13"][f"{name}_sum"]) /
                Decimal(comparison["common_success_work"]["h12"][f"{name}_sum"])
            )
        ) for name in WORK_RECEIPT_FIELDS
    }
    return scope


def iter_screen_pairs(
    stage: Path, tasks: Sequence[Task], binary: Path, binary_sha256: str,
    taskset: Path, allowed_cpus: set[int],
    cells: Mapping[CellKey, Mapping[str, Any]],
) -> Iterator[
    tuple[CellKey, dict[str, dict[str, str]], Mapping[str, Any]]
]:
    for task in tasks:
        pairs = validate_shard(
            stage, task, binary, binary_sha256, taskset, allowed_cpus, cells,
        )
        if len(pairs) != len(task.ks):
            die("validated screen shard pair count mismatch")
        for K, arms in zip(task.ks, pairs):
            key: CellKey = (
                task.construction_seed_index, task.packet_trace_root_index,
                task.schedule, K,
            )
            expected = cells.get(key)
            if expected is None:
                die("screen shard contains a cell outside the matched cohort")
            yield key, arms, expected


def aggregate_screen_pairs(
    pairs: Iterable[
        tuple[CellKey, Mapping[str, Mapping[str, str]], Mapping[str, Any]]
    ],
    expected_keys: set[CellKey],
) -> dict[str, Any]:
    scopes = {
        "full_matched": _empty_scope(), "union": _empty_scope(),
        "control": _empty_scope(),
    }
    by_stratum: dict[str, dict[str, Any]] = {}
    seen: set[CellKey] = set()
    weak_cells: list[dict[str, Any]] = []
    for key, arms, expected in pairs:
        if key in seen:
            die("screen aggregate contains a duplicate cell")
        seen.add(key)
        role = expected.get("role")
        match_id = expected.get("match_id")
        if role not in ("union", "control") or type(match_id) is not int:
            die("screen aggregate cell role/match receipt is malformed")
        _update_scope(scopes["full_matched"], arms)
        _update_scope(scopes[role], arms)
        construction_seed_index, packet_trace_root_index, schedule, K = key
        band = band_index(K)
        stratum_name = (
            f"construction{construction_seed_index}:"
            f"trace{packet_trace_root_index}:{schedule}:band{band}"
        )
        stratum = by_stratum.setdefault(stratum_name, {
            "construction_seed_index": construction_seed_index,
            "construction_seed": construction_seed_text(
                construction_seed_index),
            "packet_trace_root_index": packet_trace_root_index,
            "packet_trace_root": LOSS_ROOTS[packet_trace_root_index],
            "schedule": schedule, "band_index": band,
            "band": list(BANDS[band]), "full_matched": _empty_scope(),
            "union": _empty_scope(), "control": _empty_scope(),
        })
        _update_scope(stratum["full_matched"], arms)
        _update_scope(stratum[role], arms)
        if any(arms[arm]["result"] == "1" for arm in ARMS):
            weak_cells.append({
                **cell_key_payload(key), "match_id": match_id, "role": role,
                "outcomes": {
                    arm: arms[arm]["result_name"] for arm in ARMS
                },
                "pair_id": arms["h12"]["pair_id"],
            })
    if seen != expected_keys:
        die("screen aggregate coverage differs from the exact matched cohort")
    for scope in scopes.values():
        _finish_scope(scope)
    for stratum in by_stratum.values():
        for name in ("full_matched", "union", "control"):
            _finish_scope(stratum[name])
    return {
        "scopes": scopes,
        "descriptive_by_stratum_band": {
            name: by_stratum[name] for name in sorted(by_stratum)
        },
        "weak_cells": weak_cells,
    }


def screen_decision_payload(metrics: Mapping[str, Any]) -> dict[str, Any]:
    try:
        union = metrics["scopes"]["union"]
        full = metrics["scopes"]["full_matched"]
        union_h12 = union["arms"]["h12"]["need_more"]
        union_h13 = union["arms"]["h13"]["need_more"]
        full_h12 = full["arms"]["h12"]["need_more"]
        full_h13 = full["arms"]["h13"]["need_more"]
    except (KeyError, TypeError) as exc:
        die(f"screen metrics cannot drive the frozen decision: {exc}")
    if (any(type(value) is not int or value < 0 for value in
            (union_h12, union_h13, full_h12, full_h13)) or
            union.get("paired_cells") != EXPECTED_SOURCE_UNION or
            full.get("paired_cells") != EXPECTED_SCREEN_CELLS):
        die("screen decision scope cardinality/counts are malformed")
    union_improved = union_h13 < union_h12
    full_improved = full_h13 < full_h12
    passed = union_improved and full_improved
    return {
        "mechanical_validation_passed": True,
        "union_scope": {
            "paired_cells": EXPECTED_SOURCE_UNION,
            "h12_need_more": union_h12, "h13_need_more": union_h13,
            "strictly_improved": union_improved,
            "repairs": union["comparison"]["repairs"],
            "introductions": union["comparison"]["introductions"],
        },
        "full_matched_scope": {
            "paired_cells": EXPECTED_SCREEN_CELLS,
            "h12_need_more": full_h12, "h13_need_more": full_h13,
            "strictly_improved": full_improved,
            "repairs": full["comparison"]["repairs"],
            "introductions": full["comparison"]["introductions"],
        },
        "screen_passed": passed,
        "allk_escalation_planned": passed,
        "descriptive_statistics_are_not_gates": True,
        "single_cell_regressions_are_not_a_veto": True,
    }


def compute_screen_payload(
    result_dir: Path, tasks: Sequence[Task], contract: Mapping[str, Any],
    cohort: Mapping[str, Any], allowed_cpus: set[int], *,
    complete_required: bool = True,
) -> dict[str, Any]:
    stage = screen_stage(result_dir)
    validate_screen_inventory(stage, tasks, complete_required=complete_required)
    shards = stage / "shards"
    expected_shards = {task.stem for task in tasks}
    if (not shards.is_dir() or shards.is_symlink() or
            {entry.name for entry in shards.iterdir()} != expected_shards or
            any(entry.is_symlink() or not entry.is_dir()
                for entry in shards.iterdir())):
        die("screen shard inventory is incomplete or differs from manifest")
    cells = cohort_cell_map(cohort)
    metrics = aggregate_screen_pairs(
        iter_screen_pairs(
            stage, tasks, Path(contract["binary"]), contract["binary_sha256"],
            Path(contract["taskset"]), allowed_cpus, cells,
        ),
        set(cells),
    )
    return {
        "phase": "screen", "overhead": 0, "task_count": len(tasks),
        "paired_cell_count": metrics["scopes"]["full_matched"]["paired_cells"],
        "arm_outcome_count": 2 * metrics["scopes"]["full_matched"]["paired_cells"],
        "manifest_sha256": sha256_file(stage / "manifest.json"),
        "cohort_sha256": contract["cohort_sha256"], "metrics": metrics,
    }


def complete_screen(
    result_dir: Path, tasks: Sequence[Task], contract: Mapping[str, Any],
    cohort: Mapping[str, Any], allowed_cpus: set[int],
) -> dict[str, Any]:
    stage = screen_stage(result_dir)
    payload = compute_screen_payload(
        result_dir, tasks, contract, cohort, allowed_cpus,
        complete_required=False,
    )
    write_sealed_once(stage / "complete.json", f"{SCHEMA}.stage_complete", payload)
    stored = load_sealed(stage / "complete.json", f"{SCHEMA}.stage_complete")
    if canonical_json(stored) != canonical_json(payload):
        die("stored screen completion differs from recomputed shards")
    return payload


def verify_complete_screen(
    result_dir: Path, tasks: Sequence[Task], contract: Mapping[str, Any],
    cohort: Mapping[str, Any], allowed_cpus: set[int],
) -> dict[str, Any]:
    stored = load_sealed(
        screen_stage(result_dir) / "complete.json",
        f"{SCHEMA}.stage_complete",
    )
    recomputed = compute_screen_payload(
        result_dir, tasks, contract, cohort, allowed_cpus,
    )
    if canonical_json(stored) != canonical_json(recomputed):
        die("sealed screen completion differs from recomputed shards")
    return stored


def seal_screen_decision(
    result_dir: Path, completed: Mapping[str, Any],
) -> dict[str, Any]:
    decision = {
        **screen_decision_payload(completed["metrics"]),
        "screen_complete_sha256": sha256_file(
            screen_stage(result_dir) / "complete.json"),
    }
    path = result_dir / "screen_decision.json"
    write_sealed_once(path, f"{SCHEMA}.screen_decision", decision)
    stored = load_sealed(path, f"{SCHEMA}.screen_decision")
    if canonical_json(stored) != canonical_json(decision):
        die("stored screen decision differs from mechanical reduction")
    return decision


def verify_screen_decision(
    result_dir: Path, completed: Mapping[str, Any],
) -> dict[str, Any]:
    expected = {
        **screen_decision_payload(completed["metrics"]),
        "screen_complete_sha256": sha256_file(
            screen_stage(result_dir) / "complete.json"),
    }
    stored = load_sealed(
        result_dir / "screen_decision.json", f"{SCHEMA}.screen_decision",
    )
    if canonical_json(stored) != canonical_json(expected):
        die("sealed screen decision differs from recomputed metrics")
    return stored


def plan_allk(result_dir: Path, cohort: Mapping[str, Any]) -> str:
    stage = result_dir / "stages" / "allk"
    tasks = build_allk_tasks()
    write_sealed_once(
        stage / "manifest.json", f"{SCHEMA}.stage_manifest",
        manifest_payload("allk", tasks),
    )
    validate_allk_inventory(result_dir, cohort, complete_required=False)
    return sha256_file(stage / "manifest.json")


def screen_terminal_payload(
    result_dir: Path, decision: Mapping[str, Any],
) -> dict[str, Any]:
    passed = decision.get("screen_passed") is True
    allk = result_dir / "stages" / "allk" / "manifest.json"
    if passed != allk.exists():
        die("all-K plan presence differs from the screen decision")
    return {
        "terminal_phase": "screen",
        "screen_passed": passed,
        "allk_escalation_planned": passed,
        "screen_complete_sha256": sha256_file(
            screen_stage(result_dir) / "complete.json"),
        "screen_decision_sha256": sha256_file(
            result_dir / "screen_decision.json"),
        "allk_manifest_sha256": sha256_file(allk) if passed else None,
    }


def seal_screen_terminal(
    result_dir: Path, decision: Mapping[str, Any],
) -> dict[str, Any]:
    payload = screen_terminal_payload(result_dir, decision)
    write_sealed_once(result_dir / "terminal.json", f"{SCHEMA}.terminal", payload)
    stored = load_sealed(result_dir / "terminal.json", f"{SCHEMA}.terminal")
    if canonical_json(stored) != canonical_json(payload):
        die("stored screen terminal differs from its sealed chain")
    return stored


def verify_screen_terminal(
    result_dir: Path, decision: Mapping[str, Any], cohort: Mapping[str, Any],
) -> dict[str, Any]:
    stored = load_sealed(result_dir / "terminal.json", f"{SCHEMA}.terminal")
    expected = screen_terminal_payload(result_dir, decision)
    if canonical_json(stored) != canonical_json(expected):
        die("screen terminal differs from recomputed chain")
    if stored["screen_passed"]:
        validate_allk_inventory(result_dir, cohort, complete_required=None)
    elif (result_dir / "stages" / "allk").exists():
        die("failed screen has an all-K escalation directory")
    return stored


def seal_screen_campaign_completion(
    result_dir: Path, telemetry_end: Optional[Mapping[str, Any]],
) -> dict[str, Any]:
    null = result_dir / "telemetry_null.json"
    start = result_dir / "telemetry_start.json"
    end = result_dir / "telemetry_end.json"
    if telemetry_end is None:
        if start.exists() or end.exists():
            die("null-telemetry screen has interval receipts")
        write_sealed_once(
            null, f"{SCHEMA}.telemetry_null",
            {"external_telemetry_configured": False},
        )
        telemetry_path = null
    else:
        if null.exists() or not start.exists() or not end.exists():
            die("screen telemetry completion inventory mismatch")
        telemetry_path = end
    payload = {
        "contract_sha256": sha256_file(result_dir / "contract.json"),
        "cohort_sha256": sha256_file(result_dir / "cohort.json"),
        "controller_affinity_sha256": sha256_file(
            result_dir / "controller_affinity.json"),
        "screen_manifest_sha256": sha256_file(
            screen_stage(result_dir) / "manifest.json"),
        "screen_complete_sha256": sha256_file(
            screen_stage(result_dir) / "complete.json"),
        "screen_decision_sha256": sha256_file(
            result_dir / "screen_decision.json"),
        "terminal_sha256": sha256_file(result_dir / "terminal.json"),
        "telemetry_receipt": telemetry_path.name,
        "telemetry_receipt_sha256": sha256_file(telemetry_path),
        "external_telemetry_bound": telemetry_end is not None,
    }
    path = result_dir / "campaign_complete.json"
    write_sealed_once(path, f"{SCHEMA}.campaign_complete", payload)
    stored = load_sealed(path, f"{SCHEMA}.campaign_complete")
    if canonical_json(stored) != canonical_json(payload):
        die("screen campaign completion cross-seal mismatch")
    return stored


def verify_screen_campaign_completion(
    result_dir: Path, contract: Mapping[str, Any],
    telemetry_start_receipt: Optional[Mapping[str, Any]],
    telemetry_end_receipt: Optional[Mapping[str, Any]],
) -> dict[str, Any]:
    configured = contract["external_telemetry"]["path"] is not None
    if (configured != (telemetry_start_receipt is not None) or
            configured != (telemetry_end_receipt is not None)):
        die("screen campaign telemetry disposition mismatch")
    if telemetry_end_receipt is None:
        telemetry_path = result_dir / "telemetry_null.json"
        null = load_sealed(telemetry_path, f"{SCHEMA}.telemetry_null")
        if null != {"external_telemetry_configured": False}:
            die("screen null-telemetry receipt mismatch")
    else:
        telemetry_path = result_dir / "telemetry_end.json"
    expected = {
        "contract_sha256": sha256_file(result_dir / "contract.json"),
        "cohort_sha256": sha256_file(result_dir / "cohort.json"),
        "controller_affinity_sha256": sha256_file(
            result_dir / "controller_affinity.json"),
        "screen_manifest_sha256": sha256_file(
            screen_stage(result_dir) / "manifest.json"),
        "screen_complete_sha256": sha256_file(
            screen_stage(result_dir) / "complete.json"),
        "screen_decision_sha256": sha256_file(
            result_dir / "screen_decision.json"),
        "terminal_sha256": sha256_file(result_dir / "terminal.json"),
        "telemetry_receipt": telemetry_path.name,
        "telemetry_receipt_sha256": sha256_file(telemetry_path),
        "external_telemetry_bound": telemetry_end_receipt is not None,
    }
    stored = load_sealed(
        result_dir / "campaign_complete.json", f"{SCHEMA}.campaign_complete",
    )
    if canonical_json(stored) != canonical_json(expected):
        die("screen campaign completion cross-seal mismatch")
    return stored


def cleanup_orphan_screen_staging(stage: Path, tasks: Sequence[Task]) -> None:
    root = stage / "shards"
    if not root.exists():
        root.mkdir(parents=True)
        return
    if root.is_symlink() or not root.is_dir():
        die("screen shard root has unsafe type")
    task_ids = {task.task_id for task in tasks}
    for entry in root.iterdir():
        match = re.fullmatch(
            r"\.staging-([0-9a-f]{24})-([1-9][0-9]*)", entry.name,
        )
        if match is None:
            continue
        if match.group(1) not in task_ids:
            die("orphan screen staging entry is not bound to a manifest task")
        if entry.is_symlink() or not entry.is_dir():
            die("orphan screen staging entry has unsafe type")
        shutil.rmtree(entry)


@exclusive_campaign
def run_screen(args: argparse.Namespace, result_dir: Path) -> None:
    contract, cohort = load_contract(result_dir, require_frozen_controller=True)
    validate_root_inventory(result_dir, contract)
    tasks = load_manifest(screen_stage(result_dir), "screen", cohort)
    cpus = parse_cpu_list(args.cpus)
    workers = resolve_worker_count(args.workers, cpus)
    bind_controller_affinity(result_dir, cpus)
    allowed_cpus = set(cpus)
    terminal_path = result_dir / "terminal.json"
    if ((result_dir / "campaign_complete.json").exists() and
            not terminal_path.exists()):
        die("screen campaign completion exists without terminal")
    if (result_dir / "telemetry_end.json").exists() and not terminal_path.exists():
        die("screen telemetry end exists before terminal")
    start = telemetry_start(result_dir, contract, "screen")
    if terminal_path.exists():
        completed = verify_complete_screen(
            result_dir, tasks, contract, cohort, allowed_cpus,
        )
        decision = verify_screen_decision(result_dir, completed)
        terminal = verify_screen_terminal(result_dir, decision, cohort)
        end = telemetry_finish(result_dir, contract, start, "screen")
        if (result_dir / "campaign_complete.json").exists():
            completion = verify_screen_campaign_completion(
                result_dir, contract, start, end,
            )
        else:
            completion = seal_screen_campaign_completion(result_dir, end)
        print(canonical_json({
            **terminal, "already_terminal": True,
            "campaign_complete_sha256": sha256_file(
                result_dir / "campaign_complete.json"),
            "campaign_completion": completion,
        }))
        return
    cleanup_orphan_screen_staging(screen_stage(result_dir), tasks)
    validate_screen_inventory(
        screen_stage(result_dir), tasks, complete_required=False,
    )
    cpu_queue: queue.Queue[int] = queue.Queue()
    for cpu in cpus[:workers]:
        cpu_queue.put(cpu)
    cells = cohort_cell_map(cohort)
    dispatch_tasks(
        screen_stage(result_dir), tasks, Path(contract["binary"]),
        contract["binary_sha256"], Path(contract["taskset"]),
        contract["taskset_sha256"], cpu_queue, workers, args.timeout,
        allowed_cpus, cells,
    )
    completed = complete_screen(
        result_dir, tasks, contract, cohort, allowed_cpus,
    )
    decision = seal_screen_decision(result_dir, completed)
    if decision["screen_passed"]:
        plan_allk(result_dir, cohort)
    elif (result_dir / "stages" / "allk").exists():
        die("failed screen cannot retain an all-K escalation directory")
    terminal = seal_screen_terminal(result_dir, decision)
    end = telemetry_finish(result_dir, contract, start, "screen")
    seal_screen_campaign_completion(result_dir, end)
    verify_screen_terminal(result_dir, decision, cohort)
    verify_screen_campaign_completion(result_dir, contract, start, end)
    validate_root_inventory(result_dir, contract)
    print(canonical_json({
        **terminal,
        "campaign_complete_sha256": sha256_file(
            result_dir / "campaign_complete.json"),
        "reduce_command": [
            str(contract["controller"]), "reduce", "--result-dir",
            str(result_dir),
        ],
        "run_allk_command": (
            [str(contract["controller"]), "run-allk", "--result-dir",
             str(result_dir), "--cpus", args.cpus]
            if decision["screen_passed"] else None
        ),
    }))


def cleanup_orphan_allk_staging(stage: Path, tasks: Sequence[Task]) -> None:
    root = stage / "shards"
    if not root.exists():
        root.mkdir(parents=True)
        return
    if root.is_symlink() or not root.is_dir():
        die("all-K shard root has unsafe type")
    task_ids = {task.task_id for task in tasks}
    for entry in root.iterdir():
        match = re.fullmatch(
            r"\.staging-([0-9a-f]{24})-([1-9][0-9]*)", entry.name,
        )
        if match is None:
            continue
        if match.group(1) not in task_ids:
            die("orphan all-K staging entry is not bound to a manifest task")
        if entry.is_symlink() or not entry.is_dir():
            die("orphan all-K staging entry has unsafe type")
        shutil.rmtree(entry)


def iter_allk_pairs(
    stage: Path, tasks: Sequence[Task], binary: Path, binary_sha256: str,
    taskset: Path, allowed_cpus: set[int], oracle: StageAAllKOracle,
) -> Iterator[tuple[CellKey, dict[str, dict[str, str]]]]:
    previous = -1
    for task in tasks:
        pairs = validate_shard(
            stage, task, binary, binary_sha256, taskset, allowed_cpus, oracle,
        )
        if len(pairs) != len(task.ks):
            die("validated all-K shard pair count mismatch")
        for K, arms in zip(task.ks, pairs):
            key: CellKey = (
                task.construction_seed_index, task.packet_trace_root_index,
                task.schedule, K,
            )
            ordinal = cell_ordinal(key)
            if ordinal <= previous:
                die("all-K cells are duplicated or not in canonical order")
            previous = ordinal
            yield key, arms


def aggregate_allk_pairs(
    pairs: Iterable[
        tuple[CellKey, Mapping[str, Mapping[str, str]]]
    ],
) -> dict[str, Any]:
    overall = _empty_scope()
    by_construction_root = {
        str(index): _empty_scope() for index in range(len(CONSTRUCTION_ROOTS))
    }
    by_packet_trace_root = {
        str(index): _empty_scope() for index in range(len(LOSS_ROOTS))
    }
    by_schedule = {schedule: _empty_scope() for schedule in SCHEDULES}
    by_k_band = {
        f"band{index}:{low}-{high}": {
            "band_index": index, "band": [low, high], **_empty_scope(),
        }
        for index, (low, high) in enumerate(BANDS)
    }
    by_root_pair_schedule = {
        f"construction{construction_index}:trace{trace_index}:{schedule}":
            _empty_scope()
        for construction_index in range(len(CONSTRUCTION_ROOTS))
        for trace_index in range(len(LOSS_ROOTS))
        for schedule in SCHEDULES
    }
    failure_by_k: dict[str, dict[int, int]] = {arm: {} for arm in ARMS}
    failure_by_kc: dict[str, dict[tuple[int, int], int]] = {
        arm: {} for arm in ARMS
    }
    failure_records_kc: dict[
        str, dict[tuple[int, int], list[dict[str, Any]]]
    ] = {arm: {} for arm in ARMS}
    weak_cells: list[dict[str, Any]] = []
    count = 0
    previous_ordinal = -1
    for (construction_index, trace_index, schedule, K), arms in pairs:
        if (construction_index not in range(len(CONSTRUCTION_ROOTS)) or
                trace_index not in range(len(LOSS_ROOTS)) or
                schedule not in SCHEDULES or not K_MIN <= K <= K_MAX):
            die("all-K aggregate contains an out-of-domain cell")
        ordinal = cell_ordinal((construction_index, trace_index, schedule, K))
        if ordinal <= previous_ordinal:
            die("all-K aggregate is duplicated or not canonically ordered")
        previous_ordinal = ordinal
        count += 1
        _update_scope(overall, arms)
        _update_scope(by_construction_root[str(construction_index)], arms)
        _update_scope(by_packet_trace_root[str(trace_index)], arms)
        _update_scope(by_schedule[schedule], arms)
        band = band_index(K)
        low, high = BANDS[band]
        _update_scope(by_k_band[f"band{band}:{low}-{high}"], arms)
        _update_scope(by_root_pair_schedule[
            f"construction{construction_index}:trace{trace_index}:{schedule}"
        ], arms)
        failures = {arm: arms[arm]["result"] == "1" for arm in ARMS}
        for arm in ARMS:
            if failures[arm]:
                failure_by_k[arm][K] = failure_by_k[arm].get(K, 0) + 1
                kc = (K, construction_index)
                failure_by_kc[arm][kc] = failure_by_kc[arm].get(kc, 0) + 1
                failure_records_kc[arm].setdefault(kc, []).append({
                    "packet_trace_root_index": trace_index,
                    "packet_trace_root": LOSS_ROOTS[trace_index],
                    "schedule": schedule,
                    "pair_id": arms[arm]["pair_id"],
                })
        if any(failures.values()):
            weak_cells.append({
                **cell_key_payload((construction_index, trace_index, schedule, K)),
                "outcomes": {
                    arm: arms[arm]["result_name"] for arm in ARMS
                },
                "pair_id": arms["h12"]["pair_id"],
            })
    if count != ALLK_CELLS:
        die("all-K reduction did not contain the exact 1,727,973-cell census")
    _finish_scope(overall)
    for mapping in (
            by_construction_root, by_packet_trace_root, by_schedule,
            by_k_band, by_root_pair_schedule):
        for scope in mapping.values():
            _finish_scope(scope)
    weak_k = {}
    for arm in ARMS:
        strata_per_k = checked_product(
            len(CONSTRUCTION_ROOTS), len(LOSS_ROOTS), len(SCHEDULES),
        )
        strata_per_kc = checked_product(len(LOSS_ROOTS), len(SCHEDULES))
        kc_histogram = {
            str(value): 0 for value in range(strata_per_kc + 1)
        }
        for construction_index in range(len(CONSTRUCTION_ROOTS)):
            for K in range(K_MIN, K_MAX + 1):
                failures = failure_by_kc[arm].get(
                    (K, construction_index), 0,
                )
                if failures > strata_per_kc:
                    die("one all-K (K,construction) exceeds nine failures")
                kc_histogram[str(failures)] += 1
        k_histogram = {
            str(value): 0 for value in range(strata_per_k + 1)
        }
        for K in range(K_MIN, K_MAX + 1):
            failures = failure_by_k[arm].get(K, 0)
            if failures > strata_per_k:
                die("one all-K K exceeds the named-stratum failure count")
            k_histogram[str(failures)] += 1
        weak_k[arm] = {
            "sampled_unique_K": K_COUNT,
            "sampled_K_construction_pairs": checked_product(
                K_COUNT, len(CONSTRUCTION_ROOTS),
            ),
            "weak_K_construction_pairs": (
                checked_product(K_COUNT, len(CONSTRUCTION_ROOTS)) -
                kc_histogram["0"]
            ),
            "weak_K": K_COUNT - k_histogram["0"],
            "multi_failure_K_construction_pairs": sum(
                kc_histogram[str(value)]
                for value in range(2, strata_per_kc + 1)
            ),
            "multi_failure_K": sum(
                k_histogram[str(value)]
                for value in range(2, strata_per_k + 1)
            ),
            "maximum_failure_strata_per_K_construction": max(
                failure_by_kc[arm].values(), default=0,
            ),
            "maximum_failure_strata_per_K": max(
                failure_by_k[arm].values(), default=0,
            ),
            # Preserve the original weak-K receipt names while adding the
            # explicit per-(K,C) and per-K vocabulary used by Stage A v2.
            "maximum_failure_strata": max(
                failure_by_k[arm].values(), default=0,
            ),
            f"failure_strata_histogram_per_K_construction_0_to_{strata_per_kc}":
                kc_histogram,
            f"failure_strata_histogram_per_K_0_to_{strata_per_k}": k_histogram,
            f"failure_strata_histogram_0_to_{strata_per_k}": k_histogram,
            "weak_K_construction_records": [
                {
                    "K": K,
                    "construction_seed_index": construction_index,
                    "construction_seed": construction_seed_text(
                        construction_index),
                    "failure_strata": len(records),
                    "failures": records,
                }
                for (K, construction_index), records in sorted(
                    failure_records_kc[arm].items())
            ],
            "weak_K_records": [
                {
                    "K": K,
                    "failure_strata": failure_by_k[arm][K],
                    "construction_seed_multiplicities": [
                        {
                            "construction_seed_index": construction_index,
                            "construction_seed": construction_seed_text(
                                construction_index),
                            "failure_strata": failure_by_kc[arm].get(
                                (K, construction_index), 0),
                        }
                        for construction_index in range(
                            len(CONSTRUCTION_ROOTS))
                    ],
                }
                for K in sorted(failure_by_k[arm])
            ],
        }
    return {
        "overall": overall,
        "descriptive_scopes": {
            "by_construction_root": by_construction_root,
            "by_packet_trace_root": by_packet_trace_root,
            "by_schedule": by_schedule,
            "by_K_band": by_k_band,
            "by_root_pair_schedule": by_root_pair_schedule,
        },
        "weak_K": weak_k, "weak_cells": weak_cells,
    }


def compute_allk_payload(
    result_dir: Path, tasks: Sequence[Task], contract: Mapping[str, Any],
    allowed_cpus: set[int], oracle: StageAAllKOracle, *,
    complete_required: Optional[bool],
) -> dict[str, Any]:
    validate_allk_inventory(
        result_dir, load_cohort(result_dir, contract["cohort_sha256"]),
        complete_required=complete_required,
    )
    stage = result_dir / "stages" / "allk"
    shards = stage / "shards"
    expected = {task.stem for task in tasks}
    if (not shards.is_dir() or shards.is_symlink() or
            {entry.name for entry in shards.iterdir()} != expected or
            any(entry.is_symlink() or not entry.is_dir()
                for entry in shards.iterdir())):
        die("all-K shard inventory is incomplete or differs from manifest")
    metrics = aggregate_allk_pairs(iter_allk_pairs(
        stage, tasks, Path(contract["binary"]), contract["binary_sha256"],
        Path(contract["taskset"]), allowed_cpus, oracle,
    ))
    return {
        "phase": "allk", "overhead": 0, "task_count": len(tasks),
        "paired_cell_count": metrics["overall"]["paired_cells"],
        "arm_outcome_count": 2 * metrics["overall"]["paired_cells"],
        "manifest_sha256": sha256_file(stage / "manifest.json"),
        "screen_decision_sha256": sha256_file(
            result_dir / "screen_decision.json"),
        "metrics": metrics,
    }


def complete_allk(
    result_dir: Path, tasks: Sequence[Task], contract: Mapping[str, Any],
    allowed_cpus: set[int], oracle: StageAAllKOracle,
) -> dict[str, Any]:
    payload = compute_allk_payload(
        result_dir, tasks, contract, allowed_cpus, oracle,
        complete_required=None,
    )
    path = result_dir / "stages" / "allk" / "complete.json"
    write_sealed_once(path, f"{SCHEMA}.allk_complete", payload)
    stored = load_sealed(path, f"{SCHEMA}.allk_complete")
    if canonical_json(stored) != canonical_json(payload):
        die("stored all-K completion differs from recomputed shards")
    validate_allk_inventory(
        result_dir, load_cohort(result_dir, contract["cohort_sha256"]),
        complete_required=True,
    )
    return stored


def verify_complete_allk(
    result_dir: Path, tasks: Sequence[Task], contract: Mapping[str, Any],
    allowed_cpus: set[int], oracle: StageAAllKOracle,
) -> dict[str, Any]:
    stored = load_sealed(
        result_dir / "stages" / "allk" / "complete.json",
        f"{SCHEMA}.allk_complete",
    )
    expected = compute_allk_payload(
        result_dir, tasks, contract, allowed_cpus, oracle,
        complete_required=True,
    )
    if canonical_json(stored) != canonical_json(expected):
        die("sealed all-K completion differs from recomputed shards")
    return stored


def allk_terminal_payload(result_dir: Path) -> dict[str, Any]:
    return {
        "terminal_phase": "allk",
        "screen_decision_sha256": sha256_file(
            result_dir / "screen_decision.json"),
        "screen_campaign_complete_sha256": sha256_file(
            result_dir / "campaign_complete.json"),
        "allk_manifest_sha256": sha256_file(
            result_dir / "stages" / "allk" / "manifest.json"),
        "allk_complete_sha256": sha256_file(
            result_dir / "stages" / "allk" / "complete.json"),
    }


def seal_allk_terminal(result_dir: Path) -> dict[str, Any]:
    payload = allk_terminal_payload(result_dir)
    write_sealed_once(
        result_dir / "allk_terminal.json", f"{SCHEMA}.allk_terminal", payload,
    )
    stored = load_sealed(
        result_dir / "allk_terminal.json", f"{SCHEMA}.allk_terminal",
    )
    if canonical_json(stored) != canonical_json(payload):
        die("stored all-K terminal differs from chain")
    return stored


def verify_allk_terminal(result_dir: Path) -> dict[str, Any]:
    stored = load_sealed(
        result_dir / "allk_terminal.json", f"{SCHEMA}.allk_terminal",
    )
    expected = allk_terminal_payload(result_dir)
    if canonical_json(stored) != canonical_json(expected):
        die("all-K terminal differs from recomputed chain")
    return stored


def seal_allk_campaign_completion(
    result_dir: Path, telemetry_end: Optional[Mapping[str, Any]],
) -> dict[str, Any]:
    null = result_dir / "allk_telemetry_null.json"
    start = result_dir / "allk_telemetry_start.json"
    end = result_dir / "allk_telemetry_end.json"
    if telemetry_end is None:
        if start.exists() or end.exists():
            die("null-telemetry all-K run has interval receipts")
        write_sealed_once(
            null, f"{SCHEMA}.allk_telemetry_null",
            {"external_telemetry_configured": False},
        )
        telemetry_path = null
    else:
        if null.exists() or not start.exists() or not end.exists():
            die("all-K telemetry completion inventory mismatch")
        telemetry_path = end
    payload = {
        "contract_sha256": sha256_file(result_dir / "contract.json"),
        "screen_campaign_complete_sha256": sha256_file(
            result_dir / "campaign_complete.json"),
        "controller_affinity_sha256": sha256_file(
            result_dir / "controller_affinity.json"),
        "allk_manifest_sha256": sha256_file(
            result_dir / "stages" / "allk" / "manifest.json"),
        "allk_complete_sha256": sha256_file(
            result_dir / "stages" / "allk" / "complete.json"),
        "allk_terminal_sha256": sha256_file(
            result_dir / "allk_terminal.json"),
        "telemetry_receipt": telemetry_path.name,
        "telemetry_receipt_sha256": sha256_file(telemetry_path),
        "external_telemetry_bound": telemetry_end is not None,
    }
    path = result_dir / "allk_campaign_complete.json"
    write_sealed_once(path, f"{SCHEMA}.allk_campaign_complete", payload)
    stored = load_sealed(path, f"{SCHEMA}.allk_campaign_complete")
    if canonical_json(stored) != canonical_json(payload):
        die("all-K campaign completion cross-seal mismatch")
    return stored


def verify_allk_campaign_completion(
    result_dir: Path, contract: Mapping[str, Any],
    telemetry_start_receipt: Optional[Mapping[str, Any]],
    telemetry_end_receipt: Optional[Mapping[str, Any]],
) -> dict[str, Any]:
    configured = contract["external_telemetry"]["path"] is not None
    if (configured != (telemetry_start_receipt is not None) or
            configured != (telemetry_end_receipt is not None)):
        die("all-K telemetry disposition mismatch")
    if telemetry_end_receipt is None:
        telemetry_path = result_dir / "allk_telemetry_null.json"
        null = load_sealed(
            telemetry_path, f"{SCHEMA}.allk_telemetry_null",
        )
        if null != {"external_telemetry_configured": False}:
            die("all-K null-telemetry receipt mismatch")
    else:
        telemetry_path = result_dir / "allk_telemetry_end.json"
    expected = {
        "contract_sha256": sha256_file(result_dir / "contract.json"),
        "screen_campaign_complete_sha256": sha256_file(
            result_dir / "campaign_complete.json"),
        "controller_affinity_sha256": sha256_file(
            result_dir / "controller_affinity.json"),
        "allk_manifest_sha256": sha256_file(
            result_dir / "stages" / "allk" / "manifest.json"),
        "allk_complete_sha256": sha256_file(
            result_dir / "stages" / "allk" / "complete.json"),
        "allk_terminal_sha256": sha256_file(
            result_dir / "allk_terminal.json"),
        "telemetry_receipt": telemetry_path.name,
        "telemetry_receipt_sha256": sha256_file(telemetry_path),
        "external_telemetry_bound": telemetry_end_receipt is not None,
    }
    stored = load_sealed(
        result_dir / "allk_campaign_complete.json",
        f"{SCHEMA}.allk_campaign_complete",
    )
    if canonical_json(stored) != canonical_json(expected):
        die("all-K campaign completion cross-seal mismatch")
    return stored


def verify_passing_screen_guard(
    result_dir: Path, contract: Mapping[str, Any], cohort: Mapping[str, Any],
    allowed_cpus: set[int],
) -> tuple[dict[str, Any], dict[str, Any]]:
    screen_tasks = load_manifest(screen_stage(result_dir), "screen", cohort)
    completed = verify_complete_screen(
        result_dir, screen_tasks, contract, cohort, allowed_cpus,
    )
    decision = verify_screen_decision(result_dir, completed)
    terminal = verify_screen_terminal(result_dir, decision, cohort)
    if (decision.get("screen_passed") is not True or
            terminal.get("allk_escalation_planned") is not True):
        die("all-K execution requires an authenticated passing screen")
    screen_start = telemetry_start(result_dir, contract, "screen")
    screen_end = telemetry_finish(
        result_dir, contract, screen_start, "screen",
    )
    verify_screen_campaign_completion(
        result_dir, contract, screen_start, screen_end,
    )
    return completed, decision


@exclusive_campaign
def run_allk(args: argparse.Namespace, result_dir: Path) -> None:
    contract, cohort = load_contract(result_dir, require_frozen_controller=True)
    validate_root_inventory(result_dir, contract)
    cpus = parse_cpu_list(args.cpus)
    workers = resolve_worker_count(args.workers, cpus)
    bind_controller_affinity(result_dir, cpus)
    allowed_cpus = set(cpus)
    verify_passing_screen_guard(result_dir, contract, cohort, allowed_cpus)
    fresh_cohort, fresh_provenance = authenticate_stage_a(
        Path(contract["source_stage_a"]["result_dir"]),
    )
    if (canonical_json(fresh_cohort) != canonical_json(cohort) or
            canonical_json(fresh_provenance) !=
            canonical_json(contract["source_stage_a"])):
        die("all-K oracle replay differs from the prepared Stage-A source")
    oracle = StageAAllKOracle(
        Path(contract["source_stage_a"]["result_dir"]),
        contract["source_stage_a"],
    )
    tasks = validate_allk_inventory(
        result_dir, cohort, complete_required=None,
    )
    terminal_path = result_dir / "allk_terminal.json"
    if ((result_dir / "allk_campaign_complete.json").exists() and
            not terminal_path.exists()):
        die("all-K campaign completion exists without terminal")
    if ((result_dir / "allk_telemetry_end.json").exists() and
            not terminal_path.exists()):
        die("all-K telemetry end exists before terminal")
    if (terminal_path.exists() and
            contract["external_telemetry"]["path"] is not None and
            not (result_dir / "allk_telemetry_start.json").exists()):
        die("all-K terminal lacks its telemetry start receipt")
    start = telemetry_start(result_dir, contract, "allk")
    if terminal_path.exists():
        completed = verify_complete_allk(
            result_dir, tasks, contract, allowed_cpus, oracle,
        )
        terminal = verify_allk_terminal(result_dir)
        end = telemetry_finish(result_dir, contract, start, "allk")
        if (result_dir / "allk_campaign_complete.json").exists():
            verify_allk_campaign_completion(result_dir, contract, start, end)
        else:
            seal_allk_campaign_completion(result_dir, end)
        print(canonical_json({
            **terminal, "already_terminal": True,
            "paired_cells": completed["paired_cell_count"],
            "allk_campaign_complete_sha256": sha256_file(
                result_dir / "allk_campaign_complete.json"),
        }))
        return
    cleanup_orphan_allk_staging(result_dir / "stages" / "allk", tasks)
    validate_allk_inventory(result_dir, cohort, complete_required=None)
    cpu_queue: queue.Queue[int] = queue.Queue()
    for cpu in cpus[:workers]:
        cpu_queue.put(cpu)
    dispatch_tasks(
        result_dir / "stages" / "allk", tasks, Path(contract["binary"]),
        contract["binary_sha256"], Path(contract["taskset"]),
        contract["taskset_sha256"], cpu_queue, workers, args.timeout,
        allowed_cpus, oracle,
    )
    completed = complete_allk(
        result_dir, tasks, contract, allowed_cpus, oracle,
    )
    terminal = seal_allk_terminal(result_dir)
    end = telemetry_finish(result_dir, contract, start, "allk")
    seal_allk_campaign_completion(result_dir, end)
    verify_complete_allk(result_dir, tasks, contract, allowed_cpus, oracle)
    verify_allk_terminal(result_dir)
    verify_allk_campaign_completion(result_dir, contract, start, end)
    validate_root_inventory(result_dir, contract)
    print(canonical_json({
        **terminal, "paired_cells": completed["paired_cell_count"],
        "allk_campaign_complete_sha256": sha256_file(
            result_dir / "allk_campaign_complete.json"),
        "reduce_allk_command": [
            str(contract["controller"]), "reduce-allk", "--result-dir",
            str(result_dir),
        ],
    }))


@exclusive_campaign
def reduce_screen(args: argparse.Namespace, result_dir: Path) -> None:
    contract, cohort = load_contract(result_dir, require_frozen_controller=True)
    allowed_cpus = set(apply_sealed_controller_affinity(result_dir))
    validate_root_inventory(result_dir, contract)
    tasks = load_manifest(screen_stage(result_dir), "screen", cohort)
    completed = verify_complete_screen(
        result_dir, tasks, contract, cohort, allowed_cpus,
    )
    decision = verify_screen_decision(result_dir, completed)
    terminal = verify_screen_terminal(result_dir, decision, cohort)
    start = telemetry_start(result_dir, contract, "screen")
    end = telemetry_finish(result_dir, contract, start, "screen")
    completion = verify_screen_campaign_completion(
        result_dir, contract, start, end,
    )
    summary = {
        "schema": f"{SCHEMA}.analysis",
        "phase": "screen", "coverage": {
            "source_paired_cells_replayed": SOURCE_PAIRED_CELLS,
            "source_union_cells": EXPECTED_SOURCE_UNION,
            "matched_controls": EXPECTED_SCREEN_MATCHES,
            "screen_paired_cells": EXPECTED_SCREEN_CELLS,
            "screen_arm_outcomes": EXPECTED_SCREEN_OUTCOMES,
            "nonempty_strata": EXPECTED_NONEMPTY_STRATA,
            "zero_strata": [list(value) for value in EXPECTED_ZERO_STRATA],
            "screen_pair_native_tasks": EXPECTED_SCREEN_TASKS,
        },
        "source_stage_a": contract["source_stage_a"],
        "cohort_sha256": contract["cohort_sha256"],
        "screen": completed["metrics"], "decision": decision,
        "terminal": terminal,
        "allk_plan": ({
            "tasks": ALLK_TASKS, "paired_cells": ALLK_CELLS,
            "arm_outcomes": ALLK_OUTCOMES,
            "manifest_sha256": terminal["allk_manifest_sha256"],
            "execution_is_explicit": True,
        } if decision["screen_passed"] else None),
        "external_telemetry": {
            "policy": TELEMETRY_POLICY, "start": start, "end": end,
        },
        "campaign_completion": completion,
        "interpretation": (
            "The deterministic matched screen promotes only when H13 has "
            "strictly fewer failures on both the raw-union and full matched "
            "scopes. McNemar, per-stratum, work, timing, and individual "
            "regression records are descriptive and are not vetoes. Passing "
            "seals but does not implicitly execute the exact all-K plan."
        ),
    }
    write_sealed_once(
        result_dir / "analysis.json", f"{SCHEMA}.analysis_record", summary,
    )
    stored = load_sealed(
        result_dir / "analysis.json", f"{SCHEMA}.analysis_record",
    )
    if canonical_json(stored) != canonical_json(summary):
        die("stored screen analysis differs from verified reduction")
    validate_root_inventory(result_dir, contract)
    print(canonical_json(summary))


@exclusive_campaign
def reduce_allk(args: argparse.Namespace, result_dir: Path) -> None:
    contract, cohort = load_contract(result_dir, require_frozen_controller=True)
    allowed_cpus = set(apply_sealed_controller_affinity(result_dir))
    validate_root_inventory(result_dir, contract)
    verify_passing_screen_guard(result_dir, contract, cohort, allowed_cpus)
    fresh_cohort, fresh_provenance = authenticate_stage_a(
        Path(contract["source_stage_a"]["result_dir"]),
    )
    if (canonical_json(fresh_cohort) != canonical_json(cohort) or
            canonical_json(fresh_provenance) !=
            canonical_json(contract["source_stage_a"])):
        die("all-K reduction oracle differs from prepared Stage-A source")
    oracle = StageAAllKOracle(
        Path(contract["source_stage_a"]["result_dir"]),
        contract["source_stage_a"],
    )
    tasks = validate_allk_inventory(
        result_dir, cohort, complete_required=True,
    )
    completed = verify_complete_allk(
        result_dir, tasks, contract, allowed_cpus, oracle,
    )
    terminal = verify_allk_terminal(result_dir)
    start = telemetry_start(result_dir, contract, "allk")
    end = telemetry_finish(result_dir, contract, start, "allk")
    completion = verify_allk_campaign_completion(
        result_dir, contract, start, end,
    )
    overall = completed["metrics"]["overall"]
    summary = {
        "schema": f"{SCHEMA}.allk_analysis",
        "phase": "allk", "coverage": {
            "K_min": K_MIN, "K_max": K_MAX, "unique_K": K_COUNT,
            "construction_roots": len(CONSTRUCTION_ROOTS),
            "packet_trace_roots": len(LOSS_ROOTS),
            "schedules": len(SCHEDULES),
            "paired_cells": ALLK_CELLS, "arm_outcomes": ALLK_OUTCOMES,
            "pair_native_tasks": ALLK_TASKS,
        },
        "authenticated_stage_a_identity_stream_sha256":
            fresh_provenance["source_identity_stream_sha256"],
        "metrics": completed["metrics"],
        "comparison": {
            "h12_need_more": overall["arms"]["h12"]["need_more"],
            "h13_need_more": overall["arms"]["h13"]["need_more"],
            "h13_strictly_improved":
                overall["arms"]["h13"]["need_more"] <
                overall["arms"]["h12"]["need_more"],
            "repairs": overall["comparison"]["repairs"],
            "introductions": overall["comparison"]["introductions"],
            "exact_two_sided_mcnemar_descriptive":
                overall["comparison"]["exact_two_sided_mcnemar_descriptive"],
        },
        "terminal": terminal,
        "external_telemetry": {
            "policy": TELEMETRY_POLICY, "start": start, "end": end,
        },
        "campaign_completion": completion,
        "interpretation": (
            "Exact deterministic all-K census over all 27 named "
            "construction-root/loss-root/schedule strata, with every receipt "
            "checked against the pinned Stage-A OH0 oracle. Weak (K,C), "
            "weak-K, and per-K-band summaries describe the grouped Stage-B "
            "architecture itself. Statistical values are descriptive rather "
            "than IID population claims."
        ),
    }
    write_sealed_once(
        result_dir / "allk_analysis.json",
        f"{SCHEMA}.allk_analysis_record", summary,
    )
    stored = load_sealed(
        result_dir / "allk_analysis.json",
        f"{SCHEMA}.allk_analysis_record",
    )
    if canonical_json(stored) != canonical_json(summary):
        die("stored all-K analysis differs from verified reduction")
    validate_root_inventory(result_dir, contract)
    print(canonical_json(summary))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prepare_parser = commands.add_parser("prepare")
    prepare_parser.add_argument("--binary", type=Path, required=True)
    prepare_parser.add_argument("--stage-a-result", type=Path, required=True)
    prepare_parser.add_argument("--result-dir", type=Path, required=True)
    prepare_parser.add_argument(
        "--telemetry-log", type=Path,
        help="existing external CPU/DRAM-temperature telemetry CSV",
    )
    prepare_parser.set_defaults(handler=prepare)
    for name, handler in (("run", run_screen), ("run-allk", run_allk)):
        command = commands.add_parser(name)
        command.add_argument("--result-dir", type=Path, required=True)
        command.add_argument("--cpus", required=True)
        command.add_argument("--workers", type=int)
        command.add_argument("--timeout", type=float, default=300.0)
        command.set_defaults(handler=handler)
    reduce_parser = commands.add_parser("reduce")
    reduce_parser.add_argument("--result-dir", type=Path, required=True)
    reduce_parser.set_defaults(handler=reduce_screen)
    reduce_allk_parser = commands.add_parser("reduce-allk")
    reduce_allk_parser.add_argument("--result-dir", type=Path, required=True)
    reduce_allk_parser.set_defaults(handler=reduce_allk)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    try:
        args = build_parser().parse_args(argv)
        timeout = getattr(args, "timeout", 1.0)
        if not math.isfinite(timeout) or timeout <= 0:
            die("timeout must be finite and positive")
        args.handler(args)
        return 0
    except CampaignInterrupted as exc:
        print(f"wh2_h13_stage_b: {exc}", file=sys.stderr)
        return 128 + exc.signum
    except CampaignError as exc:
        print(f"wh2_h13_stage_b: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
