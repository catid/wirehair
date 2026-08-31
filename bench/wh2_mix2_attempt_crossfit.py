#!/usr/bin/env python3
"""Cross-fit MIX2 construction-attempt quality on already-consumed roots.

This is a bounded exploratory selector test.  It cannot promote a profile and
it never names a fresh short-screen or final-validation root.  Production
routing and state are not modified: exact uint8 attempts exercise the existing
production construction/solve code through the ``precodefail`` test hook, and
the resulting matrix is cross-fitted three ways.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

try:
    from bench import wh2_mix2_promotion_short_screen as screen
except ModuleNotFoundError as exc:
    if exc.name != "bench":
        raise
    import wh2_mix2_promotion_short_screen as screen


CONTRACT_SCHEMA = "wirehair.wh2.mix2-attempt-crossfit-contract.v1"
INPUT_SCHEMA = "wirehair.wh2.mix2-attempt-crossfit-input.v1"
RECORD_SCHEMA = "wirehair.wh2.mix2-attempt-crossfit-record.v1"
REPLAY_SCHEMA = "wirehair.wh2.mix2-attempt-crossfit-replay.v1"
SUMMARY_SCHEMA = "wirehair.wh2.mix2-attempt-crossfit-summary.v1"
BENCH_DESCRIPTION_SCHEMA = "wirehair.wh2.v2-bench-description.v1"
EXPECTED_CANDIDATE_PROFILE_SHA256 = (
    "90233b44a0893f96c1a18c19aa61ada052c935a48c6bf7d6a2813065856651f0")

K_VALUES = (
    2, 3, 4, 5, 6, 8, 16, 32, 64, 100, 101, 128, 256, 512, 513,
    1000, 1001, 2048, 4096, 5000, 5001, 8192, 10000, 10001, 16384,
    20000, 20001, 32768, 49152, 64000,
)
ROOTS = (
    "0xd1b54a32d192ed03",
    "0x94d049bb133111eb",
    "0x8538ecb5bd456ea3",
    "0xc0ac29b7c97c50dd",
    "0x3f84d5b5b5470917",
    "0x9216d5d98979fb1b",
    "0xb889883a79549774",
    "0xb5666de0987896af",
    "0x8bfca269b0bc01e0",
    "0xc4695292d9835286",
    "0x7ccd510f122fc160",
    "0x7001a960b7d9c0a4",
)
SCHEDULES = ("burst", "adversarial", "repair-only")
FOLDS = (
    (0, 1, 6, 9),
    (2, 3, 7, 10),
    (4, 5, 8, 11),
)
ATTEMPTS = tuple(range(256))
ANCHOR_ROOT_INDEX = 0
ANCHOR_SCHEDULE_INDEX = 0
BLOCK_BYTES = 2
LOSS_PPM = 500000
OVERHEAD = 0
MATRIX_CELL_COUNT = len(K_VALUES) * len(ATTEMPTS) * len(ROOTS) * len(SCHEDULES)
OOF_CELL_COUNT = len(FOLDS) * len(K_VALUES) * len(FOLDS[0]) * len(SCHEDULES)
ANCHOR_PROCESS_COUNT = len(K_VALUES) * len(ATTEMPTS)
MAX_BATCH_PROCESS_COUNT = (
    (len(ROOTS) * len(SCHEDULES) - 1) * len(ATTEMPTS))
MAX_MATRIX_PROCESS_COUNT = ANCHOR_PROCESS_COUNT + MAX_BATCH_PROCESS_COUNT
MAX_JOBS = 32
MAX_STDOUT_BYTES = 512 * 1024
MAX_STDERR_BYTES = 64 * 1024
PROCESS_TIMEOUT_SECONDS = 180
UINT64_MAX = (1 << 64) - 1
CHILD_ENVIRONMENT = {"LANG": "C", "LC_ALL": "C", "TZ": "UTC"}
ZERO_SHA256 = "0" * 64
SHA256 = re.compile(r"[0-9a-f]{64}\Z")
HEX64 = re.compile(r"0x[0-9a-f]{16}\Z")
HEX32 = re.compile(r"0x[0-9a-f]{8}\Z")
UINT = re.compile(r"(?:0|[1-9][0-9]*)\Z")
THREE_DECIMAL = re.compile(r"(?:0|[1-9][0-9]*)\.[0-9]{3}\Z")
EIGHT_DECIMAL = re.compile(r"(?:0|[1-9][0-9]*)\.[0-9]{8}\Z")

STAIRCASE_BY_K = {
    2: 2, 3: 3, 4: 3, 5: 5, 6: 6, 8: 6, 16: 13, 32: 13,
    64: 19, 100: 26, 101: 26, 128: 30, 256: 30, 512: 38,
    513: 38, 1000: 50, 1001: 50, 2048: 54, 4096: 70, 5000: 66,
    5001: 66, 8192: 78, 10000: 86, 10001: 86, 16384: 114,
    20000: 134, 20001: 134, 32768: 194, 49152: 378, 64000: 346,
}

CONTROLLER_PATH = Path(__file__).resolve()
REPO_ROOT = CONTROLLER_PATH.parent.parent

# Updated only when the preregistered contract object intentionally changes.
# Keeping the digest literal prevents a silent roster/gate edit from blessing
# itself with a newly computed hash.
EXPECTED_CONTRACT_SHA256 = (
    "cb15c941b3f303a544945934c6a31181c2109f3ae61f348153515d037dae119c")


class CrossfitError(RuntimeError):
    """The frozen exploratory matrix cannot be executed or admitted safely."""


def fail(message: str) -> None:
    raise CrossfitError(message)


def canonical_json(value: Any) -> str:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False)


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_json(value: Any) -> str:
    return sha256_bytes(canonical_json(value).encode("ascii"))


def self_hash(value: Mapping[str, Any], field: str) -> str:
    body = dict(value)
    body.pop(field, None)
    return sha256_json(body)


def _contract_body() -> Mapping[str, Any]:
    return {
        "schema": CONTRACT_SCHEMA,
        "preregistration": "wirehair-sxvz.16.1.20.38@2026-08-31T00:43Z",
        "scope": "exploratory attempt-quality cross-fit only",
        "promotion_evidence": False,
        "fresh_roots_used": False,
        "candidate": {
            "profile_sha256": EXPECTED_CANDIDATE_PROFILE_SHA256,
            "field": "GF(256)",
            "dense_anchor_layout": "two07",
            "mix_count": 2,
            "binary_dense_rows": 12,
            "gf256_heavy_rows": 12,
            "construction_seed_basis": "production-profile",
            "precode_attempt_stride":
                "0x{:016x}".format(screen.PRECODE_ATTEMPT_STRIDE),
            "packet_attempt_stride":
                "0x{:08x}".format(screen.PACKET_ATTEMPT_STRIDE),
            "profile_equivalence":
                "at B=2, the production-profile CLI basis is exactly the "
                "production-normalized-b2-v1 candidate profile",
        },
        "domain": {
            "K": list(K_VALUES),
            "attempt_range_inclusive": [0, 255],
            "attempt_count": len(ATTEMPTS),
            "roots": list(ROOTS),
            "schedules": list(SCHEDULES),
            "block_bytes": BLOCK_BYTES,
            "loss_ppm": LOSS_PPM,
            "overhead": OVERHEAD,
            "solve_semantics":
                "deterministic two-byte-wide all-zero-RHS rank/work solve; "
                "the separate payload-e2e byte reconstruction is not run",
            "staircase_by_K": {
                str(K): STAIRCASE_BY_K[K] for K in K_VALUES},
            "logical_matrix_cells": MATRIX_CELL_COUNT,
        },
        "folds": [list(fold) for fold in FOLDS],
        "fold_rule":
            "each fold contains two original roots, one v8 short root, and "
            "one v9 short root",
        "selector": {
            "priority_key": [
                "training_K_with_any_weak_cell",
                "training_weak_cells",
                "uint8_attempt",
            ],
            "per_K_rule":
                "first priority attempt passing all 24 training cells",
            "baseline": "numeric ascending first attempt passing training",
            "weak_cell":
                "any non-success, including exact construction failure",
            "heldout_cells_per_fold_K": 12,
            "out_of_fold_cells": OOF_CELL_COUNT,
        },
        "transport": {
            "anchor_root_index": ANCHOR_ROOT_INDEX,
            "anchor_schedule_index": ANCHOR_SCHEDULE_INDEX,
            "single_K_anchor_processes": ANCHOR_PROCESS_COUNT,
            "maximum_remaining_batch_processes": MAX_BATCH_PROCESS_COUNT,
            "maximum_matrix_processes": MAX_MATRIX_PROCESS_COUNT,
            "child_environment": dict(CHILD_ENVIRONMENT),
            "per_process_timeout_seconds": PROCESS_TIMEOUT_SECONDS,
            "benchmark_description_schema": BENCH_DESCRIPTION_SCHEMA,
            "configuration_failure":
                "accept only exact exit-2 stderr at the anchor and cache the "
                "root-independent failure across all 36 cells",
            "selected_out_of_fold_replay": {
                "if_all_priority_fold_K_selectable": OOF_CELL_COUNT,
                "exact_fields": [
                    "terminal_outcome", "rank_counters",
                    "inactivated_columns", "effective_seeds",
                    "block_xors", "gf256_muladds",
                ],
                "otherwise":
                    "publish zero replay rows and a scientific REJECT "
                    "because no complete selected set exists to replay",
            },
        },
        "gate": {
            "every_priority_fold_K_selectable": True,
            "every_ascending_fold_K_selectable": True,
            "priority_selector_weak_cells_max": 0,
            "strictly_fewer_weak_cells_than_ascending": True,
            "selected_replay_exact_when_priority_selectable": True,
        },
        "stop_rule":
            "a REJECT retires seed-only MIX2 repair; a PASS can license only "
            "one separately preregistered terminal v10 with fresh roots",
    }


def contract_description() -> Mapping[str, Any]:
    body = dict(_contract_body())
    digest = sha256_json(body)
    if digest != EXPECTED_CONTRACT_SHA256:
        fail("attempt-crossfit contract body differs from its frozen digest")
    body["contract_sha256"] = digest
    return body


def _validate_constants() -> None:
    if (len(K_VALUES) != 30 or tuple(sorted(set(K_VALUES))) != K_VALUES or
            K_VALUES[0] != 2 or K_VALUES[-1] != 64000):
        fail("K roster is not the frozen 30-value boundary set")
    if len(ROOTS) != 12 or len(set(ROOTS)) != len(ROOTS) or \
            any(not re.fullmatch(r"0x[0-9a-f]{16}", root) or
                int(root, 16) == 0 for root in ROOTS):
        fail("consumed-root roster is malformed")
    if ROOTS != tuple(screen.SELECTION_ROOTS) + tuple(screen.ROOTS):
        fail("cross-fit roots are not exactly the consumed v9 root roster")
    if (K_VALUES != tuple(screen.K_VALUES) or
            SCHEDULES != tuple(screen.SCHEDULES) or
            BLOCK_BYTES != screen.BLOCK_BYTES or
            LOSS_PPM != screen.LOSS_PPM or OVERHEAD != screen.OVERHEAD):
        fail("cross-fit coordinates differ from the v9 screen")
    if screen.candidate_profile_sha256() != \
            EXPECTED_CANDIDATE_PROFILE_SHA256:
        fail("cross-fit candidate differs from the frozen v9 profile")
    if screen.BENCH_DESCRIPTION_SCHEMA != BENCH_DESCRIPTION_SCHEMA:
        fail("benchmark description schema differs from the frozen preflight")
    if (set(ROOTS) & set(screen.FINAL_VALIDATION_ROOTS) or
            set(STAIRCASE_BY_K) != set(K_VALUES) or
            tuple(ATTEMPTS) != tuple(range(256))):
        fail("cross-fit domain includes an unfrozen coordinate")
    flattened = tuple(index for fold in FOLDS for index in fold)
    if tuple(sorted(flattened)) != tuple(range(len(ROOTS))):
        fail("cross-fit folds are not an exact root partition")
    for fold in FOLDS:
        if (len(fold) != 4 or
                sum(index < 6 for index in fold) != 2 or
                sum(6 <= index < 9 for index in fold) != 1 or
                sum(index >= 9 for index in fold) != 1):
            fail("cross-fit fold is not source-balanced")
    if (MATRIX_CELL_COUNT != 276480 or OOF_CELL_COUNT != 1080 or
            MAX_MATRIX_PROCESS_COUNT != 16640):
        fail("cross-fit cardinality changed")


@dataclass(frozen=True)
class Invocation:
    attempt: int
    root_index: int
    schedule_index: int
    K: Tuple[int, ...]
    phase: str

    def __post_init__(self) -> None:
        if type(self.attempt) is not int or self.attempt not in ATTEMPTS:
            fail("invocation attempt is outside the frozen uint8 domain")
        if type(self.root_index) is not int or \
                not 0 <= self.root_index < len(ROOTS) or \
                type(self.schedule_index) is not int or \
                not 0 <= self.schedule_index < len(SCHEDULES):
            fail("invocation loss coordinate is outside the frozen domain")
        if type(self.K) is not tuple or not self.K or \
                any(type(K) is not int or K not in K_VALUES for K in self.K) or \
                tuple(K for K in K_VALUES if K in self.K) != self.K:
            fail("invocation K batch is not a canonical frozen sub-roster")
        if self.phase not in ("anchor", "matrix", "replay"):
            fail("invocation phase is unknown")
        if self.phase == "anchor" and \
                (self.root_index != ANCHOR_ROOT_INDEX or
                 self.schedule_index != ANCHOR_SCHEDULE_INDEX or
                 len(self.K) != 1):
            fail("anchor invocation is not one canonical single-K probe")
        if self.phase == "matrix" and \
                self.root_index == ANCHOR_ROOT_INDEX and \
                self.schedule_index == ANCHOR_SCHEDULE_INDEX:
            fail("matrix invocation overlaps the anchor coordinate")

    @property
    def root(self) -> str:
        return ROOTS[self.root_index]

    @property
    def schedule(self) -> str:
        return SCHEDULES[self.schedule_index]

    @property
    def cell_ordinal(self) -> int:
        return self.root_index * len(SCHEDULES) + self.schedule_index

    def argv(self, bench: Path) -> List[str]:
        return [
            str(bench), "precodefail", "--N",
            ",".join(str(K) for K in self.K),
            "--bb-list", "2", "--overhead", "0", "--trials", "1",
            "--threads", "1", "--loss", "0.5", "--seed", self.root,
            "--schedule", self.schedule, "--heavy-family", "periodic",
            "--mix-count", "2", "--binary-dense-rows", "12",
            "--gf256-heavy-rows", "12", "--dense-anchors", "two07",
            "--paired-overhead-stream", "--full-payload-solve",
            "--cold-solve-wide-xor", "policy", "--exact-precode-attempt",
            str(self.attempt), "--exact-packet-attempt", str(self.attempt),
            "--construction-seed-basis", "production-profile",
        ]

    def identity(self) -> Mapping[str, Any]:
        return {
            "attempt": self.attempt,
            "root_index": self.root_index,
            "loss_root": self.root,
            "schedule_index": self.schedule_index,
            "schedule": self.schedule,
            "cell_ordinal": self.cell_ordinal,
            "K": list(self.K),
            "phase": self.phase,
        }


@dataclass(frozen=True)
class ProcessResult:
    invocation: Invocation
    command_sha256: str
    stdout_sha256: str
    returncode: int
    stdout: bytes
    stderr: bytes


def _parse_uint(text: str, maximum: int, context: str) -> int:
    if type(text) is not str or UINT.fullmatch(text) is None:
        fail("{} is not canonical unsigned decimal".format(context))
    if len(text) > len(str(maximum)):
        fail("{} exceeds its bound".format(context))
    try:
        value = int(text)
    except ValueError:
        fail("{} cannot be represented as an integer".format(context))
    if value > maximum:
        fail("{} exceeds its bound".format(context))
    return value


def _parse_integral_decimal(text: str, maximum: int, context: str) -> int:
    if type(text) is not str or THREE_DECIMAL.fullmatch(text) is None:
        fail("{} is not canonical three-place decimal".format(context))
    integer, fraction = text.split(".")
    if fraction != "000":
        fail("{} is not a nonnegative integral decimal".format(context))
    return _parse_uint(integer, maximum, context)


def _parse_nonnegative_decimal(text: str, maximum: int,
                               context: str) -> str:
    if type(text) is not str or THREE_DECIMAL.fullmatch(text) is None:
        fail("{} is not canonical three-place decimal".format(context))
    _parse_uint(text.split(".")[0], maximum, context)
    return text


def _parse_csv_line(line: str) -> List[str]:
    try:
        rows = list(csv.reader([line], strict=True))
    except csv.Error as exc:
        fail("precodefail CSV is malformed: {}".format(exc))
    if len(rows) != 1 or len(rows[0]) != len(screen.CSV_HEADER):
        fail("precodefail CSV row has the wrong field count")
    if line != ",".join(rows[0]):
        fail("precodefail CSV is not in canonical unquoted form")
    return rows[0]


def _parse_metadata(line: str, invocation: Invocation,
                    source_commit: str) -> Mapping[str, str]:
    prefix = "# precodefail: "
    if not line.startswith(prefix):
        fail("precodefail metadata line is missing")
    values: Dict[str, str] = {}
    for token in line[len(prefix):].split(" "):
        if token.count("=") != 1:
            fail("precodefail metadata token is malformed")
        key, value = token.split("=", 1)
        if not key or not value or key in values:
            fail("precodefail metadata is ambiguous")
        values[key] = value
    if set(values) != set(screen.METADATA_KEYS):
        fail("precodefail metadata fields changed")
    expected = {
        "trials": "1", "threads": "1", "loss": "0.5",
        "source_hits_override": "0", "packet_peel_seed_xor": "0x0",
        "binary_dense_rows_override": "12",
        "gf256_heavy_rows_override": "12",
        "dense_anchor_layout": "two07",
        "odd_packet_peel_seed_xor": "0x0",
        "packet_row_seed_multiplier": "0x1",
        "packet_row_seed_avalanche": "0",
        "seed_block_bytes_override": "0", "overhead_stream": "paired",
        "full_payload_solve": "1", "schedule": invocation.schedule,
        "cold_solve_wide_xor": "policy", "exact_attempt_mode": "1",
        "exact_precode_attempt": str(invocation.attempt),
        "exact_packet_attempt": str(invocation.attempt),
        "construction_seed_basis": "production-profile",
        "seed_schedule_sha256": ZERO_SHA256,
        "source_git_commit": source_commit,
    }
    for key, expected_value in expected.items():
        if values[key] != expected_value:
            fail("precodefail metadata {} changed".format(key))
    if values["seed"] != invocation.root:
        fail("precodefail metadata seed changed")
    return values


def _framed_lines(result: ProcessResult, source_commit: str) -> List[str]:
    expected_command_sha256 = sha256_json({
        "argv": result.invocation.argv(Path("wirehair_v2_bench"))[1:],
        "invocation": result.invocation.identity(),
    })
    if (SHA256.fullmatch(result.command_sha256) is None or
            result.command_sha256 != expected_command_sha256 or
            SHA256.fullmatch(result.stdout_sha256) is None or
            result.stdout_sha256 != sha256_bytes(result.stdout)):
        fail("precodefail process receipt is inconsistent")
    if (len(result.stdout) > MAX_STDOUT_BYTES or
            len(result.stderr) > MAX_STDERR_BYTES):
        fail("precodefail output exceeded its bound")
    if not result.stdout.endswith(b"\n") or b"\r" in result.stdout:
        fail("precodefail stdout framing changed")
    try:
        text = result.stdout.decode("ascii")
    except UnicodeDecodeError:
        fail("precodefail stdout is not ASCII")
    lines = text.splitlines()
    if len(lines) < 2:
        fail("precodefail omitted metadata or header")
    _parse_metadata(lines[0], result.invocation, source_commit)
    if tuple(_parse_csv_line(lines[1])) != screen.CSV_HEADER:
        fail("precodefail CSV header changed")
    return lines


def _terminal_fields(row: Mapping[str, str], K: int, attempt: int,
                     context: str) -> Mapping[str, Any]:
    expected_text = {
        "N": str(K), "bb": "2", "heavy_family": "periodic",
        "mix_count": "2", "binary_dense_rows": "12",
        "gf256_heavy_rows": "12", "dense_identity_corner": "0",
        "overhead": "0", "trials": "1", "seed_attempt": "",
        "active_packet_peel_seed_xor": "0x0",
        "precode_attempt": str(attempt), "packet_attempt": str(attempt),
        "attempt_mode": "exact",
        "construction_seed_basis": "production-profile",
        "seed_schedule_sha256": ZERO_SHA256,
    }
    for field, expected in expected_text.items():
        if row[field] != expected:
            fail("{} field {} changed".format(context, field))
    if (HEX64.fullmatch(row["effective_precode_seed"]) is None or
            HEX32.fullmatch(row["effective_packet_seed"]) is None):
        fail("{} effective seed is malformed".format(context))
    staircase = _parse_uint(row["staircase"], 400, context + " staircase")
    source_hits = _parse_uint(row["source_hits"], 8, context + " source hits")
    if staircase != STAIRCASE_BY_K[K] or \
            source_hits != (3 if K >= 10000 else 2) or \
            K + staircase + 24 > 65535:
        fail("{} realized structure changed".format(context))
    success = _parse_uint(row["success"], 1, context + " success")
    rank_fail = _parse_uint(row["rank_fail"], 1, context + " rank failure")
    error = _parse_uint(row["error"], 1, context + " error")
    if success + rank_fail + error != 1 or error != 0:
        fail("{} terminal result is invalid".format(context))
    expected_fail_rate = "0.00000000" if success else "1.00000000"
    if EIGHT_DECIMAL.fullmatch(row["fail_rate"]) is None or \
            row["fail_rate"] != expected_fail_rate:
        fail("{} failure rate disagrees with the terminal result".format(
            context))
    inactivated = _parse_integral_decimal(
        row["inact_mu"], 65535, context + " inactivation")
    if (_parse_uint(row["inact_max"], 65535,
                    context + " inactivation max") != inactivated or
            inactivated > K + staircase + 24):
        fail("{} inactivation counters disagree".format(context))
    binary_deficiency = _parse_integral_decimal(
        row["binary_def_mu"], 65535, context + " binary deficiency")
    gf256_rank_gain = _parse_integral_decimal(
        row["heavy_gain_mu"], 65535, context + " GF256 rank gain")
    if (_parse_uint(row["binary_def_max"], 65535,
                    context + " binary deficiency max") !=
            binary_deficiency or
            _parse_uint(row["heavy_gain_min"], 65535,
                        context + " GF256 rank gain min") != gf256_rank_gain or
            binary_deficiency > inactivated or gf256_rank_gain > 12 or
            gf256_rank_gain > binary_deficiency or
            bool(success) != (gf256_rank_gain == binary_deficiency) or
            bool(rank_fail) != (gf256_rank_gain < binary_deficiency)):
        fail("{} rank counters disagree with the terminal result".format(
            context))
    heavy_shortfall = _parse_uint(
        row["heavy_shortfall"], 1, context + " heavy shortfall")
    expected_shortfall = int(
        bool(rank_fail) and binary_deficiency <= 12 and
        gf256_rank_gain < binary_deficiency)
    if heavy_shortfall != expected_shortfall:
        fail("{} GF256 shortfall disagrees with the rank counters".format(
            context))
    if row["binary_def_hist"] != "{}:1".format(binary_deficiency) or \
            row["heavy_gain_hist"] != "{}:1".format(gf256_rank_gain):
        fail("{} rank histogram disagrees with the one-trial row".format(
            context))
    if row["first_rank_fail"] != ("0" if rank_fail else "-1") or \
            row["failure_trials"] != ("0" if rank_fail else ""):
        fail("{} failure diagnostics changed".format(context))
    for field in (
            "solve_ms_mu", "build_ms_mu", "peel_ms_mu", "project_ms_mu",
            "residual_ms_mu", "backsub_ms_mu"):
        _parse_nonnegative_decimal(
            row[field], UINT64_MAX, context + " " + field)
    return {
        "outcome": "success" if success else "rank_fail",
        "success": bool(success),
        "rank_fail": bool(rank_fail),
        "error": False,
        "binary_deficiency": binary_deficiency,
        "gf256_rank_gain": gf256_rank_gain,
        "inactivated_columns": inactivated,
        "heavy_shortfall": heavy_shortfall,
        "effective_precode_seed": row["effective_precode_seed"],
        "effective_packet_seed": row["effective_packet_seed"],
        "block_xors": _parse_integral_decimal(
            row["block_xors_mu"], UINT64_MAX, context + " block XORs"),
        "gf256_muladds": _parse_integral_decimal(
            row["block_muladds_mu"], UINT64_MAX,
            context + " GF256 muladds"),
    }


def parse_success_result(result: ProcessResult,
                         source_commit: str) -> Tuple[Mapping[str, Any], ...]:
    if result.returncode != 0 or result.stderr:
        fail("successful precodefail invocation did not exit cleanly")
    lines = _framed_lines(result, source_commit)
    if len(lines) != len(result.invocation.K) + 2:
        fail("precodefail emitted the wrong row count")
    parsed: List[Mapping[str, Any]] = []
    seen = set()
    realized_K = []
    for line_number, line in enumerate(lines[2:], 3):
        values = _parse_csv_line(line)
        row = dict(zip(screen.CSV_HEADER, values))
        K = _parse_uint(row["N"], 64000, "precodefail K")
        if K not in result.invocation.K or K in seen:
            fail("precodefail emitted an unexpected or duplicate K")
        seen.add(K)
        realized_K.append(K)
        terminal = _terminal_fields(
            row, K, result.invocation.attempt,
            "precodefail row {}".format(line_number))
        parsed.append({
            "K": K,
            "attempt": result.invocation.attempt,
            "root_index": result.invocation.root_index,
            "loss_root": result.invocation.root,
            "schedule_index": result.invocation.schedule_index,
            "schedule": result.invocation.schedule,
            "cell_ordinal": result.invocation.cell_ordinal,
            "trace_executed": True,
            "command_sha256": result.command_sha256,
            "stdout_sha256": result.stdout_sha256,
            "configuration_proof_sha256": None,
            **terminal,
        })
    if seen != set(result.invocation.K):
        fail("precodefail omitted a requested K")
    if tuple(realized_K) != result.invocation.K:
        fail("precodefail emitted K rows in a noncanonical order")
    return tuple(parsed)


def parse_configuration_failure(result: ProcessResult,
                                source_commit: str) -> Mapping[str, Any]:
    invocation = result.invocation
    if invocation.phase != "anchor" or len(invocation.K) != 1 or \
            invocation.cell_ordinal != 0:
        fail("configuration failure occurred outside the single-K anchor")
    lines = _framed_lines(result, source_commit)
    K = invocation.K[0]
    expected_stderr = (
        "precodefail configuration failed N={} bb=2 heavy_family=periodic "
        "mix_count=2 attempt_mode=exact precode_attempt={} "
        "packet_attempt={} result=2\n".format(
            K, invocation.attempt, invocation.attempt)).encode("ascii")
    if result.returncode != 2 or result.stderr != expected_stderr or \
            len(lines) != 2:
        fail("precodefail configuration failure is not canonical")
    return {
        "K": K,
        "attempt": invocation.attempt,
        "command_sha256": result.command_sha256,
        "stdout_sha256": result.stdout_sha256,
        "stderr_sha256": sha256_bytes(result.stderr),
    }


def _run_raw(invocation: Invocation, pinned: screen.PinnedExecutable) \
        -> ProcessResult:
    command = invocation.argv(pinned.path)
    try:
        completed = subprocess.run(
            command, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, timeout=PROCESS_TIMEOUT_SECONDS,
            check=False, shell=False, close_fds=True,
            start_new_session=True, env=CHILD_ENVIRONMENT,
            executable="/proc/self/fd/{}".format(pinned.descriptor),
            pass_fds=(pinned.descriptor,))
    except (OSError, subprocess.SubprocessError) as exc:
        fail("precodefail invocation failed to execute: {}".format(exc))
    return ProcessResult(
        invocation=invocation,
        command_sha256=sha256_json({
            "argv": command[1:], "invocation": invocation.identity()}),
        stdout_sha256=sha256_bytes(completed.stdout),
        returncode=completed.returncode,
        stdout=completed.stdout,
        stderr=completed.stderr,
    )


def _parallel(invocations: Sequence[Invocation], jobs: int,
              pinned: screen.PinnedExecutable) -> List[ProcessResult]:
    results: List[ProcessResult] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
        futures = [executor.submit(_run_raw, invocation, pinned)
                   for invocation in invocations]
        try:
            for future in concurrent.futures.as_completed(futures):
                results.append(future.result())
        except BaseException:
            for future in futures:
                future.cancel()
            raise
    return results


def _matrix_key(record: Mapping[str, Any]) -> Tuple[int, int, int]:
    return (int(record["K"]), int(record["attempt"]),
            int(record["cell_ordinal"]))


def _matrix_record_ordinal(K: int, attempt: int, cell_ordinal: int) -> int:
    return ((K_VALUES.index(K) * len(ATTEMPTS) + attempt) *
            (len(ROOTS) * len(SCHEDULES)) + cell_ordinal)


def build_matrix(pinned: screen.PinnedExecutable, source_commit: str,
                 jobs: int) -> Tuple[Dict[Tuple[int, int, int], Mapping[str, Any]],
                                     Mapping[str, int]]:
    anchor_invocations = tuple(
        Invocation(attempt, ANCHOR_ROOT_INDEX, ANCHOR_SCHEDULE_INDEX,
                   (K,), "anchor")
        for attempt in ATTEMPTS for K in K_VALUES)
    matrix: Dict[Tuple[int, int, int], Mapping[str, Any]] = {}
    constructible: Dict[int, List[int]] = {attempt: [] for attempt in ATTEMPTS}
    failures: Dict[Tuple[int, int], Mapping[str, Any]] = {}
    for result in _parallel(anchor_invocations, jobs, pinned):
        K = result.invocation.K[0]
        attempt = result.invocation.attempt
        if result.returncode == 0:
            rows = parse_success_result(result, source_commit)
            if len(rows) != 1:
                fail("single-K anchor did not emit one row")
            matrix[_matrix_key(rows[0])] = rows[0]
            constructible[attempt].append(K)
        else:
            failures[(K, attempt)] = parse_configuration_failure(
                result, source_commit)
    for attempt in ATTEMPTS:
        constructible[attempt].sort()
        if tuple(K for K in K_VALUES if (K, attempt) not in failures) != \
                tuple(constructible[attempt]):
            fail("anchor construction partition is incomplete")

    remaining_invocations = []
    for attempt in ATTEMPTS:
        values = tuple(constructible[attempt])
        if not values:
            continue
        for root_index in range(len(ROOTS)):
            for schedule_index in range(len(SCHEDULES)):
                if root_index == ANCHOR_ROOT_INDEX and \
                        schedule_index == ANCHOR_SCHEDULE_INDEX:
                    continue
                remaining_invocations.append(Invocation(
                    attempt, root_index, schedule_index, values, "matrix"))
    for result in _parallel(tuple(remaining_invocations), jobs, pinned):
        if result.returncode != 0:
            fail("constructible batched invocation failed")
        for row in parse_success_result(result, source_commit):
            key = _matrix_key(row)
            if key in matrix:
                fail("matrix contains a duplicate executed cell")
            matrix[key] = row

    for (K, attempt), proof in failures.items():
        for root_index, root in enumerate(ROOTS):
            for schedule_index, schedule in enumerate(SCHEDULES):
                cell_ordinal = root_index * len(SCHEDULES) + schedule_index
                key = (K, attempt, cell_ordinal)
                if key in matrix:
                    fail("construction failure overlaps an executed cell")
                matrix[key] = {
                    "K": K, "attempt": attempt,
                    "root_index": root_index, "loss_root": root,
                    "schedule_index": schedule_index, "schedule": schedule,
                    "cell_ordinal": cell_ordinal,
                    "trace_executed": False,
                    "outcome": "construct_failed", "success": False,
                    "rank_fail": False, "error": False,
                    "binary_deficiency": None, "gf256_rank_gain": None,
                    "inactivated_columns": None, "heavy_shortfall": None,
                    "effective_precode_seed": None,
                    "effective_packet_seed": None,
                    "block_xors": None, "gf256_muladds": None,
                    "command_sha256": proof["command_sha256"],
                    "stdout_sha256": proof["stdout_sha256"],
                    "configuration_proof_sha256": proof["stderr_sha256"],
                }
    if len(matrix) != MATRIX_CELL_COUNT:
        fail("matrix cardinality is incomplete")
    for K in K_VALUES:
        precode_bases = set()
        packet_bases = set()
        for attempt in ATTEMPTS:
            seeds = {
                (matrix[(K, attempt, cell)]["effective_precode_seed"],
                 matrix[(K, attempt, cell)]["effective_packet_seed"])
                for cell in range(len(ROOTS) * len(SCHEDULES))
                if matrix[(K, attempt, cell)]["trace_executed"]
            }
            if len(seeds) > 1:
                fail("effective construction seeds vary by loss cell")
            if seeds:
                precode_seed, packet_seed = next(iter(seeds))
                precode_bases.add(
                    (int(precode_seed, 16) -
                     attempt * screen.PRECODE_ATTEMPT_STRIDE) &
                    screen.MASK64)
                packet_bases.add(
                    (int(packet_seed, 16) -
                     attempt * screen.PACKET_ATTEMPT_STRIDE) &
                    screen.MASK32)
        if len(precode_bases) > 1 or len(packet_bases) > 1:
            fail("effective seeds do not follow the frozen attempt strides")
    return matrix, {
        "anchor_processes": len(anchor_invocations),
        "batch_processes": len(remaining_invocations),
        "configuration_failed_K_attempts": len(failures),
    }


def _training_cells(heldout_fold: int) -> Tuple[int, ...]:
    heldout = set(FOLDS[heldout_fold])
    return tuple(
        root_index * len(SCHEDULES) + schedule_index
        for root_index in range(len(ROOTS)) if root_index not in heldout
        for schedule_index in range(len(SCHEDULES)))


def _heldout_cells(heldout_fold: int) -> Tuple[int, ...]:
    return tuple(
        root_index * len(SCHEDULES) + schedule_index
        for root_index in FOLDS[heldout_fold]
        for schedule_index in range(len(SCHEDULES)))


def attempt_priority(
        matrix: Mapping[Tuple[int, int, int], Mapping[str, Any]],
        heldout_fold: int) -> Tuple[int, ...]:
    training = _training_cells(heldout_fold)
    scores = []
    for attempt in ATTEMPTS:
        weak_K = 0
        weak_cells = 0
        for K in K_VALUES:
            count = sum(
                not bool(matrix[(K, attempt, cell)]["success"])
                for cell in training)
            weak_K += int(count > 0)
            weak_cells += count
        scores.append((weak_K, weak_cells, attempt))
    return tuple(item[2] for item in sorted(scores))


def _select_attempt(
        matrix: Mapping[Tuple[int, int, int], Mapping[str, Any]], K: int,
        training: Sequence[int], priority: Sequence[int]) -> Optional[int]:
    for attempt in priority:
        if all(matrix[(K, attempt, cell)]["success"] for cell in training):
            return int(attempt)
    return None


def adjudicate(matrix: Mapping[Tuple[int, int, int], Mapping[str, Any]]) \
        -> Mapping[str, Any]:
    folds = []
    priority_selected: Dict[Tuple[int, int], int] = {}
    baseline_selected: Dict[Tuple[int, int], int] = {}
    priority_weak_total = 0
    baseline_weak_total = 0
    priority_unselectable_total = 0
    baseline_unselectable_total = 0
    for fold_index in range(len(FOLDS)):
        training = _training_cells(fold_index)
        heldout = _heldout_cells(fold_index)
        priority = attempt_priority(matrix, fold_index)
        fold_priority_weak = 0
        fold_baseline_weak = 0
        fold_priority_unselectable = []
        fold_baseline_unselectable = []
        for K in K_VALUES:
            selected = _select_attempt(matrix, K, training, priority)
            baseline = _select_attempt(matrix, K, training, ATTEMPTS)
            if selected is None:
                fold_priority_unselectable.append(K)
            else:
                priority_selected[(fold_index, K)] = selected
                fold_priority_weak += sum(
                    not matrix[(K, selected, cell)]["success"]
                    for cell in heldout)
            if baseline is None:
                fold_baseline_unselectable.append(K)
            else:
                baseline_selected[(fold_index, K)] = baseline
                fold_baseline_weak += sum(
                    not matrix[(K, baseline, cell)]["success"]
                    for cell in heldout)
        priority_weak_total += fold_priority_weak
        baseline_weak_total += fold_baseline_weak
        priority_unselectable_total += len(fold_priority_unselectable)
        baseline_unselectable_total += len(fold_baseline_unselectable)
        folds.append({
            "fold": fold_index,
            "heldout_root_indices": list(FOLDS[fold_index]),
            "priority_sha256": sha256_json(list(priority)),
            "selected_attempts_sha256": sha256_json([
                priority_selected.get((fold_index, K)) for K in K_VALUES]),
            "baseline_attempts_sha256": sha256_json([
                baseline_selected.get((fold_index, K)) for K in K_VALUES]),
            "priority_weak_cells": fold_priority_weak,
            "ascending_weak_cells": fold_baseline_weak,
            "priority_unselectable_K": fold_priority_unselectable,
            "ascending_unselectable_K": fold_baseline_unselectable,
        })
    disposition = (
        "PASS" if priority_unselectable_total == 0 and
        baseline_unselectable_total == 0 and priority_weak_total == 0 and
        priority_weak_total < baseline_weak_total else "REJECT")
    return {
        "disposition": disposition,
        "priority_weak_cells": priority_weak_total,
        "ascending_weak_cells": baseline_weak_total,
        "priority_unselectable_fold_K": priority_unselectable_total,
        "ascending_unselectable_fold_K": baseline_unselectable_total,
        "folds": folds,
        "priority_selected": priority_selected,
        "baseline_selected": baseline_selected,
    }


def replay_selected(
        matrix: Mapping[Tuple[int, int, int], Mapping[str, Any]],
        selected: Mapping[Tuple[int, int], int],
        pinned: screen.PinnedExecutable, source_commit: str, jobs: int) \
        -> Tuple[Mapping[str, Any], ...]:
    groups: Dict[Tuple[int, int, int, int], List[int]] = {}
    root_to_fold = {
        root_index: fold_index for fold_index, fold in enumerate(FOLDS)
        for root_index in fold}
    for fold_index, fold in enumerate(FOLDS):
        for K in K_VALUES:
            attempt = selected[(fold_index, K)]
            for root_index in fold:
                for schedule_index in range(len(SCHEDULES)):
                    groups.setdefault(
                        (fold_index, attempt, root_index, schedule_index),
                        []).append(K)
    invocations = tuple(
        Invocation(attempt, root_index, schedule_index, tuple(Ks), "replay")
        for (fold_index, attempt, root_index, schedule_index), Ks
        in sorted(groups.items()))
    replay = []
    seen = set()
    expected = {
        (fold_index, K,
         root_index * len(SCHEDULES) + schedule_index)
        for fold_index, fold in enumerate(FOLDS) for K in K_VALUES
        for root_index in fold
        for schedule_index in range(len(SCHEDULES))
    }
    for result in _parallel(invocations, jobs, pinned):
        if result.returncode != 0:
            fail("selected replay invocation did not exit cleanly")
        fold_index = root_to_fold[result.invocation.root_index]
        for row in parse_success_result(result, source_commit):
            K = int(row["K"])
            expected_fold = fold_index
            if selected[(expected_fold, K)] != int(row["attempt"]):
                fail("selected replay used another fold attempt")
            coordinate = (
                expected_fold, K, int(row["cell_ordinal"]))
            if coordinate in seen:
                fail("selected replay contains a duplicate coordinate")
            seen.add(coordinate)
            original = matrix[_matrix_key(row)]
            fields = (
                "outcome", "success", "rank_fail", "error",
                "binary_deficiency", "gf256_rank_gain",
                "inactivated_columns", "heavy_shortfall",
                "effective_precode_seed", "effective_packet_seed",
                "block_xors", "gf256_muladds",
            )
            if any(row[field] != original[field] for field in fields):
                fail("selected replay differs from the frozen matrix")
            replay.append({
                "schema": REPLAY_SCHEMA,
                "fold": expected_fold,
                "K": K,
                "attempt": row["attempt"],
                "cell_ordinal": row["cell_ordinal"],
                "matrix_record_ordinal": _matrix_record_ordinal(
                    K, int(row["attempt"]), int(row["cell_ordinal"])),
                "outcome": row["outcome"],
                "success": row["success"],
                "rank_fail": row["rank_fail"],
                "error": row["error"],
                "binary_deficiency": row["binary_deficiency"],
                "gf256_rank_gain": row["gf256_rank_gain"],
                "inactivated_columns": row["inactivated_columns"],
                "heavy_shortfall": row["heavy_shortfall"],
                "effective_precode_seed": row["effective_precode_seed"],
                "effective_packet_seed": row["effective_packet_seed"],
                "block_xors": row["block_xors"],
                "gf256_muladds": row["gf256_muladds"],
                "command_sha256": row["command_sha256"],
                "stdout_sha256": row["stdout_sha256"],
            })
    replay.sort(key=lambda row: (
        int(row["fold"]), K_VALUES.index(int(row["K"])),
        int(row["cell_ordinal"])))
    if len(replay) != OOF_CELL_COUNT or seen != expected:
        fail("selected replay does not cover the exact out-of-fold roster")
    return tuple(replay)


def _controller_receipt(source_commit: str) -> Mapping[str, Any]:
    if type(source_commit) is not str or \
            screen.COMMIT_TOKEN.fullmatch(source_commit) is None:
        fail("controller source commit is malformed")
    try:
        relative = CONTROLLER_PATH.relative_to(REPO_ROOT)
    except ValueError:
        fail("cross-fit controller is outside the repository")
    try:
        completed = subprocess.run(
            ["git", "status", "--porcelain=v1", "--untracked-files=all",
             "--", str(relative)], cwd=str(REPO_ROOT),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
            timeout=30)
    except (OSError, subprocess.SubprocessError) as exc:
        fail("cannot attest cross-fit controller: {}".format(exc))
    if completed.returncode != 0 or completed.stderr or completed.stdout:
        fail("cross-fit controller is not clean at the frozen commit")
    data = screen._stable_file_bytes(CONTROLLER_PATH, 4 * 1024 * 1024)
    return {
        "source_git_commit": source_commit,
        "path": str(relative),
        "sha256": sha256_bytes(data),
    }


def _canonical_file(value: Mapping[str, Any]) -> bytes:
    return (canonical_json(value) + "\n").encode("ascii")


def _write_jsonl(path: Path, rows: Iterable[Mapping[str, Any]]) -> str:
    digest = hashlib.sha256()
    try:
        descriptor = os.open(
            str(path), os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    except OSError as exc:
        fail("cannot create {}: {}".format(path, exc))
    try:
        for row in rows:
            data = (canonical_json(row) + "\n").encode("ascii")
            digest.update(data)
            offset = 0
            while offset < len(data):
                written = os.write(descriptor, data[offset:])
                if written <= 0:
                    fail("short write creating {}".format(path))
                offset += written
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    return digest.hexdigest()


def _matrix_rows(
        matrix: Mapping[Tuple[int, int, int], Mapping[str, Any]]) \
        -> Iterable[Mapping[str, Any]]:
    for K in K_VALUES:
        for attempt in ATTEMPTS:
            for cell_ordinal in range(len(ROOTS) * len(SCHEDULES)):
                row = dict(matrix[(K, attempt, cell_ordinal)])
                row["schema"] = RECORD_SCHEMA
                row["record_ordinal"] = _matrix_record_ordinal(
                    K, attempt, cell_ordinal)
                yield row


def run(bench: Path, output_dir: Path, jobs: int) -> Mapping[str, Any]:
    _validate_constants()
    contract = contract_description()
    if type(jobs) is not int or jobs < 1 or jobs > MAX_JOBS:
        fail("jobs must be in [1,{}]".format(MAX_JOBS))
    if output_dir.exists() or not output_dir.parent.is_dir():
        fail("output path is not one fresh path below an existing directory")
    base_source = screen._source_receipt()
    source_commit = base_source["source_git_commit"]
    source_receipt_sha256 = base_source["source_receipt_sha256"]
    controller = _controller_receipt(source_commit)
    temporary = Path(tempfile.mkdtemp(
        prefix=".wh2-mix2-attempt-crossfit-", dir=str(output_dir.parent)))
    published = False
    try:
        with screen._open_binary(bench, "wirehair_v2_bench") as pinned:
            screen.read_bench_description(
                pinned.path, source_commit, pinned.descriptor)
            input_value = {
                "schema": INPUT_SCHEMA,
                "contract": contract,
                "bench": pinned.receipt(),
                "source_receipt": base_source,
                "controller": controller,
                "jobs": jobs,
            }
            input_value["input_sha256"] = self_hash(
                input_value, "input_sha256")
            input_bytes = _canonical_file(input_value)
            screen._write_exclusive(
                temporary / "attempt-crossfit-input.json",
                input_bytes)
            matrix, process_counts = build_matrix(pinned, source_commit, jobs)
            matrix_sha256 = _write_jsonl(
                temporary / "attempt-crossfit-matrix.jsonl",
                _matrix_rows(matrix))
            verdict = adjudicate(matrix)
            if verdict["priority_unselectable_fold_K"] == 0:
                replay = replay_selected(
                    matrix, verdict["priority_selected"], pinned,
                    source_commit, jobs)
            else:
                replay = ()
            replay_sha256 = _write_jsonl(
                temporary / "attempt-crossfit-replay.jsonl", replay)
            if pinned.receipt() != input_value["bench"]:
                fail("benchmark executable changed during cross-fit")
            if screen._source_receipt() != base_source or \
                    _controller_receipt(source_commit) != controller:
                fail("source changed during cross-fit")
            summary = {
                "schema": SUMMARY_SCHEMA,
                "contract_sha256": contract["contract_sha256"],
                "input_sha256": input_value["input_sha256"],
                "input_file_sha256": sha256_bytes(input_bytes),
                "bench_binary_sha256": input_value["bench"]["sha256"],
                "source_git_commit": source_commit,
                "source_receipt_sha256": source_receipt_sha256,
                "controller_sha256": controller["sha256"],
                "candidate_profile_sha256":
                    EXPECTED_CANDIDATE_PROFILE_SHA256,
                "matrix_record_count": MATRIX_CELL_COUNT,
                "matrix_sha256": matrix_sha256,
                "replay_record_count": len(replay),
                "replay_sha256": replay_sha256,
                **process_counts,
                "matrix_processes":
                    process_counts["anchor_processes"] +
                    process_counts["batch_processes"],
                "replay_processes": len({
                    row["command_sha256"] for row in replay}),
                "priority_weak_cells": verdict["priority_weak_cells"],
                "ascending_weak_cells": verdict["ascending_weak_cells"],
                "priority_unselectable_fold_K":
                    verdict["priority_unselectable_fold_K"],
                "ascending_unselectable_fold_K":
                    verdict["ascending_unselectable_fold_K"],
                "priority_selected_fold_K_count":
                    len(verdict["priority_selected"]),
                "folds": verdict["folds"],
                "disposition": verdict["disposition"],
                "promotion_evidence": False,
                "fresh_roots_used": False,
            }
            summary["summary_sha256"] = self_hash(summary, "summary_sha256")
            screen._write_exclusive(
                temporary / "attempt-crossfit-summary.json",
                _canonical_file(summary))
        screen._publish_directory_noreplace(temporary, output_dir)
        published = True
        return summary
    finally:
        if not published and temporary.exists():
            shutil.rmtree(temporary)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("describe", help="print the frozen cross-fit contract")
    run_parser = subparsers.add_parser("run", help="run the consumed-root matrix")
    run_parser.add_argument("--bench", type=Path, required=True)
    run_parser.add_argument("--output-dir", type=Path, required=True)
    run_parser.add_argument("--jobs", type=int, default=12)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    arguments = _parser().parse_args(argv)
    try:
        _validate_constants()
        if arguments.command == "describe":
            print(canonical_json(contract_description()))
            return 0
        summary = run(arguments.bench, arguments.output_dir, arguments.jobs)
        print(canonical_json(summary))
        return 0 if summary["disposition"] == "PASS" else 2
    except (CrossfitError, screen.ShortScreenError) as exc:
        print("wh2 mix2 attempt cross-fit: {}".format(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
