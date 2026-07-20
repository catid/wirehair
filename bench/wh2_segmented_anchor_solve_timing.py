#!/usr/bin/env python3
"""Sealed isolated full-payload timing for raw segmented D12 anchors.

The panel compares the raw D12 two-anchor phase-zero baseline independently
against three-048 and four-0369.  Both comparisons retain mixed 10+2, P244,
q0, mix=2, and the same packet trace; no additional weak-K architecture seed
repair is applied beyond the identical baseline profile inputs.
Each codec process preflights both layouts and emits four paired ABBABAAB
cycles.  Cycle zero is discarded and cold/warm observations are sealed in the
same (K, block width, packet seed, schedule) bootstrap cluster.

Preparation freezes one exact committed source tree, native Release binary,
runner, isolation helper, thermal sampler, task manifest, build provenance,
and the completed raw all-K recovery receipt.  Run segments own the sole
CPU/DIMM sampler.  They refuse to coexist with load fillers, another I2C
reader, a cross-payload campaign, or any other wirehair benchmark process.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import fcntl
import hashlib
import importlib
import io
import json
import math
import os
from pathlib import Path
import random
import re
import shutil
import statistics
import subprocess
import sys
import time
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


sys.dont_write_bytecode = True

SCHEMA = "wirehair.wh2.segmented_anchor_solve_timing.v1"
CONTROLLER_NAME = "wh2_segmented_anchor_solve_timing.py"
BUILD_HELPER_NAME = "wh2_rank_floor_two_anchor_allk.py"
BUILD_SUPPORT_NAME = "wh2_rank_floor_two_anchor_screen.py"
ISOLATION_HELPER_NAME = "wh2_p32_dispatch_timing.py"
ISOLATION_HELPER_TEST_NAME = "wh2_p32_dispatch_timing_test.py"
SAMPLER_NAME = "wirehair_expo_thermal_sampler.py"
RECOVERY_RESULT_NAME = "wh2_segmented_anchor_allk_recovery_result.json"

P32_HELPER_SOURCE_PROVENANCE = {
    "source_commit": "97a3a0b941f3efe34a6e18609b7681d3a1110982",
    "source_tree": "ab84c084b212fecc6672a7efda5a4df6203c866c",
    "files": {
        "bench/wh2_p32_dispatch_timing.py": {
            "git_blob_sha1": "c1605efe36168f96bb45864e22a50f0e2fcbfbf8",
            "sha256": "1efbeffd05e7b4dbd9950b9eafc64016e7564f36c410a552c4e231362f6c4ba8",
        },
        "bench/wh2_p32_dispatch_timing_test.py": {
            "git_blob_sha1": "21a6c01bc196971c643809ae23fd329116447c31",
            "sha256": "22656b311658e3ed73fc5771c4ff8c12cb07e0b2f9e7c8b751e04291a5bf0fa5",
        },
        "bench/wirehair_expo_thermal_sampler.py": {
            "git_blob_sha1": "77d58141f3ecc6cf769e62eba474ffd4e4571b5e",
            "sha256": "047b6e7e70fc54631c94808bfa2aeb24470423ba3a71f376bcab9f78ab63bd56",
        },
    },
}

KS = (3200, 10000, 20000, 32000, 48466, 64000)
WIDTHS = (64, 1280, 4096)
SEEDS = (
    5141736193901540515,  # 0x475b1dd3864534a3
    7344707251978627308,  # 0x65eda454ed6be8ec
    8059702554173035811,  # 0x6fd9d0c558cad523
)
SCHEDULES = ("burst", "adversarial", "repair-only")
CACHE_STATES = ("cold", "warm")
CANDIDATES = ("three-048", "four-0369")
CLUSTER_SIZE = len(CACHE_STATES) * len(CANDIDATES)
TASK_COUNT = (len(KS) * len(WIDTHS) * len(SEEDS) * len(SCHEDULES) *
              CLUSTER_SIZE)
ORDER = "ABBABAAB"
OVERHEAD = 4
LOSS = "0.5"
PERIOD = 244
GROUPED_ROWS = 0
BUCKETS = "auto"
MAX_ENVIRONMENTAL_ATTEMPTS = 10
MAX_MINOR_FAULTS = 64
MAX_STDOUT_BYTES = 512 * 1024
MAX_STDERR_BYTES = 64 * 1024
MAX_INTERRUPTED_ARCHIVE_BYTES = 64 * 1024 * 1024
MAX_CPU_TEMP_C = 85.0
MAX_DIMM_TEMP_C = 84.0
MAX_THERMAL_GAP_S = 2.25
MAX_THERMAL_MARGIN_S = 2.25
QUIET_SAMPLE_SECONDS = 1.0
MAX_QUIET_BUSY_TICKS_PER_CPU = 1
BOOTSTRAP_REPETITIONS = 20000
PARETO_MATERIAL_SPEED_ADVANTAGE = 0.005
MALLOC_MMAP_THRESHOLD = "1073741824"
MALLOC_TRIM_THRESHOLD = "-1"

RECOVERY_VALIDATED_SUMMARY_SHA256 = (
    "ab4e8e2100829c490f45935c24593d3604c84a47993ad3a73ce220c0a595cff6")
RECOVERY_DATA_MANIFEST_SHA256 = (
    "36ea532889d1ed9fa7c8c1e536cec26c86c41fe78d5da9b80fe1b7648cac9e60")
RECOVERY_CANONICAL_SHA256 = (
    "69cad574501f6b658ba826ceab5a6b764be7e34c0500c929db193f4f393bb3a2")
RECOVERY_CONTROLLER_COMMIT = "2e7ae1ffa5ff9752a5ef820cae12f81f43eb16e9"
RECOVERY_PROTOTYPE_COMMIT = "57d79903716bfb7bb815106040e478ed6996a029"

# Exact GetDenseCount(K) values selected by the pinned baseline for this
# timing panel.  These are part of the architecture identity, not merely
# informational preamble fields.
STAIRCASE_ROWS = {
    3200: 62,
    10000: 86,
    20000: 134,
    32000: 190,
    48466: 374,
    64000: 346,
}

SPEED_POLICY = {
    "scope": "one aggregate gate per raw segmented candidate",
    "cluster_coordinates":
        "candidate,K,bb,seed_index,seed,schedule (cold/warm together)",
    "all_fixed_cells_common_success": True,
    "observed_ratio_of_sums_at_most": 1.01,
    "one_sided_95_upper_below": 1.01,
    "recovery_gate": "bound completed raw all-K evidence",
    "pareto_selection": {
        "comparison":
            "(three-048/cotimed baseline)/(four-0369/cotimed baseline)",
        "material_three_048_advantage_at_least":
            PARETO_MATERIAL_SPEED_ADVANTAGE,
        "one_sided_95_upper_below": 1.0,
        "tie_or_unsupported_advantage":
            "retain four-0369 for stronger raw recovery",
    },
}

UINT_RE = re.compile(r"0|[1-9][0-9]*\Z")
SINT_RE = re.compile(r"0|-?[1-9][0-9]*\Z")
HEX_RE = re.compile(r"0x(?:0|[1-9a-f][0-9a-f]*)\Z")
SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
TASK_DIR_RE = re.compile(r"[0-9]{4}\Z")
TXN_DIR_RE = re.compile(r"[0-9]{4}\.part\.[1-9][0-9]*\Z")
TRANSACTION_ARCHIVE_RE = re.compile(
    r"([0-9]{4}\.part\.[1-9][0-9]*)\.recovered\.[1-9][0-9]*\Z")
UNBOUND_TASK_ARCHIVE_RE = re.compile(
    r"(task-[0-9]{4}\.unbound)\.recovered\.[1-9][0-9]*\Z")
UTC_RE = re.compile(
    r"[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}\."
    r"[0-9]{3}Z\Z")

CSV_FIELDS = (
    "N", "bb", "overhead", "schedule", "seed", "loss", "cache_state",
    "cycle", "slot", "arm", "period", "grouped_rows",
    "buckets_requested", "seed_attempt", "matrix_seed", "peel_seed",
    "preflight_result", "cell_class", "common_success", "result",
    "outcome_stable", "elapsed_ns", "saturated", "cpu_before",
    "cpu_after", "cpu_migrated", "minflt_delta", "majflt_delta",
    "fault_contaminated", "inactivated", "binary_def", "heavy_gain",
    "block_xors", "block_muladds", "build_ns", "peel_ns", "project_ns",
    "residual_ns", "backsub_ns", "joint_source_xors",
    "joint_marginal_xors", "joint_marginal_copies", "joint_active_deltas",
    "joint_scratch_bytes", "dual_source_columns", "source_bytes",
    "packet_payload_bytes", "intermediate_bytes", "solve_value_arena_bytes",
    "solve_value_eager_zero_bytes", "solve_value_commit_copy_bytes",
)

PREAMBLE_FIELDS = (
    "schema", "policy", "timing_scope", "cycles", "order",
    "discard_cycle", "cycle_mode", "cycle_index", "N", "bb", "overhead",
    "loss", "seed", "schedule", "cache_state", "overhead_stream",
    "evict_bytes", "eviction_prefaulted", "control_period",
    "control_grouped_rows", "control_buckets", "control_grouped_hash_seed",
    "control_final_h_a_columns", "candidate_period",
    "candidate_grouped_rows", "candidate_buckets",
    "candidate_grouped_hash_seed", "candidate_final_h_a_columns",
    "gf256_rows", "gf16_rows", "control_dense_layout",
    "candidate_dense_layout", "dense_layout_is_only_architecture_selector",
    "control_staircase_rows", "control_dense_rows", "control_heavy_rows",
    "control_source_hits", "control_field", "control_dense_identity_corner",
    "control_dense_two_anchor_exact", "control_dense_two_anchor_phase",
    "control_segmented_dense_anchors", "control_heavy_family",
    "control_mix_count", "candidate_staircase_rows", "candidate_dense_rows",
    "candidate_heavy_rows", "candidate_source_hits", "candidate_field",
    "candidate_dense_identity_corner", "candidate_dense_two_anchor_exact",
    "candidate_dense_two_anchor_phase", "candidate_segmented_dense_anchors",
    "candidate_heavy_family", "candidate_mix_count", "control_attempt",
    "control_matrix_seed", "control_peel_seed", "candidate_attempt",
    "candidate_matrix_seed", "candidate_peel_seed", "mix", "payload",
    "payload_count", "payload_bytes", "payload_alignment",
    "payload_prefaulted", "system_build", "tls_reapply",
    "allocator_tls_state", "solve_value_storage", "solve_value_publish",
    "preflight_control_result", "preflight_candidate_result", "cell_class",
    "common_success", "trace_sha256",
)

WORK_FIELDS = (
    "result", "inactivated", "binary_def", "heavy_gain", "block_xors",
    "block_muladds", "joint_source_xors", "joint_marginal_xors",
    "joint_marginal_copies", "joint_active_deltas", "joint_scratch_bytes",
    "dual_source_columns", "source_bytes", "packet_payload_bytes",
    "intermediate_bytes", "solve_value_arena_bytes",
    "solve_value_eager_zero_bytes", "solve_value_commit_copy_bytes",
)


class TimingError(RuntimeError):
    """Substantive campaign, codec, or evidence failure."""


class ContaminationError(TimingError):
    """Retryable, explicitly receipted environmental timing contamination."""


def isolation_module():
    try:
        return importlib.import_module("wh2_p32_dispatch_timing")
    except ImportError as exc:
        raise TimingError("isolated timing support module is unavailable") from exc


def build_module():
    try:
        return importlib.import_module("wh2_rank_floor_two_anchor_allk")
    except ImportError as exc:
        raise TimingError("pinned-build support module is unavailable") from exc


def validate_isolation_policy(module) -> None:
    expected_limits = {
        "MAX_CPU_TEMP_C": MAX_CPU_TEMP_C,
        "MAX_DIMM_TEMP_C": MAX_DIMM_TEMP_C,
        "MAX_THERMAL_GAP_S": MAX_THERMAL_GAP_S,
        "MAX_THERMAL_MARGIN_S": MAX_THERMAL_MARGIN_S,
    }
    for name, expected in expected_limits.items():
        if getattr(module, name, None) != expected:
            raise TimingError("isolation helper policy differs: " + name)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(
        timespec="milliseconds").replace("+00:00", "Z")


def canonical_json(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                       allow_nan=False, ensure_ascii=True) + "\n").encode("ascii")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fsync_directory(path: Path) -> None:
    descriptor = os.open(str(path), os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def atomic_write(path: Path, data: bytes, mode: int = 0o444) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + ".part")
    if path.exists() or path.is_symlink() or part.exists() or part.is_symlink():
        raise TimingError("refusing to replace or recover output: %s" % path)
    with part.open("xb") as stream:
        stream.write(data)
        stream.flush()
        os.fsync(stream.fileno())
    os.chmod(part, mode)
    os.replace(part, path)
    fsync_directory(path.parent)


def recover_derived_atomic(path: Path, schema: str, name: str) -> str:
    """Finish or discard a crash remnant for replayable derived metadata."""
    part = path.with_name(path.name + ".part")
    for candidate in (path, part):
        if candidate.is_symlink() or \
                (candidate.exists() and not candidate.is_file()):
            raise TimingError(name + " atomic path is unsafe")
    if path.exists():
        if part.exists():
            if path.stat().st_size != part.stat().st_size or \
                    sha256_file(path) != sha256_file(part):
                raise TimingError(name + " atomic remnant disagrees with final")
            part.unlink()
            fsync_directory(path.parent)
            return "removed-identical-remnant"
        return "already-final"
    if not part.exists():
        return "absent"
    try:
        verify_sealed(load_canonical(part, name + " partial"), schema, name)
    except TimingError:
        # Campaign and summary receipts are wholly derived from immutable
        # task/segment evidence.  An incomplete pre-fsync temporary contains
        # no unique evidence and is safe to discard before exact regeneration.
        part.unlink()
        fsync_directory(path.parent)
        return "discarded-incomplete-remnant"
    os.replace(part, path)
    fsync_directory(path.parent)
    return "completed-valid-remnant"


def sealed_record(schema: str, payload: Mapping[str, object]) -> Dict[str, object]:
    if not schema or "schema" in payload or "self_sha256_excluding_field" in payload:
        raise TimingError("invalid sealed record")
    value: Dict[str, object] = {"schema": schema, **payload}
    value["self_sha256_excluding_field"] = sha256_bytes(canonical_json(value))
    return value


def verify_sealed(value: object, schema: str, name: str) -> Dict[str, object]:
    if not isinstance(value, dict) or value.get("schema") != schema:
        raise TimingError(name + " schema mismatch")
    claimed = value.get("self_sha256_excluding_field")
    if not isinstance(claimed, str) or SHA256_RE.fullmatch(claimed) is None:
        raise TimingError(name + " self hash is malformed")
    unhashed = dict(value)
    del unhashed["self_sha256_excluding_field"]
    if sha256_bytes(canonical_json(unhashed)) != claimed:
        raise TimingError(name + " self hash mismatch")
    return value


def load_canonical(path: Path, name: str) -> Dict[str, object]:
    if path.is_symlink() or not path.is_file() or \
            path.stat().st_size > MAX_INTERRUPTED_ARCHIVE_BYTES:
        raise TimingError(name + " is missing, unsafe, or oversized")
    raw = path.read_bytes()
    if len(raw) > MAX_INTERRUPTED_ARCHIVE_BYTES:
        raise TimingError(name + " exceeded its byte bound while reading")
    try:
        value = json.loads(raw.decode("ascii"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise TimingError(name + " is not ASCII JSON") from exc
    if not isinstance(value, dict) or canonical_json(value) != raw:
        raise TimingError(name + " is not canonical JSON")
    return value


def parse_uint(text: object, name: str, maximum: int = (1 << 64) - 1) -> int:
    if not isinstance(text, str) or UINT_RE.fullmatch(text) is None:
        raise TimingError(name + " is not a canonical uint")
    value = int(text)
    if value > maximum:
        raise TimingError(name + " exceeds its integer domain")
    return value


def parse_sint(text: object, name: str) -> int:
    if not isinstance(text, str) or SINT_RE.fullmatch(text) is None:
        raise TimingError(name + " is not a canonical signed integer")
    return int(text)


def validate_utc_timestamp(value: object, name: str) -> None:
    if not isinstance(value, str) or UTC_RE.fullmatch(value) is None:
        raise TimingError(name + " is not canonical millisecond UTC")
    try:
        parsed = datetime.fromisoformat(value[:-1] + "+00:00")
    except ValueError as exc:
        raise TimingError(name + " is not a valid UTC timestamp") from exc
    if parsed.tzinfo is None or parsed.utcoffset() != timezone.utc.utcoffset(parsed):
        raise TimingError(name + " is not UTC")


def _layout_receipt(layout: str) -> Dict[str, str]:
    if layout == "two-anchor":
        return {"dense_two_anchor_exact": "1", "segmented_dense_anchors": "none"}
    if layout not in CANDIDATES:
        raise TimingError("unknown raw dense layout")
    return {"dense_two_anchor_exact": "0", "segmented_dense_anchors": layout}


def generate_tasks() -> List[Dict[str, object]]:
    clusters: List[Tuple[str, List[Dict[str, object]]]] = []
    for K in KS:
        for bb in WIDTHS:
            for seed_index, seed in enumerate(SEEDS):
                for schedule in SCHEDULES:
                    cluster = {"K": K, "bb": bb, "seed_index": seed_index,
                               "seed": seed, "schedule": schedule}
                    priority = sha256_bytes(
                        b"wirehair.wh2.segmented-anchor.timing.cluster.v1\0" +
                        canonical_json(cluster))
                    first = ((CANDIDATES[0], "cold"), (CANDIDATES[1], "warm"),
                             (CANDIDATES[1], "cold"), (CANDIDATES[0], "warm"))
                    second = ((CANDIDATES[1], "cold"), (CANDIDATES[0], "warm"),
                              (CANDIDATES[0], "cold"), (CANDIDATES[1], "warm"))
                    order = first if int(priority[-1], 16) & 1 == 0 else second
                    members = [{**cluster, "candidate": candidate,
                                "cache_state": cache_state}
                               for candidate, cache_state in order]
                    clusters.append((priority, members))
    clusters.sort(key=lambda item: item[0])
    tasks: List[Dict[str, object]] = []
    for _priority, members in clusters:
        for member in members:
            coordinate = dict(member)
            task_id = sha256_bytes(
                b"wirehair.wh2.segmented-anchor.timing.task.v1\0" +
                canonical_json(coordinate))[:24]
            tasks.append({**coordinate, "job": len(tasks), "task_id": task_id})
    if len(tasks) != TASK_COUNT or TASK_COUNT != 648:
        raise TimingError("segmented timing Cartesian product changed")
    if len({str(task["task_id"]) for task in tasks}) != len(tasks):
        raise TimingError("segmented timing task identity collision")
    return tasks


def sanitized_environment(home: Path) -> Dict[str, str]:
    return {
        "HOME": str(home), "PATH": "/usr/bin:/bin", "LC_ALL": "C",
        "LANG": "C", "TZ": "UTC", "PYTHONDONTWRITEBYTECODE": "1",
        "MALLOC_MMAP_THRESHOLD_": MALLOC_MMAP_THRESHOLD,
        "MALLOC_TRIM_THRESHOLD_": MALLOC_TRIM_THRESHOLD,
    }


def command_for(design: Mapping[str, object], task: Mapping[str, object],
                *, cycle_index: Optional[int] = None,
                binary: Optional[Path] = None,
                evict_bytes: Optional[int] = None) -> List[str]:
    candidate = str(task["candidate"])
    if candidate not in CANDIDATES:
        raise TimingError("task candidate is unknown")
    root = Path(str(design["root"]))
    executable = binary or root / "frozen/wirehair_v2_bench"
    tools = design["tools"]
    command = [
        str(tools["taskset"]["path"]), "-c", str(design["core"]),
        str(tools["numactl"]["path"]),
        "--physcpubind=" + str(design["core"]),
        "--membind=" + str(design["numa_node"]), str(executable),
        "groupedtiming", "--N", str(task["K"]), "--bb", str(task["bb"]),
        "--overhead", str(OVERHEAD),
        "--control-period", str(PERIOD),
        "--control-grouped-rows", str(GROUPED_ROWS),
        "--control-buckets", BUCKETS,
        "--control-dense-layout", "two-anchor",
        "--candidate-period", str(PERIOD),
        "--candidate-grouped-rows", str(GROUPED_ROWS),
        "--candidate-buckets", BUCKETS,
        "--candidate-dense-layout", candidate,
        "--evict-bytes", str(design["evict_bytes"] if evict_bytes is None
                             else evict_bytes),
        "--cache-state", str(task["cache_state"]), "--loss", LOSS,
        "--seed", str(task["seed"]), "--schedule", str(task["schedule"]),
    ]
    if cycle_index is not None:
        command.extend(("--cycle-index", str(cycle_index)))
    return command


@dataclass(frozen=True)
class ParsedOutput:
    preamble: Mapping[str, str]
    rows: Tuple[Mapping[str, str], ...]
    stdout_sha256: str
    trace_sha256: str
    cell_class: str
    common_success: bool
    timed_control_ns: int
    timed_candidate_ns: int
    work_signatures: Mapping[str, Tuple[str, ...]]
    contaminations: Tuple[str, ...]


def parse_grouped_output(raw: bytes, task: Mapping[str, object], core: int,
                         *, replacement_cycle: Optional[int] = None,
                         expected_evict_bytes: Optional[int] = None,
                         ) -> ParsedOutput:
    if type(core) is not int or core < 0 or \
            (replacement_cycle is not None and
             (type(replacement_cycle) is not int or
              not 0 <= replacement_cycle <= 3)):
        raise TimingError("grouped timing parser domain is malformed")
    if len(raw) > MAX_STDOUT_BYTES or not raw.endswith(b"\n") or \
            b"\r" in raw or b"\0" in raw or b'"' in raw:
        raise TimingError("grouped timing output is noncanonical or oversized")
    try:
        lines = raw.decode("ascii").splitlines()
    except UnicodeDecodeError as exc:
        raise TimingError("grouped timing output is not ASCII") from exc
    expected_rows = 8 if replacement_cycle is not None else 32
    if len(lines) != expected_rows + 2:
        raise TimingError("grouped timing output line count changed")
    prefix = "# groupedtiming: "
    if not lines[0].startswith(prefix):
        raise TimingError("grouped timing preamble is missing")
    tokens = lines[0][len(prefix):].split(" ")
    if any(token.count("=") != 1 for token in tokens):
        raise TimingError("grouped timing preamble token is malformed")
    pairs = tuple(tuple(token.split("=", 1)) for token in tokens)
    if tuple(pair[0] for pair in pairs) != PREAMBLE_FIELDS:
        raise TimingError("grouped timing preamble schema/order changed")
    preamble = dict(pairs)
    candidate = str(task["candidate"])
    expected = {
        "schema": "v3", "policy": "h12-q0-grouped",
        "timing_scope": "solve",
        "cycles": "1" if replacement_cycle is not None else "4",
        "order": ORDER, "discard_cycle": "0",
        "cycle_mode": "replacement" if replacement_cycle is not None else "full",
        "cycle_index": str(replacement_cycle)
            if replacement_cycle is not None else "all",
        "N": str(task["K"]), "bb": str(task["bb"]),
        "overhead": str(OVERHEAD), "loss": LOSS, "seed": str(task["seed"]),
        "schedule": str(task["schedule"]),
        "cache_state": str(task["cache_state"]), "overhead_stream": "salted",
        "eviction_prefaulted": "1", "control_period": str(PERIOD),
        "control_grouped_rows": "0", "control_buckets": BUCKETS,
        "control_grouped_hash_seed": "0x0", "control_final_h_a_columns": "0",
        "candidate_period": str(PERIOD), "candidate_grouped_rows": "0",
        "candidate_buckets": BUCKETS, "candidate_grouped_hash_seed": "0x0",
        "candidate_final_h_a_columns": "0", "gf256_rows": "10",
        "gf16_rows": "2", "control_dense_layout": "two-anchor",
        "candidate_dense_layout": candidate,
        "dense_layout_is_only_architecture_selector": "1",
        "control_dense_rows": "12", "control_heavy_rows": "12",
        "control_field": "1", "control_dense_identity_corner": "0",
        "control_dense_two_anchor_exact": "1",
        "control_dense_two_anchor_phase": "0",
        "control_segmented_dense_anchors": "none",
        "control_heavy_family": "0", "control_mix_count": "2",
        "candidate_dense_rows": "12", "candidate_heavy_rows": "12",
        "candidate_field": "1", "candidate_dense_identity_corner": "0",
        "candidate_dense_two_anchor_phase": "0",
        "candidate_heavy_family": "0", "candidate_mix_count": "2",
        "mix": "2", "payload": "distinct-packet-zero-v1",
        "payload_count": str(int(task["K"]) + OVERHEAD),
        "payload_bytes": str((int(task["K"]) + OVERHEAD) * int(task["bb"])),
        "payload_alignment": "64", "payload_prefaulted": "1",
        "system_build": "outside-timer",
        "tls_reapply": "full-per-slot-outside-timer",
        "allocator_tls_state": "preflight-warmed",
        "solve_value_storage": "owned-noinit", "solve_value_publish": "swap",
    }
    expected.update({"candidate_" + key: value
                     for key, value in _layout_receipt(candidate).items()})
    staircase = STAIRCASE_ROWS.get(int(task["K"]))
    if staircase is None:
        raise TimingError("task K has no sealed staircase dimension")
    expected.update({
        "control_staircase_rows": str(staircase),
        "candidate_staircase_rows": str(staircase),
        "control_source_hits": "2" if int(task["K"]) < 10000 else "3",
        "candidate_source_hits": "2" if int(task["K"]) < 10000 else "3",
    })
    for key, value in expected.items():
        if preamble.get(key) != value:
            raise TimingError("preamble mismatch %s: %r != %r" %
                              (key, preamble.get(key), value))
    evict = parse_uint(preamble.get("evict_bytes"), "eviction bytes")
    if evict < 4096:
        raise TimingError("preamble eviction bound is too small")
    if expected_evict_bytes is not None and evict != expected_evict_bytes:
        raise TimingError("preamble eviction allocation changed")
    if preamble["control_staircase_rows"] != preamble["candidate_staircase_rows"] or \
            preamble["control_source_hits"] != preamble["candidate_source_hits"]:
        raise TimingError("raw layouts changed a non-layout dimension")
    for arm in ("control", "candidate"):
        parse_uint(preamble.get(arm + "_staircase_rows"), arm + " staircase")
        parse_uint(preamble.get(arm + "_source_hits"), arm + " source hits")
        parse_uint(preamble.get(arm + "_attempt"), arm + " attempt", 255)
        for suffix, maximum in (("matrix_seed", (1 << 64) - 1),
                                ("peel_seed", (1 << 32) - 1)):
            text = preamble.get(arm + "_" + suffix, "")
            if HEX_RE.fullmatch(text) is None or int(text, 16) > maximum:
                raise TimingError("malformed %s %s" % (arm, suffix))
        parse_uint(preamble.get("preflight_" + arm + "_result"),
                   arm + " preflight result", 1)
    trace = preamble.get("trace_sha256")
    if not isinstance(trace, str) or SHA256_RE.fullmatch(trace) is None:
        raise TimingError("trace hash is malformed")
    cell_class = preamble.get("cell_class", "")
    classes = {"common-success", "control-only", "candidate-only", "common-failure"}
    if cell_class not in classes:
        raise TimingError("cell class is malformed")
    control_result = parse_uint(
        preamble.get("preflight_control_result"), "control preflight result")
    candidate_result = parse_uint(
        preamble.get("preflight_candidate_result"), "candidate preflight result")
    expected_class = (
        "common-success" if control_result == 0 and candidate_result == 0 else
        "control-only" if control_result == 0 else
        "candidate-only" if candidate_result == 0 else "common-failure")
    if cell_class != expected_class:
        raise TimingError("cell class disagrees with preflight outcomes")
    common_success = expected_class == "common-success"
    if preamble.get("common_success") != ("1" if common_success else "0"):
        raise TimingError("common-success receipt disagrees with cell class")
    reader = csv.DictReader(io.StringIO("\n".join(lines[1:]) + "\n"))
    if tuple(reader.fieldnames or ()) != CSV_FIELDS:
        raise TimingError("grouped timing CSV schema changed")
    rows = tuple(dict(row) for row in reader)
    if len(rows) != expected_rows:
        raise TimingError("grouped timing row count changed")

    signatures: Dict[str, set[Tuple[str, ...]]] = {"control": set(), "candidate": set()}
    timed = {"control": 0, "candidate": 0}
    contaminations: List[str] = []
    first_cycle = replacement_cycle if replacement_cycle is not None else 0
    for index, row in enumerate(rows):
        if None in row or tuple(row) != CSV_FIELDS or \
                any(value is None for value in row.values()):
            raise TimingError("grouped timing CSV row width changed")
        cycle = first_cycle + index // 8
        slot = index % 8
        arm = "control" if ORDER[slot] == "A" else "candidate"
        exact = {
            "N": str(task["K"]), "bb": str(task["bb"]),
            "overhead": str(OVERHEAD), "schedule": str(task["schedule"]),
            "seed": str(task["seed"]), "loss": LOSS,
            "cache_state": str(task["cache_state"]), "cycle": str(cycle),
            "slot": str(slot), "arm": arm, "period": str(PERIOD),
            "grouped_rows": "0", "buckets_requested": BUCKETS,
            "seed_attempt": preamble[arm + "_attempt"],
            "matrix_seed": preamble[arm + "_matrix_seed"],
            "peel_seed": preamble[arm + "_peel_seed"],
            "preflight_result": preamble["preflight_" + arm + "_result"],
            "cell_class": cell_class,
            "common_success": "1" if common_success else "0",
            "outcome_stable": "1",
            "source_bytes": str(int(task["K"]) * int(task["bb"])),
            "packet_payload_bytes": str(
                (int(task["K"]) + OVERHEAD) * int(task["bb"])),
            "intermediate_bytes": str(
                (int(task["K"]) + staircase + 24) * int(task["bb"])),
            "solve_value_arena_bytes": str(
                (int(task["K"]) + staircase + 24) * int(task["bb"])),
            "solve_value_eager_zero_bytes": "0",
            "solve_value_commit_copy_bytes": "0",
        }
        for key, value in exact.items():
            if row.get(key) != value:
                raise TimingError("row %d mismatch %s" % (index, key))
        for field in WORK_FIELDS[1:] + (
                "elapsed_ns", "saturated", "build_ns", "peel_ns",
                "project_ns", "residual_ns", "backsub_ns"):
            parse_uint(row.get(field), "row %d %s" % (index, field))
        result = parse_uint(row.get("result"), "row result", 1)
        if result != parse_uint(
                row["preflight_result"], "row preflight result", 1):
            raise TimingError("row outcome differs from its preflight")
        elapsed = parse_uint(row["elapsed_ns"], "elapsed time")
        if elapsed == 0:
            raise TimingError("timed solve reported zero nanoseconds")
        if replacement_cycle is None and cycle != 0:
            timed[arm] += elapsed
        minflt = parse_sint(row.get("minflt_delta"), "minor faults")
        majflt = parse_sint(row.get("majflt_delta"), "major faults")
        if minflt < -1 or majflt < -1:
            raise TimingError("fault counter delta is outside the emitted domain")
        expected_fault = -1 if minflt < 0 or majflt < 0 else int(bool(minflt or majflt))
        if parse_sint(row.get("fault_contaminated"), "fault receipt") != expected_fault:
            raise TimingError("fault contamination receipt is inconsistent")
        before = parse_sint(row.get("cpu_before"), "CPU before")
        after = parse_sint(row.get("cpu_after"), "CPU after")
        migrated = parse_sint(row.get("cpu_migrated"), "CPU migration")
        expected_migration = -1 if before < 0 or after < 0 else int(before != after)
        if migrated != expected_migration:
            raise TimingError("CPU migration receipt is inconsistent")
        if parse_uint(row.get("saturated"), "saturated", 1):
            contaminations.append("row%d:saturated" % index)
        if before != core or after != core or migrated != 0:
            contaminations.append("row%d:migration:%d:%d:%d" %
                                  (index, before, after, migrated))
        if minflt < 0 or minflt > MAX_MINOR_FAULTS:
            contaminations.append("row%d:minor-fault:%d" % (index, minflt))
        if majflt != 0:
            contaminations.append("row%d:major-fault:%d" % (index, majflt))
        arena = parse_uint(row["solve_value_arena_bytes"], "arena bytes")
        intermediate = parse_uint(row["intermediate_bytes"], "intermediate bytes")
        if arena != intermediate or \
                parse_uint(row["solve_value_eager_zero_bytes"], "eager zero") != 0 or \
                parse_uint(row["solve_value_commit_copy_bytes"], "commit copy") != 0:
            raise TimingError("no-init solve arena receipt changed")
        signatures[arm].add(tuple(row[field] for field in WORK_FIELDS))
    if any(len(values) != 1 for values in signatures.values()):
        raise TimingError("repeated observations changed deterministic arm work")
    if replacement_cycle is not None:
        timed = {"control": 0, "candidate": 0}
    elif timed["control"] <= 0 or timed["candidate"] <= 0:
        raise TimingError("retained timing rows are empty")
    return ParsedOutput(
        preamble=preamble, rows=rows, stdout_sha256=sha256_bytes(raw),
        trace_sha256=trace, cell_class=cell_class,
        common_success=common_success, timed_control_ns=timed["control"],
        timed_candidate_ns=timed["candidate"],
        work_signatures={arm: next(iter(values))
                         for arm, values in signatures.items()},
        contaminations=tuple(contaminations),
    )


def task_summary(parsed: ParsedOutput) -> Dict[str, object]:
    return {
        "stdout_sha256": parsed.stdout_sha256,
        "trace_sha256": parsed.trace_sha256,
        "cell_class": parsed.cell_class,
        "common_success": parsed.common_success,
        "control_ns": parsed.timed_control_ns,
        "candidate_ns": parsed.timed_candidate_ns,
        "control_work_sha256": sha256_bytes(canonical_json(
            {"fields": list(WORK_FIELDS),
             "values": list(parsed.work_signatures["control"])})),
        "candidate_work_sha256": sha256_bytes(canonical_json(
            {"fields": list(WORK_FIELDS),
             "values": list(parsed.work_signatures["candidate"])})),
        "control_attempt": parse_uint(
            parsed.preamble["control_attempt"], "control attempt"),
        "control_matrix_seed": parsed.preamble["control_matrix_seed"],
        "control_peel_seed": parsed.preamble["control_peel_seed"],
        "candidate_attempt": parse_uint(
            parsed.preamble["candidate_attempt"], "candidate attempt"),
        "candidate_matrix_seed": parsed.preamble["candidate_matrix_seed"],
        "candidate_peel_seed": parsed.preamble["candidate_peel_seed"],
    }


def _clustered_bootstrap(values: Sequence[Mapping[str, object]],
                         candidate: str, domain: str) -> Dict[str, object]:
    if candidate not in CANDIDATES or not values:
        raise TimingError("bootstrap candidate/domain is empty")
    clusters: Dict[Tuple[object, ...], List[int]] = {}
    for value in values:
        if value.get("candidate") != candidate:
            raise TimingError("bootstrap mixed candidates")
        key = (value["K"], value["bb"], value["seed_index"],
               value["seed"], value["schedule"])
        totals = clusters.setdefault(key, [0, 0])
        control = value.get("control_ns")
        contender = value.get("candidate_ns")
        if type(control) is not int or type(contender) is not int or \
                control <= 0 or contender <= 0:
            raise TimingError("bootstrap timing value is malformed")
        totals[0] += control
        totals[1] += contender
    records = [tuple(clusters[key]) for key in sorted(clusters)]
    seed = int.from_bytes(hashlib.sha256(
        b"wirehair.wh2.segmented-anchor.bootstrap.v1\0" +
        candidate.encode("ascii") + b"\0" + domain.encode("ascii")
    ).digest()[:8], "big")
    generator = random.Random(seed)
    ratios: List[float] = []
    for _repetition in range(BOOTSTRAP_REPETITIONS):
        control_total = 0
        candidate_total = 0
        for _draw in range(len(records)):
            control, contender = records[generator.randrange(len(records))]
            control_total += control
            candidate_total += contender
        ratios.append(candidate_total / control_total)
    ratios.sort()
    return {
        "schema": "paired-cluster-ratio-of-sums-bootstrap-v1",
        "cluster_coordinates": SPEED_POLICY["cluster_coordinates"],
        "cluster_count": len(records), "seed": seed,
        "repetitions": BOOTSTRAP_REPETITIONS,
        "lower_95": format(
            ratios[int(math.floor(0.025 * (BOOTSTRAP_REPETITIONS - 1)))],
            ".12f"),
        "upper_95": format(
            ratios[int(math.ceil(0.975 * (BOOTSTRAP_REPETITIONS - 1)))],
            ".12f"),
        "upper_one_sided_95": format(
            ratios[int(math.ceil(0.95 * (BOOTSTRAP_REPETITIONS - 1)))],
            ".12f"),
    }


def _pareto_bootstrap(values: Sequence[Mapping[str, object]]) -> Dict[str, object]:
    if not values:
        raise TimingError("Pareto bootstrap domain is empty")
    clusters: Dict[Tuple[object, ...], List[int]] = {}
    fields = ("three_control_ns", "three_candidate_ns",
              "four_control_ns", "four_candidate_ns")
    for value in values:
        key = (value["K"], value["bb"], value["seed_index"],
               value["seed"], value["schedule"])
        totals = clusters.setdefault(key, [0, 0, 0, 0])
        for index, field in enumerate(fields):
            sample = value.get(field)
            if type(sample) is not int or sample <= 0:
                raise TimingError("Pareto bootstrap timing value is malformed")
            totals[index] += sample
    records = [tuple(clusters[key]) for key in sorted(clusters)]
    seed = int.from_bytes(hashlib.sha256(
        b"wirehair.wh2.segmented-anchor.pareto-bootstrap.v1"
    ).digest()[:8], "big")
    generator = random.Random(seed)
    ratios: List[float] = []
    for _repetition in range(BOOTSTRAP_REPETITIONS):
        totals = [0, 0, 0, 0]
        for _draw in range(len(records)):
            record = records[generator.randrange(len(records))]
            for index, sample in enumerate(record):
                totals[index] += sample
        three_control, three_candidate, four_control, four_candidate = totals
        ratios.append(
            (three_candidate * four_control) /
            (three_control * four_candidate))
    ratios.sort()
    return {
        "schema": "paired-cluster-baseline-normalized-bootstrap-v1",
        "cluster_coordinates":
            "K,bb,seed_index,seed,schedule (cold/warm together)",
        "ratio_direction":
            "(three-048/control)/(four-0369/control)",
        "cluster_count": len(records), "seed": seed,
        "repetitions": BOOTSTRAP_REPETITIONS,
        "lower_95": format(
            ratios[int(math.floor(0.025 * (BOOTSTRAP_REPETITIONS - 1)))],
            ".12f"),
        "upper_95": format(
            ratios[int(math.ceil(0.975 * (BOOTSTRAP_REPETITIONS - 1)))],
            ".12f"),
        "upper_one_sided_95": format(
            ratios[int(math.ceil(0.95 * (BOOTSTRAP_REPETITIONS - 1)))],
            ".12f"),
    }


def _speed_gate(summary: Mapping[str, object]) -> Dict[str, bool]:
    outcomes = summary.get("outcome_counts")
    count = summary.get("cell_count")
    ratio = summary.get("ratio_of_sums")
    bootstrap = summary.get("bootstrap")
    if not isinstance(outcomes, dict) or type(count) is not int:
        raise TimingError("speed summary is malformed")
    all_common = outcomes.get("common-success") == count and all(
        outcomes.get(name) == 0 for name in
        ("control-only", "candidate-only", "common-failure"))
    observed = isinstance(ratio, (int, float)) and not isinstance(ratio, bool) and \
        math.isfinite(float(ratio)) and float(ratio) <= 1.01
    if bootstrap is None:
        bounded = False
    elif isinstance(bootstrap, dict):
        try:
            upper = float(bootstrap["upper_one_sided_95"])
        except (KeyError, TypeError, ValueError) as exc:
            raise TimingError("speed bootstrap is malformed") from exc
        bounded = math.isfinite(upper) and upper < 1.01
    else:
        raise TimingError("speed bootstrap is malformed")
    return {
        "all_fixed_cells_common_success": all_common,
        "observed_regression_at_most_1_percent": observed,
        "one_sided_95_upper_below_1_01": bounded,
        "speed_gate_passed": all_common and observed and bounded,
    }


def candidate_pareto_comparison(
        cells: Sequence[Mapping[str, object]]) -> Dict[str, object]:
    by_coordinate: Dict[Tuple[object, ...], Dict[str, Mapping[str, object]]] = {}
    for cell in cells:
        candidate = str(cell.get("candidate"))
        if candidate not in CANDIDATES:
            raise TimingError("Pareto comparison contains an unknown candidate")
        key = (cell["K"], cell["bb"], cell["seed_index"], cell["seed"],
               cell["schedule"], cell["cache_state"])
        candidates = by_coordinate.setdefault(key, {})
        if candidate in candidates:
            raise TimingError("Pareto comparison contains a duplicate cell")
        candidates[candidate] = cell
    expected_pairs = (len(KS) * len(WIDTHS) * len(SEEDS) * len(SCHEDULES) *
                      len(CACHE_STATES))
    if len(by_coordinate) != expected_pairs or any(
            set(values) != set(CANDIDATES)
            for values in by_coordinate.values()):
        raise TimingError("Pareto comparison cell pairing is incomplete")
    common: List[Dict[str, object]] = []
    for key in sorted(by_coordinate):
        three = by_coordinate[key]["three-048"]
        four = by_coordinate[key]["four-0369"]
        if three.get("cell_class") != "common-success" or \
                four.get("cell_class") != "common-success":
            continue
        three_control = three.get("control_ns")
        three_ns = three.get("candidate_ns")
        four_control = four.get("control_ns")
        four_ns = four.get("candidate_ns")
        if any(type(value) is not int or value <= 0 for value in
               (three_control, three_ns, four_control, four_ns)):
            raise TimingError("Pareto comparison timing is malformed")
        K, bb, seed_index, seed, schedule, cache_state = key
        common.append({
            "candidate": "three-048", "K": K, "bb": bb,
            "seed_index": seed_index, "seed": seed,
            "schedule": schedule, "cache_state": cache_state,
            "three_control_ns": three_control,
            "three_candidate_ns": three_ns,
            "four_control_ns": four_control,
            "four_candidate_ns": four_ns,
        })
    three_control_total = sum(int(value["three_control_ns"])
                              for value in common)
    three_total = sum(int(value["three_candidate_ns"]) for value in common)
    four_control_total = sum(int(value["four_control_ns"])
                             for value in common)
    four_total = sum(int(value["four_candidate_ns"]) for value in common)
    raw_ratio = three_total / four_total if four_total else None
    ratio = ((three_total * four_control_total) /
             (three_control_total * four_total)) \
        if three_control_total and four_total else None
    bootstrap = _pareto_bootstrap(common) if common else None
    all_common = len(common) == expected_pairs
    materially_faster = ratio is not None and \
        ratio <= 1.0 - PARETO_MATERIAL_SPEED_ADVANTAGE
    statistically_supported = bootstrap is not None and \
        float(bootstrap["upper_one_sided_95"]) < 1.0
    return {
        "schema": "paired-candidate-pareto-v1",
        "ratio_direction": "(three-048/control)/(four-0369/control)",
        "cell_pair_count": expected_pairs,
        "paired_common_success_count": len(common),
        "all_fixed_cells_paired_common_success": all_common,
        "three_048_control_ns": three_control_total,
        "three_048_candidate_ns": three_total,
        "four_0369_control_ns": four_control_total,
        "four_0369_candidate_ns": four_total,
        "raw_candidate_ratio_of_sums": raw_ratio,
        "baseline_normalized_ratio_of_ratios": ratio,
        "bootstrap": bootstrap,
        "material_speed_advantage_required":
            PARETO_MATERIAL_SPEED_ADVANTAGE,
        "three_048_materially_faster_observed": materially_faster,
        "three_048_statistically_faster_one_sided_95":
            statistically_supported,
        "select_three_048_over_four_0369":
            all_common and materially_faster and statistically_supported,
    }


def select_pareto_architecture(
        gates: Mapping[str, Mapping[str, object]],
        pareto: Mapping[str, object]) -> Optional[str]:
    try:
        four_ok = gates["four-0369"]["speed_gate_passed"] is True
        three_ok = gates["three-048"]["speed_gate_passed"] is True
    except (KeyError, TypeError) as exc:
        raise TimingError("candidate speed-gate ledger is malformed") from exc
    if four_ok and three_ok:
        return "three-048" if pareto.get(
            "select_three_048_over_four_0369") is True else "four-0369"
    if four_ok:
        return "four-0369"
    if three_ok:
        return "three-048"
    return None


def aggregate_rows(tasks: Sequence[Mapping[str, object]],
                   parsed_by_job: Mapping[int, ParsedOutput]) -> Dict[str, object]:
    if set(parsed_by_job) != set(range(len(tasks))):
        raise TimingError("aggregate task ledger is incomplete")
    cells: List[Dict[str, object]] = []
    cache_identity: Dict[Tuple[object, ...], Tuple[object, ...]] = {}
    baseline_identity: Dict[Tuple[object, ...], Tuple[object, ...]] = {}
    work_fields = (
        "inactivated", "binary_def", "heavy_gain", "block_xors",
        "block_muladds", "joint_source_xors", "joint_marginal_xors",
        "joint_marginal_copies", "joint_active_deltas",
        "joint_scratch_bytes", "dual_source_columns", "intermediate_bytes",
    )
    for task in tasks:
        parsed = parsed_by_job[int(task["job"])]
        control_row = next(row for row in parsed.rows if row["arm"] == "control")
        candidate_row = next(row for row in parsed.rows if row["arm"] == "candidate")
        cache_key = (task["candidate"], task["K"], task["bb"],
                     task["seed_index"], task["seed"], task["schedule"])
        identity = (
            parsed.cell_class, parsed.trace_sha256,
            parsed.work_signatures["control"],
            parsed.work_signatures["candidate"],
            parsed.preamble["control_attempt"],
            parsed.preamble["control_matrix_seed"],
            parsed.preamble["control_peel_seed"],
            parsed.preamble["candidate_attempt"],
            parsed.preamble["candidate_matrix_seed"],
            parsed.preamble["candidate_peel_seed"],
        )
        previous = cache_identity.setdefault(cache_key, identity)
        if previous != identity:
            raise TimingError("cold/warm graph, trace, outcome, or work changed")
        baseline_key = (task["K"], task["bb"], task["seed_index"],
                        task["seed"], task["schedule"], task["cache_state"])
        baseline = (
            parsed.trace_sha256, parsed.work_signatures["control"],
            parsed.preamble["control_attempt"],
            parsed.preamble["control_matrix_seed"],
            parsed.preamble["control_peel_seed"],
        )
        previous_baseline = baseline_identity.setdefault(baseline_key, baseline)
        if previous_baseline != baseline:
            raise TimingError("baseline identity changed across candidate comparisons")
        cells.append({
            "job": task["job"], "candidate": task["candidate"],
            "K": task["K"], "bb": task["bb"],
            "seed_index": task["seed_index"], "seed": task["seed"],
            "schedule": task["schedule"], "cache_state": task["cache_state"],
            "cell_class": parsed.cell_class,
            "control_ns": parsed.timed_control_ns if parsed.common_success else None,
            "candidate_ns": parsed.timed_candidate_ns if parsed.common_success else None,
            "control_work": {field: int(control_row[field]) for field in work_fields},
            "candidate_work": {
                field: int(candidate_row[field]) for field in work_fields},
        })

    def summarize(values: Sequence[Mapping[str, object]],
                  candidate: str, domain: str,
                  *, bootstrap: bool) -> Dict[str, object]:
        common = [value for value in values
                  if value["cell_class"] == "common-success"]
        control_total = sum(int(value["control_ns"]) for value in common)
        candidate_total = sum(int(value["candidate_ns"]) for value in common)
        ratios = [int(value["candidate_ns"]) / int(value["control_ns"])
                  for value in common]
        outcomes = {name: sum(value["cell_class"] == name for value in values)
                    for name in ("common-success", "control-only",
                                 "candidate-only", "common-failure")}
        work_totals = {
            arm: {field: sum(int(value[arm + "_work"][field]) for value in values)
                  for field in work_fields}
            for arm in ("control", "candidate")
        }
        result: Dict[str, object] = {
            "cell_count": len(values), "outcome_counts": outcomes,
            "common_success_timing_count": len(common),
            "control_ns": control_total, "candidate_ns": candidate_total,
            "ratio_of_sums": candidate_total / control_total
                if control_total else None,
            "median_cell_ratio": statistics.median(ratios) if ratios else None,
            "candidate_faster_cells": sum(ratio < 1.0 for ratio in ratios),
            "work_totals": work_totals,
        }
        if bootstrap:
            result["bootstrap"] = _clustered_bootstrap(
                common, candidate, domain) if common else None
            result["speed_gate"] = _speed_gate(result)
        return result

    candidates: Dict[str, object] = {}
    for candidate in CANDIDATES:
        selected = [cell for cell in cells if cell["candidate"] == candidate]
        overall = summarize(selected, candidate, "overall", bootstrap=True)
        groups: Dict[str, object] = {}
        for K in KS:
            for bb in WIDTHS:
                domain = "K=%d,bb=%d" % (K, bb)
                values = [cell for cell in selected
                          if cell["K"] == K and cell["bb"] == bb]
                groups[domain] = summarize(
                    values, candidate, domain, bootstrap=True)
        for cache_state in CACHE_STATES:
            domain = "cache=%s" % cache_state
            groups[domain] = summarize(
                [cell for cell in selected
                 if cell["cache_state"] == cache_state],
                candidate, domain, bootstrap=False)
        candidates[candidate] = {
            "overall": overall, "groups": groups,
            "weak_coordinates": [{
                "K": cell["K"], "bb": cell["bb"],
                "seed_index": cell["seed_index"], "seed": cell["seed"],
                "schedule": cell["schedule"],
                "cell_class": cell["cell_class"],
            } for cell in selected if cell["cache_state"] == "cold" and
                cell["cell_class"] != "common-success"],
        }
    pareto = candidate_pareto_comparison(cells)
    return {
        "speed_policy": SPEED_POLICY, "candidates": candidates,
        "candidate_pareto_comparison": pareto,
        "baseline_identity_equal_across_candidates": True,
        "cold_warm_identity_exact": True,
    }


def checked_text(argv: Sequence[str]) -> str:
    try:
        result = subprocess.run(
            list(argv), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            timeout=30.0, check=False)
    except subprocess.TimeoutExpired as exc:
        raise TimingError("command timed out: %r" % (list(argv),)) from exc
    except OSError as exc:
        raise TimingError("command could not be executed: %r" %
                          (list(argv),)) from exc
    if len(result.stdout) > 1024 * 1024 or len(result.stderr) > 64 * 1024:
        raise TimingError("command output exceeded its preparation bound")
    if result.returncode != 0 or result.stderr:
        raise TimingError("command failed exit=%d argv=%r stderr=%r" %
                          (result.returncode, list(argv), result.stderr[:1000]))
    try:
        return result.stdout.decode("utf-8").strip()
    except UnicodeDecodeError as exc:
        raise TimingError("command output is not UTF-8") from exc


def _git_value(git: Path, repo: Path, *arguments: str) -> str:
    return checked_text((str(git), "-C", str(repo), *arguments))


def _validate_imported_helper_provenance(repo: Path, git: Path) -> None:
    provenance = P32_HELPER_SOURCE_PROVENANCE
    commit = str(provenance["source_commit"])
    tree = str(provenance["source_tree"])
    if _git_value(git, repo, "rev-parse", commit + "^{commit}") != commit or \
            _git_value(git, repo, "rev-parse", commit + "^{tree}") != tree:
        raise TimingError("imported isolation-helper source identity changed")
    files = provenance.get("files")
    if not isinstance(files, dict) or not files:
        raise TimingError("imported isolation-helper file ledger is malformed")
    for relative, receipt in files.items():
        path = repo / relative
        if not isinstance(relative, str) or not isinstance(receipt, dict) or \
                set(receipt) != {"git_blob_sha1", "sha256"} or \
                path.is_symlink() or not path.is_file() or \
                _git_value(git, repo, "rev-parse", "HEAD:" + relative) != \
                    receipt.get("git_blob_sha1") or \
                sha256_file(path) != receipt.get("sha256"):
            raise TimingError(
                "imported isolation-helper blob changed: " + str(relative))


def _validate_recovery_result(value: object) -> Dict[str, object]:
    if not isinstance(value, dict) or value.get("schema") != \
            "wirehair.wh2.segmented_anchor_allk.compact_recovery_result.v1":
        raise TimingError("raw all-K recovery receipt schema changed")
    if sha256_bytes(canonical_json(value)) != RECOVERY_CANONICAL_SHA256:
        raise TimingError("raw all-K recovery receipt content changed")
    evidence = value.get("evidence")
    architecture = value.get("architecture")
    arms = value.get("arms")
    selection = value.get("selection")
    if not all(isinstance(item, dict)
               for item in (evidence, architecture, arms, selection)):
        raise TimingError("raw all-K recovery receipt is incomplete")
    expected_evidence = {
        "prototype_commit": RECOVERY_PROTOTYPE_COMMIT,
        "controller_commit": RECOVERY_CONTROLLER_COMMIT,
        "validated_summary_sha256": RECOVERY_VALIDATED_SUMMARY_SHA256,
        "data_manifest_sha256": RECOVERY_DATA_MANIFEST_SHA256,
        "K_range": [2, 64000], "cells_per_arm": 575991,
        "total_cells": 2303964, "raw_architectures": True,
        "seed_fixes_applied": False, "timing_claim": False,
    }
    for key, expected in expected_evidence.items():
        if evidence.get(key) != expected:
            raise TimingError("raw recovery evidence changed: " + key)
    if architecture != {
            "dense_rows": 12, "gf256_rows": 10, "gf16_rows": 2,
            "period": 244, "mix_count": 2, "overhead": 0,
            "loss": "0.5", "block_bytes": 64,
            "baseline": "D12 two-anchor phase0 q0",
            "candidates": ["three-048", "three-059", "four-0369"],
    }:
        raise TimingError("raw recovery architecture changed")
    expected_failures = {"baseline": 733, "three-048": 663,
                         "four-0369": 650}
    for arm, failures in expected_failures.items():
        if not isinstance(arms.get(arm), dict) or \
                arms[arm].get("failures") != failures:
            raise TimingError("raw recovery arm result changed: " + arm)
    if selection.get("raw_recovery_winner") != "four-0369" or \
            selection.get("seed_fix_evaluation_deferred") is not True or \
            selection.get("pinned_full_payload_timing_required") is not True or \
            selection.get("timing_acceptance_limit_percent") != 1:
        raise TimingError("raw recovery selection contract changed")
    return value


def _copy_frozen(source: Path, destination: Path, mode: int) -> None:
    if source.is_symlink() or not source.is_file():
        raise TimingError("frozen source is not a regular file: %s" % source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() or destination.is_symlink():
        raise TimingError("frozen destination already exists: %s" % destination)
    shutil.copyfile(source, destination)
    os.chmod(destination, mode)


def _boundary_smoke(staging: Path, design: Mapping[str, object], binary: Path,
                    environment: Mapping[str, str]) -> Dict[str, object]:
    isolation = isolation_module()
    specs = (
        (3200, 64, "three-048"), (3200, 64, "four-0369"),
        (64000, 1280, "three-048"), (64000, 1280, "four-0369"),
    )
    records = []
    for index, (K, bb, candidate) in enumerate(specs):
        task = {
            "K": K, "bb": bb, "seed_index": 0, "seed": SEEDS[0],
            "schedule": "burst", "candidate": candidate,
            "cache_state": "warm",
        }
        argv = command_for(
            design, task, cycle_index=2, binary=binary, evict_bytes=4096)
        started = time.monotonic_ns()
        rc, stdout, stderr = isolation.run_bounded(
            argv, environment, 900.0, MAX_STDOUT_BYTES, MAX_STDERR_BYTES)
        ended = time.monotonic_ns()
        if rc != 0 or stderr:
            raise TimingError("prepare smoke failed exit=%d stderr=%r" %
                              (rc, stderr[:1000]))
        parsed = parse_grouped_output(
            stdout, task, int(design["core"]), replacement_cycle=2,
            expected_evict_bytes=4096)
        if not parsed.common_success:
            raise TimingError("prepare smoke is not common-success")
        name = "smoke.%d.%s.K%d.bb%d.csv" % (index, candidate, K, bb)
        atomic_write(staging / "provenance" / name, stdout)
        records.append({
            "task": task,
            "executed_argv": argv,
            "replay_argv": command_for(
                design, task, cycle_index=2, evict_bytes=4096),
            "executed_binary_sha256": sha256_file(binary),
            "stdout_name": name,
            "stdout_sha256": parsed.stdout_sha256,
            "trace_sha256": parsed.trace_sha256,
            "cell_class": parsed.cell_class,
            "nonpromotional_contaminations": list(parsed.contaminations),
            "started_monotonic_ns": started, "ended_monotonic_ns": ended,
            "duration_ns": ended - started,
        })
    return {
        "timing_evidence": False,
        "scope": "bounded CLI/architecture smoke under ordinary machine load",
        "records": records,
    }


def prepare_campaign(args: argparse.Namespace) -> None:
    common = build_module()
    isolation = isolation_module()
    result = Path(args.result_dir).resolve()
    repo = Path(args.repo).resolve()
    binary_hint = Path(args.binary).resolve()
    if result.exists() or result.is_symlink():
        raise TimingError("result directory already exists")
    if not repo.is_dir() or args.build_jobs <= 0 or args.core < 0 or \
            args.controller_core < 0 or args.thermal_core < 0 or \
            args.numa_node < 0 or args.evict_bytes < 268435456:
        raise TimingError("prepare arguments are outside the sealed domain")
    validate_isolation_policy(isolation)
    tool_names = ("git", "cmake", "taskset", "numactl", "sudo", "fuser",
                  "timeout", "python3", "env", "true")
    tools = {name: common._trusted_executable(name) for name in tool_names}
    try:
        sudo_probe = subprocess.run(
            (tools["sudo"]["path"], "-n", tools["true"]["path"]),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=10.0,
            check=False, env={"PATH": "/usr/bin:/bin", "LC_ALL": "C",
                              "LANG": "C", "TZ": "UTC"})
    except subprocess.TimeoutExpired as exc:
        raise TimingError("passwordless sudo preflight timed out") from exc
    except OSError as exc:
        raise TimingError("passwordless sudo preflight could not execute") from exc
    if sudo_probe.returncode != 0 or sudo_probe.stdout or sudo_probe.stderr:
        raise TimingError("passwordless sudo preflight failed")
    git = Path(tools["git"]["path"])
    head = _git_value(git, repo, "rev-parse", "HEAD^{commit}")
    tree = _git_value(git, repo, "rev-parse", "HEAD^{tree}")
    if re.fullmatch(r"[0-9a-f]{40}", head) is None or \
            re.fullmatch(r"[0-9a-f]{40}", tree) is None or \
            _git_value(git, repo, "status", "--porcelain", "--untracked-files=no"):
        raise TimingError("tracked source tree must be a clean exact commit")
    _validate_imported_helper_provenance(repo, git)
    source_dir = Path(__file__).resolve().parent
    sources = (
        Path(__file__).resolve(), source_dir / BUILD_HELPER_NAME,
        source_dir / BUILD_SUPPORT_NAME,
        source_dir / ISOLATION_HELPER_NAME,
        source_dir / ISOLATION_HELPER_TEST_NAME, source_dir / SAMPLER_NAME,
        source_dir / RECOVERY_RESULT_NAME,
    )
    if any(path.is_symlink() or not path.is_file() for path in sources):
        raise TimingError("one required committed timing source is unavailable")
    timing_topology = isolation.topology_record(args.core, args.numa_node)
    controller_topology = isolation.topology_record(
        args.controller_core, args.numa_node)
    thermal_topology = isolation.topology_record(args.thermal_core, args.numa_node)
    timing_quiet_cpus = timing_topology["llc_shared_cpus"]
    if not isinstance(timing_quiet_cpus, list) or not timing_quiet_cpus or \
            any(type(cpu) is not int or cpu < 0 for cpu in timing_quiet_cpus) or \
            timing_quiet_cpus != sorted(set(timing_quiet_cpus)) or \
            args.core not in timing_quiet_cpus:
        raise TimingError("timing LLC CPU topology is malformed")
    timing_llc = set(timing_quiet_cpus)
    if (args.controller_core in timing_llc or args.thermal_core in timing_llc or
            args.controller_core == args.thermal_core or
            args.controller_core in thermal_topology["siblings"] or
            args.thermal_core in controller_topology["siblings"]):
        raise TimingError("controller and thermal cores are not isolated")
    tasks = generate_tasks()
    manifest = b"".join(canonical_json(task) for task in tasks)
    recovery = _validate_recovery_result(
        json.loads((source_dir / RECOVERY_RESULT_NAME).read_text(encoding="utf-8")))

    with common.prepare_result_staging(result) as (staging, final):
        for directory in ("frozen", "provenance", "external", "ledger",
                          ".transactions", "interrupted", "segments"):
            (staging / directory).mkdir()
        smoke_home = staging / ".smoke-home"
        smoke_home.mkdir()
        atomic_write(staging / "controller.lock", b"", 0o644)
        atomic_write(staging / "tasks_manifest.jsonl", manifest)
        with common.fresh_pinned_benchmark_build(
                repo, head, binary_hint, tools["git"], tools["cmake"],
                tools["taskset"], build_workers=args.build_jobs,
                cpu_set=args.build_cpu_set) as fresh:
            frozen_binary = staging / "frozen/wirehair_v2_bench"
            frozen_cache = staging / "frozen/CMakeCache.txt"
            frozen_build = staging / "frozen/pinned_build.json"
            _copy_frozen(fresh.binary, frozen_binary, 0o555)
            _copy_frozen(fresh.cache, frozen_cache, 0o444)
            atomic_write(frozen_build, canonical_json(fresh.record))
            common.validate_pinned_build_record(
                json.loads(frozen_build.read_text(encoding="ascii")), frozen_cache)
            provisional_design = {
                "root": str(final), "tools": tools, "core": args.core,
                "numa_node": args.numa_node, "evict_bytes": args.evict_bytes,
            }
            prelaunch = _boundary_smoke(
                staging, provisional_design, frozen_binary,
                sanitized_environment(smoke_home))
        if smoke_home.is_symlink():
            raise TimingError("prepare smoke HOME became unsafe")
        shutil.rmtree(smoke_home)
        (staging / "runtime-home").mkdir()
        frozen_pairs = []
        for source in sources[:-1]:
            destination = staging / "frozen" / source.name
            _copy_frozen(source, destination, 0o555)
            frozen_pairs.append((source, destination))
        recovery_destination = staging / "external" / RECOVERY_RESULT_NAME
        _copy_frozen(sources[-1], recovery_destination, 0o444)
        frozen_pairs.append((sources[-1], recovery_destination))
        common.verify_frozen_sources_at_commit(repo, head, frozen_pairs, git)
        atomic_write(staging / "prelaunch_receipt.json", canonical_json(
            sealed_record(SCHEMA + ".prelaunch", {
                "created_utc": utc_now(), "boundary_smoke": prelaunch,
                "quiet_timing_launched": False,
                "recovery_receipt_sha256": sha256_file(recovery_destination),
                "recovery_validated_summary_sha256":
                    recovery["evidence"]["validated_summary_sha256"],
                "recovery_data_manifest_sha256":
                    recovery["evidence"]["data_manifest_sha256"],
            })))
        immutable: Dict[str, str] = {}
        for dirname in ("frozen", "provenance", "external"):
            for path in sorted((staging / dirname).iterdir()):
                if path.is_symlink() or not path.is_file():
                    raise TimingError("immutable prepare inventory is unsafe")
                immutable[str(path.relative_to(staging))] = sha256_file(path)
        design = sealed_record(SCHEMA + ".design", {
            "root": str(final), "head": head, "tree": tree,
            "timing_scope": "full-payload decoder precode solve",
            "raw_architecture": True, "architecture_seed_repairs": "none",
            "architecture": {
                "control": "D12 two-anchor phase0 q0",
                "candidates": list(CANDIDATES), "dense_rows": 12,
                "gf256_rows": 10, "gf16_rows": 2, "period": PERIOD,
                "grouped_rows": GROUPED_ROWS, "buckets": BUCKETS,
                "mix_count": 2,
                "dense_layout_is_only_architecture_selector": True,
            },
            "K": list(KS), "bb": list(WIDTHS), "seeds": list(SEEDS),
            "schedules": list(SCHEDULES), "cache_states": list(CACHE_STATES),
            "overhead": OVERHEAD, "loss": LOSS, "order": ORDER,
            "discard_cycle": 0, "task_count": len(tasks),
            "rows_per_task": 32, "timed_rows_per_task": 24,
            "physical_solve_rows": len(tasks) * 32,
            "retained_solve_rows": len(tasks) * 24,
            "cluster_size": CLUSTER_SIZE,
            "cluster_policy": SPEED_POLICY["cluster_coordinates"],
            "speed_policy": SPEED_POLICY,
            "max_environmental_attempts": MAX_ENVIRONMENTAL_ATTEMPTS,
            "minor_fault_max": MAX_MINOR_FAULTS,
            "output_bounds": {"stdout": MAX_STDOUT_BYTES,
                              "stderr": MAX_STDERR_BYTES},
            "allocator_environment": {
                "MALLOC_MMAP_THRESHOLD_": MALLOC_MMAP_THRESHOLD,
                "MALLOC_TRIM_THRESHOLD_": MALLOC_TRIM_THRESHOLD,
            },
            "core": args.core, "controller_core": args.controller_core,
            "thermal_core": args.thermal_core, "numa_node": args.numa_node,
            "evict_bytes": args.evict_bytes,
            "timing_quiet_cpus": timing_quiet_cpus,
            "timing_topology": timing_topology,
            "controller_topology": controller_topology,
            "thermal_topology": thermal_topology,
            "thermal_limits": {"cpu_c": MAX_CPU_TEMP_C,
                "dimm_c": MAX_DIMM_TEMP_C, "max_gap_s": MAX_THERMAL_GAP_S,
                "max_margin_s": MAX_THERMAL_MARGIN_S,
                "reject_dimm_at_or_above_limit": True},
            "thermal_policy":
                "controller-owned sole reader through post-end sample and graceful stop",
            "transaction_policy":
                "one atomic immutable task directory; bounded resumable segments",
            "imported_p32_helper_source_provenance":
                P32_HELPER_SOURCE_PROVENANCE,
            "runtime_llc_peer_policy": {
                "scope": "every timing child from immediately before launch through reap",
                "excluded_busy_cpu": args.core,
                "maximum_busy_ticks_per_peer":
                    MAX_QUIET_BUSY_TICKS_PER_CPU,
            },
            "launch_prerequisites": [
                "cross-payload v3 and every other wirehair benchmark process stopped",
                "nice-19 filler supervisor and workers stopped",
                "no pre-existing I2C reader or external thermal sampler",
                "every CPU sharing the timing LLC quiet and topology/power policy unchanged",
                "controller and owned thermal sampler outside timing LLC",
                "passwordless sudo available for sole-reader sampler lifecycle",
                "owned live thermal samples below CPU 85C and DIMM 84C with zero EDAC errors",
            ],
            "tasks_manifest_sha256": sha256_bytes(manifest),
            "controller_lock_sha256": sha256_bytes(b""),
            "prelaunch_receipt_sha256": sha256_file(
                staging / "prelaunch_receipt.json"),
            "recovery_receipt_sha256": sha256_file(recovery_destination),
            "recovery_validated_summary_sha256":
                RECOVERY_VALIDATED_SUMMARY_SHA256,
            "recovery_data_manifest_sha256": RECOVERY_DATA_MANIFEST_SHA256,
            "immutable_files": immutable, "tools": tools,
            "pinned_build_sha256": immutable["frozen/pinned_build.json"],
            "binary_sha256": immutable["frozen/wirehair_v2_bench"],
            "python": common.python_runtime_identity(),
        })
        atomic_write(staging / "design.json", canonical_json(design))
        prepare = sealed_record(SCHEMA + ".prepare", {
            "prepared_utc": utc_now(), "head": head, "tree": tree,
            "design_sha256": sha256_file(staging / "design.json"),
            "tasks_manifest_sha256": sha256_bytes(manifest),
            "prelaunch_receipt_sha256": sha256_file(
                staging / "prelaunch_receipt.json"),
            "binary_sha256": design["binary_sha256"],
            "runner_sha256": immutable["frozen/" + CONTROLLER_NAME],
            "isolation_helper_sha256": immutable[
                "frozen/" + ISOLATION_HELPER_NAME],
            "thermal_sampler_sha256": immutable["frozen/" + SAMPLER_NAME],
        })
        atomic_write(staging / "prepare_receipt.json", canonical_json(prepare))
        for dirname in ("frozen", "provenance", "external"):
            os.chmod(staging / dirname, 0o555)
        if _git_value(git, repo, "rev-parse", "HEAD^{commit}") != head or \
                _git_value(git, repo, "rev-parse", "HEAD^{tree}") != tree or \
                _git_value(git, repo, "status", "--porcelain",
                           "--untracked-files=no"):
            raise TimingError("repository changed during preparation")
    print(json.dumps({
        "result_dir": str(result), "head": head, "tree": tree,
        "task_count": len(tasks), "physical_solve_rows": len(tasks) * 32,
        "retained_solve_rows": len(tasks) * 24,
        "design_sha256": sha256_file(result / "design.json"),
        "manifest_sha256": sha256_file(result / "tasks_manifest.jsonl"),
        "binary_sha256": design["binary_sha256"],
        "quiet_timing_launched": False,
    }, sort_keys=True))


def load_design(root: Path) -> Dict[str, object]:
    common = build_module()
    validate_isolation_policy(isolation_module())
    design = verify_sealed(
        load_canonical(root / "design.json", "timing design"),
        SCHEMA + ".design", "timing design")
    if design.get("root") != str(root.resolve()):
        raise TimingError("prepared timing root moved")
    expected = {
        "timing_scope": "full-payload decoder precode solve",
        "raw_architecture": True, "architecture_seed_repairs": "none",
        "architecture": {
            "control": "D12 two-anchor phase0 q0",
            "candidates": list(CANDIDATES), "dense_rows": 12,
            "gf256_rows": 10, "gf16_rows": 2, "period": PERIOD,
            "grouped_rows": GROUPED_ROWS, "buckets": BUCKETS,
            "mix_count": 2,
            "dense_layout_is_only_architecture_selector": True,
        },
        "K": list(KS), "bb": list(WIDTHS), "seeds": list(SEEDS),
        "schedules": list(SCHEDULES), "cache_states": list(CACHE_STATES),
        "overhead": OVERHEAD, "loss": LOSS, "order": ORDER,
        "discard_cycle": 0, "task_count": TASK_COUNT, "rows_per_task": 32,
        "timed_rows_per_task": 24, "physical_solve_rows": TASK_COUNT * 32,
        "retained_solve_rows": TASK_COUNT * 24,
        "cluster_size": CLUSTER_SIZE,
        "cluster_policy": SPEED_POLICY["cluster_coordinates"],
        "speed_policy": SPEED_POLICY,
        "max_environmental_attempts": MAX_ENVIRONMENTAL_ATTEMPTS,
        "minor_fault_max": MAX_MINOR_FAULTS,
        "output_bounds": {"stdout": MAX_STDOUT_BYTES,
                          "stderr": MAX_STDERR_BYTES},
        "allocator_environment": {
            "MALLOC_MMAP_THRESHOLD_": MALLOC_MMAP_THRESHOLD,
            "MALLOC_TRIM_THRESHOLD_": MALLOC_TRIM_THRESHOLD,
        },
        "thermal_limits": {"cpu_c": MAX_CPU_TEMP_C,
            "dimm_c": MAX_DIMM_TEMP_C, "max_gap_s": MAX_THERMAL_GAP_S,
            "max_margin_s": MAX_THERMAL_MARGIN_S,
            "reject_dimm_at_or_above_limit": True},
        "thermal_policy":
            "controller-owned sole reader through post-end sample and graceful stop",
        "transaction_policy":
            "one atomic immutable task directory; bounded resumable segments",
        "imported_p32_helper_source_provenance":
            P32_HELPER_SOURCE_PROVENANCE,
        "runtime_llc_peer_policy": {
            "scope": "every timing child from immediately before launch through reap",
            "excluded_busy_cpu": design.get("core"),
            "maximum_busy_ticks_per_peer": MAX_QUIET_BUSY_TICKS_PER_CPU,
        },
        "launch_prerequisites": [
            "cross-payload v3 and every other wirehair benchmark process stopped",
            "nice-19 filler supervisor and workers stopped",
            "no pre-existing I2C reader or external thermal sampler",
            "every CPU sharing the timing LLC quiet and topology/power policy unchanged",
            "controller and owned thermal sampler outside timing LLC",
            "passwordless sudo available for sole-reader sampler lifecycle",
            "owned live thermal samples below CPU 85C and DIMM 84C with zero EDAC errors",
        ],
        "controller_lock_sha256": sha256_bytes(b""),
        "recovery_validated_summary_sha256": RECOVERY_VALIDATED_SUMMARY_SHA256,
        "recovery_data_manifest_sha256": RECOVERY_DATA_MANIFEST_SHA256,
    }
    dynamic_fields = {
        "schema", "self_sha256_excluding_field", "root", "head", "tree",
        "core", "controller_core", "thermal_core", "numa_node",
        "evict_bytes", "timing_quiet_cpus", "timing_topology", "controller_topology",
        "thermal_topology", "tasks_manifest_sha256",
        "prelaunch_receipt_sha256", "recovery_receipt_sha256",
        "immutable_files", "tools", "pinned_build_sha256", "binary_sha256",
        "python",
    }
    if set(design) != set(expected) | dynamic_fields:
        raise TimingError("frozen timing design field set changed")
    for key, expected_value in expected.items():
        if design.get(key) != expected_value:
            raise TimingError("frozen timing design changed: " + key)
    for key in ("core", "controller_core", "thermal_core", "numa_node",
                "evict_bytes"):
        if type(design.get(key)) is not int or int(design[key]) < 0:
            raise TimingError("frozen integer is malformed: " + key)
    if int(design["evict_bytes"]) < 268435456:
        raise TimingError("frozen eviction allocation is too small")
    for key in ("head", "tree"):
        if not isinstance(design.get(key), str) or \
                re.fullmatch(r"[0-9a-f]{40}", str(design[key])) is None:
            raise TimingError("frozen Git identity is malformed: " + key)
    for key in (
            "tasks_manifest_sha256", "prelaunch_receipt_sha256",
            "recovery_receipt_sha256", "pinned_build_sha256",
            "binary_sha256"):
        if not isinstance(design.get(key), str) or \
                SHA256_RE.fullmatch(str(design[key])) is None:
            raise TimingError("frozen digest is malformed: " + key)
    for key, core_key in (("timing_topology", "core"),
                          ("controller_topology", "controller_core"),
                          ("thermal_topology", "thermal_core")):
        value = design.get(key)
        if not isinstance(value, dict) or value.get("core") != design[core_key] or \
                value.get("numa_node") != design["numa_node"]:
            raise TimingError("frozen topology is malformed: " + key)
    quiet_cpus = design.get("timing_quiet_cpus")
    controller_siblings = design["controller_topology"].get("siblings")
    thermal_siblings = design["thermal_topology"].get("siblings")
    if not isinstance(quiet_cpus, list) or not quiet_cpus or \
            any(type(cpu) is not int or cpu < 0 for cpu in quiet_cpus) or \
            quiet_cpus != sorted(set(quiet_cpus)) or \
            quiet_cpus != design["timing_topology"].get("llc_shared_cpus") or \
            not isinstance(controller_siblings, list) or \
            not isinstance(thermal_siblings, list) or \
            design["core"] not in quiet_cpus or \
            design["controller_core"] in quiet_cpus or \
            design["thermal_core"] in quiet_cpus or \
            design["controller_core"] == design["thermal_core"] or \
            design["controller_core"] in thermal_siblings or \
            design["thermal_core"] in controller_siblings:
        raise TimingError("frozen timing LLC quiet CPU set is malformed")
    if design.get("python") != common.python_runtime_identity():
        raise TimingError("frozen Python runtime changed")
    tools = design.get("tools")
    required = {"git", "cmake", "taskset", "numactl", "sudo", "fuser",
                "timeout", "python3", "env", "true"}
    if not isinstance(tools, dict) or set(tools) != required:
        raise TimingError("frozen tool ledger changed")
    return design


def load_tasks(root: Path, design: Mapping[str, object]) -> List[Dict[str, object]]:
    path = root / "tasks_manifest.jsonl"
    if sha256_file(path) != design.get("tasks_manifest_sha256"):
        raise TimingError("task manifest hash mismatch")
    tasks: List[Dict[str, object]] = []
    for line_number, line in enumerate(path.read_bytes().splitlines(keepends=True), 1):
        try:
            value = json.loads(line.decode("ascii"))
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise TimingError("task manifest line is malformed") from exc
        if not isinstance(value, dict) or canonical_json(value) != line or \
                value.get("job") != len(tasks):
            raise TimingError("task manifest line %d is noncanonical" % line_number)
        tasks.append(value)
    if tasks != generate_tasks() or len(tasks) != design.get("task_count"):
        raise TimingError("task manifest Cartesian product changed")
    return tasks


def verify_immutable_inputs(root: Path, design: Mapping[str, object]) -> None:
    common = build_module()
    expected = design.get("immutable_files")
    if not isinstance(expected, dict) or not expected:
        raise TimingError("immutable file ledger is malformed")
    actual = set()
    for dirname in ("frozen", "provenance", "external"):
        directory = root / dirname
        if directory.is_symlink() or not directory.is_dir():
            raise TimingError("immutable directory is unsafe: " + dirname)
        for path in directory.iterdir():
            if path.is_symlink() or not path.is_file():
                raise TimingError("immutable directory contains a non-file")
            actual.add(str(path.relative_to(root)))
    if actual != set(expected):
        raise TimingError("immutable input inventory changed")
    for relative, digest in expected.items():
        if not isinstance(digest, str) or SHA256_RE.fullmatch(digest) is None or \
                sha256_file(root / relative) != digest:
            raise TimingError("immutable input changed: " + str(relative))
    tools = design.get("tools")
    if not isinstance(tools, dict):
        raise TimingError("runtime tool ledger is missing")
    for name, receipt in tools.items():
        if not isinstance(receipt, dict) or set(receipt) != {"path", "sha256"}:
            raise TimingError("runtime tool receipt is malformed: " + str(name))
        path = Path(str(receipt.get("path")))
        digest = receipt.get("sha256")
        if path.is_symlink() or not path.is_file() or not os.access(path, os.X_OK) or \
                not isinstance(digest, str) or SHA256_RE.fullmatch(digest) is None or \
                sha256_file(path) != digest:
            raise TimingError("runtime tool changed: " + str(name))
    if design.get("binary_sha256") != expected.get("frozen/wirehair_v2_bench") or \
            design.get("pinned_build_sha256") != expected.get(
                "frozen/pinned_build.json") or \
            design.get("recovery_receipt_sha256") != expected.get(
                "external/" + RECOVERY_RESULT_NAME):
        raise TimingError("design-to-immutable binding changed")
    record = json.loads((root / "frozen/pinned_build.json").read_text(
        encoding="ascii"))
    common.validate_pinned_build_record(record, root / "frozen/CMakeCache.txt")
    common.validate_pinned_build_contract_binding(
        {"build_policy": record["build_policy"],
         "build_command": record["build_command"],
         "source_commit": design["head"],
         "binary_sha256": design["binary_sha256"]},
        record, "segmented timing design")
    recovery = _validate_recovery_result(json.loads(
        (root / "external" / RECOVERY_RESULT_NAME).read_text(encoding="utf-8")))
    if recovery["evidence"]["validated_summary_sha256"] != \
            design["recovery_validated_summary_sha256"] or \
            recovery["evidence"]["data_manifest_sha256"] != \
            design["recovery_data_manifest_sha256"]:
        raise TimingError("recovery receipt binding changed")


def verify_active_runner(root: Path, design: Mapping[str, object]) -> None:
    script = Path(__file__).resolve()
    expected = root / "frozen" / CONTROLLER_NAME
    if script != expected or sha256_file(script) != design["immutable_files"].get(
            "frozen/" + CONTROLLER_NAME):
        raise TimingError("run/verify must use the exact frozen controller")


def validate_prepare(root: Path, design: Mapping[str, object]) -> None:
    prepare = verify_sealed(
        load_canonical(root / "prepare_receipt.json", "prepare receipt"),
        SCHEMA + ".prepare", "prepare receipt")
    prelaunch = verify_sealed(
        load_canonical(root / "prelaunch_receipt.json", "prelaunch receipt"),
        SCHEMA + ".prelaunch", "prelaunch receipt")
    if set(prepare) != {
            "schema", "self_sha256_excluding_field", "prepared_utc", "head",
            "tree", "design_sha256", "tasks_manifest_sha256",
            "prelaunch_receipt_sha256", "binary_sha256", "runner_sha256",
            "isolation_helper_sha256", "thermal_sampler_sha256"}:
        raise TimingError("prepare receipt field set changed")
    if set(prelaunch) != {
            "schema", "self_sha256_excluding_field", "created_utc",
            "boundary_smoke", "quiet_timing_launched",
            "recovery_receipt_sha256", "recovery_validated_summary_sha256",
            "recovery_data_manifest_sha256"}:
        raise TimingError("prelaunch receipt field set changed")
    validate_utc_timestamp(prepare.get("prepared_utc"), "prepare completion")
    validate_utc_timestamp(prelaunch.get("created_utc"), "prelaunch creation")
    if str(prelaunch["created_utc"]) > str(prepare["prepared_utc"]):
        raise TimingError("prepare receipt predates its prelaunch evidence")
    exact = {
        "head": design["head"], "tree": design["tree"],
        "design_sha256": sha256_file(root / "design.json"),
        "tasks_manifest_sha256": sha256_file(root / "tasks_manifest.jsonl"),
        "prelaunch_receipt_sha256": sha256_file(root / "prelaunch_receipt.json"),
        "binary_sha256": design["binary_sha256"],
        "runner_sha256": design["immutable_files"]["frozen/" + CONTROLLER_NAME],
        "isolation_helper_sha256": design["immutable_files"][
            "frozen/" + ISOLATION_HELPER_NAME],
        "thermal_sampler_sha256": design["immutable_files"][
            "frozen/" + SAMPLER_NAME],
    }
    for key, value in exact.items():
        if prepare.get(key) != value:
            raise TimingError("prepare receipt mismatch: " + key)
    prelaunch_exact = {
        "quiet_timing_launched": False,
        "recovery_receipt_sha256": design["recovery_receipt_sha256"],
        "recovery_validated_summary_sha256":
            design["recovery_validated_summary_sha256"],
        "recovery_data_manifest_sha256":
            design["recovery_data_manifest_sha256"],
    }
    if any(prelaunch.get(key) != value
           for key, value in prelaunch_exact.items()):
        raise TimingError("prelaunch recovery/launch binding changed")
    smoke = prelaunch.get("boundary_smoke")
    if not isinstance(smoke, dict) or set(smoke) != {
            "timing_evidence", "scope", "records"} or \
            smoke.get("timing_evidence") is not False or \
            smoke.get("scope") != \
                "bounded CLI/architecture smoke under ordinary machine load" or \
            not isinstance(smoke.get("records"), list) or \
            len(smoke["records"]) != 4:
        raise TimingError("prelaunch smoke receipt changed")
    expected_specs = (
        (3200, 64, "three-048"), (3200, 64, "four-0369"),
        (64000, 1280, "three-048"), (64000, 1280, "four-0369"),
    )
    for index, (record, spec) in enumerate(zip(smoke["records"], expected_specs)):
        if not isinstance(record, dict) or set(record) != {
                "task", "executed_argv", "replay_argv",
                "executed_binary_sha256", "stdout_name", "stdout_sha256",
                "trace_sha256", "cell_class", "nonpromotional_contaminations",
                "started_monotonic_ns", "ended_monotonic_ns", "duration_ns"}:
            raise TimingError("prelaunch smoke record field set changed")
        K, bb, candidate = spec
        task = {
            "K": K, "bb": bb, "seed_index": 0, "seed": SEEDS[0],
            "schedule": "burst", "candidate": candidate,
            "cache_state": "warm",
        }
        replay_argv = command_for(
            design, task, cycle_index=2, evict_bytes=4096)
        executed_argv = record.get("executed_argv")
        if not isinstance(executed_argv, list) or len(executed_argv) != \
                len(replay_argv) or len(executed_argv) <= 6 or \
                not isinstance(executed_argv[6], str):
            raise TimingError("prelaunch smoke executed command is malformed")
        executed_binary = Path(executed_argv[6])
        if record.get("task") != task or \
                record.get("replay_argv") != replay_argv or \
                executed_argv != command_for(
                    design, task, cycle_index=2, binary=executed_binary,
                    evict_bytes=4096) or \
                not executed_binary.is_absolute() or \
                executed_binary.name != "wirehair_v2_bench" or \
                executed_binary.parent.name != "frozen" or \
                record.get("executed_binary_sha256") != design["binary_sha256"]:
            raise TimingError("prelaunch smoke command changed")
        name = record.get("stdout_name")
        expected_name = "smoke.%d.%s.K%d.bb%d.csv" % (
            index, candidate, K, bb)
        if name != expected_name:
            raise TimingError("prelaunch smoke output name changed")
        raw = (root / "provenance" / expected_name).read_bytes()
        parsed = parse_grouped_output(
            raw, task, int(design["core"]), replacement_cycle=2,
            expected_evict_bytes=4096)
        start_ns = record.get("started_monotonic_ns")
        end_ns = record.get("ended_monotonic_ns")
        duration_ns = record.get("duration_ns")
        if type(start_ns) is not int or type(end_ns) is not int or \
                type(duration_ns) is not int or start_ns < 0 or end_ns <= start_ns or \
                duration_ns != end_ns - start_ns or not parsed.common_success or \
                record.get("stdout_sha256") != sha256_bytes(raw) or \
                record.get("trace_sha256") != parsed.trace_sha256 or \
                record.get("cell_class") != parsed.cell_class or \
                record.get("nonpromotional_contaminations") != \
                    list(parsed.contaminations):
            raise TimingError("prelaunch smoke output changed")


def _ancestor_pids() -> set[int]:
    result = {os.getpid()}
    pid = os.getppid()
    while pid > 1 and pid not in result:
        result.add(pid)
        try:
            raw = (Path("/proc") / str(pid) / "stat").read_bytes()
            right = raw.rfind(b")")
            fields = raw[right + 2:].split() if right >= 0 else []
            pid = int(fields[1]) if len(fields) >= 2 else 0
        except (OSError, ValueError):
            break
    return result


def conflicting_processes(*, excluded_pids: Iterable[int] = (),
                          ) -> List[Dict[str, object]]:
    ignored = _ancestor_pids() | {int(pid) for pid in excluded_pids}
    signatures = (
        b"wirehair_v2_bench", b"wirehair_load_fillers.sh",
        b"while :; do :; done", b"wirehair_expo_thermal_sampler.py",
        b"wh2_cross_payload", b"cross-payload-recovery",
        b"wh2_p32_dispatch_timing.py", b"wh2_degree",
        b"wh2_segmented_anchor_solve_timing.py",
    )
    result = []
    for path in Path("/proc").glob("[0-9]*/cmdline"):
        pid = int(path.parent.name)
        if pid in ignored:
            continue
        try:
            raw = path.read_bytes()
        except OSError:
            continue
        if any(signature in raw for signature in signatures):
            result.append({
                "pid": pid,
                "cmdline_sha256": sha256_bytes(raw),
                "cmdline_preview": raw.replace(b"\0", b" ")[:240].decode(
                    "utf-8", errors="replace"),
            })
    return sorted(result, key=lambda item: int(item["pid"]))


def cpu_ticks(cpu: int) -> Tuple[int, ...]:
    prefix = "cpu%d " % cpu
    for line in Path("/proc/stat").read_text(encoding="ascii").splitlines():
        if line.startswith(prefix):
            fields = line.split()[1:]
            if not fields or any(not field.isdigit() for field in fields):
                break
            return tuple(int(field) for field in fields)
    raise TimingError("CPU tick counter is unavailable")


def busy_ticks(before: Sequence[int], after: Sequence[int]) -> int:
    if len(before) != len(after) or len(before) < 5:
        raise TimingError("CPU tick vectors are incompatible")
    if any(type(value) is not int or value < 0
           for value in tuple(before) + tuple(after)) or \
            any(new < old for old, new in zip(before, after)):
        raise TimingError("CPU tick counter regressed or is malformed")
    idle_before = before[3] + before[4]
    idle_after = after[3] + after[4]
    # Linux reports guest and guest_nice after steal, but those two counters
    # are already included in user and nice.  Exclude them from the total so
    # a guest workload cannot be charged twice.
    accounted = min(len(before), 8)
    result = (sum(after[:accounted]) - sum(before[:accounted]) -
              (idle_after - idle_before))
    if result < 0:
        raise TimingError("CPU busy tick delta is negative")
    return result


def cpu_tick_snapshot(cpus: Sequence[int]) -> Dict[str, Tuple[int, ...]]:
    cpu_set = list(cpus)
    if not cpu_set or any(type(cpu) is not int or cpu < 0 for cpu in cpu_set) or \
            cpu_set != sorted(set(cpu_set)):
        raise TimingError("CPU tick snapshot set is malformed")
    return {str(cpu): cpu_ticks(cpu) for cpu in cpu_set}


def cpu_busy_tick_deltas(
        before: Mapping[str, Sequence[int]],
        after: Mapping[str, Sequence[int]]) -> Dict[str, int]:
    if not before or set(before) != set(after) or \
            any(not isinstance(cpu, str) or UINT_RE.fullmatch(cpu) is None
                for cpu in before):
        raise TimingError("CPU tick snapshots are incompatible")
    return {cpu: busy_ticks(before[cpu], after[cpu])
            for cpu in sorted(before, key=int)}


def measure_quiet_cpus(cpus: Sequence[int]) -> Dict[str, object]:
    cpu_set = list(cpus)
    if not cpu_set or any(type(cpu) is not int or cpu < 0 for cpu in cpu_set) or \
            cpu_set != sorted(set(cpu_set)):
        raise TimingError("quiet CPU set is malformed")
    before = cpu_tick_snapshot(cpu_set)
    started = time.monotonic()
    time.sleep(QUIET_SAMPLE_SECONDS)
    after = cpu_tick_snapshot(cpu_set)
    observed = time.monotonic() - started
    if not math.isfinite(observed) or observed < QUIET_SAMPLE_SECONDS:
        raise TimingError("quiet CPU sample interval was truncated")
    values = cpu_busy_tick_deltas(before, after)
    return {
        "measured": True,
        "cpu_set": cpu_set,
        "requested_duration_s": QUIET_SAMPLE_SECONDS,
        "observed_duration_s": observed,
        "busy_ticks": values,
        "maximum_allowed_per_cpu": MAX_QUIET_BUSY_TICKS_PER_CPU,
    }


def quiet_cpu_blocker(receipt: Mapping[str, object]) -> Optional[Dict[str, object]]:
    values = receipt.get("busy_ticks")
    if not isinstance(values, dict):
        raise TimingError("quiet CPU receipt lacks tick deltas")
    busy = {cpu: ticks for cpu, ticks in values.items()
            if type(ticks) is int and ticks > MAX_QUIET_BUSY_TICKS_PER_CPU}
    malformed = {cpu: ticks for cpu, ticks in values.items()
                 if type(ticks) is not int or ticks < 0}
    if malformed:
        raise TimingError("quiet CPU receipt has malformed tick deltas")
    if not busy:
        return None
    return {
        "name": "timing-llc-not-quiet",
        "cpu_set": receipt.get("cpu_set"),
        "busy_ticks": busy,
        "maximum_allowed_per_cpu": MAX_QUIET_BUSY_TICKS_PER_CPU,
    }


def preflight_report(root: Path, design: Mapping[str, object],
                     *, measure_quiet: bool) -> Dict[str, object]:
    isolation = isolation_module()
    blockers: List[Dict[str, object]] = []
    conflicts = conflicting_processes()
    if conflicts:
        blockers.append({"name": "conflicting-processes", "processes": conflicts})
    tools = design["tools"]
    try:
        readers = isolation.sole_i2c_readers(
            Path(tools["fuser"]["path"]), Path(tools["sudo"]["path"]),
            Path(tools["timeout"]["path"]))
    except (OSError, ValueError, RuntimeError) as exc:
        readers = ()
        blockers.append({"name": "i2c-reader-inspection-failed",
                         "detail": str(exc)})
    if readers:
        blockers.append({"name": "pre-existing-i2c-readers", "pids": list(readers)})
    topology = {}
    try:
        for core_key, receipt_key in (("core", "timing_topology"),
                                      ("controller_core", "controller_topology"),
                                      ("thermal_core", "thermal_topology")):
            current = isolation.topology_record(
                int(design[core_key]), int(design["numa_node"]))
            topology[core_key] = current
            if current != design[receipt_key]:
                blockers.append({"name": "topology-or-power-policy-changed",
                                 "core": design[core_key]})
    except (OSError, ValueError, RuntimeError) as exc:
        blockers.append({"name": "topology-unavailable", "detail": str(exc)})
    try:
        sudo = subprocess.run(
            (tools["sudo"]["path"], "-n", tools["true"]["path"]),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=10.0,
            check=False, env={"PATH": "/usr/bin:/bin", "LC_ALL": "C",
                              "LANG": "C", "TZ": "UTC"})
    except (OSError, subprocess.TimeoutExpired):
        sudo = None
    if sudo is None or sudo.returncode != 0 or sudo.stdout or sudo.stderr:
        blockers.append({"name": "passwordless-sudo-unavailable"})
    quiet_receipt: Dict[str, object] = {
        "measured": False,
        "cpu_set": design["timing_quiet_cpus"],
    }
    if measure_quiet and not conflicts:
        try:
            quiet_receipt = measure_quiet_cpus(design["timing_quiet_cpus"])
            blocker = quiet_cpu_blocker(quiet_receipt)
            if blocker is not None:
                blockers.append(blocker)
        except (OSError, ValueError, TimingError) as exc:
            quiet_receipt = {
                "measured": False,
                "cpu_set": design["timing_quiet_cpus"],
                "error": str(exc),
            }
            blockers.append({"name": "timing-llc-quiet-check-failed",
                             "detail": str(exc)})
    return {
        "ready": not blockers, "checked_utc": utc_now(),
        "blockers": blockers, "topology": topology,
        "quiet_llc_receipt": quiet_receipt,
        "owned_sampler_temperature_gate_pending": not blockers,
        "requirements": design["launch_prerequisites"],
    }


def preflight_command(args: argparse.Namespace) -> None:
    root = Path(args.result_dir).resolve()
    design = load_design(root)
    verify_root_layout(root)
    verify_active_runner(root, design)
    load_tasks(root, design)
    verify_immutable_inputs(root, design)
    validate_prepare(root, design)
    print(json.dumps(preflight_report(root, design, measure_quiet=True),
                     sort_keys=True))


def verify_root_layout(root: Path) -> None:
    required_directories = {
        ".transactions", "external", "frozen", "interrupted", "ledger",
        "provenance", "runtime-home", "segments",
    }
    required_files = {
        "controller.lock", "design.json", "prelaunch_receipt.json",
        "prepare_receipt.json", "tasks_manifest.jsonl",
    }
    optional_files = {"campaign_receipt.json", "validated_summary.json"}
    if root.is_symlink() or not root.is_dir():
        raise TimingError("campaign root is missing or unsafe")
    actual = {path.name for path in root.iterdir()}
    if not required_directories | required_files <= actual or \
            not actual <= required_directories | required_files | optional_files:
        raise TimingError("campaign root layout contains missing or unknown entries")
    for name in required_directories:
        path = root / name
        if path.is_symlink() or not path.is_dir():
            raise TimingError("campaign directory is missing or unsafe: " + name)
    for name in required_files | (optional_files & actual):
        path = root / name
        if path.is_symlink() or not path.is_file():
            raise TimingError("campaign file is missing or unsafe: " + name)
    if any((root / "runtime-home").iterdir()):
        raise TimingError("runtime HOME contains unreceipted state")


def _archive_file_receipts(directory: Path) -> Dict[str, Dict[str, object]]:
    if directory.is_symlink() or not directory.is_dir():
        raise TimingError("interrupted archive is unsafe")
    files: Dict[str, Dict[str, object]] = {}
    total = 0
    for path in sorted(directory.iterdir()):
        if path.is_symlink() or not path.is_file() or not path.name or \
                Path(path.name).name != path.name:
            raise TimingError("interrupted archive contains an unsafe entry")
        stat = path.stat()
        total += stat.st_size
        if total > MAX_INTERRUPTED_ARCHIVE_BYTES:
            raise TimingError("interrupted archive exceeds its byte bound")
        files[path.name] = {"size": stat.st_size, "sha256": sha256_file(path)}
    return files


def _archive_directory(root: Path, source: Path, stem: str,
                       *, reason: Optional[str] = None) -> Dict[str, object]:
    files = _archive_file_receipts(source)
    destination = root / "interrupted" / (
        stem + ".recovered.%d" % time.monotonic_ns())
    if destination.exists() or destination.is_symlink():
        raise TimingError("interrupted archive destination collision")
    os.chmod(source, 0o700)
    os.replace(source, destination)
    os.chmod(destination, 0o555)
    fsync_directory(source.parent)
    fsync_directory(destination.parent)
    record: Dict[str, object] = {
        "source": stem, "archive": destination.name, "files": files,
    }
    if reason is not None:
        record["reason"] = reason
    return record


def recover_interrupted_transactions(root: Path) -> List[Dict[str, object]]:
    entries = sorted((root / ".transactions").iterdir())
    if len(entries) > TASK_COUNT:
        raise TimingError("transaction staging count exceeds the task bound")
    recovered = []
    for entry in entries:
        if entry.is_symlink() or not entry.is_dir() or \
                TXN_DIR_RE.fullmatch(entry.name) is None:
            raise TimingError("unsafe interrupted transaction entry")
        recovered.append(_archive_directory(root, entry, entry.name))
    return recovered


def begin_task_transaction(root: Path, task: Mapping[str, object]) -> Path:
    job = int(task["job"])
    staging = root / ".transactions" / ("%04d.part.%d" % (job, os.getpid()))
    final = root / "ledger" / ("%04d" % job)
    if staging.exists() or staging.is_symlink() or final.exists() or final.is_symlink():
        raise TimingError("task transaction already exists")
    staging.mkdir(mode=0o700)
    atomic_write(staging / "intent.json", canonical_json(sealed_record(
        SCHEMA + ".task_intent", {
            "created_utc": utc_now(), "job": job,
            "task_id": task["task_id"],
            "task_sha256": sha256_bytes(canonical_json(task)),
        })))
    fsync_directory(staging.parent)
    return staging


def commit_task_transaction(root: Path, task: Mapping[str, object], raw: bytes,
                            stderr: bytes, receipt: Mapping[str, object],
                            staging: Path) -> Path:
    job = int(task["job"])
    final = root / "ledger" / ("%04d" % job)
    if final.exists() or final.is_symlink() or staging.parent != root / ".transactions" or \
            staging.is_symlink() or not staging.is_dir() or \
            TXN_DIR_RE.fullmatch(staging.name) is None:
        raise TimingError("task transaction identity is invalid")
    try:
        atomic_write(staging / "stdout.csv", raw)
        atomic_write(staging / "stderr.bin", stderr)
        atomic_write(staging / "receipt.json", canonical_json(receipt))
        os.replace(staging, final)
        os.chmod(final, 0o555)
        fsync_directory(final.parent)
    except BaseException:
        if staging.exists() and not staging.is_symlink():
            os.chmod(staging, 0o700)
        raise
    return final


def _task_execution(root: Path, design: Mapping[str, object],
                    task: Mapping[str, object], timeout_s: float,
                    poll_callback) -> Dict[str, object]:
    isolation = isolation_module()
    environment = sanitized_environment(root / "runtime-home")
    command = command_for(design, task)
    staging = begin_task_transaction(root, task)
    prior: List[Dict[str, object]] = []
    timing_core = int(design["core"])
    peer_cpus = [int(cpu) for cpu in design["timing_quiet_cpus"]
                 if int(cpu) != timing_core]
    if not peer_cpus:
        raise TimingError("timing LLC has no independently observable peer CPU")
    for attempt in range(MAX_ENVIRONMENTAL_ATTEMPTS):
        started_utc = utc_now()
        start_ns = time.monotonic_ns()
        peer_ticks_before = cpu_tick_snapshot(peer_cpus)
        rc, stdout, stderr = isolation.run_bounded(
            command, environment, timeout_s, MAX_STDOUT_BYTES,
            MAX_STDERR_BYTES, poll_callback=poll_callback)
        peer_ticks_after = cpu_tick_snapshot(peer_cpus)
        end_ns = time.monotonic_ns()
        peer_busy = cpu_busy_tick_deltas(peer_ticks_before, peer_ticks_after)
        if rc != 0 or stderr:
            raise TimingError("substantive task child failure exit=%d stderr=%r" %
                              (rc, stderr[:1000]))
        parsed = parse_grouped_output(
            stdout, task, int(design["core"]),
            expected_evict_bytes=int(design["evict_bytes"]))
        contaminations = list(parsed.contaminations)
        contaminations.extend(
            "llc-peer-cpu%s-busy-ticks:%d" % (cpu, ticks)
            for cpu, ticks in peer_busy.items()
            if ticks > MAX_QUIET_BUSY_TICKS_PER_CPU)
        if contaminations:
            prior.append({
                "attempt": attempt, "started_utc": started_utc,
                "start_monotonic_ns": start_ns, "end_monotonic_ns": end_ns,
                "stdout_sha256": sha256_bytes(stdout),
                "stderr_sha256": sha256_bytes(stderr),
                "llc_peer_busy_ticks": peer_busy,
                "contaminations": sorted(set(contaminations)),
            })
            if attempt + 1 == MAX_ENVIRONMENTAL_ATTEMPTS:
                raise ContaminationError("environmental retry limit exhausted")
            continue
        receipt = sealed_record(SCHEMA + ".task", {
            "job": task["job"], "task_id": task["task_id"],
            "task_sha256": sha256_bytes(canonical_json(task)),
            "argv": command, "child_environment": environment,
            "attempt": attempt, "prior_contaminations": prior,
            "llc_peer_busy_ticks": peer_busy,
            "started_utc": started_utc, "start_monotonic_ns": start_ns,
            "end_monotonic_ns": end_ns, "duration_ns": end_ns - start_ns,
            "stderr_sha256": sha256_bytes(stderr), **task_summary(parsed),
        })
        final = commit_task_transaction(
            root, task, stdout, stderr, receipt, staging)
        return {
            "job": task["job"],
            "receipt_sha256": sha256_file(final / "receipt.json"),
            "retry_count": attempt,
        }
    raise TimingError("unreachable environmental retry state")


def validate_task_entry(root: Path, design: Mapping[str, object],
                        task: Mapping[str, object]) -> Tuple[Dict[str, object], ParsedOutput]:
    directory = root / "ledger" / ("%04d" % int(task["job"]))
    if directory.is_symlink() or not directory.is_dir():
        raise TimingError("task ledger directory is missing or unsafe")
    if sorted(path.name for path in directory.iterdir()) != \
            ["intent.json", "receipt.json", "stderr.bin", "stdout.csv"] or \
            any(path.is_symlink() or not path.is_file()
                for path in directory.iterdir()):
        raise TimingError("task ledger file set changed")
    stderr = (directory / "stderr.bin").read_bytes()
    if stderr:
        raise TimingError("successful task has nonempty stderr")
    raw = (directory / "stdout.csv").read_bytes()
    parsed = parse_grouped_output(
        raw, task, int(design["core"]),
        expected_evict_bytes=int(design["evict_bytes"]))
    intent = verify_sealed(
        load_canonical(directory / "intent.json", "task intent"),
        SCHEMA + ".task_intent", "task intent")
    if set(intent) != {
            "schema", "self_sha256_excluding_field", "created_utc", "job",
            "task_id", "task_sha256"}:
        raise TimingError("task intent field set changed")
    validate_utc_timestamp(intent.get("created_utc"), "task intent creation")
    exact_identity = {
        "job": task["job"], "task_id": task["task_id"],
        "task_sha256": sha256_bytes(canonical_json(task)),
    }
    if any(intent.get(key) != value for key, value in exact_identity.items()):
        raise TimingError("task intent identity changed")
    receipt = verify_sealed(
        load_canonical(directory / "receipt.json", "task receipt"),
        SCHEMA + ".task", "task receipt")
    fields = {
        "schema", "self_sha256_excluding_field", "job", "task_id",
        "task_sha256", "argv", "child_environment", "attempt",
        "prior_contaminations", "llc_peer_busy_ticks", "started_utc",
        "start_monotonic_ns",
        "end_monotonic_ns", "duration_ns", "stderr_sha256",
        "stdout_sha256", "trace_sha256", "cell_class", "common_success",
        "control_ns", "candidate_ns", "control_work_sha256",
        "candidate_work_sha256", "control_attempt", "control_matrix_seed",
        "control_peel_seed", "candidate_attempt", "candidate_matrix_seed",
        "candidate_peel_seed",
    }
    if set(receipt) != fields:
        raise TimingError("task receipt field set changed")
    expected = {
        **exact_identity, "stderr_sha256": sha256_bytes(stderr),
        **task_summary(parsed),
    }
    for key, value in expected.items():
        if receipt.get(key) != value:
            raise TimingError("task receipt mismatch: " + key)
    if receipt.get("argv") != command_for(design, task) or \
            receipt.get("child_environment") != sanitized_environment(
                root / "runtime-home"):
        raise TimingError("task command/environment receipt changed")
    expected_peer_cpus = {
        str(cpu) for cpu in design["timing_quiet_cpus"]
        if int(cpu) != int(design["core"])}
    peer_busy = receipt.get("llc_peer_busy_ticks")
    if not expected_peer_cpus or not isinstance(peer_busy, dict) or \
            set(peer_busy) != expected_peer_cpus or \
            any(type(ticks) is not int or ticks < 0 or
                ticks > MAX_QUIET_BUSY_TICKS_PER_CPU
                for ticks in peer_busy.values()):
        raise TimingError("task timing-LLC peer receipt is malformed")
    attempt = receipt.get("attempt")
    prior = receipt.get("prior_contaminations")
    if type(attempt) is not int or not 0 <= attempt < MAX_ENVIRONMENTAL_ATTEMPTS or \
            not isinstance(prior, list) or len(prior) != attempt:
        raise TimingError("task retry ledger is malformed")
    started_ns = receipt.get("start_monotonic_ns")
    ended_ns = receipt.get("end_monotonic_ns")
    duration_ns = receipt.get("duration_ns")
    if type(started_ns) is not int or type(ended_ns) is not int or \
            type(duration_ns) is not int or started_ns < 0 or \
            ended_ns <= started_ns or duration_ns != ended_ns - started_ns:
        raise TimingError("task timing receipt is malformed")
    validate_utc_timestamp(receipt.get("started_utc"), "task start")
    if parsed.contaminations:
        raise TimingError("terminal task contains timing contamination")
    previous_end = -1
    for index, record in enumerate(prior):
        if not isinstance(record, dict) or set(record) != {
                "attempt", "started_utc", "start_monotonic_ns",
                "end_monotonic_ns", "stdout_sha256", "stderr_sha256",
                "llc_peer_busy_ticks", "contaminations"} or \
                record.get("attempt") != index:
            raise TimingError("prior contamination receipt is malformed")
        prior_start = record.get("start_monotonic_ns")
        prior_end = record.get("end_monotonic_ns")
        values = record.get("contaminations")
        prior_peer_busy = record.get("llc_peer_busy_ticks")
        if type(prior_start) is not int or type(prior_end) is not int or \
                prior_start < 0 or prior_end <= prior_start or \
                prior_start < previous_end or prior_end > started_ns or \
                SHA256_RE.fullmatch(str(record.get("stdout_sha256", ""))) is None or \
                record.get("stderr_sha256") != sha256_bytes(b"") or \
                not isinstance(prior_peer_busy, dict) or \
                set(prior_peer_busy) != expected_peer_cpus or \
                any(type(ticks) is not int or ticks < 0
                    for ticks in prior_peer_busy.values()) or \
                not isinstance(values, list) or not values or \
                values != sorted(set(values)) or \
                any(not isinstance(value, str) or not value for value in values):
            raise TimingError("prior contamination evidence is malformed")
        expected_peer_contaminations = {
            "llc-peer-cpu%s-busy-ticks:%d" % (cpu, ticks)
            for cpu, ticks in prior_peer_busy.items()
            if ticks > MAX_QUIET_BUSY_TICKS_PER_CPU}
        actual_peer_contaminations = {
            value for value in values
            if value.startswith("llc-peer-cpu")}
        if actual_peer_contaminations != expected_peer_contaminations:
            raise TimingError("prior timing-LLC contamination receipt changed")
        validate_utc_timestamp(record.get("started_utc"),
                               "prior contamination start")
        previous_end = prior_end
    return receipt, parsed


def existing_task_entries(root: Path, design: Mapping[str, object],
                          tasks: Sequence[Mapping[str, object]],
                          ) -> Dict[int, Dict[str, object]]:
    entries = sorted((root / "ledger").iterdir())
    if len(entries) > len(tasks):
        raise TimingError("task ledger exceeds manifest")
    result: Dict[int, Dict[str, object]] = {}
    for entry in entries:
        if entry.is_symlink() or not entry.is_dir() or \
                TASK_DIR_RE.fullmatch(entry.name) is None:
            raise TimingError("task ledger contains an unknown entry")
        job = int(entry.name)
        if job >= len(tasks):
            raise TimingError("task ledger job is outside the manifest")
        receipt, parsed = validate_task_entry(root, design, tasks[job])
        result[job] = {
            "receipt": receipt, "parsed": parsed,
            "receipt_sha256": sha256_file(entry / "receipt.json"),
        }
    return result


def _load_segment_terminal(path: Path) -> Dict[str, object]:
    return verify_sealed(
        load_canonical(path, "segment terminal"),
        SCHEMA + ".segment_terminal", "segment terminal")


def terminal_task_bindings(root: Path) -> Dict[int, str]:
    bindings: Dict[int, str] = {}
    terminals = sorted((root / "segments").glob("segment*.terminal.json"))
    for index, path in enumerate(terminals):
        if path.name != "segment%04d.terminal.json" % index:
            raise TimingError("segment terminal sequence has a gap")
        completed = _load_segment_terminal(path).get("completed_tasks")
        if not isinstance(completed, list):
            raise TimingError("segment completed-task binding is malformed")
        for item in completed:
            if not isinstance(item, dict) or type(item.get("job")) is not int or \
                    not isinstance(item.get("receipt_sha256"), str) or \
                    SHA256_RE.fullmatch(item["receipt_sha256"]) is None:
                raise TimingError("segment task binding is malformed")
            job = int(item["job"])
            if not 0 <= job < TASK_COUNT:
                raise TimingError("segment task binding is outside the manifest")
            if job in bindings:
                raise TimingError("task is bound by multiple thermal segments")
            bindings[job] = str(item["receipt_sha256"])
    return bindings


def recover_unbound_task_entries(root: Path,
                                 bindings: Mapping[int, str]) -> List[Dict[str, object]]:
    recovered = []
    for entry in sorted((root / "ledger").iterdir()):
        if entry.is_symlink() or not entry.is_dir() or \
                TASK_DIR_RE.fullmatch(entry.name) is None:
            raise TimingError("unsafe task ledger entry during resume")
        job = int(entry.name)
        if job >= TASK_COUNT:
            raise TimingError("unbound task ledger job is outside the manifest")
        if job in bindings:
            receipt = entry / "receipt.json"
            if receipt.is_symlink() or not receipt.is_file() or \
                    sha256_file(receipt) != bindings[job]:
                raise TimingError("terminal-bound task receipt changed")
            continue
        recovered.append(_archive_directory(
            root, entry, "task-%04d.unbound" % job,
            reason="no terminal thermal binding"))
    return recovered


def _referenced_interrupted_archives(root: Path) -> set[str]:
    referenced: set[str] = set()
    for path in sorted((root / "segments").glob("segment*.terminal.json")):
        records = _load_segment_terminal(path).get(
            "recovered_interrupted_transactions")
        if not isinstance(records, list):
            raise TimingError("segment interrupted archive ledger is malformed")
        for record in records:
            if not isinstance(record, dict) or not isinstance(
                    record.get("archive"), str):
                raise TimingError("segment interrupted archive receipt is malformed")
            name = str(record["archive"])
            if name in referenced:
                raise TimingError("interrupted archive is bound more than once")
            referenced.add(name)
    return referenced


def _archive_record(root: Path, directory: Path) -> Dict[str, object]:
    isolation = isolation_module()
    files = _archive_file_receipts(directory)
    thermal_match = isolation.THERMAL_ARCHIVE_FINAL_RE.fullmatch(directory.name)
    transaction_match = TRANSACTION_ARCHIVE_RE.fullmatch(directory.name)
    task_match = UNBOUND_TASK_ARCHIVE_RE.fullmatch(directory.name)
    if thermal_match is not None:
        failure = isolation._verify_thermal_failure_archive(directory)
        return {
            "source": "invalid-thermal-segment", "archive": directory.name,
            "files": files, "reason": str(failure["reason"]),
        }
    if transaction_match is not None:
        allowed = {
            "intent.json", "intent.json.part", "stdout.csv",
            "stdout.csv.part", "stderr.bin", "stderr.bin.part",
            "receipt.json", "receipt.json.part",
        }
        if not set(files) <= allowed:
            raise TimingError("interrupted task transaction file set changed")
        return {
            "source": transaction_match.group(1), "archive": directory.name,
            "files": files,
        }
    if task_match is not None:
        if set(files) != {
                "intent.json", "receipt.json", "stderr.bin", "stdout.csv"}:
            raise TimingError("unbound task archive file set changed")
        return {
            "source": task_match.group(1), "archive": directory.name,
            "files": files, "reason": "no terminal thermal binding",
        }
    raise TimingError("interrupted archive name is unknown")


def adopt_unbound_interrupted_archives(root: Path) -> List[Dict[str, object]]:
    referenced = _referenced_interrupted_archives(root)
    records = []
    for directory in sorted((root / "interrupted").iterdir()):
        if directory.name not in referenced:
            if directory.is_symlink() or not directory.is_dir():
                raise TimingError("unbound interrupted archive is unsafe")
            os.chmod(directory, 0o555)
            fsync_directory(directory.parent)
            records.append(_archive_record(root, directory))
    return records


def verify_interrupted_records(
    root: Path, records: Sequence[Mapping[str, object]],
    *, expected_thermal_segment: Optional[int] = None,
) -> None:
    isolation = isolation_module()
    seen: set[str] = set()
    for record in records:
        if not isinstance(record, dict) or not isinstance(record.get("archive"), str):
            raise TimingError("interrupted archive receipt is malformed")
        name = str(record["archive"])
        if not name or Path(name).name != name or name in seen:
            raise TimingError("interrupted archive identity is malformed")
        seen.add(name)
        match = isolation.THERMAL_ARCHIVE_FINAL_RE.fullmatch(name)
        if match is not None and expected_thermal_segment is not None and \
                int(match.group(1)) != expected_thermal_segment:
            raise TimingError("thermal archive is bound to the wrong segment")
        if record != _archive_record(root, root / "interrupted" / name):
            raise TimingError("interrupted archive receipt changed")


def _next_segment(root: Path) -> int:
    terminals = sorted((root / "segments").glob("segment*.terminal.json"))
    if len(terminals) >= 10000:
        raise TimingError("segment sequence exhausted its four-digit domain")
    for index, path in enumerate(terminals):
        if path.name != "segment%04d.terminal.json" % index:
            raise TimingError("segment terminal sequence has a gap")
    allowed = {
        *("segment%04d.terminal.json" % index for index in range(len(terminals))),
        *("segment%04d.thermal.csv" % index for index in range(len(terminals))),
    }
    unknown = [path.name for path in (root / "segments").iterdir()
               if path.name not in allowed]
    if unknown:
        raise TimingError("segments directory contains stale entry: " + unknown[0])
    return len(terminals)


def _sampler_public_cmdline(root: Path, design: Mapping[str, object],
                            segment: int) -> List[str]:
    command = isolation_module()._sampler_cmdline(
        root, design, segment, public=True)
    if not isinstance(command, list) or \
            any(not isinstance(value, str) or not value for value in command):
        raise TimingError("isolation helper returned a malformed sampler command")
    return command


def _validate_initial_thermal_gate(value: object) -> None:
    if not isinstance(value, dict) or set(value) != {
            "last_monotonic_s", "latest_age_s", "max_gap_s", "cpu_max_c",
            "dimm_max_c"}:
        raise TimingError("initial thermal gate receipt is malformed")
    for key in value:
        number = value[key]
        if type(number) not in (int, float) or not math.isfinite(float(number)):
            raise TimingError("initial thermal gate contains a nonfinite value")
    if float(value["latest_age_s"]) < 0.0 or \
            float(value["latest_age_s"]) > MAX_THERMAL_MARGIN_S or \
            float(value["max_gap_s"]) < 0.0 or \
            float(value["max_gap_s"]) > MAX_THERMAL_GAP_S or \
            float(value["cpu_max_c"]) >= MAX_CPU_TEMP_C or \
            float(value["dimm_max_c"]) >= MAX_DIMM_TEMP_C:
        raise TimingError("initial thermal gate exceeds the sealed limits")


def _validate_launch_preflight(value: object,
                               design: Mapping[str, object]) -> None:
    fields = {
        "ready", "checked_utc", "blockers", "topology",
        "quiet_llc_receipt", "owned_sampler_temperature_gate_pending",
        "requirements",
    }
    if not isinstance(value, dict) or set(value) != fields or \
            value.get("ready") is not True or value.get("blockers") != [] or \
            value.get("owned_sampler_temperature_gate_pending") is not True or \
            value.get("requirements") != design["launch_prerequisites"]:
        raise TimingError("segment launch preflight was not clean")
    validate_utc_timestamp(value.get("checked_utc"), "segment launch preflight")
    topology = value.get("topology")
    if topology != {
            "core": design["timing_topology"],
            "controller_core": design["controller_topology"],
            "thermal_core": design["thermal_topology"]}:
        raise TimingError("segment launch topology receipt changed")
    quiet = value.get("quiet_llc_receipt")
    expected_cpu_list = design["timing_quiet_cpus"]
    expected_cpus = {str(cpu) for cpu in expected_cpu_list}
    if not isinstance(quiet, dict) or set(quiet) != {
            "measured", "cpu_set", "requested_duration_s",
            "observed_duration_s", "busy_ticks",
            "maximum_allowed_per_cpu"} or \
            quiet.get("measured") is not True or \
            quiet.get("cpu_set") != expected_cpu_list or \
            quiet.get("requested_duration_s") != QUIET_SAMPLE_SECONDS or \
            type(quiet.get("observed_duration_s")) not in (int, float) or \
            not math.isfinite(float(quiet["observed_duration_s"])) or \
            float(quiet["observed_duration_s"]) < QUIET_SAMPLE_SECONDS or \
            quiet.get("maximum_allowed_per_cpu") != \
                MAX_QUIET_BUSY_TICKS_PER_CPU or \
            not isinstance(quiet.get("busy_ticks"), dict) or \
            set(quiet["busy_ticks"]) != expected_cpus or \
            any(type(ticks) is not int or ticks < 0 or
                ticks > MAX_QUIET_BUSY_TICKS_PER_CPU
                for ticks in quiet["busy_ticks"].values()):
        raise TimingError("segment timing-LLC quiet receipt changed")


def verify_segments(root: Path, design: Mapping[str, object],
                    entries: Mapping[int, Mapping[str, object]],
                    *, allow_unbound_archives: bool = False,
                    ) -> List[Dict[str, object]]:
    isolation = isolation_module()
    terminals = sorted((root / "segments").glob("segment*.terminal.json"))
    records: List[Dict[str, object]] = []
    completed_jobs: set[int] = set()
    referenced_archives: set[str] = set()
    saw_complete = False
    terminal_fields = {
        "schema", "self_sha256_excluding_field", "segment", "status",
        "started_utc", "ended_utc", "timing_start_monotonic_s",
        "benchmark_end_monotonic_s", "duration_s", "controller_pid",
        "controller_affinity", "design_sha256", "prepare_receipt_sha256",
        "prelaunch_receipt_sha256", "tasks_manifest_sha256",
        "launch_preflight", "initial_thermal_gate", "resumed_jobs_before",
        "completed_tasks", "recovered_interrupted_transactions", "failure",
        "sampler_identity_start", "sampler_identity_end", "graceful_stop",
        "graceful_stop_mechanism", "post_end_sample", "thermal_csv_name",
        "thermal_csv_sha256", "thermal_csv_size", "thermal_csv_device",
        "thermal_csv_inode", "thermal_csv_uid", "thermal_csv_mode",
        "thermal_csv_nlink", "thermal_summary",
    }
    for index, path in enumerate(terminals):
        if saw_complete:
            raise TimingError("a terminal segment follows campaign completion")
        if path.name != "segment%04d.terminal.json" % index:
            raise TimingError("terminal segment sequence changed")
        terminal = _load_segment_terminal(path)
        if set(terminal) != terminal_fields or terminal.get("segment") != index or \
                terminal.get("status") not in (
                    "failed", "bounded-incomplete", "complete"):
            raise TimingError("segment terminal field set/status changed")
        start = terminal.get("timing_start_monotonic_s")
        end = terminal.get("benchmark_end_monotonic_s")
        duration = terminal.get("duration_s")
        if type(start) not in (int, float) or type(end) not in (int, float) or \
                type(duration) not in (int, float) or \
                not all(math.isfinite(float(value))
                        for value in (start, end, duration)) or \
                float(start) < 0.0 or float(end) <= float(start) or \
                float(duration) != float(end) - float(start):
            raise TimingError("segment timing interval is malformed")
        validate_utc_timestamp(terminal.get("started_utc"), "segment start")
        validate_utc_timestamp(terminal.get("ended_utc"), "segment end")
        if type(terminal.get("controller_pid")) is not int or \
                int(terminal["controller_pid"]) <= 1 or \
                terminal.get("controller_affinity") != [design["controller_core"]]:
            raise TimingError("segment controller identity changed")
        for key, name in (
                ("design_sha256", "design.json"),
                ("prepare_receipt_sha256", "prepare_receipt.json"),
                ("prelaunch_receipt_sha256", "prelaunch_receipt.json"),
                ("tasks_manifest_sha256", "tasks_manifest.jsonl")):
            if terminal.get(key) != sha256_file(root / name):
                raise TimingError("segment input binding changed: " + key)
        launch = terminal.get("launch_preflight")
        _validate_launch_preflight(launch, design)
        _validate_initial_thermal_gate(terminal.get("initial_thermal_gate"))
        initial_sample_s = float(terminal["initial_thermal_gate"][
            "last_monotonic_s"])
        if abs(initial_sample_s - float(start)) > MAX_THERMAL_MARGIN_S:
            raise TimingError("initial thermal gate is not adjacent to timing start")
        start_identity = terminal.get("sampler_identity_start")
        end_identity = terminal.get("sampler_identity_end")
        identity_fields = {
            "pid", "start_tick", "process_group", "session_id", "boot_id",
            "cmdline", "cmdline_sha256", "uid", "affinity", "csv_device",
            "csv_inode",
        }
        if not isinstance(start_identity, dict) or start_identity != end_identity or \
                set(start_identity) != identity_fields or \
                start_identity.get("affinity") != [design["thermal_core"]] or \
                type(start_identity.get("pid")) is not int or \
                int(start_identity["pid"]) <= 1 or \
                type(start_identity.get("start_tick")) is not int or \
                int(start_identity["start_tick"]) <= 0 or \
                start_identity.get("process_group") != start_identity.get("session_id") or \
                type(start_identity.get("session_id")) is not int or \
                int(start_identity["session_id"]) <= 1:
            raise TimingError("sampler terminal identity is incomplete")
        if start_identity.get("uid") != 0 or \
                terminal.get("thermal_csv_uid") != 0 or \
                terminal.get("thermal_csv_mode") != 0o444 or \
                terminal.get("thermal_csv_nlink") != 1:
            raise TimingError("sampler terminal privilege/file seal changed")
        if start_identity.get("cmdline") != _sampler_public_cmdline(
                root, design, index) or \
                not isinstance(start_identity.get("cmdline_sha256"), str) or \
                SHA256_RE.fullmatch(start_identity["cmdline_sha256"]) is None or \
                start_identity.get("csv_device") != terminal.get(
                    "thermal_csv_device") or \
                start_identity.get("csv_inode") != terminal.get(
                    "thermal_csv_inode"):
            raise TimingError("sampler terminal process/CSV identity changed")
        if terminal.get("graceful_stop") is not True or \
                terminal.get("graceful_stop_mechanism") != \
                    "sudo-timeout-python-pidfd-identity-verified-sigterm" or \
                terminal.get("post_end_sample") is not True or \
                terminal.get("thermal_csv_name") != \
                    "segment%04d.thermal.csv" % index:
            raise TimingError("sampler lifecycle was not terminal and graceful")
        isolation.verify_terminal_thermal(
            root / "segments" / str(terminal["thermal_csv_name"]), terminal)
        resumed = terminal.get("resumed_jobs_before")
        completed = terminal.get("completed_tasks")
        if resumed != sorted(completed_jobs) or not isinstance(completed, list):
            raise TimingError("segment task ledger is malformed")
        pending = [job for job in range(int(design["task_count"]))
                   if job not in completed_jobs]
        completed_numbers = [item.get("job") if isinstance(item, dict) else None
                             for item in completed]
        if completed_numbers != pending[:len(completed_numbers)] or \
                len(completed_numbers) % int(design["cluster_size"]) != 0:
            raise TimingError("segment split or reordered a paired task cluster")
        for item in completed:
            if not isinstance(item, dict) or set(item) != {
                    "job", "receipt_sha256", "retry_count"} or \
                    type(item.get("job")) is not int or \
                    not isinstance(item.get("receipt_sha256"), str) or \
                    SHA256_RE.fullmatch(item["receipt_sha256"]) is None or \
                    type(item.get("retry_count")) is not int:
                raise TimingError("segment completed-task receipt is malformed")
            job = int(item["job"])
            if job in completed_jobs or job not in entries or \
                    item.get("receipt_sha256") != entries[job]["receipt_sha256"] or \
                    item.get("retry_count") != entries[job]["receipt"]["attempt"]:
                raise TimingError("segment completed-task binding changed")
            completed_jobs.add(job)
        status = terminal["status"]
        failure = terminal.get("failure")
        if status == "failed":
            if not isinstance(failure, dict) or set(failure) != {"class", "message"} or \
                    not all(isinstance(failure.get(key), str) and failure[key]
                            for key in ("class", "message")):
                raise TimingError("failed segment lacks an exact failure receipt")
        elif failure is not None:
            raise TimingError("successful segment carries a failure receipt")
        if status == "complete" and len(completed_jobs) != int(design["task_count"]):
            raise TimingError("complete segment did not finish the manifest")
        saw_complete = status == "complete"
        if status == "bounded-incomplete" and \
                (not completed or
                 len(completed_jobs) >= int(design["task_count"])):
            raise TimingError("bounded segment status contradicts the task ledger")
        recovered = terminal.get("recovered_interrupted_transactions")
        if not isinstance(recovered, list):
            raise TimingError("segment interrupted recovery ledger is malformed")
        verify_interrupted_records(
            root, recovered, expected_thermal_segment=index)
        for item in recovered:
            archive = str(item["archive"])
            if archive in referenced_archives:
                raise TimingError("interrupted archive is bound more than once")
            referenced_archives.add(archive)
        records.append(terminal)
    expected_segment_files = sorted(
        ["segment%04d.terminal.json" % index for index in range(len(terminals))] +
        ["segment%04d.thermal.csv" % index for index in range(len(terminals))])
    if sorted(path.name for path in (root / "segments").iterdir()) != \
            expected_segment_files:
        raise TimingError("segments directory contains ephemeral or unknown files")
    actual_archives = {path.name for path in (root / "interrupted").iterdir()}
    if allow_unbound_archives:
        if not referenced_archives <= actual_archives or \
                not completed_jobs <= set(entries):
            raise TimingError("bound segment/archive evidence is missing")
    elif actual_archives != referenced_archives or completed_jobs != set(entries):
        raise TimingError("task/archive ledgers are not exactly terminal-bound")
    if saw_complete != (len(completed_jobs) == int(design["task_count"])):
        raise TimingError("campaign completion lacks one exact complete terminal")
    return records


def _rollback_partial_cluster(
    root: Path, completed: List[Dict[str, object]],
) -> List[Dict[str, object]]:
    partial = len(completed) % CLUSTER_SIZE
    if partial == 0:
        return []
    records = []
    for item in completed[-partial:]:
        job = int(item["job"])
        entry = root / "ledger" / ("%04d" % job)
        records.append(_archive_directory(
            root, entry, "task-%04d.unbound" % job,
            reason="no terminal thermal binding"))
    del completed[-partial:]
    return records


def _rollback_final_cluster_after_error(
    root: Path, completed: List[Dict[str, object]],
) -> List[Dict[str, object]]:
    if len(completed) < CLUSTER_SIZE or len(completed) % CLUSTER_SIZE:
        raise TimingError("cannot reserve a complete final cluster for replay")
    records = []
    for item in completed[-CLUSTER_SIZE:]:
        job = int(item["job"])
        entry = root / "ledger" / ("%04d" % job)
        records.append(_archive_directory(
            root, entry, "task-%04d.unbound" % job,
            reason="no terminal thermal binding"))
    del completed[-CLUSTER_SIZE:]
    return records


def _rollback_unlisted_task_entries(
    root: Path, preexisting_jobs: Iterable[int],
    completed: Sequence[Mapping[str, object]],
) -> List[Dict[str, object]]:
    expected = {int(job) for job in preexisting_jobs} | {
        int(item["job"]) for item in completed}
    records = []
    actual = set()
    for entry in sorted((root / "ledger").iterdir()):
        if entry.is_symlink() or not entry.is_dir() or \
                TASK_DIR_RE.fullmatch(entry.name) is None:
            raise TimingError("task ledger contains an unsafe unlisted entry")
        job = int(entry.name)
        actual.add(job)
        if job not in expected:
            records.append(_archive_directory(
                root, entry, "task-%04d.unbound" % job,
                reason="no terminal thermal binding"))
    if not expected <= actual:
        raise TimingError("a terminal-listed task entry disappeared")
    return records


def _write_campaign_receipt(root: Path, tasks: Sequence[Mapping[str, object]]) -> None:
    task_receipts = [{
        "job": task["job"],
        "receipt_sha256": sha256_file(
            root / "ledger" / ("%04d" % int(task["job"])) / "receipt.json"),
    } for task in tasks]
    terminals = sorted((root / "segments").glob("segment*.terminal.json"))
    terminal_receipts = [{"name": path.name, "sha256": sha256_file(path)}
                         for path in terminals]
    exact = {
        "design_sha256": sha256_file(root / "design.json"),
        "prepare_receipt_sha256": sha256_file(root / "prepare_receipt.json"),
        "prelaunch_receipt_sha256": sha256_file(root / "prelaunch_receipt.json"),
        "tasks_manifest_sha256": sha256_file(root / "tasks_manifest.jsonl"),
        "task_count": len(tasks), "task_receipts": task_receipts,
        "terminal_receipts": terminal_receipts,
    }
    path = root / "campaign_receipt.json"
    if path.exists() or path.is_symlink():
        campaign = verify_sealed(
            load_canonical(path, "campaign receipt"), SCHEMA + ".campaign",
            "campaign receipt")
        if set(campaign) != {
                "schema", "self_sha256_excluding_field", "completed_utc", *exact} or \
                any(campaign.get(key) != value for key, value in exact.items()):
            raise TimingError("existing campaign receipt does not replay")
        validate_utc_timestamp(campaign.get("completed_utc"),
                               "campaign completion")
        return
    atomic_write(path, canonical_json(sealed_record(SCHEMA + ".campaign", {
        "completed_utc": utc_now(), **exact,
    })))


def _run_campaign_locked(args: argparse.Namespace, root: Path,
                         design: Mapping[str, object], termination) -> None:
    isolation = isolation_module()
    recover_derived_atomic(
        root / "campaign_receipt.json", SCHEMA + ".campaign",
        "campaign receipt")
    recover_derived_atomic(
        root / "validated_summary.json", SCHEMA + ".summary",
        "validated summary")
    verify_root_layout(root)
    tasks = load_tasks(root, design)
    verify_immutable_inputs(root, design)
    validate_prepare(root, design)
    verify_active_runner(root, design)
    if args.max_tasks <= 0 or args.max_tasks % int(design["cluster_size"]) != 0 or \
            args.timeout_seconds <= 0:
        raise TimingError("run bounds must be positive and preserve four-task clusters")
    os.sched_setaffinity(0, {int(design["controller_core"])})
    if os.sched_getaffinity(0) != {int(design["controller_core"])}:
        raise TimingError("controller affinity could not be sealed")
    for core_key, receipt_key in (
            ("core", "timing_topology"),
            ("controller_core", "controller_topology"),
            ("thermal_core", "thermal_topology")):
        current = isolation.topology_record(
            int(design[core_key]), int(design["numa_node"]))
        if current != design[receipt_key]:
            raise TimingError("runtime CPU topology/power policy changed: " + core_key)
    isolation.recover_orphaned_thermal_segments(root, design)
    bindings = terminal_task_bindings(root)
    recover_unbound_task_entries(root, bindings)
    recover_interrupted_transactions(root)
    recovered = adopt_unbound_interrupted_archives(root)
    entries = existing_task_entries(root, design, tasks)
    verify_segments(root, design, entries, allow_unbound_archives=True)
    if set(entries) != set(bindings):
        raise TimingError("resumed task ledger and thermal bindings disagree")
    expected_prefix = set(range(len(entries)))
    if set(entries) != expected_prefix or len(entries) % int(design["cluster_size"]):
        raise TimingError("resumed task ledger is not a whole-cluster prefix")
    remaining = tasks[len(entries):]
    campaign_path = root / "campaign_receipt.json"
    if remaining and (campaign_path.exists() or campaign_path.is_symlink()):
        raise TimingError("incomplete campaign has a terminal campaign receipt")
    summary_path = root / "validated_summary.json"
    if remaining and (summary_path.exists() or summary_path.is_symlink()):
        raise TimingError("incomplete campaign has a validated summary")
    if not remaining:
        if recovered:
            raise TimingError("unbound recovery evidence has no remaining task")
        if (summary_path.exists() or summary_path.is_symlink()) and \
                not campaign_path.exists():
            raise TimingError(
                "validated summary exists without its campaign receipt")
        termination.raise_if_requested()
        _write_campaign_receipt(root, tasks)
        print(json.dumps({"complete": True, "resumed": len(entries),
                          "completed": 0}, sort_keys=True))
        return
    selected = remaining[:args.max_tasks]
    if len(selected) % int(design["cluster_size"]):
        # The final bounded segment is allowed to contain fewer than max_tasks,
        # but the sealed Cartesian product itself is a whole number of clusters.
        raise TimingError("remaining manifest unexpectedly splits a paired cluster")
    launch_preflight = preflight_report(root, design, measure_quiet=True)
    if not launch_preflight["ready"]:
        raise TimingError("quiet timing launch prerequisites are not met: %s" %
                          launch_preflight["blockers"])
    segment = _next_segment(root)
    termination.raise_if_requested()
    owner = isolation.start_sampler(root, design, segment)
    timing_start_s = time.monotonic()
    started_utc = utc_now()
    completed: List[Dict[str, object]] = []
    error: Optional[BaseException] = None
    last_live_check_s = -math.inf
    initial_thermal_gate: Dict[str, object] = {}

    def monitor_owned_sampler(*, force: bool) -> None:
        nonlocal last_live_check_s, initial_thermal_gate
        termination.raise_if_requested()
        poll_time = time.monotonic()
        if not force and poll_time - last_live_check_s < 0.75:
            return
        receipt = isolation.enforce_live_thermal_safety(owner.csv_part)
        if not initial_thermal_gate:
            initial_thermal_gate = dict(receipt)
        last_live_check_s = time.monotonic()

    try:
        monitor_owned_sampler(force=True)
        for task in selected:
            completed_task = _task_execution(
                root, design, task, args.timeout_seconds,
                poll_callback=lambda: monitor_owned_sampler(force=False))
            monitor_owned_sampler(force=True)
            tools = design["tools"]
            if conflicting_processes(excluded_pids=(owner.pid, owner.process.pid)):
                raise TimingError("a conflicting benchmark/filler/sampler appeared")
            if not isolation.process_identity_matches(
                    owner.identity, int(design["thermal_core"]), owner.csv_part) or \
                    isolation.sole_i2c_readers(
                        Path(tools["fuser"]["path"]),
                        Path(tools["sudo"]["path"]),
                        Path(tools["timeout"]["path"])) != (owner.pid,):
                raise TimingError("thermal ownership changed during timing")
            completed.append(completed_task)
            if len(completed) % 12 == 0 or len(completed) == len(selected):
                print("segment=%d progress=%d/%d total=%d/%d" %
                      (segment, len(completed), len(selected),
                       len(entries) + len(completed), len(tasks)), flush=True)
    except BaseException as exc:
        error = exc
    termination_was_primary = isinstance(
        error, getattr(isolation, "ControllerTermination", ()))
    benchmark_end_s = time.monotonic()
    terminal_payload: Dict[str, object] = {}
    stop_error: Optional[BaseException] = None
    try:
        terminal_payload = isolation.stop_sampler(
            owner, root, design, timing_start_s, benchmark_end_s)
    except BaseException as exc:
        stop_error = exc
    if stop_error is not None and not terminal_payload:
        try:
            isolation.archive_invalid_thermal_segment(
                root, design, segment,
                "sampler-finalization-failure: %s" % stop_error)
        except BaseException as archive_error:
            stop_error = TimingError(
                "%s; invalid thermal evidence could not be archived: %s" %
                (stop_error, archive_error))
    if stop_error is not None:
        error = stop_error if error is None else TimingError(
            "task failure %s; sampler finalization failure %s" %
            (error, stop_error))
    if not terminal_payload and stop_error is None:
        empty_stop = TimingError("sampler stop returned no terminal evidence")
        error = empty_stop if error is None else TimingError(
            "campaign failure %s; sampler evidence failure %s" %
            (error, empty_stop))
    if terminal_payload and not initial_thermal_gate:
        missing_gate = TimingError(
            "sampler interval lacks a successful initial live thermal gate")
        try:
            isolation.archive_invalid_thermal_segment(
                root, design, segment, str(missing_gate))
        except BaseException as archive_error:
            missing_gate = TimingError(
                "%s; invalid thermal evidence could not be archived: %s" %
                (missing_gate, archive_error))
        terminal_payload = {}
        error = missing_gate if error is None else TimingError(
            "campaign failure %s; thermal gate failure %s" %
            (error, missing_gate))
    termination_error = termination.error()
    if termination_error is not None and not termination_was_primary:
        error = termination_error if error is None else TimingError(
            "campaign failure %s; deferred termination %s" %
            (error, termination_error))
    if terminal_payload:
        recovered.extend(_rollback_partial_cluster(root, completed))
        if error is not None and \
                len(entries) + len(completed) == len(tasks):
            recovered.extend(_rollback_final_cluster_after_error(
                root, completed))
        recovered.extend(_rollback_unlisted_task_entries(
            root, entries, completed))
        verify_interrupted_records(
            root, recovered, expected_thermal_segment=segment)
        status = "failed" if error is not None else \
            "complete" if len(entries) + len(completed) == len(tasks) else \
            "bounded-incomplete"
        terminal = sealed_record(SCHEMA + ".segment_terminal", {
            "segment": segment, "status": status,
            "started_utc": started_utc, "ended_utc": utc_now(),
            "timing_start_monotonic_s": timing_start_s,
            "benchmark_end_monotonic_s": benchmark_end_s,
            "duration_s": benchmark_end_s - timing_start_s,
            "controller_pid": os.getpid(),
            "controller_affinity": sorted(os.sched_getaffinity(0)),
            "design_sha256": sha256_file(root / "design.json"),
            "prepare_receipt_sha256": sha256_file(root / "prepare_receipt.json"),
            "prelaunch_receipt_sha256": sha256_file(
                root / "prelaunch_receipt.json"),
            "tasks_manifest_sha256": sha256_file(root / "tasks_manifest.jsonl"),
            "launch_preflight": launch_preflight,
            "initial_thermal_gate": initial_thermal_gate,
            "resumed_jobs_before": sorted(entries),
            "completed_tasks": completed,
            "recovered_interrupted_transactions": recovered,
            "failure": None if error is None else {
                "class": type(error).__name__, "message": str(error)},
            **terminal_payload,
        })
        terminal_path = root / "segments" / (
            "segment%04d.terminal.json" % segment)
        atomic_write(terminal_path, canonical_json(terminal))
        sealed_entries = existing_task_entries(root, design, tasks)
        verify_segments(root, design, sealed_entries)
    if error is not None:
        raise error
    if len(entries) + len(completed) == len(tasks):
        _write_campaign_receipt(root, tasks)
    print(json.dumps({
        "segment": segment, "status": status, "completed": len(completed),
        "resumed": len(entries),
        "remaining": len(tasks) - len(entries) - len(completed),
        "terminal_sha256": sha256_file(root / "segments" /
            ("segment%04d.terminal.json" % segment)),
    }, sort_keys=True))


def run_campaign(args: argparse.Namespace) -> None:
    isolation = isolation_module()
    root = Path(args.result_dir).resolve()
    design = load_design(root)
    lock_path = root / "controller.lock"
    if lock_path.is_symlink() or not lock_path.is_file() or \
            sha256_file(lock_path) != design.get("controller_lock_sha256"):
        raise TimingError("controller lock identity changed")
    with lock_path.open("rb") as stream:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise TimingError("another campaign controller owns this artifact") from exc
        with isolation.DeferredTermination() as termination:
            _run_campaign_locked(args, root, design, termination)


def verify_campaign_data(root: Path, *, require_complete: bool) -> Tuple[
        Dict[str, object], List[Dict[str, object]], Dict[int, Dict[str, object]],
        List[Dict[str, object]]]:
    design = load_design(root)
    verify_root_layout(root)
    verify_active_runner(root, design)
    tasks = load_tasks(root, design)
    verify_immutable_inputs(root, design)
    validate_prepare(root, design)
    if any((root / ".transactions").iterdir()):
        raise TimingError("unrecovered task transaction remains")
    entries = existing_task_entries(root, design, tasks)
    segments = verify_segments(root, design, entries)
    complete = len(entries) == len(tasks)
    campaign_path = root / "campaign_receipt.json"
    if require_complete and not complete:
        raise TimingError("campaign task ledger is incomplete")
    if complete:
        if campaign_path.is_symlink() or not campaign_path.is_file():
            raise TimingError("complete task ledger lacks a campaign receipt")
        campaign = verify_sealed(
            load_canonical(campaign_path, "campaign receipt"),
            SCHEMA + ".campaign", "campaign receipt")
        if set(campaign) != {
                "schema", "self_sha256_excluding_field", "completed_utc",
                "design_sha256", "prepare_receipt_sha256",
                "prelaunch_receipt_sha256", "tasks_manifest_sha256",
                "task_count", "task_receipts", "terminal_receipts"}:
            raise TimingError("campaign receipt field set changed")
        expected = {
            "design_sha256": sha256_file(root / "design.json"),
            "prepare_receipt_sha256": sha256_file(root / "prepare_receipt.json"),
            "prelaunch_receipt_sha256": sha256_file(
                root / "prelaunch_receipt.json"),
            "tasks_manifest_sha256": sha256_file(root / "tasks_manifest.jsonl"),
            "task_count": len(tasks),
            "task_receipts": [{
                "job": task["job"],
                "receipt_sha256": entries[int(task["job"])]["receipt_sha256"],
            } for task in tasks],
            "terminal_receipts": [{
                "name": "segment%04d.terminal.json" % index,
                "sha256": sha256_file(root / "segments" /
                    ("segment%04d.terminal.json" % index)),
            } for index in range(len(segments))],
        }
        if any(campaign.get(key) != value for key, value in expected.items()):
            raise TimingError("campaign receipt does not replay")
        validate_utc_timestamp(campaign.get("completed_utc"),
                               "campaign completion")
    elif campaign_path.exists() or campaign_path.is_symlink():
        raise TimingError("incomplete task ledger has a campaign receipt")
    return design, tasks, entries, segments


def _summary_payload(root: Path, tasks: Sequence[Mapping[str, object]],
                     entries: Mapping[int, Mapping[str, object]],
                     segments: Sequence[Mapping[str, object]],
                     ) -> Dict[str, object]:
    aggregates = aggregate_rows(
        tasks, {job: value["parsed"] for job, value in entries.items()})
    gates = {
        candidate: aggregates["candidates"][candidate]["overall"]["speed_gate"]
        for candidate in CANDIDATES
    }
    pareto = aggregates["candidate_pareto_comparison"]
    recovery = _validate_recovery_result(json.loads(
        (root / "external" / RECOVERY_RESULT_NAME).read_text(encoding="utf-8")))
    speed_eligible = [candidate for candidate in ("four-0369", "three-048")
                      if gates[candidate]["speed_gate_passed"] is True]
    selection = select_pareto_architecture(gates, pareto)
    if selection is None:
        selection_reason = "neither raw segmented layout passed its speed gate"
    elif selection == "three-048" and gates["four-0369"][
            "speed_gate_passed"] is True:
        selection_reason = (
            "three-048 has a material and one-sided-95-supported paired "
            "solve-speed advantage over four-0369")
    elif selection == "three-048":
        selection_reason = (
            "three-048 passed its baseline speed gate while four-0369 did not")
    elif gates["three-048"]["speed_gate_passed"] is True:
        selection_reason = (
            "both layouts passed; three-048 lacked a material, "
            "one-sided-95-supported speed advantage, so four-0369 retains "
            "the stronger raw recovery")
    else:
        selection_reason = (
            "four-0369 passed its baseline speed gate while three-048 did not")

    arm_metrics: Dict[str, object] = {}
    for arm in ("baseline", "three-048", "four-0369"):
        source = recovery["arms"][arm]
        arm_metrics[arm] = {
            "failures": source["failures"],
            "failure_reduction_vs_baseline":
                source.get("failure_reduction_vs_baseline", 0),
            "failure_reduction_percent":
                source.get("failure_reduction_percent", "0"),
            "failures_by_schedule": source["failures_by_schedule"],
            "failures_by_seed_index": source["failures_by_seed_index"],
            "field_shortfall_failures": source["field_shortfall_failures"],
            "q_gt_12_failures": source["q_gt_12_failures"],
            "weak_K": source["weak_K"],
            "maximum_weak_K_multiplicity":
                source["maximum_weak_K_multiplicity"],
            "weak_K_multiplicity_histogram":
                source["weak_K_multiplicity_histogram"],
            "seed_attempt_histogram": source["seed_attempt_histogram"],
        }
    comparison_metrics: Dict[str, object] = {}
    for candidate in CANDIDATES:
        source = recovery["comparisons_to_baseline"][candidate]
        comparison_metrics[candidate] = {
            "repairs": source["repairs"],
            "introductions": source["introductions"],
            "repair_q_gt_12": source["repair_q_gt_12"],
            "introduction_q_gt_12": source["introduction_q_gt_12"],
            "both_success": source["both_success"],
            "both_fail": source["both_fail"],
            "net_failure_change": source["net_failure_change"],
            "common_success_work_ratios":
                source["common_success_work_ratios"],
            "seed_attempt_delta_cells":
                source["seed_attempt_delta_cells"],
            "seed_attempt_candidate_minus_baseline_sum":
                source["seed_attempt_candidate_minus_baseline_sum"],
            "q_delta_cells": source["q_delta_cells"],
            "q_candidate_minus_baseline_sum":
                source["q_candidate_minus_baseline_sum"],
        }
    return {
        "design_sha256": sha256_file(root / "design.json"),
        "campaign_receipt_sha256": sha256_file(root / "campaign_receipt.json"),
        "task_count": len(tasks), "segment_count": len(segments),
        "comparison": "raw-three-048-and-four-0369-vs-raw-two-anchor",
        "timing_evidence_promotional": True,
        "architecture_seed_repairs": "none",
        "imported_p32_helper_source_provenance":
            P32_HELPER_SOURCE_PROVENANCE,
        "architecture_dimensions": {
            "dense_rows": 12, "gf256_rows": 10, "gf16_rows": 2,
            "period": PERIOD, "grouped_rows": GROUPED_ROWS,
            "buckets": BUCKETS, "mix_count": 2,
        },
        "recovery_evidence": {
            "receipt_sha256": sha256_file(
                root / "external" / RECOVERY_RESULT_NAME),
            "validated_summary_sha256":
                recovery["evidence"]["validated_summary_sha256"],
            "data_manifest_sha256":
                recovery["evidence"]["data_manifest_sha256"],
            "baseline_failures": recovery["arms"]["baseline"]["failures"],
            "three_048_failures": recovery["arms"]["three-048"]["failures"],
            "four_0369_failures": recovery["arms"]["four-0369"]["failures"],
            "three_048_failure_reduction_vs_baseline":
                recovery["arms"]["three-048"][
                    "failure_reduction_vs_baseline"],
            "four_0369_failure_reduction_vs_baseline":
                recovery["arms"]["four-0369"][
                    "failure_reduction_vs_baseline"],
            "arm_metrics": arm_metrics,
            "comparisons_to_baseline": comparison_metrics,
        },
        "candidate_speed_gates": gates,
        "candidate_pareto_comparison": pareto,
        "speed_eligible_recovery_improvements": speed_eligible,
        "architecture_promotion_ready": selection is not None,
        "architecture_promotion_selection": selection,
        "architecture_promotion_selection_reason": selection_reason,
        "architecture_promotion_selection_policy":
            "among baseline-speed-eligible layouts, select three-048 only "
            "when its cotimed-baseline-normalized solve ratio is at least "
            "0.5% lower than four-0369 and its one-sided 95% bootstrap upper bound is below "
            "parity; otherwise retain four-0369 for stronger raw recovery",
        "architecture_promotion_blocker": None if selection is not None else
            "neither raw segmented layout passed its full-payload paired speed gate",
        "aggregates": aggregates,
    }


def reduce_campaign(args: argparse.Namespace) -> None:
    root = Path(args.result_dir).resolve()
    recover_derived_atomic(
        root / "campaign_receipt.json", SCHEMA + ".campaign",
        "campaign receipt")
    recover_derived_atomic(
        root / "validated_summary.json", SCHEMA + ".summary",
        "validated summary")
    _design, tasks, entries, segments = verify_campaign_data(
        root, require_complete=True)
    path = root / "validated_summary.json"
    payload = _summary_payload(root, tasks, entries, segments)
    if path.exists() or path.is_symlink():
        summary = verify_sealed(
            load_canonical(path, "validated summary"), SCHEMA + ".summary",
            "validated summary")
        if set(summary) != {
                "schema", "self_sha256_excluding_field", "created_utc", *payload} or \
                any(summary.get(key) != value for key, value in payload.items()):
            raise TimingError("existing validated summary does not replay")
        validate_utc_timestamp(summary.get("created_utc"),
                               "validated summary creation")
    else:
        summary = sealed_record(SCHEMA + ".summary", {
            "created_utc": utc_now(), **payload,
        })
        atomic_write(path, canonical_json(summary))
    print(json.dumps({
        "validated_summary_sha256": sha256_file(path),
        "architecture_promotion_ready": summary["architecture_promotion_ready"],
        "candidate_speed_gates": summary["candidate_speed_gates"],
    }, sort_keys=True))


def verify_campaign(args: argparse.Namespace) -> None:
    root = Path(args.result_dir).resolve()
    _design, tasks, entries, segments = verify_campaign_data(
        root, require_complete=not args.allow_incomplete)
    summary_path = root / "validated_summary.json"
    if summary_path.exists() or summary_path.is_symlink():
        if len(entries) != len(tasks):
            raise TimingError("incomplete campaign has a validated summary")
        summary = verify_sealed(
            load_canonical(summary_path, "validated summary"),
            SCHEMA + ".summary", "validated summary")
        payload = _summary_payload(root, tasks, entries, segments)
        if set(summary) != {
                "schema", "self_sha256_excluding_field", "created_utc", *payload} or \
                any(summary.get(key) != value for key, value in payload.items()):
            raise TimingError("validated summary does not replay")
        validate_utc_timestamp(summary.get("created_utc"),
                               "validated summary creation")
    print(json.dumps({
        "verified": True, "complete": len(entries) == len(tasks),
        "task_count": len(tasks), "completed_tasks": len(entries),
        "segment_count": len(segments),
        "design_sha256": sha256_file(root / "design.json"),
        "campaign_receipt_sha256": sha256_file(root / "campaign_receipt.json")
            if (root / "campaign_receipt.json").exists() else None,
        "validated_summary_sha256": sha256_file(summary_path)
            if summary_path.exists() else None,
    }, sort_keys=True))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare = subparsers.add_parser("prepare")
    prepare.add_argument("--repo", required=True)
    prepare.add_argument("--binary", required=True)
    prepare.add_argument("--result-dir", required=True)
    prepare.add_argument("--core", type=int, default=8)
    prepare.add_argument("--controller-core", type=int, default=126)
    prepare.add_argument("--thermal-core", type=int, default=127)
    prepare.add_argument("--numa-node", type=int, default=0)
    prepare.add_argument("--evict-bytes", type=int, default=268435456)
    prepare.add_argument("--build-jobs", type=int, default=12)
    prepare.add_argument("--build-cpu-set", default="0-7,64-71")
    preflight = subparsers.add_parser("preflight")
    preflight.add_argument("--result-dir", required=True)
    run = subparsers.add_parser("run")
    run.add_argument("--result-dir", required=True)
    run.add_argument("--max-tasks", type=int, default=36)
    run.add_argument("--timeout-seconds", type=float, default=1800.0)
    verify = subparsers.add_parser("verify")
    verify.add_argument("--result-dir", required=True)
    verify.add_argument("--allow-incomplete", action="store_true")
    reduce_parser = subparsers.add_parser("reduce")
    reduce_parser.add_argument("--result-dir", required=True)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = build_parser().parse_args(argv)
    if args.command == "prepare":
        prepare_campaign(args)
    elif args.command == "preflight":
        preflight_command(args)
    elif args.command == "run":
        run_campaign(args)
    elif args.command == "verify":
        verify_campaign(args)
    elif args.command == "reduce":
        reduce_campaign(args)
    else:
        raise TimingError("unknown command")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("segmented-anchor timing failed: %s" % exc,
              file=sys.stderr, flush=True)
        raise
