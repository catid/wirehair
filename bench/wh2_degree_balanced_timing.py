#!/usr/bin/env python3
"""Sealed full-payload timing for the WH2 degree-balanced staircase.

The architecture commit is an exact direct child of the promoted lazy-arena
base and replays only the reviewed degree-balanced patch.  Preparation freezes
the exact base, exact architecture, and a later measurement-only binary.  A
bounded semantic gate proves the architecture's flag-off path replays the
exact base, while the measurement binary receipts the expected staircase-only
graph difference.  Promotional timing uses just that one measurement binary,
swapping the semantic baseline/balanced arms in an outer ABBABAAB order to
avoid cross-binary code-layout bias.  Every process runs four inner ABBABAAB
cycles and discards inner cycle zero, yielding 96 timed solves per architecture
in each of 252 immutable logical tasks.

``prepare`` builds, audits, smokes, and freezes the campaign but never launches
quiet timing.  ``run`` is fresh-only, automatically reuses exactly one compatible
external CPU/eight-DIMM/EDAC reader without ever signalling it, and transactionally
copies only the minimal campaign-bracketing thermal slice.  ``reduce`` replays
every raw byte and reports unconditional attempt latency and outcome strata;
common-success speed is promotional only when the entire fixed panel succeeds.
No command changes a production profile.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import io
import json
import math
import os
from pathlib import Path
import random
import re
import shutil
import stat
import subprocess
import sys
import time
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


sys.dont_write_bytecode = True

BASE_COMMIT = "07c13ecfe03c45e9195336792606c3d74a7be262"
ARCHITECTURE_COMMIT = "176c098dadeab4666bd04a79c199e3d1ec0d0bbe"
MEASUREMENT_COMMIT = "13d4e548e13a54dfe278970b25495769188d91ec"
ARCHITECTURE_PATCH_ID = "b1c882028f283d0d5ca6bb1298f6b6ae4863f716"
GROUPED_MASKOPT_COMMIT = "3d6d9f65abef106b5d44c7bbc500fb103d92982c"
RECOVERY_CONTEXT_COMMIT = "c4ba11a666256ee8a8802d1e595cd111126d6998"
RECOVERY_CONTEXT_PATH = "bench/wh2_degree_balanced_allk_result.json"
RECOVERY_CONTEXT_SHA256 = \
    "4507d0f6717cfbd866e7ca24e1639818c5f6946e22b864f37dadc2f4a5f24ffc"
RECOVERY_CONTROLLER_PATH = "bench/wh2_degree_balanced_allk_campaign.py"
RECOVERY_CONTROLLER_SHA256 = \
    "a0cd8ff51f3a62cfc5dd655c84dc5e2e9dc443fa882297dbf2743dcab0f04b06"
RECOVERY_BASE_COMMIT = "3c8a7a31e0e0b592f906d52bf7fef0bafebbf2f9"
RECOVERY_CANDIDATE_COMMIT = "bf418cabbfc70e8da89e938a6f68c439e4234948"
FROZEN_MIXED_COEFFICIENT_FINGERPRINT = "0xdcff9be773c1a37c"
KS = (3200, 9999, 10000, 20000, 32000, 48466, 64000)
TIMING_STAIRCASE_ROWS = {
    "3200": 62, "9999": 86, "10000": 86, "20000": 134,
    "32000": 190, "48466": 374, "64000": 346,
}
WIDTHS = (64, 1280, 4096)
SCHEDULE_SEEDS = (
    ("burst", 0, 0xd43616aa4e9a624b),
    ("burst", 1, 0x08527a72b7e08330),
    ("adversarial", 0, 0xa459e3d5dd440150),
    ("adversarial", 1, 0x3a7618e056091af3),
    ("repair-only", 0, 0x65a2c2234cf67452),
    ("repair-only", 1, 0x6b4071fb2a20757d),
)
CACHE_STATES = ("cold", "warm")
OUTER_ORDER = "ABBABAAB"
INNER_ORDER = "ABBABAAB"
ARCHITECTURE = {
    "period": 244, "grouped_rows": 0, "buckets": "auto",
    "geometry": "frozen",
    "gf256_rows": 10, "gf16_rows": 2, "mix_count": 2,
    "mixed_residue_skew": 0, "mixed_residue_schedule": "constant",
    "mixed_residue_hash_seed": 0,
    "mixed_independent_extension_residues": False,
    "mixed_extension_residue_seed_xor": 78,
    "packet_peel_seed_xor": 0, "packet_peel_seed_table": "none",
    "binary_dense_rows_override": 0, "gf256_heavy_rows_override": 0,
    "source_hits_override": 0, "seed_block_bytes_override": 0,
    "packet_row_seed_multiplier": 1, "packet_row_seed_avalanche": False,
    "odd_packet_peel_seed_xor": 0,
    "source_hits_policy": "certified-by-K",
    "source_hits_transition_K": 10000,
    "source_hits_below_transition": 2,
    "source_hits_at_or_above_transition": 3,
    "field": 1,
    "dense_identity_corner": False, "heavy_family": 0,
    "staircase_rows_policy": "production-GetDenseCount",
    "staircase_rows_by_timing_K": TIMING_STAIRCASE_ROWS,
    "dense_rows": 12, "heavy_rows": 12, "dense_two_anchor": True,
    "dense_two_anchor_phase": 0,
    "mixed_coefficient_fingerprint":
        FROZEN_MIXED_COEFFICIENT_FINGERPRINT,
}
BINARY_NAMES = {
    "exact_base": "wirehair_v2_bench.exact_base",
    "exact_architecture": "wirehair_v2_bench.exact_architecture",
    "measurement": "wirehair_v2_bench.measurement",
}
NEUTRAL_SMOKE_SPECS = (
    ("exact_base", "exact_base", False, "shared-x", "primary"),
    ("exact_architecture", "exact_architecture", False, "shared-x", "primary"),
    ("measurement_shared_x", "measurement", True, "shared-x", "primary"),
    ("measurement_frozen", "measurement", True, "frozen", "primary"),
    ("measurement_frozen_K10000", "measurement", True, "frozen", "boundary"),
)
OUTER_ORIENTATIONS = {"A": "forward", "B": "reverse"}
OVERHEAD = 4
LOSS_TEXT = "0.5"
MALLOC_MMAP_THRESHOLD = "1073741824"
MALLOC_TRIM_THRESHOLD = "-1"
MAX_ENVIRONMENTAL_ATTEMPTS = 10
MAX_MINOR_FAULTS = 64
TIMING_CORE = 8
TIMING_SIBLING = 72
CONTROLLER_CORE = 126
SAMPLER_CORE = 127
NUMA_NODE = 0
MAX_CPU_TEMP_C = 85.0
MAX_DIMM_TEMP_C = 90.0
MAX_THERMAL_GAP_S = 2.25
MAX_THERMAL_MARGIN_S = 2.25
BOOTSTRAP_REPS = 20000

SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
UINT_RE = re.compile(r"0|[1-9][0-9]*\Z")
SINT_RE = re.compile(r"0|-?[1-9][0-9]*\Z")
HEX_RE = re.compile(r"0x(?:0|[1-9a-f][0-9a-f]*)\Z")
UTC_RE = re.compile(
    r"[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:"
    r"[0-9]{2}\.[0-9]{3}Z\Z")

TIMING_FIELDS = (
    "N", "bb", "overhead", "schedule", "seed", "loss", "cache_state",
    "cycle", "slot", "arm", "period", "grouped_rows",
    "buckets_requested", "mixed_geometry", "degree_balanced_staircase",
    "seed_attempt",
    "matrix_seed", "peel_seed",
    "preflight_result", "cell_class", "common_success", "result",
    "outcome_stable", "elapsed_ns", "saturated", "cpu_before",
    "cpu_after", "cpu_migrated", "minflt_delta", "majflt_delta",
    "fault_contaminated", "inactivated", "binary_def", "heavy_gain",
    "block_xors", "block_muladds", "build_ns", "peel_ns", "project_ns",
    "residual_ns", "backsub_ns", "joint_source_xors",
    "joint_marginal_xors", "joint_marginal_copies", "joint_active_deltas",
    "joint_scratch_bytes", "dual_source_columns", "source_bytes",
    "packet_payload_bytes", "intermediate_bytes",
    "solve_value_arena_bytes", "solve_value_eager_zero_bytes",
    "solve_value_commit_copy_bytes", "rhs_route_expected",
    "rhs_route_actual",
)

PREAMBLE_V5 = (
    "schema", "policy", "timing_scope", "cycles", "order",
    "discard_cycle", "cycle_mode", "cycle_index", "N", "bb", "overhead",
    "loss", "seed", "schedule", "cache_state", "overhead_stream",
    "evict_bytes", "eviction_prefaulted", "control_period",
    "control_grouped_rows", "control_buckets", "control_geometry",
    "control_secondary_schedule", "control_rhs_route_expected",
    "control_preflight_rhs_route",
    "control_degree_balanced_staircase", "control_grouped_hash_seed",
    "control_final_h_a_columns", "control_staircase_fingerprint",
    "control_dense_fingerprint", "control_heavy_fingerprint",
    "control_mixed_coefficient_fingerprint",
    "candidate_period", "candidate_grouped_rows", "candidate_buckets",
    "candidate_geometry", "candidate_secondary_schedule",
    "candidate_rhs_route_expected", "candidate_preflight_rhs_route",
    "candidate_degree_balanced_staircase",
    "candidate_grouped_hash_seed", "candidate_final_h_a_columns",
    "candidate_staircase_fingerprint", "candidate_dense_fingerprint",
    "candidate_heavy_fingerprint", "candidate_mixed_coefficient_fingerprint",
    "staircase_rows", "dense_rows", "heavy_rows", "columns",
    "gf256_rows", "gf16_rows", "mixed_residue_skew",
    "mixed_residue_schedule", "mixed_residue_hash_seed",
    "mixed_independent_extension_residues",
    "mixed_extension_residue_seed_xor", "packet_peel_seed_xor",
    "packet_peel_seed_table", "binary_dense_rows_override",
    "gf256_heavy_rows_override", "source_hits_override",
    "seed_block_bytes_override", "packet_row_seed_multiplier",
    "packet_row_seed_avalanche", "odd_packet_peel_seed_xor",
    "source_hits", "field", "dense_identity_corner", "heavy_family",
    "dense_two_anchor", "dense_two_anchor_phase", "control_attempt",
    "control_matrix_seed", "control_peel_seed", "candidate_attempt",
    "candidate_matrix_seed", "candidate_peel_seed", "mix", "payload",
    "payload_count", "payload_bytes", "payload_alignment",
    "payload_prefaulted", "payload_fingerprint", "packet_trace_sha256",
    "system_build", "tls_reapply", "allocator_tls_state",
    "solve_value_storage", "solve_value_publish", "preflight_control_result",
    "preflight_candidate_result", "preflight_control_decoded_fingerprint",
    "preflight_candidate_decoded_fingerprint", "cell_class",
    "common_success", "trace_sha256",
)

LEGACY_PREAMBLE_V2 = (
    "schema", "policy", "timing_scope", "cycles", "order",
    "discard_cycle", "cycle_mode", "cycle_index", "N", "bb", "overhead",
    "loss", "seed", "schedule", "cache_state", "overhead_stream",
    "evict_bytes", "eviction_prefaulted", "control_period",
    "control_grouped_rows", "control_buckets", "control_grouped_hash_seed",
    "control_final_h_a_columns", "candidate_period",
    "candidate_grouped_rows", "candidate_buckets",
    "candidate_grouped_hash_seed", "candidate_final_h_a_columns",
    "gf256_rows", "gf16_rows", "dense_two_anchor", "control_attempt",
    "control_matrix_seed", "control_peel_seed", "candidate_attempt",
    "candidate_matrix_seed", "candidate_peel_seed", "mix", "payload",
    "payload_count", "payload_bytes", "payload_alignment",
    "payload_prefaulted", "system_build", "tls_reapply",
    "allocator_tls_state", "solve_value_storage", "solve_value_publish",
    "preflight_control_result", "preflight_candidate_result", "cell_class",
    "common_success", "trace_sha256",
)
LEGACY_FIELDS = tuple(
    field for field in TIMING_FIELDS
    if field not in {
        "mixed_geometry", "degree_balanced_staircase",
        "rhs_route_expected", "rhs_route_actual"})
FINGERPRINT_RE = re.compile(r"0x[0-9a-f]{16}\Z")

WORK_FIELDS = (
    "result", "inactivated", "binary_def", "heavy_gain", "block_xors",
    "block_muladds", "joint_source_xors", "joint_marginal_xors",
    "joint_marginal_copies", "joint_active_deltas", "joint_scratch_bytes",
    "dual_source_columns", "source_bytes", "packet_payload_bytes",
    "intermediate_bytes",
)
THERMAL_FIELDS = (
    "utc", "monotonic_s", "cpu_busy_pct", "cpu_avg_mhz", "cpu_tctl_c",
    "dimm_i2c1_50_c", "dimm_i2c1_51_c", "dimm_i2c1_52_c",
    "dimm_i2c1_53_c", "dimm_i2c2_50_c", "dimm_i2c2_51_c",
    "dimm_i2c2_52_c", "dimm_i2c2_53_c", "dimm_read_errors", "load1",
    "load5", "load15", "edac_ce", "edac_ue",
)
DIMM_FIELDS = tuple(field for field in THERMAL_FIELDS if field.startswith("dimm_i2c"))

EXECUTION_RECEIPT_FIELDS = frozenset((
    "schema", "self_sha256_excluding_field", "job", "task_id", "task_sha256",
    "outer_slot", "outer_marker", "orientation", "binary_sha256", "argv",
    "attempt", "started_utc", "start_monotonic_ns", "end_monotonic_ns",
    "duration_ns", "stderr_sha256", "prior_contamination_receipts",
    "schema_version", "stdout_sha256", "semantic_sha256", "trace_sha256",
    "work_signatures_sha256", "semantic_outcome_class",
    "base_timed_elapsed_ns", "balanced_timed_elapsed_ns", "all_elapsed_ns",
    "timed_minor_faults", "discard_minor_faults", "max_minor_faults",
    "row_count", "timed_row_count",
))
CONTAMINATION_RECEIPT_FIELDS = frozenset((
    "schema", "self_sha256_excluding_field", "name", "attempt", "argv",
    "start_monotonic_ns", "end_monotonic_ns", "stdout_name",
    "stdout_sha256", "stderr_name", "stderr_sha256", "contaminations",
))
TASK_RECEIPT_FIELDS = frozenset((
    "schema", "self_sha256_excluding_field", "job", "task_id", "task_sha256",
    "outer_order", "execution_receipts", "trace_sha256", "semantic_sha256",
    "work_signatures_sha256", "semantic_outcome_class", "common_success",
    "base_timed_elapsed_ns", "balanced_timed_elapsed_ns", "ratio",
    "forward_process_count", "reverse_process_count",
    "timed_rows_per_architecture",
))
PREPARE_RECEIPT_FIELDS = frozenset((
    "schema", "self_sha256_excluding_field", "prepared_utc", "design_sha256",
    "tasks_manifest_sha256", "immutable_files",
    "exact_base_binary_sha256", "exact_architecture_binary_sha256",
    "measurement_binary_sha256",
))
LAUNCH_RECEIPT_FIELDS = frozenset((
    "schema", "self_sha256_excluding_field", "started_utc", "ended_utc",
    "start_monotonic_s", "end_monotonic_s", "duration_s", "design_sha256",
    "prepare_receipt_sha256", "tasks_manifest_sha256", "task_count",
    "execution_count", "retry_count", "execution_receipts", "task_receipts",
    "thermal_reader", "thermal_source_device", "thermal_source_inode",
    "thermal_interval_sha256", "thermal_summary", "topology",
    "load_workers_stopped", "controller_core", "controller_affinity",
    "core_ticks_before", "core_ticks_after", "sibling_ticks_before",
    "sibling_ticks_after", "sibling_busy_ticks", "preflight_quiet_core_ticks",
    "preflight_quiet_sibling_ticks",
))


class TimingError(RuntimeError):
    """A substantive campaign or evidence error."""


class ContaminationError(TimingError):
    """A predeclared, receipted environmental condition eligible for retry."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds").replace(
        "+00:00", "Z")


def canonical_json(value: object) -> bytes:
    return (json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False,
        ensure_ascii=True,
    ) + "\n").encode("ascii")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_write(path: Path, data: bytes, mode: int = 0o444) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_name(path.name + ".part")
    if part.exists():
        raise TimingError("stale partial output: %s" % part)
    with part.open("xb") as stream:
        stream.write(data)
        stream.flush()
        os.fsync(stream.fileno())
    os.chmod(part, mode)
    os.replace(part, path)
    directory_fd = os.open(str(path.parent), os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)


def write_new(path: Path, data: bytes, mode: int = 0o444) -> None:
    if path.exists():
        raise TimingError("refusing to replace %s" % path)
    atomic_write(path, data, mode)


def sealed_record(schema: str, payload: Mapping[str, object]) -> Dict[str, object]:
    if not isinstance(schema, str) or not schema or "schema" in payload or \
            "self_sha256_excluding_field" in payload:
        raise TimingError("invalid sealed-record construction")
    value: Dict[str, object] = {"schema": schema, **payload}
    value["self_sha256_excluding_field"] = sha256_bytes(canonical_json(value))
    return value


def verify_sealed_record(value: object, schema: str, name: str) -> Dict[str, object]:
    if not isinstance(value, dict) or value.get("schema") != schema:
        raise TimingError("%s schema mismatch" % name)
    claimed = value.get("self_sha256_excluding_field")
    if not isinstance(claimed, str) or SHA256_RE.fullmatch(claimed) is None:
        raise TimingError("%s self-hash is malformed" % name)
    unhashed = dict(value)
    del unhashed["self_sha256_excluding_field"]
    if sha256_bytes(canonical_json(unhashed)) != claimed:
        raise TimingError("%s self-hash mismatch" % name)
    return value


def load_canonical(path: Path, name: str) -> Dict[str, object]:
    raw = path.read_bytes()
    try:
        value = json.loads(raw.decode("ascii"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise TimingError("%s is not ASCII JSON" % name) from exc
    if not isinstance(value, dict) or canonical_json(value) != raw:
        raise TimingError("%s is not a canonical JSON object" % name)
    return value


def parse_uint(text: object, name: str, maximum: int = (1 << 64) - 1) -> int:
    if not isinstance(text, str) or UINT_RE.fullmatch(text) is None:
        raise TimingError("%s is not a canonical uint" % name)
    value = int(text, 10)
    if value > maximum:
        raise TimingError("%s exceeds its integer domain" % name)
    return value


def parse_sint(
    text: object, name: str, minimum: int = -(1 << 63),
    maximum: int = (1 << 63) - 1,
) -> int:
    if not isinstance(text, str) or SINT_RE.fullmatch(text) is None:
        raise TimingError("%s is not a canonical signed integer" % name)
    value = int(text, 10)
    if not minimum <= value <= maximum:
        raise TimingError("%s exceeds its integer domain" % name)
    return value


def require_sha256(value: object, name: str) -> str:
    if not isinstance(value, str) or SHA256_RE.fullmatch(value) is None:
        raise TimingError("%s is not a lowercase SHA256" % name)
    return value


def stable_bytes(path: Path, attempts: int = 20) -> bytes:
    for _attempt in range(attempts):
        before = path.stat()
        raw = path.read_bytes()
        after = path.stat()
        if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) == \
                (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns) and \
                len(raw) == after.st_size:
            return raw
        time.sleep(0.02)
    raise TimingError("file did not stabilize while reading: %s" % path)


def command_result(
    argv: Sequence[str], cwd: Optional[Path] = None,
    environment: Optional[Mapping[str, str]] = None,
) -> Dict[str, object]:
    started = utc_now()
    begin = time.monotonic()
    result = subprocess.run(
        list(argv), cwd=str(cwd) if cwd is not None else None,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
        env=dict(environment) if environment is not None else None,
    )
    return {
        "argv": list(argv), "cwd": str(cwd) if cwd is not None else None,
        "environment": dict(environment) if environment is not None else None,
        "started_utc": started, "duration_s": time.monotonic() - begin,
        "returncode": result.returncode, "stdout": result.stdout,
        "stderr": result.stderr,
    }


def checked_text(argv: Sequence[str], cwd: Optional[Path] = None) -> str:
    result = subprocess.run(
        list(argv), cwd=str(cwd) if cwd is not None else None,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
    )
    if result.returncode != 0 or result.stderr:
        raise TimingError(
            "command failed: %r exit=%d stderr=%r" %
            (list(argv), result.returncode, result.stderr[:1000]))
    try:
        return result.stdout.decode("utf-8").strip()
    except UnicodeDecodeError as exc:
        raise TimingError("command output is not UTF-8: %r" % list(argv)) from exc


def resolve_tool(name: str) -> Path:
    value = shutil.which(name)
    if value is None:
        raise TimingError("required executable is unavailable: %s" % name)
    path = Path(value).resolve()
    if not path.is_file() or not os.access(path, os.X_OK):
        raise TimingError("required executable is not executable: %s" % path)
    return path


def parse_cpu_list(text: str) -> Tuple[int, ...]:
    if not text or text.strip() != text:
        raise TimingError("CPU-list text is not canonical")
    values: List[int] = []
    for token in text.split(","):
        if "-" in token:
            parts = token.split("-")
            if len(parts) != 2:
                raise TimingError("CPU-list range is malformed")
            low = parse_uint(parts[0], "CPU-list low", (1 << 31) - 1)
            high = parse_uint(parts[1], "CPU-list high", (1 << 31) - 1)
            if low >= high:
                raise TimingError("CPU-list range is not increasing")
            values.extend(range(low, high + 1))
        else:
            values.append(parse_uint(token, "CPU-list value", (1 << 31) - 1))
    if values != sorted(set(values)):
        raise TimingError("CPU-list is duplicate or nonascending")
    return tuple(values)


def topology_record(core: int, numa_node: int) -> Dict[str, object]:
    cpu = Path("/sys/devices/system/cpu") / ("cpu%d" % core)
    if cpu.is_symlink() or not cpu.is_dir():
        raise TimingError("timing CPU is unavailable")
    online = cpu / "online"
    if online.exists() and online.read_text(encoding="ascii").strip() != "1":
        raise TimingError("timing CPU is offline")
    siblings_text = (cpu / "topology/thread_siblings_list").read_text(
        encoding="ascii").strip()
    siblings = parse_cpu_list(siblings_text)
    if core not in siblings or len(siblings) != 2:
        raise TimingError("timing CPU does not have one expected SMT sibling")
    node_names = sorted(path.name for path in cpu.glob("node[0-9]*"))
    if node_names != ["node%d" % numa_node]:
        raise TimingError("timing CPU NUMA binding changed")
    caches: List[Tuple[int, int, str, Tuple[int, ...]]] = []
    for index in sorted((cpu / "cache").glob("index[0-9]*")):
        cache_type = (index / "type").read_text(encoding="ascii").strip()
        if cache_type not in ("Data", "Unified"):
            continue
        level = parse_uint(
            (index / "level").read_text(encoding="ascii").strip(),
            "cache level", 32)
        size_text = (index / "size").read_text(encoding="ascii").strip()
        match = re.fullmatch(r"([1-9][0-9]*)([KMG])", size_text)
        if match is None:
            raise TimingError("cache size is not canonical")
        scale = {"K": 1024, "M": 1024 ** 2, "G": 1024 ** 3}[match.group(2)]
        shared_text = (index / "shared_cpu_list").read_text(
            encoding="ascii").strip()
        shared = parse_cpu_list(shared_text)
        caches.append((level, int(match.group(1)) * scale, shared_text, shared))
    if not caches:
        raise TimingError("timing CPU has no data cache receipts")
    max_level = max(item[0] for item in caches)
    llc = max((item for item in caches if item[0] == max_level), key=lambda item: item[1])
    governor_path = cpu / "cpufreq/scaling_governor"
    preference_path = cpu / "cpufreq/energy_performance_preference"
    return {
        "core": core,
        "numa_node": numa_node,
        "thread_siblings_list": siblings_text,
        "sibling": next(value for value in siblings if value != core),
        "llc_level": max_level,
        "llc_bytes": llc[1],
        "llc_shared_cpu_list": llc[2],
        "llc_shared_cpus": list(llc[3]),
        "governor": governor_path.read_text(encoding="ascii").strip()
            if governor_path.exists() else "unavailable",
        "energy_performance_preference":
            preference_path.read_text(encoding="ascii").strip()
            if preference_path.exists() else "unavailable",
    }


def generate_tasks() -> Tuple[Dict[str, object], ...]:
    pending: List[Tuple[str, Dict[str, object]]] = []
    for K in KS:
        for width in WIDTHS:
            for schedule, seed_index, seed in SCHEDULE_SEEDS:
                for cache_state in CACHE_STATES:
                    task: Dict[str, object] = {
                        "K": K, "bb": width, "seed_index": seed_index,
                        "seed": seed, "schedule": schedule,
                        "cache_state": cache_state,
                    }
                    priority = sha256_bytes(
                        b"wirehair.wh2.degree-balanced.full-payload.order.v1\0" +
                        canonical_json(task))
                    pending.append((priority, task))
    pending.sort(key=lambda item: item[0])
    result: List[Dict[str, object]] = []
    for job, (_priority, task) in enumerate(pending):
        value = dict(task)
        value["job"] = job
        value["task_id"] = (
            "%03d.K%d.bb%d.seed%d.%s.%s" %
            (job, task["K"], task["bb"], task["seed_index"],
             task["schedule"], task["cache_state"]))
        result.append(value)
    if len(result) != 252:
        raise TimingError("internal timing task count changed")
    return tuple(result)


@dataclass(frozen=True)
class ParsedOutput:
    orientation: str
    preamble: Mapping[str, str]
    rows: Tuple[Mapping[str, str], ...]
    work_signatures: Mapping[str, Tuple[str, ...]]
    outcomes: Mapping[str, int]
    semantic_sha256: str
    stdout_sha256: str
    base_timed_elapsed_ns: int
    balanced_timed_elapsed_ns: int
    all_elapsed_ns: int
    timed_minor_faults: int
    discard_minor_faults: int
    max_minor_faults: int
    contaminations: Tuple[str, ...]


def _parse_preamble(
    line: str, expected_schema: str, expected_keys: Sequence[str],
) -> Dict[str, str]:
    prefix = "# groupedtiming: "
    if not line.startswith(prefix):
        raise TimingError("missing groupedtiming preamble")
    tokens = line[len(prefix):].split(" ")
    if any(token.count("=") != 1 for token in tokens):
        raise TimingError("groupedtiming preamble token is malformed")
    pairs = tuple(tuple(token.split("=", 1)) for token in tokens)
    if tuple(pair[0] for pair in pairs) != expected_keys:
        raise TimingError("groupedtiming preamble order/schema mismatch")
    if len({pair[0] for pair in pairs}) != len(pairs):
        raise TimingError("groupedtiming preamble contains a duplicate")
    result = dict(pairs)
    if result.get("schema") != expected_schema:
        raise TimingError("binary emitted the wrong groupedtiming schema")
    return result


def _certified_source_hits(K: int) -> int:
    if K < 2 or K > 64000:
        raise TimingError("source-hit block count is outside the codec domain")
    if K < int(ARCHITECTURE["source_hits_transition_K"]):
        return int(ARCHITECTURE["source_hits_below_transition"])
    return int(ARCHITECTURE["source_hits_at_or_above_transition"])


def _timing_staircase_rows(K: int) -> int:
    value = TIMING_STAIRCASE_ROWS.get(str(K))
    if value is None:
        raise TimingError("block count is outside the frozen timing K grid")
    return value


def _expected_preamble(
    task: Mapping[str, object], evict_bytes: int, orientation: str,
    geometry: str = "frozen",
) -> Dict[str, str]:
    if orientation not in ("forward", "reverse", "neutral"):
        raise TimingError("unknown timing orientation")
    if geometry not in ("frozen", "shared-x"):
        raise TimingError("unknown groupedtiming geometry")
    K = int(task["K"])
    width = int(task["bb"])
    staircase_rows = _timing_staircase_rows(K)
    control_degree = "1" if orientation == "reverse" else "0"
    candidate_degree = "1" if orientation == "forward" else "0"
    expected = {
        "policy": "h12-q0-grouped", "timing_scope": "solve",
        "cycles": "4", "order": INNER_ORDER, "discard_cycle": "0",
        "cycle_mode": "full", "cycle_index": "all", "N": str(K),
        "bb": str(width), "overhead": str(OVERHEAD), "loss": LOSS_TEXT,
        "seed": str(task["seed"]), "schedule": str(task["schedule"]),
        "cache_state": str(task["cache_state"]),
        "overhead_stream": "salted", "evict_bytes": str(evict_bytes),
        "eviction_prefaulted": "1",
        "control_period": str(ARCHITECTURE["period"]),
        "control_grouped_rows": str(ARCHITECTURE["grouped_rows"]),
        "control_buckets": str(ARCHITECTURE["buckets"]),
        "control_geometry": geometry,
        "control_secondary_schedule": "0",
        "control_rhs_route_expected": "streamed",
        "control_degree_balanced_staircase": control_degree,
        "control_grouped_hash_seed": "0x0",
        "control_final_h_a_columns": "0",
        "candidate_period": str(ARCHITECTURE["period"]),
        "candidate_grouped_rows": str(ARCHITECTURE["grouped_rows"]),
        "candidate_buckets": str(ARCHITECTURE["buckets"]),
        "candidate_geometry": geometry,
        "candidate_secondary_schedule": "0",
        "candidate_rhs_route_expected": "streamed",
        "candidate_degree_balanced_staircase": candidate_degree,
        "candidate_grouped_hash_seed": "0x0",
        "candidate_final_h_a_columns": "0", "gf256_rows": "10",
        "gf16_rows": "2", "mixed_residue_skew": "0",
        "mixed_residue_schedule": "constant",
        "mixed_residue_hash_seed": "0x0",
        "mixed_independent_extension_residues": "0",
        "mixed_extension_residue_seed_xor": "0x4e",
        "packet_peel_seed_xor": "0x0",
        "packet_peel_seed_table": "none",
        "binary_dense_rows_override": "0",
        "gf256_heavy_rows_override": "0", "source_hits_override": "0",
        "seed_block_bytes_override": "0",
        "packet_row_seed_multiplier": "0x1",
        "packet_row_seed_avalanche": "0",
        "odd_packet_peel_seed_xor": "0x0",
        "source_hits": str(_certified_source_hits(K)),
        "field": "1", "dense_identity_corner": "0", "heavy_family": "0",
        "staircase_rows": str(staircase_rows),
        "dense_rows": str(ARCHITECTURE["dense_rows"]),
        "heavy_rows": str(ARCHITECTURE["heavy_rows"]),
        "columns": str(
            K + staircase_rows + int(ARCHITECTURE["dense_rows"]) +
            int(ARCHITECTURE["heavy_rows"])),
        "dense_two_anchor": "1", "dense_two_anchor_phase": "0", "mix": "2",
        "payload": "distinct-packet-zero-v1",
        "payload_count": str(K + OVERHEAD),
        "payload_bytes": str((K + OVERHEAD) * width),
        "payload_alignment": "64", "payload_prefaulted": "1",
        "solve_value_storage": "owned-noinit", "solve_value_publish": "swap",
        "system_build": "outside-timer",
        "tls_reapply": "full-per-slot-outside-timer",
        "allocator_tls_state": "preflight-warmed",
    }
    if geometry == "frozen":
        expected.update({
            "control_mixed_coefficient_fingerprint":
                FROZEN_MIXED_COEFFICIENT_FINGERPRINT,
            "candidate_mixed_coefficient_fingerprint":
                FROZEN_MIXED_COEFFICIENT_FINGERPRINT,
        })
    return expected


def _semantic_arm(preamble: Mapping[str, str], arm: str) -> str:
    degree = preamble.get(arm + "_degree_balanced_staircase")
    if degree == "0":
        return "base"
    if degree == "1":
        return "balanced"
    raise TimingError("arm degree-balanced receipt is malformed")


def _outcome_class(control: int, candidate: int) -> str:
    control_success = control == 0
    candidate_success = candidate == 0
    if control_success and candidate_success:
        return "common-success"
    if control_success:
        return "control-only"
    if candidate_success:
        return "candidate-only"
    return "common-failure"


def _semantic_outcome_class(outcomes: Mapping[str, int]) -> str:
    base_success = outcomes["base"] == 0
    balanced_success = outcomes["balanced"] == 0
    if base_success and balanced_success:
        return "common-success"
    if base_success:
        return "base-only"
    if balanced_success:
        return "balanced-only"
    return "common-failure"


def parse_grouped_output(
    raw: bytes, orientation: str, task: Mapping[str, object], evict_bytes: int,
    expected_core: int,
) -> ParsedOutput:
    if orientation not in ("forward", "reverse"):
        raise TimingError("unknown timing orientation")
    if (not raw or not raw.endswith(b"\n") or b"\r" in raw or b"\0" in raw or
            b'"' in raw):
        raise TimingError("groupedtiming output is not canonical LF text")
    try:
        lines = raw.decode("ascii").splitlines()
    except UnicodeDecodeError as exc:
        raise TimingError("groupedtiming output is not ASCII") from exc
    if len(lines) != 34:
        raise TimingError("groupedtiming output does not have 34 lines")
    preamble = _parse_preamble(lines[0], "v5", PREAMBLE_V5)
    for key, value in _expected_preamble(task, evict_bytes, orientation).items():
        if preamble.get(key) != value:
            raise TimingError(
                "groupedtiming preamble mismatch %s: %r != %r" %
                (key, preamble.get(key), value))
    trace_sha256 = require_sha256(preamble.get("trace_sha256"), "trace hash")
    if require_sha256(
            preamble.get("packet_trace_sha256"), "packet trace hash") != trace_sha256:
        raise TimingError("packet trace aliases disagree")
    if FINGERPRINT_RE.fullmatch(str(preamble.get("payload_fingerprint"))) is None:
        raise TimingError("payload fingerprint is malformed")
    for arm in ("control", "candidate"):
        parse_uint(preamble.get(arm + "_attempt"), arm + " attempt", 255)
        for suffix, maximum in (("matrix_seed", (1 << 64) - 1),
                                ("peel_seed", (1 << 32) - 1)):
            text = preamble.get(arm + "_" + suffix, "")
            if HEX_RE.fullmatch(text) is None or int(text, 16) > maximum:
                raise TimingError("malformed %s %s receipt" % (arm, suffix))
        parse_sint(preamble.get("preflight_" + arm + "_result"),
                   "preflight result", -(1 << 31), (1 << 31) - 1)
        for suffix in (
                "staircase_fingerprint", "dense_fingerprint",
                "heavy_fingerprint", "mixed_coefficient_fingerprint"):
            if FINGERPRINT_RE.fullmatch(str(preamble.get(arm + "_" + suffix))) \
                    is None:
                raise TimingError("malformed graph fingerprint")
    for suffix in ("attempt", "matrix_seed", "peel_seed"):
        if preamble["control_" + suffix] != preamble["candidate_" + suffix]:
            raise TimingError("timing arms selected different attempts or seeds")
    if (preamble["control_dense_fingerprint"] !=
            preamble["candidate_dense_fingerprint"] or
            preamble["control_heavy_fingerprint"] !=
            preamble["candidate_heavy_fingerprint"] or
            preamble["control_mixed_coefficient_fingerprint"] !=
            preamble["candidate_mixed_coefficient_fingerprint"] or
            preamble["control_staircase_fingerprint"] ==
            preamble["candidate_staircase_fingerprint"]):
        raise TimingError("degree timing did not isolate the staircase graph")
    control_result = parse_sint(
        preamble["preflight_control_result"], "control preflight", 0, 1)
    candidate_result = parse_sint(
        preamble["preflight_candidate_result"], "candidate preflight", 0, 1)
    for arm, result in (("control", control_result),
                        ("candidate", candidate_result)):
        expected_route = "streamed" if result == 0 else "not-reached"
        if (preamble[arm + "_rhs_route_expected"] != "streamed" or
                preamble[arm + "_preflight_rhs_route"] != expected_route):
            raise TimingError("preflight RHS route receipt is inconsistent")
    if preamble.get("cell_class") != _outcome_class(
            control_result, candidate_result) or preamble.get("common_success") != (
                "1" if control_result == candidate_result == 0 else "0"):
        raise TimingError("timing outcome-class receipt is inconsistent")
    outcomes = {
        _semantic_arm(preamble, "control"): control_result,
        _semantic_arm(preamble, "candidate"): candidate_result,
    }
    if set(outcomes) != {"base", "balanced"}:
        raise TimingError("timing orientation does not cover both architectures")
    decoded: Dict[str, str] = {}
    for arm, result in (("control", control_result),
                        ("candidate", candidate_result)):
        fingerprint = preamble["preflight_" + arm + "_decoded_fingerprint"]
        if (result == 0 and FINGERPRINT_RE.fullmatch(fingerprint) is None) or \
                (result != 0 and fingerprint != "none"):
            raise TimingError("decoded-output fingerprint/result mismatch")
        decoded[_semantic_arm(preamble, arm)] = fingerprint
    if control_result == candidate_result == 0 and \
            decoded["base"] != decoded["balanced"]:
        raise TimingError("common-success decoded output differs")
    for field in ("staircase_rows", "dense_rows", "heavy_rows", "columns"):
        parse_uint(preamble.get(field), "dimension " + field, 65535)
    if int(preamble["columns"]) != int(task["K"]) + \
            int(preamble["staircase_rows"]) + int(preamble["dense_rows"]) + \
            int(preamble["heavy_rows"]):
        raise TimingError("precode dimensions do not sum")
    reader = csv.DictReader(io.StringIO("\n".join(lines[1:]) + "\n"))
    if tuple(reader.fieldnames or ()) != TIMING_FIELDS:
        raise TimingError("groupedtiming CSV schema mismatch")
    if any(line.count(",") != len(TIMING_FIELDS) - 1 for line in lines[1:]):
        raise TimingError("groupedtiming CSV field count changed")
    rows = tuple(dict(row) for row in reader)
    if len(rows) != 32 or any(set(row) != set(TIMING_FIELDS) for row in rows):
        raise TimingError("groupedtiming output does not have 32 rows")
    K = int(task["K"])
    width = int(task["bb"])
    signatures: Dict[str, set[Tuple[str, ...]]] = {
        "base": set(), "balanced": set()}
    contaminations: List[str] = []
    all_elapsed = 0
    architecture_elapsed = {"base": 0, "balanced": 0}
    timed_row_counts = {"base": 0, "balanced": 0}
    timed_minor = 0
    discard_minor = 0
    max_minor = 0
    for index, row in enumerate(rows):
        cycle, slot = divmod(index, 8)
        arm = "control" if INNER_ORDER[slot] == "A" else "candidate"
        semantic_arm = _semantic_arm(preamble, arm)
        preflight_result = str(
            control_result if arm == "control" else candidate_result)
        exact = {
            "N": str(K), "bb": str(width), "overhead": str(OVERHEAD),
            "schedule": str(task["schedule"]), "seed": str(task["seed"]),
            "loss": LOSS_TEXT, "cache_state": str(task["cache_state"]),
            "cycle": str(cycle), "slot": str(slot), "arm": arm,
            "period": str(ARCHITECTURE["period"]),
            "grouped_rows": str(ARCHITECTURE["grouped_rows"]),
            "buckets_requested": str(ARCHITECTURE["buckets"]),
            "mixed_geometry": "frozen",
            "degree_balanced_staircase": "1" if semantic_arm == "balanced" else "0",
            "seed_attempt": preamble[arm + "_attempt"],
            "matrix_seed": preamble[arm + "_matrix_seed"],
            "peel_seed": preamble[arm + "_peel_seed"],
            "preflight_result": preflight_result,
            "cell_class": preamble["cell_class"],
            "common_success": preamble["common_success"],
            "result": preflight_result, "outcome_stable": "1",
            "joint_source_xors": "0", "joint_marginal_xors": "0",
            "joint_marginal_copies": "0", "joint_active_deltas": "0",
            "joint_scratch_bytes": "0", "dual_source_columns": "0",
            "source_bytes": str(K * width),
            "packet_payload_bytes": str((K + OVERHEAD) * width),
            "rhs_route_expected": "streamed",
            "rhs_route_actual":
                "streamed" if preflight_result == "0" else "not-reached",
        }
        for key, value in exact.items():
            if row.get(key) != value:
                raise TimingError("row %d mismatch %s" % (index, key))
        for field in (
                "elapsed_ns", "saturated", "inactivated", "binary_def",
                "heavy_gain", "block_xors", "block_muladds", "build_ns",
                "peel_ns", "project_ns", "residual_ns", "backsub_ns",
                "joint_source_xors", "joint_marginal_xors",
                "joint_marginal_copies", "joint_active_deltas",
                "joint_scratch_bytes", "dual_source_columns", "source_bytes",
                "packet_payload_bytes", "intermediate_bytes"):
            parse_uint(row.get(field), "row %d %s" % (index, field))
        elapsed = parse_uint(row["elapsed_ns"], "elapsed_ns")
        if elapsed == 0:
            raise TimingError("timing row has zero elapsed time")
        all_elapsed += elapsed
        if cycle != 0:
            architecture_elapsed[semantic_arm] += elapsed
            timed_row_counts[semantic_arm] += 1
        minflt = parse_sint(row.get("minflt_delta"), "minor-fault delta")
        majflt = parse_sint(row.get("majflt_delta"), "major-fault delta")
        expected_fault = -1 if minflt < 0 or majflt < 0 else (
            1 if minflt or majflt else 0)
        if parse_sint(row.get("fault_contaminated"), "fault receipt", -1, 1) \
                != expected_fault:
            raise TimingError("fault receipt disagrees with signed counters")
        max_minor = max(max_minor, minflt)
        if cycle == 0:
            discard_minor += max(minflt, 0)
        else:
            timed_minor += max(minflt, 0)
        cpu_before = parse_sint(row.get("cpu_before"), "cpu_before", -1, 1 << 31)
        cpu_after = parse_sint(row.get("cpu_after"), "cpu_after", -1, 1 << 31)
        migrated = parse_sint(row.get("cpu_migrated"), "cpu_migrated", -1, 1)
        saturated = parse_uint(row.get("saturated"), "saturated", 1)
        if saturated:
            contaminations.append("row%d:saturated" % index)
        if cpu_before != expected_core or cpu_after != expected_core or migrated != 0:
            contaminations.append(
                "row%d:migration:%d->%d:%d" %
                (index, cpu_before, cpu_after, migrated))
        if minflt < 0 or minflt > MAX_MINOR_FAULTS:
            contaminations.append("row%d:minor-fault:%d" % (index, minflt))
        if majflt != 0:
            contaminations.append("row%d:major-fault:%d" % (index, majflt))
        intermediate = parse_uint(row["intermediate_bytes"], "intermediate_bytes")
        expected_intermediate = int(preamble["columns"]) * width
        if (intermediate % width or intermediate > expected_intermediate or
                (preflight_result == "0" and
                 intermediate != expected_intermediate)):
            raise TimingError("solve-value arena size is inconsistent")
        arena = parse_uint(row.get("solve_value_arena_bytes"), "arena bytes")
        eager = parse_uint(row.get("solve_value_eager_zero_bytes"), "eager-zero bytes")
        copied = parse_uint(row.get("solve_value_commit_copy_bytes"), "commit-copy bytes")
        if arena != expected_intermediate or eager != 0 or copied != 0:
            raise TimingError("solve-value arena receipt changed")
        signatures[semantic_arm].add(tuple(row[field] for field in WORK_FIELDS))
    if any(len(values) != 1 for values in signatures.values()) or \
            timed_row_counts != {"base": 12, "balanced": 12}:
        raise TimingError("architecture work or timed-row coverage is unstable")
    work_signatures = {
        arm: next(iter(signatures[arm])) for arm in ("base", "balanced")}
    graph = {}
    for cli_arm in ("control", "candidate"):
        semantic_arm = _semantic_arm(preamble, cli_arm)
        graph[semantic_arm] = {
            suffix: preamble[cli_arm + "_" + suffix] for suffix in (
                "attempt", "matrix_seed", "peel_seed",
                "staircase_fingerprint", "dense_fingerprint",
                "heavy_fingerprint", "mixed_coefficient_fingerprint")
        }
    architecture_receipt_fields = (
        "gf256_rows", "gf16_rows", "mixed_residue_skew",
        "mixed_residue_schedule", "mixed_residue_hash_seed",
        "mixed_independent_extension_residues",
        "mixed_extension_residue_seed_xor", "packet_peel_seed_xor",
        "packet_peel_seed_table", "binary_dense_rows_override",
        "gf256_heavy_rows_override", "source_hits_override",
        "seed_block_bytes_override", "packet_row_seed_multiplier",
        "packet_row_seed_avalanche", "odd_packet_peel_seed_xor",
        "source_hits", "field", "dense_identity_corner", "heavy_family",
        "dense_two_anchor", "dense_two_anchor_phase", "mix",
    )
    semantic = {
        "trace_sha256": trace_sha256,
        "payload_fingerprint": preamble["payload_fingerprint"],
        "dimensions": {field: preamble[field] for field in (
            "staircase_rows", "dense_rows", "heavy_rows", "columns")},
        "graph": graph, "decoded_fingerprints": decoded,
        "architecture_receipt": {
            field: preamble[field] for field in architecture_receipt_fields},
        "rhs_route": {
            "expected": "streamed",
            "base": "streamed" if outcomes["base"] == 0 else "not-reached",
            "balanced":
                "streamed" if outcomes["balanced"] == 0 else "not-reached",
        },
        "outcomes": outcomes,
        "semantic_outcome_class": _semantic_outcome_class(outcomes),
        "work_fields": list(WORK_FIELDS),
        "work_signatures": {
            arm: list(work_signatures[arm]) for arm in ("base", "balanced")},
        "row_count": len(rows), "timed_cycles": [1, 2, 3],
    }
    return ParsedOutput(
        orientation=orientation, preamble=preamble, rows=rows,
        work_signatures=work_signatures, outcomes=outcomes,
        semantic_sha256=sha256_bytes(canonical_json(semantic)),
        stdout_sha256=sha256_bytes(raw),
        base_timed_elapsed_ns=architecture_elapsed["base"],
        balanced_timed_elapsed_ns=architecture_elapsed["balanced"],
        all_elapsed_ns=all_elapsed, timed_minor_faults=timed_minor,
        discard_minor_faults=discard_minor, max_minor_faults=max_minor,
        contaminations=tuple(contaminations),
    )


def execution_name(task: Mapping[str, object], slot: int, orientation: str) -> str:
    return "%s.outer%d.%s.csv" % (task["task_id"], slot, orientation)


def command_for(
    design: Mapping[str, object], task: Mapping[str, object], orientation: str,
) -> List[str]:
    root = Path(str(design["root"]))
    tools = design["tools"]
    if orientation not in ("forward", "reverse"):
        raise TimingError("unknown timing orientation")
    binary = root / "frozen" / BINARY_NAMES["measurement"]
    core = int(design["core"])
    node = int(design["numa_node"])
    evict = int(design["evict_bytes"])
    return [
        str(tools["env"]["path"]), "-i", "LC_ALL=C", "TZ=UTC",
        "PATH=/usr/bin:/bin", "MALLOC_MMAP_THRESHOLD_=" + MALLOC_MMAP_THRESHOLD,
        "MALLOC_TRIM_THRESHOLD_=" + MALLOC_TRIM_THRESHOLD,
        str(tools["taskset"]["path"]), "-c", str(core),
        str(tools["numactl"]["path"]), "--physcpubind=" + str(core),
        "--membind=" + str(node), str(binary), "groupedtiming",
        "--N", str(task["K"]), "--bb", str(task["bb"]),
        "--overhead", str(OVERHEAD),
        "--control-period", str(ARCHITECTURE["period"]),
        "--control-grouped-rows", str(ARCHITECTURE["grouped_rows"]),
        "--control-buckets", str(ARCHITECTURE["buckets"]),
        "--control-geometry", str(ARCHITECTURE["geometry"]),
        "--control-degree-balanced-staircase",
        "0" if orientation == "forward" else "1",
        "--candidate-period", str(ARCHITECTURE["period"]),
        "--candidate-grouped-rows", str(ARCHITECTURE["grouped_rows"]),
        "--candidate-buckets", str(ARCHITECTURE["buckets"]),
        "--candidate-geometry", str(ARCHITECTURE["geometry"]),
        "--candidate-degree-balanced-staircase",
        "1" if orientation == "forward" else "0",
        "--evict-bytes", str(evict), "--cache-state", str(task["cache_state"]),
        "--loss", LOSS_TEXT, "--seed", str(task["seed"]),
        "--schedule", str(task["schedule"]),
    ]


def _receipt_summary(parsed: ParsedOutput) -> Dict[str, object]:
    return {
        "schema_version": "v5",
        "stdout_sha256": parsed.stdout_sha256,
        "semantic_sha256": parsed.semantic_sha256,
        "trace_sha256": parsed.preamble["trace_sha256"],
        "work_signatures_sha256": sha256_bytes(canonical_json({
            "fields": list(WORK_FIELDS),
            "values": {arm: list(parsed.work_signatures[arm])
                       for arm in ("base", "balanced")},
        })),
        "semantic_outcome_class": _semantic_outcome_class(parsed.outcomes),
        "base_timed_elapsed_ns": parsed.base_timed_elapsed_ns,
        "balanced_timed_elapsed_ns": parsed.balanced_timed_elapsed_ns,
        "all_elapsed_ns": parsed.all_elapsed_ns,
        "timed_minor_faults": parsed.timed_minor_faults,
        "discard_minor_faults": parsed.discard_minor_faults,
        "max_minor_faults": parsed.max_minor_faults,
        "row_count": len(parsed.rows), "timed_row_count": 24,
    }


def _register_cross_cache_identity(
    ledger: Dict[Tuple[object, ...], Dict[str, Dict[str, object]]],
    task: Mapping[str, object], parsed: ParsedOutput,
) -> None:
    """Require cold/warm cells to select exactly the same graph and work."""
    key = (
        task["K"], task["bb"], task["seed_index"], task["seed"],
        task["schedule"],
    )
    cache_state = str(task["cache_state"])
    if cache_state not in CACHE_STATES:
        raise TimingError("cross-cache ledger has an unknown cache state")
    identity: Dict[str, object] = {
        "trace_sha256": parsed.preamble["trace_sha256"],
        "payload_fingerprint": parsed.preamble["payload_fingerprint"],
        "work_signatures": {
            arm: list(parsed.work_signatures[arm])
            for arm in ("base", "balanced")},
        "outcomes": dict(parsed.outcomes),
        "semantic_sha256": parsed.semantic_sha256,
    }
    states = ledger.setdefault(key, {})
    if cache_state in states:
        raise TimingError("cross-cache ledger contains a duplicate cell")
    states[cache_state] = identity
    if len(states) == len(CACHE_STATES) and len({
            sha256_bytes(canonical_json(value)) for value in states.values()
    }) != 1:
        raise TimingError("cold/warm cells changed graph, trace, or work")


def _validate_cross_cache_ledger(
    ledger: Mapping[Tuple[object, ...], Mapping[str, Mapping[str, object]]],
) -> None:
    expected_pairs = len(KS) * len(WIDTHS) * len(SCHEDULE_SEEDS)
    if len(ledger) != expected_pairs or any(
            set(states) != set(CACHE_STATES) for states in ledger.values()):
        raise TimingError("cross-cache ledger is incomplete")


def _git_value(git: Path, repo: Path, *arguments: str) -> str:
    return checked_text((str(git), "-C", str(repo), *arguments))


def _dependency_hashes(binary: Path, ldd: Path) -> Tuple[str, Dict[str, str]]:
    output = checked_text((str(ldd), str(binary)))
    paths = set()
    for line in output.splitlines():
        match = re.search(r"=>\s+(/\S+)\s+\(", line)
        if match is None:
            match = re.match(r"\s*(/\S+)\s+\(", line)
        if match is not None:
            paths.add(str(Path(match.group(1)).resolve()))
    if not paths:
        raise TimingError("ldd did not identify any runtime dependencies")
    hashes = {path: sha256_file(Path(path)) for path in sorted(paths)}
    return output, hashes


def _build_one(
    label: str, commit: str, repo: Path, workspace: Path, frozen: Path,
    provenance_dir: Path, tools: Mapping[str, Path], jobs: int,
    c_compiler: Path, cxx_compiler: Path,
) -> Dict[str, object]:
    # Reuse the exact same checkout/build/HOME paths sequentially so __FILE__,
    # build-id inputs, command lengths, and tool caches cannot differ merely
    # because one arm is named "candidate".
    source = workspace / "source"
    build = workspace / "build"
    build_home = workspace / "home"
    if source.exists() or build.exists() or build_home.exists():
        raise TimingError("shared build workspace was not reset between binaries")
    add = command_result((
        str(tools["git"]), "-C", str(repo), "worktree", "add", "--detach",
        str(source), commit,
    ))
    if add["returncode"] != 0:
        raise TimingError("git worktree add failed: %r" % add["stderr"][:1000])
    try:
        head = _git_value(tools["git"], source, "rev-parse", "HEAD")
        tree = _git_value(tools["git"], source, "rev-parse", "HEAD^{tree}")
        if head != commit:
            raise TimingError("detached build worktree resolved the wrong commit")
        if _git_value(
                tools["git"], source, "status", "--porcelain",
                "--untracked-files=all"):
            raise TimingError("detached build worktree is not clean")
        generator = "Ninja" if "ninja" in tools else "Unix Makefiles"
        build_home.mkdir()
        build_environment = {
            "HOME": str(build_home), "PATH": "/usr/bin:/bin",
            "LC_ALL": "C", "LANG": "C", "TZ": "UTC",
        }
        configure_argv = (
            str(tools["cmake"]), "-S", str(source), "-B", str(build),
            "-G", generator, "-DCMAKE_BUILD_TYPE=Release",
            "-DBUILD_TESTS=ON", "-DBUILD_CODEC_V2=ON",
            "-DWIREHAIR_BUILD_BENCHMARKS=ON", "-DMARCH_NATIVE=ON",
            "-DWIREHAIR_STRICT_WARNINGS=ON", "-DWH_LTO=OFF",
            "-DWH_PGO_MODE=OFF", "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
            "-DCMAKE_C_COMPILER=" + str(c_compiler),
            "-DCMAKE_CXX_COMPILER=" + str(cxx_compiler),
        )
        configure = command_result(configure_argv, environment=build_environment)
        for stream_name in ("stdout", "stderr"):
            write_new(
                provenance_dir / (label + ".configure." + stream_name),
                configure[stream_name])
        if configure["returncode"] != 0:
            raise TimingError("%s configure failed" % label)
        build_argv = (
            str(tools["cmake"]), "--build", str(build), "--target",
            "wirehair_v2_bench", "--parallel", str(jobs),
        )
        built = command_result(build_argv, environment=build_environment)
        for stream_name in ("stdout", "stderr"):
            write_new(
                provenance_dir / (label + ".build." + stream_name),
                built[stream_name])
        if built["returncode"] != 0:
            raise TimingError("%s build failed" % label)
        if _git_value(
                tools["git"], source, "status", "--porcelain",
                "--untracked-files=all"):
            raise TimingError("out-of-tree build changed its source worktree")
        binary = build / "codec/wirehair_v2_bench"
        if not binary.is_file() or not os.access(binary, os.X_OK):
            raise TimingError("%s benchmark binary was not produced" % label)
        frozen_binary = frozen / BINARY_NAMES[label]
        shutil.copyfile(binary, frozen_binary)
        os.chmod(frozen_binary, 0o555)
        cache = build / "CMakeCache.txt"
        commands = build / "compile_commands.json"
        for source_path, suffix in (
                (cache, "CMakeCache.txt"), (commands, "compile_commands.json")):
            if not source_path.is_file():
                raise TimingError("%s build provenance is missing %s" %
                                  (label, suffix))
            shutil.copyfile(source_path, provenance_dir / (label + "." + suffix))
            os.chmod(provenance_dir / (label + "." + suffix), 0o444)
        ldd_text, dependencies = _dependency_hashes(frozen_binary, tools["ldd"])
        write_new(provenance_dir / (label + ".ldd.txt"),
                  (ldd_text + "\n").encode("utf-8"))
        compiler_version = checked_text((str(cxx_compiler), "--version"))
        compiler_target = checked_text((str(cxx_compiler), "-dumpmachine"))
        compiler_numeric = checked_text((
            str(cxx_compiler), "-dumpfullversion", "-dumpversion"))
        log_names = (
            label + ".configure.stdout", label + ".configure.stderr",
            label + ".build.stdout", label + ".build.stderr",
            label + ".CMakeCache.txt", label + ".compile_commands.json",
            label + ".ldd.txt",
        )
        payload: Dict[str, object] = {
            "label": label, "commit": commit, "tree": tree,
            "source_clean_before_and_after": True,
            "binary_name": frozen_binary.name,
            "binary_sha256": sha256_file(frozen_binary),
            "binary_size": frozen_binary.stat().st_size,
            "configure_argv": list(configure_argv),
            "configure_duration_s": configure["duration_s"],
            "build_argv": list(build_argv),
            "build_duration_s": built["duration_s"],
            "build_environment": build_environment,
            "generator": generator,
            "c_compiler": str(c_compiler),
            "c_compiler_sha256": sha256_file(c_compiler),
            "cxx_compiler": str(cxx_compiler),
            "cxx_compiler_sha256": sha256_file(cxx_compiler),
            "cxx_compiler_version": compiler_version,
            "cxx_compiler_target": compiler_target,
            "cxx_compiler_numeric_version": compiler_numeric,
            "runtime_dependency_sha256": dependencies,
            "evidence_files": {
                name: sha256_file(provenance_dir / name) for name in log_names
            },
        }
        provenance = sealed_record(
            "wirehair.wh2.degree_balanced_timing.build_provenance.v1", payload)
        path = provenance_dir / (label + ".json")
        write_new(path, canonical_json(provenance))
        return {"record": provenance, "path": path}
    finally:
        removal = subprocess.run((
            str(tools["git"]), "-C", str(repo), "worktree", "remove",
            "--force", str(source)), stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=False)
        if removal.returncode != 0 and source.exists():
            raise TimingError("could not remove temporary build worktree")
        if build.exists():
            shutil.rmtree(build)
        if build_home.exists():
            shutil.rmtree(build_home)


def _neutral_command(
    staging: Path, tools: Mapping[str, Path], task: Mapping[str, object],
    binary_key: str, core: int, numa_node: int, measurement: bool,
    geometry: str = "shared-x",
) -> List[str]:
    if geometry not in ("shared-x", "frozen") or \
            (not measurement and geometry != "shared-x"):
        raise TimingError("neutral command geometry is unsupported")
    binary = staging / "frozen" / BINARY_NAMES[binary_key]
    command = [
        str(tools["env"]), "-i", "LC_ALL=C", "TZ=UTC", "PATH=/usr/bin:/bin",
        "MALLOC_MMAP_THRESHOLD_=" + MALLOC_MMAP_THRESHOLD,
        "MALLOC_TRIM_THRESHOLD_=" + MALLOC_TRIM_THRESHOLD,
        str(tools["taskset"]), "-c", str(core), str(tools["numactl"]),
        "--physcpubind=" + str(core), "--membind=" + str(numa_node),
        str(binary), "groupedtiming", "--N", str(task["K"]), "--bb",
        str(task["bb"]), "--overhead", str(OVERHEAD), "--control-period",
        str(ARCHITECTURE["period"]), "--control-grouped-rows",
        str(ARCHITECTURE["grouped_rows"]), "--control-buckets",
        str(ARCHITECTURE["buckets"]),
    ]
    if measurement:
        command += [
            "--control-geometry", geometry,
            "--control-degree-balanced-staircase", "0"]
    command += [
        "--candidate-period", str(ARCHITECTURE["period"]),
        "--candidate-grouped-rows", str(ARCHITECTURE["grouped_rows"]),
        "--candidate-buckets", str(ARCHITECTURE["buckets"]),
    ]
    if measurement:
        command += [
            "--candidate-geometry", geometry,
            "--candidate-degree-balanced-staircase", "0"]
    command += [
        "--evict-bytes", "4096", "--cache-state", str(task["cache_state"]),
        "--loss", LOSS_TEXT, "--seed", str(task["seed"]), "--schedule",
        str(task["schedule"]),
    ]
    return command


def _parse_neutral_replay(
    raw: bytes, measurement: bool, task: Mapping[str, object],
    expected_geometry: str = "shared-x",
) -> Dict[str, object]:
    if expected_geometry not in ("shared-x", "frozen") or \
            (not measurement and expected_geometry != "shared-x"):
        raise TimingError("neutral replay geometry is unsupported")
    if (not raw or not raw.endswith(b"\n") or b"\r" in raw or b"\0" in raw or
            b'"' in raw):
        raise TimingError("neutral replay output is not canonical LF text")
    try:
        lines = raw.decode("ascii").splitlines()
    except UnicodeDecodeError as exc:
        raise TimingError("neutral replay output is not ASCII") from exc
    if len(lines) != 34:
        raise TimingError("neutral replay output line count changed")
    schema = "v5" if measurement else "v2"
    keys = PREAMBLE_V5 if measurement else LEGACY_PREAMBLE_V2
    fields = TIMING_FIELDS if measurement else LEGACY_FIELDS
    preamble = _parse_preamble(lines[0], schema, keys)
    expected = _expected_preamble(
        task, 4096, "neutral", geometry=expected_geometry)
    for key, value in expected.items():
        if key in preamble and preamble[key] != value:
            raise TimingError("neutral replay preamble mismatch: %s" % key)
    require_sha256(preamble.get("trace_sha256"), "neutral replay trace")
    for arm in ("control", "candidate"):
        parse_uint(preamble[arm + "_attempt"], "neutral replay attempt", 255)
        if (HEX_RE.fullmatch(preamble[arm + "_matrix_seed"]) is None or
                HEX_RE.fullmatch(preamble[arm + "_peel_seed"]) is None):
            raise TimingError("neutral replay seed receipt is malformed")
    for suffix in ("attempt", "matrix_seed", "peel_seed"):
        if preamble["control_" + suffix] != preamble["candidate_" + suffix]:
            raise TimingError("neutral replay arms selected different seeds")
    control_result = parse_sint(
        preamble["preflight_control_result"], "neutral control result", 0, 1)
    candidate_result = parse_sint(
        preamble["preflight_candidate_result"], "neutral candidate result", 0, 1)
    if (preamble["cell_class"] != _outcome_class(
            control_result, candidate_result) or
            preamble["common_success"] != (
                "1" if control_result == candidate_result == 0 else "0")):
        raise TimingError("neutral replay outcome receipt is inconsistent")
    reader = csv.DictReader(io.StringIO("\n".join(lines[1:]) + "\n"))
    if tuple(reader.fieldnames or ()) != fields:
        raise TimingError("neutral replay CSV schema changed")
    if any(line.count(",") != len(fields) - 1 for line in lines[1:]):
        raise TimingError("neutral replay CSV field count changed")
    rows = tuple(dict(row) for row in reader)
    if (len(rows) != 32 or any(set(row) != set(fields) for row in rows) or
            any(row.get("outcome_stable") != "1" for row in rows)):
        raise TimingError("neutral replay row coverage/outcome changed")
    deterministic_fields = (
        "N", "bb", "overhead", "schedule", "seed", "loss", "cache_state",
        "cycle", "slot", "arm", "period", "grouped_rows",
        "buckets_requested", "seed_attempt", "matrix_seed", "peel_seed",
        "preflight_result", "cell_class", "common_success", "result",
        "outcome_stable", "inactivated", "binary_def", "heavy_gain",
        "block_xors", "block_muladds", "joint_source_xors",
        "joint_marginal_xors", "joint_marginal_copies",
        "joint_active_deltas", "joint_scratch_bytes", "dual_source_columns",
        "source_bytes", "packet_payload_bytes", "intermediate_bytes",
        "solve_value_arena_bytes", "solve_value_eager_zero_bytes",
        "solve_value_commit_copy_bytes",
    )
    for index, row in enumerate(rows):
        cycle, slot = divmod(index, 8)
        arm = "control" if INNER_ORDER[slot] == "A" else "candidate"
        expected_result = str(
            control_result if arm == "control" else candidate_result)
        exact = {
            "N": str(task["K"]), "bb": str(task["bb"]),
            "overhead": str(OVERHEAD), "schedule": str(task["schedule"]),
            "seed": str(task["seed"]), "loss": LOSS_TEXT,
            "cache_state": str(task["cache_state"]), "cycle": str(cycle),
            "slot": str(slot), "arm": arm,
            "period": str(ARCHITECTURE["period"]),
            "grouped_rows": str(ARCHITECTURE["grouped_rows"]),
            "buckets_requested": str(ARCHITECTURE["buckets"]),
            "seed_attempt": preamble[arm + "_attempt"],
            "matrix_seed": preamble[arm + "_matrix_seed"],
            "peel_seed": preamble[arm + "_peel_seed"],
            "preflight_result": expected_result,
            "cell_class": preamble["cell_class"],
            "common_success": preamble["common_success"],
            "result": expected_result, "outcome_stable": "1",
            "source_bytes": str(int(task["K"]) * int(task["bb"])),
            "packet_payload_bytes": str(
                (int(task["K"]) + OVERHEAD) * int(task["bb"])),
        }
        if measurement:
            exact["mixed_geometry"] = expected_geometry
            exact["degree_balanced_staircase"] = "0"
            exact["rhs_route_expected"] = "streamed"
            exact["rhs_route_actual"] = \
                "streamed" if expected_result == "0" else "not-reached"
            exact.update({
                "joint_source_xors": "0", "joint_marginal_xors": "0",
                "joint_marginal_copies": "0", "joint_active_deltas": "0",
                "joint_scratch_bytes": "0", "dual_source_columns": "0",
            })
        if any(row.get(key) != value for key, value in exact.items()):
            raise TimingError("neutral replay row coordinate changed")
        for field in deterministic_fields:
            if field not in row or not isinstance(row[field], str):
                raise TimingError("neutral replay deterministic field is missing")
        for field in WORK_FIELDS[1:]:
            parse_uint(row[field], "neutral replay work " + field)
        for field in (
                "solve_value_arena_bytes", "solve_value_eager_zero_bytes",
                "solve_value_commit_copy_bytes"):
            parse_uint(row[field], "neutral replay " + field)
        intermediate = parse_uint(
            row["intermediate_bytes"], "neutral replay intermediate bytes")
        arena = parse_uint(
            row["solve_value_arena_bytes"], "neutral replay arena bytes")
        if (arena != intermediate or int(row["solve_value_eager_zero_bytes"]) != 0 or
                int(row["solve_value_commit_copy_bytes"]) != 0 or
                intermediate % int(task["bb"])):
            raise TimingError("neutral replay solve arena changed")
    signatures = {tuple(row[field] for field in WORK_FIELDS) for row in rows}
    if len(signatures) != 1:
        raise TimingError("neutral replay internal arms differ")
    shared = {
        "legacy_normalized_preamble": {
            key: preamble[key] for key in LEGACY_PREAMBLE_V2
            if key != "schema"},
        "deterministic_row_fields": list(deterministic_fields),
        "deterministic_rows": [
            [row[field] for field in deterministic_fields] for row in rows],
        "work_fields": list(WORK_FIELDS),
        "work_signature": list(next(iter(signatures))),
    }
    receipt: Dict[str, object] = {
        "schema_version": schema, "shared": shared,
        "shared_sha256": sha256_bytes(canonical_json(shared)),
    }
    if measurement:
        if FINGERPRINT_RE.fullmatch(preamble["payload_fingerprint"]) is None:
            raise TimingError("measurement-neutral payload hash is malformed")
        dimensions = {
            field: parse_uint(
                preamble[field], "measurement-neutral dimension " + field,
                65535)
            for field in ("staircase_rows", "dense_rows", "heavy_rows", "columns")
        }
        if dimensions["columns"] != int(task["K"]) + \
                dimensions["staircase_rows"] + dimensions["dense_rows"] + \
                dimensions["heavy_rows"]:
            raise TimingError("measurement-neutral dimensions do not sum")
        if (preamble["control_degree_balanced_staircase"] != "0" or
                preamble["candidate_degree_balanced_staircase"] != "0" or
                preamble["control_geometry"] != expected_geometry or
                preamble["candidate_geometry"] != expected_geometry or
                preamble["control_staircase_fingerprint"] !=
                preamble["candidate_staircase_fingerprint"] or
                preamble["control_dense_fingerprint"] !=
                preamble["candidate_dense_fingerprint"] or
                preamble["control_heavy_fingerprint"] !=
                preamble["candidate_heavy_fingerprint"] or
                preamble["control_mixed_coefficient_fingerprint"] !=
                preamble["candidate_mixed_coefficient_fingerprint"] or
                preamble["preflight_control_decoded_fingerprint"] !=
                preamble["preflight_candidate_decoded_fingerprint"] or
                preamble["packet_trace_sha256"] != preamble["trace_sha256"]):
            raise TimingError("measurement-neutral semantic receipt differs")
        for arm in ("control", "candidate"):
            expected_route = "streamed" if (
                control_result if arm == "control" else candidate_result
            ) == 0 else "not-reached"
            if (preamble[arm + "_secondary_schedule"] != "0" or
                    preamble[arm + "_rhs_route_expected"] != "streamed" or
                    preamble[arm + "_preflight_rhs_route"] != expected_route):
                raise TimingError(
                    "measurement-neutral RHS route receipt differs")
            for suffix in (
                    "staircase_fingerprint", "dense_fingerprint",
                    "heavy_fingerprint", "mixed_coefficient_fingerprint"):
                if FINGERPRINT_RE.fullmatch(preamble[arm + "_" + suffix]) is None:
                    raise TimingError("measurement-neutral graph hash is malformed")
            decoded = preamble["preflight_" + arm + "_decoded_fingerprint"]
            result = control_result if arm == "control" else candidate_result
            if ((result == 0 and FINGERPRINT_RE.fullmatch(decoded) is None) or
                    (result != 0 and decoded != "none")):
                raise TimingError("measurement-neutral decoded hash is malformed")
        receipt["measurement_receipts"] = {
            "staircase_fingerprint": preamble["control_staircase_fingerprint"],
            "dense_fingerprint": preamble["control_dense_fingerprint"],
            "heavy_fingerprint": preamble["control_heavy_fingerprint"],
            "mixed_coefficient_fingerprint":
                preamble["control_mixed_coefficient_fingerprint"],
            "payload_fingerprint": preamble["payload_fingerprint"],
            "decoded_fingerprint":
                preamble["preflight_control_decoded_fingerprint"],
            "packet_trace_sha256": preamble["packet_trace_sha256"],
            "geometry": expected_geometry,
            "rhs_route_expected": "streamed",
            "preflight_rhs_route":
                preamble["control_preflight_rhs_route"],
            "dimensions": {field: preamble[field] for field in (
                "staircase_rows", "dense_rows", "heavy_rows", "columns")},
        }
    return receipt


def _prepare_semantic_smoke(
    staging: Path, provenance_dir: Path, tools: Mapping[str, Path],
    core: int, numa_node: int,
) -> Dict[str, object]:
    """Prove exact-base replay and measurement neutrality under normal load."""
    task = next(
        value for value in generate_tasks()
        if value["K"] == 3200 and value["bb"] == 64 and
        value["schedule"] == "burst" and value["seed_index"] == 0 and
        value["cache_state"] == "warm")
    boundary_task = next(
        value for value in generate_tasks()
        if value["K"] == 10000 and value["bb"] == 64 and
        value["schedule"] == "burst" and value["seed_index"] == 0 and
        value["cache_state"] == "warm")
    smoke_tasks = {"primary": task, "boundary": boundary_task}
    records: Dict[str, object] = {}
    semantics = []
    for record_key, binary_key, measurement, geometry, task_key in \
            NEUTRAL_SMOKE_SPECS:
        smoke_task = smoke_tasks[task_key]
        command = _neutral_command(
            staging, tools, smoke_task, binary_key, core, numa_node, measurement,
            geometry=geometry)
        result = subprocess.run(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            timeout=300, check=False, start_new_session=True)
        if result.returncode != 0 or result.stderr:
            raise TimingError(
                "%s semantic smoke failed exit=%d stderr=%r" %
                (record_key, result.returncode, result.stderr[:1000]))
        parsed = _parse_neutral_replay(
            result.stdout, measurement, smoke_task,
            expected_geometry=geometry)
        stdout_name = "smoke." + record_key + ".csv"
        write_new(provenance_dir / stdout_name, result.stdout)
        records[record_key] = {
            "argv": command, "stdout_name": stdout_name,
            "stdout_sha256": sha256_bytes(result.stdout), **parsed,
        }
        if geometry == "shared-x":
            semantics.append(parsed["shared_sha256"])
    if len(semantics) != 3 or len(set(semantics)) != 1:
        raise TimingError(
            "exact base/architecture/measurement shared-x flag-off replay differs")

    tool_records = {
        name: {"path": str(path), "sha256": sha256_file(path)}
        for name, path in tools.items()}
    smoke_design = {
        "root": str(staging), "tools": tool_records, "core": core,
        "numa_node": numa_node, "evict_bytes": 4096,
    }
    orientations = {}
    for orientation in ("forward", "reverse"):
        command = command_for(smoke_design, task, orientation)
        result = subprocess.run(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            timeout=300, check=False, start_new_session=True)
        if result.returncode != 0 or result.stderr:
            raise TimingError("%s degree smoke failed" % orientation)
        parsed = parse_grouped_output(result.stdout, orientation, task, 4096, core)
        name = "smoke.degree." + orientation + ".csv"
        write_new(provenance_dir / name, result.stdout)
        orientations[orientation] = {
            "argv": command, "stdout_name": name,
            "stdout_sha256": parsed.stdout_sha256,
            "semantic_sha256": parsed.semantic_sha256,
            "nonpromotional_contaminations": list(parsed.contaminations),
        }
    if orientations["forward"]["semantic_sha256"] != \
            orientations["reverse"]["semantic_sha256"]:
        raise TimingError("forward/reverse degree smoke semantics differ")
    return {
        "scope": "bounded semantic replay under ordinary machine load",
        "timing_evidence": False, "task": task, "binaries": records,
        "orientations": orientations,
        "exact_shared_x_flag_off_replay": True,
        "measurement_frozen_neutral_replay": True,
        "source_hits_boundary": {
            "task": boundary_task, "expected_source_hits": 3,
            "record_key": "measurement_frozen_K10000",
        },
        "measurement_only_selector": True,
    }


def _commit_audit(repo: Path, git: Path) -> Dict[str, object]:
    architecture_parent = _git_value(git, repo, "rev-parse", ARCHITECTURE_COMMIT + "^")
    measurement_parent = _git_value(git, repo, "rev-parse", MEASUREMENT_COMMIT + "^")
    if architecture_parent != BASE_COMMIT or measurement_parent != ARCHITECTURE_COMMIT:
        raise TimingError("base/architecture/measurement parent chain changed")
    architecture_paths = tuple(_git_value(
        git, repo, "diff", "--name-only", BASE_COMMIT, ARCHITECTURE_COMMIT
    ).splitlines())
    expected_architecture_paths = (
        "codec/PrecodeTest.cpp", "codec/V2BenchCliTest.cmake",
        "codec/WirehairV2Bench.cpp", "codec/WirehairV2Precode.cpp",
        "codec/WirehairV2Precode.h", "codec/fuzz/V2ProfileFuzz.cpp",
    )
    measurement_paths = tuple(_git_value(
        git, repo, "diff", "--name-only", ARCHITECTURE_COMMIT,
        MEASUREMENT_COMMIT).splitlines())
    expected_measurement_paths = (
        "codec/V2BenchCliTest.cmake", "codec/WirehairV2Bench.cpp",
        "codec/WirehairV2Precode.cpp", "codec/WirehairV2Precode.h",
        "codec/WirehairV2Solve.cpp", "codec/WirehairV2Solve.h",
    )
    if architecture_paths != expected_architecture_paths or \
            measurement_paths != expected_measurement_paths:
        raise TimingError("architecture or measurement changed-path audit failed")
    architecture_diff = subprocess.run((
        str(git), "-C", str(repo), "diff", "--binary", BASE_COMMIT,
        ARCHITECTURE_COMMIT), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        check=False)
    if architecture_diff.returncode or architecture_diff.stderr:
        raise TimingError("cannot materialize architecture diff")
    patch_id = subprocess.run(
        (str(git), "patch-id", "--stable"), input=architecture_diff.stdout,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if patch_id.returncode or patch_id.stderr:
        raise TimingError("cannot compute architecture patch-id")
    tokens = patch_id.stdout.decode("ascii").split()
    if len(tokens) != 2 or tokens[0] != ARCHITECTURE_PATCH_ID:
        raise TimingError("degree-balanced architecture patch-id changed")
    measurement_diff = subprocess.run((
        str(git), "-C", str(repo), "diff", "--binary", ARCHITECTURE_COMMIT,
        MEASUREMENT_COMMIT), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        check=False)
    if measurement_diff.returncode or measurement_diff.stderr:
        raise TimingError("cannot materialize measurement diff")
    for descendant in (ARCHITECTURE_COMMIT, MEASUREMENT_COMMIT):
        ancestry = subprocess.run((
            str(git), "-C", str(repo), "merge-base", "--is-ancestor",
            GROUPED_MASKOPT_COMMIT, descendant), stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=False)
        if ancestry.returncode not in (0, 1) or ancestry.stdout or ancestry.stderr:
            raise TimingError("cannot audit grouped-maskopt ancestry")
        if ancestry.returncode == 0:
            raise TimingError("grouped mask optimization entered timing ancestry")
    return {
        "base_commit": BASE_COMMIT,
        "architecture_commit": ARCHITECTURE_COMMIT,
        "architecture_parent": architecture_parent,
        "architecture_patch_id_stable": ARCHITECTURE_PATCH_ID,
        "architecture_diff_sha256": sha256_bytes(architecture_diff.stdout),
        "architecture_changed_paths": list(architecture_paths),
        "measurement_commit": MEASUREMENT_COMMIT,
        "measurement_parent": measurement_parent,
        "measurement_diff_sha256": sha256_bytes(measurement_diff.stdout),
        "measurement_changed_paths": list(measurement_paths),
        "measurement_scope":
            "test-hook benchmark CLI, actual-route instrumentation, readback accessors, and CLI regressions only",
        "measurement_test_hook_only": True,
        "production_behavior_changed": False,
        "grouped_maskopt_commit": GROUPED_MASKOPT_COMMIT,
        "grouped_maskopt_in_ancestry": False,
    }


def _freeze_prior_recovery_context(
    repo: Path, git: Path, provenance_dir: Path,
) -> Dict[str, object]:
    result = subprocess.run((
        str(git), "-C", str(repo), "show",
        RECOVERY_CONTEXT_COMMIT + ":" + RECOVERY_CONTEXT_PATH,
    ), stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if result.returncode != 0 or result.stderr or \
            sha256_bytes(result.stdout) != RECOVERY_CONTEXT_SHA256:
        raise TimingError("prior recovery context is unavailable or changed")
    try:
        context = json.loads(result.stdout.decode("ascii"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise TimingError("prior recovery context is malformed") from exc
    controller = subprocess.run((
        str(git), "-C", str(repo), "show",
        RECOVERY_CONTEXT_COMMIT + ":" + RECOVERY_CONTROLLER_PATH,
    ), stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if (controller.returncode != 0 or controller.stderr or
            sha256_bytes(controller.stdout) != RECOVERY_CONTROLLER_SHA256):
        raise TimingError("prior recovery controller is unavailable or changed")
    recovery_parent = _git_value(
        git, repo, "rev-parse", RECOVERY_CANDIDATE_COMMIT + "^")
    if recovery_parent != RECOVERY_BASE_COMMIT:
        raise TimingError("prior recovery architecture parent changed")
    recovery_diff = subprocess.run((
        str(git), "-C", str(repo), "diff", "--binary",
        RECOVERY_BASE_COMMIT, RECOVERY_CANDIDATE_COMMIT,
    ), stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if recovery_diff.returncode or recovery_diff.stderr:
        raise TimingError("cannot materialize prior recovery architecture diff")
    recovery_patch_id = subprocess.run(
        (str(git), "patch-id", "--stable"), input=recovery_diff.stdout,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if recovery_patch_id.returncode or recovery_patch_id.stderr:
        raise TimingError("cannot identify prior recovery architecture patch")
    recovery_patch_tokens = recovery_patch_id.stdout.decode("ascii").split()
    if (len(recovery_patch_tokens) != 2 or
            recovery_patch_tokens[0] != ARCHITECTURE_PATCH_ID):
        raise TimingError("prior recovery architecture patch changed")
    try:
        architecture = context["architecture"]
        campaign = context["campaign"]
        verification = context["independent_verification"]
        paired = context["paired_comparison"]
        recovery = context["recovery"]
        valid = (
            context["schema"] ==
                "wirehair.wh2.degree_balanced_allk_holdout.v1.result_receipt" and
            architecture == {
                "base_commit":
                    RECOVERY_BASE_COMMIT,
                "candidate_commit":
                    RECOVERY_CANDIDATE_COMMIT,
                "candidate_kind":
                    "raw corrected degree-balanced staircase",
                "seed_fixes_applied": False,
            } and
            campaign["K_range"] == [2, 64000] and
            campaign["cells_per_arm"] == 575991 and
            campaign["total_cells"] == 1151982 and
            campaign["hard_loss"] == 0.5 and
            campaign["independent_holdout_seeds"] == 3 and
            campaign["schedules"] == [
                "burst", "adversarial", "repair-only"] and
            verification["status"] == "PASS" and
            verification["raw_cells_reduced"] == 1151982 and
            paired["both_fail"] == 533 and
            paired["both_success"] == 575080 and
            paired["repairs"] == 200 and
            paired["introductions"] == 178 and
            paired["net_failure_change"] == -22 and
            paired["common_success_xor_ratio"] ==
                "1.000013171921654" and
            paired["common_success_muladd_ratio"] ==
                "1.000016638686226" and
            paired["common_success_inactivation_ratio"] ==
                "1.000027334144132" and
            recovery["base"] == {
                "failure_rate_pct": "0.12725893286527046",
                "failures": 733, "maximum_weak_K_multiplicity": 8,
                "weak_K": 703,
            } and
            recovery["candidate"] == {
                "failure_rate_pct": "0.12343942874107408",
                "failures": 711, "maximum_weak_K_multiplicity": 5,
                "weak_K": 687,
            } and
            recovery["failure_schedule_base_to_candidate"] == {
                "adversarial": [216, 221], "burst": [279, 247],
                "repair-only": [238, 243],
            } and
            context["artifacts"] == {
                "data_manifest_sha256":
                    "8be8c802fbaacd9bc7337b4ea8d8f709402cbf6247ac8041fc66e77adf915c3e",
                "prelaunch_receipts_sha256":
                    "00215d89c9a1c2441359cb7a56bdbf80345755509ba6bf86cb06a591f6585f1c",
                "result_root":
                    "/tmp/wh2-degree-balanced-allk-holdout-bf418ca-v1",
                "validated_summary_sha256":
                    "bd34669900fc5fce12a2a9e8fa76e3ea34e64ea07d2d29c0d7fa6e497d05d215",
            } and
            context["verdict"] ==
                "favorable but mixed raw architecture evidence; full-payload timing and promotion decision remain")
    except (KeyError, TypeError) as exc:
        raise TimingError("prior recovery context schema changed") from exc
    if not valid:
        raise TimingError("prior recovery context values changed")
    frozen_name = "prior_degree_balance_recovery_context.json"
    controller_name = "prior_degree_balance_recovery_controller.py"
    write_new(provenance_dir / frozen_name, result.stdout)
    write_new(provenance_dir / controller_name, controller.stdout)
    return {
        "evidentiary_role":
            "sealed all-K raw-architecture recovery evidence at the BlockBytes=64 graph-seed profile; BlockBytes=1280/4096 timing cells are speed-only",
        "source_commit": RECOVERY_CONTEXT_COMMIT,
        "source_path": RECOVERY_CONTEXT_PATH,
        "source_sha256": RECOVERY_CONTEXT_SHA256,
        "controller_source_path": RECOVERY_CONTROLLER_PATH,
        "controller_source_sha256": RECOVERY_CONTROLLER_SHA256,
        "frozen_path": "provenance/" + frozen_name,
        "frozen_sha256": RECOVERY_CONTEXT_SHA256,
        "frozen_controller_path": "provenance/" + controller_name,
        "frozen_controller_sha256": RECOVERY_CONTROLLER_SHA256,
        "architecture_base_commit": RECOVERY_BASE_COMMIT,
        "architecture_candidate_commit": RECOVERY_CANDIDATE_COMMIT,
        "architecture_patch_id_stable": recovery_patch_tokens[0],
        "K_range": [2, 64000], "cells_per_arm": 575991,
        "total_cells": 1151982, "hard_loss": "0.5",
        "schedules": ["burst", "adversarial", "repair-only"],
        "independent_holdout_seeds": 3,
        "block_bytes_seed_profile": 64, "full_payload_solve": False,
        "base_failures": 733, "candidate_failures": 711,
        "base_weak_K": 703, "candidate_weak_K": 687,
        "base_maximum_weak_K_multiplicity": 8,
        "candidate_maximum_weak_K_multiplicity": 5,
        "repairs": 200, "introductions": 178, "net_failure_change": -22,
        "failure_schedule_base_to_candidate": {
            "adversarial": [216, 221], "burst": [279, 247],
            "repair-only": [238, 243],
        },
        "common_success_xor_ratio": "1.000013171921654",
        "common_success_muladd_ratio": "1.000016638686226",
        "common_success_inactivation_ratio": "1.000027334144132",
        "seed_fixes_applied": False,
        "timing_architecture_match": {
            **ARCHITECTURE,
            "degree_balanced_staircase_only_difference": True,
        },
        "payload_evidence_scope": {
            "64": "recovery-backed and timed",
            "1280": "speed-only until the separate cross-payload recovery holdout",
            "4096": "speed-only until the separate cross-payload recovery holdout",
            "reason":
                "MatrixSeedFromProfile hashes BlockBytes, so graph reliability does not transfer across widths",
        },
    }


def prepare_campaign(args: argparse.Namespace) -> None:
    result = Path(args.result_dir).resolve()
    repo = Path(args.repo).resolve()
    if result.exists():
        raise TimingError("result directory already exists")
    if not repo.is_dir():
        raise TimingError("repository directory does not exist")
    if (args.core < 0 or args.controller_core < 0 or
            args.numa_node < 0 or args.evict_bytes < 4096 or
            args.build_jobs <= 0):
        raise TimingError("prepare integer argument is outside its domain")
    if (args.core, args.controller_core, args.numa_node) != (
            TIMING_CORE, CONTROLLER_CORE, NUMA_NODE):
        raise TimingError(
            "prepare CPU/NUMA selection differs from the frozen campaign")
    tool_names = ("git", "cmake", "env", "taskset", "numactl", "ldd",
                  "sudo", "fuser", "true", "python3", "cat", "kill")
    tools = {name: resolve_tool(name) for name in tool_names}
    if Path(sys.executable).resolve() != tools["python3"]:
        raise TimingError("prepare is not running under the receipted Python")
    sudo_probe = subprocess.run(
        (str(tools["sudo"]), "-n", str(tools["true"])),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if sudo_probe.returncode != 0 or sudo_probe.stdout or sudo_probe.stderr:
        raise TimingError("passwordless sudo is unavailable for thermal verification")
    ninja = shutil.which("ninja")
    if ninja is not None:
        tools["ninja"] = Path(ninja).resolve()
    c_compiler = Path(args.c_compiler).resolve() if args.c_compiler else resolve_tool("cc")
    cxx_compiler = Path(args.cxx_compiler).resolve() if args.cxx_compiler else resolve_tool("c++")
    for compiler, name in ((c_compiler, "C compiler"),
                           (cxx_compiler, "C++ compiler")):
        if not compiler.is_file() or not os.access(compiler, os.X_OK):
            raise TimingError("%s is not executable" % name)
    for commit, name in (
            (BASE_COMMIT, "exact base"),
            (ARCHITECTURE_COMMIT, "exact architecture"),
            (MEASUREMENT_COMMIT, "measurement")):
        resolved = _git_value(tools["git"], repo, "rev-parse", commit + "^{commit}")
        if resolved != commit:
            raise TimingError("%s commit is unavailable or ambiguous" % name)
    topology = topology_record(args.core, args.numa_node)
    controller_topology = topology_record(args.controller_core, args.numa_node)
    sampler_topology = topology_record(SAMPLER_CORE, args.numa_node)
    if topology["sibling"] != TIMING_SIBLING:
        raise TimingError("timing SMT sibling differs from the frozen campaign")
    if (args.controller_core in topology["llc_shared_cpus"] or
            SAMPLER_CORE in topology["llc_shared_cpus"]):
        raise TimingError("controller or sampler CPU shares the timing LLC")
    if args.evict_bytes < max(2 * int(topology["llc_bytes"]), 256 * 1024 * 1024):
        raise TimingError("eviction allocation is smaller than the frozen LLC gate")
    staging = result.with_name(result.name + ".prepare.%d" % os.getpid())
    workspace = result.with_name(result.name + ".build.%d" % os.getpid())
    if staging.exists() or workspace.exists():
        raise TimingError("stale prepare workspace exists")
    staging.mkdir(parents=True)
    workspace.mkdir(parents=True)
    try:
        frozen = staging / "frozen"
        provenance_dir = staging / "provenance"
        frozen.mkdir()
        provenance_dir.mkdir()
        for directory in (
                "raw", "stderr", "exit", "receipts", "task_receipts",
                "contamination"):
            (staging / directory).mkdir()
        harness_source = Path(__file__).resolve()
        harness_frozen = frozen / "wh2_degree_balanced_timing.py"
        shutil.copyfile(harness_source, harness_frozen)
        os.chmod(harness_frozen, 0o555)
        builds = {}
        for label, commit in (
                ("exact_base", BASE_COMMIT),
                ("exact_architecture", ARCHITECTURE_COMMIT),
                ("measurement", MEASUREMENT_COMMIT)):
            builds[label] = _build_one(
                label, commit, repo, workspace, frozen, provenance_dir,
                tools, args.build_jobs, c_compiler, cxx_compiler)
        commit_audit = _commit_audit(repo, tools["git"])
        prior_recovery_context = _freeze_prior_recovery_context(
            repo, tools["git"], provenance_dir)
        prepare_smoke = _prepare_semantic_smoke(
            staging, provenance_dir, tools, args.core, args.numa_node)
        if sha256_file(harness_source) != sha256_file(harness_frozen):
            raise TimingError("timing harness changed during preparation")
        tasks = generate_tasks()
        manifest = b"".join(canonical_json(task) for task in tasks)
        write_new(staging / "tasks_manifest.jsonl", manifest)
        immutable_files: Dict[str, str] = {}
        for directory in (frozen, provenance_dir):
            for path in sorted(directory.iterdir()):
                if path.is_file():
                    immutable_files[str(path.relative_to(staging))] = sha256_file(path)
        tool_records = {
            name: {"path": str(path), "sha256": sha256_file(path)}
            for name, path in sorted(tools.items())
        }
        design_payload: Dict[str, object] = {
            "root": str(result), "base_commit": BASE_COMMIT,
            "architecture_commit": ARCHITECTURE_COMMIT,
            "measurement_commit": MEASUREMENT_COMMIT,
            "commit_audit": commit_audit,
            "prior_recovery_context": prior_recovery_context,
            "timing_scope": "full-payload decoder precode solve",
            "architecture": ARCHITECTURE,
            "K": list(KS), "bb": list(WIDTHS),
            "schedule_seeds": [list(item) for item in SCHEDULE_SEEDS],
            "cache_states": list(CACHE_STATES), "overhead": OVERHEAD,
            "loss": LOSS_TEXT, "outer_order": OUTER_ORDER,
            "inner_order": INNER_ORDER, "inner_cycles": 4,
            "discard_inner_cycle": 0, "task_count": len(tasks),
            "processes_per_task": len(OUTER_ORDER),
            "rows_per_process": 32, "timed_rows_per_process": 24,
            "timed_rows_per_architecture_per_task": 96,
            "measurement_binary_count": 1,
            "outer_semantics": {
                "A": "control=base,candidate=balanced",
                "B": "control=balanced,candidate=base",
            },
            "allocator_environment": {
                "MALLOC_MMAP_THRESHOLD_": MALLOC_MMAP_THRESHOLD,
                "MALLOC_TRIM_THRESHOLD_": MALLOC_TRIM_THRESHOLD,
                "purpose": "retain and prefault solve-output pages in each process while preserving timed eager-zero/copy work",
            },
            "minor_fault_policy": {
                "minimum": 0, "maximum": MAX_MINOR_FAULTS,
                "major_faults": 0, "all_cycles_receipted": True,
            },
            "max_environmental_attempts": MAX_ENVIRONMENTAL_ATTEMPTS,
            "core": args.core, "numa_node": args.numa_node,
            "topology": topology, "controller_core": args.controller_core,
            "controller_topology": controller_topology,
            "sampler_core": SAMPLER_CORE,
            "sampler_topology": sampler_topology,
            "evict_bytes": args.evict_bytes,
            "thermal_limits_c": {
                "cpu": MAX_CPU_TEMP_C, "dimm": MAX_DIMM_TEMP_C,
                "max_gap_s": MAX_THERMAL_GAP_S,
                "max_coverage_margin_s": MAX_THERMAL_MARGIN_S,
            },
            "fresh_only": True,
            "python_runtime": {
                "executable": str(tools["python3"]),
                "version": sys.version,
                "implementation": sys.implementation.name,
                "cache_tag": sys.implementation.cache_tag,
                "byteorder": sys.byteorder,
            },
            "thermal_reader_policy": {
                "ownership": "external-reused", "automatic_discovery": True,
                "compatible_reader_count": 1, "never_signal_external": True,
                "required_cpu": SAMPLER_CORE,
                "source_is_growing": True,
                "seal_minimal_bracketing_slice_transactionally": True,
            },
            "tasks_manifest_sha256": sha256_bytes(manifest),
            "immutable_files": immutable_files, "tools": tool_records,
            "build_provenance_sha256": {
                label: sha256_file(value["path"])
                for label, value in builds.items()
            },
            "prepare_smoke": prepare_smoke,
        }
        design = sealed_record(
            "wirehair.wh2.degree_balanced_timing.design.v1", design_payload)
        design_path = staging / "design.json"
        write_new(design_path, canonical_json(design))
        receipt = sealed_record(
            "wirehair.wh2.degree_balanced_timing.prepare_receipt.v1", {
                "prepared_utc": utc_now(), "design_sha256": sha256_file(design_path),
                "tasks_manifest_sha256": sha256_bytes(manifest),
                "immutable_files": immutable_files,
                "exact_base_binary_sha256": immutable_files[
                    "frozen/" + BINARY_NAMES["exact_base"]],
                "exact_architecture_binary_sha256": immutable_files[
                    "frozen/" + BINARY_NAMES["exact_architecture"]],
                "measurement_binary_sha256": immutable_files[
                    "frozen/" + BINARY_NAMES["measurement"]],
            })
        write_new(staging / "prepare_receipt.json", canonical_json(receipt))
        for directory in (frozen, provenance_dir):
            os.chmod(directory, 0o555)
        os.replace(staging, result)
        print(json.dumps({
            "result_dir": str(result), "task_count": len(tasks),
            "design_sha256": sha256_file(result / "design.json"),
            "manifest_sha256": sha256_file(result / "tasks_manifest.jsonl"),
            "exact_base_binary_sha256": receipt["exact_base_binary_sha256"],
            "exact_architecture_binary_sha256":
                receipt["exact_architecture_binary_sha256"],
            "measurement_binary_sha256": receipt["measurement_binary_sha256"],
        }, sort_keys=True))
    finally:
        if staging.exists():
            shutil.rmtree(staging)
        if workspace.exists():
            shutil.rmtree(workspace)


def _load_design(root: Path) -> Dict[str, object]:
    design = load_canonical(root / "design.json", "timing design")
    verify_sealed_record(
        design, "wirehair.wh2.degree_balanced_timing.design.v1", "timing design")
    if design.get("root") != str(root.resolve()):
        raise TimingError("timing root moved after preparation")
    if (design.get("base_commit") != BASE_COMMIT or
            design.get("architecture_commit") != ARCHITECTURE_COMMIT or
            design.get("measurement_commit") != MEASUREMENT_COMMIT):
        raise TimingError("timing commit identities changed")
    expected = {
        "timing_scope": "full-payload decoder precode solve",
        "architecture": ARCHITECTURE, "K": list(KS), "bb": list(WIDTHS),
        "schedule_seeds": [list(item) for item in SCHEDULE_SEEDS],
        "cache_states": list(CACHE_STATES), "overhead": OVERHEAD,
        "loss": LOSS_TEXT, "outer_order": OUTER_ORDER,
        "inner_order": INNER_ORDER, "inner_cycles": 4,
        "discard_inner_cycle": 0, "task_count": 252,
        "processes_per_task": 8, "rows_per_process": 32,
        "timed_rows_per_process": 24,
        "timed_rows_per_architecture_per_task": 96,
        "measurement_binary_count": 1,
        "outer_semantics": {
            "A": "control=base,candidate=balanced",
            "B": "control=balanced,candidate=base",
        },
        "allocator_environment": {
            "MALLOC_MMAP_THRESHOLD_": MALLOC_MMAP_THRESHOLD,
            "MALLOC_TRIM_THRESHOLD_": MALLOC_TRIM_THRESHOLD,
            "purpose": "retain and prefault solve-output pages in each process while preserving timed eager-zero/copy work",
        },
        "minor_fault_policy": {
            "minimum": 0, "maximum": MAX_MINOR_FAULTS,
            "major_faults": 0, "all_cycles_receipted": True,
        },
        "max_environmental_attempts": MAX_ENVIRONMENTAL_ATTEMPTS,
        "thermal_limits_c": {
            "cpu": MAX_CPU_TEMP_C, "dimm": MAX_DIMM_TEMP_C,
            "max_gap_s": MAX_THERMAL_GAP_S,
            "max_coverage_margin_s": MAX_THERMAL_MARGIN_S,
        },
        "fresh_only": True,
        "thermal_reader_policy": {
            "ownership": "external-reused", "automatic_discovery": True,
            "compatible_reader_count": 1, "never_signal_external": True,
            "required_cpu": SAMPLER_CORE,
            "source_is_growing": True,
            "seal_minimal_bracketing_slice_transactionally": True,
        },
    }
    for key, value in expected.items():
        if design.get(key) != value:
            raise TimingError("timing design policy changed: %s" % key)
    python_runtime = design.get("python_runtime")
    if python_runtime != {
            "executable": str(Path(sys.executable).resolve()),
            "version": sys.version,
            "implementation": sys.implementation.name,
            "cache_tag": sys.implementation.cache_tag,
            "byteorder": sys.byteorder,
    }:
        raise TimingError("timing Python runtime changed")
    for key in (
            "core", "controller_core", "sampler_core", "numa_node",
            "evict_bytes"):
        value = design.get(key)
        if not isinstance(value, int) or isinstance(value, bool) or value < 0:
            raise TimingError("timing design integer is malformed: %s" % key)
    if (design["core"], design["controller_core"], design["sampler_core"],
            design["numa_node"]) != (
                TIMING_CORE, CONTROLLER_CORE, SAMPLER_CORE, NUMA_NODE):
        raise TimingError("timing design CPU/NUMA selection changed")
    topology = design.get("topology")
    controller_topology = design.get("controller_topology")
    sampler_topology = design.get("sampler_topology")
    if (not isinstance(topology, dict) or
            not isinstance(controller_topology, dict) or
            not isinstance(sampler_topology, dict)):
        raise TimingError("timing design topology receipt is malformed")
    if (topology.get("core") != TIMING_CORE or
            topology.get("sibling") != TIMING_SIBLING or
            controller_topology.get("core") != CONTROLLER_CORE or
            sampler_topology.get("core") != SAMPLER_CORE or
            design["evict_bytes"] < 4096 or
            design["controller_core"] in topology.get("llc_shared_cpus", []) or
            design["sampler_core"] in topology.get("llc_shared_cpus", [])):
        raise TimingError("timing design isolation domain is malformed")
    smoke = design.get("prepare_smoke")
    expected_smoke_task = next(
        value for value in generate_tasks()
        if value["K"] == 3200 and value["bb"] == 64 and
        value["schedule"] == "burst" and value["seed_index"] == 0 and
        value["cache_state"] == "warm")
    expected_boundary_task = next(
        value for value in generate_tasks()
        if value["K"] == 10000 and value["bb"] == 64 and
        value["schedule"] == "burst" and value["seed_index"] == 0 and
        value["cache_state"] == "warm")
    if (not isinstance(smoke, dict) or smoke.get("timing_evidence") is not False or
            smoke.get("task") != expected_smoke_task or
            smoke.get("exact_shared_x_flag_off_replay") is not True or
            smoke.get("measurement_frozen_neutral_replay") is not True or
            smoke.get("source_hits_boundary") != {
                "task": expected_boundary_task, "expected_source_hits": 3,
                "record_key": "measurement_frozen_K10000",
            } or
            smoke.get("measurement_only_selector") is not True or
            set(smoke.get("binaries", {})) != {
                "exact_base", "exact_architecture", "measurement_shared_x",
                "measurement_frozen", "measurement_frozen_K10000"} or
            set(smoke.get("orientations", {})) != {"forward", "reverse"}):
        raise TimingError("prepare compatibility smoke receipt is malformed")
    audit = design.get("commit_audit")
    if (not isinstance(audit, dict) or audit.get("base_commit") != BASE_COMMIT or
            audit.get("architecture_commit") != ARCHITECTURE_COMMIT or
            audit.get("measurement_commit") != MEASUREMENT_COMMIT or
            audit.get("architecture_patch_id_stable") != ARCHITECTURE_PATCH_ID or
            audit.get("grouped_maskopt_in_ancestry") is not False or
            audit.get("measurement_test_hook_only") is not True or
            audit.get("production_behavior_changed") is not False):
        raise TimingError("commit separation audit is malformed")
    recovery = design.get("prior_recovery_context")
    expected_recovery = {
        "evidentiary_role":
            "sealed all-K raw-architecture recovery evidence at the BlockBytes=64 graph-seed profile; BlockBytes=1280/4096 timing cells are speed-only",
        "source_commit": RECOVERY_CONTEXT_COMMIT,
        "source_path": RECOVERY_CONTEXT_PATH,
        "source_sha256": RECOVERY_CONTEXT_SHA256,
        "controller_source_path": RECOVERY_CONTROLLER_PATH,
        "controller_source_sha256": RECOVERY_CONTROLLER_SHA256,
        "frozen_path":
            "provenance/prior_degree_balance_recovery_context.json",
        "frozen_sha256": RECOVERY_CONTEXT_SHA256,
        "frozen_controller_path":
            "provenance/prior_degree_balance_recovery_controller.py",
        "frozen_controller_sha256": RECOVERY_CONTROLLER_SHA256,
        "architecture_base_commit": RECOVERY_BASE_COMMIT,
        "architecture_candidate_commit": RECOVERY_CANDIDATE_COMMIT,
        "architecture_patch_id_stable": ARCHITECTURE_PATCH_ID,
        "K_range": [2, 64000], "cells_per_arm": 575991,
        "total_cells": 1151982, "hard_loss": "0.5",
        "schedules": ["burst", "adversarial", "repair-only"],
        "independent_holdout_seeds": 3,
        "block_bytes_seed_profile": 64, "full_payload_solve": False,
        "base_failures": 733, "candidate_failures": 711,
        "base_weak_K": 703, "candidate_weak_K": 687,
        "base_maximum_weak_K_multiplicity": 8,
        "candidate_maximum_weak_K_multiplicity": 5,
        "repairs": 200, "introductions": 178, "net_failure_change": -22,
        "failure_schedule_base_to_candidate": {
            "adversarial": [216, 221], "burst": [279, 247],
            "repair-only": [238, 243],
        },
        "common_success_xor_ratio": "1.000013171921654",
        "common_success_muladd_ratio": "1.000016638686226",
        "common_success_inactivation_ratio": "1.000027334144132",
        "seed_fixes_applied": False,
        "timing_architecture_match": {
            **ARCHITECTURE,
            "degree_balanced_staircase_only_difference": True,
        },
        "payload_evidence_scope": {
            "64": "recovery-backed and timed",
            "1280":
                "speed-only until the separate cross-payload recovery holdout",
            "4096":
                "speed-only until the separate cross-payload recovery holdout",
            "reason":
                "MatrixSeedFromProfile hashes BlockBytes, so graph reliability does not transfer across widths",
        },
    }
    if not isinstance(recovery, dict) or recovery != expected_recovery:
        raise TimingError("prior recovery context receipt is malformed")
    return design


def _load_tasks(root: Path, design: Mapping[str, object]) -> Tuple[Dict[str, object], ...]:
    raw = stable_bytes(root / "tasks_manifest.jsonl")
    if sha256_bytes(raw) != design.get("tasks_manifest_sha256"):
        raise TimingError("task manifest hash mismatch")
    rows: List[Dict[str, object]] = []
    for line in raw.splitlines(keepends=True):
        try:
            value = json.loads(line.decode("ascii"))
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise TimingError("task manifest line is malformed") from exc
        if not isinstance(value, dict) or canonical_json(value) != line:
            raise TimingError("task manifest line is noncanonical")
        if value.get("job") != len(rows):
            raise TimingError("task manifest jobs are not contiguous")
        rows.append(value)
    if tuple(rows) != generate_tasks() or len(rows) != design.get("task_count"):
        raise TimingError("task manifest is not the frozen Cartesian grid")
    return tuple(rows)


def _verify_immutable(root: Path, design: Mapping[str, object]) -> None:
    expected = design.get("immutable_files")
    if not isinstance(expected, dict):
        raise TimingError("immutable-file ledger is malformed")
    actual = set()
    for dirname in ("frozen", "provenance"):
        directory = root / dirname
        for path in directory.iterdir():
            if path.is_symlink() or not path.is_file():
                raise TimingError("immutable directory contains a non-regular file")
            actual.add(str(path.relative_to(root)))
    if actual != set(expected):
        raise TimingError("immutable campaign inventory changed")
    for relative, digest in expected.items():
        require_sha256(digest, "immutable-file hash")
        if sha256_file(root / relative) != digest:
            raise TimingError("immutable campaign input changed: %s" % relative)
    active_harness = Path(__file__).resolve()
    if sha256_file(active_harness) != expected.get(
            "frozen/wh2_degree_balanced_timing.py"):
        raise TimingError("active timing harness is not the frozen harness")
    for label, commit in (
            ("exact_base", BASE_COMMIT),
            ("exact_architecture", ARCHITECTURE_COMMIT),
            ("measurement", MEASUREMENT_COMMIT)):
        provenance = load_canonical(
            root / "provenance" / (label + ".json"), label + " build provenance")
        verify_sealed_record(
            provenance, "wirehair.wh2.degree_balanced_timing.build_provenance.v1",
            label + " build provenance")
        if provenance.get("label") != label or provenance.get("commit") != commit:
            raise TimingError("build provenance identity mismatch")
        binary = root / "frozen" / BINARY_NAMES[label]
        if provenance.get("binary_sha256") != sha256_file(binary):
            raise TimingError("frozen binary/provenance binding mismatch")
        evidence_files = provenance.get("evidence_files")
        if not isinstance(evidence_files, dict) or not evidence_files:
            raise TimingError("build provenance evidence ledger is missing")
        for name, digest in evidence_files.items():
            if (not isinstance(name, str) or "/" in name or
                    sha256_file(root / "provenance" / name) != digest):
                raise TimingError("build provenance evidence hash mismatch")
        dependencies = provenance.get("runtime_dependency_sha256")
        if not isinstance(dependencies, dict) or not dependencies:
            raise TimingError("runtime dependency ledger is missing")
        for path, digest in dependencies.items():
            if not isinstance(path, str) or sha256_file(Path(path)) != digest:
                raise TimingError("runtime dependency changed: %s" % path)
    tools = design.get("tools")
    if not isinstance(tools, dict):
        raise TimingError("tool ledger is malformed")
    for name, record in tools.items():
        if not isinstance(record, dict) or sha256_file(Path(str(record.get("path")))) \
                != record.get("sha256"):
            raise TimingError("campaign tool changed: %s" % name)
    if ("python3" not in tools or
            Path(sys.executable).resolve() !=
            Path(str(tools["python3"]["path"])).resolve()):
        raise TimingError("active Python is not the frozen campaign interpreter")
    smoke = design["prepare_smoke"]
    smoke_task = smoke.get("task")
    smoke_records = smoke.get("binaries")
    if not isinstance(smoke_task, dict) or not isinstance(smoke_records, dict):
        raise TimingError("prepare smoke ledger is malformed")
    smoke_tools = {
        name: Path(str(record["path"])) for name, record in tools.items()}
    boundary = smoke.get("source_hits_boundary")
    if not isinstance(boundary, dict) or not isinstance(
            boundary.get("task"), dict):
        raise TimingError("prepare source-hit boundary ledger is malformed")
    smoke_tasks = {"primary": smoke_task, "boundary": boundary["task"]}
    semantic_hashes = []
    smoke_staging_roots = set()
    for record_key, binary_key, measurement, geometry, task_key in \
            NEUTRAL_SMOKE_SPECS:
        replay_task = smoke_tasks[task_key]
        record = smoke_records.get(record_key)
        expected_name = "smoke." + record_key + ".csv"
        if not isinstance(record, dict) or record.get("stdout_name") != expected_name:
            raise TimingError("prepare smoke binary ledger is malformed")
        argv = record.get("argv")
        if not isinstance(argv, list) or any(
                not isinstance(value, str) for value in argv):
            raise TimingError("prepare smoke argv is malformed")
        binary_indices = [
            index for index, value in enumerate(argv)
            if Path(value).name == BINARY_NAMES[binary_key]]
        if len(binary_indices) != 1:
            raise TimingError("prepare smoke binary argv binding is missing")
        smoke_staging = Path(argv[binary_indices[0]]).parent.parent
        smoke_staging_roots.add(str(smoke_staging))
        if argv != _neutral_command(
                smoke_staging, smoke_tools, replay_task, binary_key,
                int(design["core"]), int(design["numa_node"]), measurement,
                geometry=geometry):
            raise TimingError("prepare smoke argv does not replay")
        raw = (root / "provenance" / expected_name).read_bytes()
        parsed = _parse_neutral_replay(
            raw, measurement, replay_task, expected_geometry=geometry)
        expected_record = {"stdout_sha256": sha256_bytes(raw), **parsed}
        if (set(record) != {"argv", "stdout_name", *expected_record} or
                any(record.get(key) != value
                    for key, value in expected_record.items())):
            raise TimingError("prepare smoke receipt does not replay")
        if geometry == "shared-x":
            semantic_hashes.append(parsed["shared_sha256"])
    if (len(semantic_hashes) != 3 or len(set(semantic_hashes)) != 1 or
            len(smoke_staging_roots) != 1):
        raise TimingError("prepare smoke cross-binary identity changed")
    orientation_records = smoke.get("orientations")
    if not isinstance(orientation_records, dict):
        raise TimingError("prepare orientation smoke ledger is malformed")
    orientation_semantics = []
    for orientation in ("forward", "reverse"):
        record = orientation_records.get(orientation)
        expected_name = "smoke.degree." + orientation + ".csv"
        if not isinstance(record, dict) or record.get("stdout_name") != expected_name:
            raise TimingError("prepare orientation smoke receipt is malformed")
        argv = record.get("argv")
        if not isinstance(argv, list) or any(
                not isinstance(value, str) for value in argv):
            raise TimingError("prepare orientation argv is malformed")
        binary_indices = [
            index for index, value in enumerate(argv)
            if Path(value).name == BINARY_NAMES["measurement"]]
        if len(binary_indices) != 1:
            raise TimingError("prepare orientation binary binding is missing")
        smoke_staging = Path(argv[binary_indices[0]]).parent.parent
        smoke_staging_roots.add(str(smoke_staging))
        smoke_design = {
            "root": str(smoke_staging), "tools": tools,
            "core": design["core"], "numa_node": design["numa_node"],
            "evict_bytes": 4096,
        }
        if argv != command_for(smoke_design, smoke_task, orientation):
            raise TimingError("prepare orientation argv does not replay")
        raw = (root / "provenance" / expected_name).read_bytes()
        parsed = parse_grouped_output(
            raw, orientation, smoke_task, 4096, int(design["core"]))
        if (set(record) != {
                "argv", "stdout_name", "stdout_sha256", "semantic_sha256",
                "nonpromotional_contaminations"} or
                record.get("stdout_sha256") != parsed.stdout_sha256 or
                record.get("semantic_sha256") != parsed.semantic_sha256 or
                record.get("nonpromotional_contaminations") !=
                list(parsed.contaminations)):
            raise TimingError("prepare orientation smoke does not replay")
        orientation_semantics.append(parsed.semantic_sha256)
    if len(set(orientation_semantics)) != 1 or len(smoke_staging_roots) != 1:
        raise TimingError("prepare orientation semantics changed")


def _validate_prepare_receipt(
    root: Path, design: Mapping[str, object],
) -> Dict[str, object]:
    prepare = load_canonical(root / "prepare_receipt.json", "prepare receipt")
    verify_sealed_record(
        prepare, "wirehair.wh2.degree_balanced_timing.prepare_receipt.v1",
        "prepare receipt")
    if set(prepare) != PREPARE_RECEIPT_FIELDS:
        raise TimingError("prepare receipt fields changed")
    if (prepare.get("design_sha256") != sha256_file(root / "design.json") or
            prepare.get("tasks_manifest_sha256") !=
            sha256_file(root / "tasks_manifest.jsonl") or
            prepare.get("immutable_files") != design.get("immutable_files") or
            prepare.get("exact_base_binary_sha256") !=
            design["immutable_files"].get(
                "frozen/" + BINARY_NAMES["exact_base"]) or
            prepare.get("exact_architecture_binary_sha256") !=
            design["immutable_files"].get(
                "frozen/" + BINARY_NAMES["exact_architecture"]) or
            prepare.get("measurement_binary_sha256") !=
            design["immutable_files"].get(
                "frozen/" + BINARY_NAMES["measurement"]) or
            not isinstance(prepare.get("prepared_utc"), str) or
            UTC_RE.fullmatch(prepare["prepared_utc"]) is None):
        raise TimingError("prepare receipt does not bind the frozen campaign")
    return prepare


def _verify_directory_inventory(root: Path) -> None:
    expected_directories = {
        "frozen", "provenance", "raw", "stderr", "exit", "receipts",
        "task_receipts", "contamination",
    }
    actual_directories = {
        str(path.relative_to(root)) for path in root.rglob("*")
        if path.is_dir() and not path.is_symlink()
    }
    if actual_directories != expected_directories or any(
            path.is_symlink() for path in root.rglob("*")):
        raise TimingError("campaign directory inventory changed")


def _fresh_output_preflight(
    root: Path, design: Mapping[str, object],
) -> None:
    _verify_directory_inventory(root)
    for dirname in (
            "raw", "stderr", "exit", "receipts", "task_receipts",
            "contamination"):
        directory = root / dirname
        entries = list(directory.iterdir())
        if entries:
            raise TimingError("fresh-only launch found %s in %s" %
                              (entries[0].name, dirname))
    for name in (
            "launch_receipt.json", "thermal_interval.csv",
            "validated_summary.json", "data_manifest.json",
            "data_manifest.sha256"):
        if (root / name).exists() or (root / (name + ".part")).exists():
            raise TimingError("fresh-only launch found stale %s" % name)
    expected_files = {
        "design.json", "prepare_receipt.json", "tasks_manifest.jsonl",
        *[str(relative) for relative in design["immutable_files"]],
    }
    actual_files = {
        str(path.relative_to(root)) for path in root.rglob("*") if path.is_file()
    }
    if actual_files != expected_files:
        raise TimingError("prepared campaign file inventory changed")


def cpu_ticks(cpu: int) -> Tuple[int, ...]:
    prefix = "cpu%d " % cpu
    for line in Path("/proc/stat").read_text(encoding="ascii").splitlines():
        if line.startswith(prefix):
            return tuple(int(value) for value in line.split()[1:])
    raise TimingError("CPU is missing from /proc/stat")


def busy_ticks(ticks: Sequence[int]) -> int:
    idle = ticks[3] + (ticks[4] if len(ticks) > 4 else 0)
    return sum(ticks) - idle


def _filler_pids() -> Tuple[int, ...]:
    result = []
    signatures = (b"wirehair_load_fillers.sh", b"while :; do :; done")
    for path in Path("/proc").glob("[0-9]*/cmdline"):
        try:
            command = path.read_bytes().replace(b"\0", b" ")
        except OSError:
            continue
        if any(signature in command for signature in signatures):
            result.append(int(path.parent.name))
    return tuple(sorted(result))


def _receipt_tool_path(
    tools: Mapping[str, object], name: str,
) -> str:
    record = tools.get(name)
    if isinstance(record, Path):
        return str(record)
    if isinstance(record, dict) and isinstance(record.get("path"), str):
        return record["path"]
    raise TimingError("tool receipt is missing: %s" % name)


def _pid_alive(pid: int, tools: Mapping[str, object]) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        result = subprocess.run(
            (_receipt_tool_path(tools, "sudo"), "-n",
             _receipt_tool_path(tools, "kill"), "-0", str(pid)),
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
        if result.returncode == 0:
            return True
        if Path("/proc/%d" % pid).exists():
            raise TimingError("cannot determine protected sampler liveness")
        return False
    return True


def _proc_bytes(
    pid: int, name: str, tools: Mapping[str, object],
) -> bytes:
    path = Path("/proc/%d/%s" % (pid, name))
    try:
        return path.read_bytes()
    except PermissionError:
        result = subprocess.run(
            (_receipt_tool_path(tools, "sudo"), "-n",
             _receipt_tool_path(tools, "cat"), str(path)),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=False)
        if result.returncode != 0 or result.stderr:
            raise TimingError("cannot inspect protected sampler process")
        return result.stdout


def _process_start_ticks(pid: int, tools: Mapping[str, object]) -> int:
    raw = _proc_bytes(pid, "stat", tools)
    try:
        return int(raw.rsplit(b") ", 1)[1].split()[19])
    except (IndexError, ValueError) as exc:
        raise TimingError("sampler process-start receipt is malformed") from exc


def _thermal_reader(
    design: Mapping[str, object], pid_file: Path,
) -> Dict[str, object]:
    if pid_file.is_symlink() or not pid_file.is_file():
        raise TimingError("thermal sampler PID file is not a regular file")
    text = pid_file.read_text(encoding="ascii")
    if not re.fullmatch(r"[1-9][0-9]*\n", text):
        raise TimingError("thermal sampler PID file is not canonical")
    pid = int(text)
    tools = design["tools"]
    if not _pid_alive(pid, tools):
        raise TimingError("thermal sampler is not alive")
    process_start_ticks = _process_start_ticks(pid, tools)
    result = subprocess.run((
        str(tools["sudo"]["path"]), "-n", str(tools["fuser"]["path"]),
        "/dev/i2c-1", "/dev/i2c-2"), stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, check=False)
    if result.returncode != 0:
        raise TimingError("cannot verify sole I2C reader with sudo fuser")
    reader_pids = sorted(set(int(value) for value in re.findall(rb"[0-9]+", result.stdout)))
    if reader_pids != [pid]:
        raise TimingError("thermal sampler is not the sole I2C reader")
    status = _proc_bytes(pid, "status", tools).decode("ascii")
    match = re.search(r"^Cpus_allowed_list:\s*(\S+)\s*$", status, re.MULTILINE)
    if match is None:
        raise TimingError("thermal sampler affinity is unavailable")
    allowed = parse_cpu_list(match.group(1))
    if allowed != (int(design["sampler_core"]),):
        raise TimingError("thermal sampler is not pinned to frozen CPU 127")
    if allowed[0] == int(design["controller_core"]):
        raise TimingError("thermal sampler shares the controller CPU")
    llc = set(int(value) for value in design["topology"]["llc_shared_cpus"])
    if set(allowed) & llc:
        raise TimingError("thermal sampler shares the timing LLC")
    raw_cmdline = _proc_bytes(pid, "cmdline", tools)
    tokens = [token.decode("utf-8") for token in raw_cmdline.split(b"\0") if token]
    script_tokens = [token for token in tokens
                     if Path(token).name == "wirehair_expo_thermal_sampler.py"]
    if (len(script_tokens) != 1 or tokens.count("--csv") != 1 or
            tokens.count("--pid-file") != 1):
        raise TimingError("I2C reader is not the expected thermal sampler")
    raw_script = Path(script_tokens[0])
    if (not raw_script.is_absolute() or raw_script.is_symlink() or
            not raw_script.is_file()):
        raise TimingError("thermal sampler script is unavailable")
    script = raw_script.resolve()
    try:
        csv_index = tokens.index("--csv")
        pid_index = tokens.index("--pid-file")
        raw_csv_path = Path(tokens[csv_index + 1])
        raw_pid_file = Path(tokens[pid_index + 1])
        if (not raw_csv_path.is_absolute() or not raw_pid_file.is_absolute() or
                raw_csv_path.is_symlink() or not raw_csv_path.is_file() or
                raw_pid_file.is_symlink() or not raw_pid_file.is_file()):
            raise TimingError("thermal sampler source path is not a regular file")
        csv_path = str(raw_csv_path.resolve())
        command_pid_file = str(raw_pid_file.resolve())
    except (ValueError, IndexError) as exc:
        raise TimingError("thermal sampler command lacks source paths") from exc
    if command_pid_file != str(pid_file.resolve()):
        raise TimingError("thermal sampler PID-file command binding changed")
    if (_process_start_ticks(pid, tools) != process_start_ticks or
            pid_file.read_text(encoding="ascii") != text or
            not _pid_alive(pid, tools)):
        raise TimingError("thermal sampler process identity changed during inspection")
    return {
        "ownership": "external-reused", "pid": pid,
        "process_start_ticks": process_start_ticks,
        "cpus_allowed_list": match.group(1),
        "cmdline_sha256": sha256_bytes(raw_cmdline),
        "script_path": str(script),
        "script_sha256": sha256_bytes(stable_bytes(script)),
        "csv_path": csv_path, "pid_file": command_pid_file,
        "unique_i2c_reader_count": len(reader_pids),
        "external_sampler_signalled": False,
    }


def _discover_external_thermal_reader(
    design: Mapping[str, object], explicit_csv: Optional[str],
    explicit_pid_file: Optional[str],
) -> Tuple[Dict[str, object], Path, Path]:
    tools = design["tools"]
    result = subprocess.run((
        str(tools["sudo"]["path"]), "-n", str(tools["fuser"]["path"]),
        "/dev/i2c-1", "/dev/i2c-2"), stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, check=False)
    if result.returncode != 0:
        raise TimingError("automatic thermal reader discovery failed")
    pids = sorted(set(int(value) for value in re.findall(rb"[0-9]+", result.stdout)))
    if len(pids) != 1:
        raise TimingError("thermal discovery requires exactly one I2C reader")
    raw = _proc_bytes(pids[0], "cmdline", tools)
    tokens = [token.decode("utf-8") for token in raw.split(b"\0") if token]
    try:
        pid_file = Path(tokens[tokens.index("--pid-file") + 1]).resolve()
        csv_path = Path(tokens[tokens.index("--csv") + 1]).resolve()
    except (ValueError, IndexError) as exc:
        raise TimingError("sole I2C reader is not a compatible sampler") from exc
    if (explicit_csv is not None and Path(explicit_csv).resolve() != csv_path) or \
            (explicit_pid_file is not None and
             Path(explicit_pid_file).resolve() != pid_file):
        raise TimingError("explicit thermal paths disagree with automatic discovery")
    reader = _thermal_reader(design, pid_file)
    if reader["pid"] != pids[0] or reader["csv_path"] != str(csv_path):
        raise TimingError("thermal discovery identity changed during preflight")
    return reader, csv_path, pid_file


def _parse_thermal(raw: bytes) -> Tuple[Tuple[Dict[str, str], ...], Tuple[bytes, ...]]:
    # The frozen sampler uses csv.DictWriter's canonical Excel line ending.
    if (not raw or not raw.endswith(b"\r\n") or b"\0" in raw or b'"' in raw):
        raise TimingError("thermal CSV is not canonical CRLF text")
    lines = tuple(raw.splitlines(keepends=True))
    if (any(not line.endswith(b"\r\n") for line in lines) or
            any(len(line[:-2].split(b",")) != len(THERMAL_FIELDS)
                for line in lines)):
        raise TimingError("thermal CSV field count changed")
    try:
        header = lines[0].decode("ascii")[:-2].split(",")
        reader = csv.DictReader(io.StringIO(raw.decode("ascii")))
    except UnicodeDecodeError as exc:
        raise TimingError("thermal CSV is not ASCII") from exc
    if tuple(header) != THERMAL_FIELDS or tuple(reader.fieldnames or ()) != THERMAL_FIELDS:
        raise TimingError("thermal CSV schema mismatch")
    rows = tuple(dict(row) for row in reader)
    if len(lines) != len(rows) + 1:
        raise TimingError("thermal CSV row accounting mismatch")
    previous = -math.inf
    for index, row in enumerate(rows):
        if re.fullmatch(
                r"[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:"
                r"[0-9]{2}\.[0-9]{3}Z", row["utc"]) is None:
            raise TimingError("thermal UTC timestamp is malformed")
        try:
            monotonic = float(row["monotonic_s"])
        except ValueError as exc:
            raise TimingError("thermal monotonic timestamp is malformed") from exc
        if not math.isfinite(monotonic) or monotonic <= previous:
            raise TimingError("thermal monotonic timestamps are not increasing")
        previous = monotonic
        for field in ("cpu_tctl_c", *DIMM_FIELDS):
            try:
                value = float(row[field])
            except ValueError as exc:
                raise TimingError("thermal temperature is missing") from exc
            if not math.isfinite(value):
                raise TimingError("thermal temperature is nonfinite")
            if ((field == "cpu_tctl_c" and not 0.0 < value < 130.0) or
                    (field != "cpu_tctl_c" and not -40.0 < value < 130.0)):
                raise TimingError("thermal temperature is implausible")
        for field in ("dimm_read_errors", "edac_ce", "edac_ue"):
            parse_uint(row[field], "thermal " + field)
    return rows, lines


def _latest_thermal_row(path: Path) -> Dict[str, str]:
    rows, _lines, _raw = _stable_thermal_snapshot(path)
    if not rows:
        raise TimingError("thermal CSV has no samples")
    return rows[-1]


def _stable_thermal_snapshot(
    path: Path, timeout_s: float = 1.0,
) -> Tuple[Tuple[Dict[str, str], ...], Tuple[bytes, ...], bytes]:
    """Read around a concurrent append without accepting malformed evidence."""
    deadline = time.monotonic() + timeout_s
    last_error: Optional[BaseException] = None
    while True:
        try:
            raw = stable_bytes(path, attempts=3)
            rows, lines = _parse_thermal(raw)
            return rows, lines, raw
        except (OSError, TimingError) as exc:
            last_error = exc
        if time.monotonic() >= deadline:
            raise TimingError(
                "thermal CSV did not yield a canonical stable snapshot: %s" %
                last_error) from last_error
        time.sleep(0.02)


def validate_sealed_thermal_interval(
    raw: bytes, start_s: float, end_s: float,
) -> Dict[str, object]:
    """Replay the minimal bracketing interval and all health/cadence gates."""
    if (not isinstance(start_s, (int, float)) or isinstance(start_s, bool) or
            not isinstance(end_s, (int, float)) or isinstance(end_s, bool) or
            not math.isfinite(float(start_s)) or not math.isfinite(float(end_s)) or
            end_s <= start_s):
        raise TimingError("thermal campaign interval is malformed")
    rows, _lines = _parse_thermal(raw)
    if not rows:
        raise TimingError("thermal interval has no samples")
    times = [float(row["monotonic_s"]) for row in rows]
    if (times[0] > start_s or times[-1] < end_s or
            (len(times) > 1 and times[1] <= start_s) or
            (len(times) > 1 and times[-2] >= end_s)):
        raise TimingError("thermal interval is not the minimal campaign bracket")
    gaps = [right - left for left, right in zip(times, times[1:])]
    start_margin = start_s - times[0]
    end_margin = times[-1] - end_s
    if (start_margin < 0 or end_margin < 0 or
            start_margin > MAX_THERMAL_MARGIN_S or
            end_margin > MAX_THERMAL_MARGIN_S or
            (gaps and max(gaps) > MAX_THERMAL_GAP_S)):
        raise TimingError("thermal interval has a coverage or cadence gap")
    if any(parse_uint(row["dimm_read_errors"], "DIMM errors") != 0
           for row in rows):
        raise TimingError("thermal interval contains a DIMM read error")
    ce = [parse_uint(row["edac_ce"], "EDAC CE") for row in rows]
    ue = [parse_uint(row["edac_ue"], "EDAC UE") for row in rows]
    if any(value != ce[0] for value in ce) or \
            any(value != ue[0] for value in ue):
        raise TimingError("EDAC counters changed during timing")
    cpu_max = max(float(row["cpu_tctl_c"]) for row in rows)
    dimm_maxima = {
        field: max(float(row[field]) for row in rows) for field in DIMM_FIELDS
    }
    if cpu_max > MAX_CPU_TEMP_C or max(dimm_maxima.values()) > MAX_DIMM_TEMP_C:
        raise TimingError("thermal interval exceeded the frozen temperature gate")
    return {
        "sample_count": len(rows), "start_margin_s": start_margin,
        "end_margin_s": end_margin, "max_gap_s": max(gaps) if gaps else 0.0,
        "cpu_max_c": cpu_max, "dimm_max_c": dimm_maxima,
        "dimm_read_errors": 0, "edac_ce_delta": ce[-1] - ce[0],
        "edac_ue_delta": ue[-1] - ue[0],
    }


def collect_thermal_interval(
    path: Path, start_s: float, end_s: float,
) -> Tuple[bytes, Dict[str, object]]:
    deadline = time.monotonic() + 10.0
    while True:
        rows, lines, _raw = _stable_thermal_snapshot(path)
        if rows and float(rows[-1]["monotonic_s"]) >= end_s:
            break
        if time.monotonic() >= deadline:
            raise TimingError("thermal sampler did not cover campaign end")
        time.sleep(0.1)
    before = [index for index, row in enumerate(rows)
              if float(row["monotonic_s"]) <= start_s]
    after = [index for index, row in enumerate(rows)
             if float(row["monotonic_s"]) >= end_s]
    if not before or not after:
        raise TimingError("thermal interval does not bracket the campaign")
    first = before[-1]
    last = after[0]
    selected_raw = lines[0] + b"".join(lines[first + 1:last + 2])
    return selected_raw, validate_sealed_thermal_interval(
        selected_raw, start_s, end_s)


def _validate_topology_again(design: Mapping[str, object]) -> Dict[str, object]:
    if (design.get("core"), design.get("controller_core"),
            design.get("sampler_core"), design.get("numa_node")) != (
                TIMING_CORE, CONTROLLER_CORE, SAMPLER_CORE, NUMA_NODE):
        raise TimingError("runtime CPU/NUMA selection differs from campaign")
    current = topology_record(int(design["core"]), int(design["numa_node"]))
    frozen = design["topology"]
    for key in (
            "core", "numa_node", "thread_siblings_list", "sibling",
            "llc_level", "llc_bytes", "llc_shared_cpu_list",
            "llc_shared_cpus", "governor", "energy_performance_preference"):
        if current.get(key) != frozen.get(key):
            raise TimingError("timing topology/power policy changed: %s" % key)
    controller = topology_record(
        int(design["controller_core"]), int(design["numa_node"]))
    if controller != design.get("controller_topology"):
        raise TimingError("controller topology/power policy changed")
    sampler = topology_record(
        int(design["sampler_core"]), int(design["numa_node"]))
    if sampler != design.get("sampler_topology"):
        raise TimingError("sampler topology/power policy changed")
    if (current.get("sibling") != TIMING_SIBLING or
            int(design["controller_core"]) in current["llc_shared_cpus"] or
            int(design["sampler_core"]) in current["llc_shared_cpus"]):
        raise TimingError("frozen CPU isolation topology changed")
    return current


def _execution_receipt(
    task: Mapping[str, object], slot: int, orientation: str,
    command: Sequence[str],
    attempt: int, started_utc: str, start_ns: int, end_ns: int,
    parsed: ParsedOutput, stderr: bytes, prior: Sequence[Mapping[str, object]],
    binary_sha256: str,
) -> Dict[str, object]:
    return sealed_record(
        "wirehair.wh2.degree_balanced_timing.execution_receipt.v1", {
            "job": task["job"], "task_id": task["task_id"],
            "task_sha256": sha256_bytes(canonical_json(task)),
            "outer_slot": slot, "outer_marker": OUTER_ORDER[slot],
            "orientation": orientation, "binary_sha256": binary_sha256,
            "argv": list(command), "attempt": attempt,
            "started_utc": started_utc, "start_monotonic_ns": start_ns,
            "end_monotonic_ns": end_ns, "duration_ns": end_ns - start_ns,
            "stderr_sha256": sha256_bytes(stderr),
            "prior_contamination_receipts": list(prior),
            **_receipt_summary(parsed),
        })


def _save_contamination(
    root: Path, name: str, attempt: int, raw: bytes, stderr: bytes,
    parsed: ParsedOutput, command: Sequence[str], start_ns: int, end_ns: int,
) -> str:
    prefix = "%s.attempt%d" % (name, attempt)
    raw_name = prefix + ".stdout"
    stderr_name = prefix + ".stderr"
    write_new(root / "contamination" / raw_name, raw)
    write_new(root / "contamination" / stderr_name, stderr)
    receipt = sealed_record(
        "wirehair.wh2.degree_balanced_timing.contamination_receipt.v1", {
            "name": name, "attempt": attempt, "argv": list(command),
            "start_monotonic_ns": start_ns, "end_monotonic_ns": end_ns,
            "stdout_name": raw_name, "stdout_sha256": sha256_bytes(raw),
            "stderr_name": stderr_name, "stderr_sha256": sha256_bytes(stderr),
            "contaminations": list(parsed.contaminations),
        })
    receipt_name = prefix + ".json"
    write_new(root / "contamination" / receipt_name, canonical_json(receipt))
    return sha256_file(root / "contamination" / receipt_name)


def run_campaign(args: argparse.Namespace) -> None:
    root = Path(args.result_dir).resolve()
    design = _load_design(root)
    controller_core = int(design["controller_core"])
    os.sched_setaffinity(0, {controller_core})
    if os.sched_getaffinity(0) != {controller_core}:
        raise TimingError("campaign controller could not pin outside the timing LLC")
    tasks = _load_tasks(root, design)
    _verify_immutable(root, design)
    _validate_prepare_receipt(root, design)
    _fresh_output_preflight(root, design)
    if _filler_pids():
        raise TimingError("load filler workers are still running")
    topology = _validate_topology_again(design)
    core = int(design["core"])
    sibling = int(topology["sibling"])
    before_core = cpu_ticks(core)
    before_sibling = cpu_ticks(sibling)
    time.sleep(0.25)
    quiet_core = busy_ticks(cpu_ticks(core)) - busy_ticks(before_core)
    quiet_sibling = busy_ticks(cpu_ticks(sibling)) - busy_ticks(before_sibling)
    if quiet_core > 1 or quiet_sibling > 1:
        raise TimingError("timing core/sibling are not quiet before launch")
    thermal_reader, thermal_csv, thermal_pid = _discover_external_thermal_reader(
        design, args.thermal_csv, args.thermal_pid_file)
    open_flags = os.O_RDONLY | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        open_flags |= os.O_NOFOLLOW
    thermal_stream = os.fdopen(
        os.open(str(thermal_csv), open_flags), "rb", buffering=0)
    thermal_fd = thermal_stream.fileno()
    thermal_stat = os.fstat(thermal_fd)
    if not stat.S_ISREG(thermal_stat.st_mode):
        thermal_stream.close()
        raise TimingError("thermal CSV source is not a regular file")
    thermal_identity = (thermal_stat.st_dev, thermal_stat.st_ino)
    # Read through the held descriptor so pathname replacement cannot switch
    # the evidence source between the identity gate and interval collection.
    thermal_fd_path = Path("/proc/self/fd/%d" % thermal_fd)
    latest = _latest_thermal_row(thermal_fd_path)
    campaign_start_s = time.monotonic()
    if campaign_start_s - float(latest["monotonic_s"]) > MAX_THERMAL_MARGIN_S:
        raise TimingError("thermal sampler has no fresh baseline")
    started_utc = utc_now()
    core_ticks_before = cpu_ticks(core)
    sibling_ticks_before = cpu_ticks(sibling)
    execution_hashes: List[Dict[str, object]] = []
    task_hashes: List[Dict[str, object]] = []
    retry_count = 0
    immutable = design["immutable_files"]
    cross_cache: Dict[
        Tuple[object, ...], Dict[str, Dict[str, object]]
    ] = {}
    for task in tasks:
        parsed_outputs: List[Tuple[str, ParsedOutput]] = []
        execution_records: List[Dict[str, object]] = []
        for slot, marker in enumerate(OUTER_ORDER):
            orientation = OUTER_ORIENTATIONS[marker]
            name = execution_name(task, slot, orientation)
            command = command_for(design, task, orientation)
            prior: List[Dict[str, object]] = []
            for attempt in range(MAX_ENVIRONMENTAL_ATTEMPTS):
                attempt_utc = utc_now()
                start_ns = time.monotonic_ns()
                try:
                    result = subprocess.run(
                        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                        timeout=args.timeout_seconds, check=False,
                        start_new_session=True)
                except subprocess.TimeoutExpired as exc:
                    raise TimingError("timing subprocess exceeded its timeout") from exc
                end_ns = time.monotonic_ns()
                if result.returncode != 0 or result.stderr:
                    raise TimingError(
                        "substantive timing command failure exit=%d stderr=%r" %
                        (result.returncode, result.stderr[:1000]))
                parsed = parse_grouped_output(
                    result.stdout, orientation, task,
                    int(design["evict_bytes"]), core)
                if parsed.contaminations:
                    digest = _save_contamination(
                        root, name, attempt, result.stdout, result.stderr,
                        parsed, command, start_ns, end_ns)
                    prior.append({"attempt": attempt, "receipt_sha256": digest})
                    retry_count += 1
                    if attempt + 1 == MAX_ENVIRONMENTAL_ATTEMPTS:
                        raise ContaminationError(
                            "environmental retry limit exhausted: %s" %
                            ",".join(parsed.contaminations))
                    continue
                raw_path = root / "raw" / name
                stderr_path = root / "stderr" / (name + ".stderr")
                exit_path = root / "exit" / (name + ".exit")
                receipt_path = root / "receipts" / (name + ".json")
                write_new(raw_path, result.stdout)
                write_new(stderr_path, result.stderr)
                write_new(exit_path, b"0\n")
                receipt = _execution_receipt(
                    task, slot, orientation, command, attempt, attempt_utc, start_ns,
                    end_ns, parsed, result.stderr, prior,
                    str(immutable["frozen/" + BINARY_NAMES["measurement"]]))
                write_new(receipt_path, canonical_json(receipt))
                receipt_hash = sha256_file(receipt_path)
                record = {"name": name, "receipt_sha256": receipt_hash}
                execution_records.append(record)
                execution_hashes.append(record)
                parsed_outputs.append((orientation, parsed))
                break
        semantic = {item.semantic_sha256 for _orientation, item in parsed_outputs}
        traces = {item.preamble["trace_sha256"]
                  for _orientation, item in parsed_outputs}
        work = {sha256_bytes(canonical_json({
            arm: list(item.work_signatures[arm])
            for arm in ("base", "balanced")}))
            for _orientation, item in parsed_outputs}
        outcome_classes = {_semantic_outcome_class(item.outcomes)
                           for _orientation, item in parsed_outputs}
        if (len(semantic) != 1 or len(traces) != 1 or len(work) != 1 or
                len(outcome_classes) != 1):
            raise TimingError("forward/reverse executions changed semantics")
        _register_cross_cache_identity(cross_cache, task, parsed_outputs[0][1])
        base_elapsed = sum(item.base_timed_elapsed_ns
                           for _orientation, item in parsed_outputs)
        balanced_elapsed = sum(item.balanced_timed_elapsed_ns
                               for _orientation, item in parsed_outputs)
        if base_elapsed <= 0 or balanced_elapsed <= 0:
            raise TimingError("task timing sum is empty")
        task_receipt = sealed_record(
            "wirehair.wh2.degree_balanced_timing.task_receipt.v1", {
                "job": task["job"], "task_id": task["task_id"],
                "task_sha256": sha256_bytes(canonical_json(task)),
                "outer_order": OUTER_ORDER,
                "execution_receipts": execution_records,
                "trace_sha256": next(iter(traces)),
                "semantic_sha256": next(iter(semantic)),
                "work_signatures_sha256": next(iter(work)),
                "semantic_outcome_class": next(iter(outcome_classes)),
                "common_success": next(iter(outcome_classes)) == "common-success",
                "base_timed_elapsed_ns": base_elapsed,
                "balanced_timed_elapsed_ns": balanced_elapsed,
                "ratio": {"numerator": balanced_elapsed,
                          "denominator": base_elapsed},
                "forward_process_count": 4, "reverse_process_count": 4,
                "timed_rows_per_architecture": 96,
            })
        task_path = root / "task_receipts" / (str(task["task_id"]) + ".json")
        write_new(task_path, canonical_json(task_receipt))
        task_hashes.append({
            "task_id": task["task_id"], "receipt_sha256": sha256_file(task_path)})
        if _filler_pids():
            raise TimingError("load filler worker appeared during timing")
        if (int(task["job"]) + 1) % 6 == 0 or int(task["job"]) + 1 == len(tasks):
            print("progress=%d/%d retries=%d" %
                  (int(task["job"]) + 1, len(tasks), retry_count), flush=True)
    _validate_cross_cache_ledger(cross_cache)
    campaign_end_s = time.monotonic()
    thermal_reader_end = _thermal_reader(design, thermal_pid)
    if thermal_reader_end != thermal_reader:
        raise TimingError("thermal reader identity/affinity changed during timing")
    thermal_stat_end = os.fstat(thermal_fd)
    if (thermal_stat_end.st_dev, thermal_stat_end.st_ino) != thermal_identity:
        raise TimingError("thermal CSV inode changed during timing")
    thermal_raw, thermal_summary = collect_thermal_interval(
        thermal_fd_path, campaign_start_s, campaign_end_s)
    thermal_stat_collected = os.fstat(thermal_fd)
    if (thermal_stat_collected.st_dev, thermal_stat_collected.st_ino) != \
            thermal_identity:
        raise TimingError("thermal CSV inode changed during interval collection")
    thermal_stream.close()
    write_new(root / "thermal_interval.csv", thermal_raw)
    core_ticks_after = cpu_ticks(core)
    sibling_ticks_after = cpu_ticks(sibling)
    sibling_busy = busy_ticks(sibling_ticks_after) - busy_ticks(sibling_ticks_before)
    if sibling_busy > 1:
        raise TimingError("SMT sibling accumulated busy ticks during timing")
    launch = sealed_record(
        "wirehair.wh2.degree_balanced_timing.launch_receipt.v1", {
            "started_utc": started_utc, "ended_utc": utc_now(),
            "start_monotonic_s": campaign_start_s,
            "end_monotonic_s": campaign_end_s,
            "duration_s": campaign_end_s - campaign_start_s,
            "design_sha256": sha256_file(root / "design.json"),
            "prepare_receipt_sha256": sha256_file(root / "prepare_receipt.json"),
            "tasks_manifest_sha256": sha256_file(root / "tasks_manifest.jsonl"),
            "task_count": len(tasks), "execution_count": len(execution_hashes),
            "retry_count": retry_count,
            "execution_receipts": execution_hashes,
            "task_receipts": task_hashes,
            "thermal_reader": thermal_reader,
            "thermal_source_device": thermal_identity[0],
            "thermal_source_inode": thermal_identity[1],
            "thermal_interval_sha256": sha256_bytes(thermal_raw),
            "thermal_summary": thermal_summary,
            "topology": topology, "load_workers_stopped": True,
            "controller_core": controller_core,
            "controller_affinity": sorted(os.sched_getaffinity(0)),
            "core_ticks_before": list(core_ticks_before),
            "core_ticks_after": list(core_ticks_after),
            "sibling_ticks_before": list(sibling_ticks_before),
            "sibling_ticks_after": list(sibling_ticks_after),
            "sibling_busy_ticks": sibling_busy,
            "preflight_quiet_core_ticks": quiet_core,
            "preflight_quiet_sibling_ticks": quiet_sibling,
        })
    write_new(root / "launch_receipt.json", canonical_json(launch))
    print(json.dumps({
        "task_count": len(tasks), "execution_count": len(execution_hashes),
        "retry_count": retry_count, "duration_s": launch["duration_s"],
        "thermal": thermal_summary,
        "launch_receipt_sha256": sha256_file(root / "launch_receipt.json"),
    }, sort_keys=True), flush=True)


def _validate_execution_receipt(
    root: Path, design: Mapping[str, object], task: Mapping[str, object],
    slot: int, not_before_ns: int, not_after_ns: int,
) -> Tuple[Dict[str, object], ParsedOutput]:
    marker = OUTER_ORDER[slot]
    orientation = OUTER_ORIENTATIONS[marker]
    name = execution_name(task, slot, orientation)
    raw = (root / "raw" / name).read_bytes()
    stderr = (root / "stderr" / (name + ".stderr")).read_bytes()
    if stderr or (root / "exit" / (name + ".exit")).read_bytes() != b"0\n":
        raise TimingError("execution stderr/exit artifact changed")
    parsed = parse_grouped_output(
        raw, orientation, task, int(design["evict_bytes"]), int(design["core"]))
    if parsed.contaminations:
        raise TimingError("accepted execution now parses as contaminated")
    path = root / "receipts" / (name + ".json")
    receipt = load_canonical(path, "execution receipt")
    verify_sealed_record(
        receipt, "wirehair.wh2.degree_balanced_timing.execution_receipt.v1",
        "execution receipt")
    if set(receipt) != EXECUTION_RECEIPT_FIELDS:
        raise TimingError("execution receipt fields changed")
    if (receipt.get("job") != task["job"] or
            receipt.get("task_id") != task["task_id"] or
            receipt.get("task_sha256") != sha256_bytes(canonical_json(task)) or
            receipt.get("outer_slot") != slot or
            receipt.get("outer_marker") != marker or
            receipt.get("orientation") != orientation or
            receipt.get("argv") != command_for(design, task, orientation) or
            receipt.get("stderr_sha256") != sha256_bytes(stderr)):
        raise TimingError("execution receipt binding mismatch")
    expected_summary = _receipt_summary(parsed)
    if any(receipt.get(key) != value for key, value in expected_summary.items()):
        raise TimingError("execution receipt parsed summary mismatch")
    if receipt.get("binary_sha256") != design["immutable_files"].get(
            "frozen/" + BINARY_NAMES["measurement"]):
        raise TimingError("execution receipt binary binding mismatch")
    start = receipt.get("start_monotonic_ns")
    end = receipt.get("end_monotonic_ns")
    attempt = receipt.get("attempt")
    prior = receipt.get("prior_contamination_receipts")
    if (not isinstance(start, int) or isinstance(start, bool) or
            not isinstance(end, int) or isinstance(end, bool) or end <= start or
            start < not_before_ns or end > not_after_ns or
            receipt.get("duration_ns") != end - start or
            not isinstance(attempt, int) or isinstance(attempt, bool) or
            not 0 <= attempt < MAX_ENVIRONMENTAL_ATTEMPTS or
            not isinstance(prior, list) or len(prior) != attempt):
        raise TimingError("execution timing/retry ledger mismatch")
    if (not isinstance(receipt.get("started_utc"), str) or
            UTC_RE.fullmatch(receipt["started_utc"]) is None):
        raise TimingError("execution UTC receipt is malformed")
    previous_contamination_end: Optional[int] = None
    for index, record in enumerate(prior):
        if (not isinstance(record, dict) or
                set(record) != {"attempt", "receipt_sha256"} or
                record.get("attempt") != index):
            raise TimingError("contamination retry order mismatch")
        require_sha256(record.get("receipt_sha256"), "contamination receipt hash")
        prefix = "%s.attempt%d" % (name, index)
        contamination_path = root / "contamination" / (prefix + ".json")
        if sha256_file(contamination_path) != record["receipt_sha256"]:
            raise TimingError("contamination receipt hash mismatch")
        contamination = load_canonical(contamination_path, "contamination receipt")
        verify_sealed_record(
            contamination,
            "wirehair.wh2.degree_balanced_timing.contamination_receipt.v1",
            "contamination receipt")
        expected_stdout_name = prefix + ".stdout"
        expected_stderr_name = prefix + ".stderr"
        if (set(contamination) != CONTAMINATION_RECEIPT_FIELDS or
                contamination.get("name") != name or
                contamination.get("attempt") != index or
                contamination.get("argv") != command_for(
                    design, task, orientation) or
                contamination.get("stdout_name") != expected_stdout_name or
                contamination.get("stderr_name") != expected_stderr_name):
            raise TimingError("contamination receipt coordinate mismatch")
        for stream in ("stdout", "stderr"):
            artifact = root / "contamination" / str(contamination[stream + "_name"])
            if sha256_file(artifact) != contamination.get(stream + "_sha256"):
                raise TimingError("contamination artifact hash mismatch")
        contaminated = parse_grouped_output(
            (root / "contamination" / str(contamination["stdout_name"])).read_bytes(),
            orientation, task, int(design["evict_bytes"]), int(design["core"]))
        contamination_stderr = (
            root / "contamination" / expected_stderr_name).read_bytes()
        contamination_start = contamination.get("start_monotonic_ns")
        contamination_end = contamination.get("end_monotonic_ns")
        if (contamination_stderr or
                not isinstance(contamination_start, int) or
                isinstance(contamination_start, bool) or
                not isinstance(contamination_end, int) or
                isinstance(contamination_end, bool) or
                contamination_end <= contamination_start or
                contamination_start < not_before_ns or
                (previous_contamination_end is not None and
                 contamination_start < previous_contamination_end) or
                contamination_end > start or
                not contaminated.contaminations or
                contaminated.semantic_sha256 != parsed.semantic_sha256 or
                contaminated.work_signatures != parsed.work_signatures or
                contaminated.preamble["trace_sha256"] !=
                parsed.preamble["trace_sha256"] or
                list(contaminated.contaminations) !=
                contamination.get("contaminations")):
            raise TimingError("contamination classification changed on replay")
        previous_contamination_end = contamination_end
    return receipt, parsed


def _fraction_record(value: Fraction) -> Dict[str, object]:
    return {
        "numerator": value.numerator, "denominator": value.denominator,
        "decimal": format(float(value), ".12f"),
    }


def _upper_median(values: Sequence[Fraction]) -> Fraction:
    if not values:
        raise TimingError("cannot summarize an empty timing group")
    return sorted(values)[len(values) // 2]


def _bootstrap(
    records: Sequence[Mapping[str, object]], domain: str,
    repetitions: int = BOOTSTRAP_REPS,
) -> Dict[str, object]:
    if not records or repetitions < 100:
        raise TimingError("bootstrap domain is empty or undersized")
    seed = int.from_bytes(hashlib.sha256(
        b"wirehair.wh2.degree-balanced.bootstrap.v1\0" + domain.encode("ascii")
    ).digest()[:8], "big")
    generator = random.Random(seed)
    clustered: Dict[Tuple[object, ...], List[Mapping[str, object]]] = {}
    for record in records:
        key = tuple(record["cluster"])
        clustered.setdefault(key, []).append(record)
    cluster_records = [(
        sum(int(record["base"]) for record in values),
        sum(int(record["candidate"]) for record in values),
    ) for _key, values in sorted(clustered.items())]
    values = []
    count = len(cluster_records)
    for _index in range(repetitions):
        base = 0
        balanced = 0
        for _draw in range(count):
            selected = cluster_records[generator.randrange(count)]
            base += selected[0]
            balanced += selected[1]
        values.append(balanced / base)
    values.sort()
    low = values[int(math.floor(0.025 * (repetitions - 1)))]
    high = values[int(math.ceil(0.975 * (repetitions - 1)))]
    upper_one_sided = values[int(math.ceil(0.95 * (repetitions - 1)))]
    return {
        "schema": "paired-task-ratio-of-sums-bootstrap-v1",
        "seed": seed, "repetitions": repetitions,
        "cluster_count": count,
        "cluster_coordinates": "K,bb,schedule,seed (cold/warm together when present)",
        "lower_95": format(low, ".12f"),
        "upper_95": format(high, ".12f"),
        "upper_one_sided_95": format(upper_one_sided, ".12f"),
    }


def _summarize(
    records: Sequence[Mapping[str, object]], domain: str,
) -> Dict[str, object]:
    if not records:
        raise TimingError("cannot summarize an empty record group")
    base = sum(int(record["base"]) for record in records)
    balanced = sum(int(record["candidate"]) for record in records)
    ratios = [Fraction(int(record["candidate"]), int(record["base"]))
              for record in records]
    faster = sum(1 for ratio in ratios if ratio < 1)
    slower = sum(1 for ratio in ratios if ratio > 1)
    ties = len(ratios) - faster - slower
    return {
        "task_count": len(records), "base_elapsed_ns": base,
        "balanced_elapsed_ns": balanced,
        "ratio_of_sums": _fraction_record(Fraction(balanced, base)),
        "upper_median_task_ratio": _fraction_record(_upper_median(ratios)),
        "balanced_faster_tasks": faster, "balanced_slower_tasks": slower,
        "tied_tasks": ties,
        "bootstrap": _bootstrap(records, domain),
    }


def _work_summary(records: Sequence[Mapping[str, object]]) -> Dict[str, object]:
    fields = (
        "inactivated", "binary_def", "heavy_gain", "block_xors",
        "block_muladds", "joint_source_xors", "joint_marginal_xors",
        "joint_marginal_copies", "joint_active_deltas",
        "joint_scratch_bytes", "dual_source_columns",
    )
    result: Dict[str, object] = {}
    for field in fields:
        base = sum(int(record["work"]["base"][field]) for record in records)
        balanced = sum(
            int(record["work"]["balanced"][field]) for record in records)
        result[field] = {
            "base_sum": base, "balanced_sum": balanced,
            "ratio": _fraction_record(Fraction(balanced, base))
                if base else ("equal-zero" if balanced == 0 else "base-zero"),
        }
    return result


def replay_campaign(root: Path) -> Tuple[Dict[str, object], set[str]]:
    design = _load_design(root)
    _verify_directory_inventory(root)
    tasks = _load_tasks(root, design)
    _verify_immutable(root, design)
    _validate_prepare_receipt(root, design)
    launch = load_canonical(root / "launch_receipt.json", "launch receipt")
    verify_sealed_record(
        launch, "wirehair.wh2.degree_balanced_timing.launch_receipt.v1",
        "launch receipt")
    if set(launch) != LAUNCH_RECEIPT_FIELDS:
        raise TimingError("launch receipt fields changed")
    if (launch.get("design_sha256") != sha256_file(root / "design.json") or
            launch.get("prepare_receipt_sha256") !=
            sha256_file(root / "prepare_receipt.json") or
            launch.get("tasks_manifest_sha256") !=
            sha256_file(root / "tasks_manifest.jsonl") or
            launch.get("task_count") != len(tasks) or
            launch.get("execution_count") != len(tasks) * len(OUTER_ORDER)):
        raise TimingError("launch receipt campaign binding mismatch")
    if launch.get("load_workers_stopped") is not True or \
            launch.get("topology") != design.get("topology") or \
            launch.get("controller_core") != design.get("controller_core") or \
            launch.get("controller_affinity") != [design.get("controller_core")]:
        raise TimingError("launch isolation receipt mismatch")
    if (not isinstance(launch.get("started_utc"), str) or
            UTC_RE.fullmatch(launch["started_utc"]) is None or
            not isinstance(launch.get("ended_utc"), str) or
            UTC_RE.fullmatch(launch["ended_utc"]) is None):
        raise TimingError("launch UTC receipt is malformed")
    thermal_reader = launch.get("thermal_reader")
    thermal_reader_fields = {
        "ownership", "pid", "process_start_ticks", "cpus_allowed_list",
        "cmdline_sha256", "script_path", "script_sha256", "csv_path",
        "pid_file", "unique_i2c_reader_count", "external_sampler_signalled",
    }
    if (not isinstance(thermal_reader, dict) or
            set(thermal_reader) != thermal_reader_fields or
            thermal_reader.get("ownership") != "external-reused" or
            thermal_reader.get("external_sampler_signalled") is not False or
            thermal_reader.get("unique_i2c_reader_count") != 1 or
            not isinstance(thermal_reader.get("pid"), int) or
            isinstance(thermal_reader.get("pid"), bool) or
            thermal_reader["pid"] <= 0 or
            not isinstance(thermal_reader.get("process_start_ticks"), int) or
            isinstance(thermal_reader.get("process_start_ticks"), bool) or
            thermal_reader["process_start_ticks"] < 0 or
            not isinstance(thermal_reader.get("cpus_allowed_list"), str) or
            len(parse_cpu_list(thermal_reader["cpus_allowed_list"])) != 1 or
            any(not isinstance(thermal_reader.get(field), str) or
                not thermal_reader[field] for field in (
                    "script_path", "csv_path", "pid_file"))):
        raise TimingError("launch thermal-reader receipt is malformed")
    require_sha256(thermal_reader.get("cmdline_sha256"), "thermal cmdline hash")
    require_sha256(thermal_reader.get("script_sha256"), "thermal script hash")
    for field in ("thermal_source_device", "thermal_source_inode"):
        value = launch.get(field)
        if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
            raise TimingError("launch thermal source identity is malformed")
    thermal_raw = (root / "thermal_interval.csv").read_bytes()
    if sha256_bytes(thermal_raw) != launch.get("thermal_interval_sha256"):
        raise TimingError("thermal interval hash mismatch")
    start_s = launch.get("start_monotonic_s")
    end_s = launch.get("end_monotonic_s")
    if (not isinstance(start_s, (int, float)) or isinstance(start_s, bool) or
            not isinstance(end_s, (int, float)) or isinstance(end_s, bool) or
            end_s <= start_s):
        raise TimingError("launch monotonic interval is malformed")
    if launch.get("duration_s") != end_s - start_s:
        raise TimingError("launch duration receipt mismatch")
    thermal_summary = validate_sealed_thermal_interval(
        thermal_raw, float(start_s), float(end_s))
    if launch.get("thermal_summary") != thermal_summary:
        raise TimingError("thermal summary does not replay")
    for prefix in ("core", "sibling"):
        before = launch.get(prefix + "_ticks_before")
        after = launch.get(prefix + "_ticks_after")
        if (not isinstance(before, list) or not isinstance(after, list) or
                len(before) < 5 or len(after) != len(before) or
                any(not isinstance(value, int) or isinstance(value, bool) or value < 0
                    for value in before + after) or
                any(right < left for left, right in zip(before, after))):
            raise TimingError("launch CPU tick receipt is malformed")
    sibling_busy = busy_ticks(launch["sibling_ticks_after"]) - \
        busy_ticks(launch["sibling_ticks_before"])
    if sibling_busy > 1 or launch.get("sibling_busy_ticks") != sibling_busy:
        raise TimingError("launch sibling-idle receipt mismatch")
    if (not isinstance(launch.get("preflight_quiet_core_ticks"), int) or
            not 0 <= launch["preflight_quiet_core_ticks"] <= 1 or
            not isinstance(launch.get("preflight_quiet_sibling_ticks"), int) or
            not 0 <= launch["preflight_quiet_sibling_ticks"] <= 1):
        raise TimingError("launch preflight quiet receipt mismatch")
    execution_ledger = []
    task_ledger = []
    result_records: List[Dict[str, object]] = []
    cross_cache: Dict[
        Tuple[object, ...], Dict[str, Dict[str, object]]
    ] = {}
    expected_files = {
        "design.json", "prepare_receipt.json", "tasks_manifest.jsonl",
        "launch_receipt.json", "thermal_interval.csv",
    }
    for relative in design["immutable_files"]:
        expected_files.add(str(relative))
    contamination_expected = set()
    replay_retry_count = 0
    previous_execution_end = int(math.floor(float(start_s) * 1_000_000_000)) - 1000
    campaign_end_ns = int(math.ceil(float(end_s) * 1_000_000_000)) + 1000
    for task in tasks:
        parsed_outputs = []
        execution_records = []
        for slot in range(len(OUTER_ORDER)):
            receipt, parsed = _validate_execution_receipt(
                root, design, task, slot, previous_execution_end,
                campaign_end_ns)
            previous_execution_end = int(receipt["end_monotonic_ns"])
            replay_retry_count += int(receipt["attempt"])
            orientation = str(receipt["orientation"])
            name = execution_name(task, slot, orientation)
            receipt_path = root / "receipts" / (name + ".json")
            record = {"name": name, "receipt_sha256": sha256_file(receipt_path)}
            execution_ledger.append(record)
            execution_records.append(record)
            parsed_outputs.append((orientation, parsed))
            expected_files.update((
                "raw/" + name, "stderr/" + name + ".stderr",
                "exit/" + name + ".exit", "receipts/" + name + ".json"))
            for prior in receipt["prior_contamination_receipts"]:
                attempt = int(prior["attempt"])
                prefix = "%s.attempt%d" % (name, attempt)
                contamination_expected.update((
                    "contamination/" + prefix + ".stdout",
                    "contamination/" + prefix + ".stderr",
                    "contamination/" + prefix + ".json"))
        semantic = {parsed.semantic_sha256
                    for _orientation, parsed in parsed_outputs}
        traces = {parsed.preamble["trace_sha256"]
                  for _orientation, parsed in parsed_outputs}
        work = {sha256_bytes(canonical_json({
            arm: list(parsed.work_signatures[arm])
            for arm in ("base", "balanced")}))
            for _orientation, parsed in parsed_outputs}
        outcomes = {_semantic_outcome_class(parsed.outcomes)
                    for _orientation, parsed in parsed_outputs}
        if (len(semantic) != 1 or len(traces) != 1 or len(work) != 1 or
                len(outcomes) != 1):
            raise TimingError("replay found forward/reverse semantic drift")
        _register_cross_cache_identity(cross_cache, task, parsed_outputs[0][1])
        base = sum(parsed.base_timed_elapsed_ns
                   for _orientation, parsed in parsed_outputs)
        candidate = sum(parsed.balanced_timed_elapsed_ns
                        for _orientation, parsed in parsed_outputs)
        task_path = root / "task_receipts" / (str(task["task_id"]) + ".json")
        task_receipt = load_canonical(task_path, "task receipt")
        verify_sealed_record(
            task_receipt, "wirehair.wh2.degree_balanced_timing.task_receipt.v1",
            "task receipt")
        if set(task_receipt) != TASK_RECEIPT_FIELDS:
            raise TimingError("task receipt fields changed")
        expected_task = {
            "job": task["job"], "task_id": task["task_id"],
            "task_sha256": sha256_bytes(canonical_json(task)),
            "outer_order": OUTER_ORDER, "execution_receipts": execution_records,
            "trace_sha256": next(iter(traces)),
            "semantic_sha256": next(iter(semantic)),
            "work_signatures_sha256": next(iter(work)),
            "semantic_outcome_class": next(iter(outcomes)),
            "common_success": next(iter(outcomes)) == "common-success",
            "base_timed_elapsed_ns": base,
            "balanced_timed_elapsed_ns": candidate,
            "ratio": {"numerator": candidate, "denominator": base},
            "forward_process_count": 4, "reverse_process_count": 4,
            "timed_rows_per_architecture": 96,
        }
        if any(task_receipt.get(key) != value for key, value in expected_task.items()):
            raise TimingError("task receipt replay mismatch")
        task_record = {
            "task_id": task["task_id"], "receipt_sha256": sha256_file(task_path)}
        task_ledger.append(task_record)
        expected_files.add("task_receipts/" + task_path.name)
        representative = parsed_outputs[0][1]
        work_by_arm = {
            arm: dict(zip(WORK_FIELDS, representative.work_signatures[arm]))
            for arm in ("base", "balanced")}
        result_records.append({
            "K": task["K"], "bb": task["bb"],
            "schedule": task["schedule"], "seed_index": task["seed_index"],
            "seed": task["seed"], "cache_state": task["cache_state"],
            "cluster": [task["K"], task["bb"], task["schedule"], task["seed"]],
            "outcome_class": next(iter(outcomes)),
            "work": work_by_arm,
            "base": base, "candidate": candidate,
        })
    _validate_cross_cache_ledger(cross_cache)
    expected_files.update(contamination_expected)
    if launch.get("execution_receipts") != execution_ledger or \
            launch.get("task_receipts") != task_ledger:
        raise TimingError("launch receipt ledger does not replay exactly")
    if launch.get("retry_count") != replay_retry_count:
        raise TimingError("launch retry count does not replay")
    actual_contamination = {
        "contamination/" + path.name for path in (root / "contamination").iterdir()
    }
    if actual_contamination != contamination_expected:
        raise TimingError("contamination artifact inventory changed")
    optional_final = {
        "validated_summary.json", "data_manifest.json", "data_manifest.sha256"}
    actual_files = {
        str(path.relative_to(root)) for path in root.rglob("*") if path.is_file()
    }
    if not expected_files <= actual_files or actual_files - expected_files - optional_final:
        raise TimingError("campaign artifact inventory has missing or extra files")
    schedule_names = tuple(dict.fromkeys(item[0] for item in SCHEDULE_SEEDS))
    overall = _summarize(result_records, "overall-unconditional")
    by_cache = {
        cache: _summarize(
            [record for record in result_records
             if record["cache_state"] == cache], "cache:" + cache)
        for cache in CACHE_STATES}
    by_K = {
        str(K): _summarize(
            [record for record in result_records if record["K"] == K],
            "K:%d" % K) for K in KS}
    by_bb = {
        str(width): _summarize(
            [record for record in result_records if record["bb"] == width],
            "bb:%d" % width) for width in WIDTHS}
    by_schedule = {
        schedule: _summarize(
            [record for record in result_records
             if record["schedule"] == schedule], "schedule:" + schedule)
        for schedule in schedule_names}
    common_success_records = [
        record for record in result_records
        if record["outcome_class"] == "common-success"]
    all_common_success = len(common_success_records) == len(result_records)
    overall_speed_bound = float(overall["bootstrap"]["upper_95"]) < 1.0
    cold_warm_speed_bounds = all(
        float(by_cache[cache]["bootstrap"]["upper_95"]) < 1.0
        for cache in CACHE_STATES)
    noninferiority_strata = {
        **{"K:" + key: value for key, value in by_K.items()},
        **{"bb:" + key: value for key, value in by_bb.items()},
        **{"schedule:" + key: value for key, value in by_schedule.items()},
    }
    noninferiority_pass = all(
        float(value["bootstrap"]["upper_one_sided_95"]) < 1.01
        for value in noninferiority_strata.values())
    timing_promotional = (
        all_common_success and overall_speed_bound and
        cold_warm_speed_bounds and noninferiority_pass)
    outcome_counts = {
        outcome: sum(record["outcome_class"] == outcome
                     for record in result_records)
        for outcome in ("common-success", "base-only", "balanced-only",
                        "common-failure")}
    summary_payload: Dict[str, object] = {
        "base_commit": BASE_COMMIT,
        "architecture_commit": ARCHITECTURE_COMMIT,
        "measurement_commit": MEASUREMENT_COMMIT,
        "architecture": ARCHITECTURE, "task_count": len(tasks),
        "execution_count": len(execution_ledger),
        "timed_rows_per_architecture": len(tasks) * 96,
        "trace_payload_dimensions_seed_receipts_exact": True,
        "cold_warm_graph_trace_work_exact": True,
        "work_may_differ_between_architectures": True,
        "outcome_counts": outcome_counts,
        "all_fixed_cells_common_success": all_common_success,
        "unconditional_attempt_latency": {
            "overall": overall, "by_cache_state": by_cache,
            "by_K": by_K, "by_bb": by_bb, "by_schedule": by_schedule,
        },
        "common_success_latency_diagnostic": _summarize(
            common_success_records, "common-success-diagnostic")
            if common_success_records else None,
        "deterministic_work": _work_summary(result_records),
        "promotion_gates": {
            "overall_paired_95_upper_below_one": overall_speed_bound,
            "cold_and_warm_95_uppers_below_one": cold_warm_speed_bounds,
            "every_K_width_schedule_one_sided_95_upper_below_1_01":
                noninferiority_pass,
            "no_post_outcome_filtering": True,
            "all_fixed_cells_common_success": all_common_success,
            "timing_promotional": timing_promotional,
            "architecture_promotion_requires_favorable_sealed_all_K_recovery": True,
            "dual_improvement_claim_from_timing_alone": False,
        },
        "timing_promotional": timing_promotional,
        "thermal_interval_sha256": sha256_bytes(thermal_raw),
        "thermal_summary": thermal_summary,
    }
    return summary_payload, expected_files


def verify_prepared_campaign(args: argparse.Namespace) -> None:
    """Replay every frozen preparation receipt without starting timing."""
    root = Path(args.result_dir).resolve()
    design = _load_design(root)
    tasks = _load_tasks(root, design)
    _verify_immutable(root, design)
    prepare = _validate_prepare_receipt(root, design)
    _fresh_output_preflight(root, design)
    print(json.dumps({
        "verified_prepared": True,
        "task_count": len(tasks),
        "design_sha256": sha256_file(root / "design.json"),
        "prepare_receipt_sha256": sha256_file(root / "prepare_receipt.json"),
        "tasks_manifest_sha256": sha256_file(root / "tasks_manifest.jsonl"),
        "exact_base_binary_sha256": prepare["exact_base_binary_sha256"],
        "exact_architecture_binary_sha256":
            prepare["exact_architecture_binary_sha256"],
        "measurement_binary_sha256": prepare["measurement_binary_sha256"],
    }, sort_keys=True))


def _data_manifest(root: Path, expected: Iterable[str]) -> List[Dict[str, object]]:
    records = []
    for relative in sorted(set(expected)):
        path = root / relative
        if path.is_symlink() or not path.is_file():
            raise TimingError("data manifest path is not a regular file: %s" % relative)
        records.append({
            "path": relative, "size": path.stat().st_size,
            "sha256": sha256_file(path),
        })
    return records


def reduce_campaign(args: argparse.Namespace) -> None:
    root = Path(args.result_dir).resolve()
    for name in ("validated_summary.json", "data_manifest.json",
                 "data_manifest.sha256"):
        if (root / name).exists():
            raise TimingError("reducer refuses to replace %s" % name)
    payload, expected = replay_campaign(root)
    summary = sealed_record(
        "wirehair.wh2.degree_balanced_timing.validated_summary.v1", payload)
    summary_path = root / "validated_summary.json"
    write_new(summary_path, canonical_json(summary))
    expected.add("validated_summary.json")
    manifest = _data_manifest(root, expected)
    manifest_path = root / "data_manifest.json"
    write_new(manifest_path, canonical_json({"files": manifest}))
    write_new(root / "data_manifest.sha256",
              (sha256_file(manifest_path) + "\n").encode("ascii"))
    print(json.dumps({
        "summary_sha256": sha256_file(summary_path),
        "data_manifest_sha256": sha256_file(manifest_path),
        "overall_ratio": summary["unconditional_attempt_latency"]["overall"]
            ["ratio_of_sums"]["decimal"],
        "cold_ratio": summary["unconditional_attempt_latency"]
            ["by_cache_state"]["cold"]["ratio_of_sums"]["decimal"],
        "warm_ratio": summary["unconditional_attempt_latency"]
            ["by_cache_state"]["warm"]["ratio_of_sums"]["decimal"],
        "timing_promotional": summary["timing_promotional"],
    }, sort_keys=True))


def verify_campaign(args: argparse.Namespace) -> None:
    root = Path(args.result_dir).resolve()
    payload, expected = replay_campaign(root)
    summary = load_canonical(root / "validated_summary.json", "validated summary")
    verify_sealed_record(
        summary, "wirehair.wh2.degree_balanced_timing.validated_summary.v1",
        "validated summary")
    expected_summary = sealed_record(
        "wirehair.wh2.degree_balanced_timing.validated_summary.v1", payload)
    if summary != expected_summary:
        raise TimingError("validated summary does not reproduce")
    expected.add("validated_summary.json")
    manifest = load_canonical(root / "data_manifest.json", "data manifest")
    if set(manifest) != {"files"} or manifest["files"] != _data_manifest(root, expected):
        raise TimingError("data manifest does not reproduce")
    sidecar = (root / "data_manifest.sha256").read_text(encoding="ascii")
    if sidecar != sha256_file(root / "data_manifest.json") + "\n":
        raise TimingError("data manifest sidecar mismatch")
    all_files = {
        str(path.relative_to(root)) for path in root.rglob("*") if path.is_file()
    }
    allowed = set(expected) | {"data_manifest.json", "data_manifest.sha256"}
    if all_files != allowed:
        raise TimingError("campaign root contains missing or extra files")
    print(json.dumps({
        "verified": True,
        "summary_sha256": sha256_file(root / "validated_summary.json"),
        "data_manifest_sha256": sha256_file(root / "data_manifest.json"),
        "file_count": len(all_files),
    }, sort_keys=True))


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare = subparsers.add_parser("prepare", help="build and freeze exact binaries")
    prepare.add_argument("--repo", required=True)
    prepare.add_argument("--result-dir", required=True)
    prepare.add_argument("--core", type=int, required=True)
    prepare.add_argument("--controller-core", type=int, required=True)
    prepare.add_argument("--numa-node", type=int, required=True)
    prepare.add_argument("--evict-bytes", type=int, default=256 * 1024 * 1024)
    prepare.add_argument("--build-jobs", type=int, default=max(1, os.cpu_count() or 1))
    prepare.add_argument("--c-compiler")
    prepare.add_argument("--cxx-compiler")
    prepare.set_defaults(function=prepare_campaign)
    run = subparsers.add_parser("run", help="run one fresh isolated campaign")
    run.add_argument("--result-dir", required=True)
    run.add_argument(
        "--thermal-csv",
        help="optional assertion; the sole external I2C reader is auto-discovered")
    run.add_argument(
        "--thermal-pid-file",
        help="optional assertion; the sole external I2C reader is auto-discovered")
    run.add_argument("--timeout-seconds", type=int, default=1800)
    run.set_defaults(function=run_campaign)
    verify_prepared = subparsers.add_parser(
        "verify-prepared", help="replay a frozen preparation without launching")
    verify_prepared.add_argument("--result-dir", required=True)
    verify_prepared.set_defaults(function=verify_prepared_campaign)
    reduce = subparsers.add_parser("reduce", help="strictly replay and summarize")
    reduce.add_argument("--result-dir", required=True)
    reduce.set_defaults(function=reduce_campaign)
    verify = subparsers.add_parser("verify", help="replay a reduced campaign")
    verify.add_argument("--result-dir", required=True)
    verify.set_defaults(function=verify_campaign)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = make_parser()
    args = parser.parse_args(argv)
    try:
        if getattr(args, "timeout_seconds", 1) <= 0:
            raise TimingError("timeout must be positive")
        args.function(args)
    except (OSError, ValueError, TimingError) as exc:
        print("degree-balanced timing failed: %s" % exc, file=sys.stderr, flush=True)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
