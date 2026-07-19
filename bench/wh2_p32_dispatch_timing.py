#!/usr/bin/env python3
"""Transactional sealed timing for WH2 P32/r7 versus prod244.

The campaign compares two raw recovery architectures in one paired benchmark
binary.  Preparation freezes an exact committed source tree, executable,
runner, thermal sampler, task ledger, build provenance, and external recovery
receipts.  ``run`` is bounded and resumable at atomic task-directory
boundaries.  Each run segment owns the sole CPU/DIMM sampler from launch
through a post-benchmark sample and graceful shutdown.  ``verify`` replays all
receipts; ``reduce`` publishes aggregate paired ratios only after verification.

The obsolete a3f5c66 v1 artifact is provenance only.  This runner neither
reads its result directories nor launches or modifies its frozen executable.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import fcntl
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import selectors
import shutil
import signal
import statistics
import subprocess
import sys
import time
from typing import Dict, List, Mapping, Optional, Sequence, Tuple


sys.dont_write_bytecode = True

KS = (3199, 3200, 4095, 4096, 6354, 6355, 8191, 8192,
      8887, 8888, 9999, 10000, 10001)
WIDTHS = (64, 256, 512, 1024, 1280, 4096)
SCHEDULES = ("burst", "adversarial", "repair-only")
SEEDS = (15111065706836454659, 10723151780598845931,
         9599682871048892067)
CACHE_STATES = ("cold", "warm")
ORDER = "ABBABAAB"
OVERHEAD = 4
LOSS = "0.5"
MAX_ENVIRONMENTAL_ATTEMPTS = 10
MAX_MINOR_FAULTS = 64
MAX_STDOUT_BYTES = 512 * 1024
MAX_STDERR_BYTES = 64 * 1024
MAX_THERMAL_GAP_S = 2.25
MAX_THERMAL_MARGIN_S = 2.25
MAX_CPU_TEMP_C = 85.0
MAX_DIMM_TEMP_C = 90.0
MALLOC_MMAP_THRESHOLD = "1073741824"
MALLOC_TRIM_THRESHOLD = "-1"

# Run as root only after opening a pidfd, then validate the process behind the
# descriptor before signaling it.  This closes the PID-reuse race inherent in
# a privileged ``kill PID`` between controller-side identity validation and
# signal delivery.
PIDFD_STOP_PROGRAM = r"""
import ctypes
import hashlib
import os
import signal
import sys

pid = int(sys.argv[1])
expected_tick = int(sys.argv[2])
expected_cmdline_sha256 = sys.argv[3]
libc = ctypes.CDLL(None, use_errno=True)
pidfd = libc.syscall(434, pid, 0)
if pidfd < 0:
    raise OSError(ctypes.get_errno(), 'pidfd_open')
try:
    stat_raw = open('/proc/%d/stat' % pid, 'rb').read()
    right = stat_raw.rfind(b')')
    fields = stat_raw[right + 2:].split() if right >= 0 else []
    cmdline_raw = open('/proc/%d/cmdline' % pid, 'rb').read()
    if (len(fields) < 20 or not fields[19].isdigit() or
            int(fields[19]) != expected_tick or
            hashlib.sha256(cmdline_raw).hexdigest() != expected_cmdline_sha256):
        raise SystemExit(73)
    if libc.syscall(424, pidfd, signal.SIGTERM, 0, 0) != 0:
        raise OSError(ctypes.get_errno(), 'pidfd_send_signal')
finally:
    os.close(pidfd)
""".strip()

SUPERSEDED_PRELAUNCH = Path(
    "/tmp/wh2-p32-dispatch-postopt-timing-a3f5c66-v1/prelaunch_receipts.json")
SUPERSEDED_PRELAUNCH_SHA256 = (
    "72904b0f1b4323222fa9dc0cff8030583f67b021ae8cfdf8f457f7f5d77dff02")
RECOVERY_ROOT = Path("/tmp/wh2-grouped-allk-finalists-0978602-v1")
RECOVERY_SUMMARY_SHA256 = (
    "c09d1e9b74a0e6cc38fea67d7038a3a77281736d7fd1be4d2867c48d6e3a2e72")
RECOVERY_DATA_MANIFEST_SHA256 = (
    "725f77bd96076507e47156c352b9d48f8e62b95dc2ddce888eefa4a6848d05a1")


def expected_external_receipts() -> Dict[str, Dict[str, str]]:
    return {
        "superseded_prelaunch_receipts.json": {
            "source": str(SUPERSEDED_PRELAUNCH),
            "sha256": SUPERSEDED_PRELAUNCH_SHA256,
        },
        "recovery_validated_summary.json": {
            "source": str(RECOVERY_ROOT / "validated_summary.json"),
            "sha256": RECOVERY_SUMMARY_SHA256,
        },
        "recovery_data_manifest.sha256": {
            "source": str(RECOVERY_ROOT / "data_manifest.sha256"),
            "sha256": RECOVERY_DATA_MANIFEST_SHA256,
        },
    }

ARMS: Dict[str, Dict[str, object]] = {
    "prod244": {
        "geometry": "frozen", "period": 244, "gf256_rows": 10,
        "gf16_rows": 2, "residue_skew": 0,
        "residue_schedule": "constant", "residue_hash_seed": "0x0",
        "residues_rotated": 0, "independent_extension_residues": 0,
        "independent_extension_seed_xor": "0x4e",
        "extension_hash_seed": "0x0", "grouped_rows": 0,
        "grouped_hash_seed": "0x0", "buckets": "auto",
        "packet_seed_multiplier": 1, "packet_seed_avalanche": 0,
        "odd_packet_peel_seed_xor": "0x0", "dense_rows_override": 0,
        "dense_two_anchor": 1, "dense_two_anchor_phase": 0,
        "mix_count": 2, "final_h_a_columns": 0,
    },
    "p32_r7": {
        "geometry": "shared-x", "period": 32, "gf256_rows": 10,
        "gf16_rows": 2, "residue_skew": 0,
        "residue_schedule": "constant", "residue_hash_seed": "0x0",
        "residues_rotated": 0, "independent_extension_residues": 0,
        "independent_extension_seed_xor": "0x4e",
        "extension_hash_seed": "0x0", "grouped_rows": 7,
        "grouped_hash_seed": "0xb7e15162", "buckets": "auto",
        "packet_seed_multiplier": 1, "packet_seed_avalanche": 0,
        "odd_packet_peel_seed_xor": "0x0", "dense_rows_override": 0,
        "dense_two_anchor": 1, "dense_two_anchor_phase": 0,
        "mix_count": 2, "final_h_a_columns": 12,
    },
}

UINT_RE = re.compile(r"0|[1-9][0-9]*\Z")
SINT_RE = re.compile(r"0|-?[1-9][0-9]*\Z")
HEX_RE = re.compile(r"0x(?:0|[1-9a-f][0-9a-f]*)\Z")
SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
TASK_DIR_RE = re.compile(r"[0-9]{4}\Z")
TXN_DIR_RE = re.compile(r"[0-9]{4}\.part(?:\.[1-9][0-9]*)?\Z")

CSV_FIELDS = (
    "N", "bb", "overhead", "schedule", "seed", "loss", "cache_state",
    "cycle", "slot", "arm", "period", "grouped_rows", "buckets_requested",
    "seed_attempt", "matrix_seed", "peel_seed", "preflight_result",
    "cell_class", "common_success", "result", "outcome_stable", "elapsed_ns",
    "saturated", "cpu_before", "cpu_after", "cpu_migrated", "minflt_delta",
    "majflt_delta", "fault_contaminated", "inactivated", "binary_def",
    "heavy_gain", "block_xors", "block_muladds", "build_ns", "peel_ns",
    "project_ns", "residual_ns", "backsub_ns", "joint_source_xors",
    "joint_marginal_xors", "joint_marginal_copies", "joint_active_deltas",
    "joint_scratch_bytes", "dual_source_columns", "source_bytes",
    "packet_payload_bytes", "intermediate_bytes", "solve_value_arena_bytes",
    "solve_value_eager_zero_bytes", "solve_value_commit_copy_bytes",
    "geometry", "gf256_rows", "gf16_rows", "residue_skew",
    "residue_schedule", "residue_hash_seed", "residues_rotated",
    "independent_extension_residues", "independent_extension_seed_xor",
    "extension_hash_seed", "grouped_hash_seed_exact", "packet_seed_multiplier",
    "packet_seed_avalanche", "odd_packet_peel_seed_xor",
    "dense_rows_override", "dense_two_anchor_exact", "dense_two_anchor_phase",
    "staircase_rows", "dense_rows", "heavy_rows", "source_hits", "field",
    "dense_identity_corner", "heavy_family", "mix_count",
    "rhs_route_expected", "rhs_route_actual",
)

THERMAL_FIELDS = (
    "utc", "monotonic_s", "cpu_busy_pct", "cpu_avg_mhz", "cpu_tctl_c",
    "dimm_i2c1_50_c", "dimm_i2c1_51_c", "dimm_i2c1_52_c",
    "dimm_i2c1_53_c", "dimm_i2c2_50_c", "dimm_i2c2_51_c",
    "dimm_i2c2_52_c", "dimm_i2c2_53_c", "dimm_read_errors", "load1",
    "load5", "load15", "edac_ce", "edac_ue",
)
DIMM_FIELDS = tuple(field for field in THERMAL_FIELDS
                    if field.startswith("dimm_i2c"))


class TimingError(RuntimeError):
    """Substantive campaign or evidence failure."""


class ContaminationError(TimingError):
    """Retryable, explicitly receipted environmental timing contamination."""


class BoundedProcessError(TimingError):
    """Child exceeded a declared output or wall-time bound."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds").replace(
        "+00:00", "Z")


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
    fd = os.open(str(path), os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


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
    try:
        os.link(part, path, follow_symlinks=False)
    except FileExistsError as exc:
        raise TimingError("refusing to replace %s" % path) from exc
    os.unlink(part)
    fsync_directory(path.parent)


def write_new(path: Path, data: bytes, mode: int = 0o444) -> None:
    if path.exists() or path.is_symlink():
        raise TimingError("refusing to replace %s" % path)
    atomic_write(path, data, mode)


def sealed_record(schema: str, payload: Mapping[str, object]) -> Dict[str, object]:
    if not schema or "schema" in payload or "self_sha256_excluding_field" in payload:
        raise TimingError("invalid sealed record")
    value: Dict[str, object] = {"schema": schema, **payload}
    value["self_sha256_excluding_field"] = sha256_bytes(canonical_json(value))
    return value


def verify_sealed(value: object, schema: str, name: str) -> Dict[str, object]:
    if not isinstance(value, dict) or value.get("schema") != schema:
        raise TimingError("%s schema mismatch" % name)
    claimed = value.get("self_sha256_excluding_field")
    if not isinstance(claimed, str) or SHA256_RE.fullmatch(claimed) is None:
        raise TimingError("%s self hash malformed" % name)
    unhashed = dict(value)
    del unhashed["self_sha256_excluding_field"]
    if sha256_bytes(canonical_json(unhashed)) != claimed:
        raise TimingError("%s self hash mismatch" % name)
    return value


def load_canonical(path: Path, name: str) -> Dict[str, object]:
    raw = path.read_bytes()
    try:
        value = json.loads(raw.decode("ascii"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise TimingError("%s is not canonical ASCII JSON" % name) from exc
    if not isinstance(value, dict) or canonical_json(value) != raw:
        raise TimingError("%s is not a canonical JSON object" % name)
    return value


def parse_uint(text: object, name: str, maximum: int = (1 << 64) - 1) -> int:
    if not isinstance(text, str) or UINT_RE.fullmatch(text) is None:
        raise TimingError("%s is not a canonical uint" % name)
    value = int(text)
    if value > maximum:
        raise TimingError("%s is outside its integer domain" % name)
    return value


def parse_sint(text: object, name: str) -> int:
    if not isinstance(text, str) or SINT_RE.fullmatch(text) is None:
        raise TimingError("%s is not a canonical signed integer" % name)
    return int(text)


def validate_utc_timestamp(value: object, name: str) -> str:
    if not isinstance(value, str) or not value.endswith("Z"):
        raise TimingError(name + " is not a UTC timestamp")
    try:
        parsed = datetime.fromisoformat(value[:-1] + "+00:00")
    except ValueError as exc:
        raise TimingError(name + " is malformed") from exc
    if parsed.tzinfo != timezone.utc:
        raise TimingError(name + " is not normalized to UTC")
    return value


def sanitized_environment(home: Path, *, allocator: bool) -> Dict[str, str]:
    result = {
        "HOME": str(home), "PATH": "/usr/bin:/bin", "LC_ALL": "C",
        "LANG": "C", "TZ": "UTC", "PYTHONDONTWRITEBYTECODE": "1",
    }
    if allocator:
        result["MALLOC_MMAP_THRESHOLD_"] = MALLOC_MMAP_THRESHOLD
        result["MALLOC_TRIM_THRESHOLD_"] = MALLOC_TRIM_THRESHOLD
    return result


def generate_tasks() -> List[Dict[str, object]]:
    pending: List[Tuple[str, Dict[str, object]]] = []
    for K in KS:
        for bb in WIDTHS:
            for seed_index, seed in enumerate(SEEDS):
                for schedule in SCHEDULES:
                    for cache_state in CACHE_STATES:
                        coordinate: Dict[str, object] = {
                            "K": K, "bb": bb, "seed_index": seed_index,
                            "seed": seed, "schedule": schedule,
                            "cache_state": cache_state, "control": "prod244",
                            "candidate": "p32_r7",
                        }
                        priority = sha256_bytes(
                            b"wirehair.wh2.p32.dispatch.timing.order.v2\0" +
                            canonical_json(coordinate))
                        pending.append((priority, coordinate))
    pending.sort(key=lambda item: item[0])
    tasks = []
    for job, (_priority, coordinate) in enumerate(pending):
        task = dict(coordinate)
        task_id = sha256_bytes(
            b"wirehair.wh2.p32.dispatch.timing.task.v2\0" +
            canonical_json(coordinate))[:24]
        task.update({"job": job, "task_id": task_id})
        tasks.append(task)
    return tasks


def expected_rhs_route(arm_name: str, K: int, bb: int) -> str:
    if arm_name == "prod244":
        return "streamed"
    if arm_name != "p32_r7":
        raise TimingError("unknown timing arm")
    if (bb >= 4096 and K >= 3200) or (bb >= 1280 and K >= 10000):
        return "joint-delta"
    return "streamed"


def command_for(design: Mapping[str, object], task: Mapping[str, object],
                *, cycle_index: Optional[int] = None,
                evict_bytes: Optional[int] = None) -> List[str]:
    root = Path(str(design["root"]))
    binary = root / "frozen/wirehair_v2_bench"
    control = ARMS["prod244"]
    candidate = ARMS["p32_r7"]
    command = [
        str(design["tools"]["taskset"]["path"]), "-c", str(design["core"]),
        str(design["tools"]["numactl"]["path"]),
        "--physcpubind=" + str(design["core"]),
        "--membind=" + str(design["numa_node"]), str(binary),
        "groupedtiming", "--N", str(task["K"]), "--bb", str(task["bb"]),
        "--overhead", str(OVERHEAD), "--control-geometry",
        str(control["geometry"]), "--control-period", str(control["period"]),
        "--control-grouped-rows", str(control["grouped_rows"]),
        "--control-buckets", str(control["buckets"]),
        "--candidate-geometry", str(candidate["geometry"]),
        "--candidate-period", str(candidate["period"]),
        "--candidate-grouped-rows", str(candidate["grouped_rows"]),
        "--candidate-buckets", str(candidate["buckets"]),
        "--evict-bytes", str(evict_bytes if evict_bytes is not None
                             else design["evict_bytes"]),
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


def _arm_preamble_expected(prefix: str, arm_name: str, task: Mapping[str, object]) -> Dict[str, str]:
    arm = ARMS[arm_name]
    result = {
        prefix + "_period": str(arm["period"]),
        prefix + "_grouped_rows": str(arm["grouped_rows"]),
        prefix + "_buckets": str(arm["buckets"]),
        prefix + "_grouped_hash_seed": str(arm["grouped_hash_seed"]),
        prefix + "_final_h_a_columns": str(arm["final_h_a_columns"]),
        prefix + "_geometry": str(arm["geometry"]),
        prefix + "_gf256_rows": str(arm["gf256_rows"]),
        prefix + "_gf16_rows": str(arm["gf16_rows"]),
        prefix + "_residue_skew": str(arm["residue_skew"]),
        prefix + "_residue_schedule": str(arm["residue_schedule"]),
        prefix + "_residue_hash_seed": str(arm["residue_hash_seed"]),
        prefix + "_residues_rotated": str(arm["residues_rotated"]),
        prefix + "_independent_extension_residues":
            str(arm["independent_extension_residues"]),
        prefix + "_independent_extension_seed_xor":
            str(arm["independent_extension_seed_xor"]),
        prefix + "_extension_hash_seed": str(arm["extension_hash_seed"]),
        prefix + "_grouped_hash_seed_exact": str(arm["grouped_hash_seed"]),
        prefix + "_packet_seed_multiplier": str(arm["packet_seed_multiplier"]),
        prefix + "_packet_seed_avalanche": str(arm["packet_seed_avalanche"]),
        prefix + "_odd_packet_peel_seed_xor":
            str(arm["odd_packet_peel_seed_xor"]),
        prefix + "_dense_rows_override": str(arm["dense_rows_override"]),
        prefix + "_dense_two_anchor_exact": str(arm["dense_two_anchor"]),
        prefix + "_dense_two_anchor_phase": str(arm["dense_two_anchor_phase"]),
        prefix + "_mix_count": str(arm["mix_count"]),
        prefix + "_rhs_route_expected": expected_rhs_route(
            arm_name, int(task["K"]), int(task["bb"])),
    }
    return result


def parse_grouped_output(raw: bytes, task: Mapping[str, object], core: int,
                         *, replacement_cycle: Optional[int] = None) -> ParsedOutput:
    if len(raw) > MAX_STDOUT_BYTES:
        raise TimingError("grouped timing output exceeds its sealed bound")
    if not raw.endswith(b"\n") or b"\r" in raw or b"\0" in raw:
        raise TimingError("grouped timing output is not canonical LF ASCII")
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
        raise TimingError("grouped timing preamble token malformed")
    preamble = dict(token.split("=", 1) for token in tokens)
    if len(preamble) != len(tokens):
        raise TimingError("grouped timing preamble has duplicate keys")
    expected = {
        "schema": "v2", "policy": "h12-q0-grouped", "timing_scope": "solve",
        "cycles": "1" if replacement_cycle is not None else "4",
        "order": ORDER, "discard_cycle": "0",
        "cycle_mode": "replacement" if replacement_cycle is not None else "full",
        "cycle_index": str(replacement_cycle) if replacement_cycle is not None else "all",
        "N": str(task["K"]), "bb": str(task["bb"]),
        "overhead": str(OVERHEAD), "loss": LOSS, "seed": str(task["seed"]),
        "schedule": str(task["schedule"]), "cache_state": str(task["cache_state"]),
        "overhead_stream": "salted", "eviction_prefaulted": "1",
        "gf256_rows": "10", "gf16_rows": "2", "dense_two_anchor": "1",
        "mix": "2", "payload": "distinct-packet-zero-v1",
        "payload_count": str(int(task["K"]) + OVERHEAD),
        "payload_bytes": str((int(task["K"]) + OVERHEAD) * int(task["bb"])),
        "payload_alignment": "64", "payload_prefaulted": "1",
        "system_build": "outside-timer",
        "tls_reapply": "full-per-slot-outside-timer",
        "allocator_tls_state": "preflight-warmed",
        "solve_value_storage": "owned-noinit", "solve_value_publish": "swap",
    }
    expected.update(_arm_preamble_expected("control", "prod244", task))
    expected.update(_arm_preamble_expected("candidate", "p32_r7", task))
    dynamic_fields = (
        "staircase_rows", "dense_rows", "heavy_rows", "source_hits", "field",
        "dense_identity_corner", "heavy_family", "mix_count",
    )
    exact_keys = set(expected) | {
        "evict_bytes", "trace_sha256", "cell_class", "common_success",
        "preflight_control_result", "preflight_candidate_result",
    }
    for prefix_name in ("control", "candidate"):
        exact_keys.update({
            prefix_name + "_attempt", prefix_name + "_matrix_seed",
            prefix_name + "_peel_seed",
            prefix_name + "_preflight_rhs_route",
        })
        exact_keys.update(prefix_name + "_" + field
                          for field in dynamic_fields)
    if set(preamble) != exact_keys:
        raise TimingError("grouped timing preamble key set changed")
    for key, value in expected.items():
        if preamble.get(key) != value:
            raise TimingError("preamble mismatch %s: %r != %r" %
                              (key, preamble.get(key), value))
    if parse_uint(preamble.get("evict_bytes"), "preamble evict bytes") < 4096:
        raise TimingError("preamble eviction allocation is too small")
    trace = preamble.get("trace_sha256", "")
    if SHA256_RE.fullmatch(trace) is None:
        raise TimingError("trace receipt is not a lowercase SHA256")
    arm_receipts: Dict[str, Dict[str, str]] = {}
    for prefix_name, arm_name in (("control", "prod244"),
                                  ("candidate", "p32_r7")):
        receipt = {
            "attempt": preamble.get(prefix_name + "_attempt", ""),
            "matrix_seed": preamble.get(prefix_name + "_matrix_seed", ""),
            "peel_seed": preamble.get(prefix_name + "_peel_seed", ""),
            "result": preamble.get("preflight_" + prefix_name + "_result", ""),
            "route": preamble.get(prefix_name + "_preflight_rhs_route", ""),
        }
        parse_uint(receipt["attempt"], prefix_name + " attempt")
        if receipt["attempt"] != "0" or HEX_RE.fullmatch(receipt["matrix_seed"]) is None or \
                HEX_RE.fullmatch(receipt["peel_seed"]) is None:
            raise TimingError(prefix_name + " graph receipt malformed or seed-fixed")
        parse_sint(receipt["result"], prefix_name + " preflight result")
        if receipt["route"] not in ("not-reached", "streamed", "dual", "joint-delta"):
            raise TimingError(prefix_name + " preflight RHS route malformed")
        if int(receipt["result"]) == 0 and receipt["route"] != \
                expected_rhs_route(arm_name, int(task["K"]), int(task["bb"])):
            raise TimingError(prefix_name + " successful preflight used wrong RHS route")
        arm_receipts[prefix_name] = receipt
    for key in ("attempt", "matrix_seed", "peel_seed"):
        if arm_receipts["control"][key] != arm_receipts["candidate"][key]:
            raise TimingError("paired arms do not share the base graph receipt")
    for field in dynamic_fields:
        left = preamble.get("control_" + field)
        right = preamble.get("candidate_" + field)
        parse_uint(left, "control " + field)
        parse_uint(right, "candidate " + field)
        if left != right:
            raise TimingError("paired arms differ in base graph field " + field)
    reader = csv.DictReader(io.StringIO("\n".join(lines[1:]) + "\n"))
    if tuple(reader.fieldnames or ()) != CSV_FIELDS:
        raise TimingError("grouped timing CSV schema mismatch")
    rows = tuple(dict(row) for row in reader)
    if len(rows) != expected_rows:
        raise TimingError("grouped timing row count changed")
    control_success = int(arm_receipts["control"]["result"]) == 0
    candidate_success = int(arm_receipts["candidate"]["result"]) == 0
    expected_class = (
        "common-success" if control_success and candidate_success else
        "control-only" if control_success else
        "candidate-only" if candidate_success else "common-failure")
    if preamble.get("cell_class") != expected_class or \
            preamble.get("common_success") != ("1" if expected_class == "common-success" else "0"):
        raise TimingError("preamble outcome classification mismatch")
    work_fields = (
        "result", "inactivated", "binary_def", "heavy_gain", "block_xors",
        "block_muladds", "joint_source_xors", "joint_marginal_xors",
        "joint_marginal_copies", "joint_active_deltas", "joint_scratch_bytes",
        "dual_source_columns", "source_bytes", "packet_payload_bytes",
        "intermediate_bytes", "rhs_route_actual",
    )
    signatures: Dict[str, set[Tuple[str, ...]]] = {"control": set(), "candidate": set()}
    elapsed = {"control": 0, "candidate": 0}
    contaminations: List[str] = []
    first_cycle = replacement_cycle if replacement_cycle is not None else 0
    for index, row in enumerate(rows):
        if None in row or tuple(row) != CSV_FIELDS or \
                any(value is None for value in row.values()):
            raise TimingError("grouped timing row width changed")
        cycle = first_cycle + index // 8
        slot = index % 8
        prefix_name = "control" if ORDER[slot] == "A" else "candidate"
        arm_name = "prod244" if prefix_name == "control" else "p32_r7"
        arm = ARMS[arm_name]
        receipt = arm_receipts[prefix_name]
        exact = {
            "N": str(task["K"]), "bb": str(task["bb"]),
            "overhead": str(OVERHEAD), "schedule": str(task["schedule"]),
            "seed": str(task["seed"]), "loss": LOSS,
            "cache_state": str(task["cache_state"]), "cycle": str(cycle),
            "slot": str(slot), "arm": prefix_name, "period": str(arm["period"]),
            "grouped_rows": str(arm["grouped_rows"]),
            "buckets_requested": str(arm["buckets"]),
            "seed_attempt": receipt["attempt"], "matrix_seed": receipt["matrix_seed"],
            "peel_seed": receipt["peel_seed"], "preflight_result": receipt["result"],
            "cell_class": expected_class,
            "common_success": "1" if expected_class == "common-success" else "0",
            "source_bytes": str(int(task["K"]) * int(task["bb"])),
            "packet_payload_bytes": str((int(task["K"]) + OVERHEAD) * int(task["bb"])),
            "geometry": str(arm["geometry"]), "gf256_rows": str(arm["gf256_rows"]),
            "gf16_rows": str(arm["gf16_rows"]),
            "residue_skew": str(arm["residue_skew"]),
            "residue_schedule": str(arm["residue_schedule"]),
            "residue_hash_seed": str(arm["residue_hash_seed"]),
            "residues_rotated": str(arm["residues_rotated"]),
            "independent_extension_residues": str(arm["independent_extension_residues"]),
            "independent_extension_seed_xor": str(arm["independent_extension_seed_xor"]),
            "extension_hash_seed": str(arm["extension_hash_seed"]),
            "grouped_hash_seed_exact": str(arm["grouped_hash_seed"]),
            "packet_seed_multiplier": str(arm["packet_seed_multiplier"]),
            "packet_seed_avalanche": str(arm["packet_seed_avalanche"]),
            "odd_packet_peel_seed_xor": str(arm["odd_packet_peel_seed_xor"]),
            "dense_rows_override": str(arm["dense_rows_override"]),
            "dense_two_anchor_exact": str(arm["dense_two_anchor"]),
            "dense_two_anchor_phase": str(arm["dense_two_anchor_phase"]),
            "mix_count": str(arm["mix_count"]),
            "rhs_route_expected": expected_rhs_route(
                arm_name, int(task["K"]), int(task["bb"])),
        }
        for field in dynamic_fields[:-1]:
            exact[field] = preamble[prefix_name + "_" + field]
        for key, value in exact.items():
            if row.get(key) != value:
                raise TimingError("row %d mismatch %s: %r != %r" %
                                  (index, key, row.get(key), value))
        if row.get("outcome_stable") != "1" or row.get("result") != receipt["result"]:
            raise TimingError("physical solve changed from arm preflight")
        if row.get("rhs_route_actual") != receipt["route"]:
            raise TimingError("physical solve RHS route changed from preflight")
        if int(receipt["result"]) == 0 and row["rhs_route_actual"] != \
                exact["rhs_route_expected"]:
            raise TimingError("successful physical solve used wrong RHS route")
        minflt = parse_sint(row.get("minflt_delta"), "minor faults")
        majflt = parse_sint(row.get("majflt_delta"), "major faults")
        expected_fault = "-1" if minflt < 0 or majflt < 0 else \
            "1" if minflt or majflt else "0"
        if row.get("fault_contaminated") != expected_fault:
            raise TimingError("signed fault contamination receipt mismatch")
        if row.get("saturated") != "0":
            contaminations.append("row%d:saturated" % index)
        if row.get("cpu_before") != str(core) or row.get("cpu_after") != str(core) or \
                row.get("cpu_migrated") != "0":
            contaminations.append("row%d:cpu-migration" % index)
        if minflt < 0 or minflt > MAX_MINOR_FAULTS or majflt != 0:
            contaminations.append("row%d:faults:%d:%d" % (index, minflt, majflt))
        if parse_uint(row.get("elapsed_ns"), "elapsed ns") == 0:
            raise TimingError("physical solve has zero elapsed time")
        for field in (
                "inactivated", "binary_def", "heavy_gain", "block_xors",
                "block_muladds", "build_ns", "peel_ns", "project_ns",
                "residual_ns", "backsub_ns", "joint_source_xors",
                "joint_marginal_xors", "joint_marginal_copies",
                "joint_active_deltas", "joint_scratch_bytes",
                "dual_source_columns", "source_bytes", "packet_payload_bytes",
                "intermediate_bytes", "solve_value_arena_bytes",
                "solve_value_eager_zero_bytes", "solve_value_commit_copy_bytes",
                "staircase_rows", "dense_rows", "heavy_rows", "source_hits",
                "field", "dense_identity_corner", "heavy_family"):
            parse_uint(row.get(field), field)
        if row["solve_value_arena_bytes"] != row["intermediate_bytes"] or \
                row["solve_value_eager_zero_bytes"] != "0" or \
                row["solve_value_commit_copy_bytes"] != "0":
            raise TimingError("owned no-init solve arena receipt changed")
        signatures[prefix_name].add(tuple(row[field] for field in work_fields))
        if cycle != 0:
            elapsed[prefix_name] += int(row["elapsed_ns"])
    if any(len(values) != 1 for values in signatures.values()):
        raise TimingError("deterministic work changed within an arm")
    return ParsedOutput(
        preamble=preamble, rows=rows, stdout_sha256=sha256_bytes(raw),
        trace_sha256=trace, cell_class=expected_class,
        common_success=expected_class == "common-success",
        timed_control_ns=elapsed["control"], timed_candidate_ns=elapsed["candidate"],
        work_signatures={key: next(iter(value)) for key, value in signatures.items()},
        contaminations=tuple(sorted(set(contaminations))),
    )


def run_bounded(argv: Sequence[str], environment: Mapping[str, str],
                timeout_s: float, stdout_limit: int = MAX_STDOUT_BYTES,
                stderr_limit: int = MAX_STDERR_BYTES,
                cwd: Optional[Path] = None) -> Tuple[int, bytes, bytes]:
    """Run a child without ever buffering beyond the declared output caps."""
    if timeout_s <= 0 or stdout_limit < 0 or stderr_limit < 0:
        raise TimingError("bounded child limits are invalid")
    process = subprocess.Popen(
        list(argv), cwd=str(cwd) if cwd is not None else None,
        env=dict(environment), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        start_new_session=True, close_fds=True)
    if process.stdout is None or process.stderr is None:
        raise TimingError("bounded child pipes are unavailable")
    os.set_blocking(process.stdout.fileno(), False)
    os.set_blocking(process.stderr.fileno(), False)
    selector = selectors.DefaultSelector()
    selector.register(process.stdout, selectors.EVENT_READ, ("stdout", stdout_limit))
    selector.register(process.stderr, selectors.EVENT_READ, ("stderr", stderr_limit))
    buffers = {"stdout": bytearray(), "stderr": bytearray()}
    deadline = time.monotonic() + timeout_s
    failure: Optional[str] = None
    try:
        while selector.get_map() or process.poll() is None:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                failure = "timeout"
                break
            if selector.get_map():
                events = selector.select(min(remaining, 0.25))
            else:
                try:
                    process.wait(timeout=min(remaining, 0.25))
                except subprocess.TimeoutExpired:
                    pass
                events = []
            for key, _mask in events:
                name, limit = key.data
                try:
                    block = os.read(key.fileobj.fileno(), 65536)
                except BlockingIOError:
                    continue
                if not block:
                    selector.unregister(key.fileobj)
                    continue
                if len(buffers[name]) + len(block) > limit:
                    allowed = max(0, limit - len(buffers[name]))
                    buffers[name].extend(block[:allowed])
                    failure = name + "-limit"
                    break
                buffers[name].extend(block)
            if failure is not None:
                break
        if failure is not None:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        returncode = process.wait(timeout=10)
    finally:
        selector.close()
        process.stdout.close()
        process.stderr.close()
    if failure is not None:
        raise BoundedProcessError("child exceeded " + failure)
    return returncode, bytes(buffers["stdout"]), bytes(buffers["stderr"])


def parse_cpu_list(text: str) -> Tuple[int, ...]:
    if not text or text.strip() != text:
        raise TimingError("CPU affinity list is not canonical")
    values: List[int] = []
    for token in text.split(","):
        if "-" in token:
            pieces = token.split("-")
            if len(pieces) != 2 or not all(piece.isdigit() for piece in pieces):
                raise TimingError("CPU affinity range malformed")
            low, high = map(int, pieces)
            if low >= high:
                raise TimingError("CPU affinity range is not increasing")
            values.extend(range(low, high + 1))
        elif token.isdigit():
            values.append(int(token))
        else:
            raise TimingError("CPU affinity value malformed")
    if values != sorted(set(values)):
        raise TimingError("CPU affinity values duplicate or reorder")
    return tuple(values)


def topology_record(core: int, numa_node: int) -> Dict[str, object]:
    cpu = Path("/sys/devices/system/cpu") / ("cpu%d" % core)
    if not cpu.is_dir():
        raise TimingError("requested CPU is unavailable")
    siblings_text = (cpu / "topology/thread_siblings_list").read_text(
        encoding="ascii").strip()
    siblings = parse_cpu_list(siblings_text)
    if core not in siblings or len(siblings) != 2:
        raise TimingError("timing topology is not one two-thread core")
    llc_dirs = sorted(cpu.glob("cache/index*"),
                      key=lambda path: int(path.name.removeprefix("index")))
    llc = None
    for cache in llc_dirs:
        if (cache / "type").read_text(encoding="ascii").strip() == "Unified":
            llc = cache
    if llc is None:
        raise TimingError("unified LLC topology is unavailable")
    llc_text = (llc / "shared_cpu_list").read_text(encoding="ascii").strip()
    llc_cpus = parse_cpu_list(llc_text)
    nodes = sorted(path.name for path in cpu.glob("node[0-9]*"))
    if nodes != ["node%d" % numa_node]:
        raise TimingError("CPU is outside requested NUMA node")
    governor_path = cpu / "cpufreq/scaling_governor"
    preference_path = cpu / "cpufreq/energy_performance_preference"
    return {
        "core": core, "numa_node": numa_node,
        "thread_siblings_list": siblings_text, "siblings": list(siblings),
        "sibling": next(value for value in siblings if value != core),
        "llc_shared_cpu_list": llc_text, "llc_shared_cpus": list(llc_cpus),
        "governor": governor_path.read_text(encoding="ascii").strip()
            if governor_path.exists() else "unavailable",
        "energy_performance_preference":
            preference_path.read_text(encoding="ascii").strip()
            if preference_path.exists() else "unavailable",
    }


def _process_start_tick(stat_raw: bytes) -> int:
    right = stat_raw.rfind(b")")
    if right < 0:
        raise TimingError("process stat comm field malformed")
    fields = stat_raw[right + 2:].split()
    if len(fields) < 20 or not fields[19].isdigit():
        raise TimingError("process start tick unavailable")
    return int(fields[19])


def capture_process_identity(pid: int, expected_cpu: int, csv_path: Path,
                             *, proc_root: Path = Path("/proc"),
                             boot_id_path: Path = Path(
                                 "/proc/sys/kernel/random/boot_id")) -> Dict[str, object]:
    if pid <= 1:
        raise TimingError("sampler PID outside its domain")
    process = proc_root / str(pid)
    stat_raw = (process / "stat").read_bytes()
    cmdline_raw = (process / "cmdline").read_bytes()
    if not cmdline_raw.endswith(b"\0"):
        raise TimingError("sampler cmdline is not NUL terminated")
    argv_raw = cmdline_raw[:-1].split(b"\0")
    try:
        argv = [item.decode("ascii") for item in argv_raw]
    except UnicodeDecodeError as exc:
        raise TimingError("sampler cmdline is not ASCII") from exc
    status = (process / "status").read_text(encoding="ascii")
    match = re.search(r"^Cpus_allowed_list:\s*(\S+)\s*$", status, re.MULTILINE)
    if match is None:
        raise TimingError("sampler affinity receipt unavailable")
    affinity = parse_cpu_list(match.group(1))
    if affinity != (expected_cpu,):
        raise TimingError("sampler is not pinned to its sole CPU")
    csv_stat = csv_path.stat()
    if not csv_path.is_file() or csv_stat.st_size == 0:
        raise TimingError("sampler CSV has no baseline")
    boot_id = boot_id_path.read_text(encoding="ascii").strip()
    if re.fullmatch(r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}",
                    boot_id) is None:
        raise TimingError("boot ID malformed")
    return {
        "pid": pid, "start_tick": _process_start_tick(stat_raw),
        "boot_id": boot_id, "cmdline": argv,
        "cmdline_sha256": sha256_bytes(cmdline_raw),
        "affinity": list(affinity), "csv_device": csv_stat.st_dev,
        "csv_inode": csv_stat.st_ino,
    }


def process_identity_matches(identity: Mapping[str, object], expected_cpu: int,
                             csv_path: Path, *, proc_root: Path = Path("/proc"),
                             boot_id_path: Path = Path(
                                 "/proc/sys/kernel/random/boot_id")) -> bool:
    try:
        current = capture_process_identity(
            int(identity["pid"]), expected_cpu, csv_path,
            proc_root=proc_root, boot_id_path=boot_id_path)
    except (OSError, TimingError, KeyError, TypeError, ValueError):
        return False
    return current == dict(identity)


def sole_i2c_readers(fuser_path: Path, sudo_path: Path) -> Tuple[int, ...]:
    result = subprocess.run(
        (str(sudo_path), "-n", str(fuser_path), "/dev/i2c-1", "/dev/i2c-2"),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
        env=sanitized_environment(Path("/tmp"), allocator=False))
    if result.returncode not in (0, 1):
        raise TimingError("cannot inspect I2C reader ownership")
    if result.returncode == 1 and result.stdout:
        raise TimingError("fuser returned contradictory no-reader output")
    return tuple(sorted(set(int(value) for value in re.findall(rb"[0-9]+", result.stdout))))


def _parse_thermal_csv(raw: bytes) -> Tuple[Dict[str, str], ...]:
    if not raw.endswith(b"\n") or b"\0" in raw:
        raise TimingError("thermal CSV is not canonical complete text")
    try:
        text = raw.decode("ascii")
    except UnicodeDecodeError as exc:
        raise TimingError("thermal CSV is not ASCII") from exc
    reader = csv.DictReader(io.StringIO(text))
    if tuple(reader.fieldnames or ()) != THERMAL_FIELDS:
        raise TimingError("thermal CSV schema mismatch")
    rows = tuple(dict(row) for row in reader)
    if not rows:
        raise TimingError("thermal CSV has no samples")
    previous = -math.inf
    previous_utc: Optional[datetime] = None
    for row in rows:
        if None in row or tuple(row) != THERMAL_FIELDS or \
                any(value is None for value in row.values()):
            raise TimingError("thermal CSV row width changed")
        try:
            monotonic = float(row["monotonic_s"])
        except ValueError as exc:
            raise TimingError("thermal monotonic timestamp malformed") from exc
        if not math.isfinite(monotonic) or monotonic <= previous:
            raise TimingError("thermal timestamps are not strictly increasing")
        previous = monotonic
        validate_utc_timestamp(row["utc"], "thermal UTC sample")
        wall_time = datetime.fromisoformat(row["utc"][:-1] + "+00:00")
        if previous_utc is not None and wall_time <= previous_utc:
            raise TimingError("thermal UTC samples are not strictly increasing")
        previous_utc = wall_time
        for field in ("cpu_busy_pct", "cpu_avg_mhz", "load1", "load5", "load15"):
            try:
                value = float(row[field])
            except ValueError as exc:
                raise TimingError("thermal utilization field missing") from exc
            if not math.isfinite(value) or value < 0.0:
                raise TimingError("thermal utilization field implausible")
            if field == "cpu_busy_pct" and value > 100.0:
                raise TimingError("thermal CPU utilization exceeds 100 percent")
            if field == "cpu_avg_mhz" and not 100.0 <= value <= 10000.0:
                raise TimingError("thermal CPU frequency implausible")
        for field in ("cpu_tctl_c", *DIMM_FIELDS):
            try:
                temperature = float(row[field])
            except ValueError as exc:
                raise TimingError("thermal temperature missing") from exc
            if not math.isfinite(temperature) or not -40.0 < temperature < 130.0:
                raise TimingError("thermal temperature implausible")
        for field in ("dimm_read_errors", "edac_ce", "edac_ue"):
            parse_uint(row[field], "thermal " + field)
    return rows


def stable_file_bytes(path: Path, attempts: int = 20) -> bytes:
    for _attempt in range(attempts):
        before = path.stat()
        raw = path.read_bytes()
        after = path.stat()
        if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) == \
                (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns) and \
                len(raw) == after.st_size:
            return raw
        time.sleep(0.02)
    raise TimingError("file changed while being sealed")


def latest_thermal_time(path: Path) -> float:
    return float(_parse_thermal_csv(stable_file_bytes(path))[-1]["monotonic_s"])


def validate_thermal_interval(raw: bytes, start_s: float,
                              benchmark_end_s: float) -> Dict[str, object]:
    if not math.isfinite(start_s) or not math.isfinite(benchmark_end_s) or \
            benchmark_end_s <= start_s:
        raise TimingError("thermal benchmark interval malformed")
    rows = _parse_thermal_csv(raw)
    times = [float(row["monotonic_s"]) for row in rows]
    if times[0] > start_s or times[-1] < benchmark_end_s:
        raise TimingError("thermal samples do not bracket benchmark end")
    if start_s - times[0] > MAX_THERMAL_MARGIN_S or \
            times[-1] - benchmark_end_s > MAX_THERMAL_MARGIN_S:
        raise TimingError("thermal coverage margin exceeds gate")
    gaps = [right - left for left, right in zip(times, times[1:])]
    if gaps and max(gaps) > MAX_THERMAL_GAP_S:
        raise TimingError("thermal cadence gap exceeds gate")
    if any(int(row["dimm_read_errors"]) != 0 for row in rows):
        raise TimingError("DIMM read error during timing")
    ce = [int(row["edac_ce"]) for row in rows]
    ue = [int(row["edac_ue"]) for row in rows]
    if any(value != ce[0] for value in ce) or \
            any(value != ue[0] for value in ue):
        raise TimingError("EDAC counters changed during timing")
    cpu_max = max(float(row["cpu_tctl_c"]) for row in rows)
    dimm_max = {field: max(float(row[field]) for row in rows)
                for field in DIMM_FIELDS}
    if cpu_max > MAX_CPU_TEMP_C or max(dimm_max.values()) > MAX_DIMM_TEMP_C:
        raise TimingError("thermal temperature gate exceeded")
    return {
        "sample_count": len(rows), "first_monotonic_s": times[0],
        "last_monotonic_s": times[-1], "post_end_margin_s":
            times[-1] - benchmark_end_s,
        "start_margin_s": start_s - times[0],
        "max_gap_s": max(gaps) if gaps else 0.0, "cpu_max_c": cpu_max,
        "dimm_max_c": dimm_max, "dimm_read_errors": 0,
        "edac_ce_delta": 0, "edac_ue_delta": 0,
    }


def verify_terminal_thermal(path: Path, receipt: Mapping[str, object]) -> bytes:
    if path.is_symlink() or not path.is_file():
        raise TimingError("terminal thermal CSV is missing or unsafe")
    raw = stable_file_bytes(path)
    stat = path.stat()
    exact = {
        "thermal_csv_sha256": sha256_bytes(raw),
        "thermal_csv_size": len(raw), "thermal_csv_device": stat.st_dev,
        "thermal_csv_inode": stat.st_ino,
    }
    for key, value in exact.items():
        if receipt.get(key) != value:
            raise TimingError("terminal thermal receipt mismatch: " + key)
    summary = validate_thermal_interval(
        raw, float(receipt["timing_start_monotonic_s"]),
        float(receipt["benchmark_end_monotonic_s"]))
    if receipt.get("thermal_summary") != summary:
        raise TimingError("terminal thermal summary does not replay")
    return raw


def recover_interrupted_transactions(root: Path) -> List[Dict[str, object]]:
    transactions = root / ".transactions"
    interrupted = root / "interrupted"
    recovered: List[Dict[str, object]] = []
    entries = sorted(transactions.iterdir())
    if len(entries) > len(generate_tasks()):
        raise TimingError("transaction staging count exceeds task bound")
    for entry in entries:
        if entry.is_symlink() or not entry.is_dir() or TXN_DIR_RE.fullmatch(entry.name) is None:
            raise TimingError("unsafe interrupted transaction entry")
        total_bytes = 0
        files = {}
        for path in sorted(entry.iterdir()):
            if path.is_symlink() or not path.is_file():
                raise TimingError("interrupted transaction contains unsafe entry")
            stat = path.stat()
            total_bytes += stat.st_size
            if total_bytes > MAX_STDOUT_BYTES + MAX_STDERR_BYTES + 1024 * 1024:
                raise TimingError("interrupted transaction exceeds recovery bound")
            files[path.name] = {"size": stat.st_size, "sha256": sha256_file(path)}
        destination = interrupted / (entry.name + ".recovered.%d" % time.monotonic_ns())
        if destination.exists():
            raise TimingError("interrupted recovery destination collision")
        os.replace(entry, destination)
        os.chmod(destination, 0o555)
        fsync_directory(transactions)
        fsync_directory(interrupted)
        recovered.append({"source": entry.name, "archive": destination.name,
                          "files": files})
    return recovered


def terminal_task_bindings(root: Path) -> Dict[int, str]:
    bindings: Dict[int, str] = {}
    terminals = sorted((root / "segments").glob("segment*.terminal.json"))
    for index, path in enumerate(terminals):
        if path.name != "segment%04d.terminal.json" % index:
            raise TimingError("terminal sequence has a gap")
        terminal = verify_sealed(
            load_canonical(path, "segment terminal"),
            "wirehair.wh2.p32_dispatch.segment_terminal.v2",
            "segment terminal")
        completed = terminal.get("completed_tasks")
        if not isinstance(completed, list):
            raise TimingError("segment completed-task binding malformed")
        for item in completed:
            if not isinstance(item, dict) or not isinstance(item.get("job"), int) or \
                    SHA256_RE.fullmatch(str(item.get("receipt_sha256", ""))) is None:
                raise TimingError("segment task binding malformed")
            job = item["job"]
            if job in bindings:
                raise TimingError("task bound by multiple terminal segments")
            bindings[job] = item["receipt_sha256"]
    return bindings


def recover_unbound_task_entries(root: Path, bindings: Mapping[int, str]) -> List[Dict[str, object]]:
    recovered: List[Dict[str, object]] = []
    for entry in sorted((root / "ledger").iterdir()):
        if entry.is_symlink() or not entry.is_dir() or TASK_DIR_RE.fullmatch(entry.name) is None:
            raise TimingError("unsafe task ledger entry during resume")
        job = int(entry.name)
        if job in bindings:
            receipt_path = entry / "receipt.json"
            if not receipt_path.is_file() or sha256_file(receipt_path) != bindings[job]:
                raise TimingError("terminal-bound task receipt changed")
            continue
        files = {}
        total = 0
        for path in sorted(entry.iterdir()):
            if path.is_symlink() or not path.is_file():
                raise TimingError("unbound task ledger contains unsafe entry")
            total += path.stat().st_size
            if total > MAX_STDOUT_BYTES + MAX_STDERR_BYTES + 1024 * 1024:
                raise TimingError("unbound task ledger exceeds recovery bound")
            files[path.name] = {"size": path.stat().st_size,
                                "sha256": sha256_file(path)}
        destination = root / "interrupted" / (
            "task-%04d.unbound.recovered.%d" % (job, time.monotonic_ns()))
        os.chmod(entry, 0o700)
        os.replace(entry, destination)
        os.chmod(destination, 0o555)
        fsync_directory(root / "ledger")
        fsync_directory(root / "interrupted")
        recovered.append({"source": entry.name, "archive": destination.name,
                          "files": files, "reason": "no terminal thermal binding"})
    return recovered


def adopt_unbound_interrupted_archives(root: Path) -> List[Dict[str, object]]:
    referenced: set[str] = set()
    for path in sorted((root / "segments").glob("segment*.terminal.json")):
        terminal = verify_sealed(
            load_canonical(path, "segment terminal"),
            "wirehair.wh2.p32_dispatch.segment_terminal.v2",
            "segment terminal")
        for item in terminal.get("recovered_interrupted_transactions", []):
            if isinstance(item, dict) and isinstance(item.get("archive"), str):
                referenced.add(item["archive"])
    records = []
    for directory in sorted((root / "interrupted").iterdir()):
        if directory.name in referenced:
            continue
        if directory.is_symlink() or not directory.is_dir():
            raise TimingError("unsafe unbound interrupted archive")
        files = {}
        total = 0
        for path in sorted(directory.iterdir()):
            if path.is_symlink() or not path.is_file():
                raise TimingError("unbound interrupted archive contains unsafe entry")
            total += path.stat().st_size
            if total > MAX_STDOUT_BYTES + MAX_STDERR_BYTES + 1024 * 1024:
                raise TimingError("unbound interrupted archive exceeds bound")
            files[path.name] = {"size": path.stat().st_size,
                                "sha256": sha256_file(path)}
        records.append({"source": "prior-controller-interruption",
                        "archive": directory.name, "files": files,
                        "reason": "archive created before terminal receipt"})
    return records


def begin_task_transaction(root: Path, task: Mapping[str, object]) -> Path:
    job = int(task["job"])
    staging = root / ".transactions" / ("%04d.part.%d" % (job, os.getpid()))
    if staging.exists() or (root / "ledger" / ("%04d" % job)).exists():
        raise TimingError("task transaction already exists")
    staging.mkdir(mode=0o700)
    write_new(staging / "intent.json", canonical_json(sealed_record(
        "wirehair.wh2.p32_dispatch.intent.v2", {
            "created_utc": utc_now(), "job": job, "task_id": task["task_id"],
            "task_sha256": sha256_bytes(canonical_json(task)),
        })))
    fsync_directory(staging.parent)
    return staging


def commit_task_transaction(root: Path, task: Mapping[str, object], raw: bytes,
                            stderr: bytes, receipt: Mapping[str, object],
                            staging: Optional[Path] = None) -> Path:
    job = int(task["job"])
    final = root / "ledger" / ("%04d" % job)
    if final.exists():
        raise TimingError("task ledger entry already exists")
    if staging is None:
        staging = begin_task_transaction(root, task)
    if staging.parent != root / ".transactions" or not staging.is_dir() or \
            TXN_DIR_RE.fullmatch(staging.name) is None:
        raise TimingError("task staging identity is invalid")
    try:
        write_new(staging / "stdout.csv", raw)
        write_new(staging / "stderr.bin", stderr)
        write_new(staging / "receipt.json", canonical_json(receipt))
        os.replace(staging, final)
        os.chmod(final, 0o555)
        fsync_directory(final.parent)
    except BaseException:
        if staging.exists():
            os.chmod(staging, 0o700)
        raise
    return final


def resolve_tool(name: str) -> Path:
    value = shutil.which(name)
    if value is None:
        raise TimingError("required tool unavailable: " + name)
    path = Path(value).resolve()
    if not path.is_file() or not os.access(path, os.X_OK):
        raise TimingError("required tool is not executable: %s" % path)
    return path


def command_capture(argv: Sequence[str], *, cwd: Optional[Path] = None,
                    environment: Optional[Mapping[str, str]] = None) -> str:
    child_environment = dict(environment) if environment is not None else \
        sanitized_environment(Path("/tmp"), allocator=False)
    returncode, stdout, stderr = run_bounded(
        argv, child_environment, 60.0, 4 * 1024 * 1024,
        1024 * 1024, cwd=cwd)
    if returncode != 0 or stderr:
        raise TimingError("command failed exit=%d argv=%r stderr=%r" %
                          (returncode, list(argv), stderr[:1000]))
    try:
        return stdout.decode("utf-8").strip()
    except UnicodeDecodeError as exc:
        raise TimingError("command output is not UTF-8") from exc


def _git_value(git: Path, repo: Path, *args: str) -> str:
    return command_capture((str(git), "-C", str(repo), *args))


def _copy_frozen(source: Path, destination: Path, mode: int) -> None:
    if destination.exists():
        raise TimingError("frozen destination already exists")
    shutil.copyfile(source, destination)
    os.chmod(destination, mode)


def _build_frozen_binary(repo: Path, head: str, workspace: Path, staging: Path,
                         tools: Mapping[str, Path], jobs: int,
                         c_compiler: Path, cxx_compiler: Path) -> Dict[str, object]:
    source = workspace / "source"
    build = workspace / "build"
    home = workspace / "home"
    home.mkdir(parents=True)
    environment = sanitized_environment(home, allocator=False)
    add = subprocess.run(
        (str(tools["git"]), "-C", str(repo), "worktree", "add", "--detach",
         str(source), head), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        check=False, env=environment)
    if add.returncode != 0:
        raise TimingError("detached build worktree failed: %r" % add.stderr[:1000])
    try:
        if _git_value(tools["git"], source, "rev-parse", "HEAD") != head or \
                _git_value(tools["git"], source, "status", "--porcelain",
                           "--untracked-files=all"):
            raise TimingError("detached build worktree identity is not exact")
        generator = "Ninja" if "ninja" in tools else "Unix Makefiles"
        configure = (
            str(tools["cmake"]), "-S", str(source), "-B", str(build),
            "-G", generator, "-DCMAKE_BUILD_TYPE=Release", "-DBUILD_TESTS=ON",
            "-DBUILD_CODEC_V2=ON", "-DWIREHAIR_BUILD_BENCHMARKS=ON",
            "-DMARCH_NATIVE=OFF", "-DWIREHAIR_STRICT_WARNINGS=ON",
            "-DWH_LTO=OFF", "-DWH_PGO_MODE=OFF",
            "-DCMAKE_EXPORT_COMPILE_COMMANDS=ON",
            "-DCMAKE_C_COMPILER=" + str(c_compiler),
            "-DCMAKE_CXX_COMPILER=" + str(cxx_compiler),
        )
        begin = time.monotonic()
        rc, stdout, stderr = run_bounded(
            configure, environment, 600, 4 * 1024 * 1024, 4 * 1024 * 1024)
        configure_duration = time.monotonic() - begin
        write_new(staging / "provenance/configure.stdout", stdout)
        write_new(staging / "provenance/configure.stderr", stderr)
        if rc != 0:
            raise TimingError("frozen configure failed")
        build_command = (
            str(tools["cmake"]), "--build", str(build), "--target",
            "wirehair_v2_bench", "--parallel", str(jobs),
        )
        begin = time.monotonic()
        rc, stdout, stderr = run_bounded(
            build_command, environment, 1800, 8 * 1024 * 1024,
            8 * 1024 * 1024)
        build_duration = time.monotonic() - begin
        write_new(staging / "provenance/build.stdout", stdout)
        write_new(staging / "provenance/build.stderr", stderr)
        if rc != 0:
            raise TimingError("frozen benchmark build failed")
        if _git_value(tools["git"], source, "status", "--porcelain",
                      "--untracked-files=all"):
            raise TimingError("out-of-tree build modified frozen source")
        binary = build / "codec/wirehair_v2_bench"
        if not binary.is_file() or not os.access(binary, os.X_OK):
            raise TimingError("frozen benchmark binary missing")
        _copy_frozen(binary, staging / "frozen/wirehair_v2_bench", 0o555)
        for source_path, destination_name in (
                (build / "CMakeCache.txt", "CMakeCache.txt"),
                (build / "compile_commands.json", "compile_commands.json")):
            if not source_path.is_file():
                raise TimingError("build provenance file missing")
            _copy_frozen(source_path, staging / "provenance" / destination_name,
                         0o444)
        return {
            "configure_argv": list(configure),
            "build_argv": list(build_command), "environment": environment,
            "configure_duration_s": configure_duration,
            "build_duration_s": build_duration, "generator": generator,
            "binary_sha256": sha256_file(
                staging / "frozen/wirehair_v2_bench"),
            "binary_size": (staging / "frozen/wirehair_v2_bench").stat().st_size,
            "c_compiler": str(c_compiler),
            "c_compiler_sha256": sha256_file(c_compiler),
            "cxx_compiler": str(cxx_compiler),
            "cxx_compiler_sha256": sha256_file(cxx_compiler),
            "cxx_compiler_version": command_capture((str(cxx_compiler), "--version")),
        }
    finally:
        remove = subprocess.run(
            (str(tools["git"]), "-C", str(repo), "worktree", "remove", "--force",
             str(source)), stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            check=False, env=environment)
        if remove.returncode != 0 and source.exists():
            raise TimingError("could not remove detached build worktree")


def _boundary_smoke(staging: Path, design: Mapping[str, object],
                    environment: Mapping[str, str]) -> Dict[str, object]:
    tasks = generate_tasks()
    coordinates = ((3199, 4096), (3200, 4096),
                   (9999, 1280), (10000, 1280))
    records = []
    for K, bb in coordinates:
        task = next(value for value in tasks if value["K"] == K and
                    value["bb"] == bb and value["schedule"] == "burst" and
                    value["cache_state"] == "warm")
        command = command_for(design, task, cycle_index=2, evict_bytes=4096)
        rc, stdout, stderr = run_bounded(command, environment, 600)
        if rc != 0 or stderr:
            raise TimingError("boundary prelaunch smoke failed")
        parsed = parse_grouped_output(
            stdout, task, int(design["core"]), replacement_cycle=2)
        name = "boundary.K%d.bb%d.csv" % (K, bb)
        write_new(staging / "provenance" / name, stdout)
        records.append({
            "K": K, "bb": bb, "argv": command, "stdout_name": name,
            "stdout_sha256": sha256_bytes(stdout),
            "trace_sha256": parsed.trace_sha256,
            "control_route": parsed.preamble["control_preflight_rhs_route"],
            "candidate_route": parsed.preamble["candidate_preflight_rhs_route"],
            "expected_candidate_route": expected_rhs_route("p32_r7", K, bb),
            "common_success": parsed.common_success,
            "nonpromotional_contaminations": list(parsed.contaminations),
        })
    if not all(record["common_success"] for record in records):
        raise TimingError("route boundary smoke was not common-success")
    return {
        "timing_evidence": False,
        "scope": "bounded architecture and automatic-route prelaunch smoke under load",
        "records": records,
    }


def prepare_campaign(args: argparse.Namespace) -> None:
    result = Path(args.result_dir).resolve()
    repo = Path(args.repo).resolve()
    if result.exists() or result.is_symlink():
        raise TimingError("result directory already exists")
    if args.build_jobs <= 0 or args.core < 0 or args.controller_core < 0 or \
            args.thermal_core < 0 or args.numa_node < 0 or args.evict_bytes < 4096:
        raise TimingError("prepare integer argument outside domain")
    tools = {name: resolve_tool(name) for name in
             ("git", "cmake", "taskset", "numactl", "sudo", "fuser",
              "python3", "env", "true")}
    ninja = shutil.which("ninja")
    if ninja is not None:
        tools["ninja"] = Path(ninja).resolve()
    sudo_probe = subprocess.run(
        (str(tools["sudo"]), "-n", str(tools["true"])),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
        env=sanitized_environment(Path("/tmp"), allocator=False))
    if sudo_probe.returncode != 0 or sudo_probe.stdout or sudo_probe.stderr:
        raise TimingError("passwordless sudo preflight failed")
    c_compiler = Path(args.c_compiler).resolve() if args.c_compiler else resolve_tool("cc")
    cxx_compiler = Path(args.cxx_compiler).resolve() if args.cxx_compiler else resolve_tool("c++")
    head = _git_value(tools["git"], repo, "rev-parse", "HEAD^{commit}")
    tree = _git_value(tools["git"], repo, "rev-parse", "HEAD^{tree}")
    if _git_value(tools["git"], repo, "status", "--porcelain",
                  "--untracked-files=no"):
        raise TimingError("tracked source tree must be clean before freezing")
    harness_source = Path(__file__).resolve()
    sampler_source = harness_source.with_name("wirehair_expo_thermal_sampler.py")
    for source_path in (harness_source, sampler_source):
        tracked = subprocess.run(
            (str(tools["git"]), "-C", str(repo), "ls-files", "--error-unmatch",
             str(source_path.relative_to(repo))), stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, check=False,
            env=sanitized_environment(Path("/tmp"), allocator=False))
        if tracked.returncode != 0:
            raise TimingError("campaign source is not committed: %s" % source_path)
    timing_topology = topology_record(args.core, args.numa_node)
    controller_topology = topology_record(args.controller_core, args.numa_node)
    thermal_topology = topology_record(args.thermal_core, args.numa_node)
    if args.controller_core in timing_topology["llc_shared_cpus"] or \
            args.thermal_core in timing_topology["llc_shared_cpus"] or \
            args.controller_core == args.thermal_core or \
            args.thermal_core in controller_topology["siblings"] or \
            args.controller_core in thermal_topology["siblings"]:
        raise TimingError("controller/thermal CPU isolation is invalid")
    staging = result.with_name(result.name + ".prepare.%d" % os.getpid())
    workspace = result.with_name(result.name + ".build.%d" % os.getpid())
    if staging.exists() or workspace.exists():
        raise TimingError("stale prepare workspace exists")
    staging.mkdir(parents=True)
    workspace.mkdir(parents=True)
    try:
        for directory in ("frozen", "provenance", "external", "ledger",
                          ".transactions", "interrupted", "segments"):
            (staging / directory).mkdir()
        _copy_frozen(harness_source,
                     staging / "frozen/wh2_p32_dispatch_timing.py", 0o555)
        _copy_frozen(sampler_source,
                     staging / "frozen/wirehair_expo_thermal_sampler.py", 0o555)
        build = _build_frozen_binary(
            repo, head, workspace, staging, tools, args.build_jobs,
            c_compiler, cxx_compiler)
        external_sources = (
            (SUPERSEDED_PRELAUNCH, SUPERSEDED_PRELAUNCH_SHA256,
             "superseded_prelaunch_receipts.json"),
            (RECOVERY_ROOT / "validated_summary.json", RECOVERY_SUMMARY_SHA256,
             "recovery_validated_summary.json"),
            (RECOVERY_ROOT / "data_manifest.sha256", RECOVERY_DATA_MANIFEST_SHA256,
             "recovery_data_manifest.sha256"),
        )
        external_receipts = {}
        for source_path, expected_hash, destination_name in external_sources:
            if sha256_file(source_path) != expected_hash:
                raise TimingError("external receipt hash mismatch: %s" % source_path)
            _copy_frozen(source_path, staging / "external" / destination_name, 0o444)
            external_receipts[destination_name] = {
                "source": str(source_path), "sha256": expected_hash,
            }
        tasks = generate_tasks()
        manifest = b"".join(canonical_json(task) for task in tasks)
        write_new(staging / "tasks_manifest.jsonl", manifest)
        write_new(staging / "controller.lock", b"")
        tool_records = {name: {"path": str(path), "sha256": sha256_file(path)}
                        for name, path in sorted(tools.items())}
        provisional_design: Dict[str, object] = {
            "root": str(staging), "head": head, "tree": tree,
            "tools": tool_records, "core": args.core,
            "controller_core": args.controller_core,
            "thermal_core": args.thermal_core, "numa_node": args.numa_node,
            "evict_bytes": args.evict_bytes,
        }
        smoke_environment = sanitized_environment(workspace / "home", allocator=True)
        prelaunch = sealed_record(
            "wirehair.wh2.p32_dispatch.prelaunch.v2", {
                "created_utc": utc_now(),
                "boundary_smoke": _boundary_smoke(
                    staging, provisional_design, smoke_environment),
                "external_receipts": external_receipts,
                "superseded_artifact_sha256": SUPERSEDED_PRELAUNCH_SHA256,
                "superseded_artifact_launched": False,
            })
        write_new(staging / "prelaunch_receipt.json", canonical_json(prelaunch))
        immutable_files: Dict[str, str] = {}
        for directory_name in ("frozen", "provenance", "external"):
            for path in sorted((staging / directory_name).iterdir()):
                if path.is_file():
                    immutable_files[str(path.relative_to(staging))] = sha256_file(path)
        design = sealed_record(
            "wirehair.wh2.p32_dispatch.design.v2", {
                "root": str(result), "head": head, "tree": tree,
                "timing_scope": "full-payload decoder precode solve",
                "raw_architecture": True, "seed_fixes": "none",
                "arms": ARMS, "K": list(KS), "bb": list(WIDTHS),
                "schedules": list(SCHEDULES), "seeds": list(SEEDS),
                "cache_states": list(CACHE_STATES), "overhead": OVERHEAD,
                "loss": LOSS, "order": ORDER, "task_count": len(tasks),
                "rows_per_task": 32, "timed_rows_per_task": 24,
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
                "timing_topology": timing_topology,
                "controller_topology": controller_topology,
                "thermal_topology": thermal_topology,
                "thermal_limits": {"cpu_c": MAX_CPU_TEMP_C,
                    "dimm_c": MAX_DIMM_TEMP_C, "max_gap_s": MAX_THERMAL_GAP_S,
                    "max_margin_s": MAX_THERMAL_MARGIN_S},
                "tasks_manifest_sha256": sha256_bytes(manifest),
                "controller_lock_sha256": sha256_bytes(b""),
                "immutable_files": immutable_files, "tools": tool_records,
                "build": build,
                "prelaunch_receipt_sha256": sha256_file(
                    staging / "prelaunch_receipt.json"),
                "external_receipts": external_receipts,
                "supersedes": {
                    "artifact": "/tmp/wh2-p32-dispatch-postopt-timing-a3f5c66-v1",
                    "prelaunch_sha256": SUPERSEDED_PRELAUNCH_SHA256,
                    "reason": "wrong prod244 geometry and nonterminal thermal PID-file dependency",
                    "launched": False,
                },
                "transaction_policy": "one atomic immutable directory per task; bounded resumable segments",
                "thermal_policy": "controller-owned sole reader through post-end sample and graceful stop",
            })
        write_new(staging / "design.json", canonical_json(design))
        prepare = sealed_record(
            "wirehair.wh2.p32_dispatch.prepare.v2", {
                "prepared_utc": utc_now(), "head": head, "tree": tree,
                "design_sha256": sha256_file(staging / "design.json"),
                "tasks_manifest_sha256": sha256_bytes(manifest),
                "prelaunch_receipt_sha256": sha256_file(
                    staging / "prelaunch_receipt.json"),
                "binary_sha256": build["binary_sha256"],
                "runner_sha256": immutable_files[
                    "frozen/wh2_p32_dispatch_timing.py"],
                "thermal_sampler_sha256": immutable_files[
                    "frozen/wirehair_expo_thermal_sampler.py"],
            })
        write_new(staging / "prepare_receipt.json", canonical_json(prepare))
        for name in ("frozen", "provenance", "external"):
            os.chmod(staging / name, 0o555)
        os.replace(staging, result)
        print(json.dumps({
            "result_dir": str(result), "head": head, "tree": tree,
            "task_count": len(tasks),
            "design_sha256": sha256_file(result / "design.json"),
            "manifest_sha256": sha256_file(result / "tasks_manifest.jsonl"),
            "binary_sha256": build["binary_sha256"],
            "prelaunch_receipt_sha256": sha256_file(
                result / "prelaunch_receipt.json"),
        }, sort_keys=True))
    finally:
        if staging.exists():
            shutil.rmtree(staging)
        if workspace.exists():
            shutil.rmtree(workspace)


def load_design(root: Path) -> Dict[str, object]:
    design = verify_sealed(load_canonical(root / "design.json", "design"),
                           "wirehair.wh2.p32_dispatch.design.v2", "design")
    if design.get("root") != str(root.resolve()):
        raise TimingError("campaign root moved after preparation")
    expected = {
        "timing_scope": "full-payload decoder precode solve",
        "raw_architecture": True, "seed_fixes": "none", "arms": ARMS,
        "K": list(KS), "bb": list(WIDTHS), "schedules": list(SCHEDULES),
        "seeds": list(SEEDS), "cache_states": list(CACHE_STATES),
        "overhead": OVERHEAD, "loss": LOSS, "order": ORDER,
        "task_count": 1404, "rows_per_task": 32, "timed_rows_per_task": 24,
        "controller_lock_sha256": sha256_bytes(b""),
        "max_environmental_attempts": MAX_ENVIRONMENTAL_ATTEMPTS,
        "minor_fault_max": MAX_MINOR_FAULTS,
        "output_bounds": {"stdout": MAX_STDOUT_BYTES,
                          "stderr": MAX_STDERR_BYTES},
        "allocator_environment": {
            "MALLOC_MMAP_THRESHOLD_": MALLOC_MMAP_THRESHOLD,
            "MALLOC_TRIM_THRESHOLD_": MALLOC_TRIM_THRESHOLD,
        },
        "transaction_policy":
            "one atomic immutable directory per task; bounded resumable segments",
        "thermal_policy":
            "controller-owned sole reader through post-end sample and graceful stop",
        "thermal_limits": {"cpu_c": MAX_CPU_TEMP_C,
                           "dimm_c": MAX_DIMM_TEMP_C,
                           "max_gap_s": MAX_THERMAL_GAP_S,
                           "max_margin_s": MAX_THERMAL_MARGIN_S},
    }
    for key, value in expected.items():
        if design.get(key) != value:
            raise TimingError("frozen design policy changed: " + key)
    supersedes = design.get("supersedes")
    expected_supersedes = {
        "artifact": "/tmp/wh2-p32-dispatch-postopt-timing-a3f5c66-v1",
        "prelaunch_sha256": SUPERSEDED_PRELAUNCH_SHA256,
        "reason":
            "wrong prod244 geometry and nonterminal thermal PID-file dependency",
        "launched": False,
    }
    if supersedes != expected_supersedes or \
            design.get("external_receipts") != expected_external_receipts():
        raise TimingError("superseded artifact receipt changed")
    for core_key, topology_key in (
            ("core", "timing_topology"),
            ("controller_core", "controller_topology"),
            ("thermal_core", "thermal_topology")):
        core = design.get(core_key)
        topology = design.get(topology_key)
        if type(core) is not int or core < 0 or not isinstance(topology, dict) or \
                topology.get("core") != core or \
                topology.get("numa_node") != design.get("numa_node"):
            raise TimingError("frozen topology receipt malformed: " + core_key)
    if type(design.get("numa_node")) is not int or design["numa_node"] < 0 or \
            type(design.get("evict_bytes")) is not int or \
            design["evict_bytes"] < 4096:
        raise TimingError("frozen NUMA/eviction policy malformed")
    required_tools = {"git", "cmake", "taskset", "numactl", "sudo", "fuser",
                      "python3", "env", "true"}
    tool_names = set(design.get("tools", {}))
    if not required_tools <= tool_names or not tool_names <= required_tools | {"ninja"}:
        raise TimingError("frozen tool ledger changed")
    return design


def load_tasks(root: Path, design: Mapping[str, object]) -> List[Dict[str, object]]:
    path = root / "tasks_manifest.jsonl"
    if sha256_file(path) != design.get("tasks_manifest_sha256"):
        raise TimingError("task manifest hash mismatch")
    tasks: List[Dict[str, object]] = []
    for line_number, line in enumerate(path.read_bytes().splitlines(keepends=True), 1):
        try:
            task = json.loads(line.decode("ascii"))
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            raise TimingError("task manifest JSON malformed") from exc
        if not isinstance(task, dict) or canonical_json(task) != line or \
                task.get("job") != len(tasks):
            raise TimingError("task manifest line %d is noncanonical" % line_number)
        tasks.append(task)
    if tasks != generate_tasks():
        raise TimingError("task manifest Cartesian product changed")
    return tasks


def verify_immutable_inputs(root: Path, design: Mapping[str, object]) -> None:
    lock_path = root / "controller.lock"
    if lock_path.is_symlink() or not lock_path.is_file() or \
            sha256_file(lock_path) != design.get("controller_lock_sha256"):
        raise TimingError("controller lock identity changed")
    immutable = design.get("immutable_files")
    if not isinstance(immutable, dict) or not immutable:
        raise TimingError("immutable input ledger missing")
    for relative, expected_hash in immutable.items():
        if not isinstance(relative, str) or not isinstance(expected_hash, str) or \
                SHA256_RE.fullmatch(expected_hash) is None:
            raise TimingError("immutable input ledger malformed")
        path = root / relative
        if path.is_symlink() or not path.is_file() or sha256_file(path) != expected_hash:
            raise TimingError("immutable input changed: " + relative)
    for name, receipt in expected_external_receipts().items():
        if immutable.get("external/" + name) != receipt["sha256"]:
            raise TimingError("external receipt is not immutably bound: " + name)
    tools = design.get("tools")
    if not isinstance(tools, dict):
        raise TimingError("tool receipt ledger missing")
    for name, receipt in tools.items():
        if not isinstance(receipt, dict):
            raise TimingError("tool receipt malformed")
        path = Path(str(receipt.get("path")))
        if not path.is_file() or sha256_file(path) != receipt.get("sha256"):
            raise TimingError("runtime tool changed: " + str(name))
    prelaunch = verify_sealed(
        load_canonical(root / "prelaunch_receipt.json", "prelaunch receipt"),
        "wirehair.wh2.p32_dispatch.prelaunch.v2", "prelaunch receipt")
    if sha256_file(root / "prelaunch_receipt.json") != \
            design.get("prelaunch_receipt_sha256") or \
            prelaunch.get("superseded_artifact_launched") is not False or \
            prelaunch.get("superseded_artifact_sha256") != \
                SUPERSEDED_PRELAUNCH_SHA256 or \
            prelaunch.get("external_receipts") != expected_external_receipts():
        raise TimingError("prelaunch binding changed")
    validate_utc_timestamp(prelaunch.get("created_utc"), "prelaunch creation")
    smoke = prelaunch.get("boundary_smoke")
    if not isinstance(smoke, dict) or set(smoke) != {
            "timing_evidence", "scope", "records"} or \
            smoke.get("timing_evidence") is not False or \
            smoke.get("scope") != \
                "bounded architecture and automatic-route prelaunch smoke under load" or \
            not isinstance(smoke.get("records"), list):
        raise TimingError("boundary prelaunch smoke receipt malformed")
    tasks = generate_tasks()
    expected_coordinates = ((3199, 4096), (3200, 4096),
                            (9999, 1280), (10000, 1280))
    if [(record.get("K"), record.get("bb"))
            for record in smoke["records"] if isinstance(record, dict)] != \
            list(expected_coordinates):
        raise TimingError("boundary prelaunch coordinate set changed")
    for record, (K, bb) in zip(smoke["records"], expected_coordinates):
        if not isinstance(record, dict) or set(record) != {
                "K", "bb", "argv", "stdout_name", "stdout_sha256",
                "trace_sha256", "control_route", "candidate_route",
                "expected_candidate_route", "common_success",
                "nonpromotional_contaminations"}:
            raise TimingError("boundary prelaunch record field set changed")
        task = next(value for value in tasks if value["K"] == K and
                    value["bb"] == bb and value["schedule"] == "burst" and
                    value["cache_state"] == "warm")
        stdout_name = "boundary.K%d.bb%d.csv" % (K, bb)
        stdout_path = root / "provenance" / stdout_name
        raw = stdout_path.read_bytes()
        parsed = parse_grouped_output(
            raw, task, int(design["core"]), replacement_cycle=2)
        expected_argv = command_for(design, task, cycle_index=2, evict_bytes=4096)
        actual_argv = record["argv"]
        if not isinstance(actual_argv, list) or len(actual_argv) != len(expected_argv):
            raise TimingError("boundary prelaunch argv malformed")
        normalized_argv = list(actual_argv)
        if Path(str(normalized_argv[6])).name != "wirehair_v2_bench":
            raise TimingError("boundary prelaunch binary identity malformed")
        normalized_argv[6] = expected_argv[6]
        exact_record = {
            "K": K, "bb": bb, "argv": actual_argv,
            "stdout_name": stdout_name, "stdout_sha256": sha256_bytes(raw),
            "trace_sha256": parsed.trace_sha256,
            "control_route": parsed.preamble["control_preflight_rhs_route"],
            "candidate_route": parsed.preamble["candidate_preflight_rhs_route"],
            "expected_candidate_route": expected_rhs_route("p32_r7", K, bb),
            "common_success": parsed.common_success,
            "nonpromotional_contaminations": list(parsed.contaminations),
        }
        if normalized_argv != expected_argv or record != exact_record or \
                record["common_success"] is not True:
            raise TimingError("boundary prelaunch receipt does not replay")
    prepare = verify_sealed(
        load_canonical(root / "prepare_receipt.json", "prepare receipt"),
        "wirehair.wh2.p32_dispatch.prepare.v2", "prepare receipt")
    if prepare.get("design_sha256") != sha256_file(root / "design.json") or \
            prepare.get("head") != design.get("head") or \
            prepare.get("tree") != design.get("tree") or \
            prepare.get("tasks_manifest_sha256") != sha256_file(
                root / "tasks_manifest.jsonl") or \
            prepare.get("prelaunch_receipt_sha256") != sha256_file(
                root / "prelaunch_receipt.json") or \
            prepare.get("binary_sha256") != immutable.get(
                "frozen/wirehair_v2_bench") or \
            prepare.get("runner_sha256") != immutable.get(
                "frozen/wh2_p32_dispatch_timing.py") or \
            prepare.get("thermal_sampler_sha256") != immutable.get(
                "frozen/wirehair_expo_thermal_sampler.py"):
        raise TimingError("prepare receipt does not bind frozen inputs")
    validate_utc_timestamp(prepare.get("prepared_utc"), "prepare completion")


def verify_root_layout(root: Path) -> None:
    required_directories = {
        ".transactions", "external", "frozen", "interrupted", "ledger",
        "provenance", "segments",
    }
    required_files = {
        "controller.lock", "design.json", "prelaunch_receipt.json",
        "prepare_receipt.json", "tasks_manifest.jsonl",
    }
    optional = {"campaign_receipt.json", "runtime-home",
                "validated_summary.json"}
    actual = {path.name for path in root.iterdir()}
    if not required_directories | required_files <= actual or \
            not actual <= required_directories | required_files | optional:
        raise TimingError("campaign root layout contains missing or unknown entries")
    for name in required_directories:
        path = root / name
        if path.is_symlink() or not path.is_dir():
            raise TimingError("campaign directory is missing or unsafe: " + name)
    for name in required_files:
        path = root / name
        if path.is_symlink() or not path.is_file():
            raise TimingError("campaign file is missing or unsafe: " + name)
    for name in optional & actual:
        path = root / name
        expected_kind = path.is_dir() if name == "runtime-home" else path.is_file()
        if path.is_symlink() or not expected_kind:
            raise TimingError("optional campaign entry is unsafe: " + name)
        if name == "runtime-home" and any(path.iterdir()):
            raise TimingError("runtime HOME contains unreceipted state")


def _task_receipt_summary(parsed: ParsedOutput) -> Dict[str, object]:
    return {
        "stdout_sha256": parsed.stdout_sha256,
        "trace_sha256": parsed.trace_sha256, "cell_class": parsed.cell_class,
        "common_success": parsed.common_success,
        "timed_control_ns": parsed.timed_control_ns,
        "timed_candidate_ns": parsed.timed_candidate_ns,
        "control_work_sha256": sha256_bytes(canonical_json(
            list(parsed.work_signatures["control"]))),
        "candidate_work_sha256": sha256_bytes(canonical_json(
            list(parsed.work_signatures["candidate"]))),
        "control_rhs_route": parsed.preamble["control_preflight_rhs_route"],
        "candidate_rhs_route": parsed.preamble["candidate_preflight_rhs_route"],
    }


def validate_task_entry(root: Path, design: Mapping[str, object],
                        task: Mapping[str, object]) -> Tuple[Dict[str, object], ParsedOutput]:
    directory = root / "ledger" / ("%04d" % int(task["job"]))
    if directory.is_symlink() or not directory.is_dir():
        raise TimingError("task ledger directory missing or unsafe")
    names = sorted(path.name for path in directory.iterdir())
    if names != ["intent.json", "receipt.json", "stderr.bin", "stdout.csv"]:
        raise TimingError("task ledger file set changed")
    for path in directory.iterdir():
        if path.is_symlink() or not path.is_file():
            raise TimingError("task ledger contains unsafe file")
    if (directory / "stderr.bin").read_bytes() != b"":
        raise TimingError("successful task has nonempty stderr")
    raw = (directory / "stdout.csv").read_bytes()
    parsed = parse_grouped_output(raw, task, int(design["core"]))
    intent = verify_sealed(
        load_canonical(directory / "intent.json", "task intent"),
        "wirehair.wh2.p32_dispatch.intent.v2", "task intent")
    if set(intent) != {"schema", "self_sha256_excluding_field", "created_utc",
                       "job", "task_id", "task_sha256"}:
        raise TimingError("task intent field set changed")
    validate_utc_timestamp(intent.get("created_utc"), "task intent creation")
    if intent.get("job") != task["job"] or intent.get("task_id") != task["task_id"] or \
            intent.get("task_sha256") != sha256_bytes(canonical_json(task)):
        raise TimingError("task intent identity mismatch")
    receipt = verify_sealed(
        load_canonical(directory / "receipt.json", "task receipt"),
        "wirehair.wh2.p32_dispatch.task.v2", "task receipt")
    receipt_fields = {
        "schema", "self_sha256_excluding_field", "job", "task_id",
        "task_sha256", "argv", "child_environment", "attempt",
        "prior_contaminations", "started_utc", "start_monotonic_ns",
        "end_monotonic_ns", "duration_ns", "stderr_sha256", "stdout_sha256",
        "trace_sha256", "cell_class", "common_success", "timed_control_ns",
        "timed_candidate_ns", "control_work_sha256", "candidate_work_sha256",
        "control_rhs_route", "candidate_rhs_route",
    }
    if set(receipt) != receipt_fields:
        raise TimingError("task receipt field set changed")
    expected = {
        "job": task["job"], "task_id": task["task_id"],
        "task_sha256": sha256_bytes(canonical_json(task)),
        "stderr_sha256": sha256_bytes(b""), **_task_receipt_summary(parsed),
    }
    for key, value in expected.items():
        if receipt.get(key) != value:
            raise TimingError("task receipt mismatch: " + key)
    command = command_for(design, task)
    if receipt.get("argv") != command or receipt.get("child_environment") != \
            sanitized_environment(root / "runtime-home", allocator=True):
        raise TimingError("task command/environment receipt mismatch")
    attempt = receipt.get("attempt")
    prior = receipt.get("prior_contaminations")
    if type(attempt) is not int or not 0 <= attempt < MAX_ENVIRONMENTAL_ATTEMPTS or \
            not isinstance(prior, list) or len(prior) != attempt:
        raise TimingError("task retry ledger malformed")
    started_ns = receipt.get("start_monotonic_ns")
    ended_ns = receipt.get("end_monotonic_ns")
    duration_ns = receipt.get("duration_ns")
    if type(started_ns) is not int or type(ended_ns) is not int or \
            type(duration_ns) is not int or started_ns < 0 or \
            ended_ns <= started_ns or duration_ns != ended_ns - started_ns:
        raise TimingError("task timing receipt malformed")
    validate_utc_timestamp(receipt.get("started_utc"), "task start")
    if parsed.contaminations:
        raise TimingError("terminal task still contains timing contamination")
    previous_end = -1
    for index, record in enumerate(prior):
        if not isinstance(record, dict) or record.get("attempt") != index or \
                set(record) != {"attempt", "started_utc", "start_monotonic_ns",
                                "end_monotonic_ns", "stdout_sha256",
                                "stderr_sha256", "contaminations"}:
            raise TimingError("prior contamination receipt malformed")
        prior_start = record["start_monotonic_ns"]
        prior_end = record["end_monotonic_ns"]
        contaminations = record["contaminations"]
        if type(prior_start) is not int or type(prior_end) is not int or \
                prior_start < 0 or prior_end <= prior_start or \
                prior_start < previous_end or prior_end > started_ns or \
                SHA256_RE.fullmatch(str(record["stdout_sha256"])) is None or \
                record["stderr_sha256"] != sha256_bytes(b"") or \
                not isinstance(contaminations, list) or not contaminations or \
                contaminations != sorted(set(contaminations)) or \
                any(not isinstance(value, str) or not value
                    for value in contaminations):
            raise TimingError("prior contamination evidence malformed")
        validate_utc_timestamp(record["started_utc"],
                               "prior contamination start")
        previous_end = prior_end
    return receipt, parsed


def existing_task_entries(root: Path, design: Mapping[str, object],
                          tasks: Sequence[Mapping[str, object]]) -> Dict[int, Dict[str, object]]:
    entries = sorted((root / "ledger").iterdir())
    if len(entries) > len(tasks):
        raise TimingError("task ledger exceeds manifest")
    result = {}
    for entry in entries:
        if entry.is_symlink() or not entry.is_dir() or TASK_DIR_RE.fullmatch(entry.name) is None:
            raise TimingError("task ledger contains unknown entry")
        job = int(entry.name)
        if job >= len(tasks):
            raise TimingError("task ledger job outside manifest")
        receipt, parsed = validate_task_entry(root, design, tasks[job])
        result[job] = {"receipt": receipt, "parsed": parsed,
                       "receipt_sha256": sha256_file(entry / "receipt.json")}
    return result


def filler_pids() -> Tuple[int, ...]:
    signatures = (b"wirehair_load_fillers.sh", b"while :; do :; done")
    result = []
    for path in Path("/proc").glob("[0-9]*/cmdline"):
        try:
            command = path.read_bytes().replace(b"\0", b" ")
        except OSError:
            continue
        if any(signature in command for signature in signatures):
            result.append(int(path.parent.name))
    return tuple(sorted(result))


def verify_runtime_topology(design: Mapping[str, object]) -> None:
    for core_key, frozen_key in (
            ("core", "timing_topology"),
            ("controller_core", "controller_topology"),
            ("thermal_core", "thermal_topology")):
        current = topology_record(int(design[core_key]), int(design["numa_node"]))
        if current != design.get(frozen_key):
            raise TimingError("runtime CPU topology/power policy changed: " + core_key)


def _public_sampler_identity(identity: Mapping[str, object]) -> Dict[str, object]:
    result = dict(identity)
    argv = list(result.get("cmdline", []))
    for index, value in enumerate(argv[:-1]):
        if value == "--pid-file":
            argv[index + 1] = "<ephemeral-bootstrap-only>"
    result["cmdline"] = argv
    return result


def _sampler_cmdline(root: Path, design: Mapping[str, object], segment: int,
                     *, public: bool) -> List[str]:
    tools = design["tools"]
    pid_file = root / "segments" / (".segment%04d.bootstrap.pid" % segment)
    return [
        str(tools["python3"]["path"]),
        str(root / "frozen/wirehair_expo_thermal_sampler.py"),
        "--csv", str(root / "segments" /
            (".segment%04d.thermal.csv.part" % segment)),
        "--pid-file", "<ephemeral-bootstrap-only>" if public else str(pid_file),
        "--interval", "1.0", "--dimm-attempts", "5",
        "--dimm-retry-delay", "0.01",
    ]


@dataclass
class SamplerOwner:
    process: subprocess.Popen[bytes]
    pid: int
    identity: Dict[str, object]
    csv_part: Path
    pid_file: Path
    segment: int


def start_sampler(root: Path, design: Mapping[str, object], segment: int) -> SamplerOwner:
    tools = design["tools"]
    sudo_path = Path(str(tools["sudo"]["path"]))
    fuser_path = Path(str(tools["fuser"]["path"]))
    if sole_i2c_readers(fuser_path, sudo_path):
        raise TimingError("refusing to signal or coexist with an existing I2C reader")
    csv_part = root / "segments" / (".segment%04d.thermal.csv.part" % segment)
    pid_file = root / "segments" / (".segment%04d.bootstrap.pid" % segment)
    for path in (csv_part, pid_file):
        if path.exists() or path.is_symlink():
            raise TimingError("stale sampler bootstrap artifact exists")
    environment = sanitized_environment(root / "runtime-home", allocator=False)
    sampler_argv = _sampler_cmdline(root, design, segment, public=False)
    command = [
        str(sudo_path), "-n", str(tools["env"]["path"]), "-i",
        "HOME=" + environment["HOME"], "PATH=" + environment["PATH"],
        "LC_ALL=C", "LANG=C", "TZ=UTC", "PYTHONDONTWRITEBYTECODE=1",
        str(tools["taskset"]["path"]), "-c", str(design["thermal_core"]),
        *sampler_argv,
    ]
    process = subprocess.Popen(
        command, env=environment, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        start_new_session=True, close_fds=True)
    deadline = time.monotonic() + 20.0
    last_error: Optional[BaseException] = None
    while time.monotonic() < deadline:
        if process.poll() is not None:
            stdout, stderr = process.communicate(timeout=1)
            raise TimingError("thermal sampler exited during bootstrap: %r %r" %
                              (stdout[:1000], stderr[:1000]))
        try:
            text = pid_file.read_text(encoding="ascii")
            if re.fullmatch(r"[1-9][0-9]*\n", text) is None:
                raise TimingError("sampler bootstrap PID is noncanonical")
            pid = int(text)
            _parse_thermal_csv(stable_file_bytes(csv_part, attempts=3))
            identity = capture_process_identity(
                pid, int(design["thermal_core"]), csv_part)
            if sole_i2c_readers(fuser_path, sudo_path) != (pid,):
                raise TimingError("sampler is not the sole I2C reader")
            if identity["cmdline"] != sampler_argv:
                raise TimingError("sampler cmdline receipt is not exact")
            return SamplerOwner(process, pid, identity, csv_part, pid_file, segment)
        except (OSError, TimingError, ValueError) as exc:
            last_error = exc
            time.sleep(0.05)
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    process.wait(timeout=10)
    raise TimingError("thermal sampler bootstrap timed out: %s" % last_error)


def stop_sampler(owner: SamplerOwner, root: Path, design: Mapping[str, object],
                 timing_start_s: float, benchmark_end_s: float) -> Dict[str, object]:
    tools = design["tools"]
    sudo_path = Path(str(tools["sudo"]["path"]))
    fuser_path = Path(str(tools["fuser"]["path"]))
    deadline = time.monotonic() + 15.0
    while True:
        if not process_identity_matches(
                owner.identity, int(design["thermal_core"]), owner.csv_part):
            raise TimingError("thermal sampler identity changed before graceful stop")
        if latest_thermal_time(owner.csv_part) >= benchmark_end_s:
            break
        if time.monotonic() >= deadline:
            raise TimingError("thermal sampler did not emit a post-end sample")
        time.sleep(0.05)
    identity_end = capture_process_identity(
        owner.pid, int(design["thermal_core"]), owner.csv_part)
    if identity_end != owner.identity or \
            sole_i2c_readers(fuser_path, sudo_path) != (owner.pid,):
        raise TimingError("thermal ownership changed before shutdown")
    environment = sanitized_environment(root / "runtime-home", allocator=False)
    stop_command = (
        str(sudo_path), "-n", str(tools["env"]["path"]), "-i",
        "HOME=" + environment["HOME"], "PATH=" + environment["PATH"],
        "LC_ALL=C", "LANG=C", "TZ=UTC", "PYTHONDONTWRITEBYTECODE=1",
        str(tools["python3"]["path"]), "-c", PIDFD_STOP_PROGRAM,
        str(owner.pid), str(owner.identity["start_tick"]),
        str(owner.identity["cmdline_sha256"]),
    )
    stop_rc, stop_stdout, stop_stderr = run_bounded(
        stop_command, environment, 10.0, 1024, 1024)
    if stop_rc != 0 or stop_stdout or stop_stderr:
        raise TimingError("graceful sampler stop command failed")
    deadline = time.monotonic() + 15.0
    while process_identity_matches(
            owner.identity, int(design["thermal_core"]), owner.csv_part):
        if time.monotonic() >= deadline:
            raise TimingError("thermal sampler ignored graceful stop")
        time.sleep(0.05)
    stdout, stderr = owner.process.communicate(timeout=15)
    if owner.process.returncode != 0 or stdout or stderr:
        raise TimingError("thermal sampler shutdown was not clean")
    if owner.pid_file.exists() or owner.pid_file.is_symlink():
        raise TimingError("sampler did not remove ephemeral bootstrap PID file")
    if sole_i2c_readers(fuser_path, sudo_path):
        raise TimingError("I2C reader remained after sampler shutdown")
    raw = stable_file_bytes(owner.csv_part)
    summary = validate_thermal_interval(raw, timing_start_s, benchmark_end_s)
    stat = owner.csv_part.stat()
    final = root / "segments" / ("segment%04d.thermal.csv" % owner.segment)
    if final.exists():
        raise TimingError("terminal thermal CSV already exists")
    os.chmod(owner.csv_part, 0o444)
    os.replace(owner.csv_part, final)
    fsync_directory(final.parent)
    return {
        "sampler_identity_start": _public_sampler_identity(owner.identity),
        "sampler_identity_end": _public_sampler_identity(identity_end),
        "graceful_stop": True,
        "graceful_stop_mechanism":
            "sudo-python-pidfd-identity-verified-sigterm",
        "post_end_sample": True,
        "thermal_csv_name": final.name,
        "thermal_csv_sha256": sha256_bytes(raw), "thermal_csv_size": len(raw),
        "thermal_csv_device": stat.st_dev, "thermal_csv_inode": stat.st_ino,
        "thermal_summary": summary,
    }


def _next_segment(root: Path) -> int:
    terminals = sorted((root / "segments").glob("segment*.terminal.json"))
    for index, path in enumerate(terminals):
        if path.name != "segment%04d.terminal.json" % index:
            raise TimingError("segment terminal sequence has a gap")
    unknown = [path for path in (root / "segments").iterdir()
               if not re.fullmatch(r"segment[0-9]{4}\.(?:terminal\.json|thermal\.csv)",
                                   path.name)]
    if unknown:
        raise TimingError("segments directory contains stale/unknown entry: %s" %
                          unknown[0].name)
    return len(terminals)


def _task_execution(root: Path, design: Mapping[str, object],
                    task: Mapping[str, object], timeout_s: float) -> Dict[str, object]:
    environment = sanitized_environment(root / "runtime-home", allocator=True)
    command = command_for(design, task)
    staging = begin_task_transaction(root, task)
    prior: List[Dict[str, object]] = []
    for attempt in range(MAX_ENVIRONMENTAL_ATTEMPTS):
        started_utc = utc_now()
        start_ns = time.monotonic_ns()
        rc, stdout, stderr = run_bounded(command, environment, timeout_s)
        end_ns = time.monotonic_ns()
        if rc != 0 or stderr:
            raise TimingError("substantive task child failure exit=%d stderr=%r" %
                              (rc, stderr[:1000]))
        parsed = parse_grouped_output(stdout, task, int(design["core"]))
        if parsed.contaminations:
            prior.append({
                "attempt": attempt, "started_utc": started_utc,
                "start_monotonic_ns": start_ns, "end_monotonic_ns": end_ns,
                "stdout_sha256": sha256_bytes(stdout),
                "stderr_sha256": sha256_bytes(stderr),
                "contaminations": list(parsed.contaminations),
            })
            if attempt + 1 == MAX_ENVIRONMENTAL_ATTEMPTS:
                raise ContaminationError("environmental retry limit exhausted")
            continue
        receipt = sealed_record(
            "wirehair.wh2.p32_dispatch.task.v2", {
                "job": task["job"], "task_id": task["task_id"],
                "task_sha256": sha256_bytes(canonical_json(task)),
                "argv": command, "child_environment": environment,
                "attempt": attempt, "prior_contaminations": prior,
                "started_utc": started_utc, "start_monotonic_ns": start_ns,
                "end_monotonic_ns": end_ns, "duration_ns": end_ns - start_ns,
                "stderr_sha256": sha256_bytes(stderr),
                **_task_receipt_summary(parsed),
            })
        final = commit_task_transaction(
            root, task, stdout, stderr, receipt, staging=staging)
        return {"job": task["job"],
                "receipt_sha256": sha256_file(final / "receipt.json"),
                "retry_count": attempt}
    raise TimingError("unreachable task retry state")


def _write_campaign_receipt(root: Path, design: Mapping[str, object],
                            tasks: Sequence[Mapping[str, object]]) -> None:
    path = root / "campaign_receipt.json"
    task_receipts = []
    for task in tasks:
        receipt_path = root / "ledger" / ("%04d" % int(task["job"])) / "receipt.json"
        task_receipts.append({"job": task["job"],
                              "receipt_sha256": sha256_file(receipt_path)})
    terminals = sorted((root / "segments").glob("segment*.terminal.json"))
    terminal_receipts = [{"name": path.name, "sha256": sha256_file(path)}
                         for path in terminals]
    exact = {
        "design_sha256": sha256_file(root / "design.json"),
        "prepare_receipt_sha256": sha256_file(root / "prepare_receipt.json"),
        "prelaunch_receipt_sha256": sha256_file(
            root / "prelaunch_receipt.json"),
        "tasks_manifest_sha256": sha256_file(root / "tasks_manifest.jsonl"),
        "task_count": len(tasks), "task_receipts": task_receipts,
        "terminal_receipts": terminal_receipts,
    }
    if path.exists() or path.is_symlink():
        campaign = verify_sealed(
            load_canonical(path, "campaign receipt"),
            "wirehair.wh2.p32_dispatch.campaign.v2", "campaign receipt")
        if set(campaign) != {
                "schema", "self_sha256_excluding_field", "completed_utc",
                *exact} or any(
                    campaign.get(key) != value for key, value in exact.items()):
            raise TimingError("existing campaign receipt does not replay")
        validate_utc_timestamp(campaign.get("completed_utc"),
                               "campaign completion")
        return
    receipt = sealed_record(
        "wirehair.wh2.p32_dispatch.campaign.v2", {
            "completed_utc": utc_now(), **exact,
        })
    write_new(path, canonical_json(receipt))


def _run_campaign_locked(args: argparse.Namespace, root: Path,
                         design: Mapping[str, object]) -> None:
    verify_root_layout(root)
    tasks = load_tasks(root, design)
    verify_immutable_inputs(root, design)
    if args.max_tasks <= 0 or args.timeout_seconds <= 0:
        raise TimingError("run bounds must be positive")
    if filler_pids():
        raise TimingError("load filler supervisor/workers are still running")
    verify_runtime_topology(design)
    os.sched_setaffinity(0, {int(design["controller_core"])})
    if os.sched_getaffinity(0) != {int(design["controller_core"])}:
        raise TimingError("controller affinity could not be sealed")
    runtime_home = root / "runtime-home"
    runtime_home.mkdir(exist_ok=True)
    if runtime_home.is_symlink() or any(runtime_home.iterdir()):
        raise TimingError("runtime HOME must be an empty real directory")
    bindings = terminal_task_bindings(root)
    preexisting = existing_task_entries(root, design, tasks)
    verify_segments(root, design, preexisting, allow_unbound_archives=True)
    recovered = adopt_unbound_interrupted_archives(root)
    recovered.extend(recover_unbound_task_entries(root, bindings))
    recovered.extend(recover_interrupted_transactions(root))
    existing = existing_task_entries(root, design, tasks)
    if set(existing) != set(bindings):
        raise TimingError("resumed task ledger and terminal bindings disagree")
    remaining = [task for task in tasks if int(task["job"]) not in existing]
    if remaining and ((root / "campaign_receipt.json").exists() or
                      (root / "campaign_receipt.json").is_symlink()):
        raise TimingError("incomplete campaign has a terminal campaign receipt")
    if not remaining:
        if recovered:
            raise TimingError(
                "unbound recovery evidence has no remaining task to rerun")
        _write_campaign_receipt(root, design, tasks)
        print(json.dumps({"complete": True, "resumed": len(existing),
                          "completed": 0}, sort_keys=True))
        return
    selected = remaining[:args.max_tasks]
    segment = _next_segment(root)
    owner = start_sampler(root, design, segment)
    timing_start_s = time.monotonic()
    started_utc = utc_now()
    completed: List[Dict[str, object]] = []
    error: Optional[BaseException] = None
    try:
        if timing_start_s - latest_thermal_time(owner.csv_part) > MAX_THERMAL_MARGIN_S:
            raise TimingError("sampler baseline is stale")
        for task in selected:
            completed.append(_task_execution(
                root, design, task, args.timeout_seconds))
            if filler_pids():
                raise TimingError("load filler appeared during quiet timing")
            tools = design["tools"]
            if not process_identity_matches(
                    owner.identity, int(design["thermal_core"]),
                    owner.csv_part) or sole_i2c_readers(
                        Path(str(tools["fuser"]["path"])),
                        Path(str(tools["sudo"]["path"]))) != (owner.pid,):
                raise TimingError("thermal ownership changed during timing")
            if len(completed) % 10 == 0 or len(completed) == len(selected):
                print("segment=%d progress=%d/%d total=%d/%d" %
                      (segment, len(completed), len(selected),
                       len(existing) + len(completed), len(tasks)), flush=True)
    except BaseException as exc:
        error = exc
    benchmark_end_s = time.monotonic()
    terminal_payload: Dict[str, object] = {}
    stop_error: Optional[BaseException] = None
    try:
        terminal_payload = stop_sampler(
            owner, root, design, timing_start_s, benchmark_end_s)
    except BaseException as exc:
        stop_error = exc
    if stop_error is not None:
        if error is None:
            error = stop_error
        else:
            error = TimingError("task failure %s; sampler finalization failure %s" %
                                (error, stop_error))
    if terminal_payload:
        status = "failed" if error is not None else \
            "complete" if len(existing) + len(completed) == len(tasks) else \
            "bounded-incomplete"
        terminal = sealed_record(
            "wirehair.wh2.p32_dispatch.segment_terminal.v2", {
                "segment": segment, "status": status,
                "started_utc": started_utc, "ended_utc": utc_now(),
                "timing_start_monotonic_s": timing_start_s,
                "benchmark_end_monotonic_s": benchmark_end_s,
                "duration_s": benchmark_end_s - timing_start_s,
                "controller_pid": os.getpid(),
                "controller_affinity": sorted(os.sched_getaffinity(0)),
                "design_sha256": sha256_file(root / "design.json"),
                "prepare_receipt_sha256": sha256_file(
                    root / "prepare_receipt.json"),
                "prelaunch_receipt_sha256": sha256_file(
                    root / "prelaunch_receipt.json"),
                "tasks_manifest_sha256": sha256_file(
                    root / "tasks_manifest.jsonl"),
                "resumed_jobs_before": sorted(existing),
                "completed_tasks": completed,
                "recovered_interrupted_transactions": recovered,
                "failure": None if error is None else {
                    "class": type(error).__name__, "message": str(error)},
                **terminal_payload,
            })
        terminal_path = root / "segments" / (
            "segment%04d.terminal.json" % segment)
        write_new(terminal_path, canonical_json(terminal))
        sealed_entries = existing_task_entries(root, design, tasks)
        verify_segments(root, design, sealed_entries)
    if error is not None:
        raise error
    if len(existing) + len(completed) == len(tasks):
        _write_campaign_receipt(root, design, tasks)
    print(json.dumps({
        "segment": segment, "status": status, "completed": len(completed),
        "resumed": len(existing), "remaining":
            len(tasks) - len(existing) - len(completed),
        "terminal_sha256": sha256_file(root / "segments" /
            ("segment%04d.terminal.json" % segment)),
    }, sort_keys=True))


def run_campaign(args: argparse.Namespace) -> None:
    root = Path(args.result_dir).resolve()
    design = load_design(root)
    lock_path = root / "controller.lock"
    if lock_path.is_symlink() or not lock_path.is_file() or \
            sha256_file(lock_path) != design.get("controller_lock_sha256"):
        raise TimingError("controller lock identity changed")
    with lock_path.open("rb") as lock_stream:
        try:
            fcntl.flock(lock_stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise TimingError("another campaign controller owns this artifact") from exc
        _run_campaign_locked(args, root, design)


def _verify_interrupted_archive(root: Path,
                                records: Sequence[Mapping[str, object]]) -> None:
    for record in records:
        archive = record.get("archive")
        files = record.get("files")
        if not isinstance(archive, str) or not isinstance(files, dict):
            raise TimingError("interrupted transaction receipt malformed")
        directory = root / "interrupted" / archive
        if directory.is_symlink() or not directory.is_dir():
            raise TimingError("interrupted transaction archive missing")
        if sorted(path.name for path in directory.iterdir()) != sorted(files):
            raise TimingError("interrupted transaction archive file set changed")
        for name, receipt in files.items():
            path = directory / name
            if path.is_symlink() or not path.is_file() or \
                    path.stat().st_size != receipt.get("size") or \
                    sha256_file(path) != receipt.get("sha256"):
                raise TimingError("interrupted transaction archive changed")


def verify_segments(root: Path, design: Mapping[str, object],
                    task_entries: Mapping[int, Mapping[str, object]],
                    *, allow_unbound_archives: bool = False) -> List[Dict[str, object]]:
    terminals = sorted((root / "segments").glob("segment*.terminal.json"))
    records: List[Dict[str, object]] = []
    completed_jobs: set[int] = set()
    referenced_archives: set[str] = set()
    for index, path in enumerate(terminals):
        if path.name != "segment%04d.terminal.json" % index:
            raise TimingError("terminal segment sequence changed")
        terminal = verify_sealed(
            load_canonical(path, "segment terminal"),
            "wirehair.wh2.p32_dispatch.segment_terminal.v2",
            "segment terminal")
        terminal_fields = {
            "schema", "self_sha256_excluding_field", "segment", "status",
            "started_utc", "ended_utc", "timing_start_monotonic_s",
            "benchmark_end_monotonic_s", "duration_s", "controller_pid",
            "controller_affinity", "design_sha256", "prepare_receipt_sha256",
            "prelaunch_receipt_sha256", "tasks_manifest_sha256",
            "resumed_jobs_before", "completed_tasks",
            "recovered_interrupted_transactions", "failure",
            "sampler_identity_start", "sampler_identity_end", "graceful_stop",
            "graceful_stop_mechanism", "post_end_sample", "thermal_csv_name",
            "thermal_csv_sha256", "thermal_csv_size", "thermal_csv_device",
            "thermal_csv_inode", "thermal_summary",
        }
        if set(terminal) != terminal_fields:
            raise TimingError("segment terminal field set changed")
        if terminal.get("segment") != index or terminal.get("status") not in \
                ("failed", "bounded-incomplete", "complete"):
            raise TimingError("segment terminal identity/status malformed")
        status = str(terminal["status"])
        timing_start = terminal.get("timing_start_monotonic_s")
        benchmark_end = terminal.get("benchmark_end_monotonic_s")
        duration = terminal.get("duration_s")
        if type(timing_start) not in (int, float) or \
                type(benchmark_end) not in (int, float) or \
                type(duration) not in (int, float) or \
                not all(math.isfinite(float(value))
                        for value in (timing_start, benchmark_end, duration)) or \
                float(timing_start) < 0.0 or \
                float(benchmark_end) <= float(timing_start) or \
                float(duration) != float(benchmark_end) - float(timing_start):
            raise TimingError("segment timing interval malformed")
        validate_utc_timestamp(terminal.get("started_utc"), "segment start")
        validate_utc_timestamp(terminal.get("ended_utc"), "segment end")
        if type(terminal.get("controller_pid")) is not int or \
                int(terminal["controller_pid"]) <= 1:
            raise TimingError("segment controller PID malformed")
        for key, file_name in (
                ("design_sha256", "design.json"),
                ("prepare_receipt_sha256", "prepare_receipt.json"),
                ("prelaunch_receipt_sha256", "prelaunch_receipt.json"),
                ("tasks_manifest_sha256", "tasks_manifest.jsonl")):
            if terminal.get(key) != sha256_file(root / file_name):
                raise TimingError("segment terminal input binding changed")
        if terminal.get("controller_affinity") != [design["controller_core"]]:
            raise TimingError("segment controller affinity changed")
        start_identity = terminal.get("sampler_identity_start")
        end_identity = terminal.get("sampler_identity_end")
        if not isinstance(start_identity, dict) or start_identity != end_identity:
            raise TimingError("sampler start/end identity changed")
        required_identity = ("pid", "start_tick", "boot_id", "cmdline",
                             "cmdline_sha256", "affinity", "csv_device", "csv_inode")
        if set(start_identity) != set(required_identity) or \
                start_identity.get("affinity") != [design["thermal_core"]]:
            raise TimingError("sampler terminal identity incomplete")
        if not isinstance(start_identity.get("pid"), int) or \
                int(start_identity["pid"]) <= 1 or \
                not isinstance(start_identity.get("start_tick"), int) or \
                int(start_identity["start_tick"]) <= 0 or \
                start_identity.get("csv_device") != terminal.get("thermal_csv_device") or \
                start_identity.get("csv_inode") != terminal.get("thermal_csv_inode"):
            raise TimingError("sampler PID/start/CSV identity malformed")
        cmdline = start_identity.get("cmdline")
        if cmdline != _sampler_cmdline(root, design, index, public=True) or \
                any("bootstrap.pid" in str(value) for value in cmdline) or \
                SHA256_RE.fullmatch(str(start_identity.get(
                    "cmdline_sha256", ""))) is None:
            raise TimingError("terminal cmdline did not redact ephemeral PID file")
        if terminal.get("graceful_stop") is not True or \
                terminal.get("graceful_stop_mechanism") != \
                    "sudo-python-pidfd-identity-verified-sigterm" or \
                terminal.get("post_end_sample") is not True:
            raise TimingError("sampler lifecycle was not terminal")
        thermal_name = terminal.get("thermal_csv_name")
        if thermal_name != "segment%04d.thermal.csv" % index:
            raise TimingError("terminal thermal file identity changed")
        verify_terminal_thermal(root / "segments" / str(thermal_name), terminal)
        completed = terminal.get("completed_tasks")
        resumed = terminal.get("resumed_jobs_before")
        if not isinstance(completed, list) or not isinstance(resumed, list) or \
                resumed != sorted(completed_jobs):
            raise TimingError("segment task ledger malformed")
        expected_pending = [job for job in range(int(design["task_count"]))
                            if job not in completed_jobs]
        completed_numbers = [item.get("job") if isinstance(item, dict) else None
                             for item in completed]
        if completed_numbers != expected_pending[:len(completed_numbers)]:
            raise TimingError("segment task completion order changed")
        for item in completed:
            if not isinstance(item, dict) or set(item) != {
                    "job", "receipt_sha256", "retry_count"} or \
                    not isinstance(item.get("job"), int):
                raise TimingError("segment completed task receipt malformed")
            job = item["job"]
            if job in completed_jobs or job not in task_entries or \
                    item.get("receipt_sha256") != task_entries[job]["receipt_sha256"] or \
                    item.get("retry_count") != task_entries[job]["receipt"]["attempt"]:
                raise TimingError("segment completed task binding changed")
            completed_jobs.add(job)
        failure = terminal.get("failure")
        if status == "failed":
            if not isinstance(failure, dict) or set(failure) != {"class", "message"} or \
                    not all(isinstance(failure[key], str) and failure[key]
                            for key in ("class", "message")):
                raise TimingError("failed segment lacks exact failure receipt")
        elif failure is not None:
            raise TimingError("successful segment carries a failure receipt")
        if status == "complete" and len(completed_jobs) != int(design["task_count"]):
            raise TimingError("complete segment did not finish the manifest")
        if status == "bounded-incomplete" and \
                len(completed_jobs) >= int(design["task_count"]):
            raise TimingError("bounded segment status contradicts task ledger")
        recovered = terminal.get("recovered_interrupted_transactions")
        if not isinstance(recovered, list):
            raise TimingError("segment interrupted recovery ledger malformed")
        _verify_interrupted_archive(root, recovered)
        for item in recovered:
            archive = str(item["archive"])
            if archive in referenced_archives:
                raise TimingError("interrupted archive bound more than once")
            referenced_archives.add(archive)
        records.append(terminal)
    segment_files = sorted(path.name for path in (root / "segments").iterdir())
    expected_files = sorted(
        ["segment%04d.terminal.json" % i for i in range(len(terminals))] +
        ["segment%04d.thermal.csv" % i for i in range(len(terminals))])
    if segment_files != expected_files:
        raise TimingError("segments directory contains ephemeral or unknown files")
    actual_archives = {path.name for path in (root / "interrupted").iterdir()}
    if not allow_unbound_archives and actual_archives != referenced_archives:
        raise TimingError("interrupted archive ledger is not exact")
    if allow_unbound_archives and not referenced_archives <= actual_archives:
        raise TimingError("terminal-referenced interrupted archive is missing")
    if allow_unbound_archives:
        if not completed_jobs <= set(task_entries):
            raise TimingError("terminal binds a missing task entry")
    elif completed_jobs != set(task_entries):
        raise TimingError("task ledger is not exactly terminal-bound")
    return records


def verify_campaign_data(root: Path, *, require_complete: bool) -> Tuple[
        Dict[str, object], List[Dict[str, object]], Dict[int, Dict[str, object]],
        List[Dict[str, object]]]:
    design = load_design(root)
    verify_root_layout(root)
    tasks = load_tasks(root, design)
    verify_immutable_inputs(root, design)
    if any((root / ".transactions").iterdir()):
        raise TimingError("unrecovered task transaction remains")
    entries = existing_task_entries(root, design, tasks)
    segments = verify_segments(root, design, entries)
    complete = len(entries) == len(tasks)
    campaign_path = root / "campaign_receipt.json"
    if require_complete and not complete:
        raise TimingError("campaign task ledger is incomplete")
    if complete:
        if not campaign_path.is_file():
            raise TimingError("complete task ledger lacks campaign receipt")
        campaign = verify_sealed(
            load_canonical(campaign_path, "campaign receipt"),
            "wirehair.wh2.p32_dispatch.campaign.v2", "campaign receipt")
        if set(campaign) != {
                "schema", "self_sha256_excluding_field", "completed_utc",
                "design_sha256", "prepare_receipt_sha256",
                "prelaunch_receipt_sha256", "tasks_manifest_sha256",
                "task_count", "task_receipts", "terminal_receipts"}:
            raise TimingError("campaign receipt field set changed")
        expected_tasks = [{"job": task["job"],
                           "receipt_sha256": entries[int(task["job"])]["receipt_sha256"]}
                          for task in tasks]
        expected_terminals = [{"name": "segment%04d.terminal.json" % index,
                               "sha256": sha256_file(root / "segments" /
                                   ("segment%04d.terminal.json" % index))}
                              for index in range(len(segments))]
        exact = {
            "design_sha256": sha256_file(root / "design.json"),
            "prepare_receipt_sha256": sha256_file(root / "prepare_receipt.json"),
            "prelaunch_receipt_sha256": sha256_file(
                root / "prelaunch_receipt.json"),
            "tasks_manifest_sha256": sha256_file(root / "tasks_manifest.jsonl"),
            "task_count": len(tasks), "task_receipts": expected_tasks,
            "terminal_receipts": expected_terminals,
        }
        for key, value in exact.items():
            if campaign.get(key) != value:
                raise TimingError("campaign receipt mismatch: " + key)
        validate_utc_timestamp(campaign.get("completed_utc"),
                               "campaign completion")
    elif campaign_path.exists():
        raise TimingError("incomplete task ledger has a campaign receipt")
    return design, tasks, entries, segments


def _aggregate_rows(tasks: Sequence[Mapping[str, object]],
                    entries: Mapping[int, Mapping[str, object]]) -> Dict[str, object]:
    cells: List[Dict[str, object]] = []
    route_counts: Dict[str, int] = {}
    cache_identity: Dict[Tuple[object, ...], Tuple[object, ...]] = {}
    work_fields = (
        "inactivated", "binary_def", "heavy_gain", "block_xors",
        "block_muladds", "joint_source_xors", "joint_marginal_xors",
        "joint_marginal_copies", "joint_active_deltas", "joint_scratch_bytes",
        "dual_source_columns", "intermediate_bytes",
    )
    for task in tasks:
        parsed: ParsedOutput = entries[int(task["job"])]["parsed"]
        route = parsed.preamble["candidate_preflight_rhs_route"]
        route_counts[route] = route_counts.get(route, 0) + 1
        control_row = next(row for row in parsed.rows if row["arm"] == "control")
        candidate_row = next(row for row in parsed.rows if row["arm"] == "candidate")
        control_success = parsed.cell_class in ("common-success", "control-only")
        candidate_success = parsed.cell_class in ("common-success", "candidate-only")
        if parsed.common_success and (parsed.timed_control_ns <= 0 or
                                      parsed.timed_candidate_ns <= 0):
            raise TimingError("common-success timing sum is empty")
        identity_key = (task["K"], task["bb"], task["seed_index"],
                        task["schedule"])
        identity = (
            parsed.cell_class, parsed.work_signatures["control"],
            parsed.work_signatures["candidate"], parsed.trace_sha256,
        )
        if identity_key in cache_identity and cache_identity[identity_key] != identity:
            raise TimingError("cold/warm recovery or work identity changed")
        cache_identity[identity_key] = identity
        cells.append({
            "job": task["job"], "K": task["K"], "bb": task["bb"],
            "seed_index": task["seed_index"], "seed": task["seed"],
            "schedule": task["schedule"], "cache_state": task["cache_state"],
            "cell_class": parsed.cell_class, "control_success": control_success,
            "candidate_success": candidate_success,
            "control_ns": parsed.timed_control_ns if parsed.common_success else None,
            "candidate_ns": parsed.timed_candidate_ns if parsed.common_success else None,
            "ratio": parsed.timed_candidate_ns / parsed.timed_control_ns
                if parsed.common_success else None,
            "control_work": {field: int(control_row[field]) for field in work_fields},
            "candidate_work": {field: int(candidate_row[field]) for field in work_fields},
        })
    def summarize(values: Sequence[Mapping[str, object]]) -> Dict[str, object]:
        common = [value for value in values
                  if value["cell_class"] == "common-success"]
        control = sum(int(value["control_ns"]) for value in common)
        candidate = sum(int(value["candidate_ns"]) for value in common)
        ratios = [float(value["ratio"]) for value in common]
        outcomes = {name: sum(value["cell_class"] == name for value in values)
                    for name in ("common-success", "control-only",
                                 "candidate-only", "common-failure")}
        work_totals = {
            arm: {field: sum(int(value[arm + "_work"][field]) for value in values)
                  for field in work_fields}
            for arm in ("control", "candidate")
        }
        return {
            "cell_count": len(values), "outcome_counts": outcomes,
            "control_success_count": sum(bool(value["control_success"])
                                         for value in values),
            "candidate_success_count": sum(bool(value["candidate_success"])
                                           for value in values),
            "recovery_advantage_count":
                outcomes["candidate-only"] - outcomes["control-only"],
            "common_success_timing_count": len(common), "control_ns": control,
            "candidate_ns": candidate,
            "ratio_of_sums": candidate / control if control else None,
            "median_cell_ratio": statistics.median(ratios) if ratios else None,
            "candidate_faster_cells": sum(value < 1.0 for value in ratios),
            "work_totals": work_totals,
        }
    groups: Dict[str, List[Mapping[str, object]]] = {}
    for cell in cells:
        keys = (
            "K=%s" % cell["K"], "bb=%s" % cell["bb"],
            "schedule=%s" % cell["schedule"],
            "seed_index=%s" % cell["seed_index"],
            "cache=%s" % cell["cache_state"],
            "K=%s,bb=%s" % (cell["K"], cell["bb"]),
            "K=%s,bb=%s,schedule=%s" %
                (cell["K"], cell["bb"], cell["schedule"]),
            "K=%s,bb=%s,cache=%s" %
                (cell["K"], cell["bb"], cell["cache_state"]),
        )
        for key in keys:
            groups.setdefault(key, []).append(cell)
    return {
        "candidate_route_counts": route_counts,
        "weak_seed_coordinates": [{
            "K": value["K"], "bb": value["bb"],
            "seed_index": value["seed_index"], "seed": value["seed"],
            "schedule": value["schedule"], "cell_class": value["cell_class"],
        } for value in cells if value["cache_state"] == "cold" and
            value["cell_class"] != "common-success"],
        "overall": summarize(cells),
        "groups": {key: summarize(value) for key, value in sorted(groups.items())},
    }


def reduce_campaign(args: argparse.Namespace) -> None:
    root = Path(args.result_dir).resolve()
    design, tasks, entries, segments = verify_campaign_data(
        root, require_complete=True)
    payload = {
        "created_utc": utc_now(), "design_sha256": sha256_file(root / "design.json"),
        "campaign_receipt_sha256": sha256_file(root / "campaign_receipt.json"),
        "task_count": len(tasks), "segment_count": len(segments),
        "comparison": "p32_r7_vs_prod244", "timing_promotional": True,
        "seed_fixes": "none", "aggregates": _aggregate_rows(tasks, entries),
    }
    summary = sealed_record("wirehair.wh2.p32_dispatch.summary.v2", payload)
    path = root / "validated_summary.json"
    write_new(path, canonical_json(summary))
    print(json.dumps({"validated_summary_sha256": sha256_file(path),
                      "overall": summary["aggregates"]["overall"]},
                     sort_keys=True))


def verify_campaign(args: argparse.Namespace) -> None:
    root = Path(args.result_dir).resolve()
    design, tasks, entries, segments = verify_campaign_data(
        root, require_complete=not args.allow_incomplete)
    summary_path = root / "validated_summary.json"
    if summary_path.exists():
        if len(entries) != len(tasks):
            raise TimingError("incomplete campaign has validated summary")
        summary = verify_sealed(
            load_canonical(summary_path, "validated summary"),
            "wirehair.wh2.p32_dispatch.summary.v2", "validated summary")
        if summary.get("design_sha256") != sha256_file(root / "design.json") or \
                summary.get("campaign_receipt_sha256") != sha256_file(
                    root / "campaign_receipt.json") or \
                summary.get("aggregates") != _aggregate_rows(tasks, entries):
            raise TimingError("validated summary does not replay exactly")
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
    prepare.add_argument("--result-dir", required=True)
    prepare.add_argument("--core", type=int, default=8)
    prepare.add_argument("--controller-core", type=int, default=126)
    prepare.add_argument("--thermal-core", type=int, default=127)
    prepare.add_argument("--numa-node", type=int, default=0)
    prepare.add_argument("--evict-bytes", type=int, default=268435456)
    prepare.add_argument("--build-jobs", type=int, default=12)
    prepare.add_argument("--c-compiler")
    prepare.add_argument("--cxx-compiler")
    run = subparsers.add_parser("run")
    run.add_argument("--result-dir", required=True)
    run.add_argument("--max-tasks", type=int, default=1404)
    run.add_argument("--timeout-seconds", type=float, default=1800.0)
    reduce_parser = subparsers.add_parser("reduce")
    reduce_parser.add_argument("--result-dir", required=True)
    verify = subparsers.add_parser("verify")
    verify.add_argument("--result-dir", required=True)
    verify.add_argument("--allow-incomplete", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> None:
    args = build_parser().parse_args(argv)
    if args.command == "prepare":
        prepare_campaign(args)
    elif args.command == "run":
        run_campaign(args)
    elif args.command == "reduce":
        reduce_campaign(args)
    elif args.command == "verify":
        verify_campaign(args)
    else:
        raise TimingError("unknown command")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("P32 timing failed: %s" % exc, file=sys.stderr, flush=True)
        raise
