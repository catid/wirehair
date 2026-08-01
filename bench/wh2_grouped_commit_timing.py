#!/usr/bin/env python3
"""Sealed cross-binary timing for the WH2 grouped decoder commit delta.

This formal campaign compares the exact schedule-cache baseline and the
production-only combined RHS-fusion rescue.  Each parent is built with the
same measurement-only ``groupedtiming`` overlay from commit d6ab1a6, then run
with the same P48/r3 graph, packet trace, payload, allocator policy, and CPU
placement.  Each cell launches the binaries
in an outer ABBABAAB order.  Every launch in turn runs the codec's four-cycle
inner ABBABAAB fixture; inner cycle zero is discarded.  Thus process-start
drift and within-process drift are both balanced without timing system
construction or payload allocation.

The ``prepare`` command checks out both exact codec parents in detached
temporary worktrees, applies and receipts the identical benchmark overlay,
freezes the binaries and build provenance, and creates an immutable 108-cell
task ledger.  ``run`` is intentionally fresh-only and requires a live
sole-reader CPU/DIMM thermal sampler.  ``reduce`` independently replays every
raw output and receipt before publishing elapsed and phase-timing summaries.
No command writes a production profile or changes the caller's worktree.
"""

from __future__ import annotations

import argparse
import ctypes
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
import select
import selectors
import shutil
import signal
import stat
import subprocess
import sys
import time
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import wh2_timing_evidence_io as evidence_io


sys.dont_write_bytecode = True

BASE_COMMIT = "48d14bc77e3f9e98605fca4d226aa218d7d03a0d"
CANDIDATE_COMMIT = "be0bc94b97d03d6ddbc23db3b7058aa7f575b5cd"
MEASUREMENT_OVERLAY_COMMIT = "d6ab1a65c9864ad97ef06c4b88c2917cb387c0be"
MEASUREMENT_OVERLAY_PARENT_COMMIT = \
    "3a659aeb132e6cf5e5e68d88094af8402cdb0e47"
MEASUREMENT_OVERLAY_FILES = (
    "codec/WirehairV2Bench.cpp",
    "codec/V2BenchCliTest.cmake",
)
MEASUREMENT_DIFF_OPTIONS = (
    "--binary", "--full-index", "--no-color", "--no-ext-diff",
    "--no-textconv", "--no-renames", "--src-prefix=a/", "--dst-prefix=b/",
    "--diff-algorithm=myers", "--no-indent-heuristic", "--unified=3",
)
KS = (3200, 10000, 20000, 32000, 48466, 64000)
WIDTHS = (64, 1280, 4096)
SCHEDULE_SEEDS = (
    ("burst", 15111065706836454659),
    ("adversarial", 10723151780598845931),
    ("repair-only", 9599682871048892067),
)
CACHE_STATES = ("cold", "warm")
OUTER_ORDER = "ABBABAAB"
INNER_ORDER = "ABBABAAB"
ARCHITECTURE = {"period": 48, "grouped_rows": 3, "buckets": "auto"}
BINARY_NAMES = {"base": "wirehair_v2_bench.A", "candidate": "wirehair_v2_bench.B"}
THERMAL_SAMPLER_NAME = "wirehair_expo_thermal_sampler.py"
PREPARED_CAMPAIGN_DIRECTORIES = (
    "frozen", "provenance", "raw", "stderr", "exit", "receipts",
    "task_receipts", "contamination", "failure",
)
OVERHEAD = 4
LOSS_TEXT = "0.5"
MALLOC_MMAP_THRESHOLD = "1073741824"
MALLOC_TRIM_THRESHOLD = "-1"
MAX_ENVIRONMENTAL_ATTEMPTS = 10
MAX_MINOR_FAULTS = 64
SIBLING_PREFLIGHT_WINDOW_NS = 5_000_000_000
SIBLING_PREFLIGHT_MAX_BUSY_TICKS = 0
SIBLING_ACCEPTED_EXECUTION_MAX_BUSY_TICKS = 0
SIBLING_SCHED_RUNTIME_MAX_PPM = 50
SIBLING_SCHED_PCOUNT_WINDOW_NS = 1_000_000_000
DEVICE_SOFTIRQ_VECTORS = frozenset((
    "HI", "NET_TX", "NET_RX", "BLOCK", "IRQ_POLL", "TASKLET",
))
EXPECTED_SOFTIRQ_VECTORS = (
    "HI", "TIMER", "NET_TX", "NET_RX", "BLOCK", "IRQ_POLL", "TASKLET",
    "SCHED", "HRTIMER", "RCU",
)
GLOBAL_HARDIRQ_VECTORS = frozenset(("ERR", "MIS"))
MANAGED_NVME_IRQ_WHITELIST = (
    (103, "IR-PCI-MSIX-0000:a9:00.0 9-edge nvme2q9",
     "nvme2q9", "8", "8"),
    (107, "IR-PCI-MSIX-0000:d1:00.0 2-edge nvme0q2",
     "nvme0q2", "5-9,69-72", "72"),
    (134, "IR-PCI-MSIX-0000:d2:00.0 14-edge nvme1q14",
     "nvme1q14", "60-63,124-127", "126"),
    (176, "IR-PCI-MSIX-0000:aa:00.0 9-edge nvme3q9",
     "nvme3q9", "8", "8"),
    (257, "IR-PCI-MSIX-0000:a9:00.0 73-edge nvme2q73",
     "nvme2q73", "72", "72"),
    (304, "IR-PCI-MSIX-0000:aa:00.0 73-edge nvme3q73",
     "nvme3q73", "72", "72"),
    (365, "IR-PCI-MSIX-0000:a9:00.0 127-edge nvme2q127",
     "nvme2q127", "126", "126"),
    (390, "IR-PCI-MSIX-0000:aa:00.0 127-edge nvme3q127",
     "nvme3q127", "126", "126"),
)
IRQ30_IDENTITY = "IOMMU-MSI 376-edge AMD-Vi0-PPR"
GUARDED_IRQ_IDENTITIES = {
    30: IRQ30_IDENTITY,
    **{
        irq: identity
        for irq, identity, _handler, _requested, _effective
        in MANAGED_NVME_IRQ_WHITELIST
    },
}
TMPFS_MAGIC = 0x01021994
SIBLING_CAMPAIGN_MAX_BUSY_PPM = 50
SIBLING_CAMPAIGN_MIN_BUSY_TICKS = 1
CLOCK_TICKS_PER_SECOND = int(os.sysconf("SC_CLK_TCK"))
SIBLING_ACCEPTED_EXECUTION_MIN_DURATION_NS = \
    ((1_000_000_000 + CLOCK_TICKS_PER_SECOND - 1) //
     CLOCK_TICKS_PER_SECOND if CLOCK_TICKS_PER_SECOND > 0 else 0)
SCHEDSTAT_VERSION = 15
SCHEDSTAT_CPU_FIELD_COUNT = 9
SCHEDSTAT_DOMAIN_FIELD_COUNT = 36
SCHEDSTAT_RUNTIME_FIELD = 7
SCHEDSTAT_PCOUNT_FIELD = 9
SIBLING_IDLE_POLICY = {
    "accounting_sources": [
        "/proc/stat", "/proc/schedstat", "/proc/interrupts",
        "/proc/softirqs"],
    "clock_ticks_per_second": CLOCK_TICKS_PER_SECOND,
    "schedstat_version": SCHEDSTAT_VERSION,
    "schedstat_runtime_field": SCHEDSTAT_RUNTIME_FIELD,
    "schedstat_pcount_field": SCHEDSTAT_PCOUNT_FIELD,
    "schedstat_sched_schedstats_sysctl_required": False,
    "preflight_window_ns": SIBLING_PREFLIGHT_WINDOW_NS,
    "preflight_max_busy_ticks": SIBLING_PREFLIGHT_MAX_BUSY_TICKS,
    "preflight_max_sched_runtime_ppm": SIBLING_SCHED_RUNTIME_MAX_PPM,
    "preflight_max_sched_pcount_per_window": 1,
    "preflight_sched_pcount_window_ns": SIBLING_SCHED_PCOUNT_WINDOW_NS,
    "accepted_execution_max_busy_ticks":
        SIBLING_ACCEPTED_EXECUTION_MAX_BUSY_TICKS,
    "accepted_execution_max_sched_runtime_ppm":
        SIBLING_SCHED_RUNTIME_MAX_PPM,
    "accepted_execution_max_sched_pcount_per_window": 1,
    "accepted_execution_sched_pcount_window_ns":
        SIBLING_SCHED_PCOUNT_WINDOW_NS,
    "accepted_execution_min_duration_ns":
        SIBLING_ACCEPTED_EXECUTION_MIN_DURATION_NS,
    "accepted_execution_hybrid_proof":
        "at-least-one-USER_HZ-tick with zero proc-stat busy ticks, "
        "at most 50-ppm schedstat runtime, and at most one schedstat "
        "pcount per started one-second window",
    "campaign_max_busy_ppm": SIBLING_CAMPAIGN_MAX_BUSY_PPM,
    "campaign_min_busy_ticks": SIBLING_CAMPAIGN_MIN_BUSY_TICKS,
    "campaign_max_sched_runtime_ppm": SIBLING_CAMPAIGN_MAX_BUSY_PPM,
    "numeric_hardirq_target_delta": 0,
    "managed_numeric_irq_whitelist": [
        {"irq": irq, "handler": handler, "requested": requested,
         "effective": effective}
        for irq, _identity, handler, requested, effective
        in MANAGED_NVME_IRQ_WHITELIST
    ],
    "device_softirq_timing_cpu_delta": 0,
    "global_hardirq_delta": 0,
    "named_hardirq_and_other_softirq_activity":
        "explicitly classified and receipted",
}
MAX_CPU_TEMP_C = 85.0
MAX_DIMM_TEMP_C = 90.0
MAX_THERMAL_GAP_S = 2.25
MAX_THERMAL_MARGIN_S = 2.25
BOOTSTRAP_REPS = 20000
MAX_GROUPED_OUTPUT_BYTES = 1024 * 1024
MAX_BENCHMARK_STDERR_BYTES = 64 * 1024
MAX_EXEC_HANDSHAKE_BYTES = 1024
MAX_EVIDENCE_FILE_BYTES = 1024 * 1024 * 1024
MAX_JSON_EVIDENCE_BYTES = 16 * 1024 * 1024
MAX_TASK_MANIFEST_BYTES = 16 * 1024 * 1024
MAX_THERMAL_CSV_BYTES = 256 * 1024 * 1024
MAX_SIGNED_COUNTER = (1 << 63) - 1
MAX_UNSIGNED_COUNTER = (1 << 64) - 1
MAX_CPU_ID = (1 << 31) - 1
MAX_CPU_LIST_CARDINALITY = 1 << 16
MAX_PROCESS_PID = (1 << 31) - 1
EXEC_BIND_TIMEOUT_S = 5.0
MAX_BENCHMARK_TIMEOUT_S = 24.0 * 60.0 * 60.0
FINAL_EXEC_MARKER = b"wirehair-final-exec-v1\n"
TIMING_C_FLAGS_RELEASE = \
    "-O3 -DNDEBUG -ffunction-sections -fdata-sections"
TIMING_CXX_FLAGS_RELEASE = (
    TIMING_C_FLAGS_RELEASE +
    " -DWIREHAIR_V2_GROUPED_TIMING_ONLY=1"
    " -DWIREHAIR_V2_DISABLE_PACKED_RESIDUAL_TEXT_SECTION=1"
    " -Wno-unused-function"
)
TIMING_EXE_LINKER_FLAGS = "-Wl,--gc-sections"
FORBIDDEN_TIMED_BINARY_SYMBOLS = (
    "CheckMixedRhsFusionOracleForTesting",
    "CheckPackedBinaryResidualOracleForTesting",
    "CmdSelfTest()", "CmdPrecodeFail(", "CmdCompare(",
    " U fork@", " U vfork@", " U clone@", " U clone3@",
    " U posix_spawn@", " U popen@", " U pthread_create@",
    " U system@", " U syscall@",
)
CHILD_EXEC_WRAPPER = r"""
import ctypes
import fcntl
import os
import signal
import sys

fd = int(sys.argv[1])
expected_parent = int(sys.argv[2])
command = sys.argv[3:]
try:
    libc = ctypes.CDLL(None, use_errno=True)
    if libc.prctl(1, signal.SIGKILL, 0, 0, 0) != 0:
        raise OSError(ctypes.get_errno(), "prctl(PR_SET_PDEATHSIG)")
    if os.getppid() != expected_parent:
        raise RuntimeError("controller parent changed before exec")
    fcntl.fcntl(fd, fcntl.F_SETFD, 0)
    os.execv(command[0], command)
except BaseException as exc:
    payload = ("pre-exec:%s:%s\n" %
               (type(exc).__name__, str(exc))).encode("utf-8", "replace")
    try:
        os.write(fd, payload[:1024])
    except OSError:
        pass
    os._exit(126)
"""
FINAL_EXEC_WRAPPER = r"""
import fcntl
import os
import sys

fd = int(sys.argv[1])
expected_parent = int(sys.argv[2])
command = sys.argv[3:]
try:
    if os.getppid() != expected_parent:
        raise RuntimeError("controller parent changed before final exec")
    os.write(fd, b"wirehair-final-exec-v1\n")
    fcntl.fcntl(fd, fcntl.F_SETFD, fcntl.FD_CLOEXEC)
    os.execv(command[0], command)
except BaseException as exc:
    payload = ("final-exec:%s:%s\n" %
               (type(exc).__name__, str(exc))).encode("utf-8", "replace")
    try:
        os.write(fd, payload[:1024])
    except OSError:
        pass
    os._exit(126)
"""
PHASE_FIELDS = (
    "build_ns", "peel_ns", "project_ns", "residual_ns", "backsub_ns",
)
PRIMARY_PHASE_FIELD = "residual_ns"
NEGATIVE_CONTROL_PHASE_FIELDS = (
    "build_ns", "peel_ns", "project_ns", "backsub_ns",
)

SHA256_RE = re.compile(r"[0-9a-f]{64}\Z")
GIT_OID_RE = re.compile(r"[0-9a-f]{40,64}\Z")
UINT_RE = re.compile(r"0|[1-9][0-9]*\Z")
SINT_RE = re.compile(r"0|-?[1-9][0-9]*\Z")
HEX_RE = re.compile(r"0x(?:0|[1-9a-f][0-9a-f]*)\Z")
DECIMAL_RE = re.compile(r"-?(?:0|[1-9][0-9]*)(?:\.[0-9]+)?\Z")
UTC_RE = re.compile(
    r"[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:"
    r"[0-9]{2}\.[0-9]{3}Z\Z")
BOOT_ID_RE = re.compile(
    r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-"
    r"[0-9a-f]{4}-[0-9a-f]{12}\Z")

GROUPED_FIELDS = (
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
    "packet_payload_bytes", "intermediate_bytes",
)

GROUPED_PREAMBLE_FIELDS = (
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
    "allocator_tls_state", "preflight_control_result",
    "preflight_candidate_result", "cell_class", "common_success",
    "trace_sha256",
)

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
THERMAL_READER_FIELDS = frozenset((
    "pid", "start_ticks", "boot_id", "cpus_allowed_list", "argv",
    "argv_sha256", "sampler_sha256", "thermal_csv", "pid_file",
))

EXECUTION_RECEIPT_FIELDS = frozenset((
    "schema", "self_sha256_excluding_field", "job", "task_id", "task_sha256",
    "outer_slot", "outer_marker", "binary_label", "binary_sha256", "argv",
    "attempt", "started_utc", "start_monotonic_ns", "end_monotonic_ns",
    "duration_ns", "stderr_sha256", "prior_contamination_receipts",
    "schema_version", "stdout_sha256", "semantic_sha256", "trace_sha256",
    "work_signature_sha256", "timed_elapsed_ns", "all_elapsed_ns",
    "timed_phase_ns", "all_phase_ns",
    "timed_minor_faults", "discard_minor_faults", "max_minor_faults",
    "row_count", "timed_row_count", "process_identity", "cleanup_action",
    "sibling_ticks_before", "sibling_ticks_after", "sibling_busy_ticks",
    "sibling_schedstat_before", "sibling_schedstat_after",
    "sibling_sched_runtime_ns", "sibling_sched_runtime_limit_ns",
    "sibling_sched_pcount", "sibling_sched_pcount_limit",
    "target_irq_snapshot_before", "target_irq_snapshot_after",
    "target_irq_delta",
))
FAILURE_RECEIPT_FIELDS = frozenset((
    "schema", "self_sha256_excluding_field", "job", "task_id",
    "task_sha256", "outer_slot", "outer_marker", "binary_label",
    "binary_sha256", "argv", "attempt", "started_utc",
    "start_monotonic_ns", "end_monotonic_ns", "duration_ns",
    "returncode", "error_type", "error_message", "stdout_name",
    "stdout_sha256", "stderr_name", "stderr_sha256", "process_identity",
    "cleanup_action", "cleanup_error",
    "sibling_ticks_before", "sibling_ticks_after", "sibling_busy_ticks",
    "sibling_schedstat_before", "sibling_schedstat_after",
    "sibling_sched_runtime_ns", "sibling_sched_runtime_limit_ns",
    "sibling_sched_pcount", "sibling_sched_pcount_limit",
    "target_irq_snapshot_before", "target_irq_snapshot_after",
    "target_irq_delta",
))
CONTAMINATION_RECEIPT_FIELDS = frozenset((
    "schema", "self_sha256_excluding_field", "name", "attempt", "argv",
    "process_identity", "cleanup_action",
    "start_monotonic_ns", "end_monotonic_ns", "duration_ns", "stdout_name",
    "stdout_sha256", "stderr_name", "stderr_sha256", "contaminations",
    "sibling_ticks_before", "sibling_ticks_after", "sibling_busy_ticks",
    "sibling_schedstat_before", "sibling_schedstat_after",
    "sibling_sched_runtime_ns", "sibling_sched_runtime_limit_ns",
    "sibling_sched_pcount", "sibling_sched_pcount_limit",
    "target_irq_snapshot_before", "target_irq_snapshot_after",
    "target_irq_delta",
))
TASK_RECEIPT_FIELDS = frozenset((
    "schema", "self_sha256_excluding_field", "job", "task_id", "task_sha256",
    "outer_order", "execution_receipts", "trace_sha256", "semantic_sha256",
    "work_signature_sha256", "base_timed_elapsed_ns",
    "candidate_timed_elapsed_ns", "ratio", "base_process_count",
    "candidate_process_count", "timed_rows_per_binary",
    "base_timed_phase_ns", "candidate_timed_phase_ns", "phase_ratios",
))
PREPARE_RECEIPT_FIELDS = frozenset((
    "schema", "self_sha256_excluding_field", "prepared_utc", "design_sha256",
    "tasks_manifest_sha256", "immutable_files", "base_binary_sha256",
    "candidate_binary_sha256",
))
LAUNCH_RECEIPT_FIELDS = frozenset((
    "schema", "self_sha256_excluding_field", "started_utc", "ended_utc",
    "start_monotonic_s", "end_monotonic_s", "duration_s",
    "start_monotonic_ns", "end_monotonic_ns", "duration_ns", "design_sha256",
    "prepare_receipt_sha256", "tasks_manifest_sha256", "task_count",
    "execution_count", "retry_count", "execution_receipts", "task_receipts",
    "thermal_reader", "thermal_source_device", "thermal_source_inode",
    "thermal_interval_sha256", "thermal_summary", "topology",
    "load_workers_stopped", "controller_core", "controller_affinity",
    "core_ticks_before", "core_ticks_after", "sibling_ticks_before",
    "sibling_ticks_after", "sibling_busy_ticks", "preflight_quiet_core_ticks",
    "preflight_quiet_sibling_ticks", "sibling_busy_limit_ticks",
    "sibling_attempt_busy_ticks", "sibling_gap_busy_ticks",
    "preflight_core_ticks_before", "preflight_core_ticks_after",
    "preflight_sibling_ticks_before", "preflight_sibling_ticks_after",
    "preflight_sibling_schedstat_before",
    "preflight_sibling_schedstat_after",
    "preflight_start_monotonic_ns", "preflight_end_monotonic_ns",
    "preflight_duration_ns", "preflight_sibling_sched_runtime_ns",
    "preflight_sibling_sched_runtime_limit_ns",
    "preflight_sibling_sched_pcount", "preflight_sibling_sched_pcount_limit",
    "sibling_schedstat_before", "sibling_schedstat_after",
    "sibling_sched_runtime_ns", "sibling_sched_runtime_limit_ns",
    "sibling_sched_pcount", "sibling_attempt_sched_runtime_ns",
    "sibling_gap_sched_runtime_ns", "sibling_attempt_sched_pcount",
    "sibling_gap_sched_pcount", "runtime_isolation_snapshot_start",
    "runtime_isolation_snapshot_start_sha256", "runtime_isolation_snapshot_end",
    "runtime_isolation_snapshot_end_sha256",
    "preflight_target_irq_snapshot_before",
    "preflight_target_irq_snapshot_after", "preflight_target_irq_delta",
    "campaign_target_irq_snapshot_before",
    "campaign_target_irq_snapshot_after", "campaign_target_irq_delta",
    "result_tmpfs_binding", "thermal_csv_tmpfs_binding",
    "thermal_pid_tmpfs_binding", "prepared_tree_tmpfs_bindings",
))

RUNTIME_ISOLATION_SNAPSHOT_FIELDS = frozenset((
    "schema", "capture_start_monotonic_ns", "capture_end_monotonic_ns",
    "capture_duration_ns", "self_cgroup", "expected_isolated_cpus",
    "kernel_isolated_cpu_list", "kernel_isolated_cpus",
    "cgroup_cpu_list", "cgroup_cpus", "cgroup_effective_cpu_list",
    "cgroup_effective_cpus", "cgroup_exclusive_cpu_list",
    "cgroup_exclusive_cpus", "cgroup_exclusive_effective_cpu_list",
    "cgroup_exclusive_effective_cpus", "cgroup_partition",
    "irq_effective_affinities", "irq30_exception",
    "managed_nvme_exceptions",
))
IRQ_EFFECTIVE_AFFINITY_FIELDS = frozenset((
    "irq", "effective_affinity_list", "effective_cpus",
))
IRQ30_EXCEPTION_FIELDS = frozenset((
    "irq", "identity", "handler_directories", "requested_affinity_list",
    "requested_cpus", "effective_affinity_list", "effective_cpus",
    "global_interrupt_count",
))
MANAGED_IRQ_EXCEPTION_FIELDS = frozenset((
    "irq", "identity", "handler_directories", "requested_affinity_list",
    "requested_cpus", "effective_affinity_list", "effective_cpus",
))
TARGET_IRQ_SNAPSHOT_FIELDS = frozenset((
    "schema", "capture_start_monotonic_ns", "capture_end_monotonic_ns",
    "capture_duration_ns", "target_cpus", "cpu_ids", "hardirq_rows",
    "hardirq_sha256", "softirq_rows", "softirq_sha256",
))
TARGET_IRQ_DELTA_FIELDS = frozenset((
    "schema", "target_cpus", "hardirq_before_sha256",
    "hardirq_after_sha256", "softirq_before_sha256",
    "softirq_after_sha256", "hardirq_deltas", "softirq_deltas",
    "classifications", "contaminations",
))
TMPFS_BINDING_FIELDS = frozenset((
    "schema", "path", "device", "inode", "mode", "uid", "gid", "nlink", "kind",
    "filesystem_magic",
))


class TimingError(RuntimeError):
    """A substantive campaign or evidence error."""


class ContaminationError(TimingError):
    """A predeclared, receipted environmental condition eligible for retry."""


@dataclass(frozen=True)
class BoundCommandResult:
    returncode: int
    stdout: bytes
    stderr: bytes
    process_identity: Mapping[str, object]
    cleanup_action: str


class BoundCommandError(TimingError):
    """A command failed after launch while retaining bounded evidence."""

    def __init__(
        self, message: str, *, stdout: bytes = b"", stderr: bytes = b"",
        returncode: Optional[int] = None,
        process_identity: Optional[Mapping[str, object]] = None,
        cleanup_action: str = "not_started",
        cleanup_error: Optional[str] = None,
    ) -> None:
        super().__init__(message)
        self.stdout = stdout
        self.stderr = stderr
        self.returncode = returncode
        self.process_identity = (
            dict(process_identity) if process_identity is not None else None)
        self.cleanup_action = cleanup_action
        self.cleanup_error = cleanup_error


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds").replace(
        "+00:00", "Z")


def canonical_json(value: object) -> bytes:
    return evidence_io.canonical_json_bytes(value, error_type=TimingError)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(
    path: Path, *, max_bytes: int = MAX_EVIDENCE_FILE_BYTES,
    require_unique: bool = False,
) -> str:
    return evidence_io.stable_file_sha256(
        path, max_bytes=max_bytes, require_unique=require_unique,
        error_type=TimingError)


def atomic_write(path: Path, data: bytes, mode: int = 0o444) -> None:
    evidence_io.publish_immutable_file(
        path, data, mode=mode, error_type=TimingError)


def write_new(path: Path, data: bytes, mode: int = 0o444) -> None:
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
    return evidence_io.load_canonical_object(
        path, max_bytes=MAX_JSON_EVIDENCE_BYTES, name=name,
        require_unique=True, error_type=TimingError)


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


def parse_decimal_float(
    text: object, name: str, minimum: float, maximum: float,
) -> float:
    if not isinstance(text, str) or DECIMAL_RE.fullmatch(text) is None:
        raise TimingError("%s is not a canonical decimal" % name)
    value = float(text)
    if not math.isfinite(value) or not minimum <= value <= maximum:
        raise TimingError("%s exceeds its numeric domain" % name)
    return value


def stable_bytes(
    path: Path, attempts: int = 20,
    max_bytes: int = MAX_EVIDENCE_FILE_BYTES,
    require_unique: bool = True,
) -> bytes:
    last_error: Optional[TimingError] = None
    for _attempt in range(attempts):
        try:
            return evidence_io.stable_file_snapshot(
                path, max_bytes=max_bytes, require_unique=require_unique,
                error_type=TimingError)
        except TimingError as exc:
            last_error = exc
            time.sleep(0.02)
    raise TimingError("file did not stabilize while reading: %s" % path) \
        from last_error


def require_pidfd_runtime() -> None:
    pidfd_open = getattr(os, "pidfd_open", None)
    if not callable(pidfd_open):
        raise TimingError(
            "campaign execution requires Python os.pidfd_open support")
    descriptor: Optional[int] = None
    try:
        descriptor = pidfd_open(os.getpid(), 0)
        if pidfd_has_exited(descriptor):
            raise TimingError("self pidfd unexpectedly reports exit")
    except (OSError, RuntimeError) as exc:
        raise TimingError(
            "campaign execution requires working kernel pidfds: %s" %
            exc) from exc
    finally:
        if descriptor is not None:
            os.close(descriptor)


def process_start_ticks(pid: int) -> Optional[int]:
    if not 1 < pid <= MAX_PROCESS_PID:
        return None
    try:
        raw = (Path("/proc") / str(pid) / "stat").read_bytes()
        return int(raw.rsplit(b") ", 1)[1].split()[19])
    except (OSError, IndexError, ValueError):
        return None


def process_tokens(pid: int) -> Optional[List[str]]:
    if not 1 < pid <= MAX_PROCESS_PID:
        return None
    try:
        raw = (Path("/proc") / str(pid) / "cmdline").read_bytes()
        return [
            token.decode("utf-8", "strict")
            for token in raw.split(b"\0") if token
        ]
    except (OSError, UnicodeDecodeError):
        return None


def process_executable_identity(pid: int) -> Optional[Dict[str, object]]:
    path = Path("/proc") / str(pid) / "exe"
    descriptor = -1
    try:
        descriptor = os.open(
            str(path), getattr(os, "O_PATH", os.O_RDONLY) |
            getattr(os, "O_CLOEXEC", 0))
        identity = os.fstat(descriptor)
        target = os.readlink(str(path))
    except OSError:
        return None
    finally:
        if descriptor >= 0:
            os.close(descriptor)
    return {
        "path": target, "device": identity.st_dev, "inode": identity.st_ino,
    }


def process_exec_snapshot(pid: int) -> Optional[Dict[str, object]]:
    executable = process_executable_identity(pid)
    tokens = process_tokens(pid)
    start_ticks = process_start_ticks(pid)
    if executable is None or tokens is None or start_ticks is None:
        return None
    return {
        "start_ticks": start_ticks, "executable": executable,
        "argv": tokens,
    }


def pidfd_has_exited(pidfd: int) -> bool:
    poller = select.poll()
    poller.register(pidfd, select.POLLIN | select.POLLHUP | select.POLLERR)
    return bool(poller.poll(0))


def _process_group_members(group: int) -> Tuple[Tuple[int, str], ...]:
    members = []
    for path in Path("/proc").glob("[0-9]*/stat"):
        try:
            raw = path.read_bytes()
            fields = raw.rsplit(b") ", 1)[1].split()
            if int(fields[2]) == group:
                members.append((
                    int(path.parent.name), fields[0].decode("ascii", "strict")))
        except (OSError, IndexError, ValueError, UnicodeDecodeError):
            continue
    return tuple(sorted(members))


def terminate_direct_child_by_pidfd(
    process: subprocess.Popen, pidfd: int,
) -> str:
    """SIGKILL and reap a private session while its leader pins PID/PGID."""
    if process.returncode is not None:
        live = [
            item for item in _process_group_members(process.pid)
            if item[1] != "Z"
        ]
        if live:
            raise TimingError(
                "reaped benchmark leader left live process-group members")
        return "already_reaped"
    leader_had_exited = pidfd_has_exited(pidfd)
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        # A zombie session leader still keeps the group addressable.  ESRCH is
        # therefore acceptable only when pidfd already proves it exited and no
        # group member remains.
        live = [
            item for item in _process_group_members(process.pid)
            if item[1] != "Z"
        ]
        if not leader_had_exited or live:
            raise TimingError(
                "owned benchmark process group vanished inconsistently")
    poller = select.poll()
    poller.register(pidfd, select.POLLIN | select.POLLHUP | select.POLLERR)
    if not poller.poll(5000):
        raise TimingError(
            "bound benchmark leader survived process-group SIGKILL")
    deadline = time.monotonic() + 5.0
    while time.monotonic() < deadline:
        live = [
            item for item in _process_group_members(process.pid)
            if item[0] != process.pid and item[1] != "Z"
        ]
        if not live:
            break
        time.sleep(0.01)
    else:
        raise TimingError("owned benchmark process group retained live members")
    try:
        process.wait(timeout=5.0)
    except subprocess.TimeoutExpired as exc:
        raise TimingError("bound benchmark leader could not be reaped") from exc
    live_after = [
        item for item in _process_group_members(process.pid)
        if item[1] != "Z"
    ]
    if live_after:
        raise TimingError(
            "owned benchmark process group survived leader reap")
    return "exited_group_swept" if leader_had_exited else "killed_group"


def terminate_direct_child_without_pidfd(
    process: subprocess.Popen,
) -> str:
    """Fallback sweep used only when pidfd acquisition itself failed."""
    leader_had_exited = process.returncode is not None
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        if _process_group_members(process.pid):
            raise TimingError(
                "unbound benchmark process group vanished inconsistently")
    try:
        process.wait(timeout=5.0)
    except subprocess.TimeoutExpired as exc:
        raise TimingError("unbound benchmark leader survived SIGKILL") from exc
    deadline = time.monotonic() + 5.0
    while time.monotonic() < deadline:
        live = [
            item for item in _process_group_members(process.pid)
            if item[1] != "Z"
        ]
        if not live:
            return (
                "exited_group_swept_without_pidfd"
                if leader_had_exited else "killed_group_without_pidfd")
        time.sleep(0.01)
    raise TimingError("unbound benchmark process group retained live members")


def _read_exec_handshake(descriptor: int, timeout_s: float) -> bytes:
    poller = select.poll()
    poller.register(
        descriptor, select.POLLIN | select.POLLHUP | select.POLLERR)
    deadline = time.monotonic() + timeout_s
    chunks: List[bytes] = []
    while True:
        remaining_ms = max(
            0, int(math.ceil((deadline - time.monotonic()) * 1000.0)))
        if remaining_ms == 0 or not poller.poll(remaining_ms):
            raise BoundCommandError("benchmark wrapper exec handshake timed out")
        try:
            chunk = os.read(
                descriptor,
                MAX_EXEC_HANDSHAKE_BYTES + 1 -
                sum(len(value) for value in chunks))
        except BlockingIOError:
            continue
        if not chunk:
            return b"".join(chunks)
        chunks.append(chunk)
        if sum(len(value) for value in chunks) > MAX_EXEC_HANDSHAKE_BYTES:
            raise BoundCommandError(
                "benchmark wrapper exec handshake exceeded its bound")


def _bind_benchmark_exec(
    process: subprocess.Popen, pidfd: int, start_ticks: int,
    binary: Path, expected_tokens: Sequence[str], deadline_s: float,
) -> Dict[str, object]:
    expected = binary.stat()
    expected_executable = {
        "path": str(binary), "device": expected.st_dev,
        "inode": expected.st_ino,
    }
    while time.monotonic() < deadline_s:
        if process_start_ticks(process.pid) != start_ticks:
            raise BoundCommandError(
                "benchmark process start identity changed")
        first = process_exec_snapshot(process.pid)
        if (first is not None and
                first["start_ticks"] == start_ticks and
                first["executable"] == expected_executable and
                first["argv"] == list(expected_tokens)):
            second = process_exec_snapshot(process.pid)
            if (second is not None and second == first and
                    not pidfd_has_exited(pidfd)):
                return {
                    "pid": process.pid, "start_ticks": start_ticks,
                    "executable": expected_executable,
                    "argv": list(expected_tokens),
                    "binding_observation": "double-proc-snapshot",
                }
        if pidfd_has_exited(pidfd):
            return {
                "pid": process.pid, "start_ticks": start_ticks,
                "executable": expected_executable,
                "argv": list(expected_tokens),
                "binding_observation": "final-exec-handshake-terminal",
            }
        time.sleep(min(0.002, max(0.0, deadline_s - time.monotonic())))
    raise BoundCommandError(
        "benchmark did not reach its exact final exec identity before timeout")


def _bounded_communicate(
    process: subprocess.Popen, pidfd: int, deadline_s: float,
) -> Tuple[bytes, bytes, str]:
    if process.stdout is None or process.stderr is None:
        raise TimingError("benchmark pipes were not captured")
    streams = {
        process.stdout.fileno(): ("stdout", MAX_GROUPED_OUTPUT_BYTES),
        process.stderr.fileno(): ("stderr", MAX_BENCHMARK_STDERR_BYTES),
    }
    chunks: Dict[str, List[bytes]] = {"stdout": [], "stderr": []}
    sizes = {"stdout": 0, "stderr": 0}
    cleanup_action: Optional[str] = None

    def failure(message: str) -> BoundCommandError:
        return BoundCommandError(
            message, stdout=b"".join(chunks["stdout"]),
            stderr=b"".join(chunks["stderr"]),
            cleanup_action=cleanup_action or "not_started")

    selector = selectors.DefaultSelector()
    try:
        for stream in (process.stdout, process.stderr):
            os.set_blocking(stream.fileno(), False)
            selector.register(stream, selectors.EVENT_READ)
        while selector.get_map():
            if cleanup_action is None and pidfd_has_exited(pidfd):
                cleanup_action = terminate_direct_child_by_pidfd(
                    process, pidfd)
            timeout = deadline_s - time.monotonic()
            if timeout <= 0:
                raise failure("timing subprocess exceeded its timeout")
            events = selector.select(min(timeout, 0.25))
            if not events and cleanup_action is not None:
                events = [
                    (key, selectors.EVENT_READ)
                    for key in list(selector.get_map().values())
                ]
            for key, _mask in events:
                try:
                    block = os.read(key.fd, 65536)
                except BlockingIOError:
                    continue
                if not block:
                    selector.unregister(key.fileobj)
                    continue
                name, maximum = streams[key.fd]
                prior_size = sizes[name]
                if prior_size + len(block) > maximum:
                    remaining = maximum - prior_size
                    if remaining > 0:
                        chunks[name].append(block[:remaining])
                    sizes[name] = maximum
                    raise failure(
                        "benchmark %s exceeded its bounded capture" % name)
                sizes[name] += len(block)
                chunks[name].append(block)
        if cleanup_action is None:
            remaining_ms = max(
                0, int(math.ceil(
                    (deadline_s - time.monotonic()) * 1000.0)))
            poller = select.poll()
            poller.register(
                pidfd, select.POLLIN | select.POLLHUP | select.POLLERR)
            if remaining_ms == 0 or not poller.poll(remaining_ms):
                raise failure("timing subprocess exceeded its timeout")
            cleanup_action = terminate_direct_child_by_pidfd(process, pidfd)
        return (
            b"".join(chunks["stdout"]), b"".join(chunks["stderr"]),
            cleanup_action)
    finally:
        selector.close()


def _drain_after_cleanup(
    process: subprocess.Popen, stdout_prefix: bytes, stderr_prefix: bytes,
) -> Tuple[bytes, bytes, Optional[str]]:
    """Best-effort bounded drain after the owned session has been stopped."""
    if process.stdout is None or process.stderr is None:
        return stdout_prefix, stderr_prefix, "benchmark capture pipes are missing"
    values = {
        "stdout": bytearray(stdout_prefix[:MAX_GROUPED_OUTPUT_BYTES]),
        "stderr": bytearray(stderr_prefix[:MAX_BENCHMARK_STDERR_BYTES]),
    }
    maxima = {
        "stdout": MAX_GROUPED_OUTPUT_BYTES,
        "stderr": MAX_BENCHMARK_STDERR_BYTES,
    }
    streams = {
        process.stdout.fileno(): "stdout",
        process.stderr.fileno(): "stderr",
    }
    selector = selectors.DefaultSelector()
    errors = []
    deadline = time.monotonic() + 1.0
    try:
        for stream in (process.stdout, process.stderr):
            try:
                os.set_blocking(stream.fileno(), False)
                selector.register(stream, selectors.EVENT_READ)
            except OSError as exc:
                errors.append("capture registration failed: %s" % exc)
        while selector.get_map() and time.monotonic() < deadline:
            events = selector.select(0.05)
            if not events:
                events = [
                    (key, selectors.EVENT_READ)
                    for key in list(selector.get_map().values())
                ]
            for key, _mask in events:
                try:
                    block = os.read(key.fd, 65536)
                except BlockingIOError:
                    continue
                except OSError as exc:
                    errors.append("capture drain failed: %s" % exc)
                    selector.unregister(key.fileobj)
                    continue
                if not block:
                    selector.unregister(key.fileobj)
                    continue
                name = streams[key.fd]
                remaining = maxima[name] - len(values[name])
                if remaining > 0:
                    values[name].extend(block[:remaining])
                if len(block) > remaining:
                    errors.append("%s capture exceeded its bound" % name)
        if selector.get_map():
            errors.append("capture pipes did not close after cleanup")
    finally:
        selector.close()
    return bytes(values["stdout"]), bytes(values["stderr"]), (
        "; ".join(errors) if errors else None)


def run_bound_command(
    command: Sequence[str], binary: Path, python: Path, timeout_s: float,
) -> BoundCommandResult:
    require_pidfd_runtime()
    if isinstance(timeout_s, bool) or not isinstance(timeout_s, (int, float)):
        raise TimingError("benchmark timeout must be a finite number")
    try:
        timeout_value = float(timeout_s)
    except (OverflowError, ValueError) as exc:
        raise TimingError("benchmark timeout must be a finite number") from exc
    if (not math.isfinite(timeout_value) or timeout_value <= 0 or
            timeout_value > MAX_BENCHMARK_TIMEOUT_S):
        raise TimingError("benchmark timeout is outside its bounded domain")
    command_list = list(command)
    if (not command_list or
            any(not isinstance(token, str) or not token or "\0" in token
                for token in command_list)):
        raise TimingError("benchmark command tokens are malformed")
    try:
        binary_index = command_list.index(str(binary))
    except ValueError as exc:
        raise TimingError("benchmark command does not contain its binary") from exc
    expected_tokens = command_list[binary_index:]
    final_wrapper = [
        str(python), "-I", "-c", FINAL_EXEC_WRAPPER,
        # The inherited handshake descriptor and controller PID are appended
        # after pipe creation below.
    ]
    read_fd, write_fd = os.pipe2(
        getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NONBLOCK", 0))
    wrapped_command = [
        *command_list[:binary_index], *final_wrapper,
        str(write_fd), str(os.getpid()), *expected_tokens,
    ]
    wrapper = [
        str(python), "-I", "-c", CHILD_EXEC_WRAPPER,
        str(write_fd), str(os.getpid()), *wrapped_command,
    ]
    process: Optional[subprocess.Popen] = None
    pidfd: Optional[int] = None
    identity: Optional[Dict[str, object]] = None
    cleanup_action = "not_started"
    deadline = time.monotonic() + timeout_value
    try:
        process = subprocess.Popen(
            wrapper, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            stdin=subprocess.DEVNULL, close_fds=True, pass_fds=(write_fd,),
            start_new_session=True,
            env={"LC_ALL": "C", "TZ": "UTC", "PATH": "/usr/bin:/bin"})
        os.close(write_fd)
        write_fd = -1
        start_ticks = process_start_ticks(process.pid)
        if start_ticks is None:
            raise BoundCommandError(
                "cannot bind benchmark direct-child start identity")
        identity = {
            "pid": process.pid, "start_ticks": start_ticks,
            "executable": None, "argv": None,
            "binding_observation": "direct-child-provisional",
            "boot_id": Path(
                "/proc/sys/kernel/random/boot_id").read_text(
                    encoding="ascii").strip(),
        }
        try:
            pidfd = os.pidfd_open(process.pid, 0)
        except OSError as exc:
            raise BoundCommandError(
                "cannot bind benchmark direct-child pidfd") from exc
        if process_start_ticks(process.pid) != start_ticks or \
                pidfd_has_exited(pidfd):
            raise BoundCommandError(
                "benchmark direct-child identity changed during pidfd binding")
        handshake = _read_exec_handshake(
            read_fd, min(EXEC_BIND_TIMEOUT_S, max(0.001, deadline - time.monotonic())))
        if handshake != FINAL_EXEC_MARKER:
            raise BoundCommandError(
                "benchmark final-exec handshake failed: %r" % handshake)
        identity = _bind_benchmark_exec(
            process, pidfd, start_ticks, binary, expected_tokens,
            min(deadline, time.monotonic() + EXEC_BIND_TIMEOUT_S))
        identity["boot_id"] = Path(
            "/proc/sys/kernel/random/boot_id").read_text(
                encoding="ascii").strip()
        stdout, stderr, cleanup_action = _bounded_communicate(
            process, pidfd, deadline)
        if (identity.get("binding_observation") ==
                "final-exec-handshake-terminal" and
                (process.returncode is None or process.returncode < 0 or
                 process.returncode == 126)):
            raise BoundCommandError(
                "terminal-only final-exec proof ended with an unsafe status",
                stdout=stdout, stderr=stderr,
                returncode=(
                    int(process.returncode)
                    if process.returncode is not None else None),
                process_identity=identity,
                cleanup_action=cleanup_action)
        return BoundCommandResult(
            returncode=int(process.returncode), stdout=stdout, stderr=stderr,
            process_identity=identity, cleanup_action=cleanup_action)
    except BoundCommandError as exc:
        cleanup_errors = []
        try:
            if exc.cleanup_action != "not_started":
                cleanup_action = exc.cleanup_action
            elif process is not None and pidfd is not None:
                cleanup_action = terminate_direct_child_by_pidfd(process, pidfd)
            elif process is not None:
                cleanup_action = terminate_direct_child_without_pidfd(process)
        except BaseException as cleanup_exc:
            cleanup_action = "cleanup_failed"
            cleanup_errors.append(
                "%s: %s" %
                (type(cleanup_exc).__name__, str(cleanup_exc)))
        if process is not None:
            try:
                stdout, stderr, drain_error = _drain_after_cleanup(
                    process, exc.stdout, exc.stderr)
                exc.stdout = stdout
                exc.stderr = stderr
                if drain_error is not None:
                    cleanup_errors.append(drain_error)
            except BaseException as drain_exc:
                cleanup_errors.append(
                    "%s: %s" %
                    (type(drain_exc).__name__, str(drain_exc)))
            if process.returncode is not None:
                exc.returncode = int(process.returncode)
        exc.process_identity = identity
        exc.cleanup_action = cleanup_action
        exc.cleanup_error = (
            "; ".join(cleanup_errors) if cleanup_errors else None)
        raise
    except BaseException as primary:
        cleanup_error = None
        try:
            if process is not None and pidfd is not None:
                terminate_direct_child_by_pidfd(process, pidfd)
            elif process is not None:
                terminate_direct_child_without_pidfd(process)
        except BaseException as cleanup_exc:
            cleanup_error = "%s: %s" % (
                type(cleanup_exc).__name__, str(cleanup_exc))
        if cleanup_error is not None:
            raise TimingError(
                "benchmark launch failed (%s: %s); cleanup also failed (%s)" %
                (type(primary).__name__, str(primary), cleanup_error)) from primary
        raise
    finally:
        # Close every local resource even if an earlier close is interrupted;
        # cleanup failures must never replace the substantive launch result.
        if write_fd >= 0:
            try:
                os.close(write_fd)
            except OSError:
                pass
        try:
            os.close(read_fd)
        except OSError:
            pass
        if pidfd is not None:
            try:
                os.close(pidfd)
            except OSError:
                pass
        if process is not None:
            if process.stdout is not None:
                try:
                    process.stdout.close()
                except OSError:
                    pass
            if process.stderr is not None:
                try:
                    process.stderr.close()
                except OSError:
                    pass


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
            if high - low + 1 > MAX_CPU_LIST_CARDINALITY - len(values):
                raise TimingError("CPU-list expands beyond its cardinality bound")
            values.extend(range(low, high + 1))
        else:
            if len(values) >= MAX_CPU_LIST_CARDINALITY:
                raise TimingError("CPU-list expands beyond its cardinality bound")
            values.append(parse_uint(token, "CPU-list value", (1 << 31) - 1))
    if values != sorted(set(values)):
        raise TimingError("CPU-list is duplicate or nonascending")
    return tuple(values)


def _bounded_kernel_bytes(path: Path, maximum: int = MAX_JSON_EVIDENCE_BYTES) -> bytes:
    """Read one procfs/sysfs pseudo-file without trusting its zero st_size."""
    if (not isinstance(maximum, int) or isinstance(maximum, bool) or
            maximum <= 0):
        raise TimingError("kernel pseudo-file byte limit is malformed")
    flags = os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW
    descriptor = -1
    try:
        descriptor = os.open(str(path), flags)
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            raise TimingError("kernel evidence path is not a plain file: %s" % path)
        blocks = []
        total = 0
        while True:
            block = os.read(descriptor, min(64 * 1024, maximum - total + 1))
            if not block:
                break
            total += len(block)
            if total > maximum:
                raise TimingError("kernel evidence exceeds byte limit: %s" % path)
            blocks.append(block)
        after = os.fstat(descriptor)
        named = os.stat(str(path), follow_symlinks=False)
        identity = lambda value: (
            value.st_dev, value.st_ino, value.st_mode, value.st_uid,
            value.st_gid)
        if identity(before) != identity(after) or identity(after) != identity(named):
            raise TimingError("kernel evidence identity changed: %s" % path)
        return b"".join(blocks)
    except OSError as exc:
        raise TimingError("cannot read kernel evidence: %s" % path) from exc
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _kernel_ascii_line(
    path: Path, context: str, *, allow_empty: bool = False,
) -> str:
    raw = _bounded_kernel_bytes(path, 4096)
    try:
        text = raw.decode("ascii", "strict")
    except UnicodeDecodeError as exc:
        raise TimingError("%s is not ASCII" % context) from exc
    if (not isinstance(allow_empty, bool) or not text.endswith("\n") or
            "\n" in text[:-1] or (not allow_empty and not text[:-1])):
        raise TimingError("%s is not one canonical line" % context)
    return text[:-1]


def _validate_tmpfs_binding_record(
    record: object, expected_path: Path, expected_kind: str, context: str,
) -> Dict[str, object]:
    expected = str(expected_path.resolve())
    if (not isinstance(record, dict) or set(record) != TMPFS_BINDING_FIELDS or
            record.get("schema") != "wirehair.wh2.tmpfs_binding.v1" or
            record.get("path") != expected or
            record.get("kind") != expected_kind or
            expected_kind not in ("directory", "regular") or
            record.get("filesystem_magic") != TMPFS_MAGIC):
        raise TimingError("%s tmpfs binding is malformed" % context)
    for field in ("device", "inode", "mode", "uid", "gid", "nlink"):
        value = record.get(field)
        if (not isinstance(value, int) or isinstance(value, bool) or value < 0):
            raise TimingError("%s tmpfs binding is malformed" % context)
    if (record["inode"] <= 0 or record["mode"] > 0o7777 or
            record["nlink"] <= 0 or
            (expected_kind == "regular" and record["nlink"] != 1)):
        raise TimingError("%s tmpfs binding is malformed" % context)
    return record


def validate_tmpfs_binding(
    record: object, expected_path: Path, expected_kind: str, context: str,
) -> Dict[str, object]:
    """Replay one binding against the currently named inode and filesystem."""
    checked = _validate_tmpfs_binding_record(
        record, expected_path, expected_kind, context)
    current = capture_tmpfs_binding(expected_path, expected_kind, context)
    if current != checked:
        raise TimingError("%s live tmpfs binding changed" % context)
    return checked


def capture_tmpfs_binding(
    path: Path, expected_kind: str, context: str,
) -> Dict[str, object]:
    """Bind one exact inode to Linux tmpfs without launching a helper."""
    resolved = path.resolve()
    if path != resolved or expected_kind not in ("directory", "regular"):
        raise TimingError("%s path is not canonical" % context)
    flags = os.O_RDONLY | os.O_CLOEXEC | os.O_NOFOLLOW
    if expected_kind == "directory":
        flags |= os.O_DIRECTORY
    descriptor = -1
    try:
        descriptor = os.open(str(resolved), flags)
        before = os.fstat(descriptor)
        expected_type = stat.S_ISDIR if expected_kind == "directory" \
            else stat.S_ISREG
        if not expected_type(before.st_mode):
            raise TimingError("%s has the wrong inode type" % context)
        if expected_kind == "regular" and before.st_nlink != 1:
            raise TimingError("%s link policy changed" % context)
        buffer = ctypes.create_string_buffer(256)
        libc = ctypes.CDLL(None, use_errno=True)
        fstatfs = libc.fstatfs
        fstatfs.argtypes = (ctypes.c_int, ctypes.c_void_p)
        fstatfs.restype = ctypes.c_int
        ctypes.set_errno(0)
        if fstatfs(descriptor, ctypes.byref(buffer)) != 0:
            raise OSError(ctypes.get_errno(), "fstatfs")
        filesystem_magic = ctypes.c_long.from_buffer(buffer).value
        after = os.fstat(descriptor)
        named = os.stat(str(resolved), follow_symlinks=False)
        identity = lambda value: (
            value.st_dev, value.st_ino, value.st_mode, value.st_uid,
            value.st_gid, value.st_nlink)
        if identity(before) != identity(after) or identity(after) != identity(named):
            raise TimingError("%s inode identity changed" % context)
        if filesystem_magic != TMPFS_MAGIC:
            raise TimingError("%s must reside on tmpfs" % context)
        record = {
            "schema": "wirehair.wh2.tmpfs_binding.v1",
            "path": str(resolved), "device": before.st_dev,
            "inode": before.st_ino, "mode": stat.S_IMODE(before.st_mode),
            "uid": before.st_uid, "gid": before.st_gid,
            "nlink": before.st_nlink,
            "kind": expected_kind, "filesystem_magic": filesystem_magic,
        }
        return _validate_tmpfs_binding_record(
            record, resolved, expected_kind, context)
    except OSError as exc:
        raise TimingError("cannot bind %s to tmpfs" % context) from exc
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _prepared_campaign_binding_paths(
    root: Path, design: Mapping[str, object],
) -> Tuple[Tuple[Path, str], ...]:
    immutable = design.get("immutable_files")
    if not isinstance(immutable, dict):
        raise TimingError("prepared tmpfs inventory ledger is malformed")
    relative_files = {
        "design.json", "prepare_receipt.json", "tasks_manifest.jsonl",
        *immutable.keys(),
    }
    if any(
            not isinstance(relative, str) or not relative or
            Path(relative).is_absolute() or ".." in Path(relative).parts or
            str(Path(relative)) != relative
            for relative in relative_files):
        raise TimingError("prepared tmpfs file path is malformed")
    entries = [(root, "directory")]
    entries.extend(
        (root / relative, "directory")
        for relative in PREPARED_CAMPAIGN_DIRECTORIES)
    entries.extend(
        (root / relative, "regular") for relative in sorted(relative_files))
    return tuple(entries)


def capture_prepared_tree_tmpfs_bindings(
    root: Path, design: Mapping[str, object],
) -> List[Dict[str, object]]:
    """Capture every prepared input and writable output-directory mount."""
    bindings = [
        capture_tmpfs_binding(path, kind, "prepared campaign tree")
        for path, kind in _prepared_campaign_binding_paths(root, design)
    ]
    if (not bindings or
            any(binding["device"] != bindings[0]["device"]
                for binding in bindings[1:])):
        raise TimingError("prepared campaign tree spans filesystems")
    return bindings


def validate_prepared_tree_tmpfs_bindings(
    records: object, root: Path, design: Mapping[str, object],
) -> List[Dict[str, object]]:
    expected = _prepared_campaign_binding_paths(root, design)
    if not isinstance(records, list) or len(records) != len(expected):
        raise TimingError("prepared campaign tmpfs ledger is malformed")
    checked = [
        validate_tmpfs_binding(record, path, kind, "prepared campaign tree")
        for record, (path, kind) in zip(records, expected)
    ]
    if (not checked or
            any(binding["device"] != checked[0]["device"]
                for binding in checked[1:])):
        raise TimingError("prepared campaign tmpfs ledger spans filesystems")
    return checked


def require_live_tmpfs_tree(
    root: Path, expected_device: int, context: str,
) -> None:
    """Reject nested mounts and non-regular objects anywhere in a result tree."""
    if (not isinstance(expected_device, int) or
            isinstance(expected_device, bool) or expected_device <= 0):
        raise TimingError("%s tmpfs device is malformed" % context)
    paths = [root, *sorted(root.rglob("*"), key=lambda path: str(path))]
    if len(paths) > MAX_CPU_LIST_CARDINALITY:
        raise TimingError("%s tmpfs tree exceeds its entry bound" % context)
    for path in paths:
        if path.is_symlink():
            raise TimingError("%s tmpfs tree contains a symlink" % context)
        if path.is_dir():
            kind = "directory"
        elif path.is_file():
            kind = "regular"
        else:
            raise TimingError("%s tmpfs tree contains a special file" % context)
        binding = capture_tmpfs_binding(
            path.resolve(), kind, "%s tmpfs tree" % context)
        if binding["device"] != expected_device:
            raise TimingError("%s tmpfs tree spans filesystems" % context)


def _parse_guarded_irq_rows(
    raw: bytes, expected_identities: Mapping[int, str],
) -> Dict[int, Dict[str, object]]:
    """Parse exact identities and global counts for guarded numeric IRQs."""
    if (not isinstance(raw, bytes) or not raw.endswith(b"\n") or
            len(raw) > MAX_JSON_EVIDENCE_BYTES or b"\r" in raw or b"\0" in raw or
            not isinstance(expected_identities, Mapping) or
            not expected_identities or
            any(not isinstance(irq, int) or isinstance(irq, bool) or
                not 0 <= irq <= MAX_CPU_ID or
                not isinstance(identity, str) or not identity
                for irq, identity in expected_identities.items())):
        raise TimingError("guarded interrupt evidence is malformed")
    lines = raw.splitlines()
    if not lines:
        raise TimingError("guarded interrupt evidence is empty")
    header = lines[0].split()
    if (not header or len(header) > MAX_CPU_LIST_CARDINALITY or
            any(re.fullmatch(rb"CPU(0|[1-9][0-9]*)", token) is None
                for token in header)):
        raise TimingError("guarded interrupt CPU header is malformed")
    cpu_ids = [parse_uint(
        token[3:].decode("ascii"), "guarded interrupt CPU", MAX_CPU_ID)
        for token in header]
    if cpu_ids != sorted(set(cpu_ids)):
        raise TimingError("guarded interrupt CPU header is unordered")
    rows: Dict[int, List[bytes]] = {}
    for line in lines[1:]:
        match = re.fullmatch(rb"[ \t]*(0|[1-9][0-9]*):[ \t]*(.*)", line)
        if match is None:
            continue
        irq = parse_uint(
            match.group(1).decode("ascii"), "guarded IRQ", MAX_CPU_ID)
        if irq in expected_identities:
            if irq in rows:
                raise TimingError("guarded interrupt row is ambiguous")
            rows[irq] = match.group(2).split()
    if set(rows) != set(expected_identities):
        raise TimingError("guarded interrupt row is missing")
    result: Dict[int, Dict[str, object]] = {}
    for irq in sorted(rows):
        tokens = rows[irq]
        if len(tokens) <= len(header):
            raise TimingError("guarded interrupt row is malformed")
        try:
            counters = [parse_uint(
                token.decode("ascii", "strict"),
                "guarded interrupt counter", MAX_UNSIGNED_COUNTER)
                for token in tokens[:len(header)]]
            identity = b" ".join(tokens[len(header):]).decode("ascii", "strict")
        except UnicodeDecodeError as exc:
            raise TimingError("guarded interrupt row is not ASCII") from exc
        total = sum(counters)
        if total > MAX_UNSIGNED_COUNTER:
            raise TimingError("guarded interrupt count exceeds its domain")
        if identity != expected_identities[irq]:
            raise TimingError("guarded interrupt identity changed")
        result[irq] = {
            "identity": identity, "total_count": total,
            "cpu_ids": cpu_ids, "counters": counters,
        }
    return result


def _target_cpu_tuple(target_cpus: Sequence[int], context: str) -> Tuple[int, ...]:
    try:
        target = tuple(target_cpus)
    except TypeError as exc:
        raise TimingError("%s target CPUs are malformed" % context) from exc
    if (len(target) != 3 or len(set(target)) != 3 or
            any(not isinstance(cpu, int) or isinstance(cpu, bool) or
                not 0 <= cpu <= MAX_CPU_ID for cpu in target)):
        raise TimingError("%s target CPUs are malformed" % context)
    return target


def _parse_target_counter_rows(
    raw: bytes, target_cpus: Sequence[int], source: str,
) -> Tuple[List[int], List[List[object]]]:
    target = _target_cpu_tuple(target_cpus, source)
    if (not isinstance(raw, bytes) or not raw.endswith(b"\n") or
            len(raw) > MAX_JSON_EVIDENCE_BYTES or b"\r" in raw or b"\0" in raw or
            source not in ("hardirq", "softirq")):
        raise TimingError("%s counter source is malformed" % source)
    lines = raw.splitlines()
    if len(lines) < 2:
        raise TimingError("%s counter source is empty" % source)
    header = lines[0].split()
    if (not header or len(header) > MAX_CPU_LIST_CARDINALITY or
            any(re.fullmatch(rb"CPU(0|[1-9][0-9]*)", token) is None
                for token in header)):
        raise TimingError("%s CPU header is malformed" % source)
    cpu_ids = [parse_uint(
        token[3:].decode("ascii"), "%s CPU" % source, MAX_CPU_ID)
        for token in header]
    if cpu_ids != sorted(set(cpu_ids)) or not set(target) <= set(cpu_ids):
        raise TimingError("%s CPU topology changed" % source)
    target_columns = [cpu_ids.index(cpu) for cpu in target]
    rows: List[List[object]] = []
    seen = set()
    numeric_phase = True
    previous_numeric = -1
    for raw_line in lines[1:]:
        match = re.fullmatch(rb"[ \t]*([^:\s]+):[ \t]*(.*)", raw_line)
        if match is None:
            raise TimingError("%s counter row is malformed" % source)
        try:
            vector = match.group(1).decode("ascii", "strict")
        except UnicodeDecodeError as exc:
            raise TimingError("%s vector is not ASCII" % source) from exc
        numeric = UINT_RE.fullmatch(vector) is not None
        if numeric:
            number = parse_uint(vector, "hardirq vector", MAX_CPU_ID)
            if not numeric_phase or number <= previous_numeric:
                raise TimingError("hardirq numeric topology is unordered")
            previous_numeric = number
            kind = "numeric"
        else:
            numeric_phase = False
            if re.fullmatch(r"[A-Z][A-Z0-9_]*", vector) is None:
                raise TimingError("%s named vector is malformed" % source)
            kind = "named"
        if vector in seen:
            raise TimingError("%s vector is duplicated" % source)
        seen.add(vector)
        tokens = match.group(2).split()
        if source == "hardirq" and vector in GLOBAL_HARDIRQ_VECTORS:
            if len(tokens) != 1:
                raise TimingError("global hardirq row width changed")
            try:
                counter = parse_uint(
                    tokens[0].decode("ascii", "strict"),
                    "global hardirq counter", MAX_UNSIGNED_COUNTER)
            except UnicodeDecodeError as exc:
                raise TimingError("hardirq counter is not ASCII") from exc
            rows.append([vector, "global", counter])
            continue
        if len(tokens) < len(header) or (
                source == "softirq" and len(tokens) != len(header)):
            raise TimingError("%s counter row width changed" % source)
        counters = []
        try:
            for token in tokens[:len(header)]:
                counters.append(parse_uint(
                    token.decode("ascii", "strict"),
                    "%s counter" % source, MAX_UNSIGNED_COUNTER))
        except UnicodeDecodeError as exc:
            raise TimingError("%s counter is not ASCII" % source) from exc
        selected = [counters[index] for index in target_columns]
        if source == "hardirq":
            if not numeric and vector in GLOBAL_HARDIRQ_VECTORS:
                raise TimingError("global hardirq kind changed")
            suffix = tokens[len(header):]
            if not suffix:
                raise TimingError("hardirq identity is missing")
            try:
                identity = b" ".join(suffix).decode("ascii", "strict")
            except UnicodeDecodeError as exc:
                raise TimingError("hardirq identity is not ASCII") from exc
            rows.append([vector, kind, identity, *selected])
        else:
            if numeric:
                raise TimingError("softirq vector must be named")
            rows.append([vector, *selected])
    if not rows:
        raise TimingError("%s counter rows are empty" % source)
    if source == "hardirq" and not GLOBAL_HARDIRQ_VECTORS <= seen:
        raise TimingError("global hardirq rows are missing")
    if source == "softirq" and [row[0] for row in rows] != \
            list(EXPECTED_SOFTIRQ_VECTORS):
        raise TimingError("softirq vector topology changed")
    return cpu_ids, rows


def _target_rows_sha256(
    source: str, target_cpus: Sequence[int], cpu_ids: Sequence[int],
    rows: Sequence[object],
) -> str:
    return sha256_bytes(canonical_json({
        "schema": "wirehair.wh2.%s_target_columns.v2" % source,
        "target_cpus": list(target_cpus), "cpu_ids": list(cpu_ids),
        "rows": list(rows),
    }))


def validate_target_irq_snapshot(
    snapshot: object, target_cpus: Sequence[int], context: str,
) -> Dict[str, object]:
    target = _target_cpu_tuple(target_cpus, context)
    if (not isinstance(snapshot, dict) or
            set(snapshot) != TARGET_IRQ_SNAPSHOT_FIELDS or
            snapshot.get("schema") != "wirehair.wh2.target_irq_snapshot.v2" or
            snapshot.get("target_cpus") != list(target)):
        raise TimingError("%s IRQ snapshot is malformed" % context)
    checked_attempt_duration_ns(
        snapshot.get("capture_start_monotonic_ns"),
        snapshot.get("capture_end_monotonic_ns"),
        snapshot.get("capture_duration_ns"), "%s IRQ capture" % context)
    for field in ("hardirq_sha256", "softirq_sha256"):
        require_sha256(snapshot.get(field), "%s %s" % (context, field))
    cpu_ids = snapshot.get("cpu_ids")
    if (not isinstance(cpu_ids, list) or not cpu_ids or
            len(cpu_ids) > MAX_CPU_LIST_CARDINALITY or
            any(not isinstance(cpu, int) or isinstance(cpu, bool) or
                not 0 <= cpu <= MAX_CPU_ID for cpu in cpu_ids) or
            cpu_ids != sorted(set(cpu_ids)) or not set(target) <= set(cpu_ids)):
        raise TimingError("%s IRQ CPU topology is malformed" % context)
    hardirq_rows = snapshot.get("hardirq_rows")
    softirq_rows = snapshot.get("softirq_rows")
    if (not isinstance(hardirq_rows, list) or not hardirq_rows or
            len(hardirq_rows) > MAX_CPU_LIST_CARDINALITY or
            not isinstance(softirq_rows, list) or not softirq_rows or
            len(softirq_rows) > MAX_CPU_LIST_CARDINALITY):
        raise TimingError("%s IRQ snapshot rows are malformed" % context)
    previous_numeric = -1
    numeric_phase = True
    hard_seen = set()
    for row in hardirq_rows:
        if not isinstance(row, list) or len(row) < 2:
            raise TimingError("%s hardirq row is malformed" % context)
        kind = row[1]
        if (((kind in ("numeric", "named") and
              len(row) != 3 + len(target)) or
             (kind == "global" and len(row) != 3) or
             kind not in ("numeric", "named", "global")) or
                not isinstance(row[0], str) or
                (kind in ("numeric", "named") and
                 (not isinstance(row[2], str) or not row[2] or
                  not row[2].isascii() or
                  " ".join(row[2].split()) != row[2])) or
                any(not isinstance(value, int) or isinstance(value, bool) or
                    not 0 <= value <= MAX_UNSIGNED_COUNTER
                    for value in row[2 if kind == "global" else 3:])):
            raise TimingError("%s hardirq row is malformed" % context)
        vector = row[0]
        numeric = UINT_RE.fullmatch(vector) is not None
        if numeric:
            number = parse_uint(
                vector, "%s hardirq vector" % context, MAX_CPU_ID)
        if ((numeric and kind != "numeric") or
                (not numeric and kind not in ("named", "global")) or
                (vector in GLOBAL_HARDIRQ_VECTORS and kind != "global") or
                (kind == "global" and vector not in GLOBAL_HARDIRQ_VECTORS) or
                vector in hard_seen):
            raise TimingError("%s hardirq identity is malformed" % context)
        hard_seen.add(vector)
        if numeric:
            if not numeric_phase or number <= previous_numeric:
                raise TimingError("%s hardirq order is malformed" % context)
            previous_numeric = number
        else:
            numeric_phase = False
            if re.fullmatch(r"[A-Z][A-Z0-9_]*", vector) is None:
                raise TimingError("%s hardirq name is malformed" % context)
    if not GLOBAL_HARDIRQ_VECTORS <= hard_seen:
        raise TimingError("%s global hardirq rows are missing" % context)
    soft_seen = set()
    for row in softirq_rows:
        if (not isinstance(row, list) or len(row) != 1 + len(target) or
                not isinstance(row[0], str) or
                re.fullmatch(r"[A-Z][A-Z0-9_]*", row[0]) is None or
                row[0] in soft_seen or
                any(not isinstance(value, int) or isinstance(value, bool) or
                    not 0 <= value <= MAX_UNSIGNED_COUNTER for value in row[1:])):
            raise TimingError("%s softirq row is malformed" % context)
        soft_seen.add(row[0])
    if [row[0] for row in softirq_rows] != list(EXPECTED_SOFTIRQ_VECTORS):
        raise TimingError("%s softirq topology is malformed" % context)
    if (snapshot["hardirq_sha256"] != _target_rows_sha256(
            "hardirq", target, cpu_ids, hardirq_rows) or
            snapshot["softirq_sha256"] != _target_rows_sha256(
                "softirq", target, cpu_ids, softirq_rows)):
        raise TimingError("%s IRQ target-column hash mismatch" % context)
    return snapshot


def parse_target_irq_snapshot(
    interrupts_raw: bytes, softirqs_raw: bytes,
    target_cpus: Sequence[int], capture_start_ns: int, capture_end_ns: int,
) -> Dict[str, object]:
    target = _target_cpu_tuple(target_cpus, "captured IRQ")
    hardirq_cpu_ids, hardirq_rows = _parse_target_counter_rows(
        interrupts_raw, target, "hardirq")
    softirq_cpu_ids, softirq_rows = _parse_target_counter_rows(
        softirqs_raw, target, "softirq")
    if hardirq_cpu_ids != softirq_cpu_ids:
        raise TimingError("hardirq/softirq CPU topology changed during capture")
    duration_ns = checked_attempt_duration_ns(
        capture_start_ns, capture_end_ns, capture_end_ns - capture_start_ns,
        "captured IRQ")
    snapshot = {
        "schema": "wirehair.wh2.target_irq_snapshot.v2",
        "capture_start_monotonic_ns": capture_start_ns,
        "capture_end_monotonic_ns": capture_end_ns,
        "capture_duration_ns": duration_ns,
        "target_cpus": list(target),
        "cpu_ids": hardirq_cpu_ids,
        "hardirq_rows": hardirq_rows,
        "hardirq_sha256": _target_rows_sha256(
            "hardirq", target, hardirq_cpu_ids, hardirq_rows),
        "softirq_rows": softirq_rows,
        "softirq_sha256": _target_rows_sha256(
            "softirq", target, softirq_cpu_ids, softirq_rows),
    }
    return validate_target_irq_snapshot(snapshot, target, "captured IRQ")


def capture_target_irq_snapshot(
    target_cpus: Sequence[int], *, proc_root: Path = Path("/proc"),
) -> Dict[str, object]:
    capture_start_ns = time.monotonic_ns()
    interrupts_raw = _bounded_kernel_bytes(proc_root / "interrupts")
    softirqs_raw = _bounded_kernel_bytes(proc_root / "softirqs")
    capture_end_ns = time.monotonic_ns()
    return parse_target_irq_snapshot(
        interrupts_raw, softirqs_raw, target_cpus,
        capture_start_ns, capture_end_ns)


def checked_target_irq_delta(
    before: object, after: object, target_cpus: Sequence[int], context: str,
) -> Dict[str, object]:
    target = _target_cpu_tuple(target_cpus, context)
    before_record = validate_target_irq_snapshot(
        before, target, "%s before" % context)
    after_record = validate_target_irq_snapshot(
        after, target, "%s after" % context)
    hard_before = before_record["hardirq_rows"]
    hard_after = after_record["hardirq_rows"]
    soft_before = before_record["softirq_rows"]
    soft_after = after_record["softirq_rows"]
    if (before_record["capture_end_monotonic_ns"] >
            after_record["capture_start_monotonic_ns"] or
            before_record["cpu_ids"] != after_record["cpu_ids"] or
            [row[:2] if row[1] == "global" else row[:3]
             for row in hard_before] !=
            [row[:2] if row[1] == "global" else row[:3]
             for row in hard_after] or
            [row[:1] for row in soft_before] != [row[:1] for row in soft_after]):
        raise TimingError("%s IRQ vector topology changed" % context)
    hardirq_deltas = []
    softirq_deltas = []
    classifications = []
    contaminations = []
    for left, right in zip(hard_before, hard_after):
        if left[1] == "global":
            if right[2] < left[2]:
                raise TimingError("%s hardirq counter reset" % context)
            delta = right[2] - left[2]
            hardirq_deltas.append([left[0], "global", delta])
            if delta:
                classification = "global-hardirq:%s:%d" % (left[0], delta)
                classifications.append(classification)
                contaminations.append(classification)
            continue
        deltas = []
        for index, cpu in enumerate(target):
            before_count = left[index + 3]
            after_count = right[index + 3]
            if after_count < before_count:
                raise TimingError("%s hardirq counter reset" % context)
            delta = after_count - before_count
            deltas.append(delta)
            if not delta:
                continue
            classification = "%s-hardirq:%s:cpu%d:%d" % (
                left[1], left[0], cpu, delta)
            classifications.append(classification)
            if left[1] == "numeric":
                contaminations.append(classification)
        hardirq_deltas.append([left[0], left[1], *deltas])
    for left, right in zip(soft_before, soft_after):
        deltas = []
        for index, cpu in enumerate(target):
            before_count = left[index + 1]
            after_count = right[index + 1]
            if after_count < before_count:
                raise TimingError("%s softirq counter reset" % context)
            delta = after_count - before_count
            deltas.append(delta)
            if not delta:
                continue
            classification = "softirq:%s:cpu%d:%d" % (left[0], cpu, delta)
            classifications.append(classification)
            if index < 2 and left[0] in DEVICE_SOFTIRQ_VECTORS:
                contaminations.append(classification)
        softirq_deltas.append([left[0], *deltas])
    return {
        "schema": "wirehair.wh2.target_irq_delta.v2",
        "target_cpus": list(target),
        "hardirq_before_sha256": before_record["hardirq_sha256"],
        "hardirq_after_sha256": after_record["hardirq_sha256"],
        "softirq_before_sha256": before_record["softirq_sha256"],
        "softirq_after_sha256": after_record["softirq_sha256"],
        "hardirq_deltas": hardirq_deltas,
        "softirq_deltas": softirq_deltas,
        "classifications": classifications,
        "contaminations": contaminations,
    }


def validate_target_irq_delta(
    delta: object, before: object, after: object,
    target_cpus: Sequence[int], context: str,
) -> Dict[str, object]:
    if not isinstance(delta, dict) or set(delta) != TARGET_IRQ_DELTA_FIELDS:
        raise TimingError("%s IRQ delta is malformed" % context)
    expected = checked_target_irq_delta(before, after, target_cpus, context)
    if delta != expected:
        raise TimingError("%s IRQ delta does not replay" % context)
    return delta


def require_target_irq_contained_interval(
    before: Mapping[str, object], after: Mapping[str, object],
    start_ns: int, end_ns: int, context: str,
) -> None:
    """Require two ordered captures to live wholly inside one interval."""
    if (not isinstance(start_ns, int) or isinstance(start_ns, bool) or
            not isinstance(end_ns, int) or isinstance(end_ns, bool) or
            not 0 <= start_ns < end_ns <= MAX_UNSIGNED_COUNTER or
            before.get("capture_start_monotonic_ns") < start_ns or
            after.get("capture_end_monotonic_ns") > end_ns):
        raise TimingError("%s IRQ capture interval is malformed" % context)


def require_target_irq_bracketing_interval(
    before: Mapping[str, object], after: Mapping[str, object],
    start_ns: int, end_ns: int, context: str,
) -> None:
    """Require captures to end before and start after the measured interval."""
    if (not isinstance(start_ns, int) or isinstance(start_ns, bool) or
            not isinstance(end_ns, int) or isinstance(end_ns, bool) or
            not 0 <= start_ns < end_ns <= MAX_UNSIGNED_COUNTER or
            before.get("capture_end_monotonic_ns") > start_ns or
            after.get("capture_start_monotonic_ns") < end_ns):
        raise TimingError("%s IRQ capture bracket is malformed" % context)


def _snapshot_cpu_list(
    value: object, parsed: object, context: str, *, allow_empty: bool = False,
) -> Tuple[int, ...]:
    if not isinstance(value, str):
        raise TimingError("%s CPU-list receipt is malformed" % context)
    if not value and allow_empty:
        cpus: Tuple[int, ...] = ()
    else:
        cpus = parse_cpu_list(value)
    if (not isinstance(parsed, list) or
            any(not isinstance(cpu, int) or isinstance(cpu, bool)
                for cpu in parsed) or list(cpus) != parsed):
        raise TimingError("%s parsed CPU-list receipt changed" % context)
    return cpus


def validate_runtime_isolation_snapshot(
    snapshot: object, expected_cpus: Sequence[int], context: str,
) -> Dict[str, object]:
    try:
        expected = tuple(sorted(expected_cpus))
    except (TypeError, ValueError) as exc:
        raise TimingError("%s expected isolated CPUs are malformed" % context) \
            from exc
    if (not expected or len(expected) != len(set(expected)) or
            any(not isinstance(cpu, int) or isinstance(cpu, bool) or
                not 0 <= cpu <= MAX_CPU_ID for cpu in expected)):
        raise TimingError("%s expected isolated CPUs are malformed" % context)
    if (not isinstance(snapshot, dict) or
            set(snapshot) != RUNTIME_ISOLATION_SNAPSHOT_FIELDS or
            snapshot.get("schema") !=
            "wirehair.wh2.runtime_isolation_snapshot.v2" or
            snapshot.get("self_cgroup") != "/wh2-timing-v4" or
            snapshot.get("expected_isolated_cpus") != list(expected) or
            snapshot.get("cgroup_partition") != "isolated"):
        raise TimingError("%s isolation snapshot is malformed" % context)
    checked_attempt_duration_ns(
        snapshot.get("capture_start_monotonic_ns"),
        snapshot.get("capture_end_monotonic_ns"),
        snapshot.get("capture_duration_ns"), "%s isolation capture" % context)
    for text_field, parsed_field in (
            ("kernel_isolated_cpu_list", "kernel_isolated_cpus"),
            ("cgroup_cpu_list", "cgroup_cpus"),
            ("cgroup_effective_cpu_list", "cgroup_effective_cpus"),
            ("cgroup_exclusive_cpu_list", "cgroup_exclusive_cpus"),
            ("cgroup_exclusive_effective_cpu_list",
             "cgroup_exclusive_effective_cpus")):
        if _snapshot_cpu_list(
                snapshot.get(text_field), snapshot.get(parsed_field),
                "%s %s" % (context, text_field)) != expected:
            raise TimingError("%s isolated CPU set changed" % context)
    affinities = snapshot.get("irq_effective_affinities")
    if not isinstance(affinities, list) or not affinities:
        raise TimingError("%s IRQ affinity inventory is malformed" % context)
    previous_irq = -1
    irq30_record: Optional[Dict[str, object]] = None
    intersecting = []
    affinity_by_irq = {}
    for record in affinities:
        if (not isinstance(record, dict) or
                set(record) != IRQ_EFFECTIVE_AFFINITY_FIELDS):
            raise TimingError("%s IRQ affinity record is malformed" % context)
        irq = record.get("irq")
        if (not isinstance(irq, int) or isinstance(irq, bool) or
                not previous_irq < irq <= MAX_CPU_ID):
            raise TimingError("%s IRQ affinity order is malformed" % context)
        previous_irq = irq
        effective = _snapshot_cpu_list(
            record.get("effective_affinity_list"),
            record.get("effective_cpus"),
            "%s IRQ %d effective" % (context, irq), allow_empty=True)
        overlap = tuple(sorted(set(effective).intersection(expected)))
        if overlap:
            intersecting.append(irq)
        affinity_by_irq[irq] = record
        if irq == 30:
            irq30_record = record
    exception = snapshot.get("irq30_exception")
    managed = snapshot.get("managed_nvme_exceptions")
    expected_intersecting = [30] + [
        irq for irq, _identity, _handler, _requested, _effective
        in MANAGED_NVME_IRQ_WHITELIST]
    if intersecting != sorted(expected_intersecting) or exception is None or \
            not isinstance(managed, list) or \
            len(managed) != len(MANAGED_NVME_IRQ_WHITELIST):
        raise TimingError("%s numeric IRQ reaches an isolated CPU" % context)
    if (not isinstance(exception, dict) or
            set(exception) != IRQ30_EXCEPTION_FIELDS or
            exception.get("irq") != 30 or
            exception.get("identity") != IRQ30_IDENTITY or
            exception.get("handler_directories") != ["AMD-Vi0-PPR"] or
            not isinstance(exception.get("global_interrupt_count"), int) or
            isinstance(exception.get("global_interrupt_count"), bool) or
            exception.get("global_interrupt_count") != 0 or
            irq30_record is None):
        raise TimingError("%s guarded IRQ30 identity/count is malformed" % context)
    requested = _snapshot_cpu_list(
        exception.get("requested_affinity_list"),
        exception.get("requested_cpus"), "%s IRQ30 requested" % context)
    effective = _snapshot_cpu_list(
        exception.get("effective_affinity_list"),
        exception.get("effective_cpus"), "%s IRQ30 effective" % context)
    if (not requested or
            exception.get("requested_affinity_list") != "0-127" or
            list(requested) != list(range(128)) or
            exception.get("effective_affinity_list") != "8" or
            list(effective) != [8] or 8 not in expected or
            exception.get("effective_affinity_list") !=
            irq30_record.get("effective_affinity_list") or
            list(effective) != irq30_record.get("effective_cpus") or
            not set(effective).intersection(expected)):
        raise TimingError("%s guarded IRQ30 affinity is malformed" % context)
    for record, frozen in zip(managed, MANAGED_NVME_IRQ_WHITELIST):
        irq, identity, handler, requested_text, effective_text = frozen
        if (not isinstance(record, dict) or
                set(record) != MANAGED_IRQ_EXCEPTION_FIELDS or
                record.get("irq") != irq or
                record.get("identity") != identity or
                record.get("handler_directories") != [handler] or
                record.get("requested_affinity_list") != requested_text or
                record.get("effective_affinity_list") != effective_text or
                irq not in affinity_by_irq):
            raise TimingError(
                "%s managed NVMe IRQ identity changed" % context)
        requested_cpus = _snapshot_cpu_list(
            record.get("requested_affinity_list"),
            record.get("requested_cpus"),
            "%s managed IRQ %d requested" % (context, irq))
        effective_cpus = _snapshot_cpu_list(
            record.get("effective_affinity_list"),
            record.get("effective_cpus"),
            "%s managed IRQ %d effective" % (context, irq))
        frozen_effective = affinity_by_irq[irq]
        if (not requested_cpus or not set(effective_cpus).intersection(expected) or
                record.get("effective_affinity_list") !=
                frozen_effective.get("effective_affinity_list") or
                list(effective_cpus) != frozen_effective.get("effective_cpus")):
            raise TimingError(
                "%s managed NVMe IRQ affinity changed" % context)
    return snapshot


def capture_runtime_isolation_snapshot(
    core: int, sibling: int, controller: int, *,
    proc_root: Path = Path("/proc"),
    cgroup_root: Path = Path("/sys/fs/cgroup"),
) -> Dict[str, object]:
    values = (core, sibling, controller)
    if (any(not isinstance(cpu, int) or isinstance(cpu, bool) or
            not 0 <= cpu <= MAX_CPU_ID for cpu in values) or
            len(set(values)) != 3):
        raise TimingError("runtime isolation CPUs are malformed or not distinct")
    expected = tuple(sorted(values))
    capture_start_ns = time.monotonic_ns()
    self_cgroup = _kernel_ascii_line(
        proc_root / "self/cgroup", "self cgroup membership")
    if self_cgroup != "0::/wh2-timing-v4":
        raise TimingError("campaign is not in cgroup /wh2-timing-v4")
    group = cgroup_root / "wh2-timing-v4"

    def cpu_list_record(
        path: Path, context: str, *, allow_empty: bool = False,
    ) -> Tuple[str, List[int]]:
        text = _kernel_ascii_line(path, context, allow_empty=allow_empty)
        return text, list(parse_cpu_list(text)) if text else []

    kernel_text, kernel_cpus = cpu_list_record(
        cgroup_root / "cpuset.cpus.isolated", "kernel isolated CPUs")
    cpu_text, cpus = cpu_list_record(group / "cpuset.cpus", "cgroup CPUs")
    effective_text, effective = cpu_list_record(
        group / "cpuset.cpus.effective", "cgroup effective CPUs")
    exclusive_text, exclusive = cpu_list_record(
        group / "cpuset.cpus.exclusive", "cgroup exclusive CPUs")
    exclusive_effective_text, exclusive_effective = cpu_list_record(
        group / "cpuset.cpus.exclusive.effective",
        "cgroup exclusive effective CPUs")
    partition = _kernel_ascii_line(
        group / "cpuset.cpus.partition", "cgroup partition")
    guarded_before = _parse_guarded_irq_rows(
        _bounded_kernel_bytes(proc_root / "interrupts"),
        GUARDED_IRQ_IDENTITIES)
    irq_root = proc_root / "irq"
    try:
        entries = list(irq_root.iterdir())
    except OSError as exc:
        raise TimingError("cannot enumerate numeric IRQs") from exc
    irq_paths = []
    for path in entries:
        if UINT_RE.fullmatch(path.name) is None:
            continue
        if path.is_symlink() or not path.is_dir():
            raise TimingError("numeric IRQ path is not a plain directory")
        irq_paths.append((
            parse_uint(path.name, "numeric IRQ", MAX_CPU_ID), path))
    irq_paths.sort(key=lambda item: item[0])
    if not irq_paths or len({irq for irq, _path in irq_paths}) != len(irq_paths):
        raise TimingError("numeric IRQ inventory is empty or duplicated")
    records = []
    exception: Optional[Dict[str, object]] = None
    managed_exceptions = []
    managed_by_irq = {
        irq: (identity, handler, requested, effective)
        for irq, identity, handler, requested, effective
        in MANAGED_NVME_IRQ_WHITELIST}
    for irq, path in irq_paths:
        affinity_text, affinity_cpus = cpu_list_record(
            path / "effective_affinity_list",
            "IRQ %d effective affinity" % irq, allow_empty=True)
        record = {
            "irq": irq, "effective_affinity_list": affinity_text,
            "effective_cpus": affinity_cpus,
        }
        records.append(record)
        if set(affinity_cpus).intersection(expected):
            if irq != 30 and irq not in managed_by_irq:
                raise TimingError(
                    "numeric IRQ %d reaches an isolated timing CPU" % irq)
            try:
                handler_directories = sorted(
                    child.name for child in path.iterdir()
                    if child.is_dir() and not child.is_symlink())
            except OSError as exc:
                raise TimingError(
                    "cannot enumerate IRQ %d handlers" % irq) from exc
            requested_text, requested_cpus = cpu_list_record(
                path / "smp_affinity_list",
                "IRQ %d requested affinity" % irq)
            if irq == 30:
                exception = {
                    "irq": 30, "identity": guarded_before[30]["identity"],
                    "handler_directories": handler_directories,
                    "requested_affinity_list": requested_text,
                    "requested_cpus": requested_cpus,
                    "effective_affinity_list": affinity_text,
                    "effective_cpus": affinity_cpus,
                    "global_interrupt_count":
                        guarded_before[30]["total_count"],
                }
            else:
                managed_exceptions.append({
                    "irq": irq,
                    "identity": guarded_before[irq]["identity"],
                    "handler_directories": handler_directories,
                    "requested_affinity_list": requested_text,
                    "requested_cpus": requested_cpus,
                    "effective_affinity_list": affinity_text,
                    "effective_cpus": affinity_cpus,
                })
    guarded_after = _parse_guarded_irq_rows(
        _bounded_kernel_bytes(proc_root / "interrupts"),
        GUARDED_IRQ_IDENTITIES)
    for irq in sorted(GUARDED_IRQ_IDENTITIES):
        before_guard = guarded_before[irq]
        after_guard = guarded_after[irq]
        if (before_guard["identity"] != after_guard["identity"] or
                before_guard["cpu_ids"] != after_guard["cpu_ids"] or
                len(before_guard["counters"]) != len(after_guard["counters"]) or
                any(right < left for left, right in zip(
                    before_guard["counters"], after_guard["counters"]))):
            raise TimingError("guarded IRQ topology/counter changed during capture")
    if (guarded_before[30]["total_count"] != 0 or
            guarded_after[30]["total_count"] != 0):
        raise TimingError("guarded IRQ30 fired during isolation capture")
    capture_end_ns = time.monotonic_ns()
    snapshot = {
        "schema": "wirehair.wh2.runtime_isolation_snapshot.v2",
        "capture_start_monotonic_ns": capture_start_ns,
        "capture_end_monotonic_ns": capture_end_ns,
        "capture_duration_ns": capture_end_ns - capture_start_ns,
        "self_cgroup": self_cgroup.split("::", 1)[1],
        "expected_isolated_cpus": list(expected),
        "kernel_isolated_cpu_list": kernel_text,
        "kernel_isolated_cpus": kernel_cpus,
        "cgroup_cpu_list": cpu_text, "cgroup_cpus": cpus,
        "cgroup_effective_cpu_list": effective_text,
        "cgroup_effective_cpus": effective,
        "cgroup_exclusive_cpu_list": exclusive_text,
        "cgroup_exclusive_cpus": exclusive,
        "cgroup_exclusive_effective_cpu_list": exclusive_effective_text,
        "cgroup_exclusive_effective_cpus": exclusive_effective,
        "cgroup_partition": partition,
        "irq_effective_affinities": records,
        "irq30_exception": exception,
        "managed_nvme_exceptions": managed_exceptions,
    }
    return validate_runtime_isolation_snapshot(
        snapshot, expected, "captured runtime")


def validate_runtime_isolation_transition(
    start: object, end: object, expected_cpus: Sequence[int],
) -> None:
    start_record = validate_runtime_isolation_snapshot(
        start, expected_cpus, "start runtime")
    end_record = validate_runtime_isolation_snapshot(
        end, expected_cpus, "end runtime")
    invariant_fields = RUNTIME_ISOLATION_SNAPSHOT_FIELDS - {
        "irq_effective_affinities", "capture_start_monotonic_ns",
        "capture_end_monotonic_ns", "capture_duration_ns"}
    for field in invariant_fields:
        if start_record[field] != end_record[field]:
            raise TimingError("runtime isolation state changed: %s" % field)


def runtime_isolation_snapshot_sha256(snapshot: object) -> str:
    if not isinstance(snapshot, dict):
        raise TimingError("runtime isolation snapshot is malformed")
    return sha256_bytes(canonical_json(snapshot))


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
            for seed_index, (schedule, seed) in enumerate(SCHEDULE_SEEDS):
                for cache_state in CACHE_STATES:
                    task: Dict[str, object] = {
                        "K": K, "bb": width, "seed_index": seed_index,
                        "seed": seed, "schedule": schedule,
                        "cache_state": cache_state,
                    }
                    priority = sha256_bytes(
                        b"wirehair.wh2.grouped-commit.cross-binary.order.v1\0" +
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
    if len(result) != 108:
        raise TimingError("internal timing task count changed")
    return tuple(result)


@dataclass(frozen=True)
class ParsedOutput:
    label: str
    schema: str
    preamble: Mapping[str, str]
    rows: Tuple[Mapping[str, str], ...]
    work_signature: Tuple[str, ...]
    semantic_sha256: str
    stdout_sha256: str
    timed_elapsed_ns: int
    all_elapsed_ns: int
    timed_phase_ns: Mapping[str, int]
    all_phase_ns: Mapping[str, int]
    timed_minor_faults: int
    discard_minor_faults: int
    max_minor_faults: int
    contaminations: Tuple[str, ...]


def _parse_preamble(line: str) -> Dict[str, str]:
    prefix = "# groupedtiming: "
    if not line.startswith(prefix):
        raise TimingError("missing groupedtiming preamble")
    tokens = line[len(prefix):].split(" ")
    if any(token.count("=") != 1 for token in tokens):
        raise TimingError("groupedtiming preamble token is malformed")
    pairs = tuple(tuple(token.split("=", 1)) for token in tokens)
    if tuple(pair[0] for pair in pairs) != GROUPED_PREAMBLE_FIELDS:
        raise TimingError("groupedtiming preamble order/schema mismatch")
    if len({pair[0] for pair in pairs}) != len(pairs):
        raise TimingError("groupedtiming preamble contains a duplicate")
    return dict(pairs)


def _expected_preamble(
    task: Mapping[str, object], evict_bytes: int,
) -> Dict[str, str]:
    K = int(task["K"])
    width = int(task["bb"])
    return {
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
        "control_grouped_hash_seed": "0xb7e15162",
        "control_final_h_a_columns": "12",
        "candidate_period": str(ARCHITECTURE["period"]),
        "candidate_grouped_rows": str(ARCHITECTURE["grouped_rows"]),
        "candidate_buckets": str(ARCHITECTURE["buckets"]),
        "candidate_grouped_hash_seed": "0xb7e15162",
        "candidate_final_h_a_columns": "12", "gf256_rows": "10",
        "gf16_rows": "2", "dense_two_anchor": "1", "mix": "2",
        "payload": "distinct-packet-zero-v1",
        "payload_count": str(K + OVERHEAD),
        "payload_bytes": str((K + OVERHEAD) * width),
        "payload_alignment": "64", "payload_prefaulted": "1",
        "system_build": "outside-timer",
        "tls_reapply": "full-per-slot-outside-timer",
        "allocator_tls_state": "preflight-warmed",
    }


def parse_grouped_output(
    raw: bytes, label: str, task: Mapping[str, object], evict_bytes: int,
    expected_core: int,
) -> ParsedOutput:
    if label not in ("base", "candidate"):
        raise TimingError("unknown binary label")
    expected_schema = "v1"
    if (not raw or len(raw) > MAX_GROUPED_OUTPUT_BYTES or
            not raw.endswith(b"\n") or b"\r" in raw or b"\0" in raw or
            b'"' in raw):
        raise TimingError("groupedtiming output is not canonical LF text")
    try:
        lines = raw.decode("ascii").splitlines()
    except UnicodeDecodeError as exc:
        raise TimingError("groupedtiming output is not ASCII") from exc
    if len(lines) != 34:
        raise TimingError("groupedtiming output does not have 34 lines")
    if any(line.count(",") != len(GROUPED_FIELDS) - 1 for line in lines[1:]):
        raise TimingError("groupedtiming CSV row arity changed")
    preamble = _parse_preamble(lines[0])
    if preamble.get("schema") != expected_schema:
        raise TimingError("binary emitted the wrong groupedtiming schema")
    for key, value in _expected_preamble(task, evict_bytes).items():
        if preamble.get(key) != value:
            raise TimingError(
                "groupedtiming preamble mismatch %s: %r != %r" %
                (key, preamble.get(key), value))
    trace_sha256 = require_sha256(preamble.get("trace_sha256"), "trace hash")
    for arm in ("control", "candidate"):
        parse_uint(preamble.get(arm + "_attempt"), arm + " attempt", 255)
        for suffix, maximum in (("matrix_seed", (1 << 64) - 1),
                                ("peel_seed", (1 << 32) - 1)):
            text = preamble.get(arm + "_" + suffix, "")
            if HEX_RE.fullmatch(text) is None or int(text, 16) > maximum:
                raise TimingError("malformed %s %s receipt" % (arm, suffix))
        parse_sint(preamble.get("preflight_" + arm + "_result"),
                   "preflight result", -(1 << 31), (1 << 31) - 1)
    for suffix in ("attempt", "matrix_seed", "peel_seed"):
        if preamble["control_" + suffix] != preamble["candidate_" + suffix]:
            raise TimingError("identical internal arms selected different graphs")
    if (preamble.get("preflight_control_result") != "0" or
            preamble.get("preflight_candidate_result") != "0" or
            preamble.get("cell_class") != "common-success" or
            preamble.get("common_success") != "1"):
        raise TimingError("timing cell is not a common-success solve")
    reader = csv.DictReader(io.StringIO("\n".join(lines[1:]) + "\n"))
    if tuple(reader.fieldnames or ()) != GROUPED_FIELDS:
        raise TimingError("groupedtiming CSV schema mismatch")
    rows = tuple(dict(row) for row in reader)
    if len(rows) != 32:
        raise TimingError("groupedtiming output does not have 32 rows")
    if any(set(row) != set(GROUPED_FIELDS) or
           any(value is None for value in row.values()) for row in rows):
        raise TimingError("groupedtiming CSV contains a malformed row")
    K = int(task["K"])
    width = int(task["bb"])
    signatures = set()
    contaminations: List[str] = []
    timed_elapsed = 0
    all_elapsed = 0
    timed_phase = {field: 0 for field in PHASE_FIELDS}
    all_phase = {field: 0 for field in PHASE_FIELDS}
    timed_minor = 0
    discard_minor = 0
    max_minor = 0
    for index, row in enumerate(rows):
        cycle, slot = divmod(index, 8)
        arm = "control" if INNER_ORDER[slot] == "A" else "candidate"
        exact = {
            "N": str(K), "bb": str(width), "overhead": str(OVERHEAD),
            "schedule": str(task["schedule"]), "seed": str(task["seed"]),
            "loss": LOSS_TEXT, "cache_state": str(task["cache_state"]),
            "cycle": str(cycle), "slot": str(slot), "arm": arm,
            "period": str(ARCHITECTURE["period"]),
            "grouped_rows": str(ARCHITECTURE["grouped_rows"]),
            "buckets_requested": str(ARCHITECTURE["buckets"]),
            "seed_attempt": preamble[arm + "_attempt"],
            "matrix_seed": preamble[arm + "_matrix_seed"],
            "peel_seed": preamble[arm + "_peel_seed"],
            "preflight_result": "0", "cell_class": "common-success",
            "common_success": "1", "result": "0", "outcome_stable": "1",
            "source_bytes": str(K * width),
            "packet_payload_bytes": str((K + OVERHEAD) * width),
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
            timed_elapsed += elapsed
        for field in PHASE_FIELDS:
            phase_ns = parse_uint(
                row[field], "row %d %s" % (index, field))
            all_phase[field] += phase_ns
            if cycle != 0:
                timed_phase[field] += phase_ns
        if sum(parse_uint(row[field], "row %d %s" % (index, field))
               for field in PHASE_FIELDS) > elapsed:
            raise TimingError("timing phases exceed elapsed time")
        inactivated = parse_uint(row["inactivated"], "inactivated")
        binary_def = parse_uint(row["binary_def"], "binary_def")
        heavy_gain = parse_uint(row["heavy_gain"], "heavy_gain")
        if binary_def > inactivated or heavy_gain != binary_def:
            raise TimingError("residual rank/work relationship changed")
        minflt = parse_sint(
            row.get("minflt_delta"), "minor-fault delta",
            -1, MAX_SIGNED_COUNTER)
        majflt = parse_sint(
            row.get("majflt_delta"), "major-fault delta",
            -1, MAX_SIGNED_COUNTER)
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
        cpu_before = parse_sint(
            row.get("cpu_before"), "cpu_before", -1, MAX_CPU_ID)
        cpu_after = parse_sint(
            row.get("cpu_after"), "cpu_after", -1, MAX_CPU_ID)
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
        if intermediate < K * width or intermediate % width:
            raise TimingError("intermediate storage size is inconsistent")
        signatures.add(tuple(row[field] for field in WORK_FIELDS))
    if len(signatures) != 1:
        raise TimingError("identical internal arms changed deterministic work")
    signature = next(iter(signatures))
    normalized_preamble = {
        key: value for key, value in preamble.items() if key != "schema"
    }
    semantic = {
        "preamble": normalized_preamble,
        "work_fields": list(WORK_FIELDS), "work_signature": list(signature),
        "row_count": len(rows), "timed_cycles": [1, 2, 3],
    }
    return ParsedOutput(
        label=label, schema=expected_schema, preamble=preamble, rows=rows,
        work_signature=signature,
        semantic_sha256=sha256_bytes(canonical_json(semantic)),
        stdout_sha256=sha256_bytes(raw), timed_elapsed_ns=timed_elapsed,
        all_elapsed_ns=all_elapsed, timed_phase_ns=timed_phase,
        all_phase_ns=all_phase, timed_minor_faults=timed_minor,
        discard_minor_faults=discard_minor, max_minor_faults=max_minor,
        contaminations=tuple(contaminations),
    )


def execution_name(task: Mapping[str, object], slot: int, label: str) -> str:
    return "%s.outer%d.%s.csv" % (task["task_id"], slot, label)


def command_for(
    design: Mapping[str, object], task: Mapping[str, object], label: str,
) -> List[str]:
    root = Path(str(design["root"]))
    tools = design["tools"]
    if label not in BINARY_NAMES:
        raise TimingError("unknown binary label")
    binary = root / "frozen" / BINARY_NAMES[label]
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
        "--candidate-period", str(ARCHITECTURE["period"]),
        "--candidate-grouped-rows", str(ARCHITECTURE["grouped_rows"]),
        "--candidate-buckets", str(ARCHITECTURE["buckets"]),
        "--evict-bytes", str(evict), "--cache-state", str(task["cache_state"]),
        "--loss", LOSS_TEXT, "--seed", str(task["seed"]),
        "--schedule", str(task["schedule"]),
    ]


def _receipt_summary(parsed: ParsedOutput) -> Dict[str, object]:
    return {
        "schema_version": parsed.schema,
        "stdout_sha256": parsed.stdout_sha256,
        "semantic_sha256": parsed.semantic_sha256,
        "trace_sha256": parsed.preamble["trace_sha256"],
        "work_signature_sha256": sha256_bytes(canonical_json(
            {"fields": list(WORK_FIELDS), "values": list(parsed.work_signature)})),
        "timed_elapsed_ns": parsed.timed_elapsed_ns,
        "all_elapsed_ns": parsed.all_elapsed_ns,
        "timed_phase_ns": dict(parsed.timed_phase_ns),
        "all_phase_ns": dict(parsed.all_phase_ns),
        "timed_minor_faults": parsed.timed_minor_faults,
        "discard_minor_faults": parsed.discard_minor_faults,
        "max_minor_faults": parsed.max_minor_faults,
        "row_count": len(parsed.rows), "timed_row_count": 24,
    }


def _ratio_record_or_none(
    numerator: int, denominator: int,
) -> Optional[Dict[str, int]]:
    if numerator < 0 or denominator < 0:
        raise TimingError("timing ratio contains a negative value")
    if denominator == 0:
        return None
    return {"numerator": numerator, "denominator": denominator}


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
        "work_signature": list(parsed.work_signature),
        "graph": {
            prefix + "_" + suffix: parsed.preamble[prefix + "_" + suffix]
            for prefix in ("control", "candidate")
            for suffix in ("attempt", "matrix_seed", "peel_seed")
        },
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


def _git_blob(git: Path, repo: Path, revision: str) -> bytes:
    result = subprocess.run(
        (str(git), "-C", str(repo), "show", revision),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if result.returncode != 0 or result.stderr:
        raise TimingError(
            "git could not read overlay blob %s: exit=%d stderr=%r" %
            (revision, result.returncode, result.stderr[:1000]))
    return result.stdout


def _stable_patch_id(git: Path, diff: bytes) -> str:
    if not diff:
        raise TimingError("measurement overlay diff is empty")
    result = subprocess.run(
        (str(git), "patch-id", "--stable"), input=diff,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    match = re.fullmatch(
        rb"([0-9a-f]{40,64}) [0-9a-f]{40,64}\n", result.stdout)
    if result.returncode != 0 or result.stderr or match is None:
        raise TimingError(
            "git patch-id rejected measurement overlay: exit=%d stderr=%r "
            "stdout=%r" %
            (result.returncode, result.stderr[:1000], result.stdout[:1000]))
    return match.group(1).decode("ascii")


def _overlay_identity(record: object) -> Dict[str, object]:
    """Validate and return the label-independent overlay provenance."""
    identity_fields = (
        "source_commit", "source_parent_commit", "source_tree", "files",
        "diff_options", "diff_sha256", "stable_patch_id",
    )
    all_fields = set(identity_fields) | {"diff_evidence_name"}
    if not isinstance(record, dict) or set(record) != all_fields:
        raise TimingError("measurement overlay provenance fields changed")
    if (record.get("source_commit") != MEASUREMENT_OVERLAY_COMMIT or
            record.get("source_parent_commit") !=
            MEASUREMENT_OVERLAY_PARENT_COMMIT):
        raise TimingError("measurement overlay lineage changed")
    source_tree = record.get("source_tree")
    if not isinstance(source_tree, str) or GIT_OID_RE.fullmatch(source_tree) is None:
        raise TimingError("measurement overlay tree is malformed")
    expected_options = list(MEASUREMENT_DIFF_OPTIONS)
    if record.get("diff_options") != expected_options:
        raise TimingError("measurement overlay diff options changed")
    diff_name = record.get("diff_evidence_name")
    if (not isinstance(diff_name, str) or "/" in diff_name or
            not diff_name.endswith(".measurement-overlay.diff")):
        raise TimingError("measurement overlay evidence name is malformed")
    require_sha256(record.get("diff_sha256"), "measurement overlay diff hash")
    patch_id = record.get("stable_patch_id")
    if not isinstance(patch_id, str) or GIT_OID_RE.fullmatch(patch_id) is None:
        raise TimingError("measurement overlay patch-id is malformed")
    files = record.get("files")
    if not isinstance(files, list) or len(files) != len(MEASUREMENT_OVERLAY_FILES):
        raise TimingError("measurement overlay file ledger is malformed")
    if [item.get("path") if isinstance(item, dict) else None
            for item in files] != list(MEASUREMENT_OVERLAY_FILES):
        raise TimingError("measurement overlay file list changed")
    expected_file_fields = {
        "path", "parent_blob_oid", "parent_sha256",
        "overlay_blob_oid", "overlay_sha256",
    }
    for item in files:
        if not isinstance(item, dict) or set(item) != expected_file_fields:
            raise TimingError("measurement overlay file receipt changed")
        for key in ("parent_blob_oid", "overlay_blob_oid"):
            value = item.get(key)
            if not isinstance(value, str) or GIT_OID_RE.fullmatch(value) is None:
                raise TimingError("measurement overlay blob ID is malformed")
        require_sha256(item.get("parent_sha256"), "parent overlay-file hash")
        require_sha256(item.get("overlay_sha256"), "overlay-file hash")
        if item["parent_blob_oid"] == item["overlay_blob_oid"] or \
                item["parent_sha256"] == item["overlay_sha256"]:
            raise TimingError("measurement overlay did not change a listed file")
    return {field: record[field] for field in identity_fields}


def _expected_overlay_status() -> set[str]:
    return {" M " + path for path in MEASUREMENT_OVERLAY_FILES}


def _git_status_lines(git: Path, source: Path) -> Tuple[str, ...]:
    result = subprocess.run(
        (str(git), "-C", str(source), "status", "--porcelain=v1",
         "--untracked-files=all"),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if result.returncode != 0 or result.stderr:
        raise TimingError(
            "git could not inspect overlay worktree: exit=%d stderr=%r" %
            (result.returncode, result.stderr[:1000]))
    try:
        text = result.stdout.decode("ascii")
    except UnicodeDecodeError as exc:
        raise TimingError("overlay worktree status is not ASCII") from exc
    if text and not text.endswith("\n"):
        raise TimingError("overlay worktree status is not canonical text")
    return tuple(text.splitlines())


def _verify_overlay_status(
    git: Path, source: Path,
    files: Optional[Sequence[Mapping[str, object]]] = None,
) -> None:
    status = _git_status_lines(git, source)
    if set(status) != _expected_overlay_status() or \
            len(status) != len(MEASUREMENT_OVERLAY_FILES):
        raise TimingError("build worktree differs from the exact measurement overlay")
    if files is not None:
        if [item.get("path") for item in files] != \
                list(MEASUREMENT_OVERLAY_FILES):
            raise TimingError("overlay content ledger does not match its paths")
        for item in files:
            path = source / str(item["path"])
            if path.is_symlink() or not path.is_file() or \
                    sha256_file(path) != item.get("overlay_sha256"):
                raise TimingError(
                    "build worktree overlay content changed: %s" % item["path"])


def _apply_measurement_overlay(
    label: str, repo: Path, source: Path, provenance_dir: Path,
    git: Path,
) -> Dict[str, object]:
    if _git_status_lines(git, source):
        raise TimingError("detached codec-parent worktree is not initially clean")
    files: List[Dict[str, object]] = []
    for relative in MEASUREMENT_OVERLAY_FILES:
        destination = source / relative
        if destination.is_symlink() or not destination.is_file():
            raise TimingError("overlay target is not a regular file: %s" % relative)
        parent_revision = "HEAD:" + relative
        overlay_revision = MEASUREMENT_OVERLAY_COMMIT + ":" + relative
        parent = _git_blob(git, source, parent_revision)
        overlay = _git_blob(git, repo, overlay_revision)
        if stable_bytes(
                destination, attempts=3, max_bytes=MAX_EVIDENCE_FILE_BYTES,
                require_unique=True) != parent:
            raise TimingError("codec-parent file disagrees with its Git blob")
        destination.write_bytes(overlay)
        files.append({
            "path": relative,
            "parent_blob_oid": _git_value(
                git, source, "rev-parse", parent_revision),
            "parent_sha256": sha256_bytes(parent),
            "overlay_blob_oid": _git_value(
                git, repo, "rev-parse", overlay_revision),
            "overlay_sha256": sha256_bytes(overlay),
        })
    _verify_overlay_status(git, source, files)
    diff_options = list(MEASUREMENT_DIFF_OPTIONS)
    diff_result = subprocess.run(
        (str(git), "-C", str(source), "diff", *diff_options, "--",
         *MEASUREMENT_OVERLAY_FILES),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if diff_result.returncode != 0 or diff_result.stderr or not diff_result.stdout:
        raise TimingError(
            "could not capture measurement overlay diff: exit=%d stderr=%r" %
            (diff_result.returncode, diff_result.stderr[:1000]))
    diff_name = label + ".measurement-overlay.diff"
    write_new(provenance_dir / diff_name, diff_result.stdout)
    record = {
        "source_commit": MEASUREMENT_OVERLAY_COMMIT,
        "source_parent_commit": _git_value(
            git, repo, "rev-parse", MEASUREMENT_OVERLAY_COMMIT + "^"),
        "source_tree": _git_value(
            git, repo, "rev-parse", MEASUREMENT_OVERLAY_COMMIT + "^{tree}"),
        "files": files, "diff_options": diff_options,
        "diff_sha256": sha256_bytes(diff_result.stdout),
        "stable_patch_id": _stable_patch_id(git, diff_result.stdout),
        "diff_evidence_name": diff_name,
    }
    _overlay_identity(record)
    return record


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
        measurement_overlay = _apply_measurement_overlay(
            label, repo, source, provenance_dir, tools["git"])
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
            "-DCMAKE_C_FLAGS_RELEASE=" + TIMING_C_FLAGS_RELEASE,
            "-DCMAKE_CXX_FLAGS_RELEASE=" + TIMING_CXX_FLAGS_RELEASE,
            "-DCMAKE_EXE_LINKER_FLAGS=" + TIMING_EXE_LINKER_FLAGS,
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
        _verify_overlay_status(
            tools["git"], source, measurement_overlay["files"])
        binary = build / "codec/wirehair_v2_bench"
        if not binary.is_file() or not os.access(binary, os.X_OK):
            raise TimingError("%s benchmark binary was not produced" % label)
        frozen_binary = frozen / BINARY_NAMES[label]
        write_new(
            frozen_binary,
            stable_bytes(binary, attempts=3),
            mode=0o555)
        cache = build / "CMakeCache.txt"
        commands = build / "compile_commands.json"
        for source_path, suffix in (
                (cache, "CMakeCache.txt"), (commands, "compile_commands.json")):
            if not source_path.is_file():
                raise TimingError("%s build provenance is missing %s" %
                                  (label, suffix))
            write_new(
                provenance_dir / (label + "." + suffix),
                stable_bytes(source_path, attempts=3))
        ldd_text, dependencies = _dependency_hashes(frozen_binary, tools["ldd"])
        write_new(provenance_dir / (label + ".ldd.txt"),
                  (ldd_text + "\n").encode("utf-8"))
        nm_text = checked_text((
            str(tools["nm"]), "-C", str(frozen_binary)))
        forbidden = [
            symbol for symbol in FORBIDDEN_TIMED_BINARY_SYMBOLS
            if symbol in nm_text
        ]
        if forbidden:
            raise TimingError(
                "%s timed binary retained forbidden test/CLI code: %s" %
                (label, ",".join(forbidden)))
        write_new(
            provenance_dir / (label + ".nm.txt"),
            (nm_text + "\n").encode("utf-8"))
        rejected_mode = command_result(
            (str(frozen_binary), "selftest"),
            environment=build_environment)
        if (rejected_mode["returncode"] != 1 or rejected_mode["stdout"] or
                rejected_mode["stderr"] != b"unknown mode: selftest\n"):
            raise TimingError(
                "%s timed binary did not reject an unrelated CLI mode" %
                label)
        write_new(
            provenance_dir / (label + ".rejected-mode.stderr"),
            rejected_mode["stderr"])
        compiler_version = checked_text((str(cxx_compiler), "--version"))
        compiler_target = checked_text((str(cxx_compiler), "-dumpmachine"))
        compiler_numeric = checked_text((
            str(cxx_compiler), "-dumpfullversion", "-dumpversion"))
        log_names = (
            label + ".configure.stdout", label + ".configure.stderr",
            label + ".build.stdout", label + ".build.stderr",
            label + ".CMakeCache.txt", label + ".compile_commands.json",
            label + ".ldd.txt", label + ".nm.txt",
            label + ".rejected-mode.stderr",
            label + ".measurement-overlay.diff",
        )
        payload: Dict[str, object] = {
            "label": label, "codec_parent_commit": commit,
            "codec_parent_tree": tree,
            "source_exact_overlay_before_and_after": True,
            "measurement_overlay": measurement_overlay,
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
            "timing_only_binary": True,
            "timing_build_flags": {
                "c_release": TIMING_C_FLAGS_RELEASE,
                "cxx_release": TIMING_CXX_FLAGS_RELEASE,
                "exe_linker": TIMING_EXE_LINKER_FLAGS,
            },
            "forbidden_symbols_absent": list(
                FORBIDDEN_TIMED_BINARY_SYMBOLS),
            "unrelated_mode_rejected": True,
            "evidence_files": {
                name: sha256_file(provenance_dir / name) for name in log_names
            },
        }
        provenance = sealed_record(
            "wirehair.wh2.grouped_commit_timing.build_provenance.v1", payload)
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


def _prepare_cross_binary_smoke(
    staging: Path, provenance_dir: Path, tools: Mapping[str, Path],
    core: int, numa_node: int,
) -> Dict[str, object]:
    """Prove both frozen v1 CLIs and their normalizer before quiet timing."""
    task = next(
        value for value in generate_tasks()
        if value["K"] == 3200 and value["bb"] == 64 and
        value["schedule"] == "burst" and value["cache_state"] == "warm")
    tool_records = {
        name: {"path": str(path), "sha256": sha256_file(path)}
        for name, path in tools.items()
    }
    design = {
        "root": str(staging), "tools": tool_records, "core": core,
        "numa_node": numa_node, "evict_bytes": 4096,
    }
    parsed = []
    records = {}
    for label in ("base", "candidate"):
        command = command_for(design, task, label)
        result = subprocess.run(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            timeout=300, check=False, start_new_session=True)
        if result.returncode != 0 or result.stderr:
            raise TimingError(
                "%s prepare smoke failed exit=%d stderr=%r" %
                (label, result.returncode, result.stderr[:1000]))
        observation = parse_grouped_output(
            result.stdout, label, task, 4096, core)
        stdout_name = "smoke." + label + ".csv"
        write_new(provenance_dir / stdout_name, result.stdout)
        records[label] = {
            "argv": command, "stdout_name": stdout_name,
            "stdout_sha256": sha256_bytes(result.stdout),
            "schema_version": observation.schema,
            "semantic_sha256": observation.semantic_sha256,
            "trace_sha256": observation.preamble["trace_sha256"],
            "work_signature_sha256": sha256_bytes(canonical_json({
                "fields": list(WORK_FIELDS),
                "values": list(observation.work_signature),
            })),
            # This bounded preparation smoke runs under the normal full-load
            # regime and is evidence of compatibility, never speed evidence.
            "nonpromotional_contaminations": list(observation.contaminations),
        }
        parsed.append(observation)
    if (len({item.semantic_sha256 for item in parsed}) != 1 or
            len({item.work_signature for item in parsed}) != 1 or
            len({item.preamble["trace_sha256"] for item in parsed}) != 1):
        raise TimingError("prepare smoke found cross-binary trace or work drift")
    return {
        "scope": "bounded compatibility smoke under ordinary machine load",
        "timing_evidence": False, "task": task, "binaries": records,
        "semantic_identity_exact": True, "work_identity_exact": True,
        "trace_identity_exact": True,
    }


def prepare_campaign(args: argparse.Namespace) -> None:
    result = Path(args.result_dir).resolve()
    repo = Path(args.repo).resolve()
    if result.exists():
        raise TimingError("result directory already exists")
    if not repo.is_dir():
        raise TimingError("repository directory does not exist")
    thermal_sampler_source = Path(args.thermal_sampler).resolve()
    thermal_sampler_raw = stable_bytes(
        thermal_sampler_source, attempts=3, max_bytes=1024 * 1024,
        require_unique=False)
    if not thermal_sampler_raw.startswith(b"#!/usr/bin/env python3\n"):
        raise TimingError("thermal sampler source has an unexpected header")
    if (args.core < 0 or args.controller_core < 0 or
            args.numa_node < 0 or args.evict_bytes < 4096 or
            args.build_jobs <= 0):
        raise TimingError("prepare integer argument is outside its domain")
    if (CLOCK_TICKS_PER_SECOND <= 0 or
            SIBLING_ACCEPTED_EXECUTION_MIN_DURATION_NS !=
            (1_000_000_000 + CLOCK_TICKS_PER_SECOND - 1) //
            CLOCK_TICKS_PER_SECOND):
        raise TimingError("SC_CLK_TCK is outside its domain")
    require_pidfd_runtime()
    tool_names = ("git", "cmake", "env", "taskset", "numactl", "ldd", "nm",
                  "sudo", "fuser", "true")
    tools = {name: resolve_tool(name) for name in tool_names}
    python = Path(sys.executable).resolve()
    if not python.is_file() or not os.access(python, os.X_OK):
        raise TimingError("active Python interpreter is not an executable file")
    tools["python"] = python
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
            (BASE_COMMIT, "base"), (CANDIDATE_COMMIT, "candidate"),
            (MEASUREMENT_OVERLAY_PARENT_COMMIT, "measurement"),
            (MEASUREMENT_OVERLAY_COMMIT, "measurement overlay")):
        resolved = _git_value(tools["git"], repo, "rev-parse", commit + "^{commit}")
        if resolved != commit:
            raise TimingError("%s commit is unavailable or ambiguous" % name)
    if _git_value(
            tools["git"], repo, "rev-parse", CANDIDATE_COMMIT + "^") != \
            BASE_COMMIT:
        raise TimingError("candidate is not the direct child of the base codec parent")
    if _git_value(
            tools["git"], repo, "rev-parse",
            MEASUREMENT_OVERLAY_COMMIT + "^") != \
            MEASUREMENT_OVERLAY_PARENT_COMMIT or \
            _git_value(
                tools["git"], repo, "rev-parse",
                MEASUREMENT_OVERLAY_PARENT_COMMIT + "^") != \
            CANDIDATE_COMMIT:
        raise TimingError(
            "measurement overlay lineage is not candidate -> measurement -> "
            "timing-only")
    topology = topology_record(args.core, args.numa_node)
    controller_topology = topology_record(args.controller_core, args.numa_node)
    if args.controller_core in topology["llc_shared_cpus"]:
        raise TimingError("controller CPU shares the timing LLC")
    if args.evict_bytes < max(2 * int(topology["llc_bytes"]), 256 * 1024 * 1024):
        raise TimingError("eviction allocation is smaller than the frozen LLC gate")
    staging = result.with_name(result.name + ".prepare.%d" % os.getpid())
    workspace = result.with_name(result.name + ".build.%d" % os.getpid())
    if staging.exists() or workspace.exists():
        raise TimingError("stale prepare workspace exists")
    staging.mkdir(mode=0o700)
    workspace.mkdir(mode=0o700)
    os.chmod(staging, 0o700)
    os.chmod(workspace, 0o700)
    staging_sealed = False
    try:
        frozen = staging / "frozen"
        provenance_dir = staging / "provenance"
        frozen.mkdir(mode=0o700)
        provenance_dir.mkdir(mode=0o700)
        for directory in (
                "raw", "stderr", "exit", "receipts", "task_receipts",
                "contamination", "failure"):
            (staging / directory).mkdir(mode=0o700)
        harness_source = Path(__file__).resolve()
        helper_source = Path(str(evidence_io.__file__)).resolve()
        if (helper_source.parent != harness_source.parent or
                helper_source.name != "wh2_timing_evidence_io.py"):
            raise TimingError("timing evidence-I/O helper is not beside the harness")
        harness_frozen = frozen / "wh2_grouped_commit_timing.py"
        helper_frozen = frozen / "wh2_timing_evidence_io.py"
        sampler_frozen = frozen / THERMAL_SAMPLER_NAME
        write_new(
            harness_frozen,
            stable_bytes(harness_source, attempts=3),
            mode=0o555)
        write_new(
            helper_frozen,
            stable_bytes(helper_source, attempts=3))
        write_new(sampler_frozen, thermal_sampler_raw, mode=0o555)
        builds = {}
        for label, commit in (("base", BASE_COMMIT),
                              ("candidate", CANDIDATE_COMMIT)):
            builds[label] = _build_one(
                label, commit, repo, workspace, frozen, provenance_dir,
                tools, args.build_jobs, c_compiler, cxx_compiler)
        base_overlay = _overlay_identity(
            builds["base"]["record"].get("measurement_overlay"))
        candidate_overlay = _overlay_identity(
            builds["candidate"]["record"].get("measurement_overlay"))
        if base_overlay != candidate_overlay:
            raise TimingError(
                "base and candidate did not receive an identical measurement overlay")
        prepare_smoke = _prepare_cross_binary_smoke(
            staging, provenance_dir, tools, args.core, args.numa_node)
        if (sha256_file(harness_source) != sha256_file(harness_frozen) or
                sha256_file(helper_source) != sha256_file(helper_frozen)):
            raise TimingError("timing controller source changed during preparation")
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
            "candidate_commit": CANDIDATE_COMMIT,
            "codec_lineage": {
                "candidate_parent_commit": BASE_COMMIT,
                "measurement_overlay_commit":
                    MEASUREMENT_OVERLAY_PARENT_COMMIT,
                "measurement_overlay_parent_commit": CANDIDATE_COMMIT,
                "timing_only_overlay_commit": MEASUREMENT_OVERLAY_COMMIT,
                "timing_only_overlay_parent_commit":
                    MEASUREMENT_OVERLAY_PARENT_COMMIT,
            },
            "measurement_overlay": base_overlay,
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
            "timed_rows_per_binary_per_task": 96,
            "allocator_environment": {
                "MALLOC_MMAP_THRESHOLD_": MALLOC_MMAP_THRESHOLD,
                "MALLOC_TRIM_THRESHOLD_": MALLOC_TRIM_THRESHOLD,
                "purpose": "retain and prefault solve-output pages in each process while preserving all timed decoder work",
            },
            "phase_timing": {
                "primary": PRIMARY_PHASE_FIELD,
                "negative_controls": list(NEGATIVE_CONTROL_PHASE_FIELDS),
                "fields": list(PHASE_FIELDS),
            },
            "timed_binary_build": {
                "grouped_timing_only": True,
                "c_release_flags": TIMING_C_FLAGS_RELEASE,
                "cxx_release_flags": TIMING_CXX_FLAGS_RELEASE,
                "exe_linker_flags": TIMING_EXE_LINKER_FLAGS,
                "forbidden_symbols_absent": list(
                    FORBIDDEN_TIMED_BINARY_SYMBOLS),
                "unrelated_modes_rejected": True,
            },
            "minor_fault_policy": {
                "minimum": 0, "maximum": MAX_MINOR_FAULTS,
                "major_faults": 0, "all_cycles_receipted": True,
            },
            "max_environmental_attempts": MAX_ENVIRONMENTAL_ATTEMPTS,
            "sibling_idle_policy": SIBLING_IDLE_POLICY,
            "core": args.core, "numa_node": args.numa_node,
            "topology": topology, "controller_core": args.controller_core,
            "controller_topology": controller_topology,
            "evict_bytes": args.evict_bytes,
            "thermal_limits_c": {
                "cpu": MAX_CPU_TEMP_C, "dimm": MAX_DIMM_TEMP_C,
                "max_gap_s": MAX_THERMAL_GAP_S,
                "max_coverage_margin_s": MAX_THERMAL_MARGIN_S,
            },
            "thermal_sampler": {
                "name": THERMAL_SAMPLER_NAME,
                "sha256": sha256_bytes(thermal_sampler_raw),
            },
            "fresh_only": True,
            "tasks_manifest_sha256": sha256_bytes(manifest),
            "immutable_files": immutable_files, "tools": tool_records,
            "build_provenance_sha256": {
                label: sha256_file(value["path"])
                for label, value in builds.items()
            },
            "prepare_smoke": prepare_smoke,
        }
        design = sealed_record(
            "wirehair.wh2.grouped_commit_timing.design.v3", design_payload)
        design_path = staging / "design.json"
        write_new(design_path, canonical_json(design))
        receipt = sealed_record(
            "wirehair.wh2.grouped_commit_timing.prepare_receipt.v1", {
                "prepared_utc": utc_now(), "design_sha256": sha256_file(design_path),
                "tasks_manifest_sha256": sha256_bytes(manifest),
                "immutable_files": immutable_files,
                "base_binary_sha256": immutable_files[
                    "frozen/" + BINARY_NAMES["base"]],
                "candidate_binary_sha256": immutable_files[
                    "frozen/" + BINARY_NAMES["candidate"]],
            })
        write_new(staging / "prepare_receipt.json", canonical_json(receipt))
        for directory in (frozen, provenance_dir):
            os.chmod(directory, 0o500)
            with evidence_io.held_private_directory(
                    directory, require_writable=False,
                    error_type=TimingError) as directory_fd:
                os.fsync(directory_fd)
        os.chmod(staging, 0o500)
        staging_sealed = True
        evidence_io.publish_directory_noreplace(
            staging, result, error_type=TimingError)
        staging_sealed = False
        os.chmod(result, 0o700)
        with evidence_io.held_private_directory(
                result, require_writable=True,
                error_type=TimingError) as result_fd:
            os.fsync(result_fd)
        with evidence_io.held_private_directory(
                result.parent, require_writable=True,
                error_type=TimingError) as parent_fd:
            os.fsync(parent_fd)
        print(json.dumps({
            "result_dir": str(result), "task_count": len(tasks),
            "design_sha256": sha256_file(result / "design.json"),
            "manifest_sha256": sha256_file(result / "tasks_manifest.jsonl"),
            "prepare_receipt_sha256":
                sha256_file(result / "prepare_receipt.json"),
            "base_binary_sha256": receipt["base_binary_sha256"],
            "candidate_binary_sha256": receipt["candidate_binary_sha256"],
        }, sort_keys=True))
    finally:
        if staging.exists() and not staging_sealed:
            shutil.rmtree(staging)
        if workspace.exists():
            shutil.rmtree(workspace)


def require_frozen_sibling_idle_policy(value: object) -> None:
    """Require exact JSON types as well as values for the proof policy."""
    if canonical_json(value) != canonical_json(SIBLING_IDLE_POLICY):
        raise TimingError("timing design policy changed: sibling_idle_policy")


def _load_design(root: Path) -> Dict[str, object]:
    design = load_canonical(root / "design.json", "timing design")
    verify_sealed_record(
        design, "wirehair.wh2.grouped_commit_timing.design.v3", "timing design")
    if design.get("root") != str(root.resolve()):
        raise TimingError("timing root moved after preparation")
    if design.get("base_commit") != BASE_COMMIT or \
            design.get("candidate_commit") != CANDIDATE_COMMIT:
        raise TimingError("timing commit identities changed")
    expected = {
        "codec_lineage": {
            "candidate_parent_commit": BASE_COMMIT,
            "measurement_overlay_commit":
                MEASUREMENT_OVERLAY_PARENT_COMMIT,
            "measurement_overlay_parent_commit": CANDIDATE_COMMIT,
            "timing_only_overlay_commit": MEASUREMENT_OVERLAY_COMMIT,
            "timing_only_overlay_parent_commit":
                MEASUREMENT_OVERLAY_PARENT_COMMIT,
        },
        "timing_scope": "full-payload decoder precode solve",
        "architecture": ARCHITECTURE, "K": list(KS), "bb": list(WIDTHS),
        "schedule_seeds": [list(item) for item in SCHEDULE_SEEDS],
        "cache_states": list(CACHE_STATES), "overhead": OVERHEAD,
        "loss": LOSS_TEXT, "outer_order": OUTER_ORDER,
        "inner_order": INNER_ORDER, "inner_cycles": 4,
        "discard_inner_cycle": 0, "task_count": 108,
        "processes_per_task": 8, "rows_per_process": 32,
        "timed_rows_per_process": 24,
        "timed_rows_per_binary_per_task": 96,
        "allocator_environment": {
            "MALLOC_MMAP_THRESHOLD_": MALLOC_MMAP_THRESHOLD,
            "MALLOC_TRIM_THRESHOLD_": MALLOC_TRIM_THRESHOLD,
            "purpose": "retain and prefault solve-output pages in each process while preserving all timed decoder work",
        },
        "phase_timing": {
            "primary": PRIMARY_PHASE_FIELD,
            "negative_controls": list(NEGATIVE_CONTROL_PHASE_FIELDS),
            "fields": list(PHASE_FIELDS),
        },
        "timed_binary_build": {
            "grouped_timing_only": True,
            "c_release_flags": TIMING_C_FLAGS_RELEASE,
            "cxx_release_flags": TIMING_CXX_FLAGS_RELEASE,
            "exe_linker_flags": TIMING_EXE_LINKER_FLAGS,
            "forbidden_symbols_absent": list(
                FORBIDDEN_TIMED_BINARY_SYMBOLS),
            "unrelated_modes_rejected": True,
        },
        "minor_fault_policy": {
            "minimum": 0, "maximum": MAX_MINOR_FAULTS,
            "major_faults": 0, "all_cycles_receipted": True,
        },
        "max_environmental_attempts": MAX_ENVIRONMENTAL_ATTEMPTS,
        "sibling_idle_policy": SIBLING_IDLE_POLICY,
        "thermal_limits_c": {
            "cpu": MAX_CPU_TEMP_C, "dimm": MAX_DIMM_TEMP_C,
            "max_gap_s": MAX_THERMAL_GAP_S,
            "max_coverage_margin_s": MAX_THERMAL_MARGIN_S,
        },
        "thermal_sampler": {
            "name": THERMAL_SAMPLER_NAME,
            "sha256": design.get("immutable_files", {}).get(
                "frozen/" + THERMAL_SAMPLER_NAME),
        },
        "fresh_only": True,
    }
    for key, value in expected.items():
        if key == "sibling_idle_policy":
            require_frozen_sibling_idle_policy(design.get(key))
        elif design.get(key) != value:
            raise TimingError("timing design policy changed: %s" % key)
    overlay = design.get("measurement_overlay")
    if _overlay_identity({
            **overlay, "diff_evidence_name": "design.measurement-overlay.diff"
    } if isinstance(overlay, dict) else overlay) != overlay:
        raise TimingError("timing design measurement-overlay receipt is malformed")
    for key in ("core", "controller_core", "numa_node", "evict_bytes"):
        value = design.get(key)
        if not isinstance(value, int) or isinstance(value, bool) or value < 0:
            raise TimingError("timing design integer is malformed: %s" % key)
    topology = design.get("topology")
    controller_topology = design.get("controller_topology")
    if not isinstance(topology, dict) or not isinstance(controller_topology, dict):
        raise TimingError("timing design topology receipt is malformed")
    if design["evict_bytes"] < 4096 or \
            design["controller_core"] in topology.get("llc_shared_cpus", []):
        raise TimingError("timing design isolation domain is malformed")
    smoke = design.get("prepare_smoke")
    expected_smoke_task = next(
        value for value in generate_tasks()
        if value["K"] == 3200 and value["bb"] == 64 and
        value["schedule"] == "burst" and value["cache_state"] == "warm")
    if (not isinstance(smoke, dict) or smoke.get("timing_evidence") is not False or
            smoke.get("task") != expected_smoke_task or
            smoke.get("semantic_identity_exact") is not True or
            smoke.get("work_identity_exact") is not True or
            smoke.get("trace_identity_exact") is not True):
        raise TimingError("prepare compatibility smoke receipt is malformed")
    return design


def _load_tasks(root: Path, design: Mapping[str, object]) -> Tuple[Dict[str, object], ...]:
    raw = stable_bytes(
        root / "tasks_manifest.jsonl",
        max_bytes=MAX_TASK_MANIFEST_BYTES)
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
            "frozen/wh2_grouped_commit_timing.py"):
        raise TimingError("active timing harness is not the frozen harness")
    active_helper = Path(str(evidence_io.__file__)).resolve()
    if (active_helper.parent != active_harness.parent or
            sha256_file(active_helper) != expected.get(
                "frozen/wh2_timing_evidence_io.py")):
        raise TimingError("active evidence-I/O helper is not the frozen helper")
    tools = design.get("tools")
    if not isinstance(tools, dict):
        raise TimingError("tool ledger is malformed")
    for name, record in tools.items():
        if not isinstance(record, dict) or sha256_file(Path(str(record.get("path")))) \
                != record.get("sha256"):
            raise TimingError("campaign tool changed: %s" % name)
    python_record = tools.get("python")
    if (not isinstance(python_record, dict) or
            Path(sys.executable).resolve() !=
            Path(str(python_record.get("path"))).resolve()):
        raise TimingError("campaign is running under a different Python interpreter")
    git_record = tools.get("git")
    if not isinstance(git_record, dict):
        raise TimingError("frozen Git tool receipt is missing")
    git = Path(str(git_record.get("path")))
    build_hashes = design.get("build_provenance_sha256")
    if not isinstance(build_hashes, dict) or set(build_hashes) != {
            "base", "candidate"}:
        raise TimingError("build-provenance hash ledger is malformed")
    for label, commit in (("base", BASE_COMMIT),
                          ("candidate", CANDIDATE_COMMIT)):
        provenance_path = root / "provenance" / (label + ".json")
        provenance = load_canonical(
            provenance_path, label + " build provenance")
        verify_sealed_record(
            provenance, "wirehair.wh2.grouped_commit_timing.build_provenance.v1",
            label + " build provenance")
        parent_tree = provenance.get("codec_parent_tree")
        if (provenance.get("label") != label or
                provenance.get("codec_parent_commit") != commit or
                not isinstance(parent_tree, str) or
                GIT_OID_RE.fullmatch(parent_tree) is None or
                provenance.get("source_exact_overlay_before_and_after") is not True or
                build_hashes.get(label) != sha256_file(provenance_path)):
            raise TimingError("build provenance identity mismatch")
        overlay = provenance.get("measurement_overlay")
        overlay_identity = _overlay_identity(overlay)
        if overlay_identity != design.get("measurement_overlay"):
            raise TimingError("build used a different measurement overlay")
        if (provenance.get("timing_only_binary") is not True or
                provenance.get("timing_build_flags") != {
                    "c_release": TIMING_C_FLAGS_RELEASE,
                    "cxx_release": TIMING_CXX_FLAGS_RELEASE,
                    "exe_linker": TIMING_EXE_LINKER_FLAGS,
                } or
                provenance.get("forbidden_symbols_absent") !=
                list(FORBIDDEN_TIMED_BINARY_SYMBOLS) or
                provenance.get("unrelated_mode_rejected") is not True):
            raise TimingError("timed binary build policy changed")
        expected_diff_name = label + ".measurement-overlay.diff"
        if overlay.get("diff_evidence_name") != expected_diff_name:
            raise TimingError("measurement overlay evidence label changed")
        overlay_diff = stable_bytes(
            root / "provenance" / expected_diff_name,
            max_bytes=MAX_EVIDENCE_FILE_BYTES)
        if (sha256_bytes(overlay_diff) != overlay.get("diff_sha256") or
                _stable_patch_id(git, overlay_diff) !=
                overlay.get("stable_patch_id")):
            raise TimingError("measurement overlay diff provenance changed")
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
        try:
            nm_text = stable_bytes(
                root / "provenance" / (label + ".nm.txt"),
                max_bytes=MAX_EVIDENCE_FILE_BYTES).decode("utf-8", "strict")
        except UnicodeDecodeError as exc:
            raise TimingError("timed binary symbol evidence is not UTF-8") \
                from exc
        if any(symbol in nm_text for symbol in FORBIDDEN_TIMED_BINARY_SYMBOLS):
            raise TimingError("timed binary symbol exclusion no longer replays")
        dependencies = provenance.get("runtime_dependency_sha256")
        if not isinstance(dependencies, dict) or not dependencies:
            raise TimingError("runtime dependency ledger is missing")
        for path, digest in dependencies.items():
            if not isinstance(path, str) or sha256_file(Path(path)) != digest:
                raise TimingError("runtime dependency changed: %s" % path)
    smoke = design["prepare_smoke"]
    smoke_task = smoke.get("task")
    smoke_records = smoke.get("binaries")
    if not isinstance(smoke_task, dict) or not isinstance(smoke_records, dict):
        raise TimingError("prepare smoke ledger is malformed")
    parsed_smoke = []
    for label in ("base", "candidate"):
        record = smoke_records.get(label)
        expected_name = "smoke." + label + ".csv"
        if not isinstance(record, dict) or record.get("stdout_name") != expected_name:
            raise TimingError("prepare smoke binary ledger is malformed")
        raw = stable_bytes(
            root / "provenance" / expected_name,
            max_bytes=MAX_GROUPED_OUTPUT_BYTES)
        parsed = parse_grouped_output(
            raw, label, smoke_task, 4096, int(design["core"]))
        expected_record = {
            "stdout_sha256": parsed.stdout_sha256,
            "schema_version": parsed.schema,
            "semantic_sha256": parsed.semantic_sha256,
            "trace_sha256": parsed.preamble["trace_sha256"],
            "work_signature_sha256": sha256_bytes(canonical_json({
                "fields": list(WORK_FIELDS),
                "values": list(parsed.work_signature),
            })),
            "nonpromotional_contaminations": list(parsed.contaminations),
        }
        if any(record.get(key) != value for key, value in expected_record.items()):
            raise TimingError("prepare smoke receipt does not replay")
        parsed_smoke.append(parsed)
    if (len({item.semantic_sha256 for item in parsed_smoke}) != 1 or
            len({item.work_signature for item in parsed_smoke}) != 1 or
            len({item.preamble["trace_sha256"] for item in parsed_smoke}) != 1):
        raise TimingError("prepare smoke cross-binary identity changed")


def _validate_prepare_receipt(
    root: Path, design: Mapping[str, object],
) -> Dict[str, object]:
    prepare = load_canonical(root / "prepare_receipt.json", "prepare receipt")
    verify_sealed_record(
        prepare, "wirehair.wh2.grouped_commit_timing.prepare_receipt.v1",
        "prepare receipt")
    if set(prepare) != PREPARE_RECEIPT_FIELDS:
        raise TimingError("prepare receipt fields changed")
    if (prepare.get("design_sha256") != sha256_file(root / "design.json") or
            prepare.get("tasks_manifest_sha256") !=
            sha256_file(root / "tasks_manifest.jsonl") or
            prepare.get("immutable_files") != design.get("immutable_files") or
            prepare.get("base_binary_sha256") != design["immutable_files"].get(
                "frozen/" + BINARY_NAMES["base"]) or
            prepare.get("candidate_binary_sha256") !=
            design["immutable_files"].get(
                "frozen/" + BINARY_NAMES["candidate"]) or
            not isinstance(prepare.get("prepared_utc"), str) or
            UTC_RE.fullmatch(prepare["prepared_utc"]) is None):
        raise TimingError("prepare receipt does not bind the frozen campaign")
    return prepare


def _require_external_prepare_anchor(
    root: Path, expected_sha256: object,
) -> str:
    expected = require_sha256(
        expected_sha256, "externally supplied prepare-receipt SHA256")
    actual = sha256_file(
        root / "prepare_receipt.json",
        max_bytes=MAX_JSON_EVIDENCE_BYTES, require_unique=True)
    if actual != expected:
        raise TimingError(
            "externally supplied prepare-receipt SHA256 does not match")
    return actual


def _verify_directory_inventory(root: Path) -> None:
    expected_directories = set(PREPARED_CAMPAIGN_DIRECTORIES)
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
            "contamination", "failure"):
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


def _schedstat_uint(token: bytes, context: str) -> int:
    if (len(token) > 20 or
            re.fullmatch(rb"0|[1-9][0-9]*", token) is None):
        raise TimingError("%s is not a canonical unsigned integer" % context)
    value = int(token)
    if value > MAX_UNSIGNED_COUNTER:
        raise TimingError("%s is outside the schedstat counter domain" % context)
    return value


def parse_schedstat_cpu(raw: bytes, cpu: int) -> Dict[str, int]:
    """Strictly parse Linux schedstat v15 and return one CPU's counters."""
    if (not isinstance(raw, bytes) or not raw or not raw.endswith(b"\n") or
            len(raw) > MAX_JSON_EVIDENCE_BYTES or b"\r" in raw or
            b"\0" in raw or not isinstance(cpu, int) or
            isinstance(cpu, bool) or not 0 <= cpu <= MAX_CPU_ID):
        raise TimingError("schedstat input is malformed")
    lines = raw.splitlines()
    if (len(lines) < 3 or lines[0] != b"version 15" or
            len(lines[1].split(b" ")) != 2 or
            lines[1].split(b" ")[0] != b"timestamp"):
        raise TimingError("schedstat version/header is malformed")
    _schedstat_uint(lines[1].split(b" ")[1], "schedstat timestamp")
    wanted: Optional[Dict[str, int]] = None
    previous_cpu = -1
    current_cpu: Optional[int] = None
    next_domain = 0
    mask_re = re.compile(rb"(?:[0-9a-f]{8},)*[0-9a-f]{8}\Z")
    for line in lines[2:]:
        fields = line.split(b" ")
        if not fields or any(field == b"" for field in fields):
            raise TimingError("schedstat line shape is malformed")
        cpu_match = re.fullmatch(rb"cpu(0|[1-9][0-9]*)", fields[0])
        if cpu_match is not None:
            if len(fields) != SCHEDSTAT_CPU_FIELD_COUNT + 1:
                raise TimingError("schedstat CPU line shape is malformed")
            cpu_id = _schedstat_uint(
                cpu_match.group(1), "schedstat CPU identifier")
            if cpu_id > MAX_CPU_ID or cpu_id <= previous_cpu:
                raise TimingError("schedstat CPU identifiers are duplicated or unordered")
            values = [
                _schedstat_uint(value, "schedstat CPU field")
                for value in fields[1:]
            ]
            # Version 15 retains fields 1 and 3-6.  Only the legacy field 2
            # is retired and emitted as zero.
            if values[1] != 0:
                raise TimingError("schedstat v15 retired CPU field is nonzero")
            if cpu_id == cpu:
                wanted = {
                    "runtime_ns": values[SCHEDSTAT_RUNTIME_FIELD - 1],
                    "pcount": values[SCHEDSTAT_PCOUNT_FIELD - 1],
                }
            previous_cpu = cpu_id
            current_cpu = cpu_id
            next_domain = 0
            continue
        domain_match = re.fullmatch(rb"domain(0|[1-9][0-9]*)", fields[0])
        if (domain_match is None or current_cpu is None or
                len(fields) != SCHEDSTAT_DOMAIN_FIELD_COUNT + 2 or
                mask_re.fullmatch(fields[1]) is None):
            raise TimingError("schedstat domain line shape is malformed")
        domain_id = _schedstat_uint(
            domain_match.group(1), "schedstat domain identifier")
        if domain_id != next_domain:
            raise TimingError("schedstat domains are duplicated or unordered")
        for value in fields[2:]:
            _schedstat_uint(value, "schedstat domain field")
        next_domain += 1
    if wanted is None:
        raise TimingError("CPU is missing from /proc/schedstat")
    return wanted


def schedstat_cpu(cpu: int) -> Dict[str, int]:
    try:
        with Path("/proc/schedstat").open("rb") as stream:
            raw = stream.read(MAX_JSON_EVIDENCE_BYTES + 1)
    except OSError as exc:
        raise TimingError("cannot read /proc/schedstat") from exc
    return parse_schedstat_cpu(raw, cpu)


def busy_ticks(ticks: Sequence[int]) -> int:
    idle = ticks[3] + (ticks[4] if len(ticks) > 4 else 0)
    return sum(ticks) - idle


def checked_busy_tick_delta(
    before: Sequence[int], after: Sequence[int], context: str,
) -> int:
    """Validate two /proc/stat CPU snapshots and return non-idle ticks."""
    if (not isinstance(before, (list, tuple)) or
            not isinstance(after, (list, tuple)) or len(before) < 5 or
            len(after) != len(before) or
            any(not isinstance(value, int) or isinstance(value, bool) or
                value < 0 for values in (before, after) for value in values) or
            any(right < left for left, right in zip(before, after))):
        raise TimingError("%s CPU tick receipt is malformed" % context)
    delta = busy_ticks(after) - busy_ticks(before)
    if delta < 0:
        raise TimingError("%s busy CPU tick delta is negative" % context)
    return delta


def require_tick_snapshot_order(
    earlier: Sequence[int], later: Sequence[int], context: str,
) -> None:
    if (not isinstance(earlier, (list, tuple)) or
            not isinstance(later, (list, tuple)) or
            len(earlier) != len(later) or
            any(right < left for left, right in zip(earlier, later))):
        raise TimingError("%s CPU tick order is malformed" % context)


def checked_schedstat_delta(
    before: Mapping[str, object], after: Mapping[str, object], context: str,
) -> Tuple[int, int]:
    """Validate schedstat snapshots and return runtime-ns/pcount deltas."""
    fields = {"runtime_ns", "pcount"}
    if (not isinstance(before, dict) or not isinstance(after, dict) or
            set(before) != fields or set(after) != fields or
            any(not isinstance(value, int) or isinstance(value, bool) or
                not 0 <= value <= MAX_UNSIGNED_COUNTER
                for snapshot in (before, after)
                for value in snapshot.values()) or
            any(after[field] < before[field] for field in fields)):
        raise TimingError("%s schedstat receipt is malformed" % context)
    return (
        int(after["runtime_ns"]) - int(before["runtime_ns"]),
        int(after["pcount"]) - int(before["pcount"]),
    )


def require_schedstat_snapshot_order(
    earlier: Mapping[str, object], later: Mapping[str, object], context: str,
) -> None:
    checked_schedstat_delta(earlier, later, context)


def require_exact_counter(value: object, expected: int, context: str) -> None:
    if (not isinstance(value, int) or isinstance(value, bool) or
            value != expected):
        raise TimingError("%s counter receipt mismatch" % context)


def _positive_duration_ns(value: object, context: str) -> int:
    if (not isinstance(value, int) or isinstance(value, bool) or
            not 0 < value <= MAX_UNSIGNED_COUNTER):
        raise TimingError("%s duration is malformed" % context)
    return value


def checked_attempt_duration_ns(
    start_ns: object, end_ns: object, recorded_duration_ns: object,
    context: str,
) -> int:
    if (not isinstance(start_ns, int) or isinstance(start_ns, bool) or
            not isinstance(end_ns, int) or isinstance(end_ns, bool) or
            not 0 <= start_ns < end_ns <= MAX_UNSIGNED_COUNTER):
        raise TimingError("%s duration is malformed" % context)
    duration_ns = end_ns - start_ns
    if (not isinstance(recorded_duration_ns, int) or
            isinstance(recorded_duration_ns, bool) or
            recorded_duration_ns != duration_ns):
        raise TimingError("%s duration receipt mismatch" % context)
    return duration_ns


def checked_campaign_interval(
    launch: Mapping[str, object],
) -> Tuple[float, float, float, int, int, int]:
    start_s = launch.get("start_monotonic_s")
    end_s = launch.get("end_monotonic_s")
    duration_s = launch.get("duration_s")
    start_ns = launch.get("start_monotonic_ns")
    end_ns = launch.get("end_monotonic_ns")
    duration_ns = launch.get("duration_ns")
    if (not isinstance(start_s, (int, float)) or isinstance(start_s, bool) or
            not isinstance(end_s, (int, float)) or isinstance(end_s, bool) or
            not isinstance(duration_s, (int, float)) or
            isinstance(duration_s, bool)):
        raise TimingError("launch monotonic interval is malformed")
    try:
        finite_seconds = all(math.isfinite(float(value)) for value in (
            start_s, end_s, duration_s))
    except (OverflowError, ValueError):
        finite_seconds = False
    if (not finite_seconds or end_s <= start_s or
            not isinstance(start_ns, int) or isinstance(start_ns, bool) or
            not isinstance(end_ns, int) or isinstance(end_ns, bool) or
            not isinstance(duration_ns, int) or isinstance(duration_ns, bool) or
            not 0 <= start_ns <= MAX_UNSIGNED_COUNTER or
            not 0 < end_ns <= MAX_UNSIGNED_COUNTER or end_ns <= start_ns or
            duration_ns != end_ns - start_ns or
            start_s != start_ns / 1_000_000_000 or
            end_s != end_ns / 1_000_000_000 or
            duration_s != end_s - start_s):
        raise TimingError("launch monotonic interval is malformed")
    return (float(start_s), float(end_s), float(duration_s),
            start_ns, end_ns, duration_ns)


def checked_preflight_limits(
    launch: Mapping[str, object], campaign_start_ns: int,
) -> Tuple[int, int, int]:
    if (not isinstance(campaign_start_ns, int) or
            isinstance(campaign_start_ns, bool) or campaign_start_ns <= 0):
        raise TimingError("launch preflight campaign boundary is malformed")
    duration_ns = checked_attempt_duration_ns(
        launch.get("preflight_start_monotonic_ns"),
        launch.get("preflight_end_monotonic_ns"),
        launch.get("preflight_duration_ns"), "launch preflight")
    if (duration_ns < SIBLING_PREFLIGHT_WINDOW_NS or
            launch["preflight_end_monotonic_ns"] > campaign_start_ns):
        raise TimingError("launch preflight interval is malformed")
    runtime_limit_ns = sibling_attempt_runtime_limit_ns(duration_ns)
    pcount_limit = sibling_attempt_pcount_limit(duration_ns)
    require_exact_counter(
        launch.get("preflight_sibling_sched_runtime_limit_ns"),
        runtime_limit_ns, "launch preflight sibling runtime limit")
    require_exact_counter(
        launch.get("preflight_sibling_sched_pcount_limit"), pcount_limit,
        "launch preflight sibling pcount limit")
    return duration_ns, runtime_limit_ns, pcount_limit


def sibling_campaign_busy_limit_ns(duration_ns: int) -> int:
    """Return the whole-tick floor of the 50-ppm /proc/stat allowance."""
    duration_ns = _positive_duration_ns(duration_ns, "sibling campaign")
    if CLOCK_TICKS_PER_SECOND <= 0:
        raise TimingError("sibling campaign tick rate is malformed")
    scaled = (duration_ns * CLOCK_TICKS_PER_SECOND *
              SIBLING_CAMPAIGN_MAX_BUSY_PPM) // (1_000_000_000 * 1_000_000)
    return max(SIBLING_CAMPAIGN_MIN_BUSY_TICKS, scaled)


def sibling_campaign_runtime_limit_ns(duration_ns: int) -> int:
    """Return the exact nanosecond floor of the 50-ppm schedstat bound."""
    duration_ns = _positive_duration_ns(duration_ns, "sibling campaign")
    return (duration_ns * SIBLING_CAMPAIGN_MAX_BUSY_PPM) // 1_000_000


def sibling_attempt_runtime_limit_ns(duration_ns: int) -> int:
    """Return the exact floor of the per-attempt 50-ppm schedstat bound."""
    duration_ns = _positive_duration_ns(duration_ns, "sibling attempt")
    return (duration_ns * SIBLING_SCHED_RUNTIME_MAX_PPM) // 1_000_000


def sibling_attempt_pcount_limit(duration_ns: int) -> int:
    """Allow one sibling schedule per started one-second attempt window."""
    duration_ns = _positive_duration_ns(duration_ns, "sibling attempt")
    return max(1, (duration_ns + SIBLING_SCHED_PCOUNT_WINDOW_NS - 1) //
               SIBLING_SCHED_PCOUNT_WINDOW_NS)


def _attempt_contaminations(
    parsed: ParsedOutput, sibling_busy: int, sibling_runtime_ns: int,
    sibling_pcount: int, duration_ns: int,
    target_irq_delta: Mapping[str, object],
) -> Tuple[str, ...]:
    for value in (sibling_busy, sibling_runtime_ns, sibling_pcount):
        if (not isinstance(value, int) or isinstance(value, bool) or value < 0):
            raise TimingError("attempt sibling counters are malformed")
    result = list(parsed.contaminations)
    if (not isinstance(target_irq_delta, dict) or
            set(target_irq_delta) != TARGET_IRQ_DELTA_FIELDS or
            not isinstance(target_irq_delta.get("contaminations"), list) or
            any(not isinstance(value, str) or not value
                for value in target_irq_delta["contaminations"])):
        raise TimingError("attempt IRQ classification is malformed")
    if target_irq_delta["contaminations"]:
        raise TimingError("target IRQ contamination is campaign-fatal")
    runtime_limit_ns = sibling_attempt_runtime_limit_ns(duration_ns)
    pcount_limit = sibling_attempt_pcount_limit(duration_ns)
    if sibling_busy > SIBLING_ACCEPTED_EXECUTION_MAX_BUSY_TICKS:
        result.append("sibling-busy:%d" % sibling_busy)
    if sibling_runtime_ns > runtime_limit_ns:
        result.append("sibling-sched-runtime-ns:%d" % sibling_runtime_ns)
    if sibling_pcount > pcount_limit:
        result.append("sibling-sched-pcount:%d" % sibling_pcount)
    return tuple(result)


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


def _thermal_pid(pid_file: Path) -> int:
    raw = stable_bytes(pid_file, attempts=3, max_bytes=64)
    if re.fullmatch(rb"[1-9][0-9]*\n", raw) is None:
        raise TimingError("thermal sampler PID file is not canonical")
    pid = int(raw[:-1])
    if not 1 < pid <= MAX_PROCESS_PID:
        raise TimingError("thermal sampler PID is outside its domain")
    return pid


def _thermal_reader(
    design: Mapping[str, object], pid_file: Path, thermal_csv: Path,
    pidfd: int, expected_pid: int, expected_start_ticks: int,
) -> Dict[str, object]:
    pid = _thermal_pid(pid_file)
    if pid != expected_pid or pidfd_has_exited(pidfd):
        raise TimingError("thermal sampler PID identity changed or exited")
    if process_start_ticks(pid) != expected_start_ticks:
        raise TimingError("thermal sampler start identity changed")
    tools = design["tools"]
    try:
        result = subprocess.run((
            str(tools["sudo"]["path"]), "-n", str(tools["fuser"]["path"]),
            "/dev/i2c-1", "/dev/i2c-2"), stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, timeout=5.0, check=False)
    except subprocess.TimeoutExpired as exc:
        raise TimingError("sudo fuser timed out during thermal verification") \
            from exc
    if result.returncode != 0:
        raise TimingError("cannot verify sole I2C reader with sudo fuser")
    reader_pids = sorted(set(int(value) for value in re.findall(rb"[0-9]+", result.stdout)))
    if reader_pids != [pid]:
        raise TimingError("thermal sampler is not the sole I2C reader")
    status = Path("/proc/%d/status" % pid).read_text(encoding="ascii")
    match = re.search(r"^Cpus_allowed_list:\s*(\S+)\s*$", status, re.MULTILINE)
    if match is None:
        raise TimingError("thermal sampler affinity is unavailable")
    allowed = parse_cpu_list(match.group(1))
    if len(allowed) != 1:
        raise TimingError("thermal sampler is not pinned to one CPU")
    if allowed[0] == int(design["controller_core"]):
        raise TimingError("thermal sampler shares the controller CPU")
    llc = set(int(value) for value in design["topology"]["llc_shared_cpus"])
    if set(allowed) & llc:
        raise TimingError("thermal sampler shares the timing LLC")
    root = Path(str(design["root"]))
    sampler = root / "frozen" / THERMAL_SAMPLER_NAME
    expected_tokens = [
        str(tools["python"]["path"]), "-I", str(sampler),
        "--csv", str(thermal_csv), "--pid-file", str(pid_file),
    ]
    tokens = process_tokens(pid)
    if tokens != expected_tokens:
        raise TimingError("thermal sampler command identity changed")
    boot_id = Path("/proc/sys/kernel/random/boot_id").read_text(
        encoding="ascii").strip()
    if BOOT_ID_RE.fullmatch(boot_id) is None:
        raise TimingError("kernel boot identity is malformed")
    if (process_start_ticks(pid) != expected_start_ticks or
            pidfd_has_exited(pidfd)):
        raise TimingError("thermal sampler changed during identity capture")
    sampler_sha256 = sha256_file(
        sampler, max_bytes=1024 * 1024, require_unique=True)
    if sampler_sha256 != design["immutable_files"].get(
            "frozen/" + THERMAL_SAMPLER_NAME):
        raise TimingError("thermal sampler source changed")
    return {
        "pid": pid, "start_ticks": expected_start_ticks,
        "boot_id": boot_id, "cpus_allowed_list": match.group(1),
        "argv": expected_tokens,
        "argv_sha256": sha256_bytes(
            b"\0".join(token.encode("utf-8") for token in expected_tokens)),
        "sampler_sha256": sampler_sha256,
        "thermal_csv": str(thermal_csv), "pid_file": str(pid_file),
    }


def _bind_thermal_reader(
    design: Mapping[str, object], pid_file: Path, thermal_csv: Path,
) -> Tuple[Dict[str, object], int]:
    pid = _thermal_pid(pid_file)
    start_ticks = process_start_ticks(pid)
    if start_ticks is None:
        raise TimingError("thermal sampler start identity is unavailable")
    try:
        pidfd = os.pidfd_open(pid, 0)
    except OSError as exc:
        raise TimingError("cannot bind thermal sampler pidfd") from exc
    try:
        if (process_start_ticks(pid) != start_ticks or
                pidfd_has_exited(pidfd)):
            raise TimingError(
                "thermal sampler changed during pidfd binding")
        return (
            _thermal_reader(
                design, pid_file, thermal_csv, pidfd, pid, start_ticks),
            pidfd,
        )
    except BaseException:
        os.close(pidfd)
        raise


def _validate_live_thermal_row(
    row: Mapping[str, str], now_s: float,
    baseline_edac_ce: int, baseline_edac_ue: int,
) -> None:
    monotonic = float(row["monotonic_s"])
    if now_s < monotonic or now_s - monotonic > MAX_THERMAL_MARGIN_S:
        raise TimingError("thermal sampler has no fresh live sample")
    if float(row["cpu_tctl_c"]) > MAX_CPU_TEMP_C:
        raise TimingError("live CPU temperature exceeded its gate")
    if max(float(row[field]) for field in DIMM_FIELDS) > MAX_DIMM_TEMP_C:
        raise TimingError("live DIMM temperature exceeded its gate")
    if int(row["dimm_read_errors"]) != 0:
        raise TimingError("live DIMM temperature read failed")
    if (int(row["edac_ce"]) != baseline_edac_ce or
            int(row["edac_ue"]) != baseline_edac_ue):
        raise TimingError("live EDAC counters changed")


def _parse_thermal(raw: bytes) -> Tuple[Tuple[Dict[str, str], ...], Tuple[bytes, ...]]:
    # The frozen sampler explicitly uses one canonical LF line ending.
    if (not raw or not raw.endswith(b"\n") or b"\r" in raw or
            b"\0" in raw or b'"' in raw):
        raise TimingError("thermal CSV is not canonical LF text")
    lines = tuple(raw.splitlines(keepends=True))
    if (any(not line.endswith(b"\n") for line in lines) or
            any(len(line[:-1].split(b",")) != len(THERMAL_FIELDS)
                for line in lines)):
        raise TimingError("thermal CSV field count changed")
    try:
        header = lines[0].decode("ascii")[:-1].split(",")
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
        if UTC_RE.fullmatch(row["utc"]) is None:
            raise TimingError("thermal UTC timestamp is malformed")
        try:
            datetime.strptime(row["utc"], "%Y-%m-%dT%H:%M:%S.%fZ")
        except ValueError as exc:
            raise TimingError("thermal UTC date is invalid") from exc
        monotonic = parse_decimal_float(
            row["monotonic_s"], "thermal monotonic timestamp", 0.0, 1e12)
        if monotonic <= previous:
            raise TimingError("thermal monotonic timestamps are not increasing")
        previous = monotonic
        parse_decimal_float(
            row["cpu_busy_pct"], "thermal CPU busy percentage", 0.0, 100.0)
        parse_decimal_float(
            row["cpu_avg_mhz"], "thermal CPU frequency", 0.001, 10000.0)
        for field in ("cpu_tctl_c", *DIMM_FIELDS):
            minimum = 0.001 if field == "cpu_tctl_c" else -39.999
            parse_decimal_float(
                row[field], "thermal " + field, minimum, 129.999)
        for field in ("load1", "load5", "load15"):
            parse_decimal_float(
                row[field], "thermal " + field, 0.0, 1000000.0)
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
            raw = stable_bytes(
                path, attempts=3, max_bytes=MAX_THERMAL_CSV_BYTES)
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
    if int(design["controller_core"]) in current["llc_shared_cpus"]:
        raise TimingError("controller CPU moved into the timing LLC")
    return current


def _execution_receipt(
    task: Mapping[str, object], slot: int, label: str, command: Sequence[str],
    attempt: int, started_utc: str, start_ns: int, end_ns: int,
    parsed: ParsedOutput, stderr: bytes, prior: Sequence[Mapping[str, object]],
    binary_sha256: str, process_identity: Mapping[str, object],
    cleanup_action: str, sibling_ticks_before: Sequence[int],
    sibling_ticks_after: Sequence[int],
    sibling_schedstat_before: Mapping[str, object],
    sibling_schedstat_after: Mapping[str, object],
    target_irq_snapshot_before: Mapping[str, object],
    target_irq_snapshot_after: Mapping[str, object],
    target_cpus: Sequence[int],
) -> Dict[str, object]:
    if (not isinstance(start_ns, int) or isinstance(start_ns, bool) or
            not isinstance(end_ns, int) or isinstance(end_ns, bool) or
            not 0 <= start_ns < end_ns <= MAX_UNSIGNED_COUNTER):
        raise TimingError("accepted execution duration is malformed")
    duration_ns = end_ns - start_ns
    if duration_ns < SIBLING_ACCEPTED_EXECUTION_MIN_DURATION_NS:
        raise TimingError("accepted execution is shorter than one USER_HZ tick")
    sibling_busy = checked_busy_tick_delta(
        sibling_ticks_before, sibling_ticks_after, "accepted execution")
    sibling_runtime_ns, sibling_pcount = checked_schedstat_delta(
        sibling_schedstat_before, sibling_schedstat_after,
        "accepted execution")
    sibling_runtime_limit_ns = sibling_attempt_runtime_limit_ns(duration_ns)
    sibling_pcount_limit = sibling_attempt_pcount_limit(duration_ns)
    target_irq_delta = checked_target_irq_delta(
        target_irq_snapshot_before, target_irq_snapshot_after, target_cpus,
        "accepted execution")
    require_target_irq_contained_interval(
        target_irq_snapshot_before, target_irq_snapshot_after,
        start_ns, end_ns, "accepted execution")
    if (sibling_busy > SIBLING_ACCEPTED_EXECUTION_MAX_BUSY_TICKS or
            sibling_runtime_ns > sibling_runtime_limit_ns or
            sibling_pcount > sibling_pcount_limit or
            target_irq_delta["contaminations"]):
        raise TimingError("accepted execution used the SMT sibling")
    receipt = sealed_record(
        "wirehair.wh2.grouped_commit_timing.execution_receipt.v4", {
            "job": task["job"], "task_id": task["task_id"],
            "task_sha256": sha256_bytes(canonical_json(task)),
            "outer_slot": slot, "outer_marker": OUTER_ORDER[slot],
            "binary_label": label, "binary_sha256": binary_sha256,
            "argv": list(command), "attempt": attempt,
            "started_utc": started_utc, "start_monotonic_ns": start_ns,
            "end_monotonic_ns": end_ns, "duration_ns": duration_ns,
            "stderr_sha256": sha256_bytes(stderr),
            "prior_contamination_receipts": list(prior),
            "process_identity": dict(process_identity),
            "cleanup_action": cleanup_action,
            "sibling_ticks_before": list(sibling_ticks_before),
            "sibling_ticks_after": list(sibling_ticks_after),
            "sibling_busy_ticks": sibling_busy,
            "sibling_schedstat_before": dict(sibling_schedstat_before),
            "sibling_schedstat_after": dict(sibling_schedstat_after),
            "sibling_sched_runtime_ns": sibling_runtime_ns,
            "sibling_sched_runtime_limit_ns": sibling_runtime_limit_ns,
            "sibling_sched_pcount": sibling_pcount,
            "sibling_sched_pcount_limit": sibling_pcount_limit,
            "target_irq_snapshot_before": dict(target_irq_snapshot_before),
            "target_irq_snapshot_after": dict(target_irq_snapshot_after),
            "target_irq_delta": target_irq_delta,
            **_receipt_summary(parsed),
        })
    if set(receipt) != EXECUTION_RECEIPT_FIELDS:
        raise TimingError("execution receipt constructor/schema drift")
    return receipt


def _save_contamination(
    root: Path, name: str, attempt: int, raw: bytes, stderr: bytes,
    parsed: ParsedOutput, command: Sequence[str],
    process_identity: Mapping[str, object], cleanup_action: str,
    start_ns: int, end_ns: int,
    sibling_ticks_before: Sequence[int], sibling_ticks_after: Sequence[int],
    sibling_schedstat_before: Mapping[str, object],
    sibling_schedstat_after: Mapping[str, object],
    target_irq_snapshot_before: Mapping[str, object],
    target_irq_snapshot_after: Mapping[str, object],
    target_cpus: Sequence[int],
) -> str:
    duration_ns = checked_attempt_duration_ns(
        start_ns, end_ns, end_ns - start_ns, "contaminated execution")
    if duration_ns < SIBLING_ACCEPTED_EXECUTION_MIN_DURATION_NS:
        raise TimingError("contaminated execution is shorter than one USER_HZ tick")
    sibling_busy = checked_busy_tick_delta(
        sibling_ticks_before, sibling_ticks_after, "contaminated execution")
    sibling_runtime_ns, sibling_pcount = checked_schedstat_delta(
        sibling_schedstat_before, sibling_schedstat_after,
        "contaminated execution")
    sibling_runtime_limit_ns = sibling_attempt_runtime_limit_ns(duration_ns)
    sibling_pcount_limit = sibling_attempt_pcount_limit(duration_ns)
    target_irq_delta = checked_target_irq_delta(
        target_irq_snapshot_before, target_irq_snapshot_after, target_cpus,
        "contaminated execution")
    require_target_irq_contained_interval(
        target_irq_snapshot_before, target_irq_snapshot_after,
        start_ns, end_ns, "contaminated execution")
    contaminations = _attempt_contaminations(
        parsed, sibling_busy, sibling_runtime_ns, sibling_pcount, duration_ns,
        target_irq_delta)
    if not contaminations:
        raise TimingError("cannot save an uncontaminated execution as contamination")
    prefix = "%s.attempt%d" % (name, attempt)
    raw_name = prefix + ".stdout"
    stderr_name = prefix + ".stderr"
    write_new(root / "contamination" / raw_name, raw)
    write_new(root / "contamination" / stderr_name, stderr)
    receipt = sealed_record(
        "wirehair.wh2.grouped_commit_timing.contamination_receipt.v4", {
            "name": name, "attempt": attempt, "argv": list(command),
            "process_identity": dict(process_identity),
            "cleanup_action": cleanup_action,
            "start_monotonic_ns": start_ns, "end_monotonic_ns": end_ns,
            "duration_ns": duration_ns,
            "stdout_name": raw_name, "stdout_sha256": sha256_bytes(raw),
            "stderr_name": stderr_name, "stderr_sha256": sha256_bytes(stderr),
            "contaminations": list(contaminations),
            "sibling_ticks_before": list(sibling_ticks_before),
            "sibling_ticks_after": list(sibling_ticks_after),
            "sibling_busy_ticks": sibling_busy,
            "sibling_schedstat_before": dict(sibling_schedstat_before),
            "sibling_schedstat_after": dict(sibling_schedstat_after),
            "sibling_sched_runtime_ns": sibling_runtime_ns,
            "sibling_sched_runtime_limit_ns": sibling_runtime_limit_ns,
            "sibling_sched_pcount": sibling_pcount,
            "sibling_sched_pcount_limit": sibling_pcount_limit,
            "target_irq_snapshot_before": dict(target_irq_snapshot_before),
            "target_irq_snapshot_after": dict(target_irq_snapshot_after),
            "target_irq_delta": target_irq_delta,
        })
    if set(receipt) != CONTAMINATION_RECEIPT_FIELDS:
        raise TimingError("contamination receipt constructor/schema drift")
    receipt_name = prefix + ".json"
    write_new(root / "contamination" / receipt_name, canonical_json(receipt))
    return sha256_file(root / "contamination" / receipt_name)


def _save_failure(
    root: Path, task: Mapping[str, object], slot: int, label: str,
    command: Sequence[str], attempt: int, started_utc: str, start_ns: int,
    end_ns: int, stdout: bytes, stderr: bytes, returncode: Optional[int],
    error: BaseException, binary_sha256: str,
    process_identity: Optional[Mapping[str, object]], cleanup_action: str,
    cleanup_error: Optional[str], sibling_ticks_before: Sequence[int],
    sibling_ticks_after: Sequence[int],
    sibling_schedstat_before: Mapping[str, object],
    sibling_schedstat_after: Mapping[str, object],
    target_irq_snapshot_before: Mapping[str, object],
    target_irq_snapshot_after: Mapping[str, object],
    target_cpus: Sequence[int],
) -> str:
    duration_ns = checked_attempt_duration_ns(
        start_ns, end_ns, end_ns - start_ns, "failed execution")
    sibling_busy = checked_busy_tick_delta(
        sibling_ticks_before, sibling_ticks_after, "failed execution")
    sibling_runtime_ns, sibling_pcount = checked_schedstat_delta(
        sibling_schedstat_before, sibling_schedstat_after,
        "failed execution")
    sibling_runtime_limit_ns = sibling_attempt_runtime_limit_ns(duration_ns)
    sibling_pcount_limit = sibling_attempt_pcount_limit(duration_ns)
    target_irq_delta = checked_target_irq_delta(
        target_irq_snapshot_before, target_irq_snapshot_after, target_cpus,
        "failed execution")
    require_target_irq_contained_interval(
        target_irq_snapshot_before, target_irq_snapshot_after,
        start_ns, end_ns, "failed execution")
    name = execution_name(task, slot, label)
    prefix = "%s.attempt%d" % (name, attempt)
    stdout_name = prefix + ".stdout"
    stderr_name = prefix + ".stderr"
    receipt_name = prefix + ".json"
    write_new(root / "failure" / stdout_name, stdout)
    write_new(root / "failure" / stderr_name, stderr)
    receipt = sealed_record(
        "wirehair.wh2.grouped_commit_timing.failure_receipt.v4", {
            "job": task["job"], "task_id": task["task_id"],
            "task_sha256": sha256_bytes(canonical_json(task)),
            "outer_slot": slot, "outer_marker": OUTER_ORDER[slot],
            "binary_label": label, "binary_sha256": binary_sha256,
            "argv": list(command), "attempt": attempt,
            "started_utc": started_utc, "start_monotonic_ns": start_ns,
            "end_monotonic_ns": end_ns, "duration_ns": duration_ns,
            "returncode": returncode, "error_type": type(error).__name__,
            "error_message": str(error), "stdout_name": stdout_name,
            "stdout_sha256": sha256_bytes(stdout),
            "stderr_name": stderr_name, "stderr_sha256": sha256_bytes(stderr),
            "process_identity": (
                dict(process_identity)
                if process_identity is not None else None),
            "cleanup_action": cleanup_action,
            "cleanup_error": cleanup_error,
            "sibling_ticks_before": list(sibling_ticks_before),
            "sibling_ticks_after": list(sibling_ticks_after),
            "sibling_busy_ticks": sibling_busy,
            "sibling_schedstat_before": dict(sibling_schedstat_before),
            "sibling_schedstat_after": dict(sibling_schedstat_after),
            "sibling_sched_runtime_ns": sibling_runtime_ns,
            "sibling_sched_runtime_limit_ns": sibling_runtime_limit_ns,
            "sibling_sched_pcount": sibling_pcount,
            "sibling_sched_pcount_limit": sibling_pcount_limit,
            "target_irq_snapshot_before": dict(target_irq_snapshot_before),
            "target_irq_snapshot_after": dict(target_irq_snapshot_after),
            "target_irq_delta": target_irq_delta,
        })
    if set(receipt) != FAILURE_RECEIPT_FIELDS:
        raise TimingError("failure receipt constructor/schema drift")
    write_new(root / "failure" / receipt_name, canonical_json(receipt))
    replay_failure_receipt(root, receipt_name, target_cpus)
    return sha256_file(root / "failure" / receipt_name)


def replay_failure_receipt(
    root: Path, receipt_name: str, target_cpus: Sequence[int],
) -> Dict[str, object]:
    """Replay a terminal failure's provenance, streams, and isolation proof."""
    if (not isinstance(receipt_name, str) or receipt_name in ("", ".", "..") or
            Path(receipt_name).name != receipt_name):
        raise TimingError("failure receipt name is malformed")
    receipt = load_canonical(
        root / "failure" / receipt_name, "failure receipt")
    verify_sealed_record(
        receipt, "wirehair.wh2.grouped_commit_timing.failure_receipt.v4",
        "failure receipt")
    if set(receipt) != FAILURE_RECEIPT_FIELDS:
        raise TimingError("failure receipt fields changed")
    design = _load_design(root)
    tasks = _load_tasks(root, design)
    expected_target_cpus = (
        int(design["core"]), int(design["topology"]["sibling"]),
        int(design["controller_core"]),
    )
    if _target_cpu_tuple(target_cpus, "failure replay") != expected_target_cpus:
        raise TimingError("failure replay target CPU scope changed")
    job = receipt.get("job")
    slot = receipt.get("outer_slot")
    attempt = receipt.get("attempt")
    if (not isinstance(job, int) or isinstance(job, bool) or
            not 0 <= job < len(tasks) or
            not isinstance(slot, int) or isinstance(slot, bool) or
            not 0 <= slot < len(OUTER_ORDER) or
            not isinstance(attempt, int) or isinstance(attempt, bool) or
            not 0 <= attempt < MAX_ENVIRONMENTAL_ATTEMPTS):
        raise TimingError("failure receipt coordinates are malformed")
    task = tasks[job]
    marker = OUTER_ORDER[slot]
    label = "base" if marker == "A" else "candidate"
    command = command_for(design, task, label)
    binary = root / "frozen" / BINARY_NAMES[label]
    prefix = "%s.attempt%d" % (execution_name(task, slot, label), attempt)
    immutable = design.get("immutable_files")
    if not isinstance(immutable, dict):
        raise TimingError("failure immutable-file ledger is malformed")
    expected_binary_sha256 = immutable.get("frozen/" + BINARY_NAMES[label])
    if (receipt_name != prefix + ".json" or
            receipt.get("task_id") != task["task_id"] or
            receipt.get("task_sha256") != sha256_bytes(canonical_json(task)) or
            receipt.get("outer_marker") != marker or
            receipt.get("binary_label") != label or
            receipt.get("binary_sha256") != expected_binary_sha256 or
            receipt.get("argv") != command or
            receipt.get("stdout_name") != prefix + ".stdout" or
            receipt.get("stderr_name") != prefix + ".stderr"):
        raise TimingError("failure receipt provenance binding mismatch")
    require_sha256(receipt.get("binary_sha256"), "failure binary hash")
    if sha256_file(binary) != expected_binary_sha256:
        raise TimingError("failure frozen binary changed")
    if (not isinstance(receipt.get("started_utc"), str) or
            UTC_RE.fullmatch(receipt["started_utc"]) is None or
            receipt.get("error_type") not in {"BoundCommandError", "TimingError"} or
            not isinstance(receipt.get("error_message"), str) or
            not 0 < len(receipt["error_message"]) <= 4096 or
            "\0" in receipt["error_message"] or
            (receipt.get("returncode") is not None and
             (not isinstance(receipt["returncode"], int) or
              isinstance(receipt["returncode"], bool) or
              not -255 <= receipt["returncode"] <= 255))):
        raise TimingError("failure receipt error metadata is malformed")
    _validate_failure_process_identity_receipt(
        receipt.get("process_identity"), receipt.get("cleanup_action"),
        receipt.get("cleanup_error"), command, binary)
    start_ns = receipt.get("start_monotonic_ns")
    end_ns = receipt.get("end_monotonic_ns")
    duration_ns = receipt.get("duration_ns")
    duration_ns = checked_attempt_duration_ns(
        start_ns, end_ns, duration_ns, "failure")
    for stream, maximum in (("stdout", MAX_GROUPED_OUTPUT_BYTES),
                            ("stderr", MAX_BENCHMARK_STDERR_BYTES)):
        name = receipt.get(stream + "_name")
        require_sha256(
            receipt.get(stream + "_sha256"), "failure %s hash" % stream)
        if (not isinstance(name, str) or
                sha256_bytes(stable_bytes(
                    root / "failure" / name, max_bytes=maximum)) !=
                receipt.get(stream + "_sha256")):
            raise TimingError("failure %s receipt does not replay" % stream)
    busy = checked_busy_tick_delta(
        receipt.get("sibling_ticks_before"),
        receipt.get("sibling_ticks_after"), "failed execution")
    require_exact_counter(
        receipt.get("sibling_busy_ticks"), busy,
        "failed execution busy ticks")
    runtime_ns, pcount = checked_schedstat_delta(
        receipt.get("sibling_schedstat_before"),
        receipt.get("sibling_schedstat_after"), "failed execution")
    require_exact_counter(
        receipt.get("sibling_sched_runtime_ns"), runtime_ns,
        "failed execution runtime")
    require_exact_counter(
        receipt.get("sibling_sched_pcount"), pcount,
        "failed execution pcount")
    require_exact_counter(
        receipt.get("sibling_sched_runtime_limit_ns"),
        sibling_attempt_runtime_limit_ns(duration_ns),
        "failed execution runtime limit")
    require_exact_counter(
        receipt.get("sibling_sched_pcount_limit"),
        sibling_attempt_pcount_limit(duration_ns),
        "failed execution pcount limit")
    validate_target_irq_delta(
        receipt.get("target_irq_delta"),
        receipt.get("target_irq_snapshot_before"),
        receipt.get("target_irq_snapshot_after"), expected_target_cpus,
        "failed execution")
    require_target_irq_contained_interval(
        receipt["target_irq_snapshot_before"],
        receipt["target_irq_snapshot_after"], start_ns, end_ns,
        "failed execution")
    return receipt


def run_campaign(args: argparse.Namespace) -> None:
    root = Path(args.result_dir).resolve()
    result_tmpfs_binding = capture_tmpfs_binding(
        root, "directory", "campaign result root")
    _require_external_prepare_anchor(root, args.expected_prepare_sha256)
    design = _load_design(root)
    controller_core = int(design["controller_core"])
    os.sched_setaffinity(0, {controller_core})
    if os.sched_getaffinity(0) != {controller_core}:
        raise TimingError("campaign controller could not pin outside the timing LLC")
    tasks = _load_tasks(root, design)
    _verify_immutable(root, design)
    _validate_prepare_receipt(root, design)
    _fresh_output_preflight(root, design)
    prepared_tree_tmpfs_bindings = capture_prepared_tree_tmpfs_bindings(
        root, design)
    if prepared_tree_tmpfs_bindings[0] != result_tmpfs_binding:
        raise TimingError("campaign result-root tmpfs binding changed at preflight")
    if _filler_pids():
        raise TimingError("load filler workers are still running")
    topology = _validate_topology_again(design)
    core = int(design["core"])
    sibling = int(topology["sibling"])
    thermal_csv = Path(args.thermal_csv).resolve()
    thermal_pid = Path(args.thermal_pid_file).resolve()
    thermal_csv_tmpfs_binding = capture_tmpfs_binding(
        thermal_csv, "regular", "thermal CSV")
    thermal_pid_tmpfs_binding = capture_tmpfs_binding(
        thermal_pid, "regular", "thermal PID file")
    if not (result_tmpfs_binding["device"] ==
            thermal_csv_tmpfs_binding["device"] ==
            thermal_pid_tmpfs_binding["device"]):
        raise TimingError("campaign artifacts and thermal evidence differ in tmpfs")
    thermal_reader, thermal_pidfd = _bind_thermal_reader(
        design, thermal_pid, thermal_csv)
    thermal_stat = thermal_csv.stat()
    thermal_identity = (thermal_stat.st_dev, thermal_stat.st_ino)
    latest = _latest_thermal_row(thermal_csv)
    baseline_edac_ce = int(latest["edac_ce"])
    baseline_edac_ue = int(latest["edac_ue"])
    _validate_live_thermal_row(
        latest, time.monotonic(), baseline_edac_ce, baseline_edac_ue)
    isolation_expected_cpus = (core, sibling, controller_core)
    target_cpus = (core, sibling, controller_core)
    preflight_start_monotonic_ns = time.monotonic_ns()
    before_core = cpu_ticks(core)
    preflight_sibling_schedstat_before = schedstat_cpu(sibling)
    before_sibling = cpu_ticks(sibling)
    preflight_target_irq_snapshot_before = capture_target_irq_snapshot(
        target_cpus)
    time.sleep(SIBLING_PREFLIGHT_WINDOW_NS / 1_000_000_000)
    preflight_target_irq_snapshot_after = capture_target_irq_snapshot(
        target_cpus)
    preflight_core_ticks_after = cpu_ticks(core)
    quiet_core = checked_busy_tick_delta(
        before_core, preflight_core_ticks_after, "preflight timing core")
    preflight_sibling_ticks_after = cpu_ticks(sibling)
    preflight_sibling_schedstat_after = schedstat_cpu(sibling)
    preflight_end_monotonic_ns = time.monotonic_ns()
    preflight_duration_ns = checked_attempt_duration_ns(
        preflight_start_monotonic_ns, preflight_end_monotonic_ns,
        preflight_end_monotonic_ns - preflight_start_monotonic_ns,
        "sibling preflight")
    if preflight_duration_ns < SIBLING_PREFLIGHT_WINDOW_NS:
        raise TimingError("sibling preflight was shorter than five seconds")
    preflight_runtime_limit_ns = sibling_attempt_runtime_limit_ns(
        preflight_duration_ns)
    preflight_pcount_limit = sibling_attempt_pcount_limit(
        preflight_duration_ns)
    preflight_target_irq_delta = checked_target_irq_delta(
        preflight_target_irq_snapshot_before,
        preflight_target_irq_snapshot_after, target_cpus, "sibling preflight")
    require_target_irq_contained_interval(
        preflight_target_irq_snapshot_before,
        preflight_target_irq_snapshot_after,
        preflight_start_monotonic_ns, preflight_end_monotonic_ns,
        "sibling preflight")
    quiet_sibling = checked_busy_tick_delta(
        before_sibling, preflight_sibling_ticks_after, "preflight sibling")
    quiet_sibling_runtime_ns, quiet_sibling_pcount = checked_schedstat_delta(
        preflight_sibling_schedstat_before,
        preflight_sibling_schedstat_after, "preflight sibling")
    if (quiet_core > SIBLING_PREFLIGHT_MAX_BUSY_TICKS or
            quiet_sibling > SIBLING_PREFLIGHT_MAX_BUSY_TICKS or
            quiet_sibling_runtime_ns > preflight_runtime_limit_ns or
            quiet_sibling_pcount > preflight_pcount_limit or
            preflight_target_irq_delta["contaminations"]):
        raise TimingError("timing core/sibling are not quiet before launch")
    runtime_isolation_snapshot_start = capture_runtime_isolation_snapshot(
        core, sibling, controller_core)
    campaign_target_irq_snapshot_before = capture_target_irq_snapshot(target_cpus)
    started_utc = utc_now()
    campaign_start_ns = time.monotonic_ns()
    campaign_start_s = campaign_start_ns / 1_000_000_000
    sibling_schedstat_before = schedstat_cpu(sibling)
    sibling_ticks_before = cpu_ticks(sibling)
    core_ticks_before = cpu_ticks(core)
    execution_hashes: List[Dict[str, object]] = []
    task_hashes: List[Dict[str, object]] = []
    retry_count = 0
    sibling_attempt_busy = 0
    sibling_attempt_runtime_ns = 0
    sibling_attempt_pcount = 0
    previous_target_irq_snapshot = campaign_target_irq_snapshot_before
    immutable = design["immutable_files"]
    cross_cache: Dict[
        Tuple[object, ...], Dict[str, Dict[str, object]]
    ] = {}
    for task in tasks:
        parsed_outputs: List[Tuple[str, ParsedOutput]] = []
        execution_records: List[Dict[str, object]] = []
        for slot, marker in enumerate(OUTER_ORDER):
            label = "base" if marker == "A" else "candidate"
            name = execution_name(task, slot, label)
            command = command_for(design, task, label)
            binary = root / "frozen" / BINARY_NAMES[label]
            python = Path(str(design["tools"]["python"]["path"]))
            binary_sha256 = str(
                immutable["frozen/" + BINARY_NAMES[label]])
            prior: List[Dict[str, object]] = []
            for attempt in range(MAX_ENVIRONMENTAL_ATTEMPTS):
                attempt_utc = utc_now()
                start_ns = time.monotonic_ns()
                attempt_sibling_schedstat_before = schedstat_cpu(sibling)
                attempt_sibling_before = cpu_ticks(sibling)
                attempt_target_irq_snapshot_before = \
                    capture_target_irq_snapshot(target_cpus)
                gap_target_irq_delta = checked_target_irq_delta(
                    previous_target_irq_snapshot,
                    attempt_target_irq_snapshot_before, target_cpus,
                    "pre-attempt campaign gap")
                if gap_target_irq_delta["contaminations"]:
                    raise TimingError(
                        "campaign gap accumulated target device IRQ activity: %s" %
                        ",".join(gap_target_irq_delta["contaminations"]))
                try:
                    result = run_bound_command(
                        command, binary, python, args.timeout_seconds)
                except BoundCommandError as exc:
                    attempt_target_irq_snapshot_after = \
                        capture_target_irq_snapshot(target_cpus)
                    attempt_sibling_after = cpu_ticks(sibling)
                    attempt_sibling_schedstat_after = schedstat_cpu(sibling)
                    end_ns = time.monotonic_ns()
                    _save_failure(
                        root, task, slot, label, command, attempt, attempt_utc,
                        start_ns, end_ns, exc.stdout, exc.stderr,
                        exc.returncode, exc, binary_sha256,
                        exc.process_identity, exc.cleanup_action,
                        exc.cleanup_error, attempt_sibling_before,
                        attempt_sibling_after,
                        attempt_sibling_schedstat_before,
                        attempt_sibling_schedstat_after,
                        attempt_target_irq_snapshot_before,
                        attempt_target_irq_snapshot_after, target_cpus)
                    raise
                attempt_target_irq_snapshot_after = \
                    capture_target_irq_snapshot(target_cpus)
                attempt_sibling_after = cpu_ticks(sibling)
                attempt_sibling_schedstat_after = schedstat_cpu(sibling)
                end_ns = time.monotonic_ns()
                attempt_sibling_delta = checked_busy_tick_delta(
                    attempt_sibling_before, attempt_sibling_after,
                    "timing attempt")
                attempt_sibling_runtime_ns, attempt_sibling_pcount = \
                    checked_schedstat_delta(
                        attempt_sibling_schedstat_before,
                        attempt_sibling_schedstat_after, "timing attempt")
                sibling_attempt_busy += attempt_sibling_delta
                sibling_attempt_runtime_ns += attempt_sibling_runtime_ns
                sibling_attempt_pcount += attempt_sibling_pcount
                attempt_target_irq_delta = checked_target_irq_delta(
                    attempt_target_irq_snapshot_before,
                    attempt_target_irq_snapshot_after, target_cpus,
                    "timing attempt")
                previous_target_irq_snapshot = attempt_target_irq_snapshot_after
                try:
                    if result.returncode != 0 or result.stderr:
                        raise TimingError(
                            "substantive timing command failure "
                            "exit=%d stderr=%r" %
                            (result.returncode, result.stderr[:1000]))
                    parsed = parse_grouped_output(
                        result.stdout, label, task,
                        int(design["evict_bytes"]), core)
                    if (end_ns - start_ns <
                            SIBLING_ACCEPTED_EXECUTION_MIN_DURATION_NS):
                        raise TimingError(
                            "execution is shorter than one USER_HZ tick")
                except TimingError as exc:
                    _save_failure(
                        root, task, slot, label, command, attempt, attempt_utc,
                        start_ns, end_ns, result.stdout, result.stderr,
                        result.returncode, exc, binary_sha256,
                        result.process_identity, result.cleanup_action, None,
                        attempt_sibling_before, attempt_sibling_after,
                        attempt_sibling_schedstat_before,
                        attempt_sibling_schedstat_after,
                        attempt_target_irq_snapshot_before,
                        attempt_target_irq_snapshot_after, target_cpus)
                    raise
                if attempt_target_irq_delta["contaminations"]:
                    exc = TimingError(
                        "timing attempt accumulated campaign-fatal target IRQ: %s" %
                        ",".join(attempt_target_irq_delta["contaminations"]))
                    _save_failure(
                        root, task, slot, label, command, attempt, attempt_utc,
                        start_ns, end_ns, result.stdout, result.stderr,
                        result.returncode, exc, binary_sha256,
                        result.process_identity, result.cleanup_action, None,
                        attempt_sibling_before, attempt_sibling_after,
                        attempt_sibling_schedstat_before,
                        attempt_sibling_schedstat_after,
                        attempt_target_irq_snapshot_before,
                        attempt_target_irq_snapshot_after, target_cpus)
                    raise exc
                contaminations = _attempt_contaminations(
                    parsed, attempt_sibling_delta,
                    attempt_sibling_runtime_ns, attempt_sibling_pcount,
                    end_ns - start_ns, attempt_target_irq_delta)
                if contaminations:
                    digest = _save_contamination(
                        root, name, attempt, result.stdout, result.stderr,
                        parsed, command, result.process_identity,
                        result.cleanup_action, start_ns, end_ns,
                        attempt_sibling_before, attempt_sibling_after,
                        attempt_sibling_schedstat_before,
                        attempt_sibling_schedstat_after,
                        attempt_target_irq_snapshot_before,
                        attempt_target_irq_snapshot_after, target_cpus)
                    prior.append({"attempt": attempt, "receipt_sha256": digest})
                    retry_count += 1
                    if attempt + 1 == MAX_ENVIRONMENTAL_ATTEMPTS:
                        raise ContaminationError(
                            "environmental retry limit exhausted: %s" %
                            ",".join(contaminations))
                    continue
                raw_path = root / "raw" / name
                stderr_path = root / "stderr" / (name + ".stderr")
                exit_path = root / "exit" / (name + ".exit")
                receipt_path = root / "receipts" / (name + ".json")
                write_new(raw_path, result.stdout)
                write_new(stderr_path, result.stderr)
                write_new(exit_path, b"0\n")
                receipt = _execution_receipt(
                    task, slot, label, command, attempt, attempt_utc, start_ns,
                    end_ns, parsed, result.stderr, prior,
                    binary_sha256, result.process_identity,
                    result.cleanup_action, attempt_sibling_before,
                    attempt_sibling_after,
                    attempt_sibling_schedstat_before,
                    attempt_sibling_schedstat_after,
                    attempt_target_irq_snapshot_before,
                    attempt_target_irq_snapshot_after, target_cpus)
                write_new(receipt_path, canonical_json(receipt))
                receipt_hash = sha256_file(receipt_path)
                record = {"name": name, "receipt_sha256": receipt_hash}
                execution_records.append(record)
                execution_hashes.append(record)
                parsed_outputs.append((label, parsed))
                break
        semantic = {item.semantic_sha256 for _label, item in parsed_outputs}
        traces = {item.preamble["trace_sha256"] for _label, item in parsed_outputs}
        work = {item.work_signature for _label, item in parsed_outputs}
        if len(semantic) != 1 or len(traces) != 1 or len(work) != 1:
            raise TimingError(
                "cross-binary executions changed trace or deterministic work")
        _register_cross_cache_identity(cross_cache, task, parsed_outputs[0][1])
        base_elapsed = sum(item.timed_elapsed_ns for label, item in parsed_outputs
                           if label == "base")
        candidate_elapsed = sum(
            item.timed_elapsed_ns for label, item in parsed_outputs
            if label == "candidate")
        if base_elapsed <= 0 or candidate_elapsed <= 0:
            raise TimingError("task timing sum is empty")
        base_phases = {
            field: sum(
                int(item.timed_phase_ns[field])
                for label, item in parsed_outputs if label == "base")
            for field in PHASE_FIELDS
        }
        candidate_phases = {
            field: sum(
                int(item.timed_phase_ns[field])
                for label, item in parsed_outputs if label == "candidate")
            for field in PHASE_FIELDS
        }
        if (base_phases[PRIMARY_PHASE_FIELD] <= 0 or
                candidate_phases[PRIMARY_PHASE_FIELD] <= 0):
            raise TimingError("task primary phase-timing sum is empty")
        phase_ratios = {
            field: _ratio_record_or_none(
                candidate_phases[field], base_phases[field])
            for field in PHASE_FIELDS
        }
        task_receipt = sealed_record(
            "wirehair.wh2.grouped_commit_timing.task_receipt.v1", {
                "job": task["job"], "task_id": task["task_id"],
                "task_sha256": sha256_bytes(canonical_json(task)),
                "outer_order": OUTER_ORDER,
                "execution_receipts": execution_records,
                "trace_sha256": next(iter(traces)),
                "semantic_sha256": next(iter(semantic)),
                "work_signature_sha256": sha256_bytes(canonical_json({
                    "fields": list(WORK_FIELDS), "values": list(next(iter(work)))})),
                "base_timed_elapsed_ns": base_elapsed,
                "candidate_timed_elapsed_ns": candidate_elapsed,
                "ratio": {"numerator": candidate_elapsed,
                          "denominator": base_elapsed},
                "base_timed_phase_ns": base_phases,
                "candidate_timed_phase_ns": candidate_phases,
                "phase_ratios": phase_ratios,
                "base_process_count": 4, "candidate_process_count": 4,
                "timed_rows_per_binary": 96,
            })
        task_path = root / "task_receipts" / (str(task["task_id"]) + ".json")
        write_new(task_path, canonical_json(task_receipt))
        task_hashes.append({
            "task_id": task["task_id"], "receipt_sha256": sha256_file(task_path)})
        if _filler_pids():
            raise TimingError("load filler worker appeared during timing")
        current_reader = _thermal_reader(
            design, thermal_pid, thermal_csv, thermal_pidfd,
            int(thermal_reader["pid"]), int(thermal_reader["start_ticks"]))
        if current_reader != thermal_reader:
            raise TimingError(
                "thermal reader identity/affinity changed during timing")
        current_thermal_stat = thermal_csv.stat()
        if (current_thermal_stat.st_dev, current_thermal_stat.st_ino) != \
                thermal_identity:
            raise TimingError("thermal CSV inode changed during timing")
        _validate_live_thermal_row(
            _latest_thermal_row(thermal_csv), time.monotonic(),
            baseline_edac_ce, baseline_edac_ue)
        if (int(task["job"]) + 1) % 6 == 0 or int(task["job"]) + 1 == len(tasks):
            print("progress=%d/%d retries=%d" %
                  (int(task["job"]) + 1, len(tasks), retry_count), flush=True)
    _validate_cross_cache_ledger(cross_cache)
    # Seal CPU isolation immediately after the final task.  Thermal identity
    # verification and interval collection may launch helpers or wait for the
    # next sampler row; neither belongs in the timed sibling-idle interval.
    core_ticks_after = cpu_ticks(core)
    sibling_ticks_after = cpu_ticks(sibling)
    sibling_schedstat_after = schedstat_cpu(sibling)
    campaign_end_ns = time.monotonic_ns()
    campaign_end_s = campaign_end_ns / 1_000_000_000
    campaign_target_irq_snapshot_after = capture_target_irq_snapshot(target_cpus)
    final_gap_target_irq_delta = checked_target_irq_delta(
        previous_target_irq_snapshot, campaign_target_irq_snapshot_after,
        target_cpus, "final campaign gap")
    sibling_busy = checked_busy_tick_delta(
        sibling_ticks_before, sibling_ticks_after, "campaign sibling")
    duration_ns = campaign_end_ns - campaign_start_ns
    duration_s = campaign_end_s - campaign_start_s
    sibling_busy_limit = sibling_campaign_busy_limit_ns(duration_ns)
    sibling_runtime_ns, sibling_pcount = checked_schedstat_delta(
        sibling_schedstat_before, sibling_schedstat_after,
        "campaign sibling")
    sibling_runtime_limit_ns = sibling_campaign_runtime_limit_ns(duration_ns)
    campaign_target_irq_delta = checked_target_irq_delta(
        campaign_target_irq_snapshot_before,
        campaign_target_irq_snapshot_after, target_cpus, "campaign")
    require_target_irq_bracketing_interval(
        campaign_target_irq_snapshot_before,
        campaign_target_irq_snapshot_after,
        campaign_start_ns, campaign_end_ns, "campaign")
    sibling_gap_busy = sibling_busy - sibling_attempt_busy
    sibling_gap_runtime_ns = sibling_runtime_ns - sibling_attempt_runtime_ns
    sibling_gap_pcount = sibling_pcount - sibling_attempt_pcount
    if (sibling_gap_busy < 0 or sibling_gap_runtime_ns < 0 or
            sibling_gap_pcount < 0):
        raise TimingError(
            "execution sibling accounting exceeds the campaign interval")
    if sibling_busy > sibling_busy_limit:
        raise TimingError(
            "SMT sibling accumulated %d busy ticks during timing (limit=%d)" %
            (sibling_busy, sibling_busy_limit))
    if sibling_runtime_ns > sibling_runtime_limit_ns:
        raise TimingError(
            "SMT sibling accumulated %d scheduled ns during timing "
            "(limit=%d)" % (sibling_runtime_ns, sibling_runtime_limit_ns))
    if (final_gap_target_irq_delta["contaminations"] or
            campaign_target_irq_delta["contaminations"]):
        raise TimingError(
            "campaign accumulated target device IRQ activity: %s" %
            ",".join(campaign_target_irq_delta["contaminations"]))
    runtime_isolation_snapshot_end = capture_runtime_isolation_snapshot(
        core, sibling, controller_core)
    validate_runtime_isolation_transition(
        runtime_isolation_snapshot_start, runtime_isolation_snapshot_end,
        isolation_expected_cpus)
    thermal_reader_end = _thermal_reader(
        design, thermal_pid, thermal_csv, thermal_pidfd,
        int(thermal_reader["pid"]), int(thermal_reader["start_ticks"]))
    if thermal_reader_end != thermal_reader:
        raise TimingError("thermal reader identity/affinity changed during timing")
    thermal_stat_end = thermal_csv.stat()
    if (thermal_stat_end.st_dev, thermal_stat_end.st_ino) != thermal_identity:
        raise TimingError("thermal CSV inode changed during timing")
    if (capture_tmpfs_binding(root, "directory", "campaign result root") !=
            result_tmpfs_binding or
            capture_tmpfs_binding(thermal_csv, "regular", "thermal CSV") !=
            thermal_csv_tmpfs_binding or
            capture_tmpfs_binding(thermal_pid, "regular", "thermal PID file") !=
            thermal_pid_tmpfs_binding):
        raise TimingError("campaign tmpfs binding changed during timing")
    validated_prepared_bindings = validate_prepared_tree_tmpfs_bindings(
        prepared_tree_tmpfs_bindings, root, design)
    if validated_prepared_bindings[0] != result_tmpfs_binding:
        raise TimingError("campaign prepared-tree tmpfs binding changed")
    thermal_raw, thermal_summary = collect_thermal_interval(
        thermal_csv, campaign_start_s, campaign_end_s)
    write_new(root / "thermal_interval.csv", thermal_raw)
    require_live_tmpfs_tree(
        root, int(result_tmpfs_binding["device"]), "completed campaign")
    launch = sealed_record(
        "wirehair.wh2.grouped_commit_timing.launch_receipt.v5", {
            "started_utc": started_utc, "ended_utc": utc_now(),
            "start_monotonic_s": campaign_start_s,
            "end_monotonic_s": campaign_end_s,
            "duration_s": duration_s,
            "start_monotonic_ns": campaign_start_ns,
            "end_monotonic_ns": campaign_end_ns,
            "duration_ns": duration_ns,
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
            "sibling_busy_limit_ticks": sibling_busy_limit,
            "sibling_attempt_busy_ticks": sibling_attempt_busy,
            "sibling_gap_busy_ticks": sibling_gap_busy,
            "preflight_quiet_core_ticks": quiet_core,
            "preflight_quiet_sibling_ticks": quiet_sibling,
            "preflight_core_ticks_before": list(before_core),
            "preflight_core_ticks_after": list(preflight_core_ticks_after),
            "preflight_sibling_ticks_before": list(before_sibling),
            "preflight_sibling_ticks_after":
                list(preflight_sibling_ticks_after),
            "preflight_sibling_schedstat_before":
                preflight_sibling_schedstat_before,
            "preflight_sibling_schedstat_after":
                preflight_sibling_schedstat_after,
            "preflight_start_monotonic_ns": preflight_start_monotonic_ns,
            "preflight_end_monotonic_ns": preflight_end_monotonic_ns,
            "preflight_duration_ns": preflight_duration_ns,
            "preflight_sibling_sched_runtime_ns": quiet_sibling_runtime_ns,
            "preflight_sibling_sched_runtime_limit_ns":
                preflight_runtime_limit_ns,
            "preflight_sibling_sched_pcount": quiet_sibling_pcount,
            "preflight_sibling_sched_pcount_limit": preflight_pcount_limit,
            "sibling_schedstat_before": sibling_schedstat_before,
            "sibling_schedstat_after": sibling_schedstat_after,
            "sibling_sched_runtime_ns": sibling_runtime_ns,
            "sibling_sched_runtime_limit_ns": sibling_runtime_limit_ns,
            "sibling_sched_pcount": sibling_pcount,
            "sibling_attempt_sched_runtime_ns":
                sibling_attempt_runtime_ns,
            "sibling_gap_sched_runtime_ns": sibling_gap_runtime_ns,
            "sibling_attempt_sched_pcount": sibling_attempt_pcount,
            "sibling_gap_sched_pcount": sibling_gap_pcount,
            "runtime_isolation_snapshot_start":
                runtime_isolation_snapshot_start,
            "runtime_isolation_snapshot_start_sha256":
                runtime_isolation_snapshot_sha256(
                    runtime_isolation_snapshot_start),
            "runtime_isolation_snapshot_end": runtime_isolation_snapshot_end,
            "runtime_isolation_snapshot_end_sha256":
                runtime_isolation_snapshot_sha256(
                    runtime_isolation_snapshot_end),
            "preflight_target_irq_snapshot_before":
                preflight_target_irq_snapshot_before,
            "preflight_target_irq_snapshot_after":
                preflight_target_irq_snapshot_after,
            "preflight_target_irq_delta": preflight_target_irq_delta,
            "campaign_target_irq_snapshot_before":
                campaign_target_irq_snapshot_before,
            "campaign_target_irq_snapshot_after":
                campaign_target_irq_snapshot_after,
            "campaign_target_irq_delta": campaign_target_irq_delta,
            "result_tmpfs_binding": result_tmpfs_binding,
            "thermal_csv_tmpfs_binding": thermal_csv_tmpfs_binding,
            "thermal_pid_tmpfs_binding": thermal_pid_tmpfs_binding,
            "prepared_tree_tmpfs_bindings": prepared_tree_tmpfs_bindings,
        })
    if set(launch) != LAUNCH_RECEIPT_FIELDS:
        raise TimingError("launch receipt constructor/schema drift")
    write_new(root / "launch_receipt.json", canonical_json(launch))
    os.close(thermal_pidfd)
    print(json.dumps({
        "task_count": len(tasks), "execution_count": len(execution_hashes),
        "retry_count": retry_count, "duration_s": launch["duration_s"],
        "thermal": thermal_summary,
        "launch_receipt_sha256": sha256_file(root / "launch_receipt.json"),
    }, sort_keys=True), flush=True)


def _validate_process_identity_receipt(
    identity: object, cleanup_action: object, command: Sequence[str],
    binary: Path, context: str,
) -> Dict[str, object]:
    binary_stat = binary.stat()
    try:
        binary_index = list(command).index(str(binary))
    except ValueError as exc:
        raise TimingError("%s command lacks the frozen binary" % context) from exc
    if (not isinstance(identity, dict) or
            set(identity) != {
                "pid", "start_ticks", "executable", "argv", "boot_id",
                "binding_observation"} or
            not isinstance(identity.get("pid"), int) or
            isinstance(identity.get("pid"), bool) or
            not 1 < identity["pid"] <= MAX_PROCESS_PID or
            not isinstance(identity.get("start_ticks"), int) or
            isinstance(identity.get("start_ticks"), bool) or
            identity["start_ticks"] <= 0 or
            identity.get("argv") != list(command)[binary_index:] or
            not isinstance(identity.get("boot_id"), str) or
            BOOT_ID_RE.fullmatch(identity["boot_id"]) is None or
            identity.get("binding_observation") not in {
                "double-proc-snapshot",
                "final-exec-handshake-terminal",
            } or cleanup_action != "exited_group_swept"):
        raise TimingError("%s process-identity receipt is malformed" % context)
    executable = identity.get("executable")
    if (not isinstance(executable, dict) or
            set(executable) != {"path", "device", "inode"} or
            executable.get("path") != str(binary) or
            executable.get("device") != binary_stat.st_dev or
            executable.get("inode") != binary_stat.st_ino):
        raise TimingError("%s executable binding mismatch" % context)
    return identity


def _validate_failure_process_identity_receipt(
    identity: object, cleanup_action: object, cleanup_error: object,
    command: Sequence[str], binary: Path,
) -> None:
    allowed_cleanup_actions = {
        "already_reaped", "exited_group_swept", "killed_group",
        "exited_group_swept_without_pidfd", "killed_group_without_pidfd",
        "cleanup_failed",
    }
    if cleanup_action not in allowed_cleanup_actions:
        raise TimingError("failure cleanup action is malformed")
    if (cleanup_error is not None and
            (not isinstance(cleanup_error, str) or
             not 0 < len(cleanup_error) <= 4096 or "\0" in cleanup_error)):
        raise TimingError("failure cleanup error is malformed")
    if cleanup_action == "cleanup_failed" and cleanup_error is None:
        raise TimingError("failure cleanup error is missing")
    if identity is None:
        if cleanup_action not in {
                "exited_group_swept_without_pidfd",
                "killed_group_without_pidfd", "cleanup_failed"}:
            raise TimingError("failure process identity is missing")
        return
    if (not isinstance(identity, dict) or
            set(identity) != {
                "pid", "start_ticks", "executable", "argv", "boot_id",
                "binding_observation"} or
            not isinstance(identity.get("pid"), int) or
            isinstance(identity.get("pid"), bool) or
            not 1 < identity["pid"] <= MAX_PROCESS_PID or
            not isinstance(identity.get("start_ticks"), int) or
            isinstance(identity.get("start_ticks"), bool) or
            identity["start_ticks"] <= 0 or
            not isinstance(identity.get("boot_id"), str) or
            BOOT_ID_RE.fullmatch(identity["boot_id"]) is None):
        raise TimingError("failure process identity is malformed")
    observation = identity.get("binding_observation")
    if observation == "direct-child-provisional":
        if identity.get("executable") is not None or identity.get("argv") is not None:
            raise TimingError("failure provisional process binding is malformed")
        return
    if observation not in {
            "double-proc-snapshot", "final-exec-handshake-terminal"}:
        raise TimingError("failure process binding observation is malformed")
    if cleanup_action in {
            "exited_group_swept_without_pidfd", "killed_group_without_pidfd"}:
        raise TimingError("failure final binding claims pidfd-less cleanup")
    try:
        binary_index = list(command).index(str(binary))
        binary_stat = binary.stat()
    except (ValueError, OSError) as exc:
        raise TimingError("failure frozen-binary binding is unavailable") from exc
    executable = identity.get("executable")
    if (identity.get("argv") != list(command)[binary_index:] or
            not isinstance(executable, dict) or
            set(executable) != {"path", "device", "inode"} or
            executable.get("path") != str(binary) or
            executable.get("device") != binary_stat.st_dev or
            executable.get("inode") != binary_stat.st_ino):
        raise TimingError("failure final process binding is malformed")


def _validate_execution_receipt(
    root: Path, design: Mapping[str, object], task: Mapping[str, object],
    slot: int, not_before_ns: int, not_after_ns: int,
) -> Tuple[
    Dict[str, object], ParsedOutput, int, int, int,
    Sequence[int], Sequence[int], Mapping[str, object], Mapping[str, object],
    Mapping[str, object], Mapping[str, object],
]:
    marker = OUTER_ORDER[slot]
    label = "base" if marker == "A" else "candidate"
    name = execution_name(task, slot, label)
    raw = stable_bytes(
        root / "raw" / name, max_bytes=MAX_GROUPED_OUTPUT_BYTES)
    stderr = stable_bytes(
        root / "stderr" / (name + ".stderr"),
        max_bytes=MAX_BENCHMARK_STDERR_BYTES)
    exit_raw = stable_bytes(
        root / "exit" / (name + ".exit"), max_bytes=32)
    if stderr or exit_raw != b"0\n":
        raise TimingError("execution stderr/exit artifact changed")
    parsed = parse_grouped_output(
        raw, label, task, int(design["evict_bytes"]), int(design["core"]))
    if parsed.contaminations:
        raise TimingError("accepted execution now parses as contaminated")
    path = root / "receipts" / (name + ".json")
    receipt = load_canonical(path, "execution receipt")
    verify_sealed_record(
        receipt, "wirehair.wh2.grouped_commit_timing.execution_receipt.v4",
        "execution receipt")
    if set(receipt) != EXECUTION_RECEIPT_FIELDS:
        raise TimingError("execution receipt fields changed")
    if (receipt.get("job") != task["job"] or
            receipt.get("task_id") != task["task_id"] or
            receipt.get("task_sha256") != sha256_bytes(canonical_json(task)) or
            receipt.get("outer_slot") != slot or
            receipt.get("outer_marker") != marker or
            receipt.get("binary_label") != label or
            receipt.get("argv") != command_for(design, task, label) or
            receipt.get("stderr_sha256") != sha256_bytes(stderr)):
        raise TimingError("execution receipt binding mismatch")
    identity = receipt.get("process_identity")
    cleanup_action = receipt.get("cleanup_action")
    binary = root / "frozen" / BINARY_NAMES[label]
    command = command_for(design, task, label)
    _validate_process_identity_receipt(
        identity, cleanup_action, command, binary, "execution")
    expected_summary = _receipt_summary(parsed)
    if any(receipt.get(key) != value for key, value in expected_summary.items()):
        raise TimingError("execution receipt parsed summary mismatch")
    if receipt.get("binary_sha256") != design["immutable_files"].get(
            "frozen/" + BINARY_NAMES[label]):
        raise TimingError("execution receipt binary binding mismatch")
    start = receipt.get("start_monotonic_ns")
    end = receipt.get("end_monotonic_ns")
    attempt = receipt.get("attempt")
    prior = receipt.get("prior_contamination_receipts")
    if (not isinstance(start, int) or isinstance(start, bool) or
            not isinstance(end, int) or isinstance(end, bool) or end <= start or
            start < not_before_ns or end > not_after_ns or
            receipt.get("duration_ns") != end - start or
            end - start < SIBLING_ACCEPTED_EXECUTION_MIN_DURATION_NS or
            not isinstance(attempt, int) or isinstance(attempt, bool) or
            not 0 <= attempt < MAX_ENVIRONMENTAL_ATTEMPTS or
            not isinstance(prior, list) or len(prior) != attempt):
        raise TimingError("execution timing/retry ledger mismatch")
    duration_ns = checked_attempt_duration_ns(
        start, end, receipt.get("duration_ns"), "accepted execution")
    if (not isinstance(receipt.get("started_utc"), str) or
            UTC_RE.fullmatch(receipt["started_utc"]) is None):
        raise TimingError("execution UTC receipt is malformed")
    accepted_sibling_before = receipt.get("sibling_ticks_before")
    accepted_sibling_after = receipt.get("sibling_ticks_after")
    accepted_sibling_busy = checked_busy_tick_delta(
        accepted_sibling_before, accepted_sibling_after,
        "accepted execution")
    recorded_accepted_sibling_busy = receipt.get("sibling_busy_ticks")
    if (not isinstance(recorded_accepted_sibling_busy, int) or
            isinstance(recorded_accepted_sibling_busy, bool) or
            accepted_sibling_busy != recorded_accepted_sibling_busy or
            accepted_sibling_busy >
            SIBLING_ACCEPTED_EXECUTION_MAX_BUSY_TICKS):
        raise TimingError("accepted execution sibling-idle receipt mismatch")
    accepted_schedstat_before = receipt.get("sibling_schedstat_before")
    accepted_schedstat_after = receipt.get("sibling_schedstat_after")
    accepted_sibling_runtime_ns, accepted_sibling_pcount = \
        checked_schedstat_delta(
            accepted_schedstat_before, accepted_schedstat_after,
            "accepted execution")
    require_exact_counter(
        receipt.get("sibling_sched_runtime_ns"),
        accepted_sibling_runtime_ns, "accepted execution runtime")
    require_exact_counter(
        receipt.get("sibling_sched_pcount"), accepted_sibling_pcount,
        "accepted execution pcount")
    accepted_runtime_limit_ns = sibling_attempt_runtime_limit_ns(duration_ns)
    accepted_pcount_limit = sibling_attempt_pcount_limit(duration_ns)
    require_exact_counter(
        receipt.get("sibling_sched_runtime_limit_ns"),
        accepted_runtime_limit_ns, "accepted execution runtime limit")
    require_exact_counter(
        receipt.get("sibling_sched_pcount_limit"), accepted_pcount_limit,
        "accepted execution pcount limit")
    if (accepted_sibling_runtime_ns > accepted_runtime_limit_ns or
            accepted_sibling_pcount > accepted_pcount_limit):
        raise TimingError("accepted execution schedstat receipt mismatch")
    target_cpus = (
        int(design["core"]), int(design["topology"]["sibling"]),
        int(design["controller_core"]))
    accepted_target_irq_delta = validate_target_irq_delta(
        receipt.get("target_irq_delta"),
        receipt.get("target_irq_snapshot_before"),
        receipt.get("target_irq_snapshot_after"), target_cpus,
        "accepted execution")
    require_target_irq_contained_interval(
        receipt["target_irq_snapshot_before"],
        receipt["target_irq_snapshot_after"], start, end,
        "accepted execution")
    if accepted_target_irq_delta["contaminations"]:
        raise TimingError("accepted execution target IRQ receipt mismatch")
    attempt_sibling_busy = accepted_sibling_busy
    attempt_sibling_runtime_ns = accepted_sibling_runtime_ns
    attempt_sibling_pcount = accepted_sibling_pcount
    previous_contamination_end: Optional[int] = None
    previous_sibling_after: Optional[Sequence[int]] = None
    first_sibling_before: Optional[Sequence[int]] = None
    previous_schedstat_after: Optional[Mapping[str, object]] = None
    first_schedstat_before: Optional[Mapping[str, object]] = None
    previous_target_irq_after: Optional[Mapping[str, object]] = None
    first_target_irq_before: Optional[Mapping[str, object]] = None
    for index, record in enumerate(prior):
        if (not isinstance(record, dict) or
                set(record) != {"attempt", "receipt_sha256"} or
                not isinstance(record.get("attempt"), int) or
                isinstance(record.get("attempt"), bool) or
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
            "wirehair.wh2.grouped_commit_timing.contamination_receipt.v4",
            "contamination receipt")
        expected_stdout_name = prefix + ".stdout"
        expected_stderr_name = prefix + ".stderr"
        if (set(contamination) != CONTAMINATION_RECEIPT_FIELDS or
                contamination.get("name") != name or
                not isinstance(contamination.get("attempt"), int) or
                isinstance(contamination.get("attempt"), bool) or
                contamination.get("attempt") != index or
                contamination.get("argv") != command_for(design, task, label) or
                contamination.get("stdout_name") != expected_stdout_name or
                contamination.get("stderr_name") != expected_stderr_name):
            raise TimingError("contamination receipt coordinate mismatch")
        _validate_process_identity_receipt(
            contamination.get("process_identity"),
            contamination.get("cleanup_action"), command, binary,
            "contaminated execution")
        for stream in ("stdout", "stderr"):
            artifact = root / "contamination" / str(contamination[stream + "_name"])
            if sha256_file(artifact) != contamination.get(stream + "_sha256"):
                raise TimingError("contamination artifact hash mismatch")
        contaminated = parse_grouped_output(
            stable_bytes(
                root / "contamination" /
                str(contamination["stdout_name"]),
                max_bytes=MAX_GROUPED_OUTPUT_BYTES),
            label, task, int(design["evict_bytes"]), int(design["core"]))
        contamination_stderr = stable_bytes(
            root / "contamination" / expected_stderr_name,
            max_bytes=MAX_BENCHMARK_STDERR_BYTES)
        contamination_start = contamination.get("start_monotonic_ns")
        contamination_end = contamination.get("end_monotonic_ns")
        contamination_duration_ns = checked_attempt_duration_ns(
            contamination_start, contamination_end,
            contamination.get("duration_ns"), "contaminated execution")
        contamination_sibling_before = contamination.get(
            "sibling_ticks_before")
        contamination_sibling_after = contamination.get(
            "sibling_ticks_after")
        contamination_sibling_busy = checked_busy_tick_delta(
            contamination_sibling_before, contamination_sibling_after,
            "contaminated execution")
        contamination_schedstat_before = contamination.get(
            "sibling_schedstat_before")
        contamination_schedstat_after = contamination.get(
            "sibling_schedstat_after")
        contamination_sibling_runtime_ns, contamination_sibling_pcount = \
            checked_schedstat_delta(
                contamination_schedstat_before, contamination_schedstat_after,
                "contaminated execution")
        contamination_target_irq_delta = validate_target_irq_delta(
            contamination.get("target_irq_delta"),
            contamination.get("target_irq_snapshot_before"),
            contamination.get("target_irq_snapshot_after"), target_cpus,
            "contaminated execution")
        require_target_irq_contained_interval(
            contamination["target_irq_snapshot_before"],
            contamination["target_irq_snapshot_after"],
            contamination_start, contamination_end,
            "contaminated execution")
        if previous_target_irq_after is not None:
            contamination_gap_delta = checked_target_irq_delta(
                previous_target_irq_after,
                contamination["target_irq_snapshot_before"], target_cpus,
                "inter-attempt contamination gap")
            if contamination_gap_delta["contaminations"]:
                raise TimingError(
                    "contamination gap contains target IRQ activity")
        expected_contaminations = _attempt_contaminations(
            contaminated, contamination_sibling_busy,
            contamination_sibling_runtime_ns,
            contamination_sibling_pcount, contamination_duration_ns,
            contamination_target_irq_delta)
        if (contamination_stderr or
                not isinstance(contamination_start, int) or
                isinstance(contamination_start, bool) or
                not isinstance(contamination_end, int) or
                isinstance(contamination_end, bool) or
                contamination_end <= contamination_start or
                contamination_duration_ns <
                SIBLING_ACCEPTED_EXECUTION_MIN_DURATION_NS or
                contamination_start < not_before_ns or
                (previous_contamination_end is not None and
                 contamination_start < previous_contamination_end) or
                contamination_end > start or
                not expected_contaminations or
                contaminated.semantic_sha256 != parsed.semantic_sha256 or
                contaminated.work_signature != parsed.work_signature or
                contaminated.preamble["trace_sha256"] !=
                parsed.preamble["trace_sha256"] or
                list(expected_contaminations) !=
                contamination.get("contaminations") or
                (not isinstance(contamination.get("sibling_busy_ticks"), int) or
                 isinstance(contamination.get("sibling_busy_ticks"), bool) or
                 contamination.get("sibling_busy_ticks") !=
                 contamination_sibling_busy) or
                (previous_sibling_after is not None and
                 (len(previous_sibling_after) !=
                  len(contamination_sibling_before) or
                  any(right < left for left, right in zip(
                      previous_sibling_after,
                      contamination_sibling_before)))) or
                (previous_schedstat_after is not None and
                 any(contamination_schedstat_before[field] <
                     previous_schedstat_after[field]
                     for field in ("runtime_ns", "pcount")))):
            raise TimingError("contamination classification changed on replay")
        require_exact_counter(
            contamination.get("sibling_sched_runtime_ns"),
            contamination_sibling_runtime_ns,
            "contaminated execution runtime")
        require_exact_counter(
            contamination.get("sibling_sched_pcount"),
            contamination_sibling_pcount, "contaminated execution pcount")
        require_exact_counter(
            contamination.get("sibling_sched_runtime_limit_ns"),
            sibling_attempt_runtime_limit_ns(contamination_duration_ns),
            "contaminated execution runtime limit")
        require_exact_counter(
            contamination.get("sibling_sched_pcount_limit"),
            sibling_attempt_pcount_limit(contamination_duration_ns),
            "contaminated execution pcount limit")
        previous_contamination_end = contamination_end
        if first_sibling_before is None:
            first_sibling_before = contamination_sibling_before
        if first_schedstat_before is None:
            first_schedstat_before = contamination_schedstat_before
        if first_target_irq_before is None:
            first_target_irq_before = contamination[
                "target_irq_snapshot_before"]
        previous_sibling_after = contamination_sibling_after
        previous_schedstat_after = contamination_schedstat_after
        previous_target_irq_after = contamination[
            "target_irq_snapshot_after"]
        attempt_sibling_busy += contamination_sibling_busy
        attempt_sibling_runtime_ns += contamination_sibling_runtime_ns
        attempt_sibling_pcount += contamination_sibling_pcount
    if previous_sibling_after is not None:
        require_tick_snapshot_order(
            previous_sibling_after, accepted_sibling_before,
            "accepted execution")
    if previous_schedstat_after is not None:
        require_schedstat_snapshot_order(
            previous_schedstat_after, accepted_schedstat_before,
            "accepted execution")
    if previous_target_irq_after is not None:
        accepted_gap_delta = checked_target_irq_delta(
            previous_target_irq_after,
            receipt["target_irq_snapshot_before"], target_cpus,
            "accepted execution gap")
        if accepted_gap_delta["contaminations"]:
            raise TimingError("accepted execution gap contains target IRQ activity")
    if first_sibling_before is None:
        first_sibling_before = accepted_sibling_before
    if first_schedstat_before is None:
        first_schedstat_before = accepted_schedstat_before
    if first_target_irq_before is None:
        first_target_irq_before = receipt["target_irq_snapshot_before"]
    return (
        receipt, parsed, attempt_sibling_busy, attempt_sibling_runtime_ns,
        attempt_sibling_pcount, first_sibling_before, accepted_sibling_after,
        first_schedstat_before, accepted_schedstat_after,
        first_target_irq_before, receipt["target_irq_snapshot_after"],
    )


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
    records: Sequence[Tuple[int, int]], domain: str,
    repetitions: int = BOOTSTRAP_REPS,
) -> Dict[str, object]:
    if (not records or repetitions < 100 or
            any(base < 0 or candidate < 0 for base, candidate in records) or
            not any(base > 0 for base, _candidate in records)):
        raise TimingError("bootstrap domain is empty or undersized")
    seed = int.from_bytes(hashlib.sha256(
        b"wirehair.wh2.grouped-commit.bootstrap.v1\0" + domain.encode("ascii")
    ).digest()[:8], "big")
    generator = random.Random(seed)
    values = []
    count = len(records)
    for _index in range(repetitions):
        for _retry in range(1000):
            base = 0
            candidate = 0
            for _draw in range(count):
                selected = records[generator.randrange(count)]
                base += selected[0]
                candidate += selected[1]
            if base > 0:
                break
        else:
            raise TimingError("bootstrap repeatedly sampled a zero denominator")
        values.append(candidate / base)
    values.sort()
    low = values[int(math.floor(0.025 * (repetitions - 1)))]
    high = values[int(math.ceil(0.975 * (repetitions - 1)))]
    return {
        "schema": "paired-task-ratio-of-sums-bootstrap-v1",
        "seed": seed, "repetitions": repetitions,
        "lower_95": format(low, ".12f"),
        "upper_95": format(high, ".12f"),
    }


def _summarize(
    records: Sequence[Mapping[str, object]], domain: str, metric: str,
) -> Dict[str, object]:
    if not records:
        raise TimingError("cannot summarize an empty record group")
    if any(int(record["base"]) < 0 or int(record["candidate"]) < 0
           for record in records):
        raise TimingError("timing summary contains a negative sample")
    base = sum(int(record["base"]) for record in records)
    candidate = sum(int(record["candidate"]) for record in records)
    if base <= 0:
        return {
            "metric": metric, "task_count": len(records),
            "base_" + metric: base, "candidate_" + metric: candidate,
            "ratio_available": False, "ratio_of_sums": None,
            "upper_median_task_ratio": None,
            "candidate_faster_tasks": 0, "candidate_slower_tasks": 0,
            "tied_tasks": 0, "ratio_defined_tasks": 0,
            "base_zero_tasks": len(records), "bootstrap": None,
        }
    ratios = [
        Fraction(int(record["candidate"]), int(record["base"]))
        for record in records if int(record["base"]) > 0
    ]
    faster = sum(1 for ratio in ratios if ratio < 1)
    slower = sum(1 for ratio in ratios if ratio > 1)
    ties = len(ratios) - faster - slower
    return {
        "metric": metric, "task_count": len(records),
        "base_" + metric: base, "candidate_" + metric: candidate,
        "ratio_available": True,
        "ratio_of_sums": _fraction_record(Fraction(candidate, base)),
        "upper_median_task_ratio": _fraction_record(_upper_median(ratios)),
        "candidate_faster_tasks": faster, "candidate_slower_tasks": slower,
        "tied_tasks": ties, "ratio_defined_tasks": len(ratios),
        "base_zero_tasks": len(records) - len(ratios),
        "bootstrap": _bootstrap(
            [(int(record["base"]), int(record["candidate"]))
             for record in records], domain),
    }


def _metric_breakdowns(
    records: Sequence[Mapping[str, object]], metric: str,
) -> Dict[str, object]:
    if metric != "elapsed_ns" and metric not in PHASE_FIELDS:
        raise TimingError("unknown timing summary metric")

    def selected(
        subset: Sequence[Mapping[str, object]],
    ) -> List[Dict[str, object]]:
        result = []
        for record in subset:
            timings = record.get("timings")
            if not isinstance(timings, dict):
                raise TimingError("result record timing ledger is malformed")
            pair = timings.get(metric)
            if not isinstance(pair, dict) or set(pair) != {"base", "candidate"}:
                raise TimingError("result record timing pair is malformed")
            result.append({
                "base": pair["base"], "candidate": pair["candidate"],
            })
        return result

    prefix = "metric:" + metric + ":"
    return {
        "overall": _summarize(
            selected(records), prefix + "overall", metric),
        "by_cache_state": {
            cache: _summarize(
                selected([
                    record for record in records
                    if record["cache_state"] == cache
                ]),
                prefix + "cache:" + cache,
                metric,
            )
            for cache in CACHE_STATES
        },
        "by_K": {
            str(K): _summarize(
                selected([record for record in records if record["K"] == K]),
                prefix + "K:%d" % K,
                metric,
            )
            for K in KS
        },
        "by_bb": {
            str(width): _summarize(
                selected([
                    record for record in records if record["bb"] == width
                ]),
                prefix + "bb:%d" % width,
                metric,
            )
            for width in WIDTHS
        },
        "by_schedule": {
            schedule: _summarize(
                selected([
                    record for record in records
                    if record["schedule"] == schedule
                ]),
                prefix + "schedule:" + schedule,
                metric,
            )
            for schedule, _seed in SCHEDULE_SEEDS
        },
    }


def replay_campaign(root: Path) -> Tuple[Dict[str, object], set[str]]:
    design = _load_design(root)
    _verify_directory_inventory(root)
    tasks = _load_tasks(root, design)
    _verify_immutable(root, design)
    _validate_prepare_receipt(root, design)
    launch = load_canonical(root / "launch_receipt.json", "launch receipt")
    verify_sealed_record(
        launch, "wirehair.wh2.grouped_commit_timing.launch_receipt.v5",
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
    result_tmpfs_binding = validate_tmpfs_binding(
        launch.get("result_tmpfs_binding"), root, "directory",
        "launch result root")
    prepared_tree_tmpfs_bindings = validate_prepared_tree_tmpfs_bindings(
        launch.get("prepared_tree_tmpfs_bindings"), root, design)
    if prepared_tree_tmpfs_bindings[0] != result_tmpfs_binding:
        raise TimingError("launch prepared-tree tmpfs binding mismatch")
    require_live_tmpfs_tree(
        root, int(result_tmpfs_binding["device"]), "replayed campaign")
    isolation_expected_cpus = (
        int(design["core"]), int(design["topology"]["sibling"]),
        int(design["controller_core"]))
    target_cpus = isolation_expected_cpus
    validate_runtime_isolation_transition(
        launch.get("runtime_isolation_snapshot_start"),
        launch.get("runtime_isolation_snapshot_end"),
        isolation_expected_cpus)
    for phase in ("start", "end"):
        hash_field = "runtime_isolation_snapshot_%s_sha256" % phase
        require_sha256(launch.get(hash_field), hash_field)
        if launch[hash_field] != runtime_isolation_snapshot_sha256(
                launch["runtime_isolation_snapshot_%s" % phase]):
            raise TimingError("runtime isolation snapshot hash mismatch")
    preflight_target_irq_delta = validate_target_irq_delta(
        launch.get("preflight_target_irq_delta"),
        launch.get("preflight_target_irq_snapshot_before"),
        launch.get("preflight_target_irq_snapshot_after"), target_cpus,
        "launch preflight")
    campaign_target_irq_delta = validate_target_irq_delta(
        launch.get("campaign_target_irq_delta"),
        launch.get("campaign_target_irq_snapshot_before"),
        launch.get("campaign_target_irq_snapshot_after"), target_cpus,
        "launch campaign")
    if (preflight_target_irq_delta["contaminations"] or
            campaign_target_irq_delta["contaminations"]):
        raise TimingError("launch target device IRQ receipt is contaminated")
    if (not isinstance(launch.get("retry_count"), int) or
            isinstance(launch.get("retry_count"), bool) or
            not 0 <= launch["retry_count"] <=
            len(tasks) * len(OUTER_ORDER) *
            (MAX_ENVIRONMENTAL_ATTEMPTS - 1)):
        raise TimingError("launch retry count is malformed")
    if (not isinstance(launch.get("started_utc"), str) or
            UTC_RE.fullmatch(launch["started_utc"]) is None or
            not isinstance(launch.get("ended_utc"), str) or
            UTC_RE.fullmatch(launch["ended_utc"]) is None):
        raise TimingError("launch UTC receipt is malformed")
    thermal_reader = launch.get("thermal_reader")
    if (not isinstance(thermal_reader, dict) or
            set(thermal_reader) != THERMAL_READER_FIELDS or
            not isinstance(thermal_reader.get("pid"), int) or
            isinstance(thermal_reader.get("pid"), bool) or
            not 1 < thermal_reader["pid"] <= MAX_PROCESS_PID or
            not isinstance(thermal_reader.get("start_ticks"), int) or
            isinstance(thermal_reader.get("start_ticks"), bool) or
            thermal_reader["start_ticks"] <= 0 or
            not isinstance(thermal_reader.get("boot_id"), str) or
            BOOT_ID_RE.fullmatch(thermal_reader["boot_id"]) is None or
            not isinstance(thermal_reader.get("cpus_allowed_list"), str) or
            len(parse_cpu_list(thermal_reader["cpus_allowed_list"])) != 1):
        raise TimingError("launch thermal-reader receipt is malformed")
    thermal_csv_tmpfs_binding = validate_tmpfs_binding(
        launch.get("thermal_csv_tmpfs_binding"),
        Path(str(thermal_reader["thermal_csv"])), "regular",
        "launch thermal CSV")
    thermal_pid_tmpfs_binding = validate_tmpfs_binding(
        launch.get("thermal_pid_tmpfs_binding"),
        Path(str(thermal_reader["pid_file"])), "regular",
        "launch thermal PID file")
    if not (result_tmpfs_binding["device"] ==
            thermal_csv_tmpfs_binding["device"] ==
            thermal_pid_tmpfs_binding["device"]):
        raise TimingError("launch tmpfs devices do not match")
    expected_sampler = (
        root / "frozen" / THERMAL_SAMPLER_NAME)
    argv = thermal_reader.get("argv")
    if (not isinstance(argv, list) or len(argv) != 7 or
            argv[:3] != [
                str(design["tools"]["python"]["path"]), "-I",
                str(expected_sampler)] or
            argv[3] != "--csv" or argv[4] != thermal_reader.get("thermal_csv") or
            argv[5] != "--pid-file" or argv[6] != thermal_reader.get("pid_file") or
            thermal_reader.get("sampler_sha256") !=
            design["immutable_files"].get(
                "frozen/" + THERMAL_SAMPLER_NAME) or
            thermal_reader.get("argv_sha256") != sha256_bytes(
                b"\0".join(str(token).encode("utf-8") for token in argv))):
        raise TimingError("launch thermal-reader command binding is malformed")
    require_sha256(
        thermal_reader.get("argv_sha256"), "thermal argv hash")
    require_sha256(
        thermal_reader.get("sampler_sha256"), "thermal sampler hash")
    for field in ("thermal_source_device", "thermal_source_inode"):
        value = launch.get(field)
        if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
            raise TimingError("launch thermal source identity is malformed")
    if (launch["thermal_source_device"] !=
            thermal_csv_tmpfs_binding["device"] or
            launch["thermal_source_inode"] !=
            thermal_csv_tmpfs_binding["inode"]):
        raise TimingError("launch thermal tmpfs/source identity mismatch")
    thermal_raw = stable_bytes(
        root / "thermal_interval.csv", max_bytes=MAX_THERMAL_CSV_BYTES)
    if sha256_bytes(thermal_raw) != launch.get("thermal_interval_sha256"):
        raise TimingError("thermal interval hash mismatch")
    (start_s, end_s, _duration_s, start_ns, end_ns, duration_ns) = \
        checked_campaign_interval(launch)
    require_target_irq_bracketing_interval(
        launch["campaign_target_irq_snapshot_before"],
        launch["campaign_target_irq_snapshot_after"],
        start_ns, end_ns, "launch campaign")
    thermal_summary = validate_sealed_thermal_interval(
        thermal_raw, float(start_s), float(end_s))
    if launch.get("thermal_summary") != thermal_summary:
        raise TimingError("thermal summary does not replay")
    checked_busy_tick_delta(
        launch.get("core_ticks_before"), launch.get("core_ticks_after"),
        "launch core")
    sibling_busy = checked_busy_tick_delta(
        launch.get("sibling_ticks_before"),
        launch.get("sibling_ticks_after"), "launch sibling")
    sibling_busy_limit = sibling_campaign_busy_limit_ns(duration_ns)
    recorded_sibling_busy = launch.get("sibling_busy_ticks")
    recorded_sibling_limit = launch.get("sibling_busy_limit_ticks")
    if (not isinstance(recorded_sibling_busy, int) or
            isinstance(recorded_sibling_busy, bool) or
            not isinstance(recorded_sibling_limit, int) or
            isinstance(recorded_sibling_limit, bool) or
            recorded_sibling_busy != sibling_busy or
            recorded_sibling_limit != sibling_busy_limit or
            sibling_busy > sibling_busy_limit):
        raise TimingError("launch sibling-idle receipt mismatch")
    sibling_schedstat_before = launch.get("sibling_schedstat_before")
    sibling_schedstat_after = launch.get("sibling_schedstat_after")
    sibling_runtime_ns, sibling_pcount = checked_schedstat_delta(
        sibling_schedstat_before, sibling_schedstat_after, "launch sibling")
    sibling_runtime_limit_ns = sibling_campaign_runtime_limit_ns(duration_ns)
    require_exact_counter(
        launch.get("sibling_sched_runtime_ns"), sibling_runtime_ns,
        "launch sibling runtime")
    require_exact_counter(
        launch.get("sibling_sched_pcount"), sibling_pcount,
        "launch sibling pcount")
    require_exact_counter(
        launch.get("sibling_sched_runtime_limit_ns"),
        sibling_runtime_limit_ns, "launch sibling runtime limit")
    if sibling_runtime_ns > sibling_runtime_limit_ns:
        raise TimingError("launch sibling schedstat receipt mismatch")
    (preflight_duration_ns, preflight_runtime_limit_ns,
     preflight_pcount_limit) = checked_preflight_limits(launch, start_ns)
    require_target_irq_contained_interval(
        launch["preflight_target_irq_snapshot_before"],
        launch["preflight_target_irq_snapshot_after"],
        launch["preflight_start_monotonic_ns"],
        launch["preflight_end_monotonic_ns"], "launch preflight")
    isolation_start = launch["runtime_isolation_snapshot_start"]
    isolation_end = launch["runtime_isolation_snapshot_end"]
    if (isolation_start["capture_start_monotonic_ns"] <
            launch["preflight_end_monotonic_ns"] or
            isolation_start["capture_end_monotonic_ns"] > start_ns or
            isolation_end["capture_start_monotonic_ns"] < end_ns):
        raise TimingError("runtime isolation capture order is malformed")
    preflight_quiet_core = checked_busy_tick_delta(
        launch.get("preflight_core_ticks_before"),
        launch.get("preflight_core_ticks_after"), "launch preflight core")
    preflight_quiet_sibling = checked_busy_tick_delta(
        launch.get("preflight_sibling_ticks_before"),
        launch.get("preflight_sibling_ticks_after"),
        "launch preflight sibling")
    require_exact_counter(
        launch.get("preflight_quiet_core_ticks"), preflight_quiet_core,
        "launch preflight core busy ticks")
    require_exact_counter(
        launch.get("preflight_quiet_sibling_ticks"), preflight_quiet_sibling,
        "launch preflight sibling busy ticks")
    if (preflight_quiet_core > SIBLING_PREFLIGHT_MAX_BUSY_TICKS or
            preflight_quiet_sibling > SIBLING_PREFLIGHT_MAX_BUSY_TICKS):
        raise TimingError("launch preflight quiet receipt mismatch")
    preflight_sibling_runtime_ns, preflight_sibling_pcount = \
        checked_schedstat_delta(
            launch.get("preflight_sibling_schedstat_before"),
            launch.get("preflight_sibling_schedstat_after"),
            "launch preflight sibling")
    require_exact_counter(
        launch.get("preflight_sibling_sched_runtime_ns"),
        preflight_sibling_runtime_ns, "launch preflight sibling runtime")
    require_exact_counter(
        launch.get("preflight_sibling_sched_pcount"),
        preflight_sibling_pcount, "launch preflight sibling pcount")
    if (preflight_sibling_runtime_ns > preflight_runtime_limit_ns or
            preflight_sibling_pcount > preflight_pcount_limit):
        raise TimingError("launch preflight schedstat receipt mismatch")
    require_schedstat_snapshot_order(
        launch["preflight_sibling_schedstat_after"],
        sibling_schedstat_before, "preflight-to-campaign sibling")
    require_tick_snapshot_order(
        launch["preflight_sibling_ticks_after"],
        launch["sibling_ticks_before"], "preflight-to-campaign sibling")
    require_tick_snapshot_order(
        launch["preflight_core_ticks_after"], launch["core_ticks_before"],
        "preflight-to-campaign timing core")
    for field in ("sibling_attempt_busy_ticks", "sibling_gap_busy_ticks"):
        value = launch.get(field)
        if (not isinstance(value, int) or isinstance(value, bool) or value < 0):
            raise TimingError("launch sibling accounting is malformed")
    for field in (
            "sibling_attempt_sched_runtime_ns",
            "sibling_gap_sched_runtime_ns", "sibling_attempt_sched_pcount",
            "sibling_gap_sched_pcount"):
        value = launch.get(field)
        if (not isinstance(value, int) or isinstance(value, bool) or value < 0):
            raise TimingError("launch sibling schedstat accounting is malformed")
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
    replay_sibling_attempt_busy = 0
    replay_sibling_attempt_runtime_ns = 0
    replay_sibling_attempt_pcount = 0
    previous_sibling_attempt_after = launch["sibling_ticks_before"]
    previous_sibling_schedstat_after = sibling_schedstat_before
    previous_target_irq_after = launch["campaign_target_irq_snapshot_before"]
    previous_execution_end = start_ns
    campaign_end_ns = end_ns
    for task in tasks:
        parsed_outputs = []
        execution_records = []
        for slot in range(len(OUTER_ORDER)):
            (receipt, parsed, attempt_sibling_busy,
             attempt_sibling_runtime_ns, attempt_sibling_pcount,
             attempt_sibling_before, attempt_sibling_after,
             attempt_sibling_schedstat_before,
             attempt_sibling_schedstat_after,
             attempt_target_irq_before,
             attempt_target_irq_after) = \
                _validate_execution_receipt(
                    root, design, task, slot, previous_execution_end,
                    campaign_end_ns)
            require_tick_snapshot_order(
                previous_sibling_attempt_after, attempt_sibling_before,
                "cross-execution sibling")
            require_schedstat_snapshot_order(
                previous_sibling_schedstat_after,
                attempt_sibling_schedstat_before,
                "cross-execution sibling")
            previous_sibling_attempt_after = attempt_sibling_after
            previous_sibling_schedstat_after = \
                attempt_sibling_schedstat_after
            target_gap_delta = checked_target_irq_delta(
                previous_target_irq_after, attempt_target_irq_before,
                target_cpus, "cross-execution target IRQ")
            if target_gap_delta["contaminations"]:
                raise TimingError(
                    "cross-execution gap contains target IRQ activity")
            previous_target_irq_after = attempt_target_irq_after
            replay_sibling_attempt_busy += attempt_sibling_busy
            replay_sibling_attempt_runtime_ns += attempt_sibling_runtime_ns
            replay_sibling_attempt_pcount += attempt_sibling_pcount
            previous_execution_end = int(receipt["end_monotonic_ns"])
            replay_retry_count += int(receipt["attempt"])
            label = str(receipt["binary_label"])
            name = execution_name(task, slot, label)
            receipt_path = root / "receipts" / (name + ".json")
            record = {"name": name, "receipt_sha256": sha256_file(receipt_path)}
            execution_ledger.append(record)
            execution_records.append(record)
            parsed_outputs.append((label, parsed))
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
        semantic = {parsed.semantic_sha256 for _label, parsed in parsed_outputs}
        traces = {parsed.preamble["trace_sha256"] for _label, parsed in parsed_outputs}
        work = {parsed.work_signature for _label, parsed in parsed_outputs}
        if len(semantic) != 1 or len(traces) != 1 or len(work) != 1:
            raise TimingError("replay found cross-binary semantic drift")
        _register_cross_cache_identity(cross_cache, task, parsed_outputs[0][1])
        base = sum(parsed.timed_elapsed_ns for label, parsed in parsed_outputs
                   if label == "base")
        candidate = sum(parsed.timed_elapsed_ns for label, parsed in parsed_outputs
                        if label == "candidate")
        base_phases = {
            field: sum(
                int(parsed.timed_phase_ns[field])
                for label, parsed in parsed_outputs if label == "base")
            for field in PHASE_FIELDS
        }
        candidate_phases = {
            field: sum(
                int(parsed.timed_phase_ns[field])
                for label, parsed in parsed_outputs if label == "candidate")
            for field in PHASE_FIELDS
        }
        if (base <= 0 or candidate <= 0 or
                base_phases[PRIMARY_PHASE_FIELD] <= 0 or
                candidate_phases[PRIMARY_PHASE_FIELD] <= 0):
            raise TimingError("replayed task contains an empty timing sum")
        phase_ratios = {
            field: _ratio_record_or_none(
                candidate_phases[field], base_phases[field])
            for field in PHASE_FIELDS
        }
        task_path = root / "task_receipts" / (str(task["task_id"]) + ".json")
        task_receipt = load_canonical(task_path, "task receipt")
        verify_sealed_record(
            task_receipt, "wirehair.wh2.grouped_commit_timing.task_receipt.v1",
            "task receipt")
        if set(task_receipt) != TASK_RECEIPT_FIELDS:
            raise TimingError("task receipt fields changed")
        expected_task = {
            "job": task["job"], "task_id": task["task_id"],
            "task_sha256": sha256_bytes(canonical_json(task)),
            "outer_order": OUTER_ORDER, "execution_receipts": execution_records,
            "trace_sha256": next(iter(traces)),
            "semantic_sha256": next(iter(semantic)),
            "work_signature_sha256": sha256_bytes(canonical_json({
                "fields": list(WORK_FIELDS), "values": list(next(iter(work)))})),
            "base_timed_elapsed_ns": base,
            "candidate_timed_elapsed_ns": candidate,
            "ratio": {"numerator": candidate, "denominator": base},
            "base_timed_phase_ns": base_phases,
            "candidate_timed_phase_ns": candidate_phases,
            "phase_ratios": phase_ratios,
            "base_process_count": 4, "candidate_process_count": 4,
            "timed_rows_per_binary": 96,
        }
        if any(task_receipt.get(key) != value for key, value in expected_task.items()):
            raise TimingError("task receipt replay mismatch")
        task_record = {
            "task_id": task["task_id"], "receipt_sha256": sha256_file(task_path)}
        task_ledger.append(task_record)
        expected_files.add("task_receipts/" + task_path.name)
        result_records.append({
            "K": task["K"], "bb": task["bb"],
            "schedule": task["schedule"], "cache_state": task["cache_state"],
            "timings": {
                "elapsed_ns": {"base": base, "candidate": candidate},
                **{
                    field: {
                        "base": base_phases[field],
                        "candidate": candidate_phases[field],
                    }
                    for field in PHASE_FIELDS
                },
            },
        })
    _validate_cross_cache_ledger(cross_cache)
    require_tick_snapshot_order(
        previous_sibling_attempt_after, launch["sibling_ticks_after"],
        "campaign-final sibling")
    require_schedstat_snapshot_order(
        previous_sibling_schedstat_after, sibling_schedstat_after,
        "campaign-final sibling")
    final_target_gap_delta = checked_target_irq_delta(
        previous_target_irq_after,
        launch["campaign_target_irq_snapshot_after"], target_cpus,
        "campaign-final target IRQ")
    if final_target_gap_delta["contaminations"]:
        raise TimingError("campaign-final gap contains target IRQ activity")
    expected_files.update(contamination_expected)
    if launch.get("execution_receipts") != execution_ledger or \
            launch.get("task_receipts") != task_ledger:
        raise TimingError("launch receipt ledger does not replay exactly")
    if launch.get("retry_count") != replay_retry_count:
        raise TimingError("launch retry count does not replay")
    replay_sibling_gap_busy = sibling_busy - replay_sibling_attempt_busy
    replay_sibling_gap_runtime_ns = \
        sibling_runtime_ns - replay_sibling_attempt_runtime_ns
    replay_sibling_gap_pcount = \
        sibling_pcount - replay_sibling_attempt_pcount
    if (replay_sibling_gap_busy < 0 or
            replay_sibling_gap_runtime_ns < 0 or
            replay_sibling_gap_pcount < 0 or
            launch.get("sibling_attempt_busy_ticks") !=
            replay_sibling_attempt_busy or
            launch.get("sibling_gap_busy_ticks") != replay_sibling_gap_busy or
            launch.get("sibling_attempt_sched_runtime_ns") !=
            replay_sibling_attempt_runtime_ns or
            launch.get("sibling_gap_sched_runtime_ns") !=
            replay_sibling_gap_runtime_ns or
            launch.get("sibling_attempt_sched_pcount") !=
            replay_sibling_attempt_pcount or
            launch.get("sibling_gap_sched_pcount") !=
            replay_sibling_gap_pcount):
        raise TimingError("launch sibling accounting does not replay")
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
    elapsed = _metric_breakdowns(result_records, "elapsed_ns")
    phase_metrics = {
        field: {
            "role": "primary" if field == PRIMARY_PHASE_FIELD
            else "negative-control",
            **_metric_breakdowns(result_records, field),
        }
        for field in PHASE_FIELDS
    }
    summary_payload: Dict[str, object] = {
        "base_commit": BASE_COMMIT, "candidate_commit": CANDIDATE_COMMIT,
        "codec_lineage": design["codec_lineage"],
        "measurement_overlay": design["measurement_overlay"],
        "architecture": ARCHITECTURE, "task_count": len(tasks),
        "execution_count": len(execution_ledger),
        "timed_rows_per_binary": len(tasks) * 96,
        "work_and_recovery_exact": True, "trace_identity_exact": True,
        "cold_warm_graph_trace_work_exact": True,
        "timing_promotional": True,
        "sibling_idle": {
            "policy": design["sibling_idle_policy"],
            "accepted_execution_dynamic_limits_replayed": True,
            "preflight": {
                "duration_ns": preflight_duration_ns,
                "timing_core_busy_ticks": preflight_quiet_core,
                "sibling_busy_ticks": preflight_quiet_sibling,
                "sibling_scheduled_runtime_ns":
                    preflight_sibling_runtime_ns,
                "sibling_scheduled_runtime_limit_ns":
                    preflight_runtime_limit_ns,
                "sibling_pcount": preflight_sibling_pcount,
                "sibling_pcount_limit": preflight_pcount_limit,
                "target_irq_classifications":
                    preflight_target_irq_delta["classifications"],
            },
            "busy_ticks": sibling_busy,
            "busy_limit_ticks": sibling_busy_limit,
            "attempt_busy_ticks": replay_sibling_attempt_busy,
            "gap_busy_ticks": replay_sibling_gap_busy,
            "scheduled_runtime_ns": sibling_runtime_ns,
            "scheduled_runtime_limit_ns": sibling_runtime_limit_ns,
            "pcount": sibling_pcount,
            "attempt_scheduled_runtime_ns":
                replay_sibling_attempt_runtime_ns,
            "gap_scheduled_runtime_ns": replay_sibling_gap_runtime_ns,
            "attempt_pcount": replay_sibling_attempt_pcount,
            "gap_pcount": replay_sibling_gap_pcount,
            "campaign_target_irq_classifications":
                campaign_target_irq_delta["classifications"],
        },
        "runtime_isolation": {
            "expected_cpus": list(sorted(isolation_expected_cpus)),
            "self_cgroup": isolation_start["self_cgroup"],
            "start_snapshot_sha256":
                launch["runtime_isolation_snapshot_start_sha256"],
            "end_snapshot_sha256":
                launch["runtime_isolation_snapshot_end_sha256"],
            "start_numeric_irq_count": len(
                isolation_start["irq_effective_affinities"]),
            "end_numeric_irq_count": len(
                isolation_end["irq_effective_affinities"]),
            "guarded_irq30_unchanged":
                isolation_start["irq30_exception"] ==
                isolation_end["irq30_exception"],
            "managed_nvme_whitelist_unchanged":
                isolation_start["managed_nvme_exceptions"] ==
                isolation_end["managed_nvme_exceptions"],
        },
        "thermal_interval_sha256": sha256_bytes(thermal_raw),
        "thermal_summary": thermal_summary,
        "overall": elapsed["overall"],
        "by_cache_state": elapsed["by_cache_state"],
        "by_K": elapsed["by_K"],
        "by_bb": elapsed["by_bb"],
        "by_schedule": elapsed["by_schedule"],
        "phase_timing": {
            "primary": PRIMARY_PHASE_FIELD,
            "negative_controls": list(NEGATIVE_CONTROL_PHASE_FIELDS),
            "metrics": phase_metrics,
        },
    }
    return summary_payload, expected_files


def _data_manifest(root: Path, expected: Iterable[str]) -> List[Dict[str, object]]:
    records = []
    for relative in sorted(set(expected)):
        path = root / relative
        if path.is_symlink() or not path.is_file():
            raise TimingError("data manifest path is not a regular file: %s" % relative)
        raw = stable_bytes(path, max_bytes=MAX_EVIDENCE_FILE_BYTES)
        records.append({
            "path": relative, "size": len(raw),
            "sha256": sha256_bytes(raw),
        })
    return records


def _publish_or_verify_reduction_artifact(path: Path, data: bytes) -> str:
    """Publish once, or accept only an exact immutable crash-recovery artifact."""
    try:
        write_new(path, data)
        return "published"
    except TimingError as publication_error:
        try:
            before = os.lstat(path)
            existing = stable_bytes(
                path, attempts=3, max_bytes=len(data), require_unique=True)
            after = os.lstat(path)
        except (OSError, TimingError):
            raise TimingError(
                "cannot publish or resume reduction artifact: %s" % path.name
            ) from publication_error
        identity_before = (
            before.st_dev, before.st_ino, before.st_mode, before.st_nlink,
            before.st_uid, before.st_gid, before.st_size, before.st_mtime_ns,
            before.st_ctime_ns,
        )
        identity_after = (
            after.st_dev, after.st_ino, after.st_mode, after.st_nlink,
            after.st_uid, after.st_gid, after.st_size, after.st_mtime_ns,
            after.st_ctime_ns,
        )
        if (not stat.S_ISREG(before.st_mode) or before.st_nlink != 1 or
                before.st_uid != os.geteuid() or
                stat.S_IMODE(before.st_mode) != 0o444 or
                identity_before != identity_after):
            raise TimingError(
                "existing reduction artifact is not immutable: %s" %
                path.name) from publication_error
        if existing != data:
            raise TimingError(
                "existing reduction artifact does not reproduce: %s" %
                path.name) from publication_error
        return "verified_existing"


def reduce_campaign(args: argparse.Namespace) -> None:
    root = Path(args.result_dir).resolve()
    _require_external_prepare_anchor(root, args.expected_prepare_sha256)
    payload, expected = replay_campaign(root)
    summary = sealed_record(
        "wirehair.wh2.grouped_commit_timing.validated_summary.v2", payload)
    summary_path = root / "validated_summary.json"
    publication = {
        "validated_summary.json": _publish_or_verify_reduction_artifact(
            summary_path, canonical_json(summary)),
    }
    expected.add("validated_summary.json")
    manifest = _data_manifest(root, expected)
    manifest_path = root / "data_manifest.json"
    publication["data_manifest.json"] = \
        _publish_or_verify_reduction_artifact(
            manifest_path, canonical_json({"files": manifest}))
    publication["data_manifest.sha256"] = \
        _publish_or_verify_reduction_artifact(
            root / "data_manifest.sha256",
            (sha256_file(manifest_path) + "\n").encode("ascii"))
    print(json.dumps({
        "summary_sha256": sha256_file(summary_path),
        "data_manifest_sha256": sha256_file(manifest_path),
        "overall_ratio": summary["overall"]["ratio_of_sums"]["decimal"],
        "residual_ratio": summary["phase_timing"]["metrics"][
            PRIMARY_PHASE_FIELD]["overall"]["ratio_of_sums"]["decimal"],
        "cold_ratio": summary["by_cache_state"]["cold"]["ratio_of_sums"]["decimal"],
        "warm_ratio": summary["by_cache_state"]["warm"]["ratio_of_sums"]["decimal"],
        "work_and_recovery_exact": True,
        "publication": publication,
    }, sort_keys=True))


def verify_campaign(args: argparse.Namespace) -> None:
    root = Path(args.result_dir).resolve()
    _require_external_prepare_anchor(root, args.expected_prepare_sha256)
    payload, expected = replay_campaign(root)
    summary = load_canonical(root / "validated_summary.json", "validated summary")
    verify_sealed_record(
        summary, "wirehair.wh2.grouped_commit_timing.validated_summary.v2",
        "validated summary")
    expected_summary = sealed_record(
        "wirehair.wh2.grouped_commit_timing.validated_summary.v2", payload)
    if summary != expected_summary:
        raise TimingError("validated summary does not reproduce")
    expected.add("validated_summary.json")
    manifest = load_canonical(root / "data_manifest.json", "data manifest")
    if set(manifest) != {"files"} or manifest["files"] != _data_manifest(root, expected):
        raise TimingError("data manifest does not reproduce")
    try:
        sidecar = stable_bytes(
            root / "data_manifest.sha256", max_bytes=128).decode(
                "ascii", "strict")
    except UnicodeDecodeError as exc:
        raise TimingError("data manifest sidecar is not ASCII") from exc
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
    prepare.add_argument("--thermal-sampler", required=True)
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
    run.add_argument("--thermal-csv", required=True)
    run.add_argument("--thermal-pid-file", required=True)
    run.add_argument("--expected-prepare-sha256", required=True)
    run.add_argument("--timeout-seconds", type=int, default=1800)
    run.set_defaults(function=run_campaign)
    reduce = subparsers.add_parser("reduce", help="strictly replay and summarize")
    reduce.add_argument("--result-dir", required=True)
    reduce.add_argument("--expected-prepare-sha256", required=True)
    reduce.set_defaults(function=reduce_campaign)
    verify = subparsers.add_parser("verify", help="replay a reduced campaign")
    verify.add_argument("--result-dir", required=True)
    verify.add_argument("--expected-prepare-sha256", required=True)
    verify.set_defaults(function=verify_campaign)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = make_parser()
    args = parser.parse_args(argv)
    try:
        if getattr(args, "timeout_seconds", 1) <= 0:
            raise TimingError("timeout must be positive")
        result = Path(args.result_dir).resolve()
        lock_name = ".wirehair-wh2-timing-%s.lock" % sha256_bytes(
            str(result).encode("utf-8"))[:24]
        with evidence_io.nonblocking_global_flock(
                result.parent / lock_name, error_type=TimingError):
            args.function(args)
    except (OSError, ValueError, TimingError) as exc:
        print("grouped commit timing failed: %s" % exc, file=sys.stderr, flush=True)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
