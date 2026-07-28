#!/usr/bin/env python3
"""Shared, fail-closed access to the native WH2 peel benchmark.

The native codec is the authority for the shipped peel distribution.  Do not
reconstruct it from an analytic soliton formula here: the production generator
has a fixed degree-one prefix through K=2048 and no degree-one prefix above it,
and its fixed-point thresholds are not the same as ``1 / K``.
"""

import csv
import copy
import hashlib
import json
import math
import os
import platform
import re
import stat
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, replace


PEEL_TABLE_SCHEMA = "wirehair-v2-peel-table/v4"
NATIVE_PMF_PROTOCOL = "wirehair-v2-bench:peelpmf:dispatch-v1:v2"
NATIVE_COMPARE_PROTOCOL = "wirehair-v2-bench:compare:wh2-target:v4"
NATIVE_PAIRED_PROTOCOL = "wirehair-v2-bench:peeltiming:dispatch-v1:v1"
SEED_DERIVATION_PROTOCOL = "blake2b64:wh2-peel-v2:split-construction-loss"
COMPLETION_REGIME_PROTOCOL = "wirehair-v2-dispatch-v1-completion/v2"
TARGET_PROFILE = "dispatch-v1"
TARGET_SEED_POLICY = "raw"
TARGET_RECEIPT_SEED_POLICY = "raw-uniform-v1"
TARGET_CONTRACT = {
    "target_profile": TARGET_PROFILE,
    "contract_id": "a98c37c23ee7feae",
    "contract_sha256":
        "a98c37c23ee7feae4171ff10627f660f705db6b7aae9268f85617ce86396583c",
    "precode_contract": 5,
    "packet_contract": 4,
    "architecture": "smallband100-d4",
    "completion": "mixed",
    "dense_rows": 4,
    "heavy_rows": 12,
    "gf256_rows": 10,
    "gf16_rows": 2,
    "period": 244,
    "geometry": "frozen",
    "residue_schedule": "constant",
    "residue_skew": 0,
    "mix_count": 3,
    "pmf_encoding": "wirehair-v2-peel-spec-v1",
}
TARGET_MEASUREMENT_POLICY = {
    "seed_policy": TARGET_RECEIPT_SEED_POLICY,
    "attempt_policy": "single",
    "seed_attempt": "0",
    "seed_tables": "none",
    "fixups": "none",
    "grouped_gf256_rows": "0",
    "independent_extension": "0",
    "heavy_family": "periodic",
    "dense_identity_corner": "0",
    "dense_two_anchor": "0",
    "packet_seed_multiplier": "1",
    "packet_seed_avalanche": "0",
    "loss": "packet-schedule-v1",
    "loss_trace": "common-id-v2",
}
# The cheap compare recovery gate and promotion-grade peeltiming path both
# bind the requested packet schedule.  ``packet-schedule-v1`` covers both the
# packet-ID ordering and its schedule-specific loss process (notably burst
# loss); the separate ``schedule`` field selects the exact algorithm.
TARGET_SCHEDULES = {
    "iid", "burst", "permutation", "systematic-first", "repair-only",
    "adversarial",
}
COMPARE_TARGET_SCHEDULES = set(TARGET_SCHEDULES)
PAIRED_ORDER = "ABBABAAB"
PAIRED_EXPERIMENTS = ("candidate-control", "identity-aa")
PAIRED_THERMAL_INTERVAL_MAX_MS = 2000
PAIRED_CACHE_STATE = "warm"
PAIRED_TIMING_SCOPE = "decoder-solve"
PAIRED_UNCERTAINTY_PROTOCOL = "paired-log-ratio-t95/v1"
PAIRED_CONFIDENCE = 0.95
PAIRED_LOSS_MODEL = "packet-schedule-v1"
PAIRED_CONSTRUCTION_SEED_DERIVATION = (
    "base_xor_d1b54a32d192ed03_times_rep_plus_1_v1"
)
PAIRED_LOSS_SEED_DERIVATION = (
    "base_xor_9e3779b97f4a7c15_times_rep_plus_1_v1"
)
PAIRED_CONTEXT_SCHEMA = "wirehair-v2-peeltiming-context/v1"
PAIRED_THERMAL_SCHEMA = "wirehair.wh2.thermal.v1"
PAIRED_CLOCK_DOMAIN = "linux-clock-monotonic"
PEELTIMING_EVICT_BYTES_MAX = 1 << 30
PEELTIMING_WORKING_SET_BYTES_MAX = 4 << 30
STOCK_CONTROL_SELECTED = "selected-paired-measurement"
STOCK_CONTROL_EMBEDDED = "embedded-paired-measurement"
STOCK_PMF_DIGEST = hashlib.sha256(b"stock").hexdigest()
PROXY_COST_MODEL = "embedded-es-cost-model/raw-calibration-unavailable"
PROXY_MEASURE_REGIME = {
    "solve_block_bytes": 0,
    "cost_model_block_bytes": 1280,
    "cost_model_verified": 0,
    "band_tracking_x": 1,
    "loss": "0.100000",
    "seed_base": 55,
    "completion": "mixed",
    "geometry": "frozen",
    "period": 244,
    "gf16_rows": 0,
}
PROXY_ORDERING_PROTOCOL = "fail-rate-then-pred-ns/v1"
SEARCH_BOX_PROTOCOL = "native-mean:[0,min(K,max(40,4*mean))]/v1"
PROXY_K_LADDER = (
    2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384,
    512, 768, 1024, 1536, 2048, 4096, 8192, 16384, 32768, 64000,
)
# These are the complete production mixed-completion settings reported by the
# native benchmark.  Search tools do not pass test-hook overrides for them:
# if a production default changes, the receipt must get a deliberate schema
# update rather than silently comparing a different codec under the old name.
PRODUCTION_COMPLETION_REGIME = {
    "precode": "1",
    "precode_cache": "0",
    "precode_profile": "mixed",
    "encoder_cache": "0",
    "decoder_cache": "0",
    "mixed_mix_count": "3",
    "mixed_period": "244",
    "mixed_gf256_rows": "10",
    "mixed_gf16_rows": "2",
    "mixed_geometry": "frozen",
    "mixed_residue_skew": "0",
    "mixed_residue_schedule": "constant",
    "mixed_residue_hash_seed": "0x0",
    "mixed_residue_hash_keyed": "0",
    "mixed_independent_extension_residues": "0",
    "mixed_grouped_gf256_rows": "0",
    "mixed_grouped_gf256_row_mask": "0x0",
    "mixed_grouped_gf256_hash_seed": "0x0",
    "mixed_grouped_final_h_a_columns": "0",
    "mixed_residue_buckets_requested": "auto",
    "packet_row_seed_multiplier": "0x1",
    "packet_row_seed_avalanche": "0",
    "precode_profile_handoff": "encoder-selected-v1",
}
# ValidatePrecodeParams() bounds K + S + D2 + H to a uint16 span.  The exact
# dispatch-v1 target fixes D2 at four rows and its completion height to the
# GF(256) plus GF(2^16) rows named in the receipt above.  Native PMF metadata
# comes from that contract, so a larger S could not have been produced by a
# valid target expansion.
_PRODUCTION_DENSE_ROWS = int(TARGET_CONTRACT["dense_rows"])
_PRODUCTION_TOTAL_SPAN_MAX = 0xffff
_PRODUCTION_MIXED_HEAVY_ROWS = (
    int(PRODUCTION_COMPLETION_REGIME["mixed_gf256_rows"]) +
    int(PRODUCTION_COMPLETION_REGIME["mixed_gf16_rows"])
)
# Exact WirehairTools::GetDenseCount() data used by dispatch-v1's
# SmallBandStaircaseCount().  Keep the legacy values here instead of trusting
# a self-reported native staircase: table documents are security boundaries,
# and a forged receipt must not be able to redefine the equation geometry.
_LEGACY_TINY_DENSE_COUNTS = (
    0, 0, 2, 3, 3, 5, 6, 6, 6, 7, 9, 10, 10, 10, 12, 14,
    13, 14, 12, 12, 15, 16, 21, 14, 14, 13, 18, 21, 22, 21, 13, 22,
    13, 24, 14, 17, 16, 24, 30, 26, 24, 18, 15, 15, 24, 18, 21, 17,
    14, 16, 21, 18, 17, 22, 25, 20, 17, 18, 21, 18, 23, 20, 19, 23,
    19,
)
_LEGACY_DENSE_POINTS = (
    (2048, 52),
    (2618, 54),
    (2826, 60),
    (3725, 62),
    (3962, 67),
    (4277, 65),
    (4547, 60),
    (5065, 64),
    (5224, 76),
    (5642, 66),
    (5909, 71),
    (6285, 76),
    (6583, 66),
    (6895, 72),
    (7448, 69),
    (7682, 76),
    (8046, 78),
    (8558, 76),
    (8963, 73),
    (9389, 81),
    (10143, 86),
    (11129, 94),
    (12593, 99),
    (12988, 105),
    (14032, 108),
    (14473, 114),
    (15397, 110),
    (16636, 113),
    (17698, 118),
    (18828, 123),
    (19420, 127),
    (20343, 136),
    (21979, 139),
    (23024, 150),
    (24119, 156),
    (25659, 162),
    (27298, 173),
    (29042, 176),
    (30898, 183),
    (31870, 190),
    (33906, 200),
    (35519, 211),
    (37208, 220),
    (38978, 234),
    (40205, 253),
    (42776, 297),
    (44122, 320),
    (45511, 336),
    (46944, 357),
    (48421, 373),
    (49177, 376),
    (50725, 380),
    (52321, 391),
    (53968, 388),
    (54811, 382),
    # This duplicate is intentional and mirrors the native graph exactly.
    (54811, 382),
    (55667, 372),
    (57419, 362),
    (58316, 356),
    (60152, 347),
    (61091, 337),
    (62045, 334),
    (63014, 340),
    (64000, 345),
)
# The native compare parser stores --trials in uint32_t and parses --bb-list
# through a positive int.  The production mixed arm additionally requires an
# even payload width for its GF(2^16) rows.
_COMPARE_TRIALS_MAX = 0xffffffff
_COMPARE_BLOCK_BYTES_MAX = 0x7fffffff
PRODUCTION_COMPARE_ARM = {
    "target_profile": TARGET_PROFILE,
    "seed_policy": TARGET_SEED_POLICY,
}
_COMPARE_BANNER_ARM = {
    "v2_profile": "base",
    "peel_candidates": "16",
    "peel_trials": "3",
    "auto_trials": "8",
    "auto_min_delta": "0.1000",
    "tune_seed": "0x9a7e11a",
    "auto_seed": "0xa570ca1",
    "dense_override": "0",
    "dense_delta": "0",
    "dense_candidate": "0",
}
RECOVERY_METRIC_SCOPE = {
    "fail": "all-trials",
    "overhead": "successful-trials-only",
    "throughput": "successful-trials-only",
}
RECOVERY_METRIC_FIELDS = {
    "construction_seed", "loss_seed", "target_receipt", "fail", "oh_mean",
    "OH_sd", "OH50", "OH95", "OH99", "OH_max", "decode_mbps",
}
PAIRED_ARM_FIELDS = {
    "trials", "fail", "oh_mean", "OH_sd", "OH50", "OH95", "OH99",
    "OH_max", "solve_mbps",
}
COMPARE_COLUMNS = (
    "codec", "bb", "trials", "fail", "N_mean", "OH_mean", "OH_sd",
    "OH50", "OH95", "OH99", "OH_max", "create_MBps", "encode_MBps",
    "decode_MBps", "recover_MBps",
)
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_SOURCE_EXTENSIONS = (
    ".c", ".cc", ".cmake", ".cpp", ".h", ".hpp", ".in", ".inc", ".map",
    ".py",
)


class MeasurementError(RuntimeError):
    """The benchmark failed or emitted an unusable receipt."""


class _ThermalCoveragePending(MeasurementError):
    """The sampler has not yet appended a complete post-run bracket."""


def valid_loss_rate(value):
    """Return whether value has the native benchmark's canonical loss form."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return False
    try:
        numeric = float(value)
    except (OverflowError, TypeError, ValueError):
        return False
    return (
        math.isfinite(numeric) and
        0.0 <= numeric <= 0.99 and
        not (numeric == 0.0 and math.copysign(1.0, numeric) < 0.0)
    )


def _is_canonical_nonnegative_float(value):
    """Accept one finite JSON float in [0,+inf), excluding signed zero."""
    return (
        type(value) is float and math.isfinite(value) and value >= 0.0 and
        not (value == 0.0 and math.copysign(1.0, value) < 0.0)
    )


def _is_canonical_scale(value):
    """Accept the exact unset sentinel or a canonical nonnegative scale."""
    return (
        type(value) is float and
        (value == -1.0 or
         (_is_canonical_nonnegative_float(value) and value <= 64000.0))
    )


def _validate_paired_degree_scale(value):
    """Reject aliases the native numeric scale parser would otherwise coerce."""
    if (value is not None and
            (not _is_canonical_nonnegative_float(value) or value > 64000.0)):
        raise ValueError(
            "paired staircase degree scale must be a canonical float")


def _is_canonical_absolute_path(value):
    """Return whether a receipt path is nonempty, absolute, and normalized."""
    if (not isinstance(value, str) or not value.startswith("/") or
            value.startswith("//") or "\0" in value or
            not os.path.isabs(value) or os.path.normpath(value) != value):
        return False
    try:
        os.fsencode(value)
    except UnicodeEncodeError:
        return False
    return True


def _canonical_staircase_scale_spec(value):
    """Return the native receipt spelling for a validated numeric scale."""
    numeric = float(value)
    if numeric == 0.0:
        numeric = 0.0
    return f"{numeric:.17g}"


def _nonnegative_llround(value):
    """Mirror C++ llround for a finite nonnegative binary64 value."""
    whole = math.floor(value)
    return whole + int(value - whole >= 0.5)


def _target_mean_spec(
        block_count, staircase, source_hits, staircase_scale):
    """Return native's exact printed mean staircase-row degree."""
    if staircase_scale == "unset":
        edge_count = block_count * min(source_hits, staircase)
    else:
        scale = float(staircase_scale)
        # MakeStaircaseDegreeMixture() uses llround(scale * S).  Mirror its
        # nonnegative halfway-away-from-zero rule without adding 0.5 first,
        # which can round upward prematurely at a binary64 boundary.
        edge_count = _nonnegative_llround(scale * staircase)
        edge_count = min(edge_count, block_count * staircase)
    return f"{edge_count / staircase:.17g}"


def _scale_realizes_legacy_edges(
        degree_scale, block_count, staircase, source_hits):
    """Whether a numeric scale constructs the exact legacy staircase budget."""
    if degree_scale is None:
        return True
    scaled_edges = min(
        _nonnegative_llround(float(degree_scale) * staircase),
        block_count * staircase)
    legacy_edges = block_count * min(source_hits, staircase)
    return scaled_edges == legacy_edges


@dataclass(frozen=True)
class StockProfile:
    """Exact dispatch-v1 peel distribution and construction contract."""

    block_count: int
    target_profile: str
    contract_id: str
    contract_sha256: str
    precode_contract: int
    packet_contract: int
    architecture: str
    staircase: int
    dense_rows: int
    heavy_rows: int
    source_hits: int
    completion: str
    gf256_rows: int
    gf16_rows: int
    period: int
    geometry: str
    residue_schedule: str
    residue_skew: int
    mix_count: int
    target_mean: float
    native_pmf_sha256: str
    pmf_encoding: str
    pmf: tuple

    def as_dict(self):
        return {
            "block_count": self.block_count,
            "target_profile": self.target_profile,
            "contract_id": self.contract_id,
            "contract_sha256": self.contract_sha256,
            "precode_contract": self.precode_contract,
            "packet_contract": self.packet_contract,
            "architecture": self.architecture,
            "staircase": self.staircase,
            "dense_rows": self.dense_rows,
            "heavy_rows": self.heavy_rows,
            "source_hits": self.source_hits,
            "completion": self.completion,
            "gf256_rows": self.gf256_rows,
            "gf16_rows": self.gf16_rows,
            "period": self.period,
            "geometry": self.geometry,
            "residue_schedule": self.residue_schedule,
            "residue_skew": self.residue_skew,
            "mix_count": self.mix_count,
            "target_mean": self.target_mean,
            "native_pmf_sha256": self.native_pmf_sha256,
            "pmf_encoding": self.pmf_encoding,
            "pmf_sha256": pmf_sha256(self.pmf),
            "pmf": list(self.pmf),
        }


@dataclass(frozen=True)
class RecoveryMetrics:
    """One exact raw-seed dispatch-v1 recovery measurement."""

    construction_seed: int
    loss_seed: int
    target_receipt: dict
    fail: int
    oh_mean: float
    oh_sd: float
    oh50: float
    oh95: float
    oh99: float
    oh_max: float
    decode_mbps: float

    def goodput(self, block_count):
        if self.fail != 0:
            return 0.0
        return self.decode_mbps * block_count / (block_count + self.oh_mean)

    def as_dict(self):
        return {
            "construction_seed": self.construction_seed,
            "loss_seed": self.loss_seed,
            "target_receipt": dict(self.target_receipt),
            "fail": self.fail,
            "oh_mean": self.oh_mean,
            "OH_sd": self.oh_sd,
            "OH50": self.oh50,
            "OH95": self.oh95,
            "OH99": self.oh99,
            "OH_max": self.oh_max,
            "decode_mbps": self.decode_mbps,
        }


PEELTIMING_MANIFEST_FIELDS = (
    "schema", "target_profile", "seed_policy", "contract_id",
    "K", "bb", "S", "D2", "H",
    "construction_seed_base", "construction_seed_derivation",
    "semantic_seed_derivation",
    "loss", "loss_seed_base", "loss_seed_derivation", "message_seed_policy",
    "schedule", "loss_model", "trace_encoding", "panels",
    "candidate_pmf_sha256", "candidate_pmf_encoding",
    "candidate_scale_requested",
    "candidate_scale_effective", "identity_pmf_sha256",
    "identity_pmf_encoding",
    "identity_scale_effective", "warmup_replicates", "replicates",
    "slots_per_panel", "panels_per_replicate", "order", "label_swap",
    "inner_reps", "max_overhead", "cache_state", "evict_bytes",
    "payload_alignment", "prefault", "cpu_affinity_policy", "payload",
    "timing_scope",
    "timing_prefix", "recovery_scope", "weak_seed_policy", "hook_path",
    "no_override_scope", "system_build", "startup_amortization", "slot_prewarm",
    "context_sha256", "uncertainty", "required_margin", "margin_rule",
    "clock_domain", "stream_hash_scope",
    "started_monotonic_ns", "expected_rows",
)
PEELTIMING_SEMANTIC_FIELDS = (
    "timed", "construction_seed", "seed_attempt", "seed_attempt_cap",
    "loss_seed", "trace_sha256", "message_sha256",
    "identity_pmf_sha256", "identity_pmf_encoding",
    "identity_scale_effective", "nohook_result", "nohook_recovery_ok",
    "nohook_overhead", "identity_result", "identity_recovery_ok",
    "identity_overhead", "nohook_system_sha256", "identity_system_sha256",
    "system_equal", "nohook_packet_rows_sha256",
    "identity_packet_rows_sha256", "packet_rows_equal",
    "nohook_payload_sha256", "identity_payload_sha256", "payload_equal",
    "nohook_direct_result", "identity_direct_result",
    "nohook_intermediate_bytes", "identity_intermediate_bytes",
    "nohook_solve_sha256", "identity_solve_sha256", "solve_equal",
    "full_recovery_equal", "pass",
)
PEELTIMING_COLUMNS = (
    "replicate", "measured", "panel", "slot", "pair", "label", "role",
    "label_swap", "construction_seed", "loss_seed", "trace_sha256",
    "common_overhead",
    "arm_overhead", "recovery_result", "recovery_ok", "timing_eligible",
    "preflight_result", "timing_result", "outcome_stable", "elapsed_ns",
    "inner_reps",
    "saturated", "cpu_before", "cpu_after", "cpu_migrated", "minflt_delta",
    "majflt_delta", "fault_contaminated", "inactivated", "binary_def",
    "heavy_gain", "block_xors", "block_muladds", "build_ns_sum",
    "peel_ns_sum", "project_ns_sum", "residual_ns_sum", "backsub_ns_sum",
    "source_bytes", "packet_payload_bytes", "intermediate_bytes",
)
_PEELTIMING_SCHEMA = "wirehair.wh2.peeltiming.v2"
_PEELTIMING_PANELS = "candidate_control+identity_aa"
_PEELTIMING_LABEL_SWAP = "alternating"
_PEELTIMING_PMF_ENCODING = "binary64-17g-comma-v1"
_PEELTIMING_IDENTITY_PMF_ENCODING = "stock-recovered-explicit-17g-v1"
_PEELTIMING_CLOCK_DOMAIN = "posix-clock-monotonic-ns-v1"
_PAIRED_BOUND_BASE_FIELDS = {
    "cpu_model", "cpu_affinity", "cpu_governors", "cache_state",
    "evict_bytes", "thermal_source", "thermal_sampling_interval_ms",
    "clock_domain",
}
_PAIRED_BOUND_CAPTURE_FIELDS = {
    "thermal_device", "thermal_inode", "thermal_prelaunch_monotonic_s",
    "thermal_prelaunch_row_sha256",
}
_PAIRED_UNSIGNED_NUMBER = re.compile(
    r"(?:0|[1-9][0-9]*)(?:[.][0-9]+)?\Z")
_PAIRED_UNSIGNED_COUNTER = re.compile(r"(?:0|[1-9][0-9]*)\Z")


@dataclass(frozen=True)
class PairedArmSummary:
    """Aggregate one logical peeltiming arm without discarding failures."""

    trials: int
    fail: int
    oh_mean: float
    oh_sd: float
    oh50: float
    oh95: float
    oh99: float
    oh_max: float
    solve_mbps: float

    def goodput(self, block_count):
        return self.solve_mbps * block_count / (block_count + self.oh_mean)

    def as_dict(self):
        return {
            "trials": self.trials,
            "fail": self.fail,
            "oh_mean": self.oh_mean,
            "OH_sd": self.oh_sd,
            "OH50": self.oh50,
            "OH95": self.oh95,
            "OH99": self.oh99,
            "OH_max": self.oh_max,
            "solve_mbps": self.solve_mbps,
        }


@dataclass(frozen=True)
class PairedMeasurement:
    """Strict aggregate of one complete two-panel native timing stream."""

    manifest: dict
    semantic: dict
    context: dict
    candidate: PairedArmSummary
    identity: PairedArmSummary
    replicate_receipts: tuple
    candidate_log_cost_mean: float
    candidate_log_cost_ci_low: float
    candidate_log_cost_ci_high: float
    aa_log_cost_mean: float
    aa_log_cost_ci_low: float
    aa_log_cost_ci_high: float
    aa_floor_log: float
    required_margin: float
    effective_margin_log: float
    timing_eligible_replicates: int
    timing_ci_available: bool
    candidate_recovery_result_counts: dict
    identity_recovery_result_counts: dict
    candidate_weak_seed_failures: int
    identity_weak_seed_failures: int
    recovery_regressions: int
    recovery_improvements: int
    both_failures: int
    overhead_regressions: int
    valid_for_promotion: bool
    stream_sha256: str
    finished_monotonic_ns: int
    final_context_sha256: str
    native_stdout: str

    def as_dict(self):
        return {
            "protocol": NATIVE_PAIRED_PROTOCOL,
            "manifest": dict(self.manifest),
            "semantic": dict(self.semantic),
            "context": {
                "schema": self.context["schema"],
                "bound": dict(self.context["bound"]),
                "thermal": dict(self.context["thermal"]),
            },
            "candidate": self.candidate.as_dict(),
            "identity": self.identity.as_dict(),
            "replicates": [dict(item) for item in self.replicate_receipts],
            "candidate_log_cost_mean": self.candidate_log_cost_mean,
            "candidate_log_cost_ci_low": self.candidate_log_cost_ci_low,
            "candidate_log_cost_ci_high": self.candidate_log_cost_ci_high,
            "aa_log_cost_mean": self.aa_log_cost_mean,
            "aa_log_cost_ci_low": self.aa_log_cost_ci_low,
            "aa_log_cost_ci_high": self.aa_log_cost_ci_high,
            "aa_floor_log": self.aa_floor_log,
            "required_margin": self.required_margin,
            "effective_margin_log": self.effective_margin_log,
            "timing_eligible_replicates": self.timing_eligible_replicates,
            "timing_ci_available": self.timing_ci_available,
            "candidate_recovery_result_counts":
                dict(self.candidate_recovery_result_counts),
            "identity_recovery_result_counts":
                dict(self.identity_recovery_result_counts),
            "candidate_weak_seed_failures":
                self.candidate_weak_seed_failures,
            "identity_weak_seed_failures":
                self.identity_weak_seed_failures,
            "recovery_regressions": self.recovery_regressions,
            "recovery_improvements": self.recovery_improvements,
            "both_failures": self.both_failures,
            "overhead_regressions": self.overhead_regressions,
            "valid_for_promotion": self.valid_for_promotion,
            "stream_sha256": self.stream_sha256,
            "evidence": {
                "context_sha256": self.manifest["context_sha256"],
                "stream_sha256": self.stream_sha256,
                "started_monotonic_ns":
                    self.manifest["started_monotonic_ns"],
                "finished_monotonic_ns": self.finished_monotonic_ns,
                "final_context_sha256": self.final_context_sha256,
            },
            "native_stdout": self.native_stdout,
        }


def _parse_ordered_kv_line(line, prefix, fields, label):
    """Parse one comma receipt while rejecting aliases and field reordering."""
    if not line.startswith(prefix):
        raise MeasurementError(f"missing {label} receipt")
    pieces = line[len(prefix):].split(",")
    if len(pieces) != len(fields):
        raise MeasurementError(f"{label} receipt has the wrong field count")
    result = {}
    names = []
    for piece in pieces:
        if piece.count("=") != 1:
            raise MeasurementError(f"malformed {label} field {piece!r}")
        name, value = piece.split("=", 1)
        if not name or not value or name in result:
            raise MeasurementError(f"duplicate or empty {label} field {name!r}")
        names.append(name)
        result[name] = value
    if tuple(names) != tuple(fields):
        raise MeasurementError(f"{label} fields are missing or reordered")
    return result


def _parse_decimal_integer(value, label, *, minimum=None, maximum=None):
    if not isinstance(value, str) or not value:
        raise MeasurementError(f"{label} is not an integer")
    try:
        parsed = int(value, 10)
    except ValueError:
        raise MeasurementError(f"{label} is not an integer")
    if str(parsed) != value:
        raise MeasurementError(f"{label} is not a canonical integer")
    if ((minimum is not None and parsed < minimum) or
            (maximum is not None and parsed > maximum)):
        raise MeasurementError(f"{label} is outside its integer domain")
    return parsed


def _parse_seed_integer(value, label):
    if not isinstance(value, str) or not value:
        raise MeasurementError(f"{label} is not a uint64")
    try:
        parsed = int(value, 0)
    except ValueError:
        raise MeasurementError(f"{label} is not a uint64")
    if not 0 <= parsed <= 0xffffffffffffffff:
        raise MeasurementError(f"{label} is outside the uint64 domain")
    canonical = str(parsed)
    canonical_hex = f"0x{parsed:x}"
    if value not in (canonical, canonical_hex):
        raise MeasurementError(f"{label} is not a canonical uint64")
    return parsed


def _parse_finite_float(value, label, *, minimum=None, maximum=None):
    try:
        parsed = float(value)
    except (OverflowError, TypeError, ValueError):
        raise MeasurementError(f"{label} is not numeric")
    if (not math.isfinite(parsed) or
            (minimum is not None and parsed < minimum) or
            (maximum is not None and parsed > maximum)):
        raise MeasurementError(f"{label} is outside its numeric domain")
    if parsed == 0.0:
        parsed = 0.0
    return parsed


def _canonical_json_sha256(value):
    """Hash one strict JSON value in the format bound by peeltiming."""
    try:
        encoded = json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode("utf-8")
    except (OverflowError, TypeError, ValueError) as error:
        raise MeasurementError(f"could not canonicalize timing context: {error}")
    return hashlib.sha256(encoded).hexdigest()


def _same_typed_json(left, right):
    """Compare JSON-shaped values without Python bool/int/float aliases."""
    if type(left) is not type(right):
        return False
    if isinstance(left, dict):
        return (
            set(left) == set(right) and
            all(_same_typed_json(left[key], right[key]) for key in left)
        )
    if isinstance(left, list):
        return (
            len(left) == len(right) and
            all(_same_typed_json(a, b) for a, b in zip(left, right))
        )
    if isinstance(left, float) and left == 0.0 and right == 0.0:
        return math.copysign(1.0, left) == math.copysign(1.0, right)
    return left == right


def _canonical_paired_pmf_spec(pmf):
    """Mirror peeltiming's parsed binary64 PMF canonicalization exactly."""
    try:
        numeric = []
        for value in pmf:
            if isinstance(value, bool):
                raise ValueError
            parsed = float(value)
            if parsed == 0.0:
                parsed = 0.0
            numeric.append(parsed)
    except (OverflowError, TypeError, ValueError):
        raise ValueError("invalid paired PMF")
    total = sum(numeric)
    if (not 1 <= len(numeric) <= 64 or
            any(value == 0.0 and math.copysign(1.0, value) < 0.0
                for value in numeric) or
            any(not math.isfinite(value) or value < 0.0 for value in numeric) or
            not math.isfinite(total) or not total > 0.0):
        raise ValueError("invalid paired PMF")
    return ",".join(f"{value:.17g}" for value in numeric)


def _paired_pmf_sha256(pmf):
    return hashlib.sha256(
        _canonical_paired_pmf_spec(pmf).encode("ascii")).hexdigest()


def _native_peel_cdf_signature(pmf, block_count):
    """Return exact uint32 sampler cutpoints after the native K/2 clamp."""
    if type(block_count) is not int or not 2 <= block_count <= 64000:
        raise ValueError("invalid paired block count")
    try:
        numeric = []
        for value in pmf:
            if isinstance(value, bool):
                raise ValueError
            parsed = float(value)
            if parsed == 0.0:
                parsed = 0.0
            numeric.append(parsed)
    except (OverflowError, TypeError, ValueError):
        raise ValueError("invalid paired PMF")
    if (not 1 <= len(numeric) <= 64 or
            any(not math.isfinite(value) or value < 0.0
                for value in numeric)):
        raise ValueError("invalid paired PMF")
    total = 0.0
    for value in numeric:
        total += value
        if not math.isfinite(total):
            raise ValueError("invalid paired PMF")
    if not total > 0.0:
        raise ValueError("invalid paired PMF")
    cdf = [1.0] * 64
    running = 0.0
    for index, value in enumerate(numeric[:63]):
        running += value
        probability = running / total
        if (not math.isfinite(probability) or
                not 0.0 <= probability <= 1.0):
            raise ValueError("invalid paired PMF")
        cdf[index] = probability
    cdf[-1] = 1.0
    realized_max_degree = min(64, block_count // 2)
    uint32_scale = 4294967296.0
    return tuple(
        max(0, min(1 << 32, math.ceil(cdf[index] * uint32_scale)))
        for index in range(realized_max_degree - 1)
    )


def _validate_paired_bound_context(
        context, *, cache_state, evict_bytes, require_capture=True):
    """Validate and hash the context subset available before native startup."""
    if (not isinstance(context, dict) or
            set(context) != {"schema", "bound", "thermal"} or
            context.get("schema") != PAIRED_CONTEXT_SCHEMA or
            not isinstance(context.get("bound"), dict) or
            (context.get("thermal") is not None and
             not isinstance(context.get("thermal"), dict))):
        raise MeasurementError("paired timing context is incomplete")
    bound = context["bound"]
    required_bound = set(_PAIRED_BOUND_BASE_FIELDS)
    if require_capture:
        required_bound.update(_PAIRED_BOUND_CAPTURE_FIELDS)
    affinity = bound.get("cpu_affinity")
    governors = bound.get("cpu_governors")
    if (set(bound) != required_bound or
            not isinstance(bound.get("cpu_model"), str) or
            not bound["cpu_model"] or
            not isinstance(affinity, list) or not affinity or
            any(type(cpu) is not int or not 0 <= cpu <= 0x7fffffff
                for cpu in affinity) or
            affinity != sorted(set(affinity)) or
            not isinstance(governors, dict) or
            set(governors) != {str(cpu) for cpu in affinity} or
            any(not isinstance(value, str) or not value or
                any(character.isspace() for character in value)
                for value in governors.values()) or
            bound.get("cache_state") != cache_state or
            type(bound.get("evict_bytes")) is not int or
            bound["evict_bytes"] != evict_bytes or
            not 4096 <= bound["evict_bytes"] <=
                PEELTIMING_EVICT_BYTES_MAX or
            not _is_canonical_absolute_path(
                bound.get("thermal_source")) or
            type(bound.get("thermal_sampling_interval_ms")) is not int or
            not 1 <= bound["thermal_sampling_interval_ms"] <=
                PAIRED_THERMAL_INTERVAL_MAX_MS or
            bound.get("clock_domain") != PAIRED_CLOCK_DOMAIN):
        raise MeasurementError("paired timing bound context is invalid")
    if require_capture and (
            type(bound.get("thermal_device")) is not int or
            not 0 <= bound["thermal_device"] <= 0xffffffffffffffff or
            type(bound.get("thermal_inode")) is not int or
            not 0 <= bound["thermal_inode"] <= 0xffffffffffffffff or
            type(bound.get("thermal_prelaunch_monotonic_s")) is not float or
            not math.isfinite(bound["thermal_prelaunch_monotonic_s"]) or
            bound["thermal_prelaunch_monotonic_s"] < 0.0 or
            (bound["thermal_prelaunch_monotonic_s"] == 0.0 and
             math.copysign(
                 1.0, bound["thermal_prelaunch_monotonic_s"]) < 0.0) or
            not _is_sha256(bound.get("thermal_prelaunch_row_sha256"))):
        raise MeasurementError("paired timing capture context is invalid")
    return _canonical_json_sha256(bound)


def _validate_paired_context(
        context, *, cache_state, evict_bytes, started_ns=None, finished_ns=None):
    """Validate CPU/cache identity plus thermal evidence covering the run."""
    bound_sha256 = _validate_paired_bound_context(
        context, cache_state=cache_state, evict_bytes=evict_bytes,
        require_capture=True)
    bound = context["bound"]
    thermal = context["thermal"]
    required_thermal = {
        "schema", "summary_sha256", "csv_sha256", "csv_window",
        "prelaunch_row",
        "rows", "valid_rows",
        "missing_read_rows", "invalid_rows", "first_monotonic_s",
        "last_monotonic_s", "cpu_tctl_max_c", "dimm_max_c",
        "edac_ce_max", "edac_ue_max",
    }
    integer_fields = {
        "rows", "valid_rows", "missing_read_rows", "invalid_rows",
        "edac_ce_max", "edac_ue_max",
    }
    float_fields = {
        "first_monotonic_s", "last_monotonic_s", "cpu_tctl_max_c",
        "dimm_max_c",
    }
    summary_payload = (
        {name: value for name, value in thermal.items()
         if name != "summary_sha256"}
        if isinstance(thermal, dict) else None
    )
    if (not isinstance(thermal, dict) or
            set(thermal) != required_thermal or
            thermal.get("schema") != PAIRED_THERMAL_SCHEMA or
            not _is_sha256(thermal.get("summary_sha256")) or
            not _is_sha256(thermal.get("csv_sha256")) or
            not isinstance(thermal.get("csv_window"), str) or
            not thermal["csv_window"] or
            not isinstance(thermal.get("prelaunch_row"), str) or
            not thermal["prelaunch_row"] or
            thermal["summary_sha256"] !=
                _canonical_json_sha256(summary_payload) or
            any(type(thermal.get(name)) is not int for name in integer_fields) or
            any(not _is_canonical_nonnegative_float(thermal.get(name))
                for name in float_fields) or
            thermal["rows"] < 2 or thermal["valid_rows"] < 2 or
            thermal["valid_rows"] > thermal["rows"] or
            thermal["missing_read_rows"] < 0 or
            thermal["missing_read_rows"] > thermal["rows"] or
            thermal["invalid_rows"] < 0 or
            thermal["valid_rows"] + thermal["invalid_rows"] !=
                thermal["rows"] or
            thermal["invalid_rows"] != thermal["missing_read_rows"] or
            bound["thermal_prelaunch_monotonic_s"] >
                thermal["first_monotonic_s"] or
            thermal["first_monotonic_s"] < 0.0 or
            thermal["last_monotonic_s"] <= thermal["first_monotonic_s"] or
            not 0.0 <= thermal["cpu_tctl_max_c"] <= 120.0 or
            not 0.0 <= thermal["dimm_max_c"] <= 100.0 or
            thermal["edac_ce_max"] != 0 or thermal["edac_ue_max"] != 0):
        raise MeasurementError("paired timing thermal context is invalid")
    try:
        csv_payload = thermal["csv_window"].encode("ascii")
    except UnicodeEncodeError:
        raise MeasurementError("paired timing thermal CSV is not ASCII")
    if hashlib.sha256(csv_payload).hexdigest() != thermal["csv_sha256"]:
        raise MeasurementError("paired timing thermal CSV digest is invalid")
    try:
        prelaunch_payload = thermal["prelaunch_row"].encode("ascii")
    except UnicodeEncodeError:
        raise MeasurementError("paired timing prelaunch row is not ASCII")
    if (not prelaunch_payload.endswith(b"\n") or
            prelaunch_payload.count(b"\n") != 1 or
            hashlib.sha256(prelaunch_payload).hexdigest() !=
                bound["thermal_prelaunch_row_sha256"]):
        raise MeasurementError(
            "paired timing prelaunch row digest is invalid")
    prelaunch = _parse_thermal_record(prelaunch_payload, -1)
    interval = bound["thermal_sampling_interval_ms"] / 1000.0
    freshness = _thermal_startup_freshness(interval)
    start_s = (
        started_ns / 1e9 if started_ns is not None
        else thermal["first_monotonic_s"])
    if (not prelaunch["valid"] or
            prelaunch["monotonic"] !=
                bound["thermal_prelaunch_monotonic_s"] or
            prelaunch["monotonic"] > thermal["first_monotonic_s"] or
            prelaunch["monotonic"] > start_s or
            start_s - prelaunch["monotonic"] > freshness):
        raise MeasurementError(
            "paired timing prelaunch row is invalid or stale")
    replay_start_ns = (
        started_ns if started_ns is not None else
        int(round(thermal["first_monotonic_s"] * 1e9))
    )
    replay_finish_ns = (
        finished_ns if finished_ns is not None else
        int(round(thermal["last_monotonic_s"] * 1e9))
    )
    replayed = _thermal_window_from_csv(
        csv_payload, replay_start_ns, replay_finish_ns,
        sampling_interval_ms=bound["thermal_sampling_interval_ms"],
        prelaunch_row=prelaunch_payload)
    if not _same_typed_json(replayed, thermal):
        raise MeasurementError(
            "paired timing thermal summary does not match embedded CSV")
    if (started_ns is not None and
            thermal["first_monotonic_s"] > started_ns / 1e9):
        raise MeasurementError("thermal evidence starts after peeltiming")
    if (finished_ns is not None and
            thermal["last_monotonic_s"] < finished_ns / 1e9):
        raise MeasurementError("thermal evidence ends before peeltiming")
    return bound_sha256


def read_paired_context(path, *, require_thermal=True):
    """Read one stable, strict external CPU/cache/thermal context artifact."""
    try:
        with open(path, "rb") as source:
            before = os.fstat(source.fileno())
            if not stat.S_ISREG(before.st_mode):
                raise MeasurementError(
                    f"paired timing context is not a regular file: {path!r}")
            payload = source.read()
            after = os.fstat(source.fileno())
    except OSError as error:
        raise MeasurementError(
            f"could not read paired timing context {path!r}: {error}")
    stable_fields = (
        "st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns",
        "st_ctime_ns")
    if (any(getattr(before, name) != getattr(after, name)
            for name in stable_fields) or len(payload) != after.st_size):
        raise MeasurementError(
            f"paired timing context changed while reading: {path!r}")
    try:
        value = strict_json_loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise MeasurementError(f"invalid paired timing context JSON: {error}")
    # Run bounds are checked after the native receipt is parsed.
    if require_thermal:
        _validate_paired_context(
            value,
            cache_state=value.get("bound", {}).get("cache_state"),
            evict_bytes=value.get("bound", {}).get("evict_bytes"))
    else:
        bound = value.get("bound", {}) if isinstance(value, dict) else {}
        require_capture = (
            isinstance(bound, dict) and
            set(bound) == _PAIRED_BOUND_BASE_FIELDS |
                _PAIRED_BOUND_CAPTURE_FIELDS
        )
        _validate_paired_bound_context(
            value, cache_state=bound.get("cache_state"),
            evict_bytes=bound.get("evict_bytes"),
            require_capture=require_capture)
    return value


def paired_context_thermal_source(path):
    """Return the canonical live sampler protected by an input context."""
    context = read_paired_context(path, require_thermal=False)
    return context["bound"]["thermal_source"]


_THERMAL_CSV_COLUMNS = (
    "utc", "monotonic_s", "cpu_busy_pct", "cpu_avg_mhz", "cpu_tctl_c",
    "dimm_i2c1_50_c", "dimm_i2c1_51_c", "dimm_i2c1_52_c",
    "dimm_i2c1_53_c", "dimm_i2c2_50_c", "dimm_i2c2_51_c",
    "dimm_i2c2_52_c", "dimm_i2c2_53_c", "dimm_read_errors", "load1",
    "load5", "load15", "edac_ce", "edac_ue",
)
_THERMAL_CAPTURE_TAIL_BYTES = 64 * 1024
_THERMAL_CAPTURE_BYTES_MAX = 16 * 1024 * 1024


def _peeltiming_run_bounds(stdout):
    """Read authenticated run bounds before final thermal parsing."""
    lines = stdout.splitlines(keepends=True)
    if len(lines) < 4 or any(not line.endswith("\n") for line in lines):
        raise MeasurementError("peeltiming output is incomplete")
    manifest = _parse_ordered_kv_line(
        lines[0][:-1], "# peeltiming,", PEELTIMING_MANIFEST_FIELDS,
        "peeltiming manifest")
    done = _parse_ordered_kv_line(
        lines[-1][:-1], "# peeltiming_done,",
        ("complete", "rows", "finished_monotonic_ns", "stream_sha256"),
        "peeltiming completion")
    done_prefix = lines[-1][:-1].rsplit("stream_sha256=", 1)[0]
    done_prefix += "stream_sha256="
    stream_sha256 = hashlib.sha256(
        ("".join(lines[:-1]) + done_prefix).encode("ascii")).hexdigest()
    if (done["complete"] != "1" or
            done["stream_sha256"] != stream_sha256):
        raise MeasurementError(
            "cannot finalize thermal evidence for an incomplete peeltiming run")
    started_ns = _parse_decimal_integer(
        manifest["started_monotonic_ns"],
        "peeltiming started_monotonic_ns", minimum=1,
        maximum=0xffffffffffffffff)
    finished_ns = _parse_decimal_integer(
        done["finished_monotonic_ns"],
        "peeltiming finished_monotonic_ns", minimum=1,
        maximum=0xffffffffffffffff)
    if finished_ns <= started_ns:
        raise MeasurementError("peeltiming completion clock did not advance")
    return started_ns, finished_ns


@dataclass(frozen=True)
class _PairedThermalCapture:
    fd: int
    path: str
    device: int
    inode: int
    initial_size: int
    tail_offset: int
    prefix: bytes
    prelaunch_line_index: int
    prelaunch_line: bytes
    prelaunch_monotonic_s: float


def _paired_number(text, label):
    if (not isinstance(text, str) or
            _PAIRED_UNSIGNED_NUMBER.fullmatch(text) is None):
        raise MeasurementError(f"{label} is not a canonical number")
    value = float(text)
    if not math.isfinite(value):
        raise MeasurementError(f"{label} is not finite")
    return value


def _paired_counter(text, label):
    if (not isinstance(text, str) or
            _PAIRED_UNSIGNED_COUNTER.fullmatch(text) is None):
        raise MeasurementError(f"{label} is not a canonical counter")
    try:
        value = int(text, 10)
    except ValueError:
        raise MeasurementError(f"{label} is outside the uint64 domain")
    if value > 0xffffffffffffffff:
        raise MeasurementError(f"{label} is outside the uint64 domain")
    return value


def _parse_thermal_record(raw, index):
    """Parse one exact sampler row without discarding tolerated read gaps."""
    if not isinstance(raw, bytes) or not raw.endswith(b"\n") or b"\r" in raw:
        return {
            "raw": raw, "monotonic": None, "valid": False,
            "tolerated_missing": False, "hard_invalid": True,
        }
    try:
        fields = raw[:-1].decode("ascii").split(",")
    except UnicodeDecodeError:
        fields = []
    if len(fields) != len(_THERMAL_CSV_COLUMNS):
        return {
            "raw": raw, "monotonic": None, "valid": False,
            "tolerated_missing": False, "hard_invalid": True,
        }
    try:
        monotonic = _paired_number(
            fields[1], f"thermal CSV row {index}.monotonic_s")
    except MeasurementError:
        return {
            "raw": raw, "monotonic": None, "valid": False,
            "tolerated_missing": False, "hard_invalid": True,
        }
    try:
        busy = _paired_number(
            fields[2], f"thermal CSV row {index}.cpu_busy_pct")
        mhz = _paired_number(
            fields[3], f"thermal CSV row {index}.cpu_avg_mhz")
        cpu = _paired_number(
            fields[4], f"thermal CSV row {index}.cpu_tctl_c")
        read_errors = _paired_counter(
            fields[13], f"thermal CSV row {index}.dimm_read_errors")
        loads = [
            _paired_number(
                fields[column],
                f"thermal CSV row {index}.{_THERMAL_CSV_COLUMNS[column]}")
            for column in (14, 15, 16)
        ]
        ce = _paired_counter(
            fields[17], f"thermal CSV row {index}.edac_ce")
        ue = _paired_counter(
            fields[18], f"thermal CSV row {index}.edac_ue")
        dimms = [
            None if fields[column] == "" else _paired_number(
                fields[column],
                f"thermal CSV row {index}.{_THERMAL_CSV_COLUMNS[column]}")
            for column in range(5, 13)
        ]
    except MeasurementError:
        return {
            "raw": raw, "monotonic": monotonic, "valid": False,
            "tolerated_missing": False, "hard_invalid": True,
        }
    missing = any(value is None for value in dimms)
    present_dimms = [value for value in dimms if value is not None]
    other_valid = (
        fields[0] != "" and monotonic >= 0.0 and
        0.0 <= busy <= 100.0 and mhz > 0.0 and
        0.0 <= cpu <= 120.0 and
        all(value >= 0.0 for value in loads) and
        all(0.0 <= value <= 100.0 for value in present_dimms)
    )
    tolerated_missing = (
        other_valid and read_errors > 0 and missing and ce == 0 and ue == 0
    )
    valid = (
        other_valid and read_errors == 0 and not missing and
        len(present_dimms) == 8 and ce == 0 and ue == 0
    )
    return {
        "raw": raw,
        "monotonic": monotonic,
        "cpu": cpu,
        "dimms": present_dimms,
        "read_errors": read_errors,
        "ce": ce,
        "ue": ue,
        "valid": valid,
        "tolerated_missing": tolerated_missing,
        "hard_invalid": not valid and not tolerated_missing,
    }


def _actual_cpu_model():
    try:
        with open("/proc/cpuinfo", encoding="ascii") as source:
            models = {
                line.split(":", 1)[1].strip()
                for line in source
                if line.startswith("model name") and ":" in line
            }
    except OSError as error:
        raise MeasurementError(f"could not read CPU model: {error}")
    if not models:
        fallback = platform.processor().strip()
        if not fallback:
            raise MeasurementError("could not identify CPU model")
        models = {fallback}
    return " | ".join(sorted(models))


def _thermal_max_gap(interval_s):
    return min(5.0, max(0.25, 2.5 * interval_s))


def _thermal_startup_freshness(interval_s):
    return min(5.0, max(1.0, 2.5 * interval_s))


def _actual_cpu_affinity():
    try:
        affinity = sorted(os.sched_getaffinity(0))
    except (AttributeError, OSError) as error:
        raise MeasurementError(f"could not read CPU affinity: {error}")
    if not affinity:
        raise MeasurementError("CPU affinity is empty")
    return affinity


def _actual_cpu_governors(affinity):
    governors = {}
    for cpu in affinity:
        path = (
            f"/sys/devices/system/cpu/cpu{cpu}/cpufreq/scaling_governor"
        )
        try:
            with open(path, encoding="ascii") as source:
                value = source.read().strip()
        except OSError as error:
            raise MeasurementError(
                f"could not read governor for CPU {cpu}: {error}")
        if not value or any(character.isspace() for character in value):
            raise MeasurementError(f"CPU {cpu} governor is malformed")
        governors[str(cpu)] = value
    return governors


def _validate_runtime_identity(bound):
    affinity = _actual_cpu_affinity()
    if (bound["cpu_model"] != _actual_cpu_model() or
            bound["cpu_affinity"] != affinity or
            bound["cpu_governors"] != _actual_cpu_governors(affinity)):
        raise MeasurementError(
            "CPU model, affinity, or governors changed during peeltiming")


def make_paired_context_config(
        thermal_source, thermal_sampling_interval_ms, *,
        cache_state, evict_bytes):
    """Capture the static host policy used to start one paired experiment."""
    affinity = _actual_cpu_affinity()
    context = {
        "schema": PAIRED_CONTEXT_SCHEMA,
        "bound": {
            "cpu_model": _actual_cpu_model(),
            "cpu_affinity": affinity,
            "cpu_governors": _actual_cpu_governors(affinity),
            "cache_state": cache_state,
            "evict_bytes": evict_bytes,
            "thermal_source": os.path.realpath(thermal_source),
            "thermal_sampling_interval_ms": thermal_sampling_interval_ms,
            "clock_domain": PAIRED_CLOCK_DOMAIN,
        },
        "thermal": None,
    }
    _validate_paired_bound_context(
        context, cache_state=cache_state, evict_bytes=evict_bytes,
        require_capture=False)
    return context


def _capture_payload_size(capture, payload):
    """Recover the absolute file size represented by one bounded snapshot."""
    header_size = len((",".join(_THERMAL_CSV_COLUMNS) + "\n").encode("ascii"))
    if capture.tail_offset == 0:
        return len(payload)
    if len(payload) < header_size:
        raise MeasurementError("thermal CSV header is incomplete")
    return capture.tail_offset + len(payload) - header_size


def _snapshot_capture(capture):
    """Read only the fixed recent tail plus bytes appended during this run."""
    header = (",".join(_THERMAL_CSV_COLUMNS) + "\n").encode("ascii")
    try:
        descriptor = os.fstat(capture.fd)
        path_stat = os.lstat(capture.path)
        if (descriptor.st_dev != capture.device or
                descriptor.st_ino != capture.inode or
                path_stat.st_dev != capture.device or
                path_stat.st_ino != capture.inode or
                not stat.S_ISREG(descriptor.st_mode) or
                not stat.S_ISREG(path_stat.st_mode)):
            raise MeasurementError(
                "thermal sampler path or descriptor identity changed")
        if descriptor.st_size < capture.initial_size:
            raise MeasurementError(
                "thermal sampler was truncated during peeltiming")
        if capture.tail_offset == 0:
            read_offset = 0
            read_size = descriptor.st_size
            if read_size > _THERMAL_CAPTURE_BYTES_MAX:
                raise MeasurementError(
                    "thermal sampler append exceeds the bounded capture")
            payload = os.pread(capture.fd, read_size, read_offset)
            expected_size = read_size
        else:
            if descriptor.st_size < capture.tail_offset:
                raise MeasurementError(
                    "thermal sampler was truncated during peeltiming")
            tail_size = descriptor.st_size - capture.tail_offset
            if tail_size > _THERMAL_CAPTURE_BYTES_MAX:
                raise MeasurementError(
                    "thermal sampler append exceeds the bounded capture")
            header_payload = os.pread(capture.fd, len(header), 0)
            tail_payload = os.pread(
                capture.fd, tail_size, capture.tail_offset)
            payload = header_payload + tail_payload
            expected_size = len(header) + tail_size
        after = os.fstat(capture.fd)
        path_after = os.lstat(capture.path)
    except OSError as error:
        raise MeasurementError(f"could not read retained thermal CSV: {error}")
    if (after.st_dev != capture.device or after.st_ino != capture.inode or
            path_after.st_dev != capture.device or
            path_after.st_ino != capture.inode or
            not stat.S_ISREG(path_after.st_mode) or
            after.st_size < descriptor.st_size or
            len(payload) != expected_size or
            not payload.startswith(capture.prefix)):
        raise MeasurementError(
            "thermal sampler path or descriptor identity changed or was "
            "modified during peeltiming")
    return payload


def _prepare_paired_context(initial_context, *, cache_state, evict_bytes):
    """Open and authenticate one fresh sampler before starting native code."""
    context = copy.deepcopy(initial_context)
    if (not isinstance(context, dict) or
            context.get("thermal") is not None or
            not isinstance(context.get("bound"), dict) or
            set(context["bound"]) != _PAIRED_BOUND_BASE_FIELDS):
        raise MeasurementError(
            "automated peeltiming requires an unfinalized context")
    _validate_paired_bound_context(
        context, cache_state=cache_state, evict_bytes=evict_bytes,
        require_capture=False)
    _validate_runtime_identity(context["bound"])
    supplied_path = context["bound"]["thermal_source"]
    real_path = os.path.realpath(supplied_path)
    fd = None
    try:
        flags = os.O_RDONLY | os.O_CLOEXEC
        if hasattr(os, "O_NOFOLLOW"):
            flags |= os.O_NOFOLLOW
        fd = os.open(real_path, flags)
        descriptor = os.fstat(fd)
        path_stat = os.lstat(real_path)
    except OSError as error:
        if fd is not None:
            os.close(fd)
        raise MeasurementError(
            f"could not open thermal sampler {supplied_path!r}: {error}")
    if (not stat.S_ISREG(descriptor.st_mode) or
            not stat.S_ISREG(path_stat.st_mode) or
            descriptor.st_dev != path_stat.st_dev or
            descriptor.st_ino != path_stat.st_ino):
        os.close(fd)
        raise MeasurementError("thermal sampler is not one canonical file")
    provisional = _PairedThermalCapture(
        fd=fd, path=real_path, device=descriptor.st_dev,
        inode=descriptor.st_ino, initial_size=descriptor.st_size,
        tail_offset=0, prefix=b"", prelaunch_line_index=-1,
        prelaunch_line=b"", prelaunch_monotonic_s=0.0)
    try:
        header = (",".join(_THERMAL_CSV_COLUMNS) + "\n").encode("ascii")
        if descriptor.st_size > len(header):
            desired = max(
                len(header),
                descriptor.st_size - _THERMAL_CAPTURE_TAIL_BYTES)
            if desired > len(header):
                alignment = os.pread(
                    fd, descriptor.st_size - desired, desired)
                newline = alignment.find(b"\n")
                if newline < 0:
                    raise MeasurementError(
                        "recent thermal sampler row exceeds the bounded "
                        "capture")
                desired += newline + 1
            provisional = replace(provisional, tail_offset=desired)
        interval = (
            context["bound"]["thermal_sampling_interval_ms"] / 1000.0)
        deadline = time.monotonic() + 2.0 + 2.0 * interval
        while True:
            payload = _snapshot_capture(provisional)
            lines = payload.splitlines(keepends=True)
            if not lines or lines[0] != header:
                if (not payload or len(payload) >= len(header) or
                        not header.startswith(payload)):
                    raise MeasurementError("thermal CSV header is invalid")
                pending_part = "header"
            else:
                if payload.endswith(b"\n"):
                    break
                # A sampler append is not atomic with respect to our
                # fstat/pread snapshot. Retry only that transient case,
                # through the retained descriptor. Any already-complete
                # malformed row still fails immediately.
                complete_lines = lines[:-1]
                complete_tail_start = max(1, len(complete_lines) - 8)
                complete_tail = [
                    _parse_thermal_record(complete_lines[index], index)
                    for index in range(
                        complete_tail_start, len(complete_lines))
                ]
                if any(row["hard_invalid"] for row in complete_tail):
                    raise MeasurementError(
                        "recent thermal sampler rows are malformed")
                pending_part = "tail"
            if time.monotonic() >= deadline:
                raise MeasurementError(
                    f"thermal CSV {pending_part} remained incomplete")
            time.sleep(0.1)
        tail_start = max(1, len(lines) - 8)
        parsed_tail = [
            _parse_thermal_record(lines[index], index)
            for index in range(tail_start, len(lines))
        ]
        if (not parsed_tail or
                any(row["hard_invalid"] for row in parsed_tail)):
            raise MeasurementError(
                "recent thermal sampler rows are malformed")
        valid = [
            (tail_start + offset, row)
            for offset, row in enumerate(parsed_tail) if row["valid"]
        ]
        if not valid:
            raise MeasurementError("thermal sampler has no recent valid row")
        selected_index, selected = valid[-1]
        newest = max(row["monotonic"] for row in parsed_tail)
        now = time.clock_gettime(time.CLOCK_MONOTONIC)
        freshness = _thermal_startup_freshness(interval)
        if (newest - selected["monotonic"] > 5.0 or
                now - selected["monotonic"] > freshness or
                selected["monotonic"] - now > 1.0):
            raise MeasurementError("thermal sampler prelaunch row is stale")
        if len(valid) >= 2:
            previous = valid[-2][1]["monotonic"]
            if (selected["monotonic"] <= previous or
                    selected["monotonic"] - previous >
                        _thermal_max_gap(interval)):
                raise MeasurementError(
                    "thermal sampler interval does not match its context")
        capture = _PairedThermalCapture(
            fd=fd, path=real_path, device=descriptor.st_dev,
            inode=descriptor.st_ino,
            initial_size=_capture_payload_size(provisional, payload),
            tail_offset=provisional.tail_offset, prefix=payload,
            prelaunch_line_index=selected_index,
            prelaunch_line=lines[selected_index],
            prelaunch_monotonic_s=float(selected["monotonic"]))
        bound = dict(context["bound"])
        bound["thermal_source"] = real_path
        bound["thermal_device"] = descriptor.st_dev
        bound["thermal_inode"] = descriptor.st_ino
        bound["thermal_prelaunch_monotonic_s"] = float(
            selected["monotonic"])
        bound["thermal_prelaunch_row_sha256"] = hashlib.sha256(
            lines[selected_index]).hexdigest()
        frozen = {
            "schema": context["schema"],
            "bound": bound,
            "thermal": None,
        }
        _validate_paired_bound_context(
            frozen, cache_state=cache_state, evict_bytes=evict_bytes,
            require_capture=True)
        return copy.deepcopy(frozen), capture
    except BaseException:
        os.close(fd)
        raise


def _thermal_window_from_csv(
        payload, start_ns, finish_ns, capture=None,
        sampling_interval_ms=1000, prelaunch_row=None):
    """Select and hash only the exact sampler window bracketing native time."""
    if not isinstance(payload, bytes):
        raise MeasurementError("thermal CSV payload is not bytes")
    if (type(sampling_interval_ms) is not int or
            not 1 <= sampling_interval_ms <=
                PAIRED_THERMAL_INTERVAL_MAX_MS):
        raise MeasurementError("thermal sampling interval is invalid")
    interval = sampling_interval_ms / 1000.0
    maximum_gap = _thermal_max_gap(interval)
    header = (",".join(_THERMAL_CSV_COLUMNS) + "\n").encode("ascii")
    lines = payload.splitlines(keepends=True)
    if not lines or lines[0] != header:
        raise MeasurementError("thermal CSV header is invalid")
    if not payload.endswith(b"\n"):
        raise _ThermalCoveragePending("thermal CSV tail is incomplete")
    first_search = 1
    if capture is not None:
        if (not payload.startswith(capture.prefix) or
                capture.prelaunch_line_index >= len(lines) or
                lines[capture.prelaunch_line_index] !=
                    capture.prelaunch_line):
            raise MeasurementError(
                "thermal CSV prelaunch evidence changed")
        first_search = capture.prelaunch_line_index
        prelaunch_row = capture.prelaunch_line
    if prelaunch_row is None:
        prelaunch_row = lines[first_search]
    if (not isinstance(prelaunch_row, bytes) or
            not prelaunch_row.endswith(b"\n")):
        raise MeasurementError("thermal prelaunch row is malformed")
    parsed = [
        _parse_thermal_record(lines[index], index)
        for index in range(first_search, len(lines))
    ]
    start_s = start_ns / 1e9
    finish_s = finish_ns / 1e9
    if not start_s < finish_s:
        raise MeasurementError("peeltiming run bounds are reversed")
    first = None
    last = None
    for offset, row in enumerate(parsed):
        if row["valid"] and row["monotonic"] <= start_s:
            first = offset
        if (first is not None and row["valid"] and
                row["monotonic"] >= finish_s):
            last = offset
            break
    if first is None:
        raise MeasurementError(
            "thermal CSV lacks its authenticated pre-run bracket")
    if last is None:
        raise _ThermalCoveragePending(
            "thermal CSV does not yet cover peeltiming")
    window = parsed[first:last + 1]
    previous = None
    for row in window:
        monotonic = row["monotonic"]
        if (row["hard_invalid"] or monotonic is None or
                (previous is not None and
                 (monotonic <= previous or
                  monotonic - previous > maximum_gap))):
            raise MeasurementError(
                "selected thermal window is malformed, reordered, or gapped")
        previous = monotonic
    valid_rows = [row for row in window if row["valid"]]
    if (len(valid_rows) < 2 or
            valid_rows[0]["monotonic"] > start_s or
            valid_rows[-1]["monotonic"] < finish_s):
        raise MeasurementError("thermal window does not bracket peeltiming")
    first_line = first_search + first
    last_line = first_search + last
    artifact = header + b"".join(lines[first_line:last_line + 1])
    summary = {
        "schema": PAIRED_THERMAL_SCHEMA,
        "csv_sha256": hashlib.sha256(artifact).hexdigest(),
        "csv_window": artifact.decode("ascii"),
        "prelaunch_row": prelaunch_row.decode("ascii"),
        "rows": len(window),
        "valid_rows": len(valid_rows),
        "missing_read_rows": sum(
            row["tolerated_missing"] for row in window),
        "invalid_rows": sum(not row["valid"] for row in window),
        "first_monotonic_s": float(valid_rows[0]["monotonic"]),
        "last_monotonic_s": float(valid_rows[-1]["monotonic"]),
        "cpu_tctl_max_c": float(max(row["cpu"] for row in window)),
        "dimm_max_c": float(max(
            value for row in window for value in row["dimms"])),
        "edac_ce_max": max(row["ce"] for row in window),
        "edac_ue_max": max(row["ue"] for row in window),
    }
    summary["summary_sha256"] = _canonical_json_sha256(summary)
    return summary


def _finalize_paired_context(frozen_context, stdout, capture):
    """Finalize thermal evidence through the same retained prelaunch file."""
    start_ns, finish_ns = _peeltiming_run_bounds(stdout)
    _validate_runtime_identity(frozen_context["bound"])
    interval = (
        frozen_context["bound"]["thermal_sampling_interval_ms"] / 1000.0)
    deadline = time.monotonic() + 2.0 + 2.0 * interval
    last_error = None
    while True:
        try:
            payload = _snapshot_capture(capture)
            thermal = _thermal_window_from_csv(
                payload, start_ns, finish_ns, capture,
                sampling_interval_ms=
                    frozen_context["bound"]["thermal_sampling_interval_ms"])
            result = {
                "schema": frozen_context["schema"],
                "bound": copy.deepcopy(frozen_context["bound"]),
                "thermal": thermal,
            }
            _validate_runtime_identity(result["bound"])
            _validate_paired_context(
                result,
                cache_state=result["bound"]["cache_state"],
                evict_bytes=result["bound"]["evict_bytes"],
                started_ns=start_ns, finished_ns=finish_ns)
            return result
        except _ThermalCoveragePending as error:
            last_error = error
        if time.monotonic() >= deadline:
            raise MeasurementError(
                f"could not finalize peeltiming thermal evidence: {last_error}")
        time.sleep(0.1)


def _beta_continued_fraction(a, b, x):
    """Numerical Recipes continued fraction for the regularized beta."""
    tiny = sys.float_info.min * 1e3
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    result = d
    for iteration in range(1, 401):
        twice = 2 * iteration
        aa = (
            iteration * (b - iteration) * x /
            ((qam + twice) * (a + twice))
        )
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        result *= d * c
        aa = -(
            (a + iteration) * (qab + iteration) * x /
            ((a + twice) * (qap + twice))
        )
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        result *= delta
        if abs(delta - 1.0) <= 4.0 * sys.float_info.epsilon:
            return result
    raise MeasurementError("Student-t critical calculation did not converge")


def _regularized_beta(a, b, x):
    if not 0.0 <= x <= 1.0:
        raise MeasurementError("invalid beta argument")
    if x in (0.0, 1.0):
        return x
    front = math.exp(
        math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) +
        a * math.log(x) + b * math.log1p(-x)
    )
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _beta_continued_fraction(a, b, x) / a
    return 1.0 - (
        front * _beta_continued_fraction(b, a, 1.0 - x) / b
    )


def _student_t_cdf(value, degrees_freedom):
    if degrees_freedom < 1 or not math.isfinite(value):
        raise MeasurementError("invalid Student-t argument")
    if value == 0.0:
        return 0.5
    x = degrees_freedom / (degrees_freedom + value * value)
    tail = 0.5 * _regularized_beta(
        degrees_freedom / 2.0, 0.5, x)
    return 1.0 - tail if value > 0.0 else tail


def _student_t_critical_95(degrees_freedom):
    """Return the deterministic two-sided 95% Student-t critical."""
    if type(degrees_freedom) is not int or degrees_freedom < 1:
        raise MeasurementError("invalid Student-t degrees of freedom")
    lo, hi = 0.0, 16.0
    for _ in range(80):
        mid = (lo + hi) / 2.0
        if _student_t_cdf(mid, degrees_freedom) < 0.975:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def _paired_t_interval(values):
    if not isinstance(values, list) or len(values) < 4:
        raise MeasurementError("paired timing needs at least four replicates")
    if any(type(value) is not float or not math.isfinite(value)
           for value in values):
        raise MeasurementError("paired timing contains invalid log costs")
    count = len(values)
    mean = sum(values) / count
    variance = sum((value - mean) ** 2 for value in values) / (count - 1)
    half = _student_t_critical_95(count - 1) * math.sqrt(variance / count)
    return mean, mean - half, mean + half


def _paired_practical_log_bounds(required_margin):
    """Return ratio-consistent log bounds for a declared speed margin."""
    if (isinstance(required_margin, bool) or
            not isinstance(required_margin, (int, float)) or
            not math.isfinite(float(required_margin)) or
            not 0.0 <= required_margin <= 1.0):
        raise MeasurementError("invalid paired practical margin")
    factor = 1.0 + float(required_margin)
    return math.log(1.0 / factor), math.log(factor)


def _percentile_integer(values, fraction):
    if not values:
        return 0.0
    ordered = sorted(values)
    index = int(fraction * (len(ordered) - 1) + 0.5)
    return float(ordered[index])


def _summarize_paired_arm(rows, replicates, block_count, role):
    """Summarize one role, counting recovery once per measured replicate."""
    overheads = []
    failures = 0
    elapsed_ns = 0
    source_bytes = 0
    for replicate in range(replicates):
        selected = [
            row for row in rows
            if row["measured"] == 1 and row["replicate_measured"] == replicate
            and row["role"] == role
        ]
        if len(selected) != 4:
            raise MeasurementError(
                f"peeltiming role {role!r} is unbalanced")
        recoveries = {
            (row["recovery_result"], row["recovery_ok"], row["arm_overhead"])
            for row in selected
        }
        if len(recoveries) != 1:
            raise MeasurementError(
                f"peeltiming role {role!r} has unstable recovery metadata")
        unused_result, recovery_ok, overhead = next(iter(recoveries))
        if recovery_ok:
            overheads.append(overhead)
            for row in selected:
                elapsed_ns += row["elapsed_ns"]
                source_bytes += row["source_bytes"] * row["inner_reps"]
        else:
            failures += 1
    if failures == replicates:
        overheads = []
    mean = sum(overheads) / len(overheads) if overheads else 0.0
    variance = (
        sum((value - mean) ** 2 for value in overheads) / len(overheads)
        if overheads else 0.0
    )
    solve_mbps = (
        source_bytes * 1000.0 / elapsed_ns if elapsed_ns > 0 else 0.0
    )
    return PairedArmSummary(
        trials=replicates,
        fail=failures,
        oh_mean=float(mean),
        oh_sd=float(math.sqrt(variance)),
        oh50=_percentile_integer(overheads, 0.50),
        oh95=_percentile_integer(overheads, 0.95),
        oh99=_percentile_integer(overheads, 0.99),
        oh_max=float(max(overheads) if overheads else 0),
        solve_mbps=float(solve_mbps),
    )


def _expected_peeltiming_roles(panel, replicate, label):
    swapped = (replicate & 1) != 0
    if panel == "candidate_control":
        if label == "A":
            return "candidate" if swapped else "identity"
        return "identity" if swapped else "candidate"
    if panel != "identity_aa":
        raise MeasurementError(f"unknown peeltiming panel {panel!r}")
    if label == "A":
        return "identity_b" if swapped else "identity_a"
    return "identity_a" if swapped else "identity_b"


def _parse_peeltiming_row(fields, label):
    if len(fields) != len(PEELTIMING_COLUMNS):
        raise MeasurementError(f"{label} has the wrong column count")
    raw = dict(zip(PEELTIMING_COLUMNS, fields))
    uint32_names = {
        "replicate", "measured", "slot", "pair", "label_swap",
        "common_overhead", "recovery_ok", "timing_eligible", "inner_reps",
        "saturated",
        "inactivated", "binary_def", "heavy_gain",
    }
    uint64_names = {
        "elapsed_ns", "block_xors", "block_muladds", "build_ns_sum",
        "peel_ns_sum", "project_ns_sum", "residual_ns_sum",
        "backsub_ns_sum", "source_bytes", "packet_payload_bytes",
        "intermediate_bytes",
    }
    int32_names = {
        "recovery_result", "preflight_result", "timing_result",
        "outcome_stable", "cpu_before", "cpu_after", "cpu_migrated",
        "fault_contaminated",
    }
    int64_names = {"minflt_delta", "majflt_delta"}
    row = {}
    for name in uint32_names:
        row[name] = _parse_decimal_integer(
            raw[name], f"{label}.{name}", minimum=0, maximum=0xffffffff)
    for name in uint64_names:
        row[name] = _parse_decimal_integer(
            raw[name], f"{label}.{name}",
            minimum=0, maximum=0xffffffffffffffff)
    for name in int32_names:
        row[name] = _parse_decimal_integer(
            raw[name], f"{label}.{name}",
            minimum=-0x80000000, maximum=0x7fffffff)
    for name in int64_names:
        row[name] = _parse_decimal_integer(
            raw[name], f"{label}.{name}",
            minimum=0, maximum=0x7fffffffffffffff)
    row["arm_overhead"] = _parse_decimal_integer(
        raw["arm_overhead"], f"{label}.arm_overhead",
        minimum=-1, maximum=0x7fffffff)
    row["construction_seed"] = _parse_decimal_integer(
        raw["construction_seed"], f"{label}.construction_seed",
        minimum=0, maximum=0xffffffffffffffff)
    row["loss_seed"] = _parse_decimal_integer(
        raw["loss_seed"], f"{label}.loss_seed",
        minimum=0, maximum=0xffffffffffffffff)
    for name in ("panel", "label", "role"):
        row[name] = raw[name]
    if not _is_sha256(raw["trace_sha256"]):
        raise MeasurementError(f"{label}.trace_sha256 is invalid")
    row["trace_sha256"] = raw["trace_sha256"]
    return row


def validate_peeltiming_dimensions(
        *, block_count, block_bytes, target_profile, seed_policy,
        construction_seed, loss, loss_seed, schedule, warmup_replicates,
        replicates, inner_reps, max_overhead, cache_state, evict_bytes,
        required_margin):
    """Mirror every native peeltiming request-domain constraint."""
    if target_profile != TARGET_PROFILE:
        raise ValueError(
            f"target_profile must be exactly {TARGET_PROFILE!r}")
    if seed_policy != TARGET_SEED_POLICY:
        raise ValueError(f"seed_policy must be exactly {TARGET_SEED_POLICY!r}")

    def working_set_bytes():
        staircase = _dispatch_staircase_count(block_count)
        width = block_count + staircase + 4 + 12
        retained_rows = min(width, 4096)
        blocks = (
            block_count + 3 * (block_count + max_overhead) +
            2 * width + 2 * retained_rows + 96
        )
        return (
            blocks * block_bytes +
            block_count * 4096 +
            evict_bytes
        )

    if (not isinstance(schedule, str) or schedule not in TARGET_SCHEDULES or
            not valid_loss_rate(loss) or
            any(
                type(seed) is not int or
                not 0 <= seed <= 0xffffffffffffffff
                for seed in (construction_seed, loss_seed)) or
            type(block_count) is not int or
            not 2 <= block_count <= 64000 or
            type(block_bytes) is not int or
            not 2 <= block_bytes <= _COMPARE_BLOCK_BYTES_MAX or
            block_bytes % 2 != 0 or
            type(warmup_replicates) is not int or
            warmup_replicates < 0 or warmup_replicates % 2 != 0 or
            type(replicates) is not int or
            replicates < 4 or replicates % 2 != 0 or
            warmup_replicates + replicates > 10000 or
            type(inner_reps) is not int or not 1 <= inner_reps <= 1024 or
            type(max_overhead) is not int or
            not 0 <= max_overhead <= min(4096, block_count + 512) or
            cache_state not in ("warm", "cold") or
            type(evict_bytes) is not int or
            not 4096 <= evict_bytes <= PEELTIMING_EVICT_BYTES_MAX or
            working_set_bytes() > PEELTIMING_WORKING_SET_BYTES_MAX or
            (cache_state == "cold" and inner_reps != 1) or
            (warmup_replicates + replicates) * 16 * inner_reps *
                (block_count + max_overhead) > 1_000_000_000_000 or
            isinstance(required_margin, bool) or
            not isinstance(required_margin, (int, float)) or
            not math.isfinite(float(required_margin)) or
            (float(required_margin) == 0.0 and
             math.copysign(1.0, float(required_margin)) < 0.0) or
            not 0.0 <= required_margin <= 1.0):
        raise ValueError("invalid peeltiming dimensions or measurement policy")


def require_distinct_paths(
        first_path, second_path, first_label="output",
        second_label="protected input"):
    """Fail closed unless two pathnames demonstrably name distinct objects."""
    try:
        if os.path.realpath(first_path) == os.path.realpath(second_path):
            raise MeasurementError(
                f"{first_label} must not replace {second_label}")
        if (os.path.lexists(first_path) and os.path.lexists(second_path) and
                os.path.samefile(first_path, second_path)):
            raise MeasurementError(
                f"{first_label} must not replace {second_label}")
    except OSError as error:
        raise MeasurementError(
            f"could not establish that {first_label} is distinct from "
            f"{second_label}: {error}")


def require_distinct_from_source_provenance(path, path_label, generator):
    """Reject a result path that aliases any source-identity input."""
    git_directories, git_files = _git_provenance_paths()
    try:
        lexical_path = os.path.abspath(path)
        lexical_root = os.path.abspath(REPO_ROOT)
        resolved_path = os.path.realpath(path)
        resolved_root = os.path.realpath(REPO_ROOT)
        lexical_is_source = (
            lexical_path.endswith(_SOURCE_EXTENSIONS) or
            os.path.basename(lexical_path) in (
                "CMakeLists.txt", ".gitignore"))
        resolved_is_source = (
            resolved_path.endswith(_SOURCE_EXTENSIONS) or
            os.path.basename(resolved_path) in (
                "CMakeLists.txt", ".gitignore"))
        creates_source_input = (
            (
                lexical_is_source and
                os.path.commonpath(
                    (lexical_path, lexical_root)) == lexical_root
            ) or (
                resolved_is_source and
                os.path.commonpath(
                    (resolved_path, resolved_root)) == resolved_root
            )
        )
    except (OSError, ValueError) as error:
        raise MeasurementError(
            f"could not establish whether {path_label} replaces a "
            f"source-provenance input: {error}")
    if creates_source_input:
        raise MeasurementError(
            f"{path_label} must not create or replace a "
            f"source-provenance input")
    for label, directory in git_directories:
        try:
            lexical_directory = os.path.abspath(directory)
            resolved_directory = os.path.realpath(directory)
            aliases_directory = (
                os.path.commonpath(
                    (lexical_path, lexical_directory)) == lexical_directory or
                os.path.commonpath(
                    (resolved_path, resolved_directory)) == resolved_directory
            )
        except (OSError, ValueError) as error:
            raise MeasurementError(
                f"could not establish that {path_label} is outside "
                f"{label}: {error}")
        if aliases_directory:
            raise MeasurementError(
                f"{path_label} must not replace {label}")
    for label, metadata_path in git_files:
        require_distinct_paths(path, metadata_path, path_label, label)

    generator_path = (
        generator if os.path.isabs(generator)
        else os.path.join(REPO_ROOT, generator)
    )
    require_distinct_paths(
        path, generator_path, path_label, "the generator source")
    generator_absolute = os.path.abspath(generator_path)
    for relative, source_path in _source_provenance_files():
        if os.path.abspath(source_path) == generator_absolute:
            continue
        require_distinct_paths(
            path, source_path, path_label,
            f"source-provenance input {relative!r}")


def parse_peeltiming_output(
        stdout, *, block_count, block_bytes, candidate_pmf, degree_scale,
        native_profile, target_profile, seed_policy, construction_seed,
        loss, loss_seed, schedule, warmup_replicates, replicates, inner_reps,
        max_overhead, cache_state, evict_bytes, context, required_margin):
    """Parse, balance, and independently aggregate one native timing stream."""
    validate_peeltiming_dimensions(
        block_count=block_count,
        block_bytes=block_bytes,
        target_profile=target_profile,
        seed_policy=seed_policy,
        construction_seed=construction_seed,
        loss=loss,
        loss_seed=loss_seed,
        schedule=schedule,
        warmup_replicates=warmup_replicates,
        replicates=replicates,
        inner_reps=inner_reps,
        max_overhead=max_overhead,
        cache_state=cache_state,
        evict_bytes=evict_bytes,
        required_margin=required_margin,
    )
    _validate_paired_degree_scale(degree_scale)
    if not isinstance(stdout, str) or not stdout:
        raise MeasurementError("peeltiming returned no output")
    try:
        stdout.encode("ascii")
    except UnicodeEncodeError:
        raise MeasurementError("peeltiming output is not ASCII")
    lines = stdout.splitlines(keepends=True)
    if (not lines or any(not line.endswith("\n") or line.endswith("\r\n")
                         for line in lines)):
        raise MeasurementError(
            "peeltiming output is incomplete or uses noncanonical newlines")
    if any(line == "\n" for line in lines):
        raise MeasurementError("peeltiming output contains a blank record")

    expected_rows = 16 * (warmup_replicates + replicates)
    if len(lines) != expected_rows + 4:
        raise MeasurementError(
            f"peeltiming emitted {len(lines) - 4} rows; expected "
            f"{expected_rows}")
    manifest_raw = _parse_ordered_kv_line(
        lines[0][:-1], "# peeltiming,", PEELTIMING_MANIFEST_FIELDS,
        "peeltiming manifest")
    semantic_raw = _parse_ordered_kv_line(
        lines[1][:-1], "# peel_semantic,", PEELTIMING_SEMANTIC_FIELDS,
        "peeltiming semantic")
    try:
        header = next(csv.reader([lines[2][:-1]], strict=True))
    except (csv.Error, StopIteration) as error:
        raise MeasurementError(f"malformed peeltiming header: {error}")
    if (lines[2][:-1] != ",".join(header) or
            tuple(header) != PEELTIMING_COLUMNS):
        raise MeasurementError("peeltiming header is missing or reordered")

    done_fields = (
        "complete", "rows", "finished_monotonic_ns", "stream_sha256")
    done = _parse_ordered_kv_line(
        lines[-1][:-1], "# peeltiming_done,", done_fields,
        "peeltiming completion")
    if (done["complete"] != "1" or
            _parse_decimal_integer(
                done["rows"], "peeltiming completion rows", minimum=0) !=
            expected_rows or
            not _is_sha256(done["stream_sha256"])):
        raise MeasurementError("peeltiming completion receipt is invalid")
    finished_ns = _parse_decimal_integer(
        done["finished_monotonic_ns"],
        "peeltiming finished_monotonic_ns", minimum=1,
        maximum=0xffffffffffffffff)
    done_prefix = lines[-1][:-1].rsplit("stream_sha256=", 1)[0]
    done_prefix += "stream_sha256="
    stream_sha256 = hashlib.sha256(
        ("".join(lines[:-1]) + done_prefix).encode("ascii")).hexdigest()
    if done["stream_sha256"] != stream_sha256:
        raise MeasurementError("peeltiming stream SHA-256 does not match")

    if (not isinstance(native_profile, StockProfile) or
            native_profile.block_count != block_count):
        raise ValueError("native_profile does not match peeltiming K")
    candidate_digest = _paired_pmf_sha256(candidate_pmf)
    identity_digest = _paired_pmf_sha256(native_profile.pmf)
    requested_scale = (
        "identity" if degree_scale is None else
        _canonical_staircase_scale_spec(degree_scale)
    )
    candidate_is_exact_identity = (
        _native_peel_cdf_signature(candidate_pmf, block_count) ==
        _native_peel_cdf_signature(native_profile.pmf, block_count) and
        _scale_realizes_legacy_edges(
            degree_scale, block_count, native_profile.staircase,
            native_profile.source_hits)
    )
    context_sha256 = _validate_paired_context(
        context, cache_state=cache_state, evict_bytes=evict_bytes)
    manifest_integer_fields = {
        "K", "bb", "S", "D2", "H", "warmup_replicates", "replicates",
        "slots_per_panel", "panels_per_replicate", "inner_reps",
        "max_overhead", "evict_bytes", "payload_alignment", "prefault",
        "expected_rows", "started_monotonic_ns",
    }
    manifest = dict(manifest_raw)
    for name in manifest_integer_fields:
        manifest[name] = _parse_decimal_integer(
            manifest_raw[name], f"peeltiming manifest.{name}", minimum=0,
            maximum=(
                0xffffffffffffffff
                if name in ("evict_bytes", "started_monotonic_ns")
                else 0xffffffff))
    manifest["construction_seed_base"] = _parse_decimal_integer(
        manifest_raw["construction_seed_base"],
        "peeltiming manifest.construction_seed_base",
        minimum=0, maximum=0xffffffffffffffff)
    manifest["loss_seed_base"] = _parse_decimal_integer(
        manifest_raw["loss_seed_base"], "peeltiming manifest.loss_seed_base",
        minimum=0, maximum=0xffffffffffffffff)
    for name in (
            "loss", "candidate_scale_effective", "identity_scale_effective",
            "required_margin"):
        manifest[name] = _parse_finite_float(
            manifest_raw[name], f"peeltiming manifest.{name}", minimum=0.0,
            maximum=1.0 if name in ("loss", "required_margin") else None)
    if finished_ns <= manifest["started_monotonic_ns"]:
        raise MeasurementError("peeltiming completion clock did not advance")
    if (manifest_raw["schema"] != _PEELTIMING_SCHEMA or
            manifest_raw["target_profile"] != target_profile or
            manifest_raw["seed_policy"] != seed_policy or
            manifest_raw["contract_id"] != TARGET_CONTRACT["contract_id"] or
            manifest["K"] != block_count or manifest["bb"] != block_bytes or
            manifest["S"] != native_profile.staircase or
            manifest["D2"] != native_profile.dense_rows or
            manifest["H"] != native_profile.heavy_rows or
            manifest["construction_seed_base"] != construction_seed or
            manifest_raw["construction_seed_derivation"] !=
                PAIRED_CONSTRUCTION_SEED_DERIVATION or
            manifest_raw["semantic_seed_derivation"] !=
                "base-plus-attempt-mod2^64-skip-measured-alias-v1" or
            manifest["loss"] != float(loss) or
            manifest_raw["loss"] != f"{float(loss):.17g}" or
            manifest["loss_seed_base"] != loss_seed or
            manifest_raw["loss_seed_derivation"] !=
                PAIRED_LOSS_SEED_DERIVATION or
            manifest_raw["message_seed_policy"] !=
                "replicate-loss-seed-v1" or
            manifest_raw["schedule"] != schedule or
            manifest_raw["loss_model"] != PAIRED_LOSS_MODEL or
            manifest_raw["trace_encoding"] !=
                "wirehair-wh2-peeltiming-loss-trace-v1" or
            manifest_raw["panels"] != _PEELTIMING_PANELS or
            manifest_raw["candidate_pmf_sha256"] != candidate_digest or
            manifest_raw["candidate_pmf_encoding"] !=
                _PEELTIMING_PMF_ENCODING or
            manifest_raw["candidate_scale_requested"] != requested_scale or
            manifest_raw["identity_pmf_sha256"] != identity_digest or
            manifest_raw["identity_pmf_encoding"] !=
                _PEELTIMING_IDENTITY_PMF_ENCODING or
            manifest["warmup_replicates"] != warmup_replicates or
            manifest["replicates"] != replicates or
            manifest["slots_per_panel"] != 8 or
            manifest["panels_per_replicate"] != 2 or
            manifest_raw["order"] != PAIRED_ORDER or
            manifest_raw["label_swap"] != _PEELTIMING_LABEL_SWAP or
            manifest["inner_reps"] != inner_reps or
            manifest["max_overhead"] != max_overhead or
            manifest_raw["cache_state"] != cache_state or
            manifest["evict_bytes"] != evict_bytes or
            manifest["payload_alignment"] != 64 or
            manifest["prefault"] != 1 or
            manifest_raw["cpu_affinity_policy"] !=
                "first-allowed-affinity-v1" or
            manifest_raw["payload"] !=
                "common-source-role-encoded-consistent-rhs-v1" or
            manifest_raw["timing_scope"] != "decoder_solve_only" or
            manifest_raw["timing_prefix"] !=
                "common-max-recovery-overhead-v1" or
            manifest_raw["recovery_scope"] !=
                "full_encode_decode_recover_memcmp" or
            manifest_raw["weak_seed_policy"] !=
                "balanced-timing-ineligible-rows-v1" or
            manifest_raw["hook_path"] !=
                "explicit-tls-peel-pmf-and-scale-v1" or
            manifest_raw["no_override_scope"] !=
                "untimed-semantic-only" or
            manifest_raw["system_build"] !=
                "per_replicate_per_scale_outside_timer" or
            manifest_raw["startup_amortization"] !=
                "per-replicate-build-and-recovery-preflight-excluded-v1" or
        manifest_raw["slot_prewarm"] !=
                "validated-plus-conditioning-matching-role-solves-"
                "same-cpu-before-cache-v1" or
            manifest_raw["context_sha256"] != context_sha256 or
            manifest_raw["uncertainty"] != PAIRED_UNCERTAINTY_PROTOCOL or
            manifest["required_margin"] != required_margin or
            manifest_raw["required_margin"] !=
                f"{float(required_margin):.17g}" or
            manifest_raw["margin_rule"] !=
                "upper-log-cost-lt-negative-required-margin-and-aa-floor-v1" or
            manifest_raw["clock_domain"] != _PEELTIMING_CLOCK_DOMAIN or
            manifest_raw["stream_hash_scope"] !=
                "body-plus-done-prefix-v1" or
            manifest["expected_rows"] != expected_rows or
            manifest["started_monotonic_ns"] < 1):
        raise MeasurementError(
            "peeltiming manifest does not match the requested experiment")
    if requested_scale == "identity":
        expected_identity_scale = (
            block_count *
            min(native_profile.source_hits, native_profile.staircase) /
            native_profile.staircase
        )
        if (manifest["identity_scale_effective"] !=
                expected_identity_scale or
                manifest["candidate_scale_effective"] !=
                expected_identity_scale or
                manifest_raw["identity_scale_effective"] !=
                    f"{expected_identity_scale:.17g}" or
                manifest_raw["candidate_scale_effective"] !=
                    f"{expected_identity_scale:.17g}"):
            raise MeasurementError(
                "identity candidate did not use the identity scale hook")
    else:
        expected_identity_scale = (
            block_count *
            min(native_profile.source_hits, native_profile.staircase) /
            native_profile.staircase
        )
        if (manifest["identity_scale_effective"] !=
                expected_identity_scale or
                manifest["candidate_scale_effective"] !=
                float(requested_scale) or
                manifest_raw["identity_scale_effective"] !=
                    f"{expected_identity_scale:.17g}" or
                manifest_raw["candidate_scale_effective"] !=
                    requested_scale):
            raise MeasurementError(
                "candidate or identity scale receipt is incorrect")
    _validate_paired_context(
        context, cache_state=cache_state, evict_bytes=evict_bytes,
        started_ns=manifest["started_monotonic_ns"], finished_ns=finished_ns)

    semantic = dict(semantic_raw)
    for name in (
            "timed", "seed_attempt", "seed_attempt_cap",
            "nohook_result", "nohook_recovery_ok",
            "nohook_overhead", "identity_result", "identity_recovery_ok",
            "identity_overhead", "system_equal", "packet_rows_equal",
            "payload_equal", "nohook_direct_result",
            "identity_direct_result", "nohook_intermediate_bytes",
            "identity_intermediate_bytes", "solve_equal",
            "full_recovery_equal", "pass"):
        semantic[name] = _parse_decimal_integer(
            semantic_raw[name], f"peeltiming semantic.{name}", minimum=0)
    semantic["construction_seed"] = _parse_decimal_integer(
        semantic_raw["construction_seed"],
        "peeltiming semantic.construction_seed",
        minimum=0, maximum=0xffffffffffffffff)
    semantic["loss_seed"] = _parse_decimal_integer(
        semantic_raw["loss_seed"], "peeltiming semantic.loss_seed",
        minimum=0, maximum=0xffffffffffffffff)
    semantic["identity_scale_effective"] = _parse_finite_float(
        semantic_raw["identity_scale_effective"],
        "peeltiming semantic.identity_scale_effective", minimum=0.0)
    semantic_hash_fields = (
        "trace_sha256", "message_sha256", "identity_pmf_sha256",
        "nohook_system_sha256", "identity_system_sha256",
        "nohook_packet_rows_sha256", "identity_packet_rows_sha256",
        "nohook_payload_sha256", "identity_payload_sha256",
        "nohook_solve_sha256", "identity_solve_sha256",
    )
    if any(not _is_sha256(semantic_raw[name])
           for name in semantic_hash_fields):
        raise MeasurementError("peeltiming semantic digest is invalid")
    if (semantic["seed_attempt_cap"] != 256 or
            semantic["seed_attempt"] >= semantic["seed_attempt_cap"]):
        raise MeasurementError("peeltiming semantic seed search is invalid")
    expected_semantic_construction_seed = (
        construction_seed + semantic["seed_attempt"]) & 0xffffffffffffffff
    expected_semantic_loss_seed = loss_seed
    semantic_intermediate_bytes = (
        block_count + native_profile.staircase +
        native_profile.dense_rows + native_profile.heavy_rows
    ) * block_bytes
    semantic_terminal_equal = (
        semantic["nohook_result"] == 0 and
        semantic["identity_result"] == 0 and
        semantic["nohook_recovery_ok"] == 1 and
        semantic["identity_recovery_ok"] == 1 and
            semantic["identity_overhead"] == semantic["nohook_overhead"] and
            semantic["identity_overhead"] <= max_overhead
    )
    if (semantic["timed"] != 0 or
            semantic["construction_seed"] !=
                expected_semantic_construction_seed or
            semantic["loss_seed"] != expected_semantic_loss_seed or
            semantic_raw["identity_pmf_sha256"] != identity_digest or
            semantic_raw["identity_pmf_encoding"] !=
                _PEELTIMING_IDENTITY_PMF_ENCODING or
            semantic["identity_scale_effective"] !=
                manifest["identity_scale_effective"] or
            semantic_raw["identity_scale_effective"] !=
                f"{expected_identity_scale:.17g}" or
            not semantic_terminal_equal or
            semantic_raw["nohook_system_sha256"] !=
                semantic_raw["identity_system_sha256"] or
            semantic["system_equal"] != 1 or
            semantic_raw["nohook_packet_rows_sha256"] !=
                semantic_raw["identity_packet_rows_sha256"] or
            semantic["packet_rows_equal"] != 1 or
            semantic_raw["nohook_payload_sha256"] !=
                semantic_raw["identity_payload_sha256"] or
            semantic["payload_equal"] != 1 or
            semantic["nohook_direct_result"] != 0 or
            semantic["identity_direct_result"] != 0 or
            semantic["nohook_intermediate_bytes"] !=
                semantic_intermediate_bytes or
            semantic["identity_intermediate_bytes"] !=
                semantic_intermediate_bytes or
            semantic_raw["nohook_solve_sha256"] !=
                semantic_raw["identity_solve_sha256"] or
            semantic["solve_equal"] != 1 or
            semantic["full_recovery_equal"] != 1 or
            semantic["pass"] != 1):
        raise MeasurementError(
            "peeltiming no-override semantic identity check failed")

    rows = []
    total_replicates = warmup_replicates + replicates
    successful_intermediate_bytes = (
        block_count + native_profile.staircase +
        native_profile.dense_rows + native_profile.heavy_rows
    ) * block_bytes
    system_width = (
        block_count + native_profile.staircase +
        native_profile.dense_rows + native_profile.heavy_rows
    )
    seen_trace = {}
    for row_index, line in enumerate(lines[3:-1]):
        try:
            fields = next(csv.reader([line[:-1]], strict=True))
        except (csv.Error, StopIteration) as error:
            raise MeasurementError(
                f"malformed peeltiming row {row_index}: {error}")
        if line[:-1] != ",".join(fields):
            raise MeasurementError(
                f"peeltiming row {row_index} is not canonical CSV")
        row = _parse_peeltiming_row(fields, f"peeltiming row {row_index}")
        expected_replicate = row_index // 16
        within = row_index % 16
        expected_panel = (
            "candidate_control" if within < 8 else "identity_aa")
        expected_slot = within % 8
        expected_label = PAIRED_ORDER[expected_slot]
        expected_measured = (
            1 if expected_replicate >= warmup_replicates else 0)
        expected_swap = expected_replicate & 1
        expected_role = _expected_peeltiming_roles(
            expected_panel, expected_replicate, expected_label)
        timing_eligible = row["timing_eligible"] == 1
        timed_contract = (
            row["recovery_result"] == 0 and
            row["preflight_result"] == row["recovery_result"] and
            row["timing_result"] == row["preflight_result"] and
            row["outcome_stable"] == 1 and row["elapsed_ns"] >= 1 and
            row["inner_reps"] == inner_reps and
            row["saturated"] == 0 and row["cpu_before"] >= 0 and
            row["cpu_after"] >= 0 and row["cpu_migrated"] == 0 and
            row["cpu_before"] == row["cpu_after"] and
            row["cpu_before"] == context["bound"]["cpu_affinity"][0] and
            row["minflt_delta"] == 0 and row["majflt_delta"] == 0 and
            row["fault_contaminated"] == 0 and
            0 <= row["heavy_gain"] <= row["binary_def"] <=
                row["inactivated"] <= system_width and
            row["heavy_gain"] <= native_profile.heavy_rows and
            (row["preflight_result"] != 0 or
             row["heavy_gain"] == row["binary_def"]) and
            (row["preflight_result"] != 1 or
             row["heavy_gain"] < row["binary_def"]) and
            row["packet_payload_bytes"] ==
                (block_count + row["common_overhead"]) * block_bytes and
            row["intermediate_bytes"] == (
                successful_intermediate_bytes
                if row["preflight_result"] == 0 else 0) and
            sum(row[name] for name in (
                "build_ns_sum", "peel_ns_sum", "project_ns_sum",
                "residual_ns_sum", "backsub_ns_sum")) <= row["elapsed_ns"]
        )
        skipped_contract = (
            row["preflight_result"] == -1 and
            row["timing_result"] == -1 and
            row["outcome_stable"] == 0 and row["elapsed_ns"] == 0 and
            row["inner_reps"] == 0 and row["saturated"] == 0 and
            row["cpu_before"] == -1 and row["cpu_after"] == -1 and
            row["cpu_migrated"] == -1 and
            row["minflt_delta"] == 0 and row["majflt_delta"] == 0 and
            row["fault_contaminated"] == 0 and
            all(row[name] == 0 for name in (
                "inactivated", "binary_def", "heavy_gain", "block_xors",
                "block_muladds", "build_ns_sum", "peel_ns_sum",
                "project_ns_sum", "residual_ns_sum", "backsub_ns_sum",
                "packet_payload_bytes", "intermediate_bytes"))
        )
        if (row["replicate"] != expected_replicate or
                row["measured"] != expected_measured or
                row["panel"] != expected_panel or
                row["slot"] != expected_slot or
                row["pair"] != expected_slot // 2 or
                row["label"] != expected_label or
                row["role"] != expected_role or
                row["label_swap"] != expected_swap or
                row["common_overhead"] > max_overhead or
                (row["recovery_ok"] == 1 and not
                 0 <= row["arm_overhead"] <= row["common_overhead"]) or
                (row["recovery_ok"] == 0 and
                 row["arm_overhead"] != -1) or
                row["recovery_result"] not in (0, 1, 3, 4, 7) or
                row["recovery_ok"] not in (0, 1) or
                row["recovery_ok"] != int(
                    row["recovery_result"] == 0) or
                row["timing_eligible"] not in (0, 1) or
                (timing_eligible and not timed_contract) or
                (not timing_eligible and not skipped_contract) or
                row["source_bytes"] != block_count * block_bytes):
            raise MeasurementError(
                f"peeltiming row {row_index} violates its balance contract")
        expected_loss_seed = (
            loss_seed ^
            (((expected_replicate + 1) * 0x9e3779b97f4a7c15) &
             0xffffffffffffffff)
        )
        expected_construction_seed = (
            construction_seed ^
            (((expected_replicate + 1) * 0xd1b54a32d192ed03) &
             0xffffffffffffffff)
        )
        if row["construction_seed"] != expected_construction_seed:
            raise MeasurementError(
                f"peeltiming row {row_index} has the wrong derived "
                "construction seed")
        if row["loss_seed"] != expected_loss_seed:
            raise MeasurementError(
                f"peeltiming row {row_index} has the wrong derived loss seed")
        trace = seen_trace.setdefault(
            expected_replicate,
            (row["construction_seed"], row["loss_seed"],
             row["trace_sha256"], row["common_overhead"]))
        if trace != (
                row["construction_seed"], row["loss_seed"],
                row["trace_sha256"], row["common_overhead"]):
            raise MeasurementError(
                f"peeltiming replicate {expected_replicate} changed trace")
        row["replicate_measured"] = expected_replicate - warmup_replicates
        rows.append(row)
    if (len(seen_trace) != total_replicates or
            len({value[0] for value in seen_trace.values()}) !=
                total_replicates or
            len({value[1] for value in seen_trace.values()}) !=
                total_replicates or
            len({value[2] for value in seen_trace.values()}) !=
                total_replicates):
        raise MeasurementError(
            "peeltiming replicate seeds or traces are missing or reused")
    if semantic["trace_sha256"] in {
            value[2] for value in seen_trace.values()}:
        raise MeasurementError(
            "peeltiming semantic trace aliases a timed replicate")
    if (semantic["construction_seed"] in {
            value[0] for value in seen_trace.values()} or
            semantic["loss_seed"] in {
                value[1] for value in seen_trace.values()}):
        raise MeasurementError(
            "peeltiming semantic seed aliases a timed replicate")
    for absolute in range(total_replicates):
        selected = [row for row in rows if row["replicate"] == absolute]
        role_rows = {
            role: [row for row in selected if row["role"] == role]
            for role in ("candidate", "identity", "identity_a", "identity_b")
        }
        if any(len(value) != 4 for value in role_rows.values()):
            raise MeasurementError(
                f"peeltiming replicate {absolute} is unbalanced")
        recovery = {}
        for role, values in role_rows.items():
            stable = {
                (row["arm_overhead"], row["recovery_result"],
                 row["recovery_ok"], row["preflight_result"],
                 row["timing_result"], row["inactivated"],
                 row["binary_def"], row["heavy_gain"], row["block_xors"],
                 row["block_muladds"])
                for row in values
            }
            if len(stable) != 1:
                raise MeasurementError(
                    f"peeltiming {role} metadata changed within replicate "
                    f"{absolute}")
            recovery[role] = next(iter(stable))
        if (recovery["identity"] != recovery["identity_a"] or
                recovery["identity"] != recovery["identity_b"]):
            raise MeasurementError(
                f"peeltiming identity arms differ in replicate {absolute}")
        if (candidate_is_exact_identity and
                recovery["candidate"] != recovery["identity"]):
            raise MeasurementError(
                f"peeltiming stock candidate/control metadata differs in "
                f"replicate {absolute}")
        eligibility = {
            row["timing_eligible"] for row in selected}
        expected_eligible = int(all(
            recovery[role][1] == 0 for role in recovery))
        if eligibility != {expected_eligible}:
            raise MeasurementError(
                f"peeltiming replicate {absolute} has inconsistent timing "
                "eligibility")
        if expected_eligible:
            for pair in range(4):
                candidate_cpu = next(
                    row["cpu_before"]
                    for row in role_rows["candidate"]
                    if row["pair"] == pair)
                identity_cpu = next(
                    row["cpu_before"]
                    for row in role_rows["identity"]
                    if row["pair"] == pair)
                aa_left_cpu = next(
                    row["cpu_before"]
                    for row in role_rows["identity_a"]
                    if row["pair"] == pair)
                aa_right_cpu = next(
                    row["cpu_before"]
                    for row in role_rows["identity_b"]
                    if row["pair"] == pair)
                if (candidate_cpu != identity_cpu or
                        aa_left_cpu != aa_right_cpu):
                    raise MeasurementError(
                        f"peeltiming pair {pair} in replicate {absolute} "
                        "ran on different CPUs")
        candidate_oh, unused_result, candidate_ok, *_ = recovery["candidate"]
        identity_oh, unused_result, identity_ok, *_ = recovery["identity"]
        expected_common = (
            max(candidate_oh, identity_oh)
            if candidate_ok and identity_ok else max_overhead
        )
        if seen_trace[absolute][3] != expected_common:
            raise MeasurementError(
                f"peeltiming replicate {absolute} has the wrong "
                "common overhead")

    replicate_receipts = []
    candidate_log_costs = []
    aa_log_costs = []
    recovery_regressions = 0
    recovery_improvements = 0
    both_failures = 0
    overhead_regressions = 0
    candidate_result_counts = {
        name: 0 for name in (
            "success", "need_more", "bad_dense_seed", "bad_peel_seed",
            "extra_insufficient")}
    identity_result_counts = dict(candidate_result_counts)
    result_names = {
        0: "success", 1: "need_more", 3: "bad_dense_seed",
        4: "bad_peel_seed", 7: "extra_insufficient",
    }
    timing_eligible_replicates = 0
    for measured_index in range(replicates):
        absolute = warmup_replicates + measured_index
        selected = [
            row for row in rows if row["replicate"] == absolute]
        role_rows = {
            role: [row for row in selected if row["role"] == role]
            for role in ("candidate", "identity", "identity_a", "identity_b")
        }
        if any(len(value) != 4 for value in role_rows.values()):
            raise MeasurementError(
                f"peeltiming measured replicate {measured_index} is unbalanced")
        role_meta = {}
        for role, values in role_rows.items():
            stable = {
                (row["arm_overhead"], row["recovery_result"],
                 row["recovery_ok"], row["preflight_result"],
                 row["timing_result"], row["inactivated"],
                 row["binary_def"], row["heavy_gain"], row["block_xors"],
                 row["block_muladds"])
                for row in values
            }
            if len(stable) != 1:
                raise MeasurementError(
                    f"peeltiming {role} metadata changed within a replicate")
            role_meta[role] = next(iter(stable))
        if role_meta["identity"] != role_meta["identity_a"] or \
                role_meta["identity"] != role_meta["identity_b"]:
            raise MeasurementError(
                "peeltiming identity arms are not equation/recovery identical")
        candidate_oh, candidate_result, candidate_ok, *_ = \
            role_meta["candidate"]
        identity_oh, identity_result, identity_ok, *_ = role_meta["identity"]
        candidate_result_counts[result_names[candidate_result]] += 1
        identity_result_counts[result_names[identity_result]] += 1
        expected_common_overhead = (
            max(candidate_oh, identity_oh)
            if candidate_ok and identity_ok else max_overhead
        )
        if seen_trace[absolute][3] != expected_common_overhead:
            raise MeasurementError(
                f"peeltiming measured replicate {measured_index} has "
                "the wrong common overhead")
        if identity_ok and not candidate_ok:
            recovery_regressions += 1
        if candidate_ok and not identity_ok:
            recovery_improvements += 1
        if not candidate_ok and not identity_ok:
            both_failures += 1
        if candidate_ok and identity_ok and candidate_oh > identity_oh:
            overhead_regressions += 1
        timing_eligible = selected[0]["timing_eligible"] == 1
        if any(
                (row["timing_eligible"] == 1) != timing_eligible
                for row in selected):
            raise MeasurementError(
                "peeltiming measured replicate changed timing eligibility")
        candidate_pair_logs = []
        aa_pair_logs = []
        if timing_eligible:
            timing_eligible_replicates += 1
            for pair in range(4):
                candidate_row = next(
                    row for row in role_rows["candidate"]
                    if row["pair"] == pair)
                identity_row = next(
                    row for row in role_rows["identity"]
                    if row["pair"] == pair)
                aa_left_row = next(
                    row for row in role_rows["identity_a"]
                    if row["pair"] == pair)
                aa_right_row = next(
                    row for row in role_rows["identity_b"]
                    if row["pair"] == pair)
                candidate_pair_logs.append(math.log(
                    candidate_row["elapsed_ns"] /
                    identity_row["elapsed_ns"]))
                aa_pair_logs.append(math.log(
                    aa_left_row["elapsed_ns"] /
                    aa_right_row["elapsed_ns"]))
        candidate_elapsed = sum(
            row["elapsed_ns"] for row in role_rows["candidate"])
        identity_elapsed = sum(
            row["elapsed_ns"] for row in role_rows["identity"])
        aa_left_elapsed = sum(
            row["elapsed_ns"] for row in role_rows["identity_a"])
        aa_right_elapsed = sum(
            row["elapsed_ns"] for row in role_rows["identity_b"])
        candidate_log_cost = (
            sum(candidate_pair_logs) / 4.0 if timing_eligible else None)
        aa_log_cost = (
            sum(aa_pair_logs) / 4.0 if timing_eligible else None)
        if timing_eligible:
            candidate_log_costs.append(float(candidate_log_cost))
            aa_log_costs.append(float(aa_log_cost))
        replicate_receipts.append({
            "replicate": measured_index,
            "construction_seed": seen_trace[absolute][0],
            "loss_seed": seen_trace[absolute][1],
            "trace_sha256": seen_trace[absolute][2],
            "common_overhead": seen_trace[absolute][3],
            "candidate_overhead": candidate_oh,
            "identity_overhead": identity_oh,
            "candidate_recovery_result": candidate_result,
            "identity_recovery_result": identity_result,
            "candidate_recovery_ok": candidate_ok,
            "identity_recovery_ok": identity_ok,
            "timing_eligible": timing_eligible,
            "candidate_elapsed_ns": candidate_elapsed,
            "identity_elapsed_ns": identity_elapsed,
            "identity_a_elapsed_ns": aa_left_elapsed,
            "identity_b_elapsed_ns": aa_right_elapsed,
            "candidate_pair_log_costs": list(candidate_pair_logs),
            "aa_pair_log_costs": list(aa_pair_logs),
            "candidate_log_cost": (
                float(candidate_log_cost) if timing_eligible else None),
            "aa_log_cost": (
                float(aa_log_cost) if timing_eligible else None),
        })

    timing_ci_available = timing_eligible_replicates >= 4
    if timing_ci_available:
        candidate_mean, candidate_low, candidate_high = _paired_t_interval(
            candidate_log_costs)
        aa_mean, aa_low, aa_high = _paired_t_interval(aa_log_costs)
        aa_floor = max(abs(aa_low), abs(aa_high))
    else:
        candidate_mean = candidate_low = candidate_high = 0.0
        aa_mean = aa_low = aa_high = aa_floor = 0.0
    required_lower_log, unused_required_upper_log = \
        _paired_practical_log_bounds(required_margin)
    required_log = (
        -required_lower_log if required_lower_log < 0.0 else 0.0)
    effective_margin = max(required_log, aa_floor)
    candidate_summary = _summarize_paired_arm(
        rows, replicates, block_count, "candidate")
    identity_summary = _summarize_paired_arm(
        rows, replicates, block_count, "identity")
    valid = (
        not candidate_is_exact_identity and
        timing_ci_available and
        aa_low <= 0.0 <= aa_high and
        recovery_regressions == 0 and
        candidate_high < -effective_margin
    )
    return PairedMeasurement(
        manifest=manifest,
        semantic=semantic,
        context={
            "schema": context["schema"],
            "bound": dict(context["bound"]),
            "thermal": dict(context["thermal"]),
        },
        candidate=candidate_summary,
        identity=identity_summary,
        replicate_receipts=tuple(replicate_receipts),
        candidate_log_cost_mean=float(candidate_mean),
        candidate_log_cost_ci_low=float(candidate_low),
        candidate_log_cost_ci_high=float(candidate_high),
        aa_log_cost_mean=float(aa_mean),
        aa_log_cost_ci_low=float(aa_low),
        aa_log_cost_ci_high=float(aa_high),
        aa_floor_log=float(aa_floor),
        required_margin=float(required_margin),
        effective_margin_log=float(effective_margin),
        timing_eligible_replicates=timing_eligible_replicates,
        timing_ci_available=timing_ci_available,
        candidate_recovery_result_counts=candidate_result_counts,
        identity_recovery_result_counts=identity_result_counts,
        candidate_weak_seed_failures=(
            candidate_result_counts["bad_dense_seed"] +
            candidate_result_counts["bad_peel_seed"]),
        identity_weak_seed_failures=(
            identity_result_counts["bad_dense_seed"] +
            identity_result_counts["bad_peel_seed"]),
        recovery_regressions=recovery_regressions,
        recovery_improvements=recovery_improvements,
        both_failures=both_failures,
        overhead_regressions=overhead_regressions,
        valid_for_promotion=valid,
        stream_sha256=stream_sha256,
        finished_monotonic_ns=finished_ns,
        final_context_sha256=_canonical_json_sha256(context),
        native_stdout=stdout,
    )


def make_search_receipt(
        metrics, *, mode, goodput, trials, block_bytes, search_kind,
        construction_seed_base, loss_seed_base, seed_domain, coordinates,
        peel_pmf, shipped_control, shipped_goodput, context=None):
    """Create a complete selected-arm receipt for a search result."""
    receipt = metrics.as_dict()
    receipt.update({
        "mode": mode,
        "goodput": goodput,
        "trials": trials,
        "block_bytes": block_bytes,
        "search_kind": search_kind,
        "construction_seed_base": construction_seed_base,
        "loss_seed_base": loss_seed_base,
        "seed_domain": seed_domain,
        "coordinates": dict(coordinates),
        "peel_pmf": list(peel_pmf),
        "peel_pmf_sha256": pmf_sha256(peel_pmf),
        "shipped_control": shipped_control.as_dict(),
        "shipped_goodput": shipped_goodput,
        "context": context or {},
    })
    return receipt


def paired_arm_goodput(arm, block_count):
    """Compute the recovery-adjusted solve rate for a paired arm summary."""
    if isinstance(arm, PairedArmSummary):
        return arm.goodput(block_count)
    if (not isinstance(arm, dict) or
            type(arm.get("solve_mbps")) is not float or
            type(arm.get("oh_mean")) is not float):
        return 0.0
    return (
        arm["solve_mbps"] * block_count /
        (block_count + arm["oh_mean"])
    )


def make_paired_search_receipt(
        measurement, *, mode, block_count, block_bytes, search_kind,
        construction_seed_base, loss_seed_base, seed_domain, coordinates,
        peel_pmf, evaluated_coordinates, evaluated_pmf,
        stock_control_measurement=None, context=None):
    """Create an authoritative v4 search receipt from one native paired run."""
    if not isinstance(measurement, PairedMeasurement):
        raise ValueError("paired search receipt requires PairedMeasurement")
    if mode not in ("trained", "scale-only", "shipped"):
        raise ValueError("invalid paired search mode")
    if mode == "shipped":
        if stock_control_measurement is not None:
            raise ValueError(
                "shipped search must use its selected paired measurement "
                "as stock control")
        stock_control_source = STOCK_CONTROL_SELECTED
        stock_control_receipt = None
    else:
        if not isinstance(stock_control_measurement, PairedMeasurement):
            raise ValueError(
                "trained search requires an embedded stock-control "
                "measurement")
        stock_control_source = STOCK_CONTROL_EMBEDDED
        stock_control_receipt = stock_control_measurement.as_dict()
    selected = (
        measurement.identity if mode == "shipped" else measurement.candidate)
    return {
        "protocol": NATIVE_PAIRED_PROTOCOL,
        "mode": mode,
        "selected_arm": "identity" if mode == "shipped" else "candidate",
        "goodput": selected.goodput(block_count),
        "shipped_goodput": measurement.identity.goodput(block_count),
        "trials": measurement.candidate.trials,
        "block_bytes": block_bytes,
        "search_kind": search_kind,
        "construction_seed_base": construction_seed_base,
        "loss_seed_base": loss_seed_base,
        "seed_domain": seed_domain,
        "coordinates": dict(coordinates),
        "peel_pmf": list(peel_pmf),
        "peel_pmf_sha256": pmf_sha256(peel_pmf),
        "evaluated_coordinates": dict(evaluated_coordinates),
        "evaluated_pmf": list(evaluated_pmf),
        "evaluated_pmf_sha256": pmf_sha256(evaluated_pmf),
        "paired_measurement": measurement.as_dict(),
        "stock_control_source": stock_control_source,
        "stock_control_measurement": stock_control_receipt,
        "context": context or {},
    }


def paired_selected_metrics(measurement, mode):
    """Return canonical top-level metrics for the selected paired arm."""
    arm = measurement.identity if mode == "shipped" else measurement.candidate
    return arm.as_dict()


def require_paired_stock_control(measurement, label="paired stock control"):
    """Require powered stock timing with no practically material bias."""
    if (not isinstance(measurement, PairedMeasurement) or
            not measurement.timing_ci_available):
        raise MeasurementError(f"{label} lacks powered timing evidence")

    # A nominal 95% interval excludes zero roughly one run in twenty under a
    # true null, so requiring both independent stock panels to contain zero
    # makes an all-K validation almost certain to abort spuriously.  Instead,
    # cross-calibrate the two panels: either panel is objectionably biased only
    # when its complete interval lies beyond both the predeclared practical
    # margin and the other panel's CI half-width.  Borrow uncertainty only,
    # never the other panel's location: including its absolute endpoints would
    # let two equally biased controls legitimize each other at any magnitude.
    # Equality is accepted because promotion itself requires a strict margin
    # improvement.
    required_lower_log, required_upper_log = \
        _paired_practical_log_bounds(measurement.required_margin)
    candidate_half_width = (
        measurement.candidate_log_cost_ci_high -
        measurement.candidate_log_cost_ci_low
    ) / 2.0
    aa_half_width = (
        measurement.aa_log_cost_ci_high -
        measurement.aa_log_cost_ci_low
    ) / 2.0
    candidate_lower_bound = min(required_lower_log, -aa_half_width)
    candidate_upper_bound = max(required_upper_log, aa_half_width)
    aa_lower_bound = min(required_lower_log, -candidate_half_width)
    aa_upper_bound = max(required_upper_log, candidate_half_width)
    candidate_biased = (
        measurement.candidate_log_cost_ci_high < candidate_lower_bound or
        measurement.candidate_log_cost_ci_low > candidate_upper_bound
    )
    aa_biased = (
        measurement.aa_log_cost_ci_high < aa_lower_bound or
        measurement.aa_log_cost_ci_low > aa_upper_bound
    )
    if candidate_biased or aa_biased:
        raise MeasurementError(
            f"{label} has a biased stock-vs-stock confidence interval")


def isolated_codec_env(overrides=None):
    """Return an environment with every inherited WH2 test hook removed."""
    env = {
        key: value for key, value in os.environ.items()
        if not key.startswith("WIREHAIR_V2_")
    }
    if overrides:
        env.update(overrides)
    return env


def _run_checked(command, env):
    try:
        result = subprocess.run(
            command, capture_output=True, text=True, env=env)
    except OSError as error:
        raise MeasurementError(
            f"could not execute benchmark {command[0]!r}: {error}")
    if result.returncode != 0:
        detail = result.stderr.strip() or result.stdout.strip() or "no output"
        raise MeasurementError(
            f"benchmark exited {result.returncode}: {' '.join(command)}: "
            f"{detail}")
    return result.stdout


def _sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as source:
        while True:
            chunk = source.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def file_sha256(path):
    """Return a file content identity as a lowercase SHA-256 digest."""
    try:
        return _sha256_file(path)
    except OSError as error:
        raise MeasurementError(f"could not hash {path!r}: {error}")


def _benchmark_identity(bench):
    """Hash a stable regular-file snapshot, detecting concurrent replacement."""
    path = os.path.realpath(bench)
    try:
        with open(path, "rb") as executable:
            before = os.fstat(executable.fileno())
            if not stat.S_ISREG(before.st_mode):
                raise OSError("benchmark is not a regular file")
            digest = hashlib.sha256()
            while True:
                chunk = executable.read(1024 * 1024)
                if not chunk:
                    break
                digest.update(chunk)
            after = os.fstat(executable.fileno())
        if not stat.S_ISREG(after.st_mode):
            raise OSError("benchmark is not a regular file")
    except OSError as error:
        raise MeasurementError(
            f"could not identify benchmark {bench!r}: {error}")
    stable_fields = (
        "st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns",
        "st_ctime_ns")
    if any(getattr(before, name) != getattr(after, name)
           for name in stable_fields):
        raise MeasurementError(
            f"benchmark changed while it was being identified: {path!r}")
    return {
        "path": path,
        "sha256": digest.hexdigest(),
        "size": before.st_size,
    }


def _source_provenance_files():
    """Enumerate the exact source paths incorporated into artifact identity."""
    try:
        source_files = subprocess.run(
            [
                "git", "-C", REPO_ROOT, "ls-files", "-z",
                "--cached", "--others", "--exclude-standard",
            ],
            capture_output=True)
        ignored_gitignores = subprocess.run(
            [
                "git", "-C", REPO_ROOT, "ls-files", "-z",
                "--others", "--ignored", "--exclude-standard", "--",
                ":(top).gitignore", ":(top,glob)**/.gitignore",
            ],
            capture_output=True)
    except OSError as error:
        raise MeasurementError(f"could not enumerate the source tree: {error}")
    if (source_files.returncode != 0 or
            ignored_gitignores.returncode != 0):
        raise MeasurementError("could not enumerate the source tree")

    relative_files = {
        raw.decode("utf-8", "surrogateescape")
        for raw in source_files.stdout.split(b"\0") if raw}
    # An untracked .gitignore can ignore itself.  The main source query still
    # consumes that file but omits it from its result, so enumerate ignored
    # ignore files separately and incorporate them into the same identity.
    relative_files.update(
        raw.decode("utf-8", "surrogateescape")
        for raw in ignored_gitignores.stdout.split(b"\0") if raw)
    selected = [
        (relative, os.path.join(REPO_ROOT, relative))
        for relative in sorted(relative_files)
        if (relative.endswith(_SOURCE_EXTENSIONS) or
            os.path.basename(relative) in ("CMakeLists.txt", ".gitignore"))
    ]
    for label, path in _git_exclude_provenance_paths():
        if os.path.isfile(path):
            selected.append((label, path))
    selected = tuple(sorted(selected))
    if not selected:
        raise MeasurementError("source identity contains no source files")
    return selected


def _git_config_introspection_environment():
    """Match ordinary Git commands when invoking the special `git config`."""
    environment = os.environ.copy()
    # The legacy variable is interpreted as `git config --file` by that
    # porcelain, but it does not redirect configuration for `git ls-files`.
    environment.pop("GIT_CONFIG", None)
    return environment


def _effective_global_git_exclude_path():
    """Return the exact global exclude path used by --exclude-standard."""
    try:
        configured = subprocess.run(
            [
                "git", "-C", REPO_ROOT, "config", "--null", "--path", "--get",
                "core.excludesFile",
            ],
            capture_output=True, env=_git_config_introspection_environment())
    except OSError as error:
        raise MeasurementError(
            f"could not locate Git's global excludes file: {error}")
    if configured.returncode == 0:
        if (configured.stderr or
                not configured.stdout.endswith(b"\0") or
                configured.stdout.count(b"\0") != 1 or
                not configured.stdout[:-1]):
            raise MeasurementError("Git's global excludes path is invalid")
        path = configured.stdout[:-1].decode(
            "utf-8", "surrogateescape")
        if not os.path.isabs(path):
            path = os.path.join(REPO_ROOT, path)
        return os.path.abspath(path)
    if (configured.returncode != 1 or configured.stdout or
            configured.stderr):
        raise MeasurementError("could not identify Git's global excludes file")

    xdg_config = os.environ.get("XDG_CONFIG_HOME")
    if not xdg_config:
        xdg_config = os.path.join(os.path.expanduser("~"), ".config")
    if not os.path.isabs(xdg_config):
        xdg_config = os.path.join(REPO_ROOT, xdg_config)
    return os.path.abspath(os.path.join(xdg_config, "git", "ignore"))


def _git_exclude_provenance_paths():
    """Locate non-.gitignore inputs consumed by --exclude-standard."""
    try:
        info = subprocess.run(
            [
                "git", "-C", REPO_ROOT, "rev-parse",
                "--path-format=absolute", "--git-path", "info/exclude",
            ],
            capture_output=True, text=True)
    except OSError as error:
        raise MeasurementError(
            f"could not locate Git's repository excludes file: {error}")
    lines = info.stdout.splitlines()
    if (info.returncode != 0 or len(lines) != 1 or
            not os.path.isabs(lines[0]) or "\0" in lines[0]):
        raise MeasurementError(
            "could not identify Git's repository excludes file")
    return (
        ("@git/info-exclude", lines[0]),
        ("@git/global-excludes", _effective_global_git_exclude_path()),
    )


def _git_configuration_provenance_paths():
    """Locate file-backed Git configuration that selects exclude inputs."""
    paths = []

    def add_path(path, base=REPO_ROOT, expand_user=False):
        if expand_user:
            path = os.path.expanduser(path)
        if "%(" in path:
            raise MeasurementError(
                "Git configuration path interpolation is unsupported")
        path = os.path.abspath(
            path if os.path.isabs(path) else os.path.join(base, path))
        if path not in paths:
            paths.append(path)

    try:
        repository_configs = subprocess.run(
            [
                "git", "-C", REPO_ROOT, "rev-parse",
                "--path-format=absolute", "--git-path", "config",
                "--git-path", "config.worktree",
            ],
            capture_output=True, text=True)
    except OSError as error:
        raise MeasurementError(
            f"could not locate repository Git configuration: {error}")
    repository_config_paths = repository_configs.stdout.splitlines()
    if (repository_configs.returncode != 0 or
            len(repository_config_paths) != 2 or
            repository_configs.stderr or
            any(not os.path.isabs(path) or "\0" in path
                for path in repository_config_paths)):
        raise MeasurementError(
            "could not identify repository Git configuration")
    for path in repository_config_paths:
        add_path(path)

    for variable in ("GIT_CONFIG_SYSTEM", "GIT_CONFIG_GLOBAL"):
        # `git var` separates multiple default paths with newlines, which is
        # ambiguous when an explicit environment pathname itself contains a
        # newline.  Each environment override is one exact pathname, so retain
        # it directly rather than reparsing Git's line-oriented rendering.
        if variable in os.environ:
            explicit_path = os.environ[variable]
            if explicit_path:
                add_path(explicit_path)
            continue
        if variable == "GIT_CONFIG_GLOBAL":
            # Git's global defaults are derived from environment pathnames,
            # which may legally contain newlines.  Preserve those candidates
            # directly instead of splitting `git var`'s line-oriented output.
            home = os.environ.get("HOME")
            if not home:
                raise MeasurementError(
                    "could not identify GIT_CONFIG_GLOBAL without HOME")
            xdg_config = os.environ.get("XDG_CONFIG_HOME")
            if xdg_config:
                add_path(os.path.join(xdg_config, "git", "config"))
            else:
                add_path(os.path.join(home, ".config", "git", "config"))
            add_path(os.path.join(home, ".gitconfig"))
            continue
        try:
            defaults = subprocess.run(
                ["git", "-C", REPO_ROOT, "var", variable],
                capture_output=True, text=True,
                env=_git_config_introspection_environment())
        except OSError as error:
            raise MeasurementError(
                f"could not locate {variable}: {error}")
        lines = defaults.stdout.splitlines()
        if defaults.returncode == 0:
            if not lines or any(not line or "\0" in line for line in lines):
                raise MeasurementError(f"{variable} paths are invalid")
            for path in lines:
                add_path(path)
        elif (defaults.returncode != 1 or defaults.stdout or
              defaults.stderr):
            raise MeasurementError(f"could not identify {variable}")

    try:
        configuration = subprocess.run(
            [
                "git", "-C", REPO_ROOT, "config", "--null",
                "--show-origin", "--list",
            ],
            capture_output=True, env=_git_config_introspection_environment())
    except OSError as error:
        raise MeasurementError(f"could not locate Git configuration: {error}")
    fields = configuration.stdout.split(b"\0")
    if fields and fields[-1] == b"":
        fields.pop()
    if configuration.returncode != 0 or len(fields) % 2 != 0:
        raise MeasurementError("could not identify Git configuration inputs")
    for raw_origin, raw_entry in zip(fields[0::2], fields[1::2]):
        origin = raw_origin.decode("utf-8", "surrogateescape")
        origin_path = None
        if origin.startswith("file:"):
            origin_path = origin[len("file:"):]
            if not origin_path or "\0" in origin_path:
                raise MeasurementError("Git configuration origin is invalid")
            add_path(origin_path)
            origin_path = os.path.abspath(
                origin_path if os.path.isabs(origin_path)
                else os.path.join(REPO_ROOT, origin_path))

        entry = raw_entry.decode("utf-8", "surrogateescape")
        name, separator, value = entry.partition("\n")
        normalized_name = name.lower()
        is_include_path = (
            normalized_name == "include.path" or
            (normalized_name.startswith("includeif.") and
             normalized_name.endswith(".path"))
        )
        if is_include_path:
            if not separator or not value or "\0" in value:
                raise MeasurementError(
                    "Git configuration include path is invalid")
            add_path(
                value,
                os.path.dirname(origin_path)
                if origin_path is not None else REPO_ROOT,
                expand_user=True)
    return tuple(paths)


def _git_provenance_paths():
    """Locate Git state read while deriving the source provenance identity."""
    try:
        locations = subprocess.run(
            [
                "git", "-C", REPO_ROOT, "rev-parse",
                "--path-format=absolute", "--git-dir", "--git-common-dir",
                "--git-path", "HEAD", "--git-path", "index",
                "--git-path", "packed-refs",
            ],
            capture_output=True, text=True)
    except OSError as error:
        raise MeasurementError(f"could not locate Git metadata: {error}")
    location_lines = locations.stdout.splitlines()
    if (locations.returncode != 0 or len(location_lines) != 5 or
            any(not os.path.isabs(path) or "\0" in path
                for path in location_lines)):
        raise MeasurementError("could not locate Git metadata")
    git_directory, common_directory, head_path, index_path, packed_refs = (
        location_lines)
    if not os.path.isdir(git_directory) or not os.path.isdir(common_directory):
        raise MeasurementError("Git metadata directories are unavailable")

    git_files = [
        ("the worktree Git locator", os.path.join(REPO_ROOT, ".git")),
        (
            "the Git common-directory locator",
            os.path.join(git_directory, "commondir"),
        ),
        ("Git HEAD", head_path),
        ("the Git index", index_path),
        ("Git packed refs", packed_refs),
    ]
    try:
        shared_index = subprocess.run(
            ["git", "-C", REPO_ROOT, "rev-parse", "--shared-index-path"],
            capture_output=True, text=True)
    except OSError as error:
        raise MeasurementError(f"could not locate Git shared index: {error}")
    shared_index_path = shared_index.stdout.rstrip("\n")
    if (shared_index.returncode != 0 or shared_index.stderr or
            "\n" in shared_index_path or "\0" in shared_index_path):
        raise MeasurementError("could not identify Git shared index")
    if shared_index_path:
        if not os.path.isabs(shared_index_path):
            shared_index_path = os.path.join(REPO_ROOT, shared_index_path)
        git_files.append(
            ("the Git shared index", os.path.abspath(shared_index_path)))
    for label, path in _git_exclude_provenance_paths():
        git_files.append((f"Git exclude input {label}", path))
    for path in _git_configuration_provenance_paths():
        git_files.append((f"Git configuration input {path}", path))
    symbolic_name = "HEAD"
    symbolic_names = {"HEAD"}
    for hop in range(256):
        try:
            symbolic = subprocess.run(
                [
                    "git", "-C", REPO_ROOT, "symbolic-ref", "-q",
                    "--no-recurse", symbolic_name,
                ],
                capture_output=True, text=True)
        except OSError as error:
            raise MeasurementError(
                f"could not locate Git symbolic ref state: {error}")
        if symbolic.returncode == 1:
            if symbolic.stdout or symbolic.stderr:
                raise MeasurementError(
                    "could not identify Git symbolic ref state")
            break
        if symbolic.returncode != 0 or symbolic.stderr:
            raise MeasurementError(
                "could not identify Git symbolic ref state")

        reference = symbolic.stdout.rstrip("\n")
        if (not reference or "\n" in reference or "\0" in reference or
                os.path.isabs(reference) or reference in symbolic_names):
            raise MeasurementError("Git symbolic ref chain is invalid")
        symbolic_names.add(reference)
        try:
            reference_path = subprocess.run(
                [
                    "git", "-C", REPO_ROOT, "rev-parse",
                    "--path-format=absolute", "--git-path", reference,
                ],
                capture_output=True, text=True)
        except OSError as error:
            raise MeasurementError(
                f"could not locate Git symbolic ref path: {error}")
        reference_lines = reference_path.stdout.splitlines()
        if (reference_path.returncode != 0 or reference_path.stderr or
                len(reference_lines) != 1 or
                not os.path.isabs(reference_lines[0]) or
                "\0" in reference_lines[0]):
            raise MeasurementError("could not locate Git symbolic ref path")
        label = (
            "the Git symbolic HEAD ref" if hop == 0 else
            f"Git symbolic ref target {reference}"
        )
        git_files.append((label, reference_lines[0]))
        symbolic_name = reference
    else:
        raise MeasurementError("Git symbolic ref chain is too deep")

    directories = []
    for label, directory in (
            ("the worktree Git directory", git_directory),
            ("the common Git directory", common_directory)):
        if all(os.path.abspath(directory) != os.path.abspath(existing)
               for _, existing in directories):
            directories.append((label, directory))
    return tuple(directories), tuple(git_files)


def _artifact_identity(bench, generator):
    """Return immutable executable and source-state identities."""
    benchmark = _benchmark_identity(bench)
    generator_path = (
        generator if os.path.isabs(generator)
        else os.path.join(REPO_ROOT, generator)
    )
    try:
        generator_sha256 = _sha256_file(generator_path)
    except OSError as error:
        raise MeasurementError(
            f"could not identify generator source: {error}")

    try:
        git = subprocess.run(
            ["git", "-C", REPO_ROOT, "rev-parse", "HEAD"],
            capture_output=True, text=True)
    except OSError as error:
        raise MeasurementError(f"could not identify source commit: {error}")
    commit = git.stdout.strip()
    if git.returncode != 0 or len(commit) not in (40, 64) or any(
            character not in "0123456789abcdef" for character in commit):
        raise MeasurementError("could not identify the source Git commit")

    source_digest = hashlib.sha256()
    source_count = 0
    for relative, path in _source_provenance_files():
        try:
            content_sha256 = _sha256_file(path)
        except OSError as error:
            raise MeasurementError(
                f"could not identify source file {relative}: {error}")
        source_digest.update(relative.encode("utf-8", "surrogateescape"))
        source_digest.update(b"\0")
        source_digest.update(content_sha256.encode("ascii"))
        source_digest.update(b"\0")
        source_count += 1
    return {
        "benchmark": benchmark,
        "source": {
            "git_commit": commit,
            "state_sha256": source_digest.hexdigest(),
            "file_count": source_count,
            "generator_sha256": generator_sha256,
        },
    }


def _python_runtime_identity():
    libc_name, libc_version = platform.libc_ver()
    return {
        "implementation": sys.implementation.name,
        "version": ".".join(str(value) for value in sys.version_info[:3]),
        "cache_tag": sys.implementation.cache_tag or "",
        "libc": f"{libc_name}-{libc_version}",
    }


def capture_artifact_identity(bench, generator):
    """Capture the executable and source state before the first measurement."""
    return _artifact_identity(bench, generator)


def verify_artifact_identity(identity, bench, generator):
    """Reject publication if executable or source state changed during a run."""
    current = _artifact_identity(bench, generator)
    if current != identity:
        raise MeasurementError(
            "benchmark, generator, or source state changed during measurement")


def verify_benchmark_identity(document, bench):
    """Refuse replay with a benchmark executable different from provenance."""
    try:
        expected = document["provenance"]["benchmark"]
        actual = _benchmark_identity(bench)
    except (KeyError, MeasurementError, TypeError) as error:
        raise MeasurementError(f"could not verify benchmark identity: {error}")
    if (actual["size"] != expected["size"] or
            actual["sha256"] != expected["sha256"]):
        raise MeasurementError(
            f"benchmark identity mismatch for {actual['path']!r}; table expects "
            f"sha256={expected['sha256']}")


def derive_seed(base_seed, *domain):
    """Derive a stable, domain-separated uint64 benchmark seed."""
    if (isinstance(base_seed, bool) or not isinstance(base_seed, int) or
            not 0 <= base_seed <= 0xffffffffffffffff):
        raise ValueError("base seed must be a uint64")
    digest = hashlib.blake2b(digest_size=8, person=b"wh2-peel-v2")
    digest.update(base_seed.to_bytes(8, "little"))
    for item in domain:
        encoded = str(item).encode("utf-8")
        digest.update(len(encoded).to_bytes(4, "little"))
        digest.update(encoded)
    return int.from_bytes(digest.digest(), "little")


_STOCK_CACHE = {}


def _production_staircase_max(block_count):
    """Largest S accepted by the production mixed precode parameter domain."""
    return (
        _PRODUCTION_TOTAL_SPAN_MAX - block_count -
        _PRODUCTION_DENSE_ROWS - _PRODUCTION_MIXED_HEAVY_ROWS
    )


def _cpp_trunc_div(numerator, denominator):
    """Divide signed integers with C++'s truncation toward zero."""
    quotient = abs(numerator) // denominator
    return -quotient if numerator < 0 else quotient


def _legacy_dense_count(block_count):
    """Mirror wirehair::GetDenseCount() exactly for K in [2, 64000]."""
    if (type(block_count) is not int or
            not 2 <= block_count <= 64000):
        raise ValueError("legacy dense-count K must be in [2,64000]")

    if block_count < len(_LEGACY_TINY_DENSE_COUNTS):
        return _LEGACY_TINY_DENSE_COUNTS[block_count]
    if block_count < 2048:
        if block_count <= 500:
            low_k, high_k, low_count, high_count = 64, 500, 26, 35
        elif block_count <= 1000:
            low_k, high_k, low_count, high_count = 500, 1000, 35, 48
        else:
            low_k, high_k, low_count, high_count = 1000, 2048, 48, 62
    else:
        # Preserve native's strict ``N > mid.N`` binary search, including its
        # duplicate K=54811 point.  Equality selects the preceding segment;
        # K=54812 selects the segment beginning at the second duplicate.
        low = 0
        high = len(_LEGACY_DENSE_POINTS) - 1
        while True:
            middle = (high + low) // 2
            if middle == low:
                break
            if block_count > _LEGACY_DENSE_POINTS[middle][0]:
                low = middle
            else:
                high = middle
        low_k, low_count = _LEGACY_DENSE_POINTS[low]
        high_k, high_count = _LEGACY_DENSE_POINTS[low + 1]

    numerator = (block_count - low_k) * (high_count - low_count)
    dense_count = (
        low_count +
        _cpp_trunc_div(numerator, high_k - low_k)
    )
    # Native rounds upward to the next D such that D mod 4 == 2.
    dense_count += (2 - dense_count % 4) % 4
    return dense_count


def _dispatch_staircase_count(block_count):
    """Mirror dispatch-v1 SmallBandStaircaseCount() exactly."""
    if (type(block_count) is not int or
            block_count < 2 or block_count > 64000):
        return 0
    inherited = _legacy_dense_count(block_count)
    if block_count > 100:
        return inherited
    staircase = 2
    while 16 * (staircase + 1) ** 2 <= 25 * block_count:
        staircase += 1
    return min(staircase, inherited)


def _dispatch_source_hits(block_count):
    """Mirror dispatch-v1 CertifiedSourceHits()."""
    return 3 if block_count >= 10000 else 2


def stock_profile(bench, block_count, *, target_profile):
    """Return the exact native PMF and dispatch contract for one K."""
    if (isinstance(block_count, bool) or
            not isinstance(block_count, int) or
            not 2 <= block_count <= 64000):
        raise ValueError("native PMF K must be in [2,64000]")
    if target_profile != TARGET_PROFILE:
        raise ValueError(
            f"target_profile must be exactly {TARGET_PROFILE!r}")
    benchmark = _benchmark_identity(bench)
    key = (
        benchmark["path"], benchmark["sha256"], target_profile,
        TARGET_CONTRACT["contract_id"], TARGET_CONTRACT["contract_sha256"],
        block_count,
    )
    cached = _STOCK_CACHE.get(key)
    if cached is not None:
        return cached

    stdout = _run_checked(
        [
            bench, "peelpmf", "--N", str(block_count),
            "--target-profile", target_profile,
        ],
        isolated_codec_env())
    if _benchmark_identity(bench) != benchmark:
        raise MeasurementError(
            "benchmark changed while reading the native peel profile")
    metadata = None
    csv_header = None
    probabilities = {}
    for line in stdout.splitlines():
        if not line.strip():
            continue
        if line.startswith("# peelpmf,"):
            if metadata is not None or csv_header is not None or probabilities:
                raise MeasurementError(
                    f"misordered or duplicate peelpmf metadata "
                    f"for K={block_count}")
            fields = line.split(",")
            try:
                values = {}
                for field in fields[1:]:
                    name, value = field.split("=", 1)
                    if name in values:
                        raise ValueError("duplicate metadata field")
                    values[name] = value
                expected_fields = {
                    "N", "target_profile", "contract_id",
                    "contract_sha256", "precode_contract",
                    "packet_contract", "architecture", "degrees",
                    "staircase", "dense_rows", "heavy_rows", "source_hits",
                    "completion", "gf256_rows", "gf16_rows", "period",
                    "geometry", "residue_schedule", "residue_skew",
                    "mix_count", "target_mean", "pmf_sha256",
                    "pmf_encoding",
                }
                if set(values) != expected_fields:
                    raise ValueError("unexpected metadata fields")
                metadata = {
                    "N": int(values["N"]),
                    "target_profile": values["target_profile"],
                    "contract_id": values["contract_id"],
                    "contract_sha256": values["contract_sha256"],
                    "precode_contract": int(values["precode_contract"]),
                    "packet_contract": int(values["packet_contract"]),
                    "architecture": values["architecture"],
                    "degrees": int(values["degrees"]),
                    "staircase": int(values["staircase"]),
                    "dense_rows": int(values["dense_rows"]),
                    "heavy_rows": int(values["heavy_rows"]),
                    "source_hits": int(values["source_hits"]),
                    "completion": values["completion"],
                    "gf256_rows": int(values["gf256_rows"]),
                    "gf16_rows": int(values["gf16_rows"]),
                    "period": int(values["period"]),
                    "geometry": values["geometry"],
                    "residue_schedule": values["residue_schedule"],
                    "residue_skew": int(values["residue_skew"]),
                    "mix_count": int(values["mix_count"]),
                    "target_mean": float(values["target_mean"]),
                    "pmf_sha256": values["pmf_sha256"],
                    "pmf_encoding": values["pmf_encoding"],
                }
            except (KeyError, TypeError, ValueError):
                raise MeasurementError(
                    f"malformed peelpmf metadata for K={block_count}: {line}")
            continue
        fields = [field.strip() for field in line.split(",")]
        if fields == ["degree", "probability"]:
            if metadata is None or csv_header is not None or probabilities:
                raise MeasurementError(
                    f"misordered or duplicate peelpmf CSV header "
                    f"for K={block_count}")
            csv_header = fields
            continue
        if csv_header is None or len(fields) != len(csv_header):
            raise MeasurementError(
                f"unexpected peelpmf output for K={block_count}: {line}")
        try:
            degree = int(fields[csv_header.index("degree")])
            probability = float(fields[csv_header.index("probability")])
        except (ValueError, IndexError):
            raise MeasurementError(
                f"malformed peelpmf row for K={block_count}: {line}")
        if degree in probabilities or degree < 1 or degree > 64:
            raise MeasurementError(
                f"invalid peelpmf degree for K={block_count}: {line}")
        probabilities[degree] = probability

    if metadata is None:
        raise MeasurementError(
            f"peelpmf returned no metadata for K={block_count}")
    expected_staircase = _dispatch_staircase_count(block_count)
    expected_source_hits = _dispatch_source_hits(block_count)
    expected_target_mean = (
        block_count * min(expected_source_hits, expected_staircase) /
        expected_staircase
    )
    contract_matches = all(
        metadata[name] == value
        for name, value in TARGET_CONTRACT.items()
        if name in metadata
    )
    if (metadata["N"] != block_count or
            metadata["target_profile"] != target_profile or
            not contract_matches or metadata["degrees"] != 64 or
            metadata["staircase"] != expected_staircase or
            metadata["staircase"] >
            _production_staircase_max(block_count) or
            metadata["source_hits"] != expected_source_hits or
            not math.isfinite(metadata["target_mean"]) or
            metadata["target_mean"] <= 0.0 or
            metadata["target_mean"] != expected_target_mean or
            metadata["pmf_sha256"] != STOCK_PMF_DIGEST):
        raise MeasurementError(
            f"invalid peelpmf metadata for K={block_count}: {metadata}")
    if set(probabilities) != set(range(1, 65)):
        raise MeasurementError(
            f"peelpmf returned {len(probabilities)} of 64 degrees "
            f"for K={block_count}")
    pmf = tuple(probabilities[degree] for degree in range(1, 65))
    if (any(value == 0.0 and math.copysign(1.0, value) < 0.0
            for value in pmf) or
            any(not math.isfinite(value) or value < 0.0 for value in pmf) or
            abs(sum(pmf) - 1.0) > 1e-12):
        raise MeasurementError(
            f"peelpmf returned a non-probability distribution "
            f"for K={block_count}")
    profile = StockProfile(
        block_count=metadata["N"],
        target_profile=metadata["target_profile"],
        contract_id=metadata["contract_id"],
        contract_sha256=metadata["contract_sha256"],
        precode_contract=metadata["precode_contract"],
        packet_contract=metadata["packet_contract"],
        architecture=metadata["architecture"],
        staircase=metadata["staircase"],
        dense_rows=metadata["dense_rows"],
        heavy_rows=metadata["heavy_rows"],
        source_hits=metadata["source_hits"],
        completion=metadata["completion"],
        gf256_rows=metadata["gf256_rows"],
        gf16_rows=metadata["gf16_rows"],
        period=metadata["period"],
        geometry=metadata["geometry"],
        residue_schedule=metadata["residue_schedule"],
        residue_skew=metadata["residue_skew"],
        mix_count=metadata["mix_count"],
        target_mean=metadata["target_mean"],
        native_pmf_sha256=metadata["pmf_sha256"],
        pmf_encoding=metadata["pmf_encoding"],
        pmf=pmf,
    )
    _STOCK_CACHE[key] = profile
    return profile


def stock_pmf(bench, block_count, *, target_profile):
    """Return degrees 1..64 of the exact PMF recovered by the native codec."""
    return list(stock_profile(
        bench, block_count, target_profile=target_profile).pmf)


def pmf_sha256(pmf):
    """Hash the exact binary64 text law passed through the benchmark hook."""
    try:
        numeric = [float(probability) for probability in pmf]
    except (OverflowError, TypeError, ValueError):
        raise ValueError("invalid PMF")
    numeric = [0.0 if value == 0.0 else value for value in numeric]
    encoded = ",".join(
        f"{value:.17g}" for value in numeric).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def family(stock, p1, tilt, dmax, absorb):
    """Apply the native search family's coordinates to a shipped PMF."""
    try:
        coordinates = tuple(float(value)
                            for value in (p1, tilt, dmax, absorb))
    except (OverflowError, TypeError, ValueError):
        return None
    if any(not math.isfinite(value) for value in coordinates):
        return None
    p1, tilt, dmax, absorb = coordinates
    dmax = max(2, min(int(dmax), len(stock)))
    weights = [0.0] * dmax
    try:
        weights[0] = (p1 / 100.0) * stock[0]
        for degree in range(2, dmax):
            weights[degree - 1] = (
                stock[degree - 1] * (degree ** (-tilt / 100.0)))
        absorption = min(max(absorb, 0), 100) / 100.0
        weights[dmax - 1] = (
            absorption * sum(stock[dmax - 1:]) +
            (1.0 - absorption) * stock[dmax - 1]
        ) * (dmax ** (-tilt / 100.0))
    except (OverflowError, TypeError, ValueError):
        return None
    total = sum(value for value in weights if value > 0.0)
    if not math.isfinite(total) or not total > 0.0:
        return None
    return [max(value, 0.0) / total for value in weights]


def compare_probe(
        bench, block_count, trials, block_bytes, peel_weights=None,
        degree_scale=None, *, native_profile, target_profile, seed_policy,
        loss, schedule, construction_seed, loss_seed):
    """Measure one exact raw dispatch-v1 arm under split deterministic seeds."""
    for name, value in (
            ("construction_seed", construction_seed),
            ("loss_seed", loss_seed)):
        if (isinstance(value, bool) or not isinstance(value, int) or
                not 0 <= value <= 0xffffffffffffffff):
            raise ValueError(f"{name} must be a uint64")
    if target_profile != TARGET_PROFILE:
        raise ValueError(
            f"target_profile must be exactly {TARGET_PROFILE!r}")
    if seed_policy != TARGET_SEED_POLICY:
        raise ValueError(f"seed_policy must be exactly {TARGET_SEED_POLICY!r}")
    if (not isinstance(schedule, str) or
            schedule not in COMPARE_TARGET_SCHEDULES):
        raise ValueError("invalid target schedule")
    if not valid_loss_rate(loss):
        raise ValueError("loss must be finite and in [0,0.99]")
    loss = float(loss)
    if (isinstance(block_count, bool) or
            not isinstance(block_count, int) or
            not 2 <= block_count <= 64000 or
            isinstance(trials, bool) or
            not isinstance(trials, int) or
            not 1 <= trials <= _COMPARE_TRIALS_MAX or
            isinstance(block_bytes, bool) or
            not isinstance(block_bytes, int) or
            not 1 <= block_bytes <= _COMPARE_BLOCK_BYTES_MAX or
            block_bytes % 2 != 0):
        raise ValueError("invalid compare K, trial count, or block size")
    if (not isinstance(native_profile, StockProfile) or
            native_profile.block_count != block_count or
            native_profile.target_profile != target_profile):
        raise ValueError("native_profile does not match compare K and target")
    _validate_native_receipt(
        native_profile.as_dict(), block_count, "compare native_profile")

    overrides = {}
    peel_spec = None
    if peel_weights is not None:
        try:
            numeric_weights = []
            for value in peel_weights:
                if isinstance(value, bool):
                    raise ValueError("Boolean peel weight")
                numeric_weights.append(float(value))
        except (OverflowError, TypeError, ValueError):
            raise ValueError("invalid peel weight vector")
        total_weight = sum(numeric_weights)
        if (not 2 <= len(numeric_weights) <= 64 or
                any(not math.isfinite(value) or value < 0.0
                    for value in numeric_weights) or
                not math.isfinite(total_weight) or not total_weight > 0.0):
            raise ValueError("invalid peel weight vector")
        peel_spec = ",".join(
            f"{value:.17g}" for value in numeric_weights)
        overrides["WIREHAIR_V2_PEEL_DEGREES"] = peel_spec
    scale_spec = "unset"
    if degree_scale is not None:
        if isinstance(degree_scale, bool):
            raise ValueError("invalid staircase degree scale")
        try:
            numeric_scale = float(degree_scale)
        except (OverflowError, TypeError, ValueError):
            raise ValueError("invalid staircase degree scale")
        if (not math.isfinite(numeric_scale) or
                not 0.0 <= numeric_scale <= 64000.0):
            raise ValueError("invalid staircase degree scale")
        scale_spec = _canonical_staircase_scale_spec(numeric_scale)
        overrides["WIREHAIR_V2_STAIRCASE_DEGREE_SCALE"] = scale_spec
    expected_pmf_digest = hashlib.sha256(
        (peel_spec if peel_spec is not None else "stock").encode(
            "ascii")).hexdigest()

    command = [
        bench, "compare",
        "--nlo", str(block_count),
        "--nhi", str(block_count),
        "--bb-list", str(block_bytes),
        "--trials", str(trials),
        "--loss", f"{loss:.17g}",
        "--loss-seed", str(loss_seed),
        "--schedule", schedule,
        # A fixed cap can silently suppress the requested row when K * bb is
        # larger than the benchmark's default 128 MiB cap.
        "--max-message-mib", "0",
        "--precode",
        "--target-profile", target_profile,
        "--seed-policy", seed_policy,
        "--construction-seed", str(construction_seed),
    ]
    stdout = _run_checked(command, isolated_codec_env(overrides))
    target_receipt = None
    compare_metadata = None
    header = None
    rows = []
    row_codecs = []
    for line in stdout.splitlines():
        if not line.strip():
            continue
        if line.startswith("# wh2_target,"):
            if (target_receipt is not None or compare_metadata is not None or
                    header is not None or row_codecs):
                raise MeasurementError(
                    f"misordered or duplicate wh2_target receipt "
                    f"for K={block_count}")
            target_receipt = {}
            for field in line[len("# wh2_target,"):].split(","):
                if "=" not in field:
                    raise MeasurementError(
                        f"malformed wh2_target receipt for K={block_count}: "
                        f"{line}")
                name, value = field.split("=", 1)
                if not name or name in target_receipt:
                    raise MeasurementError(
                        f"duplicate wh2_target field {name!r} "
                        f"for K={block_count}")
                target_receipt[name] = value
            continue
        if line.startswith("# compare:"):
            if (target_receipt is None or compare_metadata is not None or
                    header is not None or row_codecs):
                raise MeasurementError(
                    f"misordered or duplicate compare metadata "
                    f"for K={block_count}")
            compare_metadata = {}
            for field in line[len("# compare:"):].split():
                if "=" not in field:
                    raise MeasurementError(
                        f"malformed compare metadata for K={block_count}: "
                        f"{line}")
                name, value = field.split("=", 1)
                if name in compare_metadata:
                    raise MeasurementError(
                        f"duplicate compare metadata field {name!r} "
                        f"for K={block_count}")
                compare_metadata[name] = value
            continue
        if line.startswith("codec"):
            parsed_header = tuple(line.split())
            if (compare_metadata is None or header is not None or row_codecs or
                    parsed_header != COMPARE_COLUMNS):
                raise MeasurementError(
                    f"unexpected, misordered, or duplicate compare header "
                    f"for K={block_count}")
            header = list(parsed_header)
            continue
        fields = line.split()
        codec = fields[0]
        if codec not in ("baseline", "v2", "v2_target"):
            raise MeasurementError(
                f"unexpected compare output for K={block_count}: {line}")
        if (header is None or len(fields) != len(header) or
                codec in row_codecs):
            raise MeasurementError(
                f"malformed or duplicate compare row for K={block_count}: "
                f"{line}")
        row_codecs.append(codec)
        try:
            row_block_bytes = int(fields[header.index("bb")])
            row_trials = int(fields[header.index("trials")])
            row_fail = int(fields[header.index("fail")])
            row_numeric = {
                name: float(fields[header.index(name)])
                for name in (
                    "N_mean", "OH_mean", "OH_sd", "OH50", "OH95",
                    "OH99", "OH_max", "create_MBps", "encode_MBps",
                    "decode_MBps", "recover_MBps")
            }
        except (ValueError, IndexError):
            raise MeasurementError(
                f"malformed compare row for K={block_count}: {line}")
        row_n_mean = row_numeric["N_mean"]
        if (row_block_bytes != block_bytes or row_trials != trials or
                not 0 <= row_fail <= trials or
                not math.isfinite(row_n_mean) or
                row_n_mean != float(block_count) or
                any(not math.isfinite(value) or value < 0.0
                    for value in row_numeric.values()) or
                row_numeric["OH_mean"] > row_numeric["OH_max"] or
                not (row_numeric["OH50"] <= row_numeric["OH95"] <=
                     row_numeric["OH99"] <= row_numeric["OH_max"])):
            raise MeasurementError(
                f"semantically wrong compare row for K={block_count}: {line}")
        if codec == "v2_target":
            try:
                values = {
                    name: fields[header.index(name)]
                    for name in (
                        "bb", "trials", "N_mean", "fail", "OH_mean", "OH_sd",
                        "OH50", "OH95", "OH99", "OH_max", "decode_MBps")
                }
                metrics = RecoveryMetrics(
                    construction_seed=construction_seed,
                    loss_seed=loss_seed,
                    target_receipt=dict(target_receipt or {}),
                    fail=int(values["fail"]),
                    oh_mean=float(values["OH_mean"]),
                    oh_sd=float(values["OH_sd"]),
                    oh50=float(values["OH50"]),
                    oh95=float(values["OH95"]),
                    oh99=float(values["OH99"]),
                    oh_max=float(values["OH_max"]),
                    decode_mbps=float(values["decode_MBps"]),
                )
            except (KeyError, ValueError, IndexError):
                raise MeasurementError(
                    f"malformed compare row for K={block_count}: {line}")
            numeric = (
                metrics.oh_mean, metrics.oh_sd, metrics.oh50, metrics.oh95,
                metrics.oh99, metrics.oh_max, metrics.decode_mbps)
            if (not 0 <= metrics.fail <= trials or
                    any(not math.isfinite(value) or value < 0.0
                        for value in numeric) or
                    metrics.oh_mean > metrics.oh_max or
                    not (metrics.oh50 <= metrics.oh95 <= metrics.oh99 <=
                         metrics.oh_max)):
                raise MeasurementError(
                    f"invalid compare row for K={block_count}: {line}")
            rows.append(metrics)
    if target_receipt is None:
        raise MeasurementError(
            f"compare returned no wh2_target receipt for K={block_count}")
    if compare_metadata is None:
        raise MeasurementError(
            f"compare returned no metadata for K={block_count}")

    expected_target_receipt = {
        "profile": target_profile,
        "contract_id": TARGET_CONTRACT["contract_id"],
        "contract_sha256": TARGET_CONTRACT["contract_sha256"],
        "precode_contract": str(TARGET_CONTRACT["precode_contract"]),
        "packet_contract": str(TARGET_CONTRACT["packet_contract"]),
        "architecture": TARGET_CONTRACT["architecture"],
        **TARGET_MEASUREMENT_POLICY,
        "N": str(block_count),
        "bb": str(block_bytes),
        "staircase": str(native_profile.staircase),
        "dense_rows": str(native_profile.dense_rows),
        "heavy_rows": str(native_profile.heavy_rows),
        "completion": native_profile.completion,
        "gf256_rows": str(native_profile.gf256_rows),
        "gf16_rows": str(native_profile.gf16_rows),
        "period": str(native_profile.period),
        "geometry": native_profile.geometry,
        "residue_schedule": native_profile.residue_schedule,
        "residue_skew": str(native_profile.residue_skew),
        "source_hits": str(native_profile.source_hits),
        "target_mean": _target_mean_spec(
            block_count, native_profile.staircase,
            native_profile.source_hits, scale_spec),
        "mix_count": str(native_profile.mix_count),
        "packet_peel_seed": str(
            ((construction_seed & 0xffffffff) ^
             ((construction_seed >> 32) & 0xffffffff)) & 0xffffffff),
        "construction_seed": f"0x{construction_seed:x}",
        "loss_rate": f"{loss:.17g}",
        "loss_seed": f"0x{loss_seed:x}",
        "schedule": schedule,
        "pmf_sha256": expected_pmf_digest,
        "pmf_encoding": TARGET_CONTRACT["pmf_encoding"],
        "staircase_scale": scale_spec,
    }
    if target_receipt != expected_target_receipt:
        differences = sorted(
            set(target_receipt) ^ set(expected_target_receipt))
        mismatches = sorted(
            name for name in set(target_receipt) & set(expected_target_receipt)
            if target_receipt[name] != expected_target_receipt[name])
        raise MeasurementError(
            f"wh2_target receipt does not match the requested trial for "
            f"K={block_count}: fields={differences}, values={mismatches}")

    try:
        metadata_seed = int(compare_metadata["seed"], 0)
        schedule_seed = int(compare_metadata["schedule_seed"], 0)
        metadata_loss = float(compare_metadata["loss"])
        dynamic_matches = (
            compare_metadata["N"] == f"[{block_count},{block_count}]" and
            compare_metadata["trials/bb"] == str(trials) and
            compare_metadata["max_message_mib"] == "0" and
            compare_metadata["schedule"] == schedule and
            compare_metadata["loss_trace"] == "common-id-v2" and
            metadata_seed == loss_seed and schedule_seed == loss_seed and
            math.isfinite(metadata_loss) and metadata_loss == loss
        )
    except (KeyError, TypeError, ValueError):
        dynamic_matches = False
    if not dynamic_matches:
        raise MeasurementError(
            f"compare metadata does not match the requested trial for "
            f"K={block_count}")
    expected_metadata = {
        "N", "trials/bb", "loss", "seed", "max_message_mib", "schedule",
        "schedule_seed", "loss_trace",
        *PRODUCTION_COMPLETION_REGIME,
        *_COMPARE_BANNER_ARM,
    }
    if set(compare_metadata) != expected_metadata:
        raise MeasurementError(
            f"compare returned unexpected metadata fields for "
            f"K={block_count}: "
            f"{sorted(set(compare_metadata) ^ expected_metadata)}")
    completion = {
        name: compare_metadata.get(name)
        for name in PRODUCTION_COMPLETION_REGIME
    }
    if completion != PRODUCTION_COMPLETION_REGIME:
        raise MeasurementError(
            f"compare did not use production mixed completion for "
            f"K={block_count}: {completion}")
    compare_arm = {
        name: compare_metadata.get(name)
        for name in _COMPARE_BANNER_ARM
    }
    if compare_arm != _COMPARE_BANNER_ARM:
        raise MeasurementError(
            f"compare did not use the production WH2 arm for "
            f"K={block_count}: {compare_arm}")
    if row_codecs != ["baseline", "v2", "v2_target"]:
        raise MeasurementError(
            f"compare returned unexpected row sequence {row_codecs!r} "
            f"for K={block_count}")
    if len(rows) != 1:
        raise MeasurementError(
            f"compare returned {len(rows)} v2_target rows for K={block_count}")
    return rows[0]


def paired_probe(
        bench, block_count, block_bytes, candidate_pmf, degree_scale=None, *,
        native_profile, target_profile, seed_policy, loss, schedule,
        construction_seed, loss_seed, warmup_replicates, replicates,
        inner_reps, max_overhead, cache_state, evict_bytes, context,
    required_margin):
    """Run one complete same-process candidate/control plus A/A experiment."""
    validate_peeltiming_dimensions(
        block_count=block_count,
        block_bytes=block_bytes,
        target_profile=target_profile,
        seed_policy=seed_policy,
        construction_seed=construction_seed,
        loss=loss,
        loss_seed=loss_seed,
        schedule=schedule,
        warmup_replicates=warmup_replicates,
        replicates=replicates,
        inner_reps=inner_reps,
        max_overhead=max_overhead,
        cache_state=cache_state,
        evict_bytes=evict_bytes,
        required_margin=required_margin,
    )
    if (not isinstance(native_profile, StockProfile) or
            native_profile.block_count != block_count or
            native_profile.target_profile != target_profile):
        raise ValueError("native_profile does not match peeltiming target")
    try:
        numeric_pmf = []
        for value in candidate_pmf:
            if isinstance(value, bool):
                raise ValueError
            numeric_pmf.append(float(value))
    except (OverflowError, TypeError, ValueError):
        raise ValueError("invalid paired candidate PMF")
    if (not 1 <= len(numeric_pmf) <= 64 or
            any(not math.isfinite(value) or value < 0.0
                for value in numeric_pmf) or
            not math.isfinite(sum(numeric_pmf)) or
            not sum(numeric_pmf) > 0.0):
        raise ValueError("invalid paired candidate PMF")
    candidate_spec = _canonical_paired_pmf_spec(numeric_pmf)
    scale_spec = "identity"
    if degree_scale is not None:
        _validate_paired_degree_scale(degree_scale)
        scale_spec = _canonical_staircase_scale_spec(degree_scale)

    initial_context = (
        read_paired_context(context, require_thermal=False)
        if isinstance(context, str) else context)
    frozen_context, thermal_capture = _prepare_paired_context(
        initial_context, cache_state=cache_state, evict_bytes=evict_bytes)
    context_sha256 = _validate_paired_bound_context(
        frozen_context, cache_state=cache_state, evict_bytes=evict_bytes,
        require_capture=True)
    command = [
        bench, "peeltiming",
        "--N", str(block_count),
        "--bb", str(block_bytes),
        "--target-profile", target_profile,
        "--seed-policy", seed_policy,
        "--construction-seed", str(construction_seed),
        "--loss", f"{float(loss):.17g}",
        "--loss-seed", str(loss_seed),
        "--schedule", schedule,
        "--candidate-pmf", candidate_spec,
        "--candidate-scale", scale_spec,
        "--warmup-replicates", str(warmup_replicates),
        "--replicates", str(replicates),
        "--inner-reps", str(inner_reps),
        "--max-overhead", str(max_overhead),
        "--cache-state", cache_state,
        "--evict-bytes", str(evict_bytes),
        "--context-sha256", context_sha256,
        "--required-margin", f"{float(required_margin):.17g}",
    ]
    try:
        stdout = _run_checked(command, isolated_codec_env())
        final_context = _finalize_paired_context(
            frozen_context, stdout, thermal_capture)
    finally:
        os.close(thermal_capture.fd)
    return parse_peeltiming_output(
        stdout,
        block_count=block_count,
        block_bytes=block_bytes,
        candidate_pmf=numeric_pmf,
        degree_scale=degree_scale,
        native_profile=native_profile,
        target_profile=target_profile,
        seed_policy=seed_policy,
        construction_seed=construction_seed,
        loss=loss,
        loss_seed=loss_seed,
        schedule=schedule,
        warmup_replicates=warmup_replicates,
        replicates=replicates,
        inner_reps=inner_reps,
        max_overhead=max_overhead,
        cache_state=cache_state,
        evict_bytes=evict_bytes,
        context=final_context,
        required_margin=float(required_margin),
    )


def make_table_document(
        entries, *, generator, bench, construction_seed_base, loss_seed_base,
        target_profile, seed_policy, loss, schedule, settings,
        source_provenance=None, artifact_identity=None):
    """Wrap entries in the strict, replayable peel-table schema."""
    if (not isinstance(generator, str) or not generator or
            not isinstance(settings, dict) or
            target_profile != TARGET_PROFILE or
            seed_policy != TARGET_SEED_POLICY or
            not isinstance(schedule, str) or
            schedule not in TARGET_SCHEDULES or
            not valid_loss_rate(loss) or
            any(
                isinstance(seed, bool) or not isinstance(seed, int) or
                not 0 <= seed <= 0xffffffffffffffff
                for seed in (construction_seed_base, loss_seed_base))):
        raise ValueError(
            "invalid table generator, target, settings, loss, or seeds")
    if artifact_identity is None:
        identity = _artifact_identity(bench, generator)
    else:
        verify_artifact_identity(artifact_identity, bench, generator)
        identity = artifact_identity
    document = {
        "schema": PEEL_TABLE_SCHEMA,
        "provenance": {
            "generator": generator,
            "benchmark": identity["benchmark"],
            "source": identity["source"],
            "native_pmf": NATIVE_PMF_PROTOCOL,
            "native_compare": NATIVE_COMPARE_PROTOCOL,
            "native_paired": NATIVE_PAIRED_PROTOCOL,
            "seed_derivation": SEED_DERIVATION_PROTOCOL,
            "completion_regime": {
                "protocol": COMPLETION_REGIME_PROTOCOL,
                "settings": dict(PRODUCTION_COMPLETION_REGIME),
            },
            "target_contract": dict(TARGET_CONTRACT),
            "compare_arm": dict(PRODUCTION_COMPARE_ARM),
            "measurement_policy": {
                "target_profile": target_profile,
                "seed_policy": seed_policy,
                "loss": float(loss),
                "schedule": schedule,
                "loss_model": PAIRED_LOSS_MODEL,
            },
            "recovery_metric_scope": dict(RECOVERY_METRIC_SCOPE),
            "python_runtime": _python_runtime_identity(),
            "construction_seed_base": construction_seed_base,
            "loss_seed_base": loss_seed_base,
            "settings": settings,
        },
        "entries": {
            str(k): value for k, value in sorted(entries.items())
        },
    }
    if source_provenance is not None:
        document["source_provenance"] = source_provenance
    validate_table_document(document)
    return document


def _is_sha256(value):
    return (
        isinstance(value, str) and len(value) == 64 and
        all(character in "0123456789abcdef" for character in value)
    )


_TARGET_CONTRACT_INTEGER_FIELDS = {
    "precode_contract", "packet_contract", "dense_rows", "heavy_rows",
    "gf256_rows", "gf16_rows", "period", "residue_skew", "mix_count",
}
_NATIVE_INTEGER_FIELDS = {
    "block_count", "precode_contract", "packet_contract", "staircase",
    "dense_rows", "heavy_rows", "source_hits", "gf256_rows", "gf16_rows",
    "period", "residue_skew", "mix_count",
}
_PROXY_REGIME_INTEGER_FIELDS = {
    "solve_block_bytes", "cost_model_block_bytes", "cost_model_verified",
    "band_tracking_x", "seed_base", "period", "gf16_rows",
}


def _validate_target_contract_object(contract, label):
    """Require the canonical contract values and their canonical JSON types."""
    if (not isinstance(contract, dict) or
            set(contract) != set(TARGET_CONTRACT) or
            any(type(contract.get(name)) is not int
                for name in _TARGET_CONTRACT_INTEGER_FIELDS) or
            contract != TARGET_CONTRACT):
        raise MeasurementError(f"{label} is not the exact dispatch-v1 contract")


def _validate_proxy_measure_regime(regime, label):
    """Reject bool/float aliases in the fixed integer proxy receipt."""
    if (not isinstance(regime, dict) or
            set(regime) != set(PROXY_MEASURE_REGIME) or
            any(type(regime.get(name)) is not int
                for name in _PROXY_REGIME_INTEGER_FIELDS) or
            regime != PROXY_MEASURE_REGIME):
        raise MeasurementError(f"{label} is not the exact proxy regime")


_TARGET_RECEIPT_FIELDS = {
    "profile", "contract_id", "contract_sha256", "precode_contract",
    "packet_contract", "architecture", *TARGET_MEASUREMENT_POLICY,
    "N", "bb", "staircase", "dense_rows", "heavy_rows", "completion",
    "gf256_rows", "gf16_rows", "period", "geometry", "residue_schedule",
    "residue_skew", "source_hits", "target_mean", "mix_count",
    "packet_peel_seed",
    "construction_seed", "loss_rate", "loss_seed", "schedule",
    "pmf_sha256", "pmf_encoding", "staircase_scale",
}


def _validate_target_receipt(
        receipt, *, block_count, block_bytes, native, construction_seed,
        loss_seed, loss, schedule, pmf_digest, staircase_scale, label):
    """Validate every field in one native ``# wh2_target`` receipt."""
    if not isinstance(receipt, dict) or set(receipt) != _TARGET_RECEIPT_FIELDS:
        raise MeasurementError(f"{label} is incomplete")
    expected = {
        "profile": TARGET_PROFILE,
        "contract_id": TARGET_CONTRACT["contract_id"],
        "contract_sha256": TARGET_CONTRACT["contract_sha256"],
        "precode_contract": str(TARGET_CONTRACT["precode_contract"]),
        "packet_contract": str(TARGET_CONTRACT["packet_contract"]),
        "architecture": TARGET_CONTRACT["architecture"],
        **TARGET_MEASUREMENT_POLICY,
        "N": str(block_count),
        "bb": str(block_bytes),
        "staircase": str(native["staircase"]),
        "dense_rows": str(native["dense_rows"]),
        "heavy_rows": str(native["heavy_rows"]),
        "completion": native["completion"],
        "gf256_rows": str(native["gf256_rows"]),
        "gf16_rows": str(native["gf16_rows"]),
        "period": str(native["period"]),
        "geometry": native["geometry"],
        "residue_schedule": native["residue_schedule"],
        "residue_skew": str(native["residue_skew"]),
        "source_hits": str(native["source_hits"]),
        "target_mean": _target_mean_spec(
            block_count, native["staircase"], native["source_hits"],
            staircase_scale),
        "mix_count": str(native["mix_count"]),
        "packet_peel_seed": str(
            ((construction_seed & 0xffffffff) ^
             ((construction_seed >> 32) & 0xffffffff)) & 0xffffffff),
        "construction_seed": f"0x{construction_seed:x}",
        "loss_rate": f"{float(loss):.17g}",
        "loss_seed": f"0x{loss_seed:x}",
        "schedule": schedule,
        "pmf_sha256": pmf_digest,
        "pmf_encoding": TARGET_CONTRACT["pmf_encoding"],
        "staircase_scale": staircase_scale,
    }
    if receipt != expected:
        raise MeasurementError(f"{label} does not match its exact target arm")


def _validate_recovery_metrics(receipt, trials, label):
    if (not isinstance(receipt, dict) or
            not RECOVERY_METRIC_FIELDS.issubset(receipt)):
        raise MeasurementError(f"{label} is missing recovery metrics")
    fail = receipt["fail"]
    if (any(
            type(receipt[name]) is not int or
            not 0 <= receipt[name] <= 0xffffffffffffffff
            for name in ("construction_seed", "loss_seed")) or
            not isinstance(receipt["target_receipt"], dict) or
            set(receipt["target_receipt"]) != _TARGET_RECEIPT_FIELDS or
            type(fail) is not int or
            not 0 <= fail <= trials):
        raise MeasurementError(
            f"{label} has invalid construction/loss seeds, target, or fail")
    value_names = (
        "oh_mean", "OH_sd", "OH50", "OH95", "OH99", "OH_max",
        "decode_mbps")
    if any(type(receipt[name]) is not float for name in value_names):
        raise MeasurementError(f"{label} has non-canonical recovery metrics")
    values = [receipt[name] for name in value_names]
    if (any(not _is_canonical_nonnegative_float(value)
            for value in values) or
            not values[2] <= values[3] <= values[4] <= values[5] or
            values[0] > values[5]):
        raise MeasurementError(f"{label} has invalid recovery metrics")


def _validate_paired_arm_summary(receipt, trials, label):
    if (not isinstance(receipt, dict) or
            not PAIRED_ARM_FIELDS.issubset(receipt) or
            type(receipt.get("trials")) is not int or
            receipt["trials"] != trials or
            type(receipt.get("fail")) is not int or
            not 0 <= receipt["fail"] <= trials):
        raise MeasurementError(f"{label} has invalid paired counts")
    names = (
        "oh_mean", "OH_sd", "OH50", "OH95", "OH99", "OH_max",
        "solve_mbps",
    )
    if any(type(receipt.get(name)) is not float for name in names):
        raise MeasurementError(
            f"{label} has non-canonical paired summaries")
    values = [receipt[name] for name in names]
    if (any(not _is_canonical_nonnegative_float(value)
            for value in values) or
            values[0] > values[5] or
            not values[2] <= values[3] <= values[4] <= values[5] or
            values[6] <= 0.0):
        raise MeasurementError(f"{label} has invalid paired summaries")


def _canonical_goodput(metrics, block_count):
    """Recompute the exact goodput represented by recovery metrics."""
    if metrics["fail"] != 0:
        return 0.0
    return (
        metrics["decode_mbps"] * block_count /
        (block_count + metrics["oh_mean"])
    )


def _validate_search_context(
        context, block_count, search_kind, native, sampling_seed, label):
    """Validate every search-control field needed to replay one result."""
    if search_kind == "direct-real-codec":
        required = {
            "warm_start", "sampling_seed", "screen", "refine",
            "gate_trials", "gate_block_bytes", "screen_paired_replicates",
            "rank_top",
        }
        warm_start = context.get("warm_start")
        valid_warm_start = warm_start is None
        if isinstance(warm_start, list) and len(warm_start) == 4 and all(
                type(value) is int for value in warm_start):
            valid_warm_start = (
                0 <= warm_start[0] <= 400 and
                -100 <= warm_start[1] <= 1600 and
                2 <= warm_start[2] <= 64 and
                0 <= warm_start[3] <= 100
            )
        integer_bounds = {
            "screen": (1, None),
            "refine": (0, None),
            "gate_trials": (1, _COMPARE_TRIALS_MAX),
            "gate_block_bytes": (2, _COMPARE_BLOCK_BYTES_MAX),
            "screen_paired_replicates": (4, 10000),
            "rank_top": (1, None),
        }
        invalid_integer = any(
            type(context.get(name)) is not int or
            context[name] < minimum or
            (maximum is not None and context[name] > maximum)
            for name, (minimum, maximum) in integer_bounds.items()
        ) if set(context) == required else True
        if (set(context) != required or invalid_integer or
                context["gate_block_bytes"] % 2 != 0 or
                context["screen_paired_replicates"] % 2 != 0 or
                context.get("sampling_seed") != sampling_seed or
                not valid_warm_start):
            raise MeasurementError(f"{label} has invalid direct-search context")
        return
    if search_kind != "unverified-proxy-funnel":
        raise MeasurementError(f"{label} has an unknown search context")
    required = {
        "proxy_cost_model", "proxy_measure_regime", "proxy_ordering",
        "search_box", "scale_centi", "warm_start", "sampling_seed",
        "screen", "refine", "finals", "screen_cells", "gate_trials",
        "gate_block_bytes", "rank_top", "threads", "batch", "cell_base",
    }
    if set(context) != required:
        raise MeasurementError(f"{label} has incomplete funnel context")
    _validate_proxy_measure_regime(
        context.get("proxy_measure_regime"),
        f"{label}.proxy_measure_regime")
    expected_scale_max = int(math.ceil(
        100.0 * min(
            float(block_count),
            max(40.0, 4.0 * float(native["target_mean"])))
    ))
    integer_bounds = {
        "screen": (1, None),
        "refine": (0, None),
        "finals": (2, None),
        "screen_cells": (1, None),
        "gate_trials": (1, _COMPARE_TRIALS_MAX),
        "gate_block_bytes": (2, _COMPARE_BLOCK_BYTES_MAX),
        "rank_top": (1, None),
        "threads": (1, 1024),
        "batch": (1, None),
        "cell_base": (0, 0xffffffffffffffff),
    }
    invalid_integer = False
    for name, (minimum, maximum) in integer_bounds.items():
        value = context.get(name)
        if (type(value) is not int or
                value < minimum or
                (maximum is not None and value > maximum)):
            invalid_integer = True
            break
    warm_start = context.get("warm_start")
    valid_warm_start = warm_start is None
    if isinstance(warm_start, list) and len(warm_start) == 5 and all(
            type(value) is int for value in warm_start):
        valid_warm_start = (
            0 <= warm_start[0] <= expected_scale_max and
            0 <= warm_start[1] <= 400 and
            0 <= warm_start[2] <= 1700 and
            2 <= warm_start[3] <= 64 and
            0 <= warm_start[4] <= 100
        )
    scale_centi = context.get("scale_centi")
    valid_scale_centi = (
        isinstance(scale_centi, list) and
        len(scale_centi) == 2 and
        all(type(value) is int for value in scale_centi)
    )
    if (invalid_integer or
            context["cell_base"] + context["screen_cells"] > 0x100000000 or
            context["gate_block_bytes"] % 2 != 0 or
            context.get("proxy_cost_model") != PROXY_COST_MODEL or
            context.get("proxy_ordering") != PROXY_ORDERING_PROTOCOL or
            context.get("search_box") != SEARCH_BOX_PROTOCOL or
            not valid_scale_centi or
            scale_centi != [0, expected_scale_max] or
            context.get("sampling_seed") != sampling_seed or
            not valid_warm_start):
        raise MeasurementError(f"{label} has invalid funnel context")


def _validate_legacy_search_receipt(
        receipt, block_count, expected_construction_seed_base,
        expected_loss_seed_base, expected_seed_domain, expected_search_kind,
        native, measurement_policy, label):
    required = {
        "mode", "goodput", "trials", "block_bytes", "search_kind",
        "construction_seed_base", "loss_seed_base", "seed_domain",
        "coordinates", "peel_pmf", "peel_pmf_sha256", "shipped_control",
        "shipped_goodput", "context",
    }
    if (not isinstance(receipt, dict) or
            set(receipt) != required | RECOVERY_METRIC_FIELDS):
        raise MeasurementError(f"{label} is incomplete")
    trials = receipt["trials"]
    block_bytes = receipt["block_bytes"]
    if (receipt["mode"] not in ("trained", "scale-only", "shipped") or
            type(trials) is not int or
            not 1 <= trials <= _COMPARE_TRIALS_MAX or
            type(block_bytes) is not int or
            not 1 <= block_bytes <= _COMPARE_BLOCK_BYTES_MAX or
            block_bytes % 2 != 0 or
            receipt["search_kind"] != expected_search_kind or
            any(
                type(receipt[name]) is not int or
                receipt[name] != expected or
                not 0 <= receipt[name] <= 0xffffffffffffffff
                for name, expected in (
                    ("construction_seed_base",
                     expected_construction_seed_base),
                    ("loss_seed_base", expected_loss_seed_base))) or
            receipt["seed_domain"] != expected_seed_domain or
            not isinstance(receipt["coordinates"], dict) or
            not isinstance(receipt["context"], dict)):
        raise MeasurementError(f"{label} has invalid search metadata")
    expected_sampling_seed = derive_seed(
        expected_construction_seed_base, expected_seed_domain, block_count,
        "sampling")
    sampling_seed = receipt["context"].get("sampling_seed")
    if (type(sampling_seed) is not int or
            sampling_seed != expected_sampling_seed):
        raise MeasurementError(
            f"{label} has an invalid derived sampling seed")
    _validate_search_context(
        receipt["context"], block_count, receipt["search_kind"], native,
        expected_sampling_seed, f"{label}.context")
    coordinate_names = ("scale", "p1", "tilt", "dmax", "absorb")
    coordinates = receipt["coordinates"]
    if set(coordinates) != set(coordinate_names):
        raise MeasurementError(f"{label} has invalid search coordinates")
    scale_value = coordinates["scale"]
    integer_values = [
        coordinates[name] for name in ("p1", "tilt", "dmax", "absorb")
    ]
    if (type(scale_value) is not float or
            any(type(value) is not int for value in integer_values)):
        raise MeasurementError(f"{label} has non-numeric search coordinates")
    scale = scale_value
    p1 = coordinates["p1"]
    tilt = coordinates["tilt"]
    dmax = coordinates["dmax"]
    absorb = coordinates["absorb"]
    if (not _is_canonical_scale(scale) or
            not 0 <= p1 <= 400 or
            not -100 <= tilt <= 1600 or
            not 2 <= dmax <= 64 or
            not 0 <= absorb <= 100 or
            (expected_search_kind == "direct-real-codec" and
             scale != -1.0) or
            (expected_search_kind == "unverified-proxy-funnel" and
             receipt["mode"] in ("trained", "scale-only") and
             (not receipt["context"]["scale_centi"][0] <=
              int(round(scale * 100.0)) <=
              receipt["context"]["scale_centi"][1] or
              scale != int(round(scale * 100.0)) / 100.0))):
        raise MeasurementError(f"{label} has invalid search coordinates")
    peel_pmf_digest = _validate_pmf(
        receipt["peel_pmf"], f"{label}.peel_pmf")
    if receipt["peel_pmf_sha256"] != peel_pmf_digest:
        raise MeasurementError(f"{label} has an invalid peel PMF digest")
    if receipt["mode"] in ("shipped", "scale-only"):
        expected_pmf = list(native["pmf"])
        expected_scale = -1.0 if receipt["mode"] == "shipped" else scale
        if (scale != expected_scale or p1 != 100 or tilt != 0 or
                dmax != 64 or absorb != 100):
            raise MeasurementError(
                f"{label} has invalid stock-control coordinates")
        if receipt["mode"] == "scale-only" and scale < 0.0:
            raise MeasurementError(
                f"{label} has an unset scale-only coordinate")
    else:
        expected_pmf = family(native["pmf"], p1, tilt, dmax, absorb)
    if expected_pmf is None or receipt["peel_pmf"] != expected_pmf:
        raise MeasurementError(
            f"{label} peel PMF does not match its coordinates")
    # Funnel canonicalizes a native-shaped PMF with an active staircase scale
    # to the stock-hook scale-only arm. A trained receipt for the same PMF is
    # therefore not an artifact the current generators can produce.
    if (receipt["mode"] == "trained" and
            _native_peel_cdf_signature(
                receipt["peel_pmf"], block_count) ==
            _native_peel_cdf_signature(native["pmf"], block_count)):
        raise MeasurementError(
            f"{label} labels a realized stock-PMF alias as a trained arm")
    if type(receipt["goodput"]) is not float:
        raise MeasurementError(f"{label} has invalid goodput")
    goodput = receipt["goodput"]
    if not _is_canonical_nonnegative_float(goodput):
        raise MeasurementError(f"{label} has invalid goodput")
    _validate_recovery_metrics(receipt, trials, label)
    if receipt["fail"] != 0:
        raise MeasurementError(f"{label} selected a non-decoding search arm")
    _validate_recovery_metrics(
        receipt["shipped_control"], trials, f"{label}.shipped_control")
    if set(receipt["shipped_control"]) != RECOVERY_METRIC_FIELDS:
        raise MeasurementError(
            f"{label} shipped control has unexpected fields")
    expected_construction_seed = derive_seed(
        expected_construction_seed_base, receipt["seed_domain"], block_count,
        "rank", trials, block_bytes, "construction")
    expected_loss_seed = derive_seed(
        expected_loss_seed_base, receipt["seed_domain"], block_count, "rank",
        trials, block_bytes, "loss")
    if (receipt["construction_seed"] != expected_construction_seed or
            receipt["loss_seed"] != expected_loss_seed):
        raise MeasurementError(
            f"{label} has invalid derived construction/loss seeds")
    if (receipt["shipped_control"]["construction_seed"] !=
            expected_construction_seed or
            receipt["shipped_control"]["loss_seed"] != expected_loss_seed):
        raise MeasurementError(
            f"{label} control did not use the paired construction/loss seeds")
    scale_spec = (
        "unset" if scale == -1.0 else
        _canonical_staircase_scale_spec(scale))
    _validate_target_receipt(
        receipt["target_receipt"],
        block_count=block_count,
        block_bytes=block_bytes,
        native=native,
        construction_seed=expected_construction_seed,
        loss_seed=expected_loss_seed,
        loss=measurement_policy["loss"],
        schedule=measurement_policy["schedule"],
        pmf_digest=(
            STOCK_PMF_DIGEST
            if receipt["mode"] in ("shipped", "scale-only") else
            pmf_sha256(receipt["peel_pmf"])),
        staircase_scale=scale_spec,
        label=f"{label}.target_receipt",
    )
    _validate_target_receipt(
        receipt["shipped_control"]["target_receipt"],
        block_count=block_count,
        block_bytes=block_bytes,
        native=native,
        construction_seed=expected_construction_seed,
        loss_seed=expected_loss_seed,
        loss=measurement_policy["loss"],
        schedule=measurement_policy["schedule"],
        pmf_digest=STOCK_PMF_DIGEST,
        staircase_scale="unset",
        label=f"{label}.shipped_control.target_receipt",
    )
    candidate_target = dict(receipt["target_receipt"])
    control_target = dict(receipt["shipped_control"]["target_receipt"])
    for field in ("pmf_sha256", "staircase_scale", "target_mean"):
        candidate_target.pop(field)
        control_target.pop(field)
    if candidate_target != control_target:
        raise MeasurementError(
            f"{label} candidate/control changed more than the PMF or scale")
    expected_goodput = _canonical_goodput(receipt, block_count)
    if goodput != expected_goodput:
        raise MeasurementError(f"{label} has an inconsistent goodput")
    shipped = receipt["shipped_control"]
    expected_shipped_goodput = _canonical_goodput(shipped, block_count)
    if (type(receipt["shipped_goodput"]) is not float or
            receipt["shipped_goodput"] != expected_shipped_goodput):
        raise MeasurementError(
            f"{label} has an inconsistent shipped-control goodput")
    if (receipt["mode"] in ("trained", "scale-only") and
            not expected_goodput > expected_shipped_goodput):
        raise MeasurementError(
            f"{label} selected a candidate that did not beat shipped")
    if receipt["mode"] == "shipped" and any(
            receipt[name] != shipped[name]
            for name in (
                "construction_seed", "loss_seed", "target_receipt", "fail",
                "oh_mean", "OH_sd", "OH50", "OH95", "OH99", "OH_max",
                "decode_mbps")):
        raise MeasurementError(
            f"{label} shipped selection contradicts its control metrics")


def _validate_search_receipt(
        receipt, block_count, expected_construction_seed_base,
        expected_loss_seed_base, expected_seed_domain, expected_search_kind,
        native, measurement_policy, label):
    """Replay and cross-check one v4 native paired search receipt."""
    required = {
        "protocol", "mode", "selected_arm", "goodput", "shipped_goodput",
        "trials", "block_bytes", "search_kind",
        "construction_seed_base", "loss_seed_base", "seed_domain",
        "coordinates", "peel_pmf", "peel_pmf_sha256",
        "evaluated_coordinates", "evaluated_pmf",
        "evaluated_pmf_sha256", "paired_measurement",
        "stock_control_source", "stock_control_measurement", "context",
    }
    if not isinstance(receipt, dict) or set(receipt) != required:
        raise MeasurementError(f"{label} is incomplete")
    trials = receipt["trials"]
    block_bytes = receipt["block_bytes"]
    mode = receipt["mode"]
    if (receipt["protocol"] != NATIVE_PAIRED_PROTOCOL or
            mode not in ("trained", "scale-only", "shipped") or
            receipt["selected_arm"] !=
                ("identity" if mode == "shipped" else "candidate") or
            type(trials) is not int or trials < 4 or trials % 2 != 0 or
            type(block_bytes) is not int or
            not 2 <= block_bytes <= _COMPARE_BLOCK_BYTES_MAX or
            block_bytes % 2 != 0 or
            receipt["search_kind"] != expected_search_kind or
            receipt["seed_domain"] != expected_seed_domain or
            any(
                type(receipt.get(name)) is not int or
                receipt[name] != expected or
                not 0 <= receipt[name] <= 0xffffffffffffffff
                for name, expected in (
                    ("construction_seed_base",
                     expected_construction_seed_base),
                    ("loss_seed_base", expected_loss_seed_base))) or
            not isinstance(receipt["context"], dict)):
        raise MeasurementError(f"{label} has invalid paired metadata")
    if (expected_search_kind == "direct-real-codec" and
            mode == "scale-only"):
        raise MeasurementError(
            f"{label} records a scale-only arm for scale-free Direct search")
    expected_sampling_seed = derive_seed(
        expected_construction_seed_base, expected_seed_domain, block_count,
        "sampling")
    if (type(receipt["context"].get("sampling_seed")) is not int or
            receipt["context"]["sampling_seed"] != expected_sampling_seed):
        raise MeasurementError(f"{label} has an invalid sampling seed")
    _validate_search_context(
        receipt["context"], block_count, expected_search_kind, native,
        expected_sampling_seed, f"{label}.context")

    coordinate_names = {"scale", "p1", "tilt", "dmax", "absorb"}

    def validate_coordinates(coordinates, coordinate_label):
        if (not isinstance(coordinates, dict) or
                set(coordinates) != coordinate_names or
                type(coordinates.get("scale")) is not float or
                any(type(coordinates.get(name)) is not int
                    for name in ("p1", "tilt", "dmax", "absorb"))):
            raise MeasurementError(
                f"{coordinate_label} has non-canonical coordinates")
        scale = coordinates["scale"]
        if (not _is_canonical_scale(scale) or
                not 0 <= coordinates["p1"] <= 400 or
                not -100 <= coordinates["tilt"] <= 1600 or
                not 2 <= coordinates["dmax"] <= 64 or
                not 0 <= coordinates["absorb"] <= 100 or
                (expected_search_kind == "direct-real-codec" and
                 scale != -1.0)):
            raise MeasurementError(
                f"{coordinate_label} has out-of-domain coordinates")
        if (expected_search_kind == "unverified-proxy-funnel" and
                scale != -1.0):
            scale_centi = int(round(scale * 100.0))
            if (not receipt["context"]["scale_centi"][0] <= scale_centi <=
                    receipt["context"]["scale_centi"][1] or
                    scale != scale_centi / 100.0):
                raise MeasurementError(
                    f"{coordinate_label} has an off-lattice scale")
        return coordinates

    coordinates = validate_coordinates(
        receipt["coordinates"], f"{label}.coordinates")
    evaluated = validate_coordinates(
        receipt["evaluated_coordinates"],
        f"{label}.evaluated_coordinates")
    selected_pmf_digest = _validate_pmf(
        receipt["peel_pmf"], f"{label}.peel_pmf")
    evaluated_pmf_digest = _validate_pmf(
        receipt["evaluated_pmf"], f"{label}.evaluated_pmf")
    if (receipt["peel_pmf_sha256"] != selected_pmf_digest or
            receipt["evaluated_pmf_sha256"] != evaluated_pmf_digest):
        raise MeasurementError(f"{label} has an invalid PMF digest")

    stock_coordinates = {
        "scale": -1.0, "p1": 100, "tilt": 0,
        "dmax": 64, "absorb": 100,
    }
    evaluated_is_stock_shape = (
        evaluated["p1"] == 100 and evaluated["tilt"] == 0 and
        evaluated["dmax"] == 64 and evaluated["absorb"] == 100
    )
    expected_evaluated_pmf = (
        list(native["pmf"]) if evaluated_is_stock_shape else
        family(
            native["pmf"], evaluated["p1"], evaluated["tilt"],
            evaluated["dmax"], evaluated["absorb"])
    )
    if (expected_evaluated_pmf is None or
            receipt["evaluated_pmf"] != expected_evaluated_pmf):
        raise MeasurementError(
            f"{label} evaluated PMF does not match its coordinates")
    evaluated_semantics_are_stock = (
        _native_peel_cdf_signature(
            receipt["evaluated_pmf"], block_count) ==
        _native_peel_cdf_signature(native["pmf"], block_count)
    )
    if mode == "trained" and evaluated_semantics_are_stock:
        raise MeasurementError(
            f"{label} labels a realized stock-PMF alias as a trained arm")
    if mode == "shipped":
        if (coordinates != stock_coordinates or
                evaluated != stock_coordinates or
                receipt["peel_pmf"] != list(native["pmf"]) or
                receipt["evaluated_pmf"] != list(native["pmf"])):
            raise MeasurementError(
                f"{label} is not the predeclared stock fallback")
    else:
        expected_mode = (
            "scale-only" if evaluated_semantics_are_stock else "trained")
        if (mode != expected_mode or coordinates != evaluated or
                receipt["peel_pmf"] != receipt["evaluated_pmf"] or
                (expected_mode == "scale-only" and
                 not evaluated_is_stock_shape) or
                coordinates["scale"] < (
                    0.0 if expected_search_kind ==
                    "unverified-proxy-funnel" else -1.0)):
            raise MeasurementError(
                f"{label} selected configuration is inconsistent")

    expected_construction_seed = derive_seed(
        expected_construction_seed_base, expected_seed_domain, block_count,
        "rank", trials, block_bytes, "construction")
    expected_loss_seed = derive_seed(
        expected_loss_seed_base, expected_seed_domain, block_count, "rank",
        trials, block_bytes, "loss")
    measurement = _validate_paired_measurement_receipt(
        receipt["paired_measurement"],
        block_count=block_count,
        block_bytes=block_bytes,
        candidate_pmf=receipt["evaluated_pmf"],
        degree_scale=(
            None if evaluated["scale"] == -1.0 else evaluated["scale"]),
        native=native,
        measurement_policy=measurement_policy,
        construction_seed=expected_construction_seed,
        loss_seed=expected_loss_seed,
        label=f"{label}.paired_measurement",
    )
    validate_paired_stock_control_receipt(
        mode=mode,
        source=receipt["stock_control_source"],
        receipt=receipt["stock_control_measurement"],
        selected_measurement=measurement,
        block_count=block_count,
        block_bytes=block_bytes,
        native=native,
        measurement_policy=measurement_policy,
        construction_seed=expected_construction_seed,
        loss_seed=expected_loss_seed,
        label=f"{label}.stock_control",
    )
    selected = (
        measurement.identity if mode == "shipped"
        else measurement.candidate)
    if (measurement.candidate.trials != trials or
            measurement.identity.trials != trials or
            (mode != "shipped" and not measurement.valid_for_promotion) or
            receipt["goodput"] != selected.goodput(block_count) or
            receipt["shipped_goodput"] !=
                measurement.identity.goodput(block_count)):
        raise MeasurementError(
            f"{label} has an ineligible or stale paired selection")
    for name in ("goodput", "shipped_goodput"):
        if not _is_canonical_nonnegative_float(receipt[name]):
            raise MeasurementError(f"{label} has invalid {name}")
    return measurement


def _validate_native_receipt(receipt, block_count, label):
    required = {
        "block_count", "target_profile", "contract_id", "contract_sha256",
        "precode_contract", "packet_contract", "architecture", "staircase",
        "dense_rows", "heavy_rows", "source_hits", "completion",
        "gf256_rows", "gf16_rows", "period", "geometry",
        "residue_schedule", "residue_skew", "mix_count", "target_mean",
        "native_pmf_sha256", "pmf_encoding", "pmf_sha256", "pmf",
    }
    if not isinstance(receipt, dict) or set(receipt) != required:
        raise MeasurementError(f"{label} is incomplete")
    if type(block_count) is not int or not 2 <= block_count <= 64000:
        raise MeasurementError(f"{label} has invalid block count")
    if any(type(receipt.get(name)) is not int
           for name in _NATIVE_INTEGER_FIELDS):
        raise MeasurementError(f"{label} has non-canonical integer metadata")
    if type(receipt["target_mean"]) is not float:
        raise MeasurementError(f"{label} has invalid target mean")
    if (not isinstance(receipt["pmf"], list) or
            any(type(value) is not float for value in receipt["pmf"])):
        raise MeasurementError(f"{label} has a non-canonical native PMF")
    target_mean = float(receipt["target_mean"])
    pmf_digest = _validate_pmf(receipt["pmf"], f"{label}.pmf")
    expected_staircase = _dispatch_staircase_count(block_count)
    expected_source_hits = _dispatch_source_hits(block_count)
    expected_target_mean = (
        block_count * min(expected_source_hits, expected_staircase) /
        expected_staircase
    )
    contract_matches = all(
        receipt.get(name) == value
        for name, value in TARGET_CONTRACT.items())
    if (receipt["block_count"] != block_count or
            type(receipt["block_count"]) is not int or
            not contract_matches or
            type(receipt["staircase"]) is not int or
            receipt["staircase"] != expected_staircase or
            receipt["staircase"] >
            _production_staircase_max(block_count) or
            type(receipt["source_hits"]) is not int or
            receipt["source_hits"] != expected_source_hits or
            not math.isfinite(target_mean) or target_mean <= 0.0 or
            target_mean != expected_target_mean or
            len(receipt["pmf"]) != 64 or
            receipt["native_pmf_sha256"] != STOCK_PMF_DIGEST or
            receipt["pmf_sha256"] != pmf_digest):
        raise MeasurementError(f"{label} has invalid native metadata")


def _validate_pmf(pmf, label):
    if not isinstance(pmf, list) or not 2 <= len(pmf) <= 64:
        raise MeasurementError(f"{label} must contain 2..64 probabilities")
    if any(type(value) is not float for value in pmf):
        raise MeasurementError(
            f"{label} contains non-canonical probabilities")
    values = list(pmf)
    if (any(not _is_canonical_nonnegative_float(value)
            for value in values) or
            abs(sum(values) - 1.0) > 1e-12):
        raise MeasurementError(f"{label} is not a probability distribution")
    return pmf_sha256(values)


def _stock_profile_from_receipt(native):
    return StockProfile(
        block_count=native["block_count"],
        target_profile=native["target_profile"],
        contract_id=native["contract_id"],
        contract_sha256=native["contract_sha256"],
        precode_contract=native["precode_contract"],
        packet_contract=native["packet_contract"],
        architecture=native["architecture"],
        staircase=native["staircase"],
        dense_rows=native["dense_rows"],
        heavy_rows=native["heavy_rows"],
        source_hits=native["source_hits"],
        completion=native["completion"],
        gf256_rows=native["gf256_rows"],
        gf16_rows=native["gf16_rows"],
        period=native["period"],
        geometry=native["geometry"],
        residue_schedule=native["residue_schedule"],
        residue_skew=native["residue_skew"],
        mix_count=native["mix_count"],
        target_mean=native["target_mean"],
        native_pmf_sha256=native["native_pmf_sha256"],
        pmf_encoding=native["pmf_encoding"],
        pmf=tuple(native["pmf"]),
    )


def _validate_paired_measurement_receipt(
        receipt, *, block_count, block_bytes, candidate_pmf, degree_scale,
        native, measurement_policy, construction_seed, loss_seed, label):
    """Replay the strict native parser and require byte-for-byte receipt truth."""
    if not isinstance(receipt, dict):
        raise MeasurementError(f"{label} is not an object")
    manifest = receipt.get("manifest")
    context = receipt.get("context")
    if not isinstance(manifest, dict) or not isinstance(context, dict):
        raise MeasurementError(f"{label} is incomplete")
    integer_names = (
        "warmup_replicates", "replicates", "inner_reps", "max_overhead",
        "evict_bytes",
    )
    if any(type(manifest.get(name)) is not int for name in integer_names):
        raise MeasurementError(f"{label} has invalid manifest dimensions")
    required_margin = manifest.get("required_margin")
    if type(required_margin) is not float:
        raise MeasurementError(f"{label} has invalid required margin")
    try:
        parsed = parse_peeltiming_output(
            receipt.get("native_stdout"),
            block_count=block_count,
            block_bytes=block_bytes,
            candidate_pmf=candidate_pmf,
            degree_scale=degree_scale,
            native_profile=_stock_profile_from_receipt(native),
            target_profile=measurement_policy["target_profile"],
            seed_policy=measurement_policy["seed_policy"],
            construction_seed=construction_seed,
            loss=measurement_policy["loss"],
            loss_seed=loss_seed,
            schedule=measurement_policy["schedule"],
            warmup_replicates=manifest["warmup_replicates"],
            replicates=manifest["replicates"],
            inner_reps=manifest["inner_reps"],
            max_overhead=manifest["max_overhead"],
            cache_state=manifest.get("cache_state"),
            evict_bytes=manifest["evict_bytes"],
            context=context,
            required_margin=required_margin,
        )
    except (KeyError, TypeError, ValueError) as error:
        raise MeasurementError(f"{label} could not be replayed: {error}")
    if not _same_typed_json(parsed.as_dict(), receipt):
        raise MeasurementError(f"{label} aggregate does not match native rows")
    return parsed


def validate_paired_stock_control_receipt(
        *, mode, source, receipt, selected_measurement, block_count,
        block_bytes, native, measurement_policy, construction_seed,
        loss_seed, label):
    """Replay the mandatory fixed stock-vs-stock rank measurement."""
    if mode == "shipped":
        if source != STOCK_CONTROL_SELECTED or receipt is not None:
            raise MeasurementError(
                f"{label} does not reference its selected stock control")
        control = selected_measurement
    else:
        if source != STOCK_CONTROL_EMBEDDED or not isinstance(receipt, dict):
            raise MeasurementError(
                f"{label} does not embed its independent stock control")
        control = _validate_paired_measurement_receipt(
            receipt,
            block_count=block_count,
            block_bytes=block_bytes,
            candidate_pmf=native["pmf"],
            degree_scale=None,
            native=native,
            measurement_policy=measurement_policy,
            construction_seed=construction_seed,
            loss_seed=loss_seed,
            label=f"{label}.paired_measurement",
        )
        if (control.manifest["started_monotonic_ns"] <=
                selected_measurement.finished_monotonic_ns):
            raise MeasurementError(
                f"{label} did not run after the selected rank measurement")
        dimension_fields = (
            "warmup_replicates", "replicates", "inner_reps",
            "max_overhead", "cache_state", "evict_bytes",
            "required_margin",
        )
        if any(
                not _same_typed_json(
                    control.manifest.get(name),
                    selected_measurement.manifest.get(name))
                for name in dimension_fields):
            raise MeasurementError(
                f"{label} dimensions differ from the selected rank run")
        selected_bound = dict(selected_measurement.context["bound"])
        control_bound = dict(control.context["bound"])
        for name in (
                "thermal_prelaunch_monotonic_s",
                "thermal_prelaunch_row_sha256"):
            selected_bound.pop(name, None)
            control_bound.pop(name, None)
        if not _same_typed_json(control_bound, selected_bound):
            raise MeasurementError(
                f"{label} host and sampler context differs from the "
                "selected rank run")
        semantic_fields = (
            "construction_seed", "loss_seed", "trace_sha256",
        )
        replicate_fields = (
            "replicate", "construction_seed", "loss_seed", "trace_sha256",
            "identity_overhead", "identity_recovery_result",
            "identity_recovery_ok",
        )
        identity_fields = (
            "trials", "fail", "oh_mean", "OH_sd", "OH50", "OH95",
            "OH99", "OH_max",
        )
        selected_semantic = {
            name: selected_measurement.semantic[name]
            for name in semantic_fields
        }
        control_semantic = {
            name: control.semantic[name] for name in semantic_fields
        }
        selected_replicates = [
            {name: item[name] for name in replicate_fields}
            for item in selected_measurement.replicate_receipts
        ]
        control_replicates = [
            {name: item[name] for name in replicate_fields}
            for item in control.replicate_receipts
        ]
        selected_identity = {
            name: selected_measurement.identity.as_dict()[name]
            for name in identity_fields
        }
        control_identity = {
            name: control.identity.as_dict()[name]
            for name in identity_fields
        }
        if (not _same_typed_json(control_semantic, selected_semantic) or
                not _same_typed_json(
                    control_replicates, selected_replicates) or
                not _same_typed_json(control_identity, selected_identity)):
            raise MeasurementError(
                f"{label} loss trace or identity recovery differs from the "
                "selected rank run")
    require_paired_stock_control(control, label)
    return control


def _validate_legacy_validation_receipt(
        receipt, block_count, construction_seed_base, loss_seed_base, native,
        measurement_policy, label):
    required = {
        "verdict", "margin_percent", "trials", "block_bytes",
        "scale", "trained_pmf_sha256", "trained_goodput",
        "shipped_goodput", "trained", "shipped",
    }
    if not isinstance(receipt, dict) or set(receipt) != required:
        raise MeasurementError(f"{label} is incomplete")
    trials = receipt["trials"]
    block_bytes = receipt["block_bytes"]
    if (receipt["verdict"] not in ("keep", "control") or
            type(trials) is not int or
            not 1 <= trials <= _COMPARE_TRIALS_MAX or
            type(block_bytes) is not int or
            not 1 <= block_bytes <= _COMPARE_BLOCK_BYTES_MAX or
            block_bytes % 2 != 0):
        raise MeasurementError(f"{label} has invalid validation metadata")
    value_names = (
        "margin_percent", "scale", "trained_goodput", "shipped_goodput")
    if any(type(receipt[name]) is not float for name in value_names):
        raise MeasurementError(
            f"{label} has non-canonical validation values")
    values = [receipt[name] for name in value_names]
    if (not _is_canonical_nonnegative_float(values[0]) or
            not _is_canonical_scale(values[1]) or
            any(not _is_canonical_nonnegative_float(value)
                for value in values[2:]) or
            not _is_sha256(receipt.get("trained_pmf_sha256"))):
        raise MeasurementError(f"{label} has invalid validation values")
    _validate_recovery_metrics(receipt["trained"], trials, f"{label}.trained")
    _validate_recovery_metrics(receipt["shipped"], trials, f"{label}.shipped")
    if (set(receipt["trained"]) != RECOVERY_METRIC_FIELDS or
            set(receipt["shipped"]) != RECOVERY_METRIC_FIELDS):
        raise MeasurementError(f"{label} has unexpected arm fields")
    if (receipt["trained"]["construction_seed"] !=
            receipt["shipped"]["construction_seed"] or
            receipt["trained"]["loss_seed"] !=
            receipt["shipped"]["loss_seed"]):
        raise MeasurementError(
            f"{label} did not use paired construction/loss seeds")
    if receipt["trained"]["fail"] > 0 and receipt["shipped"]["fail"] > 0:
        raise MeasurementError(f"{label} records two failed recovery arms")
    expected_construction_seed = derive_seed(
        construction_seed_base, "validation", block_count, trials,
        block_bytes, "construction")
    expected_loss_seed = derive_seed(
        loss_seed_base, "validation", block_count, trials, block_bytes, "loss")
    if (receipt["trained"]["construction_seed"] !=
            expected_construction_seed or
            receipt["trained"]["loss_seed"] != expected_loss_seed):
        raise MeasurementError(
            f"{label} has invalid derived construction/loss seeds")
    scale = receipt.get("scale")
    if not _is_canonical_scale(scale):
        raise MeasurementError(f"{label} has invalid paired scale")
    for arm in ("trained", "shipped"):
        pmf_digest = (
            receipt[f"{arm}_pmf_sha256"]
            if arm == "trained" else STOCK_PMF_DIGEST)
        staircase_scale = (
            "unset" if arm == "shipped" or scale == -1.0
            else _canonical_staircase_scale_spec(scale))
        _validate_target_receipt(
            receipt[arm]["target_receipt"],
            block_count=block_count,
            block_bytes=block_bytes,
            native=native,
            construction_seed=expected_construction_seed,
            loss_seed=expected_loss_seed,
            loss=measurement_policy["loss"],
            schedule=measurement_policy["schedule"],
            pmf_digest=pmf_digest,
            staircase_scale=staircase_scale,
            label=f"{label}.{arm}.target_receipt",
        )
    trained_target = dict(receipt["trained"]["target_receipt"])
    shipped_target = dict(receipt["shipped"]["target_receipt"])
    for field in ("pmf_sha256", "staircase_scale", "target_mean"):
        trained_target.pop(field)
        shipped_target.pop(field)
    if trained_target != shipped_target:
        raise MeasurementError(
            f"{label} candidate/control changed more than the PMF or scale")
    canonical_goodputs = {}
    for arm in ("trained", "shipped"):
        metrics = receipt[arm]
        expected_goodput = _canonical_goodput(metrics, block_count)
        canonical_goodputs[arm] = expected_goodput
        if receipt[f"{arm}_goodput"] != expected_goodput:
            raise MeasurementError(
                f"{label} has inconsistent {arm} goodput")
    if (receipt["verdict"] == "keep" and
            not canonical_goodputs["trained"] >
            canonical_goodputs["shipped"] *
            (1.0 + receipt["margin_percent"] / 100.0)):
        raise MeasurementError(
            f"{label} keep verdict does not clear its margin")


def _validate_validation_receipt(
        receipt, block_count, construction_seed_base, loss_seed_base, native,
        measurement_policy, search_receipt, label):
    """Replay one paired validation experiment and verify its verdict."""
    required = {
        "protocol", "verdict", "selected_arm", "source_mode", "trials",
        "block_bytes", "scale", "trained_pmf_sha256", "trained_goodput",
        "shipped_goodput", "paired_measurement",
    }
    if not isinstance(receipt, dict) or set(receipt) != required:
        raise MeasurementError(f"{label} is incomplete")
    trials = receipt["trials"]
    block_bytes = receipt["block_bytes"]
    if (receipt["protocol"] != NATIVE_PAIRED_PROTOCOL or
            receipt["verdict"] not in ("keep", "control") or
            receipt["selected_arm"] !=
                ("candidate" if receipt["verdict"] == "keep"
                 else "identity") or
            receipt["source_mode"] != search_receipt["mode"] or
            type(trials) is not int or trials < 4 or trials % 2 != 0 or
            type(block_bytes) is not int or
            not 2 <= block_bytes <= _COMPARE_BLOCK_BYTES_MAX or
            block_bytes % 2 != 0 or
            not _is_canonical_scale(receipt["scale"]) or
            not _same_typed_json(
                receipt["scale"],
                search_receipt["coordinates"]["scale"]) or
            receipt["trained_pmf_sha256"] !=
                search_receipt["peel_pmf_sha256"]):
        raise MeasurementError(f"{label} has invalid paired metadata")
    expected_construction_seed = derive_seed(
        construction_seed_base, "validation", block_count, trials,
        block_bytes, "construction")
    expected_loss_seed = derive_seed(
        loss_seed_base, "validation", block_count, trials, block_bytes, "loss")
    measurement = _validate_paired_measurement_receipt(
        receipt["paired_measurement"],
        block_count=block_count,
        block_bytes=block_bytes,
        candidate_pmf=search_receipt["peel_pmf"],
        degree_scale=(
            None if receipt["scale"] == -1.0 else receipt["scale"]),
        native=native,
        measurement_policy=measurement_policy,
        construction_seed=expected_construction_seed,
        loss_seed=expected_loss_seed,
        label=f"{label}.paired_measurement",
    )
    if search_receipt["mode"] == "shipped":
        require_paired_stock_control(
            measurement, f"{label} shipped-source control")
    should_keep = (
        search_receipt["mode"] != "shipped" and
        measurement.valid_for_promotion
    )
    if (receipt["verdict"] == "keep") != should_keep:
        raise MeasurementError(
            f"{label} verdict contradicts paired promotion policy")
    if (receipt["trained_goodput"] !=
            measurement.candidate.goodput(block_count) or
            receipt["shipped_goodput"] !=
            measurement.identity.goodput(block_count)):
        raise MeasurementError(f"{label} has stale paired goodput")
    for name in ("trained_goodput", "shipped_goodput"):
        if not _is_canonical_nonnegative_float(receipt[name]):
            raise MeasurementError(f"{label} has invalid {name}")
    return measurement


_DIRECT_SETTINGS_FIELDS = {
    "kmin", "kmax", "screen", "refine", "gate_trials",
    "gate_block_bytes", "screen_paired_replicates",
    "paired_replicates", "rank_block_bytes", "rank_top",
    "paired_context", "paired_warmups", "paired_inner_reps",
    "max_overhead", "cache_state", "evict_bytes", "rank_margin",
    "target_profile", "seed_policy", "loss", "schedule",
}
_SWEEP_SETTINGS_FIELDS = {
    "proxy_k_ladder", "paired_replicates", "paired_context",
    "paired_warmups", "paired_inner_reps", "max_overhead", "cache_state",
    "evict_bytes", "rank_margin", "target_profile",
    "seed_policy", "loss", "schedule", "proxy_cost_model",
    "proxy_measure_regime", "proxy_ordering",
    "allow_unverified_cost_model", "search_box",
}
_ENTRY_BASE_FIELDS = {
    "K", "scale", "p1", "tilt", "dmax", "absorb",
    *PAIRED_ARM_FIELDS,
    "goodput", "native", "peel_pmf", "search_receipt",
}
_DIRECT_ENTRY_DIAGNOSTIC_FIELDS = {"seconds", "probes"}
_SWEEP_ENTRY_DIAGNOSTIC_FIELDS = {
    "S", "source_hits", "target_mean", "seconds", "screen",
    "screen_cells", "finals", "rejected", "paired_replicates",
}
_VALIDATION_ENTRY_FIELDS = {
    "validation_receipt", "verified_solve_mbps", "verified_oh",
    "shipped_solve_mbps", "solve_gain_pct",
}


def _settings_match_measurement_policy(settings, measurement_policy):
    return (
        settings.get("target_profile") ==
            measurement_policy["target_profile"] and
        settings.get("seed_policy") == measurement_policy["seed_policy"] and
        type(settings.get("loss")) is float and
        settings["loss"] == measurement_policy["loss"] and
        settings.get("schedule") == measurement_policy["schedule"] and
        measurement_policy.get("loss_model") == PAIRED_LOSS_MODEL
    )


def _validate_search_generator_settings(
        generator, settings, measurement_policy, label):
    """Validate the exact settings schema emitted by each table generator."""
    if not isinstance(settings, dict):
        raise MeasurementError(f"{label} is not an object")
    shared_integer_fields = {
        "paired_replicates", "paired_warmups", "paired_inner_reps",
        "max_overhead", "evict_bytes",
    }

    def valid_paired_settings():
        paired_context = settings.get("paired_context")
        return (
            all(type(settings.get(name)) is int
                for name in shared_integer_fields) and
            settings["paired_replicates"] >= 4 and
            settings["paired_replicates"] % 2 == 0 and
            settings["paired_warmups"] >= 0 and
            settings["paired_warmups"] % 2 == 0 and
            1 <= settings["paired_inner_reps"] <= 1024 and
            0 <= settings["max_overhead"] <= 4096 and
            4096 <= settings["evict_bytes"] <=
                PEELTIMING_EVICT_BYTES_MAX and
            settings.get("cache_state") in ("warm", "cold") and
            (settings["cache_state"] != "cold" or
             settings["paired_inner_reps"] == 1) and
            _is_canonical_absolute_path(paired_context) and
            _is_canonical_nonnegative_float(
                settings.get("rank_margin")) and
            settings["rank_margin"] <= 1.0
        )

    if generator == "tools/peel_direct.py":
        integer_fields = {
            "kmin", "kmax", "screen", "refine", "gate_trials",
            "gate_block_bytes", "screen_paired_replicates",
            "paired_replicates", "rank_block_bytes", "rank_top",
            "paired_warmups", "paired_inner_reps",
            "max_overhead", "evict_bytes",
        }
        if (set(settings) != _DIRECT_SETTINGS_FIELDS or
                any(type(settings.get(name)) is not int
                    for name in integer_fields) or
                not 2 <= settings["kmin"] <= settings["kmax"] <= 64000 or
                settings["screen"] < 1 or settings["refine"] < 0 or
                not 1 <= settings["gate_trials"] <= _COMPARE_TRIALS_MAX or
                not 2 <= settings["gate_block_bytes"] <=
                    _COMPARE_BLOCK_BYTES_MAX or
                settings["gate_block_bytes"] % 2 != 0 or
                not 2 <= settings["rank_block_bytes"] <=
                    _COMPARE_BLOCK_BYTES_MAX or
                settings["rank_block_bytes"] % 2 != 0 or
                settings["screen_paired_replicates"] < 4 or
                settings["screen_paired_replicates"] % 2 != 0 or
                settings["rank_top"] < 1 or
                not valid_paired_settings() or
                not _settings_match_measurement_policy(
                    settings, measurement_policy)):
            raise MeasurementError(
                f"{label} is not an exact direct-search settings receipt")
        try:
            # Direct uses the same native peeltiming invocation policy for
            # its discarded screen-ranking tier and its retained finalists.
            # Replaying only the final receipt would allow a table to claim
            # screen dimensions the generator/native CLI cannot execute.
            for block_count in range(
                    settings["kmin"], settings["kmax"] + 1):
                for replicates in {
                        settings["screen_paired_replicates"],
                        settings["paired_replicates"]}:
                    validate_peeltiming_dimensions(
                        block_count=block_count,
                        block_bytes=settings["rank_block_bytes"],
                        target_profile=settings["target_profile"],
                        seed_policy=settings["seed_policy"],
                        construction_seed=0,
                        loss=settings["loss"],
                        loss_seed=0,
                        schedule=settings["schedule"],
                        warmup_replicates=settings["paired_warmups"],
                        replicates=replicates,
                        inner_reps=settings["paired_inner_reps"],
                        max_overhead=settings["max_overhead"],
                        cache_state=settings["cache_state"],
                        evict_bytes=settings["evict_bytes"],
                        required_margin=settings["rank_margin"],
                    )
        except ValueError as error:
            raise MeasurementError(
                f"{label} has impossible direct peeltiming dimensions: "
                f"{error}")
        return
    if generator != "tools/peel_sweep.py":
        raise MeasurementError(
            f"peel table has unsupported search generator {generator!r}")
    ladder = settings.get("proxy_k_ladder")
    if (set(settings) != _SWEEP_SETTINGS_FIELDS or
            not isinstance(ladder, list) or not ladder or
            any(type(block_count) is not int or
                block_count not in PROXY_K_LADDER
                for block_count in ladder) or
            ladder != sorted(set(ladder)) or
            ladder != [
                block_count for block_count in PROXY_K_LADDER
                if block_count <= ladder[-1]
            ] or
            not valid_paired_settings() or
            settings.get("proxy_cost_model") != PROXY_COST_MODEL or
            settings.get("proxy_ordering") != PROXY_ORDERING_PROTOCOL or
            settings.get("allow_unverified_cost_model") is not True or
            settings.get("search_box") != SEARCH_BOX_PROTOCOL or
            not _settings_match_measurement_policy(
                settings, measurement_policy)):
        raise MeasurementError(
            f"{label} is not an exact proxy-sweep settings receipt")
    _validate_proxy_measure_regime(
        settings.get("proxy_measure_regime"),
        f"{label}.proxy_measure_regime")
    try:
        for block_count in ladder:
            validate_peeltiming_dimensions(
                block_count=block_count,
                block_bytes=4096,
                target_profile=settings["target_profile"],
                seed_policy=settings["seed_policy"],
                construction_seed=0,
                loss=settings["loss"],
                loss_seed=0,
                schedule=settings["schedule"],
                warmup_replicates=settings["paired_warmups"],
                replicates=settings["paired_replicates"],
                inner_reps=settings["paired_inner_reps"],
                max_overhead=settings["max_overhead"],
                cache_state=settings["cache_state"],
                evict_bytes=settings["evict_bytes"],
                required_margin=settings["rank_margin"],
            )
    except ValueError as error:
        raise MeasurementError(
            f"{label} has impossible sweep peeltiming dimensions: {error}")


def _sweep_budget(block_count):
    if block_count <= 1024:
        return 3000, 400, 16, 1000
    if block_count <= 8192:
        return 1200, 250, 12, 1000
    if block_count <= 32768:
        return 500, 150, 8, 600
    return 250, 100, 6, 400


def _validate_entry_generator_settings(
        generator, settings, receipt, block_count, label):
    """Cross-check one selected receipt against its generator invocation."""
    context = receipt["context"]
    if generator == "tools/peel_direct.py":
        expected_context = {
            "screen": settings["screen"],
            "refine": settings["refine"],
            "gate_trials": settings["gate_trials"],
            "gate_block_bytes": settings["gate_block_bytes"],
            "screen_paired_replicates":
                settings["screen_paired_replicates"],
            "rank_top": settings["rank_top"],
        }
        if (receipt["trials"] != settings["paired_replicates"] or
                receipt["block_bytes"] != settings["rank_block_bytes"] or
                any(context.get(name) != value
                    for name, value in expected_context.items())):
            raise MeasurementError(
                f"{label} contradicts its direct-search settings")
        manifest = receipt["paired_measurement"]["manifest"]
        if not (
                manifest.get("warmup_replicates") ==
                    settings["paired_warmups"] and
                manifest.get("replicates") ==
                    settings["paired_replicates"] and
                manifest.get("inner_reps") ==
                    settings["paired_inner_reps"] and
                manifest.get("max_overhead") ==
                    settings["max_overhead"] and
                manifest.get("cache_state") == settings["cache_state"] and
                manifest.get("evict_bytes") == settings["evict_bytes"] and
                manifest.get("required_margin") == settings["rank_margin"]):
            raise MeasurementError(
                f"{label} contradicts its paired measurement settings")
        return
    screen, refine, finals, screen_cells = _sweep_budget(block_count)
    expected_context = {
        "proxy_cost_model": settings["proxy_cost_model"],
        "proxy_measure_regime": settings["proxy_measure_regime"],
        "proxy_ordering": settings["proxy_ordering"],
        "search_box": settings["search_box"],
        "screen": screen,
        "refine": refine,
        "finals": finals,
        "screen_cells": screen_cells,
        "gate_trials": 25,
        "gate_block_bytes": 64,
        "rank_top": 3,
        "threads": 64,
        "batch": 60,
        "cell_base": 900_000_000,
    }
    if (receipt["trials"] != settings["paired_replicates"] or
            receipt["block_bytes"] != 4096 or
            any(context.get(name) != value
                for name, value in expected_context.items())):
        raise MeasurementError(
            f"{label} contradicts its proxy-sweep settings")
    manifest = receipt["paired_measurement"]["manifest"]
    paired_matches = (
        manifest.get("warmup_replicates") == settings["paired_warmups"] and
        manifest.get("replicates") == settings["paired_replicates"] and
        manifest.get("inner_reps") == settings["paired_inner_reps"] and
        manifest.get("max_overhead") == settings["max_overhead"] and
        manifest.get("cache_state") == settings["cache_state"] and
        manifest.get("evict_bytes") == settings["evict_bytes"] and
        manifest.get("required_margin") == settings["rank_margin"]
    )
    if not paired_matches:
        raise MeasurementError(
            f"{label} contradicts its paired measurement settings")


def _validate_entry_field_schema(
        entry, source_generator, native, search_receipt, *,
        validated, selected_shipped, label):
    """Require the exact top-level fields emitted for this entry state."""
    diagnostics = (
        _DIRECT_ENTRY_DIAGNOSTIC_FIELDS
        if source_generator == "tools/peel_direct.py" else
        _SWEEP_ENTRY_DIAGNOSTIC_FIELDS
    )
    expected = set(_ENTRY_BASE_FIELDS)
    if validated:
        expected.update(_VALIDATION_ENTRY_FIELDS)
        if selected_shipped:
            expected.update(
                {
                    "reverted_to_shipped",
                    "search_would_have_solve_delta_pct",
                })
        else:
            expected.update(diagnostics)
    else:
        expected.update(diagnostics)
        if search_receipt["mode"] == "shipped":
            expected.add("reverted_to_shipped")
    if set(entry) != expected:
        differences = sorted(set(entry) ^ expected)
        raise MeasurementError(
            f"{label} has unexpected/missing top-level fields {differences}")
    if ("reverted_to_shipped" in expected and
            entry["reverted_to_shipped"] is not True):
        raise MeasurementError(
            f"{label} has an invalid shipped-selection marker")

    # A validated control deliberately drops generator diagnostics when it
    # reconstructs the selected shipped arm.  All other states retain the
    # diagnostics emitted by the source generator.
    if validated and selected_shipped:
        return
    seconds = entry["seconds"]
    if not _is_canonical_nonnegative_float(seconds):
        raise MeasurementError(f"{label} has invalid elapsed-time diagnostics")
    if source_generator == "tools/peel_direct.py":
        if type(entry["probes"]) is not int or entry["probes"] < 0:
            raise MeasurementError(f"{label} has invalid probe diagnostics")
        return
    context = search_receipt["context"]
    if (type(entry["S"]) is not int or
            entry["S"] != native["staircase"] or
            type(entry["source_hits"]) is not int or
            entry["source_hits"] != native["source_hits"] or
            type(entry["target_mean"]) is not float or
            entry["target_mean"] != native["target_mean"] or
            any(type(entry[name]) is not int
                for name in (
                    "screen", "screen_cells", "finals", "rejected",
                    "paired_replicates")) or
            entry["screen"] != context["screen"] or
            entry["screen_cells"] != context["screen_cells"] or
            entry["finals"] != context["finals"] or
            not 0 <= entry["rejected"] <= 2 * entry["finals"] or
            entry["paired_replicates"] != search_receipt["trials"]):
        raise MeasurementError(f"{label} has invalid proxy-sweep diagnostics")


def _expected_source_k(
        generator, settings, source_selection, validation_settings):
    """Return the complete and selected K sets implied by provenance."""
    if generator == "tools/peel_direct.py":
        complete = list(range(settings["kmin"], settings["kmax"] + 1))
    else:
        complete = list(settings["proxy_k_ladder"])
    if source_selection is None:
        return complete, complete
    selected = [
        block_count for block_count in complete
        if block_count <= validation_settings["kmax"]
    ]
    return complete, selected


def _search_warm_start(entry, generator):
    """Reconstruct the next generator input from a preserved search result."""
    receipt = entry["search_receipt"]
    coordinates = receipt["coordinates"]
    if receipt["mode"] == "shipped":
        return None
    if generator == "tools/peel_direct.py":
        return [
            coordinates["p1"], coordinates["tilt"],
            coordinates["dmax"], coordinates["absorb"],
        ]
    if coordinates["scale"] < 0.0:
        return None
    return [
        int(round(coordinates["scale"] * 100.0)),
        coordinates["p1"], coordinates["tilt"] + 100,
        coordinates["dmax"], coordinates["absorb"],
    ]


def _validate_warm_start_chain(entries, generator, complete_k):
    """Replay the exact Direct ascent or Sweep pivot/up/down dependency."""
    predecessor = {}
    if generator == "tools/peel_direct.py":
        predecessor[complete_k[0]] = None
        predecessor.update({
            block_count: complete_k[index - 1]
            for index, block_count in enumerate(complete_k)
            if index > 0
        })
    else:
        pivot = 128 if 128 in complete_k else complete_k[len(complete_k) // 2]
        up = [block_count for block_count in complete_k if block_count >= pivot]
        down = [
            block_count for block_count in reversed(complete_k)
            if block_count < pivot
        ]
        predecessor[pivot] = None
        for path in (up, [pivot, *down]):
            predecessor.update({
                block_count: path[index - 1]
                for index, block_count in enumerate(path)
                if index > 0
            })
    for block_count, entry in entries.items():
        previous = predecessor[block_count]
        # A validated prefix may omit the original Sweep pivot/upper parent.
        # Once a predecessor is present, every preserved link is still exact.
        if previous is not None and previous not in entries:
            continue
        expected = (
            None if previous is None else
            _search_warm_start(entries[previous], generator)
        )
        if expected is not None and generator == "tools/peel_sweep.py":
            scale_min, scale_max = entry["search_receipt"]["context"][
                "scale_centi"]
            bounds = (
                (scale_min, scale_max), (0, 400), (0, 1700),
                (2, 64), (0, 100),
            )
            expected = [
                max(low, min(high, int(round(value))))
                for value, (low, high) in zip(expected, bounds)
            ]
        actual = entry["search_receipt"]["context"]["warm_start"]
        if not _same_typed_json(actual, expected):
            raise MeasurementError(
                f"peel table entry {block_count} breaks the "
                f"{generator} warm-start chain")


def _validate_table_document_impl(document):
    """Implement strict table validation behind the public numeric boundary."""
    if not isinstance(document, dict):
        raise MeasurementError("peel table root must be an object")
    if document.get("schema") != PEEL_TABLE_SCHEMA:
        raise MeasurementError(
            f"unsupported or legacy peel table schema "
            f"{document.get('schema')!r}; expected {PEEL_TABLE_SCHEMA!r}")
    generator = (
        document.get("provenance", {}).get("generator")
        if isinstance(document.get("provenance"), dict) else None)
    expected_document_fields = {"schema", "provenance", "entries"}
    if generator == "tools/peel_validate.py":
        expected_document_fields.add("source_provenance")
    if set(document) != expected_document_fields:
        raise MeasurementError("peel table has unexpected root fields")
    provenance = document.get("provenance")
    if not isinstance(provenance, dict):
        raise MeasurementError("peel table has no provenance object")
    required = {
        "generator", "benchmark", "native_pmf", "native_compare",
        "native_paired",
        "seed_derivation", "completion_regime", "target_contract",
        "compare_arm", "measurement_policy", "construction_seed_base",
        "loss_seed_base", "settings", "source", "recovery_metric_scope",
        "python_runtime",
    }
    if set(provenance) != required:
        differences = sorted(required ^ set(provenance))
        raise MeasurementError(
            f"peel table provenance has unexpected/missing {differences}")
    if provenance["native_pmf"] != NATIVE_PMF_PROTOCOL:
        raise MeasurementError("peel table used an unsupported PMF source")
    if provenance["native_compare"] != NATIVE_COMPARE_PROTOCOL:
        raise MeasurementError("peel table used an unsupported compare protocol")
    if provenance["native_paired"] != NATIVE_PAIRED_PROTOCOL:
        raise MeasurementError("peel table used an unsupported paired protocol")
    if provenance["seed_derivation"] != SEED_DERIVATION_PROTOCOL:
        raise MeasurementError(
            "peel table used an unsupported seed-derivation protocol")
    completion_regime = provenance["completion_regime"]
    if (not isinstance(completion_regime, dict) or
            set(completion_regime) != {"protocol", "settings"} or
            completion_regime.get("protocol") != COMPLETION_REGIME_PROTOCOL or
            completion_regime.get("settings") !=
            PRODUCTION_COMPLETION_REGIME):
        raise MeasurementError(
            "peel table did not use the production mixed-completion regime")
    if provenance["compare_arm"] != PRODUCTION_COMPARE_ARM:
        raise MeasurementError(
            "peel table did not use the production WH2 compare arm")
    _validate_target_contract_object(
        provenance["target_contract"], "peel table target contract")
    measurement_policy = provenance["measurement_policy"]
    if (not isinstance(measurement_policy, dict) or
            set(measurement_policy) != {
                "target_profile", "seed_policy", "loss", "schedule",
                "loss_model"} or
            measurement_policy.get("target_profile") != TARGET_PROFILE or
            measurement_policy.get("seed_policy") != TARGET_SEED_POLICY or
            measurement_policy.get("loss_model") != PAIRED_LOSS_MODEL or
            not isinstance(measurement_policy.get("schedule"), str) or
            measurement_policy.get("schedule") not in TARGET_SCHEDULES or
            type(measurement_policy.get("loss")) is not float or
            not valid_loss_rate(measurement_policy.get("loss"))):
        raise MeasurementError(
            "peel table has an invalid target measurement policy")
    if provenance["recovery_metric_scope"] != RECOVERY_METRIC_SCOPE:
        raise MeasurementError(
            "peel table has an unsupported recovery-metric scope")
    python_runtime = provenance["python_runtime"]
    if (not isinstance(python_runtime, dict) or
            set(python_runtime) !=
            {"implementation", "version", "cache_tag", "libc"} or
            any(not isinstance(value, str)
                for value in python_runtime.values()) or
            not python_runtime["implementation"] or
            not python_runtime["version"]):
        raise MeasurementError("peel table has invalid Python runtime identity")
    benchmark = provenance["benchmark"]
    source = provenance["source"]
    if (not isinstance(provenance["generator"], str) or
            not provenance["generator"] or
            not isinstance(benchmark, dict) or
            set(benchmark) != {"path", "sha256", "size"} or
            not _is_canonical_absolute_path(benchmark.get("path")) or
            not _is_sha256(benchmark.get("sha256")) or
            type(benchmark.get("size")) is not int or
            benchmark["size"] < 1 or
            not isinstance(source, dict) or
            set(source) != {
                "git_commit", "state_sha256", "file_count",
                "generator_sha256"} or
            not isinstance(source.get("git_commit"), str) or
            len(source["git_commit"]) not in (40, 64) or
            any(character not in "0123456789abcdef"
                for character in source["git_commit"]) or
            not _is_sha256(source.get("state_sha256")) or
            not _is_sha256(source.get("generator_sha256")) or
            type(source.get("file_count")) is not int or
            source["file_count"] < 1 or
            not isinstance(provenance["settings"], dict) or
            any(
                type(provenance[name]) is not int or
                not 0 <= provenance[name] <= 0xffffffffffffffff
                for name in ("construction_seed_base", "loss_seed_base"))):
        raise MeasurementError("peel table has invalid provenance fields")
    expected_search_construction_seed_base = (
        provenance["construction_seed_base"])
    expected_search_loss_seed_base = provenance["loss_seed_base"]
    source_generator = provenance["generator"]
    source_settings = provenance["settings"]
    validation_settings = None
    source_selection = None
    if provenance["generator"] == "tools/peel_validate.py":
        source_provenance = document.get("source_provenance")
        if (not isinstance(source_provenance, dict) or
                set(source_provenance) != {
                    "schema", "document_sha256", "provenance",
                    "entry_count", "selected_entry_count", "selected_K"} or
                source_provenance.get("schema") != PEEL_TABLE_SCHEMA or
                not _is_sha256(source_provenance.get("document_sha256")) or
                not isinstance(source_provenance.get("provenance"), dict) or
                type(source_provenance.get("entry_count")) is not int or
                source_provenance["entry_count"] < 1 or
                type(source_provenance.get("selected_entry_count")) is not
                    int or
                source_provenance["selected_entry_count"] < 1 or
                not isinstance(source_provenance.get("selected_K"), list) or
                not source_provenance["selected_K"] or
                any(type(k) is not int or
                    not 2 <= k <= 64000
                    for k in source_provenance["selected_K"]) or
                source_provenance["selected_K"] !=
                sorted(set(source_provenance["selected_K"]))):
            raise MeasurementError(
                "validated peel table has an invalid source-table receipt")
        source_selection = source_provenance
        source_receipt = source_provenance["provenance"]
        source_benchmark = source_receipt.get("benchmark")
        source_identity = source_receipt.get("source")
        source_runtime = source_receipt.get("python_runtime")
        _validate_target_contract_object(
            source_receipt.get("target_contract"),
            "validated peel table source target contract")
        if (set(source_receipt) != required or
                not isinstance(source_receipt.get("generator"), str) or
                not source_receipt["generator"] or
                not isinstance(source_benchmark, dict) or
                set(source_benchmark) != {"path", "sha256", "size"} or
                not _is_canonical_absolute_path(
                    source_benchmark.get("path")) or
                not _is_sha256(source_benchmark.get("sha256")) or
                type(source_benchmark.get("size")) is not int or
                source_benchmark["size"] < 1 or
                not isinstance(source_identity, dict) or
                set(source_identity) != {
                    "git_commit", "state_sha256", "file_count",
                    "generator_sha256"} or
                not isinstance(source_identity.get("git_commit"), str) or
                len(source_identity["git_commit"]) not in (40, 64) or
                any(character not in "0123456789abcdef"
                    for character in source_identity["git_commit"]) or
                not _is_sha256(source_identity.get("state_sha256")) or
                not _is_sha256(source_identity.get("generator_sha256")) or
                type(source_identity.get("file_count")) is not int or
                source_identity["file_count"] < 1 or
                not isinstance(source_receipt.get("settings"), dict) or
                not isinstance(source_runtime, dict) or
                set(source_runtime) !=
                {"implementation", "version", "cache_tag", "libc"} or
                any(not isinstance(value, str)
                    for value in source_runtime.values()) or
                not source_runtime["implementation"] or
                not source_runtime["version"] or
                source_receipt.get("native_pmf") != NATIVE_PMF_PROTOCOL or
                source_receipt.get("native_compare") !=
                NATIVE_COMPARE_PROTOCOL or
                source_receipt.get("native_paired") !=
                NATIVE_PAIRED_PROTOCOL or
                source_receipt.get("seed_derivation") !=
                SEED_DERIVATION_PROTOCOL or
                source_receipt.get("completion_regime") !=
                provenance["completion_regime"] or
                source_receipt.get("compare_arm") != PRODUCTION_COMPARE_ARM or
                not isinstance(
                    source_receipt.get("measurement_policy"), dict) or
                set(source_receipt["measurement_policy"]) != {
                    "target_profile", "seed_policy", "loss", "schedule",
                    "loss_model"} or
                source_receipt["measurement_policy"].get("target_profile") !=
                TARGET_PROFILE or
                source_receipt["measurement_policy"].get("seed_policy") !=
                TARGET_SEED_POLICY or
                source_receipt["measurement_policy"].get("loss_model") !=
                PAIRED_LOSS_MODEL or
                not isinstance(
                    source_receipt["measurement_policy"].get("schedule"),
                    str) or
                source_receipt["measurement_policy"].get("schedule") not in
                TARGET_SCHEDULES or
                type(source_receipt["measurement_policy"].get("loss")) is not
                    float or
                not valid_loss_rate(
                    source_receipt["measurement_policy"].get("loss")) or
                source_receipt.get("recovery_metric_scope") !=
                RECOVERY_METRIC_SCOPE or
                not isinstance(source_receipt.get("python_runtime"), dict) or
                any(
                    type(source_receipt.get(name)) is not int or
                    not 0 <= source_receipt[name] <= 0xffffffffffffffff
                    for name in (
                        "construction_seed_base", "loss_seed_base"))):
            raise MeasurementError(
                "validated peel table has invalid source provenance")
        if (source_benchmark["sha256"] != benchmark["sha256"] or
                source_benchmark["size"] != benchmark["size"]):
            raise MeasurementError(
                "validated peel table changed benchmark identity")
        validation_settings = provenance["settings"]
        required_validation_settings = {
            "source_table", "source_table_sha256", "source_entry_count",
            "selected_entry_count", "selected_K", "paired_replicates",
            "block_bytes", "kmax", "margin_percent", "paired_context",
            "paired_warmups", "paired_inner_reps", "max_overhead",
            "cache_state", "evict_bytes", "target_profile", "seed_policy",
            "loss", "schedule",
        }
        if (set(validation_settings) != required_validation_settings or
                not _is_canonical_absolute_path(
                    validation_settings.get("source_table")) or
                validation_settings.get("source_table_sha256") !=
                source_provenance["document_sha256"] or
                type(validation_settings.get("source_entry_count")) is not
                    int or
                validation_settings.get("source_entry_count") !=
                source_provenance["entry_count"] or
                type(validation_settings.get("selected_entry_count")) is not
                    int or
                validation_settings.get("selected_entry_count") !=
                source_provenance["selected_entry_count"] or
                not isinstance(
                    validation_settings.get("selected_K"), list) or
                any(type(k) is not int
                    for k in validation_settings["selected_K"]) or
                validation_settings.get("selected_K") !=
                source_provenance["selected_K"] or
                type(validation_settings.get("paired_replicates")) is not
                    int or
                validation_settings["paired_replicates"] < 4 or
                validation_settings["paired_replicates"] % 2 != 0 or
                type(validation_settings.get("block_bytes")) is not int or
                not 2 <= validation_settings["block_bytes"] <=
                _COMPARE_BLOCK_BYTES_MAX or
                validation_settings["block_bytes"] % 2 != 0 or
                type(validation_settings.get("kmax")) is not int or
                validation_settings["kmax"] < 2 or
                not _is_canonical_nonnegative_float(
                    validation_settings.get("margin_percent")) or
                validation_settings["margin_percent"] > 100.0 or
                not _is_canonical_absolute_path(
                    validation_settings.get("paired_context")) or
                type(validation_settings.get("paired_warmups")) is not int or
                validation_settings["paired_warmups"] < 0 or
                validation_settings["paired_warmups"] % 2 != 0 or
                type(validation_settings.get("paired_inner_reps")) is not
                    int or
                not 1 <= validation_settings["paired_inner_reps"] <= 1024 or
                type(validation_settings.get("max_overhead")) is not int or
                not 0 <= validation_settings["max_overhead"] <= 4096 or
                validation_settings.get("cache_state") not in
                    ("warm", "cold") or
                (validation_settings["cache_state"] == "cold" and
                 validation_settings["paired_inner_reps"] != 1) or
                type(validation_settings.get("evict_bytes")) is not int or
                not 4096 <= validation_settings["evict_bytes"] <=
                    PEELTIMING_EVICT_BYTES_MAX or
                type(validation_settings.get("loss")) is not float or
                validation_settings["loss"] !=
                measurement_policy["loss"] or
                validation_settings.get("target_profile") !=
                TARGET_PROFILE or
                validation_settings.get("seed_policy") !=
                TARGET_SEED_POLICY or
                validation_settings.get("schedule") !=
                measurement_policy["schedule"]):
            raise MeasurementError(
                "validated peel table settings contradict its source receipt")
        try:
            for selected_k in validation_settings["selected_K"]:
                validate_peeltiming_dimensions(
                    block_count=selected_k,
                    block_bytes=validation_settings["block_bytes"],
                    target_profile=validation_settings["target_profile"],
                    seed_policy=validation_settings["seed_policy"],
                    construction_seed=0,
                    loss=validation_settings["loss"],
                    loss_seed=0,
                    schedule=validation_settings["schedule"],
                    warmup_replicates=(
                        validation_settings["paired_warmups"]),
                    replicates=validation_settings["paired_replicates"],
                    inner_reps=validation_settings["paired_inner_reps"],
                    max_overhead=validation_settings["max_overhead"],
                    cache_state=validation_settings["cache_state"],
                    evict_bytes=validation_settings["evict_bytes"],
                    required_margin=(
                        validation_settings["margin_percent"] / 100.0),
                )
        except ValueError as error:
            raise MeasurementError(
                "validated peel table settings exceed native peeltiming "
                f"limits: {error}")
        expected_search_construction_seed_base = (
            source_receipt["construction_seed_base"])
        expected_search_loss_seed_base = source_receipt["loss_seed_base"]
        source_generator = source_receipt.get("generator")
        source_settings = source_receipt["settings"]
        source_measurement_policy = source_receipt["measurement_policy"]
    else:
        source_measurement_policy = measurement_policy
    search_protocols = {
        "tools/peel_direct.py": (
            "direct-search", "direct-real-codec"),
        "tools/peel_sweep.py": (
            "funnel-search", "unverified-proxy-funnel"),
    }
    if source_generator not in search_protocols:
        raise MeasurementError(
            f"peel table has unsupported search generator "
            f"{source_generator!r}")
    _validate_search_generator_settings(
        source_generator, source_settings, source_measurement_policy,
        "peel table search settings")
    expected_seed_domain, expected_search_kind = search_protocols[
        source_generator]
    if (source_generator == "tools/peel_sweep.py" and
            source_settings.get("allow_unverified_cost_model") is not True):
        raise MeasurementError(
            "proxy funnel table is missing its unverified-cost opt-in")
    raw_entries = document.get("entries")
    if not isinstance(raw_entries, dict) or not raw_entries:
        raise MeasurementError("peel table entries must be a non-empty object")
    entries = {}
    for key, value in raw_entries.items():
        try:
            block_count = int(key)
        except (TypeError, ValueError):
            raise MeasurementError(f"invalid peel table K key: {key!r}")
        if not 2 <= block_count <= 64000 or str(block_count) != key:
            raise MeasurementError(f"non-canonical peel table K key: {key!r}")
        if (not isinstance(value, dict) or
                type(value.get("K")) is not int or
                value.get("K") != block_count):
            raise MeasurementError(
                f"peel table entry {key} has an invalid K receipt")
        coordinate_names = ("scale", "p1", "tilt", "dmax", "absorb")
        if any(name not in value for name in coordinate_names):
            raise MeasurementError(
                f"peel table entry {key} is missing peel coordinates")
        try:
            scale = value["scale"]
        except (KeyError, TypeError):
            raise MeasurementError(
                f"peel table entry {key} has non-numeric coordinates")
        p1 = value["p1"]
        tilt = value["tilt"]
        dmax = value["dmax"]
        absorb = value["absorb"]
        if (not _is_canonical_scale(value["scale"]) or
                any(type(item) is not int
                    for item in (p1, tilt, dmax, absorb)) or
                not 0 <= p1 <= 400 or
                not -100 <= tilt <= 1600 or
                not 2 <= dmax <= 64 or
                not 0 <= absorb <= 100 or
                ("reverted_to_shipped" in value and
                 not isinstance(value["reverted_to_shipped"], bool)) or
                "reverted_to_control" in value):
            raise MeasurementError(
                f"peel table entry {key} has out-of-domain coordinates")
        if value.get("reverted_to_shipped") and (
                scale != -1.0 or p1 != 100 or tilt != 0 or
                dmax != 64 or absorb != 100):
            raise MeasurementError(
                f"peel table entry {key} has a malformed shipped receipt")
        native = value.get("native")
        _validate_native_receipt(
            native, block_count, f"peel table entry {key}.native")
        selected_pmf_digest = _validate_pmf(
            value.get("peel_pmf"),
            f"peel table entry {key}.peel_pmf")
        search_measurement = _validate_search_receipt(
            value.get("search_receipt"), block_count,
            expected_search_construction_seed_base,
            expected_search_loss_seed_base, expected_seed_domain,
            expected_search_kind, native,
            (source_receipt["measurement_policy"]
             if provenance["generator"] == "tools/peel_validate.py"
             else measurement_policy),
            f"peel table entry {key}.search_receipt")
        _validate_entry_generator_settings(
            source_generator, source_settings, value["search_receipt"],
            block_count, f"peel table entry {key}.search_receipt")
        validation_receipt = value.get("validation_receipt")
        if provenance["generator"] == "tools/peel_validate.py":
            validation_measurement = _validate_validation_receipt(
                validation_receipt, block_count,
                provenance["construction_seed_base"],
                provenance["loss_seed_base"], native, measurement_policy,
                value["search_receipt"],
                f"peel table entry {key}.validation_receipt")
            if (validation_receipt["trials"] !=
                    validation_settings["paired_replicates"] or
                    validation_receipt["block_bytes"] !=
                    validation_settings["block_bytes"] or
                    validation_receipt["scale"] !=
                    value["search_receipt"]["coordinates"]["scale"] or
                    validation_receipt["trained_pmf_sha256"] !=
                    value["search_receipt"]["peel_pmf_sha256"]):
                raise MeasurementError(
                    f"peel table entry {key} validation receipt contradicts "
                    f"its provenance settings")
            validation_manifest = validation_measurement.manifest
            if (validation_manifest["warmup_replicates"] !=
                    validation_settings["paired_warmups"] or
                    validation_manifest["replicates"] !=
                    validation_settings["paired_replicates"] or
                    validation_manifest["inner_reps"] !=
                    validation_settings["paired_inner_reps"] or
                    validation_manifest["max_overhead"] !=
                    validation_settings["max_overhead"] or
                    validation_manifest["cache_state"] !=
                    validation_settings["cache_state"] or
                    validation_manifest["evict_bytes"] !=
                    validation_settings["evict_bytes"] or
                    validation_measurement.required_margin !=
                    validation_settings["margin_percent"] / 100.0):
                raise MeasurementError(
                    f"peel table entry {key} validation measurement "
                    "contradicts its settings")
            selected_shipped = validation_receipt["verdict"] == "control"
            if (selected_shipped !=
                    bool(value.get("reverted_to_shipped"))):
                raise MeasurementError(
                    f"peel table entry {key} contradicts its validation verdict")
            expected_selected_pmf = (
                native["pmf"] if selected_shipped else
                value["search_receipt"]["peel_pmf"])
            selected_arm = (
                validation_measurement.identity if selected_shipped
                else validation_measurement.candidate)
            selected_metrics = selected_arm.as_dict()
            selected_goodput = selected_arm.goodput(block_count)
            selected_trials = validation_receipt["trials"]
            solve_delta = (
                None if validation_measurement.identity.solve_mbps == 0.0 else
                float(round(
                    100.0 * (
                        validation_measurement.candidate.solve_mbps -
                        validation_measurement.identity.solve_mbps) /
                    validation_measurement.identity.solve_mbps,
                    2))
            )
            if solve_delta == 0.0:
                solve_delta = 0.0
            expected_summaries = {
                "verified_solve_mbps": selected_arm.solve_mbps,
                "verified_oh": selected_arm.oh_mean,
                "shipped_solve_mbps":
                    validation_measurement.identity.solve_mbps,
                "solve_gain_pct":
                    0.0 if selected_shipped else solve_delta,
            }
            if selected_shipped:
                expected_summaries[
                    "search_would_have_solve_delta_pct"] = solve_delta
            if (any(
                    name not in value or
                    (expected is None and value[name] is not None) or
                    (expected is not None and (
                        type(value[name]) is not float or
                        not _same_typed_json(value[name], expected)))
                    for name, expected in expected_summaries.items()) or
                    (not selected_shipped and
                     "search_would_have_solve_delta_pct" in value)):
                raise MeasurementError(
                    f"peel table entry {key} has stale validation summary")
        elif validation_receipt is not None:
            raise MeasurementError(
                f"peel table entry {key} has an unexpected validation receipt")
        else:
            search_shipped = value["search_receipt"]["mode"] == "shipped"
            if search_shipped != bool(value.get("reverted_to_shipped")):
                raise MeasurementError(
                    f"peel table entry {key} contradicts its search receipt")
            expected_selected_pmf = value["search_receipt"]["peel_pmf"]
            selected_metrics = (
                search_measurement.identity if search_shipped
                else search_measurement.candidate).as_dict()
            selected_goodput = (
                search_measurement.identity if search_shipped
                else search_measurement.candidate).goodput(block_count)
            selected_trials = value["search_receipt"]["trials"]
        _validate_entry_field_schema(
            value, source_generator, native, value["search_receipt"],
            validated=provenance["generator"] == "tools/peel_validate.py",
            selected_shipped=(
                selected_shipped
                if provenance["generator"] == "tools/peel_validate.py"
                else search_shipped),
            label=f"peel table entry {key}",
        )
        _validate_paired_arm_summary(
            value, selected_trials, f"peel table entry {key} top-level")
        if not _is_canonical_nonnegative_float(value.get("goodput")):
            raise MeasurementError(
                f"peel table entry {key} has invalid top-level goodput")
        if (value["peel_pmf"] != expected_selected_pmf or
                selected_pmf_digest != pmf_sha256(expected_selected_pmf)):
            raise MeasurementError(
                f"peel table entry {key} selected PMF is inconsistent")
        top_coordinates = {name: value[name] for name in coordinate_names}
        expected_top_coordinates = value["search_receipt"]["coordinates"]
        if (provenance["generator"] == "tools/peel_validate.py" and
                selected_shipped):
            expected_top_coordinates = {
                "scale": -1.0,
                "p1": 100,
                "tilt": 0,
                "dmax": 64,
                "absorb": 100,
            }
        if top_coordinates != expected_top_coordinates:
            raise MeasurementError(
                f"peel table entry {key} coordinates contradict its search")
        metric_names = (
            "trials", "fail", "oh_mean", "OH_sd", "OH50", "OH95",
            "OH99", "OH_max", "solve_mbps",
        )
        if (any(not _same_typed_json(
                    value.get(name), selected_metrics[name])
                for name in metric_names) or
                not _same_typed_json(
                    value.get("goodput"), selected_goodput)):
            raise MeasurementError(
                f"peel table entry {key} top-level metrics are stale")
        entries[block_count] = value
    complete_k, expected_selected_k = _expected_source_k(
        source_generator, source_settings, source_selection,
        validation_settings)
    selected_k = sorted(entries)
    if selected_k != expected_selected_k:
        raise MeasurementError(
            "peel table K coverage contradicts its generator settings")
    _validate_warm_start_chain(
        entries, source_generator, complete_k)
    if source_selection is not None and (
            source_selection["entry_count"] != len(complete_k) or
            source_selection["selected_entry_count"] !=
                len(expected_selected_k) or
            source_selection["selected_K"] != expected_selected_k):
        raise MeasurementError(
            "validated peel table source selection contradicts its entries")
    return entries


def validate_table_document(document):
    """Validate a current table, normalizing numeric overflow to refusal."""
    try:
        return _validate_table_document_impl(document)
    except (OverflowError, ValueError) as error:
        # JSON integer tokens have arbitrary precision in Python.  Every
        # schema numeric eventually enters binary64 or a fixed native integer
        # domain; an integer too large for binary64 conversion raises
        # OverflowError, while one too large for Python's guarded decimal
        # conversion can raise ValueError while deriving a receipt seed.
        # Both are malformed input, not unexpected validator crashes.
        raise MeasurementError(
            f"peel table contains an out-of-range numeric value: {error}")


def _strict_json_object(pairs):
    """Reject duplicate object keys instead of silently accepting the last."""
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON object key {key!r}")
        result[key] = value
    return result


def _reject_json_constant(value):
    raise ValueError(f"non-standard JSON constant {value!r}")


def _strict_json_float(value):
    parsed = float(value)
    if not math.isfinite(parsed):
        raise ValueError(f"non-finite JSON number {value!r}")
    return parsed


def strict_json_loads(payload):
    """Decode strict UTF-8 JSON with unique keys and finite numbers."""
    if isinstance(payload, bytes):
        payload = payload.decode("utf-8")
    if not isinstance(payload, str):
        raise ValueError("JSON payload must be text or bytes")
    return json.loads(
        payload,
        object_pairs_hook=_strict_json_object,
        parse_float=_strict_json_float,
        parse_constant=_reject_json_constant)


def read_table_document_snapshot(path):
    """Read, hash, parse, and validate one stable regular-file snapshot."""
    try:
        with open(path, "rb") as source:
            before = os.fstat(source.fileno())
            if not stat.S_ISREG(before.st_mode):
                raise MeasurementError(
                    f"peel table is not a regular file: {path!r}")
            payload = source.read()
            after = os.fstat(source.fileno())
    except OSError as error:
        raise MeasurementError(f"could not read peel table {path!r}: {error}")
    stable_fields = (
        "st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns",
        "st_ctime_ns")
    if (not stat.S_ISREG(before.st_mode) or
            any(getattr(before, name) != getattr(after, name)
                for name in stable_fields) or
            len(payload) != after.st_size):
        raise MeasurementError(
            f"peel table changed while it was being read: {path!r}")
    try:
        document = strict_json_loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise MeasurementError(f"invalid peel table JSON: {error}")
    entries = validate_table_document(document)
    return document, entries, hashlib.sha256(payload).hexdigest()


def read_table_document(path):
    """Read a current peel table, refusing legacy unversioned results."""
    document, entries, _ = read_table_document_snapshot(path)
    return document, entries


def write_json_atomic(path, value):
    """Replace a JSON result only after the complete value is available."""
    directory = os.path.dirname(os.path.abspath(path))
    handle = tempfile.NamedTemporaryFile(
        mode="w", dir=directory, prefix=".peel-table-",
        suffix=".tmp", delete=False)
    try:
        with handle:
            json.dump(value, handle, indent=1, allow_nan=False)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(handle.name, path)
    except BaseException:
        try:
            os.unlink(handle.name)
        except FileNotFoundError:
            pass
        raise
