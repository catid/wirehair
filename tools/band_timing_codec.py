#!/usr/bin/env python3
"""Strict consumer for the native WH2 band-architecture timing protocol.

The native benchmark owns measurement and this module independently owns
validation and aggregation.  In particular, it does not trust native summary
rows, silently repair weak construction seeds, or reuse ``PairedMeasurement``:
the three timing scopes have deliberately different prefix semantics.
"""

import copy
import csv
import hashlib
import math
import os
from dataclasses import dataclass
from typing import Optional

import peel_codec


MeasurementError = peel_codec.MeasurementError

BANDTIMING_PROTOCOL_V1 = "wirehair-v2-bench:bandtiming:dispatch-v1:v1"
BANDTIMING_SCHEMA_V1 = "wirehair.wh2.bandtiming.v1"
BANDTIMING_PROTOCOL_V2 = "wirehair-v2-bench:bandtiming:dispatch-v1:v2"
BANDTIMING_SCHEMA_V2 = "wirehair.wh2.bandtiming.v2"
BANDTIMING_PROTOCOL = BANDTIMING_PROTOCOL_V2
BANDTIMING_SCHEMA = BANDTIMING_SCHEMA_V2
BANDTIMING_DIRECT_SCOPE_V1 = (
    "candidate-dispatch-pair-local-fixed-prefix-solve-v1"
)
BANDTIMING_DIRECT_SCOPE_V2 = (
    "candidate-dispatch-pair-local-fixed-prefix-encoder-"
    "intermediate-witnessed-solve-v2"
)
_BANDTIMING_SCHEMA_BY_PROTOCOL = {
    BANDTIMING_PROTOCOL_V1: BANDTIMING_SCHEMA_V1,
    BANDTIMING_PROTOCOL_V2: BANDTIMING_SCHEMA_V2,
}
_BANDTIMING_DIRECT_SCOPE_BY_SCHEMA = {
    BANDTIMING_SCHEMA_V1: BANDTIMING_DIRECT_SCOPE_V1,
    BANDTIMING_SCHEMA_V2: BANDTIMING_DIRECT_SCOPE_V2,
}
BANDTIMING_ORDER = "ABBABAAB"
BANDTIMING_LABEL_SWAP = "alternating"
BANDTIMING_SCOPE_ORDER = ("encoder", "decoder", "direct")
BANDTIMING_SCHEDULES = frozenset({
    "iid", "burst", "permutation", "systematic-first", "repair-only",
    "adversarial",
})
BANDTIMING_CONSTRUCTION_SEED_DERIVATION = (
    "base_xor_d192ed03_times_rep_plus_1_mod2^32_v1"
)
BANDTIMING_LOSS_SEED_DERIVATION = (
    "base_xor_9e3779b97f4a7c15_times_rep_plus_1_v1"
)
BANDTIMING_SEMANTIC_SEED_DERIVATION = (
    "base-plus-attempt-mod2^32-skip-measured-alias-and-shared-weak-v1"
)
BANDTIMING_UNCERTAINTY = "paired-log-ratio-t95/v1"
BANDTIMING_CLOCK_DOMAIN = "posix-clock-monotonic-ns-v1"
BANDTIMING_CONTEXT_SCHEMA = peel_codec.PAIRED_CONTEXT_SCHEMA
BANDTIMING_THERMAL_SCHEMA = peel_codec.PAIRED_THERMAL_SCHEMA
BANDTIMING_EVICT_BYTES_MAX = peel_codec.PEELTIMING_EVICT_BYTES_MAX
BANDTIMING_WORKING_SET_BYTES_MAX = (
    peel_codec.PEELTIMING_WORKING_SET_BYTES_MAX
)

# Every cross-codec contrast has its own panel.  A/A panels are separate and
# feed scope/arm-specific noise floors; they are not extra architecture
# contrasts.
BANDTIMING_PANEL_SPECS = (
    ("encoder", "encoder_candidate_dispatch", "candidate", "dispatch"),
    ("encoder", "encoder_candidate_wh1", "candidate", "wh1"),
    ("encoder", "encoder_dispatch_wh1", "dispatch", "wh1"),
    ("encoder", "encoder_candidate_aa", "candidate_a", "candidate_b"),
    ("encoder", "encoder_dispatch_aa", "dispatch_a", "dispatch_b"),
    ("encoder", "encoder_wh1_aa", "wh1_a", "wh1_b"),
    ("decoder", "decoder_candidate_dispatch", "candidate", "dispatch"),
    ("decoder", "decoder_candidate_wh1", "candidate", "wh1"),
    ("decoder", "decoder_dispatch_wh1", "dispatch", "wh1"),
    ("decoder", "decoder_candidate_aa", "candidate_a", "candidate_b"),
    ("decoder", "decoder_dispatch_aa", "dispatch_a", "dispatch_b"),
    ("decoder", "decoder_wh1_aa", "wh1_a", "wh1_b"),
    ("direct", "direct_candidate_dispatch", "candidate", "dispatch"),
    ("direct", "direct_candidate_aa", "candidate_a", "candidate_b"),
    ("direct", "direct_dispatch_aa", "dispatch_a", "dispatch_b"),
)
BANDTIMING_PANELS = tuple(spec[1] for spec in BANDTIMING_PANEL_SPECS)
BANDTIMING_CROSS_PANELS = tuple(
    spec[1] for spec in BANDTIMING_PANEL_SPECS
    if not spec[1].endswith("_aa")
)
BANDTIMING_AA_PANELS = tuple(
    spec[1] for spec in BANDTIMING_PANEL_SPECS
    if spec[1].endswith("_aa")
)
assert len(BANDTIMING_PANELS) == 15
assert len(BANDTIMING_CROSS_PANELS) == 7
assert len(BANDTIMING_AA_PANELS) == 8

_PANEL_BY_NAME = {
    name: (scope, left, right, index)
    for index, (scope, name, left, right)
    in enumerate(BANDTIMING_PANEL_SPECS)
}
_ARM_SCOPE_ORDER = (
    ("encoder", "candidate"), ("encoder", "dispatch"), ("encoder", "wh1"),
    ("decoder", "candidate"), ("decoder", "dispatch"), ("decoder", "wh1"),
    ("direct", "candidate"), ("direct", "dispatch"),
)
_ARM_SCOPE_INDEX = {
    key: index for index, key in enumerate(_ARM_SCOPE_ORDER)
}
_REPLICATE_CONSTRUCTION_FIELDS = (
    "candidate_construct_result", "candidate_construct_class",
    "dispatch_construct_result", "dispatch_construct_class",
    "wh1_construct_result", "wh1_construct_class",
)
_LEGACY_REPLICATE_FIELDS = frozenset({
    "replicate", "construction_seed", "loss_seed", "trace_sha256",
    "encoder_candidate_result_class", "encoder_candidate_overhead",
    "encoder_dispatch_result_class", "encoder_dispatch_overhead",
    "encoder_wh1_result_class", "encoder_wh1_overhead",
    "decoder_candidate_result_class", "decoder_candidate_overhead",
    "decoder_dispatch_result_class", "decoder_dispatch_overhead",
    "decoder_wh1_result_class", "decoder_wh1_overhead",
    "direct_candidate_result_class", "direct_candidate_overhead",
    "direct_dispatch_result_class", "direct_dispatch_overhead",
    "censored_panels", "fault_contaminated_panels",
})
_PRIMARY_PANEL = {
    ("encoder", "candidate"): "encoder_candidate_dispatch",
    ("encoder", "dispatch"): "encoder_candidate_dispatch",
    ("encoder", "wh1"): "encoder_candidate_wh1",
    ("decoder", "candidate"): "decoder_candidate_dispatch",
    ("decoder", "dispatch"): "decoder_candidate_dispatch",
    ("decoder", "wh1"): "decoder_candidate_wh1",
    ("direct", "candidate"): "direct_candidate_dispatch",
    ("direct", "dispatch"): "direct_candidate_dispatch",
}
_RESULT_CLASSES = frozenset({
    "success", "need_more", "weak", "error", "not_applicable",
})
_RECOVERY_CLASSES = _RESULT_CLASSES
_SOLVE_COUNT_FIELDS = (
    "inactivated", "binary_def", "heavy_gain", "block_xors",
    "block_muladds",
)
_SOLVE_TIME_FIELDS = (
    "build_ns_sum", "peel_ns_sum", "project_ns_sum",
    "residual_ns_sum", "backsub_ns_sum",
)
_SOLVE_STATS_FIELDS = _SOLVE_COUNT_FIELDS + _SOLVE_TIME_FIELDS
_BYTE_FIELDS = ("source_bytes", "packet_payload_bytes", "intermediate_bytes")


BANDTIMING_MANIFEST_FIELDS = (
    "schema", "dispatch_profile", "seed_policy", "contract_id",
    "K", "bb", "message_bytes", "completion",
    "candidate_S", "candidate_D2", "candidate_gf256", "candidate_gf16",
    "candidate_P", "candidate_x_geometry",
    "candidate_descriptor_sha256",
    "dispatch_S", "dispatch_D2", "dispatch_gf256", "dispatch_gf16",
    "dispatch_P", "dispatch_x_geometry",
    "dispatch_descriptor_sha256", "descriptor_encoding",
    "recovery_mix_count",
    "construction_seed_base", "construction_seed_derivation",
    "wh2_seed_mapping", "wh1_seed_mapping", "wh1_dense_count",
    "semantic_seed_derivation",
    "loss", "loss_seed_base", "loss_seed_derivation", "message_seed_policy",
    "schedule", "loss_model", "trace_encoding",
    "panels", "scope_order", "warmup_replicates", "replicates",
    "slots_per_panel", "panels_per_replicate", "order", "label_swap",
    "inner_reps", "max_overhead",
    "cache_state", "systematic_cache", "wh2_source_cache",
    "wh2_received_systematic_cache", "wh1_encoder_source_policy",
    "wh1_decoder_systematic_policy", "wh1_intermediate_policy",
    "evict_bytes",
    "payload_alignment", "prefault", "cpu_affinity_policy",
    "encoder_scope", "decoder_scope", "direct_scope", "weak_seed_policy",
    "hook_path", "codec_reuse", "context_sha256", "uncertainty",
    "required_margin", "margin_rule", "clock_domain",
    "stream_hash_scope", "started_monotonic_ns", "expected_rows",
)

BANDTIMING_SEMANTIC_FIELDS_V1 = (
    "timed", "construction_seed", "seed_attempt", "seed_attempt_cap",
    "canonical_S", "canonical_D2", "canonical_gf256", "canonical_gf16",
    "canonical_P", "canonical_x", "trace_sha256", "message_sha256",
    "explicit_construct_result", "dispatch_construct_result",
    "explicit_params_sha256", "dispatch_params_sha256", "params_equal",
    "explicit_coefficients_sha256", "dispatch_coefficients_sha256",
    "coefficients_equal",
    "explicit_packet_rows_sha256", "dispatch_packet_rows_sha256",
    "packet_rows_equal",
    "explicit_intermediate_sha256", "dispatch_intermediate_sha256",
    "explicit_intermediate_bytes", "dispatch_intermediate_bytes",
    "intermediate_equal",
    "explicit_payload_sha256", "dispatch_payload_sha256", "payload_equal",
    "explicit_direct_result", "dispatch_direct_result",
    "explicit_solve_sha256", "dispatch_solve_sha256", "direct_equal",
    "explicit_decode_result", "dispatch_decode_result",
    "explicit_overhead", "dispatch_overhead",
    "explicit_recovered_sha256", "dispatch_recovered_sha256",
    "recovery_equal", "message_equal", "pass",
)
BANDTIMING_SEMANTIC_FIELDS_V2 = BANDTIMING_SEMANTIC_FIELDS_V1 + (
    "explicit_direct_intermediate_sha256",
    "dispatch_direct_intermediate_sha256",
)
BANDTIMING_SEMANTIC_FIELDS = BANDTIMING_SEMANTIC_FIELDS_V2

BANDTIMING_COLUMNS_V1 = (
    "replicate", "measured", "scope", "panel", "panel_index", "slot",
    "pair", "label", "role", "label_swap", "construction_seed",
    "loss_seed", "trace_sha256", "construct_result", "result",
    "result_class", "recovery_result", "recovery_class",
    "recovery_ok",
    "encoded_symbols", "received_symbols", "arm_overhead",
    "fixed_prefix_symbols", "timing_eligible", "panel_censored",
    "censor_reason", "preflight_result", "timing_result", "outcome_stable",
    "elapsed_ns", "inner_reps", "saturated", "cpu_before", "cpu_after",
    "cpu_migrated", "minflt_delta", "majflt_delta", "fault_contaminated",
    "stats_available", *_SOLVE_STATS_FIELDS, *_BYTE_FIELDS,
)
BANDTIMING_COLUMNS_V2 = BANDTIMING_COLUMNS_V1 + (
    "intermediate_sha256",
)
BANDTIMING_COLUMNS = BANDTIMING_COLUMNS_V2


def _bandtiming_format(protocol):
    if protocol == BANDTIMING_PROTOCOL_V1:
        return (
            BANDTIMING_SCHEMA_V1,
            BANDTIMING_SEMANTIC_FIELDS_V1,
            BANDTIMING_COLUMNS_V1,
            False,
        )
    if protocol == BANDTIMING_PROTOCOL_V2:
        return (
            BANDTIMING_SCHEMA_V2,
            BANDTIMING_SEMANTIC_FIELDS_V2,
            BANDTIMING_COLUMNS_V2,
            True,
        )
    raise MeasurementError("bandtiming expected protocol is unsupported")


@dataclass(frozen=True)
class BandDescriptor:
    """One explicit mixed-completion equation descriptor."""

    staircase: int
    dense_rows: int
    gf256_rows: int
    gf16_rows: int
    period: int
    x_mode: str

    @property
    def heavy_rows(self):
        return self.gf256_rows + self.gf16_rows

    def as_dict(self):
        return {
            "S": self.staircase,
            "D2": self.dense_rows,
            "gf256": self.gf256_rows,
            "gf16": self.gf16_rows,
            "P": self.period,
            "x": self.x_mode,
        }


@dataclass(frozen=True)
class BandArmSummary:
    """Independently aggregated result for one arm in one timing scope."""

    scope: str
    arm: str
    trials: int
    successes: int
    need_more: int
    weak: int
    errors: int
    overhead_mean: float
    overhead_sd: float
    overhead50: float
    overhead95: float
    overhead99: float
    overhead_max: float
    timed_slots: int
    elapsed_ns: int
    throughput_mbps: float
    aa_eligible_replicates: int
    aa_log_cost_mean: float
    aa_log_cost_ci_low: Optional[float]
    aa_log_cost_ci_high: Optional[float]
    aa_floor_log: float

    def as_dict(self):
        return {
            "scope": self.scope,
            "arm": self.arm,
            "trials": self.trials,
            "successes": self.successes,
            "need_more": self.need_more,
            "weak": self.weak,
            "errors": self.errors,
            "overhead_mean": self.overhead_mean,
            "overhead_sd": self.overhead_sd,
            "overhead50": self.overhead50,
            "overhead95": self.overhead95,
            "overhead99": self.overhead99,
            "overhead_max": self.overhead_max,
            "timed_slots": self.timed_slots,
            "elapsed_ns": self.elapsed_ns,
            "throughput_mbps": self.throughput_mbps,
            "aa_eligible_replicates": self.aa_eligible_replicates,
            "aa_log_cost_mean": self.aa_log_cost_mean,
            "aa_log_cost_ci_low": self.aa_log_cost_ci_low,
            "aa_log_cost_ci_high": self.aa_log_cost_ci_high,
            "aa_floor_log": self.aa_floor_log,
        }


@dataclass(frozen=True)
class BandTimingContrast:
    """One of the seven real cross-architecture timing contrasts."""

    name: str
    scope: str
    left_arm: str
    right_arm: str
    eligible_replicates: int
    log_cost_mean: float
    log_cost_ci_low: Optional[float]
    log_cost_ci_high: Optional[float]
    left_aa_floor_log: float
    right_aa_floor_log: float
    effective_floor_log: float
    required_margin: float
    recovery_regressions: int
    recovery_improvements: int
    both_failures: int
    timing_ci_available: bool
    left_faster: bool

    def as_dict(self):
        return {
            "name": self.name,
            "scope": self.scope,
            "left_arm": self.left_arm,
            "right_arm": self.right_arm,
            "eligible_replicates": self.eligible_replicates,
            "log_cost_mean": self.log_cost_mean,
            "log_cost_ci_low": self.log_cost_ci_low,
            "log_cost_ci_high": self.log_cost_ci_high,
            "left_aa_floor_log": self.left_aa_floor_log,
            "right_aa_floor_log": self.right_aa_floor_log,
            "effective_floor_log": self.effective_floor_log,
            "required_margin": self.required_margin,
            "recovery_regressions": self.recovery_regressions,
            "recovery_improvements": self.recovery_improvements,
            "both_failures": self.both_failures,
            "timing_ci_available": self.timing_ci_available,
            "left_faster": self.left_faster,
        }


@dataclass(frozen=True)
class BandTimingMeasurement:
    """Complete independently checked bandtiming experiment."""

    protocol: str
    manifest: dict
    semantic: dict
    context: dict
    candidate_descriptor: BandDescriptor
    dispatch_descriptor: BandDescriptor
    arms: tuple
    contrasts: tuple
    replicate_receipts: tuple
    weak_cells: tuple
    stream_sha256: str
    finished_monotonic_ns: int
    final_context_sha256: str
    native_stdout: str

    @property
    def valid_for_promotion(self):
        return self.protocol == BANDTIMING_PROTOCOL_V2

    def as_dict(self):
        return {
            "protocol": self.protocol,
            "manifest": dict(self.manifest),
            "semantic": dict(self.semantic),
            "context": copy.deepcopy(self.context),
            "candidate_descriptor": self.candidate_descriptor.as_dict(),
            "dispatch_descriptor": self.dispatch_descriptor.as_dict(),
            "arms": [arm.as_dict() for arm in self.arms],
            "contrasts": [contrast.as_dict() for contrast in self.contrasts],
            "replicates": copy.deepcopy(list(self.replicate_receipts)),
            "weak_cells": copy.deepcopy(list(self.weak_cells)),
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


@dataclass(frozen=True)
class _ReplicateRowIndex:
    """Append-only buckets populated during validation and then read-only."""

    panel_rows: tuple
    scope_arm_rows: tuple


def dispatch_band_descriptor(block_count):
    """Return the frozen dispatch-v1 completion descriptor at ``K``."""
    if type(block_count) is not int or not 2 <= block_count <= 100:
        raise ValueError("invalid dispatch band K")
    return BandDescriptor(
        staircase=peel_codec._dispatch_staircase_count(block_count),
        dense_rows=4,
        gf256_rows=10,
        gf16_rows=2,
        period=244,
        x_mode="frozen",
    )


def _base_arm(role):
    if role.endswith("_a") or role.endswith("_b"):
        return role[:-2]
    return role


_ARM_PANEL_INDICES = {
    (scope, arm): tuple(
        panel_index
        for panel_index, (
            panel_scope, unused_panel, left, right,
        ) in enumerate(BANDTIMING_PANEL_SPECS)
        if panel_scope == scope and
        arm in (_base_arm(left), _base_arm(right))
    )
    for scope, arm in _ARM_SCOPE_ORDER
}


def _lossless_expected_tail(schedule, block_count, loss, prefix_symbols):
    """Return a deterministic lossless tail state, or ``None`` if unknown."""
    if loss != 0.0:
        return None
    if schedule in ("iid", "burst"):
        return int(prefix_symbols >= block_count)
    if schedule == "systematic-first" and prefix_symbols >= block_count:
        return 1
    if (schedule == "permutation" and
            prefix_symbols >= block_count + 512):
        # The native schedule shuffles one K+512 candidate chunk for the
        # entire accepted small-K domain, so consuming that chunk guarantees
        # that every systematic id appeared.
        return 1
    return None


def _tail_candidate_deadline(schedule, block_count):
    """Latest output rank at which systematic id K-1 can first appear."""
    if schedule in ("iid", "burst", "systematic-first"):
        return block_count
    if schedule == "permutation":
        return block_count + 512
    return None


def _validate_descriptor(descriptor, block_count, label):
    if not isinstance(descriptor, BandDescriptor):
        raise ValueError(f"{label} must be a BandDescriptor")
    if type(block_count) is not int or not 2 <= block_count <= 100:
        raise ValueError(f"invalid {label} K")
    values = (
        descriptor.staircase, descriptor.dense_rows, descriptor.gf256_rows,
        descriptor.gf16_rows, descriptor.period,
    )
    if (any(type(value) is not int for value in values) or
            not 1 <= descriptor.staircase <= 64000 or
            not 0 <= descriptor.dense_rows <= 64 or
            not 1 <= descriptor.gf256_rows <= 10 or
            not 0 <= descriptor.gf16_rows <= 4 or
            not descriptor.heavy_rows <= descriptor.period <= 244 or
            descriptor.x_mode not in ("frozen", "tracking-x") or
            block_count + descriptor.staircase + descriptor.dense_rows +
                descriptor.heavy_rows > 0xffff):
        raise ValueError(f"invalid {label}")


def validate_bandtiming_dimensions(
        *, block_count, block_bytes, candidate, dispatch_profile,
        seed_policy, construction_seed, loss, loss_seed, schedule,
        warmup_replicates, replicates, inner_reps, max_overhead,
        cache_state, systematic_cache, evict_bytes, required_margin):
    """Mirror the finite request domain of the native bandtiming command."""
    _validate_descriptor(candidate, block_count, "candidate descriptor")
    dispatch = dispatch_band_descriptor(block_count)
    _validate_descriptor(dispatch, block_count, "dispatch descriptor")

    def working_set_bytes():
        candidate_width = (
            block_count + candidate.staircase + candidate.dense_rows +
            candidate.heavy_rows
        )
        dispatch_width = (
            block_count + dispatch.staircase + dispatch.dense_rows +
            dispatch.heavy_rows
        )
        packet_slots = block_count + max_overhead
        blocks = (
            6 * packet_slots + 4 * candidate_width +
            4 * dispatch_width + 8 * block_count + 256
        )
        return blocks * block_bytes + block_count * 16384 + evict_bytes

    try:
        numeric_loss = (
            float(loss)
            if not isinstance(loss, bool) and isinstance(loss, (int, float))
            else math.nan
        )
    except OverflowError:
        numeric_loss = math.nan
    try:
        numeric_margin = (
            float(required_margin)
            if not isinstance(required_margin, bool) and
            isinstance(required_margin, (int, float))
            else math.nan
        )
    except OverflowError:
        numeric_margin = math.nan
    if (dispatch_profile != peel_codec.TARGET_PROFILE or
            seed_policy != peel_codec.TARGET_SEED_POLICY or
            type(block_count) is not int or not 2 <= block_count <= 100 or
            type(block_bytes) is not int or
            not 2 <= block_bytes <= 0x7fffffff or block_bytes % 2 != 0 or
            type(construction_seed) is not int or
            not 0 <= construction_seed <= 0xffffffff or
            type(loss_seed) is not int or
            not 0 <= loss_seed <= 0xffffffffffffffff or
            not math.isfinite(numeric_loss) or
            (numeric_loss == 0.0 and
             math.copysign(1.0, numeric_loss) < 0.0) or
            not 0.0 <= numeric_loss <= 0.99 or
            not isinstance(schedule, str) or
            schedule not in BANDTIMING_SCHEDULES or
            type(warmup_replicates) is not int or
            warmup_replicates < 0 or
            type(replicates) is not int or
            replicates < 3 or
            warmup_replicates + replicates > 10000 or
            type(inner_reps) is not int or not 1 <= inner_reps <= 1024 or
            type(max_overhead) is not int or
            not 0 <= max_overhead <= min(4096, block_count + 512) or
            cache_state not in ("warm", "cold") or
            systematic_cache not in ("off", "on") or
            type(evict_bytes) is not int or
            not 4096 <= evict_bytes <= BANDTIMING_EVICT_BYTES_MAX or
            working_set_bytes() > BANDTIMING_WORKING_SET_BYTES_MAX or
            (cache_state == "cold" and inner_reps != 1) or
            (warmup_replicates + replicates) * 120 * inner_reps *
                (block_count + max_overhead) > 1_000_000_000_000 or
            not math.isfinite(numeric_margin) or
            (numeric_margin == 0.0 and
             math.copysign(1.0, numeric_margin) < 0.0) or
            not 0.0 <= numeric_margin <= 1.0):
        raise ValueError("invalid bandtiming dimensions or policy")


def _expected_seed(base, replicate):
    return (
        base ^ (((replicate + 1) * 0xd192ed03) & 0xffffffff)
    ) & 0xffffffff


def _expected_loss_seed(base, replicate):
    return (
        base ^
        (((replicate + 1) * 0x9e3779b97f4a7c15) & 0xffffffffffffffff)
    ) & 0xffffffffffffffff


def _descriptor_spec(descriptor, block_count):
    _validate_descriptor(descriptor, block_count, "band descriptor")
    return (
        f"S={descriptor.staircase}\n"
        f"D2={descriptor.dense_rows}\n"
        f"gf256={descriptor.gf256_rows}\n"
        f"gf16={descriptor.gf16_rows}\n"
        f"P={descriptor.period}\n"
        f"x={descriptor.x_mode}\n"
    )


def descriptor_sha256(descriptor, block_count):
    """Hash the protocol's canonical architecture descriptor spelling."""
    return hashlib.sha256(
        _descriptor_spec(descriptor, block_count).encode("ascii")
    ).hexdigest()


def _expected_role(panel, replicate, label):
    try:
        unused_scope, left, right, unused_index = _PANEL_BY_NAME[panel]
    except KeyError:
        raise MeasurementError(f"unknown bandtiming panel {panel!r}")
    if label not in ("A", "B"):
        raise MeasurementError(f"invalid bandtiming label {label!r}")
    # Native deliberately reverses the logical labels on every odd replicate.
    if (replicate & 1) != 0:
        left, right = right, left
    return left if label == "A" else right


def _percentile_integer(values, probability):
    if not values:
        return 0.0
    ordered = sorted(values)
    index = max(
        0, min(len(ordered) - 1, math.ceil(probability * len(ordered)) - 1))
    return float(ordered[index])


def _parse_row(fields, label, columns):
    if len(fields) != len(columns):
        raise MeasurementError(f"{label} has the wrong column count")
    raw = dict(zip(columns, fields))
    uint32_names = {
        "replicate", "measured", "panel_index", "slot", "pair",
        "label_swap", "construction_seed", "recovery_ok",
        "encoded_symbols", "received_symbols", "timing_eligible",
        "panel_censored", "inner_reps", "saturated", "stats_available",
    }
    uint64_names = {"loss_seed", "elapsed_ns"}
    int32_names = {
        "construct_result", "result", "recovery_result", "arm_overhead",
        "fixed_prefix_symbols", "preflight_result", "timing_result",
        "outcome_stable", "cpu_before", "cpu_after", "cpu_migrated",
        "fault_contaminated",
    }
    int64_names = {
        "minflt_delta", "majflt_delta", *_SOLVE_STATS_FIELDS,
        *_BYTE_FIELDS,
    }
    row = {}
    for name in uint32_names:
        row[name] = peel_codec._parse_decimal_integer(
            raw[name], f"{label}.{name}", minimum=0, maximum=0xffffffff)
    for name in uint64_names:
        row[name] = peel_codec._parse_decimal_integer(
            raw[name], f"{label}.{name}",
            minimum=0, maximum=0xffffffffffffffff)
    for name in int32_names:
        row[name] = peel_codec._parse_decimal_integer(
            raw[name], f"{label}.{name}",
            minimum=-0x80000000, maximum=0x7fffffff)
    for name in int64_names:
        row[name] = peel_codec._parse_decimal_integer(
            raw[name], f"{label}.{name}",
            minimum=-1, maximum=0x7fffffffffffffff)
    for name in (
            "scope", "panel", "label", "role", "result_class",
            "recovery_class", "censor_reason"):
        row[name] = raw[name]
    if not peel_codec._is_sha256(raw["trace_sha256"]):
        raise MeasurementError(f"{label}.trace_sha256 is invalid")
    row["trace_sha256"] = raw["trace_sha256"]
    row["intermediate_sha256"] = raw.get(
        "intermediate_sha256", "not_applicable")
    if row["result_class"] not in _RESULT_CLASSES:
        raise MeasurementError(f"{label}.result_class is invalid")
    if row["recovery_class"] not in _RECOVERY_CLASSES:
        raise MeasurementError(f"{label}.recovery_class is invalid")
    return row


def _parse_manifest(
        raw, *, block_count, block_bytes, candidate, dispatch_profile,
        seed_policy, construction_seed, loss, loss_seed, schedule,
        warmup_replicates, replicates, inner_reps, max_overhead,
        cache_state, systematic_cache, evict_bytes, context_sha256,
        required_margin, expected_rows, expected_schema):
    integer_names = {
        "K", "bb", "message_bytes", "candidate_S", "candidate_D2",
        "candidate_gf256",
        "candidate_gf16", "candidate_P", "dispatch_S", "dispatch_D2",
        "dispatch_gf256", "dispatch_gf16", "dispatch_P",
        "recovery_mix_count", "wh1_dense_count",
        "warmup_replicates", "replicates",
        "slots_per_panel", "panels_per_replicate", "inner_reps",
        "max_overhead", "wh2_source_cache",
        "wh2_received_systematic_cache", "payload_alignment", "prefault",
        "expected_rows",
    }
    manifest = dict(raw)
    for name in integer_names:
        manifest[name] = peel_codec._parse_decimal_integer(
            raw[name], f"bandtiming manifest.{name}", minimum=0,
            maximum=0xffffffff)
    manifest["construction_seed_base"] = peel_codec._parse_decimal_integer(
        raw["construction_seed_base"],
        "bandtiming manifest.construction_seed_base",
        minimum=0, maximum=0xffffffff)
    manifest["loss_seed_base"] = peel_codec._parse_decimal_integer(
        raw["loss_seed_base"], "bandtiming manifest.loss_seed_base",
        minimum=0, maximum=0xffffffffffffffff)
    manifest["evict_bytes"] = peel_codec._parse_decimal_integer(
        raw["evict_bytes"], "bandtiming manifest.evict_bytes",
        minimum=0, maximum=0xffffffffffffffff)
    manifest["started_monotonic_ns"] = peel_codec._parse_decimal_integer(
        raw["started_monotonic_ns"],
        "bandtiming manifest.started_monotonic_ns",
        minimum=1, maximum=0xffffffffffffffff)
    for name in ("loss", "required_margin"):
        manifest[name] = peel_codec._parse_finite_float(
            raw[name], f"bandtiming manifest.{name}", minimum=0.0,
            maximum=1.0)

    dispatch = dispatch_band_descriptor(block_count)
    cache_enabled = int(systematic_cache == "on")
    if expected_schema not in _BANDTIMING_DIRECT_SCOPE_BY_SCHEMA:
        raise MeasurementError("bandtiming expected schema is unsupported")
    expected = {
        "schema": expected_schema,
        "dispatch_profile": dispatch_profile,
        "seed_policy": seed_policy,
        "contract_id": peel_codec.TARGET_CONTRACT["contract_id"],
        "K": block_count,
        "bb": block_bytes,
        "message_bytes": (block_count - 1) * block_bytes + 1,
        "completion": "mixed",
        "candidate_S": candidate.staircase,
        "candidate_D2": candidate.dense_rows,
        "candidate_gf256": candidate.gf256_rows,
        "candidate_gf16": candidate.gf16_rows,
        "candidate_P": candidate.period,
        "candidate_x_geometry": candidate.x_mode,
        "candidate_descriptor_sha256":
            descriptor_sha256(candidate, block_count),
        "dispatch_S": dispatch.staircase,
        "dispatch_D2": dispatch.dense_rows,
        "dispatch_gf256": dispatch.gf256_rows,
        "dispatch_gf16": dispatch.gf16_rows,
        "dispatch_P": dispatch.period,
        "dispatch_x_geometry": dispatch.x_mode,
        "dispatch_descriptor_sha256":
            descriptor_sha256(dispatch, block_count),
        "descriptor_encoding": "S-D2-gf256-gf16-P-x-newline-v1",
        "recovery_mix_count": 3,
        "construction_seed_base": construction_seed,
        "construction_seed_derivation":
            BANDTIMING_CONSTRUCTION_SEED_DERIVATION,
        "wh2_seed_mapping": "zero-extend-u32-precode_xor-fold-packet-v1",
        "wh1_seed_mapping": "low16-peel-high16-dense-v1",
        "wh1_dense_count": peel_codec._legacy_dense_count(block_count),
        "semantic_seed_derivation":
            BANDTIMING_SEMANTIC_SEED_DERIVATION,
        "loss": float(loss),
        "loss_seed_base": loss_seed,
        "loss_seed_derivation": BANDTIMING_LOSS_SEED_DERIVATION,
        "message_seed_policy": "replicate-loss-seed-partial-final-v1",
        "schedule": schedule,
        "loss_model": "packet-schedule-v1",
        "trace_encoding": "wirehair-wh2-bandtiming-loss-trace-v1",
        "panels": "+".join(BANDTIMING_PANELS),
        "scope_order": "+".join(BANDTIMING_SCOPE_ORDER),
        "warmup_replicates": warmup_replicates,
        "replicates": replicates,
        "slots_per_panel": 8,
        "panels_per_replicate": 15,
        "order": BANDTIMING_ORDER,
        "label_swap": BANDTIMING_LABEL_SWAP,
        "inner_reps": inner_reps,
        "max_overhead": max_overhead,
        "cache_state": cache_state,
        "systematic_cache": systematic_cache,
        "wh2_source_cache": cache_enabled,
        "wh2_received_systematic_cache": cache_enabled,
        "wh1_encoder_source_policy":
            "borrow-when-off-copy-when-on-v1",
        "wh1_decoder_systematic_policy":
            "native-staging-uncontrolled-v1",
        "wh1_intermediate_policy": "unavailable-zero-v1",
        "evict_bytes": evict_bytes,
        "payload_alignment": 64,
        "prefault": 1,
        "cpu_affinity_policy": "first-allowed-affinity-v1",
        "encoder_scope": "fresh-object-init-through-first-K-symbols-v1",
        "decoder_scope":
            "fresh-init-outside-timer-first-feed-through-own-success-v1",
        "direct_scope":
            _BANDTIMING_DIRECT_SCOPE_BY_SCHEMA[expected_schema],
        "weak_seed_policy": "panel-local-balanced-censor-v1",
        "hook_path":
            "caller-pinned-explicit-transaction-attempt-zero-v2",
        "codec_reuse": "none-fresh-object-every-inner-v1",
        "context_sha256": context_sha256,
        "uncertainty": BANDTIMING_UNCERTAINTY,
        "required_margin": float(required_margin),
        "margin_rule":
            "upper-log-cost-lt-negative-required-margin-and-arm-aa-floors-v1",
        "clock_domain": BANDTIMING_CLOCK_DOMAIN,
        "stream_hash_scope": "body-plus-done-prefix-v1",
        "expected_rows": expected_rows,
    }
    for name, value in expected.items():
        if manifest[name] != value:
            raise MeasurementError(
                f"bandtiming manifest.{name} does not match the request")
    if raw["loss"] != f"{float(loss):.17g}":
        raise MeasurementError("bandtiming loss spelling is noncanonical")
    if raw["required_margin"] != f"{float(required_margin):.17g}":
        raise MeasurementError(
            "bandtiming required margin spelling is noncanonical")
    return manifest


def _parse_semantic(
        raw, dispatch, block_count, block_bytes, construction_seed_base,
        measured_seeds, measured_traces, max_overhead, witnessed):
    integer_names = {
        "timed", "construction_seed", "seed_attempt", "seed_attempt_cap",
        "canonical_S", "canonical_D2", "canonical_gf256", "canonical_gf16",
        "canonical_P",
        "explicit_construct_result", "dispatch_construct_result",
        "params_equal", "coefficients_equal", "packet_rows_equal",
        "explicit_intermediate_bytes", "dispatch_intermediate_bytes",
        "intermediate_equal", "payload_equal",
        "explicit_direct_result", "dispatch_direct_result", "direct_equal",
        "explicit_decode_result", "dispatch_decode_result",
        "explicit_overhead", "dispatch_overhead", "recovery_equal",
        "message_equal", "pass",
    }
    semantic = dict(raw)
    for name in integer_names:
        semantic[name] = peel_codec._parse_decimal_integer(
            raw[name], f"bandtiming semantic.{name}", minimum=0,
            maximum=(
                0xffffffffffffffff
                if name in (
                    "construction_seed", "explicit_intermediate_bytes",
                    "dispatch_intermediate_bytes")
                else 0xffffffff))
    digest_names = {
        "trace_sha256", "message_sha256",
        "explicit_params_sha256", "dispatch_params_sha256",
        "explicit_coefficients_sha256",
        "dispatch_coefficients_sha256", "explicit_packet_rows_sha256",
        "dispatch_packet_rows_sha256", "explicit_intermediate_sha256",
        "dispatch_intermediate_sha256", "explicit_payload_sha256",
        "dispatch_payload_sha256", "explicit_solve_sha256",
        "dispatch_solve_sha256", "explicit_recovered_sha256",
        "dispatch_recovered_sha256",
    }
    if witnessed:
        digest_names.update({
            "explicit_direct_intermediate_sha256",
            "dispatch_direct_intermediate_sha256",
        })
    if any(not peel_codec._is_sha256(raw[name]) for name in digest_names):
        raise MeasurementError("bandtiming semantic digest is invalid")
    expected_semantic_seed = (
        construction_seed_base + semantic["seed_attempt"]) & 0xffffffff
    direct_intermediates_match = (
        not witnessed or
        (
            semantic["explicit_direct_intermediate_sha256"] ==
                semantic["explicit_intermediate_sha256"] and
            semantic["dispatch_direct_intermediate_sha256"] ==
                semantic["dispatch_intermediate_sha256"]
        )
    )
    if (semantic["timed"] != 0 or
            semantic["construction_seed"] > 0xffffffff or
            semantic["construction_seed"] != expected_semantic_seed or
            semantic["construction_seed"] in measured_seeds or
            semantic["trace_sha256"] in measured_traces or
            semantic["seed_attempt"] >= semantic["seed_attempt_cap"] or
            semantic["seed_attempt_cap"] != 256 or
            semantic["canonical_S"] != dispatch.staircase or
            semantic["canonical_D2"] != dispatch.dense_rows or
            semantic["canonical_gf256"] != dispatch.gf256_rows or
            semantic["canonical_gf16"] != dispatch.gf16_rows or
            semantic["canonical_P"] != dispatch.period or
            semantic["canonical_x"] != dispatch.x_mode or
            any(semantic[name] != 1 for name in (
                "params_equal", "coefficients_equal", "packet_rows_equal",
                "intermediate_equal", "payload_equal",
                "direct_equal", "recovery_equal", "message_equal", "pass")) or
            any(semantic[name] != 0 for name in (
                "explicit_construct_result", "dispatch_construct_result",
                "explicit_direct_result", "dispatch_direct_result",
                "explicit_decode_result", "dispatch_decode_result")) or
            semantic["explicit_params_sha256"] !=
                semantic["dispatch_params_sha256"] or
            semantic["explicit_coefficients_sha256"] !=
                semantic["dispatch_coefficients_sha256"] or
            semantic["explicit_packet_rows_sha256"] !=
                semantic["dispatch_packet_rows_sha256"] or
            semantic["explicit_intermediate_sha256"] !=
                semantic["dispatch_intermediate_sha256"] or
            semantic["explicit_intermediate_bytes"] !=
                semantic["dispatch_intermediate_bytes"] or
            semantic["explicit_intermediate_bytes"] != (
                block_count + dispatch.staircase + dispatch.dense_rows +
                dispatch.heavy_rows) * block_bytes or
            semantic["explicit_payload_sha256"] !=
                semantic["dispatch_payload_sha256"] or
            semantic["explicit_solve_sha256"] !=
                semantic["dispatch_solve_sha256"] or
            not direct_intermediates_match or
            semantic["explicit_overhead"] != semantic["dispatch_overhead"] or
            semantic["explicit_overhead"] > max_overhead or
            semantic["explicit_recovered_sha256"] !=
                semantic["dispatch_recovered_sha256"] or
            semantic["message_sha256"] !=
                semantic["explicit_recovered_sha256"]):
        raise MeasurementError(
            "bandtiming canonical explicit/dispatch semantic bridge failed")
    # The native hashes the canonical descriptor.  The independent consumer
    # binds its dimensions separately through the manifest and does not try to
    # reproduce a private C++ byte encoding here.
    _validate_descriptor(
        dispatch, block_count, "semantic dispatch descriptor")
    return semantic


def _row_succeeded(row):
    if row["construct_result"] != 0 or row["result_class"] != "success":
        return False
    if row["scope"] == "decoder":
        return (
            row["recovery_result"] == 0 and
            row["recovery_class"] == "success" and
            row["recovery_ok"] == 1
        )
    return True


def _generic_result_class(code):
    if code == -1:
        return "not_applicable"
    if code == 0:
        return "success"
    if code in (1, 7):
        return "need_more"
    if code in (3, 4):
        return "weak"
    if code in (2, 5, 6, 8, 9, 10):
        return "error"
    raise MeasurementError(f"unknown Wirehair result code {code}")


def _expected_construct_class(row):
    # Native output normalizes only probe-certified raw construction
    # singularities to an explicit weak-seed code.  A bare Error therefore
    # remains an error in every scope; the parser never trusts a role-specific
    # textual relabeling of code 8.
    return _generic_result_class(row["construct_result"])


def _failure_class(row):
    if row["construct_result"] != 0:
        return _expected_construct_class(row)
    if row["result_class"] != "success":
        return row["result_class"]
    if row["scope"] == "decoder" and row["recovery_class"] != "success":
        return row["recovery_class"]
    return "error"


def _outcome_class(row):
    return "success" if _row_succeeded(row) else _failure_class(row)


def _outcome_signature(row):
    return (
        row["construct_result"], row["result"], row["result_class"],
        row["recovery_result"], row["recovery_class"], row["recovery_ok"],
        row["encoded_symbols"], row["received_symbols"],
        row["arm_overhead"], row["fixed_prefix_symbols"],
        tuple(row[name] for name in _SOLVE_COUNT_FIELDS),
        tuple(row[name] for name in _BYTE_FIELDS),
        row["intermediate_sha256"],
    )


def _cross_panel_signature(row):
    """Bind scope receipts that native preflights once per arm/replicate."""
    signature = (
        row["construct_result"], row["result"], row["result_class"],
        row["recovery_result"], row["recovery_class"], row["recovery_ok"],
    )
    if row["scope"] in ("encoder", "decoder"):
        signature += (
            row["encoded_symbols"], row["received_symbols"],
            row["arm_overhead"], row["fixed_prefix_symbols"],
            tuple(row[name] for name in _BYTE_FIELDS),
            row["intermediate_sha256"],
        )
    elif row["scope"] == "direct":
        signature += (row["intermediate_sha256"],)
    return signature


def _validate_row_scope(
        row, *, block_count, block_bytes, candidate, dispatch,
        inner_reps, max_overhead, expected_cpu, schedule, loss,
        witnessed):
    label = (
        f"bandtiming replicate {row['replicate']} "
        f"{row['panel']} slot {row['slot']}"
    )
    base = _base_arm(row["role"])
    message_bytes = (block_count - 1) * block_bytes + 1
    if (base in ("candidate", "dispatch") and
            row["construct_result"] == 1):
        # Native normalizes a raw explicit construction NeedMore to the
        # explicit weak-seed code before any scope observation is emitted.
        raise MeasurementError(
            f"{label} bypassed explicit construction normalization")
    if (9 in (
                row["construct_result"], row["result"],
                row["recovery_result"]) or
            not 0 <= row["construct_result"] <= 10 or
            not -1 <= row["result"] <= 10 or
            not -1 <= row["recovery_result"] <= 10 or
            row["recovery_ok"] not in (0, 1) or
            row["recovery_class"] !=
                _generic_result_class(row["recovery_result"])):
        raise MeasurementError(
            f"{label} result class contradicts its code or reports OOM")
    if row["construct_result"] != 0:
        if (row["result"] != row["construct_result"] or
                row["result_class"] != _expected_construct_class(row)):
            raise MeasurementError(
                f"{label} construction failure class contradicts its code")
    elif (row["result"] == -1 or
            row["result_class"] != _generic_result_class(row["result"])):
        raise MeasurementError(f"{label} result class contradicts its code")
    if row["scope"] in ("encoder", "direct"):
        if (row["recovery_result"] != -1 or
                row["recovery_class"] != "not_applicable" or
                row["recovery_ok"] != 0):
            raise MeasurementError(
                f"{label} fabricated an inapplicable recovery result")
    elif row["result_class"] != "success":
        if (row["recovery_result"] != -1 or
                row["recovery_class"] != "not_applicable" or
                row["recovery_ok"] != 0):
            raise MeasurementError(
                f"{label} fabricated recovery after terminal failure")
    elif (row["recovery_result"] == -1 or
          row["recovery_class"] == "not_applicable" or
          ((row["recovery_class"] == "success") !=
           (row["recovery_ok"] == 1))):
        raise MeasurementError(
            f"{label} recovery class and output witness disagree")

    construction_failed = row["construct_result"] != 0
    zero_prefix = (
        row["encoded_symbols"] == 0 and
        row["received_symbols"] == 0 and
        row["arm_overhead"] == -1 and
        row["fixed_prefix_symbols"] == -1
    )
    if construction_failed:
        expected_symbols = zero_prefix
    elif row["scope"] == "encoder":
        expected_symbols = (
            row["encoded_symbols"] == block_count and
            row["received_symbols"] == 0 and
            row["arm_overhead"] == 0 and
            row["fixed_prefix_symbols"] == -1
        )
    elif row["scope"] == "decoder":
        received_in_range = (
            block_count <= row["received_symbols"] <=
                block_count + max_overhead
            if row["result_class"] == "success" else
            0 <= row["received_symbols"] <= block_count + max_overhead
        )
        expected_symbols = (
            row["encoded_symbols"] == 0 and
            received_in_range and
            row["arm_overhead"] ==
                (
                    row["received_symbols"] - block_count
                    if _row_succeeded(row) else -1
                ) and
            row["fixed_prefix_symbols"] == -1
        )
    else:
        expected_symbols = (
            row["encoded_symbols"] == 0 and
            (
                row["received_symbols"] == 0 or
                block_count <= row["received_symbols"] <=
                    block_count + max_overhead
            ) and
            row["arm_overhead"] ==
                row["received_symbols"] - block_count and
            row["fixed_prefix_symbols"] == row["received_symbols"]
        )
    if not expected_symbols:
        raise MeasurementError(
            f"{label} violates its scope-specific symbol counts")
    if (row["scope"] == "decoder" and
            row["construct_result"] == 0 and
            row["result"] == 1 and
            row["received_symbols"] not in (
                0, block_count + max_overhead)):
        # DecodeResult returns NeedMore after each nonterminal packet.  The
        # native loop can therefore receipt NeedMore only before it starts
        # (a shared-payload failure) or after it exhausts the complete
        # delivered schedule.
        raise MeasurementError(
            f"{label} stopped early while still reporting need-more")

    expected_stats = (
        row["timing_eligible"] == 1 and
        row["scope"] in ("decoder", "direct") and base != "wh1"
    )
    if row["stats_available"] != int(expected_stats):
        raise MeasurementError(f"{label} has dishonest stats availability")
    stats = tuple(row[name] for name in _SOLVE_STATS_FIELDS)
    if ((expected_stats and any(value < 0 for value in stats)) or
            (not expected_stats and any(value != -1 for value in stats))):
        raise MeasurementError(
            f"{label} fabricated or omitted solver statistics")

    descriptor = candidate if base == "candidate" else dispatch
    if expected_stats:
        inactivated = row["inactivated"]
        binary_def = row["binary_def"]
        heavy_gain = row["heavy_gain"]
        system_width = (
            block_count + descriptor.staircase + descriptor.dense_rows +
            descriptor.heavy_rows
        )
        phase_sum = sum(row[name] for name in _SOLVE_TIME_FIELDS)
        if (not 0 <= heavy_gain <= binary_def <= inactivated <= system_width or
                heavy_gain > descriptor.heavy_rows or
                (_row_succeeded(row) and heavy_gain != binary_def) or
                phase_sum > row["elapsed_ns"]):
            raise MeasurementError(
                f"{label} has impossible solver statistics")
    expected_intermediate = (
        0 if base == "wh1" else
        (
            block_count + descriptor.staircase + descriptor.dense_rows +
            descriptor.heavy_rows
        ) * block_bytes
    )
    if (row["source_bytes"] != message_bytes or
            row["intermediate_bytes"] != expected_intermediate):
        raise MeasurementError(
            f"{label} byte provenance contradicts its architecture")
    if witnessed:
        witness_expected = (
            base != "wh1" and _row_succeeded(row)
        )
        witness_valid = peel_codec._is_sha256(
            row["intermediate_sha256"])
        if witness_valid != witness_expected or (
                not witness_expected and
                row["intermediate_sha256"] != "not_applicable"):
            raise MeasurementError(
                f"{label} has an invalid intermediate witness")
    if zero_prefix:
        payload_ok = row["packet_payload_bytes"] == 0
    elif row["scope"] == "encoder":
        payload_ok = row["packet_payload_bytes"] == message_bytes
    elif row["scope"] == "direct" and not _row_succeeded(row):
        possible_remainders = (
            (0,) if schedule in ("repair-only", "adversarial") else (0, 1)
        )
        payload_ok = (
            0 <= row["packet_payload_bytes"] <=
                row["received_symbols"] * block_bytes and
            row["packet_payload_bytes"] % block_bytes in possible_remainders
        )
    else:
        # The only short packet is systematic id K-1, whose one-byte tail
        # reduces a complete prefix by exactly block_bytes - 1.  Native
        # payload construction cannot produce an arbitrary size in between.
        full_prefix_bytes = row["received_symbols"] * block_bytes
        possible_prefix_bytes = (
            (full_prefix_bytes,)
            if schedule in ("repair-only", "adversarial") else
            (
                full_prefix_bytes,
                full_prefix_bytes - (block_bytes - 1),
            )
        )
        payload_ok = row["packet_payload_bytes"] in possible_prefix_bytes
        expected_tail = _lossless_expected_tail(
            schedule, block_count, loss, row["received_symbols"])
        if payload_ok and expected_tail is not None:
            # With no loss these schedules necessarily deliver every
            # systematic id within the first K packets.  A complete prefix
            # of at least K packets must therefore include the one-byte
            # systematic tail at id K-1.
            payload_ok = row["packet_payload_bytes"] == (
                full_prefix_bytes -
                expected_tail * (block_bytes - 1))
    if not payload_ok:
        raise MeasurementError(
            f"{label} packet payload byte count is impossible")

    if row["timing_eligible"]:
        if (row["panel_censored"] != 0 or
                row["censor_reason"] != "none" or
                row["elapsed_ns"] <= 0 or
                row["inner_reps"] != inner_reps or
                row["saturated"] != 0 or
                row["cpu_before"] != expected_cpu or
                row["cpu_after"] != row["cpu_before"] or
                row["cpu_migrated"] != 0 or
                row["minflt_delta"] < 0 or row["majflt_delta"] < 0 or
                row["fault_contaminated"] != int(
                    row["minflt_delta"] > 0 or row["majflt_delta"] > 0) or
                row["preflight_result"] != row["result"] or
                row["timing_result"] != row["result"] or
                row["outcome_stable"] != 1):
            raise MeasurementError(
                f"{label} is contaminated or falsely timing-eligible")
    else:
        if (row["panel_censored"] != 1 or
                row["censor_reason"] == "none" or
                row["elapsed_ns"] != 0 or row["inner_reps"] != 0 or
                row["saturated"] != 0 or
                row["cpu_before"] != -1 or row["cpu_after"] != -1 or
                row["cpu_migrated"] != 0 or
                row["minflt_delta"] != 0 or
                row["majflt_delta"] != 0 or
                row["fault_contaminated"] != 0 or
                row["preflight_result"] != row["result"] or
                row["timing_result"] != -1 or
                row["outcome_stable"] != 0):
            raise MeasurementError(
                f"{label} did work or claimed stability while censored")


def _expected_censor_reason(panel_rows):
    panel = panel_rows[0]["panel"]
    unused_scope, left, right, unused_index = _PANEL_BY_NAME[panel]
    failures = []
    seen = set()
    for role in (left, right):
        base = _base_arm(role)
        if base in seen:
            continue
        seen.add(base)
        row = next(
            item for item in panel_rows
            if _base_arm(item["role"]) == base)
        if not _row_succeeded(row):
            failures.append(f"{base}_{_failure_class(row)}")
    return "none" if not failures else "+".join(failures)


def _aa_statistics(row_index, warmup_replicates, replicates, scope, arm):
    panel_name = f"{scope}_{arm}_aa"
    panel_index = _PANEL_BY_NAME[panel_name][3]
    replicate_logs = []
    for measured in range(replicates):
        absolute = warmup_replicates + measured
        selected = row_index[absolute].panel_rows[panel_index]
        if (not selected[0]["timing_eligible"] or
                any(row["fault_contaminated"] for row in selected)):
            continue
        pair_logs = []
        for pair in range(4):
            left = next(
                row for row in selected
                if row["pair"] == pair and row["role"] == f"{arm}_a")
            right = next(
                row for row in selected
                if row["pair"] == pair and row["role"] == f"{arm}_b")
            pair_logs.append(math.log(left["elapsed_ns"] / right["elapsed_ns"]))
        replicate_logs.append(sum(pair_logs) / 4.0)
    if len(replicate_logs) >= 4:
        mean, low, high = peel_codec._paired_t_interval(replicate_logs)
        floor = max(abs(low), abs(high))
    else:
        mean = (
            sum(replicate_logs) / len(replicate_logs)
            if replicate_logs else 0.0
        )
        low = high = None
        floor = 0.0
    return (
        len(replicate_logs), float(mean),
        None if low is None else float(low),
        None if high is None else float(high),
        float(floor),
    )


def _arm_summary(
        row_index, warmup_replicates, replicates, scope, arm,
        message_bytes):
    outcomes = []
    timed = []
    scope_arm_index = _ARM_SCOPE_INDEX[(scope, arm)]
    panel_indices = _ARM_PANEL_INDICES[(scope, arm)]
    primary_panel_index = _PANEL_BY_NAME[_PRIMARY_PANEL[(scope, arm)]][3]
    for measured in range(replicates):
        absolute = warmup_replicates + measured
        replicate_index = row_index[absolute]
        all_selected = replicate_index.scope_arm_rows[scope_arm_index]
        core_signatures = {
            (
                row["construct_result"], row["result"],
                row["result_class"], row["recovery_result"],
                row["recovery_class"], row["recovery_ok"],
            )
            for row in all_selected
        }
        if len(core_signatures) != 1:
            raise MeasurementError(
                f"bandtiming {scope}/{arm} outcome changed across panels "
                f"in replicate {absolute}")
        primary = [
            row
            for row in replicate_index.panel_rows[primary_panel_index]
            if _base_arm(row["role"]) == arm
        ]
        if not primary:
            raise MeasurementError(
                f"bandtiming {scope}/{arm} primary panel is missing")
        outcomes.append(primary[0])
        for panel_index in panel_indices:
            panel_rows = replicate_index.panel_rows[panel_index]
            if any(row["fault_contaminated"] for row in panel_rows):
                continue
            timed.extend(
                row for row in panel_rows
                if _base_arm(row["role"]) == arm and
                row["timing_eligible"])
    result_counts = {
        name: sum(_outcome_class(row) == name for row in outcomes)
        for name in _RESULT_CLASSES
    }
    overheads = [
        row["arm_overhead"] for row in outcomes if _row_succeeded(row)
    ]
    if overheads:
        mean = sum(overheads) / len(overheads)
        variance = (
            sum((value - mean) ** 2 for value in overheads) /
            len(overheads)
        )
    else:
        mean = variance = 0.0
    elapsed = sum(row["elapsed_ns"] for row in timed)
    work_bytes = sum(row["inner_reps"] * message_bytes for row in timed)
    throughput = work_bytes * 1000.0 / elapsed if elapsed else 0.0
    aa = _aa_statistics(
        row_index, warmup_replicates, replicates, scope, arm)
    return BandArmSummary(
        scope=scope,
        arm=arm,
        trials=replicates,
        successes=result_counts["success"],
        need_more=result_counts["need_more"],
        weak=result_counts["weak"],
        errors=result_counts["error"],
        overhead_mean=float(mean),
        overhead_sd=float(math.sqrt(variance)),
        overhead50=_percentile_integer(overheads, 0.50),
        overhead95=_percentile_integer(overheads, 0.95),
        overhead99=_percentile_integer(overheads, 0.99),
        overhead_max=float(max(overheads) if overheads else 0),
        timed_slots=len(timed),
        elapsed_ns=elapsed,
        throughput_mbps=float(throughput),
        aa_eligible_replicates=aa[0],
        aa_log_cost_mean=aa[1],
        aa_log_cost_ci_low=aa[2],
        aa_log_cost_ci_high=aa[3],
        aa_floor_log=aa[4],
    )


def _contrast(
        row_index, summaries, warmup_replicates, replicates, panel_name,
        required_margin):
    scope, left, right, panel_index = _PANEL_BY_NAME[panel_name]
    replicate_logs = []
    regressions = improvements = both_failures = 0
    for measured in range(replicates):
        absolute = warmup_replicates + measured
        selected = row_index[absolute].panel_rows[panel_index]
        left_row = next(row for row in selected if row["role"] == left)
        right_row = next(row for row in selected if row["role"] == right)
        left_ok = _row_succeeded(left_row)
        right_ok = _row_succeeded(right_row)
        if not left_ok and right_ok:
            regressions += 1
        elif left_ok and not right_ok:
            improvements += 1
        elif not left_ok and not right_ok:
            both_failures += 1
        if (not selected[0]["timing_eligible"] or
                any(row["fault_contaminated"] for row in selected)):
            continue
        pair_logs = []
        for pair in range(4):
            left_pair = next(
                row for row in selected
                if row["pair"] == pair and row["role"] == left)
            right_pair = next(
                row for row in selected
                if row["pair"] == pair and row["role"] == right)
            pair_logs.append(
                math.log(left_pair["elapsed_ns"] / right_pair["elapsed_ns"]))
        replicate_logs.append(sum(pair_logs) / 4.0)
    if len(replicate_logs) >= 4:
        mean, low, high = peel_codec._paired_t_interval(replicate_logs)
        available = True
    else:
        mean = (
            sum(replicate_logs) / len(replicate_logs)
            if replicate_logs else 0.0
        )
        low = high = None
        available = False
    left_summary = summaries[(scope, _base_arm(left))]
    right_summary = summaries[(scope, _base_arm(right))]
    lower, unused_upper = peel_codec._paired_practical_log_bounds(
        required_margin)
    required_log = -lower if lower < 0.0 else 0.0
    effective = max(
        required_log, left_summary.aa_floor_log,
        right_summary.aa_floor_log)
    return BandTimingContrast(
        name=panel_name,
        scope=scope,
        left_arm=_base_arm(left),
        right_arm=_base_arm(right),
        eligible_replicates=len(replicate_logs),
        log_cost_mean=float(mean),
        log_cost_ci_low=None if low is None else float(low),
        log_cost_ci_high=None if high is None else float(high),
        left_aa_floor_log=left_summary.aa_floor_log,
        right_aa_floor_log=right_summary.aa_floor_log,
        effective_floor_log=float(effective),
        required_margin=float(required_margin),
        recovery_regressions=regressions,
        recovery_improvements=improvements,
        both_failures=both_failures,
        timing_ci_available=available,
        left_faster=(
            available and
            left_summary.aa_eligible_replicates >= 4 and
            right_summary.aa_eligible_replicates >= 4 and
            regressions == 0 and high < -effective
        ),
    )


def _parse_stream_envelope(
        stdout, expected_rows, semantic_fields, columns):
    if not isinstance(stdout, str) or not stdout:
        raise MeasurementError("bandtiming returned no output")
    try:
        stdout.encode("ascii")
    except UnicodeEncodeError:
        raise MeasurementError("bandtiming output is not ASCII")
    lines = stdout.splitlines(keepends=True)
    if (not lines or
            any(not line.endswith("\n") or line.endswith("\r\n")
                for line in lines) or
            any(line == "\n" for line in lines)):
        raise MeasurementError(
            "bandtiming output is incomplete or noncanonical")
    if len(lines) != expected_rows + 4:
        raise MeasurementError(
            f"bandtiming emitted {len(lines) - 4} rows; "
            f"expected {expected_rows}")
    manifest_raw = peel_codec._parse_ordered_kv_line(
        lines[0][:-1], "# bandtiming,", BANDTIMING_MANIFEST_FIELDS,
        "bandtiming manifest")
    semantic_raw = peel_codec._parse_ordered_kv_line(
        lines[1][:-1], "# band_semantic,", semantic_fields,
        "bandtiming semantic")
    try:
        header = next(csv.reader([lines[2][:-1]], strict=True))
    except (csv.Error, StopIteration) as error:
        raise MeasurementError(f"malformed bandtiming header: {error}")
    if (lines[2][:-1] != ",".join(header) or
            tuple(header) != columns):
        raise MeasurementError("bandtiming header is missing or reordered")
    done = peel_codec._parse_ordered_kv_line(
        lines[-1][:-1], "# bandtiming_done,",
        ("complete", "rows", "finished_monotonic_ns", "stream_sha256"),
        "bandtiming completion")
    finished_ns = peel_codec._parse_decimal_integer(
        done["finished_monotonic_ns"],
        "bandtiming completion.finished_monotonic_ns",
        minimum=1, maximum=0xffffffffffffffff)
    if (done["complete"] != "1" or
            peel_codec._parse_decimal_integer(
                done["rows"], "bandtiming completion.rows",
                minimum=0, maximum=0xffffffff) != expected_rows or
            not peel_codec._is_sha256(done["stream_sha256"])):
        raise MeasurementError("bandtiming completion receipt is invalid")
    done_prefix = lines[-1][:-1].rsplit("stream_sha256=", 1)[0]
    done_prefix += "stream_sha256="
    stream_sha256 = hashlib.sha256(
        ("".join(lines[:-1]) + done_prefix).encode("ascii")).hexdigest()
    if done["stream_sha256"] != stream_sha256:
        raise MeasurementError("bandtiming stream SHA-256 does not match")
    return (
        lines, manifest_raw, semantic_raw, finished_ns, stream_sha256,
    )


def parse_bandtiming_output(
        stdout, *, block_count, block_bytes, candidate,
        dispatch_profile=peel_codec.TARGET_PROFILE,
        seed_policy=peel_codec.TARGET_SEED_POLICY,
        construction_seed, loss, loss_seed, schedule,
        warmup_replicates, replicates, inner_reps, max_overhead,
        cache_state, systematic_cache, evict_bytes, context,
        required_margin, _protocol=BANDTIMING_PROTOCOL_V2):
    """Validate and independently aggregate one native 15-panel stream."""
    expected_schema, semantic_fields, columns, witnessed = \
        _bandtiming_format(_protocol)
    validate_bandtiming_dimensions(
        block_count=block_count,
        block_bytes=block_bytes,
        candidate=candidate,
        dispatch_profile=dispatch_profile,
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
        systematic_cache=systematic_cache,
        evict_bytes=evict_bytes,
        required_margin=required_margin,
    )
    total_replicates = warmup_replicates + replicates
    expected_rows = total_replicates * 120
    (
        lines, manifest_raw, semantic_raw, finished_ns, stream_sha256,
    ) = _parse_stream_envelope(
        stdout, expected_rows, semantic_fields, columns)
    started_ns = peel_codec._parse_decimal_integer(
        manifest_raw["started_monotonic_ns"],
        "bandtiming manifest.started_monotonic_ns",
        minimum=1, maximum=0xffffffffffffffff)
    if finished_ns <= started_ns:
        raise MeasurementError("bandtiming completion clock did not advance")
    context_sha256 = peel_codec._validate_paired_context(
        context, cache_state=cache_state, evict_bytes=evict_bytes,
        started_ns=started_ns, finished_ns=finished_ns)
    manifest = _parse_manifest(
        manifest_raw,
        block_count=block_count,
        block_bytes=block_bytes,
        candidate=candidate,
        dispatch_profile=dispatch_profile,
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
        systematic_cache=systematic_cache,
        evict_bytes=evict_bytes,
        context_sha256=context_sha256,
        required_margin=float(required_margin),
        expected_rows=expected_rows,
        expected_schema=expected_schema,
    )
    if manifest["started_monotonic_ns"] != started_ns:
        raise MeasurementError("bandtiming start clock changed while parsing")
    dispatch = dispatch_band_descriptor(block_count)

    indexed_rows = tuple(
        _ReplicateRowIndex(
            panel_rows=tuple([] for unused in BANDTIMING_PANELS),
            scope_arm_rows=tuple([] for unused in _ARM_SCOPE_ORDER),
        )
        for unused_replicate in range(total_replicates)
    )
    seen_trace = {}
    total_timed_ns = 0
    for row_index, line in enumerate(lines[3:-1]):
        try:
            fields = next(csv.reader([line[:-1]], strict=True))
        except (csv.Error, StopIteration) as error:
            raise MeasurementError(
                f"malformed bandtiming row {row_index}: {error}")
        if line[:-1] != ",".join(fields):
            raise MeasurementError(
                f"bandtiming row {row_index} uses noncanonical CSV")
        row = _parse_row(
            fields, f"bandtiming row {row_index}", columns)
        expected_replicate = row_index // 120
        within_replicate = row_index % 120
        expected_panel_index = within_replicate // 8
        expected_slot = within_replicate % 8
        scope, panel, unused_left, unused_right = \
            BANDTIMING_PANEL_SPECS[expected_panel_index]
        expected_label = BANDTIMING_ORDER[expected_slot]
        expected_role = _expected_role(
            panel, expected_replicate, expected_label)
        if (row["replicate"] != expected_replicate or
                row["measured"] !=
                    int(expected_replicate >= warmup_replicates) or
                row["scope"] != scope or row["panel"] != panel or
                row["panel_index"] != expected_panel_index or
                row["slot"] != expected_slot or
                row["pair"] != expected_slot // 2 or
                row["label"] != expected_label or
                row["role"] != expected_role or
                row["label_swap"] != (expected_replicate & 1)):
            raise MeasurementError(
                f"bandtiming row {row_index} is reordered or mislabeled")
        expected_construction = _expected_seed(
            construction_seed, expected_replicate)
        expected_loss = _expected_loss_seed(loss_seed, expected_replicate)
        if row["construction_seed"] != expected_construction:
            raise MeasurementError(
                f"bandtiming row {row_index} has the wrong construction seed")
        if row["loss_seed"] != expected_loss:
            raise MeasurementError(
                f"bandtiming row {row_index} has the wrong loss seed")
        trace = seen_trace.setdefault(
            expected_replicate,
            (
                row["construction_seed"], row["loss_seed"],
                row["trace_sha256"],
            ))
        if trace != (
                row["construction_seed"], row["loss_seed"],
                row["trace_sha256"]):
            raise MeasurementError(
                f"bandtiming replicate {expected_replicate} changed trace")
        _validate_row_scope(
            row,
            block_count=block_count,
            block_bytes=block_bytes,
            candidate=candidate,
            dispatch=dispatch,
            inner_reps=inner_reps,
            max_overhead=max_overhead,
            expected_cpu=context["bound"]["cpu_affinity"][0],
            schedule=schedule,
            loss=float(loss),
            witnessed=witnessed,
        )
        replicate_index = indexed_rows[expected_replicate]
        replicate_index.panel_rows[expected_panel_index].append(row)
        replicate_index.scope_arm_rows[
            _ARM_SCOPE_INDEX[(scope, _base_arm(expected_role))]
        ].append(row)
        total_timed_ns += row["elapsed_ns"]
    if total_timed_ns > finished_ns - started_ns:
        raise MeasurementError(
            "bandtiming aggregate timed interval exceeds the authenticated "
            "sequential run interval")
    if (len(seen_trace) != total_replicates or
            len({item[0] for item in seen_trace.values()}) !=
                total_replicates or
            len({item[1] for item in seen_trace.values()}) !=
                total_replicates or
            len({item[2] for item in seen_trace.values()}) !=
                total_replicates):
        raise MeasurementError(
            "bandtiming seeds or traces are missing, reused, or aliased")

    for replicate in range(total_replicates):
        replicate_index = indexed_rows[replicate]
        base_counts = {
            arm: sum(
                len(replicate_index.scope_arm_rows[
                    _ARM_SCOPE_INDEX[(scope, arm)]
                ])
                for scope in (
                    ("encoder", "decoder", "direct")
                    if arm != "wh1" else ("encoder", "decoder")
                )
            )
            for arm in ("candidate", "dispatch", "wh1")
        }
        if base_counts != {"candidate": 44, "dispatch": 44, "wh1": 32}:
            raise MeasurementError(
                f"bandtiming replicate {replicate} has wrong arm cardinality")
        for panel_index, panel in enumerate(BANDTIMING_PANELS):
            panel_scope = BANDTIMING_PANEL_SPECS[panel_index][0]
            selected = replicate_index.panel_rows[panel_index]
            if len(selected) != 8:
                raise MeasurementError(
                    f"bandtiming replicate {replicate} panel {panel} "
                    "is incomplete")
            per_role = {}
            for row in selected:
                per_role.setdefault(row["role"], []).append(row)
            if set(per_role) != {
                    BANDTIMING_PANEL_SPECS[panel_index][2],
                    BANDTIMING_PANEL_SPECS[panel_index][3]} or \
                    any(len(items) != 4 for items in per_role.values()) or \
                    any(len({_outcome_signature(row) for row in items}) != 1
                        for items in per_role.values()):
                raise MeasurementError(
                    f"bandtiming replicate {replicate} panel {panel} "
                    "changed arm semantics")
            if panel.endswith("_aa"):
                signatures = {
                    _outcome_signature(items[0])
                    for items in per_role.values()
                }
                if len(signatures) != 1:
                    raise MeasurementError(
                        f"bandtiming {panel} A/A arms differ semantically")
            expected_reason = _expected_censor_reason(selected)
            expected_eligible = int(expected_reason == "none")
            if ({row["timing_eligible"] for row in selected} !=
                    {expected_eligible} or
                    {row["panel_censored"] for row in selected} !=
                    {1 - expected_eligible} or
                    {row["censor_reason"] for row in selected} !=
                    {expected_reason}):
                raise MeasurementError(
                    f"bandtiming replicate {replicate} panel {panel} "
                    "has a forged censor receipt")
            if panel_scope == "direct" and expected_eligible:
                for pair in range(4):
                    prefixes = {
                        row["fixed_prefix_symbols"] for row in selected
                        if row["pair"] == pair
                    }
                    if len(prefixes) != 1:
                        raise MeasurementError(
                            f"bandtiming direct panel {panel} pair {pair} "
                            "did not share one fixed prefix")
            if panel_scope == "direct":
                panel_bases = {
                    _base_arm(row["role"]) for row in selected
                }
                decoder_rows = {
                    arm: replicate_index.scope_arm_rows[
                        _ARM_SCOPE_INDEX[("decoder", arm)]
                    ][0]
                    for arm in panel_bases
                }
                expected_prefix = (
                    max(row["received_symbols"]
                        for row in decoder_rows.values())
                    if all(_row_succeeded(row)
                           for row in decoder_rows.values()) else
                    block_count + max_overhead
                )
                constructed = [
                    row for row in selected
                    if row["construct_result"] == 0
                ]
                if (any(
                        row["received_symbols"] != expected_prefix or
                        row["fixed_prefix_symbols"] != expected_prefix or
                        row["arm_overhead"] !=
                            expected_prefix - block_count
                        for row in constructed) or
                        (expected_eligible and len({
                            row["packet_payload_bytes"]
                            for row in constructed
                        }) > 1)):
                    raise MeasurementError(
                        f"bandtiming direct panel {panel} did not use its "
                        "decoder-bound fixed prefix and payload")
        for scope, arm in _ARM_SCOPE_ORDER:
            core_signatures = {
                _cross_panel_signature(row)
                for row in replicate_index.scope_arm_rows[
                    _ARM_SCOPE_INDEX[(scope, arm)]
                ]
            }
            if len(core_signatures) != 1:
                raise MeasurementError(
                    f"bandtiming {scope}/{arm} outcome changed across panels "
                    f"in replicate {replicate}")

        decoder_by_arm = {
            arm: replicate_index.scope_arm_rows[
                _ARM_SCOPE_INDEX[("decoder", arm)]
            ][0]
            for arm in ("candidate", "dispatch", "wh1")
        }
        for arm in ("candidate", "dispatch", "wh1"):
            shared_constructor_rows = [
                row
                for scope in (
                    ("encoder", "decoder", "direct")
                    if arm != "wh1" else ("encoder", "decoder")
                )
                for row in replicate_index.scope_arm_rows[
                    _ARM_SCOPE_INDEX[(scope, arm)]
                ]
            ]
            constructor_codes = {
                row["construct_result"] for row in shared_constructor_rows
            }
            if (len(constructor_codes) != 1 or
                    (
                        next(iter(constructor_codes)) != 0 and
                        {
                            row["result"]
                            for row in shared_constructor_rows
                        } != constructor_codes
                    )):
                raise MeasurementError(
                    f"bandtiming {arm} constructor code changed across "
                    f"shared payload scopes in replicate {replicate}")
        if witnessed:
            for arm in ("candidate", "dispatch"):
                witnesses_by_scope = {
                    scope: {
                        row["intermediate_sha256"]
                        for row in replicate_index.scope_arm_rows[
                            _ARM_SCOPE_INDEX[(scope, arm)]
                        ]
                        if row["intermediate_sha256"] != "not_applicable"
                    }
                    for scope in BANDTIMING_SCOPE_ORDER
                }
                witnesses = set().union(*witnesses_by_scope.values())
                if len(witnesses) > 1 or (
                        any(witnesses_by_scope[scope]
                            for scope in ("decoder", "direct")) and
                        not witnesses_by_scope["encoder"]):
                    raise MeasurementError(
                        f"bandtiming {arm} intermediate witness changed "
                        f"across scopes in replicate {replicate}")
        complete_prefix_rows = [
            row for row in decoder_by_arm.values()
            if row["construct_result"] == 0 and
            row["received_symbols"] > 0
        ]
        for arm in ("candidate", "dispatch"):
            decoder_row = decoder_by_arm[arm]
            direct_rows = replicate_index.scope_arm_rows[
                _ARM_SCOPE_INDEX[("direct", arm)]
            ]
            # A constructed decoder that consumed no packet can only be the
            # early shared-payload failure path.  Direct preflight sees the
            # same payload object and must copy that exact failure and byte
            # prefix into every panel.
            if (decoder_row["construct_result"] == 0 and
                    decoder_row["received_symbols"] == 0):
                failure_signatures = {
                    (
                        row["construct_result"], row["result"],
                        row["result_class"], row["received_symbols"],
                        row["arm_overhead"], row["fixed_prefix_symbols"],
                        row["packet_payload_bytes"],
                    )
                    for row in direct_rows
                }
                if (len(failure_signatures) != 1 or
                        next(iter(failure_signatures))[:3] != (
                            decoder_row["construct_result"],
                            decoder_row["result"],
                            decoder_row["result_class"],
                        )):
                    raise MeasurementError(
                        f"bandtiming {arm} shared payload failure changed "
                        f"between decoder and direct scopes in replicate "
                        f"{replicate}")
                if any(
                        row["received_symbols"] <= 0 or
                        row["packet_payload_bytes"] >
                            (row["received_symbols"] - 1) * block_bytes
                        for row in direct_rows):
                    raise MeasurementError(
                        f"bandtiming {arm} shared payload failure claimed "
                        f"a complete direct prefix in replicate {replicate}")
            # Once the decoder consumed a packet, its payload was completely
            # constructed.  A successful direct solve proves the same fact
            # independently.  Such rows must use a complete schedule prefix.
            complete_prefix_rows.extend(
                row for row in direct_rows
                if row["construct_result"] == 0 and
                (
                    decoder_row["received_symbols"] > 0 or
                    _row_succeeded(row)
                )
            )

            direct_counts_by_prefix = {}
            for row in direct_rows:
                if row["stats_available"] != 1:
                    continue
                direct_counts_by_prefix.setdefault(
                    row["fixed_prefix_symbols"], set()).add(
                        tuple(row[name] for name in _SOLVE_COUNT_FIELDS))
            if any(
                    len(counts) > 1
                    for counts in direct_counts_by_prefix.values()):
                raise MeasurementError(
                    f"bandtiming direct/{arm} solver counts changed for one "
                    f"fixed prefix in replicate {replicate}")

        # Every codec arm reuses the same id schedule.  For a complete prefix,
        # the payload is either n full blocks or exactly one short systematic
        # tail.  Equal prefixes must agree, and after that tail first appears
        # it remains included in every longer prefix.
        short_by_prefix = {}
        for row in complete_prefix_rows:
            full_bytes = row["received_symbols"] * block_bytes
            if row["packet_payload_bytes"] == full_bytes:
                short = 0
            elif row["packet_payload_bytes"] == (
                    full_bytes - (block_bytes - 1)):
                short = 1
            else:
                raise MeasurementError(
                    f"bandtiming shared payload prefix is incomplete in "
                    f"replicate {replicate}")
            short_by_prefix.setdefault(
                row["received_symbols"], set()).add(short)

        # A shared-payload failure leaves the failed and subsequent
        # DataBytes entries at zero.  Its direct receipt still authenticates
        # the exact number of packets completed before the failure and
        # whether that prefix already contained the one-byte systematic
        # tail.  Fold those partial prefixes into the same common-schedule
        # history instead of discarding their information.
        for arm in ("candidate", "dispatch"):
            decoder_row = decoder_by_arm[arm]
            if (decoder_row["construct_result"] != 0 or
                    decoder_row["received_symbols"] != 0):
                continue
            direct_rows = replicate_index.scope_arm_rows[
                _ARM_SCOPE_INDEX[("direct", arm)]
            ]
            for row in direct_rows:
                payload_bytes = row["packet_payload_bytes"]
                remainder = payload_bytes % block_bytes
                if remainder == 0:
                    completed = payload_bytes // block_bytes
                    short = 0
                elif remainder == 1:
                    completed = (
                        payload_bytes + block_bytes - 1) // block_bytes
                    short = 1
                else:
                    # Row validation rejects this first; keep this branch
                    # fail-closed if the validation order ever changes.
                    raise MeasurementError(
                        f"bandtiming {arm} shared payload failure has an "
                        f"impossible partial prefix in replicate {replicate}")
                if completed >= row["received_symbols"]:
                    raise MeasurementError(
                        f"bandtiming {arm} shared payload failure claimed "
                        f"a complete direct prefix in replicate {replicate}")
                expected_tail = _lossless_expected_tail(
                    schedule, block_count, loss, completed)
                if expected_tail is not None and short != expected_tail:
                    raise MeasurementError(
                        f"bandtiming {arm} shared payload failure omitted "
                        f"or fabricated the deterministic systematic tail "
                        f"in replicate "
                        f"{replicate}")
                short_by_prefix.setdefault(completed, set()).add(short)
        ordered_prefixes = sorted(short_by_prefix)
        deadline = _tail_candidate_deadline(schedule, block_count)
        if (any(len(values) != 1 for values in short_by_prefix.values()) or
                any(
                    (
                        next(iter(short_by_prefix[left])) >
                            next(iter(short_by_prefix[right]))
                    ) or (
                        deadline is not None and left >= deadline and
                        next(iter(short_by_prefix[left])) == 0 and
                        next(iter(short_by_prefix[right])) == 1
                    )
                    for left, right in zip(
                        ordered_prefixes, ordered_prefixes[1:]))):
            raise MeasurementError(
                f"bandtiming shared payload prefix changed packet-byte "
                f"history in replicate {replicate}")

        for arm in ("candidate", "dispatch"):
            decoder_counts = {
                tuple(row[name] for name in _SOLVE_COUNT_FIELDS)
                for row in replicate_index.scope_arm_rows[
                    _ARM_SCOPE_INDEX[("decoder", arm)]
                ]
                if row["stats_available"] == 1
            }
            if len(decoder_counts) > 1:
                raise MeasurementError(
                    f"bandtiming decoder/{arm} solver counts changed across "
                    f"panels in replicate {replicate}")

    semantic = _parse_semantic(
        semantic_raw,
        dispatch,
        block_count,
        block_bytes,
        construction_seed,
        {item[0] for item in seen_trace.values()},
        {item[2] for item in seen_trace.values()},
        max_overhead,
        witnessed,
    )
    message_bytes = (block_count - 1) * block_bytes + 1
    summaries = {
        (scope, arm): _arm_summary(
            indexed_rows, warmup_replicates, replicates, scope, arm,
            message_bytes)
        for scope, arm in _ARM_SCOPE_ORDER
    }
    contrasts = tuple(
        _contrast(
            indexed_rows, summaries, warmup_replicates, replicates, panel,
            float(required_margin))
        for panel in BANDTIMING_CROSS_PANELS
    )
    weak_cells = []
    replicate_receipts = []
    for measured in range(replicates):
        absolute = warmup_replicates + measured
        replicate_index = indexed_rows[absolute]
        trace = seen_trace[absolute]
        outcomes = {}
        construction = {}
        censored_panels = []
        contaminated_panels = []
        for arm in ("candidate", "dispatch", "wh1"):
            selected = replicate_index.scope_arm_rows[
                _ARM_SCOPE_INDEX[("encoder", arm)]
            ][0]
            construct_result = selected["construct_result"]
            construction[f"{arm}_construct_result"] = construct_result
            construction[f"{arm}_construct_class"] = \
                _generic_result_class(construct_result)
        for scope, arm in _ARM_SCOPE_ORDER:
            selected = replicate_index.scope_arm_rows[
                _ARM_SCOPE_INDEX[(scope, arm)]
            ][0]
            outcomes[f"{scope}_{arm}_result_class"] = \
                _outcome_class(selected)
            outcomes[f"{scope}_{arm}_overhead"] = selected["arm_overhead"]
        for panel_index, panel in enumerate(BANDTIMING_PANELS):
            selected = replicate_index.panel_rows[panel_index]
            if selected[0]["panel_censored"]:
                censored_panels.append(panel)
            if any(row["fault_contaminated"] for row in selected):
                contaminated_panels.append(panel)
        for scope in BANDTIMING_SCOPE_ORDER:
            arms = (
                ("candidate", "dispatch", "wh1")
                if scope != "direct" else ("candidate", "dispatch")
            )
            weak = tuple(
                arm for arm in arms
                if outcomes[f"{scope}_{arm}_result_class"] == "weak"
            )
            if weak:
                if "candidate" in weak and len(weak) > 1:
                    category = "shared"
                elif weak == ("candidate",):
                    category = "candidate-only"
                else:
                    category = "control-only"
                weak_cells.append({
                    "replicate": measured,
                    "scope": scope,
                    "construction_seed": trace[0],
                    "loss_seed": trace[1],
                    "trace_sha256": trace[2],
                    "category": category,
                    "weak_arms": list(weak),
                })
        replicate_receipts.append({
            "replicate": measured,
            "construction_seed": trace[0],
            "loss_seed": trace[1],
            "trace_sha256": trace[2],
            **construction,
            **outcomes,
            "censored_panels": list(censored_panels),
            "fault_contaminated_panels": list(contaminated_panels),
        })
    return BandTimingMeasurement(
        protocol=_protocol,
        manifest=manifest,
        semantic=semantic,
        context=copy.deepcopy(context),
        candidate_descriptor=candidate,
        dispatch_descriptor=dispatch,
        arms=tuple(summaries[key] for key in _ARM_SCOPE_ORDER),
        contrasts=contrasts,
        replicate_receipts=tuple(replicate_receipts),
        weak_cells=tuple(weak_cells),
        stream_sha256=stream_sha256,
        finished_monotonic_ns=finished_ns,
        final_context_sha256=peel_codec._canonical_json_sha256(context),
        native_stdout=stdout,
    )


def _bandtiming_run_bounds(stdout):
    """Extract authenticated run bounds for the shared thermal finalizer."""
    if not isinstance(stdout, str):
        raise MeasurementError("bandtiming output is not text")
    try:
        stdout.encode("ascii")
    except UnicodeEncodeError:
        raise MeasurementError("bandtiming output is not ASCII")
    lines = stdout.splitlines(keepends=True)
    if (len(lines) < 4 or
            any(not line.endswith("\n") or line.endswith("\r\n")
                for line in lines)):
        raise MeasurementError("bandtiming output is incomplete")
    manifest = peel_codec._parse_ordered_kv_line(
        lines[0][:-1], "# bandtiming,", BANDTIMING_MANIFEST_FIELDS,
        "bandtiming manifest")
    done = peel_codec._parse_ordered_kv_line(
        lines[-1][:-1], "# bandtiming_done,",
        ("complete", "rows", "finished_monotonic_ns", "stream_sha256"),
        "bandtiming completion")
    expected_rows = peel_codec._parse_decimal_integer(
        manifest["expected_rows"], "bandtiming manifest.expected_rows",
        minimum=1, maximum=1_200_000)
    rows = peel_codec._parse_decimal_integer(
        done["rows"], "bandtiming completion.rows",
        minimum=1, maximum=1_200_000)
    done_prefix = lines[-1][:-1].rsplit("stream_sha256=", 1)[0]
    done_prefix += "stream_sha256="
    stream_sha256 = hashlib.sha256(
        ("".join(lines[:-1]) + done_prefix).encode("ascii")).hexdigest()
    if (manifest["schema"] != BANDTIMING_SCHEMA_V2 or
            done["complete"] != "1" or rows != expected_rows or
            len(lines) != rows + 4 or
            done["stream_sha256"] != stream_sha256):
        raise MeasurementError(
            "cannot finalize thermal evidence for incomplete bandtiming")
    started_ns = peel_codec._parse_decimal_integer(
        manifest["started_monotonic_ns"],
        "bandtiming manifest.started_monotonic_ns",
        minimum=1, maximum=0xffffffffffffffff)
    finished_ns = peel_codec._parse_decimal_integer(
        done["finished_monotonic_ns"],
        "bandtiming completion.finished_monotonic_ns",
        minimum=1, maximum=0xffffffffffffffff)
    if finished_ns <= started_ns:
        raise MeasurementError("bandtiming completion clock did not advance")
    return started_ns, finished_ns


def bandtiming_probe(
        bench, block_count, block_bytes, candidate, *,
        dispatch_profile=peel_codec.TARGET_PROFILE,
        seed_policy=peel_codec.TARGET_SEED_POLICY,
        construction_seed, loss, loss_seed, schedule,
        warmup_replicates, replicates, inner_reps, max_overhead,
        cache_state, systematic_cache, evict_bytes, context,
        required_margin):
    """Run one native bandtiming process with retained thermal evidence."""
    validate_bandtiming_dimensions(
        block_count=block_count,
        block_bytes=block_bytes,
        candidate=candidate,
        dispatch_profile=dispatch_profile,
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
        systematic_cache=systematic_cache,
        evict_bytes=evict_bytes,
        required_margin=required_margin,
    )
    initial_context = (
        peel_codec.read_paired_context(context, require_thermal=False)
        if isinstance(context, str) else context
    )
    frozen_context, capture = peel_codec._prepare_paired_context(
        initial_context, cache_state=cache_state, evict_bytes=evict_bytes)
    try:
        context_sha256 = peel_codec._validate_paired_bound_context(
            frozen_context, cache_state=cache_state, evict_bytes=evict_bytes,
            require_capture=True)
        command = [
            bench, "bandtiming",
            "--N", str(block_count),
            "--bb", str(block_bytes),
            "--dispatch-profile", dispatch_profile,
            "--seed-policy", seed_policy,
            "--construction-seed", str(construction_seed),
            "--loss", f"{float(loss):.17g}",
            "--loss-seed", str(loss_seed),
            "--schedule", schedule,
            "--candidate-staircase", str(candidate.staircase),
            "--candidate-dense-rows", str(candidate.dense_rows),
            "--candidate-gf256-rows", str(candidate.gf256_rows),
            "--candidate-gf16-rows", str(candidate.gf16_rows),
            "--candidate-period", str(candidate.period),
            "--candidate-x-geometry", candidate.x_mode,
            "--warmup-replicates", str(warmup_replicates),
            "--replicates", str(replicates),
            "--inner-reps", str(inner_reps),
            "--max-overhead", str(max_overhead),
            "--cache-state", cache_state,
            "--systematic-cache", systematic_cache,
            "--evict-bytes", str(evict_bytes),
            "--context-sha256", context_sha256,
            "--required-margin", f"{float(required_margin):.17g}",
        ]
        stdout = peel_codec._run_checked(
            command, peel_codec.isolated_codec_env())
        final_context = peel_codec._finalize_paired_context(
            frozen_context, stdout, capture,
            run_bounds=_bandtiming_run_bounds)
    finally:
        os.close(capture.fd)
    return parse_bandtiming_output(
        stdout,
        block_count=block_count,
        block_bytes=block_bytes,
        candidate=candidate,
        dispatch_profile=dispatch_profile,
        seed_policy=seed_policy,
        construction_seed=construction_seed,
        loss=float(loss),
        loss_seed=loss_seed,
        schedule=schedule,
        warmup_replicates=warmup_replicates,
        replicates=replicates,
        inner_reps=inner_reps,
        max_overhead=max_overhead,
        cache_state=cache_state,
        systematic_cache=systematic_cache,
        evict_bytes=evict_bytes,
        context=final_context,
        required_margin=float(required_margin),
    )


def _descriptor_from_receipt(value, label):
    if (not isinstance(value, dict) or
            set(value) != {"S", "D2", "gf256", "gf16", "P", "x"} or
            any(type(value[name]) is not int
                for name in ("S", "D2", "gf256", "gf16", "P")) or
            not isinstance(value["x"], str)):
        raise MeasurementError(f"{label} is malformed")
    return BandDescriptor(
        staircase=value["S"],
        dense_rows=value["D2"],
        gf256_rows=value["gf256"],
        gf16_rows=value["gf16"],
        period=value["P"],
        x_mode=value["x"],
    )


def replay_bandtiming_receipt(receipt, *, expected_request=None):
    """Reparse native evidence and reject stale or forged derived summaries."""
    fields = {
        "protocol", "manifest", "semantic", "context",
        "candidate_descriptor", "dispatch_descriptor", "arms", "contrasts",
        "replicates", "weak_cells", "stream_sha256", "evidence",
        "native_stdout",
    }
    if not isinstance(receipt, dict) or set(receipt) != fields:
        raise MeasurementError("bandtiming receipt schema is incomplete")
    protocol = receipt.get("protocol")
    if (not isinstance(protocol, str) or
            protocol not in _BANDTIMING_SCHEMA_BY_PROTOCOL):
        raise MeasurementError("bandtiming receipt protocol is invalid")
    manifest = receipt.get("manifest")
    if not isinstance(manifest, dict):
        raise MeasurementError("bandtiming receipt manifest is invalid")
    if manifest.get("schema") != _BANDTIMING_SCHEMA_BY_PROTOCOL[protocol]:
        raise MeasurementError(
            "bandtiming receipt protocol and native schema disagree")
    candidate = _descriptor_from_receipt(
        receipt.get("candidate_descriptor"),
        "bandtiming candidate descriptor")
    try:
        expected_dispatch = dispatch_band_descriptor(manifest.get("K"))
    except ValueError as error:
        raise MeasurementError(
            f"bandtiming receipt manifest is invalid: {error}")
    if _descriptor_from_receipt(
            receipt.get("dispatch_descriptor"),
            "bandtiming dispatch descriptor") != expected_dispatch:
        raise MeasurementError("bandtiming dispatch descriptor is forged")
    request = {
        "block_count": manifest.get("K"),
        "block_bytes": manifest.get("bb"),
        "candidate": candidate,
        "dispatch_profile": manifest.get("dispatch_profile"),
        "seed_policy": manifest.get("seed_policy"),
        "construction_seed": manifest.get("construction_seed_base"),
        "loss": manifest.get("loss"),
        "loss_seed": manifest.get("loss_seed_base"),
        "schedule": manifest.get("schedule"),
        "warmup_replicates": manifest.get("warmup_replicates"),
        "replicates": manifest.get("replicates"),
        "inner_reps": manifest.get("inner_reps"),
        "max_overhead": manifest.get("max_overhead"),
        "cache_state": manifest.get("cache_state"),
        "systematic_cache": manifest.get("systematic_cache"),
        "evict_bytes": manifest.get("evict_bytes"),
        "context": receipt.get("context"),
        "required_margin": manifest.get("required_margin"),
    }
    if expected_request is not None:
        if not isinstance(expected_request, dict):
            raise MeasurementError("expected bandtiming request is not a dict")
        for name, expected in expected_request.items():
            if name == "candidate":
                same = (
                    type(expected) is BandDescriptor
                    and peel_codec._same_typed_json(
                        request[name].as_dict(), expected.as_dict())
                )
            else:
                same = (
                    name in request
                    and peel_codec._same_typed_json(
                        request[name], expected)
                )
            if name not in request or not same:
                raise MeasurementError(
                    f"bandtiming replay request changed {name}")
    try:
        parsed = parse_bandtiming_output(
            receipt.get("native_stdout"), _protocol=protocol, **request)
    except ValueError as error:
        raise MeasurementError(
            f"bandtiming receipt request is invalid: {error}")
    replayed = parsed.as_dict()
    if peel_codec._same_typed_json(replayed, receipt):
        return parsed

    # Pre-enrichment v1 receipts are durable evidence already archived by the
    # benchmark campaign.  Accept only their exact replicate schema, and only
    # when the entire receipt matches a fresh replay after removing precisely
    # the six constructor fields that did not exist at the time.  Partially
    # enriched or otherwise incomplete receipts remain invalid.
    legacy_replicates = receipt.get("replicates")
    if (
        protocol == BANDTIMING_PROTOCOL_V1
        and isinstance(legacy_replicates, list)
        and legacy_replicates
        and all(
            isinstance(replicate, dict)
            and set(replicate) == _LEGACY_REPLICATE_FIELDS
            for replicate in legacy_replicates
        )
    ):
        legacy_replayed = copy.deepcopy(replayed)
        for replicate in legacy_replayed["replicates"]:
            for field in _REPLICATE_CONSTRUCTION_FIELDS:
                del replicate[field]
        if peel_codec._same_typed_json(legacy_replayed, receipt):
            return parsed

    raise MeasurementError(
        "bandtiming receipt has stale or forged derived fields")
