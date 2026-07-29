#!/usr/bin/env python3
"""Strict independent consumer for WH2 repair-v1 v3 timing evidence.

The native benchmark owns measurement.  This module owns the wire-format
contract, rejects malformed or internally impossible evidence, and exposes a
lossless normalized receipt.  It deliberately performs no campaign outcome
aggregation and does not alter the durable bandtiming v1/v2 consumers.
"""

import copy
import csv
import ctypes
import functools
import hashlib
import io
import json
import math
import os
import selectors
import signal
import stat
import struct
import subprocess
import sys
import time
import zlib
from dataclasses import dataclass
from typing import Optional

import peel_codec


MeasurementError = peel_codec.MeasurementError

REPAIRTIMING_PROTOCOL_V3 = (
    "wirehair-v2-bench:repairtiming:repair-v1:v3"
)
REPAIRTIMING_SCHEMA_V3 = "wirehair.wh2.repairtiming.v3"
REPAIRTIMING_PROTOCOL = REPAIRTIMING_PROTOCOL_V3
REPAIRTIMING_SCHEMA = REPAIRTIMING_SCHEMA_V3

REPAIR_V1_POLICY_NAME = "repair-v1"
REPAIR_V1_POLICY_SHA256 = (
    "5e67a150d1f909d6ed80468185fa2dd0"
    "e82eb2fc3486c0fa662e213cf3100b42"
)
REPAIR_V1_ATTEMPT_CAP = 8
REPAIR_V1_DESCRIPTOR_MAGIC = "W2R3"
REPAIR_V1_DESCRIPTOR_ENCODING = (
    "ascii-W2R3-u64le-contract-id-u32le-root-u8-attempt-v1"
)
REPAIR_V1_DESCRIPTOR_BYTES = 17
REPAIR_V1_DESCRIPTOR_POLICY = (
    "exact-length-two-known-ids-attempt-lt8-transactional-v1"
)

REPAIRTIMING_CONSTRUCTION_SEED_DERIVED = (
    "base_xor_d192ed03_times_rep_plus_1_mod2^32_v1"
)
REPAIRTIMING_CONSTRUCTION_SEED_FIXED = "fixed_base_u32_v1"
REPAIRTIMING_CONSTRUCTION_SEED_DERIVATIONS = frozenset({
    REPAIRTIMING_CONSTRUCTION_SEED_DERIVED,
    REPAIRTIMING_CONSTRUCTION_SEED_FIXED,
})
REPAIRTIMING_LOSS_SEED_DERIVATION = (
    "base_xor_9e3779b97f4a7c15_times_rep_plus_1_v1"
)
REPAIRTIMING_SCHEDULES = frozenset({
    "iid", "burst", "permutation", "systematic-first", "repair-only",
    "adversarial",
})
REPAIRTIMING_ORDER = "ABBABAAB"
REPAIRTIMING_CLOCK_DOMAIN = "posix-clock-monotonic-ns-v1"
REPAIRTIMING_UNCERTAINTY = "paired-log-ratio-t95/v1"


@dataclass(frozen=True)
class RepairArm:
    """One exact provisional repair-v1 equation contract."""

    name: str
    provisional_id: int
    contract_sha256: str
    gf256_rows: int


REPAIR_V1_ARMS = (
    RepairArm(
        name="pure8_s0m1_d3_repair_v1",
        provisional_id=0x19cccf775ce0bf09,
        contract_sha256=(
            "19cccf775ce0bf098c9a425cb349714c"
            "4c4a880e7cf136c3bc365e13c05089a5"
        ),
        gf256_rows=8,
    ),
    RepairArm(
        name="pure9_s0m1_d3_repair_v1",
        provisional_id=0xa530f9105beaa450,
        contract_sha256=(
            "a530f9105beaa450dee70ad9b2a5cc54"
            "c944d3cd47f0aa6534630b8971608541"
        ),
        gf256_rows=9,
    ),
)
_REPAIR_ARM_BY_ID = {arm.provisional_id: arm for arm in REPAIR_V1_ARMS}
_REPAIR_ARM_BY_NAME = {arm.name: arm for arm in REPAIR_V1_ARMS}

REPAIRTIMING_PANEL_SPECS = (
    (
        "encoder_selector_forced", "encoder_selector",
        "repair_selector_encoder", "repair_forced_encoder",
    ),
    (
        "encoder_selector_aa", "encoder_selector",
        "repair_selector_encoder", "repair_selector_encoder",
    ),
    (
        "encoder_forced_aa", "encoder_forced",
        "repair_forced_encoder", "repair_forced_encoder",
    ),
    (
        "encoder_selector_wh1", "encoder_selector",
        "repair_selector_encoder", "wh1_encoder",
    ),
    (
        "encoder_wh1_aa", "encoder_selector",
        "wh1_encoder", "wh1_encoder",
    ),
    (
        "decoder_feed_wh1", "decoder_feed",
        "repair_decoder_feed", "wh1_decoder_feed",
    ),
    (
        "decoder_feed_aa", "decoder_feed",
        "repair_decoder_feed", "repair_decoder_feed",
    ),
    (
        "decoder_feed_wh1_aa", "decoder_feed",
        "wh1_decoder_feed", "wh1_decoder_feed",
    ),
    (
        "decoder_full_wh1", "decoder_full",
        "repair_decoder_full", "wh1_decoder_full",
    ),
    (
        "decoder_full_aa", "decoder_full",
        "repair_decoder_full", "repair_decoder_full",
    ),
    (
        "decoder_full_wh1_aa", "decoder_full",
        "wh1_decoder_full", "wh1_decoder_full",
    ),
    (
        "direct_dispatch", "direct",
        "repair_direct", "dispatch_direct",
    ),
    (
        "direct_aa", "direct",
        "repair_direct", "repair_direct",
    ),
    (
        "direct_dispatch_aa", "direct",
        "dispatch_direct", "dispatch_direct",
    ),
)
REPAIRTIMING_PANELS = tuple(spec[0] for spec in REPAIRTIMING_PANEL_SPECS)
REPAIRTIMING_SELECTOR_CALL_CAP = 1 + REPAIR_V1_ATTEMPT_CAP + 1
assert REPAIRTIMING_SELECTOR_CALL_CAP == 10
REPAIRTIMING_ROWS_PER_CELL = 1 + REPAIR_V1_ATTEMPT_CAP + (
    len(REPAIRTIMING_PANEL_SPECS) * len(REPAIRTIMING_ORDER)
)
assert REPAIRTIMING_ROWS_PER_CELL == 121

# Every stable field exposed by PrecodeSolveStats in a test-hook build.
REPAIRTIMING_STATS_FIELDS = (
    "packet_rows",
    "peeled_columns",
    "inactivated_columns",
    "residual_rows",
    "residual_rank",
    "binary_residual_rank",
    "binary_row_references",
    "binary_row_storage_bytes",
    "binary_adjacency_storage_bytes",
    "binary_row_storage_allocations",
    "binary_adjacency_storage_allocations",
    "block_xors",
    "block_muladds",
    "block_copies",
    "block_zero_fills",
    "block_add_sets",
    "block_add_set_sources",
    "peel_adjacency_visits",
    "peel_row_scan_steps",
    "peel_heap_operations",
    "projection_word_xors",
    "residual_coeff_word_xors",
    "residual_coeff_byte_ops",
    "build_ns",
    "peel_ns",
    "project_ns",
    "residual_ns",
    "backsub_ns",
    "packet_seed_attempt",
    "tiny_mixed_fastpath_acceptances",
    "mixed_joint_source_xors",
    "mixed_joint_marginal_xors",
    "mixed_joint_marginal_copies",
    "mixed_joint_scratch_bytes",
    "mixed_joint_active_deltas",
    "mixed_dual_source_columns",
)
assert len(REPAIRTIMING_STATS_FIELDS) == 36
REPAIRTIMING_STATS_TIME_FIELDS = (
    "build_ns", "peel_ns", "project_ns", "residual_ns", "backsub_ns",
)
REPAIRTIMING_STRUCTURAL_STATS_FIELDS = tuple(
    name for name in REPAIRTIMING_STATS_FIELDS
    if name not in REPAIRTIMING_STATS_TIME_FIELDS
)
assert len(REPAIRTIMING_STRUCTURAL_STATS_FIELDS) == 31

REPAIRTIMING_COMMON_FIELDS = (
    "row_kind",
    "cell_index",
    "measured",
    "construction_root",
    "loss_seed",
    "trace_sha256",
    "message_sha256",
)
REPAIRTIMING_SELECTOR_FIELDS = (
    "selector_result",
    "selected_attempt",
    "attempts_executed",
    "calls_executed",
    "real_configuration_calls",
    "structural_probe_calls",
    "cap_exhausted",
    "fatal_attempt_zero_mismatch",
    "oom",
    "committed",
    "descriptor_hex",
    "descriptor_sha256",
    "selector_structural_sha256",
)
REPAIRTIMING_ROLES = ("raw", "repaired", "dispatch", "wh1")
REPAIRTIMING_ROLE_SUFFIXES = (
    "encode_result",
    "decode_construct_result",
    "feed_result",
    "recover_result",
    "recovery_ok",
    "encoded_symbols",
    "received_symbols",
    "overhead",
    "payload_sha256",
    "recovered_sha256",
)
REPAIRTIMING_ROLE_FIELDS = tuple(
    f"{role}_{suffix}"
    for role in REPAIRTIMING_ROLES
    for suffix in REPAIRTIMING_ROLE_SUFFIXES
)
REPAIRTIMING_CONTROL_FIELDS = (
    "forced_a_result",
    "forced_b_result",
    "forced_a_payload_sha256",
    "forced_b_payload_sha256",
    "forced_equal",
    "direct_fixed_prefix_symbols",
    "repair_intermediate_bytes",
    "repair_intermediate_sha256",
    "repair_direct_executed",
    "repair_direct_result",
    "repair_direct_intermediate_bytes",
    "repair_direct_intermediate_sha256",
    "repair_direct_witness_equal",
    "dispatch_intermediate_bytes",
    "dispatch_intermediate_sha256",
    "dispatch_direct_executed",
    "dispatch_direct_result",
    "dispatch_direct_intermediate_bytes",
    "dispatch_direct_intermediate_sha256",
    "dispatch_direct_witness_equal",
)
REPAIRTIMING_ATTEMPT_FIELDS = (
    "attempt",
    "probe_executed",
    "probe_result",
    "probe_stats_available",
    *(f"probe_{name}" for name in REPAIRTIMING_STATS_FIELDS),
    "real_executed",
    "real_result",
    "real_stats_available",
    *(f"real_{name}" for name in REPAIRTIMING_STATS_FIELDS),
)
REPAIRTIMING_TIMING_FIELDS = (
    "timing_panel",
    "timing_panel_index",
    "timing_slot",
    "timing_pair",
    "timing_label",
    "timing_role",
    "timing_scope",
    "timing_eligible",
    "timing_executed",
    "timing_censor_reason",
    "timing_construct_result",
    "timing_result",
    "timing_recover_result",
    "timing_recovery_ok",
    "timing_stable",
    "timing_elapsed_ns",
    "timing_inner_reps",
    "timing_saturated",
    "timing_cpu_before",
    "timing_cpu_after",
    "timing_cpu_migrated",
    "timing_minflt",
    "timing_majflt",
    "timing_fault_contaminated",
    "timing_stats_available",
    *(f"timing_{name}" for name in REPAIRTIMING_STATS_FIELDS),
    "timing_fixed_prefix_symbols",
    "timing_intermediate_bytes",
    "timing_intermediate_sha256",
)
REPAIRTIMING_SUMMARY_COLUMNS = (
    REPAIRTIMING_COMMON_FIELDS +
    REPAIRTIMING_SELECTOR_FIELDS +
    REPAIRTIMING_ROLE_FIELDS +
    REPAIRTIMING_CONTROL_FIELDS
)
REPAIRTIMING_ATTEMPT_COLUMNS = (
    REPAIRTIMING_COMMON_FIELDS + REPAIRTIMING_ATTEMPT_FIELDS
)
REPAIRTIMING_TIMING_COLUMNS = (
    REPAIRTIMING_COMMON_FIELDS + REPAIRTIMING_TIMING_FIELDS
)
assert len(REPAIRTIMING_SUMMARY_COLUMNS) == 80
assert len(REPAIRTIMING_ATTEMPT_COLUMNS) == 86
assert len(REPAIRTIMING_TIMING_COLUMNS) == 71
REPAIRTIMING_COLUMNS = tuple(dict.fromkeys(
    REPAIRTIMING_SUMMARY_COLUMNS +
    REPAIRTIMING_ATTEMPT_COLUMNS +
    REPAIRTIMING_TIMING_COLUMNS
))
assert len(REPAIRTIMING_COLUMNS) == len(set(REPAIRTIMING_COLUMNS))

REPAIRTIMING_MANIFEST_FIELDS = (
    "protocol",
    "schema",
    "repair_arm",
    "repair_contract_id",
    "repair_contract_sha256",
    "repair_policy_sha256",
    "repair_attempt_cap",
    "descriptor_magic",
    "descriptor_encoding",
    "descriptor_bytes",
    "descriptor_policy",
    "K",
    "bb",
    "message_bytes",
    "message_tail_bytes",
    "dispatch_profile",
    "dispatch_contract_id",
    "construction_seed_base",
    "construction_seed_derivation",
    "loss",
    "loss_seed_base",
    "loss_seed_derivation",
    "schedule",
    "loss_model",
    "trace_encoding",
    "message_seed_policy",
    "warmup_replicates",
    "replicates",
    "cells",
    "summary_rows_per_cell",
    "attempt_rows_per_cell",
    "timing_panels",
    "timing_panels_per_cell",
    "timing_slots_per_panel",
    "timing_rows_per_cell",
    "timing_order",
    "timing_label_policy",
    "inner_reps",
    "max_overhead",
    "cache_state",
    "systematic_cache",
    "evict_bytes",
    "payload_alignment",
    "prefault",
    "cpu_affinity_policy",
    "encoder_selector_scope",
    "encoder_forced_scope",
    "decoder_feed_scope",
    "decoder_full_scope",
    "direct_scope",
    "codec_reuse",
    "scratch_reporting",
    "timing_stats_policy",
    "selector_timing_stats_policy",
    "structural_digest_encoding",
    "context_sha256",
    "uncertainty",
    "required_margin",
    "margin_rule",
    "clock_domain",
    "stream_hash_scope",
    "started_monotonic_ns",
    "expected_summary_rows",
    "expected_attempt_rows",
    "expected_timing_rows",
    "expected_rows",
)

_RESULT_CODES = frozenset(range(11))
_DIGEST_FIELDS = frozenset({
    "trace_sha256",
    "message_sha256",
    "descriptor_sha256",
    "selector_structural_sha256",
    *(f"{role}_payload_sha256" for role in REPAIRTIMING_ROLES),
    *(f"{role}_recovered_sha256" for role in REPAIRTIMING_ROLES),
    "forced_a_payload_sha256",
    "forced_b_payload_sha256",
    "repair_intermediate_sha256",
    "repair_direct_intermediate_sha256",
    "dispatch_intermediate_sha256",
    "dispatch_direct_intermediate_sha256",
    "timing_intermediate_sha256",
})
_STRING_FIELDS = frozenset({
    "row_kind",
    "descriptor_hex",
    "timing_panel",
    "timing_label",
    "timing_role",
    "timing_scope",
    "timing_censor_reason",
    *_DIGEST_FIELDS,
})
_NUMERIC_FIELDS = tuple(
    name for name in REPAIRTIMING_COLUMNS if name not in _STRING_FIELDS
)

# Native streams are bounded before splitting or CSV parsing.  These maxima
# are derived from the exact union schema and canonical field spellings rather
# than copied from an older timing protocol.
REPAIRTIMING_MAX_CELLS = 768
REPAIRTIMING_HEADER_BYTE_LIMIT = 64 * 1024
REPAIRTIMING_MANIFEST_BYTE_LIMIT = 64 * 1024


def _wire_field_byte_limit(name):
    if name in _NUMERIC_FIELDS:
        # uint64 max is 20 digits; the only negative wire value is "-1".
        return 20
    if name in _DIGEST_FIELDS:
        return 64
    if name == "descriptor_hex":
        return 2 * REPAIR_V1_DESCRIPTOR_BYTES
    if name == "row_kind":
        return len("summary")
    if name == "timing_panel":
        return max(len(value) for value in REPAIRTIMING_PANELS)
    if name == "timing_label":
        return len("not_applicable")
    if name == "timing_role":
        return max(
            len("not_applicable"),
            *(len(spec[2]) for spec in REPAIRTIMING_PANEL_SPECS),
            *(len(spec[3]) for spec in REPAIRTIMING_PANEL_SPECS),
        )
    if name == "timing_scope":
        return max(
            len("not_applicable"),
            *(len(spec[1]) for spec in REPAIRTIMING_PANEL_SPECS),
        )
    if name == "timing_censor_reason":
        return max(
            len(value) for value in (
                "not_applicable",
                "none",
                "selector_uncommitted",
                "role_preflight_failed",
                "cpu_migrated",
                "fault_contaminated",
                "saturated",
                "outcome_unstable",
            )
        )
    raise AssertionError(f"unbounded repairtiming field {name}")


REPAIRTIMING_FIELD_BYTE_LIMITS = {
    name: _wire_field_byte_limit(name)
    for name in REPAIRTIMING_COLUMNS
}
REPAIRTIMING_FIELD_BYTE_LIMIT = max(
    REPAIRTIMING_FIELD_BYTE_LIMITS.values())
def _columns_byte_limit(columns):
    return (
        sum(REPAIRTIMING_FIELD_BYTE_LIMITS[name] for name in columns) +
        len(columns)  # commas plus final LF
    )


REPAIRTIMING_SUMMARY_ROW_BYTE_LIMIT = _columns_byte_limit(
    REPAIRTIMING_SUMMARY_COLUMNS)
REPAIRTIMING_ATTEMPT_ROW_BYTE_LIMIT = _columns_byte_limit(
    REPAIRTIMING_ATTEMPT_COLUMNS)
REPAIRTIMING_TIMING_ROW_BYTE_LIMIT = _columns_byte_limit(
    REPAIRTIMING_TIMING_COLUMNS)
REPAIRTIMING_ROW_BYTE_LIMIT = max(
    REPAIRTIMING_SUMMARY_ROW_BYTE_LIMIT,
    REPAIRTIMING_ATTEMPT_ROW_BYTE_LIMIT,
    REPAIRTIMING_TIMING_ROW_BYTE_LIMIT,
)
assert REPAIRTIMING_ROW_BYTE_LIMIT < 32 * 1024
REPAIRTIMING_FIXED_TEXT_BYTE_LIMIT = (
    REPAIRTIMING_MANIFEST_BYTE_LIMIT +
    REPAIRTIMING_HEADER_BYTE_LIMIT +
    1024
)
REPAIRTIMING_CELL_TEXT_BYTE_LIMIT = (
    REPAIRTIMING_SUMMARY_ROW_BYTE_LIMIT +
    REPAIR_V1_ATTEMPT_CAP * REPAIRTIMING_ATTEMPT_ROW_BYTE_LIMIT +
    len(REPAIRTIMING_PANELS) * len(REPAIRTIMING_ORDER) *
        REPAIRTIMING_TIMING_ROW_BYTE_LIMIT
)
REPAIRTIMING_NATIVE_STDOUT_BYTE_LIMIT = (
    REPAIRTIMING_FIXED_TEXT_BYTE_LIMIT +
    REPAIRTIMING_MAX_CELLS * REPAIRTIMING_CELL_TEXT_BYTE_LIMIT
)
REPAIRTIMING_STDOUT_BYTE_LIMIT = REPAIRTIMING_NATIVE_STDOUT_BYTE_LIMIT
REPAIRTIMING_RECEIPT_CELL_CANONICAL_BYTE_LIMIT = (
    REPAIRTIMING_CELL_TEXT_BYTE_LIMIT +
    REPAIRTIMING_ROWS_PER_CELL +  # JSON escapes each LF as two bytes
    4096  # Per-cell canonical-object punctuation/headroom.
)
assert REPAIRTIMING_RECEIPT_CELL_CANONICAL_BYTE_LIMIT <= 256 * 1024
REPAIRTIMING_CONTEXT_CANONICAL_BYTE_LIMIT = 1 * 1024 * 1024
REPAIRTIMING_RECEIPT_FIXED_CANONICAL_BYTE_LIMIT = (
    REPAIRTIMING_CONTEXT_CANONICAL_BYTE_LIMIT +
    REPAIRTIMING_FIXED_TEXT_BYTE_LIMIT +
    512 * 1024
)
assert REPAIRTIMING_RECEIPT_FIXED_CANONICAL_BYTE_LIMIT <= 2 * 1024 * 1024
REPAIRTIMING_RECEIPT_CANONICAL_BYTE_LIMIT = (
    REPAIRTIMING_RECEIPT_FIXED_CANONICAL_BYTE_LIMIT +
    REPAIRTIMING_MAX_CELLS *
        REPAIRTIMING_RECEIPT_CELL_CANONICAL_BYTE_LIMIT
)
REPAIRTIMING_GZIP_COMPRESSED_BYTE_LIMIT = 512 * 1024 * 1024
REPAIRTIMING_GZIP_DECOMPRESSED_BYTE_LIMIT = (
    REPAIRTIMING_NATIVE_STDOUT_BYTE_LIMIT
)
REPAIRTIMING_STDERR_BYTE_LIMIT = 64 * 1024
REPAIRTIMING_DIAGNOSTIC_BYTE_LIMIT = 4096
REPAIRTIMING_CAPTURE_CHUNK_BYTES = 64 * 1024
REPAIRTIMING_REAP_TIMEOUT_SECONDS = 2.0
REPAIRTIMING_TERMINATE_GRACE_SECONDS = 0.1
REPAIRTIMING_GROUP_CONFIRM_SECONDS = 2.0
REPAIRTIMING_SELECT_POLL_SECONDS = 0.05
REPAIRTIMING_POST_EXIT_DRAIN_SECONDS = 0.5
_REPAIRTIMING_PR_SET_PDEATHSIG = 1

if sys.platform.startswith("linux"):
    _REPAIRTIMING_LIBC = ctypes.CDLL(None, use_errno=True)
    _REPAIRTIMING_PRCTL = getattr(_REPAIRTIMING_LIBC, "prctl", None)
    if _REPAIRTIMING_PRCTL is not None:
        _REPAIRTIMING_PRCTL.argtypes = (
            ctypes.c_int,
            ctypes.c_ulong,
            ctypes.c_ulong,
            ctypes.c_ulong,
            ctypes.c_ulong,
        )
        _REPAIRTIMING_PRCTL.restype = ctypes.c_int
else:
    _REPAIRTIMING_PRCTL = None


def repair_arm(value):
    """Resolve an exact provisional arm by object, name, or integer id."""
    if type(value) is RepairArm and value in REPAIR_V1_ARMS:
        return value
    if type(value) is str and value in _REPAIR_ARM_BY_NAME:
        return _REPAIR_ARM_BY_NAME[value]
    if type(value) is int and value in _REPAIR_ARM_BY_ID:
        return _REPAIR_ARM_BY_ID[value]
    raise ValueError("unknown repair-v1 arm")


def _expected_construction_root(base, replicate, derivation):
    if derivation == REPAIRTIMING_CONSTRUCTION_SEED_FIXED:
        return base
    if derivation == REPAIRTIMING_CONSTRUCTION_SEED_DERIVED:
        return (
            base ^ (((replicate + 1) * 0xd192ed03) & 0xffffffff)
        ) & 0xffffffff
    raise ValueError("unknown construction seed derivation")


def _expected_loss_seed(base, replicate):
    return (
        base ^
        (((replicate + 1) * 0x9e3779b97f4a7c15) &
         0xffffffffffffffff)
    ) & 0xffffffffffffffff


def _message_bytes(block_count, block_bytes):
    return (block_count - 1) * block_bytes + 1


def _candidate_staircase(block_count):
    staircase = peel_codec._dispatch_staircase_count(block_count) - 1
    if staircase < 1:
        raise ValueError("repair-v1 staircase is invalid")
    return staircase


def _candidate_intermediate_bytes(arm, block_count, block_bytes):
    width = (
        block_count + _candidate_staircase(block_count) +
        3 + arm.gf256_rows
    )
    return width * block_bytes


def _dispatch_intermediate_bytes(block_count, block_bytes):
    dispatch = peel_codec.TARGET_CONTRACT
    width = (
        block_count + peel_codec._dispatch_staircase_count(block_count) +
        dispatch["dense_rows"] + dispatch["heavy_rows"]
    )
    return width * block_bytes


def validate_repairtiming_dimensions(
        *, block_count, block_bytes, repair_arm,
        repair_policy=REPAIR_V1_POLICY_NAME,
        dispatch_profile=peel_codec.TARGET_PROFILE,
        construction_seed, construction_seed_derivation,
        loss, loss_seed, schedule, warmup_replicates, replicates, inner_reps,
        max_overhead, cache_state, systematic_cache, evict_bytes,
        required_margin):
    """Mirror the finite request domain of native ``repairtiming``."""
    arm = globals()["repair_arm"](repair_arm)
    del arm
    try:
        numeric_loss = (
            float(loss)
            if not isinstance(loss, bool) and isinstance(loss, (int, float))
            else math.nan
        )
        numeric_margin = (
            float(required_margin)
            if not isinstance(required_margin, bool) and
            isinstance(required_margin, (int, float))
            else math.nan
        )
    except OverflowError:
        numeric_loss = numeric_margin = math.nan
    cells = (
        warmup_replicates + replicates
        if type(warmup_replicates) is int and type(replicates) is int
        else -1
    )
    if (
        repair_policy != REPAIR_V1_POLICY_NAME or
        dispatch_profile != peel_codec.TARGET_PROFILE or
        type(block_count) is not int or not 2 <= block_count <= 100 or
        type(block_bytes) is not int or
        not 2 <= block_bytes <= 4096 or block_bytes % 2 != 0 or
        type(construction_seed) is not int or
        not 0 <= construction_seed <= 0xffffffff or
        type(construction_seed_derivation) is not str or
        construction_seed_derivation not in
            REPAIRTIMING_CONSTRUCTION_SEED_DERIVATIONS or
        type(loss_seed) is not int or
        not 0 <= loss_seed <= 0xffffffffffffffff or
        not math.isfinite(numeric_loss) or
        (numeric_loss == 0.0 and
         math.copysign(1.0, numeric_loss) < 0.0) or
        not 0.0 <= numeric_loss <= 0.99 or
        type(schedule) is not str or
        schedule not in REPAIRTIMING_SCHEDULES or
        type(warmup_replicates) is not int or warmup_replicates < 0 or
        type(replicates) is not int or not 3 <= replicates <= 768 or
        not 1 <= cells <= REPAIRTIMING_MAX_CELLS or
        type(inner_reps) is not int or not 1 <= inner_reps <= 1024 or
        type(max_overhead) is not int or
        max_overhead != 64 or
        cache_state not in ("warm", "cold") or
        systematic_cache not in ("off", "on") or
        type(evict_bytes) is not int or
        not 4096 <= evict_bytes <= peel_codec.PEELTIMING_EVICT_BYTES_MAX or
        (cache_state == "cold" and inner_reps != 1) or
        not math.isfinite(numeric_margin) or
        (numeric_margin == 0.0 and
         math.copysign(1.0, numeric_margin) < 0.0) or
        not 0.0 <= numeric_margin <= 1.0
    ):
        raise ValueError("invalid repairtiming dimensions or policy")


def _parse_decimal(text, label, *, minimum=-1,
                   maximum=0xffffffffffffffff):
    return peel_codec._parse_decimal_integer(
        text, label, minimum=minimum, maximum=maximum)


def _is_digest_or_na(value):
    return value == "not_applicable" or peel_codec._is_sha256(value)


def _canonical_csv_fields(line, label, *, columns=None):
    try:
        fields = next(csv.reader([line], strict=True))
    except (csv.Error, StopIteration) as error:
        raise MeasurementError(f"malformed {label}: {error}")
    if line != ",".join(fields):
        raise MeasurementError(f"{label} is not canonical CSV")
    if columns is None:
        if any(len(field.encode("ascii")) > REPAIRTIMING_FIELD_BYTE_LIMIT
               for field in fields):
            raise MeasurementError(f"{label} contains an oversized field")
    elif any(
            len(field.encode("ascii")) >
                REPAIRTIMING_FIELD_BYTE_LIMITS[name]
            for name, field in zip(columns, fields)
    ):
        raise MeasurementError(f"{label} contains an oversized field")
    return fields


def _parse_typed_row(line, row_index, columns):
    fields = _canonical_csv_fields(
        line, f"repairtiming row {row_index}",
        columns=columns)
    if len(fields) != len(columns):
        raise MeasurementError(
            f"repairtiming row {row_index} has the wrong column count")
    raw = dict(zip(columns, fields))
    row = {}
    for name in columns:
        if name in _STRING_FIELDS:
            row[name] = raw[name]
        else:
            row[name] = _parse_decimal(
                raw[name], f"repairtiming row {row_index}.{name}")
    for name in _DIGEST_FIELDS.intersection(columns):
        if not _is_digest_or_na(row[name]):
            raise MeasurementError(
                f"repairtiming row {row_index}.{name} is invalid")
    return row


def _stats_from_row(row, prefix):
    return {
        name: row[f"{prefix}_{name}"]
        for name in REPAIRTIMING_STATS_FIELDS
    }


def _validate_stats(row, prefix, available, label, *,
                    elapsed_ns: Optional[int] = None,
                    acceptance_max: Optional[int] = None,
                    call_limit: int = 1):
    stats = _stats_from_row(row, prefix)
    if (
        type(call_limit) is not int or
        not 1 <= call_limit <= REPAIRTIMING_SELECTOR_CALL_CAP
    ):
        raise AssertionError("invalid internal repairtiming call limit")
    if available not in (0, 1):
        raise MeasurementError(f"{label} has invalid stats availability")
    if not available:
        if any(value != -1 for value in stats.values()):
            raise MeasurementError(f"{label} fabricated unavailable stats")
        return stats
    if any(value < 0 for value in stats.values()):
        raise MeasurementError(f"{label} omitted available stats")
    if (
        stats["residual_rank"] > stats["residual_rows"] or
        stats["binary_residual_rank"] > stats["residual_rank"] or
        (acceptance_max is not None and
         stats["tiny_mixed_fastpath_acceptances"] > acceptance_max) or
        stats["inactivated_columns"] > 0xffff * call_limit or
        stats["packet_seed_attempt"] > 255 * call_limit or
        stats["mixed_joint_scratch_bytes"] >
            peel_codec.PEELTIMING_WORKING_SET_BYTES_MAX * call_limit or
        stats["binary_row_storage_bytes"] >
            peel_codec.PEELTIMING_WORKING_SET_BYTES_MAX * call_limit or
        stats["binary_adjacency_storage_bytes"] >
            peel_codec.PEELTIMING_WORKING_SET_BYTES_MAX * call_limit or
        stats["mixed_joint_active_deltas"] >
            stats["inactivated_columns"] or
        stats["binary_row_storage_bytes"] == 0 and
            stats["binary_row_storage_allocations"] != 0 or
        stats["binary_adjacency_storage_bytes"] == 0 and
            stats["binary_adjacency_storage_allocations"] != 0 or
        stats["binary_row_storage_allocations"] >
            stats["binary_row_storage_bytes"] or
        stats["binary_adjacency_storage_allocations"] >
            stats["binary_adjacency_storage_bytes"] or
        stats["block_add_sets"] > stats["block_add_set_sources"]
    ):
        raise MeasurementError(f"{label} has impossible solver counters")
    phase_ns = sum(stats[name] for name in (
        "build_ns", "peel_ns", "project_ns", "residual_ns", "backsub_ns"))
    if elapsed_ns is not None and phase_ns > elapsed_ns:
        raise MeasurementError(f"{label} phase time exceeds elapsed time")
    return stats


def _descriptor_bytes(arm, construction_root, selected_attempt):
    return (
        REPAIR_V1_DESCRIPTOR_MAGIC.encode("ascii") +
        struct.pack("<QI", arm.provisional_id, construction_root) +
        bytes((selected_attempt,))
    )


def _result_code(value, label, *, allow_na=False):
    if value == -1 and allow_na:
        return value
    if value not in _RESULT_CODES:
        raise MeasurementError(f"{label} has an invalid Wirehair result")
    return value


def _result_class(value):
    if value == -1:
        return "not_applicable"
    if value == 0:
        return "success"
    if value in (1, 7):
        return "need_more"
    if value in (3, 4):
        return "weak"
    if value in _RESULT_CODES:
        return "error"
    raise MeasurementError(f"unknown Wirehair result code {value}")


def _json_safe_copy(value):
    return copy.deepcopy(value)


@dataclass(frozen=True)
class RepairTimingCell:
    """Lossless normalized evidence for one physical recovery cell."""

    replicate: int
    construction_root: int
    loss_seed: int
    selector_key: tuple
    selector_witness: dict
    real_witness: dict

    def as_dict(self):
        return {
            "replicate": self.replicate,
            "construction_root": self.construction_root,
            "loss_seed": self.loss_seed,
            "selector_key": list(self.selector_key),
            "selector_witness": _json_safe_copy(self.selector_witness),
            "real_witness": _json_safe_copy(self.real_witness),
        }


@dataclass(frozen=True)
class RepairTimingEvidence:
    """A fully validated native v3 stream and its normalized evidence."""

    protocol: str
    manifest: dict
    context: dict
    cells: tuple
    stream_sha256: str
    evidence: dict
    native_stdout: str

    def as_dict(self):
        return {
            "protocol": self.protocol,
            "manifest": _json_safe_copy(self.manifest),
            "context": _json_safe_copy(self.context),
            "stream_sha256": self.stream_sha256,
            "evidence": _json_safe_copy(self.evidence),
            "native_stdout": self.native_stdout,
        }


def selector_projection(cell):
    """Return the exact payload-independent selector-dedup projection."""
    if not isinstance(cell, RepairTimingCell):
        raise TypeError("selector projection requires RepairTimingCell")
    return _json_safe_copy(cell.selector_witness)


def real_projection(cell):
    """Return the exact loss/payload/timing projection for aggregation."""
    if not isinstance(cell, RepairTimingCell):
        raise TypeError("real projection requires RepairTimingCell")
    return _json_safe_copy(cell.real_witness)


def _parse_manifest(
        raw, *, arm, block_count, block_bytes, dispatch_profile,
        construction_seed, construction_seed_derivation, loss, loss_seed,
        schedule, warmup_replicates, replicates, inner_reps, max_overhead,
        cache_state, systematic_cache, evict_bytes, context_sha256,
        required_margin):
    integer_fields = {
        "repair_attempt_cap",
        "descriptor_bytes",
        "K",
        "bb",
        "message_bytes",
        "message_tail_bytes",
        "construction_seed_base",
        "loss_seed_base",
        "warmup_replicates",
        "replicates",
        "cells",
        "summary_rows_per_cell",
        "attempt_rows_per_cell",
        "timing_panels_per_cell",
        "timing_slots_per_panel",
        "timing_rows_per_cell",
        "inner_reps",
        "max_overhead",
        "evict_bytes",
        "payload_alignment",
        "prefault",
        "started_monotonic_ns",
        "expected_summary_rows",
        "expected_attempt_rows",
        "expected_timing_rows",
        "expected_rows",
    }
    manifest = dict(raw)
    for name in integer_fields:
        manifest[name] = _parse_decimal(
            raw[name], f"repairtiming manifest.{name}", minimum=0)
    manifest["loss"] = peel_codec._parse_finite_float(
        raw["loss"], "repairtiming manifest.loss",
        minimum=0.0, maximum=0.99)
    manifest["required_margin"] = peel_codec._parse_finite_float(
        raw["required_margin"], "repairtiming manifest.required_margin",
        minimum=0.0, maximum=1.0)
    cells = warmup_replicates + replicates
    expected = {
        "protocol": REPAIRTIMING_PROTOCOL_V3,
        "schema": REPAIRTIMING_SCHEMA_V3,
        "repair_arm": arm.name,
        "repair_contract_id": f"{arm.provisional_id:016x}",
        "repair_contract_sha256": arm.contract_sha256,
        "repair_policy_sha256": REPAIR_V1_POLICY_SHA256,
        "repair_attempt_cap": REPAIR_V1_ATTEMPT_CAP,
        "descriptor_magic": REPAIR_V1_DESCRIPTOR_MAGIC,
        "descriptor_encoding": REPAIR_V1_DESCRIPTOR_ENCODING,
        "descriptor_bytes": REPAIR_V1_DESCRIPTOR_BYTES,
        "descriptor_policy": REPAIR_V1_DESCRIPTOR_POLICY,
        "K": block_count,
        "bb": block_bytes,
        "message_bytes": _message_bytes(block_count, block_bytes),
        "message_tail_bytes": 1,
        "dispatch_profile": dispatch_profile,
        "dispatch_contract_id":
            peel_codec.TARGET_CONTRACT["contract_id"],
        "construction_seed_base": construction_seed,
        "construction_seed_derivation": construction_seed_derivation,
        "loss": float(loss),
        "loss_seed_base": loss_seed,
        "loss_seed_derivation": REPAIRTIMING_LOSS_SEED_DERIVATION,
        "schedule": schedule,
        "loss_model": "packet-schedule-v1",
        "trace_encoding":
            "wirehair-wh2-repairtiming-loss-trace-v1",
        "message_seed_policy":
            "replicate-loss-seed-partial-final-v1",
        "warmup_replicates": warmup_replicates,
        "replicates": replicates,
        "cells": cells,
        "summary_rows_per_cell": 1,
        "attempt_rows_per_cell": REPAIR_V1_ATTEMPT_CAP,
        "timing_panels": "+".join(REPAIRTIMING_PANELS),
        "timing_panels_per_cell": len(REPAIRTIMING_PANELS),
        "timing_slots_per_panel": len(REPAIRTIMING_ORDER),
        "timing_rows_per_cell":
            len(REPAIRTIMING_PANELS) * len(REPAIRTIMING_ORDER),
        "timing_order": REPAIRTIMING_ORDER,
        "timing_label_policy": "fixed-logical-role-v1",
        "inner_reps": inner_reps,
        "max_overhead": max_overhead,
        "cache_state": cache_state,
        "systematic_cache": systematic_cache,
        "evict_bytes": evict_bytes,
        "payload_alignment": 64,
        "prefault": 1,
        "cpu_affinity_policy": "first-allowed-affinity-v1",
        "encoder_selector_scope":
            "lazy-selector-init-first-k-encode-w2r3-serialize-v1",
        "encoder_forced_scope":
            "forced-selected-init-first-k-encode-w2r3-serialize-v1",
        "decoder_feed_scope":
            "exact-selected-descriptor-bound-feed-through-first-success-v1",
        "decoder_full_scope":
            "w2r3-parse-create-feed-through-first-success-recover-v1",
        "direct_scope":
            "selected-attempt-pair-local-fixed-prefix-encoder-"
            "intermediate-witnessed-solve-v1",
        "codec_reuse": "none-fresh-object-every-inner-v1",
        "scratch_reporting":
            "complete-precode-solve-stats-per-call-v1",
        "timing_stats_policy":
            "first-inner-structural-counters-summed-phase-ns-v1",
        "selector_timing_stats_policy":
            "all-executed-selector-subcalls-summed-per-inner-v1",
        "structural_digest_encoding":
            "wirehair-wh2-repair-v1-selector-structural-v1",
        "context_sha256": context_sha256,
        "uncertainty": REPAIRTIMING_UNCERTAINTY,
        "required_margin": float(required_margin),
        "margin_rule":
            "upper-log-cost-lt-negative-required-margin-and-aa-floors-v1",
        "clock_domain": REPAIRTIMING_CLOCK_DOMAIN,
        "stream_hash_scope": "body-plus-done-prefix-v1",
        "expected_summary_rows": cells,
        "expected_attempt_rows": cells * REPAIR_V1_ATTEMPT_CAP,
        "expected_timing_rows":
            cells * len(REPAIRTIMING_PANELS) *
            len(REPAIRTIMING_ORDER),
        "expected_rows": cells * REPAIRTIMING_ROWS_PER_CELL,
    }
    for name, expected_value in expected.items():
        if manifest[name] != expected_value:
            raise MeasurementError(
                f"repairtiming manifest.{name} does not match the request")
    if raw["loss"] != f"{float(loss):.17g}":
        raise MeasurementError("repairtiming loss spelling is noncanonical")
    if raw["required_margin"] != f"{float(required_margin):.17g}":
        raise MeasurementError(
            "repairtiming required margin spelling is noncanonical")
    if (manifest["loss"] > 0.99 or
            manifest["started_monotonic_ns"] < 1):
        raise MeasurementError("repairtiming manifest domain is invalid")
    return manifest


def _read_line(stream, label, *, byte_limit):
    line = stream.readline(byte_limit + 2)
    if not line:
        raise MeasurementError(f"repairtiming is truncated before {label}")
    if (len(line) > byte_limit or not line.endswith("\n") or
            line.endswith("\r\n") or line == "\n" or "\0" in line):
        raise MeasurementError(
            f"repairtiming {label} is oversized or noncanonical")
    return line


def _validate_common(
        row, *, row_kind, cell_index, measured, construction_root,
        loss_seed, common_hashes, label):
    if (
        row["row_kind"] != row_kind or
        row["cell_index"] != cell_index or
        row["measured"] != measured or
        row["construction_root"] != construction_root or
        row["loss_seed"] != loss_seed or
        not peel_codec._is_sha256(row["trace_sha256"]) or
        not peel_codec._is_sha256(row["message_sha256"])
    ):
        raise MeasurementError(f"{label} has invalid common fields")
    hashes = (row["trace_sha256"], row["message_sha256"])
    if common_hashes is not None and hashes != common_hashes:
        raise MeasurementError(f"{label} changed trace or message identity")
    return hashes


def _summary_numeric_fields():
    return tuple(
        name for name in (
            REPAIRTIMING_SELECTOR_FIELDS +
            REPAIRTIMING_ROLE_FIELDS +
            REPAIRTIMING_CONTROL_FIELDS
        )
        if name not in _STRING_FIELDS
    )


def _summary_string_fields():
    return tuple(
        name for name in (
            REPAIRTIMING_SELECTOR_FIELDS +
            REPAIRTIMING_ROLE_FIELDS +
            REPAIRTIMING_CONTROL_FIELDS
        )
        if name in _STRING_FIELDS
    )


_SUMMARY_NUMERIC_FIELDS = _summary_numeric_fields()
_SUMMARY_STRING_FIELDS = _summary_string_fields()
_TIMING_NUMERIC_FIELDS = tuple(
    name for name in REPAIRTIMING_TIMING_FIELDS
    if name not in _STRING_FIELDS
)
_TIMING_STRING_FIELDS = tuple(
    name for name in REPAIRTIMING_TIMING_FIELDS
    if name in _STRING_FIELDS
)


def _parse_role_summary(row, role, *, block_count, max_overhead,
                        message_sha256, label):
    values = {
        suffix: row[f"{role}_{suffix}"]
        for suffix in REPAIRTIMING_ROLE_SUFFIXES
    }
    encode = _result_code(
        values["encode_result"], f"{label}.{role}_encode_result")
    decode_construct = _result_code(
        values["decode_construct_result"],
        f"{label}.{role}_decode_construct_result", allow_na=True)
    feed = _result_code(
        values["feed_result"], f"{label}.{role}_feed_result",
        allow_na=True)
    recover = _result_code(
        values["recover_result"], f"{label}.{role}_recover_result",
        allow_na=True)
    if (
        values["recovery_ok"] not in (0, 1) or
        values["encoded_symbols"] < 0 or
        values["received_symbols"] < 0 or
        not _is_digest_or_na(values["payload_sha256"]) or
        not _is_digest_or_na(values["recovered_sha256"])
    ):
        raise MeasurementError(f"{label}.{role} recovery receipt is invalid")

    if encode != 0:
        valid = (
            decode_construct == -1 and feed == -1 and recover == -1 and
            values["recovery_ok"] == 0 and
            values["encoded_symbols"] == 0 and
            values["received_symbols"] == 0 and
            values["overhead"] == -1 and
            values["payload_sha256"] == "not_applicable" and
            values["recovered_sha256"] == "not_applicable"
        )
    else:
        recovered = (
            decode_construct == 0 and feed == 0 and recover == 0 and
            values["recovery_ok"] == 1
        )
        valid = (
            values["encoded_symbols"] == block_count + max_overhead and
            peel_codec._is_sha256(values["payload_sha256"]) and
            decode_construct != -1 and
            (
                (
                    recovered and
                    block_count <= values["received_symbols"] <=
                        block_count + max_overhead and
                    values["overhead"] ==
                        values["received_symbols"] - block_count and
                    values["recovered_sha256"] == message_sha256
                ) or
                (
                    not recovered and
                    values["recovery_ok"] == 0 and
                    values["overhead"] == -1 and
                    values["recovered_sha256"] == "not_applicable" and
                    0 <= values["received_symbols"] <=
                        block_count + max_overhead
                )
            )
        )
        if decode_construct != 0:
            valid = (
                valid and feed == -1 and recover == -1 and
                values["received_symbols"] == 0
            )
        elif feed != 0:
            valid = (
                valid and feed != -1 and recover == -1 and
                values["received_symbols"] >= 1
            )
        else:
            valid = (
                valid and recover != -1 and
                values["received_symbols"] >= block_count and
                (recover == 0) == (values["recovery_ok"] == 1)
            )
    if not valid:
        raise MeasurementError(
            f"{label}.{role} has an impossible recovery sequence")
    values["encode_class"] = _result_class(encode)
    values["feed_class"] = _result_class(feed)
    values["recover_class"] = _result_class(recover)
    outcome = next(
        (
            result for result in (
                encode, decode_construct, feed, recover)
            if result not in (-1, 0)
        ),
        0 if values["recovery_ok"] else 8,
    )
    values["outcome_class"] = _result_class(outcome)
    return values


def _validate_summary_row(
        row, *, arm, block_count, block_bytes, max_overhead, label):
    selector_result = _result_code(
        row["selector_result"], f"{label}.selector_result")
    selected = row["selected_attempt"]
    if selected != -1 and not 0 <= selected < REPAIR_V1_ATTEMPT_CAP:
        raise MeasurementError(f"{label}.selected_attempt is invalid")
    for name in (
            "cap_exhausted", "fatal_attempt_zero_mismatch", "oom",
            "committed", "forced_equal", "repair_direct_witness_equal",
            "dispatch_direct_witness_equal"):
        if row[name] not in (0, 1):
            raise MeasurementError(f"{label}.{name} is not boolean")
    if (
        not 1 <= row["attempts_executed"] <= REPAIR_V1_ATTEMPT_CAP or
        not 1 <= row["calls_executed"] <=
            REPAIRTIMING_SELECTOR_CALL_CAP or
        not 1 <= row["real_configuration_calls"] <= 2 or
        not 0 <= row["structural_probe_calls"] <=
            REPAIR_V1_ATTEMPT_CAP or
        row["calls_executed"] !=
            row["real_configuration_calls"] +
            row["structural_probe_calls"] or
        row["oom"] != int(selector_result == 9) or
        row["committed"] != int(selector_result == 0) or
        not peel_codec._is_sha256(row["selector_structural_sha256"])
    ):
        raise MeasurementError(f"{label} selector accounting is invalid")

    if row["committed"]:
        if selected < 0:
            raise MeasurementError(
                f"{label} committed without a selected attempt")
        expected_descriptor = _descriptor_bytes(
            arm, row["construction_root"], selected)
        if (
            row["descriptor_hex"] != expected_descriptor.hex() or
            row["descriptor_sha256"] !=
                hashlib.sha256(expected_descriptor).hexdigest()
        ):
            raise MeasurementError(f"{label} W2R3 descriptor is invalid")
    elif (
        row["descriptor_hex"] != "not_applicable" or
        row["descriptor_sha256"] != "not_applicable"
    ):
        raise MeasurementError(
            f"{label} serialized an uncommitted W2R3 descriptor")

    roles = {
        role: _parse_role_summary(
            row, role, block_count=block_count,
            max_overhead=max_overhead,
            message_sha256=row["message_sha256"], label=label)
        for role in REPAIRTIMING_ROLES
    }
    if roles["repaired"]["encode_result"] != selector_result:
        raise MeasurementError(
            f"{label} repaired result disagrees with selector")
    if (
        row["committed"] and selected == 0 and
        any(
            roles["raw"][suffix] != roles["repaired"][suffix]
            for suffix in REPAIRTIMING_ROLE_SUFFIXES
        )
    ):
        raise MeasurementError(
            f"{label} attempt-zero selector changed raw work")
    if row["committed"]:
        if (
            row["forced_a_result"] != 0 or
            row["forced_b_result"] != 0 or
            row["forced_equal"] != 1
        ):
            raise MeasurementError(
                f"{label} candidate/forced encoder results disagree")
        payload = roles["repaired"]["payload_sha256"]
        if (
            not peel_codec._is_sha256(payload) or
            row["forced_a_payload_sha256"] != payload or
            row["forced_b_payload_sha256"] != payload
        ):
            raise MeasurementError(
                f"{label} forced-selected payload witness differs")
    else:
        if (
            row["forced_a_result"] != -1 or
            row["forced_b_result"] != -1 or
            row["forced_equal"] != 0 or
            row["forced_a_payload_sha256"] != "not_applicable" or
            row["forced_b_payload_sha256"] != "not_applicable"
        ):
            raise MeasurementError(
                f"{label} forced an uncommitted selector")

    expected_prefix = (
        max(
            roles["repaired"]["received_symbols"],
            roles["dispatch"]["received_symbols"],
        )
        if (
            roles["repaired"]["recovery_ok"] == 1 and
            roles["dispatch"]["recovery_ok"] == 1
        )
        else block_count + max_overhead
    )
    if (
        row["direct_fixed_prefix_symbols"] != expected_prefix or
        not block_count <= expected_prefix <=
            block_count + max_overhead or
        expected_prefix > block_count + max_overhead
    ):
        raise MeasurementError(
            f"{label} direct fixed-prefix witness is invalid")

    for prefix, expected_bytes in (
        (
            "repair",
            _candidate_intermediate_bytes(arm, block_count, block_bytes),
        ),
        (
            "dispatch",
            _dispatch_intermediate_bytes(block_count, block_bytes),
        ),
    ):
        intermediate_bytes = row[f"{prefix}_intermediate_bytes"]
        intermediate_sha = row[f"{prefix}_intermediate_sha256"]
        direct_executed = row[f"{prefix}_direct_executed"]
        if direct_executed not in (0, 1):
            raise MeasurementError(
                f"{label}.{prefix}_direct_executed is invalid")
        direct_result = _result_code(
            row[f"{prefix}_direct_result"],
            f"{label}.{prefix}_direct_result",
            allow_na=not direct_executed)
        direct_bytes = row[f"{prefix}_direct_intermediate_bytes"]
        direct_sha = row[f"{prefix}_direct_intermediate_sha256"]
        witness_equal = row[f"{prefix}_direct_witness_equal"]
        role_success = (
            roles["repaired" if prefix == "repair" else "dispatch"]
            ["encode_result"] == 0
        )
        if role_success:
            if (
                intermediate_bytes != expected_bytes or
                not peel_codec._is_sha256(intermediate_sha)
            ):
                raise MeasurementError(
                    f"{label}.{prefix} encoder intermediate is invalid")
        elif (
            intermediate_bytes != -1 or
            intermediate_sha != "not_applicable"
        ):
            raise MeasurementError(
                f"{label}.{prefix} fabricated encoder intermediate")
        if not direct_executed:
            if (
                direct_result != -1 or direct_bytes != -1 or
                direct_sha != "not_applicable" or witness_equal != 0
            ):
                raise MeasurementError(
                    f"{label}.{prefix} fabricated unexecuted direct work")
        elif direct_result == 0:
            if (
                not role_success or direct_bytes != expected_bytes or
                direct_sha != intermediate_sha or witness_equal != 1
            ):
                raise MeasurementError(
                    f"{label}.{prefix} direct witness drifted")
        elif (
            direct_bytes != -1 or
            direct_sha != "not_applicable" or
            witness_equal != 0
        ):
            raise MeasurementError(
                f"{label}.{prefix} failed direct witness is invalid")
    if (
        row["repair_direct_executed"] != row["committed"] or
        row["dispatch_direct_executed"] !=
            int(roles["dispatch"]["encode_result"] == 0)
    ):
        raise MeasurementError(
            f"{label} direct execution does not match available encoders")
    return roles


def _parse_attempt_row(row, *, expected_attempt, label):
    if row["attempt"] != expected_attempt:
        raise MeasurementError(f"{label} is missing or reordered")
    for kind in ("probe", "real"):
        executed = row[f"{kind}_executed"]
        result = row[f"{kind}_result"]
        available = row[f"{kind}_stats_available"]
        if executed not in (0, 1) or available not in (0, 1):
            raise MeasurementError(f"{label}.{kind} flags are invalid")
        if not executed:
            if result != -1 or available != 0:
                raise MeasurementError(
                    f"{label}.{kind} fabricated an unexecuted call")
        else:
            _result_code(result, f"{label}.{kind}_result")
        stats = _validate_stats(
            row, kind, available, f"{label}.{kind}",
            acceptance_max=1)
        if available and not executed:
            raise MeasurementError(
                f"{label}.{kind} exposed stats without a call")
        row[f"_{kind}_stats"] = stats
    return {
        "attempt": expected_attempt,
        "probe": {
            "executed": row["probe_executed"],
            "result": row["probe_result"],
            "stats_available": row["probe_stats_available"],
            "stats": dict(row["_probe_stats"]),
        },
        "real": {
            "executed": row["real_executed"],
            "result": row["real_result"],
            "stats_available": row["real_stats_available"],
            "stats": dict(row["_real_stats"]),
        },
    }


def _validate_policy(summary, attempts, roles, label):
    if len(attempts) != REPAIR_V1_ATTEMPT_CAP:
        raise MeasurementError(f"{label} does not contain exactly 8 attempts")
    initial = attempts[0]["real"]
    if not initial["executed"]:
        raise MeasurementError(f"{label} omitted attempt-zero real work")
    if roles["raw"]["encode_result"] != initial["result"]:
        raise MeasurementError(
            f"{label} raw receipt disagrees with attempt-zero real work")

    expected_selected = -1
    expected_result = initial["result"]
    expected_cap = 0
    expected_fatal = 0
    expected_probe_count = 0
    expected_real_indices = {0}
    considered = {0}

    if initial["result"] == 0:
        expected_selected = 0
    elif initial["result"] in (1, 8):
        terminal = False
        for attempt in range(REPAIR_V1_ATTEMPT_CAP):
            probe = attempts[attempt]["probe"]
            if not probe["executed"]:
                raise MeasurementError(
                    f"{label} stopped structural probes before a terminal "
                    "policy result")
            expected_probe_count += 1
            considered.add(attempt)
            result = probe["result"]
            if result == 1:
                continue
            terminal = True
            if result == 0:
                if attempt == 0:
                    expected_result = 8
                    expected_fatal = 1
                else:
                    expected_selected = attempt
                    expected_real_indices.add(attempt)
                    selected_real = attempts[attempt]["real"]
                    if not selected_real["executed"]:
                        raise MeasurementError(
                            f"{label} omitted selected real solve")
                    expected_result = selected_real["result"]
            else:
                expected_result = result
            break
        if not terminal:
            expected_result = 4
            expected_cap = 1
    # Every result other than NeedMore/Error is terminal at attempt-zero.

    executed_probes = {
        attempt["attempt"] for attempt in attempts
        if attempt["probe"]["executed"]
    }
    executed_reals = {
        attempt["attempt"] for attempt in attempts
        if attempt["real"]["executed"]
    }
    if (
        executed_probes != set(range(expected_probe_count)) or
        executed_reals != expected_real_indices
    ):
        raise MeasurementError(f"{label} violates lazy selector execution")
    expected_attempts = max(considered) + 1
    expected_real_count = len(expected_real_indices)
    expected_calls = expected_real_count + expected_probe_count
    if (
        summary["selector_result"] != expected_result or
        summary["selected_attempt"] != expected_selected or
        summary["attempts_executed"] != expected_attempts or
        summary["calls_executed"] != expected_calls or
        summary["real_configuration_calls"] != expected_real_count or
        summary["structural_probe_calls"] != expected_probe_count or
        summary["cap_exhausted"] != expected_cap or
        summary["fatal_attempt_zero_mismatch"] != expected_fatal or
        summary["oom"] != int(expected_result == 9) or
        summary["committed"] != int(expected_result == 0)
    ):
        raise MeasurementError(f"{label} contradicts repair-v1 policy")


def _expected_encoder_timing_stats(summary, attempts, label):
    """Derive deterministic selector/forced structural timing counters."""
    if not summary["committed"]:
        return {}
    executed = [
        call
        for attempt in attempts
        for call in (attempt["probe"], attempt["real"])
        if call["executed"]
    ]
    if (
        len(executed) != summary["calls_executed"] or
        any(not call["stats_available"] for call in executed)
    ):
        raise MeasurementError(
            f"{label} committed without complete per-call stats")
    totals = {}
    for name in REPAIRTIMING_STATS_FIELDS:
        total = sum(call["stats"][name] for call in executed)
        if total > 0xffffffffffffffff:
            raise MeasurementError(
                f"{label} selector stats aggregate overflows uint64")
        totals[name] = total
    selected = summary["selected_attempt"]
    selected_real = attempts[selected]["real"]
    if (
        not selected_real["executed"] or
        not selected_real["stats_available"] or
        selected_real["result"] != 0
    ):
        raise MeasurementError(
            f"{label} selected real stats are incomplete")
    return {
        "repair_selector_encoder": {
            name: totals[name]
            for name in REPAIRTIMING_STRUCTURAL_STATS_FIELDS
        },
        "repair_forced_encoder": {
            name: selected_real["stats"][name]
            for name in REPAIRTIMING_STRUCTURAL_STATS_FIELDS
        },
    }


def _structural_digest_text(arm, block_count, construction_root,
                            summary, attempts):
    lines = [
        "wirehair-wh2-repair-v1-selector-structural-v1",
        f"arm_id={arm.provisional_id:016x}",
        f"K={block_count}",
        f"construction_root={construction_root}",
        f"repair_policy_sha256={REPAIR_V1_POLICY_SHA256}",
    ]
    for name in (
            "selector_result",
            "selected_attempt",
            "attempts_executed",
            "calls_executed",
            "real_configuration_calls",
            "structural_probe_calls",
            "cap_exhausted",
            "fatal_attempt_zero_mismatch",
            "oom",
            "committed"):
        lines.append(f"{name}={summary[name]}")
    lines.append(f"descriptor_sha256={summary['descriptor_sha256']}")
    for attempt in attempts:
        probe = attempt["probe"]
        pieces = [
            f"attempt={attempt['attempt']}",
            f"probe_executed={probe['executed']}",
            f"probe_result={probe['result']}",
            f"probe_stats_available={probe['stats_available']}",
        ]
        pieces.extend(
            f"probe_{name}={probe['stats'][name]}"
            for name in REPAIRTIMING_STRUCTURAL_STATS_FIELDS
        )
        lines.append(",".join(pieces))
    return "\n".join(lines) + "\n"


def _make_selector_witness(arm, block_count, construction_root,
                           summary, attempts):
    structural_attempts = []
    for attempt in attempts:
        probe = attempt["probe"]
        structural_attempts.append({
            "attempt": attempt["attempt"],
            "probe_executed": probe["executed"],
            "probe_result": probe["result"],
            "probe_stats_available": probe["stats_available"],
            "probe_stats": {
                name: probe["stats"][name]
                for name in REPAIRTIMING_STRUCTURAL_STATS_FIELDS
            },
        })
    return {
        "schema": "wirehair.wh2.repair-v1.selector-structural.v1",
        "arm_id": arm.provisional_id,
        "K": block_count,
        "construction_root": construction_root,
        "repair_policy_sha256": REPAIR_V1_POLICY_SHA256,
        "selector_result": summary["selector_result"],
        "selected_attempt": summary["selected_attempt"],
        "attempts_executed": summary["attempts_executed"],
        "calls_executed": summary["calls_executed"],
        "real_configuration_calls":
            summary["real_configuration_calls"],
        "structural_probe_calls": summary["structural_probe_calls"],
        "cap_exhausted": summary["cap_exhausted"],
        "fatal_attempt_zero_mismatch":
            summary["fatal_attempt_zero_mismatch"],
        "oom": summary["oom"],
        "committed": summary["committed"],
        "descriptor_sha256": summary["descriptor_sha256"],
        "attempts": structural_attempts,
        "selector_structural_sha256":
            summary["selector_structural_sha256"],
    }


def _timing_scope_for_role(role):
    if role == "repair_forced_encoder":
        return "encoder_forced"
    if role == "repair_selector_encoder":
        return "encoder_selector"
    if role == "wh1_encoder":
        return "encoder_wh1"
    if role in ("repair_decoder_feed", "wh1_decoder_feed"):
        return "decoder_feed"
    if role in ("repair_decoder_full", "wh1_decoder_full"):
        return "decoder_full"
    if role in ("repair_direct", "dispatch_direct"):
        return "direct"
    raise MeasurementError(f"unknown repairtiming role {role!r}")


def _timing_role_summary(role, roles, summary):
    if role.startswith("repair_"):
        recovered = roles["repaired"]
    elif role.startswith("wh1_"):
        recovered = roles["wh1"]
    elif role == "dispatch_direct":
        recovered = roles["dispatch"]
    else:
        raise MeasurementError(f"unknown repairtiming role {role!r}")

    scope = _timing_scope_for_role(role)
    if scope in ("encoder_selector", "encoder_forced", "encoder_wh1"):
        return {
            "construct_result": recovered["encode_result"],
            "result": recovered["encode_result"],
            "recover_result": -1,
            "recovery_ok": 0,
            "success": recovered["encode_result"] == 0,
        }
    if scope == "decoder_feed":
        return {
            "construct_result": recovered["decode_construct_result"],
            "result": recovered["feed_result"],
            "recover_result": -1,
            "recovery_ok": 0,
            "success": (
                recovered["decode_construct_result"] == 0 and
                recovered["feed_result"] == 0
            ),
        }
    if scope == "decoder_full":
        return {
            "construct_result": recovered["decode_construct_result"],
            "result": recovered["feed_result"],
            "recover_result": recovered["recover_result"],
            "recovery_ok": recovered["recovery_ok"],
            "success": (
                recovered["decode_construct_result"] == 0 and
                recovered["feed_result"] == 0 and
                recovered["recover_result"] == 0 and
                recovered["recovery_ok"] == 1
            ),
        }
    prefix = "repair" if role == "repair_direct" else "dispatch"
    return {
        "construct_result": recovered["encode_result"],
        "result": summary[f"{prefix}_direct_result"],
        "recover_result": -1,
        "recovery_ok": 0,
        "success": (
            recovered["encode_result"] == 0 and
            summary[f"{prefix}_direct_result"] == 0
        ),
    }


def _timing_intermediate(role, arm, block_count, block_bytes, summary):
    if (
        role.startswith("wh1_") or
        role in ("repair_decoder_feed", "repair_decoder_full")
    ):
        return -1, "not_applicable"
    if role == "dispatch_direct":
        return (
            _dispatch_intermediate_bytes(block_count, block_bytes),
            summary["dispatch_direct_intermediate_sha256"],
        )
    return (
        _candidate_intermediate_bytes(arm, block_count, block_bytes),
        (
            summary["repair_direct_intermediate_sha256"]
            if role == "repair_direct"
            else summary["repair_intermediate_sha256"]
        ),
    )


def _parse_timing_row(
        row, *, arm, block_count, block_bytes, inner_reps,
        expected_cpu, panel_index, slot, summary, roles, label):
    panel, unused_scope, left, right = \
        REPAIRTIMING_PANEL_SPECS[panel_index]
    expected_label = REPAIRTIMING_ORDER[slot]
    expected_role = left if expected_label == "A" else right
    expected_scope = _timing_scope_for_role(expected_role)
    if (
        row["timing_panel"] != panel or
        row["timing_panel_index"] != panel_index or
        row["timing_slot"] != slot or
        row["timing_pair"] != slot // 2 or
        row["timing_label"] != expected_label or
        row["timing_role"] != expected_role or
        row["timing_scope"] != expected_scope
    ):
        raise MeasurementError(f"{label} is reordered or mislabeled")

    expected = _timing_role_summary(expected_role, roles, summary)
    panel_roles = (
        _timing_role_summary(left, roles, summary),
        _timing_role_summary(right, roles, summary),
    )
    requires_repair = (
        left.startswith("repair_") or right.startswith("repair_"))
    selector_ready = (
        not requires_repair or summary["committed"] == 1)
    panel_success = all(item["success"] for item in panel_roles)
    should_execute = selector_ready and panel_success
    expected_prefix = (
        summary["direct_fixed_prefix_symbols"]
        if should_execute and expected_scope == "direct" else -1
    )
    if row["timing_fixed_prefix_symbols"] != expected_prefix:
        raise MeasurementError(f"{label} changed the direct fixed prefix")
    if not should_execute:
        expected_reason = (
            "selector_uncommitted" if not selector_ready
            else "role_preflight_failed")
        if (
            row["timing_executed"] != 0 or
            row["timing_eligible"] != 0 or
            row["timing_censor_reason"] != expected_reason
        ):
            raise MeasurementError(
                f"{label} has dishonest preflight censoring")
        for name in (
                "timing_construct_result", "timing_result",
                "timing_recover_result", "timing_recovery_ok",
                "timing_stable", "timing_elapsed_ns",
                "timing_inner_reps", "timing_saturated",
                "timing_cpu_before", "timing_cpu_after",
                "timing_cpu_migrated", "timing_minflt",
                "timing_majflt", "timing_fault_contaminated",
                "timing_stats_available"):
            if row[name] != -1:
                raise MeasurementError(f"{label} did work while censored")
        _validate_stats(row, "timing", 0, label)
        if (
            row["timing_intermediate_bytes"] != -1 or
            row["timing_intermediate_sha256"] != "not_applicable"
        ):
            raise MeasurementError(
                f"{label} fabricated a censored intermediate witness")
        return {
            name[len("timing_"):]: row[name]
            for name in REPAIRTIMING_TIMING_FIELDS
        }

    if (
        row["timing_executed"] != 1 or
        row["timing_eligible"] not in (0, 1) or
        row["timing_censor_reason"] not in (
            "none", "cpu_migrated", "fault_contaminated",
            "saturated", "outcome_unstable")
    ):
        raise MeasurementError(f"{label} has invalid execution flags")
    for name in (
            "timing_construct_result", "timing_result",
            "timing_recover_result"):
        _result_code(
            row[name], f"{label}.{name}", allow_na=(name.endswith("recover_result")))
    if (
        row["timing_recovery_ok"] not in (0, 1) or
        row["timing_stable"] not in (0, 1) or
        row["timing_elapsed_ns"] < 1 or
        row["timing_inner_reps"] != inner_reps or
        row["timing_saturated"] not in (0, 1) or
        row["timing_cpu_before"] < 0 or
        row["timing_cpu_after"] < 0 or
        row["timing_cpu_migrated"] != int(
            row["timing_cpu_before"] != row["timing_cpu_after"]) or
        row["timing_minflt"] < 0 or row["timing_majflt"] < 0 or
        row["timing_fault_contaminated"] != int(
            row["timing_minflt"] > 0 or row["timing_majflt"] > 0)
    ):
        raise MeasurementError(f"{label} timing receipt is impossible")

    construct_result = row["timing_construct_result"]
    result = row["timing_result"]
    recover_result = row["timing_recover_result"]
    recovery_ok = row["timing_recovery_ok"]
    if construct_result != 0 and result != construct_result:
        raise MeasurementError(
            f"{label} continued after failed construction")
    if expected_scope == "direct" and construct_result != 0:
        raise MeasurementError(
            f"{label} has a failed direct precondition")
    if expected_scope == "decoder_full":
        if result == 0:
            if recover_result == -1 or (
                    (recover_result == 0) != (recovery_ok == 1)):
                raise MeasurementError(
                    f"{label} has an impossible recovery stage")
        elif recover_result != -1 or recovery_ok != 0:
            raise MeasurementError(
                f"{label} recovered after a failed feed")
    elif recover_result != -1 or recovery_ok != 0:
        raise MeasurementError(
            f"{label} populated an inapplicable recovery stage")

    expected_stats_available = int(
        expected_scope == "direct" or
        (
            expected_role.startswith("repair_") and
            result == 0
        )
    )
    if row["timing_stats_available"] != expected_stats_available:
        raise MeasurementError(f"{label} has dishonest stats availability")
    call_limit = (
        REPAIRTIMING_SELECTOR_CALL_CAP
        if expected_role == "repair_selector_encoder"
        else 1
    )
    _validate_stats(
        row, "timing", row["timing_stats_available"], label,
        elapsed_ns=row["timing_elapsed_ns"],
        acceptance_max=call_limit,
        call_limit=call_limit)
    success_bytes, success_sha = _timing_intermediate(
        expected_role, arm, block_count, block_bytes, summary)
    witness_role = expected_role in (
        "repair_selector_encoder",
        "repair_forced_encoder",
        "repair_direct",
        "dispatch_direct",
    )
    expected_bytes, expected_sha = (
        (success_bytes, success_sha)
        if witness_role and row["timing_result"] == 0
        else (-1, "not_applicable")
    )
    observed_intermediate = (
        row["timing_intermediate_bytes"],
        row["timing_intermediate_sha256"],
    )
    if observed_intermediate != (expected_bytes, expected_sha):
        raise MeasurementError(f"{label} intermediate witness drifted")
    observable_matches = (
        row["timing_construct_result"] == expected["construct_result"] and
        row["timing_result"] == expected["result"] and
        row["timing_recover_result"] == expected["recover_result"] and
        row["timing_recovery_ok"] == expected["recovery_ok"] and
        observed_intermediate == (expected_bytes, expected_sha)
    )
    if row["timing_stable"] == 1 and not observable_matches:
        raise MeasurementError(
            f"{label} claims stability after an outcome change")
    outcome_unstable = (
        row["timing_stable"] != 1 or not observable_matches
    )
    local_reason = (
        "cpu_migrated" if row["timing_cpu_migrated"] else
        "fault_contaminated" if row["timing_fault_contaminated"] else
        "saturated" if row["timing_saturated"] else
        "outcome_unstable" if outcome_unstable else
        "none"
    )
    normalized = {
        name[len("timing_"):]: row[name]
        for name in REPAIRTIMING_TIMING_FIELDS
    }
    normalized["_local_censor_reason"] = local_reason
    normalized["_expected_cpu"] = expected_cpu
    return normalized


def _validate_stable_timing_signature(
        row, signatures, expected_encoder_stats, label):
    """Cross-check deterministic fields for every stable role execution."""
    if row["timing_executed"] != 1 or row["timing_stable"] != 1:
        return
    role = row["timing_role"]
    expected_stats = expected_encoder_stats.get(role)
    if expected_stats is not None and any(
            row[f"timing_{name}"] != expected_stats[name]
            for name in REPAIRTIMING_STRUCTURAL_STATS_FIELDS
    ):
        raise MeasurementError(
            f"{label} changed {role} aggregate structural stats")
    signature = (
        row["timing_construct_result"],
        row["timing_result"],
        row["timing_recover_result"],
        row["timing_recovery_ok"],
        row["timing_stats_available"],
        *(row[f"timing_{name}"]
          for name in REPAIRTIMING_STRUCTURAL_STATS_FIELDS),
        row["timing_fixed_prefix_symbols"],
        row["timing_intermediate_bytes"],
        row["timing_intermediate_sha256"],
    )
    prior = signatures.setdefault(role, signature)
    if signature != prior:
        raise MeasurementError(
            f"{label} changed stable structural work for role {role}")


def _validate_timing_panel(rows, label):
    if len(rows) != len(REPAIRTIMING_ORDER):
        raise MeasurementError(f"{label} is incomplete")
    if not rows[0]["executed"]:
        expected_reason = rows[0]["censor_reason"]
        if any(
                row["executed"] != 0 or row["eligible"] != 0 or
                row["censor_reason"] != expected_reason
                for row in rows):
            raise MeasurementError(f"{label} preflight censoring changed")
        return
    if any(row["executed"] != 1 for row in rows):
        raise MeasurementError(f"{label} partially executed")
    local_reasons = {row["_local_censor_reason"] for row in rows}
    expected_reason = next(
        (
            reason for reason in (
                "cpu_migrated", "fault_contaminated", "saturated",
                "outcome_unstable")
            if reason in local_reasons
        ),
        "none",
    )
    expected_eligible = int(expected_reason == "none")
    if any(
            row["eligible"] != expected_eligible or
            row["censor_reason"] != expected_reason
            for row in rows):
        raise MeasurementError(
            f"{label} has the wrong panel-wide censor reason")
    if expected_eligible and any(
            row["cpu_before"] != row["_expected_cpu"] or
            row["cpu_after"] != row["_expected_cpu"] or
            not row["stats_available"] and
                not row["role"].startswith("wh1_")
            for row in rows):
        raise MeasurementError(
            f"{label} marked incomplete measurements timing-eligible")


_REPAIRTIMING_DONE_FIELDS = (
    "complete",
    "cells",
    "summary_rows",
    "attempt_rows",
    "timing_rows",
    "rows",
    "finished_monotonic_ns",
    "stream_sha256",
)


def _attempt_values(attempt):
    values = [
        attempt["attempt"],
        attempt["probe"]["executed"],
        attempt["probe"]["result"],
        attempt["probe"]["stats_available"],
    ]
    values.extend(
        attempt["probe"]["stats"][name]
        for name in REPAIRTIMING_STATS_FIELDS)
    values.extend((
        attempt["real"]["executed"],
        attempt["real"]["result"],
        attempt["real"]["stats_available"],
    ))
    values.extend(
        attempt["real"]["stats"][name]
        for name in REPAIRTIMING_STATS_FIELDS)
    if len(values) != len(REPAIRTIMING_ATTEMPT_FIELDS):
        raise AssertionError("normalized attempt schema drifted")
    return values


def _normalized_timing_values(row):
    values = [row[name] for name in REPAIRTIMING_TIMING_FIELDS]
    if len(values) != len(REPAIRTIMING_TIMING_FIELDS):
        raise AssertionError("normalized timing schema drifted")
    return values


def parse_repairtiming_output(
        stdout, *, block_count, block_bytes, repair_arm,
        repair_policy=REPAIR_V1_POLICY_NAME,
        dispatch_profile=peel_codec.TARGET_PROFILE,
        construction_seed, construction_seed_derivation, loss, loss_seed,
        schedule, warmup_replicates, replicates, inner_reps, max_overhead,
        cache_state, systematic_cache, evict_bytes, context,
        required_margin):
    """Validate one native v3 stream without aggregating campaign outcomes."""
    arm = globals()["repair_arm"](repair_arm)
    validate_repairtiming_dimensions(
        block_count=block_count,
        block_bytes=block_bytes,
        repair_arm=arm,
        repair_policy=repair_policy,
        dispatch_profile=dispatch_profile,
        construction_seed=construction_seed,
        construction_seed_derivation=construction_seed_derivation,
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
    if not isinstance(stdout, str) or not stdout:
        raise MeasurementError("repairtiming returned no output")
    if len(stdout) > REPAIRTIMING_NATIVE_STDOUT_BYTE_LIMIT:
        raise MeasurementError("repairtiming output exceeds its byte limit")
    if not stdout.isascii():
        raise MeasurementError("repairtiming output is not ASCII")
    if not isinstance(context, dict):
        raise MeasurementError("repairtiming context is not an object")
    try:
        context_canonical = json.dumps(
            context, sort_keys=True, separators=(",", ":"),
            ensure_ascii=True, allow_nan=False).encode("ascii")
    except (TypeError, ValueError):
        raise MeasurementError("repairtiming context is not canonical JSON")
    if len(context_canonical) > REPAIRTIMING_CONTEXT_CANONICAL_BYTE_LIMIT:
        raise MeasurementError(
            "repairtiming context exceeds its canonical byte limit")
    context_sha256 = peel_codec._validate_paired_bound_context(
        context, cache_state=cache_state, evict_bytes=evict_bytes,
        require_capture=True)
    affinity = context["bound"]["cpu_affinity"]
    if len(affinity) != 1:
        raise MeasurementError(
            "repairtiming context must bind one assigned CPU")
    expected_cpu = affinity[0]

    cells_count = warmup_replicates + replicates
    expected_rows = cells_count * REPAIRTIMING_ROWS_PER_CELL
    stream = io.StringIO(stdout)
    stream_digest = hashlib.sha256()

    manifest_line = _read_line(
        stream, "manifest", byte_limit=REPAIRTIMING_MANIFEST_BYTE_LIMIT)
    stream_digest.update(manifest_line.encode("ascii"))
    manifest_raw = peel_codec._parse_ordered_kv_line(
        manifest_line[:-1],
        "# repairtiming,",
        REPAIRTIMING_MANIFEST_FIELDS,
        "repairtiming manifest",
    )
    manifest = _parse_manifest(
        manifest_raw,
        arm=arm,
        block_count=block_count,
        block_bytes=block_bytes,
        dispatch_profile=dispatch_profile,
        construction_seed=construction_seed,
        construction_seed_derivation=construction_seed_derivation,
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
        context_sha256=context_sha256,
        required_margin=float(required_margin),
    )

    header_specs = (
        (
            "# repairtiming_summary_columns,",
            REPAIRTIMING_SUMMARY_COLUMNS,
        ),
        (
            "# repairtiming_attempt_columns,",
            REPAIRTIMING_ATTEMPT_COLUMNS,
        ),
        (
            "# repairtiming_timing_columns,",
            REPAIRTIMING_TIMING_COLUMNS,
        ),
    )
    for prefix, columns in header_specs:
        line = _read_line(
            stream, "typed column header",
            byte_limit=REPAIRTIMING_HEADER_BYTE_LIMIT)
        stream_digest.update(line.encode("ascii"))
        if line != prefix + ",".join(columns) + "\n":
            raise MeasurementError(
                "repairtiming typed columns are missing or reordered")

    cells = []
    seen_selectors = {}
    seen_loss_seeds = set()
    seen_traces = set()
    seen_roots = set()
    total_elapsed_ns = 0
    row_index = 0
    for cell_index in range(cells_count):
        measured = int(cell_index >= warmup_replicates)
        construction_root = _expected_construction_root(
            construction_seed, cell_index,
            construction_seed_derivation)
        cell_loss_seed = _expected_loss_seed(loss_seed, cell_index)

        line = _read_line(
            stream, f"summary row for cell {cell_index}",
            byte_limit=REPAIRTIMING_SUMMARY_ROW_BYTE_LIMIT)
        stream_digest.update(line.encode("ascii"))
        summary = _parse_typed_row(
            line[:-1], row_index, REPAIRTIMING_SUMMARY_COLUMNS)
        row_index += 1
        common_hashes = _validate_common(
            summary,
            row_kind="summary",
            cell_index=cell_index,
            measured=measured,
            construction_root=construction_root,
            loss_seed=cell_loss_seed,
            common_hashes=None,
            label=f"repairtiming cell {cell_index} summary",
        )
        roles = _validate_summary_row(
            summary,
            arm=arm,
            block_count=block_count,
            block_bytes=block_bytes,
            max_overhead=max_overhead,
            label=f"repairtiming cell {cell_index} summary",
        )

        attempts = []
        for attempt_index in range(REPAIR_V1_ATTEMPT_CAP):
            line = _read_line(
                stream,
                f"attempt {attempt_index} row for cell {cell_index}",
                byte_limit=REPAIRTIMING_ATTEMPT_ROW_BYTE_LIMIT)
            stream_digest.update(line.encode("ascii"))
            attempt_row = _parse_typed_row(
                line[:-1], row_index, REPAIRTIMING_ATTEMPT_COLUMNS)
            row_index += 1
            _validate_common(
                attempt_row,
                row_kind="attempt",
                cell_index=cell_index,
                measured=measured,
                construction_root=construction_root,
                loss_seed=cell_loss_seed,
                common_hashes=common_hashes,
                label=(
                    f"repairtiming cell {cell_index} "
                    f"attempt {attempt_index}"
                ),
            )
            attempts.append(_parse_attempt_row(
                attempt_row,
                expected_attempt=attempt_index,
                label=(
                    f"repairtiming cell {cell_index} "
                    f"attempt {attempt_index}"
                ),
            ))
        _validate_policy(
            summary, attempts, roles,
            f"repairtiming cell {cell_index}")
        expected_encoder_stats = _expected_encoder_timing_stats(
            summary, attempts, f"repairtiming cell {cell_index}")
        structural_text = _structural_digest_text(
            arm, block_count, construction_root, summary, attempts)
        structural_sha256 = hashlib.sha256(
            structural_text.encode("ascii")).hexdigest()
        if summary["selector_structural_sha256"] != structural_sha256:
            raise MeasurementError(
                f"repairtiming cell {cell_index} structural digest differs")

        normalized_timing_rows = []
        stable_timing_signatures = {}
        for panel_index, panel_spec in enumerate(
                REPAIRTIMING_PANEL_SPECS):
            normalized_panel = []
            for slot in range(len(REPAIRTIMING_ORDER)):
                line = _read_line(
                    stream,
                    (
                        f"timing panel {panel_index} slot {slot} "
                        f"for cell {cell_index}"
                    ),
                    byte_limit=REPAIRTIMING_TIMING_ROW_BYTE_LIMIT)
                stream_digest.update(line.encode("ascii"))
                timing_row = _parse_typed_row(
                    line[:-1], row_index, REPAIRTIMING_TIMING_COLUMNS)
                row_index += 1
                _validate_common(
                    timing_row,
                    row_kind="timing",
                    cell_index=cell_index,
                    measured=measured,
                    construction_root=construction_root,
                    loss_seed=cell_loss_seed,
                    common_hashes=common_hashes,
                    label=(
                        f"repairtiming cell {cell_index} "
                        f"timing panel {panel_index} slot {slot}"
                    ),
                )
                normalized = _parse_timing_row(
                    timing_row,
                    arm=arm,
                    block_count=block_count,
                    block_bytes=block_bytes,
                    inner_reps=inner_reps,
                    expected_cpu=expected_cpu,
                    panel_index=panel_index,
                    slot=slot,
                    summary=summary,
                    roles=roles,
                    label=(
                        f"repairtiming cell {cell_index} "
                        f"timing panel {panel_index} slot {slot}"
                    ),
                )
                _validate_stable_timing_signature(
                    timing_row,
                    stable_timing_signatures,
                    expected_encoder_stats,
                    (
                        f"repairtiming cell {cell_index} "
                        f"timing panel {panel_index} slot {slot}"
                    ),
                )
                normalized_panel.append(normalized)
                normalized_timing_rows.append(
                    _normalized_timing_values(timing_row))
                if timing_row["timing_executed"] == 1:
                    total_elapsed_ns += timing_row["timing_elapsed_ns"]
            _validate_timing_panel(
                normalized_panel,
                (
                    f"repairtiming cell {cell_index} "
                    f"panel {panel_spec[0]}"
                ),
            )

        selector_witness = _make_selector_witness(
            arm, block_count, construction_root, summary, attempts)
        selector_key = (
            f"{arm.provisional_id:016x}",
            block_count,
            construction_root,
        )
        prior = seen_selectors.setdefault(
            selector_key, selector_witness)
        if not peel_codec._same_typed_json(prior, selector_witness):
            raise MeasurementError(
                "repairtiming repeated selector key changed structure")

        controls = {
            name: summary[name] for name in REPAIRTIMING_CONTROL_FIELDS
        }
        real_witness = {
            "schema": "wirehair.wh2.repairtiming.real-witness.v3",
            "trace_sha256": summary["trace_sha256"],
            "message_sha256": summary["message_sha256"],
            "descriptor_hex": summary["descriptor_hex"],
            "roles": _json_safe_copy(roles),
            "controls": controls,
            "attempt_columns": list(REPAIRTIMING_ATTEMPT_FIELDS),
            "attempt_rows": [
                _attempt_values(attempt) for attempt in attempts
            ],
            "timing_columns": list(REPAIRTIMING_TIMING_FIELDS),
            "timing_rows": normalized_timing_rows,
        }
        cells.append(RepairTimingCell(
            replicate=cell_index,
            construction_root=construction_root,
            loss_seed=cell_loss_seed,
            selector_key=selector_key,
            selector_witness=selector_witness,
            real_witness=real_witness,
        ))
        if cell_loss_seed in seen_loss_seeds:
            raise MeasurementError("repairtiming reused a loss seed")
        if summary["trace_sha256"] in seen_traces:
            raise MeasurementError("repairtiming reused a trace")
        seen_loss_seeds.add(cell_loss_seed)
        seen_traces.add(summary["trace_sha256"])
        seen_roots.add(construction_root)

    if row_index != expected_rows:
        raise MeasurementError("repairtiming parsed the wrong row count")

    done_line = _read_line(
        stream, "completion", byte_limit=4096)
    done = peel_codec._parse_ordered_kv_line(
        done_line[:-1],
        "# repairtiming_done,",
        _REPAIRTIMING_DONE_FIELDS,
        "repairtiming completion",
    )
    done_integer_fields = (
        "cells", "summary_rows", "attempt_rows", "timing_rows", "rows",
        "finished_monotonic_ns",
    )
    done_values = {
        name: _parse_decimal(
            done[name], f"repairtiming completion.{name}", minimum=0)
        for name in done_integer_fields
    }
    expected_done = {
        "cells": cells_count,
        "summary_rows": cells_count,
        "attempt_rows": cells_count * REPAIR_V1_ATTEMPT_CAP,
        "timing_rows": (
            cells_count * len(REPAIRTIMING_PANELS) *
            len(REPAIRTIMING_ORDER)
        ),
        "rows": expected_rows,
    }
    if (
        done["complete"] != "1" or
        any(done_values[name] != value
            for name, value in expected_done.items()) or
        done_values["finished_monotonic_ns"] < 1 or
        not peel_codec._is_sha256(done["stream_sha256"])
    ):
        raise MeasurementError("repairtiming completion is invalid")
    done_prefix = done_line[:-1].rsplit("stream_sha256=", 1)[0]
    done_prefix += "stream_sha256="
    stream_digest.update(done_prefix.encode("ascii"))
    stream_sha256 = stream_digest.hexdigest()
    if done["stream_sha256"] != stream_sha256:
        raise MeasurementError("repairtiming stream SHA-256 does not match")
    if stream.read(1):
        raise MeasurementError("repairtiming has trailing output")

    started_ns = manifest["started_monotonic_ns"]
    finished_ns = done_values["finished_monotonic_ns"]
    if finished_ns <= started_ns:
        raise MeasurementError("repairtiming completion clock did not advance")
    final_context_sha256 = peel_codec._validate_paired_context(
        context,
        cache_state=cache_state,
        evict_bytes=evict_bytes,
        started_ns=started_ns,
        finished_ns=finished_ns,
    )
    if final_context_sha256 != context_sha256:
        raise MeasurementError("repairtiming context changed during parsing")
    if total_elapsed_ns > finished_ns - started_ns:
        raise MeasurementError(
            "repairtiming elapsed work exceeds its sequential run interval")
    expected_unique_roots = (
        1 if construction_seed_derivation ==
            REPAIRTIMING_CONSTRUCTION_SEED_FIXED
        else cells_count
    )
    if len(seen_roots) != expected_unique_roots:
        raise MeasurementError(
            "repairtiming construction roots are missing or aliased")

    return RepairTimingEvidence(
        protocol=REPAIRTIMING_PROTOCOL_V3,
        manifest=manifest,
        context=_json_safe_copy(context),
        cells=tuple(cells),
        stream_sha256=stream_sha256,
        evidence={
            "schema": "wirehair.wh2.repairtiming.evidence.v3",
            "started_monotonic_ns": started_ns,
            "finished_monotonic_ns": finished_ns,
            "summary_rows": cells_count,
            "attempt_rows": cells_count * REPAIR_V1_ATTEMPT_CAP,
            "timing_rows": (
                cells_count * len(REPAIRTIMING_PANELS) *
                len(REPAIRTIMING_ORDER)
            ),
            "rows": expected_rows,
            "stdout_bytes": len(stdout),
            "context_sha256": final_context_sha256,
        },
        native_stdout=stdout,
    )


def _repairtiming_run_bounds(stdout):
    """Extract authenticated run clocks without splitting the large stream."""
    if (
        not isinstance(stdout, str) or not stdout or
        len(stdout) > REPAIRTIMING_NATIVE_STDOUT_BYTE_LIMIT or
        not stdout.isascii()
    ):
        raise MeasurementError("repairtiming output is missing or oversized")
    first_end = stdout.find("\n")
    last_start = stdout.rfind("\n", 0, len(stdout) - 1)
    if (
        first_end < 0 or last_start < 0 or
        not stdout.endswith("\n") or last_start <= first_end
    ):
        raise MeasurementError("repairtiming output is incomplete")
    manifest = peel_codec._parse_ordered_kv_line(
        stdout[:first_end],
        "# repairtiming,",
        REPAIRTIMING_MANIFEST_FIELDS,
        "repairtiming manifest",
    )
    done = peel_codec._parse_ordered_kv_line(
        stdout[last_start + 1:-1],
        "# repairtiming_done,",
        _REPAIRTIMING_DONE_FIELDS,
        "repairtiming completion",
    )
    started_ns = _parse_decimal(
        manifest["started_monotonic_ns"],
        "repairtiming manifest.started_monotonic_ns",
        minimum=1)
    finished_ns = _parse_decimal(
        done["finished_monotonic_ns"],
        "repairtiming completion.finished_monotonic_ns",
        minimum=1)
    if finished_ns <= started_ns:
        raise MeasurementError("repairtiming completion clock did not advance")
    return started_ns, finished_ns


def _arm_repairtiming_parent_death(expected_parent_pid):
    """Arm Linux parent-death SIGKILL, closing the fork-to-prctl race."""
    if (
        _REPAIRTIMING_PRCTL is None or
        _REPAIRTIMING_PRCTL(
            _REPAIRTIMING_PR_SET_PDEATHSIG,
            signal.SIGKILL,
            0,
            0,
            0,
        ) != 0
    ):
        os._exit(127)
    # PR_SET_PDEATHSIG is not retroactive.  If the campaign worker died
    # between fork() and prctl(), terminate before exec rather than orphaning.
    if os.getppid() != expected_parent_pid:
        os.kill(os.getpid(), signal.SIGKILL)
        os._exit(127)


def _repairtiming_process_group_exists(process_group):
    try:
        os.killpg(process_group, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _signal_repairtiming_process_group(process_group, signal_number):
    try:
        os.killpg(process_group, signal_number)
    except ProcessLookupError:
        return False
    return True


def _terminate_repairtiming_process_group(process):
    """TERM, KILL, reap the leader, and confirm its private group is gone."""
    process_group = process.pid
    if _repairtiming_process_group_exists(process_group):
        _signal_repairtiming_process_group(
            process_group, signal.SIGTERM)
        deadline = time.monotonic() + REPAIRTIMING_TERMINATE_GRACE_SECONDS
        while (
            _repairtiming_process_group_exists(process_group) and
            time.monotonic() < deadline
        ):
            process.poll()
            time.sleep(0.005)
        if _repairtiming_process_group_exists(process_group):
            _signal_repairtiming_process_group(
                process_group, signal.SIGKILL)
    elif process.poll() is None:
        # This is only reachable if child setup failed before setsid().
        process.kill()

    try:
        process.wait(timeout=REPAIRTIMING_REAP_TIMEOUT_SECONDS)
    except subprocess.TimeoutExpired:
        if _repairtiming_process_group_exists(process_group):
            _signal_repairtiming_process_group(
                process_group, signal.SIGKILL)
        else:
            process.kill()
        try:
            process.wait(timeout=REPAIRTIMING_REAP_TIMEOUT_SECONDS)
        except subprocess.TimeoutExpired as error:
            raise MeasurementError(
                "could not promptly reap repairtiming benchmark leader"
            ) from error

    deadline = time.monotonic() + REPAIRTIMING_GROUP_CONFIRM_SECONDS
    while (
        _repairtiming_process_group_exists(process_group) and
        time.monotonic() < deadline
    ):
        time.sleep(0.005)
    if _repairtiming_process_group_exists(process_group):
        raise MeasurementError(
            "repairtiming benchmark process group survived SIGKILL")


def _repairtiming_diagnostic(stdout, stderr, limit):
    """Return one bounded, escaped diagnostic from already bounded pipes."""
    label, payload = (
        ("stderr", stderr) if stderr else
        ("stdout", stdout) if stdout else
        ("output", "")
    )
    # Slice before repr(): escaping an otherwise valid 150-MiB stdout must not
    # transiently allocate a many-times-larger failure message.
    sample = payload[:limit]
    rendered = repr(sample)
    truncated = len(sample) != len(payload) or len(rendered) > limit
    rendered = rendered[:limit]
    if truncated:
        rendered += "... [truncated]"
    return f"{label}={rendered}"


def _repairtiming_command_label(command):
    """Render only a bounded executable label, never attacker-sized argv."""
    try:
        executable = os.fsdecode(os.fspath(command[0]))
    except (IndexError, TypeError, ValueError):
        return "<invalid-command>"
    sample = executable[:256]
    rendered = repr(sample)
    truncated = len(sample) != len(executable) or len(rendered) > 256
    rendered = rendered[:256]
    if truncated:
        rendered += "... [truncated]"
    return rendered


def _run_repairtiming_checked(
        command, env, *, stdout_limit=None, stderr_limit=None,
        diagnostic_limit=None):
    """Run repairtiming while bounding and concurrently draining both pipes."""
    if stdout_limit is None:
        stdout_limit = REPAIRTIMING_NATIVE_STDOUT_BYTE_LIMIT
    if stderr_limit is None:
        stderr_limit = REPAIRTIMING_STDERR_BYTE_LIMIT
    if diagnostic_limit is None:
        diagnostic_limit = REPAIRTIMING_DIAGNOSTIC_BYTE_LIMIT
    for label, value in (
        ("stdout", stdout_limit),
        ("stderr", stderr_limit),
        ("diagnostic", diagnostic_limit),
    ):
        if type(value) is not int or value < 1:
            raise MeasurementError(
                f"repairtiming {label} byte limit must be a positive integer")
    if _REPAIRTIMING_PRCTL is None:
        raise MeasurementError(
            "repairtiming subprocess isolation requires Linux prctl")

    command_label = _repairtiming_command_label(command)
    expected_parent_pid = os.getpid()
    try:
        # A private session lets this runner clean up every benchmark
        # descendant.  PDEATHSIG covers abrupt campaign-worker death despite
        # that isolation; the child-side PPID check closes the preexec race.
        process = subprocess.Popen(
            command,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=env,
            bufsize=0,
            start_new_session=True,
            preexec_fn=functools.partial(
                _arm_repairtiming_parent_death, expected_parent_pid),
        )
    except (OSError, subprocess.SubprocessError, ValueError) as error:
        raw_detail = str(error)
        sample = raw_detail[:diagnostic_limit]
        rendered = repr(sample)
        detail = rendered[:diagnostic_limit]
        if (
            len(sample) != len(raw_detail) or
            len(rendered) > diagnostic_limit
        ):
            detail += "... [truncated]"
        raise MeasurementError(
            f"could not execute repairtiming benchmark {command_label}: "
            f"{detail}") from error

    selector = None
    streams = {}
    buffers = {}
    limits = {}
    returncode = None
    overflow = None
    pipe_holder = False
    unexpected_descendants = False
    post_exit_deadline = None
    try:
        selector = selectors.DefaultSelector()
        streams = {
            "stdout": process.stdout,
            "stderr": process.stderr,
        }
        buffers = {
            "stdout": bytearray(),
            "stderr": bytearray(),
        }
        limits = {
            "stdout": stdout_limit,
            "stderr": stderr_limit,
        }
        for name, stream in streams.items():
            if stream is None:
                raise MeasurementError(
                    f"repairtiming benchmark has no {name} pipe")
            os.set_blocking(stream.fileno(), False)
            selector.register(stream, selectors.EVENT_READ, name)

        while selector.get_map():
            try:
                events = selector.select(
                    timeout=REPAIRTIMING_SELECT_POLL_SECONDS)
            except InterruptedError:
                continue
            now = time.monotonic()
            if process.poll() is not None and post_exit_deadline is None:
                post_exit_deadline = (
                    now + REPAIRTIMING_POST_EXIT_DRAIN_SECONDS)
            if not events:
                if post_exit_deadline is not None:
                    pipe_holder = True
                    break
                continue
            for key, unused_mask in events:
                name = key.data
                target = buffers[name]
                remaining = limits[name] - len(target)
                read_size = min(
                    REPAIRTIMING_CAPTURE_CHUNK_BYTES, remaining + 1)
                try:
                    chunk = os.read(key.fd, read_size)
                except BlockingIOError:
                    continue
                if not chunk:
                    selector.unregister(key.fileobj)
                    key.fileobj.close()
                    continue
                target.extend(chunk)
                if len(target) > limits[name]:
                    overflow = name
                    break
            if overflow is not None:
                break
            if process.poll() is not None and post_exit_deadline is None:
                post_exit_deadline = (
                    time.monotonic() +
                    REPAIRTIMING_POST_EXIT_DRAIN_SECONDS
                )
            if (
                post_exit_deadline is not None and
                selector.get_map() and
                time.monotonic() >= post_exit_deadline
            ):
                pipe_holder = True
                break

        if overflow is not None:
            raise MeasurementError(
                f"repairtiming benchmark {overflow} exceeds its byte limit "
                f"({limits[overflow]} bytes)")
        if pipe_holder:
            raise MeasurementError(
                "repairtiming benchmark leader exited while a descendant "
                "held an output pipe open")
        try:
            returncode = process.wait(
                timeout=REPAIRTIMING_REAP_TIMEOUT_SECONDS)
        except subprocess.TimeoutExpired as error:
            raise MeasurementError(
                "repairtiming benchmark closed its output pipes before exit"
            ) from error
        unexpected_descendants = _repairtiming_process_group_exists(
            process.pid)
    finally:
        try:
            _terminate_repairtiming_process_group(process)
        finally:
            try:
                if selector is not None:
                    selector.close()
            finally:
                for stream in (process.stdout, process.stderr):
                    if stream is not None:
                        stream.close()

    if unexpected_descendants:
        raise MeasurementError(
            "repairtiming benchmark left an unexpected descendant")
    decoded = {}
    for name, payload in buffers.items():
        try:
            decoded[name] = payload.decode("ascii")
        except UnicodeDecodeError as error:
            raise MeasurementError(
                f"repairtiming benchmark emitted non-ASCII {name}") from error
    stdout = decoded["stdout"]
    stderr = decoded["stderr"]
    diagnostic = _repairtiming_diagnostic(
        stdout, stderr, diagnostic_limit)
    if returncode != 0:
        raise MeasurementError(
            f"repairtiming benchmark {command_label} exited {returncode}: "
            f"{diagnostic}")
    if stderr:
        raise MeasurementError(
            "repairtiming benchmark wrote unexpected stderr: "
            f"{diagnostic}")
    return stdout


def repairtiming_probe(
        bench, block_count, block_bytes, repair_arm, *,
        repair_policy=REPAIR_V1_POLICY_NAME,
        dispatch_profile=peel_codec.TARGET_PROFILE,
        construction_seed, construction_seed_derivation, loss, loss_seed,
        schedule, warmup_replicates, replicates, inner_reps, max_overhead,
        cache_state, systematic_cache, evict_bytes, context,
        required_margin):
    """Run native ``repairtiming`` while retaining complete thermal evidence."""
    arm = globals()["repair_arm"](repair_arm)
    validate_repairtiming_dimensions(
        block_count=block_count,
        block_bytes=block_bytes,
        repair_arm=arm,
        repair_policy=repair_policy,
        dispatch_profile=dispatch_profile,
        construction_seed=construction_seed,
        construction_seed_derivation=construction_seed_derivation,
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
            frozen_context,
            cache_state=cache_state,
            evict_bytes=evict_bytes,
            require_capture=True,
        )
        derivation_cli = (
            "fixed"
            if construction_seed_derivation ==
                REPAIRTIMING_CONSTRUCTION_SEED_FIXED
            else "derived"
        )
        command = [
            bench,
            "repairtiming",
            "--N", str(block_count),
            "--bb", str(block_bytes),
            "--repair-arm", arm.name,
            "--repair-policy", repair_policy,
            "--dispatch-profile", dispatch_profile,
            "--construction-seed", str(construction_seed),
            "--construction-seed-derivation", derivation_cli,
            "--loss", f"{float(loss):.17g}",
            "--loss-seed", str(loss_seed),
            "--schedule", schedule,
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
        stdout = _run_repairtiming_checked(
            command, peel_codec.isolated_codec_env())
        if len(stdout) > REPAIRTIMING_NATIVE_STDOUT_BYTE_LIMIT:
            raise MeasurementError(
                "repairtiming output exceeds its byte limit")
        final_context = peel_codec._finalize_paired_context(
            frozen_context,
            stdout,
            capture,
            run_bounds=_repairtiming_run_bounds,
        )
    finally:
        os.close(capture.fd)
    return parse_repairtiming_output(
        stdout,
        block_count=block_count,
        block_bytes=block_bytes,
        repair_arm=arm,
        repair_policy=repair_policy,
        dispatch_profile=dispatch_profile,
        construction_seed=construction_seed,
        construction_seed_derivation=construction_seed_derivation,
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


def replay_repairtiming_receipt(receipt, *, expected_request=None):
    """Reparse native text and reject forged or stale normalized fields."""
    receipt_fields = {
        "protocol",
        "manifest",
        "context",
        "stream_sha256",
        "evidence",
        "native_stdout",
    }
    if not isinstance(receipt, dict) or set(receipt) != receipt_fields:
        raise MeasurementError("repairtiming receipt schema is incomplete")
    if receipt.get("protocol") != REPAIRTIMING_PROTOCOL_V3:
        raise MeasurementError("repairtiming receipt protocol is invalid")
    manifest = receipt.get("manifest")
    if not isinstance(manifest, dict):
        raise MeasurementError("repairtiming receipt manifest is invalid")
    try:
        arm = repair_arm(manifest.get("repair_arm"))
    except ValueError as error:
        raise MeasurementError(f"repairtiming receipt arm is invalid: {error}")
    request = {
        "block_count": manifest.get("K"),
        "block_bytes": manifest.get("bb"),
        "repair_arm": arm,
        "repair_policy": REPAIR_V1_POLICY_NAME,
        "dispatch_profile": manifest.get("dispatch_profile"),
        "construction_seed": manifest.get("construction_seed_base"),
        "construction_seed_derivation":
            manifest.get("construction_seed_derivation"),
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
            raise MeasurementError(
                "expected repairtiming request is not a dict")
        for name, expected in expected_request.items():
            if name not in request:
                raise MeasurementError(
                    f"repairtiming replay request has no {name}")
            actual = request[name]
            if name == "repair_arm":
                try:
                    same = repair_arm(expected) == actual
                except ValueError:
                    same = False
            else:
                same = peel_codec._same_typed_json(actual, expected)
            if not same:
                raise MeasurementError(
                    f"repairtiming replay request changed {name}")
    try:
        parsed = parse_repairtiming_output(
            receipt.get("native_stdout"), **request)
    except ValueError as error:
        raise MeasurementError(
            f"repairtiming receipt request is invalid: {error}")
    if not peel_codec._same_typed_json(parsed.as_dict(), receipt):
        raise MeasurementError(
            "repairtiming receipt has forged or stale derived fields")
    return parsed


class _ByteLimitExceeded(Exception):
    pass


class _CappedHashReader:
    def __init__(self, source, limit):
        self.source = source
        self.limit = limit
        self.count = 0
        self.digest = hashlib.sha256()

    def read(self, size):
        remaining = self.limit - self.count
        if remaining < 0:
            raise _ByteLimitExceeded
        request = min(size, remaining + 1)
        data = self.source.read(request)
        self.count += len(data)
        if self.count > self.limit:
            raise _ByteLimitExceeded
        self.digest.update(data)
        return data


def read_repairtiming_gzip_text(
        path, *,
        compressed_limit=REPAIRTIMING_GZIP_COMPRESSED_BYTE_LIMIT,
        decompressed_limit=REPAIRTIMING_GZIP_DECOMPRESSED_BYTE_LIMIT):
    """Read exactly one stable gzip member under hard byte limits."""
    if (
        type(compressed_limit) is not int or compressed_limit < 1 or
        type(decompressed_limit) is not int or decompressed_limit < 1
    ):
        raise ValueError("repairtiming gzip byte limits are invalid")
    stable_fields = (
        "st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns",
        "st_ctime_ns",
    )
    reader = None
    try:
        nofollow = getattr(os, "O_NOFOLLOW", 0)
        if not nofollow:
            raise MeasurementError(
                "repairtiming gzip requires no-follow file opens")

        def open_nofollow(raw_path, flags):
            return os.open(
                raw_path,
                flags | nofollow | getattr(os, "O_CLOEXEC", 0),
            )

        with open(path, "rb", opener=open_nofollow) as source:
            before = os.fstat(source.fileno())
            if not stat.S_ISREG(before.st_mode):
                raise MeasurementError(
                    f"repairtiming gzip is not a regular file: {path!r}")
            if before.st_size > compressed_limit:
                raise MeasurementError(
                    "repairtiming gzip exceeds its compressed byte limit")
            reader = _CappedHashReader(source, compressed_limit)
            archive = zlib.decompressobj(16 + zlib.MAX_WBITS)
            chunks = []
            total = 0
            while not archive.eof:
                compressed = reader.read(64 * 1024)
                if not compressed:
                    raise MeasurementError("repairtiming gzip is truncated")
                pending = compressed
                while pending and not archive.eof:
                    remaining_with_probe = decompressed_limit - total + 1
                    chunk = archive.decompress(
                        pending,
                        min(64 * 1024, remaining_with_probe),
                    )
                    pending = archive.unconsumed_tail
                    total += len(chunk)
                    if total > decompressed_limit:
                        raise MeasurementError(
                            "repairtiming gzip exceeds its decompressed "
                            "byte limit")
                    chunks.append(chunk)
                if archive.unused_data:
                    raise MeasurementError(
                        "repairtiming gzip has trailing data")
            if reader.read(1):
                raise MeasurementError(
                    "repairtiming gzip has trailing data")
            after = os.fstat(source.fileno())
        current = os.stat(path, follow_symlinks=False)
        if (
            not stat.S_ISREG(current.st_mode) or
            any(getattr(before, name) != getattr(after, name)
                for name in stable_fields) or
            any(getattr(after, name) != getattr(current, name)
                for name in stable_fields) or
            reader.count != after.st_size
        ):
            raise MeasurementError(
                "repairtiming gzip changed while being read")
        payload = b"".join(chunks)
    except _ByteLimitExceeded:
        raise MeasurementError(
            "repairtiming gzip exceeds its compressed byte limit")
    except (OSError, zlib.error) as error:
        raise MeasurementError(f"could not read repairtiming gzip: {error}")
    try:
        text = payload.decode("ascii")
    except UnicodeDecodeError:
        raise MeasurementError("repairtiming gzip payload is not ASCII")
    return text, {
        "compressed_sha256": reader.digest.hexdigest(),
        "decompressed_sha256": hashlib.sha256(payload).hexdigest(),
        "compressed_bytes": reader.count,
        "decompressed_bytes": len(payload),
    }
