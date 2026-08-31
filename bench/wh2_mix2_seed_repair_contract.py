#!/usr/bin/env python3
"""Freeze, derive, validate, and export the WH2 Two07/mix2 repair map.

This is a separate v9 bench contract.  It does not reinterpret the rejected
v5 contract, the uniform-raw descriptors used there, or the closed MIX2
decision.
The native worker performs exact byte recovery.  Derivation searches attempts
only offline; validation sends exactly one already-frozen uint8 attempt for
each K and has no retry command.
"""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
import os
from pathlib import Path
import re
import selectors
import signal
import stat
import subprocess
import sys
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Any, Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Tuple, Union

try:
    from bench import wh2_mix2_promotion_short_screen as short_screen_controller
except ModuleNotFoundError as exc:  # Direct execution via python3 bench/<script>.
    if exc.name != "bench":
        raise
    import wh2_mix2_promotion_short_screen as short_screen_controller


CONTRACT_SCHEMA = "wirehair.wh2.mix2-seed-repair-contract.v9"
PROFILE_SCHEMA = "wirehair.wh2.mix2-production-profile.v1"
DESCRIPTION_SCHEMA = "wirehair.wh2.mix2-seed-repair-worker-description.v3"
WORKER_SCHEMA = "wirehair.wh2.mix2-seed-repair-worker.v3"
DERIVATION_RECORD_SCHEMA = \
    "wirehair.wh2.mix2-seed-repair-derivation-record.v3"
VALIDATION_RECORD_SCHEMA = \
    "wirehair.wh2.mix2-seed-repair-validation-record.v3"
MANIFEST_SCHEMA = "wirehair.wh2.mix2-seed-repair-manifest.v2"
FREEZE_SCHEMA = "wirehair.wh2.mix2-seed-repair-freeze.v2"
MAP_SCHEMA = "wirehair.wh2.mix2-repair-map.v2"
DERIVATION_SUMMARY_SCHEMA = \
    "wirehair.wh2.mix2-seed-repair-derivation-summary.v2"
VALIDATION_SUMMARY_SCHEMA = \
    "wirehair.wh2.mix2-seed-repair-validation-summary.v2"
MAP_EXPORT_SCHEMA = "wirehair.wh2.mix2-repair-map-cxx-include.v2"
SEMANTIC_REPLAY_SCHEMA = \
    "wirehair.wh2.mix2-repair-map-semantic-replay-roster.v2"

SHORT_SCREEN_CONTRACT_SCHEMA = \
    "wirehair.wh2.mix2-promotion-short-screen.v2"
SHORT_SCREEN_INPUT_SCHEMA = \
    "wirehair.wh2.mix2-promotion-short-screen-input.v2"
SHORT_SCREEN_RESULT_SCHEMA = \
    "wirehair.wh2.mix2-promotion-short-screen-record.v2"
SHORT_SCREEN_SUMMARY_SCHEMA = \
    "wirehair.wh2.mix2-promotion-short-screen-summary.v2"
SHORT_SCREEN_ATTEMPT_STREAM_SCHEMA = \
    "wirehair.wh2.mix2-promotion-short-screen-attempt-stream.v2"
# Frozen from bench/wh2_mix2_promotion_short_screen.py before any all-K run.
SHORT_SCREEN_CONTRACT_SHA256 = \
    "28e6affc80680377df2e8099825b1b475b65e51d04116a8c2fd6465566041cdb"
SHORT_SCREEN_ROOTS = (
    "0xc4695292d9835286",
    "0x7ccd510f122fc160",
    "0x7001a960b7d9c0a4",
)
SHORT_SCREEN_ARMS = (
    "current_disabled_mix3",
    "candidate_two07_mix2",
)
SHORT_SCREEN_TIMING_ORDERS = ("ABBA", "BAAB")
SHORT_SCREEN_OBSERVATIONS_PER_ARM = 2
SHORT_SCREEN_OFFICIAL_SCOPE = \
    "v9 production-basis bounded promotion screen"

DEFAULT_CONTRACT = Path(__file__).with_name(
    "wh2_mix2_seed_repair_contract_v9.json")
CONTROLLER_PATH = Path(__file__).resolve()
REPO_ROOT = CONTROLLER_PATH.parent.parent
K_MIN = 2
K_MAX = 64000
MAP_BYTES = K_MAX - K_MIN + 1
SELECTION_CELL_COUNT = 27
VALIDATION_CELL_COUNT = 9
SHORT_SCREEN_ATTEMPT_COUNT = 30
SHORT_SCREEN_CELL_COUNT = 270
SHORT_SCREEN_RECORD_COUNT = (
    SHORT_SCREEN_CELL_COUNT * len(SHORT_SCREEN_ARMS) *
    SHORT_SCREEN_OBSERVATIONS_PER_ARM)
MAX_JOBS = 64
MAX_RECORD_BYTES = 256 * 1024
MAX_AUDIT_BYTES = MAP_BYTES * MAX_RECORD_BYTES
MAX_STDERR_BYTES = 64 * 1024
MAX_DESCRIPTION_BYTES = 64 * 1024
MAX_DEADLINE_SECONDS = 7 * 24 * 60 * 60
SCIENTIFIC_REJECT_EXIT_CODE = 2
SEMANTIC_REPLAY_DEADLINE_SECONDS = 30 * 60
MASK32 = (1 << 32) - 1
MASK64 = (1 << 64) - 1
SHA256 = re.compile(r"[0-9a-f]{64}\Z")
COMMIT = re.compile(r"[0-9a-f]{40}\Z")
SHORT_SCREEN_DECIMAL = re.compile(r"(?:0|[1-9][0-9]*)\.[0-9]{3}\Z")

PRODUCTION_SEED_BASIS = "production-normalized-b2-v1"
PRODUCTION_SEED_SCHEDULE_CANONICAL = (
    "wirehair:wh2:two07-mix2:graph-b2:seed-schedule-v1;"
    "profile=SelectSeedProfile(K,2);"
    "precode=MatrixSeedFromProfile(profile,0,0x763263707265636f);"
    "packet=PacketPeelSeedFromProfile(profile,0x76327265636f7665);"
    "precode_stride=0x9e3779b97f4a7c15;"
    "packet_stride=0x9e3779b9"
)
PRODUCTION_SEED_SCHEDULE_SHA256 = (
    "fe8101781024b4a30797d66cb1512640"
    "e6a939586ad5e32f73c5b7ff6f411294"
)
PRECODE_ATTEMPT_STRIDE = 0x9E3779B97F4A7C15
PACKET_ATTEMPT_STRIDE = 0x9E3779B9
SELECTION_ROOTS = (
    "0xd1b54a32d192ed03",
    "0x94d049bb133111eb",
    "0x8538ecb5bd456ea3",
    "0xc0ac29b7c97c50dd",
    "0x3f84d5b5b5470917",
    "0x9216d5d98979fb1b",
    "0xb889883a79549774",
    "0xb5666de0987896af",
    "0x8bfca269b0bc01e0",
)
VALIDATION_ROOTS = (
    "0xefd20c982041a46b",
    "0x8827bc36ed906555",
    "0x86029f23d6132efa",
)
VALIDATION_ROOT_NAMESPACE = \
    "wirehair2-two07-mix2-graph-b2-all-k-v9:holdout-root:"
VALIDATION_ROOT_FULL_SHA256 = (
    "efd20c982041a46b00dc0db6b4117e52e6cdc1265dc05e51b4d4c79705327196",
    "8827bc36ed90655550b226632b470105958ff2ca85ad65cd9b55ba729c4f0a94",
    "86029f23d6132efa80a6cf4b2484f8e789eae1fc76a4f4a6e4144b824b958081",
)
SCHEDULES = ("burst", "adversarial", "repair-only")
VALIDATION_ROSTER_SCHEMA = \
    "wirehair.wh2.mix2-seed-repair-validation-roster.v1"
VALIDATION_ROSTER_CELL_ORDER = "root-major-then-schedule"
VALIDATION_ROSTER_SHA256 = \
    "27bbf2b02f51da35527d47c6313a8a0ea4f4b1c1483c636341ce607ef5fb8581"
EXPECTED_VALIDATION_ROSTER = {
    "schema": VALIDATION_ROSTER_SCHEMA,
    "cell_order": VALIDATION_ROSTER_CELL_ORDER,
    "root_count": len(VALIDATION_ROOTS),
    "schedule_count": len(SCHEDULES),
    "cell_count": VALIDATION_CELL_COUNT,
    "sha256": VALIDATION_ROSTER_SHA256,
}
SEMANTIC_REPLAY_K = (
    2, 3, 4, 5, 6, 8, 16, 32, 64, 100, 101, 128, 256, 512, 513,
    1000, 1001, 2048, 4096, 5000, 5001, 8192, 10000, 10001, 16384,
    20000, 20001, 32768, 49152, 64000,
)
SEMANTIC_REPLAY_ROSTER_SHA256 = (
    "1cbff6985440e83668b18ab1cb8970ba"
    "998ef1ba0522cf6bb72df89d73a74ae1"
)
SELECTION_RULE = (
    "for each K, test uint8 construction attempts in ascending order; reject "
    "an attempt at its first failing cell in frozen cell order; select the "
    "first attempt that byte-recovers at overhead 0 on all twenty-seven cells"
)
LOWER_ATTEMPT_EVIDENCE = (
    "retain exactly one first-failure witness for every lower attempt and all "
    "twenty-seven success receipts for the selected attempt"
)
VALIDATION_RULE = (
    "replay exactly the frozen map byte indexed by K-2 on all nine disjoint "
    "cells; runtime search, retries, substitutions, and post-hoc map edits "
    "are forbidden"
)
VALIDATION_GATE = (
    "zero overhead-0 weak K values and zero unsupported, construction, fatal, "
    "internal, or byte-mismatch outcomes"
)
MAP_FREEZE_RULE = (
    "SHA-256 the exact 63999 bytes and complete canonical derivation audit "
    "before any validation worker starts"
)
BUNDLE_CHAIN_RULE = (
    "the derivation manifest binds an authenticated v9 short-screen PASS; the "
    "validation manifest binds the repair-map receipt; that receipt binds the "
    "short-screen hashes, derivation manifest, freeze, summary, audit, and "
    "map; every self-hash and semantic relationship is rechecked before "
    "validation"
)
EXPECTED_APPLICABILITY = {
    "equation_scope":
        "only the exact candidate_profile object, including its K-indexed "
        "SelectSeedProfile(K,2) graph basis, before uint8 attempt stepping",
    "forbidden_inference":
        "the repair map is not evidence for the legacy width-sensitive "
        "profile, another dense-anchor layout, or another packet mix count",
}
EXPECTED_VALIDATION_ROOT_DERIVATION = {
    "namespace": VALIDATION_ROOT_NAMESPACE,
    "indices": [0, 1, 2],
    "hash": "SHA-256",
    "extraction":
        "first 8 digest bytes interpreted as unsigned big-endian uint64",
    "full_sha256": list(VALIDATION_ROOT_FULL_SHA256),
    "collision_rule":
        "fail if any root is zero, duplicated, or appears in selection roots; "
        "do not retry or substitute",
}
EXPECTED_PACKET_TRACE = {
    "loss_rate_conversion": "double(loss_ppm) / 1000000.0",
    "trace_seed":
        "loss_seed xor (K * 0x9e3779b97f4a7c15) xor (block_bytes * "
        "0xbf58476d1ce4e5b9) mod 2^64",
    "rng":
        "SplitMix64; Unit=(Next() >> 11) / 2^53; drop iff Unit < loss_rate",
    "burst":
        "candidate packet IDs 0,1,...; RNG seed trace_seed xor 0x10fade; "
        "eight-drop bursts start on an idle candidate with probability "
        "loss_rate/(8-7*loss_rate)",
    "adversarial":
        "candidate packet IDs UINT32_MAX-2*i modulo 2^32; RNG seed "
        "trace_seed xor 0x10fade and IID drops",
    "repair-only":
        "candidate packet IDs K+i modulo 2^32; RNG seed trace_seed xor "
        "0x10fade and IID drops",
    "nested_overhead":
        "generate one K+4 delivered-ID prefix; score OH0 from exactly its "
        "first K IDs",
    "candidate_limit": "256*(K+4)+65536 attempted candidates",
    "trace_sha256":
        "SHA-256 of all K+4 delivered IDs encoded consecutively as unsigned "
        "32-bit little-endian words",
}
OUTCOMES = frozenset((
    "success", "need_more_at_oh0", "need_more_at_cap",
    "construct_failed", "unsupported", "fatal",
))

DERIVATION_RECORD_FIELDS = frozenset((
    "K", "base_packet_seed", "base_precode_seed",
    "candidate_profile_sha256", "effective_packet_seed",
    "effective_precode_seed", "lower_attempt_failure_witnesses", "mode",
    "ordinal", "schema", "selected_attempt", "selected_successes",
    "source_sha256", "worker_binary_sha256",
))
VALIDATION_RECORD_FIELDS = frozenset((
    "K", "all_oh0_success", "base_packet_seed", "base_precode_seed",
    "candidate_profile_sha256", "cells", "construction_attempt",
    "effective_packet_seed",
    "effective_precode_seed", "mode", "ordinal", "schema",
    "source_sha256", "worker_binary_sha256",
))
CELL_FIELDS = frozenset((
    "attempted_candidates", "cell_ordinal", "construction_attempt",
    "decoded_extra", "loss_ppm", "loss_root", "outcome", "result_code",
    "root_index", "schedule", "trace_sha256",
))
DESCRIPTION_FIELDS = frozenset((
    "binary_sha256", "candidate_profile", "candidate_profile_sha256",
    "contract_schema", "derivation_schema", "protocol", "schema",
    "source_git_commit", "validation_roster_schema",
    "validation_roster_sha256", "validation_schema", "worker_schema",
))
SHORT_SCREEN_BINDING_FIELDS = (
    "short_screen_contract_sha256",
    "short_screen_input_sha256",
    "short_screen_attempt_stream_sha256",
    "short_screen_result_stream_sha256",
    "short_screen_summary_sha256",
    "short_screen_map_sha256",
)
SHORT_SCREEN_INPUT_FIELDS = frozenset((
    "schema", "contract", "contract_sha256", "source_receipt",
    "source_receipt_sha256", "bench_binary", "repair_worker_binary",
    "repair_worker_description", "repair_worker_description_sha256",
    "attempt_selection_worker_command_sha256",
    "attempt_selection_stream_schema", "attempt_selection_stream_sha256",
    "attempt_selection_record_count", "selected_attempts",
    "worker_validation_commands_issued",
    "worker_final_validation_stream_present", "invocations", "input_sha256",
))
SHORT_SCREEN_SUMMARY_FIELDS = frozenset((
    "schema", "contract_sha256", "input_sha256", "input_file_sha256",
    "attempt_selection_stream_schema", "attempt_selection_stream_sha256",
    "attempt_selection_record_count", "worker_validation_commands_issued",
    "worker_final_validation_stream_present", "result_stream_sha256",
    "record_count", "invocation_count", "bench_binary_sha256",
    "repair_worker_binary_sha256", "source_git_commit", "source_receipt",
    "source_receipt_sha256", "candidate_profile_sha256",
    "architecture_selection_performed", "offline_attempt_derivation_performed",
    "official_scope", "all_K_recovery_claimed", "wirehair1_timing_claimed",
    "disposition", "gates", "candidate_weak_observations",
    "control_weak_observations", "candidate_weak_coordinates",
    "control_weak_coordinates", "common_success_cells",
    "common_success_timing_pairs", "aggregates", "ratios",
    "complexity_receipt", "summary_sha256",
))
SHORT_SCREEN_BINARY_RECEIPT_FIELDS = frozenset((
    "path", "sha256", "size",
))
SHORT_SCREEN_SOURCE_RECEIPT_FIELDS = frozenset((
    "source_git_commit", "tracked_and_untracked_linked_sources_clean",
    "clean_status_scope", "source_files", "source_receipt_sha256",
))
SHORT_SCREEN_SOURCE_FILES = (
    "CMakeLists.txt",
    "bench/wh2_mix2_promotion_short_screen.py",
    "bench/Wh2Mix2SeedRepairWorker.cpp",
    "bench/Wh2FrozenTrace.cpp",
    "bench/Wh2FrozenTrace.h",
    "bench/Wh2NativeCodec.cpp",
    "bench/Wh2NativeCodec.h",
    "codec/WirehairV2Bench.cpp",
    "codec/WirehairV2Precode.cpp",
    "codec/WirehairV2Precode.h",
    "codec/WirehairV2PrecodeEncode.cpp",
    "codec/WirehairV2PrecodeEncode.h",
    "codec/WirehairV2Profile.cpp",
    "codec/WirehairV2Seeds.h",
    "include/wirehair/wirehair.h",
    "WirehairDenseFixups.inc",
    "WirehairPeelFixups.inc",
)
SHORT_SCREEN_CLEAN_STATUS_SCOPE = (
    "CMakeLists.txt", "cmake", "codec", "include", "tables",
    "wirehair.cpp", "gf256.cpp", "gf256.h", "WirehairCodec.cpp",
    "WirehairCodec.h", "WirehairCompiler.h", "WirehairEnvironment.h",
    "WirehairHeavy.h", "WirehairTools.cpp", "WirehairTools.h",
    "WirehairDenseFixups.inc", "WirehairPeelFixups.inc",
    "bench/Wh2Mix2SeedRepairWorker.cpp", "bench/Wh2FrozenTrace.cpp",
    "bench/Wh2FrozenTrace.h", "bench/Wh2NativeCodec.cpp",
    "bench/Wh2NativeCodec.h",
    "bench/wh2_mix2_promotion_short_screen.py",
)
SHORT_SCREEN_GATE_FIELDS = frozenset((
    "candidate_benchmark_weak_observations_zero",
    "candidate_benchmark_weak_observations_not_above_control",
    "repeated_deterministic_work_equal",
    "block_xor_ratio_at_most_0.9829453304",
    "gf256_muladd_ratio_at_most_1",
    "inactivation_ratio_at_most_1",
    "aggregate_b2_solve_time_ratio_at_most_1",
    "aggregate_b2_build_time_ratio_at_most_1.05",
))
SHORT_SCREEN_VERDICT_FIELDS = frozenset((
    "disposition", "gates", "candidate_weak_observations",
    "control_weak_observations", "candidate_weak_coordinates",
    "control_weak_coordinates", "common_success_cells",
    "common_success_timing_pairs", "aggregates", "ratios",
    "complexity_receipt",
))
SHORT_SCREEN_INVOCATION_FIELDS = frozenset((
    "ordinal", "coordinate_ordinal", "K", "root_index", "loss_root",
    "schedule", "cell_ordinal", "arm", "construction_attempt",
    "timing_order", "timing_slot", "observation_index", "argv",
    "command_sha256",
))
SHORT_SCREEN_RESULT_FIELDS = frozenset((
    "schema", "ordinal", "coordinate_ordinal", "arm", "K", "block_bytes",
    "loss_ppm", "overhead", "root_index", "loss_root", "schedule",
    "cell_ordinal", "timing_order", "timing_slot", "observation_index",
    "attempt_selection_stream_sha256",
    "attempt_selection_cell_receipts_used_as_promotion_evidence",
    "benchmark_loss_trace_hash_recorded",
    "full_payload_byte_recovery_verified", "candidate_profile_sha256",
    "construction_attempt", "attempt_mode", "effective_precode_seed",
    "effective_packet_seed", "staircase", "binary_dense_rows",
    "gf256_heavy_rows", "source_hits", "dense_anchor_layout", "mix_count",
    "success", "rank_fail", "error", "weak", "inactivated_columns",
    "block_xors", "gf256_muladds", "solve_ms", "build_ms", "peel_ms",
    "project_ms", "residual_ms", "backsub_ms",
    "extra_dense_basis_capacity_entries",
    "extra_dense_basis_capacity_bytes", "command_sha256",
    "bench_stdout_sha256", "bench_binary_sha256",
    "bench_source_git_commit", "source_receipt_sha256",
))
MAP_RECEIPT_FIELDS = frozenset((
    "schema", "contract_sha256", "candidate_arm",
    "candidate_profile_sha256", "construction_seed_basis",
    "seed_schedule_sha256", "entry_kind", "map_bytes", "map_sha256",
    "derivation_audit_records",
    "derivation_audit_sha256", "worker_binary_sha256", "source_git_commit",
    "controller_sha256", "selection_cells_per_selected_attempt",
    "selection_rule", "derivation_manifest_sha256",
    "derivation_freeze_sha256", "derivation_summary_sha256",
    "receipt_sha256",
)) | frozenset(SHORT_SCREEN_BINDING_FIELDS)
MANIFEST_FIELDS = frozenset((
    "schema", "action", "contract_sha256", "candidate_arm",
    "candidate_profile_sha256", "construction_seed_basis",
    "seed_schedule_sha256", "controller_sha256", "worker_binary_sha256",
    "source_git_commit", "K_min", "K_max", "map_bytes",
    "selection_roots", "validation_roots", "schedules",
    "selection_cell_count", "validation_cell_count",
    "validation_roster_sha256",
    "repair_map_sha256", "map_receipt_sha256", "manifest_sha256",
)) | frozenset(SHORT_SCREEN_BINDING_FIELDS)
FREEZE_FIELDS = frozenset((
    "schema", "action", "manifest_sha256", "manifest", "freeze_sha256",
))
DERIVATION_SUMMARY_FIELDS = frozenset((
    "schema", "contract_sha256", "candidate_profile_sha256",
    "manifest_sha256", "audit_sha256", "audit_records", "map_sha256",
    "map_bytes", "maximum_selected_attempt", "selected_attempt_histogram",
    "runtime_search", "summary_sha256",
))
VALIDATION_SUMMARY_FIELDS = frozenset((
    "schema", "contract_sha256", "candidate_profile_sha256",
    "manifest_sha256", "repair_map_sha256", "map_receipt_sha256",
    "audit_sha256", "audit_records", "weak_K_count", "weak_K",
    "runtime_search", "disposition", "summary_sha256",
))

DERIVATION_FREEZE_NAME = "derivation-freeze.json"
DERIVATION_AUDIT_NAME = "derivation-audit.jsonl"
MAP_NAME = "repair-map.bin"
MAP_RECEIPT_NAME = "repair-map.json"
DERIVATION_SUMMARY_NAME = "derivation-summary.json"
VALIDATION_FREEZE_NAME = "validation-freeze.json"
VALIDATION_AUDIT_NAME = "validation-audit.jsonl"
VALIDATION_SUMMARY_NAME = "validation-summary.json"
SHORT_SCREEN_INPUT_NAME = "promotion-short-screen-input.json"
SHORT_SCREEN_ATTEMPT_NAME = "promotion-short-screen-attempts.jsonl"
SHORT_SCREEN_RESULT_NAME = "promotion-short-screen-results.jsonl"
SHORT_SCREEN_SUMMARY_NAME = "promotion-short-screen-summary.json"


class ContractError(RuntimeError):
    """The v9 contract or one of its bound artifacts is invalid."""


def fail(message: str) -> None:
    raise ContractError(message)


def canonical_json(value: Any) -> str:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False)


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def sha256_json(value: Any) -> str:
    return sha256_bytes(canonical_json(value).encode("ascii"))


def validation_roster_sha256(
        roots: Sequence[str] = VALIDATION_ROOTS,
        schedules: Sequence[str] = SCHEDULES) -> str:
    """Digest the exact root-major V roster before any holdout is probed."""
    body = {
        "schema": VALIDATION_ROSTER_SCHEMA,
        "cell_order": VALIDATION_ROSTER_CELL_ORDER,
        "root_count": len(roots),
        "roots": list(roots),
        "schedule_count": len(schedules),
        "schedules": list(schedules),
        "cell_count": len(roots) * len(schedules),
    }
    return sha256_json(body)


def _same_frozen_json(value: Any, expected: Any) -> bool:
    """Compare JSON values without Python's bool/int/float aliases."""
    return type(value) is type(expected) and \
        canonical_json(value) == canonical_json(expected)


def _semantic_replay_roster_sha256() -> str:
    if (not SEMANTIC_REPLAY_K or
            tuple(sorted(set(SEMANTIC_REPLAY_K))) != SEMANTIC_REPLAY_K or
            any(K < K_MIN or K > K_MAX for K in SEMANTIC_REPLAY_K)):
        fail("semantic replay roster is not the frozen 30-K boundary set")
    digest = sha256_json({
        "schema": SEMANTIC_REPLAY_SCHEMA,
        "K": list(SEMANTIC_REPLAY_K),
    })
    if digest != SEMANTIC_REPLAY_ROSTER_SHA256:
        fail("semantic replay roster differs from its frozen digest")
    return digest


def _shared_coordinate_identity(value: Mapping[str, Any]) -> bytes:
    """Canonical attempt-independent identity shared by both all-K audits."""
    identity = {
        "K": value["K"],
        "base_packet_seed": value["base_packet_seed"],
        "base_precode_seed": value["base_precode_seed"],
        "source_sha256": value["source_sha256"],
    }
    return (canonical_json(identity) + "\n").encode("ascii")


def _object_without_self_hash(value: Mapping[str, Any], field: str) \
        -> Mapping[str, Any]:
    result = dict(value)
    result.pop(field, None)
    return result


def _reject_duplicate_pairs(pairs: Sequence[Tuple[str, Any]]) \
        -> Mapping[str, Any]:
    result: Dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            fail("JSON contains duplicate key {!r}".format(key))
        result[key] = value
    return result


def _parse_json_bytes(data: bytes, context: str) -> Any:
    try:
        text = data.decode("utf-8")
        return json.loads(
            text, object_pairs_hook=_reject_duplicate_pairs,
            parse_constant=lambda token: fail(
                "{} contains {}".format(context, token)))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        fail("{} is not strict UTF-8 JSON: {}".format(context, exc))
    return None


def _read_regular(path: Path, cap: int, context: str) -> bytes:
    descriptor = -1
    try:
        nofollow = getattr(os, "O_NOFOLLOW", 0)
        if nofollow == 0:
            fail("{} cannot be opened fail-closed".format(context))
        descriptor = os.open(
            str(path), os.O_RDONLY | os.O_CLOEXEC | nofollow)
        before = os.fstat(descriptor)
        if (not stat.S_ISREG(before.st_mode) or before.st_size < 0 or
                before.st_size > cap):
            fail("{} has an invalid type or size".format(context))
        chunks = []
        remaining = before.st_size
        while remaining:
            chunk = os.read(descriptor, min(1024 * 1024, remaining))
            if not chunk:
                fail("{} was truncated".format(context))
            chunks.append(chunk)
            remaining -= len(chunk)
        if os.read(descriptor, 1):
            fail("{} grew while reading".format(context))
        after = os.fstat(descriptor)
        identity = lambda value: (
            value.st_dev, value.st_ino, value.st_size,
            value.st_mtime_ns, value.st_ctime_ns)
        if identity(before) != identity(after):
            fail("{} changed while reading".format(context))
        return b"".join(chunks)
    except OSError as exc:
        fail("cannot read {}: {}".format(context, exc))
    finally:
        if descriptor >= 0:
            os.close(descriptor)
    return b""


def _load_canonical_object(path: Path, cap: int, context: str) \
        -> Mapping[str, Any]:
    data = _read_regular(path, cap, context)
    value = _parse_json_bytes(data, context)
    if not isinstance(value, dict) or \
            data != (canonical_json(value) + "\n").encode("ascii"):
        fail("{} is not one canonical JSON object".format(context))
    return value


def load_contract(path: Path = DEFAULT_CONTRACT) -> Mapping[str, Any]:
    data = _read_regular(path, 64 * 1024, "mix2 v9 contract")
    value = _parse_json_bytes(data, "mix2 v9 contract")
    if not isinstance(value, dict) or value.get("schema") != CONTRACT_SCHEMA:
        fail("contract schema is not the closed v9 identity")
    expected_top = frozenset((
        "schema", "contract_id", "candidate_profile", "applicability",
        "domain", "selection", "packet_trace", "short_screen",
        "validation", "attempt_schedule", "repair_map",
    ))
    if set(value) != expected_top or \
            value.get("contract_id") != \
                "wirehair2-two07-mix2-graph-b2-all-k-v9":
        fail("contract identity or top-level fields changed")
    profile = value.get("candidate_profile")
    expected_profile = {
        "schema": PROFILE_SCHEMA,
        "arm": "wirehair2_two07_mix2_graph_b2_v1",
        "codec": "wirehair2_candidate",
        "binary_dense_rows": 12,
        "gf256_heavy_rows": 12,
        "graph_seed_block_bytes": 2,
        "heavy_family": "periodic-cauchy",
        "source_hits": "certified-by-K",
        "dense_identity_corner": False,
        "dense_anchor_layout": "two07",
        "field": "GF(256)",
        "mix_count": 2,
        "construction_seed_basis": PRODUCTION_SEED_BASIS,
        "precode_seed_salt": "0x763263707265636f",
        "packet_seed_salt": "0x76327265636f7665",
        "seed_schedule_sha256": PRODUCTION_SEED_SCHEDULE_SHA256,
    }
    if not _same_frozen_json(profile, expected_profile):
        fail("contract candidate is not normalized-production Two07/mix2")
    if not _same_frozen_json(
            value.get("applicability"), EXPECTED_APPLICABILITY):
        fail("contract permits an invalid cross-profile inference")
    domain = value.get("domain")
    if not _same_frozen_json(domain, {
            "K_min": K_MIN, "K_max": K_MAX, "map_entries": MAP_BYTES,
            "block_bytes": 2, "overhead": 0, "loss_ppm": 500000,
            "source_seed_derivation":
                "MakeDeterministicSource(K,2,0x6d69783273656564); state is "
                "seed xor (K*0xd6e8feb86659fd93) xor "
                "(block_bytes*0xa0761d6478bd642f), followed by SplitMix64 "
                "words serialized little-endian",
            "cell_order": "phase-root-major-then-schedule",
            "schedules": list(SCHEDULES)}):
        fail("contract all-K OH0 domain changed")
    selection = value.get("selection")
    short_screen = value.get("short_screen")
    validation = value.get("validation")
    if (not _same_frozen_json(selection, {
            "roots": list(SELECTION_ROOTS),
            "cell_count": SELECTION_CELL_COUNT,
            "selection_rule": SELECTION_RULE,
            "lower_attempt_evidence": LOWER_ATTEMPT_EVIDENCE}) or
            not _same_frozen_json(short_screen, {
                "contract_schema": SHORT_SCREEN_CONTRACT_SCHEMA,
                "contract_sha256": SHORT_SCREEN_CONTRACT_SHA256,
                "input_schema": SHORT_SCREEN_INPUT_SCHEMA,
                "attempt_stream_schema":
                    SHORT_SCREEN_ATTEMPT_STREAM_SCHEMA,
                "result_schema": SHORT_SCREEN_RESULT_SCHEMA,
                "summary_schema": SHORT_SCREEN_SUMMARY_SCHEMA,
                "K": list(SEMANTIC_REPLAY_K),
                "attempt_records": SHORT_SCREEN_ATTEMPT_COUNT,
                "attempt_bytes": SHORT_SCREEN_ATTEMPT_COUNT,
                "selection_cells_per_attempt": SELECTION_CELL_COUNT,
                "required_disposition": "PASS",
                "final_holdout_access": False,
                "full_derivation_reproduction":
                    "exact canonical 30-K derivation records and selected "
                    "attempt bytes before validation"}) or
            not _same_frozen_json(validation, {
                "roots": list(VALIDATION_ROOTS),
                "cell_count": VALIDATION_CELL_COUNT,
                "root_derivation": EXPECTED_VALIDATION_ROOT_DERIVATION,
                "roster": EXPECTED_VALIDATION_ROSTER,
                "rule": VALIDATION_RULE,
                "gate": VALIDATION_GATE}) or
            set(SELECTION_ROOTS) & set(VALIDATION_ROOTS)):
        fail("selection, short-screen, or validation identity changed")
    if validation_roster_sha256() != VALIDATION_ROSTER_SHA256:
        fail("validation roster digest differs from its canonical identity")
    derived_validation_roots = []
    for index, expected_digest in enumerate(VALIDATION_ROOT_FULL_SHA256):
        digest = hashlib.sha256(
            (VALIDATION_ROOT_NAMESPACE + str(index)).encode("ascii")
        ).hexdigest()
        if digest != expected_digest:
            fail("validation root full digest changed")
        derived_validation_roots.append("0x" + digest[:16])
    if (tuple(derived_validation_roots) != VALIDATION_ROOTS or
            len(set(SELECTION_ROOTS + VALIDATION_ROOTS)) !=
                len(SELECTION_ROOTS) + len(VALIDATION_ROOTS) or
            any(int(root, 16) == 0
                for root in SELECTION_ROOTS + VALIDATION_ROOTS)):
        fail("validation roots are not their frozen disjoint derivation")
    if not _same_frozen_json(
            value.get("packet_trace"), EXPECTED_PACKET_TRACE):
        fail("packet trace law differs from the frozen v9 contract")
    attempt = value.get("attempt_schedule")
    if not _same_frozen_json(attempt, {
            "attempt_count": 256,
            "graph_seed_profile": "SelectSeedProfile(K,2)",
            "precode_base_seed":
                "MatrixSeedFromProfile(graph_seed_profile,0,"
                "0x763263707265636f)",
            "precode_attempt_stride": "0x9e3779b97f4a7c15",
            "packet_base_seed":
                "PacketPeelSeedFromProfile(graph_seed_profile,"
                "0x76327265636f7665)",
            "packet_attempt_stride": "0x9e3779b9"}):
        fail("normalized-production uint8 attempt schedule changed")
    if hashlib.sha256(
            PRODUCTION_SEED_SCHEDULE_CANONICAL.encode("ascii")).hexdigest() != \
            PRODUCTION_SEED_SCHEDULE_SHA256:
        fail("normalized-production seed schedule digest is inconsistent")
    repair = value.get("repair_map")
    if not _same_frozen_json(repair, {
            "schema": MAP_SCHEMA,
            "encoding": "63999 uint8 attempt values indexed by K-2",
            "freeze_rule": MAP_FREEZE_RULE,
            "bundle_chain": BUNDLE_CHAIN_RULE,
            "runtime_search": False}):
        fail("repair-map identity or no-runtime-search rule changed")
    return value


def contract_sha256(contract: Mapping[str, Any]) -> str:
    return sha256_json(contract)


def candidate_profile_sha256(contract: Mapping[str, Any]) -> str:
    return sha256_json(contract["candidate_profile"])


def _require_hex_seed(value: Any, digits: int, context: str) -> int:
    if (not isinstance(value, str) or len(value) != digits + 2 or
            not value.startswith("0x") or
            re.fullmatch(r"[0-9a-f]+", value[2:]) is None):
        fail("{} is not a fixed-width lowercase seed".format(context))
    return int(value[2:], 16)


def _effective_precode_seed(base: int, attempt: int) -> str:
    return "0x{:016x}".format(
        (base + attempt * PRECODE_ATTEMPT_STRIDE) & MASK64)


def _effective_packet_seed(base: int, attempt: int) -> str:
    return "0x{:08x}".format(
        (base + attempt * PACKET_ATTEMPT_STRIDE) & MASK32)


def _controller_sha256() -> str:
    return sha256_bytes(_read_regular(
        CONTROLLER_PATH, 4 * 1024 * 1024, "v9 controller"))


def _sha256_open_fd(descriptor: int, context: str) -> str:
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            fail("{} is no longer a regular file".format(context))
        digest = hashlib.sha256()
        offset = 0
        while offset < before.st_size:
            chunk = os.pread(
                descriptor, min(1024 * 1024, before.st_size - offset), offset)
            if not chunk:
                fail("{} was truncated while hashing".format(context))
            digest.update(chunk)
            offset += len(chunk)
        if os.pread(descriptor, 1, offset):
            fail("{} grew while hashing".format(context))
        after = os.fstat(descriptor)
        identity = lambda value: (
            value.st_dev, value.st_ino, value.st_size,
            value.st_mtime_ns, value.st_ctime_ns)
        if identity(before) != identity(after):
            fail("{} changed while hashing".format(context))
        return digest.hexdigest()
    except OSError as exc:
        fail("cannot hash {}: {}".format(context, exc))
    return ""


def _require_clean_source(expected_commit: str) -> None:
    commands = (
        ["git", "rev-parse", "--verify", "HEAD^{commit}"],
        ["git", "status", "--porcelain=v1", "--untracked-files=all", "--",
         "CMakeLists.txt", "cmake", "codec", "include", "tables",
         "wirehair.cpp", "gf256.cpp", "gf256.h", "WirehairCodec.cpp",
         "WirehairCodec.h", "WirehairCompiler.h", "WirehairEnvironment.h",
         "WirehairHeavy.h", "WirehairTools.cpp", "WirehairTools.h",
         "WirehairDenseFixups.inc", "WirehairPeelFixups.inc",
         "bench/Wh2Mix2SeedRepairWorker.cpp", "bench/Wh2FrozenTrace.cpp",
         "bench/Wh2FrozenTrace.h", "bench/Wh2NativeCodec.cpp",
         "bench/Wh2NativeCodec.h",
         "bench/wh2_mix2_promotion_short_screen.py",
         "bench/wh2_mix2_seed_repair_contract.py",
         "bench/wh2_mix2_seed_repair_contract_v9.json"],
        ["git", "rev-parse", "--verify", "HEAD^{commit}"],
    )
    outputs = []
    try:
        for command in commands:
            completed = subprocess.run(
                command, cwd=str(REPO_ROOT), stdin=subprocess.DEVNULL,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                timeout=30.0, check=False)
            if completed.returncode != 0:
                fail("cannot bind clean source provenance")
            outputs.append(completed.stdout)
    except (OSError, subprocess.TimeoutExpired) as exc:
        fail("cannot bind clean source provenance: {}".format(exc))
    first = outputs[0].decode("ascii", "strict").strip()
    second = outputs[2].decode("ascii", "strict").strip()
    if first != expected_commit or second != expected_commit or first != second:
        fail("source commit changed or differs from the worker receipt")
    if outputs[1]:
        fail("provenance operation requires a clean relevant source worktree")


def _require_sha(value: Any, context: str) -> str:
    if not isinstance(value, str) or SHA256.fullmatch(value) is None:
        fail("{} is not a lowercase SHA-256".format(context))
    return value


def _require_commit(value: Any, context: str) -> str:
    if not isinstance(value, str) or COMMIT.fullmatch(value) is None:
        fail("{} is not a lowercase 40-hex commit".format(context))
    return value


def _require_int(value: Any, low: int, high: int, context: str) -> int:
    if type(value) is not int or not low <= value <= high:
        fail("{} is outside [{},{}]".format(context, low, high))
    return value


def _expected_cell(cell_ordinal: int, validation: bool) -> Tuple[int, str, str]:
    roots = VALIDATION_ROOTS if validation else SELECTION_ROOTS
    cell_count = VALIDATION_CELL_COUNT if validation else SELECTION_CELL_COUNT
    _require_int(cell_ordinal, 0, cell_count - 1, "cell ordinal")
    root_index = cell_ordinal // len(SCHEDULES)
    schedule = SCHEDULES[cell_ordinal % len(SCHEDULES)]
    return root_index, roots[root_index], schedule


def _validate_cell(
        value: Any, K: int, attempt: int, validation: bool,
        require_success: bool) -> bool:
    if not isinstance(value, dict) or set(value) != CELL_FIELDS:
        fail("cell receipt schema changed")
    cell_count = VALIDATION_CELL_COUNT if validation else SELECTION_CELL_COUNT
    roots = VALIDATION_ROOTS if validation else SELECTION_ROOTS
    ordinal = _require_int(
        value.get("cell_ordinal"), 0, cell_count - 1, "cell ordinal")
    root_index, root, schedule = _expected_cell(ordinal, validation)
    cell_attempt = _require_int(
        value.get("construction_attempt"), 0, 255,
        "cell construction attempt")
    loss_ppm = _require_int(
        value.get("loss_ppm"), 500000, 500000, "cell loss ppm")
    if (_require_int(
            value.get("root_index"), 0, len(roots) - 1, "root index") !=
            root_index or value.get("loss_root") != root or
            value.get("schedule") != schedule or
            loss_ppm != 500000 or
            cell_attempt != attempt):
        fail("cell receipt differs from its frozen coordinate")
    _require_sha(value.get("trace_sha256"), "cell trace")
    maximum_candidates = (K + 4) * 256 + 65536
    _require_int(
        value.get("attempted_candidates"), K + 4, maximum_candidates,
        "attempted candidates")
    _require_int(
        value.get("result_code"), -0x7fffffff, 0x7fffffff, "result code")
    outcome = value.get("outcome")
    if outcome not in OUTCOMES:
        fail("cell outcome is unknown")
    decoded = value.get("decoded_extra")
    result_code = value["result_code"]
    success = outcome == "success"
    if success:
        _require_int(decoded, 0, 0, "OH0 decoded overhead")
        if result_code != 0:
            fail("OH0 success must report decoded_extra=0")
    elif outcome == "need_more_at_oh0":
        _require_int(decoded, 1, 4, "later decoded overhead")
        if result_code != 0:
            fail("later nested success must retain Wirehair_Success")
    elif outcome == "need_more_at_cap":
        if decoded is not None or result_code != 1:
            fail("need-more-at-cap result receipt is inconsistent")
    elif outcome == "construct_failed":
        if decoded is not None or result_code != 4:
            fail("construction-failure result receipt is inconsistent")
    elif decoded is not None:
        fail("non-success cell must use decoded_extra=null")
    if outcome in ("unsupported", "fatal"):
        fail("unsupported/fatal recovery outcomes abort the contract")
    if require_success and not success:
        fail("selected attempt does not pass every OH0 cell")
    return success


def validate_derivation_record(
        value: Any, contract: Mapping[str, Any],
        worker_binary_sha256: str,
        expected_K: Optional[int] = None) -> int:
    if (not isinstance(value, dict) or
            set(value) != DERIVATION_RECORD_FIELDS or
            value.get("schema") != DERIVATION_RECORD_SCHEMA or
            value.get("mode") != "derive"):
        fail("derivation record schema or mode changed")
    K = _require_int(value.get("K"), K_MIN, K_MAX, "derivation K")
    if expected_K is not None and K != expected_K:
        fail("derivation worker returned another K")
    if _require_int(
            value.get("ordinal"), 0, MAP_BYTES - 1,
            "derivation ordinal") != K - K_MIN:
        fail("derivation ordinal is not K-2")
    if (value.get("candidate_profile_sha256") !=
            candidate_profile_sha256(contract) or
            value.get("worker_binary_sha256") != worker_binary_sha256):
        fail("derivation record substitutes profile or worker")
    _require_sha(value.get("source_sha256"), "derivation source")
    selected = _require_int(
        value.get("selected_attempt"), 0, 255, "selected attempt")
    base_precode = _require_hex_seed(
        value.get("base_precode_seed"), 16, "derivation base precode seed")
    base_packet = _require_hex_seed(
        value.get("base_packet_seed"), 8, "derivation base packet seed")
    if (value.get("effective_precode_seed") !=
            _effective_precode_seed(base_precode, selected) or
            value.get("effective_packet_seed") !=
            _effective_packet_seed(base_packet, selected)):
        fail("selected attempt seed receipt is invalid")
    lower = value.get("lower_attempt_failure_witnesses")
    if not isinstance(lower, list) or len(lower) != selected:
        fail("derivation must retain one witness for every lower attempt")
    for attempt, witness in enumerate(lower):
        if _validate_cell(
                witness, K, attempt, validation=False,
                require_success=False):
            fail("lower-attempt witness is not a failure")
    successes = value.get("selected_successes")
    if not isinstance(successes, list) or \
            len(successes) != SELECTION_CELL_COUNT:
        fail("selected attempt must retain exactly twenty-seven successes")
    ordinals = []
    for receipt in successes:
        _validate_cell(
            receipt, K, selected, validation=False, require_success=True)
        ordinals.append(receipt["cell_ordinal"])
    if ordinals != list(range(SELECTION_CELL_COUNT)):
        fail("selected successes are not in frozen cell order")
    for witness in lower:
        selected_receipt = successes[witness["cell_ordinal"]]
        if (witness["trace_sha256"] !=
                selected_receipt["trace_sha256"] or
                witness["attempted_candidates"] !=
                selected_receipt["attempted_candidates"]):
            fail("attempt-independent selection trace identity changed")
    return selected


def validate_validation_record(
        value: Any, contract: Mapping[str, Any],
        worker_binary_sha256: str, expected_K: int,
        expected_attempt: int) -> bool:
    if (not isinstance(value, dict) or
            set(value) != VALIDATION_RECORD_FIELDS or
            value.get("schema") != VALIDATION_RECORD_SCHEMA or
            value.get("mode") != "validate"):
        fail("validation record schema or mode changed")
    K = _require_int(value.get("K"), K_MIN, K_MAX, "validation K")
    if (K != expected_K or
            _require_int(
                value.get("ordinal"), 0, MAP_BYTES - 1,
                "validation ordinal") != K - K_MIN):
        fail("validation worker returned another coordinate")
    attempt = _require_int(
        value.get("construction_attempt"), 0, 255,
        "validation construction attempt")
    if attempt != expected_attempt:
        fail("validation performed runtime search or substituted the map")
    base_precode = _require_hex_seed(
        value.get("base_precode_seed"), 16, "validation base precode seed")
    base_packet = _require_hex_seed(
        value.get("base_packet_seed"), 8, "validation base packet seed")
    if (value.get("candidate_profile_sha256") !=
            candidate_profile_sha256(contract) or
        value.get("worker_binary_sha256") != worker_binary_sha256 or
        value.get("effective_precode_seed") !=
                _effective_precode_seed(base_precode, attempt) or
        value.get("effective_packet_seed") !=
                _effective_packet_seed(base_packet, attempt)):
        fail("validation construction identity changed")
    _require_sha(value.get("source_sha256"), "validation source")
    cells = value.get("cells")
    if not isinstance(cells, list) or \
            len(cells) != VALIDATION_CELL_COUNT:
        fail("validation must retain all nine disjoint cells")
    results = []
    ordinals = []
    for cell in cells:
        results.append(_validate_cell(
            cell, K, attempt, validation=True, require_success=False))
        ordinals.append(cell["cell_ordinal"])
    if ordinals != list(range(VALIDATION_CELL_COUNT)):
        fail("validation cells are not in frozen order")
    passed = all(results)
    if type(value.get("all_oh0_success")) is not bool or \
            value["all_oh0_success"] != passed:
        fail("validation aggregate disagrees with its cells")
    return passed


@dataclass(frozen=True)
class VerifiedWorker:
    path: Path
    descriptor: int
    sha256: str
    source_commit: str
    description: Mapping[str, Any]
    validation_roster_sha256: str


@dataclass(frozen=True)
class VerifiedShortScreen:
    """Authenticated v9 bounded-screen PASS and its frozen 30 attempts."""

    binding: Mapping[str, str]
    attempt_map: bytes
    derivation_records: Tuple[Tuple[int, bytes], ...]
    summary: Mapping[str, Any]


@dataclass(frozen=True)
class VerifiedRepairBundle:
    """Self-consistent derivation artifacts plus their shared identity."""

    repair_map: bytes
    receipt: Mapping[str, Any]
    shared_identity_sha256: str
    semantic_replay_records: Tuple[Tuple[int, bytes], ...]


@dataclass(frozen=True)
class VerifiedValidationBundle:
    """Verified validation artifacts bound to one repair-map receipt."""

    freeze: Mapping[str, Any]
    manifest: Mapping[str, Any]
    summary: Mapping[str, Any]
    shared_identity_sha256: str
    semantic_replay_records: Tuple[Tuple[int, bytes], ...]


def _canonical_jsonl_objects(
        data: bytes, expected_count: int, context: str) \
        -> Tuple[Tuple[Mapping[str, Any], bytes], ...]:
    lines = data.splitlines(keepends=True)
    if len(lines) != expected_count:
        fail("{} has the wrong record count".format(context))
    result = []
    for index, raw in enumerate(lines):
        if (len(raw) > MAX_RECORD_BYTES or not raw.endswith(b"\n") or
                b"\r" in raw):
            fail("{} record {} has invalid framing".format(context, index))
        body = raw[:-1]
        value = _parse_json_bytes(
            body, "{} record {}".format(context, index))
        if (not isinstance(value, dict) or
                body != canonical_json(value).encode("ascii")):
            fail("{} record {} is not canonical".format(context, index))
        result.append((value, body))
    return tuple(result)


def _validate_short_screen_binary_receipt(
        value: Any, context: str, expected_sha256: Optional[str] = None) \
        -> Tuple[str, str]:
    if type(value) is not dict or \
            set(value) != SHORT_SCREEN_BINARY_RECEIPT_FIELDS:
        fail("{} is not an exact binary receipt".format(context))
    path = value.get("path")
    if (type(path) is not str or not path or "\x00" in path or
            not Path(path).is_absolute()):
        fail("{} path is not one absolute path".format(context))
    digest = _require_sha(value.get("sha256"), context + " hash")
    _require_int(value.get("size"), 1, (1 << 63) - 1, context + " size")
    if expected_sha256 is not None and digest != expected_sha256:
        fail("{} hash changed".format(context))
    return path, digest


def _validate_short_screen_source_receipt(
        value: Any, worker: VerifiedWorker) -> str:
    if type(value) is not dict or \
            set(value) != SHORT_SCREEN_SOURCE_RECEIPT_FIELDS:
        fail("short-screen source receipt is not an exact receipt")
    receipt_sha256 = _require_sha(
        value.get("source_receipt_sha256"), "short-screen source receipt")
    source_files = value.get("source_files")
    if (value.get("source_git_commit") != worker.source_commit or
            value.get("tracked_and_untracked_linked_sources_clean") is not
                True or
            not _same_frozen_json(
                value.get("clean_status_scope"),
                list(SHORT_SCREEN_CLEAN_STATUS_SCOPE)) or
            type(source_files) is not dict or
            set(source_files) != set(SHORT_SCREEN_SOURCE_FILES) or
            any(type(path) is not str or type(digest) is not str or
                not SHA256.fullmatch(digest)
                for path, digest in source_files.items()) or
            sha256_json(_object_without_self_hash(
                value, "source_receipt_sha256")) != receipt_sha256):
        fail("short-screen source receipt differs from the clean worker source")
    return receipt_sha256


def _validate_short_screen_result(
        record: Mapping[str, Any], expected: Mapping[str, Any],
        attempt_sha256: str, profile_sha256: str, bench_sha256: str,
        source_commit: str, source_receipt_sha256: str,
        derivation: Mapping[str, Any]) -> bool:
    ordinal = expected["ordinal"]
    context = "short-screen result {}".format(ordinal)
    if type(record) is not dict or set(record) != SHORT_SCREEN_RESULT_FIELDS:
        fail("{} fields differ from the frozen v2 schema".format(context))
    expected_identity = {
        field: expected[field] for field in (
            "ordinal", "coordinate_ordinal", "arm", "K", "root_index",
            "loss_root", "schedule", "cell_ordinal", "timing_order",
            "timing_slot", "observation_index")
    }
    actual_identity = {field: record.get(field) for field in expected_identity}
    candidate = expected["arm"] == "candidate_two07_mix2"
    fixed = {
        "schema": SHORT_SCREEN_RESULT_SCHEMA,
        "block_bytes": 2,
        "loss_ppm": 500000,
        "overhead": 0,
        "attempt_selection_stream_sha256": attempt_sha256,
        "attempt_selection_cell_receipts_used_as_promotion_evidence": False,
        "benchmark_loss_trace_hash_recorded": False,
        "candidate_profile_sha256": profile_sha256,
        "binary_dense_rows": 12,
        "gf256_heavy_rows": 12,
        "dense_anchor_layout": "two07" if candidate else "disabled",
        "mix_count": 2 if candidate else 3,
        "command_sha256": expected["command_sha256"],
        "bench_binary_sha256": bench_sha256,
        "bench_source_git_commit": source_commit,
        "source_receipt_sha256": source_receipt_sha256,
    }
    actual_fixed = {field: record.get(field) for field in fixed}
    if (not _same_frozen_json(actual_identity, expected_identity) or
            not _same_frozen_json(actual_fixed, fixed)):
        fail("{} is outside the fresh frozen benchmark roster".format(context))
    weak = record.get("weak")
    recovery = record.get("full_payload_byte_recovery_verified")
    if type(weak) is not bool or type(recovery) is not bool or \
            recovery != (not weak):
        fail("{} has inconsistent recovery evidence".format(context))
    success = _require_int(
        record.get("success"), 0, 1, context + " success count")
    rank_fail = _require_int(
        record.get("rank_fail"), 0, 1, context + " rank-failure count")
    error = _require_int(
        record.get("error"), 0, 1, context + " error count")
    if success + rank_fail + error != 1 or weak != (success != 1):
        fail("{} outcome counts disagree with weakness".format(context))
    _require_int(record.get("staircase"), 1, 4096,
                 context + " staircase")
    _require_int(record.get("source_hits"), 1, 8,
                 context + " source hits")
    for field in (
            "inactivated_columns", "block_xors", "gf256_muladds",
            "extra_dense_basis_capacity_entries",
            "extra_dense_basis_capacity_bytes"):
        _require_int(record.get(field), 0, (1 << 63) - 1,
                     context + " " + field.replace("_", " "))
    for field in (
            "solve_ms", "build_ms", "peel_ms", "project_ms",
            "residual_ms", "backsub_ms"):
        value = record.get(field)
        if (type(value) is not str or len(value) > 64 or
                not SHORT_SCREEN_DECIMAL.fullmatch(value)):
            fail("{} {} is not a bounded canonical decimal".format(
                context, field.replace("_", " ")))
    construction_attempt = _require_int(
        record.get("construction_attempt"), 0, 255,
        context + " construction attempt")
    if candidate:
        if (record.get("attempt_mode") != "exact" or
                construction_attempt != expected["construction_attempt"] or
                record.get("effective_precode_seed") !=
                    derivation["effective_precode_seed"] or
                record.get("effective_packet_seed") !=
                    derivation["effective_packet_seed"]):
            fail("{} did not use the frozen derived attempt".format(context))
    elif record.get("attempt_mode") != "selected":
        fail("{} does not identify production selection".format(context))
    _require_hex_seed(
        record.get("effective_precode_seed"), 16,
        context + " effective precode seed")
    _require_hex_seed(
        record.get("effective_packet_seed"), 8,
        context + " effective packet seed")
    _require_sha(record.get("bench_stdout_sha256"), context + " stdout")
    return weak


def load_short_screen_bundle(
        contract: Mapping[str, Any], directory: Path, worker: VerifiedWorker,
        expected_summary_sha256: str) -> VerifiedShortScreen:
    """Authenticate the prerequisite v9 bounded-screen PASS fail-closed."""
    expected_summary_sha256 = _require_sha(
        expected_summary_sha256, "externally recorded short-screen summary")
    if directory.is_symlink() or not directory.is_dir():
        fail("short-screen bundle must be a real directory")
    input_path = directory / SHORT_SCREEN_INPUT_NAME
    attempt_path = directory / SHORT_SCREEN_ATTEMPT_NAME
    result_path = directory / SHORT_SCREEN_RESULT_NAME
    summary_path = directory / SHORT_SCREEN_SUMMARY_NAME
    input_bytes = _read_regular(
        input_path, 16 * 1024 * 1024, "short-screen input")
    summary_bytes = _read_regular(
        summary_path, 4 * 1024 * 1024, "short-screen summary")
    input_receipt = _parse_json_bytes(input_bytes, "short-screen input")
    summary = _parse_json_bytes(summary_bytes, "short-screen summary")
    if (not isinstance(input_receipt, dict) or
            set(input_receipt) != SHORT_SCREEN_INPUT_FIELDS or
            input_bytes !=
                (canonical_json(input_receipt) + "\n").encode("ascii")):
        fail("short-screen input is not its exact canonical v2 object")
    if (not isinstance(summary, dict) or
            set(summary) != SHORT_SCREEN_SUMMARY_FIELDS or
            summary_bytes != (canonical_json(summary) + "\n").encode("ascii")):
        fail("short-screen summary is not its exact canonical v2 object")

    screen_contract = input_receipt.get("contract")
    try:
        expected_screen_contract = short_screen_controller.contract_description()
    except Exception as exc:
        fail("cannot reconstruct the frozen short-screen contract: {}".format(
            exc))
    screen_attempt_selection = (
        screen_contract.get("attempt_selection")
        if isinstance(screen_contract, dict) else None)
    screen_timing = (
        screen_contract.get("timing_protocol")
        if isinstance(screen_contract, dict) else None)
    if (not isinstance(screen_contract, dict) or
            not _same_frozen_json(
                screen_contract, expected_screen_contract) or
            screen_contract.get("schema") != SHORT_SCREEN_CONTRACT_SCHEMA or
            screen_contract.get("contract_sha256") !=
                SHORT_SCREEN_CONTRACT_SHA256 or
            sha256_json(_object_without_self_hash(
                screen_contract, "contract_sha256")) !=
                SHORT_SCREEN_CONTRACT_SHA256 or
            screen_contract.get("K") != list(SEMANTIC_REPLAY_K) or
            screen_contract.get("selection_roots") !=
                list(SELECTION_ROOTS) or
            screen_contract.get("benchmark_roots") !=
                list(SHORT_SCREEN_ROOTS) or
            screen_contract.get("schedules") != list(SCHEDULES) or
            screen_contract.get("arms") != list(SHORT_SCREEN_ARMS) or
            type(screen_attempt_selection) is not dict or
            screen_attempt_selection.get("commands") !=
                SHORT_SCREEN_ATTEMPT_COUNT or
            screen_attempt_selection.get("cells_per_K") !=
                SELECTION_CELL_COUNT or
            type(screen_timing) is not dict or
            screen_timing.get("orders") !=
                list(SHORT_SCREEN_TIMING_ORDERS) or
            screen_timing.get("observations_per_arm_per_coordinate") !=
                SHORT_SCREEN_OBSERVATIONS_PER_ARM or
            set(SHORT_SCREEN_ROOTS) &
                (set(SELECTION_ROOTS) | set(VALIDATION_ROOTS))):
        fail("short-screen contract is not the frozen v9 prerequisite")
    input_sha256 = _require_sha(
        input_receipt.get("input_sha256"), "short-screen input identity")
    if (input_receipt.get("schema") != SHORT_SCREEN_INPUT_SCHEMA or
            input_receipt.get("contract_sha256") !=
                SHORT_SCREEN_CONTRACT_SHA256 or
            sha256_json(_object_without_self_hash(
                input_receipt, "input_sha256")) != input_sha256):
        fail("short-screen input identity or self-hash is invalid")

    source_receipt = input_receipt.get("source_receipt")
    source_receipt_sha256 = _validate_short_screen_source_receipt(
        source_receipt, worker)
    if input_receipt.get("source_receipt_sha256") != source_receipt_sha256:
        fail("short-screen input does not bind its source receipt")
    worker_binary = input_receipt.get("repair_worker_binary")
    bench_binary = input_receipt.get("bench_binary")
    _validate_short_screen_binary_receipt(
        worker_binary, "short-screen worker binary", worker.sha256)
    bench_path, bench_binary_sha256 = _validate_short_screen_binary_receipt(
        bench_binary, "short-screen bench binary")
    description = input_receipt.get("repair_worker_description")
    expected_worker_commands = "".join(
        "D {} {}\n".format(K - K_MIN, K) for K in SEMANTIC_REPLAY_K) + \
        "Q\n"
    expected_worker_command_sha256 = sha256_bytes(
        expected_worker_commands.encode("ascii"))
    if (not _same_frozen_json(description, worker.description) or
            _require_sha(
                input_receipt.get("repair_worker_description_sha256"),
                "short-screen worker description") !=
                sha256_json(description) or
            _require_sha(
                input_receipt.get("attempt_selection_worker_command_sha256"),
                "short-screen worker command") !=
                expected_worker_command_sha256):
        fail("short-screen worker identity or command receipt changed")

    selected_attempts = input_receipt.get("selected_attempts")
    invocations = input_receipt.get("invocations")
    if (input_receipt.get("attempt_selection_stream_schema") !=
            SHORT_SCREEN_ATTEMPT_STREAM_SCHEMA or
            _require_int(
                input_receipt.get("attempt_selection_record_count"),
                SHORT_SCREEN_ATTEMPT_COUNT, SHORT_SCREEN_ATTEMPT_COUNT,
                "short-screen attempt records") !=
                    SHORT_SCREEN_ATTEMPT_COUNT or
            not isinstance(selected_attempts, list) or
            len(selected_attempts) != SHORT_SCREEN_ATTEMPT_COUNT or
            any(type(value) is not int or not 0 <= value <= 255
                for value in selected_attempts) or
            type(invocations) is not list or
            len(invocations) != SHORT_SCREEN_RECORD_COUNT or
            _require_int(
                input_receipt.get("worker_validation_commands_issued"), 0, 0,
                "short-screen validation command count") != 0 or
            input_receipt.get("worker_final_validation_stream_present") is not
                False):
        fail("short-screen input did not freeze exactly 30 offline attempts")
    try:
        expected_invocations = tuple(
            invocation.identity(Path(bench_path))
            for invocation in short_screen_controller.make_invocations(dict(
                zip(SEMANTIC_REPLAY_K, selected_attempts))))
    except Exception as exc:
        fail("cannot reconstruct the frozen short-screen invocation roster: "
             "{}".format(exc))
    if len(expected_invocations) != SHORT_SCREEN_RECORD_COUNT:
        fail("short-screen controller produced the wrong invocation count")
    for ordinal, (actual, expected) in enumerate(zip(
            invocations, expected_invocations)):
        if (type(actual) is not dict or
                set(actual) != SHORT_SCREEN_INVOCATION_FIELDS or
                not _same_frozen_json(actual, expected)):
            fail("short-screen invocation {} differs from the frozen command "
                 "roster".format(ordinal))

    attempt_bytes = _read_regular(
        attempt_path, SHORT_SCREEN_ATTEMPT_COUNT * MAX_RECORD_BYTES,
        "short-screen attempt stream")
    attempt_sha256 = sha256_bytes(attempt_bytes)
    if _require_sha(
            input_receipt.get("attempt_selection_stream_sha256"),
            "short-screen attempt stream") != attempt_sha256:
        fail("short-screen attempt stream differs from its input receipt")
    parsed_attempts = _canonical_jsonl_objects(
        attempt_bytes, SHORT_SCREEN_ATTEMPT_COUNT,
        "short-screen attempt stream")
    retained_records = []
    derivation_by_K: Dict[int, Mapping[str, Any]] = {}
    for index, ((record, raw), K, expected_attempt) in enumerate(zip(
            parsed_attempts, SEMANTIC_REPLAY_K, selected_attempts)):
        selected = validate_derivation_record(
            record, contract, worker.sha256, expected_K=K)
        if selected != expected_attempt:
            fail("short-screen attempt map differs at index {}".format(index))
        retained_records.append((K, raw))
        derivation_by_K[K] = record
    short_map = bytes(selected_attempts)

    result_bytes = _read_regular(
        result_path, 64 * 1024 * 1024, "short-screen result stream")
    result_sha256 = sha256_bytes(result_bytes)
    parsed_results = _canonical_jsonl_objects(
        result_bytes, SHORT_SCREEN_RECORD_COUNT, "short-screen result stream")
    for (record, unused), expected in zip(
            parsed_results, expected_invocations):
        _validate_short_screen_result(
            record, expected, attempt_sha256,
            candidate_profile_sha256(contract), bench_binary_sha256,
            worker.source_commit, source_receipt_sha256,
            derivation_by_K[expected["K"]])
    try:
        replayed_verdict = short_screen_controller.adjudicate(
            [record for record, unused in parsed_results], attempt_sha256)
    except Exception as exc:
        fail("short-screen result replay failed: {}".format(exc))
    if (type(replayed_verdict) is not dict or
            set(replayed_verdict) != SHORT_SCREEN_VERDICT_FIELDS or
            any(not _same_frozen_json(summary.get(field), value)
                for field, value in replayed_verdict.items())):
        fail("short-screen summary differs from independently replayed "
             "adjudication")

    gates = summary.get("gates")
    if (summary.get("schema") != SHORT_SCREEN_SUMMARY_SCHEMA or
            _require_sha(summary.get("summary_sha256"),
                         "short-screen summary") !=
                expected_summary_sha256 or
            sha256_json(_object_without_self_hash(
                summary, "summary_sha256")) != expected_summary_sha256 or
            summary.get("contract_sha256") != SHORT_SCREEN_CONTRACT_SHA256 or
            summary.get("input_sha256") != input_sha256 or
            summary.get("input_file_sha256") != sha256_bytes(input_bytes) or
            summary.get("attempt_selection_stream_schema") !=
                SHORT_SCREEN_ATTEMPT_STREAM_SCHEMA or
            summary.get("attempt_selection_stream_sha256") != attempt_sha256 or
            _require_int(
                summary.get("attempt_selection_record_count"),
                SHORT_SCREEN_ATTEMPT_COUNT, SHORT_SCREEN_ATTEMPT_COUNT,
                "short-screen summary attempt records") !=
                    SHORT_SCREEN_ATTEMPT_COUNT or
            summary.get("result_stream_sha256") != result_sha256 or
            _require_int(
                summary.get("record_count"), SHORT_SCREEN_RECORD_COUNT,
                SHORT_SCREEN_RECORD_COUNT, "short-screen result records") !=
                    SHORT_SCREEN_RECORD_COUNT or
            _require_int(
                summary.get("invocation_count"), SHORT_SCREEN_RECORD_COUNT,
                SHORT_SCREEN_RECORD_COUNT, "short-screen invocations") !=
                    SHORT_SCREEN_RECORD_COUNT or
            summary.get("bench_binary_sha256") != bench_binary_sha256 or
            summary.get("repair_worker_binary_sha256") != worker.sha256 or
            summary.get("source_git_commit") != worker.source_commit or
            not _same_frozen_json(
                summary.get("source_receipt"), source_receipt) or
            summary.get("source_receipt_sha256") != source_receipt_sha256 or
            summary.get("candidate_profile_sha256") !=
                candidate_profile_sha256(contract) or
            _require_int(
                summary.get("worker_validation_commands_issued"), 0, 0,
                "short-screen summary validation command count") != 0 or
            summary.get("worker_final_validation_stream_present") is not
                False or
            summary.get("architecture_selection_performed") is not False or
            summary.get("offline_attempt_derivation_performed") is not True or
            summary.get("official_scope") != SHORT_SCREEN_OFFICIAL_SCOPE or
            summary.get("all_K_recovery_claimed") is not False or
            summary.get("wirehair1_timing_claimed") is not False or
            summary.get("disposition") != "PASS" or
            type(gates) is not dict or
            set(gates) != SHORT_SCREEN_GATE_FIELDS or
            any(value is not True for value in gates.values()) or
            _require_int(
                summary.get("candidate_weak_observations"), 0, 0,
                "short-screen candidate weak observations") !=
                    replayed_verdict["candidate_weak_observations"] or
            _require_int(
                summary.get("candidate_weak_coordinates"), 0, 0,
                "short-screen candidate weak coordinates") !=
                    replayed_verdict["candidate_weak_coordinates"]):
        fail("short-screen summary is not an authenticated frozen PASS")

    binding = {
        "short_screen_contract_sha256": SHORT_SCREEN_CONTRACT_SHA256,
        "short_screen_input_sha256": input_sha256,
        "short_screen_attempt_stream_sha256": attempt_sha256,
        "short_screen_result_stream_sha256": result_sha256,
        "short_screen_summary_sha256": expected_summary_sha256,
        "short_screen_map_sha256": sha256_bytes(short_map),
    }
    return VerifiedShortScreen(
        binding, short_map, tuple(retained_records), summary)


def _open_worker(path: Path) -> Tuple[Path, int, str]:
    descriptor = -1
    try:
        original = os.lstat(str(path))
        if stat.S_ISLNK(original.st_mode):
            fail("worker path must not be a symlink")
        resolved = path.resolve(strict=True)
        descriptor = os.open(
            str(resolved), os.O_RDONLY | os.O_CLOEXEC |
            getattr(os, "O_NOFOLLOW", 0))
        opened = os.fstat(descriptor)
        current = os.stat(str(resolved), follow_symlinks=False)
        if (not stat.S_ISREG(opened.st_mode) or not opened.st_mode & 0o111 or
                (opened.st_dev, opened.st_ino) !=
                    (current.st_dev, current.st_ino)):
            fail("worker is not one stable executable regular file")
        return resolved, descriptor, _sha256_open_fd(
            descriptor, "mix2 repair worker")
    except BaseException:
        if descriptor >= 0:
            os.close(descriptor)
        raise


def _run_description(path: Path, descriptor: int) -> Mapping[str, Any]:
    execution = "/proc/self/fd/{}".format(descriptor)
    try:
        completed = subprocess.run(
            [execution, "--describe"], stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            pass_fds=(descriptor,), timeout=30.0, check=False,
            start_new_session=True)
    except (OSError, subprocess.TimeoutExpired) as exc:
        fail("cannot describe mix2 repair worker: {}".format(exc))
    if (completed.returncode != 0 or completed.stderr or
            len(completed.stdout) > MAX_DESCRIPTION_BYTES or
            not completed.stdout.endswith(b"\n") or
            completed.stdout.count(b"\n") != 1):
        fail("worker description command failed or emitted extra output")
    value = _parse_json_bytes(
        completed.stdout[:-1], "worker description")
    if (not isinstance(value, dict) or set(value) != DESCRIPTION_FIELDS or
            completed.stdout !=
                (canonical_json(value) + "\n").encode("ascii")):
        fail("worker description is not one canonical closed object")
    return value


def _validate_worker_description_identity(
        description: Mapping[str, Any], actual_sha256: str,
        expected_source_commit: str, contract: Mapping[str, Any]) -> None:
    if (description.get("binary_sha256") != actual_sha256 or
            description.get("source_git_commit") != expected_source_commit or
            description.get("schema") != DESCRIPTION_SCHEMA or
            description.get("worker_schema") != WORKER_SCHEMA or
            description.get("contract_schema") != CONTRACT_SCHEMA or
            description.get("derivation_schema") !=
                DERIVATION_RECORD_SCHEMA or
            description.get("validation_schema") !=
                VALIDATION_RECORD_SCHEMA or
            description.get("protocol") !=
                "D ordinal K | V ordinal K attempt | Q" or
            description.get("validation_roster_schema") !=
                VALIDATION_ROSTER_SCHEMA or
            description.get("validation_roster_sha256") !=
                VALIDATION_ROSTER_SHA256 or
            not _same_frozen_json(
                description.get("candidate_profile"),
                contract["candidate_profile"]) or
            description.get("candidate_profile_sha256") !=
                candidate_profile_sha256(contract)):
        fail("worker binary, source, protocol, profile, or validation roster "
             "identity changed")


def verify_worker(
        path: Path, expected_sha256: str, expected_source_commit: str,
        contract: Mapping[str, Any]) -> VerifiedWorker:
    expected_sha256 = _require_sha(expected_sha256, "expected worker")
    expected_source_commit = _require_commit(
        expected_source_commit, "expected source commit")
    resolved, descriptor, actual_sha256 = _open_worker(path)
    try:
        description = _run_description(resolved, descriptor)
        if actual_sha256 != expected_sha256:
            fail("worker binary differs from its expected identity")
        _validate_worker_description_identity(
            description, actual_sha256, expected_source_commit, contract)
        return VerifiedWorker(
            resolved, descriptor, actual_sha256, expected_source_commit,
            description, description["validation_roster_sha256"])
    except BaseException:
        os.close(descriptor)
        raise


def _short_screen_binding(value: Mapping[str, Any]) -> Mapping[str, str]:
    if not isinstance(value, Mapping):
        fail("short-screen binding is not an object")
    binding = {
        field: _require_sha(value.get(field), field.replace("_", " "))
        for field in SHORT_SCREEN_BINDING_FIELDS
    }
    if binding["short_screen_contract_sha256"] != \
            SHORT_SCREEN_CONTRACT_SHA256:
        fail("short-screen binding does not name the frozen v9 contract")
    return binding


def _manifest_unsigned(
        action: str, contract: Mapping[str, Any], worker: VerifiedWorker,
        short_screen_binding: Mapping[str, Any],
        repair_map_sha256: Optional[str] = None,
        map_receipt_sha256: Optional[str] = None) -> Mapping[str, Any]:
    if worker.validation_roster_sha256 != VALIDATION_ROSTER_SHA256:
        fail("worker validation roster was not preflight-authenticated")
    if action not in ("derive", "validate"):
        fail("manifest action is invalid")
    if action == "derive":
        repair_map_sha256 = "0" * 64
        map_receipt_sha256 = "0" * 64
    else:
        _require_sha(repair_map_sha256, "validation repair map")
        _require_sha(map_receipt_sha256, "validation map receipt")
    return {
        "schema": MANIFEST_SCHEMA,
        "action": action,
        "contract_sha256": contract_sha256(contract),
        "candidate_arm": contract["candidate_profile"]["arm"],
        "candidate_profile_sha256": candidate_profile_sha256(contract),
        "construction_seed_basis": PRODUCTION_SEED_BASIS,
        "seed_schedule_sha256": PRODUCTION_SEED_SCHEDULE_SHA256,
        "controller_sha256": _controller_sha256(),
        "worker_binary_sha256": worker.sha256,
        "source_git_commit": worker.source_commit,
        "K_min": K_MIN,
        "K_max": K_MAX,
        "map_bytes": MAP_BYTES,
        "selection_roots": list(SELECTION_ROOTS),
        "validation_roots": list(VALIDATION_ROOTS),
        "schedules": list(SCHEDULES),
        "selection_cell_count": SELECTION_CELL_COUNT,
        "validation_cell_count": VALIDATION_CELL_COUNT,
        "validation_roster_sha256": worker.validation_roster_sha256,
        "repair_map_sha256": repair_map_sha256,
        "map_receipt_sha256": map_receipt_sha256,
        **_short_screen_binding(short_screen_binding),
    }


def create_manifest(
        action: str, contract: Mapping[str, Any], worker: VerifiedWorker,
        short_screen_binding: Mapping[str, Any],
        repair_map_sha256: Optional[str] = None,
        map_receipt_sha256: Optional[str] = None) -> Mapping[str, Any]:
    value = dict(_manifest_unsigned(
        action, contract, worker, short_screen_binding,
        repair_map_sha256, map_receipt_sha256))
    value["manifest_sha256"] = sha256_json(value)
    return value


def _validate_manifest_identity(
        manifest: Mapping[str, Any], expected_sha256: str) -> None:
    expected_sha256 = _require_sha(expected_sha256, "manifest")
    actual = manifest.get("manifest_sha256")
    if (actual != expected_sha256 or
            sha256_json(_object_without_self_hash(
                manifest, "manifest_sha256")) != actual):
        fail("fresh manifest disagrees with preregistered identity")


class _Registry:
    def __init__(self) -> None:
        self.lock = threading.Lock()
        self.processes: Dict[int, subprocess.Popen] = {}
        self.cancelled = False

    def add(self, process: subprocess.Popen) -> None:
        with self.lock:
            if self.cancelled:
                _kill_process(process)
                fail("campaign was cancelled")
            self.processes[process.pid] = process

    def remove(self, process: subprocess.Popen) -> None:
        with self.lock:
            self.processes.pop(process.pid, None)

    def cancel(self) -> None:
        with self.lock:
            self.cancelled = True
            processes = list(self.processes.values())
        for process in processes:
            _kill_process(process)


def _kill_process(process: subprocess.Popen) -> None:
    if process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except OSError:
        try:
            process.kill()
        except OSError:
            pass
    try:
        process.wait(timeout=2.0)
    except (OSError, subprocess.TimeoutExpired):
        pass


def _read_response_line(
        process: subprocess.Popen, deadline: float) -> bytes:
    if process.stdout is None:
        fail("worker stdout pipe is unavailable")
    selector = selectors.DefaultSelector()
    data = bytearray()
    try:
        selector.register(process.stdout.fileno(), selectors.EVENT_READ)
        while True:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                fail("campaign deadline expired")
            events = selector.select(remaining)
            if not events:
                fail("campaign deadline expired")
            chunk = os.read(process.stdout.fileno(), 65536)
            if not chunk:
                fail("worker exited before a complete response")
            data.extend(chunk)
            if len(data) > MAX_RECORD_BYTES:
                fail("worker response exceeds its byte cap")
            newline = data.find(b"\n")
            if newline >= 0:
                if newline != len(data) - 1:
                    fail("worker emitted multiple responses for one command")
                return bytes(data[:-1])
    finally:
        selector.close()


def _drain_worker_stdout(
        process: subprocess.Popen, deadline: float, context: str) -> None:
    """Drain through EOF after Q and reject every trailing output byte."""
    if process.stdout is None:
        fail("{} stdout pipe is unavailable".format(context))
    selector = selectors.DefaultSelector()
    extra_bytes = 0
    try:
        selector.register(process.stdout.fileno(), selectors.EVENT_READ)
        while True:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                fail("{} deadline expired while draining stdout".format(
                    context))
            events = selector.select(remaining)
            if not events:
                fail("{} deadline expired while draining stdout".format(
                    context))
            chunk = os.read(process.stdout.fileno(), 65536)
            if not chunk:
                break
            extra_bytes += len(chunk)
            if extra_bytes > MAX_RECORD_BYTES:
                fail("{} emitted excessive trailing stdout".format(context))
    finally:
        selector.close()
    if extra_bytes:
        fail("{} emitted trailing stdout after its final response".format(
            context))


def _read_stderr_file(stream: Any, context: str) -> bytes:
    size = os.fstat(stream.fileno()).st_size
    if size > MAX_STDERR_BYTES:
        fail("{} stderr exceeds its byte cap".format(context))
    stream.seek(0)
    data = stream.read(MAX_STDERR_BYTES + 1)
    if len(data) != size:
        fail("{} stderr changed while reading".format(context))
    return data


def _split_ranges(jobs: int) -> List[Tuple[int, int]]:
    jobs = min(jobs, MAP_BYTES)
    quotient, remainder = divmod(MAP_BYTES, jobs)
    result = []
    start = K_MIN
    for index in range(jobs):
        count = quotient + (1 if index < remainder else 0)
        result.append((start, start + count - 1))
        start += count
    if start != K_MAX + 1:
        fail("internal K partition is incomplete")
    return result


@dataclass(frozen=True)
class _ShardResult:
    start_K: int
    end_K: int
    path: Path
    map_bytes: bytes
    weak_K: Tuple[int, ...]


def _run_shard(
        action: str, start_K: int, end_K: int, shard_path: Path,
        attempts: Optional[bytes], contract: Mapping[str, Any],
        worker: VerifiedWorker, deadline: float, registry: _Registry) \
        -> _ShardResult:
    if worker.validation_roster_sha256 != VALIDATION_ROSTER_SHA256:
        fail("worker validation roster changed before shard execution")
    execution = "/proc/self/fd/{}".format(worker.descriptor)
    stderr_file = tempfile.TemporaryFile(mode="w+b")
    try:
        process = subprocess.Popen(
            [execution, "--worker"], stdin=subprocess.PIPE,
            stdout=subprocess.PIPE, stderr=stderr_file,
            pass_fds=(worker.descriptor,), start_new_session=True)
    except OSError as exc:
        stderr_file.close()
        registry.cancel()
        fail("cannot start mix2 repair worker: {}".format(exc))
    registry.add(process)
    selected = bytearray()
    weak = []
    try:
        if process.stdin is None:
            fail("worker stdin pipe is unavailable")
        with shard_path.open("xb") as output:
            for K in range(start_K, end_K + 1):
                ordinal = K - K_MIN
                if action == "derive":
                    command = "D {} {}\n".format(ordinal, K)
                else:
                    if attempts is None or len(attempts) != MAP_BYTES:
                        fail("validation repair map is incomplete")
                    command = "V {} {} {}\n".format(
                        ordinal, K, attempts[ordinal])
                process.stdin.write(command.encode("ascii"))
                process.stdin.flush()
                raw = _read_response_line(process, deadline)
                value = _parse_json_bytes(raw, "worker response")
                if raw != canonical_json(value).encode("ascii"):
                    fail("worker response is not canonical JSON")
                if action == "derive":
                    attempt = validate_derivation_record(
                        value, contract, worker.sha256, expected_K=K)
                    selected.append(attempt)
                else:
                    attempt = attempts[ordinal]
                    if not validate_validation_record(
                            value, contract, worker.sha256, K, attempt):
                        weak.append(K)
                output.write(raw)
                output.write(b"\n")
            output.flush()
            os.fsync(output.fileno())
        process.stdin.write(b"Q\n")
        process.stdin.flush()
        process.stdin.close()
        _drain_worker_stdout(process, deadline, "campaign worker")
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            fail("campaign deadline expired")
        try:
            process.wait(timeout=remaining)
        except subprocess.TimeoutExpired:
            fail("campaign deadline expired")
        stderr = _read_stderr_file(stderr_file, "campaign worker")
        if process.returncode != 0 or stderr:
            fail("worker did not shut down cleanly: {}".format(
                stderr[:MAX_STDERR_BYTES].decode("utf-8", "replace").strip()))
        expected_count = end_K - start_K + 1
        if action == "derive" and len(selected) != expected_count:
            fail("derivation shard omitted map bytes")
        return _ShardResult(
            start_K, end_K, shard_path, bytes(selected), tuple(weak))
    except BaseException:
        registry.cancel()
        _kill_process(process)
        raise
    finally:
        registry.remove(process)
        for stream in (process.stdin, process.stdout, process.stderr):
            if stream is not None and not stream.closed:
                stream.close()
        stderr_file.close()


ArtifactSource = Union[bytes, Path]


def _publish_sources(
        directory: Path, sources: Mapping[str, ArtifactSource]) -> None:
    if not directory.is_dir() or directory.is_symlink():
        fail("output directory must be a real directory")
    token = "{:x}-{:x}".format(os.getpid(), time.monotonic_ns())
    staged: Dict[str, Path] = {}
    published: List[Path] = []
    try:
        for name, source in sources.items():
            target = directory / name
            if target.exists() or target.is_symlink():
                fail("refusing to replace existing artifact {}".format(name))
            stage = directory / (".{}.{}.stage".format(name, token))
            descriptor = os.open(
                str(stage), os.O_WRONLY | os.O_CREAT | os.O_EXCL |
                os.O_CLOEXEC, 0o600)
            # Register immediately so every failure after creation, including
            # source reads, writes, and fsync, removes the private stage.
            staged[name] = stage
            try:
                if isinstance(source, bytes):
                    chunks: Iterable[bytes] = (source,)
                else:
                    def source_chunks() -> Iterator[bytes]:
                        with source.open("rb") as stream:
                            while True:
                                chunk = stream.read(1024 * 1024)
                                if not chunk:
                                    return
                                yield chunk
                    chunks = source_chunks()
                for chunk in chunks:
                    view = memoryview(chunk)
                    while view:
                        written = os.write(descriptor, view)
                        if written <= 0:
                            fail("short artifact write")
                        view = view[written:]
                os.fsync(descriptor)
            finally:
                os.close(descriptor)
        for name, stage in staged.items():
            target = directory / name
            os.link(str(stage), str(target), follow_symlinks=False)
            published.append(target)
        dir_fd = os.open(str(directory), os.O_RDONLY | os.O_DIRECTORY)
        try:
            for stage in staged.values():
                os.unlink(stage)
            os.fsync(dir_fd)
        finally:
            os.close(dir_fd)
    except BaseException:
        for target in reversed(published):
            try:
                os.unlink(target)
            except OSError:
                pass
        for stage in staged.values():
            try:
                os.unlink(stage)
            except OSError:
                pass
        raise


def _create_output_dir(path: Path) -> Path:
    try:
        path.mkdir(parents=True, exist_ok=False)
    except OSError as exc:
        fail("cannot create fresh output directory: {}".format(exc))
    if path.is_symlink() or not path.is_dir():
        fail("output path is not a real directory")
    return path.resolve(strict=True)


def _freeze_object(manifest: Mapping[str, Any]) -> Mapping[str, Any]:
    freeze = {
        "schema": FREEZE_SCHEMA,
        "action": manifest["action"],
        "manifest_sha256": manifest["manifest_sha256"],
        "manifest": manifest,
    }
    freeze["freeze_sha256"] = sha256_json(freeze)
    return freeze


def _merge_shards(
        shards: Sequence[_ShardResult], output: Path) -> Tuple[str, int]:
    digest = hashlib.sha256()
    records = 0
    with output.open("xb") as destination:
        expected_K = K_MIN
        for shard in sorted(shards, key=lambda item: item.start_K):
            if shard.start_K != expected_K:
                fail("campaign shard partition has a gap")
            with shard.path.open("rb") as source:
                while True:
                    chunk = source.read(1024 * 1024)
                    if not chunk:
                        break
                    destination.write(chunk)
                    digest.update(chunk)
                    records += chunk.count(b"\n")
            expected_K = shard.end_K + 1
        if expected_K != K_MAX + 1:
            fail("campaign shard partition does not end at K=64000")
        destination.flush()
        os.fsync(destination.fileno())
    if records != MAP_BYTES:
        fail("campaign audit does not contain exactly 63,999 records")
    return digest.hexdigest(), records


def _verify_short_screen_reproduction(
        audit_path: Path, repair_map: bytes,
        short_screen: VerifiedShortScreen) -> None:
    """Require the full derivation to reproduce the frozen 30-K screen."""
    if (type(repair_map) is not bytes or len(repair_map) != MAP_BYTES or
            len(short_screen.attempt_map) != SHORT_SCREEN_ATTEMPT_COUNT or
            tuple(K for K, unused in short_screen.derivation_records) !=
                SEMANTIC_REPLAY_K):
        fail("short-screen reproduction inputs are incomplete")
    expected_records = dict(short_screen.derivation_records)
    screen_index = 0
    rows = iter(_iter_canonical_jsonl(
        audit_path, "full derivation short-screen reproduction"))
    for K in range(K_MIN, K_MAX + 1):
        try:
            record = next(rows)
        except StopIteration:
            fail("full derivation ended before short-screen reproduction")
        if K not in expected_records:
            continue
        canonical_record = canonical_json(record).encode("ascii")
        if canonical_record != expected_records[K]:
            fail("full derivation differs from short screen at K={}".format(K))
        selected = repair_map[K - K_MIN]
        if selected != short_screen.attempt_map[screen_index] or \
                record.get("selected_attempt") != selected:
            fail("full derivation attempt differs from short screen at K={}".format(
                K))
        screen_index += 1
    try:
        next(rows)
        fail("full derivation contains records beyond K=64000")
    except StopIteration:
        pass
    if screen_index != SHORT_SCREEN_ATTEMPT_COUNT:
        fail("full derivation omitted short-screen coordinates")


def _run_all_k(
        action: str, contract: Mapping[str, Any], worker: VerifiedWorker,
        output_dir: Path, jobs: int, deadline_seconds: float,
        expected_manifest_sha256: str,
        attempts: Optional[bytes] = None,
        repair_identity: Optional[Mapping[str, Any]] = None,
        repair_shared_identity_sha256: Optional[str] = None,
        short_screen: Optional[VerifiedShortScreen] = None) \
        -> Mapping[str, Any]:
    if type(jobs) is not int or not 1 <= jobs <= MAX_JOBS:
        fail("jobs must be in [1,{}]".format(MAX_JOBS))
    if (not math.isfinite(deadline_seconds) or deadline_seconds <= 0 or
            deadline_seconds > MAX_DEADLINE_SECONDS):
        fail("deadline must be finite and within seven days")
    if action == "derive":
        if short_screen is None:
            fail("derivation requires an authenticated short-screen PASS")
        short_binding: Mapping[str, Any] = short_screen.binding
    else:
        if repair_identity is None:
            fail("validation requires an authenticated repair identity")
        short_binding = repair_identity
    manifest = create_manifest(
        action, contract, worker, short_binding,
        None if repair_identity is None else repair_identity["map_sha256"],
        None if repair_identity is None else
            repair_identity["receipt_sha256"])
    _validate_manifest_identity(manifest, expected_manifest_sha256)
    _require_clean_source(worker.source_commit)
    output_dir = _create_output_dir(output_dir)
    freeze_name = DERIVATION_FREEZE_NAME if action == "derive" else \
        VALIDATION_FREEZE_NAME
    freeze = _freeze_object(manifest)
    _publish_sources(output_dir, {
        freeze_name: (canonical_json(freeze) + "\n").encode("ascii")})

    deadline = time.monotonic() + deadline_seconds
    registry = _Registry()
    with tempfile.TemporaryDirectory(prefix="wh2-mix2-v9-") as temporary:
        temporary_path = Path(temporary)
        ranges = _split_ranges(jobs)
        futures = []
        shards: List[_ShardResult] = []
        with ThreadPoolExecutor(max_workers=len(ranges)) as executor:
            for index, (start_K, end_K) in enumerate(ranges):
                futures.append(executor.submit(
                    _run_shard, action, start_K, end_K,
                    temporary_path / "shard-{:03d}.jsonl".format(index),
                    attempts, contract, worker, deadline, registry))
            try:
                for future in as_completed(futures):
                    shards.append(future.result())
            except BaseException:
                registry.cancel()
                for future in futures:
                    future.cancel()
                raise
        if len(shards) != len(ranges):
            fail("campaign did not complete every shard")
        if (_controller_sha256() != manifest["controller_sha256"] or
                _sha256_open_fd(worker.descriptor, "mix2 repair worker") !=
                    worker.sha256):
            fail("controller or worker changed during the all-K campaign")
        _require_clean_source(worker.source_commit)
        merged = temporary_path / "merged.jsonl"
        audit_sha256, record_count = _merge_shards(shards, merged)

        if action == "validate":
            expected_shared_identity = _require_sha(
                repair_shared_identity_sha256,
                "repair shared-coordinate identity")
            if _audit_shared_identity_sha256(merged) != \
                    expected_shared_identity:
                fail("derivation and validation source/base identities differ")

        if action == "derive":
            ordered = sorted(shards, key=lambda item: item.start_K)
            repair_map = b"".join(item.map_bytes for item in ordered)
            if len(repair_map) != MAP_BYTES:
                fail("derived repair map is not exactly 63,999 bytes")
            _verify_short_screen_reproduction(
                merged, repair_map, short_screen)
            map_sha256 = sha256_bytes(repair_map)
            histogram = Counter(repair_map)
            unsigned_summary = {
                "schema": DERIVATION_SUMMARY_SCHEMA,
                "contract_sha256": contract_sha256(contract),
                "candidate_profile_sha256":
                    candidate_profile_sha256(contract),
                "manifest_sha256": manifest["manifest_sha256"],
                "audit_sha256": audit_sha256,
                "audit_records": record_count,
                "map_sha256": map_sha256,
                "map_bytes": len(repair_map),
                "maximum_selected_attempt": max(repair_map),
                "selected_attempt_histogram": {
                    str(key): histogram[key] for key in sorted(histogram)},
                "runtime_search": False,
            }
            summary = dict(unsigned_summary)
            summary["summary_sha256"] = sha256_json(unsigned_summary)
            unsigned_receipt = {
                "schema": MAP_SCHEMA,
                "contract_sha256": contract_sha256(contract),
                "candidate_arm": contract["candidate_profile"]["arm"],
                "candidate_profile_sha256":
                    candidate_profile_sha256(contract),
                "construction_seed_basis": PRODUCTION_SEED_BASIS,
                "seed_schedule_sha256": PRODUCTION_SEED_SCHEDULE_SHA256,
                "entry_kind": "uint8_attempt_indexed_by_K_minus_2",
                "map_bytes": MAP_BYTES,
                "map_sha256": map_sha256,
                "derivation_audit_records": record_count,
                "derivation_audit_sha256": audit_sha256,
                "worker_binary_sha256": worker.sha256,
                "source_git_commit": worker.source_commit,
                "controller_sha256": _controller_sha256(),
                "derivation_manifest_sha256": manifest["manifest_sha256"],
                "derivation_freeze_sha256": freeze["freeze_sha256"],
                "derivation_summary_sha256": summary["summary_sha256"],
                "selection_cells_per_selected_attempt":
                    SELECTION_CELL_COUNT,
                "selection_rule":
                    contract["selection"]["selection_rule"],
                **_short_screen_binding(short_binding),
            }
            receipt = dict(unsigned_receipt)
            receipt["receipt_sha256"] = sha256_json(unsigned_receipt)
            _publish_sources(output_dir, {
                DERIVATION_AUDIT_NAME: merged,
                MAP_NAME: repair_map,
                MAP_RECEIPT_NAME:
                    (canonical_json(receipt) + "\n").encode("ascii"),
                DERIVATION_SUMMARY_NAME:
                    (canonical_json(summary) + "\n").encode("ascii"),
            })
            return summary

        weak_K = sorted(K for shard in shards for K in shard.weak_K)
        unsigned_summary = {
            "schema": VALIDATION_SUMMARY_SCHEMA,
            "contract_sha256": contract_sha256(contract),
            "candidate_profile_sha256": candidate_profile_sha256(contract),
            "manifest_sha256": manifest["manifest_sha256"],
            "repair_map_sha256": repair_identity["map_sha256"],
            "map_receipt_sha256": repair_identity["receipt_sha256"],
            "audit_sha256": audit_sha256,
            "audit_records": record_count,
            "weak_K_count": len(weak_K),
            "weak_K": weak_K,
            "runtime_search": False,
            "disposition": "PASS" if not weak_K else "REJECT",
        }
        summary = dict(unsigned_summary)
        summary["summary_sha256"] = sha256_json(unsigned_summary)
        _publish_sources(output_dir, {
            VALIDATION_AUDIT_NAME: merged,
            VALIDATION_SUMMARY_NAME:
                (canonical_json(summary) + "\n").encode("ascii"),
        })
        return summary


def _iter_canonical_jsonl(path: Path, context: str) \
        -> Iterator[Mapping[str, Any]]:
    descriptor = -1
    try:
        nofollow = getattr(os, "O_NOFOLLOW", 0)
        if nofollow == 0:
            fail("{} cannot be opened fail-closed".format(context))
        descriptor = os.open(
            str(path), os.O_RDONLY | os.O_CLOEXEC | nofollow)
        before = os.fstat(descriptor)
        if (not stat.S_ISREG(before.st_mode) or before.st_size < 0 or
                before.st_size > MAX_AUDIT_BYTES):
            fail("{} has an invalid type or size".format(context))
        total = 0
        with os.fdopen(os.dup(descriptor), "rb") as stream:
            for line_number, raw in enumerate(stream, 1):
                total += len(raw)
                if total > before.st_size:
                    fail("{} grew while reading".format(context))
                if len(raw) > MAX_RECORD_BYTES:
                    fail("{} line {} exceeds its byte cap".format(
                        context, line_number))
                if not raw.endswith(b"\n") or b"\r" in raw:
                    fail("{} line {} has invalid framing".format(
                        context, line_number))
                value = _parse_json_bytes(
                    raw[:-1], "{} line {}".format(context, line_number))
                if (not isinstance(value, dict) or
                        raw != (canonical_json(value) + "\n").encode("ascii")):
                    fail("{} line {} is not canonical".format(
                        context, line_number))
                yield value
        after = os.fstat(descriptor)
        identity = lambda value: (
            value.st_dev, value.st_ino, value.st_size,
            value.st_mtime_ns, value.st_ctime_ns)
        if total != before.st_size or identity(before) != identity(after):
            fail("{} changed while reading".format(context))
    except OSError as exc:
        fail("cannot read {}: {}".format(context, exc))
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _audit_shared_identity_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    rows = iter(_iter_canonical_jsonl(path, "shared-coordinate audit"))
    for K in range(K_MIN, K_MAX + 1):
        try:
            record = next(rows)
        except StopIteration:
            fail("shared-coordinate audit ended before K=64000")
        if record.get("K") != K:
            fail("shared-coordinate audit is not ordered by K")
        digest.update(_shared_coordinate_identity(record))
    try:
        next(rows)
        fail("shared-coordinate audit contains records beyond K=64000")
    except StopIteration:
        pass
    return digest.hexdigest()


def load_repair_bundle(
        contract: Mapping[str, Any], directory: Path) \
        -> VerifiedRepairBundle:
    if directory.is_symlink() or not directory.is_dir():
        fail("repair bundle must be a real directory")
    receipt = _load_canonical_object(
        directory / MAP_RECEIPT_NAME, 1024 * 1024, "repair-map receipt")
    _require_int(
        receipt.get("map_bytes"), MAP_BYTES, MAP_BYTES,
        "repair-map receipt bytes")
    _require_int(
        receipt.get("derivation_audit_records"), MAP_BYTES, MAP_BYTES,
        "repair-map derivation records")
    _require_int(
        receipt.get("selection_cells_per_selected_attempt"),
        SELECTION_CELL_COUNT, SELECTION_CELL_COUNT,
        "repair-map selection cell count")
    if (set(receipt) != MAP_RECEIPT_FIELDS or
            receipt.get("schema") != MAP_SCHEMA or
            receipt.get("contract_sha256") != contract_sha256(contract) or
            receipt.get("candidate_arm") !=
                contract["candidate_profile"]["arm"] or
            receipt.get("candidate_profile_sha256") !=
                candidate_profile_sha256(contract) or
            receipt.get("construction_seed_basis") !=
                PRODUCTION_SEED_BASIS or
            receipt.get("seed_schedule_sha256") !=
                PRODUCTION_SEED_SCHEDULE_SHA256 or
            receipt.get("entry_kind") !=
                "uint8_attempt_indexed_by_K_minus_2" or
            receipt.get("selection_rule") !=
                contract["selection"]["selection_rule"] or
            receipt.get("controller_sha256") != _controller_sha256() or
            sha256_json(_object_without_self_hash(
                receipt, "receipt_sha256")) !=
                receipt.get("receipt_sha256")):
        fail("repair-map receipt identity or self-hash is invalid")
    _require_sha(receipt.get("worker_binary_sha256"), "map worker")
    _require_commit(receipt.get("source_git_commit"), "map source commit")
    _require_sha(receipt.get("map_sha256"), "map bytes")
    _require_sha(receipt.get("derivation_audit_sha256"), "derivation audit")
    derivation_manifest_sha256 = _require_sha(
        receipt.get("derivation_manifest_sha256"), "derivation manifest")
    derivation_freeze_sha256 = _require_sha(
        receipt.get("derivation_freeze_sha256"), "derivation freeze")
    derivation_summary_sha256 = _require_sha(
        receipt.get("derivation_summary_sha256"), "derivation summary")

    freeze = _load_canonical_object(
        directory / DERIVATION_FREEZE_NAME,
        1024 * 1024, "derivation freeze")
    if (set(freeze) != FREEZE_FIELDS or
            freeze.get("schema") != FREEZE_SCHEMA or
            freeze.get("action") != "derive" or
            freeze.get("manifest_sha256") !=
                derivation_manifest_sha256 or
            freeze.get("freeze_sha256") != derivation_freeze_sha256 or
            sha256_json(_object_without_self_hash(
                freeze, "freeze_sha256")) != freeze.get("freeze_sha256")):
        fail("derivation freeze identity or self-hash is invalid")
    manifest = freeze.get("manifest")
    if not isinstance(manifest, dict) or set(manifest) != MANIFEST_FIELDS:
        fail("derivation manifest fields changed")
    _require_int(manifest.get("K_min"), K_MIN, K_MIN, "derivation K_min")
    _require_int(manifest.get("K_max"), K_MAX, K_MAX, "derivation K_max")
    _require_int(
        manifest.get("map_bytes"), MAP_BYTES, MAP_BYTES,
        "derivation map bytes")
    if (
            manifest.get("schema") != MANIFEST_SCHEMA or
            manifest.get("action") != "derive" or
            manifest.get("manifest_sha256") !=
                derivation_manifest_sha256 or
            sha256_json(_object_without_self_hash(
                manifest, "manifest_sha256")) !=
                derivation_manifest_sha256 or
            manifest.get("contract_sha256") != contract_sha256(contract) or
            manifest.get("candidate_arm") !=
                contract["candidate_profile"]["arm"] or
            manifest.get("candidate_profile_sha256") !=
                candidate_profile_sha256(contract) or
            manifest.get("construction_seed_basis") !=
                PRODUCTION_SEED_BASIS or
            manifest.get("seed_schedule_sha256") !=
                PRODUCTION_SEED_SCHEDULE_SHA256 or
            manifest.get("controller_sha256") !=
                receipt["controller_sha256"] or
            manifest.get("worker_binary_sha256") !=
                receipt["worker_binary_sha256"] or
            manifest.get("source_git_commit") !=
                receipt["source_git_commit"] or
            manifest.get("selection_roots") != list(SELECTION_ROOTS) or
            manifest.get("validation_roots") != list(VALIDATION_ROOTS) or
            manifest.get("schedules") != list(SCHEDULES) or
            _require_int(
                manifest.get("selection_cell_count"),
                SELECTION_CELL_COUNT, SELECTION_CELL_COUNT,
                "derivation selection cell count") != SELECTION_CELL_COUNT or
            _require_int(
                manifest.get("validation_cell_count"),
                VALIDATION_CELL_COUNT, VALIDATION_CELL_COUNT,
                "derivation validation cell count") != VALIDATION_CELL_COUNT or
            manifest.get("validation_roster_sha256") !=
                VALIDATION_ROSTER_SHA256 or
            _short_screen_binding(manifest) !=
                _short_screen_binding(receipt) or
            manifest.get("repair_map_sha256") != "0" * 64 or
            manifest.get("map_receipt_sha256") != "0" * 64):
        fail("derivation manifest differs from its frozen campaign")

    summary = _load_canonical_object(
        directory / DERIVATION_SUMMARY_NAME,
        4 * 1024 * 1024, "derivation summary")
    _require_int(
        summary.get("audit_records"), MAP_BYTES, MAP_BYTES,
        "derivation summary records")
    _require_int(
        summary.get("map_bytes"), MAP_BYTES, MAP_BYTES,
        "derivation summary map bytes")
    maximum_selected_attempt = _require_int(
        summary.get("maximum_selected_attempt"), 0, 255,
        "maximum selected attempt")
    if (set(summary) != DERIVATION_SUMMARY_FIELDS or
            summary.get("schema") != DERIVATION_SUMMARY_SCHEMA or
            summary.get("contract_sha256") != contract_sha256(contract) or
            summary.get("candidate_profile_sha256") !=
                candidate_profile_sha256(contract) or
            summary.get("manifest_sha256") !=
                derivation_manifest_sha256 or
            summary.get("audit_sha256") !=
                receipt["derivation_audit_sha256"] or
            summary.get("map_sha256") != receipt["map_sha256"] or
            summary.get("runtime_search") is not False or
            summary.get("summary_sha256") != derivation_summary_sha256 or
            sha256_json(_object_without_self_hash(
                summary, "summary_sha256")) != summary.get("summary_sha256")):
        fail("derivation summary identity or self-hash is invalid")
    repair_map = _read_regular(
        directory / MAP_NAME, MAP_BYTES, "repair-map bytes")
    if len(repair_map) != MAP_BYTES or \
            sha256_bytes(repair_map) != receipt["map_sha256"]:
        fail("repair map is not its authenticated 63,999-byte stream")
    histogram = Counter(repair_map)
    expected_histogram = {
        str(key): histogram[key] for key in sorted(histogram)}
    received_histogram = summary.get("selected_attempt_histogram")
    if (not isinstance(received_histogram, dict) or
            set(received_histogram) != set(expected_histogram) or
            any(type(value) is not int or value < 1
                for value in received_histogram.values()) or
            received_histogram != expected_histogram or
            maximum_selected_attempt != max(repair_map)):
        fail("derivation summary disagrees with the repair map")
    digest = hashlib.sha256()
    shared_identity_digest = hashlib.sha256()
    semantic_replay_records: List[Tuple[int, bytes]] = []
    semantic_replay_set = frozenset(SEMANTIC_REPLAY_K)
    count = 0
    rows = iter(_iter_canonical_jsonl(
        directory / DERIVATION_AUDIT_NAME, "derivation audit"))
    for K in range(K_MIN, K_MAX + 1):
        try:
            record = next(rows)
        except StopIteration:
            fail("derivation audit ended before K=64000")
        canonical_record = canonical_json(record).encode("ascii")
        digest.update(canonical_record + b"\n")
        selected = validate_derivation_record(
            record, contract, receipt["worker_binary_sha256"], expected_K=K)
        if repair_map[K - K_MIN] != selected:
            fail("repair map differs from its derivation audit")
        shared_identity_digest.update(_shared_coordinate_identity(record))
        if K in semantic_replay_set:
            semantic_replay_records.append((K, canonical_record))
        count += 1
    try:
        next(rows)
        fail("derivation audit contains records beyond K=64000")
    except StopIteration:
        pass
    if (count != MAP_BYTES or digest.hexdigest() !=
            receipt["derivation_audit_sha256"]):
        fail("derivation audit cardinality or SHA-256 is invalid")
    if tuple(K for K, unused in semantic_replay_records) != \
            SEMANTIC_REPLAY_K:
        fail("derivation audit omitted the semantic replay roster")
    return VerifiedRepairBundle(
        repair_map, receipt, shared_identity_digest.hexdigest(),
        tuple(semantic_replay_records))


def _require_derivation_manifest_anchor(
        repair_bundle: VerifiedRepairBundle,
        expected_manifest_sha256: str) -> str:
    """Authenticate a repair bundle with its pre-derivation external anchor."""
    if not isinstance(repair_bundle, VerifiedRepairBundle):
        fail("derivation manifest anchor requires one verified repair bundle")
    expected = _require_sha(
        expected_manifest_sha256, "externally recorded derivation manifest")
    receipt = repair_bundle.receipt
    if (not isinstance(receipt, Mapping) or
            _require_sha(
                receipt.get("derivation_manifest_sha256"),
                "repair-bundle derivation manifest") != expected):
        fail("repair bundle differs from its preregistered derivation manifest")
    return expected


def load_validation_bundle(
        contract: Mapping[str, Any], directory: Path,
        repair_bundle: VerifiedRepairBundle, worker_binary_sha256: str,
        source_git_commit: str) -> VerifiedValidationBundle:
    """Verify a saved PASS bundle against a verified repair bundle.

    The repair bundle's essential identity and self-hash are rechecked here so
    a caller cannot accidentally cross-wire two bundle chains.
    """
    if not isinstance(repair_bundle, VerifiedRepairBundle):
        fail("validation requires one verified repair bundle")
    repair_map = repair_bundle.repair_map
    repair_receipt = repair_bundle.receipt
    repair_shared_identity_sha256 = _require_sha(
        repair_bundle.shared_identity_sha256,
        "repair shared-coordinate identity")
    worker_binary_sha256 = _require_sha(
        worker_binary_sha256, "validation worker")
    source_git_commit = _require_commit(
        source_git_commit, "validation source commit")
    if directory.is_symlink() or not directory.is_dir():
        fail("validation bundle must be a real directory")
    if (type(repair_map) is not bytes or len(repair_map) != MAP_BYTES or
            not isinstance(repair_receipt, dict) or
            set(repair_receipt) != MAP_RECEIPT_FIELDS or
            repair_receipt.get("schema") != MAP_SCHEMA or
            repair_receipt.get("contract_sha256") !=
                contract_sha256(contract) or
            repair_receipt.get("candidate_arm") !=
                contract["candidate_profile"]["arm"] or
            repair_receipt.get("candidate_profile_sha256") !=
                candidate_profile_sha256(contract) or
            repair_receipt.get("construction_seed_basis") !=
                PRODUCTION_SEED_BASIS or
            repair_receipt.get("seed_schedule_sha256") !=
                PRODUCTION_SEED_SCHEDULE_SHA256 or
            repair_receipt.get("map_bytes") != MAP_BYTES or
            sha256_bytes(repair_map) !=
                repair_receipt.get("map_sha256") or
            repair_receipt.get("worker_binary_sha256") !=
                worker_binary_sha256 or
            repair_receipt.get("source_git_commit") != source_git_commit or
            sha256_json(_object_without_self_hash(
                repair_receipt, "receipt_sha256")) !=
                repair_receipt.get("receipt_sha256")):
        fail("validation repair bundle identity is invalid")

    freeze = _load_canonical_object(
        directory / VALIDATION_FREEZE_NAME,
        1024 * 1024, "validation freeze")
    if (set(freeze) != FREEZE_FIELDS or
            freeze.get("schema") != FREEZE_SCHEMA or
            freeze.get("action") != "validate" or
            _require_sha(
                freeze.get("manifest_sha256"),
                "validation manifest") !=
                freeze.get("manifest_sha256") or
            _require_sha(
                freeze.get("freeze_sha256"),
                "validation freeze") != freeze.get("freeze_sha256") or
            sha256_json(_object_without_self_hash(
                freeze, "freeze_sha256")) != freeze.get("freeze_sha256")):
        fail("validation freeze identity or self-hash is invalid")
    manifest = freeze.get("manifest")
    if not isinstance(manifest, dict) or set(manifest) != MANIFEST_FIELDS:
        fail("validation manifest fields changed")
    if (_require_int(
            manifest.get("K_min"), K_MIN, K_MIN,
            "validation K_min") != K_MIN or
            _require_int(
                manifest.get("K_max"), K_MAX, K_MAX,
                "validation K_max") != K_MAX or
            _require_int(
                manifest.get("map_bytes"), MAP_BYTES, MAP_BYTES,
                "validation map bytes") != MAP_BYTES):
        fail("validation manifest domain changed")
    manifest_sha256 = freeze["manifest_sha256"]
    if (manifest.get("schema") != MANIFEST_SCHEMA or
            manifest.get("action") != "validate" or
            manifest.get("manifest_sha256") != manifest_sha256 or
            sha256_json(_object_without_self_hash(
                manifest, "manifest_sha256")) != manifest_sha256 or
            manifest.get("contract_sha256") != contract_sha256(contract) or
            manifest.get("candidate_arm") !=
                contract["candidate_profile"]["arm"] or
            manifest.get("candidate_profile_sha256") !=
                candidate_profile_sha256(contract) or
            manifest.get("construction_seed_basis") !=
                PRODUCTION_SEED_BASIS or
            manifest.get("seed_schedule_sha256") !=
                PRODUCTION_SEED_SCHEDULE_SHA256 or
            manifest.get("controller_sha256") != _controller_sha256() or
            manifest.get("controller_sha256") !=
                repair_receipt.get("controller_sha256") or
            manifest.get("worker_binary_sha256") != worker_binary_sha256 or
            manifest.get("source_git_commit") != source_git_commit or
            manifest.get("selection_roots") != list(SELECTION_ROOTS) or
            manifest.get("validation_roots") != list(VALIDATION_ROOTS) or
            manifest.get("schedules") != list(SCHEDULES) or
            _require_int(
                manifest.get("selection_cell_count"),
                SELECTION_CELL_COUNT, SELECTION_CELL_COUNT,
                "validation selection cell count") != SELECTION_CELL_COUNT or
            _require_int(
                manifest.get("validation_cell_count"),
                VALIDATION_CELL_COUNT, VALIDATION_CELL_COUNT,
                "validation cell count") != VALIDATION_CELL_COUNT or
            manifest.get("validation_roster_sha256") !=
                VALIDATION_ROSTER_SHA256 or
            _short_screen_binding(manifest) !=
                _short_screen_binding(repair_receipt) or
            manifest.get("repair_map_sha256") !=
                repair_receipt.get("map_sha256") or
            manifest.get("map_receipt_sha256") !=
                repair_receipt.get("receipt_sha256")):
        fail("validation manifest differs from its frozen campaign")

    summary = _load_canonical_object(
        directory / VALIDATION_SUMMARY_NAME,
        4 * 1024 * 1024, "validation summary")
    if not isinstance(summary, dict) or \
            set(summary) != VALIDATION_SUMMARY_FIELDS:
        fail("validation summary fields changed")
    audit_sha256 = _require_sha(
        summary.get("audit_sha256"), "validation audit")
    if (_require_int(
            summary.get("audit_records"), MAP_BYTES, MAP_BYTES,
            "validation audit records") != MAP_BYTES or
            _require_int(
                summary.get("weak_K_count"), 0, 0,
                "validation weak K count") != 0):
        fail("validation summary cardinality changed")
    if (summary.get("schema") != VALIDATION_SUMMARY_SCHEMA or
            summary.get("contract_sha256") != contract_sha256(contract) or
            summary.get("candidate_profile_sha256") !=
                candidate_profile_sha256(contract) or
            summary.get("manifest_sha256") != manifest_sha256 or
            summary.get("repair_map_sha256") !=
                repair_receipt.get("map_sha256") or
            summary.get("map_receipt_sha256") !=
                repair_receipt.get("receipt_sha256") or
            summary.get("weak_K") != [] or
            summary.get("runtime_search") is not False or
            summary.get("disposition") != "PASS" or
            _require_sha(
                summary.get("summary_sha256"),
                "validation summary") != summary.get("summary_sha256") or
            sha256_json(_object_without_self_hash(
                summary, "summary_sha256")) !=
                summary.get("summary_sha256")):
        fail("validation summary is not an authenticated zero-weak-K PASS")

    _semantic_replay_roster_sha256()
    digest = hashlib.sha256()
    shared_identity_digest = hashlib.sha256()
    semantic_replay_records: List[Tuple[int, bytes]] = []
    semantic_replay_set = frozenset(SEMANTIC_REPLAY_K)
    count = 0
    rows = iter(_iter_canonical_jsonl(
        directory / VALIDATION_AUDIT_NAME, "validation audit"))
    for K in range(K_MIN, K_MAX + 1):
        try:
            record = next(rows)
        except StopIteration:
            fail("validation audit ended before K=64000")
        canonical_record = canonical_json(record).encode("ascii")
        digest.update(canonical_record + b"\n")
        attempt = repair_map[K - K_MIN]
        if not validate_validation_record(
                record, contract, worker_binary_sha256, K, attempt):
            fail("validation audit contains an OH0-weak K")
        shared_identity_digest.update(_shared_coordinate_identity(record))
        if K in semantic_replay_set:
            semantic_replay_records.append((K, canonical_record))
        count += 1
    try:
        next(rows)
        fail("validation audit contains records beyond K=64000")
    except StopIteration:
        pass
    if count != MAP_BYTES or digest.hexdigest() != audit_sha256:
        fail("validation audit cardinality or SHA-256 is invalid")
    shared_identity_sha256 = shared_identity_digest.hexdigest()
    if shared_identity_sha256 != repair_shared_identity_sha256:
        fail("derivation and validation source/base identities differ")
    if tuple(K for K, unused in semantic_replay_records) != \
            SEMANTIC_REPLAY_K:
        fail("validation audit omitted the semantic replay roster")
    return VerifiedValidationBundle(
        freeze, manifest, summary, shared_identity_sha256,
        tuple(semantic_replay_records))


def _render_repair_map_cpp_include(
        contract: Mapping[str, Any], repair_bundle: VerifiedRepairBundle,
        validation: VerifiedValidationBundle) -> bytes:
    """Render initializer bytes plus deterministic authenticated provenance."""
    repair_map = repair_bundle.repair_map
    repair_receipt = repair_bundle.receipt
    summary = validation.summary
    lines = [
        "// schema={}".format(MAP_EXPORT_SCHEMA),
        "// contract_sha256={}".format(contract_sha256(contract)),
        "// candidate_profile_sha256={}".format(
            candidate_profile_sha256(contract)),
        "// repair_map_sha256={}".format(repair_receipt["map_sha256"]),
        "// repair_receipt_sha256={}".format(
            repair_receipt["receipt_sha256"]),
        "// derivation_manifest_sha256={}".format(
            repair_receipt["derivation_manifest_sha256"]),
        "// derivation_freeze_sha256={}".format(
            repair_receipt["derivation_freeze_sha256"]),
        "// derivation_audit_sha256={}".format(
            repair_receipt["derivation_audit_sha256"]),
        "// derivation_summary_sha256={}".format(
            repair_receipt["derivation_summary_sha256"]),
        "// short_screen_contract_sha256={}".format(
            repair_receipt["short_screen_contract_sha256"]),
        "// short_screen_input_sha256={}".format(
            repair_receipt["short_screen_input_sha256"]),
        "// short_screen_attempt_stream_sha256={}".format(
            repair_receipt["short_screen_attempt_stream_sha256"]),
        "// short_screen_result_stream_sha256={}".format(
            repair_receipt["short_screen_result_stream_sha256"]),
        "// short_screen_summary_sha256={}".format(
            repair_receipt["short_screen_summary_sha256"]),
        "// short_screen_map_sha256={}".format(
            repair_receipt["short_screen_map_sha256"]),
        "// validation_manifest_sha256={}".format(
            validation.manifest["manifest_sha256"]),
        "// validation_freeze_sha256={}".format(
            validation.freeze["freeze_sha256"]),
        "// validation_audit_sha256={}".format(summary["audit_sha256"]),
        "// validation_summary_sha256={}".format(
            summary["summary_sha256"]),
        "// validation_roster_sha256={}".format(
            VALIDATION_ROSTER_SHA256),
        "// shared_coordinate_identity_sha256={}".format(
            validation.shared_identity_sha256),
        "// semantic_derivation_validation_replay_roster_sha256={}".format(
            _semantic_replay_roster_sha256()),
        "// semantic_derivation_replay_records={}".format(
            len(repair_bundle.semantic_replay_records)),
        "// semantic_validation_replay_records={}".format(
            len(validation.semantic_replay_records)),
        "// validation_disposition={}".format(summary["disposition"]),
        "// runtime_search=false",
        "// worker_binary_sha256={}".format(
            repair_receipt["worker_binary_sha256"]),
        "// source_git_commit={}".format(
            repair_receipt["source_git_commit"]),
        "// map_index=K-{}; K_min={}; K_max={}; map_bytes={}".format(
            K_MIN, K_MIN, K_MAX, MAP_BYTES),
        "// clang-format off",
    ]
    for offset in range(0, len(repair_map), 16):
        lines.append("    " + ", ".join(
            "0x{:02x}".format(value)
            for value in repair_map[offset:offset + 16]) + ",")
    lines.append("// clang-format on")
    return ("\n".join(lines) + "\n").encode("ascii")


def _replay_semantic_sample(
        worker: VerifiedWorker, repair_map: bytes,
        derivation_records: Sequence[Tuple[int, bytes]],
        validation_records: Sequence[Tuple[int, bytes]]) -> str:
    """Re-run frozen 30-K derivation and validation records exactly."""
    if worker.validation_roster_sha256 != VALIDATION_ROSTER_SHA256:
        fail("worker validation roster changed before semantic replay")
    roster_sha256 = _semantic_replay_roster_sha256()
    if (type(repair_map) is not bytes or len(repair_map) != MAP_BYTES or
            tuple(K for K, unused in derivation_records) !=
                SEMANTIC_REPLAY_K or
            tuple(K for K, unused in validation_records) !=
                SEMANTIC_REPLAY_K):
        fail("semantic replay inputs differ from the frozen roster")
    expected_derivation = dict(derivation_records)
    expected_validation = dict(validation_records)
    if (len(expected_derivation) != len(SEMANTIC_REPLAY_K) or
            len(expected_validation) != len(SEMANTIC_REPLAY_K) or
            any(type(raw) is not bytes or not raw
                for raw in expected_derivation.values()) or
            any(type(raw) is not bytes or not raw
                for raw in expected_validation.values())):
        fail("semantic replay records are incomplete")

    execution = "/proc/self/fd/{}".format(worker.descriptor)
    process: Optional[subprocess.Popen] = None
    stderr_file = tempfile.TemporaryFile(mode="w+b")
    deadline = time.monotonic() + SEMANTIC_REPLAY_DEADLINE_SECONDS
    try:
        process = subprocess.Popen(
            [execution, "--worker"], stdin=subprocess.PIPE,
            stdout=subprocess.PIPE, stderr=stderr_file,
            pass_fds=(worker.descriptor,), start_new_session=True)
        if process.stdin is None:
            fail("semantic replay worker stdin is unavailable")
        for K in SEMANTIC_REPLAY_K:
            attempt = repair_map[K - K_MIN]
            process.stdin.write(
                "D {} {}\n".format(K - K_MIN, K).encode("ascii"))
            process.stdin.flush()
            raw = _read_response_line(process, deadline)
            if raw != expected_derivation[K]:
                fail("semantic derivation replay differs at K={}".format(K))
            process.stdin.write(
                "V {} {} {}\n".format(K - K_MIN, K, attempt).encode(
                    "ascii"))
            process.stdin.flush()
            raw = _read_response_line(process, deadline)
            if raw != expected_validation[K]:
                fail("semantic validation replay differs at K={}".format(K))
        process.stdin.write(b"Q\n")
        process.stdin.flush()
        process.stdin.close()
        _drain_worker_stdout(process, deadline, "semantic replay worker")
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            fail("semantic replay deadline expired")
        try:
            process.wait(timeout=remaining)
        except subprocess.TimeoutExpired:
            fail("semantic replay deadline expired")
        stderr = _read_stderr_file(stderr_file, "semantic replay worker")
        if process.returncode != 0 or stderr:
            fail("semantic replay worker did not shut down cleanly: {}".format(
                stderr[:MAX_STDERR_BYTES].decode(
                    "utf-8", "replace").strip()))
        return roster_sha256
    finally:
        if process is not None:
            if process.poll() is None:
                _kill_process(process)
            for stream in (process.stdin, process.stdout, process.stderr):
                if stream is not None and not stream.closed:
                    stream.close()
        stderr_file.close()


def export_validated_repair_map_include(
        contract: Mapping[str, Any], repair_directory: Path,
        validation_directory: Path, worker: VerifiedWorker,
        expected_derivation_manifest_sha256: str,
        expected_validation_manifest_sha256: str,
        expected_derivation_audit_sha256: str,
        expected_validation_audit_sha256: str,
        output_path: Path) -> Mapping[str, Any]:
    """Publish a deterministic C++ initializer include after full validation."""
    # These caller-supplied digests are the external trust anchors for both
    # preregistered manifests and complete audits.  The bundle's unkeyed
    # self-hashes detect corruption and cross-wiring but, by themselves, cannot
    # authenticate a fully resealed artifact.  Record all four digests outside
    # the bundle at their prescribed freeze points before invoking export.
    expected_derivation_manifest_sha256 = _require_sha(
        expected_derivation_manifest_sha256,
        "externally recorded derivation manifest")
    expected_validation_manifest_sha256 = _require_sha(
        expected_validation_manifest_sha256,
        "externally recorded validation manifest")
    expected_derivation_audit_sha256 = _require_sha(
        expected_derivation_audit_sha256,
        "externally recorded derivation audit")
    expected_validation_audit_sha256 = _require_sha(
        expected_validation_audit_sha256,
        "externally recorded validation audit")
    _require_clean_source(worker.source_commit)
    if _sha256_open_fd(worker.descriptor, "mix2 repair worker") != \
            worker.sha256:
        fail("worker changed before repair-map export")
    repair_bundle = load_repair_bundle(contract, repair_directory)
    _require_derivation_manifest_anchor(
        repair_bundle, expected_derivation_manifest_sha256)
    repair_map = repair_bundle.repair_map
    repair_receipt = repair_bundle.receipt
    validation = load_validation_bundle(
        contract, validation_directory, repair_bundle,
        worker.sha256, worker.source_commit)
    if (validation.manifest["manifest_sha256"] !=
            expected_validation_manifest_sha256 or
            repair_receipt["derivation_audit_sha256"] !=
            expected_derivation_audit_sha256 or
            validation.summary["audit_sha256"] !=
            expected_validation_audit_sha256):
        fail("bundle differs from its externally recorded manifest or audit "
             "digests")
    semantic_replay_roster_sha256 = _replay_semantic_sample(
        worker, repair_map, repair_bundle.semantic_replay_records,
        validation.semantic_replay_records)
    payload = _render_repair_map_cpp_include(
        contract, repair_bundle, validation)
    if not output_path.name:
        fail("map include output must name one file")
    _require_clean_source(worker.source_commit)
    if (_controller_sha256() != repair_receipt["controller_sha256"] or
            _sha256_open_fd(worker.descriptor, "mix2 repair worker") !=
                worker.sha256):
        fail("controller, worker, or source changed during repair-map export")
    _publish_sources(output_path.parent, {output_path.name: payload})
    return {
        "schema": MAP_EXPORT_SCHEMA,
        "contract_sha256": contract_sha256(contract),
        "candidate_profile_sha256": candidate_profile_sha256(contract),
        "repair_map_sha256": repair_receipt["map_sha256"],
        "derivation_manifest_sha256":
            expected_derivation_manifest_sha256,
        "validation_manifest_sha256":
            expected_validation_manifest_sha256,
        "derivation_audit_sha256": expected_derivation_audit_sha256,
        "validation_audit_sha256": expected_validation_audit_sha256,
        "validation_summary_sha256":
            validation.summary["summary_sha256"],
        "semantic_replay_roster_sha256":
            semantic_replay_roster_sha256,
        "output_bytes": len(payload),
        "output_sha256": sha256_bytes(payload),
    }


def _identity_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--worker", type=Path, required=True)
    parser.add_argument("--worker-sha256", required=True)
    parser.add_argument("--source-commit", required=True)


def _short_screen_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--short-screen-dir", type=Path, required=True)
    parser.add_argument("--short-screen-summary-sha256", required=True)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    describe = commands.add_parser("describe")
    describe.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    manifest = commands.add_parser("manifest")
    _identity_args(manifest)
    _short_screen_args(manifest)
    derive = commands.add_parser("derive")
    _identity_args(derive)
    _short_screen_args(derive)
    derive.add_argument("--manifest-sha256", required=True)
    derive.add_argument("--output-dir", type=Path, required=True)
    derive.add_argument("--jobs", type=int, default=8)
    derive.add_argument("--deadline-seconds", type=float, default=172800.0)
    validation_manifest = commands.add_parser("validation-manifest")
    _identity_args(validation_manifest)
    validation_manifest.add_argument("--repair-dir", type=Path, required=True)
    validation_manifest.add_argument(
        "--derivation-manifest-sha256", required=True)
    validate = commands.add_parser("validate")
    _identity_args(validate)
    validate.add_argument("--repair-dir", type=Path, required=True)
    validate.add_argument("--derivation-manifest-sha256", required=True)
    validate.add_argument("--manifest-sha256", required=True)
    validate.add_argument("--output-dir", type=Path, required=True)
    validate.add_argument("--jobs", type=int, default=8)
    validate.add_argument("--deadline-seconds", type=float, default=172800.0)
    export_map = commands.add_parser("export-map-include")
    _identity_args(export_map)
    export_map.add_argument("--repair-dir", type=Path, required=True)
    export_map.add_argument("--validation-dir", type=Path, required=True)
    export_map.add_argument("--derivation-manifest-sha256", required=True)
    export_map.add_argument("--validation-manifest-sha256", required=True)
    export_map.add_argument("--derivation-audit-sha256", required=True)
    export_map.add_argument("--validation-audit-sha256", required=True)
    export_map.add_argument("--output", type=Path, required=True)
    return parser


def _result_exit_code(command: str, result: Mapping[str, Any]) -> int:
    if command != "validate":
        return 0
    disposition = result.get("disposition")
    if disposition == "PASS":
        return 0
    if disposition == "REJECT":
        return SCIENTIFIC_REJECT_EXIT_CODE
    fail("validation result omitted its PASS/REJECT disposition")
    return 1


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    worker: Optional[VerifiedWorker] = None
    try:
        contract = load_contract(args.contract)
        if args.command == "describe":
            result = {
                "schema": contract["schema"],
                "contract_id": contract["contract_id"],
                "contract_sha256": contract_sha256(contract),
                "candidate_profile_sha256": candidate_profile_sha256(contract),
                "map_bytes": MAP_BYTES,
                "selected_attempt_selection_cells":
                    MAP_BYTES * SELECTION_CELL_COUNT,
                "lower_attempt_witness_cells": "data-dependent",
                "validation_cells": MAP_BYTES * VALIDATION_CELL_COUNT,
                "validation_roster_sha256": VALIDATION_ROSTER_SHA256,
                "runtime_search": False,
            }
        else:
            worker = verify_worker(
                args.worker, args.worker_sha256, args.source_commit, contract)
            if args.command == "manifest":
                short_screen = load_short_screen_bundle(
                    contract, args.short_screen_dir, worker,
                    args.short_screen_summary_sha256)
                result = create_manifest(
                    "derive", contract, worker, short_screen.binding)
            elif args.command == "derive":
                short_screen = load_short_screen_bundle(
                    contract, args.short_screen_dir, worker,
                    args.short_screen_summary_sha256)
                result = _run_all_k(
                    "derive", contract, worker, args.output_dir, args.jobs,
                    args.deadline_seconds, args.manifest_sha256,
                    short_screen=short_screen)
            elif args.command == "export-map-include":
                result = export_validated_repair_map_include(
                    contract, args.repair_dir, args.validation_dir,
                    worker, args.derivation_manifest_sha256,
                    args.validation_manifest_sha256,
                    args.derivation_audit_sha256,
                    args.validation_audit_sha256, args.output)
            else:
                repair_bundle = load_repair_bundle(
                    contract, args.repair_dir)
                _require_derivation_manifest_anchor(
                    repair_bundle, args.derivation_manifest_sha256)
                repair_map = repair_bundle.repair_map
                receipt = repair_bundle.receipt
                if (receipt["worker_binary_sha256"] != worker.sha256 or
                        receipt["source_git_commit"] != worker.source_commit):
                    fail("validation worker differs from derivation worker")
                repair_identity = {
                    "map_sha256": receipt["map_sha256"],
                    "receipt_sha256": receipt["receipt_sha256"],
                    **_short_screen_binding(receipt),
                }
                if args.command == "validation-manifest":
                    result = create_manifest(
                        "validate", contract, worker, repair_identity,
                        repair_identity["map_sha256"],
                        repair_identity["receipt_sha256"])
                else:
                    result = _run_all_k(
                        "validate", contract, worker, args.output_dir,
                        args.jobs, args.deadline_seconds,
                        args.manifest_sha256, repair_map, repair_identity,
                        repair_bundle.shared_identity_sha256)
        print(canonical_json(result))
        return _result_exit_code(args.command, result)
    except (ContractError, OSError, ValueError, KeyError) as exc:
        invalid = {
            "schema": "wirehair.wh2.mix2-seed-repair-invalid.v2",
            "disposition": "INVALID",
            "diagnostic": str(exc),
        }
        print(canonical_json(invalid), file=sys.stderr)
        return 1
    finally:
        if worker is not None:
            os.close(worker.descriptor)


if __name__ == "__main__":
    raise SystemExit(main())
