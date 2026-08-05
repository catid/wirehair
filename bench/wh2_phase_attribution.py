#!/usr/bin/env python3
"""Strictly validate and analyze the bounded WH2 n16 phase screen.

This diagnostic is deliberately outside the canonical timing selector.  It
accepts exactly the 24 retry-zero P coordinates for the n16 profile, retains
weak cells as complete null ledgers, and reports phase attribution only.  It
contains no process orchestration or publication path: a separate sealed
runner must bind and terminally revalidate the completed recovery campaign's
same-source, same-binary timing-proxy witness before using these consumers.
"""

from __future__ import annotations

import hashlib
import math
from pathlib import Path
import re
import sys
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence
from typing import Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_benchmark_contract as contract_api
import wh2_native_short_screen as native_api
import wh2_run_native_short_screen as runner_api


PHASE_TRACE_SCHEMA = \
    "wirehair.wh2.native-phase-attribution-trace-record.v1"
PHASE_CELL_SCHEMA = "wirehair.wh2.phase-attribution-cell.v1"
PHASE_RECORD_SCHEMA = "wirehair.wh2.native-phase-attribution-record.v1"
PHASE_ANALYSIS_SCHEMA = "wirehair.wh2.phase-attribution-analysis.v1"
PHASE_DECISION_SCHEMA = "wirehair.wh2.phase-attribution-decision.v1"
NATIVE_WORK_SCHEMA = "wirehair.wh2.native-work.v1"
TRACE_MANIFEST_DOMAIN = b"wirehair.wh2.phase-attribution-trace-manifest.v1\0"
SOURCE_DOMAIN = b"wirehair.wh2.source.v1\0"
ZERO_SHA256 = "0" * 64

PROFILE_NAME = "n16"
PROFILE_ORDINAL = 0
INVOCATIONS_PER_SLOT = 16
PAIR_COUNT = 8
PHASE_COORDINATE_COUNT = 24
WORKER_COUNT = 8
FIXED_RECEIVED_OVERHEAD = 4
RECEIVE_OVERHEAD_CAP = 256
T95_DF7 = 2.3646242510102993
MAX_TRACE_BYTES = 256 * 1024
MAX_NATIVE_BYTES = 4 * 1024 * 1024
MAX_TRACE_LINE_BYTES = 8 * 1024
MAX_NATIVE_LINE_BYTES = 128 * 1024
UINT32_MAX = (1 << 32) - 1
UINT64_MAX = (1 << 64) - 1
INT63_MAX = (1 << 63) - 1

EXPECTED_K = (
    8, 32, 100, 128, 512, 1000, 2048, 5000, 8192, 20000, 32768, 64000,
)
EXPECTED_WIDTHS = (64, 1280)
EXPECTED_ARMS = (
    ("wirehair2_head", "wirehair2_certified",
     "4cafe27a8fb388ca9a4249b2c279b1406e7a0a86bcf14e98246988c7c503fa7a"),
    ("wirehair1", "wirehair1",
     "d5a24d404e69efeb439907cd8271eba98d6af86b58efe159a820fb7aea08883d"),
    ("wirehair2_dense_two07_basis_v1", "wirehair2_experiment",
     "9527f200ad38c7eec6502b2f768fdd67b92787fb227eed3d7616274ffc2df388"),
)
LEFT_ARM = EXPECTED_ARMS[2][0]
RIGHT_ARM = EXPECTED_ARMS[0][0]

TRACE_FIELDS = frozenset((
    "candidate_count", "coordinate_ordinal", "packet_count",
    "phase_cell_sha256_by_profile", "schema",
    "source_base_cell_sha256", "trace_qualified_timing_cell_sha256",
    "trace_sha256",
))
NATIVE_RECORD_FIELDS = frozenset((
    "cpu", "finished_monotonic_ns", "message_sha256", "ordinal", "payload",
    "schema", "started_monotonic_ns", "work_sha256",
    "worker_binary_sha256", "worker_pid", "worker_process_start_ticks",
))
PHASE_PAYLOAD_FIELDS = frozenset((
    "K", "band", "base_loss_seed", "binary_sha256", "block_bytes",
    "block_muladds_semantics", "cell_sha256", "construction_attempt",
    "coordinate_ordinal", "fixed_received_overhead", "interleave_policy",
    "invocations_per_slot", "left_arm", "left_arm_descriptor_sha256",
    "left_non_timing_counters", "left_realized_construction_sha256",
    "left_repair_map_sha256", "loss_ppm", "loss_retry_offset", "loss_seed",
    "measured_observations", "order", "panel_comparable", "panel_kind",
    "profile", "profile_ordinal", "receive_overhead_cap", "replicate",
    "right_arm", "right_arm_descriptor_sha256",
    "right_non_timing_counters", "right_realized_construction_sha256",
    "right_repair_map_sha256", "schedule", "scope", "slot_sums",
    "source_base_cell_sha256", "trace_qualified_timing_cell_sha256",
    "trace_sha256", "warmups", "weak_ledger",
))
COUNTER_FIELDS = frozenset((
    "binary_adjacency_storage_allocations",
    "binary_adjacency_storage_bytes", "binary_residual_rank",
    "binary_row_references", "binary_row_storage_allocations",
    "binary_row_storage_bytes", "block_xors",
    "full_block_gf256_multiply_add_divide_normalize_ops",
    "inactivated_columns", "packet_rows", "packet_seed_attempt",
    "peeled_columns", "residual_rank", "residual_rows",
))
COUNTER_UINT32_FIELDS = frozenset((
    "binary_adjacency_storage_allocations", "binary_residual_rank",
    "binary_row_storage_allocations", "inactivated_columns", "packet_rows",
    "packet_seed_attempt", "peeled_columns", "residual_rank",
    "residual_rows",
))
COUNTER_UINT64_FIELDS = COUNTER_FIELDS - COUNTER_UINT32_FIELDS
OBSERVATION_FIELDS = frozenset((
    "arm", "back_sub_ns", "build_ns", "bytes_verified", "counter_sha256",
    "outcome", "outer_ns", "peel_ns", "project_ns", "residual_ns",
    "wirehair_result",
))
MEASURED_FIELDS = frozenset(("block", "observation", "repeat", "slot"))
SLOT_FIELDS = frozenset((
    "back_sub_ns", "build_ns", "outer_ns", "peel_ns", "project_ns",
    "residual_ns", "slot",
))
PHASE_KEYS = (
    ("outer", "outer_ns"),
    ("build", "build_ns"),
    ("peel", "peel_ns"),
    ("project", "project_ns"),
    ("residual", "residual_ns"),
    ("back_sub", "back_sub_ns"),
)
TIMING_KEYS = tuple(value for _, value in PHASE_KEYS)
DESCRIPTION_FIELDS = frozenset((
    "schema", "source_git_commit", "binary_sha256", "arms",
))
DESCRIPTION_ARM_FIELDS = frozenset((
    "arm", "codec", "arm_descriptor_sha256",
))
GIT_COMMIT = re.compile(r"[0-9a-f]{40}\Z")


class PhaseRunnerError(runner_api.RunnerError):
    """The phase-attribution evidence cannot continue safely."""


def fail(message: str) -> None:
    raise PhaseRunnerError(message)


def _is_sha256(value: Any) -> bool:
    return isinstance(value, str) and \
        contract_api.SHA256.fullmatch(value) is not None


def _is_uint32(value: Any) -> bool:
    return type(value) is int and 0 <= value <= UINT32_MAX


def _is_uint64(value: Any) -> bool:
    return type(value) is int and 0 <= value <= UINT64_MAX


def _is_int63(value: Any) -> bool:
    return type(value) is int and 0 <= value <= INT63_MAX


def _exact_mapping(value: Any, fields: Iterable[str], context: str) \
        -> Mapping[str, Any]:
    if not isinstance(value, dict) or set(value) != set(fields):
        fail("{} has an unexpected schema".format(context))
    return value


def _parse_jsonl_bytes(
        data: bytes, context: str, expected_count: int,
        max_file_bytes: int, max_line_bytes: int) -> List[Mapping[str, Any]]:
    if (type(data) is not bytes or not data or len(data) > max_file_bytes or
            not data.endswith(b"\n") or
            data.count(b"\n") != expected_count):
        fail("{} has invalid bounded cardinality".format(context))
    rows = []
    for index, line in enumerate(data.splitlines(keepends=True)):
        if len(line) > max_line_bytes:
            fail("{} line {} exceeds its bound".format(context, index + 1))
        try:
            rows.append(runner_api._parse_canonical_line(
                line, "{} line {}".format(context, index + 1)))
        except runner_api.RunnerError as exc:
            fail(str(exc))
    return rows


def _phase_coordinates(contract: Mapping[str, Any]) \
        -> List[Mapping[str, Any]]:
    bases = list(contract_api.iter_timing_base_cells(contract, "development"))
    if len(bases) != 2304:
        fail("phase projection requires the exact development timing domain")
    coordinates = []
    for width_index, width in enumerate(EXPECTED_WIDTHS):
        for k_index, K in enumerate(EXPECTED_K):
            source_ordinal = width_index * 96 + k_index * 8
            base = bases[source_ordinal]
            if (base.get("phase") != "development" or base.get("K") != K or
                    base.get("block_bytes") != width or
                    base.get("replicate") != 0 or
                    base.get("base_seed_attempt") != 0 or
                    base.get("loss_ppm") != 100000 or
                    base.get("schedule") != "iid" or
                    base.get("fixed_received_overhead") !=
                        FIXED_RECEIVED_OVERHEAD or
                    base.get("receive_overhead_cap") !=
                        RECEIVE_OVERHEAD_CAP or
                    base.get("interleave_policy") !=
                        "self-counterbalanced-repeat-major-v1"):
                fail("phase coordinate differs from its natural retry-zero projection")
            qualified = dict(base)
            qualified["base_cell_sha256"] = contract_api.sha256_json(base)
            qualified["loss_retry_offset"] = 0
            qualified["loss_seed"] = base["base_loss_seed"]
            coordinates.append({
                "coordinate_ordinal": len(coordinates),
                "source_ordinal": source_ordinal,
                "base": base,
                "qualified": qualified,
            })
    if len(coordinates) != PHASE_COORDINATE_COUNT:
        fail("phase projection has the wrong coordinate count")
    return coordinates


def _splitmix64_next(state: int) -> Tuple[int, int]:
    mask = (1 << 64) - 1
    state = (state + 0x9e3779b97f4a7c15) & mask
    value = state
    value = ((value ^ (value >> 30)) * 0xbf58476d1ce4e5b9) & mask
    value = ((value ^ (value >> 27)) * 0x94d049bb133111eb) & mask
    return state, (value ^ (value >> 31)) & mask


def _deterministic_message_sha256(base: Mapping[str, Any]) -> str:
    base_json = contract_api.canonical_json(base).encode("ascii")
    source = hashlib.sha256(SOURCE_DOMAIN + base_json).digest()
    source_seed = int.from_bytes(source[:8], "big")
    K = base["K"]
    block_bytes = base["block_bytes"]
    mask = (1 << 64) - 1
    state = source_seed ^ ((K * 0xd6e8feb86659fd93) & mask) ^ \
        ((block_bytes * 0xa0761d6478bd642f) & mask)
    remaining = K * block_bytes
    digest = hashlib.sha256()
    chunk = bytearray()
    while remaining:
        state, word = _splitmix64_next(state)
        encoded = word.to_bytes(8, "little")
        take = min(remaining, 8)
        chunk.extend(encoded[:take])
        remaining -= take
        if len(chunk) >= 1024 * 1024 or remaining == 0:
            digest.update(chunk)
            chunk.clear()
    return digest.hexdigest()


def _validate_description(value: Mapping[str, Any],
                          binary_sha256: Optional[str] = None) \
        -> Mapping[str, Any]:
    _exact_mapping(value, DESCRIPTION_FIELDS, "phase worker description")
    source_git_commit = value.get("source_git_commit")
    if (value.get("schema") != runner_api.DESCRIPTION_SCHEMA or
            not isinstance(source_git_commit, str) or
            GIT_COMMIT.fullmatch(source_git_commit) is None or
            not _is_sha256(value.get("binary_sha256")) or
            (binary_sha256 is not None and
             value.get("binary_sha256") != binary_sha256)):
        fail("phase worker description has invalid source/binary identity")
    arms = value.get("arms")
    if not isinstance(arms, list) or len(arms) != len(EXPECTED_ARMS):
        fail("phase worker description has the wrong exact arm roster")
    for index, expected in enumerate(EXPECTED_ARMS):
        arm = _exact_mapping(
            arms[index], DESCRIPTION_ARM_FIELDS,
            "phase worker arm {}".format(index))
        if (arm.get("arm"), arm.get("codec"),
                arm.get("arm_descriptor_sha256")) != expected:
            fail("phase worker arm {} differs from its frozen descriptor".format(
                index))
    return value


def _realized_sha256(
        arm: Mapping[str, Any], K: int, block_bytes: int) -> str:
    if arm["codec"] == "wirehair2_certified":
        return contract_api.realized_construction_sha256(
            arm["arm_descriptor_sha256"], K, block_bytes, 0)
    return contract_api.sha256_json({
        "K": K,
        "arm_descriptor_sha256": arm["arm_descriptor_sha256"],
        "block_bytes": block_bytes,
        "codec": arm["codec"],
        "construction_attempt": 0,
        "schema": "wirehair.wh2.realized-construction.v1",
    })


def _phase_cell(
        coordinate: Mapping[str, Any], profile_ordinal: int,
        description: Mapping[str, Any], trace_sha256: str) \
        -> Mapping[str, Any]:
    if type(profile_ordinal) is not int or profile_ordinal not in (0, 1):
        fail("phase cell profile is outside the closed worker domain")
    qualified = coordinate["qualified"]
    profile = ("n16", "n24")[profile_ordinal]
    n = (16, 24)[profile_ordinal]
    left = description["arms"][2]
    right = description["arms"][0]
    return {
        "K": qualified["K"],
        "band": qualified["band"],
        "base_loss_seed": qualified["base_loss_seed"],
        "block_bytes": qualified["block_bytes"],
        "construction_attempt": 0,
        "coordinate_ordinal": coordinate["coordinate_ordinal"],
        "diagnostic_phase": "development",
        "fixed_received_overhead": FIXED_RECEIVED_OVERHEAD,
        "interleave_policy": "self-counterbalanced-repeat-major-v1",
        "invocations_per_slot": n,
        "left_arm": left["arm"],
        "left_arm_descriptor_sha256": left["arm_descriptor_sha256"],
        "left_realized_construction_sha256": _realized_sha256(
            left, qualified["K"], qualified["block_bytes"]),
        "left_repair_map_sha256": ZERO_SHA256,
        "loss_ppm": qualified["loss_ppm"],
        "loss_retry_offset": 0,
        "loss_seed": qualified["loss_seed"],
        "order": "ABBA",
        "panel_kind": "ab",
        "profile": profile,
        "profile_ordinal": profile_ordinal,
        "receive_overhead_cap": RECEIVE_OVERHEAD_CAP,
        "replicate": 0,
        "right_arm": right["arm"],
        "right_arm_descriptor_sha256": right["arm_descriptor_sha256"],
        "right_realized_construction_sha256": _realized_sha256(
            right, qualified["K"], qualified["block_bytes"]),
        "right_repair_map_sha256": ZERO_SHA256,
        "schedule": "iid",
        "schema": PHASE_CELL_SCHEMA,
        "scope": "decoder_solve",
        "source_base_cell_sha256": qualified["base_cell_sha256"],
        "trace_qualified_timing_cell_sha256":
            contract_api.sha256_json(qualified),
        "trace_sha256": trace_sha256,
    }


def validate_trace_manifest(
        contract: Mapping[str, Any], description: Mapping[str, Any],
        data: bytes) -> Tuple[List[Mapping[str, Any]], str]:
    """Validate all 24 trace rows and both worker-domain profile cell hashes."""
    _validate_description(description)
    rows = _parse_jsonl_bytes(
        data, "phase trace manifest", PHASE_COORDINATE_COUNT,
        MAX_TRACE_BYTES, MAX_TRACE_LINE_BYTES)
    coordinates = _phase_coordinates(contract)
    for ordinal, (row, coordinate) in enumerate(zip(rows, coordinates)):
        _exact_mapping(row, TRACE_FIELDS, "phase trace record")
        qualified = coordinate["qualified"]
        profile_hashes = row.get("phase_cell_sha256_by_profile")
        packet_count = row.get("packet_count")
        candidate_count = row.get("candidate_count")
        if (row.get("schema") != PHASE_TRACE_SCHEMA or
                not _is_uint32(row.get("coordinate_ordinal")) or
                row.get("coordinate_ordinal") != ordinal or
                not _is_uint32(packet_count) or
                packet_count != qualified["K"] + RECEIVE_OVERHEAD_CAP or
                not _is_uint64(candidate_count) or
                not row["packet_count"] <= row["candidate_count"] <=
                    256 * row["packet_count"] + 65536 or
                row.get("source_base_cell_sha256") !=
                    qualified["base_cell_sha256"] or
                row.get("trace_qualified_timing_cell_sha256") !=
                    contract_api.sha256_json(qualified) or
                not _is_sha256(row.get("trace_sha256")) or
                not isinstance(profile_hashes, list) or
                len(profile_hashes) != 2):
            fail("phase trace record {} differs from its frozen cell".format(
                ordinal))
        expected_hashes = [
            contract_api.sha256_json(_phase_cell(
                coordinate, profile, description, row["trace_sha256"]))
            for profile in range(2)
        ]
        if profile_hashes != expected_hashes:
            fail("phase trace cell hashes do not bind both profiles")
    all_cells = [
        digest for row in rows
        for digest in row["phase_cell_sha256_by_profile"]
    ]
    if len(set(all_cells)) != 48:
        fail("phase trace manifest reuses a profile cell hash")
    digest = hashlib.sha256(TRACE_MANIFEST_DOMAIN + data).hexdigest()
    return rows, digest


def _validate_counters(value: Any, context: str) -> Mapping[str, Any]:
    counters = _exact_mapping(value, COUNTER_FIELDS, context)
    if (any(not _is_uint32(counters[field])
            for field in COUNTER_UINT32_FIELDS) or
            any(not _is_uint64(counters[field])
                for field in COUNTER_UINT64_FIELDS)):
        fail("{} contains a noncanonical counter".format(context))
    return counters


def _expected_observation_arm(block: int, slot: int) -> str:
    primary_left = slot % 4 in (0, 3)
    left = primary_left if block == 0 else not primary_left
    return LEFT_ARM if left else RIGHT_ARM


def _validate_observation(
        value: Any, counters_by_arm: Mapping[str, Mapping[str, Any]],
        timing_visible: bool, expected_arm: str, context: str) \
        -> Mapping[str, Any]:
    observation = _exact_mapping(value, OBSERVATION_FIELDS, context)
    if observation.get("arm") != expected_arm:
        fail("{} has the wrong chronological arm".format(context))
    result = observation.get("wirehair_result")
    outcome = observation.get("outcome")
    valid_outcome = type(result) is int and (
        (result == 0 and outcome == "success") or
        (result == 1 and outcome == "need_more_at_cap") or
        (result in (3, 4) and outcome == "construct_failed")
    )
    counters = counters_by_arm.get(expected_arm)
    if (not valid_outcome or counters is None or
            observation.get("counter_sha256") !=
                contract_api.sha256_json(counters) or
            type(observation.get("bytes_verified")) is not bool):
        fail("{} has invalid outcome/counter identity".format(context))
    if result == 0 and not observation["bytes_verified"]:
        fail("{} successful solve did not verify bytes".format(context))
    if result != 0 and observation["bytes_verified"]:
        fail("{} failed solve claims verified bytes".format(context))
    timings = [observation.get(key) for key in TIMING_KEYS]
    if timing_visible:
        if (result != 0 or any(not _is_int63(item) for item in timings) or
                observation["outer_ns"] <= 0 or
                sum(observation[key] for _, key in PHASE_KEYS[1:]) >
                    observation["outer_ns"]):
            fail("{} has invalid visible phase timing".format(context))
    elif any(item is not None for item in timings):
        fail("{} exposes timing from a weak panel".format(context))
    return observation


def _stable_observation_identity(value: Mapping[str, Any]) \
        -> Tuple[Any, ...]:
    return (
        value["arm"], value["wirehair_result"], value["outcome"],
        value["bytes_verified"], value["counter_sha256"],
    )


def _validate_payload(
        payload: Any, coordinate: Mapping[str, Any],
        trace: Mapping[str, Any], description: Mapping[str, Any]) \
        -> Mapping[str, Any]:
    payload = _exact_mapping(
        payload, PHASE_PAYLOAD_FIELDS, "phase result payload")
    qualified = coordinate["qualified"]
    expected_cell = _phase_cell(
        coordinate, PROFILE_ORDINAL, description, trace["trace_sha256"])
    expected_scalars = {
        "K": qualified["K"],
        "band": qualified["band"],
        "base_loss_seed": qualified["base_loss_seed"],
        "binary_sha256": description["binary_sha256"],
        "block_bytes": qualified["block_bytes"],
        "block_muladds_semantics":
            "full-block-gf256-multiply-add-divide-normalize-operations",
        "cell_sha256": contract_api.sha256_json(expected_cell),
        "construction_attempt": 0,
        "coordinate_ordinal": coordinate["coordinate_ordinal"],
        "fixed_received_overhead": FIXED_RECEIVED_OVERHEAD,
        "interleave_policy": "self-counterbalanced-repeat-major-v1",
        "invocations_per_slot": INVOCATIONS_PER_SLOT,
        "left_arm": expected_cell["left_arm"],
        "left_arm_descriptor_sha256":
            expected_cell["left_arm_descriptor_sha256"],
        "left_realized_construction_sha256":
            expected_cell["left_realized_construction_sha256"],
        "left_repair_map_sha256": ZERO_SHA256,
        "loss_ppm": qualified["loss_ppm"],
        "loss_retry_offset": 0,
        "loss_seed": qualified["loss_seed"],
        "order": "ABBA",
        "panel_kind": "ab",
        "profile": PROFILE_NAME,
        "profile_ordinal": PROFILE_ORDINAL,
        "receive_overhead_cap": RECEIVE_OVERHEAD_CAP,
        "replicate": 0,
        "right_arm": expected_cell["right_arm"],
        "right_arm_descriptor_sha256":
            expected_cell["right_arm_descriptor_sha256"],
        "right_realized_construction_sha256":
            expected_cell["right_realized_construction_sha256"],
        "right_repair_map_sha256": ZERO_SHA256,
        "schedule": "iid",
        "scope": "decoder_solve",
        "source_base_cell_sha256": qualified["base_cell_sha256"],
        "trace_qualified_timing_cell_sha256":
            contract_api.sha256_json(qualified),
        "trace_sha256": trace["trace_sha256"],
    }
    if any(
            payload.get(field) != expected or
            (type(expected) is int and type(payload.get(field)) is not int)
            for field, expected in expected_scalars.items()):
        fail("phase payload differs from its exact retry-zero cell")
    if (trace["phase_cell_sha256_by_profile"][PROFILE_ORDINAL] !=
            payload["cell_sha256"]):
        fail("phase payload cell hash differs from the trace manifest")

    left_counters = _validate_counters(
        payload["left_non_timing_counters"], "left non-timing counters")
    right_counters = _validate_counters(
        payload["right_non_timing_counters"], "right non-timing counters")
    counters_by_arm = {
        LEFT_ARM: left_counters,
        RIGHT_ARM: right_counters,
    }
    comparable = payload.get("panel_comparable")
    weak_ledger = payload.get("weak_ledger")
    if type(comparable) is not bool or type(weak_ledger) is not bool or \
            weak_ledger == comparable:
        fail("phase panel comparability/weak-ledger flags are inconsistent")

    warmups = _exact_mapping(
        payload.get("warmups"), ("left", "right"), "phase warmups")
    left_warm = _validate_observation(
        warmups["left"], counters_by_arm, comparable, LEFT_ARM, "left warmup")
    right_warm = _validate_observation(
        warmups["right"], counters_by_arm, comparable, RIGHT_ARM,
        "right warmup")
    if comparable != (
            left_warm["wirehair_result"] == 0 and
            right_warm["wirehair_result"] == 0):
        fail("phase panel comparable flag differs from its warmups")
    warm_by_arm = {
        LEFT_ARM: _stable_observation_identity(left_warm),
        RIGHT_ARM: _stable_observation_identity(right_warm),
    }

    measured = payload.get("measured_observations")
    if not isinstance(measured, list) or len(measured) != 4 * \
            INVOCATIONS_PER_SLOT:
        fail("phase payload lacks the exact n16 measured chronology")
    expected_tags = [
        (block, repeat, block * 4 + block_slot)
        for block in range(2)
        for repeat in range(PAIR_COUNT)
        for block_slot in range(4)
    ]
    observations_by_slot: Dict[int, List[Mapping[str, Any]]] = {
        slot: [] for slot in range(8)
    }
    for index, (raw, tag) in enumerate(zip(measured, expected_tags)):
        row = _exact_mapping(
            raw, MEASURED_FIELDS, "phase measured row {}".format(index))
        actual_tag = (
            row.get("block"), row.get("repeat"), row.get("slot"))
        if (any(type(item) is not int for item in actual_tag) or
                actual_tag != tag):
            fail("phase measured chronology is not exact ABBA/BAAB n16")
        expected_arm = _expected_observation_arm(tag[0], tag[2])
        observation = _validate_observation(
            row["observation"], counters_by_arm, comparable, expected_arm,
            "phase measured observation {}".format(index))
        if _stable_observation_identity(observation) != \
                warm_by_arm[expected_arm]:
            fail("phase measured outcome/counter identity drifted")
        observations_by_slot[tag[2]].append(observation)

    slot_sums = payload.get("slot_sums")
    if not isinstance(slot_sums, list) or len(slot_sums) != 8:
        fail("phase payload lacks exactly eight slot ledgers")
    for slot, raw in enumerate(slot_sums):
        totals = _exact_mapping(
            raw, SLOT_FIELDS, "phase slot {}".format(slot))
        if (not _is_uint32(totals.get("slot")) or
                totals.get("slot") != slot):
            fail("phase slot ledger is out of order")
        values = [totals.get(key) for key in TIMING_KEYS]
        if comparable:
            if (any(not _is_int63(value) for value in values) or
                    totals["outer_ns"] <= 0 or
                    sum(totals[key] for _, key in PHASE_KEYS[1:]) >
                        totals["outer_ns"]):
                fail("phase slot has invalid visible timing totals")
            for key in TIMING_KEYS:
                if totals[key] != sum(
                        observation[key]
                        for observation in observations_by_slot[slot]):
                    fail("phase slot totals differ from measured observations")
        elif any(value is not None for value in values):
            fail("weak phase slot ledger exposes timing")
    return payload


def _expected_work_sha256(cell_sha256: str, ordinal: int) -> str:
    return contract_api.sha256_json({
        "cell_sha256": cell_sha256,
        "evidence_kind": "phase_attribution",
        "ordinal": ordinal,
        "phase": "development",
        "schema": NATIVE_WORK_SCHEMA,
    })


def validate_native_records(
        contract: Mapping[str, Any], description: Mapping[str, Any],
        trace_data: bytes, data: bytes,
        coordinate_cpus: Sequence[int],
        runtime_workers: Mapping[int, Tuple[int, int]],
        window_start_ns: int, window_end_ns: int,
        ) -> Tuple[List[Mapping[str, Any]], Mapping[str, Any]]:
    """Strictly validate and order the exact 24 n16 native envelopes."""
    _validate_description(description)
    traces, trace_manifest_sha256 = validate_trace_manifest(
        contract, description, trace_data)
    if not isinstance(coordinate_cpus, (list, tuple)):
        fail("phase coordinate CPU schedule must be a concrete sequence")
    worker_cpus = tuple(coordinate_cpus[:WORKER_COUNT])
    if (len(coordinate_cpus) != PHASE_COORDINATE_COUNT or
            len(worker_cpus) != WORKER_COUNT or
            any(not _is_uint32(cpu) for cpu in coordinate_cpus) or
            len(set(worker_cpus)) != WORKER_COUNT or
            tuple(coordinate_cpus) != worker_cpus * 3 or
            not _is_uint64(window_start_ns) or window_start_ns == 0 or
            not _is_uint64(window_end_ns) or
            window_end_ns <= window_start_ns or
            not isinstance(runtime_workers, Mapping) or
            set(runtime_workers) != set(worker_cpus)):
        fail("phase execution geometry is not the exact 24-job/eight-worker domain")
    runtime_pids = []
    for cpu in worker_cpus:
        identity = runtime_workers[cpu]
        if (not isinstance(identity, tuple) or len(identity) != 2 or
                not _is_uint32(identity[0]) or identity[0] == 0 or
                not _is_uint64(identity[1]) or identity[1] == 0):
            fail("phase live-worker roster has invalid PID/start identity")
        runtime_pids.append(identity[0])
    if len(set(runtime_pids)) != WORKER_COUNT:
        fail("phase live-worker roster reuses a PID")
    rows = _parse_jsonl_bytes(
        data, "phase native result stream", PHASE_COORDINATE_COUNT,
        MAX_NATIVE_BYTES, MAX_NATIVE_LINE_BYTES)
    coordinates = _phase_coordinates(contract)
    expected_messages = [
        _deterministic_message_sha256(coordinate["base"])
        for coordinate in coordinates
    ]

    ordered: List[Optional[Mapping[str, Any]]] = \
        [None] * PHASE_COORDINATE_COUNT
    cpu_to_pid: Dict[int, int] = {}
    pid_to_cpu: Dict[int, int] = {}
    pid_to_ticks: Dict[int, int] = {}
    worker_start: Optional[int] = None
    worker_end: Optional[int] = None
    work_hashes: List[Optional[str]] = [None] * PHASE_COORDINATE_COUNT
    message_hashes: List[Optional[str]] = [None] * PHASE_COORDINATE_COUNT
    for row in rows:
        _exact_mapping(row, NATIVE_RECORD_FIELDS, "phase native envelope")
        if (row.get("schema") != PHASE_RECORD_SCHEMA or
                not _is_uint32(row.get("cpu")) or
                not _is_uint32(row.get("ordinal")) or
                not _is_uint32(row.get("worker_pid")) or
                not _is_uint64(row.get("worker_process_start_ticks")) or
                not _is_uint64(row.get("started_monotonic_ns")) or
                not _is_uint64(row.get("finished_monotonic_ns"))):
            fail("phase native envelope has invalid schema/integer fields")
        ordinal = row["ordinal"]
        if ordinal < 0 or ordinal % 2 != PROFILE_ORDINAL:
            fail("phase native ordinal is not an n16 worker-domain ordinal")
        coordinate_ordinal = ordinal // 2
        if (coordinate_ordinal >= PHASE_COORDINATE_COUNT or
                ordered[coordinate_ordinal] is not None):
            fail("phase native stream has a duplicate/out-of-domain ordinal")
        cpu = row["cpu"]
        pid = row["worker_pid"]
        ticks = row["worker_process_start_ticks"]
        start = row["started_monotonic_ns"]
        end = row["finished_monotonic_ns"]
        if (cpu != coordinate_cpus[coordinate_ordinal] or pid <= 0 or
                ticks <= 0 or not window_start_ns <= start < end or
                end > window_end_ns):
            fail("phase native envelope has invalid CPU/PID/time provenance")
        try:
            native_api._require_process_predates_window(
                ticks, window_start_ns, "phase native worker")
        except native_api.NativeEvidenceError as exc:
            fail(str(exc))
        if (cpu in cpu_to_pid and cpu_to_pid[cpu] != pid) or \
                (pid in pid_to_cpu and pid_to_cpu[pid] != cpu):
            fail("phase persistent worker CPU/PID identity drifted")
        cpu_to_pid[cpu] = pid
        pid_to_cpu[pid] = cpu
        if pid_to_ticks.setdefault(pid, ticks) != ticks:
            fail("phase persistent worker start ticks drifted")
        if runtime_workers.get(cpu) != (pid, ticks):
            fail("phase native envelope differs from its live worker")
        for field in (
                "worker_binary_sha256", "message_sha256", "work_sha256"):
            if not _is_sha256(row.get(field)):
                fail("phase native envelope has an invalid {} hash".format(
                    field))
        if (row["worker_binary_sha256"] != description["binary_sha256"] or
                row["message_sha256"] !=
                    expected_messages[coordinate_ordinal]):
            fail("phase native envelope differs from its source/binary identity")
        payload = _validate_payload(
            row.get("payload"), coordinates[coordinate_ordinal],
            traces[coordinate_ordinal], description)
        expected_ordinal = coordinate_ordinal * 2 + PROFILE_ORDINAL
        expected_work = _expected_work_sha256(
            payload["cell_sha256"], expected_ordinal)
        if ordinal != expected_ordinal or row["work_sha256"] != expected_work:
            fail("phase native work hash differs from its exact command")
        ordered[coordinate_ordinal] = payload
        work_hashes[coordinate_ordinal] = row["work_sha256"]
        message_hashes[coordinate_ordinal] = row["message_sha256"]
        worker_start = start if worker_start is None else min(
            worker_start, start)
        worker_end = end if worker_end is None else max(worker_end, end)
    if (any(value is None for value in ordered) or
            set(cpu_to_pid) != set(coordinate_cpus)):
        fail("phase native stream is incomplete or did not use every worker")
    metadata = {
        "coordinate_cpus": list(coordinate_cpus),
        "message_set_sha256": contract_api.sha256_json(message_hashes),
        "native_stream_sha256": hashlib.sha256(data).hexdigest(),
        "record_count": PHASE_COORDINATE_COUNT,
        "trace_manifest_sha256": trace_manifest_sha256,
        "worker_cpus": list(worker_cpus),
        "worker_end_monotonic_ns": worker_end,
        "worker_start_monotonic_ns": worker_start,
        "workers": [
            {
                "cpu": cpu,
                "pid": cpu_to_pid[cpu],
                "process_start_ticks": pid_to_ticks[cpu_to_pid[cpu]],
            }
            for cpu in worker_cpus
        ],
        "work_set_sha256": contract_api.sha256_json(work_hashes),
    }
    return [value for value in ordered if value is not None], metadata


def _cell_phase_units(payload: Mapping[str, Any]) -> Mapping[str, Any]:
    left_counters = payload["left_non_timing_counters"]
    right_counters = payload["right_non_timing_counters"]
    warmups = payload["warmups"]
    identity = {
        "K": payload["K"],
        "band": payload["band"],
        "block_bytes": payload["block_bytes"],
        "cell_sha256": payload["cell_sha256"],
        "coordinate_ordinal": payload["coordinate_ordinal"],
        "left_non_timing_counters": dict(left_counters),
        "left_minus_right_non_timing_counter_deltas": {
            field: left_counters[field] - right_counters[field]
            for field in sorted(COUNTER_FIELDS)
        },
        "arm_outcome_ledger": {
            "left": {
                "arm": warmups["left"]["arm"],
                "bytes_verified": warmups["left"]["bytes_verified"],
                "counter_sha256": warmups["left"]["counter_sha256"],
                "outcome": warmups["left"]["outcome"],
                "wirehair_result": warmups["left"]["wirehair_result"],
            },
            "right": {
                "arm": warmups["right"]["arm"],
                "bytes_verified": warmups["right"]["bytes_verified"],
                "counter_sha256": warmups["right"]["counter_sha256"],
                "outcome": warmups["right"]["outcome"],
                "wirehair_result": warmups["right"]["wirehair_result"],
            },
        },
        "panel_comparable": payload["panel_comparable"],
        "profile": PROFILE_NAME,
        "profile_ordinal": PROFILE_ORDINAL,
        "right_non_timing_counters": dict(right_counters),
        "source_base_cell_sha256": payload["source_base_cell_sha256"],
        "trace_qualified_timing_cell_sha256":
            payload["trace_qualified_timing_cell_sha256"],
        "trace_sha256": payload["trace_sha256"],
        "weak_ledger": payload["weak_ledger"],
    }
    by_tag = {
        (row["block"], row["repeat"], row["slot"]): row["observation"]
        for row in payload["measured_observations"]
    }
    phases = {}
    for name, key in PHASE_KEYS:
        if not payload["panel_comparable"]:
            phases[name] = {
                "mean_log_ratio": None,
                "pair_log_ratios": None,
                "point_ratio": None,
                "status": "non_comparable_weak",
            }
            continue
        pair_logs = []
        unresolved = False
        for repeat in range(PAIR_COUNT):
            primary_left = sum(
                by_tag[(0, repeat, slot)][key] for slot in (0, 3))
            primary_right = sum(
                by_tag[(0, repeat, slot)][key] for slot in (1, 2))
            opposite_left = sum(
                by_tag[(1, repeat, slot)][key] for slot in (5, 6))
            opposite_right = sum(
                by_tag[(1, repeat, slot)][key] for slot in (4, 7))
            if min(primary_left, primary_right,
                   opposite_left, opposite_right) <= 0:
                unresolved = True
                break
            pair_logs.append(0.5 * (
                math.log(primary_left / primary_right) +
                math.log(opposite_left / opposite_right)))
        if unresolved:
            phases[name] = {
                "mean_log_ratio": None,
                "pair_log_ratios": None,
                "point_ratio": None,
                "status": "unresolved_zero_phase",
            }
        else:
            mean = math.fsum(pair_logs) / PAIR_COUNT
            phases[name] = {
                "mean_log_ratio": mean,
                "pair_log_ratios": pair_logs,
                "point_ratio": math.exp(mean),
                "status": "comparable",
            }
    identity["phases"] = phases
    return identity


def _confidence_interval(values: Sequence[float]) -> Mapping[str, float]:
    if (len(values) != PAIR_COUNT or
            any(type(value) is not float or not math.isfinite(value)
                for value in values)):
        fail("phase confidence interval requires eight finite pair means")
    mean = math.fsum(values) / PAIR_COUNT
    variance = math.fsum(
        (value - mean) * (value - mean) for value in values) / \
        (PAIR_COUNT - 1)
    half_width = T95_DF7 * math.sqrt(variance / PAIR_COUNT)
    lower = mean - half_width
    upper = mean + half_width
    return {
        "geometric_mean_ratio": math.exp(mean),
        "lower_95_log_ratio": lower,
        "lower_95_ratio": math.exp(lower),
        "mean_log_ratio": mean,
        "upper_95_log_ratio": upper,
        "upper_95_ratio": math.exp(upper),
    }


def build_phase_analysis(
        payloads: Sequence[Mapping[str, Any]]) -> Mapping[str, Any]:
    """Analyze eight repeat-pair units without treating cells as iid."""
    if (len(payloads) != PHASE_COORDINATE_COUNT or
            [payload.get("coordinate_ordinal") for payload in payloads] !=
                list(range(PHASE_COORDINATE_COUNT)) or
            any(payload.get("profile_ordinal") != PROFILE_ORDINAL
                for payload in payloads)):
        fail("phase analysis requires the exact ordered n16 payload domain")
    cells = [_cell_phase_units(payload) for payload in payloads]
    groups: List[Tuple[str, Any, List[Mapping[str, Any]]]] = [
        ("global", "all", cells),
    ]
    groups.extend((
        "width", width,
        [cell for cell in cells if cell["block_bytes"] == width])
        for width in EXPECTED_WIDTHS)
    groups.extend((
        "K", K, [cell for cell in cells if cell["K"] == K])
        for K in EXPECTED_K)
    bands = []
    for cell in cells:
        if cell["band"] not in bands:
            bands.append(cell["band"])
    groups.extend((
        "band", band, [cell for cell in cells if cell["band"] == band])
        for band in bands)

    aggregates = []
    for scope, key, selected in groups:
        phase_results = {}
        for phase_name, _ in PHASE_KEYS:
            statuses = [
                cell["phases"][phase_name]["status"] for cell in selected
            ]
            if "non_comparable_weak" in statuses:
                status = "non_comparable_weak"
                statistics: Mapping[str, Any] = {
                    "geometric_mean_ratio": None,
                    "lower_95_log_ratio": None,
                    "lower_95_ratio": None,
                    "mean_log_ratio": None,
                    "upper_95_log_ratio": None,
                    "upper_95_ratio": None,
                }
            elif "unresolved_zero_phase" in statuses:
                status = "unresolved_zero_phase"
                statistics = {
                    "geometric_mean_ratio": None,
                    "lower_95_log_ratio": None,
                    "lower_95_ratio": None,
                    "mean_log_ratio": None,
                    "upper_95_log_ratio": None,
                    "upper_95_ratio": None,
                }
            elif not selected or any(value != "comparable" for value in statuses):
                fail("phase aggregate has an unknown cell status")
            else:
                lane_means = [
                    math.fsum(
                        cell["phases"][phase_name]["pair_log_ratios"][lane]
                        for cell in selected) / len(selected)
                    for lane in range(PAIR_COUNT)
                ]
                status = "comparable"
                statistics = _confidence_interval(lane_means)
            phase_results[phase_name] = {
                "cell_count": len(selected),
                "independent_pair_count": PAIR_COUNT,
                "status": status,
                **statistics,
            }
        aggregates.append({
            "group_key": key,
            "profile": PROFILE_NAME,
            "profile_ordinal": PROFILE_ORDINAL,
            "scope": scope,
            "phases": phase_results,
        })
    threshold_log_ratio = math.log1p(0.02)
    global_width = [
        row["phases"]["outer"] for row in aggregates
        if row["scope"] in ("global", "width")
    ]
    local = [
        row["phases"]["outer"] for row in aggregates
        if row["scope"] in ("K", "band")
    ]
    if len(global_width) != 3 or not local:
        fail("phase decision aggregate geometry is incomplete")
    decision_resolved = all(
        row["status"] == "comparable" for row in global_width + local)
    upper_within = None if not decision_resolved else all(
        row["upper_95_log_ratio"] <= threshold_log_ratio
        for row in global_width)
    no_local_lower_above = None if not decision_resolved else all(
        row["lower_95_log_ratio"] <= threshold_log_ratio for row in local)
    sufficient = upper_within is True and no_local_lower_above is True
    decision = {
        "global_and_width_upper95_within_threshold": upper_within,
        "metric": "outer",
        "no_K_or_band_lower95_above_threshold": no_local_lower_above,
        "outcome": "n16_sufficient" if sufficient else "inconclusive",
        "permitted_inconclusive_next_step":
            "stop_or_one_separately_sealed_n24_screen",
        "profile": PROFILE_NAME,
        "profile_ordinal": PROFILE_ORDINAL,
        "ratio_direction": "left_over_right",
        "schema": PHASE_DECISION_SCHEMA,
        "threshold_log_ratio": threshold_log_ratio,
        "threshold_ratio": 1.02,
    }
    return {
        "aggregates": aggregates,
        "cells": cells,
        "decision": decision,
        "independent_unit":
            "counterbalanced-primary/opposite repeat pair",
        "left_arm": LEFT_ARM,
        "profile": PROFILE_NAME,
        "profile_ordinal": PROFILE_ORDINAL,
        "right_arm": RIGHT_ARM,
        "schema": PHASE_ANALYSIS_SCHEMA,
        "student_t_critical_df7": T95_DF7,
    }
