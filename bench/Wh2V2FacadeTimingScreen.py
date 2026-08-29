#!/usr/bin/env python3
"""Sealed controller and adjudicator for the public V2 facade timing falsifier.

The scientific run is deliberately a one-shot, role-separated experiment.  A
current worker exposes C/E/D and an exact-parent worker exposes P/L.  Both stay
alive, pinned to the same singleton CPU, while this controller interleaves the
individual public-API invocations in counterbalanced blocks.  The controller
does not time work; it only validates and sums the worker clocks.

This module is Python 3.8 compatible.  ``--selftest`` and ``--replay`` never
start a scientific workload.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
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
import time
import types
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


CAMPAIGN = "wh2-v2-facade-default-parent-falsifier-r0"
PARENT_COMMIT = "101584b7e5c30326b1791429221c331c82a00807"
CURRENT_IMPLEMENTATION_COMMIT = "d6ddb35e046956174202584fb5a26c9a79679ea8"
FROZEN_WORKER_SOURCE_SHA256 = \
    "fbafe52dbd677eee9c16bca5f8255f289f9e09bb62650a5372b12caa28d4d3c6"

WORKER_CONFIG_SCHEMA = "wirehair.wh2.v2-facade-timing-worker.config.v1"
WORKER_INVOCATION_SCHEMA = "wirehair.wh2.v2-facade-timing-worker.invocation.v1"
WORKER_TERMINAL_SCHEMA = "wirehair.wh2.v2-facade-timing-worker.terminal.v1"
SCREEN_CONFIG_SCHEMA = "wirehair.wh2.v2-facade-timing-screen.config.v1"
PANEL_SCHEMA = "wirehair.wh2.v2-facade-timing-screen.panel.v1"
RAW_TERMINAL_SCHEMA = "wirehair.wh2.v2-facade-timing-screen.terminal.v1"
SUMMARY_SCHEMA = "wirehair.wh2.v2-facade-timing-screen.summary.v1"
PROVENANCE_SCHEMA = "wirehair.wh2.v2-facade-timing-screen.provenance.v1"
BUILD_RECEIPT_SCHEMA = "wirehair.wh2.v2-facade-timing-build-receipt.v1"
COMPLETE_SCHEMA = "wirehair.wh2.v2-facade-timing-screen.complete.v1"

CELLS: Tuple[Tuple[int, int], ...] = tuple(
    (k, block_bytes)
    for k in (8, 128, 512, 5000, 64000)
    for block_bytes in (64, 1280)
)
SCOPES = (
    "prebuilt-systematic",
    "fresh-systematic",
    "fresh-repair",
    "prebuilt-repair",
)
SYSTEMATIC_SCOPES = frozenset(("prebuilt-systematic", "fresh-systematic"))
REPAIR_SCOPES = frozenset(("fresh-repair", "prebuilt-repair"))
ROLE_ARMS: Mapping[str, Tuple[str, ...]] = {
    "current": ("C", "E", "D"),
    "parent": ("P", "L"),
}
ARM_ROLE = {
    arm: role for role, arms in ROLE_ARMS.items() for arm in arms
}
COMPARISONS: Tuple[Tuple[str, str, str], ...] = (
    ("D/D", "D", "D"),
    ("E/E", "E", "E"),
    ("C/C", "C", "C"),
    ("P/P", "P", "P"),
    ("L/L", "L", "L"),
    ("E/D", "E", "D"),
    ("C/E", "C", "E"),
    ("D/P", "D", "P"),
    ("C/L", "C", "L"),
)
COMPARISON_ARMS = {name: (left, right) for name, left, right in COMPARISONS}
AA_COMPARISONS = frozenset(("D/D", "E/E", "C/C", "P/P", "L/L"))
EQUIVALENCE_COMPARISONS = frozenset(("E/D", "D/P"))
SYSTEMATIC_IMPROVEMENT_COMPARISONS = frozenset(("C/E", "C/L"))

REPLICATES = 12
INVOCATION_BUDGET = 65536
MINIMUM_INVOCATIONS = 24
INTERNAL_DEADLINE_SECONDS = 840
EXTERNAL_DEADLINE_SECONDS = 900
HEALTH_FINALIZATION_RESERVE_SECONDS = 20
T11_975 = 2.200985160082949
MAX_UINT64 = (1 << 64) - 1
MAX_INT63 = (1 << 63) - 1
MAX_JSON_LINE_BYTES = 64 * 1024
MAX_RAW_LINE_BYTES = 1024 * 1024
MAX_WORKER_STDOUT_BYTES = 4 * 1024 * 1024 * 1024
MAX_WORKER_STDERR_BYTES = 1024 * 1024
MAX_RAW_BYTES = 32 * 1024 * 1024
MAX_RECEIPT_BYTES = 8 * 1024 * 1024
MAX_FAILURES = 512
EXPECTED_PYTHON_IMAGE_NLINK = 2
EXPECTED_PANEL_COUNT = REPLICATES * len(CELLS) * len(SCOPES) * len(COMPARISONS)
EXPECTED_OUTPUT_NAMES = (
    "raw.jsonl",
    "summary.json",
    "provenance.json",
    "current.stderr",
    "parent.stderr",
)
HEALTH_ADAPTER_PATH = "bench/Wh2DirectSystematicComplementScreen.py"
SEALED_ENV = {
    "LANG": "C.UTF-8",
    "LC_ALL": "C.UTF-8",
    "PATH": "/usr/bin:/bin",
    "TZ": "UTC",
}
LOWER40 = re.compile(r"^[0-9a-f]{40}$")
LOWER64 = re.compile(r"^[0-9a-f]{64}$")
HEX_U64 = re.compile(r"^0x[0-9a-f]{16}$")


class ValidationError(RuntimeError):
    """The evidence is invalid, as distinct from a valid scientific reject."""


def fail(message: str) -> None:
    raise ValidationError(message)


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False,
    ).encode("ascii")


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def unique_object(pairs: Sequence[Tuple[str, Any]]) -> Dict[str, Any]:
    result: Dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            fail("duplicate JSON key: {}".format(key))
        result[key] = value
    return result


def reject_constant(value: str) -> None:
    fail("non-finite JSON number: {}".format(value))


def parse_canonical_line(
    data: bytes, where: str, maximum: int = MAX_JSON_LINE_BYTES,
) -> Dict[str, Any]:
    if (
        not data or len(data) > maximum or
        not data.endswith(b"\n") or data.count(b"\n") != 1 or b"\r" in data
    ):
        fail("{} is not one bounded LF-terminated line".format(where))
    try:
        value = json.loads(
            data[:-1].decode("ascii"), object_pairs_hook=unique_object,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        fail("{} is malformed: {}".format(where, exc))
    if type(value) is not dict or canonical_bytes(value) + b"\n" != data:
        fail("{} is not a canonical JSON object".format(where))
    return value


def exact_keys(value: Mapping[str, Any], expected: Iterable[str], where: str) -> None:
    wanted = set(expected)
    actual = set(value.keys())
    if actual != wanted:
        fail("{} keys differ: expected {}, got {}".format(
            where, sorted(wanted), sorted(actual)))


def exact_int(value: Any, lower: int, upper: int, where: str) -> int:
    if type(value) is not int or value < lower or value > upper:
        fail("{} is not an integer in [{},{}]".format(where, lower, upper))
    return value


def exact_string(value: Any, expected: str, where: str) -> str:
    if type(value) is not str or value != expected:
        fail("{} differs".format(where))
    return value


def bounded_string(value: Any, where: str, maximum: int = 4096) -> str:
    if type(value) is not str or not value or len(value) > maximum:
        fail("{} is not a bounded nonempty string".format(where))
    return value


def lower_hash(value: Any, where: str) -> str:
    if type(value) is not str or LOWER64.fullmatch(value) is None:
        fail("{} is not a lowercase SHA-256".format(where))
    return value


def lower_commit(value: Any, where: str) -> str:
    if type(value) is not str or LOWER40.fullmatch(value) is None:
        fail("{} is not a lowercase Git object ID".format(where))
    return value


def total_invocations(k: int) -> int:
    return max(MINIMUM_INVOCATIONS, (INVOCATION_BUDGET + k - 1) // k)


def invocations_by_replicate(k: int) -> List[int]:
    quotient, remainder = divmod(total_invocations(k), REPLICATES)
    return [
        quotient + (1 if replicate < remainder else 0)
        for replicate in range(REPLICATES)
    ]


def band_for(k: int) -> str:
    if k in (8, 128):
        return "small"
    if k in (512, 5000):
        return "medium"
    if k == 64000:
        return "large"
    raise AssertionError("unsealed K")


def panel_key_sha256(
    replicate: int, k: int, block_bytes: int, scope: str, comparison: str,
) -> str:
    return sha256_bytes(canonical_bytes({
        "K": k,
        "block_bytes": block_bytes,
        "comparison": comparison,
        "replicate": replicate,
        "scope": scope,
    }))


def panel_order(
    replicate: int, k: int, block_bytes: int, scope: str, comparison: str,
) -> str:
    # Replicate parity counterbalances the stable panel-specific phase bit.
    digest = panel_key_sha256(0, k, block_bytes, scope, comparison)
    phase = bytes.fromhex(digest)[-1] & 1
    return "ABBA" if (replicate & 1) == phase else "BAAB"


def sides_for(order: str) -> Tuple[str, ...]:
    if order == "ABBA":
        first = ("left", "right", "right", "left")
        second = ("right", "left", "left", "right")
    elif order == "BAAB":
        first = ("right", "left", "left", "right")
        second = ("left", "right", "right", "left")
    else:
        fail("unknown panel order")
        raise AssertionError
    return first + second


def lane_contrast(values: Sequence[int], order: str) -> float:
    if len(values) != 8 or any(type(value) is not int or value <= 0 for value in values):
        fail("lane contrast requires eight positive integer sums")
    logs = [math.log(value) for value in values]

    def contrast(offset: int, block_order: str) -> float:
        if block_order == "ABBA":
            return ((logs[offset] - logs[offset + 1]) +
                    (logs[offset + 3] - logs[offset + 2])) / 2.0
        return ((logs[offset + 1] - logs[offset]) +
                (logs[offset + 2] - logs[offset + 3])) / 2.0

    opposite = "BAAB" if order == "ABBA" else "ABBA"
    return 0.5 * (contrast(0, order) + contrast(4, opposite))


def sample_summary(values: Sequence[float]) -> Dict[str, Any]:
    if len(values) != REPLICATES or any(
        type(value) not in (int, float) or not math.isfinite(value)
        for value in values
    ):
        fail("statistical sample is incomplete or non-finite")
    mean = math.fsum(values) / REPLICATES
    variance = math.fsum((value - mean) ** 2 for value in values) / (REPLICATES - 1)
    if variance < 0.0 or not math.isfinite(variance):
        fail("statistical variance is invalid")
    standard_error = math.sqrt(variance / REPLICATES)
    lower_log = mean - T11_975 * standard_error
    upper_log = mean + T11_975 * standard_error
    return {
        "geometric_mean": math.exp(mean),
        "log_mean": mean,
        "log_standard_error": standard_error,
        "lower95": math.exp(lower_log),
        "lower95_log": lower_log,
        "n": REPLICATES,
        "upper95": math.exp(upper_log),
        "upper95_log": upper_log,
    }


WORKER_CONFIG_KEYS = {
    "arm_descriptors", "campaign", "cells", "implementation_git_commit",
    "internal_deadline_seconds", "invocation_budget", "minimum_invocations",
    "panel_replicates", "role", "schema", "supported_arms", "target_cpu",
    "worker_git_commit",
}
WORKER_CELL_KEYS = {
    "K", "block_bytes", "invocations_by_replicate", "message_bytes",
    "oracles", "source_seed", "source_sha256",
}
ORACLE_KEYS = {
    "arm", "arm_descriptor_sha256", "borrowed_systematic_ids",
    "equation_configuration_sha256", "first_repair_sha256", "high_id_sha256",
    "public_state_receipt_sha256", "repair_sha256",
    "roundtrip_sha256", "roundtrip_verified",
    "serialized_profile_bytes", "serialized_profile_hex",
    "serialized_profile_sha256",
    "systematic_sha256",
}
ARM_DESCRIPTOR_ENTRY_KEYS = {"arm", "descriptor", "descriptor_sha256"}
ARM_DESCRIPTOR_KEYS = {
    "api", "arm", "codec", "equation_profile", "schema", "source_policy",
}
ARM_DESCRIPTOR_SCHEMA = "wirehair.wh2.v2-facade-timing-worker.arm.v1"
ARM_DESCRIPTOR_VALUES: Mapping[str, Mapping[str, str]] = {
    "C": {
        "api": "wirehair_v2_encoder_create_with_options",
        "codec": "wirehair2",
        "equation_profile": "certified-2026-07",
        "source_policy": "borrowed-immutable",
    },
    "E": {
        "api": "wirehair_v2_encoder_create_with_options",
        "codec": "wirehair2",
        "equation_profile": "certified-2026-07",
        "source_policy": "independent",
    },
    "D": {
        "api": "wirehair_v2_encoder_create",
        "codec": "wirehair2",
        "equation_profile": "certified-2026-07",
        "source_policy": "default-independent",
    },
    "P": {
        "api": "wirehair_v2_encoder_create",
        "codec": "wirehair2",
        "equation_profile": "certified-2026-07",
        "source_policy": "pre-feature-independent",
    },
    "L": {
        "api": "wirehair_encoder_create_ex",
        "codec": "wirehair1",
        "equation_profile": "legacy-current",
        "source_policy": "borrowed",
    },
}
INVOCATION_KEYS = {
    "K", "arm", "block_bytes", "constructor_ns", "correctness_sha256",
    "descriptor_sha256", "borrowed_systematic_invocations", "elapsed_ns",
    "init_first_repair_ns", "result", "role", "schema", "scope", "seq",
    "systematic_ns", "target_cpu",
}
TERMINAL_KEYS = {"invocation_count", "role", "schema", "status", "target_cpu"}


@dataclass(frozen=True)
class WorkerContract:
    role: str
    cpu: int
    implementation_commit: str
    worker_commit: str
    config: Dict[str, Any]
    oracles: Mapping[Tuple[int, int, str], Dict[str, Any]]


def validate_serialized_profile(oracle: Mapping[str, Any], arm: str, where: str) -> None:
    profile_bytes = oracle["serialized_profile_bytes"]
    profile = oracle["serialized_profile_hex"]
    digest = oracle["serialized_profile_sha256"]
    if arm == "L":
        if profile_bytes != 0 or profile is not None or digest is not None:
            fail("{} WH1 profile fields must be null".format(where))
        return
    exact_int(profile_bytes, 32, 32, where + ".serialized_profile_bytes")
    if type(profile) is not str or len(profile) != 64:
        fail("{} serialized profile is not bounded lowercase hex".format(where))
    if len(profile) & 1 or re.fullmatch(r"[0-9a-f]+", profile) is None:
        fail("{} serialized profile is not lowercase hex".format(where))
    expected = lower_hash(digest, "{}.serialized_profile_sha256".format(where))
    if sha256_bytes(bytes.fromhex(profile)) != expected:
        fail("{} serialized profile hash differs".format(where))


def validate_worker_config(
    config: Dict[str, Any], role: str, cpu: int, worker_commit: str,
    implementation_commit: str,
) -> WorkerContract:
    where = "{} worker config".format(role)
    exact_keys(config, WORKER_CONFIG_KEYS, where)
    exact_string(config["schema"], WORKER_CONFIG_SCHEMA, where + ".schema")
    exact_string(config["campaign"], CAMPAIGN, where + ".campaign")
    exact_string(config["role"], role, where + ".role")
    exact_int(config["target_cpu"], cpu, cpu, where + ".target_cpu")
    exact_int(
        config["internal_deadline_seconds"], INTERNAL_DEADLINE_SECONDS,
        INTERNAL_DEADLINE_SECONDS, where + ".internal_deadline_seconds")
    exact_int(config["invocation_budget"], INVOCATION_BUDGET, INVOCATION_BUDGET,
              where + ".invocation_budget")
    exact_int(config["minimum_invocations"], MINIMUM_INVOCATIONS,
              MINIMUM_INVOCATIONS, where + ".minimum_invocations")
    exact_int(config["panel_replicates"], REPLICATES, REPLICATES,
              where + ".panel_replicates")
    exact_string(config["worker_git_commit"], worker_commit,
                 where + ".worker_git_commit")
    exact_string(config["implementation_git_commit"], implementation_commit,
                 where + ".implementation_git_commit")
    if config["supported_arms"] != list(ROLE_ARMS[role]):
        fail("{} supported arm roster differs".format(where))
    descriptors = config["arm_descriptors"]
    if type(descriptors) is not list or len(descriptors) != len(ROLE_ARMS[role]):
        fail("{} arm descriptor roster differs".format(where))
    descriptor_hashes: Dict[str, str] = {}
    for descriptor_index, (arm, entry) in enumerate(
        zip(ROLE_ARMS[role], descriptors)
    ):
        descriptor_where = "{}.arm_descriptors[{}]".format(
            where, descriptor_index)
        if type(entry) is not dict:
            fail("{} is not an object".format(descriptor_where))
        exact_keys(entry, ARM_DESCRIPTOR_ENTRY_KEYS, descriptor_where)
        exact_string(entry["arm"], arm, descriptor_where + ".arm")
        descriptor = entry["descriptor"]
        if type(descriptor) is not dict:
            fail("{}.descriptor is not an object".format(descriptor_where))
        exact_keys(descriptor, ARM_DESCRIPTOR_KEYS, descriptor_where + ".descriptor")
        exact_string(descriptor["schema"], ARM_DESCRIPTOR_SCHEMA,
                     descriptor_where + ".descriptor.schema")
        exact_string(descriptor["arm"], arm,
                     descriptor_where + ".descriptor.arm")
        for field_name, expected_value in ARM_DESCRIPTOR_VALUES[arm].items():
            exact_string(descriptor[field_name], expected_value,
                         descriptor_where + ".descriptor." + field_name)
        digest = lower_hash(entry["descriptor_sha256"],
                            descriptor_where + ".descriptor_sha256")
        if sha256_bytes(canonical_bytes(descriptor)) != digest:
            fail("{} descriptor hash differs".format(descriptor_where))
        descriptor_hashes[arm] = digest

    cells = config["cells"]
    if type(cells) is not list or len(cells) != len(CELLS):
        fail("{} cell roster differs".format(where))
    oracles: Dict[Tuple[int, int, str], Dict[str, Any]] = {}
    for cell_index, ((k, block_bytes), cell) in enumerate(zip(CELLS, cells)):
        cell_where = "{}.cells[{}]".format(where, cell_index)
        if type(cell) is not dict:
            fail("{} is not an object".format(cell_where))
        exact_keys(cell, WORKER_CELL_KEYS, cell_where)
        exact_int(cell["K"], k, k, cell_where + ".K")
        exact_int(cell["block_bytes"], block_bytes, block_bytes,
                  cell_where + ".block_bytes")
        exact_int(cell["message_bytes"], k * block_bytes, k * block_bytes,
                  cell_where + ".message_bytes")
        expected_n = invocations_by_replicate(k)
        if cell["invocations_by_replicate"] != expected_n:
            fail("{} invocation allocation differs".format(cell_where))
        seed = cell["source_seed"]
        if type(seed) is not str or HEX_U64.fullmatch(seed) is None:
            fail("{} source seed is not canonical hex".format(cell_where))
        lower_hash(cell["source_sha256"], cell_where + ".source_sha256")
        arm_oracles = cell["oracles"]
        if type(arm_oracles) is not list or len(arm_oracles) != len(ROLE_ARMS[role]):
            fail("{} oracle roster differs".format(cell_where))
        for oracle_index, (arm, oracle) in enumerate(zip(ROLE_ARMS[role], arm_oracles)):
            oracle_where = "{}.oracles[{}]".format(cell_where, oracle_index)
            if type(oracle) is not dict:
                fail("{} is not an object".format(oracle_where))
            exact_keys(oracle, ORACLE_KEYS, oracle_where)
            exact_string(oracle["arm"], arm, oracle_where + ".arm")
            for field in ORACLE_KEYS - {
                "arm", "borrowed_systematic_ids", "roundtrip_verified",
                "serialized_profile_bytes", "serialized_profile_hex",
                "serialized_profile_sha256",
            }:
                lower_hash(oracle[field], oracle_where + "." + field)
            expected_hits = k if arm == "C" else 0
            exact_int(oracle["borrowed_systematic_ids"], expected_hits,
                      expected_hits, oracle_where + ".borrowed_systematic_ids")
            if oracle["roundtrip_verified"] is not True:
                fail("{} lacks public decode roundtrip".format(oracle_where))
            exact_string(oracle["roundtrip_sha256"], cell["source_sha256"],
                         oracle_where + ".roundtrip_sha256")
            exact_string(oracle["systematic_sha256"], cell["source_sha256"],
                         oracle_where + ".systematic_sha256")
            validate_serialized_profile(oracle, arm, oracle_where)
            exact_string(oracle["arm_descriptor_sha256"], descriptor_hashes[arm],
                         oracle_where + ".arm_descriptor_sha256")
            oracles[(k, block_bytes, arm)] = dict(oracle)
    return WorkerContract(
        role=role, cpu=cpu, implementation_commit=implementation_commit,
        worker_commit=worker_commit, config=config, oracles=oracles,
    )


def validate_cross_role_contracts(
    current: WorkerContract, parent: WorkerContract,
) -> None:
    if current.cpu != parent.cpu or current.worker_commit != parent.worker_commit:
        fail("worker roles do not share the exact harness and CPU")
    v2_arms = ("C", "E", "D", "P")
    v2_fields = (
        "equation_configuration_sha256", "first_repair_sha256",
        "high_id_sha256", "public_state_receipt_sha256", "repair_sha256",
        "serialized_profile_bytes", "serialized_profile_hex",
        "serialized_profile_sha256",
        "systematic_sha256",
    )
    for cell_index, (k, block_bytes) in enumerate(CELLS):
        current_cell = current.config["cells"][cell_index]
        parent_cell = parent.config["cells"][cell_index]
        if (
            current_cell["source_seed"] != parent_cell["source_seed"] or
            current_cell["source_sha256"] != parent_cell["source_sha256"] or
            current_cell["message_bytes"] != parent_cell["message_bytes"] or
            current_cell["invocations_by_replicate"] !=
                parent_cell["invocations_by_replicate"]
        ):
            fail("K{} B{} source/work allocation differs across roles".format(
                k, block_bytes))
        reference = current.oracles[(k, block_bytes, "C")]
        for arm in v2_arms:
            contract = current if arm in ROLE_ARMS["current"] else parent
            oracle = contract.oracles[(k, block_bytes, arm)]
            for field in v2_fields:
                if oracle[field] != reference[field]:
                    fail("K{} B{} WH2 {} differs at {}".format(
                        k, block_bytes, arm, field))
        # WH1 is intentionally a different equation system.  Only the source
        # systematic bytes are cross-bound; its repair/high-ID receipts remain
        # exact, stable L-local evidence.
        l_oracle = parent.oracles[(k, block_bytes, "L")]
        if l_oracle["systematic_sha256"] != reference["systematic_sha256"]:
            fail("K{} B{} WH1 source systematic oracle differs".format(
                k, block_bytes))
        descriptors = [
            (current if arm in ROLE_ARMS["current"] else parent).oracles[
                (k, block_bytes, arm)]["arm_descriptor_sha256"]
            for arm in ("C", "E", "D", "P", "L")
        ]
        if len(set(descriptors)) != len(descriptors):
            fail("K{} B{} arm descriptors are not unique".format(k, block_bytes))


def expected_correctness_hash(
    contract: WorkerContract, k: int, block_bytes: int, arm: str, scope: str,
) -> str:
    oracle = contract.oracles[(k, block_bytes, arm)]
    return oracle[
        "systematic_sha256" if scope in SYSTEMATIC_SCOPES else "repair_sha256"
    ]


def validate_invocation(
    value: Dict[str, Any], contract: WorkerContract, seq: int, arm: str,
    scope: str, k: int, block_bytes: int,
) -> Dict[str, Any]:
    where = "{} invocation {}".format(contract.role, seq)
    exact_keys(value, INVOCATION_KEYS, where)
    exact_string(value["schema"], WORKER_INVOCATION_SCHEMA, where + ".schema")
    exact_int(value["seq"], seq, seq, where + ".seq")
    exact_string(value["role"], contract.role, where + ".role")
    exact_string(value["arm"], arm, where + ".arm")
    exact_string(value["scope"], scope, where + ".scope")
    exact_int(value["K"], k, k, where + ".K")
    exact_int(value["block_bytes"], block_bytes, block_bytes,
              where + ".block_bytes")
    exact_int(value["target_cpu"], contract.cpu, contract.cpu,
              where + ".target_cpu")
    exact_int(value["result"], 0, 0, where + ".result")
    elapsed = exact_int(value["elapsed_ns"], 1, MAX_INT63, where + ".elapsed_ns")
    oracle = contract.oracles[(k, block_bytes, arm)]
    exact_string(value["descriptor_sha256"], oracle["arm_descriptor_sha256"],
                 where + ".descriptor_sha256")
    exact_string(
        value["correctness_sha256"],
        expected_correctness_hash(contract, k, block_bytes, arm, scope),
        where + ".correctness_sha256")
    expected_hits = k if arm == "C" and scope in SYSTEMATIC_SCOPES else 0
    exact_int(value["borrowed_systematic_invocations"], expected_hits,
              expected_hits, where + ".borrowed_systematic_invocations")
    if scope == "fresh-systematic":
        constructor = exact_int(value["constructor_ns"], 1, elapsed,
                                where + ".constructor_ns")
        systematic = exact_int(value["systematic_ns"], 1, elapsed,
                               where + ".systematic_ns")
        if constructor + systematic > elapsed:
            fail("{} nested clocks exceed total".format(where))
        if value["init_first_repair_ns"] is not None:
            fail("{} has repair clock in systematic scope".format(where))
    elif scope == "fresh-repair":
        exact_int(value["init_first_repair_ns"], 1, elapsed,
                  where + ".init_first_repair_ns")
        if value["constructor_ns"] is not None or value["systematic_ns"] is not None:
            fail("{} has systematic clocks in repair scope".format(where))
    elif (
        value["constructor_ns"] is not None or
        value["systematic_ns"] is not None or
        value["init_first_repair_ns"] is not None
    ):
        fail("{} has nested clocks in a prebuilt scope".format(where))
    return value


def validate_worker_terminal(
    value: Dict[str, Any], contract: WorkerContract, invocation_count: int,
) -> None:
    where = "{} worker terminal".format(contract.role)
    exact_keys(value, TERMINAL_KEYS, where)
    exact_string(value["schema"], WORKER_TERMINAL_SCHEMA, where + ".schema")
    exact_string(value["role"], contract.role, where + ".role")
    exact_string(value["status"], "complete", where + ".status")
    exact_int(value["target_cpu"], contract.cpu, contract.cpu, where + ".target_cpu")
    exact_int(value["invocation_count"], invocation_count, invocation_count,
              where + ".invocation_count")


SCREEN_CONFIG_KEYS = {
    "campaign", "cells", "comparisons", "current_worker",
    "expected_panel_count", "internal_deadline_seconds", "invocation_budget",
    "minimum_invocations", "panel_replicates", "parent_worker", "schema",
    "scopes", "target_cpu",
}
PANEL_KEYS = {
    "K", "block_bytes", "comparison", "first_sequence",
    "invocations_per_slot", "last_sequence", "left_arm", "order",
    "panel_key_sha256", "replicate", "right_arm", "schema", "scope",
    "slots", "warmup_response_stream_sha256", "warmup_sequences",
}
SLOT_KEYS = {
    "arm", "borrowed_systematic_invocations", "constructor_ns", "elapsed_ns",
    "first_sequence", "init_first_repair_ns", "invocation_count",
    "last_sequence", "response_stream_sha256", "side", "systematic_ns",
}
RAW_TERMINAL_KEYS = {
    "current_invocations", "panel_count", "parent_invocations",
    "raw_config_sha256", "schema", "status",
}


def make_screen_config(
    current: WorkerContract, parent: WorkerContract,
) -> Dict[str, Any]:
    return {
        "campaign": CAMPAIGN,
        "cells": [{"K": k, "block_bytes": block_bytes}
                  for k, block_bytes in CELLS],
        "comparisons": [
            {"comparison": name, "left_arm": left, "right_arm": right}
            for name, left, right in COMPARISONS
        ],
        "current_worker": current.config,
        "expected_panel_count": EXPECTED_PANEL_COUNT,
        "internal_deadline_seconds": INTERNAL_DEADLINE_SECONDS,
        "invocation_budget": INVOCATION_BUDGET,
        "minimum_invocations": MINIMUM_INVOCATIONS,
        "panel_replicates": REPLICATES,
        "parent_worker": parent.config,
        "schema": SCREEN_CONFIG_SCHEMA,
        "scopes": list(SCOPES),
        "target_cpu": current.cpu,
    }


def validate_screen_config(
    value: Dict[str, Any], cpu: int, worker_commit: str,
    current_commit: str,
) -> Tuple[WorkerContract, WorkerContract]:
    exact_keys(value, SCREEN_CONFIG_KEYS, "screen config")
    exact_string(value["schema"], SCREEN_CONFIG_SCHEMA, "screen config.schema")
    exact_string(value["campaign"], CAMPAIGN, "screen config.campaign")
    exact_int(value["target_cpu"], cpu, cpu, "screen config.target_cpu")
    exact_int(value["internal_deadline_seconds"], INTERNAL_DEADLINE_SECONDS,
              INTERNAL_DEADLINE_SECONDS, "screen config deadline")
    exact_int(value["invocation_budget"], INVOCATION_BUDGET, INVOCATION_BUDGET,
              "screen config invocation budget")
    exact_int(value["minimum_invocations"], MINIMUM_INVOCATIONS,
              MINIMUM_INVOCATIONS, "screen config minimum invocations")
    exact_int(value["panel_replicates"], REPLICATES, REPLICATES,
              "screen config replicates")
    exact_int(value["expected_panel_count"], EXPECTED_PANEL_COUNT,
              EXPECTED_PANEL_COUNT, "screen config panel count")
    expected_cells = [{"K": k, "block_bytes": block_bytes}
                      for k, block_bytes in CELLS]
    if value["cells"] != expected_cells or value["scopes"] != list(SCOPES):
        fail("screen config cell/scope roster differs")
    expected_comparisons = [
        {"comparison": name, "left_arm": left, "right_arm": right}
        for name, left, right in COMPARISONS
    ]
    if value["comparisons"] != expected_comparisons:
        fail("screen config comparison roster differs")
    current_value = value["current_worker"]
    parent_value = value["parent_worker"]
    if type(current_value) is not dict or type(parent_value) is not dict:
        fail("screen config worker contracts are not objects")
    current = validate_worker_config(
        current_value, "current", cpu, worker_commit, current_commit)
    parent = validate_worker_config(
        parent_value, "parent", cpu, worker_commit, PARENT_COMMIT)
    validate_cross_role_contracts(current, parent)
    return current, parent


def expected_panel_slot_sequence(
    first_sequence: int, n: int, slot_index: int,
) -> Tuple[int, int]:
    if slot_index < 4:
        first = first_sequence + 2 + slot_index
    else:
        first = first_sequence + 2 + 4 * n + (slot_index - 4)
    return first, first + 4 * (n - 1)


def validate_panel(
    panel: Dict[str, Any], replicate: int, k: int, block_bytes: int,
    scope: str, comparison: str, next_sequence: int,
) -> Tuple[float, int]:
    where = "panel r{} K{} B{} {} {}".format(
        replicate, k, block_bytes, scope, comparison)
    exact_keys(panel, PANEL_KEYS, where)
    exact_string(panel["schema"], PANEL_SCHEMA, where + ".schema")
    exact_int(panel["K"], k, k, where + ".K")
    exact_int(panel["block_bytes"], block_bytes, block_bytes,
              where + ".block_bytes")
    exact_string(panel["scope"], scope, where + ".scope")
    exact_string(panel["comparison"], comparison, where + ".comparison")
    exact_int(panel["replicate"], replicate, replicate, where + ".replicate")
    left_arm, right_arm = COMPARISON_ARMS[comparison]
    exact_string(panel["left_arm"], left_arm, where + ".left_arm")
    exact_string(panel["right_arm"], right_arm, where + ".right_arm")
    expected_key = panel_key_sha256(
        replicate, k, block_bytes, scope, comparison)
    exact_string(panel["panel_key_sha256"], expected_key,
                 where + ".panel_key_sha256")
    order = panel_order(replicate, k, block_bytes, scope, comparison)
    exact_string(panel["order"], order, where + ".order")
    n = invocations_by_replicate(k)[replicate]
    exact_int(panel["invocations_per_slot"], n, n,
              where + ".invocations_per_slot")
    first_sequence = exact_int(panel["first_sequence"], next_sequence,
                               next_sequence, where + ".first_sequence")
    last_sequence = first_sequence + 2 + 8 * n - 1
    exact_int(panel["last_sequence"], last_sequence, last_sequence,
              where + ".last_sequence")
    if panel["warmup_sequences"] != [first_sequence, first_sequence + 1]:
        fail("{} warmup sequence roster differs".format(where))
    lower_hash(panel["warmup_response_stream_sha256"],
               where + ".warmup_response_stream_sha256")
    slots = panel["slots"]
    sides = sides_for(order)
    if type(slots) is not list or len(slots) != 8:
        fail("{} does not contain eight slot sums".format(where))
    elapsed: List[int] = []
    for slot_index, (slot, side) in enumerate(zip(slots, sides)):
        slot_where = "{}.slots[{}]".format(where, slot_index)
        if type(slot) is not dict:
            fail("{} is not an object".format(slot_where))
        exact_keys(slot, SLOT_KEYS, slot_where)
        exact_string(slot["side"], side, slot_where + ".side")
        arm = left_arm if side == "left" else right_arm
        exact_string(slot["arm"], arm, slot_where + ".arm")
        exact_int(slot["invocation_count"], n, n,
                  slot_where + ".invocation_count")
        sequence_first, sequence_last = expected_panel_slot_sequence(
            first_sequence, n, slot_index)
        exact_int(slot["first_sequence"], sequence_first, sequence_first,
                  slot_where + ".first_sequence")
        exact_int(slot["last_sequence"], sequence_last, sequence_last,
                  slot_where + ".last_sequence")
        lower_hash(slot["response_stream_sha256"],
                   slot_where + ".response_stream_sha256")
        total = exact_int(slot["elapsed_ns"], n, MAX_INT63,
                          slot_where + ".elapsed_ns")
        elapsed.append(total)
        expected_hits = n * k if arm == "C" and scope in SYSTEMATIC_SCOPES else 0
        exact_int(slot["borrowed_systematic_invocations"], expected_hits,
                  expected_hits,
                  slot_where + ".borrowed_systematic_invocations")
        if scope == "fresh-systematic":
            constructor = exact_int(slot["constructor_ns"], n, total,
                                    slot_where + ".constructor_ns")
            systematic = exact_int(slot["systematic_ns"], n, total,
                                   slot_where + ".systematic_ns")
            if constructor + systematic > total:
                fail("{} nested sums exceed total".format(slot_where))
            exact_int(slot["init_first_repair_ns"], 0, 0,
                      slot_where + ".init_first_repair_ns")
        elif scope == "fresh-repair":
            exact_int(slot["constructor_ns"], 0, 0,
                      slot_where + ".constructor_ns")
            exact_int(slot["systematic_ns"], 0, 0,
                      slot_where + ".systematic_ns")
            exact_int(slot["init_first_repair_ns"], n, total,
                      slot_where + ".init_first_repair_ns")
        else:
            for field_name in (
                "constructor_ns", "systematic_ns", "init_first_repair_ns"
            ):
                exact_int(slot[field_name], 0, 0,
                          slot_where + "." + field_name)
    return lane_contrast(elapsed, order), last_sequence + 1


def aggregate_logs(
    logs: Mapping[Tuple[str, int, int, str], List[float]], comparison: str,
    scope: str, cells: Sequence[Tuple[int, int]],
) -> List[float]:
    return [
        math.fsum(logs[(comparison, k, block_bytes, scope)][replicate]
                  for k, block_bytes in cells) / len(cells)
        for replicate in range(REPLICATES)
    ]


def adjudicate(
    logs: Mapping[Tuple[str, int, int, str], List[float]],
) -> Tuple[Dict[str, Any], List[str]]:
    failures: List[str] = []
    cell_stats: Dict[str, Any] = {}
    equivalence_limit_log = math.log1p(0.02)
    improvement_point_log = math.log(1.02)
    improvement_group_log = math.log(0.99)
    repair_limit_log = math.log(1.02)

    for comparison, _, _ in COMPARISONS:
        for k, block_bytes in CELLS:
            for scope in SCOPES:
                key = (comparison, k, block_bytes, scope)
                label = "{}:{}:K{}:B{}".format(
                    comparison, scope, k, block_bytes)
                summary = sample_summary(logs[key])
                cell_stats[label] = summary
                if comparison in AA_COMPARISONS:
                    if not (
                        summary["lower95_log"] > -equivalence_limit_log and
                        summary["upper95_log"] < equivalence_limit_log
                    ):
                        failures.append("aa_ci:" + label)
                elif comparison in EQUIVALENCE_COMPARISONS:
                    if not (
                        summary["lower95_log"] > -equivalence_limit_log and
                        summary["upper95_log"] < equivalence_limit_log
                    ):
                        failures.append("equivalence_ci:" + label)
                elif (
                    comparison in SYSTEMATIC_IMPROVEMENT_COMPARISONS and
                    scope in SYSTEMATIC_SCOPES
                ):
                    if summary["log_mean"] > improvement_point_log:
                        failures.append("systematic_cell_point:" + label)
                elif comparison == "C/E" and scope in REPAIR_SCOPES:
                    if summary["log_mean"] > repair_limit_log:
                        failures.append("repair_cell_point:" + label)
                    if summary["upper95_log"] > repair_limit_log:
                        failures.append("repair_cell_upper95:" + label)

    group_roster: List[Tuple[str, Tuple[Tuple[int, int], ...]]] = [
        ("equal-cell", CELLS),
        ("width-B64", tuple(cell for cell in CELLS if cell[1] == 64)),
        ("width-B1280", tuple(cell for cell in CELLS if cell[1] == 1280)),
    ]
    group_roster.extend(
        ("band-" + band,
         tuple(cell for cell in CELLS if band_for(cell[0]) == band))
        for band in ("small", "medium", "large")
    )
    aggregate_stats: Dict[str, Any] = {}
    for comparison in ("C/E", "C/L"):
        comparison_groups: Dict[str, Any] = {}
        for scope in SCOPES:
            scope_groups: Dict[str, Any] = {}
            for group_name, cells in group_roster:
                summary = sample_summary(
                    aggregate_logs(logs, comparison, scope, cells))
                scope_groups[group_name] = summary
                if scope in SYSTEMATIC_SCOPES:
                    if not summary["upper95_log"] < improvement_group_log:
                        failures.append("systematic_group_upper95:{}:{}:{}".format(
                            comparison, scope, group_name))
                elif comparison == "C/E":
                    if summary["upper95_log"] > repair_limit_log:
                        failures.append("repair_group_upper95:{}:{}:{}".format(
                            comparison, scope, group_name))
            comparison_groups[scope] = scope_groups
        aggregate_stats[comparison] = comparison_groups
    if len(failures) > MAX_FAILURES:
        fail("failed-gate roster exceeds its hard cap")
    statistics = {
        "aggregates": aggregate_stats,
        "cells": cell_stats,
        "failed_gates": failures,
        "gate_policy": {
            "aa_and_equivalence_two_sided_ratio": 1.02,
            "repair_upper_ratio_inclusive": 1.02,
            "systematic_group_upper_ratio_strict": 0.99,
            "systematic_point_upper_ratio_inclusive": 1.02,
            "wh1_repair_is_descriptive_only": True,
        },
        "student_t_975_df11": T11_975,
    }
    return statistics, failures


def classify_outcome(
    infrastructure: Sequence[str], raw_complete: bool,
    gate_failures: Sequence[str],
) -> str:
    """Apply the preregistered invalid-versus-reject disposition law."""
    if (
        infrastructure or not raw_complete
        or any(item.startswith("aa_ci:") for item in gate_failures)
    ):
        return "invalid"
    return "reject" if gate_failures else "pass"


def parse_raw_transcript(
    raw: bytes, cpu: int, worker_commit: str, current_commit: str,
) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any], List[str]]:
    if (
        not raw or len(raw) > MAX_RAW_BYTES or not raw.endswith(b"\n") or
        b"\r" in raw
    ):
        fail("raw transcript is not bounded LF-terminated JSONL")
    lines = raw.splitlines(keepends=True)
    if len(lines) != EXPECTED_PANEL_COUNT + 2:
        fail("raw transcript record count differs")
    records: List[Dict[str, Any]] = []
    for index, line in enumerate(lines):
        if len(line) > MAX_RAW_LINE_BYTES:
            fail("raw transcript line {} exceeds its cap".format(index))
        records.append(parse_canonical_line(
            line, "raw line {}".format(index), MAX_RAW_LINE_BYTES))
    config = records[0]
    validate_screen_config(config, cpu, worker_commit, current_commit)
    logs: Dict[Tuple[str, int, int, str], List[float]] = {}
    next_sequence = 0
    record_index = 1
    for replicate in range(REPLICATES):
        for k, block_bytes in CELLS:
            for scope in SCOPES:
                for comparison, _, _ in COMPARISONS:
                    value, next_sequence = validate_panel(
                        records[record_index], replicate, k, block_bytes,
                        scope, comparison, next_sequence)
                    logs.setdefault(
                        (comparison, k, block_bytes, scope), []).append(value)
                    record_index += 1
    terminal = records[record_index]
    exact_keys(terminal, RAW_TERMINAL_KEYS, "raw terminal")
    exact_string(terminal["schema"], RAW_TERMINAL_SCHEMA,
                 "raw terminal.schema")
    exact_string(terminal["status"], "complete", "raw terminal.status")
    exact_int(terminal["panel_count"], EXPECTED_PANEL_COUNT,
              EXPECTED_PANEL_COUNT, "raw terminal.panel_count")
    exact_string(terminal["raw_config_sha256"], sha256_bytes(canonical_bytes(config)),
                 "raw terminal.raw_config_sha256")
    expected_role_counts = {"current": 0, "parent": 0}
    # Two warmups plus eight timed lanes per panel, attributed by arm.
    for replicate in range(REPLICATES):
        for k, _ in CELLS:
            n = invocations_by_replicate(k)[replicate]
            for _scope in SCOPES:
                for _comparison, left, right in COMPARISONS:
                    expected_role_counts[ARM_ROLE[left]] += 1 + 4 * n
                    expected_role_counts[ARM_ROLE[right]] += 1 + 4 * n
    exact_int(terminal["current_invocations"], expected_role_counts["current"],
              expected_role_counts["current"], "raw terminal current count")
    exact_int(terminal["parent_invocations"], expected_role_counts["parent"],
              expected_role_counts["parent"], "raw terminal parent count")
    if sum(expected_role_counts.values()) != next_sequence:
        fail("raw sequence count and role counts differ")
    statistics, failures = adjudicate(logs)
    return config, terminal, statistics, failures


def _synthetic_hash(label: str) -> str:
    return sha256_bytes(label.encode("ascii"))


def synthetic_worker_config(
    role: str, cpu: int, worker_commit: str, implementation_commit: str,
) -> Dict[str, Any]:
    cells: List[Dict[str, Any]] = []
    for k, block_bytes in CELLS:
        cell_label = "K{}-B{}".format(k, block_bytes)
        source = _synthetic_hash("source-" + cell_label)
        v2_profile = _synthetic_hash("v2-profile-" + cell_label)
        v2_fields = {
            "equation_configuration_sha256": _synthetic_hash(
                "v2-equation-" + cell_label),
            "first_repair_sha256": _synthetic_hash(
                "v2-first-repair-" + cell_label),
            "high_id_sha256": _synthetic_hash("v2-high-" + cell_label),
            "public_state_receipt_sha256": _synthetic_hash(
                "v2-state-" + cell_label),
            "repair_sha256": _synthetic_hash("v2-repair-" + cell_label),
            "serialized_profile_bytes": 32,
            "serialized_profile_hex": v2_profile,
            "serialized_profile_sha256": sha256_bytes(bytes.fromhex(v2_profile)),
            "systematic_sha256": source,
        }
        oracles: List[Dict[str, Any]] = []
        for arm in ROLE_ARMS[role]:
            fields = dict(v2_fields)
            if arm == "L":
                fields.update({
                    "equation_configuration_sha256": _synthetic_hash(
                        "l-equation-" + cell_label),
                    "first_repair_sha256": _synthetic_hash(
                        "l-first-repair-" + cell_label),
                    "high_id_sha256": _synthetic_hash("l-high-" + cell_label),
                    "public_state_receipt_sha256": _synthetic_hash(
                        "l-state-" + cell_label),
                    "repair_sha256": _synthetic_hash("l-repair-" + cell_label),
                    "serialized_profile_bytes": 0,
                    "serialized_profile_hex": None,
                    "serialized_profile_sha256": None,
                })
            oracle = {
                "arm": arm,
                "arm_descriptor_sha256": _synthetic_hash(
                    "descriptor-{}-{}".format(arm, cell_label)),
                "borrowed_systematic_ids": k if arm == "C" else 0,
                "roundtrip_sha256": source,
                "roundtrip_verified": True,
            }
            oracle.update(fields)
            oracles.append(oracle)
        cells.append({
            "K": k,
            "block_bytes": block_bytes,
            "invocations_by_replicate": invocations_by_replicate(k),
            "message_bytes": k * block_bytes,
            "oracles": oracles,
            "source_seed": "0x{:016x}".format(
                0x5470000000000000 ^ (k << 16) ^ block_bytes),
            "source_sha256": source,
        })
    arm_descriptors: List[Dict[str, Any]] = []
    for arm in ROLE_ARMS[role]:
        descriptor = {
            "api": ARM_DESCRIPTOR_VALUES[arm]["api"],
            "arm": arm,
            "codec": ARM_DESCRIPTOR_VALUES[arm]["codec"],
            "equation_profile": ARM_DESCRIPTOR_VALUES[arm]["equation_profile"],
            "schema": ARM_DESCRIPTOR_SCHEMA,
            "source_policy": ARM_DESCRIPTOR_VALUES[arm]["source_policy"],
        }
        arm_descriptors.append({
            "arm": arm,
            "descriptor": descriptor,
            "descriptor_sha256": sha256_bytes(canonical_bytes(descriptor)),
        })
    # Rebind the cell-local receipt to the canonical top-level descriptor.
    descriptor_by_arm = {
        entry["arm"]: entry["descriptor_sha256"] for entry in arm_descriptors
    }
    for cell in cells:
        for oracle in cell["oracles"]:
            oracle["arm_descriptor_sha256"] = descriptor_by_arm[oracle["arm"]]
    return {
        "arm_descriptors": arm_descriptors,
        "campaign": CAMPAIGN,
        "cells": cells,
        "implementation_git_commit": implementation_commit,
        "internal_deadline_seconds": INTERNAL_DEADLINE_SECONDS,
        "invocation_budget": INVOCATION_BUDGET,
        "minimum_invocations": MINIMUM_INVOCATIONS,
        "panel_replicates": REPLICATES,
        "role": role,
        "schema": WORKER_CONFIG_SCHEMA,
        "supported_arms": list(ROLE_ARMS[role]),
        "target_cpu": cpu,
        "worker_git_commit": worker_commit,
    }


def synthetic_contracts(
    cpu: int, worker_commit: str, current_commit: str,
) -> Tuple[WorkerContract, WorkerContract]:
    current_config = synthetic_worker_config(
        "current", cpu, worker_commit, current_commit)
    parent_config = synthetic_worker_config(
        "parent", cpu, worker_commit, PARENT_COMMIT)
    current = validate_worker_config(
        current_config, "current", cpu, worker_commit, current_commit)
    parent = validate_worker_config(
        parent_config, "parent", cpu, worker_commit, PARENT_COMMIT)
    validate_cross_role_contracts(current, parent)
    return current, parent


def synthetic_ratio(comparison: str, scope: str) -> float:
    if comparison in AA_COMPARISONS or comparison in EQUIVALENCE_COMPARISONS:
        return 1.0
    if comparison == "C/E":
        return 0.95 if scope in SYSTEMATIC_SCOPES else 1.0
    if comparison == "C/L":
        return 0.94 if scope in SYSTEMATIC_SCOPES else 1.10
    raise AssertionError("unsealed comparison")


def synthetic_panel(
    replicate: int, k: int, block_bytes: int, scope: str,
    comparison: str, first_sequence: int,
) -> Dict[str, Any]:
    left_arm, right_arm = COMPARISON_ARMS[comparison]
    order = panel_order(replicate, k, block_bytes, scope, comparison)
    n = invocations_by_replicate(k)[replicate]
    ratio = synthetic_ratio(comparison, scope)
    base = n * 10_000_000
    left_elapsed = max(n, int(round(base * ratio)))
    right_elapsed = base
    slots: List[Dict[str, Any]] = []
    for slot_index, side in enumerate(sides_for(order)):
        arm = left_arm if side == "left" else right_arm
        elapsed = left_elapsed if side == "left" else right_elapsed
        sequence_first, sequence_last = expected_panel_slot_sequence(
            first_sequence, n, slot_index)
        slot = {
            "arm": arm,
            "borrowed_systematic_invocations": (
                n * k if arm == "C" and scope in SYSTEMATIC_SCOPES else 0),
            "constructor_ns": elapsed // 3 if scope == "fresh-systematic" else 0,
            "elapsed_ns": elapsed,
            "first_sequence": sequence_first,
            "init_first_repair_ns": elapsed // 2 if scope == "fresh-repair" else 0,
            "invocation_count": n,
            "last_sequence": sequence_last,
            "response_stream_sha256": _synthetic_hash(
                "responses-{}-{}-{}-{}-{}-{}".format(
                    replicate, k, block_bytes, scope, comparison, slot_index)),
            "side": side,
            "systematic_ns": elapsed // 2 if scope == "fresh-systematic" else 0,
        }
        slots.append(slot)
    return {
        "K": k,
        "block_bytes": block_bytes,
        "comparison": comparison,
        "first_sequence": first_sequence,
        "invocations_per_slot": n,
        "last_sequence": first_sequence + 2 + 8 * n - 1,
        "left_arm": left_arm,
        "order": order,
        "panel_key_sha256": panel_key_sha256(
            replicate, k, block_bytes, scope, comparison),
        "replicate": replicate,
        "right_arm": right_arm,
        "schema": PANEL_SCHEMA,
        "scope": scope,
        "slots": slots,
        "warmup_response_stream_sha256": _synthetic_hash(
            "warmup-{}-{}-{}-{}-{}".format(
                replicate, k, block_bytes, scope, comparison)),
        "warmup_sequences": [first_sequence, first_sequence + 1],
    }


def synthetic_transcript(
    cpu: int, worker_commit: str, current_commit: str,
) -> bytes:
    current, parent = synthetic_contracts(cpu, worker_commit, current_commit)
    config = make_screen_config(current, parent)
    records: List[Dict[str, Any]] = [config]
    next_sequence = 0
    role_counts = {"current": 0, "parent": 0}
    for replicate in range(REPLICATES):
        for k, block_bytes in CELLS:
            for scope in SCOPES:
                for comparison, left, right in COMPARISONS:
                    panel = synthetic_panel(
                        replicate, k, block_bytes, scope, comparison,
                        next_sequence)
                    records.append(panel)
                    n = invocations_by_replicate(k)[replicate]
                    role_counts[ARM_ROLE[left]] += 1 + 4 * n
                    role_counts[ARM_ROLE[right]] += 1 + 4 * n
                    next_sequence = panel["last_sequence"] + 1
    records.append({
        "current_invocations": role_counts["current"],
        "panel_count": EXPECTED_PANEL_COUNT,
        "parent_invocations": role_counts["parent"],
        "raw_config_sha256": sha256_bytes(canonical_bytes(config)),
        "schema": RAW_TERMINAL_SCHEMA,
        "status": "complete",
    })
    return b"".join(canonical_bytes(record) + b"\n" for record in records)


def synthetic_logs(default_ratio: float = 1.0) -> Dict[Tuple[str, int, int, str], List[float]]:
    return {
        (comparison, k, block_bytes, scope): [math.log(default_ratio)] * REPLICATES
        for comparison, _, _ in COMPARISONS
        for k, block_bytes in CELLS
        for scope in SCOPES
    }


def expect_invalid(callback: Any, label: str) -> None:
    try:
        callback()
    except ValidationError:
        return
    fail("selftest mutation was accepted: {}".format(label))


def selftest() -> int:
    cpu = 7
    worker_commit = "1" * 40
    current_commit = "2" * 40
    if EXPECTED_PANEL_COUNT != 4320:
        fail("selftest panel roster constant differs")
    if invocations_by_replicate(8) != [683] * 8 + [682] * 4:
        fail("selftest K=8 replicate allocation differs")
    if invocations_by_replicate(64000) != [2] * 12:
        fail("selftest K=64000 replicate allocation differs")
    left = 950000
    right = 1000000
    for order in ("ABBA", "BAAB"):
        values = [left if side == "left" else right for side in sides_for(order)]
        if abs(math.exp(lane_contrast(values, order)) - 0.95) > 1e-12:
            fail("selftest lane contrast differs for {}".format(order))

    raw = synthetic_transcript(cpu, worker_commit, current_commit)
    config, terminal, statistics, failures = parse_raw_transcript(
        raw, cpu, worker_commit, current_commit)
    if failures or terminal["panel_count"] != EXPECTED_PANEL_COUNT:
        fail("synthetic passing replay did not pass")
    replay = parse_raw_transcript(raw, cpu, worker_commit, current_commit)
    if canonical_bytes(replay[2]) != canonical_bytes(statistics):
        fail("synthetic replay is not deterministic")
    if config["current_worker"]["cells"][0]["oracles"][0][
        "borrowed_systematic_ids"] != 8:
        fail("synthetic borrowed-source receipt differs")

    expect_invalid(
        lambda: parse_canonical_line(b'{"a":1,"a":2}\n', "duplicate"),
        "duplicate JSON key")
    expect_invalid(
        lambda: parse_canonical_line(b'{"a": 1}\n', "spacing"),
        "noncanonical JSON spacing")
    current, parent = synthetic_contracts(cpu, worker_commit, current_commit)
    bad_parent = json.loads(canonical_bytes(parent.config).decode("ascii"))
    bad_parent["cells"][0]["oracles"][0]["repair_sha256"] = "f" * 64
    expect_invalid(
        lambda: validate_cross_role_contracts(
            current, validate_worker_config(
                bad_parent, "parent", cpu, worker_commit, PARENT_COMMIT)),
        "WH2 cross-role repair mismatch")

    # Exercise invalid-versus-reject adjudication without constructing another
    # 4,320-record transcript.  The preregistration makes A/A drift invalid
    # evidence, while an otherwise sound failed candidate gate is a rejection.
    logs = synthetic_logs()
    for key in list(logs.keys()):
        comparison, _, _, scope = key
        if comparison == "C/E":
            logs[key] = [math.log(0.95 if scope in SYSTEMATIC_SCOPES else 1.0)] * REPLICATES
        elif comparison == "C/L":
            logs[key] = [math.log(0.94 if scope in SYSTEMATIC_SCOPES else 1.50)] * REPLICATES
    pass_stats, pass_failures = adjudicate(logs)
    if pass_failures or not pass_stats["gate_policy"]["wh1_repair_is_descriptive_only"]:
        fail("selftest passing statistics differ")
    logs[("D/D", 8, 64, "prebuilt-systematic")] = [math.log(1.03)] * REPLICATES
    _, invalid_failures = adjudicate(logs)
    if (
        not any(item.startswith("aa_ci:D/D") for item in invalid_failures)
        or classify_outcome([], True, invalid_failures) != "invalid"
    ):
        fail("selftest invalid A/A evidence was not detected")
    logs = synthetic_logs()
    for key in list(logs.keys()):
        comparison, _, _, scope = key
        if comparison == "C/E":
            logs[key] = [math.log(
                1.03 if scope in SYSTEMATIC_SCOPES else 1.0
            )] * REPLICATES
        elif comparison == "C/L":
            logs[key] = [math.log(
                0.94 if scope in SYSTEMATIC_SCOPES else 1.50
            )] * REPLICATES
    _, reject_failures = adjudicate(logs)
    if (
        not reject_failures
        or any(item.startswith("aa_ci:") for item in reject_failures)
        or classify_outcome([], True, reject_failures) != "reject"
    ):
        fail("selftest valid scientific rejection was not detected")

    current, _parent = synthetic_contracts(cpu, worker_commit, current_commit)
    for seq, scope in enumerate(SCOPES):
        oracle = current.oracles[(8, 64, "C")]
        invocation = {
            "K": 8,
            "arm": "C",
            "block_bytes": 64,
            "borrowed_systematic_invocations": 8 if scope in SYSTEMATIC_SCOPES else 0,
            "constructor_ns": 30 if scope == "fresh-systematic" else None,
            "correctness_sha256": expected_correctness_hash(
                current, 8, 64, "C", scope),
            "descriptor_sha256": oracle["arm_descriptor_sha256"],
            "elapsed_ns": 100,
            "init_first_repair_ns": 50 if scope == "fresh-repair" else None,
            "result": 0,
            "role": "current",
            "schema": WORKER_INVOCATION_SCHEMA,
            "scope": scope,
            "seq": seq,
            "systematic_ns": 60 if scope == "fresh-systematic" else None,
            "target_cpu": cpu,
        }
        validate_invocation(invocation, current, seq, "C", scope, 8, 64)
    bad_invocation = dict(invocation)
    bad_invocation["borrowed_systematic_invocations"] = 1
    expect_invalid(
        lambda: validate_invocation(
            bad_invocation, current, 3, "C", "prebuilt-repair", 8, 64),
        "borrowed eligibility mutation")
    validate_worker_terminal({
        "invocation_count": 4,
        "role": "current",
        "schema": WORKER_TERMINAL_SCHEMA,
        "status": "complete",
        "target_cpu": cpu,
    }, current, 4)

    compile(GUARDIAN_CODE, "<guardian-selftest>", "exec")
    compile(WORKER_STUB_CODE, "<worker-stub-selftest>", "exec")
    selftest_build_authority()
    with tempfile.TemporaryDirectory(prefix="wh2-facade-screen-selftest-") as root:
        artifact_path = Path(root) / "python-image"
        artifact_alias = Path(root) / "python-image-alias"
        artifact_data = b"synthetic-python-image\n"
        artifact_path.write_bytes(artifact_data)
        artifact_path.chmod(0o755)
        os.link(str(artifact_path), str(artifact_alias))
        linked = open_artifact(
            "linked Python image selftest", str(artifact_path),
            sha256_bytes(artifact_data), 4096, executable=True,
            expected_uid=os.getuid(), expected_gid=os.getgid(),
            expected_mode=0o755, expected_nlink=2,
        )
        try:
            if linked.before.st_nlink != 2:
                fail("selftest linked Python image receipt differs")
        finally:
            linked.close()
        multiply_linked = open_artifact(
            "system dynamic link-count selftest", str(artifact_path),
            sha256_bytes(artifact_data), 4096, executable=True,
            expected_uid=os.getuid(), expected_gid=os.getgid(),
            expected_mode=0o755, expected_nlink=None,
        )
        try:
            if multiply_linked.before.st_nlink != 2:
                fail("selftest system dynamic link-count receipt differs")
        finally:
            multiply_linked.close()
        expect_invalid(
            lambda: open_artifact(
                "wrong-link Python image selftest", str(artifact_path),
                sha256_bytes(artifact_data), 4096, executable=True,
                expected_uid=os.getuid(), expected_gid=os.getgid(),
                expected_mode=0o755,
            ),
            "Python image link-count mismatch",
        )
        output_path = str(Path(root) / CAMPAIGN)
        expected_parent = output_parent_identity(
            os.stat(root, follow_symlinks=False))
        bundle = OutputBundle(output_path, expected_parent)
        try:
            bundle.append_raw({"schema": "synthetic-output-selftest"})
            bundle.write_once(
                "provenance.json", {"schema": "synthetic-provenance"}, 4096)
            bundle.write_once(
                "summary.json", {"schema": "synthetic-summary"}, 4096)
            complete = bundle.seal("pass")
            if (
                complete["outcome"] != "pass" or
                stat.S_IMODE(os.stat(output_path).st_mode) != 0o500 or
                sorted(os.listdir(output_path)) !=
                    sorted(EXPECTED_OUTPUT_NAMES + ("COMPLETE",))
            ):
                fail("selftest immutable output seal differs")
            expect_invalid(lambda: OutputBundle(output_path, expected_parent),
                           "exclusive output reservation")
        finally:
            bundle.close()
            os.chmod(output_path, 0o700)

    print("WH2 V2 facade timing screen selftest passed")
    return 0


ARTIFACT_KEYS = {"bytes", "path", "sha256"}
BUILD_RECEIPT_KEYS = {
    "archive_closure", "build_argv", "build_type", "compile_argv",
    "compiler_path", "compiler_sha256", "configure_argv",
    "dt_needed", "dynamic_dependencies", "generator", "git_path", "git_sha256",
    "harness_git_commit", "harness_source_root", "harness_tree_oid",
    "harness_tree_listing_sha256",
    "implementation_git_commit", "implementation_library_path",
    "implementation_library_sha256", "link_argv", "object_closure",
    "implementation_source_root", "implementation_tree_oid",
    "implementation_tree_listing_sha256",
    "receipt_preimage_sha256", "role", "role_library_path",
    "role_library_sha256", "runpath", "schema", "source_manifest", "worker_path",
    "worker_sha256",
}
SOURCE_ARTIFACT_KEYS = {
    "authority", "bytes", "git_blob_oid", "git_mode", "path",
    "repository_relative_path", "sha256",
}
LOWER40_OR_64 = re.compile(r"^(?:[0-9a-f]{40}|[0-9a-f]{64})$")
REQUIRED_HARNESS_SOURCES = frozenset((
    "CMakeLists.txt",
    "bench/Wh2FrozenTrace.cpp",
    "bench/Wh2FrozenTrace.h",
    "bench/Wh2V2FacadeTimingScreen.py",
    "bench/Wh2V2FacadeTimingWorker.cpp",
    "bench/public_facade_timing/CMakeLists.txt",
))
REQUIRED_IMPLEMENTATION_SOURCES = frozenset((
    "CMakeLists.txt",
    "include/wirehair/wirehair.h",
))


def exact_absolute_path(value: Any, where: str) -> str:
    if (
        type(value) is not str or not value.startswith("/") or
        len(value) > 4096 or os.path.normpath(value) != value or "\x00" in value
    ):
        fail("{} is not a canonical absolute path".format(where))
    return value


def exact_argv(value: Any, where: str) -> List[str]:
    if (
        type(value) is not list or not value or len(value) > 4096 or
        any(type(item) is not str or not item or len(item) > 4096
            for item in value)
    ):
        fail("{} is not a bounded nonempty argv".format(where))
    return value


def validate_artifact_roster(
    value: Any, where: str, *, nonempty: bool = True,
) -> List[Dict[str, Any]]:
    if (
        type(value) is not list or len(value) > 65536 or
        (nonempty and not value)
    ):
        fail("{} artifact roster shape differs".format(where))
    paths: List[str] = []
    result: List[Dict[str, Any]] = []
    for index, artifact in enumerate(value):
        item_where = "{}[{}]".format(where, index)
        if type(artifact) is not dict:
            fail("{} is not an object".format(item_where))
        exact_keys(artifact, ARTIFACT_KEYS, item_where)
        path = exact_absolute_path(artifact["path"], item_where + ".path")
        exact_int(artifact["bytes"], 1, MAX_INT63, item_where + ".bytes")
        lower_hash(artifact["sha256"], item_where + ".sha256")
        paths.append(path)
        result.append(artifact)
    if paths != sorted(set(paths)):
        fail("{} paths are not sorted and unique".format(where))
    return result


def validate_source_artifact_roster(value: Any, where: str) -> List[Dict[str, Any]]:
    if type(value) is not list or not value or len(value) > 65536:
        fail("{} source roster shape differs".format(where))
    paths: List[str] = []
    authorities: Dict[str, set] = {"harness": set(), "implementation": set()}
    result: List[Dict[str, Any]] = []
    for index, entry in enumerate(value):
        item_where = "{}[{}]".format(where, index)
        if type(entry) is not dict:
            fail("{} is not an object".format(item_where))
        exact_keys(entry, SOURCE_ARTIFACT_KEYS, item_where)
        if entry["authority"] not in authorities:
            fail("{} authority differs".format(item_where))
        path = exact_absolute_path(entry["path"], item_where + ".path")
        relative = entry["repository_relative_path"]
        if (
            type(relative) is not str or not relative or relative.startswith("/") or
            os.path.normpath(relative) != relative or relative.startswith("../") or
            len(relative) > 4096
        ):
            fail("{} repository-relative path differs".format(item_where))
        oid = entry["git_blob_oid"]
        if type(oid) is not str or LOWER40_OR_64.fullmatch(oid) is None:
            fail("{} Git blob OID differs".format(item_where))
        if entry["git_mode"] not in ("100644", "100755"):
            fail("{} Git mode differs".format(item_where))
        exact_int(entry["bytes"], 1, MAX_INT63, item_where + ".bytes")
        lower_hash(entry["sha256"], item_where + ".sha256")
        paths.append(path)
        authorities[entry["authority"]].add(relative)
        result.append(entry)
    if paths != sorted(set(paths)):
        fail("{} paths are not sorted and unique".format(where))
    if not REQUIRED_HARNESS_SOURCES.issubset(authorities["harness"]):
        fail("{} lacks required controller/worker/trace/CMake sources".format(where))
    if not REQUIRED_IMPLEMENTATION_SOURCES.issubset(
        authorities["implementation"]
    ):
        fail("{} lacks required public header/build sources".format(where))
    # A role library built from only the two minimum entries is not a transitive
    # implementation closure.  Require a substantive codec/core source roster.
    implementation_code = [
        relative for relative in authorities["implementation"]
        if relative.startswith(("codec/", "src/")) and
        relative.endswith((".c", ".cc", ".cpp", ".h", ".hpp"))
    ]
    if len(implementation_code) < 8:
        fail("{} lacks the facade/core transitive source closure".format(where))
    return result


def parse_canonical_document(data: bytes, where: str, maximum: int) -> Dict[str, Any]:
    if not data or len(data) > maximum:
        fail("{} exceeds its bounded document cap".format(where))
    return parse_canonical_line(data, where, maximum)


def validate_build_receipt(
    receipt: Dict[str, Any], role: str, harness_commit: str,
    implementation_commit: str, worker_path: str, worker_sha256: str,
    library_path: str, library_sha256: str, implementation_library_path: str,
    implementation_library_sha256: str, harness_source_root: str,
    controller_path: str, controller_sha256: str,
) -> None:
    where = "{} build receipt".format(role)
    exact_keys(receipt, BUILD_RECEIPT_KEYS, where)
    exact_string(receipt["schema"], BUILD_RECEIPT_SCHEMA, where + ".schema")
    exact_string(receipt["role"], role, where + ".role")
    exact_string(receipt["harness_git_commit"], harness_commit,
                 where + ".harness_git_commit")
    exact_string(receipt["implementation_git_commit"], implementation_commit,
                 where + ".implementation_git_commit")
    exact_string(receipt["worker_path"], worker_path, where + ".worker_path")
    exact_string(receipt["worker_sha256"], worker_sha256,
                 where + ".worker_sha256")
    exact_string(receipt["role_library_path"], library_path,
                 where + ".role_library_path")
    exact_string(receipt["role_library_sha256"], library_sha256,
                 where + ".role_library_sha256")
    exact_string(receipt["implementation_library_path"],
                 implementation_library_path,
                 where + ".implementation_library_path")
    exact_string(receipt["implementation_library_sha256"],
                 implementation_library_sha256,
                 where + ".implementation_library_sha256")
    if library_sha256 != implementation_library_sha256:
        fail("{} runtime and implementation libraries differ in bytes".format(where))
    exact_string(receipt["generator"], "Ninja", where + ".generator")
    exact_string(receipt["build_type"], "Release", where + ".build_type")
    compiler = exact_absolute_path(receipt["compiler_path"],
                                   where + ".compiler_path")
    if compiler != "/usr/bin/x86_64-linux-gnu-g++-13":
        fail("{} compiler path differs".format(where))
    lower_hash(receipt["compiler_sha256"], where + ".compiler_sha256")
    git_path = exact_absolute_path(receipt["git_path"], where + ".git_path")
    if git_path != "/usr/bin/git":
        fail("{} Git executable path differs".format(where))
    lower_hash(receipt["git_sha256"], where + ".git_sha256")
    for authority in ("harness", "implementation"):
        tree_oid = receipt[authority + "_tree_oid"]
        if type(tree_oid) is not str or LOWER40_OR_64.fullmatch(tree_oid) is None:
            fail("{} {} tree OID differs".format(where, authority))
        lower_hash(receipt[authority + "_tree_listing_sha256"],
                   "{} {} tree listing SHA-256".format(where, authority))
    actual_harness_root = exact_absolute_path(
        receipt["harness_source_root"], where + ".harness_source_root")
    actual_implementation_root = exact_absolute_path(
        receipt["implementation_source_root"],
        where + ".implementation_source_root")
    exact_string(actual_harness_root, harness_source_root,
                 where + ".harness_source_root")
    for root_name, root_path in (
        ("harness", actual_harness_root),
        ("implementation", actual_implementation_root),
    ):
        try:
            root_info = os.stat(root_path, follow_symlinks=False)
        except OSError as exc:
            fail("{} {} source root cannot be inspected: {}".format(
                where, root_name, exc))
        if (
            not stat.S_ISDIR(root_info.st_mode) or
            os.path.realpath(root_path) != root_path
        ):
            fail("{} {} source root is not a canonical real directory".format(
                where, root_name))
    configure_argv = exact_argv(receipt["configure_argv"],
                                where + ".configure_argv")
    build_argv = exact_argv(receipt["build_argv"], where + ".build_argv")
    compile_argv = exact_argv(receipt["compile_argv"], where + ".compile_argv")
    link_argv = exact_argv(receipt["link_argv"], where + ".link_argv")
    if not any(item.endswith("cmake") for item in configure_argv[:1]):
        fail("{} configure command does not use CMake".format(where))
    if not any(item.endswith("cmake") for item in build_argv[:1]):
        fail("{} build command does not use CMake".format(where))
    for required in ("-O3", "-DNDEBUG"):
        if required not in compile_argv:
            fail("{} compile command lacks {}".format(where, required))
    if library_path not in link_argv:
        fail("{} link command lacks the exact role library".format(where))
    needed = receipt["dt_needed"]
    if (
        type(needed) is not list or not needed or len(needed) > 64 or
        any(type(item) is not str or not item or len(item) > 256 or "/" in item
            for item in needed) or
        len(set(needed)) != len(needed) or "libwirehair.so.2" not in needed
    ):
        fail("{} DT_NEEDED closure differs".format(where))
    exact_string(receipt["runpath"], "$ORIGIN", where + ".runpath")
    source_manifest = validate_source_artifact_roster(
        receipt["source_manifest"], where + ".source_manifest")
    validate_artifact_roster(receipt["object_closure"], where + ".object_closure")
    validate_artifact_roster(
        receipt["archive_closure"], where + ".archive_closure", nonempty=False)
    dynamic = validate_artifact_roster(
        receipt["dynamic_dependencies"], where + ".dynamic_dependencies")
    dynamic_basenames = [Path(item["path"]).name for item in dynamic]
    matched_needed = []
    for basename in dynamic_basenames:
        matches = [
            soname for soname in needed
            if basename == soname or basename.startswith(soname + ".")
        ]
        if len(matches) != 1:
            fail("{} resolved dynamic dependency name is ambiguous".format(where))
        matched_needed.append(matches[0])
    if (
        len(set(dynamic_basenames)) != len(dynamic_basenames) or
        sorted(matched_needed) != sorted(needed)
    ):
        fail("{} resolved dynamic closure differs from DT_NEEDED".format(where))
    source_by_authority_relative: Dict[Tuple[str, str], Dict[str, Any]] = {}
    roots = {
        "harness": actual_harness_root,
        "implementation": actual_implementation_root,
    }
    for entry in source_manifest:
        key = (entry["authority"], entry["repository_relative_path"])
        if key in source_by_authority_relative:
            fail("{} source manifest repeats {}:{}".format(where, *key))
        source_by_authority_relative[key] = entry
        expected_path = os.path.join(roots[entry["authority"]], key[1])
        if entry["path"] != expected_path or os.path.realpath(expected_path) != expected_path:
            fail("{} source {}:{} is outside its exact root".format(
                where, *key))
    controller_entry = source_by_authority_relative.get(
        ("harness", "bench/Wh2V2FacadeTimingScreen.py"))
    if (
        controller_entry is None or controller_entry["path"] != controller_path or
        controller_entry["sha256"] != controller_sha256
    ):
        fail("{} source manifest does not crossbind the controller".format(where))
    worker_source = source_by_authority_relative.get(
        ("harness", "bench/Wh2V2FacadeTimingWorker.cpp"))
    if (
        worker_source is None or
        worker_source["sha256"] != FROZEN_WORKER_SOURCE_SHA256
    ):
        fail("{} source manifest does not crossbind the frozen worker".format(where))
    if library_path not in {item["path"] for item in dynamic}:
        fail("{} dynamic closure lacks the exact role library".format(where))
    preimage_hash = lower_hash(
        receipt["receipt_preimage_sha256"], where + ".receipt_preimage_sha256")
    preimage = dict(receipt)
    preimage["receipt_preimage_sha256"] = None
    if sha256_bytes(canonical_bytes(preimage)) != preimage_hash:
        fail("{} self-preimage differs".format(where))


def file_sha256_fd(fd: int, size: int) -> str:
    digest = hashlib.sha256()
    offset = 0
    while offset < size:
        block = os.pread(fd, min(1024 * 1024, size - offset), offset)
        if not block:
            fail("short read while hashing an artifact")
        digest.update(block)
        offset += len(block)
    return digest.hexdigest()


def stat_receipt(value: os.stat_result) -> Dict[str, Any]:
    return {
        "ctime_ns": value.st_ctime_ns,
        "device": value.st_dev,
        "gid": value.st_gid,
        "inode": value.st_ino,
        "mode": stat.S_IMODE(value.st_mode),
        "mtime_ns": value.st_mtime_ns,
        "nlink": value.st_nlink,
        "size": value.st_size,
        "uid": value.st_uid,
    }


@dataclass
class OpenArtifact:
    name: str
    path: str
    expected_sha256: str
    fd: int
    before: os.stat_result
    sha256: str

    def receipt(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "path": self.path,
            "sha256": self.sha256,
            "stat": stat_receipt(self.before),
        }

    def verify_unchanged(self) -> None:
        after = os.fstat(self.fd)
        try:
            path_after = os.stat(self.path, follow_symlinks=False)
        except OSError as exc:
            fail("{} path recheck failed: {}".format(self.name, exc))
        if (
            stat_receipt(after) != stat_receipt(self.before) or
            stat_receipt(path_after) != stat_receipt(self.before) or
            file_sha256_fd(self.fd, self.before.st_size) != self.sha256
        ):
            fail("{} changed after sealing".format(self.name))

    def close(self) -> None:
        if self.fd >= 0:
            os.close(self.fd)
            self.fd = -1


def open_artifact(
    name: str, path: str, expected_sha256: str, maximum_bytes: int,
    *, executable: bool = False, expected_uid: Optional[int] = None,
    expected_gid: Optional[int] = None, expected_mode: Optional[int] = None,
    expected_nlink: Optional[int] = 1,
) -> OpenArtifact:
    exact_absolute_path(path, name + " path")
    lower_hash(expected_sha256, name + " expected SHA-256")
    if expected_nlink is not None:
        exact_int(expected_nlink, 1, MAX_INT63, name + " expected link count")
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        fd = os.open(path, flags)
    except OSError as exc:
        fail("could not open {}: {}".format(name, exc))
    try:
        info = os.fstat(fd)
        mode = stat.S_IMODE(info.st_mode)
        if (
            not stat.S_ISREG(info.st_mode) or
            (expected_nlink is not None and info.st_nlink != expected_nlink) or
            info.st_size < 1 or info.st_size > maximum_bytes or mode & 0o022
        ):
            fail("{} file policy differs".format(name))
        if expected_uid is not None and info.st_uid != expected_uid:
            fail("{} owner UID differs".format(name))
        if expected_gid is not None and info.st_gid != expected_gid:
            fail("{} owner GID differs".format(name))
        if expected_mode is not None and mode != expected_mode:
            fail("{} mode differs".format(name))
        if executable and mode & 0o111 == 0:
            fail("{} is not executable".format(name))
        digest = file_sha256_fd(fd, info.st_size)
        if digest != expected_sha256:
            fail("{} SHA-256 differs".format(name))
        path_info = os.stat(path, follow_symlinks=False)
        if stat_receipt(path_info) != stat_receipt(info):
            fail("{} path and descriptor differ".format(name))
        return OpenArtifact(name, path, expected_sha256, fd, info, digest)
    except BaseException:
        os.close(fd)
        raise


def read_artifact(artifact: OpenArtifact, maximum: int) -> bytes:
    if artifact.before.st_size > maximum:
        fail("{} exceeds its read cap".format(artifact.name))
    data = bytearray()
    offset = 0
    while offset < artifact.before.st_size:
        block = os.pread(
            artifact.fd, min(1024 * 1024, artifact.before.st_size - offset), offset)
        if not block:
            fail("short read from {}".format(artifact.name))
        data.extend(block)
        offset += len(block)
    return bytes(data)


def git_blob_oid_fd(fd: int, size: int, hexadecimal_length: int) -> str:
    if hexadecimal_length == 40:
        digest = hashlib.sha1()
    elif hexadecimal_length == 64:
        digest = hashlib.sha256()
    else:
        fail("unsupported Git object identifier width")
    digest.update("blob {}\0".format(size).encode("ascii"))
    offset = 0
    while offset < size:
        block = os.pread(fd, min(1024 * 1024, size - offset), offset)
        if not block:
            fail("short read while computing Git blob identity")
        digest.update(block)
        offset += len(block)
    return digest.hexdigest()


def run_held_git(
    git: OpenArtifact, root: str, arguments: Sequence[str], maximum: int,
) -> bytes:
    if maximum < 1 or maximum > 64 * 1024 * 1024:
        fail("sealed Git output cap differs")
    git.verify_unchanged()
    argv = [
        git.path,
        "-c", "core.fsmonitor=false",
        "-c", "core.filemode=false",
        "-c", "core.checkStat=default",
        "-c", "core.trustctime=true",
        "-c", "safe.directory=" + root,
        "-C", root,
    ] + list(arguments)
    if any(type(item) is not str or not item or len(item) > 4096 for item in argv):
        fail("sealed Git argv differs")
    executable = "/proc/self/fd/{}".format(git.fd)
    with tempfile.TemporaryFile() as stdout_file, tempfile.TemporaryFile() as stderr_file:
        process = subprocess.Popen(
            argv, executable=executable, stdin=subprocess.DEVNULL,
            stdout=stdout_file, stderr=stderr_file, close_fds=True,
            pass_fds=(git.fd,), start_new_session=True,
            env={
                "LANG": "C.UTF-8", "LC_ALL": "C.UTF-8",
                "PATH": "/usr/bin:/bin", "TZ": "UTC",
                "GIT_CONFIG_GLOBAL": "/dev/null",
                "GIT_CONFIG_NOSYSTEM": "1",
                "GIT_NO_LAZY_FETCH": "1",
                "GIT_NO_REPLACE_OBJECTS": "1",
                "GIT_OPTIONAL_LOCKS": "0",
                "GIT_TERMINAL_PROMPT": "0",
            },
        )
        try:
            code = process.wait(timeout=15.0)
        except subprocess.TimeoutExpired:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except OSError:
                process.kill()
            process.wait(timeout=5.0)
            fail("sealed Git command exceeded its deadline")
        stdout_size = os.fstat(stdout_file.fileno()).st_size
        stderr_size = os.fstat(stderr_file.fileno()).st_size
        if stdout_size > maximum or stderr_size > 64 * 1024:
            fail("sealed Git command exceeded its output cap")
        stdout_file.seek(0)
        stderr_file.seek(0)
        output = stdout_file.read()
        error = stderr_file.read()
    git.verify_unchanged()
    if code != 0 or error:
        fail("sealed Git command failed rc={} stderr_sha256={}".format(
            code, sha256_bytes(error)))
    return output


def parse_git_tree_listing(data: bytes, where: str) -> List[Dict[str, str]]:
    if not data or not data.endswith(b"\0"):
        fail("{} is not a nonempty NUL-terminated Git tree".format(where))
    result: List[Dict[str, str]] = []
    paths: set = set()
    for record in data[:-1].split(b"\0"):
        try:
            header, path_bytes = record.split(b"\t", 1)
            mode, kind, oid = header.split(b" ")
            path = path_bytes.decode("utf-8")
            mode_text = mode.decode("ascii")
            oid_text = oid.decode("ascii")
        except (ValueError, UnicodeDecodeError):
            fail("{} record is malformed".format(where))
        if (
            kind != b"blob" or mode_text not in ("100644", "100755") or
            LOWER40_OR_64.fullmatch(oid_text) is None or not path or
            path.startswith("/") or path.startswith("../") or
            os.path.normpath(path) != path or "\x00" in path or path in paths
        ):
            fail("{} contains a nonregular or noncanonical entry".format(where))
        paths.add(path)
        result.append({"git_blob_oid": oid_text, "git_mode": mode_text,
                       "repository_relative_path": path})
    if not result or len(result) > 65536:
        fail("{} entry count differs".format(where))
    return result


def parse_git_stage_listing(data: bytes, where: str) -> List[Dict[str, str]]:
    if not data or not data.endswith(b"\0"):
        fail("{} is not a nonempty NUL-terminated Git index".format(where))
    result: List[Dict[str, str]] = []
    paths: set = set()
    for record in data[:-1].split(b"\0"):
        try:
            header, path_bytes = record.split(b"\t", 1)
            mode, oid, stage = header.split(b" ")
            path = path_bytes.decode("utf-8")
            mode_text = mode.decode("ascii")
            oid_text = oid.decode("ascii")
        except (ValueError, UnicodeDecodeError):
            fail("{} record is malformed".format(where))
        if (
            stage != b"0" or mode_text not in ("100644", "100755") or
            LOWER40_OR_64.fullmatch(oid_text) is None or not path or
            path.startswith("/") or path.startswith("../") or
            os.path.normpath(path) != path or path in paths
        ):
            fail("{} contains a non-stage-0 or noncanonical entry".format(where))
        paths.add(path)
        result.append({"git_blob_oid": oid_text, "git_mode": mode_text,
                       "repository_relative_path": path})
    if not result or len(result) > 65536:
        fail("{} entry count differs".format(where))
    return result


def verify_tracked_worktree(
    root: str, entries: Sequence[Mapping[str, str]], authority: str,
    expected_uid: int = 1000, expected_gid: int = 1000,
    expected_directory_uid: int = 0, expected_directory_gid: int = 0,
    expected_directory_mode: int = 0o555,
) -> Dict[str, Any]:
    total_bytes = 0
    manifest_digest = hashlib.sha256()
    directories = {root}
    for entry in entries:
        parent = os.path.dirname(os.path.join(
            root, entry["repository_relative_path"]))
        while parent != root:
            if not parent.startswith(root + os.sep):
                fail("{} tracked directory escapes its root".format(authority))
            directories.add(parent)
            parent = os.path.dirname(parent)
    directory_digest = hashlib.sha256()
    for directory in sorted(directories):
        info = os.stat(directory, follow_symlinks=False)
        if (
            not stat.S_ISDIR(info.st_mode) or
            info.st_uid != expected_directory_uid or
            info.st_gid != expected_directory_gid or
            stat.S_IMODE(info.st_mode) != expected_directory_mode or
            os.path.realpath(directory) != directory
        ):
            fail("{} tracked directory {} seal differs".format(
                authority, directory))
        directory_digest.update(canonical_bytes({
            "path": os.path.relpath(directory, root),
            "stat": stat_receipt(info),
        }) + b"\n")
    for entry in entries:
        relative = entry["repository_relative_path"]
        path = os.path.join(root, relative)
        if os.path.realpath(path) != path:
            fail("{} tracked path {} traverses a link".format(authority, relative))
        flags = (os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
                 getattr(os, "O_NOFOLLOW", 0))
        try:
            fd = os.open(path, flags)
        except OSError as exc:
            fail("{} tracked path {} cannot be opened: {}".format(
                authority, relative, exc))
        try:
            before = os.fstat(fd)
            path_before = os.stat(path, follow_symlinks=False)
            if (
                not stat.S_ISREG(before.st_mode) or before.st_nlink != 1 or
                before.st_uid != expected_uid or before.st_gid != expected_gid or
                stat.S_IMODE(before.st_mode) != 0o444 or before.st_size < 1 or
                before.st_size > 1024 * 1024 * 1024 or
                stat_receipt(before) != stat_receipt(path_before)
            ):
                fail("{} tracked path {} owner/type/mode seal differs".format(
                    authority, relative))
            total_bytes += before.st_size
            if total_bytes > 4 * 1024 * 1024 * 1024:
                fail("{} tracked worktree exceeds its total byte cap".format(
                    authority))
            blob_oid = git_blob_oid_fd(
                fd, before.st_size, len(entry["git_blob_oid"]))
            if blob_oid != entry["git_blob_oid"]:
                fail("{} tracked path {} content differs from HEAD".format(
                    authority, relative))
            content_sha256 = file_sha256_fd(fd, before.st_size)
            after = os.fstat(fd)
            path_after = os.stat(path, follow_symlinks=False)
            if (
                stat_receipt(after) != stat_receipt(before) or
                stat_receipt(path_after) != stat_receipt(before)
            ):
                fail("{} tracked path {} changed while hashing".format(
                    authority, relative))
            manifest_digest.update(canonical_bytes({
                "bytes": before.st_size,
                "git_blob_oid": blob_oid,
                "git_mode": entry["git_mode"],
                "repository_relative_path": relative,
                "sha256": content_sha256,
            }) + b"\n")
        finally:
            os.close(fd)
    return {
        "bytes": total_bytes,
        "entries": len(entries),
        "directory_manifest_sha256": directory_digest.hexdigest(),
        "worktree_manifest_sha256": manifest_digest.hexdigest(),
    }


def verify_git_repository(
    git: OpenArtifact, root: str, commit: str, authority: str,
    expected_tree_oid: str, expected_tree_listing_sha256: str,
    source_entries: Sequence[Mapping[str, Any]],
) -> Dict[str, Any]:
    top = run_held_git(git, root, ("rev-parse", "--show-toplevel"), 8192)
    if top != root.encode("utf-8") + b"\n":
        fail("{} Git top-level differs from its canonical root".format(authority))
    head = run_held_git(git, root, ("rev-parse", "--verify", "HEAD"), 8192)
    if head != commit.encode("ascii") + b"\n":
        fail("{} source root is not at its exact commit".format(authority))
    branch = run_held_git(
        git, root, ("rev-parse", "--abbrev-ref", "HEAD"), 8192)
    if branch != b"HEAD\n":
        fail("{} source root is not detached".format(authority))
    status = run_held_git(
        git, root,
        ("status", "--porcelain=v1", "-z", "--untracked-files=all"),
        32 * 1024 * 1024)
    if status:
        fail("{} source root is not full-tree clean (status sha256 {})".format(
            authority, sha256_bytes(status)))
    untracked = run_held_git(
        git, root, ("ls-files", "--others", "--exclude-standard", "-z"),
        32 * 1024 * 1024)
    ignored = run_held_git(
        git, root,
        ("ls-files", "--others", "--ignored", "--exclude-standard", "-z"),
        32 * 1024 * 1024)
    if untracked or ignored:
        fail("{} source root contains untracked/ignored files ({}/{})".format(
            authority, sha256_bytes(untracked), sha256_bytes(ignored)))
    tree_oid = run_held_git(
        git, root, ("rev-parse", "--verify", commit + "^{tree}"), 8192)
    if tree_oid != expected_tree_oid.encode("ascii") + b"\n":
        fail("{} root tree OID differs".format(authority))
    tree = run_held_git(
        git, root,
        ("ls-tree", "-r", "-z", "--full-tree", commit),
        32 * 1024 * 1024)
    if not tree or not tree.endswith(b"\0"):
        fail("{} committed tree listing differs".format(authority))
    tree_entries = parse_git_tree_listing(tree, authority + " HEAD tree")
    tree_by_path = {
        entry["repository_relative_path"]: entry for entry in tree_entries
    }
    if len(source_entries) != len(tree_entries):
        fail("{} source manifest is not the complete tracked tree".format(authority))
    for source_entry in source_entries:
        relative = source_entry["repository_relative_path"]
        tree_entry = tree_by_path.get(relative)
        if tree_entry is None or tree_entry != {
            "git_blob_oid": source_entry["git_blob_oid"],
            "git_mode": source_entry["git_mode"],
            "repository_relative_path": relative,
        }:
            fail("{} manifest source {} lacks exact commit membership".format(
                authority, relative))
    if {
        entry["repository_relative_path"] for entry in source_entries
    } != set(tree_by_path):
        fail("{} source manifest roster differs from the tracked tree".format(
            authority))
    stage = run_held_git(
        git, root, ("ls-files", "--stage", "-z"), 32 * 1024 * 1024)
    stage_entries = parse_git_stage_listing(stage, authority + " stage-0 index")
    if stage_entries != tree_entries:
        fail("{} index does not exactly match the HEAD tree".format(authority))
    worktree = verify_tracked_worktree(root, tree_entries, authority)
    tree_listing_sha256 = sha256_bytes(tree)
    if tree_listing_sha256 != expected_tree_listing_sha256:
        fail("{} full-tree listing SHA-256 differs".format(authority))
    return {
        "authority": authority,
        "commit": commit,
        "root": root,
        "status_sha256": sha256_bytes(status),
        "tree_bytes": len(tree),
        "tree_oid": expected_tree_oid,
        "tree_sha256": tree_listing_sha256,
        "worktree": worktree,
    }


def verify_build_receipt_closure(
    receipt: Mapping[str, Any], role: str, git: OpenArtifact,
) -> str:
    """Rehash every receipt roster entry from a no-follow descriptor."""
    entries: Dict[str, Mapping[str, Any]] = {}
    for roster_name in (
        "source_manifest", "object_closure", "archive_closure",
        "dynamic_dependencies",
    ):
        for entry in receipt[roster_name]:
            path = entry["path"]
            previous = entries.get(path)
            if previous is not None and previous != entry:
                fail("{} build closure gives {} conflicting receipts".format(
                    role, path))
            entries[path] = entry
    compiler_entry = {
        "bytes": os.stat(receipt["compiler_path"], follow_symlinks=False).st_size,
        "path": receipt["compiler_path"],
        "sha256": receipt["compiler_sha256"],
    }
    entries[compiler_entry["path"]] = compiler_entry
    verified: List[Dict[str, Any]] = []
    source_paths = {entry["path"]: entry for entry in receipt["source_manifest"]}
    object_paths = {entry["path"] for entry in receipt["object_closure"]}
    archive_paths = {entry["path"] for entry in receipt["archive_closure"]}
    dynamic_paths = {entry["path"] for entry in receipt["dynamic_dependencies"]}
    for path in sorted(entries):
        entry = entries[path]
        artifact = open_artifact(
            "{} closure {}".format(role, path), path, entry["sha256"],
            1024 * 1024 * 1024,
            executable=(path == receipt["compiler_path"]),
            expected_nlink=(
                None
                if path in dynamic_paths
                and path != receipt["role_library_path"]
                else 1
            ),
        )
        try:
            exact_int(artifact.before.st_size, entry["bytes"], entry["bytes"],
                      "{} closure {} bytes".format(role, path))
            source_entry = source_paths.get(path)
            if source_entry is not None:
                if (
                    artifact.before.st_uid != 1000 or
                    artifact.before.st_gid != 1000 or
                    stat.S_IMODE(artifact.before.st_mode) != 0o444
                ):
                    fail("{} source {} owner/mode seal differs".format(role, path))
                oid = git_blob_oid_fd(
                    artifact.fd, artifact.before.st_size,
                    len(source_entry["git_blob_oid"]))
                if oid != source_entry["git_blob_oid"]:
                    fail("{} source {} live Git blob differs".format(role, path))
            elif path == receipt["compiler_path"]:
                if (
                    artifact.before.st_uid != 0 or artifact.before.st_gid != 0 or
                    stat.S_IMODE(artifact.before.st_mode) != 0o755
                ):
                    fail("{} compiler owner/mode seal differs".format(role))
            elif path in object_paths or path in archive_paths:
                if (
                    artifact.before.st_uid != 0 or artifact.before.st_gid != 0 or
                    stat.S_IMODE(artifact.before.st_mode) != 0o444
                ):
                    fail("{} build input {} owner/mode seal differs".format(
                        role, path))
            elif path in dynamic_paths:
                mode = stat.S_IMODE(artifact.before.st_mode)
                if (
                    artifact.before.st_uid != 0 or artifact.before.st_gid != 0 or
                    mode not in (0o444, 0o555, 0o644, 0o755) or
                    (path == receipt["role_library_path"] and mode != 0o444)
                ):
                    fail("{} dynamic input {} owner/mode seal differs".format(
                        role, path))
            verified.append({
                "bytes": artifact.before.st_size,
                "path": path,
                "sha256": artifact.sha256,
                "stat": stat_receipt(artifact.before),
            })
        finally:
            artifact.close()
    repository_receipts: List[Dict[str, Any]] = []
    authority_contracts = (
        ("harness", receipt["harness_source_root"],
         receipt["harness_git_commit"], receipt["harness_tree_oid"],
         receipt["harness_tree_listing_sha256"]),
        ("implementation", receipt["implementation_source_root"],
         receipt["implementation_git_commit"],
         receipt["implementation_tree_oid"],
         receipt["implementation_tree_listing_sha256"]),
    )
    for authority, root, commit, tree_oid, tree_listing_sha256 in authority_contracts:
        authority_sources = [
            entry for entry in receipt["source_manifest"]
            if entry["authority"] == authority
        ]
        repository_receipts.append(
            verify_git_repository(git, root, commit,
                                  "{} {}".format(role, authority),
                                  tree_oid, tree_listing_sha256,
                                  authority_sources))
    git.verify_unchanged()
    closure = {
        "artifacts": verified,
        "git": git.receipt(),
        "repositories": repository_receipts,
    }
    return sha256_bytes(canonical_bytes(closure))


def selftest_build_authority() -> None:
    def write_file(path: Path, data: bytes) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        fd = os.open(str(path), os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        try:
            if os.write(fd, data) != len(data):
                fail("synthetic authority fixture write was short")
        finally:
            os.close(fd)
        os.chmod(str(path), 0o444)

    with tempfile.TemporaryDirectory(
        prefix="wh2-facade-build-authority-selftest-"
    ) as temporary:
        base = Path(temporary)
        harness = base / "harness"
        implementation = base / "implementation"
        harness_relatives = sorted(REQUIRED_HARNESS_SOURCES)
        implementation_relatives = sorted(
            set(REQUIRED_IMPLEMENTATION_SOURCES) |
            {"codec/Synthetic{}.cpp".format(index) for index in range(8)})
        controller_relative = "bench/Wh2V2FacadeTimingScreen.py"
        controller_data = b"synthetic-controller\n"
        for relative in harness_relatives:
            write_file(
                harness / relative,
                controller_data if relative == controller_relative else b"harness\n")
        for relative in implementation_relatives:
            write_file(implementation / relative, b"implementation\n")
        worker_path = str(base / "worker")
        role_library = str(base / "libwirehair.so.2")
        implementation_library = str(base / "implementation-libwirehair.so.2")
        object_path = str(base / "worker.o")
        archive_path = str(base / "libfixture.a")
        for path in (worker_path, role_library, implementation_library,
                     object_path, archive_path):
            write_file(Path(path), b"fixture\n")

        source_manifest: List[Dict[str, Any]] = []
        for authority, root, relatives in (
            ("harness", harness, harness_relatives),
            ("implementation", implementation, implementation_relatives),
        ):
            for relative in relatives:
                data = (root / relative).read_bytes()
                digest = sha256_bytes(data)
                if authority == "harness" and relative == (
                    "bench/Wh2V2FacadeTimingWorker.cpp"
                ):
                    digest = FROZEN_WORKER_SOURCE_SHA256
                source_manifest.append({
                    "authority": authority,
                    "bytes": len(data),
                    "git_blob_oid": "1" * 40,
                    "git_mode": "100644",
                    "path": str(root / relative),
                    "repository_relative_path": relative,
                    "sha256": digest,
                })
        source_manifest.sort(key=lambda item: item["path"])
        artifact = lambda path: {
            "bytes": os.stat(path).st_size,
            "path": path,
            "sha256": sha256_bytes(Path(path).read_bytes()),
        }
        receipt: Dict[str, Any] = {
            "archive_closure": [artifact(archive_path)],
            "build_argv": ["/usr/bin/cmake", "--build", str(base)],
            "build_type": "Release",
            "compile_argv": ["/usr/bin/x86_64-linux-gnu-g++-13", "-O3",
                             "-DNDEBUG", "-c", object_path],
            "compiler_path": "/usr/bin/x86_64-linux-gnu-g++-13",
            "compiler_sha256": "2" * 64,
            "configure_argv": ["/usr/bin/cmake", "-G", "Ninja"],
            "dt_needed": ["libwirehair.so.2"],
            "dynamic_dependencies": [artifact(role_library)],
            "generator": "Ninja",
            "git_path": "/usr/bin/git",
            "git_sha256": "3" * 64,
            "harness_git_commit": "4" * 40,
            "harness_source_root": str(harness),
            "harness_tree_listing_sha256": "5" * 64,
            "harness_tree_oid": "6" * 40,
            "implementation_git_commit": CURRENT_IMPLEMENTATION_COMMIT,
            "implementation_library_path": implementation_library,
            "implementation_library_sha256": "7" * 64,
            "implementation_source_root": str(implementation),
            "implementation_tree_listing_sha256": "8" * 64,
            "implementation_tree_oid": "9" * 40,
            "link_argv": ["/usr/bin/x86_64-linux-gnu-g++-13", role_library],
            "object_closure": [artifact(object_path)],
            "receipt_preimage_sha256": None,
            "role": "current",
            "role_library_path": role_library,
            "role_library_sha256": "7" * 64,
            "runpath": "$ORIGIN",
            "schema": BUILD_RECEIPT_SCHEMA,
            "source_manifest": source_manifest,
            "worker_path": worker_path,
            "worker_sha256": "a" * 64,
        }
        receipt["receipt_preimage_sha256"] = sha256_bytes(canonical_bytes(receipt))
        validate_build_receipt(
            receipt, "current", "4" * 40, CURRENT_IMPLEMENTATION_COMMIT,
            worker_path, "a" * 64, role_library, "7" * 64,
            implementation_library, "7" * 64, str(harness),
            str(harness / controller_relative), sha256_bytes(controller_data))
        mutated = json.loads(canonical_bytes(receipt).decode("ascii"))
        for item in mutated["source_manifest"]:
            if item["repository_relative_path"] == (
                "bench/Wh2V2FacadeTimingWorker.cpp"
            ):
                item["sha256"] = "f" * 64
                break
        mutated["receipt_preimage_sha256"] = None
        mutated["receipt_preimage_sha256"] = sha256_bytes(canonical_bytes(mutated))
        expect_invalid(lambda: validate_build_receipt(
            mutated, "current", "4" * 40, CURRENT_IMPLEMENTATION_COMMIT,
            worker_path, "a" * 64, role_library, "7" * 64,
            implementation_library, "7" * 64, str(harness),
            str(harness / controller_relative), sha256_bytes(controller_data)),
            "build receipt worker-source mutation")

        tracked = base / "tracked-fixture"
        original = b"same-size-a\n"
        write_file(tracked, original)
        info = os.stat(str(tracked), follow_symlinks=False)
        entry = {
            "git_blob_oid": hashlib.sha1(
                "blob {}\0".format(len(original)).encode("ascii") + original
            ).hexdigest(),
            "git_mode": "100644",
            "repository_relative_path": tracked.name,
        }
        verify_tracked_worktree(
            str(base), [entry], "synthetic tracked tree",
            expected_uid=os.getuid(), expected_gid=os.getgid(),
            expected_directory_uid=os.getuid(),
            expected_directory_gid=os.getgid(),
            expected_directory_mode=stat.S_IMODE(os.stat(str(base)).st_mode))
        os.chmod(str(tracked), 0o644)
        fd = os.open(str(tracked), os.O_WRONLY | os.O_TRUNC)
        try:
            replacement = b"same-size-b\n"
            if os.write(fd, replacement) != len(replacement):
                fail("synthetic tracked mutation write was short")
        finally:
            os.close(fd)
        os.chmod(str(tracked), 0o444)
        os.utime(str(tracked), ns=(info.st_atime_ns, info.st_mtime_ns))
        expect_invalid(lambda: verify_tracked_worktree(
            str(base), [entry], "synthetic tracked tree",
            expected_uid=os.getuid(), expected_gid=os.getgid(),
            expected_directory_uid=os.getuid(),
            expected_directory_gid=os.getgid(),
            expected_directory_mode=stat.S_IMODE(os.stat(str(base)).st_mode)),
            "same-size restored-mtime tracked mutation")


def exception_text(exc: BaseException) -> str:
    text = "{}: {}".format(type(exc).__name__, exc)
    return text if len(text) <= 4096 else text[:4093] + "..."


def remaining(deadline: float, where: str) -> float:
    value = deadline - time.monotonic()
    if value <= 0.0:
        fail("internal deadline expired while {}".format(where))
    return value


@dataclass
class WorkerProcess:
    role: str
    path: str
    cpu: int
    process: subprocess.Popen
    stdout_buffer: bytearray = field(default_factory=bytearray)
    stdout_queue: List[bytes] = field(default_factory=list)
    stdout_bytes: int = 0
    stderr: bytearray = field(default_factory=bytearray)
    stderr_overflow: bool = False
    stdout_eof: bool = False
    stderr_eof: bool = False
    invocation_count: int = 0
    contract: Optional[WorkerContract] = None
    terminal_received: bool = False
    release_fd: int = -1
    stub_argv: List[str] = field(default_factory=list)

    @property
    def stdin_fd(self) -> int:
        if self.process.stdin is None:
            fail("{} worker stdin is absent".format(self.role))
        return self.process.stdin.fileno()

    def send(self, command: bytes) -> None:
        if (
            not command or len(command) > 4096 or not command.endswith(b"\n") or
            command.count(b"\n") != 1 or b"\r" in command
        ):
            fail("{} worker command is not one bounded line".format(self.role))
        if self.process.poll() is not None:
            fail("{} worker exited before command".format(self.role))
        offset = 0
        while offset < len(command):
            try:
                written = os.write(self.stdin_fd, command[offset:])
            except BrokenPipeError:
                fail("{} worker closed stdin".format(self.role))
            if written <= 0:
                fail("{} worker command write made no progress".format(self.role))
            offset += written

    def release(self) -> None:
        if self.release_fd < 0:
            fail("{} worker release gate is absent".format(self.role))
        try:
            if os.write(self.release_fd, b"G") != 1:
                fail("{} worker release made no progress".format(self.role))
        finally:
            os.close(self.release_fd)
            self.release_fd = -1

    def terminate(self) -> None:
        if self.release_fd >= 0:
            try:
                os.close(self.release_fd)
            except OSError:
                pass
            self.release_fd = -1
        if self.process.poll() is None:
            try:
                os.killpg(self.process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            except OSError:
                try:
                    self.process.kill()
                except OSError:
                    pass


class WorkerMux:
    def __init__(
        self, workers: Sequence[WorkerProcess], sampler_poll: Optional[Any] = None,
        guardian_poll: Optional[Any] = None,
    ) -> None:
        self.workers = list(workers)
        self.sampler_poll = sampler_poll
        self.guardian_poll = guardian_poll
        self.sampler_event = "none"
        self.guardian_event = "none"
        self.next_authority_poll = 0.0
        self.selector = selectors.DefaultSelector()
        for worker in self.workers:
            if worker.process.stdout is None or worker.process.stderr is None:
                fail("{} worker capture pipes are absent".format(worker.role))
            stdout_fd = worker.process.stdout.fileno()
            stderr_fd = worker.process.stderr.fileno()
            os.set_blocking(stdout_fd, False)
            os.set_blocking(stderr_fd, False)
            self.selector.register(stdout_fd, selectors.EVENT_READ,
                                   (worker, "stdout"))
            self.selector.register(stderr_fd, selectors.EVENT_READ,
                                   (worker, "stderr"))

    def close(self) -> None:
        self.selector.close()

    def _unregister(self, fd: int) -> None:
        try:
            self.selector.unregister(fd)
        except (KeyError, ValueError):
            pass

    def _consume_stdout(self, worker: WorkerProcess, data: bytes) -> None:
        worker.stdout_bytes += len(data)
        if worker.stdout_bytes > MAX_WORKER_STDOUT_BYTES:
            fail("{} worker stdout exceeds its total cap".format(worker.role))
        worker.stdout_buffer.extend(data)
        while True:
            newline = worker.stdout_buffer.find(b"\n")
            if newline < 0:
                if len(worker.stdout_buffer) > MAX_JSON_LINE_BYTES:
                    fail("{} worker stdout line exceeds its cap".format(worker.role))
                break
            line = bytes(worker.stdout_buffer[:newline + 1])
            del worker.stdout_buffer[:newline + 1]
            if len(line) > MAX_JSON_LINE_BYTES or b"\r" in line:
                fail("{} worker emitted a noncanonical line".format(worker.role))
            worker.stdout_queue.append(line)
            if len(worker.stdout_queue) > 4:
                fail("{} worker emitted unsolicited records".format(worker.role))

    def _read_ready(self, worker: WorkerProcess, stream: str, fd: int) -> None:
        try:
            data = os.read(fd, 65536)
        except BlockingIOError:
            return
        if not data:
            self._unregister(fd)
            if stream == "stdout":
                worker.stdout_eof = True
                if worker.stdout_buffer:
                    fail("{} worker stdout ended mid-line".format(worker.role))
            else:
                worker.stderr_eof = True
            return
        if stream == "stdout":
            self._consume_stdout(worker, data)
        else:
            available = MAX_WORKER_STDERR_BYTES - len(worker.stderr)
            if available > 0:
                worker.stderr.extend(data[:available])
            if len(data) > available:
                worker.stderr_overflow = True
                fail("{} worker stderr exceeds its cap".format(worker.role))

    def check_authorities(self, force: bool = False) -> None:
        if self.guardian_event != "none":
            fail("outer guardian supervision: " + self.guardian_event)
        if self.sampler_event != "none":
            fail("sampler supervision: " + self.sampler_event)
        now = time.monotonic()
        if not force and now < self.next_authority_poll:
            return
        self.next_authority_poll = now + 0.05
        event = "none"
        if self.guardian_poll is not None:
            try:
                event = self.guardian_poll()
            except BaseException as exc:
                event = "guardian-monitor-invalid:" + exception_text(exc)
            if event != "none":
                self.guardian_event = event
        if event == "none" and self.sampler_poll is not None:
            try:
                event = self.sampler_poll()
            except BaseException as exc:
                event = "sampler-monitor-invalid:" + exception_text(exc)
            if event != "none":
                self.sampler_event = event
        if event != "none":
            for worker in self.workers:
                worker.terminate()
            if self.guardian_event != "none":
                fail("outer guardian supervision: " + event)
            fail("sampler supervision: " + event)

    # Keep the protocol-facing name so every invocation/terminal check uses the
    # combined sampler + independent-guardian authority poll.
    def check_sampler(self, force: bool = False) -> None:
        self.check_authorities(force)

    def pump(self, deadline: float, where: str) -> None:
        self.check_authorities()
        poll_slice = 0.05 if (
            self.sampler_poll is not None or self.guardian_poll is not None
        ) else 1.0
        events = self.selector.select(min(poll_slice, remaining(deadline, where)))
        self.check_authorities()
        if not events:
            for worker in self.workers:
                code = worker.process.poll()
                if (
                    code is not None and not worker.stdout_queue and
                    not worker.terminal_received
                ):
                    fail("{} worker exited early with {}".format(worker.role, code))
            return
        for key, _mask in events:
            worker, stream = key.data
            self._read_ready(worker, stream, key.fd)

    def read_line(
        self, worker: WorkerProcess, deadline: float, where: str,
    ) -> bytes:
        while not worker.stdout_queue:
            if worker.stdout_eof:
                fail("{} worker stdout closed while {}".format(worker.role, where))
            self.pump(deadline, where)
        self.check_authorities()
        return worker.stdout_queue.pop(0)

    def require_quiet(self) -> None:
        self.check_authorities()
        for worker in self.workers:
            if worker.stdout_queue:
                fail("{} worker emitted an unsolicited record".format(worker.role))

    def drain_and_reap(self, deadline: float) -> None:
        while (
            any(worker.process.poll() is None for worker in self.workers) or
            bool(self.selector.get_map())
        ):
            self.pump(deadline, "reaping workers")
        for worker in self.workers:
            if worker.process.returncode != 0:
                fail("{} worker exited with {}".format(
                    worker.role, worker.process.returncode))
            if worker.stdout_queue or worker.stdout_buffer:
                fail("{} worker left extra stdout".format(worker.role))

    def contain_and_drain(self, timeout: float = 5.0) -> List[str]:
        """Kill all worker PGs, retain stderr through EOF, and reap boundedly."""
        errors: List[str] = []
        deadline = time.monotonic() + timeout
        for worker in self.workers:
            worker.terminate()
            if worker.process.stdin is not None:
                try:
                    worker.process.stdin.close()
                except OSError as exc:
                    errors.append("{} stdin close: {}".format(
                        worker.role, exception_text(exc)))
        while self.selector.get_map() and time.monotonic() < deadline:
            events = self.selector.select(min(0.05, deadline - time.monotonic()))
            for key, _mask in events:
                worker, stream = key.data
                try:
                    data = os.read(key.fd, 65536)
                except BlockingIOError:
                    continue
                except OSError as exc:
                    errors.append("{} {} drain: {}".format(
                        worker.role, stream, exception_text(exc)))
                    self._unregister(key.fd)
                    continue
                if not data:
                    self._unregister(key.fd)
                    if stream == "stdout":
                        worker.stdout_eof = True
                    else:
                        worker.stderr_eof = True
                    continue
                if stream == "stderr":
                    available = MAX_WORKER_STDERR_BYTES - len(worker.stderr)
                    if available > 0:
                        worker.stderr.extend(data[:available])
                    if len(data) > available:
                        worker.stderr_overflow = True
                else:
                    worker.stdout_bytes += len(data)
        if self.selector.get_map():
            errors.append("worker containment pipes did not reach EOF")
        for worker in self.workers:
            remaining_wait = max(0.0, deadline - time.monotonic())
            try:
                worker.process.wait(timeout=remaining_wait)
            except BaseException as exc:
                errors.append("{} containment reap: {}".format(
                    worker.role, exception_text(exc)))
        return errors


WORKER_STUB_CODE = r'''import os,sys
fd=int(sys.argv[1]); cpu=int(sys.argv[2]); target=sys.argv[3]; argv=sys.argv[3:]
os.sched_setaffinity(0,{cpu})
token=os.read(fd,1); os.close(fd)
if token != b"G": sys.exit(4)
env={"LANG":"C.UTF-8","LC_ALL":"C.UTF-8","PATH":"/usr/bin:/bin","TZ":"UTC"}
os.execve(target,argv,env)
'''


def launch_worker(role: str, path: str, cpu: int) -> WorkerProcess:
    worker_argv = [path, "--serve", "--role", role, "--cpu", str(cpu)]
    environment = dict(SEALED_ENV)
    release_read, release_write = os.pipe2(getattr(os, "O_CLOEXEC", 0))
    argv = [
        "/usr/bin/python3", "-I", "-S", "-c", WORKER_STUB_CODE,
        str(release_read), str(cpu),
    ] + worker_argv
    try:
        process = subprocess.Popen(
            argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, close_fds=True, start_new_session=True,
            env=environment, bufsize=0, pass_fds=(release_read,),
        )
    except BaseException as exc:
        os.close(release_read)
        os.close(release_write)
        fail("could not launch {} worker: {}".format(role, exc))
    os.close(release_read)
    return WorkerProcess(
        role=role, path=path, cpu=cpu, process=process,
        release_fd=release_write, stub_argv=argv)


def validate_worker_affinity(worker: WorkerProcess) -> None:
    try:
        affinity = os.sched_getaffinity(worker.process.pid)
    except (AttributeError, OSError) as exc:
        fail("could not read {} worker affinity: {}".format(worker.role, exc))
    if affinity != {worker.cpu}:
        fail("{} worker affinity is not singleton CPU {}".format(
            worker.role, worker.cpu))


def invoke_worker(
    mux: WorkerMux, worker: WorkerProcess, seq: int, arm: str, scope: str,
    k: int, block_bytes: int, deadline: float,
) -> Tuple[Dict[str, Any], bytes]:
    if worker.contract is None:
        fail("{} worker has no validated contract".format(worker.role))
    if ARM_ROLE.get(arm) != worker.role:
        fail("arm {} was routed to the wrong role".format(arm))
    mux.check_sampler()
    mux.require_quiet()
    command = "invoke {} {} {} {} {}\n".format(
        seq, arm, scope, k, block_bytes).encode("ascii")
    worker.send(command)
    line = mux.read_line(worker, deadline, "reading invocation {}".format(seq))
    mux.require_quiet()
    value = parse_canonical_line(line, "{} invocation {}".format(worker.role, seq))
    validate_invocation(
        value, worker.contract, seq, arm, scope, k, block_bytes)
    worker.invocation_count += 1
    return value, line


@dataclass
class SlotAccumulator:
    side: str
    arm: str
    first_sequence: int
    last_sequence: int
    invocation_count: int = 0
    elapsed_ns: int = 0
    constructor_ns: int = 0
    systematic_ns: int = 0
    init_first_repair_ns: int = 0
    borrowed_systematic_invocations: int = 0
    digest: Any = field(default_factory=hashlib.sha256)

    def add(self, value: Mapping[str, Any], line: bytes) -> None:
        self.invocation_count += 1
        self.elapsed_ns += value["elapsed_ns"]
        self.constructor_ns += value["constructor_ns"] or 0
        self.systematic_ns += value["systematic_ns"] or 0
        self.init_first_repair_ns += value["init_first_repair_ns"] or 0
        self.borrowed_systematic_invocations += value[
            "borrowed_systematic_invocations"]
        for name, total in (
            ("elapsed", self.elapsed_ns),
            ("constructor", self.constructor_ns),
            ("systematic", self.systematic_ns),
            ("init-first-repair", self.init_first_repair_ns),
            ("borrowed-systematic", self.borrowed_systematic_invocations),
        ):
            if total < 0 or total > MAX_INT63:
                fail("slot {} sum overflows".format(name))
        self.digest.update(line)

    def record(self) -> Dict[str, Any]:
        return {
            "arm": self.arm,
            "borrowed_systematic_invocations":
                self.borrowed_systematic_invocations,
            "constructor_ns": self.constructor_ns,
            "elapsed_ns": self.elapsed_ns,
            "first_sequence": self.first_sequence,
            "init_first_repair_ns": self.init_first_repair_ns,
            "invocation_count": self.invocation_count,
            "last_sequence": self.last_sequence,
            "response_stream_sha256": self.digest.hexdigest(),
            "side": self.side,
            "systematic_ns": self.systematic_ns,
        }


def run_panel(
    mux: WorkerMux, workers: Mapping[str, WorkerProcess], seq: int,
    replicate: int, k: int, block_bytes: int, scope: str, comparison: str,
    deadline: float,
) -> Tuple[Dict[str, Any], int]:
    left_arm, right_arm = COMPARISON_ARMS[comparison]
    order = panel_order(replicate, k, block_bytes, scope, comparison)
    n = invocations_by_replicate(k)[replicate]
    first_sequence = seq
    warmup_digest = hashlib.sha256()
    for arm in (left_arm, right_arm):
        worker = workers[ARM_ROLE[arm]]
        _value, line = invoke_worker(
            mux, worker, seq, arm, scope, k, block_bytes, deadline)
        warmup_digest.update(line)
        seq += 1

    sides = sides_for(order)
    slots: List[SlotAccumulator] = []
    for slot_index, side in enumerate(sides):
        arm = left_arm if side == "left" else right_arm
        sequence_first, sequence_last = expected_panel_slot_sequence(
            first_sequence, n, slot_index)
        slots.append(SlotAccumulator(
            side=side, arm=arm, first_sequence=sequence_first,
            last_sequence=sequence_last))
    for block_offset in (0, 4):
        for _repeat in range(n):
            for slot_index in range(block_offset, block_offset + 4):
                slot = slots[slot_index]
                worker = workers[ARM_ROLE[slot.arm]]
                value, line = invoke_worker(
                    mux, worker, seq, slot.arm, scope, k, block_bytes,
                    deadline)
                slot.add(value, line)
                seq += 1
    record = {
        "K": k,
        "block_bytes": block_bytes,
        "comparison": comparison,
        "first_sequence": first_sequence,
        "invocations_per_slot": n,
        "last_sequence": seq - 1,
        "left_arm": left_arm,
        "order": order,
        "panel_key_sha256": panel_key_sha256(
            replicate, k, block_bytes, scope, comparison),
        "replicate": replicate,
        "right_arm": right_arm,
        "schema": PANEL_SCHEMA,
        "scope": scope,
        "slots": [slot.record() for slot in slots],
        "warmup_response_stream_sha256": warmup_digest.hexdigest(),
        "warmup_sequences": [first_sequence, first_sequence + 1],
    }
    validate_panel(
        record, replicate, k, block_bytes, scope, comparison, first_sequence)
    return record, seq


@dataclass
class OutputFile:
    name: str
    fd: int
    bytes_written: int = 0
    digest: Any = field(default_factory=hashlib.sha256)
    closed: bool = False

    def append(self, data: bytes, maximum: int) -> None:
        if self.closed or self.bytes_written + len(data) > maximum:
            fail("output {} exceeds its write policy".format(self.name))
        offset = 0
        while offset < len(data):
            written = os.write(self.fd, data[offset:])
            if written <= 0:
                fail("output {} write made no progress".format(self.name))
            offset += written
        self.bytes_written += len(data)
        self.digest.update(data)

    def receipt(self) -> Dict[str, Any]:
        return {
            "bytes": self.bytes_written,
            "name": self.name,
            "sha256": self.digest.hexdigest(),
        }


OUTPUT_PARENT_KEYS = {"device", "gid", "inode", "mode", "uid"}
EMERGENCY_OUTPUT_SCHEMA = \
    "wirehair.wh2.v2-facade-timing-screen.emergency-invalid.v1"


def output_parent_identity(info: os.stat_result) -> Dict[str, int]:
    return {
        "device": info.st_dev,
        "gid": info.st_gid,
        "inode": info.st_ino,
        "mode": stat.S_IMODE(info.st_mode),
        "uid": info.st_uid,
    }


def expected_output_parent(args: argparse.Namespace) -> Dict[str, int]:
    value = {
        "device": args.expected_output_parent_device,
        "gid": args.expected_output_parent_gid,
        "inode": args.expected_output_parent_inode,
        "mode": args.expected_output_parent_mode,
        "uid": args.expected_output_parent_uid,
    }
    exact_keys(value, OUTPUT_PARENT_KEYS, "expected output parent")
    for name in ("device", "gid", "inode", "uid"):
        lower = 1 if name == "inode" else 0
        exact_int(value[name], lower, MAX_UINT64,
                  "expected output parent " + name)
    if value["uid"] != 1000 or value["gid"] != 1000:
        fail("expected output parent owner is not UID/GID1000")
    mode = exact_int(value["mode"], 0, 0o777,
                     "expected output parent mode")
    if mode & 0o700 != 0o700 or mode & 0o022:
        fail("expected output parent mode is not protected")
    return value


class OutputBundle:
    def __init__(self, path: str, expected_parent: Mapping[str, Any]) -> None:
        exact_absolute_path(path, "output directory")
        if Path(path).name != CAMPAIGN:
            fail("output directory basename differs from the fixed namespace")
        parent = str(Path(path).parent)
        exact_keys(expected_parent, OUTPUT_PARENT_KEYS, "expected output parent")
        for name in ("device", "gid", "inode", "uid"):
            exact_int(expected_parent[name], 0, MAX_UINT64,
                      "expected output parent " + name)
        expected_mode = exact_int(
            expected_parent["mode"], 0, 0o777,
            "expected output parent mode")
        if expected_mode & 0o700 != 0o700 or expected_mode & 0o022:
            fail("expected output parent is not owner-writable protected storage")
        self.path = path
        self.parent_path = parent
        self.parent_fd = -1
        self.directory_fd = -1
        self.files: Dict[str, OutputFile] = {}
        self.sealed = False
        self.created = False
        self.parent_receipt: Dict[str, Any] = {}
        self.directory_receipt: Dict[str, Any] = {}
        try:
            self.parent_fd = os.open(
                parent, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) |
                getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0))
            parent_info = os.fstat(self.parent_fd)
            parent_path_info = os.stat(parent, follow_symlinks=False)
            if (
                not stat.S_ISDIR(parent_info.st_mode) or
                output_parent_identity(parent_info) != dict(expected_parent) or
                output_parent_identity(parent_path_info) != dict(expected_parent)
            ):
                fail("output parent does not match its external identity")
            self.parent_receipt = output_parent_identity(parent_info)
            try:
                os.mkdir(CAMPAIGN, 0o700, dir_fd=self.parent_fd)
            except OSError as exc:
                fail("could not exclusively reserve output directory: {}".format(exc))
            self.created = True
            self.directory_fd = os.open(
                CAMPAIGN, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) |
                getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0),
                dir_fd=self.parent_fd)
            directory_info = os.fstat(self.directory_fd)
            directory_path_info = os.stat(
                CAMPAIGN, dir_fd=self.parent_fd, follow_symlinks=False)
            if (
                not stat.S_ISDIR(directory_info.st_mode) or
                stat.S_IMODE(directory_info.st_mode) != 0o700 or
                stat_receipt(directory_info) != stat_receipt(directory_path_info)
            ):
                fail("reserved output directory policy differs")
            self.directory_receipt = {
                "device": directory_info.st_dev,
                "gid": directory_info.st_gid,
                "inode": directory_info.st_ino,
                "mode": stat.S_IMODE(directory_info.st_mode),
                "uid": directory_info.st_uid,
            }
            for name in EXPECTED_OUTPUT_NAMES:
                fd = os.open(
                    name, os.O_RDWR | os.O_CREAT | os.O_EXCL |
                    getattr(os, "O_CLOEXEC", 0) |
                    getattr(os, "O_NOFOLLOW", 0),
                    0o600, dir_fd=self.directory_fd)
                self.files[name] = OutputFile(name=name, fd=fd)
            os.fsync(self.directory_fd)
            self.verify_parent()
        except BaseException as exc:
            if self.created:
                try:
                    if self.directory_fd < 0 and self.parent_fd >= 0:
                        self.directory_fd = os.open(
                            CAMPAIGN,
                            os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) |
                            getattr(os, "O_CLOEXEC", 0) |
                            getattr(os, "O_NOFOLLOW", 0),
                            dir_fd=self.parent_fd)
                    self.emergency_invalid(
                        "output reservation failed: " + exception_text(exc))
                except BaseException:
                    # The root launcher owns the externally reserved parent and
                    # treats a consumed namespace without a valid COMPLETE as
                    # invalid if even emergency terminalization is impossible.
                    pass
            self.close()
            raise

    def verify_parent(self) -> None:
        if self.parent_fd < 0:
            fail("output parent descriptor is closed")
        descriptor = os.fstat(self.parent_fd)
        path_info = os.stat(self.parent_path, follow_symlinks=False)
        if (
            output_parent_identity(descriptor) != self.parent_receipt or
            output_parent_identity(path_info) != self.parent_receipt
        ):
            fail("output parent identity changed")

    def receipt(self) -> Dict[str, Any]:
        self.verify_parent()
        directory = os.fstat(self.directory_fd)
        if (
            directory.st_dev != self.directory_receipt["device"] or
            directory.st_ino != self.directory_receipt["inode"] or
            directory.st_uid != self.directory_receipt["uid"] or
            directory.st_gid != self.directory_receipt["gid"]
        ):
            fail("reserved output directory identity changed")
        return {
            "directory": dict(self.directory_receipt),
            "parent": dict(self.parent_receipt),
            "path": self.path,
        }

    def append_raw(self, record: Mapping[str, Any]) -> None:
        self.files["raw.jsonl"].append(canonical_bytes(record) + b"\n", MAX_RAW_BYTES)

    def write_once(self, name: str, value: Mapping[str, Any], maximum: int) -> None:
        output = self.files[name]
        if output.bytes_written != 0:
            fail("output {} was already written".format(name))
        output.append(canonical_bytes(value) + b"\n", maximum)

    def raw_bytes(self) -> bytes:
        output = self.files["raw.jsonl"]
        os.fsync(output.fd)
        data = bytearray()
        offset = 0
        while offset < output.bytes_written:
            block = os.pread(
                output.fd, min(1024 * 1024, output.bytes_written - offset), offset)
            if not block:
                fail("short read while replaying raw output")
            data.extend(block)
            offset += len(block)
        return bytes(data)

    def _verify_named_output(
        self, output: OutputFile, expected_mode: int,
    ) -> Dict[str, Any]:
        descriptor = os.fstat(output.fd)
        named = os.stat(
            output.name, dir_fd=self.directory_fd, follow_symlinks=False)
        readback_sha256 = file_sha256_fd(output.fd, descriptor.st_size)
        descriptor_after = os.fstat(output.fd)
        named_after = os.stat(
            output.name, dir_fd=self.directory_fd, follow_symlinks=False)
        if (
            not stat.S_ISREG(descriptor.st_mode) or descriptor.st_nlink != 1 or
            descriptor.st_uid != 1000 or descriptor.st_gid != 1000 or
            stat.S_IMODE(descriptor.st_mode) != expected_mode or
            stat_receipt(descriptor) != stat_receipt(named) or
            descriptor.st_size != output.bytes_written or
            readback_sha256 != output.digest.hexdigest() or
            stat_receipt(descriptor_after) != stat_receipt(descriptor) or
            stat_receipt(named_after) != stat_receipt(descriptor)
        ):
            fail("output {} named-file/readback seal differs".format(output.name))
        return output.receipt()

    def seal(self, outcome: str) -> Dict[str, Any]:
        if self.sealed:
            fail("output bundle was already sealed")
        self.verify_parent()
        receipts: List[Dict[str, Any]] = []
        for name in EXPECTED_OUTPUT_NAMES:
            output = self.files[name]
            os.fsync(output.fd)
            os.fchmod(output.fd, 0o400)
            receipts.append(self._verify_named_output(output, 0o400))
        self.verify_parent()
        if sorted(os.listdir(self.directory_fd)) != sorted(EXPECTED_OUTPUT_NAMES):
            fail("output namespace differs before COMPLETE")
        complete: Dict[str, Any] = {
            "campaign": CAMPAIGN,
            "files": receipts,
            "outcome": outcome,
            "preimage_sha256": None,
            "schema": COMPLETE_SCHEMA,
        }
        complete["preimage_sha256"] = sha256_bytes(canonical_bytes(complete))
        complete_data = canonical_bytes(complete) + b"\n"
        complete_fd = os.open(
            "COMPLETE", os.O_RDWR | os.O_CREAT | os.O_EXCL |
            getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0),
            0o400, dir_fd=self.directory_fd)
        try:
            offset = 0
            while offset < len(complete_data):
                written = os.write(complete_fd, complete_data[offset:])
                if written <= 0:
                    fail("COMPLETE marker write made no progress")
                offset += written
            os.fsync(complete_fd)
            if stat.S_IMODE(os.fstat(complete_fd).st_mode) != 0o400:
                fail("COMPLETE marker mode differs")
            complete_info = os.fstat(complete_fd)
            complete_named = os.stat(
                "COMPLETE", dir_fd=self.directory_fd, follow_symlinks=False)
            if (
                complete_info.st_nlink != 1 or complete_info.st_uid != 1000 or
                complete_info.st_gid != 1000 or
                stat_receipt(complete_info) != stat_receipt(complete_named) or
                complete_info.st_size != len(complete_data) or
                file_sha256_fd(complete_fd, complete_info.st_size) !=
                    sha256_bytes(complete_data)
            ):
                fail("COMPLETE named-file/readback seal differs")
        finally:
            os.close(complete_fd)
        os.fsync(self.directory_fd)
        os.fchmod(self.directory_fd, 0o500)
        os.fsync(self.directory_fd)
        self.sealed = True
        return complete

    def _actual_file_receipt(self, name: str, fd: int) -> Dict[str, Any]:
        info = os.fstat(fd)
        if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
            fail("emergency output {} policy differs".format(name))
        return {
            "bytes": info.st_size,
            "name": name,
            "sha256": file_sha256_fd(fd, info.st_size),
        }

    def emergency_invalid(self, error: str) -> None:
        """Best-effort terminalization after this namespace was consumed."""
        if self.sealed or self.directory_fd < 0:
            return
        if type(error) is not str or not error:
            error = "unspecified emergency publication failure"
        error = error[:4096]
        for name in EXPECTED_OUTPUT_NAMES:
            if name not in self.files:
                flags = (os.O_RDWR | getattr(os, "O_CLOEXEC", 0) |
                         getattr(os, "O_NOFOLLOW", 0))
                try:
                    fd = os.open(name, flags, dir_fd=self.directory_fd)
                except FileNotFoundError:
                    fd = os.open(
                        name, flags | os.O_CREAT | os.O_EXCL, 0o600,
                        dir_fd=self.directory_fd)
                self.files[name] = OutputFile(name=name, fd=fd)
            output = self.files[name]
            info = os.fstat(output.fd)
            if info.st_size == 0 and name in (
                "raw.jsonl", "provenance.json", "summary.json",
            ):
                output.append(canonical_bytes({
                    "campaign": CAMPAIGN,
                    "error": error,
                    "outcome": "invalid",
                    "schema": EMERGENCY_OUTPUT_SCHEMA,
                }) + b"\n", MAX_RECEIPT_BYTES)
            os.fsync(output.fd)
            os.fchmod(output.fd, 0o400)
        receipts = []
        for name in EXPECTED_OUTPUT_NAMES:
            output = self.files[name]
            actual = self._actual_file_receipt(name, output.fd)
            output.bytes_written = actual["bytes"]
            output.digest = hashlib.sha256()
            offset = 0
            while offset < actual["bytes"]:
                block = os.pread(
                    output.fd, min(1024 * 1024, actual["bytes"] - offset), offset)
                if not block:
                    fail("emergency output {} readback is short".format(name))
                output.digest.update(block)
                offset += len(block)
            receipts.append(self._verify_named_output(output, 0o400))
        self.verify_parent()
        if sorted(os.listdir(self.directory_fd)) != sorted(EXPECTED_OUTPUT_NAMES):
            fail("emergency output namespace differs before COMPLETE")
        complete: Dict[str, Any] = {
            "campaign": CAMPAIGN,
            "files": receipts,
            "outcome": "invalid",
            "preimage_sha256": None,
            "schema": COMPLETE_SCHEMA,
        }
        complete["preimage_sha256"] = sha256_bytes(canonical_bytes(complete))
        data = canonical_bytes(complete) + b"\n"
        fd = os.open(
            "COMPLETE", os.O_RDWR | os.O_CREAT | os.O_EXCL |
            getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0),
            0o400, dir_fd=self.directory_fd)
        try:
            offset = 0
            while offset < len(data):
                written = os.write(fd, data[offset:])
                if written <= 0:
                    fail("emergency COMPLETE write made no progress")
                offset += written
            os.fsync(fd)
        finally:
            os.close(fd)
        os.fsync(self.directory_fd)
        os.fchmod(self.directory_fd, 0o500)
        os.fsync(self.directory_fd)
        self.sealed = True

    def close(self) -> None:
        for output in self.files.values():
            if not output.closed:
                try:
                    os.close(output.fd)
                except OSError:
                    pass
                output.closed = True
        if self.directory_fd >= 0:
            try:
                os.close(self.directory_fd)
            except OSError:
                pass
            self.directory_fd = -1
        if self.parent_fd >= 0:
            try:
                os.close(self.parent_fd)
            except OSError:
                pass
            self.parent_fd = -1


def make_self_hashed(value: Dict[str, Any], field_name: str) -> Dict[str, Any]:
    if field_name not in value or value[field_name] is not None:
        fail("self-hash placeholder differs")
    value[field_name] = sha256_bytes(canonical_bytes(value))
    return value


def make_summary(
    outcome: str, raw: OutputFile, raw_complete: bool,
    statistics: Optional[Mapping[str, Any]], failures: Sequence[str],
    provenance_sha256: str, worker_terminals: Optional[Mapping[str, Any]],
    health: Optional[Mapping[str, Any]],
) -> Dict[str, Any]:
    if outcome not in ("pass", "reject", "invalid"):
        fail("summary outcome differs")
    value: Dict[str, Any] = {
        "campaign": CAMPAIGN,
        "failures": list(failures)[:MAX_FAILURES],
        "health": health,
        "outcome": outcome,
        "provenance_sha256": provenance_sha256,
        "raw_bytes": raw.bytes_written,
        "raw_complete": raw_complete,
        "raw_sha256": raw.digest.hexdigest(),
        "schema": SUMMARY_SCHEMA,
        "statistics": statistics,
        "summary_preimage_sha256": None,
        "valid_evidence": outcome != "invalid",
        "worker_terminals": worker_terminals,
    }
    return make_self_hashed(value, "summary_preimage_sha256")


def make_provenance(
    args: argparse.Namespace, artifacts: Mapping[str, OpenArtifact],
    build_receipts: Mapping[str, Mapping[str, Any]], start_ns: int,
    end_ns: int, error: str, closure_sha256: Mapping[str, str],
    health_adapter: HealthAdapter, health: Optional[Mapping[str, Any]],
    guardian: Optional[Mapping[str, Any]],
    worker_processes: Optional[Mapping[str, Any]],
    output_reservation: Mapping[str, Any],
) -> Dict[str, Any]:
    value: Dict[str, Any] = {
        "artifacts": [artifacts[name].receipt() for name in sorted(artifacts)],
        "build_receipts": {
            role: build_receipts[role] for role in ("current", "parent")
        },
        "campaign": CAMPAIGN,
        "closure_verification_sha256": dict(closure_sha256),
        "controller_end_monotonic_ns": end_ns,
        "controller_security_boundary": args.controller_security_boundary,
        "controller_start_monotonic_ns": start_ns,
        "error": error,
        "expected_current_implementation_commit":
            args.expected_current_implementation_commit,
        "expected_controller_sha256": args.expected_controller_sha256,
        "expected_harness_commit": args.expected_harness_commit,
        "expected_parent_implementation_commit":
            args.expected_parent_implementation_commit,
        "internal_deadline_seconds": INTERNAL_DEADLINE_SECONDS,
        "health": health,
        "health_adapter_sha256": health_adapter.adapter_sha256,
        "health_module_loader": health_adapter.modules.receipt,
        "health_source_manifest": health_adapter.source_manifest,
        "health_source_git_receipt": health_adapter.git_receipt,
        "guardian": guardian,
        "launch_argv": {
            "current": [args.current_worker, "--serve", "--role", "current",
                        "--cpu", str(args.cpu)],
            "parent": [args.parent_worker, "--serve", "--role", "parent",
                       "--cpu", str(args.cpu)],
        },
        "preimage_sha256": None,
        "schema": PROVENANCE_SCHEMA,
        "sealed_environment": dict(SEALED_ENV),
        "target_cpu": args.cpu,
        "trusted_security_assumption": (
            "root-owned-0555 source directories and no concurrent hostile "
            "UID1000 process; root launcher independently journals outputs"
        ),
        "worker_processes": worker_processes,
        "outer_guardian": {
            "code_sha256": sha256_bytes(GUARDIAN_CODE.encode("ascii")),
            "deadline_seconds": EXTERNAL_DEADLINE_SECONDS,
            "python": "/usr/bin/python3",
        },
        "output_reservation": dict(output_reservation),
    }
    return make_self_hashed(value, "preimage_sha256")


@dataclass(frozen=True)
class HealthAdapter:
    module: Any
    modules: Any
    source_manifest: Mapping[str, Any]
    adapter_sha256: str
    git_receipt: Mapping[str, Any]


def load_health_adapter(
    source_root: Path, expected_manifest_sha256: str,
    expected_adapter_sha256: str,
    expected_harness_commit: str, expected_git_sha256: str,
) -> HealthAdapter:
    if not sys.flags.isolated or not sys.dont_write_bytecode or sys.flags.optimize != 0:
        fail("scientific run requires unoptimized Python -I -B")
    adapter_path = source_root / HEALTH_ADAPTER_PATH
    data = adapter_path.read_bytes()
    lower_hash(expected_adapter_sha256, "expected health adapter SHA-256")
    if sha256_bytes(data) != expected_adapter_sha256:
        fail("health adapter bytes differ before execution")
    module_name = "_wh2_v2_facade_sealed_health_adapter"
    if module_name in sys.modules:
        fail("health adapter was already loaded")
    module = types.ModuleType(module_name)
    module.__file__ = str(adapter_path)
    module.__package__ = ""
    sys.modules[module_name] = module
    try:
        code = compile(data, str(adapter_path), "exec", dont_inherit=True, optimize=0)
        exec(code, module.__dict__)
        manifest = module.source_manifest(source_root)
        if manifest["sha256"] != expected_manifest_sha256:
            fail("health source manifest differs from its presealed value")
        modules = module.load_sealed_health_modules(source_root, manifest)
        git = module.git_receipt(
            source_root, expected_harness_commit, expected_git_sha256)
        return HealthAdapter(
            module=module, modules=modules, source_manifest=manifest,
            adapter_sha256=sha256_bytes(data), git_receipt=git)
    except BaseException:
        if sys.modules.get(module_name) is module:
            del sys.modules[module_name]
        raise


GUARDIAN_CODE = r'''import os,select,signal,sys,time
fd=int(sys.argv[1]); ready=int(sys.argv[2]); controller=int(sys.argv[3]); deadline_ns=int(sys.argv[4]); cpu=int(sys.argv[5])
os.sched_setaffinity(0,{cpu}); anchor_ns=time.monotonic_ns(); os.write(ready,(str(anchor_ns)+" "+str(deadline_ns)+"\n").encode("ascii")); os.close(ready)
pids=[]; pending=b""; done=False
while not done:
    left=(deadline_ns-time.monotonic_ns())/1000000000.0
    if left <= 0: break
    ready,_,_=select.select([fd],[],[],left)
    if not ready: break
    data=os.read(fd,4096)
    if not data: break
    pending += data
    while b"\n" in pending:
        line,pending=pending.split(b"\n",1)
        parts=line.split()
        if parts == [b"done"]: done=True; break
        if len(parts)==2 and parts[0]==b"pid" and parts[1].isdigit():
            pid=int(parts[1])
            if pid > 1 and pid not in pids: pids.append(pid)
        else: pending=b"invalid"; break
    if pending == b"invalid": break
if done: sys.exit(0)
for pid in pids:
    try: os.killpg(pid, signal.SIGKILL)
    except ProcessLookupError: pass
    except OSError: pass
try: os.kill(controller, signal.SIGKILL)
except ProcessLookupError: pass
except OSError: pass
sys.exit(3)
'''


def process_start_ticks(pid: int) -> int:
    try:
        data = Path("/proc/{}/stat".format(pid)).read_bytes()
    except OSError as exc:
        fail("cannot read process {} start identity: {}".format(pid, exc))
    marker = data.rfind(b") ")
    if marker < 0:
        fail("process {} stat comm field is malformed".format(pid))
    fields = data[marker + 2:].split()
    if len(fields) < 20 or not fields[19].isdigit():
        fail("process {} stat start field is malformed".format(pid))
    return exact_int(int(fields[19]), 1, MAX_UINT64,
                     "process {} start ticks".format(pid))


def bounded_proc_file(pid: int, name: str, maximum: int) -> bytes:
    try:
        with open("/proc/{}/{}".format(pid, name), "rb", buffering=0) as source:
            data = source.read(maximum + 1)
    except OSError as exc:
        fail("cannot read process {} {}: {}".format(pid, name, exc))
    if len(data) > maximum:
        fail("process {} {} exceeds its cap".format(pid, name))
    return data


def exact_nul_vector(data: bytes, where: str) -> List[bytes]:
    if not data or not data.endswith(b"\0") or b"\0\0" in data:
        fail("{} is not a canonical NUL vector".format(where))
    return data[:-1].split(b"\0")


def process_security_receipt(
    pid: int, expected_affinity: Sequence[int],
    expected_executable_stat: Mapping[str, Any],
    expected_cmdline: Optional[Sequence[str]] = None,
) -> Dict[str, Any]:
    status_data = bounded_proc_file(pid, "status", 1024 * 1024)
    status_fields: Dict[bytes, bytes] = {}
    wanted = {
        b"Uid", b"Gid", b"Groups", b"CapInh", b"CapPrm", b"CapEff",
        b"CapBnd", b"CapAmb", b"NoNewPrivs",
    }
    for line in status_data.splitlines():
        if b":" not in line:
            continue
        key, value = line.split(b":", 1)
        if key in wanted:
            if key in status_fields:
                fail("process {} status repeats {}".format(
                    pid, key.decode("ascii")))
            status_fields[key] = value.strip()
    if set(status_fields) != wanted:
        fail("process {} security status is incomplete".format(pid))
    if status_fields[b"Uid"].split() != [b"1000"] * 4:
        fail("process {} UID boundary differs".format(pid))
    if status_fields[b"Gid"].split() != [b"1000"] * 4:
        fail("process {} GID boundary differs".format(pid))
    if status_fields[b"Groups"].split():
        fail("process {} supplementary groups are not empty".format(pid))
    for name in (b"CapInh", b"CapPrm", b"CapEff", b"CapBnd", b"CapAmb"):
        try:
            capability = int(status_fields[name], 16)
        except ValueError:
            fail("process {} {} is malformed".format(pid, name.decode("ascii")))
        if capability != 0:
            fail("process {} {} is nonzero".format(pid, name.decode("ascii")))
    if status_fields[b"NoNewPrivs"] != b"1":
        fail("process {} lacks no-new-privileges".format(pid))
    environment = bounded_proc_file(pid, "environ", 64 * 1024)
    expected_environment = b"".join(
        (key + "=" + value + "\0").encode("ascii")
        for key, value in SEALED_ENV.items())
    if environment != expected_environment:
        fail("process {} environment differs from the exact four keys".format(pid))
    cmdline_data = bounded_proc_file(pid, "cmdline", 1024 * 1024)
    cmdline = [item.decode("utf-8") for item in exact_nul_vector(
        cmdline_data, "process {} cmdline".format(pid))]
    if expected_cmdline is not None and cmdline != list(expected_cmdline):
        fail("process {} command line differs".format(pid))
    affinity = sorted(os.sched_getaffinity(pid))
    if affinity != list(expected_affinity):
        fail("process {} affinity differs".format(pid))
    executable_info = os.stat("/proc/{}/exe".format(pid))
    executable_receipt = stat_receipt(executable_info)
    if executable_receipt != dict(expected_executable_stat):
        fail("process {} executable image differs".format(pid))
    return {
        "affinity": affinity,
        "caps": {name.decode("ascii"): status_fields[name].decode("ascii")
                 for name in (b"CapInh", b"CapPrm", b"CapEff", b"CapBnd", b"CapAmb")},
        "cmdline_sha256": sha256_bytes(cmdline_data),
        "environment_sha256": sha256_bytes(environment),
        "executable_stat": executable_receipt,
        "gid": 1000,
        "groups": [],
        "no_new_privs": 1,
        "pid": pid,
        "process_start_ticks": process_start_ticks(pid),
        "uid": 1000,
    }


def validate_controller_boundary(args: argparse.Namespace) -> Dict[str, Any]:
    if (
        not hasattr(os, "getresuid") or os.getresuid() != (1000, 1000, 1000) or
        not hasattr(os, "getresgid") or os.getresgid() != (1000, 1000, 1000) or
        os.getgroups() != []
    ):
        fail("controller UID/GID/supplementary-group boundary differs")
    if dict(os.environ) != SEALED_ENV:
        fail("controller environment is not the exact sealed four-key map")
    if (
        not sys.flags.isolated or not sys.dont_write_bytecode or
        sys.flags.optimize != 0
    ):
        fail("controller requires /usr/bin/python3 -I -B without optimization")
    controller_path = str(Path(__file__).resolve(strict=True))
    expected_cmdline = [
        "/usr/bin/python3", "-I", "-B", controller_path,
    ] + list(sys.argv[1:])
    python_path = exact_absolute_path(
        args.expected_python_image_path, "expected Python image path")
    if os.path.realpath("/usr/bin/python3") != python_path:
        fail("/usr/bin/python3 does not resolve to the externally sealed image")
    lower_hash(args.expected_python_image_sha256,
               "expected Python image SHA-256")
    python_info = os.stat(python_path, follow_symlinks=False)
    if not stat.S_ISREG(python_info.st_mode):
        fail("expected Python image is not regular")
    return process_security_receipt(
        os.getpid(), (120, 121, 122, 123), stat_receipt(python_info),
        expected_cmdline)


class OuterGuardian:
    def __init__(self, deadline_monotonic_ns: int) -> None:
        exact_int(deadline_monotonic_ns, time.monotonic_ns() + 1, MAX_INT63,
                  "outer guardian absolute deadline")
        read_fd, write_fd = os.pipe2(getattr(os, "O_CLOEXEC", 0))
        ready_read, ready_write = os.pipe2(getattr(os, "O_CLOEXEC", 0))
        self.started_monotonic_ns = time.monotonic_ns()
        argv = [
            "/usr/bin/python3", "-I", "-S", "-c", GUARDIAN_CODE,
            str(read_fd), str(ready_write), str(os.getpid()),
            str(deadline_monotonic_ns), "123",
        ]
        self.argv = list(argv)
        try:
            self.process = subprocess.Popen(
                argv, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL, close_fds=True,
                pass_fds=(read_fd, ready_write), start_new_session=True,
                env=dict(SEALED_ENV),
            )
        except BaseException:
            os.close(read_fd)
            os.close(write_fd)
            os.close(ready_read)
            os.close(ready_write)
            raise
        os.close(read_fd)
        os.close(ready_write)
        self.write_fd = write_fd
        self.completed = False
        try:
            ready, _, _ = __import__("select").select([ready_read], [], [], 5.0)
            token = os.read(ready_read, 256) if ready else b""
        finally:
            os.close(ready_read)
        parts = token.split()
        if (
            len(parts) != 2 or any(not part.isdigit() for part in parts) or
            int(parts[1]) != deadline_monotonic_ns or
            self.process.poll() is not None
        ):
            self.close_descriptor()
            self.process.kill()
            self.process.wait(timeout=5.0)
            fail("outer guardian did not reach its pinned ready state")
        affinity = sorted(os.sched_getaffinity(self.process.pid))
        if affinity != [123]:
            self.close_descriptor()
            self.process.kill()
            self.process.wait(timeout=5.0)
            fail("outer guardian affinity differs from CPU123")
        self.affinity = affinity
        self.process_start_ticks = process_start_ticks(self.process.pid)
        self.armed_monotonic_ns = int(parts[0])
        self.deadline_monotonic_ns = deadline_monotonic_ns
        if self.armed_monotonic_ns >= self.deadline_monotonic_ns:
            self.close_descriptor()
            self.process.kill()
            self.process.wait(timeout=5.0)
            fail("outer guardian armed at or after its deadline")

    def poll(self) -> str:
        """Return a bounded terminal/identity event without raising."""
        code = self.process.poll()
        if code is not None:
            return "guardian-exited:{}".format(code)
        try:
            if process_start_ticks(self.process.pid) != self.process_start_ticks:
                return "guardian-start-identity-changed"
            if sorted(os.sched_getaffinity(self.process.pid)) != self.affinity:
                return "guardian-affinity-changed"
            if os.getpgid(self.process.pid) != self.process.pid:
                return "guardian-process-group-changed"
        except BaseException as exc:
            return "guardian-monitor-invalid:" + exception_text(exc)
        return "none"

    def receipt(self) -> Dict[str, Any]:
        return {
            "affinity": list(self.affinity),
            "code_sha256": sha256_bytes(GUARDIAN_CODE.encode("ascii")),
            "armed_monotonic_ns": self.armed_monotonic_ns,
            "deadline_monotonic_ns": self.deadline_monotonic_ns,
            "deadline_seconds": EXTERNAL_DEADLINE_SECONDS,
            "pid": self.process.pid,
            "process_start_ticks": self.process_start_ticks,
            "started_monotonic_ns": self.started_monotonic_ns,
            "worker_stub_sha256": sha256_bytes(WORKER_STUB_CODE.encode("ascii")),
        }

    def _write(self, data: bytes) -> None:
        if self.completed or self.write_fd < 0:
            fail("outer guardian control pipe is closed")
        offset = 0
        while offset < len(data):
            written = os.write(self.write_fd, data[offset:])
            if written <= 0:
                fail("outer guardian control write made no progress")
            offset += written

    def register(self, worker: WorkerProcess) -> None:
        event = self.poll()
        if event != "none":
            fail("cannot register worker with unhealthy guardian: " + event)
        self._write("pid {}\n".format(worker.process.pid).encode("ascii"))
        event = self.poll()
        if event != "none":
            fail("guardian failed while registering worker: " + event)

    def complete(self) -> None:
        if self.completed:
            return
        event = self.poll()
        if event != "none":
            fail("outer guardian was not live at completion: " + event)
        self._write(b"done\n")
        os.close(self.write_fd)
        self.write_fd = -1
        try:
            code = self.process.wait(timeout=5.0)
        except subprocess.TimeoutExpired:
            self.process.kill()
            self.process.wait(timeout=5.0)
            fail("outer guardian did not acknowledge completion")
        self.completed = True
        if code != 0:
            fail("outer guardian exited {}".format(code))

    def close_descriptor(self) -> None:
        if self.write_fd >= 0:
            os.close(self.write_fd)
            self.write_fd = -1


def wait_process_security(
    process: subprocess.Popen, affinity: Sequence[int],
    executable_stat: Mapping[str, Any], cmdline: Sequence[str],
    deadline: float, authority_poll: Any, where: str,
) -> Dict[str, Any]:
    local_deadline = min(deadline, time.monotonic() + 2.0)
    last_error = "process did not enter its sealed state"
    while time.monotonic() < local_deadline:
        event = authority_poll()
        if event != "none":
            fail("{} lost supervision: {}".format(where, event))
        if process.poll() is not None:
            fail("{} exited {} before security validation".format(
                where, process.returncode))
        try:
            return process_security_receipt(
                process.pid, affinity, executable_stat, cmdline)
        except ValidationError as exc:
            last_error = str(exc)
        time.sleep(0.005)
    fail("{} did not reach its sealed security state: {}".format(
        where, last_error))


@dataclass
class PreparedInputs:
    artifacts: Dict[str, OpenArtifact]
    build_receipts: Dict[str, Dict[str, Any]]
    closure_sha256: Dict[str, str]

    def close(self) -> None:
        for artifact in self.artifacts.values():
            artifact.close()


def prepare_inputs(args: argparse.Namespace) -> PreparedInputs:
    exact_int(args.cpu, 120, 120, "scientific target CPU")
    exact_int(args.controller_cpu, 121, 121, "scientific controller CPU")
    exact_int(args.sampler_cpu, 122, 122, "scientific sampler CPU")
    exact_string(args.expected_current_implementation_commit,
                 CURRENT_IMPLEMENTATION_COMMIT,
                 "expected current implementation commit")
    exact_string(args.expected_parent_implementation_commit, PARENT_COMMIT,
                 "expected parent implementation commit")
    lower_commit(args.expected_harness_commit, "expected harness commit")
    lower_hash(args.expected_controller_sha256,
               "expected controller SHA-256")
    exact_absolute_path(args.output_dir, "output directory")
    expected_output_parent(args)

    path_specs = (
        ("current_worker", args.current_worker,
         args.expected_current_worker_sha256, 256 * 1024 * 1024, True,
         0, 0, 0o555),
        ("parent_worker", args.parent_worker,
         args.expected_parent_worker_sha256, 256 * 1024 * 1024, True,
         0, 0, 0o555),
        ("current_library", args.current_library,
         args.expected_current_library_sha256, 1024 * 1024 * 1024, False,
         0, 0, 0o444),
        ("parent_library", args.parent_library,
         args.expected_parent_library_sha256, 1024 * 1024 * 1024, False,
         0, 0, 0o444),
        ("current_implementation_library", args.current_implementation_library,
         args.expected_current_implementation_library_sha256,
         1024 * 1024 * 1024, False, 0, 0, 0o444),
        ("parent_implementation_library", args.parent_implementation_library,
         args.expected_parent_implementation_library_sha256,
         1024 * 1024 * 1024, False, 0, 0, 0o444),
        ("current_build_receipt", args.current_build_receipt,
         args.expected_current_build_receipt_sha256, MAX_RECEIPT_BYTES, False,
         0, 0, 0o444),
        ("parent_build_receipt", args.parent_build_receipt,
         args.expected_parent_build_receipt_sha256, MAX_RECEIPT_BYTES, False,
         0, 0, 0o444),
        ("python_image", args.expected_python_image_path,
         args.expected_python_image_sha256, 256 * 1024 * 1024, True,
         0, 0, 0o755),
    )
    artifacts: Dict[str, OpenArtifact] = {}
    try:
        for (name, path, digest, maximum, executable, uid, gid,
             mode) in path_specs:
            artifacts[name] = open_artifact(
                name, path, digest, maximum, executable=executable,
                expected_uid=uid, expected_gid=gid, expected_mode=mode,
                expected_nlink=(
                    EXPECTED_PYTHON_IMAGE_NLINK
                    if name == "python_image" else 1
                ))
        controller_path = str(Path(__file__).resolve(strict=True))
        artifacts["controller"] = open_artifact(
            "controller", controller_path, args.expected_controller_sha256,
            8 * 1024 * 1024, expected_uid=1000, expected_gid=1000,
            expected_mode=0o444)
        if stat_receipt(os.stat("/proc/self/exe")) != stat_receipt(
            artifacts["python_image"].before
        ):
            fail("running Python image differs from the held external image")

        if Path(args.current_worker).name != (
            "wirehair_wh2_v2_facade_timing_worker_current"
        ) or Path(args.parent_worker).name != (
            "wirehair_wh2_v2_facade_timing_worker_parent"
        ):
            fail("worker basenames differ from their sealed roles")
        for role in ("current", "parent"):
            worker_path = getattr(args, role + "_worker")
            library_path = getattr(args, role + "_library")
            if (
                Path(library_path).name != "libwirehair.so.2" or
                Path(worker_path).parent != Path(library_path).parent
            ):
                fail("{} runtime library is not the adjacent SONAME image".format(role))
            if getattr(args, "expected_" + role + "_library_sha256") != getattr(
                args, "expected_" + role + "_implementation_library_sha256"
            ):
                fail("{} runtime/implementation library hashes differ".format(role))
        if (
            artifacts["current_worker"].before.st_ino ==
                artifacts["parent_worker"].before.st_ino and
            artifacts["current_worker"].before.st_dev ==
                artifacts["parent_worker"].before.st_dev
        ):
            fail("role workers alias the same file")

        build_receipts: Dict[str, Dict[str, Any]] = {}
        controller_root = str(Path(controller_path).parents[1])
        for role, implementation_commit in (
            ("current", CURRENT_IMPLEMENTATION_COMMIT),
            ("parent", PARENT_COMMIT),
        ):
            receipt_artifact = artifacts[role + "_build_receipt"]
            receipt = parse_canonical_document(
                read_artifact(receipt_artifact, MAX_RECEIPT_BYTES),
                role + " build receipt", MAX_RECEIPT_BYTES)
            validate_build_receipt(
                receipt, role, args.expected_harness_commit,
                implementation_commit, getattr(args, role + "_worker"),
                getattr(args, "expected_" + role + "_worker_sha256"),
                getattr(args, role + "_library"),
                getattr(args, "expected_" + role + "_library_sha256"),
                getattr(args, role + "_implementation_library"),
                getattr(args,
                        "expected_" + role + "_implementation_library_sha256"),
                controller_root, controller_path,
                args.expected_controller_sha256,
            )
            build_receipts[role] = receipt
        for field_name in (
            "git_path", "git_sha256", "harness_source_root",
            "harness_tree_oid", "harness_tree_listing_sha256",
        ):
            if build_receipts["current"][field_name] != build_receipts["parent"][
                field_name
            ]:
                fail("role build receipts disagree on " + field_name)
        current_harness = [
            item for item in build_receipts["current"]["source_manifest"]
            if item["authority"] == "harness"
        ]
        parent_harness = [
            item for item in build_receipts["parent"]["source_manifest"]
            if item["authority"] == "harness"
        ]
        if current_harness != parent_harness:
            fail("role build receipts disagree on the harness source roster")
        artifacts["git"] = open_artifact(
            "Git executable", build_receipts["current"]["git_path"],
            build_receipts["current"]["git_sha256"],
            256 * 1024 * 1024, executable=True, expected_uid=0,
            expected_gid=0, expected_mode=0o755)
        closure_sha256: Dict[str, str] = {}
        for role in ("current", "parent"):
            closure_sha256[role] = verify_build_receipt_closure(
                build_receipts[role], role, artifacts["git"])
        return PreparedInputs(
            artifacts=artifacts, build_receipts=build_receipts,
            closure_sha256=closure_sha256)
    except BaseException:
        for artifact in artifacts.values():
            artifact.close()
        raise


@dataclass
class ScientificResult:
    raw_complete: bool = False
    statistics: Optional[Dict[str, Any]] = None
    gate_failures: List[str] = field(default_factory=list)
    infrastructure_failures: List[str] = field(default_factory=list)
    worker_terminals: Optional[Dict[str, Any]] = None
    health: Optional[Dict[str, Any]] = None
    health_state: Optional[Dict[str, Any]] = None
    child_start_ns: Optional[int] = None
    child_reap_ns: Optional[int] = None
    guardian: Optional[Dict[str, Any]] = None
    worker_processes: Optional[Dict[str, Any]] = None


def execute_scientific_run(
    args: argparse.Namespace, bundle: OutputBundle,
    health_adapter: HealthAdapter, deadline: float,
    external_deadline_monotonic_ns: int,
    artifacts: Mapping[str, OpenArtifact],
) -> ScientificResult:
    result = ScientificResult()
    workers_list: List[WorkerProcess] = []
    mux: Optional[WorkerMux] = None
    guardian: Optional[OuterGuardian] = None
    health_state: Optional[Dict[str, Any]] = None
    terminals: Dict[str, Any] = {}
    worker_deadline = deadline - HEALTH_FINALIZATION_RESERVE_SECONDS
    try:
        initial_affinity = sorted(os.sched_getaffinity(0))
        args.controller_initial_affinity = initial_affinity
        if initial_affinity != [120, 121, 122, 123]:
            fail("controller initial affinity is not exactly CPUs120-123")
        os.sched_setaffinity(0, {args.controller_cpu})
        if os.sched_getaffinity(0) != {args.controller_cpu}:
            fail("controller did not retain singleton CPU121 affinity")

        guardian = OuterGuardian(external_deadline_monotonic_ns)
        result.guardian = guardian.receipt()
        result.guardian["security"] = wait_process_security(
            guardian.process, (123,), stat_receipt(artifacts["python_image"].before),
            guardian.argv, deadline, guardian.poll, "outer guardian")
        health_state = health_adapter.module.prepare_health(
            args, args.cpu, deadline, health_adapter.modules)
        result.health_state = health_state
        if not health_state.get("ready"):
            fail("hardened health admission did not license workers")
        sampler_handles = health_state.get("sampler_handles")
        if not sampler_handles:
            fail("health admission lacks live sampler handles")
        prelaunch_event = health_adapter.module.poll_sampler_supervision(
            sampler_handles)
        if prelaunch_event != "none":
            health_state["sampler_supervision_event"] = prelaunch_event
            fail("sampler terminal state preceded workers: " + prelaunch_event)

        def prelaunch_authority_poll() -> str:
            guardian_event = guardian.poll()
            if guardian_event != "none":
                return guardian_event
            try:
                return health_adapter.module.poll_sampler_supervision(
                    sampler_handles)
            except BaseException as exc:
                return "sampler-monitor-invalid:" + exception_text(exc)

        result.child_start_ns = time.monotonic_ns()
        blocked_receipts: Dict[str, Any] = {}
        current_worker = launch_worker("current", args.current_worker, args.cpu)
        workers_list.append(current_worker)
        guardian.register(current_worker)
        blocked_receipts["current"] = wait_process_security(
            current_worker.process, (args.cpu,),
            stat_receipt(artifacts["python_image"].before),
            current_worker.stub_argv, worker_deadline,
            prelaunch_authority_poll, "blocked current worker stub")
        parent_worker = launch_worker("parent", args.parent_worker, args.cpu)
        workers_list.append(parent_worker)
        guardian.register(parent_worker)
        blocked_receipts["parent"] = wait_process_security(
            parent_worker.process, (args.cpu,),
            stat_receipt(artifacts["python_image"].before),
            parent_worker.stub_argv, worker_deadline,
            prelaunch_authority_poll, "blocked parent worker stub")
        guardian_event = guardian.poll()
        if guardian_event != "none":
            fail("outer guardian failed before worker release: " + guardian_event)
        current_worker.release()
        parent_worker.release()
        workers = {"current": current_worker, "parent": parent_worker}
        mux = WorkerMux(
            workers_list,
            sampler_poll=lambda: health_adapter.module.poll_sampler_supervision(
                sampler_handles),
            guardian_poll=guardian.poll,
        )
        current_line = mux.read_line(
            current_worker, worker_deadline, "reading current startup config")
        parent_line = mux.read_line(
            parent_worker, worker_deadline, "reading parent startup config")
        mux.require_quiet()
        current_config = parse_canonical_line(current_line, "current startup config")
        parent_config = parse_canonical_line(parent_line, "parent startup config")
        current_contract = validate_worker_config(
            current_config, "current", args.cpu, args.expected_harness_commit,
            args.expected_current_implementation_commit)
        parent_contract = validate_worker_config(
            parent_config, "parent", args.cpu, args.expected_harness_commit,
            args.expected_parent_implementation_commit)
        validate_cross_role_contracts(current_contract, parent_contract)
        current_worker.contract = current_contract
        parent_worker.contract = parent_contract
        validate_worker_affinity(current_worker)
        validate_worker_affinity(parent_worker)
        process_receipts: Dict[str, Any] = {}
        for worker in workers_list:
            try:
                pgid = os.getpgid(worker.process.pid)
                executable_info = os.stat(
                    "/proc/{}/exe".format(worker.process.pid))
            except OSError as exc:
                fail("{} worker process receipt failed: {}".format(
                    worker.role, exc))
            if pgid != worker.process.pid:
                fail("{} worker lacks its isolated process group".format(worker.role))
            expected_worker = artifacts[worker.role + "_worker"]
            released = process_security_receipt(
                worker.process.pid, (args.cpu,), stat_receipt(expected_worker.before),
                [worker.path, "--serve", "--role", worker.role,
                 "--cpu", str(args.cpu)])
            if (
                released["process_start_ticks"] !=
                blocked_receipts[worker.role]["process_start_ticks"]
            ):
                fail("{} worker PID identity changed across exec".format(worker.role))
            process_receipts[worker.role] = {
                "artifact": expected_worker.receipt(),
                "blocked_stub": blocked_receipts[worker.role],
                "pgid": pgid,
                "released_worker": released,
                "terminal_boundary": None,
            }
        result.worker_processes = process_receipts
        if result.guardian is not None:
            result.guardian["registered_worker_pgids"] = [
                process_receipts[role]["pgid"] for role in ("current", "parent")
            ]
        screen_config = make_screen_config(current_contract, parent_contract)
        bundle.append_raw(screen_config)

        seq = 0
        panel_count = 0
        for replicate in range(REPLICATES):
            for k, block_bytes in CELLS:
                for scope in SCOPES:
                    for comparison, _, _ in COMPARISONS:
                        panel, seq = run_panel(
                            mux, workers, seq, replicate, k, block_bytes,
                            scope, comparison, worker_deadline)
                        bundle.append_raw(panel)
                        panel_count += 1
        if panel_count != EXPECTED_PANEL_COUNT:
            fail("scientific panel roster is incomplete")

        mux.check_sampler(force=True)
        for worker in workers_list:
            expected_worker = artifacts[worker.role + "_worker"]
            terminal_boundary = process_security_receipt(
                worker.process.pid, (args.cpu,), stat_receipt(expected_worker.before),
                [worker.path, "--serve", "--role", worker.role,
                 "--cpu", str(args.cpu)])
            if terminal_boundary["process_start_ticks"] != process_receipts[
                worker.role
            ]["released_worker"]["process_start_ticks"]:
                fail("{} worker PID identity changed before terminal".format(
                    worker.role))
            process_receipts[worker.role]["terminal_boundary"] = terminal_boundary
        current_worker.send(b"quit\n")
        parent_worker.send(b"quit\n")
        for worker in workers_list:
            line = mux.read_line(
                worker, worker_deadline,
                "reading {} terminal".format(worker.role))
            value = parse_canonical_line(
                line, "{} terminal".format(worker.role))
            if worker.contract is None:
                fail("{} terminal lacks contract".format(worker.role))
            validate_worker_terminal(value, worker.contract, worker.invocation_count)
            worker.terminal_received = True
            terminals[worker.role] = value
        mux.require_quiet()
        mux.drain_and_reap(worker_deadline)
        result.child_reap_ns = time.monotonic_ns()
        raw_terminal = {
            "current_invocations": current_worker.invocation_count,
            "panel_count": panel_count,
            "parent_invocations": parent_worker.invocation_count,
            "raw_config_sha256": sha256_bytes(canonical_bytes(screen_config)),
            "schema": RAW_TERMINAL_SCHEMA,
            "status": "complete",
        }
        bundle.append_raw(raw_terminal)
        raw = bundle.raw_bytes()
        _config, _terminal, statistics, failures = parse_raw_transcript(
            raw, args.cpu, args.expected_harness_commit,
            args.expected_current_implementation_commit)
        result.statistics = statistics
        result.gate_failures = failures
        result.worker_terminals = terminals
        result.raw_complete = True
    except BaseException as exc:
        result.infrastructure_failures.append(exception_text(exc))
    finally:
        if mux is not None and health_state is not None:
            health_state["sampler_supervision_event"] = mux.sampler_event
        if mux is not None:
            result.infrastructure_failures.extend(mux.contain_and_drain())
        else:
            for worker in workers_list:
                worker.terminate()
            for worker in workers_list:
                try:
                    _stdout, stderr = worker.process.communicate(timeout=5.0)
                    if stderr:
                        available = MAX_WORKER_STDERR_BYTES - len(worker.stderr)
                        worker.stderr.extend(stderr[:available])
                        if len(stderr) > available:
                            worker.stderr_overflow = True
                except BaseException as exc:
                    result.infrastructure_failures.append(
                        "{} containment: {}".format(
                            worker.role, exception_text(exc)))
        if workers_list and result.child_reap_ns is None:
            result.child_reap_ns = time.monotonic_ns()
        if mux is not None:
            try:
                mux.close()
            except BaseException as exc:
                result.infrastructure_failures.append(
                    "worker mux close: " + exception_text(exc))
        if health_state is not None:
            try:
                health, health_errors = health_adapter.module.finish_health(
                    args, args.cpu, result.child_start_ns, result.child_reap_ns,
                    health_state, deadline, health_adapter.modules)
                result.health = health
                result.infrastructure_failures.extend(health_errors)
            except BaseException as exc:
                result.infrastructure_failures.append(
                    "health finish: " + exception_text(exc))
                receipt = health_state.get("receipt")
                if type(receipt) is dict:
                    result.health = dict(receipt)
        if guardian is not None:
            try:
                guardian_event = guardian.poll()
                if guardian_event != "none":
                    fail("outer guardian failed before health closure: " +
                         guardian_event)
                guardian.complete()
            except BaseException as exc:
                result.infrastructure_failures.append(
                    "outer guardian: " + exception_text(exc))
                guardian.close_descriptor()
        for worker in workers_list:
            stderr_name = worker.role + ".stderr"
            try:
                bundle.files[stderr_name].append(
                    bytes(worker.stderr), MAX_WORKER_STDERR_BYTES)
            except BaseException as exc:
                result.infrastructure_failures.append(
                    "{} publication: {}".format(stderr_name, exception_text(exc)))
            if worker.stderr:
                result.infrastructure_failures.append(
                    "{} worker emitted stderr".format(worker.role))
            if worker.stderr_overflow:
                result.infrastructure_failures.append(
                    "{} worker stderr overflowed".format(worker.role))
    return result


def bounded_unique_failures(values: Iterable[str]) -> List[str]:
    result: List[str] = []
    for value in values:
        if type(value) is not str or not value:
            value = "invalid empty failure diagnostic"
        if len(value) > 4096:
            value = value[:4000] + "...[sha256:{}]".format(
                sha256_bytes(value.encode("utf-8")))
        if value not in result:
            result.append(value)
        if len(result) >= MAX_FAILURES:
            break
    return result


def run_once(args: argparse.Namespace) -> int:
    prepared: Optional[PreparedInputs] = None
    bundle: Optional[OutputBundle] = None
    try:
        args.controller_security_boundary = validate_controller_boundary(args)
        prepared = prepare_inputs(args)
        source_root = Path(__file__).resolve(strict=True).parents[1]
        health_adapter = load_health_adapter(
            source_root, args.expected_health_source_manifest_sha256,
            args.expected_health_adapter_sha256,
            args.expected_harness_commit,
            prepared.build_receipts["current"]["git_sha256"])
        health_adapter.module.validate_sampler_source_authority(
            args, health_adapter.source_manifest, source_root)
        bundle = OutputBundle(args.output_dir, expected_output_parent(args))
    except BaseException as exc:
        if prepared is not None:
            prepared.close()
        print("preflight invalid: {}".format(exception_text(exc)), file=sys.stderr)
        return 1

    start_ns = time.monotonic_ns()
    deadline_ns = start_ns + INTERNAL_DEADLINE_SECONDS * 1_000_000_000
    external_deadline_ns = start_ns + EXTERNAL_DEADLINE_SECONDS * 1_000_000_000
    deadline = deadline_ns / 1_000_000_000.0
    result = execute_scientific_run(
        args, bundle, health_adapter, deadline, external_deadline_ns,
        prepared.artifacts)
    infrastructure = list(result.infrastructure_failures)
    try:
        if result.worker_processes is not None:
            for role in ("current", "parent"):
                if result.worker_processes[role]["released_worker"][
                    "executable_stat"
                ] != stat_receipt(prepared.artifacts[role + "_worker"].before):
                    fail("{} running executable differs from held worker".format(role))
        for artifact in prepared.artifacts.values():
            artifact.verify_unchanged()
        for role in ("current", "parent"):
            after = verify_build_receipt_closure(
                prepared.build_receipts[role], role, prepared.artifacts["git"])
            if after != prepared.closure_sha256[role]:
                fail("{} build closure changed during execution".format(role))
        adapter_data = (
            Path(__file__).resolve(strict=True).parents[1] / HEALTH_ADAPTER_PATH
        ).read_bytes()
        if sha256_bytes(adapter_data) != args.expected_health_adapter_sha256:
            fail("health adapter changed during execution")
        source_root = Path(__file__).resolve(strict=True).parents[1]
        health_manifest_after = health_adapter.module.source_manifest(source_root)
        if canonical_bytes(health_manifest_after) != canonical_bytes(
            health_adapter.source_manifest
        ):
            fail("health source manifest changed during execution")
        health_git_after = health_adapter.module.git_receipt(
            source_root, args.expected_harness_commit,
            prepared.build_receipts["current"]["git_sha256"])
        if canonical_bytes(health_git_after) != canonical_bytes(
            health_adapter.git_receipt
        ):
            fail("health source Git receipt changed during execution")
    except BaseException as exc:
        infrastructure.append("postflight: " + exception_text(exc))
    end_ns = time.monotonic_ns()
    if end_ns >= start_ns + INTERNAL_DEADLINE_SECONDS * 1_000_000_000:
        infrastructure.append("controller reached the inclusive 840-second deadline")
    infrastructure = bounded_unique_failures(infrastructure)
    outcome = classify_outcome(
        infrastructure, result.raw_complete, result.gate_failures,
    )
    error = "; ".join(infrastructure) if infrastructure else "none"

    publication_failed = False
    try:
        provenance = make_provenance(
            args, prepared.artifacts, prepared.build_receipts, start_ns, end_ns,
            error, prepared.closure_sha256, health_adapter, result.health,
            result.guardian, result.worker_processes, bundle.receipt())
        bundle.write_once("provenance.json", provenance, MAX_RECEIPT_BYTES)
        provenance_sha = bundle.files["provenance.json"].digest.hexdigest()
        aa_invalid = any(
            item.startswith("aa_ci:") for item in result.gate_failures
        )
        summary_failures = (
            bounded_unique_failures([
                *infrastructure,
                *(result.gate_failures if aa_invalid else ()),
            ])
            if outcome == "invalid" else result.gate_failures
        )
        summary = make_summary(
            outcome, bundle.files["raw.jsonl"], result.raw_complete,
            result.statistics if result.raw_complete else None,
            summary_failures, provenance_sha, result.worker_terminals,
            result.health)
        bundle.write_once("summary.json", summary, MAX_RECEIPT_BYTES)
        bundle.seal(outcome)
    except BaseException as exc:
        publication_failed = True
        print("publication failed: {}".format(exception_text(exc)), file=sys.stderr)
        try:
            bundle.emergency_invalid("publication failed: " + exception_text(exc))
        except BaseException as terminal_exc:
            print("emergency terminalization failed: {}".format(
                exception_text(terminal_exc)), file=sys.stderr)
    finally:
        bundle.close()
        prepared.close()
    if publication_failed or outcome == "invalid":
        return 1
    return 0 if outcome == "pass" else 2


REPLAY_SCHEMA = "wirehair.wh2.v2-facade-timing-screen.replay.v1"


def replay(args: argparse.Namespace) -> int:
    exact_int(args.cpu, 120, 120, "replay target CPU")
    lower_commit(args.expected_harness_commit, "replay harness commit")
    exact_string(args.expected_current_implementation_commit,
                 CURRENT_IMPLEMENTATION_COMMIT,
                 "replay current implementation commit")
    path = exact_absolute_path(args.replay, "replay raw path")
    data = Path(path).read_bytes()
    _config, terminal, statistics, failures = parse_raw_transcript(
        data, args.cpu, args.expected_harness_commit,
        args.expected_current_implementation_commit)
    outcome = classify_outcome([], True, failures)
    value = {
        "campaign": CAMPAIGN,
        "failed_gates": failures,
        "outcome": outcome,
        "raw_bytes": len(data),
        "raw_sha256": sha256_bytes(data),
        "raw_terminal": terminal,
        "schema": REPLAY_SCHEMA,
        "statistics": statistics,
    }
    sys.stdout.buffer.write(canonical_bytes(value) + b"\n")
    return {"pass": 0, "invalid": 1, "reject": 2}[outcome]


def canonical_uint(text: str) -> int:
    if re.fullmatch(r"0|[1-9][0-9]*", text) is None:
        raise argparse.ArgumentTypeError("expected a canonical unsigned integer")
    value = int(text)
    if value > MAX_UINT64:
        raise argparse.ArgumentTypeError("unsigned integer exceeds uint64")
    return value


def add_run_arguments(parser: argparse.ArgumentParser) -> None:
    for option in (
        "current-worker", "parent-worker", "current-library", "parent-library",
        "current-implementation-library", "parent-implementation-library",
        "current-build-receipt", "parent-build-receipt", "output-dir",
        "expected-current-worker-sha256", "expected-parent-worker-sha256",
        "expected-current-library-sha256", "expected-parent-library-sha256",
        "expected-current-implementation-library-sha256",
        "expected-parent-implementation-library-sha256",
        "expected-current-build-receipt-sha256",
        "expected-parent-build-receipt-sha256", "expected-harness-commit",
        "expected-controller-sha256",
        "expected-python-image-path", "expected-python-image-sha256",
        "expected-current-implementation-commit",
        "expected-parent-implementation-commit",
        "expected-health-source-manifest-sha256",
        "expected-health-adapter-sha256", "sampler-script", "sampler-csv",
        "sampler-pid-file", "sampler-validation-jsonl", "sampler-receipt",
        "expected-sampler-script-sha256", "expected-sampler-cmdline-sha256",
        "expected-sampler-environ-sha256", "expected-sampler-executable-sha256",
    ):
        parser.add_argument("--" + option)
    for option in (
        "cpu", "controller-cpu", "sampler-cpu", "sampler-pid",
        "expected-sampler-process-start-ticks", "expected-sampler-csv-device",
        "expected-sampler-csv-inode", "expected-sampler-pid-file-device",
        "expected-sampler-pid-file-inode", "expected-sampler-validation-device",
        "expected-sampler-validation-inode", "expected-sampler-receipt-device",
        "expected-sampler-receipt-inode", "expected-sampler-uid",
        "expected-sampler-gid", "expected-sampler-i2c-gid",
        "expected-output-parent-device", "expected-output-parent-inode",
        "expected-output-parent-uid", "expected-output-parent-gid",
        "expected-output-parent-mode",
    ):
        parser.add_argument("--" + option, type=canonical_uint)


RUN_REQUIRED = (
    "current_worker", "parent_worker", "current_library", "parent_library",
    "current_implementation_library", "parent_implementation_library",
    "current_build_receipt", "parent_build_receipt", "output_dir",
    "expected_current_worker_sha256", "expected_parent_worker_sha256",
    "expected_current_library_sha256", "expected_parent_library_sha256",
    "expected_current_implementation_library_sha256",
    "expected_parent_implementation_library_sha256",
    "expected_current_build_receipt_sha256",
    "expected_parent_build_receipt_sha256", "expected_harness_commit",
    "expected_controller_sha256",
    "expected_python_image_path", "expected_python_image_sha256",
    "expected_current_implementation_commit",
    "expected_parent_implementation_commit",
    "expected_health_source_manifest_sha256", "expected_health_adapter_sha256",
    "cpu", "controller_cpu", "sampler_cpu", "sampler_pid", "sampler_script",
    "sampler_csv", "sampler_pid_file", "sampler_validation_jsonl",
    "sampler_receipt", "expected_sampler_process_start_ticks",
    "expected_sampler_script_sha256", "expected_sampler_csv_device",
    "expected_sampler_csv_inode", "expected_sampler_pid_file_device",
    "expected_sampler_pid_file_inode", "expected_sampler_validation_device",
    "expected_sampler_validation_inode", "expected_sampler_receipt_device",
    "expected_sampler_receipt_inode", "expected_sampler_cmdline_sha256",
    "expected_sampler_environ_sha256", "expected_sampler_executable_sha256",
    "expected_sampler_uid", "expected_sampler_gid", "expected_sampler_i2c_gid",
    "expected_output_parent_device", "expected_output_parent_inode",
    "expected_output_parent_uid", "expected_output_parent_gid",
    "expected_output_parent_mode",
)


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Strict public V2 facade timing falsifier")
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--selftest", action="store_true")
    modes.add_argument("--run-once", action="store_true")
    modes.add_argument("--replay")
    add_run_arguments(parser)
    args = parser.parse_args(argv)
    if args.run_once:
        missing = [name for name in RUN_REQUIRED if getattr(args, name) is None]
        if missing:
            parser.error("--run-once lacks: " + ", ".join(missing))
        for name in (
            "current_worker", "parent_worker", "current_library", "parent_library",
            "current_implementation_library", "parent_implementation_library",
            "current_build_receipt", "parent_build_receipt", "output_dir",
            "expected_python_image_path",
            "sampler_script", "sampler_csv", "sampler_pid_file",
            "sampler_validation_jsonl", "sampler_receipt",
        ):
            exact_absolute_path(getattr(args, name), "--" + name.replace("_", "-"))
        args.sampler_script = Path(args.sampler_script)
        args.sampler_csv = Path(args.sampler_csv)
        args.sampler_pid_file = Path(args.sampler_pid_file)
        args.sampler_validation_jsonl = Path(args.sampler_validation_jsonl)
        args.sampler_receipt = Path(args.sampler_receipt)
    elif args.replay is not None:
        required = (
            "cpu", "expected_harness_commit",
            "expected_current_implementation_commit",
        )
        missing = [name for name in required if getattr(args, name) is None]
        if missing:
            parser.error("--replay lacks: " + ", ".join(missing))
        supplied = [
            name for name in RUN_REQUIRED
            if name not in required and getattr(args, name) is not None
        ]
        if supplied:
            parser.error("--replay received run-only options: " + ", ".join(supplied))
    else:
        supplied = [name for name in RUN_REQUIRED if getattr(args, name) is not None]
        if supplied:
            parser.error("--selftest received run-only options: " + ", ".join(supplied))
    return args


def main(argv: Sequence[str]) -> int:
    try:
        args = parse_args(argv)
        if args.selftest:
            return selftest()
        if args.replay is not None:
            return replay(args)
        return run_once(args)
    except ValidationError as exc:
        print("invalid: {}".format(exc), file=sys.stderr)
        return 1
    except (OSError, ValueError) as exc:
        print("invalid: {}".format(exception_text(exc)), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
