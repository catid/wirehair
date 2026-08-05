#!/usr/bin/env python3
"""Run the bounded D12-reference versus Two07 WH2 recovery screen.

``run`` freezes exactly four same-binary arms: descriptive production WH2 and
Wirehair1 controls, the uniform-raw D12/H12 recovery reference, and the
uniform-raw Two07 candidate.  It binds the D12-to-production timing-proxy
witness and the two-arm work/rank sidecar before emitting recovery results.
The proxy witness is structure-only: at development attempt 0 it proves that
the disabled-anchor D12 structure under the production timing seed policy is
the production head, while separately binding the distinct uniform-raw
recovery seeds.  It does not equate their graph, peel, rank, or solve work.
The compatibility ``combine`` command may re-open that sole campaign as a
logical bundle, but is not itself a native execution receipt.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import os
from pathlib import Path
import re
import stat
import subprocess
import sys
import tempfile
import time
from types import SimpleNamespace
from typing import (Any, Callable, Dict, List, Mapping, Optional, Sequence, Set,
                    Tuple)

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_benchmark_contract as contract_api
import wh2_native_short_screen as native_api
import wh2_precodefail_work_screen as work_api
import wh2_run_native_short_screen as runner_api


CANDIDATE_SPECS = (
    ("two07", "wirehair2_dense_two07_basis_v1",
     "9527f200ad38c7eec6502b2f768fdd67b92787fb227eed3d7616274ffc2df388"),
)
CANDIDATE_BY_ID = {
    candidate_id: (arm, descriptor)
    for candidate_id, arm, descriptor in CANDIDATE_SPECS
}
CONTROL_ARMS = (
    ("wirehair2_head", "wirehair2_certified"),
    ("wirehair1", "wirehair1"),
)
RAW_CONTROL_ARM = (
    "wirehair2_raw_d12_h12_periodic", "wirehair2_experiment")
RAW_CONTROL_DESCRIPTOR_SHA256 = \
    "739092a7824449e6168f08b46661dfbe8ad5495ea4166b36073c79cd3bacdd11"
RAW_SEED_BASIS = "uniform-raw-v1"
RAW_SEED_SCHEDULE_SHA256 = \
    "90a98a3db207852dabdf5fb27573ef48bce52e0228cee4e291d96fa44ed509a7"
CANDIDATE_DESCRIPTION_SCHEMA = "wirehair.wh2.native-worker-description.v3"
CANDIDATE_DESCRIPTION_ARM_FIELDS = frozenset((
    "arm", "codec", "arm_descriptor_sha256",
    "construction_seed_basis", "dense_anchor_layout",
    "seed_schedule_sha256",
))
CONTROL_DESCRIPTOR_SHA256S = (
    "4cafe27a8fb388ca9a4249b2c279b1406e7a0a86bcf14e98246988c7c503fa7a",
    "d5a24d404e69efeb439907cd8271eba98d6af86b58efe159a820fb7aea08883d",
)

CAMPAIGN_SUMMARY_SCHEMA = "wirehair.wh2.native-recovery-screen-run.v3"
COMBINATION_SUMMARY_SCHEMA = \
    "wirehair.wh2.logical-recovery-combination.v3"
CAMPAIGN_SUMMARY_FIELDS = frozenset((
    "schema", "status", "output_dir", "candidate_id", "candidate_arm",
    "source_git_commit", "contract_sha256", "domain_sha256",
    "trace_manifest_sha256", "worker_binary_sha256", "controller_cpu",
    "worker_cpus", "recovery_records", "recovery_freeze_sha256",
    "recovery_result_sha256", "recovery_execution_receipt_sha256",
    "thermal_samples", "cpu_tctl_max_millic", "dimm_max_millic",
    "construction_seed_basis", "seed_schedule_sha256", "summary_sha256",
    "timing_proxy_witness_sha256", "work_rank_summary_sha256",
    "work_rank_result_stream_sha256", "work_rank_domain_sha256",
    "raw_identity_join_count", "raw_identity_join_sha256",
))
COMBINATION_SUMMARY_FIELDS = frozenset((
    "schema", "status", "artifact_kind", "is_execution_receipt", "phase",
    "contract_sha256", "domain_sha256", "trace_manifest_sha256",
    "trace_file_sha256", "source_git_commit", "worker_binary_sha256",
    "candidate_roster", "arm_roster", "campaigns",
    "logical_freeze_manifest_sha256", "logical_result_stream_sha256",
    "construction_seed_basis", "seed_schedule_sha256", "validator_summary",
    "validator_summary_sha256", "combination_sha256",
    "work_rank_summary_sha256", "work_rank_result_stream_sha256",
    "work_rank_domain_sha256", "raw_identity_join_count",
    "raw_identity_join_sha256",
))
CAMPAIGN_BINDING_FIELDS = frozenset((
    "candidate_id", "candidate_arm", "run_summary_sha256",
    "freeze_manifest_sha256", "result_stream_sha256",
    "execution_receipt_sha256",
))
ZERO_SHA256 = "0" * 64
RECOVERY_RECORDS = 1440
LOGICAL_RECORDS = 1440
RECOVERY_WORKER_COUNT = 8
MAX_RECOVERY_RECORD_BYTES = 4096
MAX_RECOVERY_METADATA_BYTES = 1024 * 1024
RECOVERY_NATIVE_STREAM_BYTE_CAP = \
    RECOVERY_RECORDS * MAX_RECOVERY_RECORD_BYTES
RECOVERY_RESULT_STREAM_BYTE_CAP = \
    RECOVERY_RECORDS * MAX_RECOVERY_RECORD_BYTES
LOGICAL_RESULT_STREAM_BYTE_CAP = \
    LOGICAL_RECORDS * MAX_RECOVERY_RECORD_BYTES
RAW_IDENTITY_JOIN_COUNT = 720
TIMING_PROXY_WITNESS_SCHEMA = \
    "wirehair.wh2.native-timing-proxy-witness.v2"
TIMING_PROXY_DOMAIN = b"wirehair.wh2.timing-proxy-domain.v2\0"
TIMING_PROXY_FIELDS = frozenset((
    "schema", "source_git_commit", "binary_sha256", "proof_scope",
    "evidence_phase", "construction_attempts", "applicability",
    "production_timing_proxy_arm",
    "production_timing_proxy_arm_descriptor_sha256",
    "raw_recovery_reference_arm",
    "raw_recovery_reference_arm_descriptor_sha256",
    "raw_recovery_seed_basis", "raw_recovery_seed_schedule_sha256",
    "timing_candidate_arm", "timing_candidate_arm_descriptor_sha256",
    "timing_seed_basis", "timing_seed_schedule_sha256",
    "timing_seed_policy_arms", "seed_relationship",
    "witness_domain_sha256", "cells",
))
TIMING_PROXY_CELL_FIELDS = frozenset((
    "K", "block_bytes", "construction_attempt",
    "normalized_structure_sha256",
    "production_timing_configuration_sha256",
    "production_timing_equation_system_sha256",
    "production_timing_precode_seed", "production_timing_packet_seed",
    "raw_recovery_precode_seed", "raw_recovery_packet_seed",
    "seeds_differ",
))
RECOVERY_COMPLETED_ARTIFACT_BYTE_LIMITS = {
    "run-summary.json": MAX_RECOVERY_METADATA_BYTES,
    "recovery-freeze.json": MAX_RECOVERY_METADATA_BYTES,
    "recovery-traces.jsonl":
        360 * MAX_RECOVERY_RECORD_BYTES,
    "recovery-native-results.jsonl": RECOVERY_NATIVE_STREAM_BYTE_CAP,
    "recovery-results.jsonl": RECOVERY_RESULT_STREAM_BYTE_CAP,
    "recovery-execution.json": MAX_RECOVERY_METADATA_BYTES,
    "timing-proxy-witness.json": MAX_RECOVERY_METADATA_BYTES,
    "sampler-attestation.json": MAX_RECOVERY_METADATA_BYTES,
}
MAX_COMPLETED_ARTIFACT_BYTES = max(
    RECOVERY_COMPLETED_ARTIFACT_BYTE_LIMITS.values())
MAX_COMPLETED_CAMPAIGN_BYTES = sum(
    RECOVERY_COMPLETED_ARTIFACT_BYTE_LIMITS.values())
RECOVERY_COMPLETED_SOURCE_NAMES = {
    "summary": "run-summary.json",
    "freeze": "recovery-freeze.json",
    "trace": "recovery-traces.jsonl",
    "native": "recovery-native-results.jsonl",
    "result": "recovery-results.jsonl",
    "receipt": "recovery-execution.json",
    "witness": "timing-proxy-witness.json",
    "sampler": "sampler-attestation.json",
}
RECOVERY_COMPLETED_DEPENDENCY_NAMES = {
    key: name for key, name in RECOVERY_COMPLETED_SOURCE_NAMES.items()
    if key != "summary"
}
RECOVERY_COMPLETED_CONTEXTS = {
    "summary": "campaign run summary",
    "freeze": "campaign recovery freeze",
    "trace": "campaign trace manifest",
    "native": "campaign native recovery stream",
    "result": "campaign recovery ledger",
    "receipt": "campaign recovery execution receipt",
    "witness": "campaign timing proxy witness",
    "sampler": "campaign sampler attestation",
}


class RecoveryRunnerError(RuntimeError):
    """The recovery-only campaign or logical combination is invalid."""


def fail(message: str) -> None:
    raise RecoveryRunnerError(message)


def _close_recovery_descriptors(
        descriptors: Sequence[int], context: str,
        primary: Optional[BaseException] = None) -> None:
    """Exhaust closes, ignoring OSError but preserving unexpected control flow."""
    unexpected: List[BaseException] = []
    for descriptor in descriptors:
        if type(descriptor) is not int or descriptor < 0:
            continue
        try:
            os.close(descriptor)
        except OSError:
            # The inode/name transaction is already decided; a close(2) error
            # cannot be repaired and must not revoke a published marker.
            pass
        except BaseException as exc:
            unexpected.append(exc)
    if not unexpected:
        return
    messages = []
    control_flow: Optional[BaseException] = None
    if primary is not None and not isinstance(primary, Exception):
        control_flow = primary
    elif primary is not None:
        messages.append("primary failure: {}".format(primary))
    messages.append("{} raised {} unexpected close failure(s): {}".format(
        context, len(unexpected), "; ".join(
            "{}: {}".format(type(exc).__name__, exc)
            for exc in unexpected)))
    details = RecoveryRunnerError("; ".join(messages))
    if control_flow is None:
        control_flow = next(
            (exc for exc in unexpected if not isinstance(exc, Exception)), None)
    if control_flow is not None:
        raise control_flow from details
    if primary is not None:
        raise details from primary
    raise details from unexpected[0]


def _candidate(candidate_id: Any) -> Tuple[str, str]:
    if not isinstance(candidate_id, str) or candidate_id not in CANDIDATE_BY_ID:
        fail("candidate ID is outside the closed recovery-only roster")
    return CANDIDATE_BY_ID[candidate_id]


def _expected_arms(candidate_id: str) \
        -> Tuple[Tuple[str, str], ...]:
    candidate_arm, _ = _candidate(candidate_id)
    return CONTROL_ARMS + (RAW_CONTROL_ARM,) + (
        (candidate_arm, "wirehair2_experiment"),)


def _is_sha256(value: Any) -> bool:
    return isinstance(value, str) and \
        contract_api.SHA256.fullmatch(value) is not None


def _hash_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _hash_jsonl(values: Sequence[Mapping[str, Any]]) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(contract_api.canonical_json(value).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def _self_hash(value: Mapping[str, Any], field: str) -> str:
    return contract_api.sha256_json({
        key: item for key, item in value.items() if key != field
    })


def describe_candidate_worker(
        worker_path: Path, candidate_id: str, deadline: float,
        ) -> Mapping[str, Any]:
    """Describe one exact four-arm worker configuration and bind its file."""
    candidate_arm, expected_descriptor = _candidate(candidate_id)
    try:
        resolved = worker_path.resolve(strict=True)
        info = resolved.stat()
    except OSError as exc:
        fail("cannot resolve native worker {}: {}".format(worker_path, exc))
    if not stat.S_ISREG(info.st_mode) or not os.access(str(resolved), os.X_OK):
        fail("native worker must be an executable regular file")
    binary_sha256 = runner_api._hash_file(resolved)
    description = runner_api._parse_canonical_line(
        runner_api._run_command(
            [str(resolved), "--describe-recovery-candidate", candidate_id],
            deadline, "recovery candidate worker description"),
        "recovery candidate worker description")
    if (set(description) != runner_api.DESCRIPTION_FIELDS or
            description.get("schema") != CANDIDATE_DESCRIPTION_SCHEMA or
            not isinstance(description.get("source_git_commit"), str) or
            re.fullmatch(r"[0-9a-f]{40}",
                         description["source_git_commit"]) is None or
            description.get("binary_sha256") != binary_sha256):
        fail("recovery candidate description does not bind its executable")
    arms = description.get("arms")
    expected = _expected_arms(candidate_id)
    if not isinstance(arms, list) or len(arms) != len(expected):
        fail("recovery candidate description is not an exact four-arm roster")
    for index, expected_arm in enumerate(expected):
        arm = arms[index]
        if (not isinstance(arm, dict) or
                set(arm) != CANDIDATE_DESCRIPTION_ARM_FIELDS or
                arm.get("arm") != expected_arm[0] or
                arm.get("codec") != expected_arm[1] or
                not _is_sha256(arm.get("arm_descriptor_sha256"))):
            fail("recovery candidate description has an invalid arm at index {}"
                 .format(index))
    expected_policies = (
        (contract_api.PRODUCTION_CONSTRUCTION_SEED_BASIS, ZERO_SHA256),
        (contract_api.NOT_APPLICABLE_CONSTRUCTION_SEED_BASIS, ZERO_SHA256),
        (RAW_SEED_BASIS, RAW_SEED_SCHEDULE_SHA256),
        (RAW_SEED_BASIS, RAW_SEED_SCHEDULE_SHA256),
    )
    if tuple((arm["construction_seed_basis"],
              arm["seed_schedule_sha256"]) for arm in arms) != \
            expected_policies:
        fail("recovery candidate description has an invalid seed-policy roster")
    if tuple(arm["dense_anchor_layout"] for arm in arms) != (
            "disabled", "not-applicable", "disabled", "two07"):
        fail("recovery candidate description has an invalid layout roster")
    if tuple(arms[index]["arm_descriptor_sha256"] for index in (0, 1)) != \
            CONTROL_DESCRIPTOR_SHA256S:
        fail("recovery candidate description substitutes a control descriptor")
    if (arms[2]["arm"] != RAW_CONTROL_ARM[0] or
            arms[2]["arm_descriptor_sha256"] !=
            RAW_CONTROL_DESCRIPTOR_SHA256):
        fail("recovery candidate description substitutes its raw control")
    if arms[3]["arm"] != candidate_arm or \
            arms[3]["arm_descriptor_sha256"] != expected_descriptor:
        fail("candidate descriptor does not match its closed candidate ID")
    if len({arm["arm_descriptor_sha256"] for arm in arms}) != len(arms):
        fail("recovery candidate description reuses an arm descriptor")
    result = dict(description)
    result["resolved_path"] = str(resolved)
    result["candidate_id"] = candidate_id
    return result


def load_timing_proxy_witness(
        description: Mapping[str, Any], source_git_commit: str,
        deadline: float) -> Mapping[str, Any]:
    """Validate the attempt-0 structure-only production timing witness."""
    worker = description.get("resolved_path")
    if not isinstance(worker, str) or not worker:
        fail("worker description lacks its resolved executable path")
    try:
        witness = runner_api._parse_canonical_line(
            runner_api._run_command(
                [worker, "--emit-timing-proxy-witness"], deadline,
                "D12 timing proxy witness"),
            "D12 timing proxy witness")
    except runner_api.RunnerError as exc:
        fail(str(exc))
    return _validate_timing_proxy_witness(
        witness, description.get("binary_sha256"), source_git_commit)


def _validate_timing_proxy_witness(
        witness: Mapping[str, Any], binary_sha256: Any,
        source_git_commit: Any) -> Mapping[str, Any]:
    """Pure validator used both at creation and immutable campaign reopen."""
    if (set(witness) != TIMING_PROXY_FIELDS or
            witness.get("schema") != TIMING_PROXY_WITNESS_SCHEMA or
            witness.get("source_git_commit") != source_git_commit or
            witness.get("binary_sha256") != binary_sha256 or
            witness.get("proof_scope") !=
                "d12-disabled-structure-under-production-timing-seed-policy-v1" or
            witness.get("evidence_phase") != "development" or
            not isinstance(witness.get("construction_attempts"), list) or
            len(witness["construction_attempts"]) != 1 or
            type(witness["construction_attempts"][0]) is not int or
            witness["construction_attempts"][0] != 0 or
            witness.get("applicability") !=
                "development-attempt-0-only-new-witness-required-for-other-attempt-semantics" or
            witness.get("raw_recovery_reference_arm") !=
                RAW_CONTROL_ARM[0] or
            witness.get("raw_recovery_reference_arm_descriptor_sha256") !=
                RAW_CONTROL_DESCRIPTOR_SHA256 or
            witness.get("raw_recovery_seed_basis") != RAW_SEED_BASIS or
            witness.get("raw_recovery_seed_schedule_sha256") !=
                RAW_SEED_SCHEDULE_SHA256 or
            witness.get("production_timing_proxy_arm") !=
                "wirehair2_head" or
            witness.get(
                "production_timing_proxy_arm_descriptor_sha256") !=
                CONTROL_DESCRIPTOR_SHA256S[0] or
            witness.get("timing_candidate_arm") != CANDIDATE_SPECS[0][1] or
            witness.get("timing_candidate_arm_descriptor_sha256") !=
                CANDIDATE_SPECS[0][2] or
            witness.get("timing_seed_basis") !=
                contract_api.PRODUCTION_CONSTRUCTION_SEED_BASIS or
            witness.get("timing_seed_schedule_sha256") != ZERO_SHA256 or
            witness.get("timing_seed_policy_arms") != [
                "wirehair2_head", CANDIDATE_SPECS[0][1]] or
            witness.get("seed_relationship") !=
                "raw-recovery-precode-and-packet-seeds-differ-from-production-timing"):
        fail("D12 structure-only timing witness substitutes its proof scope")
    expected_coordinates = [
        {"K": K, "block_bytes": block_bytes,
         "construction_attempt": 0}
        for block_bytes in (64, 1280)
        for K in contract_api.EXPECTED_TIMING_SHORT_K
    ]
    expected_domain_sha256 = hashlib.sha256(
        TIMING_PROXY_DOMAIN + contract_api.canonical_json(
            expected_coordinates).encode("utf-8")).hexdigest()
    cells = witness.get("cells")
    if (not isinstance(cells, list) or len(cells) != 24 or
            witness.get("witness_domain_sha256") !=
                expected_domain_sha256):
        fail("D12 timing proxy witness has the wrong frozen domain")
    observed_coordinates = []
    for index, cell in enumerate(cells):
        if (not isinstance(cell, dict) or
                set(cell) != TIMING_PROXY_CELL_FIELDS or
                any(not _is_sha256(cell.get(field)) for field in (
                    "normalized_structure_sha256",
                    "production_timing_configuration_sha256",
                    "production_timing_equation_system_sha256")) or
                not isinstance(cell.get("K"), int) or
                isinstance(cell.get("K"), bool) or
                not isinstance(cell.get("block_bytes"), int) or
                isinstance(cell.get("block_bytes"), bool) or
                type(cell.get("construction_attempt")) is not int or
                cell["construction_attempt"] != 0 or
                cell.get("seeds_differ") is not True or
                not isinstance(cell.get("production_timing_precode_seed"), str) or
                contract_api.HEX64.fullmatch(
                    cell["production_timing_precode_seed"]) is None or
                not isinstance(cell.get("raw_recovery_precode_seed"), str) or
                contract_api.HEX64.fullmatch(
                    cell["raw_recovery_precode_seed"]) is None or
                not isinstance(cell.get("production_timing_packet_seed"), str) or
                contract_api.HEX32.fullmatch(
                    cell["production_timing_packet_seed"]) is None or
                not isinstance(cell.get("raw_recovery_packet_seed"), str) or
                contract_api.HEX32.fullmatch(
                    cell["raw_recovery_packet_seed"]) is None or
                cell["raw_recovery_precode_seed"] !=
                    contract_api._effective_raw_precode_seed(0) or
                cell["raw_recovery_packet_seed"] !=
                    contract_api._effective_raw_packet_seed(0) or
                cell["production_timing_precode_seed"] ==
                    cell["raw_recovery_precode_seed"] or
                cell["production_timing_packet_seed"] ==
                    cell["raw_recovery_packet_seed"]):
            fail("D12 timing proxy witness cell {} is malformed".format(
                index))
        observed_coordinates.append({
            "K": cell["K"], "block_bytes": cell["block_bytes"],
            "construction_attempt": cell["construction_attempt"],
        })
    if observed_coordinates != expected_coordinates:
        fail("D12 timing proxy witness omits or reorders a frozen cell")
    return witness


def _candidate_commands(
        description: Mapping[str, Any], candidate_id: str,
        cpus: Sequence[int]) -> List[List[str]]:
    _candidate(candidate_id)
    worker = description.get("resolved_path")
    if not isinstance(worker, str) or not worker:
        fail("worker description lacks its resolved executable path")
    return [
        [worker, "--describe-recovery-candidate", candidate_id],
        [worker, "--emit-traces", "recovery"],
    ] + [
        [worker, "--recovery-candidate-worker", candidate_id, str(cpu)]
        for cpu in cpus
    ] + [[worker, "--emit-timing-proxy-witness"]]


def _candidate_freeze(
        contract: Mapping[str, Any], description: Mapping[str, Any],
        candidate_id: str, cpus: Sequence[int], controller_cpu: int,
        source_commit: str, trace_sha256: str,
        timing_proxy_witness_sha256: str,
        work_binding: Mapping[str, Any]) -> Mapping[str, Any]:
    expected = _expected_arms(candidate_id)
    cpu_list = list(cpus)
    if (len(cpu_list) != RECOVERY_WORKER_COUNT or
            any(type(cpu) is not int or cpu < 0 for cpu in cpu_list) or
            cpu_list != sorted(set(cpu_list))):
        fail("recovery-only screen requires eight sorted unique worker CPUs")
    roster = [arm for arm, _ in expected]
    arms = []
    for index, value in enumerate(description["arms"]):
        if value["arm"] != roster[index]:
            fail("worker description changed before recovery freeze")
        arms.append({
            "arm": value["arm"],
            "codec": value["codec"],
            "binary_sha256": description["binary_sha256"],
            "arm_descriptor_sha256": value["arm_descriptor_sha256"],
            "construction_policy": "not_applicable"
                if value["codec"] == "wirehair1" else "raw_base",
            "repair_map_sha256": ZERO_SHA256,
            "construction_seed_basis": value["construction_seed_basis"],
            "dense_anchor_layout": value["dense_anchor_layout"],
            "seed_schedule_sha256": value["seed_schedule_sha256"],
        })
    return {
        "schema": contract_api.RAW_FREEZE_SCHEMA,
        "contract_sha256": contract_api.contract_sha256(contract),
        "evidence_kind": "recovery",
        "phase": "development",
        "domain_sha256": contract["recovery"]["domains"]["development"][
            "domain_sha256"],
        "source_git_commit": source_commit,
        "arm_roster": roster,
        "arm_roster_sha256": contract_api.arm_roster_sha256(roster),
        "trace_manifest_sha256": trace_sha256,
        "repair_training_trace_manifest_sha256": ZERO_SHA256,
        "commands": _candidate_commands(
            description, candidate_id, cpu_list),
        "cpu_affinity": cpu_list,
        "host_identity": runner_api._host_identity(controller_cpu),
        "architecture_roles": dict(
            contract_api.EXPECTED_RAW_ARCHITECTURE_ROLES),
        "timing_proxy_witness_sha256": timing_proxy_witness_sha256,
        "work_rank_summary_sha256":
            work_binding["work_rank_summary_sha256"],
        "work_rank_result_stream_sha256":
            work_binding["work_rank_result_stream_sha256"],
        "work_rank_domain_sha256":
            work_binding["work_rank_domain_sha256"],
        "arms": arms,
    }


def _validate_candidate_freeze(
        freeze: Mapping[str, Any], candidate_id: str,
        worker_path: Optional[str] = None) -> None:
    candidate_arm, expected_descriptor = _candidate(candidate_id)
    expected_roster = [arm for arm, _ in _expected_arms(candidate_id)]
    if freeze.get("arm_roster") != expected_roster or \
            freeze.get("evidence_kind") != "recovery" or \
            freeze.get("phase") != "development" or \
            freeze.get("schema") != contract_api.RAW_FREEZE_SCHEMA or \
            freeze.get("architecture_roles") != \
                contract_api.EXPECTED_RAW_ARCHITECTURE_ROLES:
        fail("candidate freeze does not bind the exact four-arm recovery roster")
    arms = freeze.get("arms")
    if (not isinstance(arms, list) or len(arms) != 4 or
            any(not isinstance(arm, dict) for arm in arms) or
            tuple((arm.get("arm"), arm.get("codec")) for arm in arms) !=
            _expected_arms(candidate_id) or
            arms[2].get("arm") != RAW_CONTROL_ARM[0] or
            arms[2].get("arm_descriptor_sha256") !=
            RAW_CONTROL_DESCRIPTOR_SHA256 or
            arms[3].get("arm") != candidate_arm or
            arms[3].get("arm_descriptor_sha256") != expected_descriptor):
        fail("candidate freeze substitutes its candidate descriptor")
    if tuple(arms[index].get("arm_descriptor_sha256")
             for index in (0, 1)) != CONTROL_DESCRIPTOR_SHA256S:
        fail("candidate freeze substitutes a control descriptor")
    if tuple(arm.get("dense_anchor_layout") for arm in arms) != (
            "disabled", "not-applicable", "disabled", "two07"):
        fail("candidate freeze substitutes a dense-anchor layout")
    for field in (
            "timing_proxy_witness_sha256", "work_rank_summary_sha256",
            "work_rank_result_stream_sha256", "work_rank_domain_sha256"):
        if not _is_sha256(freeze.get(field)) or \
                freeze.get(field) == ZERO_SHA256:
            fail("candidate freeze lacks its nonzero {} binding".format(
                field))
    binaries = {arm.get("binary_sha256") for arm in arms}
    if len(binaries) != 1 or not _is_sha256(next(iter(binaries))):
        fail("candidate freeze does not use one bound worker binary")
    commands = freeze.get("commands")
    cpus = freeze.get("cpu_affinity")
    if (not isinstance(commands, list) or not isinstance(cpus, list) or
            len(cpus) != RECOVERY_WORKER_COUNT):
        fail("candidate freeze lacks its exact eight-worker argv roster")
    if worker_path is None:
        try:
            worker_path = commands[0][0]
        except (IndexError, TypeError):
            fail("candidate freeze lacks its description argv")
    description_stub = {"resolved_path": worker_path}
    if commands != _candidate_commands(description_stub, candidate_id, cpus):
        fail("candidate freeze commands differ from the exact candidate argv")


def _raw_policy_from_freeze(
        freeze: Mapping[str, Any]) -> Tuple[str, str]:
    """Return the one worker-attested raw policy bound by a v3 freeze."""
    if freeze.get("schema") != contract_api.RAW_FREEZE_SCHEMA:
        fail("candidate evidence does not use the raw v3 freeze schema")
    arms = freeze.get("arms")
    if not isinstance(arms, list):
        fail("candidate freeze arm records are malformed")
    policies = {
        (arm.get("construction_seed_basis"),
         arm.get("seed_schedule_sha256"))
        for arm in arms
        if isinstance(arm, dict) and
        arm.get("construction_seed_basis") == RAW_SEED_BASIS
    }
    if policies != {(RAW_SEED_BASIS, RAW_SEED_SCHEDULE_SHA256)}:
        fail("candidate freeze does not bind one uniform raw seed policy")
    return next(iter(policies))


def write_recovery_freeze(
        contract: Mapping[str, Any], description: Mapping[str, Any],
        candidate_id: str, cpus: Sequence[int], controller_cpu: int,
        source_commit: str, trace_sha256: str, output_dir: Path,
        timing_proxy_witness_sha256: str,
        work_binding: Mapping[str, Any],
        ) -> Mapping[str, Any]:
    path = output_dir / "recovery-freeze.json"
    runner_api._atomic_write_object(path, _candidate_freeze(
        contract, description, candidate_id, cpus, controller_cpu,
        source_commit, trace_sha256, timing_proxy_witness_sha256,
        work_binding))
    try:
        loaded = contract_api.load_freeze_manifest(
            contract, "development", path, "recovery")
    except contract_api.ContractError as exc:
        fail(str(exc))
    _validate_candidate_freeze(
        loaded, candidate_id, description["resolved_path"])
    return loaded


def spawn_candidate_workers(
        description: Mapping[str, Any], candidate_id: str,
        cpus: Sequence[int], deadline: float,
        ) -> List[runner_api.PersistentWorker]:
    """Spawn workers in the recovery-only mode; that mode rejects every T job."""
    _candidate(candidate_id)
    workers: List[runner_api.PersistentWorker] = []
    provisional_worker: Optional[Any] = None
    try:
        for cpu in cpus:
            argv = [description["resolved_path"],
                    "--recovery-candidate-worker", candidate_id, str(cpu)]
            # Allocate the owner before Popen.  Once Popen returns, retain the
            # child here until PersistentWorker construction and list handoff
            # both finish, so every observable interruption has an exhaustive
            # cleanup owner.
            provisional_worker = SimpleNamespace(
                cpu=cpu, process=None, start_ticks=0, pending=None, buffer=b"")
            try:
                provisional_worker.process = subprocess.Popen(
                    argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE, bufsize=0)
            except OSError as exc:
                fail("cannot spawn recovery worker on CPU {}: {}".format(
                    cpu, exc))
            process = provisional_worker.process
            worker = runner_api.PersistentWorker(cpu, process, 0)
            workers.append(worker)
            provisional_worker = None
        pending = set(range(len(workers)))
        ready_deadline = min(deadline, time.monotonic() + 10.0)
        while pending:
            if time.monotonic() >= ready_deadline:
                fail("recovery workers did not establish singleton affinity")
            for index in list(pending):
                worker = workers[index]
                if worker.process.poll() is not None:
                    fail("recovery worker on CPU {} exited during startup: {}"
                         .format(worker.cpu,
                                 runner_api._worker_stderr(worker)))
                try:
                    ticks = native_api._parse_proc_start_ticks(
                        worker.process.pid)
                    affinity = os.sched_getaffinity(worker.process.pid)
                except (native_api.NativeEvidenceError,
                        AttributeError, OSError):
                    continue
                if affinity == {worker.cpu}:
                    worker.start_ticks = ticks
                    os.set_blocking(worker.process.stdout.fileno(), False)
                    os.set_blocking(worker.process.stderr.fileno(), False)
                    pending.remove(index)
            if pending:
                time.sleep(0.01)
        return workers
    except BaseException:
        cleanup_workers = list(workers)
        if (provisional_worker is not None and
                provisional_worker.process is not None and
                all(worker.process is not provisional_worker.process
                    for worker in cleanup_workers)):
            cleanup_workers.append(provisional_worker)
        _finish_recovery_cleanup(
            cleanup_workers, False, False, set(), sys.exc_info()[1])
        raise


def _candidate_recovery_jobs(
        contract: Mapping[str, Any]) -> List[runner_api.Job]:
    """Build the exact four-arm job ledger without weakening legacy gates."""
    cells = contract["recovery"]["domains"]["development"][
        "expected_cells_per_arm"]
    if type(cells) is not int or cells != 360:
        fail("development recovery domain is not exactly 360 frozen cells")
    jobs = [
        runner_api.Job("recovery", cell, arm, cell * 4 + arm)
        for cell in range(cells) for arm in range(4)
    ]
    if len(jobs) != RECOVERY_RECORDS or any(
            job.ordinal != index for index, job in enumerate(jobs)):
        fail("candidate recovery jobs differ from the four-arm frozen roster")
    return jobs


def _run_recovery_jobs(
        contract: Mapping[str, Any], freeze: Mapping[str, Any],
        description: Mapping[str, Any],
        workers: Sequence[runner_api.PersistentWorker], output_dir: Path,
        window_start_ns: int, deadline: float) -> Tuple[Path, int]:
    if len(workers) != RECOVERY_WORKER_COUNT:
        fail("recovery-only screen requires exactly eight workers")
    jobs = _candidate_recovery_jobs(contract)
    if len(jobs) != RECOVERY_RECORDS or \
            any(not job.command().startswith(b"R ") for job in jobs):
        fail("recovery-only job roster is not exactly 1440 R commands")
    path = output_dir / "recovery-native-results.jsonl"
    strict_validator = runner_api._strict_response_validator(
        contract, freeze, "recovery", description, window_start_ns)
    sink = runner_api.AtomicLineSink(
        path, maximum_bytes=RECOVERY_NATIVE_STREAM_BYTE_CAP)

    def bounded_validator(
            value: Mapping[str, Any], line: bytes,
            worker: runner_api.PersistentWorker, job: runner_api.Job) -> int:
        if len(line) > MAX_RECOVERY_RECORD_BYTES:
            fail("native recovery response exceeds the 4096-byte record cap")
        return strict_validator(value, line, worker, job)

    try:
        maximum_end, used = runner_api.run_job_batch(
            workers, jobs, 0, sink, deadline,
            bounded_validator,
            maximum_response_bytes=MAX_RECOVERY_RECORD_BYTES)
        if used != {worker.cpu for worker in workers}:
            fail("recovery campaign did not exercise every frozen CPU")
        sink.publish()
        return path, maximum_end
    finally:
        sink.abort()


def _finish_recovery_cleanup(
        workers: Sequence[runner_api.PersistentWorker], clean_shutdown: bool,
        controller_pinned: bool, original_affinity: Set[int],
        primary: Optional[BaseException]) -> None:
    """Attempt every cleanup and preserve ordinary and control-flow failures."""
    cleanup_failures: List[Tuple[str, BaseException]] = []
    if workers and not clean_shutdown:
        try:
            runner_api.terminate_workers(workers)
        except BaseException as cleanup:
            cleanup_failures.append((
                "recovery worker cleanup failed", cleanup))
    if controller_pinned:
        try:
            runner_api._restore_controller_affinity(original_affinity)
        except BaseException as cleanup:
            cleanup_failures.append((
                "controller affinity cleanup failed", cleanup))
    if not cleanup_failures:
        return

    cleanup_message = "; ".join(
        "{}: {}".format(label, failure)
        for label, failure in cleanup_failures)
    control_flow: Optional[BaseException]
    if primary is not None and not isinstance(primary, Exception):
        control_flow = primary
    else:
        control_flow = next((
            failure for _label, failure in cleanup_failures
            if not isinstance(failure, Exception)), None)
    if control_flow is not None:
        details = []
        if primary is not None and primary is not control_flow:
            details.append("primary failure: {}".format(primary))
        details.append(cleanup_message)
        raise control_flow from RecoveryRunnerError("; ".join(details))
    if primary is not None:
        raise RecoveryRunnerError("{}; {}".format(
            primary, cleanup_message)) from primary
    raise RecoveryRunnerError(cleanup_message) from cleanup_failures[0][1]


def run_recovery_screen(
        args: argparse.Namespace,
        finish_hard_wall: Optional[Callable[[], None]] = None,
        ) -> Mapping[str, Any]:
    candidate_id = args.candidate
    candidate_arm, _ = _candidate(candidate_id)
    if (not math.isfinite(args.deadline_seconds) or
            not 0.0 < args.deadline_seconds <= runner_api.MAX_WALL_SECONDS):
        fail("--deadline-seconds must be in (0,7200]")
    deadline = time.monotonic() + args.deadline_seconds
    contract = contract_api.load_contract(args.contract)
    try:
        work_screen = work_api.load_completed_work_screen(
            contract, args.work_rank_dir)
    except work_api.WorkScreenError as exc:
        fail(str(exc))
    work_summary = work_screen["summary"]
    work_binding = {
        "work_rank_summary_sha256": work_summary["summary_sha256"],
        "work_rank_result_stream_sha256":
            work_summary["result_stream_sha256"],
        "work_rank_domain_sha256": work_summary["work_domain_sha256"],
    }
    description = describe_candidate_worker(args.worker, candidate_id, deadline)
    try:
        original_affinity = set(os.sched_getaffinity(0))
    except (AttributeError, OSError) as exc:
        fail("cannot inspect initial controller affinity: {}".format(exc))
    explicit = runner_api.parse_cpu_list(args.cpus) \
        if args.cpus is not None else None
    cpus, controller_cpu = runner_api.select_cpu_layout(
        args.sampler_cpu, explicit, args.controller_cpu,
        affinity=original_affinity)
    runner_api._preflight_sampler(
        args.sampler_pid, args.sampler_cpu,
        args.sampler_script, args.sampler_csv)
    source_commit = runner_api._git_head(deadline)
    runner_api._require_worker_source_commit(description, source_commit)
    if work_screen.get("source_git_commit") != source_commit:
        fail("native and work/rank evidence use different source commits")
    timing_proxy_witness = load_timing_proxy_witness(
        description, source_commit, deadline)
    timing_proxy_witness_sha256 = contract_api.sha256_json(
        timing_proxy_witness)
    output_dir = runner_api._create_output_dir(args.output_dir)
    runner_api._atomic_write_object(
        output_dir / "timing-proxy-witness.json", timing_proxy_witness)

    trace_path, trace_sha256 = runner_api._emit_and_assemble_trace(
        contract, "recovery", description, output_dir, deadline)
    if runner_api._git_head(deadline) != source_commit:
        fail("codec source commit changed before the recovery freeze")
    freeze = write_recovery_freeze(
        contract, description, candidate_id, cpus, controller_cpu,
        source_commit, trace_sha256, output_dir,
        timing_proxy_witness_sha256, work_binding)
    freeze_path = output_dir / "recovery-freeze.json"

    workers: List[runner_api.PersistentWorker] = []
    clean_shutdown = False
    controller_pinned = False
    completed_summary: Optional[Mapping[str, Any]] = None
    try:
        workers = spawn_candidate_workers(
            description, candidate_id, cpus, deadline)
        controller_pinned = True
        runner_api._pin_controller(controller_cpu)
        window_start_ns = runner_api.choose_new_sampler_start(
            args.sampler_csv, deadline)
        native_path, maximum_worker_end = _run_recovery_jobs(
            contract, freeze, description, workers, output_dir,
            window_start_ns, deadline)
        window_end_ns = runner_api._wait_for_sampler_sample(
            args.sampler_csv, deadline,
            at_or_after_ns=maximum_worker_end)
        if runner_api._git_head(deadline) != source_commit:
            fail("codec source commit changed during the recovery screen")
        sampler_path = output_dir / "sampler-attestation.json"
        try:
            native_api.write_sampler_attestation(
                args.sampler_pid, args.sampler_cpu,
                args.sampler_script, args.sampler_csv,
                window_start_ns, window_end_ns, sampler_path)
        except native_api.NativeEvidenceError as exc:
            fail(str(exc))

        result_path = output_dir / "recovery-results.jsonl"
        receipt_path = output_dir / "recovery-execution.json"
        try:
            assembled = native_api.assemble_results(
                contract, "recovery", "development", freeze_path,
                trace_path, native_path, sampler_path,
                result_path, receipt_path, verify_live_sampler=True)
            runner_api._remaining(deadline, "assembling recovery results")
            validated = native_api.validate_execution_receipt(
                contract, "recovery", "development", freeze_path,
                trace_path, native_path, result_path, receipt_path,
                verify_live_sampler=True, sampler_path=sampler_path)
            runner_api._remaining(deadline, "validating recovery execution")
        except (native_api.NativeEvidenceError,
                contract_api.ContractError) as exc:
            fail(str(exc))
        if assembled != validated:
            fail("recovery assembly and terminal validation disagree")
        runner_api.quit_workers(workers, deadline)
        clean_shutdown = True

        receipt = assembled["execution_receipt"]
        if receipt.get("record_count") != RECOVERY_RECORDS:
            fail("recovery execution receipt does not contain 1440 records")
        result_bytes = _read_regular_bytes(
            result_path, "campaign recovery ledger")
        if _hash_bytes(result_bytes) != receipt["result_stream_sha256"]:
            fail("campaign recovery ledger changed after terminal validation")
        result_rows = _parse_exact_jsonl(
            result_bytes, "campaign recovery ledger")
        joined_work = _bind_work_rank_identities(
            result_rows, work_screen, source_commit)
        if any(joined_work[field] != value
               for field, value in work_binding.items()):
            fail("joined work/rank evidence differs from its pre-result freeze")
        thermal = receipt["thermal"]
        raw_basis, raw_schedule_sha256 = _raw_policy_from_freeze(freeze)
        summary: Dict[str, Any] = {
            "schema": CAMPAIGN_SUMMARY_SCHEMA,
            "status": "complete",
            "output_dir": str(output_dir),
            "candidate_id": candidate_id,
            "candidate_arm": candidate_arm,
            "source_git_commit": source_commit,
            "contract_sha256": contract_api.contract_sha256(contract),
            "domain_sha256": contract["recovery"]["domains"]
                ["development"]["domain_sha256"],
            "trace_manifest_sha256": trace_sha256,
            "worker_binary_sha256": description["binary_sha256"],
            "controller_cpu": controller_cpu,
            "worker_cpus": list(cpus),
            "recovery_records": receipt["record_count"],
            "recovery_freeze_sha256": receipt["freeze_manifest_sha256"],
            "recovery_result_sha256": receipt["result_stream_sha256"],
            "recovery_execution_receipt_sha256": receipt["receipt_sha256"],
            "thermal_samples": thermal["sample_count"],
            "cpu_tctl_max_millic": thermal["cpu_tctl_max_millic"],
            "dimm_max_millic": thermal["dimm_max_millic"],
            "construction_seed_basis": raw_basis,
            "seed_schedule_sha256": raw_schedule_sha256,
            "timing_proxy_witness_sha256":
                timing_proxy_witness_sha256,
            **joined_work,
        }
        summary["summary_sha256"] = contract_api.sha256_json(summary)
        completed_summary = summary
    finally:
        _finish_recovery_cleanup(
            workers, clean_shutdown, controller_pinned, original_affinity,
            sys.exc_info()[1])
    if completed_summary is None:
        fail("recovery screen reached cleanup without a completed summary")
    if runner_api._git_head(deadline) != source_commit:
        fail("codec source commit changed before recovery completion")
    _commit_completed_recovery_screen(
        contract, output_dir, completed_summary, freeze, assembled,
        timing_proxy_witness, finish_hard_wall)
    # The successful hard-link is the final fallible action.  No work may be
    # added here without extending the pre-publication transaction around it.
    return completed_summary


def _load_canonical_object(path: Path, context: str) -> Mapping[str, Any]:
    return _parse_canonical_object(_read_regular_bytes(path, context), context)


def _parse_canonical_object(data: bytes, context: str) -> Mapping[str, Any]:
    try:
        value = contract_api._load_json_bytes(data, context)
    except contract_api.ContractError as exc:
        fail(str(exc))
    if not isinstance(value, dict):
        fail("{} must be a JSON object".format(context))
    logical = contract_api.canonical_json(value).encode("utf-8")
    if data != logical + b"\n":
        fail("{} must be exact canonical JSON followed by LF".format(context))
    return value


def _revalidate_terminal_timing_proxy_witness(
        path: Path, description: Mapping[str, Any], source_git_commit: str,
        expected_sha256: str) -> Mapping[str, Any]:
    """Reopen the published witness immediately before terminal status."""
    witness = _parse_canonical_object(
        _read_regular_bytes(path, "terminal timing proxy witness"),
        "terminal timing proxy witness")
    _validate_timing_proxy_witness(
        witness, description["binary_sha256"], source_git_commit)
    if contract_api.sha256_json(witness) != expected_sha256:
        fail("timing proxy witness changed before terminal summary publication")
    return witness


def _revalidate_terminal_campaign_execution(
        contract: Mapping[str, Any], candidate_id: str,
        expected_freeze: Mapping[str, Any],
        expected_validation: Mapping[str, Any], freeze_path: Path,
        trace_path: Path, native_path: Path, result_path: Path,
        receipt_path: Path, sampler_path: Path) -> Mapping[str, Any]:
    """Reopen every execution dependency immediately before completion."""
    try:
        validated = native_api.validate_execution_receipt(
            contract, "recovery", "development", freeze_path, trace_path,
            native_path, result_path, receipt_path,
            verify_live_sampler=False, sampler_path=sampler_path)
        reopened_freeze = contract_api.load_freeze_manifest(
            contract, "development", freeze_path, "recovery")
    except (native_api.NativeEvidenceError,
            contract_api.ContractError) as exc:
        fail(str(exc))
    _validate_candidate_freeze(reopened_freeze, candidate_id)
    if contract_api.canonical_json(reopened_freeze) != \
            contract_api.canonical_json(expected_freeze):
        fail("campaign freeze changed after native validation")
    if contract_api.canonical_json(validated) != \
            contract_api.canonical_json(expected_validation):
        fail("campaign execution changed after native validation")
    return validated


def _read_regular_bytes(
        path: Path, context: str, directory_fd: Optional[int] = None) -> bytes:
    """Read one bounded, already-opened regular artifact without symlink races."""
    descriptor = -1
    try:
        nofollow = getattr(os, "O_NOFOLLOW", 0)
        if nofollow == 0:
            fail("{} cannot be opened fail-closed without O_NOFOLLOW".format(
                context))
        flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | nofollow | \
            getattr(os, "O_NONBLOCK", 0)
        descriptor = os.open(str(path), flags, dir_fd=directory_fd)
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            fail("{} must be a regular non-symlink file".format(context))
        if before.st_size > MAX_COMPLETED_ARTIFACT_BYTES:
            fail("{} exceeds the bounded artifact size".format(context))
        chunks = []
        size = 0
        while True:
            block = os.read(descriptor, 1024 * 1024)
            if not block:
                break
            size += len(block)
            if size > MAX_COMPLETED_ARTIFACT_BYTES:
                fail("{} exceeds the bounded artifact size".format(context))
            chunks.append(block)
        after = os.fstat(descriptor)
        stable_before = (
            before.st_dev, before.st_ino, before.st_size,
            getattr(before, "st_mtime_ns", None),
            getattr(before, "st_ctime_ns", None))
        stable_after = (
            after.st_dev, after.st_ino, after.st_size,
            getattr(after, "st_mtime_ns", None),
            getattr(after, "st_ctime_ns", None))
        data = b"".join(chunks)
        if stable_before != stable_after or len(data) != before.st_size:
            fail("{} changed while it was being read".format(context))
        return data
    except OSError as exc:
        fail("cannot read {} {}: {}".format(context, path, exc))
    finally:
        if descriptor >= 0:
            os.close(descriptor)
    return b""


def _read_terminal_combination_artifacts(
        output_dir: Path) -> Mapping[str, bytes]:
    """Read the three named logical outputs through one pinned directory."""
    directory_fd = -1
    try:
        nofollow = getattr(os, "O_NOFOLLOW", 0)
        directory_flag = getattr(os, "O_DIRECTORY", 0)
        if nofollow == 0 or directory_flag == 0:
            fail("logical output directory cannot be opened fail-closed")
        directory_fd = os.open(
            str(output_dir), os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
            nofollow | directory_flag)
        opened_info = os.fstat(directory_fd)
        if not stat.S_ISDIR(opened_info.st_mode):
            fail("logical output path must be a real directory")
        resolved_info = output_dir.resolve(strict=True).stat()
        if (opened_info.st_dev, opened_info.st_ino) != \
                (resolved_info.st_dev, resolved_info.st_ino):
            fail("logical output directory identity changed")
        return {
            "trace": _read_regular_bytes(
                Path("logical-recovery-traces.jsonl"),
                "terminal logical recovery trace", directory_fd),
            "freeze": _read_regular_bytes(
                Path("logical-recovery-freeze.json"),
                "terminal logical recovery freeze", directory_fd),
            "result": _read_regular_bytes(
                Path("logical-recovery-results.jsonl"),
                "terminal logical recovery result", directory_fd),
        }
    except OSError as exc:
        fail("cannot reopen logical output directory: {}".format(exc))
    finally:
        if directory_fd >= 0:
            os.close(directory_fd)


def _parse_exact_jsonl(
        data: bytes, context: str) -> List[Mapping[str, Any]]:
    _require_bounded_jsonl_records(data, context)
    rows = []
    # JSONL is delimited by LF only.  bytes.splitlines() also treats CR, VT,
    # and several other control bytes as boundaries, which can fragment one
    # oversized physical record into apparently bounded pieces.
    for index, payload in enumerate(data[:-1].split(b"\n"), 1):
        line = payload + b"\n"
        try:
            rows.append(runner_api._parse_canonical_line(
                line, "{} line {}".format(context, index)))
        except runner_api.RunnerError as exc:
            fail(str(exc))
    return rows


def _require_bounded_jsonl_records(data: bytes, context: str) -> None:
    """Reject malformed framing and records above the recovery wire cap."""
    if not data or not data.endswith(b"\n"):
        fail("{} must be a nonempty newline-terminated JSONL stream".format(
            context))
    for index, payload in enumerate(data[:-1].split(b"\n"), 1):
        if len(payload) + 1 > MAX_RECOVERY_RECORD_BYTES:
            fail("{} line {} exceeds the 4096-byte record cap".format(
                context, index))


def _validate_bound_recovery_trace_bytes(
        contract: Mapping[str, Any], freeze: Mapping[str, Any], data: bytes,
        ) -> None:
    rows = _parse_exact_jsonl(data, "campaign trace manifest")
    cells = list(contract_api.iter_recovery_cells(contract, "development"))
    expected_count = contract["recovery"]["domains"]["development"][
        "expected_cells_per_arm"]
    if len(cells) != expected_count or len(rows) != expected_count:
        fail("campaign trace manifest has the wrong frozen cardinality")
    digest = contract_api._trace_manifest_hasher(
        contract, "recovery", "development")
    for ordinal, (row, cell) in enumerate(zip(rows, cells)):
        if (set(row) != contract_api.TRACE_FIELDS or
                type(row.get("ordinal")) is not int or
                row.get("ordinal") != ordinal or
                not _is_sha256(row.get("cell_sha256")) or
                not _is_sha256(row.get("trace_sha256")) or
                row.get("cell_sha256") != contract_api.sha256_json(cell)):
            fail("campaign trace manifest row {} differs from its frozen cell"
                 .format(ordinal))
        contract_api._hash_trace_manifest_row(digest, row)
    if digest.hexdigest() != freeze.get("trace_manifest_sha256"):
        fail("campaign trace bytes differ from the frozen trace hash")


def _write_snapshot(path: Path, data: bytes) -> None:
    try:
        with path.open("xb") as output:
            output.write(data)
            output.flush()
            os.fsync(output.fileno())
    except OSError as exc:
        fail("cannot write private campaign snapshot {}: {}".format(path, exc))


def _validate_completed_campaign_snapshots(
        contract: Mapping[str, Any], resolved: Path,
        directory_identity: Tuple[int, int],
        snapshots: Mapping[str, bytes],
        ) -> Mapping[str, Any]:
    """Semantically validate one exact private eight-artifact snapshot."""
    if (set(snapshots) != set(RECOVERY_COMPLETED_SOURCE_NAMES) or
            any(not isinstance(value, bytes) for value in snapshots.values())):
        fail("completed recovery snapshot is not the exact eight-artifact bundle")
    total_bytes = 0
    for key, name in RECOVERY_COMPLETED_SOURCE_NAMES.items():
        data = snapshots[key]
        if len(data) > RECOVERY_COMPLETED_ARTIFACT_BYTE_LIMITS[name]:
            fail("{} exceeds its completed-artifact byte cap".format(
                RECOVERY_COMPLETED_CONTEXTS[key]))
        total_bytes += len(data)
    if total_bytes > MAX_COMPLETED_CAMPAIGN_BYTES:
        fail("completed recovery campaign exceeds its aggregate cap")
    for key in ("trace", "native", "result"):
        _require_bounded_jsonl_records(
            snapshots[key], RECOVERY_COMPLETED_CONTEXTS[key])
    source_names = RECOVERY_COMPLETED_SOURCE_NAMES
    with tempfile.TemporaryDirectory(prefix="wh2-recovery-snapshot-") as raw:
        snapshot_root = Path(raw)
        paths = {
            key: snapshot_root / name
            for key, name in source_names.items()
        }
        for key, path in paths.items():
            _write_snapshot(path, snapshots[key])

        summary = _parse_canonical_object(
            snapshots["summary"], "campaign run summary")
        if set(summary) != CAMPAIGN_SUMMARY_FIELDS or \
                summary.get("schema") != CAMPAIGN_SUMMARY_SCHEMA or \
                summary.get("status") != "complete" or \
                summary.get("summary_sha256") != _self_hash(
                    summary, "summary_sha256"):
            fail("campaign run summary is incomplete or has an invalid self-hash")
        candidate_id = summary.get("candidate_id")
        candidate_arm, _ = _candidate(candidate_id)
        if summary.get("candidate_arm") != candidate_arm or \
                summary.get("output_dir") != str(resolved):
            fail("campaign run summary is bound to another candidate or directory")
        try:
            validated = native_api.validate_execution_receipt(
                contract, "recovery", "development", paths["freeze"],
                paths["trace"], paths["native"], paths["result"],
                paths["receipt"], verify_live_sampler=False,
                sampler_path=paths["sampler"])
            freeze = contract_api.load_freeze_manifest(
                contract, "development", paths["freeze"], "recovery")
        except (native_api.NativeEvidenceError,
                contract_api.ContractError) as exc:
            fail(str(exc))
        _validate_candidate_freeze(freeze, candidate_id)
        receipt = validated["execution_receipt"]
        if contract_api.freeze_manifest_sha256(freeze) != \
                receipt.get("freeze_manifest_sha256"):
            fail("reopened candidate freeze differs from its execution receipt")
        raw_basis, raw_schedule_sha256 = _raw_policy_from_freeze(freeze)
        witness = _parse_canonical_object(
            snapshots["witness"], "campaign timing proxy witness")
        _validate_timing_proxy_witness(
            witness, freeze["arms"][0]["binary_sha256"],
            freeze["source_git_commit"])
        witness_sha256 = contract_api.sha256_json(witness)
        if (witness_sha256 != freeze["timing_proxy_witness_sha256"] or
                witness_sha256 != summary["timing_proxy_witness_sha256"]):
            fail("campaign timing proxy witness differs from freeze or summary")
        for field in (
                "work_rank_summary_sha256",
                "work_rank_result_stream_sha256",
                "work_rank_domain_sha256"):
            if summary.get(field) != freeze.get(field):
                fail("campaign work/rank binding differs from its freeze")
        expected_summary = {
            "source_git_commit": freeze["source_git_commit"],
            "contract_sha256": contract_api.contract_sha256(contract),
            "domain_sha256": freeze["domain_sha256"],
            "trace_manifest_sha256": freeze["trace_manifest_sha256"],
            "worker_binary_sha256": freeze["arms"][0]["binary_sha256"],
            "recovery_records": receipt["record_count"],
            "recovery_freeze_sha256": receipt["freeze_manifest_sha256"],
            "recovery_result_sha256": receipt["result_stream_sha256"],
            "recovery_execution_receipt_sha256": receipt["receipt_sha256"],
            "worker_cpus": receipt["worker_cpus"],
            "controller_cpu": freeze.get("host_identity", {}).get(
                "controller_cpu"),
            "thermal_samples": receipt["thermal"]["sample_count"],
            "cpu_tctl_max_millic": receipt["thermal"]["cpu_tctl_max_millic"],
            "dimm_max_millic": receipt["thermal"]["dimm_max_millic"],
            "construction_seed_basis": raw_basis,
            "seed_schedule_sha256": raw_schedule_sha256,
            "timing_proxy_witness_sha256": witness_sha256,
            "work_rank_summary_sha256":
                freeze["work_rank_summary_sha256"],
            "work_rank_result_stream_sha256":
                freeze["work_rank_result_stream_sha256"],
            "work_rank_domain_sha256":
                freeze["work_rank_domain_sha256"],
        }
        if (type(expected_summary["controller_cpu"]) is not int or
                expected_summary["controller_cpu"] < 0 or
                any(contract_api.canonical_json(summary.get(field)) !=
                    contract_api.canonical_json(value)
                    for field, value in expected_summary.items()) or
                contract_api.canonical_json(
                    summary.get("recovery_records")) !=
                contract_api.canonical_json(RECOVERY_RECORDS)):
            fail("campaign run summary differs from its validated execution")
        rows = _parse_exact_jsonl(
            snapshots["result"], "campaign recovery ledger")
        if len(rows) != RECOVERY_RECORDS or \
                _hash_bytes(snapshots["result"]) != \
                receipt["result_stream_sha256"]:
            fail("campaign recovery ledger differs from its receipted result hash")
        raw_join = _native_raw_identity_join_binding(rows)
        if any(contract_api.canonical_json(summary.get(field)) !=
               contract_api.canonical_json(value)
               for field, value in raw_join.items()):
            fail("campaign raw identity join receipt differs from native rows")
        _validate_bound_recovery_trace_bytes(
            contract, freeze, snapshots["trace"])
        return {
            "directory": str(resolved),
            "directory_identity": directory_identity,
            "candidate_id": candidate_id,
            "candidate_arm": candidate_arm,
            "summary": summary,
            "freeze": freeze,
            "receipt": receipt,
            "rows": rows,
            "trace_bytes": snapshots["trace"],
            "timing_proxy_witness": witness,
        }


def _open_pinned_completed_campaign_bundle(
        source_names: Mapping[str, str], directory_fd: int,
        descriptors: Dict[str, int],
        ) -> Tuple[Dict[str, bytes],
                   Dict[str, runner_api.CompletedFingerprint]]:
    """Retain every recovery artifact inode under recovery-specific caps."""
    if (not isinstance(descriptors, dict) or descriptors or
            not source_names or
            not set(source_names).issubset(RECOVERY_COMPLETED_SOURCE_NAMES) or
            any(RECOVERY_COMPLETED_SOURCE_NAMES[key] != name
                for key, name in source_names.items())):
        fail("pinned completed recovery bundle is not an exact known subset")
    snapshots: Dict[str, bytes] = {}
    fingerprints: Dict[str, runner_api.CompletedFingerprint] = {}
    total_bytes = 0
    try:
        for key, name in source_names.items():
            limit = RECOVERY_COMPLETED_ARTIFACT_BYTE_LIMITS.get(name)
            if type(limit) is not int or limit <= 0:
                fail("completed recovery bundle names an unbounded artifact")
            try:
                data, fingerprint = \
                    runner_api._open_completed_regular_snapshot(
                        Path(name), RECOVERY_COMPLETED_CONTEXTS[key],
                        directory_fd, limit, descriptors, key)
            except runner_api.RunnerError as exc:
                fail(str(exc))
            snapshots[key] = data
            fingerprints[key] = fingerprint
            total_bytes += len(data)
            if total_bytes > MAX_COMPLETED_CAMPAIGN_BYTES:
                fail("completed recovery campaign exceeds its aggregate cap")
        return snapshots, fingerprints
    except BaseException:
        try:
            _close_recovery_descriptors(
                list(descriptors.values()),
                "pinned recovery bundle failure cleanup",
                sys.exc_info()[1])
        finally:
            descriptors.clear()
        raise


def _reread_pinned_completed_campaign_bundle(
        source_names: Mapping[str, str], directory_fd: int,
        descriptors: Mapping[str, int],
        expected_fingerprints:
            Mapping[str, runner_api.CompletedFingerprint],
        ) -> Tuple[Dict[str, bytes],
                   Dict[str, runner_api.CompletedFingerprint]]:
    """Reread retained inodes and prove all public names still select them."""
    if (set(descriptors) != set(source_names) or
            set(expected_fingerprints) != set(source_names)):
        fail("pinned completed recovery bundle is incomplete")
    snapshots: Dict[str, bytes] = {}
    fingerprints: Dict[str, runner_api.CompletedFingerprint] = {}
    total_bytes = 0
    for key, name in source_names.items():
        limit = RECOVERY_COMPLETED_ARTIFACT_BYTE_LIMITS[name]
        try:
            data, fingerprint = runner_api._read_completed_descriptor_bytes(
                descriptors[key], RECOVERY_COMPLETED_CONTEXTS[key], limit)
            named = os.stat(
                name, dir_fd=directory_fd, follow_symlinks=False)
        except runner_api.RunnerError as exc:
            fail(str(exc))
        except OSError as exc:
            fail("cannot terminally inspect {}: {}".format(
                RECOVERY_COMPLETED_CONTEXTS[key], exc))
        if (not stat.S_ISREG(named.st_mode) or
                fingerprint != expected_fingerprints[key] or
                runner_api._completed_fingerprint(named) != fingerprint):
            fail("{} changed during semantic validation".format(
                RECOVERY_COMPLETED_CONTEXTS[key]))
        snapshots[key] = data
        fingerprints[key] = fingerprint
        total_bytes += len(data)
        if total_bytes > MAX_COMPLETED_CAMPAIGN_BYTES:
            fail("completed recovery campaign exceeds its aggregate cap")
    return snapshots, fingerprints


def _require_pinned_completed_campaign_unchanged(
        source_names: Mapping[str, str], directory_fd: int,
        descriptors: Mapping[str, int],
        expected_fingerprints:
            Mapping[str, runner_api.CompletedFingerprint]) -> None:
    """Make the last cheap retained-inode/name check before publication."""
    if (set(descriptors) != set(source_names) or
            set(expected_fingerprints) != set(source_names)):
        fail("pinned completed recovery bundle is incomplete")
    for key, name in source_names.items():
        try:
            retained = os.fstat(descriptors[key])
            named = os.stat(
                name, dir_fd=directory_fd, follow_symlinks=False)
        except OSError as exc:
            fail("cannot terminally inspect {}: {}".format(
                RECOVERY_COMPLETED_CONTEXTS[key], exc))
        expected = expected_fingerprints[key]
        if (not stat.S_ISREG(retained.st_mode) or
                not stat.S_ISREG(named.st_mode) or
                runner_api._completed_fingerprint(retained) != expected or
                runner_api._completed_fingerprint(named) != expected):
            fail("{} changed before completion publication".format(
                RECOVERY_COMPLETED_CONTEXTS[key]))


def _read_completed_campaign_bundle(
        directory_fd: int,
        ) -> Tuple[Dict[str, bytes],
                   Dict[str, runner_api.CompletedFingerprint]]:
    """Read one full bundle while retaining all eight selected inodes."""
    descriptors: Dict[str, int] = {}
    try:
        snapshots, fingerprints = _open_pinned_completed_campaign_bundle(
            RECOVERY_COMPLETED_SOURCE_NAMES, directory_fd, descriptors)
        _require_pinned_completed_campaign_unchanged(
            RECOVERY_COMPLETED_SOURCE_NAMES, directory_fd, descriptors,
            fingerprints)
        return snapshots, fingerprints
    finally:
        _close_recovery_descriptors(
            list(descriptors.values()), "completed recovery bundle cleanup",
            sys.exc_info()[1])


def _open_completed_campaign_directory(
        campaign_dir: Path,
        expected_identity: Optional[Tuple[int, int]] = None,
        descriptor_holder: Optional[List[int]] = None,
        ) -> Tuple[int, Path, Tuple[int, int]]:
    """Open and pin a real campaign directory, rejecting path substitution."""
    if (descriptor_holder is not None and
            (not isinstance(descriptor_holder, list) or
             descriptor_holder != [-1])):
        fail("campaign directory descriptor holder is invalid")
    directory_fd = -1
    try:
        nofollow = getattr(os, "O_NOFOLLOW", 0)
        directory_flag = getattr(os, "O_DIRECTORY", 0)
        if nofollow == 0 or directory_flag == 0:
            fail("campaign directory cannot be opened fail-closed")
        directory_fd = os.open(
            str(campaign_dir), os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) |
            nofollow | directory_flag)
        if descriptor_holder is not None:
            descriptor_holder[0] = directory_fd
        opened_info = os.fstat(directory_fd)
        if not stat.S_ISDIR(opened_info.st_mode):
            fail("campaign path must be a real directory, not a symlink")
        identity = (opened_info.st_dev, opened_info.st_ino)
        if expected_identity is not None and identity != expected_identity:
            fail("campaign directory changed before terminal reread")
        resolved = campaign_dir.resolve(strict=True)
        resolved_info = os.stat(str(resolved), follow_symlinks=False)
        if (not stat.S_ISDIR(resolved_info.st_mode) or
                (resolved_info.st_dev, resolved_info.st_ino) != identity):
            fail("campaign directory identity changed while opening")
        return directory_fd, resolved, identity
    except OSError as exc:
        if directory_fd >= 0:
            if descriptor_holder is not None:
                descriptor_holder[0] = -1
            _close_recovery_descriptors(
                [directory_fd], "failed recovery directory-open cleanup",
                sys.exc_info()[1])
        fail("cannot open campaign directory {}: {}".format(campaign_dir, exc))
    except BaseException:
        if directory_fd >= 0:
            if descriptor_holder is not None:
                descriptor_holder[0] = -1
            _close_recovery_descriptors(
                [directory_fd], "interrupted recovery directory-open cleanup",
                sys.exc_info()[1])
        raise
    raise AssertionError("unreachable")


def load_completed_campaign(
        contract: Mapping[str, Any], campaign_dir: Path,
        ) -> Mapping[str, Any]:
    """Twice snapshot and revalidate one terminal four-arm campaign."""
    initial_fd = -1
    initial_fd_holder = [-1]
    try:
        opened_fd, resolved, directory_identity = \
            _open_completed_campaign_directory(
                campaign_dir, descriptor_holder=initial_fd_holder)
        initial_fd = initial_fd_holder[0]
        if opened_fd != initial_fd:
            fail("initial campaign directory ownership handoff failed")
        initial_snapshots, initial_fingerprints = \
            _read_completed_campaign_bundle(initial_fd)
        runner_api._verify_completed_directory_path(
            resolved, directory_identity)
    except runner_api.RunnerError as exc:
        fail(str(exc))
    finally:
        owned_initial_fd = initial_fd_holder[0]
        initial_fd_holder[0] = -1
        initial_fd = -1
        if owned_initial_fd >= 0:
            _close_recovery_descriptors(
                [owned_initial_fd], "initial recovery directory cleanup",
                sys.exc_info()[1])

    validated = _validate_completed_campaign_snapshots(
        contract, resolved, directory_identity, initial_snapshots)

    terminal_fd = -1
    terminal_fd_holder = [-1]
    try:
        opened_fd, terminal_resolved, terminal_identity = \
            _open_completed_campaign_directory(
                resolved, directory_identity, terminal_fd_holder)
        terminal_fd = terminal_fd_holder[0]
        if opened_fd != terminal_fd:
            fail("terminal campaign directory ownership handoff failed")
        terminal_snapshots, terminal_fingerprints = \
            _read_completed_campaign_bundle(terminal_fd)
        runner_api._verify_completed_directory_path(
            terminal_resolved, directory_identity)
    except runner_api.RunnerError as exc:
        fail(str(exc))
    finally:
        owned_terminal_fd = terminal_fd_holder[0]
        terminal_fd_holder[0] = -1
        terminal_fd = -1
        if owned_terminal_fd >= 0:
            _close_recovery_descriptors(
                [owned_terminal_fd], "terminal recovery directory cleanup",
                sys.exc_info()[1])
    if (terminal_identity != directory_identity or
            terminal_resolved != resolved or
            terminal_snapshots != initial_snapshots or
            terminal_fingerprints != initial_fingerprints):
        fail("completed recovery evidence changed during semantic validation")
    return validated


def _commit_completed_recovery_screen(
        contract: Mapping[str, Any], output_dir: Path,
        summary: Mapping[str, Any], freeze: Mapping[str, Any],
        expected_validation: Mapping[str, Any],
        expected_witness: Mapping[str, Any],
        finish_hard_wall: Optional[Callable[[], None]] = None) -> None:
    """Validate an unnamed summary, finish the wall, then hard-link once."""
    summary_path = output_dir / RECOVERY_COMPLETED_SOURCE_NAMES["summary"]
    summary_bytes = (
        contract_api.canonical_json(summary) + "\n").encode("utf-8")
    if len(summary_bytes) > RECOVERY_COMPLETED_ARTIFACT_BYTE_LIMITS[
            "run-summary.json"]:
        fail("recovery completion marker exceeds its byte cap")
    parent_fd = -1
    parent_fd_holder = [-1]
    summary_fd_holder = [-1]
    dependency_fds: Dict[str, int] = {}
    parent_identity: Optional[Tuple[int, int]] = None
    try:
        try:
            opened_parent_fd, parent_identity = \
                runner_api._open_completion_parent(
                    summary_path, parent_fd_holder)
            parent_fd = parent_fd_holder[0]
            if opened_parent_fd != parent_fd:
                fail("recovery completion parent ownership handoff failed")
            runner_api._require_completion_parent_descriptor(
                parent_fd, parent_identity)
            runner_api._verify_completed_directory_path(
                output_dir, parent_identity)
            runner_api._require_completion_marker_absent(
                summary_path, parent_fd)
            summary_fingerprint = runner_api._open_unnamed_completion_marker(
                parent_fd, summary_bytes, summary_fd_holder)
        except runner_api.RunnerError as exc:
            fail(str(exc))
        summary_fd = summary_fd_holder[0]

        dependency_snapshots, dependency_fingerprints = \
            _open_pinned_completed_campaign_bundle(
                RECOVERY_COMPLETED_DEPENDENCY_NAMES, parent_fd,
                dependency_fds)
        if len(summary_bytes) + sum(
                len(value) for value in dependency_snapshots.values()) > \
                MAX_COMPLETED_CAMPAIGN_BYTES:
            fail("completed recovery campaign exceeds its aggregate cap")
        prospective_snapshots = dict(dependency_snapshots)
        prospective_snapshots["summary"] = summary_bytes
        resolved = output_dir.resolve(strict=True)
        resolved_info = os.stat(str(resolved), follow_symlinks=False)
        if (not stat.S_ISDIR(resolved_info.st_mode) or
                (resolved_info.st_dev, resolved_info.st_ino) !=
                parent_identity):
            fail("recovery campaign directory changed before validation")
        validated = _validate_completed_campaign_snapshots(
            contract, resolved, parent_identity, prospective_snapshots)
        expected_receipt = expected_validation.get("execution_receipt") \
            if isinstance(expected_validation, Mapping) else None
        if (not isinstance(validated, Mapping) or
                validated.get("directory_identity") != parent_identity or
                validated.get("candidate_id") != summary.get("candidate_id") or
                validated.get("candidate_arm") != summary.get("candidate_arm") or
                contract_api.canonical_json(validated.get("summary")) !=
                contract_api.canonical_json(summary) or
                contract_api.canonical_json(validated.get("freeze")) !=
                contract_api.canonical_json(freeze) or
                contract_api.canonical_json(validated.get("receipt")) !=
                contract_api.canonical_json(expected_receipt) or
                contract_api.canonical_json(
                    validated.get("timing_proxy_witness")) !=
                contract_api.canonical_json(expected_witness)):
            fail("prospective recovery validation differs from terminal evidence")

        try:
            final_summary_bytes, final_summary_fingerprint = \
                runner_api._read_unnamed_completion_marker(
                    summary_fd,
                    RECOVERY_COMPLETED_ARTIFACT_BYTE_LIMITS[
                        "run-summary.json"])
        except runner_api.RunnerError as exc:
            fail(str(exc))
        final_dependencies, final_dependency_fingerprints = \
            _reread_pinned_completed_campaign_bundle(
                RECOVERY_COMPLETED_DEPENDENCY_NAMES, parent_fd,
                dependency_fds, dependency_fingerprints)
        if (final_summary_bytes != summary_bytes or
                final_summary_fingerprint != summary_fingerprint):
            fail("prospective recovery completion marker changed")
        if (final_dependencies != dependency_snapshots or
                final_dependency_fingerprints != dependency_fingerprints):
            fail("completed recovery evidence changed during validation")
        try:
            runner_api._require_completion_parent_descriptor(
                parent_fd, parent_identity)
            runner_api._verify_completed_directory_path(
                output_dir, parent_identity)
            runner_api._require_completion_marker_absent(
                summary_path, parent_fd)
        except runner_api.RunnerError as exc:
            fail(str(exc))

        if finish_hard_wall is not None:
            finish_hard_wall()
        try:
            post_finish_summary = os.fstat(summary_fd)
        except OSError as exc:
            fail("cannot recheck unnamed recovery completion marker: {}".
                 format(exc))
        if runner_api._completed_fingerprint(post_finish_summary) != \
                summary_fingerprint:
            fail("prospective recovery completion marker changed before publication")
        _require_pinned_completed_campaign_unchanged(
            RECOVERY_COMPLETED_DEPENDENCY_NAMES, parent_fd, dependency_fds,
            dependency_fingerprints)
        try:
            runner_api._require_completion_parent_descriptor(
                parent_fd, parent_identity)
            runner_api._verify_completed_directory_path(
                output_dir, parent_identity)
            runner_api._require_completion_marker_absent(
                summary_path, parent_fd)
            runner_api._link_unnamed_completion_marker(
                summary_path, summary_fd, parent_fd, parent_identity)
        except runner_api.RunnerError as exc:
            fail(str(exc))
    finally:
        _close_recovery_descriptors(
            [*dependency_fds.values(), summary_fd_holder[0],
             parent_fd_holder[0]],
            "recovery completion descriptor cleanup", sys.exc_info()[1])


def _combine_loaded_campaigns(
        contract: Mapping[str, Any], campaigns: Sequence[Mapping[str, Any]],
        ) -> Tuple[Mapping[str, Any], List[Mapping[str, Any]],
                   List[Mapping[str, Any]], bytes]:
    """Re-open the sole four-arm campaign as a logical evidence bundle."""
    if len(campaigns) != len(CANDIDATE_SPECS):
        fail("combination requires exactly one completed campaign directory")
    by_id: Dict[str, Mapping[str, Any]] = {}
    identities = set()
    for campaign in campaigns:
        candidate_id = campaign.get("candidate_id")
        _candidate(candidate_id)
        if candidate_id in by_id:
            fail("combination contains a duplicate candidate campaign")
        identity = campaign.get("directory_identity")
        if identity is not None and identity in identities:
            fail("combination reuses the same campaign directory")
        identities.add(identity)
        by_id[candidate_id] = campaign
    expected_ids = [value[0] for value in CANDIDATE_SPECS]
    if set(by_id) != set(expected_ids):
        fail("combination does not cover the closed candidate roster")
    ordered = [by_id[candidate_id] for candidate_id in expected_ids]
    first = ordered[0]

    common_fields = (
        "contract_sha256", "domain_sha256", "trace_manifest_sha256",
        "source_git_commit",
    )
    for campaign in ordered[1:]:
        for field in common_fields:
            if campaign["freeze"].get(field) != first["freeze"].get(field):
                fail("campaigns differ in their {}".format(field))
        if campaign.get("trace_bytes") != first.get("trace_bytes"):
            fail("campaign trace manifests are not byte-identical")
    expected_contract = contract_api.contract_sha256(contract)
    expected_domain = contract["recovery"]["domains"]["development"][
        "domain_sha256"]
    if first["freeze"].get("contract_sha256") != expected_contract or \
            first["freeze"].get("domain_sha256") != expected_domain:
        fail("campaigns do not bind the loaded recovery contract/domain")

    binaries = set()
    for campaign in ordered:
        if contract_api.freeze_manifest_sha256(campaign["freeze"]) != \
                campaign["receipt"].get("freeze_manifest_sha256"):
            fail("candidate campaign freeze differs from its receipt hash")
        freeze_arms = campaign["freeze"].get("arms")
        if not isinstance(freeze_arms, list):
            fail("campaign freeze arm records are malformed")
        binaries.update(arm.get("binary_sha256") for arm in freeze_arms)
    if len(binaries) != 1 or not _is_sha256(next(iter(binaries))):
        fail("campaigns were not produced by one identical worker binary")
    first_common = first["freeze"]["arms"][:3]
    for campaign in ordered[1:]:
        if campaign["freeze"]["arms"][:3] != first_common:
            fail("campaigns substitute a common control arm descriptor")
    cell_count = contract["recovery"]["domains"]["development"][
        "expected_cells_per_arm"]
    if cell_count * 4 != RECOVERY_RECORDS:
        fail("development recovery cardinality is not 360 cells per arm")
    for campaign in ordered:
        rows = campaign.get("rows")
        if not isinstance(rows, list) or len(rows) != RECOVERY_RECORDS:
            fail("candidate campaign does not contain exactly 1440 payload rows")
        if _hash_jsonl(rows) != campaign["receipt"].get(
                "result_stream_sha256"):
            fail("candidate campaign rows differ from their result-stream hash")
        expected_roster = campaign["freeze"]["arm_roster"]
        for cell in range(cell_count):
            for arm_index, arm in enumerate(expected_roster):
                if rows[cell * 4 + arm_index].get("arm") != arm:
                    fail("candidate ledger rows are not in frozen ordinal order")
    for campaign in ordered[1:]:
        for cell in range(cell_count):
            for common_index in (0, 1, 2):
                if campaign["rows"][cell * 4 + common_index] != \
                        first["rows"][cell * 4 + common_index]:
                    fail("campaign common payloads are not cell-identical")

    roster = [value[0] for value in CONTROL_ARMS] + [RAW_CONTROL_ARM[0]] + [
        value[1] for value in CANDIDATE_SPECS
    ]
    logical_arms = [dict(value) for value in first_common] + [
        dict(campaign["freeze"]["arms"][3]) for campaign in ordered
    ]
    commands: List[List[str]] = []
    cpu_affinity: Set[int] = set()
    input_freezes = []
    logical_rows: List[Mapping[str, Any]] = []
    for campaign in ordered:
        commands.extend(campaign["freeze"]["commands"])
        cpu_affinity.update(campaign["freeze"]["cpu_affinity"])
        input_freezes.append(contract_api.freeze_manifest_sha256(
            campaign["freeze"]))
    for cell in range(cell_count):
        logical_rows.extend((first["rows"][cell * 4],
                             first["rows"][cell * 4 + 1],
                             first["rows"][cell * 4 + 2]))
        logical_rows.extend(
            campaign["rows"][cell * 4 + 3] for campaign in ordered)
    if len(logical_rows) != LOGICAL_RECORDS:
        fail("logical recovery ledger is not exactly four arms by 360 cells")

    logical_freeze = {
        "schema": contract_api.RAW_FREEZE_SCHEMA,
        "contract_sha256": expected_contract,
        "evidence_kind": "recovery",
        "phase": "development",
        "domain_sha256": expected_domain,
        "source_git_commit": first["freeze"]["source_git_commit"],
        "arm_roster": roster,
        "arm_roster_sha256": contract_api.arm_roster_sha256(roster),
        "trace_manifest_sha256": first["freeze"]["trace_manifest_sha256"],
        "repair_training_trace_manifest_sha256": ZERO_SHA256,
        "commands": commands,
        "cpu_affinity": sorted(cpu_affinity),
        "host_identity": {
            "artifact_kind": "logical_recovery_combination",
            "input_freeze_manifest_sha256s": input_freezes,
        },
        "architecture_roles": dict(
            first["freeze"]["architecture_roles"]),
        "timing_proxy_witness_sha256":
            first["freeze"]["timing_proxy_witness_sha256"],
        "work_rank_summary_sha256":
            first["freeze"]["work_rank_summary_sha256"],
        "work_rank_result_stream_sha256":
            first["freeze"]["work_rank_result_stream_sha256"],
        "work_rank_domain_sha256":
            first["freeze"]["work_rank_domain_sha256"],
        "arms": logical_arms,
    }
    try:
        # Validate the synthesized structure before allowing publication.
        contract_api._exact_keys(
            logical_freeze, contract_api.RAW_V3_FREEZE_FIELDS,
            "logical recovery freeze")
    except contract_api.ContractError as exc:
        fail(str(exc))
    bindings = [{
        "candidate_id": campaign["candidate_id"],
        "candidate_arm": campaign["candidate_arm"],
        "run_summary_sha256": campaign["summary"]["summary_sha256"],
        "freeze_manifest_sha256": campaign["receipt"]
            ["freeze_manifest_sha256"],
        "result_stream_sha256": campaign["receipt"]["result_stream_sha256"],
        "execution_receipt_sha256": campaign["receipt"]["receipt_sha256"],
    } for campaign in ordered]
    return logical_freeze, logical_rows, bindings, first["trace_bytes"]


def _native_raw_identity_join_binding(
        rows: Sequence[Mapping[str, Any]]) -> Mapping[str, Any]:
    """Recompute the sidecar-join receipt projection from native rows only."""
    joined = []
    seen = set()
    for row in rows:
        if row.get("construction_seed_basis") != RAW_SEED_BASIS:
            continue
        if (row.get("base_seed_attempt") != row.get("trial") or
                row.get("construction_attempt") !=
                    row.get("precode_attempt") or
                row.get("precode_attempt") != row.get("packet_attempt")):
            fail("native raw row has inconsistent construction attempts")
        key = (
            row.get("arm"), row.get("K"), row.get("block_bytes"),
            row.get("trial"), row.get("schedule"), row.get("loss_ppm"),
            row.get("loss_seed"),
        )
        if key in seen:
            fail("native logical ledger repeats a raw construction identity")
        seen.add(key)
        try:
            joined.append({
                "K": row["K"],
                "arm": row["arm"],
                "block_bytes": row["block_bytes"],
                "cell_sha256": row["cell_sha256"],
                "dense_anchor_layout": row["dense_anchor_layout"],
                "loss_ppm": row["loss_ppm"],
                "loss_seed": row["loss_seed"],
                "realized_construction_sha256":
                    row["realized_construction_sha256"],
                "schedule": row["schedule"],
                "trace_sha256": row["trace_sha256"],
                "trial": row["trial"],
            })
        except KeyError as exc:
            fail("native raw identity row lacks {}".format(exc.args[0]))
    if len(joined) != RAW_IDENTITY_JOIN_COUNT:
        fail("native raw identity receipt has the wrong cardinality")
    return {
        "raw_identity_join_count": len(joined),
        "raw_identity_join_sha256": contract_api.sha256_json(joined),
    }


def _bind_work_rank_identities(
        logical_rows: Sequence[Mapping[str, Any]],
        work_screen: Mapping[str, Any], source_git_commit: str,
        ) -> Mapping[str, Any]:
    """Prove every native raw construction equals all four sidecar prefixes."""
    summary = work_screen.get("summary")
    work_rows = work_screen.get("rows")
    if (not isinstance(summary, Mapping) or
            not isinstance(work_rows, (list, tuple))):
        fail("work/rank artifact loader returned a malformed result")
    provenance = summary.get("source_provenance")
    if not isinstance(provenance, Mapping) or \
            provenance.get("source_git_commit") != source_git_commit:
        fail("native and work/rank evidence use different source commits")
    if (summary.get("construction_seed_basis") != RAW_SEED_BASIS or
            summary.get("seed_schedule_sha256") !=
            RAW_SEED_SCHEDULE_SHA256):
        fail("work/rank summary does not bind the frozen raw seed policy")

    grouped: Dict[Tuple[Any, ...], List[Mapping[str, Any]]] = {}
    for row in work_rows:
        if not isinstance(row, Mapping):
            fail("work/rank result contains a malformed row")
        key = (
            row.get("arm"), row.get("K"), row.get("block_bytes"),
            row.get("trial"), row.get("schedule"), row.get("loss_ppm"),
            row.get("loss_seed"),
        )
        grouped.setdefault(key, []).append(row)
    if len(work_rows) != work_api.EXPECTED_RECORDS or \
            len(grouped) != RAW_IDENTITY_JOIN_COUNT:
        fail("work/rank evidence has the wrong raw identity cardinality")

    common_fields = (
        "arm_descriptor_sha256", "construction_seed_basis",
        "seed_schedule_sha256", "precode_attempt", "packet_attempt",
        "effective_precode_seed", "effective_packet_seed", "staircase",
        "binary_dense_rows", "gf256_heavy_rows", "source_hits",
        "dense_anchor_layout", "dense_identity_corner", "heavy_family",
        "mix_count",
        "realized_construction_sha256",
    )
    seen = set()
    for native in logical_rows:
        if native.get("construction_seed_basis") != RAW_SEED_BASIS:
            continue
        if native.get("base_seed_attempt") != native.get("trial") or \
                native.get("construction_attempt") != \
                native.get("precode_attempt") or \
                native.get("precode_attempt") != native.get("packet_attempt"):
            fail("native raw row has inconsistent construction attempts")
        key = (
            native.get("arm"), native.get("K"), native.get("block_bytes"),
            native.get("trial"), native.get("schedule"),
            native.get("loss_ppm"), native.get("loss_seed"),
        )
        if key in seen:
            fail("native logical ledger repeats a raw construction identity")
        seen.add(key)
        matches = grouped.get(key)
        if matches is None or \
                {row.get("overhead") for row in matches} != {0, 1, 2, 4} or \
                len(matches) != 4:
            fail("native raw row lacks its exact four sidecar prefixes")
        for sidecar in matches:
            if any(native.get(field) != sidecar.get(field)
                   for field in common_fields):
                fail("native and work/rank effective constructions differ")
            if (native.get("cell_sha256") != sidecar.get("cell_sha256") or
                    native.get("trace_sha256") !=
                    sidecar.get("frozen_trace_sha256")):
                fail("native and work/rank frozen trace identities differ")
    if len(seen) != RAW_IDENTITY_JOIN_COUNT or seen != set(grouped):
        fail("native and work/rank raw identity domains differ")
    for field in ("summary_sha256", "result_stream_sha256",
                  "work_domain_sha256"):
        if not _is_sha256(summary.get(field)):
            fail("work/rank summary {} is not a SHA-256".format(field))
    return {
        "work_rank_summary_sha256": summary["summary_sha256"],
        "work_rank_result_stream_sha256": summary["result_stream_sha256"],
        "work_rank_domain_sha256": summary["work_domain_sha256"],
        **_native_raw_identity_join_binding(logical_rows),
    }


def combine_recovery_screens(args: argparse.Namespace) -> Mapping[str, Any]:
    contract = contract_api.load_contract(args.contract)
    campaigns = [load_completed_campaign(contract, path)
                 for path in args.campaign_dir]
    try:
        work_screen = work_api.load_completed_work_screen(
            contract, args.work_rank_dir)
    except work_api.WorkScreenError as exc:
        fail(str(exc))
    logical_freeze, logical_rows, bindings, trace_bytes = \
        _combine_loaded_campaigns(contract, campaigns)
    work_binding = _bind_work_rank_identities(
        logical_rows, work_screen, logical_freeze["source_git_commit"])
    for field in (
            "work_rank_summary_sha256",
            "work_rank_result_stream_sha256",
            "work_rank_domain_sha256"):
        if contract_api.canonical_json(work_binding.get(field)) != \
                contract_api.canonical_json(logical_freeze.get(field)):
            fail("loaded work/rank artifact differs from campaign freeze {}".
                 format(field))
    output_dir = runner_api._create_output_dir(args.output_dir)
    trace_path = output_dir / "logical-recovery-traces.jsonl"
    freeze_path = output_dir / "logical-recovery-freeze.json"
    result_path = output_dir / "logical-recovery-results.jsonl"
    runner_api._atomic_write_bytes(trace_path, trace_bytes)
    runner_api._atomic_write_object(freeze_path, logical_freeze)
    sink = runner_api.AtomicLineSink(
        result_path, maximum_bytes=LOGICAL_RESULT_STREAM_BYTE_CAP)
    try:
        for row in logical_rows:
            line = (contract_api.canonical_json(row) + "\n").encode("utf-8")
            if len(line) > MAX_RECOVERY_RECORD_BYTES:
                fail("logical recovery row exceeds the 4096-byte record cap")
            sink.write(line)
        sink.publish()
    finally:
        sink.abort()
    try:
        loaded_freeze = contract_api.load_freeze_manifest(
            contract, "development", freeze_path, "recovery")
        validator_summary = contract_api.validate_ledger(
            contract, "development", result_path, freeze_path, trace_path)
    except contract_api.ContractError as exc:
        fail(str(exc))
    terminal = _read_terminal_combination_artifacts(output_dir)
    terminal_freeze = _parse_canonical_object(
        terminal["freeze"], "terminal logical recovery freeze")
    terminal_freeze_sha256 = contract_api.freeze_manifest_sha256(
        terminal_freeze)
    if (contract_api.canonical_json(terminal_freeze) !=
            contract_api.canonical_json(logical_freeze) or
            terminal_freeze_sha256 !=
                contract_api.freeze_manifest_sha256(loaded_freeze) or
            terminal_freeze_sha256 !=
                validator_summary.get("freeze_manifest_sha256")):
        fail("terminal logical freeze differs from validated evidence")
    if terminal["trace"] != trace_bytes:
        fail("terminal logical trace differs from validated evidence")
    _validate_bound_recovery_trace_bytes(
        contract, terminal_freeze, terminal["trace"])
    terminal_rows = _parse_exact_jsonl(
        terminal["result"], "terminal logical recovery result")
    result_hash = _hash_jsonl(logical_rows)
    if (contract_api.canonical_json(terminal_rows) !=
            contract_api.canonical_json(logical_rows) or
            _hash_bytes(terminal["result"]) != result_hash):
        fail("terminal logical result differs from validated evidence")
    raw_basis, raw_schedule_sha256 = _raw_policy_from_freeze(terminal_freeze)
    summary: Dict[str, Any] = {
        "schema": COMBINATION_SUMMARY_SCHEMA,
        "status": "complete",
        "artifact_kind": "logical_recovery_combination",
        "is_execution_receipt": False,
        "phase": "development",
        "contract_sha256": contract_api.contract_sha256(contract),
        "domain_sha256": logical_freeze["domain_sha256"],
        "trace_manifest_sha256": logical_freeze["trace_manifest_sha256"],
        "trace_file_sha256": _hash_bytes(terminal["trace"]),
        "source_git_commit": logical_freeze["source_git_commit"],
        "worker_binary_sha256": logical_freeze["arms"][0]["binary_sha256"],
        "candidate_roster": [value[0] for value in CANDIDATE_SPECS],
        "arm_roster": logical_freeze["arm_roster"],
        "campaigns": bindings,
        "construction_seed_basis": raw_basis,
        "seed_schedule_sha256": raw_schedule_sha256,
        "logical_freeze_manifest_sha256":
            terminal_freeze_sha256,
        "logical_result_stream_sha256": result_hash,
        "validator_summary": validator_summary,
        "validator_summary_sha256":
            contract_api.sha256_json(validator_summary),
        **work_binding,
    }
    summary["combination_sha256"] = contract_api.sha256_json(summary)
    if set(summary) != COMBINATION_SUMMARY_FIELDS or \
            any(set(value) != CAMPAIGN_BINDING_FIELDS
                for value in summary["campaigns"]):
        fail("internal logical-combination summary schema mismatch")
    runner_api._atomic_write_object(
        output_dir / "logical-recovery-summary.json", summary)
    return summary


def _add_common_run_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--contract", type=Path, default=contract_api.DEFAULT_CONTRACT)
    parser.add_argument(
        "--worker", type=Path,
        default=Path("build/wirehair_wh2_contract_worker"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--candidate", choices=tuple(CANDIDATE_BY_ID),
                        required=True)
    parser.add_argument(
        "--work-rank-dir", type=Path, required=True,
        help="completed raw-v3 D12/Two07 work/rank artifact directory")
    parser.add_argument("--sampler-pid", type=int, required=True)
    parser.add_argument("--sampler-cpu", type=int, required=True)
    parser.add_argument("--sampler-script", type=Path, required=True)
    parser.add_argument("--sampler-csv", type=Path, required=True)
    parser.add_argument(
        "--cpus", help="strictly increasing comma-separated logical CPUs")
    parser.add_argument(
        "--controller-cpu", type=int,
        help="dedicated non-worker, non-sampler physical core")
    parser.add_argument(
        "--deadline-seconds", type=float,
        default=runner_api.MAX_WALL_SECONDS)


def main(argv: Sequence[str] = ()) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser(
        "run", help="run one exact four-arm native recovery campaign")
    _add_common_run_arguments(run)
    combine = commands.add_parser(
        "combine", help="combine the sole completed Two07 recovery campaign")
    combine.add_argument(
        "--contract", type=Path, default=contract_api.DEFAULT_CONTRACT)
    combine.add_argument(
        "--campaign-dir", type=Path, action="append", required=True,
        help="supply exactly once for the sole closed Two07 candidate")
    combine.add_argument(
        "--work-rank-dir", type=Path, required=True,
        help="completed raw-v3 D12/Two07 work/rank artifact directory")
    combine.add_argument("--output-dir", type=Path, required=True)
    combine.add_argument(
        "--deadline-seconds", type=float,
        default=runner_api.MAX_WALL_SECONDS)
    args = parser.parse_args(argv if argv else None)
    try:
        if args.command == "run":
            if (not math.isfinite(args.deadline_seconds) or
                    not 0.0 < args.deadline_seconds <=
                    runner_api.MAX_WALL_SECONDS):
                fail("--deadline-seconds must be in (0,7200]")
            with runner_api._hard_wall(
                    args.deadline_seconds) as finish_hard_wall:
                summary = run_recovery_screen(args, finish_hard_wall)
        else:
            if len(args.campaign_dir) != len(CANDIDATE_SPECS):
                fail("--campaign-dir must be supplied exactly one time")
            if (not math.isfinite(args.deadline_seconds) or
                    not 0.0 < args.deadline_seconds <=
                    runner_api.MAX_WALL_SECONDS):
                fail("--deadline-seconds must be in (0,7200]")
            with runner_api._hard_wall(args.deadline_seconds):
                summary = combine_recovery_screens(args)
    except (RecoveryRunnerError, runner_api.RunnerError,
            native_api.NativeEvidenceError, contract_api.ContractError,
            OSError, UnicodeError) as exc:
        print("wh2 native recovery screen: {}".format(exc), file=sys.stderr)
        return 1
    print(contract_api.canonical_json(summary))
    return 0


if __name__ == "__main__":
    sys.exit(main())
