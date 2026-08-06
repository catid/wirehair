#!/usr/bin/env python3
"""Tests for fail-closed native WH2 short-screen assembly."""

from __future__ import annotations

import copy
import csv
import errno
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Any, Dict, List, Mapping, NamedTuple, Optional, Sequence, Tuple
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_benchmark_contract as contract_api
import wh2_native_short_screen as subject


ARMS = ("wirehair2_head", "wirehair1", "candidate")
RAW_RECOVERY_ARMS = (
    "wirehair2_head", "wirehair1",
    "wirehair2_raw_d12_h12_periodic",
    "wirehair2_dense_two07_basis_v1",
)
RAW_RECOVERY_DESCRIPTORS = {
    "wirehair2_head":
        "4cafe27a8fb388ca9a4249b2c279b1406e7a0a86bcf14e98246988c7c503fa7a",
    "wirehair1":
        "d5a24d404e69efeb439907cd8271eba98d6af86b58efe159a820fb7aea08883d",
    "wirehair2_raw_d12_h12_periodic":
        "739092a7824449e6168f08b46661dfbe8ad5495ea4166b36073c79cd3bacdd11",
    "wirehair2_dense_two07_basis_v1":
        "9527f200ad38c7eec6502b2f768fdd67b92787fb227eed3d7616274ffc2df388",
}
WORKER_CPUS = tuple(range(8))
QUALIFICATION_CPUS = tuple(range(9))
WORKER_BINARY = hashlib.sha256(b"one native worker").hexdigest()
SOURCE_COMMIT = "1" * 40


class BuiltEvidence(NamedTuple):
    freeze: Path
    traces: Path
    native: Path
    records: List[Dict[str, Any]]
    qualification: Optional[contract_api.TimingQualification]
    qualification_native: Optional[Path]
    qualification_map: Optional[Path]
    qualification_audit: Optional[Path]
    qualification_map_sha256: Optional[str]
    qualification_sampler: Optional[Path]
    sampler: Optional[Path]
    qualification_execution_receipt: Optional[Path]
    qualification_execution_receipt_sha256: Optional[str]


class BuiltQualification(NamedTuple):
    value: contract_api.TimingQualification
    native: Path
    audit: Path
    map: Path
    map_sha256: str
    records: List[Dict[str, Any]]
    metadata: Mapping[str, Any]
    execution_receipt: Path
    execution_receipt_sha256: str


def timing_evidence_kwargs(value: BuiltEvidence) -> Dict[str, Any]:
    if value.qualification is None:
        return {}
    freeze = json.loads(value.freeze.read_text(encoding="utf-8"))
    worker_cpus = sorted({row["cpu"] for row in value.records})
    controller_cpu = freeze["host_identity"]["controller_cpu"]
    protected_cpus = sorted(worker_cpus + [controller_cpu])
    first = min(row["started_monotonic_ns"] for row in value.records)
    last = max(row["finished_monotonic_ns"] for row in value.records)
    if first < 3 or len(value.records) % len(worker_cpus) != 0:
        raise AssertionError("timing isolation fixture has invalid waves")
    checks = [
        {"ordinal": 0, "stage": "before-worker-spawn",
         "monotonic_ns": first - 3, "foreign_task_count": 0},
        {"ordinal": 1, "stage": "before-timing-window",
         "monotonic_ns": first - 2, "foreign_task_count": 0},
    ]
    chronological = sorted(
        value.records, key=lambda row: row["started_monotonic_ns"])
    wave_rows = [
        chronological[offset:offset + len(worker_cpus)]
        for offset in range(0, len(chronological), len(worker_cpus))
    ]
    for wave_ordinal, rows in enumerate(wave_rows):
        wave_start = min(row["started_monotonic_ns"] for row in rows)
        wave_end = max(row["finished_monotonic_ns"] for row in rows)
        checks.extend((
            {"ordinal": len(checks),
             "stage": "before-timing-wave-{}".format(wave_ordinal),
             "monotonic_ns": wave_start - 1, "foreign_task_count": 0},
            {"ordinal": len(checks) + 1,
             "stage": "after-timing-wave-{}".format(wave_ordinal),
             "monotonic_ns": wave_end + 1, "foreign_task_count": 0},
        ))
    checks.append({
        "ordinal": len(checks), "stage": "after-timing-interval",
        "monotonic_ns": last + 2, "foreign_task_count": 0,
    })
    isolation = {
        "schema": subject.SIBLING_ISOLATION_SCHEMA,
        "policy": "exact-affinity-wave-boundary-v1",
        "worker_cpus": worker_cpus,
        "controller_cpu": controller_cpu,
        "protected_cpus": protected_cpus,
        "sibling_cpus": subject.timing_sibling_cpus(protected_cpus),
        "paused_processes": [],
        "checks": checks,
        "check_count": len(checks),
        "checks_sha256": contract_api.sha256_json(checks),
        "first_check_monotonic_ns": checks[0]["monotonic_ns"],
        "last_check_monotonic_ns": checks[-1]["monotonic_ns"],
        "terminal_status": "clean",
    }
    return {
        "timing_qualification_map_path": value.qualification_map,
        "timing_qualification_audit_path": value.qualification_audit,
        "timing_qualification_map_sha256":
            value.qualification_map_sha256,
        "timing_qualification_native_path": value.qualification_native,
        "timing_qualification_sampler_path": value.qualification_sampler,
        "timing_qualification_execution_receipt_path":
            value.qualification_execution_receipt,
        "timing_qualification_execution_receipt_sha256":
            value.qualification_execution_receipt_sha256,
        "sibling_isolation": isolation,
    }


def digest(label: str) -> str:
    return hashlib.sha256(label.encode("ascii")).hexdigest()


def canonical_file(path: Path, value: Mapping[str, Any]) -> None:
    path.write_bytes(
        (contract_api.canonical_json(value) + "\n").encode("utf-8"))


def canonical_jsonl(path: Path, values: Sequence[Mapping[str, Any]]) -> None:
    path.write_bytes("".join(
        contract_api.canonical_json(value) + "\n" for value in values
    ).encode("utf-8"))


def arm_codec(arm: str) -> str:
    if arm == "wirehair2_head":
        return "wirehair2_certified"
    if arm == "wirehair1":
        return "wirehair1"
    return "wirehair2_experiment"


def descriptor(arm: str) -> str:
    return digest("descriptor:" + arm)


def trace_records(
        contract: Mapping[str, Any], kind: str, phase: str,
        qualification: Optional[contract_api.TimingQualification] = None,
        ) -> Tuple[List[Dict[str, Any]], List[str]]:
    cells = contract_api.iter_recovery_cells(contract, phase) if \
        kind == "recovery" else contract_api.iter_timing_cells(
            contract, phase, qualification)
    rows = []
    hashes = []
    for ordinal, cell in enumerate(cells):
        trace = digest("{}:{}:{}".format(
            kind, ordinal, contract_api.sha256_json(cell)))
        hashes.append(trace)
        rows.append({
            "schema": subject.TRACE_RECORD_SCHEMA,
            "ordinal": ordinal,
            "cell_sha256": contract_api.sha256_json(cell),
            "trace_sha256": trace,
            "packet_count": cell["K"] + (
                4 if kind == "recovery" else
                contract["timing"]["receive_overhead_cap"]),
            "candidate_count": cell["K"] + (
                4 if kind == "recovery" else
                contract["timing"]["receive_overhead_cap"]),
        })
    return rows, hashes


def timing_freeze_manifest(
        contract: Mapping[str, Any], phase: str, trace_sha256: str,
        qualification: contract_api.TimingQualification,
        ) -> Dict[str, Any]:
    records = []
    for arm in ARMS:
        codec = arm_codec(arm)
        records.append({
            "arm": arm,
            "codec": codec,
            "binary_sha256": WORKER_BINARY,
            "arm_descriptor_sha256": descriptor(arm),
            "construction_policy":
                "not_applicable" if codec == "wirehair1" else "raw_base",
            "repair_map_sha256": "0" * 64,
        })
    return {
        "schema": contract_api.FREEZE_SCHEMA,
        "contract_sha256": contract_api.contract_sha256(contract),
        "evidence_kind": "timing",
        "phase": phase,
        "domain_sha256": qualification.qualified_domain_sha256,
        "source_git_commit": SOURCE_COMMIT,
        "arm_roster": list(ARMS),
        "arm_roster_sha256": contract_api.arm_roster_sha256(ARMS),
        "trace_manifest_sha256": trace_sha256,
        "repair_training_trace_manifest_sha256": "0" * 64,
        "commands": [["wirehair_v2_contract_worker", "timing", phase]],
        "cpu_affinity": list(WORKER_CPUS),
        "host_identity": {
            "name": "native-fixture", "controller_cpu": 8,
        },
        "arms": records,
    }


def qualification_controls() -> List[Dict[str, Any]]:
    return [{
        "arm": arm,
        "scope": scope,
        "binary_sha256": WORKER_BINARY,
        "arm_descriptor_sha256": descriptor(arm),
        "construction_policy":
            "not_applicable" if arm == "wirehair1" else "raw_base",
        "repair_map_sha256": "0" * 64,
    } for arm, scope in (
        ("wirehair2_head", "decoder_solve"),
        ("wirehair1", "receive_to_success"),
    )]


def qualification_records(
        contract: Mapping[str, Any],
        retry_offsets: Optional[Sequence[int]] = None,
        cpus: Sequence[int] = QUALIFICATION_CPUS,
        ) -> List[Dict[str, Any]]:
    base_cells = list(contract_api.iter_timing_base_cells(
        contract, "development"))
    selected_offsets = [0] * len(base_cells) if retry_offsets is None else \
        list(retry_offsets)
    if len(selected_offsets) != len(base_cells):
        raise AssertionError("qualification fixture offset cardinality")
    rows: List[Dict[str, Any]] = []
    for base_ordinal, (base_cell, selected) in enumerate(zip(
            base_cells, selected_offsets)):
        if type(selected) is not int or not 0 <= selected < 256:
            raise AssertionError("qualification fixture retry range")
        base_sha256 = contract_api.sha256_json(base_cell)
        message_sha256 = subject._timing_qualification_message_sha256(
            base_cell)
        for retry_offset in range(selected + 1):
            terminal = retry_offset == selected
            qualified_cell = dict(base_cell)
            qualified_cell.update({
                "base_cell_sha256": base_sha256,
                "loss_retry_offset": retry_offset,
                "loss_seed": contract_api._qualified_timing_loss_seed(
                    base_cell["base_loss_seed"], retry_offset),
            })
            cell_sha256 = contract_api.sha256_json(qualified_cell)
            payload = {
                "ordinal": base_ordinal,
                "base_cell_sha256": base_sha256,
                "loss_retry_offset": retry_offset,
                "loss_seed": qualified_cell["loss_seed"],
                "cell_sha256": cell_sha256,
                "trace_sha256": digest(
                    "qualification-trace:{}:{}".format(
                        base_ordinal, retry_offset)),
                "packet_count": base_cell["K"] +
                    contract["timing"]["receive_overhead_cap"],
                "candidate_count": base_cell["K"] +
                    contract["timing"]["receive_overhead_cap"],
                "wirehair2_head_outcome":
                    "success" if terminal else "need_more_at_bound",
                "wirehair2_head_decoded_extra": 4 if terminal else None,
                "wirehair1_outcome": "success",
                "wirehair1_decoded_extra": 0,
            }
            outer_ordinal = base_ordinal * 256 + retry_offset
            cpu = cpus[outer_ordinal % len(cpus)]
            start = 1200000000 + outer_ordinal * 2
            rows.append({
                "schema": subject.TIMING_QUALIFICATION_RECORD_SCHEMA,
                "ordinal": outer_ordinal,
                "cpu": cpu,
                "worker_pid": 2000 + cpu,
                "worker_process_start_ticks": 1,
                "started_monotonic_ns": start,
                "finished_monotonic_ns": start + 1,
                "worker_binary_sha256": WORKER_BINARY,
                "message_sha256": message_sha256,
                "work_sha256":
                    subject._expected_timing_qualification_work_sha256(
                        outer_ordinal, cell_sha256),
                "payload": payload,
            })
    return rows


def build_qualification(
        root: Path, contract: Mapping[str, Any],
        retry_offsets: Optional[Sequence[int]] = None,
        records: Optional[Sequence[Mapping[str, Any]]] = None,
        ) -> BuiltQualification:
    qualification_native = root / "timing-qualification-native.jsonl"
    audit_path = root / "timing-qualification-audit.jsonl"
    map_path = root / "timing-qualification-map.json"
    execution_receipt_path = root / \
        "timing-qualification-execution.json"
    native_records = qualification_records(
        contract, retry_offsets) if records is None else [
            copy.deepcopy(value) for value in records]
    # Qualification workers may finish in any order.  The assembler must
    # canonicalize the stream by base ordinal and retry offset.
    canonical_jsonl(qualification_native, list(reversed(native_records)))
    value, metadata, execution_receipt_sha256 = \
        subject.assemble_timing_qualification(
        contract, "development", qualification_native, audit_path, map_path,
        execution_receipt_path, SOURCE_COMMIT, qualification_controls(),
        QUALIFICATION_CPUS,
        verify_live_workers=False)
    return BuiltQualification(
        value, qualification_native, audit_path, map_path,
        value.map_sha256, native_records, metadata,
        execution_receipt_path, execution_receipt_sha256)


def realized(arm: str, K: int, block_bytes: int, attempt: int) -> str:
    return contract_api.generic_realized_construction_sha256(
        arm_codec(arm), descriptor(arm), K, block_bytes, attempt)


def envelope(ordinal: int, message_ordinal: int,
             payload: Mapping[str, Any], schema: str,
             cpu: Optional[int] = None,
             started_monotonic_ns: Optional[int] = None,
             finished_monotonic_ns: Optional[int] = None) -> Dict[str, Any]:
    kind = "recovery" if schema in (
        subject.RECOVERY_RECORD_SCHEMA,
        subject.LEGACY_RAW_RECOVERY_RECORD_SCHEMA,
        subject.RAW_RECOVERY_RECORD_SCHEMA,
    ) else "timing"
    selected_cpu = WORKER_CPUS[ordinal % len(WORKER_CPUS)] \
        if cpu is None else cpu
    start = 2000000000 + ordinal * 2 \
        if started_monotonic_ns is None else started_monotonic_ns
    finish = start + 1 \
        if finished_monotonic_ns is None else finished_monotonic_ns
    return {
        "schema": schema,
        "ordinal": ordinal,
        "cpu": selected_cpu,
        "worker_pid": 1000 + selected_cpu,
        "worker_process_start_ticks": 1,
        "started_monotonic_ns": start,
        "finished_monotonic_ns": finish,
        "worker_binary_sha256": WORKER_BINARY,
        "message_sha256": digest("message:{}".format(message_ordinal)),
        "work_sha256": subject._expected_work_sha256(
            kind, "development", ordinal, payload),
        "payload": dict(payload),
    }


def raw_v3_recovery_fixture(
        root: Path, contract: Mapping[str, Any]) -> BuiltEvidence:
    """Build one exact mixed control/raw native recovery execution."""
    native_traces, trace_hashes = trace_records(
        contract, "recovery", "development")
    native_trace_path = root / "raw-v3-native-traces.jsonl"
    trace_path = root / "raw-v3-traces.jsonl"
    canonical_jsonl(native_trace_path, list(reversed(native_traces)))
    trace_sha256 = subject.assemble_trace_manifest(
        contract, "recovery", "development", native_trace_path, trace_path)

    arms = []
    for arm in RAW_RECOVERY_ARMS:
        codec = arm_codec(arm)
        raw = arm in RAW_RECOVERY_ARMS[2:]
        arms.append({
            "arm": arm,
            "codec": codec,
            "binary_sha256": WORKER_BINARY,
            "arm_descriptor_sha256": RAW_RECOVERY_DESCRIPTORS[arm],
            "construction_policy":
                "not_applicable" if codec == "wirehair1" else "raw_base",
            "repair_map_sha256": "0" * 64,
            "construction_seed_basis":
                contract_api.RAW_CONSTRUCTION_SEED_BASIS if raw else
                (contract_api.NOT_APPLICABLE_CONSTRUCTION_SEED_BASIS
                 if codec == "wirehair1" else
                 contract_api.PRODUCTION_CONSTRUCTION_SEED_BASIS),
            "seed_schedule_sha256":
                contract_api.RAW_SEED_SCHEDULE_SHA256 if raw else "0" * 64,
            "dense_anchor_layout":
                "two07" if arm == RAW_RECOVERY_ARMS[3] else
                ("not-applicable" if codec == "wirehair1" else "disabled"),
        })
    freeze = {
        "schema": contract_api.RAW_FREEZE_SCHEMA,
        "contract_sha256": contract_api.contract_sha256(contract),
        "evidence_kind": "recovery",
        "phase": "development",
        "domain_sha256": contract["recovery"]["domains"]["development"]
            ["domain_sha256"],
        "source_git_commit": SOURCE_COMMIT,
        "arm_roster": list(RAW_RECOVERY_ARMS),
        "arm_roster_sha256": contract_api.arm_roster_sha256(
            RAW_RECOVERY_ARMS),
        "trace_manifest_sha256": trace_sha256,
        "repair_training_trace_manifest_sha256": "0" * 64,
        "commands": [["wirehair_wh2_contract_worker", "raw-v3-recovery"]],
        "cpu_affinity": list(WORKER_CPUS),
        "host_identity": {
            "name": "raw-v3-native-fixture", "controller_cpu": 8,
        },
        "architecture_roles": copy.deepcopy(
            contract_api.EXPECTED_RAW_ARCHITECTURE_ROLES),
        "timing_proxy_witness_sha256": "1" * 64,
        "work_rank_summary_sha256": "2" * 64,
        "work_rank_result_stream_sha256": "3" * 64,
        "work_rank_domain_sha256": "4" * 64,
        "arms": arms,
    }
    freeze_path = root / "raw-v3-freeze.json"
    canonical_file(freeze_path, freeze)

    records = []
    for cell_ordinal, cell in enumerate(contract_api.iter_recovery_cells(
            contract, "development")):
        for arm_index, arm in enumerate(RAW_RECOVERY_ARMS):
            codec = arm_codec(arm)
            attempt = 0 if codec == "wirehair1" else \
                cell["base_seed_attempt"]
            payload = {
                "arm": arm,
                **cell,
                "outcome": "success",
                "decoded_extra": 0,
                "cell_sha256": contract_api.sha256_json(cell),
                "trace_sha256": trace_hashes[cell_ordinal],
                "binary_sha256": WORKER_BINARY,
                "arm_descriptor_sha256": RAW_RECOVERY_DESCRIPTORS[arm],
                "construction_attempt": attempt,
                "repair_map_sha256": "0" * 64,
            }
            schema = subject.RECOVERY_RECORD_SCHEMA
            if arm in RAW_RECOVERY_ARMS[2:]:
                raw_fields = {
                    "construction_seed_basis":
                        contract_api.RAW_CONSTRUCTION_SEED_BASIS,
                    "seed_schedule_sha256":
                        contract_api.RAW_SEED_SCHEDULE_SHA256,
                    "precode_attempt": attempt,
                    "packet_attempt": attempt,
                    "effective_precode_seed":
                        contract_api._effective_raw_precode_seed(attempt),
                    "effective_packet_seed":
                        contract_api._effective_raw_packet_seed(attempt),
                    "staircase": contract_api._raw_staircase_for_K(cell["K"]),
                    "binary_dense_rows": 12,
                    "gf256_heavy_rows": 12,
                    "source_hits": 3 if cell["K"] >= 10000 else 2,
                    "dense_anchor_layout":
                        "two07" if arm == RAW_RECOVERY_ARMS[3] else
                        "disabled",
                    "dense_identity_corner": False,
                    "heavy_family": "periodic-cauchy",
                    "mix_count": 3,
                }
                payload.update(raw_fields)
                payload["realized_construction_sha256"] = \
                    contract_api.raw_realized_construction_sha256(
                        codec, arm, RAW_RECOVERY_DESCRIPTORS[arm],
                        cell["K"], cell["block_bytes"], raw_fields)
                schema = subject.RAW_RECOVERY_RECORD_SCHEMA
            else:
                payload["realized_construction_sha256"] = \
                    contract_api.generic_realized_construction_sha256(
                        codec, RAW_RECOVERY_DESCRIPTORS[arm], cell["K"],
                        cell["block_bytes"], attempt)
            ordinal = cell_ordinal * len(RAW_RECOVERY_ARMS) + arm_index
            records.append(envelope(
                ordinal, cell_ordinal, payload, schema))
    native_path = root / "raw-v3-native-results.jsonl"
    canonical_jsonl(native_path, list(reversed(records)))
    return BuiltEvidence(
        freeze_path, trace_path, native_path, records,
        None, None, None, None, None, None, None, None, None)


def timing_records(
        contract: Mapping[str, Any], hashes: Sequence[str],
        qualification: contract_api.TimingQualification,
        ) -> List[Dict[str, Any]]:
    rows = []
    panels = contract_api.timing_panels(contract, ARMS)
    cells = list(contract_api.iter_timing_cells(
        contract, "development", qualification))
    cell_fields = contract["timing"]["cell_key"]
    stable_fields = [
        field for field in cell_fields
        if field not in (
            "replicate", "base_loss_seed", "base_cell_sha256",
            "loss_retry_offset", "loss_seed")
    ]
    identities = [contract_api.canonical_json({
        field: cell[field] for field in stable_fields
    }) for cell in cells]
    identity_order = []
    for cell, identity in zip(cells, identities):
        if cell["replicate"] == 0:
            identity_order.append(identity)
    identity_indexes = {
        identity: index for index, identity in enumerate(identity_order)
    }
    repetitions = contract["timing"]["domains"]["development"][
        "paired_repetitions"]
    jobs_per_wave = contract["timing"]["execution_geometry"][
        "jobs_per_wave"]
    independent_rounds = contract["timing"]["domains"]["development"][
        "independent_rounds"]
    if independent_rounds * jobs_per_wave != repetitions:
        raise AssertionError("fixture timing rounds do not partition repetitions")
    cohort_count = len(identity_order) * len(panels)
    for cell_ordinal, cell in enumerate(cells):
        for panel_index, panel in enumerate(panels):
            order = contract_api.timing_order(panel, cell["replicate"])
            payload: Dict[str, Any] = {
                **cell,
                **panel,
                "order": order,
                "left_outcome": "success",
                "right_outcome": "success",
                "left_decoded_extra": None if panel["scope"] ==
                    "encoder_init_plus_first_K_symbols" else 4,
                "right_decoded_extra": None if panel["scope"] ==
                    "encoder_init_plus_first_K_symbols" else 4,
                "elapsed_ns": [100000] * 8,
                "cell_sha256": contract_api.sha256_json(cell),
                "trace_sha256": hashes[cell_ordinal],
            }
            for side in ("left", "right"):
                arm = panel[side + "_arm"]
                codec = arm_codec(arm)
                attempt = 0 if codec == "wirehair1" else \
                    cell["base_seed_attempt"]
                payload[side + "_binary_sha256"] = WORKER_BINARY
                payload[side + "_arm_descriptor_sha256"] = descriptor(arm)
                payload[side + "_construction_attempt"] = attempt
                payload[side + "_realized_construction_sha256"] = realized(
                    arm, cell["K"], cell["block_bytes"], attempt)
                payload[side + "_repair_map_sha256"] = "0" * 64
            ordinal = cell_ordinal * len(panels) + panel_index
            cohort_index = identity_indexes[identities[cell_ordinal]] * \
                len(panels) + panel_index
            round_index = cell["replicate"] // jobs_per_wave
            global_wave = round_index * cohort_count + cohort_index
            position = cell["replicate"] % jobs_per_wave
            cpu = WORKER_CPUS[contract_api.timing_worker_slot(
                contract, "development", ARMS, cohort_index,
                cell["replicate"])]
            start = 3000000000 + global_wave * 1000000 + position * 10
            rows.append(envelope(
                ordinal, cell_ordinal, payload, subject.TIMING_RECORD_SCHEMA,
                cpu=cpu, started_monotonic_ns=start,
                finished_monotonic_ns=start + 500000))
    return rows


def timing_wave_rows(
        contract: Mapping[str, Any], records: Sequence[Mapping[str, Any]],
        qualification: contract_api.TimingQualification,
        ) -> Mapping[Tuple[int, int], List[int]]:
    panels = contract_api.timing_panels(contract, ARMS)
    cells = list(contract_api.iter_timing_cells(
        contract, "development", qualification))
    stable_fields = [
        field for field in contract["timing"]["cell_key"]
        if field not in (
            "replicate", "base_loss_seed", "base_cell_sha256",
            "loss_retry_offset", "loss_seed")
    ]
    identities = [contract_api.canonical_json({
        field: cell[field] for field in stable_fields
    }) for cell in cells]
    identity_order = [
        identity for cell, identity in zip(cells, identities)
        if cell["replicate"] == 0
    ]
    identity_indexes = {
        identity: index for index, identity in enumerate(identity_order)
    }
    jobs_per_wave = contract["timing"]["execution_geometry"][
        "jobs_per_wave"]
    result: Dict[Tuple[int, int], List[int]] = {}
    for row_index, row in enumerate(records):
        cell_ordinal, panel_index = divmod(row["ordinal"], len(panels))
        cell = cells[cell_ordinal]
        cohort_index = identity_indexes[identities[cell_ordinal]] * \
            len(panels) + panel_index
        wave_index = cell["replicate"] // jobs_per_wave
        result.setdefault((cohort_index, wave_index), []).append(row_index)
    return result


def write_cpu_topology(
        root: Path, values: Mapping[int, Tuple[int, int]]) -> Path:
    for cpu, (package, core) in values.items():
        topology = root / "cpu{}".format(cpu) / "topology"
        topology.mkdir(parents=True)
        (topology / "physical_package_id").write_text(
            str(package) + "\n", encoding="ascii")
        (topology / "core_id").write_text(
            str(core) + "\n", encoding="ascii")
    return root


def sampler_fixture(root: Path, *, edac: int = 0,
                    cpu_c: float = 61.0, prefix: str = "",
                    window_start_ns: int = 1000000000,
                    window_end_ns: int = 10000000000) -> Path:
    script_path = root / (prefix + "sampler.py")
    script_path.write_text("# frozen sampler fixture\n", encoding="ascii")
    csv_path = root / (prefix + "thermal.csv")
    with csv_path.open("w", encoding="ascii", newline="") as output:
        writer = csv.writer(output, lineterminator="\n")
        writer.writerow(subject.THERMAL_HEADER)
        for second in range(1, 11):
            writer.writerow([
                "2026-08-02T00:00:{:02d}.000Z".format(second),
                "{}.0".format(second), "99.0", "4000.0", str(cpu_c),
                "40.0", "40.25", "40.5", "40.75",
                "41.0", "41.25", "41.5", "41.75",
                "0", "1.0", "1.0", "1.0", str(edac), "0",
            ])
    info = csv_path.stat()
    attestation = {
        "schema": subject.SAMPLER_SCHEMA,
        "pid": 4242,
        "cpu": 127,
        "process_start_ticks": 1,
        "script_path": str(script_path),
        "script_sha256": subject._sha256_file(script_path),
        "csv_path": str(csv_path),
        "csv_device": info.st_dev,
        "csv_inode": info.st_ino,
        "window_start_monotonic_ns": window_start_ns,
        "window_end_monotonic_ns": window_end_ns,
        "terminal_status": "ok",
    }
    path = root / (prefix + "sampler.json")
    canonical_file(path, attestation)
    return path


class NativeShortScreenTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.topology_temporary = tempfile.TemporaryDirectory()
        cls.sysfs_root = write_cpu_topology(
            Path(cls.topology_temporary.name) / "sysfs",
            {cpu: (0, cpu) for cpu in list(range(9)) + [127]})
        cls.topology_patch = mock.patch.object(
            subject, "CPU_SYSFS_ROOT", cls.sysfs_root)
        cls.topology_patch.start()
        cls.contract = contract_api.load_contract()

    @classmethod
    def tearDownClass(cls) -> None:
        cls.topology_patch.stop()
        cls.topology_temporary.cleanup()

    def build_kind(self, root: Path, kind: str) -> BuiltEvidence:
        if kind == "recovery":
            return raw_v3_recovery_fixture(root, self.contract)
        if kind != "timing":
            raise AssertionError("fixture evidence kind must be recovery or timing")
        phase = "development"
        trace_path = root / "timing-traces.jsonl"
        qualification_fixture = build_qualification(root, self.contract)
        qualification = qualification_fixture.value
        hashes = list(qualification.selected_trace_sha256s)
        trace_sha = subject.publish_timing_trace_manifest(
            self.contract, phase, qualification, trace_path)
        freeze = timing_freeze_manifest(
            self.contract, phase, trace_sha, qualification)
        freeze_path = root / "timing-freeze.json"
        canonical_file(freeze_path, freeze)
        records = timing_records(self.contract, hashes, qualification)
        native_path = root / "timing-native-results.jsonl"
        canonical_jsonl(native_path, list(reversed(records)))
        sampler = sampler_fixture(
            root, window_start_ns=3000000000,
            window_end_ns=10000000000)
        qualification_sampler = root / "qualification-sampler.json"
        qualification_sampler_value = json.loads(
            sampler.read_text(encoding="utf-8"))
        qualification_sampler_value["window_start_monotonic_ns"] = 1000000000
        qualification_sampler_value["window_end_monotonic_ns"] = 2000000000
        canonical_file(qualification_sampler, qualification_sampler_value)
        return BuiltEvidence(
            freeze_path, trace_path, native_path, records,
            qualification_fixture.value, qualification_fixture.native,
            qualification_fixture.map, qualification_fixture.audit,
            qualification_fixture.map_sha256, qualification_sampler,
            sampler, qualification_fixture.execution_receipt,
            qualification_fixture.execution_receipt_sha256)

    def test_full_native_recovery_and_one_candidate_timing_publish(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampler_path = sampler_fixture(root)
            outputs = {}
            timing_expected = self.contract["timing"]["domains"][
                "development"]["expected_cells"] * len(
                    contract_api.timing_panels(self.contract, ARMS))
            for kind, expected in (
                    ("recovery", 1440), ("timing", timing_expected)):
                built = self.build_kind(root, kind)
                output = root / (kind + "-output.jsonl")
                receipt = root / (kind + "-execution.json")
                result = subject.assemble_results(
                    self.contract, kind, "development", built.freeze,
                    built.traces, built.native,
                    built.sampler or sampler_path, output, receipt,
                    verify_live_sampler=False,
                    **timing_evidence_kwargs(built))
                self.assertEqual(
                    result["execution_receipt"]["record_count"], expected)
                self.assertEqual(
                    result["execution_receipt"]["worker_cpus"],
                    list(WORKER_CPUS))
                self.assertEqual(
                    result["execution_receipt"]["thermal"]["edac_ce_max"], 0)
                self.assertTrue(output.exists())
                self.assertTrue(receipt.exists())
                revalidated = subject.validate_execution_receipt(
                    self.contract, kind, "development", built.freeze,
                    built.traces, built.native, output, receipt,
                    verify_live_sampler=False,
                    sampler_path=built.sampler or sampler_path,
                    **timing_evidence_kwargs(built))
                self.assertEqual(
                    revalidated["execution_receipt"]["record_count"], expected)
                outputs[kind] = (built, output)

            interpreters = [
                path for name in ("python3.8", "python3.12")
                for path in [shutil.which(name)] if path is not None
            ]
            self.assertTrue(interpreters, "at least one contract interpreter")
            for interpreter in interpreters:
                for kind in ("recovery", "timing"):
                    built, output = outputs[kind]
                    command = "validate-ledger" if kind == "recovery" else \
                        "validate-timing"
                    arguments = [
                        interpreter,
                        str(Path(contract_api.__file__).resolve()), command,
                        "--phase", "development", "--freeze-manifest",
                        str(built.freeze), "--trace-manifest",
                        str(built.traces),
                    ]
                    if kind == "timing":
                        arguments.extend([
                            "--timing-qualification-map",
                            str(built.qualification_map),
                            "--timing-qualification-audit",
                            str(built.qualification_audit),
                            "--timing-qualification-map-sha256",
                            str(built.qualification_map_sha256),
                        ])
                    arguments.append(str(output))
                    completed = subprocess.run(
                       arguments, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                       check=False)
                    self.assertEqual(
                        completed.returncode, 0,
                        completed.stderr.decode("utf-8", "replace"))

    def test_native_record_schema_versions_and_timing_batch_formula(self) \
            -> None:
        self.assertEqual(
            subject.RECOVERY_RECORD_SCHEMA,
            "wirehair.wh2.native-recovery-record.v1")
        self.assertEqual(
            subject.TIMING_RECORD_SCHEMA,
            "wirehair.wh2.native-timing-record.v4")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            recovery = self.build_kind(root, "recovery").records
            self.assertEqual(
                sum(row["schema"] == subject.RECOVERY_RECORD_SCHEMA
                    for row in recovery), 720)
            self.assertEqual(
                sum(row["schema"] == subject.RAW_RECOVERY_RECORD_SCHEMA
                    for row in recovery), 720)
            self.assertTrue(all(
                "invocations_per_slot" not in row["payload"]
                for row in recovery))

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            built = self.build_kind(root, "timing")
            timing = built.records
            freeze = contract_api.load_freeze_manifest(
                self.contract, "development", built.freeze, "timing",
                built.qualification)
            self.assertTrue(all(
                row["schema"] == subject.TIMING_RECORD_SCHEMA
                for row in timing))
            for row in timing:
                payload = row["payload"]
                expected = max(2, (65536 + payload["K"] - 1) // payload["K"])
                self.assertEqual(payload["invocations_per_slot"], expected)
                self.assertEqual(
                    payload["interleave_policy"],
                    "self-counterbalanced-repeat-major-v1")
            legacy = copy.deepcopy(timing)
            legacy[0]["schema"] = "wirehair.wh2.native-timing-record.v2"
            with self.assertRaises(subject.NativeEvidenceError):
                subject._validate_native_records(
                    self.contract, freeze, "timing", "development", legacy,
                    built.qualification)

    def test_timing_receipt_binds_qualification_and_separate_thermal_windows(
            self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            built = self.build_kind(root, "timing")
            result_path = root / "timing-results.jsonl"
            receipt_path = root / "timing-execution.json"
            assembled = subject.assemble_results(
                self.contract, "timing", "development", built.freeze,
                built.traces, built.native, built.sampler, result_path,
                receipt_path, verify_live_sampler=False,
                **timing_evidence_kwargs(built))
            receipt = assembled["execution_receipt"]
            self.assertEqual(set(receipt), subject.TIMING_EXECUTION_FIELDS)
            self.assertEqual(receipt["schema"], subject.TIMING_EXECUTION_SCHEMA)
            self.assertEqual(
                receipt["timing_base_domain_sha256"],
                built.qualification.base_domain_sha256)
            self.assertEqual(
                receipt["timing_qualified_domain_sha256"],
                built.qualification.qualified_domain_sha256)
            self.assertEqual(
                receipt["timing_qualification_map_sha256"],
                built.qualification_map_sha256)
            self.assertEqual(
                receipt["timing_qualification_audit_sha256"],
                built.qualification.qualification_audit_sha256)
            self.assertEqual(
                receipt[
                    "timing_qualification_execution_receipt_sha256"],
                built.qualification_execution_receipt_sha256)
            self.assertEqual(
                receipt["qualification_worker_cpus"],
                list(QUALIFICATION_CPUS))
            self.assertEqual(
                receipt["qualification_thermal"][
                    "window_end_monotonic_ns"], 2000000000)
            self.assertEqual(
                receipt["thermal"]["window_start_monotonic_ns"],
                3000000000)
            self.assertLess(
                receipt["qualification_thermal"][
                    "window_end_monotonic_ns"],
                receipt["thermal"]["window_start_monotonic_ns"])
            self.assertEqual(
                receipt["sibling_isolation"],
                timing_evidence_kwargs(built)["sibling_isolation"])

            for missing in (
                    "timing_qualification_execution_receipt_path",
                    "timing_qualification_execution_receipt_sha256",
                    "sibling_isolation"):
                kwargs = timing_evidence_kwargs(built)
                kwargs[missing] = None
                missing_result = root / (missing + "-result.jsonl")
                missing_receipt = root / (missing + "-receipt.json")
                with self.subTest(missing=missing), self.assertRaises(
                        subject.NativeEvidenceError):
                    subject.assemble_results(
                        self.contract, "timing", "development", built.freeze,
                        built.traces, built.native, built.sampler,
                        missing_result, missing_receipt,
                        verify_live_sampler=False, **kwargs)
                self.assertFalse(missing_result.exists())
                self.assertFalse(missing_receipt.exists())

            changed = copy.deepcopy(receipt)
            changed["sibling_isolation"]["checks"][0][
                "foreign_task_count"] = 1
            changed["sibling_isolation"]["checks_sha256"] = \
                contract_api.sha256_json(
                    changed["sibling_isolation"]["checks"])
            changed["receipt_sha256"] = contract_api.sha256_json({
                key: item for key, item in changed.items()
                if key != "receipt_sha256"
            })
            mutant = root / "sibling-isolation-receipt.json"
            canonical_file(mutant, changed)
            with self.assertRaisesRegex(
                    subject.NativeEvidenceError, "sibling-isolation check"):
                subject.validate_execution_receipt(
                    self.contract, "timing", "development",
                    built.freeze, built.traces, built.native, result_path,
                    mutant, verify_live_sampler=False,
                    sampler_path=built.sampler,
                    **timing_evidence_kwargs(built))

            incomplete = copy.deepcopy(receipt)
            del incomplete["sibling_isolation"]["checks"][2]
            incomplete["sibling_isolation"]["check_count"] -= 1
            incomplete["sibling_isolation"]["checks_sha256"] = \
                contract_api.sha256_json(
                    incomplete["sibling_isolation"]["checks"])
            incomplete["receipt_sha256"] = contract_api.sha256_json({
                key: item for key, item in incomplete.items()
                if key != "receipt_sha256"
            })
            incomplete_path = root / "incomplete-sibling-wave-receipt.json"
            canonical_file(incomplete_path, incomplete)
            with self.assertRaisesRegex(
                    subject.NativeEvidenceError, "wave checks are incomplete"):
                subject.validate_execution_receipt(
                    self.contract, "timing", "development",
                    built.freeze, built.traces, built.native, result_path,
                    incomplete_path, verify_live_sampler=False,
                    sampler_path=built.sampler,
                    **timing_evidence_kwargs(built))

            missed_wave = copy.deepcopy(receipt)
            first_before = missed_wave["sibling_isolation"]["checks"][2]
            first_after = missed_wave["sibling_isolation"]["checks"][3]
            first_after["monotonic_ns"] = first_before["monotonic_ns"] + 1
            missed_wave["sibling_isolation"]["checks_sha256"] = \
                contract_api.sha256_json(
                    missed_wave["sibling_isolation"]["checks"])
            missed_wave["receipt_sha256"] = contract_api.sha256_json({
                key: item for key, item in missed_wave.items()
                if key != "receipt_sha256"
            })
            missed_wave_path = root / "missed-native-wave-receipt.json"
            canonical_file(missed_wave_path, missed_wave)
            with self.assertRaisesRegex(
                    subject.NativeEvidenceError, "misses its native wave"):
                subject.validate_execution_receipt(
                    self.contract, "timing", "development",
                    built.freeze, built.traces, built.native, result_path,
                    missed_wave_path, verify_live_sampler=False,
                    sampler_path=built.sampler,
                    **timing_evidence_kwargs(built))

            mutations = {
                "base_domain": ("timing_base_domain_sha256", digest("bad")),
                "qualified_domain": (
                    "timing_qualified_domain_sha256", digest("bad")),
                "map": ("timing_qualification_map_sha256", digest("bad")),
                "audit": (
                    "timing_qualification_audit_sha256", digest("bad")),
                "native_stream": (
                    "timing_qualification_native_stream_sha256",
                    digest("bad")),
                "qualification_execution": (
                    "timing_qualification_execution_receipt_sha256",
                    digest("bad")),
                "cpu_roster": (
                    "qualification_worker_cpus",
                    list(reversed(QUALIFICATION_CPUS))),
            }
            for name, (field, value) in mutations.items():
                changed = copy.deepcopy(receipt)
                changed[field] = value
                changed["receipt_sha256"] = contract_api.sha256_json({
                    key: item for key, item in changed.items()
                    if key != "receipt_sha256"
                })
                mutant = root / (name + "-receipt.json")
                canonical_file(mutant, changed)
                with self.subTest(name=name), self.assertRaises(
                        subject.NativeEvidenceError):
                    subject.validate_execution_receipt(
                        self.contract, "timing", "development",
                        built.freeze, built.traces, built.native, result_path,
                        mutant, verify_live_sampler=False,
                        sampler_path=built.sampler,
                        **timing_evidence_kwargs(built))

    def test_timing_receipt_rejects_sampler_splice_and_touching_windows(
            self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            built = self.build_kind(root, "timing")
            result_path = root / "timing-results.jsonl"
            receipt_path = root / "timing-execution.json"
            receipt = subject.assemble_results(
                self.contract, "timing", "development", built.freeze,
                built.traces, built.native, built.sampler, result_path,
                receipt_path, verify_live_sampler=False,
                **timing_evidence_kwargs(built))["execution_receipt"]

            cases = []
            spliced = sampler_fixture(
                root, prefix="spliced-", window_start_ns=1000000000,
                window_end_ns=2000000000)
            spliced_value = json.loads(spliced.read_text(encoding="utf-8"))
            spliced_thermal = subject._thermal_window(
                spliced_value, 1200000000, 1300000000,
                QUALIFICATION_CPUS, False)
            cases.append(("spliced_identity", spliced, spliced_thermal))

            touching = root / "touching-qualification-sampler.json"
            touching_value = json.loads(
                built.qualification_sampler.read_text(encoding="utf-8"))
            touching_value["window_end_monotonic_ns"] = 3000000000
            canonical_file(touching, touching_value)
            touching_thermal = subject._thermal_window(
                touching_value, 1200000000, 1300000000,
                QUALIFICATION_CPUS, False)
            cases.append(("touching_windows", touching, touching_thermal))

            for name, qualification_sampler, qualification_thermal in cases:
                changed = copy.deepcopy(receipt)
                changed["qualification_thermal"] = qualification_thermal
                changed["receipt_sha256"] = contract_api.sha256_json({
                    key: value for key, value in changed.items()
                    if key != "receipt_sha256"
                })
                mutant = root / (name + "-receipt.json")
                canonical_file(mutant, changed)
                kwargs = timing_evidence_kwargs(built)
                kwargs["timing_qualification_sampler_path"] = \
                    qualification_sampler
                with self.subTest(name=name), self.assertRaises(
                        subject.NativeEvidenceError):
                    subject.validate_execution_receipt(
                        self.contract, "timing", "development",
                        built.freeze, built.traces, built.native, result_path,
                        mutant, verify_live_sampler=False,
                        sampler_path=built.sampler, **kwargs)

    def test_timing_receipt_rejects_qualification_and_trace_substitution(
            self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            built = self.build_kind(root, "timing")
            result_path = root / "timing-results.jsonl"
            receipt_path = root / "timing-execution.json"
            subject.assemble_results(
                self.contract, "timing", "development", built.freeze,
                built.traces, built.native, built.sampler, result_path,
                receipt_path, verify_live_sampler=False,
                **timing_evidence_kwargs(built))

            artifacts = (
                ("map", built.qualification_map),
                ("audit", built.qualification_audit),
                ("trace", built.traces),
                ("native", built.qualification_native),
                ("qualification_receipt",
                 built.qualification_execution_receipt),
            )
            for name, path in artifacts:
                original = path.read_bytes()
                if name in ("map", "qualification_receipt"):
                    value = json.loads(original.decode("utf-8"))
                    if name == "map":
                        value["source_git_commit"] = "2" * 40
                    else:
                        value["qualification_attempt_count"] += 1
                        value["receipt_sha256"] = contract_api.sha256_json({
                            key: item for key, item in value.items()
                            if key != "receipt_sha256"
                        })
                    canonical_file(path, value)
                else:
                    rows = list(contract_api._parse_canonical_jsonl(
                        path, name + " mutation fixture"))
                    if name in ("audit", "trace"):
                        rows[0] = dict(rows[0])
                        rows[0]["trace_sha256"] = digest(
                            name + " substitution")
                    else:
                        rows.reverse()
                    canonical_jsonl(path, rows)
                try:
                    with self.subTest(name=name), self.assertRaises(
                            subject.NativeEvidenceError):
                        subject.validate_execution_receipt(
                            self.contract, "timing", "development",
                            built.freeze, built.traces, built.native,
                            result_path, receipt_path,
                            verify_live_sampler=False,
                            sampler_path=built.sampler,
                            **timing_evidence_kwargs(built))
                finally:
                    path.write_bytes(original)

            substitute_root = root / "substitute"
            substitute_root.mkdir()
            offsets = [0] * len(built.qualification.retry_offsets)
            offsets[0] = 1
            substitute = build_qualification(
                substitute_root, self.contract, offsets)
            kwargs = timing_evidence_kwargs(built)
            kwargs.update({
                "timing_qualification_map_path": substitute.map,
                "timing_qualification_audit_path": substitute.audit,
                "timing_qualification_map_sha256": substitute.map_sha256,
                "timing_qualification_native_path": substitute.native,
                "timing_qualification_execution_receipt_path":
                    substitute.execution_receipt,
                "timing_qualification_execution_receipt_sha256":
                    substitute.execution_receipt_sha256,
            })
            with self.assertRaises(subject.NativeEvidenceError):
                subject.validate_execution_receipt(
                    self.contract, "timing", "development", built.freeze,
                    built.traces, built.native, result_path, receipt_path,
                    verify_live_sampler=False, sampler_path=built.sampler,
                    **kwargs)

    def test_timing_terminal_reopen_rejects_qualification_evidence_races(
            self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            built = self.build_kind(root, "timing")
            sampler_value = json.loads(
                built.qualification_sampler.read_text(encoding="utf-8"))
            csv_path = Path(sampler_value["csv_path"])
            originals = {
                "receipt": built.qualification_execution_receipt.read_bytes(),
                "sampler": built.qualification_sampler.read_bytes(),
                "csv": csv_path.read_bytes(),
            }

            def mutate(name: str) -> None:
                if name == "receipt":
                    value = json.loads(
                        built.qualification_execution_receipt.read_text(
                            encoding="utf-8"))
                    value["qualification_attempt_count"] += 1
                    value["receipt_sha256"] = contract_api.sha256_json({
                        key: item for key, item in value.items()
                        if key != "receipt_sha256"
                    })
                    canonical_file(
                        built.qualification_execution_receipt, value)
                elif name == "sampler":
                    value = json.loads(
                        built.qualification_sampler.read_text(
                            encoding="utf-8"))
                    value["window_start_monotonic_ns"] = 0
                    canonical_file(built.qualification_sampler, value)
                else:
                    lines = csv_path.read_text(encoding="ascii").splitlines()
                    fields = lines[2].split(",")
                    fields[4] = "63.0"
                    lines[2] = ",".join(fields)
                    csv_path.write_text(
                        "\n".join(lines) + "\n", encoding="ascii")

            original_validate = contract_api.validate_timing_receipt
            for name in ("receipt", "sampler", "csv"):
                output = root / (name + "-timing-result.jsonl")
                receipt = root / (name + "-timing-execution.json")

                def validate_then_mutate(*args, _name=name, **kwargs):
                    summary = original_validate(*args, **kwargs)
                    mutate(_name)
                    return summary

                try:
                    with self.subTest(name=name), mock.patch.object(
                            contract_api, "validate_timing_receipt",
                            side_effect=validate_then_mutate), \
                            self.assertRaises(subject.NativeEvidenceError):
                        subject.assemble_results(
                            self.contract, "timing", "development",
                            built.freeze, built.traces, built.native,
                            built.sampler, output, receipt,
                            verify_live_sampler=False,
                            **timing_evidence_kwargs(built))
                    self.assertFalse(output.exists())
                    self.assertFalse(receipt.exists())
                finally:
                    built.qualification_execution_receipt.write_bytes(
                        originals["receipt"])
                    built.qualification_sampler.write_bytes(
                        originals["sampler"])
                    csv_path.write_bytes(originals["csv"])

    def test_qualification_canonicalizes_retries_and_exercises_all_cpus(
            self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            offsets = [0] * self.contract["timing"]["domains"][
                "development"]["expected_cells"]
            offsets[0] = 1
            offsets[1] = 2
            built = build_qualification(root, self.contract, offsets)
            self.assertEqual(list(built.value.retry_offsets), offsets)
            self.assertEqual(
                built.metadata["qualification_worker_cpus"],
                list(QUALIFICATION_CPUS))
            self.assertEqual(
                [value["cpu"] for value in
                 built.metadata["qualification_workers"]],
                list(QUALIFICATION_CPUS))
            execution_receipt = json.loads(
                built.execution_receipt.read_text(encoding="utf-8"))
            self.assertEqual(
                set(execution_receipt),
                subject.TIMING_QUALIFICATION_EXECUTION_FIELDS)
            self.assertEqual(
                execution_receipt["schema"],
                subject.TIMING_QUALIFICATION_EXECUTION_SCHEMA)
            self.assertEqual(
                execution_receipt["receipt_sha256"],
                built.execution_receipt_sha256)
            self.assertEqual(
                execution_receipt["qualification_attempt_count"],
                len(built.records))
            self.assertEqual(
                execution_receipt["qualification_allowed_cpus"],
                list(QUALIFICATION_CPUS))
            reopened, reopened_metadata = \
                subject.load_timing_qualification_execution_receipt(
                    self.contract, "development", built.value, built.native,
                    built.execution_receipt,
                    expected_receipt_sha256=
                        built.execution_receipt_sha256,
                    expected_cpus=QUALIFICATION_CPUS)
            self.assertEqual(reopened, execution_receipt)
            self.assertEqual(reopened_metadata, built.metadata)

            audit = list(contract_api._parse_canonical_jsonl(
                built.audit, "qualification audit fixture"))
            self.assertEqual(
                [(row["ordinal"], row["loss_retry_offset"])
                 for row in audit[:5]],
                [(0, 0), (0, 1), (1, 0), (1, 1), (1, 2)])
            first_attempts = sorted(
                (row for row in built.records
                 if row["payload"]["ordinal"] == 0),
                key=lambda row: row["payload"]["loss_retry_offset"])
            self.assertEqual(len(first_attempts), 2)
            self.assertEqual(
                {row["message_sha256"] for row in first_attempts},
                {first_attempts[0]["message_sha256"]})
            self.assertEqual(
                len({row["payload"]["cell_sha256"]
                     for row in first_attempts}), 2)
            self.assertEqual(
                len({row["payload"]["trace_sha256"]
                     for row in first_attempts}), 2)
            self.assertEqual(
                len({row["work_sha256"] for row in first_attempts}), 2)
            self.assertTrue(all(
                set(row["payload"]) ==
                subject.TIMING_QUALIFICATION_PAYLOAD_FIELDS
                for row in built.records))

    def test_qualification_rejects_gaps_duplicates_speculation_and_fatal_rows(
            self) -> None:
        count = self.contract["timing"]["domains"]["development"][
            "expected_cells"]
        base_offsets = [0] * count
        cases: Dict[str, List[Dict[str, Any]]] = {}

        gap_offsets = list(base_offsets)
        gap_offsets[0] = 2
        cases["gap"] = [
            row for row in qualification_records(self.contract, gap_offsets)
            if row["ordinal"] != 1]

        duplicate = qualification_records(self.contract)
        duplicate.append(copy.deepcopy(duplicate[0]))
        cases["duplicate"] = duplicate

        speculative_offsets = list(base_offsets)
        speculative_offsets[0] = 1
        speculative = qualification_records(
            self.contract, speculative_offsets)
        speculative[0]["payload"]["wirehair2_head_outcome"] = "success"
        speculative[0]["payload"]["wirehair2_head_decoded_extra"] = 4
        cases["speculative"] = speculative

        fatal = qualification_records(self.contract)
        fatal[0]["payload"]["wirehair2_head_outcome"] = "fatal"
        fatal[0]["payload"]["wirehair2_head_decoded_extra"] = None
        cases["fatal"] = fatal

        exhausted_offsets = list(base_offsets)
        exhausted_offsets[0] = 255
        exhausted = qualification_records(self.contract, exhausted_offsets)
        exhausted[255]["payload"]["wirehair2_head_outcome"] = \
            "need_more_at_bound"
        exhausted[255]["payload"]["wirehair2_head_decoded_extra"] = None
        cases["retry255"] = exhausted

        malformed = qualification_records(self.contract)
        malformed[0]["worker_pid"] = True
        cases["malformed_provenance"] = malformed

        cpu_subset = qualification_records(self.contract)
        for row in cpu_subset:
            if row["cpu"] == QUALIFICATION_CPUS[-1]:
                row["cpu"] = QUALIFICATION_CPUS[0]
                row["worker_pid"] = 2000 + QUALIFICATION_CPUS[0]
        cases["cpu_subset"] = cpu_subset

        for name, rows in cases.items():
            with self.subTest(name=name), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                native = root / "qualification-native.jsonl"
                audit = root / "qualification-audit.jsonl"
                map_path = root / "qualification-map.json"
                execution_receipt = root / "qualification-execution.json"
                canonical_jsonl(native, list(reversed(rows)))
                with self.assertRaises(subject.NativeEvidenceError):
                    subject.assemble_timing_qualification(
                        self.contract, "development", native, audit,
                        map_path, execution_receipt, SOURCE_COMMIT,
                        qualification_controls(), QUALIFICATION_CPUS,
                        verify_live_workers=False)
                self.assertFalse(audit.exists())
                self.assertFalse(map_path.exists())
                self.assertFalse(execution_receipt.exists())

    def test_qualification_execution_receipt_rejects_tamper_and_post_q_mutation(
            self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            built = build_qualification(root, self.contract)
            receipt = json.loads(
                built.execution_receipt.read_text(encoding="utf-8"))
            mutations = {}
            self_hash = copy.deepcopy(receipt)
            self_hash["receipt_sha256"] = digest("bad self hash")
            mutations["self_hash"] = self_hash
            for name, field, value in (
                    ("allowed_cpu_subset", "qualification_allowed_cpus",
                     list(QUALIFICATION_CPUS[:-1])),
                    ("observed_cpu_subset", "qualification_worker_cpus",
                     list(QUALIFICATION_CPUS[:-1])),
                    ("map_hash", "timing_qualification_map_sha256",
                     digest("bad map")),
                    ("native_hash",
                     "timing_qualification_native_stream_sha256",
                     digest("bad native"))):
                changed = copy.deepcopy(receipt)
                changed[field] = value
                changed["receipt_sha256"] = contract_api.sha256_json({
                    key: item for key, item in changed.items()
                    if key != "receipt_sha256"
                })
                mutations[name] = changed
            worker = copy.deepcopy(receipt)
            worker["qualification_workers"][0]["pid"] += 1
            worker["receipt_sha256"] = contract_api.sha256_json({
                key: value for key, value in worker.items()
                if key != "receipt_sha256"
            })
            mutations["worker_roster"] = worker

            def add_rehashed(name: str, changed: Dict[str, Any]) -> None:
                changed["receipt_sha256"] = contract_api.sha256_json({
                    key: value for key, value in changed.items()
                    if key != "receipt_sha256"
                })
                mutations[name] = changed

            allowed_cpu_alias = copy.deepcopy(receipt)
            allowed_cpu_alias["qualification_allowed_cpus"][0] = 0.0
            add_rehashed("allowed_cpu_numeric_alias", allowed_cpu_alias)
            for name, value in (
                    ("allowed_cpu_list_container", []),
                    ("allowed_cpu_object_container", {})):
                changed = copy.deepcopy(receipt)
                changed["qualification_allowed_cpus"][0] = value
                add_rehashed(name, changed)
            worker_cpu_alias = copy.deepcopy(receipt)
            worker_cpu_alias["qualification_worker_cpus"][0] = 0.0
            add_rehashed("worker_cpu_numeric_alias", worker_cpu_alias)
            worker_field_aliases = (
                ("worker_cpu_bool_alias", "cpu", False),
                ("worker_pid_numeric_alias", "pid", float(
                    receipt["qualification_workers"][0]["pid"])),
                ("worker_start_ticks_numeric_alias", "process_start_ticks",
                 float(receipt["qualification_workers"][0]
                       ["process_start_ticks"])),
            )
            for name, field, value in worker_field_aliases:
                changed = copy.deepcopy(receipt)
                changed["qualification_workers"][0][field] = value
                add_rehashed(name, changed)
            numeric_aliases = (
                ("attempt_count_numeric_alias",
                 "qualification_attempt_count"),
                ("worker_start_numeric_alias",
                 "qualification_worker_start_monotonic_ns"),
                ("worker_end_numeric_alias",
                 "qualification_worker_end_monotonic_ns"),
            )
            for name, field in numeric_aliases:
                changed = copy.deepcopy(receipt)
                changed[field] = float(changed[field])
                add_rehashed(name, changed)
            binary_container = copy.deepcopy(receipt)
            binary_container["qualification_worker_binary_sha256s"][0] = {}
            add_rehashed("worker_binary_object_container", binary_container)

            for name, value in mutations.items():
                mutant = root / (name + "-qualification-execution.json")
                canonical_file(mutant, value)
                with self.subTest(name=name), self.assertRaises(
                        subject.NativeEvidenceError):
                    subject.load_timing_qualification_execution_receipt(
                        self.contract, "development", built.value,
                        built.native, mutant,
                        expected_receipt_sha256=value["receipt_sha256"])

            for name, value in (("list", []), ("object", {})):
                expected_cpus = [value, *QUALIFICATION_CPUS[1:]]
                with self.subTest(expected_cpu_container=name), \
                        self.assertRaises(subject.NativeEvidenceError):
                    subject.load_timing_qualification_execution_receipt(
                        self.contract, "development", built.value,
                        built.native, built.execution_receipt,
                        expected_receipt_sha256=
                            built.execution_receipt_sha256,
                        expected_cpus=expected_cpus)

            original_native = built.native.read_bytes()
            rows = list(contract_api._parse_canonical_jsonl(
                built.native, "post-Q mutation fixture"))
            canonical_jsonl(built.native, list(reversed(rows)))
            try:
                with self.assertRaises(subject.NativeEvidenceError):
                    subject.load_timing_qualification_execution_receipt(
                        self.contract, "development", built.value,
                        built.native, built.execution_receipt,
                        expected_receipt_sha256=
                            built.execution_receipt_sha256)
            finally:
                built.native.write_bytes(original_native)

    def test_qualification_receipt_commit_rejects_dependency_mutation(
            self) -> None:
        for dependency in ("map", "audit"):
            with self.subTest(dependency=dependency), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                native = root / "qualification-native.jsonl"
                audit = root / "qualification-audit.jsonl"
                map_path = root / "qualification-map.json"
                execution_receipt = root / "qualification-execution.json"
                canonical_jsonl(native, qualification_records(self.contract))
                real_publish = subject._publish_new
                publish_calls = 0

                def mutate_before_receipt_publish(
                        staged, destination, expected_identity=None):
                    nonlocal publish_calls
                    publish_calls += 1
                    if publish_calls == 3:
                        self.assertEqual(destination, execution_receipt)
                        if dependency == "map":
                            value = json.loads(
                                map_path.read_text(encoding="utf-8"))
                            value["source_git_commit"] = "2" * 40
                            canonical_file(map_path, value)
                        else:
                            rows = list(contract_api._parse_canonical_jsonl(
                                audit, "qualification race fixture"))
                            rows[0] = dict(rows[0])
                            rows[0]["trace_sha256"] = digest(
                                "mutated qualification audit")
                            canonical_jsonl(audit, rows)
                    real_publish(staged, destination, expected_identity)

                with mock.patch.object(
                        subject, "_publish_new",
                        side_effect=mutate_before_receipt_publish), \
                        self.assertRaises((
                            subject.NativeEvidenceError,
                            contract_api.ContractError)):
                    subject.assemble_timing_qualification(
                        self.contract, "development", native, audit,
                        map_path, execution_receipt, SOURCE_COMMIT,
                        qualification_controls(), QUALIFICATION_CPUS,
                        verify_live_workers=False)
                self.assertEqual(publish_calls, 3)
                self.assertTrue(execution_receipt.exists())
                receipt = json.loads(
                    execution_receipt.read_text(encoding="utf-8"))
                with self.assertRaises(contract_api.ContractError):
                    contract_api.load_timing_qualification_map(
                        self.contract, "development", map_path, audit,
                        receipt["timing_qualification_map_sha256"])

    def test_recovery_schema_dispatch_rejects_non_string_arm(self) -> None:
        freeze = {
            "schema": contract_api.RAW_FREEZE_SCHEMA,
            "arms_by_name": {
                "wirehair2_raw_fixture": {
                    "construction_seed_basis":
                        contract_api.RAW_CONSTRUCTION_SEED_BASIS,
                },
            },
        }
        with self.assertRaisesRegex(
                subject.NativeEvidenceError, "arm must be a string"):
            subject.recovery_record_schema_fields(freeze, {"arm": []})
        schema, fields = subject.recovery_record_schema_fields(
            freeze, {"arm": "wirehair2_raw_fixture"})
        self.assertEqual(schema, subject.RAW_RECOVERY_RECORD_SCHEMA)
        self.assertEqual(fields, contract_api.RAW_RECOVERY_RECORD_FIELDS)

    def test_raw_v3_mixed_stream_receipt_downgrade_and_terminal_alias(
            self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            built = raw_v3_recovery_fixture(root, self.contract)
            sampler = sampler_fixture(root)
            result = root / "raw-v3-result.jsonl"
            receipt_path = root / "raw-v3-execution.json"
            assembled = subject.assemble_results(
                self.contract, "recovery", "development", built.freeze,
                built.traces, built.native, sampler, result, receipt_path,
                verify_live_sampler=False)
            receipt = assembled["execution_receipt"]
            self.assertEqual(receipt["schema"], subject.RAW_EXECUTION_SCHEMA)
            self.assertEqual(receipt["record_count"], 1440)
            self.assertEqual(
                {row["schema"] for row in built.records},
                {subject.RECOVERY_RECORD_SCHEMA,
                 subject.RAW_RECOVERY_RECORD_SCHEMA})
            subject.validate_execution_receipt(
                self.contract, "recovery", "development", built.freeze,
                built.traces, built.native, result, receipt_path,
                verify_live_sampler=False, sampler_path=sampler)

            downgraded_records = copy.deepcopy(built.records)
            downgraded_records[2]["schema"] = \
                subject.LEGACY_RAW_RECOVERY_RECORD_SCHEMA
            downgraded_native = root / "downgraded-raw-native.jsonl"
            canonical_jsonl(downgraded_native, downgraded_records)
            with self.assertRaisesRegex(
                    subject.NativeEvidenceError, "unknown schema"):
                subject.assemble_results(
                    self.contract, "recovery", "development", built.freeze,
                    built.traces, downgraded_native, sampler,
                    root / "downgraded-result.jsonl",
                    root / "downgraded-execution.json",
                    verify_live_sampler=False)

            downgraded_receipt = copy.deepcopy(receipt)
            downgraded_receipt["schema"] = subject.LEGACY_RAW_EXECUTION_SCHEMA
            downgraded_receipt["receipt_sha256"] = contract_api.sha256_json({
                key: value for key, value in downgraded_receipt.items()
                if key != "receipt_sha256"
            })
            downgraded_receipt_path = root / "downgraded-receipt.json"
            canonical_file(downgraded_receipt_path, downgraded_receipt)
            with self.assertRaisesRegex(
                    subject.NativeEvidenceError, "unknown schema"):
                subject.validate_execution_receipt(
                    self.contract, "recovery", "development", built.freeze,
                    built.traces, built.native, result,
                    downgraded_receipt_path, verify_live_sampler=False,
                    sampler_path=sampler)

            aliased_receipt = copy.deepcopy(receipt)
            aliased_receipt["record_count"] = float(
                aliased_receipt["record_count"])
            aliased_receipt["receipt_sha256"] = contract_api.sha256_json({
                key: value for key, value in aliased_receipt.items()
                if key != "receipt_sha256"
            })
            aliased_path = root / "aliased-receipt.json"
            canonical_file(aliased_path, aliased_receipt)
            with self.assertRaisesRegex(
                    subject.NativeEvidenceError,
                    "differs from native evidence"):
                subject.validate_execution_receipt(
                    self.contract, "recovery", "development", built.freeze,
                    built.traces, built.native, result, aliased_path,
                    verify_live_sampler=False, sampler_path=sampler)

            original_validate = contract_api.validate_ledger

            def validate_then_alias(*args, **kwargs):
                summary = original_validate(*args, **kwargs)
                value = json.loads(
                    receipt_path.read_text(encoding="utf-8"))
                value["record_count"] = float(value["record_count"])
                canonical_file(receipt_path, value)
                return summary

            with mock.patch.object(
                    contract_api, "validate_ledger",
                    side_effect=validate_then_alias), \
                    self.assertRaisesRegex(
                        subject.NativeEvidenceError,
                        "execution receipt changed"):
                subject.validate_execution_receipt(
                    self.contract, "recovery", "development", built.freeze,
                    built.traces, built.native, result, receipt_path,
                    verify_live_sampler=False, sampler_path=sampler)

    def test_terminal_timing_geometry_rejects_placement_and_barrier_drift(
            self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            built = self.build_kind(root, "timing")
            records = built.records
            freeze = contract_api.load_freeze_manifest(
                self.contract, "development", built.freeze, "timing",
                built.qualification)
            waves = timing_wave_rows(
                self.contract, records, built.qualification)

            def validate(
                    values: Sequence[Mapping[str, Any]],
                    frozen: Mapping[str, Any] = freeze,
                    sysfs_root: Optional[Path] = None) -> Any:
                return subject._validate_native_records(
                    self.contract, frozen, "timing", "development", values,
                    built.qualification, sysfs_root=sysfs_root)
            rounds = self.contract["timing"]["domains"]["development"][
                "independent_rounds"]
            expected_cells = self.contract["timing"]["domains"][
                "development"]["expected_cells"]
            self.assertEqual(len(records), expected_cells * 11)
            self.assertEqual(len(waves), 264 * rounds)
            self.assertTrue(all(len(indexes) == 8
                                for indexes in waves.values()))
            # The v4 domain has 25,344 rows.  Rebuilding the 2,304-cell index
            # for every row is both unnecessary and pathologically expensive.
            with mock.patch.object(
                    contract_api, "_timing_cell_indexes",
                    wraps=contract_api._timing_cell_indexes) as cell_indexes:
                validate(records)
            self.assertEqual(cell_indexes.call_count, 1)

            first_wave = waves[(0, 0)]
            second_wave = waves[(0, 1)]
            first_by_replicate = {
                records[index]["payload"]["replicate"]: index
                for index in first_wave
            }
            second_by_replicate = {
                records[index]["payload"]["replicate"]: index
                for index in second_wave
            }
            self.assertEqual(
                [records[first_by_replicate[replicate]]["cpu"]
                 for replicate in range(8)], list(WORKER_CPUS))
            self.assertEqual(
                [records[second_by_replicate[replicate]]["cpu"]
                 for replicate in range(8, 16)],
                [1, 2, 3, 4, 5, 6, 7, 0])

            cpu_remap = copy.deepcopy(records)
            for row in cpu_remap:
                if row["cpu"] == 0:
                    row["cpu"] = 1
                elif row["cpu"] == 1:
                    row["cpu"] = 0
            with self.assertRaisesRegex(
                    subject.NativeEvidenceError, "worker-slot placement"):
                validate(cpu_remap)

            swapped = copy.deepcopy(records)
            panels = contract_api.timing_panels(self.contract, ARMS)
            cells = list(contract_api.iter_timing_cells(
                self.contract, "development", built.qualification))
            first_identity = {
                field: cells[0][field]
                for field in self.contract["timing"]["cell_key"]
                if field not in (
                    "replicate", "base_loss_seed", "base_cell_sha256",
                    "loss_retry_offset", "loss_seed")
            }
            matching_cells = [
                index for index, cell in enumerate(cells)
                if all(cell[field] == value
                       for field, value in first_identity.items())
            ]
            cell0, cell1 = matching_cells[:2]
            identity_fields = (
                "ordinal", "payload", "message_sha256", "work_sha256",
            )
            for panel_index in range(len(panels)):
                left = cell0 * len(panels) + panel_index
                right = cell1 * len(panels) + panel_index
                for field in identity_fields:
                    swapped[left][field], swapped[right][field] = \
                        swapped[right][field], swapped[left][field]
            with self.assertRaisesRegex(
                    subject.NativeEvidenceError, "worker-slot placement"):
                validate(swapped)

            ordered_wave_keys = [
                (cohort, wave)
                for wave in range(rounds) for cohort in range(264)
            ]
            prior_indexes = waves[ordered_wave_keys[0]]
            next_indexes = waves[ordered_wave_keys[1]]
            prior_finish = max(
                records[index]["finished_monotonic_ns"]
                for index in prior_indexes)
            next_min_index = min(
                next_indexes,
                key=lambda index: records[index]["started_monotonic_ns"])
            equality = copy.deepcopy(records)
            equality[next_min_index]["started_monotonic_ns"] = prior_finish
            validate(equality)
            overlap = copy.deepcopy(equality)
            overlap[next_min_index]["started_monotonic_ns"] = prior_finish - 1
            with self.assertRaisesRegex(
                    subject.NativeEvidenceError, "overlap"):
                validate(overlap)

            # The superseded v3 scheduler completed every round for one
            # cohort before moving to the next cohort.  Even if each old wave
            # is internally well formed, its chronology is not a v4 proof.
            cohort_major = copy.deepcopy(records)
            for (cohort_index, round_index), indexes in waves.items():
                old_wave = cohort_index * rounds + round_index
                for index in indexes:
                    position = cohort_major[index]["payload"]["replicate"] % 8
                    start = 2000000000 + old_wave * 1000000 + position * 10
                    cohort_major[index]["started_monotonic_ns"] = start
                    cohort_major[index]["finished_monotonic_ns"] = \
                        start + 500000
            with self.assertRaisesRegex(
                    subject.NativeEvidenceError, "overlap"):
                validate(cohort_major)

            with self.assertRaises(subject.NativeEvidenceError):
                validate([row for index, row in enumerate(records)
                          if index not in set(first_wave)])
            with self.assertRaises(subject.NativeEvidenceError):
                validate(records + [copy.deepcopy(records[index])
                                    for index in first_wave])

            for roster in (
                    list(range(7)), list(range(9)),
                    [1, 0, 2, 3, 4, 5, 6, 7]):
                changed_freeze = copy.deepcopy(freeze)
                changed_freeze["cpu_affinity"] = roster
                with self.subTest(roster=roster), self.assertRaisesRegex(
                        subject.NativeEvidenceError,
                        "exactly eight sorted worker CPUs"):
                    validate(records, changed_freeze)

            sibling_root = write_cpu_topology(
                root / "sibling-sysfs",
                {**{cpu: (0, cpu) for cpu in range(9)},
                 64: (0, 0), 127: (0, 127)})
            sibling_records = copy.deepcopy(records)
            for row in sibling_records:
                if row["cpu"] == 7:
                    row["cpu"] = 64
                    row["worker_pid"] = 1064
            sibling_freeze = copy.deepcopy(freeze)
            sibling_freeze["cpu_affinity"] = list(range(7)) + [64]
            with self.assertRaisesRegex(
                    subject.NativeEvidenceError, "unique physical cores"):
                validate(sibling_records, sibling_freeze, sibling_root)

            controller_freeze = copy.deepcopy(freeze)
            controller_freeze["host_identity"]["controller_cpu"] = 64
            with self.assertRaisesRegex(
                    subject.NativeEvidenceError, "controller shares"):
                validate(records, controller_freeze, sibling_root)

            sampler_case = root / "sampler-physical-core"
            sampler_case.mkdir()
            sampler_path = sampler_fixture(sampler_case)
            sampler = json.loads(sampler_path.read_text(encoding="utf-8"))
            sampler["cpu"] = 64
            with self.assertRaisesRegex(
                    subject.NativeEvidenceError, "worker physical core"):
                subject._thermal_window(
                    sampler, 2000000000, 2500000000, WORKER_CPUS, False,
                    controller_cpu=8, sysfs_root=sibling_root)

            controller_sampler_root = write_cpu_topology(
                root / "controller-sampler-sysfs",
                {**{cpu: (0, cpu) for cpu in range(9)},
                 64: (0, 8), 127: (0, 127)})
            with self.assertRaisesRegex(
                    subject.NativeEvidenceError, "controller physical core"):
                subject._thermal_window(
                    sampler, 2000000000, 2500000000, WORKER_CPUS, False,
                    controller_cpu=8, sysfs_root=controller_sampler_root)

    def test_timing_cell_scalar_mutations_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            built = self.build_kind(root, "timing")
            records = built.records
            for mutation in (
                    "missing", "wrong", "noncanonical",
                    "receive_cap_missing", "receive_cap_wrong",
                    "receive_cap_noncanonical"):
                case = root / mutation
                case.mkdir()
                changed = copy.deepcopy(records)
                payload = changed[0]["payload"]
                if mutation == "missing":
                    del payload["invocations_per_slot"]
                elif mutation == "wrong":
                    payload["invocations_per_slot"] += 1
                elif mutation == "noncanonical":
                    # bool is deliberately chosen because it is an int
                    # subclass: the authoritative cell validator must use
                    # exact scalar types rather than isinstance(value, int).
                    payload["invocations_per_slot"] = True
                elif mutation == "receive_cap_missing":
                    del payload["receive_overhead_cap"]
                elif mutation == "receive_cap_wrong":
                    payload["receive_overhead_cap"] -= 1
                else:
                    payload["receive_overhead_cap"] = True
                native = case / "mutant-native.jsonl"
                canonical_jsonl(native, changed)
                output = case / "output.jsonl"
                receipt = case / "receipt.json"
                expected_error = {
                    "missing": "native result payload has an unexpected schema",
                    "wrong": (
                        "timing row is outside the frozen development domain"),
                    "noncanonical": (
                        "timing cell key uses a noncanonical scalar type"),
                    "receive_cap_missing": (
                        "native result payload has an unexpected schema"),
                    "receive_cap_wrong": (
                        "timing row is outside the frozen development domain"),
                    "receive_cap_noncanonical": (
                        "timing cell key uses a noncanonical scalar type"),
                }[mutation]
                with self.subTest(mutation=mutation), \
                        self.assertRaisesRegex(
                            subject.NativeEvidenceError, expected_error):
                    subject.assemble_results(
                        self.contract, "timing", "development", built.freeze,
                        built.traces, native, built.sampler, output,
                        receipt, verify_live_sampler=False,
                        **timing_evidence_kwargs(built))
                self.assertFalse(output.exists())
                self.assertFalse(receipt.exists())

    def test_trace_stream_rejects_missing_duplicate_and_bad_identity(self) -> None:
        mutations = ("missing", "duplicate", "identity", "candidate_count")
        for mutation in mutations:
            with self.subTest(mutation=mutation), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                records, _ = trace_records(
                    self.contract, "recovery", "development")
                if mutation == "missing":
                    records.pop()
                elif mutation == "duplicate":
                    records[-1] = copy.deepcopy(records[0])
                elif mutation == "identity":
                    records[0]["cell_sha256"] = digest("wrong cell")
                else:
                    records[0]["candidate_count"] = 1
                native = root / "native.jsonl"
                canonical_jsonl(native, records)
                with self.assertRaises(subject.NativeEvidenceError):
                    subject.assemble_trace_manifest(
                        self.contract, "recovery", "development", native,
                        root / "traces.jsonl")

    def test_native_envelopes_reject_core_ordinal_binary_and_message_drift(
            self) -> None:
        mutations = (
            "core", "ordinal", "binary", "message", "work", "start_ticks",
            "missing",
        )
        for mutation in mutations:
            with self.subTest(mutation=mutation), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                built = self.build_kind(root, "recovery")
                freeze = contract_api.load_freeze_manifest(
                    self.contract, "development", built.freeze, "recovery")
                changed = copy.deepcopy(built.records)
                if mutation == "core":
                    changed[0]["cpu"] = 9
                elif mutation == "ordinal":
                    changed[0]["ordinal"] = 1
                elif mutation == "binary":
                    changed[0]["worker_binary_sha256"] = digest("substitute")
                elif mutation == "message":
                    changed[1]["message_sha256"] = digest("substitute")
                elif mutation == "work":
                    for row in changed:
                        row["work_sha256"] = changed[0]["work_sha256"]
                elif mutation == "start_ticks":
                    changed[0]["worker_process_start_ticks"] = 999999999999
                else:
                    changed.pop()
                with self.assertRaises(subject.NativeEvidenceError):
                    subject._validate_native_records(
                        self.contract, freeze, "recovery", "development",
                        changed)

    def test_authoritative_validator_rejects_descriptor_trace_construction_and_relabel(
            self) -> None:
        mutations = ("descriptor", "trace", "construction", "relabel")
        for mutation in mutations:
            with self.subTest(mutation=mutation), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                built = self.build_kind(root, "timing")
                changed = copy.deepcopy(built.records)
                payload = changed[0]["payload"]
                if mutation == "descriptor":
                    payload["left_arm_descriptor_sha256"] = digest("wrong")
                    payload["right_arm_descriptor_sha256"] = digest("wrong")
                elif mutation == "trace":
                    payload["trace_sha256"] = digest("wrong")
                elif mutation == "construction":
                    payload["left_construction_attempt"] = 17
                    payload["right_construction_attempt"] = 17
                else:
                    payload["scope"] = "receive_to_success"
                native = root / "mutant-native.jsonl"
                canonical_jsonl(native, changed)
                with self.assertRaises(subject.NativeEvidenceError):
                    subject.assemble_results(
                        self.contract, "timing", "development", built.freeze,
                        built.traces, native, built.sampler,
                        root / "output.jsonl", root / "receipt.json",
                        verify_live_sampler=False,
                        **timing_evidence_kwargs(built))
                self.assertFalse((root / "output.jsonl").exists())
                self.assertFalse((root / "receipt.json").exists())

    def test_thermal_and_edac_mutations_fail_before_publication(self) -> None:
        mutations = (
            "hot", "edac", "short", "overlap", "inode", "late_process",
        )
        for mutation in mutations:
            with self.subTest(mutation=mutation), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                sampler_path = sampler_fixture(
                    root, edac=1 if mutation == "edac" else 0,
                    cpu_c=121.0 if mutation == "hot" else 61.0)
                sampler = json.loads(sampler_path.read_text(encoding="utf-8"))
                if mutation == "short":
                    sampler["window_end_monotonic_ns"] = 2000000000
                elif mutation == "overlap":
                    sampler["cpu"] = 0
                elif mutation == "inode":
                    sampler["csv_inode"] += 1
                elif mutation == "late_process":
                    sampler["process_start_ticks"] = 999999999999
                canonical_file(root / "changed.json", sampler)
                with self.assertRaises(subject.NativeEvidenceError):
                    subject._thermal_window(
                        sampler, 2000000000, 3000000000, WORKER_CPUS, False)

    def test_live_sampler_rejects_python_c_decoy(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            script_path = root / "sampler.py"
            script_path.write_text("# never executed\n", encoding="ascii")
            csv_path = root / "thermal.csv"
            csv_path.write_text("unused\n", encoding="ascii")
            cpu = min(os.sched_getaffinity(0))
            process = subprocess.Popen([
                sys.executable, "-c", "import time; time.sleep(30)",
                str(script_path), "--csv", str(csv_path),
            ])
            try:
                os.sched_setaffinity(process.pid, {cpu})
                start_ticks = subject._parse_proc_start_ticks(process.pid)
                with self.assertRaises(subject.NativeEvidenceError):
                    subject._verify_live_sampler_process(
                        process.pid, cpu, start_ticks, script_path, csv_path)
            finally:
                process.terminate()
                process.wait(timeout=10)

    def test_live_sampler_rejects_relative_script_and_csv_arguments(self) \
            -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            launched = root / "launched"
            attested = root / "attested"
            launched.mkdir()
            attested.mkdir()
            script = "import time\ntime.sleep(30)\n"
            for directory in (launched, attested):
                (directory / "sampler.py").write_text(
                    script, encoding="ascii")
                (directory / "thermal.csv").write_text(
                    "unused\n", encoding="ascii")
            cpu = min(os.sched_getaffinity(0))
            process = subprocess.Popen([
                sys.executable, "sampler.py", "--csv", "thermal.csv",
            ], cwd=str(launched))
            original_cwd = Path.cwd()
            try:
                os.sched_setaffinity(process.pid, {cpu})
                start_ticks = subject._parse_proc_start_ticks(process.pid)
                # Before the absolute-argv rule, verifier-relative resolution
                # incorrectly bound these unrelated same-named files.
                os.chdir(str(attested))
                with self.assertRaises(subject.NativeEvidenceError):
                    subject._verify_live_sampler_process(
                        process.pid,
                        cpu,
                        start_ticks,
                        attested / "sampler.py",
                        attested / "thermal.csv")
            finally:
                os.chdir(str(original_cwd))
                process.terminate()
                process.wait(timeout=10)

    def test_freeze_mutation_during_validation_cannot_publish(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            built = self.build_kind(root, "recovery")
            output = root / "result.jsonl"
            receipt = root / "execution.json"
            original_validate = contract_api.validate_ledger

            def mutate_then_validate(*args, **kwargs):
                summary = original_validate(*args, **kwargs)
                value = json.loads(built.freeze.read_text(encoding="utf-8"))
                value["source_git_commit"] = "2" * 40
                canonical_file(built.freeze, value)
                return summary

            with mock.patch.object(
                    contract_api, "validate_ledger",
                    side_effect=mutate_then_validate):
                with self.assertRaises(subject.NativeEvidenceError):
                    subject.assemble_results(
                        self.contract,
                        "recovery",
                        "development",
                        built.freeze,
                        built.traces,
                        built.native,
                        sampler_fixture(root),
                        output,
                        receipt,
                        verify_live_sampler=False)
            self.assertFalse(output.exists())
            self.assertFalse(receipt.exists())

    def test_terminal_reopen_rejects_native_trace_and_sampler_races(self) \
            -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            built = self.build_kind(root, "recovery")
            sampler_path = sampler_fixture(root)
            sampler_value = json.loads(
                sampler_path.read_text(encoding="utf-8"))
            csv_path = Path(sampler_value["csv_path"])
            originals = {
                "native": built.native.read_bytes(),
                "trace": built.traces.read_bytes(),
                "sampler": sampler_path.read_bytes(),
                "csv": csv_path.read_bytes(),
            }

            def mutate(name: str) -> None:
                if name == "native":
                    rows = list(contract_api._parse_canonical_jsonl(
                        built.native, "native race fixture"))
                    canonical_jsonl(built.native, list(reversed(rows)))
                elif name == "trace":
                    rows = list(contract_api._parse_canonical_jsonl(
                        built.traces, "trace race fixture"))
                    rows[0] = dict(rows[0])
                    rows[0]["trace_sha256"] = digest("raced trace")
                    canonical_jsonl(built.traces, rows)
                elif name == "sampler":
                    value = json.loads(
                        sampler_path.read_text(encoding="utf-8"))
                    value["window_end_monotonic_ns"] = 9000000000
                    canonical_file(sampler_path, value)
                else:
                    lines = csv_path.read_text(encoding="ascii").splitlines()
                    fields = lines[5].split(",")
                    fields[4] = "62.0"
                    lines[5] = ",".join(fields)
                    csv_path.write_text(
                        "\n".join(lines) + "\n", encoding="ascii")

            original_validate = contract_api.validate_ledger
            for name in ("native", "trace", "sampler", "csv"):
                output = root / (name + "-result.jsonl")
                receipt = root / (name + "-execution.json")

                def validate_then_mutate(*args, _name=name, **kwargs):
                    summary = original_validate(*args, **kwargs)
                    mutate(_name)
                    return summary

                try:
                    with self.subTest(name=name), mock.patch.object(
                            contract_api, "validate_ledger",
                            side_effect=validate_then_mutate), \
                            self.assertRaises(subject.NativeEvidenceError):
                        subject.assemble_results(
                            self.contract, "recovery", "development",
                            built.freeze, built.traces, built.native,
                            sampler_path, output, receipt,
                            verify_live_sampler=False)
                    self.assertFalse(output.exists())
                    self.assertFalse(receipt.exists())
                finally:
                    built.native.write_bytes(originals["native"])
                    built.traces.write_bytes(originals["trace"])
                    sampler_path.write_bytes(originals["sampler"])
                    csv_path.write_bytes(originals["csv"])

    def test_execution_receipt_mutation_is_detected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            built = self.build_kind(root, "recovery")
            result = root / "result.jsonl"
            receipt = root / "execution.json"
            sampler_path = sampler_fixture(root)
            subject.assemble_results(
                self.contract, "recovery", "development", built.freeze,
                built.traces, built.native, sampler_path, result, receipt,
                verify_live_sampler=False)
            sampler_value = json.loads(
                (root / "sampler.json").read_text(encoding="utf-8"))
            # The real sampler remains live.  Rows appended after the sealed
            # endpoint must not invalidate the exact prior window.
            with Path(sampler_value["csv_path"]).open(
                    "a", encoding="ascii", newline="") as output:
                output.write(
                    "2026-08-02T00:00:05.000Z,5.0,99.0,4000.0,61.0,"
                    "40.0,40.25,40.5,40.75,41.0,41.25,41.5,41.75,"
                    "0,1.0,1.0,1.0,0,0\n")
            subject.validate_execution_receipt(
                self.contract, "recovery", "development", built.freeze,
                built.traces, built.native, result, receipt,
                verify_live_sampler=False, sampler_path=sampler_path)
            authentic = json.loads(receipt.read_text(encoding="utf-8"))
            alias_mutations = {
                "record-count-float": lambda value: value.__setitem__(
                    "record_count", float(value["record_count"])),
                "worker-start-float": lambda value: value.__setitem__(
                    "worker_start_monotonic_ns",
                    float(value["worker_start_monotonic_ns"])),
                "worker-cpu-bool": lambda value: value["worker_cpus"].
                    __setitem__(0, False),
                "worker-record-cpu-bool": lambda value: value["workers"][0].
                    __setitem__("cpu", False),
                "thermal-sample-count-float": lambda value: value[
                    "thermal"].__setitem__(
                        "sample_count", float(
                            value["thermal"]["sample_count"])),
            }
            for name, mutate in alias_mutations.items():
                with self.subTest(exact_json_alias=name):
                    aliased = copy.deepcopy(authentic)
                    mutate(aliased)
                    aliased["receipt_sha256"] = contract_api.sha256_json({
                        key: value for key, value in aliased.items()
                        if key != "receipt_sha256"
                    })
                    aliased_path = root / (name + "-execution.json")
                    canonical_file(aliased_path, aliased)
                    with self.assertRaises(subject.NativeEvidenceError):
                        subject.validate_execution_receipt(
                            self.contract, "recovery", "development",
                            built.freeze, built.traces, built.native, result,
                            aliased_path, verify_live_sampler=False,
                            sampler_path=sampler_path)
            changed = json.loads(receipt.read_text(encoding="utf-8"))
            changed["thermal"]["edac_ce_max"] = 1
            unsigned = {key: value for key, value in changed.items()
                        if key != "receipt_sha256"}
            changed["receipt_sha256"] = contract_api.sha256_json(unsigned)
            mutant = root / "mutant-execution.json"
            canonical_file(mutant, changed)
            with self.assertRaises(subject.NativeEvidenceError):
                subject.validate_execution_receipt(
                    self.contract, "recovery", "development", built.freeze,
                    built.traces, built.native, result, mutant,
                    verify_live_sampler=False, sampler_path=sampler_path)

            original_receipt = receipt.read_bytes()
            original_validate = contract_api.validate_ledger

            def validate_then_mutate_receipt(*args, **kwargs):
                summary = original_validate(*args, **kwargs)
                value = json.loads(receipt.read_text(encoding="utf-8"))
                # Retain the old self-hash while exploiting Python's
                # 1440 == 1440.0 structural alias.  The terminal comparison
                # must authenticate exact JSON bytes/types, not dict equality.
                value["record_count"] = float(value["record_count"])
                canonical_file(receipt, value)
                return summary

            try:
                with mock.patch.object(
                        contract_api, "validate_ledger",
                        side_effect=validate_then_mutate_receipt), \
                        self.assertRaises(subject.NativeEvidenceError):
                    subject.validate_execution_receipt(
                        self.contract, "recovery", "development",
                        built.freeze, built.traces, built.native, result,
                        receipt, verify_live_sampler=False,
                        sampler_path=sampler_path)
            finally:
                receipt.write_bytes(original_receipt)

    def test_timing_terminal_receipt_rejects_stale_hash_numeric_alias_race(
            self) -> None:
        """Exercise the terminal exact-JSON check on an accepted v1 timing freeze."""
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            built = self.build_kind(root, "timing")
            result = root / "timing-result.jsonl"
            receipt = root / "timing-execution.json"
            evidence = timing_evidence_kwargs(built)
            subject.assemble_results(
                self.contract, "timing", "development", built.freeze,
                built.traces, built.native, built.sampler, result, receipt,
                verify_live_sampler=False, **evidence)
            original_validate = contract_api.validate_timing_receipt

            def validate_then_alias(*args, **kwargs):
                summary = original_validate(*args, **kwargs)
                value = json.loads(receipt.read_text(encoding="utf-8"))
                # Preserve the stale hash: dict equality aliases this float
                # with the authentic integer, canonical JSON equality does not.
                value["record_count"] = float(value["record_count"])
                canonical_file(receipt, value)
                return summary

            with mock.patch.object(
                    contract_api, "validate_timing_receipt",
                    side_effect=validate_then_alias), \
                    self.assertRaisesRegex(
                        subject.NativeEvidenceError,
                        "execution receipt changed"):
                subject.validate_execution_receipt(
                    self.contract, "timing", "development", built.freeze,
                    built.traces, built.native, result, receipt,
                    verify_live_sampler=False, sampler_path=built.sampler,
                    **evidence)

    def test_existing_destination_cannot_leave_a_partial_artifact_pair(
            self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            built = self.build_kind(root, "recovery")
            sampler = sampler_fixture(root)
            result = root / "result.jsonl"
            receipt = root / "execution.json"
            sentinel = b"preexisting\n"

            receipt.write_bytes(sentinel)
            with self.assertRaises(subject.NativeEvidenceError):
                subject.assemble_results(
                    self.contract, "recovery", "development", built.freeze,
                    built.traces, built.native, sampler, result, receipt,
                    verify_live_sampler=False)
            self.assertFalse(result.exists())
            self.assertEqual(receipt.read_bytes(), sentinel)

            receipt.unlink()
            result.write_bytes(sentinel)
            with self.assertRaises(subject.NativeEvidenceError):
                subject.assemble_results(
                    self.contract, "recovery", "development", built.freeze,
                    built.traces, built.native, sampler, result, receipt,
                    verify_live_sampler=False)
            self.assertEqual(result.read_bytes(), sentinel)
            self.assertFalse(receipt.exists())

    def test_publish_new_preserves_visible_commit_after_async_exception(
            self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            staged = root / "staged"
            destination = root / "destination"
            staged.write_bytes(b"new artifact\n")
            real_link = os.link

            class InjectedInterrupt(BaseException):
                pass

            def link_then_interrupt(
                    source: str, target: str, *args: Any,
                    **kwargs: Any) -> None:
                real_link(source, target, *args, **kwargs)
                raise InjectedInterrupt()

            with mock.patch.object(os, "link", side_effect=link_then_interrupt):
                with self.assertRaises(InjectedInterrupt):
                    subject._publish_new(staged, destination)
            self.assertTrue(staged.exists())
            self.assertEqual(destination.read_bytes(), b"new artifact\n")

            destination.write_bytes(b"preexisting winner\n")
            staged = root / "second-staged"
            staged.write_bytes(b"second artifact\n")
            with self.assertRaises(subject.NativeEvidenceError):
                subject._publish_new(staged, destination)
            self.assertEqual(destination.read_bytes(), b"preexisting winner\n")
            self.assertEqual(staged.read_bytes(), b"second artifact\n")

    def test_native_publishers_preserve_commits_after_successful_link(
            self) -> None:
        class InjectedInterrupt(BaseException):
            pass

        for fail_on in (1, 2, 3):
            with self.subTest(publisher="qualification", fail_on=fail_on), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                native = root / "qualification-native.jsonl"
                audit = root / "qualification-audit.jsonl"
                map_path = root / "qualification-map.json"
                execution_receipt = root / "qualification-execution.json"
                canonical_jsonl(native, qualification_records(self.contract))
                real_publish = subject._publish_new
                calls = 0

                def publish_then_interrupt(
                        staged, destination, expected_identity=None):
                    nonlocal calls
                    calls += 1
                    real_publish(staged, destination, expected_identity)
                    if calls == fail_on:
                        raise InjectedInterrupt()

                with mock.patch.object(
                        subject, "_publish_new",
                        side_effect=publish_then_interrupt):
                    with self.assertRaises(InjectedInterrupt):
                        subject.assemble_timing_qualification(
                            self.contract, "development", native, audit,
                            map_path, execution_receipt, SOURCE_COMMIT,
                            qualification_controls(), QUALIFICATION_CPUS,
                            verify_live_workers=False)
                self.assertEqual(audit.exists(), fail_on >= 1)
                self.assertEqual(map_path.exists(), fail_on >= 2)
                self.assertEqual(execution_receipt.exists(), fail_on >= 3)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            qualification = build_qualification(root, self.contract).value
            trace_path = root / "timing-traces.jsonl"
            real_publish = subject._publish_new

            def trace_publish_then_interrupt(
                    staged, destination, expected_identity=None):
                real_publish(staged, destination, expected_identity)
                raise InjectedInterrupt()

            with mock.patch.object(
                    subject, "_publish_new",
                    side_effect=trace_publish_then_interrupt):
                with self.assertRaises(InjectedInterrupt):
                    subject.publish_timing_trace_manifest(
                        self.contract, "development", qualification,
                        trace_path)
            self.assertTrue(trace_path.exists())

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            qualification = build_qualification(root, self.contract).value
            trace_path = root / "replacement-timing-traces.jsonl"
            real_publish = subject._publish_new

            def replace_then_interrupt(
                    staged, destination, expected_identity=None):
                real_publish(staged, destination, expected_identity)
                destination.unlink()
                destination.write_bytes(b"replacement winner\n")
                raise InjectedInterrupt()

            with mock.patch.object(
                    subject, "_publish_new",
                    side_effect=replace_then_interrupt):
                with self.assertRaises(InjectedInterrupt):
                    subject.publish_timing_trace_manifest(
                        self.contract, "development", qualification,
                        trace_path)
            self.assertEqual(trace_path.read_bytes(), b"replacement winner\n")

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = sampler_fixture(
                root, prefix="source-", window_start_ns=1000000000,
                window_end_ns=2000000000)
            source_value = json.loads(source.read_text(encoding="utf-8"))
            output = root / "published-sampler.json"
            real_publish = subject._publish_new

            def sampler_publish_then_interrupt(
                    staged, destination, expected_identity=None):
                real_publish(staged, destination, expected_identity)
                raise InjectedInterrupt()

            with mock.patch.object(
                    subject, "_parse_proc_start_ticks", return_value=1), \
                    mock.patch.object(
                        subject, "_verify_live_sampler_process"), \
                    mock.patch.object(
                        subject, "_publish_new",
                        side_effect=sampler_publish_then_interrupt):
                with self.assertRaises(InjectedInterrupt):
                    subject.write_sampler_attestation(
                        4242, 127, Path(source_value["script_path"]),
                        Path(source_value["csv_path"]), 1000000000,
                        2000000000, output)
            self.assertTrue(output.exists())

    def test_native_publishers_never_remove_preexisting_destinations(
            self) -> None:
        sentinel = b"preexisting winner\n"
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            native = root / "qualification-native.jsonl"
            audit = root / "qualification-audit.jsonl"
            map_path = root / "qualification-map.json"
            execution_receipt = root / "qualification-execution.json"
            canonical_jsonl(native, qualification_records(self.contract))
            audit.write_bytes(sentinel)
            with self.assertRaises(subject.NativeEvidenceError):
                subject.assemble_timing_qualification(
                    self.contract, "development", native, audit, map_path,
                    execution_receipt, SOURCE_COMMIT, qualification_controls(),
                    QUALIFICATION_CPUS, verify_live_workers=False)
            self.assertEqual(audit.read_bytes(), sentinel)
            self.assertFalse(map_path.exists())
            self.assertFalse(execution_receipt.exists())

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            qualification = build_qualification(root, self.contract).value
            trace_path = root / "timing-traces.jsonl"
            trace_path.write_bytes(sentinel)
            with self.assertRaises(subject.NativeEvidenceError):
                subject.publish_timing_trace_manifest(
                    self.contract, "development", qualification, trace_path)
            self.assertEqual(trace_path.read_bytes(), sentinel)

    def test_result_receipt_publication_is_receipt_last(self) -> None:
        for failure_point in ("result", "receipt"):
            with self.subTest(failure_point=failure_point), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                built = self.build_kind(root, "recovery")
                sampler = sampler_fixture(root)
                result = root / "result.jsonl"
                receipt = root / "execution.json"
                target = result if failure_point == "result" else receipt
                original_publish = subject._publish_new

                def publish_then_fail(
                        staged, destination, expected_identity=None):
                    original_publish(staged, destination, expected_identity)
                    if destination == target:
                        raise OSError("injected post-link failure")

                with mock.patch.object(
                        subject, "_publish_new",
                        side_effect=publish_then_fail):
                    with self.assertRaises(OSError):
                        subject.assemble_results(
                            self.contract,
                            "recovery",
                            "development",
                            built.freeze,
                            built.traces,
                            built.native,
                            sampler,
                            result,
                            receipt,
                            verify_live_sampler=False)
                self.assertTrue(result.exists())
                self.assertEqual(receipt.exists(), failure_point == "receipt")

    def test_published_result_validation_failure_never_commits_receipt(
            self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            built = self.build_kind(root, "recovery")
            sampler = sampler_fixture(root)
            result = root / "result.jsonl"
            receipt = root / "execution.json"
            original_validate = contract_api.validate_ledger
            calls = 0

            def reject_published_dependency(*args, **kwargs):
                nonlocal calls
                calls += 1
                if calls == 2:
                    raise contract_api.ContractError(
                        "injected published-result rejection")
                return original_validate(*args, **kwargs)

            with mock.patch.object(
                    contract_api, "validate_ledger",
                    side_effect=reject_published_dependency), \
                    self.assertRaisesRegex(
                        subject.NativeEvidenceError,
                        "injected published-result rejection"):
                subject.assemble_results(
                    self.contract, "recovery", "development", built.freeze,
                    built.traces, built.native, sampler, result, receipt,
                    verify_live_sampler=False)
            self.assertEqual(calls, 2)
            self.assertTrue(result.exists())
            self.assertFalse(receipt.exists())

    def test_result_mutation_during_revalidation_is_detected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            built = self.build_kind(root, "recovery")
            result = root / "result.jsonl"
            receipt = root / "execution.json"
            sampler_path = sampler_fixture(root)
            subject.assemble_results(
                self.contract,
                "recovery",
                "development",
                built.freeze,
                built.traces,
                built.native,
                sampler_path,
                result,
                receipt,
                verify_live_sampler=False)
            original_validate = contract_api.validate_ledger

            def reorder_then_validate(*args, **kwargs):
                summary = original_validate(*args, **kwargs)
                rows = list(contract_api._parse_canonical_jsonl(
                    result, "mutation fixture"))
                canonical_jsonl(result, list(reversed(rows)))
                return summary

            with mock.patch.object(
                    contract_api, "validate_ledger",
                    side_effect=reorder_then_validate):
                with self.assertRaises(subject.NativeEvidenceError):
                    subject.validate_execution_receipt(
                        self.contract,
                        "recovery",
                        "development",
                        built.freeze,
                        built.traces,
                        built.native,
                        result,
                        receipt,
                        verify_live_sampler=False,
                        sampler_path=sampler_path)

    def test_thermal_window_ignores_only_prewindow_history(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampler_path = sampler_fixture(root)
            sampler = json.loads(sampler_path.read_text(encoding="utf-8"))
            csv_path = Path(sampler["csv_path"])
            lines = csv_path.read_text(encoding="ascii").splitlines(True)
            csv_path.write_text(
                lines[0] + "malformed historical row\n" + "".join(lines[1:]),
                encoding="ascii")
            summary = subject._thermal_window(
                sampler, 2000000000, 3000000000, WORKER_CPUS, False)
            self.assertEqual(summary["sample_count"], 10)

            lines = csv_path.read_text(encoding="ascii").splitlines(True)
            # Header, ignored historical row, first valid row, then a row
            # inside the selected interval: the latter must fail closed.
            lines[3] = "malformed in-window row\n"
            csv_path.write_text("".join(lines), encoding="ascii")
            with self.assertRaises(subject.NativeEvidenceError):
                subject._thermal_window(
                    sampler, 2000000000, 3000000000, WORKER_CPUS, False)

    def test_thermal_fdopen_failure_closes_duplicated_descriptor(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampler_path = sampler_fixture(root)
            sampler = json.loads(sampler_path.read_text(encoding="utf-8"))
            descriptors = []

            def reject_fdopen(descriptor: int, _mode: str):
                descriptors.append(descriptor)
                raise OSError("injected thermal fdopen failure")

            with mock.patch.object(
                    subject.os, "fdopen", side_effect=reject_fdopen), \
                    self.assertRaisesRegex(
                        subject.NativeEvidenceError,
                        "injected thermal fdopen failure"):
                subject._thermal_window(
                    sampler, 2000000000, 3000000000, WORKER_CPUS, False)
            self.assertEqual(len(descriptors), 1)
            with self.assertRaises(OSError) as raised:
                os.fstat(descriptors[0])
            self.assertEqual(raised.exception.errno, errno.EBADF)

    def test_thermal_endpoint_must_be_newline_terminated(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampler_path = sampler_fixture(root)
            sampler = json.loads(sampler_path.read_text(encoding="utf-8"))
            csv_path = Path(sampler["csv_path"])
            contents = csv_path.read_bytes()
            self.assertTrue(contents.endswith(b"\n"))
            csv_path.write_bytes(contents[:-1])
            with self.assertRaises(subject.NativeEvidenceError):
                subject._thermal_window(
                    sampler, 2000000000, 3000000000, WORKER_CPUS, False)


if __name__ == "__main__":
    unittest.main()
