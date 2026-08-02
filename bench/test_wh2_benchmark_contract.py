#!/usr/bin/env python3
"""Fail-closed tests for the minimal Wirehair2 benchmark contract."""

from __future__ import annotations

import copy
import hashlib
import io
import json
import math
from pathlib import Path
import sys
import tempfile
from typing import Any, Callable, Dict, List, Mapping, Optional, Sequence, Tuple
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_benchmark_contract as subject


EXPECTED_RECOVERY_HASHES = {
    "development": "f97f28c211428cd77aed97160073b192d93014cb4a61a844bc7d76375ac61b77",
    "final_raw": "c33160c0379675dac9de3162463e3100041bc8eeb55c002dcf6c3473973e7749",
    "final_repaired": "e247b0ff51f47383f415b1747a783e0c4b4a91e97862e9e72aa660241ff0aeee",
    "final_validation": "a63611ad300396a1b8977281fcd387fde4bd207a3ca82043c5e234e8fc39d4de",
    "cross_width_validation": "80fb4c19b31f88101ef4b3c480e6d5ff4a3bba7d45c96f02dc3518cd7c1c0399",
}
EXPECTED_TIMING_BASE_HASHES = {
    "development": "eab25fe96642d8dd12d4b64e91fe00b533d464f1bd7487c84e670ce888e2d164",
    "final": "5b53ef3780f9b5b3cab83b411da0c8d0d06abaa3c76a9dc51d96425cd9de6c94",
}
EXPECTED_CONTRACT_SHA256 = \
    "64aa60eea7a13b143420349d878b1ad334d314eb6f996ed6018833a73ef89e3e"


def digest(label: str) -> str:
    return hashlib.sha256(label.encode("ascii")).hexdigest()


Outcome = Tuple[str, Any]
OutcomeFunction = Callable[[str, int, Mapping[str, Any]], Outcome]


def successful(_arm: str, _ordinal: int, _cell: Mapping[str, Any]) -> Outcome:
    return "success", 0


def ledger_rows(
        contract: Mapping[str, Any],
        arms: Sequence[str] = ("wirehair2_head", "wirehair1", "candidate"),
        outcome: OutcomeFunction = successful,
        trace_hashes: Sequence[str] = (),
        phase: str = "development",
        repair_attempts: Optional[Mapping[str, Mapping[int, int]]] = None,
        repair_map_hashes: Optional[Mapping[str, str]] = None,
) -> List[Dict[str, Any]]:
    if repair_attempts is None:
        repair_attempts = {}
    if repair_map_hashes is None:
        repair_map_hashes = {}
    rows: List[Dict[str, Any]] = []
    for ordinal, cell in enumerate(subject.iter_recovery_cells(contract, phase)):
        trace = trace_hashes[ordinal] if trace_hashes else digest(
            "trace:" + subject.sha256_json(cell))
        for arm in arms:
            result, extra = outcome(arm, ordinal, cell)
            descriptor = \
                "739092a7824449e6168f08b46661dfbe8ad5495ea4166b36073c79cd3bacdd11" \
                if arm == RAW_TEST_ARM else digest("descriptor:" + arm)
            codec = "wirehair2_certified" if arm == "wirehair2_head" else \
                "wirehair1" if arm == "wirehair1" else "wirehair2_experiment"
            attempt = repair_attempts.get(arm, {}).get(
                cell["K"],
                0 if codec == "wirehair1" else cell["base_seed_attempt"])
            rows.append({
                "arm": arm,
                **cell,
                "outcome": result,
                "decoded_extra": extra,
                "cell_sha256": subject.sha256_json(cell),
                "trace_sha256": trace,
                "binary_sha256": digest("binary:" + arm),
                "arm_descriptor_sha256": descriptor,
                "construction_attempt": attempt,
                "realized_construction_sha256":
                    subject.generic_realized_construction_sha256(
                        codec, descriptor, cell["K"], cell["block_bytes"],
                        attempt),
                "repair_map_sha256": repair_map_hashes.get(arm, "0" * 64),
            })
    return rows


def write_ledger(root: Path, rows: Sequence[Mapping[str, Any]],
                 *, terminal_newline: bool = True) -> Path:
    path = root / "ledger.jsonl"
    data = "\n".join(subject.canonical_json(row) for row in rows)
    if terminal_newline:
        data += "\n"
    path.write_bytes(data.encode("utf-8"))
    return path


def trace_rows(contract: Mapping[str, Any],
               trace_hashes: Sequence[str] = (),
               phase: str = "development") -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for ordinal, cell in enumerate(subject.iter_recovery_cells(
            contract, phase)):
        trace = trace_hashes[ordinal] if trace_hashes else digest(
            "trace:" + subject.sha256_json(cell))
        rows.append({
            "ordinal": ordinal,
            "cell_sha256": subject.sha256_json(cell),
            "trace_sha256": trace,
        })
    return rows


def write_trace_manifest(
        root: Path, contract: Mapping[str, Any],
        trace_hashes: Sequence[str] = (),
        phase: str = "development") -> Path:
    path = root / "traces.jsonl"
    data = "".join(
        subject.canonical_json(row) + "\n"
        for row in trace_rows(contract, trace_hashes, phase))
    path.write_bytes(data.encode("utf-8"))
    return path


def write_freeze_manifest(
        root: Path, contract: Mapping[str, Any], trace_path: Path,
        arms: Sequence[str], phase: str = "development",
        repair_map_hashes: Optional[Mapping[str, str]] = None,
        training_trace_sha256: str = "0" * 64) -> Path:
    if repair_map_hashes is None:
        repair_map_hashes = {}
    arm_records = []
    for arm in arms:
        codec = "wirehair2_certified" if arm == "wirehair2_head" else \
            "wirehair1" if arm == "wirehair1" else "wirehair2_experiment"
        mapped = arm in repair_map_hashes
        arm_records.append({
            "arm": arm,
            "codec": codec,
            "binary_sha256": digest("binary:" + arm),
            "arm_descriptor_sha256":
                "739092a7824449e6168f08b46661dfbe8ad5495ea4166b36073c79cd3bacdd11"
                if arm == RAW_TEST_ARM else digest("descriptor:" + arm),
            "construction_policy":
                "not_applicable" if codec == "wirehair1" else
                "repair_map" if mapped else "raw_base",
            "repair_map_sha256": repair_map_hashes.get(arm, "0" * 64),
        })
    manifest = {
        "schema": subject.FREEZE_SCHEMA,
        "contract_sha256": subject.contract_sha256(contract),
        "evidence_kind": "recovery",
        "phase": phase,
        "domain_sha256": contract["recovery"]["domains"][phase][
            "domain_sha256"],
        "source_git_commit": "1" * 40,
        "arm_roster": list(arms),
        "arm_roster_sha256": subject.arm_roster_sha256(arms),
        "trace_manifest_sha256": subject.trace_manifest_sha256(
            contract, "recovery", phase, trace_path),
        "repair_training_trace_manifest_sha256": training_trace_sha256,
        "commands": [["wirehair_v2_bench", "raw-cell", "--frozen"]],
        "cpu_affinity": [0, 1],
        "host_identity": {"name": "unit-test-host"},
        "arms": arm_records,
    }
    path = root / "freeze.json"
    path.write_bytes((subject.canonical_json(manifest) + "\n").encode("utf-8"))
    return path


RAW_TEST_ARM = "wirehair2_raw_d12_h12_periodic"


def raw_construction_fields(attempt: int, K: int = 2) -> Dict[str, Any]:
    return {
        "construction_seed_basis": subject.RAW_CONSTRUCTION_SEED_BASIS,
        "seed_schedule_sha256": subject.RAW_SEED_SCHEDULE_SHA256,
        "precode_attempt": attempt,
        "packet_attempt": attempt,
        "effective_precode_seed":
            subject._effective_raw_precode_seed(attempt),
        "effective_packet_seed":
            subject._effective_raw_packet_seed(attempt),
        "staircase": subject._raw_staircase_for_K(K),
        "binary_dense_rows": 12,
        "gf256_heavy_rows": 12,
        "source_hits": 3 if K >= 10000 else 2,
        "dense_identity_corner": False,
        "heavy_family": "periodic-cauchy",
        "mix_count": 3,
    }


def raw_ledger_rows(
        contract: Mapping[str, Any], trace_hashes: Sequence[str] = (),
        phase: str = "development") -> List[Dict[str, Any]]:
    arms = ("wirehair2_head", "wirehair1", RAW_TEST_ARM)
    rows = ledger_rows(
        contract, arms=arms, trace_hashes=trace_hashes, phase=phase)
    for row in rows:
        if row["arm"] != RAW_TEST_ARM:
            continue
        fields = raw_construction_fields(
            row["construction_attempt"], row["K"])
        row.update(fields)
        row["realized_construction_sha256"] = \
            subject.raw_realized_construction_sha256(
                "wirehair2_experiment", RAW_TEST_ARM,
                row["arm_descriptor_sha256"], row["K"],
                row["block_bytes"], fields)
    return rows


def write_raw_freeze_manifest(
        root: Path, contract: Mapping[str, Any], trace_path: Path,
        phase: str = "development") -> Path:
    arms = ("wirehair2_head", "wirehair1", RAW_TEST_ARM)
    arm_records = []
    for arm in arms:
        codec = arm_codec(arm)
        if arm == "wirehair2_head":
            basis = subject.PRODUCTION_CONSTRUCTION_SEED_BASIS
            schedule = "0" * 64
        elif arm == "wirehair1":
            basis = subject.NOT_APPLICABLE_CONSTRUCTION_SEED_BASIS
            schedule = "0" * 64
        else:
            basis = subject.RAW_CONSTRUCTION_SEED_BASIS
            schedule = subject.RAW_SEED_SCHEDULE_SHA256
        arm_records.append({
            "arm": arm,
            "codec": codec,
            "binary_sha256": digest("binary:" + arm),
            "arm_descriptor_sha256":
                "739092a7824449e6168f08b46661dfbe8ad5495ea4166b36073c79cd3bacdd11"
                if arm == RAW_TEST_ARM else digest("descriptor:" + arm),
            "construction_policy":
                "not_applicable" if codec == "wirehair1" else "raw_base",
            "repair_map_sha256": "0" * 64,
            "construction_seed_basis": basis,
            "seed_schedule_sha256": schedule,
        })
    manifest = {
        "schema": subject.RAW_FREEZE_SCHEMA,
        "contract_sha256": subject.contract_sha256(contract),
        "evidence_kind": "recovery",
        "phase": phase,
        "domain_sha256": contract["recovery"]["domains"][phase][
            "domain_sha256"],
        "source_git_commit": "1" * 40,
        "arm_roster": list(arms),
        "arm_roster_sha256": subject.arm_roster_sha256(arms),
        "trace_manifest_sha256": subject.trace_manifest_sha256(
            contract, "recovery", phase, trace_path),
        "repair_training_trace_manifest_sha256": "0" * 64,
        "commands": [["wirehair_v2_bench", "raw-cell", "--frozen"]],
        "cpu_affinity": [0, 1],
        "host_identity": {"name": "unit-test-host"},
        "arms": arm_records,
    }
    path = root / "raw-freeze.json"
    path.write_bytes((subject.canonical_json(manifest) + "\n").encode("utf-8"))
    return path


def write_timing_qualification(
        root: Path, contract: Mapping[str, Any], phase: str = "development",
        retry_offsets: Optional[Sequence[int]] = None,
        ) -> Tuple[subject.TimingQualification, Path, Path, str, List[str]]:
    base_cells = list(subject.iter_timing_base_cells(contract, phase))
    if retry_offsets is None:
        retry_offsets = [0] * len(base_cells)
    offsets = list(retry_offsets)
    if len(offsets) != len(base_cells):
        raise AssertionError("test qualification offset count")
    audit_rows = []
    selected_traces = []
    for ordinal, (base_cell, selected) in enumerate(zip(base_cells, offsets)):
        for retry in range(selected + 1):
            terminal = retry == selected
            trace = digest(
                "timing-qualified-trace:{}:{}:{}".format(
                    phase, ordinal, retry))
            audit_rows.append({
                "ordinal": ordinal,
                "base_cell_sha256": subject.sha256_json(base_cell),
                "loss_retry_offset": retry,
                "loss_seed": subject._qualified_timing_loss_seed(
                    base_cell["base_loss_seed"], retry),
                "trace_sha256": trace,
                "wirehair2_head_outcome":
                    "success" if terminal else "need_more_at_bound",
                "wirehair2_head_decoded_extra": 4 if terminal else None,
                "wirehair1_outcome": "success",
                "wirehair1_decoded_extra": 0,
            })
            if terminal:
                selected_traces.append(trace)
    audit_path = root / ("timing-qualification-" + phase + ".jsonl")
    audit_path.write_bytes("".join(
        subject.canonical_json(row) + "\n" for row in audit_rows
    ).encode("utf-8"))
    controls = []
    for arm, scope in (
            ("wirehair2_head", "decoder_solve"),
            ("wirehair1", "receive_to_success")):
        mapped = phase == "final" and arm == "wirehair2_head"
        controls.append({
            "arm": arm,
            "scope": scope,
            "binary_sha256": digest("binary:" + arm),
            "arm_descriptor_sha256": digest("descriptor:" + arm),
            "construction_policy":
                "not_applicable" if arm == "wirehair1" else
                "repair_map" if mapped else "raw_base",
            "repair_map_sha256":
                digest("production map:" + arm) if mapped else "0" * 64,
        })
    value = {
        "schema": subject.TIMING_QUALIFICATION_MAP_SCHEMA,
        "contract_sha256": subject.contract_sha256(contract),
        "phase": phase,
        "source_git_commit": "1" * 40,
        "base_domain_sha256": contract["timing"]["domains"][phase][
            "base_domain_sha256"],
        "qualified_domain_sha256":
            subject._timing_domain_sha256_from_offsets(
                contract, phase, offsets),
        "entry_kind": subject.TIMING_QUALIFICATION_ENTRY_KIND,
        "controls": controls,
        "qualification_audit_sha256":
            subject.timing_qualification_audit_sha256(
                contract, phase, audit_path),
        "selected_trace_roster_sha256":
            subject.timing_selected_trace_roster_sha256(selected_traces),
        "retry_offsets": offsets,
    }
    map_path = root / ("timing-qualification-" + phase + ".json")
    map_path.write_bytes(
        (subject.canonical_json(value) + "\n").encode("utf-8"))
    map_sha256 = subject.timing_qualification_map_sha256(value)
    qualification = subject.load_timing_qualification_map(
        contract, phase, map_path, audit_path, map_sha256)
    return qualification, map_path, audit_path, map_sha256, selected_traces


def timing_trace_hashes(
        contract: Mapping[str, Any],
        qualification: subject.TimingQualification) -> List[str]:
    if qualification.phase != "development":
        raise AssertionError("test helper expects development qualification")
    return list(qualification.selected_trace_sha256s)


def write_timing_trace_manifest(
        root: Path, contract: Mapping[str, Any],
        hashes: Sequence[str],
        qualification: subject.TimingQualification) -> Path:
    rows = []
    for ordinal, cell in enumerate(subject.iter_timing_cells(
            contract, qualification.phase, qualification)):
        rows.append({
            "ordinal": ordinal,
            "cell_sha256": subject.sha256_json(cell),
            "trace_sha256": hashes[ordinal],
        })
    path = root / "timing-traces.jsonl"
    path.write_bytes("".join(
        subject.canonical_json(row) + "\n" for row in rows).encode("utf-8"))
    return path


def arm_codec(arm: str) -> str:
    if arm == "wirehair2_head":
        return "wirehair2_certified"
    if arm == "wirehair1":
        return "wirehair1"
    return "wirehair2_experiment"


def write_timing_freeze_manifest(
        root: Path, contract: Mapping[str, Any], trace_path: Path,
        arms: Sequence[str],
        qualification: subject.TimingQualification) -> Path:
    arm_records = []
    for arm in arms:
        codec = arm_codec(arm)
        arm_records.append({
            "arm": arm,
            "codec": codec,
            "binary_sha256": digest("binary:" + arm),
            "arm_descriptor_sha256": digest("descriptor:" + arm),
            "construction_policy":
                "not_applicable" if codec == "wirehair1" else "raw_base",
            "repair_map_sha256": "0" * 64,
        })
    manifest = {
        "schema": subject.FREEZE_SCHEMA,
        "contract_sha256": subject.contract_sha256(contract),
        "evidence_kind": "timing",
        "phase": "development",
        "domain_sha256": qualification.qualified_domain_sha256,
        "source_git_commit": "1" * 40,
        "arm_roster": list(arms),
        "arm_roster_sha256": subject.arm_roster_sha256(arms),
        "trace_manifest_sha256": subject.trace_manifest_sha256(
            contract, "timing", "development", trace_path, qualification),
        "repair_training_trace_manifest_sha256": "0" * 64,
        "commands": [["wirehair_v2_bench", "timing-cell", "--frozen"]],
        "cpu_affinity": [0, 1],
        "host_identity": {"name": "unit-test-host"},
        "arms": arm_records,
    }
    path = root / "timing-freeze.json"
    path.write_bytes((subject.canonical_json(manifest) + "\n").encode("utf-8"))
    return path


def write_final_freeze_manifests(
        root: Path, contract: Mapping[str, Any],
        timing_qualification: subject.TimingQualification,
        arms: Sequence[str] = ("wirehair2_head", "wirehair1", "candidate"),
        ) -> Dict[Tuple[str, str], Path]:
    training_trace = digest("final repaired training traces")
    map_hashes = {
        arm: digest("production map:" + arm)
        for arm in arms if arm != "wirehair1"
    }
    timing_trace_path = write_timing_trace_manifest(
        root, contract, timing_qualification.selected_trace_sha256s,
        timing_qualification)
    timing_trace_sha256 = subject.trace_manifest_sha256(
        contract, "timing", "final", timing_trace_path,
        timing_qualification)
    paths: Dict[Tuple[str, str], Path] = {}
    keys = (
        ("recovery", "final_raw"),
        ("recovery", "final_repaired"),
        ("recovery", "final_validation"),
        ("recovery", "cross_width_validation"),
        ("timing", "final"),
    )
    for evidence_kind, phase in keys:
        raw = phase == "final_raw"
        records = []
        for arm in arms:
            codec = arm_codec(arm)
            records.append({
                "arm": arm,
                "codec": codec,
                "binary_sha256": digest("binary:" + arm),
                "arm_descriptor_sha256": digest("descriptor:" + arm),
                "construction_policy":
                    "not_applicable" if codec == "wirehair1" else
                    "raw_base" if raw else "repair_map",
                "repair_map_sha256":
                    "0" * 64 if raw or codec == "wirehair1" else
                    map_hashes[arm],
            })
        trace_sha256 = timing_trace_sha256 if evidence_kind == "timing" else \
            training_trace if phase == "final_repaired" else \
            digest("trace:{}:{}".format(evidence_kind, phase))
        manifest = {
            "schema": subject.FREEZE_SCHEMA,
            "contract_sha256": subject.contract_sha256(contract),
            "evidence_kind": evidence_kind,
            "phase": phase,
            "domain_sha256":
                timing_qualification.qualified_domain_sha256
                if evidence_kind == "timing" else
                contract[evidence_kind]["domains"][phase]["domain_sha256"],
            "source_git_commit": "1" * 40,
            "arm_roster": list(arms),
            "arm_roster_sha256": subject.arm_roster_sha256(arms),
            "trace_manifest_sha256": trace_sha256,
            "repair_training_trace_manifest_sha256":
                "0" * 64 if raw else training_trace,
            "commands": [["wirehair_v2_bench", evidence_kind, phase]],
            "cpu_affinity": [0, 1],
            "host_identity": {"name": "unit-test-host"},
            "arms": records,
        }
        path = root / (evidence_kind + "-" + phase + "-freeze.json")
        path.write_bytes(
            (subject.canonical_json(manifest) + "\n").encode("utf-8"))
        paths[(evidence_kind, phase)] = path
    return paths


def architecture_selection_receipt(
        contract: Mapping[str, Any],
        timing_qualification: subject.TimingQualification,
        selected: str = "candidate",
        ) -> Dict[str, Any]:
    selected_artifact = {
        "codec": "wirehair2_experiment",
        "arm_descriptor_sha256": digest("descriptor:" + selected),
    }
    result = {
        "schema": subject.SCHEMA + ".architecture-selection.v1",
        "contract_sha256": subject.contract_sha256(contract),
        "recovery_domain_sha256": contract["recovery"]["domains"]
            ["development"]["domain_sha256"],
        "timing_base_domain_sha256": contract["timing"]["domains"]
            ["development"]["base_domain_sha256"],
        "timing_domain_sha256":
            timing_qualification.qualified_domain_sha256,
        "timing_qualification_map_sha256": timing_qualification.map_sha256,
        "recovery_freeze_manifest_sha256": digest("recovery freeze"),
        "timing_freeze_manifest_sha256": digest("timing freeze"),
        "architecture_artifact_sha256": digest("architecture artifacts"),
        "recovery_cells_per_arm": 360,
        "timing_rows": 25344,
        "candidate_roster": [selected],
        "eligible_candidates": [selected],
        "eligible_overhead0_failures": {selected: 0},
        "minimum_overhead0_failures": 0,
        "recovery_equivalence_allowance": 1,
        "recovery_equivalent_candidates": [selected],
        "ranking": [{
            "arm": selected,
            "decoder_solve_mean_log_ratio": -0.1,
            "failures_overhead0": 0,
            "failures_overhead1": 0,
            "failures_overhead2": 0,
            "failures_overhead4": 0,
        }],
        "selected_arm": selected,
        "selected_codec": selected_artifact["codec"],
        "selected_arm_descriptor_sha256":
            selected_artifact["arm_descriptor_sha256"],
        "selected_architecture_sha256":
            subject.selected_architecture_sha256(selected, selected_artifact),
    }
    result["selection_sha256"] = subject.sha256_json(result)
    return result


def timing_receipt_rows(
        contract: Mapping[str, Any], hashes: Sequence[str],
        qualification: subject.TimingQualification,
        arms: Sequence[str] = ("wirehair2_head", "wirehair1", "candidate"),
        candidate_scale: float = 0.8) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    panels = subject.timing_panels(contract, arms)
    for ordinal, cell in enumerate(subject.iter_timing_cells(
            contract, "development", qualification)):
        for panel in panels:
            left = panel["left_arm"]
            right = panel["right_arm"]
            order = subject.timing_order(panel, cell["replicate"])
            if panel["panel_kind"] == "AA":
                left_ns = right_ns = 100000
            else:
                left_ns = int(100000 * candidate_scale)
                right_ns = 100000
            primary = [left_ns, right_ns, right_ns, left_ns] \
                if order == "ABBA" else \
                [right_ns, left_ns, left_ns, right_ns]
            opposite = [right_ns, left_ns, left_ns, right_ns] \
                if order == "ABBA" else \
                [left_ns, right_ns, right_ns, left_ns]
            elapsed = primary + opposite
            row: Dict[str, Any] = {
                **cell,
                **panel,
                "order": order,
                "left_outcome": "success",
                "right_outcome": "success",
                "left_decoded_extra": None if panel["scope"] ==
                    "encoder_init_plus_first_K_symbols" else
                    cell["receive_overhead_cap"] if panel["scope"] ==
                    "receive_to_success" else cell["fixed_received_overhead"],
                "right_decoded_extra": None if panel["scope"] ==
                    "encoder_init_plus_first_K_symbols" else
                    cell["receive_overhead_cap"] if panel["scope"] ==
                    "receive_to_success" else cell["fixed_received_overhead"],
                "elapsed_ns": elapsed,
                "cell_sha256": subject.sha256_json(cell),
                "trace_sha256": hashes[ordinal],
            }
            for side, arm in (("left", left), ("right", right)):
                codec = arm_codec(arm)
                descriptor = digest("descriptor:" + arm)
                attempt = 0 if codec == "wirehair1" else \
                    cell["base_seed_attempt"]
                row[side + "_binary_sha256"] = digest("binary:" + arm)
                row[side + "_arm_descriptor_sha256"] = descriptor
                row[side + "_construction_attempt"] = attempt
                row[side + "_realized_construction_sha256"] = \
                    subject.generic_realized_construction_sha256(
                        codec, descriptor, cell["K"], cell["block_bytes"],
                        attempt)
                row[side + "_repair_map_sha256"] = "0" * 64
            rows.append(row)
    return rows


def write_timing_receipt(root: Path,
                         rows: Sequence[Mapping[str, Any]]) -> Path:
    path = root / "timing.jsonl"
    path.write_bytes("".join(
        subject.canonical_json(row) + "\n" for row in rows).encode("utf-8"))
    return path


def write_repair_map(
        root: Path, contract: Mapping[str, Any], arm: str,
        training_trace_sha256: str,
        retry_overrides: Optional[Mapping[int, int]] = None,
        ) -> Tuple[Path, str, Dict[int, int]]:
    if retry_overrides is None:
        retry_overrides = {}
    offsets = [retry_overrides.get(K, 0) for K in range(2, 64001)]
    value = {
        "schema": subject.REPAIR_MAP_SCHEMA,
        "contract_sha256": subject.contract_sha256(contract),
        "training_domain_sha256": contract["recovery"]["domains"][
            "final_repaired"]["domain_sha256"],
        "training_trace_manifest_sha256": training_trace_sha256,
        "arm": arm,
        "source_git_commit": "1" * 40,
        "binary_sha256": digest("binary:" + arm),
        "arm_descriptor_sha256": digest("descriptor:" + arm),
        "production_base_seed_attempt": 0,
        "entry_kind": "retry_offset_indexed_by_K_minus_2",
        "attempt_derivation": subject.EXPECTED_ATTEMPT_DERIVATION,
        "repair_rule": subject.EXPECTED_REPAIR_RULE,
        "retry_offsets": offsets,
    }
    path = root / (arm + "-repair-map.json")
    path.write_bytes((subject.canonical_json(value) + "\n").encode("utf-8"))
    return path, subject.repair_map_sha256(value), {
        K: (contract["seeds"]["production_base_seed_attempt"] + offset) & 0xff
        for K, offset in zip(range(2, 64001), offsets)
    }


class ContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = subject.load_contract()
        cls.qualification_temp = tempfile.TemporaryDirectory()
        root = Path(cls.qualification_temp.name)
        cls.development_qualification = write_timing_qualification(
            root, cls.contract, "development")[0]
        cls.final_qualification = write_timing_qualification(
            root, cls.contract, "final")[0]

    @classmethod
    def tearDownClass(cls) -> None:
        cls.qualification_temp.cleanup()

    def test_v6_timing_identity_and_default_contract_path_are_exact(self) -> None:
        self.assertEqual(subject.SCHEMA,
                         "wirehair.wh2.benchmark-contract.v3")
        self.assertEqual(subject.DEFAULT_CONTRACT.name,
                         "wh2_benchmark_contract_v3.json")
        self.assertEqual(self.contract["schema"], subject.SCHEMA)
        self.assertEqual(self.contract["contract_id"],
                         "wh2-pure-gf256-1pct-v6")
        self.assertEqual(subject.contract_sha256(self.contract),
                         EXPECTED_CONTRACT_SHA256)

    def test_exact_recovery_domains_and_hashes(self) -> None:
        expected_counts = {
            "development": 360,
            "final_raw": 575991,
            "final_repaired": 575991,
            "final_validation": 575991,
            "cross_width_validation": 1440,
        }
        for phase, expected_hash in EXPECTED_RECOVERY_HASHES.items():
            with self.subTest(phase=phase):
                domain = self.contract["recovery"]["domains"][phase]
                self.assertEqual(domain["expected_cells_per_arm"], expected_counts[phase])
                self.assertEqual(domain["domain_sha256"], expected_hash)
                self.assertEqual(subject.recovery_domain_sha256(
                    self.contract, phase), expected_hash)

    def test_exact_timing_domains_and_hashes(self) -> None:
        for phase, expected_hash in EXPECTED_TIMING_BASE_HASHES.items():
            with self.subTest(phase=phase):
                domain = self.contract["timing"]["domains"][phase]
                self.assertEqual(subject.timing_base_domain_sha256(
                    self.contract, phase), expected_hash)
                self.assertEqual(domain["base_domain_sha256"], expected_hash)
                qualification = self.development_qualification \
                    if phase == "development" else self.final_qualification
                self.assertEqual(subject.timing_domain_sha256(
                    self.contract, phase, qualification),
                    qualification.qualified_domain_sha256)
        self.assertEqual(
            self.contract["timing"]["domains"]["development"]["expected_cells"],
            2304,
        )
        self.assertEqual(
            self.contract["timing"]["domains"]["final"]["expected_cells"],
            14400,
        )

    def test_timing_qualification_selects_lowest_retry_and_preserves_source(
            self) -> None:
        offsets = [0] * 2304
        offsets[0] = 2
        with tempfile.TemporaryDirectory() as temporary:
            qualification, _map_path, _audit_path, _map_sha, _traces = \
                write_timing_qualification(
                    Path(temporary), self.contract, "development", offsets)
            cells = list(subject.iter_timing_cells(
                self.contract, "development", qualification))
        base = next(subject.iter_timing_base_cells(
            self.contract, "development"))
        self.assertEqual(cells[0]["loss_retry_offset"], 2)
        self.assertEqual(
            cells[0]["loss_seed"],
            subject._qualified_timing_loss_seed(base["base_loss_seed"], 2))
        self.assertEqual(
            cells[0]["base_cell_sha256"], subject.sha256_json(base))
        self.assertEqual(
            {key: cells[0][key]
             for key in self.contract["timing"]["base_cell_key"]}, base)
        self.assertNotEqual(
            subject._qualified_timing_loss_seed(base["base_loss_seed"], 1),
            cells[0]["loss_seed"])

    def test_timing_qualification_map_and_audit_mutants_fail_closed(
            self) -> None:
        offsets = [0] * 2304
        offsets[0] = 1
        mutations = (
            "nonlowest", "terminal", "seed", "candidate_audit",
            "candidate_map", "candidate_control", "qualified_domain",
            "selected_trace_roster", "base_cell", "gap", "offset_type",
            "source_commit",
        )
        for mutation in mutations:
            with self.subTest(mutation=mutation), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                _qualification, map_path, audit_path, _map_sha, _traces = \
                    write_timing_qualification(
                        root, self.contract, "development", offsets)
                value = json.loads(map_path.read_text(encoding="utf-8"))
                rows = [json.loads(line) for line in
                        audit_path.read_text(encoding="utf-8").splitlines()]
                if mutation == "nonlowest":
                    rows[0]["wirehair2_head_outcome"] = "success"
                    rows[0]["wirehair2_head_decoded_extra"] = 4
                elif mutation == "terminal":
                    rows[1]["wirehair1_outcome"] = "need_more_at_bound"
                    rows[1]["wirehair1_decoded_extra"] = None
                elif mutation == "seed":
                    rows[0]["loss_seed"] = "0x0000000000000000"
                elif mutation == "candidate_audit":
                    rows[0]["candidate_outcome"] = "success"
                elif mutation == "candidate_map":
                    value["candidate_timing_ns"] = 1
                elif mutation == "candidate_control":
                    candidate = copy.deepcopy(value["controls"][0])
                    candidate["arm"] = "candidate"
                    value["controls"].append(candidate)
                elif mutation == "qualified_domain":
                    value["qualified_domain_sha256"] = digest(
                        "forged qualified domain")
                elif mutation == "selected_trace_roster":
                    value["selected_trace_roster_sha256"] = digest(
                        "forged selected trace roster")
                elif mutation == "base_cell":
                    rows[0]["base_cell_sha256"] = digest(
                        "substituted source cell")
                elif mutation == "gap":
                    del rows[0]
                elif mutation == "offset_type":
                    value["retry_offsets"][0] = True
                else:
                    value["source_git_commit"] = "candidate"
                audit_path.write_bytes("".join(
                    subject.canonical_json(row) + "\n" for row in rows
                ).encode("utf-8"))
                value["qualification_audit_sha256"] = \
                    subject.timing_qualification_audit_sha256(
                        self.contract, "development", audit_path)
                map_path.write_bytes(
                    (subject.canonical_json(value) + "\n").encode("utf-8"))
                expected = subject.timing_qualification_map_sha256(value)
                with self.assertRaises(subject.ContractError):
                    subject.load_timing_qualification_map(
                        self.contract, "development", map_path, audit_path,
                        expected)

    def test_direct_timing_qualification_forgery_is_rejected(self) -> None:
        with self.assertRaises(subject.ContractError):
            subject.TimingQualification(
                subject.contract_sha256(self.contract), "development",
                digest("map"), EXPECTED_TIMING_BASE_HASHES["development"],
                digest("domain"), "1" * 40, (), digest("audit"),
                tuple([0] * 2304), digest("trace roster"),
                tuple([digest("trace")] * 2304))

        forged = object.__new__(subject.TimingQualification)
        original = self.development_qualification
        for key in subject.TIMING_QUALIFICATION_FIELDS:
            object.__setattr__(forged, key, getattr(original, key))
        forged_offsets = list(original.retry_offsets)
        forged_offsets[0] = 1
        values = {
            "map_sha256": digest("syntactically valid forged map"),
            "qualified_domain_sha256":
                subject._timing_domain_sha256_from_offsets(
                    self.contract, "development", forged_offsets),
            "retry_offsets": tuple(forged_offsets),
        }
        for key, value in values.items():
            object.__setattr__(forged, key, value)
        with self.assertRaises(subject.ContractError):
            list(subject.iter_timing_cells(
                self.contract, "development", forged))

        trace_forgery = object.__new__(subject.TimingQualification)
        for key in subject.TIMING_QUALIFICATION_FIELDS:
            object.__setattr__(trace_forgery, key, getattr(original, key))
        substituted_traces = list(original.selected_trace_sha256s)
        substituted_traces[0] = digest("substituted selected trace")
        object.__setattr__(
            trace_forgery, "selected_trace_sha256s",
            tuple(substituted_traces))
        with self.assertRaises(subject.ContractError):
            list(subject.iter_timing_cells(
                self.contract, "development", trace_forgery))

        # Reproduce the stronger attack: every public field and digest is
        # self-consistent, but no strict audit loader created this identity.
        self_consistent = object.__new__(subject.TimingQualification)
        forged_offsets = tuple([1] * 2304)
        forged_traces = tuple([digest("forged terminal trace")] * 2304)
        forged_domain = subject._timing_domain_sha256_from_offsets(
            self.contract, "development", forged_offsets)
        forged_roster = subject.timing_selected_trace_roster_sha256(
            forged_traces)
        control_records = [dict(zip(
            subject.TIMING_QUALIFICATION_CONTROL_ORDER, control))
            for control in original.controls]
        forged_map = {
            "schema": subject.TIMING_QUALIFICATION_MAP_SCHEMA,
            "contract_sha256": subject.contract_sha256(self.contract),
            "phase": "development",
            "source_git_commit": original.source_git_commit,
            "base_domain_sha256": original.base_domain_sha256,
            "qualified_domain_sha256": forged_domain,
            "entry_kind": subject.TIMING_QUALIFICATION_ENTRY_KIND,
            "controls": control_records,
            "qualification_audit_sha256": "a" * 64,
            "selected_trace_roster_sha256": forged_roster,
            "retry_offsets": list(forged_offsets),
        }
        forged_values = {
            "contract_sha256": subject.contract_sha256(self.contract),
            "phase": "development",
            "map_sha256": subject.timing_qualification_map_sha256(
                forged_map),
            "base_domain_sha256": original.base_domain_sha256,
            "qualified_domain_sha256": forged_domain,
            "source_git_commit": original.source_git_commit,
            "controls": original.controls,
            "qualification_audit_sha256": "a" * 64,
            "retry_offsets": forged_offsets,
            "selected_trace_roster_sha256": forged_roster,
            "selected_trace_sha256s": forged_traces,
        }
        for key in subject.TIMING_QUALIFICATION_FIELDS:
            object.__setattr__(self_consistent, key, forged_values[key])
        with self.assertRaises(subject.ContractError):
            next(subject.iter_timing_cells(
                self.contract, "development", self_consistent))

        # Even object.__setattr__ cannot mutate an identity after the loader
        # placed its complete immutable snapshot in the provenance registry.
        with tempfile.TemporaryDirectory() as temporary:
            loaded = write_timing_qualification(
                Path(temporary), self.contract, "development")[0]
            object.__setattr__(loaded, "map_sha256", digest("mutated map"))
            with self.assertRaises(subject.ContractError):
                next(subject.iter_timing_cells(
                    self.contract, "development", loaded))

    def test_timing_trace_manifest_must_replay_terminal_audit_rows(
            self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            qualification, _map, _audit, _map_sha, traces = \
                write_timing_qualification(root, self.contract)
            traces[0] = digest("substituted terminal timing trace")
            path = write_timing_trace_manifest(
                root, self.contract, traces, qualification)
            manifest_sha256 = subject.trace_manifest_sha256(
                self.contract, "timing", "development", path,
                qualification)
            with self.assertRaises(subject.ContractError):
                subject.load_timing_trace_manifest(
                    self.contract, "development", path, manifest_sha256,
                    qualification)

    def test_timing_batch_formula_and_production_training_seed_are_exact(
            self) -> None:
        protocol = self.contract["timing"]["panel_protocol"]
        geometry = self.contract["timing"]["execution_geometry"]
        self.assertEqual(protocol["measured_batch_block_target"], 65536)
        self.assertEqual(
            protocol["invocations_per_slot"],
            "max(2,ceil(measured_batch_block_target/K))")
        self.assertEqual(protocol["interleave_policy"],
                         subject.TIMING_INTERLEAVE_POLICY)
        self.assertEqual(protocol["slots_per_panel"], 8)
        self.assertEqual(
            geometry, subject.EXPECTED_TIMING_EXECUTION_GEOMETRY)
        self.assertEqual(
            geometry["cohort_identity"],
            "all timing cell fields except replicate, base_loss_seed, "
            "base_cell_sha256, loss_retry_offset, and loss_seed, plus panel")
        self.assertEqual(geometry["timing_worker_count"], 8)
        self.assertEqual(geometry["jobs_per_wave"], 8)
        for domain in self.contract["timing"]["domains"].values():
            self.assertEqual(domain["paired_repetitions"], 96)
            self.assertEqual(domain["independent_rounds"], 12)
            self.assertEqual(
                domain["paired_repetitions"],
                domain["independent_rounds"] *
                geometry["timing_worker_count"])
        self.assertEqual(
            self.contract["timing"]["statistics"][
                "t_critical_by_independent_rounds"]["12"],
            2.200985160082949)
        arms = ("wirehair2_head", "wirehair1", "candidate")
        multi_arms = arms + ("candidate2",)
        self.assertEqual(subject.timing_cohort_count(
            self.contract, "development", arms), 264)
        self.assertEqual(subject.timing_cohort_count(
            self.contract, "development", multi_arms), 432)
        self.assertEqual(subject.timing_cohort_count(
            self.contract, "final", arms), 1650)
        for cohort in range(264):
            for independent_round in range(12):
                slots = {
                    subject.timing_worker_slot(
                        self.contract, "development", arms, cohort,
                        independent_round * 8 + lane)
                    for lane in range(8)
                }
                self.assertEqual(slots, set(range(8)))
        for replicate in range(96):
            counts = [0] * 8
            for cohort in range(264):
                counts[subject.timing_worker_slot(
                    self.contract, "development", arms,
                    cohort, replicate)] += 1
            self.assertEqual(counts, [33] * 8)
        self.assertEqual(
            subject.timing_worker_slot(
                self.contract, "development", arms, 263, 15),
            (15 % 8 + 263 + 1) % 8)
        self.assertEqual(
            subject.timing_worker_slot(
                self.contract, "final", arms, 1649, 95),
            (95 % 8 + 1649 + 11) % 8)
        for phase, cohort, replicate in (
                ("missing", 0, 0), (True, 0, 0),
                ("development", -1, 0),
                ("development", True, 0), ("development", 0, True),
                ("development", 264, 0), ("development", 0, 96),
                ("final", 1650, 0), ("final", 0, 96)):
            with self.assertRaises(subject.ContractError):
                subject.timing_worker_slot(
                    self.contract, phase, arms, cohort, replicate)
        self.assertEqual(
            subject.timing_worker_slot(
                self.contract, "development", multi_arms, 431, 95),
            (95 % 8 + 431 + 11) % 8)
        with self.assertRaises(subject.ContractError):
            subject.timing_worker_slot(
                self.contract, "development", multi_arms, 432, 0)
        for invalid_arms in (
                (), ("wirehair2_head", "candidate"),
                ("wirehair2_head", "wirehair1", "candidate", "candidate"),
                "wirehair2_head"):
            with self.assertRaises(subject.ContractError):
                subject.timing_cohort_count(
                    self.contract, "development", invalid_arms)
        expected = {
            8: 8192, 32: 2048, 100: 656, 128: 512,
            512: 128, 1000: 66, 2048: 32, 5000: 14,
            8192: 8, 20000: 4, 32768: 2, 64000: 2,
            65536: 2, 65537: 2,
        }
        for K, invocations in expected.items():
            with self.subTest(K=K):
                self.assertEqual(subject.timing_invocations_per_slot(
                    self.contract, K), invocations)
                self.assertEqual(
                    subject.timing_invocations_for_elapsed_slot(
                        self.contract, K, 0),
                    (invocations + 1) // 2)
                self.assertEqual(
                    subject.timing_invocations_for_elapsed_slot(
                        self.contract, K, 7),
                    invocations // 2)
        self.assertEqual(subject.timing_invocations_per_slot(
            self.contract, 6), 10923)
        self.assertEqual([
            subject.timing_invocations_for_elapsed_slot(
                self.contract, 6, slot) for slot in range(8)
        ], [5462] * 4 + [5461] * 4)
        for K in (0, -1, True, 1.0):
            with self.subTest(invalid_K=K), \
                    self.assertRaises(subject.ContractError):
                subject.timing_invocations_per_slot(self.contract, K)
        for slot in (-1, 8, True, 1.0):
            with self.subTest(invalid_slot=slot), \
                    self.assertRaises(subject.ContractError):
                subject.timing_invocations_for_elapsed_slot(
                    self.contract, 8, slot)

        development = list(subject.iter_timing_cells(
            self.contract, "development", self.development_qualification))
        self.assertEqual(self.contract["timing"]["domains"]["development"]
                         ["seed_mode"], "production_training")
        self.assertEqual({cell["base_seed_attempt"] for cell in development},
                         {self.contract["seeds"]
                          ["production_base_seed_attempt"]})
        for cell in development:
            self.assertEqual(cell["fixed_received_overhead"], 4)
            self.assertEqual(cell["receive_overhead_cap"], 256)
            base_cell = {
                key: cell[key]
                for key in self.contract["timing"]["base_cell_key"]
            }
            self.assertEqual(
                cell["base_cell_sha256"], subject.sha256_json(base_cell))
            self.assertEqual(
                cell["invocations_per_slot"],
                max(2, (65536 + cell["K"] - 1) // cell["K"]))
            self.assertEqual(cell["interleave_policy"],
                             subject.TIMING_INTERLEAVE_POLICY)
        first_by_replicate = {}
        for cell in development:
            first_by_replicate.setdefault(cell["replicate"], cell)
        mask = (1 << 64) - 1
        roots = self.contract["seeds"]["training_loss_roots"]
        for replicate, cell in first_by_replicate.items():
            salt = (replicate * 0x9e3779b97f4a7c15) & mask
            expected_seed = subject._splitmix64(
                int(roots[replicate % len(roots)], 16) ^ salt)
            self.assertEqual(cell["base_loss_seed"],
                             "0x{:016x}".format(expected_seed))
            self.assertEqual(cell["loss_seed"], cell["base_loss_seed"])
        cells_per_round = len(self.contract["k_sets"]["timing_short"]) * 2 * 8
        self.assertEqual(
            [cell["replicate"] for cell in development[:16]],
            list(range(8)) + list(range(8)))
        self.assertEqual(
            {cell["replicate"] // 8
             for cell in development[:cells_per_round]}, {0})
        self.assertEqual(
            {cell["replicate"] // 8
             for cell in development[cells_per_round:2 * cells_per_round]},
            {1})

    def test_timing_batch_contract_mutations_fail_closed(self) -> None:
        mutations = []
        wrong_target = copy.deepcopy(self.contract)
        wrong_target["timing"]["panel_protocol"][
            "measured_batch_block_target"] = 65535
        mutations.append(("target", wrong_target))
        bool_target = copy.deepcopy(self.contract)
        bool_target["timing"]["panel_protocol"][
            "measured_batch_block_target"] = True
        mutations.append(("target_type", bool_target))
        wrong_formula = copy.deepcopy(self.contract)
        wrong_formula["timing"]["panel_protocol"][
            "invocations_per_slot"] = "max(2,65536//K)"
        mutations.append(("formula", wrong_formula))
        missing_key = copy.deepcopy(self.contract)
        missing_key["timing"]["cell_key"].remove("invocations_per_slot")
        mutations.append(("cell_key", missing_key))
        wrong_receive_cap = copy.deepcopy(self.contract)
        wrong_receive_cap["timing"]["receive_overhead_cap"] = 255
        mutations.append(("receive_overhead_cap", wrong_receive_cap))
        missing_receive_cap_key = copy.deepcopy(self.contract)
        missing_receive_cap_key["timing"]["cell_key"].remove(
            "receive_overhead_cap")
        mutations.append(("receive_overhead_cap_cell_key",
                          missing_receive_cap_key))
        missing_interleave_key = copy.deepcopy(self.contract)
        missing_interleave_key["timing"]["cell_key"].remove(
            "interleave_policy")
        mutations.append(("interleave_cell_key", missing_interleave_key))
        forged_interleave = copy.deepcopy(self.contract)
        forged_interleave["timing"]["panel_protocol"][
            "interleave_policy"] = "forged"
        mutations.append(("interleave_policy", forged_interleave))
        wrong_round_count = copy.deepcopy(self.contract)
        wrong_round_count["timing"]["domains"]["development"][
            "independent_rounds"] = 96
        mutations.append(("independent_rounds", wrong_round_count))
        raw_timing = copy.deepcopy(self.contract)
        raw_timing["timing"]["domains"]["development"][
            "seed_mode"] = "raw_paired_training"
        mutations.append(("development_seed_mode", raw_timing))
        wrong_workers = copy.deepcopy(self.contract)
        wrong_workers["timing"]["execution_geometry"][
            "timing_worker_count"] = 7
        mutations.append(("worker_count", wrong_workers))
        bool_workers = copy.deepcopy(self.contract)
        bool_workers["timing"]["execution_geometry"][
            "timing_worker_count"] = True
        mutations.append(("worker_count_type", bool_workers))
        wrong_wave = copy.deepcopy(self.contract)
        wrong_wave["timing"]["execution_geometry"]["jobs_per_wave"] = 7
        mutations.append(("wave_size", wrong_wave))
        wrong_cohort = copy.deepcopy(self.contract)
        wrong_cohort["timing"]["execution_geometry"][
            "cohort_identity"] = "K and panel"
        mutations.append(("cohort_identity", wrong_cohort))
        wrong_rotation = copy.deepcopy(self.contract)
        wrong_rotation["timing"]["execution_geometry"][
            "worker_slot"] = "replicate%8"
        mutations.append(("worker_rotation", wrong_rotation))
        wrong_barrier = copy.deepcopy(self.contract)
        wrong_barrier["timing"]["execution_geometry"][
            "barrier"] = "after_each_cohort"
        mutations.append(("wave_barrier", wrong_barrier))
        wrong_cohort_formula = copy.deepcopy(self.contract)
        wrong_cohort_formula["timing"]["execution_geometry"][
            "cohort_count_formula"] = "24*11"
        mutations.append(("cohort_count_formula", wrong_cohort_formula))
        for name, changed in mutations:
            with self.subTest(name=name), \
                    self.assertRaises(subject.ContractError):
                subject._validate_structure(changed, check_domain_hashes=False)

    def test_aa_repeatability_rejects_exact_margin_equality(self) -> None:
        floor = math.log1p(
            self.contract["timing"]["practical_margin_ppm"] / 1000000.0)
        ab = {
            "selectable": True, "mean_log_ratio": 0.0,
            "lower_95": 0.0, "upper_95": 0.0,
        }
        stable_aa = dict(ab)
        for boundary_aa in (
                {**stable_aa, "lower_95": -floor},
                {**stable_aa, "upper_95": floor}):
            decision = subject._timing_group_decision(
                self.contract, ab, boundary_aa, stable_aa)
            self.assertFalse(decision["selectable"])
            self.assertEqual(
                decision["reason"], "aa_repeatability_threshold")
        inside = {
            **stable_aa,
            "lower_95": -floor + 1e-15,
            "upper_95": floor - 1e-15,
        }
        self.assertTrue(subject._timing_group_decision(
            self.contract, ab, inside, stable_aa)["selectable"])

    def test_every_k_hard_domain_is_575991_not_a_cross_product(self) -> None:
        cells = self.contract["recovery"]["domains"]["final_raw"][
            "expected_cells_per_arm"]
        self.assertEqual(cells, 63999 * 3 * 3)
        first = next(subject.iter_recovery_cells(self.contract, "final_raw"))
        self.assertEqual(first["K"], 2)
        self.assertEqual(first["base_seed_attempt"],
                         self.contract["seeds"]["raw_base_seed_attempts"][0])
        self.assertEqual(first["loss_seed"],
                         self.contract["seeds"]["training_loss_roots"][0])

    def test_v3_bands_and_short_cohorts_are_immutable(self) -> None:
        for mutation in ("band", "short", "timing"):
            changed = copy.deepcopy(self.contract)
            if mutation == "band":
                changed["k_bands"][0]["last"] = 99
                changed["k_bands"][1]["first"] = 100
            elif mutation == "short":
                changed["k_sets"]["short"].insert(9, 99)
            else:
                changed["k_sets"]["timing_short"][0] = 7
            with self.subTest(mutation=mutation), \
                    self.assertRaises(subject.ContractError):
                subject._validate_structure(changed, check_domain_hashes=False)

    def test_mutated_domain_hash_rejects(self) -> None:
        for evidence_kind in ("recovery", "timing"):
            changed = copy.deepcopy(self.contract)
            field = "domain_sha256" if evidence_kind == "recovery" else \
                "base_domain_sha256"
            changed[evidence_kind]["domains"]["development"][field] = \
                "0" * 64
            with self.subTest(evidence_kind=evidence_kind), \
                    self.assertRaises(subject.ContractError):
                subject._validate_structure(changed, check_domain_hashes=True)

    def test_v3_semantics_cannot_be_resealed_to_an_easier_contract(self) -> None:
        changes = []
        easy_loss = copy.deepcopy(self.contract)
        easy_loss["strata_sets"]["development"][1] = {
            "schedule": "easy", "loss_ppm": 1}
        easy_loss["strata_sets"]["hard"][0] = {
            "schedule": "easy", "loss_ppm": 1}
        changes.append(("loss", easy_loss))
        short_final = copy.deepcopy(self.contract)
        short_final["recovery"]["domains"]["final_raw"]["k_set"] = "short"
        short_final["recovery"]["domains"]["final_raw"][
            "expected_cells_per_arm"] = 270
        changes.append(("domain", short_final))
        easy_timing = copy.deepcopy(self.contract)
        easy_timing["timing"]["schedule"] = "made-up"
        easy_timing["timing"]["loss_ppm"] = -7
        changes.append(("timing", easy_timing))
        for name, changed in changes:
            with self.subTest(name=name), self.assertRaises(subject.ContractError):
                subject._validate_structure(changed, check_domain_hashes=False)

    def test_bool_float_numeric_aliases_reject_as_contract_errors(self) -> None:
        mutations = (
            ("goal", ("goal", "primary_overhead"), False),
            ("range", ("k_sets", "all", "first"), 2.0),
            ("count", ("recovery", "domains", "development",
                       "expected_cells_per_arm"), 360.0),
        )
        for name, path, replacement in mutations:
            changed = copy.deepcopy(self.contract)
            target = changed
            for key in path[:-1]:
                target = target[key]
            target[path[-1]] = replacement
            with self.subTest(name=name), self.assertRaises(subject.ContractError):
                subject._validate_structure(changed, check_domain_hashes=False)

    def test_malformed_roots_and_extreme_json_numbers_fail_closed(self) -> None:
        changed = copy.deepcopy(self.contract)
        changed["seeds"]["training_loss_roots"][0] = {}
        with self.assertRaises(subject.ContractError):
            subject._validate_structure(changed, check_domain_hashes=False)
        for payload in (b"1e400", b"1" + b"0" * 5000):
            with self.subTest(payload_bytes=len(payload)), \
                    self.assertRaises(subject.ContractError):
                subject._load_json_bytes(payload, "mutant")

    def test_loss_roots_are_disjoint_and_raw_attempts_are_explicit(self) -> None:
        seeds = self.contract["seeds"]
        self.assertIn(seeds["production_base_seed_attempt"],
                      seeds["raw_base_seed_attempts"])
        self.assertTrue(set(seeds["training_loss_roots"]).isdisjoint(
            seeds["validation_loss_roots"]))
        self.assertEqual(self.contract["recovery"]["raw_construction_attempts"], 1)
        self.assertEqual(
            self.contract["recovery"]["production_seed_fixups_in_raw_phase"], 0)

    def test_certified_realized_profile_receipt_is_canonical_little_endian(
            self) -> None:
        profile = subject.serialized_wh2_profile(2, 64, 7)
        self.assertEqual(len(profile), 32)
        self.assertEqual(profile[:8], b"WHV2\x01\x00\x20\x00")
        self.assertEqual(profile[8:16], bytes.fromhex("c9f9f447bb5b294b"))
        self.assertEqual(profile[16:24], (128).to_bytes(8, "little"))
        self.assertEqual(profile[24:29], b"\x40\x00\x00\x00\x07")
        self.assertEqual(profile[29:], b"\0\0\0")
        self.assertEqual(
            subject.realized_construction_sha256("0" * 64, 3, 64, 2),
            "de138eba17c4ee043e13ef0e545ff5c1a0ca631233058646ba414fad2126d12c")
        self.assertEqual(
            subject.generic_realized_construction_sha256(
                "wirehair2_certified", "0" * 64, 3, 64, 2),
            "de138eba17c4ee043e13ef0e545ff5c1a0ca631233058646ba414fad2126d12c")
        self.assertEqual(subject.freeze_manifest_sha256({
            "schema": subject.FREEZE_SCHEMA, "x": 1,
        }), "604f625718e3ec58a18d363eee71c73e0e3117ce60d4f5ec9f0a2583a30ec1a6")

    def test_uniform_raw_seed_schedule_goldens_are_K_independent(self) -> None:
        expected = {
            0: ("0x487468302aad7105", "0x4ec72102"),
            1: ("0xe6abe1e9a9f7ed1a", "0xecfe9abb"),
            2: ("0x84e35ba32942692f", "0x8b361474"),
            255: ("0xe1b6a7f5f5df09f0", "0xe8096049"),
        }
        for K in (3, 4, 1001, 3061, 5550, 7533):
            for attempt, seeds in expected.items():
                with self.subTest(K=K, attempt=attempt):
                    fields = raw_construction_fields(attempt)
                    self.assertEqual(fields["effective_precode_seed"], seeds[0])
                    self.assertEqual(fields["effective_packet_seed"], seeds[1])
                    # Validation and hashing accept exactly the same schedule
                    # for small, large, and production-fixup-active K values.
                    self.assertRegex(subject.raw_realized_construction_sha256(
                        "wirehair2_experiment", RAW_TEST_ARM,
                        digest("descriptor:" + RAW_TEST_ARM), K, 2, fields),
                        r"^[0-9a-f]{64}$")
        staircase_bytes = b"".join(
            subject._raw_staircase_for_K(K).to_bytes(2, "little")
            for K in range(2, 64001))
        self.assertEqual(
            hashlib.sha256(staircase_bytes).hexdigest(),
            "9ec9f5431db5342f01d8bff0dbbf943e80ecfcec7ac5f5eb078d72fad3a58dca")

    def test_raw_v2_freeze_binds_policy_and_rejects_legacy_raw_receipts(
            self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            trace_path = write_trace_manifest(root, self.contract)
            path = write_raw_freeze_manifest(root, self.contract, trace_path)
            loaded = subject.load_freeze_manifest(
                self.contract, "development", path)
            self.assertEqual(loaded["schema"], subject.RAW_FREEZE_SCHEMA)
            self.assertEqual(
                loaded["arms_by_name"][RAW_TEST_ARM]
                    ["seed_schedule_sha256"],
                subject.RAW_SEED_SCHEDULE_SHA256)

            original = json.loads(path.read_text(encoding="utf-8"))
            mutations = []
            missing = copy.deepcopy(original)
            del missing["arms"][2]["construction_seed_basis"]
            mutations.append(("missing policy", missing))
            forged_basis = copy.deepcopy(original)
            forged_basis["arms"][2]["construction_seed_basis"] = \
                subject.PRODUCTION_CONSTRUCTION_SEED_BASIS
            mutations.append(("forged basis", forged_basis))
            forged_schedule = copy.deepcopy(original)
            forged_schedule["arms"][2]["seed_schedule_sha256"] = digest(
                "forged schedule")
            mutations.append(("forged schedule", forged_schedule))
            forged_head = copy.deepcopy(original)
            forged_head["arms"][0]["construction_seed_basis"] = \
                subject.RAW_CONSTRUCTION_SEED_BASIS
            forged_head["arms"][0]["seed_schedule_sha256"] = \
                subject.RAW_SEED_SCHEDULE_SHA256
            mutations.append(("raw certified control", forged_head))
            nonraw = copy.deepcopy(original)
            nonraw["arm_roster"][2] = "candidate"
            nonraw["arm_roster_sha256"] = subject.arm_roster_sha256(
                nonraw["arm_roster"])
            nonraw["arms"][2]["arm"] = "candidate"
            mutations.append(("nonraw candidate", nonraw))
            for name, value in mutations:
                with self.subTest(name=name):
                    path.write_bytes((subject.canonical_json(value) + "\n").encode(
                        "utf-8"))
                    with self.assertRaises(subject.ContractError):
                        subject.load_freeze_manifest(
                            self.contract, "development", path)
            legacy = copy.deepcopy(original)
            legacy["schema"] = subject.FREEZE_SCHEMA
            for arm in legacy["arms"]:
                del arm["construction_seed_basis"]
                del arm["seed_schedule_sha256"]
            path.write_bytes((subject.canonical_json(legacy) + "\n").encode(
                "utf-8"))
            with self.assertRaises(subject.ContractError):
                subject.load_freeze_manifest(
                    self.contract, "development", path)

            # A legacy raw receipt cannot evade the schema boundary merely by
            # relabeling the experimental arm without the conventional raw
            # name prefix.
            legacy["arm_roster"][2] = "candidate"
            legacy["arm_roster_sha256"] = subject.arm_roster_sha256(
                legacy["arm_roster"])
            legacy["arms"][2]["arm"] = "candidate"
            legacy_descriptors = (
                "0550e0ed0c62d5491ff6915652fd96ed25f3c7782462da8c551636ec2e0294dd",
                "80f57b83eac9b8e3a19d8235cc067e39990980510e46588ddefa16f9561e1c38",
                "0eb3aef0602b5e7de15c822de84a5dbfc5dfdd99b76fbfd41538f7a13248c3a5",
                "2dc244661b3b073569319377ee3e55333a82ddad7bd328e1b0fef67395174614",
                "91d7c1a558e1cf93b002fcf2062b7657d301faca03972215495bdf2429499e90",
                "739092a7824449e6168f08b46661dfbe8ad5495ea4166b36073c79cd3bacdd11",
                "7c7889747a97ac160726b807fb03349344d49d4bec84c9e8220aa4689b00d2ca",
                "c70e0f57bb8d7783fa29b0decbed5da5058a8eb532d57d540f72108e114f091a",
            )
            for descriptor_sha256 in legacy_descriptors:
                with self.subTest(legacy_descriptor=descriptor_sha256):
                    legacy["arms"][2]["arm_descriptor_sha256"] = \
                        descriptor_sha256
                    path.write_bytes(
                        (subject.canonical_json(legacy) + "\n").encode(
                            "utf-8"))
                    with self.assertRaises(subject.ContractError):
                        subject.load_freeze_manifest(
                            self.contract, "development", path)
            malformed_descriptor = copy.deepcopy(legacy)
            malformed_descriptor["arms"][2]["arm_descriptor_sha256"] = []
            path.write_bytes(
                (subject.canonical_json(malformed_descriptor) + "\n").encode(
                    "utf-8"))
            with self.assertRaisesRegex(subject.ContractError, "not a SHA-256"):
                subject.load_freeze_manifest(
                    self.contract, "development", path)

    def test_v1_production_freeze_may_retain_selected_raw_arm_name(self) -> None:
        phase = "cross_width_validation"
        arms = ("wirehair2_head", "wirehair1", RAW_TEST_ARM)
        training_trace_sha256 = digest("production training traces")
        map_hashes = {
            arm: digest("production map:" + arm)
            for arm in arms if arm != "wirehair1"
        }
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            trace_path = write_trace_manifest(
                root, self.contract, phase=phase)
            path = write_freeze_manifest(
                root, self.contract, trace_path, arms, phase,
                map_hashes, training_trace_sha256)
            loaded = subject.load_freeze_manifest(
                self.contract, phase, path)
            self.assertEqual(loaded["schema"], subject.FREEZE_SCHEMA)
            self.assertEqual(
                loaded["arms_by_name"][RAW_TEST_ARM]["construction_policy"],
                "repair_map")

    def test_raw_v2_is_recovery_only_not_development_timing(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            hashes = timing_trace_hashes(
                self.contract, self.development_qualification)
            trace_path = write_timing_trace_manifest(
                root, self.contract, hashes, self.development_qualification)
            timing_path = write_timing_freeze_manifest(
                root, self.contract, trace_path,
                ("wirehair2_head", "wirehair1", "candidate"),
                self.development_qualification)
            # The existing timing-development v1 receipt remains valid, but
            # relabeling it as raw v2 is rejected before any result is scored.
            subject.load_freeze_manifest(
                self.contract, "development", timing_path, "timing",
                self.development_qualification)
            value = json.loads(timing_path.read_text(encoding="utf-8"))
            value["schema"] = subject.RAW_FREEZE_SCHEMA
            for arm in value["arms"]:
                arm["construction_seed_basis"] = \
                    subject.NOT_APPLICABLE_CONSTRUCTION_SEED_BASIS
                arm["seed_schedule_sha256"] = "0" * 64
            timing_path.write_bytes(
                (subject.canonical_json(value) + "\n").encode("utf-8"))
            with self.assertRaises(subject.ContractError):
                subject.load_freeze_manifest(
                    self.contract, "development", timing_path, "timing",
                    self.development_qualification)

    def test_final_freezes_bind_one_architecture_and_production_map(self) -> None:
        mutations = ("binary", "map", "training_trace")
        for mutation in mutations:
            with self.subTest(mutation=mutation), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                paths = write_final_freeze_manifests(
                    root, self.contract, self.final_qualification)
                summary = subject.validate_final_freeze_continuity(
                    self.contract, paths,
                    architecture_selection_receipt(
                        self.contract, self.development_qualification),
                    self.final_qualification, root / "timing-traces.jsonl")
                self.assertEqual(summary["selected_candidate"], "candidate")
                path = paths[("recovery", "final_validation")]
                value = json.loads(path.read_text(encoding="utf-8"))
                candidate = next(
                    arm for arm in value["arms"]
                    if arm["arm"] == "candidate")
                if mutation == "binary":
                    candidate["binary_sha256"] = digest("substitute binary")
                elif mutation == "map":
                    candidate["repair_map_sha256"] = digest("substitute map")
                else:
                    value["repair_training_trace_manifest_sha256"] = \
                        digest("substitute training trace")
                path.write_bytes(
                    (subject.canonical_json(value) + "\n").encode("utf-8"))
                with self.assertRaises(subject.ContractError):
                    subject.validate_final_freeze_continuity(
                        self.contract, paths,
                        architecture_selection_receipt(
                            self.contract, self.development_qualification),
                        self.final_qualification,
                        root / "timing-traces.jsonl")

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            paths = write_final_freeze_manifests(
                root, self.contract, self.final_qualification)
            with self.assertRaises(subject.ContractError):
                subject.validate_final_freeze_continuity(
                    self.contract, paths,
                    architecture_selection_receipt(
                        self.contract, self.development_qualification,
                        "substitute"), self.final_qualification,
                    root / "timing-traces.jsonl")
            substituted = architecture_selection_receipt(
                self.contract, self.development_qualification)
            substituted["selected_arm_descriptor_sha256"] = \
                digest("different selected equations")
            substituted["selected_architecture_sha256"] = \
                subject.selected_architecture_sha256("candidate", {
                    "codec": substituted["selected_codec"],
                    "arm_descriptor_sha256":
                        substituted["selected_arm_descriptor_sha256"],
                })
            unsigned = {key: value for key, value in substituted.items()
                        if key != "selection_sha256"}
            substituted["selection_sha256"] = subject.sha256_json(unsigned)
            with self.assertRaises(subject.ContractError):
                subject.validate_final_freeze_continuity(
                    self.contract, paths, substituted,
                    self.final_qualification, root / "timing-traces.jsonl")

    def test_final_continuity_rejects_valid_qualification_map_substitution(
            self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            original, map_path, audit_path, _map_sha, _traces = \
                write_timing_qualification(root, self.contract, "final")
            paths = write_final_freeze_manifests(
                root, self.contract, original)
            selection = architecture_selection_receipt(
                self.contract, self.development_qualification)
            trace_path = root / "timing-traces.jsonl"
            subject.validate_final_freeze_continuity(
                self.contract, paths, selection, original, trace_path)

            # Create a second strictly valid map with the same source,
            # controls, offsets, and qualified domain but a different terminal
            # trace roster.  The final timing freeze was made with the first
            # map and must not accept this substitution.
            rows = [json.loads(line) for line in
                    audit_path.read_text(encoding="utf-8").splitlines()]
            rows[0]["trace_sha256"] = digest("substituted final trace")
            audit_path.write_bytes("".join(
                subject.canonical_json(row) + "\n" for row in rows
            ).encode("utf-8"))
            map_value = json.loads(map_path.read_text(encoding="utf-8"))
            map_value["qualification_audit_sha256"] = \
                subject.timing_qualification_audit_sha256(
                    self.contract, "final", audit_path)
            map_value["selected_trace_roster_sha256"] = \
                subject.timing_selected_trace_roster_sha256(
                    [row["trace_sha256"] for row in rows])
            map_path.write_bytes(
                (subject.canonical_json(map_value) + "\n").encode("utf-8"))
            substitute = subject.load_timing_qualification_map(
                self.contract, "final", map_path, audit_path,
                subject.timing_qualification_map_sha256(map_value))
            self.assertEqual(
                substitute.qualified_domain_sha256,
                original.qualified_domain_sha256)
            self.assertNotEqual(substitute.map_sha256, original.map_sha256)
            with self.assertRaises(subject.ContractError):
                subject.validate_final_freeze_continuity(
                    self.contract, paths, selection, substitute, trace_path)

    def test_reserved_codec_kinds_are_bijective_in_final_freezes(self) -> None:
        for codec, policy, map_hash in (
                ("wirehair1", "not_applicable", "0" * 64),
                ("wirehair2_certified", "repair_map",
                 digest("production map:candidate"))):
            with self.subTest(codec=codec), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                paths = write_final_freeze_manifests(
                    root, self.contract, self.final_qualification)
                path = paths[("timing", "final")]
                value = json.loads(path.read_text(encoding="utf-8"))
                candidate = next(
                    arm for arm in value["arms"]
                    if arm["arm"] == "candidate")
                candidate["codec"] = codec
                candidate["construction_policy"] = policy
                candidate["repair_map_sha256"] = map_hash
                path.write_bytes(
                    (subject.canonical_json(value) + "\n").encode("utf-8"))
                with self.assertRaises(subject.ContractError):
                    subject.load_freeze_manifest(
                        self.contract, "final", path, "timing",
                        self.final_qualification)


class LedgerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = subject.load_contract()
        cls.trace_hashes = [
            digest("trace:" + subject.sha256_json(cell))
            for cell in subject.iter_recovery_cells(cls.contract, "development")
        ]
        cls.valid_rows = ledger_rows(
            cls.contract, trace_hashes=cls.trace_hashes)

    def validate(self, rows: Sequence[Mapping[str, Any]],
                 arms: Sequence[str] = (
                     "wirehair2_head", "wirehair1", "candidate"),
                 trace_hashes: Sequence[str] = ()) -> Mapping[str, Any]:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = write_ledger(root, rows)
            trace_path = write_trace_manifest(
                root, self.contract, trace_hashes or self.trace_hashes)
            freeze_path = write_freeze_manifest(
                root, self.contract, trace_path, arms)
            return subject.validate_ledger(
                self.contract, "development", path, freeze_path, trace_path)

    def assert_rejects(self, rows: Sequence[Mapping[str, Any]],
                       arms: Sequence[str] = (
                           "wirehair2_head", "wirehair1", "candidate")) -> None:
        with self.assertRaises(subject.ContractError):
            self.validate(rows, arms)

    def test_valid_common_ledger_reports_repairs_introductions_and_tails(self) -> None:
        def outcomes(arm: str, ordinal: int, _cell: Mapping[str, Any]) -> Outcome:
            if arm == "wirehair2_head" and ordinal in (0, 1):
                return "need_more_at_cap", None
            if arm == "candidate" and ordinal in (1, 2):
                return "construct_failed", None
            if ordinal == 3:
                return "success", 2
            return "success", 0

        summary = self.validate(ledger_rows(self.contract, outcome=outcomes))
        self.assertEqual(summary["schema"],
                         subject.SCHEMA + ".ledger-summary.v1")
        self.assertNotIn("freeze_schema", summary)
        self.assertEqual(summary["excluded_cells"], 0)
        self.assertEqual(summary["arms"]["wirehair2_head"]["failure_by_overhead"],
                         {"0": 3, "1": 3, "2": 2, "4": 2})
        self.assertEqual(summary["arms"]["candidate"]["failure_by_overhead"],
                         {"0": 3, "1": 3, "2": 2, "4": 2})
        comparison = summary["comparisons"]["candidate"]["controls"][
            "wirehair2_head"]
        self.assertEqual(comparison["overhead0_repairs"], 1)
        self.assertEqual(comparison["overhead0_introductions"], 1)
        self.assertEqual(comparison["overhead0_shared_failures"], 2)
        self.assertTrue(comparison["all_noninferiority_gates_pass"])
        self.assertIn("wirehair1",
                      summary["comparisons"]["candidate"]["controls"])

    def test_explicit_unsupported_rows_stay_in_denominator_and_charge_cap(self) -> None:
        def outcomes(_arm: str, ordinal: int, _cell: Mapping[str, Any]) -> Outcome:
            return ("unsupported", None) if ordinal < 40 else ("success", 0)

        summary = self.validate(ledger_rows(self.contract, outcome=outcomes))
        for arm in ("wirehair2_head", "wirehair1", "candidate"):
            self.assertEqual(summary["arms"][arm]["cells"], 360)
            self.assertEqual(summary["arms"][arm]["status_counts"]["unsupported"], 40)
            self.assertEqual(summary["arms"][arm]["capped_overhead_sum"], 160)
        self.assertFalse(
            summary["comparisons"]["candidate"]["all_controls_noninferior"])

    def test_mandatory_control_unsupported_makes_candidates_ineligible(self) -> None:
        def outcomes(arm: str, ordinal: int, _cell: Mapping[str, Any]) -> Outcome:
            if arm == "wirehair1" and ordinal == 0:
                return "unsupported", None
            return "success", 0

        summary = self.validate(ledger_rows(self.contract, outcome=outcomes))
        self.assertFalse(summary["mandatory_controls_supported"])
        self.assertFalse(summary["comparisons"]["candidate"][
            "architecture_eligible"])

        def construct_failure(
                arm: str, ordinal: int, _cell: Mapping[str, Any]) -> Outcome:
            if arm == "wirehair1" and ordinal == 0:
                return "construct_failed", None
            return "success", 0

        scored = self.validate(ledger_rows(
            self.contract, outcome=construct_failure))
        self.assertTrue(scored["mandatory_controls_supported"])
        self.assertTrue(scored["comparisons"]["candidate"][
            "architecture_eligible"])

    def test_missing_duplicate_extra_seed_swap_and_partial_roster_reject(self) -> None:
        cases: List[Tuple[str, List[Dict[str, Any]], Sequence[str]]] = []
        cases.append(("missing", copy.deepcopy(self.valid_rows[1:]),
                      ("wirehair2_head", "wirehair1", "candidate")))
        cases.append(("duplicate", copy.deepcopy(self.valid_rows + [self.valid_rows[0]]),
                      ("wirehair2_head", "wirehair1", "candidate")))
        extra = copy.deepcopy(self.valid_rows)
        extra[0]["K"] = 7
        cases.append(("extra K", extra,
                      ("wirehair2_head", "wirehair1", "candidate")))
        swapped = copy.deepcopy(self.valid_rows)
        swapped[0]["loss_seed"] = \
            self.contract["seeds"]["training_loss_roots"][1]
        cases.append(("seed swap", swapped,
                      ("wirehair2_head", "wirehair1", "candidate")))
        cases.append(("partial roster", copy.deepcopy(self.valid_rows),
                      ("wirehair2_head", "wirehair1", "candidate", "third")))
        for name, rows, arms in cases:
            with self.subTest(name=name):
                self.assert_rejects(rows, arms)

    def test_wrong_band_trial_loss_and_cell_hash_reject(self) -> None:
        mutations = (
            ("band", "101-1000"),
            ("trial", 2),
            ("loss_ppm", 350000),
            ("cell_sha256", "0" * 64),
        )
        for field, value in mutations:
            rows = copy.deepcopy(self.valid_rows)
            rows[0][field] = value
            with self.subTest(field=field):
                self.assert_rejects(rows)

    def test_trace_arm_binary_and_descriptor_drift_reject(self) -> None:
        for field in ("trace_sha256", "binary_sha256", "arm_descriptor_sha256",
                      "realized_construction_sha256", "repair_map_sha256"):
            rows = copy.deepcopy(self.valid_rows)
            # Row 3 is the Wirehair2 control's next cell.  This makes each field drift
            # either between common-cell arms or within one frozen arm.
            rows[3][field] = digest("mutated:" + field)
            with self.subTest(field=field):
                self.assert_rejects(rows)

    def test_duplicate_frozen_trace_contents_are_allowed(self) -> None:
        duplicate_traces = [digest("same legitimate trace")] * 360
        rows = ledger_rows(self.contract, trace_hashes=duplicate_traces)
        self.validate(rows, trace_hashes=duplicate_traces)

    def test_forged_trace_receipt_rejects_against_frozen_manifest(self) -> None:
        rows = copy.deepcopy(self.valid_rows)
        rows[0]["trace_sha256"] = digest("forged trace")
        self.assert_rejects(rows)

    def test_trace_manifest_is_hashed_and_loaded_from_one_open(self) -> None:
        first_rows = trace_rows(self.contract, self.trace_hashes)
        second_rows = copy.deepcopy(first_rows)
        second_rows[0]["trace_sha256"] = digest("unfrozen replacement")
        first_data = "".join(
            subject.canonical_json(row) + "\n" for row in first_rows
        ).encode("utf-8")
        second_data = "".join(
            subject.canonical_json(row) + "\n" for row in second_rows
        ).encode("utf-8")

        class ReplacedPath:
            def __init__(self) -> None:
                self.opens = 0

            def open(self, _mode: str) -> io.BytesIO:
                self.opens += 1
                return io.BytesIO(first_data if self.opens == 1 else second_data)

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            stable = root / "stable-traces.jsonl"
            stable.write_bytes(first_data)
            expected = subject.trace_manifest_sha256(
                self.contract, "recovery", "development", stable)
        replaced = ReplacedPath()
        loaded = subject.load_trace_manifest(
            self.contract, "development", replaced, expected)  # type: ignore[arg-type]
        self.assertEqual(replaced.opens, 1)
        self.assertEqual(loaded[0], first_rows[0]["trace_sha256"])

    def test_raw_attempt_or_repair_map_rejects(self) -> None:
        for field, value in (
                ("construction_attempt", 1),
                ("construction_attempt", True),
                ("repair_map_sha256", digest("unexpected repair map"))):
            rows = copy.deepcopy(self.valid_rows)
            rows[0][field] = value
            with self.subTest(field=field, value=value):
                self.assert_rejects(rows)

    def test_repaired_validation_replays_exact_authenticated_map(self) -> None:
        phase = "cross_width_validation"
        trace_hashes = [
            digest("cross-trace:" + subject.sha256_json(cell))
            for cell in subject.iter_recovery_cells(self.contract, phase)
        ]
        training_trace_sha = digest("frozen final_repaired trace manifest")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            map_paths: Dict[str, Path] = {}
            map_hashes: Dict[str, str] = {}
            attempts: Dict[str, Mapping[int, int]] = {}
            for arm, offset in (("wirehair2_head", 7), ("candidate", 8)):
                map_path, map_hash, arm_attempts = write_repair_map(
                    root, self.contract, arm, training_trace_sha, {2: offset})
                map_paths[arm] = map_path
                map_hashes[arm] = map_hash
                attempts[arm] = arm_attempts
            trace_path = write_trace_manifest(
                root, self.contract, trace_hashes, phase)
            freeze_path = write_freeze_manifest(
                root, self.contract, trace_path,
                ("wirehair2_head", "wirehair1", "candidate"), phase,
                map_hashes, training_trace_sha)
            rows = ledger_rows(
                self.contract, trace_hashes=trace_hashes, phase=phase,
                repair_attempts=attempts, repair_map_hashes=map_hashes)
            ledger_path = write_ledger(root, rows)
            summary = subject.validate_ledger(
                self.contract, phase, ledger_path, freeze_path, trace_path,
                map_paths)
            self.assertTrue(summary["arms"]["candidate"][
                "phase_recovery_gate_pass"])

            mutated = copy.deepcopy(rows)
            candidate_row = next(
                row for row in mutated
                if row["arm"] == "candidate" and row["K"] == 2)
            candidate_row["construction_attempt"] = 9
            candidate_row["realized_construction_sha256"] = \
                subject.generic_realized_construction_sha256(
                    "wirehair2_experiment",
                    candidate_row["arm_descriptor_sha256"], 2,
                    candidate_row["block_bytes"], 9)
            ledger_path = write_ledger(root, mutated)
            with self.assertRaises(subject.ContractError):
                subject.validate_ledger(
                    self.contract, phase, ledger_path, freeze_path, trace_path,
                    map_paths)

    def test_missing_declared_control_rejects_before_scoring(self) -> None:
        rows = ledger_rows(
            self.contract, arms=("wirehair2_head", "candidate"))
        self.assert_rejects(rows, ("wirehair2_head", "candidate"))

    def test_fatal_unknown_and_contradictory_outcomes_reject(self) -> None:
        mutations = (
            ("fatal", None),
            ("success", None),
            ("unsupported", 4),
            ("success", 5),
        )
        for outcome, extra in mutations:
            rows = copy.deepcopy(self.valid_rows)
            rows[0]["outcome"] = outcome
            rows[0]["decoded_extra"] = extra
            with self.subTest(outcome=outcome, extra=extra):
                self.assert_rejects(rows)

    def test_noncanonical_types_reject(self) -> None:
        for field, value in (("trial", True), ("K", 2.0), ("arm", ["control"])):
            rows = copy.deepcopy(self.valid_rows)
            rows[0][field] = value
            with self.subTest(field=field):
                self.assert_rejects(rows)

    def test_duplicate_json_key_and_partial_line_reject(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            duplicate = root / "duplicate.jsonl"
            duplicate.write_bytes(b'{"arm":"control","arm":"candidate"}\n')
            trace_path = write_trace_manifest(root, self.contract, self.trace_hashes)
            freeze_path = write_freeze_manifest(
                root, self.contract, trace_path,
                ("wirehair2_head", "wirehair1", "candidate"))
            with self.assertRaises(subject.ContractError):
                subject.validate_ledger(
                    self.contract, "development", duplicate,
                    freeze_path, trace_path)
            partial = write_ledger(root, self.valid_rows, terminal_newline=False)
            with self.assertRaises(subject.ContractError):
                subject.validate_ledger(
                    self.contract, "development", partial,
                    freeze_path, trace_path)

    def test_canonical_crlf_evidence_is_portable(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            trace_path = write_trace_manifest(
                root, self.contract, self.trace_hashes)
            freeze_path = write_freeze_manifest(
                root, self.contract, trace_path,
                ("wirehair2_head", "wirehair1", "candidate"))
            ledger_path = write_ledger(root, self.valid_rows)
            for path in (trace_path, freeze_path, ledger_path):
                path.write_bytes(path.read_bytes().replace(b"\n", b"\r\n"))
            subject.validate_ledger(
                self.contract, "development", ledger_path,
                freeze_path, trace_path)


class RawLedgerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = subject.load_contract()
        cls.trace_hashes = [
            digest("trace:" + subject.sha256_json(cell))
            for cell in subject.iter_recovery_cells(
                cls.contract, "development")
        ]
        cls.valid_rows = raw_ledger_rows(cls.contract, cls.trace_hashes)

    def validate(self, rows: Sequence[Mapping[str, Any]]) \
            -> Mapping[str, Any]:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            trace_path = write_trace_manifest(
                root, self.contract, self.trace_hashes)
            freeze_path = write_raw_freeze_manifest(
                root, self.contract, trace_path)
            ledger_path = write_ledger(root, rows)
            return subject.validate_ledger(
                self.contract, "development", ledger_path,
                freeze_path, trace_path)

    def assert_rejects(self, rows: Sequence[Mapping[str, Any]]) -> None:
        with self.assertRaises(subject.ContractError):
            self.validate(rows)

    def test_valid_raw_v2_ledger_has_separate_summary_schema(self) -> None:
        summary = self.validate(self.valid_rows)
        self.assertEqual(
            summary["schema"], subject.SCHEMA + ".ledger-summary.v2")
        self.assertEqual(summary["freeze_schema"], subject.RAW_FREEZE_SCHEMA)
        self.assertEqual(summary["arms"][RAW_TEST_ARM]["cells"], 360)

    def test_raw_rows_require_exact_policy_and_only_raw_rows_gain_fields(
            self) -> None:
        missing = copy.deepcopy(self.valid_rows)
        target = next(row for row in missing if row["arm"] == RAW_TEST_ARM)
        del target["construction_seed_basis"]
        self.assert_rejects(missing)

        extra_on_control = copy.deepcopy(self.valid_rows)
        control = next(
            row for row in extra_on_control if row["arm"] == "wirehair2_head")
        control.update(raw_construction_fields(control["construction_attempt"]))
        self.assert_rejects(extra_on_control)

        for field, value in (
                ("construction_seed_basis",
                 subject.PRODUCTION_CONSTRUCTION_SEED_BASIS),
                ("seed_schedule_sha256", digest("forged schedule")),
                ("effective_precode_seed", "0x0000000000000000"),
                ("effective_packet_seed", "0x00000000"),
                ("precode_attempt", True),
                ("staircase", 0),
                ("binary_dense_rows", 1.0),
                ("gf256_heavy_rows", -1),
                ("source_hits", 9),
                ("dense_identity_corner", 1),
                ("heavy_family", "unknown"),
                ("mix_count", 0)):
            rows = copy.deepcopy(self.valid_rows)
            raw = next(row for row in rows if row["arm"] == RAW_TEST_ARM)
            raw[field] = value
            with self.subTest(field=field):
                self.assert_rejects(rows)

    def test_raw_paired_attempt_cannot_be_reinterpreted(self) -> None:
        rows = copy.deepcopy(self.valid_rows)
        raw = next(row for row in rows if row["arm"] == RAW_TEST_ARM)
        raw["precode_attempt"] = 1
        raw["effective_precode_seed"] = \
            subject._effective_raw_precode_seed(1)
        fields = {field: raw[field] for field in subject.RAW_CONSTRUCTION_FIELDS}
        raw["realized_construction_sha256"] = \
            subject.raw_realized_construction_sha256(
                "wirehair2_experiment", RAW_TEST_ARM,
                raw["arm_descriptor_sha256"], raw["K"],
                raw["block_bytes"], fields)
        self.assert_rejects(rows)

    def test_raw_descriptor_binds_arm_name_and_fixed_geometry(self) -> None:
        for field, value in (
                ("binary_dense_rows", 13),
                ("gf256_heavy_rows", 11),
                ("mix_count", 4),
                ("staircase", 3),
                ("source_hits", 3),
                ("dense_identity_corner", True)):
            rows = copy.deepcopy(self.valid_rows)
            raw = next(row for row in rows if row["arm"] == RAW_TEST_ARM)
            raw[field] = value
            fields = {
                item: raw[item] for item in subject.RAW_CONSTRUCTION_FIELDS
            }
            raw["realized_construction_sha256"] = \
                subject.raw_realized_construction_sha256(
                    "wirehair2_experiment", RAW_TEST_ARM,
                    raw["arm_descriptor_sha256"], raw["K"],
                    raw["block_bytes"], fields)
            with self.subTest(field=field):
                self.assert_rejects(rows)

        for label, mutate in (
                ("retired_descriptor", lambda manifest, arm: arm.__setitem__(
                    "arm_descriptor_sha256",
                    "0550e0ed0c62d5491ff6915652fd96ed25f3c7782462da8c551636ec2e0294dd")),
                ("arm_alias", lambda manifest, arm: (
                    arm.__setitem__("arm", "wirehair2_raw_fake"),
                    manifest["arm_roster"].__setitem__(
                        manifest["arm_roster"].index(RAW_TEST_ARM),
                        "wirehair2_raw_fake"),
                    manifest.__setitem__(
                        "arm_roster_sha256", subject.arm_roster_sha256(
                            manifest["arm_roster"]))))):
            with self.subTest(freeze_mutation=label), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                trace_path = write_trace_manifest(
                    root, self.contract, self.trace_hashes)
                freeze_path = write_raw_freeze_manifest(
                    root, self.contract, trace_path)
                manifest = json.loads(freeze_path.read_bytes())
                raw_arm = next(
                    arm for arm in manifest["arms"]
                    if arm["arm"] == RAW_TEST_ARM)
                mutate(manifest, raw_arm)
                freeze_path.write_bytes(
                    (subject.canonical_json(manifest) + "\n").encode("utf-8"))
                ledger_path = write_ledger(root, self.valid_rows)
                with self.assertRaises(subject.ContractError):
                    subject.validate_ledger(
                        self.contract, "development", ledger_path,
                        freeze_path, trace_path)

    def test_raw_realized_hash_is_sensitive_to_every_variable_field(self) -> None:
        descriptor = digest("descriptor:" + RAW_TEST_ARM)
        base = raw_construction_fields(0)
        baseline = subject.raw_realized_construction_sha256(
            "wirehair2_experiment", RAW_TEST_ARM, descriptor, 3, 2, base)
        variants = []
        precode = copy.deepcopy(base)
        precode["precode_attempt"] = 1
        precode["effective_precode_seed"] = \
            subject._effective_raw_precode_seed(1)
        variants.append(("precode_attempt", precode))
        packet = copy.deepcopy(base)
        packet["packet_attempt"] = 1
        packet["effective_packet_seed"] = \
            subject._effective_raw_packet_seed(1)
        variants.append(("packet_attempt", packet))
        for field, value in (
                ("staircase", 11), ("binary_dense_rows", 13),
                ("gf256_heavy_rows", 13), ("source_hits", 3),
                ("dense_identity_corner", True), ("mix_count", 2)):
            changed = copy.deepcopy(base)
            changed[field] = value
            variants.append((field, changed))
        for field, changed in variants:
            with self.subTest(field=field):
                self.assertNotEqual(
                    subject.raw_realized_construction_sha256(
                        "wirehair2_experiment", RAW_TEST_ARM, descriptor,
                        3, 2, changed), baseline)
        for label, args in (
                ("arm", ("wirehair2_experiment",
                         "wirehair2_raw_other", descriptor, 3, 2, base)),
                ("descriptor", ("wirehair2_experiment", RAW_TEST_ARM,
                                digest("other descriptor"), 3, 2, base)),
                ("K", ("wirehair2_experiment", RAW_TEST_ARM,
                       descriptor, 4, 2, base)),
                ("block_bytes", ("wirehair2_experiment", RAW_TEST_ARM,
                                 descriptor, 3, 3, base))):
            with self.subTest(field=label):
                self.assertNotEqual(
                    subject.raw_realized_construction_sha256(*args), baseline)
        for field, value in (
                ("construction_seed_basis", "forged"),
                ("seed_schedule_sha256", digest("other schedule")),
                ("heavy_family", "other")):
            changed = copy.deepcopy(base)
            changed[field] = value
            with self.subTest(rejected_fixed_field=field), \
                    self.assertRaises(subject.ContractError):
                subject.raw_realized_construction_sha256(
                    "wirehair2_experiment", RAW_TEST_ARM, descriptor,
                    3, 2, changed)


class TimingReceiptTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = subject.load_contract()
        cls.arms = ("wirehair2_head", "wirehair1", "candidate")
        cls.qualification_temp = tempfile.TemporaryDirectory()
        (cls.qualification, cls.qualification_map_path,
         cls.qualification_audit_path, cls.qualification_map_sha256,
         _selected_traces) = write_timing_qualification(
             Path(cls.qualification_temp.name), cls.contract)
        cls.hashes = timing_trace_hashes(cls.contract, cls.qualification)
        cls.valid_rows = timing_receipt_rows(
            cls.contract, cls.hashes, cls.qualification)

    @classmethod
    def tearDownClass(cls) -> None:
        cls.qualification_temp.cleanup()

    def validate(self, rows: Sequence[Mapping[str, Any]]) -> Mapping[str, Any]:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            trace_path = write_timing_trace_manifest(
                root, self.contract, self.hashes, self.qualification)
            freeze_path = write_timing_freeze_manifest(
                root, self.contract, trace_path, self.arms,
                self.qualification)
            receipt_path = write_timing_receipt(root, rows)
            return subject.validate_timing_receipt(
                self.contract, "development", receipt_path,
                freeze_path, trace_path, self.qualification_map_path,
                self.qualification_audit_path,
                self.qualification_map_sha256)

    def assert_rejects(self, rows: Sequence[Mapping[str, Any]]) -> None:
        with self.assertRaises(subject.ContractError):
            self.validate(rows)

    def test_exact_panels_counterbalancing_and_student_t_gate(self) -> None:
        summary = self.validate(self.valid_rows)
        self.assertEqual(summary["rows"], 2304 * 11)
        self.assertTrue(summary["candidates"]["candidate"][
            "phase_speed_gate_pass"])
        decisions = summary["candidates"]["candidate"]["panels"]
        for panel in decisions.values():
            for decision in panel["by_band_width"].values():
                self.assertEqual(decision["status"], "faster")
                self.assertLess(decision["upper_95"],
                                -decision["effective_floor"])
                self.assertAlmostEqual(
                    decision["block_phase_difference"]["mean_log_ratio"],
                    0.0)

    def test_opposite_order_blocks_are_averaged_and_diagnosed(self) -> None:
        rows = copy.deepcopy(self.valid_rows)
        for row in rows:
            if row["panel_kind"] != "AB":
                continue
            left_primary = 64000
            right_primary = 100000
            primary = [left_primary, right_primary,
                       right_primary, left_primary] \
                if row["order"] == "ABBA" else \
                [right_primary, left_primary,
                 left_primary, right_primary]
            # The opposite block is an exact tie.  Its left/right-oriented
            # contrast is therefore zero regardless of the opposite order.
            row["elapsed_ns"] = primary + [100000] * 4
        summary = self.validate(rows)
        panel = {
            "panel_kind": "AB", "scope": "decoder_solve",
            "left_arm": "candidate", "right_arm": "wirehair2_head",
        }
        stats = summary["panel_statistics"][subject.canonical_json(panel)][
            "by_band_width"]["2-100|64"]
        self.assertAlmostEqual(stats["mean_log_ratio"], math.log(0.8),
                               places=12)
        block_phase = stats["block_phase_difference"]
        self.assertAlmostEqual(block_phase["mean_log_ratio"], math.log(0.64),
                               places=12)

    def test_student_t_uses_twelve_round_means_not_ninety_six_lanes(
            self) -> None:
        rows = copy.deepcopy(self.valid_rows)
        realized_round_logs: Dict[int, float] = {}
        for row in rows:
            if row["panel_kind"] != "AB":
                continue
            independent_round = row["replicate"] // 8
            lane = row["replicate"] % 8
            target_log = -0.10 + independent_round * 0.004 + \
                (lane - 3.5) * 0.0002
            right_ns = 1000000
            left_ns = int(round(right_ns * math.exp(target_log)))
            realized_log = math.log(left_ns) - math.log(right_ns)
            realized_round_logs.setdefault(independent_round, 0.0)
            # Every frozen K has this same lane value, so equal-K weighting
            # leaves the lane mean unchanged.
            if row["K"] == 8 and row["block_bytes"] == 64 and \
                    row["scope"] == "decoder_solve" and \
                    row["right_arm"] == "wirehair2_head":
                realized_round_logs[independent_round] += realized_log / 8.0
            primary = [left_ns, right_ns, right_ns, left_ns] \
                if row["order"] == "ABBA" else \
                [right_ns, left_ns, left_ns, right_ns]
            opposite = [right_ns, left_ns, left_ns, right_ns] \
                if row["order"] == "ABBA" else \
                [left_ns, right_ns, right_ns, left_ns]
            row["elapsed_ns"] = primary + opposite
        summary = self.validate(rows)
        panel = {
            "panel_kind": "AB", "scope": "decoder_solve",
            "left_arm": "candidate", "right_arm": "wirehair2_head",
        }
        stats = summary["panel_statistics"][subject.canonical_json(panel)][
            "by_band_width"]["2-100|64"]
        expected = subject._timing_confidence_interval(
            self.contract,
            [realized_round_logs[round_index] for round_index in range(12)],
            12)
        for field in ("mean_log_ratio", "lower_95", "upper_95"):
            self.assertAlmostEqual(stats[field], expected[field], places=12)

    def test_aggregate_equal_weights_every_k_and_width_cell(self) -> None:
        rows = copy.deepcopy(self.valid_rows)
        realized_by_cell: Dict[Tuple[int, int], float] = {}
        k_indexes = {
            K: index for index, K in enumerate(
                subject.EXPECTED_TIMING_SHORT_K)
        }
        for row in rows:
            if (row["panel_kind"] != "AB" or
                    row["scope"] != "decoder_solve" or
                    row["right_arm"] != "wirehair2_head"):
                continue
            target_log = -0.08 + 0.002 * k_indexes[row["K"]] + \
                (0.007 if row["block_bytes"] == 1280 else 0.0)
            right_ns = 1000000000
            left_ns = int(round(right_ns * math.exp(target_log)))
            realized_log = math.log(left_ns) - math.log(right_ns)
            realized_by_cell[(row["K"], row["block_bytes"])] = realized_log
            primary = [left_ns, right_ns, right_ns, left_ns] \
                if row["order"] == "ABBA" else \
                [right_ns, left_ns, left_ns, right_ns]
            opposite = [right_ns, left_ns, left_ns, right_ns] \
                if row["order"] == "ABBA" else \
                [left_ns, right_ns, right_ns, left_ns]
            row["elapsed_ns"] = primary + opposite
        self.assertEqual(len(realized_by_cell), 24)
        summary = self.validate(rows)
        panel = {
            "panel_kind": "AB", "scope": "decoder_solve",
            "left_arm": "candidate", "right_arm": "wirehair2_head",
        }
        aggregate = summary["panel_statistics"][
            subject.canonical_json(panel)]["aggregate"]
        expected = math.fsum(realized_by_cell.values()) / \
            len(realized_by_cell)
        self.assertAlmostEqual(
            aggregate["mean_log_ratio"], expected, places=12)

    def test_resolved_slowdown_fails_development_gate(self) -> None:
        rows = timing_receipt_rows(
            self.contract, self.hashes, self.qualification,
            candidate_scale=1.2)
        summary = self.validate(rows)
        self.assertFalse(summary["candidates"]["candidate"][
            "phase_speed_gate_pass"])
        statuses = {
            decision["status"]
            for panel in summary["candidates"]["candidate"]["panels"].values()
            for decision in panel["by_band_width"].values()
        }
        self.assertEqual(statuses, {"resolved_slower"})

    def test_aa_noise_cannot_expand_the_allowed_regression(self) -> None:
        rows = timing_receipt_rows(
            self.contract, self.hashes, self.qualification,
            candidate_scale=1.2)
        for row in rows:
            if row["panel_kind"] != "AA":
                continue
            primary = [200000, 100000, 100000, 200000] \
                if row["order"] == "ABBA" else \
                [100000, 200000, 200000, 100000]
            opposite = [100000, 200000, 200000, 100000] \
                if row["order"] == "ABBA" else \
                [200000, 100000, 100000, 200000]
            row["elapsed_ns"] = primary + opposite
        summary = self.validate(rows)
        candidate = summary["candidates"]["candidate"]
        self.assertFalse(candidate["phase_speed_gate_pass"])
        reasons = {
            decision.get("reason")
            for panel in candidate["panels"].values()
            for decision in panel["by_band_width"].values()
        }
        self.assertEqual(reasons, {"aa_repeatability_threshold"})

    def test_final_speed_gate_requires_all_four_declared_wins(self) -> None:
        faster = {"by_band_width": {"band|64": {"status": "faster"}},
                  "aggregate": {"status": "faster"}}
        decisions = {
            "decoder_solve|wirehair2_head": copy.deepcopy(faster),
            "receive_to_success|wirehair1": copy.deepcopy(faster),
            "encoder_init_plus_first_K_symbols|wirehair2_head":
                copy.deepcopy(faster),
            "encoder_init_plus_first_K_symbols|wirehair1":
                copy.deepcopy(faster),
        }
        self.assertTrue(subject._timing_speed_gate("final", decisions))
        mutations = (
            ("decoder_solve|wirehair2_head", "band", "noninferior"),
            ("receive_to_success|wirehair1", "band", "noninferior"),
            ("encoder_init_plus_first_K_symbols|wirehair2_head",
             "band", "unresolved"),
            ("encoder_init_plus_first_K_symbols|wirehair2_head",
             "aggregate", "noninferior"),
            ("encoder_init_plus_first_K_symbols|wirehair1",
             "aggregate", "noninferior"),
        )
        for panel, scope, status in mutations:
            changed = copy.deepcopy(decisions)
            if scope == "band":
                changed[panel]["by_band_width"]["band|64"]["status"] = status
            else:
                changed[panel]["aggregate"]["status"] = status
            with self.subTest(panel=panel, scope=scope):
                self.assertFalse(subject._timing_speed_gate("final", changed))

    def test_missing_order_trace_timer_and_identity_mutants_reject(self) -> None:
        cases: List[Tuple[str, List[Dict[str, Any]]]] = []
        cases.append(("missing", copy.deepcopy(self.valid_rows[1:])))
        wrong_order = copy.deepcopy(self.valid_rows)
        wrong_order[0]["order"] = \
            "BAAB" if wrong_order[0]["order"] == "ABBA" else "ABBA"
        cases.append(("order", wrong_order))
        wrong_trace = copy.deepcopy(self.valid_rows)
        wrong_trace[0]["trace_sha256"] = digest("wrong timing trace")
        cases.append(("trace", wrong_trace))
        missing_interleave = copy.deepcopy(self.valid_rows)
        del missing_interleave[0]["interleave_policy"]
        cases.append(("missing interleave policy", missing_interleave))
        forged_interleave = copy.deepcopy(self.valid_rows)
        forged_interleave[0]["interleave_policy"] = "forged"
        forged_cell = {
            key: forged_interleave[0][key]
            for key in self.contract["timing"]["cell_key"]
        }
        forged_interleave[0]["cell_sha256"] = subject.sha256_json(forged_cell)
        cases.append(("forged interleave policy", forged_interleave))
        zero_timer = copy.deepcopy(self.valid_rows)
        zero_timer[0]["elapsed_ns"][0] = 0
        cases.append(("timer", zero_timer))
        short_timer = copy.deepcopy(self.valid_rows)
        short_timer[0]["elapsed_ns"] = short_timer[0]["elapsed_ns"][:4]
        cases.append(("eight-slot cardinality", short_timer))
        wrong_batch = copy.deepcopy(self.valid_rows)
        wrong_batch[0]["invocations_per_slot"] += 1
        wrong_batch_cell = {
            key: wrong_batch[0][key]
            for key in self.contract["timing"]["cell_key"]
        }
        wrong_batch[0]["cell_sha256"] = subject.sha256_json(wrong_batch_cell)
        cases.append(("batch", wrong_batch))
        wrong_receive_cap = copy.deepcopy(self.valid_rows)
        wrong_receive_cap[0]["receive_overhead_cap"] = 255
        wrong_receive_cap_cell = {
            key: wrong_receive_cap[0][key]
            for key in self.contract["timing"]["cell_key"]
        }
        wrong_receive_cap[0]["cell_sha256"] = subject.sha256_json(
            wrong_receive_cap_cell)
        cases.append(("receive cap", wrong_receive_cap))
        bool_batch = copy.deepcopy(self.valid_rows)
        bool_batch[0]["invocations_per_slot"] = True
        cases.append(("batch type", bool_batch))
        float_batch = copy.deepcopy(self.valid_rows)
        float_batch[0]["invocations_per_slot"] = float(
            float_batch[0]["invocations_per_slot"])
        cases.append(("batch float", float_batch))
        descriptor_drift = copy.deepcopy(self.valid_rows)
        descriptor_drift[0]["left_arm_descriptor_sha256"] = digest("drift")
        cases.append(("identity", descriptor_drift))
        for name, rows in cases:
            with self.subTest(name=name):
                self.assert_rejects(rows)

    def test_receive_decoded_extra_uses_256_cap_while_solve_stays_at_4(
            self) -> None:
        receive_rows = [
            row for row in self.valid_rows
            if row["scope"] == "receive_to_success"
        ]
        solve_rows = [
            row for row in self.valid_rows if row["scope"] == "decoder_solve"
        ]
        self.assertTrue(receive_rows)
        self.assertTrue(solve_rows)
        self.assertEqual(
            {row[side + "_decoded_extra"]
             for row in receive_rows for side in ("left", "right")},
            {256})
        self.assertEqual(
            {row[side + "_decoded_extra"]
             for row in solve_rows for side in ("left", "right")},
            {4})

        above_receive_cap = copy.deepcopy(self.valid_rows)
        receive = next(
            row for row in above_receive_cap
            if row["scope"] == "receive_to_success")
        receive["left_decoded_extra"] = 257
        self.assert_rejects(above_receive_cap)

        solve_at_receive_cap = copy.deepcopy(self.valid_rows)
        solve = next(
            row for row in solve_at_receive_cap
            if row["scope"] == "decoder_solve")
        solve["left_decoded_extra"] = 256
        self.assert_rejects(solve_at_receive_cap)

    def test_explicit_failure_is_not_dropped_and_makes_panel_nonselectable(
            self) -> None:
        rows = copy.deepcopy(self.valid_rows)
        target = rows[0]
        target_cell = {key: target[key] for key in self.contract["timing"][
            "cell_key"]}
        target_scope = target["scope"]
        target_arm = target["left_arm"]
        for row in rows:
            if (all(row[key] == value for key, value in target_cell.items()) and
                    row["scope"] == target_scope):
                changed = False
                for side in ("left", "right"):
                    if row[side + "_arm"] == target_arm:
                        row[side + "_outcome"] = "need_more_at_cap"
                        row[side + "_decoded_extra"] = None
                        changed = True
                if changed:
                    row["elapsed_ns"] = [None] * 8
        summary = self.validate(rows)
        self.assertFalse(summary["candidates"]["candidate"][
            "phase_speed_gate_pass"])


class SelectionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = subject.load_contract()
        cls.qualification_temp = tempfile.TemporaryDirectory()
        cls.qualification = write_timing_qualification(
            Path(cls.qualification_temp.name), cls.contract)[0]

    @classmethod
    def tearDownClass(cls) -> None:
        cls.qualification_temp.cleanup()

    def test_recovery_equivalent_set_precedes_speed_ranking(self) -> None:
        artifact = digest("one development architecture roster")
        candidates = {
            "fast_but_worse": (7, -0.40),
            "reliable": (5, -0.05),
            "selected": (6, -0.20),
        }
        controls = ("wirehair2_head", "wirehair1")
        arm_failures = {
            **{arm: 0 for arm in controls},
            **{arm: value[0] for arm, value in candidates.items()},
        }
        arm_artifacts = {
            arm: {
                "codec": arm_codec(arm),
                "binary_sha256": digest("binary:" + arm),
                "arm_descriptor_sha256": digest("descriptor:" + arm),
            }
            for arm in arm_failures
        }
        recovery: Dict[str, Any] = {
            "schema": subject.SCHEMA + ".ledger-summary.v2",
            "phase": "development",
            "freeze_schema": subject.RAW_FREEZE_SCHEMA,
            "contract_sha256": subject.contract_sha256(self.contract),
            "domain_sha256": self.contract["recovery"]["domains"]
                ["development"]["domain_sha256"],
            "freeze_manifest_sha256": digest("recovery freeze"),
            "architecture_artifact_sha256": artifact,
            "arm_artifacts": copy.deepcopy(arm_artifacts),
            "excluded_cells": 0,
            "mandatory_controls_supported": True,
            "comparisons": {
                arm: {
                    "architecture_eligible": True,
                    "controls": {control: {} for control in controls},
                }
                for arm in candidates
            },
            "arms": {
                arm: {
                    "cells": 360,
                    "failure_by_overhead": {
                        "0": failure0, "1": max(0, failure0 - 2),
                        "2": max(0, failure0 - 4), "4": 0,
                    },
                }
                for arm, failure0 in arm_failures.items()
            },
        }
        timing: Dict[str, Any] = {
            "schema": subject.SCHEMA + ".timing-summary.v1",
            "phase": "development",
            "contract_sha256": subject.contract_sha256(self.contract),
            "base_domain_sha256": self.contract["timing"]["domains"]
                ["development"]["base_domain_sha256"],
            "domain_sha256": self.qualification.qualified_domain_sha256,
            "timing_qualification_map_sha256": self.qualification.map_sha256,
            "freeze_manifest_sha256": digest("timing freeze"),
            "architecture_artifact_sha256": artifact,
            "arm_artifacts": copy.deepcopy(arm_artifacts),
            "excluded_cells": 0,
            "rows": 2304 * (4 + 7 * len(candidates)),
            "candidates": {
                arm: {
                    "phase_speed_gate_pass": True,
                    "panels": {
                        "decoder_solve|wirehair2_head": {
                            "aggregate": {
                                "selectable": True,
                                "mean_log_ratio": mean,
                            },
                        },
                    },
                }
                for arm, (_failure0, mean) in candidates.items()
            },
        }
        result = subject.select_development_architecture(
            self.contract, recovery, timing, self.qualification)
        subject.validate_selection_receipt(self.contract, result)
        self.assertEqual(result["recovery_equivalence_allowance"], 1)
        self.assertEqual(
            result["recovery_equivalent_candidates"],
            ["reliable", "selected"])
        self.assertEqual(result["selected_arm"], "selected")
        self.assertEqual(result["candidate_roster"], sorted(candidates))
        tampered_receipt = copy.deepcopy(result)
        tampered_receipt["recovery_freeze_manifest_sha256"] = \
            digest("substituted after selection")
        with self.assertRaises(subject.ContractError):
            subject.validate_selection_receipt(
                self.contract, tampered_receipt)
        forged_minimum = copy.deepcopy(result)
        forged_minimum["minimum_overhead0_failures"] = 0
        unsigned = {key: value for key, value in forged_minimum.items()
                    if key != "selection_sha256"}
        forged_minimum["selection_sha256"] = subject.sha256_json(unsigned)
        with self.assertRaises(subject.ContractError):
            subject.validate_selection_receipt(
                self.contract, forged_minimum)
        with self.assertRaises(subject.ContractError):
            subject.validate_selection_receipt(
                self.contract,
                architecture_selection_receipt(
                    self.contract, self.qualification, "wirehair1"))

        mutations = []
        changed = copy.deepcopy(recovery)
        changed["schema"] = subject.SCHEMA + ".ledger-summary.v1"
        del changed["freeze_schema"]
        mutations.append(("legacy raw summary", changed, timing))
        changed = copy.deepcopy(recovery)
        changed["domain_sha256"] = digest("other recovery domain")
        mutations.append(("recovery domain", changed, timing))
        changed = copy.deepcopy(recovery)
        changed["excluded_cells"] = 1
        mutations.append(("recovery exclusion", changed, timing))
        changed = copy.deepcopy(recovery)
        changed["arms"]["selected"]["cells"] = 1
        mutations.append(("recovery cardinality", changed, timing))
        changed = copy.deepcopy(recovery)
        changed["arms"]["selected"]["failure_by_overhead"]["1"] = 7
        mutations.append(("nonnested tail", changed, timing))
        changed_timing = copy.deepcopy(timing)
        changed_timing["domain_sha256"] = digest("other timing domain")
        mutations.append(("timing domain", recovery, changed_timing))
        changed_timing = copy.deepcopy(timing)
        changed_timing["rows"] = 1
        mutations.append(("timing cardinality", recovery, changed_timing))
        changed_timing = copy.deepcopy(timing)
        changed_timing["freeze_manifest_sha256"] = "bad"
        mutations.append(("timing freeze", recovery, changed_timing))
        changed_timing = copy.deepcopy(timing)
        changed_timing["architecture_artifact_sha256"] = digest("other binary")
        mutations.append(("architecture", recovery, changed_timing))
        for name, recovery_value, timing_value in mutations:
            with self.subTest(mutation=name), \
                    self.assertRaises(subject.ContractError):
                subject.select_development_architecture(
                    self.contract, recovery_value, timing_value,
                    self.qualification)


if __name__ == "__main__":
    unittest.main()
