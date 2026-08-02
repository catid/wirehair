#!/usr/bin/env python3
"""Tests for fail-closed native WH2 short-screen assembly."""

from __future__ import annotations

import copy
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_benchmark_contract as contract_api
import wh2_native_short_screen as subject


ARMS = ("wirehair2_head", "wirehair1", "candidate")
WORKER_CPUS = tuple(range(8))
WORKER_BINARY = hashlib.sha256(b"one native worker").hexdigest()


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
        contract: Mapping[str, Any], kind: str, phase: str
        ) -> Tuple[List[Dict[str, Any]], List[str]]:
    cells = contract_api.iter_recovery_cells(contract, phase) if \
        kind == "recovery" else contract_api.iter_timing_cells(contract, phase)
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
            "packet_count": cell["K"] + 4,
            "candidate_count": cell["K"] + 4,
        })
    return rows, hashes


def freeze_manifest(
        contract: Mapping[str, Any], kind: str, phase: str,
        trace_sha256: str) -> Dict[str, Any]:
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
        "evidence_kind": kind,
        "phase": phase,
        "domain_sha256": contract[kind]["domains"][phase]["domain_sha256"],
        "source_git_commit": "1" * 40,
        "arm_roster": list(ARMS),
        "arm_roster_sha256": contract_api.arm_roster_sha256(ARMS),
        "trace_manifest_sha256": trace_sha256,
        "repair_training_trace_manifest_sha256": "0" * 64,
        "commands": [["wirehair_v2_contract_worker", kind, phase]],
        "cpu_affinity": list(WORKER_CPUS),
        "host_identity": {
            "name": "native-fixture", "controller_cpu": 8,
        },
        "arms": records,
    }


def realized(arm: str, K: int, block_bytes: int, attempt: int) -> str:
    return contract_api.generic_realized_construction_sha256(
        arm_codec(arm), descriptor(arm), K, block_bytes, attempt)


def envelope(ordinal: int, message_ordinal: int,
             payload: Mapping[str, Any], schema: str,
             cpu: Optional[int] = None,
             started_monotonic_ns: Optional[int] = None,
             finished_monotonic_ns: Optional[int] = None) -> Dict[str, Any]:
    kind = "recovery" if schema == subject.RECOVERY_RECORD_SCHEMA else "timing"
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


def recovery_records(
        contract: Mapping[str, Any], hashes: Sequence[str]
        ) -> List[Dict[str, Any]]:
    rows = []
    for cell_ordinal, cell in enumerate(contract_api.iter_recovery_cells(
            contract, "development")):
        for arm_index, arm in enumerate(ARMS):
            codec = arm_codec(arm)
            attempt = 0 if codec == "wirehair1" else \
                cell["base_seed_attempt"]
            payload = {
                "arm": arm,
                **cell,
                "outcome": "success",
                "decoded_extra": 0,
                "cell_sha256": contract_api.sha256_json(cell),
                "trace_sha256": hashes[cell_ordinal],
                "binary_sha256": WORKER_BINARY,
                "arm_descriptor_sha256": descriptor(arm),
                "construction_attempt": attempt,
                "realized_construction_sha256": realized(
                    arm, cell["K"], cell["block_bytes"], attempt),
                "repair_map_sha256": "0" * 64,
            }
            ordinal = cell_ordinal * len(ARMS) + arm_index
            rows.append(envelope(
                ordinal, cell_ordinal, payload,
                subject.RECOVERY_RECORD_SCHEMA))
    return rows


def timing_records(
        contract: Mapping[str, Any], hashes: Sequence[str]
        ) -> List[Dict[str, Any]]:
    rows = []
    panels = contract_api.timing_panels(contract, ARMS)
    cells = list(contract_api.iter_timing_cells(contract, "development"))
    cell_fields = contract["timing"]["cell_key"]
    stable_fields = [
        field for field in cell_fields
        if field not in ("replicate", "loss_seed")
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
    waves_per_cohort = repetitions // jobs_per_wave
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
                "elapsed_ns": [100000, 100000, 100000, 100000],
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
            wave_index = cell["replicate"] // jobs_per_wave
            global_wave = cohort_index * waves_per_cohort + wave_index
            position = cell["replicate"] % jobs_per_wave
            cpu = WORKER_CPUS[contract_api.timing_worker_slot(
                contract, "development", ARMS, cohort_index,
                cell["replicate"])]
            start = 2000000000 + global_wave * 1000000 + position * 10
            rows.append(envelope(
                ordinal, cell_ordinal, payload, subject.TIMING_RECORD_SCHEMA,
                cpu=cpu, started_monotonic_ns=start,
                finished_monotonic_ns=start + 500000))
    return rows


def timing_wave_rows(
        contract: Mapping[str, Any], records: Sequence[Mapping[str, Any]],
        ) -> Mapping[Tuple[int, int], List[int]]:
    panels = contract_api.timing_panels(contract, ARMS)
    cells = list(contract_api.iter_timing_cells(contract, "development"))
    stable_fields = [
        field for field in contract["timing"]["cell_key"]
        if field not in ("replicate", "loss_seed")
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
                    cpu_c: float = 61.0) -> Path:
    script_path = root / "sampler.py"
    script_path.write_text("# frozen sampler fixture\n", encoding="ascii")
    csv_path = root / "thermal.csv"
    with csv_path.open("w", encoding="ascii", newline="") as output:
        writer = csv.writer(output, lineterminator="\n")
        writer.writerow(subject.THERMAL_HEADER)
        for second in range(1, 5):
            writer.writerow([
                "2026-08-02T00:00:0{}.000Z".format(second),
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
        "window_start_monotonic_ns": 1000000000,
        "window_end_monotonic_ns": 4000000000,
        "terminal_status": "ok",
    }
    path = root / "sampler.json"
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

    def build_kind(self, root: Path, kind: str) -> Tuple[
            Path, Path, Path, List[Dict[str, Any]]]:
        phase = "development"
        native_traces, hashes = trace_records(self.contract, kind, phase)
        native_trace_path = root / (kind + "-native-traces.jsonl")
        trace_path = root / (kind + "-traces.jsonl")
        canonical_jsonl(native_trace_path, list(reversed(native_traces)))
        trace_sha = subject.assemble_trace_manifest(
            self.contract, kind, phase, native_trace_path, trace_path)
        freeze = freeze_manifest(self.contract, kind, phase, trace_sha)
        freeze_path = root / (kind + "-freeze.json")
        canonical_file(freeze_path, freeze)
        records = recovery_records(self.contract, hashes) if \
            kind == "recovery" else timing_records(self.contract, hashes)
        native_path = root / (kind + "-native-results.jsonl")
        canonical_jsonl(native_path, list(reversed(records)))
        return freeze_path, trace_path, native_path, records

    def test_full_native_recovery_and_one_candidate_timing_publish(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sampler_path = sampler_fixture(root)
            outputs = {}
            for kind, expected in (("recovery", 1080), ("timing", 4224)):
                freeze, traces, native, _ = self.build_kind(root, kind)
                output = root / (kind + "-output.jsonl")
                receipt = root / (kind + "-execution.json")
                result = subject.assemble_results(
                    self.contract, kind, "development", freeze, traces,
                    native, sampler_path, output, receipt,
                    verify_live_sampler=False)
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
                    self.contract, kind, "development", freeze, traces,
                    native, output, receipt, verify_live_sampler=False)
                self.assertEqual(
                    revalidated["execution_receipt"]["record_count"], expected)
                outputs[kind] = (freeze, traces, output)

            interpreters = [
                path for name in ("python3.8", "python3.12")
                for path in [shutil.which(name)] if path is not None
            ]
            self.assertTrue(interpreters, "at least one contract interpreter")
            for interpreter in interpreters:
                for kind in ("recovery", "timing"):
                    freeze, traces, output = outputs[kind]
                    command = "validate-ledger" if kind == "recovery" else \
                        "validate-timing"
                    completed = subprocess.run([
                        interpreter,
                        str(Path(contract_api.__file__).resolve()), command,
                        "--phase", "development", "--freeze-manifest",
                        str(freeze), "--trace-manifest", str(traces),
                        str(output),
                    ], stdout=subprocess.PIPE, stderr=subprocess.PIPE,
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
            "wirehair.wh2.native-timing-record.v2")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            _, _, _, recovery = self.build_kind(root, "recovery")
            self.assertTrue(all(
                row["schema"] == subject.RECOVERY_RECORD_SCHEMA
                for row in recovery))
            self.assertTrue(all(
                "invocations_per_slot" not in row["payload"]
                for row in recovery))

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            freeze_path, _, _, timing = self.build_kind(root, "timing")
            freeze = contract_api.load_freeze_manifest(
                self.contract, "development", freeze_path, "timing")
            self.assertTrue(all(
                row["schema"] == subject.TIMING_RECORD_SCHEMA
                for row in timing))
            for row in timing:
                payload = row["payload"]
                expected = max(1, (65536 + payload["K"] - 1) // payload["K"])
                self.assertEqual(payload["invocations_per_slot"], expected)
            legacy = copy.deepcopy(timing)
            legacy[0]["schema"] = "wirehair.wh2.native-timing-record.v1"
            with self.assertRaises(subject.NativeEvidenceError):
                subject._validate_native_records(
                    self.contract, freeze, "timing", "development", legacy)

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

    def test_terminal_timing_geometry_rejects_placement_and_barrier_drift(
            self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            freeze_path, _, _, records = self.build_kind(root, "timing")
            freeze = contract_api.load_freeze_manifest(
                self.contract, "development", freeze_path, "timing")
            waves = timing_wave_rows(self.contract, records)
            self.assertEqual(len(records), 4224)
            self.assertEqual(len(waves), 528)
            self.assertTrue(all(len(indexes) == 8
                                for indexes in waves.values()))
            subject._validate_native_records(
                self.contract, freeze, "timing", "development", records)

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
                subject._validate_native_records(
                    self.contract, freeze, "timing", "development",
                    cpu_remap)

            swapped = copy.deepcopy(records)
            panels = contract_api.timing_panels(self.contract, ARMS)
            cells = list(contract_api.iter_timing_cells(
                self.contract, "development"))
            first_identity = {
                field: cells[0][field]
                for field in self.contract["timing"]["cell_key"]
                if field not in ("replicate", "loss_seed")
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
                subject._validate_native_records(
                    self.contract, freeze, "timing", "development", swapped)

            ordered_wave_keys = [
                (cohort, wave)
                for cohort in range(264) for wave in range(2)
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
            subject._validate_native_records(
                self.contract, freeze, "timing", "development", equality)
            overlap = copy.deepcopy(equality)
            overlap[next_min_index]["started_monotonic_ns"] = prior_finish - 1
            with self.assertRaisesRegex(
                    subject.NativeEvidenceError, "overlap"):
                subject._validate_native_records(
                    self.contract, freeze, "timing", "development", overlap)

            with self.assertRaises(subject.NativeEvidenceError):
                subject._validate_native_records(
                    self.contract, freeze, "timing", "development",
                    [row for index, row in enumerate(records)
                     if index not in set(first_wave)])
            with self.assertRaises(subject.NativeEvidenceError):
                subject._validate_native_records(
                    self.contract, freeze, "timing", "development",
                    records + [copy.deepcopy(records[index])
                               for index in first_wave])

            for roster in (
                    list(range(7)), list(range(9)),
                    [1, 0, 2, 3, 4, 5, 6, 7]):
                changed_freeze = copy.deepcopy(freeze)
                changed_freeze["cpu_affinity"] = roster
                with self.subTest(roster=roster), self.assertRaisesRegex(
                        subject.NativeEvidenceError,
                        "exactly eight sorted worker CPUs"):
                    subject._validate_native_records(
                        self.contract, changed_freeze, "timing",
                        "development", records)

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
                subject._validate_native_records(
                    self.contract, sibling_freeze, "timing", "development",
                    sibling_records, sysfs_root=sibling_root)

            controller_freeze = copy.deepcopy(freeze)
            controller_freeze["host_identity"]["controller_cpu"] = 64
            with self.assertRaisesRegex(
                    subject.NativeEvidenceError, "controller shares"):
                subject._validate_native_records(
                    self.contract, controller_freeze, "timing", "development",
                    records, sysfs_root=sibling_root)

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

    def test_timing_batch_count_mutations_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            freeze, traces, _, records = self.build_kind(root, "timing")
            for mutation in ("missing", "wrong", "noncanonical"):
                case = root / mutation
                case.mkdir()
                changed = copy.deepcopy(records)
                payload = changed[0]["payload"]
                if mutation == "missing":
                    del payload["invocations_per_slot"]
                elif mutation == "wrong":
                    payload["invocations_per_slot"] += 1
                else:
                    # bool is deliberately chosen because it is an int
                    # subclass: the authoritative cell validator must use
                    # exact scalar types rather than isinstance(value, int).
                    payload["invocations_per_slot"] = True
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
                }[mutation]
                with self.subTest(mutation=mutation), \
                        self.assertRaisesRegex(
                            subject.NativeEvidenceError, expected_error):
                    subject.assemble_results(
                        self.contract, "timing", "development", freeze,
                        traces, native, sampler_fixture(case), output,
                        receipt, verify_live_sampler=False)
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
                freeze_path, _, _, records = self.build_kind(root, "recovery")
                freeze = contract_api.load_freeze_manifest(
                    self.contract, "development", freeze_path, "recovery")
                changed = copy.deepcopy(records)
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
                freeze, traces, _, records = self.build_kind(root, "timing")
                changed = copy.deepcopy(records)
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
                        self.contract, "timing", "development", freeze,
                        traces, native, sampler_fixture(root),
                        root / "output.jsonl", root / "receipt.json",
                        verify_live_sampler=False)
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
            freeze, traces, native, _ = self.build_kind(root, "recovery")
            output = root / "result.jsonl"
            receipt = root / "execution.json"
            original_validate = contract_api.validate_ledger

            def mutate_then_validate(*args, **kwargs):
                value = json.loads(freeze.read_text(encoding="utf-8"))
                value["source_git_commit"] = "2" * 40
                canonical_file(freeze, value)
                return original_validate(*args, **kwargs)

            with mock.patch.object(
                    contract_api, "validate_ledger",
                    side_effect=mutate_then_validate):
                with self.assertRaises(subject.NativeEvidenceError):
                    subject.assemble_results(
                        self.contract,
                        "recovery",
                        "development",
                        freeze,
                        traces,
                        native,
                        sampler_fixture(root),
                        output,
                        receipt,
                        verify_live_sampler=False)
            self.assertFalse(output.exists())
            self.assertFalse(receipt.exists())

    def test_execution_receipt_mutation_is_detected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            freeze, traces, native, _ = self.build_kind(root, "recovery")
            result = root / "result.jsonl"
            receipt = root / "execution.json"
            subject.assemble_results(
                self.contract, "recovery", "development", freeze, traces,
                native, sampler_fixture(root), result, receipt,
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
                self.contract, "recovery", "development", freeze, traces,
                native, result, receipt, verify_live_sampler=False)
            changed = json.loads(receipt.read_text(encoding="utf-8"))
            changed["thermal"]["edac_ce_max"] = 1
            unsigned = {key: value for key, value in changed.items()
                        if key != "receipt_sha256"}
            changed["receipt_sha256"] = contract_api.sha256_json(unsigned)
            mutant = root / "mutant-execution.json"
            canonical_file(mutant, changed)
            with self.assertRaises(subject.NativeEvidenceError):
                subject.validate_execution_receipt(
                    self.contract, "recovery", "development", freeze,
                    traces, native, result, mutant,
                    verify_live_sampler=False)

    def test_existing_destination_cannot_leave_a_partial_artifact_pair(
            self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            freeze, traces, native, _ = self.build_kind(root, "recovery")
            sampler = sampler_fixture(root)
            result = root / "result.jsonl"
            receipt = root / "execution.json"
            sentinel = b"preexisting\n"

            receipt.write_bytes(sentinel)
            with self.assertRaises(subject.NativeEvidenceError):
                subject.assemble_results(
                    self.contract, "recovery", "development", freeze,
                    traces, native, sampler, result, receipt,
                    verify_live_sampler=False)
            self.assertFalse(result.exists())
            self.assertEqual(receipt.read_bytes(), sentinel)

            receipt.unlink()
            result.write_bytes(sentinel)
            with self.assertRaises(subject.NativeEvidenceError):
                subject.assemble_results(
                    self.contract, "recovery", "development", freeze,
                    traces, native, sampler, result, receipt,
                    verify_live_sampler=False)
            self.assertEqual(result.read_bytes(), sentinel)
            self.assertFalse(receipt.exists())

    def test_post_link_failure_rolls_back_both_artifacts(self) -> None:
        for failure_point in ("result", "receipt"):
            with self.subTest(failure_point=failure_point), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                freeze, traces, native, _ = self.build_kind(root, "recovery")
                sampler = sampler_fixture(root)
                result = root / "result.jsonl"
                receipt = root / "execution.json"
                target = result if failure_point == "result" else receipt
                original_publish = subject._publish_new

                def publish_then_fail(staged, destination):
                    original_publish(staged, destination)
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
                            freeze,
                            traces,
                            native,
                            sampler,
                            result,
                            receipt,
                            verify_live_sampler=False)
                self.assertFalse(result.exists())
                self.assertFalse(receipt.exists())

    def test_result_mutation_during_revalidation_is_detected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            freeze, traces, native, _ = self.build_kind(root, "recovery")
            result = root / "result.jsonl"
            receipt = root / "execution.json"
            subject.assemble_results(
                self.contract,
                "recovery",
                "development",
                freeze,
                traces,
                native,
                sampler_fixture(root),
                result,
                receipt,
                verify_live_sampler=False)
            original_validate = contract_api.validate_ledger

            def reorder_then_validate(*args, **kwargs):
                rows = list(contract_api._parse_canonical_jsonl(
                    result, "mutation fixture"))
                canonical_jsonl(result, list(reversed(rows)))
                return original_validate(*args, **kwargs)

            with mock.patch.object(
                    contract_api, "validate_ledger",
                    side_effect=reorder_then_validate):
                with self.assertRaises(subject.NativeEvidenceError):
                    subject.validate_execution_receipt(
                        self.contract,
                        "recovery",
                        "development",
                        freeze,
                        traces,
                        native,
                        result,
                        receipt,
                        verify_live_sampler=False)

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
            self.assertEqual(summary["sample_count"], 4)

            lines = csv_path.read_text(encoding="ascii").splitlines(True)
            # Header, ignored historical row, first valid row, then a row
            # inside the selected interval: the latter must fail closed.
            lines[3] = "malformed in-window row\n"
            csv_path.write_text("".join(lines), encoding="ascii")
            with self.assertRaises(subject.NativeEvidenceError):
                subject._thermal_window(
                    sampler, 2000000000, 3000000000, WORKER_CPUS, False)

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
