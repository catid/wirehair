#!/usr/bin/env python3
"""Focused tests for the bounded native WH2 recovery-only controller."""

from __future__ import annotations

import copy
from contextlib import ExitStack, redirect_stderr
import hashlib
from io import StringIO
import os
from pathlib import Path
from types import SimpleNamespace
import sys
import tempfile
import time
from typing import Any, Dict, Mapping
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_benchmark_contract as contract_api
import wh2_native_short_screen as native_api
import wh2_run_native_recovery_screen as subject


FAKE_WORKER = r'''#!/usr/bin/env python3
import hashlib
import json
from pathlib import Path
import sys

def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"))

candidate = {
    "d12-h11-periodic": (
        "wirehair2_raw_d12_h11_periodic",
        "91d7c1a558e1cf93b002fcf2062b7657d301faca03972215495bdf2429499e90"),
    "d12-h13-periodic": (
        "wirehair2_raw_d12_h13_periodic",
        "7c7889747a97ac160726b807fb03349344d49d4bec84c9e8220aa4689b00d2ca"),
    "d13-h12-periodic": (
        "wirehair2_raw_d13_h12_periodic",
        "c70e0f57bb8d7783fa29b0decbed5da5058a8eb532d57d540f72108e114f091a"),
}
if len(sys.argv) == 3 and sys.argv[1] == "--describe-recovery-candidate":
    arm, descriptor = candidate.get(sys.argv[2], (None, None))
    if arm is None:
        raise SystemExit(2)
    arms = [
        {"arm": "wirehair2_head", "codec": "wirehair2_certified",
         "arm_descriptor_sha256":
             "4cafe27a8fb388ca9a4249b2c279b1406e7a0a86bcf14e98246988c7c503fa7a"},
        {"arm": "wirehair1", "codec": "wirehair1",
         "arm_descriptor_sha256":
             "d5a24d404e69efeb439907cd8271eba98d6af86b58efe159a820fb7aea08883d"},
        {"arm": "wirehair2_raw_d12_h12_periodic",
         "codec": "wirehair2_experiment",
         "arm_descriptor_sha256":
             "739092a7824449e6168f08b46661dfbe8ad5495ea4166b36073c79cd3bacdd11"},
        {"arm": arm, "codec": "wirehair2_experiment",
         "arm_descriptor_sha256": descriptor},
    ]
    for index, value in enumerate(arms):
        if index == 0:
            value["construction_seed_basis"] = "production-profile-v1"
            value["seed_schedule_sha256"] = "0" * 64
        elif index == 1:
            value["construction_seed_basis"] = "not-applicable"
            value["seed_schedule_sha256"] = "0" * 64
        else:
            value["construction_seed_basis"] = "uniform-raw-v1"
            value["seed_schedule_sha256"] = (
                "90a98a3db207852dabdf5fb27573ef48b"
                "ce52e0228cee4e291d96fa44ed509a7")
    print(canonical({
        "arms": arms,
        "binary_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "schema": "wirehair.wh2.native-worker-description.v2",
        "source_git_commit": "1" * 40,
    }))
    raise SystemExit(0)
raise SystemExit(2)
'''


class RecoveryOnlyRunnerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = contract_api.load_contract()

    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def _fake_worker(self) -> Path:
        path = self.root / "fake_recovery_worker.py"
        path.write_text(FAKE_WORKER, encoding="utf-8")
        path.chmod(0o755)
        return path

    def _description(self, candidate_id: str) -> Mapping[str, Any]:
        return subject.describe_candidate_worker(
            self._fake_worker(), candidate_id, time.monotonic() + 5.0)

    def _campaign(
            self, candidate_id: str, ordinal: int,
            source: str = "1" * 40, binary: str = "a" * 64,
            trace: bytes = b"identical trace\n") -> Mapping[str, Any]:
        candidate_arm, candidate_descriptor = \
            subject.CANDIDATE_BY_ID[candidate_id]
        controls = [
            {
                "arm": "wirehair2_head", "codec": "wirehair2_certified",
                "binary_sha256": binary,
                "arm_descriptor_sha256":
                    subject.CONTROL_DESCRIPTOR_SHA256S[0],
                "construction_policy": "raw_base",
                "repair_map_sha256": subject.ZERO_SHA256,
                "construction_seed_basis":
                    contract_api.PRODUCTION_CONSTRUCTION_SEED_BASIS,
                "seed_schedule_sha256": subject.ZERO_SHA256,
            },
            {
                "arm": "wirehair1", "codec": "wirehair1",
                "binary_sha256": binary,
                "arm_descriptor_sha256":
                    subject.CONTROL_DESCRIPTOR_SHA256S[1],
                "construction_policy": "not_applicable",
                "repair_map_sha256": subject.ZERO_SHA256,
                "construction_seed_basis":
                    contract_api.NOT_APPLICABLE_CONSTRUCTION_SEED_BASIS,
                "seed_schedule_sha256": subject.ZERO_SHA256,
            },
        ]
        raw_control = {
            "arm": subject.RAW_CONTROL_ARM[0],
            "codec": subject.RAW_CONTROL_ARM[1],
            "binary_sha256": binary,
            "arm_descriptor_sha256":
                subject.RAW_CONTROL_DESCRIPTOR_SHA256,
            "construction_policy": "raw_base",
            "repair_map_sha256": subject.ZERO_SHA256,
            "construction_seed_basis": subject.RAW_SEED_BASIS,
            "seed_schedule_sha256": subject.RAW_SEED_SCHEDULE_SHA256,
        }
        candidate = {
            "arm": candidate_arm, "codec": "wirehair2_experiment",
            "binary_sha256": binary,
            "arm_descriptor_sha256": candidate_descriptor,
            "construction_policy": "raw_base",
            "repair_map_sha256": subject.ZERO_SHA256,
            "construction_seed_basis": subject.RAW_SEED_BASIS,
            "seed_schedule_sha256": subject.RAW_SEED_SCHEDULE_SHA256,
        }
        arms = controls + [raw_control, candidate]
        roster = [value["arm"] for value in arms]
        worker = "/fixture/worker"
        description = {"resolved_path": worker}
        freeze = {
            "schema": contract_api.RAW_FREEZE_SCHEMA,
            "contract_sha256": contract_api.contract_sha256(self.contract),
            "evidence_kind": "recovery",
            "phase": "development",
            "domain_sha256": self.contract["recovery"]["domains"]
                ["development"]["domain_sha256"],
            "source_git_commit": source,
            "arm_roster": roster,
            "arm_roster_sha256": contract_api.arm_roster_sha256(roster),
            "trace_manifest_sha256": "d" * 64,
            "repair_training_trace_manifest_sha256": subject.ZERO_SHA256,
            "commands": subject._candidate_commands(
                description, candidate_id, list(range(8))),
            "cpu_affinity": list(range(8)),
            "host_identity": {"controller_cpu": 8},
            "arms": arms,
        }
        rows = []
        cells = self.contract["recovery"]["domains"]["development"][
            "expected_cells_per_arm"]
        for cell in range(cells):
            rows.extend((
                {"arm": "wirehair2_head", "cell": cell, "score": cell % 5},
                {"arm": "wirehair1", "cell": cell, "score": (cell + 1) % 5},
                {"arm": subject.RAW_CONTROL_ARM[0], "cell": cell,
                 "score": (cell + 2) % 5},
                {"arm": candidate_arm, "cell": cell,
                 "score": (cell + ordinal) % 5},
            ))
        run_hash = hashlib.sha256(candidate_id.encode("ascii")).hexdigest()
        result_hash = subject._hash_jsonl(rows)
        return {
            "directory": "/fixture/{}".format(candidate_id),
            "directory_identity": (1, ordinal + 1),
            "candidate_id": candidate_id,
            "candidate_arm": candidate_arm,
            "summary": {"summary_sha256": run_hash},
            "freeze": freeze,
            "receipt": {
                "freeze_manifest_sha256":
                    contract_api.freeze_manifest_sha256(freeze),
                "result_stream_sha256": result_hash,
                "receipt_sha256": "f" * 64,
            },
            "rows": rows,
            "trace_bytes": trace,
        }

    def _campaigns(self) -> list:
        return [
            self._campaign(candidate_id, ordinal)
            for ordinal, (candidate_id, _, _) in enumerate(
                subject.CANDIDATE_SPECS)
        ]

    def test_raw_seed_schedule_identity_is_frozen(self) -> None:
        self.assertEqual(subject.RAW_SEED_BASIS, "uniform-raw-v1")
        self.assertEqual(
            subject.RAW_SEED_SCHEDULE_SHA256,
            "90a98a3db207852dabdf5fb27573ef48bce52e0228cee4e291d96fa44ed509a7")

    def test_closed_candidate_descriptions_bind_exact_binary_and_roster(self) \
            -> None:
        worker = self._fake_worker()
        binary = hashlib.sha256(worker.read_bytes()).hexdigest()
        for candidate_id, candidate_arm, descriptor in \
                subject.CANDIDATE_SPECS:
            with self.subTest(candidate_id=candidate_id):
                description = subject.describe_candidate_worker(
                    worker, candidate_id, time.monotonic() + 5.0)
                self.assertEqual(description["binary_sha256"], binary)
                self.assertEqual(
                    [value["arm"] for value in description["arms"]],
                    ["wirehair2_head", "wirehair1",
                     subject.RAW_CONTROL_ARM[0], candidate_arm])
                self.assertEqual(
                    description["arms"][3]["arm_descriptor_sha256"],
                    descriptor)
                self.assertEqual(
                    description["arms"][2]["arm_descriptor_sha256"],
                    subject.RAW_CONTROL_DESCRIPTOR_SHA256)
        for invalid in ("", "D12-H11-PERIODIC", "d12-h12-periodic",
                        "d12-h11-periodic ", None):
            with self.subTest(invalid=invalid), \
                    self.assertRaises(subject.RecoveryRunnerError):
                subject.describe_candidate_worker(
                    worker, invalid, time.monotonic() + 5.0)

    def test_description_rejects_candidate_descriptor_substitution(self) \
            -> None:
        worker = self._fake_worker()
        original = subject.runner_api._run_command

        def tamper(argv, deadline, context):
            data = original(argv, deadline, context)
            value = subject.runner_api._parse_canonical_line(data, "fixture")
            value["arms"][3]["arm_descriptor_sha256"] = "0" * 64
            return (contract_api.canonical_json(value) + "\n").encode("utf-8")

        with mock.patch.object(subject.runner_api, "_run_command",
                               side_effect=tamper):
            with self.assertRaisesRegex(
                    subject.RecoveryRunnerError, "closed candidate ID"):
                subject.describe_candidate_worker(
                    worker, "d12-h11-periodic", time.monotonic() + 5.0)

    def test_description_rejects_raw_control_descriptor_substitution(self) \
            -> None:
        worker = self._fake_worker()
        original = subject.runner_api._run_command

        def tamper(argv, deadline, context):
            data = original(argv, deadline, context)
            value = subject.runner_api._parse_canonical_line(data, "fixture")
            value["arms"][2]["arm_descriptor_sha256"] = "0" * 64
            return (contract_api.canonical_json(value) + "\n").encode("utf-8")

        with mock.patch.object(subject.runner_api, "_run_command",
                               side_effect=tamper):
            with self.assertRaisesRegex(
                    subject.RecoveryRunnerError, "raw control"):
                subject.describe_candidate_worker(
                    worker, "d12-h11-periodic", time.monotonic() + 5.0)

    def test_description_rejects_control_descriptor_substitution(self) \
            -> None:
        worker = self._fake_worker()
        original = subject.runner_api._run_command

        def tamper(argv, deadline, context):
            data = original(argv, deadline, context)
            value = subject.runner_api._parse_canonical_line(data, "fixture")
            value["arms"][0]["arm_descriptor_sha256"] = "0" * 64
            return (contract_api.canonical_json(value) + "\n").encode("utf-8")

        with mock.patch.object(subject.runner_api, "_run_command",
                               side_effect=tamper):
            with self.assertRaisesRegex(
                    subject.RecoveryRunnerError, "control descriptor"):
                subject.describe_candidate_worker(
                    worker, "d12-h11-periodic", time.monotonic() + 5.0)

    def test_freeze_contains_only_exact_recovery_candidate_argv(self) -> None:
        candidate_id = "d12-h11-periodic"
        description = self._description(candidate_id)
        output = self.root / "freeze"
        output.mkdir()
        freeze = subject.write_recovery_freeze(
            self.contract, description, candidate_id, list(range(8)), 9,
            "1" * 40, "d" * 64, output)
        self.assertEqual(freeze["arm_roster"], [
            "wirehair2_head", "wirehair1",
            "wirehair2_raw_d12_h12_periodic",
            "wirehair2_raw_d12_h11_periodic"])
        self.assertEqual(freeze["commands"], [
            [description["resolved_path"],
             "--describe-recovery-candidate", candidate_id],
            [description["resolved_path"], "--emit-traces", "recovery"],
        ] + [[description["resolved_path"],
              "--recovery-candidate-worker", candidate_id, str(cpu)]
             for cpu in range(8)])
        flattened = " ".join(
            argument for command in freeze["commands"] for argument in command)
        self.assertNotIn("timing", flattened.lower())
        self.assertNotIn("--worker ", flattened)

        mutated = copy.deepcopy(freeze)
        mutated["commands"][-1][2] = "d12-h13-periodic"
        with self.assertRaisesRegex(
                subject.RecoveryRunnerError, "exact candidate argv"):
            subject._validate_candidate_freeze(mutated, candidate_id)
        mutated = copy.deepcopy(freeze)
        mutated["arms"][3]["codec"] = "routed_composite"
        with self.assertRaisesRegex(
                subject.RecoveryRunnerError, "candidate descriptor"):
            subject._validate_candidate_freeze(mutated, candidate_id)

    def test_job_domain_is_exactly_1440_recovery_commands(self) -> None:
        jobs = subject._candidate_recovery_jobs(self.contract)
        self.assertEqual(len(jobs), subject.RECOVERY_RECORDS)
        self.assertEqual(
            {job.ordinal for job in jobs}, set(range(subject.RECOVERY_RECORDS)))
        self.assertTrue(all(job.command().startswith(b"R ") for job in jobs))
        self.assertFalse(any(job.command().startswith(b"T ") for job in jobs))

    def test_logical_combination_has_six_arms_and_preserves_envelopes(self) \
            -> None:
        campaigns = self._campaigns()
        freeze, rows, bindings, trace = subject._combine_loaded_campaigns(
            self.contract, list(reversed(campaigns)))
        self.assertEqual(trace, b"identical trace\n")
        self.assertEqual(freeze["arm_roster"], [
            "wirehair2_head", "wirehair1",
            "wirehair2_raw_d12_h12_periodic",
            "wirehair2_raw_d12_h11_periodic",
            "wirehair2_raw_d12_h13_periodic",
            "wirehair2_raw_d13_h12_periodic",
        ])
        self.assertEqual(len(rows), subject.LOGICAL_RECORDS)
        self.assertEqual([row["arm"] for row in rows[:6]],
                         freeze["arm_roster"])
        self.assertEqual(
            [value["candidate_id"] for value in bindings],
            [value[0] for value in subject.CANDIDATE_SPECS])
        # Logical rows are payloads only.  Native envelope provenance remains
        # bound through each input receipt and is never relabeled six-arm work.
        self.assertTrue(all("work_sha256" not in row for row in rows))
        self.assertEqual(len(bindings), 3)

    def test_work_rank_binding_joins_every_effective_raw_construction(self) \
            -> None:
        raw_arms = (
            (subject.RAW_CONTROL_ARM[0],
             subject.RAW_CONTROL_DESCRIPTOR_SHA256, 12, 12),
            (subject.CANDIDATE_SPECS[0][1],
             subject.CANDIDATE_SPECS[0][2], 12, 11),
            (subject.CANDIDATE_SPECS[1][1],
             subject.CANDIDATE_SPECS[1][2], 12, 13),
            (subject.CANDIDATE_SPECS[2][1],
             subject.CANDIDATE_SPECS[2][2], 13, 12),
        )
        native_rows = []
        work_rows = []
        for arm, descriptor, dense_rows, heavy_rows in raw_arms:
            for cell in contract_api.iter_recovery_cells(
                    self.contract, "development"):
                attempt = cell["base_seed_attempt"]
                raw = {
                    "construction_seed_basis": subject.RAW_SEED_BASIS,
                    "seed_schedule_sha256":
                        subject.RAW_SEED_SCHEDULE_SHA256,
                    "precode_attempt": attempt,
                    "packet_attempt": attempt,
                    "effective_precode_seed":
                        contract_api._effective_raw_precode_seed(attempt),
                    "effective_packet_seed":
                        contract_api._effective_raw_packet_seed(attempt),
                    "staircase": 1 + cell["K"] % 7,
                    "binary_dense_rows": dense_rows,
                    "gf256_heavy_rows": heavy_rows,
                    "source_hits": 3,
                    "dense_identity_corner": False,
                    "heavy_family": "periodic-cauchy",
                    "mix_count": 3,
                }
                realized = contract_api.raw_realized_construction_sha256(
                    "wirehair2_experiment", arm, descriptor, cell["K"],
                    cell["block_bytes"], raw)
                cell_sha256 = contract_api.sha256_json(cell)
                trace_sha256 = hashlib.sha256(
                    ("{}:{}".format(arm, cell_sha256)).encode("ascii")
                ).hexdigest()
                native = {
                    **cell,
                    "arm": arm,
                    "arm_descriptor_sha256": descriptor,
                    "cell_sha256": cell_sha256,
                    "construction_attempt": attempt,
                    "realized_construction_sha256": realized,
                    "trace_sha256": trace_sha256,
                    **raw,
                }
                native_rows.append(native)
                for overhead in (0, 1, 2, 4):
                    work_rows.append({
                        **native,
                        "frozen_trace_sha256": trace_sha256,
                        "overhead": overhead,
                    })
        summary = {
            "source_provenance": {"source_git_commit": "1" * 40},
            "construction_seed_basis": subject.RAW_SEED_BASIS,
            "seed_schedule_sha256": subject.RAW_SEED_SCHEDULE_SHA256,
            "summary_sha256": "a" * 64,
            "result_stream_sha256": "b" * 64,
            "work_domain_sha256": "c" * 64,
        }
        binding = subject._bind_work_rank_identities(
            native_rows, {"summary": summary, "rows": work_rows},
            "1" * 40)
        self.assertEqual(
            binding["raw_identity_join_count"],
            subject.RAW_IDENTITY_JOIN_COUNT)
        self.assertRegex(binding["raw_identity_join_sha256"],
                         contract_api.SHA256)
        work_rows[0]["effective_packet_seed"] = "0x00000000"
        with self.assertRaisesRegex(
                subject.RecoveryRunnerError, "effective constructions differ"):
            subject._bind_work_rank_identities(
                native_rows, {"summary": summary, "rows": work_rows},
                "1" * 40)
        work_rows[0]["effective_packet_seed"] = \
            native_rows[0]["effective_packet_seed"]
        work_rows[0]["frozen_trace_sha256"] = "0" * 64
        with self.assertRaisesRegex(
                subject.RecoveryRunnerError, "trace identities differ"):
            subject._bind_work_rank_identities(
                native_rows, {"summary": summary, "rows": work_rows},
                "1" * 40)

    def test_combiner_rejects_identity_and_control_drift(self) -> None:
        mutations: Dict[str, Any] = {}
        controls = self._campaigns()
        controls[1]["rows"][0]["score"] += 1
        controls[1]["receipt"]["result_stream_sha256"] = \
            subject._hash_jsonl(controls[1]["rows"])
        mutations["cell-identical"] = controls
        source = self._campaigns()
        source[1]["freeze"]["source_git_commit"] = "2" * 40
        source[1]["receipt"]["freeze_manifest_sha256"] = \
            contract_api.freeze_manifest_sha256(source[1]["freeze"])
        mutations["source_git_commit"] = source
        binary = self._campaigns()
        binary[2]["freeze"]["arms"][2]["binary_sha256"] = "9" * 64
        binary[2]["receipt"]["freeze_manifest_sha256"] = \
            contract_api.freeze_manifest_sha256(binary[2]["freeze"])
        mutations["identical worker binary"] = binary
        trace = self._campaigns()
        trace[1]["trace_bytes"] = b"different trace\n"
        mutations["byte-identical"] = trace
        duplicate = self._campaigns()
        duplicate[2] = copy.deepcopy(duplicate[1])
        duplicate[2]["directory_identity"] = (1, 99)
        mutations["duplicate candidate"] = duplicate
        receipted = self._campaigns()
        receipted[2]["rows"][2]["score"] += 1
        mutations["result-stream hash"] = receipted
        frozen = self._campaigns()
        frozen[0]["freeze"]["host_identity"]["tampered"] = True
        mutations["freeze differs"] = frozen
        for diagnostic, campaigns in mutations.items():
            with self.subTest(diagnostic=diagnostic), \
                    self.assertRaisesRegex(
                        subject.RecoveryRunnerError, diagnostic):
                subject._combine_loaded_campaigns(self.contract, campaigns)

    def test_completed_campaign_binds_controller_to_frozen_host(self) -> None:
        campaign = self._campaigns()[0]
        directory = self.root / "campaign"
        directory.mkdir()
        freeze = campaign["freeze"]
        receipt = dict(campaign["receipt"])
        receipt.update({
            "record_count": subject.RECOVERY_RECORDS,
            "worker_cpus": list(range(8)),
            "thermal": {
                "sample_count": 2,
                "cpu_tctl_max_millic": 65000,
                "dimm_max_millic": 42000,
            },
        })
        summary = {
            "schema": subject.CAMPAIGN_SUMMARY_SCHEMA,
            "status": "complete",
            "output_dir": str(directory.resolve()),
            "candidate_id": campaign["candidate_id"],
            "candidate_arm": campaign["candidate_arm"],
            "source_git_commit": freeze["source_git_commit"],
            "contract_sha256": freeze["contract_sha256"],
            "domain_sha256": freeze["domain_sha256"],
            "trace_manifest_sha256": freeze["trace_manifest_sha256"],
            "worker_binary_sha256": freeze["arms"][0]["binary_sha256"],
            "controller_cpu": 8,
            "worker_cpus": receipt["worker_cpus"],
            "recovery_records": receipt["record_count"],
            "recovery_freeze_sha256": receipt["freeze_manifest_sha256"],
            "recovery_result_sha256": receipt["result_stream_sha256"],
            "recovery_execution_receipt_sha256": receipt["receipt_sha256"],
            "thermal_samples": receipt["thermal"]["sample_count"],
            "cpu_tctl_max_millic": receipt["thermal"]
                ["cpu_tctl_max_millic"],
            "dimm_max_millic": receipt["thermal"]["dimm_max_millic"],
            "construction_seed_basis": subject.RAW_SEED_BASIS,
            "seed_schedule_sha256": subject.RAW_SEED_SCHEDULE_SHA256,
        }
        summary["summary_sha256"] = contract_api.sha256_json(summary)

        def publish(value):
            (directory / "run-summary.json").write_text(
                contract_api.canonical_json(value) + "\n", encoding="utf-8")

        publish(summary)
        (directory / "recovery-traces.jsonl").write_bytes(
            campaign["trace_bytes"])
        for name in ("recovery-freeze.json",
                     "recovery-native-results.jsonl",
                     "recovery-execution.json"):
            (directory / name).write_bytes(b"fixture\n")

        def publish_rows(rows):
            with (directory / "recovery-results.jsonl").open(
                    "w", encoding="utf-8") as output:
                for row in rows:
                    output.write(contract_api.canonical_json(row) + "\n")

        publish_rows(campaign["rows"])
        validated = {"summary": {}, "execution_receipt": receipt}
        with mock.patch.object(
                subject.native_api, "validate_execution_receipt",
                return_value=validated), \
                mock.patch.object(
                    subject.contract_api, "load_freeze_manifest",
                    return_value=freeze), \
                mock.patch.object(
                    subject, "_validate_bound_recovery_trace_bytes"):
            loaded = subject.load_completed_campaign(self.contract, directory)
            self.assertEqual(loaded["candidate_id"], campaign["candidate_id"])
            self.assertEqual(
                loaded["summary"]["construction_seed_basis"],
                subject.RAW_SEED_BASIS)
            self.assertEqual(
                loaded["summary"]["seed_schedule_sha256"],
                subject.RAW_SEED_SCHEDULE_SHA256)
            freeze["host_identity"]["tampered_after_validation"] = True
            with self.assertRaisesRegex(
                    subject.RecoveryRunnerError, "reopened candidate freeze"):
                subject.load_completed_campaign(self.contract, directory)
            freeze["host_identity"].pop("tampered_after_validation")
            tampered_rows = copy.deepcopy(campaign["rows"])
            tampered_rows[2]["score"] += 1
            publish_rows(tampered_rows)
            with self.assertRaisesRegex(
                    subject.RecoveryRunnerError, "receipted result hash"):
                subject.load_completed_campaign(self.contract, directory)
            publish_rows(campaign["rows"])
            mutated = dict(summary)
            mutated["construction_seed_basis"] = "production-profile"
            mutated["summary_sha256"] = contract_api.sha256_json({
                key: value for key, value in mutated.items()
                if key != "summary_sha256"
            })
            publish(mutated)
            with self.assertRaisesRegex(
                    subject.RecoveryRunnerError, "validated execution"):
                subject.load_completed_campaign(self.contract, directory)
            mutated = dict(summary)
            mutated["controller_cpu"] = 9
            mutated["summary_sha256"] = contract_api.sha256_json({
                key: value for key, value in mutated.items()
                if key != "summary_sha256"
            })
            publish(mutated)
            with self.assertRaisesRegex(
                    subject.RecoveryRunnerError, "validated execution"):
                subject.load_completed_campaign(self.contract, directory)

    def test_reopened_trace_bytes_are_semantically_bound_to_freeze(self) \
            -> None:
        campaign = self._campaigns()[0]
        freeze = copy.deepcopy(campaign["freeze"])
        digest = contract_api._trace_manifest_hasher(
            self.contract, "recovery", "development")
        rows = []
        for ordinal, cell in enumerate(contract_api.iter_recovery_cells(
                self.contract, "development")):
            row = {
                "ordinal": ordinal,
                "cell_sha256": contract_api.sha256_json(cell),
                "trace_sha256": hashlib.sha256(
                    str(ordinal).encode("ascii")).hexdigest(),
            }
            contract_api._hash_trace_manifest_row(digest, row)
            rows.append(row)
        freeze["trace_manifest_sha256"] = digest.hexdigest()
        path = self.root / "trace.jsonl"

        def publish(values):
            with path.open("w", encoding="utf-8") as output:
                for value in values:
                    output.write(contract_api.canonical_json(value) + "\n")

        publish(rows)
        expected = path.read_bytes()
        subject._validate_bound_recovery_trace_bytes(
            self.contract, freeze, expected)
        rows[0]["trace_sha256"] = "0" * 64
        publish(rows)
        with self.assertRaisesRegex(
                subject.RecoveryRunnerError, "frozen trace hash"):
            subject._validate_bound_recovery_trace_bytes(
                self.contract, freeze, path.read_bytes())

    def test_restore_failure_cannot_publish_complete_run_summary(self) -> None:
        candidate_id = subject.CANDIDATE_SPECS[0][0]
        description = {
            "resolved_path": "/fixture/worker",
            "binary_sha256": "a" * 64,
            "source_git_commit": "1" * 40,
            "arms": [],
        }
        receipt = {
            "record_count": subject.RECOVERY_RECORDS,
            "freeze_manifest_sha256": "b" * 64,
            "result_stream_sha256": "c" * 64,
            "receipt_sha256": "d" * 64,
            "thermal": {
                "sample_count": 2,
                "cpu_tctl_max_millic": 65000,
                "dimm_max_millic": 42000,
            },
        }
        assembled = {"summary": {}, "execution_receipt": receipt}
        args = SimpleNamespace(
            candidate=candidate_id, deadline_seconds=60.0,
            contract=contract_api.DEFAULT_CONTRACT,
            worker=Path("worker"), cpus=None, controller_cpu=None,
            sampler_pid=1, sampler_cpu=127,
            sampler_script=Path("sampler.py"),
            sampler_csv=Path("thermal.csv"),
            output_dir=self.root / "never-published")
        cpus = list(range(8))
        source = "1" * 40
        with ExitStack() as stack:
            stack.enter_context(mock.patch.object(
                subject.contract_api, "load_contract",
                return_value=self.contract))
            stack.enter_context(mock.patch.object(
                subject, "describe_candidate_worker",
                return_value=description))
            stack.enter_context(mock.patch.object(
                subject.os, "sched_getaffinity", return_value=set(range(10))))
            stack.enter_context(mock.patch.object(
                subject.runner_api, "select_cpu_layout",
                return_value=(cpus, 8)))
            stack.enter_context(mock.patch.object(
                subject.runner_api, "_preflight_sampler"))
            stack.enter_context(mock.patch.object(
                subject.runner_api, "_git_head",
                side_effect=(source, source, source)))
            stack.enter_context(mock.patch.object(
                subject.runner_api, "_require_worker_source_commit"))
            stack.enter_context(mock.patch.object(
                subject.runner_api, "_create_output_dir",
                return_value=args.output_dir))
            stack.enter_context(mock.patch.object(
                subject.runner_api, "_emit_and_assemble_trace",
                return_value=(Path("trace.jsonl"), "e" * 64)))
            stack.enter_context(mock.patch.object(
                subject, "write_recovery_freeze", return_value={}))
            stack.enter_context(mock.patch.object(
                subject, "spawn_candidate_workers", return_value=[object()]))
            stack.enter_context(mock.patch.object(
                subject.runner_api, "_pin_controller"))
            stack.enter_context(mock.patch.object(
                subject.runner_api, "choose_new_sampler_start",
                return_value=100))
            stack.enter_context(mock.patch.object(
                subject, "_run_recovery_jobs",
                return_value=(Path("native.jsonl"), 200)))
            stack.enter_context(mock.patch.object(
                subject.runner_api, "_wait_for_sampler_sample",
                return_value=300))
            stack.enter_context(mock.patch.object(
                subject.native_api, "write_sampler_attestation"))
            stack.enter_context(mock.patch.object(
                subject.native_api, "assemble_results",
                return_value=assembled))
            stack.enter_context(mock.patch.object(
                subject.native_api, "validate_execution_receipt",
                return_value=assembled))
            stack.enter_context(mock.patch.object(
                subject.runner_api, "quit_workers"))
            stack.enter_context(mock.patch.object(
                subject.runner_api, "_restore_controller_affinity",
                side_effect=subject.runner_api.RunnerError(
                    "injected restore failure")))
            publish = stack.enter_context(mock.patch.object(
                subject.runner_api, "_atomic_write_object"))
            with self.assertRaisesRegex(
                    subject.runner_api.RunnerError, "restore failure"):
                subject.run_recovery_screen(args)
            publish.assert_not_called()

    def test_combination_summary_is_hash_bound_and_not_a_receipt(self) -> None:
        campaigns = self._campaigns()
        output = self.root / "combined"
        args = SimpleNamespace(
            contract=contract_api.DEFAULT_CONTRACT,
            campaign_dir=[Path("one"), Path("two"), Path("three")],
            work_rank_dir=Path("work-rank"),
            output_dir=output)
        work_binding = {
            "work_rank_summary_sha256": "1" * 64,
            "work_rank_result_stream_sha256": "2" * 64,
            "work_rank_domain_sha256": "3" * 64,
            "raw_identity_join_count": subject.RAW_IDENTITY_JOIN_COUNT,
            "raw_identity_join_sha256": "4" * 64,
        }

        def validate_ledger(contract, phase, result_path, freeze_path,
                            trace_path):
            freeze = contract_api.load_freeze_manifest(
                contract, phase, freeze_path, "recovery")
            return {
                "schema": "fixture.logical-summary.v1",
                "freeze_manifest_sha256":
                    contract_api.freeze_manifest_sha256(freeze),
                "rows": len(result_path.read_text().splitlines()),
            }

        with mock.patch.object(
                subject, "load_completed_campaign",
                side_effect=campaigns), \
                mock.patch.object(
                    subject.work_api, "load_completed_work_screen",
                    return_value={"summary": {}, "rows": []}), \
                mock.patch.object(
                    subject, "_bind_work_rank_identities",
                    return_value=work_binding), \
                mock.patch.object(
                    subject.contract_api, "validate_ledger",
                    side_effect=validate_ledger):
            summary = subject.combine_recovery_screens(args)
        self.assertEqual(set(summary), subject.COMBINATION_SUMMARY_FIELDS)
        self.assertEqual(summary["schema"],
                         subject.COMBINATION_SUMMARY_SCHEMA)
        self.assertNotEqual(summary["schema"], native_api.EXECUTION_SCHEMA)
        self.assertIs(summary["is_execution_receipt"], False)
        self.assertEqual(summary["artifact_kind"],
                         "logical_recovery_combination")
        self.assertEqual(
            summary["construction_seed_basis"], subject.RAW_SEED_BASIS)
        self.assertEqual(summary["seed_schedule_sha256"],
                         subject.RAW_SEED_SCHEDULE_SHA256)
        self.assertEqual(
            summary["candidate_roster"],
            [value[0] for value in subject.CANDIDATE_SPECS])
        self.assertEqual(summary["validator_summary"]["rows"],
                         subject.LOGICAL_RECORDS)
        self.assertEqual(
            summary["combination_sha256"],
            subject._self_hash(summary, "combination_sha256"))
        stored = subject._load_canonical_object(
            output / "logical-recovery-summary.json", "stored summary")
        self.assertEqual(stored, summary)
        self.assertFalse((output / "recovery-execution.json").exists())

    def test_main_requires_exactly_three_combination_directories(self) -> None:
        stderr = StringIO()
        with redirect_stderr(stderr):
            result = subject.main([
                "combine", "--campaign-dir", "one",
                "--campaign-dir", "two", "--work-rank-dir", "work",
                "--output-dir", "out"])
        self.assertEqual(result, 1)
        self.assertIn("exactly three", stderr.getvalue())

    def test_completed_artifact_preflight_rejects_fifo_and_oversize(self) \
            -> None:
        fifo = self.root / "artifact.fifo"
        os.mkfifo(str(fifo))
        started = time.monotonic()
        with self.assertRaisesRegex(
                subject.RecoveryRunnerError, "regular non-symlink"):
            subject._read_regular_bytes(fifo, "fixture FIFO")
        self.assertLess(time.monotonic() - started, 1.0)
        oversized = self.root / "oversized.bin"
        with oversized.open("wb") as output:
            output.truncate(subject.MAX_COMPLETED_ARTIFACT_BYTES + 1)
        with self.assertRaisesRegex(
                subject.RecoveryRunnerError, "bounded artifact size"):
            subject._read_regular_bytes(oversized, "oversized fixture")

        source = self.root / "pinned-source"
        source.mkdir()
        (source / "artifact").write_bytes(b"original")
        directory_fd = os.open(
            str(source), os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
        try:
            source.rename(self.root / "held-source")
            source.mkdir()
            (source / "artifact").write_bytes(b"substitute")
            self.assertEqual(
                subject._read_regular_bytes(
                    Path("artifact"), "pinned fixture", directory_fd),
                b"original")
        finally:
            os.close(directory_fd)

    def test_combine_hard_wall_interrupts_synchronous_validation(self) -> None:
        stderr = StringIO()
        started = time.monotonic()
        argv = [
            "combine", "--campaign-dir", "one", "--campaign-dir", "two",
            "--campaign-dir", "three", "--work-rank-dir", "work",
            "--output-dir", "out",
            "--deadline-seconds", "0.05",
        ]
        with redirect_stderr(stderr), mock.patch.object(
                subject, "combine_recovery_screens",
                side_effect=lambda _args: time.sleep(2.0)):
            result = subject.main(argv)
        self.assertEqual(result, 1)
        self.assertLess(time.monotonic() - started, 1.0)
        self.assertIn("hard wall expired", stderr.getvalue())


if __name__ == "__main__":
    unittest.main()
