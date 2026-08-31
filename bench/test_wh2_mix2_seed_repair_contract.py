#!/usr/bin/env python3
"""Focused pure tests for the isolated Two07/mix2 v9 repair contract."""

from __future__ import annotations

import copy
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest
from unittest import mock

from bench import wh2_mix2_seed_repair_contract as subject
from bench import wh2_mix2_promotion_short_screen as short_subject


def cell(K, attempt, ordinal, validation=False, success=True):
    root_index, root, schedule = subject._expected_cell(ordinal, validation)
    trace_sha256 = hashlib.sha256(
        "{}:{}:{}".format(K, int(validation), ordinal).encode(
            "ascii")).hexdigest()
    return {
        "attempted_candidates": K + 4 + ordinal,
        "cell_ordinal": ordinal,
        "construction_attempt": attempt,
        "decoded_extra": 0 if success else None,
        "loss_ppm": 500000,
        "loss_root": root,
        "outcome": "success" if success else "need_more_at_cap",
        "result_code": 0 if success else 1,
        "root_index": root_index,
        "schedule": schedule,
        "trace_sha256": trace_sha256,
    }


def derivation_record(contract, K=8, selected=2, worker="b" * 64):
    base_precode = 0x0123456789ABCDEF
    base_packet = 0x89ABCDEF
    return {
        "K": K,
        "base_packet_seed": "0x{:08x}".format(base_packet),
        "base_precode_seed": "0x{:016x}".format(base_precode),
        "candidate_profile_sha256":
            subject.candidate_profile_sha256(contract),
        "effective_packet_seed":
            subject._effective_packet_seed(base_packet, selected),
        "effective_precode_seed":
            subject._effective_precode_seed(base_precode, selected),
        "lower_attempt_failure_witnesses": [
            cell(K, attempt, attempt % subject.SELECTION_CELL_COUNT,
                 success=False)
            for attempt in range(selected)
        ],
        "mode": "derive",
        "ordinal": K - subject.K_MIN,
        "schema": subject.DERIVATION_RECORD_SCHEMA,
        "selected_attempt": selected,
        "selected_successes": [
            cell(K, selected, ordinal)
            for ordinal in range(subject.SELECTION_CELL_COUNT)
        ],
        "source_sha256": "c" * 64,
        "worker_binary_sha256": worker,
    }


def validation_record(
        contract, K=8, attempt=2, worker="b" * 64, failed=None):
    base_precode = 0x0123456789ABCDEF
    base_packet = 0x89ABCDEF
    cells = [
        cell(K, attempt, ordinal, validation=True,
             success=ordinal != failed)
        for ordinal in range(subject.VALIDATION_CELL_COUNT)
    ]
    return {
        "K": K,
        "all_oh0_success": failed is None,
        "base_packet_seed": "0x{:08x}".format(base_packet),
        "base_precode_seed": "0x{:016x}".format(base_precode),
        "candidate_profile_sha256":
            subject.candidate_profile_sha256(contract),
        "cells": cells,
        "construction_attempt": attempt,
        "effective_packet_seed":
            subject._effective_packet_seed(base_packet, attempt),
        "effective_precode_seed":
            subject._effective_precode_seed(base_precode, attempt),
        "mode": "validate",
        "ordinal": K - subject.K_MIN,
        "schema": subject.VALIDATION_RECORD_SCHEMA,
        "source_sha256": "c" * 64,
        "worker_binary_sha256": worker,
    }


def canonical_bytes(value):
    return (subject.canonical_json(value) + "\n").encode("ascii")


def reseal(value, field):
    unsigned = dict(value)
    unsigned.pop(field, None)
    value[field] = subject.sha256_json(unsigned)


class ContractIdentityTests(unittest.TestCase):
    def setUp(self):
        self.contract = subject.load_contract()

    def test_closed_v9_identity_and_cardinality(self):
        self.assertEqual(self.contract["schema"], subject.CONTRACT_SCHEMA)
        self.assertEqual(subject.MAP_BYTES, 63999)
        self.assertEqual(
            subject.candidate_profile_sha256(self.contract),
            "90233b44a0893f96c1a18c19aa61ada052c935a48c6bf7d6a2813065856651f0")
        self.assertEqual(self.contract["candidate_profile"]["mix_count"], 2)
        self.assertEqual(
            self.contract["candidate_profile"]["graph_seed_block_bytes"], 2)
        self.assertEqual(
            self.contract["candidate_profile"]["construction_seed_basis"],
            subject.PRODUCTION_SEED_BASIS)
        self.assertEqual(
            self.contract["candidate_profile"]["seed_schedule_sha256"],
            subject.PRODUCTION_SEED_SCHEDULE_SHA256)
        self.assertEqual(
            hashlib.sha256(
                subject.PRODUCTION_SEED_SCHEDULE_CANONICAL.encode(
                    "ascii")).hexdigest(),
            subject.PRODUCTION_SEED_SCHEDULE_SHA256)
        self.assertEqual(
            self.contract["candidate_profile"]["dense_anchor_layout"],
            "two07")
        self.assertEqual(
            self.contract["repair_map"]["bundle_chain"],
            subject.BUNDLE_CHAIN_RULE)
        self.assertFalse(self.contract["repair_map"]["runtime_search"])
        self.assertTrue(
            set(subject.SELECTION_ROOTS).isdisjoint(subject.VALIDATION_ROOTS))
        self.assertEqual(subject.SELECTION_CELL_COUNT, 27)
        self.assertEqual(subject.VALIDATION_CELL_COUNT, 9)

    def test_v9_final_roots_recompute_from_full_frozen_digests(self):
        derived = []
        for index, expected in enumerate(subject.VALIDATION_ROOT_FULL_SHA256):
            digest = hashlib.sha256(
                (subject.VALIDATION_ROOT_NAMESPACE + str(index)).encode(
                    "ascii")).hexdigest()
            self.assertEqual(digest, expected)
            derived.append("0x" + digest[:16])
        self.assertEqual(tuple(derived), subject.VALIDATION_ROOTS)
        self.assertEqual(
            subject.VALIDATION_ROOT_NAMESPACE,
            "wirehair2-two07-mix2-graph-b2-all-k-v9:holdout-root:")
        self.assertEqual(len(set(subject.SELECTION_ROOTS)), 9)
        self.assertEqual(len(set(subject.VALIDATION_ROOTS)), 3)
        self.assertTrue(
            set(subject.SELECTION_ROOTS).isdisjoint(subject.VALIDATION_ROOTS))
        self.assertNotIn(0, [int(value, 16)
                             for value in subject.VALIDATION_ROOTS])

    def test_validation_roster_digest_binds_order_and_cardinality(self):
        self.assertEqual(
            subject.validation_roster_sha256(),
            subject.VALIDATION_ROSTER_SHA256)
        self.assertEqual(
            self.contract["validation"]["roster"],
            subject.EXPECTED_VALIDATION_ROSTER)
        permuted_roots = list(subject.VALIDATION_ROOTS)
        permuted_roots[0], permuted_roots[1] = (
            permuted_roots[1], permuted_roots[0])
        permuted_schedules = list(subject.SCHEDULES)
        permuted_schedules[0], permuted_schedules[1] = (
            permuted_schedules[1], permuted_schedules[0])
        self.assertNotEqual(
            subject.validation_roster_sha256(permuted_roots),
            subject.VALIDATION_ROSTER_SHA256)
        self.assertNotEqual(
            subject.validation_roster_sha256(
                subject.VALIDATION_ROOTS, permuted_schedules),
            subject.VALIDATION_ROSTER_SHA256)
        self.assertNotEqual(
            subject.validation_roster_sha256(
                subject.VALIDATION_ROOTS[:-1], subject.SCHEDULES),
            subject.VALIDATION_ROSTER_SHA256)

    def test_worker_description_roster_tamper_fails_preflight(self):
        worker_sha256 = "b" * 64
        source_commit = "d" * 40
        description = {
            "binary_sha256": worker_sha256,
            "candidate_profile": self.contract["candidate_profile"],
            "candidate_profile_sha256":
                subject.candidate_profile_sha256(self.contract),
            "contract_schema": subject.CONTRACT_SCHEMA,
            "derivation_schema": subject.DERIVATION_RECORD_SCHEMA,
            "protocol": "D ordinal K | V ordinal K attempt | Q",
            "schema": subject.DESCRIPTION_SCHEMA,
            "source_git_commit": source_commit,
            "validation_roster_schema": subject.VALIDATION_ROSTER_SCHEMA,
            "validation_roster_sha256": subject.VALIDATION_ROSTER_SHA256,
            "validation_schema": subject.VALIDATION_RECORD_SCHEMA,
            "worker_schema": subject.WORKER_SCHEMA,
        }
        subject._validate_worker_description_identity(
            description, worker_sha256, source_commit, self.contract)
        for field, value in (
                ("validation_roster_schema", "retired-roster.v0"),
                ("validation_roster_sha256", "f" * 64)):
            with self.subTest(field=field):
                broken = copy.deepcopy(description)
                broken[field] = value
                with self.assertRaises(subject.ContractError):
                    subject._validate_worker_description_identity(
                        broken, worker_sha256, source_commit, self.contract)

    def test_worker_description_candidate_profile_types_are_strict(self):
        worker_sha256 = "b" * 64
        source_commit = "d" * 40
        description = {
            "binary_sha256": worker_sha256,
            "candidate_profile": copy.deepcopy(
                self.contract["candidate_profile"]),
            "candidate_profile_sha256":
                subject.candidate_profile_sha256(self.contract),
            "contract_schema": subject.CONTRACT_SCHEMA,
            "derivation_schema": subject.DERIVATION_RECORD_SCHEMA,
            "protocol": "D ordinal K | V ordinal K attempt | Q",
            "schema": subject.DESCRIPTION_SCHEMA,
            "source_git_commit": source_commit,
            "validation_roster_schema": subject.VALIDATION_ROSTER_SCHEMA,
            "validation_roster_sha256": subject.VALIDATION_ROSTER_SHA256,
            "validation_schema": subject.VALIDATION_RECORD_SCHEMA,
            "worker_schema": subject.WORKER_SCHEMA,
        }
        for field, alias in (
                ("dense_identity_corner", 0),
                ("binary_dense_rows", 12.0)):
            with self.subTest(field=field):
                broken = copy.deepcopy(description)
                broken["candidate_profile"][field] = alias
                with self.assertRaises(subject.ContractError):
                    subject._validate_worker_description_identity(
                        broken, worker_sha256, source_commit, self.contract)

    def test_retired_contracts_and_artifact_schemas_fail_closed(self):
        for name in ("wh2_mix2_seed_repair_contract_v5.json",
                     "wh2_mix2_seed_repair_contract_v7.json",
                     "wh2_mix2_seed_repair_contract_v8.json"):
            with self.subTest(name=name), self.assertRaises(
                    subject.ContractError):
                subject.load_contract(Path(subject.__file__).with_name(name))
        self.assertTrue(subject.CONTRACT_SCHEMA.endswith(".v9"))
        for schema in (
                subject.DESCRIPTION_SCHEMA, subject.WORKER_SCHEMA,
                subject.DERIVATION_RECORD_SCHEMA,
                subject.VALIDATION_RECORD_SCHEMA):
            self.assertTrue(schema.endswith(".v3"))
        for schema in (
                subject.MANIFEST_SCHEMA, subject.FREEZE_SCHEMA,
                subject.MAP_SCHEMA, subject.DERIVATION_SUMMARY_SCHEMA,
                subject.VALIDATION_SUMMARY_SCHEMA,
                subject.MAP_EXPORT_SCHEMA,
                subject.SEMANTIC_REPLAY_SCHEMA):
            self.assertTrue(schema.endswith(".v2"))

    def test_short_screen_contract_is_the_frozen_v9_prerequisite(self):
        screen = short_subject.contract_description()
        self.assertEqual(
            screen["contract_sha256"], subject.SHORT_SCREEN_CONTRACT_SHA256)
        self.assertEqual(
            self.contract["short_screen"]["contract_sha256"],
            subject.SHORT_SCREEN_CONTRACT_SHA256)
        self.assertFalse(
            self.contract["short_screen"]["final_holdout_access"])

    def test_contract_rejects_mix3_or_validation_root_reuse(self):
        for mutation in (
                "mix", "root", "validation_rule", "applicability",
                "extra_repair_field", "bool_overhead", "float_attempts",
                "integer_runtime_search", "validation_roster_digest",
                "validation_roster_count"):
            value = copy.deepcopy(self.contract)
            if mutation == "mix":
                value["candidate_profile"]["mix_count"] = 3
            elif mutation in ("root", "validation_rule"):
                if mutation == "root":
                    value["validation"]["roots"][0] = \
                        subject.SELECTION_ROOTS[0]
                else:
                    value["validation"]["rule"] = "runtime retry is allowed"
            elif mutation == "applicability":
                value["applicability"]["forbidden_inference"] = \
                    "all profiles are equivalent"
            elif mutation == "bool_overhead":
                value["domain"]["overhead"] = False
            elif mutation == "float_attempts":
                value["attempt_schedule"]["attempt_count"] = 256.0
            elif mutation == "integer_runtime_search":
                value["repair_map"]["runtime_search"] = 0
            elif mutation == "validation_roster_digest":
                value["validation"]["roster"]["sha256"] = "f" * 64
            elif mutation == "validation_roster_count":
                value["validation"]["roster"]["cell_count"] = 8
            else:
                value["repair_map"]["runtime_retry"] = False
            with self.subTest(mutation=mutation), \
                    tempfile.TemporaryDirectory() as temporary:
                path = Path(temporary) / "contract.json"
                path.write_text(json.dumps(value), encoding="utf-8")
                with self.assertRaises(subject.ContractError):
                    subject.load_contract(path)

    def test_K_specific_attempt_step_goldens_include_uint8_endpoint(self):
        base_precode = 0x0123456789ABCDEF
        base_packet = 0x89ABCDEF
        self.assertEqual(
            [subject._effective_precode_seed(base_precode, value)
             for value in (0, 1, 2, 255)],
            ["0x0123456789abcdef", "0x9f5abf2108f64a04",
             "0x3d9238da8840c619", "0x9a65852d54dd66da"])
        self.assertEqual(
            [subject._effective_packet_seed(base_packet, value)
             for value in (0, 1, 2, 255)],
            ["0x89abcdef", "0x27e347a8", "0xc61ac161", "0x22ee0d36"])


class DerivationRecordTests(unittest.TestCase):
    def setUp(self):
        self.contract = subject.load_contract()
        self.worker = "b" * 64

    def test_lowest_attempt_proof_accepts_exact_witnesses_and_27_passes(self):
        value = derivation_record(self.contract, selected=2, worker=self.worker)
        self.assertEqual(subject.validate_derivation_record(
            value, self.contract, self.worker, expected_K=8), 2)

    def test_missing_or_successful_lower_witness_is_invalid(self):
        value = derivation_record(self.contract, selected=2, worker=self.worker)
        value["lower_attempt_failure_witnesses"].pop()
        with self.assertRaises(subject.ContractError):
            subject.validate_derivation_record(value, self.contract, self.worker)
        value = derivation_record(self.contract, selected=2, worker=self.worker)
        value["lower_attempt_failure_witnesses"][0] = cell(8, 0, 0)
        with self.assertRaises(subject.ContractError):
            subject.validate_derivation_record(value, self.contract, self.worker)

    def test_selected_attempt_requires_all_27_frozen_cells(self):
        value = derivation_record(self.contract, selected=1, worker=self.worker)
        value["selected_successes"][26] = cell(8, 1, 26, success=False)
        with self.assertRaises(subject.ContractError):
            subject.validate_derivation_record(value, self.contract, self.worker)

    def test_lower_witness_reuses_attempt_independent_trace_identity(self):
        value = derivation_record(self.contract, selected=2, worker=self.worker)
        value["lower_attempt_failure_witnesses"][0]["trace_sha256"] = \
            "e" * 64
        with self.assertRaises(subject.ContractError):
            subject.validate_derivation_record(
                value, self.contract, self.worker)
        value = derivation_record(self.contract, selected=2, worker=self.worker)
        value["lower_attempt_failure_witnesses"][1][
            "attempted_candidates"] += 1
        with self.assertRaises(subject.ContractError):
            subject.validate_derivation_record(
                value, self.contract, self.worker)
        value = derivation_record(self.contract, selected=1, worker=self.worker)
        value["selected_successes"][4]["loss_root"] = subject.VALIDATION_ROOTS[0]
        with self.assertRaises(subject.ContractError):
            subject.validate_derivation_record(value, self.contract, self.worker)

    def test_boolean_attempt_alias_is_invalid(self):
        value = derivation_record(self.contract, selected=0, worker=self.worker)
        value["selected_attempt"] = False
        with self.assertRaises(subject.ContractError):
            subject.validate_derivation_record(value, self.contract, self.worker)

    def test_loss_ppm_rejects_bool_and_float_aliases(self):
        for replacement in (False, 500000.0):
            value = derivation_record(
                self.contract, selected=0, worker=self.worker)
            value["selected_successes"][0]["loss_ppm"] = replacement
            with self.subTest(replacement=replacement), \
                    self.assertRaises(subject.ContractError):
                subject.validate_derivation_record(
                    value, self.contract, self.worker)
        value = derivation_record(self.contract, K=2, selected=0,
                                  worker=self.worker)
        value["ordinal"] = False
        with self.assertRaises(subject.ContractError):
            subject.validate_derivation_record(value, self.contract, self.worker)
        value = derivation_record(self.contract, selected=0, worker=self.worker)
        value["selected_successes"][0]["construction_attempt"] = False
        with self.assertRaises(subject.ContractError):
            subject.validate_derivation_record(value, self.contract, self.worker)


class ValidationRecordTests(unittest.TestCase):
    def setUp(self):
        self.contract = subject.load_contract()
        self.worker = "b" * 64

    def test_disjoint_validation_pass_and_scored_failure(self):
        passed = validation_record(self.contract, worker=self.worker)
        self.assertTrue(subject.validate_validation_record(
            passed, self.contract, self.worker, 8, 2))
        failed = validation_record(
            self.contract, worker=self.worker, failed=4)
        self.assertFalse(subject.validate_validation_record(
            failed, self.contract, self.worker, 8, 2))

    def test_attempt_substitution_and_training_root_are_invalid(self):
        value = validation_record(self.contract, worker=self.worker)
        with self.assertRaises(subject.ContractError):
            subject.validate_validation_record(
                value, self.contract, self.worker, 8, 1)
        value = validation_record(self.contract, worker=self.worker)
        value["cells"][0]["loss_root"] = subject.SELECTION_ROOTS[0]
        with self.assertRaises(subject.ContractError):
            subject.validate_validation_record(
                value, self.contract, self.worker, 8, 2)

    def test_aggregate_cannot_hide_cell_failure(self):
        value = validation_record(
            self.contract, worker=self.worker, failed=3)
        value["all_oh0_success"] = True
        with self.assertRaises(subject.ContractError):
            subject.validate_validation_record(
                value, self.contract, self.worker, 8, 2)


class ShortScreenBundleTests(unittest.TestCase):
    def setUp(self):
        self.contract = subject.load_contract()
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.worker_sha256 = "b" * 64
        self.source_commit = "d" * 40
        self.description = {
            "binary_sha256": self.worker_sha256,
            "candidate_profile": self.contract["candidate_profile"],
            "candidate_profile_sha256":
                subject.candidate_profile_sha256(self.contract),
            "contract_schema": subject.CONTRACT_SCHEMA,
            "derivation_schema": subject.DERIVATION_RECORD_SCHEMA,
            "protocol": "D ordinal K | V ordinal K attempt | Q",
            "schema": subject.DESCRIPTION_SCHEMA,
            "source_git_commit": self.source_commit,
            "validation_roster_schema": subject.VALIDATION_ROSTER_SCHEMA,
            "validation_roster_sha256": subject.VALIDATION_ROSTER_SHA256,
            "validation_schema": subject.VALIDATION_RECORD_SCHEMA,
            "worker_schema": subject.WORKER_SCHEMA,
        }
        self.worker = subject.VerifiedWorker(
            Path("/unused"), -1, self.worker_sha256,
            self.source_commit, self.description,
            subject.VALIDATION_ROSTER_SHA256)
        self.source_receipt = {
            "source_git_commit": self.source_commit,
            "tracked_and_untracked_linked_sources_clean": True,
            "clean_status_scope":
                list(subject.SHORT_SCREEN_CLEAN_STATUS_SCOPE),
            "source_files": {
                path: "c" * 64 for path in subject.SHORT_SCREEN_SOURCE_FILES
            },
        }
        reseal(self.source_receipt, "source_receipt_sha256")
        self.screen_contract = short_subject.contract_description()
        self.attempt_records = [
            derivation_record(
                self.contract, K=K, selected=0, worker=self.worker_sha256)
            for K in subject.SEMANTIC_REPLAY_K
        ]
        self.attempt_bytes = b"".join(
            canonical_bytes(record) for record in self.attempt_records)
        self.result_records = self._make_result_records()
        self.result_bytes = b"".join(
            canonical_bytes(record) for record in self.result_records)
        self._make_objects()
        self._write_bundle()

    def _make_result_records(self):
        derivation_by_K = {
            record["K"]: record for record in self.attempt_records
        }
        invocations = short_subject.make_invocations({
            K: 0 for K in subject.SEMANTIC_REPLAY_K
        })
        records = []
        for invocation in invocations:
            candidate = invocation.arm == "candidate_two07_mix2"
            derivation = derivation_by_K[invocation.K]
            staircase = 346 if invocation.K == 64000 else 2
            entries = short_subject.extra_dense_basis_capacity_entries(
                invocation.K, staircase)
            identity = invocation.identity(Path("/unused/bench"))
            records.append({
                "schema": subject.SHORT_SCREEN_RESULT_SCHEMA,
                "ordinal": invocation.ordinal,
                "coordinate_ordinal": invocation.coordinate_ordinal,
                "arm": invocation.arm,
                "K": invocation.K,
                "block_bytes": short_subject.BLOCK_BYTES,
                "loss_ppm": short_subject.LOSS_PPM,
                "overhead": short_subject.OVERHEAD,
                "root_index": invocation.root_index,
                "loss_root": invocation.root,
                "schedule": invocation.schedule,
                "cell_ordinal": invocation.cell_ordinal,
                "timing_order": invocation.timing_order,
                "timing_slot": invocation.timing_slot,
                "observation_index": invocation.observation_index,
                "attempt_selection_stream_sha256":
                    hashlib.sha256(self.attempt_bytes).hexdigest(),
                "attempt_selection_cell_receipts_used_as_promotion_evidence":
                    False,
                "benchmark_loss_trace_hash_recorded": False,
                "full_payload_byte_recovery_verified": True,
                "candidate_profile_sha256":
                    subject.candidate_profile_sha256(self.contract),
                "construction_attempt": 0,
                "attempt_mode": "exact" if candidate else "selected",
                "effective_precode_seed":
                    derivation["effective_precode_seed"],
                "effective_packet_seed":
                    derivation["effective_packet_seed"],
                "staircase": staircase,
                "binary_dense_rows": 12,
                "gf256_heavy_rows": 12,
                "source_hits": 2,
                "dense_anchor_layout": "two07" if candidate else "disabled",
                "mix_count": 2 if candidate else 3,
                "success": 1,
                "rank_fail": 0,
                "error": 0,
                "weak": False,
                "inactivated_columns": 100,
                "block_xors": 98 if candidate else 100,
                "gf256_muladds": 100,
                "solve_ms": "1.000",
                "build_ms": "1.000",
                "peel_ms": "1.000",
                "project_ms": "1.000",
                "residual_ms": "1.000",
                "backsub_ms": "1.000",
                "extra_dense_basis_capacity_entries":
                    entries if candidate else 0,
                "extra_dense_basis_capacity_bytes":
                    entries * 4 if candidate else 0,
                "command_sha256": identity["command_sha256"],
                "bench_stdout_sha256": hashlib.sha256(
                    str(invocation.ordinal).encode("ascii")).hexdigest(),
                "bench_binary_sha256": "a" * 64,
                "bench_source_git_commit": self.source_commit,
                "source_receipt_sha256":
                    self.source_receipt["source_receipt_sha256"],
            })
        return records

    def _make_objects(self):
        worker_binary = {
            "path": "/unused/worker",
            "sha256": self.worker_sha256,
            "size": 1,
        }
        bench_binary = {
            "path": "/unused/bench",
            "sha256": "a" * 64,
            "size": 1,
        }
        self.input_receipt = {
            "schema": subject.SHORT_SCREEN_INPUT_SCHEMA,
            "contract": self.screen_contract,
            "contract_sha256": subject.SHORT_SCREEN_CONTRACT_SHA256,
            "source_receipt": self.source_receipt,
            "source_receipt_sha256":
                self.source_receipt["source_receipt_sha256"],
            "bench_binary": bench_binary,
            "repair_worker_binary": worker_binary,
            "repair_worker_description": self.description,
            "repair_worker_description_sha256":
                subject.sha256_json(self.description),
            "attempt_selection_worker_command_sha256": hashlib.sha256(
                ("".join(
                    "D {} {}\n".format(K - subject.K_MIN, K)
                    for K in subject.SEMANTIC_REPLAY_K) + "Q\n").encode(
                        "ascii")).hexdigest(),
            "attempt_selection_stream_schema":
                subject.SHORT_SCREEN_ATTEMPT_STREAM_SCHEMA,
            "attempt_selection_stream_sha256":
                hashlib.sha256(self.attempt_bytes).hexdigest(),
            "attempt_selection_record_count": 30,
            "selected_attempts": [0] * 30,
            "worker_validation_commands_issued": 0,
            "worker_final_validation_stream_present": False,
            "invocations": [
                invocation.identity(Path("/unused/bench"))
                for invocation in short_subject.make_invocations({
                    K: 0 for K in subject.SEMANTIC_REPLAY_K
                })
            ],
        }
        reseal(self.input_receipt, "input_sha256")
        input_bytes = canonical_bytes(self.input_receipt)
        self.summary = {
            "schema": subject.SHORT_SCREEN_SUMMARY_SCHEMA,
            "contract_sha256": subject.SHORT_SCREEN_CONTRACT_SHA256,
            "input_sha256": self.input_receipt["input_sha256"],
            "input_file_sha256": hashlib.sha256(input_bytes).hexdigest(),
            "attempt_selection_stream_schema":
                subject.SHORT_SCREEN_ATTEMPT_STREAM_SCHEMA,
            "attempt_selection_stream_sha256":
                hashlib.sha256(self.attempt_bytes).hexdigest(),
            "attempt_selection_record_count": 30,
            "worker_validation_commands_issued": 0,
            "worker_final_validation_stream_present": False,
            "result_stream_sha256":
                hashlib.sha256(self.result_bytes).hexdigest(),
            "record_count": 1080,
            "invocation_count": 1080,
            "bench_binary_sha256": "a" * 64,
            "repair_worker_binary_sha256": self.worker_sha256,
            "source_git_commit": self.source_commit,
            "source_receipt": self.source_receipt,
            "source_receipt_sha256":
                self.source_receipt["source_receipt_sha256"],
            "candidate_profile_sha256":
                subject.candidate_profile_sha256(self.contract),
            "architecture_selection_performed": False,
            "offline_attempt_derivation_performed": True,
            "official_scope": subject.SHORT_SCREEN_OFFICIAL_SCOPE,
            "all_K_recovery_claimed": False,
            "wirehair1_timing_claimed": False,
            **short_subject.adjudicate(
                self.result_records,
                hashlib.sha256(self.attempt_bytes).hexdigest()),
        }
        reseal(self.summary, "summary_sha256")

    def _write_bundle(self):
        (self.root / subject.SHORT_SCREEN_INPUT_NAME).write_bytes(
            canonical_bytes(self.input_receipt))
        (self.root / subject.SHORT_SCREEN_ATTEMPT_NAME).write_bytes(
            self.attempt_bytes)
        (self.root / subject.SHORT_SCREEN_RESULT_NAME).write_bytes(
            self.result_bytes)
        (self.root / subject.SHORT_SCREEN_SUMMARY_NAME).write_bytes(
            canonical_bytes(self.summary))

    def _reseal_input_and_summary(self):
        reseal(self.input_receipt, "input_sha256")
        self.summary["input_sha256"] = self.input_receipt["input_sha256"]
        self.summary["input_file_sha256"] = hashlib.sha256(
            canonical_bytes(self.input_receipt)).hexdigest()
        self.summary["attempt_selection_stream_sha256"] = hashlib.sha256(
            self.attempt_bytes).hexdigest()
        reseal(self.summary, "summary_sha256")
        self._write_bundle()

    def _reseal_results_and_summary(self):
        self.result_bytes = b"".join(
            canonical_bytes(record) for record in self.result_records)
        self.summary["result_stream_sha256"] = hashlib.sha256(
            self.result_bytes).hexdigest()
        reseal(self.summary, "summary_sha256")
        self._write_bundle()

    def _load(self, expected=None):
        return subject.load_short_screen_bundle(
            self.contract, self.root, self.worker,
            self.summary["summary_sha256"] if expected is None else expected)

    def test_complete_short_screen_pass_authenticates_30_attempts(self):
        verified = self._load()
        self.assertEqual(verified.attempt_map, bytes(30))
        self.assertEqual(
            tuple(K for K, unused in verified.derivation_records),
            subject.SEMANTIC_REPLAY_K)
        self.assertEqual(
            verified.binding["short_screen_map_sha256"],
            hashlib.sha256(bytes(30)).hexdigest())
        self.assertEqual(
            verified.binding["short_screen_summary_sha256"],
            self.summary["summary_sha256"])

    def test_external_summary_hash_and_pass_gates_are_mandatory(self):
        with self.assertRaises(subject.ContractError):
            self._load("f" * 64)
        gate = next(iter(subject.SHORT_SCREEN_GATE_FIELDS))
        self.summary["gates"][gate] = False
        reseal(self.summary, "summary_sha256")
        self._write_bundle()
        with self.assertRaises(subject.ContractError):
            self._load()

    def test_resealed_attempt_reordering_is_rejected_semantically(self):
        lines = self.attempt_bytes.splitlines(keepends=True)
        lines[0], lines[1] = lines[1], lines[0]
        self.attempt_bytes = b"".join(lines)
        self.input_receipt["attempt_selection_stream_sha256"] = \
            hashlib.sha256(self.attempt_bytes).hexdigest()
        self._reseal_input_and_summary()
        with self.assertRaises(subject.ContractError):
            self._load()

    def test_bool_alias_cannot_claim_zero_validation_commands(self):
        self.input_receipt["worker_validation_commands_issued"] = False
        self._reseal_input_and_summary()
        with self.assertRaises(subject.ContractError):
            self._load()

    def test_resealed_wrong_derivation_command_hash_is_rejected(self):
        self.input_receipt["attempt_selection_worker_command_sha256"] = \
            "f" * 64
        self._reseal_input_and_summary()
        with self.assertRaises(subject.ContractError):
            self._load()

    def test_resealed_partial_all_true_gate_set_is_rejected(self):
        self.summary["gates"].pop(next(iter(self.summary["gates"])))
        reseal(self.summary, "summary_sha256")
        self._write_bundle()
        with self.assertRaises(subject.ContractError):
            self._load()

    def test_binary_and_source_receipts_have_closed_strict_types(self):
        self.input_receipt["repair_worker_binary"]["size"] = True
        self._reseal_input_and_summary()
        with self.assertRaises(subject.ContractError):
            self._load()

    def test_resealed_worker_description_roster_tamper_is_rejected(self):
        embedded_description = copy.deepcopy(self.description)
        embedded_description["validation_roster_sha256"] = "f" * 64
        self.input_receipt["repair_worker_description"] = \
            embedded_description
        self.input_receipt["repair_worker_description_sha256"] = \
            subject.sha256_json(embedded_description)
        self._reseal_input_and_summary()
        with self.assertRaises(subject.ContractError):
            self._load()

    def test_resealed_worker_description_profile_type_alias_is_rejected(self):
        embedded_description = copy.deepcopy(self.description)
        embedded_description["candidate_profile"][
            "dense_identity_corner"] = 0
        self.input_receipt["repair_worker_description"] = \
            embedded_description
        self.input_receipt["repair_worker_description_sha256"] = \
            subject.sha256_json(embedded_description)
        self._reseal_input_and_summary()
        with self.assertRaises(subject.ContractError):
            self._load()

    def test_source_file_hash_bool_alias_is_rejected_cleanly(self):
        path = subject.SHORT_SCREEN_SOURCE_FILES[0]
        self.source_receipt["source_files"][path] = True
        reseal(self.source_receipt, "source_receipt_sha256")
        self.input_receipt["source_receipt_sha256"] = \
            self.source_receipt["source_receipt_sha256"]
        self.summary["source_receipt_sha256"] = \
            self.source_receipt["source_receipt_sha256"]
        self._reseal_input_and_summary()
        with self.assertRaises(subject.ContractError):
            self._load()

    def test_source_receipt_field_set_is_closed(self):
        self.source_receipt["unexpected"] = "resealed"
        reseal(self.source_receipt, "source_receipt_sha256")
        self.input_receipt["source_receipt_sha256"] = \
            self.source_receipt["source_receipt_sha256"]
        self.summary["source_receipt_sha256"] = \
            self.source_receipt["source_receipt_sha256"]
        self._reseal_input_and_summary()
        with self.assertRaises(subject.ContractError):
            self._load()

    def test_resealed_result_root_and_attempt_binding_tamper_is_rejected(self):
        for field, value in (
                ("loss_root", subject.SELECTION_ROOTS[0]),
                ("attempt_selection_stream_sha256", "f" * 64)):
            with self.subTest(field=field):
                original = self.result_records[0][field]
                self.result_records[0][field] = value
                self._reseal_results_and_summary()
                with self.assertRaises(subject.ContractError):
                    self._load()
                self.result_records[0][field] = original
                self._reseal_results_and_summary()

    def test_fabricated_all_true_performance_verdict_is_rejected(self):
        candidate = next(
            record for record in self.result_records
            if record["arm"] == "candidate_two07_mix2")
        candidate["block_xors"] = 10000
        pair = next(
            record for record in self.result_records
            if record["arm"] == candidate["arm"] and
            record["coordinate_ordinal"] == candidate["coordinate_ordinal"] and
            record is not candidate)
        pair["block_xors"] = 10000
        self._reseal_results_and_summary()
        with self.assertRaises(subject.ContractError):
            self._load()

    def test_result_outcome_work_and_timing_types_are_strict(self):
        record = self.result_records[0]
        for field, value in (
                ("success", 0),
                ("block_xors", -1),
                ("gf256_muladds", True),
                ("solve_ms", "-1.000")):
            with self.subTest(field=field):
                original = record[field]
                record[field] = value
                self._reseal_results_and_summary()
                with self.assertRaises(subject.ContractError):
                    self._load()
                record[field] = original
                self._reseal_results_and_summary()

    def test_symlinked_short_screen_artifact_is_rejected(self):
        path = self.root / subject.SHORT_SCREEN_ATTEMPT_NAME
        target = self.root / "saved-attempts.jsonl"
        os.replace(path, target)
        path.symlink_to(target.name)
        with self.assertRaises(subject.ContractError):
            self._load()


class FreezeAndPublicationTests(unittest.TestCase):
    def test_K_partition_is_complete_and_ordered(self):
        ranges = subject._split_ranges(7)
        flattened = [K for low, high in ranges for K in range(low, high + 1)]
        self.assertEqual(flattened[0], subject.K_MIN)
        self.assertEqual(flattened[-1], subject.K_MAX)
        self.assertEqual(len(flattened), subject.MAP_BYTES)
        self.assertEqual(len(set(flattened)), subject.MAP_BYTES)

    def test_bundle_publication_is_no_clobber(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            subject._publish_sources(root, {"a": b"first\n", "b": b"two\n"})
            with self.assertRaises(subject.ContractError):
                subject._publish_sources(root, {"a": b"replace\n"})
            self.assertEqual((root / "a").read_bytes(), b"first\n")

    def test_publication_write_failure_removes_private_stage(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with mock.patch.object(
                    subject.os, "write", side_effect=OSError("injected")):
                with self.assertRaises(OSError):
                    subject._publish_sources(root, {"artifact": b"value\n"})
            self.assertEqual(list(root.iterdir()), [])

    def test_clean_source_scope_rejects_dirty_fixup_tables(self):
        commit = "d" * 40
        for path in ("WirehairDenseFixups.inc", "WirehairPeelFixups.inc"):
            completed = (
                subprocess.CompletedProcess([], 0, (commit + "\n").encode(), b""),
                subprocess.CompletedProcess(
                    [], 0, (" M " + path + "\n").encode(), b""),
                subprocess.CompletedProcess([], 0, (commit + "\n").encode(), b""),
            )
            with self.subTest(path=path), \
                    mock.patch.object(
                        subject.subprocess, "run", side_effect=completed) \
                    as run, self.assertRaises(subject.ContractError):
                subject._require_clean_source(commit)
            status_command = run.call_args_list[1].args[0]
            self.assertIn(path, status_command)

    def test_worker_stderr_uses_a_nonblocking_bounded_file(self):
        with tempfile.TemporaryFile(mode="w+b") as stream:
            stream.write(b"diagnostic")
            stream.flush()
            self.assertEqual(
                subject._read_stderr_file(stream, "test"), b"diagnostic")
        with tempfile.TemporaryFile(mode="w+b") as stream:
            stream.write(b"x" * (subject.MAX_STDERR_BYTES + 1))
            stream.flush()
            with self.assertRaises(subject.ContractError):
                subject._read_stderr_file(stream, "test")

    def test_stdout_drain_rejects_more_than_pipe_capacity_without_deadlock(self):
        process = subprocess.Popen(
            [sys.executable, "-c",
             "import sys; sys.stdout.buffer.write(b'x' * 1048576); "
             "sys.stdout.buffer.flush()"],
            stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE, start_new_session=True)
        try:
            with self.assertRaises(subject.ContractError):
                subject._drain_worker_stdout(
                    process, time.monotonic() + 5, "malicious worker")
        finally:
            subject._kill_process(process)
            for stream in (process.stdout, process.stderr):
                if stream is not None:
                    stream.close()

    def test_bad_validation_roster_cannot_spawn_any_worker(self):
        contract = subject.load_contract()
        worker = subject.VerifiedWorker(
            Path("/unused"), -1, "b" * 64, "d" * 40, {}, "f" * 64)
        with tempfile.TemporaryDirectory() as temporary, \
                mock.patch.object(subject.subprocess, "Popen") as spawn:
            with self.assertRaises(subject.ContractError):
                subject._run_shard(
                    "validate", 2, 2, Path(temporary) / "shard.jsonl",
                    bytes(subject.MAP_BYTES), contract, worker,
                    time.monotonic() + 5, subject._Registry())
            with self.assertRaises(subject.ContractError):
                subject._replay_semantic_sample(worker, b"", (), ())
        spawn.assert_not_called()


class SemanticReplayTests(unittest.TestCase):
    class FakeProcess:
        def __init__(self, trailing_stdout=b""):
            self.stdin = io.BytesIO()
            read_fd, write_fd = os.pipe()
            if trailing_stdout:
                os.write(write_fd, trailing_stdout)
            os.close(write_fd)
            self.stdout = os.fdopen(read_fd, "rb")
            self.stderr = io.BytesIO()
            self.returncode = 0
            self.pid = 12345

        def wait(self, timeout=None):
            del timeout
            return self.returncode

        def poll(self):
            return self.returncode

    def setUp(self):
        self.contract = subject.load_contract()
        self.worker = subject.VerifiedWorker(
            Path("/unused"), -1, "b" * 64, "d" * 40, {},
            subject.VALIDATION_ROSTER_SHA256)
        self.derivation = derivation_record(
            self.contract, K=2, selected=0, worker=self.worker.sha256)
        self.row = validation_record(
            self.contract, K=2, attempt=0, worker=self.worker.sha256)
        self.expected_derivation = subject.canonical_json(
            self.derivation).encode("ascii")
        self.expected = subject.canonical_json(self.row).encode("ascii")
        replay_sha256 = subject.sha256_json({
            "schema": subject.SEMANTIC_REPLAY_SCHEMA,
            "K": [2],
        })
        self.constants = mock.patch.multiple(
            subject, K_MAX=2, MAP_BYTES=1, SEMANTIC_REPLAY_K=(2,),
            SEMANTIC_REPLAY_ROSTER_SHA256=replay_sha256)
        self.constants.start()
        self.addCleanup(self.constants.stop)

    def _run(self, returned_derivation, returned_validation,
             trailing_stdout=b""):
        with mock.patch.object(
                subject.subprocess, "Popen",
                side_effect=lambda *args, **kwargs:
                    self.FakeProcess(trailing_stdout)), \
                mock.patch.object(
                    subject, "_read_response_line",
                    side_effect=(returned_derivation, returned_validation)):
            return subject._replay_semantic_sample(
                self.worker, b"\x00", ((2, self.expected_derivation),),
                ((2, self.expected),))

    def test_exact_saved_record_replay_passes(self):
        self.assertEqual(
            self._run(self.expected_derivation, self.expected),
            subject._semantic_replay_roster_sha256())

    def test_replay_rejects_changed_derivation_record(self):
        returned = copy.deepcopy(self.derivation)
        returned["selected_successes"][0]["trace_sha256"] = "e" * 64
        with self.assertRaises(subject.ContractError):
            self._run(
                subject.canonical_json(returned).encode("ascii"),
                self.expected)

    def test_replay_rejects_every_semantic_receipt_class(self):
        mutations = (
            ("source", lambda row: row.__setitem__("source_sha256", "e" * 64)),
            ("base", lambda row: row.__setitem__(
                "base_packet_seed", "0x00000000")),
            ("trace", lambda row: row["cells"][0].__setitem__(
                "trace_sha256", "e" * 64)),
            ("candidates", lambda row: row["cells"][0].__setitem__(
                "attempted_candidates",
                row["cells"][0]["attempted_candidates"] + 1)),
            ("outcome", lambda row: row["cells"][0].__setitem__(
                "outcome", "need_more_at_oh0")),
        )
        for label, mutate in mutations:
            with self.subTest(label=label):
                returned = copy.deepcopy(self.row)
                mutate(returned)
                with self.assertRaises(subject.ContractError):
                    self._run(
                        self.expected_derivation,
                        subject.canonical_json(returned).encode("ascii"))

    def test_replay_rejects_trailing_worker_stdout(self):
        with self.assertRaises(subject.ContractError):
            self._run(
                self.expected_derivation, self.expected,
                trailing_stdout=b"unexpected\n")

    def test_campaign_shard_rejects_trailing_worker_stdout(self):
        process = self.FakeProcess(b"unexpected\n")
        registry = subject._Registry()
        with tempfile.TemporaryDirectory() as temporary, \
                mock.patch.object(
                    subject.subprocess, "Popen", return_value=process), \
                mock.patch.object(
                    subject, "_read_response_line",
                    return_value=self.expected_derivation):
            with self.assertRaises(subject.ContractError):
                subject._run_shard(
                    "derive", 2, 2, Path(temporary) / "shard.jsonl", None,
                    self.contract, self.worker, time.monotonic() + 5, registry)
        self.assertTrue(process.stdin.closed)
        self.assertTrue(process.stdout.closed)

    def test_replay_failure_kills_and_closes_live_worker(self):
        process = self.FakeProcess()
        process.poll = lambda: None
        with mock.patch.object(
                subject.subprocess, "Popen", return_value=process), \
                mock.patch.object(
                    subject, "_read_response_line",
                    side_effect=subject.ContractError("injected")), \
                mock.patch.object(subject, "_kill_process") as kill:
            with self.assertRaises(subject.ContractError):
                subject._replay_semantic_sample(
                    self.worker, b"\x00", ((2, self.expected_derivation),),
                    ((2, self.expected),))
        kill.assert_called_once_with(process)
        self.assertTrue(process.stdin.closed)
        self.assertTrue(process.stdout.closed)
        self.assertTrue(process.stderr.closed)


class RepairBundleChainTests(unittest.TestCase):
    def setUp(self):
        self.contract = subject.load_contract()
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        one_K_replay_sha256 = subject.sha256_json({
            "schema": subject.SEMANTIC_REPLAY_SCHEMA,
            "K": [2],
        })
        self.constants = mock.patch.multiple(
            subject, K_MAX=2, MAP_BYTES=1,
            MAX_AUDIT_BYTES=subject.MAX_RECORD_BYTES,
            SEMANTIC_REPLAY_K=(2,),
            SEMANTIC_REPLAY_ROSTER_SHA256=one_K_replay_sha256)
        self.constants.start()
        self.addCleanup(self.constants.stop)
        self.worker_sha256 = "b" * 64
        self.source_commit = "d" * 40
        self.worker = subject.VerifiedWorker(
            Path("/unused"), -1, self.worker_sha256,
            self.source_commit, {}, subject.VALIDATION_ROSTER_SHA256)
        self.short_binding = {
            field: hashlib.sha256(field.encode("ascii")).hexdigest()
            for field in subject.SHORT_SCREEN_BINDING_FIELDS
        }
        self.short_binding["short_screen_contract_sha256"] = \
            subject.SHORT_SCREEN_CONTRACT_SHA256
        self.manifest = dict(subject.create_manifest(
            "derive", self.contract, self.worker, self.short_binding))
        self.derivation_manifest_anchor = self.manifest["manifest_sha256"]
        self.freeze = dict(subject._freeze_object(self.manifest))
        self.record = derivation_record(
            self.contract, K=2, selected=0, worker=self.worker_sha256)
        self.repair_map = b"\x00"
        self._write_complete_bundle()
        self.validation_root = self.root / "validation"
        self.validation_root.mkdir()
        self.validation_manifest = dict(subject.create_manifest(
            "validate", self.contract, self.worker, self.short_binding,
            self.receipt["map_sha256"], self.receipt["receipt_sha256"]))
        self.validation_manifest_anchor = \
            self.validation_manifest["manifest_sha256"]
        self.validation_freeze = dict(subject._freeze_object(
            self.validation_manifest))
        self.validation_row = validation_record(
            self.contract, K=2, attempt=self.repair_map[0],
            worker=self.worker_sha256)
        self._write_complete_validation_bundle()

    def _write_object(self, name, value):
        (self.root / name).write_bytes(canonical_bytes(value))

    def _audit_bytes(self):
        return canonical_bytes(self.record)

    def _make_summary(self):
        audit = self._audit_bytes()
        histogram = {str(value): self.repair_map.count(value)
                     for value in sorted(set(self.repair_map))}
        summary = {
            "schema": subject.DERIVATION_SUMMARY_SCHEMA,
            "contract_sha256": subject.contract_sha256(self.contract),
            "candidate_profile_sha256":
                subject.candidate_profile_sha256(self.contract),
            "manifest_sha256": self.manifest["manifest_sha256"],
            "audit_sha256": hashlib.sha256(audit).hexdigest(),
            "audit_records": 1,
            "map_sha256": hashlib.sha256(self.repair_map).hexdigest(),
            "map_bytes": 1,
            "maximum_selected_attempt": max(self.repair_map),
            "selected_attempt_histogram": histogram,
            "runtime_search": False,
        }
        reseal(summary, "summary_sha256")
        return summary

    def _make_receipt(self):
        receipt = {
            "schema": subject.MAP_SCHEMA,
            "contract_sha256": subject.contract_sha256(self.contract),
            "candidate_arm": self.contract["candidate_profile"]["arm"],
            "candidate_profile_sha256":
                subject.candidate_profile_sha256(self.contract),
            "construction_seed_basis": subject.PRODUCTION_SEED_BASIS,
            "seed_schedule_sha256":
                subject.PRODUCTION_SEED_SCHEDULE_SHA256,
            "entry_kind": "uint8_attempt_indexed_by_K_minus_2",
            "map_bytes": 1,
            "map_sha256": hashlib.sha256(self.repair_map).hexdigest(),
            "derivation_audit_records": 1,
            "derivation_audit_sha256":
                hashlib.sha256(self._audit_bytes()).hexdigest(),
            "worker_binary_sha256": self.worker_sha256,
            "source_git_commit": self.source_commit,
            "controller_sha256": subject._controller_sha256(),
            "selection_cells_per_selected_attempt":
                subject.SELECTION_CELL_COUNT,
            "selection_rule": self.contract["selection"]["selection_rule"],
            "derivation_manifest_sha256":
                self.manifest["manifest_sha256"],
            "derivation_freeze_sha256": self.freeze["freeze_sha256"],
            "derivation_summary_sha256": self.summary["summary_sha256"],
            **self.short_binding,
        }
        reseal(receipt, "receipt_sha256")
        return receipt

    def _write_complete_bundle(self):
        self.summary = self._make_summary()
        self.receipt = self._make_receipt()
        self._write_object(subject.DERIVATION_FREEZE_NAME, self.freeze)
        (self.root / subject.DERIVATION_AUDIT_NAME).write_bytes(
            self._audit_bytes())
        (self.root / subject.MAP_NAME).write_bytes(self.repair_map)
        self._write_object(subject.MAP_RECEIPT_NAME, self.receipt)
        self._write_object(subject.DERIVATION_SUMMARY_NAME, self.summary)

    def _reseal_audit_summary_receipt(self):
        audit = self._audit_bytes()
        (self.root / subject.DERIVATION_AUDIT_NAME).write_bytes(audit)
        self.summary["audit_sha256"] = hashlib.sha256(audit).hexdigest()
        reseal(self.summary, "summary_sha256")
        self.receipt["derivation_audit_sha256"] = \
            self.summary["audit_sha256"]
        self.receipt["derivation_summary_sha256"] = \
            self.summary["summary_sha256"]
        reseal(self.receipt, "receipt_sha256")
        self._write_object(subject.DERIVATION_SUMMARY_NAME, self.summary)
        self._write_object(subject.MAP_RECEIPT_NAME, self.receipt)

    def _reseal_summary_receipt(self):
        reseal(self.summary, "summary_sha256")
        self.receipt["derivation_summary_sha256"] = \
            self.summary["summary_sha256"]
        reseal(self.receipt, "receipt_sha256")
        self._write_object(subject.DERIVATION_SUMMARY_NAME, self.summary)
        self._write_object(subject.MAP_RECEIPT_NAME, self.receipt)

    def _write_validation_object(self, name, value):
        (self.validation_root / name).write_bytes(canonical_bytes(value))

    def _validation_audit_bytes(self):
        return canonical_bytes(self.validation_row)

    def _make_validation_summary(self):
        summary = {
            "schema": subject.VALIDATION_SUMMARY_SCHEMA,
            "contract_sha256": subject.contract_sha256(self.contract),
            "candidate_profile_sha256":
                subject.candidate_profile_sha256(self.contract),
            "manifest_sha256":
                self.validation_manifest["manifest_sha256"],
            "repair_map_sha256": self.receipt["map_sha256"],
            "map_receipt_sha256": self.receipt["receipt_sha256"],
            "audit_sha256":
                hashlib.sha256(self._validation_audit_bytes()).hexdigest(),
            "audit_records": 1,
            "weak_K_count": 0,
            "weak_K": [],
            "runtime_search": False,
            "disposition": "PASS",
        }
        reseal(summary, "summary_sha256")
        return summary

    def _write_complete_validation_bundle(self):
        self.validation_summary = self._make_validation_summary()
        self._write_validation_object(
            subject.VALIDATION_FREEZE_NAME, self.validation_freeze)
        (self.validation_root / subject.VALIDATION_AUDIT_NAME).write_bytes(
            self._validation_audit_bytes())
        self._write_validation_object(
            subject.VALIDATION_SUMMARY_NAME, self.validation_summary)

    def _reseal_validation_audit_summary(self):
        audit = self._validation_audit_bytes()
        (self.validation_root / subject.VALIDATION_AUDIT_NAME).write_bytes(
            audit)
        self.validation_summary["audit_sha256"] = \
            hashlib.sha256(audit).hexdigest()
        reseal(self.validation_summary, "summary_sha256")
        self._write_validation_object(
            subject.VALIDATION_SUMMARY_NAME, self.validation_summary)

    def _reseal_short_prerequisite_chain(self):
        field = "short_screen_input_sha256"
        self.short_binding[field] = "f" * 64
        self.manifest[field] = self.short_binding[field]
        reseal(self.manifest, "manifest_sha256")
        self.freeze = dict(subject._freeze_object(self.manifest))
        self._write_complete_bundle()
        self.validation_manifest = dict(subject.create_manifest(
            "validate", self.contract, self.worker, self.short_binding,
            self.receipt["map_sha256"], self.receipt["receipt_sha256"]))
        self.validation_freeze = dict(subject._freeze_object(
            self.validation_manifest))
        self._write_complete_validation_bundle()

    def _load_validation(self):
        repair_bundle = subject.load_repair_bundle(
            self.contract, self.root)
        return subject.load_validation_bundle(
            self.contract, self.validation_root, repair_bundle,
            self.worker_sha256, self.source_commit)

    def _export(self, output):
        with mock.patch.object(subject, "_require_clean_source"), \
                mock.patch.object(
                    subject, "_sha256_open_fd",
                    return_value=self.worker_sha256), \
                mock.patch.object(
                    subject, "_replay_semantic_sample",
                    return_value=subject._semantic_replay_roster_sha256()):
            return subject.export_validated_repair_map_include(
                self.contract, self.root, self.validation_root,
                self.worker, self.derivation_manifest_anchor,
                self.validation_manifest_anchor,
                self.receipt["derivation_audit_sha256"],
                self.validation_summary["audit_sha256"], output)

    def test_complete_bundle_loads_and_validation_manifest_binds_receipt(self):
        repair_bundle = subject.load_repair_bundle(
            self.contract, self.root)
        self.assertEqual(repair_bundle.repair_map, b"\x00")
        receipt = repair_bundle.receipt
        worker = subject.VerifiedWorker(
            Path("/unused"), -1, self.worker_sha256,
            self.source_commit, {}, subject.VALIDATION_ROSTER_SHA256)
        manifest = subject.create_manifest(
            "validate", self.contract, worker, self.short_binding,
            receipt["map_sha256"], receipt["receipt_sha256"])
        self.assertEqual(
            manifest["map_receipt_sha256"], receipt["receipt_sha256"])

    def test_fully_resealed_short_chain_needs_external_manifest_anchors(self):
        derivation_audit = self._audit_bytes()
        validation_audit = self._validation_audit_bytes()
        self._reseal_short_prerequisite_chain()
        self.assertEqual(self._audit_bytes(), derivation_audit)
        self.assertEqual(self._validation_audit_bytes(), validation_audit)

        repair_bundle = subject.load_repair_bundle(self.contract, self.root)
        with self.assertRaises(subject.ContractError):
            subject._require_derivation_manifest_anchor(
                repair_bundle, self.derivation_manifest_anchor)
        verified_validation = subject.load_validation_bundle(
            self.contract, self.validation_root, repair_bundle,
            self.worker_sha256, self.source_commit)
        self.assertEqual(verified_validation.summary["disposition"], "PASS")

        output = self.root / "forged-prerequisite-map.inc"
        with mock.patch.object(subject, "_require_clean_source"), \
                mock.patch.object(
                    subject, "_sha256_open_fd",
                    return_value=self.worker_sha256):
            with self.assertRaises(subject.ContractError):
                subject.export_validated_repair_map_include(
                    self.contract, self.root, self.validation_root,
                    self.worker, self.manifest["manifest_sha256"],
                    self.validation_manifest_anchor,
                    hashlib.sha256(derivation_audit).hexdigest(),
                    hashlib.sha256(validation_audit).hexdigest(), output)
        self.assertFalse(output.exists())

    def test_freeze_tamper_rejected_after_local_self_hash_reseal(self):
        self.freeze["manifest"]["K_max"] = 3
        reseal(self.freeze, "freeze_sha256")
        self._write_object(subject.DERIVATION_FREEZE_NAME, self.freeze)
        with self.assertRaises(subject.ContractError):
            subject.load_repair_bundle(self.contract, self.root)

    def test_forged_manifest_rejected_after_full_chain_reseal(self):
        self.manifest["K_max"] = 3
        reseal(self.manifest, "manifest_sha256")
        self.freeze["manifest_sha256"] = self.manifest["manifest_sha256"]
        self.freeze["manifest"] = self.manifest
        reseal(self.freeze, "freeze_sha256")
        self.summary["manifest_sha256"] = self.manifest["manifest_sha256"]
        reseal(self.summary, "summary_sha256")
        self.receipt["derivation_manifest_sha256"] = \
            self.manifest["manifest_sha256"]
        self.receipt["derivation_freeze_sha256"] = \
            self.freeze["freeze_sha256"]
        self.receipt["derivation_summary_sha256"] = \
            self.summary["summary_sha256"]
        reseal(self.receipt, "receipt_sha256")
        self._write_object(subject.DERIVATION_FREEZE_NAME, self.freeze)
        self._write_object(subject.DERIVATION_SUMMARY_NAME, self.summary)
        self._write_object(subject.MAP_RECEIPT_NAME, self.receipt)
        with self.assertRaises(subject.ContractError):
            subject.load_repair_bundle(self.contract, self.root)

    def test_resealed_derivation_manifest_roster_tamper_is_rejected(self):
        self.manifest["validation_roster_sha256"] = "f" * 64
        reseal(self.manifest, "manifest_sha256")
        self.freeze["manifest_sha256"] = self.manifest["manifest_sha256"]
        self.freeze["manifest"] = self.manifest
        reseal(self.freeze, "freeze_sha256")
        self.summary["manifest_sha256"] = self.manifest["manifest_sha256"]
        reseal(self.summary, "summary_sha256")
        self.receipt["derivation_manifest_sha256"] = \
            self.manifest["manifest_sha256"]
        self.receipt["derivation_freeze_sha256"] = \
            self.freeze["freeze_sha256"]
        self.receipt["derivation_summary_sha256"] = \
            self.summary["summary_sha256"]
        reseal(self.receipt, "receipt_sha256")
        self._write_object(subject.DERIVATION_FREEZE_NAME, self.freeze)
        self._write_object(subject.DERIVATION_SUMMARY_NAME, self.summary)
        self._write_object(subject.MAP_RECEIPT_NAME, self.receipt)
        with self.assertRaises(subject.ContractError):
            subject.load_repair_bundle(self.contract, self.root)

    def test_summary_tamper_rejected_after_local_self_hash_reseal(self):
        self.summary["maximum_selected_attempt"] = 1
        reseal(self.summary, "summary_sha256")
        self._write_object(subject.DERIVATION_SUMMARY_NAME, self.summary)
        with self.assertRaises(subject.ContractError):
            subject.load_repair_bundle(self.contract, self.root)

    def test_numeric_metadata_rejects_bool_and_float_aliases(self):
        self.summary["maximum_selected_attempt"] = False
        self._reseal_summary_receipt()
        with self.assertRaises(subject.ContractError):
            subject.load_repair_bundle(self.contract, self.root)

        self._write_complete_bundle()
        self.summary["selected_attempt_histogram"] = {"0": 1.0}
        self._reseal_summary_receipt()
        with self.assertRaises(subject.ContractError):
            subject.load_repair_bundle(self.contract, self.root)

        self._write_complete_bundle()
        self.receipt["map_bytes"] = True
        reseal(self.receipt, "receipt_sha256")
        self._write_object(subject.MAP_RECEIPT_NAME, self.receipt)
        with self.assertRaises(subject.ContractError):
            subject.load_repair_bundle(self.contract, self.root)

    def test_map_audit_disagreement_rejected_after_full_chain_reseal(self):
        self.repair_map = b"\x01"
        (self.root / subject.MAP_NAME).write_bytes(self.repair_map)
        self.summary["map_sha256"] = hashlib.sha256(
            self.repair_map).hexdigest()
        self.summary["maximum_selected_attempt"] = 1
        self.summary["selected_attempt_histogram"] = {"1": 1}
        reseal(self.summary, "summary_sha256")
        self.receipt["map_sha256"] = self.summary["map_sha256"]
        self.receipt["derivation_summary_sha256"] = \
            self.summary["summary_sha256"]
        reseal(self.receipt, "receipt_sha256")
        self._write_object(subject.DERIVATION_SUMMARY_NAME, self.summary)
        self._write_object(subject.MAP_RECEIPT_NAME, self.receipt)
        with self.assertRaises(subject.ContractError):
            subject.load_repair_bundle(self.contract, self.root)

    def test_audit_forgery_rejected_after_full_chain_reseal(self):
        self.record["ordinal"] = 1
        self._reseal_audit_summary_receipt()
        with self.assertRaises(subject.ContractError):
            subject.load_repair_bundle(self.contract, self.root)

    def test_symlinked_audit_is_rejected(self):
        audit = self.root / subject.DERIVATION_AUDIT_NAME
        target = self.root / "saved-audit.jsonl"
        os.replace(audit, target)
        audit.symlink_to(target.name)
        with self.assertRaises(subject.ContractError):
            subject.load_repair_bundle(self.contract, self.root)

    def test_receipt_self_hash_tamper_is_rejected(self):
        self.receipt["receipt_sha256"] = "0" * 64
        self._write_object(subject.MAP_RECEIPT_NAME, self.receipt)
        with self.assertRaises(subject.ContractError):
            subject.load_repair_bundle(self.contract, self.root)

    def test_complete_validation_bundle_authenticates_zero_weak_K_pass(self):
        verified = self._load_validation()
        self.assertEqual(
            verified.freeze["freeze_sha256"],
            self.validation_freeze["freeze_sha256"])
        self.assertEqual(verified.summary["disposition"], "PASS")
        self.assertEqual(verified.summary["weak_K"], [])

    def test_validation_manifest_tamper_rejected_after_self_hash_reseal(self):
        self.validation_manifest["repair_map_sha256"] = "e" * 64
        reseal(self.validation_manifest, "manifest_sha256")
        self.validation_freeze["manifest_sha256"] = \
            self.validation_manifest["manifest_sha256"]
        self.validation_freeze["manifest"] = self.validation_manifest
        reseal(self.validation_freeze, "freeze_sha256")
        self._write_validation_object(
            subject.VALIDATION_FREEZE_NAME, self.validation_freeze)
        with self.assertRaises(subject.ContractError):
            self._load_validation()

    def test_resealed_validation_manifest_roster_tamper_is_rejected(self):
        self.validation_manifest["validation_roster_sha256"] = "f" * 64
        reseal(self.validation_manifest, "manifest_sha256")
        self.validation_freeze["manifest_sha256"] = \
            self.validation_manifest["manifest_sha256"]
        self.validation_freeze["manifest"] = self.validation_manifest
        reseal(self.validation_freeze, "freeze_sha256")
        self.validation_summary["manifest_sha256"] = \
            self.validation_manifest["manifest_sha256"]
        reseal(self.validation_summary, "summary_sha256")
        self._write_validation_object(
            subject.VALIDATION_FREEZE_NAME, self.validation_freeze)
        self._write_validation_object(
            subject.VALIDATION_SUMMARY_NAME, self.validation_summary)
        with self.assertRaises(subject.ContractError):
            self._load_validation()

    def test_validation_summary_gate_tamper_rejected_after_reseal(self):
        mutations = (
            {"runtime_search": True},
            {"disposition": "REJECT"},
            {"weak_K_count": 1, "weak_K": [2]},
            {"repair_map_sha256": "e" * 64},
            {"map_receipt_sha256": "e" * 64},
        )
        original = copy.deepcopy(self.validation_summary)
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                self.validation_summary = copy.deepcopy(original)
                self.validation_summary.update(mutation)
                reseal(self.validation_summary, "summary_sha256")
                self._write_validation_object(
                    subject.VALIDATION_SUMMARY_NAME,
                    self.validation_summary)
                with self.assertRaises(subject.ContractError):
                    self._load_validation()

    def test_validation_audit_tamper_rejected_after_summary_reseal(self):
        self.validation_row["ordinal"] = 1
        audit = self._validation_audit_bytes()
        (self.validation_root / subject.VALIDATION_AUDIT_NAME).write_bytes(
            audit)
        self.validation_summary["audit_sha256"] = \
            hashlib.sha256(audit).hexdigest()
        reseal(self.validation_summary, "summary_sha256")
        self._write_validation_object(
            subject.VALIDATION_SUMMARY_NAME, self.validation_summary)
        with self.assertRaises(subject.ContractError):
            self._load_validation()

    def test_validation_source_or_base_drift_rejected_after_full_reseal(self):
        original = copy.deepcopy(self.validation_row)
        mutations = (
            ("source_sha256", "e" * 64),
            ("base_packet_seed", "0x00000000"),
            ("base_precode_seed", "0x0000000000000000"),
        )
        for field, replacement in mutations:
            with self.subTest(field=field):
                self.validation_row = copy.deepcopy(original)
                self.validation_row[field] = replacement
                attempt = self.validation_row["construction_attempt"]
                if field == "base_packet_seed":
                    self.validation_row["effective_packet_seed"] = \
                        subject._effective_packet_seed(0, attempt)
                elif field == "base_precode_seed":
                    self.validation_row["effective_precode_seed"] = \
                        subject._effective_precode_seed(0, attempt)
                self._reseal_validation_audit_summary()
                with self.assertRaises(subject.ContractError):
                    self._load_validation()

    def test_validation_worker_or_source_substitution_is_rejected(self):
        repair_bundle = subject.load_repair_bundle(
            self.contract, self.root)
        for worker, source in (
                ("e" * 64, self.source_commit),
                (self.worker_sha256, "e" * 40)):
            with self.subTest(worker=worker, source=source), \
                    self.assertRaises(subject.ContractError):
                subject.load_validation_bundle(
                    self.contract, self.validation_root, repair_bundle,
                    worker, source)

    def test_symlinked_validation_artifacts_are_rejected(self):
        for name in (
                subject.VALIDATION_FREEZE_NAME,
                subject.VALIDATION_AUDIT_NAME,
                subject.VALIDATION_SUMMARY_NAME):
            with self.subTest(name=name):
                path = self.validation_root / name
                target = self.validation_root / ("saved-" + name)
                os.replace(path, target)
                path.symlink_to(target.name)
                with self.assertRaises(subject.ContractError):
                    self._load_validation()
                path.unlink()
                os.replace(target, path)

    def test_validated_map_include_is_deterministic_and_no_clobber(self):
        first = self.root / "repair-map-a.inc"
        second = self.root / "repair-map-b.inc"
        first_result = self._export(first)
        second_result = self._export(second)
        self.assertEqual(first.read_bytes(), second.read_bytes())
        self.assertEqual(
            first_result["output_sha256"], second_result["output_sha256"])
        text = first.read_text(encoding="ascii")
        self.assertIn(
            "// repair_map_sha256={}".format(self.receipt["map_sha256"]),
            text)
        self.assertIn(
            "// validation_summary_sha256={}".format(
                self.validation_summary["summary_sha256"]), text)
        self.assertIn(
            "// shared_coordinate_identity_sha256={}".format(
                self._load_validation().shared_identity_sha256), text)
        self.assertIn(
            "// semantic_derivation_validation_replay_roster_sha256={}".format(
                subject._semantic_replay_roster_sha256()), text)
        self.assertIn("// semantic_derivation_replay_records=1", text)
        self.assertIn("// semantic_validation_replay_records=1", text)
        self.assertIn("    0x00,\n", text)
        with self.assertRaises(subject.ContractError):
            self._export(first)

    def test_export_never_publishes_on_dirty_source_or_replay_failure(self):
        output = self.root / "blocked-map.inc"
        with mock.patch.object(
                subject, "_require_clean_source",
                side_effect=subject.ContractError("dirty")):
            with self.assertRaises(subject.ContractError):
                subject.export_validated_repair_map_include(
                    self.contract, self.root, self.validation_root,
                    self.worker, self.derivation_manifest_anchor,
                    self.validation_manifest_anchor,
                    self.receipt["derivation_audit_sha256"],
                    self.validation_summary["audit_sha256"], output)
        self.assertFalse(output.exists())

    def test_export_second_provenance_gate_is_fail_closed(self):
        output = self.root / "late-blocked-map.inc"
        derivation_audit = self.receipt["derivation_audit_sha256"]
        validation_audit = self.validation_summary["audit_sha256"]
        replay_sha256 = subject._semantic_replay_roster_sha256()

        with mock.patch.object(
                subject, "_require_clean_source",
                side_effect=(None, subject.ContractError("late dirty"))) \
                as clean, \
                mock.patch.object(
                    subject, "_sha256_open_fd",
                    return_value=self.worker_sha256), \
                mock.patch.object(
                    subject, "_replay_semantic_sample",
                    return_value=replay_sha256):
            with self.assertRaises(subject.ContractError):
                subject.export_validated_repair_map_include(
                    self.contract, self.root, self.validation_root,
                    self.worker, self.derivation_manifest_anchor,
                    self.validation_manifest_anchor,
                    derivation_audit, validation_audit, output)
        self.assertEqual(clean.call_count, 2)
        self.assertFalse(output.exists())

        with mock.patch.object(subject, "_require_clean_source"), \
                mock.patch.object(
                    subject, "_sha256_open_fd",
                    side_effect=(self.worker_sha256, "e" * 64)), \
                mock.patch.object(
                    subject, "_replay_semantic_sample",
                    return_value=replay_sha256):
            with self.assertRaises(subject.ContractError):
                subject.export_validated_repair_map_include(
                    self.contract, self.root, self.validation_root,
                    self.worker, self.derivation_manifest_anchor,
                    self.validation_manifest_anchor,
                    derivation_audit, validation_audit, output)
        self.assertFalse(output.exists())

        with mock.patch.object(subject, "_require_clean_source"), \
                mock.patch.object(
                    subject, "_sha256_open_fd",
                    return_value=self.worker_sha256):
            with self.assertRaises(subject.ContractError):
                subject.export_validated_repair_map_include(
                    self.contract, self.root, self.validation_root,
                    self.worker, self.derivation_manifest_anchor,
                    self.validation_manifest_anchor, "e" * 64,
                    self.validation_summary["audit_sha256"], output)
        self.assertFalse(output.exists())

        with mock.patch.object(subject, "_require_clean_source"), \
                mock.patch.object(
                    subject, "_sha256_open_fd",
                    return_value=self.worker_sha256), \
                mock.patch.object(
                    subject, "_replay_semantic_sample",
                    side_effect=subject.ContractError("replay drift")):
            with self.assertRaises(subject.ContractError):
                subject.export_validated_repair_map_include(
                    self.contract, self.root, self.validation_root,
                    self.worker, self.derivation_manifest_anchor,
                    self.validation_manifest_anchor,
                    self.receipt["derivation_audit_sha256"],
                    self.validation_summary["audit_sha256"], output)
        self.assertFalse(output.exists())

    def test_export_map_include_command_is_registered(self):
        args = subject._parser().parse_args([
            "export-map-include", "--worker", "/unused",
            "--worker-sha256", self.worker_sha256,
            "--source-commit", self.source_commit,
            "--repair-dir", str(self.root),
            "--validation-dir", str(self.validation_root),
            "--derivation-manifest-sha256",
            self.derivation_manifest_anchor,
            "--validation-manifest-sha256",
            self.validation_manifest_anchor,
            "--derivation-audit-sha256",
            self.receipt["derivation_audit_sha256"],
            "--validation-audit-sha256",
            self.validation_summary["audit_sha256"],
            "--output", str(self.root / "repair-map.inc"),
        ])
        self.assertEqual(args.command, "export-map-include")


class ExitDispositionTests(unittest.TestCase):
    def test_validation_pass_and_reject_have_distinct_exit_codes(self):
        self.assertEqual(
            subject._result_exit_code("validate", {"disposition": "PASS"}), 0)
        self.assertEqual(
            subject._result_exit_code("validate", {"disposition": "REJECT"}),
            subject.SCIENTIFIC_REJECT_EXIT_CODE)
        self.assertEqual(subject.SCIENTIFIC_REJECT_EXIT_CODE, 2)

    def test_invalid_validation_disposition_fails_closed(self):
        for result in ({}, {"disposition": "INVALID"}):
            with self.subTest(result=result), \
                    self.assertRaises(subject.ContractError):
                subject._result_exit_code("validate", result)
        self.assertEqual(subject._result_exit_code("derive", {}), 0)


if __name__ == "__main__":
    unittest.main()
