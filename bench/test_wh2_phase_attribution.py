#!/usr/bin/env python3
"""Synthetic fail-closed tests for the WH2 n16 phase consumer."""

from __future__ import annotations

import copy
import hashlib
import math
from pathlib import Path
import sys
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_benchmark_contract as contract_api
import wh2_native_short_screen as native_api
import wh2_phase_attribution as subject


def digest(label: str) -> str:
    return hashlib.sha256(label.encode("ascii")).hexdigest()


def canonical_jsonl(values) -> bytes:
    return "".join(
        contract_api.canonical_json(value) + "\n" for value in values
    ).encode("ascii")


def fake_message_sha256(base) -> str:
    return digest("message:{}:{}".format(base["K"], base["block_bytes"]))


def counters(tag: int):
    result = {}
    for index, field in enumerate(sorted(subject.COUNTER_FIELDS)):
        result[field] = tag * 100 + index
    return result


def success_timing(arm: str, block: int, zero_phase=None):
    if arm == subject.LEFT_ARM:
        scale = 2 if block == 0 else 0.5
    else:
        scale = 1
    values = {
        "outer_ns": int(1000 * scale),
        "build_ns": int(100 * scale),
        "peel_ns": int(100 * scale),
        "project_ns": int(100 * scale),
        "residual_ns": int(100 * scale),
        "back_sub_ns": int(100 * scale),
    }
    if zero_phase is not None:
        values[zero_phase + "_ns"] = 0
    return values


def observation(arm: str, arm_counters, block: int, comparable: bool,
                zero_phase=None):
    if comparable:
        result = {
            "arm": arm,
            "bytes_verified": True,
            "counter_sha256": contract_api.sha256_json(arm_counters),
            "outcome": "success",
            "wirehair_result": 0,
        }
        result.update(success_timing(arm, block, zero_phase))
        return result
    wirehair_result = 1 if arm == subject.LEFT_ARM else 0
    return {
        "arm": arm,
        "back_sub_ns": None,
        "build_ns": None,
        "bytes_verified": wirehair_result == 0,
        "counter_sha256": contract_api.sha256_json(arm_counters),
        "outcome": "need_more_at_cap" if wirehair_result == 1 else "success",
        "outer_ns": None,
        "peel_ns": None,
        "project_ns": None,
        "residual_ns": None,
        "wirehair_result": wirehair_result,
    }


class PhaseAttributionTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.contract = contract_api.load_contract()
        cls.description = {
            "arms": [{
                "arm": arm,
                "arm_descriptor_sha256": descriptor,
                "codec": codec,
            } for arm, codec, descriptor in subject.EXPECTED_ARMS],
            "binary_sha256": digest("phase-worker"),
            "schema": subject.runner_api.DESCRIPTION_SCHEMA,
            "source_git_commit": "1" * 40,
        }
        cls.coordinates = subject._phase_coordinates(cls.contract)
        cls.traces = []
        for ordinal, coordinate in enumerate(cls.coordinates):
            trace_sha256 = digest("trace:{}".format(ordinal))
            packet_count = coordinate["qualified"]["K"] + \
                subject.RECEIVE_OVERHEAD_CAP
            cls.traces.append({
                "candidate_count": packet_count,
                "coordinate_ordinal": ordinal,
                "packet_count": packet_count,
                "phase_cell_sha256_by_profile": [
                    contract_api.sha256_json(subject._phase_cell(
                        coordinate, profile, cls.description, trace_sha256))
                    for profile in range(2)
                ],
                "schema": subject.PHASE_TRACE_SCHEMA,
                "source_base_cell_sha256":
                    coordinate["qualified"]["base_cell_sha256"],
                "trace_qualified_timing_cell_sha256":
                    contract_api.sha256_json(coordinate["qualified"]),
                "trace_sha256": trace_sha256,
            })
        cls.trace_data = canonical_jsonl(cls.traces)
        cls.coordinate_cpus = tuple(range(subject.WORKER_COUNT)) * 3
        cls.runtime_workers = {
            cpu: (1000 + cpu, 1 + cpu)
            for cpu in range(subject.WORKER_COUNT)
        }
        cls.window_start_ns = 1_000_000_000_000
        cls.window_end_ns = cls.window_start_ns + 1_000_000

    @classmethod
    def build_payload(cls, ordinal: int, weak: bool = False,
                      zero_phase=None):
        coordinate = cls.coordinates[ordinal]
        qualified = coordinate["qualified"]
        trace = cls.traces[ordinal]
        cell = subject._phase_cell(
            coordinate, subject.PROFILE_ORDINAL, cls.description,
            trace["trace_sha256"])
        left_counters = counters(1 + ordinal * 2)
        right_counters = counters(2 + ordinal * 2)
        by_arm = {
            subject.LEFT_ARM: left_counters,
            subject.RIGHT_ARM: right_counters,
        }
        comparable = not weak
        measured = []
        observations_by_slot = {slot: [] for slot in range(8)}
        for block in range(2):
            for repeat in range(subject.PAIR_COUNT):
                for block_slot in range(4):
                    slot = block * 4 + block_slot
                    arm = subject._expected_observation_arm(block, slot)
                    value = observation(
                        arm, by_arm[arm], block, comparable, zero_phase)
                    measured.append({
                        "block": block,
                        "observation": value,
                        "repeat": repeat,
                        "slot": slot,
                    })
                    observations_by_slot[slot].append(value)
        slots = []
        for slot in range(8):
            row = {"slot": slot}
            for key in subject.TIMING_KEYS:
                row[key] = sum(
                    item[key] for item in observations_by_slot[slot]
                ) if comparable else None
            slots.append(row)
        return {
            "K": qualified["K"],
            "band": qualified["band"],
            "base_loss_seed": qualified["base_loss_seed"],
            "binary_sha256": cls.description["binary_sha256"],
            "block_bytes": qualified["block_bytes"],
            "block_muladds_semantics":
                "full-block-gf256-multiply-add-divide-normalize-operations",
            "cell_sha256": contract_api.sha256_json(cell),
            "construction_attempt": 0,
            "coordinate_ordinal": ordinal,
            "fixed_received_overhead": subject.FIXED_RECEIVED_OVERHEAD,
            "interleave_policy":
                "self-counterbalanced-repeat-major-v1",
            "invocations_per_slot": subject.INVOCATIONS_PER_SLOT,
            "left_arm": cell["left_arm"],
            "left_arm_descriptor_sha256":
                cell["left_arm_descriptor_sha256"],
            "left_non_timing_counters": left_counters,
            "left_realized_construction_sha256":
                cell["left_realized_construction_sha256"],
            "left_repair_map_sha256": subject.ZERO_SHA256,
            "loss_ppm": qualified["loss_ppm"],
            "loss_retry_offset": 0,
            "loss_seed": qualified["loss_seed"],
            "measured_observations": measured,
            "order": "ABBA",
            "panel_comparable": comparable,
            "panel_kind": "ab",
            "profile": subject.PROFILE_NAME,
            "profile_ordinal": subject.PROFILE_ORDINAL,
            "receive_overhead_cap": subject.RECEIVE_OVERHEAD_CAP,
            "replicate": 0,
            "right_arm": cell["right_arm"],
            "right_arm_descriptor_sha256":
                cell["right_arm_descriptor_sha256"],
            "right_non_timing_counters": right_counters,
            "right_realized_construction_sha256":
                cell["right_realized_construction_sha256"],
            "right_repair_map_sha256": subject.ZERO_SHA256,
            "schedule": "iid",
            "scope": "decoder_solve",
            "slot_sums": slots,
            "source_base_cell_sha256": qualified["base_cell_sha256"],
            "trace_qualified_timing_cell_sha256":
                contract_api.sha256_json(qualified),
            "trace_sha256": trace["trace_sha256"],
            "warmups": {
                "left": observation(
                    subject.LEFT_ARM, left_counters, 0, comparable,
                    zero_phase),
                "right": observation(
                    subject.RIGHT_ARM, right_counters, 0, comparable,
                    zero_phase),
            },
            "weak_ledger": weak,
        }

    @classmethod
    def build_records(cls, weak=frozenset(), zero_phase=None):
        zero_phase = {} if zero_phase is None else zero_phase
        rows = []
        for ordinal, coordinate in enumerate(cls.coordinates):
            payload = cls.build_payload(
                ordinal, ordinal in weak, zero_phase.get(ordinal))
            worker_ordinal = ordinal * 2
            cpu = cls.coordinate_cpus[ordinal]
            started = cls.window_start_ns + 1000 + ordinal * 100
            rows.append({
                "cpu": cpu,
                "finished_monotonic_ns": started + 50,
                "message_sha256": fake_message_sha256(coordinate["base"]),
                "ordinal": worker_ordinal,
                "payload": payload,
                "schema": subject.PHASE_RECORD_SCHEMA,
                "started_monotonic_ns": started,
                "work_sha256": subject._expected_work_sha256(
                    payload["cell_sha256"], worker_ordinal),
                "worker_binary_sha256": cls.description["binary_sha256"],
                "worker_pid": cls.runtime_workers[cpu][0],
                "worker_process_start_ticks":
                    cls.runtime_workers[cpu][1],
            })
        return list(reversed(rows))

    def validate(self, records, **overrides):
        arguments = {
            "contract": self.contract,
            "description": self.description,
            "trace_data": self.trace_data,
            "data": canonical_jsonl(records),
            "coordinate_cpus": self.coordinate_cpus,
            "runtime_workers": self.runtime_workers,
            "window_start_ns": self.window_start_ns,
            "window_end_ns": self.window_end_ns,
        }
        arguments.update(overrides)
        with mock.patch.object(
                subject, "_deterministic_message_sha256",
                side_effect=fake_message_sha256):
            return subject.validate_native_records(**arguments)

    @staticmethod
    def retime_left_outer(payload, numerator: int, denominator: int):
        warm = payload["warmups"]["left"]
        warm["outer_ns"] = warm["outer_ns"] * numerator // denominator
        for row in payload["measured_observations"]:
            value = row["observation"]
            if value["arm"] == subject.LEFT_ARM:
                value["outer_ns"] = \
                    value["outer_ns"] * numerator // denominator
        for slot in payload["slot_sums"]:
            slot["outer_ns"] = sum(
                row["observation"]["outer_ns"]
                for row in payload["measured_observations"]
                if row["slot"] == slot["slot"])

    @staticmethod
    def aggregate(analysis, scope, key):
        matches = [
            row for row in analysis["aggregates"]
            if row["scope"] == scope and row["group_key"] == key
        ]
        if len(matches) != 1:
            raise AssertionError("aggregate lookup is not unique")
        return matches[0]

    def test_exact_trace_manifest_binds_both_profiles(self):
        rows, manifest_sha256 = subject.validate_trace_manifest(
            self.contract, self.description, self.trace_data)
        self.assertEqual(len(rows), 24)
        self.assertEqual(len({
            value for row in rows
            for value in row["phase_cell_sha256_by_profile"]
        }), 48)
        self.assertEqual(
            manifest_sha256,
            hashlib.sha256(
                subject.TRACE_MANIFEST_DOMAIN + self.trace_data).hexdigest())

        changed = copy.deepcopy(self.traces)
        changed[0]["phase_cell_sha256_by_profile"][1] = digest("forged")
        with self.assertRaises(subject.PhaseRunnerError):
            subject.validate_trace_manifest(
                self.contract, self.description, canonical_jsonl(changed))

        malformed_description = copy.deepcopy(self.description)
        malformed_description["source_git_commit"] = 1
        with self.assertRaises(subject.PhaseRunnerError):
            subject.validate_trace_manifest(
                self.contract, malformed_description, self.trace_data)

    def test_physical_record_cap_precedes_semantic_parsing(self):
        limit = 4096
        fragmented = b"x" * (limit - 1) + b"\r" + \
            b"y" * (limit - 1) + b"\n"
        parser = mock.Mock(
            side_effect=AssertionError("oversized record reached parser"))
        with mock.patch.object(
                subject.runner_api, "_parse_canonical_line", parser), \
                self.assertRaisesRegex(
                    subject.PhaseRunnerError, "line 1 exceeds its bound"):
            subject._parse_jsonl_bytes(
                fragmented, "fragmented fixture", 1, len(fragmented), limit)
        parser.assert_not_called()

    def test_n16_records_and_counterbalanced_statistics(self):
        payloads, metadata = self.validate(self.build_records())
        self.assertEqual(
            [row["coordinate_ordinal"] for row in payloads], list(range(24)))
        self.assertEqual(metadata["record_count"], 24)
        self.assertEqual(metadata["worker_cpus"], list(range(8)))
        self.assertEqual(
            metadata["coordinate_cpus"], list(self.coordinate_cpus))
        analysis = subject.build_phase_analysis(payloads)
        global_outer = self.aggregate(analysis, "global", "all")[
            "phases"]["outer"]
        self.assertEqual(global_outer["independent_pair_count"], 8)
        self.assertEqual(global_outer["status"], "comparable")
        self.assertAlmostEqual(global_outer["geometric_mean_ratio"], 1.0)
        self.assertAlmostEqual(global_outer["lower_95_ratio"], 1.0)
        self.assertAlmostEqual(global_outer["upper_95_ratio"], 1.0)
        self.assertEqual(analysis["student_t_critical_df7"],
                         subject.T95_DF7)
        self.assertEqual(analysis["decision"]["outcome"], "n16_sufficient")
        self.assertIs(
            analysis["decision"][
                "global_and_width_upper95_within_threshold"], True)
        self.assertEqual(len([
            row for row in analysis["aggregates"] if row["scope"] == "K"
        ]), 12)
        self.assertEqual(
            self.aggregate(analysis, "width", 64)["phases"]["outer"][
                "cell_count"], 12)
        first = analysis["cells"][0]
        self.assertEqual(
            first["arm_outcome_ledger"]["left"]["outcome"], "success")
        self.assertEqual(
            first["left_minus_right_non_timing_counter_deltas"][
                "packet_rows"], -100)
        self.assertEqual(len(first["phases"]["outer"]["pair_log_ratios"]), 8)
        self.assertTrue(all(math.isclose(value, 0.0, abs_tol=1e-15)
                            for value in
                            first["phases"]["outer"]["pair_log_ratios"]))

    def test_weak_cell_nulls_affected_aggregates_and_keeps_ledger(self):
        payloads, _ = self.validate(self.build_records(weak={0}))
        analysis = subject.build_phase_analysis(payloads)
        first = analysis["cells"][0]
        self.assertTrue(first["weak_ledger"])
        self.assertEqual(
            first["arm_outcome_ledger"]["left"]["wirehair_result"], 1)
        self.assertEqual(
            first["arm_outcome_ledger"]["right"]["wirehair_result"], 0)
        self.assertEqual(first["phases"]["outer"]["status"],
                         "non_comparable_weak")
        for scope, key in (
                ("global", "all"), ("width", 64), ("K", 8),
                ("band", "2-100")):
            phase = self.aggregate(analysis, scope, key)["phases"]["outer"]
            self.assertEqual(phase["status"], "non_comparable_weak")
            self.assertIsNone(phase["geometric_mean_ratio"])
        self.assertEqual(
            self.aggregate(analysis, "width", 1280)["phases"]["outer"][
                "status"], "comparable")
        self.assertEqual(analysis["decision"]["outcome"], "inconclusive")

    def test_zero_component_is_unresolved_without_poisoning_outer(self):
        payloads, _ = self.validate(self.build_records(
            zero_phase={0: "project"}))
        analysis = subject.build_phase_analysis(payloads)
        self.assertEqual(
            analysis["cells"][0]["phases"]["project"]["status"],
            "unresolved_zero_phase")
        self.assertEqual(
            self.aggregate(analysis, "global", "all")["phases"][
                "project"]["status"], "unresolved_zero_phase")
        self.assertEqual(
            self.aggregate(analysis, "global", "all")["phases"][
                "outer"]["status"], "comparable")
        self.assertEqual(analysis["decision"]["outcome"], "n16_sufficient")

    def test_visible_success_rejects_zero_internal_phase_total(self):
        records = self.build_records()
        observation = records[-1]["payload"]["measured_observations"][0][
            "observation"]
        for _, key in subject.PHASE_KEYS[1:]:
            observation[key] = 0
        with self.assertRaisesRegex(
                subject.PhaseRunnerError, "invalid visible phase timing"):
            self.validate(records)

    def test_visible_slot_rejects_zero_internal_phase_total(self):
        records = self.build_records()
        slot = records[-1]["payload"]["slot_sums"][0]
        for _, key in subject.PHASE_KEYS[1:]:
            slot[key] = 0
        with self.assertRaisesRegex(
                subject.PhaseRunnerError,
                "phase slot has invalid visible timing totals"):
            self.validate(records)

    def test_frozen_two_percent_decision_is_conservative(self):
        records = self.build_records()
        for row in records:
            self.retime_left_outer(row["payload"], 103, 100)
        payloads, _ = self.validate(records)
        analysis = subject.build_phase_analysis(payloads)
        decision = analysis["decision"]
        self.assertEqual(decision["outcome"], "inconclusive")
        self.assertIs(
            decision["global_and_width_upper95_within_threshold"], False)
        self.assertIs(
            decision["no_K_or_band_lower95_above_threshold"], False)
        self.assertAlmostEqual(
            self.aggregate(analysis, "global", "all")["phases"]["outer"][
                "geometric_mean_ratio"], 1.03)

    def test_counter_native_widths_are_exact(self):
        value = {field: 0 for field in subject.COUNTER_FIELDS}
        for field in subject.COUNTER_UINT32_FIELDS:
            value[field] = subject.UINT32_MAX
        for field in subject.COUNTER_UINT64_FIELDS:
            value[field] = subject.UINT64_MAX
        subject._validate_counters(value, "boundary counters")
        for field, invalid in (
                (next(iter(subject.COUNTER_UINT32_FIELDS)), 1 << 32),
                (next(iter(subject.COUNTER_UINT64_FIELDS)), 1 << 64)):
            changed = dict(value)
            changed[field] = invalid
            with self.assertRaises(subject.PhaseRunnerError):
                subject._validate_counters(changed, "overflow counters")
        changed = dict(value)
        changed["packet_rows"] = True
        with self.assertRaises(subject.PhaseRunnerError):
            subject._validate_counters(changed, "boolean counter")

    def test_corrupt_native_provenance_and_chronology_fail_closed(self):
        mutators = {}

        def mutate(name, callback):
            records = self.build_records()
            by_coordinate = {
                row["payload"]["coordinate_ordinal"]: row for row in records
            }
            callback(by_coordinate)
            mutators[name] = records

        mutate("n24-ordinal", lambda rows: rows[0].__setitem__("ordinal", 1))
        mutate("cpu-schedule", lambda rows: rows[8].__setitem__("cpu", 1))
        mutate("pid-drift", lambda rows: rows[8].__setitem__(
            "worker_pid", rows[8]["worker_pid"] + 100))
        mutate("message", lambda rows: rows[0].__setitem__(
            "message_sha256", digest("substitute-source")))
        mutate("work", lambda rows: rows[0].__setitem__(
            "work_sha256", digest("substitute-work")))
        mutate("worker-overlap", lambda rows: (
            rows[8].__setitem__(
                "started_monotonic_ns",
                rows[0]["finished_monotonic_ns"] - 1),
            rows[8].__setitem__(
                "finished_monotonic_ns",
                rows[0]["finished_monotonic_ns"] + 1)))
        mutate("worker-reverse-order", lambda rows: (
            rows[8].__setitem__(
                "started_monotonic_ns",
                rows[0]["started_monotonic_ns"] - 2),
            rows[8].__setitem__(
                "finished_monotonic_ns",
                rows[0]["started_monotonic_ns"] - 1)))
        mutate("cross-lane-wave-overlap", lambda rows: (
            rows[8].__setitem__(
                "started_monotonic_ns",
                rows[0]["finished_monotonic_ns"] + 1),
            rows[8].__setitem__(
                "finished_monotonic_ns",
                rows[0]["finished_monotonic_ns"] + 2)))
        mutate("chronology", lambda rows: rows[0]["payload"][
            "measured_observations"][0].__setitem__("slot", 1))
        mutate("slot-sum", lambda rows: rows[0]["payload"][
            "slot_sums"][0].__setitem__("outer_ns", 1))
        mutate("bool-result", lambda rows: rows[0]["payload"]["warmups"][
            "left"].__setitem__("wirehair_result", False))
        mutate("int63-timing", lambda rows: rows[0]["payload"][
            "measured_observations"][0]["observation"].__setitem__(
                "outer_ns", subject.INT63_MAX + 1))
        for name, records in mutators.items():
            with self.subTest(name=name):
                with self.assertRaises(subject.PhaseRunnerError):
                    self.validate(records)

        weak = self.build_records(weak={0})
        by_coordinate = {
            row["payload"]["coordinate_ordinal"]: row for row in weak
        }
        by_coordinate[0]["payload"]["measured_observations"][0][
            "observation"]["outer_ns"] = 1
        with self.assertRaises(subject.PhaseRunnerError):
            self.validate(weak)

    def test_execution_geometry_and_canonical_boundaries_fail_closed(self):
        records = self.build_records()
        with self.assertRaises(subject.PhaseRunnerError):
            self.validate(records, coordinate_cpus=None)
        with self.assertRaises(subject.PhaseRunnerError):
            self.validate(records, coordinate_cpus=tuple(range(8)) * 2 +
                          tuple(reversed(range(8))))
        runtime = dict(self.runtime_workers)
        del runtime[7]
        with self.assertRaises(subject.PhaseRunnerError):
            self.validate(records, runtime_workers=runtime)
        with self.assertRaises(subject.PhaseRunnerError):
            self.validate(records, window_start_ns=0)
        with self.assertRaises(subject.PhaseRunnerError):
            self.validate(
                records, window_end_ns=self.window_start_ns)
        zero_duration = self.build_records()
        zero_duration[0]["finished_monotonic_ns"] = \
            zero_duration[0]["started_monotonic_ns"]
        with self.assertRaises(subject.PhaseRunnerError):
            self.validate(zero_duration)
        noncanonical = canonical_jsonl(records).replace(b"{", b"{ ", 1)
        with self.assertRaises(subject.PhaseRunnerError):
            self.validate(records, data=noncanonical)

    def test_canonical_selector_categorically_rejects_phase_evidence(self):
        self.assertNotEqual(subject.PHASE_RECORD_SCHEMA,
                            native_api.TIMING_RECORD_SCHEMA)
        self.assertNotEqual(subject.PHASE_TRACE_SCHEMA,
                            native_api.TRACE_RECORD_SCHEMA)
        self.assertNotEqual(
            subject.PHASE_ANALYSIS_SCHEMA,
            contract_api.SCHEMA + ".timing-summary.v1")
        self.assertNotEqual(
            subject.PHASE_DECISION_SCHEMA,
            contract_api.SCHEMA + ".architecture-selection.v1")
        with self.assertRaises(native_api.NativeEvidenceError):
            native_api._trace_cells(
                self.contract, "phase_attribution", "development")
        with self.assertRaises(contract_api.ContractError):
            contract_api.validate_selection_receipt(self.contract, {
                "schema": subject.PHASE_ANALYSIS_SCHEMA,
            })


if __name__ == "__main__":
    unittest.main()
