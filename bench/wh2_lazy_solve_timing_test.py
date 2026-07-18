#!/usr/bin/env python3
"""Bounded unit tests for the WH2 lazy-arena cross-binary timing harness."""

from __future__ import annotations

import csv
import io
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import wh2_lazy_solve_timing as timing


def one_task(**changes):
    value = {
        "K": 3200, "bb": 64, "seed_index": 0,
        "seed": timing.SCHEDULE_SEEDS[0][1], "schedule": "burst",
        "cache_state": "warm", "job": 0,
        "task_id": "000.K3200.bb64.seed0.burst.warm",
    }
    value.update(changes)
    return value


def grouped_stdout(label="base", task=None, core=8, evict_bytes=4096):
    task = task or one_task()
    schema = "v1" if label == "base" else "v2"
    preamble = timing._expected_preamble(task, evict_bytes)
    preamble.update({
        "schema": schema,
        "control_attempt": "0", "control_matrix_seed": "0x123",
        "control_peel_seed": "0x456", "candidate_attempt": "0",
        "candidate_matrix_seed": "0x123", "candidate_peel_seed": "0x456",
        "preflight_control_result": "0", "preflight_candidate_result": "0",
        "cell_class": "common-success", "common_success": "1",
        "trace_sha256": "a" * 64,
    })
    if label == "candidate":
        preamble.update({
            "solve_value_storage": "owned-noinit",
            "solve_value_publish": "swap",
        })
    keys = timing.PREAMBLE_V1 if label == "base" else timing.PREAMBLE_V2
    lines = [
        "# groupedtiming: " + " ".join(
            "%s=%s" % (key, preamble[key]) for key in keys),
        ",".join(timing.BASE_FIELDS if label == "base" else
                 timing.BASE_FIELDS + timing.LAZY_ONLY_FIELDS),
    ]
    for cycle in range(4):
        for slot, marker in enumerate(timing.INNER_ORDER):
            arm = "control" if marker == "A" else "candidate"
            values = {
                "N": str(task["K"]), "bb": str(task["bb"]),
                "overhead": str(timing.OVERHEAD),
                "schedule": str(task["schedule"]), "seed": str(task["seed"]),
                "loss": timing.LOSS_TEXT,
                "cache_state": str(task["cache_state"]),
                "cycle": str(cycle), "slot": str(slot), "arm": arm,
                "period": str(timing.ARCHITECTURE["period"]),
                "grouped_rows": str(timing.ARCHITECTURE["grouped_rows"]),
                "buckets_requested": str(timing.ARCHITECTURE["buckets"]),
                "seed_attempt": "0", "matrix_seed": "0x123",
                "peel_seed": "0x456", "preflight_result": "0",
                "cell_class": "common-success", "common_success": "1",
                "result": "0", "outcome_stable": "1",
                "elapsed_ns": str(1000 + cycle * 10 + slot),
                "saturated": "0", "cpu_before": str(core),
                "cpu_after": str(core), "cpu_migrated": "0",
                "minflt_delta": "0", "majflt_delta": "0",
                "fault_contaminated": "0", "inactivated": "109",
                "binary_def": "8", "heavy_gain": "8",
                "block_xors": "74095", "block_muladds": "1880",
                "build_ns": "10", "peel_ns": "20", "project_ns": "30",
                "residual_ns": "40", "backsub_ns": "50",
                "joint_source_xors": "0", "joint_marginal_xors": "0",
                "joint_marginal_copies": "0", "joint_active_deltas": "0",
                "joint_scratch_bytes": "0", "dual_source_columns": "0",
                "source_bytes": str(task["K"] * task["bb"]),
                "packet_payload_bytes": str(
                    (task["K"] + timing.OVERHEAD) * task["bb"]),
                "intermediate_bytes": "210304",
                "solve_value_arena_bytes": "210304",
                "solve_value_eager_zero_bytes": "0",
                "solve_value_commit_copy_bytes": "0",
            }
            fields = timing.BASE_FIELDS if label == "base" else \
                timing.BASE_FIELDS + timing.LAZY_ONLY_FIELDS
            lines.append(",".join(values[field] for field in fields))
    return ("\n".join(lines) + "\n").encode("ascii")


def mutate_csv(raw, row_index, field, value):
    lines = raw.decode("ascii").splitlines()
    header = lines[1].split(",")
    fields = lines[row_index + 2].split(",")
    fields[header.index(field)] = str(value)
    lines[row_index + 2] = ",".join(fields)
    return ("\n".join(lines) + "\n").encode("ascii")


def mutate_preamble(raw, key, value):
    lines = raw.decode("ascii").splitlines()
    tokens = lines[0].split(" ")
    for index, token in enumerate(tokens):
        if token.startswith(key + "="):
            tokens[index] = key + "=" + value
            break
    else:
        raise AssertionError("preamble key missing")
    lines[0] = " ".join(tokens)
    return ("\n".join(lines) + "\n").encode("ascii")


def thermal_csv(times=(9.0, 10.0, 11.0, 12.0)):
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=timing.THERMAL_FIELDS,
                            lineterminator="\r\n")
    writer.writeheader()
    for index, monotonic in enumerate(times):
        row = {
            "utc": "2026-07-18T00:00:%02d.000Z" % index,
            "monotonic_s": str(monotonic), "cpu_busy_pct": "1.0",
            "cpu_avg_mhz": "3000.0", "cpu_tctl_c": "55.0",
            "dimm_read_errors": "0", "load1": "0", "load5": "0",
            "load15": "0", "edac_ce": "4", "edac_ue": "0",
        }
        for dimm in timing.DIMM_FIELDS:
            row[dimm] = "48.0"
        writer.writerow(row)
    return output.getvalue().encode("ascii")


class TaskLedgerTests(unittest.TestCase):
    def test_exact_grid_is_deterministic_and_complete(self):
        first = timing.generate_tasks()
        second = timing.generate_tasks()
        self.assertEqual(first, second)
        self.assertEqual(len(first), 108)
        coordinates = {
            (row["K"], row["bb"], row["schedule"], row["seed"],
             row["cache_state"]) for row in first
        }
        expected = {
            (K, width, schedule, seed, cache)
            for K in timing.KS for width in timing.WIDTHS
            for schedule, seed in timing.SCHEDULE_SEEDS
            for cache in timing.CACHE_STATES
        }
        self.assertEqual(coordinates, expected)
        self.assertEqual([row["job"] for row in first], list(range(108)))

    def test_outer_order_balances_four_processes_per_binary(self):
        self.assertEqual(timing.OUTER_ORDER, "ABBABAAB")
        self.assertEqual(timing.OUTER_ORDER.count("A"), 4)
        self.assertEqual(timing.OUTER_ORDER.count("B"), 4)
        self.assertEqual(timing.INNER_ORDER, timing.OUTER_ORDER)
        self.assertEqual(len(timing.BINARY_NAMES["base"]),
                         len(timing.BINARY_NAMES["candidate"]))

    def test_cross_cache_identity_is_exact_and_complete(self):
        ledger = {}
        parsed = timing.parse_grouped_output(
            grouped_stdout("base"), "base", one_task(), 4096, 8)
        for task in timing.generate_tasks():
            timing._register_cross_cache_identity(ledger, task, parsed)
        timing._validate_cross_cache_ledger(ledger)

        mismatched = {}
        cold = one_task(cache_state="cold")
        warm = one_task(cache_state="warm")
        timing._register_cross_cache_identity(mismatched, cold, parsed)
        changed = timing.parse_grouped_output(
            mutate_preamble(grouped_stdout("base", warm),
                            "trace_sha256", "b" * 64),
            "base", warm, 4096, 8)
        with self.assertRaisesRegex(timing.TimingError, "cold/warm"):
            timing._register_cross_cache_identity(mismatched, warm, changed)


class OutputParserTests(unittest.TestCase):
    def test_v1_v2_normalize_to_identical_trace_and_work(self):
        task = one_task()
        base = timing.parse_grouped_output(
            grouped_stdout("base", task), "base", task, 4096, 8)
        candidate = timing.parse_grouped_output(
            grouped_stdout("candidate", task), "candidate", task, 4096, 8)
        self.assertEqual(base.schema, "v1")
        self.assertEqual(candidate.schema, "v2")
        self.assertEqual(base.semantic_sha256, candidate.semantic_sha256)
        self.assertEqual(base.work_signature, candidate.work_signature)
        self.assertEqual(base.preamble["trace_sha256"],
                         candidate.preamble["trace_sha256"])
        self.assertEqual(base.timed_elapsed_ns, candidate.timed_elapsed_ns)
        self.assertEqual(base.contaminations, ())
        self.assertEqual(candidate.contaminations, ())

    def test_lazy_storage_receipts_are_mandatory(self):
        task = one_task()
        raw = mutate_csv(grouped_stdout("candidate", task), 7,
                         "solve_value_eager_zero_bytes", 1)
        with self.assertRaisesRegex(timing.TimingError, "arena receipt changed"):
            timing.parse_grouped_output(raw, "candidate", task, 4096, 8)
        raw = mutate_preamble(grouped_stdout("candidate", task),
                              "solve_value_publish", "copy")
        with self.assertRaisesRegex(timing.TimingError, "owned no-init"):
            timing.parse_grouped_output(raw, "candidate", task, 4096, 8)

    def test_work_or_trace_drift_is_substantive(self):
        task = one_task()
        raw = mutate_csv(grouped_stdout("base", task), 8, "block_xors", 74096)
        with self.assertRaisesRegex(timing.TimingError, "deterministic work"):
            timing.parse_grouped_output(raw, "base", task, 4096, 8)
        raw = mutate_preamble(grouped_stdout("base", task), "trace_sha256", "z" * 64)
        with self.assertRaisesRegex(timing.TimingError, "trace hash"):
            timing.parse_grouped_output(raw, "base", task, 4096, 8)

    def test_migration_saturation_and_faults_are_receipted_contamination(self):
        task = one_task()
        raw = mutate_csv(grouped_stdout("base", task), 9, "cpu_after", 9)
        raw = mutate_csv(raw, 10, "saturated", 1)
        raw = mutate_csv(raw, 11, "minflt_delta", 65)
        raw = mutate_csv(raw, 11, "fault_contaminated", 1)
        raw = mutate_csv(raw, 12, "majflt_delta", 1)
        raw = mutate_csv(raw, 12, "fault_contaminated", 1)
        parsed = timing.parse_grouped_output(raw, "base", task, 4096, 8)
        self.assertTrue(any("migration" in value for value in parsed.contaminations))
        self.assertTrue(any("saturated" in value for value in parsed.contaminations))
        self.assertTrue(any("minor-fault" in value for value in parsed.contaminations))
        self.assertTrue(any("major-fault" in value for value in parsed.contaminations))

    def test_signed_fault_receipt_cannot_be_forged(self):
        task = one_task()
        raw = mutate_csv(grouped_stdout("base", task), 0, "minflt_delta", -1)
        with self.assertRaisesRegex(timing.TimingError, "fault receipt disagrees"):
            timing.parse_grouped_output(raw, "base", task, 4096, 8)
        raw = mutate_csv(raw, 0, "fault_contaminated", -1)
        parsed = timing.parse_grouped_output(raw, "base", task, 4096, 8)
        self.assertIn("row0:minor-fault:-1", parsed.contaminations)

    def test_schema_and_preamble_order_are_exact(self):
        task = one_task()
        raw = mutate_preamble(grouped_stdout("base", task), "schema", "v2")
        with self.assertRaisesRegex(timing.TimingError, "wrong groupedtiming schema"):
            timing.parse_grouped_output(raw, "base", task, 4096, 8)
        lines = grouped_stdout("base", task).decode("ascii").splitlines()
        tokens = lines[0].split(" ")
        tokens[-1], tokens[-2] = tokens[-2], tokens[-1]
        lines[0] = " ".join(tokens)
        with self.assertRaisesRegex(timing.TimingError, "order/schema"):
            timing.parse_grouped_output(
                ("\n".join(lines) + "\n").encode("ascii"),
                "base", task, 4096, 8)


class ReceiptAndThermalTests(unittest.TestCase):
    def test_sealed_record_detects_mutation(self):
        value = timing.sealed_record("example.v1", {"number": 4})
        timing.verify_sealed_record(value, "example.v1", "fixture")
        value["number"] = 5
        with self.assertRaisesRegex(timing.TimingError, "self-hash mismatch"):
            timing.verify_sealed_record(value, "example.v1", "fixture")

    def test_execution_receipt_schema_matches_constructor(self):
        task = one_task()
        parsed = timing.parse_grouped_output(
            grouped_stdout("base", task), "base", task, 4096, 8)
        receipt = timing._execution_receipt(
            task, 0, "base", ["example"], 0,
            "2026-07-18T00:00:00.000Z", 10, 20, parsed, b"", [],
            "c" * 64)
        self.assertEqual(set(receipt), timing.EXECUTION_RECEIPT_FIELDS)

    def test_thermal_interval_is_exactly_bracketed_and_bounded(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "thermal.csv"
            path.write_bytes(thermal_csv())
            raw, summary = timing.collect_thermal_interval(path, 10.25, 11.75)
        rows, _lines = timing._parse_thermal(raw)
        self.assertEqual([float(row["monotonic_s"]) for row in rows],
                         [10.0, 11.0, 12.0])
        self.assertEqual(summary["sample_count"], 3)
        self.assertEqual(summary["dimm_read_errors"], 0)
        self.assertEqual(summary["edac_ce_delta"], 0)

    def test_thermal_snapshot_retries_a_concurrent_partial_append(self):
        complete = thermal_csv()
        with mock.patch.object(
                timing, "stable_bytes",
                side_effect=[complete[:-3], complete]), \
                mock.patch.object(timing.time, "sleep"):
            rows, lines, raw = timing._stable_thermal_snapshot(
                Path("unused"), timeout_s=1.0)
        self.assertEqual(raw, complete)
        self.assertEqual(len(rows) + 1, len(lines))

    def test_thermal_gap_and_edac_change_fail_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "thermal.csv"
            path.write_bytes(thermal_csv((9.0, 10.0, 12.5)))
            with self.assertRaisesRegex(timing.TimingError, "coverage or cadence"):
                timing.collect_thermal_interval(path, 9.5, 12.0)
            text = thermal_csv().decode("ascii").splitlines()
            header = text[0].split(",")
            fields = text[-1].split(",")
            fields[header.index("edac_ce")] = "5"
            text[-1] = ",".join(fields)
            path.write_bytes(("\r\n".join(text) + "\r\n").encode("ascii"))
            with self.assertRaisesRegex(timing.TimingError, "EDAC counters"):
                timing.collect_thermal_interval(path, 10.25, 11.75)

    def test_bootstrap_is_deterministic(self):
        rows = [(100, 90), (200, 205), (300, 270), (400, 390)]
        first = timing._bootstrap(rows, "unit-test", repetitions=200)
        second = timing._bootstrap(rows, "unit-test", repetitions=200)
        self.assertEqual(first, second)
        self.assertLess(float(first["lower_95"]), float(first["upper_95"]))


if __name__ == "__main__":
    unittest.main()
