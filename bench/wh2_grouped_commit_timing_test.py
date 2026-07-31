#!/usr/bin/env python3
"""Bounded unit tests for the WH2 grouped-commit cross-binary timing harness."""

from __future__ import annotations

import csv
import io
import os
from pathlib import Path
import tempfile
import time
import unittest
from unittest import mock

import wh2_grouped_commit_timing as timing


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
    if label not in ("base", "candidate"):
        raise AssertionError("unknown fixture label")
    schema = "v1"
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
    keys = timing.GROUPED_PREAMBLE_FIELDS
    lines = [
        "# groupedtiming: " + " ".join(
            "%s=%s" % (key, preamble[key]) for key in keys),
        ",".join(timing.GROUPED_FIELDS),
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
            }
            lines.append(",".join(values[field] for field in timing.GROUPED_FIELDS))
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
                            lineterminator="\n")
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
    def test_v1_v1_normalize_to_identical_trace_work_and_timing(self):
        task = one_task()
        base = timing.parse_grouped_output(
            grouped_stdout("base", task), "base", task, 4096, 8)
        candidate = timing.parse_grouped_output(
            grouped_stdout("candidate", task), "candidate", task, 4096, 8)
        self.assertEqual(base.schema, "v1")
        self.assertEqual(candidate.schema, "v1")
        self.assertEqual(base.semantic_sha256, candidate.semantic_sha256)
        self.assertEqual(base.work_signature, candidate.work_signature)
        self.assertEqual(base.preamble["trace_sha256"],
                         candidate.preamble["trace_sha256"])
        self.assertEqual(base.timed_elapsed_ns, candidate.timed_elapsed_ns)
        self.assertEqual(base.timed_phase_ns, {
            "build_ns": 240, "peel_ns": 480, "project_ns": 720,
            "residual_ns": 960, "backsub_ns": 1200,
        })
        self.assertEqual(base.all_phase_ns, {
            "build_ns": 320, "peel_ns": 640, "project_ns": 960,
            "residual_ns": 1280, "backsub_ns": 1600,
        })
        self.assertEqual(base.timed_phase_ns, candidate.timed_phase_ns)
        self.assertEqual(base.contaminations, ())
        self.assertEqual(candidate.contaminations, ())

    def test_both_binary_labels_require_the_same_v1_schema(self):
        task = one_task()
        raw = mutate_preamble(grouped_stdout("candidate", task),
                              "schema", "v2")
        with self.assertRaisesRegex(timing.TimingError,
                                    "wrong groupedtiming schema"):
            timing.parse_grouped_output(raw, "candidate", task, 4096, 8)

    def test_phase_timing_changes_are_measured_not_semantic_drift(self):
        task = one_task()
        base = timing.parse_grouped_output(
            grouped_stdout("base", task), "base", task, 4096, 8)
        changed = timing.parse_grouped_output(
            mutate_csv(grouped_stdout("candidate", task), 8,
                       "residual_ns", 41),
            "candidate", task, 4096, 8)
        self.assertEqual(base.semantic_sha256, changed.semantic_sha256)
        self.assertEqual(base.work_signature, changed.work_signature)
        self.assertEqual(changed.timed_phase_ns["residual_ns"],
                         base.timed_phase_ns["residual_ns"] + 1)

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

    def test_csv_quotes_extra_cells_and_oversize_are_rejected(self):
        task = one_task()
        quoted = grouped_stdout("base", task).replace(
            b"burst", b'"burst"', 1)
        with self.assertRaisesRegex(timing.TimingError, "canonical LF"):
            timing.parse_grouped_output(quoted, "base", task, 4096, 8)
        lines = grouped_stdout("base", task).decode("ascii").splitlines()
        lines[2] += ",surplus"
        with self.assertRaisesRegex(timing.TimingError, "arity"):
            timing.parse_grouped_output(
                ("\n".join(lines) + "\n").encode("ascii"),
                "base", task, 4096, 8)
        with self.assertRaisesRegex(timing.TimingError, "canonical LF"):
            timing.parse_grouped_output(
                b"x" * timing.MAX_GROUPED_OUTPUT_BYTES + b"\n",
                "base", task, 4096, 8)

    def test_impossible_domains_and_phase_overrun_are_rejected(self):
        task = one_task()
        for field, value, error in (
                ("cpu_before", timing.MAX_CPU_ID + 1, "integer domain"),
                ("minflt_delta", -2, "integer domain"),
                ("majflt_delta", -2, "integer domain"),
                ("heavy_gain", 7, "rank/work")):
            raw = mutate_csv(
                grouped_stdout("base", task), 0, field, value)
            with self.assertRaisesRegex(timing.TimingError, error):
                timing.parse_grouped_output(raw, "base", task, 4096, 8)
        raw = mutate_csv(
            grouped_stdout("base", task), 0, "residual_ns", 1000)
        with self.assertRaisesRegex(timing.TimingError, "phases exceed"):
            timing.parse_grouped_output(raw, "base", task, 4096, 8)

    def test_zero_negative_control_phase_is_valid(self):
        task = one_task()
        raw = grouped_stdout("base", task)
        for row_index in range(32):
            raw = mutate_csv(raw, row_index, "build_ns", 0)
        parsed = timing.parse_grouped_output(
            raw, "base", task, 4096, 8)
        self.assertEqual(parsed.all_phase_ns["build_ns"], 0)
        self.assertEqual(
            timing._ratio_record_or_none(0, 0), None)


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
            "c" * 64, {
                "pid": 123, "start_ticks": 456,
                "executable": {
                    "path": "/example", "device": 7, "inode": 8,
                },
                "argv": ["/example"], "boot_id":
                    "1788608a-7aa1-48de-8f7c-848485be3cc3",
                "binding_observation": "double-proc-snapshot",
            }, "exited_group_swept")
        self.assertEqual(set(receipt), timing.EXECUTION_RECEIPT_FIELDS)
        self.assertEqual(receipt["timed_phase_ns"], parsed.timed_phase_ns)
        self.assertEqual(receipt["all_phase_ns"], parsed.all_phase_ns)

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
            path.write_bytes(("\n".join(text) + "\n").encode("ascii"))
            with self.assertRaisesRegex(timing.TimingError, "EDAC counters"):
                timing.collect_thermal_interval(path, 10.25, 11.75)

    def test_thermal_numeric_dates_and_transient_edac_are_strict(self):
        def mutate(raw, row_index, field, value):
            text = raw.decode("ascii").splitlines()
            header = text[0].split(",")
            fields = text[row_index + 1].split(",")
            fields[header.index(field)] = value
            text[row_index + 1] = ",".join(fields)
            return ("\n".join(text) + "\n").encode("ascii")

        raw = mutate(thermal_csv(), 0, "utc",
                     "2026-02-31T00:00:00.000Z")
        with self.assertRaisesRegex(timing.TimingError, "UTC date"):
            timing._parse_thermal(raw)
        raw = mutate(thermal_csv(), 0, "cpu_busy_pct", "1e2")
        with self.assertRaisesRegex(timing.TimingError, "canonical decimal"):
            timing._parse_thermal(raw)
        raw = mutate(thermal_csv(), 1, "edac_ce", "5")
        with self.assertRaisesRegex(timing.TimingError, "EDAC counters"):
            timing.validate_sealed_thermal_interval(raw, 9.5, 11.5)

    def test_bootstrap_is_deterministic(self):
        rows = [(100, 90), (200, 205), (300, 270), (400, 390)]
        first = timing._bootstrap(rows, "unit-test", repetitions=200)
        second = timing._bootstrap(rows, "unit-test", repetitions=200)
        self.assertEqual(first, second)
        self.assertLess(float(first["lower_95"]), float(first["upper_95"]))

    def test_zero_denominator_phase_summary_is_explicit(self):
        unavailable = timing._summarize(
            [{"base": 0, "candidate": 0}], "zero", "build_ns")
        self.assertFalse(unavailable["ratio_available"])
        self.assertIsNone(unavailable["ratio_of_sums"])
        available = timing._summarize(
            [{"base": 0, "candidate": 4},
             {"base": 10, "candidate": 9}],
            "mixed-zero", "build_ns")
        self.assertTrue(available["ratio_available"])
        self.assertEqual(available["base_zero_tasks"], 1)

    def test_phase_breakdowns_cover_the_frozen_grid(self):
        records = []
        for task in timing.generate_tasks():
            records.append({
                **task,
                "timings": {
                    field: {"base": 100 + task["job"],
                            "candidate": 101 + task["job"]}
                    for field in ("elapsed_ns", *timing.PHASE_FIELDS)
                },
            })

        def summarize(rows, domain, metric):
            return {
                "task_count": len(rows), "domain": domain, "metric": metric,
            }

        with mock.patch.object(timing, "_summarize",
                               side_effect=summarize):
            result = timing._metric_breakdowns(records, "residual_ns")
        self.assertEqual(result["overall"]["task_count"], 108)
        self.assertEqual(
            {value["task_count"]
             for value in result["by_cache_state"].values()}, {54})
        self.assertEqual(
            {value["task_count"] for value in result["by_K"].values()}, {18})
        self.assertEqual(
            {value["task_count"] for value in result["by_bb"].values()}, {36})
        self.assertEqual(
            {value["task_count"]
             for value in result["by_schedule"].values()}, {36})
        self.assertTrue(
            result["overall"]["domain"].startswith("metric:residual_ns:"))
        self.assertEqual(result["overall"]["metric"], "residual_ns")


class BoundCommandTests(unittest.TestCase):
    python = Path("/usr/bin/python3.12")
    bash = Path("/usr/bin/bash")
    sleep = Path("/usr/bin/sleep")
    env = Path("/usr/bin/env")
    true = Path("/usr/bin/true")

    @classmethod
    def setUpClass(cls):
        if not all(path.is_file() for path in (
                cls.python, cls.bash, cls.sleep, cls.env, cls.true)):
            raise unittest.SkipTest("Linux lifecycle test tools are unavailable")
        try:
            timing.require_pidfd_runtime()
        except timing.TimingError as exc:
            raise unittest.SkipTest(str(exc))

    def test_exact_final_exec_is_bound_and_output_is_bounded(self):
        command = [
            str(self.bash), "-c",
            "printf stdout; printf stderr >&2; sleep 0.02",
        ]
        result = timing.run_bound_command(
            command, self.bash, self.python, 2.0)
        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout, b"stdout")
        self.assertEqual(result.stderr, b"stderr")
        self.assertEqual(result.cleanup_action, "exited_group_swept")
        self.assertEqual(result.process_identity["argv"], command)
        self.assertEqual(
            (result.process_identity["executable"]["device"],
             result.process_identity["executable"]["inode"]),
            (self.bash.stat().st_dev, self.bash.stat().st_ino))

    def test_prefixed_fast_final_exec_is_proven_repeatedly(self):
        command = [str(self.env), "-i", str(self.true)]
        for _attempt in range(20):
            result = timing.run_bound_command(
                command, self.true, self.python, 1.0)
            self.assertEqual(result.returncode, 0)
            self.assertEqual(result.stdout, b"")
            self.assertEqual(result.stderr, b"")
            self.assertIn(
                result.process_identity["binding_observation"],
                ("double-proc-snapshot",
                 "final-exec-handshake-terminal"))

    def test_prefixed_fast_nonzero_output_is_retained(self):
        command = [
            str(self.env), "-i", str(self.bash), "-c",
            "printf quick-out; printf quick-err >&2; exit 7",
        ]
        result = timing.run_bound_command(
            command, self.bash, self.python, 1.0)
        self.assertEqual(result.returncode, 7)
        self.assertEqual(result.stdout, b"quick-out")
        self.assertEqual(result.stderr, b"quick-err")

    def test_marker_then_wrapper_signal_is_not_accepted_as_exec(self):
        dying_wrapper = r"""
import os
import signal
import sys
fd = int(sys.argv[1])
os.write(fd, b"wirehair-final-exec-v1\n")
os.kill(os.getpid(), signal.SIGKILL)
"""
        with mock.patch.object(
                timing, "FINAL_EXEC_WRAPPER", dying_wrapper):
            with self.assertRaisesRegex(
                    timing.BoundCommandError,
                    "terminal-only final-exec proof") as raised:
                timing.run_bound_command(
                    [str(self.true)], self.true, self.python, 1.0)
        self.assertEqual(raised.exception.returncode, -9)
        self.assertEqual(
            raised.exception.process_identity["binding_observation"],
            "final-exec-handshake-terminal")

    def test_prefix_failure_retains_bounded_stderr_and_returncode(self):
        command = [
            str(self.env), "--definitely-not-an-env-option", str(self.true),
        ]
        with self.assertRaisesRegex(
                timing.BoundCommandError,
                "final-exec handshake failed") as raised:
            timing.run_bound_command(
                command, self.true, self.python, 1.0)
        error = raised.exception
        self.assertNotEqual(error.stderr, b"")
        self.assertIsNotNone(error.returncode)
        self.assertIsNotNone(error.process_identity)
        self.assertEqual(
            error.process_identity["binding_observation"],
            "direct-child-provisional")

    def test_timeout_kills_the_bound_direct_child(self):
        started = time.monotonic()
        with self.assertRaisesRegex(
                timing.BoundCommandError, "exceeded its timeout") as raised:
            timing.run_bound_command(
                [str(self.sleep), "30"], self.sleep, self.python, 0.05)
        self.assertLess(time.monotonic() - started, 1.0)
        error = raised.exception
        self.assertIn(
            error.cleanup_action, ("terminated_group", "killed_group"))
        self.assertIsNotNone(error.process_identity)
        pid = int(error.process_identity["pid"])
        self.assertIsNone(timing.process_start_ticks(pid))

    def test_stdout_overflow_kills_child_without_unbounded_capture(self):
        command = [
            str(self.bash), "-c",
            "sleep 0.02; head -c %d /dev/zero; sleep 30" %
            (timing.MAX_GROUPED_OUTPUT_BYTES + 1),
        ]
        with self.assertRaisesRegex(
                timing.BoundCommandError, "exceeded its bounded capture") as raised:
            timing.run_bound_command(command, self.bash, self.python, 2.0)
        error = raised.exception
        self.assertLessEqual(
            len(error.stdout), timing.MAX_GROUPED_OUTPUT_BYTES)
        self.assertLessEqual(
            len(error.stderr), timing.MAX_BENCHMARK_STDERR_BYTES)
        self.assertIn(
            error.cleanup_action, ("terminated_group", "killed_group"))
        self.assertIsNone(
            timing.process_start_ticks(int(error.process_identity["pid"])))

    def test_stdout_overflow_evidence_is_the_exact_stream_prefix(self):
        prefix_size = timing.MAX_GROUPED_OUTPUT_BYTES - 10
        program = (
            "import os,time;"
            "os.write(1,b'A'*%d);time.sleep(.05);"
            "os.write(1,b'B'*20);time.sleep(30)" % prefix_size
        )
        command = [str(self.python), "-I", "-c", program]
        with self.assertRaisesRegex(
                timing.BoundCommandError,
                "exceeded its bounded capture") as raised:
            timing.run_bound_command(
                command, self.python, self.python, 2.0)
        self.assertEqual(
            raised.exception.stdout,
            b"A" * prefix_size + b"B" * 10)

    def test_timeout_kills_descendants_in_the_owned_session(self):
        command = [
            str(self.bash), "-c",
            "sleep 30 & child=$!; printf '%s\\n' \"$child\"; sleep 30",
        ]
        with self.assertRaisesRegex(
                timing.BoundCommandError, "exceeded its timeout") as raised:
            timing.run_bound_command(command, self.bash, self.python, 0.10)
        error = raised.exception
        child = int(error.stdout.strip())
        self.assertIn(
            error.cleanup_action, ("terminated_group", "killed_group"))
        self.assertIsNone(timing.process_start_ticks(child))

    def test_success_sweeps_background_descendant_with_closed_stdio(self):
        command = [
            str(self.bash), "-c",
            "sleep 30 </dev/null >/dev/null 2>&1 & printf '%s\\n' \"$!\"",
        ]
        result = timing.run_bound_command(
            command, self.bash, self.python, 1.0)
        child = int(result.stdout.strip())
        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.cleanup_action, "exited_group_swept")
        live = [
            item for item in timing._process_group_members(
                int(result.process_identity["pid"]))
            if item[1] != "Z"
        ]
        self.assertEqual(live, [])
        if timing.process_start_ticks(child) is not None:
            state = (Path("/proc") / str(child) / "stat").read_bytes(
            ).rsplit(b") ", 1)[1].split()[0]
            self.assertEqual(state, b"Z")

    def test_timeout_overflow_is_rejected_without_fd_leak(self):
        before = len(os.listdir("/proc/self/fd"))
        for _attempt in range(8):
            with self.assertRaisesRegex(
                    timing.TimingError, "bounded domain|finite number"):
                timing.run_bound_command(
                    [str(self.sleep), "1"], self.sleep, self.python, 10 ** 400)
        self.assertEqual(len(os.listdir("/proc/self/fd")), before)

    def test_cleanup_failure_is_attached_without_masking_primary(self):
        original = timing.terminate_direct_child_by_pidfd

        def cleanup_then_fail(process, pidfd):
            original(process, pidfd)
            raise timing.TimingError("injected cleanup report")

        with mock.patch.object(
                timing, "terminate_direct_child_by_pidfd",
                side_effect=cleanup_then_fail):
            with self.assertRaisesRegex(
                    timing.BoundCommandError,
                    "exceeded its timeout") as raised:
                timing.run_bound_command(
                    [str(self.sleep), "30"],
                    self.sleep, self.python, 0.05)
        self.assertEqual(raised.exception.cleanup_action, "cleanup_failed")
        self.assertIn(
            "injected cleanup report",
            str(raised.exception.cleanup_error))

    def test_drain_failure_is_attached_without_masking_primary(self):
        with mock.patch.object(
                timing, "_drain_after_cleanup",
                side_effect=timing.TimingError("injected drain report")):
            with self.assertRaisesRegex(
                    timing.BoundCommandError,
                    "exceeded its timeout") as raised:
                timing.run_bound_command(
                    [str(self.sleep), "30"],
                    self.sleep, self.python, 0.05)
        self.assertIn(
            "injected drain report",
            str(raised.exception.cleanup_error))


class ReductionPublicationTests(unittest.TestCase):
    def test_exact_partial_reduction_is_resumable(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            summary = root / "validated_summary.json"
            manifest = root / "data_manifest.json"
            sidecar = root / "data_manifest.sha256"
            self.assertEqual(
                timing._publish_or_verify_reduction_artifact(
                    summary, b'{"summary":1}\n'),
                "published")

            # Model a controller failure before the remaining publications.
            self.assertEqual(
                timing._publish_or_verify_reduction_artifact(
                    summary, b'{"summary":1}\n'),
                "verified_existing")
            self.assertEqual(
                timing._publish_or_verify_reduction_artifact(
                    manifest, b'{"files":[]}\n'),
                "published")
            self.assertEqual(
                timing._publish_or_verify_reduction_artifact(
                    sidecar, b"a" * 64 + b"\n"),
                "published")

    def test_partial_reduction_must_be_exact_and_immutable(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / "validated_summary.json"
            timing._publish_or_verify_reduction_artifact(target, b"first\n")
            with self.assertRaisesRegex(
                    timing.TimingError, "does not reproduce"):
                timing._publish_or_verify_reduction_artifact(
                    target, b"second\n")
            target.chmod(0o644)
            with self.assertRaisesRegex(
                    timing.TimingError, "not immutable"):
                timing._publish_or_verify_reduction_artifact(
                    target, b"first\n")


class OverlayProvenanceTests(unittest.TestCase):
    def overlay_record(self):
        return {
            "source_commit": timing.MEASUREMENT_OVERLAY_COMMIT,
            "source_parent_commit":
                timing.MEASUREMENT_OVERLAY_PARENT_COMMIT,
            "source_tree": "1" * 40,
            "files": [
                {
                    "path": path,
                    "parent_blob_oid": str(index + 2) * 40,
                    "parent_sha256": str(index + 4) * 64,
                    "overlay_blob_oid": str(index + 6) * 40,
                    "overlay_sha256": str(index + 8) * 64,
                }
                for index, path in enumerate(timing.MEASUREMENT_OVERLAY_FILES)
            ],
            "diff_options": list(timing.MEASUREMENT_DIFF_OPTIONS),
            "diff_sha256": "a" * 64,
            "stable_patch_id": "b" * 40,
            "diff_evidence_name": "base.measurement-overlay.diff",
        }

    def test_exact_commit_pair_and_overlay_are_frozen(self):
        self.assertEqual(
            timing.BASE_COMMIT,
            "48d14bc77e3f9e98605fca4d226aa218d7d03a0d")
        self.assertEqual(
            timing.CANDIDATE_COMMIT,
            "c7203519b4ef42a3d5b7bd5073152a04f89eb9d3")
        self.assertEqual(
            timing.MEASUREMENT_OVERLAY_COMMIT,
            "3eb1eaf41ded5031393ac84200a62e3c0a0b5456")
        self.assertEqual(
            timing.MEASUREMENT_OVERLAY_PARENT_COMMIT,
            "243f8ed86b7bf102fa1cb7156a481c170935e57b")
        self.assertIn(
            "WIREHAIR_V2_GROUPED_TIMING_ONLY",
            timing.TIMING_CXX_FLAGS_RELEASE)
        self.assertIn(
            "WIREHAIR_V2_DISABLE_PACKED_RESIDUAL_TEXT_SECTION",
            timing.TIMING_CXX_FLAGS_RELEASE)
        identity = timing._overlay_identity(self.overlay_record())
        self.assertNotIn("diff_evidence_name", identity)
        self.assertEqual(
            [item["path"] for item in identity["files"]],
            list(timing.MEASUREMENT_OVERLAY_FILES))

    def test_overlay_identity_rejects_noop_or_file_drift(self):
        record = self.overlay_record()
        record["files"][0]["overlay_sha256"] = \
            record["files"][0]["parent_sha256"]
        with self.assertRaisesRegex(timing.TimingError,
                                    "did not change a listed file"):
            timing._overlay_identity(record)
        record = self.overlay_record()
        record["files"].reverse()
        with self.assertRaisesRegex(timing.TimingError, "file list changed"):
            timing._overlay_identity(record)

    def test_stable_patch_id_is_parsed_strictly(self):
        completed = mock.Mock(
            returncode=0, stderr=b"",
            stdout=(b"c" * 40) + b" " + (b"0" * 40) + b"\n")
        with mock.patch.object(timing.subprocess, "run",
                               return_value=completed):
            self.assertEqual(
                timing._stable_patch_id(Path("/usr/bin/git"), b"diff\n"),
                "c" * 40)

    def test_porcelain_status_preserves_leading_index_column(self):
        completed = mock.Mock(
            returncode=0, stderr=b"",
            stdout=b" M codec/V2BenchCliTest.cmake\n"
                   b" M codec/WirehairV2Bench.cpp\n")
        with mock.patch.object(timing.subprocess, "run",
                               return_value=completed):
            rows = timing._git_status_lines(
                Path("/usr/bin/git"), Path("/unused"))
        self.assertEqual(set(rows), timing._expected_overlay_status())


if __name__ == "__main__":
    unittest.main()
