#!/usr/bin/env python3
"""Bounded unit tests for the WH2 grouped-commit cross-binary timing harness."""

from __future__ import annotations

import csv
import copy
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


def schedstat_text(records=((0, 1000, 10), (1, 2000, 20))):
    lines = ["version 15", "timestamp 123456"]
    for cpu, runtime_ns, pcount in records:
        # Fields 1 and 3-6 are legitimate nonzero v15 counters.  Field 2 is
        # the sole retired field and must remain zero.
        values = (7, 0, 8, 9, 10, 11, runtime_ns, 12, pcount)
        lines.append("cpu%d %s" % (
            cpu, " ".join(str(value) for value in values)))
        lines.append("domain0 00000001 " + " ".join(
            "0" for _index in range(timing.SCHEDSTAT_DOMAIN_FIELD_COUNT)))
    return ("\n".join(lines) + "\n").encode("ascii")


def isolation_snapshot(start_ns=100, end_ns=120):
    expected = [6, 70, 80]
    affinity_records = [
        {"irq": timing.ZERO_COUNT_GUARDED_IRQ,
         "effective_affinity_list": "6", "effective_cpus": [6]},
        {"irq": 31, "effective_affinity_list": "0-5",
         "effective_cpus": list(range(6))},
    ]
    managed = []
    for irq, identity, handler, requested, effective in \
            timing.MANAGED_NVME_IRQ_WHITELIST:
        affinity_records.append({
            "irq": irq, "effective_affinity_list": effective,
            "effective_cpus": list(timing.parse_cpu_list(effective)),
        })
        managed.append({
            "irq": irq, "identity": identity,
            "handler_directories": [handler],
            "requested_affinity_list": requested,
            "requested_cpus": list(timing.parse_cpu_list(requested)),
            "effective_affinity_list": effective,
            "effective_cpus": list(timing.parse_cpu_list(effective)),
        })
    affinity_records.sort(key=lambda value: value["irq"])
    return {
        "schema": "wirehair.wh2.runtime_isolation_snapshot.v3",
        "capture_start_monotonic_ns": start_ns,
        "capture_end_monotonic_ns": end_ns,
        "capture_duration_ns": end_ns - start_ns,
        "self_cgroup": "/wh2-timing-v4",
        "expected_isolated_cpus": expected,
        "kernel_isolated_cpu_list": "6,70,80",
        "kernel_isolated_cpus": expected,
        "cgroup_cpu_list": "6,70,80", "cgroup_cpus": expected,
        "cgroup_effective_cpu_list": "6,70,80",
        "cgroup_effective_cpus": expected,
        "cgroup_exclusive_cpu_list": "6,70,80",
        "cgroup_exclusive_cpus": expected,
        "cgroup_exclusive_effective_cpu_list": "6,70,80",
        "cgroup_exclusive_effective_cpus": expected,
        "cgroup_partition": "isolated",
        "irq_effective_affinities": affinity_records,
        "zero_count_irq_exception": {
            "irq": timing.ZERO_COUNT_GUARDED_IRQ,
            "identity": timing.ZERO_COUNT_GUARDED_IRQ_IDENTITY,
            "handler_directories": [timing.ZERO_COUNT_GUARDED_IRQ_HANDLER],
            "requested_affinity_list": "0-127",
            "requested_cpus": list(range(128)),
            "effective_affinity_list": "6", "effective_cpus": [6],
            "global_interrupt_count": 0,
        },
        "managed_nvme_exceptions": managed,
    }


def target_irq_snapshot(
        hard=None, soft=None, global_hard=None,
        start_ns=100, end_ns=110, cpu_ids=(8, 72, 126)):
    hard = hard or {}
    soft = soft or {}
    global_hard = global_hard or {}
    hard_vectors = ((str(timing.ZERO_COUNT_GUARDED_IRQ),
                     timing.ZERO_COUNT_GUARDED_IRQ_IDENTITY),
                    ("142", timing.GUARDED_IRQ_IDENTITIES[142]),
                    ("NMI", "Non-maskable interrupts"),
                    ("LOC", "Local timer interrupts"))
    interrupts = [" " + " ".join("CPU%d" % cpu for cpu in cpu_ids)]
    for vector, suffix in hard_vectors:
        counts = hard.get(vector, (0, 0, 0))
        by_cpu = dict(zip((8, 72, 126), counts))
        interrupts.append(" %s: %s %s" % (
            vector, " ".join(str(by_cpu.get(cpu, 0)) for cpu in cpu_ids),
            suffix))
    for vector in ("ERR", "MIS"):
        interrupts.append(" %s: %d" % (vector, global_hard.get(vector, 0)))
    softirqs = [" " + " ".join("CPU%d" % cpu for cpu in cpu_ids)]
    for vector in timing.EXPECTED_SOFTIRQ_VECTORS:
        counts = soft.get(vector, (0, 0, 0))
        by_cpu = dict(zip((8, 72, 126), counts))
        softirqs.append(" %s: %s" % (
            vector, " ".join(str(by_cpu.get(cpu, 0)) for cpu in cpu_ids)))
    return timing.parse_target_irq_snapshot(
        ("\n".join(interrupts) + "\n").encode("ascii"),
        ("\n".join(softirqs) + "\n").encode("ascii"),
        (8, 72, 126), start_ns, end_ns)


def zero_irq_arguments(start_ns=10, end_ns=None):
    if end_ns is None:
        end_ns = start_ns + timing.SIBLING_ACCEPTED_EXECUTION_MIN_DURATION_NS
    step = (end_ns - start_ns) // 5
    if step <= 0:
        raise AssertionError("fixture IRQ interval is too short")
    before = target_irq_snapshot(
        start_ns=start_ns + step, end_ns=start_ns + 2 * step)
    after = target_irq_snapshot(
        start_ns=start_ns + 3 * step, end_ns=start_ns + 4 * step)
    return before, after, (8, 72, 126)


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

    def test_tmpfs_binding_is_live_exact_and_tamper_evident(self):
        with tempfile.TemporaryDirectory(dir="/dev/shm") as temporary:
            root = Path(temporary).resolve()
            sample = root / "sample"
            sample.write_bytes(b"fixture")
            directory = timing.capture_tmpfs_binding(
                root, "directory", "fixture directory")
            regular = timing.capture_tmpfs_binding(
                sample, "regular", "fixture file")
            self.assertEqual(directory["filesystem_magic"], timing.TMPFS_MAGIC)
            self.assertEqual(directory["device"], regular["device"])
            timing.validate_tmpfs_binding(
                directory, root, "directory", "fixture directory")
            for field, value in (
                    ("filesystem_magic", 0),
                    ("device", regular["device"] + 1),
                    ("inode", regular["inode"] + 1),
                    ("mode", regular["mode"] ^ 0o100),
                    ("uid", regular["uid"] + 1),
                    ("gid", regular["gid"] + 1),
                    ("nlink", 2)):
                tampered = {**regular, field: value}
                with self.subTest(tmpfs_field=field):
                    with self.assertRaises(timing.TimingError):
                        timing.validate_tmpfs_binding(
                            tampered, sample, "regular", "fixture file")
        with tempfile.TemporaryDirectory(dir="/tmp") as temporary:
            with self.assertRaisesRegex(
                    timing.TimingError, "must reside on tmpfs"):
                timing.capture_tmpfs_binding(
                    Path(temporary).resolve(), "directory",
                    "non-tmpfs fixture")

    def test_prepared_tmpfs_ledger_covers_every_input_and_output_directory(self):
        with tempfile.TemporaryDirectory(dir="/dev/shm") as temporary:
            root = Path(temporary).resolve()
            for relative in timing.PREPARED_CAMPAIGN_DIRECTORIES:
                (root / relative).mkdir(mode=0o700)
            immutable_name = "frozen/fixture"
            for relative in (
                    "design.json", "prepare_receipt.json",
                    "tasks_manifest.jsonl", immutable_name):
                path = root / relative
                path.write_bytes(b"fixture")
            design = {"immutable_files": {immutable_name: "0" * 64}}
            bindings = timing.capture_prepared_tree_tmpfs_bindings(root, design)
            self.assertEqual(
                len(bindings),
                1 + len(timing.PREPARED_CAMPAIGN_DIRECTORIES) + 4)
            self.assertEqual(
                timing.validate_prepared_tree_tmpfs_bindings(
                    bindings, root, design), bindings)
            timing.require_live_tmpfs_tree(
                root, int(bindings[0]["device"]), "fixture")

            real_capture = timing.capture_tmpfs_binding

            def forged_device(path, kind, context):
                binding = real_capture(path, kind, context)
                if Path(path).name == "raw":
                    return {**binding, "device": int(binding["device"]) + 1}
                return binding

            with mock.patch.object(
                    timing, "capture_tmpfs_binding",
                    side_effect=forged_device):
                with self.assertRaisesRegex(
                        timing.TimingError, "spans filesystems"):
                    timing.capture_prepared_tree_tmpfs_bindings(root, design)
            link = root / "nested-link"
            link.symlink_to(root / "design.json")
            with self.assertRaisesRegex(timing.TimingError, "symlink"):
                timing.require_live_tmpfs_tree(
                    root, int(bindings[0]["device"]), "fixture")

    def test_target_irq_capture_intervals_are_ordered_and_bound(self):
        before = target_irq_snapshot(start_ns=100, end_ns=110)
        after = target_irq_snapshot(start_ns=120, end_ns=130)
        timing.checked_target_irq_delta(
            before, after, (8, 72, 126), "fixture")
        timing.require_target_irq_contained_interval(
            before, after, 90, 140, "fixture")
        with self.assertRaisesRegex(timing.TimingError, "capture interval"):
            timing.require_target_irq_contained_interval(
                before, after, 105, 140, "fixture")
        with self.assertRaisesRegex(timing.TimingError, "topology changed"):
            timing.checked_target_irq_delta(
                before,
                target_irq_snapshot(start_ns=105, end_ns=115),
                (8, 72, 126), "fixture")
        outer_before = target_irq_snapshot(start_ns=70, end_ns=80)
        outer_after = target_irq_snapshot(start_ns=150, end_ns=160)
        timing.require_target_irq_bracketing_interval(
            outer_before, outer_after, 90, 140, "fixture")
        with self.assertRaisesRegex(timing.TimingError, "capture bracket"):
            timing.require_target_irq_bracketing_interval(
                before, outer_after, 105, 140, "fixture")

    def test_target_irq_contamination_is_never_retryable(self):
        task = one_task()
        parsed = timing.parse_grouped_output(
            grouped_stdout("base", task), "base", task, 4096, 8)
        before = target_irq_snapshot(start_ns=100, end_ns=110)
        after = target_irq_snapshot(
            hard={"142": (1, 0, 0)}, start_ns=120, end_ns=130)
        delta = timing.checked_target_irq_delta(
            before, after, (8, 72, 126), "fixture")
        with self.assertRaisesRegex(timing.TimingError, "campaign-fatal"):
            timing._attempt_contaminations(
                parsed, 0, 0, 0, 1_000_000, delta)

    def test_process_identity_validator_binds_cleanup_and_frozen_inode(self):
        with tempfile.TemporaryDirectory() as temporary:
            binary = Path(temporary) / "binary"
            binary.write_bytes(b"fixture")
            command = ["/usr/bin/env", str(binary), "arg"]
            identity = {
                "pid": 123, "start_ticks": 456,
                "executable": {
                    "path": str(binary), "device": binary.stat().st_dev,
                    "inode": binary.stat().st_ino,
                },
                "argv": [str(binary), "arg"],
                "boot_id": "1788608a-7aa1-48de-8f7c-848485be3cc3",
                "binding_observation": "double-proc-snapshot",
            }
            timing._validate_process_identity_receipt(
                identity, "exited_group_swept", command, binary, "fixture")
            with self.assertRaisesRegex(timing.TimingError, "process-identity"):
                timing._validate_process_identity_receipt(
                    identity, "terminated_group", command, binary, "fixture")

    def test_execution_receipt_schema_matches_constructor(self):
        task = one_task()
        parsed = timing.parse_grouped_output(
            grouped_stdout("base", task), "base", task, 4096, 8)
        start_ns = 10
        end_ns = start_ns + timing.SIBLING_ACCEPTED_EXECUTION_MIN_DURATION_NS
        receipt = timing._execution_receipt(
            task, 0, "base", ["example"], 0,
            "2026-07-18T00:00:00.000Z", start_ns, end_ns, parsed, b"", [],
            "c" * 64, {
                "pid": 123, "start_ticks": 456,
                "executable": {
                    "path": "/example", "device": 7, "inode": 8,
                },
                "argv": ["/example"], "boot_id":
                    "1788608a-7aa1-48de-8f7c-848485be3cc3",
                "binding_observation": "double-proc-snapshot",
            }, "exited_group_swept",
            (0, 0, 0, 10, 0), (0, 0, 0, 11, 0),
            {"runtime_ns": 100, "pcount": 5},
            {"runtime_ns": 100, "pcount": 5}, *zero_irq_arguments())
        self.assertEqual(set(receipt), timing.EXECUTION_RECEIPT_FIELDS)
        self.assertEqual(receipt["timed_phase_ns"], parsed.timed_phase_ns)
        self.assertEqual(receipt["all_phase_ns"], parsed.all_phase_ns)
        self.assertEqual(receipt["sibling_busy_ticks"], 0)
        self.assertEqual(receipt["sibling_sched_runtime_ns"], 0)
        self.assertEqual(receipt["sibling_sched_pcount"], 0)
        self.assertEqual(
            receipt["sibling_sched_runtime_limit_ns"],
            timing.sibling_attempt_runtime_limit_ns(end_ns - start_ns))
        self.assertEqual(receipt["sibling_sched_pcount_limit"], 1)

    def test_execution_receipt_rejects_any_sibling_busy_tick(self):
        task = one_task()
        parsed = timing.parse_grouped_output(
            grouped_stdout("base", task), "base", task, 4096, 8)
        end_ns = 10 + timing.SIBLING_ACCEPTED_EXECUTION_MIN_DURATION_NS
        with self.assertRaisesRegex(
                timing.TimingError, "accepted execution used"):
            timing._execution_receipt(
                task, 0, "base", ["example"], 0,
                "2026-07-18T00:00:00.000Z", 10, end_ns, parsed, b"", [],
                "c" * 64, {
                    "pid": 123, "start_ticks": 456,
                    "executable": {
                        "path": "/example", "device": 7, "inode": 8,
                    },
                    "argv": ["/example"], "boot_id":
                        "1788608a-7aa1-48de-8f7c-848485be3cc3",
                    "binding_observation": "double-proc-snapshot",
                }, "exited_group_swept",
                (0, 0, 0, 10, 0), (1, 0, 0, 11, 0),
                {"runtime_ns": 100, "pcount": 5},
                {"runtime_ns": 100, "pcount": 5}, *zero_irq_arguments())

    def test_execution_receipt_uses_exact_schedstat_limits_and_rejects_short_run(self):
        task = one_task()
        parsed = timing.parse_grouped_output(
            grouped_stdout("base", task), "base", task, 4096, 8)
        duration_ns = timing.SIBLING_ACCEPTED_EXECUTION_MIN_DURATION_NS
        runtime_limit = timing.sibling_attempt_runtime_limit_ns(duration_ns)
        arguments = (
            task, 0, "base", ["example"], 0,
            "2026-07-18T00:00:00.000Z", 10,
            10 + duration_ns,
            parsed, b"", [], "c" * 64, {}, "exited_group_swept",
            (0, 0, 0, 10, 0), (0, 0, 0, 11, 0),
            {"runtime_ns": 100, "pcount": 5},
        )
        accepted = timing._execution_receipt(
            *arguments,
            {"runtime_ns": 100 + runtime_limit, "pcount": 6},
            *zero_irq_arguments())
        self.assertEqual(
            accepted["sibling_sched_runtime_limit_ns"], runtime_limit)
        self.assertEqual(accepted["sibling_sched_pcount_limit"], 1)
        with self.assertRaisesRegex(timing.TimingError, "used the SMT sibling"):
            timing._execution_receipt(
                *arguments,
                {"runtime_ns": 101 + runtime_limit, "pcount": 6},
                *zero_irq_arguments())
        with self.assertRaisesRegex(timing.TimingError, "used the SMT sibling"):
            timing._execution_receipt(
                *arguments,
                {"runtime_ns": 100 + runtime_limit, "pcount": 7},
                *zero_irq_arguments())
        with self.assertRaisesRegex(timing.TimingError, "shorter"):
            timing._execution_receipt(
                task, 0, "base", ["example"], 0,
                "2026-07-18T00:00:00.000Z", 10,
                9 + timing.SIBLING_ACCEPTED_EXECUTION_MIN_DURATION_NS,
                parsed, b"", [], "c" * 64, {}, "exited_group_swept",
                (0, 0, 0, 10, 0), (0, 0, 0, 11, 0),
                {"runtime_ns": 100, "pcount": 5},
                {"runtime_ns": 100, "pcount": 5}, *zero_irq_arguments())

    def test_execution_receipt_replay_rejects_resealed_schedstat_tamper(self):
        task = one_task()
        raw = grouped_stdout("base", task)
        parsed = timing.parse_grouped_output(raw, "base", task, 4096, 8)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for directory in ("frozen", "raw", "stderr", "exit",
                              "receipts", "contamination"):
                (root / directory).mkdir(mode=0o700)
            binary = root / "frozen" / timing.BINARY_NAMES["base"]
            binary.write_bytes(b"fixture")
            design = {
                "root": str(root), "core": 8, "numa_node": 0,
                "evict_bytes": 4096, "controller_core": 126,
                "topology": {"sibling": 72},
                "tools": {
                    name: {"path": "/usr/bin/" + name}
                    for name in ("env", "taskset", "numactl")
                },
                "immutable_files": {
                    "frozen/" + timing.BINARY_NAMES["base"]:
                        timing.sha256_file(binary),
                },
            }
            command = timing.command_for(design, task, "base")
            name = timing.execution_name(task, 0, "base")
            contamination_start = 10
            contamination_end = contamination_start + \
                timing.SIBLING_ACCEPTED_EXECUTION_MIN_DURATION_NS
            contamination_digest = timing._save_contamination(
                root, name, 0, raw, b"", parsed, command,
                {
                    "pid": 122, "start_ticks": 455,
                    "executable": {
                        "path": str(binary), "device": binary.stat().st_dev,
                        "inode": binary.stat().st_ino,
                    },
                    "argv": command[command.index(str(binary)):],
                    "boot_id": "1788608a-7aa1-48de-8f7c-848485be3cc3",
                    "binding_observation": "double-proc-snapshot",
                }, "exited_group_swept",
                contamination_start, contamination_end,
                (0, 0, 0, 10, 0), (3, 0, 0, 11, 0),
                {"runtime_ns": 100, "pcount": 5},
                {"runtime_ns": 140, "pcount": 6},
                *zero_irq_arguments(contamination_start, contamination_end))
            binary_index = command.index(str(binary))
            binary_stat = binary.stat()
            start_ns = contamination_end + 100
            end_ns = start_ns + \
                timing.SIBLING_ACCEPTED_EXECUTION_MIN_DURATION_NS
            receipt = timing._execution_receipt(
                task, 0, "base", command, 1,
                "2026-07-18T00:00:00.000Z", start_ns, end_ns,
                parsed, b"", [{
                    "attempt": 0,
                    "receipt_sha256": contamination_digest,
                }], timing.sha256_file(binary), {
                    "pid": 123, "start_ticks": 456,
                    "executable": {
                        "path": str(binary), "device": binary_stat.st_dev,
                        "inode": binary_stat.st_ino,
                    },
                    "argv": command[binary_index:],
                    "boot_id": "1788608a-7aa1-48de-8f7c-848485be3cc3",
                    "binding_observation": "double-proc-snapshot",
                }, "exited_group_swept",
                (3, 0, 0, 20, 0), (3, 0, 0, 21, 0),
                {"runtime_ns": 200, "pcount": 10},
                {"runtime_ns": 200, "pcount": 10},
                *zero_irq_arguments(start_ns, end_ns))
            (root / "raw" / name).write_bytes(raw)
            (root / "stderr" / (name + ".stderr")).write_bytes(b"")
            (root / "exit" / (name + ".exit")).write_bytes(b"0\n")
            receipt_path = root / "receipts" / (name + ".json")
            receipt_path.write_bytes(timing.canonical_json(receipt))
            validated = timing._validate_execution_receipt(
                root, design, task, 0, 0, end_ns)
            self.assertEqual(validated[2:5], (3, 40, 1))

            payload = {
                key: value for key, value in receipt.items()
                if key not in ("schema", "self_sha256_excluding_field")
            }
            payload["sibling_schedstat_after"] = {
                "runtime_ns": 201, "pcount": 10}
            tampered = timing.sealed_record(receipt["schema"], payload)
            receipt_path.write_bytes(timing.canonical_json(tampered))
            with self.assertRaisesRegex(timing.TimingError, "counter receipt"):
                timing._validate_execution_receipt(
                    root, design, task, 0, 0, end_ns)

    def test_sibling_tick_delta_and_50ppm_campaign_boundary(self):
        before = (10, 20, 30, 100, 5, 6, 7, 0, 0, 0)
        after = (11, 20, 31, 200, 6, 6, 8, 0, 0, 0)
        self.assertEqual(
            timing.checked_busy_tick_delta(before, after, "fixture"), 3)
        with self.assertRaisesRegex(timing.TimingError, "malformed"):
            timing.checked_busy_tick_delta(before, after[:-1], "fixture")
        with self.assertRaisesRegex(timing.TimingError, "malformed"):
            timing.checked_busy_tick_delta(before, (*after[:-1], True),
                                           "fixture")
        with mock.patch.object(timing, "CLOCK_TICKS_PER_SECOND", 100):
            self.assertEqual(
                timing.sibling_campaign_busy_limit_ns(180_000_000_000), 1)
            self.assertEqual(
                timing.sibling_campaign_busy_limit_ns(200_000_000_000), 1)
            self.assertEqual(
                timing.sibling_campaign_busy_limit_ns(200_000_001_000), 1)
            self.assertEqual(
                timing.sibling_campaign_busy_limit_ns(399_999_999_999), 1)
            self.assertEqual(
                timing.sibling_campaign_busy_limit_ns(400_000_000_000), 2)
            self.assertEqual(
                timing.sibling_campaign_busy_limit_ns(1_500_000_000_000), 7)
        self.assertEqual(timing.sibling_campaign_runtime_limit_ns(19_999), 0)
        self.assertEqual(timing.sibling_campaign_runtime_limit_ns(20_000), 1)
        self.assertEqual(
            timing.sibling_campaign_runtime_limit_ns(200_000_000_000),
            10_000_000)
        self.assertEqual(timing.sibling_attempt_runtime_limit_ns(19_999), 0)
        self.assertEqual(timing.sibling_attempt_runtime_limit_ns(20_000), 1)
        self.assertEqual(timing.sibling_attempt_pcount_limit(1), 1)
        self.assertEqual(
            timing.sibling_attempt_pcount_limit(1_000_000_000), 1)
        self.assertEqual(
            timing.sibling_attempt_pcount_limit(1_000_000_001), 2)
        self.assertEqual(
            timing.sibling_attempt_pcount_limit(2_000_000_000), 2)
        for malformed in (0, -1, True, 1.5):
            with self.subTest(malformed=malformed):
                with self.assertRaisesRegex(timing.TimingError, "duration"):
                    timing.sibling_attempt_runtime_limit_ns(malformed)
                with self.assertRaisesRegex(timing.TimingError, "duration"):
                    timing.sibling_attempt_pcount_limit(malformed)

    def test_schedstat_v15_parser_is_strict_and_accepts_live_fields(self):
        raw = schedstat_text()
        self.assertEqual(timing.parse_schedstat_cpu(raw, 1), {
            "runtime_ns": 2000, "pcount": 20})
        malformed = (
            raw.replace(b"version 15", b"version 14", 1),
            raw.replace(b"cpu1 ", b"cpu0 ", 1),
            raw.replace(b"cpu0 7 0 ", b"cpu0 7 1 ", 1),
            raw.replace(b"domain0 ", b"domain1 ", 1),
            raw.replace(b" 1000 12 10\n", b" -1 12 10\n", 1),
            raw.replace(b" 1000 12 10\n", b" " + b"9" * 21 +
                        b" 12 10\n", 1),
            raw.replace(b"cpu0 7 0 8 9 10 11 1000 12 10",
                        b"cpu0 7 0 8 9 10 11 1000 12"),
        )
        for value in malformed:
            with self.subTest(value=value[:40]):
                with self.assertRaises(timing.TimingError):
                    timing.parse_schedstat_cpu(value, 0)
        with self.assertRaisesRegex(timing.TimingError, "malformed"):
            timing.parse_schedstat_cpu(raw, True)
        with self.assertRaisesRegex(timing.TimingError, "missing"):
            timing.parse_schedstat_cpu(raw, 2)

    def test_schedstat_delta_order_and_typed_tamper_checks(self):
        before = {"runtime_ns": 100, "pcount": 5}
        after = {"runtime_ns": 140, "pcount": 6}
        self.assertEqual(
            timing.checked_schedstat_delta(before, after, "fixture"),
            (40, 1))
        with self.assertRaisesRegex(timing.TimingError, "malformed"):
            timing.checked_schedstat_delta(
                before, {"runtime_ns": 99, "pcount": 6}, "fixture")
        with self.assertRaisesRegex(timing.TimingError, "malformed"):
            timing.checked_schedstat_delta(
                before, {"runtime_ns": 140, "pcount": True}, "fixture")
        with self.assertRaisesRegex(timing.TimingError, "counter receipt"):
            timing.require_exact_counter(True, 1, "tampered pcount")

    def test_target_irq_parser_delta_classification_reset_and_topology(self):
        before = target_irq_snapshot()
        after = target_irq_snapshot(
            hard={"LOC": (2, 0, 0)},
            soft={"TIMER": (0, 1, 0), "RCU": (0, 0, 3)},
            start_ns=120, end_ns=130)
        delta = timing.checked_target_irq_delta(
            before, after, (8, 72, 126), "fixture")
        self.assertEqual(delta["contaminations"], [])
        self.assertEqual(delta["classifications"], [
            "named-hardirq:LOC:cpu8:2",
            "softirq:TIMER:cpu72:1",
            "softirq:RCU:cpu126:3",
        ])
        self.assertEqual(
            timing.validate_target_irq_delta(
                delta, before, after, (8, 72, 126), "fixture"), delta)

        numeric = timing.checked_target_irq_delta(
            before, target_irq_snapshot(
                hard={"142": (1, 0, 0)}, start_ns=120, end_ns=130),
            (8, 72, 126), "fixture")
        self.assertEqual(
            numeric["contaminations"], ["numeric-hardirq:142:cpu8:1"])
        device_soft = timing.checked_target_irq_delta(
            before,
            target_irq_snapshot(
                soft={"NET_RX": (0, 1, 1)}, start_ns=120, end_ns=130),
            (8, 72, 126), "fixture")
        self.assertEqual(
            device_soft["contaminations"], ["softirq:NET_RX:cpu72:1"])
        self.assertIn(
            "softirq:NET_RX:cpu126:1", device_soft["classifications"])
        global_hard = timing.checked_target_irq_delta(
            before, target_irq_snapshot(
                global_hard={"ERR": 1}, start_ns=120, end_ns=130),
            (8, 72, 126), "fixture")
        self.assertEqual(
            global_hard["contaminations"], ["global-hardirq:ERR:1"])
        self.assertIn(
            ["ERR", "global", 1], global_hard["hardirq_deltas"])

        with self.assertRaisesRegex(timing.TimingError, "counter reset"):
            timing.checked_target_irq_delta(
                target_irq_snapshot(
                    hard={"142": (1, 0, 0)}, start_ns=80, end_ns=90),
                before,
                (8, 72, 126), "fixture")
        with self.assertRaisesRegex(timing.TimingError, "counter reset"):
            timing.checked_target_irq_delta(
                target_irq_snapshot(
                    global_hard={"MIS": 1}, start_ns=80, end_ns=90), before,
                (8, 72, 126), "fixture")
        changed_topology = copy.deepcopy(after)
        changed_topology["hardirq_rows"].pop()
        changed_topology["hardirq_sha256"] = timing._target_rows_sha256(
            "hardirq", (8, 72, 126), changed_topology["cpu_ids"],
            changed_topology["hardirq_rows"])
        with self.assertRaises(timing.TimingError):
            timing.checked_target_irq_delta(
                before, changed_topology, (8, 72, 126), "fixture")
        malformed_header = (
            b" CPU8 CPU72\n 30: 0 0 AMD-Vi0-PPR\n")
        with self.assertRaisesRegex(timing.TimingError, "topology changed"):
            timing.parse_target_irq_snapshot(
                malformed_header,
                b" CPU8 CPU72\n TIMER: 0 0\n", (8, 72, 126), 100, 110)
        tampered = copy.deepcopy(before)
        tampered["hardirq_rows"][0][2] = True
        with self.assertRaisesRegex(timing.TimingError, "hardirq row"):
            timing.validate_target_irq_snapshot(
                tampered, (8, 72, 126), "tampered")
        non_ascii_identity = copy.deepcopy(before)
        non_ascii_identity["hardirq_rows"][0][2] += "\N{SNOWMAN}"
        non_ascii_identity["hardirq_sha256"] = timing._target_rows_sha256(
            "hardirq", (8, 72, 126), non_ascii_identity["cpu_ids"],
            non_ascii_identity["hardirq_rows"])
        with self.assertRaisesRegex(timing.TimingError, "hardirq row"):
            timing.validate_target_irq_snapshot(
                non_ascii_identity, (8, 72, 126), "tampered")
        for global_tamper in (
                ["ERR", "global", 0, 0],
                ["ERR", "named", "Error counters", 0, 0, 0],
                ["142", "global", 0]):
            tampered = copy.deepcopy(before)
            index = next(
                index for index, row in enumerate(tampered["hardirq_rows"])
                if row[0] == "ERR")
            tampered["hardirq_rows"][index] = global_tamper
            tampered["hardirq_sha256"] = timing._target_rows_sha256(
                "hardirq", (8, 72, 126), tampered["cpu_ids"],
                tampered["hardirq_rows"])
            with self.subTest(global_tamper=global_tamper):
                with self.assertRaises(timing.TimingError):
                    timing.validate_target_irq_snapshot(
                        tampered, (8, 72, 126), "tampered")
        malformed_global = (
            b" CPU8 CPU72 CPU126\n ERR: 0 0\n")
        with self.assertRaisesRegex(timing.TimingError, "width changed"):
            timing.parse_target_irq_snapshot(
                malformed_global,
                b" CPU8 CPU72 CPU126\n TIMER: 0 0 0\n",
                (8, 72, 126), 100, 110)
        full_width_global = b" CPU8 CPU72 CPU126\n ERR: 0 0 0\n"
        full_softirq = (
            " CPU8 CPU72 CPU126\n" + "".join(
                " %s: 0 0 0\n" % vector
                for vector in timing.EXPECTED_SOFTIRQ_VECTORS)
        ).encode("ascii")
        with self.assertRaisesRegex(timing.TimingError, "global hardirq"):
            timing.parse_target_irq_snapshot(
                full_width_global, full_softirq,
                (8, 72, 126), 100, 110)

        missing_global = copy.deepcopy(before)
        missing_global["hardirq_rows"] = [
            row for row in missing_global["hardirq_rows"] if row[0] != "MIS"]
        missing_global["hardirq_sha256"] = timing._target_rows_sha256(
            "hardirq", (8, 72, 126), missing_global["cpu_ids"],
            missing_global["hardirq_rows"])
        with self.assertRaisesRegex(timing.TimingError, "global hardirq"):
            timing.validate_target_irq_snapshot(
                missing_global, (8, 72, 126), "tampered")

        high_priority = timing.checked_target_irq_delta(
            before, target_irq_snapshot(
                soft={"HI": (1, 1, 0)}, start_ns=120, end_ns=130),
            (8, 72, 126), "fixture")
        self.assertEqual(high_priority["contaminations"], [
            "softirq:HI:cpu8:1", "softirq:HI:cpu72:1"])

        changed_header = target_irq_snapshot(
            start_ns=120, end_ns=130, cpu_ids=(1, 8, 72, 126))
        with self.assertRaisesRegex(timing.TimingError, "topology changed"):
            timing.checked_target_irq_delta(
                before, changed_header, (8, 72, 126), "fixture")

        changed_identity = copy.deepcopy(after)
        numeric_index = next(
            index for index, row in enumerate(changed_identity["hardirq_rows"])
            if row[0] == "142")
        changed_identity["hardirq_rows"][numeric_index][2] = "foreign-device"
        changed_identity["hardirq_sha256"] = timing._target_rows_sha256(
            "hardirq", (8, 72, 126), changed_identity["cpu_ids"],
            changed_identity["hardirq_rows"])
        with self.assertRaisesRegex(timing.TimingError, "topology changed"):
            timing.checked_target_irq_delta(
                before, changed_identity, (8, 72, 126), "fixture")

        missing_softirq = copy.deepcopy(before)
        missing_softirq["softirq_rows"].pop(0)
        missing_softirq["softirq_sha256"] = timing._target_rows_sha256(
            "softirq", (8, 72, 126), missing_softirq["cpu_ids"],
            missing_softirq["softirq_rows"])
        with self.assertRaisesRegex(timing.TimingError, "softirq topology"):
            timing.validate_target_irq_snapshot(
                missing_softirq, (8, 72, 126), "tampered")

    def test_sibling_policy_and_launch_interval_reject_bool_or_nonfinite_aliases(self):
        timing.require_frozen_sibling_idle_policy(timing.SIBLING_IDLE_POLICY)
        changed = {
            **timing.SIBLING_IDLE_POLICY,
            "campaign_min_busy_ticks": True,
        }
        with self.assertRaisesRegex(timing.TimingError, "sibling_idle_policy"):
            timing.require_frozen_sibling_idle_policy(changed)
        changed = {
            **timing.SIBLING_IDLE_POLICY,
            "schedstat_sched_schedstats_sysctl_required": 0,
        }
        with self.assertRaisesRegex(timing.TimingError, "sibling_idle_policy"):
            timing.require_frozen_sibling_idle_policy(changed)
        valid = {
            "start_monotonic_s": 10.0, "end_monotonic_s": 11.0,
            "duration_s": 1.0, "start_monotonic_ns": 10_000_000_000,
            "end_monotonic_ns": 11_000_000_000,
            "duration_ns": 1_000_000_000,
        }
        self.assertEqual(timing.checked_campaign_interval(valid)[-1],
                         1_000_000_000)
        for field, value in (("duration_s", True),
                             ("duration_s", float("nan")),
                             ("start_monotonic_s", float("inf")),
                             ("start_monotonic_s", 10 ** 400),
                             ("duration_ns", True)):
            with self.subTest(field=field, value=value):
                with self.assertRaisesRegex(timing.TimingError, "interval"):
                    timing.checked_campaign_interval({**valid, field: value})

    def test_five_second_preflight_limits_are_exact_and_tamper_evident(self):
        duration_ns = timing.SIBLING_PREFLIGHT_WINDOW_NS
        launch = {
            "preflight_start_monotonic_ns": 100,
            "preflight_end_monotonic_ns": 100 + duration_ns,
            "preflight_duration_ns": duration_ns,
            "preflight_sibling_sched_runtime_limit_ns":
                timing.sibling_attempt_runtime_limit_ns(duration_ns),
            "preflight_sibling_sched_pcount_limit":
                timing.sibling_attempt_pcount_limit(duration_ns),
        }
        self.assertEqual(
            timing.checked_preflight_limits(
                launch, launch["preflight_end_monotonic_ns"]),
            (duration_ns, 250_000, 5))
        for field, value in (
                ("preflight_duration_ns", duration_ns - 1),
                ("preflight_sibling_sched_runtime_limit_ns", 250_001),
                ("preflight_sibling_sched_pcount_limit", True)):
            with self.subTest(field=field):
                with self.assertRaises(timing.TimingError):
                    timing.checked_preflight_limits(
                        {**launch, field: value},
                        launch["preflight_end_monotonic_ns"])
        with self.assertRaisesRegex(timing.TimingError, "interval"):
            timing.checked_preflight_limits(
                launch, launch["preflight_end_monotonic_ns"] - 1)

    def test_runtime_isolation_snapshot_transition_and_tamper_checks(self):
        start = isolation_snapshot(100, 120)
        end = isolation_snapshot(10_000, 10_030)
        timing.validate_runtime_isolation_transition(
            start, end, (6, 70, 80))
        digest = timing.runtime_isolation_snapshot_sha256(start)
        self.assertRegex(digest, r"^[0-9a-f]{64}$")
        changed = copy.deepcopy(start)
        changed["capture_end_monotonic_ns"] += 1
        changed["capture_duration_ns"] += 1
        self.assertNotEqual(
            digest, timing.runtime_isolation_snapshot_sha256(changed))

        malformed_cases = []
        malformed = copy.deepcopy(start)
        malformed["self_cgroup"] = "/"
        malformed_cases.append(malformed)
        malformed = copy.deepcopy(start)
        del malformed["cgroup_partition"]
        malformed_cases.append(malformed)
        malformed = copy.deepcopy(start)
        malformed["irq_effective_affinities"][1] = {
            "irq": 31, "effective_affinity_list": "70",
            "effective_cpus": [70],
        }
        malformed_cases.append(malformed)
        malformed = copy.deepcopy(start)
        malformed["zero_count_irq_exception"][
            "global_interrupt_count"] = True
        malformed_cases.append(malformed)
        malformed = copy.deepcopy(start)
        malformed["zero_count_irq_exception"]["identity"] = \
            "OTHER 999-edge AMDI0010:05"
        malformed_cases.append(malformed)
        malformed = copy.deepcopy(start)
        malformed["zero_count_irq_exception"][
            "requested_affinity_list"] = "0-126"
        malformed["zero_count_irq_exception"][
            "requested_cpus"] = list(range(127))
        malformed_cases.append(malformed)
        malformed = copy.deepcopy(start)
        malformed["zero_count_irq_exception"][
            "effective_affinity_list"] = "70"
        malformed["zero_count_irq_exception"]["effective_cpus"] = [70]
        malformed_cases.append(malformed)
        malformed = copy.deepcopy(start)
        malformed["zero_count_irq_exception"][
            "handler_directories"] = ["AMDI0010:04"]
        malformed_cases.append(malformed)
        malformed = copy.deepcopy(start)
        malformed["managed_nvme_exceptions"][0][
            "handler_directories"] = ["nvme2q10"]
        malformed_cases.append(malformed)
        malformed = copy.deepcopy(start)
        malformed["managed_nvme_exceptions"][0]["identity"] = \
            "IR-PCI-MSIX-ffff:ff:ff.f 9-edge nvme2q9"
        malformed_cases.append(malformed)
        malformed = copy.deepcopy(start)
        malformed["managed_nvme_exceptions"][0][
            "requested_affinity_list"] = "8-9"
        malformed["managed_nvme_exceptions"][0]["requested_cpus"] = [8, 9]
        malformed_cases.append(malformed)
        malformed = copy.deepcopy(start)
        malformed["managed_nvme_exceptions"][0][
            "effective_affinity_list"] = "72"
        malformed["managed_nvme_exceptions"][0]["effective_cpus"] = [72]
        malformed_cases.append(malformed)
        malformed = copy.deepcopy(start)
        malformed["managed_nvme_exceptions"].reverse()
        malformed_cases.append(malformed)
        malformed = copy.deepcopy(start)
        malformed["capture_duration_ns"] += 1
        malformed_cases.append(malformed)
        for value in malformed_cases:
            with self.subTest(value=value):
                with self.assertRaises(timing.TimingError):
                    timing.validate_runtime_isolation_snapshot(
                        value, (6, 70, 80), "tampered")
        with self.assertRaisesRegex(timing.TimingError, "cardinality"):
            timing.parse_cpu_list("0-2147483647")

    def test_runtime_isolation_capture_audits_cgroup_and_all_numeric_irqs(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            proc = root / "proc"
            cgroup = root / "cgroup"
            (proc / "self").mkdir(parents=True)
            (proc / "self/cgroup").write_text(
                "0::/wh2-timing-v4\n", encoding="ascii")
            zero_irq = proc / ("irq/%d" % timing.ZERO_COUNT_GUARDED_IRQ)
            (zero_irq / timing.ZERO_COUNT_GUARDED_IRQ_HANDLER).mkdir(
                parents=True)
            (zero_irq / "effective_affinity_list").write_text(
                "6\n", encoding="ascii")
            (zero_irq / "smp_affinity_list").write_text(
                "0-127\n", encoding="ascii")
            old_zero_irq = proc / "irq/30"
            (old_zero_irq / "AMD-Vi0-PPR").mkdir(parents=True)
            (old_zero_irq / "effective_affinity_list").write_text(
                "8\n", encoding="ascii")
            for irq, _identity, handler, requested, effective in \
                    timing.MANAGED_NVME_IRQ_WHITELIST:
                path = proc / ("irq/%d" % irq)
                (path / handler).mkdir(parents=True)
                (path / "effective_affinity_list").write_text(
                    effective + "\n", encoding="ascii")
                (path / "smp_affinity_list").write_text(
                    requested + "\n", encoding="ascii")
            irq2 = proc / "irq/2"
            irq2.mkdir(parents=True)
            (irq2 / "effective_affinity_list").write_text(
                "\n", encoding="ascii")
            irq31 = proc / "irq/31"
            irq31.mkdir(parents=True)
            (irq31 / "effective_affinity_list").write_text(
                "0-5\n", encoding="ascii")
            interrupt_lines = [" CPU0 CPU1"] + [
                " %d: 0 0 %s" % (irq, identity)
                for irq, identity in sorted(
                    timing.GUARDED_IRQ_IDENTITIES.items())
            ]
            (proc / "interrupts").write_text(
                "\n".join(interrupt_lines) + "\n", encoding="ascii")
            group = cgroup / "wh2-timing-v4"
            group.mkdir(parents=True)
            (cgroup / "cpuset.cpus.isolated").write_text(
                "6,70,80\n", encoding="ascii")
            for name in (
                    "cpuset.cpus", "cpuset.cpus.effective",
                    "cpuset.cpus.exclusive",
                    "cpuset.cpus.exclusive.effective"):
                (group / name).write_text("6,70,80\n", encoding="ascii")
            (group / "cpuset.cpus.partition").write_text(
                "isolated\n", encoding="ascii")
            with mock.patch.object(
                    timing.time, "monotonic_ns", side_effect=[100, 200]):
                captured = timing.capture_runtime_isolation_snapshot(
                    6, 70, 80, proc_root=proc, cgroup_root=cgroup)
            self.assertEqual(captured["capture_duration_ns"], 100)
            self.assertEqual(
                captured["zero_count_irq_exception"][
                    "global_interrupt_count"], 0)
            self.assertEqual(
                [item["irq"] for item in
                 captured["irq_effective_affinities"]],
                sorted([2, timing.ZERO_COUNT_GUARDED_IRQ, 30, 31] + [
                    item[0] for item in timing.MANAGED_NVME_IRQ_WHITELIST]))
            irq30_record = next(
                item for item in captured["irq_effective_affinities"]
                if item["irq"] == 30)
            self.assertEqual(irq30_record["effective_cpus"], [8])
            self.assertEqual(
                captured["irq_effective_affinities"][0]["effective_cpus"], [])

            (irq31 / "effective_affinity_list").write_text(
                "70\n", encoding="ascii")
            with self.assertRaisesRegex(timing.TimingError, "IRQ 31"):
                timing.capture_runtime_isolation_snapshot(
                    6, 70, 80, proc_root=proc, cgroup_root=cgroup)
            (irq31 / "effective_affinity_list").write_text(
                "0-5\n", encoding="ascii")
            (proc / "interrupts").write_text(
                " CPU0 CPU1\n"
                " 23: 0 1 IR-IO-APIC 23-edge AMDI0010:05\n",
                encoding="ascii")
            with self.assertRaises(timing.TimingError):
                timing.capture_runtime_isolation_snapshot(
                    6, 70, 80, proc_root=proc, cgroup_root=cgroup)

    def test_guarded_irq_parser_rejects_nonzero_identity_or_malformed_rows(self):
        zero = (
            b" CPU0 CPU1\n"
            b" 23: 0 0 IR-IO-APIC 23-edge AMDI0010:05\n")
        expected = {
            timing.ZERO_COUNT_GUARDED_IRQ:
                timing.ZERO_COUNT_GUARDED_IRQ_IDENTITY}
        self.assertEqual(
            timing._parse_guarded_irq_rows(
                zero, expected)[timing.ZERO_COUNT_GUARDED_IRQ]["total_count"],
            0)
        self.assertEqual(timing._parse_guarded_irq_rows(
            zero.replace(b"23: 0 0", b"23: 0 1"), expected
        )[timing.ZERO_COUNT_GUARDED_IRQ]["total_count"], 1)
        for malformed in (
                zero[:-1],
                zero.replace(b"CPU1", b"CPU0"),
                zero + zero.splitlines(keepends=True)[1],
                zero.replace(b"AMDI0010:05", b"other")):
            with self.subTest(malformed=malformed):
                with self.assertRaises(timing.TimingError):
                    timing._parse_guarded_irq_rows(malformed, expected)

    def test_sibling_only_contamination_is_receipted_for_retry(self):
        task = one_task()
        parsed = timing.parse_grouped_output(
            grouped_stdout("base", task), "base", task, 4096, 8)
        self.assertEqual(parsed.contaminations, ())
        duration_ns = timing.SIBLING_ACCEPTED_EXECUTION_MIN_DURATION_NS
        runtime_limit = timing.sibling_attempt_runtime_limit_ns(duration_ns)
        irq_before, irq_after, target_cpus = zero_irq_arguments()
        irq_delta = timing.checked_target_irq_delta(
            irq_before, irq_after, target_cpus, "fixture")
        self.assertEqual(
            timing._attempt_contaminations(
                parsed, 0, runtime_limit, 1, duration_ns, irq_delta), ())
        self.assertEqual(
            timing._attempt_contaminations(
                parsed, 3, runtime_limit + 1, 2, duration_ns, irq_delta),
            ("sibling-busy:3",
             "sibling-sched-runtime-ns:%d" % (runtime_limit + 1),
             "sibling-sched-pcount:2"))
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "contamination").mkdir(mode=0o700)
            digest = timing._save_contamination(
                root, "sample", 0, grouped_stdout("base", task), b"",
                parsed, ["example"], {
                    "pid": 123, "start_ticks": 456,
                    "executable": {"path": "/example", "device": 7,
                                   "inode": 8},
                    "argv": ["/example"],
                    "boot_id": "1788608a-7aa1-48de-8f7c-848485be3cc3",
                    "binding_observation": "double-proc-snapshot",
                }, "exited_group_swept", 10, 10 + duration_ns,
                (0, 0, 0, 10, 0), (3, 0, 0, 11, 0),
                {"runtime_ns": 100, "pcount": 5},
                {"runtime_ns": 101 + runtime_limit, "pcount": 7},
                *zero_irq_arguments())
            receipt = timing.load_canonical(
                root / "contamination" / "sample.attempt0.json",
                "test contamination")
        self.assertRegex(digest, r"^[0-9a-f]{64}$")
        self.assertEqual(set(receipt), timing.CONTAMINATION_RECEIPT_FIELDS)
        self.assertEqual(receipt["schema"],
                         "wirehair.wh2.grouped_commit_timing."
                         "contamination_receipt.v4")
        self.assertEqual(receipt["contaminations"], [
            "sibling-busy:3",
            "sibling-sched-runtime-ns:%d" % (runtime_limit + 1),
            "sibling-sched-pcount:2"])
        self.assertEqual(receipt["sibling_busy_ticks"], 3)
        self.assertEqual(receipt["duration_ns"], duration_ns)
        self.assertEqual(
            receipt["sibling_sched_runtime_ns"], runtime_limit + 1)
        self.assertEqual(receipt["sibling_sched_runtime_limit_ns"], runtime_limit)
        self.assertEqual(receipt["sibling_sched_pcount"], 2)
        self.assertEqual(receipt["sibling_sched_pcount_limit"], 1)

    def test_failure_receipt_records_exact_schedstat_accounting(self):
        task = one_task()
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "failure").mkdir(mode=0o700)
            frozen = root / "frozen"
            frozen.mkdir(mode=0o700)
            binary = frozen / timing.BINARY_NAMES["base"]
            binary.write_bytes(b"fixture")
            binary.chmod(0o555)
            binary_sha256 = timing.sha256_file(binary)
            command = [str(binary), "example"]
            identity = {
                "pid": 123, "start_ticks": 456,
                "executable": {
                    "path": str(binary), "device": binary.stat().st_dev,
                    "inode": binary.stat().st_ino,
                },
                "argv": command,
                "boot_id": "1788608a-7aa1-48de-8f7c-848485be3cc3",
                "binding_observation": "double-proc-snapshot",
            }
            design = {
                "immutable_files": {
                    "frozen/" + timing.BINARY_NAMES["base"]: binary_sha256,
                },
                "core": 8, "topology": {"sibling": 72},
                "controller_core": 126,
            }
            with mock.patch.object(
                    timing, "_load_design", return_value=design), \
                    mock.patch.object(
                        timing, "_load_tasks", return_value=(task,)), \
                    mock.patch.object(
                        timing, "command_for", return_value=command):
                digest = timing._save_failure(
                    root, task, 0, "base", command, 0,
                    "2026-07-18T00:00:00.000Z", 10, 20,
                    b"out", b"err", 7, timing.TimingError("failure"),
                    binary_sha256, identity, "exited_group_swept", None,
                    (0, 0, 0, 10, 0), (2, 0, 0, 11, 0),
                    {"runtime_ns": 100, "pcount": 5},
                    {"runtime_ns": 125, "pcount": 6},
                    *zero_irq_arguments(10, 20))
                name = (
                    timing.execution_name(task, 0, "base") +
                    ".attempt0.json")
                receipt_path = root / "failure" / name
                receipt = timing.load_canonical(receipt_path, "test failure")
                with self.assertRaises(timing.TimingError):
                    timing.replay_failure_receipt(root, name, (1, 2, 3))

                for field, forged in (
                        ("cleanup_action", "impossible"),
                        ("process_identity", {"forged": True}),
                        ("argv", ["forged"]),
                        ("binary_sha256", "0" * 64)):
                    tampered_payload = {
                        key: value for key, value in receipt.items()
                        if key not in ("schema", "self_sha256_excluding_field")
                    }
                    tampered_payload[field] = forged
                    tampered = timing.sealed_record(
                        receipt["schema"], tampered_payload)
                    receipt_path.chmod(0o600)
                    receipt_path.write_bytes(timing.canonical_json(tampered))
                    receipt_path.chmod(0o444)
                    with self.subTest(forged_field=field):
                        with self.assertRaises(timing.TimingError):
                            timing.replay_failure_receipt(
                                root, name, (8, 72, 126))

                tampered_payload = {
                    key: value for key, value in receipt.items()
                    if key not in ("schema", "self_sha256_excluding_field")
                }
                tampered_payload["sibling_sched_pcount"] = True
                tampered = timing.sealed_record(
                    receipt["schema"], tampered_payload)
                receipt_path.chmod(0o600)
                receipt_path.write_bytes(timing.canonical_json(tampered))
                receipt_path.chmod(0o444)
                with self.assertRaisesRegex(
                        timing.TimingError, "counter receipt"):
                    timing.replay_failure_receipt(
                        root, name, (8, 72, 126))
                tampered_payload["sibling_sched_pcount"] = 1
                tampered_payload["sibling_sched_runtime_limit_ns"] = 1
                tampered = timing.sealed_record(
                    receipt["schema"], tampered_payload)
                receipt_path.chmod(0o600)
                receipt_path.write_bytes(timing.canonical_json(tampered))
                receipt_path.chmod(0o444)
                with self.assertRaisesRegex(timing.TimingError, "runtime limit"):
                    timing.replay_failure_receipt(
                        root, name, (8, 72, 126))
        self.assertRegex(digest, r"^[0-9a-f]{64}$")
        self.assertEqual(set(receipt), timing.FAILURE_RECEIPT_FIELDS)
        self.assertEqual(receipt["sibling_sched_runtime_ns"], 25)
        self.assertEqual(receipt["sibling_sched_runtime_limit_ns"], 0)
        self.assertEqual(receipt["sibling_sched_pcount"], 1)
        self.assertEqual(receipt["sibling_sched_pcount_limit"], 1)

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
            # The trailing builtin prevents bash from replacing itself with the
            # final sleep before the controller can bind its exact executable.
            "sleep 30 & child=$!; printf '%s\\n' \"$child\"; sleep 30; :",
        ]
        with self.assertRaisesRegex(
                timing.BoundCommandError, "exceeded its timeout") as raised:
            # Leave ample launch/binding slack under the intentionally saturated
            # host; this deadline is testing post-bind process-group cleanup.
            timing.run_bound_command(command, self.bash, self.python, 1.0)
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
            "be0bc94b97d03d6ddbc23db3b7058aa7f575b5cd")
        self.assertEqual(
            timing.MEASUREMENT_OVERLAY_COMMIT,
            "d6ab1a65c9864ad97ef06c4b88c2917cb387c0be")
        self.assertEqual(
            timing.MEASUREMENT_OVERLAY_PARENT_COMMIT,
            "3a659aeb132e6cf5e5e68d88094af8402cdb0e47")
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
