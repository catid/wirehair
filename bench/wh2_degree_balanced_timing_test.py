#!/usr/bin/env python3
"""Bounded parser, receipt, and tamper tests for degree-balanced timing."""

from __future__ import annotations

import csv
import io
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import wh2_degree_balanced_timing as timing


def one_task(**changes):
    value = {
        "K": 3200, "bb": 64, "seed_index": 0,
        "seed": timing.SCHEDULE_SEEDS[0][2], "schedule": "burst",
        "cache_state": "warm", "job": 0,
        "task_id": "000.K3200.bb64.seed0.burst.warm",
    }
    value.update(changes)
    return value


def grouped_stdout(
    orientation="forward", task=None, core=8, evict_bytes=4096,
    base_result=0, balanced_result=0, geometry="frozen",
):
    task = task or one_task()
    preamble = timing._expected_preamble(
        task, evict_bytes, orientation, geometry=geometry)
    if orientation == "forward":
        cli_semantics = {"control": "base", "candidate": "balanced"}
    elif orientation == "reverse":
        cli_semantics = {"control": "balanced", "candidate": "base"}
    elif orientation == "neutral":
        cli_semantics = {"control": "base", "candidate": "base"}
    else:
        raise AssertionError("unknown fixture orientation")
    results = {"base": base_result, "balanced": balanced_result}
    outcome = timing._outcome_class(
        results[cli_semantics["control"]], results[cli_semantics["candidate"]])
    common = results["base"] == results["balanced"] == 0
    staircase_rows = timing._timing_staircase_rows(int(task["K"]))
    preamble.update({
        "schema": "v5", "staircase_rows": str(staircase_rows),
        "dense_rows": "12", "heavy_rows": "12",
        "columns": str(int(task["K"]) + staircase_rows + 24),
        "control_attempt": "0", "control_matrix_seed": "0x123",
        "control_peel_seed": "0x456", "candidate_attempt": "0",
        "candidate_matrix_seed": "0x123", "candidate_peel_seed": "0x456",
        "payload_fingerprint": "0x5555555555555555",
        "packet_trace_sha256": "a" * 64, "cell_class": outcome,
        "common_success": "1" if common else "0", "trace_sha256": "a" * 64,
    })
    for arm, semantic in cli_semantics.items():
        result = results[semantic]
        preamble[arm + "_staircase_fingerprint"] = (
            "0x1111111111111111" if semantic == "base" else
            "0x2222222222222222")
        preamble[arm + "_dense_fingerprint"] = "0x3333333333333333"
        preamble[arm + "_heavy_fingerprint"] = "0x4444444444444444"
        preamble[arm + "_mixed_coefficient_fingerprint"] = (
            timing.FROZEN_MIXED_COEFFICIENT_FINGERPRINT
            if geometry == "frozen" else "0x7777777777777777")
        preamble["preflight_" + arm + "_result"] = str(result)
        preamble[arm + "_preflight_rhs_route"] = \
            "streamed" if result == 0 else "not-reached"
        preamble["preflight_" + arm + "_decoded_fingerprint"] = (
            "0x6666666666666666" if result == 0 else "none")
    missing = [key for key in timing.PREAMBLE_V5 if key not in preamble]
    if missing:
        raise AssertionError("fixture preamble missing %r" % missing)
    lines = [
        "# groupedtiming: " + " ".join(
            "%s=%s" % (key, preamble[key]) for key in timing.PREAMBLE_V5),
        ",".join(timing.TIMING_FIELDS),
    ]
    expected_arena = int(preamble["columns"]) * int(task["bb"])
    for cycle in range(4):
        for slot, marker in enumerate(timing.INNER_ORDER):
            arm = "control" if marker == "A" else "candidate"
            semantic = cli_semantics[arm]
            result = results[semantic]
            balanced = semantic == "balanced"
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
                "mixed_geometry": geometry,
                "degree_balanced_staircase": "1" if balanced else "0",
                "seed_attempt": "0", "matrix_seed": "0x123",
                "peel_seed": "0x456", "preflight_result": str(result),
                "cell_class": outcome, "common_success": "1" if common else "0",
                "result": str(result), "outcome_stable": "1",
                "elapsed_ns": str((900 if balanced else 1000) + cycle * 10 + slot),
                "saturated": "0", "cpu_before": str(core),
                "cpu_after": str(core), "cpu_migrated": "0",
                "minflt_delta": "0", "majflt_delta": "0",
                "fault_contaminated": "0",
                "inactivated": "107" if balanced else "109",
                "binary_def": "7" if balanced else "8", "heavy_gain": "8",
                "block_xors": "73000" if balanced else "74095",
                "block_muladds": "1880", "build_ns": "10",
                "peel_ns": "20", "project_ns": "30", "residual_ns": "40",
                "backsub_ns": "50", "joint_source_xors": "0",
                "joint_marginal_xors": "0", "joint_marginal_copies": "0",
                "joint_active_deltas": "0", "joint_scratch_bytes": "0",
                "dual_source_columns": "0",
                "source_bytes": str(int(task["K"]) * int(task["bb"])),
                "packet_payload_bytes": str(
                    (int(task["K"]) + timing.OVERHEAD) * int(task["bb"])),
                "intermediate_bytes": str(expected_arena),
                "solve_value_arena_bytes": str(expected_arena),
                "solve_value_eager_zero_bytes": "0",
                "solve_value_commit_copy_bytes": "0",
                "rhs_route_expected": "streamed",
                "rhs_route_actual":
                    "streamed" if result == 0 else "not-reached",
            }
            lines.append(",".join(values[field] for field in timing.TIMING_FIELDS))
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


def legacy_stdout_from_neutral_measurement(raw):
    lines = raw.decode("ascii").splitlines()
    pairs = dict(token.split("=", 1) for token in lines[0].split(" ")[2:])
    pairs["schema"] = "v2"
    output = [
        "# groupedtiming: " + " ".join(
            "%s=%s" % (key, pairs[key]) for key in timing.LEGACY_PREAMBLE_V2)
    ]
    header = lines[1].split(",")
    output.append(",".join(timing.LEGACY_FIELDS))
    for line in lines[2:]:
        row = dict(zip(header, line.split(",")))
        output.append(",".join(row[field] for field in timing.LEGACY_FIELDS))
    return ("\n".join(output) + "\n").encode("ascii")


def thermal_csv(times=(9.0, 10.0, 11.0, 12.0)):
    output = io.StringIO()
    writer = csv.DictWriter(
        output, fieldnames=timing.THERMAL_FIELDS, lineterminator="\r\n")
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


class TaskAndCommandTests(unittest.TestCase):
    def test_command_result_returns_complete_receipt(self):
        completed = mock.Mock(returncode=0, stdout=b"out", stderr=b"")
        with mock.patch.object(
                timing.subprocess, "run", return_value=completed) as run:
            receipt = timing.command_result(
                ("tool", "argument"), cwd=Path("/tmp"),
                environment={"LC_ALL": "C"})
        self.assertEqual(receipt["returncode"], 0)
        self.assertEqual(receipt["stdout"], b"out")
        self.assertEqual(receipt["stderr"], b"")
        self.assertEqual(receipt["argv"], ["tool", "argument"])
        self.assertEqual(receipt["cwd"], "/tmp")
        self.assertEqual(receipt["environment"], {"LC_ALL": "C"})
        self.assertGreaterEqual(receipt["duration_s"], 0.0)
        run.assert_called_once()

    def test_exact_grid_is_deterministic_and_complete(self):
        first = timing.generate_tasks()
        self.assertEqual(first, timing.generate_tasks())
        self.assertEqual(len(first), 252)
        coordinates = {
            (row["K"], row["bb"], row["schedule"], row["seed_index"],
             row["seed"], row["cache_state"]) for row in first
        }
        expected = {
            (K, width, schedule, seed_index, seed, cache)
            for K in timing.KS for width in timing.WIDTHS
            for schedule, seed_index, seed in timing.SCHEDULE_SEEDS
            for cache in timing.CACHE_STATES
        }
        self.assertEqual(coordinates, expected)
        self.assertEqual([row["job"] for row in first], list(range(252)))

    def test_prepare_rejects_nonfrozen_cpu_layout(self):
        with tempfile.TemporaryDirectory() as temporary:
            args = mock.Mock(
                result_dir=str(Path(temporary) / "result"), repo=temporary,
                core=9, controller_core=timing.CONTROLLER_CORE,
                numa_node=timing.NUMA_NODE, evict_bytes=256 * 1024 * 1024,
                build_jobs=1)
            with self.assertRaisesRegex(timing.TimingError, "frozen campaign"):
                timing.prepare_campaign(args)
        with self.assertRaisesRegex(timing.TimingError, "runtime CPU"):
            timing._validate_topology_again({
                "core": 9, "controller_core": timing.CONTROLLER_CORE,
                "sampler_core": timing.SAMPLER_CORE,
                "numa_node": timing.NUMA_NODE,
            })

    def test_production_source_hit_boundary_is_receipted(self):
        self.assertEqual(timing._certified_source_hits(9999), 2)
        self.assertEqual(timing._certified_source_hits(10000), 3)
        self.assertEqual(timing._certified_source_hits(64000), 3)
        with self.assertRaisesRegex(timing.TimingError, "codec domain"):
            timing._certified_source_hits(1)
        for K, expected in ((9999, "2"), (10000, "3")):
            task = one_task(K=K)
            preamble = timing._expected_preamble(task, 4096, "forward")
            self.assertEqual(preamble["source_hits"], expected)
            parsed = timing.parse_grouped_output(
                grouped_stdout(task=task), "forward", task, 4096, 8)
            self.assertEqual(parsed.preamble["source_hits"], expected)
        task = one_task(K=10000)
        tampered = mutate_preamble(
            grouped_stdout(task=task), "source_hits", "2")
        with self.assertRaisesRegex(timing.TimingError, "source_hits"):
            timing.parse_grouped_output(
                tampered, "forward", task, 4096, 8)

    def test_exact_timing_staircase_row_table_is_receipted(self):
        expected = {
            3200: 62, 9999: 86, 10000: 86, 20000: 134,
            32000: 190, 48466: 374, 64000: 346,
        }
        self.assertEqual(
            {K: timing._timing_staircase_rows(K) for K in timing.KS},
            expected)
        with self.assertRaisesRegex(timing.TimingError, "timing K grid"):
            timing._timing_staircase_rows(3199)
        task = one_task(K=32000)
        parsed = timing.parse_grouped_output(
            grouped_stdout(task=task), "forward", task, 4096, 8)
        self.assertEqual(parsed.preamble["staircase_rows"], "190")
        staircase = mutate_preamble(
            grouped_stdout(task=task), "staircase_rows", "191")
        staircase = mutate_preamble(staircase, "columns", "32215")
        with self.assertRaisesRegex(timing.TimingError, "staircase_rows"):
            timing.parse_grouped_output(
                staircase, "forward", task, 4096, 8)

    def test_neutral_smoke_replay_table_covers_every_record(self):
        self.assertEqual(timing.NEUTRAL_SMOKE_SPECS, (
            ("exact_base", "exact_base", False, "shared-x", "primary"),
            ("exact_architecture", "exact_architecture", False,
             "shared-x", "primary"),
            ("measurement_shared_x", "measurement", True,
             "shared-x", "primary"),
            ("measurement_frozen", "measurement", True,
             "frozen", "primary"),
            ("measurement_frozen_K10000", "measurement", True,
             "frozen", "boundary"),
        ))
        self.assertEqual(
            {record_key for record_key, *_rest in timing.NEUTRAL_SMOKE_SPECS},
            {"exact_base", "exact_architecture", "measurement_shared_x",
             "measurement_frozen", "measurement_frozen_K10000"})

    def test_outer_order_swaps_one_measurement_binary(self):
        design = {
            "root": "/tmp/frozen-campaign", "core": 8, "numa_node": 0,
            "evict_bytes": 4096,
            "tools": {name: {"path": "/usr/bin/" + name}
                      for name in ("env", "taskset", "numactl")},
        }
        forward = timing.command_for(design, one_task(), "forward")
        reverse = timing.command_for(design, one_task(), "reverse")
        self.assertEqual(timing.OUTER_ORDER.count("A"), 4)
        self.assertEqual(timing.OUTER_ORDER.count("B"), 4)
        self.assertIn(timing.BINARY_NAMES["measurement"],
                      [Path(value).name for value in forward])
        self.assertIn(timing.BINARY_NAMES["measurement"],
                      [Path(value).name for value in reverse])
        self.assertEqual(
            forward[forward.index("--control-degree-balanced-staircase") + 1],
            "0")
        self.assertEqual(
            reverse[reverse.index("--control-degree-balanced-staircase") + 1],
            "1")
        self.assertEqual(
            forward[forward.index("--candidate-degree-balanced-staircase") + 1],
            "1")
        self.assertEqual(
            reverse[reverse.index("--candidate-degree-balanced-staircase") + 1],
            "0")

    def test_cross_cache_identity_is_exact(self):
        ledger = {}
        parsed = timing.parse_grouped_output(
            grouped_stdout(), "forward", one_task(), 4096, 8)
        for task in timing.generate_tasks():
            timing._register_cross_cache_identity(ledger, task, parsed)
        timing._validate_cross_cache_ledger(ledger)
        mismatched = {}
        cold = one_task(cache_state="cold")
        warm = one_task(cache_state="warm")
        cold_parsed = timing.parse_grouped_output(
            grouped_stdout(task=cold), "forward", cold, 4096, 8)
        warm_raw = mutate_preamble(
            grouped_stdout(task=warm), "trace_sha256", "b" * 64)
        warm_raw = mutate_preamble(warm_raw, "packet_trace_sha256", "b" * 64)
        warm_parsed = timing.parse_grouped_output(
            warm_raw, "forward", warm, 4096, 8)
        timing._register_cross_cache_identity(mismatched, cold, cold_parsed)
        with self.assertRaisesRegex(timing.TimingError, "cold/warm"):
            timing._register_cross_cache_identity(mismatched, warm, warm_parsed)


class OutputParserTests(unittest.TestCase):
    def test_exact_legacy_and_measurement_neutral_normalize(self):
        task = one_task()
        measurement_raw = grouped_stdout("neutral", task, geometry="shared-x")
        legacy_raw = legacy_stdout_from_neutral_measurement(measurement_raw)
        measurement = timing._parse_neutral_replay(
            measurement_raw, True, task, expected_geometry="shared-x")
        legacy = timing._parse_neutral_replay(legacy_raw, False, task)
        self.assertEqual(measurement["shared_sha256"], legacy["shared_sha256"])
        frozen = timing._parse_neutral_replay(
            grouped_stdout("neutral", task, geometry="frozen"), True, task,
            expected_geometry="frozen")
        self.assertEqual(
            frozen["measurement_receipts"]["geometry"], "frozen")
        self.assertEqual(
            frozen["measurement_receipts"]["rhs_route_expected"], "streamed")
        tampered = mutate_csv(measurement_raw, 3, "block_xors", 74096)
        with self.assertRaisesRegex(timing.TimingError, "internal arms differ"):
            timing._parse_neutral_replay(
                tampered, True, task, expected_geometry="shared-x")
        selector = mutate_csv(
            measurement_raw, 3, "degree_balanced_staircase", "1")
        with self.assertRaisesRegex(timing.TimingError, "coordinate changed"):
            timing._parse_neutral_replay(
                selector, True, task, expected_geometry="shared-x")

    def test_forward_reverse_normalize_to_one_semantic_record(self):
        task = one_task()
        forward = timing.parse_grouped_output(
            grouped_stdout("forward", task), "forward", task, 4096, 8)
        reverse = timing.parse_grouped_output(
            grouped_stdout("reverse", task), "reverse", task, 4096, 8)
        self.assertEqual(forward.semantic_sha256, reverse.semantic_sha256)
        self.assertEqual(forward.work_signatures, reverse.work_signatures)
        self.assertEqual(forward.outcomes, {"base": 0, "balanced": 0})
        self.assertEqual(forward.contaminations, ())
        self.assertGreater(forward.base_timed_elapsed_ns,
                           forward.balanced_timed_elapsed_ns)

    def test_only_staircase_graph_may_differ(self):
        task = one_task()
        raw = mutate_preamble(
            grouped_stdout(task=task), "candidate_dense_fingerprint",
            "0x7777777777777777")
        with self.assertRaisesRegex(timing.TimingError, "isolate the staircase"):
            timing.parse_grouped_output(raw, "forward", task, 4096, 8)

    def test_coefficient_and_actual_route_receipts_are_exact(self):
        task = one_task()
        coefficient = mutate_preamble(
            grouped_stdout(task=task),
            "candidate_mixed_coefficient_fingerprint",
            "0x8888888888888888")
        with self.assertRaisesRegex(timing.TimingError, "preamble mismatch"):
            timing.parse_grouped_output(
                coefficient, "forward", task, 4096, 8)
        dense_rows = mutate_preamble(
            grouped_stdout(task=task), "dense_rows", "13")
        with self.assertRaisesRegex(timing.TimingError, "preamble mismatch"):
            timing.parse_grouped_output(
                dense_rows, "forward", task, 4096, 8)
        preflight_route = mutate_preamble(
            grouped_stdout(task=task), "candidate_preflight_rhs_route",
            "dual")
        with self.assertRaisesRegex(timing.TimingError, "RHS route"):
            timing.parse_grouped_output(
                preflight_route, "forward", task, 4096, 8)
        row_route = mutate_csv(
            grouped_stdout(task=task), 1, "rhs_route_actual", "dual")
        with self.assertRaisesRegex(timing.TimingError, "row 1 mismatch"):
            timing.parse_grouped_output(row_route, "forward", task, 4096, 8)
        joint_work = mutate_csv(
            grouped_stdout(task=task), 1, "joint_source_xors", 1)
        with self.assertRaisesRegex(timing.TimingError, "row 1 mismatch"):
            timing.parse_grouped_output(joint_work, "forward", task, 4096, 8)
        raw = mutate_preamble(
            grouped_stdout(task=task), "candidate_staircase_fingerprint",
            "0x1111111111111111")
        with self.assertRaisesRegex(timing.TimingError, "isolate the staircase"):
            timing.parse_grouped_output(raw, "forward", task, 4096, 8)

    def test_work_drift_within_an_arm_is_rejected(self):
        task = one_task()
        raw = mutate_csv(grouped_stdout(task=task), 8, "block_xors", 74096)
        with self.assertRaisesRegex(timing.TimingError, "work.*unstable"):
            timing.parse_grouped_output(raw, "forward", task, 4096, 8)

    def test_non_common_success_is_retained_unconditionally(self):
        task = one_task()
        parsed = timing.parse_grouped_output(
            grouped_stdout(task=task, balanced_result=1),
            "forward", task, 4096, 8)
        self.assertEqual(parsed.outcomes, {"base": 0, "balanced": 1})
        self.assertEqual(timing._semantic_outcome_class(parsed.outcomes),
                         "base-only")
        self.assertGreater(parsed.base_timed_elapsed_ns, 0)
        self.assertGreater(parsed.balanced_timed_elapsed_ns, 0)

    def test_common_success_decoded_output_must_match(self):
        task = one_task()
        raw = mutate_preamble(
            grouped_stdout(task=task),
            "preflight_candidate_decoded_fingerprint",
            "0x7777777777777777")
        with self.assertRaisesRegex(timing.TimingError, "decoded output differs"):
            timing.parse_grouped_output(raw, "forward", task, 4096, 8)

    def test_migration_saturation_and_faults_are_contamination(self):
        task = one_task()
        raw = mutate_csv(grouped_stdout(task=task), 9, "cpu_after", 9)
        raw = mutate_csv(raw, 10, "saturated", 1)
        raw = mutate_csv(raw, 11, "minflt_delta", 65)
        raw = mutate_csv(raw, 11, "fault_contaminated", 1)
        raw = mutate_csv(raw, 12, "majflt_delta", 1)
        raw = mutate_csv(raw, 12, "fault_contaminated", 1)
        parsed = timing.parse_grouped_output(raw, "forward", task, 4096, 8)
        self.assertTrue(any("migration" in item for item in parsed.contaminations))
        self.assertTrue(any("saturated" in item for item in parsed.contaminations))
        self.assertTrue(any("minor-fault" in item for item in parsed.contaminations))
        self.assertTrue(any("major-fault" in item for item in parsed.contaminations))

    def test_signed_fault_receipt_cannot_be_forged(self):
        task = one_task()
        raw = mutate_csv(grouped_stdout(task=task), 0, "minflt_delta", -1)
        with self.assertRaisesRegex(timing.TimingError, "fault receipt disagrees"):
            timing.parse_grouped_output(raw, "forward", task, 4096, 8)

    def test_schema_and_preamble_order_are_exact(self):
        task = one_task()
        raw = mutate_preamble(grouped_stdout(task=task), "schema", "v2")
        with self.assertRaisesRegex(timing.TimingError, "wrong groupedtiming schema"):
            timing.parse_grouped_output(raw, "forward", task, 4096, 8)
        lines = grouped_stdout(task=task).decode("ascii").splitlines()
        tokens = lines[0].split(" ")
        tokens[-1], tokens[-2] = tokens[-2], tokens[-1]
        lines[0] = " ".join(tokens)
        with self.assertRaisesRegex(timing.TimingError, "order/schema"):
            timing.parse_grouped_output(
                ("\n".join(lines) + "\n").encode("ascii"),
                "forward", task, 4096, 8)

    def test_quoted_or_extra_csv_fields_are_noncanonical(self):
        task = one_task()
        quoted = grouped_stdout(task=task).replace(b",burst,", b',"burst",', 1)
        with self.assertRaisesRegex(timing.TimingError, "canonical LF text"):
            timing.parse_grouped_output(quoted, "forward", task, 4096, 8)
        lines = grouped_stdout(task=task).decode("ascii").splitlines()
        lines[2] += ",ignored"
        extra = ("\n".join(lines) + "\n").encode("ascii")
        with self.assertRaisesRegex(timing.TimingError, "field count"):
            timing.parse_grouped_output(extra, "forward", task, 4096, 8)


class ReceiptStatisticsAndThermalTests(unittest.TestCase):
    def test_sealed_record_detects_mutation(self):
        value = timing.sealed_record("example.v1", {"number": 4})
        timing.verify_sealed_record(value, "example.v1", "fixture")
        value["number"] = 5
        with self.assertRaisesRegex(timing.TimingError, "self-hash mismatch"):
            timing.verify_sealed_record(value, "example.v1", "fixture")

    def test_execution_receipt_schema_matches_constructor(self):
        task = one_task()
        parsed = timing.parse_grouped_output(
            grouped_stdout(task=task), "forward", task, 4096, 8)
        receipt = timing._execution_receipt(
            task, 0, "forward", ["example"], 0,
            "2026-07-18T00:00:00.000Z", 10, 20, parsed, b"", [], "c" * 64)
        self.assertEqual(set(receipt), timing.EXECUTION_RECEIPT_FIELDS)

    def test_bootstrap_clusters_cold_and_warm_together(self):
        records = [
            {"base": 100, "candidate": 90, "cluster": [3200, 64, "burst", 1]},
            {"base": 110, "candidate": 100, "cluster": [3200, 64, "burst", 1]},
            {"base": 200, "candidate": 210, "cluster": [3200, 64, "burst", 2]},
            {"base": 210, "candidate": 215, "cluster": [3200, 64, "burst", 2]},
        ]
        first = timing._bootstrap(records, "unit-test", repetitions=200)
        second = timing._bootstrap(records, "unit-test", repetitions=200)
        self.assertEqual(first, second)
        self.assertEqual(first["cluster_count"], 2)
        self.assertLessEqual(float(first["lower_95"]),
                             float(first["upper_95"]))

    def test_protected_sampler_fallbacks_use_frozen_tools(self):
        tools = {
            "sudo": {"path": "/frozen/sudo"},
            "kill": {"path": "/frozen/kill"},
            "cat": {"path": "/frozen/cat"},
        }
        completed = mock.Mock(returncode=0, stdout=b"sample", stderr=b"")
        with mock.patch.object(timing.os, "kill", side_effect=PermissionError), \
                mock.patch.object(timing.subprocess, "run", return_value=completed) \
                as run:
            self.assertTrue(timing._pid_alive(123, tools))
        self.assertEqual(
            run.call_args.args[0],
            ("/frozen/sudo", "-n", "/frozen/kill", "-0", "123"))
        with mock.patch.object(Path, "read_bytes", side_effect=PermissionError), \
                mock.patch.object(timing.subprocess, "run", return_value=completed) \
                as run:
            self.assertEqual(timing._proc_bytes(123, "status", tools), b"sample")
        self.assertEqual(
            run.call_args.args[0],
            ("/frozen/sudo", "-n", "/frozen/cat", "/proc/123/status"))

    def test_sampler_metadata_rejects_pid_reuse_splice(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            script = root / "wirehair_expo_thermal_sampler.py"
            source = root / "thermal.csv"
            pid_file = root / "thermal.pid"
            script.write_text("# fixture\n", encoding="ascii")
            source.write_text("fixture\n", encoding="ascii")
            pid_file.write_text("123\n", encoding="ascii")
            cmdline = b"\0".join((
                b"/usr/bin/python3", str(script).encode(), b"--csv",
                str(source).encode(), b"--pid-file", str(pid_file).encode(),
                b""))
            design = {
                "controller_core": 126,
                "sampler_core": 127,
                "topology": {"llc_shared_cpus": [8, 9]},
                "tools": {
                    "sudo": {"path": "/frozen/sudo"},
                    "fuser": {"path": "/frozen/fuser"},
                    "kill": {"path": "/frozen/kill"},
                    "cat": {"path": "/frozen/cat"},
                },
            }
            fuser = mock.Mock(returncode=0, stdout=b"123\n", stderr=b"")
            proc_values = {
                "status": b"Name:\tpython3\nCpus_allowed_list:\t127\n",
                "cmdline": cmdline,
            }
            with mock.patch.object(timing, "_pid_alive", return_value=True), \
                    mock.patch.object(
                        timing, "_process_start_ticks", side_effect=[10, 11]), \
                    mock.patch.object(
                        timing, "_proc_bytes",
                        side_effect=lambda _pid, name, _tools: proc_values[name]), \
                    mock.patch.object(
                        timing.subprocess, "run", return_value=fuser):
                with self.assertRaisesRegex(
                        timing.TimingError, "identity changed"):
                    timing._thermal_reader(design, pid_file)

    def test_sampler_affinity_must_be_exact_cpu127(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            pid_file = root / "thermal.pid"
            pid_file.write_text("123\n", encoding="ascii")
            design = {
                "controller_core": timing.CONTROLLER_CORE,
                "sampler_core": timing.SAMPLER_CORE,
                "topology": {"llc_shared_cpus": [8, 72]},
                "tools": {
                    "sudo": {"path": "/frozen/sudo"},
                    "fuser": {"path": "/frozen/fuser"},
                    "kill": {"path": "/frozen/kill"},
                    "cat": {"path": "/frozen/cat"},
                },
            }
            fuser = mock.Mock(returncode=0, stdout=b"123\n", stderr=b"")
            with mock.patch.object(timing, "_pid_alive", return_value=True), \
                    mock.patch.object(
                        timing, "_process_start_ticks", return_value=10), \
                    mock.patch.object(
                        timing, "_proc_bytes",
                        return_value=(
                            b"Name:\tpython3\nCpus_allowed_list:\t120\n")), \
                    mock.patch.object(
                        timing.subprocess, "run", return_value=fuser):
                with self.assertRaisesRegex(
                        timing.TimingError, "frozen CPU 127"):
                    timing._thermal_reader(design, pid_file)

    def test_thermal_interval_is_exactly_bracketed_and_bounded(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "thermal.csv"
            path.write_bytes(thermal_csv())
            raw, summary = timing.collect_thermal_interval(path, 10.25, 11.75)
        rows, _lines = timing._parse_thermal(raw)
        self.assertEqual([float(row["monotonic_s"]) for row in rows],
                         [10.0, 11.0, 12.0])
        self.assertEqual(summary["sample_count"], 3)
        self.assertEqual(summary["edac_ce_delta"], 0)

    def test_thermal_snapshot_retries_partial_append(self):
        complete = thermal_csv()
        with mock.patch.object(
                timing, "stable_bytes", side_effect=[complete[:-3], complete]), \
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

    def test_mid_interval_edac_increment_then_reset_fails_closed(self):
        text = thermal_csv().decode("ascii").splitlines()
        header = text[0].split(",")
        fields = text[2].split(",")
        fields[header.index("edac_ce")] = "5"
        text[2] = ",".join(fields)
        raw = ("\r\n".join(text) + "\r\n").encode("ascii")
        with self.assertRaisesRegex(timing.TimingError, "EDAC counters"):
            timing.validate_sealed_thermal_interval(raw, 9.25, 11.75)


if __name__ == "__main__":
    unittest.main()
