#!/usr/bin/env python3
"""Unit tests for the sealed segmented-anchor timing controller."""

from __future__ import annotations

from dataclasses import replace
import csv
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock


sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_segmented_anchor_solve_timing as timing  # noqa: E402


CORE = 8


def fixture(task: dict[str, object], *, candidate_elapsed: int = 99,
            control_elapsed: int = 100) -> bytes:
    candidate = str(task["candidate"])
    K = int(task["K"])
    bb = int(task["bb"])
    staircase = timing.STAIRCASE_ROWS[K]
    source_hits = 2 if K < 10000 else 3
    intermediate = (K + staircase + 12 + 12) * bb
    preamble = {
        "schema": "v3", "policy": "h12-q0-grouped", "timing_scope": "solve",
        "cycles": "4", "order": timing.ORDER, "discard_cycle": "0",
        "cycle_mode": "full", "cycle_index": "all", "N": str(K),
        "bb": str(bb), "overhead": str(timing.OVERHEAD), "loss": timing.LOSS,
        "seed": str(task["seed"]), "schedule": str(task["schedule"]),
        "cache_state": str(task["cache_state"]), "overhead_stream": "salted",
        "evict_bytes": "268435456", "eviction_prefaulted": "1",
        "control_period": str(timing.PERIOD), "control_grouped_rows": "0",
        "control_buckets": timing.BUCKETS, "control_grouped_hash_seed": "0x0",
        "control_final_h_a_columns": "0", "candidate_period": str(timing.PERIOD),
        "candidate_grouped_rows": "0", "candidate_buckets": timing.BUCKETS,
        "candidate_grouped_hash_seed": "0x0", "candidate_final_h_a_columns": "0",
        "gf256_rows": "10", "gf16_rows": "2",
        "control_dense_layout": "two-anchor", "candidate_dense_layout": candidate,
        "dense_layout_is_only_architecture_selector": "1",
        "control_staircase_rows": str(staircase), "control_dense_rows": "12",
        "control_heavy_rows": "12", "control_source_hits": str(source_hits),
        "control_field": "1", "control_dense_identity_corner": "0",
        "control_dense_two_anchor_exact": "1",
        "control_dense_two_anchor_phase": "0",
        "control_segmented_dense_anchors": "none", "control_heavy_family": "0",
        "control_mix_count": "2", "candidate_staircase_rows": str(staircase),
        "candidate_dense_rows": "12", "candidate_heavy_rows": "12",
        "candidate_source_hits": str(source_hits), "candidate_field": "1",
        "candidate_dense_identity_corner": "0",
        "candidate_dense_two_anchor_exact": "0",
        "candidate_dense_two_anchor_phase": "0",
        "candidate_segmented_dense_anchors": candidate,
        "candidate_heavy_family": "0", "candidate_mix_count": "2",
        "control_attempt": "0", "control_matrix_seed": "0x1234",
        "control_peel_seed": "0x5678", "candidate_attempt": "0",
        "candidate_matrix_seed": "0x1234", "candidate_peel_seed": "0x5678",
        "mix": "2", "payload": "distinct-packet-zero-v1",
        "payload_count": str(K + timing.OVERHEAD),
        "payload_bytes": str((K + timing.OVERHEAD) * bb),
        "payload_alignment": "64", "payload_prefaulted": "1",
        "system_build": "outside-timer",
        "tls_reapply": "full-per-slot-outside-timer",
        "allocator_tls_state": "preflight-warmed",
        "solve_value_storage": "owned-noinit", "solve_value_publish": "swap",
        "preflight_control_result": "0", "preflight_candidate_result": "0",
        "cell_class": "common-success", "common_success": "1",
        "trace_sha256": "1" * 64,
    }
    if tuple(preamble) != timing.PREAMBLE_FIELDS:
        raise AssertionError("fixture preamble order drifted")
    output = io.StringIO(newline="")
    output.write("# groupedtiming: " + " ".join(
        "%s=%s" % item for item in preamble.items()) + "\n")
    writer = csv.DictWriter(output, fieldnames=timing.CSV_FIELDS,
                            lineterminator="\n")
    writer.writeheader()
    for cycle in range(4):
        for slot, letter in enumerate(timing.ORDER):
            arm = "control" if letter == "A" else "candidate"
            elapsed = control_elapsed if arm == "control" else candidate_elapsed
            row = {field: "0" for field in timing.CSV_FIELDS}
            row.update({
                "N": str(K), "bb": str(bb), "overhead": str(timing.OVERHEAD),
                "schedule": str(task["schedule"]), "seed": str(task["seed"]),
                "loss": timing.LOSS, "cache_state": str(task["cache_state"]),
                "cycle": str(cycle), "slot": str(slot), "arm": arm,
                "period": str(timing.PERIOD), "grouped_rows": "0",
                "buckets_requested": timing.BUCKETS, "seed_attempt": "0",
                "matrix_seed": "0x1234", "peel_seed": "0x5678",
                "preflight_result": "0", "cell_class": "common-success",
                "common_success": "1", "result": "0", "outcome_stable": "1",
                "elapsed_ns": str(elapsed), "saturated": "0",
                "cpu_before": str(CORE), "cpu_after": str(CORE),
                "cpu_migrated": "0", "minflt_delta": "0", "majflt_delta": "0",
                "fault_contaminated": "0", "inactivated": "20",
                "binary_def": "0", "heavy_gain": "12",
                "block_xors": "1000" if arm == "control" else "1001",
                "block_muladds": "200", "build_ns": "10", "peel_ns": "20",
                "project_ns": "30", "residual_ns": "40", "backsub_ns": "50",
                "joint_source_xors": "10", "joint_marginal_xors": "11",
                "joint_marginal_copies": "12", "joint_active_deltas": "13",
                "joint_scratch_bytes": "64", "dual_source_columns": "0",
                "source_bytes": str(K * bb),
                "packet_payload_bytes": str((K + timing.OVERHEAD) * bb),
                "intermediate_bytes": str(intermediate),
                "solve_value_arena_bytes": str(intermediate),
                "solve_value_eager_zero_bytes": "0",
                "solve_value_commit_copy_bytes": "0",
            })
            writer.writerow(row)
    return output.getvalue().encode("ascii")


def mutate_csv(raw: bytes, callback) -> bytes:
    lines = raw.decode("ascii").splitlines()
    reader = csv.DictReader(io.StringIO("\n".join(lines[1:]) + "\n"))
    rows = [dict(row) for row in reader]
    callback(rows)
    output = io.StringIO(newline="")
    output.write(lines[0] + "\n")
    writer = csv.DictWriter(output, fieldnames=timing.CSV_FIELDS,
                            lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue().encode("ascii")


def mutate_preamble(raw: bytes, replacements: dict[str, str]) -> bytes:
    lines = raw.decode("ascii").splitlines()
    prefix = "# groupedtiming: "
    values = [token.split("=", 1) for token in lines[0][len(prefix):].split(" ")]
    for pair in values:
        if pair[0] in replacements:
            pair[1] = replacements[pair[0]]
    lines[0] = prefix + " ".join("=".join(pair) for pair in values)
    return ("\n".join(lines) + "\n").encode("ascii")


class SegmentedAnchorTimingTest(unittest.TestCase):
    def setUp(self) -> None:
        self.task = {
            "K": 3200, "bb": 64, "seed_index": 0,
            "seed": timing.SEEDS[0], "schedule": "burst",
            "candidate": "four-0369", "cache_state": "warm",
        }
        self.raw = fixture(self.task)

    def parse(self, raw: bytes | None = None) -> timing.ParsedOutput:
        return timing.parse_grouped_output(
            self.raw if raw is None else raw, self.task, CORE,
            expected_evict_bytes=268435456)

    def test_cartesian_manifest_and_cluster_adjacency(self) -> None:
        tasks = timing.generate_tasks()
        self.assertEqual(len(tasks), 648)
        self.assertEqual([task["job"] for task in tasks], list(range(648)))
        self.assertEqual(len({task["task_id"] for task in tasks}), 648)
        for offset in range(0, len(tasks), 4):
            cluster = tasks[offset:offset + 4]
            coordinates = {(item["K"], item["bb"], item["seed_index"],
                            item["seed"], item["schedule"]) for item in cluster}
            members = {(item["candidate"], item["cache_state"])
                       for item in cluster}
            self.assertEqual(len(coordinates), 1)
            self.assertEqual(members, set(
                (candidate, cache) for candidate in timing.CANDIDATES
                for cache in timing.CACHE_STATES))

    def test_exact_parser_and_retained_cycles(self) -> None:
        parsed = self.parse()
        self.assertTrue(parsed.common_success)
        self.assertEqual(parsed.cell_class, "common-success")
        self.assertEqual(parsed.timed_control_ns, 12 * 100)
        self.assertEqual(parsed.timed_candidate_ns, 12 * 99)
        self.assertEqual(parsed.preamble["candidate_dense_layout"], "four-0369")
        self.assertEqual(parsed.contaminations, ())

    def test_parser_rejects_layout_eviction_and_width_drift(self) -> None:
        bad_layout = mutate_preamble(
            self.raw, {"candidate_dense_layout": "three-048"})
        with self.assertRaisesRegex(timing.TimingError, "preamble mismatch"):
            self.parse(bad_layout)
        with self.assertRaisesRegex(timing.TimingError, "eviction allocation"):
            self.parse(mutate_preamble(self.raw, {"evict_bytes": "4096"}))
        with self.assertRaisesRegex(timing.TimingError, "preamble mismatch"):
            self.parse(mutate_preamble(
                self.raw, {"control_staircase_rows": "63",
                           "candidate_staircase_rows": "63"}))
        bad_arena = mutate_csv(self.raw, lambda rows: rows[0].update({
            "intermediate_bytes": "1", "solve_value_arena_bytes": "1"}))
        with self.assertRaisesRegex(timing.TimingError, "row 0 mismatch"):
            self.parse(bad_arena)
        lines = self.raw.decode("ascii").splitlines()
        lines[2] += ",unexpected"
        with self.assertRaisesRegex(timing.TimingError, "row width"):
            self.parse(("\n".join(lines) + "\n").encode("ascii"))
        with self.assertRaisesRegex(timing.TimingError, "noncanonical"):
            self.parse(self.raw.replace(b"\n3200,", b'\n"3200",', 1))
        with self.assertRaisesRegex(timing.TimingError, "parser domain"):
            timing.parse_grouped_output(
                self.raw, self.task, CORE, replacement_cycle=4)

    def test_parser_derives_outcome_and_restricts_codec_results(self) -> None:
        inconsistent = mutate_preamble(
            self.raw, {"preflight_control_result": "1"})
        with self.assertRaisesRegex(timing.TimingError, "cell class disagrees"):
            self.parse(inconsistent)
        failed = mutate_preamble(self.raw, {
            "preflight_control_result": "1",
            "preflight_candidate_result": "1",
            "cell_class": "common-failure", "common_success": "0",
        })
        failed = mutate_csv(failed, lambda rows: [row.update({
            "preflight_result": "1", "result": "1",
            "cell_class": "common-failure", "common_success": "0",
        }) for row in rows])
        parsed = self.parse(failed)
        self.assertFalse(parsed.common_success)
        self.assertEqual(parsed.cell_class, "common-failure")
        invalid = mutate_preamble(self.raw, {
            "preflight_control_result": "-1",
            "preflight_candidate_result": "-1",
            "cell_class": "common-failure", "common_success": "0",
        })
        with self.assertRaisesRegex(timing.TimingError, "canonical uint"):
            self.parse(invalid)

    def test_contamination_receipts(self) -> None:
        migrated = mutate_csv(self.raw, lambda rows: rows[7].update({
            "cpu_after": "9", "cpu_migrated": "1"}))
        parsed = self.parse(migrated)
        self.assertIn("row7:migration:8:9:1", parsed.contaminations)
        inconsistent_migration = mutate_csv(
            self.raw, lambda rows: rows[7].update({"cpu_migrated": "1"}))
        with self.assertRaisesRegex(timing.TimingError, "migration receipt"):
            self.parse(inconsistent_migration)
        bad_fault = mutate_csv(
            self.raw, lambda rows: rows[0].update({"fault_contaminated": "1"}))
        with self.assertRaisesRegex(timing.TimingError, "fault contamination"):
            self.parse(bad_fault)
        negative_fault = mutate_csv(self.raw, lambda rows: rows[0].update({
            "minflt_delta": "-2", "fault_contaminated": "-1"}))
        with self.assertRaisesRegex(timing.TimingError, "fault counter delta"):
            self.parse(negative_fault)

    def test_aggregate_pairing_and_speed_gates(self) -> None:
        tasks = timing.generate_tasks()
        parsed = {}
        by_candidate = {}
        for candidate in timing.CANDIDATES:
            task = {**self.task, "candidate": candidate}
            by_candidate[candidate] = timing.parse_grouped_output(
                fixture(task), task, CORE, expected_evict_bytes=268435456)
        for task in tasks:
            parsed[int(task["job"])] = by_candidate[str(task["candidate"])]
        original = timing.BOOTSTRAP_REPETITIONS
        timing.BOOTSTRAP_REPETITIONS = 200
        try:
            result = timing.aggregate_rows(tasks, parsed)
        finally:
            timing.BOOTSTRAP_REPETITIONS = original
        for candidate in timing.CANDIDATES:
            overall = result["candidates"][candidate]["overall"]
            self.assertEqual(overall["cell_count"], 324)
            self.assertAlmostEqual(overall["ratio_of_sums"], 0.99)
            self.assertTrue(overall["speed_gate"]["speed_gate_passed"])
            self.assertEqual(overall["bootstrap"]["cluster_count"], 162)
        pareto = result["candidate_pareto_comparison"]
        self.assertEqual(pareto["baseline_normalized_ratio_of_ratios"], 1.0)
        self.assertFalse(pareto["select_three_048_over_four_0369"])
        gates = {candidate: result["candidates"][candidate]["overall"][
            "speed_gate"] for candidate in timing.CANDIDATES}
        self.assertEqual(
            timing.select_pareto_architecture(gates, pareto), "four-0369")

        fast_three = {}
        elapsed = {"three-048": 98, "four-0369": 100}
        for candidate in timing.CANDIDATES:
            task = {**self.task, "candidate": candidate}
            fast_three[candidate] = timing.parse_grouped_output(
                fixture(task, candidate_elapsed=elapsed[candidate]),
                task, CORE, expected_evict_bytes=268435456)
        timing.BOOTSTRAP_REPETITIONS = 200
        try:
            fast_result = timing.aggregate_rows(
                tasks, {int(task["job"]): fast_three[str(task["candidate"])]
                        for task in tasks})
        finally:
            timing.BOOTSTRAP_REPETITIONS = original
        fast_pareto = fast_result["candidate_pareto_comparison"]
        self.assertAlmostEqual(
            fast_pareto["baseline_normalized_ratio_of_ratios"], 0.98)
        self.assertTrue(fast_pareto["select_three_048_over_four_0369"])
        fast_gates = {candidate: fast_result["candidates"][candidate][
            "overall"]["speed_gate"] for candidate in timing.CANDIDATES}
        self.assertEqual(timing.select_pareto_architecture(
            fast_gates, fast_pareto), "three-048")

        # The cross-architecture decision normalizes each candidate by its
        # cotimed identical baseline, so a deliberate process-scale drift
        # cannot reverse the architectural comparison.
        drifted = {}
        for candidate in timing.CANDIDATES:
            task = {**self.task, "candidate": candidate}
            control_ns, candidate_ns = (
                (200, 198) if candidate == "three-048" else (100, 100))
            drifted[candidate] = timing.parse_grouped_output(
                fixture(task, candidate_elapsed=candidate_ns,
                        control_elapsed=control_ns),
                task, CORE, expected_evict_bytes=268435456)
        timing.BOOTSTRAP_REPETITIONS = 200
        try:
            drifted_result = timing.aggregate_rows(
                tasks, {int(task["job"]): drifted[str(task["candidate"])]
                        for task in tasks})
        finally:
            timing.BOOTSTRAP_REPETITIONS = original
        drifted_pareto = drifted_result["candidate_pareto_comparison"]
        self.assertAlmostEqual(
            drifted_pareto["raw_candidate_ratio_of_sums"], 1.98)
        self.assertAlmostEqual(
            drifted_pareto["baseline_normalized_ratio_of_ratios"], 0.99)
        self.assertTrue(
            drifted_pareto["select_three_048_over_four_0369"])
        broken = dict(parsed)
        job = next(job for job, value in enumerate(tasks)
                   if value["candidate"] == "four-0369" and
                   value["cache_state"] == "warm")
        broken[job] = replace(broken[job], trace_sha256="2" * 64)
        with self.assertRaisesRegex(timing.TimingError, "cold/warm"):
            timing.aggregate_rows(tasks, broken)

    def test_recovery_receipt_is_exactly_bound(self) -> None:
        path = Path(__file__).with_name(timing.RECOVERY_RESULT_NAME)
        value = json.loads(path.read_text(encoding="utf-8"))
        validated = timing._validate_recovery_result(value)
        self.assertEqual(validated["arms"]["four-0369"]["failures"], 650)
        corrupted = json.loads(json.dumps(value))
        corrupted["arms"]["four-0369"]["failures"] = 651
        with self.assertRaisesRegex(timing.TimingError, "content changed"):
            timing._validate_recovery_result(corrupted)

    def test_imported_helper_provenance_and_codec_isolation(self) -> None:
        repo = Path(__file__).resolve().parent.parent
        provenance = timing.P32_HELPER_SOURCE_PROVENANCE
        self.assertEqual(provenance["source_commit"],
                         "97a3a0b941f3efe34a6e18609b7681d3a1110982")
        self.assertEqual(provenance["source_tree"],
                         "ab84c084b212fecc6672a7efda5a4df6203c866c")
        for relative, receipt in provenance["files"].items():
            raw = (repo / relative).read_bytes()
            git_object = b"blob " + str(len(raw)).encode("ascii") + b"\0" + raw
            self.assertEqual(hashlib.sha1(git_object).hexdigest(),
                             receipt["git_blob_sha1"])
            self.assertEqual(hashlib.sha256(raw).hexdigest(),
                             receipt["sha256"])

        # The imported branch also contained unrelated P32 codec work.  The
        # timing-only segmented branch must retain these four core files
        # byte-for-byte from its 9df95d6 recovery checkpoint and must not
        # acquire the P32 groupedtiming geometry/RHS interface.
        unchanged_codec_sha256 = {
            "codec/WirehairV2Precode.cpp":
                "2c671d21665c7b3cd711cec0f424b5563e52e3bf50b9c56188073439403463e9",
            "codec/WirehairV2Precode.h":
                "3614827cf7aa4733be4fb010f16663798d66b98264834b9b88ef2cd8bf2549ee",
            "codec/WirehairV2Solve.cpp":
                "f30329c60a95f27eb21dd56e1a86f5b77f2f2da2f117bd2dfd6d66dcff7572f3",
            "codec/WirehairV2Solve.h":
                "9ca07438d5f192bd41d4db9490ab9fb24953cd0fe12a647930e2166db96c1f49",
        }
        for relative, digest in unchanged_codec_sha256.items():
            self.assertEqual(hashlib.sha256(
                (repo / relative).read_bytes()).hexdigest(), digest)
        grouped_sources = b"\n".join(
            (repo / relative).read_bytes() for relative in (
                "codec/WirehairV2Bench.cpp", "codec/V2BenchCliTest.cmake"))
        for token in (b"--control-geometry", b"--candidate-geometry",
                      b"GroupedTimingRhsRouteName",
                      b"GroupedTimingExpectedRhsRoute"):
            self.assertNotIn(token, grouped_sources)

    def test_sealed_records_and_busy_ticks(self) -> None:
        value = timing.sealed_record("example.v1", {"answer": 42})
        self.assertEqual(timing.verify_sealed(value, "example.v1", "example"), value)
        corrupted = dict(value)
        corrupted["answer"] = 43
        with self.assertRaisesRegex(timing.TimingError, "self hash mismatch"):
            timing.verify_sealed(corrupted, "example.v1", "example")
        self.assertEqual(
            timing.busy_ticks((10, 0, 5, 100, 20), (14, 0, 8, 110, 22)), 7)
        self.assertEqual(timing.cpu_busy_tick_deltas(
            {"9": (10, 0, 5, 100, 20)},
            {"9": (12, 0, 5, 110, 20)}), {"9": 2})
        with self.assertRaisesRegex(timing.TimingError, "regressed"):
            timing.busy_ticks(
                (10, 0, 5, 100, 20), (9, 0, 8, 110, 22))

    def test_non_sibling_timing_llc_peer_activity_is_rejected(self) -> None:
        cpus = [8, 9, 72, 73]
        samples = {
            8: [(10, 0, 5, 100, 20), (10, 0, 5, 110, 20)],
            9: [(10, 0, 5, 100, 20), (12, 0, 5, 110, 20)],
            72: [(10, 0, 5, 100, 20), (10, 0, 5, 110, 20)],
            73: [(10, 0, 5, 100, 20), (10, 0, 5, 110, 20)],
        }

        def read_ticks(cpu: int) -> tuple[int, ...]:
            return samples[cpu].pop(0)

        with mock.patch.object(timing, "cpu_ticks", side_effect=read_ticks), \
                mock.patch.object(timing.time, "sleep"), \
                mock.patch.object(
                    timing.time, "monotonic", side_effect=(100.0, 101.0)):
            receipt = timing.measure_quiet_cpus(cpus)
        blocker = timing.quiet_cpu_blocker(receipt)
        self.assertIsNotNone(blocker)
        self.assertEqual(blocker["name"], "timing-llc-not-quiet")
        self.assertEqual(blocker["cpu_set"], cpus)
        self.assertEqual(blocker["busy_ticks"], {"9": 2})

    def test_partial_cluster_and_unlisted_task_rollback(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "ledger").mkdir()
            (root / "interrupted").mkdir()
            for job in range(5):
                directory = root / "ledger" / ("%04d" % job)
                directory.mkdir()
                (directory / "receipt.json").write_bytes(str(job).encode("ascii"))
            completed = [{"job": job} for job in (1, 2, 3)]
            partial = timing._rollback_partial_cluster(root, completed)
            self.assertEqual(len(partial), 3)
            self.assertFalse((root / "ledger/0001").exists())
            unlisted = timing._rollback_unlisted_task_entries(
                root, (0,), ())
            self.assertEqual([record["source"] for record in unlisted],
                             ["task-0004.unbound"])
            self.assertTrue((root / "ledger/0000").is_dir())
            self.assertEqual(
                {path.name for path in (root / "interrupted").iterdir()},
                {record["archive"] for record in partial + unlisted})

    def test_interrupted_transaction_recovery(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / ".transactions").mkdir()
            (root / "interrupted").mkdir()
            transaction = root / ".transactions/0007.part.123"
            transaction.mkdir()
            (transaction / "intent.json.part").write_bytes(b"partial")
            records = timing.recover_interrupted_transactions(root)
            self.assertEqual(len(records), 1)
            self.assertEqual(records[0]["source"], "0007.part.123")
            self.assertEqual(records[0]["files"]["intent.json.part"]["size"], 7)
            self.assertEqual(list((root / ".transactions").iterdir()), [])
            self.assertTrue((root / "interrupted" /
                             records[0]["archive"]).is_dir())

    def test_terminal_error_reserves_final_cluster_for_replay(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "ledger").mkdir()
            (root / "interrupted").mkdir()
            completed = []
            for job in range(timing.CLUSTER_SIZE):
                directory = root / "ledger" / ("%04d" % job)
                directory.mkdir()
                (directory / "receipt.json").write_bytes(
                    str(job).encode("ascii"))
                completed.append({"job": job})
            records = timing._rollback_final_cluster_after_error(
                root, completed)
            self.assertEqual(completed, [])
            self.assertEqual(len(records), timing.CLUSTER_SIZE)
            self.assertEqual(list((root / "ledger").iterdir()), [])
            self.assertTrue(all(record["reason"] ==
                                "no terminal thermal binding"
                                for record in records))

    def test_derived_atomic_receipt_recovery(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / "campaign_receipt.json"
            valid = timing.sealed_record("campaign.v1", {"complete": True})
            (root / "campaign_receipt.json.part").write_bytes(
                timing.canonical_json(valid))
            self.assertEqual(timing.recover_derived_atomic(
                path, "campaign.v1", "campaign"),
                "completed-valid-remnant")
            self.assertEqual(timing.load_canonical(path, "campaign"), valid)
            path.unlink()
            (root / "campaign_receipt.json.part").write_bytes(b"truncated")
            self.assertEqual(timing.recover_derived_atomic(
                path, "campaign.v1", "campaign"),
                "discarded-incomplete-remnant")
            self.assertFalse(path.exists())

    def test_speed_gate_reports_empty_common_success_domain(self) -> None:
        gate = timing._speed_gate({
            "cell_count": 2,
            "outcome_counts": {"common-success": 0, "control-only": 1,
                               "candidate-only": 1, "common-failure": 0},
            "ratio_of_sums": None, "bootstrap": None,
        })
        self.assertEqual(gate, {
            "all_fixed_cells_common_success": False,
            "observed_regression_at_most_1_percent": False,
            "one_sided_95_upper_below_1_01": False,
            "speed_gate_passed": False,
        })


if __name__ == "__main__":
    unittest.main()
