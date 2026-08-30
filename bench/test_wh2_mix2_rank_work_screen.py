#!/usr/bin/env python3
"""Focused pure tests for the closed mix2 rank/work controller."""

import sys
sys.dont_write_bytecode = True

import copy
import hashlib
import os
from pathlib import Path
import signal
import tempfile
import threading
import time
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parent))
import wh2_mix2_rank_work_screen as subject


FAKE_BINARY = {
    "path": "/frozen/wh2-bench",
    "sha256": "a" * 64,
    "build_id": "b" * 40,
    "embedded_source_commit": "c" * 40,
    "test_hooks_verified": True,
}
FAKE_CONTROLLER = {
    "path": "/frozen/wh2_mix2_rank_work_screen.py",
    "sha256": "d" * 64,
}


def manifest():
    return subject.build_manifest(FAKE_BINARY, FAKE_CONTROLLER)


def _hash(label):
    return hashlib.sha256(label.encode("ascii")).hexdigest()


def _set_status(record, failed):
    record["success"] = 0 if failed else 1
    record["rank_fail"] = 1 if failed else 0
    record["error"] = 0
    record["binary_deficiency"] = 2
    record["gf256_rank_gain"] = 1 if failed else 2
    record["heavy_shortfall"] = 1 if failed else 0
    record["fail_rate"] = "1.00000000" if failed else "0.00000000"
    record["failure_trials"] = "0" if failed else ""
    record["first_rank_fail"] = 0 if failed else -1


def synthetic_records(with_repair=True):
    receipt = manifest()
    rows = []
    for invocation in subject.make_invocations():
        stdout_hash = _hash("stdout-{}".format(invocation.ordinal))
        metadata_hash = _hash("metadata-{}".format(invocation.ordinal))
        argv_hash = receipt["invocations"][invocation.ordinal]["argv_sha256"]
        row_ordinal = 0
        for K in subject.K_VALUES:
            for mix_count in subject.MIX_COUNTS:
                arm = subject.ARM_BY_LAYOUT_MIX[(invocation.layout, mix_count)]
                for overhead in subject.OVERHEADS:
                    row = {
                        "active_packet_peel_seed_xor": "0x0",
                        "arm": arm,
                        "argv_sha256": argv_hash,
                        "attempt_mode": "exact",
                        "binary_build_id": FAKE_BINARY["build_id"],
                        "binary_deficiency": 2,
                        "binary_dense_rows": 12,
                        "binary_sha256": FAKE_BINARY["sha256"],
                        "block_bytes": 1280,
                        "block_xors": 95 if arm == "D" else 100,
                        "construction_attempt": invocation.attempt,
                        "construction_seed_basis": subject.RAW_SEED_BASIS,
                        "controller_sha256": FAKE_CONTROLLER["sha256"],
                        "csv_row_sha256": _hash(
                            "row-{}-{}".format(invocation.ordinal, row_ordinal)),
                        "dense_anchor_layout": invocation.layout,
                        "dense_identity_corner": False,
                        "diagnostic_timings": {
                            field: "1.000" for field in subject.TIMING_FIELDS
                        },
                        "effective_packet_seed":
                            subject._effective_packet_seed(invocation.attempt),
                        "effective_precode_seed":
                            subject._effective_precode_seed(invocation.attempt),
                        "error": 0,
                        "fail_rate": "0.00000000",
                        "failure_trials": "",
                        "first_rank_fail": -1,
                        "gf256_heavy_rows": 12,
                        "gf256_muladds": 99 if arm == "D" else 100,
                        "gf256_rank_gain": 2,
                        "heavy_family": "periodic",
                        "heavy_shortfall": 0,
                        "inactivated_columns": 100,
                        "invocation_ordinal": invocation.ordinal,
                        "K": K,
                        "legacy_seed_attempt": "",
                        "loss_ppm": invocation.loss_ppm,
                        "loss_root": invocation.root,
                        "metadata_sha256": metadata_hash,
                        "mix_count": mix_count,
                        "overhead": overhead,
                        "packet_attempt": invocation.attempt,
                        "precode_attempt": invocation.attempt,
                        "rank_fail": 0,
                        "record_ordinal":
                            invocation.ordinal * subject.ROWS_PER_INVOCATION +
                            row_ordinal,
                        "root_index": invocation.root_index,
                        "row_ordinal": row_ordinal,
                        "schedule": invocation.schedule,
                        "schema": subject.RECORD_SCHEMA,
                        "seed_schedule_sha256":
                            subject.RAW_SEED_SCHEDULE_SHA256,
                        "source_commit": FAKE_BINARY["embedded_source_commit"],
                        "source_hits": 3,
                        "staircase": subject.STAIRCASE_BY_K[K],
                        "stdout_sha256": stdout_hash,
                        "success": 1,
                        "trials": 1,
                    }
                    if (with_repair and arm == "A" and
                            invocation.attempt == 0 and
                            invocation.root_index == 0 and
                            invocation.schedule == "iid" and K == 10000 and
                            overhead == 0):
                        _set_status(row, True)
                    rows.append(row)
                    row_ordinal += 1
    return rows


def _find(rows, arm, attempt=0, root_index=0, schedule="iid", K=10000,
          overhead=0):
    matches = [row for row in rows if
               row["arm"] == arm and
               row["construction_attempt"] == attempt and
               row["root_index"] == root_index and
               row["schedule"] == schedule and row["K"] == K and
               row["overhead"] == overhead]
    assert len(matches) == 1
    return matches[0]


def metadata_line(invocation):
    values = {
        "trials": "1", "threads": "1",
        "loss": "0.10000000000000001" if invocation.loss_ppm == 100000
                else "0.5",
        "seed": "0x{:x}".format(int(invocation.root, 16)),
        "source_hits_override": "0", "packet_peel_seed_xor": "0x0",
        "binary_dense_rows_override": "12",
        "gf256_heavy_rows_override": "12",
        "dense_anchor_layout": invocation.layout,
        "odd_packet_peel_seed_xor": "0x0",
        "packet_row_seed_multiplier": "0x1",
        "packet_row_seed_avalanche": "0",
        "seed_block_bytes_override": "0", "overhead_stream": "paired",
        "full_payload_solve": "0", "schedule": invocation.schedule,
        "cold_solve_wide_xor": "policy", "exact_attempt_mode": "1",
        "exact_precode_attempt": str(invocation.attempt),
        "exact_packet_attempt": str(invocation.attempt),
        "construction_seed_basis": subject.RAW_SEED_BASIS,
        "seed_schedule_sha256": subject.RAW_SEED_SCHEDULE_SHA256,
        "source_git_commit": FAKE_BINARY["embedded_source_commit"],
    }
    return "# precodefail: " + " ".join(
        "{}={}".format(key, values[key]) for key in subject.METADATA_KEYS)


def csv_row(invocation, K, mix_count, overhead):
    values = {field: "" for field in subject.CSV_HEADER}
    values.update({
        "N": str(K), "bb": "1280", "heavy_family": "periodic",
        "mix_count": str(mix_count),
        "staircase": str(subject.STAIRCASE_BY_K[K]),
        "binary_dense_rows": "12", "gf256_heavy_rows": "12",
        "source_hits": "3", "dense_identity_corner": "0",
        "overhead": str(overhead), "trials": "1", "success": "1",
        "rank_fail": "0", "error": "0", "fail_rate": "0.00000000",
        "inact_mu": "100.000", "inact_max": "100",
        "binary_def_mu": "2.000", "binary_def_max": "2",
        "heavy_gain_mu": "2.000", "heavy_gain_min": "2",
        "heavy_shortfall": "0", "solve_ms_mu": "1.000",
        "build_ms_mu": "1.000", "peel_ms_mu": "1.000",
        "project_ms_mu": "1.000", "residual_ms_mu": "1.000",
        "backsub_ms_mu": "1.000", "seed_attempt": "",
        "block_xors_mu": "100.000", "block_muladds_mu": "100.000",
        "first_rank_fail": "-1", "binary_def_hist": "2:1",
        "heavy_gain_hist": "2:1", "failure_trials": "",
        "active_packet_peel_seed_xor": "0x0",
        "precode_attempt": str(invocation.attempt),
        "packet_attempt": str(invocation.attempt), "attempt_mode": "exact",
        "construction_seed_basis": subject.RAW_SEED_BASIS,
        "seed_schedule_sha256": subject.RAW_SEED_SCHEDULE_SHA256,
        "effective_precode_seed":
            subject._effective_precode_seed(invocation.attempt),
        "effective_packet_seed":
            subject._effective_packet_seed(invocation.attempt),
    })
    return ",".join(values[field] for field in subject.CSV_HEADER)


def invocation_stdout(invocation):
    lines = [metadata_line(invocation), ",".join(subject.CSV_HEADER)]
    lines.extend(csv_row(invocation, K, mix_count, overhead)
                 for K in subject.K_VALUES
                 for mix_count in subject.MIX_COUNTS
                 for overhead in subject.OVERHEADS)
    return ("\n".join(lines) + "\n").encode("ascii")


def mutate_first_csv(stdout, field, value):
    lines = stdout.decode("ascii").splitlines()
    row = lines[2].split(",")
    row[subject.CSV_HEADER.index(field)] = value
    lines[2] = ",".join(row)
    return ("\n".join(lines) + "\n").encode("ascii")


def set_failure_block(rows, count_a, count_d):
    keys = []
    for row in rows:
        if row["arm"] == "A" and row["overhead"] == 0 and len(keys) < count_a:
            keys.append((row["construction_attempt"], row["root_index"],
                         row["schedule"], row["K"]))
    assert len(keys) == count_a
    for index, (attempt, root, schedule, K) in enumerate(keys):
        _set_status(_find(rows, "A", attempt, root, schedule, K, 0), True)
        if index < count_d:
            _set_status(_find(rows, "C", attempt, root, schedule, K, 0), True)
            _set_status(_find(rows, "D", attempt, root, schedule, K, 0), True)


class ManifestTests(unittest.TestCase):
    def test_exact_manifest_count_and_allowlisted_argv(self):
        receipt = manifest()
        self.assertEqual(receipt["expected_invocations"], 96)
        self.assertEqual(receipt["expected_records"], 6144)
        self.assertEqual(len(receipt["invocations"]), 96)
        self.assertEqual(
            receipt["policy"]["child_environment"],
            subject.CHILD_ENVIRONMENT)
        self.assertEqual(receipt["manifest_sha256"], subject.sha256_json({
            key: value for key, value in receipt.items()
            if key != "manifest_sha256"
        }))
        forbidden = {
            "--payload-e2e", "--full-payload-solve", "--source-hits",
            "--packet-peel-seed-xor", "--odd-packet-peel-seed-xor",
            "--packet-row-seed-multiplier",
            "--packet-row-seed-avalanche", "--seed-block-bytes",
            "--fail-thread-launch-after",
        }
        for entry in receipt["invocations"]:
            argv = entry["argv"]
            self.assertIn("--binary-dense-rows", argv)
            self.assertEqual(argv[argv.index("--binary-dense-rows") + 1], "12")
            self.assertIn("--gf256-heavy-rows", argv)
            self.assertEqual(argv[argv.index("--gf256-heavy-rows") + 1], "12")
            self.assertEqual(argv[argv.index("--mix-count") + 1], "3,2")
            self.assertEqual(
                argv[argv.index("--cold-solve-wide-xor") + 1], "policy")
            self.assertTrue(forbidden.isdisjoint(argv))

    def test_resigned_manifest_experiment_tamper_is_invalid(self):
        receipt = copy.deepcopy(manifest())
        receipt["experiment"]["block_bytes"] = 2
        receipt = subject._resign_manifest(receipt)
        with self.assertRaises(subject.ScreenError):
            subject.reduce_records(synthetic_records(), receipt)

    def test_binary_test_hook_markers_are_required_prelaunch(self):
        commit = FAKE_BINARY["embedded_source_commit"]
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "worker"
            data = b"\x7fELF" + commit.encode("ascii")
            path.write_bytes(data)
            path.chmod(0o700)
            with mock.patch.object(subject, "_gnu_build_id",
                                   return_value=FAKE_BINARY["build_id"]):
                with self.assertRaisesRegex(subject.ScreenError,
                                            "test-hook markers"):
                    subject.verify_binary(
                        path, subject.sha256_bytes(data),
                        FAKE_BINARY["build_id"], commit)
            data += b"\0" + b"\0".join(
                marker.encode("ascii")
                for marker in subject.REQUIRED_TEST_HOOK_MARKERS)
            path.write_bytes(data)
            with mock.patch.object(subject, "_gnu_build_id",
                                   return_value=FAKE_BINARY["build_id"]):
                verified = subject.verify_binary(
                    path, subject.sha256_bytes(data), FAKE_BINARY["build_id"],
                    commit)
            try:
                self.assertTrue(verified.test_hooks_verified)
            finally:
                os.close(verified.descriptor)

    def test_run_manifest_mismatch_stops_before_first_child(self):
        with tempfile.TemporaryDirectory() as temporary:
            descriptor = os.open("/dev/null", os.O_RDONLY)
            verified = subject.VerifiedBinary(
                Path(FAKE_BINARY["path"]), descriptor, 1, 2, 3,
                FAKE_BINARY["sha256"], FAKE_BINARY["build_id"],
                FAKE_BINARY["embedded_source_commit"], True)
            with mock.patch.object(subject, "verify_binary",
                                   return_value=verified), \
                    mock.patch.object(subject, "_manifest_for_verified",
                                      return_value=manifest()), \
                    mock.patch.object(subject, "_execute_invocation") as run:
                with self.assertRaisesRegex(
                        subject.ScreenError, "preregistered identity"):
                    subject.run_campaign(
                        Path(FAKE_BINARY["path"]), Path(temporary),
                        FAKE_BINARY["sha256"], FAKE_BINARY["build_id"],
                        FAKE_BINARY["embedded_source_commit"], "e" * 64)
                run.assert_not_called()


class ParseTests(unittest.TestCase):
    def setUp(self):
        self.invocation = subject.make_invocations()[0]
        self.command = self.invocation.argv(Path(FAKE_BINARY["path"]))
        self.identity = {
            "build_id": FAKE_BINARY["build_id"],
            "embedded_source_commit": FAKE_BINARY["embedded_source_commit"],
            "sha256": FAKE_BINARY["sha256"],
        }

    def parse(self, stdout):
        return subject.parse_invocation_output(
            self.invocation, self.command, stdout, b"", self.identity,
            FAKE_CONTROLLER["sha256"])

    def test_accepts_exact_metadata_header_and_64_rows(self):
        result = self.parse(invocation_stdout(self.invocation))
        self.assertEqual(len(result.rows), 64)
        self.assertEqual([row["arm"] for row in result.rows[:5]],
                         ["A", "A", "A", "A", "B"])
        self.assertTrue(all(row["legacy_seed_attempt"] == ""
                            for row in result.rows))

    def test_missing_and_duplicate_metadata_are_invalid(self):
        original = invocation_stdout(self.invocation).decode("ascii")
        missing = original.replace(" threads=1", "", 1).encode("ascii")
        duplicate = original.replace(
            "\n", " trials=1\n", 1).encode("ascii")
        for value in (missing, duplicate):
            with self.subTest(value=value[:100]):
                with self.assertRaises(subject.ScreenError):
                    self.parse(value)

    def test_missing_and_duplicate_rows_are_invalid(self):
        lines = invocation_stdout(self.invocation).decode("ascii").splitlines()
        missing = ("\n".join(lines[:-1]) + "\n").encode("ascii")
        duplicate_lines = list(lines)
        duplicate_lines[-1] = duplicate_lines[-2]
        duplicate = ("\n".join(duplicate_lines) + "\n").encode("ascii")
        for value in (missing, duplicate):
            with self.subTest(kind=len(value)):
                with self.assertRaises(subject.ScreenError):
                    self.parse(value)

    def test_payload_blind_argv_and_seed_identity_mutations_are_invalid(self):
        stdout = invocation_stdout(self.invocation)
        with self.assertRaises(subject.ScreenError):
            subject.parse_invocation_output(
                self.invocation, self.command + ["--payload-e2e"], stdout,
                b"", self.identity, FAKE_CONTROLLER["sha256"])
        for field, value in (
                ("seed_attempt", "0"),
                ("precode_attempt", "1"),
                ("packet_attempt", "1"),
                ("effective_precode_seed", "0x0000000000000001"),
                ("effective_packet_seed", "0x00000001")):
            with self.subTest(field=field):
                with self.assertRaises(subject.ScreenError):
                    self.parse(mutate_first_csv(stdout, field, value))

    def test_inactivation_cap_accepts_1024_and_rejects_1025(self):
        stdout = invocation_stdout(self.invocation)
        at_cap = mutate_first_csv(stdout, "inact_mu", "1024.000")
        at_cap = mutate_first_csv(at_cap, "inact_max", "1024")
        self.assertEqual(self.parse(at_cap).rows[0]["inactivated_columns"], 1024)
        over = mutate_first_csv(stdout, "inact_mu", "1025.000")
        over = mutate_first_csv(over, "inact_max", "1025")
        with self.assertRaises(subject.ScreenError):
            self.parse(over)


class ReducerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.receipt = manifest()

    def test_pass_at_ratio_boundaries_and_one_repair(self):
        decision = subject.reduce_records(
            synthetic_records(), self.receipt)
        self.assertEqual(decision["disposition"], "PASS")
        self.assertEqual(decision["failure_totals"]["A"], 1)
        self.assertEqual(decision["failure_totals"]["D"], 0)
        self.assertEqual(decision["repair_count_vs_A"], 1)
        global_work = decision["common_success_work"]["global"]
        self.assertTrue(all(metric["pass"] for metric in global_work.values()))
        self.assertEqual(global_work["block_xors"]["limit_numerator"], 95)
        self.assertEqual(global_work["gf256_muladds"]["limit_numerator"], 99)
        self.assertEqual(
            global_work["inactivated_columns"]["limit_numerator"], 100)

    def test_no_control_failures_is_inconclusive(self):
        decision = subject.reduce_records(
            synthetic_records(with_repair=False), self.receipt)
        self.assertEqual(decision["disposition"], "INCONCLUSIVE")

    def test_failure_introduction_rejects(self):
        rows = synthetic_records()
        _set_status(_find(rows, "D", K=10001), True)
        decision = subject.reduce_records(rows, self.receipt)
        self.assertEqual(decision["disposition"], "REJECT")
        self.assertEqual(decision["introduction_counts"],
                         {"vs_A": 1, "vs_C": 1})

    def test_introductions_vs_A_and_C_are_independent(self):
        rows = synthetic_records()
        _set_status(_find(rows, "C", K=10001), True)
        _set_status(_find(rows, "D", K=10001), True)
        decision = subject.reduce_records(rows, self.receipt)
        self.assertEqual(decision["introduction_counts"],
                         {"vs_A": 1, "vs_C": 0})
        rows = synthetic_records()
        _set_status(_find(rows, "A", K=10001), True)
        _set_status(_find(rows, "D", K=10001), True)
        decision = subject.reduce_records(rows, self.receipt)
        self.assertEqual(decision["introduction_counts"],
                         {"vs_A": 0, "vs_C": 1})

    def test_overhead_nesting_rejects(self):
        rows = synthetic_records()
        for arm in ("A", "C", "D"):
            _set_status(_find(rows, arm, K=10001, overhead=4), True)
        decision = subject.reduce_records(rows, self.receipt)
        self.assertEqual(decision["disposition"], "REJECT")
        self.assertFalse(decision["gates"]["overhead_nesting"])

    def test_exact_integer_five_percent_recovery_boundary_passes(self):
        rows = synthetic_records(with_repair=False)
        set_failure_block(rows, 20, 19)
        decision = subject.reduce_records(rows, self.receipt)
        self.assertEqual(decision["failure_totals"]["A"], 20)
        self.assertEqual(decision["failure_totals"]["D"], 19)
        self.assertTrue(decision["gates"]["five_percent_reduction_vs_A"])
        self.assertEqual(decision["disposition"], "PASS")

    def test_integer_five_percent_one_unit_miss_rejects(self):
        rows = synthetic_records(with_repair=False)
        set_failure_block(rows, 21, 20)
        decision = subject.reduce_records(rows, self.receipt)
        self.assertEqual(decision["failure_totals"]["A"], 21)
        self.assertEqual(decision["failure_totals"]["D"], 20)
        self.assertFalse(decision["gates"]["five_percent_reduction_vs_A"])
        self.assertEqual(decision["disposition"], "REJECT")

    def test_work_ratio_failure_and_zero_denominator_rules(self):
        rows = synthetic_records()
        for row in rows:
            if row["arm"] == "D":
                row["block_xors"] = 96
        decision = subject.reduce_records(rows, self.receipt)
        self.assertEqual(decision["disposition"], "REJECT")
        self.assertFalse(decision["gates"]["common_success_global"])

        rows = synthetic_records()
        for row in rows:
            if row["arm"] in ("A", "D"):
                row["block_xors"] = 0
        zero = subject.reduce_records(rows, self.receipt)
        self.assertEqual(zero["disposition"], "PASS")
        self.assertEqual(
            zero["common_success_work"]["global"]["block_xors"]
                ["zero_denominator_rule"], "both-zero-pass")
        _find(rows, "D", K=10001)["block_xors"] = 1
        positive = subject.reduce_records(rows, self.receipt)
        self.assertEqual(positive["disposition"], "REJECT")

    def test_empty_K_and_schedule_common_success_reject(self):
        rows = synthetic_records()
        for row in rows:
            if (row["K"] == 10001 or row["schedule"] == "burst") and \
                    row["arm"] in ("A", "C", "D"):
                _set_status(row, True)
        decision = subject.reduce_records(rows, self.receipt)
        self.assertFalse(decision["gates"]["common_success_per_K_schedule"])
        k_entry = next(value for value in
                       decision["common_success_work"]["per_stratum"]["K"]
                       if value["value"] == 10001)
        schedule_entry = next(value for value in
                              decision["common_success_work"]["per_stratum"]
                                      ["schedule"]
                              if value["value"] == "burst")
        self.assertFalse(k_entry["pass"])
        self.assertFalse(schedule_entry["pass"])

    def test_per_stratum_101_boundary_and_one_unit_miss(self):
        for dimension, selected in (("K", 10001), ("schedule", "burst")):
            rows = synthetic_records()
            pairs = [row for row in rows
                     if row[dimension] == selected and row["arm"] in ("A", "D")]
            for row in pairs:
                row["block_xors"] = 0
            control = next(row for row in pairs if row["arm"] == "A")
            candidate = _find(
                rows, "D", control["construction_attempt"],
                control["root_index"], control["schedule"], control["K"],
                control["overhead"])
            control["block_xors"] = 100
            candidate["block_xors"] = 101
            boundary = subject.reduce_records(rows, self.receipt)
            entry = next(value for value in
                         boundary["common_success_work"]["per_stratum"]
                                 [dimension]
                         if value["value"] == selected)
            self.assertTrue(entry["metrics"]["block_xors"]["pass"])
            candidate["block_xors"] = 102
            miss = subject.reduce_records(rows, self.receipt)
            entry = next(value for value in
                         miss["common_success_work"]["per_stratum"][dimension]
                         if value["value"] == selected)
            self.assertFalse(entry["metrics"]["block_xors"]["pass"])

    def test_timing_mutation_preserves_decision_hash_and_disposition(self):
        rows = synthetic_records()
        before = subject.reduce_records(rows, self.receipt)
        changed = list(rows)
        changed[1] = dict(changed[1])
        changed[1]["diagnostic_timings"] = dict(
            changed[1]["diagnostic_timings"])
        changed[1]["diagnostic_timings"]["solve_ms_mu"] = "999.999"
        after = subject.reduce_records(changed, self.receipt)
        self.assertEqual(before["decision_sha256"], after["decision_sha256"])
        self.assertEqual(before["disposition"], after["disposition"])

    def test_missing_duplicate_and_tampered_invocation_identity_invalid(self):
        rows = synthetic_records()
        with self.assertRaises(subject.ScreenError):
            subject.reduce_records(rows[:-1], self.receipt)
        duplicate = list(rows)
        duplicate[-1] = duplicate[-2]
        with self.assertRaises(subject.ScreenError):
            subject.reduce_records(duplicate, self.receipt)
        for field, value in (
                ("dense_anchor_layout", "four0369"),
                ("loss_root", subject.ROOTS[1]),
                ("schedule", "burst"),
                ("construction_attempt", 1)):
            mutated = list(rows)
            mutated[0] = dict(mutated[0])
            mutated[0][field] = value
            with self.subTest(field=field):
                with self.assertRaises(subject.ScreenError):
                    subject.reduce_records(mutated, self.receipt)

    def test_row_stdout_receipt_swap_is_invalid(self):
        rows = synthetic_records()
        rows[0] = dict(rows[0])
        rows[0]["stdout_sha256"] = rows[64]["stdout_sha256"]
        with self.assertRaises(subject.ScreenError):
            subject.reduce_records(rows, self.receipt)

    def test_boolean_integer_identity_is_invalid(self):
        rows = synthetic_records()
        rows[0] = dict(rows[0])
        rows[0]["construction_attempt"] = False
        with self.assertRaises(subject.ScreenError):
            subject.reduce_records(rows, self.receipt)


class PublicationTests(unittest.TestCase):
    def test_pair_success_and_no_clobber(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            subject._publish_pair(root, b"result\n", b"summary\n")
            self.assertEqual((root / subject.RESULT_NAME).read_bytes(), b"result\n")
            self.assertEqual((root / subject.SUMMARY_NAME).read_bytes(), b"summary\n")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / subject.RESULT_NAME
            target.write_bytes(b"keep\n")
            with self.assertRaises(subject.ScreenError):
                subject._publish_pair(root, b"replace\n", b"summary\n")
            self.assertEqual(target.read_bytes(), b"keep\n")
            self.assertFalse((root / subject.SUMMARY_NAME).exists())

    def test_output_directory_swap_and_publication_interrupt_roll_back(self):
        with tempfile.TemporaryDirectory() as temporary:
            parent = Path(temporary)
            root = parent / "out"
            root.mkdir()
            descriptor, pinned, identity = subject._open_output_directory(root)
            moved = parent / "moved"
            root.rename(moved)
            root.mkdir()
            try:
                with self.assertRaises(subject.ScreenError):
                    subject._publish_pair_at(
                        descriptor, pinned, identity, b"r\n", b"s\n")
            finally:
                os.close(descriptor)
            self.assertFalse((root / subject.RESULT_NAME).exists())
            self.assertFalse((moved / subject.RESULT_NAME).exists())
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            real_link = os.link
            calls = [0]

            def interrupted_link(*args, **kwargs):
                calls[0] += 1
                if calls[0] == 2:
                    raise KeyboardInterrupt()
                return real_link(*args, **kwargs)

            with mock.patch.object(subject.os, "link", side_effect=interrupted_link):
                with self.assertRaises(KeyboardInterrupt):
                    subject._publish_pair(root, b"r\n", b"s\n")
            self.assertFalse((root / subject.RESULT_NAME).exists())
            self.assertFalse((root / subject.SUMMARY_NAME).exists())

    def test_stage_write_and_fsync_failures_leave_no_temporary(self):
        for operation in ("write", "fsync"):
            with self.subTest(operation=operation), \
                    tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                descriptor, _path, _identity = \
                    subject._open_output_directory(root)
                try:
                    with mock.patch.object(
                            subject.os, operation,
                            side_effect=OSError("injected {}".format(
                                operation))):
                        with self.assertRaises(OSError):
                            subject._stage_file(
                                descriptor, subject.RESULT_NAME, b"data\n")
                finally:
                    os.close(descriptor)
                self.assertEqual(list(root.iterdir()), [])

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            real_link = os.link

            def link_then_interrupt(*args, **kwargs):
                real_link(*args, **kwargs)
                raise KeyboardInterrupt()

            with mock.patch.object(
                    subject.os, "link", side_effect=link_then_interrupt):
                with self.assertRaises(KeyboardInterrupt):
                    subject._publish_pair(root, b"r\n", b"s\n")
            self.assertFalse((root / subject.RESULT_NAME).exists())
            self.assertFalse((root / subject.SUMMARY_NAME).exists())


class ProcessLifecycleTests(unittest.TestCase):
    def test_sigterm_kills_active_process_group(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            pid_path = root / "child.pid"
            worker = root / "worker.py"
            worker.write_text(
                "#!/usr/bin/python3\n"
                "import os\n"
                "import time\n"
                "assert os.environ['LANG'] == 'C'\n"
                "assert os.environ['LC_ALL'] == 'C'\n"
                "assert os.environ['TZ'] == 'UTC'\n"
                "assert not any(key.startswith('LD_') for key in os.environ)\n"
                "with open({!r}, 'w') as stream:\n"
                "    stream.write(str(os.getpid()))\n"
                "time.sleep(60)\n".format(str(pid_path)),
                encoding="ascii")
            worker.chmod(0o700)
            descriptor = os.open(str(worker), os.O_RDONLY)
            verified = subject.VerifiedBinary(
                worker, descriptor, 1, 2, 3, "a" * 64, "b" * 40,
                "c" * 40, True)
            sender_errors = []

            def terminate_when_started():
                deadline = time.monotonic() + 5.0
                while time.monotonic() < deadline and not pid_path.exists():
                    time.sleep(0.01)
                if not pid_path.exists():
                    sender_errors.append("child did not start")
                    return
                os.kill(os.getpid(), signal.SIGTERM)

            previous = subject._install_termination_handlers()
            sender = threading.Thread(target=terminate_when_started)
            sender.start()
            try:
                with self.assertRaisesRegex(subject.ScreenError, "SIGTERM"):
                    subject._execute_invocation(
                        verified, subject.make_invocations()[0], 10.0)
            finally:
                subject._restore_termination_handlers(previous)
                sender.join(timeout=5.0)
                os.close(descriptor)
            self.assertFalse(sender.is_alive())
            self.assertEqual(sender_errors, [])
            child_pid = int(pid_path.read_text(encoding="ascii"))
            with self.assertRaises(ProcessLookupError):
                os.kill(child_pid, 0)


if __name__ == "__main__":
    unittest.main()
