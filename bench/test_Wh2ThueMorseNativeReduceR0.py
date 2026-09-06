#!/usr/bin/env python3
"""Synthetic full-panel tests; no actual tables, trace roots or native codec."""
import copy
import hashlib
import importlib.util
from pathlib import Path
import unittest
from unittest import mock


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("tm_native_reducer_test", HERE / "Wh2ThueMorseNativeReduceR0.py")
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


def panels():
    result = []
    for rep, B, scope, comparison, condition in M.panel_roster():
        base = (M.WIDTHS.index(B) * 2048 + 1) * 4096
        addresses = [base + index * 256 * 4096 for index in range(6)]
        sides = (0, 1, 1, 0, 1, 0, 0, 1) if condition % 2 == 0 else (1, 0, 0, 1, 0, 1, 1, 0)
        times = [900 if comparison >= 3 and side == 0 else 1000 for side in sides]
        receipt = dict(native_status=0, left_code=0, right_code=0,
                       handles_count=10 if scope in (2, 3, 4) else 640,
                       handles_sha256="1" * 64, handles_first=128, handles_last=256,
                       handles_min=128, handles_max=256)
        extra = 0 if scope in (5, 6) else None
        result.append([rep, B, scope, comparison, condition, 0, extra, extra, times, addresses, receipt])
    return result


class PanelTests(unittest.TestCase):
    def test_full_counterbalanced_roster_and_clean_controls(self):
        rows = panels()
        samples, unavailable = M.validate_panels(rows)
        self.assertEqual(len(samples), 5040)
        self.assertFalse(unavailable)
        summaries, failures = M.statistics(samples)
        self.assertFalse(failures)
        self.assertEqual(len(summaries), 105)
        for row in summaries:
            self.assertAlmostEqual(row["matched"]["geometric_mean"],
                                   1.0 if row["comparison"] in M.COMPARISONS[:3] else 0.9)

    def test_conditional_bias_cannot_be_pooled_away(self):
        samples, _ = M.validate_panels(panels())
        for rep in range(12):
            samples[(rep, 2, 0, 0, 0)] = 0.03
            samples[(rep, 2, 0, 0, 1)] = -0.03
        summary, failed = M.statistics(samples)
        self.assertIn([2, 0, 0, 0], failed)
        self.assertIn([2, 0, 0, 1], failed)
        self.assertAlmostEqual(summary[0]["matched"]["geometric_mean"], 1)

    def test_unavailable_panel_keeps_unrelated_complete_statistics(self):
        rows = panels()
        rows[4][5], rows[4][8] = 1, None
        rows[4][10].update(left_code=3, handles_count=0, handles_first=None,
                            handles_last=None, handles_min=None, handles_max=None)
        samples, unavailable = M.validate_panels(rows)
        self.assertEqual(len(unavailable), 1)
        summary, _ = M.statistics(samples)
        self.assertIsNone(summary[1]["conditions"][0])
        self.assertIsNone(summary[1]["matched"])
        self.assertIsNotNone(summary[-1]["matched"])

    def test_forgeries_rejected(self):
        mutations = (
            (0, [0], True), (0, [1], 3), (0, [4], 2), (0, [5], 2),
            (0, [8, 0], 0), (0, [8, 0], True), (0, [8, 0], 10**20),
            (0, [9, 0], 4097), (0, [9, 1], 4096),
            (0, [10, "handles_count"], 639), (0, [10, "handles_sha256"], "x" * 64),
            (0, [10, "left_code"], 1),
        )
        for index, path, replacement in mutations:
            rows = panels()
            target = rows[index]
            for key in path[:-1]:
                target = target[key]
            target[path[-1]] = replacement
            with self.subTest(path=path), self.assertRaises(ValueError):
                M.validate_panels(rows)

    def test_same_receive_trace_cannot_change_success_point(self):
        rows = panels()
        row = next(row for row in rows if row[:5] == [1, 2, 5, 0, 1])
        row[6] = 1
        with self.assertRaisesRegex(ValueError, "same receive trace"):
            M.validate_panels(rows)

    def test_physical_ranges_not_just_addresses(self):
        rows = panels()
        row = next(row for row in rows if row[1] == 1280)
        row[9][1] = row[9][0] + 4096
        with self.assertRaisesRegex(ValueError, "spans"):
            M.validate_panels(rows)


class CorrectnessTests(unittest.TestCase):
    def fixture(self):
        data = dict(traces=[dict(B=2, ids=list(range(10)), ranks=[5, 6, 6, 6, 6])],
                    history=[dict(ids=list(range(7)), rank=6)],
                    fixtures=[dict(ids=list(range(6)), rank=5)], history_cases=[[0, 2, 2]],
                    fixture_cases=[[0, 2, 2]], partial_cases=[[0, 2, 1]])
        records = []
        for coord, ids, rank, first, exact_first in M.expected_cases(data):
            extra = first if exact_first else 0
            arm = [0, extra, 0, 0, 0, 0, 0, 6 + extra] if rank == 6 else [1, None, 0, 0, 0, 1, None, len(ids)]
            records.append(coord + [list(arm) for _ in range(3)])
        total = sum(len(case[1]) for case in M.expected_cases(data))
        return data, dict(records=records, cases=4, status=0, candidate_packets=total,
                          control_packets=2 * total, recovered_messages=9)

    def test_rank_and_payload_denominators(self):
        data, raw = self.fixture()
        summary = M.validate_correctness(raw, data)
        self.assertEqual(summary["status_counts"], [[3, 1, 0, 0, 0, 0, 0, 0]] * 3)
        self.assertEqual(summary["fresh_first_success"], [[0, 1, 0, 0, 0, 0]] * 3)

    def test_count_status_types_and_early_success_forgeries(self):
        changes = ((["cases"], 4.0), (["status"], False), (["candidate_packets"], 0),
                   (["control_packets"], 0), (["recovered_messages"], 8),
                   (["records", 0, 4, 1], 0), (["records", 0, 4, 7], 6),
                   (["records", 0, 4, 0], 7), (["records", 0, 4, 2], True))
        for path, value in changes:
            data, raw = self.fixture()
            target = raw
            for key in path[:-1]:
                target = target[key]
            target[path[-1]] = value
            with self.subTest(path=path), self.assertRaises(ValueError):
                M.validate_correctness(raw, data)


class FullRawTests(unittest.TestCase):
    def fixture(self):
        data, correctness = CorrectnessTests().fixture()
        identity_bytes = b"neutral-target-only"
        identity = dict(after_cpu=50, before_cpu=50, requested_cpu=50,
                        canonical_bytes=len(identity_bytes), canonical_hex=identity_bytes.hex(),
                        canonical_sha256=hashlib.sha256(identity_bytes).hexdigest(),
                        capabilities=dict.fromkeys(M.H.TARGET_CAPABILITY_KEYS, 0),
                        derived=dict.fromkeys(M.H.TARGET_DERIVED_KEYS, 0),
                        raw_capture_complete=True, semantic_validation_passed=True,
                        singleton_affinity_verified=True)
        identity["derived"].update(threads_per_core=1, logical_processors_per_package=1)
        measured = dict(callbacks=40320, encoder_creates=368640, decoder_creates=737280,
                        encode_calls=6635520, feed_calls=4423680, recover_calls=737280)
        raw = dict(schema=M.PROTOCOL + ".raw.v1", protocol=M.PROTOCOL, target_cpu=50,
                   gf_runtime=dict(polynomial=0x14d, address=256, ssse3=1, avx2=1, gfni=0, avx512=0),
                   visit_sha256=M.VISIT_SHA, data_sha256=hashlib.sha256(M.H.canonical_bytes(data)).hexdigest(),
                   target_identity_before=identity, target_identity_after=copy.deepcopy(identity),
                   config=M.expected_config(), correctness=correctness, panels=panels(),
                   timed_scope_counts=dict(measured=measured, warm={key: value // 4 for key, value in measured.items()}),
                   elapsed_seconds=5.0, outcome="COMPLETE", diagnostic="")
        return data, raw

    def reduce(self, raw, data):
        with mock.patch.object(M, "sibling", side_effect=AssertionError("no real artifact or codec in neutral tests")):
            return M.reduce_raw(raw, data=data)

    def test_complete_synthetic_result_replays_with_exact_scope_counts(self):
        data, raw = self.fixture()
        result = self.reduce(raw, data)
        self.assertEqual(result["outcome"], "PASS")
        self.assertTrue(result["WH1_compared"])
        self.assertFalse(result["promotion_claimed"])
        self.assertFalse(result["all_K_claimed"])

    def test_control_failure_cannot_become_promotion(self):
        data, raw = self.fixture()
        for row in raw["panels"]:
            if row[1:5] == [2, 0, 0, 0]:
                sides = (0, 1, 1, 0, 1, 0, 0, 1)
                row[8] = [1030 if side == 0 else 1000 for side in sides]
        result = self.reduce(raw, data)
        self.assertEqual(result["outcome"], "CONTROL_FAIL")
        self.assertIn([2, 0, 0, 0], result["failed_controls"])

    def test_configuration_outcome_counts_and_identity_forgeries(self):
        mutations = ((["outcome"], "NONCOMPARABLE"), (["config", "batch_messages"], 16),
                     (["config", "receive_ids", 9], 16), (["data_sha256"], "0" * 64),
                     (["visit_sha256"], "0" * 64), (["target_cpu"], 50.0),
                     (["gf_runtime", "polynomial"], 0x11d),
                     (["target_identity_after", "after_cpu"], 114),
                     (["timed_scope_counts", "measured", "feed_calls"], 4423679),
                     (["timed_scope_counts", "warm", "callbacks"], 10080.0),
                     (["elapsed_seconds"], 0.00001), (["elapsed_seconds"], float("nan")))
        for path, replacement in mutations:
            data, raw = self.fixture()
            target = raw
            for key in path[:-1]:
                target = target[key]
            target[path[-1]] = replacement
            with self.subTest(path=path), self.assertRaises((ValueError, M.H.ValidationError)):
                self.reduce(raw, data)

    def test_one_unavailable_control_cannot_hide_other_call_counts(self):
        data, raw = self.fixture()
        row = raw["panels"][4]
        row[5], row[8], row[10]["left_code"] = 1, None, 3
        row[10].update(handles_count=0, handles_first=None, handles_last=None,
                       handles_min=None, handles_max=None)
        raw["outcome"] = "NONCOMPARABLE"
        self.assertEqual(self.reduce(raw, data)["outcome"], "NONCOMPARABLE")
        raw["timed_scope_counts"]["measured"]["encode_calls"] = 0
        with self.assertRaises(ValueError):
            self.reduce(raw, data)

    def test_impossible_control_error_records_rejected(self):
        invalid = ([2, 0, 5, 0, 0, 0, 0, 0], [3, None, 0, 5, None, None, None, 0],
                   [4, None, 0, 0, 5, None, None, 0], [5, 0, 0, 0, 0, 5, None, 6],
                   [6, None, 0, 0, 0, 0, 5, 6])
        for value in invalid:
            with self.subTest(value=value), self.assertRaises(ValueError):
                M.validate_arm(value, 10)


if __name__ == "__main__":
    unittest.main()
